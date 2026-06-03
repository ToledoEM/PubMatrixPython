"""Tests for pubmatrix.core."""

import xml.etree.ElementTree as ET
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import requests

from pubmatrix.core import (
    _build_base_url,
    _extract_count,
    _fetch_count,
    _validate_daterange,
    pubmatrix,
    pubmatrix_from_file,
)
from tests.conftest import (
    MALFORMED_XML,
    MISSING_COUNT_XML,
    NON_NUMERIC_COUNT_XML,
    VALID_XML,
    ZERO_COUNT_XML,
)


# ---------------------------------------------------------------------------
# _extract_count
# ---------------------------------------------------------------------------

class TestExtractCount:
    def test_valid_xml(self):
        assert _extract_count(VALID_XML) == 42

    def test_zero_count(self):
        assert _extract_count(ZERO_COUNT_XML) == 0

    def test_malformed_xml_raises(self):
        with pytest.raises(ValueError, match="Could not parse"):
            _extract_count(MALFORMED_XML)

    def test_missing_count_element_raises(self):
        with pytest.raises(ValueError, match="missing <Count>"):
            _extract_count(MISSING_COUNT_XML)

    def test_non_numeric_count_raises(self):
        with pytest.raises(ValueError, match="not numeric"):
            _extract_count(NON_NUMERIC_COUNT_XML)


# ---------------------------------------------------------------------------
# _validate_daterange
# ---------------------------------------------------------------------------

class TestValidateDaterange:
    def test_none_returns_none(self):
        assert _validate_daterange(None) is None

    def test_valid_range(self):
        assert _validate_daterange([2010, 2020]) == (2010, 2020)

    def test_same_year(self):
        assert _validate_daterange([2015, 2015]) == (2015, 2015)

    def test_floats_rounded(self):
        assert _validate_daterange([2010.4, 2020.6]) == (2010, 2021)

    def test_wrong_length_raises(self):
        with pytest.raises(ValueError, match="exactly 2"):
            _validate_daterange([2010])

    def test_start_after_end_raises(self):
        with pytest.raises(ValueError, match="must be <="):
            _validate_daterange([2020, 2010])


# ---------------------------------------------------------------------------
# _fetch_count
# ---------------------------------------------------------------------------

class TestFetchCount:
    def _mock_response(self, xml_text, status=200):
        mock = MagicMock()
        mock.text = xml_text
        mock.raise_for_status = MagicMock()
        if status != 200:
            mock.raise_for_status.side_effect = requests.HTTPError(f"HTTP {status}")
        return mock

    def test_success_first_attempt(self):
        with patch("pubmatrix.core.requests.get") as mock_get:
            mock_get.return_value = self._mock_response(VALID_XML)
            result = _fetch_count("http://fake", "term", n_tries=2)
        assert result == 42
        assert mock_get.call_count == 1

    def test_retries_on_network_error(self):
        with patch("pubmatrix.core.requests.get") as mock_get, \
             patch("pubmatrix.core.time.sleep"):
            mock_get.side_effect = [
                requests.ConnectionError("timeout"),
                self._mock_response(VALID_XML),
            ]
            result = _fetch_count("http://fake", "term", n_tries=2)
        assert result == 42
        assert mock_get.call_count == 2

    def test_all_retries_exhausted_raises(self):
        with patch("pubmatrix.core.requests.get") as mock_get, \
             patch("pubmatrix.core.time.sleep"):
            mock_get.side_effect = requests.ConnectionError("down")
            with pytest.raises(RuntimeError, match="Failed to fetch count"):
                _fetch_count("http://fake", "term", n_tries=2)

    def test_cache_hit_skips_network(self, tmp_path):
        import json
        cache_file = tmp_path / (
            __import__("hashlib").md5(b"http://fake&term=x&usehistory=y").hexdigest() + ".json"
        )
        cache_file.write_text(json.dumps({"count": 99, "url": "http://fake&term=x&usehistory=y"}))

        with patch("pubmatrix.core.requests.get") as mock_get:
            result = _fetch_count("http://fake", "x", cache_dir=tmp_path)
        assert result == 99
        mock_get.assert_not_called()

    def test_cache_miss_writes_file(self, tmp_path):
        with patch("pubmatrix.core.requests.get") as mock_get, \
             patch("pubmatrix.core.time.sleep"):
            mock_get.return_value = self._mock_response(VALID_XML)
            result = _fetch_count("http://fake", "term", cache_dir=tmp_path)

        assert result == 42
        cache_files = list(tmp_path.iterdir())
        assert len(cache_files) == 1


# ---------------------------------------------------------------------------
# pubmatrix() — validation
# ---------------------------------------------------------------------------

class TestPubmatrixValidation:
    def test_invalid_database(self):
        with pytest.raises(ValueError, match="database must be"):
            pubmatrix(["a"], ["b"], database="notadb")

    def test_invalid_export_format(self):
        with pytest.raises(ValueError, match="export_format must be"):
            pubmatrix(["a"], ["b"], export_format="xlsx")

    def test_export_format_without_outfile(self):
        with pytest.raises(ValueError, match="outfile must be"):
            pubmatrix(["a"], ["b"], export_format="csv")

    def test_empty_A_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            pubmatrix([], ["b"])

    def test_empty_B_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            pubmatrix(["a"], [])

    def test_whitespace_terms_A(self):
        with pytest.raises(ValueError, match="empty or whitespace"):
            pubmatrix(["   "], ["b"])

    def test_n_tries_zero_raises(self):
        with pytest.raises(ValueError, match="n_tries"):
            pubmatrix(["a"], ["b"], n_tries=0)

    def test_n_workers_zero_raises(self):
        with pytest.raises(ValueError, match="n_workers"):
            pubmatrix(["a"], ["b"], n_workers=0)


# ---------------------------------------------------------------------------
# pubmatrix() — DataFrame structure (mocked network)
# ---------------------------------------------------------------------------

class TestPubmatrixDataFrame:
    def _mock_fetch(self, count=5):
        with patch("pubmatrix.core._fetch_count", return_value=count), \
             patch("pubmatrix.core.time.sleep"):
            return pubmatrix(["A1", "A2"], ["B1", "B2"])

    def test_returns_dataframe(self):
        df = self._mock_fetch()
        assert isinstance(df, pd.DataFrame)

    def test_shape(self):
        df = self._mock_fetch()
        assert df.shape == (2, 2)

    def test_columns_match_A(self):
        df = self._mock_fetch()
        assert list(df.columns) == ["A1", "A2"]

    def test_index_matches_B(self):
        df = self._mock_fetch()
        assert list(df.index) == ["B1", "B2"]

    def test_values_populated(self):
        df = self._mock_fetch(count=7)
        assert (df.values == 7).all()


# ---------------------------------------------------------------------------
# pubmatrix_from_file
# ---------------------------------------------------------------------------

class TestPubmatrixFromFile:
    def test_basic_parsing(self, tmp_path):
        f = tmp_path / "terms.txt"
        f.write_text("A1\nA2\n#\nB1\nB2\n")
        with patch("pubmatrix.core._fetch_count", return_value=3), \
             patch("pubmatrix.core.time.sleep"):
            df = pubmatrix_from_file(str(f))
        assert df.shape == (2, 2)

    def test_missing_separator_raises(self, tmp_path):
        f = tmp_path / "bad.txt"
        f.write_text("A1\nA2\nB1\n")
        with pytest.raises(ValueError, match="'#' line"):
            pubmatrix_from_file(str(f))

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            pubmatrix_from_file(str(tmp_path / "missing.txt"))

    def test_empty_section_raises(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("#\nB1\n")
        with pytest.raises(ValueError, match="terms both before"):
            pubmatrix_from_file(str(f))
