"""Tests for pubmatrix.core."""

import urllib.parse
import xml.etree.ElementTree as ET
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import requests

from pubmatrix.core import (
    _build_base_url,
    _build_hyperlink_url,
    _export_csv,
    _export_ods,
    _extract_count,
    _fetch_count,
    _make_hyperlink_formula,
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

    def test_infinite_value_raises(self):
        with pytest.raises(ValueError, match="finite"):
            _validate_daterange([2010, float("inf")])

    def test_nan_value_raises(self):
        with pytest.raises(ValueError, match="finite"):
            _validate_daterange([float("nan"), 2020])


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

    def test_url_has_retmax_zero_and_no_usehistory(self):
        with patch("pubmatrix.core.requests.get") as mock_get:
            mock_get.return_value = self._mock_response(VALID_XML)
            _fetch_count("http://fake", "term")
        url = mock_get.call_args[0][0]
        assert "retmax=0" in url
        assert "usehistory" not in url

    def test_backoff_is_exponential(self):
        with patch("pubmatrix.core.requests.get") as mock_get, \
             patch("pubmatrix.core.time.sleep") as mock_sleep:
            mock_get.side_effect = requests.ConnectionError("down")
            with pytest.raises(RuntimeError):
                _fetch_count("http://fake", "term", n_tries=3)
        assert [c.args[0] for c in mock_sleep.call_args_list] == [0.25, 0.5]

    def test_cache_hit_skips_network(self, tmp_path):
        import json
        url = "http://fake&term=x&retmax=0"
        cache_file = tmp_path / (
            __import__("hashlib").md5(url.encode()).hexdigest() + ".json"
        )
        cache_file.write_text(json.dumps({"count": 99, "url": url}))

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

    def test_whitespace_terms_B(self):
        with pytest.raises(ValueError, match="B contains empty or whitespace"):
            pubmatrix(["a"], ["   "])

    def test_n_tries_zero_raises(self):
        with pytest.raises(ValueError, match="n_tries"):
            pubmatrix(["a"], ["b"], n_tries=0)

    def test_n_workers_zero_raises(self):
        with pytest.raises(ValueError, match="n_workers"):
            pubmatrix(["a"], ["b"], n_workers=0)

    def test_duplicate_A_raises(self):
        with pytest.raises(ValueError, match=r"A must not contain duplicate"):
            pubmatrix(["a", "a"], ["b"])

    def test_duplicate_B_raises(self):
        with pytest.raises(ValueError, match=r"B must not contain duplicate"):
            pubmatrix(["a"], ["b", "b"])

    def test_duplicate_error_names_term(self):
        with pytest.raises(ValueError, match="'a'"):
            pubmatrix(["a", "a"], ["b"])

    def test_duplicates_detected_after_stripping(self):
        with pytest.raises(ValueError, match="duplicate"):
            pubmatrix(["a", " a "], ["b"])


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

    def test_query_terms_use_A_and_B_order(self):
        with patch("pubmatrix.core._fetch_count", return_value=5) as mock_fetch, \
             patch("pubmatrix.core.time.sleep"):
            pubmatrix(["A1", "A2"], ["B1", "B2"])
        queried = [
            urllib.parse.unquote(call.args[1]) for call in mock_fetch.call_args_list
        ]
        assert queried == ["A1 AND B1", "A2 AND B1", "A1 AND B2", "A2 AND B2"]


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

    def test_blank_lines_ignored(self, tmp_path):
        f = tmp_path / "blanks.txt"
        f.write_text("A1\n\n  \nA2\n#\n\nB1\nB2\n\n")
        with patch("pubmatrix.core._fetch_count", return_value=3), \
             patch("pubmatrix.core.time.sleep"):
            df = pubmatrix_from_file(str(f))
        assert list(df.columns) == ["A1", "A2"]
        assert list(df.index) == ["B1", "B2"]

    def test_missing_separator_raises(self, tmp_path):
        f = tmp_path / "bad.txt"
        f.write_text("A1\nA2\nB1\n")
        with pytest.raises(ValueError, match="bad.txt"):
            pubmatrix_from_file(str(f))

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            pubmatrix_from_file(str(tmp_path / "missing.txt"))

    def test_no_terms_before_separator_raises(self, tmp_path):
        f = tmp_path / "empty_a.txt"
        f.write_text("#\nB1\n")
        with pytest.raises(ValueError, match="before the '#' separator in:.*empty_a.txt"):
            pubmatrix_from_file(str(f))

    def test_no_terms_after_separator_raises(self, tmp_path):
        f = tmp_path / "empty_b.txt"
        f.write_text("A1\n#\n")
        with pytest.raises(ValueError, match="after the '#' separator in:.*empty_b.txt"):
            pubmatrix_from_file(str(f))


# ---------------------------------------------------------------------------
# URL builders
# ---------------------------------------------------------------------------

class TestBuildBaseUrl:
    def test_minimal_url(self):
        url = _build_base_url("pubmed", None, None)
        assert "db=pubmed" in url
        assert "api_key" not in url
        assert "mindate" not in url

    def test_api_key_included(self):
        url = _build_base_url("pubmed", "SECRET", None)
        assert "api_key=SECRET" in url

    def test_daterange_included(self):
        url = _build_base_url("pmc", None, (2010, 2020))
        assert "mindate=2010" in url
        assert "maxdate=2020" in url
        assert "datetype=pdat" in url

    def test_api_key_and_daterange_together(self):
        url = _build_base_url("pubmed", "KEY", (2015, 2016))
        assert "api_key=KEY" in url and "mindate=2015" in url


class TestBuildHyperlinkUrl:
    def test_database_and_term_encoded(self):
        url = _build_hyperlink_url("pubmed", "WNT1 AND obesity")
        assert url.startswith("https://www.ncbi.nlm.nih.gov/pubmed/?term=")
        assert "WNT1%20AND%20obesity" in url

    def test_pmc_database(self):
        assert "/pmc/" in _build_hyperlink_url("pmc", "x")


class TestMakeHyperlinkFormula:
    def test_formula_shape(self):
        assert _make_hyperlink_formula("http://x", 7) == '=HYPERLINK("http://x","7")'


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def export_df():
    return pd.DataFrame(
        [[10, 5], [3, 8]], index=["B1", "B2"], columns=["A1", "A2"]
    )


class TestExportCsv:
    def test_writes_csv_next_to_outfile(self, export_df, tmp_path):
        out = tmp_path / "results"
        _export_csv(export_df, str(out), "pubmed")
        assert (tmp_path / "results.csv").exists()

    def test_suffix_not_doubled(self, export_df, tmp_path):
        out = tmp_path / "results.csv"
        _export_csv(export_df, str(out), "pubmed")
        assert (tmp_path / "results.csv").exists()
        assert not (tmp_path / "results.csv.csv").exists()

    def test_cells_contain_hyperlink_formula(self, export_df, tmp_path):
        out = tmp_path / "r"
        _export_csv(export_df, str(out), "pubmed")
        text = (tmp_path / "r.csv").read_text()
        assert "=HYPERLINK(" in text
        assert "WNT" not in text  # sanity: only our terms appear

    def test_hyperlink_term_is_A_and_B_order(self, export_df, tmp_path):
        out = tmp_path / "r"
        _export_csv(export_df, str(out), "pubmed")
        text = (tmp_path / "r.csv").read_text()
        assert "A1%20AND%20B1" in text
        assert "B1%20AND%20A1" not in text

    def test_counts_preserved_in_formula(self, export_df, tmp_path):
        out = tmp_path / "r"
        _export_csv(export_df, str(out), "pubmed")
        text = (tmp_path / "r.csv").read_text()
        # to_csv doubles the quotes inside each formula
        for count in ("10", "5", "3", "8"):
            assert f'""{count}"")' in text

    def test_row_and_column_labels_written(self, export_df, tmp_path):
        out = tmp_path / "r"
        _export_csv(export_df, str(out), "pubmed")
        text = (tmp_path / "r.csv").read_text()
        assert "A1" in text and "B1" in text

    def test_pmc_database_in_urls(self, export_df, tmp_path):
        out = tmp_path / "r"
        _export_csv(export_df, str(out), "pmc")
        assert "/pmc/?term=" in (tmp_path / "r.csv").read_text()


class TestExportOds:
    def test_writes_ods_file(self, export_df, tmp_path):
        pytest.importorskip("odf")
        out = tmp_path / "results"
        _export_ods(export_df, str(out), "pubmed")
        assert (tmp_path / "results.ods").exists()

    def test_suffix_not_doubled(self, export_df, tmp_path):
        pytest.importorskip("odf")
        out = tmp_path / "results.ods"
        _export_ods(export_df, str(out), "pubmed")
        assert (tmp_path / "results.ods").exists()
        assert not (tmp_path / "results.ods.ods").exists()

    def test_cells_carry_links_and_counts(self, export_df, tmp_path):
        pytest.importorskip("odf")
        from odf.opendocument import load
        from odf.table import Table, TableRow
        from odf.text import A as OdfA

        out = tmp_path / "r"
        _export_ods(export_df, str(out), "pubmed")

        doc = load(str(tmp_path / "r.ods"))
        links = doc.getElementsByType(OdfA)
        assert len(links) == export_df.size

        hrefs = [
            link.getAttribute("href") for link in links
        ]
        assert any("A1%20AND%20B1" in h for h in hrefs)

        # header row plus one row per B term
        table = doc.getElementsByType(Table)[0]
        assert len(table.getElementsByType(TableRow)) == len(export_df.index) + 1

    def test_pmc_database_in_urls(self, export_df, tmp_path):
        pytest.importorskip("odf")
        from odf.opendocument import load
        from odf.text import A as OdfA

        out = tmp_path / "r"
        _export_ods(export_df, str(out), "pmc")
        doc = load(str(tmp_path / "r.ods"))
        hrefs = [
            link.getAttribute("href") for link in doc.getElementsByType(OdfA)
        ]
        assert all("/pmc/?term=" in h for h in hrefs)


# ---------------------------------------------------------------------------
# pubmatrix() — export wiring and concurrent path
# ---------------------------------------------------------------------------

class TestPubmatrixExportWiring:
    def test_csv_export_written(self, tmp_path):
        out = tmp_path / "wired"
        with patch("pubmatrix.core._fetch_count", return_value=4), \
             patch("pubmatrix.core.time.sleep"):
            pubmatrix(["A1"], ["B1"], outfile=str(out), export_format="csv")
        assert (tmp_path / "wired.csv").exists()

    def test_ods_export_written(self, tmp_path):
        pytest.importorskip("odf")
        out = tmp_path / "wired"
        with patch("pubmatrix.core._fetch_count", return_value=4), \
             patch("pubmatrix.core.time.sleep"):
            pubmatrix(["A1"], ["B1"], outfile=str(out), export_format="ods")
        assert (tmp_path / "wired.ods").exists()

    def test_no_export_when_format_none(self, tmp_path):
        out = tmp_path / "none"
        with patch("pubmatrix.core._fetch_count", return_value=4), \
             patch("pubmatrix.core.time.sleep"):
            pubmatrix(["A1"], ["B1"], outfile=str(out))
        assert list(tmp_path.iterdir()) == []


class TestPubmatrixConcurrent:
    def test_concurrent_path_matches_serial(self):
        counts = {"A1 AND B1": 1, "A2 AND B1": 2, "A1 AND B2": 3, "A2 AND B2": 4}

        def fake_fetch(base_url, encoded, *args, **kwargs):
            return counts[urllib.parse.unquote(encoded)]

        with patch("pubmatrix.core._fetch_count", side_effect=fake_fetch), \
             patch("pubmatrix.core.time.sleep"):
            serial = pubmatrix(["A1", "A2"], ["B1", "B2"], n_workers=1)
            concurrent = pubmatrix(["A1", "A2"], ["B1", "B2"], n_workers=4)

        pd.testing.assert_frame_equal(serial, concurrent)

    def test_concurrent_results_land_in_right_cells(self):
        counts = {"A1 AND B1": 1, "A2 AND B1": 2, "A1 AND B2": 3, "A2 AND B2": 4}

        def fake_fetch(base_url, encoded, *args, **kwargs):
            return counts[urllib.parse.unquote(encoded)]

        with patch("pubmatrix.core._fetch_count", side_effect=fake_fetch), \
             patch("pubmatrix.core.time.sleep"):
            df = pubmatrix(["A1", "A2"], ["B1", "B2"], n_workers=4)

        assert df.loc["B1", "A1"] == 1
        assert df.loc["B1", "A2"] == 2
        assert df.loc["B2", "A1"] == 3
        assert df.loc["B2", "A2"] == 4

    def test_concurrent_creates_cache_dir(self, tmp_path):
        cache = tmp_path / "nested" / "cache"
        with patch("pubmatrix.core._fetch_count", return_value=1), \
             patch("pubmatrix.core.time.sleep"):
            pubmatrix(["A1"], ["B1"], n_workers=2, cache_dir=str(cache))
        assert cache.is_dir()

    def test_worker_exception_propagates(self):
        with patch("pubmatrix.core._fetch_count",
                   side_effect=RuntimeError("boom")), \
             patch("pubmatrix.core.time.sleep"):
            with pytest.raises(RuntimeError, match="boom"):
                pubmatrix(["A1", "A2"], ["B1"], n_workers=2)
