"""Tests for pubmatrix.heatmap."""

import warnings

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")  # non-interactive backend for tests

from matplotlib.figure import Figure
from pubmatrix.heatmap import (
    _auto_font_size,
    _clustered_order,
    _handle_na,
    _overlap_percentage,
    _to_numeric_matrix,
    plot_pubmatrix_heatmap,
    pubmatrix_heatmap,
)


# ---------------------------------------------------------------------------
# _to_numeric_matrix
# ---------------------------------------------------------------------------

class TestToNumericMatrix:
    def test_dataframe_input(self, sample_matrix):
        arr = _to_numeric_matrix(sample_matrix)
        assert isinstance(arr, np.ndarray)
        assert arr.shape == (3, 3)

    def test_ndarray_input(self):
        arr = _to_numeric_matrix(np.array([[1, 2], [3, 4]]))
        assert arr.dtype == float

    def test_list_of_lists(self):
        arr = _to_numeric_matrix([[1, 2], [3, 4]])
        assert arr.shape == (2, 2)

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            _to_numeric_matrix(np.array([[]]))

    def test_1d_raises(self):
        with pytest.raises(ValueError):
            _to_numeric_matrix(np.array([1, 2, 3]))


# ---------------------------------------------------------------------------
# _handle_na
# ---------------------------------------------------------------------------

class TestHandleNa:
    def test_no_nans_unchanged(self):
        arr = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = _handle_na(arr)
        np.testing.assert_array_equal(result, arr)

    def test_nans_replaced_with_zero(self):
        arr = np.array([[1.0, np.nan], [np.nan, 4.0]])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = _handle_na(arr)
        assert not np.isnan(result).any()
        assert result[0, 1] == 0.0
        assert len(w) == 1

    def test_all_zero_no_warning(self, all_zero_matrix):
        arr = _to_numeric_matrix(all_zero_matrix)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _handle_na(arr)
        assert len(w) == 0


# ---------------------------------------------------------------------------
# _overlap_percentage
# ---------------------------------------------------------------------------

class TestOverlapPercentage:
    def test_shape_preserved(self, sample_matrix):
        arr = _to_numeric_matrix(sample_matrix)
        pct = _overlap_percentage(arr)
        assert pct.shape == arr.shape

    def test_zero_matrix_returns_zeros(self, all_zero_matrix):
        arr = _to_numeric_matrix(all_zero_matrix)
        pct = _overlap_percentage(arr)
        assert (pct == 0.0).all()

    def test_values_between_0_and_100(self, sample_matrix):
        arr = _to_numeric_matrix(sample_matrix)
        pct = _overlap_percentage(arr)
        assert (pct >= 0).all() and (pct <= 100).all()

    def test_single_cell(self, single_cell_matrix):
        arr = _to_numeric_matrix(single_cell_matrix)
        pct = _overlap_percentage(arr)
        assert pct[0, 0] == 100.0  # 5 / (5+5-5) * 100 = 100


# ---------------------------------------------------------------------------
# _clustered_order
# ---------------------------------------------------------------------------

class TestClusteredOrder:
    def test_single_row_returns_identity(self):
        arr = np.array([[1.0, 2.0, 3.0]])
        assert _clustered_order(arr) == [0]

    def test_uniform_rows_returns_identity(self):
        arr = np.ones((4, 3))
        assert _clustered_order(arr) == [0, 1, 2, 3]

    def test_returns_permutation(self, sample_matrix):
        arr = _to_numeric_matrix(sample_matrix)
        order = _clustered_order(arr)
        assert sorted(order) == list(range(arr.shape[0]))

    def test_float_near_uniform_handled(self):
        # Rows differ only by floating-point noise — should not crash
        arr = np.ones((3, 3)) + np.random.default_rng(42).random((3, 3)) * 1e-15
        order = _clustered_order(arr)
        assert sorted(order) == [0, 1, 2]


# ---------------------------------------------------------------------------
# _auto_font_size
# ---------------------------------------------------------------------------

class TestAutoFontSize:
    def test_small_matrix(self):
        assert _auto_font_size(3, 4) == 10.0

    def test_medium_matrix(self):
        assert _auto_font_size(8, 10) == 8.0

    def test_large_matrix(self):
        assert _auto_font_size(15, 20) == 6.0

    def test_very_large_matrix(self):
        assert _auto_font_size(50, 50) == 4.0


# ---------------------------------------------------------------------------
# plot_pubmatrix_heatmap
# ---------------------------------------------------------------------------

class TestPlotPubmatrixHeatmap:
    def test_returns_fig_and_ax(self, sample_matrix):
        fig, ax = plot_pubmatrix_heatmap(sample_matrix)
        assert isinstance(fig, Figure)
        import matplotlib.axes
        assert isinstance(ax, matplotlib.axes.Axes)

    def test_saves_file(self, sample_matrix, tmp_path):
        out = str(tmp_path / "test.png")
        fig, ax = plot_pubmatrix_heatmap(sample_matrix, filename=out)
        assert (tmp_path / "test.png").exists()

    def test_no_show_by_default(self, sample_matrix):
        import matplotlib.pyplot as plt
        with pytest.MonkeyPatch().context() as mp:
            shown = []
            mp.setattr(plt, "show", lambda: shown.append(True))
            plot_pubmatrix_heatmap(sample_matrix)
        assert shown == []

    def test_show_true_calls_show(self, sample_matrix):
        import matplotlib.pyplot as plt
        with pytest.MonkeyPatch().context() as mp:
            shown = []
            mp.setattr(plt, "show", lambda: shown.append(True))
            plot_pubmatrix_heatmap(sample_matrix, show=True)
        assert shown == [True]

    def test_no_clustering(self, sample_matrix):
        fig, ax = plot_pubmatrix_heatmap(sample_matrix, cluster_rows=False, cluster_cols=False)
        assert fig is not None

    def test_single_cell(self, single_cell_matrix):
        fig, ax = plot_pubmatrix_heatmap(single_cell_matrix)
        assert fig is not None


# ---------------------------------------------------------------------------
# pubmatrix_heatmap (wrapper)
# ---------------------------------------------------------------------------

class TestPubmatrixHeatmap:
    def test_returns_fig_and_ax(self, sample_matrix):
        fig, ax = pubmatrix_heatmap(sample_matrix)
        assert isinstance(fig, Figure)
