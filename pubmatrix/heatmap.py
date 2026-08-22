"""
PubMatrix heatmap visualisation — mirrors heatmap_functions.R from PubMatrixR.

Plots co-occurrence counts, with optional percentage views and hierarchical clustering.
"""

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

logger = logging.getLogger(__name__)

# Default red gradient matching R pheatmap palette
_RED_GRADIENT = ["#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#99000d"]

VALID_VALUES = ("raw", "row_pct", "relative")

_CBAR_LABELS = {"raw": "Count", "row_pct": "Row %", "relative": "Overlap %"}


def _to_numeric_matrix(matrix) -> np.ndarray:
    """Coerce input to a 2-D numeric numpy array."""
    if isinstance(matrix, pd.DataFrame):
        arr = matrix.values.astype(float)
    elif isinstance(matrix, np.ndarray):
        arr = matrix.astype(float)
    else:
        arr = np.array(matrix, dtype=float)

    if arr.ndim != 2 or arr.shape[0] == 0 or arr.shape[1] == 0:
        raise ValueError("matrix must be a non-empty 2-D array or DataFrame")

    return arr


def _handle_na(arr: np.ndarray) -> np.ndarray:
    """Replace NaN with 0, emitting a warning if any were found."""
    nan_mask = np.isnan(arr)
    if nan_mask.any():
        positions = list(zip(*np.where(nan_mask)))
        warnings.warn(
            f"NA values found at positions {positions[:5]}{'...' if len(positions) > 5 else ''}. "
            "Converting to 0.",
            UserWarning,
            stacklevel=3,
        )
        arr = arr.copy()
        arr[nan_mask] = 0.0
    return arr


def _overlap_percentage(arr: np.ndarray) -> np.ndarray:
    """
    Compute the "relative" overlap percentage for each cell.

    overlap[i, j] = count / (row_total[i] + col_total[j] - count) * 100

    This is **not** a Jaccard index and is **not comparable across runs**. The
    totals are sums over whichever partner terms happen to be in the supplied
    matrix, not the marginal publication counts for each term on its own, so
    adding an unrelated term changes the value reported for every existing cell.
    A true Jaccard index would need the single-term counts, which pubmatrix()
    does not fetch. Use it only to compare cells within one fixed matrix.
    """
    row_totals = arr.sum(axis=1, keepdims=True)   # sum across columns per row
    col_totals = arr.sum(axis=0, keepdims=True)   # sum across rows per column
    union = row_totals + col_totals - arr
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(union > 0, arr / union * 100, 0.0)
    return pct


def _row_percentage(arr: np.ndarray) -> np.ndarray:
    """Express each cell as a percentage of its row total (empty rows stay 0)."""
    row_totals = arr.sum(axis=1, keepdims=True)
    denom = np.where(row_totals == 0, 1.0, row_totals)
    return np.round(arr / denom * 100, 1)


def _display_values(arr: np.ndarray, values: str) -> np.ndarray:
    """Build the matrix that is actually plotted, per the `values` selector."""
    if values == "raw":
        return arr
    if values == "row_pct":
        return _row_percentage(arr)
    return np.round(_overlap_percentage(arr), 1)


def _clustered_order(arr: np.ndarray) -> list[int]:
    """Return row indices reordered by Euclidean distance / average linkage."""
    if arr.shape[0] < 2:
        return list(range(arr.shape[0]))
    if np.allclose(arr, arr[0]):  # no variation — skip clustering
        return list(range(arr.shape[0]))
    dist = pdist(arr, metric="euclidean")
    Z = linkage(dist, method="average")
    dend = dendrogram(Z, no_plot=True)
    return dend["leaves"]


def _auto_font_size(n_rows: int, n_cols: int) -> float:
    """Scale annotation font size based on matrix dimensions."""
    max_dim = max(n_rows, n_cols)
    if max_dim <= 5:
        return 10.0
    elif max_dim <= 10:
        return 8.0
    elif max_dim <= 20:
        return 6.0
    else:
        return 4.0


def plot_pubmatrix_heatmap(
    matrix,
    values: str = "raw",
    title: str = "PubMatrix Co-occurrence Heatmap",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    show_numbers: bool = True,
    color_palette: list[str] | None = None,
    filename: str | None = None,
    width: float = 10,
    height: float = 8,
    scale_font: bool = True,
    show: bool = False,
) -> tuple[Figure, plt.Axes]:
    """
    Create a publication-ready heatmap of PubMatrix co-occurrence results.

    Cells show the raw co-occurrence counts by default. See `values` for the
    percentage views.

    Parameters
    ----------
    matrix : DataFrame or array-like
        PubMatrix result (rows = B terms, cols = A terms, values = counts).
    values : str
        What each cell shows:

        - ``"raw"`` (default) — the co-occurrence counts themselves.
        - ``"row_pct"`` — each count as a percentage of its row total, useful
          for comparing how one row's attention is distributed across columns.
        - ``"relative"`` — ``count / (row_total + col_total - count) * 100``.
          Not a Jaccard index and not comparable across runs; its totals are
          sums over the supplied matrix, so adding an unrelated term changes
          the value reported for every existing cell.
    title : str
        Heatmap title.
    cluster_rows, cluster_cols : bool
        Apply Euclidean distance / average-linkage clustering.
    show_numbers : bool
        Annotate cells with the plotted values.
    color_palette : list of str, optional
        Custom hex color list for gradient. Defaults to red gradient.
    filename : str, optional
        Save to this path (PNG). If None, the figure is not saved automatically.
    width, height : float
        Figure size in inches.
    scale_font : bool
        Auto-scale annotation font size based on matrix dimensions.
    show : bool
        Call plt.show() after plotting (default False). Useful in interactive sessions.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
    """
    # --- Input handling ---
    if values not in VALID_VALUES:
        raise ValueError(f"values must be one of {VALID_VALUES}, got {values!r}")

    row_labels = list(matrix.index) if isinstance(matrix, pd.DataFrame) else None
    col_labels = list(matrix.columns) if isinstance(matrix, pd.DataFrame) else None

    arr = _to_numeric_matrix(matrix)
    arr = _handle_na(arr)
    display = _display_values(arr, values)

    n_rows, n_cols = arr.shape

    # --- Clustering ---
    row_order = _clustered_order(display) if cluster_rows else list(range(n_rows))
    col_order = _clustered_order(display.T) if cluster_cols else list(range(n_cols))

    display_ordered = display[np.ix_(row_order, col_order)]
    row_labels_ordered = [row_labels[i] for i in row_order] if row_labels else row_order
    col_labels_ordered = [col_labels[i] for i in col_order] if col_labels else col_order

    # --- Color map ---
    colors = color_palette or _RED_GRADIENT
    cmap = LinearSegmentedColormap.from_list("pubmatrix", colors)

    # A matrix where every value is the same is valid input, but the colour
    # scale degenerates — bracket it explicitly so one flat colour is drawn.
    low, high = float(display_ordered.min()), float(display_ordered.max())
    if low == high:
        warnings.warn(
            f"All plotted values are identical ({low}). "
            "Clustering has no effect and the heatmap will be a single colour.",
            UserWarning,
            stacklevel=2,
        )
        vmin, vmax = low - 0.5, high + 0.5
    else:
        vmin, vmax = low, high

    # --- Font size ---
    annot_kws = {}
    if scale_font:
        annot_kws["size"] = _auto_font_size(n_rows, n_cols)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        display_ordered,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        annot=show_numbers,
        # Counts are integers; percentages read better with one decimal.
        fmt=".0f" if values == "raw" else ".1f",
        annot_kws=annot_kws or None,
        xticklabels=col_labels_ordered,
        yticklabels=row_labels_ordered,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": _CBAR_LABELS[values]},
    )

    ax.set_title(title, pad=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()

    if filename:
        fig.savefig(filename, dpi=150, bbox_inches="tight")
        logger.info("Saved heatmap to %s", filename)

    if show:
        plt.show()

    return fig, ax


def pubmatrix_heatmap(
    matrix,
    title: str = "PubMatrix Results",
    values: str = "raw",
) -> tuple[Figure, plt.Axes]:
    """
    Convenience wrapper for plot_pubmatrix_heatmap() with default parameters.

    Parameters
    ----------
    matrix : DataFrame or array-like
        PubMatrix result matrix.
    title : str
        Heatmap title.
    values : str
        Passed through to plot_pubmatrix_heatmap(). Defaults to ``"raw"``
        (co-occurrence counts).

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
    """
    return plot_pubmatrix_heatmap(matrix, values=values, title=title)
