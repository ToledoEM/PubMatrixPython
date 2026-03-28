"""
PubMatrix heatmap visualisation — mirrors heatmap_functions.R from PubMatrixR.

Provides overlap-percentage heatmaps with optional hierarchical clustering.
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


# Default red gradient matching R pheatmap palette
_RED_GRADIENT = ["#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#99000d"]


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
    Compute Jaccard-style overlap percentage for each cell.

    overlap[i, j] = intersection / union * 100
    where union = row_total[i] + col_total[j] - intersection
    """
    row_totals = arr.sum(axis=1, keepdims=True)   # sum across columns per row
    col_totals = arr.sum(axis=0, keepdims=True)   # sum across rows per column
    union = row_totals + col_totals - arr
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(union > 0, arr / union * 100, 0.0)
    return pct


def _clustered_order(arr: np.ndarray) -> list[int]:
    """Return row indices reordered by Euclidean distance / average linkage."""
    if arr.shape[0] < 2:
        return list(range(arr.shape[0]))
    if np.all(arr == arr[0]):  # no variation — skip clustering
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
    title: str = "PubMatrix Co-occurrence Heatmap",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    show_numbers: bool = True,
    color_palette: list[str] | None = None,
    filename: str | None = None,
    width: float = 10,
    height: float = 8,
    scale_font: bool = True,
) -> plt.Axes:
    """
    Create a publication-ready heatmap of PubMatrix co-occurrence results.

    Cell values show overlap percentage: (intersection / union) × 100,
    where union = row_total + col_total - intersection.

    Parameters
    ----------
    matrix : DataFrame or array-like
        PubMatrix result (rows = B terms, cols = A terms, values = counts).
    title : str
        Heatmap title.
    cluster_rows, cluster_cols : bool
        Apply Euclidean distance / average-linkage clustering.
    show_numbers : bool
        Annotate cells with overlap percentage values.
    color_palette : list of str, optional
        Custom hex color list for gradient. Defaults to red gradient.
    filename : str, optional
        Save to this path (PNG). If None, display inline.
    width, height : float
        Figure size in inches.
    scale_font : bool
        Auto-scale annotation font size based on matrix dimensions.

    Returns
    -------
    matplotlib.axes.Axes
    """
    # --- Input handling ---
    row_labels = list(matrix.index) if isinstance(matrix, pd.DataFrame) else None
    col_labels = list(matrix.columns) if isinstance(matrix, pd.DataFrame) else None

    arr = _to_numeric_matrix(matrix)
    arr = _handle_na(arr)
    pct = _overlap_percentage(arr)

    n_rows, n_cols = arr.shape

    # --- Clustering ---
    row_order = _clustered_order(pct) if cluster_rows else list(range(n_rows))
    col_order = _clustered_order(pct.T) if cluster_cols else list(range(n_cols))

    pct_ordered = pct[np.ix_(row_order, col_order)]
    row_labels_ordered = [row_labels[i] for i in row_order] if row_labels else row_order
    col_labels_ordered = [col_labels[i] for i in col_order] if col_labels else col_order

    # --- Color map ---
    colors = color_palette or _RED_GRADIENT
    cmap = LinearSegmentedColormap.from_list("pubmatrix", colors)

    # --- Font size ---
    annot_kws = {}
    if scale_font:
        annot_kws["size"] = _auto_font_size(n_rows, n_cols)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        pct_ordered,
        ax=ax,
        cmap=cmap,
        annot=show_numbers,
        fmt=".1f",
        annot_kws=annot_kws or None,
        xticklabels=col_labels_ordered,
        yticklabels=row_labels_ordered,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Overlap %"},
    )

    ax.set_title(title, pad=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()

    if filename:
        fig.savefig(filename, dpi=150, bbox_inches="tight")
        print(f"Saved heatmap to {filename}")
    else:
        plt.show()

    return ax


def pubmatrix_heatmap(matrix, title: str = "PubMatrix Results") -> plt.Axes:
    """
    Convenience wrapper for plot_pubmatrix_heatmap() with default parameters.

    Parameters
    ----------
    matrix : DataFrame or array-like
        PubMatrix result matrix.
    title : str
        Heatmap title.

    Returns
    -------
    matplotlib.axes.Axes
    """
    return plot_pubmatrix_heatmap(matrix, title=title)
