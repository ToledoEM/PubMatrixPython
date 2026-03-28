# PubMatrixPython v0.1

<img src="https://toledoem.github.io/img/LogoPubmatrixP.png" align="right" width="150"/>

Python port of the [PubMatrixR](https://github.com/ToledoEM/PubMatrixR-v2) R package.

For every pair of search terms `(A, B)`, it counts how many PubMed or PMC publications mention both. Good for mapping relationships between genes, diseases, and pathways across the literature.

Based on: Becker et al. (2003) *PubMatrix: a tool for multiplex literature mining*. BMC Bioinformatics 4:61. https://doi.org/10.1186/1471-2105-4-61

---

## Setup

Requires [uv](https://docs.astral.sh/uv/). Install it with:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Clone and install dependencies:

```bash
git clone <repo-url>
cd PubMatrixPython
uv sync --all-groups
```

---

## Running the notebooks

All `uv` commands must be run from the **project root** (`PubMatrixPython/`), where `pyproject.toml` lives.

```bash
cd /path/to/PubMatrixPython
uv run jupyter lab
```

Then open any notebook from the `notebooks/` folder in the browser.

| Notebook | What it covers |
|----------|---------------|
| `01_pubmatrix.ipynb` | Basic queries, date filtering, PMC database, file input, CSV export, heatmap visualisation |
| `02_example_wnt.ipynb` | Full worked example: WNT genes × obesity genes |

---

## Quick start (script or REPL)

### Interactive REPL

```bash
uv run python
```

```python
from pubmatrix import pubmatrix, plot_pubmatrix_heatmap

A = ["WNT1", "WNT2", "CTNNB1"]
B = ["obesity", "diabetes", "cancer"]

result = pubmatrix(A=A, B=B)
print(result)

plot_pubmatrix_heatmap(result, title="WNT × Disease")
```

### Running a script

Create a file `my_analysis.py`:

```python
from pubmatrix import pubmatrix, plot_pubmatrix_heatmap

A = ["WNT1", "WNT2", "WNT3A", "WNT5A", "CTNNB1"]
B = ["obesity", "diabetes", "cancer", "inflammation"]

result = pubmatrix(
    A=A,
    B=B,
    database="pubmed",
    daterange=[2010, 2024],   # optional date filter
    outfile="results",
    export_format="csv",      # saves results.csv with PubMed hyperlinks
)

print(result)

plot_pubmatrix_heatmap(
    result,
    title="WNT Genes × Disease",
    filename="heatmap.png",   # saves to file instead of displaying
)
```

Run it with:

```bash
uv run python my_analysis.py
```

### Loading terms from a file

Create `terms.txt`:

```
WNT1
WNT2
CTNNB1
#
obesity
diabetes
cancer
```

```python
from pubmatrix import pubmatrix_from_file

result = pubmatrix_from_file("terms.txt")
print(result)
```

```bash
uv run python my_analysis.py
```

---

## API reference

### `pubmatrix(A, B, ...)`

Query PubMed and return a `pandas.DataFrame` (rows = B, cols = A).

```python
pubmatrix(
    A,                        # list of str — column terms
    B,                        # list of str — row terms
    api_key=None,             # NCBI API key (10 req/s vs 3 req/s default)
    database="pubmed",        # "pubmed" or "pmc"
    daterange=None,           # e.g. [2015, 2024]
    outfile=None,             # base filename for export
    export_format=None,       # None | "csv" | "ods"
    n_tries=2,                # retries on network failure
)
```

### `pubmatrix_from_file(filepath, ...)`

Load terms from a plain-text file and run `pubmatrix()`.

File format:
```
WNT1
WNT2
#
obesity
diabetes
```

```python
result = pubmatrix_from_file("terms.txt", database="pubmed")
```

### `plot_pubmatrix_heatmap(matrix, ...)`

Heatmap of overlap percentages with optional hierarchical clustering.

```python
plot_pubmatrix_heatmap(
    matrix,                                        # DataFrame from pubmatrix()
    title="PubMatrix Co-occurrence Heatmap",
    cluster_rows=True,
    cluster_cols=True,
    show_numbers=True,
    color_palette=None,                            # list of hex colours
    filename=None,                                 # save to PNG if set
    width=10, height=8,
    scale_font=True,
)
```

### `pubmatrix_heatmap(matrix, title=...)`

Quick wrapper around `plot_pubmatrix_heatmap()` with all defaults.

---

## NCBI API key

Without a key: 3 requests/second. With a key: 10 requests/second.
Get one at https://account.ncbi.nlm.nih.gov/

```python
result = pubmatrix(A=A, B=B, api_key="YOUR_KEY_HERE")
```
