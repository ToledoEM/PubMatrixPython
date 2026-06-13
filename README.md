# PubMatrixPython v0.2

<img src="https://toledoem.github.io/img/LogoPubmatrixP.png" align="right" width="150"/>

![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Tests](https://img.shields.io/badge/tests-60%20passed-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

Python port of the [PubMatrixR](https://github.com/ToledoEM/PubMatrixR-v2) R package.

For every pair of search terms `(A, B)`, it counts how many PubMed or PMC publications mention both. Good for mapping relationships between genes, diseases, and pathways across the literature.

Based on: Becker et al. (2003) *PubMatrix: a tool for multiplex literature mining*. BMC Bioinformatics 4:61. https://doi.org/10.1186/1471-2105-4-61

---

## Key features

- **Pairwise literature search** — automatically searches every combination of terms from two lists
- **PubMed or PMC** — query MEDLINE abstracts or PMC full text via NCBI E-utilities
- **Heatmap visualisation** — overlap-percentage heatmaps with optional hierarchical clustering
- **Export to CSV or ODS** — results include clickable hyperlinks to the matching PubMed search
- **Date filtering** — restrict searches to a publication year range
- **Flexible input** — pass term lists directly, or load them from a text file
- **Concurrency** — `n_workers` for parallel queries, respecting NCBI rate limits
- **Disk caching** — `cache_dir` persists query results between runs
- **Progress tracking** — built-in progress bar for long searches

## Use cases

- **Gene–disease association studies** — explore literature connections between genes and diseases
- **Pathway analysis** — investigate co-occurrence of genes within or across biological pathways
- **Drug–target research** — analyse relationships between compounds and potential targets
- **Systematic literature reviews** — quantify research coverage across multiple topics
- **Knowledge gap identification** — find under-researched combinations of terms
- **Bibliometric analysis** — measure research activity in a domain over time

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
    export_format="csv",      # saves results_result.csv with PubMed hyperlinks
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
    A,                    # list of str — column terms
    B,                    # list of str — row terms
    api_key=None,         # NCBI API key (10 req/s vs 3 req/s default)
    database="pubmed",    # "pubmed" or "pmc"
    daterange=None,       # e.g. [2015, 2024]
    outfile=None,         # base filename for export
    export_format=None,   # None | "csv" | "ods"
    n_tries=2,            # retries on network failure
    n_workers=1,          # parallel workers for concurrent queries
    timeout=30,           # HTTP request timeout in seconds
    cache_dir=None,       # directory to cache query results on disk
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

Heatmap of overlap percentages with optional hierarchical clustering. Returns `(fig, ax)`.

```python
fig, ax = plot_pubmatrix_heatmap(
    matrix,                                        # DataFrame from pubmatrix()
    title="PubMatrix Co-occurrence Heatmap",
    cluster_rows=True,
    cluster_cols=True,
    show_numbers=True,
    color_palette=None,                            # list of hex colours
    filename=None,                                 # save to PNG if set
    width=10, height=8,
    scale_font=True,
    show=False,                                    # call plt.show() after plotting
)
```

### `pubmatrix_heatmap(matrix, title=...)`

Quick wrapper around `plot_pubmatrix_heatmap()` with all defaults. Returns `(fig, ax)`.

---

## Output files

When `outfile` and `export_format` are set, results are written to
`{outfile}_result.{extension}` (`.csv` or `.ods`). Each cell contains the
publication count and a hyperlink to the matching PubMed search. Row names
come from `B`, column names from `A`.

ODS export requires the optional `odfpy` dependency:

```bash
pip install pubmatrixpython[ods]
```

---

## NCBI API key

Without a key: 3 requests/second. With a key: 10 requests/second.
Get one at https://account.ncbi.nlm.nih.gov/

```python
result = pubmatrix(A=A, B=B, api_key="YOUR_KEY_HERE")
```

---

## More documentation

- [Performance notes](docs/performance.md) — rate limits, caching, concurrency
- [Troubleshooting](docs/troubleshooting.md) — empty results, rate limiting, slow searches
- [Full reference notebook](https://toledoem.github.io/pubmatrixp/) — every parameter and feature, with output

---

## License & citation

This project is licensed under the MIT License — see [`LICENSE.md`](LICENSE.md).

If you use PubMatrixPython in your research, please cite:

> Becker KG, Hosack DA, Dennis G Jr, Lempicki RA, Bright TJ, Cheadle C, Engel J.
> *PubMatrix: a tool for multiplex literature mining.*
> BMC Bioinformatics. 2003 Dec 10;4:61. https://doi.org/10.1186/1471-2105-4-61

**Developers:**
- Tyler Laird (Author, original PubMatrixR)
- Enrique Toledo (Author, maintainer)
