# PubMatrixPython v0.3

<img src="https://toledoem.github.io/img/LogoPubmatrixP.png" align="right" width="150"/>

![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Tests](https://img.shields.io/badge/tests-83%20passed-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)
[![PyPI](https://img.shields.io/pypi/v/pubmatrixpython)](https://pypi.org/project/pubmatrixpython/)
[![codecov](https://codecov.io/gh/ToledoEM/PubMatrixPython/graph/badge.svg?token=FI57AHLE8V)](https://codecov.io/gh/ToledoEM/PubMatrixPython)

Python port of the [PubMatrixR](https://github.com/ToledoEM/PubMatrixR-v2) R package.

Give it two lists of search terms. For every pair, it asks PubMed or PMC how many publications mention both, and hands back the counts as a table. Useful when you want to know which gene/disease combinations the literature has actually covered, and which nobody has looked at.

Based on: Becker et al. (2003) *PubMatrix: a tool for multiplex literature mining*. BMC Bioinformatics 4:61. https://doi.org/10.1186/1471-2105-4-61

---

## What it does

Queries every combination of terms from your two lists against MEDLINE abstracts (`pubmed`) or full text (`pmc`) through the NCBI E-utilities. Results come back as a `pandas.DataFrame` and can be exported to CSV or ODS, where each cell links to the PubMed search that produced it.

You can restrict searches to a publication year range, pass terms directly or read them from a text file, and plot the result as a heatmap with optional clustering. Long runs show a progress bar.

Two things worth knowing for larger matrices: `n_workers` runs queries in parallel while staying under the NCBI rate limit, and `cache_dir` keeps results on disk so a re-run does not re-query terms it has already seen.

---

## Installation

### pip

```bash
pip install pubmatrixpython
```

### uv

```bash
uv add pubmatrixpython
```

### pixi

```bash
pixi add --pypi pubmatrixpython
```

For ODS export you also need `odfpy`:

```bash
pip install pubmatrixpython[ods]
```

---

## Development setup

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

Run every `uv` command from the project root, where `pyproject.toml` lives. It
will not resolve the environment from anywhere else.

```bash
cd /path/to/PubMatrixPython
uv run jupyter lab
```

The notebooks are in `notebooks/`.

| Notebook | What it covers |
|----------|---------------|
| `01_pubmatrix.ipynb` | Basic queries, date filtering, PMC database, file input, CSV export, heatmap visualisation |
| `02_example_wnt.ipynb` | Full worked example: WNT genes × obesity genes |

---

## Quick start

### In the REPL

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

Blank lines are ignored. If the `#` is missing, or either side of it is empty,
you get an error naming the file.

---

## API reference

### `pubmatrix(A, B, ...)`

Queries PubMed and returns a `pandas.DataFrame` with rows from `B` and columns from `A`.

```python
pubmatrix(
    A,                    # list of str — column terms
    B,                    # list of str — row terms
    api_key=None,         # NCBI API key (10 req/s vs 3 req/s default)
    database="pubmed",    # "pubmed" or "pmc"
    daterange=None,       # e.g. [2015, 2024]
    outfile=None,         # base filename for export
    export_format=None,   # None | "csv" | "ods"
    n_tries=3,            # attempts per query, exponential backoff between them
    n_workers=1,          # parallel workers for concurrent queries
    timeout=30,           # HTTP request timeout in seconds
    cache_dir=None,       # directory to cache query results on disk
)
```

### `pubmatrix_from_file(filepath, ...)`

Reads terms from a plain-text file, then hands them to `pubmatrix()` along with
any other arguments you pass.

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

Plots the co-occurrence counts, clustering rows and columns unless you turn that
off. Returns `(fig, ax)`.

```python
fig, ax = plot_pubmatrix_heatmap(
    matrix,                                        # DataFrame from pubmatrix()
    values="raw",                                  # "raw" | "row_pct" | "relative"
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

`values` selects what each cell shows:

| Value | Cell contents |
|-------|---------------|
| `"raw"` (default) | The co-occurrence counts themselves |
| `"row_pct"` | Each count as a percentage of its row total |
| `"relative"` | `count / (row_total + col_total - count) × 100` |

A warning about `"relative"`: it is not a Jaccard index, and its numbers do not
carry between runs. The totals in that formula are sums over whichever partner
terms happen to be in your matrix, not the publication count for each term on
its own, so adding one unrelated term shifts the value in every existing cell. A
real Jaccard index would need the single-term counts, and `pubmatrix()` never
fetches those. Compare cells within one fixed matrix and nothing beyond that.

### `pubmatrix_heatmap(matrix, title=..., values="raw")`

Quick wrapper around `plot_pubmatrix_heatmap()` with all defaults. Returns `(fig, ax)`.

---

## Output files

Set `outfile` and `export_format` and the results are written to `{outfile}.csv`
or `{outfile}.ods`. Every cell holds the publication count as a hyperlink to the
search that produced it. Rows are named from `B`, columns from `A`.

ODS export needs the optional `odfpy` dependency. See [Installation](#installation).

---

## NCBI API key

Without a key: 3 requests/second. With a key: 10 requests/second.
Get one at https://account.ncbi.nlm.nih.gov/

```python
result = pubmatrix(A=A, B=B, api_key="YOUR_KEY_HERE")
```

---

## More documentation

- [Performance notes](docs/performance.md) covers rate limits, retries, caching and concurrency
- [Troubleshooting](docs/troubleshooting.md) covers empty results, rate limiting and slow searches
- [Full reference notebook](https://toledoem.github.io/pubmatrixp/) walks through every parameter with output

---

## License & citation

MIT licensed. See [`LICENSE.md`](LICENSE.md).

If you use PubMatrixPython in your research, please cite:

> Becker KG, Hosack DA, Dennis G Jr, Lempicki RA, Bright TJ, Cheadle C, Engel J.
> *PubMatrix: a tool for multiplex literature mining.*
> BMC Bioinformatics. 2003 Dec 10;4:61. https://doi.org/10.1186/1471-2105-4-61

**Developers:**
- Tyler Laird (Author, original PubMatrixR)
- Enrique Toledo (Author, maintainer)
