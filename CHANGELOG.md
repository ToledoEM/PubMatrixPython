# Changelog

## [0.3.0] — 2026-08-22

Parity with [PubMatrixR 1.0.1](https://github.com/ToledoEM/PubMatrixR-v2/blob/main/NEWS.md).

### Breaking changes

- `plot_pubmatrix_heatmap()` and `pubmatrix_heatmap()` now plot raw co-occurrence counts by
  default, and take a new `values` argument to choose otherwise. Earlier versions always plotted
  a percentage the docs called `intersection / union * 100`, but the totals in that formula were
  row and column sums of whatever matrix you passed in, not the publication counts for each term
  on its own. Adding an unrelated term changed the number shown in every other cell. The old
  calculation is still there as `values="relative"`, now documented with that caveat.
  `values="row_pct"` divides each count by its row total.
- `pubmatrix()` rejects duplicate terms within `A` or within `B` instead of silently producing an
  ambiguous matrix.
- Queries drop `usehistory=y` and set `retmax=0`, so responses no longer carry PMIDs the package
  never looks at. **This changes the cache key** — entries written by earlier versions are ignored
  and will be re-fetched once.

### Bug fixes

- Query terms are now built as `"A AND B"`, matching the hyperlink written on export. Counts were
  unaffected (the query is symmetric), but the exported link was not the URL that had been queried.
- Term file parse errors name the offending file, and distinguish "no terms before `#`" from
  "no terms after `#`".
- A matrix where every plotted value is identical now warns and draws a single-colour heatmap with
  an explicit colour range.

### Other improvements

- Failed queries back off exponentially between retries, and the default `n_tries` is now 3.
- Heatmap cells are formatted as integers for `values="raw"` and to one decimal for the percentage
  views; the colourbar label follows the selected `values`.
- Clustering now runs on the matrix that is actually plotted, not always on the percentages.
- Removed an unreachable count-length check in `pubmatrix()`. It was ported from R, where the
  parallel map could return a short list; the Python version builds its result by index, so a
  missing count already raises before the check could fire.
- Test suite grew from 60 to 111 tests, covering the CSV and ODS export paths and the concurrent
  `n_workers` branch for the first time. Coverage is at 100%.

## [0.2.0] — 2026-06-03

- `n_workers` parameter for concurrent NCBI queries (respects rate limits automatically)
- `cache_dir` parameter to cache query results to disk and skip redundant requests
- `timeout` parameter exposed on `pubmatrix()` (was hardcoded at 30 s)
- `plot_pubmatrix_heatmap()` now returns `(fig, ax)` tuple instead of `ax` alone
- Added `show` parameter to `plot_pubmatrix_heatmap()`; `plt.show()` no longer called automatically
- Replaced `print()` with `logging` throughout — output can now be suppressed or redirected
- Narrowed exception handling in `_fetch_count()` to `requests.RequestException`
- Fixed float-unsafe clustering check (`np.allclose` instead of `==`)
- `n_tries` and `n_workers` now validated with clear error messages
- `odfpy` moved to optional extra: `pip install pubmatrixpython[ods]`
- Requires Python ≥ 3.10 (relaxed from 3.13)
- Full test suite: 60 tests covering core queries, XML parsing, caching, retry logic, and heatmap rendering

## [0.1.0] — 2026-03-28

Initial release. Python port of [PubMatrixR](https://github.com/ToledoEM/PubMatrixR-v2).

- `pubmatrix()` — pairwise PubMed/PMC co-occurrence queries with progress bar
- `pubmatrix_from_file()` — load term lists from a plain-text file
- `plot_pubmatrix_heatmap()` — heatmap with optional clustering, custom colours, PNG export
- `pubmatrix_heatmap()` — quick wrapper with defaults
- Date range filtering via `daterange`
- CSV and ODS export with PubMed hyperlinks
- NCBI API key support for higher rate limits
