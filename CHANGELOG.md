# Changelog

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
