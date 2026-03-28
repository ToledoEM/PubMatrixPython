# Changelog

## [0.1.0] — 2026-03-28

Initial release. Python port of [PubMatrixR](https://github.com/ToledoEM/PubMatrixR-v2).

- `pubmatrix()` — pairwise PubMed/PMC co-occurrence queries with progress bar
- `pubmatrix_from_file()` — load term lists from a plain-text file
- `plot_pubmatrix_heatmap()` — heatmap with optional clustering, custom colours, PNG export
- `pubmatrix_heatmap()` — quick wrapper with defaults
- Date range filtering via `daterange`
- CSV and ODS export with PubMed hyperlinks
- NCBI API key support for higher rate limits
