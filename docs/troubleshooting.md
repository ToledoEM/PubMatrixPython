# Troubleshooting

## Empty results

If many searches return 0 results, try:

- Using broader search terms
- Expanding the date range
- Checking spelling of scientific terms
- Using alternative gene names or synonyms

## Rate limiting errors (HTTP 429)

If `pubmatrix()` raises a `RuntimeError` after repeated 429 responses:

- Obtain and use an NCBI API key (`api_key=...`)
- Reduce `n_workers` or the size of your search matrix
- Increase `n_tries` to allow more retries with backoff

## Long search times

For large matrices:

- Break the search into smaller sub-matrices
- Use a more specific `daterange`
- Set `cache_dir` so repeated runs skip already-fetched pairs
- Increase `n_workers`, within the NCBI rate limit
