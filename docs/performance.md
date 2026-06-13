# Performance notes

## Rate limiting

NCBI allows 3 requests/second without an API key, and 10/second with one.
Pass `api_key="YOUR_KEY"` to `pubmatrix()` to use the higher limit. Get a key
at https://account.ncbi.nlm.nih.gov/.

## Search time

Total search time depends on:

- **Matrix size** — `len(A) * len(B)` queries are made
- **`n_workers`** — parallel queries, bounded by the NCBI rate limit
- **Network speed** — each query is a round trip to `eutils.ncbi.nlm.nih.gov`

## Caching

Set `cache_dir` to write each query result to disk as JSON. Repeated calls
with the same terms and parameters load from disk instead of re-querying
NCBI:

```python
result = pubmatrix(A=A, B=B, cache_dir="cache", api_key=API_KEY)
```

## Memory usage

Results are held in memory as a `pandas.DataFrame`. Very large matrices
(hundreds of terms on each axis) may require substantial RAM.
