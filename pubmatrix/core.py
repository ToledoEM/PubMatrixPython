"""
PubMatrix core — systematic literature co-occurrence analysis via NCBI E-utilities.

Mirrors the R PubMatrixR package (https://github.com/ToledoEM/PubMatrixR-v2).
Reference: Becker et al. (2003) BMC Bioinformatics 4:61. doi:10.1186/1471-2105-4-61
"""

import math
import time
import urllib.parse
import xml.etree.ElementTree as ET
from itertools import product
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_SEARCH_BASE = "https://www.ncbi.nlm.nih.gov/{db}/?term={term}"

VALID_DATABASES = {"pubmed", "pmc"}
VALID_EXPORT_FORMATS = {None, "csv", "ods"}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_count(xml_text: str) -> int:
    """Parse publication count from NCBI esearch XML response."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        raise ValueError(f"Could not parse NCBI XML response: {e}") from e

    count_el = root.find(".//Count")
    if count_el is None:
        raise ValueError("NCBI XML response missing <Count> element")

    text = (count_el.text or "").strip()
    if not text.isdigit():
        raise ValueError(f"<Count> value is not numeric: {text!r}")

    return int(text)


def _fetch_count(base_url: str, encoded_term: str, n_tries: int = 2) -> int:
    """Fetch publication count for a single search term with retry logic."""
    url = f"{base_url}&term={encoded_term}&usehistory=y"
    last_error = None

    for attempt in range(n_tries):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return _extract_count(response.text)
        except Exception as e:
            last_error = e
            if attempt < n_tries - 1:
                time.sleep(0.25 * (attempt + 1))

    raise RuntimeError(
        f"Failed to fetch count after {n_tries} attempts for term {encoded_term!r}: {last_error}"
    )


def _validate_daterange(daterange):
    """Validate and normalise daterange parameter. Returns (start, end) tuple or None."""
    if daterange is None:
        return None

    if len(daterange) != 2:
        raise ValueError("daterange must have exactly 2 elements: [start_year, end_year]")

    start, end = daterange
    if not (math.isfinite(start) and math.isfinite(end)):
        raise ValueError("daterange values must be finite numbers")

    start, end = int(round(start)), int(round(end))
    if start > end:
        raise ValueError(f"daterange start ({start}) must be <= end ({end})")

    return (start, end)


def _build_base_url(database: str, api_key: str | None, daterange) -> str:
    """Construct the NCBI esearch base URL with optional API key and date range."""
    params = [f"db={database}", "rettype=count", "retmode=xml"]

    if api_key:
        params.append(f"api_key={api_key}")

    if daterange is not None:
        start, end = daterange
        params.append(f"mindate={start}&maxdate={end}&datetype=pdat")

    return f"{NCBI_BASE}?{'&'.join(params)}"


def _build_hyperlink_url(database: str, term: str) -> str:
    """Build a PubMed/PMC search URL for a given term."""
    encoded = urllib.parse.quote(term)
    return PUBMED_SEARCH_BASE.format(db=database, term=encoded)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def pubmatrix(
    A: list[str],
    B: list[str],
    api_key: str | None = None,
    database: str = "pubmed",
    daterange=None,
    outfile: str | None = None,
    export_format: str | None = None,
    n_tries: int = 2,
) -> pd.DataFrame:
    """
    Query PubMed/PMC and build a pairwise co-occurrence matrix.

    For each pair (a, b) in A × B, counts publications matching 'a AND b'.

    Parameters
    ----------
    A : list of str
        Search terms for matrix columns.
    B : list of str
        Search terms for matrix rows.
    api_key : str, optional
        NCBI API key (allows 10 req/s instead of 3 req/s).
    database : str
        'pubmed' (default) or 'pmc'.
    daterange : list or tuple of 2 ints, optional
        [start_year, end_year] to filter by publication date.
    outfile : str, optional
        Base filename for export (required if export_format is set).
    export_format : str, optional
        None (no export), 'csv', or 'ods'.
    n_tries : int
        Number of retry attempts for failed requests (default 2).

    Returns
    -------
    pandas.DataFrame
        Rows = B terms, columns = A terms, values = publication counts.
    """
    # --- Validation (fail-fast, same order as R package) ---
    if export_format not in VALID_EXPORT_FORMATS:
        raise ValueError(f"export_format must be one of {VALID_EXPORT_FORMATS}, got {export_format!r}")

    if export_format is not None and outfile is None:
        raise ValueError("outfile must be specified when export_format is set")

    if database not in VALID_DATABASES:
        raise ValueError(f"database must be one of {VALID_DATABASES}, got {database!r}")

    daterange = _validate_daterange(daterange)

    if not A or not B:
        raise ValueError("A and B must be non-empty lists")

    A = [str(t).strip() for t in A]
    B = [str(t).strip() for t in B]

    if any(not t for t in A):
        raise ValueError("A contains empty or whitespace-only terms")
    if any(not t for t in B):
        raise ValueError("B contains empty or whitespace-only terms")

    # --- Build queries ---
    pairs = list(product(B, A))  # rows × cols, matches R expand.grid(B, A)
    encoded_terms = [
        urllib.parse.quote(f"{b} AND {a}") for b, a in pairs
    ]

    base_url = _build_base_url(database, api_key, daterange)

    # --- Fetch counts ---
    counts = []
    for encoded in tqdm(encoded_terms, desc="Querying NCBI", unit="query"):
        counts.append(_fetch_count(base_url, encoded, n_tries=n_tries))

    if len(counts) != len(B) * len(A):
        raise RuntimeError(
            f"Expected {len(B) * len(A)} counts, got {len(counts)}"
        )

    # --- Assemble matrix (rows=B, cols=A) ---
    data = {}
    for j, a in enumerate(A):
        data[a] = [counts[i * len(A) + j] for i in range(len(B))]

    df = pd.DataFrame(data, index=B)
    df.index.name = None

    # --- Optional export ---
    if export_format == "csv":
        _export_csv(df, outfile, database)
    elif export_format == "ods":
        _export_ods(df, outfile, database)

    return df


def pubmatrix_from_file(filepath: str, **kwargs) -> pd.DataFrame:
    """
    Load search terms from a file and run pubmatrix().

    File format:
        term_A1
        term_A2
        #
        term_B1
        term_B2

    Parameters
    ----------
    filepath : str
        Path to a plain-text file with A terms, a '#' separator, then B terms.
    **kwargs
        Passed directly to pubmatrix().

    Returns
    -------
    pandas.DataFrame
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    lines = [ln.strip() for ln in path.read_text().splitlines()]
    lines = [ln for ln in lines if ln]  # drop blank lines

    if "#" not in lines:
        raise ValueError("File must contain a '#' line separating A and B terms")

    sep = lines.index("#")
    A = lines[:sep]
    B = lines[sep + 1:]

    if not A or not B:
        raise ValueError("File must contain terms both before and after the '#' separator")

    return pubmatrix(A=A, B=B, **kwargs)


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def _make_hyperlink_formula(url: str, value: int) -> str:
    """Excel-compatible HYPERLINK formula."""
    return f'=HYPERLINK("{url}","{value}")'


def _export_csv(df: pd.DataFrame, outfile: str, database: str) -> None:
    """Export matrix to CSV with Excel HYPERLINK formulas."""
    path = Path(outfile).with_suffix(".csv")
    rows = []

    for b_term in df.index:
        row = {}
        for a_term in df.columns:
            term = f"{a_term} AND {b_term}"
            url = _build_hyperlink_url(database, term)
            count = df.loc[b_term, a_term]
            row[a_term] = _make_hyperlink_formula(url, count)
        rows.append(row)

    export_df = pd.DataFrame(rows, index=df.index)
    export_df.to_csv(path)
    print(f"Saved CSV to {path}")


def _export_ods(df: pd.DataFrame, outfile: str, database: str) -> None:
    """Export matrix to ODS with hyperlinks."""
    from odf.opendocument import OpenDocumentSpreadsheet
    from odf.style import Style, TextProperties
    from odf.table import Table, TableRow, TableCell
    from odf.text import A as OdfA, P

    path = Path(outfile).with_suffix(".ods")
    doc = OpenDocumentSpreadsheet()

    link_style = Style(name="LinkStyle", family="text")
    link_style.addElement(TextProperties(color="#0000EE", textunderlinestyle="solid"))
    doc.styles.addElement(link_style)

    table = Table(name="PubMatrix")

    # Header row
    header_row = TableRow()
    header_row.addElement(TableCell(valuetype="string"))  # empty corner
    for a_term in df.columns:
        cell = TableCell(valuetype="string")
        cell.addElement(P(text=a_term))
        header_row.addElement(cell)
    table.addElement(header_row)

    # Data rows
    for b_term in df.index:
        row = TableRow()
        label_cell = TableCell(valuetype="string")
        label_cell.addElement(P(text=b_term))
        row.addElement(label_cell)

        for a_term in df.columns:
            count = int(df.loc[b_term, a_term])
            term = f"{a_term} AND {b_term}"
            url = _build_hyperlink_url(database, term)

            cell = TableCell(valuetype="string")
            p = P()
            link = OdfA(href=url, text=str(count))
            p.addElement(link)
            cell.addElement(p)
            row.addElement(cell)

        table.addElement(row)

    doc.spreadsheet.addElement(table)
    doc.save(str(path))
    print(f"Saved ODS to {path}")
