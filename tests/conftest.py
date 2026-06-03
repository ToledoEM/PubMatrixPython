"""Shared fixtures for PubMatrix tests."""

import numpy as np
import pandas as pd
import pytest


VALID_XML = """<?xml version="1.0" encoding="UTF-8"?>
<eSearchResult>
  <Count>42</Count>
</eSearchResult>"""

ZERO_COUNT_XML = """<?xml version="1.0" encoding="UTF-8"?>
<eSearchResult>
  <Count>0</Count>
</eSearchResult>"""

MALFORMED_XML = "this is not xml <<<"

MISSING_COUNT_XML = """<?xml version="1.0" encoding="UTF-8"?>
<eSearchResult>
  <RetMax>0</RetMax>
</eSearchResult>"""

NON_NUMERIC_COUNT_XML = """<?xml version="1.0" encoding="UTF-8"?>
<eSearchResult>
  <Count>not_a_number</Count>
</eSearchResult>"""


@pytest.fixture
def sample_matrix():
    """3×3 DataFrame mimicking pubmatrix() output."""
    return pd.DataFrame(
        [[10, 5, 2], [3, 8, 1], [0, 4, 6]],
        index=["B1", "B2", "B3"],
        columns=["A1", "A2", "A3"],
    )


@pytest.fixture
def all_zero_matrix():
    return pd.DataFrame(
        np.zeros((3, 3), dtype=int),
        index=["B1", "B2", "B3"],
        columns=["A1", "A2", "A3"],
    )


@pytest.fixture
def single_cell_matrix():
    return pd.DataFrame([[5]], index=["B1"], columns=["A1"])
