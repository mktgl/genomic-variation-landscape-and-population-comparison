# tests/test_filter.py
# ─────────────────────────────────────────────────────────────
# PURPOSE: Verify that our filtering functions work correctly
# Run with: pytest tests/
# ─────────────────────────────────────────────────────────────

import pandas as pd
import sys
sys.path.append('../src')
from filter import filter_pass_only, drop_missing_af, apply_maf_threshold


def make_test_df():
    """Creates a small fake variant table for testing."""
    return pd.DataFrame({
        'pos'          : [100, 200, 300, 400, 500],
        'filter_status': ['PASS', 'FAIL', 'PASS', 'PASS', 'FAIL'],
        'af_global'    : [0.05, 0.10, None, 0.001, 0.50],
        'variant_type' : ['SNP', 'SNP', 'INDEL', 'SNP', 'INDEL']
    })


def test_pass_filter_reduces_rows():
    df = make_test_df()
    result = filter_pass_only(df)
    # Only rows with PASS should remain (rows 0, 2, 3)
    assert len(result) == 3, f'Expected 3 rows, got {len(result)}'
    assert (result['filter_status'] == 'PASS').all()


def test_pass_filter_does_not_modify_original():
    df = make_test_df()
    _ = filter_pass_only(df)
    assert len(df) == 5  # original should be untouched


def test_drop_missing_af():
    df = make_test_df()
    result = drop_missing_af(df)
    assert result['af_global'].isna().sum() == 0
    assert len(result) == 4  # row 2 (None) should be removed


def test_maf_threshold():
    df = make_test_df()
    result = apply_maf_threshold(df, min_maf=0.01)
    # Only rows with af_global >= 0.01 should remain
    assert (result['af_global'] >= 0.01).all()
    assert len(result) < len(df)


def test_maf_threshold_strict_less_than_lenient():
    df = make_test_df()
    lenient = apply_maf_threshold(df, min_maf=0.001)
    strict  = apply_maf_threshold(df, min_maf=0.01)
    assert len(strict) <= len(lenient), 'Strict filter should keep same or fewer rows'
