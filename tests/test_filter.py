# tests/test_filter.py
import pytest
import pandas as pd
from src.filter import (
    filter_pass_only,
    drop_missing_af,
    apply_maf_threshold,
    full_filter_pipeline
)

@pytest.fixture
def raw_variant_df():
    """Sample raw DataFrame with various filtering conditions"""
    data = {
        'filter_status': ['PASS', 'LowQual', 'PASS', 'PASS', 'FAIL'],
        'af_global': [0.05, 0.002, None, 0.15, 0.008],
        'pos': [1000, 2000, 3000, 4000, 5000],
        'variant_type': ['SNP', 'SNP', 'INDEL', 'SNP', 'INDEL']
    }
    return pd.DataFrame(data)


def test_filter_pass_only(raw_variant_df):
    result = filter_pass_only(raw_variant_df)
    
    assert len(result) == 3
    assert all(result['filter_status'] == 'PASS')


def test_drop_missing_af(raw_variant_df):
    result = drop_missing_af(raw_variant_df)
    
    assert len(result) == 4
    assert result['af_global'].isna().sum() == 0


def test_apply_maf_threshold_default(raw_variant_df):
    result = apply_maf_threshold(raw_variant_df, min_maf=0.01)
    
    assert len(result) == 2  # Only 0.05 and 0.15 are >= 0.01
    assert all(result['af_global'] >= 0.01)


def test_apply_maf_threshold_custom_threshold(raw_variant_df):
    result = apply_maf_threshold(raw_variant_df, min_maf=0.1)
    
    assert len(result) == 1  # Only 0.15
    assert result.iloc[0]['af_global'] == 0.15


def test_full_filter_pipeline(raw_variant_df):
    result = full_filter_pipeline(raw_variant_df, min_maf=0.01)
    
    assert len(result) == 2
    assert all(result['filter_status'] == 'PASS')
    assert result['af_global'].isna().sum() == 0
    assert all(result['af_global'] >= 0.01)


def test_full_filter_pipeline_print_statement(capsys, raw_variant_df):
    """Test that the pipeline prints the filtering summary"""
    full_filter_pipeline(raw_variant_df, min_maf=0.01)
    captured = capsys.readouterr()
    assert "Filtering:" in captured.out
    assert "→" in captured.out
