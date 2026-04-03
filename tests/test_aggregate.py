# tests/test_aggregate.py
import pandas as pd
import sys
sys.path.append('../src')
from aggregate import compute_window_density, melt_population_af


def make_test_df():
    return pd.DataFrame({
        'pos'          : [50000, 150000, 160000, 250000, 350000],
        'variant_type' : ['SNP', 'SNP', 'INDEL', 'SNP', 'INDEL'],
        'af_global'    : [0.05, 0.10, 0.20, 0.30, 0.40],
        'af_afr'       : [0.06, 0.12, 0.18, None, 0.42],
        'af_eur'       : [0.04, 0.08, 0.22, 0.28, 0.38],
        'af_eas'       : [None, 0.09, 0.21, 0.31, 0.39],
        'af_amr'       : [0.05, 0.11, 0.19, 0.29, 0.41],
        'af_sas'       : [0.07, 0.10, 0.20, 0.32, 0.40],
    })


def test_window_counts_sum_to_total():
    df = make_test_df()
    windows = compute_window_density(df, window_size=100_000)
    # total variant count across all windows should equal len(df)
    assert windows['variant_count'].sum() == len(df)


def test_snp_indel_counts_match_variant_count():
    df = make_test_df()
    windows = compute_window_density(df, window_size=100_000)
    # snp_count + indel_count should equal variant_count for each window
    assert (windows['snp_count'] + windows['indel_count'] == windows['variant_count']).all()


def test_melt_population_af_drops_missing():
    df = make_test_df()
    result = melt_population_af(df)
    assert result['allele_frequency'].isna().sum() == 0


def test_melt_population_af_has_correct_populations():
    df = make_test_df()
    result = melt_population_af(df)
    expected_pops = {'African', 'European', 'East Asian', 'American', 'South Asian'}
    assert set(result['population'].unique()) == expected_pops
