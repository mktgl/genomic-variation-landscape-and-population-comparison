# tests/test_ingest.py
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import pandas as pd
from pathlib import Path
import tempfile
import os
from src.ingest import parse_vcf_to_dataframe

# Note: For proper VCF testing, we create a minimal valid VCF in a temp file

@pytest.fixture
def sample_vcf_file():
    """Create a minimal VCF file for testing"""
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AF_afr,Number=A,Type=Float,Description="AF in African">
##INFO=<ID=AF_nfe,Number=A,Type=Float,Description="AF in European (non-Finnish)">
##INFO=<ID=AF_eas,Number=A,Type=Float,Description="AF in East Asian">
##INFO=<ID=AF_amr,Number=A,Type=Float,Description="AF in American">
##INFO=<ID=AF_sas,Number=A,Type=Float,Description="AF in South Asian">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele Number">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100000	.	A	T	.	PASS	AF=0.05;AF_afr=0.1;AF_nfe=0.02;AF_eas=0.03;AF_amr=0.04;AF_sas=0.06;AC=5;AN=100
chr1	200000	.	AT	ATT	.	LowQual	AF=0.12;AF_afr=0.15;AF_nfe=0.1;AF_eas=0.08;AF_amr=0.14;AF_sas=0.11;AC=12;AN=100
chr1	300000	.	G	C	.	PASS	AF=0.003;AF_afr=0.01;AF_nfe=0;AF_eas=0.005;AF_amr=.;AF_sas=0.002;AC=3;AN=1000
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        temp_path = f.name
    
    yield temp_path
    
    # Cleanup
    if os.path.exists(temp_path):
        os.unlink(temp_path)


def test_parse_vcf_to_dataframe(sample_vcf_file):
    df = parse_vcf_to_dataframe(sample_vcf_file)
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3
    assert list(df.columns) == [
        'chrom', 'pos', 'ref', 'alt', 'qual', 'filter_status',
        'af_global', 'af_afr', 'af_eur', 'af_eas', 'af_amr', 'af_sas',
        'ac', 'an', 'variant_type'
    ]
    
    # Check data types
    assert df['pos'].dtype == 'int64'
    assert df['af_global'].dtype == 'float64'
    
    # Check variant type classification
    assert df.iloc[0]['variant_type'] == 'SNP'   # A->T
    assert df.iloc[1]['variant_type'] == 'INDEL' # AT->ATT


def test_parse_vcf_to_dataframe_max_rows(sample_vcf_file):
    df = parse_vcf_to_dataframe(sample_vcf_file, max_rows=2)
    assert len(df) == 2


def test_parse_vcf_to_dataframe_variant_type_logic(sample_vcf_file):
    df = parse_vcf_to_dataframe(sample_vcf_file)
    
    assert df.iloc[0]['variant_type'] == 'SNP'
    assert df.iloc[1]['variant_type'] == 'INDEL'
