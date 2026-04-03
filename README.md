# 🧬 Genomic Variant Landscape & Population Comparison Dashboard
> DS3294 — Data Science Practice Project Team #6

---

## 📌 Overview

This project analyses genomic variant patterns using real human genetic variation data from the **gnomAD (Genome Aggregation Database)**. The goal is to explore how variant density and allele frequencies vary across genomic regions and populations — and to demonstrate principled handling of large, high-dimensional biological data.

In simple terms: we take a massive file of human DNA differences, clean it up, and produce charts and an interactive dashboard that answers questions like:
- *Where on the Y chromosome are variants most dense?*
- *How do allele frequencies differ across world populations?*
- *Does the choice of quality filter change our conclusions?*

---

## 📂 Data Source

**Dataset:** [gnomAD v4.1](https://gnomad.broadinstitute.org/)

| Field | Detail |
|---|---|
| File | `gnomad.genomes.v4.1.sites.chrY.vcf.bgz` |
| Chromosome | Y chromosome |
| Format | Compressed VCF (`.vcf.bgz`) |
| Access | Free public download from gnomad.broadinstitute.org |

**What the file contains:**
- Genomic coordinates (chromosome + position)
- Reference and alternate alleles (what changed in the DNA)
- Variant filter status (`PASS` = high quality)
- Allele counts and global allele frequencies
- Population-specific allele frequency breakdowns (AFR, AMR, EAS, EUR, SAS, etc.)

---

## 🗂️ Repository Structure

```
genomic-variant-dashboard/
│
├── data/
│   ├── raw/            ← Original downloaded VCF file goes here (do not modify)
│   ├── processed/      ← Cleaned, filtered Parquet files created by our pipeline
│   └── sample/         ← Small test subset (first 5000 rows) for fast testing
│
├── notebooks/
│   ├── 01_data_ingestion.ipynb       ← Load VCF → save as Parquet
│   ├── 02_cleaning_filtering.ipynb   ← Remove bad rows, apply quality thresholds
│   ├── 03_aggregation.ipynb          ← Count variants per region and population
│   ├── 04_visualisation.ipynb        ← All charts and graphs
│   └── 05_scalability_analysis.ipynb ← Memory benchmarks, VCF vs Parquet comparison
│
├── src/
│   ├── ingest.py       ← Functions to read and parse the VCF file
│   ├── filter.py       ← Functions to apply quality and frequency filters
│   ├── aggregate.py    ← Functions to count and summarise variants by region
│   └── plot.py         ← Reusable chart-drawing functions
│
├── tests/
│   ├── test_ingest.py      ← Tests that loading works correctly
│   ├── test_filter.py      ← Tests that filtering removes the right rows
│   └── test_aggregate.py   ← Tests that aggregation counts are correct
│
├── dashboard/
│   └── app.py          ← Interactive Streamlit dashboard (final product)
│
├── requirements.txt    ← All Python libraries needed
└── README.md           ← This file
```

---

## ⚙️ Pre-Processing Plan

### Step 1 — Ingest (`notebooks/01_data_ingestion.ipynb`)
- Read `.vcf.bgz` using `cyvcf2` Python library
- Parse key fields: `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, `FILTER`
- Unpack the `INFO` column (it contains 50+ annotations as a single string) into separate columns
- Save result as `.parquet` file for fast reloading

### Step 2 — Clean & Filter (`notebooks/02_cleaning_filtering.ipynb`)
- Keep only `FILTER == PASS` variants (removes low-quality calls)
- Drop rows where allele frequency (`AF`) is missing
- Apply MAF threshold: keep variants with frequency > 0.001
- **Sensitivity check:** re-run with MAF > 0.01 and document how results change

### Step 3 — Aggregate (`notebooks/03_aggregation.ipynb`)
- Divide the Y chromosome into 100kb windows
- Count variants per window → produces density map
- Group by population label → compare AF distributions
- Calculate SNP vs INDEL ratio across windows

### Step 4 — Visualise (`notebooks/04_visualisation.ipynb`)
| Chart | What it shows |
|---|---|
| Density bar plot | Variant count per 100kb genomic window |
| Violin plot | AF distribution by population |
| Heatmap | Variant type breakdown across genomic regions |
| Scatter plot | AF in one population vs another |

### Step 5 — Scalability Analysis (`notebooks/05_scalability_analysis.ipynb`)
- Compare load time and memory usage: raw VCF vs Parquet
- Document peak RAM usage per operation
- Propose chunked-loading strategy for larger chromosomes

---

## 🔄 Data Representation Comparison

| Strategy | Format | Read Speed | Memory Use |
|---|---|---|---|
| Raw compressed VCF | `.vcf.bgz` | Slow (sequential parse) | High |
| Parsed tabular | `.parquet` | Fast (columnar) | ~10× lower |

**Conclusion:** After ingestion, we work exclusively with Parquet files for speed and reproducibility.

---

## 🧪 Testing & Debugging Plan

- **Unit tests** using `pytest` in the `tests/` folder
- Each test loads the sample data (not the full file) so tests run in seconds
- Tests check:
  - Schema: do columns have the right types?
  - Filtering: does output have fewer rows than input?
  - Aggregation: do window counts sum to total variant count?
- Notebooks include `.describe()` and `.isna().sum()` checks at each step
- Chromosome Y results will be cross-checked against published gnomAD statistics

---

## 📈 Scalability & Memory

- chrY file is manageable (~hundreds of MB) — suitable for a single machine
- INFO field unpacking is the most memory-intensive step → done in chunks
- Future extension to other chromosomes will use Dask for out-of-core processing

---

## 🚀 Current Status

| Task | Status |
|---|---|
| Data source identified (gnomAD chrY) | ✅ Done |
| Repository structure created | ✅ Done |
| Schema designed | ✅ Done |
| Data ingestion notebook | 🔄 In Progress |
| Cleaning & filtering | 📋 Planned |
| Aggregation | 📋 Planned |
| Visualisations | 📋 Planned |
| Streamlit dashboard | 📋 Planned |
| Tests | 📋 Planned |
| Scalability benchmarks | 📋 Planned |

---

## 📦 Setup & Run

```bash
# 1. Clone the repo
git clone https://github.com/<your-username>/genomic-variant-dashboard.git
cd genomic-variant-dashboard

# 2. Install dependencies
pip install -r requirements.txt

# 3. Place your gnomAD file in data/raw/

# 4. Run notebooks in order (01 → 02 → 03 → 04 → 05)
jupyter notebook

# 5. Launch the dashboard
streamlit run dashboard/app.py
```

---

## 🔑 Key Dependencies

```
cyvcf2        # Read VCF files
pandas        # Data manipulation
pyarrow       # Parquet read/write
matplotlib    # Plotting
seaborn       # Statistical plots
plotly        # Interactive charts
streamlit     # Dashboard app
pytest        # Testing
```
