# 🧬 Genomic Variant Landscape & Population Comparison

> **DS3294 — Data Science Practice · Team \#6**

[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Open%20Access-brightgreen?style=flat-square)](https://gnomad.broadinstitute.org/)
[![Data](https://img.shields.io/badge/Data-gnomAD%20v4.1-blueviolet?style=flat-square)](https://gnomad.broadinstitute.org/)
[![Dashboard](https://img.shields.io/badge/Dashboard-Streamlit-FF4B4B?style=flat-square&logo=streamlit&logoColor=white)](https://streamlit.io/)

---

## 📌 Overview

This project analyses human Y-chromosome genomic variation using data from the **gnomAD v4.1 (Genome Aggregation Database)**. Starting from a 571 MB compressed VCF containing **1,169,063** raw variant records, a five-stage pipeline produces cleaned data, publication-quality figures, and an interactive dashboard.

**Core questions answered:**
- Where on the Y chromosome are high-quality variants most concentrated?
- How do allele frequencies differ across AFR, EUR, EAS, AMR, and SAS populations?
- Does the MAF threshold choice materially change conclusions?
- How much faster is Parquet than raw VCF for iterative analysis?

---

## 👥 Authors

| Name | Roll No. |
|------|----------|
| Dharavath Anil | 20211225 |
| Jatin Raghuwanshi | 20221126 |
| Mukta Gajanan Londhe | 20221163 |
| Sattwik Maji | 20221239 |

---

## 📂 Data Source

**Dataset:** [gnomAD v4.1](https://gnomad.broadinstitute.org/) · Chromosome Y

| Field | Detail |
|-------|--------|
| File | `gnomad.genomes.v4.1.sites.chrY.vcf.bgz` |
| Genome Build | GRCh38/hg38 |
| Format | Block-gzip VCF (`.vcf.bgz`) |
| Raw File Size | 571.4 MB |
| Raw Variants | 1,169,063 |
| Retained after filtering (MAF > 0.001) | **37,089** |
| Populations | AFR, AMR, EAS, EUR, SAS |
| Access | Free public download |

---

## 🗂️ Repository Structure

```
genomic-variant-dashboard/
│
├── data/
│   ├── raw/                          ← Original VCF + raw Parquet (do not modify)
│   ├── processed/                    ← variants_filtered.parquet (pipeline output)
│   └── sample/                       ← First 5,000 rows for fast testing
│
├── notebooks/
│   ├── 01_data_ingestion.ipynb       ← Parse VCF → save as Parquet (1,169,063 rows)
│   ├── 02_cleaning_filtering.ipynb   ← Quality filters → 37,089 variants retained
│   ├── 03_aggregation.ipynb          ← Window density, population summaries
│   ├── 04_visualisation.ipynb        ← Five publication-quality figures
│   └── 05_scalability_analysis.ipynb ← VCF vs Parquet benchmarks
│
├── src/
│   ├── ingest.py                     ← parse_vcf_to_dataframe()
│   ├── filter.py                     ← full_filter_pipeline()
│   ├── aggregate.py                  ← windowed counts, population summaries
│   └── plot.py                       ← five reusable charting functions
│
├── tests/
│   ├── test_ingest.py
│   ├── test_filter.py
│   └── test_aggregate.py
│
├── results/
│   └── figures/                      ← Five PNG outputs (see Visualisations below)
│
├── dashboard/
│   └── app.py                        ← Interactive Streamlit dashboard
│
├── requirements.txt
└── README.md
```

---

## ⚙️ Pipeline

The analysis runs in five sequential stages, each backed by a Jupyter notebook and a dedicated source module.

```
Raw VCF  →  [NB 01 Ingest]  →  [NB 02 Clean & Filter]  →  [NB 03 Aggregate]  →  [NB 04 Visualise]  →  Dashboard
                                                                              ↘
                                                                         [NB 05 Scalability]
```

### Stage 1 — Ingest (`01_data_ingestion.ipynb` · `src/ingest.py`)
- Stream `.vcf.bgz` using `cyvcf2` — no full decompression needed
- Extract: `CHROM`, `POS`, `REF`, `ALT`, `FILTER`, `AF`, population AFs, `AC`, `AN`, variant type
- Output: `data/raw/variants.parquet` — **1,169,063 rows × 16 columns**

### Stage 2 — Clean & Filter (`02_cleaning_filtering.ipynb` · `src/filter.py`)

Three sequential filters applied:

| Step | Filter | Variants remaining |
|------|--------|--------------------|
| Raw | — | 1,169,063 |
| 1 | `FILTER == PASS` only | ~110,000 |
| 2 | Drop rows where `AF` is `NaN` | ~100,000 |
| 3 | MAF > 0.001 | **37,089** |

> **96.8% removed.** The non-recombining MSY has inherently low diversity; most variants are singletons or fail quality thresholds.

A sensitivity check with MAF > 0.01 confirmed qualitatively consistent results.

### Stage 3 — Aggregate (`03_aggregation.ipynb` · `src/aggregate.py`)
- Partition chrY into **100 kb non-overlapping windows** → variant density map
- Group by population → median, IQR, and distributional shape of allele frequencies
- Partition into **1 Mb segments** → SNP vs INDEL counts per region (heatmap input)

### Stage 4 — Visualise (`04_visualisation.ipynb` · `src/plot.py`)

| Figure | File | Description |
|--------|------|-------------|
| Variant density | `variant_density.png` | Variant count per 100 kb genomic window |
| AF violin | `af_violin.png` | Allele frequency distribution by population |
| Type heatmap | `variant_heatmap.png` | SNP vs INDEL counts across 1 Mb regions |
| EUR vs AFR scatter | `af_scatter_EUR_AFR.png` | Per-variant AF comparison between two populations |
| Faceted histograms | `af_histograms_faceted.png` | AF distribution across all five populations |

### Stage 5 — Scalability (`05_scalability_analysis.ipynb`)

| Format | File Size | Load Time | RAM Used | Python-native? |
|--------|-----------|-----------|----------|----------------|
| Raw VCF `.vcf.bgz` | 571.4 MB | 4.88 s *(50k rows only)* | 18 MB | No (`cyvcf2`) |
| Parquet `.parquet` | **18.8 MB** | **0.75 s** *(full dataset)* | 351 MB | Yes (`pandas`) |

**30.4× smaller on disk · ~6.5× faster to load.** All downstream analysis uses Parquet exclusively.

---

## 📊 Visualisations

All figures are produced by `src/plot.py` and saved to `results/figures/`.

Key findings:
- **Density:** highly non-uniform along chrY; euchromatic arms show peaks, heterochromatic core is sparse
- **Population AFs:** significant differentiation across groups, consistent with Y-haplogroup structure
- **EUR vs AFR scatter:** most variants lie off-diagonal (population-private), few shared variants along the diagonal
- **SNP:INDEL ratio:** exceeds 10:1 across all regions; modest INDEL enrichment in ampliconic segments

---

## 🖥️ Interactive Dashboard

Built with Streamlit (`dashboard/app.py`):

- Sidebar controls: MAF threshold, window size, population selector
- All five chart types regenerated dynamically on parameter change
- Data table view with column filtering
- CSV download of the filtered dataset

```bash
streamlit run dashboard/app.py
```

---

## 🧪 Testing

Unit tests in `tests/` use `pytest` and run on the 5,000-row sample subset (sub-second execution).

Tests cover:
- **Schema integrity** — correct column names and dtypes
- **Filter monotonicity** — output row count < input row count at each step
- **Aggregation correctness** — window counts sum to total filtered variant count
- **Edge cases** — empty input, all-NaN AF column, single-row DataFrame

```bash
pytest tests/ -v
```

---

## 🚀 Setup & Run

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/genomic-variant-dashboard.git
cd genomic-variant-dashboard

# 2. Install dependencies
pip install -r requirements.txt

# 3. Place the gnomAD VCF in data/raw/
#    Download from: https://gnomad.broadinstitute.org/

# 4. Run notebooks in order
jupyter nbconvert --to notebook --execute notebooks/01_data_ingestion.ipynb
jupyter nbconvert --to notebook --execute notebooks/02_cleaning_filtering.ipynb
jupyter nbconvert --to notebook --execute notebooks/03_aggregation.ipynb
jupyter nbconvert --to notebook --execute notebooks/04_visualisation.ipynb
jupyter nbconvert --to notebook --execute notebooks/05_scalability_analysis.ipynb

# 5. Run the test suite
pytest tests/ -v

# 6. Launch the dashboard
streamlit run dashboard/app.py
```

---

## 📦 Key Dependencies

```
cyvcf2>=0.30.18     # Stream VCF files efficiently
pandas>=2.0.0       # Data manipulation
pyarrow>=14.0.0     # Parquet read/write
matplotlib>=3.7.0   # Base plotting
seaborn>=0.12.0     # Statistical plots
plotly>=5.18.0      # Interactive charts
streamlit>=1.30.0   # Dashboard
pytest>=7.4.0       # Unit testing
psutil>=5.9.0       # Memory profiling
```

Full list in `requirements.txt`.

---

## 📈 Scalability Notes

The chrY VCF (571 MB) is tractable on a single machine. Scaling to full-genome gnomAD data would require:

- **Chunked VCF parsing** — stream records in batches
- **Dask DataFrames** — out-of-core drop-in for pandas
- **Partitioned Parquet** — partition by genomic region for predicate pushdown

---

## 📄 License

Data sourced from gnomAD under its open-access terms. Code in this repository is released for academic use under the course guidelines of DS3294.
