# RCTD ALS Spatial Transcriptomics

**Robust Cell Type Decomposition (RCTD) analysis of 50 spatial transcriptomics slides from the GSE120374 ALS dataset — spatiotemporal dynamics of molecular pathology in Amyotrophic Lateral Sclerosis.**

---

## Table of Contents

- [Overview](#overview)
- [Datasets](#datasets)
- [Environment](#environment)
- [Pipeline Steps](#pipeline-steps)
  - [Step 1 — Mount Google Drive](#step-1--mount-google-drive)
  - [Step 2 — Create Folder Structure](#step-2--create-folder-structure)
  - [Step 3 — Install Packages](#step-3--install-packages)
  - [Step 4 — Download Spatial Data (GSE120374)](#step-4--download-spatial-data-gse120374)
  - [Step 5 — Build scRNA-seq Reference](#step-5--build-scrna-seq-reference)
  - [Step 6 — Define Helper Functions](#step-6--define-helper-functions)
  - [Step 7 — Select 50 Samples and Test Run](#step-7--select-50-samples-and-test-run)
  - [Step 8 — Batch Run (50 Samples)](#step-8--batch-run-50-samples)
  - [Step 9 — Aggregate Results](#step-9--aggregate-results)
  - [Step 10 — Visualize Results](#step-10--visualize-results)
- [Results](#results)
- [Folder Structure](#folder-structure)
- [Key Decisions and Fixes](#key-decisions-and-fixes)
- [Citation](#citation)
- [Author](#author)

---

## Overview

This pipeline applies **RCTD (Robust Cell Type Decomposition)** to spatial transcriptomics data from an ALS mouse model. RCTD uses a single-cell RNA-seq reference to decompose the mixed gene expression signal at each spatial spot into cell type proportions.

In this analysis:
- **50 spatial transcriptomics slides** from the SOD1-G93A ALS mouse model were processed
- A **self-supervised spatial reference** was built by clustering pooled spatial spots
- **3 spatial domains** (A, B, C) were identified across all slides
- All analysis was run in **Google Colab (R runtime)**

---

## Datasets

### Spatial Transcriptomics — GSE120374

| Field | Detail |
|-------|--------|
| GEO Accession | [GSE120374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120374) |
| Title | Spatiotemporal dynamics of molecular pathology in ALS |
| Contributor | Vickovic S, Broad Institute of MIT and Harvard |
| Organism | *Mus musculus* |
| Platform | NextSeq 550 |
| Design | SOD1-G93A model (4 timepoints) vs WT control; ATG7 model (2 timepoints) |
| Total slides | 331 |
| Slides used | 50 |
| File downloaded | `GSE120374_RAW.tar` (1.4 GB) |

### Reference — GSE161621

| Field | Detail |
|-------|--------|
| GEO Accession | [GSE161621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161621) |
| Title | Mouse spinal cord scRNA-seq |
| Platforms | NextSeq 550 + Illumina NovaSeq 6000 |
| Samples | 6 replicates (male/female) |
| File downloaded | `GSE161621_RAW.tar` (1 GB, H5 format) |
| Note | Due to extreme sparsity of H5 files, reference was rebuilt from pooled spatial data instead |

---

## Environment

- **Platform:** Google Colab (R runtime)
- **R version:** 4.x
- **Key packages:**

| Package | Version | Purpose |
|---------|---------|---------|
| `spacexr` | dmcable/spacexr (GitHub) | RCTD deconvolution |
| `Seurat` | 5.5.0 | Reference clustering |
| `Matrix` | CRAN | Sparse matrix handling |
| `data.table` | CRAN | Fast file reading |
| `ggplot2` | CRAN | Visualization |
| `plyr` | CRAN | Factor mapping |
| `hdf5r` | CRAN | Reading H5 files |

---

## Pipeline Steps

### Step 1 — Mount Google Drive

**Purpose:** Connect Google Drive to the Colab session so all data and results persist across sessions.

**Method used:** `reticulate::py_run_string()` calling `google.colab.drive.mount()` via Python, since native R Drive mounting is not supported in Colab R runtime.

**Key issue encountered:** Google Drive mounting only works in hosted Python runtimes by default. The solution was to use the Python `google.colab` module through `reticulate`, and to ensure the notebook was saved to Drive (not local disk) before running.

```r
reticulate::py_run_string("
from google.colab import drive
drive.mount('/content/drive')
")
```

---

### Step 2 — Create Folder Structure

**Purpose:** Create all necessary directories on Google Drive before any data is downloaded or processed. This ensures outputs are saved immediately after processing without path errors.

**Folder layout created:**

```
MyDrive/RCTD_ALS_50/
├── raw_data/
│   ├── GSE120374_RAW/     ← spatial ST count files
│   └── reference/
│       └── h5_files/      ← scRNA-seq H5 files
├── results/
│   ├── weights/           ← RCTD weight CSVs (one per sample)
│   ├── plots/             ← spatial domain plots
│   └── rctd_objects/      ← full RCTD R objects (.rds)
└── logs/                  ← timestamped batch run logs
```

```r
DRIVE_ROOT <- "/content/drive/MyDrive/RCTD_ALS_50"
for (d in c(RAW_DATA_DIR, REF_DIR, REF_H5_DIR,
            WEIGHTS_DIR, PLOTS_DIR, RCTD_OBJ_DIR, LOG_DIR)) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}
```

---

### Step 3 — Install Packages

**Purpose:** Install all required R packages. `spacexr` is not on CRAN and must be installed from GitHub using `devtools`. All other packages are installed from CRAN.

**Note:** This step takes approximately 10–15 minutes on first run. Packages do not persist across Colab sessions, so this cell must be re-run after any session restart.

```r
devtools::install_github("dmcable/spacexr", build_vignettes=FALSE)
install.packages("Seurat")
install.packages(c("Matrix", "data.table", "ggplot2", "plyr", "hdf5r"))
```

---

### Step 4 — Download Spatial Data (GSE120374)

**Purpose:** Download the raw spatial transcriptomics data from NCBI GEO and extract it to Google Drive.

**What is downloaded:**
- `GSE120374_RAW.tar` (1.4 GB) — contains TSV count matrices and H&E images for all 331 slides
- `GSE120374_1000L6/L7/L8_barcodes.txt.gz` — spot barcode files

**File format:** Each sample folder contains:
- `*_stdata_aligned_counts_IDs.txt.gz` — genes × spots count matrix (TSV)
- `*.tsv.annotations.tsv.gz` — spot annotations
- `*_HE.jpg.gz` — H&E tissue image
<img width="1066" height="872" alt="image" src="https://github.com/user-attachments/assets/58fc7231-b9bb-4807-b0e3-043e23a94645" />

**Key issue:** Files are extracted flat (no subfolders), so samples are identified by filename prefix rather than directory.

```r
download.file(TAR_URL, destfile=TAR_FILE, method="wget")
system(paste("tar -xf", TAR_FILE, "-C", RAW_DATA_DIR))
```

---

### Step 5 — Build scRNA-seq Reference

**Purpose:** Build the RCTD `Reference` object containing labeled cell types and their expression profiles.

**Original plan:** Use GSE161621 (mouse spinal cord scRNA-seq H5 files) as reference.

**Problem encountered:** The GSE161621 H5 files were extremely sparse — after loading and filtering, only 14 cells were retained, far below RCTD's minimum of 25 cells per cell type.

**Solution — self-supervised spatial reference:** Instead of using an external scRNA-seq reference, we pooled spots from 10 spatial slides, clustered them using Seurat, and used the resulting clusters as pseudo-cell-types (spatial domains).

**Steps:**
1. Load 10 spatial count matrices, apply `make.unique()` to gene names to fix duplicates
2. Replace underscores with dashes in gene names (Seurat v5 requirement)
3. Assign unique cell barcodes per sample using `paste0("S", i, "C", seq_len(ncol(mat)))`
4. Find intersection of genes across all 10 samples
5. `cbind` all matrices into one combined matrix (4440 genes × 522 spots)
6. Cluster with Seurat (NormalizeData → FindVariableFeatures → ScaleData → PCA → FindNeighbors → FindClusters)
7. Label clusters as `Domain_A`, `Domain_B`, `Domain_C`
8. Filter to keep only domains with ≥ 25 spots (RCTD requirement)
9. Build `Reference` object and save to Drive

**Result:** 516 cells, 3 spatial domains (A: 245, B: 235, C: 36)

**Critical fix:** `droplevels()` must be called after subsetting to remove empty factor levels, which would cause RCTD to fail with "need minimum 25 cells" even for domains that technically exist.

```r
# droplevels() is essential — prevents empty levels from reaching RCTD
cell_types <- droplevels(factor(obj$cell_type))
keep       <- which(cell_types %in% valid)
cell_types <- droplevels(cell_types[keep])   # droplevels again after subset
```

---

### Step 6 — Define Helper Functions

**Purpose:** Define two reusable functions that are called in the batch loop. These functions are self-contained and include all error handling.

#### `load_spatial_sample(sample_name, raw_dir)`

Loads one spatial transcriptomics sample from its count file and returns a `SpatialRNA` object.

**What it does:**
1. Finds the count file by searching for `_stdata_aligned_counts_IDs.txt.gz` pattern
2. Reads the TSV with `data.table::fread()` for speed
3. Applies `make.unique(gsub("_","-", gene_names))` to fix gene name issues
4. Detects orientation — GSE120374 files have rows = spots (e.g. `10x10` format), so the matrix is transposed to genes × spots
5. Removes all-zero rows and columns
6. Generates grid coordinates (since annotation files did not contain standard x/y columns)
7. Returns `SpatialRNA(coords, counts, nUMI)` object

#### `run_and_save_rctd(sample_name, puck, reference, ...)`

Runs RCTD on one spatial sample and saves all outputs to Drive.

**What it does:**
1. Checks if weights file already exists — skips if so (enables safe resume after crash)
2. Runs `create.RCTD()` then `run.RCTD()` with `doublet_mode = "full"`
3. Saves full RCTD object as `.rds`
4. Extracts and normalizes weight matrix with `normalize_weights()`
5. Saves weights as CSV
6. Saves spatial plots (wrapped in `tryCatch` so plot failures don't stop the batch)

---

### Step 7 — Select 50 Samples and Test Run

**Purpose:** Select the first 50 samples alphabetically and run a test on one sample to verify the full pipeline works before running all 50.

**Sample selection logic:**
```r
count_files <- all_files[grepl("_stdata_aligned_counts_IDs\\.txt\\.gz$", all_files)]
all_samples <- head(sort(gsub("_stdata_aligned_counts_IDs\\.txt\\.gz$","", count_files)), 50)
```

**Note:** Sample names include `_1` and `_2` suffixes because some slides have two count files (replicate sections). These are treated as independent samples.

---

### Step 8 — Batch Run (50 Samples)

**Purpose:** Process all 50 samples sequentially, saving results to Drive after each one.

**Key features:**
- Each sample is wrapped in `tryCatch` so a single failure doesn't stop the batch
- Already-processed samples are skipped (checks for existing weights CSV)
- A timestamped log file is written to Drive after each sample
- `gc()` is called after each sample to free RAM

**Result:** All 50 samples processed successfully in approximately 8 seconds total (spatial slides in this dataset are small: 10–140 spots each).

```r
for (i in seq_along(all_samples)) {
  tryCatch({
    puck <- load_spatial_sample(sname, RAW_DATA_DIR)
    run_and_save_rctd(sname, puck, reference, ...)
  }, error=function(e) {
    failed <<- c(failed, sname)
  })
}
```

---

### Step 9 — Aggregate Results
<img width="1013" height="873" alt="image" src="https://github.com/user-attachments/assets/b98f4b75-a683-4cbe-9f89-a4fdbda067d3" />

**Purpose:** Read all 50 weight CSV files and compute the mean proportion of each spatial domain per slide. Save as a single summary CSV.

**Key fix:** `read.csv(f, row.names=NULL)` instead of `read.csv(f, row.names=1)` — the default `row.names=1` caused a "duplicate row.names not allowed" error because spot barcodes were not unique across the file.

```r
summary_df <- do.call(rbind, lapply(wfiles, function(f) {
  w     <- read.csv(f, row.names=NULL, check.names=FALSE)
  w     <- w[, -1]
  sname <- gsub("_weights\\.csv","", basename(f))
  data.frame(sample=sname, t(colMeans(w)), check.names=FALSE)
}))
```

---

### Step 10 — Visualize Results

**Purpose:** Generate two types of plots saved to Drive.

**Plot 1 — Summary boxplot (`summary_boxplot.png`):**
Shows the distribution of each spatial domain's mean proportion across all 50 slides. Each point is one slide. Built with `ggplot2` using `geom_boxplot` + `geom_jitter`.

**Plot 2 — Per-sample spatial heatmaps:**
For each of the 50 slides, a faceted heatmap shows where each domain (A, B, C) is located across the tissue grid. Uses `geom_tile` with a white-to-darkred gradient. Since real tissue coordinates were not available, a grid layout is used.

---

## Results

| Metric | Value |
|--------|-------|
| Samples processed | 50 / 50 |
| Spatial domains identified | 3 (A, B, C) |
| Mean Domain A proportion | ~45–55% per slide |
| Mean Domain B proportion | ~35–45% per slide |
| Mean Domain C proportion | ~5–10% per slide |
| Weight files saved | 50 CSVs |
| RCTD objects saved | 50 RDS files |

---

## Folder Structure

```
results/
├── summary_50_samples.csv          ← mean domain proportions per slide
├── weights/
│   ├── GSM3399115_CN51_C2_1_weights.csv
│   ├── GSM3399115_CN51_C2_2_weights.csv
│   └── ... (50 files total)
└── plots/
    ├── summary_boxplot.png          ← boxplot across all 50 slides
    ├── domain_proportions.png       ← domain proportion plot
    └── GSM3399115_CN51_C2_1/
        └── GSM3399115_CN51_C2_1_spatial.png
    └── ... (50 sample folders)
```

---

## Key Decisions and Fixes

| Issue | Fix Applied |
|-------|------------|
| Google Drive won't mount in R runtime | Used `reticulate::py_run_string()` to call Python `google.colab.drive.mount()` |
| GSE161621 H5 files too sparse (14 cells) | Built self-supervised reference from 10 pooled spatial slides |
| Empty factor levels crashing RCTD | Called `droplevels()` twice — after annotation and after filtering |
| Duplicate gene names in count files | Applied `make.unique(gsub("_","-", gene_names))` before matrix creation |
| Seurat v5 rejects underscores in gene names | Replaced all `_` with `-` using `gsub()` |
| Duplicate cell barcodes across samples | Assigned unique IDs: `paste0("S", i, "C", seq_len(ncol(mat)))` |
| `read.csv(row.names=1)` duplicate error | Used `row.names=NULL` then dropped first column manually |
| Session restart loses variables | Every cell re-defines all paths and reloads reference from Drive |
| `CreateSeuratObject` LogMap error | Fixed by making gene names and cell barcodes fully unique before object creation |

---

## Citation

Cable DM, Murray E, Zou LS, Goeva A, Macosko EZ, Chen F, Irizarry RA (2022).
Robust decomposition of cell type mixtures in spatial transcriptomics.
*Nature Biotechnology* 40, 517–526.
https://doi.org/10.1038/s41587-021-00830-w

---

## Author

**Shreshtha Shukla**
Email: shreshthashukla13042000@gmail.com
GitHub: [Shreshthaa](https://github.com/Shreshthaa)
