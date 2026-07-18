# Triticale Marker & GWAS Ideogram Pipeline

Physical ideograms for triticale (wheat A/B/B genome + rye R genome), built
from a linkage-map marker table and GAPIT GWAS results. Each marker/SNP is
placed on its actual chromosome using BLAST-derived physical coordinates,
rather than genetic map order, so gaps, clustering, and non-recombining
regions are directly visible.

## What this pipeline does

1. **`make_marker_ideogram.py`** - draws every marker from the linkage map
   onto its physical position on the wheat and rye chromosomes. This is the
   background/reference ideogram.
2. **`make_gwas_snp_ideogram.py`** - draws the significant SNPs from a single
   GAPIT MLM run on top of that same background, colored by significance.
3. **`run_all_gwas_ideograms.py`** - batch-runs script 2 over every trait/run
   in a GAPIT results folder (`rest_1`, `rest_2`, ...), producing one ideogram
   per run.
4. **`make_combined_gwas_ideogram.py`** - draws the significant SNPs from
   *all* runs on one ideogram, color-coded by run, with SNPs that are
   significant in 2+ runs highlighted separately as likely real/replicated
   hits.

All four scripts write both an HTML file (viewable in any browser) and a
standalone SVG file (for figures/publications).

## Requirements

```
pip install pandas openpyxl
```

Python 3.8+. No other dependencies - the ideograms are drawn as raw SVG, no
plotting library required.

## Input data

### 1. Marker/BLAST location table (`location_with_alleles_blast_filter_out.xlsx`)

One row per marker, with (at least) these columns:

| column | meaning |
|---|---|
| `Nr` | marker ID - **this is what GAPIT's `SNP` column is matched against** |
| `Chromosom` | chromosome as called on the genetic/linkage map (e.g. `2B`) |
| `wheat_Chr` | wheat BLAST hit accession, wrapped like `ref\|NC_057794.1\|` |
| `wheat_Pos` | marker's physical position on that wheat chromosome |
| `rye_Chr` | rye BLAST hit accession, wrapped like `emb\|OZ284383.1\|` |
| `rye_Pos` | marker's physical position on that rye chromosome |

`wheat_Pos`/`rye_Pos` are expected to already be the midpoint of the BLAST
hit (so they're safe to use even when the hit is on the minus strand and
`Start > End`).

Markers with no BLAST hit on a named pseudochromosome (e.g. only hitting an
unplaced `NW_*` scaffold) can't be placed on the ideogram and are dropped;
every script prints how many rows were dropped so nothing disappears
silently.

### 2. Chromosome name & length lookup table (`chr_NCBI_name_v1.xlsx`)

One row per chromosome, with columns:

| column | meaning |
|---|---|
| `Chromosome` | display name used on the ideogram (e.g. `2B`, `4R`) |
| `GenBank` | plain NCBI accession, e.g. `NC_057794.1` (matches the accession embedded in `wheat_Chr`/`rye_Chr` above, once the `ref\|...\|` / `emb\|...\|` wrapper is stripped) |
| `plant` | `wheat` or `rye` |
| `end` | chromosome length, used to scale the ideogram |

Wheat D-genome chromosomes (`1D`-`7D`) are intentionally excluded from the
wheat panel in the GWAS scripts (`make_gwas_snp_ideogram.py`,
`make_combined_gwas_ideogram.py`) - only A and B genome chromosomes are
drawn. `make_marker_ideogram.py` draws all of A, B, and D.

### 3. GAPIT results directory

Expected layout (used by `run_all_gwas_ideograms.py` and
`make_combined_gwas_ideogram.py`):

```
GAPIT_results/
├── rest_1/
│   └── MLM/
│       └── Significant_SNPs_p0.01.csv
├── rest_2/
│   └── MLM/
│       └── Significant_SNPs_p0.01.csv
...
```

Each `Significant_SNPs_p0.01.csv` must have (at least) `SNP` (marker ID,
matches `Nr` in the location table), `Chr`, and `P.value` columns - this is
the standard GAPIT MLM output format. The file-matching pattern
(`Significant_SNPs_p0*01.csv`) tolerates both `p0.01` and `p0_01` naming.

## Usage

### Base marker ideogram (run once)

```bash
python make_marker_ideogram.py \
    location_with_alleles_blast_filter_out.xlsx \
    chr_NCBI_name_v1.xlsx \
    triticale_marker_ideogram.html \
    triticale_marker_ideogram.svg
```

### Single GWAS run

```bash
python make_gwas_snp_ideogram.py \
    GAPIT_results/rest_1/MLM/Significant_SNPs_p0.01.csv \
    location_with_alleles_blast_filter_out.xlsx \
    chr_NCBI_name_v1.xlsx \
    rest_1_ideogram.html rest_1_ideogram.svg
```

### All GWAS runs at once (batch)

```bash
python run_all_gwas_ideograms.py \
    GAPIT_results \
    location_with_alleles_blast_filter_out.xlsx \
    chr_NCBI_name_v1.xlsx \
    ideogram_out
```

Writes `ideogram_out/rest_N_ideogram.{html,svg}` for every `rest_N` folder
found. Missing or unmatched folders are reported and skipped, not treated as
errors.

### All GWAS runs combined into one ideogram

```bash
python make_combined_gwas_ideogram.py \
    GAPIT_results \
    location_with_alleles_blast_filter_out.xlsx \
    chr_NCBI_name_v1.xlsx \
    combined_ideogram.html combined_ideogram.svg
```

## Reading the plots

- **Grey ticks**: every marker from the linkage map, for context (marker
  density/coverage along the chromosome).
- **Colored ticks** (single-run ideograms): significant SNPs, colored by
  -log10(P) - darker red is more significant.
- **Colored ticks** (combined ideogram): one fixed color per GWAS
  run/trait (see legend).
- **Black ticks** (combined ideogram only): a SNP that is significant in
  2 or more runs/traits - these are the most likely to be real, replicated
  associations rather than run-specific noise.
- Hovering over any tick in the HTML output shows the SNP ID, P-value, and
  (for the combined ideogram) which run(s) it was significant in.

## Known limitations

- Only markers/SNPs with a physical BLAST position on a named
  pseudochromosome can be placed - anything landing only on an unplaced
  scaffold is dropped (and reported).
- Chromosome scaling within each panel uses the longest chromosome
  *actually drawn* in that panel (so excluding the D genome doesn't distort
  the A/B genome scale).
- `run_all_gwas_ideograms.py` calls `make_gwas_snp_ideogram.py` as a
  subprocess, so both files need to stay in the same directory.
