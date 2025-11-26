# Triticale Restitution GWAS: Complete Pipeline & Workflow Documentation

## Overview

This repository provides a comprehensive workflow for conducting Genome-Wide Association Studies (GWAS) on restitution traits in triticale. It covers all steps from initial data preparation, GWAS analysis, postprocessing of significant marker results, functional annotation, and BLAST mapping of markers to reference genomes. The workflow leverages **Python** for data manipulation, sequence annotation, and BLAST, and **R** for GWAS and linkage disequilibrium analysis.

---

## 1. Data Preparation for GAPIT

**Script:** `preparing_for_second_GAPIT.py`

**Steps:**
- Reads marker data (`Tabela_S1.xlsx`) and phenotype data (`20250304_pheno_2n_gametes_ABDR.csv`).
- Harmonizes sample IDs, cleans and filters for missingness, and removes monomorphic SNPs:
  - Excludes markers with >20% missing genotypes.
  - Excludes samples with >20% missing data.
  - Removes monomorphic markers.
- Encodes genotypes: `a` → `0`, `b` → `2`, and imputes missing values.
- Prepares metadata for GWAS.

**Outputs:**
- `GAPIT_genotype.csv` (Taxa x SNPs)
- `GAPIT_phenotype.csv`
- `GAPIT_snp_metadata.csv` (SNP info: chromosome, Kosambi cM, etc.)

---

## 2. GWAS Analysis in R with GAPIT

**Script:** `GWAS_kosmbi_full_pipeline_all_models.R`

**Steps:**
- Loads cleaned genotype, phenotype, and SNP metadata.
- Orders chromosomes: `1A–7A, 1B–7B, 1R–7R`.
- Intersects and synchronizes marker IDs across all files; imputes missing genotypes if needed.
- Iterates over all **traits** (`rest_1`–`rest_5`) and **models** (`GLM`, `MLM`, `CMLM`, etc).
- For each trait/model:
  - Generates Manhattan plots (`Manhattan.pdf`).
  - Outputs lists of significant SNPs (`Significant_SNPs_p001.csv`, `Significant_SNPs_p0005.csv`).
- All results are saved under `GAPIT_results/`.

---

## 3. Post-GWAS: Merging Significant SNP Results

**Script:** `merge_significant_snps.py`

**Steps:**
- Merges significant SNPs across all traits and models into consolidated summary tables.

---

## 4. Annotating Markers with Allele Sequences

**Script:** `add_allele_seq.py`

**Steps:**
- Adds allele sequences to `location.xlsx` for downstream BLAST analysis.

---

## 5. Mapping Marker Sequences via BLAST

**Reference Genomes:**
- **Wheat:** `/ref_wheat/ncbi_dataset/data/GCF_018294505.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna`
- **Rye:** `/ref_rye/ncbi_dataset/data/GCA_965641915.1/GCA_965641915.1_lpSecCere.Lo7.IPK.v3_genomic.fna`

### 5.1 Creating BLAST Databases

**Commands:**
```sh
makeblastdb -in GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna \
  -dbtype nucl \
  -parse_seqids \
  -out wheat_db
makeblastdb -in GCA_902687465.1_Rye_Lo7_2018_v1p1p1_genomic.fna \
  -dbtype nucl \
  -parse_seqids \
  -out rye_db
```

### 5.2 BLAST: Relaxed and Strict Filtering

**Script:** `blast_allele_seq.py`  
- **Relaxed search:**  
  - Tool: `blastn-short`
  - E-value: 1e-3
  - Word size: 7
  - No complexity filter
- **Strict filtering:**  
  - ≥95% identity
  - Alignment length ≥60 bp
  - ≤2 mismatches, no gaps
  - E-value ≤1e-10

**Output:**  
`relaxed_results.xlsx` – best hit per CloneID/marker per genome

### 5.3 Standardizing Chromosome Names

**Script:** `rename_chr_NCBI.py`

**Purpose:**  
- Extracts GenBank accessions from BLAST hits and maps them to chromosome names (e.g., `NC_057794.1` → `2B`).

**Output:**  
`relaxed_results_with_chr_names.xlsx`

---

## 6. Integrating BLAST Hits with Significant SNPs

**Script:** `Significant_SNPs_with_relaxed_hits.py`

**Purpose:**
- Loads `All_Significant_SNPs_p001_merged_python.csv` and `relaxed_results.xlsx`.
- Matches GWAS SNPs to BLAST markers (CloneID).
- Adds genomic coordinates.

**Output:**  
`Significant_SNPs_with_relaxed_hits.xlsx`  
(Main table connecting restitution GWAS hits to wheat/rye genomic positions.)

---

## 7. Linkage Disequilibrium (LD) Analysis

**Scripts:**  
- `v2_LD.R` (genome-wide LD)
- `LD_visualization.R` (LD heatmaps per chromosome)
- `LD_region_haploview_fixed.R` (local LD around markers)
- `LD_blocks_threshold0.8.R` (extract LD blocks)

### 7.1 Genome-wide LD
```sh
nohup Rscript v2_LD.R &
```
### 7.2 Chromosome LD Heatmaps
```sh
Rscript LD_visualization.R
```
### 7.3 LD Around a Marker
```r
source("LD_region_haploview_fixed.R")
plot_ld_region("100008483", window = 30)
```
### 7.4 LD Blocks Extraction
```sh
Rscript LD_blocks_threshold0.8.R
```

**Outputs:**  
- LD matrices  
- LD blocks per chromosome  
- Regional LD structures supporting candidate gene prioritization

---

## 8. Functional Annotation Integration

### 8.1 Wheat Annotation

**Files:**  
- High-confidence gene models: `iwgsc_refseqv2.1_annotation_200916_HC.gff3`
- Functional annotation: `iwgsc_refseqv2.1_functional_annotation.csv`

**Script:**  
`enrich_gff_with_functional_annotation.py`

**Purpose:**  
- Adds GO, Pfam, InterPro, and functional descriptions to the GFF3 file.

### 8.2 Rye Annotation

**Files:**
- Raw annotation: `Secale_cereale.Rye_Lo7_2018_v1p1p1.62.gff3`
- High-confidence + GO: `Secale_cereale_Lo7_2018v1p1p1.pgsb.Feb2019.HC.gff3`

---

## 9. Extracting Candidate Genes Near Significant Loci

**Script:**  
`get_gwas_region_genes.py --chr Chr2B --pos 84439437 --window 500000`

**Output:**  
List of genes within ±500 kb around the SNP, with functional annotations for wheat and rye.

---

## 10. GO Enrichment Analysis for Candidate Regions

**Script:**  
`gwas_chr2B_GO_enrichment.R`

**Outputs:**
- List of genes in the region (with annotation)
- GO enrichment table: `GO_enrichment_Chr2B_84439437_w1000000.tsv`
- Candidate gene descriptions: `candidate_genes_Chr2B_84439437_w1000000.tsv`

---

## 11. Integrating Functional Annotation with Significant SNP Table

**Script:**  
`add_genes_with_GO_to_SNP_table.R`

**Purpose:**  
Adds nearby genes and functional details directly into the GWAS SNP table.  
Final output used for restitution loci interpretation.

---

## 12. Summary Workflow

Here’s a streamlined flow summarizing the full pipeline:

### Workflow Diagram

```
Raw marker + phenotype data
      │
preprocessing (preparing_for_GAPIT.py)
      │
GWAS (all GAPIT models)
      ├─ per-trait/model results
      └─ merged significant SNPs
      │
Add allele sequences (add_allele_seq.py)
      │
BLAST mapping (wheat + rye, blast_allele_seq.py)
      └─ strict filtering & normalization
      │
Chromosome naming (rename_chr_NCBI.py)
      │
Merge BLAST hits with significant SNPs
      │
LD analysis (genome-wide & local, R scripts)
      │
Functional annotation enrichment (Python)
      │
Candidate gene extraction near hits
      │
GO enrichment analysis (R)
      │
Integrated output: SNP + BLAST + LD + functional genes
      │
Ready for interpretation & manuscript prep
```

---


**Contact:**  
For questions, contact shima.mahmoudi@msn.com
