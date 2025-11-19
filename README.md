# Triticale Restitution GWAS: Full Pipeline and Workflow Documentation

## Overview

This repository contains the complete workflow for Genome-Wide Association Study (GWAS) on restitution traits in triticale. It covers data preparation, GWAS analysis, postprocessing significant marker results, annotation with sequences, and BLAST for marker mapping to reference genomes. The process utilizes **Python** for data manipulation & BLAST, and **R** for GWAS and linkage disequilibrium calculations.

---

## 1. Data Preparation for GAPIT

**Script:** `preparing_for_second_GAPIT.py`

**Workflow:**
- Reads marker data (`Tabela_S1.xlsx`) and phenotype data (`20250304_pheno_2n_gametes_ABDR.csv`).
- Cleans, harmonizes genotype and phenotype sample IDs, filters by missingness, and drops monomorphic SNPs.
- Encodes genotypes (`a` → `0`, `b` → `2`) and prepares metadata.
- Outputs:
  - `GAPIT_genotype.csv` (GD matrix, Taxa x SNPs)
  - `GAPIT_phenotype.csv` (phenotype matrix)
  - `GAPIT_snp_metadata.csv` (SNP info, CHR/Kosambi cM, etc.)

---

## 2. Running GWAS in R with GAPIT

**Script:** `GWAS_kosmbi_full_pipeline_all_models.R`

**Workflow:**
- Loads cleaned genotype, phenotype, and SNP metadata tables.
- Orders chromosomes for triticale (`1A–7A, 1B–7B, 1R–7R`).
- Intersects and matches marker IDs in all files; imputes missing genotypes.
- Loops through set of **traits** (`rest_1`–`rest_5`) and **models** (`GLM`, `MLM`, `CMLM`, ...).
- Outputs for each trait/model:
  - Manhattan plots: `Manhattan.pdf`
  - Significant SNPs (`p < 0.01`, `< 0.005`): `Significant_SNPs_p001.csv`, `Significant_SNPs_p0005.csv`
- Saves all results under `GAPIT_results/`.

---

## 3. Post-GWAS: Merging Significant SNP Results

**Script:** `merge_significant_snps.py`

**Workflow:**
- Merges significant SNPs across traits/models into comprehensive summary sheets.

---

## 4. Annotating Markers: Adding Allele Sequence Information

**Script:** `add_allele_seq.py`

**Workflow:**
- Adds allele sequences to `location.xlsx` for downstream BLAST.

---

## 5. BLAST Marker Sequences to Reference Genomes

**Reference Genomes:**
- **Wheat:** `/ref_wheat/ncbi_dataset/data/GCF_018294505.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna`
- **Rye:** `/ref_rye/ncbi_dataset/data/GCA_965641915.1/GCA_965641915.1_lpSecCere.Lo7.IPK.v3_genomic.fna`

**BLAST Commands:**
- Create BLAST DB for Wheat:
  ```sh
  makeblastdb \
    -in GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna \
    -dbtype nucl \
    -parse_seqids \
    -out wheat_db
  ```
- Create BLAST DB for Rye:
  ```sh
  makeblastdb \
    -in GCA_965641915.1_lpSecCere.Lo7.IPK.v3_genomic.fna \
    -dbtype nucl \
    -parse_seqids \
    -out rye_db
  ```
- Run BLAST:
  ```
  python3 blast_allele_seq.py
  ```
- For relaxed criteria (e.g., mismapped hits):
  ```
  python3 unmapped.py
  ```

---

## 6. Prepare Data for PLINK & LD Calculation

**Script:** `v2_LD.R`

**Workflow:**
- Reads marker table (`Tabela_S1.csv`), detects individual columns, and retains only mapped markers (with Kosambi cM).
- Recodes A/B/- genotype coding to numeric 0/1/NA. Converts genotype matrix to SnpMatrix (from `snpStats`).
- Calculates LD (R²) matrix for all markers, saves to `triticale_3083markers_r2_matrix.rds`.

---

## Folder Structure

```
├── preparing_for_second_GAPIT.py             # Data preparation for GAPIT
├── GAPIT_genotype.csv                        # Output: Genotype matrix
├── GAPIT_phenotype.csv                       # Output: Phenotype matrix
├── GAPIT_snp_metadata.csv                    # Output: SNP metadata
├── GWAS_kosmbi_full_pipeline_all_models.R    # GAPIT GWAS full pipeline
├── GAPIT_results/                            # All GWAS result subdirectories
│    └── [trait]/[model]/...                  # e.g., rest_5/FarmCPU/Manhattan.pdf
├── merge_significant_snps.py                 # Merges significant SNP results
├── add_allele_seq.py                         # Adds sequences for BLAST/annotation
├── blast_allele_seq.py                       # BLASTs against reference genomes
├── unmapped.py                               # BLAST with relaxed criteria
├── v2_LD.R                                   # Prepares genotype for PLINK/LD
├── ref_wheat/                                # Wheat reference genome files
├── ref_rye/                                  # Rye reference genome files
└── triticale_3083markers_r2_matrix.rds       # Saved LD matrix
```

---

## Dependencies

- **Python 3**: `pandas`, `numpy`
- **R**: `GAPIT`, `ggplot2`, `data.table`, `snpStats`
- **BLAST+**: `makeblastdb`, `blastn`

---

## Running the Pipeline

1. **Data preparation:**
   ```
   python3 preparing_for_second_GAPIT.py
   ```
2. **GWAS run:**
   ```
   Rscript GWAS_kosmbi_full_pipeline_all_models.R
   ```
3. **Merge results:**
   ```
   python3 merge_significant_snps.py
   ```
4. **Add sequence annotation:**
   ```
   python3 add_allele_seq.py
   ```
5. **BLAST mapping:**
   ```
   python3 blast_allele_seq.py
   python3 unmapped.py      # for relaxed
   ```
6. **LD calculation / PLINK prep:**
   ```
   Rscript v2_LD.R
   ```

---

## Notes

- **Chromosome ordering:** Ensures triticale chromosomes match `1A–7A, 1B–7B, 1R–7R`.
- **Genotype encoding:** Always verify `a/b` or `A/B` coding before numeric conversion.
- **Missing data:** Filtering by sample/SNP missingness is crucial for valid GWAS results.
- **BLAST usage:** Marker position validation on reference sequence; use stringent/relaxed criteria as needed.
- **LD analysis:** Essential for exploring marker correlations and GWAS peak validation.

---

## Contact

For questions or troubleshooting, please raise an issue or contact the repository maintainer.
