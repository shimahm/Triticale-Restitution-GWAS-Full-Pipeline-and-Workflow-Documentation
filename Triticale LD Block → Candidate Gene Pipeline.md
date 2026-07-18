# Triticale LD Block → Candidate Gene Pipeline

This repository documents the workflow used to go from **genetic-map-based LD
blocks** in a triticale (wheat × rye) mapping population, to their
**estimated physical location** in the wheat and rye reference genomes, to
the **GWAS significant SNPs** that fall inside them, to the **annotated
genes** in those regions, and finally to a **candidate-gene shortlist**
based on known Arabidopsis meiosis genes.

## Why this pipeline exists

LD blocks in this dataset were originally defined purely from the **genetic
map** (Kosambi cM), by breaking the marker order wherever adjacent-marker r²
dropped below a threshold. That's useful for describing recombination
structure, but a geneticist/breeder ultimately needs to know **where a block
sits on the physical genome** (bp) so it can be matched to gene annotation.
The physical position isn't directly available for every marker (BLAST hits
against the reference are missing for many), and even where hits exist,
wheat and rye's polyploid/duplicated genome structure means a single marker
can BLAST equally well to more than one locus. Each stage below exists to
solve one part of that gap.

## Pipeline overview

```
LD_blocks_threshold0_9.csv  ─┐
                              ├─▶ estimate_block_locations.py ─▶ LD_blocks_estimated_locations.csv
location_with_alleles_blast.xlsx ┘        (naive min/max over all member markers)
                                                    │
                                                    ▼
                                     real_block_location.py
                              (majority-cluster physical position, reduces
                               spurious overlap from mismapped/homoeologous hits)
                                                    │
                                                    ▼
                          LD_blocks_threshold0_9_real_location.csv
                                                    │
                                    add_chrom_names.py (+ chr_NCBI_name_v1.xlsx)
                                                    │
                                                    ▼
                    LD_blocks_threshold0_9_real_location_named.csv
                                                    │
                                       (manual review / cleanup by user)
                                                    │
                                                    ▼
        manual_edited_LD_blocks_threshold0_9_real_location_named.csv
                                                    │
                    ┌───────────────────────────────┴───────────────────────────────┐
                    ▼                                                               ▼
     assign_snps_to_blocks_v2.py                                    find_genes_in_blocks_v2.py
  (+ Significant_SNPs_p0_01.csv,                              (+ wheat_annotated.gff3, rye_annotated.gff3)
     location_with_alleles_blast.xlsx)
                    │                                                               │
                    ▼                                                               ▼
significant_SNPs_with_LD_block_v2.csv                                  LD_block_genes_full.csv
                    │                                                               │
                    └───────────────────────────────┬───────────────────────────────┘
                                                      ▼
                                         find_meiosis_genes.py
                                      (+ At_Meiosis_geneID.txt)
                                                      ▼
                                       LD_block_meiosis_genes.csv
```

A separate exploratory branch (`build_physical_blocks.py`) was also tried,
which rebuilds the blocks from scratch sorted by physical position instead
of genetic position — see [Alternative approach considered](#alternative-approach-considered-physical-position-based-blocks)
below for why the cM-based blocks + majority-cluster location were used
going forward instead.

---

## Step-by-step

### 1. Starting inputs

| File | What it is |
|---|---|
| `LD_blocks_threshold0_9.csv` | LD blocks as originally defined: markers on each chromosome ordered by genetic position (Kosambi cM), broken into blocks wherever adjacent-marker r² < 0.9. Columns: `chrom`, `block_id`, `start_marker`, `end_marker`, `n_snps`, `start_pos_cM`, `end_pos_cM`. |
| `LD_blocks_threshold0_8.R` | The original R script that built the blocks above (kept for reference/reproducibility). |
| `location_with_alleles_blast.xlsx` | Per-marker map + BLAST annotation: `CloneID`, `Nr`, `Chromosom`, `Kosambi` (genetic map info) plus `wheat_Chr/Start/End/Pos` and `rye_Chr/Start/End/Pos` (physical BLAST hit against each reference genome, with hit quality columns `wheat_pident`, `wheat_evalue`, etc.). |
| `triticale_3083markers_r2_matrix.rds` | Full pairwise marker × marker LD (r²) matrix underlying the block-building script. |
| `Tabela_S1.xlsx` | Underlying map/genotype table (same core marker/map columns as the location file, used as the R script's direct input). |

### 2. Estimate each block's physical location from *all* its member markers

**Script:** `estimate_block_locations.py`
**Input:** `LD_blocks_threshold0_9.csv`, `location_with_alleles_blast.xlsx`
**Output:** `LD_blocks_estimated_locations.csv`

**Why:** the LD block table only stores the two *boundary* markers
(`start_marker`, `end_marker`) — not the markers in between — and many
boundary markers have no BLAST hit at all. Since blocks are just a
contiguous run of markers ordered by cM, every marker on that chromosome
whose `Kosambi` falls between `start_pos_cM` and `end_pos_cM` is a genuine
member of the block, whether or not it happens to be a stored boundary
marker. This step reconstructs the full membership list and estimates each
block's physical span as the min/max BLAST position across **all**
members, not just the two boundary ones. (Sanity check: reconstructed
membership count matched `n_snps` exactly for all 169 blocks.)

### 3. Diagnose and reduce block-to-block overlap

**Problem found:** many adjacent blocks' estimated physical ranges
overlapped (103/166 wheat pairs, 83/153 rye pairs). Root cause: min/max is
extremely sensitive to a single stray BLAST hit, and wheat/rye's duplicated
(homoeologous) genome structure means some markers genuinely BLAST with
high identity to more than one locus — only one of which is "real."

**Script:** `real_block_location.py`
**Input:** `LD_blocks_threshold0_9.csv`, `location_with_alleles_blast.xlsx`
**Output:** `LD_blocks_threshold0_9_real_location.csv`

**Approach:** for each block, cluster its member markers' physical
positions by proximity (within 20 Mb = same cluster) and report only the
**largest agreeing cluster**, discarding smaller/stray clusters as likely
mismapped or homoeologous. This cut overlaps roughly 3-fold (103→34 wheat
pairs, 83→33 rye pairs). Remaining overlap reflects genuine genetic
map / physical assembly discordance, not a data artifact — it can't be
filtered away further without additional information.

*(An exploratory alternative, `estimate_block_locations.py`'s min/max and a
MAD-based outlier trim, were tested first; the gap-based clustering above
gave the best overlap reduction and is the version carried forward.)*

### 4. Translate BLAST accession codes to friendly chromosome names

**Script:** `add_chrom_names.py`
**Input:** `LD_blocks_threshold0_9_real_location.csv`, `chr_NCBI_name_v1.xlsx`
**Output:** `LD_blocks_threshold0_9_real_location_named.csv`

**Why:** the BLAST hit columns store raw NCBI/ENA accessions
(`ref|NC_057794.1|` for wheat, `emb|OZ284382.1|` for rye). `chr_NCBI_name_v1.xlsx`
maps each bare accession (`NC_057794.1`, `OZ284382.1`, ...) to its
conventional chromosome name (`1A`, `1R`, ...). This step strips the
`ref|`/`emb|` wrapper and joins on the bare accession to add
`wheat_chrom_name` / `rye_chrom_name` columns.

### 5. Manual review

**File:** `manual_edited_LD_blocks_threshold0_9_real_location_named.csv`

The named output above was manually reviewed and cleaned up (simplified to
one `real_start`/`real_end` pair per block, since each block's `chrom`
label already indicates which genome — `A`/`B` = wheat, `R` = rye — so a
separate wheat/rye column pair per row is redundant). **This manually
finalized file is the block-location reference table used by every
downstream step.**

### 6. Assign GWAS-significant SNPs to their LD block

**Script:** `assign_snps_to_blocks_v2.py`
**Input:** `Significant_SNPs_p0_01.csv`, `location_with_alleles_blast.xlsx`, `manual_edited_LD_blocks_threshold0_9_real_location_named.csv`
**Output:** `significant_SNPs_with_LD_block_v2.csv`

**Why:** the GWAS output identifies significant markers by an internal `SNP`
ID and a `CHR` chromosome label, but doesn't know which LD block each one
falls in. This step:
1. Maps `SNP` → `Nr` in the location file to get each SNP's `CloneID` and
   genetic (`Kosambi`) position.
2. Cross-checks the GWAS file's own `CHR` against the location file's
   `Chromosom` for consistency (flags any disagreement — none found).
3. Restricts the block search to blocks on that same `CHR`, then finds the
   one block whose `[start_pos_cM, end_pos_cM]` range contains the SNP's
   position.

All 29 significant SNPs matched cleanly to exactly one block.

### 7. Find genes overlapping each block's physical location

**Scripts:** `find_genes_in_blocks.py` (first pass, key columns only) →
`find_genes_in_blocks_v2.py` (**final** — keeps every column)
**Input:** `significant_SNPs_with_LD_block_v2.csv`, `wheat_annotated.gff3`, `rye_annotated.gff3`
**Output:** `LD_block_genes_full.csv`

**Why/how:** for each unique block (deduplicated across SNPs that share a
block), the chromosome suffix decides which genome's annotation to search
(`A`/`B` → `wheat_annotated.gff3`, `R` → `rye_annotated.gff3`), with the
chromosome label converted to match each file's own naming convention
(`5A` → `Chr5A` for wheat, `4R` → `chr4R` for rye). All `gene`-type
features whose coordinates overlap `block_real_start`–`block_real_end` are
returned. Every fixed GFF3 column (`seqid`, `source`, `type`, `start`,
`end`, `score`, `strand`, `phase`) and every individual attribute key
(`ID`, `Name`, `Ortholog`, `Symbol`, `GO`, `RBH_pident`, etc. for wheat;
`primary_confidence_class`, `Ontology_term`, etc. for rye) is kept — not
just a hand-picked subset — because the rye file's `description` field
lives on the child `mRNA` line rather than the `gene` line, so that's
pulled up and merged in too (`mRNA_*` columns).

### 8. Shortlist candidate genes using known Arabidopsis meiosis genes

**Script:** `find_meiosis_genes.py`
**Input:** `LD_block_genes_full.csv`, `At_Meiosis_geneID.txt`
**Output:** `LD_block_meiosis_genes.csv`

**Why:** `At_Meiosis_geneID.txt` is a curated list of 171 Arabidopsis genes
known to function in meiosis. The wheat/rye `Ortholog` column stores each
gene's closest Arabidopsis match as a transcript ID (e.g. `AT1G67370.1`) —
the `.1` transcript suffix is stripped before matching against the gene-ID
list. This step flags any block gene whose stripped ortholog ID appears in
the meiosis list, giving a direct answer to "does this LD block contain a
known meiosis gene." Result: 27 matching gene records, covering 23 of the
171 meiosis genes, across 8 distinct blocks.

**Caveat carried into interpretation:** 20 of those 27 hits fall in a
single very large block (5A block 3, spanning ~405–700 Mb) — a block that
size is statistically likely to contain *some* meiosis gene by chance
simply from covering so much of the chromosome, so those hits are weaker
evidence than the ones in the smaller, tighter blocks (2B, 3B, 6B, 4R, 5R,
7R).

---

## Alternative approach considered: physical-position-based blocks

**Script:** `build_physical_blocks.py` (parameterized for `GENOME = 'wheat'`
or `'rye'`)
**Input:** `triticale_3083markers_r2_matrix.rds` (via the `rdata` Python
package — `pyreadr` cannot parse a bare R matrix object), `location_with_alleles_blast.xlsx`
**Output:** `LD_blocks_physical_wheat.csv`, `LD_blocks_physical_rye.csv`

Instead of keeping the cM-defined blocks, this rebuilds LD blocks from
scratch by sorting markers by **physical** position per chromosome
accession and applying the same adjacent-r² break rule used in the
original R script. This guarantees zero overlap **by construction** (each
block is a contiguous physical interval), but:
- only markers with a physical hit on that genome can be included (1702/3083
  for wheat, 1085/3083 for rye),
- it produces two separate, mutually incompatible block sets (one per
  genome), and
- many resulting blocks are singletons, since markers that were near each
  other genetically often turn out to be far apart physically.

This was evaluated but **not** used for the final gene-mapping steps —
the cM-based blocks (Step 3 onward) were kept because they preserve the
original LD-block definition the analysis is built around, and the
majority-cluster location estimate handles the overlap problem well enough
for the candidate-gene use case.

---

## File index

| File | Type | Produced by |
|---|---|---|
| `LD_blocks_threshold0_9.csv` | input | (external R analysis, `LD_blocks_threshold0_8.R`) |
| `location_with_alleles_blast.xlsx` | input | (external BLAST/map pipeline) |
| `triticale_3083markers_r2_matrix.rds` | input | (external LD analysis) |
| `Tabela_S1.xlsx` | input | (external map/genotype table) |
| `chr_NCBI_name_v1.xlsx` | input | (manually compiled accession→name lookup) |
| `Significant_SNPs_p0_01.csv` | input | (external GWAS analysis) |
| `wheat_annotated.gff3`, `rye_annotated.gff3` | input | (external genome annotation) |
| `At_Meiosis_geneID.txt` | input | (curated Arabidopsis meiosis gene list) |
| `LD_blocks_estimated_locations.csv` | intermediate | `estimate_block_locations.py` |
| `LD_blocks_threshold0_9_real_location.csv` | intermediate | `real_block_location.py` |
| `LD_blocks_threshold0_9_real_location_named.csv` | intermediate | `add_chrom_names.py` |
| `manual_edited_LD_blocks_threshold0_9_real_location_named.csv` | **reference table** | manual review |
| `significant_SNPs_with_LD_block_v2.csv` | output | `assign_snps_to_blocks_v2.py` |
| `LD_block_genes_full.csv` | output | `find_genes_in_blocks_v2.py` |
| `LD_block_meiosis_genes.csv` | **final output** | `find_meiosis_genes.py` |
| `LD_blocks_physical_wheat.csv`, `LD_blocks_physical_rye.csv` | exploratory (not used downstream) | `build_physical_blocks.py` |

## Requirements

```
pandas
numpy
openpyxl      # reading .xlsx files
rdata         # reading the .rds LD matrix (pyreadr cannot parse a bare matrix)
```

## Known limitations

- Physical block boundaries are **estimates**, not exact — they depend on
  which member markers happen to have a BLAST hit, and on the majority-
  cluster heuristic used to reject mismapped/homoeologous hits.
- Some overlap between adjacent blocks' physical ranges remains even after
  clustering — this reflects genuine disagreement between the genetic map
  and the physical assembly in those regions, not a bug.
- Gene lists for very large blocks (spanning tens to hundreds of Mb) should
  be treated as low-resolution/weak evidence, not a tight candidate-gene
  set.
