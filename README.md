# GENOME-TO-PROTEIN
R Script which downloads a genome, translates all 6 reading frames into ORF peptides, BLASTs reference protein/exon queries against them, extracts matching ORFs, and reconstructs full-length proteins with coverage, substitutions, and exon-source outputs.
# Genome_To_Protein R Script

## Citation required

If you use this script, workflow, or a modified version of it in published research, reports, teaching material, or derivative analyses, please cite it as:

**van den Hurk, Y. (2026). _Genome_To_Protein: an R workflow for genome-to-protein ORF extraction, exon matching, and protein reconstruction_. GitHub repository.**

A `CITATION.cff` file is provided in this repository. Please use that file for the most up-to-date citation information. If a DOI-linked release is available, please cite the archived release DOI.

---

## Overview

`Genome_To_Protein` is an R-based workflow for identifying and reconstructing protein sequences from genome assemblies. It was designed for cases where a target genome does not have a suitable annotated proteome, but where reference protein sequences are available.

The workflow does three main things:

1. **Translates a genome assembly into six-frame ORF peptides**
2. **Searches reference proteins or exon-level reference sequences against those ORFs using BLASTP**
3. **Reconstructs full-length proteins by placing the best-matching ORF evidence onto reference protein sequences**

The final reconstructed proteins preserve the reference protein length. Residues supported by the target genome are written in **uppercase**, while unsupported positions filled from the reference are written in **lowercase**.

This makes the output useful for checking which regions of a protein are directly supported by genome-derived ORFs and which regions remain reference-based.

---

## Software requirements

Before running the scripts, install the following software.

### Required software

- **R**
- **RStudio** is recommended but not required
- **NCBI BLAST+**
- Internet access for the NCBI genome download step, if using the NCBI assembly downloader

### Required R packages

The scripts use the following R packages:

- `Biostrings`
- `jsonlite`
- `data.table`

Some scripts install missing packages automatically. However, it is safest to install them manually first:

```r
install.packages("jsonlite")
install.packages("data.table")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings")
```

### NCBI BLAST+

The Part A scripts require the following BLAST+ programs:

- `blastp`
- `makeblastdb`
- `blastdbcmd`

These must either be available in your system `PATH`, or you can set the BLAST+ folder manually in R before running the script:

```r
Sys.setenv(NCBI_BLAST_BIN = "C:/path/to/blast-2.xx.x/bin")
```

For example, on Windows this might look like:

```r
Sys.setenv(NCBI_BLAST_BIN = "C:/Users/yourname/Documents/blast-2.17.0/bin")
```

The GitHub-ready versions of the scripts avoid hard-coded personal paths and instead use the current working directory.

---

## Required input files

The reconstruction workflow requires reference protein FASTA files. These are used as the guide for finding and reconstructing proteins from the target genome.

The required files are:

```text
FullProteins.fasta
SequenceProteins.fasta
```

Both should be placed in the same working directory as the script and the generated ORF peptide files.

---

## Reference FASTA files

### 1. `FullProteins.fasta`

This file contains the full-length reference proteins.

These can be downloaded from UniProt or another reference source. The standard UniProt FASTA header format can be kept.

Example:

```text
>sp|P12345|PROTEIN_NAME Species name
MKWVTFISLLFLFSSAYSRGVFR...
```

The script uses these full-length proteins as the final reference framework for reconstruction.

The final reconstructed protein will have the same overall length as the corresponding sequence in `FullProteins.fasta`.

---

### 2. `SequenceProteins.fasta`

This file contains the exon-level or protein-window reference sequences.

These should correspond to smaller sections of the full proteins. Each entry must include the associated protein ID and the exact amino acid coordinates within the full protein.

The expected header format is:

```text
>anything|ProteinID|start-end
SEQUENCE
```

Example:

```text
>sp|P12345|25-63
GEPGPPGPPGPPGLGGNFAPQLSYGYDEKSTGISVPGPM
```

In this example, the sequence represents amino acids 25-63 of protein `P12345`.

This file is very important because it tells the script where each reference exon or protein window belongs within the full-length protein.

The script expects the protein IDs in `SequenceProteins.fasta` to match those in `FullProteins.fasta`.

---

## General workflow

The workflow has three main steps.

---

## Step 1: Genome to ORF peptide translation

The Step 1 script takes a genome assembly and translates it in all six reading frames.

For NCBI genomes, the script can automatically resolve an NCBI assembly accession or NCBI Datasets URL.

Example input:

```r
ASSEMBLY_INPUT <- "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_041834305.1/"
```

or:

```r
ASSEMBLY_INPUT <- "GCA_041834305.1"
```

The script then:

1. Resolves the correct NCBI FTP directory
2. Downloads the `*_genomic.fna.gz` file if needed
3. Reads the genome FASTA file
4. Translates all six reading frames
5. Extracts stop-to-stop ORFs
6. Writes the ORFs to a peptide FASTA file

The resulting file has a name similar to:

```text
ASSEMBLY_NAME_ORFs_min12.pep
```

The ORF FASTA headers are designed to be compatible with BLAST database creation using `makeblastdb -parse_seqids`.

---

## Important Step 1 settings

The main settings are near the top of the Step 1 script.

### `ASSEMBLY_INPUT`

The NCBI accession or NCBI Datasets URL for the target genome.

```r
ASSEMBLY_INPUT <- "GCA_041834305.1"
```

or:

```r
ASSEMBLY_INPUT <- "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_041834305.1/"
```

---

### `MIN_AA`

Minimum ORF length in amino acids.

```r
MIN_AA <- 12
```

A lower value keeps shorter ORFs but creates a larger peptide database. A higher value removes short ORFs but may miss short exons.

---

### `REQUIRE_START_M`

Controls whether ORFs must start with methionine.

```r
REQUIRE_START_M <- FALSE
```

The default is `FALSE`, meaning the script extracts stop-to-stop ORFs. This is often better for genome-derived exon recovery because real exon fragments may not begin with methionine.

---

### `MIN_CONTIG_NT`

Minimum contig length to process.

```r
MIN_CONTIG_NT <- 1000L
```

Shorter contigs are skipped. Set this to `0L` to process all contigs.

---

### `WRITE_GZ_OUTPUT`

Controls whether the ORF peptide output is compressed.

```r
WRITE_GZ_OUTPUT <- FALSE
```

---

## Step 2: Part A - locating matching ORFs/exons

Step 2 searches the ORF peptide database using the reference sequences from `SequenceProteins.fasta`.

There are two similar versions of Step 2:

```text
partA_NCBI_github_ready.R
partA_DNAZOO_CNGB_github_ready.R
```

Both versions perform the same general task, but they differ in how they deal with FASTA headers and BLAST ID extraction.

---

## Step 2 version 1: NCBI genomes

Use the NCBI Part A script for ORF peptide files generated from NCBI genome assemblies, especially those produced by Step 1.

This version assumes that the ORF headers are already pipeline-safe and compatible with BLAST database indexing.

It does the following:

1. Finds `.pep` files in the working directory
2. Builds or verifies BLAST protein databases
3. Splits reference queries into short and long sequences
4. Runs BLASTP for longer queries
5. Runs `blastp-short` for short exon-level queries
6. Merges all BLAST hits
7. Deduplicates the BLAST table
8. Extracts matching ORF sequences using `blastdbcmd`
9. Writes cleaned FASTA headers for Part B

Main output files:

```text
blast_hits_merged.txt
exon_hits_fullcontigs.fasta
exon_hits_fullcontigs_clean.fasta
query_id_map.tsv
```

---

## Step 2 version 2: DNAZOO and CNGB genomes

Use the DNAZOO/CNGB Part A script for genome assemblies or ORF files where FASTA headers may be less predictable.

This version is more defensive and attempts to recover the correct BLAST sequence IDs from multiple fields.

It is useful when BLAST reports truncated or altered IDs.

It does the following:

1. Finds `.pep` files in the working directory
2. Builds or verifies BLAST databases
3. Runs BLASTP and `blastp-short`
4. Reads the BLAST hit table
5. Tests whether sequence IDs can be retrieved using:
   - `sseqid`
   - `sacc`
   - the first token of `stitle`
6. Automatically rebuilds BLAST databases with `-parse_seqids` if needed
7. Extracts the matching ORF sequences
8. Writes a cleaned FASTA file while preserving full scaffold-style IDs

This version is especially useful for assemblies with headers such as:

```text
HiC_scaffold_...
```

or other non-NCBI scaffold naming systems.

Main output files:

```text
blast_hits_merged.txt
exon_hits_fullcontigs.fasta
exon_hits_fullcontigs_clean.fasta
```

---

## Important Step 2 settings

### `data_dir`

The GitHub-ready scripts use the current working directory:

```r
data_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
```

This means you should set your working directory to the folder containing the input files before running the script.

In R:

```r
setwd("path/to/your/project/folder")
```

---

### `query_fasta`

The default query file is:

```r
query_fasta <- file.path(data_dir, "SequenceProteins.fasta")
```

This means BLAST searches are performed using the exon-level or protein-window reference sequences.

There is also usually an optional line in the script for searching full proteins instead:

```r
# query_fasta <- file.path(data_dir, "FullProteins.fasta")
```

For most reconstruction workflows, `SequenceProteins.fasta` should be used.

---

### Short and long query handling

The script splits queries into:

```text
short queries: <=60 amino acids
long queries:  >60 amino acids
```

Long queries are searched with regular BLASTP settings.

Short queries are searched with `blastp-short`, which is more suitable for short exon-level reference sequences.

---

## Step 3: Part B - protein reconstruction

Part B reconstructs full-length proteins from the ORF/exon evidence found in Part A.

It requires:

```text
FullProteins.fasta
SequenceProteins.fasta
exon_hits_fullcontigs_clean.fasta
blast_hits_merged.txt
```

The BLAST hit table is optional in newer versions of the script. If it is missing, the script can attempt reconstruction using pseudo-hits from the extracted exon FASTA file, but the preferred workflow is to keep the BLAST hit table.

Part B does the following:

1. Reads the full-length reference proteins
2. Reads the exon/window reference sequences
3. Reads the ORF sequences extracted in Part A
4. Reads the BLAST hit table, if available
5. Groups candidate ORFs by contig/scaffold and strand
6. Selects the best-supported cluster of ORFs
7. Paints matching ORF evidence onto the full reference protein
8. Uses reference residues where target evidence is missing
9. Reports coverage, substitutions, omitted exons, and ORFs used

---

## Interpretation of reconstructed sequences

The reconstructed FASTA sequences use uppercase and lowercase letters.

### Uppercase residues

Uppercase residues are supported by ORF evidence from the target genome.

Example:

```text
GEPGPPGPPGPPGLGGNFA
```

### Lowercase residues

Lowercase residues are copied from the reference protein because no suitable target ORF evidence was found.

Example:

```text
gepgppgppgppglggnfa
```

### Mixed sequence

A reconstructed sequence may look like this:

```text
MKWVTFISLLflfssaysrGVFRRDAHKSEVAHRFKDLGE
```

This means some parts are supported by the target genome, while other parts remain reference-based.

---

## Important Part B settings

Part B contains many tunable parameters. Most users should start with the default settings.

### `block_min`

Minimum contiguous exon block required for primary cluster selection.

```r
block_min <- 2L
```

Higher values make cluster selection stricter.

---

### `max_reuse_per_sseq`

Maximum number of times a single ORF sequence can be reused for one protein.

```r
max_reuse_per_sseq <- 2L
```

This prevents one ORF from being reused too often across many exons.

---

### `allow_immediate_orf_reuse`

Allows an ORF that was used for one exon to continue into the immediately following exon.

```r
allow_immediate_orf_reuse <- TRUE
```

This is useful when one ORF spans multiple adjacent reference exons.

---

### `boundary_forgive_ends`

Allows mismatches at exon boundaries to be treated more leniently during matching.

```r
boundary_forgive_ends <- 1L
```

---

### `boundary_fill_max`

Allows very small internal gaps to be filled using the reference sequence.

```r
boundary_fill_max <- 2L
```

---

### `nofill_exon_max_len`

Very short reference exons are not filled or rescued.

```r
nofill_exon_max_len <- 5L
```

---

### `restrict_to_single_cluster`

Restricts reconstruction to a single genomic cluster when possible.

```r
restrict_to_single_cluster <- TRUE
```

This helps avoid incorrectly combining ORFs from unrelated regions.

---

### `try_secondary_cluster_if_poor`

Allows the script to test an additional cluster if the primary cluster gives poor coverage.

```r
try_secondary_cluster_if_poor <- TRUE
```

---

### `min_cov_for_primary_accept`

Minimum coverage required to accept the primary cluster without trying alternatives.

```r
min_cov_for_primary_accept <- 0.70
```

---

### `subs_to_retry_threshold`

If the number of substitutions is higher than this threshold, the script may test alternative clusters.

```r
subs_to_retry_threshold <- 20L
```

---

## Forced cluster option

Part B includes an optional forced-cluster setting:

```r
forced_cluster_specs_raw <- c(
  # "ProteinID" = "ORF_ID"
)
```

This is useful after visual inspection of the BLAST results or ORF locations.

For example, if you know that a specific protein should be reconstructed from a particular genomic region, you can force the script to use that region.

Example:

```r
forced_cluster_specs_raw <- c(
  "P12345" = "scaffold_001_pF2_o000123_456789"
)
```

You can also provide two forced clusters for one protein:

```r
forced_cluster_specs_raw <- c(
  "P12345" = "scaffold_001_pF2_o000123_456789 ; scaffold_009_mF1_o000456_987654"
)
```

When two clusters are provided, the script can test a split reconstruction, using one cluster for one part of the protein and the second cluster for another part.

This option should be used carefully. It is intended for cases where the automatic cluster choice has been checked visually and a specific locus is known or strongly suspected to be correct.

---

## Output files

Part B writes several output files. The filenames are based on the `.pep` file used in the working directory.

For example, if the ORF peptide database is:

```text
GCA_000000000.1_ORFs_min12.pep
```

then the main outputs will be:

```text
GCA_000000000.1_ORFs_min12.pep_reconstructed.fasta
GCA_000000000.1_ORFs_min12.pep_reconstruction_coverage.tsv
GCA_000000000.1_ORFs_min12.pep_reconstructed_with_exons.fasta
```

---

### `*_reconstructed.fasta`

This file contains the reconstructed full-length proteins.

Each FASTA header includes summary information such as:

```text
>ProteinID | 12 exons | 94.2% covered | subs=3 | omit(1-15)=0, omit(16-30)=1, omit(31+)=0 | clusters=1 | OK
```

The sequence below the header is the reconstructed protein.

Uppercase residues are supported by ORF evidence.

Lowercase residues are copied from the reference.

---

### `*_reconstruction_coverage.tsv`

This is a tab-separated summary table.

It includes columns such as:

```text
protein
n_exons
covered_percent
subs
omit_1_15
omit_16_30
omit_31p
n_clusters
status
exons_used
```

This file is useful for quickly comparing reconstruction quality across proteins.

---

### `*_reconstructed_with_exons.fasta`

This file contains:

1. The reconstructed full-length protein
2. A list of ORFs/exons used
3. The full sequences of the ORFs/exons used for reconstruction
4. Extra copies of ORFs that introduced substitutions

This file is useful for checking exactly which genome-derived ORFs contributed to each reconstructed protein.

---

## Status labels

The reconstruction header contains a status label.

### `OK`

The reconstruction has an acceptable number of substitutions according to the current settings.

### `FAIL`

The reconstruction has more substitutions than the chosen threshold.

A `FAIL` label does not always mean the result is useless. It means the reconstruction should be checked manually.

---

## Recommended folder structure

A simple working directory should look like this:

```text
project_folder/
├── FullProteins.fasta
├── SequenceProteins.fasta
├── genome_to_protein_step1.R
├── partA_NCBI_github_ready.R
├── partA_DNAZOO_CNGB_github_ready.R
├── partB_protein_reconstruction.R
└── GCA_000000000.1_ORFs_min12.pep
```

After running the full workflow, the folder will also contain files such as:

```text
blast_hits_merged.txt
exon_hits_fullcontigs.fasta
exon_hits_fullcontigs_clean.fasta
*_reconstructed.fasta
*_reconstruction_coverage.tsv
*_reconstructed_with_exons.fasta
```

---

## Which Step 2 version should I use?

Use:

```text
partA_NCBI_github_ready.R
```

for NCBI genome assemblies and ORF peptide files generated by the Step 1 script.

Use:

```text
partA_DNAZOO_CNGB_github_ready.R
```

for DNAZOO, CNGB, or other genome assemblies where scaffold headers may be non-standard, truncated, or difficult for BLAST to retrieve.

The two versions are very similar. The main difference is that the DNAZOO/CNGB version includes extra ID-checking and self-healing steps for difficult FASTA headers.

---

## Notes and cautions

This workflow is intended to assist protein reconstruction from genome assemblies, but the results should be checked carefully.

Important things to inspect include:

- Whether the correct reference proteins were used
- Whether `FullProteins.fasta` and `SequenceProteins.fasta` use matching protein IDs
- Whether exon/window coordinates in `SequenceProteins.fasta` are correct
- Whether the chosen genomic cluster makes biological sense
- Whether reconstructed substitutions are real or caused by assembly errors, paralogy, contamination, or incorrect ORF placement
- Whether forced clusters were used appropriately

For publication-quality analyses, reconstructed proteins should be manually inspected and validated.

---

## License

This software is released under the MIT License.

See the `LICENSE` file for details.

Please note that the license describes the legal terms for reuse of the software. The citation request above describes the expected academic practice when the workflow is used in research.
