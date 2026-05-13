#!/usr/bin/env Rscript
# ============================================================
# PART A – NCBI GENOMES – Locating matching exons
# For genome ORF .pep files produced by Step 1
# GitHub-ready version: no hard-coded C:/Users/... paths.
#
# How to use:
#   1) Put this script in the same folder as:
#        - SequenceProteins.fasta
#        - one or more .pep ORF databases, e.g. *_ORFs_min12.pep
#   2) Make sure NCBI BLAST+ is installed.
#      Either add BLAST+ to your PATH, or set BLAST_BIN below / env var NCBI_BLAST_BIN.
#   3) In R: setwd("path/to/your/project/folder")
#   4) source("partA_NCBI_github_ready.R")
#
# Outputs:
#   - blast_hits_merged.txt
#   - exon_hits_fullcontigs.fasta
#   - exon_hits_fullcontigs_clean.fasta
#   - query_id_map.tsv
# ============================================================

# ---- 1) Portable paths ----
# Project/data folder. Default is current working directory.
data_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# BLAST+ location.
# Recommended: leave as "" if blastp/makeblastdb/blastdbcmd are on PATH.
# Alternatively, set it here, e.g. BLAST_BIN <- "C:/Program Files/NCBI/blast-2.17.0+/bin"
# Or set environment variable NCBI_BLAST_BIN before running R.
BLAST_BIN <- Sys.getenv("NCBI_BLAST_BIN", unset = "")

find_blast_exe <- function(exe_name, blast_bin = BLAST_BIN) {
  exe_file <- if (.Platform$OS.type == "windows" && !grepl("\\.exe$", exe_name, ignore.case = TRUE)) {
    paste0(exe_name, ".exe")
  } else {
    exe_name
  }

  if (nzchar(blast_bin)) {
    candidate <- file.path(blast_bin, exe_file)
    if (file.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    stop("Could not find ", exe_file, " in BLAST_BIN / NCBI_BLAST_BIN:\n  ", blast_bin,
         "\nEither correct BLAST_BIN or add BLAST+ to your PATH.", call. = FALSE)
  }

  found <- Sys.which(exe_file)
  if (nzchar(found)) return(normalizePath(found, winslash = "/", mustWork = TRUE))

  stop("Could not find ", exe_file, " on PATH.\n",
       "Install NCBI BLAST+ and add it to PATH, or set BLAST_BIN / NCBI_BLAST_BIN.", call. = FALSE)
}

blastp      <- find_blast_exe("blastp")
makeblastdb <- find_blast_exe("makeblastdb")
blastdbcmd  <- find_blast_exe("blastdbcmd")

# Queries
query_fasta  <- file.path(data_dir, "SequenceProteins.fasta")
# query_fasta <- file.path(data_dir, "FullProteins.fasta")

merged_hits   <- file.path(data_dir, "blast_hits_merged.txt")
merged_fasta  <- file.path(data_dir, "exon_hits_fullcontigs.fasta")
clean_fasta   <- file.path(data_dir, "exon_hits_fullcontigs_clean.fasta")
query_map_tsv <- file.path(data_dir, "query_id_map.tsv")

# ---- 2) Helpers ----
norm_path <- function(p) {
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

run_cmd <- function(prog, args, label = NULL, stop_on_error = TRUE) {
  if (!is.null(label)) cat(label, "\n")
  out <- system2(prog, args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  if (status != 0L) {
    cat("❌ Command failed (exit status ", status, ")\n", sep = "")
    cat("Program: ", prog, "\nArgs:\n  ", paste(args, collapse = " "), "\n\n", sep = "")
    cat("--- Output (stdout+stderr) ---\n")
    if (length(out)) {
      cat(paste(out, collapse = "\n"), "\n")
    } else {
      cat("(no stdout/stderr captured)\n")
    }
    if (isTRUE(stop_on_error)) {
      stop("Error: Stopping due to failed command above.", call. = FALSE)
    }
  }

  invisible(list(status = status, output = out))
}

detect_db_peps <- function(dirp) {
  parts <- list.files(dirp, pattern = "^part_\\d+[a-z]?\\.pep$", full.names = TRUE, ignore.case = TRUE)
  if (length(parts)) return(parts)

  orfs <- list.files(dirp, pattern = "_ORFs_min[0-9]+.*\\.pep$", full.names = TRUE, ignore.case = TRUE)
  if (length(orfs)) return(orfs)

  peps <- list.files(dirp, pattern = "\\.pep$", full.names = TRUE, ignore.case = TRUE)
  if (length(peps)) return(peps)

  character(0)
}

db_is_indexed <- function(db_prefix_path) {
  dirn <- dirname(db_prefix_path)
  base <- basename(db_prefix_path)

  esc <- gsub("([.\\+\\-\\^\\$\\(\\)\\[\\]\\{\\}\\|\\?\\*])", "\\\\\\1", base)
  pat_vol <- paste0("^", esc, "(\\.[0-9]{2})?\\.(pin|phr|psq|pog|psd|psi)$")

  idx <- list.files(dirn, pattern = pat_vol, full.names = TRUE, ignore.case = TRUE)
  pal <- file.path(dirn, paste0(base, ".pal"))
  length(idx) > 0L || file.exists(pal)
}

get_db_version <- function(db_prefix, blastdbcmd) {
  out <- system2(
    blastdbcmd,
    c("-db", shQuote(db_prefix, type = "cmd"), "-info"),
    stdout = TRUE, stderr = TRUE
  )
  txt <- paste(out, collapse = "\n")
  m <- regexpr("BLASTDB Version:\\s*([0-9]+)", txt, perl = TRUE)
  if (m[1] == -1L) return(NA_integer_)
  ver <- sub(".*BLASTDB Version:\\s*([0-9]+).*", "\\1", txt, perl = TRUE)
  suppressWarnings(as.integer(ver))
}

remove_existing_db_indexes <- function(db_prefix) {
  dirn <- dirname(db_prefix)
  base <- basename(db_prefix)

  esc <- gsub("([.\\+\\-\\^\\$\\(\\)\\[\\]\\{\\}\\|\\?\\*])", "\\\\\\1", base)
  pat <- paste0("^", esc, "(\\.[0-9]{2})?\\.(pin|phr|psq|pog|psd|psi|pto|ptf|pdb|pjs|pni|pnd|pos|pot)$")

  idx_files <- list.files(dirn, pattern = pat, full.names = TRUE, ignore.case = TRUE)
  pal_file <- file.path(dirn, paste0(base, ".pal"))
  to_remove <- c(idx_files, pal_file[file.exists(pal_file)])

  if (length(to_remove)) {
    cat("🗑️ Removing old DB index files for:", base, "\n")
    ok <- file.remove(to_remove)
    if (!all(ok)) warning("Some old DB index files could not be removed.")
  }
}

ensure_v4_db <- function(db_pep, makeblastdb, blastdbcmd) {
  db_prefix <- norm_path(sub("\\.pep$", "", db_pep, ignore.case = TRUE))

  indexed <- db_is_indexed(db_prefix)
  version <- if (indexed) get_db_version(db_prefix, blastdbcmd) else NA_integer_

  if (indexed && isTRUE(version == 4L)) {
    cat("✅ DB already indexed as VERSION 4:", basename(db_pep), "\n")
    return(invisible(db_prefix))
  }

  if (indexed) {
    cat("⚠️ DB exists but is not VERSION 4. Rebuilding:", basename(db_pep), "\n")
    remove_existing_db_indexes(db_prefix)
  } else {
    cat("🔧 Building BLAST database VERSION 4 for:", basename(db_pep), "\n")
  }

  run_cmd(
    makeblastdb,
    c("-in", shQuote(db_pep, type = "cmd"),
      "-dbtype", "prot",
      "-parse_seqids",
      "-blastdb_version", "4",
      "-out", shQuote(db_prefix, type = "cmd")),
    stop_on_error = TRUE
  )

  if (!db_is_indexed(db_prefix)) {
    stop("❌ makeblastdb did not produce usable index for: ", basename(db_pep), call. = FALSE)
  }

  new_ver <- get_db_version(db_prefix, blastdbcmd)
  cat("✅ DB ready:", basename(db_pep), " | version=", new_ver, "\n", sep = "")
  invisible(db_prefix)
}

# ---- 3) Detect DB inputs ----
db_parts <- detect_db_peps(data_dir)

if (!length(db_parts)) {
  stop("❌ No .pep databases found in:\n   ", data_dir,
       "\nPut your *_ORFs_min*.pep or part_*.pep there and re-run.", call. = FALSE)
}

db_parts <- norm_path(db_parts)

cat("📂 Data directory:\n   ", data_dir, "\n", sep = "")
cat("📂 Detected", length(db_parts), "pep DB file(s):\n")
cat(paste0("   - ", basename(db_parts)), sep = "\n")
cat("\n")

# ---- 4) Verify tools & inputs ----
cat("🔍 Checking executables and query file...\n")
cat("   blastp:      ", blastp, "\n", sep = "")
cat("   makeblastdb: ", makeblastdb, "\n", sep = "")
cat("   blastdbcmd:  ", blastdbcmd, "\n", sep = "")

for (exe in c(blastp, makeblastdb, blastdbcmd)) {
  if (!file.exists(exe)) stop("❌ Missing executable: ", exe, call. = FALSE)
}
if (!file.exists(query_fasta)) stop("❌ Missing query FASTA: ", query_fasta, call. = FALSE)

blastp       <- norm_path(blastp)
makeblastdb  <- norm_path(makeblastdb)
blastdbcmd   <- norm_path(blastdbcmd)
query_fasta  <- norm_path(query_fasta)
merged_hits  <- norm_path(merged_hits)
merged_fasta <- norm_path(merged_fasta)
clean_fasta  <- norm_path(clean_fasta)
query_map_tsv <- norm_path(query_map_tsv)
data_dir     <- norm_path(data_dir)

cat("✅ All executables and query verified.\n\n")

# ---- 5) Ensure BLAST DB VERSION 4 ----
cat("🧱 Verifying BLAST DB indexes...\n")
for (db_pep in db_parts) {
  ensure_v4_db(db_pep, makeblastdb, blastdbcmd)
}
cat("\n")

# ---- 6) Split queries by length, retaining safe IDs internally ----
suppressPackageStartupMessages(library(Biostrings))

qs <- readAAStringSet(query_fasta)
qlen <- width(qs)

orig_names <- names(qs)
if (is.null(orig_names)) orig_names <- rep("", length(qs))

bad <- is.na(orig_names) | !nzchar(trimws(orig_names))
if (any(bad)) {
  orig_names[bad] <- paste0("unnamed_query_", seq_len(sum(bad)))
}

safe_names <- sprintf("Q%06d", seq_along(qs))
names(qs) <- safe_names

query_map <- data.frame(
  safe_qseqid     = safe_names,
  original_qseqid = orig_names,
  length_aa       = as.integer(qlen),
  stringsAsFactors = FALSE
)
write.table(query_map, query_map_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

short_idx <- which(qlen <= 60)
long_idx  <- which(qlen > 60)

q_short_fa <- norm_path(file.path(data_dir, "queries_short_le60_SAFEIDS.faa"))
q_long_fa  <- norm_path(file.path(data_dir, "queries_long_gt60_SAFEIDS.faa"))

if (file.exists(q_short_fa)) file.remove(q_short_fa)
if (file.exists(q_long_fa))  file.remove(q_long_fa)

if (length(short_idx)) writeXStringSet(qs[short_idx], q_short_fa, width = 60)
if (length(long_idx))  writeXStringSet(qs[long_idx],  q_long_fa,  width = 60)

cat("🧪 Query split:\n")
cat("   short (<=60 aa): ", length(short_idx), "\n", sep = "")
cat("   long  (>60 aa):  ", length(long_idx),  "\n", sep = "")
if (length(short_idx)) cat("   short FASTA: ", q_short_fa, "\n", sep = "")
if (length(long_idx))  cat("   long  FASTA: ", q_long_fa,  "\n", sep = "")
cat("\n")

# ---- 7) Run BLASTP: normal long pass + short-query rescue pass ----
cat("🚀 Running BLASTP across database(s)...\n")
if (file.exists(merged_hits)) file.remove(merged_hits)

outfmt_arg <- "6 qseqid sseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"

n_threads <- parallel::detectCores(logical = TRUE)
if (is.na(n_threads) || n_threads < 1L) n_threads <- 4L
cat("⚙️  Using BLAST threads:", n_threads, "\n\n")

for (db_pep in db_parts) {
  db_prefix <- norm_path(sub("\\.pep$", "", db_pep, ignore.case = TRUE))

  if (length(long_idx)) {
    out_long <- norm_path(paste0(db_prefix, "_hits_long.txt"))
    if (file.exists(out_long)) file.remove(out_long)

    cat("🔬 [NORMAL] Processing:", basename(db_pep), "on", length(long_idx), "long queries (>60 aa)\n")

    blast_args_long <- c(
      "-query", shQuote(q_long_fa, type = "cmd"),
      "-db",    shQuote(db_prefix, type = "cmd"),
      "-out",   shQuote(out_long, type = "cmd"),
      "-evalue", "1e-3",
      "-word_size", "2",
      "-seg", "yes",
      "-comp_based_stats", "0",
      "-matrix", "BLOSUM62",
      "-gapopen", "11",
      "-gapextend", "1",
      "-outfmt", shQuote(outfmt_arg, type = "cmd"),
      "-max_target_seqs", "500",
      "-num_threads", as.character(n_threads)
    )

    run_cmd(blastp, blast_args_long, stop_on_error = TRUE)

    if (file.exists(out_long) && file.info(out_long)$size > 0) {
      file.append(merged_hits, out_long)
      cat("   ✅ added normal hits\n")
    } else {
      cat("   ⚠️ no normal hits for this DB\n")
    }
  }

  if (length(short_idx)) {
    out_short <- norm_path(paste0(db_prefix, "_hits_short.txt"))
    if (file.exists(out_short)) file.remove(out_short)

    cat("🛟 [SHORT]  Processing:", basename(db_pep), "on", length(short_idx), "short queries (<=60 aa)\n")

    blast_args_short <- c(
      "-query", shQuote(q_short_fa, type = "cmd"),
      "-db",    shQuote(db_prefix, type = "cmd"),
      "-out",   shQuote(out_short, type = "cmd"),
      "-task", "blastp-short",
      "-seg",  "no",
      "-comp_based_stats", "0",
      "-matrix", "PAM30",
      "-evalue", "50",
      "-word_size", "2",
      "-gapopen", "9",
      "-gapextend", "1",
      "-outfmt", shQuote(outfmt_arg, type = "cmd"),
      "-max_target_seqs", "20000",
      "-num_threads", as.character(n_threads)
    )

    run_cmd(blastp, blast_args_short, stop_on_error = TRUE)

    if (file.exists(out_short) && file.info(out_short)$size > 0) {
      file.append(merged_hits, out_short)
      cat("   ✅ added short-exon rescue hits\n")
    } else {
      cat("   ⚠️ no short-exon hits for this DB\n")
    }
  }
}

cat("\n✅ All BLAST runs completed.\nMerged hits saved at:\n", merged_hits, "\n", sep = "")

# ---- 8) Deduplicate merged hits and write Part-B-compatible 14-column table ----
if (!file.exists(merged_hits) || file.info(merged_hits)$size == 0) {
  stop("❌ No BLAST hits found. Nothing to extract.", call. = FALSE)
}

hits <- read.table(merged_hits, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", fill = TRUE)

if (ncol(hits) < 14L) {
  stop("❌ blast_hits_merged.txt has fewer than 14 columns (found ", ncol(hits), ").", call. = FALSE)
}
if (ncol(hits) > 14L) {
  hits <- hits[, seq_len(14), drop = FALSE]
}

colnames(hits) <- c("qseqid","sseqid","sacc","stitle","pident","length","mismatch","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore")

hits$sseqid <- sub("([|\\s].*)$", "", hits$sseqid)

qmap <- setNames(query_map$original_qseqid, query_map$safe_qseqid)
hits$qseqid <- unname(qmap[hits$qseqid])

bad_map <- is.na(hits$qseqid) | !nzchar(hits$qseqid)
if (any(bad_map)) {
  stop("❌ Some BLAST qseqid values could not be mapped back to original query names.", call. = FALSE)
}

suppressWarnings({
  hits$pident   <- as.numeric(hits$pident)
  hits$length   <- as.integer(hits$length)
  hits$mismatch <- as.integer(hits$mismatch)
  hits$gapopen  <- as.integer(hits$gapopen)
  hits$qstart   <- as.integer(hits$qstart)
  hits$qend     <- as.integer(hits$qend)
  hits$sstart   <- as.integer(hits$sstart)
  hits$send     <- as.integer(hits$send)
  hits$evalue   <- as.numeric(hits$evalue)
  hits$bitscore <- as.numeric(hits$bitscore)
})

hits <- unique(hits)

hits_out <- hits[, c("qseqid","sseqid","sacc","stitle","pident","length","mismatch","gapopen",
                     "qstart","qend","sstart","send","evalue","bitscore")]

write.table(hits_out, file = merged_hits, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE, na = "")

cat("🧼 De-duplicated merged hits. Unique rows:", nrow(hits_out), "\n")
cat("🧩 Wrote Part-B-compatible 14-column BLAST table.\n\n")

# ---- 9) Extract sequences for BLAST hits ----
cat("📦 Extracting hit sequences (batch mode)...\n")
hit_ids <- unique(hits_out$sseqid)
cat("🧩 Unique hits to extract:", length(hit_ids), "\n")

if (file.exists(merged_fasta)) file.remove(merged_fasta)

id_file <- norm_path(file.path(data_dir, "blast_hit_ids.txt"))
writeLines(hit_ids, id_file)

for (db_pep in db_parts) {
  db_prefix <- norm_path(sub("\\.pep$", "", db_pep, ignore.case = TRUE))
  tmp_out <- norm_path(file.path(data_dir, paste0(basename(db_prefix), "_extract_tmp.fasta")))

  if (file.exists(tmp_out)) file.remove(tmp_out)

  out <- system2(
    blastdbcmd,
    c("-db", shQuote(db_prefix, type = "cmd"),
      "-entry_batch", shQuote(id_file, type = "cmd"),
      "-out", shQuote(tmp_out, type = "cmd")),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  if (status != 0L) {
    cat("⚠️ blastdbcmd failed for DB:", basename(db_pep), "\n")
    if (length(out)) cat(paste(out, collapse = "\n"), "\n")
    cat("   ↪ continuing to next DB\n")
    next
  }

  if (file.exists(tmp_out) && file.info(tmp_out)$size > 0) {
    file.append(merged_fasta, tmp_out)
    cat("✅ extracted from", basename(db_pep), "\n")
  } else {
    cat("⚠️ nothing extracted from", basename(db_pep), "\n")
  }
}

if (!file.exists(merged_fasta) || file.info(merged_fasta)$size == 0) {
  stop("❌ blastdbcmd produced no sequences. Check DB indexing and that sseqids exist in the DB.", call. = FALSE)
}

cat("\n🎉 All sequences extracted.\nRaw FASTA saved at:\n", merged_fasta, "\n", sep = "")

# ---- 10) Clean headers ----
gx_raw <- readAAStringSet(merged_fasta)
names(gx_raw) <- sub("^([^\\s]+).*", "\\1", names(gx_raw))
writeXStringSet(gx_raw, clean_fasta, width = 60)

cat("🧽 Clean FASTA written:\n", clean_fasta, "\n", sep = "")
cat("\n✅ Part A finished successfully.\n")
