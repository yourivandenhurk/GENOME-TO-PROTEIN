#!/usr/bin/env Rscript
# ============================================================
# PART A – DNAZOO & CNGB GENOMES – Locating matching exons
# Robust + self-healing version for genome ORF .pep files
# GitHub-ready version: no hard-coded C:/Users/... paths.
#
# Fixes/features:
#   - Probes blastdbcmd using IDs from sseqid, sacc, and stitle token 1
#   - If none match, auto-rebuilds BLAST DBs with -parse_seqids and retries
#   - Never stops early if probe fails; attempts extraction and only then stops
#   - Rebuilds exon_hits_fullcontigs_clean.fasta correctly, preserving full IDs
#
# How to use:
#   1) Put this script in the same folder as:
#        - SequenceProteins.fasta
#        - one or more .pep ORF databases, e.g. *_ORFs_min12.pep
#   2) Make sure NCBI BLAST+ is installed.
#      Either add BLAST+ to your PATH, or set BLAST_BIN below / env var NCBI_BLAST_BIN.
#   3) In R: setwd("path/to/your/project/folder")
#   4) source("partA_DNAZOO_CNGB_github_ready.R")
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

query_fasta  <- file.path(data_dir, "SequenceProteins.fasta")
# query_fasta <- file.path(data_dir, "FullProteins.fasta")

merged_hits  <- file.path(data_dir, "blast_hits_merged.txt")
merged_fasta <- file.path(data_dir, "exon_hits_fullcontigs.fasta")
clean_fasta  <- file.path(data_dir, "exon_hits_fullcontigs_clean.fasta")

# If probe fails, automatically rebuild existing DB indices with -parse_seqids
# This is often useful for DNAZOO/CNGB-style headers.
auto_rebuild_db_on_probe_fail <- TRUE

norm_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

# ---- 2) Detect DB inputs (.pep) ----
detect_db_peps <- function(dirp) {
  parts <- list.files(dirp, pattern = "^part_\\d+[a-z]?\\.pep$", full.names = TRUE, ignore.case = TRUE)
  if (length(parts)) return(parts)

  orfs <- list.files(dirp, pattern = "_ORFs_min[0-9]+.*\\.pep$", full.names = TRUE, ignore.case = TRUE)
  if (length(orfs)) return(orfs)

  peps <- list.files(dirp, pattern = "\\.pep$", full.names = TRUE, ignore.case = TRUE)
  if (length(peps)) return(peps)

  character(0)
}

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

# ---- 3) Verify tools & inputs ----
cat("🔍 Checking executables and query file...\n")
cat("   blastp:      ", blastp, "\n", sep = "")
cat("   makeblastdb: ", makeblastdb, "\n", sep = "")
cat("   blastdbcmd:  ", blastdbcmd, "\n", sep = "")

for (exe in c(blastp, makeblastdb, blastdbcmd)) {
  if (!file.exists(exe)) stop("❌ Missing executable: ", exe, call. = FALSE)
}
if (!file.exists(query_fasta)) stop("❌ Missing query FASTA: ", query_fasta, call. = FALSE)
cat("✅ All executables and query verified.\n\n")

# ---- 4) DB indexing helpers ----
db_is_indexed <- function(db_prefix_path) {
  dirn <- dirname(db_prefix_path)
  base <- basename(db_prefix_path)
  esc <- gsub("([.\\+\\-\\^\\$\\(\\)\\[\\]\\{\\}\\|\\?\\*])", "\\\\\\1", base)
  pat_vol <- paste0("^", esc, "(\\.[0-9]{2})?\\.(pin|phr|psq)$")
  idx <- list.files(dirn, pattern = pat_vol, full.names = TRUE, ignore.case = TRUE)
  pal <- file.path(dirn, paste0(base, ".pal"))
  length(idx) > 0 || file.exists(pal)
}

run_cmd <- function(prog, args, label = NULL, stop_on_error = TRUE) {
  if (!is.null(label)) cat(label, "\n")
  out <- suppressWarnings(system2(prog, args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat("❌ Command failed (exit status ", status, ")\n", sep = "")
    cat("Program: ", prog, "\nArgs:\n  ", paste(args, collapse = " "), "\n\n", sep = "")
    cat("--- Output (stdout+stderr) ---\n")
    cat(paste(out, collapse = "\n"), "\n")
    if (stop_on_error) stop("Error: Stopping due to failed command above.", call. = FALSE)
  }
  invisible(out)
}

rebuild_db <- function(db_pep) {
  db_prefix <- sub("\\.pep$", "", db_pep, ignore.case = TRUE)
  cat("🔧 Rebuilding BLAST DB with -parse_seqids for:", basename(db_pep), "\n")
  run_cmd(
    makeblastdb,
    c("-in", shQuote(db_pep, type = "cmd"),
      "-dbtype", "prot",
      "-parse_seqids",
      "-out", shQuote(db_prefix, type = "cmd")),
    stop_on_error = TRUE
  )
  if (!db_is_indexed(db_prefix)) stop("❌ makeblastdb rebuild failed for: ", basename(db_pep), call. = FALSE)
  cat("✅ Rebuilt:", basename(db_pep), "\n\n")
}

# Initial build, only if missing
cat("🧱 Building / verifying BLAST DB indexes...\n")
for (db_pep in db_parts) {
  db_prefix <- sub("\\.pep$", "", db_pep, ignore.case = TRUE)
  if (db_is_indexed(db_prefix)) {
    cat("✅ DB already indexed:", basename(db_pep), "\n")
  } else {
    rebuild_db(db_pep)
  }
}
cat("\n")

# ---- 5) Split queries by length ----
suppressPackageStartupMessages(library(Biostrings))
qs <- readAAStringSet(query_fasta)
qlen <- width(qs)

short_idx <- which(qlen <= 60)
long_idx  <- which(qlen > 60)

q_short_fa <- file.path(tempdir(), "queries_short_le60.faa")
q_long_fa  <- file.path(tempdir(), "queries_long_gt60.faa")

if (length(short_idx)) writeXStringSet(qs[short_idx], q_short_fa, width = 60)
if (length(long_idx))  writeXStringSet(qs[long_idx],  q_long_fa,  width = 60)

cat("🧪 Query split:\n")
cat("   short (<=60 aa): ", length(short_idx), "\n", sep = "")
cat("   long  (>60 aa):  ", length(long_idx),  "\n\n", sep = "")

# ---- 6) Run BLASTP, 14 columns including stitle ----
cat("🚀 Running BLASTP across database(s)...\n")
if (file.exists(merged_hits)) file.remove(merged_hits)

outfmt_arg <- "6 qseqid sseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"

n_threads <- parallel::detectCores(logical = TRUE)
if (is.na(n_threads) || n_threads < 1) n_threads <- 4
cat("⚙️  Using BLAST threads:", n_threads, "\n\n")

for (db_pep in db_parts) {
  db_prefix <- sub("\\.pep$", "", db_pep, ignore.case = TRUE)

  if (length(long_idx)) {
    out_long <- paste0(db_prefix, "_hits_long.txt")
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
    run_cmd(blastp, blast_args_long, stop_on_error = FALSE)
    if (file.exists(out_long) && file.info(out_long)$size > 0) {
      file.append(merged_hits, out_long)
      cat("   ✅ added normal hits\n")
    } else cat("   ⚠️ no normal hits for this DB\n")
  }

  if (length(short_idx)) {
    out_short <- paste0(db_prefix, "_hits_short.txt")
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
    run_cmd(blastp, blast_args_short, stop_on_error = FALSE)
    if (file.exists(out_short) && file.info(out_short)$size > 0) {
      file.append(merged_hits, out_short)
      cat("   ✅ added short-exon rescue hits\n")
    } else cat("   ⚠️ no short-exon hits for this DB\n")
  }
}

cat("\n✅ All BLAST runs completed.\nMerged hits saved at:\n", merged_hits, "\n", sep = "")

# ---- 7) Read + deduplicate merged hits ----
if (!file.exists(merged_hits) || file.info(merged_hits)$size == 0)
  stop("❌ No BLAST hits found. Nothing to extract.", call. = FALSE)

hits <- read.table(merged_hits, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
if (ncol(hits) != 14) {
  stop("❌ BLAST merged hits is not 14 columns. It has ", ncol(hits),
       " columns. Check outfmt_arg and rerun.", call. = FALSE)
}
colnames(hits) <- c("qseqid","sseqid","sacc","stitle","pident","length","mismatch","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore")

# ID cleaner: first token up to whitespace; never trims underscores.
clean_db_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^>", "", x)

  tok <- sub("^([^[:space:]]+).*", "\\1", x)

  tok <- sub("^lcl\\|", "", tok)
  idx <- grepl("^gnl\\|", tok)
  if (any(idx, na.rm = TRUE)) tok[idx] <- sub("^.*\\|", "", tok[idx])

  tok
}

hits$sseqid_raw  <- hits$sseqid
hits$sacc_raw    <- hits$sacc
hits$stitle_raw  <- hits$stitle

hits$sseqid  <- clean_db_id(hits$sseqid_raw)
hits$sacc    <- clean_db_id(hits$sacc_raw)
hits$stitle_id <- clean_db_id(hits$stitle_raw)

hits <- unique(hits)

write.table(hits[, c("qseqid","sseqid","sacc","stitle","pident","length","mismatch","gapopen",
                     "qstart","qend","sstart","send","evalue","bitscore")],
            file = merged_hits, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("🧼 De-duplicated merged hits. Unique rows:", nrow(hits), "\n")
cat("🔎 Example IDs (raw -> cleaned):\n")
for (i in seq_len(min(5, nrow(hits)))) {
  cat("   sseqid:", hits$sseqid_raw[i], " -> ", hits$sseqid[i], "\n", sep = "")
}
cat("\n")

# ---- 8) Pick the right ID source by probing blastdbcmd ----
db_prefixes <- sub("\\.pep$", "", db_parts, ignore.case = TRUE)

probe_id_exists <- function(id) {
  id <- clean_db_id(id)
  if (is.na(id) || !nzchar(id)) return(FALSE)
  if (nchar(id) < 8L) return(FALSE)

  for (dbp in db_prefixes) {
    out <- suppressWarnings(system2(
      blastdbcmd,
      c("-db", shQuote(dbp, type = "cmd"),
        "-entry", shQuote(id, type = "cmd"),
        "-outfmt", "%i"),
      stdout = TRUE, stderr = TRUE
    ))
    status <- attr(out, "status")
    if (is.null(status) || status == 0) {
      if (length(out) && nzchar(out[1])) return(TRUE)
    }
  }
  FALSE
}

score_column <- function(vec, n_test = 30L) {
  u <- unique(clean_db_id(vec))
  u <- u[!is.na(u) & nzchar(u) & nchar(u) >= 8L]
  if (!length(u)) return(0L)
  u <- head(u, n_test)
  sum(vapply(u, probe_id_exists, logical(1)))
}

do_probe <- function() {
  cat("🧪 Probing which field matches blastdbcmd IDs...\n")
  sseqid_score <- score_column(hits$sseqid,    n_test = 30L)
  sacc_score   <- score_column(hits$sacc,      n_test = 30L)
  stitle_score <- score_column(hits$stitle_id, n_test = 30L)
  cat("   matches (sseqid):    ", sseqid_score, " / 30\n", sep = "")
  cat("   matches (sacc):      ", sacc_score,   " / 30\n", sep = "")
  cat("   matches (stitle_id): ", stitle_score, " / 30\n", sep = "")
  list(sseqid = sseqid_score, sacc = sacc_score, stitle_id = stitle_score)
}

scores <- do_probe()

if (max(unlist(scores)) == 0 && isTRUE(auto_rebuild_db_on_probe_fail)) {
  cat("\n⚠️  No probe matches. Auto-rebuilding DB(s) with -parse_seqids and retrying probe...\n\n")
  for (db_pep in db_parts) rebuild_db(db_pep)
  scores <- do_probe()
}

id_field_order <- c("stitle_id", "sseqid", "sacc")
best_field <- names(which.max(unlist(scores)))
if (!nzchar(best_field)) best_field <- "stitle_id"

cat("\n✅ Preferred ID field for extraction:", best_field, "\n\n")

# ---- 9) Extract sequences, trying best field first and falling back if needed ----
extract_with_field <- function(field) {
  ids <- unique(clean_db_id(hits[[field]]))
  ids <- ids[!is.na(ids) & nzchar(ids) & nchar(ids) >= 8L]
  if (!length(ids)) return(FALSE)

  id_file <- tempfile(fileext = ".txt")
  writeLines(ids, id_file, useBytes = TRUE)

  if (file.exists(merged_fasta)) file.remove(merged_fasta)
  any_extracted <- FALSE

  for (dbp in db_prefixes) {
    tmp_out <- tempfile(fileext = ".fasta")
    out <- suppressWarnings(system2(
      blastdbcmd,
      c("-db", shQuote(dbp, type = "cmd"),
        "-entry_batch", shQuote(id_file, type = "cmd"),
        "-out", shQuote(tmp_out, type = "cmd")),
      stdout = TRUE, stderr = TRUE
    ))
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) {
      cat("⚠️ blastdbcmd nonzero for DB:", basename(dbp), " (field=", field, ")\n", sep = "")
      next
    }
    if (file.exists(tmp_out) && file.info(tmp_out)$size > 0) {
      file.append(merged_fasta, tmp_out)
      any_extracted <- TRUE
    }
  }

  isTRUE(any_extracted) && file.exists(merged_fasta) && file.info(merged_fasta)$size > 0
}

cat("📦 Extracting hit sequences...\n")

try_fields <- unique(c(best_field, id_field_order))
ok <- FALSE
for (fld in try_fields) {
  cat("   → trying field:", fld, "\n")
  ok <- extract_with_field(fld)
  if (ok) {
    cat("   ✅ extraction succeeded using:", fld, "\n\n")
    best_field <- fld
    break
  } else {
    cat("   ⚠️ extraction produced nothing with:", fld, "\n")
  }
}

if (!ok) {
  stop(
    "❌ blastdbcmd produced no sequences using any candidate ID field (stitle_id/sseqid/sacc).\n",
    "This means BLAST is not reporting retrievable DB IDs, OR the DB headers are not parseable.\n\n",
    "NEXT CHECK:\n",
    "  1) Open your .pep and confirm the FIRST token after '>' is the full ID you want, e.g. HiC_scaffold_...\n",
    "  2) Ensure no earlier script rewrote the .pep headers into HiC_ / dup names.\n",
    call. = FALSE
  )
}

cat("🎉 All sequences extracted.\nRaw FASTA saved at:\n", merged_fasta, "\n", sep = "")

# ---- 10) Clean headers properly ----
gx_raw <- readAAStringSet(merged_fasta)

new_ids <- clean_db_id(names(gx_raw))
new_ids <- make.unique(new_ids, sep = "_dup")
names(gx_raw) <- new_ids

writeXStringSet(gx_raw, clean_fasta, width = 60)

cat("🧽 Clean FASTA written (fixed headers):\n", clean_fasta, "\n", sep = "")

cat("\n🔍 Clean FASTA header examples:\n")
cat(paste0("   - ", head(names(gx_raw), 6)), sep = "\n")
cat("\n")
cat("\n✅ Part A finished successfully.\n")
