#!/usr/bin/env Rscript
# pantera_dfam_export.R -- write a Dfam submission (Stockholm seed alignments)
# from a finished panteraGA run. Reads the run folder only; changes nothing.
#
# For every element that passes panteraGA's checks:
#   * the alignment in alignments/ is trimmed to the TE (the flanks svfind
#     adds are removed): the consensus is mapped onto the alignment columns
#     with mafft --add --keeplength --mapout and everything outside the first
#     and last mapped column is dropped. Inner columns are all kept.
#   * row names become Smitten V2 identifiers, accession:seqid:start-end_strand
#     (1-based, closed), recomputed for the trimmed piece of each copy.
#   * #=GC RF holds the panteraGA consensus on its columns. Dfam calls its own
#     consensus; if dfam-curator's `stk` is available it is run here
#     (stk edit --update-consensus, then stk lint).
#
# Requires: R >= 4.1, data.table, getopt, mafft (all in the pantera env).

suppressPackageStartupMessages({
  library(data.table)
  library(getopt)
})

spec <- matrix(c(
  "input",        "i", 1, "character", "panteraGA output folder [required]",
  "lib_name",     "b", 1, "character", "Library name (-b of panteraGA) [autodetected]",
  "output",       "o", 1, "character", "Output Stockholm file [<input>/<lib>-dfam.stk]",
  "authors",      "u", 1, "character", "Authors, 'First Last; First Last', optional ORCID:xxxx-xxxx-xxxx-xxxx prefix [required]",
  "taxa",         "t", 1, "character", "NCBI scientific name(s) for OC, ';'-separated [required]",
  "accessions",   "a", 1, "character", "TSV: genome file <tab> assembly accession [derived from GCA_/GCF_ file names]",
  "description",  "e", 1, "character", "DE template; {class} and {taxon} are replaced [\"{class} element identified by panteraGA in {taxon}\"]",
  "doi",          "d", 1, "character", "DOI to cite (RN/RD) [none]",
  "min_tsd_conf", "c", 1, "double",    "Min TSD_confidence to write TD [0.67]",
  "stk",          "k", 1, "character", "Path to dfam-curator's stk [stk in PATH, if any]",
  "no_stk",       "n", 0, "logical",   "Do not run stk even if found",
  "threads",      "T", 1, "integer",   "Threads for mafft [1]",
  "help",         "h", 0, "logical",   "This help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)
usage <- function(msg = NULL) {
  if (!is.null(msg)) cat("Error:", msg, "\n\n", file = stderr())
  cat(getopt(spec, usage = TRUE, command = "pantera_dfam_export.R"), file = stderr())
  quit(status = if (is.null(msg)) 0 else 1)
}
if (!is.null(opt$help)) usage()
for (req in c("input", "authors", "taxa")) if (is.null(opt[[req]])) usage(paste0("--", req, " is required"))
if (is.null(opt$min_tsd_conf)) opt$min_tsd_conf <- 0.67
if (is.null(opt$threads)) opt$threads <- 1L
if (is.null(opt$description)) opt$description <- "{class} element identified by panteraGA in {taxon}"
msg <- function(...) message("[dfam-export] ", ...)
die <- function(...) { message("[dfam-export] ERROR: ", ...); quit(status = 1) }

VAREXT <- 100L  # svfind -x, hard-coded in panteraGA

# Aligned FASTA -> data.table(header, seq)
read_afa <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- x[nzchar(x)]
  h <- startsWith(x, ">")
  if (!any(h)) stop("No FASTA records in ", file)
  grp <- cumsum(h)
  dt <- data.table(grp = grp[!h], line = x[!h])[, .(seq = paste(line, collapse = "")), by = grp]
  dt <- merge(data.table(grp = seq_len(sum(h)), header = sub("^>", "", x[h])), dt, by = "grp", all.x = TRUE)
  if (anyNA(dt$seq)) stop("Header without sequence in ", file)
  if (length(unique(nchar(dt$seq))) != 1) stop("Rows of different length in ", file)
  dt[, grp := NULL][]
}

parse_row_names <- function(h) {
  rev <- startsWith(h, "_R_")
  body <- sub("^_R_", "", h)
  m <- regmatches(body, regexec("^(.*)#(.+):([0-9]+)-([0-9]+)$", body))
  bad <- lengths(m) != 5
  if (any(bad)) stop("Unexpected row name(s): ", paste(h[bad], collapse = ", "))
  m <- do.call(rbind, m)
  data.table(genome = m[, 2], chr = m[, 3],
             start0 = as.numeric(m[, 4]), end0 = as.numeric(m[, 5]),
             strand = ifelse(rev, "-", "+"))
}

# Majority non-gap base per column (the subject the consensus is aligned to)
column_profile <- function(mat) {
  bases <- c("A", "C", "G", "T")
  counts <- sapply(bases, function(b) colSums(mat == b))
  if (is.null(dim(counts))) counts <- matrix(counts, nrow = 1)
  prof <- bases[max.col(counts, ties.method = "first")]
  prof[rowSums(counts) == 0] <- "N"
  paste(prof, collapse = "")
}

# For each consensus base, the alignment column it maps to (NA = the consensus
# has a base where the alignment has no column: an insertion). Uses
# mafft --add --keeplength --mapout, which profile-aligns the consensus to the
# existing alignment without changing its columns.
map_consensus <- function(cons, mat, mafft = Sys.which("mafft"), threads = 1L) {
  if (!nzchar(mafft)) stop("mafft not found in PATH")
  td <- tempfile("dfam_map"); dir.create(td); on.exit(unlink(td, recursive = TRUE))
  cf <- file.path(td, "cons.fa"); af <- file.path(td, "aln.fa")
  writeLines(c(">cons", cons), cf)
  writeLines(c(rbind(paste0(">s", seq_len(nrow(mat))),
                     apply(mat, 1, paste, collapse = ""))), af)
  st <- system2(mafft, c("--quiet", "--thread", threads, "--add", shQuote(cf),
                         "--keeplength", "--mapout", shQuote(af)),
                stdout = FALSE, stderr = FALSE)
  mf <- paste0(cf, ".map")
  if (st != 0 || !file.exists(mf)) stop("mafft --mapout failed")
  m <- fread(mf, skip = 2, header = FALSE, sep = ",", strip.white = TRUE,
             colClasses = "character")
  if (nrow(m) != nchar(cons)) stop("mafft map has ", nrow(m), " rows for a ",
                                   nchar(cons), " bp consensus")
  col <- suppressWarnings(as.integer(m[[3]]))   # "-" -> NA
  n_unplaced <- sum(is.na(col))
  col <- repair_map(col, cons, mat)
  prof <- strsplit(column_profile(mat), "")[[1]]
  cb <- strsplit(cons, "")[[1]]
  ok <- !is.na(col)
  list(col = col, identity = 100 * mean(cb[ok] == prof[col[ok]]),
       n_repaired = n_unplaced - sum(is.na(col)))
}

# mafft --add can leave a consensus base unplaced inside a homopolymer even
# when a free column is available. For each unplaced base, take the stretch
# between the nearest placed neighbours (widened to the whole homopolymer) and
# re-place its bases in order on the columns in between, maximising the number
# of rows that carry the same base. Bases that still have no room stay NA and
# get an all-gap column later.
repair_map <- function(col, cons, mat) {
  cb <- strsplit(cons, "")[[1]]; n <- length(cb)
  for (pass in 1:3) {
    na <- which(is.na(col)); if (!length(na)) break
    i <- na[1]
    j1 <- i; while (j1 > 1 && (is.na(col[j1 - 1]) || cb[j1 - 1] == cb[i])) j1 <- j1 - 1
    j2 <- i; while (j2 < n && (is.na(col[j2 + 1]) || cb[j2 + 1] == cb[i])) j2 <- j2 + 1
    a <- if (j1 > 1) col[j1 - 1] else 0L
    z <- if (j2 < n) col[j2 + 1] else ncol(mat) + 1L
    cols <- seq_len(z - a - 1L) + a
    k <- j2 - j1 + 1L; m <- length(cols)
    if (m < k) { col[j1:j2][is.na(col[j1:j2])] <- -1L; next }  # no room
    # score[b, c] = rows with consensus base b at column c
    sc <- sapply(cols, function(c) colSums(outer(mat[, c], cb[j1:j2], "==")))
    sc <- matrix(sc, nrow = k)
    # order-preserving assignment of k bases to m columns (DP)
    D <- matrix(-Inf, k + 1, m + 1); D[1, ] <- 0
    for (x in 1:k) for (y in x:m)
      D[x + 1, y + 1] <- max(D[x + 1, y], D[x, y] + sc[x, y])
    y <- m; pick <- integer(k)
    for (x in k:1) {
      while (D[x + 1, y + 1] == D[x + 1, y]) y <- y - 1L
      pick[x] <- cols[y]; y <- y - 1L
    }
    col[j1:j2] <- pick
  }
  col[!is.na(col) & col < 0] <- NA_integer_
  col
}

#' Trim an alignment to the consensus span.
#' Returns the trimmed matrix (upper case, "." gaps), the RF line (consensus
#' bases on their columns, "." elsewhere), the Smitten names, and diagnostics.
trim_to_consensus <- function(afa, cons, threads = 1L) {
  cons <- toupper(cons)
  rows <- parse_row_names(afa$header)
  mat <- toupper(do.call(rbind, strsplit(afa$seq, "", fixed = TRUE)))
  mat[mat %in% c("-", ".", "_", "~")] <- "-"

  mp <- map_consensus(cons, mat, threads = threads)
  col <- mp$col
  cbase <- strsplit(cons, "")[[1]]

  # Insertions: give each unmapped consensus base its own all-gap column,
  # placed right after the previous mapped column.
  n_ins <- sum(is.na(col))
  if (n_ins > 0) {
    newmat <- NULL; newcol <- integer(length(col)); prev <- 0L; out <- 0L
    for (k in seq_along(col)) {
      if (!is.na(col[k])) {
        take <- (prev + 1L):col[k]
        newmat <- cbind(newmat, mat[, take, drop = FALSE])
        out <- out + length(take); prev <- col[k]
      } else {
        newmat <- cbind(newmat, matrix("-", nrow(mat), 1)); out <- out + 1L
      }
      newcol[k] <- out
    }
    if (prev < ncol(mat)) newmat <- cbind(newmat, mat[, (prev + 1L):ncol(mat), drop = FALSE])
    mat <- newmat; col <- newcol
  }

  c1 <- col[1]; c2 <- col[length(col)]
  if (c1 > c2) stop("Consensus maps in reverse order; check orientation")

  # Bases kept per row, in row order
  isb <- mat != "-"
  cum <- t(apply(isb, 1, cumsum))
  if (nrow(mat) == 1) cum <- matrix(cum, nrow = 1)
  L  <- cum[, ncol(mat)]
  q1 <- (if (c1 > 1) cum[, c1 - 1] else 0) + 1   # first kept base
  q2 <- cum[, c2]                                 # last kept base

  # Flank on the genomic left of the segment. 200 in total means 100 + 100;
  # anything else means a contig edge, which the .maf alone cannot resolve.
  seglen <- rows$end0 - rows$start0
  flank_ok <- (L - seglen) == 2 * VAREXT
  lx <- VAREXT

  fwd <- rows$strand == "+"
  gstart <- ifelse(fwd, rows$start0 - lx + q1, rows$start0 - lx + L - q2 + 1)
  gend   <- ifelse(fwd, rows$start0 - lx + q2, rows$start0 - lx + L - q1 + 1)

  keep <- (q2 >= q1) & flank_ok
  tmat <- mat[keep, c1:c2, drop = FALSE]
  tmat[tmat == "-"] <- "."

  rf <- rep(".", ncol(mat)); rf[col] <- cbase
  rf <- paste(rf[c1:c2], collapse = "")

  info <- data.table(row = afa$header, kept = keep, flank_ok = flank_ok,
                     n_bases = pmax(q2 - q1 + 1, 0),
                     # offsets of the consensus edges from the svfind segment
                     # edges, in row orientation (0 = exactly at the breakpoint)
                     left_off = q1 - (VAREXT + 1), right_off = q2 - (L - VAREXT))
  list(mat = tmat, rf = rf,
       genome = rows$genome[keep], chr = rows$chr[keep],
       start = gstart[keep], end = gend[keep], strand = rows$strand[keep],
       cols = c(c1, c2), aln_width = ncol(mat), n_insert_cols = n_ins,
       map_identity = mp$identity, n_repaired = mp$n_repaired, rows = info)
}

# ---------------------------------------------------------------- inputs ---
inp <- normalizePath(opt$input, mustWork = FALSE)
if (!dir.exists(inp)) die("input folder not found: ", inp)
if (is.null(opt$lib_name)) {
  fas <- list.files(inp, "-pantera-final\\.fa$")
  if (length(fas) != 1) die("cannot autodetect the library name (", length(fas), " *-pantera-final.fa files); use -b")
  opt$lib_name <- sub("-pantera-final\\.fa$", "", fas)
}
lib <- opt$lib_name
fa_file <- file.path(inp, paste0(lib, "-pantera-final.fa"))
st_file <- list.files(inp, paste0("^", gsub("([.+])", "\\\\\\1", lib), "-pantera-final[._]stats\\.tsv$"), full.names = TRUE)
if (!file.exists(fa_file)) die("not found: ", fa_file)
if (length(st_file) != 1) die("stats table not found for library '", lib, "'")
aln_dir <- file.path(inp, "alignments")
if (!dir.exists(aln_dir)) die("alignments/ folder not found (run panteraGA keeping the alignments)")
if (!nzchar(Sys.which("mafft"))) die("mafft not found in PATH")
out_file <- if (is.null(opt$output)) file.path(inp, paste0(lib, "-dfam.stk")) else opt$output

# consensus sequences
read_fasta <- function(f) {
  x <- readLines(f, warn = FALSE); h <- startsWith(x, ">")
  s <- tapply(x[!h], cumsum(h)[!h], paste, collapse = "")
  setNames(toupper(as.character(s)), sub("^>", "", x[h])[as.integer(names(s))])
}
cons <- read_fasta(fa_file)

stats <- fread(st_file, sep = "\t", colClasses = list(character = c("name", "TSD_motif", "homology_class")))
stats[, name := sub("^>", "", name)]
stats[, pass := as.logical(pass)]
sel <- stats[pass == TRUE]
msg(nrow(stats), " elements in ", basename(st_file), ", ", nrow(sel), " pass")
if (!nrow(sel)) die("no element passes")

# authors, taxa, DE
authors <- trimws(strsplit(opt$authors, ";")[[1]]); authors <- authors[nzchar(authors)]
if (any(grepl(",", authors))) die("separate authors with ';', not ','")
bad_au <- authors[!grepl("^(ORCID:[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X] )?[^ .]{2,}( [^ ]+)+$", authors)]
if (length(bad_au)) die("authors must be 'First Last' with the given name spelled out: ", paste(bad_au, collapse = "; "))
taxa <- trimws(strsplit(opt$taxa, ";")[[1]]); taxa <- taxa[nzchar(taxa)]

# genome file -> accession
acc_tab <- if (!is.null(opt$accessions)) {
  a <- fread(opt$accessions, header = FALSE, sep = "\t", col.names = c("genome", "acc"))
  setNames(a$acc, a$genome)
} else character(0)
accession_for <- function(g) {
  if (g %in% names(acc_tab)) return(acc_tab[[g]])
  if (basename(g) %in% names(acc_tab)) return(acc_tab[[basename(g)]])
  m <- regmatches(basename(g), regexec("^(GC[AF]_[0-9]+\\.[0-9]+)", basename(g)))[[1]]
  if (length(m) == 2) m[2] else NA_character_
}

version <- {
  lg <- file.path(inp, "pantera.log")
  v <- if (file.exists(lg)) regmatches(readLines(lg, n = 50, warn = FALSE),
                                        regexpr("panteraGA [0-9][0-9.]*", readLines(lg, n = 50, warn = FALSE))) else character(0)
  if (length(v)) v[1] else "panteraGA"
}

# alignment file for an element: make.names() as panteraGA does, else by
# separator-insensitive match
aln_files <- list.files(aln_dir, "\\.maf$", full.names = TRUE)
norm <- function(x) gsub("[^A-Za-z0-9]", "_", x)
aln_for <- function(nm) {
  f <- file.path(aln_dir, paste0(make.names(nm), ".maf"))
  if (file.exists(f)) return(f)
  hit <- aln_files[norm(sub("\\.maf$", "", basename(aln_files))) == norm(nm)]
  if (length(hit) == 1) hit else NA_character_
}

iupac_ok <- function(s) grepl("^[ACGTRYSWKMBDHVN]+$", s)
genome_cols <- if ("Ns" %in% names(stats)) names(stats)[(match("Ns", names(stats)) + 1):ncol(stats)] else character(0)

# --------------------------------------------------------------- records ---
records <- character(0); report <- list(); missing_acc <- character(0)
for (k in seq_len(nrow(sel))) {
  s <- sel[k]; nm <- s$name
  if (is.na(cons[nm])) { msg("WARNING ", nm, ": not in ", basename(fa_file), ", skipped"); next }
  f <- aln_for(nm)
  if (is.na(f)) { msg("WARNING ", nm, ": no alignment, skipped"); report[[nm]] <- data.table(name = nm, status = "no alignment"); next }
  r <- tryCatch(trim_to_consensus(read_afa(f), cons[[nm]], threads = opt$threads),
                error = function(e) { msg("WARNING ", nm, ": ", conditionMessage(e), ", skipped"); NULL })
  if (is.null(r)) { report[[nm]] <- data.table(name = nm, status = "error"); next }
  if (!length(r$chr)) { msg("WARNING ", nm, ": no copy left after trimming, skipped"); next }

  acc <- vapply(r$genome, accession_for, "")
  missing_acc <- union(missing_acc, r$genome[is.na(acc)])
  ids <- sprintf("%s:%s:%.0f-%.0f_%s", acc, r$chr, r$start, r$end, r$strand)
  if (anyDuplicated(ids)) { keep <- !duplicated(ids); ids <- ids[keep]; r$mat <- r$mat[keep, , drop = FALSE] }

  cls <- sub("^[^#]*#", "", nm)
  if (!grepl("/", cls) && cls %in% c("PLE")) msg("WARNING ", nm, ": class '", cls, "' has no Dfam equivalent without a subtype")
  de <- gsub("{class}", cls, gsub("{taxon}", taxa[1], opt$description, fixed = TRUE), fixed = TRUE)
  if (nchar(de) > 80) de <- sub(" in .*$", "", de)
  if (nchar(de) > 80) die("DE longer than 80 characters: ", de)

  cc <- sprintf("panteraGA element %s, from a comparison of %d genomes.", nm, max(1, length(genome_cols)))
  if (!is.na(s$homology_class) && s$homology_class != "none" && s$homology_class != cls)
    cc <- c(cc, sprintf("Closest protein homology: %s (probability %s).", s$homology_class, s$homology_prob))
  if (length(genome_cols)) {
    cp <- unlist(s[, ..genome_cols])
    lab <- vapply(genome_cols, function(g) { a <- accession_for(g); if (is.na(a)) g else a }, "")
    cc <- c(cc, paste0("Copies per genome: ", paste(lab, cp, sep = " ", collapse = ", "), "."))
  }

  gf <- function(tag, val) sprintf("#=GF %-4s  %s", tag, val)
  rec <- c("# STOCKHOLM 1.0",
           gf("DE", de),
           gf("AU", paste(authors, collapse = "; ")),
           gf("TP", cls),
           gf("OC", taxa),
           gf("SQ", length(ids)),
           gf("SE", version),
           gf("BM", paste(version, "(https://github.com/piosierra/panteraGA); flanks trimmed at the consensus ends")))
  tsd <- toupper(s$TSD_motif)
  if (!is.na(tsd) && nzchar(tsd) && iupac_ok(tsd) && !is.na(s$TSD_confidence) && s$TSD_confidence >= opt$min_tsd_conf)
    rec <- c(rec, gf("TD", tsd))
  if (!is.null(opt$doi)) rec <- c(rec, gf("RN", "[1]"), gf("RD", opt$doi))
  rec <- c(rec, gf("CC", cc))
  w <- max(nchar(c(ids, "#=GC RF"))) + 2
  rec <- c(rec, sprintf("%-*s%s", w, "#=GC RF", r$rf),
           sprintf("%-*s%s", w, ids, apply(r$mat, 1, paste, collapse = "")), "//")
  records <- c(records, rec)
  report[[nm]] <- data.table(name = nm, status = "ok", copies = length(ids),
                             dropped = nrow(r$rows) - length(ids), width = nchar(r$rf),
                             trimmed_cols = r$aln_width - nchar(r$rf))
}

if (length(missing_acc))
  die("no assembly accession for: ", paste(missing_acc, collapse = ", "),
      "\n  Give a TSV with -a (genome file <tab> accession, e.g. GCA_000002195.1).")
if (!length(records)) die("nothing to write")

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
writeLines(records, out_file)
rep <- rbindlist(report, fill = TRUE)
print(rep, row.names = FALSE)
msg("wrote ", sum(rep$status == "ok"), " families to ", out_file)

# ------------------------------------------------- dfam-curator (optional) ---
stk <- if (!is.null(opt$stk)) opt$stk else Sys.which("stk")
if (is.null(opt$no_stk) && nzchar(stk)) {
  tmp <- paste0(out_file, ".tmp")
  if (system2(stk, c("edit", "--update-consensus", "-o", shQuote(tmp), shQuote(out_file))) == 0) {
    file.rename(tmp, out_file)
    msg("#=GC RF replaced by Dfam's consensus (stk edit --update-consensus)")
  } else msg("WARNING stk edit failed; RF keeps the panteraGA consensus")
  system2(stk, c("lint", "--min-severity", "warn", shQuote(out_file)))
} else {
  msg("stk (dfam-curator) not run: RF holds the panteraGA consensus. Before submitting, run\n",
      "  stk edit --update-consensus -o ", out_file, ".tmp ", out_file, " && mv ", out_file, ".tmp ", out_file, "\n",
      "  stk lint --genome <assembly> ", out_file)
}
