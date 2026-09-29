# =============================================================================
# maf_to_dfam_stockholm.R
#
# Converts an aligned-FASTA file (headers like:
#   >GCA_964199755.2_bLarMic1.hap1.2_genomic.fna#OZ118746.2:159792819-159797672
#   >_R_GCA_964199725.2_bLarMic1.hap2.2_genomic.fna#OZ287837.1:11165779-11170623
# ) into a Dfam-compliant Stockholm 1.0 "seed" file, following the required
# fields described in:
#   https://github.com/Dfam-consortium/dfam-curator/blob/main/Dfam_Seeds.md
#
# Required #=GF fields written: DE, AU, TP, OC, SQ
# Required per-column annotation written: #=GC RF (consensus)
#
# NOTE: the input is *not* UCSC/multiz MAF format (no "a score" / "s" lines) -
# it is an aligned FASTA where each record's header encodes assembly
# accession, sequence id, and 1-based closed coordinates, exactly like the
# example given. Sequence identifiers are converted to the Dfam "Smitten"
# format: ASSEMBLY:SEQID:START-END_STRAND
#
# A leading "_R_" in the header (as added by `mafft --adjustdirection`)
# is interpreted as "this copy was reverse-complemented to align it" and
# is encoded as strand "-" in the Smitten identifier.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- IUPAC ambiguity lookup for consensus calling --------------------------
.iupac_lookup <- c(
  "A" = "A", "C" = "C", "G" = "G", "T" = "T",
  "AG" = "R", "CT" = "Y", "CG" = "S", "AT" = "W", "GT" = "K", "AC" = "M",
  "CGT" = "B", "AGT" = "D", "ACT" = "H", "ACG" = "V",
  "ACGT" = "N"
)

.call_consensus_column <- function(countA, countC, countG, countT, ambig_threshold) {
  counts <- c(A = countA, C = countC, G = countG, T = countT)
  total <- sum(counts)
  if (total == 0L) return(".")
  maxc <- max(counts)
  selected <- names(counts)[counts > 0 & counts >= ambig_threshold * maxc]
  key <- paste0(sort(selected), collapse = "")
  out <- .iupac_lookup[[key]]
  if (is.null(out)) "N" else out
}

# ---- Header parser -----------------------------------------------------
# Extracts assembly/source label, sequence id, start, end and strand from a
# header of the form:
#   [_R_]<assembly-or-label>[_<free text>]#<seqid>:<start>-<end>
#
# Two label styles are supported:
#   1. NCBI-style accessions, optionally followed by free text that is
#      discarded:      GCA_964199755.2_bLarMic1.hap1.2_genomic.fna#...
#   2. Any other bare label with no free text to strip:
#      Dana#scaffold_51:14762-18274
.parse_header <- function(header) {
  strand <- "+"
  body <- header
  if (startsWith(body, "_R_")) {
    strand <- "-"
    body <- sub("^_R_", "", body)
  }
  
  hash_pos <- regexpr("#", body, fixed = TRUE)
  if (hash_pos < 0) {
    stop("Header does not match expected format (no '#' found): '", header, "'")
  }
  before <- substr(body, 1, hash_pos - 1)
  after  <- substr(body, hash_pos + 1, nchar(body))
  
  # If `before` looks like an NCBI accession with trailing free text,
  # keep only the accession; otherwise use `before` verbatim as the label.
  acc_pattern <- "^(GC[AF]_[0-9]+\\.[0-9]+)(?:_.*)?$"
  m_acc <- regmatches(before, regexec(acc_pattern, before))[[1]]
  assembly <- if (length(m_acc) == 2) m_acc[2] else before
  
  coord_pattern <- "^([^:]+):([0-9]+)-([0-9]+)$"
  m_coord <- regmatches(after, regexec(coord_pattern, after))[[1]]
  if (length(m_coord) != 4) {
    stop("Header does not match expected format (bad seqid:start-end part): '", header, "'")
  }
  
  list(
    strand   = strand,
    assembly = assembly,
    seqid    = m_coord[2],
    start    = m_coord[3],
    end      = m_coord[4]
  )
}

#' Convert an aligned-FASTA ("MAF") file to a Dfam Stockholm seed file
#'
#' @param maf_file        Path to the input aligned-FASTA file.
#' @param output_dir      Directory the .stk file will be written into
#'                         (created if it doesn't exist).
#' @param de              #=GF DE  Short description of the family (<=80 chars).
#' @param au               #=GF AU  Author name(s), "First Last". A character
#'                         vector produces one AU line per element; a single
#'                         string with ';'-separated names is also accepted.
#' @param tp              #=GF TP  RepeatMasker classification path (or shorthand
#'                         alias, e.g. "LTR/ERV1").
#' @param oc              #=GF OC  Taxonomic scope. A character vector produces
#'                         one OC line per element (NCBI scientific names).
#' @param output_filename Optional output file name (default: input basename
#'                         with a .stk extension).
#' @param gap_char        Gap character to use in the output alignment
#'                         (Dfam convention: "."). Any of "-", "_", "~", "."
#'                         found in the input are normalized to this.
#' @param ambig_threshold Fraction of the top column count (0-1) at which an
#'                         additional base is folded into an IUPAC ambiguity
#'                         code for the RF consensus line. 1 = strict
#'                         majority-only consensus (ties only); lower values
#'                         admit more ambiguity codes. Default 0.5.
#' @param overwrite       Overwrite an existing output file? Default FALSE.
#'
#' @return (invisibly) the path to the Stockholm file that was written.
maf_to_dfam_stockholm <- function(maf_file,
                                  output_dir,
                                  de,
                                  au,
                                  tp,
                                  oc,
                                  output_filename = NULL,
                                  gap_char = ".",
                                  ambig_threshold = 0.5,
                                  overwrite = FALSE) {
  
  stopifnot(file.exists(maf_file))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (nchar(de) > 80) stop("`de` must be <= 80 characters (Dfam DE field limit).")
  
  # ---- 1. Read and group lines into (header, sequence) records via data.table
  raw <- readLines(maf_file, warn = FALSE)
  raw <- raw[nzchar(raw)]  # drop blank lines (block format is not used here)
  
  dt <- data.table(line = raw)
  dt[, is_header := startsWith(line, ">")]
  dt[, grp := cumsum(is_header)]
  
  headers <- dt[is_header == TRUE, .(grp, header = sub("^>", "", line))]
  seqs    <- dt[is_header == FALSE, .(seq = paste0(line, collapse = "")), by = grp]
  
  records <- merge(headers, seqs, by = "grp", all.x = TRUE, sort = TRUE)
  setorder(records, grp)
  
  if (nrow(records) == 0) stop("No FASTA records found in ", maf_file)
  if (any(is.na(records$seq)) || any(!nzchar(records$seq))) {
    stop("One or more headers has no associated sequence.")
  }
  
  # ---- 2. Parse headers into Smitten identifiers -----------------------
  parsed <- rbindlist(lapply(records$header, .parse_header))
  records <- cbind(records, parsed)
  records[, smitten_name := paste0(assembly, ":", seqid, ":", start, "-", end, "_", strand)]
  
  if (anyDuplicated(records$smitten_name)) {
    stop("Duplicate sequence identifiers after conversion: ",
         paste(records$smitten_name[duplicated(records$smitten_name)], collapse = ", "))
  }
  
  # ---- 3. Validate equal alignment length & build the alignment matrix --
  aln_len <- unique(nchar(records$seq))
  if (length(aln_len) != 1) {
    stop("Sequences are not all the same length (unaligned/block-format input?).")
  }
  
  seq_mat <- do.call(rbind, strsplit(records$seq, "", fixed = TRUE))
  rownames(seq_mat) <- records$smitten_name
  
  # normalize all recognized gap characters to `gap_char`
  seq_mat[seq_mat %in% c("-", "_", "~", ".")] <- gap_char
  
  # ---- 4. Compute the #=GC RF consensus line (vectorized column tallies) --
  upper_mat <- toupper(seq_mat)
  countA <- colSums(upper_mat == "A")
  countC <- colSums(upper_mat == "C")
  countG <- colSums(upper_mat == "G")
  countT <- colSums(upper_mat == "T")
  
  rf <- mapply(.call_consensus_column, countA, countC, countG, countT,
               MoreArgs = list(ambig_threshold = ambig_threshold))
  rf_line <- paste0(rf, collapse = "")
  
  # ---- 5. Assemble output lines -----------------------------------------
  seq_rows <- vapply(seq_len(nrow(seq_mat)),
                     function(i) paste0(seq_mat[i, ], collapse = ""),
                     character(1))
  names(seq_rows) <- rownames(seq_mat)
  
  au <- unlist(strsplit(au, ";\\s*(?=[A-Z])", perl = TRUE))
  au <- trimws(au)
  
  label_width <- max(nchar(c(names(seq_rows), "#=GC RF"))) + 4
  
  gf_lines <- c(
    sprintf("#=GF DE%s%s", strrep(" ", max(4, label_width - 7)), de),
    sprintf("#=GF AU%s%s", strrep(" ", max(4, label_width - 7)), au),
    sprintf("#=GF TP%s%s", strrep(" ", max(4, label_width - 7)), tp),
    sprintf("#=GF OC%s%s", strrep(" ", max(4, label_width - 7)), oc),
    sprintf("#=GF SQ%s%d", strrep(" ", max(4, label_width - 7)), nrow(seq_mat))
  )
  
  rf_out <- sprintf("%-*s%s", label_width, "#=GC RF", rf_line)
  seq_out <- sprintf("%-*s%s", label_width, names(seq_rows), seq_rows)
  
  out_lines <- c("# STOCKHOLM 1.0", gf_lines, rf_out, seq_out, "//")
  
  # ---- 6. Write the file --------------------------------------------------
  if (is.null(output_filename)) {
    output_filename <- paste0(tools::file_path_sans_ext(basename(maf_file)), ".stk")
  }
  out_path <- file.path(output_dir, output_filename)
  if (file.exists(out_path) && !overwrite) {
    stop("Output file already exists (set overwrite = TRUE to replace it): ", out_path)
  }
  
  writeLines(out_lines, con = out_path)
  message("Wrote Dfam Stockholm seed file: ", out_path)
  invisible(out_path)
}

# =============================================================================
# Example usage:
#
# maf_to_dfam_stockholm(
#   maf_file    = "family1.maf",
#   output_dir  = "stockholm_out",
#   de          = "Repetitive element family identified in Lariidae",
#   au          = "Jane Doe",
#   tp          = "Interspersed_Repeat;Unknown",
#   oc          = c("Larus michahellis", "Laridae")
# )
# =============================================================================