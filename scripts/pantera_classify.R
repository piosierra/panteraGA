#!/usr/bin/env Rscript
### pantera_classify.R: homology classifier of TE sequences with the pantera HMM libraries.
###
### R port of the unified caller of the panteraGA HMM classifier project (step 12: classify_library.sh,
### 1202-1204, and the call rules of steps 5, 7, 8 and 9 it uses). Same output, byte for byte.
### Loaded by panteraGA with source(), or run alone:
###   Rscript pantera_classify.R -i library.fa -o outdir -L libdir [-T threads] [-C chunk] [-k]
###
### Pieces, per chunk of sequences (one process per chunk, one thread per tool). Results do not depend on the
### number of threads or the chunk size: fixed -Z in hmmsearch and bit-score thresholds everywhere.
###   protein : miniprot (one run per sequence) of the reference proteins plus six-frame ORFs; hmmsearch --cut_ga
###             of the protein domain models; element rule (label = best sum of domain classes)
###   tir     : nhmmer of the TIR end models on the end windows (first and last 150 bp)
###   hel     : nhmmer of the Helitron end models on the end windows
###   sine    : nhmmer of the SINE models; cmscan of the Rfam-derived head models and a cmalign check that the
###             head is not a whole RNA gene
###   flags   : TRF (tandem), dustmasker (low complexity), cmscan of 42 eukaryotic ncRNA families, DIAMOND
###             against Swiss-Prot host proteins (only when libdir/flags/sp_euk_host.dmnd exists)
### Then, over the whole library: end partners (BLASTN of the ends of sequences without a protein call against
### the ends of autonomous DNA transposons of the same library and of a curated set), and the final call.
###
### Output: outdir/classification.tsv (one row per sequence) and outdir/summary.tsv. Columns: see hc_cols below
### and the README.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# Settings. Values are those of the step-12 caller (tables/settings.tsv and the defaults of the step scripts).
hmm_cfg <- list(
  chunk = 500L,            # max sequences per chunk (one process each); results do not depend on the chunks
  chunk_min = 50L,         # min sequences per chunk by default: each chunk loads the models again (the TIR end
                           #   search takes about 14 s to start), so smaller chunks waste time
  min_aa = 80L,            # ORFs: min length (codons) of a stop-to-stop stretch kept as a protein
  mp_opts = "-N 100 -p 0.3 --outs=0.3",  # miniprot options
  mp_overlap = 0.2,        # miniprot: an alignment overlapping a kept one by more than this fraction of the
                           #   shorter (+ 1 bp) is dropped (one protein per locus)
  orf_in_mp = 0.5,         # ORFs: dropped when this fraction or more lies in a kept miniprot locus (same strand)
  hmm_z = 69302L,          # hmmsearch -Z and --domZ (fixed: scores do not depend on the input size)
  dom_min_cov = 0.2,       # protein rule: a domain hit counts if it covers at least this fraction of its model
  dom_exempt = 50,         #   or scores at least this many bits
  lin_margin = 20,         # sub label: min lead over the next sub label (bits)
  low_score = 100,         # protein calls on one domain class below this score are 'low'
  host_like = c("RH", "YR", "CHD", "HEL"),  # protein calls on one of these classes alone are 'low'
  end_win = 150L,          # TIR, Helitron and partner end windows (bp)
  end_hit = 10,            # TIR and Helitron ends: min bits of a hit (also nhmmer -T)
  end_anchor = 30L,        #   an end hit must lie within this many bp of the window end
  tir_sum = 30,            # TIR call: both ends with sL + sR >= this,
  tir_single = 60,         #   or one end alone >= this
  tir_high = 40,           #   'high' when both ends and sL + sR >= this
  tir_lead = 10,           #   superfamily: min lead over the second superfamily (else 'DNA')
  tir_sub_lead = 10,       #   sub label: min lead over the second, and within this of the superfamily score
  hel_sum = 30,            # Helitron call: both ends with s5 + s3 >= this,
  hel_single = 80,         #   or one end alone >= this
  hel_high = 40,           #   'high' when both ends and s5 + s3 >= this
  hel_glead = 10,          #   group (Helitron-1 / Helentron): min lead over the other (else 'RC')
  sine_hit = 10,           # SINE: min bits of a hit (also nhmmer -T)
  sine_head_window = 60L,  #   a HEAD hit must start within this many bp of the element's 5' end
  sine_min_total = 18,     #   min HEAD + BODY score of the called label
  sine_lead = 10,          #   min lead over the next label (else 'ambiguous')
  sine_max_len = 800L,     #   longer elements get 'SINE_part' instead of a SINE label
  sine_high = 25,          #   'high' with a HEAD hit and a score >= this
  sine_rna_frac = 0.8,     #   'RNA_gene-like' when the head fills >= this fraction of its Rfam model
  rna_min = 12,            # RNA check: min score of the Rfam head hit (also cmscan -T),
  rna_span = 0.8,          #   min fraction of the model spanned by the hit,
  rna_rest = 20L,          #   and max bp of the element after the hit
  trf_opts = "2 7 7 80 10 50 500 -h -ngs",  # TRF parameters
  host_evalue = "1e-10",   # DIAMOND e-value for host proteins
  partner_win = 150L,      # end partners: end windows (bp)
  partner_anchor = 20L,    #   a hit must lie within this many bp of the query and partner termini
  partner_lead = 10,       #   the partner's label must lead the best other label by this many bits
  partner_min_len = 50L,   #   min length of a query
  partner_block = 2000L,   #   queries per BLASTN run
  partner_evalue = "1e-5", #   BLASTN e-value
  # Final call (step 12, tables/settings.tsv)
  weak_single = c("LTR:INT", "LTR:RH", "LTR:RT:low", "DIRS:YR", "PLE:GIY"),
                           # protein calls whose only domain class is this one (TYPE:class; ':low' = only when
                           #   the protein confidence is low) get confidence 'weak'
  alt = c("LTR:INT" = "DNA/Ginger-1,DNA/Ginger-2,DNA/IS481EU,POLINTON",  # labels sharing that domain
          "LTR:RH" = "LINE",
          "LTR:RT" = "DIRS/Viper,DIRS/TATE",
          "DIRS:YR" = "CRYPTON",
          "PLE:GIY" = "other GIY-YIG nucleases"),
  close_lead = 50,         # note close_runner_up:<label> when the call leads its runner-up by less (bits)
  chimera_runner_bits = 100, # chimera: runner-up label of another TYPE with at least this score
  chimera_runner_frac = 0.5, #   and at least this fraction of the call score
  tandem_frac = 0.5,       # flag tandem: TRF arrays cover at least this fraction of the sequence
  lowcomp_frac = 0.3,      # flag low_complexity: dustmasker covers at least this fraction
  ncrna_frac = 0.3,        # flag ncRNA: Rfam hits cover at least this fraction (or any rRNA hit)
  rrna = c("SSU_rRNA_eukarya", "LSU_rRNA_eukarya", "5_8S_rRNA", "SSU_rRNA_microsporidia"),
  host_pid = 50,           # flag host_protein: a protein without a pantera hit matches a Swiss-Prot host
  host_cov = 0.5,          #   protein at >= host_pid % identity over >= host_cov of the Swiss-Prot protein,
  host_aln = 100           #   and the alignment is >= host_aln aa or the sequence has no high protein call
)

# Columns of classification.tsv
hc_cols <- c("input_id", "id", "length", "label", "confidence", "evidence", "sub", "sub_source", "sub_models",
             "domains", "fits_type", "prot_score", "prot_lead", "runner_up", "tir", "partner", "helitron", "sine",
             "alternatives", "flags", "notes")


### Small helpers

# Log message: panteraGA's lx() when loaded by panteraGA, otherwise stderr
hc_log <- function(...) {
  m <- paste(...)
  if (exists("lx", mode = "function")) lx(m) else message(format(Sys.time(), "%H:%M:%S "), m)
}

# Runs a shell command; stops naming the command if it fails (ok_status: exit codes that are not failures)
hc_sys <- function(cmd, ok_status = 0L) {
  st <- system(cmd)
  if (!(st %in% ok_status)) stop("command failed (exit ", st, "): ", cmd)
  invisible(st)
}

# Formats numbers as Python's f'{x:.Nf}' (C rounding of the binary value, as sprintf)
f1 <- function(x) sprintf("%.1f", x)
f0 <- function(x) sprintf("%.0f", x)

# Sum with Neumaier compensation, in the given order, as Python's sum() of floats (Python >= 3.12). The protein
# rule sums domain scores per label; a plain sum can differ from Python's in the last bit, which changes the
# rounding of sums that are exactly x.5 (runner_up and lineages are written with 0 decimals) and the order of
# exactly tied labels.
hc_pysum <- function(x) {
  f <- 0; c <- 0
  for (v in x) {
    t <- f + v
    if (abs(f) >= abs(v)) c <- c + ((f - t) + v) else c <- c + ((v - t) + f)
    f <- t
  }
  if (c != 0 && is.finite(c)) f <- f + c
  f
}

# Whitespace-split fields of each line (as Python's str.split()); lines starting with '#' dropped
hc_fields <- function(f) {
  if (!file.exists(f) || file.size(f) == 0) return(list())
  l <- readLines(f)
  l <- l[!startsWith(l, "#")]
  strsplit(trimws(l), "[[:space:]]+")
}

# Reads a TSV with a header, all columns as text (values are passed on as written)
hc_tsv <- function(f) fread(f, sep = "\t", colClasses = "character", na.strings = NULL, quote = "")

# Writes a data.table as TSV, as Python's '\t'.join (no quoting); all columns must be text
hc_write <- function(dt, f) {
  con <- file(f, "w")
  writeLines(paste(names(dt), collapse = "\t"), con)
  if (nrow(dt) > 0) writeLines(do.call(paste, c(as.list(dt), sep = "\t")), con)
  close(con)
}

# Reads a FASTA file: id (first word of the header) and sequence (lines joined, whitespace removed)
hc_fasta <- function(f) {
  l <- readLines(f)
  h <- which(startsWith(l, ">"))
  if (length(h) == 0) return(data.table(id = character(), seq = character()))
  grp <- findInterval(seq_along(l), h)
  ids <- vapply(strsplit(trimws(substring(l[h], 2)), "[[:space:]]+"), function(w) if (length(w)) w[1] else "", "")
  body <- l[-h]
  s <- vapply(split(body, factor(grp[-h], levels = seq_along(h))), function(x) paste(trimws(x), collapse = ""), "")
  data.table(id = ids, seq = unname(s))
}

hc_write_fasta <- function(ids, seqs, f) writeLines(paste0(">", ids, "\n", seqs), f)

# Reverse complement (ACGTN; N stays N)
hc_rc <- function(s) stringi::stri_reverse(chartr("ACGTN", "TGCAN", s))

# Length of the union of closed intervals [a, b] (as 1204's union())
hc_union <- function(a, b) {
  if (length(a) == 0) return(0)
  o <- order(a, b)
  a <- a[o]; b <- b[o]
  n <- 0; e <- 0
  for (i in seq_along(a)) {
    if (b[i] <= e) next
    n <- n + b[i] - max(a[i], e + 1) + 1
    e <- b[i]
  }
  n
}

# Decompresses a gzipped library file once into dir (tools read it once per chunk)
hc_unzip <- function(f, dir) {
  out <- file.path(dir, sub("\\.gz$", "", basename(f)))
  if (!file.exists(out)) {
    tmp <- paste0(out, ".tmp")
    hc_sys(paste("gzip -dc", shQuote(f), ">", shQuote(tmp)))
    file.rename(tmp, out)
  }
  out
}


### Input

# Splits the input FASTA into chunks. Ids as the step-9 protein piece writes them: first word of the header, cut
# at '#', '/' and '|' replaced by '_'. Sequences in upper case, anything but ACGTN becomes N.
hc_split <- function(input, chunk) {
  l <- readLines(input)
  l <- l[startsWith(l, ">") | nzchar(trimws(l))]           # blank lines are skipped
  h <- startsWith(l, ">")
  if (!h[1]) stop("input is not FASTA: ", input)
  w <- vapply(strsplit(trimws(substring(l[h], 2)), "[[:space:]]+"), function(x) if (length(x)) x[1] else "", "")
  if (any(w == "")) stop("empty sequence id in ", input)
  ids <- gsub("[/|]", "_", sub("#.*", "", w))
  bad <- ids == "" | duplicated(ids)
  if (any(bad)) stop("empty or duplicated sequence id after cleaning: ", w[bad][1], " -> ", ids[bad][1])
  grp <- cumsum(h)
  body <- l[!h]
  body <- gsub("[^ACGTN]", "N", toupper(trimws(body)))
  body_grp <- grp[!h]
  keep <- body != ""
  seqs <- vapply(split(body[keep], factor(body_grp[keep], levels = seq_along(ids))), paste, "", collapse = "")
  dt <- data.table(id = ids, input_id = w, seq = unname(seqs))
  dt[, chunk := sprintf("c%04d", (seq_len(.N) - 1L) %/% chunk)]
  dt
}


### Protein piece

# Codon table (standard code); codons with N translate to X
hc_code <- local({
  b <- c("T", "C", "A", "G")
  aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
  setNames(aa, paste0(rep(b, each = 16), rep(rep(b, each = 4), 4), rep(b, 16)))
})

# Translation of s from its first base, codon by codon (as 902's translate)
hc_translate <- function(s) {
  n <- nchar(s)
  if (n < 3) return(character())
  st <- seq(1L, n - 2L, by = 3L)
  aa <- hc_code[substring(s, st, st + 2L)]
  aa[is.na(aa)] <- "X"
  unname(aa)
}

# Reads miniprot --gff --trans output: one row per alignment (##PAF), with its translation (##STA) and CDS parts
hc_read_miniprot <- function(f) {
  l <- readLines(f)
  out <- list(); cur <- NULL
  for (x in l) {
    if (startsWith(x, "##PAF")) {
      if (!is.null(cur)) out[[length(out) + 1]] <- cur
      p <- strsplit(x, "\t", fixed = TRUE)[[1]]
      tg <- p[-(1:13)]                                  # SAM-like tags (AS:i:...) after the 13 PAF columns
      tg <- tg[lengths(regmatches(tg, gregexpr(":", tg))) >= 2]
      as_tag <- tg[startsWith(tg, "AS:")]
      cur <- list(ref = p[2], entry = p[7], strand = p[6], a = as.integer(p[9]) + 1L, b = as.integer(p[10]),
                  score = if (length(as_tag)) as.integer(sub("^AS:[^:]*:", "", as_tag[length(as_tag)])) else 0L,
                  cds_a = integer(), cds_b = integer(), sta = "")
    } else if (startsWith(x, "##STA") && !is.null(cur)) {
      cur$sta <- strsplit(x, "\t", fixed = TRUE)[[1]][2]
    } else if (!startsWith(x, "#") && !is.null(cur)) {
      p <- strsplit(x, "\t", fixed = TRUE)[[1]]
      if (length(p) > 8 && p[3] == "CDS") {
        cur$cds_a <- c(cur$cds_a, as.integer(p[4])); cur$cds_b <- c(cur$cds_b, as.integer(p[5]))
      }
    }
  }
  if (!is.null(cur)) out[[length(out) + 1]] <- cur
  out
}

hc_ov <- function(a1, b1, a2, b2) pmax(0L, pmin(b1, b2) - pmax(a1, a2) + 1L)

# Proteins of one sequence (as 902_proteins.py): miniprot translations of the reference proteins (one per locus,
# best score first) plus every stop-to-stop ORF of >= min_aa codons in the six frames that does not lie mostly in a
# kept miniprot locus on the same strand. Returns data.table(prot_id, location, aa_len, product, protein).
hc_proteins_one <- function(e, s, aln) {
  rows <- list()
  kept <- list()
  if (length(aln)) {
    o <- order(-vapply(aln, `[[`, 0L, "score"), vapply(aln, `[[`, 0L, "a"), vapply(aln, `[[`, "", "ref"),
               method = "radix")
    for (a in aln[o]) {
      if (!nzchar(a$sta) || length(a$cds_a) == 0) next
      clash <- FALSE
      for (k in kept) {
        if (hc_ov(a$a, a$b, k$a, k$b) > hmm_cfg$mp_overlap * min(a$b - a$a, k$b - k$a) + 1) { clash <- TRUE; break }
      }
      if (!clash) kept[[length(kept) + 1]] <- a
    }
  }
  if (length(kept)) {
    ka <- vapply(kept, `[[`, 0L, "a")
    kept <- kept[order(ka, method = "radix")]
    for (i in seq_along(kept)) {
      a <- kept[[i]]
      oc <- order(a$cds_a, a$cds_b, method = "radix")
      loc <- paste0(a$cds_a[oc], "..", a$cds_b[oc], collapse = ",")
      if (length(oc) > 1) loc <- paste0("join(", loc, ")")
      if (a$strand == "-") loc <- paste0("complement(", loc, ")")
      pr <- gsub("*", "X", a$sta, fixed = TRUE)
      rows[[length(rows) + 1]] <- list(paste0(e, "|mp", i, "|", a$a, "-", a$b), length(oc), loc, nchar(pr),
                                       paste0("miniprot:", a$ref), pr)
    }
  }
  n <- nchar(s)
  for (strand in c("+", "-")) {
    t <- if (strand == "+") s else hc_rc(s)
    for (f in 0:2) {
      aa <- hc_translate(substring(t, f + 1L))
      if (!length(aa)) next
      stops <- which(aa == "*")
      st <- c(1L, stops + 1L)
      en <- c(stops - 1L, length(aa))
      L <- en - st + 1L
      for (j in which(L >= hmm_cfg$min_aa)) {
        na <- f + 3L * (st[j] - 1L) + 1L
        nb <- na + 3L * L[j] - 1L
        if (strand == "+") { a1 <- na; b1 <- nb } else { a1 <- n - nb + 1L; b1 <- n - na + 1L }
        inside <- FALSE
        for (k in kept) {
          if (k$strand == strand && hc_ov(a1, b1, k$a, k$b) >= hmm_cfg$orf_in_mp * (b1 - a1 + 1)) { inside <- TRUE; break }
        }
        if (inside) next
        loc <- if (strand == "+") paste0(a1, "..", b1) else paste0("complement(", a1, "..", b1, ")")
        rows[[length(rows) + 1]] <- list(paste0(e, "|", strand, f + 1L, "|", a1, "-", b1), 1L, loc, L[j], "orf",
                                         paste(aa[st[j]:en[j]], collapse = ""))
      }
    }
  }
  if (!length(rows)) return(NULL)
  r <- rbindlist(lapply(rows, function(x) setNames(x, c("prot_id", "cds_n", "location", "aa_len", "product", "protein"))))
  r[, entry := e]
  r
}

# Element rule of the protein piece (as 66_element_calls.py): per sequence, a domain hit counts if it covers
# dom_min_cov of its model or scores dom_exempt bits; per label (TYPE/subtype) the best score of each domain class
# is summed and the highest sum is the call. Sub label (lineage) when the best lineage sum leads the next by
# lin_margin. Confidence 'low' for a call on one class with < low_score bits or on a host_like class alone.
hc_prot_calls <- function(dom_file, ids) {
  empty <- data.table(entry = ids, call = "", call_score = "", call_classes = "", confidence = "", lineage = "",
                      lineage_score = "", lineage_margin = "", runner_up = "", lineages = "")
  p <- hc_fields(dom_file)
  if (!length(p)) return(empty)
  d <- data.table(e = sub("\\|.*", "", vapply(p, `[`, "", 1L)), model = vapply(p, `[`, "", 4L),
                  qlen = as.integer(vapply(p, `[`, "", 6L)), sc = as.numeric(vapply(p, `[`, "", 14L)),
                  hf = as.integer(vapply(p, `[`, "", 16L)), ht = as.integer(vapply(p, `[`, "", 17L)))
  d[, ord := .I]
  d <- d[!((ht - hf + 1) / qlen < hmm_cfg$dom_min_cov & sc < hmm_cfg$dom_exempt) & sc > 0]
  # best score per sequence and model; the order of first appearance is kept (Python dict order)
  b <- d[, .(sc = max(sc), ord = min(ord)), by = .(e, model)][order(ord)]
  b[, full := sub(":.*", "", model)]
  b[, cls := sub(":.*", "", sub("^[^:]*:", "", model))]
  parts <- strsplit(b$full, "/", fixed = TRUE)
  b[, lab := vapply(parts, function(x) paste(x[seq_len(min(2L, length(x)))], collapse = "/"), "")]
  b[, lin := vapply(parts, function(x) if (length(x) > 2) x[3] else NA_character_, "")]
  res <- lapply(ids, function(x) {
    bx <- b[e == x]
    if (!nrow(bx)) return(NULL)
    cl <- bx[, .(s = max(sc), ord = min(ord)), by = .(lab, cls)][order(ord)]
    tot <- cl[, .(s = hc_pysum(s)), by = lab]                # classes in order of first appearance
    setorder(tot, -s, -lab)
    call <- tot$lab[1]; s <- tot$s[1]
    cs <- sort(cl[lab == call]$cls, method = "radix")
    conf <- if (length(cs) == 1 && (s < hmm_cfg$low_score || cs[1] %in% hmm_cfg$host_like)) "low" else "high"
    ln <- bx[lab == call & !is.na(lin)][, .(s = max(sc), ord = min(ord)), by = .(lin, cls)][order(ord)]
    lt <- ln[, .(s = hc_pysum(s)), by = lin]
    setorder(lt, -s, -lin)
    lin_call <- ""; ls <- ""; lm <- ""
    if (nrow(lt)) {
      lsv <- lt$s[1]; lmv <- lsv - (if (nrow(lt) > 1) lt$s[2] else 0)
      if (lmv >= hmm_cfg$lin_margin) lin_call <- lt$lin[1]
      ls <- f1(lsv); lm <- f1(lmv)
    }
    data.table(entry = x, call = call, call_score = f1(s), call_classes = paste(cs, collapse = ","),
               confidence = conf, lineage = lin_call, lineage_score = ls, lineage_margin = lm,
               runner_up = if (nrow(tot) > 1) paste0(tot$lab[2], ":", f0(tot$s[2])) else "",
               lineages = if (nrow(lt)) paste0(lt$lin, ":", f0(lt$s), collapse = ";") else "")
  })
  res <- rbindlist(res)
  if (!nrow(res)) return(empty)
  out <- res[empty[, .(entry)], on = "entry"]
  for (c in names(out)) set(out, which(is.na(out[[c]])), c, "")
  out
}

# Protein piece of one chunk: proteins, hmmsearch, element calls. Writes prot/{proteins.faa, proteins.tsv,
# prot.dom, prot_calls.tsv}.
hc_piece_prot <- function(ch, dir, lib) {
  o <- file.path(dir, "prot"); dir.create(o, showWarnings = FALSE)
  ids <- ch$id
  aln <- list()
  for (i in seq_len(nrow(ch))) {
    one <- file.path(o, "one.fna")
    hc_write_fasta(ch$id[i], ch$seq[i], one)
    gff <- file.path(o, "one.gff")
    hc_sys(paste("miniprot -t 1", hmm_cfg$mp_opts, "--gff --trans", shQuote(one), shQuote(lib$ref),
                 ">", shQuote(gff), "2>>", shQuote(file.path(o, "miniprot.log"))))
    aln[[ch$id[i]]] <- hc_read_miniprot(gff)
  }
  unlink(file.path(o, c("one.fna", "one.gff")))
  so <- order(ids, method = "radix")
  prots <- rbindlist(lapply(so, function(i) {
    a <- aln[[ids[i]]]
    a <- a[vapply(a, function(x) x$entry == ids[i], TRUE)]
    hc_proteins_one(ids[i], ch$seq[i], a)
  }))
  faa <- file.path(o, "proteins.faa")
  if (nrow(prots)) {
    wrap <- vapply(prots$protein, function(p) {
      n <- nchar(p); st <- seq(1L, n, by = 80L)
      paste(substring(p, st, pmin(st + 79L, n)), collapse = "\n")
    }, "")
    writeLines(paste0(">", prots$prot_id, "\n", wrap), faa)
  } else file.create(faa)
  tsv <- if (nrow(prots)) prots[, .(prot_id, entry, class_kw = "", pantera_label = "", label_from = "dfam",
                                    cds_n = as.character(cds_n), location, aa_len = as.character(aa_len), product,
                                    note = "", species = "", lineage = "")] else
    data.table(prot_id = character(), entry = character(), class_kw = character(), pantera_label = character(),
               label_from = character(), cds_n = character(), location = character(), aa_len = character(),
               product = character(), note = character(), species = character(), lineage = character())
  hc_write(tsv, file.path(o, "proteins.tsv"))
  dom <- file.path(o, "prot.dom")
  if (nrow(prots)) {
    hc_sys(paste("hmmsearch --cpu 1 --noali -o /dev/null --domtblout", shQuote(dom), "-Z", hmm_cfg$hmm_z,
                 "--domZ", hmm_cfg$hmm_z, "--cut_ga", shQuote(lib$prot_hmm), shQuote(faa)))
  } else file.create(dom)
  calls <- hc_prot_calls(dom, ids)
  hc_write(calls, file.path(o, "prot_calls.tsv"))
  calls
}


### End windows (TIR and Helitron pieces)

hc_windows <- function(seqs, w) ifelse(nchar(seqs) <= 2L * w, seqs, paste0(substr(seqs, 1L, w), substring(seqs, nchar(seqs) - w + 1L)))

# nhmmer table: target, query, alignment from/to, sequence length, strand, score
hc_nhmmer_tbl <- function(f) {
  p <- hc_fields(f)
  p <- p[lengths(p) >= 15]
  if (!length(p)) return(data.table(t = character(), q = character(), a1 = integer(), a2 = integer(),
                                    sl = integer(), st = character(), sc = numeric(), ord = integer()))
  data.table(t = vapply(p, `[`, "", 1L), q = vapply(p, `[`, "", 3L), a1 = as.integer(vapply(p, `[`, "", 7L)),
             a2 = as.integer(vapply(p, `[`, "", 8L)), sl = as.integer(vapply(p, `[`, "", 11L)),
             st = vapply(p, `[`, "", 12L), sc = as.numeric(vapply(p, `[`, "", 14L)))[, ord := .I][]
}

# Best hit per group, first in file order among equal scores (Python's 'if sc > best' updates)
hc_best <- function(d, by) {
  if (!nrow(d)) return(d)
  d <- d[order(-sc, ord)]
  d[!duplicated(d, by = by)]
}


### TIR piece

# TIR end scores (805) and calls (812, tircall.py), with the sub label (911, tirsub.py).
# Left end: + strand hit starting within end_anchor bp of the window start; right end: - strand hit ending within
# end_anchor bp of the window end. Superfamily score: its best left + best right hits ('mixed' models count for none).
hc_piece_tir <- function(ch, dir, lib) {
  o <- file.path(dir, "tir"); dir.create(o, showWarnings = FALSE)
  win <- hc_windows(ch$seq, hmm_cfg$end_win)
  wl <- setNames(nchar(win), ch$id)
  hc_write_fasta(ch$id, win, file.path(o, "windows.fna"))
  hc_sys(paste("nhmmer --cpu 1 --dna -T", hmm_cfg$end_hit, "--noali --tblout", shQuote(file.path(o, "hits.tbl")),
               "-o /dev/null", shQuote(lib$tir_hmm), shQuote(file.path(o, "windows.fna"))))
  h <- hc_nhmmer_tbl(file.path(o, "hits.tbl"))[sc >= hmm_cfg$end_hit]
  h[, L := wl[t]]
  h[, side := fifelse(st == "+" & pmin(a1, a2) <= hmm_cfg$end_anchor, "L",
                      fifelse(st == "-" & pmax(a1, a2) > L - hmm_cfg$end_anchor, "R", NA_character_))]
  h <- h[!is.na(side)]
  h[, sf := lib$tir_sf[q]]
  h[, sf_sub := lib$sub_sf[q]]
  h[, sub := lib$sub_sub[q]]
  out <- data.table(id = ch$id, length = as.character(nchar(ch$seq)), sL = "0", sR = "0", ends = "", call = "none",
                    superfamily = "", confidence = "", second = "", second_score = "", modelL = "", modelR = "",
                    sub = "", sub_score = "", sub_second = "")
  if (nrow(h)) {
    bs <- hc_best(h, c("t", "side"))                       # best per side, any model
    bsf <- hc_best(h, c("t", "side", "sf"))                # best per side and superfamily (805)
    for (x in unique(h$t)) {
      i <- match(x, out$id)
      sL <- bs[t == x & side == "L"]; sR <- bs[t == x & side == "R"]
      vL <- if (nrow(sL)) sL$sc else 0; vR <- if (nrow(sR)) sR$sc else 0
      g <- bsf[t == x & sf != "mixed", .(s = sum(sc)), by = sf]
      setorder(g, -s, -sf)
      s1 <- if (nrow(g) >= 1) g$s[1] else 0; s2 <- if (nrow(g) >= 2) g$s[2] else 0
      n1 <- if (nrow(g) >= 1) g$sf[1] else ""; n2 <- if (nrow(g) >= 2) g$sf[2] else ""
      both <- vL > 0 && vR > 0 && vL + vR >= hmm_cfg$tir_sum
      set(out, i, c("sL", "sR"), list(f1(vL), f1(vR)))
      if (!(both || max(vL, vR) >= hmm_cfg$tir_single)) next
      sfc <- if (nzchar(n1) && s1 > 0 && s1 - s2 >= hmm_cfg$tir_lead) n1 else "DNA"
      conf <- if (both && vL + vR >= hmm_cfg$tir_high) "high" else "low"
      set(out, i, c("ends", "call", "superfamily", "confidence", "second", "second_score", "modelL", "modelR"),
          list(if (both) "both" else "one", "TIR", sfc, conf, n2, f1(s2),
               if (nrow(sL)) sL$q else "", if (nrow(sR)) sR$q else ""))
      if (!sfc %in% c("", "DNA")) {
        # sub label (tirsub.sub_call): best sub of the superfamily, leading the second by tir_sub_lead and within
        # tir_sub_lead of the superfamily score (superfamilies here from the sub table, as 911)
        hs <- h[t == x]
        bsub <- hc_best(hs[!is.na(sub) & nzchar(sub) & sf_sub == sfc], c("side", "sub"))
        u <- bsub[, .(s = sum(sc)), by = sub]
        setorder(u, -s, -sub)
        ss <- if (nrow(u) >= 1) u$s[1] else 0; s2b <- if (nrow(u) >= 2) u$s[2] else 0
        sfs <- hc_best(hs[sf_sub == sfc], "side")[, sum(sc)]
        ok <- ss > 0 && ss - s2b >= hmm_cfg$tir_sub_lead && sfs - ss <= hmm_cfg$tir_sub_lead
        sb <- if (ok) u$sub[1] else ""
        if (nzchar(sb) || ss != 0) set(out, i, c("sub", "sub_score", "sub_second"), list(sb, f1(ss), f1(s2b)))
      }
    }
  }
  hc_write(out, file.path(o, "tir_calls.tsv"))
  out
}


### Helitron piece

# Helitron end scores (705) and calls (712, helcall.py). Per strand of the element: best anchored 5' and 3' end
# hits; the strand with the higher s5 + s3 is kept (the first one in the table when equal).
hc_piece_hel <- function(ch, dir, lib) {
  o <- file.path(dir, "hel"); dir.create(o, showWarnings = FALSE)
  win <- hc_windows(ch$seq, hmm_cfg$end_win)
  wl <- setNames(nchar(win), ch$id)
  hc_write_fasta(ch$id, win, file.path(o, "windows.fna"))
  hc_sys(paste("nhmmer --cpu 1 --dna -T", hmm_cfg$end_hit, "--noali --tblout", shQuote(file.path(o, "nhmmer.tbl")),
               "-o /dev/null", shQuote(lib$hel_hmm), shQuote(file.path(o, "windows.fna"))))
  h <- hc_nhmmer_tbl(file.path(o, "nhmmer.tbl"))[sc >= hmm_cfg$end_hit]
  h[, L := wl[t]]
  h[, end := lib$hel_end[q]]
  h[, grp := lib$hel_group[q]]
  A <- hmm_cfg$end_anchor
  h <- h[fifelse(st == "+", fifelse(end == "5", a1 <= A, a2 > L - A), fifelse(end == "5", a1 > L - A, a2 <= A))]
  out <- data.table(id = ch$id, length = as.character(nchar(ch$seq)), strand = "", s5 = "0", s3 = "0", ends = "",
                    call = "none", group = "", confidence = "", model5 = "", model3 = "")
  if (nrow(h)) {
    first <- h[, .(ord = min(ord)), by = .(t, st)][order(ord)]       # (sequence, strand) in table order
    be <- hc_best(h, c("t", "st", "end"))
    bg <- hc_best(h[grp %in% c("Helitron-1", "Helentron")], c("t", "st", "end", "grp"))
    for (x in unique(first$t)) {
      best <- NULL
      for (s in first[t == x]$st) {
        e5 <- be[t == x & st == s & end == "5"]; e3 <- be[t == x & st == s & end == "3"]
        v5 <- if (nrow(e5)) e5$sc else 0; v3 <- if (nrow(e3)) e3$sc else 0
        gs <- vapply(c("Helitron-1", "Helentron"), function(g) bg[t == x & st == s & grp == g, sum(sc)], 0)
        r <- list(strand = s, s5 = v5, s3 = v3, score = v5 + v3, m5 = if (nrow(e5)) e5$q else "",
                  m3 = if (nrow(e3)) e3$q else "", g = gs)
        if (is.null(best) || r$score > best$score) best <- r
      }
      i <- match(x, out$id)
      set(out, i, c("strand", "s5", "s3", "model5", "model3"), list(best$strand, f1(best$s5), f1(best$s3), best$m5, best$m3))
      both <- best$s5 > 0 && best$s3 > 0 && best$s5 + best$s3 >= hmm_cfg$hel_sum
      if (!(both || max(best$s5, best$s3) >= hmm_cfg$hel_single)) next
      g <- best$g
      ab <- names(g)[order(-g, method = "radix")]               # stable: Helitron-1 first when equal
      grp_call <- if (g[ab[1]] > 0 && g[ab[1]] - g[ab[2]] >= hmm_cfg$hel_glead) ab[1] else "RC"
      conf <- if (both && best$s5 + best$s3 >= hmm_cfg$hel_high) "high" else "low"
      set(out, i, c("ends", "call", "group", "confidence"), list(if (both) "both" else "one", "Helitron", grp_call, conf))
    }
  }
  hc_write(out, file.path(o, "hel_calls.tsv"))
  out
}


### SINE piece

# SINE calls (76_calls.py). Per sequence and strand, per label: best HEAD hit starting within sine_head_window bp
# of the 5' end plus best BODY hit. The call is the label with the highest sum over both strands, if >=
# sine_min_total and leading the next label by sine_lead ('ambiguous' otherwise). rna: output of hc_rna_check.
hc_sine_calls <- function(tbl, models, rna = NULL) {
  h <- hc_nhmmer_tbl(tbl)[sc >= hmm_cfg$sine_hit]
  h[, kind := models$kind[match(q, models$model)]]
  h[, lab := models$label[match(q, models$model)]]
  h <- h[!(kind == "HEAD" & fifelse(st == "+", a1, sl - a1 + 1L) > hmm_cfg$sine_head_window)]
  h[, `:=`(lo = pmin(a1, a2), hi = pmax(a1, a2))]
  if (!nrow(h)) return(data.table(id = character(), call = character(), confidence = character(),
                                  best_label = character(), score = character(), lead = character(),
                                  coverage = character(), strand = character(), head_model = character(),
                                  body_model = character(), top_labels = character(), rfam_best = character()))
  b <- hc_best(h, c("t", "st", "lab", "kind"))
  so <- h[, .(ord = min(ord)), by = .(t, st)]
  b[so, on = c("t", "st"), sord := i.ord]
  res <- lapply(sort(unique(h$t), method = "radix"), function(x) {
    bx <- b[t == x]
    elen <- h[t == x][order(ord)]$sl[1]
    # per strand (table order), per label: HEAD + BODY. Labels of one strand in order of name (Python iterates a
    # set here, in hash order; this only matters for exact ties, see README)
    tot <- NULL
    for (s in unique(bx[order(sord)]$st)) {
      bs <- bx[st == s]
      for (lb in sort(unique(bs$lab), method = "radix")) {
        hd <- bs[lab == lb & kind == "HEAD"]; bd <- bs[lab == lb & kind == "BODY"]
        v <- (if (nrow(hd)) hd$sc else 0) + (if (nrow(bd)) bd$sc else 0)
        hs <- rbind(hd, bd)
        cov <- (max(hs$hi) - min(hs$lo) + 1) / hs$sl[1]
        k <- if (is.null(tot)) NA_integer_ else match(lb, tot$lab)
        if (is.na(k)) {
          tot <- rbind(tot, data.table(lab = lb, s = v, strand = s, hm = if (nrow(hd)) hd$q else "",
                                       bm = if (nrow(bd)) bd$q else "", cov = cov))
        } else if (v > tot$s[k]) {
          set(tot, k, c("s", "strand", "hm", "bm", "cov"), list(v, s, if (nrow(hd)) hd$q else "", if (nrow(bd)) bd$q else "", cov))
        }
      }
    }
    r <- tot[order(-s)]                                     # stable, as Python's sorted(key=-score)
    lead <- r$s[1] - (if (nrow(r) > 1) r$s[2] else 0)
    call <- if (r$s[1] < hmm_cfg$sine_min_total) "none" else if (lead < hmm_cfg$sine_lead) "ambiguous" else r$lab[1]
    if (startsWith(call, "SINE") && elen > hmm_cfg$sine_max_len) call <- "SINE_part"
    conf <- if (startsWith(call, "SINE/") && nzchar(r$hm[1]) && r$s[1] >= hmm_cfg$sine_high) "high" else
      if (startsWith(call, "SINE/")) "low" else ""
    rg <- if (!is.null(rna)) rna[id == x] else NULL
    if (call != "none" && !is.null(rg) && nrow(rg) && as.numeric(rg$match_frac) >= hmm_cfg$sine_rna_frac) call <- "RNA_gene-like"
    top <- r[seq_len(min(3L, nrow(r)))]
    data.table(id = x, call = call, confidence = conf, best_label = r$lab[1], score = f1(r$s[1]), lead = f1(lead),
               coverage = sprintf("%.2f", r$cov[1]), strand = r$strand[1], head_model = r$hm[1], body_model = r$bm[1],
               top_labels = paste0(top$lab, ":", f0(top$s), collapse = ";"),
               rfam_best = if (!is.null(rg) && nrow(rg)) paste0(rg$rfam_model, ":", rg$score, ":fill", rg$match_frac) else "")
  })
  rbindlist(res)
}

# RNA-gene check (76a_rna_check.py): candidates are sequences whose best Rfam head hit (>= rna_min bits) spans >=
# rna_span of its model and ends < rna_rest bp before the sequence end. The hit region is aligned to its CM
# (cmalign, global) and match_frac = occupied consensus columns / model length.
hc_rna_check <- function(tbl, seqs, cm, dir) {
  p <- hc_fields(tbl)
  p <- p[lengths(p) >= 28]
  empty <- data.table(id = character(), rfam_model = character(), score = character(), span_frac = character(),
                      rest = character(), match_frac = character())
  if (!length(p)) return(empty)
  d <- data.table(model = vapply(p, `[`, "", 2L), q = vapply(p, `[`, "", 4L),
                  mf = as.integer(vapply(p, `[`, "", 8L)), mt = as.integer(vapply(p, `[`, "", 9L)),
                  sf = as.integer(vapply(p, `[`, "", 10L)), st = as.integer(vapply(p, `[`, "", 11L)),
                  strand = vapply(p, `[`, "", 12L), sc = as.numeric(vapply(p, `[`, "", 17L)),
                  ml = as.integer(vapply(p, `[`, "", 27L)), sl = as.integer(vapply(p, `[`, "", 28L)))
  d[, ord := .I]
  d <- d[sc >= hmm_cfg$rna_min]
  if (!nrow(d)) return(empty)
  first <- d[, .(first = min(ord)), by = q]               # Python dict order: first hit of each sequence
  b <- hc_best(d, "q")[first, on = "q"][order(first)]
  b[, span := (mt - mf + 1) / ml]
  b[, rest := fifelse(strand == "+", sl - st, st - 1L)]
  cand <- b[span >= hmm_cfg$rna_span & rest < hmm_cfg$rna_rest]
  if (!nrow(cand)) return(empty)
  cand[, s := seqs[q]]
  cand[, sub := fifelse(strand == "+", substr(s, sf, st), hc_rc(substr(s, st, sf)))]
  cand[, frac := NA_real_]
  for (m in sort(unique(cand$model), method = "radix")) {
    it <- which(cand$model == m)
    cmf <- file.path(dir, "m.cm"); fa <- file.path(dir, "q.fna"); sto <- file.path(dir, "a.sto")
    hc_sys(paste("cmfetch", shQuote(cm), shQuote(m), ">", shQuote(cmf)))
    hc_write_fasta(paste0("s", seq_along(it) - 1L), cand$sub[it], fa)
    hc_sys(paste("cmalign --cpu 1 --outformat Pfam -o", shQuote(sto), shQuote(cmf), shQuote(fa), "> /dev/null 2>&1"))
    l <- readLines(sto)
    rf <- sub(".*[[:space:]]", "", l[startsWith(l, "#=GC RF")])
    al <- l[nzchar(trimws(l)) & !startsWith(l, "#") & !startsWith(l, "//")]
    al <- strsplit(trimws(al), "[[:space:]]+")
    an <- setNames(vapply(al, `[`, "", 2L), vapply(al, `[`, "", 1L))
    cols <- which(!strsplit(rf, "")[[1]] %in% c(".", "-", "~"))
    for (k in seq_along(it)) {
      a <- strsplit(an[[paste0("s", k - 1L)]], "")[[1]]
      set(cand, it[k], "frac", sum(!a[cols] %in% c(".", "-", "~")) / length(cols))
    }
    unlink(c(cmf, fa, sto))
  }
  cand <- cand[order(q, method = "radix")]
  cand[, .(id = q, rfam_model = model, score = f1(sc), span_frac = sprintf("%.2f", span), rest = as.character(rest),
           match_frac = sprintf("%.3f", frac))]
}

# SINE piece of one chunk (as step-10 1005): nhmmer, first calls, head CMs on the called sequences, RNA check,
# final calls. Writes sine/{all.tbl, calls0.tsv, rna_check.tsv, calls.tsv}.
hc_piece_sine <- function(ch, fna, dir, lib) {
  o <- file.path(dir, "sine"); dir.create(o, showWarnings = FALSE)
  tbl <- file.path(o, "all.tbl")
  hc_sys(paste("nhmmer --cpu 1 --dna -T", hmm_cfg$sine_hit, "--noali --tblout", shQuote(tbl), "-o /dev/null",
               shQuote(lib$sine_hmm), shQuote(fna)))
  c0 <- hc_sine_calls(tbl, lib$sine_models)
  hc_write(c0, file.path(o, "calls0.tsv"))
  called <- c0[call != "none"]$id
  rna <- data.table(id = character(), rfam_model = character(), score = character(), span_frac = character(),
                    rest = character(), match_frac = character())
  if (length(called)) {
    sub <- ch[id %in% called]
    cf <- file.path(o, "called.fna")
    hc_write_fasta(sub$id, sub$seq, cf)
    ct <- file.path(o, "cmscan.tbl")
    hc_sys(paste("cmscan --cpu 1 --nohmmonly --notrunc -T", hmm_cfg$rna_min, "--tblout", shQuote(ct),
                 "--fmt 2 -o /dev/null", shQuote(lib$heads_cm), shQuote(cf)))
    rna <- hc_rna_check(ct, setNames(sub$seq, sub$id), lib$heads_cm, o)
  }
  hc_write(rna, file.path(o, "rna_check.tsv"))
  calls <- hc_sine_calls(tbl, lib$sine_models, rna)
  hc_write(calls, file.path(o, "calls.tsv"))
  calls
}


### Flags piece

# Tandem (TRF), low complexity (dustmasker), ncRNA (cmscan of Rfam families) and host proteins (DIAMOND, optional)
# of one chunk. Returns per sequence the covered lengths and hits used by hc_combine.
hc_piece_flags <- function(ch, fna, dir, lib, prot_dir) {
  o <- file.path(dir, "flags"); dir.create(o, showWarnings = FALSE)
  trf <- file.path(o, "trf.ngs")
  hc_sys(paste("cd", shQuote(o), "&& trf", shQuote(normalizePath(fna)), hmm_cfg$trf_opts, ">", shQuote(basename(trf)),
               "2> /dev/null"), ok_status = 0:255)               # TRF exits non-zero by design
  hc_sys(paste("dustmasker -in", shQuote(fna), "-outfmt interval -out", shQuote(file.path(o, "dust.txt"))))
  hc_sys(paste("cmscan --cpu 1 --rfam --cut_ga --nohmmonly --noali --fmt 2 --tblout", shQuote(file.path(o, "rfam.tbl")),
               "-o /dev/null", shQuote(lib$ncrna_cm), shQuote(fna)))
  faa <- file.path(prot_dir, "proteins.faa")
  m8 <- file.path(o, "host.m8")
  if (!is.null(lib$host_dmnd) && file.exists(faa) && file.size(faa) > 0) {
    hc_sys(paste("diamond blastp --sensitive --threads 1 -q", shQuote(faa), "-d", shQuote(sub("\\.dmnd$", "", lib$host_dmnd)),
                 "--evalue", hmm_cfg$host_evalue, "--max-target-seqs 5 --outfmt 6 qseqid sseqid pident length qlen",
                 "slen qstart qend sstart send evalue bitscore -o", shQuote(m8), "--quiet"))
  } else file.create(m8)
  invisible(NULL)
}

# Flags of one chunk from its tool outputs (as 1204): tandem and low-complexity intervals, Rfam hits, best
# host-protein hit per sequence
hc_read_flags <- function(dir) {
  o <- file.path(dir, "flags")
  # TRF -ngs: '@name' lines, then one line per array (start end ...)
  l <- readLines(file.path(o, "trf.ngs"))
  hd <- startsWith(l, "@")
  nm <- cumsum(hd)
  names_ <- vapply(strsplit(substring(l[hd], 2), "[[:space:]]+"), `[`, "", 1L)
  p <- strsplit(trimws(l[!hd]), "[[:space:]]+")
  keep <- lengths(p) > 2
  tandem <- data.table(id = names_[nm[!hd]][keep], a = as.integer(vapply(p[keep], `[`, "", 1L)),
                       b = as.integer(vapply(p[keep], `[`, "", 2L)))
  # dustmasker intervals (0-based): '>name' lines, then 'a - b'
  l <- readLines(file.path(o, "dust.txt"))
  hd <- startsWith(l, ">")
  nm <- cumsum(hd)
  names_ <- vapply(strsplit(substring(l[hd], 2), "[[:space:]]+"), `[`, "", 1L)
  p <- strsplit(trimws(l[!hd]), "[[:space:]]+")
  lowc <- data.table(id = names_[nm[!hd]], a = as.integer(vapply(p, `[`, "", 1L)) + 1L,
                     b = as.integer(vapply(p, `[`, "", 3L)) + 1L)
  # Rfam hits (cmscan --fmt 2); hits marked '=' in the olp column (overlap with a better hit) are skipped
  p <- hc_fields(file.path(o, "rfam.tbl"))
  p <- p[vapply(p, function(x) length(x) >= 20 && x[20] != "=", TRUE)]
  rfam <- if (length(p)) data.table(id = vapply(p, `[`, "", 4L), fam = vapply(p, `[`, "", 2L),
                                    x1 = as.integer(vapply(p, `[`, "", 10L)), x2 = as.integer(vapply(p, `[`, "", 11L))) else
    data.table(id = character(), fam = character(), x1 = integer(), x2 = integer())
  rfam[, `:=`(a = pmin(x1, x2), b = pmax(x1, x2))]
  list(tandem = tandem, lowc = lowc, rfam = rfam)
}

# Best host-protein hit per sequence (as 1204): proteins with a pantera domain hit are skipped; a hit counts at
# >= host_pid % identity over >= host_cov of the Swiss-Prot protein; best by (bitscore, subject, identity,
# coverage, length), as Python's tuple comparison
hc_read_host <- function(dir) {
  m8 <- file.path(dir, "flags", "host.m8")
  if (!file.exists(m8) || file.size(m8) == 0) return(NULL)
  h <- fread(m8, sep = "\t", header = FALSE, colClasses = "character",
             col.names = c("q", "s", "pid", "ln", "ql", "sl", "qs", "qe", "ss", "se", "ev", "bs"))
  dom <- file.path(dir, "prot", "prot.dom")
  tehit <- unique(vapply(hc_fields(dom), `[`, "", 1L))
  h <- h[!q %in% tehit]
  h[, cov := (as.integer(se) - as.integer(ss) + 1) / as.integer(sl)]
  h <- h[as.numeric(pid) >= hmm_cfg$host_pid & cov >= hmm_cfg$host_cov]
  if (!nrow(h)) return(NULL)
  pe <- hc_tsv(file.path(dir, "prot", "proteins.tsv"))
  h[, entry := pe$entry[match(q, pe$prot_id)]]
  h[, `:=`(bsn = as.numeric(bs), pidn = as.numeric(pid), lnn = as.integer(ln))]
  setorder(h, -bsn, -s, -pidn, -cov, -lnn)                 # data.table sorts text in C locale, as Python
  h[!duplicated(entry), .(entry, s, lnn)]
}


### End partners (whole library)

# Partners (1203): sequences without a protein call (>= partner_min_len bp) are searched by their ends against the
# ends of (a) sequences of the same library with a high-confidence DNA protein call that includes a TPase domain
# and (b) curated autonomous DNA transposons. A hit counts only when anchored at both outer ends (within
# partner_anchor bp of the query and partner termini). Call: best partner with both ends, if its label leads the
# best partner of another label by partner_lead bits ('ambiguous' otherwise).
hc_partners <- function(all, prot, dir, lib, threads) {
  o <- file.path(dir, "partners"); dir.create(o, showWarnings = FALSE)
  W <- hmm_cfg$partner_win; A <- hmm_cfg$partner_anchor
  seq <- setNames(gsub("[^ACGT]", "N", all$seq), all$id)
  win <- function(s) if (nchar(s) < 2 * W) { k <- (nchar(s) + 1) %/% 2; c(substr(s, 1, k), substring(s, k + 1)) } else
    c(substr(s, 1, W), substring(s, nchar(s) - W + 1))
  P <- prot[match(all$id, entry)]
  queries <- all$id[(is.na(P$call) | !nzchar(P$call)) & nchar(all$seq) >= hmm_cfg$partner_min_len]
  lp <- prot[startsWith(call, "DNA/") & confidence == "high" &
               vapply(strsplit(call_classes, ",", fixed = TRUE), function(x) "TPase" %in% x, TRUE)]
  lp <- lp[order(entry, method = "radix")]
  cur <- hc_fasta(lib$partners)
  hp <- strsplit(cur$id, "|", fixed = TRUE)
  cur[, `:=`(e = vapply(hp, `[`, "", 2L), end = vapply(hp, `[`, "", 3L), lab = vapply(hp, `[`, "", 4L),
             sb = vapply(hp, `[`, "", 5L))]
  cur[sb == "-", sb := ""]
  plab <- rbind(data.table(src = "v", p = lp$entry, L = lp$call, S = lp$lineage),
                unique(cur[, .(src = "c", p = e, L = lab, S = sb)], by = "p", fromLast = TRUE))
  db <- file.path(o, "db.fna")
  con <- file(db, "w")
  for (x in lp$entry) { w <- win(seq[[x]]); writeLines(c(paste0(">v|", x, "|5"), w[1], paste0(">v|", x, "|3"), w[2]), con) }
  writeLines(paste0(">c|", cur$e, "|", cur$end, "\n", cur$seq), con)
  close(con)
  hc_sys(paste("makeblastdb -in", shQuote(db), "-dbtype nucl -out", shQuote(file.path(o, "db")), "> /dev/null"))
  blocks <- split(queries, (seq_along(queries) - 1L) %/% hmm_cfg$partner_block)
  fmt <- "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"
  res <- mclapply(seq_along(blocks), function(k) {
    q <- file.path(o, sprintf("q%04d.fna", k - 1L))
    con <- file(q, "w")
    for (x in blocks[[k]]) { w <- win(seq[[x]]); writeLines(c(paste0(">", x, "|5"), w[1], paste0(">", x, "|3"), w[2]), con) }
    close(con)
    m8 <- paste0(q, ".m8")
    hc_sys(paste("blastn -task blastn -query", shQuote(q), "-db", shQuote(file.path(o, "db")), "-evalue",
                 hmm_cfg$partner_evalue, "-max_target_seqs 1000 -num_threads 1 -outfmt", shQuote(fmt), ">", shQuote(m8)))
    m8
  }, mc.cores = threads, mc.preschedule = FALSE)
  if (any(vapply(res, inherits, TRUE, "try-error"))) stop("partner search failed: ", res[vapply(res, inherits, TRUE, "try-error")][[1]])
  hits <- rbindlist(lapply(res, function(f) if (file.size(f) > 0) fread(f, sep = "\t", header = FALSE,
                                                                      colClasses = "character") else NULL))
  cols <- c("id", "ends", "partner", "partner_source", "partner_label", "partner_sub", "score", "score_5", "score_3",
            "lead", "second_label", "call", "library_partner", "library_label", "library_sub", "library_score")
  out <- data.table(matrix("", nrow = length(queries), ncol = length(cols), dimnames = list(NULL, cols)))
  out[, id := queries]
  if (nrow(hits)) {
    setnames(hits, c("q", "s", "pid", "ln", "qs", "qe", "ss", "se", "ev", "bs", "ql", "sl"))
    hits[, `:=`(x = sub("\\|[^|]*$", "", q), qend = sub(".*\\|", "", q))]
    sp <- strsplit(hits$s, "|", fixed = TRUE)
    hits[, `:=`(src = vapply(sp, `[`, "", 1L), p = vapply(sp, `[`, "", 2L), pe = vapply(sp, `[`, "", 3L))]
    for (c in c("qs", "qe", "ss", "se", "ql", "sl")) set(hits, j = c, value = as.integer(hits[[c]]))
    hits[, bsn := as.numeric(bs)]
    hits <- hits[!(qend == "5" & qs > A) & !(qend == "3" & qe < ql - A)]
    hits[, `:=`(lo = pmin(ss, se), hi = pmax(ss, se))]
    hits <- hits[!(pe == "5" & lo > A) & !(pe == "3" & hi < sl - A)]
    b <- hits[bsn > 0, .(v = max(bsn)), by = .(x, src, p, qend)]
  }
  if (nrow(hits) && nrow(b)) {
    w <- dcast(b, x + src + p ~ qend, value.var = "v")
    if (!"5" %in% names(w)) w[, `5` := NA_real_]
    if (!"3" %in% names(w)) w[, `3` := NA_real_]
    w[, both := !is.na(`5`) & !is.na(`3`)]
    w[, s := fifelse(is.na(`5`), 0, `5`) + fifelse(is.na(`3`), 0, `3`)]
    w[plab, on = c("src", "p"), `:=`(L = i.L, S = i.S)]
    setorder(w, x, -both, -s, -src, -p)
    for (xq in unique(w$x)) {
      i <- match(xq, out$id)
      wx <- w[x == xq]
      t1 <- wx[1]
      set(out, i, c("ends", "partner", "partner_source", "partner_label", "partner_sub", "score", "score_5", "score_3"),
          list(if (t1$both) "both" else "one", t1$p, if (t1$src == "v") "library" else "curated", t1$L, t1$S, f1(t1$s),
               f1(if (is.na(t1$`5`)) 0 else t1$`5`), f1(if (is.na(t1$`3`)) 0 else t1$`3`)))
      if (t1$both) {
        oth <- wx[both & !(L == t1$L & S == t1$S)]
        s2 <- if (nrow(oth)) oth$s[1] else 0
        set(out, i, c("second_label", "lead", "call"),
            list(if (nrow(oth)) paste(c(oth$L[1], oth$S[1])[nzchar(c(oth$L[1], oth$S[1]))], collapse = "/") else "",
                 f1(t1$s - s2), if (t1$s - s2 >= hmm_cfg$partner_lead) t1$L else "ambiguous"))
        sg <- wx[both & src == "v"]
        if (nrow(sg)) set(out, i, c("library_partner", "library_label", "library_sub", "library_score"),
                          list(sg$p[1], sg$L[1], sg$S[1], f1(sg$s[1])))
      }
    }
  }
  hc_write(out, file.path(dir, "partners.tsv"))
  out
}


### Final call

# Final call per sequence (1204). Precedence: protein > TIR end models and/or end partner > Helitron ends > SINE.
# Confidence: protein as the element rule, 'weak' for calls on one shared domain class (weak_single; the
# alternatives column lists the labels sharing that domain); TIR + partner of the same superfamily 'high'; TIR alone
# its own confidence; partner alone 'low'; TIR and partner disagreeing: the TIR label, 'low'; Helitron ends alone
# 'weak'; SINE its own confidence. Flags are warnings: tandem, low_complexity, ncRNA, host_protein, chimera.
hc_combine <- function(all, P, Ti, H, S, PA, fl, host) {
  Pm <- P[match(all$id, entry)]; Tm <- Ti[match(all$id, id)]; Hm <- H[match(all$id, id)]
  Sm <- S[match(all$id, id)]; PAm <- PA[match(all$id, id)]
  nz <- function(x) !is.na(x) & nzchar(x)
  g <- function(x) fifelse(is.na(x), "", x)
  fit <- function(type, cls) switch(type,
    LINE = "RT" %in% cls, LTR = "RT" %in% cls && ("INT" %in% cls || "RH" %in% cls),
    DIRS = "YR" %in% cls && "RT" %in% cls, PLE = "RT" %in% cls, DNA = "TPase" %in% cls,
    RC = "HEL" %in% cls || "REP" %in% cls, POLINTON = sum(c("ATPase", "POLB", "INT") %in% cls) >= 2,
    CRYPTON = "YR" %in% cls, TRUE)
  tir_label <- function(sf) if (sf %in% c("", "DNA")) "DNA" else if (sf == "Maverick") "POLINTON/Maverick" else paste0("DNA/", sf)
  tand <- split(fl$tandem, by = "id"); lowc <- split(fl$lowc, by = "id"); rfam <- split(fl$rfam, by = "id")
  rows <- vector("list", nrow(all))
  summ <- list()
  for (k in seq_len(nrow(all))) {
    x <- all$id[k]; n <- nchar(all$seq[k])
    o <- setNames(as.list(rep("", length(hc_cols))), hc_cols)
    o$id <- x; o$input_id <- all$input_id[k]; o$length <- as.character(n)
    notes <- character()
    pc <- g(Pm$call[k]); cls <- if (nzchar(pc)) strsplit(Pm$call_classes[k], ",", fixed = TRUE)[[1]] else character()
    ptype <- sub("/.*", "", pc)
    ru <- g(Pm$runner_up[k])
    if (nzchar(pc)) {
      o$domains <- Pm$call_classes[k]; o$prot_score <- Pm$call_score[k]; o$runner_up <- ru
      lead <- as.numeric(Pm$call_score[k]) - (if (nzchar(ru)) as.numeric(sub(".*:", "", ru)) else 0)
      o$prot_lead <- f1(lead)
      o$fits_type <- if (fit(ptype, cls)) "yes" else "no"
      if (nzchar(ru) && as.numeric(o$prot_lead) < hmm_cfg$close_lead) notes <- c(notes, paste0("close_runner_up:", sub(":[^:]*$", "", ru)))
    }
    tir <- identical(Tm$call[k], "TIR")
    tl <- if (tir) tir_label(g(Tm$superfamily[k])) else ""
    if (tir) o$tir <- paste0(tl, ":", Tm$ends[k], ":", Tm$confidence[k], if (nz(Tm$sub[k])) paste0(":", Tm$sub[k]) else "")
    pcall <- g(PAm$call[k]); if (pcall == "ambiguous") pcall <- ""
    if (nz(PAm$ends[k])) {
      o$partner <- paste0(if (nz(PAm$call[k])) PAm$call[k] else "-", ":", PAm$partner[k], ":", PAm$partner_source[k], ":",
                          PAm$ends[k], ":", PAm$score[k], if (nz(PAm$library_partner[k])) paste0(":library=", PAm$library_label[k]) else "")
    }
    hel <- identical(Hm$call[k], "Helitron")
    hl <- if (hel) (if (Hm$group[k] %in% c("Helitron-1", "Helentron")) paste0("RC/", Hm$group[k]) else "RC") else ""
    if (hel) o$helitron <- paste0(hl, ":", Hm$ends[k], ":", Hm$confidence[k])
    sc <- g(Sm$call[k]); if (sc == "none") sc <- ""
    if (nzchar(sc)) o$sine <- paste0(sc, ":", g(Sm$confidence[k]), ":", if (nz(Sm$head_model[k])) Sm$head_model[k] else "-")
    if (nzchar(pc)) {
      key <- if (length(cls) == 1) paste0(ptype, ":", cls[1]) else ""
      lab <- pc; ev <- "protein"
      weak <- key %in% hmm_cfg$weak_single || (paste0(key, ":low") %in% hmm_cfg$weak_single && Pm$confidence[k] == "low")
      conf <- if (weak) "weak" else Pm$confidence[k]
      if (weak) { o$alternatives <- g(hmm_cfg$alt[key]); notes <- c(notes, paste0("single_domain:", sub(".*:", "", key))) }
      if (nz(Pm$lineage[k])) {
        o$sub <- Pm$lineage[k]; o$sub_source <- "protein"
        o$sub_models <- as.character(sum(nzchar(strsplit(Pm$lineages[k], ";", fixed = TRUE)[[1]])))
      }
    } else if (tir || nzchar(pcall)) {
      if (tir && nzchar(pcall)) {
        if (tl %in% c(pcall, "DNA")) { lab <- pcall; ev <- "TIR+partner"; conf <- "high" }
        else { lab <- tl; ev <- "TIR"; conf <- "low"; notes <- c(notes, paste0("partner_conflict:", pcall)) }
      } else if (tir) { lab <- tl; ev <- "TIR"; conf <- Tm$confidence[k] }
      else { lab <- pcall; ev <- "partner"; conf <- "low" }
      sb <- if (tir && nz(Tm$sub[k])) Tm$sub[k] else if (nzchar(pcall) && lab == pcall) g(PAm$partner_sub[k]) else ""
      if (nzchar(sb)) { o$sub <- sb; o$sub_source <- if (tir && nz(Tm$sub[k])) "TIR" else "partner" }
    } else if (hel) { lab <- hl; ev <- "Helitron"; conf <- "weak" }
    else if (startsWith(sc, "SINE/")) { lab <- sc; ev <- "SINE"; conf <- if (nz(Sm$confidence[k])) Sm$confidence[k] else "low" }
    else { lab <- ""; ev <- ""; conf <- "none" }
    o$label <- lab; o$evidence <- ev; o$confidence <- conf
    f <- character()
    tx <- tand[[x]]; if (!is.null(tx) && hc_union(tx$a, tx$b) / n >= hmm_cfg$tandem_frac) f <- c(f, "tandem")
    lx_ <- lowc[[x]]; if (!is.null(lx_) && hc_union(lx_$a, lx_$b) / n >= hmm_cfg$lowcomp_frac) f <- c(f, "low_complexity")
    rh <- rfam[[x]]
    if (!is.null(rh) && nrow(rh)) {
      fams <- sort(unique(rh$fam), method = "radix")
      if (hc_union(rh$a, rh$b) / n >= hmm_cfg$ncrna_frac || any(fams %in% hmm_cfg$rrna)) f <- c(f, paste0("ncRNA:", paste(fams, collapse = "+")))
    }
    hp <- if (!is.null(host)) host[entry == x] else NULL
    phigh <- ev == "protein" && identical(Pm$confidence[k], "high")
    if (!is.null(hp) && nrow(hp) && (hp$lnn >= hmm_cfg$host_aln || !phigh)) f <- c(f, paste0("host_protein:", hp$s))
    if (phigh) {
      if (tir && identical(Tm$confidence[k], "high") && ptype != "DNA") f <- c(f, "chimera:protein+TIR")
      if (hel && ptype != "RC") f <- c(f, "chimera:protein+Helitron")
      if (nzchar(ru)) {
        rl <- sub(":[^:]*$", "", ru); rs <- as.numeric(sub(".*:", "", ru))
        if (sub("/.*", "", rl) != ptype && rs >= hmm_cfg$chimera_runner_bits &&
            rs >= hmm_cfg$chimera_runner_frac * as.numeric(Pm$call_score[k])) f <- c(f, paste0("chimera:runner_up=", rl))
      }
    }
    o$flags <- paste(f, collapse = ","); o$notes <- paste(notes, collapse = ",")
    rows[[k]] <- o
    summ[[k]] <- c(if (nzchar(lab)) lab else "-", conf, if (nzchar(ev)) ev else "-")
  }
  cl <- rbindlist(rows)
  sm <- as.data.table(do.call(rbind, summ)); setnames(sm, c("label", "confidence", "evidence"))
  sm <- sm[, .(n = as.character(.N)), by = .(label, confidence, evidence)]
  setorder(sm, label, confidence, evidence)
  list(classification = cl, summary = sm)
}


### Driver

# Library files (libdir as built by the classifier project's 1201_libdir.sh); gzipped models are decompressed
# once into work
hc_lib <- function(libdir, work) {
  f <- function(...) { p <- file.path(libdir, ...); if (!file.exists(p)) stop("missing library file: ", p); p }
  dir.create(work, showWarnings = FALSE, recursive = TRUE)
  tm <- hc_tsv(f("tir", "pantera_tir_ends.models.tsv"))
  ts <- hc_tsv(f("tir", "pantera_tir_ends.subs.tsv"))
  hm <- hc_tsv(f("hel", "pantera_helitron_ends.models.tsv"))
  host <- file.path(libdir, "flags", "sp_euk_host.dmnd")
  list(prot_hmm = hc_unzip(f("prot", "pantera_v5.cal.hmm.gz"), work),
       ref = hc_unzip(f("prot", "ref.c70.faa.gz"), work),
       tir_hmm = hc_unzip(f("tir", "pantera_tir_ends.hmm.gz"), work),
       tir_sf = setNames(tm$superfamily, tm$model),
       sub_sf = setNames(ts$superfamily, ts$model), sub_sub = setNames(ts$sub, ts$model),
       hel_hmm = hc_unzip(f("hel", "pantera_helitron_ends.hmm.gz"), work),
       hel_end = setNames(hm$end, hm$model), hel_group = setNames(hm$group, hm$model),
       sine_hmm = f("sine", "pantera_sine.hmm"), sine_models = hc_tsv(f("sine", "pantera_sine.models.tsv")),
       heads_cm = f("sine", "heads", "heads.cm"), ncrna_cm = f("flags", "euk_ncrna.cm"),
       host_dmnd = if (file.exists(host)) host else NULL, partners = f("partners", "partners_curated.fna"))
}

# One chunk: every per-sequence piece, resumable ('done' marker)
hc_chunk <- function(ch, outdir, lib) {
  c <- ch$chunk[1]
  dir <- file.path(outdir, "chunks", c)
  if (file.exists(file.path(dir, "done"))) return(c)
  unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
  fna <- file.path(outdir, "chunks", paste0(c, ".fna"))
  hc_write_fasta(ch$id, ch$seq, fna)
  hc_piece_prot(ch, dir, lib)
  hc_piece_tir(ch, dir, lib)
  hc_piece_hel(ch, dir, lib)
  hc_piece_sine(ch, fna, dir, lib)
  hc_piece_flags(ch, fna, dir, lib, file.path(dir, "prot"))
  unlink(file.path(dir, c("tir/windows.fna", "hel/windows.fna")))
  file.create(file.path(dir, "done"))
  c
}

# Classifies the sequences of `input` (FASTA). Writes outdir/classification.tsv and outdir/summary.tsv and
# returns the classification table. Intermediate files stay in outdir/chunks unless keep = FALSE.
hmm_classify <- function(input, outdir, libdir, threads = 1L, chunk = NULL, keep = TRUE) {
  dir.create(file.path(outdir, "chunks"), recursive = TRUE, showWarnings = FALSE)
  if (is.null(chunk)) {        # default: one chunk per thread, between hmm_cfg$chunk_min and hmm_cfg$chunk
    n <- sum(startsWith(readLines(input), ">"))
    chunk <- min(hmm_cfg$chunk, max(hmm_cfg$chunk_min, as.integer(ceiling(n / threads))))
  }
  all <- hc_split(input, chunk)
  hc_log("HMM classifier:", nrow(all), "sequences,", uniqueN(all$chunk), "chunks")
  lib <- hc_lib(libdir, file.path(outdir, "lib"))
  chs <- split(all, by = "chunk")
  res <- mclapply(chs, hc_chunk, outdir = outdir, lib = lib, mc.cores = threads, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, TRUE, "try-error")
  if (any(bad)) stop("HMM classifier failed in chunk ", names(chs)[bad][1], ": ", res[bad][[1]])
  rd <- function(p) rbindlist(lapply(names(chs), function(c) hc_tsv(file.path(outdir, "chunks", c, p))))
  P <- rd("prot/prot_calls.tsv"); Ti <- rd("tir/tir_calls.tsv"); H <- rd("hel/hel_calls.tsv"); S <- rd("sine/calls.tsv")
  hc_log("HMM classifier: end partners")
  PA <- hc_partners(all, P, outdir, lib, threads)
  fl <- lapply(names(chs), function(c) hc_read_flags(file.path(outdir, "chunks", c)))
  fl <- list(tandem = rbindlist(lapply(fl, `[[`, "tandem")), lowc = rbindlist(lapply(fl, `[[`, "lowc")),
             rfam = rbindlist(lapply(fl, `[[`, "rfam")))
  host <- rbindlist(lapply(names(chs), function(c) hc_read_host(file.path(outdir, "chunks", c))))
  if (!nrow(host)) host <- NULL
  r <- hc_combine(all, P, Ti, H, S, PA, fl, host)
  hc_write(r$classification, file.path(outdir, "classification.tsv"))
  hc_write(r$summary, file.path(outdir, "summary.tsv"))
  hc_write(P, file.path(outdir, "prot_calls.tsv")); hc_write(Ti, file.path(outdir, "tir_calls.tsv"))
  hc_write(H, file.path(outdir, "hel_calls.tsv")); hc_write(S, file.path(outdir, "sine_calls.tsv"))
  hc_write(all[, .(id, input_id)], file.path(outdir, "ids.tsv"))
  if (!keep) unlink(file.path(outdir, c("chunks", "partners", "lib")), recursive = TRUE)
  cn <- r$classification[, .N, by = confidence]
  hc_log("HMM classifier:", paste(cn$confidence, cn$N, collapse = ", "))
  r$classification
}


### Command line (when run with Rscript, not when loaded with source())
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(getopt))
  spec <- matrix(c("input", "i", 1, "character", "output", "o", 1, "character", "libdir", "L", 1, "character",
                   "threads", "T", 1, "integer", "chunk", "C", 1, "integer", "clean", "k", 0, "logical",
                   "help", "h", 0, "logical"), byrow = TRUE, ncol = 4)
  opt <- getopt(spec)
  if (!is.null(opt$help) || is.null(opt$input) || is.null(opt$output) || is.null(opt$libdir)) {
    cat("Usage: pantera_classify.R -i library.fa -o outdir -L libdir [-T threads (1)] [-C chunk (500)]",
        "[-k: remove intermediate files]\n")
    q(status = if (is.null(opt$help)) 1 else 0)
  }
  hmm_classify(opt$input, opt$output, opt$libdir, threads = if (is.null(opt$threads)) 1L else opt$threads,
               chunk = opt$chunk, keep = is.null(opt$clean))
}
