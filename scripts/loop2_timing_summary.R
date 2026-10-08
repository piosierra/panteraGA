#!/usr/bin/env Rscript
# Summarises the loop 2 timing files that panteraGA writes when run with -d.
#   Rscript scripts/loop2_timing_summary.R <output folder> [threads used]
# Loop 2 has two steps: clustering per length window (loop2_window_timing.tsv),
# then alignment and consensus per cluster (loop2_mafft_timing.tsv).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
if (length(args) < 1) stop("Usage: loop2_timing_summary.R <output folder> [threads used]")
w <- fread(file.path(args[1], "loop2_window_timing.tsv"))
m <- fread(file.path(args[1], "loop2_mafft_timing.tsv"))
threads <- if (length(args) > 1) as.numeric(args[2]) else NA

iv <- rbind(w[, .(t_start, t_end)], m[, .(t_start, t_end)])
wall <- max(iv$t_end) - min(iv$t_start)
busy <- sum(iv$t_end - iv$t_start)
cat(sprintf("Loop 2 wall time: %.0f s (%.1f min)\n", wall, wall / 60))
cat(sprintf("  clustering step: %.0f s; alignment step: %.0f s\n",
            max(w$t_end) - min(w$t_start), max(m$t_end) - min(m$t_start)))
cat(sprintf("Sum of task times: %.0f s; average busy workers: %.1f%s\n",
            busy, busy / wall, if (is.na(threads)) "" else sprintf(" of %g", threads)))
cat(sprintf("Of the task time: mafft %.0f%%, cd-hit %.0f%%, rest (consensus, TSD...) %.0f%%\n",
            100 * sum(m$seconds, na.rm = TRUE) / busy, 100 * sum(w$cdhit_s) / busy,
            100 * (busy - sum(m$seconds, na.rm = TRUE) - sum(w$cdhit_s)) / busy))

# How many tasks were running over time, and the time spent with few running
ev <- rbind(data.table(t = iv$t_start, d = 1), data.table(t = iv$t_end, d = -1))[order(t, -d)]
ev[, active := cumsum(d)]
ev[, dt := shift(t, type = "lead") - t]
ev <- ev[!is.na(dt)]
for (k in c(1, 2, 4)) {
  cat(sprintf("Time with at most %d task(s) running: %.0f s (%.0f%% of loop 2)\n",
              k, ev[active <= k, sum(dt)], 100 * ev[active <= k, sum(dt)] / wall))
}

cat(sprintf("\nAlignments: %d; mafft total %.0f s\n", nrow(m), sum(m$seconds, na.rm = TRUE)))
setorder(m, -seconds)
top <- head(m, 10)
cat(sprintf("The 10 slowest take %.0f%% of all mafft time:\n",
            100 * sum(top$seconds, na.rm = TRUE) / sum(m$seconds, na.rm = TRUE)))
print(top[, .(win_start, n_seqs, mean_len, aln_len, seconds, total_s)], row.names = FALSE)

cat("\nSlowest windows to cluster:\n")
setorder(w, -total_s)
print(head(w[, .(win_start, n_segments, n_clusters, cdhit_s, total_s)], 5), row.names = FALSE)
