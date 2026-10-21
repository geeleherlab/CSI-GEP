# Tutorial step 7: compare the GEPs CSI-GEP found with the simulation's 41 true programs.
# For each GEP, take its top N genes (highest community score) and test overlap with each true
# program's gene set (hypergeometric, Bonferroni over all GEP x program pairs); report the best match.
# Usage: Rscript check_ground_truth.R <name>_community_geps.RDS ground_truth_programs.tsv <out.tsv> [top_n]
args <- commandArgs(trailingOnly = TRUE)
geps <- readRDS(args[1])                       # genes x GEPs
truth <- read.delim(args[2], stringsAsFactors = FALSE)
out <- args[3]
top_n <- if (length(args) > 3) as.integer(args[4]) else 50

universe <- rownames(geps)
truth <- truth[truth$gene %in% universe, ]
sets <- split(truth$gene, truth$program)
n_tests <- ncol(geps) * length(sets)

res <- do.call(rbind, lapply(seq_len(ncol(geps)), function(i) {
  top <- names(sort(geps[, i], decreasing = TRUE))[seq_len(top_n)]
  ov <- sapply(sets, function(s) sum(top %in% s))
  p <- sapply(names(sets), function(nm)
    phyper(unname(ov[nm]) - 1, length(sets[[nm]]), length(universe) - length(sets[[nm]]), top_n, lower.tail = FALSE))
  best <- names(which.max(ov))
  data.frame(GEP = i, best_program = best, overlap = unname(ov[best]), top_n = top_n,
             p_bonf = min(1, unname(p[best]) * n_tests))
}))
res$significant <- res$p_bonf < 0.05
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)

matched <- unique(res$best_program[res$significant])
cat(sprintf("GEPs found: %d | significantly matched to a true program: %d | true programs recovered: %d of %d\n",
            ncol(geps), sum(res$significant), length(matched), length(sets)))
dup <- table(res$best_program[res$significant])
if (any(dup > 1)) cat("Programs matched by more than one GEP:", paste(names(dup)[dup > 1], collapse = ", "), "\n")
missed <- setdiff(names(sets), matched)
if (length(missed)) cat("True programs not recovered:", paste(missed, collapse = ", "), "\n")
