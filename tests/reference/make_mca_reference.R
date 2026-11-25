# make_mca_reference.R (2026-11-25)
# Purpose: reference p-values from R jaccard::jaccard.test.mca (ncchung/jaccard@ad1d8fc, as used by CSI-GEP) for the
#          Python port. Cases = what CSI-GEP produces: x, y binary over m genes with nx, ny ones and c shared.
#          Equal sizes (nx = ny = JL-1 for the JL grid) at every overlap, plus a few unequal-size cases.
# Output: tests/reference/mca_reference.tsv (m, nx, ny, c, statistics, pvalue, expectation)
# Run   : Rscript --vanilla tests/reference/make_mca_reference.R   (inside csi-gep_r:2026.11; see make_mca_reference.bsub)
suppressPackageStartupMessages(library(jaccard))
mk <- function(m, nx, ny, c) { x <- y <- integer(m); x[seq_len(nx)] <- 1L; y[seq_len(c)] <- 1L
  if (ny > c) y[nx + seq_len(ny - c)] <- 1L; list(x = x, y = y) }
cases <- list()
for (m in c(2000, 5100)) for (n in c(2:9, seq(19, 99, 10))) for (c in 0:n) cases[[length(cases) + 1]] <- c(m, n, n, c)
for (v in list(c(2000, 9, 8, 3), c(2000, 19, 21, 5), c(5100, 49, 50, 10), c(5100, 99, 98, 2), c(500, 30, 40, 12), c(100, 10, 10, 4)))
  cases[[length(cases) + 1]] <- v
res <- do.call(rbind, lapply(cases, function(v) { xy <- mk(v[1], v[2], v[3], v[4])
  r <- suppressWarnings(jaccard.test.mca(xy$x, xy$y))
  data.frame(m = v[1], nx = v[2], ny = v[3], c = v[4], statistics = sprintf("%.17g", r$statistics),
             pvalue = sprintf("%.17g", r$pvalue), expectation = sprintf("%.17g", r$expectation)) }))
out <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "mca_reference.tsv")
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("STATUS: OK | cases=%d | pvalue==0: %d | out=%s\n", nrow(res), sum(as.numeric(res$pvalue) == 0), out))
