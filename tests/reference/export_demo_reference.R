# export_demo_reference.R (2026-11-25; generic run since 2026-11-25)
# Purpose: R reference for the Python port, from a demo run made with the R scripts (tutorial_run_img2026.11):
#          Jaccard step (significant pairs, community number per K x JL; all p-values at one K), selection intermediates
#          (loess curve, elbow, slope per JL, re-running GEP_analysis.R's own code), final outputs (rescue 0, and rescue 1
#          via a fresh GEP_analysis.R run).
# Usage : Rscript --vanilla export_demo_reference.R <run_dir/> <out_dir/> <K csv> <K for p-values, or 0> [name] [fresh 0/1]
#         fresh=1: both rescue 0 and 1 are fresh GEP_analysis.R runs on <K csv> (else rescue0 = <run_dir>GEP_results/)
# Output: <out_dir>/{sig_pairs.tsv, communities.tsv, pvalues_K<K>.tsv, selection.tsv, curves.tsv,
#         rescue<0|1>/{params.tsv, gepindex.tsv, community_geps.tsv}}
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(elbow)})
a <- commandArgs(TRUE); run <- a[1]; out <- a[2]; ranks <- scan(text = a[3], sep = ",", quiet = TRUE); kp <- as.numeric(a[4])
nm <- if (length(a) >= 5) a[5] else "demo"; fresh <- length(a) >= 6 && a[6] == "1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
jls <- seq(10, 100, 10); sp <- list(); cm <- list(); pv <- list()
for (k in ranks) { R <- readRDS(sprintf("%sresults/jaccardtest_results_K%d.RDS", run, k))
  for (jl in jls) { r <- R[[as.character(jl)]]; s <- as.data.frame(r$sigFwerJaccard)
    if (nrow(s)) sp[[length(sp) + 1]] <- data.frame(K = k, JL = jl, row = s$row, col = s$col, ord = seq_len(nrow(s)))
    cm[[length(cm) + 1]] <- data.frame(K = k, JL = jl, community_number = r$jaccard.com, n_sig = nrow(s))
    if (k == kp) pv[[length(pv) + 1]] <- data.frame(JL = jl, i = rep(seq_len(k), each = k), j = rep(seq_len(k), k),
                                                    p = sprintf("%.17g", r$jaccardIndPs)) } }
w <- function(x, f) write.table(x, paste0(out, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(do.call(rbind, sp), "sig_pairs.tsv"); w(do.call(rbind, cm), "communities.tsv"); if (length(pv)) w(do.call(rbind, pv), sprintf("pvalues_K%d.tsv", kp))
com.df <- do.call(rbind, cm)[, c("community_number", "K", "JL")]; colnames(com.df) <- c("Community_Number", "Rank", "Jaccard_Length")
sel <- list(); cur <- list()
for (x in jls) { jl.df <- com.df %>% filter(Jaccard_Length == x)          # GEP_analysis.R selection block
  g <- ggplot_build(jl.df %>% ggplot(aes(Rank, Community_Number)) + geom_point() + geom_smooth(method = "loess", span = 0.3, method.args = list(degree = 1)))
  e <- data.frame(X = g$data[[2]]$x, Y = g$data[[2]]$y); ip <- elbow(data = e, plot = FALSE)
  f <- ip$data %>% filter(X >= ip$X_selected); reg <- lm(f$benefits ~ f$X)
  sel[[length(sel) + 1]] <- data.frame(JL = x, Rank = ranks[which.min(abs(ranks - ip$X_selected))], Slope = sprintf("%.17g", reg$coefficients[2]),
                                       X_selected = sprintf("%.17g", ip$X_selected))
  cur[[length(cur) + 1]] <- data.frame(JL = x, X = sprintf("%.17g", e$X), Y = sprintf("%.17g", e$Y), benefits = sprintf("%.17g", ip$data$benefits)) }
w(do.call(rbind, sel), "selection.tsv"); w(do.call(rbind, cur), "curves.tsv")
for (rs in 0:1) { d <- sprintf("%srescue%d/", out, rs); dir.create(d, showWarnings = FALSE)
  src <- if (rs == 0 && !fresh) paste0(run, "GEP_results/") else d
  if (src == d) system(sprintf("Rscript --vanilla /app/GEP_analysis.R %s %s %s %s %d > %sGEP_analysis.log 2>&1", run, nm, d, paste(ranks, collapse = ","), rs, d))
  p <- readRDS(paste0(src, nm, "_optimized_parameters.RDS")); w(p, sprintf("rescue%d/params.tsv", rs))
  gi <- readRDS(paste0(src, nm, "_gepindex.RDS"))
  w(do.call(rbind, lapply(seq_along(gi), function(i) data.frame(GEP = i, split = c(rep(0, length(gi[[i]]$Split0)), rep(1, length(gi[[i]]$Split1))),
      program = c(gi[[i]]$Split0, gi[[i]]$Split1), pos = c(seq_along(gi[[i]]$Split0), seq_along(gi[[i]]$Split1))))), sprintf("rescue%d/gepindex.tsv", rs))
  m <- readRDS(paste0(src, nm, "_community_geps.RDS")); m <- data.frame(gene = rownames(m), apply(m, 2, function(v) sprintf("%.17g", v)))
  w(m, sprintf("rescue%d/community_geps.tsv", rs))
  if (src != d) file.copy(paste0(src, "output_summary.txt"), sprintf("%srescue%d/output_summary.txt", out, rs), overwrite = TRUE) }
cat(sprintf("STATUS: OK | K=%s | out=%s\n", a[3], out))
