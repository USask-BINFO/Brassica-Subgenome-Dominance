#!/usr/bin/env Rscript
# Is homoeolog expression bias INHERITED from the diploid progenitors, or CREATED in the hybrid?
#
# This is the question the restructured story turns on, and the paper currently answers it only
# as a single percentage. The pairs_*.csv files already carry both numbers per gene pair:
#   allo_lfc = log2(copy1/copy2) measured WITHIN the allotetraploid
#   par_lfc  = log2(progenitor1/progenitor2) for the same two genes in the diploids
# so the inheritance question is a direct per-gene comparison, no new quantification needed.
#
# Reported per species-tissue:
#   1. Spearman and Pearson of allo_lfc vs par_lfc  - how much bias tracks the parents
#   2. slope of allo_lfc ~ par_lfc                  - <1 means bias is DAMPENED in the hybrid
#   3. the Rapp/Udall-style 4-way classification    - inherited / novel / lost / reversed
#   4. variance partition                           - share of bias variance explained by parents
suppressMessages({library(stats)})
IN <- "02_heb/theta2_a05"; TH <- 1          # 2-fold on log2 scale
f <- list.files(IN, "^pairs_.*csv$", full.names=TRUE)
out <- list()
cat(sprintf("%-10s %-9s %7s | %6s %6s %6s | %7s %7s %7s %7s | %6s\n",
            "species","tissue","n","rho","r","R2","inherit","novel","lost","revers","slope"))
cat(strrep("-", 104), "\n")
for (p in f) {
  nm <- sub("^pairs_","",sub("\\.csv$","",basename(p)))
  sp <- sub("_.*$","",nm); ti <- sub("^[^_]*_","",nm)
  d <- read.csv(p, stringsAsFactors=FALSE)
  d <- d[is.finite(d$allo_lfc) & is.finite(d$par_lfc), ]
  if (nrow(d) < 50) next
  rho <- suppressWarnings(cor(d$allo_lfc, d$par_lfc, method="spearman"))
  r   <- cor(d$allo_lfc, d$par_lfc)
  fit <- lm(allo_lfc ~ par_lfc, data=d)
  sl  <- unname(coef(fit)[2]); R2 <- summary(fit)$r.squared
  aB <- abs(d$allo_lfc) >= TH; pB <- abs(d$par_lfc) >= TH
  same <- sign(d$allo_lfc) == sign(d$par_lfc)
  inherit <- sum(aB & pB & same); revers <- sum(aB & pB & !same)
  novel   <- sum(aB & !pB);       lost   <- sum(!aB & pB)
  tot <- nrow(d)
  cat(sprintf("%-10s %-9s %7d | %6.3f %6.3f %6.3f | %6.1f%% %6.1f%% %6.1f%% %6.1f%% | %6.3f\n",
      sp, ti, tot, rho, r, R2,
      100*inherit/tot, 100*novel/tot, 100*lost/tot, 100*revers/tot, sl))
  out[[nm]] <- data.frame(species=sp, tissue=ti, n=tot, rho=rho, r=r, R2=R2, slope=sl,
                          pct_inherited=100*inherit/tot, pct_novel=100*novel/tot,
                          pct_lost=100*lost/tot, pct_reversed=100*revers/tot)
}
res <- do.call(rbind, out)
write.csv(res, "05_validation/inherited_vs_novel.csv", row.names=FALSE)
cat("\nInterpretation guide\n")
cat("  R2            share of allotetraploid bias variance explained by the parental difference\n")
cat("  slope < 1     parental bias is DAMPENED in the hybrid; > 1 amplified\n")
cat("  inherited     biased in the allotetraploid AND in the parents, same direction\n")
cat("  novel         biased in the allotetraploid but NOT in the parents (created by hybridization)\n")
cat("  lost          biased in the parents but NOT in the allotetraploid (erased by hybridization)\n")
cat("  reversed      biased in both but in OPPOSITE directions\n")
cat("\nWrote 05_validation/inherited_vs_novel.csv\n")
