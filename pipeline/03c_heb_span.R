# HEB per SPEC.md sections 2 and 3.
#  - within-allotetraploid contrast: both copies come from ONE library, so library size is
#    common and cancels -> size factors held at 1, paired design ~ replicate + copy
#  - progenitor contrast spans species/libraries -> median-of-ratios size factors
#  - equivalence testing: biased / balanced / indeterminate, never "balanced by default"
suppressMessages(library(DESeq2))
args   <- commandArgs(trailingOnly = TRUE)
THETA  <- as.numeric(args[1]); ALPHA <- as.numeric(args[2]); MINC <- as.numeric(args[3])
OUT    <- args[4]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

g <- read.delim("01_counts/gref_counts_span.tsv", check.names = FALSE)
key <- paste(g$AT_gene, g$br_sub, sep = "|")
cc  <- grep("\\|", colnames(g), value = TRUE)
pref <- sub("\\|.*", "", cc); lib <- sub(".*\\|", "", cc)

SAMPLES <- list(
  list(id="napus_pollen",     a1="Bna_A", a2="Bna_C", p1="Bra_A", p2="Bol_C", at="Pol",     pt="Pol"),
  list(id="napus_stigma",     a1="Bna_A", a2="Bna_C", p1="Bra_A", p2="Bol_C", at="Stig",    pt="Stig"),
  list(id="carinata_pollen",  a1="Bca_B", a2="Bca_C", p1="Bni_B", p2="Bol_C", at="Pol",     pt="Pol"),
  list(id="carinata_stigmaE", a1="Bca_B", a2="Bca_C", p1="Bni_B", p2="Bol_C", at="E_Stig",  pt="Stig"),
  list(id="carinata_stigmaL", a1="Bca_B", a2="Bca_C", p1="Bni_B", p2="Bol_C", at="L_Stig",  pt="Stig"),
  list(id="juncea_pollen",    a1="Bju_A", a2="Bju_B", p1="Bra_A", p2="Bni_B", at="Pol",     pt="Pol"),
  list(id="juncea_stigma",    a1="Bju_A", a2="Bju_B", p1="Bra_A", p2="Bni_B", at="Stig",    pt="Stig"))

pick <- function(p, tis) {
  i <- which(pref == p & grepl(paste0(tis, "[ABC]$"), lib))
  i[order(lib[i])]
}
# classify one contrast into biased-1 / biased-2 / balanced / indeterminate
classify <- function(m, grp, equal_sf) {
  cd <- data.frame(grp = factor(grp, levels = unique(grp)),
                   rep = factor(rep(seq_len(ncol(m)/2), 2)))
  d <- DESeqDataSetFromMatrix(round(m), cd, ~ rep + grp)
  if (equal_sf) sizeFactors(d) <- rep(1, ncol(m)) else d <- estimateSizeFactors(d)
  d <- DESeq(d, quiet = TRUE, fitType = "parametric")
  lv <- levels(cd$grp)
  up <- results(d, contrast=c("grp",lv[1],lv[2]), lfcThreshold=THETA, altHypothesis="greaterAbs", alpha=ALPHA)
  eq <- results(d, contrast=c("grp",lv[1],lv[2]), lfcThreshold=THETA, altHypothesis="lessAbs",    alpha=ALPHA)
  st <- rep("indeterminate", nrow(m))
  sig <- !is.na(up$padj) & up$padj <= ALPHA
  st[sig & up$log2FoldChange > 0] <- "toward1"
  st[sig & up$log2FoldChange < 0] <- "toward2"
  st[!sig & !is.na(eq$padj) & eq$padj <= ALPHA] <- "balanced"
  list(status = st, lfc = up$log2FoldChange)
}

res <- list()
for (s in SAMPLES) {
  ia1 <- pick(s$a1, s$at); ia2 <- pick(s$a2, s$at)
  ip1 <- pick(s$p1, s$pt); ip2 <- pick(s$p2, s$pt)
  stopifnot(length(ia1)==3, length(ia2)==3, length(ip1)==3, length(ip2)==3)
  A1 <- as.matrix(g[, cc[ia1]]); A2 <- as.matrix(g[, cc[ia2]])
  P1 <- as.matrix(g[, cc[ip1]]); P2 <- as.matrix(g[, cc[ip2]])
  keep <- complete.cases(A1,A2,P1,P2) & (rowSums(A1)+rowSums(A2) >= MINC) & (rowSums(P1)+rowSums(P2) >= MINC)
  allo <- classify(cbind(A1,A2)[keep,], rep(c("s1","s2"), each=3), TRUE)
  par  <- classify(cbind(P1,P2)[keep,], rep(c("p1","p2"), each=3), FALSE)
  df <- data.frame(key = key[keep], allo = allo$status, allo_lfc = allo$lfc,
                   par = par$status, par_lfc = par$lfc)
  write.csv(df, file.path(OUT, paste0("pairs_", s$id, ".csv")), row.names = FALSE)
  a <- table(factor(df$allo, levels=c("toward1","toward2","balanced","indeterminate")))
  r <- a["toward1"]/a["toward2"]
  cat(sprintf("%-18s tested %6d | toward %s %5d : toward %s %5d  ratio %5.2f | balanced %5d | indeterminate %5d\n",
      s$id, sum(keep), s$a1, a["toward1"], s$a2, a["toward2"], r, a["balanced"], a["indeterminate"]))
  res[[s$id]] <- data.frame(sample=s$id, tested=sum(keep), toward1=a["toward1"], toward2=a["toward2"],
                            ratio=round(r,3), balanced=a["balanced"], indeterminate=a["indeterminate"])
}
write.csv(do.call(rbind,res), file.path(OUT,"heb_summary.csv"), row.names=FALSE)
cat("\nwritten to", OUT, "\n")
