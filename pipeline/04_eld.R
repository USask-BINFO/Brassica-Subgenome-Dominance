# ELD per SPEC.md section 4 (Rapp/Udall/Wendel framework), with equivalence testing.
#
#  total pair expression in the allotetraploid vs each diploid progenitor.
#  Step 1  parents: PEC (P1 equals P2) or PED (P1 differs from P2)
#  Step 2  PED -> equals P1 and differs from P2 = ELD toward P1
#               equals P2 and differs from P1 = ELD toward P2
#               differs from both, lies between = additive
#               outside the parental range        = transgressive
#          PEC -> equals both = unchanged; above both / below both = transgressive
#
# "equals" is an EQUIVALENCE TEST (lfcThreshold + lessAbs), not a failed difference test.
# That is what makes additive reachable; defining equality as non-significance makes it
# unreachable at a large threshold, which is the flaw the manuscript already identified.
# Mid-parent is NOT used as a pseudo-library: the parental comparison is done directly and
# deterministically, so there is no sample() randomness.
suppressMessages(library(DESeq2))
args  <- commandArgs(trailingOnly = TRUE)
THETA <- as.numeric(args[1]); ALPHA <- as.numeric(args[2]); MINC <- as.numeric(args[3])
OUT   <- args[4]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

g   <- read.delim("01_counts/gref_counts.tsv", check.names = FALSE)
key <- paste(g$AT_gene, g$br_sub, sep = "|")
cc  <- grep("\\|", colnames(g), value = TRUE)
pref <- sub("\\|.*", "", cc); lib <- sub(".*\\|", "", cc)
pick <- function(p, tis) { i <- which(pref == p & grepl(paste0(tis, "[ABC]$"), lib)); i[order(lib[i])] }

SAMPLES <- list(
  list(id="napus_pollen",     a1="Bna_A",a2="Bna_C",p1="Bra_A",p2="Bol_C",at="Pol",   pt="Pol"),
  list(id="napus_stigma",     a1="Bna_A",a2="Bna_C",p1="Bra_A",p2="Bol_C",at="Stig",  pt="Stig"),
  list(id="carinata_pollen",  a1="Bca_B",a2="Bca_C",p1="Bni_B",p2="Bol_C",at="Pol",   pt="Pol"),
  list(id="carinata_stigmaE", a1="Bca_B",a2="Bca_C",p1="Bni_B",p2="Bol_C",at="E_Stig",pt="Stig"),
  list(id="carinata_stigmaL", a1="Bca_B",a2="Bca_C",p1="Bni_B",p2="Bol_C",at="L_Stig",pt="Stig"),
  list(id="juncea_pollen",    a1="Bju_A",a2="Bju_B",p1="Bra_A",p2="Bni_B",at="Pol",   pt="Pol"),
  list(id="juncea_stigma",    a1="Bju_A",a2="Bju_B",p1="Bra_A",p2="Bni_B",at="Stig",  pt="Stig"))

# one contrast on a shared fit: returns "up"/"down"/"equal"/"indet" per gene
contrast3 <- function(dds, a, b) {
  up <- results(dds, contrast=c("grp",a,b), lfcThreshold=THETA, altHypothesis="greaterAbs", alpha=ALPHA)
  eq <- results(dds, contrast=c("grp",a,b), lfcThreshold=THETA, altHypothesis="lessAbs",    alpha=ALPHA)
  s <- rep("indet", nrow(up))
  sg <- !is.na(up$padj) & up$padj <= ALPHA
  s[sg & up$log2FoldChange > 0] <- "up"
  s[sg & up$log2FoldChange < 0] <- "down"
  s[!sg & !is.na(eq$padj) & eq$padj <= ALPHA] <- "equal"
  s
}

res <- list()
for (s in SAMPLES) {
  ia1<-pick(s$a1,s$at); ia2<-pick(s$a2,s$at); ip1<-pick(s$p1,s$pt); ip2<-pick(s$p2,s$pt)
  A1<-as.matrix(g[,cc[ia1]]); A2<-as.matrix(g[,cc[ia2]])
  P1<-as.matrix(g[,cc[ip1]]); P2<-as.matrix(g[,cc[ip2]])
  keep <- complete.cases(A1,A2,P1,P2)
  F1 <- A1[keep,] + A2[keep,]                      # allotetraploid TOTAL for the pair
  p1 <- P1[keep,]; p2 <- P2[keep,]
  keep2 <- rowSums(F1)+rowSums(p1)+rowSums(p2) >= MINC
  F1<-F1[keep2,]; p1<-p1[keep2,]; p2<-p2[keep2,]
  M <- cbind(F1,p1,p2); colnames(M) <- paste0(rep(c("F1","P1","P2"),each=3),"_",1:3)
  cd <- data.frame(grp=factor(rep(c("F1","P1","P2"),each=3), levels=c("F1","P1","P2")),
                   row.names=colnames(M))
  d <- DESeq(DESeqDataSetFromMatrix(round(M), cd, ~ grp), quiet=TRUE)
  vp  <- contrast3(d,"P1","P2")   # parents differ?
  v1  <- contrast3(d,"F1","P1")
  v2  <- contrast3(d,"F1","P2")
  cls <- rep("unclassified", length(vp))
  pec <- vp=="equal"; ped <- vp %in% c("up","down")
  cls[pec & v1=="equal" & v2=="equal"] <- "PEC_unchanged"
  cls[pec & v1=="up"    & v2=="up"   ] <- "PEC_transgressive_up"
  cls[pec & v1=="down"  & v2=="down" ] <- "PEC_transgressive_down"
  cls[ped & v1=="equal" & v2!="equal"] <- "PED_ELD_P1"
  cls[ped & v2=="equal" & v1!="equal"] <- "PED_ELD_P2"
  cls[ped & v1=="up"    & v2=="up"   ] <- "PED_transgressive_up"
  cls[ped & v1=="down"  & v2=="down" ] <- "PED_transgressive_down"
  cls[ped & ((v1=="up" & v2=="down")|(v1=="down" & v2=="up"))] <- "PED_additive"
  df <- data.frame(key=key[keep][keep2], parents=vp, f1_vs_p1=v1, f1_vs_p2=v2, class=cls)
  write.csv(df, file.path(OUT,paste0("eld_",s$id,".csv")), row.names=FALSE)
  t <- table(factor(cls, levels=c("PED_ELD_P1","PED_ELD_P2","PED_additive","PED_transgressive_up",
       "PED_transgressive_down","PEC_unchanged","PEC_transgressive_up","PEC_transgressive_down","unclassified")))
  pedn <- sum(t[1:5])
  cat(sprintf("%-18s n=%6d | PED %5d: ELD_P1 %4d (%4.1f%%) ELD_P2 %4d (%4.1f%%) additive %4d TRE %4d | PEC %5d | unclassified %5d\n",
      s$id, nrow(M), pedn, t[1], 100*t[1]/max(pedn,1), t[2], 100*t[2]/max(pedn,1), t[3], t[4]+t[5], sum(t[6:8]), t[9]))
  res[[s$id]] <- data.frame(sample=s$id, tested=nrow(M), PED=pedn, ELD_P1=t[1], ELD_P2=t[2],
      additive=t[3], PED_TRE=t[4]+t[5], PEC=sum(t[6:8]), unclassified=t[9],
      ELD_P1_pct=round(100*t[1]/max(pedn,1),2), ELD_P2_pct=round(100*t[2]/max(pedn,1),2))
}
write.csv(do.call(rbind,res), file.path(OUT,"eld_summary.csv"), row.names=FALSE)
cat("\nwritten to", OUT, "\n")
