# SPEC section 6, corrected specification.
#   global bias     : ~ tissue + rep + copy            -> `copy` is the AVERAGE effect over tissues
#   tissue-specific : LRT, full ~ ... + copy:tissue vs reduced ~ ... , per gene
# With an interaction present the `copy` coefficient is the effect at the REFERENCE tissue,
# not the average, so the two questions need two fits.
suppressMessages(library(DESeq2))
g   <- read.delim("01_counts/gref_counts.tsv", check.names=FALSE)
KEY <- paste(g$AT_gene, g$br_sub, sep="|")
cc  <- grep("\\|", colnames(g), value=TRUE)
pref<- sub("\\|.*","",cc); lib <- sub(".*\\|","",cc)
pick<- function(p,t){ i<-which(pref==p & grepl(paste0(t,"[ABC]$"),lib)); i[order(lib[i])] }
MINC <- 20
SP <- list(list(id="napus",a1="Bna_A",a2="Bna_C",t1="Pol",t2="Stig"),
           list(id="juncea",a1="Bju_A",a2="Bju_B",t1="Pol",t2="Stig"),
           list(id="carinata",a1="Bca_B",a2="Bca_C",t1="Pol",t2="L_Stig"))
dir.create("04_model", showWarnings=FALSE); out <- list()
cat(sprintf("%-10s %7s | %-34s | %s\n","species","pairs","GLOBAL bias (average over tissues)","bias differs by tissue"))
for (s in SP) {
  A <- as.matrix(g[,cc[c(pick(s$a1,s$t1),pick(s$a1,s$t2))]])
  B <- as.matrix(g[,cc[c(pick(s$a2,s$t1),pick(s$a2,s$t2))]])
  keep <- complete.cases(A,B) & (rowSums(A)+rowSums(B) >= MINC)
  M <- cbind(A[keep,],B[keep,]); kk <- KEY[keep]
  cd <- data.frame(copy=factor(rep(c("c1","c2"),each=6),levels=c("c2","c1")),
                   tissue=factor(rep(rep(c("pollen","stigma"),each=3),2),levels=c("pollen","stigma")),
                   rep=factor(rep(rep(1:3,2),2)))
  rownames(cd)<-colnames(M)<-paste0(cd$copy,"_",cd$tissue,"_",cd$rep)
  # 1. global: no interaction
  d1 <- DESeqDataSetFromMatrix(round(M), cd, ~ tissue + rep + copy)
  sizeFactors(d1) <- rep(1, ncol(M)); d1 <- DESeq(d1, quiet=TRUE)
  r1 <- results(d1, name="copy_c1_vs_c2", alpha=0.05)
  up<-sum(!is.na(r1$padj)&r1$padj<=0.05&r1$log2FoldChange>0)
  dn<-sum(!is.na(r1$padj)&r1$padj<=0.05&r1$log2FoldChange<0)
  # 2. does the bias depend on tissue? LRT
  d2 <- DESeqDataSetFromMatrix(round(M), cd, ~ tissue + rep + copy + copy:tissue)
  sizeFactors(d2) <- rep(1, ncol(M))
  d2 <- DESeq(d2, test="LRT", reduced = ~ tissue + rep + copy, quiet=TRUE)
  r2 <- results(d2, alpha=0.05)
  ni <- sum(!is.na(r2$padj) & r2$padj <= 0.05)
  cat(sprintf("%-10s %7d | %s %5d : %s %5d  ratio %5.2f  median lfc %+.3f | %5d (%4.1f%%)\n",
      s$id, nrow(M), s$a1, up, s$a2, dn, up/dn, median(r1$log2FoldChange,na.rm=TRUE), ni, 100*ni/nrow(M)))
  write.csv(data.frame(key=kk, global_lfc=r1$log2FoldChange, global_padj=r1$padj,
                       tissue_lrt_padj=r2$padj),
            file.path("04_model", paste0("joint_", s$id, ".csv")), row.names=FALSE)
  out[[s$id]] <- data.frame(species=s$id, pairs=nrow(M), copy1=s$a1, copy2=s$a2,
      global_toward1=up, global_toward2=dn, global_ratio=round(up/dn,3),
      median_lfc=round(median(r1$log2FoldChange,na.rm=TRUE),4),
      tissue_dependent=ni, tissue_dependent_pct=round(100*ni/nrow(M),1))
}
write.csv(do.call(rbind,out), "04_model/joint_summary.csv", row.names=FALSE)
