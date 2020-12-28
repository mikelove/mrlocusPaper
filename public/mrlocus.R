cmd_args=commandArgs(TRUE)

dir <- cmd_args[1]
tsv.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
out.filename <- cmd_args[4]

set.seed(1)

devtools::load_all("../../mrlocus")
source("common.R") # common function for mrlocus and ecaviar-mrlocus
out <- getTrimmedSumStats(dir, tsv.filename, ld.filename, r2_threshold=0.05)
sum_stat <- out$sum_stat
ld_mat <- out$ld_mat

out1 <- collapseHighCorSNPs(sum_stat, ld_mat, plot=FALSE)
a2_plink <- "Major_plink"

out2 <- flipAllelesAndGather(out1$sum_stat, out1$ld_mat,
                             out1$ld_mat2,
                             a="eQTL", b="GWAS",
                             ref="Ref", eff="Effect",
                             beta="beta", se="se",
                             a2_plink=a2_plink,
                             snp_id="SNP", sep="_",
                             ab_last=TRUE, plot=FALSE)

png(file=file.path(dir,"initial.png"))
plotInitEstimates(out2)
dev.off()

library(Matrix)
Sigma_npd <- out2$Sigma
for (j in seq_along(out2$Sigma)) {
  Sigma_npd[[j]] <- as.matrix(nearPD(out2$Sigma[[j]])$mat)
}
Sigma_npd2 <- Sigma_npd

nsnp <- lengths(out2$beta_hat_a)

options(mc.cores=4)
beta_hat_a <- list()
beta_hat_b <- list()
for (j in seq_along(nsnp)) {
  print(j)
  if (nsnp[j] > 1) {
    cap.out <- capture.output({ 
      fit <- fitBetaColoc(beta_hat_a=out2$beta_hat_a[[j]],
                          beta_hat_b=out2$beta_hat_b[[j]],
                          se_a=out2$se_a[[j]],
                          se_b=out2$se_b[[j]],
                          Sigma_a=Sigma_npd[[j]],
                          Sigma_b=Sigma_npd2[[j]],
                          verbose=FALSE,
                          open_progress=FALSE,
                          show_messages=FALSE,
                          refresh=-1)
    })
    beta_hat_a[[j]] <- fit$beta_hat_a
    beta_hat_b[[j]] <- fit$beta_hat_b
  } else {
    beta_hat_a[[j]] <- unname(out2$beta_hat_a[[j]])
    beta_hat_b[[j]] <- unname(out2$beta_hat_b[[j]])
  }
}

# make a results list for slope fitting
res <- list(beta_hat_a=beta_hat_a,
            beta_hat_b=beta_hat_b,
            sd_a=out2$se_a,
            sd_b=out2$se_b,
            alleles=out2$alleles)

res <- extractForSlope(res, plot=FALSE)

# remove any clusters that were below p threshold
z.thr <- qnorm(.001/2, lower.tail=FALSE)
z.idx <- res$beta_hat_a/res$sd_a > z.thr
res <- lapply(res, function(x) {
  if (is.numeric(x)) {
    x[z.idx]
  } else {
    x[z.idx,]
  }
})

# second round trimming based on candidate SNPs
load("LDmatrix_allSNPs.rda")
r2 <-  get(paste0("LD_", dir))
idx <- sapply(res$alleles$id, function(snp) grep(snp, colnames(r2)))
r2 <- r2[idx,idx]
trim_clusters <- trimClusters(r2, r2_threshold=0.05)
trim2.filename <- file.path(dir, paste0(dir, "_mrlocus.trim2"))
if (length(trim_clusters) > 0) {
  res <- lapply(res, function(x) {
    if (is.numeric(x)) {
      x[-trim_clusters]
    } else {
      x[-trim_clusters,,drop=FALSE]
    }
  })
  write(trim_clusters, file=trim2.filename, ncolumns=length(trim_clusters))
} else {
  write(NA, file=trim2.filename)
}

write.table(res$alleles, file.path(dir,"mrlocus_alleles.txt"), row.names=FALSE, quote=FALSE)

res <- fitSlope(res, iter=10000)

save(res, file=out.filename)

sessionInfo()
