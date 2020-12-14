cmd_args=commandArgs(TRUE)

dir <- cmd_args[1]
tsv.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
ecav.filename <- cmd_args[4]
out.filename <- cmd_args[5]

set.seed(1)

source("common.R") # common function for mrlocus and ecaviar-mrlocus
out <- getTrimmedSumStats(dir, tsv.filename,ld.filename)
sum_stat <- out$sum_stat
ld_mat <- out$ld_mat
trim.clumps <- out$trim.clumps

# ecaviar
load(ecav.filename)
if (length(trim.clumps) > 0) {
  ecav.coloc <- ecav.coloc[-trim.clumps]
}

nclumps <- length(sum_stat)
beta_hat_a <- numeric(nclumps)
beta_hat_b <- numeric(nclumps)
sd_a <- numeric(nclumps)
sd_b <- numeric(nclumps)
alleles <- data.frame(id=character(nclumps), ref=character(nclumps), eff=character(nclumps))
for (j in seq_len(nclumps)) {
  idx <- ecav.coloc[[j]]$SNP_ID[ which.max(ecav.coloc[[j]]$CLPP) ]
  row <- sum_stat[[j]][ sum_stat[[j]]$SNP == idx,]
  beta_hat_a[j] <- abs(row$beta_eQTL)
  ref <- ifelse(row$beta_eQTL > 0, row$Ref_eQTL, row$Effect_eQTL)
  eff <- ifelse(row$beta_eQTL > 0, row$Effect_eQTL, row$Ref_eQTL)
  # flip sign so that eQTL effect size is > 0
  beta_hat_b[j] <- sign(row$beta_eQTL) * row$beta_GWAS
  # flip sign to deal with the reference allele being different
  if (row$Ref_eQTL != row$Ref_GWAS) {
    beta_hat_b[j] <- -1 * beta_hat_b[j]
  }
  sd_a[j] <- row$se_eQTL
  sd_b[j] <- row$se_GWAS
  alleles$id[j] <- row$SNP
  alleles$ref[j] <- ref
  alleles$eff[j] <- eff
}
res <- list(beta_hat_a=beta_hat_a, beta_hat_b=beta_hat_b, sd_a=sd_a, sd_b=sd_b)
write.table(alleles, file.path(dir,"ecav_alleles.txt"), row.names=FALSE, quote=FALSE)

# mrlocus
devtools::load_all("../../mrlocus")
res <- fitSlope(res, iter=10000)

save(res, file=out.filename)
