cmd_args=commandArgs(TRUE)

dir <- cmd_args[1]
tsv.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
ecav.filename <- cmd_args[4]
out.filename <- cmd_args[5]

dir <- "Artery_MRAS_CAD"
tsv.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.tsv"
ld.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.ld"
ecav.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.ecav"

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
for (j in seq_len(nclumps)) {
  idx <- ecav.coloc[[j]]$SNP_ID[ which.max(ecav.coloc[[j]]$CLPP) ]
  row <- sum_stat[[j]][ sum_stat[[j]]$SNP == idx,]
  beta_hat_a[j] <- abs(row$beta_eQTL)
  beta_hat_b[j] <- sign(row$beta_eQTL) * row$beta_GWAS
  sd_a[j] <- row$se_eQTL
  sd_b[j] <- row$se_GWAS
}
res <- list(beta_hat_a=beta_hat_a, beta_hat_b=beta_hat_b, sd_a=sd_a, sd_b=sd_b)

# mrlocus
devtools::load_all("../../mrlocus")
res <- fitSlope(res, iter=10000)

save(res, file=out.filename)
