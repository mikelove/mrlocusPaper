cmd_args=commandArgs(TRUE)

dir <- cmd_args[1]
tsv.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
ecav.filename <- cmd_args[4]
out.filename <- cmd_args[5]

if (FALSE) {
  dir <- "Artery_MRAS_CAD"
  tsv.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.tsv"
  ld.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.ld"
  ecav.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.ecav"
}

set.seed(1)

devtools::load_all("../../mrlocus")
source("common.R") # common function for mrlocus and ecaviar-mrlocus
out <- getTrimmedSumStats(dir, tsv.filename, ld.filename, r2_threshold=0.05)
sum_stat <- out$sum_stat
ld_mat <- out$ld_mat
trim_clusters <- out$trim_clusters

# ecaviar
load(ecav.filename)
if (length(trim_clusters) > 0) {
  ecav.coloc <- ecav.coloc[-trim_clusters]
}

nclusters <- length(sum_stat)
beta_hat_a <- numeric(nclusters)
beta_hat_b <- numeric(nclusters)
sd_a <- numeric(nclusters)
sd_b <- numeric(nclusters)
alleles <- data.frame(id=character(nclusters), ref=character(nclusters), eff=character(nclusters))
z.thr <- qnorm(.001/2, lower.tail=FALSE)
for (j in seq_len(nclusters)) {
  m <- merge(sum_stat[[j]], ecav.coloc[[j]], by.x="SNP", by.y="SNP_ID")
  # remove SNPs below p threshold
  m <- m[ m$abs.z > z.thr,]
  if (nrow(m) == 0) next
  row <- m[ which.max(m$CLPP), ]
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

# remove any clusters that were below p threshold
alleles <- alleles[res$beta_hat_a != 0,]
res <- lapply(res, function(x) x[x != 0])

# second round trimming based on candidate SNPs
load("LDmatrix_allSNPs.rda")
r2 <-  get(paste0("LD_", dir))
idx <- sapply(alleles$id, function(snp) grep(snp, colnames(r2)))
r2 <- r2[idx,idx]
trim_clusters <- trimClusters(r2, r2_threshold=0.05)
trim2.filename <- file.path(dir, paste0(dir, "_ecav.trim2"))
if (length(trim_clusters) > 0) {
  alleles <- alleles[-trim_clusters,,drop=FALSE]
  res <- lapply(res, function(x) x[-trim_clusters])
  write(trim_clusters, file=trim2.filename, ncolumns=length(trim_clusters))
} else {
  write(NA, file=trim2.filename)
}
write.table(alleles, file.path(dir,"ecav_alleles.txt"), row.names=FALSE, quote=FALSE)

# mrlocus
res <- fitSlope(res, iter=10000)

save(res, file=out.filename)
