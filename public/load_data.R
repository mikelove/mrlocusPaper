dir <- file.path("Artery_Tibial","PHACTR1")
cond.files <- sub(".tsv","",list.files(dir, ".tsv"))
ncond <- length(cond.files)
ld.mat <- list()
sum.stats <- list()
for (j in 1:ncond) {
  filename <- file.path(dir,paste0(cond.files[j], ".ld"))
  ld.mat[[j]] <- as.matrix(read.delim(filename,header=FALSE))
  filename <- file.path(dir,paste0(cond.files[j], ".tsv"))
  sum.stats[[j]] <- read.delim(filename)
  sum.stats[[j]] <- sum.stats[[j]][order(sum.stats[[j]]$pos),]
}

ld.mat <- ld.mat[2:4]
sum.stats <- sum.stats[2:4]

sapply(sum.stats, nrow)
# 74 113 157
# 74 75 12

library(pheatmap)
for (j in 1:3) {
  pheatmap(ld.mat[[j]])
  hc <- hclust(as.dist(2-ld.mat[[j]]))
  thresh <- .05
  plot(hc); abline(h=1+thresh,col="red")
  clusters <- cutree(hc, h=1+thresh)
  nclust <- max(clusters)
  idx <- match(seq_len(nclust), clusters)
  ld.mat[[j]] <- ld.mat[[j]][idx,idx]
  sum.stats[[j]] <- sum.stats[[j]][idx,]
}

pheatmap(ld.mat[[3]])

sapply(sum.stats, nrow)

nsnp <- 50
ncond <- 3
Sigma_a <- array(NA, dim=c(nsnp,nsnp,ncond))
beta_hat_a <- matrix(NA, nrow=nsnp, ncol=ncond)
beta_hat_b <- matrix(NA, nrow=nsnp, ncol=ncond)
se_a <- matrix(NA, nrow=nsnp, ncol=ncond)
se_b <- matrix(NA, nrow=nsnp, ncol=ncond)
indices <- list(25:74, 51:100, 1:50)
for (j in 1:3) {
  idx <- with(sum.stats[[j]], which.max(abs(beta_eQTL/se_eQTL)))
  print(idx)
  with(sum.stats[[j]], stopifnot(Ref_GWAS == Ref_eQTL | Ref_GWAS == Effect_eQTL))
  with(sum.stats[[j]], table(Ref_GWAS == Ref_eQTL))
  flip <- with(sum.stats[[j]], which(Ref_GWAS != Ref_eQTL))
  # flip GWAS so that same allele is described for GWAS as for eQTL
  sum.stats[[j]]$beta_GWAS_flipped <- sum.stats[[j]]$beta_GWAS
  sum.stats[[j]]$beta_GWAS_flipped[flip] <- -1 * sum.stats[[j]]$beta_GWAS[flip]
  ld.sign <- sign(ld.mat[[j]][,idx])
  beta_a <- with(sum.stats[[j]], ifelse(Major_plink == Ref_eQTL, 1, -1) * ld.sign * beta_eQTL)
  beta_b <- with(sum.stats[[j]], ifelse(Major_plink == Ref_eQTL, 1, -1) * ld.sign * beta_GWAS_flipped)
  # flip signs so that effect size is in terms of expression increasing allele in eQTL
  if (beta_a[idx] < 0) {
    beta_a <- beta_a * -1
    beta_b <- beta_b * -1
  }
  #plot(beta_a, beta_b)
  beta_hat_a[,j] <- beta_a[indices[[j]]]
  beta_hat_b[,j] <- beta_b[indices[[j]]]
  se_a[,j] <- sum.stats[[j]]$se_eQTL[indices[[j]]]
  se_b[,j] <- sum.stats[[j]]$se_GWAS[indices[[j]]]
  # flip LD matrix so that alleles are positively correlated with eQTL index
  ld.sign <- sign(ld.mat[[j]][,idx])
  ld.flipped <- t(t(ld.mat[[j]]) * ld.sign) * ld.sign
  Sigma_a[,,j] <- ld.flipped[indices[[j]],indices[[j]]]
}

plot(beta_hat_a, beta_hat_b, col=rep(1:3,each=50))

library(Matrix)
Sigma_np <- Sigma_a
for (j in 1:3) {
  Sigma_np[,,j] <- as.matrix(nearPD(Sigma_a[,,j])$mat)
}
pheatmap::pheatmap(Sigma_a[,,1], cluster_rows=FALSE, cluster_cols=FALSE)
pheatmap::pheatmap(Sigma_np[,,1], cluster_rows=FALSE, cluster_cols=FALSE)
Sigma_b <- Sigma_a <- Sigma_np

data <- list(nsnp=nsnp,
             ncond=ncond,
             beta_hat_a=beta_hat_a,
             beta_hat_b=beta_hat_b,
             se_a=se_a,
             se_b=se_b,
             Sigma_a=Sigma_a,
             Sigma_b=Sigma_b)
load_all("../../mrlocus")
options(mc.cores=2)
fit1 <- fitBetaEcaviar(data)
