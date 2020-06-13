dir <- file.path("Artery_Tibial","MRAS")
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

nsnp <- sapply(sum.stats, nrow)
Sigma_a <- list()
beta_hat_a <- list()
beta_hat_b <- list()
se_a <- list()
se_b <- list()
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
  beta_hat_a[[j]] <- beta_a
  beta_hat_b[[j]] <- beta_b
  se_a[[j]] <- sum.stats[[j]]$se_eQTL
  se_b[[j]] <- sum.stats[[j]]$se_GWAS
  # flip LD matrix so that alleles are positively correlated with eQTL index
  ld.sign <- sign(ld.mat[[j]][,idx])
  ld.flipped <- t(t(ld.mat[[j]]) * ld.sign) * ld.sign
  Sigma_a[[j]] <- ld.flipped
}

plot(unlist(beta_hat_a), unlist(beta_hat_b), col=rep(1:3,nsnp))

library(Matrix)
Sigma_np <- Sigma_a
for (j in 1:3) {
  Sigma_np[[j]] <- as.matrix(nearPD(Sigma_a[[j]])$mat)
}

pheatmap::pheatmap(Sigma_a[[1]], cluster_rows=FALSE, cluster_cols=FALSE)
pheatmap::pheatmap(Sigma_np[[1]], cluster_rows=FALSE, cluster_cols=FALSE)
Sigma_b <- Sigma_a <- Sigma_np

load_all("../../mrlocus")
options(mc.cores=2)
j <- 3
fit1 <- fitBetaEcaviar(nsnp=nsnp[j],
                       beta_hat_a=beta_hat_a[[j]],
                       beta_hat_b=beta_hat_b[[j]],
                       se_a=se_a[[j]],
                       se_b=se_b[[j]],
                       Sigma_a=Sigma_a[[j]],
                       Sigma_b=Sigma_b[[j]])

rstan::stan_plot(fit1, pars=paste0("beta_a[",1:nsnp[j],"]"))
rstan::stan_plot(fit1, pars=paste0("beta_b[",1:nsnp[j],"]"))
print(fit1, pars=paste0("beta_a[",1:nsnp[3],"]"), digits=3)
coefs1 <- rstan::extract(fit1)
par(mfrow=c(1,2))
beta_clean_a <- colMeans(coefs1$beta_a)
plot(beta_hat_a[[j]], beta_clean_a, asp=1)
abline(0,1); abline(h=0, lty=2)
beta_clean_b <- colMeans(coefs1$beta_b)
plot(beta_hat_b[[j]], beta_clean_b, asp=1)
abline(0,1); abline(h=0, lty=2)
