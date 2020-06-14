dir <- file.path("Artery_Tibial","PHACTR1")
cond.files <- sub(".tsv","",list.files(dir, ".tsv"))
nclust <- length(cond.files)
ld_mat <- list()
sum_stat <- list()
for (j in 1:nclust) {
  filename <- file.path(dir,paste0(cond.files[j], ".ld"))
  ld_mat[[j]] <- as.matrix(read.delim(filename,header=FALSE))
  filename <- file.path(dir,paste0(cond.files[j], ".tsv"))
  sum_stat[[j]] <- read.delim(filename)
  sum_stat[[j]] <- sum_stat[[j]][order(sum_stat[[j]]$pos),]
}

ld_mat <- ld_mat[2:4]
sum_stat <- sum_stat[2:4]

sapply(sum_stat, nrow)

library(pheatmap)
library(gridExtra)
load_all("../../mrlocus")
out <- collapseHighCorSNPs(sum_stat, ld_mat)
sum_stat <- out$sum_stat
ld_mat <- out$ld_mat

flipAlleles <- function(sum_stat, ld_mat) {
  Sigma <- list()
  beta_hat_a <- list()
  beta_hat_b <- list()
  se_a <- list()
  se_b <- list()
  for (j in seq_along(sum_stat)) {
    # the index SNP in the signal cluster for A (eQTL)
    idx <- with(sum_stat[[j]], which.max(abs(beta_eQTL/se_eQTL)))
    # reference B (GWAS) allele must be either reference or effect allele in A (eQTL)
    with(sum_stat[[j]], stopifnot(Ref_GWAS == Ref_eQTL | Ref_GWAS == Effect_eQTL))
    # flip B (GWAS) so that same allele is described for B (GWAS) as for A (eQTL)
    flip <- with(sum_stat[[j]], which(Ref_GWAS != Ref_eQTL))
    sum_stat[[j]]$beta_GWAS_flipped <- sum_stat[[j]]$beta_GWAS
    sum_stat[[j]]$beta_GWAS_flipped[flip] <- -1 * sum_stat[[j]]$beta_GWAS[flip]
    # flip alleles other than the index so they have positive correlation (LD) with index
    ld.sign <- sign(ld_mat[[j]][,idx])
    beta_a <- with(sum_stat[[j]], ifelse(Major_plink == Ref_eQTL, 1, -1) * ld.sign * beta_eQTL)
    beta_b <- with(sum_stat[[j]], ifelse(Major_plink == Ref_eQTL, 1, -1) * ld.sign * beta_GWAS_flipped)
    ld.flipped <- t(t(ld_mat[[j]]) * ld.sign) * ld.sign
    Sigma[[j]] <- ld.flipped
    # flip alleles so that effect size is positive for index SNP for A (eQTL)
    if (beta_a[idx] < 0) {
      beta_a <- beta_a * -1
      beta_b <- beta_b * -1
    }
    beta_hat_a[[j]] <- beta_a
    beta_hat_b[[j]] <- beta_b
    se_a[[j]] <- sum_stat[[j]]$se_eQTL
    se_b[[j]] <- sum_stat[[j]]$se_GWAS
  }
}

nsnp <- sapply(sum_stat, nrow)
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
