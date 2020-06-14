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

nsnp <- sapply(sum_stat, nrow)
out <- flipAllelesAndGather(sum_stat, ld_mat,
                            a="eQTL", b="GWAS",
                            ref="Ref", eff="Effect",
                            beta="beta", se="se",
                            major_plink="Major_plink",
                            sep="_", ab_last=TRUE)
plot(unlist(out$beta_hat_a), unlist(out$beta_hat_b), col=rep(1:3,nsnp))

library(Matrix)
Sigma_np <- out$Sigma
for (j in 1:3) {
  Sigma_npd[[j]] <- as.matrix(nearPD(Sigma_a[[j]])$mat)
}
Sigma_b <- Sigma_a <- Sigma_npd

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
