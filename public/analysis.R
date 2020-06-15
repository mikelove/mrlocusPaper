dir <- file.path("Artery_Tibial","PHACTR1")
dir <- file.path("Liver","CETP")
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

sapply(sum_stat, nrow)

library(pheatmap)
library(gridExtra)
load_all("../../mrlocus")
out1 <- collapseHighCorSNPs(sum_stat, ld_mat)
sum_stat <- out1$sum_stat
ld_mat <- out1$ld_mat

nsnp <- sapply(sum_stat, nrow)
out2 <- flipAllelesAndGather(sum_stat, ld_mat,
                            a="eQTL", b="GWAS",
                            ref="Ref", eff="Effect",
                            beta="beta", se="se",
                            major_plink="Major_plink",
                            sep="_", ab_last=TRUE)

library(Matrix)
Sigma_npd <- out2$Sigma
for (j in seq_along(out2$Sigma)) {
  Sigma_npd[[j]] <- as.matrix(nearPD(out2$Sigma[[j]])$mat)
}
Sigma_b <- Sigma_a <- Sigma_npd


plotInitEstimates(out2)

load_all("../../mrlocus")
options(mc.cores=2)
j <- 2
fit1 <- fitBetaEcaviar(nsnp=nsnp[j],
                       beta_hat_a=out2$beta_hat_a[[j]],
                       beta_hat_b=out2$beta_hat_b[[j]],
                       se_a=out2$se_a[[j]],
                       se_b=out2$se_b[[j]],
                       Sigma_a=Sigma_a[[j]],
                       Sigma_b=Sigma_b[[j]],
                       iter=2000)

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
