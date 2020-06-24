files <- sub(".clumped","",list.files("out",pattern=".clumped"))
i <- 3
clumped <- read.table(paste0("out/",files[i],".clumped"), strings=FALSE, header=TRUE)
big_sum_stat <- read.delim(paste0("out/",files[i],".scan.tsv"), strings=FALSE)
big_ld_mat <- as.matrix(read.table(paste0("out/",files[i],".ld")))

noparen <- function(z) sub("\\(1\\)","",z)
clumps <- lapply(strsplit(clumped$SP2,split=","), noparen)
# add index SNP to clumps
clumps <- lapply(seq_along(clumps), function(j) c(clumped$SNP[j], clumps[[j]]))

ld_mat <- lapply(clumps, function(x) {
  idx <- big_sum_stat$snp %in% x
  unname(big_ld_mat[idx,idx])
})
sapply(ld_mat, nrow)
sum_stat <- lapply(clumps, function(x) {
  z <- big_sum_stat[big_sum_stat$snp %in% x,]
  z$abs.z <- abs(z$eqtl.beta/z$eqtl.se)
  z
})

devtools::load_all("../../mrlocus")

library(pheatmap)
sapply(sum_stat, function(x) any(x$eqtl.true != 0))
eqtl.true <- sapply(sum_stat, function(x) if (any(x$eqtl.true != 0))
                                            x$eqtl.true[x$eqtl.true != 0] else NA)
out1 <- collapseHighCorSNPs(sum_stat, ld_mat, score="abs.z", plot=FALSE)
sapply(out1$sum_stat, function(x) any(x$eqtl.true != 0))

out2 <- flipAllelesAndGather(out1$sum_stat, out1$ld_mat,
                             a="eqtl", b="gwas",
                             ref="a0", eff="a1",
                             beta="beta", se="se",
                             a2_plink="a0",
                             snp_id="snp", sep=".",
                             ab_last=FALSE, alleles_same=TRUE)

plotInitEstimates(out2)

library(Matrix)
Sigma <- lapply(out2$Sigma, function(x) as.matrix(nearPD(x)$mat))

library(matrixStats)
options(mc.cores=1)
nsnp <- lengths(out2$beta_hat_a)
beta_hat_a <- list()
beta_hat_b <- list()
se_a <- list()
se_b <- list()
#save(out1, out2, Sigma, file="mrlocus_input.rda")

for (j in seq_along(nsnp)) {
  if (nsnp[j] > 1) {
    system.time({
      fit1 <- fitBetaEcaviar(nsnp=nsnp[j],
                             beta_hat_a=out2$beta_hat_a[[j]],
                             beta_hat_b=out2$beta_hat_b[[j]],
                             se_a=out2$se_a[[j]],
                             se_b=out2$se_b[[j]],
                             Sigma_a=Sigma[[j]],
                             Sigma_b=Sigma[[j]],
                             verbose=FALSE,
                             open_progress=FALSE,
                             show_messages=FALSE,
                             refresh=-1, iter=10000)
    })
    rstan::stan_plot(fit1, pars=paste0("beta_a[",1:nsnp[j],"]"))
    #rstan::stan_plot(fit1, pars=paste0("beta_b[",1:nsnp[j],"]"))
    coefs1 <- rstan::extract(fit1)
    beta_hat_a[[j]] <- colMeans(coefs1$beta_a)
    beta_hat_b[[j]] <- colMeans(coefs1$beta_b)
    se_a[[j]] <- colSds(coefs1$beta_a)
    se_b[[j]] <- colSds(coefs1$beta_b)
  } else {
    beta_hat_a[[j]] <- out2$beta_hat_a[[j]]
    beta_hat_b[[j]] <- out2$beta_hat_b[[j]]
    se_a[[j]] <- out2$se_a[[j]]
    se_b[[j]] <- out2$se_b[[j]]
  }
}

save(beta_hat_a, beta_hat_b, se_a, se_b, file="mrlocus_part1.rda")

plotInitEstimates(out2)

nsnp <- lengths(beta_hat_a)
plot(unlist(beta_hat_a), unlist(beta_hat_b), col=rep(1:5,nsnp))
beta_hat_a <- unlist(beta_hat_a)
beta_hat_b <- unlist(beta_hat_b)
sd_a <- unlist(se_a)
sd_b <- unlist(se_b)
idx <- beta_hat_a > .05
n <- sum(idx)

fit2 <- fitBetaNonzero(beta_hat_a[idx],
                       beta_hat_b[idx],
                       sd_a[idx], sd_b[idx],
                       iter=4000)

rstan::stan_plot(fit2, pars=c("alpha","sigma"))
print(fit2, pars=c("alpha","sigma"), digits=3)

coefs2 <- rstan::extract(fit2)

plot(beta_hat_a[idx], beta_hat_b[idx],
     xlim=c(0,max(beta_hat_a[idx])),
     ylim=c(0,max(beta_hat_b[idx])))
abline(0, mean(coefs2$alpha), lwd=2)
abline(mean(coefs2$sigma), mean(coefs2$alpha), col="blue")
abline(-mean(coefs2$sigma), mean(coefs2$alpha), col="blue")

