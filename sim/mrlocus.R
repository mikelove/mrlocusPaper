cmd_args=commandArgs(TRUE)

clumped.filename <- cmd_args[1]
scan.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
out.filename <- cmd_args[4]

if (FALSE) {
  files <- list.files("out", pattern="clumped")
  file <- sub(".clumped","",files[1])
  clumped.filename <- paste0("out/",file,".clumped")
  scan.filename <- paste0("out/",file,".scan.tsv")
  ld.filename <- paste0("out/",file,".ld")
}

clumped <- read.table(clumped.filename, strings=FALSE, header=TRUE)
big_sum_stat <- read.delim(scan.filename, strings=FALSE)
big_ld_mat <- as.matrix(read.table(ld.filename))

noparen <- function(z) sub("\\(1\\)","",z)
clumps <- lapply(strsplit(clumped$SP2,split=","), noparen)
# add index SNP to clumps
clumps <- lapply(seq_along(clumps), function(j) c(clumped$SNP[j], clumps[[j]]))

ld_mat <- lapply(clumps, function(x) {
  idx <- big_sum_stat$snp %in% x
  unname(big_ld_mat[idx,idx,drop=FALSE])
})
sapply(ld_mat, nrow)
sum_stat <- lapply(clumps, function(x) {
  z <- big_sum_stat[big_sum_stat$snp %in% x,]
  z$abs.z <- abs(z$eqtl.beta/z$eqtl.se)
  z
})

devtools::load_all("../../mrlocus")

sapply(sum_stat, function(x) any(x$eqtl.true != 0))
out1 <- collapseHighCorSNPs(sum_stat, ld_mat, score="abs.z", plot=FALSE)
sapply(out1$sum_stat, function(x) any(x$eqtl.true != 0))
eqtl.true <- sapply(out1$sum_stat, function(x) if (any(x$eqtl.true != 0))
                                                 x$eqtl.true[x$eqtl.true != 0] else NA)

out2 <- flipAllelesAndGather(out1$sum_stat, out1$ld_mat,
                             a="eqtl", b="gwas",
                             ref="a0", eff="a1",
                             beta="beta", se="se",
                             a2_plink="a0",
                             snp_id="snp", sep=".",
                             ab_last=FALSE, alleles_same=TRUE,
                             plot=FALSE)

#plotInitEstimates(out2)

library(Matrix)
Sigma <- lapply(out2$Sigma, function(x) as.matrix(nearPD(x)$mat))

library(matrixStats)
#options(mc.cores=3)
nsnp <- lengths(out2$beta_hat_a)
beta_hat_a <- list()
beta_hat_b <- list()
se_a <- list()
se_b <- list()

for (j in seq_along(nsnp)) {
  print(j)
  if (nsnp[j] > 1) {
    fit1 <- fitBetaColoc(beta_hat_a=out2$beta_hat_a[[j]],
                         beta_hat_b=out2$beta_hat_b[[j]],
                         se_a=out2$se_a[[j]],
                         se_b=out2$se_b[[j]],
                         Sigma_a=Sigma[[j]],
                         Sigma_b=Sigma[[j]],
                         verbose=FALSE,
                         open_progress=FALSE,
                         show_messages=FALSE,
                         refresh=-1)
    #rstan::stan_plot(fit1$stan, pars=paste0("beta_a[",1:nsnp[j],"]"))
    #rstan::stan_plot(fit1$stan, pars=paste0("beta_b[",1:nsnp[j],"]"))
    coefs1 <- rstan::extract(fit1$stan)
    beta_hat_a[[j]] <- colMeans(coefs1$beta_a)
    beta_hat_b[[j]] <- colMeans(coefs1$beta_b) / fit1$scale_b
    se_a[[j]] <- colSds(coefs1$beta_a)
    se_b[[j]] <- colSds(coefs1$beta_b) / fit1$scale_b
  } else {
    beta_hat_a[[j]] <- out2$beta_hat_a[[j]]
    beta_hat_b[[j]] <- out2$beta_hat_b[[j]]
    se_a[[j]] <- out2$se_a[[j]]
    se_b[[j]] <- out2$se_b[[j]]
  }
}

#nsnp <- lengths(out2$beta_hat_a)
#plot(unlist(beta_hat_a), unlist(beta_hat_b),
#     col=rep(seq_along(nsnp),each=nsnp),
#     pch=rep(seq_along(nsnp),each=nsnp))

out3 <- extractForSlope(beta_hat_a,
                        beta_hat_b,
                        sd_a=out2$se_a,
                        sd_b=out2$se_b,
                        plot=FALSE)

fit2 <- fitSlope(out3$beta_hat_a,
                 out3$beta_hat_b,
                 out3$sd_a,
                 out3$sd_b,
                 iter=10000)

#rstan::stan_plot(fit2, pars=c("alpha","sigma"))
#print(fit2, pars=c("alpha","sigma"), digits=3)

if (FALSE) {
  coefs2 <- rstan::extract(fit2)
  x <- out3$beta_hat_a
  y <- out3$beta_hat_b
  sd_x <- out3$sd_a
  sd_y <- out3$sd_b
  plot(x,y,xlim=c(0,1.2*max(x)),
       ylim=c(-1.2*max(abs(y)),1.2*max(abs(y))))
  arrows(x-sd_x,y,x+sd_x,y,angle=90,length=.05,code=3)
  arrows(x,y-sd_y,x,y+sd_y,angle=90,length=.05,code=3)
  abline(0, mean(coefs2$alpha), lwd=2)
  abline(mean(coefs2$sigma), mean(coefs2$alpha), col="blue")
  abline(-mean(coefs2$sigma), mean(coefs2$alpha), col="blue")
}

if (is(fit2, "stanfit")) {
  s <- summary(fit2, pars="alpha")$summary
  s[,c("n_eff","Rhat")]
  s[,c("mean","sd")]
  mrlocus.out <- s[,c("mean","sd")]
} else {
  mrlocus.out <- fit2
}

write(mrlocus.out, file=out.filename)
