cmd_args=commandArgs(TRUE)

clumped.filename <- cmd_args[1]
scan.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
out.filename <- cmd_args[4]

if (FALSE) {
  dir <- file.path("out")
  files <- list.files(dir, pattern="clumped")
  file <- sub(".clumped","",files[1])
  dir <- "out"
  clumped.filename <- paste0(dir,"/",file,".clumped")
  scan.filename <- paste0(dir,"/",file,".scan.tsv")
  ld.filename <- paste0(dir,"/",file,".ld")
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
for (j in seq_along(nsnp)) {
  print(j)
  if (nsnp[j] > 1) {
    fit <- fitBetaColoc(beta_hat_a=out2$beta_hat_a[[j]],
                        beta_hat_b=out2$beta_hat_b[[j]],
                        se_a=out2$se_a[[j]],
                        se_b=out2$se_b[[j]],
                        Sigma_a=Sigma[[j]],
                        Sigma_b=Sigma[[j]],
                        verbose=FALSE,
                        open_progress=FALSE,
                        show_messages=FALSE,
                        refresh=-1)
    ## rstan::stan_plot(fit$stan,
    ##                  pars=c(paste0("beta_a[",1:nsnp[j],"]"),
    ##                         paste0("beta_b[",1:nsnp[j],"]")))
    beta_hat_a[[j]] <- fit$beta_hat_a
    beta_hat_b[[j]] <- fit$beta_hat_b
  } else {
    beta_hat_a[[j]] <- out2$beta_hat_a[[j]]
    beta_hat_b[[j]] <- out2$beta_hat_b[[j]]
  }
}

# make a results list for slope fitting
res <- list(beta_hat_a=beta_hat_a,
            beta_hat_b=beta_hat_b,
            sd_a=out2$se_a,
            sd_b=out2$se_b)
res <- extractForSlope(res, plot=FALSE)
res <- fitSlope(res, iter=10000)

#rstan::stan_plot(res$stanfit, pars=c("alpha","sigma"))
if ("stanfit" %in% res) {
  print(res$stanfit, pars=c("alpha","sigma"), digits=3)
  summary(res$stanfit, pars="alpha", probs=c(.1,.9))$summary
}

if (FALSE) {
  coefs <- rstan::extract(res$stanfit)
  x <- res$beta_hat_a
  y <- res$beta_hat_b
  sd_x <- res$sd_a
  sd_y <- res$sd_b
  plot(x,y,pch=19,
       xlim=c(0,1.2*max(x)),
       ylim=c(-1.2*max(abs(y)),1.2*max(abs(y))))
  arrows(x-sd_x,y,x+sd_x,y,angle=90,length=.05,code=3)
  arrows(x,y-sd_y,x,y+sd_y,angle=90,length=.05,code=3)
  abline(v=0, h=0, lty=2)
  abline(0, mean(coefs$alpha), lwd=2)
  abline(mean(coefs$sigma), mean(coefs$alpha), col="blue")
  abline(-mean(coefs$sigma), mean(coefs$alpha), col="blue")
}

if ("stanfit" %in% names(res)) {
  s <- summary(res$stanfit, pars="alpha", probs=c(.025,.05,.1,.9,.95,.975))$summary
  mrlocus.out <- matrix(s[,c("mean","sd","10%","90%","5%","95%","2.5%","97.5%")],
                        ncol=2, byrow=TRUE)
} else {
  mrlocus.out <- matrix(c(res$est,
                          res$est[1] + qnorm(.1) * res$est[2],
                          res$est[1] + qnorm(.9) * res$est[2],
                          res$est[1] + qnorm(.05) * res$est[2],
                          res$est[1] + qnorm(.95) * res$est[2],
                          res$est[1] + qnorm(.025) * res$est[2],
                          res$est[1] + qnorm(.975) * res$est[2]),
                          ncol=2, byrow=TRUE)
}

write.table(mrlocus.out, file=out.filename, row.names=FALSE, col.names=FALSE)
