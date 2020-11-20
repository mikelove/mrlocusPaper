cmd_args=commandArgs(TRUE)

clumped.filename <- cmd_args[1]
scan.filename <- cmd_args[2]
keep.filename <- cmd_args[3]
ecav.filename <- cmd_args[4]
out.filename <- cmd_args[5]

if (FALSE) {
  i <- "high_n"
  dir <- file.path("out",i)
  files <- list.files(dir, pattern="ve.clumped")
  file <- sub(".clumped","",files[24])
  clumped.filename <- paste0(dir,"/",file,".clumped")
  scan.filename <- paste0(dir,"/",file,".scan.tsv")
  keep.filename <- paste0(dir,"/",file,".mrl_keep")
  ecav.filename <- paste0(dir,"/",file,".ecav")
}

clumped <- read.table(clumped.filename, strings=FALSE, header=TRUE)
big_sum_stat <- read.delim(scan.filename, strings=FALSE)

noparen <- function(z) sub("\\(1\\)","",z)
clumps <- lapply(strsplit(clumped$SP2,split=","), noparen)
# add index SNP to clumps
clumps <- lapply(seq_along(clumps), function(j) c(clumped$SNP[j], clumps[[j]]))

sum_stat <- lapply(clumps, function(x) {
  out <- big_sum_stat[big_sum_stat$snp %in% x,]
  out$z <- out$eqtl.beta/out$eqtl.se
  out$abs.z <- abs(out$z)
  out$z2 <- out$gwas.beta/out$gwas.se # for testing ecaviar
  out
})

# read in the clumps to keep
keep.clumps <- scan(keep.filename)

# trim clumps
sum_stat <- sum_stat[keep.clumps]

# ecaviar
load(ecav.filename)
ecav.coloc <- ecav.coloc[keep.clumps]

nclumps <- length(sum_stat)
beta_hat_a <- numeric(nclumps)
beta_hat_b <- numeric(nclumps)
sd_a <- numeric(nclumps)
sd_b <- numeric(nclumps)
for (j in seq_len(nclumps)) {
  idx <- ecav.coloc[[j]]$SNP_ID[ which.max(ecav.coloc[[j]]$CLPP) ]
  row <- sum_stat[[j]][ sum_stat[[j]]$snp == idx,]
  beta_hat_a[j] <- abs(row$eqtl.beta)
  beta_hat_b[j] <- sign(row$eqtl.beta) * row$gwas.beta
  sd_a[j] <- row$eqtl.se
  sd_b[j] <- row$gwas.se
}
res <- list(beta_hat_a=beta_hat_a, beta_hat_b=beta_hat_b, sd_a=sd_a, sd_b=sd_b)

# mrlocus
devtools::load_all("../../mrlocus")
res <- fitSlope(res, iter=10000)

#rstan::stan_plot(res$stanfit, pars=c("alpha","sigma"))
if ("stanfit" %in% res) {
  print(res$stanfit, pars=c("alpha","sigma"), digits=3)
  summary(res$stanfit, pars="alpha", probs=c(.1,.9))$summary
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
