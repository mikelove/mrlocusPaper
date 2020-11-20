cmd_args=commandArgs(TRUE)

extra_methods <- FALSE

if (extra_methods) {
  scan.filename <- cmd_args[1]
  twmr.filename <- cmd_args[2]
  twmr2.filename <- cmd_args[3]
  ptwas.filename <- cmd_args[4]
  ptwas2.filename <- cmd_args[5]
  mrlocus.filename <- cmd_args[6]
  mrlocus2.filename <- cmd_args[7]
  ecav.mrlocus.filename <- cmd_args[8]
  out.filename <- cmd_args[9]
} else {
  scan.filename <- cmd_args[1]
  twmr.filename <- cmd_args[2]
  ptwas.filename <- cmd_args[3]
  mrlocus.filename <- cmd_args[4]
  ecav.mrlocus.filename <- cmd_args[5]
  out.filename <- cmd_args[6]
}

ivw.fixed.effects <- function(eqtl.beta, gwas.beta, eqtl.se, gwas.se) {
  res <- TwoSampleMR::mr_ivw_fe(eqtl.beta, gwas.beta, eqtl.se, gwas.se)
  c("Estimate"=res$b, "Std. Error"=res$se)
}

estimate <- function(filename) {
  x <- read.delim(filename)
  mycoef <- function(fit) suppressWarnings(summary(fit)$coefficients[1,1:2])
  out <- cbind(true=mycoef(lm(x$gwas.true ~ x$eqtl.true + 0)),
               caus=with(x[x$eqtl.true != 0,],mycoef(lm(gwas.beta ~ eqtl.beta + 0))),
               caus.wt=with(x[x$eqtl.true != 0,],ivw.fixed.effects(eqtl.beta, gwas.beta, eqtl.se, gwas.se)),
               all=mycoef(lm(x$gwas.beta ~ x$eqtl.beta + 0)),
               all.wt=with(x, ivw.fixed.effects(eqtl.beta, gwas.beta, eqtl.se, gwas.se)))
  out[2,1] <- NA
  out
}

twmr <- as.numeric(unname(read.table(twmr.filename, header=TRUE)[1,2:3]))

ptwas <- scan(ptwas.filename, what="char", sep="\n")
ptwas <- ptwas[length(ptwas)]
ptwas <- as.numeric(trimws(strsplit(ptwas, "\t")[[1]][5:6]))

mrlocus <- unname(as.matrix(read.table(mrlocus.filename, header=FALSE))[1,])

ecav.mrlocus <- unname(as.matrix(read.table(ecav.mrlocus.filename, header=FALSE))[1,])

if (extra_methods) {
  # twmr at an alternative threshold
  if (file.info(twmr2.filename)$size > 0) {
    twmr2 <- as.numeric(unname(read.table(twmr2.filename, header=TRUE)[1,2:3]))
  } else {
    twmr2 <- c(0,0)
  }
  
  # ptwas at an alternative threshold
  ptwas2 <- scan(ptwas2.filename, what="char", sep="\n")
  ptwas2 <- ptwas2[length(ptwas2)]
  ptwas2 <- as.numeric(trimws(strsplit(ptwas2, "\t")[[1]][5:6]))
  
  # mrlocus at an alternative threshold
  if (file.info(mrlocus2.filename)$size > 0) {
    mrlocus2 <- unname(as.matrix(read.table(mrlocus2.filename, header=FALSE))[1,])
  } else {
    mrlocus2 <- c(0,0)
  }
}

out <- estimate(scan.filename)
out <- t(out)

if (extra_methods) {
  dimn <- list(c("twmr","twmr2","ptwas","ptwas2","mrlocus","mrlocus2","ecav-mrlocus"),
               c("Estimate","Std. Error"))
  out <- rbind(out,
               matrix(c(twmr, twmr2, ptwas, ptwas2, mrlocus, mrlocus2, ecav.mrlocus),
                      byrow=TRUE, ncol=2, dimnames=dimn))
} else {
  dimn <- list(c("twmr","ptwas","mrlocus","ecav-mrlocus"),
               c("Estimate","Std. Error"))
  out <- rbind(out,
               matrix(c(twmr, ptwas, mrlocus, ecav.mrlocus),
                      byrow=TRUE, ncol=2, dimnames=dimn))
}

write.table(format(out, digits=2), file=out.filename, quote=FALSE, sep="\t")
