cmd_args=commandArgs(TRUE)

scan.filename <- cmd_args[1]
twmr.filename <- cmd_args[2]
ptwas.filename <- cmd_args[3]
mrlocus.filename <- cmd_args[4]
out.filename <- cmd_args[5]

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

out <- estimate(scan.filename)
out <- t(out)

dimn <- list(c("twmr","ptwas","mrlocus"), c("Estimate","Std. Error"))
out <- rbind(out,
             matrix(c(twmr, ptwas, mrlocus),
                    byrow=TRUE, ncol=2, dimnames=dimn))
write.table(format(out, digits=2), file=out.filename, quote=FALSE, sep="\t")
