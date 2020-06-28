cmd_args=commandArgs(TRUE)

scan.filename <- cmd_args[1]
twmr.filename <- cmd_args[2]
ptwas.filename <- cmd_args[3]
mrlocus.filename <- cmd_args[4]
out.filename <- cmd_args[5]

estimate <- function(filename) {
  x <- read.delim(filename)
  #with(x, plot(gwas.true ~ eqtl.true))
  mycoef <- function(fit) suppressWarnings(summary(fit)$coefficients[1,1:2])
  out <- cbind(true=mycoef(lm(x$gwas.true ~ x$eqtl.true + 0)),
               caus=with(x[x$eqtl.true != 0,],mycoef(lm(gwas.beta ~ eqtl.beta + 0))),
               caus.wt=with(x[x$eqtl.true != 0,],mycoef(lm(gwas.beta ~ eqtl.beta + 0, weight=1/eqtl.se^2))),
               all=mycoef(lm(x$gwas.beta ~ x$eqtl.beta + 0)),
               all.wt=mycoef(lm(x$gwas.beta ~ x$eqtl.beta + 0, weight=1/x$eqtl.se^2)))
  out[2,1] <- NA
  out
}

twmr <- as.numeric(unname(read.table(twmr.filename, header=TRUE)[1,2:3]))

ptwas <- scan(ptwas.filename, what="char", sep="\n")
ptwas <- ptwas[length(ptwas)]
ptwas <- as.numeric(trimws(strsplit(ptwas, "\t")[[1]][5:6]))

mrlocus <- scan(mrlocus.filename)

out <- estimate(scan.filename)
out <- t(out)

dimn <- list(c("twmr","ptwas","mrlocus"), c("Estimate","Std. Error"))
out <- rbind(out,
             matrix(c(twmr, ptwas, mrlocus),
                    byrow=TRUE, ncol=2, dimnames=dimn))
write.table(format(out, digits=2), file=out.filename, quote=FALSE, sep="\t")
