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
