estimate <- function(filename) {
  x <- read.delim(filename)
  with(x, plot(gwas.true ~ eqtl.true))
  c(true=unname(coef(lm(x$gwas.true ~ x$eqtl.true)))[2],
    caus=with(x[x$eqtl.true != 0,],unname(coef(lm(gwas.beta ~ eqtl.beta)))[2]),
    caus.wt=with(x[x$eqtl.true != 0,],unname(coef(lm(gwas.beta ~ eqtl.beta, weight=1/eqtl.se^2)))[2]),
    all=unname(coef(lm(x$gwas.beta ~ x$eqtl.beta)))[2],
    all.wt=unname(coef(lm(x$gwas.beta ~ x$eqtl.beta, weight=1/x$eqtl.se^2)))[2])
}
