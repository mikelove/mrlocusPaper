estimate <- function(filename) {
  x <- read.delim(filename)
  with(x, plot(gwas.true ~ eqtl.true))
  c(true=unname(coef(lm(x$gwas.true ~ x$eqtl.true + 0)))[1],
    caus=with(x[x$eqtl.true != 0,],unname(coef(lm(gwas.beta ~ eqtl.beta + 0)))[1]),
    caus.wt=with(x[x$eqtl.true != 0,],unname(coef(lm(gwas.beta ~ eqtl.beta + 0, weight=1/eqtl.se^2)))[1]),
    all=unname(coef(lm(x$gwas.beta ~ x$eqtl.beta + 0)))[1],
    all.wt=unname(coef(lm(x$gwas.beta ~ x$eqtl.beta + 0, weight=1/x$eqtl.se^2)))[1])
}
