cmd_args=commandArgs(TRUE)

ss.filename <- cmd_args[1] # is scan.tsv
ld.filename <- cmd_args[2] # is LD matrix

sum.stat <- read.delim(ss.filename, strings=FALSE)
ld <- as.matrix(read.table(ld.filename))

sigma <- as.matrix(Matrix::nearPD(ld)$mat)

library(PMR)
Z1 <- with(sum.stat, eqtl.beta/eqtl.se)
Z2 <- with(sum.stat, gwas.beta/gwas.se)

res <- tryCatch({
  PMR_summary_Egger(Zscore_1=Z1, Zscore_2=Z2,
                    Sigma1sin=sigma, Sigma2sin=sigma,
                    samplen1=500, samplen2=100000,
                    lambda=.15)
  }, error=function(e) NULL)

save(res, file=cmd_args[3])
