sim <- read.delim("out/1_1pctmodel_0.1h2g_0.001ve.scan.tsv", strings=FALSE)
idx <- which(sim$eqtl.true != 0)
with(sim, plot(eqtl.beta/eqtl.se, xlim=c(0,50),
               col=ifelse(eqtl.true == 0, 1, 2),
               pch=ifelse(eqtl.true == 0, 1, 19)))
snps <- sim$snp[1:50]
write(snps, file="snp.list")
# /proj/milovelab/bin/plink/plink --bfile out/1_1pctmodel_0.1h2g_0.001ve --r2 --ld-snp-list snp.list --allow-no-sex --ld-window-r2 0 --ld-window 1000000 --out snp --memory 2048

ld <- read.table("snp.ld", header=TRUE)
ld <- ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]

ld.mat <- matrix(nrow=50, ncol=50, dimnames=list(snps, snps))
