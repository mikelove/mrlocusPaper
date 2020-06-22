system("/proj/milovelab/love/bin/dap/dap_src/dap-g -d_est eqtfile -d_ld test.ld -d_n 1000 -d_syy 1000 -o foo --output_all")
system("grep '((' foo > foo2")
x <- read.table("foo2")
xx <- x[x$V5 > 0,]
y <- read.delim("out/2_1pctmodel_0.1h2g_0.01ve.scan.tsv")
yy <- y[y$eqtl.true != 0,]
table(xx$V2 %in% yy$snp)
y$true <- y$eqtl.true != 0
y$dap <- factor(x$V5[match(y$snp, x$V2)])
plot(y$eqtl.beta/y$eqtl.se, col=y$dap, pch=19)
points(y$eqtl.beta/y$eqtl.se, col=ifelse(y$true, "dodgerblue", NA), cex=2, lwd=2)
