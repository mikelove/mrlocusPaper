a <- sub(".clumped","",list.files("out/", pattern="clumped"))
par(mfrow=c(3,1),mar=c(3,3,1,1))
for (i in 1:3) {
  x <- read.table(paste0("out/",a[i],".clumped"), header=TRUE, strings=FALSE)
  y <- read.delim(paste0("out/",a[i],".scan.tsv"))
  yy <- y[y$eqtl.true != 0,]
  table(x$SNP %in% yy$snp)
  noparen <- function(z) sub("\\(1\\)","",z)
  clumps <- lapply(strsplit(x$SP2,split=","), noparen)
  y$true <- y$eqtl.true != 0
  y$index <-y$snp %in% x$SNP
  y$clump <- 1
  for (j in seq_along(clumps)) {
    y$clump[y$snp %in% clumps[[j]] | y$snp == x$SNP[j]] <- j + 1
  }
  plot(y$eqtl.beta/y$eqtl.se, col=y$clump, pch=19)
  points(y$eqtl.beta/y$eqtl.se, col=ifelse(y$true, "dodgerblue", NA), cex=2, lwd=2)
  points(y$eqtl.beta/y$eqtl.se, col=ifelse(y$index, "red", NA), pch=5, cex=3, lwd=2)
}
