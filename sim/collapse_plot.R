pre <- as.numeric(scan("collapse_pre.txt", what="char"))
post <- as.numeric(scan("collapse_post.txt", what="char"))

dat <- data.frame(type=factor(rep(c("pre","post"),each=length(pre)),c("pre","post")),
                  number=c(pre,post))

library(ggplot2)
png(file="~/Desktop/snps-per-clump.png", width=800, height=800, res=125)
ggplot(dat, aes(x=number, fill=type)) +
  geom_histogram(position="identity",
                 alpha=.25, color="grey50",
                 breaks=0:30*5) +
  theme_bw() + xlab("number of SNPs per clump")
dev.off()
