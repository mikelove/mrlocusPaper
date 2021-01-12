pre <- as.numeric(scan("collapse_pre.txt", what="char"))
post <- as.numeric(scan("collapse_post.txt", what="char"))

dat <- data.frame(type=factor(rep(c("pre","post"),
                                  c(length(pre),length(post))),
                              c("pre","post")),
                  number=c(pre,post))

library(ggplot2)
pdf(file="../supp/figs/snps_per_clump.pdf", height=4)
ggplot(dat, aes(number, after_stat(density), color=type)) +
  geom_histogram(aes(fill=type), position="identity", alpha=.2, color="grey50", breaks=0:30*5) +
  geom_density(size=1) + 
  xlab("number of SNPs per clump") +
  xlim(0,150)
dev.off()

library(dplyr)
dat %>%
  group_by(type) %>%
  summarize(mean=mean(number),median=median(number),n=n())
# mean 33.8 -> 15.9
# median 26 -> 13
# n 1517 -> 1250
