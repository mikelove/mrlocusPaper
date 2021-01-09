sum_stat <- read.delim("out/hp/hp_2r_500nq_1pctm_0.1h2_0.01ve.scan.tsv")
library(dplyr)
library(ggplot2)
dat <- sum_stat %>%
  mutate(gwas.z=abs(gwas.beta/gwas.se),
         eqtl.z=abs(eqtl.beta/eqtl.se),
         status=ifelse(gwas.true != 0, "gene exp. mediated",
                ifelse(gwas.true.hp != 0, "trait-only", "non-causal")))
table(dat$status)
dat2 <- dat %>% filter(status != "non-causal")
pdf(file="../supp/figs/hp_example.pdf", height=3, width=7)
ggplot(dat, aes(pos, gwas.z, col=status)) +
  geom_point() + xlab("position") + ylab("GWAS |z-score|") + 
  geom_hline(yintercept=0) +
  geom_segment(data=dat2, aes(xend=pos, yend=0), lineend="butt", show.legend=FALSE)
dev.off()
