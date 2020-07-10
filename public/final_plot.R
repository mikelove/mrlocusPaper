dat <- read.delim("estimates.tsv")
library(ggplot2)
q <- qnorm(.9)
ggplot(dat, aes(x=alpha,y=method,xmin=alpha-q*se,xmax=alpha+q*se)) +
  geom_pointrange() +
  facet_wrap(~tissue * gene, scales="free") +
  geom_vline(xintercept = 0, lty=2)
