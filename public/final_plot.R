dat <- read.delim("estimates.tsv")
library(ggplot2)
library(dplyr)
q <- qnorm(.9)
dat <- dat %>% mutate(xmin=alpha-q*se,xmax=alpha+q*se)

# replace SD based intervals with quantile based
ints <- c(0.04757, 0.13474,
          -0.3938, -0.2832,
          0.02920, 0.06893,
          -0.9084,  0.0652,
          -0.15772,  0.04035,
          -0.04712, -0.04113)
ints <- matrix(ints, ncol=2, byrow=TRUE)
dat[13:18,6:7] <- ints

png(file="~/Desktop/realdata.png", res=125, width=1200, height=600)
ggplot(dat, aes(x=alpha,y=method,xmin=xmin,xmax=xmax)) +
  geom_pointrange() +
  facet_wrap(~tissue * gene, scales="free") +
  geom_vline(xintercept = 0, lty=2) +
  xlab("gene-to-trait estimate")
dev.off()
