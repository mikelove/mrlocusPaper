files <- list.files(pattern="*.txt")
dat <- data.frame(t(unname(sapply(files, scan, quiet=TRUE))))
names(dat) <- c("n","r2","FPR")
dat$n <- factor(dat$n, c(8,6,4))
dat$r2 <- factor(dat$r2)
library(dplyr)
dat <- dat %>% mutate(se=sqrt(FPR*(1-FPR)/200),ymin=FPR-1.96*se, ymax=FPR+1.96*se)
library(ggplot2)
ggplot(dat, aes(r2, FPR, col=n, group=n, ymin=ymin, ymax=ymax)) +
  geom_point() + 
  geom_errorbar(width=.2) +
  geom_line() +
  geom_hline(yintercept=0.2, lty=2)
  
