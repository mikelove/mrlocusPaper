files <- list.files(pattern="*.txt")
dat <- data.frame(t(unname(sapply(files, scan, quiet=TRUE))))
names(dat) <- c("n","r2","FPR")
dat$n <- factor(dat$n)
dat$r2 <- factor(dat$r2)
library(dplyr)
niter <- 400
dat <- dat %>% mutate(se=sqrt(FPR*(1-FPR)/niter),ymin=FPR-1.96*se, ymax=FPR+1.96*se)
save(dat, file="corr_instr_sim.rda")
library(ggplot2)
pdf(file="../supp/figs/corr_instr_sim.pdf", width=5, height=5)
ggplot(dat, aes(r2, FPR, col=n, group=n, ymin=ymin, ymax=ymax)) +
  geom_point() + 
  geom_errorbar(width=.2) +
  geom_line() +
  geom_hline(yintercept=0.2, lty=2) +
  scale_y_continuous(breaks=c(1:5/10), limits=c(0,.45))
dev.off()
