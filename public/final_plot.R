dat <- read.delim("estimates.tsv")
library(ggplot2)
library(dplyr)
q <- qnorm(.9)
dat <- dat %>% mutate(xmin=alpha-q*se,xmax=alpha+q*se)

# replace SD based intervals with quantile based
mrlocusFiles <- list.files(pattern=".txt")
ints <- lapply(mrlocusFiles, function(f) read.table(f, header=TRUE)[1,3:4])
ints <- do.call(rbind, ints)
dat[dat$method == "mrlocus",c("xmin","xmax")] <- ints

dat$label <- paste0(dat$gene," (",dat$tissue,")", " → ", dat$trait)
dat$label <- factor(dat$label, levels=dat$label[1:5])

png(file="../supp/figs/forest.png", width=1200, height=600, res=150)
ggplot(dat, aes(x=alpha,y=method,xmin=xmin,xmax=xmax)) +
  geom_pointrange() +
  facet_wrap(~label, scales="free") +
  geom_vline(xintercept = 0, lty=2) +
  xlab("gene-to-trait estimate")
dev.off()
