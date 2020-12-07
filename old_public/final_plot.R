library(ggplot2)
library(dplyr)

dat <- read.delim("estimates.tsv")
q <- qnorm(.9)
dat <- dat %>% mutate(xmin=alpha-q*se,xmax=alpha+q*se)

# replace SD based intervals with quantile based
mrlocusFiles <- list.files(pattern=".txt")
ints <- lapply(mrlocusFiles, function(f) read.table(f, header=TRUE)[1,3:4])
ints <- do.call(rbind, ints)
dat[dat$method == "mrlocus",c("xmin","xmax")] <- ints

dat$label <- paste0(dat$gene," (",dat$tissue,")", " → ", dat$trait)

dat <- rbind(dat, data.frame(tissue=NA,gene=NA,trait=NA,method="twmr",
                             alpha=NA,se=NA,I2=NA,sigma=NA,xmin=NA,xmax=NA,
                             label="heterogeneity"))

lvls <- c(dat$label[1:5],"heterogeneity")
dat$label <- factor(dat$label, levels=lvls)

library(ggpmisc)
tab <- cbind(dat[dat$method == "ptwas",c("gene","I2")], dat[dat$method == "mrlocus","sigma",drop=FALSE])
names(tab) <- c("gene","I^2", "sigma")
data.tb <- tibble(x=0, y="twmr",
                  label=factor("heterogeneity", levels=lvls),
                  tb=list(tab))

png(file="../supp/figs/forest.png", width=1200, height=600, res=150)
ggplot(dat, aes(alpha,method,xmin=xmin,xmax=xmax)) +
  geom_pointrange() +
  facet_wrap(~label, scales="free") +
  geom_vline(xintercept = 0, lty=2) +
  xlab("gene-to-trait estimate") +
  geom_table(data=data.tb, aes(x, y, label=tb),
             table.theme = ttheme_gtlight,
             stat="fmt_tb", size=4)
dev.off()
