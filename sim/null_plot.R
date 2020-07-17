i <- "null"

files <- sub(".final","",list.files("out", pattern="final"))
files <- grep(paste0("^",i),files,value=TRUE)
files

final <- list()
for (k in seq_along(files)) {
  final[[k]] <- read.table(paste0("out/",files[k],".final"),skip=1)
}

mrlocus <- list()
for (k in seq_along(files)) {
  mrlocus[[k]] <- read.table(paste0("out/",files[k],".mrlocus"))
}
mrlocus10 <- sapply(mrlocus, function(x) x[2,1])
mrlocus90 <- sapply(mrlocus, function(x) x[2,2])

h2 <- 100*as.numeric(sub(".*_(.*)h2_.*","\\1",files))

idx <- c(3,5,6:8)
dat <- data.frame(rep=rep(rep(1:20,each=5),3),
                  h2=h2,
                  true=rep(sapply(final, function(x) x[1,2]), each=5),
                  method=rep(c("causal","all","twmr","ptwas","mrlocus"),times=20),
                  estimate=as.vector(sapply(final, function(x) x$V2[idx])),
                  se=as.vector(sapply(final, function(x) x$V3[idx])))
dat$method <- factor(dat$method, c("causal","all","twmr","ptwas","mrlocus"))
dat$est_nozero <- ifelse(dat$est == 0, NA, dat$est)
q <- qnorm(.9)
dat$min <- dat$estimate - q * dat$se
dat$max <- dat$estimate + q * dat$se
dat$min[dat$method == "mrlocus"] <- mrlocus10
dat$max[dat$method == "mrlocus"] <- mrlocus90
dat$h2 <- factor(dat$h2)
levels(dat$h2) <- paste0("h2: ",levels(dat$h2),"%")

dat$contain <- dat$true > dat$min & dat$true < dat$max & !is.na(dat$est_nozero)

library(dplyr)
tab <- dat %>% group_by(h2, method) %>%
  summarize(cov=paste0(100*mean(contain),"%"))
tab
tab$x <- "left"
tab$y <- "top"

library(ggplot2)
library(ggpmisc)
png(file="~/Desktop/nullplot.png", res=125, width=1200, height=800)
ggplot(dat, aes(estimate,rep,xmin=min,xmax=max,color=contain)) +
  geom_pointrange(shape="square",size=.5) + facet_grid(h2 ~ method) +
  geom_vline(xintercept=0) +
  scale_color_manual(values=c(2,1)) +
  geom_text_npc(data=tab, aes(npcx=x, npcy=y, label=cov))
dev.off()

