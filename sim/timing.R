myscan <- function(x) as.numeric(scan(x, what="char", quiet=TRUE)[10])
mylistfiles <- function(x) {
  c(list.files("bench/1", x, full.names=TRUE),
    list.files("bench/high_n", x, full.names=TRUE))
  }

mrl <- mylistfiles("_mrlocus")
mrl.time <- sapply(mrl, myscan, USE.NAMES=FALSE)
length(mrl.time)
hist(mrl.time)

ecav <- mylistfiles("_ecav.bench")
ecav.time <- sapply(ecav, myscan, USE.NAMES=FALSE)
length(ecav.time)
hist(ecav.time)

emrl <- mylistfiles("ecav-mrlocus")
emrl.time <- sapply(emrl, myscan, USE.NAMES=FALSE)
length(emrl.time)
hist(emrl.time)

ptw <- mylistfiles("ptwas")
ptw.time <- sapply(ptw, myscan, USE.NAMES=FALSE)
length(ptw.time)
hist(ptw.time)

twmr <- mylistfiles("twmr")
twmr.time <- sapply(twmr, myscan, USE.NAMES=FALSE)
length(twmr.time)
hist(twmr.time)

#save(twmr.time,ptw.time,mrl.time, file="timing.rda")
load("timing.rda")

library(ggplot2)
mths <- c("twmr","ptwas","ecaviar","ecaviar-mrlocus","mrlocus")
dat <- data.frame(seconds=c(twmr.time,ptw.time,
                            ecav.time,emrl.time,mrl.time),
                  method=factor(rep(mths,each=80),levels=mths))

library(dplyr)
ms <- dat %>% group_by(method) %>%
  summarize(mean=mean(seconds)) %>%
  pull(mean)

tab <- data.frame(
  method=mths,
  seconds=c(10,10,10,60,300),
  label=paste("mean:",round(ms,2),"s"))

#pdf(file="../supp/figs/runtime.pdf", width=4, height=4)

ggplot(dat, aes(method,seconds)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1,300), breaks=c(1,2,5,10,30,60,120,300)) +
  ylab("runtime (seconds)") +
  geom_text(aes(method,seconds,label=label),
            data.frame(method=rep("twmr",3),
                       seconds=c(60,120,300),
                       label=c("1 min","2 min","5 min")),
            fontface="bold") +
  geom_label(aes(method,seconds,label=label),
             fill="cornsilk",
             tab)

#dev.off()

mean(mrl.time)
mean(twmr.time)
mean(ptw.time)
mean(mrl.time)/mean(twmr.time)
mean(mrl.time)/mean(ptw.time)
