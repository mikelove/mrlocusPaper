myscan <- function(x) as.numeric(scan(x, what="char", quiet=TRUE)[10])

mrl <- list.files("bench","mrlocus", full.names=TRUE)
mrl.time <- sapply(mrl, myscan)
hist(log10(mrl.time))

ptw <- list.files("bench","ptwas", full.names=TRUE)
ptw.time <- sapply(ptw, myscan)
hist(log10(ptw.time))

twmr <- list.files("bench","twmr", full.names=TRUE)
twmr.time <- sapply(twmr, myscan)
hist(log10(twmr.time))

#save(twmr.time,ptw.time,mrl.time, file="timing.rda")
load("timing.rda")

library(ggplot2)
mths <- c("twmr","ptwas","mrlocus")
dat <- data.frame(seconds=c(twmr.time,ptw.time,mrl.time),
                  method=factor(rep(mths,each=240),levels=mths))

pdf(file="../supp/figs/runtime.pdf", width=4, height=4)

ggplot(dat, aes(method,seconds)) +
  geom_boxplot() +
  scale_y_log10(breaks=c(2,5,10,30,60,120,300,600)) +
  ylab("runtime (seconds)") +
  geom_text(aes(method,seconds,label=label),
            data.frame(method=rep("twmr",4),
                       seconds=c(60,120,300,600),
                       label=c("1 min","2 min","5 min","10 min")),
            fontface="bold") +
  geom_label(aes(method,seconds,label=label),
             fill="cornsilk",
             data.frame(
               method=c("twmr","ptwas","mrlocus"),
               seconds=c(15,15,300),
               label=paste("mean:",
                           round(sapply(list(twmr.time,ptw.time,mrl.time), mean),2),
                           "s")))

dev.off()

mean(mrl.time)
mean(twmr.time)
mean(ptw.time)
mean(mrl.time)/mean(twmr.time)
mean(mrl.time)/mean(ptw.time)
