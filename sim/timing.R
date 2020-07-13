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

boxplot(list(twmr.time,ptw.time,mrl.time),
        names=c("TWMR","PTWAS","MRLocus"),
        log="y",ylim=c(1,1000),range=0, ylab="seconds")
abline(h=c(1,5,10,30,60,120,300,600),lty=2,col=rgb(0,0,0,.5))
text(c(1,1,1,1,1),c(30,60,120,300,600),
     c("30 sec","1 min","2 min","5 min","10 min"),pos=3)

mean(mrl.time)
mean(twmr.time)
mean(ptw.time)
mean(mrl.time)/mean(twmr.time)
mean(mrl.time)/mean(ptw.time)
