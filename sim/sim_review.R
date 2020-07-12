i2 <- lapply(1:12, function(i) {
  if (i < 10) {
    dir <- file.path("out",i)
    files <- list.files(dir, pattern="ptwas")
  } else {
    ii <- i - 9
    dir <- file.path("out","null")
    files <- list.files(dir, pattern=paste0("null",ii,".*ptwas"))
  }
  sapply(files, USE.NAMES=FALSE, function(f) {
    x <- scan(file.path(dir,f), what="char", quiet=TRUE)
    if (all(x[(length(x)-5):length(x)] == 0)) {
      NA
    } else {
      as.numeric(x[length(x)])
    }
  })
})
names(i2) <- c(1:9,paste0("null",1:3))

causal <- lapply(1:12, function(i) {
  if (i < 10) {
    dir <- file.path("out",i)
    files <- list.files(dir, pattern="scan.tsv")
  } else {
    ii <- i - 9
    dir <- file.path("out","null")
    files <- list.files(dir, pattern=paste0("null",ii,".*scan.tsv"))
  }
  sapply(files, USE.NAMES=FALSE, function(f) {
    x <- read.delim(file.path(dir,f))
    sum(x$eqtl.true != 0)
  })
})
names(causal) <- c(1:9,paste0("null",1:3))

clumps <- lapply(1:12, function(i) {
  if (i < 10) {
    dir <- file.path("out",i)
    files <- list.files(dir, pattern="clumped")
  } else {
    ii <- i - 9
    dir <- file.path("out","null")
    files <- list.files(dir, pattern=paste0("null",ii,".*clumped"))
  }
  sapply(files, USE.NAMES=FALSE, function(f) {
    length(readLines(file.path(dir, f))) - 1
  })
})
names(clumps) <- c(1:9,paste0("null",1:3))

dap <- lapply(1:12, function(i) {
  if (i < 10) {
    dir <- file.path("out",i)
    files <- list.files(dir, pattern="dap$")
  } else {
    ii <- i - 9
    dir <- file.path("out","null")
    files <- list.files(dir, pattern=paste0("null",ii,".*dap$"))
  }
  sapply(files, USE.NAMES=FALSE, function(f) {
    x <- scan(file.path(dir, f), what="char", sep="\n", quiet=TRUE)
    if ("Independent association signal clusters" %in% x) {
      idx <- which(x == "Independent association signal clusters")
      x <- x[(idx+2):length(x)]
      length(x)
    } else {
      NA
    }
  })
})
names(dap) <- c(1:9,paste0("null",1:3))
#save(i2, causal, clumps, dap, file="sim_review.rda")

load("sim_review.rda")

.1 .01 - 1
.05 .01 - 2 
.2 .01 - 3
.1 .005 - 4 
.1 .001 - 5 
.2 .005 - 6 
.2 .001 - 7 
.05 .005 - 8 
.05 .001 - 9

h2 <- rep(c(.1,.2,.05),4)
ve <- c(rep(c(.01,.005,.001),each=3),c(0,0,0))
sim <- c(LETTERS[1:9], paste0("Null-",c(.1,.2,.05)))
key <- data.frame(value=c(h2,ve),
                  sim=factor(rep(sim,2),levels=sim),
                  id=rep(c(1,3,2, 4,6,8, 5,7,9, "null1","null2","null3"), 2),
                  type=rep(c("eQTL h2","variance explained"),each=12),
                  group=factor(rep(rep(1:4,each=3),2)))

library(ggplot2)
cols <- palette.colors(4, palette="Set 2")
ggplot(key, aes(sim, value, fill=group)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  facet_wrap(~type, nrow=2, scales="free") +
  scale_fill_manual(values=cols)

key2 <- key[1:12,c("sim","id","group")]

library(ggbeeswarm)
plotit <- function(x, dot=TRUE, bee=FALSE) {
  dat <- data.frame(number=unlist(x),
                    id=rep(names(x), lengths(x)))

  idx <- match(dat$id, key2$id)
  dat$sim <- key2$sim[idx]
  dat$group <- key2$group[idx]
  dat$sim <- factor(dat$sim, levels=key2$sim)

  g <- ggplot(dat, aes(sim, number, fill=group)) +
    geom_violin(show.legend=FALSE) +
    theme_bw() +
    scale_fill_manual(values=cols)
  
  if (dot) {
    g <- g + geom_dotplot(binaxis="y",
                          stackdir="center",
                          dotsize=.5,
                          stackratio=1.25,
                          fill="black",
                          show.legend=FALSE)
  }
  if (bee) {
    g <- g + geom_quasirandom(show.legend=FALSE)
  }
  g
}

plotit(i2, dot=FALSE, bee=TRUE)
plotit(causal) + ylim(0,15)
plotit(clumps) + ylim(0,20)
plotit(dap, dot=FALSE, bee=TRUE) + ylim(0,40)
