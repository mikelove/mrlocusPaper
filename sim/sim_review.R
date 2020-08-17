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
  nlines <- sapply(files, USE.NAMES=FALSE, function(f) {
    x <- readLines(file.path(dir, f))
    x <- x[x != ""]
    length(x)
  })
  nlines - 1 # remove header
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

## h2 - v.e. - sim number - sim letter 
## -----------------------------------
## .1  - .01  - 1 - A
## .05 - .01  - 2 - C
## .2  - .01  - 3 - B

## .1  - .005 - 4 - D
## .1  - .001 - 5 - G

## .2  - .005 - 6 - E
## .2  - .001 - 7 - H

## .05 - .005 - 8 - F
## .05 - .001 - 9 - I

## .1  - 0    - null1
## .2  - 0    - null2
## .05 - 0    - null3

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

pdf(file="../supp/figs/sim_types.pdf", height=5)
ggplot(key, aes(sim, value, fill=group)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  facet_wrap(~type, nrow=2, scales="free") +
  scale_fill_manual(values=cols) + ylab("")
dev.off()

key2 <- key[1:12,c("sim","id","group")]

library(ggbeeswarm)
plotit <- function(x, dot=TRUE, bee=FALSE) {
  dat <- data.frame(number=unlist(x),
                    id=rep(names(x), lengths(x)))

  idx <- match(dat$id, key2$id)
  dat$sim <- key2$sim[idx]
  dat$group <- key2$group[idx]
  dat$sim <- factor(dat$sim, levels=key2$sim)

  # shorten sim names
  levels(dat$sim) <- c(LETTERS[1:9], c("N.1","N.2","N.05"))
  
  g <- ggplot(dat, aes(sim, number, fill=group)) +
    geom_violin(show.legend=FALSE) +
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

plotit(i2, dot=FALSE, bee=TRUE) + ggtitle("I2 values for PTWAS")

library(patchwork)
pdf(file="../supp/figs/sim_details.pdf", width=14, height=5)
p1 <- plotit(causal) + ylim(0,15) + ggtitle("number of true eQTL SNPs per simuation")
p2 <- plotit(clumps) + ylim(0,15) + ggtitle("number of PLINK clumps per simulation")
p3 <- plotit(dap, dot=FALSE, bee=TRUE) + ylim(0,40) + ggtitle("number of DAP signal clusters per simulation")
p1 + p2 + p3 + plot_annotation(tag_levels="A")
dev.off()
