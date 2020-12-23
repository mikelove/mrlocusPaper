filesvec <- function(i, ending) {
  if (i < 10) {
    dir <- file.path("out",i)
  } else if (i < 13) {
    ii <- i - 9
    dir <- file.path("out",paste0("null",ii))
  } else {
    dir <- file.path("out","high_n")
  }
  files <- list.files(dir, pattern=ending)
  list(dir=dir, files=files)
}

idx <- c(1,13,2:12)
nms <- c(1,"high_n",2:9,paste0("null",1:3))
names(idx) <- nms

i2 <- lapply(idx, function(i) {
  out <- filesvec(i, "*\\.ptwas")
  dir <- out$dir; files <- out$files
  sapply(files, USE.NAMES=FALSE, function(f) {
    x <- scan(file.path(dir,f), what="char", quiet=TRUE)
    if (all(x[(length(x)-5):length(x)] == 0)) {
      NA
    } else {
      as.numeric(x[length(x)])
    }
  })
})

causal <- lapply(idx, function(i) {
  out <- filesvec(i, "*\\.scan.tsv")
  dir <- out$dir; files <- out$files
  sapply(files, USE.NAMES=FALSE, function(f) {
    x <- read.delim(file.path(dir,f))
    sum(x$eqtl.true != 0)
  })
})

clumps <- lapply(idx, function(i) {
  out <- filesvec(i, "*\\.clumped")
  dir <- out$dir; files <- out$files
  files <- grep("p1e-4", files, value=TRUE, invert=TRUE)
  nlines <- sapply(files, USE.NAMES=FALSE, function(f) {
    x <- readLines(file.path(dir, f))
    x <- x[x != ""]
    length(x)
  })
  nlines - 1 # remove header
})

dap <- lapply(idx, function(i) {
  out <- filesvec(i, "*\\.dap$")
  dir <- out$dir; files <- out$files
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

# number of clumps kept in round one
kept <- lapply(idx, function(i) {
  out <- filesvec(i, "*\\.mrl_keep$")
  dir <- out$dir; files <- out$files
  unname(sapply(files, function(f) length(scan(file.path(dir, f), quiet=TRUE))))
})
lengths(kept)
sapply(kept, function(x) sum(x[1:20] > 1)) # 17 or less => more reps
sapply(kept, function(x) sum(x > 1))

# number of clumps kept in round two
mrl.kept <- lapply(idx, function(i) {
  out <- filesvec(i, "*\\.mrl_keep2$")
  dir <- out$dir; files <- out$files
  unname(sapply(files, function(f) length(scan(file.path(dir, f), quiet=TRUE))))
})
lengths(mrl.kept)

# number of clumps kept in round two
ecav.kept <- lapply(idx, function(i) {
  out <- filesvec(i, "ecav-mrl_keep2$")
  dir <- out$dir; files <- out$files
  unname(sapply(files, function(f) length(scan(file.path(dir, f), quiet=TRUE))))
})
lengths(ecav.kept)

#save(i2, causal, clumps, dap, kept, mrl.kept, ecav.kept, file="sim_review.rda")

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

h2 <- c(c(.1,.1,.2,.05),rep(c(.1,.2,.05),3))
ve <- c(rep(c(.01,.005,.001),c(4,3,3)),c(0,0,0))
sim <- c("A","High-N",LETTERS[2:9],paste0("Null-",c(.1,.2,.05)))
key <- data.frame(value=c(h2,ve),
                  sim=factor(rep(sim,2),levels=sim),
                  id=rep(c(1,"high_n",3,2, 4,6,8, 5,7,9, "null1","null2","null3"),2),
                  type=rep(c("eQTL h2","variance explained"),each=13),
                  group=factor(rep(rep(1:4,c(4,3,3,3)),2)))
write.table(key, file="sim_review.tsv", row.names=FALSE, sep="\t")
library(ggplot2)
cols <- palette.colors(4, palette="Set 2")

pdf(file="../supp/figs/sim_types.pdf", height=5, width=8)
ggplot(key, aes(sim, value, fill=group)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  facet_wrap(~type, nrow=2, scales="free") +
  scale_fill_manual(values=cols) + ylab("")
dev.off()

key2 <- key[1:13,c("sim","id","group")]

library(ggbeeswarm)
plotit <- function(x, dot=FALSE, bee=TRUE) {
  dat <- data.frame(number=unlist(x),
                    id=rep(names(x), lengths(x)))
  dat$number <- dat$number + runif(nrow(dat),-.25,.25)
  idx <- match(dat$id, key2$id)
  dat$sim <- key2$sim[idx]
  dat$group <- key2$group[idx]
  dat$sim <- factor(dat$sim, levels=key2$sim)
  # shorten sim names
  levels(dat$sim) <- c("A", "Hi-N", LETTERS[2:9], c("N.1","N.2","N.05"))
  g <- ggplot(dat, aes(sim, number, fill=group)) +
    geom_violin(show.legend=FALSE) +
    scale_fill_manual(values=cols)
  if (dot) {
    g <- g + geom_dotplot(binaxis="y", stackdir="center", dotsize=.5,
                          fill="black", show.legend=FALSE)
  }
  if (bee) {
    g <- g + geom_quasirandom(show.legend=FALSE, groupOnX=TRUE)
  }
  g
}

plotit(i2, dot=FALSE, bee=TRUE) + ggtitle("I2 values for PTWAS")

scl <- scale_y_continuous(limits=c(0,10), breaks=0:5*2)

library(patchwork)
pdf(file="../supp/figs/sim_details.pdf", width=14, height=14)
p1 <- plotit(causal) + ylim(0,15) + ggtitle("# of true eQTL SNPs")
p2 <- plotit(clumps) + ylim(0,15) + ggtitle("# of PLINK clumps")
p3 <- plotit(lapply(dap, na.omit)) + ylim(0,40) + ggtitle("# of DAP signal clusters")
p4 <- plotit(kept) + scl + ggtitle("# of clusters passing 1st round trim")
p5 <- plotit(mrl.kept) + scl + ggtitle("# of clusters passing 2nd round trim - MRLocus")
p6 <- plotit(ecav.kept) + scl + ggtitle("# of clusters passing 2nd round trim - eCAVAIR")
(p1 | p2) / (p3 | p4) / (p5 + p6) + plot_annotation(tag_levels="A")
dev.off()
