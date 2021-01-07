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

# r2 of kept clumps
r2 <- lapply(idx, function(i) {
  out <- filesvec(i, "*\\.mrl_r2")
  dir <- out$dir; files <- out$files
  unlist(lapply(files, function(f) {
    scan(file.path(dir,f), quiet=TRUE)
  }))
})

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

#save(i2, causal, clumps, dap, kept, r2, mrl.kept, ecav.kept, file="sim_review.rda")

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
#write.table(key, file="sim_review.tsv", row.names=FALSE, sep="\t")
library(ggplot2)
cols <- palette.colors(4, palette="Set 2")

# sim types plot
pdf(file="../supp/figs/sim_types.pdf", height=5, width=8)
ggplot(key, aes(sim, value, fill=group)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  facet_wrap(~type, nrow=2, scales="free") +
  scale_fill_manual(values=cols) + ylab("")
dev.off()

# subset
key2 <- key[1:13,c("sim","id","group")]

library(ggbeeswarm)
plotit <- function(x, dot=FALSE, bee=TRUE, noise=TRUE) {
  dat <- data.frame(number=unlist(x),
                    id=rep(names(x), lengths(x)))
  if (noise) {
    dat$number <- dat$number + runif(nrow(dat),-.25,.25)
  }
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

#plotit(i2, dot=FALSE, bee=TRUE) + ggtitle("I2 values for PTWAS")

# r2 of kept clusters plot
pdf(file="../supp/figs/sim_cluster_r2.pdf", height=5, width=8)
plotit(r2, noise=FALSE) + ggtitle("Pairwise r2 of clusters provided to eCAVIAR/MRLocus colocalization") + ylab("r2") + scale_y_log10() +
  geom_hline(yintercept=0.05, color="red", lty=2)
dev.off()

scl <- scale_y_continuous(limits=c(0,10), breaks=0:5*2)
twomore <- function(x) x[x >= 2]

library(patchwork)
pdf(file="../supp/figs/sim_details.pdf", width=14, height=14)
p1 <- plotit(causal) + ylim(0,15) + ggtitle("# of true eQTL SNPs")
p2 <- plotit(clumps) + ylim(0,15) + ggtitle("# of PLINK clumps")
p3 <- plotit(lapply(dap, na.omit)) + ylim(0,40) + ggtitle("# of DAP signal clusters")
p4 <- plotit(lapply(kept,twomore)) + scl + ggtitle("# of clusters passing 1st round LD-based trimming")
p5 <- plotit(lapply(mrl.kept,twomore)) + scl + ggtitle("# of clusters passing 2nd round LD-based trimming - MRLocus")
p6 <- plotit(lapply(ecav.kept,twomore)) + scl + ggtitle("# of clusters passing 2nd round LD-based trimming - eCAVAIR")
(p1 | p2) / (p3 | p4) / (p5 + p6) + plot_annotation(tag_levels="A")
dev.off()

###

# how much are eQTL beta over-estimated?
overest <- lapply(idx, function(i) {
  out <- filesvec(i, "scan.tsv$")
  dir <- out$dir; files <- out$files
  z.thr <- qnorm(.001/2, lower.tail=FALSE)
  oe <- lapply(files, function(f) {
    x <- read.table(file.path(dir, f), header=TRUE)
    x <- x[which.max(abs(x$eqtl.true)),]
    x$abs.z <- with(x, abs(eqtl.beta/eqtl.se))
    x <- x[x$abs.z > z.thr,]
    if (nrow(x) > 0) {
      pmax(x$eqtl.beta/x$eqtl.true, 0)
    } else {
      NULL
    }
  })
  unlist(oe)
})

library(ggplot2)
key3 <- key[1:13,]
dat <- data.frame(ratio=unlist(overest),
                  id=factor(rep(names(overest), lengths(overest))))
dat$h2g <- key3$value[match(dat$id, key2$id)]
pdf(file="../supp/figs/sim_overest.pdf", width=10, height=5)
ggplot(dat, aes(ratio)) +
  geom_histogram(fill="dodgerblue", alpha=0.5, bins=20) +
  geom_vline(xintercept=1, lty=2, col="red") +
  scale_x_log10() +
  xlab("eQTL coefficient estimate:true ratio") +
  facet_wrap(~h2g)
dev.off()
