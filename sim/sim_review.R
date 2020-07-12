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
save(i2, causal, clumps, dap, file="sim_review.rda")

library(ggplot2)
library(ggbeeswarm)
plotit <- function(x, dot=TRUE, bee=FALSE) {
  dat <- data.frame(number=unlist(x),
                    sim=rep(names(x), lengths(x)))
  g <- ggplot(dat, aes(sim, number)) +
    geom_violin() +
    theme_bw()
  if (dot) {
    g <- g + geom_dotplot(binaxis="y",
                          stackdir="center",
                          dotsize=.5,
                          stackratio=1.25)
  }
  if (bee) {
    g <- g + geom_quasirandom()
  }
  g
}

plotit(i2, dot=FALSE, bee=TRUE)
plotit(causal) + ylim(0,15)
plotit(clumps) + ylim(0,20)
plotit(dap, dot=FALSE, bee=TRUE) + ylim(0,40)
