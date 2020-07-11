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
stripchart(i2, ylab="I^2", xlab="simulation", main="PTWAS I^2", method="jitter", pch=1, vertical=TRUE)
abline(v=1:12, col=rgb(0,0,0,.3), lty=2)
