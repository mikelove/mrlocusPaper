key <- read.table("sim_review.tsv", header=TRUE)
key <- key[1:13,2:3]
key <- rbind(key[1,,drop=FALSE], data.frame(sim="HP", id="hp"), key[2:13,])

library(pbapply)

key$mrlocus <- NA
key$ecaviar <- NA
for (type in key$id) {
  print(paste("---",type,"---"))
  dir <- file.path("out",type)
  files <- list.files(dir, pattern="ve.clumped")
  out <- pblapply(seq_along(files), function(i) {
    file <- sub(".clumped","",files[i])
    clumped.filename <- paste0(dir,"/",file,".clumped")
    scan.filename <- paste0(dir,"/",file,".scan.tsv")
    clumped <- read.table(clumped.filename, strings=FALSE, header=TRUE)
    big_sum_stat <- read.delim(scan.filename, strings=FALSE)
    noparen <- function(z) sub("\\(1\\)","",z)
    clumps <- lapply(strsplit(clumped$SP2,split=","), noparen)
    # add index SNP to clumps
    clumps <- lapply(seq_along(clumps), function(j) c(clumped$SNP[j], clumps[[j]]))
    sum_stat <- lapply(clumps, function(x) {
      out <- big_sum_stat[big_sum_stat$snp %in% x,]
      out$z <- out$eqtl.beta/out$eqtl.se
      out$abs.z <- abs(out$z)
      out$z2 <- out$gwas.beta/out$gwas.se # for testing ecaviar
      out
    })
    # load coloc results
    load(paste0(dir,"/",file,".mrl_coloc")) # MRLocus coloc results - `res`
    load(paste0(dir,"/",file,".ecav")) # eCAVIAR coloc results - `ecav.coloc`
    # trim `sum_stat` and `ecav.coloc` to remove high r2 clumps
    # (these are already trimmed from MRLocus `res`)
    keep <- scan(paste0(dir,"/",file,".mrl_keep"), quiet=TRUE)
    sum_stat <- sum_stat[keep]
    ecav.coloc <- ecav.coloc[keep]
    # any true eSNPs in the clumps
    in.clump <- sapply(sum_stat, function(x) any(x$eqtl.true != 0))
    if (any(in.clump)) {
      mrlocus <- sapply(which(in.clump), function(j) {
        eqtl.true <- sum_stat[[j]]$eqtl.true != 0
        eqtl.nms <- sum_stat[[j]]$snp[eqtl.true]
        snps <- res$alleles[[j]]$collapsed[ which.max(res$beta_hat_a[[j]]) ]
        any(sapply(eqtl.nms, function(p) grepl(p, snps)))
      })
      ecaviar <- sapply(which(in.clump), function(j) {
        eqtl.true <- sum_stat[[j]]$eqtl.true != 0
        eqtl.nms <- sum_stat[[j]]$snp[eqtl.true]
        snp <- ecav.coloc[[j]]$SNP_ID[which.max(ecav.coloc[[j]]$CLPP)]
        # get all correlated SNP IDs
        snps <- grep(snp, res$alleles[[j]]$collapsed, value=TRUE)
        any(sapply(eqtl.nms, function(p) grepl(p, snps)))
      })
      data.frame(mrlocus=mrlocus, ecaviar=ecaviar)
    } else {
      NULL
    }
  })
  dat <- do.call(rbind, out)
  key[key$id == type,c("mrlocus","ecaviar")] <- colMeans(dat)
}
#save(key, file="assess_coloc.rda")
load("assess_coloc.rda")

library(reshape)
dat <- melt(key, id=c("sim","id"), variable_name="method")
names(dat)[4] <- "accuracy"
dat$sim <- factor(dat$sim, levels=unique(dat$sim))
library(ggplot2)
pdf(file="../supp/figs/assess_coloc.pdf", width=10, height=5)
ggplot(dat, aes(sim, accuracy, fill=method)) +
  geom_bar(stat="identity", position="dodge", width=.66) +
  ylim(0,1) + ylab("accuracy in colocalization")
dev.off()
