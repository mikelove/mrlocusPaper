type <- "high_n"
dir <- file.path("out",type)
files <- list.files(dir, pattern="ve.clumped")

library(pbapply)
out <- pblapply(1:20, function(i) {
  file <- sub(".clumped","",files[i])
  clumped.filename <- paste0(dir,"/",file,".clumped")
  scan.filename <- paste0(dir,"/",file,".scan.tsv")
  ld.filename <- paste0(dir,"/",file,".ld")
  clumped <- read.table(clumped.filename, strings=FALSE, header=TRUE)
  big_sum_stat <- read.delim(scan.filename, strings=FALSE)
  big_ld_mat <- as.matrix(read.table(ld.filename))
  noparen <- function(z) sub("\\(1\\)","",z)
  clumps <- lapply(strsplit(clumped$SP2,split=","), noparen)
  # add index SNP to clumps
  clumps <- lapply(seq_along(clumps), function(j) c(clumped$SNP[j], clumps[[j]]))
  ld_mat <- lapply(clumps, function(x) {
    idx <- big_sum_stat$snp %in% x
    unname(big_ld_mat[idx,idx,drop=FALSE])
  })
  sapply(ld_mat, nrow)
  sum_stat <- lapply(clumps, function(x) {
    out <- big_sum_stat[big_sum_stat$snp %in% x,]
    out$z <- out$eqtl.beta/out$eqtl.se
    out$abs.z <- abs(out$z)
    out$z2 <- out$gwas.beta/out$gwas.se # for testing ecaviar
    out
  })
  load(paste0(dir,"/",file,".mrl_coloc")) # res
  load(paste0(dir,"/",file,".ecav")) # ecav.coloc
  in.clump <- sapply(sum_stat, function(x) any(x$eqtl.true != 0))
  if (any(in.clump)) {
    mrlocus <- sapply(which(in.clump), function(j) {
      eqtl.true <- sum_stat[[j]]$eqtl.true != 0
      eqtl.nms <- sum_stat[[j]]$snp[eqtl.true]
      snps <- res$alleles[[j]]$collapsed[ which.max(res$beta_hat_a[[j]]) ]
      any(sapply(eqtl.nms, function(p) grepl(p, snps)))
    })
    eCAVIAR <- sapply(which(in.clump), function(j) {
      eqtl.true <- sum_stat[[j]]$eqtl.true != 0
      eqtl.nms <- sum_stat[[j]]$snp[eqtl.true]
      snp <- ecav.coloc[[j]]$SNP_ID[which.max(ecav.coloc[[j]]$CLPP)]
      # get all correlated SNP IDs
      snps <- grep(snp, res$alleles[[j]]$collapsed, value=TRUE)
      any(sapply(eqtl.nms, function(p) grepl(p, snps)))
    })
    data.frame(mrlocus=mrlocus, eCAVIAR=eCAVIAR)
  } else {
    NULL
  }
})

dat <- do.call(rbind, out)
colMeans(dat)
