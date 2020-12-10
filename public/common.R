getTrimmedSumStats <- function(dir, tsv.filename, ld.filename) {

  # read in the summary stats and LD matrices
  tsv_files <- scan(tsv.filename, what="char")
  ld_files <- scan(ld.filename, what="char")
  stopifnot(length(tsv_files) == length(ld_files))
  nclumps <- length(tsv_files)
  ld_mat <- list()
  sum_stat <- list()
  for (j in 1:nclumps) {
    filename <- file.path(dir,tsv_files[j])
    out <- read.table(filename, header=TRUE)
    out$z <- out$beta_eQTL/out$se_eQTL
    out$abs.z <- abs(out$z)
    sum_stat[[j]] <- out
    filename <- file.path(dir,ld_files[j])
    ld_mat[[j]] <- as.matrix(read.table(filename))
  }

  # mapping from TSV files to the indexLD matrix
  info <- read.table(list.files(dir, pattern="indexinfo", full=TRUE), header=TRUE)

  # matrix of LD for the index SNPs
  ld <- as.matrix(read.table(list.files(dir, pattern="indexLD", full=TRUE)))
  ld <- ld[info$idx, info$idx] # reorder according to indexinfo
  r2 <- ld^2

  # check order of Z stat (decreasing)
  z_stat <- sapply(seq_along(sum_stat), function(i) sum_stat[[i]][sum_stat[[i]]$SNP == info$idxSNP[i],"abs.z"])
  stopifnot(all(z_stat == sort(z_stat, decreasing=TRUE)))

  # find clumps that have pairwise correlation with other clumps above a threshold
  r2.threshold <- 0.01
  trim.clumps <- c()
  if (nclumps > 1) {
    diag(r2) <- 0 # useful for logic below
    if (nclumps > 1 & any(r2 > r2.threshold)) {
      for (j in nclumps:2) {
        if (any( (r2[j,] > r2.threshold)[!1:nclumps %in% trim.clumps] )) {
          trim.clumps <- c(trim.clumps, j)
        }
      }
    }
    trim.clumps <- rev(trim.clumps)
  }

  if (length(trim.clumps) > 0) {
    # trim those clumps
    sum_stat <- sum_stat[-trim.clumps]
    ld_mat <- ld_mat[-trim.clumps]
  }

  list(sum_stat=sum_stat, ld_mat=ld_mat, trim.clumps=trim.clumps)
}
