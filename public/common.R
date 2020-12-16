getTrimmedSumStats <- function(dir, tsv.filename, ld.filename, r2_threshold) {

  # read in the summary stats and LD matrices
  tsv_files <- scan(tsv.filename, what="char")
  ld_files <- scan(ld.filename, what="char")
  stopifnot(length(tsv_files) == length(ld_files))
  nclumps <- length(tsv_files)
  ld_mat <- list()
  sum_stat <- list()
  cc <- c("character", "character", "numeric",
          "character", "character", "numeric", "numeric", # GWAS
          "character", "character", "numeric", "numeric", # eQTL
          "character", "character") # plink
  for (j in 1:nclumps) {
    filename <- file.path(dir,tsv_files[j])
    out <- read.delim(filename, header=TRUE, colClasses=cc)
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

  trim_clusters <- clusterTrimmer(r2, r2_threshold=r2_threshold)

  if (length(trim_clusters) > 0) {
    # trim those clumps
    sum_stat <- sum_stat[-trim_clusters]
    ld_mat <- ld_mat[-trim_clusters]
  }

  list(sum_stat=sum_stat, ld_mat=ld_mat, trim_clusters=trim_clusters)
}
