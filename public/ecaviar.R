cmd_args=commandArgs(TRUE)

ecav.bin <- cmd_args[1]
dir <- cmd_args[2]
tsv.filename <- cmd_args[3]
ld.filename <- cmd_args[4]
out.filename <- cmd_args[5]

ecav.bin <- "/proj/milovelab/love/bin/caviar/caviar/eCAVIAR"
dir <- "Artery_MRAS_CAD"
tsv.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.tsv"
ld.filename <- "Artery_MRAS_CAD/Artery_MRAS_CAD.ld"

tsv_files <- scan(tsv.filename, what="char")
ld_files <- scan(ld.filename, what="char")
stopifnot(length(tsv_files) == length(ld_files))
info <- read.table(list.files(dir, pattern="indexinfo", full=TRUE), header=TRUE)
ld <- as.matrix(read.table(list.files(dir, pattern="indexLD", full=TRUE)))
ld <- ld[info$idx, info$idx] # reorder according to indexinfo

# TODO

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

# write out files for testing ecaviar
nclust <- length(ld_mat)
ecav.bin <- "/proj/milovelab/love/bin/caviar/caviar/eCAVIAR"
ecav.coloc <- list()
for (j in seq_len(nclust)) {
  tmp <- sub("ecav", paste0("ecav_out_j",j), out.filename)
  ld.tmp <- sub("ecav_out","ecav_ld",tmp)
  z.tmp <- sub("ecav_out","ecav_z",tmp)
  z2.tmp <- sub("ecav_out","ecav_z2",tmp)
  write.table(ld_mat[[j]], file=ld.tmp,
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(sum_stat[[j]][,c("snp","z")], file=z.tmp,
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(sum_stat[[j]][,c("snp","z2")], file=z2.tmp,
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  system(paste(ecav.bin, "-o", tmp,
               "-l", ld.tmp, "-l", ld.tmp,
               "-z", z.tmp, "-z", z2.tmp, "-c 1"))
  ecav.coloc[[j]] <- read.table(paste0(tmp,"_col"), header=TRUE)
  file.remove(ld.tmp,z.tmp,z2.tmp) # remove ecaviar input
  tmp.out <- list.files(dirname(tmp),basename(tmp), full.names=TRUE)
  file.remove(tmp.out) # remove ecaviar output
}

save(ecav.coloc, file=out.filename)
