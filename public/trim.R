cmd_args=commandArgs(TRUE)

dir <- cmd_args[1]
tsv.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
out.filename <- cmd_args[4]

set.seed(1)

source("common.R") # common function for mrlocus and ecaviar-mrlocus
out <- getTrimmedSumStats(dir, tsv.filename,ld.filename)
sum_stat <- out$sum_stat
ld_mat <- out$ld_mat
trim_clusters <- out$trim_clusters

if (length(trim_clusters) > 0) {
  write(trim_clusters, file=out.filename, ncolumns=length(trim_clusters))
} else {
  write(NA, file=out.filename)
}
