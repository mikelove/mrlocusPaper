cmd_args=commandArgs(TRUE)

clumped.filename <- cmd_args[1]
scan.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
keep.filename <- cmd_args[4]
ecav.filename <- cmd_args[5]
out.filename <- cmd_args[6]

if (FALSE) {
  i <- "1"
  dir <- file.path("out",i)
  files <- list.files(dir, pattern="ve.clumped")
  file <- sub(".clumped","",files[5])
  clumped.filename <- paste0(dir,"/",file,".clumped")
  scan.filename <- paste0(dir,"/",file,".scan.tsv")
  ld.filename <- paste0(dir,"/",file,".ld")
  keep.filename <- paste0(dir,"/",file,".mrl_keep")
  ecav.filename <- paste0(dir,"/",file,".ecav")
}

clumped <- read.table(clumped.filename, strings=FALSE, header=TRUE)
big_sum_stat <- read.delim(scan.filename, strings=FALSE)
big_ld_mat <- as.matrix(read.table(ld.filename))

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

# read in the clusters to keep
keep.clusters <- scan(keep.filename)

# trim clusters
sum_stat <- sum_stat[keep.clusters]

# ecaviar
load(ecav.filename)
ecav.coloc <- ecav.coloc[keep.clusters]

nclusters <- length(sum_stat)
beta_hat_a <- numeric(nclusters)
beta_hat_b <- numeric(nclusters)
sd_a <- numeric(nclusters)
sd_b <- numeric(nclusters)
snp <- character(nclusters)
for (j in seq_len(nclusters)) {
  idx <- ecav.coloc[[j]]$SNP_ID[ which.max(ecav.coloc[[j]]$CLPP) ]
  row <- sum_stat[[j]][ sum_stat[[j]]$snp == idx,]
  beta_hat_a[j] <- abs(row$eqtl.beta)
  beta_hat_b[j] <- sign(row$eqtl.beta) * row$gwas.beta
  sd_a[j] <- row$eqtl.se
  sd_b[j] <- row$gwas.se
  snp[j] <- row$snp
}
res <- list(beta_hat_a=beta_hat_a, beta_hat_b=beta_hat_b, sd_a=sd_a, sd_b=sd_b, snp=snp)

# remove clusters that were below p threshold
z.thr <- qnorm(.001/2, lower.tail=FALSE)
abs.z <- with(res, abs(beta_hat_a/sd_a))
z.idx <- abs.z > z.thr
if (!any(z.idx)) {
  # such a sim won't contribute to sim eval,
  # bc it doesn't have allelic heterogeneity.
  # just run slope fitting with one SNP to produce output
  z.idx <- abs.z == max(abs.z) 
}
res <- lapply(res, `[`, z.idx)

# construct LD matrix of eCAVIAR candidate SNPs
ld.idx <- match(res$snp, big_sum_stat$snp)
r2 <- big_ld_mat[ld.idx, ld.idx]^2

# second round of trimming using eCAVIAR candidate SNPs
devtools::load_all("../../mrlocus")

nclusters <- length(sum_stat)
if (nclusters > 1) {
  # find clusters that have pairwise r2 with other clumps above a threshold
  trim_clusters <- trimClusters(r2, r2_threshold=0.05)
  
  if (length(trim_clusters) > 0) {
    res <- lapply(res, function(x) x[-trim_clusters])
  }  
} else {
  trim_clusters <- numeric()
}

# write out the second round clusters that are kept
keep_clusters <- 1:nclusters
if (length(trim_clusters) > 0) {
  keep_clusters <- keep_clusters[-trim_clusters]
}
write(keep_clusters, file=sub("ecav-mrlocus","ecav-mrl_keep2",out.filename), ncolumns=length(keep_clusters))

# mrlocus slop fitting
res <- fitSlope(res, iter=10000)

#rstan::stan_plot(res$stanfit, pars=c("alpha","sigma"))
if ("stanfit" %in% res) {
  print(res$stanfit, pars=c("alpha","sigma"), digits=3)
  summary(res$stanfit, pars="alpha", probs=c(.1,.9))$summary
}

if ("stanfit" %in% names(res)) {
  s <- summary(res$stanfit, pars="alpha", probs=c(.025,.05,.1,.9,.95,.975))$summary
  mrlocus.out <- matrix(s[,c("mean","sd","10%","90%","5%","95%","2.5%","97.5%")],
                        ncol=2, byrow=TRUE)
} else {
  mrlocus.out <- matrix(c(res$est,
                          res$est[1] + qnorm(.1) * res$est[2],
                          res$est[1] + qnorm(.9) * res$est[2],
                          res$est[1] + qnorm(.05) * res$est[2],
                          res$est[1] + qnorm(.95) * res$est[2],
                          res$est[1] + qnorm(.025) * res$est[2],
                          res$est[1] + qnorm(.975) * res$est[2]),
                          ncol=2, byrow=TRUE)
}

write.table(mrlocus.out, file=out.filename, row.names=FALSE, col.names=FALSE)
