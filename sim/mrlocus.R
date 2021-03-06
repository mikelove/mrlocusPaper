cmd_args=commandArgs(TRUE)

clumped.filename <- cmd_args[1]
scan.filename <- cmd_args[2]
ld.filename <- cmd_args[3]
out.filename <- cmd_args[4]

if (FALSE) {
  i <- "high_n"
  dir <- file.path("out",i)
  files <- list.files(dir, pattern="ve.clumped")
  file <- sub(".clumped","",files[5])
  clumped.filename <- paste0(dir,"/",file,".clumped")
  scan.filename <- paste0(dir,"/",file,".scan.tsv")
  ld.filename <- paste0(dir,"/",file,".ld")
  out.filename <- paste0(dir,"/",file,".mrlocus")
}

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

# construct LD matrix of index eSNPs (largest abs Z-score)
ld.idx <- sapply(sum_stat, function(x) x$snp[which.max(x$abs.z)])
ld.idx <- match(ld.idx, big_sum_stat$snp)
r2 <- big_ld_mat[ld.idx, ld.idx]^2

# trim clusters using index eSNP
devtools::load_all("../../mrlocus")

nclusters <- length(sum_stat)
if (nclusters > 1) {
  # find clusters that have pairwise r2 with other clumps above a threshold
  trim_clusters <- trimClusters(r2, r2_threshold=0.05)

  # write out pairwise r2 of index eSNPs that remain
  if (length(trim_clusters) > 0) {
    r2.out <- r2[-trim_clusters,-trim_clusters]
  } else {
    r2.out <- r2
  }
  r2.lower <- r2.out[lower.tri(r2.out)]
  # if trim_clusters is n-1 then there is no pairwise r2 to write out...
  if (length(trim_clusters) < nclusters-1) {
    write(r2.lower, file=sub("mrlocus","mrl_r2",out.filename), ncolumns=length(r2.lower))
  }  
} else {
  trim_clusters <- numeric(0)
}

# write out the clusters to keep
keep_clusters <- 1:nclusters
if (length(trim_clusters) > 0) {
  keep_clusters <- keep_clusters[-trim_clusters]
}
write(keep_clusters, file=sub("mrlocus","mrl_keep",out.filename), ncolumns=length(keep_clusters))

# trim clusters
if (length(trim_clusters) > 0) {
  sum_stat <- sum_stat[-trim_clusters]
  ld_mat <- ld_mat[-trim_clusters]
}

# mrlocus pre-processing steps

sapply(sum_stat, function(x) any(x$eqtl.true != 0))
out1 <- collapseHighCorSNPs(sum_stat, ld_mat, score="abs.z", plot=FALSE, snp_id="snp")
# naive look to see if true eQTLs are still present
sapply(out1$sum_stat, function(x) any(x$eqtl.true != 0))
# taking into account collapsed SNPs, the true eQTL should survive
sapply(seq_along(sum_stat), function(j) {
  eqtl.true <- sum_stat[[j]]$eqtl.true != 0
  if (any(eqtl.true)) {
    eqtl.nms <- sum_stat[[j]]$snp[eqtl.true]
    any(sapply(eqtl.nms, function(p) grepl(p, out1$sum_stat[[j]]$collapsed)))
  } else {
    FALSE
  }
})

out2 <- flipAllelesAndGather(out1$sum_stat, out1$ld_mat,
                             a="eqtl", b="gwas",
                             ref="a0", eff="a1",
                             beta="beta", se="se",
                             a2_plink="a0",
                             snp_id="snp", sep=".",
                             ab_last=FALSE, alleles_same=TRUE,
                             plot=FALSE)

#plotInitEstimates(out2)

library(Matrix)
Sigma <- lapply(out2$Sigma, function(x) as.matrix(nearPD(x)$mat))

# save runtime by rerunning slope fitting from after colocalization
coloc.filename <- sub("mrlocus","mrl_coloc",out.filename)
if (file.exists(coloc.filename)) {

  load(coloc.filename)
  
} else {

  # mrlocus colocalization
  
  options(mc.cores=4)
  nsnp <- lengths(out2$beta_hat_a)
  beta_hat_a <- list()
  beta_hat_b <- list()
  for (j in seq_along(nsnp)) {
    print(j)
    if (nsnp[j] > 1) {
      fit <- fitBetaColoc(beta_hat_a=out2$beta_hat_a[[j]],
                          beta_hat_b=out2$beta_hat_b[[j]],
                          se_a=out2$se_a[[j]],
                          se_b=out2$se_b[[j]],
                          Sigma_a=Sigma[[j]],
                          Sigma_b=Sigma[[j]],
                          verbose=FALSE,
                          open_progress=FALSE,
                          show_messages=FALSE,
                          refresh=-1)
      ## rstan::stan_plot(fit$stan,
      ##                  pars=c(paste0("beta_a[",1:nsnp[j],"]"),
      ##                         paste0("beta_b[",1:nsnp[j],"]")))
      beta_hat_a[[j]] <- fit$beta_hat_a
      beta_hat_b[[j]] <- fit$beta_hat_b
    } else {
      beta_hat_a[[j]] <- out2$beta_hat_a[[j]]
      beta_hat_b[[j]] <- out2$beta_hat_b[[j]]
    }
  }
  # make a results list for slope fitting
  res <- list(beta_hat_a=beta_hat_a,
              beta_hat_b=beta_hat_b,
              sd_a=out2$se_a,
              sd_b=out2$se_b,
              alleles=out2$alleles)
  # save the colocalization output
  save(res, file=sub("mrlocus","mrl_coloc",out.filename))

}

res <- extractForSlope(res, plot=FALSE)

nclusters <- nrow(res$alleles)
if (nclusters > 1) {
  # second round of trimming
  ld.idx <- match(res$alleles$id, big_sum_stat$snp)
  r2 <- big_ld_mat[ld.idx, ld.idx]^2
  if (nclusters > 1) {
    trim_clusters <- trimClusters(r2, r2_threshold=0.05)
    if (length(trim_clusters) > 0) {
      res <- lapply(res, function(x) {
        if (is.numeric(x)) {
          x[-trim_clusters]
        } else {
          x[-trim_clusters,,drop=FALSE]
        }
      })
    }
  }
} else {
  trim_clusters <- numeric()
}

# write out the second round clusters that are kept
keep_clusters <- 1:nclusters
if (length(trim_clusters) > 0) {
  keep_clusters <- keep_clusters[-trim_clusters]
}
write(keep_clusters, file=sub("mrlocus","mrl_keep2",out.filename), ncolumns=length(keep_clusters))

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
