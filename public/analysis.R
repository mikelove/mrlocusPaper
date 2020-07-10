genes <- list(Artery_Tibial=c("MRAS","PHACTR1"),
              Blood="LIPC",
              Liver=c("CETP","LIPC","SORT1"))

for (tissue in names(genes)) {
  for (gene in genes[[tissue]]) {

    set.seed(1)
    print(paste("---",tissue,"-",gene,"---"))
    ptm <- proc.time() 
    
    dir <- file.path(tissue, gene)
    cond.files <- sub(".tsv","",list.files(dir, ".tsv"))
    nclust <- length(cond.files)
    ld_mat <- list()
    sum_stat <- list()
    for (j in 1:nclust) {
      filename <- file.path(dir,paste0(cond.files[j], ".ld"))
      ld_mat[[j]] <- as.matrix(read.delim(filename,header=FALSE))
      filename <- file.path(dir,paste0(cond.files[j], ".tsv"))
      sum_stat[[j]] <- read.delim(filename)
      sum_stat[[j]] <- sum_stat[[j]][order(sum_stat[[j]]$pos),]
    }

    sapply(sum_stat, nrow)

    library(pheatmap)
    library(gridExtra)
    load_all("../../mrlocus")
    out1 <- collapseHighCorSNPs(sum_stat, ld_mat, plot=FALSE)
    out2 <- flipAllelesAndGather(out1$sum_stat, out1$ld_mat,
                                 a="eQTL", b="GWAS",
                                 ref="Ref", eff="Effect",
                                 beta="beta", se="se",
                                 a2_plink="Major_plink",
                                 snp_id="SNP", sep="_",
                                 ab_last=TRUE, plot=FALSE)

    library(Matrix)
    Sigma_npd <- out2$Sigma
    for (j in seq_along(out2$Sigma)) {
      Sigma_npd[[j]] <- as.matrix(nearPD(out2$Sigma[[j]])$mat)
    }

    nsnp <- lengths(out2$beta_hat_a)

    # plotInitEstimates(out2)

    load_all("../../mrlocus")
    options(mc.cores=2)
    beta_hat_a <- list()
    beta_hat_b <- list()
    for (j in seq_along(nsnp)) {
      print(j)
      if (nsnp[j] > 1) {
        cap.out <- capture.output({ 
          fit <- fitBetaColoc(beta_hat_a=out2$beta_hat_a[[j]],
                              beta_hat_b=out2$beta_hat_b[[j]],
                              se_a=out2$se_a[[j]],
                              se_b=out2$se_b[[j]],
                              Sigma_a=Sigma_npd[[j]],
                              Sigma_b=Sigma_npd[[j]],
                              verbose=FALSE,
                              open_progress=FALSE,
                              show_messages=FALSE,
                              refresh=-1)
        })
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
                sd_b=out2$se_b)

    res <- extractForSlope(res, plot=FALSE)
    res <- fitSlope(res, iter=10000)

    print(proc.time() - ptm)
    
    pdf(file=paste0(tissue,"-",gene,".pdf"))
    plotMrlocus(res, main=paste(tissue,"-",gene))
    dev.off()

    out <- summary(res$stanfit, pars="alpha", probs=c(.1,.9))$summary[,c("mean","sd","10%","90%"),drop=FALSE]
    write.table(format(out, digits=4), file=paste0(tissue,"-",gene,".txt"), quote=FALSE, row.names=FALSE)

  }
}

sessionInfo()
