genes <- list(Artery_Tibial=c("MRAS","PHACTR1"),
              Liver=c("CETP","LIPC","SORT1"))

library(pheatmap)
library(gridExtra)

for (tissue in names(genes)) {
  two.ancestry <- grepl("BBJ", tissue)
  for (gene in genes[[tissue]]) {

    set.seed(1)
    print(paste("---",tissue,"-",gene,"---"))
    ptm <- proc.time() 
    
    dir <- file.path(tissue, gene)
    cond.files <- sub(".tsv","",list.files(dir, ".tsv"))
    nclust <- length(cond.files)
    ld_mat <- list()
    ld_mat2 <- list()
    sum_stat <- list()
    for (j in 1:nclust) {
      if (two.ancestry) {
        filename <- file.path(dir,paste0(cond.files[j], "_EUR.ld"))
        ld_mat[[j]] <- as.matrix(read.delim(filename,header=FALSE))
        gwas_filename <- file.path(dir,paste0(cond.files[j], "_EAS.ld"))
        ld_mat2[[j]] <- as.matrix(read.delim(gwas_filename,header=FALSE))
      } else {
        filename <- file.path(dir,paste0(cond.files[j], ".ld"))
        ld_mat[[j]] <- as.matrix(read.delim(filename,header=FALSE))
      }
      filename <- file.path(dir,paste0(cond.files[j], ".tsv"))
      sum_stat[[j]] <- read.delim(filename)
      sum_stat[[j]] <- sum_stat[[j]][order(sum_stat[[j]]$pos),]
    }

    sapply(sum_stat, nrow)

    devtools::load_all("../../mrlocus")
    out1 <- collapseHighCorSNPs(sum_stat, ld_mat, plot=FALSE)
    #out1 <- collapseHighCorSNPs(sum_stat, ld_mat, ld_mat2, plot=FALSE)
    if (two.ancestry) {
      a2_plink <- "Major_plink_EUR"
      a2_plink_mat2 <- "Major_plink_EAS"
    } else {
      a2_plink <- "Major_plink"
      a2_plink_mat2 <- NULL
    }
    out2 <- flipAllelesAndGather(out1$sum_stat, out1$ld_mat,
                                 out1$ld_mat2,
                                 a="eQTL", b="GWAS",
                                 ref="Ref", eff="Effect",
                                 beta="beta", se="se",
                                 a2_plink=a2_plink,
                                 a2_plink_mat2=a2_plink_mat2,
                                 snp_id="SNP", sep="_",
                                 ab_last=TRUE, plot=FALSE)
    
    library(Matrix)
    Sigma_npd <- out2$Sigma
    for (j in seq_along(out2$Sigma)) {
      Sigma_npd[[j]] <- as.matrix(nearPD(out2$Sigma[[j]])$mat)
    }
    if (two.ancestry) {
      Sigma_npd2 <- out2$Sigma2
      for (j in seq_along(out2$Sigma2)) {
        Sigma_npd2[[j]] <- as.matrix(nearPD(out2$Sigma2[[j]])$mat)
      }
    } else {
      Sigma_npd2 <- Sigma_npd
    }

    nsnp <- lengths(out2$beta_hat_a)

    plotInitEstimates(out2)

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
                              Sigma_b=Sigma_npd2[[j]],
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
    
    save(res, file=paste0(tissue,"-",gene,".rda"))
  }
}

sessionInfo()

if (FALSE) {
  devtools::load_all("../../mrlocus")
  genes <- list(Artery_Tibial=c("MRAS","PHACTR1"),
                Liver=c("CETP","LIPC"))
  trait <- c(MRAS="CAD",PHACTR1="CAD",CETP="HDL",LIPC="HDL",SORT1="LDL")

  png(file="../supp/figs/realloci.png", width=1000, height=1000, res=120)
  #png(file="../supp/figs/sort1.png", width=800, height=800, res=150)
  par(mfrow=c(2,2), mar=c(5,5,2,1))
  for (tissue in names(genes)) {
    for (gene in genes[[tissue]]) {
      load(paste0(tissue,"-",gene,".rda"))
      main <- paste0("SNPs → ",gene," (",sub("_"," ",tissue),") → ",trait[gene])
      if (gene == "CETP") {
        plotMrlocus(res, main=main, ylim=c(-1,1))
      } else if (gene == "LIPC") {
        plotMrlocus(res, main=main, ylim=c(-.25,.25))
      } else {
        plotMrlocus(res, main=main)
      }
    }
  }
  dev.off()

  genes <- list(Artery_Tibial=c("MRAS","PHACTR1"),
                Liver=c("CETP","LIPC","SORT1"))
  trait <- c(MRAS="CAD",PHACTR1="CAD",CETP="HDL",LIPC="HDL",SORT1="LDL")
  for (tissue in names(genes)) {
    for (gene in genes[[tissue]]) {
      load(paste0(tissue,"-",gene,".rda"))
      out <- rstan::summary(res$stanfit, pars="alpha", probs=c(.1,.9))$summary[,c("mean","sd","10%","90%"),drop=FALSE]
      write.table(format(out, digits=4), file=paste0(tissue,"-",gene,".txt"), quote=FALSE, row.names=FALSE)
    }
  }
  
}
