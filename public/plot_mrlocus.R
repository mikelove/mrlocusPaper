devtools::load_all("../../mrlocus")
tissues <- list(SORT1="Liver",MRAS="Artery",PHACTR1="Artery",
                CETP="Liver",LIPC=c("Liver","Blood"))          
traits <- c(MRAS="CAD",PHACTR1="CAD",CETP="HDL",LIPC="HDL",SORT1="LDL")

xmax <- c(MRAS=.7, PHACTR1=.5, CETP=.3, LIPC=.8, SORT1=4)
ymax <- c(MRAS=.15, PHACTR1=.2, CETP=1, LIPC=.3, SORT1=.15)

gene <- "SORT1"; tissue <- "Liver"
method <- "mrlocus"
method <- "ecav-mrlocus"
png(file=paste0("../supp/figs/real_loci_",method,".png"), width=1800, height=1400, res=175)
#png(file="../supp/figs/sort1.png", width=1200, height=1200, res=200)
par(mfrow=c(2,3), mar=c(5,5,2,1))
for (gene in names(tissues)) {
  for (tissue in tissues[[gene]]) {
    trait <- traits[gene]
    dir <- paste(tissue, gene, trait, sep="_")
    load(file.path(dir, paste0(dir, ".", method)))
    main <- paste0("SNPs → ",gene," (",sub("_"," ",tissue),") → ",trait)
    plotMrlocus(res, main=main, xlim=c(0,xmax[gene]),
                ylim=c(-ymax[gene],ymax[gene]))
  }
}
dev.off()

# this chunk before final plot
genes <- list(Artery=c("MRAS","PHACTR1"),
              Blood="LIPC",
              Liver=c("CETP","LIPC","SORT1"))
traits <- c(MRAS="CAD",PHACTR1="CAD",CETP="HDL",LIPC="HDL",SORT1="LDL")
for (tissue in names(genes)) {
  for (gene in genes[[tissue]]) {
    trait <- traits[gene]
    dir <- paste(tissue, gene, trait, sep="_")
    for (method in c("mrlocus","ecav-mrlocus")) {
      load(file.path(dir, paste0(dir, ".", method)))
      if (!is.list(res)) {
        out <- data.frame(mean=0,sd=0,"10%"=0,"90%"=0)
      } else if ("stanfit" %in% names(res)) {
        out <- rstan::summary(
                        res$stanfit, pars="alpha",
                        probs=c(.1,.9)
                      )$summary[,c("mean","sd","10%","90%"),drop=FALSE]
      } else {
        out <- data.frame(mean=res$est[1], sd=res$est[2])
        out$"10%" <- out$mean + qnorm(.1) * out$sd
        out$"90%" <- out$mean + qnorm(.9) * out$sd
      }
      write.table(format(out, digits=4),
                  file=file.path(dir, paste0(dir, "_", method, ".txt")),
                  quote=FALSE, row.names=FALSE)
    }
  }
}

# dispersion
method <- "mrlocus"
dat <- data.frame(mean_beta_a=numeric(),
                  alpha=numeric(),
                  sigma=numeric())
for (gene in names(tissues)) {
  for (tissue in tissues[[gene]]) {
    trait <- traits[gene]
    dir <- paste(tissue, gene, trait, sep="_")
    load(file.path(dir, paste0(dir, ".", method)))
    mean_beta_a <- mean(res$beta_hat_a)
    alpha <- rstan::summary(res$stanfit, pars="alpha")$summary[1]
    sigma <- rstan::summary(res$stanfit, pars="sigma")$summary[1]
    dat <- rbind(dat, data.frame(mean_beta_a=mean_beta_a, alpha=alpha, sigma=sigma))
  }
}
rownames(dat) <- paste(names(unlist(tissues)), unlist(tissues))
dat$mean_mediated <- with(dat, alpha * mean_beta_a)
dat$sigma_over_mm <- with(dat, sigma / abs(mean_mediated))
write.csv(format(dat,digits=3), file="mrlocus_disp_table.txt", quote=FALSE)
