devtools::load_all("../../mrlocus")
genes <- list(Artery=c("MRAS","PHACTR1"),
              Blood="LIPC",
              Liver=c("CETP","LIPC","SORT1"))
traits <- c(MRAS="CAD",PHACTR1="CAD",CETP="HDL",LIPC="HDL",SORT1="LDL")

#png(file="../supp/figs/realloci.png", width=1000, height=1000, res=120)
#png(file="../supp/figs/sort1.png", width=1200, height=1200, res=200)
method <- "mrlocus"
method <- "ecav-mrlocus"
xmax <- c(MRAS=.7, PHACTR1=.5, CETP=.3, LIPC=.8, SORT1=3.5)
ymax <- c(MRAS=.15, PHACTR1=.2, CETP=1, LIPC=.3, SORT1=.15)
par(mfrow=c(2,3), mar=c(5,5,2,1))
for (tissue in names(genes)) {
  for (gene in genes[[tissue]]) {
    trait <- traits[gene]
    dir <- paste(tissue, gene, trait, sep="_")
    load(file.path(dir, paste0(dir, ".", method)))
    main <- paste0("SNPs → ",gene," (",sub("_"," ",tissue),") → ",trait)
    plotMrlocus(res, main=main, xlim=c(0,xmax[gene]), ylim=c(-ymax[gene],ymax[gene]))
  }
}
#dev.off()

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
      out <- rstan::summary(
                      res$stanfit, pars="alpha",
                      probs=c(.1,.9)
                    )$summary[,c("mean","sd","10%","90%"),drop=FALSE]
      write.table(format(out, digits=4),
                  file=file.path(dir, paste0(dir, "_", method, ".txt")),
                  quote=FALSE, row.names=FALSE)
    }
  }
}
