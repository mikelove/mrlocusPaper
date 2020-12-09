devtools::load_all("../../mrlocus")
genes <- list(Artery=c("MRAS","PHACTR1"),
              Blood="LIPC",
              Liver=c("CETP","LIPC","SORT1"))
trait <- c(MRAS="CAD",PHACTR1="CAD",CETP="HDL",LIPC="HDL",SORT1="LDL")

png(file="../supp/figs/realloci.png", width=1000, height=1000, res=120)
#png(file="../supp/figs/sort1.png", width=1200, height=1200, res=200)
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

for (tissue in names(genes)) {
  for (gene in genes[[tissue]]) {
    load(paste0(tissue,"-",gene,".rda"))
    out <- rstan::summary(res$stanfit, pars="alpha", probs=c(.1,.9))$summary[,c("mean","sd","10%","90%"),drop=FALSE]
    write.table(format(out, digits=4), file=paste0(tissue,"-",gene,".txt"), quote=FALSE, row.names=FALSE)
  }
}
