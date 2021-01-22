library(ggplot2)
library(dplyr)

genes <- list(Artery=c("MRAS","PHACTR1"),
              Blood="LIPC",
              Liver=c("CETP","LIPC","SORT1"))
traits <- c(MRAS="CAD",PHACTR1="CAD",CETP="HDL",LIPC="HDL",SORT1="LDL")
dat <- list()
for (tissue in names(genes)) {
  for (gene in genes[[tissue]]) {
    trait <- unname( traits[gene] )
    dir <- paste(tissue, gene, trait, sep="_")
    twmr <- read.table(file.path(dir,paste0(dir,".twmr")), header=TRUE)
    mrlocus <- read.table(file.path(dir,paste0(dir,"_mrlocus.txt")), header=TRUE)
    ecav.mrlocus <- read.table(file.path(dir,paste0(dir,"_ecav-mrlocus.txt")), header=TRUE)
    dat[[dir]] <- data.frame(tissue=tissue,gene=gene,trait=trait,
                             method=c("twmr","mrlocus","ecaviar-mrlocus"),
                             alpha=c(twmr$alpha,mrlocus$mean,ecav.mrlocus$mean),
                             se=c(twmr$SE,mrlocus$sd,ecav.mrlocus$sd),
                             xmin=c(NA,mrlocus$X10.,ecav.mrlocus$X10.),
                             xmax=c(NA,mrlocus$X90.,ecav.mrlocus$X90.))
  }
}
dat <- do.call(rbind, unname(dat))

# add PTWAS
ptwas <- read.table("ptwas.txt", header=TRUE)
ptwas$method <- "ptwas"
ptwas$xmin <- NA
ptwas$xmax <- NA

dat <- rbind(dat, ptwas[,c("tissue","gene","trait","method","alpha","se","xmin","xmax")])

q <- qnorm(.9)
dat <- dat %>% mutate(xmin=ifelse(is.na(xmin), alpha-q*se, xmin),
                      xmax=ifelse(is.na(xmax), alpha+q*se, xmax))

dat <- dat %>% mutate(label=paste0(gene," (",tissue,")", " → ", trait))
lvls <- unique(dat$label)[c(6,1,2,4,5,3)]
dat$label <- factor(dat$label, levels=lvls)

png(file="../supp/figs/forest.png", width=1200, height=600, res=150)
ggplot(dat, aes(alpha,method,xmin=xmin,xmax=xmax)) +
  geom_pointrange(size=.3) +
  facet_wrap(~label, scales="free") +
  geom_vline(xintercept = 0, lty=2) +
  xlab("gene-to-trait estimate")
dev.off()
