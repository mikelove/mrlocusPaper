i <- 1
  
files <- sub(".final","",list.files("out", pattern="final"))
files <- grep(paste0("^",i,"_"),files,value=TRUE)
files
final <- list()
for (k in seq_along(files)) {
  final[[k]] <- read.table(paste0("out/",files[k],".final"),skip=1)
}

mrlocus <- list()
for (k in seq_along(files)) {
  mrlocus[[k]] <- read.table(paste0("out/",files[k],".mrlocus"))
}
mrlocus10 <- sapply(mrlocus, function(x) x[2,1])
mrlocus90 <- sapply(mrlocus, function(x) x[2,2])

h2 <- as.numeric(sub(".*_(.*)h2_.*","\\1",files[1]))
ve <- as.numeric(sub(".*_(.*)ve$","\\1",files[1]))
ttl <- paste0("Simulation: ",100*h2,"% h2g, ",100*ve,"% var. exp.")
ttl

idx <- c(3,5,6:8)
dat <- data.frame(rep=rep(1:20,each=5),
                  true=rep(sapply(final, function(x) x[1,2]), each=5),
                  method=rep(c("causal","all","twmr","ptwas","mrlocus"),times=20),
                  estimate=as.vector(sapply(final, function(x) x$V2[idx])),
                  se=as.vector(sapply(final, function(x) x$V3[idx])))
dat$method <- factor(dat$method, c("causal","all","twmr","ptwas","mrlocus"))
dat$est_nozero <- ifelse(dat$est == 0, NA, dat$est)
q <- qnorm(.9)
dat$min <- dat$estimate - q * dat$se
dat$max <- dat$estimate + q * dat$se
dat$min[dat$method == "mrlocus"] <- mrlocus10
dat$max[dat$method == "mrlocus"] <- mrlocus90

if (FALSE) {
  # add additional methods
  est.lda <- unname(sapply(files, function(f)
    scan(paste0("ldamregger/",f,".ldamregger"),quiet=TRUE)[1]))
  est.pmr <- unname(sapply(files, function(f) {
    load(paste0("pmr/",f,".pmr"));
    if (is.null(res)) NA else res$causal_effect
  }))
  ests <- c(est.lda, est.pmr)
  ests[is.na(ests)] <- 0
  dat2 <- data.frame(rep=rep(1:20,2),
                     true=rep(dat$true[dat$method=="causal"],2),
                     method=rep(c("lda-mr-egger","pmr-sum-egger"),each=20),
                     estimate=ests, se=rep(1,40),
                     est_nozero=ests,
                     min=ests, max=ests)
  dat <- rbind(dat, dat2)
}

library(dplyr)
tab <- dat %>% group_by(method) %>%
  summarize(
    RMAE=mean(abs(est_nozero-true)/abs(true),na.rm=TRUE),
    MAE=mean(abs(est_nozero-true),na.rm=TRUE)
    )
tab
mx <- max(abs(dat$true))
lex <- 1.2 # limits expansion
data.tb <- tibble(x=-lex*mx, y=lex*mx, tb=list(tab))

library(ggplot2)
library(ggpmisc)
nl <-  nlevels(dat$method)
cols <- unname(palette.colors( nl+2 ))[-c(1,5)]
shps <- c(24,25,17,15,16, head(7:14,nl-5) )
png(file=paste0("../supp/figs/sim",i,".png"), res=150, width=800, height=800)
#png(file=paste0("../supp/figs/sim",i,"extra.png"), res=150, width=800, height=800)
ggplot(dat, aes(true,estimate,color=method,shape=method)) +
  geom_point(size=2) +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values=cols) +
  scale_shape_manual(values=shps) +
  geom_table(data=data.tb, aes(x, y, label=tb),
             table.theme = ttheme_gtlight,
             stat="fmt_tb") +
  xlim(-lex*mx,lex*mx) + ylim(-lex*mx,lex*mx) +
  ggtitle(ttl)
dev.off()

###

dat$contain <- dat$true > dat$min & dat$true < dat$max & !is.na(dat$est_nozero)

tab <- dat %>% group_by(method) %>%
  summarize(cov=paste0("cov: ",100*mean(contain),"%"))
tab
mx <- max(abs(dat$true))
tab$x <- "left"
tab$y <- "top"

#png(file=paste0("../supp/figs/cover",i,".png"), res=170, width=1200, height=800)
p2 <- ggplot(dat, aes(true,estimate,ymin=min,ymax=max,color=contain)) +
  geom_pointrange(shape="square", size=.5) + facet_wrap(~method) +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values=c(2,1)) +
  geom_text_npc(data=tab, aes(npcx=x, npcy=y, label=cov)) + 
  xlim(-1.2*mx,1.2*mx) + ylim(-1.2*mx,1.2*mx) +
  ggtitle(ttl)
#dev.off()

library(patchwork)
png(file="../supp/figs/fig2.png", res=170, width=2000, height=800)
p1 + p2 + plot_annotation(tag_levels = "A")
dev.off()

