i <- "1"

#extra_methods <- (i %in% c("1","high_n"))
extra_methods <- FALSE

files <- sub(".final","",list.files(paste0("out/",i), pattern="final"))
files <- grep(paste0("^",i,"_"),files,value=TRUE)
files
final <- list()
for (k in seq_along(files)) {
  final[[k]] <- read.table(paste0("out/",i,"/",files[k],".final"),skip=1)
}

mrlocus <- list()
for (k in seq_along(files)) {
  mrlocus[[k]] <- read.table(paste0("out/",i,"/",files[k],".mrlocus"))
}
mrlocus10 <- sapply(mrlocus, function(x) x[2,1])
mrlocus90 <- sapply(mrlocus, function(x) x[2,2])

if (extra_methods) {
  # alternative threshold
  mrlocus2 <- list()
  for (k in seq_along(files)) {
    mrl.file <- paste0("out/",i,"/",files[k],".mrlocus_p1e-4")
    if (file.info(mrl.file)$size > 0) {
      mrlocus2[[k]] <- read.table(mrl.file)
    } else {
      mrlocus2[[k]] <- matrix(c(0,0,0,0),ncol=2)
    }
  }
  mrlocus2.10 <- sapply(mrlocus2, function(x) x[2,1])
  mrlocus2.90 <- sapply(mrlocus2, function(x) x[2,2])
}

# ecaviar-mrlocus
ecav.mrlocus <- list()
for (k in seq_along(files)) {
  ecav.mrlocus[[k]] <- read.table(paste0("out/",i,"/",files[k],".ecav-mrlocus"))
}
ecav.mrlocus10 <- sapply(ecav.mrlocus, function(x) x[2,1])
ecav.mrlocus90 <- sapply(ecav.mrlocus, function(x) x[2,2])


h2 <- as.numeric(sub(".*_(.*)h2_.*","\\1",files[1]))
ve <- as.numeric(sub(".*_(.*)ve$","\\1",files[1]))
ttl <- paste0("Simulation: ",100*h2,"% h2g, ",100*ve,"% var. exp.")
ttl

idx <- c(3,5,6:nrow(final[[1]]))
if (extra_methods) {
  meths <- c("causal","all","twmr","twmr_p1e-4","ptwas",
             "ptwas_t0.1","mrlocus","mrlocus_p1e-4","ecaviar-mrlocus")
} else {
  meths <- c("causal","all","twmr","ptwas","mrlocus","ecaviar-mrlocus")
}
meths.big <- c("causal","all","twmr","twmr_p1e-4","ptwas",
               "ptwas_t0.1","mrlocus","mrlocus_p1e-4","ecaviar-mrlocus")

nsim <- length(files)
dat <- data.frame(rep=rep(1:nsim, each=length(idx)),
                  true=rep(sapply(final, function(x) x[1,2]), each=length(idx)),
                  method=rep(meths,times=nsim),
                  estimate=as.vector(sapply(final, function(x) x$V2[idx])),
                  se=as.vector(sapply(final, function(x) x$V3[idx])))
dat$method <- factor(dat$method, meths)
dat$est_nozero <- ifelse(dat$est == 0, NA, dat$est)
q <- qnorm(.9)
dat$min <- dat$estimate - q * dat$se
dat$max <- dat$estimate + q * dat$se
dat$min[dat$method == "mrlocus"] <- mrlocus10
dat$max[dat$method == "mrlocus"] <- mrlocus90
dat$min[dat$method == "ecaviar-mrlocus"] <- ecav.mrlocus10
dat$max[dat$method == "ecaviar-mrlocus"] <- ecav.mrlocus90
if (extra_methods) {
  dat$min[dat$method == "mrlocus_p1e-4"] <- mrlocus2.10
  dat$max[dat$method == "mrlocus_p1e-4"] <- mrlocus2.90
}

if (FALSE) {
  # add even more methods (one plot only)
  est.lda <- unname(sapply(files, function(f)
    scan(paste0("ldamregger/",f,".ldamregger"),quiet=TRUE)[1]))
  est.pmr <- unname(sapply(files, function(f) {
    load(paste0("pmr/",f,".pmr"));
    if (is.null(res)) NA else res$causal_effect
  }))
  ests <- c(est.lda, est.pmr)
  ests[is.na(ests)] <- 0
  dat2 <- data.frame(rep=rep(1:nsim,2), true=rep(dat$true[dat$method=="causal"],2),
                     method=rep(c("lda-mr-egger","pmr-sum-egger"),each=nsim),
                     estimate=ests, se=rep(1,2*nsim), est_nozero=ests,
                     min=ests, max=ests)
  dat <- rbind(dat, dat2)
  meths <- c(meths, c("lda-mr-egger", "pmr-sum-egger"))
}

# check number of instruments
num_instr <- unname(sapply(files, function(f) {
  length(scan(paste0("out/",i,"/",f,".mrl_keep"),quiet=TRUE))
}))
table(num_instr)
dat$two_plus_instr <- factor(rep( ifelse(num_instr > 1, "yes", "no"), each=length(meths) ))

library(dplyr)
tab <- dat %>% filter(two_plus_instr == "yes") %>% group_by(method) %>%
  summarize(
    RMAE=mean(abs(est_nozero-true)/abs(true),na.rm=TRUE),
    MAE=mean(abs(est_nozero-true),na.rm=TRUE)
    )
tab
mx <- max(abs(dat$true))
lex <- 1.2 # limits expansion
data.tb <- tibble(x=0, y=lex*mx, tb=list(tab))

dat2 <- dat %>% mutate(estimate = sign(true) * estimate, true = abs(true))
library(ggplot2)
library(ggpmisc)
cols <- unname(palette.colors())[-c(1,5)]
cols <- cols[c(1:3,3,4,4,5,5,5)]
names(cols) <- meths.big
shps <- c(24,25,17,18,15,7,16,13,10)
names(shps) <- meths.big
if (length(meths) == 11) { # for the plot with LDA and PMR
  cols <- unname(palette.colors())[-c(1,5)]
  cols <- cols[c(1:3,3,4,4,5,5,5,6,7)]
  names(cols) <- meths
  shps <- c(24,25,17,18,15,7,16,13,10,3,8)
  names(shps) <- meths
  mx <- max(abs(dat$estimate))
}
#png(file=paste0("../supp/figs/sim",i,".png"), res=150, width=800, height=800)
#png(file=paste0("../supp/figs/sim",i,"extra.png"), res=150, width=1200, height=800)
p1 <- ggplot(dat2, aes(true,estimate,color=method,shape=method)) +
  geom_point(size=2) +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values=cols) +
  scale_shape_manual(values=shps) +
  geom_table(data=data.tb, aes(x, y, label=tb),
             table.theme = ttheme_gtlight,
             stat="fmt_tb") +
#  xlim(0,lex*mx) + ylim(-.25*lex*mx,lex*mx) +
  xlim(0,lex*mx) + ylim(-.6,2.1) +
  ggtitle(ttl)
p1
#dev.off()

###

dat$contain <- dat$true > dat$min & dat$true < dat$max & !is.na(dat$est_nozero)

if (FALSE) {
  mrl_cover <- dat %>% filter(method=="mrlocus") %>% pull(contain)
  addmargins(table(mrl_cover, num_instr), 1)
}

tab <- dat %>% filter(two_plus_instr == "yes") %>% group_by(method) %>%
  summarize(cov=paste0("cov: ",100*round(mean(contain),2),"%"))
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
p2
#dev.off()

library(patchwork)
#png(file="../supp/figs/fig2.png", res=170, width=2000, height=800)
p1 + p2 + plot_annotation(tag_levels = "A")
#dev.off()

