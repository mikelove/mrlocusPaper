i <- 4

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
ttl <- paste0("Sim ",i,": ",100*h2,"% h2g, ",100*ve,"% var. exp.")
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

library(dplyr)
tab <- dat %>% group_by(method) %>%
  summarize(
    RMAE=mean(abs(est_nozero-true)/abs(true),na.rm=TRUE),
    MAE=mean(abs(est_nozero-true),na.rm=TRUE)
    )
tab
mx <- max(abs(dat$true))
data.tb <- tibble(x=-1.2*mx, y=1.2*mx, tb=list(tab))

library(ggplot2)
library(ggpmisc)
cols <- unname(palette.colors(7))[-c(1,5)]
ggplot(dat, aes(true,estimate,col=method)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values=cols) +
  geom_table(data=data.tb, aes(x, y, label=tb),
             table.theme = ttheme_gtlight,
             stat="fmt_tb") +
  xlim(-1.2*mx,1.2*mx) + ylim(-1.2*mx,1.2*mx) +
  ggtitle(ttl)

###

dat$contain <- dat$true > dat$min & dat$true < dat$max & !is.na(dat$est_nozero)

tab <- dat %>% group_by(method) %>%
  summarize(cov=paste0("cov: ",100*mean(contain),"%"))
tab
mx <- max(abs(dat$true))
tab$x <- "left"
tab$y <- "top"

ggplot(dat, aes(true,estimate,ymin=min,ymax=max,color=contain)) +
  geom_pointrange() + facet_wrap(~method) +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values=c(2,1)) +
  geom_text_npc(data=tab, aes(npcx=x, npcy=y, label=cov)) + 
  xlim(-1.2*mx,1.2*mx) + ylim(-1.2*mx,1.2*mx) +
  ggtitle(ttl)
