files <- lapply(c("null1","null2","null3"), function(i) {
  sub(".final","",list.files(paste0("out/",i), pattern="final"))
  })
files <- unlist(files)

final <- list()
dir <- substr(files, 1, 5)
for (k in seq_along(files)) {
  final[[k]] <- read.table(paste0("out/",dir[k],"/",files[k],".final"),skip=1)
}

mrlocus <- list()
for (k in seq_along(files)) {
  mrlocus[[k]] <- read.table(paste0("out/",dir[k],"/",files[k],".mrlocus"))
}
mrlocus10 <- sapply(mrlocus, function(x) x[2,1])
mrlocus90 <- sapply(mrlocus, function(x) x[2,2])

# ecaviar-mrlocus
ecav.mrlocus <- list()
for (k in seq_along(files)) {
  ecav.mrlocus[[k]] <- read.table(paste0("out/",dir[k],"/",files[k],".ecav-mrlocus"))
}
ecav.mrlocus10 <- sapply(ecav.mrlocus, function(x) x[2,1])
ecav.mrlocus90 <- sapply(ecav.mrlocus, function(x) x[2,2])

h2 <- 100*as.numeric(sub(".*_(.*)h2_.*","\\1",files))

idx <- c(3,5,6:nrow(final[[1]]))
meths <- c("causal","all","twmr","ptwas","mrlocus","ecaviar-mrlocus")
nsim <- unname(table(dir))

dat <- data.frame(rep=unlist(lapply(1:3, function(i) rep(1:nsim[i],each=length(idx)))),
                  h2=rep(h2,each=length(idx)),
                  true=0,
                  method=rep(meths,sum(nsim)),
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
dat$h2 <- factor(dat$h2)
levels(dat$h2) <- paste0("h2: ",levels(dat$h2),"%")

# check number of instruments
num_instr <- unname(sapply(seq_along(files), function(k) {
  length(scan(paste0("out/",dir[k],"/",files[k],".mrl_keep"),quiet=TRUE))
}))
sapply(split(num_instr, dir), table)
dat$two_plus_instr <- factor(rep( ifelse(num_instr > 1, "yes", "no"), each=length(meths) ))

dat$contain <- dat$true > dat$min & dat$true < dat$max & !is.na(dat$est_nozero)

library(dplyr)
tab <- dat %>% filter(two_plus_instr=="yes") %>%
  group_by(h2, method) %>%
  summarize(cov=paste0(100*round(mean(contain),2),"%"))
tab
tab$x <- "left"
tab$y <- "top"

# also show MAE
tab2 <- dat %>% filter(two_plus_instr=="yes") %>%
  group_by(h2, method) %>%
  summarize(MAE=round(mean(abs(estimate)),3))
tab2
tab2$x <- "left"
tab2$y <- "middle"

dat2 <- dat %>% filter(two_plus_instr=="yes")

ord10 <- dat2 %>% filter(h2 == "h2: 10%" & method == "causal") %>%
  arrange(estimate) %>% pull(rep)
ord20 <- dat2 %>% filter(h2 == "h2: 20%" & method == "causal") %>%
  arrange(estimate) %>% pull(rep)
ord5 <- dat2 %>% filter(h2 == "h2: 5%" & method == "causal") %>%
  arrange(estimate) %>% pull(rep)

dat2$replicate <- NA
dat2$replicate[dat2$h2 == "h2: 10%"] <- match(dat2$rep[dat2$h2 == "h2: 10%"], ord10)
dat2$replicate[dat2$h2 == "h2: 20%"] <- match(dat2$rep[dat2$h2 == "h2: 20%"], ord20)
dat2$replicate[dat2$h2 == "h2: 5%"] <- match(dat2$rep[dat2$h2 == "h2: 5%"], ord5)

library(ggplot2)
library(ggpmisc)
png(file="../supp/figs/nullplot.png", res=125, width=1200, height=800)
ggplot(dat2, aes(estimate,replicate,xmin=min,xmax=max,color=contain)) +
  geom_pointrange(shape="square",size=.35) +
  facet_grid(h2 ~ method, scales="free_y") +
  coord_cartesian(xlim=c(-1,1)) + 
  geom_vline(xintercept=0) +
  scale_color_manual(values=c(2,1)) +
  geom_text_npc(data=tab, aes(npcx=x, npcy=y, label=cov)) +
  geom_text_npc(data=tab2, aes(npcx=x, npcy=y, label=MAE))
dev.off()
