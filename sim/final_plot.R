files <- sub(".final","",list.files("out", pattern="final"))
files <- grep("^3_",files,value=TRUE)
files
final <- list()
for (i in seq_along(files)) {
  final[[i]] <- read.table(paste0("out/",files[i],".final"),skip=1)
}

dat <- data.frame(rep=rep(1:20,each=5),
                  true=rep(sapply(final, function(x) x[1,2]), each=5),
                  method=rep(c("causal","all","twmr","ptwas","mrlocus"),times=20),
                  estimate=as.vector(sapply(final, function(x) x$V2[c(3,5,6:8)])))
dat$method <- factor(dat$method, c("causal","all","twmr","ptwas","mrlocus"))
dat$est_nozero <- ifelse(dat$est == 0, NA, dat$est)

library(dplyr)
tab <- dat %>% group_by(method) %>%
  summarize(MAE=mean(abs(est_nozero-true),na.rm=TRUE),
            RMSE=sqrt(mean((est_nozero-true)^2,na.rm=TRUE)))
tab
data.tb <- tibble(x=-.7, y=.7, tb=list(tab))
data.tb <- tibble(x=-.9, y=.9, tb=list(tab))
data.tb <- tibble(x=-.4, y=.4, tb=list(tab))

library(ggplot2)
library(ggpmisc)
cols <- unname(palette.colors(6))[-5]
ggplot(dat, aes(true,estimate,col=method)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values=cols) +
  geom_table(data=data.tb, aes(x, y, label=tb),
             table.theme = ttheme_gtlight,
             stat="fmt_tb") +
#  xlim(-.75,.75) + ylim(-.75,.75) + 
#  ggtitle("Sim 1: 10% h2g, 1% var. exp.")
#  xlim(-1,1) + ylim(-1,1) + 
#  ggtitle("Sim 2: 10% h2g, 2% var. exp.")
  xlim(-.5,.5) + ylim(-.5,.5) + 
  ggtitle("Sim 3: 20% h2g, 1% var. exp.")
