files <- sub(".final","",list.files("out", pattern="final"))
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

library(ggplot2)
library(ggpmisc)
cols <- unname(palette.colors(6))[-5]
ggplot(dat, aes(true,estimate,col=method)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values=cols) +
  xlim(-.75,.75) + ylim(-.75,.75) +
  geom_table(data=data.tb, aes(x, y, label=tb),
             table.theme = ttheme_gtlight,
             stat="fmt_tb")
