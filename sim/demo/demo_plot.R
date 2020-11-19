df$est <- factor(df$est, levels=c("summary","mrlocus","true"))
alt <- df$SNP[df$est == "true" & df$coef > 0 & df$data == "eQTL"]
df$color <- factor(ifelse(df$SNP == alt, 1, 2))
library(ggplot2)
ggplot(df, aes(x=SNP, y=coef, col=color)) +
  geom_point(size=2, show.legend=FALSE) +
  geom_errorbar(aes(ymin=coef-1.96*se,ymax=coef+1.96*se), show.legend=FALSE) +
  geom_hline(yintercept=0) +
  facet_grid(data ~ est, scales="free_y") +
  scale_color_manual(values=c("1"="blue","2"="black")) +
  scale_x_continuous(breaks=seq_len(max(df$SNP)))
