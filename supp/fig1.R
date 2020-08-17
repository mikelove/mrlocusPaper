devtools::load_all("~/proj/mrlocus")

png(file="~/Downloads/fig1_sim.png", width=2400, height=900, res=300)
par(mfrow=c(1,3))
load("~/Downloads/sim_lowdisp.rda")
plotMrlocus(res, main="SNPs → gene → trait\n(mediation with low dispersion)",
            label="Effect size of", legend=FALSE, pointers=TRUE, ylim=c(-2.5,2.5))
load("~/Downloads/sim_highdisp.rda")
plotMrlocus(res, main="SNPs → gene → trait\n(mediation with high dispersion)",
            label="Effect size of", legend=FALSE, ylim=c(-2.5,2.5))
load("~/Downloads/sim_alpha0_highdisp.rda")
plotMrlocus(res, main="SNPs → gene → trait\n(no mediation)",
            label="Effect size of", legend=FALSE, ylim=c(-2.5,2.5))
dev.off()
