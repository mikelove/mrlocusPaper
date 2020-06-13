sim <- read.delim("out/1_1pctmodel_0.1h2g_0.001ve.scan.tsv", strings=FALSE)
idx <- which(sim$eqtl.true != 0)
with(sim, plot(eqtl.beta/eqtl.se, 
               col=ifelse(eqtl.true == 0, 1, 2),
               pch=ifelse(eqtl.true == 0, 1, 19)))
with(sim, plot(gwas.beta/gwas.se, 
               col=ifelse(gwas.true == 0, 1, 2),
               pch=ifelse(gwas.true == 0, 1, 19)))
with(sim, plot(eqtl.beta/eqtl.se, 
               col=ifelse(eqtl.true == 0, 1, 2),
               pch=ifelse(eqtl.true == 0, 1, 19)))


snps <- sim$snp[1:50]
write(snps, file="snp.list")
# /proj/milovelab/bin/plink/plink --bfile out/1_1pctmodel_0.1h2g_0.001ve --r2 --ld-snp-list snp.list --allow-no-sex --ld-window-r2 0 --ld-window 1000000 --out snp --memory 2048

ld <- read.table("snp.ld", header=TRUE)
ld <- ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]

ld.mat <- matrix(nrow=50, ncol=50, dimnames=list(snps, snps))
ld.mat[cbind(ld$SNP_A,ld$SNP_B)] <- sqrt(ld$R2)
round(ld.mat[1:10,1:10], 2)

library(Matrix)
Sigma_a <- ld.mat
Sigma_np <- as.matrix(nearPD(Sigma_a)$mat)

load_all("../../mrlocus")
options(mc.cores=2)
idx <- 1:50
fit1 <- fitBetaEcaviar(nsnp=50,
                       beta_hat_a=-1 * sim$eqtl.beta[idx],
                       beta_hat_b=-1 * sim$gwas.beta[idx],
                       se_a=sim$eqtl.se[idx],
                       se_b=sim$gwas.se[idx],
                       Sigma_a=Sigma_np,
                       Sigma_b=Sigma_np)

nsnp <- 50
rstan::stan_plot(fit1, pars=paste0("beta_a[",1:nsnp,"]"))
rstan::stan_plot(fit1, pars=paste0("beta_b[",1:nsnp,"]"))
print(fit1, pars=paste0("beta_a[",1:nsnp,"]"), digits=3)
coefs1 <- rstan::extract(fit1)
par(mfrow=c(1,2), mar=c(4.5,4,3,1))
col <- ifelse(sim$eqtl.true[idx] == 0, 1, 2)
pch <- ifelse(sim$eqtl.true[idx] == 0, 1, 19)
beta_clean_a <- colMeans(coefs1$beta_a)
plot(-1 * sim$eqtl.beta[idx], beta_clean_a, asp=1, main="eQTL", xlab="summary stat", ylab="mrlocus", col=col, pch=pch)
abline(0,1); abline(h=0, lty=2)
beta_clean_b <- colMeans(coefs1$beta_b)
plot(-1 * sim$gwas.beta[idx], beta_clean_b, asp=1, main="GWAS", xlab="summary stat", ylab="mrlocus", col=col, pch=pch)
abline(0,1); abline(h=0, lty=2)
