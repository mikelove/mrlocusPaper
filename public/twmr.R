cmd_args=commandArgs(TRUE)

out <- c("gene","alpha","SE","P","Nsnps","Ngene")

# arg 1 is directory
# arg 2 is TSV files
# arg 3 is output

dir <- cmd_args[1]
tsv_files <- scan(cmd_args[2], what="char")
info <- read.table(list.files(dir, pattern="indexinfo", full=TRUE), header=TRUE)
ld <- as.matrix(read.table(list.files(dir, pattern="indexLD", full=TRUE)))
ld <- ld[info$idx, info$idx] # reorder according to indexinfo

# check that the TSV files in same order as 'c_i' in indexinfo file
stopifnot(all(sub(".*_(.*).eQTLBase.intersect.tsv","\\1",tsv_files) == as.character(info$c_i)))

# get the tissue and trait
tissue <- sub("_.*","",dir)
trait <- sub(".*_","",dir)

# hard-code sample sizes
CAD <- c(60801, 123504) # cases / controls
eqtl_samp_size <- c(Artery=663, Blood=31684, Liver=588)
gwas_samp_size <- c(CAD=4/(1/CAD[1] + 1/CAD[2]), HDL=315133, LDL=343621)

# set from hard-coded values
N_eQTLs <- eqtl_samp_size[tissue]
Ngwas <- gwas_samp_size[trait]

sum_stat <- t(sapply(seq_along(tsv_files), function(i) {
  x <- read.table(file.path(dir,tsv_files[i]),header=TRUE)
  row <- x[x$SNP == info$idxSNP[i],]
  swap <- row[1,"Ref_GWAS"] != row[1,"Ref_eQTL"]
  beta_eqtl <- row[1,"beta_eQTL"]
  beta_gwas <- row[1,"beta_GWAS"]
  if (swap) {
    beta_gwas <- -1 * beta_gwas
  }
  c(beta_eqtl=beta_eqtl, beta_gwas=beta_gwas)
}))

sum_stat <- as.data.frame(sum_stat)

# beta is a matrix of gene effect sizes
beta <- as.matrix(sum_stat$beta_eqtl, ncol=1)

# gamma is a matrix of GWAS effect sizes
gamma <- as.matrix(sum_stat$beta_gwas, ncol=1)

# C is the LD matrix
C <- ld

S<-t(beta)%*%solve(C)%*%beta
H<-(1-1/sqrt(3781))*S+(1/sqrt(3781))*diag(length(S[,1]))
alpha<-solve(H)%*%(t(beta)%*%solve(C)%*%gamma)

alpha<-as.vector(alpha)

C_inv <- solve(C)
GCG_inv <- t(beta) %*% solve(C) %*% beta
GCG_inv<-(1-1/sqrt(3781))*GCG_inv+(1/sqrt(3781))*diag(length(GCG_inv[,1]))
GCG_inv<-solve(GCG_inv)

   
df_dg <- GCG_inv %*% t(beta) %*% C_inv
df_dG <- (GCG_inv %x% (t(gamma) %*% C_inv %*% ((beta %*% GCG_inv %*% t(beta)) %*% C_inv + diag(nrow(beta))))) + ((-t(gamma) %*% C_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% C_inv))
J <- cbind(df_dG, df_dg)

SEs<-c(rep(1/sqrt(N_eQTLs),length(beta[1,])*length(beta[,1])),rep(1/sqrt(Ngwas),length(gamma[,1])))
R<-diag(length(beta[1,])+1)
Sigma <- (SEs %*% t(SEs)) * (C %x% R)   
V <- J %*% Sigma %*% t(J)
se<- sqrt(V[1,1])

N=length(beta[,1])
Ngene=length(beta[1,])
Z<-alpha[1]/se
pval<-2*pnorm(abs(Z),lower.tail=FALSE)
line<-c("gene",alpha[1],se,pval,N,Ngene)
out<-rbind(out,line)

write.table(out,file=cmd_args[3],quote=F,col.names=F,row.names=F)
