cmd_args=commandArgs(TRUE)

Ngwas <- 100000
N_eQTLs <- 1000
out <- c("gene","alpha","SE","P","Nsnps","Ngene")

# arg 1 is clumped
# arg 2 is scan.tsv
# arg 3 is LD matrix
# arg 4 is output

clumped <- read.table(cmd_args[1], header=TRUE)[,-12]
sum_stat <- read.table(cmd_args[2], header=TRUE)
ld_mat <- as.matrix(read.table(cmd_args[3]))

idx <- sum_stat$snp %in% clumped$SNP

# beta is a matrix of gene effect sizes
beta <- as.matrix(sum_stat$eqtl.beta[idx], ncol=1)

# gamma is a matrix of GWAS effect sizes
gamma <- as.matrix(sum_stat$gwas.beta[idx], ncol=1)

# C is the LD matrix
C <- ld_mat[idx,idx]

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

write.table(out,file=cmd_args[4],quote=F,col.names=F,row.names=F)
