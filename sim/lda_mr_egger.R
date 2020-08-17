cmd_args=commandArgs(TRUE)

ss.filename <- cmd_args[1] # is scan.tsv
ld.filename <- cmd_args[2] # is LD matrix

sum.stat <- read.delim(ss.filename, strings=FALSE)
ld <- as.matrix(read.table(ld.filename))

# https://rbarfield.github.io/Barfield_website/pages/Rcode.html
LDA.MREgger<-function(X,Y,W){
  bX<-cbind(1,X)
  bread<-solve(crossprod(bX,W)%*%bX)
  theEsts<-bread%*%crossprod(bX,W%*%Y)
  theresid<-c(Y-theEsts[1]-X*theEsts[2])
  Sig.Est<-c(crossprod(theresid,W%*%theresid))/(length(X)-2)
  finresults<- cbind(theEsts,diag(bread)*Sig.Est)
  TestStat<-theEsts/sqrt(finresults[,2])
  Pvals<-2*pt(abs(TestStat),df = nrow(bX)-2,lower.tail = F)
  return(cbind(finresults,TestStat,Pvals))
}

sigma <- as.matrix(Matrix::nearPD(ld)$mat)
se_G <- matrix(sum.stat$gwas.se,ncol=1)

X <- sum.stat$eqtl.beta
Y <- sum.stat$gwas.beta
eps <- diag(rep(1e-12,length(se_G)))

W <- sigma %*% solve(tcrossprod(se_G) + eps)
res <- LDA.MREgger(X,Y,W)
write(unname(c(res[2,1], sqrt(res[2,2]))), file=cmd_args[3])
