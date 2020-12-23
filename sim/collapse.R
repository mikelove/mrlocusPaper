devtools::load_all("../../mrlocus")
noparen <- function(z) sub("\\(1\\)","",z)
sims <- c(1:9,paste0("null",1:3),"high_n")
for (s in sims) {
  print(paste("---",s,"---"))
  dir <- file.path("out",s)
  files <- list.files(dir, pattern="ve\\.clumped$")
  for (f in files) {
    file <- sub(".clumped","",f)
    clumped.filename <- paste0(dir,"/",file,".clumped")
    scan.filename <- paste0(dir,"/",file,".scan.tsv")
    mrl.coloc.filename <- paste0(dir,"/",file,".mrl_coloc")
    clumped <- read.table(clumped.filename, strings=FALSE, header=TRUE)
    big_sum_stat <- read.delim(scan.filename, strings=FALSE)
    clumps <- lapply(strsplit(clumped$SP2,split=","), noparen)
    # add index SNP to clumps
    clumps <- lapply(seq_along(clumps), function(j) c(clumped$SNP[j], clumps[[j]]))
    sum_stat <- lapply(clumps, function(x) {
      z <- big_sum_stat[big_sum_stat$snp %in% x,]
      z$abs.z <- abs(z$eqtl.beta/z$eqtl.se)
      z
    })
    pre <- sapply(sum_stat, nrow)
    load(mrl.coloc.filename)
    post <- sapply(res$alleles, nrow)
    write(pre, file="collapse_pre.txt", ncolumns=length(sum_stat), append=TRUE)
    write(post, file="collapse_post.txt", ncolumns=length(sum_stat), append=TRUE)
  }
}
