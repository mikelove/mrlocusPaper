devtools::load_all("../../mrlocus")
noparen <- function(z) sub("\\(1\\)","",z)

sims <- c(1:9,"null")
for (s in sims) {
  print(paste("---",s,"---"))
  dir <- file.path("out",s)
  files <- list.files(dir, pattern="clumped")
  for (f in files) {
    print(sub("(.*)_500nq.*","\\1",f))
    file <- sub(".clumped","",f)
    clumped.filename <- paste0(dir,"/",file,".clumped")
    scan.filename <- paste0(dir,"/",file,".scan.tsv")
    ld.filename <- paste0(dir,"/",file,".ld")

    clumped <- read.table(clumped.filename, strings=FALSE, header=TRUE)
    big_sum_stat <- read.delim(scan.filename, strings=FALSE)
    big_ld_mat <- as.matrix(read.table(ld.filename))
    
    clumps <- lapply(strsplit(clumped$SP2,split=","), noparen)
    # add index SNP to clumps
    clumps <- lapply(seq_along(clumps), function(j) c(clumped$SNP[j], clumps[[j]]))

    ld_mat <- lapply(clumps, function(x) {
      idx <- big_sum_stat$snp %in% x
      unname(big_ld_mat[idx,idx,drop=FALSE])
    })
    sapply(ld_mat, nrow)
    sum_stat <- lapply(clumps, function(x) {
      z <- big_sum_stat[big_sum_stat$snp %in% x,]
      z$abs.z <- abs(z$eqtl.beta/z$eqtl.se)
      z
    })
    
    out1 <- collapseHighCorSNPs(sum_stat, ld_mat, score="abs.z", plot=FALSE)
    
    pre <- sapply(sum_stat, nrow)
    post <- sapply(out1$sum_stat, nrow)

    # filter out too large clumps
    pre <- pre[post < 100]
    post <- post[post < 100]
    
    write(pre, file="collapse_pre.txt", ncolumns=length(sum_stat), append=TRUE)
    write(post, file="collapse_post.txt", ncolumns=length(sum_stat), append=TRUE)
  }
}
