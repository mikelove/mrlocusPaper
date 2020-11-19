devtools::load_all("../../mrlocus")
library(pbapply)

niter <- 100
out <- pbsapply(1:niter, function(i) {
  n <- 6
  r2 <- 0.9
  sdev <- .5
  corr <- sqrt(r2)
  Sigma <- matrix(corr, nrow=n, ncol=n)
  diag(Sigma) <- 1
  beta_hat_a <- MASS::mvrnorm(1, mu=rep(5,n), Sigma)  
  beta_hat_b <- MASS::mvrnorm(1, mu=rep(0,n), Sigma)
  sd_b <- sd_a <- rep(sdev,n)
  res <- list(beta_hat_a=beta_hat_a, sd_a=sd_a, beta_hat_b=beta_hat_b, sd_b=sd_b)
  co <- capture.output({
    suppressWarnings({ 
      fit <- fitSlope(res, verbose=FALSE, open_progress=FALSE, show_messages=FALSE, refresh=-1)
    })
  })
  bnds <- rstan::summary(fit$stanfit, pars="alpha", probs=c(.1,.9))$summary[1,c("10%","90%")]
  unname(bnds)
})
mean(apply(out, 2, function(x) sign(x[1]) == sign(x[2])))
