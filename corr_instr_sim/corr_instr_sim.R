n <- 8
r2 <- .9
sd <- .1
corr <- sqrt(r2)
Sigma <- matrix(corr, nrow=n, ncol=n)
diag(Sigma) <- 1
beta_hat <- MASS::mvrnorm(2, mu=rep(0,n), Sigma)
rownames(beta_hat) <- c("a","b")
beta_hat["b",] <- sign(beta_hat["a",]) * beta_hat["b",]
beta_hat["a",] <- abs(beta_hat["a",])
sd_b <- sd_a <- rep(sd,n)
#devtools::load_all("../../mrlocus")
res <- list(beta_hat_a=beta_hat["a",], sd_a=sd_a, beta_hat_b=beta_hat["b",], sd_b=sd_b)
suppressWarnings({ fit <- fitSlope(res, verbose=FALSE, open_progress=FALSE, show_messages=FALSE, refresh=-1) })
rstan::summary(fit$stanfit, pars="alpha", probs=c(.1,.9))$summary
plotMrlocus(fit)
