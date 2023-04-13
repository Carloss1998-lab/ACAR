paths <- function(n,P,n_cat,theta, eta_0 = 0.5){
  K <- n_cat-1
  x <- x_path(n,P)
  y <- matrix(data = 0, ncol = n_cat, nrow = n)
  y[1,1] <- 1
  ## vérifier la taille de theta
  omega <- theta[1:K]
  gamma <- theta[(K+1):(K+P)]
  alpha <- theta[(K+P+1):(2*K + P)]
  beta <- theta[(2*K + P+1):(3*K + P)]

  eta <- eta_0 + numeric(length = K)

  for(needle_t in 2:n){
    eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
    for (needle_k in 1:K) {
      numerator <- exp(cumsum(eta))
      pi <- numerator  /(1+sum( numerator) )
    }
    y[needle_t, ] <- rmultinom(1, 1, c(max(0,1-sum(pi)), pi))[,1]
  }
  return(list(y = y, x = x))
}
