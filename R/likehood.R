#debut de code
likehood <- function(y, x, theta){
  # print("needle")
  # print("yes")
  K = ncol(y)-1
  n = nrow(y)
  P = ncol(x)
  omega <- theta[1:K]
  gamma <- theta[(K+1):(K+P)]
  alpha <- theta[(K+P+1):(2*K + P)]
  beta <- theta[(2*K + P+1):(3*K + P)]
  likehood_value <- 0
  eta <- 0.5 + numeric(length = K)
  for(needle_t in 2:n){
    eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
    cum_eta <- cumsum(eta)
    likehood_value <- likehood_value  - sum(y[needle_t, 2:(K+1)]*cum_eta)+log(1+sum(exp(cum_eta)))
  }
  # print(likehood_value)
  # print("fin")
  return(likehood_value)
}
