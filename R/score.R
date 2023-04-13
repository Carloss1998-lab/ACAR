
score <- function(y, x, theta){
  K = ncol(y)-1
  n = nrow(y)
  P = ncol(x)
  number_of_parameters <- length(theta)
  omega <- theta[1:K]
  gamma <- theta[(K+1):(K+P)]
  alpha <- theta[(K+P+1):(2*K + P)]
  beta <- theta[(2*K + P+1):(3*K + P)]
  #partial_eta_initial <-  numeric(length = number_of_parameters)
  matrix_partial_eta <- matrix(0,nrow = number_of_parameters, ncol = K)
  eta <- 0.5 + numeric(length = K)
  score_value <- numeric(number_of_parameters)
  for (needle_t in 2:n) {
    eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
    for (needle_k in 1:K) {
      matrix_partial_eta[1:K,needle_k] <- (1:K == needle_k)*1 + beta * (1:K == needle_k) * as.vector(matrix_partial_eta[1:K,needle_k])
      matrix_partial_eta[(K+1):(K+P),needle_k] <- x[needle_t-1, ] + beta[needle_k] * as.vector(matrix_partial_eta[(K+1):(K+P),needle_k])
      matrix_partial_eta[(K+P+1):(2*K+P),needle_k] <- y[needle_t-1, 2:(K+1)] + beta[needle_k]  * as.vector(matrix_partial_eta[(K+P+1):(2*K+P),needle_k])
      matrix_partial_eta[(2*K+P+1):(3*K+P), needle_k] <- eta * (1:K == needle_k) + beta[needle_k] * as.vector(matrix_partial_eta[(2*K+P+1):(3*K+P),needle_k])
    }
    inverse_cum_y <- rev(cumsum(rev(y[needle_t, 2:(K+1)])))
    cumulate_eta <- cumsum(eta)
    err <- inverse_cum_y - rev(cumsum(rev(cumulate_eta)))
    score_value <- score_value + apply(-err * t(matrix_partial_eta), 2, sum)
  }
  score_value
}
