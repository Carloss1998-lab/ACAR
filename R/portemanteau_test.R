
#' Title
#'
#' @param y
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
portemanteau_test <- function(y, x, theta){
  K = ncol(y)-1
  n = nrow(y)
  P = ncol(x)
  number_of_parameters <- length(theta)
  omega <- theta[1:K]
  gamma <- theta[(K+1):(K+P)]
  alpha <- theta[(K+P+1):(2*K + P)]
  beta <- theta[(2*K + P+1):(3*K + P)]

  matrix_partial_eta <- matrix(0,nrow = number_of_parameters, ncol = K)
  eta <- 0.5 + numeric(length = K)
  Sensibility <-  sensibility_matrix(y, x,  theta = theta)
  Information <- information_matrix(y, x,  theta = theta)
  c_10 <- numeric(number_of_parameters)
  c_20 <- numeric(number_of_parameters)
  c_3plus0 <- matrix(0,nrow = number_of_parameters, ncol = K-2)
  overline_mu2 <- matrix(0, nrow = K, ncol = K)
  G <- matrix(0,ncol = K, nrow  = number_of_parameters)
  rho <- numeric(K)
  err_moins1 <- rev(cumsum(rev(exp(cumsum(rep(1/2, K))))))/(1+sum(exp(cumsum(rep(1/2, K)))))
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
    err <- inverse_cum_y - rev(cumsum(rev(exp(cumulate_eta))))/(1+sum(exp(cumulate_eta)))

    if(needle_t>10){

      for(needle_k in 1:K){
        G[,needle_k] <- G[,needle_k] + err_moins1[needle_k] * err[needle_k] * apply(err* t(matrix_partial_eta) %*% t(solve(Sensibility)), 2, sum)
      }

      #### calcul de h_un
      right_none <- numeric(number_of_parameters)
      for(needle_ell in 1:K){
        right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
      }
      c_10 <- c_10 -  err[1]*(1/((1+sum(exp(cumulate_eta)))^2))* right_none


      #### calcul de h_deux
      right_none <- numeric(number_of_parameters)
      for(needle_ell in 2:K){
        right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
      }
      c_20 <- c_20 -  err[2]*(1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[2:K]) * matrix_partial_eta[,1]  + (1 + exp(cumulate_eta)[1]) * right_none)

      if(K>2){
        for (needle_k in 3:K) {
          left_none <- matrix_partial_eta[,1]
          for(needle_ell in 2:(needle_k-1)){
            left_none <- left_none + (1+sum(exp(cumulate_eta)[1:(needle_ell-1)])) * matrix_partial_eta[,needle_ell]
          }
          right_none <- numeric(number_of_parameters)
          for(needle_ell in needle_k:K){
            right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
          }
          c_3plus0[, needle_k-2] <- c_3plus0[, needle_k-2] - err[needle_k]*(1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[needle_k:K]) * left_none + (1 + sum(exp(cumulate_eta)[1:(needle_k -1)])) * right_none)
        }
      }

      rho_cur <- err_moins1 * err
      rho <- rho + rho_cur
      overline_mu2 <- overline_mu2 + rho_cur %*% t(rho_cur)
    }

    err_moins1 <- err
  }
  matrix_C <- cbind(c_10, c_20, c_3plus0)/(n-10)
  G <- G/(n-10)
  rho <- rho/(n-10)
  overline_mu2 <- overline_mu2/(n-10)

  matrix_W <- overline_mu2 + t(matrix_C) %*% solve(Sensibility)%*%Information%*%t(solve(Sensibility)) %*% matrix_C + t(matrix_C) %*% G + t(G) %*% matrix_C
  statistic <- (n-10) * rho %*% solve(matrix_W) %*% rho
  p_value <- pchisq(statistic, df = K, lower.tail = FALSE)
  list(statistic = statistic, p_value = p_value)
}
