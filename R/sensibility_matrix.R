
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
sensibility_matrix <- function(y, x, theta){
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
  sensibility <- matrix(0, ncol = number_of_parameters, nrow = number_of_parameters)
  eta <- 0.5 + numeric(length = K)
  for (needle_t in 2:n) {
    eta  <- omega  + sum(gamma * x[needle_t-1, ]) + sum(alpha * y[needle_t-1, 2:(K+1)]) + beta * eta
    for (needle_k in 1:K) {
      matrix_partial_eta[1:K,needle_k] <- (1:K == needle_k)*1 + beta * (1:K == needle_k) * as.vector(matrix_partial_eta[1:K,needle_k])
      matrix_partial_eta[(K+1):(K+P),needle_k] <- x[needle_t-1, ] + beta[needle_k] * as.vector(matrix_partial_eta[(K+1):(K+P),needle_k])
      matrix_partial_eta[(K+P+1):(2*K+P),needle_k] <- y[needle_t-1, 2:(K+1)] + beta[needle_k]  * as.vector(matrix_partial_eta[(K+P+1):(2*K+P),needle_k])
      matrix_partial_eta[(2*K+P+1):(3*K+P), needle_k] <- eta * (1:K == needle_k) + beta[needle_k] * as.vector(matrix_partial_eta[(2*K+P+1):(3*K+P),needle_k])
    }
    cumulate_eta <- cumsum(eta)

    sensibility_t <- matrix(0, ncol = number_of_parameters, nrow = number_of_parameters)


    if (needle_t > 10){

      #### calcul de h_un
      right_none <- numeric(number_of_parameters)
      for(needle_ell in 1:K){
        right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
      }
      h_un <-  (1/((1+sum(exp(cumulate_eta)))^2))* right_none  %*% t(matrix_partial_eta[,1])


      #### calcul de h_deux
      right_none <- numeric(number_of_parameters)
      for(needle_ell in 2:K){
        right_none <- right_none +   sum(exp(cumulate_eta)[needle_ell:K]) * matrix_partial_eta[,needle_ell]
      }
      h_deux <- (1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[2:K]) * matrix_partial_eta[,1]  + (1 + exp(cumulate_eta)[1]) * right_none) %*% t(matrix_partial_eta[,2])

      sensibility_t <- sensibility_t + h_un + h_deux
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
          sensibility_t <- sensibility_t + (1/((1+sum(exp(cumulate_eta)))^2)) * (sum(exp(cumulate_eta)[needle_k:K]) * left_none + (1 + sum(exp(cumulate_eta)[1:(needle_k -1)])) * right_none) %*% t(matrix_partial_eta[,needle_k])
        }
      }
    }
    sensibility <- sensibility + sensibility_t

  }
  sensibility/(n-10)
}
