
#' Title
#'
#' @param n
#' @param P
#' @param n_cat
#' @param theta
#' @param theta_other
#' @param eta_0
#'
#' @return
#' @export
#'
#' @examples
bipaths_concorde <- function(n,P,n_cat,theta, theta_other,  eta_0 = 0.5){
  K <- n_cat-1
  x1 <- x_path(n,P)
  x2 <-  x_path(n,P)
  y1 <- matrix(data = 0, ncol = n_cat, nrow = n)
  y1[1,1] <- 1
  y2 <- matrix(data = 0, ncol = n_cat, nrow = n)
  y2[1,1] <- 1
  ## vérifier la taille de theta
  omega <- theta[1:K]
  gamma <- theta[(K+1):(K+P)]
  alpha <- theta[(K+P+1):(2*K + P)]
  beta <- theta[(2*K + P+1):(3*K + P)]

  omega_other <- theta_other[1:K]
  gamma_other <- theta_other[(K+1):(K+P)]
  alpha_other <- theta_other[(K+P+1):(2*K + P)]
  beta_other <- theta_other[(2*K + P+1):(3*K + P)]

  eta1 <- eta_0 + numeric(length = K)
  eta2 <- eta_0 + numeric(length = K)

  for(needle_t in 2:n){
    eta1  <- omega  + sum(gamma * x1[needle_t-1, ]) + sum(alpha * y1[needle_t-1, 2:(K+1)]) + beta * eta1
    eta2  <- omega_other  + sum(gamma_other * x2[needle_t-1, ]) + sum(alpha_other * y2[needle_t-1, 2:(K+1)]) + beta_other * eta2

    for (needle_k in 1:K) {
      numerator1 <- exp(cumsum(eta1))
      pi1 <- numerator1  /(1+sum( numerator1) )
      numerator2 <- exp(cumsum(eta2))
      pi2 <- numerator2  /(1+sum( numerator2) )
    }
    u <- runif(1)
    y1[needle_t, ] <- (which.max(u-cumsum(c(1-sum(pi1), pi1))<0)==1:n_cat) * 1
    y2[needle_t, ] <- (which.max(u-cumsum(c(1-sum(pi2), pi2))<0)==1:n_cat) * 1
  }
  return(list(y1 = y1, x1 = x1, y2 = y2, x2 = x2))
}

