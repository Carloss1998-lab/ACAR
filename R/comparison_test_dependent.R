#' Title
#'
#' @param model_poly_autoregressive1
#' @param model_poly_autoregressive2
#'
#' @return
#' @export
#'
#' @examples
comparison_test_dependent <- function(model_poly_autoregressive1, model_poly_autoregressive2){
  nsample <- model_poly_autoregressive1$nsample ## vérifier que les deux ont même taille
  hat_theta_1 <- model_poly_autoregressive1$parameters
  hat_theta_2 <- model_poly_autoregressive2$parameters
  inv_sensibility_1 <- solve(model_poly_autoregressive1$Sensibility)
  inv_sensibility_2 <- solve(model_poly_autoregressive2$Sensibility)
  score_1 <- sequence_score(model_poly_autoregressive1$y, model_poly_autoregressive1$x, model_poly_autoregressive1$parameters)
  score_2 <- sequence_score(model_poly_autoregressive2$y, model_poly_autoregressive2$x, model_poly_autoregressive2$parameters)
  S <- matrix(0, ncol = length(hat_theta_2), nrow = length(hat_theta_2))
  for (needle_t in 1:nrow(score_1)) {
    S <- S + score_1[needle_t,] %*% t(score_2[needle_t,])
  }
  S <- S/nrow(score_1)
  matrix_v <- (nsample-1)* model_poly_autoregressive1$parameters_covariance +  (nsample-1)* model_poly_autoregressive2$parameters_covariance  +  inv_sensibility_1 %*% S %*% t(inv_sensibility_2) + inv_sensibility_2 %*% t(S) %*% t(inv_sensibility_1)
  global_statistic <- (nsample-1) * (hat_theta_1 - hat_theta_2) %*% solve(matrix_v) %*% (hat_theta_1 - hat_theta_2)
  p_value <- pchisq(global_statistic, df = length(hat_theta_2), lower.tail = FALSE)
  list(global_statistic = global_statistic, p_value = p_value, test_variable =abs(sqrt(nsample-1) * (hat_theta_1 - hat_theta_2) / sqrt(diag(matrix_v))))
}


