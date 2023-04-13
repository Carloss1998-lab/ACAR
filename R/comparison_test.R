
#' Title
#'
#' @param model_poly_autoregressive1
#' @param model_poly_autoregressive2
#'
#' @return
#' @export
#'
#' @examples
comparison_test <- function(model_poly_autoregressive1, model_poly_autoregressive2){
  nsample <- model_poly_autoregressive1$nsample ## vérifier que les deux ont même taille
  hat_theta_1 <- model_poly_autoregressive1$parameters
  hat_theta_2 <- model_poly_autoregressive2$parameters
  matrix_v <- model_poly_autoregressive1$nsample * model_poly_autoregressive1$parameters_covariance +  model_poly_autoregressive2$nsample * model_poly_autoregressive2$parameters_covariance
  global_statistic <- (nsample-1) * (hat_theta_1 - hat_theta_2) %*% solve(matrix_v) %*% (hat_theta_1 - hat_theta_2)
  p_value <- pchisq(global_statistic, df = length(hat_theta_2), lower.tail = FALSE)
  list(global_statistic = global_statistic, p_value = p_value)
}
