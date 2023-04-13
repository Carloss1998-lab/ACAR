#' Title
#'
#' @param y
#' @param x
#'
#' @return
#' @export
#'
#' @examples


poly_autoregression <- function(y, x,max_iteration= 10000 ,repeat_opt = 20){
  likehood_data <- function(theta) likehood(y=as.matrix(y), x = as.matrix(x), theta = theta)
  K = ncol(y)-1
  P = ncol(x)

  list_values <- numeric(repeat_opt)*NA
  list_code <- numeric(repeat_opt)*NA
  mat_param <- matrix(NA, ncol = 3*K+P , nrow = repeat_opt)
  for(needle_rep in 1:repeat_opt){
    parauto <- rexp(2*K)
    theta_init <- c(rnorm(K+P), parauto/(2*sum(parauto)))
    res_opt <- tryCatch({
      optim(theta_init, fn = likehood_data, control = list(maxit = max_iteration), method = "L-BFGS-B", 
            lower = c(rep(-Inf, 2*K+P), (-1+1e-5) * rep(1, K)),
            upper = c(rep(Inf, 2*K+P), (1-1e-5) * rep(1, K)))
    }, error = function(e) {
      return(NA)
    })
    
    print("optim function finished running")
    if (length(res_opt) > 1) {
      list_values[needle_rep] <- res_opt$value
      list_code[needle_rep] <- res_opt$convergence
      mat_param[needle_rep,] <- res_opt$par
    } else {
      print("optim function encountered an error")
    }
    if(length(res_opt)>1){
      list_values[needle_rep] <- res_opt$value
      mat_param[needle_rep, ] <- res_opt$par
    }
  }

  parameters <- mat_param[which.min(list_values),]
  loss_value <- min(list_values, na.rm = TRUE)
  Information <-  information_matrix(y, x, parameters)
  Sensibility <-  sensibility_matrix(y, x, parameters)
  parameters_covariance <- solve(Sensibility)%*%Information%*%t(solve(Sensibility))/(nrow(y)-1)
  standard_deviation <- sqrt(diag(parameters_covariance))
  AIC <- 2*loss_value + 2 * (3 * K + P)
  portemanteau = portemanteau_test(y = y, x = x, theta = parameters)
  list(parameters = parameters, standard_deviation = standard_deviation, AIC =  AIC,
       loss_value =  loss_value, Information = Information, Sensibility = Sensibility,
       parameters_covariance = parameters_covariance, nsample = nrow(y),
       number_of_parameters = (3 * K + P), portemanteau = portemanteau, y=y, x = x, convergence = list_code[which.min(loss_value)]) #
}


