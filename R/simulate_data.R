#'
#' #' Title
#' #'
#' #' @param n
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' simulate_data <- function(n) {
#'   set.seed(123) ## pour avoir des resultats reproductibles
#'
#'   ## Simulation des données
#'   X1 <- rnorm(n)
#'   X2 <- rbinom(n, 1, 0.5)
#'   X3 <- rpois(n, 1)
#'   X4 <- runif(n)
#'   X5 <- rexp(n)
#'   X6 <- rnorm(n)
#'   X <- data.frame(X1, X2, X3, X4, X5, X6)
#'
#'   ## Calcul de Y
#'   Y_linear <- X1 + 2*X2 + 3*X3 + 4*X4 + 5*X5 + 6*X6
#'   breaks <- quantile(Y_linear, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))
#'   labels <- c("A", "B", "C", "D", "E")
#'
#'   labels <- c("Very Unsatisfied", "Unsatisfied", "Neutral", "Satisfied", "Very Satisfied")
#'   Y <- cut(Y_linear, breaks, labels = labels, include.lowest = TRUE)
#'
#'   ## Création d'un data.frame avec les variables X et Y
#'   data <- data.frame(X, Y)
#'
#'   return(data)
#' }
#'
#' acar_df_sim = simulate_data(100)
