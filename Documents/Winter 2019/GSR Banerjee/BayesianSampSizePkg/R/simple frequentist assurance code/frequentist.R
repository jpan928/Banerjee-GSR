#' Frequentist Power Computation
#'
#' Takes in a set of parameters and returns the power
#' @param n sample size (either vector or numerical value)
#' @param sig_sq Fixed variance
#' @param theta_0 Prior parameter provided by client
#' @param theta_1 Prior parameter under design and analysis stage
#' @param alpha significance level
#' @return Frequentist Power
#' @export
pwr_freq <- function(n, sig_sq, alpha, theta_0, theta_1){
  delta <- theta_1 - theta_0
  z_alpha <- qnorm(alpha)
  pwr_val <- ((sqrt(n) * delta) / sqrt(sig_sq)) + z_alpha
  pwr <- pnorm(pwr_val)
  return(pwr)
}

















