#' Bayesian Assurance Computation
#'
#' Takes in a set of parameters and returns the bayesian assurance
#' @param theta_0 Prior parameter provided by client
#' @param theta_1 Prior parameter under design and analysis stage
#' @param sig_sq Fixed variance
#' @param n sample size (either vector or numerical value)
#' @param n_a sample size at analysis stage
#' @param n_d sample size at design stage
#' @param alpha significance level
#' @return Bayesian assurance
#' @export
assurance_nd_na <- function(theta_0, theta_1, sig_sq, n, n_a, n_d, alpha){
  delta <- theta_1 - theta_0
  z_alpha <- qnorm(alpha)
  phi_val <- sqrt((n_d / (n + n_d)) * (1 + n_a / n)) * (((sqrt(n) * delta) / sqrt(sig_sq)) + z_alpha)
  b_assurance <- pnorm(phi_val)
  return(b_assurance)
}


