#' Frequentist and Bayesian Curve
#'
#' Takes in a set of parameters and returns the power
#' @param n sample size (either vector or numerical value)
#' @param y1 vector of Frequentist powers computed from a vector of n using pwr_freq function
#' @param y2 vector of Bayesian assurance computed from a vector of n using assurance_nd_na function
#' @return power curves outtputed using ggplot2
#' @export
pwr_curve <- function(n, y1, y2, freq, bayes){
  library(ggplot2)

  df <- as.data.frame(cbind(n, y1))
  df2 <- as.data.frame(cbind(n, y2))

  if(freq == TRUE & bayes == FALSE){
    p <- ggplot(df, aes(x = n, y = y1, color="Frequentist")) + geom_line(lwd=1.2)
    pwr_curves <- p + xlab("Sample Size (n)") + ylab("Power") +
      ggtitle("Frequentist Power Curve")
  }else if(freq == FALSE & bayes == TRUE){
    p <- ggplot(df2, aes(x = n, y = y2, color = "Bayesian")) + geom_line(lwd=1.2)
    pwr_curves <- p + xlab("Sample Size (n)") + ylab("Assurance") +
      ggtitle("Bayesian Assurance Curve")
  }else{
    p <- ggplot(df, aes(x = n, y = y1, color="Frequentist")) + geom_line(lwd=1.2)
    p1 <- p + geom_line(data = df2, aes(x = n, y = y2, color="Bayesian"),lwd=1.2)
    p2 <- p1 + xlab("Sample Size (n)") + ylab("Power/Assurance") +
      ggtitle("Power Curves of Frequentist and Bayesian Methods")
    # pwr_curves <- p2 + scale_color_discrete(name = "Power Curves",
    #                                         labels = c("Bayesian", "Frequentist"))
    pwr_curves <- p2 + theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"))
  }

  return(pwr_curves)
}


