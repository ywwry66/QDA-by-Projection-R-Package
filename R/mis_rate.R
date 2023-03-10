##' Classification Error for Heteroscedastic Gaussian Models after 1-D
##' Projection
##'
##' This function computes the population classification
##' error/misclassification rate for heteroscedastic Gaussian models
##' after 1-D projection, given means and covariances of the model and
##' the direction of projection.
##'
##' This function only handles two-class heteroscedastic Gaussian models.
##'
##' @param a A vector representing the direction of projection.
##' @param mu0 A vector representing the mean of class 0
##' @param mu1 A vector representing the mean of class 1
##' @param sigma0 An intertible matrix representing the covariance
##'   matrix of class 0
##' @param sigma1 An invertible matrix representing the covariance
##'   matrix of class 1
##' @return The return value is the classification
##'   error/misclassification rate of the two-class heteroscedastic
##'   Gaussian model, whose class means are 'mu0' and 'mu1' and class
##'   covariance matrices are 'sigma0' and 'sigma1', after 1-D
##'   projection onto the direction 'a'.
##' @export

mis_rate <- function(a, mu0, mu1, sigma0, sigma1, p0, p1) {
  m0 <- t(a) %*% mu0
  m1 <- t(a) %*% mu1
  s0 <- sqrt(t(a) %*% sigma0 %*% a)
  s1 <- sqrt(t(a) %*% sigma1 %*% a)
  dm <- m1 - m0
  ds2 <- s0^2 - s1^2
  delta <- sqrt((dm)^2 + ds2 * log(s0^2 / s1^2))
  temp <- p1 * pnorm((s1 * dm - s0 * delta) / ds2) -
    p1 * pnorm((s1 * dm + s0 * delta) / ds2) +
    p0 * pnorm((s0 * dm + s1 * delta) / ds2) -
    p0 * pnorm((s0 * dm - s1 * delta) / ds2)
  if (identical(s0, s1)) {
    s <- s0
    return(pnorm(-abs(dm) / (2 * s)))
  } else if (s0 > s1) {
    return(p1 + temp)
  } else {
    return(p0 + temp)
  }
}
