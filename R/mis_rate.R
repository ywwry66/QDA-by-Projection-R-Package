##' This is the classification error/misclassification error function
##'
##' Place holder
##' @title mis_rate function
##' @param a direction
##' @param mu0 mean of class 0
##' @param mu1 mean of class 1
##' @param sigma0 covariance matrix of class 0
##' @param sigma1 covariance matrix of class 1
##' @return classification error
##' @author Ruiyang Wu


mis_rate <- function(a, mu0, mu1, sigma0, sigma1) {
    m0 <- t(a) %*% mu0; m1 <- t(a) %*% mu1
    s0 <- sqrt(t(a) %*% sigma0 %*% a); s1 <- sqrt(t(a) %*% sigma1 %*% a)
    delta <- sqrt((m0 - m1)^2 + (s0^2 - s1^2) * log(s0^2 / s1^2))
    if (s0 == s1) {
        s <- s0
        return(pnorm(-abs(m0 - m1) / (2 * s)))
    } else
        return((1 / 2) *
               (pnorm((s1 * (m1 - m0) - s0 * delta) / (s0^2 - s1^2)) -
                pnorm((s1 * (m1 - m0) + s0 * delta) / (s0^2 - s1^2)) +
                pnorm((s0 * (m1 - m0) + s1 * delta) / (s0^2 - s1^2)) -
                pnorm((s0 * (m1 - m0) - s1 * delta) / (s0^2 - s1^2)) +
                1))
}
