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
##' @export

mis_rate <- function(a, mu0, mu1, sigma0, sigma1, p0, p1) {
    m0 <- t(a) %*% mu0; m1 <- t(a) %*% mu1
    s0 <- sqrt(t(a) %*% sigma0 %*% a); s1 <- sqrt(t(a) %*% sigma1 %*% a)
    dm <- m1 - m0
    ds2 <- s0 ^ 2 - s1 ^ 2
    delta <- sqrt((dm)^2 + ds2 * log(s0^2 / s1^2))
    temp <- p1 * pnorm((s1 * dm - s0 * delta) / ds2) -
        p1 * pnorm((s1 * dm + s0 * delta) / ds2) +
        p0 * pnorm((s0 * dm + s1 * delta) / ds2) -
        p0 * pnorm((s0 * dm - s1 * delta) / ds2)
    if (identical(s0, s1)) {
        s <- s0
        return(pnorm(-abs(dm) / (2 * s)))
    } else if (s0 > s1)
        return(p1 + temp)
    else
        return(p0 + temp)
}
