k <- function(x, y) (x %*% y) ^ 2

##' This function executes the QDA by Projection method with kernel trick
##'
##' Place holder
##' @title k_qdap
##' @export
k_qdap <- function(x, y, xnew, lambda = 0, iter = 1, method = "Penalization") {
    x <- data.matrix(x)
    xnew <- data.matrix(xnew)
    x0 <- x[which(y == 0), ]
    x1 <- x[which(y == 1), ]
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    k0 <- k(x, t(x0))
    k1 <- k(x, t(x1))
    mu0 <- rowMeans(k0)
    mu1 <- rowMeans(k1)
    sigma0 <- k0 %*% (diag(n0) - outer(rep(1, n0), rep(1, n0)) / n0) %*% t(k0) /
        (n0 - 1) + 0.001 * diag(n0 + n1)
    sigma1 <- k1 %*% (diag(n1) - outer(rep(1, n1), rep(1, n1)) / n1) %*% t(k1) /
        (n1 - 1) + 0.001 * diag(n0 + n1)
    sigma <- ((n0 - 1) * sigma0 + (n1 - 1) * sigma1) /
        (n0 + n1 - 1)
    drt <- drt(mu0, mu1, sigma0, sigma1, sigma, iter, lambda, method)
    a <- drt$a
    ## print(a)
    conv <- drt$conv
    qda.fit <- qda(t(t(a) %*% k(x, t(x))), y)
    ynew <- predict(qda.fit, t(t(a) %*% k(x, t(xnew))))$class
    return(list(class = ynew, drt = a, conv = conv))
}
