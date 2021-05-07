##' This function executes the QDA by Projection method
##'
##' Place holder
##' @title qdap
##' @param x training predictors
##' @param y training labels
##' @param xnew test predictors
##' @param lambda 
##' @param iter number of iterations
##' @param method 
##' @return a list of the predicted labels, optimal direction and convergence status
##' @author Ruiyang Wu
##' @export

qdap <- function(x, y, xnew, lambda = 0, iter = 1,
                 method = "Penalization", optim = "BFGS") {
    x <- data.matrix(x)
    xnew <- data.matrix(xnew)
    x0 <- x[which(y == 0), ]
    x1 <- x[which(y == 1), ]
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    p0 <- n0 / (n0 + n1)
    p1 <- n1 / (n0 + n1)
    p <- ncol(x)
    mu0 <- colMeans(x0)
    mu1 <- colMeans(x1)
    sigma0 <- cov(x0)
    sigma1 <- cov(x1)
    if (!is_invertible(sigma0))
        sigma0 <- sigma0 + diag(10e-7, p)
    if (!is_invertible(sigma1))
        sigma1 <- sigma1 + diag(10e-7, p)
    sigma <- ((n0 - 1) * sigma0 + (n1 - 1) * sigma1) /
        (n0 + n1 - 2)
    drt <- drt(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
               iter, lambda, method, optim)
    a <- drt$a
    ## print(a)
    conv <- drt$conv
    ## 1D qda
    x <- as.vector(x %*% a)
    x0 <- x[which(y == 0)]
    x1 <- x[which(y == 1)]
    mu0 <- mean(x0)
    mu1 <- mean(x1)
    sigma0 <- var(x0)
    sigma1 <- var(x1)
    c2 <- 1 / sigma0 - 1 / sigma1
    c1 <- -2 * (mu0 / sigma0 - mu1 / sigma1)
    c0 <- mu0 ^ 2 / sigma0 - mu1 ^ 2 / sigma1 + log(sigma0 / sigma1) -
        2 * log(p0 / p1)
    qda_rule <- function(xnew) as.integer(c2 * xnew ^ 2 + c1 * xnew + c0 > 0)
    ynew <- sapply(as.vector(xnew %*% a), FUN = qda_rule)
    return(list(class = ynew, drt = a, conv = conv))
}

## mylda <- function(x, y, xnew) {
##     x <- data.matrix(x)
##     xnew <- data.matrix(xnew)
##     x0 <- as.matrix(x[which(y == 0), ])
##     x1 <- as.matrix(x[which(y == 1), ])
##     p <- ncol(x)
##     n0 <- nrow(x0)
##     n1 <- nrow(x1)
##     mu0 <- colMeans(x0)
##     mu1 <- colMeans(x1)
##     sigma0 <- cov(x0)
##     sigma1 <- cov(x1)
##     sigma <- ((n0 - 1) * sigma0 + (n1 - 1) * sigma1) /
##         (n0 + n1 - 2)
##     a0 <- solve(sigma + diag(0.0000001, p)) %*% (mu1 - mu0)
##     if (p == 1) a0 <- solve(sigma) %*% (mu1 - mu0)
##     pred <- function(xnew) xnew %*% a0 > (mu0 + mu1) %*% a0 / 2 - log(n1 / n0)
##     ynew <- apply(xnew, MARGIN = 1, FUN = pred)
##     return(list(class = ynew))
## }
