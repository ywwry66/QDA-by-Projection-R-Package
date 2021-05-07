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
    qda.fit <- qda(x %*% a, y)
    ynew <- predict(qda.fit, xnew %*% a)$class
    return(list(class = ynew, drt = a, conv = conv))
}
