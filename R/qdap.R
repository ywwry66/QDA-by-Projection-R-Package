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

qdap <- function(x, y, xnew, lambda = 0, iter = 1, method = "Penalization") {
    x <- data.matrix(x)
    xnew <- data.matrix(xnew)
    x0 <- x[which(y == 0), ]
    x1 <- x[which(y == 1), ]
    p <- ncol(x)
    mu0 <- colMeans(x0)
    mu1 <- colMeans(x1)
    sigma0 <- cov(x0)
    sigma1 <- cov(x1)
    sigma <- ((nrow(x0) - 1) * sigma0 + (nrow(x1) - 1) * sigma1) /
        (nrow(x0) + nrow(x1) - 1)
    drt <- drt(mu0, mu1, sigma0, sigma1, sigma, iter, lambda, method)
    if (!is_invertible(sigma0))
        sigma0 <- sigma0 + diag(10e-7, p)
    if (!is_invertible(sigma1))
        sigma1 <- sigma1 + diag(10e-7, p)
    a <- drt$a
    ## print(a)
    conv <- drt$conv
    qda.fit <- qda(x %*% a, y)
    ynew <- predict(qda.fit, xnew %*% a)$class
    return(list(class = ynew, drt = a, conv = conv))
}
