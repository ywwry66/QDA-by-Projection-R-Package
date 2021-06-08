drt_1iter <- function(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
                      lambda = 0, method = "Penalization",
                      optim = "codesc") {
    ## initialization
    if (method != "Penalization" & method != "Thresholding")
        stop("Wrong method")
    p <- length(mu0)
    a <- rep(0, p)
    if (is_invertible(sigma))
        par0 <- solve(sigma) %*% (mu0 - mu1)
    if (is_invertible(sigma0) & is_invertible(sigma1)) {
        m0 <- maxquadratio(sigma0, sigma1); m1 <- maxquadratio(sigma1, sigma0)
        par1 <- if (m0$value > m1$value) m0$arg else m1$arg
    }
    if (mis_rate(par1, mu0, mu1, sigma0, sigma1, p0, p1) <
        mis_rate(par0, mu0, mu1, sigma0, sigma1, p0, p1))
        par0 <- par1
    par0 <- par0 / sqrt(sum(par0^2))
    ## optimization wrt lda direction as initial
    if (lambda == 0 | method == "Thresholding") {
        if (optim == "BFGS")
            op <- optim(par = par0, fn = mis_rate, method = "BFGS",
                        mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1,
                        p0 = p0, p1 = p1, control = list(maxit = 1000))
        else if (optim == "codesc")
            op <- codesc(par = par0, fun = mis_rate, max_iter = 1000,
                         mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1,
                         p0 = p0, p1 = p1, step_size = 0.1)
        else if (optim == "frank_wolfe")
            op <- frank_wolfe(par = par0, fun = mis_rate, max_iter = 500,
                              mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1,
                              p0 = p0, p1 = p1)
        else
            stop("Wrong optimization method")
        ## print(op)
        d <- op$par
        conv <- op$convergence
    } else{
        const <- function(a) sum(a^2) - 1
        target <- function(a) mis_rate(a, mu0, mu1, sigma0, sigma1) +
                                  lambda * sum(abs(a))
        op <- nloptr::nloptr(x0 = par0, eval_f = target, eval_g_eq = const,
                     opts = list("algorithm" = "NLOPT_GN_ISRES",
                                 "maxeval" = 300000),
                     lb = rep(-1, p), ub = rep(1, p))
        ## print(op)
        d <- op$solution
        conv <- op$status
    }
    a <- d / sqrt(sum(d^2))
    if (method == "Thresholding" & lambda != 0) {
        a[which(abs(a) < lambda)] <- 0
        a[which(a >= lambda)] <- a[which(a >= lambda)] - lambda
        a[which(a <= -lambda)] <- a[which(a <= -lambda)] + lambda
        a <- a / sqrt(sum(a^2))
    }
    return(list(a = a, conv = conv))
}

drt <- function(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
                iter = 1, lambda = 0, method = "Penalization",
                optim = "codesc") {
    if (iter > 1 & lambda == 0) {
        p <- length(mu0)
        a <- matrix(0, p, iter)
        conv <- rep(0, iter)
        drt <- drt_1iter(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
                         lambda, method, optim)
        a[, 1] <- drt$a
        conv[1] <- drt$conv
        for (j in 2:iter) {
            Q <- qr.Q(qr(a[, 1:(j - 1)]), complete = TRUE)[, j:p]
            drt <- drt_1iter(mu0 = t(Q) %*% mu0, mu1 = t(Q) %*% mu1,
                             sigma0 = t(Q) %*% sigma0 %*% Q,
                             sigma1 = t(Q) %*% sigma1 %*% Q,
                             p0, p1, lambda, method, optim)
            d <- Q %*% drt$a
            a[, j] <- d / sqrt(sum(d^2))
            conv[j] <- drt$conv
        }
        return(list(a = a, conv = conv))
    } else
        return(drt_1iter(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
                         lambda, method, optim))
}

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

qdap <- function(x, y, xnew = NULL, lambda = 0, iter = 1,
                 method = "Penalization", optim = "codesc") {
    x <- data.matrix(x)
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
        sigma0 <- sigma0 + diag(1e-6, p)
    if (!is_invertible(sigma1))
        sigma1 <- sigma1 + diag(1e-6, p)
    sigma <- ((n0 - 1) * sigma0 + (n1 - 1) * sigma1) /
        (n0 + n1 - 2)
    drt <- drt(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
               iter, lambda, method, optim)
    a <- drt$a
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
    qdap_rule <- function(xnew)
        as.integer(c2 * (xnew %*% a) ^ 2 + c1 * (xnew %*% a) + c0 > 0)
    if (!is.null(xnew)) {
        xnew <- data.matrix(xnew)
        ynew <- apply(xnew, MARGIN = 1, FUN = qdap_rule)
        return(list(class = ynew, qdap_rule = qdap_rule, drt = a, conv = conv))
    } else
        return(list(qdap_rule = qdap_rule, drt = a, conv = conv))
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
