##' Maximizing the ratio of two full rank quadratic forms
##'
##' Find the direction `a' that maximizes (t(a) %*% sigma0 %*% a)/(t(a) %*% sigma1 %*% a)
##' @title Maximizing the ratio of two full rank quadratic forms
##' @param sigma0 The matrix which determines the top full rank quadratic form
##' @param sigma1 The matrix which determines the bottom full rank quadratic form
##' @return
##' @author Ruiyang Wu
maxquadratio <- function(sigma0, sigma1) {
    ei <- eigen(solve(sigma1) %*% sigma0, symmetric = TRUE)
    arg <- as.vector(ei$vectors[, 1])
    value <- ei$values[1]
    return(list(arg = arg, value = value))
}

is_invertible <- function(m) class(try(solve(m), silent = T))[1] == "matrix"

drt_1iter <- function(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
                      lambda = 0, method = "Penalization",
                      optim = "BFGS") {
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
        else
            op <- frank_wolfe(par = par0, fun = mis_rate, max_iter = 500,
                              mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1,
                              p0 = p0, p1 = p1)
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
                optim = "BFGS") {
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
