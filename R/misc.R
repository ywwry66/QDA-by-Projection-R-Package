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

drt_1iter <- function(mu0, mu1, sigma0, sigma1, sigma,
                      lambda = 0, method = "Penalization") {
    ## initialization
    if (method != "Penalization" & method != "Thresholding")
        stop("Wrong method")
    p <- length(mu0)
    a <- rep(0, p)
    a0 <- ginv(sigma) %*% (mu0 - mu1)
    a0 <- a0 / sqrt(sum(a0^2))
    const <- function(a) sum(a^2) - 1
    target <- function(a) mis_rate(a, mu0, mu1, sigma0, sigma1) +
                              lambda * sum(abs(a))
    ## optimization wrt lda direction as initial
    if (lambda == 0 | method == "Thresholding") {
        op <- optim(par = a0, fn = mis_rate, method = "BFGS",
                    mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1)
        d <- op$par
        conv <- op$convergence
    } else{
        op <- nloptr(x0 = a0, eval_f = target, eval_g_eq = const,
                     opts = list("algorithm" = "NLOPT_GN_ISRES",
                                 "maxeval" = 300000),
                     lb = rep(-1, p), ub = rep(1, p))
        d <- op$solution
        print(op)
        conv <- op$status
    }
    ## optimization wrt another possible initial direction
    if (is_invertible(sigma0) & is_invertible(sigma1)) {
        m0 <- maxquadratio(sigma0, sigma1); m1 <- maxquadratio(sigma1, sigma0)
        if (m0$value > m1$value) a0 <- m0$arg else a0 <- m1$arg
        if (lambda == 0 | method == "Thresholding") {
            op1 <- optim(par = a0, fn = mis_rate, method = "BFGS",
                         mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1)
            if (op1$value < op$value) {
                d <- op1$par
                conv <- op1$convergence
            }
        } else {
            op1 <- nloptr(x0 = a0, eval_f = target, eval_g_eq = const,
                          opts = list("algorithm" = "NLOPT_GN_ISRES",
                                      "maxeval" = 300000),
                          lb = rep(-1, p), ub = rep(1, p))
            if (op1$objective < op$objective) {
                d <- op1$solution
                conv <- op1$status
            }
        }
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

drt <- function(mu0, mu1, sigma0, sigma1, sigma,
                iter = 1, lambda = 0, method = "Penalization") {
    if (iter > 1 & lambda == 0) {
        p <- length(mu0)
        a <- matrix(0, p, iter)
        conv <- rep(0, iter)
        drt <- drt_1iter(mu0, mu1, sigma0, sigma1, sigma, lambda, method)
        a[, 1] <- drt$a
        conv[1] <- drt$conv
        for (j in 2:iter) {
            Q <- qr.Q(qr(a[, 1:(j - 1)]), complete = TRUE)[, j:p]
            drt <- drt_1iter(mu0 = t(Q) %*% mu0, mu1 = t(Q) %*% mu1,
                             sigma0 = t(Q) %*% sigma0 %*% Q,
                             sigma1 = t(Q) %*% sigma1 %*% Q, lambda, method)
            d <- Q %*% drt$a
            a[, j] <- d / sqrt(sum(d^2))
            conv[j] <- drt$conv
        }
        return(list(a = a, conv = conv))
    } else return(drt_1iter(mu0, mu1, sigma0, sigma1, sigma, lambda, method))
}
