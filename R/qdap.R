drt_1iter <- function(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
                      lambda = 0, method = "Penalization",
                      optim = "codesc") {
    ## initialization
    if (method != "Penalization" & method != "Thresholding")
        stop("Wrong sparse method")
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
            stop("Wrong optimization method")
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

##' Qdadratic Discriminant Analysis by Projection
##'
##' This function runs the Quadratic Discriminant Analysis by Projection (QDAP) method.
##'
##' This function only handles two-class classification problems. It tries to find the direction that minimizes the sample classification error under the assumption that both class follows normal distribution. It then projects the data onto the optimal direction, and performs 1-D regular QDA.
##'
##' The initial direction for the optimization subroutine is either the LDA direction, or the direction that maximizes the ratio of two quadratic forms induced by the covariance matrices.
##'
##' If the covariance matrices are singular, a tiny scalar matrix is added to them for numerical stability.
##'
##' Multiple rounds of applications of QDAP is implemented (still in beta status). See the argument 'iter'.
##'
##' Penalization/Thresholding is implemented to find sparse optimal direction (still in beta status). See the arguments 'lambda' and 'method'.
##'
##' @param x A matrix containing the predictors of the training data.
##' @param y A 0-1 vector containing the class labels of the training data.
##' @param xnew A matrix containing the predictors of the test data.
##' @param lambda The tuning parameter used for either the "Penalization" method or the "Thresholding" method. (Beta)
##' @param iter Number of iterations to apply QDAP. If greater than 1, this will keep searching the optimal direction in the orthogonal complement of the previous optimal subspace. (Beta)
##' @param method A method to get sparse optimal direction. "Penalization" for penalizing over the l2 norm of the optimal direction, or "Thresholding" for thresholding over each entry of the optimal direction. (Beta)
##' @param optim The optimization method used towards the classification error function. "BFGS" for Broyden–Fletcher–Goldfarb–Shanno algorithm, or "codesc" for coordinate descent algorithm.
##' @return If 'xnew' is not supplied, the return value is a list containing 'qdap_rule', 'drt' and 'conv'. If 'xnew' is supplied, in addition to these, it also contains 'class':
##' \item{class}{A 0-1 vector containing the predicted class label of the test data 'xnew'.}
##' \item{qdap_rule}{The QDAP classification rule, a function that takes a vector of the same dimension as each training sample, and returns the predicted class label.}
##' \item{drt}{a vector representing the optimal direction to project data onto.}
##' \item{conv}{An integer code. ‘0’ indicates successful completion of the optimization method.}
##' @examples
##' Iris <- iris[-which(iris$Species == "virginica"), ] # use the first two species only
##'
##' ## Set up training and test data
##' set.seed(2021)
##' n <- nrow(Iris)
##' train <- sample(1:n, n/2)
##' x <- Iris[train, 1:4]
##' y <- as.integer(Iris[train, 5] == "setosa")
##' xnew <- Iris[-train, 1:4]
##' ynew <- as.integer(Iris[-train, 5] == "setosa")
##'
##' ## Calculate classification error
##' fit <- qdap(x, y, xnew)
##' sum(fit$class != ynew)/length(ynew)
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
