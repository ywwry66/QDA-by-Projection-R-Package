drt_1iter <- function(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
                      lambda = 0, method = "Penalization",
                      par = NULL, max_iter = 1000, optim = "codesc") {
  ## initialization
  if (method != "Penalization" & method != "Thresholding") {
    stop("Wrong sparse method")
  }
  p <- length(mu0)
  a <- rep(0, p)
  if (is.null(par)) {
    ## equal covariance
    par <- solve(sigma) %*% (mu0 - mu1)
    ## equal mean
    m0 <- maxquadratio(sigma0, sigma1)
    m1 <- maxquadratio(sigma1, sigma0)
    par1 <- if (m0$value > m1$value) m0$arg else m1$arg
    if (mis_rate(par1, mu0, mu1, sigma0, sigma1, p0, p1) <
          mis_rate(par, mu0, mu1, sigma0, sigma1, p0, p1)) {
      par <- par1
    }
  }
  target <- function(a, mu0, mu1, sigma0, sigma1, p0, p1, lambda) {
    mis_rate(a, mu0, mu1, sigma0, sigma1, p0, p1) +
      lambda * sum(abs(a)) / sqrt(sum(a^2))
  }
  par <- par / sqrt(sum(par^2))
  ## optimization using specified method
  if (optim == "BFGS") {
    op <- optim(
      par = par, fn = target, method = "BFGS",
      mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1,
      p0 = p0, p1 = p1, lambda = lambda,
      control = list(maxit = max_iter)
    )
  } else if (optim == "codesc") {
    op <- codesc(
      par = par, fun = target, max_iter = max_iter,
      mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1,
      p0 = p0, p1 = p1, lambda = lambda
    )
  } else {
    stop("Wrong optimization method")
  }
  d <- op$par
  conv <- op$convergence
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
                par = NULL, max_iter = 1000, optim = "codesc") {
  if (iter > 1 & lambda == 0) {
    p <- length(mu0)
    a <- matrix(0, p, iter)
    conv <- rep(0, iter)
    drt <- drt_1iter(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
      lambda = lambda, method = method,
      par = par, max_iter = max_iter, optim = optim
    )
    a[, 1] <- drt$a
    conv[1] <- drt$conv
    for (j in 2:iter) {
      Q <- qr.Q(qr(a[, 1:(j - 1)]), complete = TRUE)[, j:p]
      drt <- drt_1iter(
        mu0 = t(Q) %*% mu0, mu1 = t(Q) %*% mu1,
        sigma0 = t(Q) %*% sigma0 %*% Q,
        sigma1 = t(Q) %*% sigma1 %*% Q,
        sigma = t(Q) %*% sigma %*% Q,
        p0, p1, lambda = lambda, method = method,
        par = par, max_iter = max_iter, optim = optim
      )
      d <- Q %*% drt$a
      a[, j] <- d / sqrt(sum(d^2))
      conv[j] <- drt$conv
    }
    return(list(a = a, conv = conv))
  } else {
    return(drt_1iter(mu0, mu1, sigma0, sigma1, sigma, p0, p1,
             lambda = lambda, method = method,
             par = par, max_iter = max_iter, optim = optim
           ))
  }
}

##' Qdadratic Discriminant Analysis by Projection
##'
##' This function runs the Quadratic Discriminant Analysis by
##' Projection (QDAP) method.
##'
##' This function only handles two-class classification problems. It
##' tries to find the direction that minimizes the sample
##' classification error under the heteroscedastic Gaussian
##' assumption. It then projects the data onto the optimal direction,
##' and performs 1-D regular QDA.
##'
##' If not specified by the user, the initial direction for the
##' optimization subroutine is either the LDA direction, or the
##' direction that maximizes the ratio of two quadratic forms induced
##' by the covariance matrices.
##'
##' If the covariance matrices are singular, a tiny scalar matrix is
##' added to them for numerical stability.
##'
##' Multiple rounds of applications of QDAP is implemented (still in
##' beta status). See the argument 'iter'.
##'
##' Penalization/Thresholding is implemented to find sparse optimal
##' direction (still in beta status). See the arguments 'lambda' and
##' 'method'.
##'
##' @param x A matrix containing the predictors of the training data.
##' @param y A 0-1 vector containing the class labels of the training
##'   data.
##' @param xnew A matrix containing the predictors of the test data.
##' @param lambda The tuning parameter used for either the
##'   "Penalization" method or the "Thresholding" method. (Beta)
##' @param iter Number of iterations to apply QDAP. If greater than 1,
##'   this will keep searching the optimal direction in the orthogonal
##'   complement of the previous optimal subspace. (Beta)
##' @param method A method to get sparse optimal
##'   direction. "Penalization" for penalizing over the l2 norm of the
##'   optimal direction, or "Thresholding" for thresholding over each
##'   entry of the optimal direction. (Beta)
##' @param par A vector representing the initial direction for the
##'   optimization subroutine.
##' @param max_iter Maximum number of iterations for the optimization
##'   subroutine to run.
##' @param optim The optimization method used towards the
##'   classification error function. "BFGS" for
##'   Broyden–Fletcher–Goldfarb–Shanno algorithm, or "codesc" for
##'   coordinate descent algorithm.
##' @return If 'xnew' is not supplied, the return value is a list
##'   containing 'qdap_rule', 'drt' and 'conv'. If 'xnew' is supplied,
##'   in addition to these, it also contains 'class': \item{class}{A
##'   0-1 vector containing the predicted class label of the test data
##'   'xnew'.}  \item{qdap_rule}{The QDAP classification rule, a
##'   function that takes a vector of the same dimension as each
##'   training sample, and returns the predicted class label.}
##'   \item{drt}{a vector representing the optimal direction to
##'   project data onto.}  \item{conv}{An integer code. ‘0’ indicates
##'   successful completion of the optimization method.}
##' @examples Iris <- iris[-which(iris$Species == "virginica"), ] #
##'   use the first two species only
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
                 method = "Penalization", par = NULL,
                 max_iter = 1000, optim = "codesc") {
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
  if (!is_invertible(sigma0)) {
    sigma0 <- sigma0 + diag(1e-6, p)
  }
  if (!is_invertible(sigma1)) {
    sigma1 <- sigma1 + diag(1e-6, p)
  }
  sigma <- ((n0 - 1) * sigma0 + (n1 - 1) * sigma1) /
    (n0 + n1 - 2)
  drt <- drt(
    mu0, mu1, sigma0, sigma1, sigma, p0, p1,
    iter, lambda, method, par, max_iter, optim
  )
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
  c0 <- mu0^2 / sigma0 - mu1^2 / sigma1 + log(sigma0 / sigma1) -
    2 * log(p0 / p1)
  qdap_rule <- function(xnew) {
    as.integer(c2 * (xnew %*% a)^2 + c1 * (xnew %*% a) + c0 > 0)
  }
  if (!is.null(xnew)) {
    xnew <- data.matrix(xnew)
    ynew <- apply(xnew, MARGIN = 1, FUN = qdap_rule)
    return(list(class = ynew, qdap_rule = qdap_rule, drt = a, conv = conv))
  } else {
    return(list(qdap_rule = qdap_rule, drt = a, conv = conv))
  }
}

qdap_cv <- function(x, y, xnew = NULL, lambda = c(1, 2, 4, 8, 16),
                    method = "Penalization", max_iter = 1000, optim = "codesc",
                    folds = 5, seed = 2020) {
  if (!is.null(seed)) {
    ## reinstate system seed after simulation
    sys_seed <- .GlobalEnv$.Random.seed
    on.exit({
      if (!is.null(sys_seed)) {
        .GlobalEnv$.Random.seed <- sys_seed
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    })
    set.seed(seed)
  }
  x <- data.matrix(x)
  x0 <- x[which(y == 0), ]
  x1 <- x[which(y == 1), ]
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  x0 <- x0[sample(n0), ]
  x1 <- x1[sample(n1), ]
  folds0 <- cut(seq_len(n0), breaks = folds, labels = FALSE)
  folds1 <- cut(seq_len(n1), breaks = folds, labels = FALSE)
  pred_err <- rep(0, length(lambda))
  par <- rep(list(NULL), folds)
  for (j in seq_along(lambda)) {
    for (i in seq_len(folds)) {
      test_ind0 <- which(folds0 == i)
      test_ind1 <- which(folds1 == i)
      test_n0 <- length(test_ind0)
      test_n1 <- length(test_ind1)
      fit <-
        qdap(
          x = rbind(x0[-test_ind0, ], x1[-test_ind1, ]),
          y = c(
            rep(0, n0 - test_n0),
            rep(1, n1 - test_n1)
          ),
          xnew = rbind(x0[test_ind0, ], x1[test_ind1, ]),
          lambda = lambda[j],
          method = method,
          par = par[[i]],
          max_iter = max_iter,
          optim = optim
        )
      par[[i]] <- fit$drt
      ypred <- fit$class
      pred_err[j] <- pred_err[j] +
        sum(ypred != c(rep(0, test_n0), rep(1, test_n1)))
    }
  }
  lambda_best <- lambda[which.min(pred_err)]
  out <- qdap(x, y, xnew,
    lambda = lambda_best,
    method = method,
    max_iter = max_iter,
    optim = optim
  )
  return(c(out, list(lambda = lambda_best)))
}
