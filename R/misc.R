## Maximize the ratio of two full rank quadratic forms
## r(a) := (t(a) %*% sigma0 %*% a)/(t(a) %*% sigma1 %*% a).
## Return a list with the optimal direction and the maximum ratio.
maxquadratio <- function(sigma0, sigma1) {
    ei <- eigen(solve(sigma1) %*% sigma0, symmetric = TRUE)
    arg <- as.vector(ei$vectors[, 1])
    value <- ei$values[1]
    return(list(arg = arg, value = value))
}

## Test whether a matrix is invertible.
## If so, return TRUE; otherwise, return FALSE.
is_invertible <- function(m) class(try(solve(m), silent = T))[1] == "matrix"

## Calculate the precision matrix (generalized inverse of covariance)
## from centered data.
prec <- function(x, df = nrow(x) - 1, tol = 1e-4) {
    x_svd <- svd(x, nu = 0)
    rank <- sum(x_svd$d > tol)
    prec <- x_svd$v[, 1:rank] %*% diag(1 / x_svd$d[1:rank]) %*%
        t(x_svd$v[, 1:rank] %*% diag(1 / x_svd$d[1:rank])) * df
    return(prec)
}
