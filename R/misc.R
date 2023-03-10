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
is_invertible <- function(m) class(try(solve(m), silent = TRUE))[1] == "matrix"
