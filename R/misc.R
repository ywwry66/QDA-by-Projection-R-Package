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
