codesc_1d <- function(fun, par, hess_size = 1e-5, step_size = 0.1) {
    y0 <- fun(par - hess_size)
    y1 <- fun(par)
    y2 <- fun(par + hess_size)
    fst_coef <- (y2 - y0) / (2 * hess_size) -
        (y0 + y2 - 2 * y1) * par / hess_size ^ 2
    sec_coef <- (y0 + y2 - 2 * y1) / (2 * hess_size ^ 2)
    if (sec_coef < 0 || sec_coef == 0) {
        return(par - sign(fst_coef + 2 * sec_coef * par) * step_size)
        }
    else
        return(-fst_coef / (2 * sec_coef))
}

codesc_1iter <- function(fun, par, step_size = 0.1) {
    for (i in seq_along(par)) {
        f <- function(x) {
            par[i] <- x
            fun(par)
        }
        par[i] <- codesc_1d(f, par[i], step_size = step_size)
    }
    return(par / sqrt(sum(par ^ 2)))
}

codesc <- function(fun, par, max_iter = 500,
                   reltol = 1e-7, abstol = 1e-9, step_size = 0.1, ...) {
    fun1 <- function(a) fun(a, ...)
    convergence <- 1
    par <- par / sqrt(sum(par ^ 2))
    value <- numeric(max_iter + 1)
    value[1] <- fun1(par)
    if (value[1] < abstol)
        return(list(par = par, value = value[1], iters = 0, convergence = 0))
    for (i in seq_len(max_iter)) {
        par <- codesc_1iter(fun1, par, step_size)
        value[i + 1] <- fun1(par)
        if (abs(value[i] - value[i + 1]) <
            max(abstol, reltol * abs(value[i]))) {
            convergence <- 0
            break
        }
    }
    return(list(par = par, value = value[i + 1],
                iters = i, convergence = convergence))
}
