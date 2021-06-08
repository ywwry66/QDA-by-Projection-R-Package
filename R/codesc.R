codesc_1d <- function(fun, par, hess_size = 1e-5, step_size = 0.1) {
    y0 <- fun(par - hess_size)
    y1 <- fun(par)
    y2 <- fun(par + hess_size)
    fst_coef <- (y2 - y0) / (2 * hess_size) -
        (y0 + y2 - 2 * y1) * par / hess_size ^ 2
    sec_coef <- (y0 + y2 - 2 * y1) / (2 * hess_size ^ 2)
    if (sec_coef < 0 || sec_coef == 0) {
        ## print("negative")
        ## return(- fst_coef / (2 * sec_coef))
        return(par - sign(fst_coef + 2 * sec_coef * par) * step_size)
        ## return(par)
        }
    else
        return(-fst_coef / (2 * sec_coef))
}

codesc_1iter <- function(fun, par, step_size = 0.1) {
    for (i in seq_along(par)) {
        ## value_old <- fun(par)
        f <- function(x) {
            par[i] <- x
            fun(par)
        }
        par[i] <- codesc_1d(f, par[i], step_size = step_size)
        ## value_new <- fun(par)
        ## if (value_new > value_old) print(1)
    }
    ## return(par)
    return(par / sqrt(sum(par ^ 2)))
}

##' @export
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

frank_wolfe_1iter <- function(fun, par, step_size, grad_size = 1e-7) {
    grad <- numeric(length(par))
    for (i in seq_along(par)) {
        temp <- par
        temp[i] <- temp[i] + grad_size
        grad[i] <- (fun(temp) - fun(par)) / grad_size
    }
    v <- -grad / sqrt(sum(grad ^ 2))
    par1 <- par + step_size * (v - par)
    return(par1 / sqrt(sum(par1 ^ 2)))
}

frank_wolfe <- function(fun, par, max_iter = 500,
                        reltol = 1e-7, abstol = 1e-9, ...) {
    fun1 <- function(a) fun(a, ...)
    convergence <- 1
    par <- par / sqrt(sum(par ^ 2))
    value <- numeric(max_iter + 1)
    value[1] <- fun1(par)
    if (value[1] < abstol)
        return(list(par = par, value = value[1], iters = 0, convergence = 0))
    for (i in seq_len(max_iter)) {
        par <- frank_wolfe_1iter(fun1, par, 2 / (i + 2))
        value[i + 1] <- fun1(par)
        if (abs(value[i] - value[i + 1]) < max(abstol, reltol * abs(value[i]))) {
            convergence <- 0
            break
        }
    }
    return(list(par = par, value = value[i + 1], iters = i, convergence = convergence))
}

##' @export
codesc_test <- function(par, sigma0, sigma1, mu0, mu1) {
    fun <- function(a) mis_rate(a, mu0, mu1, sigma0, sigma1)
    rs_BFGS$time <- system.time(rs_BFGS <- optim(par = par, fn = fun, method = "BFGS",
                                                 control = list(maxit = 2000)))
    print("step size 0.1")
    rs_codesc_1$time <- system.time(rs_codesc_1 <- codesc(fun = fun, par = par, max_iter = 1000))
    print("step size 1")
    rs_codesc_2$time <- system.time(rs_codesc_2 <- codesc(fun = fun, par = par, max_iter = 1000, step_size = 1))
    print("step size 10")
    rs_codesc_3$time <- system.time(rs_codesc_3 <- codesc(fun = fun, par = par, max_iter = 1000, step_size = 10))
    rs_frank_wolfe$time <- system.time(rs_frank_wolfe <- frank_wolfe(fun = fun, par = par, max_iter = 2000))
    list(BFGS = rs_BFGS,
         codesc_1 = rs_codesc_1,
         codesc_2 = rs_codesc_2,
         codesc_3 = rs_codesc_3,
         frank_wolfe = rs_frank_wolfe)
}

##' @export
codesc_test_summary <- function(num_runs = 10, seed = 1) {
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
    p <- 50
    s1 <- matrix(rnorm(p ^ 2), nrow = p)
    s2 <- matrix(rnorm(p ^ 2), nrow = p)
    out <- data.frame()
    for (i in seq_len(num_runs)) {
        temp <- codesc_test(rnorm(p), s1 %*% t(s1), s2 %*% t(s2), rep(0, p), rep(1, p))
        out <- rbind(
            out, data.frame(BFGS_time = unname(temp$BFGS$time[3]),
                            BFGS_value = temp$BFGS$value,
                            codesc_1_time = unname(temp$codesc_1$time[3]),
                            codesc_1_value = temp$codesc_1$value,
                            codesc_2_time = unname(temp$codesc_2$time[3]),
                            codesc_2_value = temp$codesc_2$value,
                            codesc_3_time = unname(temp$codesc_3$time[3]),
                            codesc_3_value = temp$codesc_3$value,
                            frank_wolfe_time = unname(temp$frank_wolfe$time[3]),
                            frank_wolfe_value = temp$frank_wolfe$value))
    }
    return(out)
}
