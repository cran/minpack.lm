####

.First.lib <- function(lib, pkg)
    library.dynam("minpack.lm", pkg, lib)

nls.lm <- function(par, fn, jac = NULL, control = list(), ...)
{
    fn1  <- function(par) fn(par, ...)
    jac1 <- if (!is.null(jac))
        function(par) jac(par, ...)

    con <- list(ftol = sqrt(.Machine$double.eps),
                ptol = sqrt(.Machine$double.eps),
                gtol = 0,
                diag = numeric(),
                epsfcn = 0,
                factor = 100,
                maxfev = 100*(length(par) + 1),
                nprint = 0)
    con[names(control)] <- control
    out <- .Call("nls_lm", par, fn1, jac1, con, new.env(), PACKAGE = "minpack.lm")

    out$hessian <- matrix(out$hessian, nrow = length(par))

    names(out$par)        <-
    rownames(out$hessian) <-
    colnames(out$hessian) <-
    names(out$diag)       <- names(par)

    out
}
