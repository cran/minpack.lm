nlsLM <- function (formula, data = parent.frame(), start, jac = NULL, 
                  algorithm = "LM", control = nls.lm.control(), 
                  lower = NULL, upper = NULL, trace = FALSE, subset, weights, na.action, 
                  model = FALSE, ...) 
{
  formula <- as.formula(formula)
  if (!is.list(data) && !is.environment(data)) stop("'data' must be a list or an environment")
  mf <- match.call()
  varNames <- all.vars(formula)
  
  if (length(formula) == 2L) {
    formula[[3L]] <- formula[[2L]]
    formula[[2L]] <- 0
  }

  form2 <- formula
  form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)
  mWeights <- missing(weights)
  
  ## if trace = TRUE, set nls.lm.control$nprint = 1
  if (trace) control$nprint <- 1  
  
  pnames <- if (missing(start)) {
    if (!is.null(attr(data, "parameters"))) {
      names(attr(data, "parameters"))
    }
    else {
      cll <- formula[[length(formula)]]
      func <- get(as.character(cll[[1L]]))
      if (!is.null(pn <- attr(func, "pnames"))) 
        as.character(as.list(match.call(func, call = cll))[-1L][pn])
    }
  } else names(start)

  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  if (length(pnames)) varNames <- varNames[is.na(match(varNames, pnames))]
  lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)), error = function(e) -1)
  
  if (length(varNames)) {
    n <- sapply(varNames, lenVar)
    if (any(not.there <- n == -1)) {
      nnn <- names(n[not.there])
      if (missing(start)) {
        warning("No starting values specified for some parameters.\n", 
                "Initializing ", paste(sQuote(nnn), collapse = ", "), 
                " to '1.'.\n", "Consider specifying 'start' or using a selfStart model")
        start <- as.list(rep(1, length(nnn)))
        names(start) <- nnn
        varNames <- varNames[i <- is.na(match(varNames, nnn))]
        n <- n[i]
      }
      else stop("parameters without starting value in 'data': ", paste(nnn, collapse = ", "))
    }
  } else {
    if (length(pnames) && any((np <- sapply(pnames, lenVar)) == -1)) {
      message("fitting parameters ", paste(sQuote(pnames[np == -1]), collapse = ", "), " without any variables")
      n <- integer()
    } else stop("no parameters to fit")
  }

  respLength <- length(eval(formula[[2L]], data, env))

  if (length(n) > 0L) {
    varIndex <- n%%respLength == 0
    if (is.list(data) && diff(range(n[names(n) %in% names(data)])) > 0) {
      mf <- data
      if (!missing(subset)) warning("argument 'subset' will be ignored")
      if (!missing(na.action)) warning("argument 'na.action' will be ignored")
      if (missing(start)) start <- getInitial(formula, mf)
      startEnv <- new.env(hash = FALSE, parent = environment(formula))
      for (i in names(start)) assign(i, start[[i]], envir = startEnv)
      rhs <- eval(formula[[3L]], data, startEnv)
      n <- NROW(rhs)
      wts <- if (mWeights) rep(1, n) else eval(substitute(weights), data, environment(formula))
    } else {
      mf$formula <- as.formula(paste("~", paste(varNames[varIndex], collapse = "+")), env = environment(formula))
      mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
      mf$lower <- mf$upper <- NULL
      mf[[1L]] <- as.name("model.frame")
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights) model.weights(mf) else rep(1, n)
    }
    if (any(wts < 0 | is.na(wts))) stop("missing or negative weights not allowed")
  } else {
    varIndex <- logical()
    mf <- list(0)
    wts <- numeric()
  }
  
  if (missing(start)) start <- getInitial(formula, mf)
  for (var in varNames[!varIndex]) mf[[var]] <- eval(as.name(var), data, env)
  varNamesRHS <- varNamesRHS[varNamesRHS %in% varNames[varIndex]]

  ## added nls.lm from 'minpack.lm' package
  mf <- c(mf, start)
  lhs <- eval(formula[[2L]], envir = mf)
  m <- match(names(start), names(mf)) 
  .swts <- if (!missing(wts) && length(wts)) sqrt(wts)
  
  FCT <- function(par) {    
    mf[m] <- par
    rhs <- eval(formula[[3L]], envir = mf)
    res <- lhs - rhs
    res <- .swts * res    
    res
  }
  
  NLS <- nls.lm(par = start, fn = FCT, jac = jac, control = control, lower = lower, upper = upper, ...)
  ## previous versions before boundaries were included
  # NLS <- nls.lm(par = start, fn = FCT, jac = jac, control = control, ...)
  start <- NLS$par
        
  ##pass optimized parameters to 'nlsModel'
  m <- stats:::nlsModel(formula, mf, start, wts)
    
  ## => internal 'nls' iterations by 'C_nls_iter' switched off
  ## use 'nls.lm' output for convergence info
  if (NLS$info %in% c(1, 2, 3, 4)) isConv <- TRUE else isConv <- FALSE
  finIter <- NLS$niter   
  finTol <- nls.lm.control()$ftol
  convInfo <- list(isConv = isConv, finIter = finIter, finTol = finTol, 
                    stopCode = NLS$info, stopMessage = NLS$message)
  nls.out <- list(m = m, convInfo = convInfo, data = substitute(data), call = match.call())

  nls.out$call$algorithm <- algorithm
  ## need to use '$tol' parameter from nls.control to make 'predict.nls' work
  nls.out$call$control <- nls.control()
  nls.out$call$trace <- FALSE
  nls.out$na.action <- attr(mf, "na.action")
  nls.out$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  if (model) 
    nls.out$model <- mf
  if (!mWeights) 
    nls.out$weights <- wts
  nls.out$control <- control
  class(nls.out) <- "nls"
  nls.out
}
