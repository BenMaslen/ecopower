check_coeffs = function (coeffs) {
  if (is.null(coeffs)) {
    stop("coeffs is null, specify coeffs")
  }
}

check_equivtest_args = function (coeffs, object0) {
  check_coeffs(coeffs)
  
  if (!is.null(object0) && class(object0) != "cord") {
    stop("object0 is not of class cord")
  }
}

check_effect_args = function(object, term) {
  if (!class(object$data[,term]) %in% c("integer", "numeric", "factor", "character")) {
    stop("This function does not recognise this type of term")
  }

  if (attr(object$terms, "intercept") == 0) {
    stop("Model without intercept is not supported")
  }

  if (length(unlist(strsplit(term, ":"))) > 1) {
    stop("Interaction term is not supported")
  }
}

get_ncores = function(ncores) {
  if ("_R_CHECK_LIMIT_CORES_" %in% names(s <- Sys.getenv())) {
    ncores = 2
  } else {
    ncores = ncores
  }
  return (ncores)
}

print_time = function (elapsed, show.time) {
  if (show.time == TRUE) {
    hour = elapsed%/%3600
    minute = elapsed%%3600%/%60
    second = round(elapsed%%3600%%60)
    cat("Time elapsed:", hour, "hr", minute, "min", second, "sec\n")
  }
}

#' @importFrom stats formula
get_topnote = function(object, ..., resamp) {
  objects <- list(object, ...)
  dots <- list(...)
  nModels = length(objects)
  nRows <- nrow(object$y)
  nVars <- ncol(object$y)
  nParam <- ncol(object$x)
  cor.type = object$cor.type

  if (cor.type == "R") {
    corrnum <- 0
    if ( nVars > nRows ) # p>N
       warning("number of variables is greater than number of parameters so R cannot be estimated reliably -- suggest using cor.type='shrink'.")
  }
  else if (cor.type == "I") corrnum <- 1
  else if (cor.type == "shrink") corrnum <- 2
  else stop("'cor.type' not defined. Choose one of 'I', 'R', 'shrink'")

  if (resamp=="case") resampnum <- 0  #case
  # To exclude case resampling
  #if (resamp=="case") stop("Sorry, case resampling is not yet available.")
  else if (substr(resamp,1,4)=="resi") resampnum <- 1  # residual
  else if (resamp=="score") resampnum <- 2  # score
  else if (substr(resamp,1,4) =="perm") resampnum <- 3 # permuation
#    else if (substr(resamp,1,1) =="f") resampnum <- 4 # free permuation
  else if (substr(resamp,1,4) ==  "mont") resampnum <- 5 # montecarlo
  else if (substr(resamp,1,3) ==  "pit") resampnum <- 8 # PIT residual bootstrap
  else stop("'resamp' not defined. Choose one of 'case', 'resid', 'score', 'perm.resid', 'montecarlo', 'pit.trap'")

  if (nModels==1) {
    # test the significance of each model terms
    X <- object$x
    varseq <- object$assign
    resdev <- resdf <- NULL
    tl <- attr(object$terms, "term.labels")
    # if intercept is included
    if (attr(object$terms,"intercept")==0) {
        minterm = 1
        nterms = max(1, varseq)
    } else {
        minterm = 0
        nterms <- max(0, varseq)+1
        tl <- c("(Intercept)", tl)
    }
    tl <- tl[1 + unique(object$assign)] # attempt to deal with bug
    if ( nParam==1 )
        stop("An intercept model is comparing to itself. Stopped")

    XvarIn <- matrix(ncol=nParam, nrow=nterms, 1)
    for ( i in 0:(nterms-2)) { # exclude object itself
        XvarIn[nterms-i, varseq>i+minterm] <- 0 # in reversed order
        ncoef <- nParam-length(varseq[varseq>i+minterm])
        resdf <- c(resdf, nRows-ncoef)
    }

    resdf <- c(resdf, object$df.residual)
    # get the shrinkage estimates
    tX <- matrix(1, nrow=nRows, ncol=1)
    if (corrnum==2 | resampnum==5){ # shrinkage or montecarlo bootstrap
        # shrink.param <- c(rep(NA, nterms))
        # use a single shrinkage parameter for all models
        if (object$cor.type == "shrink") {
            # shrink.param[1] <- object$shrink.param
            shrink.param <- rep(object$shrink.param,nterms)
        } else  {
            # shrink.param[1] <- ridgeParamEst(dat=object$residuals, X=tX,
            # only.ridge=TRUE)$ridgeParam
            lambda <- ridgeParamEst(dat=object$residuals, X=tX,
                        only.ridge=TRUE)$ridgeParam
            shrink.param <- rep(lambda, nterms)
        }
        # for ( i in 0:(nterms-2)){ # exclude object itself
        #     fit <- .Call("RtoGlm", modelParam, Y, X[,varseq<=i+minterm,drop=FALSE],
        #        PACKAGE="mvabund")
        #     shrink.param[nterms-i] <- ridgeParamEst(dat=fit$residuals,
        #           X=tX, only.ridge=TRUE)$ridgeParam # in reversed order
        #  }
    } 
    else if (corrnum == 0) shrink.param <- c(rep(1, nterms))
    else if (corrnum == 1) shrink.param <- c(rep(0, nterms))
    # resdev <- c(resdev, object$deviance)
    nModels <- nterms
    ord <- (nterms-1):1
#        topnote <- paste("Model:", deparse(object$call))
    topnote <- paste("Model:", deparse(formula(object), width.cutoff=500),
                     collapse = "\n")
  } else {
    targs <- match.call(expand.dots = FALSE)
    if (targs[[1]] == "example" || any(class(object) == "traitglm"))
        modelnamelist <- paste("Model", format(1:nModels))
    else
        modelnamelist <- as.character(c(targs[[2]], targs[[3]]))

    resdf   <- as.numeric(sapply(objects, function(x) x$df.residual))
    ####### check input arguments #######
    # check the order of models, so that each model is tested against the next smaller one
    ord <- order(resdf, decreasing=TRUE)
    objects <- objects[ord]
    resdf <- resdf[ord]
    modelnamelist <- modelnamelist[ord]

    # get the shrinkage estimates
    if (corrnum == 2 | resampnum == 5) { # shrinkage or parametric bootstrap
        shrink.param <- c(rep(NA,nModels))
            tX <- matrix(1, nrow=nRows, ncol=1)
        for ( i in 1:nModels ) {
            if (objects[[i]]$cor.type == "shrink")
                    shrink.param[i] <- objects[[i]]$shrink.param
            else shrink.param[i] <- ridgeParamEst(dat=objects[[i]]$residuals, X=tX, only.ridge=TRUE)$ridgeParam
        }
    } else if (corrnum == 0) {
        shrink.param <- c(rep(1,nModels))
    } else if (corrnum == 1) {
        shrink.param <- c(rep(0,nModels))
    }
    # Test if the input models are nested, construct the full matrix
    Xnull <- as.matrix(objects[[1]]$x, "numeric")
    nx <- dim(Xnull)[2]
    for ( i in 2:nModels ) {
        XAlt  <- as.matrix(objects[[i]]$x, "numeric")
        Xarg  <- cbind(XAlt, Xnull)
        tmp <- qr(Xarg)
        Xplus <- qr(XAlt)
        if ( tmp$rank == Xplus$rank ) {
            Beta <- qr.coef(Xplus, Xnull) 
            # equivalent to (XAlt\XNull) in matlab
            # The following gets the left null space of beta, ie.LT=null(t(beta));
            # note that LT is an orthogonal complement of Beta, and [Beta, LT] together forms the orthogonal basis that span the column space of XAlt
            # For some reason, it must be null(beta) instead of null(t(beta)) in R to get the same answer in matlab.
            # In case of redundant terms in the model
            RowToOmit <- which(rowSums(is.na(Beta))>0)
            Beta <- na.omit(Beta)
            tmp <- qr(Beta)
            set <- if(tmp$rank == 0) 1:ncol(Beta) else  - (1:tmp$rank)
            LT <- qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
            # Update the next null matrix with the expanded one
            if ( length(RowToOmit)>0 )
                Xnull <- cbind(Xnull, XAlt[,-RowToOmit]%*%LT)
            else
                Xnull <- cbind(Xnull, XAlt%*%LT)

            # Update the Xnull dimension
            nx <- rbind(nx, dim(Xnull)[2])
        }
    }
    X <- Xnull  # full matrix
    # Construct XvarIn
    nParam <- dim(X)[2]
    XvarIn <- matrix(nrow=nModels, ncol=nParam, as.integer(0))
    for ( i in 1:nModels ) {
        XvarIn[nModels+1-i, 1:nx[i]] <- as.integer(1)
    }
    Xnames <- list()   # formula of each model
    Xnames <- lapply(objects, function(x) paste(deparse(formula(x),
                     width.cutoff=500), collapse = "\n"))
    topnote <- paste(modelnamelist, ": ", Xnames, sep = "", collapse = "\n")
    tl <- modelnamelist
    if (tl[1]==tl[2]) {
        warning(paste("Two identical models. Second model's name changed to ", tl[2], "_2", sep=""))
        tl[2] <- paste(tl[2], "_2", sep="")
    }
    ord <- (nModels-1):1
  }
  return (topnote)
}
