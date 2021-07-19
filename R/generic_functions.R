#' @describeIn effect_alt Specify multivariate effect sizes
#' @export effect_alt
effect_alt = function(object, effect_size, increasers, decreasers, term, K=NULL) {  

  UseMethod("effect_alt")
}

#' @describeIn effect_null Specify null effects for multivariate abundance data
#' @export effect_null
effect_null = function(object, term) {

  UseMethod("effect_null")
}

#' @describeIn equivtest Multivariate equivalence testing
#' @export equivtest
equivtest = function(object, coeffs, term=NULL, object0=NULL, 
  stats=NULL, test="LR", nsim=999, ncores=detectCores()-1, show.time=TRUE) {
  
  UseMethod("equivtest")
}

#' @describeIn extend Simulate or extend multivariate abundance data
#' @export extend
extend = function(object, N=nrow(object$obj$data), coeffs=coef(object$obj),
  newdata=NULL, n_replicate=NULL, do.fit=FALSE, seed=NULL) {

  UseMethod("extend")
}

#' @describeIn powersim Provide power estimates for multivariate abundance models
#' @export powersim
powersim = function(object, coeffs, term, N=nrow(object$obj$data),
  coeffs0=effect_null(object$obj, term), nsim=999, test="score",
  alpha=0.05, newdata=NULL, n_replicate=NULL,
  ncores=detectCores()-1, show.time=TRUE) {
  
  UseMethod("powersim")
}

#' @export
print.powersim.cord = function(x, ...) {
  default.print.powersim.cord(x, ...)
  return (invisible())
}

default.print.powersim.cord = function(x, ...) {
  power = x$power
  cat("Power:", power)
  cat("\n")
  invisible (powersim)
}
