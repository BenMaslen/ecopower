#' Specify null effects for multivariate abundance data
#'
#' @description
#' \code{effect_null} returns a coefficient matrix to be parsed to \code{\link{powersim}} by default
#' to specify a null effect.
#'
#' @details
#' \code{effect_null} produces a coefficient matrix with a null effect that is specified by setting the parameter
#' estimates of a predictor of interest \code{term} to 0. This function is used by default in \code{\link{powersim}}.
#' Note that intercept values are parameterised as in \code{\link{effect_alt}}.
#' @param object objects of class \code{manyglm}, typically the result of a call to \code{\link[mvabund]{manyglm}}.
#' @param term Name of predictor of interest in quotes.
#' @return A coefficient matrix with the null effect.
#' @seealso \code{\link{effect_alt}}, \code{\link{powersim}}
#' @import mvabund
#' @importFrom stats update
#' @importFrom stats drop.terms
#' @rdname effect_null
#' @examples
#' library(mvabund)
#' data(spider)
#' spiddat = mvabund(spider$abund)
#' X = data.frame(spider$x)
#'
#' # Find null effect size for continuous predictor
#' spid.glm = manyglm(spiddat~soil.dry, family="negative.binomial", data=X)
#' coeffs0 = effect_null(spid.glm, term="soil.dry")
#' @export

effect_null.manyglm = function(object, term) {

  check_effect_args(object, term)
  coeff = object$coefficients

  object0 = get_object0(object, term, NULL)
  coeff0 = object0$coefficients
  if (formula(object0)[-2] == ~1) {
    rownames(coeff0) = "(Intercept)"
  }

  whichRow = which(!rownames(coeff) %in% rownames(coeff0))
  coeff[-whichRow,] = coeff0
  coeff[whichRow,] = 0

  return(coeff)
}
