#' Specify multivariate effect sizes
#'
#' @description
#' \code{effect_alt} returns a coefficient matrix to be parsed to \code{\link{extend}}, \code{\link{powersim}}
#' and \code{\link{equivtest}} to specify an effect size of interest.
#'
#' @details
#' \code{effect_alt} helps users to create interpretable multivariate effect sizes to be parsed into \code{\link{extend}},
#' \code{\link{powersim}} and \code{\link{equivtest}}, so that researchers can investigate the relationship between effect size, power and
#' sample size in a complicated multivariate abundance setting.
#'
#' \code{effect_alt} creates an effect of size \code{log(effect_size)} for a predictor of interest (\code{term}), for responses
#' who have been specified to increase (\code{increasers}) and \code{-log(effect_size)} for responses who have been specified
#' to decrease (\code{decreasers}). Responses that have not been specified in the \code{increasers} or \code{decreasers} vectors
#' are specified to have no effect with a coefficient of 0. The effect has been logged to make the effect size
#' interpretable within the coefficient matrix.
#'
#' For poisson regression \code{family=poisson()} and negative binomial regression \code{family="negative.binomial"} the effect
#' size is interpreted for a categorical variable as the multiplicative change in mean abundance in the treatment
#' group relative to the control group, whilst for a continuous variable it is interpreted as the multiplicative
#' change in abundance for a 1 unit increase in the predictor of interest.
#'
#' For logit regression \code{family=binomial("logit")} the effect size is interpreted as an odds ratio. For a categorical
#' variable this is the change in odds of obtaining outcome \code{1} when being in the treatment group relative to the
#' control group. Whilst for continuous variables, this is interpreted as the change in odds of obtaining outcome \code{1}
#' with a 1 unit increase in the predictor of interest.
#'
#' For cloglog regression \code{family=binomial("cloglog")} the effect size is interpreted similarly to poisson and
#' negative binomial regression. For a categorical variable it is interpreted as the multiplicative change in
#' the mean of the underlying count in the treatment group relative to the control. Whilst for a continuous
#' variable it is interpreted as the multiplicative change in the mean of the underlying count for a 1 unit
#' increase in the predictor of interest.
#'
#' For categorical variables, the intercept is also changed to be the group mean intercept by taking the
#' intercept of a model without the categorical predictor of interest. This is done to avoid messy comparisons
#' of null control groups.
#'
#' For categorical variables with more than two levels, effect size is changed to \code{effect_size^K[i]} where K defaults
#' to be \code{c(1,2,...,nlevels - 1)}, where \code{nlevels} are the number of levels of the categorical variable and is
#' specified along the order of the levels. To change this, specify a vector \code{K} with length of \code{nlevels - 1}. 
#' To change the control group, this must be done prior to specifying the \code{\link[mvabund]{manyglm}} object
#' using \code{relevel} (which can also change the order of the levels).
#'
#' Note that if the predictor of interest is a categorical variable it must be classed either as a factor or
#' character otherwise results may be misleading.
#'
#' @param object objects of class \code{manyglm}, typically the result of a call to \code{\link[mvabund]{manyglm}}.
#' @param effect_size An effect size of interest, see details for interpretation.
#' @param increasers A vector list of responses which increase relative to the control group/intercept.
#' @param decreasers A vector list of responses which decrease relative to the control group/intercept.
#' @param term Name of predictor of interest in quotes.
#' @param K A vector of length \code{nlevels - 1}. If \code{NULL}, the effect size will increase by its exponent according to
#' the order of factor variables. Alternatively, specify a vector \code{K} that corresponds to the exponent of the
#' \code{effect_size} for each level of a factor variable. Defaults to \code{NULL}, see details.
#' @return A coefficient matrix with the specified effect size.
#' @seealso \code{\link{extend}}, \code{\link{equivtest}}, \code{\link{powersim}}
#' @import mvabund
#' @rdname effect_alt
#' @examples
#' library(mvabund)
#' data(spider)
#' spiddat = mvabund(spider$abund)
#' X = data.frame(spider$x)
#'
#' # Specify increasers and decreasers
#' increasers = c("Alopacce", "Arctlute", "Arctperi", "Pardnigr", "Pardpull")
#' decreasers = c("Alopcune", "Alopfabr", "Zoraspin")
#'
#' # Obtain an effect matrix of effect_size=3
#' spid.glm = manyglm(spiddat~soil.dry, family="negative.binomial", data=X)
#' effect_mat = effect_alt(spid.glm, effect_size=3,
#'          increasers, decreasers, term="soil.dry")
#'
#' # Obtain an effect matrix of effect_size=1.5
#' X$Treatment = rep(c("A","B","C","D"),each=7)
#' spid.glm = manyglm(spiddat~Treatment, family="negative.binomial", data=X)
#' effect_mat = effect_alt(spid.glm, effect_size=1.5,
#'          increasers, decreasers, term="Treatment")
#'
#' # Change effect size parameterisation
#' effect_mat = effect_alt(spid.glm, effect_size=1.5,
#'                          increasers, decreasers, term="Treatment",
#'                          K=c(3,1,2))
#' @export

effect_alt.manyglm = function(object, effect_size, increasers, decreasers, term, K=NULL) {

  coeff = effect_null(object, term)
  rowNames = row.names(coeff)[startsWith(row.names(coeff), term)]
  nRow = length(rowNames)
  nCol = ncol(coeff)
 
  #specify term based on class
  if (class(object$data[,term]) %in% c("integer", "numeric")) {
    coeff[term,][increasers] = log(effect_size)
    coeff[term,][decreasers] = -log(effect_size)
  }

  if (class(object$data[,term]) %in% c("factor", "character")) {
    K = get_K(K, nRow)

    for (iRow in 1:nRow) {
      coeff[rowNames[iRow],][increasers] = log(effect_size^K[iRow])
      coeff[rowNames[iRow],][decreasers] = -log(effect_size^K[iRow])
    }
  }

  return(coeff)
}

get_K = function(K, nRow) {
  if (is.null(K)) {
    K = seq(nRow)
  } else {
    if (length(K) == nRow) {
      K = K
    } else {
      stop("Length of K is not equal to the number of term levels - 1")
    }
  }
  return (K)
}
