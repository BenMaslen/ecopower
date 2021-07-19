#' Multivariate equivalence testing
#'
#' @description
#' \code{equivtest} takes in a copula model fitted to data and a matrix of effect sizes to execute a
#'  a multivariate equivalence test.
#'
#' @details
#' \code{equivtest} takes a \code{\link[ecoCopula]{cord}} object and a coefficient matrix \code{coeffs} which specifies an effect size of
#' interest to perform an equivalence test.
#' 
#' First, marginal parameters of the data are obtained from a \code{\link[mvabund]{manyglm}} object. Next, a copula model is fitted
#' using \code{\link[ecoCopula]{cord}} to estimate the factor analytic covariance structure of the data. The \code{\link[ecoCopula]{cord}} function uses two
#' factors by default. The p-value is then obtained by parsing the \code{\link[ecoCopula]{cord}} object into \code{\link{extend}},
#' \code{nsim} times with an effect size specified by \code{coeffs}.
#'
#' The test statistics are simulated under the hypothesis that the effect size equals a certain threshold.
#' The p-value is computed as the proportion of times the simulated test statistics are less than the observed
#' statistic. Equivalence is declared if the estimated effect is less than the threshold.
#'
#' \code{equivtest} can handle any user-defined null hypothesis, so only the fitted null model (\code{object0}) or the predictor of
#' interest (\code{term}) needs to be specified. If both \code{object0} and \code{term} are \code{NULL}, \code{equivtest} will
#' automatically set the predictor of interest as the last term in the fitted \code{object} model or drop the only term in the model
#' to obtain the intercept model.
#'
#' Simulations are computed in parallel using the "socket" approach, which uses all available cores minus 1 for clustering
#' to improve computation efficiency. Using 1 less than the number of available cores for your
#' machine (\code{detectCores()-1}) is recommended to leave one core available for other computer processes.
#' @param object objects of class \code{cord}, typically the result of a call to \code{\link[ecoCopula]{cord}}.
#' @param coeffs Coefficient matrix for a \code{\link[mvabund]{manyglm}} object that characterises the size of effects to be simulated.
#' See \code{\link{effect_alt}} for help in producing this matrix.
#' @param term Name of predictor of interest in quotes. Defaults to \code{NULL}, see details.
#' @param object0 object of class \code{cord} that specifies the null hypothesis. Defaults to \code{NULL}, see details.
#' @param stats Statistics simulated under the null hypothesis. Optional, defaults to \code{NULL}. If not \code{NULL}, \code{equivtest} will not
#' simulate test statistics and use the \code{stats} specified. 
#' @param test Test statistic for computing p-value. Defaults to \code{"LR"}.
#' @param nsim Number of simulations for p-value estimate to be based upon. Defaults to \code{999}.
#' @param ncores Number of cores for parallel computing. Defaults to the total number of cores available on the
#' machine minus 1.
#' @param show.time Logical. Displays time elapsed. Defaults to \code{TRUE}.
#' @return Equivalence test results, and;
#' \item{\code{p}}{p-value;}
#' \item{\code{stat_obs}}{observed statistic;}
#' \item{\code{stats}}{simulated statistics.}
#' @seealso \code{\link{effect_alt}}
#' @import ecoCopula
#' @import mvabund
#' @import parallel
#' @importFrom stats anova
#' @importFrom stats update
#' @importFrom stats drop.terms
#' @importFrom stats formula
#' @importFrom stats na.omit
#' @importFrom stats terms
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @rdname equivtest
#' @examples
#' \dontrun{
#' library(ecoCopula)
#' library(mvabund)
#' data(spider)
#' spiddat = mvabund(spider$abund)
#' X = data.frame(spider$x)
#'
#' # Specify increasers and decreasers
#' increasers = c("Alopacce", "Arctlute", "Arctperi", "Pardnigr", "Pardpull")
#' decreasers = c("Alopcune", "Alopfabr", "Zoraspin")
#'
#' # Equivalence test for continuous predictor at effect_size=1.5
#' fit.glm = manyglm(spiddat~bare.sand, family="negative.binomial", data=X)
#' threshold = effect_alt(fit.glm, effect_size=1.5,
#'        increasers, decreasers, term="bare.sand")
#' fit.cord = cord(fit.glm)
#' equivtest(fit.cord, coeffs=threshold, term="bare.sand", nsim=99)
#'
#' # Equivalence test for categorical predictor with 4 levels at effect_size=1.5
#' X$Treatment = rep(c("A","B","C","D"),each=7)
#' fit_factors.glm = manyglm(spiddat~Treatment, family="negative.binomial", data=X)
#' threshold = effect_alt(fit_factors.glm, effect_size=1.5,
#'        increasers, decreasers, term="Treatment")
#' fit_factors.cord = cord(fit_factors.glm)
#' equivtest(fit_factors.cord, coeffs=threshold, term="Treatment", nsim=99)
#'
#' # Specify object0
#' object0.glm = manyglm(spiddat~1, family="negative.binomial")
#' object0.cord = cord(object0.glm)
#' equivtest(fit_factors.cord, coeffs=threshold, object0=object0.cord, nsim=99)
#' }
#' @export

equivtest.cord = function(object, coeffs, term=NULL, object0=NULL, 
  stats=NULL, test="LR", nsim=999, ncores=detectCores()-1, show.time=TRUE) {
  
  ptm = proc.time()

  check_equivtest_args(coeffs, object0)
  term = get_term(object$obj, term, object0$obj)
  object0 = get_object0(object$obj, term, object0$obj)
  stats = get_stats(object, coeffs, term, object0, stats, test, nsim, ncores)
  object = object$obj
  anova_obj = anova(object0, object, test=test, nBoot=1, show.time="none")
  stat_obs = anova_obj$table[2,ncol(anova_obj$table)-1]
  p = get_pval(stat_obs, stats, nsim)

  # inputs for default.print.anova.manyglm
  topnote = get_topnote(object0, object, resamp=anova_obj$resamp)
  attr(anova_obj$table, "heading") = c("Equivalence Test Table\n", topnote)
  anova_obj$table[2,ncol(anova_obj$table)] = p
  anova_obj$resamp = "cord"

  out = list(
    p = p,
    stat_obs = stat_obs,
    coefficients = coeffs,
    term = term,
    stats = stats, 
    test = anova_obj$test,
    cor.type = anova_obj$cor.type,
    resamp = anova_obj$resamp,
    nsim = nsim,
    n.bootsdone = nsim,
    p.uni = anova_obj$p.uni,
    table = anova_obj$table
  )

  elapsed = proc.time()[3] - ptm[3]
  print_time(elapsed, show.time)

  class(out) <- c("equivtest.cord", "anova.manyglm")
  return(out)
}

get_term = function(object, term, object0) {
  if (is.null(object0)) {
    if (is.null(term)) {
      if (length(labels(terms(object))) > 1) {
        term = labels(terms(object, keep.order = TRUE))[-1]
      } else {
        term = labels(terms(object, keep.order = TRUE))
      }
    } else {
      term = term
    }
  } else {
    term = term
    if (!is.null(term)) {
      stop("Both object0 and term are not null, specify either object0 or term.")
    }
  }
  return (term)
}

get_object0 = function(object, term, object0) {
  if (is.null(object0)) {
    if (length(labels(terms(object))) > 1) {
      term_position = grep(term, attr(object$terms, "term.labels"))
      object0 = update(
        object, 
        formula = drop.terms(
          object$terms, 
          term_position, 
          keep.response = TRUE
        )
      )
    } else {
      object0 = update(
        object,
        formula = update(
          object$formula, ~1
        )
      )
    }
  } else {
    object0 = object0
  }
  return (object0)
}

get_stats = function(object, coeffs, term, object0, stats, test, nsim, ncores) {
  if (is.null(stats)) {
    stats = rep(NA,nsim)
    N = nrow(object$obj$data)
    coeffs = coeffs
    newdata = NULL
    n_replicate = NULL
    do.fit = TRUE

    ncores = get_ncores(ncores)
    cl = makeCluster(ncores)
    clusterExport(cl, objects(envir = .GlobalEnv), envir = .GlobalEnv)
    clusterExport(cl, objects(envir = environment()), envir = environment())

    clusterEvalQ(cl, {
      library(mvabund)
      library(ecoCopula)
    })

    if (is.null(term)) {
      stats = parSapply(cl, stats, MVApowerstatObj)
    } else {
      stats = parSapply(cl, stats, MVApowerstat, coeffs=coeffs)
    }

    stopCluster(cl)
  } else {
    stats = stats
  }
  return (stats)
}

get_pval = function(stat_obs, stats, nsim) {
  p = rep(NA, length=nsim)
  p = stats < stat_obs + 1e-8
  p = (sum(p)+1) / (nsim + 1)
  return (p)
}
