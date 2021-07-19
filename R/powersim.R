#' Provide power estimates for multivariate abundance models
#'
#' @description
#' \code{powersim} returns a power estimate for a \code{\link[ecoCopula]{cord}} object for a given sample size \code{N}
#' and effect size of interest.
#'
#' @details
#' \code{powersim} takes a \code{\link[ecoCopula]{cord}} object, sample size \code{N} and coefficient matrix \code{coeffs} which
#' specifies an effect size of interest and returns a power estimate.
#'
#' The power estimate is obtained by first parsing the \code{\link[ecoCopula]{cord}} object into \code{\link{extend}},
#' \code{nsim} times with an effect size specified by \code{coeffs}. Next, the \code{\link[ecoCopula]{cord}} object is parsed into
#' \code{\link{extend}} an additional \code{nsim} times with a null effect, which is defined by default by
#' \code{\link{effect_null}}. This effectively simulates \code{nsim} \code{manyglm} models under both the null
#' and alternative hypothesis.
#'
#' For each simulated \code{\link[mvabund]{manyglm}} object a test statistic \code{test} is obtained. A critical test statistic
#' is then obtained as the upper 1 - \code{alpha} quantile of simulated test statistics under the null
#' hypothesis. Power is then estimated as the proportion of times the test statistics simulated under
#' the alternative hypothesis exceed the critical test statistic under the null.
#'
#' To improve computation time, simulations are computed in parallel using the "socket" approach, which
#' by default uses all available cores minus 1 for clustering. Using 1 less than the number of available cores for your
#' machine (\code{detectCores()-1}) is recommended to better avoid error relating to clustering or nodes.
#' @param object objects of class \code{cord}, typically the result of a call to \code{\link[ecoCopula]{cord}}.
#' @param coeffs Coefficient matrix for a \code{\link[mvabund]{manyglm}} object that characterises the size of effects to be simulated.
#' See \code{\link{effect_alt}} for help in producing this matrix.
#' @param term Name of predictor of interest in quotes.
#' @param N Number of samples for power estimate. Defaults to the number of observations in the original sample.
#' @param coeffs0 Coefficient matrix under the null hypothesis. Defaults to being specified by \code{\link{effect_null}}.
#' @param nsim Number of simulations for power estimate to be based upon. Defaults to \code{999}.
#' @param test Test statistic for power estimate to based upon. Defaults to \code{"score"}, however \code{"wald"} is also allowed.
#' @param alpha Type I error rate for power estimate, defaults to \code{0.05}.
#' @param newdata Data frame of the same size as the original data frame from the \code{\link[ecoCopula]{cord}} object
#' (\code{object$obj$data}), that specifies a different design of interest.
#' @param n_replicate Number of unique replicates of the original data frame. Defaults to \code{NULL}, overwrites \code{N} if specified.
#' @param ncores Number of cores for parallel computing. Defaults to the total number of cores available on the
#' machine minus 1.
#' @param show.time Logical. Displays time elapsed. Defaults to \code{TRUE}.
#' @return Power estimate result, and;
#' \item{\code{power}}{power.}
#' @seealso \code{\link{effect_alt}}, \code{\link{effect_null}}, \code{\link{extend}}
#' @import ecoCopula
#' @import mvabund
#' @import parallel
#' @importFrom stats anova
#' @importFrom stats quantile
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @rdname powersim
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
#' # Find power for continuous predictor at effect_size=1.5
#' fit.glm = manyglm(spiddat~bare.sand, family="negative.binomial", data=X)
#' effect_mat = effect_alt(fit.glm, effect_size=1.5,
#'        increasers, decreasers, term="bare.sand")
#' fit.cord = cord(fit.glm)
#' powersim(fit.cord, coeffs=effect_mat, term="bare.sand", nsim=99)
#'
#' # Find power for categorical predictor with 4 levels at effect_size=1.5
#' X$Treatment = rep(c("A","B","C","D"),each=7)
#' fit_factors.glm = manyglm(spiddat~Treatment, family="negative.binomial", data=X)
#' effect_mat = effect_alt(fit_factors.glm, effect_size=1.5,
#'        increasers, decreasers, term="Treatment")
#' fit_factors.cord = cord(fit_factors.glm)
#' powersim(fit_factors.cord, coeffs=effect_mat, term="Treatment", nsim=99)
#'
#' # Change effect size parameterisation
#' effect_mat = effect_alt(fit_factors.glm, effect_size=1.5,
#'                          increasers, decreasers, term="Treatment",
#'                          K=c(3,1,2))
#' powersim(fit_factors.cord, coeffs=effect_mat, term="Treatment", nsim=99)
#' }
#' @export

powersim.cord = function(object, coeffs, term, N=nrow(object$obj$data),
   coeffs0=effect_null(object$obj, term), nsim=999, test="score",
   alpha=0.05, newdata=NULL, n_replicate=NULL,
   ncores=detectCores()-1, show.time=TRUE) {

  check_coeffs(coeffs)
  stats.null = stats = rep(NA,nsim)
  do.fit = TRUE

  ptm = proc.time()

  ncores = get_ncores(ncores)
  cl = makeCluster(ncores)
  clusterExport(cl, objects(envir = .GlobalEnv), envir = .GlobalEnv)
  clusterExport(cl, objects(envir = environment()), envir = environment())
  libraries = clusterEvalQ(cl, {
    library(mvabund)
    library(parallel)
  })

  # Find critical value under the null
  stats.null = parSapply(cl, stats.null, MVApowerstat, coeffs=coeffs0)
  criticalStat = quantile(
    unlist(stats.null[!is.na(unlist(stats.null))]), 1-alpha, na.rm=TRUE
  )
  # criticalStat_lower = quantile(
  #   unlist(stats.null[!is.na(unlist(stats.null))]), 1-alpha+0.025, na.rm=TRUE
  # )
  # criticalStat_upper = quantile(
  #   unlist(stats.null[!is.na(unlist(stats.null))]), 1-alpha-0.025, na.rm=TRUE
  # )
  
  # Observe the proportion of times our test statistics exceed this value
  stats = unlist(parSapply(cl, stats, MVApowerstat, coeffs=coeffs))
  
  stopCluster(cl)
  
  #print(c(mean(stats[!is.na(stats)]>criticalStat),mean(stats[!is.na(stats)]>criticalStat_lower),mean(stats[!is.na(stats)]>criticalStat_upper),new))
  power = get_power(criticalStat, stats, nsim)

  out = list(power = power)

  elapsed = proc.time()[3] - ptm[3]
  print_time(elapsed, show.time)

  #print(binom.confint(table((stats[!is.na(stats)]>criticalStat))[2],n=nsim,method="wilson")[c(5,6)])
  class(out) = "powersim.cord" 
  return (out)
}

get_power = function(criticalStat, stats, nsim) {
  p = rep(NA, length=nsim)
  p = stats + 1e-8 > criticalStat
  p = (sum(p)+1) / (nsim + 1)
  return (p)
}
