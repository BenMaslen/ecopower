#' Crayweed dataset
#'
#' Dataset of fish abundances recorded at crayweed reference and restored sites.
#'
#' @name crayweed
#' @docType data
#'
#' @usage data(crayweed)
#'
#' @format An object of class \code{"list"} containing:
#' \describe{
#'   \item{abund}{A matrix with 27 observations of abundance of 34 fish species.}
#'   \item{X}{A data frame with treatment and time variables.}
#' }
#' @details
#' The matrix \code{abund} has the following species abundances:
#' \itemize{
#'   \item Abudefduf.sp.
#'   \item Acanthopagrus.australis
#'   \item Acanthurus.nigrofuscus
#'   \item Achoerodus.viridis
#'   \item Aplodactylus.lophodon
#'   \item Atypichthys.strigatus
#'   \item Cheilodactylus.fuscus
#'   \item Chromis.hypsilepis
#'   \item Crinodus.lophodon
#'   \item Girella.elevata
#'   \item Girella.tricuspidata
#'   \item Hypoplectrodes.maccullochi
#'   \item Needle.fish.unidentified
#'   \item Notolabrus.gymnogenis
#'   \item Odax.cyanomelas
#'   \item Olisthops.cyanomelas
#'   \item Ophthalmolepis.lineolatus
#'   \item Parma.microlepis
#'   \item Parma.unifasciata
#'   \item Pempheris.compressa
#'   \item Pempheris.multiradiata
#'   \item Pictilabrus.laticlavius
#'   \item Prionurus.microlepidotus
#'   \item Pseudocaranx.dentex
#'   \item Pseudojuloides.elongatus
#'   \item Pseudolabrus.gymnogenis
#'   \item Sardinops.neopilchardus
#'   \item Scorpis.lineolatus
#'   \item Seriola.lalandi
#'   \item Sphyraena.obtusata
#'   \item Tetractenos.hamiltoni
#'   \item Trachinops.taeniatus
#'   \item Trachurus.novaezelandiae
#'   \item Upeneichthyes.lineatus
#' }
#' The data frame \code{X} has the following variables:
#' \itemize{
#'   \item treatment - reference / restored
#'   \item time - sample period with seven time points 
#' }
#' 
#' @references Data attributed to the crayweed restoration project (\url{http://www.operationcrayweed.com/}).
#' @keywords datasets
#' @examples
#'
#' data(crayweed)
#' head(crayweed$abund)
#' head(crayweed$X)
#'
NULL
