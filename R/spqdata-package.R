#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

#' @docType package
#' @name spqdata-package
#' @rdname spqdata-package
#'
#' @details Functionality included in \pkg{spqdata} are:
#'
#' @section Datasets: \pkg{spqdata} includes two different datasets: spain and Fastfood. These sets are used to illustrate the capabilities of different functions. Briefly, their main characteristics are the following \cr
#'   \itemize{
#'     \item The \emph{FastFood} An object of class sf with georeferenced points of a selection of fast food restaurants in Toronto, Canada.
#'      \item The \emph{spain.sp} An object of class sf with the boundaries of Spanish provinces and several socio-economic variables.
#'    }
#'
#' @references Breusch T, Pagan A (1980). The Lagrange multiplier test and its applications to model specification in econometrics. \emph{Review of Economic Studies} 47: 239-254.
#' @references LeSage, J., and Pace, R. K. (2009). \emph{Introduction to spatial econometrics}. Chapman and Hall/CRC.
#' @references López, F.A., Mur, J., and Angulo, A. (2014). Spatial model selection strategies in a SUR framework. The case of regional productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#' @references López, F.A., Martínez-Ortiz, P.J., and Cegarra-Navarro, J.G. (2017). Spatial spillovers in public expenditure on a municipal level in Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#' @references Mur, J., López, F., and Herrera, M. (2010). Testing for spatial effects in seemingly unrelated regressions. \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'
#' @importFrom dplyr group_by rename select
#' @importFrom gtools combinations permutations
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar
#' @importFrom ggplot2 aes element_text labs theme
#' @importFrom magrittr %>%
#' @importFrom Matrix solve
#' @importFrom methods as
#' @importFrom rsample bootstraps analysis
#' @importFrom sf st_coordinates st_distance
#' @importFrom spdep knearneigh knn2nb poly2nb
NULL
