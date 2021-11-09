#' Spain by provinces
#'
#' Data set with areas with isolated areas
#'
#'
#' @format a sf object with 50 rows and 7 variables:
#' \describe{
#'   \item{Province}{Province name.}
#'   \item{Population}{Population.}
#'   \item{Older65}{Factor 3 categorias en funcion del porcen mayores.}
#'   \item{AverageAge}{Ageed average.}
#'   \item{MenWoman}{Factor 2 categorias men si porcen hombres mayor mujeres.}
#'   \item{MassTransitSystems}{Has masive transport sistem.}
#'   \item{Coast}{Has coast.}
#' }
#'
#' @usage data(Spain)
#'
#' @source Páez et al. (2020) \url{https://onlinelibrary.wiley.com/doi/full/10.1111/gean.12241}
#'
#' @references
#'   \itemize{
#'     \item Paez, A., Lopez, F. A., Menezes, T., Cavalcanti, R., & Pitta, M. (2020). \emph{A Spatio‐Temporal Analysis of
#'      the Environmental Correlates of COVID‐19 Incidence in Spain.}. Geographical Analysis.
#'   }
"spain.sf"

#' Distribution Fast-Food restaurants in Toronto
#'
#' @docType data
#'
#' @usage data(FastFood)
#'
#' @format A sf object with 877 rows and 4 variables:
#'
#' \describe{
#'   \item{ID}{Identification}
#'   \item{Lat}{Latitude.}
#'   \item{Lon}{Longitude.}
#'   \item{Type}{Factor with 3 types of restaurant, "H" Hamburguer; "P" Pizza; "S" Shanwich.}
#' }
#'
#' @source Ruiz et al. (2010)
#' \url{https://link.springer.com/article/10.1007/s10109-009-0100-1}
#'
#' @references
#'   \itemize{
#'     \item Ruiz M, López FA, A Páez. (2010). \emph{Testing for spatial association of qualitative
#'     data using symbolic dynamics}. Journal of Geographical Systems. 12 (3) 281-309
#'   }
"FastFood.sf"
