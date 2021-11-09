#' Ethnicity, age, marital status of residents in Newark according to 1880 US Census.
#'
#' A dataset containing various attributes of residents in Newark. A historical dataset developed and maintained by Spatial Structures in the Social Sciences, Brown University (contact: Prof. John Logan; john_logan@brown.edu)
#'
#'@docType data
#'
#'@usage data(Newark)
#'
#' @format A data frame with 21520 rows and 10 variables:
#' \describe{
#'   \item{ID}{Unique identifier for observation}
#'   \item{X}{Geocode: Longitude of the observation}
#'   \item{Y}{Geocode: Latitude of the observation}
#'   \item{X_UTM2}{Geocode: X coordinate of the observation in Universal Traverse Mercator, false origin, and jiggled to create unique coordinates for observations}
#'   \item{Y_UTM2}{Geocode: Y coordinate of the observation in Universal Traverse Mercator, false origin, and jiggled to create unique coordinates for observations}
#'   \item{NW}{Ethnicity indicator: American}
#'   \item{IRISH}{Ethnicity indicator: Irish}
#'   \item{GERMAN}{Ethnicity indicator: German}
#'   \item{under30}{Age of respondent is under 30}
#'   \item{mar}{Indicator variable for marital status = married}
#'   \item{usborn}{Respondent was born in the USA}
#'   ...
#' }
#' @source \url{https://s4.ad.brown.edu/Projects/UTP/index.htm}
"Newark"
