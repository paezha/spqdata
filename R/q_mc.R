#' @title A funcion to calculate Q using Monte Carlo
#'
#' @description This function calculates Q0, a measure of spatial association based on symbolic entropy.
#' @param Y a factor of the same length as the coordinates x
#' @param x coordenadas asociadas a la observación de Y
#' @param m amplitud de las m-historias
#' @param s grado de solapamiento
#' @param nsim number of permutations
#' @usage q_mc(Y, x, m, nsim = 999)
#' @keywords spatial association, qualitative variable, symbolic entropy, symbols
#' @details Aquí Antonio escribe una linda historia
#' @return decir que cosas son las que devuelve
#'   \tabular{ll}{
#'     \code{Q0s_p} \tab value of the statistic of the observed distribution. Symbols p\cr
#'     \code{p.value.p} \tab  the pseudo p-value of the Q0s_p test \cr
#'     \code{efp_symb} \tab frecuencia empírica de los simbolos "p" de cada permutación\cr
#'     \code{Q0s_c} \tab value of the statistic of the observed distribution. Symbols c\cr
#'     \code{p.value.c} \tab  the pseudo p-value of the Q0s_c test \cr
#'     \code{efc_symb} \tab frecuencia empírica de los simbolos "c" de cada permutación\cr
#'     \code{nsim} \tab  nsim simulated values of statistic \cr
#'     }
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz  \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     Testing for spatial association of qualitative data using symbolic dynamics.
#'       \emph{Journal of Geographical Systems}, 12(3), 281-309.
#'     \item López, F., and Páez, A. (2012).
#'     Distribution-free inference for Q(m) based on permutational bootstrapping: an application
#'     to the spatial co-location pattern of firms in Madrid.
#'       \emph{Estadística Española}, 177, 135-156.
#'   }
#' @seealso
#' \code{\link{dgp_spq}}, \code{\link{m_surr_no}},\code{\link{q_symb}}
#'
#' @examples
#' # Example 1
#' N <- 1000
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' listw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(cbind(cx,cy), k=4)))
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' fx <- dgp_spq(x = x, p = p, listw = listw, rho = rho)
#' Q0s <- q_mc(fx = fx, x = x, m = 3, nsim = 199)
#'
#' # Example 2
#'
#' # Load dataset
#' data("FastFood")
#' # Define coordinates
#' x <- cbind(FastFood.sf$Lon,FastFood.sf$Lat)
#' m <- 3
#' Q0s <- q_mc(fx = FastFood.sf$Type, x = x, m = 3, nsim = 199)
q_mc <- function(fx, x, m, nsim = 999,
                    seedinit = 123,
                    distance = "Euclidean") {
  Y <- fx
  if (length(Y) != dim(x)[1])
    stop("La longitud e Y no coincide con la dimensión de las coordenadas")
  k <- nlevels(Y)
  N <- length(Y)
  mdtfull <- sf::st_distance(sf::st_as_sf(x),
                             which = distance)
  # full distance matrix
  ms <- mdtms <- matrix(0, nrow = nrow(mdtfull),
                        ncol = m)
  ms[, 1] <- 1:N
  rownames(mdtms) <- ms[, 1]
  colnames(mdtms) <- NULL
  for (i in 1:N) {
    mdti <- mdtfull[i, ]
    mdtms[i, 1] <- mdti[i]
    # maximum distance ith row
    max_dt_mdti <- mdti[which.max(mdti)]
    mdti[i] <- mdti[i] + max_dt_mdti
    # distance with the same point is always zero...
    for (j in 2:m) {
      indx_mdti <- which.min(mdti)
      ms[i, j] <- indx_mdti
      mdtms[i, j] <- mdti[indx_mdti]
      mdti[indx_mdti] <- mdti[indx_mdti] + max_dt_mdti
    }
  }
  symb <- cr_symb(k, m)
  Q0 <- q_symb_A2(Y, ms, symb)
  set.seed(seedinit)
  mcsamp <- rsample::bootstraps(as.data.frame(as.factor(Y)),
                                  times = nsim)
  Qfull_mc <- purrr::map(mcsamp$splits,
                           q_symb_A2, ms, symb)
  Qfull_stat <- unlist(Qfull_mc)
  Qpmc <- Qfull_stat[names(Qfull_stat) == "qp"]
  Qcmc <- Qfull_stat[names(Qfull_stat) == "qc"]
  mefp_symb <- matrix(0, nrow = nrow(symb$p_symb),
                      ncol = nsim)
  mefc_symb <- matrix(0, nrow = nrow(symb$c_symb),
                      ncol = nsim)
  rownames(mefp_symb) <- names(Q0$efp_symb)
  rownames(mefc_symb) <- names(Q0$efc_symb)
  colnames(mefp_symb) <- paste("sim", 1:nsim, sep = "")
  colnames(mefc_symb) <- paste("sim", 1:nsim, sep = "")
  for (i in 1:nsim) {
    mefp_symb[,i] <- Qfull_mc[[i]]$efp_symb
    mefc_symb[,i] <- Qfull_mc[[i]]$efc_symb
  }
  pvaluemc_p <- sum(Qpmc > Q0$qp) / (nsim + 1)
  pvaluemc_c <- sum(Qcmc > Q0$qc) / (nsim + 1)
  results <- list(Q0$qp, pvaluemc_p,
                  Q0$qc, pvaluemc_c,
                  Q0$qp_symb, Q0$qc_symb,
                  Q0$PSymb, Q0$CSymb,
                  Q0$efp_symb, Q0$efc_symb,
                  Qpmc, Qcmc, mefp_symb, mefc_symb,
                  ms, mdtms, symb, distance)
  names(results) <- c("qp", "pvaluemc_qp",
                      "qc", "pvaluemc_qc",
                      "qp_symb", "qc_symb",
                      "PSymb", "CSymb",
                      "efp_symb", "efc_symb",
                      "qpmc", "qcmc",
                      "efp_symb_mc", "efc_symb_mc",
                      "ms", "mdtms", "symb", "distance")
  return(results)
}
