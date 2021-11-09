#' @name Q.test
#' @rdname Q.test
#' @title A function to compute Q test for spatial qualitative data using
#'        asymptotic distribution
#' @usage Q.test(formula = NULL, data = NULL, na.action,
#' m = 3, r = 1, distr = "asymptotic",
#' fx = NULL, coor = NULL, control = list())
#' @param formula a symbolic description of the factor(s).
#' @param data an (optional) data frame or a sf object
#'             with points/multipolygons geometry containing the
#'             variable(s) to be tested.
#' @param fx  a factor  or a matrix
#'            of factors in columns
#' @param m length of m-surrounding (default = 3).
#' @param r maximum overlapping between any two m-surroundings (default = 1).
#' @param distr character. Distribution type "asymptotic" (default) or "mc".
#' @param coor (optional) a 2xN vector with spatial coordinates.
#'             Used when *data* is not a spatial object
#' @param control Optional argument. See Control Argument section.
#' @description A function to compute Q test for spatial qualitative data
#' @usage Q.test(formula = NULL, data = NULL, na.action, fx = NULL, coor = NULL, m = 3,
#' r = 1, distr = "asymptotic", control = list())
#' @details Aquí Antonio escribe una linda historia ....
#' @return An list of two object of the class \code{htest}. Each element of the list return the:
#'   \tabular{ll}{
#'     \code{data.name} \tab a character string giving the names of the data.\cr
#'     \code{statistic} \tab Value of the Q test \cr
#'     \code{N} \tab total number of observations.\cr
#'     \code{R} \tab total number of simbolized observations.\cr
#'     \code{m} \tab length m-surrounding.\cr
#'     \code{r} \tab degree of overlapping.\cr
#'     \code{df} \tab degree of freedom.\cr
#'     \code{distr} \tab type of distribution used to get the significance of the Q test.\cr
#'     \code{type} \tab type of simbols.\cr
#'     }
#' @keywords spatial qualitative data, spatial dependence
#' @section Control arguments:
#' \describe{
#' \item{distance}{character to select the type of distance.
#' Default = "Euclidean" for Cartesian coordinates only: one of Euclidean,
#' Hausdorff or Frechet; for geodetic coordinates,
#' great circle distances are computed (see sf::st_distance())}
#' \item{dtmaxabs}{Delete degenerate surrounding based on the absolute distance between observations.}
#' \item{dtmaxpc}{A value between 0 and 1. Delete degenerate surrounding based on the distance. Delete m-surrounding when the maximum distance between observation
#' is upper than k percentage of maximum distance between anywhere observation.}
#' \item{dtmaxknn}{A integer value 'k'. Delete degenerate surrounding based on the near neighborhood criteria. Delete m-surrounding
#' is a element of the m-surrounding is not include in the set of k near neighborhood of the first element}
#' \item{nsim}{number of simulations for get the Monte Carlo distribution.
#'  Default = 999}
#' \item{seedinit}{seed to select the initial element to star
#' the algorithm to get compute the m-surroundings or to start the simulation}
#' }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Ruiz M, López FA, A Páez. (2010). \emph{Testing for spatial association of qualitative
#'     data using symbolic dynamics}. Journal of Geographical Systems. 12 (3) 281-309
#'   }
#' @export
#' @examples
#'
#' # Case 1: With coordinates
#' rm(list = ls())
#' N <- 100
#' cx <- runif(N)
#' cy <- runif(N)
#' coor <- cbind(cx,cy)
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' listw <- spdep::nb2listw(knn2nb(knearneigh(cbind(cx,cy), k = 4)))
#' fx <- dgp.spq(list = listw, p = p, rho = rho)
#' q.test <- Q.test(fx = fx, coor = coor, m = 3, r = 1)
#' summary(q.test)
#' plot(q.test)
#' print(q.test)
#' q.test.mc <- Q.test(fx = fx, coor = coor, m = 3, r = 1, distr = "mc")
#' summary(q.test.mc)
#' plot(q.test.mc)
#' print(q.test.mc)
#'
#' # Case 2: With a sf object
#' rm(list = ls())
#' data("FastFood")
#' f1 <- ~ Type
#' q.test <- Q.test(formula = f1, data = FastFood.sf, m = c(3, 4),
#' r = c(1, 2, 3), control = list(distance ="Euclidean"))
#' summary(q.test)
#' plot(q.test)
#' print(q.test)

#' # Case 3: With a sf object with isolated areas
#' rm(list = ls())
#' data("Spain")
#' f1 <- ~ Older65 + MenWoman
#' q.test <- Q.test(formula = f1,
#' data = spain.sf, m = 3, r = 1, control = list(seedinit = 1111))
#' summary(q.test)
#' print(q.test)
#' plot(q.test)
#' q.test.mc <- Q.test(formula = f1, data = spain.sf, m = 4, r = 3, distr = "mc", control = list(seedinit = 1111))
#' summary(q.test.mc)
#' print(q.test.mc)
#' plot(q.test.mc)

#' # Case 4: Examples with multipolygons
#' rm(list = ls())
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' qb79 <- quantile(nc$BIR79)
#' nc$QBIR79 <- (nc$BIR79 > qb79[2]) + (nc$BIR79 > qb79[3]) +
#' (nc$BIR79 >= qb79[4]) + 1
#' nc$QBIR79 <- as.factor(nc$QBIR79)
#' plot(nc["QBIR79"], pal = c("#FFFEDE","#FFDFA2", "#FFA93F", "#D5610D"),
#'      main = "BIR79 (Quartiles)")
#' sid79 <- quantile(nc$SID79)
#' nc$QSID79 <- (nc$SID79 > sid79[2]) + (nc$SID79 > sid79[3]) +
#' (nc$SID79 >= sid79[4]) + 1
#' nc$QSID79 <- as.factor(nc$QSID79)
#' plot(nc["QSID79"], pal = c("#FFFEDE","#FFDFA2", "#FFA93F", "#D5610D"),
#'      main = "SID79 (Quartiles)")
#' f1 <- ~ QSID79 + QBIR79
#' lq1nc <- Q.test(formula = f1, data = nc, m = 5, r = 2,
#' control = list(seedinit = 1111, dtmaxpc = 0.5, distance = "Euclidean") )
#' print(lq1nc)
#'
#' lq2nc <- Q.test(formula = f1, data = nc, m = 5, r = 2, control = list(dtmaxpc = 0.2) )
#' print(lq2nc)
#'
#' lq3nc <- Q.test(formula = f1, data = nc, m = 5, r = 2, control = list(dtmaxknn = 5) )
#' print(lq3nc)
#'
#' # Case 5: Examples with points and matrix of variables
#' fx <- matrix(c(nc$QBIR79, nc$QSID79), ncol = 2, byrow = TRUE)
#' mctr <- suppressWarnings(sf::st_centroid(sf::st_geometry(nc)))
#' mcoor <- sf::st_coordinates(mctr)[,c("X","Y")]
#' q.test <- Q.test(fx = fx, coor = mcoor, m = 5, r = 2,control = list(seedinit = 1111, dtmaxpc = 0.5))
#' print(q.test)
#' plot(q.test)
#'
#'
#'
Q.test <- function(formula = NULL, data = NULL, na.action,
                 fx = NULL, coor = NULL,
                 m = 3, r = 1, distr = "asymptotic",
                 control = list()) {
  con <- list(distance = "Euclidean",
              nsim = 999, seedinit = 1111,
              dtmaxabs = 0, dtmaxpc = 0, dtmaxknn = 0)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  distance <- con$distance
  nsim <- con$nsim
  seedinit <- con$seedinit
  dtmaxabs <- con$dtmaxabs
  dtmaxpc <- con$dtmaxpc
  dtmaxknn <- con$dtmaxknn
  cl <- match.call()
   # Lectura Datos.
  if (!is.null(formula) && !is.null(data)) {
    if (inherits(data, "Spatial")) data <- as(data, "sf")
    mfx <- model.frame(formula, data, na.action = na.action)
    #mfx <- get_all_vars(formula, data)
  } else if (!is.null(fx) && !is.null(coor)) {
    mfx <- fx
    if (!is.matrix(mfx) && !is.data.frame(mfx)) mfx <- as.matrix(mfx, ncol = 1)
    mfx <- as.data.frame(mfx)
    if (is.matrix(coor)) coor <- sp::SpatialPoints(coor)
    if (inherits(coor, "Spatial")) coor <- as(coor, "sf")
    data <- coor #sf object
  } else stop("input data wrong")
  for (i in 1:ncol(mfx)) {
    if (!is.factor(mfx[,i])) mfx[,i] <- as.factor(mfx[,i])
  }
  mfx <- as.matrix(mfx)
  lres <- NULL # Initialize the list of results
  # Same random initial observation for every m-sorrounding
  if (!is.null(seedinit)) set.seed(seedinit)
  initobs <- sample(1:nrow(mfx), 1)
  for (i in 1:ncol(mfx)) {
    xfi <- mfx[, i]
    if (!is.factor(xfi)) xfi <- as.factor(xfi)
    N <- length(xfi)
    ki <- length(levels(xfi))
    lms <- vector(mode = "list", length = length(m)) #List of m-sorroundings
    for (j in seq_len(length(m))) {
      mj <- m[j]
      symbi <- cr_symb(ki, mj)
      if (distr == "asymptotic") {
        rcut <- r[!(r > (mj - 1))]
        lms[[j]] <- vector(mode = "list",
                           length = length(rcut))
        for (h in seq_len(length(rcut))) {
          rh <- rcut[h]
           ### M-SURROUNDINGS...
          lms[[j]][[h]] <- m.surround(x = data, m = mj, r = rh,
                                      distance = distance,
                                      control = list(initobs = initobs,
                                                     dtmaxabs = dtmaxabs,
                                                     dtmaxpc = dtmaxpc,
                                                     dtmaxknn = dtmaxknn))
          # if (typems == "no") {
          R <- nrow(lms[[j]][[h]]$ms)
          #   lms[[j]][[h]] <- m.surround(x = data, m = mj, r = rh,
          #                        control = list(seedinit = seedinit,
          #                                       dtmaxabs = dtmaxabs,
          #                                       dtmaxpc = dtmaxpc))
          #   #lms__jkpr <- m_surr_no(x = sf::st_coordinates(data), m = mj, s = rh)
          # } else if (typems == "cdt") {
          #   lms[[j]][[h]] <- m_surr_cdt(x = data, m = mj, r = rh,
          #                        control = list(dtmaxabs = dtmaxabs,
          #                                       dtmaxpc = dtmaxpc))
          # } else if (typems == "cbl") {
          #   lms[[j]][[h]] <- m_surr_cbl(x = data, m = mj, r = rh,
          #                        control = list(dtmaxabs = dtmaxabs,
          #                                       dtmaxpc = dtmaxpc))
          # } else stop("Value for argument 'typems' not valid.")
          qi <- q_symb_A2(xfi, lms[[j]][[h]]$ms, symbi)
          qi$qp_df <- nrow(symbi$p_symb) - 1
          qi$qc_df <- nrow(symbi$c_symb) - 1
          statistic <-  c(qi$qp, qi$qc)
          names(statistic) <- c("Qp", "Qc")
          method <- c("Qp (asymptotic distrib.) for standard symbolization based on permutations",
                      "Qc (asymptotic distrib.) for equivalent symbolization based on combinations")
          if (!is.null(colnames(mfx)[i])) {
            data.name <- rep(colnames(mfx)[i], 2)
          } else data.name <- c(NULL, NULL)
          parameter <- c(qi$qp_df, qi$qc_df)
          names(parameter) <- rep("df", 2)
           p.value <- rep(0, length(parameter))
          for (l in 1:length(p.value)) {
            p.value[l] <- pchisq(statistic[l], df = parameter[l],
                                 lower.tail = FALSE)
          }
          lres_i <- vector(mode = "list", length = length(statistic))
          for (t in 1:length(statistic)) {
            lres_i[[t]]<- list(statistic = statistic[t],
                               parameter = parameter[t],
                               p.value = p.value[t],
                               method = method[t],
                               data.name = paste("Variable",
                                                 data.name[t],
                                                 " m = ",mj," r = ", rh),
                               var.name = data.name[t],
                               N = N,
                               R = R,
                               m = mj, r = rh,
                               k = ki, df = parameter[t])
            lres_i[[t]]$distr <- distr
            lres_i[[t]]$ms <- lms[[j]][[h]]$ms
            lres_i[[t]]$mdtms <- lms[[j]][[h]]$mdtms
            lres_i[[t]]$initobs <- lms[[j]][[h]]$initobs
            lres_i[[t]]$distance <- lms[[j]][[h]]$distance
            if(names(statistic)[t] == "Qp") {
              lres_i[[t]]$type <- "standard-permutations"
              lres_i[[t]]$symb <- symbi$p_symb
              lres_i[[t]]$efp_symb <- qi$efp_symb
              lres_i[[t]]$qp_symb <- qi$qp_symb
              lres_i[[t]]$PSymb <- qi$PSymb
            } else if(names(statistic)[t] == "Qc") {
              lres_i[[t]]$type <- "equivalent-combinations"
              lres_i[[t]]$symb <- symbi$c_symb
              lres_i[[t]]$efc_symb <- qi$efc_symb
              lres_i[[t]]$qc_symb <- qi$qc_symb
              lres_i[[t]]$CSymb <- qi$CSymb
            } else { stop("statistic neither Qp or Qc") }
            lres_i[[t]]$n <- nrow(lres_i[[t]]$symb)
            class(lres_i[[t]]) <- c("htest")
          }
          lres <- c(lres, lres_i)
        } # end for (h in seq_len(length(rcut)))
      } else if (distr == "mc") {
        qi <- q_mc(fx = xfi, x = data,
                      m = mj, nsim = nsim,
                      distance = distance,
                      seedinit = seedinit)
        statistic <-  c(qi$qp, qi$qc)
        names(statistic) <- c("Qp", "Qc")
        method <- c("Qp (mc distrib.) for symbolization based on permutations",
                    "Qc (mc distrib.) for symbolization based on combinations")
        if (!is.null(colnames(mfx)[i])) {
          data.name <- rep(colnames(mfx)[i], 2)
        } else data.name <- c(NULL, NULL)
        parameter <- c("NA", "NA")
        names(parameter) <- rep("df", 2)
        p.value <- c(qi$pvaluemc_qp,
                     qi$pvaluemc_qc)
        lres_i <- vector(mode = "list",
                         length = length(statistic))
        for (t in 1:length(statistic)) {
          lres_i[[t]]<- list(statistic = statistic[t],
                             parameter = parameter[t],
                             p.value = p.value[t],
                             method = method[t],
                             data.name = paste("Variable",
                                               data.name[t],
                                               " m = ",mj,
                                               " r = ", mj-1),
                             var.name = data.name[t],
                             N = N,
                             R = nrow(qi$ms),
                             m = mj, r = mj-1,
                             k = ki, df = parameter[t])
          lres_i[[t]]$distr <- distr
          lres_i[[t]]$ms <- qi$ms
          lres_i[[t]]$mdtms <- qi$mdtms
          lres_i[[t]]$distance <- qi$distance
          if(names(statistic)[t] == "Qp") {
            lres_i[[t]]$type <- "standard-permutations"
            lres_i[[t]]$symb <- symbi$p_symb
            lres_i[[t]]$efp_symb <- qi$efp_symb
            lres_i[[t]]$qp_symb <- qi$qp_symb
            lres_i[[t]]$PSymb <- qi$PSymb
            lres_i[[t]]$qp_mc <- qi$qpmc
            lres_i[[t]]$efp_symb_mc <- qi$efp_symb_mc
          } else if(names(statistic)[t] == "Qc") {
            lres_i[[t]]$type <- "equivalent-combinations"
            lres_i[[t]]$symb <- symbi$c_symb
            lres_i[[t]]$efc_symb <- qi$efc_symb
            lres_i[[t]]$qc_symb <- qi$qc_symb
            lres_i[[t]]$CSymb <- qi$CSymb
            lres_i[[t]]$qc_mc <- qi$qcmc
            lres_i[[t]]$efc_symb_mc <- qi$efc_symb_mc
          } else { stop("statistic neither Qp or Qc") }
          lres_i[[t]]$n <- nrow(lres_i[[t]]$symb)
          class(lres_i[[t]]) <- c("htest")
        }
        lres <- c(lres, lres_i)
      }
    } # end for (j in seq_len(length(m)))
  } # end for (i in 1:ncol(mfx))
  class(lres) <- c("spqtest", class(lres))
  return(lres)
}
