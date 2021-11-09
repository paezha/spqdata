#' @name methods_scantest
#' @title Methods for class scan
#' @description Cambiar bla-bla por localsrq The \code{plot()} function allows the user to plot both beta and spatial
#'  coefficients for all equations of the spsur model. The argument
#'  \code{viewplot} is used to choose between interactive or non-interactive
#'  plots. The \code{print()} function is used to print short tables including the values of beta and
#'  spatial coefficients as well as p-values of significance test for each
#'  coefficient.
#'  \code{\link{summary.localsrq}} when a table (gt format) is needed.
#'  The rest of methods works in the usual way.
#'
#' @param object a \code{scan} object created by \code{\link{scan.test}}.
#' @param x similar to \code{object} argument for \code{print()}
#'  and \code{plot} functions.
#' @param viewplot logical value to show interactively the plots.
#'   Default = \code{TRUE}
#' @param ... further arguments passed to or from other methods.
#' @examples
#' rm(list = ls()) # Clean memory
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     A test for global and local homogeneity of categorical data based on spatial runs.
#'       \emph{Geographical Analysis}.
#'   }
#'

NULL

#' @name plot.scantest
#' @rdname methods_scantest
#' @export
#'
plot.scantest <- function(object,  sf = NULL, coor = NULL, viewplot = TRUE,...){
  z <- object
  if (!inherits(z, "scantest")) stop("Argument must be a scantest object")
  # if (is.null(z$coor)) stop("Include the sf object to generate the plot")
  a <- matrix(0,ncol = 1,nrow = z$N)
  a[z$MLC] <- 1

  #####################
  ### Plot
  #####################
  if (!is.null(z$coor)){
    coor <- as.data.frame(z$coor)
    sf <- st_as_sf(coor,coords = names(coor))
    mysize = 4
  }
  if (!is.null(sf)){
    if (class(st_geometry(sf))[1]=="sfc_MULTIPOLYGON") mysize = .2
    if (class(st_geometry(sf))[1]=="sfc_POINT") mysize = 4
  }

  sf$levels <- addNA(a)
  if (z$distr=="multinomial"){
    levels(sf$levels)[levels(sf$levels)=="0"] <- "non-sig"
    levels(sf$levels)[levels(sf$levels)=="1"] <- "sig"
    cols <- c("non-sig" = "white", "sig" = "red")

  }
  if (z$distr=="bernoulli"){
    browser()
    levels(sf$levels)[levels(sf$levels)=="0"] <- "non-sig"
    if (z$cases.expect > z$cases.observ)
      levels(sf$levels)[levels(sf$levels)=="1"] <- "sig High"
    if (z$cases.expect < z$cases.observ)
      levels(sf$levels)[levels(sf$levels)=="2"] <- "sig Low"
    cols <- c("non-sig" = "white", "sig High" = "red", "sig Low" = "blue")
  }
  lplot_runs <- ggplot(sf) +
    geom_sf(aes(fill = levels), color = "black", shape = 21, size = mysize) +
    theme_bw() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
    xlab(paste0("Significance p-value = ", z$p.value)) +
    scale_fill_manual(values = cols, na.value ="orange",drop = FALSE)

  lplot_runs
}

#' @name summary.scantest
#' @rdname methods_scantest
#'
#' @export
summary.scantest <- function(object, ...) {
  z <- object

  # Print output
  if (z$distr=="bernoulli"){
  cat("\nSummary of data:\n")
  cat(paste0("Distribution....................: ",z$distr,"\n"))
  cat(paste0("Number of locations.............: ",z$N,"\n"))
  cat(paste0("Total number of cases...........: ",z$N,"\n"))
  cat("\nScan statistic:\n")
  cat(paste0("Type of cluster (alternative)...: ",z$alternative,"\n"))
  cat(paste0("Value of statisitic (loglik)....: ",round(z$statistic, digits = 2) ,"\n"))
  cat(paste0("p-value.........................: ",round(z$p.value, digits = 3),"\n"))
  cat("\nIDs of cluster detect:\n")
  cat("Location IDs included..................: ",z$MLC)
  }
  if (z$distr=="multinomial"){
    cat("\nSummary of data:\n")
    cat(paste0("Distribution....................: ",z$distr,"\n"))
    cat(paste0("Number of locations.............: ",z$N,"\n"))
    cat(paste0("Total number of cases...........: ",z$N,"\n"))
    cat("\nScan statistic:\n")
    cat("Observed cases in the MLC...: ")
    cat(round(z$cases.expect, digits = 2))
    cat("\n")
    cat("Expected cases in the MLC...: ")
    cat(z$cases.observ)
    cat("\n")
    cat(paste0("Value of statisitic (loglik)....: ",round(z$statistic, digits = 2) ,"\n"))
    cat(paste0("p-value.........................: ",round(z$p.value, digits = 3),"\n"))
    cat("\nIDs of cluster detect:\n")
    cat("Location IDs included..................: ",z$MLC)
  }
}
