#' @name methods_localsrq
#' @title Methods for class localsrq
#' @description Cambiar bla-bla por localsrq The \code{plot()} function allows the user to plot both beta and spatial
#'  coefficients for all equations of the spsur model. The argument
#'  \code{viewplot} is used to choose between interactive or non-interactive
#'  plots. The \code{print()} function is used to print short tables including the values of beta and
#'  spatial coefficients as well as p-values of significance test for each
#'  coefficient.
#'  \code{\link{summary.localsrq}} when a table (gt format) is needed.
#'  The rest of methods works in the usual way.
#'
#' @param object a \code{localsrq} object created by \code{\link{qtest}}.
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

#' @name print.localsrq
#' @rdname methods_localsrq
#' @export
print.localsrq <- function(object,...) {
  if (!inherits(object, "localsrq")) stop("Argument must be a localsrq object")
print(object$local.SRQ)
invisible(object)
  }

#' @name plot.localsrq
#' @rdname methods_localsrq
#' @export
plot.localsrq <- function(object,  sf = NULL, coor = NULL,  sig = 0.05, viewplot = TRUE,...){
  if (!inherits(object, "localsrq")) stop("Argument must be a localsrq object")
  lsrq <- object
  if (!is.null(lsrq$nsim)){
    a <- as.factor((lsrq$local.SRQ$pseudo.value < sig)*(lsrq$local.SRQ$zseudo.value > 0)*1 +
                     (lsrq$local.SRQ$pseudo.value < sig)*(lsrq$local.SRQ$zseudo.value < 0)*2)
  }
  else {
    a <- as.factor((lsrq$local.SRQ$p.value < sig)*(lsrq$local.SRQ$z.value > 0)*1 +
                     (lsrq$local.SRQ$p.value < sig)*(lsrq$local.SRQ$z.value < 0)*2)
  }
  #####################
  ### Plot SR Local
  #####################
  if (is.null(sf)){
    if (is.null(coor) && (class(lsrq$listw)[1] == "knn")){
    coor <- as.data.frame(lsrq$listw$x)
    }
    if (!is.null(coor) && (class(lsrq$listw)[1] != "knn")){
      coor <- as.data.frame(coor)
    }
    sf <- st_as_sf(coor,coords = names(coor))
    mysize = 4
  }
  if (!is.null(sf)){
    if (class(st_geometry(sf))[1]=="sfc_MULTIPOLYGON") mysize = .2
    if (class(st_geometry(sf))[1]=="sfc_POINT") mysize = 4
  }
     sf$levels <- addNA(a)
      levels(sf$levels)[levels(sf$levels)=="0"] <- "non-sig"
      levels(sf$levels)[levels(sf$levels)=="1"] <- "sig +"
      levels(sf$levels)[levels(sf$levels)=="2"] <- "sig -"
      cols <- c("non-sig" = "grey77", "sig +" = "red", "sig -" = "blue")
      lplot_runs <- ggplot(sf) +
        geom_sf(aes(fill = levels), color = "black", shape = 21, size = mysize) +
        theme_bw() +
        theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
        xlab(paste0("Significance p-value = ", sig)) +
        scale_fill_manual(values = cols, na.value ="grey",drop = FALSE)

      lplot_runs
}
