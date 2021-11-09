#' @name summary.spqtest
#' @rdname summary.spqtest
#'
#' @title Summary of estimated objects of class \emph{spqtest}.
#'
#' @description  This function summarizes estimated \emph{spqtest} objects.
#'   The tables in the output include basic information for each test.
#'   blablabla...
#'
#' @param object An \emph{spqtest} object including a list of \emph{htest}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \emph{summary.spqtest}
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   }
#'
#' @seealso
#'   \code{\link{print.summary.spqtest}}
#' @examples
#'
#' @export
summary.spqtest <- function(object, ...) {
  z <- object
  stopifnot(inherits(z, "spqtest"))
  ## Build a tibble with the results...
  table_res <- NULL
  for (i in 1:length(z)) {
    table_res <- rbind(table_res, broom::tidy(z[[i]]))
  }
  table_res <- dplyr::rename(table_res, df = parameter)
  if (all(table_res$df == "NA")) table_res$df <- NULL
  table_res <- dplyr::rename(table_res, Q = statistic)
  table_res$method <- NULL
  for (i in 1:length(z)) {
    z[[i]]$Q <- z[[i]]$statistic <- NULL
    z[[i]]$df <- z[[i]]$parameter <- NULL
    z[[i]]$p.value <- z[[i]]$method <- NULL
    z[[i]]$data.name <- z[[i]]$symb <- z[[i]]$ms <- NULL
    z[[i]]$mdtms <- z[[i]]$efp_symb <- z[[i]]$efc_symb <- NULL
    z[[i]]$PSymb <- z[[i]]$CSymb <- NULL
    z[[i]]$qp_symb <- z[[i]]$qc_symb <- NULL
    z[[i]]$efp_symb <- z[[i]]$efc_symb <- NULL
    z[[i]]$qp_mc <- z[[i]]$qc_mc <- NULL
    z[[i]]$efp_symb_mc <- z[[i]]$efc_symb_mc <- NULL
    class(z[[i]]) <- c("list")
  }
  ltable <- t(list2DF(z))
  for (i in 1:ncol(ltable)) {
    newcol <- unlist(ltable[, i])
    oldcolnames <- colnames(table_res)
    table_res <- cbind(table_res, newcol)
    colnames(table_res) <- c(oldcolnames,
                             names(z[[1]])[i])
  }
  table_res <- tidyr::as_tibble(table_res)
  table_res$`R/n` <- table_res$R / table_res$n
  table_res$`5k^m` <- 5*(table_res$k^table_res$m)
  distance <- unique(table_res$distance)
  distribution <- unique(table_res$distr)
  if (distribution == "asymptotic") {
    tbl <- table_res %>%
      dplyr::select(var.name, type, Q, df, p.value,
             k, N, m, r, R,
             n, `R/n`, `5k^m`) %>%
      dplyr::group_by(var.name, type)
  } else if (distribution == "mc") {
    tbl <- table_res %>%
    dplyr::select(var.name, type, Q, p.value,
           k, N, m, R,
           n, `R/n`, `5k^m`) %>%
    dplyr::group_by(var.name, type)
  } else stop("distribution must be asymptotic or mc")
  gt_tbl <- gt::gt(tbl) %>%
    gt::tab_header(
      title = "Qualitative Dependence Test (Q)",
      subtitle = paste("Distribution: ",distribution,
                       ". Distance: ", distance,
                       sep = "") ) %>%
    gt::fmt_number(
      columns = c("Q","R/n"),
      decimals = 2) %>%
    gt::fmt_number(
      columns = c("p.value"),
      decimals = 5)
  class(gt_tbl) <- c("summary.spqtest", "gt_tbl")
  return(gt_tbl)
}


