na_przemian <- function(x, y) {
  temp <- c()
  for (i in seq_len(length(x))) {
    temp <- c(temp, x[i], y[i])
  }
  return(temp)
}

#' Print "multifanova" object
#'
#' Prints the summary of the global and multiple contrasts testing for functional data.
#'
#' @param object a "multifanova" object.
#' @param ... integer indicating the number of decimal places to be used to present the numerical results.
#' It can be named \code{digits} as in the \code{round()} function (see examples).
#'
#' @details The function prints out the information about the number of samples,
#' number of observations in each sample, number of design time points, contrasts used,
#' test statistics, critical values, \eqn{p}-values of tests performed by the
#' \code{multiFANOVA()} function. It also gives the decisions.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' # Some of the examples may run some time.
#'
#' # Canadian weather data set
#' # There are three samples of mean temperatures for
#' # fifteen weather stations in Eastern Canada,
#' # another fifteen in Western Canada, and
#' # the remaining five in Northern Canada.
#' library(fda)
#' data_set <- t(CanadianWeather$dailyAv[,, "Temperature.C"])
#' k <- 3
#' gr_label <- rep(c(1, 2, 3), c(15, 15, 5))
#' \donttest{
#' # Tukey's contrast matrix
#' h_tukey <- GFDmcv::contr_mat(k, type = "Tukey")
#' # testing without parallel computing
#' res <- multiFANOVA(data_set, gr_label, h_tukey)
#' summary(res, digits = 3)}
#' \dontshow{
#' h_tukey <- GFDmcv::contr_mat(k, type = "Tukey")
#' res <- multiFANOVA(data_set[, 1:10], gr_label, h_tukey, n_boot = 3)
#' summary(res, digits = 3)}
#'
#' @import fda
#'
#' @method summary multifanova
#' @export

# decisions based on p-values
summary.multifanova <- function(object, ...) {
  res_global_round <- round(object$res_global, ...)
  temp_global <- data.frame(statistic = c(res_global_round[1, 1], res_global_round[1, 3]),
                            p.value = c(res_global_round[1, 2], res_global_round[1, 4]),
                            decision = c(ifelse(object$res_global[1, 2] <= object$alpha, "H1", "H0"),
                                         ifelse(object$res_global[1, 4] <= object$alpha, "H1", "H0")))
  rownames(temp_global) <- c("GPH", "mGPH")
  res_multi_round <- round(object$res_multi, ...)
  kontrasty <- rep("", 2 * nrow(object$h))
  statystyki <- rep("", 2 * nrow(object$h))
  for (i in seq_len(2 * nrow(object$h))) {
    if (i %% 2 == 1) {
      # kontrasty[i] <- paste("(", paste(object$h[0.5 * i + 0.5, ], collapse = ", "), ")", sep = "")
      kontrasty[i] <- paste("Constrast", 0.5 * i + 0.5, sep = "_")
      statystyki[i] <- res_multi_round[0.5 * i + 0.5, 1]
    }
  }
  temp_multi <- data.frame(contrast = kontrasty,
                           statistic = statystyki,
                           test = rep(c("GPH", "mGPH"), times = nrow(object$h)),
                           critical_value = na_przemian(res_multi_round[, 2], res_multi_round[, 4]),
                           p.value = na_przemian(res_multi_round[, 3], res_multi_round[, 5]),
                           decision = ifelse(na_przemian(res_multi_round[, 3], res_multi_round[, 5]) <= object$alpha, "H1", "H0"))
  temp_contrasts <- as.data.frame(object$h,
                                  row.names = paste("Contrast", seq_len(nrow(object$h)), sep = "_"))
  colnames(temp_contrasts) <- paste("Group", seq_len(object$k), sep = "_")

  cat("#--- Multiple contrast tests for functional data ------------------------#", "\n", "\n")
  cat("- Number of samples:", object$k, "\n")
  cat("- Number of observations in samples:", object$n, "\n")
  cat("- Number of design time points:", object$j, "\n")
  cat("- Significance level:", object$alpha)
  cat("\n", "\n")
  cat("#--- Constrasts ---------------------------------------------------------#", "\n")
  print(temp_contrasts)
  cat("\n")
  cat("#--- Overall results ----------------------------------------------------#", "\n")
  print(temp_global)
  cat("\n")
  cat("#--- Multiple contrast testing results ----------------------------------#", "\n")
  print(temp_multi)
  cat("#------------------------------------------------------------------------#", "\n")
}
