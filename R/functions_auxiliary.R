
# Wald-type tests with all bootstrap methods, but integral only

# multiple contrast tests

# For critical values calculation from multi_critical_values_Marc_ls.R
# half partition search, quicker for higher degree of dependence
# This is crit_values2, but without unecessary calculations. So it is faster.
crit_values_ls <- function(data_mat, alpha) {
  # the input is a matrix data_mat of dimenson dim_obs x n_obs containing the B (resampling)
  # observations X_b of dimension dim_obs
  # for our specific purpose X_i is the vector of the test statistics based on the different
  # contrast vectors, i.e. dim_obs = r in the pdf
  # each column represent one observation (i.e. one resampling iteration)
  # the output are the critical values for each dimension (large X_b lead to a rejection)
  # such that the overall level is alpha
  n <- length(data_mat[1, ])
  dimension <- length(data_mat[, 1])

  # First forget the connection of the coordinates, and sort the values per coordinate
  data_order <- t(apply(data_mat, 1, sort))

  # worst case1
  #for each point outside the box only one coordinate leads to a rejection.
  # Thus, for each coordinate (alpha/dim) * number of obs are outside the box  (Bonferoni correction)
  j_low <- ceiling((alpha / dimension) * n) - 1
  # A <- data_mat / data_order[, n - j_low]
  # alpha_low <- mean(apply(A, 2, max) > 1) # count the points being outside the box (where the box borders are given by the critical values)

  # worst case1
  # something like totally linear dependence (in the dimension two)
  j_high <- ceiling(alpha * n)
  # A <- data_mat / data_order[, n - j_low] # ???
  # alpha_high <- mean(apply(A, 2, max) > 1) # count the points being outside the box (where the box borders are given by the critical values)

  # now we search for values j_low and j_high = j_low + 1, such that alpha_low <= alpha and alpha_high > alpha
  while (j_high - j_low > 1) {
    # approx. middle between j_low and j_high
    j_mid <- ceiling(j_low + (j_high - j_low) / 2)
    A <- data_mat / data_order[, n - j_mid]
    alpha_sim <- mean(apply(A, 2, max) > 1)
    ifelse(alpha_sim <= alpha, j_low <- j_mid, j_high <- j_mid)
  }
  # critical values
  return(data_order[, n - j_low])
}

# group mean and covariance calculation
# x - matrix of observations n times j (n = n_1 + ... + n_k)
# gr_label - a vector with group labels
mean_cov <- function(x, gr_label) {
  n <- nrow(x)
  j <- ncol(x)
  k <- length(unique(gr_label))
  n_vec <- numeric(k)
  for (i in seq_len(k)) n_vec[i] <- sum(gr_label == i)
  if (n != sum(n_vec)) stop("different number of observations in x and gr_label")

  # estimators of means and covariance functions in groups
  gr_means <- matrix(0, nrow = k, ncol = j)
  gr_cov <- vector("list", k)
  for (i in seq_len(k)) {
    x_i <- x[gr_label == i, ]
    mean_i <- colMeans(x_i)
    gr_means[i, ] <- mean_i
    z_i <- x_i - matrix(1, nrow = n_vec[i], ncol = 1) %*% mean_i
    gr_cov[[i]] <- t(z_i) %*% z_i / (n_vec[i] - 1)
  }

  return(list(n = n, j = j, k = k, n_vec = n_vec,
              gr_means = gr_means, gr_cov = gr_cov))
}

# globalizing pointwise Hotelling's T^2 -test (GPH) statistic
# H - contrast matrix
# mc - the result of mean_cov() function
gph_f <- function(H, mc) {
  # pointwise gph
  gph <- numeric(mc$j)
  H_gr_means <- H %*% mc$gr_means
  H_gr_means_t <- t(H_gr_means)
  H_t <- t(H)
  for (i in seq_len(mc$j)) {
    gamma_hat_t <- numeric(mc$k)
    for (i_k in seq_len(mc$k)) {
      gamma_hat_t[i_k] <- mc$gr_cov[[i_k]][i, i]
    }
    gph[i] <- mc$n * H_gr_means_t[i, ] %*% MASS::ginv(H %*% diag(mc$n * gamma_hat_t / mc$n_vec) %*% H_t) %*%
      H_gr_means[, i]
  }
  return(sum(gph))
}

# # ??? with c(t)
# gph_f <- function(H, mc, cc = 0) {
#   # pointwise gph
#   gph <- numeric(mc$j)
#   H_gr_means <- H %*% mc$gr_means - cc
#   H_gr_means_t <- t(H_gr_means)
#   H_t <- t(H)
#   for (i in seq_len(mc$j)) {
#     gamma_hat_t <- numeric(mc$k)
#     for (i_k in seq_len(mc$k)) {
#       gamma_hat_t[i_k] <- mc$gr_cov[[i_k]][i, i]
#     }
#     gph[i] <- mc$n * H_gr_means_t[i, ] %*% MASS::ginv(H %*% diag(mc$n * gamma_hat_t / mc$n_vec) %*% H_t) %*%
#       H_gr_means[, i]
#   }
#   return(sum(gph))
# }

#' Pointwise Hotelling's \eqn{T^2}-test statistic
#'
#' The function \code{ph_test_statistic()} calculates the pointwise Hotelling's \eqn{T^2}-test statistic.
#'
#' @param x matrix of observations \eqn{n\times j} (\eqn{n = n_1 + ... + n_k}, \eqn{j} is a number of design time points).
#' @param gr_label a vector with group labels; the integer labels (from 1 to a number of groups) should be used.
#' @param h contrast matrix. For Dunnett’s and Tukey’s contrasts, it can be created by
#' the \code{contr_mat()} function from the package \code{GFDmcv} (see examples).
#'
#' @details For details, see the documentation of the \code{multiFANOVA()} function or
#' the paper Munko et al. (2023).
#'
#' @return A vector of values of the pointwise Hotelling's \eqn{T^2}-test statistic.
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
#' # trajectories of mean temperatures
#' matplot(t(data_set), type = "l", col = gr_label, lty = 1,
#'         xlab = "Day", ylab = "Temperature (C)",
#'         main = "Canadian weather data set")
#' legend("bottom", legend = c("Eastern Canada", "Western Canada", "Northern Canada"),
#'        col = 1:3, lty = 1)
#' \donttest{
#' # Tukey's contrast matrix
#' h_tukey <- GFDmcv::contr_mat(k, type = "Tukey")
#' # testing without parallel computing
#' res <- multiFANOVA(data_set, gr_label, h_tukey)
#' summary(res, digits = 3)
#' # plots for pointwise Hotelling's T^2-test statistics
#' oldpar <- par(mfrow = c(2, 2), mar = c(2, 2, 2, 0.1))
#' plot(ph_test_statistic(data_set, gr_label, h_tukey), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Global hypothesis")
#' plot(ph_test_statistic(data_set, gr_label, matrix(h_tukey[1, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Contrast 1")
#' plot(ph_test_statistic(data_set, gr_label, matrix(h_tukey[2, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Contrast 2")
#' plot(ph_test_statistic(data_set, gr_label, matrix(h_tukey[3, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Contrast 3")
#' par(oldpar)}
#' \dontshow{
#' h_tukey <- GFDmcv::contr_mat(k, type = "Tukey")
#' res <- multiFANOVA(data_set[, 1:10], gr_label, h_tukey, n_boot = 3)
#' summary(res, digits = 3)
#' # plots for pointwise Hotelling's T^2-test statistics
#' oldpar <- par(mfrow = c(2, 2), mar = c(2, 2, 2, 0.1))
#' plot(ph_test_statistic(data_set, gr_label, h_tukey), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Global hypothesis")
#' plot(ph_test_statistic(data_set, gr_label, matrix(h_tukey[1, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Contrast 1")
#' plot(ph_test_statistic(data_set, gr_label, matrix(h_tukey[2, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Contrast 2")
#' plot(ph_test_statistic(data_set, gr_label, matrix(h_tukey[3, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, gr_label, h_tukey))),
#'      main = "Contrast 3")
#' par(oldpar)}
#'
#' @references Dunnett C. (1955) A multiple comparison procedure for comparing several treatments
#' with a control. Journal of the American Statistical Association 50, 1096-1121.
#'
#' Munko M., Ditzhaus M., Pauly M., Smaga L., Zhang J.T. (2023) General multiple tests for functional data. Preprint https://arxiv.org/abs/2306.15259
#'
#' Tukey J.W. (1953) The problem of multiple comparisons. Princeton University.
#'
#' @import MASS
#' @import fda
#'
#' @export
ph_test_statistic <- function(x, gr_label, h) {
  mc <- mean_cov(x, gr_label)
  # pointwise gph
  gph <- numeric(mc$j)
  H_gr_means <- h %*% mc$gr_means
  H_gr_means_t <- t(H_gr_means)
  H_t <- t(h)
  for (i in seq_len(mc$j)) {
    gamma_hat_t <- numeric(mc$k)
    for (i_k in seq_len(mc$k)) {
      gamma_hat_t[i_k] <- mc$gr_cov[[i_k]][i, i]
    }
    gph[i] <- mc$n * H_gr_means_t[i, ] %*% MASS::ginv(h %*% diag(mc$n * gamma_hat_t / mc$n_vec) %*% H_t) %*%
      H_gr_means[, i]
  }
  return(gph)
}

# # ??? with c(t)
# ph_test_statistic <- function(x, gr_label, h, cc = 0) {
#   mc <- mean_cov(x, gr_label)
#   # pointwise gph
#   gph <- numeric(mc$j)
#   H_gr_means <- h %*% mc$gr_means - cc
#   H_gr_means_t <- t(H_gr_means)
#   H_t <- t(h)
#   for (i in seq_len(mc$j)) {
#     gamma_hat_t <- numeric(mc$k)
#     for (i_k in seq_len(mc$k)) {
#       gamma_hat_t[i_k] <- mc$gr_cov[[i_k]][i, i]
#     }
#     gph[i] <- mc$n * H_gr_means_t[i, ] %*% MASS::ginv(h %*% diag(mc$n * gamma_hat_t / mc$n_vec) %*% H_t) %*%
#       H_gr_means[, i]
#   }
#   return(gph)
# }
