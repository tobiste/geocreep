# min_max_std <- function(a, b) {
#   x <- range(a, b)
#   std <- diff(x) / 2
#   mean <- mean(x)
#   c(mean, std)
# }

#' @importFrom units set_units as_units
set_units_if <- function(x, unit) {
  if (is.null(unit)) {
    x
  } else {
    units::set_units(x, units::as_units(unit), mode = "standard")
  }
}

#' Estimates from Monte Carlo Simulation
#'
#' @param x numeric vector. The values from n Monte Carlo Simulations
#' @param unit (optional) object of class `units` or `symbolic_units`, or in the case of `set_units` expression with symbols.
#'
#' @returns `'MC_sim'` object, i.e. a list.
#' \describe{
#' \item{`median`}{median of the Monte Carlo simulations}
#' \item{`mean`}{geometric mean the Monte Carlo simulations}
#' \item{`stderr_log`}{standard error of `log(samples)`}
#' \item{`ir_95`}{the 95% and 68% interpercentile range}
#' \item{`ir_68`}{the 68% interpercentile range}
#' \item{`samples`}{the Monte Carlo simulation}
#' }
#'
#' @details Equations of the form \eqn{X = A b^{n \pm \sigma}} create non-normal, left-skewed distributions.
#' Thus, it is recommended to report median and percentiles instead of mean, standard deviation and confidence intervals.
#'
#' @importFrom stats median quantile t.test
#'
#' @export
#'
#' @examples
#' mc_stats(rnorm(100), "Pa")
mc_stats <- function(x, unit = NULL) {
  median_s <- stats::median(x)
  ir_95 <- stats::quantile(x, c(0.025, 0.975))
  ir_68 <- stats::quantile(x, c(0.16, 0.84))

  log_s <- log10(x)

  log_s_tt <- t.test(log_s)

  mean_log_s <- log_s_tt$estimate
  # CI95_log_s <- log_s_tt$conf.int
  stderr_log <- log_s_tt$stderr

  mean_s <- 10^mean_log_s

  # sd_log_s   <- sd(log_s)
  # se_log_s   <- sd_log_s / sqrt(length(x))
  # sde_range <- 10^c(mean_log_s - se_log_s, mean_log_s + se_log_s)
  # conf.int <- 10^CI95_log_s

  out <- list(
    median = median_s |> set_units_if(unit),
    mean = mean_s |> set_units_if(unit),
    # sde = sde_range |> |> set_units_if(unit),
    ir_95 = ir_95 |> set_units_if(unit),
    ir_68 = ir_68 |> set_units_if(unit),
    # conf.int = conf.int |> |> set_units_if(unit),
    stderr_log = stderr_log,
    samples = x |> set_units_if(unit)
  )
  class(out) <- append(class(out), "MC_sim")

  return(out)
}
