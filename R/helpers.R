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
#' @param object numeric vector of class `"MCS"` or `"MCS_log"`. The values
#' from n Monte Carlo Simulations
#' @param unit (optional) object of class `units` or `symbolic_units`, or in
#' the case of `set_units` expression with symbols.
#' @param ... additional arguments affecting the summary produced.
#'
#' @returns a list.
#' If class of `object` is `"MCS"`, the list contains the following elements:
#' \describe{
#' \item{`median`}{median of the Monte Carlo simulations}
#' \item{`ir.95`}{the 95% and 68% interpercentile range}
#' \item{`ir.68`}{the 68% interpercentile range}
#' \item{`mean`}{arithmetic mean the Monte Carlo simulations}
#' \item{`sd`}{1\eqn{\sigma} standard deviation}
#' \item{`conf.int`}{95% confidence intverval about the mean}
#' \item{`stderr`}{standard error}
#' \item{`t.test`}{Statistic and p-value of the Student's t-test}
#' \item{`n`}{Number of samples}
#' }
#'
#' If class of `object` is `"MCS_log"`, the list contains the following elements:
#' \describe{
#' \item{`median`}{median of the Monte Carlo simulations}
#' \item{`ir.95`}{the 95% and 68% interpercentile range}
#' \item{`ir.68`}{the 68% interpercentile range}
#' \item{`mean`}{geometric mean the Monte Carlo simulations}
#' \item{`sd`}{1\eqn{\sigma} range about the mean}
#' \item{`sd2`}{2\eqn{\sigma} range about the mean}
#' \item{`conf.int`}{95% confidence intverval about the mean}
#' \item{`stderr.log`}{standard error of `log(samples)`}
#' \item{`t.test`}{Statistic and p-value of the Student's t-test of `log(samples)`}
#' \item{`n`}{Number of samples}
#' }
#' Values will be in the unit specified by parameter `unit` or be equal to the
#' unit of `x` if `x` is a `units` object.
#'
#' @details Equations of the form \eqn{X = A b^{n \pm \sigma}} create
#' non-normal, left-skewed distributions (e.g. flow laws, and grain-size
#' piezometers).
#' Thus, it is recommended to report median and percentiles instead of mean,
#' standard deviation and confidence intervals.
#'
#' @importFrom stats median quantile t.test
#' @importFrom units drop_units
#'
#' @name summary-MCS
#'
#' @examples
#' set.seed(20250411)
#' MC_res <- grainsize_piezometry(12.2)
#' summary(MC_res)
#'
#' n <- 100
#' temperature <- units::set_units(rnorm(n, 300, 25), degC)
#' pressure <- units::set_units(rnorm(n, 400, 50), MPa)
#' MC_res2 <- ps_fugacity(pressure, temperature) # 37 MPa
#' summary(MC_res2)
NULL

#' @rdname summary-MCS
#' @exportS3Method base::summary
summary.MCS_log <- function(object, unit = NULL, ...) {
  if (inherits(object, "units")) {
    unit <- units(object)
    object <- units::drop_units(object)
  }
  x <- as.numeric(na.omit(object))


  median_s <- stats::median(x)
  ir_95 <- stats::quantile(x, c(0.025, 0.975))
  ir_68 <- stats::quantile(x, c(0.16, 0.84))

  log_s <- log10(x)

  log_s_tt <- t.test(log_s)

  mean_log_s <- unname(log_s_tt$estimate)
  CI95_log_s <- log_s_tt$conf.int
  stderr_log <- log_s_tt$stderr

  mean_s <- 10^mean_log_s

  sd_log_s <- sd(log_s)
  sd_range <- 10^c(mean_log_s - sd_log_s, mean_log_s + sd_log_s)
  sd2_range <- 10^c(mean_log_s - 2 * sd_log_s, mean_log_s + 2 * sd_log_s)

  se_log_s <- sd_log_s / sqrt(length(x))
  sde_range <- 10^c(mean_log_s - se_log_s, mean_log_s + se_log_s)
  conf.int <- 10^c(CI95_log_s) |>
    setNames(c("2.5%", "97.5%"))

  out <- list(
    # quantiles
    median = median_s |> set_units_if(unit),
    ir.68 = ir_68 |> set_units_if(unit),
    ir.95 = ir_95 |> set_units_if(unit),

    # log-normal stats
    mean = mean_s |> set_units_if(unit),
    sd.int = sd_range |> set_units_if(unit),
    sd2.int = sd2_range |> set_units_if(unit),
    sde.int = sde_range |> set_units_if(unit),
    conf.int = conf.int |> set_units_if(unit),
    stderr.log = stderr_log,
    t.test = c(log_s_tt$statistic, log_s_tt$p.value),
    n = length(object)
  )

  normal <- log_s_tt$p.value <= 0.05
  normal_text <- if (normal) {
    "p<0.05"
  } else {
    "not significant"
  }

  unit_text <- if(length(unit$numerator)>0 & length(unit$denominator)>0){
    paste(unit$numerator, '/', unit$denominator)
  } else if (length(unit$numerator)==0 & length(unit$denominator)>0){
    paste0('/', unit$denominator)
  } else {
    paste(unit$numerator, unit$denominator)
  }


  message(
    "Statistical summary of ", out$n, " Monte Carlo simulations\n\n",
    "Median:                      ", signif(out$median, 2), " ", unit_text, "\n",
    "95% interpercentile range:   ", signif(out$ir.95[1], 2), " - ", signif(out$ir.95[2], 2), " ", unit_text, "\n",
    "Standard error in log-space: ", signif(out$stderr.log), "\n",
    "Student's t-Test:            ", normal_text
  )

  return(
    invisible(out)
  )
}

#' @rdname summary-MCS
#' @exportS3Method base::summary
summary.MCS <- function(object, unit = NULL, ...) {
  if (inherits(object, "units")) {
    unit <- units(object)
    object <- units::drop_units(object)
  }
  x <- as.numeric(na.omit(object))

  median_s <- stats::median(x)
  ir_95 <- stats::quantile(x, c(0.025, 0.975))
  ir_68 <- stats::quantile(x, c(0.16, 0.84))


  s_tt <- t.test(x)

  mean_s <- unname(s_tt$estimate)
  CI95_s <- s_tt$conf.int
  stderr_s <- s_tt$stderr

  sd_s <- sd(x)
  # se_s <- sd_s / sqrt(length(x))

  out <- list(
    # quantiles
    median = median_s |> set_units_if(unit),
    ir.68 = ir_68 |> set_units_if(unit),
    ir.95 = ir_95 |> set_units_if(unit),

    mean = mean_s |> set_units_if(unit),
    sd = sd_s |> set_units_if(unit),
    # se = se_s |> set_units_if(unit),
    conf.int = CI95_s |> set_units_if(unit),
    stderr = stderr_s,
    t.test = c(s_tt$statistic, s_tt$p.value),
    n = length(object)
  )

  normal <- s_tt$p.value <= 0.05
  normal_text <- if (normal) {
    "p<0.05"
  } else {
    "not significant"
  }

  unit_text <- if(length(unit$numerator)>0 & length(unit$denominator)>0){
    paste(unit$numerator, '/', unit$denominator)
  } else if (length(unit$numerator)==0 & length(unit$denominator)>0){
    paste0('/', unit$denominator)
  } else {
    paste(unit$numerator, unit$denominator)
  }


  message(
    "Statistical summary of ", out$n, " Monte Carlo simulations\n\n",
    "Mean:                    ", signif(out$mean, 2), " ", unit_text, "\n",
    "95% confidence interval: ", signif(out$conf.int[1], 2), " - ", signif(out$conf.int[2], 2), " ", unit_text, "\n",
    "Standard error:          ", signif(out$stderr), "\n",
    "Student's t-Test:        ", normal_text
  )

  return(
    invisible(out)
  )
}


gas_const <- function() {
  units::set_units(8.31446261815324, J * K^-1 * mol^-1)
}

negative_strainrate <- function(x) {
  x[x < 0] <- 0
  return(x)
}
