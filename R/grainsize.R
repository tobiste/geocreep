#' Grain-Size Piezometry
#'
#' Calculates differential stress from grain size.
#' Uses Monte Carlo sampling for propagating parameter uncertainties in flow model.
#'
#' @param d numeric. Grain size in &mu;m or `units` object
#' @param sd (optional) Standard deviation of `d`
#' @param model character. One of
#' \describe{
#' \item{`"Stipp-reg2-3"`}{Piezometer for dislocation creep regime 2 and 3 (for deviatoric stress <368 MPa) after
#' Stipp and Tullis (2003)}
#' \item{`"Stripp-reg1"`}{Piezometer for dislocation creep regime 1 (deviatoric stress &ge;368 MPa) after
#' Stipp and Tullis (2003)}
#' \item{`"Cross-1"`}{Piezometer after Cross et al. (2017) for 1 &mu;m
#' step size resolution in EBSD data}
#' \item{`"Cross-sliding"`}{Sliding resolution piezometer after Cross et al. 2017.
#' According to authors, more accurately estimates stress in fine-grained (<10 &mu;m) samples}
#' }
#' @param sim non-negative integer. Number of Monte Carlo simulations
#'
#' @returns list. Differential stress in MPa. If Monte Carlo Simulation was used,
#' and object of class `"MCS"` is returned (see [summary()] for detailed description of output).
#' The piezometer produce log-normal distributed estimates considering the
#' uncertainties in the equation parameter. Hence it is recommended to report
#' the median (or geometric mean), and the interpercentile range.
#'
#' @details
#' General formula for grain-size piezometer is:
#' \deqn{\sigma = \left(\frac{d}{k}\right)^\frac{1}{n}}
#' where \eqn{\sigma} is the differential stress, \eqn{d} is the grain size, and
#' \eqn{k} and \eqn{n} are empirical parameters.
#'
#'
#' @references
#' Stipp, M., \& Tullis, J. (2003). The recrystallized grain size piezometer for
#' quartz. Geophysical Research Letters, 30(21), 1–5. https://doi.org/10.1029/2003GL018444
#'
#' Cross, A. J., Prior, D. J., Stipp, M., \& Kidder, S. (2017). The
#' recrystallized grain size piezometer for quartz: An EBSD-based calibration.
#' Geophysical Research Letters, 44(13), 6667–6674. https://doi.org/10.1002/2017GL073836
#'
#' @seealso [subgrainsize_piezometry()], [units::set_units()]
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' grainsize_piezometry(12.2) |> summary() # 92
#' grainsize_piezometry(31) |> summary() # 44
grainsize_piezometry <- function(d, sd = NULL, model = c("Stipp-reg2-3", "Stripp-reg1", "Cross-1", "Cross-sliding"), sim = 1e6) {
  # d in micrometre
  if (!is.null(sd)) d <- rnorm(sim, d, sd)

  d <- units::set_units(d, "um") |>
    as.numeric()

  model <- match.arg(model)

  if (model == "Stipp-reg2-3") {
    log_k <- 3.56
    log_k_ci95 <- 0.27
    n <- -1.26
    n_ci95 <- 0.13
  } else if (model == "Stripp-reg1") {
    log_k <- 1.89
    log_k_ci95 <- 0.11
    n <- -0.61
    n_ci95 <- 0.04
  } else if (model == "Cross-1") {
    log_k <- 3.91
    log_k_ci95 <- 0.41
    n <- -1.41
    n_ci95 <- 0.21
  } else if (model == "Cross-sliding") {
    log_k <- 4.22
    log_k_ci95 <- 0.51
    n <- -1.59
    n_ci95 <- 0.26
  }

  log_k_sd <- log_k_ci95 / 1.96 # convert 95% CI to SD
  n_sd <- n_ci95 / 1.96 # convert 95% CI to SD

  # k is the grain size exponent; n is the stress exponent
  k_samples <- 10^rnorm(sim, mean = log_k, sd = log_k_sd)
  n_samples <- rnorm(sim, mean = n, sd = n_sd)

  # solve D = k * sigma ^ n
  # stress_samples <- (k_samples / d)^(-1 / n_samples) # in MPa
  stress_samples <- (d / k_samples)^(1 / n_samples)

  stress <- units::set_units(stress_samples, "MPa")
  class(stress) <- append('MCS_log', class(stress))
  return(stress)
}

#' Subgrain‐Size Piezometer Calibrated for EBSD
#'
#' Calculates differential stress using the Subgrain‐Size Piezometer of Goddard et al. (2020).
#' Uses Monte Carlo sampling for propagating parameter uncertainties in flow model.
#'
#' @param lambda mean line intercept length in &mu;m or `units` object.
#' @param sd (optional) Standard deviation of `lambda`
#' @param calibrated logical. Whether the calibration of Holyoke and Kronenberg (2010) is considered or not.
#' @param min character. The Mineral uses. one of `"q"` for quartz, `"fo90"` for Olive with 90&#37; Forsterite, or `"fo50"` for Olivine with 50&#37; Forsterite.
#' @param sim non-negative integer. Number of Monte Carlo simulations
#'
#' @returns list. Differential stress in MPa. If Monte Carlo Simulation was used,
#' and object of class `"MCS"` is returned (see [summary()] for detailed description of output).
#' The piezometer produce log-normal distributed estimates considering the
#' uncertainties in the equation parameter. Hence it is recommended to report
#' the median (or geometric mean), and the interpercentile range.
#'
#' @details The sub-grain size piezometer is
#' \deqn{\frac{\lambda}{b} = 10^a \left(\frac{\sigma}{\mu}\right)^b}
#' where \eqn{\lambda} is the mean line intercept length, \eqn{b} is the Burgers vector,
#' \eqn{\sigma} is the differential stress, \eqn{\mu} is the shear modulus, and
#' \eqn{a} and \eqn{b} are the empirical exponents.
#'
#' @seealso [grainsize_piezometry()], [units::set_units()]
#'
#' @importFrom stats rnorm
#' @importFrom units set_units
#'
#' @references Goddard, R. M., Hansen, L. N., Wallis, D., Stipp, M., Holyoke,
#' C. W., Kumamoto, K. M., \& Kohlstedt, D. L. (2020). A Subgrain‐Size Piezometer
#' Calibrated for EBSD. Geophysical Research Letters, 47(23).
#' https://doi.org/10.1029/2020GL090056
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' subgrainsize_piezometry(9, min = "fo50") |> summary() # 420 MPa
#' subgrainsize_piezometry(18, min = "q") |> summary() # 240 MPa
subgrainsize_piezometry <- function(lambda, sd = NULL, calibrated = TRUE, min = c("q", "fo90", "fo50"), sim = 1e6) {
  min <- match.arg(min)

  if (!is.null(sd)) lambda <- stats::rnorm(sim, lambda, sd)

  lambda <- units::set_units(lambda, "um") |>
    units::set_units("m") |>
    as.numeric()

  b_q <- 5.1e-4
  b_ol <- 5.0e-4
  b <- c(q = b_q, "fo90" = b_ol, "fo50" = b_ol) |>
    units::set_units("um") |>
    units::set_units("m")
  burgers <- as.numeric(b[min])

  s <- c(q = 42.0, "fo90" = 77.8, "fo50" = 62.6) |>
    units::set_units("GPa") |>
    units::set_units("MPa")
  shear_m <- as.numeric(s[min])

  if (isTRUE(calibrated)) {
    a <- stats::rnorm(sim, 0.6, sd = 0.7 / 1.96)
    b <- stats::rnorm(sim, -1.2, sd = 0.3 / 1.96)
  } else {
    a <- stats::rnorm(sim, 1.2, sd = 1 / 1.96)
    b <- stats::rnorm(sim, -1.0, sd = 0.4 / 1.96)
  }

  stress <- shear_m * (lambda / (burgers * 10^a))^(1 / b)
  stress <- stress * 10 |> # to match the values as published in Goddard
    units::set_units("MPa")
  class(stress) <- append('MCS_log', class(stress))
  return(stress)
}
