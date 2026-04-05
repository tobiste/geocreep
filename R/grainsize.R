#' Grain size piezometry
#'
#' @param d numeric. Grain size in micrometer or units object
#' @param sd (optional) numeric. Standard deviation of `d`
#' @param method character. One of
#' \describe{
#' \item{Stipp-reg2-3"}{Piezometer for dislocation creep regime 2 and 3 after Stipp and Tullis (2003)}
#' \item{"Stripp-reg1"}{Piezometer for dislocation creep regime 1 after Stipp and Tullis (2003)}
#' \item{"Cross-1"}{Piezmeter after Cross et al. (2017) for 1 um step size resolution in EBSD data}
#' \item{"Cross-sliding"}{Sliding resolution piezometer after Cross et al. 2017. According to authors, more accurately estimates stress in fine-grained (<10 μm) samples}
#' }
#' @param sim non-negative integer. Number of Monte Carlo simulations
#'
#' @returns list. Stress in MPa given as median of the Monte Carlo simulations `median`, mean (`mean`), standard error range (`sde`), 95% and 68% interpercentile range (`ir_95` and `ir_68`), and the Monte Carlo simulation (`samples`).
#'
#' @references
#' Stipp, M., & Tullis, J. (2003). The recrystallized grain size piezometer for
#' quartz. Geophysical Research Letters, 30(21), 1–5. https://doi.org/10.1029/2003GL018444
#'
#' Cross, A. J., Prior, D. J., Stipp, M., & Kidder, S. (2017). The
#' recrystallized grain size piezometer for quartz: An EBSD-based calibration.
#' Geophysical Research Letters, 44(13), 6667–6674. https://doi.org/10.1002/2017GL073836
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' grainsize_piezometry(12.2) # 92
#' grainsize_piezometry(31) # 44
grainsize_piezometry <- function(d, sd = NULL, method = c("Stipp-reg2-3", "Stripp-reg1", "Cross-1", "Cross-sliding"), sim = 1e6) {
  # d in micrometre
  if(!is.null(sd)) d <- rnorm(sim, d, sd)

  d <- units::set_units(d, um) |> as.numeric()

  method <- match.arg(method)

  if (method == "Stipp-reg2-3") {
    log_k <- 3.56
    log_k_sd <- 0.27
    n <- -1.26
    n_sd <- 0.13
  } else if (method == "Stripp-reg1") {
    log_k <- 1.89
    log_k_sd = 0.11
    n <- -0.61
    n_std = 0.04
  } else if (method == "Cross-1") {
    log_k <- 3.91
    log_k_sd <- 0.41
    n <- -1.41
    n_sd <- 0.21
  } else if (method == "Cross-sliding") {
    log_k <- 4.22
    log_k_std <- 0.51
    n <- -1.59
    n_sd <- 0.26
  }

  # k is the grain size exponent; n is the stress exponent
  k_samples <- 10^rnorm(sim, mean = log_k, sd = log_k_sd)
  n_samples <- rnorm(sim, mean = n, sd = n_sd)

  # solve D = k * sigma ^ n
  #stress_samples <- (k_samples / d)^(-1 / n_samples) # in MPa
  stress_samples <- (d / k_samples)^(1 / n_samples)

  mc_stats(stress_samples, "MPa")
}
