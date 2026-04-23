#' @noRd
#' @importFrom stats rnorm runif
#' @importFrom units set_units
flow_models <- function(model) {
  ci2sd <- 1 / 1.96 # convert 95% CI to SD

  #x <- flow_model_params(model)

  list(
    Kronenberg1984 = function(stress, temperature, fugacity = NULL, grainsize, pressure = NULL, sim, propagate_err = TRUE) {
      H_min <- 120
      H_max <- 150
      n_min <- 2.9
      n_max <- 3.2
      m <- 0.18

      if (isTRUE(propagate_err)) {
        n <- runif(sim, n_min, n_max)
        H <- runif(sim, H_min, H_max) |> set_units("kJ mol-1")
      } else {
        n <- mean(c(n_min, n_max))
        H <- mean(c(H_min, H_max)) |> set_units("kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      stress^n * grainsize^m * arrhenius
    },
    Paterson1990 = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim, propagate_err = FALSE) {
      A <- 6.5e-8
      H <- 135 |> set_units("kJ mol-1")
      n <- 3.1

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * arrhenius
    },
    Luan1992 = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s
      A <- 4e-10
      H <- 152
      H_std <- 71
      n <- 4
      n_std <- 0.8
      # A * stress^n * grainsize^m * fugacity^r * exp(-H / (R * temperature))

      if (isTRUE(propagate_err)) {
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * arrhenius
    },
    Gleason1995 = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s
      A_factor <- 1.1
      A_mean <- -4
      A_std <- 2
      H <- 223
      H_std <- 56
      n <- 4
      n_std <- 0.9

      if (isTRUE(propagate_err)) {
        A <- A_factor * 10^rnorm(sim, mean = A_mean, sd = A_std)
        H <- rnorm(sim, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        A <- A_factor * (10^A_mean)
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * arrhenius
    },
    Gleason1995_melt = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s
      A_factor <- 1.8
      A_mean <- -8
      A_std <- 2
      H <- 137
      H_std <- 34
      n <- 4
      n_std <- 0.9

      if (isTRUE(propagate_err)) {
        A <- A_factor * 10^rnorm(sim, mean = A_mean, sd = A_std)
        H <- rnorm(sim,mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        A <- A_factor * (10^A_mean)
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * arrhenius
    },
    Rutter2004 = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      log_A <- -4.93
      log_A_std <- 0.34
      H <- 242
      H_std <- 24 # 1s
      n <- 2.97
      n_std <- 0.29 # 1s
      r <- 1

      if (isTRUE(propagate_err)) {
        A <- 10^rnorm(sim, log_A, log_A_std)
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        A <- 10^log_A
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    },
    Fukuda2018_LT = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s
      log_A <- -2.97
      log_A_std <- 0.23
      H <- 129
      H_std <- 33 # std given as 1s?
      n_min <- 2.9
      n_max <- 5.2
      r <- 1

      if (isTRUE(propagate_err)) {
        A <- 10^rnorm(sim, log_A, log_A_std)
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- runif(sim, n_min, n_max)
      } else {
        A <- 10^log_A
        H <- set_units(H, "kJ mol-1")
        n <- mean(c(n_min, n_max))
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    },
    Fukuda2018_HT = function(stress, temperature, fugacity, grainsize, pressure = NULL, sim, propagate_err = TRUE) {
      log_A <- -2.97
      log_A_std <- 0.23
      H <- 183.0
      H_std <- 25 # std given as 1s?
      n <- 1.7
      n_std <- 0.2 # std given as 1s?
      m <- -0.51
      m_std <- 0.13 # std given as 1s?
      r <- 1.0
      r_std <- 0.2 # std given as 1s?

      if (isTRUE(propagate_err)) {
        A <- 10^rnorm(sim, log_A, log_A_std)
        H <- rnorm(sim, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
        m <- rnorm(sim, m, m_std)
        r <- rnorm(sim, r, r_std)
      } else {
        A <- 10^log_A
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * grainsize^m * fugacity^r * arrhenius
    },
    Richter2018 = function(stress, temperature, fugacity = NULL, grainsize, pressure = NULL, sim, propagate_err = TRUE) {
      A <- 3.1e-4
      # H_min <- 168
      # H_max <- 170
      H <- 170
      H_std <- 72
      n <- 1.9
      n_std <- 0.6
      m <- 1.08

      if (isTRUE(propagate_err)) {
        # std given as 1s?
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * grainsize^m * arrhenius
    },
    Hirth2001 = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s

      log_A <- -11.2
      log_A_std <- 0.6 # prefactor in MPa^{-n} / s
      H <- 135 # Activation enthalpy in kJ/mol
      H_std <- 15
      n <- 4 # stress exponent; adopted from Gleason and Tullis (1995) and Luan and Paterson (1992)
      n_std <- 0.85
      r <- 1 # water fugacity exponent

      if (isTRUE(propagate_err)) {
        A <- 10^rnorm(sim, log_A, log_A_std)
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        A <- 10^log_A
        H <- set_units(H, "kJ mol-1")
      }
      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    },
    Lu2019 = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s
      A <- 6
      A_std <- 5
      A_k <- -15
      H <- 132
      H_std <- 5
      # V <- 35.3 |> set_units("cm3 mol-1") # cm3/mol
      n <- 4 # adopted from Gleason and Tullis (1995) and Luan and Paterson (1992)
      n_std <- 0.85 # mean std from Gleason and Tullis (1995) and Luan and Paterson (1992)
      r <- 2.7

      if (isTRUE(propagate_err)) {
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        A <- truncnorm::rtruncnorm(sim, a = 0, mean = A, sd = A_std) * 10^(A_k)
        n <- rnorm(sim, n, n_std)
      } else {
        A <- A* 10^(A_k)
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    },
    Tokle2019_HT = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s
      A <- 8e-12
      H <- 140
      H_std <- 15
      n <- 4
      n_std <- 0.3
      r <- 1

      if (isTRUE(propagate_err)) {
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    },
    Tokle2019_LT = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim, propagate_err = TRUE) {
      # uncertainties not specified; assuming 1s
      A <- 5.4e-12
      H <- 105
      H_std <- 15
      n <- 2.7
      n_std <- 0.3
      r <- 1.1

      if (isTRUE(propagate_err)) {
        H <- truncnorm::rtruncnorm(sim, a = 0, mean = H, sd = H_std) |> set_units("kJ mol-1")
        n <- rnorm(sim, n, n_std)
      } else {
        H <- set_units(H, "kJ mol-1")
      }

      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    },
    Lusk2021_LP = function(stress, temperature, fugacity, grainsize = NULL, pressure, sim, propagate_err = TRUE) {
      # "low-pressure" <560 MPa

      # std given as 1sd
      log_A <- -9.3 # MPa^(-n-r) s^-1
      log_A_std <- 0.66
      n <- 3.5
      n_std <- 0.2
      r <- 0.49
      r_std <- 0.13
      Q <- 118
      Q_std <- 5 # kJ mol-1
      V <- 2.59
      V_std <- 2.45 # cm3 mol-1

      if (isTRUE(propagate_err)) {
        A <- 10^rnorm(sim, log_A, log_A_std)
        n <- rnorm(sim, n, n_std)
        r <- rnorm(sim, r, r_std)
        Q <- truncnorm::rtruncnorm(sim, a = 0, mean = Q, sd = Q_std) |> set_units("kJ mol-1")
        V <- truncnorm::rtruncnorm(sim, a = 0, mean = V, sd = V_std) |> set_units("cm3 mol-1")
      } else {
        A <- 10^log_A
        Q <- set_units(Q, "kJ mol-1")
        V <- set_units(V, "cm3 mol-1")
      }
      H <- Q + V * pressure
      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    },
    Lusk2021_HP = function(stress, temperature, fugacity, grainsize = NULL, pressure, sim, propagate_err = TRUE) {
      # “high-pressure” 700–1600 MPa

      # std given as 1sd
      log_A <- -7.90
      log_A_std <- 0.34 # MPa−n−r s−1;
      n <- 2.0
      n_std <- 0.1
      r <- 0.49
      r_std <- 0.13
      Q <- 77
      Q_std <- 8 # kJ mol−1;
      V <- 2.59
      V_std <- 2.45 # cm3 mol−1

      if (isTRUE(propagate_err)) {
        A <- 10^rnorm(sim, log_A, log_A_std)
        n <- rnorm(sim, n, n_std)
        r <- rnorm(sim, r, r_std)
        Q <- truncnorm::rtruncnorm(sim, a = 0, mean = Q, sd = Q_std) |> set_units("kJ mol-1")
        V <- truncnorm::rtruncnorm(sim, a = 0, mean = V, sd = V_std) |> set_units("cm3 mol-1")
      } else {
        A <- 10^log_A
        Q <- set_units(Q, "kJ mol-1")
        V <- set_units(V, "cm3 mol-1")
      }

      H <- Q + V * pressure
      RT <- gas_const() * temperature

      term <- -H / RT
      stopifnot(units(term) == units::unitless)
      arrhenius <- exp(as.numeric(term))

      A * stress^n * fugacity^r * arrhenius
    }
  )
}


#' Strain Rates for Creep in Quartz
#'
#' Calculates strain rates of deforming quartz from stress, temperature, and grain size from
#' experimentally determined creep law parameters.
#' Monte Carlo simulation is used for propagating parameter uncertainties in to creep estimate.
#'
#' @param stress Differential stress in MPa or `units` object
#' @param temperature Temperature in Kelvin or `units` object
#' @param fugacity Water fugacity in MPa or `units` object
#' @param grainsize Grain size in \eqn{\mu}m or `units` object
#' @param pressure Pressure in MPa or `units` object
#' @param sim non-negative number. Number of Monte Carlo simulations
#' @param propagate_err logical. Whether errors of the flow law parameters
#' should be propagated. `TRUE` by default.
#' @param model character specifying the flow law to be used:
#'  \describe{
#' \item{`"Hirth2001"`}{Hirth and Tullis (2001), dislocation creep}
#' \item{`"Paterson1990"`}{Paterson and Luan (1990): dislocation creep; axial compression}
#' \item{`"Kronenberg1984"`}{Kronenberg and Tullis (1984): deformation mechanism: dislocation creep and grain-size sensitive creep; strain geometry: axial compression}
#' \item{`"Luan1992"`}{Luan and Paterson (1990): dislocation creep; axial compression}
#' \item{`"Gleason1995"`}{Gleason and Tullis (1995): dislocation creep; axial compression}
#' \item{`"Gleason1995_melt"`}{Gleason and Tullis (1995): dislocation creep; axial compression; 1&ndash;2&#37; melt}
#' \item{`"Rutter2004"`}{Rutter and Brodie (2004b): dislocation creep; axial compression}
#' \item{`"Fukuda2018_LT"`}{Fukada et al. (2018): dislocation creep; axial compression; low temperatures (600&ndash;750 &deg;C)}
#' \item{`"Fukuda2018_HT"`}{Fukada et al. (2018): dislocation creep and grain-size sensitive creep; axial compression; high temperatures (800&ndash;950 &deg;C)}
#' \item{`"Richter2018"`}{Richter et al. (2018): dislocation creep and grain-size sensitive creep; general shear; 800&ndash;1000 &deg;C}
#' \item{`"Lu2019"`}{Lu and Jiang (2019): dislocation creep}
#' \item{`"Tokle2019_HT"`}{Tokle et al. (2019): dislocation creep and grain-size sensitive creep; high temperatures/low stress}
#' \item{`"Tokle2019_LT"`}{Tokle et al. (2019): dislocation creep and grain-size sensitive creep; low temperature/high stress}
#' \item{`"Lusk2021_LP"`}{Lusk et al. (2021): dislocation-dominated creep in wet quartz, for low pressures (&le;560 MPa)}
#' \item{`"Lusk2021_HP"`}{Lusk et al. (2021): dislocation-dominated creep in wet quartz, for high pressures (700&ndash;1600 MPa)}
#' }
#'
#' @details General flow law giving the strain rate is  \deqn{\dot{\varepsilon} = A \sigma^n d^m f_{H_2O}^r \, e^{\left({\frac{-H}{RT}}\right)}}
#'
#' where \eqn{\sigma} is the differential stress, \eqn{d} is the grain size,
#' \eqn{f_{H_2O}} is the water fugacity, \eqn{T} is the temperature,
#' \eqn{H} is the enthalpy, and \eqn{R} is the ideal gas constant. The flow
#' parameters are the prefactor \eqn{A}, and the exponents \eqn{n}, \eqn{m}, and \eqn{r}.
#'
#' To propagate the uncertainties of the flow parameters Monte Carlo simulation is used here.
#' \itemize{
#' \item{If the flow law parameters are given by a mean value and a marginal error (\eqn{\mu \pm z}),
#' the Monte Carlo simulation assumes a normal distribution given by
#' \eqn{X = N\left(\mu, \sigma\right)}, where \eqn{\mu} is the mean and \eqn{\sigma}
#'  is the standard deviation of the mean (\eqn{\sigma = \text{z}/1.96}).}
#' \item{If the parameter is given by a range of possible values \eqn{\left[x_\text{min}, x_\text{max}\right]},
#'  the Monte Carlo simulation assumes an uniform distribution given by \eqn{X = U\left(x_\text{min}, x_\text{max}\right)}.}
#' }
#'
#' @returns list. Strain rate in 1/s. If Monte Carlo Simulation was used, and
#' object of class `"MCS"` is returned (see [summary()] for detailed description of output).
#' The flow laws produce log-normal distributed estimates considering the
#' uncertainties in the parameter. Hence it is recommended to report the median
#' (or geometric mean), and the interpercentile range.
#'
#' @references
#' Fukuda, J., Holyoke, C. W., & Kronenberg, A. K. (2018). Deformation of Fine‐Grained Quartz Aggregates by Mixed Diffusion and Dislocation Creep. Journal of Geophysical Research: Solid Earth, 123(6), 4676-4696. \doi{10.1029/2017JB015133}
#'
#' Gleason, G. C., & Tullis, J. (1995). A flow law for dislocation creep of quartz aggregates determined with the molten salt cell. Tectonophysics, 247(1-4), 1-23. \doi{10.1016/0040-1951(95)00011-B}
#'
#' Hirth, G., Teyssier, C., & Dunlap, W. J. (2001). An evaluation of quartzite flow laws based on comparisons between experimentally and naturally deformed rocks. International Journal of Earth Sciences, 90(1), 77-87. \doi{10.1007/s005310000152}
#'
#' Kronenberg, A. K., & Tullis, J. (1984). Flow strengths of quartz aggregates: Grain size and pressure effects due to hydrolytic weakening. Journal of Geophysical Research: Solid Earth, 89(B6), 4281–4297. \doi{10.1029/JB089iB06p04281}
#'
#' Lu, L. X., & Jiang, D. (2019). Quartz Flow Law Revisited: The Significance of Pressure Dependence of the Activation Enthalpy. Journal of Geophysical Research: Solid Earth, 124(1), 241–256. \doi{10.1029/2018JB016226}
#'
#' Luan, F. C., & Paterson, M. S. (1992). Preparation and Deformation of Synthetic Aggregates of Quartz. Journal of Geophysical Research, 97, 301–320. https://doi.org/10.1029/91JB01748
#'
#' Lusk, A. D. J., Platt, J. P., & Platt, J. A. (2021). Natural and Experimental Constraints on a Flow Law for Dislocation‐Dominated Creep in Wet Quartz. Journal of Geophysical Research: Solid Earth, 126(5), 1-25. \doi{10.1029/2020JB021302}
#'
#' Paterson, M. S., & Luan, F. C. (1990). Quartzite rheology under geological conditions. Geological Society, London, Special Publications, 54(1), 299–307. \doi{10.1144/GSL.SP.1990.054.01.26}
#'
#' Richter, B., Stünitz, H., & Heilbronner, R. (2018). The brittle-to-viscous transition in polycrystalline quartz: An experimental study. Journal of Structural Geology, 114(September 2017), 1-21. \doi{10.1016/j.jsg.2018.06.005}
#'
#' Rutter, E. H., & Brodie, K. H. (2004a). Experimental grain size-sensitive flow of hot-pressed Brazilian quartz aggregates. Journal of Structural Geology, 26(11), 2011–2023. \doi{10.1016/j.jsg.2004.04.006}
#'
#' Rutter, E. ., & Brodie, K. . (2004b). Experimental intracrystalline plastic flow in hot-pressed synthetic quartzite prepared from Brazilian quartz crystals. Journal of Structural Geology, 26(2), 259–270. \doi{10.1016/S0191-8141(03)00096-8}
#'
#' Tokle, L., Hirth, G., & Behr, W. M. (2019). Flow laws and fabric transitions in wet quartzite. Earth and Planetary Science Letters, 505, 152-161. \doi{10.1016/j.epsl.2018.10.017}
#' @export
#'
#' @seealso [units::set_units()] to set up `units` objects; [summary.MCS_log()] for statistical parameters of Monte Carlo samples;
#' [creep_quartz_analytic()] for an analytical solution
#'
#' @importFrom stats rnorm
#' @importFrom units set_units unitless
#' @importFrom truncnorm rtruncnorm
#'
#' @examples
#' set.seed(20250411)
#' stress <- units::set_units(100, MPa)
#' temperature <- units::set_units(300, degC)
#' pressure <- units::set_units(400, MPa)
#' fugacity <- ps_fugacity(pressure, temperature)
#'
#' creep_quartz(
#'   stress = stress, temperature = temperature, fugacity = fugacity,
#'   model = "Hirth2001") |>
#'  summary()
creep_quartz <- function(stress, temperature, fugacity = NULL,
                         grainsize = NULL, pressure = NULL,
                         model = c("Hirth2001", "Paterson1990", "Kronenberg1984", "Luan1992", "Gleason1995", "Gleason1995_melt", "Rutter2004", "Fukuda2018_LT", "Fukuda2018_HT", "Richter2018", "Lu2019", "Tokle2019_LT", "Tokle2019_HT", "Lusk2021_LP", "Lusk2021_HP"),
                         propagate_err = TRUE,
                         sim = 1e6) {
  model <- match.arg(model)

  # stress, pressure and fugacity in Mega-Pascal, temperature in Kelvins
  temperature <- units::set_units(temperature, "K") #|> as.numeric()
  stress <- units::set_units(stress, "MPa") |> as.numeric()

  args <- list(stress = stress, temperature = temperature, sim = sim, propagate_err = propagate_err)

  if (!is.null(fugacity)) {
    args$fugacity <- units::set_units(fugacity, "MPa") |> as.numeric()
  } else {
    args$fugacity <- 1
  }

  if (!is.null(pressure)) args$pressure <- units::set_units(pressure, "MPa") #|> as.numeric()

  if (!is.null(grainsize)) {
    args$grainsize <- units::set_units(grainsize, "um") |> as.numeric()
  } else {
    args$grainsize <- 1
  }

  # Validate model name
  fm <- flow_models()
  if (!model %in% names(fm)) {
    stop(paste("Unknown model '", model, "'. Available models: ",
      paste(names(fm), collapse = ", "),
      sep = ""
    ))
  }

  # Load model parameters and function
  fm_FUN <- fm[[model]]

  edot <- do.call(fm_FUN, args = args) |>
    negative_strainrate() |>
    units::set_units("1/s")

  if (length(edot) > 1) {
    class(edot) <- append("MCS_log", class(edot))
  }
  return(edot)
}
