#' @noRd
#' @importFrom stats rnorm runif
#' @importFrom units set_units
flow_models <- function() {
  R <- gas_const() |>
    # set_units(MPa*cm3*K^-1*mol^-1) |>
    units::set_units(kJ * K^-1 * mol^-1) |>
    as.numeric()

  ci2sd <- 1 / 1.96 # convert 95% CI to SD

  list(
    Kronenberg1984 = function(stress, temperature, fugacity = NULL, grainsize, pressure = NULL, sim) {
      H_min <- 120
      H_max <- 150
      n_min <- 2.9
      n_max <- 3.2
      m <- 0.18
      n <- runif(sim, n_min, n_max)
      H <- runif(sim, H_min, H_max)

      stress^n * grainsize^m * exp(-(H) / (R * temperature))
    },
    Paterson1990 = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim) {
      A <- 6.5e-8
      H <- 135
      n <- 3

      A * stress^n * exp(-H / (R * temperature))
    },
    Luan1992 = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim) {
      A <- 4e-10
      H <- 152
      H_std <- 71
      n <- 4
      n_std <- 0.8
      # A * stress^n * grainsize^m * fugacity^r * exp(-H / (R * temperature))
      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)

      A * stress^n * exp(-H / (R * temperature))
    },
    Gleason1995 = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim) {
      A_factor <- 1.1
      A_mean <- -4
      A_sd <- 2
      H <- 223
      H_std <- 56
      n <- 4
      n_std <- 0.9

      A <- A_factor * 10^rnorm(sim, A_mean, A_sd * ci2sd)
      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)

      A * stress^n * exp(-H / (R * temperature))
    },
    Gleason1995_melt = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim) {
      A_factor <- 1.8
      A_mean <- -8
      A_sd <- 2
      H <- 137
      H_std <- 34
      n <- 4
      n_std <- 0.9

      A <- A_factor * 10^rnorm(sim, A_mean, A_sd * ci2sd)
      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)

      A * stress^n * exp(-H / (R * temperature))
    },
    Rutter2004 = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim) {
      A <- 1.2e-5
      H <- 242
      H_std <- 24
      n <- 2.97
      n_std <- 0.29
      r <- 1

      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)

      A * stress^n * fugacity^r * exp(-H / (R * temperature))
    },
    Fukuda2018_LT = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim) {
      H <- 129
      H_std <- 33
      n_min <- 2.9
      n_max <- 5.2
      r <- 1
      H <- rnorm(sim, H, H_std * ci2sd)
      n <- runif(sim, n_min, n_max)

      stress^n * fugacity^r * exp(-H / (R * temperature))
    },
    Fukuda2018_HT = function(stress, temperature, fugacity, grainsize, pressure = NULL, sim) {
      A <- -2.97
      A_std <- 0.23
      H <- 183.0
      H_std <- 25
      n <- 1.7
      n_std <- 0.2
      m <- -0.51
      m_std <- 0.13
      r <- 1.0
      r_std <- 0.2

      A <- 10^rnorm(sim, A, A_std * ci2sd)
      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)
      m <- rnorm(sim, m, m_std * ci2sd)
      r <- rnorm(sim, r, r_std * ci2sd)

      A * stress^n * grainsize^m * fugacity^r * exp(-H / (R * temperature))
    },
    Richter2018 = function(stress, temperature, fugacity = NULL, grainsize = NULL, pressure = NULL, sim) {
      A <- 3.1e-4
      # H_min <- 168
      # H_max <- 170
      H <- 170
      H_std <- 72
      n <- 1.9
      n_std <- 0.6
      m <- 1.08

      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)

      A * stress^n * grainsize^m * exp(-H / (R * temperature))
    },
    Hirth2001 = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim) {
      log_A <- -11.2
      log_A_std <- 0.6 # prefactor in MPa^{-n} / s
      H <- 135 # Activation enthalpy in kJ/mol
      H_std <- 15
      n <- 4 # stress exponent
      r <- 1 # water fugacity exponent
      A <- 10^rnorm(sim, log_A, log_A_std * ci2sd)
      H <- rnorm(sim, H, H_std * ci2sd)

      A * stress^n * fugacity^r * exp(-H / (R * temperature))
    },
    Lu2019 = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim) {
      A <- 6e-15
      A_std <- 5e-15
      H <- 132
      H_std <- 5
      V <- 35.3 # cm3/mol
      n <- 4 # adopted from Gleason and Tullis (1995) and Luan and Paterson (1992)
      n_std <- 0.85 # mean std from Gleason and Tullis (1995) and Luan and Paterson (1992)
      r <- 2.7

      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)
      A <- rnorm(sim, A, A_std * ci2sd)

      A * stress^n * fugacity^r * exp(-H / (R * temperature))
    },
    Tokle2019_HT = function(stress, temperature, fugacity, grainsize = NULL, pressure = NULL, sim) {
      A <- 8e-12
      H <- 140
      H_std <- 15
      n <- 4
      n_std <- 0.3
      r <- 1

      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)

      A * stress^n * fugacity^r * exp(-H / (R * temperature))
    },
    Tokle2019_LT = function(stress, temperature, fugacity, sim) {
      A <- 5.4e-12
      H <- 105
      H_std <- 15
      n <- 2.7
      n_std <- 0.3
      r <- 1.1

      H <- rnorm(sim, H, H_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)

      A * stress^n * fugacity^r * exp(-H / (R * temperature))
    },
    Lusk2021 = function(stress, temperature, fugacity, grainsize = NULL, pressure, sim) {
      # "low-pressure" <560 MPa
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

      A <- 10^rnorm(sim, log_A, log_A_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)
      r <- rnorm(sim, r, r_std * ci2sd)
      Q <- rnorm(sim, Q, Q_std * ci2sd)
      V <- rnorm(sim, V, V_std * ci2sd)
      H <- Q + V * pressure

      A * stress^n * grainsize^m * fugacity^r * exp(-H / (R * temperature))
    },
    Lusk2021_HP = function(stress, temperature, fugacity, grainsize = NULL, pressure, sim) {
      # “high-pressure” 700–1600 MPa
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

      A <- 10^rnorm(sim, log_A, log_A_std * ci2sd)
      n <- rnorm(sim, n, n_std * ci2sd)
      r <- rnorm(sim, r, r_std * ci2sd)
      Q <- rnorm(sim, Q, Q_std * ci2sd)
      V <- rnorm(sim, V, V_std * ci2sd)
      H <- Q + V * pressure

      A * stress^n * grainsize^m * fugacity^r * exp(-H / (R * temperature))
    }
  )
}


#' Strain Rates for Creep in Quartz
#'
#' Calculates strain rates of deforming quartz from stress, temperature, and grain size from
#' experimentally determined creep law parameters.
#' Monte Carlo sampling is used for propagating parameter uncertainties in to creep estimate.
#'
#' @param stress Differential stress in MPa or `units` object
#' @param temperature Temperature in Kelvin or `units` object
#' @param fugacity Water fugacity in MPa or `units` object
#' @param grainsize Grain size in \eqn{\mu}m or `units` object
#' @param pressure Pressure in MPa or `units` object
#' @param sim non-negative number. Number of Monte Carlo simulations
#' @param model character specifying the flow law to be used:
#'  \describe{
#' \item{`"Hirth2001"`}{Hirth and Tullis (2001), dislocation creep}
#' \item{`"Paterson1990"`}{Paterson and Luan (1990): dislocation creep; axial compression}
#' \item{`"Kronenberg1984"`}{Kronenberg and Tullis (1984): deformation mechanism: dislocation creep and grain-size sensitive creep; strain geometry: axial compression}
#' \item{`"Luan1992"`}{Luan and Paterson (1990): dislocation creep; axial compression}
#' \item{`"Gleason1995"`}{Gleason and Tullis (1995): dislocation creep; axial compression}
#' \item{`"Gleason1995_melt"`}{Gleason and Tullis (1995): dislocation creep; axial compression; 1-2\% melt}
#' \item{`"Rutter2004"`}{Rutter and Brodie (2004): dislocation creep; axial compression}
#' \item{`"Fukuda2018_LT"`}{Fukada et al. (2018): dislocation creep; axial compression}
#' \item{`"Fukuda2018_HT"`}{Fukada et al. (2018): dislocation creep and grain-size sensitive creep; axial compression; HT-fit}
#' \item{`"Richter2018"`}{Richter et al. (2018): dislocation creep and grain-size sensitive creep; general shear; between 800 and 1000 degree Celsius}
#' \item{`"Lu2019"`}{Lu and Jiang (2019): dislocation creep}
#' \item{`"Tokle2019_HT"`}{Tokle et al. (2019): dislocation creep and grain-size sensitive creep; high temperatures/low stress}
#' \item{`"Tokle2019_LT"`}{Tokle et al. (2019): dislocation creep and grain-size sensitive creep; low temperature/high stress}
#' \item{`"Lusk2021"`}{Lusk et al. (2021): dislocation-dominated creep in wet quartz, for low pressures (less than 560 MPa)}
#' \item{`"Lusk2021_HP"`}{Lusk et al. (2021): dislocation-dominated creep in wet quartz, for high pressures (700-1600 MPa)}
#' }
#'
#' @details General flow law giving the strain rate is  \deqn{\dot{\epsilon} = A \sigma^n d^m f_{H_2O}^r \, e^{\left({\frac{-H}{RT}}\right)}}
#'
#' where \eqn{\sigma} is the differential stress, \eqn{d} is the grain size,
#' \eqn{f_{H_2O}} is the water fugacity, \eqn{T} is the temperature,
#' \eqn{H} is the enthalpy, and \eqn{R} is the ideal gas constant. The flow
#' parameters are the prefactor \eqn{A}, and the exponents \eqn{n}, \eqn{m}, and \eqn{r}.
#'
#' To propagate the uncertainties of the flow parameters Monte Carlo simulation is used here.
#' \itemize{
#' \item{ If the flow law parameters are given by a mean value and a marginal error (\eqn{\mu \pm z}),
#' the Monte Carlo simulation assumes a normal distribution given by
#' \eqn{X = N\left(\mu, \sigma\right)}, where \eqn{\mu} is the mean and \eqn{\sigma}
#'  is the standard deviation of the mean (\eqn{\sigma = \text{z}/1.96}).}
#' \item{ If the parameter is given by a range of possible values \eqn{\left[x_\text{min}, x_\text{max}\right]},
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
#' Hirth, G., Teyssier, C., \& Dunlap, W. J. (2001). An evaluation of quartzite flow laws based on comparisons between experimentally and naturally deformed rocks. International Journal of Earth Sciences, 90(1), 77-87. \doi{10.1007/s005310000152}
#'
#' Tokle, L., Hirth, G., \& Behr, W. M. (2019). Flow laws and fabric transitions in wet quartzite. Earth and Planetary Science Letters, 505, 152-161. \doi{10.1016/j.epsl.2018.10.017}
#'
#' Lu, L. X., \& Jiang, D. (2019). Quartz Flow Law Revisited: The Significance of Pressure Dependence of the Activation Enthalpy. Journal of Geophysical Research: Solid Earth, 124(1), 241–256. \doi{10.1029/2018JB016226}
#'
#' Lusk, A. D. J., Platt, J. P., \& Platt, J. A. (2021). Natural and Experimental Constraints on a Flow Law for Dislocation‐Dominated Creep in Wet Quartz. Journal of Geophysical Research: Solid Earth, 126(5), 1-25. \doi{10.1029/2020JB021302}
#'
#' Richter, B., Stünitz, H., \& Heilbronner, R. (2018). The brittle-to-viscous transition in polycrystalline quartz: An experimental study. Journal of Structural Geology, 114(September 2017), 1-21. \doi{10.1016/j.jsg.2018.06.005}
#'
#' @export
#'
#' @seealso [units::set_units()]
#'
#' @importFrom stats rnorm
#' @importFrom units set_units
#'
#' @examples
#' set.seed(20250411)
#' stress <- units::set_units(100, MPa)
#' temperature <- units::set_units(300, degC)
#' pressure <- units::set_units(400, MPa)
#' fugacity <- ps_fugacity(pressure, temperature)
#'
#' creep_quartz(stress = stress, temperature = temperature, model = "Paterson1990")
#' creep_quartz(stress = stress, temperature = temperature, fugacity = fugacity, model = "Hirth2001")
#' creep_quartz(stress = stress, temperature = temperature, fugacity = fugacity, model = "Rutter2004")
#' creep_quartz(stress = stress, temperature = temperature, fugacity = fugacity, model = "Lu2019")
creep_quartz <- function(stress, temperature, fugacity = NULL,
                         grainsize = NULL, pressure = NULL,
                         sim = 1e6,
                         model = c("Hirth2001", "Paterson1990", "Kronenberg1984", "Luan1992", "Gleason1995", "Gleason1995_melt",  "Rutter2004", "Fukuda2018_LT","Fukuda2018_HT", "Richter2018", "Lu2019", "Tokle2019_LT", "Tokle2019_HT", "Lusk2021", "Lusk2021_HP")) {
  model <- match.arg(model)

  # stress, pressure and fugacity in Mega-Pascal, temperature in Kelvins
  temperature <- units::set_units(temperature, "K") |> as.numeric()
  stress <- units::set_units(stress, "MPa") |> as.numeric()

  args <- list(stress = stress, temperature = temperature, sim = sim)

  if (!is.null(fugacity)) args$fugacity <- units::set_units(fugacity, "MPa") |> as.numeric()

  if (!is.null(pressure)) args$pressure <- units::set_units(pressure, "MPa") |> as.numeric()

  if (!is.null(grainsize)) args$grainsize <- units::set_units(grainsize, "um") |> as.numeric()

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
    units::set_units("1/s")

  if (length(edot) > 1) {
    class(edot) <- append("MCS_log", class(edot))
  }
  return(edot)
}
