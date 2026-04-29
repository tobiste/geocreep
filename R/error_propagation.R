#' Algebraic error propagation of flow laws
#'
#' Uses 1st order Taylor expansion (linearization) of the natural log-space flow law
#'
#' @inheritParams creep_quartz
#'
#' @returns list.\describe{
#' \item{`"e_best"`}{Strain rate (1/s)}
#' \item{`"sd_e_range"`}{Range of 1s standard deviation  (1/s)}
#' \item{`"var_log_e_total"`}{Variance in natural log-space}
#' \item{`"sd_log_e_total`}{Standard deviation in natural log-space}
#' \item{`"var_log_e"`}{Individual variance components in natural log-space}
#' }
#' @export
#'
#' @seealso [creep_quartz()]
#'
#' @examples
#' stress <- units::set_units(100, MPa)
#' temperature <- units::set_units(300, degC)
#' pressure <- units::set_units(400, MPa)
#' fugacity <- ps_fugacity(pressure, temperature)
#'
#' creep_quartz_analytic(
#'   stress = stress, temperature = temperature, fugacity = fugacity, pressure = pressure,
#'   model = "Hirth2001"
#' )
#'
#' creep_quartz_analytic(
#'   stress = stress, temperature = temperature, fugacity = fugacity, pressure = pressure,
#'   model = "Lusk2021_LP"
#' )
#'
#' creep_quartz_analytic(
#'   stress = stress, temperature = temperature, fugacity = fugacity, pressure = pressure,
#'   model = "Kronenberg1984"
#' )
creep_quartz_analytic <- function(stress, temperature, fugacity = NULL,
                                  grainsize = NULL, pressure = NULL,
                                  model = c("Hirth2001", "Paterson1990", "Kronenberg1984", "Luan1992", "Gleason1995", "Gleason1995_melt", "Rutter2004", "Fukuda2018_LT", "Fukuda2018_HT", "Richter2018", "Lu2019", "Tokle2019_LT", "Tokle2019_HT", "Lusk2021_LP", "Lusk2021_HP")) {
  model <- match.arg(model)

  # stress, pressure and fugacity in Mega-Pascal, temperature in Kelvins
  temperature <- units::set_units(temperature, "K") #|> as.numeric()
  stress <- units::set_units(stress, "MPa") #|> as.numeric()

  if (!is.null(fugacity)) {
    fugacity <- units::set_units(fugacity, "MPa") #|> as.numeric()
  } else {
    fugacity <- units::set_units(1, "MPa")
  }

  if (!is.null(pressure)) {
    pressure <- units::set_units(pressure, "MPa") #|> as.numeric()
  } else {
    pressure <- units::set_units(1, "MPa")
  }

  if (!is.null(grainsize)) {
    grainsize <- units::set_units(grainsize, "um") #|> as.numeric()
  } else {
    grainsize <- units::set_units(1, "um") # |> as.numeric()
  }

  # assuming normal distribution to calculate sd; if there is only one value, sd() returns NA, and the propagated error will also be NA as it has no effect
  temperature_sd <- sd(temperature) |> replace_na_with_zero()
  units(temperature_sd) <- units(temperature)
  stress_sd <- sd(stress) |> replace_na_with_zero()
  units(stress_sd) <- units(stress)
  grainsize_sd <- sd(grainsize) |> replace_na_with_zero()
  units(grainsize_sd) <- units(grainsize)
  pressure_sd <- sd(pressure) |> replace_na_with_zero()
  units(pressure_sd) <- units(pressure)
  fugacity_sd <- sd(fugacity) |> replace_na_with_zero()
  units(fugacity_sd) <- units(fugacity)

  temperature_mean <- mean(temperature)
  stress_mean <- mean(stress)
  grainsize_mean <- mean(grainsize)
  pressure_mean <- mean(pressure)
  fugacity_mean <- mean(fugacity)

  # load model parameters:
  x <- flow_model_params2(model)
  R <- gas_const() #|> as.numeric()
  RT <- R * temperature_mean

  # 1st order Taylor expansion:
  prefactor_error <- (log(10) * x$log_A_sd)^2
  stress_exponent_error <- (log(stress_mean) * x$n_sd)^2
  fugacity_exponent_error <- (log(fugacity_mean) * x$r_sd)^2
  grainsize_exponent_error <- (log(grainsize_mean) * x$m_sd)^2
  stress_error <- (stress_mean * x$n  / stress_sd)^2
  fugacity_error <- (fugacity_mean * x$r / fugacity_sd)^2
  grainsize_error <- (grainsize_mean * x$m / grainsize_sd)^2

  # replace NA, NaN and Inf with zero
  stress_error <- replace_with_zero(stress_error)
  fugacity_error <- replace_with_zero(fugacity_error)
  grainsize_error <- replace_with_zero(grainsize_error)

  if (is.null(x$Q)) {
    H <- units::set_units(x$H, "kJ mol-1") #|> as.numeric()
    H_sd <- units::set_units(x$H_sd, "kJ mol-1") #|> as.numeric()
    enthalpy_error <- H_sd / (R * temperature_mean^2)
  } else {
    Q <- units::set_units(x$Q, "kJ mol-1") # |> as.numeric()
    Q_sd <- units::set_units(x$Q_sd, "kJ mol-1") #|> as.numeric()
    V <- units::set_units(x$V, "cm3 mol-1") # |> as.numeric()
    V_sd <- units::set_units(x$V_sd, "cm3 mol-1") # |> as.numeric()

    H <- Q + V * pressure_mean
    #enthalpy_error <- 1 / (RT)^2 * (Q_sd^2 + (pressure_mean * V_sd)^2 + (V * pressure_sd)^2)
    enthalpy_error <-  (Q_sd/RT)^2 + (pressure_mean * V_sd/RT)^2 + (V * pressure_sd/RT)^2
  }

  temperature_error <- (H * temperature_sd / (RT*temperature_mean))^2

  # merge all contributors
  errors <- c("prefactor" = prefactor_error, "stress_exponent" = stress_exponent_error, "fugacity_exponent" = fugacity_exponent_error, "grainsize_exponent" = grainsize_exponent_error, "stress" = stress_error, "fugacity" = fugacity_error, "grainsize" = grainsize_error, "temperature" = temperature_error, "enthalpy" = enthalpy_error)


  # sum of all contributors is the log-variance
  var_log_e <- sum(errors, na.rm = TRUE)
  sd_log_e <- sqrt(var_log_e)

  # calculate 'best-fit' strain rate:
  term <- -H / RT
  stopifnot(units(term) == units::unitless)
  arrhenius <- exp(as.numeric(term))
  epsilon <- 10^replace_with_zero(x$log_A) * as.numeric(stress_mean)^x$n * as.numeric(fugacity_mean)^x$r * as.numeric(grainsize_mean)^x$m * arrhenius
  e_best <- units::set_units(epsilon, "s-1")
  sd_e_range <- e_best * exp(c(-sd_log_e, sd_log_e))

  # return
  list(e_best = e_best, sd_e_range = sd_e_range, var_log_e_total = var_log_e, sd_log_e_total = sd_log_e, var_log_e = errors)
}


flow_model_params2 <- function(model) {
  x <- flow_model_params(model)

  if (!is.null(x$n_min)) {
    x$n <- mean(c(x$n_min, x$n_max))
    x$n_sd <- (x$n_max - x$n_min)^2 / 12
  }
  if (!is.null(x$H_min)) {
    x$H <- mean(c(x$H_min, x$H_max))
    x$H_sd <- (x$H_max - x$H_min)^2 / 12
  }
  if (!is.null(x$A_factor)) {
    x$log_A <- x$A_mean + log10(x$A_factor)
    x$log_A_sd <- x$A_sd
  }
  if (!is.null(x$A) & is.null(x$A_k)) {
    x$log_A <- log10(x$A)
    x$log_A_sd <- 0
  }
  if (!is.null(x$A_k)) {
    x$log_A <- x$A_k + log10(x$A)
    x$log_A_sd <- x$A_sd / (x$A * log(10))
  }
  x
}

#' @noRd
flow_model_params <- function(model) {
  list(
    Kronenberg1984 = list(
      log_A = NA,
      log_A_sd = 0,
      H_min = 120,
      H_max = 150,
      n_min = 2.9,
      n_max = 3.2,
      m = 0.18,
      m_sd = 0,
      r = 0, r_sd = 0
    ),
    Paterson1990 = list(
      A = 6.5e-8,
      A_sd = 0,
      H = 135,
      H_sd = 0,
      n = 3.1,
      n_sd = 0,
      m = 0, m_sd = 0,
      r = 0, r_sd = 0
    ),
    Luan1992 = list(
      # uncertainties not specified; assuming 1s
      A = 4e-10,
      A_sd = 0,
      H = 152,
      H_sd = 71,
      n = 4,
      n_sd = 0.8,
      r = 0, r_sd = 0,
      m = 0, m_sd = 0
    ),
    Gleason1995 = list(
      # uncertainties not specified; assuming 1s
      A_factor = 1.1,
      A_mean = -4,
      A_sd = 2,
      H = 223,
      H_sd = 56,
      n = 4,
      n_sd = 0.9,
      r = 0, r_sd = 0,
      m = 0, m_sd = 0
    ),
    Gleason1995_melt = list(
      # uncertainties not specified; assuming 1s
      A_factor = 1.8,
      A_mean = -8,
      A_sd = 2,
      H = 137,
      H_sd = 34,
      n = 4,
      n_sd = 0.9,
      r = 0, r_sd = 0,
      m = 0, m_sd = 0
    ),
    Rutter2004 = list(
      log_A = -4.93,
      log_A_sd = 0.34,
      H = 242,
      H_sd = 24, # 1s
      n = 2.97,
      n_sd = 0.29, # 1s
      r = 1,
      r_sd = 0,
      m = 0, m_sd = 0
    ),
    Fukuda2018_LT = list(
      # uncertainties not specified; assuming 1s
      log_A = -2.97,
      log_A_sd = 0.23,
      H = 129,
      H_sd = 33, # std given as 1s?
      n_min = 2.9,
      n_max = 5.2,
      r = 1,
      r_sd = 0,
      m = 0, m_sd = 0
    ),
    Fukuda2018_HT = list(
      log_A = -2.97,
      log_A_sd = 0.23,
      H = 183.0,
      H_sd = 25, # std given as 1s?
      n = 1.7,
      n_sd = 0.2, # std given as 1s?
      m = -0.51,
      m_sd = 0.13, # std given as 1s?
      r = 1.0,
      r_sd = 0.2 # std given as 1s?
    ),
    Richter2018 = list(
      A = 3.1e-4,
      A_sd = 0,
      # H_min =168
      # H_max =170
      H = 170,
      H_sd = 72,
      n = 1.9,
      n_sd = 0.6,
      m = 1.08,
      m_sd = 0,
      r = 0, r_sd = 0
    ),
    Hirth2001 = list(
      # uncertainties not specified; assuming 1s
      log_A = -11.2,
      log_A_sd = 0.6, # prefactor in MPa^{-n} / s
      H = 135, # Activation enthalpy in kJ/mol
      H_sd = 15,
      n = 4, # stress exponent; adopted from Gleason and Tullis (1995) and Luan and Paterson (1992)
      n_sd = 0.85,
      r = 1, # water fugacity exponent
      r_sd = 0,
      m = 0, m_sd = 0
    ),
    Lu2019 = list(
      # uncertainties not specified; assuming 1s
      A = 6,
      A_sd = 5,
      A_k = -15,
      H = 132,
      H_sd = 5,
      V = 35.3, # |> set_units("cm3 mol-1") # cm3/mol
      n = 4, # adopted from Gleason and Tullis (1995) and Luan and Paterson (1992)
      n_sd = 0.85, # mean std from Gleason and Tullis (1995) and Luan and Paterson (1992)
      r = 2.7,
      r_sd = 0,
      m = 0, m_sd = 0
    ),
    Tokle2019_HT = list(
      # uncertainties not specified; assuming 1s
      A = 8e-12,
      A_sd = 0,
      H = 140,
      H_sd = 15,
      n = 4,
      n_sd = 0.3,
      r = 1,
      r_sd = 0,
      m = 0, m_sd = 0
    ),
    Tokle2019_LT = list(
      # uncertainties not specified; assuming 1s
      A = 5.4e-12,
      A_sd = 0,
      H = 105,
      H_sd = 15,
      n = 2.7,
      n_sd = 0.3,
      r = 1.1,
      r_sd = 0,
      m = 0, m_sd = 0
    ),
    Lusk2021_LP = list(
      # "low-pressure" <560 MPa

      # std given as 1sd
      log_A = -9.3, # MPa^(-n-r) s^-1
      log_A_sd = 0.66,
      n = 3.5,
      n_sd = 0.2,
      r = 0.49,
      r_sd = 0.13,
      Q = 118,
      Q_sd = 5, # kJ mol-1
      V = 2.59,
      V_sd = 2.45, # cm3 mol-1
      m = 0, m_sd = 0
    ),
    Lusk2021_HP = list(
      # “high-pressure” 700–1600 MPa

      # std given as 1sd
      log_A = -7.90,
      log_A_sd = 0.34, # MPa−n−r s−1;
      n = 2.0,
      n_sd = 0.1,
      r = 0.49,
      r_sd = 0.13,
      Q = 77,
      Q_sd = 8, # kJ mol−1;
      V = 2.59,
      V_sd = 2.45, # cm3 mol−1
      m = 0, m_sd = 0
    )
  )[[model]]
}
