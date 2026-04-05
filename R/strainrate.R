flow_models <- function() {
  ci2sd <- 1 / 1.96  # convert 95% CI to SD

  list(
    Kronenberg1984 = list(
      H = min_max_std(120, 150), H_std = min_max_std(120, 150)[2], n = min_max_std(2.9, 3.2)[1], n_std = min_max_std(2.9, 3.2)[2], m = 0.18
    ),
    Paterson1990 = list(
      A = 6.5e-8, H = 135, n = 3
    ),
    Luan1992 = list(
      A = 4e-10, H = 152, H_std = 71 * ci2sd, n = 4, n_std = 0.8 * ci2sd
    ),
    Gleason2004 = list(
      A = 1.1e-4, H = 223, H_std = 56 * ci2sd, n = 4, n_std = 0.9 * ci2sd
    ),
    Rutter2004 = list(
      A = 1.2e-5, H = 242, H_std = 24 * ci2sd, n = 2.97, n_std = 0.29 * ci2sd, r = 1
    ),
    Fukuda2018 = list(
      H = 129, H_std = 33 * ci2sd, n = min_max_std(2.9, 5.2)[1], n_std = min_max_std(2.9, 5.2)[2], r = 1
    ),
    Fukuda2018_2 = list(
      A = exp(-2.97), H = 183, H_std = 25 * ci2sd, n = 1.7, n_std = 0.2 * ci2sd, m = 0.51, m_std = 0.51 * ci2sd, r = 1, r_std = 0.2 * ci2sd
    ),
    Richter2018 = list(
      A = 3.1e-4, H = min_max_std(168, 170)[1], H_std = min_max_std(168, 170)[2], n = 1.9, n_std = 0.6 * ci2sd, m = 1.08
    ),
    Hirth2001 = list(
      log_A = -11.2, log_A_std = 0.6 * ci2sd, # prefactor in MPa^{-n} / s
      H = 135, # Activation enthalpy in kJ/mol
      H_std = 15 * ci2sd,
      n = 4, # stress exponent
      r = 1 # water fugacity exponent
    ),
    Lu2019 = list(
      A = 6e-15, H = 132, H_std = 5 * ci2sd, V = 35.3, n = 4, r = 2.7
    ),
    Tokle2019 = list(
      A = 8e-12, H = 140, H_std = 15 * ci2sd, n = 4, n_std = 0.3 * ci2sd, r = 1
    ),
    Tokle2019_2 = list(
      A = 5.e-12, H = 105, H_std = 15 * ci2sd, n = 2.7, n_std = 0.3 * ci2sd, r = 1.1
    ),
    Lusk2021 = list( # "low-pressure" <560 MPa
      log_A = -9.3, # MPa^(-n-r) s^-1
      log_A_std = 0.66 * ci2sd,
      n = 3.5, n_std = 0.2 * ci2sd,
      r = 0.49, r_std = 0.13,
      Q = 118, Q_std = 5 * ci2sd, # kJ mol-1
      V = 2.59, V_std = 2.45 * ci2sd # cm3 mol-1
    ),
    Lusk2021_HP = list( # “high-pressure” 700–1600 MPa
      log_A = -7.90, log_A_std = 0.34 * ci2sd, # MPa−n−r s−1;
      n = 2.0, n_std = 0.1 * ci2sd,
      r = 0.49, r_std = 0.13 * ci2sd,
      Q = 77, Q_std = 8 * ci2sd, # kJ mol−1;
      V = 2.59, V_std = 2.45 * ci2sd # cm3 mol−1
    )
  )
}
#' Strain rate
#'
#' Calculates strain rate from stress, temperature, and grain size.
#' Uses Monte Carlo sampling if flow model parameters contain uncertainties.
#'
#' @param stress Deviatoric stress in MPa or units object
#' @param temperature Temperature in Kelvin or units object
#' @param fugacity Water fugacity in MPa or units object
#' @param grainsize Grainsize in cm  or units object
#' @param pressure Pressure in MPa or units object
#' @param sim non-negative number. Number of Monte Carlo simulations
#' @param model character. Flow law.
#'
#' @details General flow law: \deqn{\dot{\epsilon} = A \sigma^n d^m f_{H_2O}^r \, e^{\left({\frac{-H}{RT}}\right)}}
#'
#' @returns list. Strain rate in 1/s given as median of the Monte Carlo simulations `median`, mean (`mean`), standard error of `log(samples)` (`stderr_log`), the 95% confidence interval of the mean (`conf.int`), and the Monte Carlo simulation (`samples`)
#'
#' @references
#'  Hirth, G., Teyssier, C., & Dunlap, W. J. (2001). An evaluation of quartzite flow laws based on comparisons between experimentally and naturally deformed rocks. International Journal of Earth Sciences, 90(1), 77–87. https://doi.org/10.1007/s005310000152
#'
#'  Tokle, L., Hirth, G., & Behr, W. M. (2019). Flow laws and fabric transitions in wet quartzite. Earth and Planetary Science Letters, 505, 152–161. https://doi.org/10.1016/j.epsl.2018.10.017
#'
#'  Lusk, A. D. J., Platt, J. P., & Platt, J. A. (2021). Natural and Experimental Constraints on a Flow Law for Dislocation‐Dominated Creep in Wet Quartz. Journal of Geophysical Research: Solid Earth, 126(5), 1–25. https://doi.org/10.1029/2020JB021302
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' stress <- units::set_units(100, MPa)
#' temperature <- units::set_units(300, degC)
#' pressure <- units::set_units(400, MPa)
#' fugacity <- ps_fugacity(pressure, temperature)
#' strain_rate(stress = stress, temperature = temperature, model = "Paterson1990")
#' strain_rate(stress = stress, temperature = temperature, fugacity = fugacity, model = "Hirth2001")
#' strain_rate(stress = stress, temperature = temperature, fugacity = fugacity, model = "Rutter2004")
#' strain_rate(stress = stress, temperature = temperature, fugacity = fugacity, model = "Lu2019")
strain_rate <- function(stress, temperature, fugacity = NULL,
                        grainsize = NULL, pressure = NULL,
                        sim = 1e6,
                        model = c("Hirth2001", "Paterson1990", "Kronenberg1984", "Luan1992", "Gleason2004", "Rutter2004", "Fukuda2018", "Richter2018", "Lu2019", "Tokle2019", "Lusk2021")) {
  model <- match.arg(model)

  # tress, pressure and fugacity in Mega-Pascal, temperature in Kelvins
  temperature <- units::set_units(temperature, K) |> as.numeric()
  stress <- units::set_units(stress, MPa) |> as.numeric()
  if (!is.null(fugacity)) {
    fugacity <- units::set_units(fugacity, MPa) |> as.numeric()
  } else {
    fugacity <- 1
  }
  if (!is.null(pressure)) pressure <- units::set_units(pressure, MPa) |> as.numeric()

  if (!is.null(grainsize)) {
    grainsize <- units::set_units(grainsize, cm) |> as.numeric()
  } else {
    grainsize <- 1
  }

  # ideal gas constant
  R <- set_units(8.31446261815324, J * K^-1 * mol^-1) |>
    # set_units(MPa*cm3*K^-1*mol^-1) |>
    set_units(kJ * K^-1 * mol^-1) |>
    as.numeric()

  # Validate model name
  fm <- flow_models()
  if (!model %in% names(fm)) {
    stop(paste("Unknown model '", model, "'. Available models: ",
      paste(names(fm), collapse = ", "),
      sep = ""
    ))
  }

  # Load model parameters
  p <- fm[[model]]

  # prefactor (log MPa^-n s-1)
  A <- if ("A" %in% names(p)) {
    if ("A_std" %in% names(p)) {
      rnorm(sim, mean = p$A, sd = p$A_std)
    } else {
      p$A
    }
  } else if ("log_A" %in% names(p)) {
    log_A <- p$log_A

    if ("log_A_std" %in% names(p)) {
      log_A_std <- p$log_A_std
      log_A <- rnorm(sim, mean = log_A, sd = log_A_std)
    }
    10^log_A
  } else {
    1
  }

  # Enthalpy (kJ mol-1)
  H <- if ("H" %in% names(p)) {
    if ("H_std" %in% names(p)) {
      rnorm(sim, mean = p$H, sd = p$H_std)
    } else {
      p$H
    }
  } else {
    1
  }

  # Stress exponent
  n <- if ("n" %in% names(p)) {
    if ("n_std" %in% names(p)) {
      rnorm(sim, mean = p$n, sd = p$n_std)
    } else {
      p$n
    }
  } else {
    0
  }

  # Grain size expopnent
  r <- if ("r" %in% names(p)) {
    if ("r_std" %in% names(p)) {
      rnorm(sim, mean = p$r, sd = p$r_std)
    } else {
      p$r
    }
  } else {
    0
  }

  # Fugacity exponent
  m <- if ("m" %in% names(p)) {
    if ("m_std" %in% names(p)) {
      rnorm(sim, mean = p$m, sd = p$m_std)
    } else {
      p$m
    }
  } else {
    0
  }

  if (all(c("Q", "V") %in% names(p))) {
    Q <- if ("Q_std" %in% names(p)) {
      rnorm(sim, mean = p$Q, sd = p$Q_std)
    } else {
      p$Q
    }

    V <- if ("V_std" %in% names(p)) {
      rnorm(sim, mean = p$V, sd = p$V_std)
    } else {
      p$V
    }

    H <- Q + V * pressure
  }

  edot <- A * stress^n * grainsize^m * fugacity^r * exp(-H / (R * temperature))

  if (length(edot) == 1) {
    set_units(edot, "1/s")
  } else {
    mc_stats(edot, "1/s")
  }
}
