# Functions to calculate H2O molar volume (ps_volume) and fugacity (ps_fugacity)
# using the Pitzer and Sterner equation of state.
# Pitzer, K.S. and Sterner, S.M., 1994. Equations of state valid
# continuously from zero to extreme pressures for H2O and CO2.
# Journal of Chemical Physics. 101: 3111-3116.

# Coefficient matrix (10 rows x 6 columns)
ps_coeff <- function() {
  matrix(c(
    0, 0, 0.24657688e6, 0.51359951e2, 0, 0,
    0, 0, 0.58638965e0, -0.28646939e-2, 0.31375577e-4, 0,
    0, 0, -0.62783840e1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7,
    0, 0, 0, -0.42719875e0, -0.16325155e-4, 0,
    0, 0, 0.56654978e4, -0.16580167e2, 0.76560762e-1, 0,
    0, 0, 0, 0.10917883e0, 0, 0,
    0.38878656e13, -0.13494878e9, 0.30916564e6, 0.75591105e1, 0, 0,
    0, 0, -0.65537898e5, 0.18810675e3, 0, 0,
    -0.14182435e14, 0.18165390e9, -0.19769068e6, -0.23530318e2, 0, 0,
    0, 0, 0.92093375e5, 0.12246777e3, 0, 0
  ), nrow = 10, ncol = 6, byrow = TRUE)
}


ps_eos <- function(volume, temperature, targetP) {
  # volume in cm3/mol, temperature in Kelvins, targetP in bars
  R_const <- 8314462.61815324 # ideal gas constant: Pa*cm3/K/mol
  den <- 1 / volume
  cv <- numeric(10)
  coeff <- ps_coeff()
  for (i in 1:10) {
    cv[i] <- coeff[i, 1] * temperature^-4 +
      coeff[i, 2] * temperature^-2 +
      coeff[i, 3] * temperature^-1 +
      coeff[i, 4] +
      coeff[i, 5] * temperature +
      coeff[i, 6] * temperature^2
  }
  pressure <- (den + cv[1] * den^2 - den^2 * ((cv[3] + 2 * cv[4] * den + 3 * cv[5] * den^2
                                               + 4 * cv[6] * den^3) / (cv[2] + cv[3] * den + cv[4] * den^2 + cv[5] * den^3
                                                                       + cv[6] * den^4)^2) + cv[7] * den^2 * exp(-cv[8] * den)
               + cv[9] * den^2 * exp(-cv[10] * den)) * R_const * temperature / 1e5
  return(pressure - targetP) # bars
}


#' Water fugacity and H\eqn{_2}O molar volume
#'
#' calculate H\eqn{_2}O molar volume and fugacity using the Pitzer and Sterner (1994)
#' equation of state.
#'
#' @param pressure numeric. Pressure either in bar or as `units` object
#' @param temperature numeric. Temperature either in Kelvin or as `units` object
#'
#' @references Pitzer, K.S. and Sterner, S.M., 1994. Equations of state valid
#' continuously from zero to extreme pressures for H2O and CO2.
#' Journal of Chemical Physics. 101: 3111-3116.
#'
#' @source This is modified version of the `fugacity.py` script: https://github.com/forsterite/fugacity/tree/master
#' This R version uses a precise gas constant and accepts inputs in different units.
#'
#' @returns units object
#'
#' @importFrom nleqslv nleqslv
#' @importFrom units set_units
#' @name pitzer
#'
#' @seealso [units::set_units()]
#'
#' @examplesIf require(units)
#' pressure <- units::set_units(1, atm)
#' temperature <- units::set_units(25, degC)
#' ps_volume(pressure, temperature) # 18.7231
#' ps_fugacity(pressure, temperature) # 0.04946
#'
#' temperature2 <- units::set_units(300, degC)
#' pressure2 <- units::set_units(400, MPa)
#' ps_fugacity(pressure2, temperature2) # 37 MPa
#' ps_volume(pressure2, temperature2)
NULL

#' @rdname pitzer
#' @export
ps_volume <- function(pressure, temperature) {
  # pressure in bars, temperature in Kelvins
  pressure <- units::set_units(pressure, bar) |> as.numeric()
  temperature <- units::set_units(temperature, K) |> as.numeric()

  result <- nleqslv::nleqslv(10, ps_eos, temperature = temperature, targetP = pressure)
  units::set_units(result$x, cm3 / mol)
}

#' @rdname pitzer
#' @export
ps_fugacity <- function(pressure, temperature) {
  # pressure in bars, temperature in Kelvins
  pressure <- units::set_units(pressure, bar) |> as.numeric()
  temperature <- units::set_units(temperature, K) |> as.numeric()

  R_const <- 8314462.61815324 # ideal gas constant: Pa*cm3/K/mol

  coeff <- ps_coeff()
  cv <- numeric(10)
  for (i in 1:10) {
    cv[i] <- coeff[i, 1] * temperature^-4 +
      coeff[i, 2] * temperature^-2 +
      coeff[i, 3] * temperature^-1 +
      coeff[i, 4] +
      coeff[i, 5] * temperature +
      coeff[i, 6] * temperature^2
  }
  volume <- ps_volume(pressure, temperature) |> as.numeric()
  den <- 1 / volume
  fug <- exp(log(den) + cv[1] * den + (1 / (cv[2] + cv[3] * den + cv[4] * den^2
                                            + cv[5] * den^3 + cv[6] * den^4) - 1 / cv[2])
             - cv[7] / cv[8] * (exp(-cv[8] * den) - 1)
             - cv[9] / cv[10] * (exp(-cv[10] * den) - 1)
             + pressure * 1e5 / (den * R_const * temperature)
             + log(R_const * temperature) - 1) / 1e5
  units::set_units(fug, bar) # bars
}
