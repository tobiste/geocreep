#' Peierls Creep
#'
#' Strain rate for low-temperature, plastic deformation mechanism that corresponds to
#' dislocation glide at high stress (differential stress larger than 200 MPa).
#'
#' @param stress Differential stress in MPa or `units` object
#' @param temperature Temperature in Kelvin or `units` object
#' @param model character. Model to be used, one of `"Goetze1979"`
#' (Goetze and Evans, 1979), `"Demouchy2013"` (Demouchy et al, 2013), and
#' `"Mei2010"` (Mei et al., 2010)
#' @param sim integer. Number of Monte Carlo simulations.
#'
#' @returns Strain rate in 1/s
#' @export
#'
#' @examples
#' peierls_creep(300, 500, model = "Goetze1979")
#' peierls_creep(300, 500, model = "Demouchy2013")
#' peierls_creep(300, 500, model = "Mei2010")
peierls_creep <- function(stress, temperature, model = c("Goetze1979", "Demouchy2013", "Mei2010"), sim = 1e6) {
  model <- match.arg(model)

  R <- gas_const() |>  # in SI units
    as.numeric()

  # stress, pressure and fugacity in Mega-Pascal, temperature in Kelvins
  TK <- units::set_units(temperature, "K") |>
    as.numeric()
  sigma_d <- units::set_units(stress, "MPa") |>
    units::set_units('Pa') |>
    as.numeric()

  if (model == "Mei2010") {
    A <- units::set_units(1.4e-7, "MPa-2 s-1") |>
      units::set_units('Pa-2 s-1') |>
      as.numeric()
    sigma_p <- units::set_units(5.9, "GPa") |>
      units::set_units('Pa') |>
      as.numeric()
    H <- rnorm(sim, 320, 50 / 1.96) |>
      units::set_units("kJ mol-1") |>
      units::set_units('J mol-1') |>
      as.numeric()

    edot <- A * sigma_d^2 * exp(-H / (R * TK) * (1 - sqrt(sigma_d / sigma_p)))
  } else {
    if (model == "Goetze1979") {
      sigma_p <- units::set_units(8.5, "GPa") |>
        units::set_units('Pa') |>
        as.numeric()
      H <- units::set_units(536, "kJ mol-1") |>
        units::set_units('J mol-1') |>
        as.numeric()
      A <- units::set_units(5.7e11, "s-1") |>
        as.numeric()
      q <- 1
    } else {
      sigma_p <- units::set_units(15, "GPa") |>
        units::set_units('Pa') |>
        as.numeric()
      H <- units::set_units(450, "kJ mol-1") |>
        units::set_units('J mol-1') |>
        as.numeric()
      A <- units::set_units(10e6, "s-1") |>
        as.numeric()
      q <- 1 / 2
    }

    edot <- A * exp(-H / (R * TK) * (1 - (sigma_d / sigma_p)^q)^2)
  }

  edot <- units::set_units(edot, 's-1')
  if(length(edot)>1){
    class(edot) <- append(class(edot), 'MCS')
  }
    return(edot)

}
