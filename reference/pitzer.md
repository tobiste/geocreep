# Fugacity and molar volume for H\\\_2\\O or CO\\\_2\\

calculate molar volume and fugacity for H\\\_2\\O or CO\\\_2\\ using the
Pitzer and Sterner (1994) equation of state.

## Usage

``` r
ps_volume(pressure, temperature, phase = c("H2O", "CO2"))

ps_fugacity(pressure, temperature, phase = c("H2O", "CO2"), ...)
```

## Source

This is modified version of the `fugacity.py` script:
https://github.com/forsterite/fugacity/tree/master This R version uses a
precise gas constant, accepts inputs in different units, and adds the
CO\\\_2\\ parameters.

## Arguments

- pressure:

  numeric. Pressure either in bar or as `units` object

- temperature:

  numeric. Temperature either in Kelvin or as `units` object

- phase:

  character. Fluid phase for which fugacity or volume should be
  calculated for; one of `"H2O"` (the default) and `"CO2"`.

- ...:

  optional arguments passed to
  [`future.apply::future_mapply()`](https://future.apply.futureverse.org/reference/future_mapply.html)

## Value

units object

## References

Pitzer, K.S. & Sterner, S.M., 1994. Equations of state valid
continuously from zero to extreme pressures for H2O and CO2. Journal of
Chemical Physics. 101: 3111-3116.

## See also

[`units::set_units()`](https://r-quantities.github.io/units/reference/units.html)

## Examples

``` r
pressure <- units::set_units(1, atm)
temperature <- units::set_units(25, degC)
ps_volume(pressure, temperature) # 18.7231
#> 18.72311 [cm^3/mol]
ps_fugacity(pressure, temperature) # 0.04946
#> 0.04946146 [bar]

temperature2 <- units::set_units(300, degC)
pressure2 <- units::set_units(400, MPa)
ps_fugacity(pressure2, temperature2) # 37 MPa
#> 371.3371 [bar]
ps_volume(pressure2, temperature2)
#> 18.69728 [cm^3/mol]
```
