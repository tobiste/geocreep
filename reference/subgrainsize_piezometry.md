# Subgrain‐Size Piezometer Calibrated for EBSD

Calculates differential stress using the Subgrain‐Size Piezometer of
Goddard et al. (2020). Uses Monte Carlo sampling for propagating
parameter uncertainties in flow model.

## Usage

``` r
subgrainsize_piezometry(
  lambda,
  sd = NULL,
  calibrated = TRUE,
  min = c("q", "fo90", "fo50"),
  sim = 1e+06,
  propagate_err = TRUE
)
```

## Arguments

- lambda:

  mean line intercept length in μm or `units` object.

- sd:

  (optional) Standard deviation of `lambda`

- calibrated:

  logical. Whether the calibration of Holyoke and Kronenberg (2010) is
  considered or not.

- min:

  character. The Mineral uses. one of `"q"` for quartz, `"fo90"` for
  Olive with 90% Forsterite, or `"fo50"` for Olivine with 50%
  Forsterite.

- sim:

  non-negative number. Number of Monte Carlo simulations

- propagate_err:

  logical. Whether errors of the flow law parameters should be
  propagated. `TRUE` by default.

## Value

list. Differential stress in MPa. If Monte Carlo Simulation was used,
and object of class `"MCS"` is returned (see
[`summary()`](https://rdrr.io/r/base/summary.html) for detailed
description of output). The piezometer produce log-normal distributed
estimates considering the uncertainties in the equation parameter. Hence
it is recommended to report the median (or geometric mean), and the
interpercentile range.

## Details

The sub-grain size piezometer is \$\$\frac{\lambda}{b} = 10^a
\left(\frac{\sigma}{\mu}\right)^b\$\$ where \\\lambda\\ is the mean line
intercept length, \\b\\ is the Burgers vector, \\\sigma\\ is the
differential stress, \\\mu\\ is the shear modulus, and \\a\\ and \\b\\
are the empirical exponents.

## References

Goddard, R. M., Hansen, L. N., Wallis, D., Stipp, M., Holyoke, C. W.,
Kumamoto, K. M., & Kohlstedt, D. L. (2020). A Subgrain‐Size Piezometer
Calibrated for EBSD. Geophysical Research Letters, 47(23).
https://doi.org/10.1029/2020GL090056

## See also

[`grainsize_piezometry()`](https://tobiste.github.io/geocreep/reference/grainsize_piezometry.md),
[`units::set_units()`](https://r-quantities.github.io/units/reference/units.html)

## Examples

``` r
set.seed(20250411)
subgrainsize_piezometry(9, min = "fo50") |> summary() # 420 MPa
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      560 MPa 
#> 95% interpercentile range:   34 - 3600 MPa 
#> Standard error in log-space: 0.000516207
#> Student's t-Test:            p<0.05
subgrainsize_piezometry(18, min = "q") |> summary() # 240 MPa
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      220 MPa 
#> 95% interpercentile range:   11 - 1500 MPa 
#> Standard error in log-space: 0.00054436
#> Student's t-Test:            p<0.05
```
