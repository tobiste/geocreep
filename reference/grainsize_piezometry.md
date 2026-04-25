# Grain-Size Piezometry

Calculates differential stress from grain size. Uses Monte Carlo
sampling for propagating parameter uncertainties in flow model.

## Usage

``` r
grainsize_piezometry(
  d,
  sd = NULL,
  model = c("Stipp-reg2-3", "Stripp-reg1", "Cross-1", "Cross-sliding"),
  sim = 1e+06,
  propagate_err = TRUE
)
```

## Arguments

- d:

  numeric. Grain size in μm or `units` object

- sd:

  (optional) Standard deviation of `d`

- model:

  character. One of

  `"Stipp-reg2-3"`

  :   Piezometer for dislocation creep regime 2 and 3 (for deviatoric
      stress \<368 MPa) after Stipp and Tullis (2003)

  `"Stripp-reg1"`

  :   Piezometer for dislocation creep regime 1 (deviatoric stress ≥368
      MPa) after Stipp and Tullis (2003)

  `"Cross-1"`

  :   Piezometer after Cross et al. (2017) for 1 μm step size resolution
      in EBSD data

  `"Cross-sliding"`

  :   Sliding resolution piezometer after Cross et al. 2017. According
      to authors, more accurately estimates stress in fine-grained (\<10
      μm) samples

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

General formula for grain-size piezometer is: \$\$\sigma =
\left(\frac{d}{k}\right)^\frac{1}{n}\$\$ where \\\sigma\\ is the
differential stress, \\d\\ is the grain size, and \\k\\ and \\n\\ are
empirical parameters.

## References

Stipp, M., & Tullis, J. (2003). The recrystallized grain size piezometer
for quartz. Geophysical Research Letters, 30(21), 1–5.
https://doi.org/10.1029/2003GL018444

Cross, A. J., Prior, D. J., Stipp, M., & Kidder, S. (2017). The
recrystallized grain size piezometer for quartz: An EBSD-based
calibration. Geophysical Research Letters, 44(13), 6667–6674.
https://doi.org/10.1002/2017GL073836

## See also

[`subgrainsize_piezometry()`](https://tobiste.github.io/geocreep/reference/subgrainsize_piezometry.md),
[`units::set_units()`](https://r-quantities.github.io/units/reference/units.html)

## Examples

``` r
set.seed(20250411)
grainsize_piezometry(12.2) |> summary() # 92
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      92 MPa 
#> 95% interpercentile range:   49 - 190 MPa 
#> Log-variance:                0.0229959
#> Student's t-Test:            p<0.05
grainsize_piezometry(31) |> summary() # 44
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      44 MPa 
#> 95% interpercentile range:   24 - 86 MPa 
#> Log-variance:                0.019715
#> Student's t-Test:            p<0.05
```
