# Peierls Creep

Strain rate for low-temperature, plastic deformation mechanism that
corresponds to dislocation glide at high stress (differential stress
larger than 200 MPa).

## Usage

``` r
peierls_creep(
  stress,
  temperature,
  model = c("Goetze1979", "Demouchy2013", "Mei2010"),
  sim = 1e+06
)
```

## Arguments

- stress:

  Differential stress in MPa or `units` object

- temperature:

  Temperature in Kelvin or `units` object

- model:

  character. Model to be used, one of `"Goetze1979"` (Goetze and Evans,
  1979), `"Demouchy2013"` (Demouchy et al, 2013), and `"Mei2010"` (Mei
  et al., 2010)

- sim:

  integer. Number of Monte Carlo simulations.

## Value

Strain rate in 1/s

## Examples

``` r
peierls_creep(300, 500, model = "Goetze1979")
#> 4.407903e-41 [1/s]
peierls_creep(300, 500, model = "Demouchy2013")
#> 2.21843e-28 [1/s]
peierls_creep(300, 500, model = "Mei2010") |> summary()
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 0.000e+00 0.000e+00 1.037e-23 0.000e+00 7.382e-19 
```
