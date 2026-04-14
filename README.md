
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geocreep

<!-- badges: start -->

[![R-CMD-check](https://github.com/tobiste/geocreep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tobiste/geocreep/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of geocreep is to quick and easy calculate strain rates,
fugacities and differential stress from grain sizes using R. The codes
use Monte Carlo simulations to propagate uncertainties of the flow
parameters into the final estimate.

## Installation

You can install the development version of geocreep from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("tobiste/geocreep")
```

## Example

``` r
library(geocreep)
library(units) # load this package to define your units

# Set seed for reproducibility
set.seed(20250411)

# define a number of Monte Carlo simulations
nmc <- 1e3 # here a small number has been chosen just for demonstration. In practice you would use a number >>1e5
```

### Water fugacity

Calculate water fugacity from temperature and pressures using the Pitzer
and Sterner (1994) equation:

``` r
#Define temperature and pressure
temperature <- units::set_units(rnorm(nmc, 300, 50/1.96), degC)
pressure <- units::set_units(rnorm(nmc, 400, 10/1.96), MPa)

# Calculate fugacity
fugacity <- ps_fugacity(pressure, temperature)

summary(fugacity)
#> Statistical summary of 1000 Monte Carlo simulations
#> 
#> Median:                      370 bar 
#> 95% interpercentile range:   200 - 590 bar 
#> Standard error in log-space: 0.00375502
#> Student's t-Test:            p<0.05
```

By default, the functions create 100,000 samples for each parameter with
uncertainties, creating 1e6 results. Here, a small number has been
chosen just for demonstration purpose. In practice, however, you would
use a number \>\>1e5.

In general, the Monte Carlo simulations results do not follow a normal
(symmetric) distribution because some parameters are exponents in the
power laws. Therefore, the median and interpercentile range provide the
best estimators for average and dispersion, rather than mean and
standard deviation.

### Grain-size piezometry

Calculating differential stress from grain size (e.g. using the Stipp
and Tullis, 2003, piezometer):

``` r
# Define grain size
grainsize <- units::set_units(11, um)

# Calculate equivalent differential stress
stress <- grainsize_piezometry(grainsize, model = "Stipp-reg2-3", sim = nmc)

# Summary stats of Monte Carlo simulation
summary(stress)
#> Statistical summary of 1000 Monte Carlo simulations
#> 
#> Median:                      100 MPa 
#> 95% interpercentile range:   52 - 200 MPa 
#> Standard error in log-space: 0.00481395
#> Student's t-Test:            p<0.05
```

Note: There is also a subgrain-size piezometer (Goddard et al. 2021):
`subgrainsize_piezometer()`

### Flow laws

Calculate strain rates using a defined flow law for quartz from
differential stress, temperature, and fugacity:

``` r
# Calculate strain rates using temperature, the fugacity, and the MC estimates for differential stress calculated before
edot <- creep_quartz(stress = stress, temperature = temperature, fugacity = fugacity, model = "Hirth2001", sim = nmc)

# Summary stats of Monte Carlo simulation
summary(edot)
#> Statistical summary of 1000 Monte Carlo simulations
#> 
#> Median:                      1.1e-14 /s
#> 95% interpercentile range:   3.8e-17 - 2.5e-12 /s
#> Standard error in log-space: 0.0380114
#> Student's t-Test:            p<0.05
```

## Author

Tobias Stephan (<tstephan@lakeheadu.ca>)

## Feedback, issues, and contributions

I welcome feedback, suggestions, issues, and contributions! If you have
found a bug, please file
[here](https://github.com/tobiste/geocreep/issues) with minimal code to
reproduce the issue.

## License

GPL-3.0 License
