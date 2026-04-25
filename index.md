# geocreep

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
```

### Water fugacity

Calculate water fugacity from temperature and pressures using the Pitzer
and Sterner (1994) equation:

``` r
nmc <- 1e3 # create 1,000 samples for temperature and pressure using a normal distribution

# Define temperature and pressure
temperature <- units::set_units(rnorm(nmc, 300, 50 / 1.96), degC)
pressure <- units::set_units(rnorm(nmc, 400, 100 / 1.96), MPa)

# Calculate fugacity
fugacity <- ps_fugacity(pressure, temperature)

fugacity_stats <- summary(fugacity)
#> Statistical summary of 1000 Monte Carlo simulations
#> 
#> Mean:                    380 bar 
#> 95% confidence interval: 370 - 390 bar 
#> Variance:                17054
#> Student's t-Test:        p<0.05
```

> Here, the Monte Carlo simulation produces 1,000 normal-distributed
> estimates for the fugacity.

### Grain-size piezometry

Calculating differential stress from grain size (e.g. using the Stipp
and Tullis, 2003, piezometer):

``` r
# Define grain size
grainsize <- units::set_units(11, um)

# Calculate equivalent differential stress
stress <- grainsize_piezometry(grainsize, model = "Stipp-reg2-3")

# Summary stats of Monte Carlo simulation
summary(stress)
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      100 MPa 
#> 95% interpercentile range:   53 - 210 MPa 
#> Log-variance:                0.0233963
#> Student's t-Test:            p<0.05
```

> By default, the functions create 1,000,000 samples for each parameter
> with uncertainties, creating 1e6 results. Results of Monte Carlo
> simulations based on flow law equations do not follow a normal
> (symmetric) distribution because some parameters are exponents of the
> power laws. Therefore, the median and interpercentile range provide
> the best estimators for average and dispersion, rather than mean and
> standard deviation.

Note: There is also a subgrain-size piezometer (Goddard et al. 2021):
`subgrainsize_piezometer()`

### Flow laws

Calculate strain rates using a defined flow law for dislocation creep in
quartz from differential stress, temperature, and fugacity:

``` r
fug_distr <- rnorm(1e6, fugacity_stats$mean, fugacity_stats$sd)
temp_distr <- units::set_units(rnorm(1e6, 300, 50 / 1.96), degC)

# Calculate strain rates using temperature, fugacity, and differential stress
edot <- creep_quartz(
  stress = stress,
  temperature = temp_distr,
  fugacity = fug_distr,
  model = "Hirth2001"
)

# Summary stats of Monte Carlo simulation
summary(edot)
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      1.2e-13 /s
#> 95% interpercentile range:   2.4e-18 - 7.2e-09 /s
#> Log-variance:                5.83321
#> Student's t-Test:            p<0.05
```

An theoretical approach to calculate the uncertainties in strain rates
uses standard error propagation and is provided through the
[`creep_quartz_analytic()`](https://tobiste.github.io/geocreep/reference/creep_quartz_analytic.md)
function:

``` r
creep_quartz_analytic(
  stress = stress,
  temperature = temp_distr,
  fugacity = fug_distr,
  model = "Hirth2001"
)
#> $e_best
#> 1.621933e-13 [1/s]
#> 
#> $sd_e_range
#> Units: [1/s]
#> [1] 2.379690e-23 1.105466e-03
#> 
#> $var_log_e_total
#> [1] 512.6828
#> 
#> $sd_log_e_total
#> [1] 22.6425
#> 
#> $var_log_e
#>          prefactor    stress_exponent  fugacity_exponent grainsize_exponent 
#>       1.908683e+00       1.582356e+01       0.000000e+00       0.000000e+00 
#>             stress           fugacity          grainsize        temperature 
#>       1.105766e+02       3.827764e+02       0.000000e+00       1.592087e+00 
#>           enthalpy 
#>       5.491489e-03
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
