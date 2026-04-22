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
#> Standard error:          4.12965
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
#> Standard error in log-space: 0.000152958
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
# Calculate strain rates using temperature, fugacity, and differential stress
edot <- creep_quartz(
  stress = stress,
  temperature = units::set_units(rnorm(1e6, 300, 50 / 1.96), degC),
  fugacity = rnorm(1e6, fugacity_stats$mean, fugacity_stats$sd),
  model = "Hirth2001"
)

# Summary stats of Monte Carlo simulation
summary(edot)
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      1.1e-13 /s
#> 95% interpercentile range:   1.8e-18 - 6.8e-09 /s
#> Standard error in log-space: 0.00242432
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
