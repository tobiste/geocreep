
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geocreep

<!-- badges: start -->

<!-- badges: end -->

The goal of geocreep is to quick and easy calculate strain rates, water
fugacities and deviatoric stress from grain sizes using R. The codes use
Monte Carlo simulations to propagate uncertainites of the flow
parameters into the final estimate.

## Installation

You can install the development version of geocreep from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("tobiste/geocreep")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(geocreep)

# Set a sedd for reproducibility
set.seed(20250411)
```

### Water fugacity

Calculate water fugacity from temperature and pressures using the Pitzer
and Sterner (1994) equation:

``` r
temperature <- units::set_units(300, degC)
pressure <- units::set_units(400, MPa)
fugacity <- ps_fugacity(pressure, temperature)
print(fugacity)
#> 371.3371 [bar]
```

### Grain-size piezometry

Calculating deviatoric stress from grain size (using the Stipp and
Tullis perometer):

``` r
grainsize <- units::set_units(11, um)
stress <- grainsize_piezometry(grainsize, method = "Stipp-reg2-3")
print(stress$median)
#> 99.79369 [MPa]
```

### Flow laws

Calculate strain rate using a defined flow law from stress, temperature,
and fugacity:

``` r
edot <- strain_rate(stress = stress$median, temperature = temperature, fugacity = fugacity, model = "Hirth2001")
print(edot$median)
#> 1.152363e-14 [1/s]
```
