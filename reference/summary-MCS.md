# Estimates from Monte Carlo Simulation

Estimates from Monte Carlo Simulation

## Usage

``` r
# S3 method for class 'MCS_log'
summary(object, unit = NULL, ...)

# S3 method for class 'MCS'
summary(object, unit = NULL, ...)
```

## Arguments

- object:

  numeric vector of class `"MCS"` or `"MCS_log"`. The values from n
  Monte Carlo Simulations

- unit:

  (optional) object of class `units` or `symbolic_units`, or in the case
  of `set_units` expression with symbols.

- ...:

  additional arguments affecting the summary produced.

## Value

a list. If class of `object` is `"MCS"`, the list contains the following
elements:

- `median`:

  median of the Monte Carlo simulations

- `ir.95`:

  the 95% and 68% interpercentile range

- `ir.68`:

  the 68% interpercentile range

- `mean`:

  arithmetic mean the Monte Carlo simulations

- `sd`:

  1\\\sigma\\ standard deviation about the mean

- `conf.int`:

  95% confidence intverval about the mean

- `"var"`:

  Variance

- `stderr`:

  standard error

- `t.test`:

  Statistic and p-value of the Student's t-test

- `n`:

  Number of samples

If class of `object` is `"MCS_log"`, the list contains the following
elements:

- `median`:

  median of the Monte Carlo simulations

- `ir.95`:

  the 95% and 68% interpercentile range

- `ir.68`:

  the 68% interpercentile range

- `mean`:

  geometric mean the Monte Carlo simulations

- `sd`:

  1\\\sigma\\ range about the mean

- `sd2`:

  2\\\sigma\\ range about the mean

- `conf.int`:

  95% confidence intverval about the mean

- `var.log`:

  Log-variance

- `stderr.log`:

  standard error of `log(samples)`

- `t.test`:

  Statistic and p-value of the Student's t-test of `log(samples)`

- `n`:

  Number of samples

Values will be in the unit specified by parameter `unit` or be equal to
the unit of `x` if `x` is a `units` object.

## Details

Equations of the form \\X = A b^{n \pm \sigma}\\ create non-normal,
left-skewed distributions (e.g. flow laws, and grain-size piezometers).
Thus, it is recommended to report median and percentiles instead of
mean, standard deviation and confidence intervals.

## Examples

``` r
set.seed(20250411)
MC_res <- grainsize_piezometry(12.2)
summary(MC_res)
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      92 MPa 
#> 95% interpercentile range:   49 - 190 MPa 
#> Log-variance:                0.0229959
#> Student's t-Test:            p<0.05

n <- 100
temperature <- units::set_units(rnorm(n, 300, 25), degC)
pressure <- units::set_units(rnorm(n, 400, 50), MPa)
MC_res2 <- ps_fugacity(pressure, temperature) # 37 MPa
summary(MC_res2)
#> Statistical summary of 100 Monte Carlo simulations
#> 
#> Mean:                    390 bar 
#> 95% confidence interval: 370 - 420 bar 
#> Variance:                14272.4
#> Student's t-Test:        p<0.05
```
