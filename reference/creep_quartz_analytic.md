# Algebraic error propagation of flow laws

Uses 1st order Taylor expansion (linearization) of the natural log-space
flow law

## Usage

``` r
creep_quartz_analytic(
  stress,
  temperature,
  fugacity = NULL,
  grainsize = NULL,
  pressure = NULL,
  model = c("Hirth2001", "Paterson1990", "Kronenberg1984", "Luan1992", "Gleason1995",
    "Gleason1995_melt", "Rutter2004", "Fukuda2018_LT", "Fukuda2018_HT", "Richter2018",
    "Lu2019", "Tokle2019_LT", "Tokle2019_HT", "Lusk2021_LP", "Lusk2021_HP")
)
```

## Arguments

- stress:

  Differential stress in MPa or `units` object

- temperature:

  Temperature in Kelvin or `units` object

- fugacity:

  Water fugacity in MPa or `units` object

- grainsize:

  Grain size in \\\mu\\m or `units` object

- pressure:

  Pressure in MPa or `units` object

- model:

  character specifying the flow law to be used:

  `"Hirth2001"`

  :   Hirth and Tullis (2001), dislocation creep

  `"Paterson1990"`

  :   Paterson and Luan (1990): dislocation creep; axial compression

  `"Kronenberg1984"`

  :   Kronenberg and Tullis (1984): deformation mechanism: dislocation
      creep and grain-size sensitive creep; strain geometry: axial
      compression

  `"Luan1992"`

  :   Luan and Paterson (1990): dislocation creep; axial compression

  `"Gleason1995"`

  :   Gleason and Tullis (1995): dislocation creep; axial compression

  `"Gleason1995_melt"`

  :   Gleason and Tullis (1995): dislocation creep; axial compression;
      1–2% melt

  `"Rutter2004"`

  :   Rutter and Brodie (2004b): dislocation creep; axial compression

  `"Fukuda2018_LT"`

  :   Fukada et al. (2018): dislocation creep; axial compression; low
      temperatures (600–750 °C)

  `"Fukuda2018_HT"`

  :   Fukada et al. (2018): dislocation creep and grain-size sensitive
      creep; axial compression; high temperatures (800–950 °C)

  `"Richter2018"`

  :   Richter et al. (2018): dislocation creep and grain-size sensitive
      creep; general shear; 800–1000 °C

  `"Lu2019"`

  :   Lu and Jiang (2019): dislocation creep

  `"Tokle2019_HT"`

  :   Tokle et al. (2019): dislocation creep and grain-size sensitive
      creep; high temperatures/low stress

  `"Tokle2019_LT"`

  :   Tokle et al. (2019): dislocation creep and grain-size sensitive
      creep; low temperature/high stress

  `"Lusk2021_LP"`

  :   Lusk et al. (2021): dislocation-dominated creep in wet quartz, for
      low pressures (≤560 MPa)

  `"Lusk2021_HP"`

  :   Lusk et al. (2021): dislocation-dominated creep in wet quartz, for
      high pressures (700–1600 MPa)

## Value

list.

- `"e_best"`:

  Strain rate (1/s)

- `"sd_e_range"`:

  Range of 1s standard deviation (1/s)

- `"var_log_e_total"`:

  Variance in natural log-space

- `"sd_log_e_total`:

  Standard deviation in natural log-space

- `"var_log_e"`:

  Individual variance components in natural log-space

## See also

[`creep_quartz()`](https://tobiste.github.io/geocreep/reference/creep_quartz.md)

## Examples

``` r
stress <- units::set_units(100, MPa)
temperature <- units::set_units(300, degC)
pressure <- units::set_units(400, MPa)
fugacity <- ps_fugacity(pressure, temperature)

creep_quartz_analytic(
  stress = stress, temperature = temperature, fugacity = fugacity, pressure = pressure,
  model = "Hirth2001"
)
#> $e_best
#> 1.165839e-14 [1/s]
#> 
#> $sd_e_range
#> Units: [1/s]
#> [1] 1.834746e-16 7.408012e-13
#> 
#> $var_log_e_total
#> [1] 17.23666
#> 
#> $sd_log_e_total
#> [1] 4.151706
#> 
#> $var_log_e
#>          prefactor    stress_exponent  fugacity_exponent grainsize_exponent 
#>        1.908683320       15.322485539        0.000000000        0.000000000 
#>             stress           fugacity          grainsize        temperature 
#>        0.000000000        0.000000000        0.000000000        0.000000000 
#>           enthalpy 
#>        0.005491873 
#> 

creep_quartz_analytic(
  stress = stress, temperature = temperature, fugacity = fugacity, pressure = pressure,
  model = "Lusk2021_LP"
)
#> $e_best
#> 4.17758e-13 [1/s]
#> 
#> $sd_e_range
#> Units: [1/s]
#> [1] 4.982208e-14 3.502899e-12
#> 
#> $var_log_e_total
#> [1] 4.521764
#> 
#> $sd_log_e_total
#> [1] 2.126444
#> 
#> $var_log_e
#>          prefactor    stress_exponent  fugacity_exponent grainsize_exponent 
#>          2.3095068          0.8483037          0.2207950          0.0000000 
#>             stress           fugacity          grainsize        temperature 
#>          0.0000000          0.0000000          0.0000000          0.0000000 
#>           enthalpy 
#>          1.1431584 
#> 

creep_quartz_analytic(
  stress = stress, temperature = temperature, fugacity = fugacity, pressure = pressure,
  model = "Kronenberg1984"
)
#> $e_best
#> 6.264269e-07 [1/s]
#> 
#> $sd_e_range
#> Units: [1/s]
#> [1] 5.288804e-07 7.419650e-07
#> 
#> $var_log_e_total
#> [1] 0.02865229
#> 
#> $sd_log_e_total
#> [1] 0.1692699
#> 
#> $var_log_e
#>          prefactor    stress_exponent  fugacity_exponent grainsize_exponent 
#>        0.000000000        0.001192927        0.000000000        0.000000000 
#>             stress           fugacity          grainsize        temperature 
#>        0.000000000        0.000000000        0.000000000        0.000000000 
#>           enthalpy 
#>        0.027459365 
#> 
```
