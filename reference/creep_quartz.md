# Strain Rates for Creep in Quartz

Calculates strain rates of deforming quartz from stress, temperature,
and grain size from experimentally determined creep law parameters.
Monte Carlo simulation is used for propagating parameter uncertainties
in to creep estimate.

## Usage

``` r
creep_quartz(
  stress,
  temperature,
  fugacity = NULL,
  grainsize = NULL,
  pressure = NULL,
  model = c("Hirth2001", "Paterson1990", "Kronenberg1984", "Luan1992", "Gleason1995",
    "Gleason1995_melt", "Rutter2004", "Fukuda2018_LT", "Fukuda2018_HT", "Richter2018",
    "Lu2019", "Tokle2019_LT", "Tokle2019_HT", "Lusk2021_LP", "Lusk2021_HP"),
  propagate_err = TRUE,
  sim = 1e+06
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

- propagate_err:

  logical. Whether errors of the flow law parameters should be
  propagated. `TRUE` by default.

- sim:

  non-negative number. Number of Monte Carlo simulations

## Value

list. Strain rate in 1/s. If Monte Carlo Simulation was used, and object
of class `"MCS"` is returned (see
[`summary()`](https://rdrr.io/r/base/summary.html) for detailed
description of output). The flow laws produce log-normal distributed
estimates considering the uncertainties in the parameter. Hence it is
recommended to report the median (or geometric mean), and the
interpercentile range.

## Details

General flow law giving the strain rate is \$\$\dot{\varepsilon} = A
\sigma^n d^m f\_{H_2O}^r \\ e^{\left({\frac{-H}{RT}}\right)}\$\$

where \\\sigma\\ is the differential stress, \\d\\ is the grain size,
f_H_2O is the water fugacity, \\T\\ is the temperature, \\H\\ is the
enthalpy, and \\R\\ is the ideal gas constant. The flow parameters are
the prefactor \\A\\, and the exponents \\n\\, \\m\\, and \\r\\.

To propagate the uncertainties of the flow parameters Monte Carlo
simulation is used here.

- If the flow law parameters are given by a mean value and a marginal
  error (\\\mu \pm z\\), the Monte Carlo simulation assumes a normal
  distribution given by \\X = N\left(\mu, \sigma\right)\\, where \\\mu\\
  is the mean and \\\sigma\\ is the standard deviation of the mean
  (\\\sigma = \text{z}/1.96\\).

- If the parameter is given by a range of possible values
  \\\left\[x\_\text{min}, x\_\text{max}\right\]\\, the Monte Carlo
  simulation assumes an uniform distribution given by \\X =
  U\left(x\_\text{min}, x\_\text{max}\right)\\.

## References

Fukuda, J., Holyoke, C. W., & Kronenberg, A. K. (2018). Deformation of
Fine‐Grained Quartz Aggregates by Mixed Diffusion and Dislocation Creep.
Journal of Geophysical Research: Solid Earth, 123(6), 4676-4696.
[doi:10.1029/2017JB015133](https://doi.org/10.1029/2017JB015133)

Gleason, G. C., & Tullis, J. (1995). A flow law for dislocation creep of
quartz aggregates determined with the molten salt cell. Tectonophysics,
247(1-4), 1-23.
[doi:10.1016/0040-1951(95)00011-B](https://doi.org/10.1016/0040-1951%2895%2900011-B)

Hirth, G., Teyssier, C., & Dunlap, W. J. (2001). An evaluation of
quartzite flow laws based on comparisons between experimentally and
naturally deformed rocks. International Journal of Earth Sciences,
90(1), 77-87.
[doi:10.1007/s005310000152](https://doi.org/10.1007/s005310000152)

Kronenberg, A. K., & Tullis, J. (1984). Flow strengths of quartz
aggregates: Grain size and pressure effects due to hydrolytic weakening.
Journal of Geophysical Research: Solid Earth, 89(B6), 4281–4297.
[doi:10.1029/JB089iB06p04281](https://doi.org/10.1029/JB089iB06p04281)

Lu, L. X., & Jiang, D. (2019). Quartz Flow Law Revisited: The
Significance of Pressure Dependence of the Activation Enthalpy. Journal
of Geophysical Research: Solid Earth, 124(1), 241–256.
[doi:10.1029/2018JB016226](https://doi.org/10.1029/2018JB016226)

Lusk, A. D. J., Platt, J. P., & Platt, J. A. (2021). Natural and
Experimental Constraints on a Flow Law for Dislocation‐Dominated Creep
in Wet Quartz. Journal of Geophysical Research: Solid Earth, 126(5),
1-25. [doi:10.1029/2020JB021302](https://doi.org/10.1029/2020JB021302)

Paterson, M. S., & Luan, F. C. (1990). Quartzite rheology under
geological conditions. Geological Society, London, Special Publications,
54(1), 299–307.
[doi:10.1144/GSL.SP.1990.054.01.26](https://doi.org/10.1144/GSL.SP.1990.054.01.26)

Richter, B., Stünitz, H., & Heilbronner, R. (2018). The
brittle-to-viscous transition in polycrystalline quartz: An experimental
study. Journal of Structural Geology, 114(September 2017), 1-21.
[doi:10.1016/j.jsg.2018.06.005](https://doi.org/10.1016/j.jsg.2018.06.005)

Rutter, E. H., & Brodie, K. H. (2004a). Experimental grain
size-sensitive flow of hot-pressed Brazilian quartz aggregates. Journal
of Structural Geology, 26(11), 2011–2023.
[doi:10.1016/j.jsg.2004.04.006](https://doi.org/10.1016/j.jsg.2004.04.006)

Rutter, E. ., & Brodie, K. . (2004b). Experimental intracrystalline
plastic flow in hot-pressed synthetic quartzite prepared from Brazilian
quartz crystals. Journal of Structural Geology, 26(2), 259–270.
[doi:10.1016/S0191-8141(03)00096-8](https://doi.org/10.1016/S0191-8141%2803%2900096-8)

Tokle, L., Hirth, G., & Behr, W. M. (2019). Flow laws and fabric
transitions in wet quartzite. Earth and Planetary Science Letters, 505,
152-161.
[doi:10.1016/j.epsl.2018.10.017](https://doi.org/10.1016/j.epsl.2018.10.017)

## See also

[`units::set_units()`](https://r-quantities.github.io/units/reference/units.html)
to set up `units` objects;
[`summary.MCS_log()`](https://tobiste.github.io/geocreep/reference/summary-MCS.md)
for statistical parameters of Monte Carlo samples;
[`creep_quartz_analytic()`](https://tobiste.github.io/geocreep/reference/creep_quartz_analytic.md)
for an analytical solution

## Examples

``` r
set.seed(20250411)
stress <- units::set_units(100, MPa)
temperature <- units::set_units(300, degC)
pressure <- units::set_units(400, MPa)
fugacity <- ps_fugacity(pressure, temperature)

creep_quartz(
  stress = stress, temperature = temperature, fugacity = fugacity,
  model = "Hirth2001") |>
 summary()
#> Statistical summary of 1000000 Monte Carlo simulations
#> 
#> Median:                      1.2e-14 1 / s
#> 95% interpercentile range:   4.3e-19 - 3.1e-10 1 / s
#> Standard error in log-space: 0.00226041
#> Student's t-Test:            p<0.05
```
