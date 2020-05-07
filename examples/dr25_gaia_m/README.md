# Required Files:
To reproduce a calculation from Hsu et al. (2019), you will need the following data files (under the 'data' sub-directory of the repository):
* "DR25topwinfuncs.jld": Summary of window function data taken from the DR25 Completeness Products
* "dr25m_osds.jld": One-sigma depth functions for the filtered catalog of Kepler planet search targets.
* "KeplerMAST_TargetProperties.csv": Summary of key properties for Kepler Targets that are not included in the other catalogs.
* "q1q17_dr25_gaia_m.jld": Filtered DR25 stellar catalog of M targets (see Hsu et al. 2020 for explanation of cuts).
* "q1_q17_dr25_koi.csv": DR25 KOI catalog with transit depths replaced with median depths from the MCMC posterior chains.

# Running Calculation:
To perform the actual calculations, you run:
```
> julia abc_run.jl
```

This defaults to computing a planet candidate occurrence rate and rate density for 7 period-radius bins spanning 0.25-4 R_Earth and 8-16 days.
For other ranges, you'd edit the values of p_bin_lim and r_bin_lim in param.in.