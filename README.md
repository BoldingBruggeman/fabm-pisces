# FABM-PISCES

This is a [FABM](https://fabm.net) port of [the PISCES model](https://doi.org/10.5194/gmd-8-2465-2015). It is based on the PISCES code that comes with the r4.0-HEAD.r13720 version of [NEMO](https://www.nemo-ocean.eu/).

## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=pisces -DFABM_PISCES_BASE=<PISCESDIR>`

Here, `<PISCESDIR>` is the directory that with the FABM-PISCES code (the same directory that contains this readme file). Note that `-DFABM_INSTITUTES=pisces` will make FABM compile PISCES as the *only* avaialblke biogoechemical model. If you additionally want to have access to other biogeochemcial modls included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES=pisces;ersem`.

For instance, to use PISCES with the [General Ocean Turbulence Model (GOTM)](https://gotm.net/), follow [these instructions](https://github.com/fabm-model/fabm/wiki/GOTM), but in the cmake step, replace `cmake <GOTMDIR>` with:

```
cmake <GOTMDIR> -DFABM_INSTITUTES=pisces -DFABM_PISCES_BASE=<PISCESDIR>
```

## To do

* hook up POC lability parameterization (`p4zpoc.F90`) to consumption and produciton terms (currently zero)
* source of iron due to sea ice melt (`p4zsed.F90`, `ln_ironice`)
* silicate, phosphate, iron input due to dust (`p4zsed.F90`, `ln_dust`) [NB river inputs `ln_river`, nitrogen deposition `ln_ndepo` to be handled at TOP-FABM level]
* iron input from hydrothermal vents - implement in TOP-FABM or in PISCES code? (`p4zsed.F90`, `ln_hydrofe`)
* sediment denitrification (`p4zsed.F90`)
* bottom fluxes for sinking matter (`p4zsed.F90`)
* The annual maximum silicate concentration at the surface is currently set in `fabm.yaml`. FABM will need minor changes to compute it on the fly as in `p4zint.F90`.
* The turbocline depth (mixing layer depth) is currently hardcoded at 10 m in src/turbocline.F90. It should be computed dynamically, but for that, FABM and GOTM need tweaks to pass the turbulent diffusivity.
* much of the chemistry code (`p4zche.F90` in the original code) uses in-situ temperature. For the moment we substitute the native temperature provided by the host (hydrodynamic model), which often is potential or absolute temperature.
* support for an iron ligand tracer (`lk_ligand`, `p4z_ligand.F90`) is currently not implemented.
* improved initial estimate for H+ with ahini_for_at (`p4zche.F90`)
* check light fields: Are they for horizontal average of entire grid cell, or ice-free section only? What do the various processes expect?
* correct silicate dissolution (use vertical integral as in `p4zrem.F90`)

## Differences from the published PISCES description

The code refers to the equations in the [the PISCES-v2 paper](https://doi.org/10.5194/gmd-8-2465-2015) where possible. However, the PISCES version it is based on (NEMO r4.0-HEAD.r13720) has already evolved beyond the description in the paper. Thus, there are several differences between the code and the paper. We will attempt to track these here:

* diazotroph temperature dependence and NO3/NH4 limitation (`p4zsed.F90`) appear to have functionally different forms from those given in the paper
* threshold for diazotroph limitation is 0.8 in paper, 0.9 in code (`p4zsed.F90`)
* ...

## Questions to PISCES authors

* The nitrogen fixation section in `p4zsed.F90` contains non-conservative DOC-related production of PO4 that seems unrelated to nitrogen fixation, depends on a half-saturation from diatoms, and is not documented in the PISCES paper. What does this represent?
```
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + concdnh4 / ( concdnh4 + trb(ji,jj,jk,jppo4) ) &
                  &                     * 0.001 * trb(ji,jj,jk,jpdoc) * xstep
```

* Likewise, there is the following block in `p4zsed.F90`. What does this represent?
```
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
```
* Bacterial uptake of iron in `p4zrem.F90` uses a hardcoded maximum phytoplankton growth rate of 0.6 d-1 (at 0 degrees Celsius). That matches the paper, but not the latest PISCES code, in which this constant has been replaced by 0.8 d-1. Should the bacterial uptake constant not be replaced too?
*  What is the unit of zpres in `p4zche.F90`?  L248 seems to convert from dbar to bar, but comments on L345-354 suggest dbar is needed. (Asked Olivier Aumont 19 July 2021)
* why the different treatment of POC and GOC disaggregation in the mixed layer (`p4zpoc.F90`)? Is that because the sinking rate of POC is constant, whereas the sinking of GOC is (potentially) depth-dependent?