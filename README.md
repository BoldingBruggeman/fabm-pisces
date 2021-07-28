# FABM-PISCES

This is a [FABM](https://fabm.net) port of [the PISCES model](https://doi.org/10.5194/gmd-8-2465-2015). It is based on the PISCES code that comes with the r4.0-HEAD.r13720 version of [NEMO](https://www.nemo-ocean.eu/).

## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=pisces -DFABM_PISCES_BASE=<PISCESDIR>`

Here, `<PISCESDIR>` is the directory that with the FABM-PISCES code (the same directory that contains this readme file). Note that `-DFABM_INSTITUTES=pisces` will make FABM compile PISCES as the *only* available biogeochemical model. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="pisces;ersem"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

For instance, to use PISCES with the latest stable release of the [General Ocean Turbulence Model (GOTM)](https://gotm.net/), do the following:

```
git clone --recurse-submodules -b v6.0 https://github.com/gotm-model/code.git gotm
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/BoldingBruggeman/fabm-pisces.git
mkdir build
cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_INSTITUTES=pisces -DFABM_PISCES_BASE=../fabm-pisces
make install
```

This will install the GOTM executable with support for PISCES at `~/local/gotm/bin/gotm`.

## How to run a FABM-PISCES simulation

A `fabm.yaml` file with the PISCES configuration is provided under `<PISCESDIR>/testcases`. You can drop this file in the working directory of a FABM-compatible model such as GOTM to use it during simulation. Note that for GOTM, you will also need to ensure that `fabm/use` is set to `true` in `gotm.yaml`. Otherwise GOTM would run with physics only.

## To do

* hook up POC lability parameterization (`p4zpoc.F90`) to consumption and production terms (currently zero)
* correct silicate dissolution (use vertical integral as in `p4zrem.F90`)
* use improved initial estimate for H+ with ahini_for_at (`p4zche.F90`)
* add source of iron due to sea ice melt (`p4zsed.F90`, `ln_ironice`)
* add iron input from hydrothermal vents - implement in TOP-FABM or in PISCES code? (`p4zsed.F90`, `ln_hydrofe`)
* The annual maximum silicate concentration at the surface is currently set in `fabm.yaml`. FABM will need minor changes to compute it on the fly as in `p4zint.F90`.
* much of the chemistry code (`p4zche.F90` in the original code) uses in-situ temperature. For the moment we substitute the native temperature provided by the host (hydrodynamic model), which often is potential or conservative temperature.
* support for an iron ligand tracer (`lk_ligand`, `p4z_ligand.F90`) is currently not implemented.
* check light fields: Are they for horizontal average of entire grid cell, or ice-free section only? What do the various processes expect?
* check dust inputs for top layer (why no dissolution at the center in the original PISCES code?)

## Differences from the published PISCES description

The code refers to the equations in the [the PISCES-v2 paper](https://doi.org/10.5194/gmd-8-2465-2015) where possible. However, the PISCES version it is based on (NEMO r4.0-HEAD.r13720) has already evolved beyond the description in the paper. Thus, there are several differences between the code and the paper. We will attempt to track these here:

* diazotroph temperature dependence and NO3/NH4 limitation (`p4zsed.F90`) appear to have functionally different forms from those given in the paper
* threshold for diazotroph limitation is 0.8 in paper, 0.9 in code (`p4zsed.F90`)
* phytoplankton maximum growth rate [its reference value at 0 degrees Celsius] is 0.6 d-1 in paper, 0.8 d-1 in code.
* particulate organic mater (POC, GOC) remineralisation is now depth-dependent and calculated using a configurable number of lability classes (`p4zpoc.F90`). This has replaced Eq 38.
* exponent for silicate dissolution (`p4zrem.F90`) is now 9.25 (9 in Eq 52 of paper)
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
* in `p4zsed.F90`, what does the 270 represent in `zwdust = 0.03 * rday / ( wdust * 55.85 ) / ( 270. * rday )`?
* dust dissolution (`p4zsed.F90`) has two parts: instantaneous dissolution upon deposition at the surface, and (slower) dissolution throughout the column. But why does the latter not apply to the (center of) the top layer?
* in `p4zsed.F90`, calcite dissolution in sediment depends on the calcite saturation state of the overlying water (Eq 91 in paper): 
```
               zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
               zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
```
First question: does `zfactcal` represent the preserved fraction, as hinted at by the paper and matching the fact that in increases with increasing saturation [= decreasing `zfactcal`]. Or does it represent the dissolved fraction, as it is treated inthe original PISCES code? In the latter case, however, `excess` takes negative values when the water is supersaturated. This leads to `zfactcal` being positive, and dissolution occuring. Is that intentional? It leads to further oversaturation and even faster dissolution - a positive feedback.