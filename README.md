# FABM-PISCES

This is a [FABM](https://fabm.net) port of [the PISCES model](https://doi.org/10.5194/gmd-8-2465-2015). It is based on the PISCES code that comes with the r4.0-HEAD.r13720 version of [NEMO](https://www.nemo-ocean.eu/).

The code has been modularised to the point where it is straightforward to create configurations with any number
of phytoplankton and zooplankton types and any number of particulate organic matter classes - all by adjusting
the runtime configuration (`fabm.yaml`), no code change or recompilation needed. However, some more work would
be needed to fully support such configurations. Specifically, the zooplankton code would need to be changed to
handle a runtime-configurable number of prey types.

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

* add source of iron due to sea ice melt (`p4zsed.F90`, `ln_ironice`)
* add iron input from hydrothermal vents - implement in TOP-FABM or in PISCES code? (`p4zsed.F90`, `ln_hydrofe`)
* enable iron input fom sediments throughout the column using coast and island mask (now restricted to bottom only)
* much of the chemistry code (`p4zche.F90` in the original code) uses in-situ temperature. For the moment we substitute the native temperature provided by the host (hydrodynamic model), which often is potential or conservative temperature.
* support for an iron ligand tracer (`lk_ligand`, `p4z_ligand.F90`) is currently not implemented [also not used in default PISCES configuration]
* check light fields: Are they for horizontal average of entire grid cell, or ice-free section only? What do the various processes expect?
* check dust inputs for top layer (why no dissolution at the center in the original PISCES code?)
* depth-dependent sinking velocity of large POM not implemented yet [but also not used in default parameterization]

## Differences from the published PISCES description

The code refers to the equations in the [the PISCES-v2 paper](https://doi.org/10.5194/gmd-8-2465-2015) where possible. However, the PISCES version it is based on (NEMO r4.0-HEAD.r13720) has already evolved beyond the description in the paper. Thus, there are several differences between the code and the paper. We will attempt to track these here:

* phytoplankton maximum growth rate [its reference value at 0 degrees Celsius] is 0.6 d-1 in paper, 0.8 d-1 in code.
* particulate organic mater (POC, GOC) remineralisation is now depth-dependent and calculated using a configurable number of lability classes (`p4zpoc.F90`). This has replaced Eq 38.
* exponent for silicate dissolution (`p4zrem.F90`) is 9 in Eq 52 of paper, but 9.25 in code. The latter matches the [Ridgwell paper](https://doi.org/10.1029/2002GB001877).
* quadratic mortality of diatoms depends in a different way on nutrient limitation than described in the paper (Eq 13)
* quadratic and linear mortality of nanophytoplankton is now modified by cell size
* a threshold concentration `xmort` of diatoms and nanophytoplankton is now protected from mortality
* the P-I slope of diatoms is now the average of separate slopes of small and large cells (and the small : large ratio is dynamic, inferred from model state)
* the P-I slope of nanophytoplankton is now optionally temperature-dependent through the `xadap` parameter (but this is off by default as `xadap` defaults to 0)
* only the new formulation of the P-I slope is supported in the code. In effect the `ln_newprod` switch mentioned in the paper is always on, Eq 2a is always used, Eq 2b never.
* The time spent below the euphotic zone is now incorporated in a reduced effective day length. This makes Eqs 3b-3d (and any use of f2, as e.g. for chlorophyll synthesis) in paper redundant. Also, only one term (the equivalent of f1(L_day)) is now used in the photosynthesis equation, while in the paper f1(L_day) is used for maximum photosynthesis, but L_day within the light limitation term.
* maximum chl : C ratio is now temperature dependent
* chlorophyll synthesis (Eq 15b) now lacks division by nutrient limitation in light limitation term. This is similar to the difference between Eq 2a and 2b (has chlorophyll been "upgraded" to become like Eq 2a?)
* the upper limit of 5.4 in equation 22 for silicate uptake seem to have been removed
* latitude dependence of silicate limitation differs (paper: threshold at equator, above limiter is 1; code: threshold at 30 degrees South, above limiter is 1 + Si3/(Si3+c))
* zooplankton: a threshold concentration of 1e-9 mol C/L is protected from mortality (linear and quadratic)
* quadratic zooplankton mortality and grazing are now depressed when oxygen becomes limiting
* zooplankton: new expressions for decrease in gross growth efficiency due to low food quality or high food quantity [latter off by default]
* dom/bacteria: the specific ammonification rate is now clipped to equal or exceed 2.74e-4 d-1 (Eq 33a, 33b)
* Brownian aggregation of DOC is controlled by a rate constant a4=5095 in the paper, but in the code this constant is now 2.4.
* Wherever DOC features in agggregation processes, it is now multiplied with 0.3
* Values for aggregation parameters a8 and a9 seem to have been swapped
* GOC-DOC Brownian aggregation is new (in code but not in Eq 36c)
* Values for the two oxygen thresholds (1 and 6 umol/L) in the oxygen limitation factor denitrification (Eq 57) have been swapped, but that seems to have been a typo in the paper since it states denitrifcation should stop at 6 umol/L (ok in code)
* diazotroph temperature dependence and NO3/NH4 limitation (`p4zsed.F90`) appear to have functionally different forms from those given in the paper
* threshold for diazotroph limitation is 0.8 in paper, 0.9 in code (`p4zsed.F90`)
* iron scavenging due to aggregation of colloids now misses the Brrownion POC-colloid interaction 9the term with a4 in Eq 61a)
* iron scavenging due to precipitation is now modelled as function of the concentration of free inorganic iron and its solubility. This replaces Eq 62
* all mortality of the calcifying fraction of nanophytoplankton now produces calcite - in the paper that was half (0.5 in Eq 76)
* the calcifying fraction (rain ratio) of nanophytoplankton is now constrained between 0.02 and 0.8
* dust dissolution is now a source of PO4 throughout the column (surface only in the paper)
* the dust concentration now decays exponentially with depth. Likely that is to account for dissolution while it sinks.
* ...

## Questions to PISCES authors

In many cases, answers to these questions would help decide how to implement particular functionality in the FABM port of PISCES.
In some cases, the questions could potentially point to minor issues in the PISCES code. It is more likely, however, that they reflect
our lack of understandng of what that original code does.

*  What is the unit of `zpres` in `p4zche.F90`?  L248 seems to convert from dbar to bar, but comments on L345-354 suggest dbar is needed.
* in `p4zsed.F90`, calcite dissolution in sediment depends on the calcite saturation state of the overlying water (Eq 91 in paper): 
```
               zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
               zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
```
Does `zfactcal` represent the preserved fraction, as hinted at by the paper and matching the fact that it increases with increasing saturation [= decreasing `zfactcal`]? Or does it represent the dissolved fraction, as it is treated in the original PISCES code? In the latter case, however, `excess` takes negative values when the water is supersaturated. This leads to `zfactcal` being positive and dissolution occuring. Is that intentional? It leads to further oversaturation and even faster dissolution - a positive feedback.
* The nitrogen fixation section in `p4zsed.F90` contains non-conservative DOC-related production of PO4 that seems unrelated to nitrogen fixation, depends on a half-saturation from diatoms, and is not documented in the PISCES paper. What does this represent?
```
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + concdnh4 / ( concdnh4 + trb(ji,jj,jk,jppo4) ) &
                  &                     * 0.001 * trb(ji,jj,jk,jpdoc) * xstep
```

* There is non-conservative iron production in `p4zsed.F90`, linked to diazotrophs via `zsoufer`. What process does this represent?
```
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
```
* Bacterial uptake of iron in `p4zrem.F90` uses a hardcoded maximum phytoplankton growth rate of 0.6 d-1 (at 0 degrees Celsius). That matches the paper, but not the latest PISCES code, in which this constant has been replaced by 0.8 d-1. Should the bacterial uptake constant not be replaced too?
* why the different treatment of POC and GOC remineralisation (constant for GOC, but dynamically calculated for POC) in the mixed layer (`p4zpoc.F90`)? 
* about pom remineralisation (`p4zpoc.F90`): why is the POC production due to GOC disaggregation proportional to `zorem3(ji,jj,jk) * alphag(ji,jj,jk,jn)` instead of `reminp(jn) * alphag(ji,jj,jk,jn)`? (in other words: why does it depend on the *mean* remineraliation rate instead of on the lability-class-specific rate?). And why is this GOC->POC conversion not taken into acccount within the mixed layer as additional POC production term? Why is the specific consumption rate in the mixed layer computed by first calculating the specific rate and then depth-integrating it, instead of computing the integral of consumption and dividing that by the integral of POC concentration?
* dust dissolution (`p4zsed.F90`) has two parts: instantaneous dissolution upon deposition at the surface, and (slower) dissolution throughout the column. But why does the latter not apply to the (center of) the top layer? (the responsible loop starts at `jk=2`) That introduces a grid scale dependence of the model solution (it depends on the thickness of the top layer)
* Assuming `dust` is in g m-2 s-1 and `wdust` in m d-1, shouldn't `wdust` be divided by `rday` instead of multiplied by it in `CALL iom_put( "pdust"  , dust(:,:) / ( wdust * rday )  * tmask(:,:,1) ) ` (in `p4zsed.F90`)
* quadratic mortality losses of mesozooplankton are ascribed to higher trophic levels. Accordingly a fraction of this loss is respired, a fraction ends up in fecal pellets. The remainder supposedly inside those HTLs - but since they are not tracked by PISCES, that is a non-conservative loss terms. What's the idea behind this? (where does that mass ultimately go? harvested?)
