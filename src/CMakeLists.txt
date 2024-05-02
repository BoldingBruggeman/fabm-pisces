add_library(fabm_models_pisces OBJECT
            pisces_model_library.F90
            tracer.F90
            phytoplankton.F90
            zooplankton.F90
            daylength.F90
            turbocline.F90
            shear.F90
            optics.F90
            carbonate_chemistry.F90
            oxygen.F90
            shared.F90
            dom_remineralization.F90
            pom.F90
            silica_dissolution.F90
            calcite_dissolution.F90
            aggregation.F90
            nitrification.F90
            nitrogen_fixation.F90
            iron.F90
            dust.F90
            sediment.F90
           )

target_link_libraries(fabm_models_pisces PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2021 - BB
#-----------------------------------------------------------------------
