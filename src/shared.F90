#include "fabm_driver.h"

module pisces_shared

   use fabm_types

   implicit none

   public

   real(rk), parameter :: rno3    =  16._rk / 122._rk
   real(rk), parameter :: po4r    =   1._rk / 122._rk
   real(rk), parameter :: o2nit   =  32._rk / 122._rk
   real(rk), parameter :: o2ut    = 133._rk / 122._rk
   real(rk), parameter :: rdenit  =  ( ( o2ut + o2nit ) * 0.80_rk - rno3 - rno3 * 0.60_rk ) / rno3

   real(rk), parameter :: rtrn    = 0.5_rk * EPSILON( 1.e0_rk )    !: truncation value
   real(rk), parameter :: rday    = 24._rk*60._rk*60._rk      !: day                                [s]
   real(rk), parameter :: r1_rday = 1._rk / rday

   real(rk), parameter :: xstep = r1_rday

   real(rk), parameter :: nyear_len = 365._rk  ! Jorn - TODO should be made dynamic to account for leap years, different calendars, etc.

   real(rk) :: maxdt = 1800._rk ! Jorn: maximum time step in seconds, used to limit benthic-pelagic fluxes in order to prevent negative values

   type (type_universal_standard_variable), parameter :: total_chlorophyll = type_universal_standard_variable(name='total_chlorophyll', units='mg m-3', aggregate_variable=.true.)
   type (type_interior_standard_variable), parameter :: shear_rate = type_interior_standard_variable(name='shear_rate', units='s-1')
   type (type_surface_standard_variable), parameter :: mixed_layer_thickness_defined_by_vertical_tracer_diffusivity = type_surface_standard_variable(name='mixed_layer_thickness_defined_by_vertical_tracer_diffusivity', units='m')
   type (type_interior_standard_variable), parameter :: calcite_saturation_state = type_interior_standard_variable(name='calcite_saturation_state', units='-')

   type (type_universal_standard_variable), parameter :: calcite_production = type_universal_standard_variable(name='calcite_production', units='mol m-3 s-1', aggregate_variable=.true.)
end module