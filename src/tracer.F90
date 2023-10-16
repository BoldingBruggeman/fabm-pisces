#include "fabm_driver.h"

module pisces_tracer
   use fabm_types
   use fabm_particle

   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_tracer
      type (type_state_variable_id) :: id_c, id_fe, id_si
   contains
      procedure :: initialize
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_tracer), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      logical :: has_carbon, has_nitrogen, has_phosphorus, has_silicon, has_iron

      call self%register_implemented_routines()

      call self%get_parameter(has_carbon,     'has_carbon',     '', 'tracer contains carbon',     default=.false.)
      call self%get_parameter(has_nitrogen,   'has_nitrogen',   '', 'tracer contains nitrogen',   default=.false.)
      call self%get_parameter(has_phosphorus, 'has_phosphorus', '', 'tracer contains phosphorus', default=.false.)
      call self%get_parameter(has_silicon,    'has_silicon',    '', 'tracer contains silicon',    default=.false.)
      call self%get_parameter(has_iron,       'has_iron',       '', 'tracer contains iron',       default=.false.)

      if (has_carbon .or. has_nitrogen .or. has_phosphorus) then
         call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'concentration in carbon units', minimum=0.0_rk)
         if (has_carbon)     call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c, scale_factor=1e6_rk)
         if (has_nitrogen)   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c, scale_factor=rno3 * 1e6_rk)
         if (has_phosphorus) call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      end if
      if (has_iron) then
         call self%register_state_variable(self%id_fe, 'fe', 'mol Fe L-1', 'concentration in iron units', minimum=0.0_rk)
         call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)
      end if
      if (has_silicon) then
         call self%register_state_variable(self%id_si, 'si', 'mol Si L-1', 'concentration in silicon units', minimum=0.0_rk)
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_si, scale_factor=1e6_rk)
      end if
   end subroutine initialize

end module
