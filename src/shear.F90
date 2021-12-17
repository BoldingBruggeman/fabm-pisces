#include "fabm_driver.h"

module pisces_shear
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_shear
      type (type_diagnostic_variable_id) :: id_xdiss
      type (type_dependency_id)          :: id_gdept_n
      type (type_surface_dependency_id)  :: id_hmld
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_shear), intent(inout), target :: self
      integer,                   intent(in)            :: configunit

      call self%register_implemented_routines((/source_do/))

      call self%register_diagnostic_variable(self%id_xdiss, 'xdiss', 's-1', 'shear rate', standard_variable=shear_rate)
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_shear), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: hmld, gdept_n, xdiss

      _LOOP_BEGIN_
         _GET_SURFACE_(self%id_hmld, hmld)
         _GET_(self%id_gdept_n, gdept_n)

         ! Jorn: from p4zbio.F90, documented at end of section 4.1.1
         xdiss = 1._rk
         if (gdept_n > hmld) xdiss = 0.01_rk

         _SET_DIAGNOSTIC_(self%id_xdiss, xdiss)
      _LOOP_END_
   end subroutine

end module