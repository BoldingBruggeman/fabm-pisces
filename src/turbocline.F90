#include "fabm_driver.h"

! Based on NEMO's OCE/ZDF/zdfmxl.F90

module pisces_turbocline
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_turbocline
      type (type_surface_diagnostic_variable_id) :: id_hmld
      type (type_dependency_id)                  :: id_gdept_n, id_avt
      type (type_bottom_dependency_id)           :: id_bdepth
      real(rk) :: avt_c, minh
   contains
      procedure :: initialize
      procedure :: do_column
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_turbocline), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      call self%get_parameter(self%avt_c, 'avt_c', 'm2 s-1', 'critical vertical diffusivity', default=5.e-4_rk)
      call self%get_parameter(self%minh, 'minh', 'm', 'minimum thickness', default=10._rk)

      call self%register_diagnostic_variable(self%id_hmld, 'hmld', 'm', 'turbocline depth', &
         standard_variable=mixed_layer_thickness_defined_by_vertical_tracer_diffusivity, source=source_do_column)
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_bdepth, standard_variables%bottom_depth)
      call self%register_dependency(self%id_avt, type_interior_standard_variable(name='vertical_tracer_diffusivity', units='m2 s-1'))
   end subroutine initialize

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_pisces_turbocline), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: avt, gdept_n, h

      _GET_BOTTOM_(self%id_bdepth, h)
      _UPWARD_LOOP_BEGIN_
         _GET_(self%id_avt, avt)
         _GET_(self%id_gdept_n, gdept_n)
         if (avt < self%avt_c .and. gdept_n > self%minh) h = gdept_n
      _UPWARD_LOOP_END_
      _SET_SURFACE_DIAGNOSTIC_(self%id_hmld, h)
   end subroutine

end module