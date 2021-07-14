#include "fabm_driver.h"

module pisces_turbocline
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_turbocline
      type (type_surface_diagnostic_variable_id) :: id_hmld
   contains
      procedure :: initialize
      procedure :: do_column
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_turbocline), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      call self%register_diagnostic_variable(self%id_hmld, 'hmld', 'm', 'turbocline depth', &
         standard_variable=turbocline_depth, source=source_do_column)
   end subroutine initialize

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_pisces_turbocline), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      _UPWARD_LOOP_BEGIN_
         _SET_SURFACE_DIAGNOSTIC_(self%id_hmld, 10._rk)
      _UPWARD_LOOP_END_
   end subroutine

end module