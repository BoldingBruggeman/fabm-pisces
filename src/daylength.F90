#include "fabm_driver.h"

module pisces_daylength
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_daylength
      type (type_global_dependency_id)           :: id_nday_year
      type (type_surface_dependency_id)          :: id_gphit
      type (type_surface_diagnostic_variable_id) :: id_zstrn
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_daylength), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      call self%register_dependency(self%id_nday_year, standard_variables%number_of_days_since_start_of_the_year)
      call self%register_dependency(self%id_gphit, standard_variables%latitude)
      call self%register_diagnostic_variable(self%id_zstrn, 'zstrn', 'h', 'day length')
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_pisces_daylength), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
      
      real(rk)            :: nday_year, zrum, zcodel, gphit, zstrn, zargu
      real(rk), parameter :: rpi      = 3.141592653589793_rk             !: pi
      real(rk), parameter :: rad      = 3.141592653589793_rk / 180._rk   !: conversion from degre into radian

      ! compute the day length depending on latitude and the day
      _GET_GLOBAL_(self%id_nday_year, nday_year)
      zrum = CEILING( nday_year - 80._rk ) / nyear_len
      zcodel = ASIN(  SIN( zrum * rpi * 2._rk ) * SIN( rad * 23.5_rk )  )

      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_gphit, gphit)
         ! day length in hours
         zstrn = 0._rk
         zargu = TAN( zcodel ) * TAN( gphit * rad )
         zargu = MAX( -1._rk, MIN(  1._rk, zargu ) )
         zstrn = MAX( 0.0_rk, 24._rk - 2._rk * ACOS( zargu ) / rad / 15._rk )
         _SET_SURFACE_DIAGNOSTIC_(self%id_zstrn, zstrn)
      _SURFACE_LOOP_END_
   end subroutine

end module