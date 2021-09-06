#include "fabm_driver.h"

module pisces_nitrification
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_nitrification
      type (type_state_variable_id) :: id_no3, id_nh4, id_oxy, id_tal
      type (type_dependency_id) :: id_nitrfac, id_etot, id_gdepw_n
      type (type_surface_dependency_id) :: id_fr_i, id_emoy, id_hmld
      type (type_diagnostic_variable_id) :: id_nit, id_denit
      real(rk) :: nitrif
   contains
      procedure :: initialize
      procedure :: do
   end type

   real(rk), parameter :: rdenita =   3._rk /  5._rk

contains

   subroutine initialize(self, configunit)
      class (type_pisces_nitrification), intent(inout), target :: self
      integer,                           intent(in)            :: configunit

      call self%get_parameter(self%nitrif, 'nitrif', 'd-1', 'NH4 nitrification rate', default=0.05_rk)

      call self%register_diagnostic_variable(self%id_nit, 'nit', 'mol N L-1 s-1', 'nitrification')
      call self%register_diagnostic_variable(self%id_denit, 'denit', 'mol N L-1 s-1', 'NH4 consumed in denitrification???')

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)

      call self%register_dependency(self%id_nitrfac, 'nitrfac', '1', 'denitrification factor')
      call self%register_dependency(self%id_etot, 'etot', 'W m-2', 'instantaneous PAR averaged')
      call self%register_dependency(self%id_emoy, 'emoy', 'W m-2', 'instantaneous PAR averaged over mixing layer')
      call self%register_dependency(self%id_fr_i, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_gdepw_n, standard_variables%depth)
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_nitrification), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: no3, nh4, nitrfac, emoy, fr_i, hmld, gdepw_n
      real(rk) :: zonitr, zdenitnh4

      real(rk), parameter :: xstep = r1_rday

      _LOOP_BEGIN_
         _GET_(self%id_no3, no3)            ! nitrate (in carbon units! mol C L-1)
         _GET_(self%id_nh4, nh4)            ! ammonium (in carbon units! mol C L-1)
         _GET_(self%id_nitrfac, nitrfac)    ! oxygen limitation factor (dimensionless)
         _GET_(self%id_etot, emoy)          ! local instantaneous PAR (W m-2)
         _GET_(self%id_gdepw_n, gdepw_n)    ! depth (m)
         _GET_SURFACE_(self%id_hmld, hmld)  ! depth of mixing layer (m)
         _GET_SURFACE_(self%id_fr_i, fr_i)  ! sea ice area fraction (1)

         IF (gdepw_n <= hmld) THEN
            ! We are inside the mixing layer - replace local instantaneous PAR by PAR averaged over mixing layer
            _GET_SURFACE_(self%id_emoy, emoy)
         END IF

         ! Jorn: from p4zrem.F90
         ! NH4 nitrification to NO3. Ceased for oxygen concentrations
         ! below 2 umol/L. Inhibited at strong light 
         ! ----------------------------------------------------------
         zonitr  = self%nitrif * xstep * nh4 * ( 1.- nitrfac )  &    ! Jorn: Eq 56
         &         / ( 1.+ emoy ) * ( 1. + fr_i * emoy ) 
         zdenitnh4 = self%nitrif * xstep * nh4 * nitrfac
         zdenitnh4 = MIN(  ( no3 - rtrn ) / rdenita, zdenitnh4 ) 
         ! Update of the tracers trends
         ! ----------------------------
         _ADD_SOURCE_(self%id_nh4, - zonitr - zdenitnh4)
         _ADD_SOURCE_(self%id_no3, + zonitr - rdenita * zdenitnh4)
         _ADD_SOURCE_(self%id_oxy, - o2nit * zonitr)
         _ADD_SOURCE_(self%id_tal, - 2 * rno3 * zonitr + rno3 * ( rdenita - 1. ) * zdenitnh4)
         _SET_DIAGNOSTIC_(self%id_nit, zonitr * rno3)
         _SET_DIAGNOSTIC_(self%id_denit, zdenitnh4 * rno3)

      _LOOP_END_
   end subroutine
end module