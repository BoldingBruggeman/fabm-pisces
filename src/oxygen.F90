#include "fabm_driver.h"

module pisces_oxygen
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_oxygen
      type (type_dependency_id) :: id_tempis, id_salinprac
      type (type_surface_dependency_id) :: id_wndm, id_fr_i, id_patm
      type (type_state_variable_id) :: id_oxy
      type (type_diagnostic_variable_id) :: id_chemo2, id_nitrfac
      type (type_surface_diagnostic_variable_id) :: id_Oflx, id_Dpo2
      real(rk) :: oxymin
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type

   real(rk), parameter ::   atcox  = 0.20946_rk         ! units atm

contains

   subroutine initialize(self, configunit)
      class (type_pisces_oxygen), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      call self%register_implemented_routines((/source_do, source_do_surface/))

      call self%get_parameter(self%oxymin, 'oxymin', 'mol O2 L-1', 'half-saturation constant for anoxia', default=1.E-6_rk)

      call self%register_state_variable(self%id_oxy, 'O2', 'mol O2 L-1', 'concentration', initial_value=177.6e-6_rk, minimum=0.0_rk)
      call self%register_diagnostic_variable(self%id_chemo2, 'chemo2', 'mol O2 (L atm)-1', 'solubility')
      call self%register_diagnostic_variable(self%id_nitrfac, 'nitrfac', '1', 'denitrication factor')
      call self%register_diagnostic_variable(self%id_Oflx, 'Oflx', 'mol m-2 s-1', 'air-sea O2 flux')
      call self%register_diagnostic_variable(self%id_Dpo2, 'Dpo2', 'uatm', 'delta pO2')

      call self%register_dependency(self%id_tempis, standard_variables%temperature) ! should be in-situ temperature (as opposed to conservative/potential)
      call self%register_dependency(self%id_salinprac, standard_variables%practical_salinity)
      call self%register_dependency(self%id_wndm, standard_variables%wind_speed)
      call self%register_dependency(self%id_fr_i, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_patm, standard_variables%surface_air_pressure)
   end subroutine initialize

   ! Jorn: from p4zche.F90
   elemental function solubility(tempis, salinprac) result(chemo2)
      real(rk), intent(in) :: tempis, salinprac
      real(rk) :: chemo2

      real(rk), parameter ::   o2atm  = 1._rk / ( 1000._rk * 0.20946_rk )
      real(rk), parameter ::   oxyco  = 1._rk / 22.4144_rk  ! converts from liters of an ideal gas to moles
                                                            ! coeff. for seawater pressure correction : millero 95

      real(rk) :: ztkel, zsal, zsal2, ztgg, ztgg2, ztgg3, ztgg4, ztgg5, zoxy

      ztkel = tempis + 273.15_rk
      zsal  = salinprac ! salinprac(ji,jj,jk) + ( 1.- tmask(ji,jj,jk) ) * 35.
      zsal2 = zsal * zsal
      ztgg  = LOG( ( 298.15_rk - tempis ) / ztkel )  ! Set the GORDON & GARCIA scaled temperature
      ztgg2 = ztgg  * ztgg
      ztgg3 = ztgg2 * ztgg
      ztgg4 = ztgg3 * ztgg
      ztgg5 = ztgg4 * ztgg

      zoxy  = 2.00856_rk + 3.22400_rk * ztgg + 3.99063_rk * ztgg2 + 4.80299_rk * ztgg3    &
      &       + 9.78188e-1_rk * ztgg4 + 1.71069_rk * ztgg5 + zsal * ( -6.24097e-3_rk   &
      &       - 6.93498e-3_rk * ztgg - 6.90358e-3_rk * ztgg2 - 4.29155e-3_rk * ztgg3 )   &
      &       - 3.11680e-7_rk * zsal2
      chemo2 = ( EXP( zoxy ) * o2atm ) * oxyco * atcox     ! mol/(L atm)
   end function

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: tempis, salinprac, oxy
      real(rk) :: nitrfac

      _LOOP_BEGIN_
         _GET_(self%id_salinprac, salinprac)
         _GET_(self%id_tempis, tempis)
         _GET_(self%id_oxy, oxy)
         _SET_DIAGNOSTIC_(self%id_chemo2, solubility(tempis, salinprac))

         ! denitrification factor computed from O2 levels - Jorn: from pzlim.F90, Eq 57, O2_min1 hardcoded to 6 umol/L (but paper gives 1 umol/L)
         nitrfac = MAX(  0.e0_rk, 0.4_rk * ( 6.e-6_rk  - oxy )    &
            &                                / ( self%oxymin + oxy )  )
         nitrfac = MIN( 1._rk, nitrfac )
         _SET_DIAGNOSTIC_(self%id_nitrfac, nitrfac)
      _LOOP_END_
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_pisces_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk), parameter ::   xconv  = 0.01_rk / 3600._rk   !: coefficients for conversion
      real(rk), parameter ::   atm_per_pa = 1._rk / 101325._rk

      real(rk) :: oxy, tempis, salinprac, wndm, fr_i, patm
      real(rk) :: chemo2, ztc, ztc2, ztc3, ztc4, zsch_o2, zws, zkgwan, zkgo2, zfld16, zflu16, zoflx

      ! Jorn: from p4zflx.F90
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_oxy, oxy)                ! oxygen concentration (mol O2/L)
         _GET_(self%id_tempis, tempis)          ! temperature (degrees Celsius)
         _GET_(self%id_salinprac, salinprac)    ! practical salinity
         _GET_SURFACE_(self%id_wndm, wndm)      ! wind speed (m/s)
         _GET_SURFACE_(self%id_fr_i, fr_i)      ! ice area fraction (1)
         _GET_SURFACE_(self%id_patm, patm)      ! atmospheric pressure (Pa)

         chemo2 = solubility(tempis, salinprac)

         ztc  = MIN( 35._rk, tempis )
         ztc2 = ztc * ztc
         ztc3 = ztc * ztc2
         ztc4 = ztc2 * ztc2

         ! Compute the schmidt Number
         zsch_o2  = 1920.4_rk - 135.6_rk  * ztc + 5.2122_rk * ztc2 - 0.109390_rk * ztc3 + 0.0009377_rk * ztc4

         !  wind speed
         zws  = wndm * wndm

         ! Compute the piston velocity for O2 and CO2
         zkgwan = 0.251_rk * zws
         zkgwan = zkgwan * xconv * ( 1._rk - fr_i )   ! oxygen equivalent of Eq 82

         zkgo2 = zkgwan * SQRT( 660._rk/ zsch_o2 )

         zfld16 = atm_per_pa * patm * chemo2 * zkgo2          ! (mol/L) * (m/s)
         zflu16 = oxy * zkgo2
         zoflx = ( zfld16 - zflu16 )

         _ADD_SURFACE_FLUX_(self%id_oxy, zoflx)
         _SET_SURFACE_DIAGNOSTIC_(self%id_Oflx, zoflx * 1000._rk)
         _SET_SURFACE_DIAGNOSTIC_(self%id_Dpo2, atcox * atm_per_pa * patm - atcox * oxy / ( chemo2 + rtrn ) )
      _SURFACE_LOOP_END_
   end subroutine

end module
