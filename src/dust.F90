#include "fabm_driver.h"

module pisces_dust
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_dust
      type (type_surface_dependency_id)          :: id_dustdep, id_solub
      type (type_dependency_id)                  :: id_gdept_n
      type (type_diagnostic_variable_id)         :: id_zdust
      type (type_surface_diagnostic_variable_id) :: id_zirondep, id_pdust
      type (type_state_variable_id)              :: id_sil, id_po4, id_fer
      logical :: ln_solub
      real(rk) :: solub, mfrac, wdust
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_dust), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%register_implemented_routines((/source_do, source_do_surface/))

      call self%get_parameter(self%ln_solub, 'ln_solub', '', 'variable solubility of iron in dust', default=.true.)
      if (.not. self%ln_solub) call self%get_parameter(self%solub, 'solub', '1', 'solubility of iron in dust', default=0.02_rk)
      call self%get_parameter(self%mfrac, 'mfrac', '1', 'Fe mineral fraction of dust', default=0.035_rk)
      call self%get_parameter(self%wdust, 'wdust', 'm d-1', 'sinking speed of dust', default=2._rk)

      call self%register_dependency(self%id_dustdep, 'dustdep', 'g m-2 s-1', 'dust deposition')
      if (self%ln_solub) call self%register_dependency(self%id_solub, 'solub', '1', 'solubility of iron in dust')
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)

      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')

      call self%register_diagnostic_variable(self%id_zirondep, 'zirondep', 'mol m-2 s-1', 'iron deposition')
      call self%register_diagnostic_variable(self%id_pdust, 'pdust', 'g m-3', 'concentration at the surface')
      call self%register_diagnostic_variable(self%id_zdust, 'zdust', 'g m-3', 'concentration')
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_pisces_dust), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: dust, solub
      real(rk) :: zirondep, zsidep, zpdep

      real(rk), parameter :: ryyss    = nyear_len * rday    ! number of seconds per year, Jorn: nyear_len should account for leap years, calendars, etc.
      real(rk), parameter :: r1_ryyss = 1. / ryyss

      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_dustdep, dust)
         !                                              ! Iron and Si deposition at the surface
         IF( self%ln_solub ) THEN
            _GET_SURFACE_(self%id_solub, solub)
         ELSE
            solub = self%solub
         ENDIF
         zirondep = solub  * dust * self%mfrac / 55.85 + 3.e-10 * r1_ryyss   ! Jorn: dropped division by e3t_n because FABM wants flux per m2, dropped multiplication with rfact2 [time step in seconds]
         zsidep  = 8.8 * 0.075 * dust * self%mfrac / 28.1      ! Jorn: 28.1 is atomic mass of Si
         zpdep = 0.1 * 0.021 * dust * self%mfrac / 31. / po4r  ! Jorn: solubility of P in dust is assumed 10%, 31 is atomic mass of P. 0.021 * default mfrac (0.035) approximates the 750 ppm P content of dust given in the paper
         _ADD_SURFACE_FLUX_(self%id_sil, zsidep)
         _ADD_SURFACE_FLUX_(self%id_po4, zpdep)
         _ADD_SURFACE_FLUX_(self%id_fer, zirondep)
         _SET_SURFACE_DIAGNOSTIC_(self%id_zirondep, zirondep * 1.e+3)
         _SET_SURFACE_DIAGNOSTIC_(self%id_pdust, dust / ( self%wdust / rday ))
      _SURFACE_LOOP_END_
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_dust), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: dust, gdept_n, zdust, zirondep, zpdep

      _LOOP_BEGIN_
         _GET_SURFACE_(self%id_dustdep, dust)
         _GET_(self%id_gdept_n, gdept_n)

         !                                              ! Iron solubilization of particles in the water column
         !                                              ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/j
         zdust  = dust / ( self%wdust / rday ) &    ! Jorn: Eq 84, note if .not. ln_dust, dust is (coupled to) 0
         &  * EXP( -gdept_n / 540. )

         _SET_DIAGNOSTIC_(self%id_zdust, zdust)

         zirondep = zdust * self%mfrac * 0.03 / 55.85 / (270.* rday)  ! Jorn : 0.03 / 270 is presumably the dissolution rate in d-1 (it approximates the value of 0.01 % per day given in the paper)
         zpdep    = zirondep * 0.023
         !                                              ! Iron solubilization of particles in the water column
         _ADD_SOURCE_(self%id_po4, + zpdep   )          ! Jorn: TODO this should exclude the top layer
         _ADD_SOURCE_(self%id_fer, + zirondep)          ! Jorn: TODO this should exclude the top layer
      _LOOP_END_
   end subroutine

end module