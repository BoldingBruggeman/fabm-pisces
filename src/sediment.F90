#include "fabm_driver.h"

module pisces_sediment
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_sediment
      type (type_bottom_diagnostic_variable_id) :: id_SedCal, id_SedSi, id_SedC, id_Sdenit
      type (type_state_variable_id)             :: id_sil, id_dic, id_tal, id_oxy, id_no3, id_nh4, id_po4, id_doc
      type (type_bottom_dependency_id)          :: id_bdepth, id_cflux, id_siflux, id_calflux
      type (type_dependency_id)                 :: id_zomegaca, id_nitrfac
      type (type_bottom_diagnostic_variable_id) :: id_bc, id_bsi, id_bcal, id_bfe
      real(rk) :: sedsilfrac = 0.03_rk
      real(rk) :: sedcalfrac = 0.6_rk
      real(rk) :: maxdt = 1800._rk ! Jorn: maximum time step in seconds, used to limit benthic-pelagic fluxes in order to prevent negative values
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_sediment), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      call self%register_diagnostic_variable(self%id_SedCal, 'SedCal', 'mol m-2 s-1',    'calcite burial')
      call self%register_diagnostic_variable(self%id_SedSi,  'SedSi',  'mol Si m-2 s-1', 'silicon burial')
      call self%register_diagnostic_variable(self%id_SedC,   'SedC',   'mol C m-2 s-1',  'organic carbon burial')
      call self%register_diagnostic_variable(self%id_Sdenit, 'Sdenit', 'mol N m-2 s-1',  'nitrate reduction')

      call self%register_dependency(self%id_cflux,  'cflux',   'mmol C m-2 s-1',  'bottom carbon flux')
      call self%register_dependency(self%id_siflux, 'siflux',  'mmol Si m-2 s-1', 'bottom silicate flux')
      call self%register_dependency(self%id_calflux,'calflux', 'mmol m-2 s-1',    'bottom calcite flux')

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')

      call self%register_dependency(self%id_bdepth, standard_variables%bottom_depth)
      call self%register_dependency(self%id_zomegaca, calcite_saturation_state)
      call self%register_dependency(self%id_nitrfac, 'nitrfac', '1', 'denitrification factor computed from O2 levels')

      ! Fake benthic pools that the pelagic modules can couple to and deposit their sinking material in.
      ! The resulting downward fluxes are picked up by this module to compute sediment processes (see request_coupling calls below).
      call self%register_diagnostic_variable(self%id_bc,   'c',   'mmol C m-2',  'carbon',   source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%register_diagnostic_variable(self%id_bsi,  'si',  'mmol Si m-2', 'silicate', source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%register_diagnostic_variable(self%id_bcal, 'cal', 'mmol m-2',    'calcite',  source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%register_diagnostic_variable(self%id_bfe,  'fe',  'mmol Fe m-2', 'iron',     source=source_constant, output=output_none, act_as_state_variable=.true.)

      call self%request_coupling(self%id_cflux,'./c_sms_tot')
      call self%request_coupling(self%id_siflux, './si_sms_tot')
      call self%request_coupling(self%id_calflux, './cal_sms_tot')
   end subroutine

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_pisces_sediment), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: gdepw_n, cflux, oxy, no3, nitrfac, zflx, zo2, zno3, zdep, zdenit2d, zbureff
      real(rk) :: zsiloss, zcaloss, zrivsil, zomegaca, excess, zfactcal, zrivalk
      real(rk) :: zrivno3, zwstpoc, zpdenit, z1pdenit, zolimit

      _BOTTOM_LOOP_BEGIN_
         _GET_BOTTOM_(self%id_bdepth, gdepw_n)
         _GET_BOTTOM_(self%id_cflux, cflux)
         _GET_(self%id_oxy, oxy)
         _GET_(self%id_no3, no3)
         _GET_(self%id_nitrfac, nitrfac)

         zflx = cflux * 1E3 / 1E4 * rday   ! Jorn: flux should be umol cm-2 d-1
         zflx  = LOG10( MAX( 1E-3, zflx ) )
         zo2   = LOG10( MAX( 10. , oxy * 1E6 ) )
         zno3  = LOG10( MAX( 1.  , no3 * 1E6 * rno3 ) )
         zdep  = LOG10( gdepw_n )
         zdenit2d = -2.2567 - 1.185 * zflx - 0.221 * zflx**2 - 0.3995 * zno3 * zo2 + 1.25 * zno3    &
            &                + 0.4721 * zo2 - 0.0996 * zdep + 0.4256 * zflx * zo2                     ! Jorn: Eq 89
         zdenit2d = 10.0**( zdenit2d )
            !
         zflx = cflux * rday   ! Jorn: flux should be mmol m-2 d-1
         zbureff = 0.013 + 0.53 * zflx**2 / ( 7.0 + zflx )**2   ! Jorn: Eq 88

         _GET_BOTTOM_(self%id_siflux, zsiloss)   ! Jorn: this is the bottom flux in mmol m-2 s-1; in PISCES it is multiplied with time step and layer thickness
         _GET_BOTTOM_(self%id_calflux, zcaloss)
         zsiloss = zsiloss * 1.E-6_rk
         zcaloss = zcaloss * 1.E-6_rk
         zrivsil = 1._rk - self%sedsilfrac
         _ADD_BOTTOM_FLUX_(self%id_sil, + zsiloss * zrivsil)
         _GET_(self%id_zomegaca, zomegaca)
         !
         excess = 1._rk - zomegaca
         zfactcal = MIN( excess, 0.2 )
         zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )   ! Jorn: Eq 91
         zrivalk  = self%sedcalfrac * zfactcal
         _ADD_BOTTOM_FLUX_(self%id_tal, + zcaloss * zrivalk * 2.0)
         _ADD_BOTTOM_FLUX_(self%id_dic, + zcaloss * zrivalk)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_SedCal, (1.0 - zrivalk) * zcaloss * 1.e+3_rk)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_SedSi, (1.0 - zrivsil) * zsiloss * 1.e+3_rk)

         zrivno3 = 1. - zbureff
         zwstpoc = cflux * 1E-6_rk
         zpdenit  = MIN( 0.5 * ( no3 - rtrn ) / rdenit / self%maxdt, zdenit2d * zwstpoc * zrivno3 )   ! Jorn: Eq 90a, added maxdt because the enitre expression is now a rate rather than an increment
         z1pdenit = zwstpoc * zrivno3 - zpdenit   ! Jorn: Eq 90b except for the o2ut scale factor
         zolimit = MIN( ( oxy - rtrn ) / o2ut / self%maxdt, z1pdenit * ( 1.- nitrfac ) )   ! Jorn: added maxdt because the enitre expression is now a rate rather than an increment
         _ADD_BOTTOM_FLUX_(self%id_doc, + z1pdenit - zolimit)
         _ADD_BOTTOM_FLUX_(self%id_po4, + zpdenit + zolimit)
         _ADD_BOTTOM_FLUX_(self%id_nh4, + zpdenit + zolimit)
         _ADD_BOTTOM_FLUX_(self%id_no3, - rdenit * zpdenit)
         _ADD_BOTTOM_FLUX_(self%id_oxy, - zolimit * o2ut)
         _ADD_BOTTOM_FLUX_(self%id_tal, + rno3 * (zolimit + (1.+rdenit) * zpdenit ))
         _ADD_BOTTOM_FLUX_(self%id_dic, + zpdenit + zolimit )
         _SET_BOTTOM_DIAGNOSTIC_(self%id_Sdenit,  rdenit * zpdenit * 1.e+3_rk * rno3)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_SedC, (1. - zrivno3) * zwstpoc * 1.e+3_rk)
      _BOTTOM_LOOP_END_
   end subroutine

end module