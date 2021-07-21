#include "fabm_driver.h"

module pisces_nitrogen_fixation
   use fabm_types
   use fabm_particle
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_nitrogen_fixation
      type (type_state_variable_id) :: id_no3, id_nh4, id_po4, id_biron, id_tal, id_doc, id_poc, id_goc, id_sfe, id_bfe, id_oxy
      type (type_dependency_id) :: id_tem, id_etot_ndcy
      type (type_surface_dependency_id) :: id_fr_i
      type (type_diagnostic_variable_id) :: id_Nfix
      real(rk) :: nitrfix, diazolight, concfediaz
      real(rk) :: concnnh4, concnno3, concdnh4
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_nitrogen_fixation), intent(inout), target :: self
      integer,                               intent(in)            :: configunit

      call self%get_parameter(self%nitrfix, 'nitrfix', 'mol C L-1 d-1', 'maximum nitrogen fixation rate', default=1.e-7_rk)
      call self%get_parameter(self%diazolight, 'diazolight', 'W m-2', 'diazotroph photosynthetic parameter', default=50._rk)
      call self%get_parameter(self%concfediaz, 'concfediaz', 'mol Fe L-1', 'diazotroph half-saturation constant for iron', default=1.e-10_rk)
      call self%get_parameter(self%concnnh4, 'concnnh4', 'mol C L-1', 'ammonium half-saturation constant for phytoplankton', default=1.E-7_rk)
      call self%get_parameter(self%concnno3, 'concnno3', 'mol C L-1', 'nitrate half-saturation constant for phytoplankton', default=1.e-6_rk)
      call self%get_parameter(self%concdnh4, 'concdnh4', 'mol C L-1', 'ammonium half-saturation constant for diatoms', default=3.E-7_rk)

      call self%register_diagnostic_variable(self%id_Nfix, 'Nfix', 'mol N L-1 s-1', 'nitrogen fixation')

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phopshate')
      call self%register_state_dependency(self%id_biron, 'biron', 'mol Fe L-1', 'bioavailable iron')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_state_dependency(self%id_goc, 'goc', 'mol C L-1', 'large particulate organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')

      call self%request_coupling_to_model(self%id_doc, 'dom', 'c')
      call self%request_coupling_to_model(self%id_poc, 'pom', 'c')
      call self%request_coupling_to_model(self%id_goc, 'gom', 'c')
      call self%request_coupling_to_model(self%id_sfe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_bfe, 'gom', 'fe')

      call self%register_dependency(self%id_etot_ndcy, 'etot_ndcy', 'W m-2', 'daily mean PAR')
      call self%register_dependency(self%id_fr_i, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_tem, standard_variables%temperature)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_nitrogen_fixation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: nh4, no3, po4, biron, ztemp, etot_ndcy, fr_i, doc
      real(rk) :: zlight, zsoufer, zmudia, xdianh4, xdiano3, zlim, zfact, ztrfer, ztrpo4, ztrdp, nitrpot

      _LOOP_BEGIN_
         _GET_(self%id_nh4, nh4)
         _GET_(self%id_no3, no3)
         _GET_(self%id_po4, po4)
         _GET_(self%id_biron, biron)
         _GET_(self%id_tem, ztemp)
         _GET_(self%id_doc, doc)
         _GET_(self%id_etot_ndcy, etot_ndcy)
         _GET_SURFACE_(self%id_fr_i, fr_i)

         zlight  =  ( 1.- EXP( -etot_ndcy / self%diazolight ) ) * ( 1. - fr_i )  ! Jorn: Eq 58b, should etot_ndcy be for ice-free water?
         zsoufer = zlight * 2E-11 / ( 2E-11 + biron )

         !                      ! Potential nitrogen fixation dependant on temperature and iron
         zmudia = MAX( 0.,-0.001096*ztemp**2 + 0.057*ztemp -0.637 ) * 7.625
         !       Potential nitrogen fixation dependant on temperature and iron
         xdianh4 = nh4 / ( self%concnnh4 + nh4 )
         xdiano3 = no3 / ( self%concnno3 + no3 ) * (1. - xdianh4)
         zlim = ( 1.- xdiano3 - xdianh4 )
         IF( zlim <= 0.1 )   zlim = 0.01
         zfact = zlim   ! Jorn : dropped multiplication with rfact2 [time step in seconds]
         ztrfer = biron / ( self%concfediaz + biron )
         ztrpo4 = po4 / ( 1E-6 + po4 )
         ztrdp = ztrpo4
         nitrpot =  zmudia * r1_rday * zfact * MIN( ztrfer, ztrdp ) * zlight

         zfact = nitrpot * self%nitrfix
         _ADD_SOURCE_(self%id_nh4, + zfact / 3.0)
         _ADD_SOURCE_(self%id_tal, + rno3 * zfact / 3.0)
         _ADD_SOURCE_(self%id_po4, - zfact * 2.0 / 3.0)
         _ADD_SOURCE_(self%id_doc, + zfact * 1.0 / 3.0)
         _ADD_SOURCE_(self%id_poc, + zfact * 1.0 / 3.0 * 2.0 / 3.0)
         _ADD_SOURCE_(self%id_goc, + zfact * 1.0 / 3.0 * 1.0 / 3.0)
         _ADD_SOURCE_(self%id_oxy, + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0)
         _ADD_SOURCE_(self%id_biron, - 30E-6 * zfact * 1.0 / 3.0)
         _ADD_SOURCE_(self%id_sfe, + 30E-6 * zfact * 1.0 / 3.0 * 2.0 / 3.0)
         _ADD_SOURCE_(self%id_bfe, + 30E-6 * zfact * 1.0 / 3.0 * 1.0 / 3.0)
         _ADD_SOURCE_(self%id_biron, + 0.002 * 4E-10 * zsoufer / rday)   ! Jorn : dropped multiplication with rfact2 [time step in seconds]
         _ADD_SOURCE_(self%id_po4, + self%concdnh4 / ( self%concdnh4 + po4 )  * 0.001 * doc * xstep)   ! Jorn: ???? seems to have nothing to do with N2 fixation? adsorption to DOC? why diatom half sat?

         _SET_DIAGNOSTIC_(self%id_Nfix, zfact * rno3)
      _LOOP_END_
   end subroutine
end module