#include "fabm_driver.h"

module pisces_dom_remineralization
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_dom_remineralization
      type (type_state_variable_id) :: id_no3, id_nh4, id_po4, id_fer, id_doc, id_oxy, id_dic, id_tal, id_sfe, id_bfe
      type (type_dependency_id) :: id_zoo, id_mes, id_gdept_n, id_tem
      type (type_surface_dependency_id) :: id_hmld, id_heup
      type (type_dependency_id) :: id_zdepbac, id_nitrfac
      type (type_diagnostic_variable_id) :: id_remin, id_denit
      type (type_diagnostic_variable_id) :: id_febact, id_blim

      real(rk) :: xremik, xkdoc, concno3, concnh4, concfe
      real(rk) :: feratb, xkferb
   contains
      procedure :: initialize
      procedure :: do
   end type

   type, extends(type_base_model), public :: type_pisces_bacteria
      type (type_dependency_id) :: id_gdept_n, id_zoo, id_mes
      type (type_surface_dependency_id) :: id_hmld, id_heup
      type (type_diagnostic_variable_id) :: id_bact
   contains
      procedure :: initialize => bacteria_initialize
      procedure :: do_column  => bacteria_do_column
   end type

   real(rk), parameter :: rdenit  =  ( ( o2ut + o2nit ) * 0.80 - rno3 - rno3 * 0.60 ) / rno3

contains

   subroutine initialize(self, configunit)
      class (type_pisces_dom_remineralization), intent(inout), target :: self
      integer,                                  intent(in)            :: configunit

      class (type_pisces_bacteria), pointer :: pbacteria

      call self%get_parameter(self%xremik, 'xremik', 'd-1', 'remineralization rate', default=0.3_rk)
      call self%get_parameter(self%xkdoc, 'xkdoc', 'mol C L-1', 'DOC half-saturation constant', default=417.E-6_rk)
      call self%get_parameter(self%concno3, 'concno3', 'mol C L-1', 'nitrate half-saturation constant', default=2.E-7_rk)   ! ???? 0.03 umol N/L in paper
      call self%get_parameter(self%concnh4, 'concnh4', 'mol C L-1', 'ammonium/phosphate half-saturation constant', default=2.E-8_rk)! ???? 0.003 umol N or P/L in paper
      call self%get_parameter(self%concfe, 'concfe', 'mol Fe L-1', 'iron half-saturation constant', default=1.E-11_rk)

      call self%get_parameter(self%feratb, 'feratb', 'mol Fe (mol C)-1', 'Fe/C quota in bacteria', default=10.E-6_rk)
      call self%get_parameter(self%xkferb, 'xkferb', 'mol Fe L-1', 'half-saturation constant for bacteria Fe/C', default=3.E-10_rk)

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)

      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_heup, 'heup', 'm', 'euphotic depth')
      call self%register_dependency(self%id_hmld, turbocline_depth)
      call self%register_dependency(self%id_zoo, 'zoo', 'mol C L-1', 'microzooplankton')
      call self%register_dependency(self%id_mes, 'mes', 'mol C L-1', 'mesozooplankton')
      call self%register_dependency(self%id_tem, standard_variables%temperature)

      allocate(pbacteria)
      call self%add_child(pbacteria, 'bacteria', configunit=-1)
      call pbacteria%request_coupling(pbacteria%id_heup, '../heup')
      call pbacteria%request_coupling(pbacteria%id_zoo, '../zoo')
      call pbacteria%request_coupling(pbacteria%id_mes, '../mes')

      call self%register_dependency(self%id_zdepbac, 'zdepbac', 'mol C L-1', 'bacterial biomass proxy')
      call self%request_coupling(self%id_zdepbac, './bacteria/bact')
      call self%register_dependency(self%id_nitrfac, 'nitrfac', '1', 'denitrification factor computed from O2 levels')

      call self%register_diagnostic_variable(self%id_remin, 'remin', 'mol N L-1 s-1', 'remineralization')
      call self%register_diagnostic_variable(self%id_denit, 'denit', 'mol N L-1 s-1', 'denitrification')
      call self%register_diagnostic_variable(self%id_febact, 'febact', 'mol Fe L-1 s-1', 'bacterial uptake of Fe')
      call self%register_diagnostic_variable(self%id_blim, 'blim', 'umol C L-1', 'effective bacterial biomass for iron uptake')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_dom_remineralization), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: no3, nh4, po4, fer, doc, oxy
      real(rk) :: hmld, heup, gdept_n, tem
      real(rk) :: zdenom, xnanono3, xnanonh4
      real(rk) :: zlim1, zlim2, zlim3, zlim4, xlimbacl, xlimbac
      real(rk) :: zdepbac, nitrfac, nitrfac2
      real(rk) :: zremik, zolimit, zolimi, zammonic, denitr, zoxyremc
      real(rk) :: tgfunc, zdep, zdepmin, zdepprod, zdepeff, zbactfer, blim

      real(rk), parameter :: xstep = r1_rday

      _LOOP_BEGIN_
         _GET_(self%id_no3, no3)
         _GET_(self%id_nh4, nh4)
         _GET_(self%id_po4, po4)
         _GET_(self%id_fer, fer)
         _GET_(self%id_doc, doc)
         _GET_(self%id_oxy, oxy)
         _GET_(self%id_zdepbac, zdepbac)
         _GET_(self%id_nitrfac, nitrfac)
         _GET_SURFACE_(self%id_hmld, hmld)
         _GET_SURFACE_(self%id_heup, heup)
         _GET_(self%id_gdept_n, gdept_n)
         _GET_(self%id_tem, tem)                  ! temperature (degrees Celsius)

         ! Jorn: from p4zlim.F90
         zdenom = 1._rk /  ( self%concno3 * self%concnh4 + self%concnh4 * no3 + self%concno3 * nh4 )  ! Jorn: denominator in Eq 34g,h
         xnanono3 = no3 * self%concnh4 * zdenom                              ! Jorn: Eq 34h
         xnanonh4 = nh4 * self%concno3 * zdenom                              ! Jorn: Eq 34g
         !
         zlim1    = xnanono3 + xnanonh4
         zlim2    = po4 / ( po4 + self%concnh4 )         ! Jorn: Eq 34e
         zlim3    = fer / ( self%concfe + fer )          ! Jorn: Eq 34d
         zlim4    = doc / ( self%xkdoc   + doc )         ! Jorn: Eq 34b
         xlimbacl = MIN( zlim1, zlim2, zlim3 )           ! Jorn: Eq 34c
         xlimbac  = MIN( zlim1, zlim2, zlim3 ) * zlim4   ! Jorn: Eq 34a

         ! denitrification factor computed from NO3 levels
         nitrfac2 = MAX( 0.e0_rk,       ( 1.E-6_rk - no3 )  &
            &                                / ( 1.E-6_rk + no3 ) )
         nitrfac2 = MIN( 1., nitrfac2 )

         ! Jorn: from p4zrem.F90
         ! DOC ammonification. Depends on depth, phytoplankton biomass
         ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 

         zremik = self%xremik * xstep / 1.e-6_rk * xlimbac * zdepbac   ! Jorn: 1.e-6_rk likely takes role of Bact_ref
         zremik = MAX( zremik, 2.74e-4_rk * xstep )
         ! Ammonification in oxic waters with oxygen consumption
         ! -----------------------------------------------------
         zolimit = zremik * ( 1._rk- nitrfac ) * doc 
         zolimi = MIN( ( oxy - rtrn ) / o2ut, zolimit )        ! Jorn: Eq 33a
         ! Ammonification in suboxic waters with denitrification
         ! -------------------------------------------------------
         zammonic = zremik * nitrfac * doc
         denitr  = zammonic * ( 1._rk - nitrfac2 )
         denitr  = MIN( ( no3 - rtrn ) / rdenit, denitr )      ! Jorn: Eq 33b
         zoxyremc          = zammonic - denitr
         !
         zolimi = MAX( 0.e0_rk, zolimi )
         denitr = MAX( 0.e0_rk, denitr )
         zoxyremc  = MAX( 0.e0_rk, zoxyremc )

         _ADD_SOURCE_(self%id_po4, + zolimi + denitr + zoxyremc)
         _ADD_SOURCE_(self%id_nh4, + zolimi + denitr + zoxyremc)
         _ADD_SOURCE_(self%id_no3, - denitr * rdenit)
         _ADD_SOURCE_(self%id_doc, - zolimi - denitr - zoxyremc)
         _ADD_SOURCE_(self%id_oxy, - zolimi * o2ut)
         _ADD_SOURCE_(self%id_dic, + zolimi + denitr + zoxyremc)
         _ADD_SOURCE_(self%id_tal, + rno3 * ( zolimi + zoxyremc    &
         &                     + ( rdenit + 1._rk) * denitr ))

         _SET_DIAGNOSTIC_(self%id_remin, zolimi)
         _SET_DIAGNOSTIC_(self%id_denit, denitr * rdenit * rno3)

         zdep = MAX( hmld, heup )
         zdepmin = MIN( 1._rk, zdep / gdept_n )
         zdepprod = zdepmin**0.273_rk
         zdepeff  = 0.3_rk * zdepmin**0.3_rk
         tgfunc = EXP( 0.063913_rk * tem )  ! Jorn: Eq 4a in PISCES-v2 paper, NB EXP(0.063913) = 1.066 = b_P

         ! Bacterial uptake of iron. No iron is available in DOC. So
         ! Bacteries are obliged to take up iron from the water. Some
         ! studies (especially at Papa) have shown this uptake to be significant
         ! ----------------------------------------------------------
         zbactfer = self%feratb * 0.6_rk / rday * tgfunc * xlimbacl     & ! Jorn: Eq 63, mu_P hardcoded to 0.6 (as in paper), but in latest code has evolved to 0.8. Dropped multiplication with rfact2 [time step in seconds]
            &              * fer / ( self%xkferb + fer )    &
            &              * zdepprod * zdepeff * zdepbac
         _ADD_SOURCE_(self%id_fer, - zbactfer*0.33_rk)
         _ADD_SOURCE_(self%id_sfe, + zbactfer*0.25_rk)
         _ADD_SOURCE_(self%id_bfe, + zbactfer*0.08_rk)
         blim      = xlimbacl  * zdepbac / 1.e-6_rk * zdepprod
         _SET_DIAGNOSTIC_(self%id_febact, zbactfer * 0.33_rk)
         _SET_DIAGNOSTIC_(self%id_blim, blim)

      _LOOP_END_
   end subroutine

   subroutine bacteria_initialize(self, configunit)
      class (type_pisces_bacteria), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      call self%register_diagnostic_variable(self%id_bact, 'bact', 'mol C L-1', 'biomass proxy', source=source_do_column)

      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_heup, 'heup', 'm', 'euphotic depth')
      call self%register_dependency(self%id_hmld, turbocline_depth)
      call self%register_dependency(self%id_zoo, 'zoo', 'mol C L-1', 'microzooplankton')
      call self%register_dependency(self%id_mes, 'mes', 'mol C L-1', 'mesozooplankton')
   end subroutine

   subroutine bacteria_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_pisces_bacteria), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: gdept_n, hmld, heup, zoo, mes
      real(rk) :: zdep, zdepbac, ztempbac, zdepmin

      _DOWNWARD_LOOP_BEGIN_
         ! Jorn from p4zrem.F90
         ! Computation of the mean phytoplankton concentration as
         ! a crude estimate of the bacterial biomass
         ! this parameterization has been deduced from a model version
         ! that was modeling explicitely bacteria
         ! -------------------------------------------------------
         _GET_(self%id_gdept_n, gdept_n)
         _GET_SURFACE_(self%id_hmld, hmld)
         _GET_SURFACE_(self%id_heup, heup)
         _GET_(self%id_zoo, zoo)
         _GET_(self%id_mes, mes)
         zdep = MAX( hmld, heup )      ! Jorn: Eq 35a
         IF( gdept_n < zdep ) THEN
            zdepbac = MIN( 0.7_rk * ( zoo + 2._rk* mes ), 4.e-6_rk )    ! Jorn: Eq 35b
            ztempbac   = zdepbac     ! Jorn: this saves the last value that is then used in the deep (clause below)
         ELSE
            zdepmin = MIN( 1._rk, zdep / gdept_n )
            zdepbac = zdepmin**0.683_rk * ztempbac    ! Jorn: Eq 35b
         ENDIF
         _SET_DIAGNOSTIC_(self%id_bact, zdepbac)
      _DOWNWARD_LOOP_END_
   end subroutine

end module