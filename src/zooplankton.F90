#include "fabm_driver.h"

module pisces_zooplankton

   use fabm_types
   use fabm_particle
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_zooplankton
      type (type_state_variable_id) :: id_c
      type (type_state_variable_id) :: id_don, id_dop, id_oxy, id_fer
      type (type_state_variable_id) :: id_dia, id_dfe, id_dsi, id_dch, id_phy, id_nfe, id_nch, id_zoo, id_poc, id_sfe, id_goc, id_bfe, id_gsi
      type (type_state_variable_id) :: id_po4, id_no3, id_nh4, id_doc, id_dic, id_tal, id_poc_waste, id_pof_waste, id_pos_waste, id_cal
      type (type_dependency_id)     :: id_tem, id_nitrfac, id_quotan, id_quotad, id_xfracal, id_wspoc, id_wsgoc
      type (type_diagnostic_variable_id) :: id_zfezoo, id_zgrazing, id_zfrac, id_pcal

      real(rk) :: grazrat, resrat, xkmort, mzrat, xthresh, xkgraz, ferat, epsher, epshermin, unass, sigma, part, grazflux
      real(rk) :: xthreshdia, xthreshphy, xthreshzoo, xthreshpoc
      real(rk) :: xprefn, xprefz, xprefd, xprefc, xsizedia, xdismort, phlim
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_zooplankton), intent(inout), target :: self
      integer,                         intent(in)            :: configunit

      call self%get_parameter(self%epsher, 'epsher', '1', 'maximum growth efficiency', minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%epshermin, 'epshermin', '1', 'minimum growth efficiency', minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%unass, 'unass', '1', 'non-assimilated fraction', default=0.3_rk, minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%sigma, 'sigma', '1', 'excretion as dissolved organic matter', default=0.6_rk, minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%grazrat, 'grazrat', 'd-1', 'maximum grazing rate', minimum=0._rk)
      call self%get_parameter(self%grazflux, 'grazflux', '(m mol C L-1)-1', 'flux-feeding rate', minimum=0._rk)
      call self%get_parameter(self%xkgraz, 'xkgraz', 'mol C L-1', 'half-saturation constant for grazing', default=20.e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xprefn, 'xprefn', '-', 'preference for nanophytoplankton', minimum=0._rk)
      call self%get_parameter(self%xprefd, 'xprefd', '-', 'preference for diatoms', minimum=0._rk)
      call self%get_parameter(self%xprefz, 'xprefz', '-', 'preference for microzooplankton', minimum=0._rk)
      call self%get_parameter(self%xprefc, 'xprefc', '-', 'preference for particulate organic carbon', minimum=0._rk)
      call self%get_parameter(self%xthresh, 'xthresh', 'mol C L-1', 'food threshold', default=0.3e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshdia, 'xthreshdia', 'mol C L-1', 'food threshold for diatoms', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshphy, 'xthreshphy', 'mol C L-1', 'food threshold for nanophytoplankton', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshzoo, 'xthreshzoo', 'mol C L-1', 'food threshold for microzooplankton', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshpoc, 'xthreshpoc', 'mol C L-1', 'food threshold for particulate organic carbon', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%mzrat, 'mzrat', '(umol C L-1)-1 d-1', 'quadratic mortality', minimum=0._rk)
      call self%get_parameter(self%resrat, 'resrat', 'd-1', 'linear mortality', minimum=0._rk)
      call self%get_parameter(self%xkmort, 'xkmort', 'mol C L-1', 'half saturation constant for mortality', default=0.2e-6_rk, minimum=0._rk)
      call self%get_parameter(self%part, 'part', '1', 'fraction of calcite that does not dissolve in guts', minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%ferat, 'ferat', 'mol Fe mol C-1', 'Fe / C ratio', default=10.e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xsizedia, 'xsizedia', 'mol C L-1', 'maximum accessible diatom biomass (above this threshold cells are too large)', default=1e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xdismort, 'xdismort', '1', 'fraction of quadratic mortality directed to nutrient pools', default=( 1._rk - self%epsher - self%unass ) /( 1._rk - self%epsher ), minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%phlim, 'phlim', '1', 'relative grazing on nanophytoplankton if cells are small', minimum=0._rk, maximum=1._rk)

      call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'carbon')
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_c, scale_factor=self%ferat * 1e9_rk)

      call self%register_diagnostic_variable(self%id_zfezoo, 'fezoo', 'nmol Fe m-3 s-1', 'iron recycling rate')
      call self%register_diagnostic_variable(self%id_zgrazing, 'graz', 'mol C m-3 s-1', 'grazing')
      call self%register_diagnostic_variable(self%id_zfrac, 'zfrac', 'mol C m-3 s-1', 'fractionation of large POM')
      call self%register_diagnostic_variable(self%id_pcal, 'pcal', 'mol m-3 s-1', 'calcite production')
      call self%add_to_aggregate_variable(calcite_production, self%id_pcal)

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_state_dependency(self%id_poc_waste, 'poc_waste', 'mol C L-1', 'particulate carbon waste')
      call self%register_state_dependency(self%id_pof_waste, 'pof_waste', 'mol Fe L-1', 'particulate iron waste')
      call self%register_state_dependency(self%id_pos_waste, 'pos_waste', 'mol Si L-1', 'particulate silicon waste')
      call self%register_state_dependency(self%id_cal, 'cal', 'mol C L-1', 'calcite')

      ! Prey
      call self%register_state_dependency(self%id_dia, 'diac', 'mol C L-1', 'diatoms')
      call self%register_state_dependency(self%id_dfe, 'dfe', 'mol Fe L-1', 'diatom iron')
      call self%register_state_dependency(self%id_dsi, 'dsi', 'mol Si L-1', 'diatom silicon')
      call self%register_state_dependency(self%id_dch, 'dch', 'g L-1', 'diatom chlorophyll')
      call self%register_state_dependency(self%id_phy, 'phyc', 'mol C L-1', 'nanophytoplankton')
      call self%register_state_dependency(self%id_nfe, 'nfe', 'mol Fe L-1', 'nanophytoplankton iron')
      call self%register_state_dependency(self%id_nch, 'nch', 'g L-1', 'nanophytoplankton chlorophyll')
      call self%register_state_dependency(self%id_zoo, 'zoo', 'mol C L-1', 'microzooplankton')
      call self%register_state_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_goc, 'goc', 'mol C L-1', 'large particulate organic carbon')
      call self%register_state_dependency(self%id_gsi, 'gsi', 'mol Si L-1', 'large particulate organic silicon')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')
      call self%register_dependency(self%id_quotan, 'quotan', '-', 'proxy for N quota of nanophytoplankton')
      call self%register_dependency(self%id_quotad, 'quotad', '-', 'proxy for N quota of diatoms')
      call self%register_dependency(self%id_xfracal, 'xfracal', '1', 'calcifying fraction of nanophytoplankton')
      call self%register_dependency(self%id_wspoc, 'wspoc', 'm d-1', 'sinking velocity of small particulate organic matter')
      call self%register_dependency(self%id_wsgoc, 'wsgoc', 'm d-1', 'sinking velocity of large particulate organic matter')

      ! Allow wholesale coupling of prey, which then automatically sets up the constituent coupling
      call self%request_coupling_to_model(self%id_dia, 'dia', 'c')
      call self%request_coupling_to_model(self%id_dfe, 'dia', 'fe')
      call self%request_coupling_to_model(self%id_dch, 'dia', 'ch')
      call self%request_coupling_to_model(self%id_dsi, 'dia', 'si')
      call self%request_coupling_to_model(self%id_quotad, 'dia', 'quota')
      call self%request_coupling_to_model(self%id_phy, 'phy', 'c')
      call self%request_coupling_to_model(self%id_nfe, 'phy', 'fe')
      call self%request_coupling_to_model(self%id_nch, 'phy', 'ch')
      call self%request_coupling_to_model(self%id_quotan, 'phy', 'quota')
      call self%request_coupling_to_model(self%id_xfracal, 'phy', 'xfracal')
      call self%request_coupling_to_model(self%id_poc, 'pom', 'c')
      call self%request_coupling_to_model(self%id_sfe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_wspoc, 'pom', 'ws')
      call self%request_coupling_to_model(self%id_goc, 'gom', 'c')
      call self%request_coupling_to_model(self%id_bfe, 'gom', 'fe')
      call self%request_coupling_to_model(self%id_gsi, 'gom', 'si')
      call self%request_coupling_to_model(self%id_cal, 'gom', 'cal')
      call self%request_coupling_to_model(self%id_wsgoc, 'gom', 'ws')

      call self%register_dependency(self%id_tem, standard_variables%temperature)
      call self%register_dependency(self%id_nitrfac, 'nitrfac', '1', 'denitrification factor computed from O2 levels')
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: c, tem, nitrfac
      real(rk) :: tgfunc2, zcompa, zfact, zrespz, ztortz
      real(rk) :: zcompadi, zcompaph, zcompapoc, zcompaz
      real(rk) :: zfood, zfoodlim, zdenom, zdenom2, zgraze
      real(rk) :: dia, dfe, dsi, dch, phy, nfe, nch, zoo, poc, sfe, goc, bfe, gsi, quotan, quotad, xfracal, wspoc, wsgoc
      real(rk) :: zgrazp, zgrazpoc, zgrazd, zgrazz, zgrazpf, zgrazpof, zgrazsf, zgraztotc, zgraztotf, zgraztotn
      real(rk) :: zgrazffeg, zgrazfffg, zgrazffep, zgrazfffp, zproport, zratio, zratio2, zfrac, zfracfe
      real(rk) :: zgrasrat, zgrasratn, zepshert, zbeta, zepsherf, zepsherq, zepsherv, zgrafer, zgrarem, zgrarsig, zmortz, zmortzgoc
      real(rk) :: zprcaca, zfracal, zgrazcal, cal

      _LOOP_BEGIN_
         _GET_(self%id_c, c)
         _GET_(self%id_tem, tem)
         _GET_(self%id_nitrfac, nitrfac)

         tgfunc2 = EXP( 0.07608_rk  * tem )         ! Jorn: from p4zint.F90, equivalent to Eq 25b as NB EXP(0.07608) = 1.079 = b_Z

         zcompa = MAX( ( c - 1.e-9_rk ), 0.e0_rk )
         zfact   = xstep * tgfunc2 * zcompa

         !   Michaelis-Menten mortality rates of microzooplankton - Jorn: last (4th) term in Eq 24
         !   -----------------------------------------------------
         zrespz = self%resrat * zfact * ( c / ( self%xkmort + c )  &
         &        + 3._rk * nitrfac )

         !   Zooplankton mortality. A square function has been selected with
         !   no real reason except that it seems to be more stable and may mimic predation.
         !   Jorn: 3rd term in Eq 24, except that a (1._rk - nitrfac) factor has been added here, and zcompaz is used to introduce a threshold
         !   ------------------------------------------------------------------------------
         ztortz = self%mzrat * 1.e6_rk * zfact * c * (1._rk - nitrfac)

         ! Prey
         _GET_(self%id_dia, dia)
         _GET_(self%id_dfe, dfe)
         _GET_(self%id_dsi, dsi)
         _GET_(self%id_dch, dch)
         _GET_(self%id_phy, phy)
         _GET_(self%id_nfe, nfe)
         _GET_(self%id_nch, nch)
         _GET_(self%id_zoo, zoo)
         _GET_(self%id_poc, poc)
         _GET_(self%id_sfe, sfe)
         _GET_(self%id_goc, goc)
         _GET_(self%id_gsi, gsi)
         _GET_(self%id_bfe, bfe)
         _GET_(self%id_cal, cal)
         _GET_(self%id_wspoc, wspoc)
         _GET_(self%id_wsgoc, wsgoc)
         _GET_(self%id_quotan, quotan)
         _GET_(self%id_quotad, quotad)
         _GET_(self%id_xfracal, xfracal)
         zcompadi  = MIN( MAX( ( dia - self%xthreshdia ), 0.e0_rk ), self%xsizedia )  ! Jorn: xsizedia = 0 for mesozoo
         zcompaph  = MAX( ( phy - self%xthreshphy ), 0.e0_rk )
         zcompapoc = MAX( ( poc - self%xthreshpoc ), 0.e0_rk )
         zcompaz   = MAX( ( zoo - self%xthreshzoo ), 0.e0_rk )

         ! Jorn - mesozo only:
         ! Size effect of nanophytoplankton on grazing : the smaller it is, the less prone
         ! it is to predation by mesozooplankton
         ! -------------------------------------------------------------------------------
         zcompaph  = zcompaph &
            &      * MIN(1._rk, MAX( self%phlim, ( quotan - 0.2_rk) / 0.3_rk ) )

         ! Grazing
         ! ----------------------------------
         zfood     = self%xprefn * zcompaph + self%xprefc * zcompapoc + self%xprefd * zcompadi + self%xprefz * zcompaz   ! Jorn: 1st term in Eq 26a
         zfoodlim  = MAX( 0. , zfood - min(self%xthresh, 0.5_rk * zfood) )                                               ! Jorn: 2nd term in Eq 26a
         zdenom    = zfoodlim / ( self%xkgraz + zfoodlim )
         zdenom2   = zdenom / ( zfood + rtrn )
         zgraze    = self%grazrat * xstep * tgfunc2 * c * (1. - nitrfac)    ! Jorn: compared to paper (Eq 26a), (1._rk - nitrfac) factor seems to have been added

         zgrazp    = zgraze  * self%xprefn * zcompaph  * zdenom2     ! Jorn: ingestion of nanophytoplankton carbon
         zgrazpoc  = zgraze  * self%xprefc * zcompapoc * zdenom2     ! Jorn: ingestion of POC
         zgrazd    = zgraze  * self%xprefd * zcompadi  * zdenom2     ! Jorn: ingestion of diatom carbon
         zgrazz    = zgraze  * self%xprefd * zcompaz   * zdenom2     ! Jorn: ingestion of microzooplankton carbon

         ! Jorn: compute specific loss rates for prey carbon, and apply those to prey iron too
         zgrazpf   = zgrazp    * nfe / (phy + rtrn)   ! Jorn: ingestion of nanophytoplankton Fe
         zgrazpof  = zgrazpoc  * sfe / (poc + rtrn)   ! Jorn: ingestion of POFe
         zgrazsf   = zgrazd    * dfe / (dia + rtrn)   ! Jorn: ingestion of diatom Fe

         ! Flux feeding on GOC
         ! ----------------------------------
         zgrazffeg = self%grazflux  * xstep * wsgoc      &  ! Jorn: Eq 29b but with added factor (1._rk - nitrfac). NB sinking velocity wsgoc should be m d-1!
         &           * tgfunc2 * goc * c &
         &           * (1._rk - nitrfac)
         zgrazfffg = zgrazffeg * bfe / (goc + rtrn)
         zgrazffep = self%grazflux  * xstep *  wspoc     &  ! Jorn: Eq 29a but with added factor (1._rk - nitrfac). NB sinking velocity wsgoc should be m d-1!
         &           * tgfunc2 * poc * c &
         &           * (1._rk - nitrfac)
         zgrazfffp = zgrazffep * sfe / (poc + rtrn)
         !
         zgraztotc = zgrazd + zgrazp + zgrazz + zgrazpoc + zgrazffep + zgrazffeg   ! Jorn: this seems to be potential ingestion, since zgrazffep and zgrazffeg will later be scaled with proportion of filter feeders, zproport
         ! Compute the proportion of filter feeders
         zproport  = (zgrazffep + zgrazffeg)/(rtrn + zgraztotc)
         ! Compute fractionation of aggregates. It is assumed that 
         ! diatoms based aggregates are more prone to fractionation
         ! since they are more porous (marine snow instead of fecal pellets)
         zratio    = gsi / ( goc + rtrn )
         zratio2   = zratio * zratio
         zfrac     = zproport * self%grazflux  * xstep * wsgoc      &
         &          * goc * c          &
         &          * ( 0.2_rk + 3.8_rk * zratio2 / ( 1._rk**2 + zratio2 ) )
         zfracfe   = zfrac * bfe / (goc + rtrn)

         ! Jorn: scale potential ingestion due to flux feeding with the proprotion of predators that is filter-feeding
         zgrazffep = zproport * zgrazffep
         zgrazffeg = zproport * zgrazffeg
         zgrazfffp = zproport * zgrazfffp
         zgrazfffg = zproport * zgrazfffg

         zgraztotc = zgrazp  + zgrazd  + zgrazpoc + zgrazz              + zgrazffep + zgrazffeg    ! Jorn: total ingestion of carbon
         zgraztotf = zgrazpf + zgrazsf + zgrazpof + zgrazz * self%ferat + zgrazfffp + zgrazfffg    ! Jorn: total ingestion of iron, note use of ferat assumes microzoo prey and its predator have same Fe : C
         zgraztotn = zgrazp * quotan + zgrazd * quotad + zgrazpoc + zgrazz + zgrazffep + zgrazffeg ! Jorn: total ingestion of nitrogen, already expressed in carbon units (i.e. relative to rno3)

         ! Grazing by microzooplankton
         _SET_DIAGNOSTIC_(self%id_zgrazing, zgraztotc * 1e3_rk)

         ! Microzooplankton efficiency. 
         ! We adopt a formulation proposed by Mitra et al. (2007)
         ! The gross growth efficiency is controled by the most limiting nutrient.
         ! Growth is also further decreased when the food quality is poor. This is currently
         ! hard coded : it can be decreased by up to 50% (zepsherq)
         ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and 
         ! Fulton, 2012)
         ! -----------------------------------------------------------------------------
         zgrasrat  = ( zgraztotf + rtrn ) / ( zgraztotc + rtrn )  ! Jorn: Fe : C ratio in ingested prey
         !write (*,*) zgraztotf, zgraztotc, rtrn
         zgrasratn = ( zgraztotn + rtrn ) / ( zgraztotc + rtrn )  ! Jorn: N : C ratio in ingested prey, but N already expressed in C units
         zepshert  =  MIN( 1._rk, zgrasratn, zgrasrat / self%ferat)  ! Jorn: Eq 27a, maximum rate of biomass production, derived from incoming C, N, Fe
         zbeta     = MAX(0._rk, (self%epsher - self%epshermin) )
         zepsherf  = self%epshermin + zbeta / ( 1.0_rk + 0.04E6_rk * 12._rk * zfood * zbeta )
         zepsherq  = 0.5_rk + (1.0_rk - 0.5_rk) * zepshert * ( 1.0_rk + 1.0_rk ) / ( zepshert + 1.0_rk )
         zepsherv  = zepsherf * zepshert * zepsherq   ! Jorn: gross growth efficiency 

!         zgrafer   = zgraztotc * MAX( 0._rk , ( 1._rk - self%unass ) * zgrasrat - self%ferat * zepsherv ) 
         zgrafer   = ( 1._rk - self%unass ) * zgraztotf - self%ferat * zepsherv * zgraztotc &  ! Jorn: total dissolved Fe waste (organic + inorganic)
         &         + self%ferat * self%xdismort * ztortz ! Jorn: Eq30b. this line for mesozoo only (xdismort=0 otherwise), ztortz is quadratic mortality
         zgrarem   = zgraztotc * ( 1._rk - zepsherv - self%unass ) &   ! Jorn: total dissolved C/N/P waste (organic + inorganic)
         &         + self%xdismort * ztortz ! Jorn: Eq30b. this line for mesozoo only (xdismort=0 otherwise), ztortz is quadratic mortality

         !  Update of the TRA arrays
         !  ------------------------
         zgrarsig  = zgrarem * self%sigma
         _ADD_SOURCE_(self%id_po4, zgrarsig)
         _ADD_SOURCE_(self%id_nh4, zgrarsig)
         _ADD_SOURCE_(self%id_doc, zgrarem - zgrarsig)
         !
         !IF( ln_ligand ) THEN
         !   tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + (zgrarem - zgrarsig) * ldocz
         !   zzligprod(ji,jj,jk) = (zgrarem - zgrarsig) * ldocz
         !ENDIF
         !
         _ADD_SOURCE_(self%id_oxy, - o2ut * zgrarsig)
         _ADD_SOURCE_(self%id_fer, + zgrafer)
         _SET_DIAGNOSTIC_(self%id_zfezoo, zgrafer * 1e12_rk)
         !prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + zgrapoc
         _ADD_SOURCE_(self%id_dic, zgrarsig)
         _ADD_SOURCE_(self%id_tal, rno3 * zgrarsig)

         !   Update the arrays TRA which contain the biological sources and sinks
         !   --------------------------------------------------------------------
         zmortz = ztortz + zrespz
         zmortzgoc = (1._rk - self%xdismort) * ztortz + zrespz  ! Jorn: Eq 30a: mortality directed to POM pool (as opposed to dissolved pools). TODO! this currently also includes the biomass that goes into the HTLs!
         ! Jorn: separated own biomass, prey, waste, fractionation
         _ADD_SOURCE_(self%id_c, - zmortz + zepsherv * zgraztotc )

         ! Prey losses
         _ADD_SOURCE_(self%id_phy, - zgrazp)
         _ADD_SOURCE_(self%id_dia, - zgrazd)
         _ADD_SOURCE_(self%id_zoo, - zgrazz)
         _ADD_SOURCE_(self%id_nch, - zgrazp  * nch / (phy + rtrn))
         _ADD_SOURCE_(self%id_dch, - zgrazd * dch / (dia + rtrn))
         _ADD_SOURCE_(self%id_dsi, - zgrazd * dsi / (dia + rtrn))
         _ADD_SOURCE_(self%id_nfe, - zgrazpf)
         _ADD_SOURCE_(self%id_dfe, - zgrazsf)
         _ADD_SOURCE_(self%id_poc, - zgrazpoc - zgrazffep)
         _ADD_SOURCE_(self%id_sfe, - zgrazpof - zgrazfffp)
         _ADD_SOURCE_(self%id_goc,            - zgrazffeg)
         _ADD_SOURCE_(self%id_bfe,            - zgrazfffg)
         !prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortz
         !conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazm

         ! Particulate waste from feeding and mortality
         _ADD_SOURCE_(self%id_poc_waste, + zgraztotc * self%unass + zmortzgoc)
         _ADD_SOURCE_(self%id_pos_waste, + zgrazd * dsi / (dia + rtrn))
         _ADD_SOURCE_(self%id_pof_waste, + zgraztotf * self%unass + self%ferat * zmortzgoc)

         ! Fractionation (break-up of large POM into small POM)
         _ADD_SOURCE_(self%id_poc, + zfrac)
         _ADD_SOURCE_(self%id_sfe, + zfracfe)
         _ADD_SOURCE_(self%id_goc, - zfrac)
         _ADD_SOURCE_(self%id_bfe, - zfracfe)
         _SET_DIAGNOSTIC_(self%id_zfrac, + zfrac * 1e3_rk)
         !
         ! Jorn: calcite is consumed as part of flux feeding (uptake proportional to cal/poc * poc_uptake).
         ! The fraction that dissolves in the gut (1-part) is removed from the ambient calcite pool
         zfracal = cal / (poc + goc + rtrn )
         zgrazcal = (zgrazffeg + zgrazpoc) * (1._rk - self%part) * zfracal

         ! calcite production
         zprcaca = xfracal * zgrazp
         !prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         zprcaca = self%part * zprcaca
         _ADD_SOURCE_(self%id_dic, + zgrazcal - zprcaca)
         _ADD_SOURCE_(self%id_tal, 2._rk * (zgrazcal - zprcaca))
         _ADD_SOURCE_(self%id_cal, - zgrazcal + zprcaca)
         _SET_DIAGNOSTIC_(self%id_pcal, zprcaca * 1e3_rk)
      _LOOP_END_
   end subroutine

end module