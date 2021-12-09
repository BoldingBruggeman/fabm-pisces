#include "fabm_driver.h"

module pisces_phytoplankton

   use fabm_types
   use fabm_particle
   use fabm_expressions
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_phytoplankton
      type (type_state_variable_id)     :: id_c, id_ch, id_fe, id_si
      type (type_state_variable_id)     :: id_no3, id_nh4, id_po4, id_sil, id_biron, id_doc, id_dic, id_tal, id_oxy
      type (type_state_variable_id)     :: id_poc, id_sfe, id_goc, id_gsi, id_bfe, id_cal, id_prodpoc, id_prodgoc
      type (type_dependency_id)         :: id_tem, id_gdept_n, id_xdiss
      type (type_dependency_id)         :: id_pe1, id_pe2, id_pe3, id_etot_ndcy, id_etot_w
      type (type_surface_dependency_id) :: id_zstrn, id_hmld, id_heup_01, id_etot_wm
      type (type_surface_dependency_id) :: id_gphit, id_fr_i, id_xksi_
      type (type_horizontal_dependency_id) :: id_xksi
      type (type_diagnostic_variable_id) :: id_quota, id_xfracal, id_pcal
      type (type_diagnostic_variable_id) :: id_zlim1, id_zlim2, id_zlim3, id_zlim4, id_xlim
      type (type_diagnostic_variable_id) :: id_PPPHY, id_PPNEW, id_PBSi, id_PFe
      type (type_diagnostic_variable_id) :: id_zprmax, id_Mu, id_Llight

      logical :: diatom
      logical :: calcify
      real(rk) :: mumax0
      real(rk) :: logbp
      real(rk) :: fpday
      real(rk) :: bresp
      real(rk) :: pislope_s
      real(rk) :: pislope_l
      real(rk) :: xadap
      real(rk) :: excret
      !real(rk) :: concpo4
      real(rk) :: concnh4
      real(rk) :: concno3
      real(rk) :: xksi1
      real(rk) :: xksi2
      real(rk) :: concfer
      real(rk) :: xsizer
      real(rk) :: grosip
      real(rk) :: qfelim
      real(rk) :: fecm
      real(rk) :: chlcm
      real(rk) :: chlcmin
      real(rk) :: xsize
      real(rk) :: texcret
      real(rk) :: caco3r
      real(rk) :: xmort, wchl, wchlm, mprat, xkmort
   contains
      procedure :: initialize
      procedure :: do
   end type

   type, extends(type_base_model) :: type_par
      type (type_surface_dependency_id) :: id_hmld, id_heup_01
      type (type_dependency_id) :: id_pe1, id_pe2, id_pe3, id_e3t_n
      type (type_diagnostic_variable_id) :: id_etot
      type (type_surface_diagnostic_variable_id) :: id_etotm
      real(rk) :: beta1, beta2, beta3
   contains
      procedure :: initialize => par_initialize
      procedure :: do_column  => par_do_column
   end type

   type, extends(type_base_model) :: type_silicate_half_saturation
      type (type_dependency_id)                  :: id_sil
      type (type_surface_diagnostic_variable_id) :: id_xksi
      real(rk) :: concsil
      real(rk) :: xksilim
   contains
      procedure :: initialize => silicate_half_saturation_initialize
      procedure :: do_surface => silicate_half_saturation_do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_phytoplankton), intent(inout), target :: self
      integer,                           intent(in)            :: configunit

      class (type_par),                      pointer :: par_model
      class (type_silicate_half_saturation), pointer :: silicate_half_saturation

      allocate(par_model)

      call self%get_parameter(self%diatom, 'diatom', '', 'use silicate', default=.false.)
      call self%get_parameter(self%calcify, 'calcify', '', 'calcify', default=.false.)
      call self%get_parameter(self%mumax0, 'mumax0', 'd-1', 'maximum growth rate at 0 degrees Celsius', default=0.8_rk)    ! default=0.8 in NEMO-PISCES 4 code, confirmed by Olivier Aumont 2021-04-21
      call self%get_parameter(self%logbp, 'logbp', '-', '(log of) Temperature sensitivity of growth', default=0.063913_rk) ! AC -07.12.2021 - Aumont et al. 2015, Table 1, Eq. 4a, Using log(b_p) instead of b_p. 
      call self%get_parameter(self%fpday, 'fpday', '-', 'day-length factor for growth', default=1.5_rk) ! AC -07.12.2021 - Aumont et al. 2015, Eq 3a
      call self%get_parameter(self%bresp, 'bresp', 'd-1', 'basal respiration', default=0.033_rk)
      call self%get_parameter(self%pislope_s, 'pislope_s', 'g C (g Chl)-1 (W m-2)-1 d-1', 'P-I slope for small cells', default=2._rk)
      call self%get_parameter(self%pislope_l, 'pislope_l', 'g C (g Chl)-1 (W m-2)-1 d-1', 'P-I slope for large cells', default=self%pislope_s) ! AC -07.12.2021 - changed default value to pislope_s to ensure common perturbation.
      call self%get_parameter(self%xadap, 'xadap', '-', 'acclimation factor to low light', default=0._rk)
      call self%get_parameter(self%excret, 'excret', '1', 'fraction of production that is excreted', default=0.05_rk)
      call self%get_parameter(par_model%beta1, 'beta1', '1', 'absorption in blue part of the light')
      call self%get_parameter(par_model%beta2, 'beta2', '1', 'absorption in green part of the light')
      call self%get_parameter(par_model%beta3, 'beta3', '1', 'absorption in red part of the light')
      !call self%get_parameter(self%concpo4, 'concpo4', 'mol C L-1', 'minimum half-saturation constant for phosphate (phosphorus units multiplied with C:P ratio of biomass)')  ! ??? code seems to use concnh4
      call self%get_parameter(self%concnh4, 'concnh4', 'mol C L-1', 'minimum half-saturation constant for ammonium (nitrogen units multiplied with C:N ratio of biomass)')
      call self%get_parameter(self%concno3, 'concno3', 'mol C L-1', 'minimum half-saturation constant for nitrate (nitrogen units multiplied with C:N ratio of biomass)')
      if (self%diatom) then
         allocate(silicate_half_saturation)
         call self%add_child(silicate_half_saturation, 'silicate_half_saturation')
         call self%get_parameter(silicate_half_saturation%concsil, 'concsil', 'mol Si L-1', 'minimum half-saturation constant for silicate uptake', default=1e-6_rk)
         call self%get_parameter(silicate_half_saturation%xksilim, 'xksilim', 'mol Si L-1', 'parameter for the half-saturation constant for silicate uptake', default=16.5e-6_rk)  ! default=16.6 in paper
         call self%get_parameter(self%xksi1, 'xksi1', 'mol Si L-1', 'parameter 1 for Si / C', default=2.e-6_rk)
         call self%get_parameter(self%xksi2, 'xksi2', 'mol Si L-1', 'parameter 2 for Si / C', default=20.e-6_rk)
      end if
      call self%get_parameter(self%concfer, 'concfer', 'mol Fe L-1', 'minimum half-saturation constant for iron uptake')
      call self%get_parameter(self%xsizer, 'xsizer', '-', 'size ratio', default=3._rk)
      if (self%diatom) call self%get_parameter(self%grosip, 'grosip', 'mol Si/mol C', 'optimal Si / C uptake ratio', default=0.159_rk)
      call self%get_parameter(self%qfelim, 'qfelim', 'mol Fe (mol C)-1', 'optimal iron quota', default=7.e-6_rk)
      call self%get_parameter(self%fecm, 'fecm', 'mol Fe (mol C)-1', 'maximum iron quota', default=40.e-6_rk)
      call self%get_parameter(self%chlcm, 'chlcm', 'g Chl (g C)-1', 'maximum Chl / C ratio')
      call self%get_parameter(self%chlcmin, 'chlcmin', 'g Chl (g C)-1', 'minimum Chl / C ratio', default=0.0033_rk)  ! 0.004 in namelist_pisces_ref
      call self%get_parameter(self%xsize, 'xsize', 'mol C L-1', 'threshold concentration for size dependency (biomass above this threshold consists of large cells)', default=1.e-6_rk)
      call self%get_parameter(self%xmort, 'xmort', 'mol C L-1', 'threshold concentration for mortality')
      call self%get_parameter(self%wchl, 'wchl', 'd-1 (umol C L-1)-1', 'quadratic mortality', default=0.01_rk)
      call self%get_parameter(self%wchlm, 'wchlm', 'd-1 (umol C L-1)-1', 'maximum additional quadratic mortality due to nutrient limitation', default=0._rk)
      call self%get_parameter(self%mprat, 'mprat', 'd-1', 'mortality', default=0.01_rk)
      call self%get_parameter(self%xkmort, 'xkmort', 'mol C L-1', 'half saturation constant for mortality', default=2.e-7_rk)
      !!!call self%get_parameter(self%lthet, 'lthet', '-', 'proportional loss of ligands due to Fe uptake', default=1.0_rk)
      self%texcret = 1._rk - self%excret
      if (self%calcify) call self%get_parameter(self%caco3r, 'caco3r', '1', 'mean rain ratio', default=0.3_rk)

      ! Set up submodel for computing available light
      call self%add_child(par_model, 'par')
      call par_model%request_coupling(par_model%id_heup_01, '../heup_01')
      call par_model%request_coupling(par_model%id_pe1, '../pe1')
      call par_model%request_coupling(par_model%id_pe2, '../pe2')
      call par_model%request_coupling(par_model%id_pe3, '../pe3')

      call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'carbon')
      call self%register_state_variable(self%id_ch, 'ch', 'g Chl L-1', 'chlorophyll')
      call self%register_state_variable(self%id_fe, 'fe', 'mol Fe L-1', 'iron')
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)
      call self%add_to_aggregate_variable(total_chlorophyll, self%id_ch, scale_factor=1e6_rk)

      if (self%diatom) then
         call self%register_state_variable(self%id_si, 'si', 'mol Si L-1', 'silicon')
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_si, scale_factor=1e6_rk)
      end if

      call self%register_state_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_prodpoc, 'prodpoc', 'mol C L-1', 'produced small particulate organic carbon')

      call self%register_state_dependency(self%id_goc, 'goc', 'mol C L-1', 'large particulate organic carbon')
      call self%register_state_dependency(self%id_gsi, 'gsi', 'mol Si L-1', 'large particulate organic silicate')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')
      call self%register_state_dependency(self%id_prodgoc, 'prodgoc', 'mol C L-1', 'produced large particulate organic carbon')
      call self%register_state_dependency(self%id_cal, 'cal', 'mol C L-1', 'calcite')

      call self%request_coupling_to_model(self%id_poc, 'pom', 'c')
      call self%request_coupling_to_model(self%id_sfe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_prodpoc, 'pom', 'prod')
      call self%request_coupling_to_model(self%id_goc, 'gom', 'c')
      call self%request_coupling_to_model(self%id_bfe, 'gom', 'fe')
      call self%request_coupling_to_model(self%id_gsi, 'gom', 'si')
      call self%request_coupling_to_model(self%id_prodgoc, 'gom', 'prod')
      if (self%calcify) then
         call self%request_coupling_to_model(self%id_cal, 'gom', 'cal')
      else
         call self%request_coupling(self%id_cal, 'zero')
      end if

      call self%register_dependency(self%id_tem, standard_variables%temperature)
      call self%register_dependency(self%id_zstrn, 'zstrn', 'h', 'day length')
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      call self%register_dependency(self%id_xdiss, shear_rate)
      call self%register_dependency(self%id_heup_01, 'heup_01', 'm', 'euphotic layer depth (PAR > 0.5 W m-2)')
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_pe1, 'pe1', 'W m-2', 'daily mean PAR in blue waveband')
      call self%register_dependency(self%id_pe2, 'pe2', 'W m-2', 'daily mean PAR in green waveband')
      call self%register_dependency(self%id_pe3, 'pe3', 'W m-2', 'daily mean PAR in red waveband')
      call self%register_dependency(self%id_etot_ndcy, 'etot_ndcy', 'W m-2', 'daily mean PAR')
      call self%register_dependency(self%id_etot_w, 'etot_w', 'W m-2', 'daily mean PAR weighted by wavelength-specific absorption coefficients')
      call self%register_dependency(self%id_etot_wm, 'etot_wm', 'W m-2', 'daily mean PAR weighted by wavelength-specific absorption coefficients, averaged over euphotic layer')
      call self%register_dependency(self%id_fr_i, standard_variables%ice_area_fraction)
      call self%request_coupling(self%id_etot_w, 'par/etot')
      call self%request_coupling(self%id_etot_wm, 'par/etotm')

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_biron, 'biron', 'mol Fe L-1', 'bioavailable iron')
      call self%register_state_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')

      if (self%diatom) then
         call self%register_dependency(self%id_gphit, standard_variables%latitude)
         call self%register_dependency(self%id_xksi_, 'xksi', 'mol Si L-1', 'instantaneous silicate half-saturation constant')
         call silicate_half_saturation%request_coupling(silicate_half_saturation%id_sil, '../sil')
         call self%request_coupling(self%id_xksi_, 'silicate_half_saturation/xksi')
         call self%register_dependency(self%id_xksi, temporal_maximum(self%id_xksi_, period=nyear_len * rday, resolution=nyear_len * rday, missing_value=2.e-6_rk))
      else
         call self%request_coupling(self%id_sil, 'zero')
      end if

      call self%register_diagnostic_variable(self%id_quota, 'quota', 'mol N (mol C)-1', 'proxy of the N/C ratio')
      call self%register_diagnostic_variable(self%id_zlim1, 'LN', '-', 'nitrogen limitation term')
      call self%register_diagnostic_variable(self%id_zlim2, 'LP', '-', 'phosphorus limitation term')
      call self%register_diagnostic_variable(self%id_zlim3, 'LSi', '-', 'silicate limitation term')
      call self%register_diagnostic_variable(self%id_zlim4, 'LFe', '-', 'iron limitation term')
      call self%register_diagnostic_variable(self%id_xlim, 'Lnut', '-', 'nutrient limitation term')
      call self%register_diagnostic_variable(self%id_PPPHY, 'PPPHY', 'mol C m-3 s-1', 'primary production')
      call self%register_diagnostic_variable(self%id_PPNEW, 'PPNEW', 'mol C m-3 s-1', 'new primary production')
      if (self%diatom) call self%register_diagnostic_variable(self%id_PBSi, 'PBSi', 'mol Si m-3 s-1', 'biogenic silica production')
      call self%register_diagnostic_variable(self%id_PFe, 'PFe', 'mol Fe m-3 s-1', 'biogenic iron production')
      call self%register_diagnostic_variable(self%id_zprmax, 'mumax', 's-1', 'maximum growth rate after temperature correction')
      call self%register_diagnostic_variable(self%id_Mu, 'Mu', 's-1', 'realized growth rate')
      call self%register_diagnostic_variable(self%id_Llight, 'Llight', '1', 'light limitation term')
      if (self%calcify) call self%register_diagnostic_variable(self%id_xfracal, 'xfracal', '1', 'calcifying fraction')
      call self%register_diagnostic_variable(self%id_pcal, 'pcal', 'mol m-3 s-1', 'calcite production')
      call self%add_to_aggregate_variable(calcite_production, self%id_pcal)

   end subroutine initialize

   subroutine par_initialize(self, configunit)
      class (type_par), intent(inout), target :: self
      integer,          intent(in)            :: configunit

      call self%register_diagnostic_variable(self%id_etot, 'etot', 'W m-2', 'daily mean PAR (sum over all wavebands, weighted by absorption coefficients)', source=source_do_column)
      call self%register_diagnostic_variable(self%id_etotm, 'etotm', 'W m-2', 'daily mean PAR averaged over mixed layer (sum over all wavebands, weighted by absorption coefficients)', source=source_do_column)
      call self%register_dependency(self%id_pe1, 'pe1', 'W m-2', 'daily mean PAR in blue waveband')
      call self%register_dependency(self%id_pe2, 'pe2', 'W m-2', 'daily mean PAR in green waveband')
      call self%register_dependency(self%id_pe3, 'pe3', 'W m-2', 'daily mean PAR in red waveband')
      call self%register_dependency(self%id_e3t_n, standard_variables%cell_thickness)
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      call self%register_dependency(self%id_heup_01, 'heup_01', 'm', 'euphotic layer depth (PAR > 0.5 W m-2)')
   end subroutine

   subroutine par_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_par), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: hmld, heup_01, pe1, pe2, pe3, gdepw_n, e3t_n, zetmp1, zdepmoy, etot, etotm

      gdepw_n = 0._rk
      zetmp1 = 0._rk
      zdepmoy = 0._rk
      _GET_SURFACE_(self%id_hmld, hmld)
      _GET_SURFACE_(self%id_heup_01, heup_01)
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_pe1, pe1)     ! daily mean PAR in blue waveband
         _GET_(self%id_pe2, pe2)     ! daily mean PAR in green waveband
         _GET_(self%id_pe3, pe3)     ! daily mean PAR in red waveband
         _GET_(self%id_e3t_n, e3t_n) ! cell thickness
         etot = self%beta1 * pe1 + self%beta2 * pe2 + self%beta2 * pe2   ! Eq 5b, daily mean PAR weighted by waveband-specific absorption
         _SET_DIAGNOSTIC_(self%id_etot, etot)
         gdepw_n = gdepw_n + e3t_n
         IF (gdepw_n <= MIN(hmld, heup_01)) THEN
            zetmp1 = zetmp1 + etot * e3t_n
            zdepmoy = gdepw_n
         END IF
      _DOWNWARD_LOOP_END_

      etotm = zetmp1 / ( zdepmoy + rtrn )
      _SET_SURFACE_DIAGNOSTIC_(self%id_etotm, etotm)
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_phytoplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      ! Coefficient for iron limitation - Jorn Eq 20
      real(rk), parameter ::  xcoef1   = 0.0016  / 55.85  
      real(rk), parameter ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
      real(rk), parameter ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5 

      real(rk) :: c, ch, fe, si
      real(rk) :: nh4, no3, po4, biron, sil
      real(rk) :: tem, gdept_n, zstrn, hmld, heup_01, etot_ndcy, etot_w, etot_wm, gphit, fr_i
      real(rk) :: tgfunc, zconc, zconc2, z1_trb, concfe, zconcno3, zconcnh4, zdenom, xno3, xnh4, xpo4, zlim1, zlim2, xlim, xlimsi
      real(rk) :: xksi, zlim3
      real(rk) :: zratio, zironmin, zlim4, xlimfe
      real(rk) :: zprmax, zval, zmxl_chl, zmxl_fac, zpr
      real(rk) :: ztn, zadap, pislope, zpislopead, zpislope, zprch
      real(rk) :: zlim, zsilim, zsilfac, zsiborn, zsilfac2, zysopt
      real(rk) :: zprorca, zpronew, zproreg, zdocprod, zmax, zprofe
      real(rk) :: ztot, zprod, zprochl, chlcm_n, quota
      real(rk) :: ztem1, ztem2, zetot1, zetot2, xfracal
      real(rk) :: xfraresp, xfratort, zcompa, zsizerat, zrespp, ztortp, zmortp, zfactfe, zfactch, zfactsi, zprcaca, xdiss

      _LOOP_BEGIN_
         _GET_(self%id_c, c)                      ! carbon (mol/L)
         _GET_(self%id_ch, ch)                    ! chlorophyll (g/L)
         _GET_(self%id_fe,  fe)                   ! internal iron pool (mol/L)

         _GET_(self%id_nh4, nh4)                  ! ammonium (mol C/L - nitrogen units multiplied with C:N ratio of biomass)
         _GET_(self%id_no3, no3)                  ! nitrate (mol C/L - nitrogen units multiplied with C:N ratio of biomass)
         _GET_(self%id_po4, po4)                  ! phosphate (mol C/L - phosphorus units multiplied with C:P ratio of biomass)
         _GET_(self%id_biron, biron)              ! bioavailable iron (mol Fe/L)
         _GET_(self%id_sil, sil)                  ! ambient silicate concentration (mol Si/L)

         _GET_(self%id_tem, tem)                  ! temperature (degrees Celsius)
         _GET_(self%id_gdept_n, gdept_n)          ! depth (m)
         _GET_(self%id_xdiss, xdiss)              ! shear rate (s-1)
         _GET_SURFACE_(self%id_zstrn, zstrn)      ! day length (h)
         _GET_SURFACE_(self%id_hmld, hmld)        ! mixing layer depth (m)
         _GET_SURFACE_(self%id_heup_01, heup_01)  ! euphotic layer depth (m) based on PAR > 0.5 W m-2
         _GET_(self%id_etot_ndcy, etot_ndcy)      ! daily mean Photosynthetically Available Radiation (W m-2), equal weights for all wavelengths
         _GET_(self%id_etot_w, etot_w)            ! daily mean Photosynthetically Available Radiation (W m-2) weighted by absorption per waveband
         _GET_SURFACE_(self%id_etot_wm, etot_wm)  ! daily mean Photosynthetically Available Radiation (W m-2) weighted by absorption per waveband, averaged over mixed/euphotic layer
         _GET_SURFACE_(self%id_fr_i, fr_i)        ! sea ice area fraction (1)

         IF (gdept_n > hmld) etot_wm = etot_w     ! below turbocline: the experienced daily mean PAR is the actual in-situ value (not the average over the mixing layer)

         ! ======================================================================================
         ! Jorn: From p4zint
         !tgfunc = EXP( 0.063913 * tem )  ! Jorn: Eq 4a in PISCES-v2 paper, NB EXP(0.063913) = 1.066 = b_P
         tgfunc = EXP( self%logbp * tem )  ! AC 07.12.2021 - replaced hard-coded value with parameter for bp. Log(bp) considered instead of bp for perturbation purposes).
         ! ======================================================================================

         ! ======================================================================================
         ! Jorn: From p4zlim

         zconc   = MAX( 0._rk , c - self%xsize ) ! Jorn: Eq 7b, biomass (mol C/L) in large cells
         zconc2  = c - zconc                     ! Jorn: equivalent to Eq 7a, carbon biomass (mol C/L) in small cells
         z1_trb   = 1._rk / ( c + rtrn )         ! 1 / carbon biomass

         concfe   = MAX( self%concfer, ( zconc2 * self%concfer + self%concfer * self%xsizer * zconc ) * z1_trb ) ! Jorn: Eq 7c (18b), size-weighted iron half saturation
         zconcno3 = MAX( self%concno3, ( zconc2 * self%concno3 + self%concno3 * self%xsizer * zconc ) * z1_trb ) ! Jorn: Eq 7c, size-weighted nitrate half saturation
         zconcnh4 = MAX( self%concnh4, ( zconc2 * self%concnh4 + self%concnh4 * self%xsizer * zconc ) * z1_trb ) ! Jorn: Eq 7c, size-weighted ammonium half saturation

         ! Michaelis-Menten Limitation term for nutrients - Jorn Eqs 6d, 6e
         ! -----------------------------------------------
         zdenom = 1._rk /  ( zconcno3 * zconcnh4 + zconcnh4 * no3 + zconcno3 * nh4 )
         xno3 = no3 * zconcnh4 * zdenom
         xnh4 = nh4 * zconcno3 * zdenom
         !
         zlim1    = xno3 + xnh4                ! Jorn: Eq 6c, nitrogen limitation (dimensionless)
         zlim2    = po4 / ( po4 + zconcnh4 )   ! Jorn: Eq 6b: phosphorus limitation (dimensionless), note it uses NH4 half saturation (same units - concentration of carbon equivalents!)
         if (self%diatom) then
            ! Jorn: From p4zint
            _GET_SURFACE_(self%id_xksi, xksi)
            zlim3    = sil / ( sil + xksi )    ! Eq 11b
         else
            zlim3 = 1._rk
         end if
         zratio   = fe* z1_trb  ! Jorn: iron quota (based on internal iron pool, not ambient concentration!), mol Fe/mol C
         zironmin = xcoef1 * ch * z1_trb + xcoef2 * zlim1 + xcoef3 * xno3  ! Eq 20, mol Fe/mol C
         zlim4    = MAX( 0.,( zratio - zironmin ) / self%qfelim )   ! Jorn Eq 6f, iron limitation (dimensionless)
         xpo4 = zlim2
         xlimfe = MIN( 1., zlim4 )
         xlim = MIN( zlim1, zlim2, zlim3, zlim4 )   ! Jorn: Eq 11a, combined N, P, (Si), Fe limitation factor (dimensionless)
         xlimsi = MIN( zlim1, zlim2, zlim4 )        ! Jorn: Eq 6a (also part of Eq 23a): combined limitation by all nutrients (N, P, Fe) EXCEPT Si (dimensionless)
         ! ======================================================================================

         _SET_DIAGNOSTIC_(self%id_zlim1, zlim1)
         _SET_DIAGNOSTIC_(self%id_zlim2, zlim2)
         _SET_DIAGNOSTIC_(self%id_zlim3, zlim3)
         _SET_DIAGNOSTIC_(self%id_zlim4, zlim4)
         _SET_DIAGNOSTIC_(self%id_xlim, xlim)

         ! Computation of the optimal production - Jorn note conversion to 1/s
         zprmax = self%mumax0 * r1_rday * tgfunc   ! Jorn: Eq 4b in PISCES-v2 paper; NEMO-PISCES uses a hardcoded mumax0 = 0.8, which evolved from the value of 0.6 in the paper (Olivier Aumont 2021-04-21)

         ! Impact of the day duration and light intermittency on phytoplankton growth
         IF( etot_ndcy > 1.E-3 ) THEN
            zval = MAX( 1., zstrn )   ! Jorn: clip day length to a minimum of 1 hour
            IF( gdept_n <= hmld ) THEN
               zval = zval * MIN(1._rk, heup_01 / ( hmld + rtrn ))   ! Jorn: when in mixing layer, multiply with fraction of time spent in euphotic depth; it seems to be an easier-to-understand substitute for Eq 3b-3d
            ENDIF
            zmxl_chl = zval / 24._rk  ! Jorn: from number of hours to fraction of day
            !zmxl_fac = 1.5_rk * zval / ( 12._rk + zval )  ! Jorn: Eq 3a in PISCES-v2 paper - but note that time spent in euphotic layer has already been incorporated in zval! Eqs 3b-3d are not used
            zmxl_fac = self%fpday * zval / ( 12._rk + zval )  ! AC 07.12.2021 - fpday as parameter 

            zpr = zprmax * zmxl_fac  ! Jorn: product of muP and f1*f2 (sort of - those have been changed) in Eq 2a, units are 1/s

            ! Maximum light intensity
            !WHERE( zstrn(:,:) < 1.e0 ) zstrn(:,:) = 24. ! Jorn: seems unused???

            ! Computation of the P-I slope for nanos and diatoms (Jorn: product of alpha and chl:C)
            ztn         = MAX( 0._rk, tem - 15._rk )
            zadap       = self%xadap * ztn / ( 2._rk + ztn )
            pislope = (self%pislope_s * zconc2 + self%pislope_l * zconc) * z1_trb   ! Weighted mean of PI slope (g C (g Chl)-1 d-1 (W m-2)-1) for small and large cells
            zpislopead = pislope * ( 1._rk + zadap  * EXP( -0.25_rk * etot_w ) )  & ! PI slope multiplied with Chl/C, resulting units are d-1 (W m-2)-1. NB xadap = zadap = 0
            &                   *ch / ( c * 12._rk + rtrn)                          ! g chl:g C

            ! Computation of production function for Carbon - Jorn: Eq 2a in PISCES-v2 paper
            !  ---------------------------------------------
            zpislope = zpislopead / ( ( r1_rday + self%bresp * r1_rday ) &  ! Jorn: prefactor mu_ref=1 for first r1_rday is omitted. Resulting units are (W m-2)-1
            &            * zmxl_fac * rday + rtrn)
            zpr = zpr * ( 1._rk- EXP( -zpislope * etot_w )  )               ! Units are 1/s

            !  Computation of production function for Chlorophyll - Jorn Eq15b
            !--------------------------------------------------
            zpislope = zpislopead / ( zprmax * zmxl_chl * rday + rtrn )   ! note zprmax is in 1/s and multiplied with rday to convert to d-1. Resulting units are (W m-2)-1. No divison by nutrient limitation (similar to Eq 2a replacing 2b)
            zprch = zprmax * ( 1._rk - EXP( -zpislope * etot_wm ) )       ! Units 1/s, note this uses mean PAR experienced in the euphotic layer (or 0 if below ML) - units are 1/d

            !  Computation of a proxy of the N/C ratio - Jorn: this is dimensionless and thus relative to rno3 (it is not the absolute N:C!)
            !  ---------------------------------------
            zval = MIN( xpo4, ( xnh4 + xno3 ) )   &    ! Jorn: nitrogen and phosphate limitation, divided by light limitation
            &      * zprmax / ( zpr + rtrn )
            quota = MIN( 1., 0.2_rk + 0.8_rk * zval )

            IF( self%diatom ) THEN
               _GET_SURFACE_(self%id_gphit, gphit)
                  !    Si/C of diatoms - Jorn: section 4.1.5 in PISCES-v2 paper
                  !    ------------------------
                  !    Si/C increases with iron stress and silicate availability
                  !    Si/C is arbitrariliy increased for very high Si concentrations
                  !    to mimic the very high ratios observed in the Southern Ocean (silpot2)
               zlim  = sil / ( sil + self%xksi1 )                                                                             ! Jorn Eq 23c
               zsilim = MIN( zpr / ( zprmax + rtrn ), xlimsi )                                                                ! Jorn Eq 23a
               zsilfac = 4.4_rk * EXP( -4.23_rk * zsilim ) * MAX( 0._rk, MIN( 1._rk, 2.2_rk * ( zlim - 0.5_rk ) )  ) + 1._rk  ! Jorn Eq 22, 23b
               zsiborn = sil * sil * sil
               IF (gphit < -30._rk ) THEN   ! threshold is 0 degrees in paper
                  zsilfac2 = 1._rk + 2._rk * zsiborn / ( zsiborn + self%xksi2**3 )  ! Eq 22 (last part), 23d
               ELSE
                  zsilfac2 = 1._rk +         zsiborn / ( zsiborn + self%xksi2**3 )  ! ??? 23d suggests this term is 1 for lat > 0
               ENDIF
               zysopt = self%grosip * zlim * zsilfac * zsilfac2  ! Eq 22 but it seems to miss the MIN(5.4, ...) clipping used there, realized Si / C uptake ratio
            ELSE
               zysopt = 0._rk
            ENDIF

            !  Mixed-layer effect on production 
            !  Sea-ice effect on production
            zpr = zpr * ( 1._rk - fr_i )   ! Jorn: production was computed using PAR in ice-free water; here that is converted into production over the entire cell (ice-free + ice-covered fractions)

            ! Computation of the various production terms 
            zprorca = zpr  * xlim* c                           ! total production (mol C/L/s), dropped multiplication with rfact2 [time step in seconds]
            zpronew  = zprorca* xno3 / ( xno3 + xnh4 + rtrn )  ! Eq 8, new production (mol C/L/s)
            !
            zratio = fe / ( c * self%fecm + rtrn )   ! Jorn: internal iron pool relative to maximum value (dimensionless)
            zmax   = MAX( 0., ( 1._rk - zratio ) / ABS( 1.05_rk - zratio ) )            ! ratio in Eq 17 (dimensionless)
            zprofe = self%fecm * zprmax * ( 1.0 - fr_i )  &                             ! Increase in internal iron pool in mol Fe/L/s
            &             * ( 4._rk - 4.5_rk * xlimfe / ( xlimfe + 0.5_rk ) )    &      ! Eq 19, note xlimfe is based on internal iron quota (not ambient concentration)
            &             * biron / ( biron + concfe )  &                               ! Eq 18a
            &             * zmax * c                                                    ! Jorn: dropped multiplication with rfact2 [time step in seconds]

            ! Computation of the chlorophyll production terms
            !  production terms for nanophyto. ( chlorophyll )
            ztot = etot_wm / ( zmxl_chl + rtrn )         ! Jorn: PAR/L_day in Eq 15a, the divison of 24h mean PAR by day length makes it equivalent to mean PAR experienced *during the daytime*
            zprod    = rday * zprorca * zprch * xlim     ! Eq15b Note zprorca was the increment in carbon (mol C/L) over a single time step, but we divide by timestep and thus have the rate of production
            zprochl = self%chlcmin * 12._rk * zprorca       ! Jorn: first part of Eq14, increase in Chl associated with increase in carbon (using carbon production and minimum Chl:C), units are g Chl/L/s
            chlcm_n   = MIN ( self%chlcm, ( self%chlcm / (1. - 1.14 / 43.4 *tem)) * (1. - 1.14 / 43.4 * 20.))  ! temperature correction of max Chl?
            zprochl = zprochl + (chlcm_n-self%chlcmin) * 12. * zprod / &
                                    & (  zpislopead * ztot +rtrn)
         ELSE
            zprochl = 0._rk
            zprorca = 0._rk
            zpronew = 0._rk
            zprofe = 0._rk
            zysopt = 0._rk
            zpr = 0._rk
            quota = 1._rk
         ENDIF

         _ADD_SOURCE_(self%id_ch, zprochl * self%texcret)

         !   Update the arrays TRA which contain the biological sources and sinks
         zproreg  = zprorca - zpronew    ! Regenerated production, equivalent to Eq 8 (part 2)
         zdocprod = self%excret * zprorca
         _ADD_SOURCE_(self%id_po4, - zprorca)
         _ADD_SOURCE_(self%id_no3, - zpronew)
         _ADD_SOURCE_(self%id_nh4, - zproreg)
         _ADD_SOURCE_(self%id_c, zprorca * self%texcret)
         _ADD_SOURCE_(self%id_fe, zprofe * self%texcret)
         _ADD_SOURCE_(self%id_doc, zdocprod)
         _ADD_SOURCE_(self%id_oxy, o2ut * zproreg + ( o2ut + o2nit ) * zpronew)
         !
         _ADD_SOURCE_(self%id_biron, -self%texcret * zprofe)
         if (self%diatom) then
            _ADD_SOURCE_(self%id_si, zprorca * zysopt * self%texcret)
            _ADD_SOURCE_(self%id_sil, -self%texcret * zprorca * zysopt)
         end if
         _ADD_SOURCE_(self%id_dic, -zprorca)
         _ADD_SOURCE_(self%id_tal, rno3 * zpronew - rno3 * zproreg)

         !
         !IF( ln_ligand ) THEN
         !   IF( etot_ndcy > 1.E-3 ) THEN
         !      zdocprod = excret * zprorca
         !      zfeup    = texcret * zprofe
         !      _ADD_SOURCE_(self%id_lgw, zdocprod * ldocp - zfeup * plig * self%lthet)
         !      zpligprod1 = zdocprod * ldocp
         !      zpligprod2 = zfeup * plig * self%lthet
         !   ENDIF
         !ENDIF

         _SET_DIAGNOSTIC_(self%id_quota, quota)

       !! Total primary production per year
       !IF( iom_use( "tintpp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
       !     & tpp = glob_sum( 'p4zprod', ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )
       !
       !IF( lk_iomput .AND.  knt == nrdttrc ) THEN
       !   zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
       !   !
         _SET_DIAGNOSTIC_(self%id_PPPHY, zprorca * 1.e+3) ! primary production
         _SET_DIAGNOSTIC_(self%id_PPNEW, zpronew * 1.e+3) ! new primary production
         if (self%diatom) _SET_DIAGNOSTIC_(self%id_PBSi, zprorca * 1.e+3 * zysopt) ! biogenic silica production
         _SET_DIAGNOSTIC_(self%id_PFe, zprofe * 1.e+3) ! biogenic iron production
       !   IF( ln_ligand ) THEN
       !     CALL iom_put( "LPRODP"  , zpligprod1(:,:,:) * 1e9 * zfact * tmask(:,:,:) )
       !     CALL iom_put( "LDETP"   , zpligprod2(:,:,:) * 1e9 * zfact * tmask(:,:,:) )
       !   ENDIF
         _SET_DIAGNOSTIC_(self%id_zprmax, zprmax) ! Maximum growth rate
         _SET_DIAGNOSTIC_(self%id_Mu, zpr * xlim) ! Realized growth rate
         _SET_DIAGNOSTIC_(self%id_Llight, zpr / (zprmax + rtrn))  ! light limitation term
       !   CALL iom_put( "TPP"     , ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * zfact * tmask(:,:,:)  )  ! total primary production
       !   CALL iom_put( "TPNEW"   , ( zpronewn(:,:,:) + zpronewd(:,:,:) ) * zfact * tmask(:,:,:)  ) ! total new production
       !   CALL iom_put( "TPBFE"   , ( zprofen(:,:,:) + zprofed(:,:,:) ) * zfact * tmask(:,:,:)  )  ! total biogenic iron production
       !   CALL iom_put( "tintpp"  , tpp * zfact )  !  global total integrated primary production molC/s
       !ENDIF   

         if (self%calcify) then
            ! Jorn: from p4zlim.F90, Eq 77
            zlim1 =  ( no3 * self%concnh4 + nh4 * self%concno3 )    &
               &   / ( self%concno3 * self%concnh4 + self%concnh4 * no3 + self%concno3 * nh4 ) 
            zlim2  = po4 / ( po4 + self%concnh4 )
            zlim3  = fe / ( fe +  5.E-11_rk   )   ! iron half saturation hardcoded at 0.05 nmol/L (much lower than 1-3 nmol/L used for nanophytoplankotn and diatoms!)
            ztem1  = MAX( 0._rk, tem )
            ztem2  = tem - 10._rk
            zetot1 = MAX( 0._rk, etot_ndcy - 1._rk) / ( 4._rk + etot_ndcy ) 
            zetot2 = 30._rk / ( 30._rk + etot_ndcy )

            xfracal = self%caco3r * MIN( zlim1, zlim2, zlim3 )                  &
               &                  * ztem1 / ( 0.1_rk + ztem1 )                     &
               &                  * MAX( 1._rk, c * 1.e6_rk / 2._rk )  &
               &                  * zetot1 * zetot2               &
               &                  * ( 1._rk + EXP(-ztem2 * ztem2 / 25._rk ) )         &
               &                  * MIN( 1., 50._rk / ( hmld + rtrn ) )
            xfracal = MIN( 0.8_rk , xfracal )
            xfracal = MAX( 0.02_rk, xfracal )
            _SET_DIAGNOSTIC_(self%id_xfracal, xfracal)
            xfraresp = 0.5_rk * xfracal   ! Jorn: fraction of material produced by linear mortality that is sent to large detritus pool
            xfratort = 0.5_rk * xfracal   ! Jorn: fraction of material produced by quadratic mortality that is sent to large detritus pool
         else
            xfracal = 0._rk
            xfraresp = 1.0_rk   ! Jorn: fraction of material produced by linear mortality that is sent to large detritus pool
            xfratort = 0.5_rk   ! Jorn: fraction of material produced by quadratic mortality that is sent to large detritus pool
         end if

         ! Jorn: from p4zmort.F90
         zcompa = MAX( ( c - self%xmort ), 0.e0_rk )

         !    Aggregation term for diatoms is increased in case of nutrient
         !    stress as observed in reality. The stressed cells become more
         !    sticky and coagulate to sink quickly out of the euphotic zone
         !     ------------------------------------------------------------
         ! Jorn: for non-diatoms this will not have any effect because wchldm will be 0 (see zrespp)
         zlim2   = xlim * xlim
         zlim1   = 0.25_rk * ( 1._rk - zlim2 ) / ( 0.25_rk + zlim2 )    ! Jorn: seems to have replaced Eq 13

         !     When highly limited by macronutrients, very small cells 
         !     dominate the community. As a consequence, aggregation
         !     due to turbulence is negligible. Mortality is also set
         !     to 0
         if (self%calcify) then
            zsizerat = MIN(1._rk, MAX( 0._rk, (quota - 0.2_rk) / 0.3_rk) ) * c
         else
            zsizerat = c
         end if

         !     Squared mortality of Phyto similar to a sedimentation term during
         !     blooms (Doney et al. 1996)
         zrespp = (self%wchl + self%wchlm * zlim1) * 1.e6_rk * xstep * xdiss * zcompa * zsizerat    ! Jorn: 3rd term in Eqs 1,9, except for wchlm contribution, zsizerat, minimum threshold in zcompa

         !     Phytoplankton mortality. This mortality loss is slightly
         !     increased when nutrients are limiting phytoplankton growth
         !     as observed for instance in case of iron limitation.
         ztortp = self%mprat * xstep * zcompa / ( self%xkmort + c ) * zsizerat    ! Jorn: hyperbolic part of 5th term in Eq 37, except for zsizerat, minimum threshold in zcompa

         zmortp = zrespp + ztortp

         !   Update the arrays TRA which contains the biological sources and sinks

         zfactfe = fe/(c+rtrn)
         zfactch = ch/(c+rtrn)
         _ADD_SOURCE_(self%id_c, - zmortp)
         _ADD_SOURCE_(self%id_ch, - zmortp * zfactch)
         _ADD_SOURCE_(self%id_fe, - zmortp * zfactfe)
         if (self%diatom) then
            _GET_(self%id_si, si)
            zfactsi = si/(c+rtrn)
            _ADD_SOURCE_(self%id_si, - zmortp * zfactsi)
            _ADD_SOURCE_(self%id_gsi,  zmortp * zfactsi)
         end if

         ! Jorn: if not self%calcify, xfracal will be 0
         zprcaca = xfracal * zmortp
         _ADD_SOURCE_(self%id_dic, - zprcaca)
         _ADD_SOURCE_(self%id_tal, - 2._rk * zprcaca)
         _ADD_SOURCE_(self%id_cal, + zprcaca)
         _SET_DIAGNOSTIC_(self%id_pcal, zprcaca * 1e3_rk)

         _ADD_SOURCE_(self%id_poc, + ( 1._rk - xfraresp ) * zrespp + ( 1._rk - xfratort ) * ztortp)
         _ADD_SOURCE_(self%id_goc, + xfraresp * zrespp + xfratort * ztortp)
         _ADD_SOURCE_(self%id_prodpoc, + ( 1._rk - xfraresp ) * zrespp + ( 1._rk - xfratort ) * ztortp)
         _ADD_SOURCE_(self%id_prodgoc, + xfraresp * zrespp + xfratort * ztortp)
         _ADD_SOURCE_(self%id_sfe, + (( 1._rk - xfraresp ) * zrespp + ( 1._rk - xfratort ) * ztortp) * zfactfe)
         _ADD_SOURCE_(self%id_bfe, + (xfraresp * zrespp + xfratort * ztortp) * zfactfe)

      _LOOP_END_
   end subroutine

   subroutine silicate_half_saturation_initialize(self, configunit)
      class (type_silicate_half_saturation), intent(inout), target :: self
      integer,                               intent(in)            :: configunit

      call self%register_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')
      call self%register_diagnostic_variable(self%id_xksi, 'xksi', 'mol Si L-1', 'silicate half saturation')
   end subroutine

   subroutine silicate_half_saturation_do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_silicate_half_saturation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: sil, zvar, xksi

      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_sil, sil)
         zvar = sil * sil
         xksi = self%concsil + 7.e-6_rk * zvar / ( self%xksilim * self%xksilim + zvar )    ! Eq 12, note self%concsil=1e-6 is hardcoded in NEMO-PISCES, p4zint.F90
         _SET_SURFACE_DIAGNOSTIC_(self%id_xksi, xksi)
      _SURFACE_LOOP_END_
   end subroutine

end module