#include "fabm_driver.h"

! Note: as in the original PISCES code, we provide radiative fluxes (etot_ndcy, etot, pe1, pe2, pe3, emoy, zpar)
! for the *ice-free part* of the water.
! That is, downwelling irradiance just below the surface, which is averaged over the grid cell, is divided by (1 - ice_area_fraction)
! in order to get downwelling irradiance in the ice-free part of the grid cell (NB this assumes no light penetrates the ice!)
! NB the only PISCES variable that is already grid cell averaged (and not for the ice-free part) is etot3, but that variable is not
! used within PISCES, but only by NEMO for the heating of the water column due to light absorption.

module pisces_optics
   use fabm_types
   use fabm_expressions
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_optics
      type (type_dependency_id)                  :: id_chltot, id_e3t_n
      type (type_surface_dependency_id)          :: id_qsr, id_par_varsw, id_hmld, id_fr_i
      type (type_horizontal_dependency_id)       :: id_qsr_mean
      type (type_diagnostic_variable_id)         :: id_pe1, id_pe2, id_pe3, id_etot_ndcy, id_etot
      type (type_surface_diagnostic_variable_id) :: id_heup, id_heup_01, id_emoy, id_zpar

      logical :: ln_varpar
      real(rk) :: xparsw
      real(rk) :: xsi0r
   contains
      procedure :: initialize
      procedure :: do_column
   end type

   REAL(rk), DIMENSION(3,61) :: rkrgb

contains

   subroutine initialize(self, configunit)
      class (type_pisces_optics), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      real(rk) :: parlux, rn_si0

      call self%register_implemented_routines((/source_do_column/))

      call self%get_parameter(self%ln_varpar, 'ln_varpar', '', 'use variable PAR : SWR ratio', default=.true.)
      if (.not. self%ln_varpar) then
         call self%get_parameter(parlux, 'parlux', '-', 'PAR : SWR ratio', default=0.43_rk)
         self%xparsw = parlux / 3.0_rk   ! Divide PAR equally over R-G-B wavebands (division by 3 in Eq 5a)
      end if
      call self%get_parameter(rn_si0, 'rn_si0', 'm', 'extinction depth for non-visible radiation (SWR - PAR)', default=0.35_rk)  ! this comes from NEMO, traqsr
      self%xsi0r = 1._rk / rn_si0

      call trc_oce_rgb(self, rkrgb)

      call self%register_dependency(self%id_chltot, total_chlorophyll)
      call self%register_dependency(self%id_e3t_n, standard_variables%cell_thickness)
      call self%register_dependency(self%id_qsr, standard_variables%surface_downwelling_shortwave_flux)
      if (self%ln_varpar) call self%register_dependency(self%id_par_varsw, 'par_varsw', '1', 'PAR : SWR ratio')
      call self%register_dependency(self%id_qsr_mean, temporal_mean(self%id_qsr, period=86400._rk,resolution=3600._rk))
      call self%register_dependency(self%id_fr_i, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)

      call self%register_diagnostic_variable(self%id_pe1, 'pe1', 'W m-2', 'daily mean PAR in blue band [in ice-free water]', source=source_do_column)
      call self%register_diagnostic_variable(self%id_pe2, 'pe2', 'W m-2', 'daily mean PAR in green band [in ice-free water]', source=source_do_column)
      call self%register_diagnostic_variable(self%id_pe3, 'pe3', 'W m-2', 'daily mean PAR in red band [in ice-free water]', source=source_do_column)
      call self%register_diagnostic_variable(self%id_etot_ndcy, 'etot_ndcy', 'W m-2', 'daily mean PAR [in ice-free water]', source=source_do_column)
      call self%register_diagnostic_variable(self%id_etot, 'etot', 'W m-2', 'instantaneous PAR [in ice-free water]', source=source_do_column)
      call self%register_diagnostic_variable(self%id_heup, 'heup', 'm', 'euphotic layer depth', source=source_do_column)
      call self%register_diagnostic_variable(self%id_heup_01, 'heup_01', 'm', 'depth where daily mean PAR [in ice-free water] equals 0.5 W m-2', source=source_do_column)
      call self%register_diagnostic_variable(self%id_emoy, 'emoy', 'W m-2', 'instantaneous PAR averaged over mixing layer [in ice-free water]', source=source_do_column)
      call self%register_diagnostic_variable(self%id_zpar, 'zpar', 'W m-2', 'daily mean PAR averaged over mixing layer [in ice-free water]', source=source_do_column)
   end subroutine initialize

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_pisces_optics), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: par_varsw, qsr, qsr_mean, zqsr, zqsr_mean, hmld, pqsr100, fr_i
      real(rk) :: ekb, ekg, ekr
      real(rk) :: zchl3d, zchl, e3t_n, gdepw_n
      integer :: irgb
      real(rk) :: f1, f2, f3, pe1, pe2, pe3, heup, heup_01, etot_ndcy, etot
      real(rk) :: zetmp1, zetmp2, zdepmoy
      logical :: first

      _GET_SURFACE_(self%id_qsr, qsr)            ! instantaneous shortwave radiation (W m-2)
      _GET_SURFACE_(self%id_qsr_mean, qsr_mean)  ! daily mean shortwave radiation (W m-2)
      _GET_SURFACE_(self%id_hmld, hmld)          ! mixing layer depth (m)
      _GET_SURFACE_(self%id_fr_i, fr_i)          ! sea ice area fraction (1)

      ! Convert to below-surface irradiances in ice-free part of the grid cell
      qsr      = qsr      / (1. - fr_i + rtrn)
      qsr_mean = qsr_mean / (1. - fr_i + rtrn)

      !  Real shortwave - Jorn: Eq 5a, zqsr is PAR for a single waveband, total par is split equally over R-G-B
      ! Note that SWR is assumed to be just below the water surface:
      ! impacts of reflection and notably ice are assumed to have been already accounted for.
      ! Ice-covered areas are assumed to be completely dark.
      IF( self%ln_varpar ) THEN
         ! PAR is variable fraction of SWR
         _GET_SURFACE_(self%id_par_varsw, par_varsw)
         zqsr      = par_varsw   * qsr      / 3.0_rk   ! Jorn: Eq 5a, distribute equally over three wavebands
         zqsr_mean = par_varsw   * qsr_mean / 3.0_rk   ! Jorn: Eq 5a, distribute equally over three wavebands
      ELSE
         ! PAR is constant fraction of SWR.
         ! Note: xparsw already includes 1/3 scale factor for distribution over three wavebands (Eq 5a)
         zqsr      = self%xparsw * qsr
         zqsr_mean = self%xparsw * qsr_mean
      ENDIF

      !  Light at the euphotic depth 
      pqsr100 = 0.01_rk * 3._rk * zqsr_mean

      ekb = 0._rk
      ekg = 0._rk
      ekr = 0._rk

      heup    = 0._rk
      heup_01 = 0._rk
      gdepw_n = 0._rk
      first = .true.
      zetmp1 = 0._rk
      zetmp2 = 0._rk
      zdepmoy = 0._rk

      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_chltot, zchl3d)
         _GET_(self%id_e3t_n, e3t_n)

         zchl = ( zchl3d + rtrn )                         ! dropped multiplication with 1.e6 as standard variable is already in mg/m3 = ug/L
         zchl = MIN(  10._rk , MAX( 0.05_rk, zchl )  )
         irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )    ! determine index into R-G-B specifc attenuation coefficients based on total chlorophyll

         ! Jorn: ensure index is valid even if zchl is non-finite to prevent crashes due to out-of-bounds rkrgb access
         ! This in principle makes the MIN( 10. , ...) above redundant
         irgb = MIN( SIZE( rkrgb, 2 ) , MAX( 1, irgb )  )

         ! Note ekb, ekg, ekr are depth-integrated attenuation coefficients,
         ! multiplied by 2 since we increment with a whole e3t_n when moving half a grid cell down.
         ! Thus, they need to be multiplied by 0.5 when taking the final EXP.
         ekb = ekb + rkrgb(1,irgb) * e3t_n
         ekg = ekg + rkrgb(2,irgb) * e3t_n
         ekr = ekr + rkrgb(3,irgb) * e3t_n

         ! Note that the irradiances below are all horizontally averaged over the entire grid cell and thus
         ! consider shading by ice. Compared to the original PISCES implementation, this is similar to etot3,
         ! but different from etot_ndcy and etot, which in PISCES are representative for the ice-free part
         ! of the grid cell only, and thus higher. We consider this interpretation specific to phytoplankton,
         ! and a similar ice correction is therefore instead done there.
         f1 = EXP( -0.5_rk * ekb )
         f2 = EXP( -0.5_rk * ekg )
         f3 = EXP( -0.5_rk * ekr )
         pe1 = zqsr_mean * f1
         pe2 = zqsr_mean * f2
         pe3 = zqsr_mean * f3
         etot = zqsr * (f1 + f2 + f3)
         etot_ndcy = pe1 + pe2 + pe3
         _SET_DIAGNOSTIC_(self%id_pe1, pe1)
         _SET_DIAGNOSTIC_(self%id_pe2, pe2)
         _SET_DIAGNOSTIC_(self%id_pe3, pe3)
         _SET_DIAGNOSTIC_(self%id_etot_ndcy, etot_ndcy)
         _SET_DIAGNOSTIC_(self%id_etot, etot)

         gdepw_n = gdepw_n + e3t_n
         IF (first .or. etot_ndcy >= pqsr100) heup    = gdepw_n  ! Euphotic layer depth
         IF (first .or. etot_ndcy >= 0.50)    heup_01 = gdepw_n  ! Euphotic layer depth (light level definition)

         IF (gdepw_n <= hmld) THEN
            zetmp1 = zetmp1 + etot      * e3t_n   ! integrating instantaneous PAR
            zetmp2 = zetmp2 + etot_ndcy * e3t_n   ! integrating daily mean PAR
            zdepmoy = gdepw_n
         END IF

         ekb = ekb + rkrgb(1,irgb) * e3t_n
         ekg = ekg + rkrgb(2,irgb) * e3t_n
         ekr = ekr + rkrgb(3,irgb) * e3t_n
         first = .false.

         !f1 = EXP( -0.5_rk * ekb )
         !f2 = EXP( -0.5_rk * ekg )
         !f3 = EXP( -0.5_rk * ekr )
         !etot_bot = zqsr * (f1 * f1 + f2 + f3)
      _DOWNWARD_LOOP_END_

      heup    = MIN( 300._rk, heup    )
      heup_01 = MIN( 300._rk, heup_01 )

      _SET_SURFACE_DIAGNOSTIC_(self%id_heup, heup)
      _SET_SURFACE_DIAGNOSTIC_(self%id_heup_01, heup_01)
      _SET_SURFACE_DIAGNOSTIC_(self%id_emoy, zetmp1 / (zdepmoy + rtrn))
      _SET_SURFACE_DIAGNOSTIC_(self%id_zpar, zetmp2 / (zdepmoy + rtrn))
   end subroutine

   ! From NEMO, trc_oce.F90
   SUBROUTINE trc_oce_rgb( self, prgb )
      class (type_pisces_optics), intent(in) :: self
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   Set a look up table for the optical coefficients
      !!                i.e. the attenuation coefficient for R-G-B light 
      !!                tabulated in Chlorophyll class (from JM Andre)
      !!
      !! ** Action  :   prgb(3,61) tabulated R-G-B attenuation coef. 
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      REAL(rk), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !
      INTEGER  ::   jc     ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(rk) ::   zchl   ! temporary scalar
      REAL(rk), DIMENSION(4,61) ::   zrgb   ! tabulated attenuation coefficient (formerly read in 'kRGB61.txt')
      !!----------------------------------------------------------------------
      !
      !  Chlorophyll        !     Blue attenuation     !     Green attenuation    !     Red attenuation      !
      zrgb(1, 1) =  0.010   ;   zrgb(2, 1) = 0.01618   ;   zrgb(3, 1) = 0.07464   ;   zrgb(4, 1) = 0.37807
      zrgb(1, 2) =  0.011   ;   zrgb(2, 2) = 0.01654   ;   zrgb(3, 2) = 0.07480   ;   zrgb(4, 2) = 0.37823
      zrgb(1, 3) =  0.013   ;   zrgb(2, 3) = 0.01693   ;   zrgb(3, 3) = 0.07499   ;   zrgb(4, 3) = 0.37840
      zrgb(1, 4) =  0.014   ;   zrgb(2, 4) = 0.01736   ;   zrgb(3, 4) = 0.07518   ;   zrgb(4, 4) = 0.37859
      zrgb(1, 5) =  0.016   ;   zrgb(2, 5) = 0.01782   ;   zrgb(3, 5) = 0.07539   ;   zrgb(4, 5) = 0.37879
      zrgb(1, 6) =  0.018   ;   zrgb(2, 6) = 0.01831   ;   zrgb(3, 6) = 0.07562   ;   zrgb(4, 6) = 0.37900
      zrgb(1, 7) =  0.020   ;   zrgb(2, 7) = 0.01885   ;   zrgb(3, 7) = 0.07586   ;   zrgb(4, 7) = 0.37923
      zrgb(1, 8) =  0.022   ;   zrgb(2, 8) = 0.01943   ;   zrgb(3, 8) = 0.07613   ;   zrgb(4, 8) = 0.37948
      zrgb(1, 9) =  0.025   ;   zrgb(2, 9) = 0.02005   ;   zrgb(3, 9) = 0.07641   ;   zrgb(4, 9) = 0.37976
      zrgb(1,10) =  0.028   ;   zrgb(2,10) = 0.02073   ;   zrgb(3,10) = 0.07672   ;   zrgb(4,10) = 0.38005
      zrgb(1,11) =  0.032   ;   zrgb(2,11) = 0.02146   ;   zrgb(3,11) = 0.07705   ;   zrgb(4,11) = 0.38036
      zrgb(1,12) =  0.035   ;   zrgb(2,12) = 0.02224   ;   zrgb(3,12) = 0.07741   ;   zrgb(4,12) = 0.38070
      zrgb(1,13) =  0.040   ;   zrgb(2,13) = 0.02310   ;   zrgb(3,13) = 0.07780   ;   zrgb(4,13) = 0.38107
      zrgb(1,14) =  0.045   ;   zrgb(2,14) = 0.02402   ;   zrgb(3,14) = 0.07821   ;   zrgb(4,14) = 0.38146
      zrgb(1,15) =  0.050   ;   zrgb(2,15) = 0.02501   ;   zrgb(3,15) = 0.07866   ;   zrgb(4,15) = 0.38189
      zrgb(1,16) =  0.056   ;   zrgb(2,16) = 0.02608   ;   zrgb(3,16) = 0.07914   ;   zrgb(4,16) = 0.38235
      zrgb(1,17) =  0.063   ;   zrgb(2,17) = 0.02724   ;   zrgb(3,17) = 0.07967   ;   zrgb(4,17) = 0.38285
      zrgb(1,18) =  0.071   ;   zrgb(2,18) = 0.02849   ;   zrgb(3,18) = 0.08023   ;   zrgb(4,18) = 0.38338
      zrgb(1,19) =  0.079   ;   zrgb(2,19) = 0.02984   ;   zrgb(3,19) = 0.08083   ;   zrgb(4,19) = 0.38396
      zrgb(1,20) =  0.089   ;   zrgb(2,20) = 0.03131   ;   zrgb(3,20) = 0.08149   ;   zrgb(4,20) = 0.38458
      zrgb(1,21) =  0.100   ;   zrgb(2,21) = 0.03288   ;   zrgb(3,21) = 0.08219   ;   zrgb(4,21) = 0.38526
      zrgb(1,22) =  0.112   ;   zrgb(2,22) = 0.03459   ;   zrgb(3,22) = 0.08295   ;   zrgb(4,22) = 0.38598
      zrgb(1,23) =  0.126   ;   zrgb(2,23) = 0.03643   ;   zrgb(3,23) = 0.08377   ;   zrgb(4,23) = 0.38676
      zrgb(1,24) =  0.141   ;   zrgb(2,24) = 0.03842   ;   zrgb(3,24) = 0.08466   ;   zrgb(4,24) = 0.38761
      zrgb(1,25) =  0.158   ;   zrgb(2,25) = 0.04057   ;   zrgb(3,25) = 0.08561   ;   zrgb(4,25) = 0.38852
      zrgb(1,26) =  0.178   ;   zrgb(2,26) = 0.04289   ;   zrgb(3,26) = 0.08664   ;   zrgb(4,26) = 0.38950
      zrgb(1,27) =  0.200   ;   zrgb(2,27) = 0.04540   ;   zrgb(3,27) = 0.08775   ;   zrgb(4,27) = 0.39056
      zrgb(1,28) =  0.224   ;   zrgb(2,28) = 0.04811   ;   zrgb(3,28) = 0.08894   ;   zrgb(4,28) = 0.39171
      zrgb(1,29) =  0.251   ;   zrgb(2,29) = 0.05103   ;   zrgb(3,29) = 0.09023   ;   zrgb(4,29) = 0.39294
      zrgb(1,30) =  0.282   ;   zrgb(2,30) = 0.05420   ;   zrgb(3,30) = 0.09162   ;   zrgb(4,30) = 0.39428
      zrgb(1,31) =  0.316   ;   zrgb(2,31) = 0.05761   ;   zrgb(3,31) = 0.09312   ;   zrgb(4,31) = 0.39572
      zrgb(1,32) =  0.355   ;   zrgb(2,32) = 0.06130   ;   zrgb(3,32) = 0.09474   ;   zrgb(4,32) = 0.39727
      zrgb(1,33) =  0.398   ;   zrgb(2,33) = 0.06529   ;   zrgb(3,33) = 0.09649   ;   zrgb(4,33) = 0.39894
      zrgb(1,34) =  0.447   ;   zrgb(2,34) = 0.06959   ;   zrgb(3,34) = 0.09837   ;   zrgb(4,34) = 0.40075
      zrgb(1,35) =  0.501   ;   zrgb(2,35) = 0.07424   ;   zrgb(3,35) = 0.10040   ;   zrgb(4,35) = 0.40270
      zrgb(1,36) =  0.562   ;   zrgb(2,36) = 0.07927   ;   zrgb(3,36) = 0.10259   ;   zrgb(4,36) = 0.40480
      zrgb(1,37) =  0.631   ;   zrgb(2,37) = 0.08470   ;   zrgb(3,37) = 0.10495   ;   zrgb(4,37) = 0.40707
      zrgb(1,38) =  0.708   ;   zrgb(2,38) = 0.09056   ;   zrgb(3,38) = 0.10749   ;   zrgb(4,38) = 0.40952
      zrgb(1,39) =  0.794   ;   zrgb(2,39) = 0.09690   ;   zrgb(3,39) = 0.11024   ;   zrgb(4,39) = 0.41216
      zrgb(1,40) =  0.891   ;   zrgb(2,40) = 0.10374   ;   zrgb(3,40) = 0.11320   ;   zrgb(4,40) = 0.41502
      zrgb(1,41) =  1.000   ;   zrgb(2,41) = 0.11114   ;   zrgb(3,41) = 0.11639   ;   zrgb(4,41) = 0.41809
      zrgb(1,42) =  1.122   ;   zrgb(2,42) = 0.11912   ;   zrgb(3,42) = 0.11984   ;   zrgb(4,42) = 0.42142
      zrgb(1,43) =  1.259   ;   zrgb(2,43) = 0.12775   ;   zrgb(3,43) = 0.12356   ;   zrgb(4,43) = 0.42500
      zrgb(1,44) =  1.413   ;   zrgb(2,44) = 0.13707   ;   zrgb(3,44) = 0.12757   ;   zrgb(4,44) = 0.42887
      zrgb(1,45) =  1.585   ;   zrgb(2,45) = 0.14715   ;   zrgb(3,45) = 0.13189   ;   zrgb(4,45) = 0.43304
      zrgb(1,46) =  1.778   ;   zrgb(2,46) = 0.15803   ;   zrgb(3,46) = 0.13655   ;   zrgb(4,46) = 0.43754
      zrgb(1,47) =  1.995   ;   zrgb(2,47) = 0.16978   ;   zrgb(3,47) = 0.14158   ;   zrgb(4,47) = 0.44240
      zrgb(1,48) =  2.239   ;   zrgb(2,48) = 0.18248   ;   zrgb(3,48) = 0.14701   ;   zrgb(4,48) = 0.44765
      zrgb(1,49) =  2.512   ;   zrgb(2,49) = 0.19620   ;   zrgb(3,49) = 0.15286   ;   zrgb(4,49) = 0.45331
      zrgb(1,50) =  2.818   ;   zrgb(2,50) = 0.21102   ;   zrgb(3,50) = 0.15918   ;   zrgb(4,50) = 0.45942
      zrgb(1,51) =  3.162   ;   zrgb(2,51) = 0.22703   ;   zrgb(3,51) = 0.16599   ;   zrgb(4,51) = 0.46601
      zrgb(1,52) =  3.548   ;   zrgb(2,52) = 0.24433   ;   zrgb(3,52) = 0.17334   ;   zrgb(4,52) = 0.47313
      zrgb(1,53) =  3.981   ;   zrgb(2,53) = 0.26301   ;   zrgb(3,53) = 0.18126   ;   zrgb(4,53) = 0.48080
      zrgb(1,54) =  4.467   ;   zrgb(2,54) = 0.28320   ;   zrgb(3,54) = 0.18981   ;   zrgb(4,54) = 0.48909
      zrgb(1,55) =  5.012   ;   zrgb(2,55) = 0.30502   ;   zrgb(3,55) = 0.19903   ;   zrgb(4,55) = 0.49803
      zrgb(1,56) =  5.623   ;   zrgb(2,56) = 0.32858   ;   zrgb(3,56) = 0.20898   ;   zrgb(4,56) = 0.50768
      zrgb(1,57) =  6.310   ;   zrgb(2,57) = 0.35404   ;   zrgb(3,57) = 0.21971   ;   zrgb(4,57) = 0.51810
      zrgb(1,58) =  7.079   ;   zrgb(2,58) = 0.38154   ;   zrgb(3,58) = 0.23129   ;   zrgb(4,58) = 0.52934
      zrgb(1,59) =  7.943   ;   zrgb(2,59) = 0.41125   ;   zrgb(3,59) = 0.24378   ;   zrgb(4,59) = 0.54147
      zrgb(1,60) =  8.912   ;   zrgb(2,60) = 0.44336   ;   zrgb(3,60) = 0.25725   ;   zrgb(4,60) = 0.55457
      zrgb(1,61) = 10.000   ;   zrgb(2,61) = 0.47804   ;   zrgb(3,61) = 0.27178   ;   zrgb(4,61) = 0.56870
      !
      prgb(:,:) = zrgb(2:4,:)
      !
      !r_si2 = 1.e0 / zrgb(2, 1)        ! blue with the smallest chlorophyll concentration)
      !IF(lwp) WRITE(numout,*) '      RGB longest depth of extinction    r_si2 = ', r_si2
      !
      DO jc = 1, 61                         ! check
         zchl = zrgb(1,jc)
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )
         IF( irgb /= jc ) THEN
            WRITE(*,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL self%fatal_error( 'trc_oce_rgb' ,'inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      !
   END SUBROUTINE trc_oce_rgb

end module