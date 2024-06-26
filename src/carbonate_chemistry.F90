#include "fabm_driver.h"

module pisces_carbonate_chemistry
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_carbonate_chemistry
      type (type_state_variable_id) :: id_dic
      type (type_state_variable_id) :: id_tal
      type (type_dependency_id) :: id_tempis, id_salinprac, id_rhop, id_sil, id_po4, id_zpres, id_h2co3, id_hi_old
      type (type_surface_dependency_id) :: id_wndm, id_fr_i, id_patm, id_satmco2
      type (type_diagnostic_variable_id) :: id_ph, id_hi, id_CO3, id_CO3sat, id_zomegaca, id_zh2co3
      type (type_surface_diagnostic_variable_id) :: id_Cflx, id_Kg, id_Dpco2, id_pCO2sea
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type

   REAL(rk), PARAMETER :: rgas   = 83.14472_rk      ! universal gas constants

   REAL(rk), PARAMETER ::   calcon = 1.03E-2   ! mean calcite concentration [Ca2+] in sea water [mole/kg solution]

   REAL(rk), PARAMETER :: devk10  = -25.5_rk
   REAL(rk), PARAMETER :: devk11  = -15.82_rk
   REAL(rk), PARAMETER :: devk12  = -29.48_rk
   REAL(rk), PARAMETER :: devk13  = -20.02_rk
   REAL(rk), PARAMETER :: devk14  = -18.03_rk
   REAL(rk), PARAMETER :: devk15  = -9.78_rk
   REAL(rk), PARAMETER :: devk16  = -48.76_rk
   REAL(rk), PARAMETER :: devk17  = -14.51_rk
   REAL(rk), PARAMETER :: devk18  = -23.12_rk
   REAL(rk), PARAMETER :: devk19  = -26.57_rk
   REAL(rk), PARAMETER :: devk110  = -29.48_rk
   !
   REAL(rk), PARAMETER :: devk20  = 0.1271_rk
   REAL(rk), PARAMETER :: devk21  = -0.0219_rk
   REAL(rk), PARAMETER :: devk22  = 0.1622_rk
   REAL(rk), PARAMETER :: devk23  = 0.1119_rk
   REAL(rk), PARAMETER :: devk24  = 0.0466_rk
   REAL(rk), PARAMETER :: devk25  = -0.0090_rk
   REAL(rk), PARAMETER :: devk26  = 0.5304_rk
   REAL(rk), PARAMETER :: devk27  = 0.1211_rk
   REAL(rk), PARAMETER :: devk28  = 0.1758_rk
   REAL(rk), PARAMETER :: devk29  = 0.2020_rk
   REAL(rk), PARAMETER :: devk210  = 0.1622_rk
   !
   REAL(rk), PARAMETER :: devk30  = 0._rk
   REAL(rk), PARAMETER :: devk31  = 0._rk
   REAL(rk), PARAMETER :: devk32  = 2.608E-3_rk
   REAL(rk), PARAMETER :: devk33  = -1.409e-3_rk
   REAL(rk), PARAMETER :: devk34  = 0.316e-3_rk
   REAL(rk), PARAMETER :: devk35  = -0.942e-3_rk
   REAL(rk), PARAMETER :: devk36  = 0._rk
   REAL(rk), PARAMETER :: devk37  = -0.321e-3_rk
   REAL(rk), PARAMETER :: devk38  = -2.647e-3_rk
   REAL(rk), PARAMETER :: devk39  = -3.042e-3_rk
   REAL(rk), PARAMETER :: devk310  = -2.6080e-3_rk
   !
   REAL(rk), PARAMETER :: devk40  = -3.08E-3_rk
   REAL(rk), PARAMETER :: devk41  = 1.13E-3_rk
   REAL(rk), PARAMETER :: devk42  = -2.84E-3_rk
   REAL(rk), PARAMETER :: devk43  = -5.13E-3_rk
   REAL(rk), PARAMETER :: devk44  = -4.53e-3_rk
   REAL(rk), PARAMETER :: devk45  = -3.91e-3_rk
   REAL(rk), PARAMETER :: devk46  = -11.76e-3_rk
   REAL(rk), PARAMETER :: devk47  = -2.67e-3_rk
   REAL(rk), PARAMETER :: devk48  = -5.15e-3_rk
   REAL(rk), PARAMETER :: devk49  = -4.08e-3_rk
   REAL(rk), PARAMETER :: devk410  = -2.84e-3_rk
   !
   REAL(rk), PARAMETER :: devk50  = 0.0877E-3_rk
   REAL(rk), PARAMETER :: devk51  = -0.1475E-3_rk
   REAL(rk), PARAMETER :: devk52  = 0._rk
   REAL(rk), PARAMETER :: devk53  = 0.0794E-3_rk
   REAL(rk), PARAMETER :: devk54  = 0.09e-3_rk
   REAL(rk), PARAMETER :: devk55  = 0.054e-3_rk
   REAL(rk), PARAMETER :: devk56  = 0.3692E-3_rk
   REAL(rk), PARAMETER :: devk57  = 0.0427e-3_rk
   REAL(rk), PARAMETER :: devk58  = 0.09e-3_rk
   REAL(rk), PARAMETER :: devk59  = 0.0714e-3_rk
   REAL(rk), PARAMETER :: devk510  = 0.0_rk

contains

   subroutine initialize(self, configunit)
      class (type_pisces_carbonate_chemistry), intent(inout), target :: self
      integer,                                 intent(in)            :: configunit

      call self%register_implemented_routines((/source_do, source_do_surface/))

      call self%register_state_variable(self%id_dic, 'DIC', 'mol C L-1', 'dissolved inorganic carbon concentration', &
         initial_value=1.99e-3_rk, minimum=0.0_rk, &
         standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_dic, scale_factor=1e6_rk)
      call self%register_state_variable(self%id_tal, 'Alkalini', 'mol L-1', 'alkalinity concentration', &
         initial_value=2.31e-3_rk, minimum=0.0_rk, &
         standard_variable=standard_variables%alkalinity_expressed_as_mole_equivalent)

      call self%register_diagnostic_variable(self%id_PH, 'PH', '1', 'pH', standard_variable=standard_variables%ph_reported_on_total_scale)
      call self%register_diagnostic_variable(self%id_hi, 'hi', 'mol L-1', 'hydrogen ion concentration')
      call self%register_diagnostic_variable(self%id_CO3, 'CO3', 'mol m-3', 'CO3 concentration')
      call self%register_diagnostic_variable(self%id_CO3sat, 'CO3sat', 'mol m-3', 'CO3 concentration at saturation')
      call self%register_diagnostic_variable(self%id_zomegaca, 'zomegaca', '1', 'calcite saturation state', standard_variable=calcite_saturation_state)
      call self%register_diagnostic_variable(self%id_zh2co3, 'zh2co3', 'mol L-1', 'carbonic acid concentration')

      call self%register_diagnostic_variable(self%id_Cflx, 'Cflx', 'mol m-2 s-1', 'air-sea CO2 flux')
      call self%register_diagnostic_variable(self%id_Kg, 'Kg', 'm s-1', 'gas transfer velocity')
      call self%register_diagnostic_variable(self%id_Dpco2, 'Dpco2', 'uatm', 'delta pCO2')
      call self%register_diagnostic_variable(self%id_pCO2sea, 'pCO2sea', 'uatm', 'surface ocean pCO2')

      call self%register_dependency(self%id_tempis, standard_variables%temperature) ! TODO should be in-situ temperature (as opposed to conservative/potential)
      call self%register_dependency(self%id_salinprac, standard_variables%practical_salinity)
      call self%register_dependency(self%id_wndm, standard_variables%wind_speed)
      call self%register_dependency(self%id_fr_i, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_patm, standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_rhop, standard_variables%density)
      call self%register_dependency(self%id_zpres, standard_variables%pressure)
      call self%register_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')
      call self%register_dependency(self%id_satmco2,standard_variables%mole_fraction_of_carbon_dioxide_in_air)
      call self%register_dependency(self%id_h2co3, 'zh2co3', 'mol L-1', 'carbonic acid concentration')
      call self%register_dependency(self%id_hi_old, 'hi', 'mol L-1', 'previous hydrogen ion concentration')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_carbonate_chemistry), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: tempis, salinprac, rhop, dic, tal, po4, sil
      real(rk) :: ztkel, zsal
      real(rk) :: zpres, zsqrt, zsal15, zlogt, ztr, zis, zis2, zisqrt, ztc
      real(rk) :: zcl, zst, zft, zcks, zckf, zckb
      real(rk) :: zck1, zck2, zckw, zck1p, zck2p, zck3p, zcksi, zaksp0
      real(rk) :: total2free, free2SWS, total2SWS, SWS2total
      real(rk) :: zak1, zak2, zakb, zakw, zaksp1, zak1p, zak2p, zak3p, zaksi
      real(rk) :: zcpexp, zcpexp2
      real(rk) :: zbuf1, zbuf2, ak13, ak23, akb3, akw3, aks3, akf3, ak1p3, ak2p3, ak3p3, aksi3
      real(rk) :: aksp, borat, sulfat, fluorid
      real(rk) :: p_hini, zhi, hi
      real(rk) :: zph, zh2co3, zco3, zcalcon, zfact, zomegaca, zco3sat

      _LOOP_BEGIN_
         _GET_(self%id_tempis, tempis)        ! in-situ temperature (TODO! currently this is the "native" temperature from the physical model, which could be conservative or potential)
         _GET_(self%id_salinprac, salinprac)  ! practical salinity (PSU)
         _GET_(self%id_rhop, rhop)            ! density (kg m-3
         _GET_(self%id_zpres, zpres)          ! pressure (dbar)
         _GET_(self%id_dic, dic)              ! total dissolved inorganic carbon (mol C L-1)
         _GET_(self%id_tal, tal)              ! alkalinity (mol L-1)
         _GET_(self%id_sil, sil)              ! silicate (mol Si L-1)
         _GET_(self%id_po4, po4)              ! phosphate (in carbon units! mol C L-1)

         !                             ! SET ABSOLUTE TEMPERATURE
         ztkel = tempis + 273.15_rk
         zsal  = salinprac !(ji,jj,1) + ( 1.- tmask(ji,jj,1) ) * 35.

         zpres = zpres / 10.0 - 1.   ! Jorn: dbar -> bar from p4zchem.F90, added -1 because approximate pressure there equals 0 at 0 depth

         zsqrt  = SQRT( zsal )
         zsal15  = zsqrt * zsal
         zlogt  = LOG( ztkel )
         ztr    = 1. / ztkel
         zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
         zis2   = zis * zis
         zisqrt = SQRT( zis )
         ztc     = tempis !(ji,jj,jk) + ( 1.- tmask(ji,jj,jk) ) * 20.

         ! CHLORINITY (WOOSTER ET AL., 1969)
         zcl     = zsal / 1.80655

         ! TOTAL SULFATE CONCENTR. [MOLES/kg soln]
         zst     = 0.14 * zcl /96.062

         ! TOTAL FLUORIDE CONCENTR. [MOLES/kg soln]
         zft     = 0.000067 * zcl /18.9984

         ! DISSOCIATION CONSTANT FOR SULFATES on free H scale (Dickson 1990)
         zcks    = EXP(-4276.1 * ztr + 141.328 - 23.093 * zlogt         &
         &         + (-13856. * ztr + 324.57 - 47.986 * zlogt) * zisqrt &
         &         + (35474. * ztr - 771.54 + 114.723 * zlogt) * zis    &
         &         - 2698. * ztr * zis**1.5 + 1776.* ztr * zis2         &
         &         + LOG(1.0 - 0.001005 * zsal))

         ! DISSOCIATION CONSTANT FOR FLUORIDES on free H scale (Dickson and Riley 79)
         zckf    = EXP( 1590.2*ztr - 12.641 + 1.525*zisqrt   &
         &         + LOG(1.0d0 - 0.001005d0*zsal)            &
         &         + LOG(1.0d0 + zst/zcks))

         ! DISSOCIATION CONSTANT FOR CARBONATE AND BORATE
         zckb=  (-8966.90 - 2890.53*zsqrt - 77.942*zsal        &
         &      + 1.728*zsal15 - 0.0996*zsal*zsal)*ztr         &
         &      + (148.0248 + 137.1942*zsqrt + 1.62142*zsal)   &
         &      + (-24.4344 - 25.085*zsqrt - 0.2474*zsal)      & 
         &      * zlogt + 0.053105*zsqrt*ztkel

         ! DISSOCIATION COEFFICIENT FOR CARBONATE ACCORDING TO 
         ! MEHRBACH (1973) REFIT BY MILLERO (1995), seawater scale
         zck1    = -1.0*(3633.86*ztr - 61.2172 + 9.6777*zlogt  &
            - 0.011555*zsal + 0.0001152*zsal*zsal)
         zck2    = -1.0*(471.78*ztr + 25.9290 - 3.16967*zlogt      &
            - 0.01781*zsal + 0.0001122*zsal*zsal)

         ! PKW (H2O) (MILLERO, 1995) from composite data
         zckw    = -13847.26 * ztr + 148.9652 - 23.6521 * zlogt + ( 118.67 * ztr    &
                     - 5.977 + 1.0495 * zlogt ) * zsqrt - 0.01615 * zsal

         ! CONSTANTS FOR PHOSPHATE (MILLERO, 1995)
         zck1p    = -4576.752*ztr + 115.540 - 18.453*zlogt   &
         &          + (-106.736*ztr + 0.69171) * zsqrt       &
         &          + (-0.65643*ztr - 0.01844) * zsal

         zck2p    = -8814.715*ztr + 172.1033 - 27.927*zlogt  &
         &          + (-160.340*ztr + 1.3566)*zsqrt          &
         &          + (0.37335*ztr - 0.05778)*zsal

         zck3p    = -3070.75*ztr - 18.126                    &
         &          + (17.27039*ztr + 2.81197) * zsqrt       &
         &          + (-44.99486*ztr - 0.09984) * zsal 

         ! CONSTANT FOR SILICATE, MILLERO (1995)
         zcksi    = -8904.2*ztr  + 117.400 - 19.334*zlogt   &
         &          + (-458.79*ztr + 3.5913) * zisqrt       &
         &          + (188.74*ztr - 1.5998) * zis           &
         &          + (-12.1652*ztr + 0.07871) * zis2       &
         &          + LOG(1.0 - 0.001005*zsal)

         ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
         !       (S=27-43, T=2-25 DEG C) at pres =0 (atmos. pressure) (MUCCI 1983)
         zaksp0  = -171.9065 -0.077993*ztkel + 2839.319*ztr + 71.595*LOG10( ztkel )   &
            &      + (-0.77712 + 0.00284263*ztkel + 178.34*ztr) * zsqrt  &
            &      - 0.07711*zsal + 0.0041249*zsal15

         ! CONVERT FROM DIFFERENT PH SCALES
         total2free  = 1.0/(1.0 + zst/zcks)
         free2SWS    = 1. + zst/zcks + zft/(zckf*total2free)
         total2SWS   = total2free * free2SWS
         SWS2total   = 1.0 / total2SWS

         ! K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
         zak1    = 10**(zck1) * total2SWS
         zak2    = 10**(zck2) * total2SWS
         zakb    = EXP( zckb ) * total2SWS
         zakw    = EXP( zckw )
         zaksp1  = 10**(zaksp0)
         zak1p   = exp( zck1p )
         zak2p   = exp( zck2p )
         zak3p   = exp( zck3p )
         zaksi   = exp( zcksi )
         zckf    = zckf * total2SWS

         ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
         !        (REFERENCE TO CULBERSON & PYTKOQICZ (1968) AS MADE
         !        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
         !        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
         !        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
         !        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
         !        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
         !        & GIESKES (1970), P. 1285-1286 (THE SMALL
         !        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
         !        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON P. 1285))
         zcpexp  = zpres / (rgas*ztkel)
         zcpexp2 = zpres * zcpexp

         ! KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
         !        CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
         !        (CF. BROECKER ET AL., 1982)

         zbuf1  = -     ( devk10 + devk20 * ztc + devk30 * ztc * ztc )
         zbuf2  = 0.5 * ( devk40 + devk50 * ztc )
         ak13 = zak1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk11 + devk21 * ztc + devk31 * ztc * ztc )
         zbuf2  = 0.5 * ( devk41 + devk51 * ztc )
         ak23 = zak2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk12 + devk22 * ztc + devk32 * ztc * ztc )
         zbuf2  = 0.5 * ( devk42 + devk52 * ztc )
         akb3 = zakb * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk13 + devk23 * ztc + devk33 * ztc * ztc )
         zbuf2  = 0.5 * ( devk43 + devk53 * ztc )
         akw3 = zakw * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk14 + devk24 * ztc + devk34 * ztc * ztc )
         zbuf2  = 0.5 * ( devk44 + devk54 * ztc )
         aks3 = zcks * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk15 + devk25 * ztc + devk35 * ztc * ztc )
         zbuf2  = 0.5 * ( devk45 + devk55 * ztc )
         akf3 = zckf * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk17 + devk27 * ztc + devk37 * ztc * ztc )
         zbuf2  = 0.5 * ( devk47 + devk57 * ztc )
         ak1p3 = zak1p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk18 + devk28 * ztc + devk38 * ztc * ztc )
         zbuf2  = 0.5 * ( devk48 + devk58 * ztc )
         ak2p3 = zak2p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk19 + devk29 * ztc + devk39 * ztc * ztc )
         zbuf2  = 0.5 * ( devk49 + devk59 * ztc )
         ak3p3 = zak3p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk110 + devk210 * ztc + devk310 * ztc * ztc )
         zbuf2  = 0.5 * ( devk410 + devk510 * ztc )
         aksi3 = zaksi * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         ! CONVERT FROM DIFFERENT PH SCALES
         total2free  = 1.0/(1.0 + zst/aks3)
         free2SWS    = 1. + zst/aks3 + zft/akf3
         total2SWS   = total2free * free2SWS
         SWS2total   = 1.0 / total2SWS

         ! Convert to total scale
         ak13  = ak13  * SWS2total
         ak23  = ak23  * SWS2total
         akb3  = akb3  * SWS2total
         akw3  = akw3  * SWS2total
         ak1p3 = ak1p3 * SWS2total
         ak2p3 = ak2p3 * SWS2total
         ak3p3 = ak3p3 * SWS2total
         aksi3 = aksi3 * SWS2total
         akf3  = akf3  / total2free

         ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE 
         !        AS FUNCTION OF PRESSURE FOLLOWING MILLERO
         !        (P. 1285) AND BERNER (1976)
         zbuf1  =     - ( devk16 + devk26 * ztc + devk36 * ztc * ztc )
         zbuf2  = 0.5 * ( devk46 + devk56 * ztc )
         aksp = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         ! TOTAL F, S, and BORATE CONCENTR. [MOLES/L]
         borat = 0.0002414 * zcl / 10.811
         sulfat = zst
         fluorid = zft 

         _GET_(self%id_hi_old, p_hini)
         zhi = solve_at_general( rhop, dic, tal, po4, sil, borat, sulfat, fluorid, &
            ak13, ak23, akb3, akw3, aks3, akf3, ak1p3, ak2p3, ak3p3, aksi3, p_hini )
         hi = zhi * rhop / 1000.    ! from p4zlys.F90
         _SET_DIAGNOSTIC_(self%id_hi, hi)
         _SET_DIAGNOSTIC_(self%id_ph, -1. * LOG10( MAX( hi, rtrn ) ))    ! from p4zlys.F90

         ! from p4zlys.F90
         zco3 = dic * ak13 * ak23 / (zhi**2   &
            &             + ak13 * zhi + ak13 * ak23 + rtrn )

         zfact    = rhop / 1000._rk

         ! from p4zflx.F90
         zph   = MAX( hi, 1.e-10_rk ) / zfact
         zh2co3 = dic/(1. + ak13/zph + ak13*ak23/zph**2)
         _SET_DIAGNOSTIC_(self%id_zh2co3, zh2co3)

         ! DEVIATION OF [CO3--] FROM SATURATION VALUE
         ! Salinity dependance in zomegaca and divide by rhop/1000 to have good units
         ! Jorn: from p4zlys.F90, we do this here so we do not need to export zco3 (or its constituents) and aksp
         zcalcon  = calcon * ( salinprac / 35._rk )
         zomegaca = ( zcalcon * zco3 ) / ( aksp * zfact + rtrn )
         zco3sat = aksp * zfact / ( zcalcon + rtrn )
         _SET_DIAGNOSTIC_(self%id_CO3, zco3 * 1.e3_rk)          ! from p4zlys.F90
         _SET_DIAGNOSTIC_(self%id_CO3sat, zco3sat * 1.e3_rk)    ! from p4zlys.F90
         _SET_DIAGNOSTIC_(self%id_zomegaca, zomegaca)
      _LOOP_END_
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_pisces_carbonate_chemistry), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: dic, tal, tempis, salinprac, rhop, wndm, fr_i, patm, satmco2, zh2co3
      real(rk) :: ztkel, zt, zsal, zcek1, chemc(3), ztc, ztc2, ztc3, ztc4, zsch_co2, zws, zkgwan, zkgco2
      real(rk) :: zvapsw, zpco2atm, zxc2, zfugcoeff, zfco2, zfld, zflu, oce_co2

      real(rk), parameter ::   xconv  = 0.01_rk / 3600._rk   !: coefficients for conversion 
      real(rk), parameter ::   atm_per_pa = 1._rk / 101325._rk

      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_dic, dic)
         _GET_(self%id_tal, tal)
         _GET_(self%id_tempis, tempis)          ! temperature (degrees Celsius)
         _GET_(self%id_salinprac, salinprac)    ! practical salinity
         _GET_(self%id_rhop, rhop)              ! density (kg m-3)
         _GET_(self%id_h2co3, zh2co3)           ! carbonic acid concentration (mol L-1)
         _GET_SURFACE_(self%id_wndm, wndm)      ! wind speed (m/s)
         _GET_SURFACE_(self%id_fr_i, fr_i)      ! ice area fraction (1)
         _GET_SURFACE_(self%id_patm, patm)      ! atmospheric pressure (Pa)
         _GET_SURFACE_(self%id_satmco2, satmco2)  ! atmospheric pCO2 (ppm)
         patm = patm * atm_per_pa                 ! convert atmospheric pressure from Pa to atm

         ztkel = tempis + 273.15_rk
         zt    = ztkel * 0.01
         zsal  = salinprac !(ji,jj,1) + ( 1.- tmask(ji,jj,1) ) * 35.
         !                             ! LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1980)
         !                             !     AND FOR THE ATMOSPHERE FOR NON IDEAL GAS
         zcek1 = 9345.17/ztkel - 60.2409 + 23.3585 * LOG(zt) + zsal*(0.023517 - 0.00023656*ztkel    &
         &       + 0.0047036e-4*ztkel**2)

         chemc(1) = EXP( zcek1 ) * 1E-6 * rhop / 1000. ! mol/(L atm)
         chemc(2) = -1636.75 + 12.0408*ztkel - 0.0327957*ztkel**2 + 0.0000316528*ztkel**3
         chemc(3) = 57.7 - 0.118*ztkel

         ztc  = MIN( 35._rk, tempis )
         ztc2 = ztc * ztc
         ztc3 = ztc * ztc2 
         ztc4 = ztc2 * ztc2

         ! Compute the schmidt Number
         zsch_co2 = 2116.8_rk - 136.25_rk * ztc + 4.7353_rk * ztc2 - 0.092307_rk * ztc3 + 0.0007555_rk * ztc4

         !  wind speed 
         zws  = wndm * wndm

         ! Compute the piston velocity for O2 and CO2
         zkgwan = 0.251_rk * zws
         zkgwan = zkgwan * xconv * ( 1._rk - fr_i )   ! Jorn: Eq 82

         zkgco2 = zkgwan * SQRT( 660._rk / zsch_co2 )

         zvapsw    = EXP(24.4543 - 67.4509*(100.0/ztkel) - 4.8489*LOG(ztkel/100) - 0.000544*zsal)   ! Jorn: water vapour pressure (atm)
         zpco2atm = satmco2 * ( patm - zvapsw )
         zxc2      = ( 1.0 - zpco2atm * 1E-6 )**2
         zfugcoeff = EXP( patm * (chemc(2) + 2.0 * zxc2 * chemc(3) )   &
         &           / ( 82.05736 * ztkel ))
         zfco2 = zpco2atm * zfugcoeff

         ! Compute CO2 flux for the sea and air
         zfld = zfco2 * chemc(1) * zkgco2  ! (mol/L) * (m/s)
         zflu = zh2co3 * zkgco2                                   ! (mol/L) (m/s) ?
         oce_co2 = ( zfld - zflu ) !* tmask(ji,jj,1) 
         ! compute the trend
         _ADD_SURFACE_FLUX_(self%id_dic, oce_co2)

         _SET_SURFACE_DIAGNOSTIC_(self%id_Cflx, oce_co2 * 1000._rk)
         _SET_SURFACE_DIAGNOSTIC_(self%id_Kg, zkgco2)
         _SET_SURFACE_DIAGNOSTIC_(self%id_Dpco2, zpco2atm - zh2co3 / (chemc(1) + rtrn))
         _SET_SURFACE_DIAGNOSTIC_(self%id_pCO2sea, zh2co3 / (chemc(1) + rtrn))
      _SURFACE_LOOP_END_
   end subroutine

end module