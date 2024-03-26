#include "fabm_driver.h"

module pisces_shared

   use fabm_types

   implicit none

   public

   real(rk), parameter :: rno3    =  16._rk / 122._rk
   real(rk), parameter :: po4r    =   1._rk / 122._rk
   real(rk), parameter :: o2nit   =  32._rk / 122._rk
   real(rk), parameter :: o2ut    = 133._rk / 122._rk
   real(rk), parameter :: rdenit  =  ( ( o2ut + o2nit ) * 0.80_rk - rno3 - rno3 * 0.60_rk ) / rno3

   real(rk), parameter :: rtrn    = 0.5_rk * EPSILON( 1.e0_rk )    !: truncation value
   real(rk), parameter :: rday    = 24._rk*60._rk*60._rk      !: day                                [s]
   real(rk), parameter :: r1_rday = 1._rk / rday

   real(rk), parameter :: xstep = r1_rday

   real(rk), parameter :: nyear_len = 365._rk  ! Jorn - TODO should be made dynamic to account for leap years, different calendars, etc.

   real(rk) :: maxdt = 1800._rk ! Jorn: maximum time step in seconds, used to limit benthic-pelagic fluxes in order to prevent negative values

   type (type_universal_standard_variable), parameter :: total_chlorophyll = type_universal_standard_variable(name='total_chlorophyll', units='mg m-3', aggregate_variable=.true.)
   type (type_interior_standard_variable), parameter :: shear_rate = type_interior_standard_variable(name='shear_rate', units='s-1')
   type (type_surface_standard_variable), parameter :: mixed_layer_thickness_defined_by_vertical_tracer_diffusivity = type_surface_standard_variable(name='mixed_layer_thickness_defined_by_vertical_tracer_diffusivity', units='m')
   type (type_interior_standard_variable), parameter :: calcite_saturation_state = type_interior_standard_variable(name='calcite_saturation_state', units='-')

   type (type_universal_standard_variable), parameter :: calcite_production = type_universal_standard_variable(name='calcite_production', units='mol m-3 s-1', aggregate_variable=.true.)
   type (type_universal_standard_variable), parameter :: zooplankton_production = type_universal_standard_variable(name='gross_mole_production_of_biomass_expressed_as_carbon_by_zooplankton', units='mol m-3 s-1', aggregate_variable=.true.)

   ! Maximum number of iterations for each method
   INTEGER, PARAMETER :: jp_maxniter_atgen    = 20
   REAL(rk), PARAMETER :: pp_rdel_ah_target = 1.E-4_rk

contains

   PURE FUNCTION solve_at_general( rhop, dic, tal, po4, sil, borat, sulfat, fluorid, &
      ak13, ak23, akb3, akw3, aks3, akf3, ak1p3, ak2p3, ak3p3, aksi3, p_hini) result(zhi) bind(c)
!$omp declare simd(solve_at_general)
   ! Universal pH solver that converges from any given initial value,
   ! determines upper an lower bounds for the solution if required

   ! Argument variables
   !--------------------
   REAL(rk), INTENT(IN), VALUE   :: rhop, dic, tal, po4, sil, borat, sulfat, fluorid
   REAL(rk), INTENT(IN), VALUE   :: ak13, ak23, akb3, akw3, aks3, akf3, ak1p3, ak2p3, ak3p3, aksi3
   REAL(rk), INTENT(IN), VALUE   :: p_hini
   REAL(rk)  :: zhi

   ! Local variables
   !-----------------
   INTEGER   ::  jn
   REAL(rk)  ::  zh_ini, zh, zh_prev, zh_lnfactor
   REAL(rk)  ::  zdelta, zh_delta
   REAL(rk)  ::  zeqn, zdeqndh, zalka
   REAL(rk)  ::  aphscale
   REAL(rk)  ::  znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
   REAL(rk)  ::  znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
   REAL(rk)  ::  znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
   REAL(rk)  ::  znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
   REAL(rk)  ::  znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
   REAL(rk)  ::  znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
   REAL(rk)  ::  zalk_wat, zdalk_wat
   REAL(rk)  ::  zfact, p_alktot, zdic, zbot, zpt, zst, zft, zsit
   LOGICAL   ::  l_exitnow
   REAL(rk), PARAMETER :: pz_exp_threshold = 1.0
   REAL(rk) :: zalknw_inf, zalknw_sup, rmask, zh_min, zh_max, zeqn_absmin

   zalknw_inf =  -po4 * 1000. / (rhop + rtrn) - sulfat  &
   &              - fluorid
   zalknw_sup =   (2. * dic + 2. * po4 + sil )    &
   &               * 1000. / (rhop + rtrn) + borat

   zhi   = 0.

   ! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
   p_alktot = tal * 1000. / (rhop + rtrn)
   aphscale = 1. + sulfat/aks3
   zh_ini = p_hini
   if (zh_ini < 0._rk) zh_ini = ahini_for_at(rhop, dic, tal, borat, ak13, ak23, akb3)

   zdelta = (p_alktot-zalknw_inf)**2 + 4.*akw3/aphscale

   IF(p_alktot >= zalknw_inf) THEN
      zh_min = 2.*akw3 /( p_alktot-zalknw_inf + SQRT(zdelta) )
   ELSE
      zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2.
   ENDIF

   zdelta = (p_alktot-zalknw_sup)**2 + 4.*akw3/aphscale

   IF(p_alktot <= zalknw_sup) THEN
      zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2.
   ELSE
      zh_max = 2.*akw3 /( p_alktot-zalknw_sup + SQRT(zdelta) )
   ENDIF

   zhi = MAX(MIN(zh_max, zh_ini), zh_min)

   zeqn_absmin = HUGE(1._rk)

   rmask = 1.

   DO jn = 1, jp_maxniter_atgen 
      IF (rmask == 1.) THEN 
         zfact = rhop / 1000. + rtrn
         p_alktot = tal / zfact
         zdic  = dic / zfact
         zbot  = borat
         zpt = po4 / zfact * po4r
         zsit = sil / zfact
         zst = sulfat
         zft = fluorid
         aphscale = 1. + sulfat/aks3
         zh = zhi
         zh_prev = zh

         ! H2CO3 - HCO3 - CO3 : n=2, m=0
         znumer_dic = 2.*ak13*ak23 + zh*ak13
         zdenom_dic = ak13*ak23 + zh*(ak13 + zh)
         zalk_dic   = zdic * (znumer_dic/zdenom_dic)
         zdnumer_dic = ak13*ak13*ak23 + zh     &
                        *(4.*ak13*ak23 + zh*ak13)
         zdalk_dic   = -zdic*(zdnumer_dic/zdenom_dic**2)


         ! B(OH)3 - B(OH)4 : n=1, m=0
         znumer_bor = akb3
         zdenom_bor = akb3 + zh
         zalk_bor   = zbot * (znumer_bor/zdenom_bor)
         zdnumer_bor = akb3
         zdalk_bor   = -zbot*(zdnumer_bor/zdenom_bor**2)


         ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
         znumer_po4 = 3.*ak1p3*ak2p3*ak3p3  &
         &            + zh*(2.*ak1p3*ak2p3 + zh* ak1p3)
         zdenom_po4 = ak1p3*ak2p3*ak3p3     &
         &            + zh*( ak1p3*ak2p3 + zh*(ak1p3 + zh))
         zalk_po4   = zpt * (znumer_po4/zdenom_po4 - 1.) ! Zero level of H3PO4 = 1
         zdnumer_po4 = ak1p3*ak2p3*ak1p3*ak2p3*ak3p3  &
         &             + zh*(4.*ak1p3*ak1p3*ak2p3*ak3p3         &
         &             + zh*(9.*ak1p3*ak2p3*ak3p3                         &
         &             + ak1p3*ak1p3*ak2p3                                &
         &             + zh*(4.*ak1p3*ak2p3 + zh * ak1p3 ) ) )
         zdalk_po4   = -zpt * (zdnumer_po4/zdenom_po4**2)

         ! H4SiO4 - H3SiO4 : n=1, m=0
         znumer_sil = aksi3
         zdenom_sil = aksi3 + zh
         zalk_sil   = zsit * (znumer_sil/zdenom_sil)
         zdnumer_sil = aksi3
         zdalk_sil   = -zsit * (zdnumer_sil/zdenom_sil**2)

         ! HSO4 - SO4 : n=1, m=1
         aphscale = 1.0 + zst/aks3
         znumer_so4 = aks3 * aphscale
         zdenom_so4 = aks3 * aphscale + zh
         zalk_so4   = zst * (znumer_so4/zdenom_so4 - 1.)
         zdnumer_so4 = aks3
         zdalk_so4   = -zst * (zdnumer_so4/zdenom_so4**2)

         ! HF - F : n=1, m=1
         znumer_flu =  akf3
         zdenom_flu =  akf3 + zh
         zalk_flu   =  zft * (znumer_flu/zdenom_flu - 1.)
         zdnumer_flu = akf3
         zdalk_flu   = -zft * (zdnumer_flu/zdenom_flu**2)

         ! H2O - OH
         aphscale = 1.0 + zst/aks3
         zalk_wat   = akw3/zh - zh/aphscale
         zdalk_wat  = -akw3/zh**2 - 1./aphscale

         ! CALCULATE [ALK]([CO3--], [HCO3-])
         zeqn = zalk_dic + zalk_bor + zalk_po4 + zalk_sil   &
         &      + zalk_so4 + zalk_flu                       &
         &      + zalk_wat - p_alktot

         zalka = p_alktot - (zalk_bor + zalk_po4 + zalk_sil   &
         &       + zalk_so4 + zalk_flu + zalk_wat)

         zdeqndh = zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
         &         + zdalk_so4 + zdalk_flu + zdalk_wat

         ! Adapt bracketing interval
         IF(zeqn > 0._rk) THEN
            zh_min = zh_prev
         ELSEIF(zeqn < 0._rk) THEN
            zh_max = zh_prev
         ENDIF

         IF(ABS(zeqn) >= 0.5_rk*zeqn_absmin) THEN
         ! if the function evaluation at the current point is
         ! not decreasing faster than with a bisection step (at least linearly)
         ! in absolute value take one bisection step on [ph_min, ph_max]
         ! ph_new = (ph_min + ph_max)/2d0
         !
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_min + ph_max)/2d0)
         !         = SQRT(10**(-(ph_min + phmax)))
         !         = SQRT(zh_max * zh_min)
            zh = SQRT(zh_max * zh_min)
            zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
         ELSE
         ! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
         !           = -zdeqndh * LOG(10) * [H]
         ! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))
         !
         ! pH_new = pH_old + \deltapH
         !
         ! [H]_new = 10**(-pH_new)
         !         = 10**(-pH_old - \Delta pH)
         !         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
         !         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
         !         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

            zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

            IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
               zh          = zh_prev*EXP(zh_lnfactor)
            ELSE
               zh_delta    = zh_lnfactor*zh_prev
               zh          = zh_prev + zh_delta
            ENDIF

            IF( zh < zh_min ) THEN
               ! if [H]_new < [H]_min
               ! i.e., if ph_new > ph_max then
               ! take one bisection step on [ph_prev, ph_max]
               ! ph_new = (ph_prev + ph_max)/2d0
               ! In terms of [H]_new:
               ! [H]_new = 10**(-ph_new)
               !         = 10**(-(ph_prev + ph_max)/2d0)
               !         = SQRT(10**(-(ph_prev + phmax)))
               !         = SQRT([H]_old*10**(-ph_max))
               !         = SQRT([H]_old * zh_min)
               zh                = SQRT(zh_prev * zh_min)
               zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
            ENDIF

            IF( zh > zh_max ) THEN
               ! if [H]_new > [H]_max
               ! i.e., if ph_new < ph_min, then
               ! take one bisection step on [ph_min, ph_prev]
               ! ph_new = (ph_prev + ph_min)/2d0
               ! In terms of [H]_new:
               ! [H]_new = 10**(-ph_new)
               !         = 10**(-(ph_prev + ph_min)/2d0)
               !         = SQRT(10**(-(ph_prev + ph_min)))
               !         = SQRT([H]_old*10**(-ph_min))
               !         = SQRT([H]_old * zhmax)
               zh                = SQRT(zh_prev * zh_max)
               zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
            ENDIF
         ENDIF

         zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)

         ! Stop iterations once |\delta{[H]}/[H]| < rdel
         ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
         ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

         ! Alternatively:
         ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
         !             ~ 1/LOG(10) * |\Delta [H]|/[H]
         !             < 1/LOG(10) * rdel

         ! Hence |zeqn/(zdeqndh*zh)| < rdel

         ! rdel <-- pp_rdel_ah_target
         l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

         IF(l_exitnow) THEN 
            rmask = 0.
         ENDIF

         zhi =  zh

         IF(jn >= jp_maxniter_atgen) THEN
            zhi = -1._rk
         ENDIF
      ENDIF
   END DO
   !

   END FUNCTION solve_at_general

   PURE FUNCTION ahini_for_at(rhop, dic, tal, borat, ak13, ak23, akb3) RESULT(p_hini) BIND(C)
!$omp declare simd(ahini_for_at)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE ahini_for_at  ***
      !!
      !! Subroutine returns the root for the 2nd order approximation of the
      !! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic 
      !! polynomial) around the local minimum, if it exists.
      !! Returns * 1E-03_wp if p_alkcb <= 0
      !!         * 1E-10_wp if p_alkcb >= 2*p_dictot + p_bortot
      !!         * 1E-07_wp if 0 < p_alkcb < 2*p_dictot + p_bortot
      !!                    and the 2nd order approximation does not have 
      !!                    a solution
      !!---------------------------------------------------------------------
      REAL(rk), INTENT(IN), VALUE   ::  rhop, dic, tal, borat, ak13, ak23, akb3
      REAL(rk)  ::  p_hini
      REAL(rk)  ::  zca1, zba1
      REAL(rk)  ::  zd, zsqrtd, zhmin
      REAL(rk)  ::  za2, za1, za0
      REAL(rk)  ::  p_dictot, p_bortot, p_alkcb 
      !!---------------------------------------------------------------------

      p_alkcb  = tal * 1000. / (rhop + rtrn)
      p_dictot = dic * 1000. / (rhop + rtrn)
      p_bortot = borat
      IF (p_alkcb <= 0.) THEN
         p_hini = 1.e-3
      ELSEIF (p_alkcb >= (2.*p_dictot + p_bortot)) THEN
         p_hini = 1.e-10_rk
      ELSE
         zca1 = p_dictot/( p_alkcb + rtrn )
         zba1 = p_bortot/ (p_alkcb + rtrn )
   ! Coefficients of the cubic polynomial
         za2 = aKb3*(1. - zba1) + ak13*(1.-zca1)
         za1 = ak13*akb3*(1. - zba1 - zca1)    &
         &     + ak13*ak23*(1. - (zca1+zca1))
         za0 = ak13*ak23*akb3*(1. - zba1 - (zca1+zca1))
                                 ! Taylor expansion around the minimum
         zd = za2*za2 - 3.*za1   ! Discriminant of the quadratic equation
                                 ! for the minimum close to the root

         IF(zd > 0.) THEN        ! If the discriminant is positive
            zsqrtd = SQRT(zd)
            IF(za2 < 0) THEN
               zhmin = (-za2 + zsqrtd)/3.
            ELSE
               zhmin = -za1/(za2 + zsqrtd)
            ENDIF
            p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
         ELSE
            p_hini = 1.e-7
         ENDIF
      !
      ENDIF
   END FUNCTION ahini_for_at

end module