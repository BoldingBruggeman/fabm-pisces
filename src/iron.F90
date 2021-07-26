#include "fabm_driver.h"

module pisces_iron
   use fabm_types
   use fabm_particle
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_iron
      type (type_state_variable_id) :: id_fer, id_sfe, id_bfe
      type (type_dependency_id) :: id_tempis, id_salinprac, id_xdiss, id_doc, id_poc, id_goc, id_cal, id_gsi, id_hi, id_oxy, id_etot, id_gdept_n, id_zdust
      type (type_surface_dependency_id) :: id_gphit
      type (type_diagnostic_variable_id) :: id_scav, id_coll, id_Fe3, id_FeL1, id_zTL1
      real(rk) :: ligand, xlam1, xlamdust, kfep, wdust
   contains
      procedure :: initialize
      procedure :: do
   end type

   logical, parameter :: ln_ligvar = .false.
   logical, parameter :: ln_ligand = .false.

contains

   subroutine initialize(self, configunit)
      class (type_pisces_iron), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%get_parameter(self%ligand, 'ligand', 'mol L-1', 'total concentration of iron ligands', default=0.7E-9_rk)
      call self%get_parameter(self%xlam1, 'xlam1', 'd-1 umol-1 L', 'scavenging rate', default=0.005_rk)
      call self%get_parameter(self%xlamdust, 'xlamdust', 'd-1 mg-1 L', 'scavenging rate of dust', default=150.0_rk)
      call self%get_parameter(self%kfep, 'kfep', 'd-1', 'nanoparticle formation rate constant', default=0.01_rk)

      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')

      call self%register_diagnostic_variable(self%id_scav, 'scav', 'mol Fe L-1 s-1', 'scavenging')
      call self%register_diagnostic_variable(self%id_coll, 'coll', 'mol Fe L-1 s-1', 'colloidal pumping of FeL')
      call self%register_diagnostic_variable(self%id_Fe3, 'Fe3', 'nmol Fe L-1', 'iron III concentration')
      call self%register_diagnostic_variable(self%id_FeL1, 'FeL1', 'nmol L-1', 'complexed Iron concentration with L1')
      call self%register_diagnostic_variable(self%id_zTL1, 'zTL1', 'nmol L-1', 'total L1 concentration')

      call self%register_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_dependency(self%id_goc, 'goc', 'mol C L-1', 'large particulate organic carbon')
      call self%register_dependency(self%id_cal, 'cal', 'mol C L-1', 'calcite')
      call self%register_dependency(self%id_gsi, 'gsi', 'mol Si L-1', 'large particulate organic silicon')

      call self%request_coupling_to_model(self%id_doc, 'dom', 'c')
      call self%request_coupling_to_model(self%id_poc, 'pom', 'c')
      call self%request_coupling_to_model(self%id_goc, 'gom', 'c')
      call self%request_coupling_to_model(self%id_sfe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_bfe, 'gom', 'fe')
      call self%request_coupling_to_model(self%id_cal, 'gom', 'cal')
      call self%request_coupling_to_model(self%id_gsi, 'gom', 'si')
      call self%register_dependency(self%id_xdiss, shear_rate)
      call self%register_dependency(self%id_hi, 'hi', 'mol L-1', 'hydrogen ion concentration')
      call self%register_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_dependency(self%id_etot, 'etot', 'W m-2', 'instantaneous PAR')
      call self%register_dependency(self%id_gphit, standard_variables%latitude)
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_tempis, standard_variables%temperature) ! TODO should be in-situ temperature (as opposed to conservative/potential)
      call self%register_dependency(self%id_salinprac, standard_variables%practical_salinity)
      call self%register_dependency(self%id_zdust, 'zdust', 'g m-2', 'dust concentration')
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_iron), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: fer, doc, poc, goc, cal, gsi, hi, oxy, etot, xdiss, gphit, tempis, salinprac, gdept_n
      real(rk) :: ztkel, zsal, zis, fekeq, ztkel1, fesol(5)
      real(rk) :: ztotlig, zTL1, zkeq, zfesatur, ztfe, zFe3, zFeL1, zdust, zhplus, fe3sol, zfeequi, zfecoll, precip, ztrc
      real(rk) :: zxlam, zlam1a, zlam1b, zscave, zdenom1, zdenom2, zlamfac, zdep, zcoag, zaggdfea, zaggdfeb

      _LOOP_BEGIN_
         _GET_(self%id_fer, fer)
         _GET_(self%id_doc, doc)
         _GET_(self%id_poc, poc)
         _GET_(self%id_goc, goc)
         _GET_(self%id_cal, cal)
         _GET_(self%id_gsi, gsi)
         _GET_(self%id_hi, hi)
         _GET_(self%id_oxy, oxy)
         _GET_(self%id_etot, etot)
         _GET_(self%id_xdiss, xdiss)
         _GET_SURFACE_(self%id_gphit, gphit)
         _GET_(self%id_tempis, tempis)
         _GET_(self%id_salinprac, salinprac)
         _GET_(self%id_gdept_n, gdept_n)
         _GET_(self%id_zdust, zdust)

         ztkel = tempis + 273.15_rk
         zsal  = salinprac !(ji,jj,1) + ( 1.- tmask(ji,jj,1) ) * 35.
         zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )

         fekeq  = 10**( 17.27 - 1565.7 / ztkel )

         ! Liu and Millero (1999) only valid 5 - 50 degC
         ztkel1 = MAX( 5. , tempis ) + 273.16
         fesol(1) = 10**(-13.486 - 0.1856* zis**0.5 + 0.3073*zis + 5254.0/ztkel1)
         fesol(2) = 10**(2.517 - 0.8885*zis**0.5 + 0.2139 * zis - 1320.0/ztkel1 )
         fesol(3) = 10**(0.4511 - 0.3305*zis**0.5 - 1996.0/ztkel1 )
         fesol(4) = 10**(-0.2965 - 0.7881*zis**0.5 - 4086.0/ztkel1 )
         fesol(5) = 10**(4.4466 - 0.8505*zis**0.5 - 7980.0/ztkel1 )

         ! Total ligand concentration : Ligands can be chosen to be constant or variable
         ! Parameterization from Tagliabue and Voelker (2011)
         ! -------------------------------------------------
         IF( ln_ligvar ) THEN
            ztotlig =  0.09 * doc * 1E6 + self%ligand * 1E9    ! Jorn: Eq 67 (note max operator in that Eq has no effect)
            ztotlig =  MIN( ztotlig, 10. )
         ELSE
           IF( ln_ligand ) THEN  ;   !ztotlig = lgw * 1E9      ! Jonr: TODO
           ELSE                  ;   ztotlig = self%ligand * 1E9
           ENDIF
         ENDIF

         ! ------------------------------------------------------------
         !  from Aumont and Bopp (2006)
         ! This model is based on one ligand and Fe' 
         ! Chemistry is supposed to be fast enough to be at equilibrium
         ! ------------------------------------------------------------
         zTL1  = ztotlig
         zkeq            = fekeq
         zfesatur        = zTL1 * 1E-9
         ztfe            = fer
         ! Fe' is the root of a 2nd order polynom
         zFe3 = ( -( 1. + zfesatur * zkeq - zkeq * ztfe )               &
            &              + SQRT( ( 1. + zfesatur * zkeq - zkeq * ztfe )**2       &
            &              + 4. * ztfe * zkeq) ) / ( 2. * zkeq )    ! Jorn: Eq65
         zFe3 = zFe3 * 1E9                    ! Jorn: free inorganic iron
         zFeL1 = MAX( 0., fer * 1E9 - zFe3 )  ! Jorn: "complexed" iron
         _SET_DIAGNOSTIC_(self%id_Fe3, zFe3)
         _SET_DIAGNOSTIC_(self%id_FeL1, zFeL1)
         _SET_DIAGNOSTIC_(self%id_zTL1, zTL1)

         ! Scavenging rate of iron. This scavenging rate depends on the load of particles of sea water. 
         ! This parameterization assumes a simple second order kinetics (k[Particles][Fe]).
         ! Scavenging onto dust is also included as evidenced from the DUNE experiments.
         ! --------------------------------------------------------------------------------------
         zhplus  = max( rtrn, hi )
         fe3sol  = fesol(1) * ( zhplus**3 + fesol(2) * zhplus**2  &
         &         + fesol(3) * zhplus + fesol(4)     &
         &         + fesol(5) / zhplus )
         !
         zfeequi = zFe3 * 1E-9
         zfecoll = 0.5 * zFeL1 * 1E-9   ! Jorn: Eq 66
         ! precipitation of Fe3+, creation of nanoparticles
         precip = MAX( 0., ( zFe3 * 1E-9 - fe3sol ) ) * self%kfep * xstep   ! Jorn: replaces Eq 62?
         !
         ztrc   = ( poc + goc + cal + gsi ) * 1.e6 
         IF (ln_ligand) THEN
            zxlam  = self%xlam1 * MAX( 1.E-3, EXP(-2 * etot / 10. ) * (1. - EXP(-2 * oxy / 100.E-6 ) ))
         ELSE
            zxlam  = self%xlam1 * 1.0
         ENDIF
         zlam1b = 3.e-5 + self%xlamdust * zdust + zxlam * ztrc    ! Jorn: Eq 50a, note lambda_fe^in is here hardcoded to 3.e-5 (same value as in paper)
         zscave = zfeequi * zlam1b * xstep                        ! Jorn: Eq 50b

         ! Compute the different ratios for scavenging of iron
         ! to later allocate scavenged iron to the different organic pools
         ! ---------------------------------------------------------
         zdenom1 = zxlam * poc / zlam1b
         zdenom2 = zxlam * goc / zlam1b

         !  Increased scavenging for very high iron concentrations found near the coasts 
         !  due to increased lithogenic particles and let say it is unknown processes (precipitation, ...)
         !  -----------------------------------------------------------
         zlamfac = MAX( 0.e0, ( gphit + 55.) / 30. )
         zlamfac = MIN( 1.  , zlamfac )
         zdep    = MIN( 1., 1000. / gdept_n )
         zcoag   = 1E-4 * ( 1. - zlamfac ) * zdep * xstep * fer

         !  Compute the coagulation of colloidal iron. This parameterization 
         !  could be thought as an equivalent of colloidal pumping.
         !  It requires certainly some more work as it is very poorly constrained.
         !  ----------------------------------------------------------------
         zlam1a   = ( 0.369  * 0.3 * doc + 102.4  * poc ) * xdiss    &
               &      + ( 114.   * 0.3 * doc )
         zaggdfea = zlam1a * xstep * zfecoll    ! Jorn: Eq 61a
         !
         zlam1b   = 3.53E3 * goc * xdiss
         zaggdfeb = zlam1b * xstep * zfecoll    ! Jorn: Eq 61b
         !
         _ADD_SOURCE_(self%id_fer, - zscave - zaggdfea - zaggdfeb - zcoag - precip)
         _ADD_SOURCE_(self%id_sfe, + zscave * zdenom1 + zaggdfea)
         _ADD_SOURCE_(self%id_bfe, + zscave * zdenom2 + zaggdfeb)
         _SET_DIAGNOSTIC_(self%id_scav, zscave)
         _SET_DIAGNOSTIC_(self%id_coll, zaggdfea + zaggdfeb)
         !
         !
         !  Define the bioavailable fraction of iron
         !  ----------------------------------------
         !biron = fer   ! Jorn: unused????

         !IF( ln_ligand ) THEN   ! Jorn: TODO, ln_ligand not implemented yet
         !   !
         !   zlam1a   = ( 0.369  * 0.3 * doc + 102.4  * poc ) * xdiss    &
         !         &    + ( 114.   * 0.3 * doc )
         !   !
         !   zlam1b   = 3.53E3 *   goc * xdiss
         !   zligco   = 0.5 * lgw
         !   zaggliga = zlam1a * xstep * zligco
         !   zaggligb = zlam1b * xstep * zligco
         !   _ADD_SOURCE_(self%id_lgw, - zaggliga - zaggligb)
         !   !zlcoll3d(ji,jj,jk)  = zaggliga + zaggligb
         !   !
         !   plig =  MAX( 0., ( ( zFeL1 * 1E-9 ) / ( fer +rtrn ) ) )
         !   !
         !ENDIF

      _LOOP_END_
   end subroutine

end module
