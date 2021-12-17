#include "fabm_driver.h"

module pisces_aggregation
   use fabm_types
   use fabm_particle
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_aggregation
      type (type_state_variable_id) :: id_doc, id_poc, id_goc, id_sfe, id_bfe, id_conspoc, id_prodgoc
      type (type_dependency_id) :: id_xdiss
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_aggregation), intent(inout), target :: self
      integer,                         intent(in)            :: configunit

      call self%register_implemented_routines((/source_do/))

      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_state_dependency(self%id_goc, 'goc', 'mol C L-1', 'large particulate organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')

      call self%register_state_dependency(self%id_conspoc, 'conspoc', 'mol C L-1', 'consumed small particulate organic carbon')
      call self%register_state_dependency(self%id_prodgoc, 'prodgoc', 'mol C L-1', 'produced large particulate organic carbon')

      call self%request_coupling_to_model(self%id_doc, 'dom', 'c')
      call self%request_coupling_to_model(self%id_poc, 'pom', 'c')
      call self%request_coupling_to_model(self%id_goc, 'gom', 'c')
      call self%request_coupling_to_model(self%id_sfe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_bfe, 'gom', 'fe')
      call self%request_coupling_to_model(self%id_conspoc, 'pom', 'cons')
      call self%request_coupling_to_model(self%id_prodgoc, 'gom', 'prod')
      call self%register_dependency(self%id_xdiss, shear_rate)
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_aggregation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: doc, poc, goc, sfe, xdiss
      real(rk) :: zfact, zagg1, zagg2, zagg3, zagg4, zagg, zaggfe, zaggdoc, zaggdoc2, zaggdoc3

      _LOOP_BEGIN_

         _GET_(self%id_doc, doc)      ! dissolved organic carbon (mol C L-1)
         _GET_(self%id_poc, poc)      ! small particulate organic carbon (mol C L-1)
         _GET_(self%id_sfe, sfe)      ! small particulate organic iron (mol Fe L-1)
         _GET_(self%id_goc, goc)      ! large particulate organic carbon (mol C L-1)
         _GET_(self%id_xdiss, xdiss)  ! shear rate (s-1)

         zfact = xstep * xdiss
         !  Part I : Coagulation dependent on turbulence
         zagg1 = 25.9  * zfact * poc * poc    ! Jorn Eq 39, POC^2 shear term, value of 25.9 matches a6 in paper
         zagg2 = 4452. * zfact * poc * goc    ! Jorn Eq 39, POC*GOC shear term, value of 4452. matches a7 in paper

         ! Part II : Differential settling

         !  Aggregation of small into large particles
         zagg3 =  47.1 * xstep * poc * goc    ! Jorn Eq 39, value of 47.1 matches a9 in paper, but a8 expected per Eq 39
         zagg4 =  3.3  * xstep * poc * poc    ! Jorn Eq 39, value of 3.3 matches a8 in paper, but a9 expected per Eq 39

         zagg   = zagg1 + zagg2 + zagg3 + zagg4
         zaggfe = zagg * sfe / ( poc + rtrn )

         ! Aggregation of DOC to POC : 
         ! 1st term is shear aggregation of DOC-DOC
         ! 2nd term is shear aggregation of DOC-POC
         ! 3rd term is differential settling of DOC-POC
         zaggdoc  = ( ( 0.369 * 0.3 * doc + 102.4 * poc ) * zfact       &   ! Jorn: Eq 36a, values of 0.369 and 102.4 match a1=0.37 and a2=102 in paper, respectively
         &            + 2.4 * xstep * poc ) * 0.3 * doc                     ! Jorn: presumably POC part of Eq 36c, but value of 2.4 does not match a4=5095 given in paper
         ! transfer of DOC to GOC : 
         ! 1st term is shear aggregation
         ! 2nd term is differential settling 
         zaggdoc2 = ( 3.53E3 * zfact + 0.1 * xstep ) * goc * 0.3 * doc  ! Jorn: Eq 36b? value of 3.53E3 matches a3 in paper. GOC*DOC Brownian part with scale factor 0.1 not in paper
         ! tranfer of DOC to POC due to brownian motion
         zaggdoc3 =  114. * 0.3 * doc *xstep * 0.3 * doc                ! Jorn: Eq 36c, DOC part only, value of 114. matches a5 in paper

         !  Update the trends
         _ADD_SOURCE_(self%id_poc, - zagg + zaggdoc + zaggdoc3)
         _ADD_SOURCE_(self%id_goc, + zagg + zaggdoc2)
         _ADD_SOURCE_(self%id_sfe, - zaggfe)
         _ADD_SOURCE_(self%id_bfe, + zaggfe)
         _ADD_SOURCE_(self%id_doc, - zaggdoc - zaggdoc2 - zaggdoc3)
         !
         _ADD_SOURCE_(self%id_conspoc, - zagg + zaggdoc + zaggdoc3)
         _ADD_SOURCE_(self%id_prodgoc, + zagg + zaggdoc2)
      _LOOP_END_
   end subroutine

end module
