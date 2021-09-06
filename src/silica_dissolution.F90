#include "fabm_driver.h"

module pisces_silica_dissolution
   use fabm_types
   use pisces_shared

   implicit none

   private

   public type_pisces_silica_dissolution

   type, extends(type_base_model) :: type_pisces_silica_dissolution
      type (type_state_variable_id)      :: id_gsi, id_sil
      type (type_diagnostic_variable_id) :: id_remin
      type (type_dependency_id)          :: id_zsiremin   ! dissoluton rate computed by type_dissolution_rate module below
   contains
      procedure :: initialize
      procedure :: do
   end type

   type, extends(type_base_model) :: type_dissolution_rate
      type (type_diagnostic_variable_id) :: id_zsiremin
      type (type_dependency_id)          :: id_tem, id_gdept_n, id_ws, id_sil, id_e3t_n
      type (type_surface_dependency_id)  :: id_hmld, id_heup_01
      real(rk) :: xsilab, xsiremlab, xsirem
   contains
      procedure :: initialize => dissolution_rate_initialize
      procedure :: do_column  => dissolution_rate_do_column
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_silica_dissolution), intent(inout), target :: self
      integer,                                intent(in)            :: configunit

      class (type_dissolution_rate), pointer :: dissolution_rate

      allocate(dissolution_rate)
      call self%get_parameter(dissolution_rate%xsirem, 'xsirem', 'd-1', 'remineralization rate of Si', default=0.003_rk)
      call self%get_parameter(dissolution_rate%xsiremlab, 'xsiremlab', 'd-1', 'fast remineralization rate of Si', default=0.03_rk)
      call self%get_parameter(dissolution_rate%xsilab, 'xsilab', '1', 'labile fraction of biogenic silica', default=0.5_rk)
      call self%add_child(dissolution_rate, 'rate')

      call self%register_diagnostic_variable(self%id_remin, 'remin', 'mol Si L-1 s-1', 'rate')

      call self%register_state_dependency(self%id_gsi, 'gsi', 'mol Si L-1', 'particulate organic silicon')
      call self%register_state_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')
      call self%register_dependency(self%id_zsiremin, 'zsiremin', 's-1', 'rate')
      call dissolution_rate%request_coupling(dissolution_rate%id_sil, '../sil')
      call self%request_coupling(self%id_zsiremin, 'rate/zsiremin')
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_silica_dissolution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: gsi, zsiremin, zosil

      _LOOP_BEGIN_
         _GET_(self%id_gsi, gsi)
         _GET_(self%id_zsiremin, zsiremin)

         zosil    = zsiremin * gsi
         !
         _ADD_SOURCE_(self%id_gsi, - zosil)
         _ADD_SOURCE_(self%id_sil, + zosil)
         _SET_DIAGNOSTIC_(self%id_remin, zosil)
      _LOOP_END_
   end subroutine

   subroutine dissolution_rate_initialize(self, configunit)
      class (type_dissolution_rate), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      call self%register_diagnostic_variable(self%id_zsiremin, 'zsiremin', 's-1', 'specific rate', source=source_do_column)

      call self%register_dependency(self%id_ws, 'ws', 'm d-1', 'sinking velocity of particulate organic silicon')
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_heup_01, 'heup_01', 'm', 'depth where daily mean PAR equals 0.5 W m-2')
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      call self%register_dependency(self%id_tem, standard_variables%temperature)    ! Jorn: TODO should be in-situ temperature
      call self%register_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')
      call self%register_dependency(self%id_e3t_n, standard_variables%cell_thickness)
   end subroutine

   subroutine dissolution_rate_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_dissolution_rate), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: heup_01, hmld, zdep
      real(rk) :: gdept_n, tem, ws, sil, e3t_n
      real(rk) :: ztkel, sio3eq, zsatur, zsatur2, znusil, zfacsib, zfacsi, zsiremin

      _GET_SURFACE_(self%id_heup_01, heup_01)
      _GET_SURFACE_(self%id_hmld, hmld)
      zdep = MAX( hmld, heup_01 )

      zfacsi = self%xsilab                           ! Jorn: labile fraction of frustule (1)
      zfacsib = self%xsilab / ( 1.0 - self%xsilab )  ! Jorn: ratio of labile : non-labile fractions of frustule (-)
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_gdept_n, gdept_n)
         _GET_(self%id_tem, tem)
         _GET_(self%id_ws, ws)
         _GET_(self%id_sil, sil)
         _GET_(self%id_e3t_n, e3t_n)
         ztkel = tem + 273.15_rk    ! Jorn: TODO this is in-situ temperature in the original PISCES code
         sio3eq = EXP(  LOG( 10.) * ( 6.44 - 968. / ztkel )  ) * 1.e-6
         zsatur   = MAX( rtrn, ( sio3eq - sil ) / ( sio3eq + rtrn ) )
         zsatur2  = ( 1. + tem / 400.)**37
         znusil   = 0.225  * ( 1. + tem / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25   ! Jorn Eq 52 except for lambda_PSi scale factor. Outer exponent now 9.25 (9 in paper)
         ! Remineralization rate of BSi depedant on T and saturation
         ! ---------------------------------------------------------
         IF ( gdept_n > zdep ) THEN
            zfacsib = zfacsib * EXP( -0.5 * ( self%xsiremlab - self%xsirem )  &
            &                   * znusil * e3t_n / ws )
            zfacsi  = zfacsib / ( 1.0 + zfacsib )
            zfacsib = zfacsib * EXP( -0.5 * ( self%xsiremlab - self%xsirem )    &
            &                   * znusil * e3t_n / ws )
         ENDIF
         zsiremin = ( self%xsiremlab * zfacsi + self%xsirem * ( 1. - zfacsi ) ) * xstep * znusil
         _SET_DIAGNOSTIC_(self%id_zsiremin, zsiremin)
      _DOWNWARD_LOOP_END_
   end subroutine

end module