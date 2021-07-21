#include "fabm_driver.h"

module pisces_tracer
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_tracer
      type (type_state_variable_id) :: id_c, id_fe, id_si, id_cal
      type (type_diagnostic_variable_id) :: id_ws
      real(rk) :: ws, wsmax, wsscale
   contains
      procedure :: initialize
   end type

   type, extends(type_base_model) :: type_pisces_calcite_dissolution
      type (type_state_variable_id)      :: id_cal, id_dic, id_tal
      type (type_diagnostic_variable_id) :: id_dcal
      type (type_dependency_id)          :: id_zomegaca
      real(rk) :: kdca, nca
   contains
      procedure :: initialize => calcite_dissolution_initialize
      procedure :: do         => calcite_dissolution_do
   end type

   type, extends(type_base_model) :: type_pisces_silicate_dissolution
      type (type_state_variable_id)      :: id_gsi, id_sil
      type (type_diagnostic_variable_id) :: id_remin
      type (type_dependency_id)          :: id_tem, id_gdept_n, id_ws
      type (type_surface_dependency_id)  :: id_hmld, id_heup_01
      real(rk) :: xsilab, xsiremlab, xsirem
   contains
      procedure :: initialize => silicate_dissolution_initialize
      procedure :: do         => silicate_dissolution_do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_tracer), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      character(len=64) :: name
      class (type_pisces_calcite_dissolution),  pointer :: calcite_dissolution
      class (type_pisces_silicate_dissolution), pointer :: silicate_dissolution

      ! Jorn: default initial values taken from trcice_pisces.F90, global prescribed concentrations
      call self%get_parameter(name, 'type', '', 'type of tracer (nitrate, ammonium, phosphate)')
      select case (name)
      case ('nitrate')
         call self%register_state_variable(self%id_c, 'c', 'mol C/L', 'concentration in carbon units', initial_value=5.79e-6_rk / rno3)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
      case ('ammonium')
         call self%register_state_variable(self%id_c, 'c', 'mol C/L', 'concentration in carbon units', initial_value=3.22e-7_rk / rno3)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
      case ('phosphate')
         call self%register_state_variable(self%id_c, 'c', 'mol C/L', 'concentration in carbon units', initial_value=5.77e-7_rk / po4r)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      case ('silicate')
         call self%register_state_variable(self%id_si, 'si', 'mol Si/L', 'concentration', initial_value=1.07e-7_rk)
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_si, scale_factor=1e6_rk)
      case ('iron')
         call self%register_state_variable(self%id_fe, 'fe', 'mol Fe/L', 'concentration', initial_value=4.06e-10_rk)
         call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)
      case ('dom')
         call self%register_state_variable(self%id_c, 'c', 'mol C/L', 'carbon concentration', initial_value=2.04e-5_rk)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      case ('pom')
         call self%get_parameter(self%ws, 'ws', 'm d-1', 'sinking velocity', default=2._rk)
         call self%register_diagnostic_variable(self%id_ws, 'ws', 'm d-1', 'sinking velocity', missing_value=self%ws, source=source_constant)

         call self%register_state_variable(self%id_c, 'c', 'mol C/L', 'carbon concentration', initial_value=1.27e-6_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
         call self%register_state_variable(self%id_fe, 'fe', 'mol Fe/L', 'iron concentration', initial_value=2.51e-11_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)
      case ('gom')
         call self%get_parameter(self%ws, 'ws', 'm d-1', 'minimum sinking velocity', default=50._rk)
         call self%get_parameter(self%wsmax, 'wsmax', 'm d-1', 'maximum sinking velocity', default=50._rk)
         call self%get_parameter(self%wsscale, 'wsscale', 'm', 'maximum sinking velocity', default=5000._rk)
         call self%register_diagnostic_variable(self%id_ws, 'ws', 'm d-1', 'sinking velocity', missing_value=self%ws, source=source_constant)

         call self%register_state_variable(self%id_c, 'c', 'mol C/L', 'carbon concentration', initial_value=5.23e-8_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
         call self%register_state_variable(self%id_fe, 'fe', 'mol Fe/L', 'iron concentration', initial_value=9.84e-13_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)
         call self%register_state_variable(self%id_si, 'si', 'mol Si/L', 'silicon concentration', initial_value=1.53e-8_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_si, scale_factor=1e6_rk)
         call self%register_state_variable(self%id_cal, 'cal', 'mol C/L', 'calcite concentration', initial_value=1.04e-8_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_cal, scale_factor=1e6_rk)

         allocate(calcite_dissolution)
         call self%add_child(calcite_dissolution, 'calcite_dissolution', configunit=-1)
         call calcite_dissolution%request_coupling(calcite_dissolution%id_cal, '../cal')

         allocate(silicate_dissolution)
         call self%add_child(silicate_dissolution, 'silicate_dissolution', configunit=-1)
         call silicate_dissolution%request_coupling(silicate_dissolution%id_gsi, '../si')
         call silicate_dissolution%request_coupling(silicate_dissolution%id_ws, '../ws')
         call silicate_dissolution%request_coupling(silicate_dissolution%id_sil, '../../sil/si')  ! Jorn: hack, TODO fix!!
         call silicate_dissolution%request_coupling(silicate_dissolution%id_heup_01, '../../optics/heup_01')  ! Jorn: hack, TODO fix!!
      case default
         call self%fatal_error('initialize', 'Unknown tracer type "' // trim(name) // '" requested.')
      end select
   end subroutine initialize

   subroutine calcite_dissolution_initialize(self, configunit)
      class (type_pisces_calcite_dissolution), intent(inout), target :: self
      integer,                                 intent(in)            :: configunit

      character(len=64) :: name

      call self%get_parameter(self%kdca, 'kdca', 'month-1', 'dissolution rate constant', default=6._rk)   ! 0.197 d-1 in paper
      call self%get_parameter(self%nca, 'nca', '-', 'exponent in the dissolution rate', default=1._rk)

      call self%register_diagnostic_variable(self%id_dcal, 'DCAL', 'mol C L-1 s-1', 'rate')

      call self%register_state_dependency(self%id_cal, 'cal', 'mol C L-1', 'calcite')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_dependency(self%id_zomegaca, calcite_saturation_state)
   end subroutine

   subroutine calcite_dissolution_do(self, _ARGUMENTS_DO_)
      class (type_pisces_calcite_dissolution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: cal, zomegaca, excess, zexcess0, zexcess, zdispot, zcaldiss
      real(rk), parameter :: nyear_len    =  365._rk     !: number of days in a year - Jorn: TODO this should depend on calendar, leap year
      real(rk), parameter :: raamo    =  12._rk          !: number of months in one year
      real(rk), parameter :: rmtss = nyear_len * rday / raamo

      _LOOP_BEGIN_
         _GET_(self%id_cal, cal)
         _GET_(self%id_zomegaca, zomegaca)

         ! Jorn: from p4zlys.F90, Eq 78, 79
         ! SET DEGREE OF UNDER-/SUPERSATURATION
         excess = 1._rk - zomegaca
         zexcess0 = MAX( 0., excess )
         zexcess  = zexcess0**self%nca

         ! AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
         !       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
         !       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
         zdispot = self%kdca * zexcess * cal
         !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
         !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
         zcaldiss  = zdispot / rmtss ! calcite dissolution   Jorn: dropped multiplication with rfact2 [time step in seconds]
         !
         _ADD_SOURCE_(self%id_tal, + 2. * zcaldiss)
         _ADD_SOURCE_(self%id_cal, -      zcaldiss)
         _ADD_SOURCE_(self%id_dic, +      zcaldiss)
         _SET_DIAGNOSTIC_(self%id_dcal, zcaldiss)
      _LOOP_END_
   end subroutine


   subroutine silicate_dissolution_initialize(self, configunit)
      class (type_pisces_silicate_dissolution), intent(inout), target :: self
      integer,                                  intent(in)            :: configunit

      character(len=64) :: name

      call self%get_parameter(self%xsirem, 'xsirem', 'd-1', 'remineralization rate of Si', default=0.003_rk)
      call self%get_parameter(self%xsiremlab, 'xsiremlab', 'd-1', 'fast remineralization rate of Si', default=0.03_rk)
      call self%get_parameter(self%xsilab, 'xsilab', '1', 'labile fraction of biogenic silica', default=0.5_rk)

      call self%register_diagnostic_variable(self%id_remin, 'remin', 'mol Si L-1 s-1', 'rate')

      call self%register_state_dependency(self%id_gsi, 'gsi', 'mol Si L-1', 'particulate organic silicon')
      call self%register_state_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')
      call self%register_dependency(self%id_ws, 'ws', 'm d-1', 'sinking velocity of particulate organic silicon')
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_heup_01, 'heup_01', 'm', 'depth where daily mean PAR equals 0.5 W m-2')
      call self%register_dependency(self%id_hmld, turbocline_depth)
      call self%register_dependency(self%id_tem, standard_variables%temperature)    ! Jorn: TODO should be in-situ temperature
   end subroutine

   subroutine silicate_dissolution_do(self, _ARGUMENTS_DO_)
      class (type_pisces_silicate_dissolution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sil, gsi, tem, ws, gdept_n, heup_01, hmld, wsbio4
      real(rk) :: ztkel, sio3eq, zdep, zsatur, zsatur2, znusil, zfacsib, zfacsi, zsiremin, zosil

      _LOOP_BEGIN_
         _GET_(self%id_sil, sil)
         _GET_(self%id_gsi, gsi)
         _GET_(self%id_gdept_n, gdept_n)
         _GET_(self%id_tem, tem)
         _GET_(self%id_ws, ws)
         _GET_SURFACE_(self%id_heup_01, heup_01)
         _GET_SURFACE_(self%id_hmld, hmld)

         ztkel = tem + 273.15_rk    ! Jorn: TODO this is in-situ temperature in the original PISCES code
         sio3eq = EXP(  LOG( 10.) * ( 6.44 - 968. / ztkel )  ) * 1.e-6
         zdep     = MAX( hmld, heup_01 )
         zsatur   = MAX( rtrn, ( sio3eq - sil ) / ( sio3eq + rtrn ) )
         zsatur2  = ( 1. + tem / 400.)**37
         znusil   = 0.225  * ( 1. + tem / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25   ! Jorn Eq 52 except for lambda_PSi scale factor. Outer exponent now 9.25 (9 in paper)
         ! Remineralization rate of BSi depedant on T and saturation
         ! ---------------------------------------------------------
         zfacsib = self%xsilab / ( 1.0 - self%xsilab ) * EXP( -(self%xsiremlab - self%xsirem )  &  ! Jorn: labile : non-labile ratio, but otherwise Eq 53a. TODO: use depth integration as in original code because znusil is depth dependent
         &                   * znusil * MAX(0._rk, gdept_n - zdep) / ws )
         zfacsi  = zfacsib / ( 1.0 + zfacsib )                                                     ! convert to labile : total ratio
         zsiremin = ( self%xsiremlab * zfacsi + self%xsirem * ( 1. - zfacsi ) ) * xstep * znusil   ! Jorn: Eq 53b
         zosil    = zsiremin * gsi
         !
         _ADD_SOURCE_(self%id_gsi, - zosil)
         _ADD_SOURCE_(self%id_sil, + zosil)
         _SET_DIAGNOSTIC_(self%id_remin, zosil)
      _LOOP_END_
   end subroutine

end module