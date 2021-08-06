#include "fabm_driver.h"

module pisces_tracer
   use fabm_types
   use fabm_particle

   use pisces_shared
   use pisces_calcite_dissolution
   use pisces_silica_dissolution

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_tracer
      type (type_state_variable_id) :: id_c, id_fe, id_si, id_cal
      type (type_bottom_state_variable_id) :: id_bc, id_bfe, id_bsi, id_bcal
      type (type_diagnostic_variable_id) :: id_ws, id_prod, id_cons
      type (type_dependency_id) :: id_e3t_n
      real(rk) :: ws, wsmax, wsscale
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_tracer), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      character(len=64) :: name
      class (type_pisces_calcite_dissolution), pointer :: calcite_dissolution
      class (type_pisces_silica_dissolution),  pointer :: silica_dissolution

      ! Jorn: default initial values taken from trcice_pisces.F90, global prescribed concentrations
      call self%get_parameter(name, 'type', '', 'type of tracer (nitrate, ammonium, phosphate)')
      select case (name)
      case ('nitrate')
         call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'concentration in carbon units', initial_value=5.79e-6_rk / rno3)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
      case ('ammonium')
         call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'concentration in carbon units', initial_value=3.22e-7_rk / rno3)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
      case ('phosphate')
         call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'concentration in carbon units', initial_value=5.77e-7_rk / po4r)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      case ('silicate')
         call self%register_state_variable(self%id_si, 'si', 'mol Si/L', 'concentration', initial_value=1.07e-7_rk)
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_si, scale_factor=1e6_rk)
      case ('iron')
         call self%register_state_variable(self%id_fe, 'fe', 'mol Fe L-1', 'concentration', initial_value=4.06e-10_rk)
         call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)
      case ('dom')
         call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'carbon concentration', initial_value=2.04e-5_rk)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      case ('pom')
         call self%get_parameter(self%ws, 'ws', 'm d-1', 'sinking velocity', default=2._rk)
         call self%register_diagnostic_variable(self%id_ws, 'ws', 'm d-1', 'sinking velocity', missing_value=self%ws, source=source_constant)

         call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'carbon concentration', initial_value=1.27e-6_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
         call self%register_state_variable(self%id_fe, 'fe', 'mol Fe/L', 'iron concentration', initial_value=2.51e-11_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)

         call self%register_diagnostic_variable(self%id_prod, 'prod', 'mol C L-1', 'produced particulate organic carbon', &
            act_as_state_variable=.true., source=source_constant, output=output_none)
         call self%register_diagnostic_variable(self%id_cons, 'cons', 'mol C L-1', 'consumed particulate organic carbon', &
            act_as_state_variable=.true., source=source_constant, output=output_none)
      case ('gom')
         call self%get_parameter(self%ws, 'ws', 'm d-1', 'minimum sinking velocity', default=50._rk)
         call self%get_parameter(self%wsmax, 'wsmax', 'm d-1', 'maximum sinking velocity', default=50._rk)
         call self%get_parameter(self%wsscale, 'wsscale', 'm', 'maximum sinking velocity', default=5000._rk)
         call self%register_diagnostic_variable(self%id_ws, 'ws', 'm d-1', 'sinking velocity', missing_value=self%ws, source=source_constant)

         call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'carbon concentration', initial_value=5.23e-8_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
         call self%register_state_variable(self%id_fe, 'fe', 'mol Fe L-1', 'iron concentration', initial_value=9.84e-13_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_fe, scale_factor=1e9_rk)
         call self%register_state_variable(self%id_si, 'si', 'mol Si L-1', 'silicon concentration', initial_value=1.53e-8_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_si, scale_factor=1e6_rk)
         call self%register_state_variable(self%id_cal, 'cal', 'mol C L-1', 'calcite concentration', initial_value=1.04e-8_rk, vertical_movement=-self%ws * r1_rday)
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_cal, scale_factor=1e6_rk)

         allocate(calcite_dissolution)
         call self%add_child(calcite_dissolution, 'calcite_dissolution')
         call calcite_dissolution%request_coupling(calcite_dissolution%id_cal, '../cal')

         allocate(silica_dissolution)
         call self%add_child(silica_dissolution, 'silica_dissolution')
         call silica_dissolution%request_coupling(silica_dissolution%id_gsi, '../si')
         call silica_dissolution%request_coupling('rate/ws', '../ws')
         call silica_dissolution%request_coupling(silica_dissolution%id_sil, '../../sil/si')  ! Jorn: hack, TODO fix!!
         call silica_dissolution%request_coupling('rate/heup_01', '../../optics/heup_01')  ! Jorn: hack, TODO fix!!

         call self%register_diagnostic_variable(self%id_prod, 'prod', 'mol C L-1', 'produced particulate organic carbon', &
            act_as_state_variable=.true., source=source_constant, output=output_none)
         call self%register_diagnostic_variable(self%id_cons, 'cons', 'mol C L-1', 'consumed particulate organic carbon', &
            act_as_state_variable=.true., source=source_constant, output=output_none)
      case default
         call self%fatal_error('initialize', 'Unknown tracer type "' // trim(name) // '" requested.')
      end select

      if (self%ws /= 0._rk) then
         ! Couple to benthic pools to deposit sinking material in.
         ! Note that these already have FABM default units so that ultimately we can couple to other benthic models too.
         call self%register_state_dependency(self%id_bc, 'bc', 'mmol C m-2', 'target pool for organic carbon depositon')
         call self%register_state_dependency(self%id_bsi, 'bsi', 'mmol Si m-2', 'target pool for silicate depositon')
         call self%register_state_dependency(self%id_bfe, 'bfe', 'umol Fe m-2', 'target pool for iron depositon')
         call self%register_state_dependency(self%id_bcal, 'bcal', 'mmol C m-2', 'target pool for calcite depositon')
         call self%request_coupling_to_model(self%id_bc, 'sed', 'c')
         call self%request_coupling_to_model(self%id_bsi, 'sed', 'si')
         call self%request_coupling_to_model(self%id_bfe, 'sed', 'fe')
         call self%request_coupling_to_model(self%id_bcal, 'sed', 'cal')
         call self%register_dependency(self%id_e3t_n, standard_variables%cell_thickness)
      end if
   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_pisces_tracer), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: ws, e3t_n, c, fe, si, cal

      if (self%ws == 0._rk) return

      _BOTTOM_LOOP_BEGIN_
         _GET_(self%id_e3t_n, e3t_n)
         ws = min(0.99_rk * e3t_n / maxdt, self%ws * xstep)   ! Jorn: protect against CFL violation as in p4zsed.F90
         _GET_(self%id_c, c)
         _ADD_BOTTOM_FLUX_(self%id_c, -ws * c)
         _ADD_BOTTOM_SOURCE_(self%id_bc, ws * c * 1.e6_rk)
         _GET_(self%id_fe, fe)
         _ADD_BOTTOM_FLUX_(self%id_fe, -ws * fe)
         _ADD_BOTTOM_SOURCE_(self%id_bfe, ws * fe * 1.e9_rk)
         if (_VARIABLE_REGISTERED_(self%id_si)) then
            _GET_(self%id_si, si)
            _ADD_BOTTOM_FLUX_(self%id_si, -ws * si)
            _ADD_BOTTOM_SOURCE_(self%id_bsi, ws * si * 1.e6_rk)
         end if
         if (_VARIABLE_REGISTERED_(self%id_cal)) then
            _GET_(self%id_cal, cal)
            _ADD_BOTTOM_FLUX_(self%id_cal, -ws * cal)
            _ADD_BOTTOM_SOURCE_(self%id_bcal, ws * cal * 1.e6_rk)
         end if
      _BOTTOM_LOOP_END_
   end subroutine

end module