#include "fabm_driver.h"

module pisces_calcite_dissolution
   use fabm_types
   use pisces_shared

   implicit none

   private

   public type_pisces_calcite_dissolution

   type, extends(type_base_model) :: type_pisces_calcite_dissolution
      type (type_state_variable_id)      :: id_cal, id_dic, id_tal
      type (type_diagnostic_variable_id) :: id_dcal
      type (type_dependency_id)          :: id_zomegaca
      real(rk) :: kdca, nca
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_calcite_dissolution), intent(inout), target :: self
      integer,                                 intent(in)            :: configunit

      call self%register_implemented_routines((/source_do/))

      call self%get_parameter(self%kdca, 'kdca', 'month-1', 'dissolution rate constant', default=6._rk)   ! 0.197 d-1 in paper
      call self%get_parameter(self%nca, 'nca', '-', 'exponent in the dissolution rate', default=1._rk)

      call self%register_diagnostic_variable(self%id_dcal, 'DCAL', 'mol C m-3 s-1', 'rate')

      call self%register_state_dependency(self%id_cal, 'cal', 'mol C L-1', 'calcite')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_dependency(self%id_zomegaca, calcite_saturation_state)
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_calcite_dissolution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: cal, zomegaca, excess, zexcess0, zexcess, zdispot, zcaldiss
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
         _SET_DIAGNOSTIC_(self%id_dcal, zcaldiss * 1.e3_rk)
      _LOOP_END_
   end subroutine

end module