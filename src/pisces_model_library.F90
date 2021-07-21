module pisces_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use pisces_tracer
   use pisces_phytoplankton
   use pisces_zooplankton
   use pisces_optics
   use pisces_carbonate_chemistry
   use pisces_oxygen
   use pisces_daylength
   use pisces_turbocline
   use pisces_shear
   use pisces_remineralization
   use pisces_aggregation
   use pisces_iron
   use pisces_nitrification
   use pisces_nitrogen_fixation

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: pisces_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('tracer');              allocate(type_pisces_tracer::model)
         case ('phytoplankton');       allocate(type_pisces_phytoplankton::model)
         case ('zooplankton');         allocate(type_pisces_zooplankton::model)
         case ('optics');              allocate(type_pisces_optics::model)
         case ('carbonate_chemistry'); allocate(type_pisces_carbonate_chemistry::model)
         case ('oxygen');              allocate(type_pisces_oxygen::model)
         case ('daylength');           allocate(type_pisces_daylength::model)
         case ('turbocline');          allocate(type_pisces_turbocline::model)
         case ('shear');               allocate(type_pisces_shear::model)
         case ('remineralization');    allocate(type_pisces_remineralization::model)
         case ('aggregation');         allocate(type_pisces_aggregation::model)
         case ('nitrification');       allocate(type_pisces_nitrification::model)
         case ('nitrogen_fixation');   allocate(type_pisces_nitrogen_fixation::model)
         case ('iron');                allocate(type_pisces_iron::model)
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name, model)
      end select
   end subroutine create

end module
