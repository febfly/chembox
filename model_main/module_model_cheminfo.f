!Data module that stores information of chemical mechansim
      module module_model_cheminfo
      use module_model_parameter, only : DP,MAX_NSPEC,MAX_NRXN,MAX_STR1,
     +                                   MAX_STR2,MAX_NPROD,MAX_NREAC
      implicit none

      integer       :: nspec_all, nspec_active, nspec_inactive
      character(len=MAX_STR1),dimension(MAX_NSPEC) :: spec_name
      character(len=1)       ,dimension(MAX_NSPEC) :: spec_flag
      real(kind=DP)          ,dimension(MAX_NSPEC) :: spec_default_conc
      
      integer       :: nrxn_all, nrxn_photo
      integer, dimension(MAX_NRXN)  :: rxn_nreac, rxn_nprod 
      integer, dimension(MAX_NREAC, MAX_NRXN) :: rxn_reac
      integer, dimension(MAX_NPROD, MAX_NRXN) :: rxn_prod
      real(kind=DP), dimension(MAX_NREAC, MAX_NRXN) :: rxn_reac_coef 
      real(kind=DP), dimension(MAX_NPROD, MAX_NRXN) :: rxn_prod_coef 
      character(len=1),dimension(MAX_NRXN)          :: rxn_flag

      
      endmodule
