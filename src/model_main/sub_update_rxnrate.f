      subroutine update_rxnrate
      use module_model_option,only: option_chemmech
      use module_met_common,only: temperature

      use module_geos_rxntype,only: geos_rxn_rate=>rxn_rate
      use module_ream_rxntype,only: ream_rxn_rate=>rxn_rate
      implicit none

      !
      endsubroutine update_rxnrate
