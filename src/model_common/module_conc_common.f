      module module_conc_common
      use module_model_parameter,only : DP, NIJK, MAX_NSPEC
      
      real(kind=DP) :: gas_conc(NIJK,MAX_NSPEC)
     
      endmodule module_conc_common
