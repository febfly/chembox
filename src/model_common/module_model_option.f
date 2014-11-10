      module module_model_option
      implicit none
      integer, parameter :: option_solver=1

      !option_chemmech
      ! 1: Geos-chem
      ! 2: REAM
      ! 3: MCM
      integer, parameter :: option_chemmech=3

      character(len=255) :: path_input
      character(len=255) :: path_output

      endmodule module_model_option
