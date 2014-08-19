      module module_model_parameter
      implicit none
      !Numerical parameters
      integer, parameter :: DP = 8

      !Chemical mechanism parameters
      integer, parameter :: MAX_NSPEC   = 300
      integer, parameter :: MAX_STR1    = 15  !For names,e.g., species
      integer, parameter :: MAX_STR2    = 50  !For long names
      !integer, parameter :: MAX_STR3    = 200 !For Formula Expression
      !integer, parameter :: MAX_NRXNTYPE= 30
      integer, parameter :: MAX_NREAC   = 4
      integer, parameter :: MAX_NPROD   = 20

      integer, parameter :: MAX_NRXN    = 500
      integer, parameter :: MAX_NPHOTO  = 90
 

      endmodule module_model_parameter