      module module_model_parameter
      implicit none
      !Numerical parameters
      integer, parameter :: DP = 8

      !Grid parameters
      integer, parameter :: IX = 1
      integer, parameter :: JX = 1
      integer, parameter :: KX = 1
      integer, parameter :: NIJ  = IX*JX
      integer, parameter :: NIJK = IX*JX*KX 

      !Chemical mechanism parameters

      !For REAM or GeosChem
      !integer, parameter :: MAX_NSPEC   = 300
      !integer, parameter :: MAX_STR1    = 15  !For names,e.g., species
      !integer, parameter :: MAX_STR2    = 50  !For long names
      !integer, parameter :: MAX_NREAC   = 4
      !integer, parameter :: MAX_NPROD   = 20

      !integer, parameter :: MAX_NRXN    = 500
      !integer, parameter :: MAX_NPHOTO  = 90

      !For MCM
      integer, parameter :: MAX_NSPEC   = 2000
      integer, parameter :: MAX_STR1    = 15  !For names,e.g., species
      integer, parameter :: MAX_STR2    = 50  !For long names
      integer, parameter :: MAX_NREAC   = 3
      integer, parameter :: MAX_NPROD   = 6

      integer, parameter :: MAX_NRXN    = 5000
      integer, parameter :: MAX_NPHOTO  = 35

      !Solver paramters
      integer, parameter :: MAX_BLK = 15

      !Time step
      real(kind=DP),parameter :: TS_MODEL_MIN = 1
      real(kind=DP),parameter :: TS_CHEM_MIN = 1
      real(kind=DP),parameter :: TS_OUTPUT_MIN = 2.
 

      endmodule module_model_parameter
