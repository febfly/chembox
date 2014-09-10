      module module_ream_parameter
      use module_model_parameter,only : MDP => DP,MMAX_NSPEC=>MAX_NSPEC,
     +                         MMAX_NRXN=>MAX_NRXN, MMAX_STR1=>MAX_STR1,
     +                      MMAX_NREAC=>MAX_NREAC,MMAX_NPROD=>MAX_NPROD
      implicit none
      integer,parameter :: DP        = MDP
      integer,parameter :: MAX_NSPEC = MMAX_NSPEC
      integer,parameter :: MAX_NRXN  = MMAX_NRXN
      integer,parameter :: MAX_STR1  = MMAX_STR1
      integer,parameter :: MAX_NREAC = MMAX_NREAC
      integer,parameter :: MAX_NPROD = MMAX_NPROD
!      integer,parameter :: DP        = 8
!      integer,parameter :: MAX_NSPEC = 300
!      integer,parameter :: MAX_NRXN  = 300
!      integer,parameter :: MAX_STR1  = 50
!      integer,parameter :: MAX_NREAC = 4
!      integer,parameter :: MAX_NPROD = 20
      
      endmodule module_ream_parameter
