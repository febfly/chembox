      module module_solver_driver
      use 
      implicit none
      public :: solver_init
      public :: solver_solve
      contains

      subroutine solver_init
      use module_chemmech_common,only:nspec,nrxn,nactive,spec_name,
     +    reac_id,prod_id,prod_coefs,nphotorxn
      use module_chemmech_common,only:MAX_NREAC, MAX_NPROD, MAX_STR1
use module_time_common,only:chem_tstep!second
      use module_smvgear_interface,only:smvgear_setup
      if (opt_solver.eq.1) then !smvgear
         call smvgear_setup(1,chem_tstep,nspec,nactive, MAX_STR1,
     +        spec_name, nrxn, nphotorxn, MAX_NREAC,MAX_NPROD,
     +        reac_id,prod_id,...,prod_coefs)
      endif
      endsubroutine

      subroutine solver_solve
      use module_chemmech_common,only:nspec,nrxn
      use module_conc_common,only:conc,
      logical,save ::firsttime=.true.
      
      

      call smvgear_solve(1, ifsun??,n_grid??,
      endsubroutine

      endmodule module_solver_driver
