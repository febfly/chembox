      module module_solver_driver
      use module_ream_cheminfo,only:

      use 
      implicit none
      public :: solver_init
      public :: solver_solve
      contains

      subroutine solver_init
      if (opt_solver.eq.1) then !smvgear

      endif
      endsubroutine

      subroutine solver_solve

      endsubroutine

      endmodule module_solver_driver
