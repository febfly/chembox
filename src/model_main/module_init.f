       module module_init
       implicit none
       public :: do_initialize
       contains
       subroutine do_initialize
       use module_domain_common,only:calc_domain_index, NIJK
       use module_time_common,only:time_step_init
       use module_model_option,only:option_chemmech,option_solver
       use module_geoschem_io,only:read_geoschem
       use module_ream_io,only:read_ream
       use module_chemmech_common,only:spec_defconc
       use module_conc_common,only:conc
       character(len=255) :: filename
       integer            :: ig

       call calc_domain_index
       call time_step_init

       if (option_chemmech.eq.1) then
          filename='globchem.dat'
          call read_geoschem(trim(filename))
       elseif (option_chemmech.eq.2) then
          filename='chem.dat'
          call read_ream(trim(filename))
       endif

       do ig = 1, NIJK
          conc(:,ig)  = spec_defconc(:)*2.5d19
       enddo
       
       endsubroutine do_initialize
       endmodule module_init
