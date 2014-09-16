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
       use module_chemmech_common,only:spec_defconc,spec_getid
       use module_conc_common,only:gas_conc
       use module_met_common,only:airdensity
       character(len=255) :: filename
       integer            :: ig,i,j,k,ido3,idno2,idisop

       !Read configuration
       !call read_config

       !Calculate index correspondence in 3D, 2D representation 
       !and 1D representation  
       !Not necessary for 1D model
       call calc_domain_index

       !Initialize time vairables
       call time_step_init(2011,10,10,0,0,0.,
     +                     2011,10,11,0,0,0.)

       !Initialize chemical mechanism
       if (option_chemmech.eq.1) then
          filename='globchem.dat'
          call read_geoschem(trim(filename))
       elseif (option_chemmech.eq.2) then
          filename='chem.dat'
          call read_ream(trim(filename))
       endif

       !Chemical initial condition
       ido3=spec_getid('O3')
       idno2=spec_getid('NO2')
       idisop=spec_getid('ISOP')
       do ig = 1, NIJK
          i=grid_1_3_i(ig)
          j=grid_1_3_j(ig)
          k=grid_1_3_k(ig)
          gas_conc(ig,:)  = spec_defconc(:)*airdensity(i,j,k)
          gas_conc(ig,ido3) = 40d-9*airdensity(i,j,k)
          gas_conc(ig,idno2) = 5d-9*airdensity(i,j,k)
          gas_conc(ig,idisop) = 1d-9*airdensity(i,j,k)
       enddo
       
       endsubroutine do_initialize
       endmodule module_init