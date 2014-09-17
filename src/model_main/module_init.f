       module module_init
       implicit none
       public :: do_initialize
       contains

       subroutine do_initialize
       use module_model_parameter,only: DP, TS_CHEM_MIN
       use module_domain_common,only:calc_domain_index, NIJK
       use module_time_common,only:time_step_init
       use module_model_option,only:option_chemmech,option_solver
       use module_geoschem_io,only:read_geoschem
       use module_ream_io,only:read_ream
       use module_chemmech_common,only:spec_defconc,spec_getid
       use module_chemmech_common,only:nspec,nrxn,nactive,spec_name,
     +     nreac,reac_id,prod_id,prod_coefs,nphotorxn
       use module_chemmech_common,only:MAX_NREAC, MAX_NPROD, MAX_STR1
       use module_smvgear_interface,only:smvgear_setup       
       use module_conc_common,only:gas_conc
       use module_met_common,only:airdensity
       character(len=255) :: filename
       integer            :: i,j,k,ido3,idno2,idisop
       integer            :: ig, irxn
       real(kind=DP), dimension(MAX_NREAC, MAX_NRXN) :: reac_coefs
       real(kind=DP)  :: chem_tstep

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
       
       !Initialize chemical solver
       if (opt_solver.eq.1) then !smvgear
          reac_coefs=0d0
          do irxn=1,nrxn
             reac_coefs(1:nreac(irxn),irxn)=1d0
          enddo
          chem_tstep = TS_CHEM_MIN * 60d0 !unit: min to second
          call smvgear_setup(1,chem_tstep,nspec,nactive, MAX_STR1,
     +         spec_name, nrxn, nphotorxn, MAX_NREAC,MAX_NPROD,
     +         reac_id,prod_id,reac_coefs,prod_coefs)
       endif       
       endsubroutine do_initialize
       endmodule module_init
