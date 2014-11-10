       module module_init
       implicit none
       public :: do_initialize
       contains

       subroutine do_initialize
       use module_model_parameter,only: DP, TS_CHEM_MIN, MAX_NREAC,
     +     MAX_NPROD, MAX_STR1, MAX_NRXN
       use module_domain_common,only:calc_domain_index, NIJK,
     +     grid_1_3_i,grid_1_3_j,grid_1_3_k
       use module_time_common,only:time_step_init
       use module_model_option,only:option_chemmech,option_solver
       use module_geoschem_io,only:geos_read
       use module_ream_io,only:ream_read
       use module_chemmech_common,only:spec_defconc,spec_getid
       use module_chemmech_common,only:nspec,nrxn,nactive,spec_name,
     +                   nreac,reac_id,prod_id,prod_coefs,nphotorxn
       use mod_smvgear_interface,only:smvgear_setup 
       use module_conc_common,only:gas_conc
       use module_met_common,only:update_met,airdensity
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
       call time_step_init(2011,10,10,0,0,0d0,
     +                     2011,10,11,0,0,30d0)

       !Initialize chemical mechanism
       if (option_chemmech.eq.1) then
          filename='globchem.dat'
          call geos_read(trim(filename))
       elseif (option_chemmech.eq.2) then
          filename='chem.dat'
          call ream_read(trim(filename))
       endif

       !Initialize met
       call update_met

       !Chemical initial condition
       ido3=spec_getid('O3')
       idno2=spec_getid('NO2')
       idisop=spec_getid('ISOP')
       do ig = 1, NIJK
          i=grid_1_3_i(ig)
          j=grid_1_3_j(ig)
          k=grid_1_3_k(ig)
          gas_conc(ig,:)  = spec_defconc(:)*airdensity(i,j,k)
          if (ido3.gt.0) gas_conc(ig,ido3) = 40d-9*airdensity(i,j,k)
          if (idno2.gt.0) gas_conc(ig,idno2) = 5d-9*airdensity(i,j,k)
          if (idisop.gt.0) gas_conc(ig,idisop) = 1d-9*airdensity(i,j,k)
       enddo
       !Initialize chemical solver
       if (option_solver.eq.1) then !smvgear
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
