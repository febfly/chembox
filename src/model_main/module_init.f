       module module_init
       implicit none
       public :: do_initialize
       contains

       subroutine do_initialize
       use module_model_parameter,only: DP, TS_CHEM, MAX_NREAC,
     +     MAX_NPROD, MAX_STR1, MAX_NRXN
       use module_domain_common,only:calc_domain_index, NIJK,
     +     grid_1_3_i,grid_1_3_j,grid_1_3_k
       use module_time_common,only:time_init
       use module_model_option,only:option_chemmech,option_solver
       use module_geoschem_io,only:geos_read
       use module_ream_io,only:ream_read
       use module_mcm_io,only:mcm_read

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
       call time_init

       !Initialize chemical mechanism
       if (option_chemmech.eq.1) then
          filename='globchem.dat'
          call geos_read(trim(filename))
       elseif (option_chemmech.eq.2) then
          filename='chem.dat'
          call ream_read(trim(filename))
       elseif (option_chemmech.eq.3) then
          filename='mcm_chem.txt'
          call mcm_read(trim(filename))
       endif

       !Initialize met
       call update_met

       !Chemical initial condition
       ido3=spec_getid('O3')
       do ig = 1, NIJK
          i=grid_1_3_i(ig)
          j=grid_1_3_j(ig)
          k=grid_1_3_k(ig)
          gas_conc(ig,:)  = spec_defconc(:)*airdensity(i,j,k)
          if (ido3.gt.0) gas_conc(ig,ido3) = 40d-9*airdensity(i,j,k)
       enddo

       call read_init_conc

       !Initialize chemical solver
       if (option_solver.eq.1) then !smvgear
          reac_coefs=0d0
          do irxn=1,nrxn
             reac_coefs(1:nreac(irxn),irxn)=1d0
          enddo
          chem_tstep = TS_CHEM !unit: second
          call smvgear_setup(1,chem_tstep,nspec,nactive, MAX_STR1,
     +         spec_name, nrxn, nphotorxn, MAX_NREAC,MAX_NPROD,
     +         reac_id,prod_id,reac_coefs,prod_coefs)
       endif       
       endsubroutine do_initialize

       subroutine read_init_conc
       use module_model_parameter, only: DP, MAX_STR1
       use module_conc_common,     only:gas_conc
       use module_chemmech_common, only:spec_getid
       use module_met_common,      only:airdensity
       use module_domain_common,   only:NIJK,
     +     grid_1_3_i,grid_1_3_j,grid_1_3_k
       character(len=255)      :: text
       character(len=MAX_STR1) :: spec
       real(kind=DP)           :: init_conc
       integer                 :: iostatus, id, ig, i, j, k
       open(unit=12,file='mcm_init.txt')
       read(12,*) text
       iostatus=0
       do while(iostatus.eq.0) 
          read(12,*) spec, init_conc
          id=spec_getid(spec)
          if (id.eq.0) then
             print*,'Warning: cannot find '//trim(spec)//' mcm_init.txt'
          else
             do ig = 1, NIJK
                i = grid_1_3_i(ig)
                j = grid_1_3_j(ig)
                k = grid_1_3_k(ig)
                gas_conc(ig,id) = init_conc*1d-9*airdensity(i,j,k)
             enddo
          endif
       enddo
       close(12)
       endsubroutine
       endmodule module_init
