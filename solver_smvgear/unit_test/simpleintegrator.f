C
C File:   simpleintegrator.f
C Author: Yuzhong
C
C Created on July 26, 2014, 11:30 PM
C
      program simpleintegrator
          use omp_lib
          implicit none
          integer, parameter :: n_spec=47
          integer, parameter :: n_active=43
          integer, parameter :: n_rxn=139
          integer, parameter :: n_photo=34
          integer, parameter :: n_mxreac=3
          integer, parameter :: n_mxprod=12
          integer, parameter :: specname_len=14
          integer, parameter :: n_grid=8
          integer :: ncs, ifsun, flag
          character(len=14), dimension(n_spec) :: specname
          real(kind=8), dimension(8,n_spec) :: conc0, conc1
          real(kind=8), dimension(n_spec,8) :: xxconc0, xxconc1
          real(kind=8), dimension(8,n_rxn)  :: rate_const
          real(kind=8), dimension(n_rxn,8)  :: xxrate_const
          integer, dimension(n_mxreac,n_rxn) :: reac
          integer, dimension(n_mxprod,n_rxn) :: prod
          real(kind=8), dimension(n_mxreac,n_rxn) :: reac_coef
          real(kind=8), dimension(n_mxprod,n_rxn) :: prod_coef
          character(len=14) :: rxntype
          character(len=255) :: filename
          integer :: i,j,iloop
          real(kind=8) :: diffre
          integer :: res1,res2
          
          call getdata(n_spec, n_rxn, specname, reac, prod, reac_coef, prod_coef,
     +                 xxconc0, xxconc1, xxrate_const)
          
          ncs = 1
          call smvgear_init

          call smvgear_setpara(ncs,900d0)

!          call smvgear_setblock(ncs) 

          rxntype="TEST"
          call smvgear_setchem(ncs, n_spec, n_active, specname_len, specname, 
     +   n_rxn,  n_photo, n_mxreac, n_mxprod, reac, prod, reac_coef, prod_coef, rxntype)

          call jsparse(ncs)
          !call omp_set_num_treads(2)
          res1=time()
!$OMP PARALLEL DO         
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( iloop, i, ifsun, conc0,rate_const,conc1,flag)
!$OMP+SCHEDULE( STATIC,5)
          do iloop=1,10
              !flag=OMP_get_thread_num()
              !print*,flag,iloop
              ifsun=1
              do i=1,8
                 conc0(i,:)=xxconc0(:,i)
                 rate_const(i,:)=xxrate_const(:,i)
              enddo
             ! print*,'xxx'

             call smvgear_solve (ncs, ifsun, n_grid, n_spec, n_rxn, conc0, rate_const, conc1, flag)
c             write(filename,'(a7,i2.2,a4)')'compar2',i,'.txt'
c             open(unit=19,file=trim(filename))
c             do j=1,n_spec
c                 if (xxconc1(j,i).ne.0.) then
c                     diffre=(conc1(1,j)-xxconc1(j,i))/xxconc1(j,i)
c                 else
c                     diffre=0
c                 endif
c                 write(19,*) specname(j),conc0(1,j),conc1(1,j),xxconc1(j,i), diffre
c             enddo
          enddo
!$OMP END PARALLEL DO
          res2=time()
          print*,res2,res1,res2-res1
      
      endprogram simpleintegrator
