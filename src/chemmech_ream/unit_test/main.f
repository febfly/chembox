      program main
      use moduel_chemmech_common
      use module_ream_common !DP,MAX_NSPEC
      use module_ream_io,only: ream_read
      use module_ream_rxntype,only: rxn_rate
      use module_ream_cheminfo,only:nr,r_type,paras,nphoto
      use module_ream_cheminfo,only:ninactrxn,inactrxn,
     +                  nemisrxn,emisrxn,ndeprxn,deprxn
      use module_ream_cheminfo,only:ns,specname,spec_getid,nactive
      implicit none
      character(len=255) :: filename
      real(kind=DP) :: Na, R
      real(kind=DP) :: temp, pres, o2,n2, denair, h2o
      real(kind=DP) :: aer_area,aer_radius
      real(kind=DP),allocatable,dimension(:) :: ratek
      real(kind=DP),dimension(MAX_NSPEC)     :: conc
      integer       :: i

      Na=6.02d23 !Constant
      R =8.314   !Constant

      temp=280!K
      pres=1e5!Pa
      denair=pres*Na/R/temp/1e6 !#/cm3
      denair=2.5d19
      o2=0.21 *denair !#/cm3
      n2=0.78 *denair !#/cm3
      h2o=0.001*denair

      aer_area  =12e-7 !cm2/cm3
      aer_radius=1e-4  !cm

      filename="chem_test1.dat"


      call ream_read(filename)

      allocate(ratek(nr))
      conc(:)=0d0
      conc(nactive+1:ns)    = 0.22060057
      conc(spec_getid('H2'))=1.1d13
      conc(spec_getid('CH4'))=4.3d13
      conc(spec_getid('MOH'))=4.4d10
      conc(spec_getid('M'))  =2.5d19
      conc(spec_getid('O2')) =0.21*2.5d19
      conc(spec_getid('N2')) =0.78*2.5d19

      call rxn_rate(temp,pres,o2,n2,denair,
     +     h2o,aer_area,aer_radius,ratek)

      do i=1,nemisrxn
         ratek(emisrxn(i))=1d0
      enddo
       print*,'DEP'
      do i=1,ndeprxn
         print*,i,deprxn(i)
         ratek(deprxn(i))=2d0
      enddo
      do i=1,ninactrxn
         print*,specname(inactrxn(2,i))
         ratek(inactrxn(1,i))=ratek(inactrxn(1,i))*conc(inactrxn(2,i))
      enddo
      print*,'Reaction rate constant:'
      do i=1,nr
         print*,i,ratek(i)
      enddo
      endprogram
