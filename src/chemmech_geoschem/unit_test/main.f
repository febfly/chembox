      program main
      use module_geoschem_parameter !DP
      use module_geoschem_io,only: geos_read
      use module_geoschem_rxntype,only: rxn_rate
      use module_geoschem_cheminfo,only:nr,r_type,paras
      implicit none
      character(len=255) :: filename
      real(kind=DP) :: Na, R
      real(kind=DP) :: temp, pres, o2,n2, denair, h2o
      real(kind=DP) :: aer_area,aer_radius
      real(kind=DP),allocatable,dimension(:) :: ratek
      integer       :: i

      Na=6.02d23 !Constant
      R =8.314   !Constant

      temp=280!K
      pres=1e5!Pa
      denair=pres*Na/R/temp !#/cm3
      o2=0.21 *denair !#/cm3
      n2=0.78 *denair !#/cm3
      h2o=0.001*denair

      aer_area  =12e-7 !cm2/cm3
      aer_radius=1e-4  !cm

      filename="chem_test1.dat"
      call geos_read(filename)

      allocate(ratek(nr))

      call rxn_rate(nr,r_type,paras,temp,pres,o2,n2,denair,
     +     h2o,aer_area,aer_radius,ratek)

      print*,'Reaction rate constant:'
      do i=1,nr
         print*,i,ratek(i)
      enddo
      endprogram
