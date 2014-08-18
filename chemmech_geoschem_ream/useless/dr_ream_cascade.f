      program dr_ream_cascade
      use module_ream,only:ream_read
      use module_analyze,only:analyze_do,analyze_cascade
      implicit none
      integer,parameter   :: ne = 14
      character(len=1024) :: pream
      character(len=15)   :: tc
      character(len=15),dimension(ne) :: ec
      character(len=15),dimension(500) :: oc
      integer             :: no,i
      pream="/data12/yzhang425/REAM/REAMSource/v2011_daq_12/Run/"
     +    //"chem.dat"
      tc="ALK4"
      ec(1)="NO"
      ec(2)='NO2'
      ec(3)='OH'
      ec(4)='HO2'
      ec(5)='M'
      ec(6)='O2'
      ec(7)='CO2'
      ec(8)='H2O'
      ec(9)='H2'
      ec(10)='H2O2'
      ec(11)='O3'
      ec(12)='N2O5'
      ec(13)='NO3'
      ec(14)='HNO3'
      call ream_read(trim(pream))
      call analyze_do
      call analyze_cascade(tc,ec,ne,oc,no)

      do i=1,no
         print*,i,oc(i)
      enddo

      endprogram
