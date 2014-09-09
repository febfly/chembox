      module module_geoschem_io
      use module_model_parameter, only: DP
      implicit none

      public : geos_read
      
      private : read_rxn

      contains
      
      subroutine geos_read
      use module_geoschem_cheminfo,only: spec_add, spec_finish_add, rxn_add
      use module_geoschem_rxntype,only:  MAX_NPARA
      character(len=*),intent(in)        :: filename
      integer                            :: u, tpid_pre, ifok, tpid
      character(len=5)                   :: head
      character(len=1)                   :: spec_stat
      character(len=14)                  :: spec_name
      real(kind=DP)                      :: spec_conc
      real(kind=DP),dimension(MAX_NPROD) :: coef
      real(kind=DP),dimension(5)         :: pinp
      integer                            :: nreac, nprod
      integer,dimension(MAX_NREAC)       :: reac_id
      integer,dimension(MAX_NPROD)       :: prod_id
      real(kind=DP),dimension(MAX_NPARA) :: paralist

      integer                            :: i,j,id

      u = 90 !I/O unit
 
      !open input file chem.dat
      open(unit=u,file=trim(filename),status='old')

      head=""
      do while (trim(head).ne.'BEGIN')
         read(u,*) head
      enddo
      read(u,*) head
      read(u,*) head

      !species
      read(u,10) spec_stat,spec_name,(pinp(j),j=1,5)
      read(u,*)  head
      do while(trim(spec_name).ne.'END')
         id=spec_add(spec_name,spec_stat,pinp(4))
         read(u,10) spec_stat,spec_name,(pinp(j),j=1,5)
         read(u,*)  head
      enddo
      call spec_finish_add
      print*,'done reading species'

      !read reaction
      head=""
      do while (trim(head).ne.'BEGIN')
         read(u,*) head
      enddo

      !read non-photolysis reactions
      tpid_pre = 0
      call read_rxn()
      do while(ifok.ne.2) !END encountered
         if (ifok.lt.0) then
            stop
         elseif (ifok.eq.0) then

         endif
         call read_rxn()
      enddo

      head=""
      do while (trim(head).ne.'BEGIN')
         read(u,*) head
      enddo

      !read photolysis reactions
      call read_photorxn()
      do while(ifok.ne.2) 
         if (ifok.lt.0) then
            stop
         elseif (ifok.eq.0) then
            call rxn_add()
         endif
         call read_photorxn()
      enddo
      print*,'done reading reactions'

  10  format(A1,1X,A14,3X,0PF6.2,4(1PE10.3))
      endsubroutine geos_read

      subroutine read_rxn()

      endsubroutine read_rxn

!===============================================================================
! subroutine : read_line
!              Read a line in GEOS-Chem format input file
!===============================================================================
      subroutine read_line(u,ifphoto,reac,reacid,nreac,prod,prodid,nprod,coef,
     +                     stat, ord,para,flag,flagid)
      use module_geoschem_cheminfo,only:spec_getid
      use module_geoschem_rxntype, only:rxntype_id
      integer                        :: u,ord,ifphoto
      character(len=14),dimension(4) :: reac
      integer,dimension(4)           :: reacid
      real(kind=DP),dimension(20)    :: coef
      character(len=14),dimension(20):: prod
      integer,dimension(20)          :: prodid
      integer                        :: nreac, nprod
      character(len=1)               :: stat,flag
      real(kind=DP),dimension(21)    :: para
      integer                        :: flagid
      character(len=20)              :: comment

      real(kind=DP)                  :: a,b,e,f,g
      integer                        :: c,d
      character(len=1),dimension(16) :: dummychar
      real,             dimension(4) :: dummyreal
      integer                        :: i,j,id

      !Read parameters
      para(:)=0d0
      read(u,51) stat,ord,a,b,c,d,tp,e,f,g,comment

      para(1)=a
      para(2)=b
      para(3)=float(c)
      para(4)=float(d)
      para(5)=e
      para(6)=f
      para(7)=g
      flagid =rxntype_id(flag,ifphoto)

      !Continue to read if this type of special reaction requires more
      !than 1 line of information
      if (nrxnline(flagid).gt.1) then 
         do i=2,nrxnline(flagid)
            read(u,51) stat,ord,a,b,c,d,tp,e,f,g,comment
            para((i-1)*7+1) = a
            para((i-1)*7+2) = b
            para((i-1)*7+3) = float(c)
            para((i-1)*7+4) = float(d)
            para((i-1)*7+5) = e
            para((i-1)*7+6) = f
            para((i-1)*7+7) = g
         enddo 
      endif

      !Read reactants
      read(u,52) (dummychar(i),dummyreal(i),reac(i),i=1,4)

      !END encountered
      if (ord.eq.9999) then
         ifok=2
         return
      endif

      nreac=0
      do j = 1,4
          if (reac(j).ne." ")then
             nreac = nreac + 1
             id = spec_getid(reac(j))
             if (id.eq.0) then
                print*,'cannot find reactions in species list:geos-chem'
                print*,reac(j)
                stop
             endif
             reacid(nreac) = id
         endif
      enddo

      !Read products
       read(u,53)(dummychar(i),coef(i),prod(i),i=1,16)
       nprod=0
       do j=  1,16
          if (prod(j).ne." ".and.coef(j).ne.0.) then
             nprod = nprod + 1
             id = spec_getid(prod(j))
             if (id.eq.0) then
                print*,'cannot find product in species list:geos'
                stop
             endif
             prod(nprod) = id
          endif
       enddo

       ifok=0 !Successfully read in a new reaction
  51   format(A1,1X,I4,1X,ES8.2,1X,ES8.1,1X,I6,1X,I1,1X,A1,1X,F6.2,1X,
     1       2(F6.0,1X),A20)
  52   format(4(A1,0PF5.3,A14))
  53   format(4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14)/
     1        4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14))

      endsubroutine read_line

      endmodule module_geoschem_io
