!======================================================================
! module   : module_ream_io
!            provides method to 
!              1. read in the REAM chemical mechanism input file
!              2. write the information in module_ream_cheminfo in the
!                 format of the REAM chemical mechanism input file
! Written by Yuzhong Zhang, 9/8/2014
!======================================================================
      module module_ream_io
      use module_ream_parameter,only: DP,MAX_NREAC,MAX_NPROD
      implicit none
      !public functions
      public  :: ream_read
      !public  :: ream_write

      !private functions
      private :: read_line
      private :: read_rxn
      private :: read_photorxn
      !private :: write_line
      !private :: write_rxn
      !private :: write_photo

      contains

!======================================================================
! subroutine : ream_read
!              reads a REAM chemical mechanism input file
!              stores relevant information in module_ream_cheminfo
!======================================================================
      subroutine ream_read(filename)
      use module_ream_cheminfo,only : cheminfo_init,spec_add,
     +                               spec_finish_add,rxn_add
      use module_ream_rxntype,only  : MAX_NPARA,rxntype_init
      character(len=*),intent(in)        :: filename
      integer                            :: u, tpid_pre, ifok, tpid
      character(len=5)                   :: head
      character(len=1)                   :: spec_stat
      character(len=8)                   :: spec_name
      real(kind=DP)                      :: spec_conc
      real(kind=DP),dimension(MAX_NPROD) :: coef
      integer                            :: nreac, nprod
      integer,dimension(MAX_NREAC)       :: reac_id
      integer,dimension(MAX_NPROD)       :: prod_id
      real(kind=DP),dimension(MAX_NPARA) :: paralist

      u = 90 !I/O unit
      call cheminfo_init
      call rxntype_init

      !open input file chem.dat
      open(unit=u,file=trim(filename),status='old')

      !read heading for species list
      head = ''
      do while(head.ne.'BEGIN')
         read(u,'(a5)') head
      enddo

      !read species list
      spec_name=''
      read(u,11) spec_stat, spec_name,spec_conc
      do while(spec_name.ne.'END')
         call spec_add(spec_name, spec_stat, spec_conc)
         read(u,11) spec_stat, spec_name,spec_conc
      enddo
      call spec_finish_add
      print*,'done read species'

      !read heading for reaction list
      head = ''
      do while(head.ne.'BEGIN')
         read(u,'(a5)') head
      enddo

      !read non-photolysis reactions
      tpid_pre = 0
      call read_rxn(u, tpid_pre, nreac, nprod, reac_id,prod_id,tpid,
     +           coef,paralist, ifok)

      do while (ifok.ne.2) !END encountered
         if (ifok.lt.0) then
            !If error occurs
            stop
         elseif (ifok.eq.0) then 
            !If an active reaction, add
            call rxn_add(0,nreac,nprod, reac_id, prod_id,coef,
     +           tpid, paralist)
            tpid_pre = tpid
         endif
         call read_rxn(u, tpid_pre, nreac, nprod, reac_id,prod_id,tpid,
     +           coef,paralist, ifok)
      enddo
      print*,'done read non-photolysis reactions'

      !read photolysis reactions
      head=""
      do while (trim(head).ne.'BEGIN')
         read(u,*) head
      enddo

      call read_photorxn(u,nreac, nprod, reac_id,prod_id,tpid,
     +           coef,paralist, ifok)
      do while (ifok.ne.2) !END encountered
         if (ifok.lt.0) then 
            !if error occurs
            stop
         elseif (ifok.eq.0) then
             call rxn_add(1,nreac,nprod,reac_id,prod_id,coef,
     +                   tpid,paralist)
         endif
         call read_photorxn(u,nreac, nprod, reac_id,prod_id,tpid,
     +           coef,paralist, ifok)
      enddo
      print*,'done read photolysis reactions'

 11   format(a1, 1x, a8, 1x, e10.3)
      endsubroutine ream_read

!==============================================================
! subroutine : read_rxn
!              Read a reaction and associated parameters in chem.dat
!==============================================================
      subroutine read_rxn(u, tpid_pre, nreac0, nprod0, reac_id0, prod_id0,
     +           tpid0, coef0, paralist, ifok)
      use module_ream_rxntype,only:nrxnline, preceeding_type,
     +           succeeding_type, MAX_NPARA

      integer, intent(in)  :: u, tpid_pre
      integer              :: nreac0, nreac,nprod0, nprod,tpid0, tpid
      integer                               :: iord
      character(len=1)                      :: dinp, flag,flagtmp
      character(len=8),dimension(MAX_NREAC) :: reac_name,reac_name2
      character(len=8),dimension(MAX_NPROD) :: prod_name,prod_name2
      real(kind=DP),dimension(MAX_NPROD)    :: coef0, coef, coef2
      integer,dimension(MAX_NREAC)          :: reac_id, reac_id0
      integer,dimension(MAX_NPROD)          :: prod_id, prod_id0
      
      integer, intent(out):: ifok
      real(kind=DP),dimension(3)         :: para
      real(kind=DP),dimension(MAX_NPARA),intent(out) :: paralist
      integer :: succ, prec, il

      paralist(:) = 0.

      !Read a line
      call read_line(u,reac_name,reac_id0,nreac0,prod_name,prod_id0,nprod0,
     +               coef0, dinp,iord,para,flag,tpid0)
      !if it is the last line of reaction list
      if (reac_name(1).eq.'END') then
         ifok=2
         return
      endif
      !if it is a dead reaction
      if (dinp.eq.'D') then
         ifok=1
         return
      endif
         
      paralist(1:3) = para(1:3)

      !Check if tpid and tpid_pre satisfy defined preceeding_type
      !and/or succeeding_type
      if (tpid_pre.gt.0) then
         succ = succeeding_type(tpid_pre)
         if (succ.gt.0) then
            if (tpid0.ne.succ) then
               print*,'Error: chem.dat'
               print*,'***Requires certain type of special reaction***'
               print*,'***Check module_ream_rxntype.f'
               ifok=-1 
               return
            endif
         endif
      endif

      !Check if tpid and tpid_pre satisfy defined succeeding_type
      !and/or succeeding_type
      if (tpid0.gt.0) then
         prec = preceeding_type(tpid0)
         if (prec.gt.0) then
            if (tpid_pre.ne.prec) then
               print*,'Error: chem.dat'
               print*,'***Requires certain type of special reaction***'
               print*,'***Check module_ream_rxntype.f'
               ifok=-2
              return
           endif
         endif
      endif

      !Continue to read if this type of special reaction requires more
      !than 1 line of information
      if (nrxnline(tpid0).gt.1) then
         !Used for checking if the next line is the same reaction or not
         reac_name2=reac_name
         prod_name2=prod_name
         coef2     =coef0
         !Read multiple lines
         do il = 2, nrxnline(tpid)
            para(:) = 0.
            call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,
     +           nprod, coef, dinp, iord, para, flagtmp, tpid, 
     +           reac_name2, prod_name2,coef2)
            paralist((il-1)*3+1:il*3) = para(1:3)
         enddo
      endif
      ifok = 0 !successful reading a reaction
      return
      endsubroutine read_rxn

!===============================================================================
! subroutine : read_line
!              Read a line in chem.dat
!===============================================================================
      subroutine read_line(u,reac,reacid,nreac,prod,prodid,nprod,coef,
     &                     stat,ord,para,flag,flagid,reac0,prod0,coef0)
      use module_ream_cheminfo,only:spec_getid
      use module_ream_rxntype ,only:rxntype_id
      integer                       :: u,ord
      character(len=8),dimension(4) :: reac
      character(len=8),dimension(4),optional :: reac0
      integer,dimension(4)          :: reacid
      real(kind=DP),dimension(20)            :: coef
      real(kind=DP),dimension(20),optional   :: coef0
      character(len=8),dimension(20):: prod
      character(len=8),dimension(20),optional:: prod0
      integer,dimension(20)         :: prodid
      integer                       :: nreac, nprod
      character(len=1)              :: stat,flag
      real(kind=DP),dimension(3)             :: para
      integer                       :: flagid

      real(kind=DP)                 :: a,b,c
      integer                       :: c0  , i, id
      character(len=10)             :: comment
      !read reactants
      read(u,12) (reac(i), i=1,4)
      !read parameters
      read(u,13) ord,stat,nprod,a,b,c0,flag, comment
      !read products
      do i = 1, nprod
          read(u,14) coef(i), prod(i)
      enddo

      !reaction parameters 
      c = float(c0)
      para(1) = a
      para(2) = b
      para(3) = c
      
      !reaction type
      flagid = rxntype_id(flag,0)
      if (flagid.eq.0 ) then
         print*,'cannot find reaction type '//flag
         stop
      endif

      nreac=0
      reacid = 0
      !count number of reactants and find their species id
      do i=1,4
         if (reac(i).ne." ")then
            nreac = nreac + 1
            id = spec_getid(reac(i))
            if (id.eq.0) then
               if (reac(i).ne.'END') then
                  print*,'cannot find reactants in species list'
                  stop
               else
                  return
               endif
            endif
            reacid(nreac) = id
         endif
      enddo

      !find product's species id
      prodid = 0
      do i=1,nprod
         if (prod(i).ne." ")then
            id = spec_getid(prod(i))
            if (id.eq.0) then
               print*,'cannot find product in species list'
               stop
            endif
            prodid(i) = id
         endif
      enddo

      !Additional checking for special reactions
      !check if reactant name matches
      if (present(reac0)) then
         do i = 1,nreac-1 ! 
            if (reac(i).ne.reac0(i)) then
               print*,'reactants not match last reaction, check!'
               stop
            endif
         enddo
      endif

      !check if products name and coef match
      if (present(prod0).and.present(coef0)) then
         do i = 1, nprod
            if (prod(i).ne.prod0(i).or.coef(i).ne.coef0(i)) then
               print*,'products not match last reaction, check!'
               stop
            endif
         enddo
      endif

 12   format(1x,4(a8,1x))
 13   format(i4,1x,a1,1x,i2,1x,es9.2,0p,1x,es9.2,0p,1x,i6,1x,a1,1x,a9)
 14   format(2x, f5.3, 1x, a8)
      endsubroutine read_line

!==============================================================
! subroutine : read_photoline
!              read a photolysis reaction
!==============================================================
      subroutine read_photorxn(u,nreac, nprod, reacid,prodid,
     &           flagid,coef,para,ifok)
      use module_ream_cheminfo,only:spec_getid
      use module_ream_rxntype ,only:rxntype_id
      integer                       :: u,ord
      character(len=8),dimension(4) :: reac
      integer,dimension(4)          :: reacid
      real(kind=DP),dimension(20)   :: coef
      character(len=8),dimension(20):: prod
      integer,dimension(20)         :: prodid
      integer                       :: nreac, nprod
      character(len=1)              :: stat,flag
      real(kind=DP),dimension(3)    :: para
      integer                       :: flagid, ifok

      integer                       :: i,id
      integer                       :: ifphoto

      !Read the line
      read(u,15)reac(1),ord,stat,nprod,para(1),flag
      !END encountered
      if (reac(1).eq.'END') then
         ifok=2 
         return
      endif

      do i = 1, nprod
         read(u,16) coef(i),prod(i)
      enddo
      !Dead reaction encountered
      if (stat.ne.'A') then
         ifok=1
         return
      endif

      !Check reaction type
      ifphoto=1
      flagid = rxntype_id(flag,ifphoto)
      if (flagid.eq.0 ) then
         print*,'cannot find reaction type '//flag
         stop
      endif

      !Find reactant
      nreac = 1
      reacid = 0
      id = spec_getid(reac(1))
      if (reac(1).eq.'END') then 
         ifok=2
         return
      endif
      if (id.eq.0) then
          print*,'cannot find reactants'//reac(1)// 'in species list'
          stop
      endif

      !Find products
      do i=1,nprod
         if (prod(i).ne." ")then
            id = spec_getid(prod(i))
            if (id.eq.0) then
               print*,'cannot find product in species list'
               stop
            endif
            prodid(i) = id
         endif
      enddo
      ifok = 0

 15   format(1x,a8,1x,i4,1x,a1,1x,i2,1x,f4.1,1x,a1)
 16   format(2x,f5.3,1x,a8)

      endsubroutine read_photorxn

      endmodule module_ream_io
