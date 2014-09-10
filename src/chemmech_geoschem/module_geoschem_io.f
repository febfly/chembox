      module module_geoschem_io
      use module_geoschem_parameter, only: DP,MAX_NREAC,MAX_NPROD
      implicit none

      public :: geos_read
      
      private :: read_rxn

      contains
      
      subroutine geos_read(filename)
      use module_geoschem_cheminfo,only: spec_add, spec_finish_add, 
     +                                   cheminfo_init,rxn_add
      use module_geoschem_rxntype,only:  MAX_NPARA, rxntype_init

      !common variables used for printing only,..
      use module_geoschem_cheminfo,only: specname, ns, status,def_conc,
     +                                   nr, reacs, prods, r_type, 
     +                                   prod_coefs,nreac_list=>nreac,
     +                                   nprod_list=>nprod
      use module_geoschem_rxntype,only:  symbol

      character(len=*),intent(in)        :: filename
      integer                            :: u,  ifok, tpid, tpid_pre
      character(len=5)                   :: head
      character(len=1)                   :: spec_stat
      character(len=14)                  :: spec_name
      real(kind=DP)                      :: spec_conc
      real(kind=DP),dimension(MAX_NPROD) :: coef
      real(kind=DP),dimension(5)         :: pinp
      integer                            :: nreac, nprod,ord
      integer,dimension(MAX_NREAC)       :: reac_id
      integer,dimension(MAX_NPROD)       :: prod_id
      real(kind=DP),dimension(MAX_NPARA) :: paralist

      integer                            :: i,j,id

      u = 90 !I/O unit
      call cheminfo_init
      call rxntype_init
    
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
         call spec_add(spec_name,spec_stat,pinp(4))
         read(u,10) spec_stat,spec_name,(pinp(j),j=1,5)
         read(u,*)  head
      enddo
      call spec_finish_add

      !print species
      print*,'Species List:'
      do i=1,ns
         print 11,i, specname(i),status(i),def_conc(i)
      enddo
      print*,''

      !read reaction
      head=""
      do while (trim(head).ne.'BEGIN')
         read(u,*) head
      enddo

      !read non-photolysis reactions
      tpid_pre = 0
      call read_rxn(u,tpid_pre,0,reac_id,nreac,prod_id,nprod,
     +              coef,ord, paralist,tpid,ifok)
      do while(ifok.ne.2) !END encountered
         if (ifok.lt.0) then
            stop
         elseif (ifok.eq.0) then
            call rxn_add(0,nreac,nprod,reac_id,prod_id,coef,
     +                  tpid,paralist)
            tpid_pre=tpid
         endif
         call read_rxn(u,tpid_pre,0,reac_id,nreac,prod_id,nprod,
     +              coef,ord, paralist,tpid,ifok)
      enddo

      head=""
      do while (trim(head).ne.'BEGIN')
         read(u,*) head
      enddo

      !read photolysis reactions
      tpid_pre = 0
      call read_rxn(u,tpid_pre,1,reac_id,nreac,prod_id,nprod,coef,ord,
     +              paralist,tpid,ifok)
      do while(ifok.ne.2) 
         if (ifok.lt.0) then
            stop
         elseif (ifok.eq.0) then
            call rxn_add(1,nreac,nprod,reac_id,prod_id,coef,
     +                  tpid,paralist)
         endif
         call read_rxn(u,tpid_pre,1,reac_id,nreac,prod_id,nprod,coef,
     +              ord, paralist,tpid,ifok)
      enddo

      !print reactions
      print*,'Reaction list:'
      do i=1,nr
          print 12,i,nreac_list(i),nprod_list(i),symbol(r_type(i))
      enddo

  10  format(A1,1X,A14,3X,0PF6.2,4(1PE10.3))
  11  format(I4,X,A15,X,A1,E10.2)
  12  format(I3,X,I2,I3,X,A2)
      endsubroutine geos_read


!===============================================================================
! subroutine : read_rxn
!              Read a line in GEOS-Chem format input file
!===============================================================================
      subroutine read_rxn(u,flagid_pre,ifphoto,reacid,nreac,prodid,
     +                   nprod,coef,ord,para,flagid,ifok)
      use module_geoschem_cheminfo,only:spec_getid
      use module_geoschem_rxntype, only:rxntype_id,MAX_NPARA,SYMLEN,
     +                             preceeding_type,succeeding_type
      integer                                :: u,ord,ifphoto
      character(len=14),dimension(MAX_NREAC) :: reac
      integer,dimension(MAX_NREAC)           :: reacid
      real(kind=DP),dimension(MAX_NPROD)     :: coef
      character(len=14),dimension(MAX_NPROD) :: prod
      integer,dimension(MAX_NPROD)           :: prodid
      integer                                :: nreac, nprod
      character(len=1)                       :: stat,stat1
      character(len=SYMLEN)                  :: flag
      real(kind=DP),dimension(MAX_NPARA)     :: para
      integer                                :: flagid,flagid_pre
      character(len=20)                      :: comment
      integer                                :: ifok
      integer                                :: nline

      real(kind=DP)                  :: a,b,e,f,g
      integer                        :: c,d, succ,prec
      character(len=1),dimension(16) :: dummychar
      real,             dimension(4) :: dummyreal
      integer                        :: i,j,id

      !Read parameters
      para(:)=0d0
      read(u,51) stat,ord,a,b,c,d,flag,e,f,g,comment

      para(1)=a
      para(2)=b
      para(3)=float(c)
      nline  =float(d)
      para(4)=e
      para(5)=f
      para(6)=g
      flagid =rxntype_id(flag,ifphoto)

      !Continue to read if this type of special reaction requires more
      !than 1 line of information
      if (nline.ge.1) then 
         do i=1,nline
            read(u,51) stat1,ord,a,b,c,d,flag,e,f,g,comment
            para(i*6+1) = a
            para(i*6+2) = b
            para(i*6+3) = float(c)
            para(i*6+4) = e
            para(i*6+5) = f
            para(i*6+6) = g
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
      reacid = 0
      do j = 1,4
          if (reac(j).ne." ")then
             nreac = nreac + 1
             id = spec_getid(reac(j))
             if (id.eq.0) then
                print*,'cannot find reactants in species list:geos-chem'
                print*,reac(j)
                stop
             endif
             reacid(nreac) = id
         endif
      enddo

      !Read products
       read(u,53)(dummychar(i),coef(i),prod(i),i=1,16)
       nprod=0
       prodid=0
       do j=  1,16
          if (prod(j).ne." ".and.coef(j).ne.0.) then
             nprod = nprod + 1
             id = spec_getid(prod(j))
             if (id.eq.0) then
                print*,'cannot find products in species list:geos'
                stop
             endif
             prodid(nprod) = id
          endif
       enddo

       !Dead reactions
       if (stat.eq.'D') then 
          ifok=1
          return
       endif

       !Check flagid satisfy defined succeeding type of flagid_pre
       if(flagid_pre.gt.0) then
          succ = succeeding_type(flagid_pre)
          if (succ.gt.0) then 
             if (flagid.ne.succ) then
                print*,'Error: read in chem mechanism'
                ifok=-1
                return
             endif
          endif
       endif

       !Check if flagid_pre satisfies defined preceeding type of flagid
       if (flagid.gt.0) then
          prec = preceeding_type(flagid)
          if (prec.gt.0) then
             if (flagid_pre.ne.prec) then
                print*,'Error: read in chem mechanism'
                ifok=-2
                return
             endif
          endif
       endif
       
       ifok=0 !Successfully read in a new reaction

  51   format(A1,1X,I4,1X,ES8.2,1X,ES8.1,1X,I6,1X,I1,1X,A2,F6.2,1X,
     1       2(F6.0,1X),A20)
  52   format(4(A1,0PF5.3,A14))
  53   format(4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14)/
     1        4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14))

      endsubroutine read_rxn

      endmodule module_geoschem_io
