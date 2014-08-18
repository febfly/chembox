      module module_ream
      implicit none
      contains
!========================================================================
      subroutine ream_read(filename)
      use module_spec, only : spec_init, spec_add, spec_finish_add,ns
      use module_reaction, only: rxn_init,rxn_add,reaction,nr
      use module_photolysis,only:photolysis_index

      character(len=*),intent(in) :: filename

      character(len=5)            :: head
      character(len=1)            :: spec_stat
      character(len=8)            :: spec_name
      real                        :: spec_conc

      character(len=8),dimension(4) :: reac_name, reac_name2
      integer,dimension(4)          :: reac_id

      integer                       :: iord
      character(len=1)              :: dinp, flag, flag0, flagtmp

      real,dimension(20)            :: coefs, coefs2
      character(len=8),dimension(20):: prod_name, prod_name2
      integer,dimension(20)         :: prod_id
      real,dimension(3)             :: para
      type(reaction)                :: rtemp

      integer                     :: u, i, id, nreac, nprod, typeid
      u = 90

      call spec_init
      call rxn_init

      open(unit=u,file=trim(filename),status='old')

      !read heading
      head = ''
      do while(head.ne.'BEGIN')
         read(u,'a5') head
      enddo

      !read spec_list
      spec_name=''
      read(u,11) spec_stat, spec_name,spec_conc
      do while(spec_name.ne.'END')
         id=spec_add(spec_name, spec_stat, spec_conc)
         read(u,11) spec_stat, spec_name,spec_conc
      enddo
      call spec_finish_add
      print*,'done read species'
      !read heading
      head = ''
      do while(head.ne.'BEGIN')
         read(u,'a5') head
      enddo

      !read reaction
      do while(reac_name(1).ne.'END')

         call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &                     dinp,iord,para,flag,typeid)
         if (reac_name(1).eq.'END') exit
         if (dinp.eq.'D') cycle
         rtemp%nreac = nreac
         rtemp%reacs(1:nreac) = reac_id(1:nreac)
         rtemp%nprod = nprod
         rtemp%prods(1:nprod) = prod_id(1:nprod)
         rtemp%prod_coefs(1:nprod) = coefs(1:nprod)
         rtemp%status         = dinp
         rtemp%sn             = iord
         rtemp%r_type         = typeid
         rtemp%paras(1:3)     = para(1:3)
         rtemp%rindex         = 0
c         flag0
         !special types
         reac_name2 = reac_name
         coefs2 = coefs
         prod_name2 = prod_name
         if (flag.eq.'E') then !equilibrium
             if (flag0.ne.'P') stop
         elseif (flag .eq. "P") then !pressure dependent
             if (reac_name(2).ne.'M'.and.reac_name(3).ne.'M') then
                print*,'Pressure dependent reaction,',reac_name
                stop
             endif

             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(4:6) = para(1:3)
            
             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(7)  = para(1) 

         elseif (flag.eq.'V') then !MCO3+MO2
             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(4:6) = para(1:3) 
         elseif (flag.eq.'B') then
             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(4:6) = para(1:3) 
         elseif (flag.eq.'A') then
             if (flag0.ne.'B') stop
             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(4) = para(1) 
         elseif (flag.eq.'X') then
             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(4:6) = para(1:3) 
             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(7:9) = para(1:3) 
         elseif (flag.eq.'Z') then
             call read_line(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flagtmp,typeid,reac_name2, prod_name2,coefs2)
             rtemp%paras(4:6) = para(1:3)             
         endif
         flag0=flag
         call rxn_add(rtemp,.false.)
      enddo

      print*,'done read reactions'
      !read photolysis reactions
      head=""
      do while (trim(head).ne.'BEGIN')
         read(u,*) head
      enddo

      do while (reac_name(1).ne.'END')
         call read_photoline(u,reac_name,reac_id,nreac,prod_name,prod_id,nprod,coefs,
     &            dinp,iord,para,flag,typeid)
         if (reac_name(1).eq.'END') exit
         if (dinp.ne.'A') cycle
         rtemp%nreac = nreac
         rtemp%reacs(1:nreac) = reac_id(1:nreac)
         rtemp%prods(1:nprod)= prod_id(1:nprod)
         rtemp%prod_coefs(1:nprod) = coefs(1:nprod)
         rtemp%status = dinp
         rtemp%sn     = iord
         rtemp%r_type = typeid
         rtemp%paras(1) = para(1)
         rtemp%rindex = photolysis_index(reac_name(1),dble(para(1)))
         call rxn_add(rtemp,.true.)
      enddo
      print*,'done read photolysis'
      close(u)
      print*,'number of species',ns
      print*,'number of reactions',nr
 11   format(a1, 1x, a8, 1x, e10.3)
      endsubroutine ream_read
!===================================================================================
      subroutine read_line(u,reac,reacid,nreac,prod,prodid,nprod,coef,
     &                     stat,ord,para,flag,flagid,reac0,prod0,coef0,
     &                     _ifphoto)
      use module_spec,only:spec_getid
      use module_reaction,only:rxn_findtype
      integer                       :: u,ord
      logical,optional              :: _ifphoto
      logical                       :: ifphoto
      character(len=8),dimension(4) :: reac 
      character(len=8),dimension(4),optional :: reac0
      integer,dimension(4)          :: reacid
      real,dimension(20)            :: coef
      real,dimension(20),optional   :: coef0
      character(len=8),dimension(20):: prod
      character(len=8),dimension(20),optional:: prod0
      integer,dimension(20)         :: prodid
      integer                       :: nreac, nprod
      character(len=1)              :: stat,flag
      real,dimension(3)             :: para
      integer                       :: flagid
     
      real                          :: a,b,c
      integer                       :: c0  , i, id
      character(len=10)             :: comment
      read(u,12) (reac(i), i=1,4)      
      read(u,13) ord,stat,nprod,a,b,c0,flag, comment
      do i = 1, nprod
          read(u,14) coef(i), prod(i)
      enddo

      c = float(c0)  
      para(1) = a
      para(2) = b
      para(3) = c
      ifphoto=.false.
      if (present(_ifphoto)) ifphoto=_ifphoto
      flagid = rxn_findtype(flag,1,ifphoto)
      if (flagid.eq.0 ) then 
         print*,'cannot find reaction type '//flag
         stop
      endif

      nreac=0
      reacid = 0
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

      if (present(reac0)) then
         do i = 1,nreac-1 
            if (reac(i).ne.reac0(i)) then 
               print*,'reactants not match last reaction, check!'
               stop
            endif
         enddo
      endif

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
      endsubroutine

!==============================================================
      subroutine read_photoline(u,reac,reacid,nreac,prod,prodid,nprod,
     &           coef,stat,ord,para,flag,flagid)
      use module_spec,only:spec_getid
      use module_reaction,only:rxn_findtype
      integer                       :: u,ord
      character(len=8),dimension(4) :: reac 
      integer,dimension(4)          :: reacid
      real,dimension(20)            :: coef
      character(len=8),dimension(20):: prod
      integer,dimension(20)         :: prodid
      integer                       :: nreac, nprod
      character(len=1)              :: stat,flag
      real,dimension(3)             :: para
      integer                       :: flagid
      
      integer                       :: i,id
      logical                       :: ifphoto
      read(u,15)reac(1),ord,stat,nprod,para(1),flag
      do i = 1, nprod
         read(u,16) coef(i),prod(i)
      enddo

      ifphoto=.true.
      flagid = rxn_findtype(flag,1,ifphoto)
      if (flagid.eq.0 ) then 
         print*,'cannot find reaction type '//flag
         stop
      endif

      nreac = 1
      reacid = 0
      id = spec_getid(reac(1))
      if (reac(1).eq.'END') return
      if (id.eq.0) then 
          print*,'cannot find reactants'//reac(1)// 'in species list'
          stop
      endif

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


 15   format(1x,a8,1x,i4,1x,a1,1x,i2,1x,f4.1,1x,a1)
 16   format(2x,f5.3,1x,a8)
 
      endsubroutine

!===============================================================
      subroutine ream_write
      use module_spec,only: ns, spec
      use module_reaction,only:nr,nphoto,rxn
      integer           :: u
      integer           :: i, j

      u = 98
      open(unit=u,file='chem.dat')
      write(u,*) ' *** chem.dat'
      write(u,*) ' *** created by zyz'
      write(u,*) ' *** species list ***'
      write(u,*) 'BEGIN'
      do i = 1, ns
         write(u,'(a1,x,a8,x,e10.3)')
     $   spec(i)%status, spec(i)%name, spec(i)%default_conc
      enddo
      write(u,*) '  END'
      write(u,*) ''
      write(u,*) ' *** reaction list ***'
      write(u,*) 'BEGIN'
      !print*,rxn(2)%reacs
      do i = 1,nr-nphoto
         call write_rxn(u,rxn(i)%nreac,rxn(i)%reacs,rxn(i)%nprod,
     $    rxn(i)%prods,rxn(i)%prod_coefs,rxn(i)%status,
     $    rxn(i)%sn, rxn(i)%r_type, rxn(i)%paras)
      enddo
      write(u,*) ' END'
      write(u,*) ''
      write(u,*) ''
      write(u,*) '*** photolysis reaction list ***'
      write(u,*) 'BEGIN'
      do j=1,nphoto
         i=j+nr-nphoto
         call write_photo(u,rxn(i)%reacs(1),rxn(i)%nprod,
     $    rxn(i)%prods,rxn(i)%prod_coefs,rxn(i)%sn,rxn(i)%status,
     $    rxn(i)%paras(1),rxn(i)%r_type)
      enddo
      write(u,*) ' END         0'
      close(u)
      endsubroutine ream_write
!=========================================================================
      subroutine write_rxn(u,nreact,react,nprod,prod,coef,stat,ord,flag,para)
      use module_reaction,only:reac_max,prod_max,para_max,rxn_type
      use module_spec,only:spec_getid
      integer :: u,nprod, ord,nreact
      integer        ,dimension(reac_max) :: react
      integer        ,dimension(prod_max):: prod
      real*8          ,dimension(prod_max):: coef
      real*8          ,dimension(para_max)  :: para
      integer                       :: flag
      character(len=1)              :: stat,sym
!      character(len=9)              :: comment
!      comment=" "
      sym=rxn_type(flag)%sym1
      call write_line(u,nreact,react,nprod,prod,coef,ord,stat,
     &         para(1),para(2),para(3),flag)

      if (sym.eq.'P') then
         call write_line(u,2,react,nprod,prod,coef,ord,stat,
     &          para(4),para(5),para(6),flag)
         call write_line(u,2,react,nprod,prod,coef,ord,stat,
     &         para(7),0d0,0d0,0)
      endif

      if (sym.eq.'B'.or.sym.eq.'A'.or.sym.eq.'V'.or.sym.eq.'Z') then
         call write_line(u,nreact,react,nprod,prod,coef,ord,stat,
     &         para(4),para(5),para(6),0)
      endif

      if (sym.eq.'X') then
         call write_line(u,nreact,react,nprod,prod,coef,ord,stat,
     &         para(4),para(5),para(6),flag)
         react(3)=spec_getid('M')
         call write_line(u,3,react,nprod,prod,coef,ord,stat,
     &         para(7),para(8),para(9),0)
      endif

      endsubroutine write_rxn
!=========================================================================
      subroutine write_line(u,nreact,react,nprod,prod,coef,ord,stat,
     &           a,b,c, flag)
      use module_spec,only: spec
      use module_reaction,only:reac_max,prod_max,para_max,rxn_type
      integer  :: u, nprod, nreact,ord
      integer        ,dimension(reac_max) :: react
      integer        ,dimension(prod_max):: prod
      real*8          ,dimension(prod_max):: coef
      character(len=1)              :: stat,sym
      integer                       :: flag
      character(len=9)              :: comment
      integer                       :: i
      real*8                        :: a,b,c
      comment=" "
      if (flag.eq.0) then
        sym=" "
      else
        sym=rxn_type(flag)%sym1
      endif
      
      if (nreact.eq.3) then
      write(u,12) spec(react(1:nreact))%name
      write(6,12) spec(react(1:nreact))%name
      elseif (nreact.eq.2) then
      write(u,11) spec(react(1:nreact))%name
      write(6,11) spec(react(1:nreact))%name
      elseif (nreact.eq.1) then
      write(u,10) spec(react(1))%name
      write(6,10) spec(react(1))%name
      endif

      write(u,13) ord,stat,nprod,a,b,floor(c+0.1d0),sym,comment
      do i=1,nprod
         write(u,14) coef(i),spec(prod(i))%name
      enddo

 12   format(1x,3(a8,1x))
 11   format(1x,2(a8,1x))
 10   format(1x,a8,1x)
 13   format(i4,1x,a1,1x,i2,1x,es9.2,0p,1x,es9.2,0p,1x,i6,1x,a1,1x,a9)
 14   format(2x, f5.3, 1x, a8)
      endsubroutine write_line
!=====================================================================
      subroutine write_photo(u,react,nprod,prod,coef,ord,stat,
     &         para,flag)
      use module_spec,only: spec
      use module_reaction,only:reac_max,prod_max,para_max,rxn_type
      integer  :: u, nprod,ord
      integer  :: react
      integer, dimension(prod_max) :: prod
      real*8 , dimension(prod_max) :: coef
      character(len=1)             :: stat
      integer                      :: flag
      real*8                       :: para
      character(len=9)             :: comment
      integer                      :: i
      real*8                       :: defj

      comment=" "
      write(u,15) spec(react)%name,ord,stat,nprod,para,
     &                rxn_type(flag)%sym1,comment
      do i = 1, nprod
         write(u,14) coef(i),spec(prod(i))%name
      enddo
 14   format(2x, f5.3, 1x, a8)
 15   format(1X,A8,1X,I4,1X,A1,1X,I2,1X,F4.1,1X,A1,1X,A9)
      endsubroutine write_photo

      endmodule module_ream
