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
      do i = 1,nr-nphoto
         call write_rxn(u,rxn(i)%reacs,rxn(i)%nprod,
     $    rxn(i)%prods,rxn(i)%prod_coefs,rxn(i)%status,
     $    rxn(i)%sn, rxn_type(rxn(i)%r_type)%sym1, rxn(i)%paras)
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
     $    rxn(i)%paras(1),rxn_type(rxn(i)%r_type)%sym1)
      enddo
      write(u,*) ' END         0'
      close(u)
      endsubroutine ream_write      

!=========================================================================
      subroutine write_rxn(u,react,nprod,prdct,coef,stat,ord,flag,para)
      use module_reaction,only:reac_max,prod_max,para_max

      integer :: u,nprod, ord
      integer        ,dimension(reac_max) :: react
      integer        ,dimension(prod_max):: prdct
      real*8          ,dimension(prod_max):: coef
      real*8          ,dimension(para_max)  :: para
      character(len=1)              :: stat,flag
      character(len=9)              :: comment
      comment=" "
      call write_line(u,4,react,nprod,prod,coef,ord,stat,
     &         para(1),para(2),para(3),flag,comment)
      
      if (flag.eq.'P') then
         call write_line(u,2,react,nprod,prod,coef,ord,stat,
               para(4),para(5),para(6),flag,comment)
         call write_line(u,2,react,nprod,prod,coef,ord,stat,
     &         para(7),0d0,0d0,' ',comment) 
      endif

      if (flag.eq.'B'.or.flag.eq.'A'.or.flag.eq.'V'.or.flag.eq.'Z') then
         call write_line(u,4,react,nprod,prod,coef,ord,stat,
     &         para(4),para(5),para(6),' ',comment)
      endif

      if (flag.eq.'X') then
         call write_line(u,2,react,nprod,prod,coef,ord,stat,
               para(4),para(5),para(6),flag,comment)
         react(3)='M'
         call write_line(u,4,react,nprod,prod,coef,ord,stat,
     &         para(7),para(8),para(9),' ',comment)
      endif

      endsubroutine write_rxn

!=========================================================================
      subroutine write_line(u,nreact,react,nprod,prod,coef,ord,stat,
     &           a,b,c, flag)
      use module_spec,only: spec
      use module_reaction,only:reac_max,prod_max,para_max
      integer  :: u, nprod, nreact,ord
      integer        ,dimension(reac_max) :: react
      integer        ,dimension(prod_max):: prod
      real*8          ,dimension(prod_max):: coef
      character(len=1)              :: stat,flag
      character(len=9)              :: comment
      integer                       :: i
      real*8                        :: a,b,c
      comment=" "
      if (nreact.eq.2) then 
         write(u,11) spec(react(1:2))%name 
      else 
         write(u,12) spec(react)%name
      endif
      write(u,13) ord,stat,nprod,a,b,c,flag,comment
      do i=1,nprod
         write(u,14) coef(i),spec(prod(i))%name
      enddo

 12   format(1x,3(a8,1x))
 11   format(1x,2(a8,1x))
 13   format(i4,1x,a1,1x,i2,1x,es9.2,0p,1x,es9.2,0p,1x,i6,1x,a1,1x,a9)
 14   format(2x, f5.3, 1x, a8)
      endsubroutine write_line

!=====================================================================
      subroutine write_photo(u,react,nprod,prod,coef,ord,stat,
     &         para,flag)
      use module_spec,only: spec
      use module_reaction,only:reac_max,prod_max,para_max
      integer  :: u, nprod
      integer  :: react
      integer, dimension(prod_max) :: prod
      real*8 , dimension(prod_max) :: coef
      character(len=1)             :: stat, flag
      real*8                       :: para
      character(len=9)             :: comment
      integer                      :: i
      real*8                       :: defj
      
      comment=" "
      write(u,15) spec(react)%name,ord,stat,nprod,para,flag,comment
      do i = 1, nprod
         write(u,14) coef(i),spec(prod(i))%name
      enddo
 14   format(2x, f5.3, 1x, a8)
 15   format(1X,A8,1X,I4,1X,A1,1X,I2,1X,F4.1,1X,A1,1X,A9)
      endsubroutine write_photo
