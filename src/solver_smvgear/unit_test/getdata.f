C
C File:   getdata.f
C Author: Yuzhong
C
C Created on July 27, 2014, 1:32 PM
C
      subroutine getdata(n_spec, n_rxn, specname, reac, prod, reac_coef, prod_coef,
     +                 conc0, conc1, rate_const)
          use mod_comode
          implicit none
          integer :: n_spec, n_rxn
          character(len=14), dimension(n_spec) :: specname
          integer, dimension(3, n_rxn) :: reac
          integer, dimension(12, n_rxn) :: prod
          real(kind=8), dimension(3, n_rxn) :: reac_coef
          real(kind=8), dimension(12, n_rxn) :: prod_coef
          real(kind=8), dimension(n_spec,8) :: conc0, conc1
          real(kind=8), dimension(n_rxn,8) :: rate_const
                  
          integer :: xns, xnts, xnr, xntr, xnar
          character(len=14), dimension(0:139) :: xname
          integer, dimension(  16,352) :: xirm
          real,dimension(16, 352) :: xcoef
          integer, dimension(352) ::xrxnmap
          real,dimension(139,8)   :: xc0,xc1
          real,dimension(352,8)   :: xrt, xrt0
          integer,dimension(NMTRATE*2) :: xnewfold
          integer a,i, j, xk, xs, st,ed
          character(len=14) :: blk
          character(len=6*8) :: str
          character(len=20) :: filename
          blk="              "
          open(unit=19,file="data/CHEMDAT",form='unformatted')
          read(19) xnts
          print*,xnts
          read(19) xns
          read(19) xntr
          read(19) xnr
          read(19) xnar
          read(19) xname
          read(19) xirm
          read(19) xcoef
          read(19) xrxnmap
          close(19)
          
          print*,xnts,xns,n_spec
          print*,xntr,xnr,xnar,n_rxn
          open(unit=19,file="data/CHEM-NEWFOLD",form='unformatted')
          read(19) xnewfold
          close(19)
          
          specname(:)=xname(1:n_spec)
          do i=1, xnar                 
               xk=xrxnmap(i)
               reac(:,i)=xirm(1:3,xk)
               reac_coef(:,i)=xcoef(1:3,xk)
               prod(:,i)=xirm(5:16,xk)
               prod_coef(:,i)=xcoef(5:16,xk)
          enddo!i
          
          
          do i=1,8 
              write(filename,'(a2,i2.2)') 'CR',i
              open(unit=19,file="data/"//trim(filename),form='unformatted')
              read(19) xc0(:,i)
              read(19) xrt0(:,i)
              read(19) xrt(:,i)
              read(19) xc1(:,i)           
              close(19)
          enddo
          
          do a=1,8
              do i=1, xnar                 
                 xk=xrxnmap(i)
                 rate_const(i,a)=xrt(xnewfold(xk),a)
              enddo!i
          enddo!a
          conc0(:,:) =xc0(1:n_spec,:)
          conc1(:,:) =xc1(1:n_spec,:)
          
      endsubroutine getdata
