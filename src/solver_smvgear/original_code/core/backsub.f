      subroutine backsub(ncsp)
!
      use mod_comode
      implicit none
      !input variables
      integer, intent(in) :: ncsp
      
      !local variables
      integer :: ij,i,kzt,kl5,kh5,kl4,kh4,kl3,kh3,kl2,kh2,kl1,kh1,kc
      integer :: j0,ij0,ij1,ij2,ij3,ij4,j1,j2,j3,j4,k,mzt,ml5,mh5,ml4,mh4
      integer :: ml3,mh3,ml2,mh2,ml1,mh1,mc
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!     BBBBBBB     A     CCCCCCC  K     K  SSSSSSS  U      U  BBBBBBB
!     B     B    A A    C        K   K    S        U      U  B     B
!     BBBBBBB   A   A   C        K K      SSSSSSS  U      U  BBBBBBB
!     B     B  AAAAAAA  C        K   K          S  U      U  B     B
!     BBBBBBB A       A CCCCCCC  K     K  SSSSSSS  UUUUUUUU  BBBBBBB
!
! *********************************************************************
! ******* perform back-substitutions on the decomposed matrix   ******* 
! *********************************************************************

! *********************************************************************
! * this subroutine solves the linear set of equations ax = b for x,  *
! * the correction vector, where "a" is the l-u decompostion of the   *
! * original matrix,                                                  * 
! *                                                                   *
! *                       p = i - h x bo x j,                         * 
! *                                                                   *
! * i = identity matrix, h = time-step, bo = a coefficient that       * 
! * depends on the order of the integration method, and j is the      * 
! * matrix of partial derivatives. b is sent from smvgear as a        * 
! * corrected value of the first derivatives of the ordinary differ-  * 
! * ential equations. subroutine decomp.f solved for "a", the         * 
! * decomposed matrix. see press et al. (1992) numerical recipes.     *
! * cambridge university press, for a better description of the back- *
! * substitution process.                                             *
! *                                                                   *
! * this back-substitution process uses sparse-matrix techniques,     *
! * vectorizes around the grid-cell dimension, and uses no partial    *
! * pivoting. tests by sherman & hindmarsh (1980) lawrence livermore  * 
! * rep. ucrl-84102 and by us have confirmed that the removal of      *
! * partial pivoting has little effect on results.                    *
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call backsub.f from smvgear.f with                               * 
! *     ncs  = 1..ncsgas for gas chemistry                            *
! *     ncsp = ncs        for daytime   gas chem                      *  
! *     ncsp = ncs   +ics for nighttime gas chem                      *  
! *********************************************************************
!
! *********************************************************************
! *                        backsub loop # 1                           *
! * first, adjust right side of ax = b using lower triangular matrix  * 
! *********************************************************************
! sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!
! ktloop   = number of grid-cells in a grid-block
! ischan   = order of matrix
! cc2      = array holding values of decomposed matrix. 
! gloss    = array initially holding right side of equation. these
!            values are converted to the solution during back-substitution.
! kzeroa,..= arrays identifying terms in gloss array 
!
      !print*,'backsub 1'
      ij            = 1 
      do kzt    = kztlo(ncsp), kzthi(ncsp)
       i            = ikztot(kzt)
       kl5          = kbl5(  kzt)
       kh5          = kbh5(  kzt)
       kl4          = kbl4(  kzt)
       kh4          = kbh4(  kzt)
       kl3          = kbl3(  kzt)
       kh3          = kbh3(  kzt)
       kl2          = kbl2(  kzt)
       kh2          = kbh2(  kzt)
       kl1          = kbl1(  kzt)
       kh1          = kbh1(  kzt)
!
! ***********************  sum 5 terms at a time  ********************* 
!
       do kc    = kl5, kh5
        ij0         = ij
        ij1         = ij + 1
        ij2         = ij + 2 
        ij3         = ij + 3 
        ij4         = ij + 4 
        ij          = ij + 5      
        j0          = kzeroa(kc)
        j1          = kzerob(kc) 
        j2          = kzeroc(kc)
        j3          = kzerod(kc)
        j4          = kzeroe(kc)
        do k    = 1, ktloop
         gloss(k,i) = gloss(k,i)
     1              - cc2(k,ij0) * gloss(k,j0)
     2              - cc2(k,ij1) * gloss(k,j1)
     3              - cc2(k,ij2) * gloss(k,j2)
     4              - cc2(k,ij3) * gloss(k,j3)
     5              - cc2(k,ij4) * gloss(k,j4)
        enddo!k
       enddo!kc
!
! ***********************  sum 4 terms at a time  ********************* 
!
       do kc    = kl4, kh4 
        ij0         = ij
        ij1         = ij + 1
        ij2         = ij + 2 
        ij3         = ij + 3 
        ij          = ij + 4      
        j0          = kzeroa(kc)
        j1          = kzerob(kc) 
        j2          = kzeroc(kc)
        j3          = kzerod(kc)
        do k    = 1, ktloop
         gloss(k,i) = gloss(k,i)
     1              - cc2(k,ij0) * gloss(k,j0)
     2              - cc2(k,ij1) * gloss(k,j1)
     3              - cc2(k,ij2) * gloss(k,j2)
     4              - cc2(k,ij3) * gloss(k,j3)
        enddo!k
       enddo!kc
!
! ***********************  sum 3 terms at a time  ********************* 
!
       do kc    = kl3, kh3 
        ij0         = ij
        ij1         = ij + 1
        ij2         = ij + 2 
        ij          = ij + 3       
        j0          = kzeroa(kc)
        j1          = kzerob(kc) 
        j2          = kzeroc(kc)
        do  k    = 1, ktloop
         gloss(k,i) = gloss(k,i)
     1              - cc2(k,ij0) * gloss(k,j0)
     2              - cc2(k,ij1) * gloss(k,j1)
     3              - cc2(k,ij2) * gloss(k,j2)
        enddo!k
       enddo!kc
!
! ***********************  sum 2 terms at a time  ********************* 
!
       do kc    = kl2, kh2 
        ij0         = ij
        ij1         = ij + 1
        ij          = ij + 2        
        j0          = kzeroa(kc)
        j1          = kzerob(kc) 
        do k    = 1, ktloop
         gloss(k,i) = gloss(k,i)
     1              - cc2(k,ij0) * gloss(k,j0)
     2              - cc2(k,ij1) * gloss(k,j1)
        enddo!k
       enddo!kc
!
! ***********************  sum 1 term at a time  ********************** 
!
       do kc    = kl1, kh1 
        ij0         = ij
        ij          = ij + 1         
        j0          = kzeroa(kc)
        do k    = 1, ktloop
         gloss(k,i) = gloss(k,i)
     1              - cc2(k,ij0) * gloss(k,j0)
        enddo!k
       enddo!kc
      enddo!kzt
!
! *********************************************************************
! *                         backsub loop # 2                          * 
! * backsubstite with upper triangular matrix to find solution        *
! *********************************************************************
! again, sum up several terms at a time to improve vectorization.
! vdiag  = diagonal term from l-u decompostion.
! gloss  = solution on output
!
      do i       = ischan, 1, -1
       mzt           = imztot(i,ncsp)  
       if (mzt.gt.0) then
        ml5          = mbl5(  mzt)
        mh5          = mbh5(  mzt)
        ml4          = mbl4(  mzt)
        mh4          = mbh4(  mzt)
        ml3          = mbl3(  mzt)
        mh3          = mbh3(  mzt)
        ml2          = mbl2(  mzt)
        mh2          = mbh2(  mzt)
        ml1          = mbl1(  mzt)
        mh1          = mbh1(  mzt)
!
! ***********************  sum 5 terms at a time  ********************* 
!
        do mc    = ml5, mh5
         ij0         = ij
         ij1         = ij + 1
         ij2         = ij + 2 
         ij3         = ij + 3 
         ij4         = ij + 4 
         ij          = ij + 5      
         j0          = mzeroa(mc)
         j1          = mzerob(mc) 
         j2          = mzeroc(mc)
         j3          = mzerod(mc)
         j4          = mzeroe(mc)
         do k    = 1, ktloop
          gloss(k,i) = gloss(k,i)
     1               - cc2(k,ij0) * gloss(k,j0)
     2               - cc2(k,ij1) * gloss(k,j1)
     3               - cc2(k,ij2) * gloss(k,j2)
     4               - cc2(k,ij3) * gloss(k,j3)
     5               - cc2(k,ij4) * gloss(k,j4)
         enddo!k
        enddo!mc
!
! ***********************  sum 4 terms at a time  ********************* 
!
        do mc    = ml4, mh4 
         ij0         = ij
         ij1         = ij + 1
         ij2         = ij + 2 
         ij3         = ij + 3 
         ij          = ij + 4      
         j0          = mzeroa(mc)
         j1          = mzerob(mc) 
         j2          = mzeroc(mc)
         j3          = mzerod(mc)
         do k    = 1, ktloop
          gloss(k,i) = gloss(k,i)
     1               - cc2(k,ij0) * gloss(k,j0)
     2               - cc2(k,ij1) * gloss(k,j1)
     3               - cc2(k,ij2) * gloss(k,j2)
     4               - cc2(k,ij3) * gloss(k,j3)
         enddo!k
        enddo!mc
!
! ***********************  sum 3 terms at a time  ********************* 
!
        do mc    = ml3, mh3 
         ij0         = ij
         ij1         = ij + 1
         ij2         = ij + 2 
         ij          = ij + 3      
         j0          = mzeroa(mc)
         j1          = mzerob(mc) 
         j2          = mzeroc(mc)
         do k    = 1, ktloop
          gloss(k,i) = gloss(k,i)
     1               - cc2(k,ij0) * gloss(k,j0)
     2               - cc2(k,ij1) * gloss(k,j1)
     3               - cc2(k,ij2) * gloss(k,j2)
         enddo!k
        enddo!mc
!
! ***********************  sum 2 terms at a time  ********************* 
!
        do mc    = ml2, mh2 
         ij0         = ij
         ij1         = ij + 1
         ij          = ij + 2      
         j0          = mzeroa(mc)
         j1          = mzerob(mc) 
         do k    = 1, ktloop
          gloss(k,i) = gloss(k,i)
     1               - cc2(k,ij0) * gloss(k,j0)
     2               - cc2(k,ij1) * gloss(k,j1)
         enddo!k
        enddo!mc
!
! ***********************  sum 1 term at a time  ********************** 
!
        do mc    = ml1, mh1 
         ij0         = ij
         ij          = ij + 1       
         j0          = mzeroa(mc)
         do k    = 1, ktloop
          gloss(k,i) = gloss(k,i)
     1               - cc2(k,ij0) * gloss(k,j0)
         enddo!k
        enddo!mc
       endif
!      endif mzt.gt.0
!
! ***************  adjust gloss with diagonal element  **************** 
!
       do k    = 1, ktloop
        gloss(k,i) = gloss(k,i) * vdiag(k,i)
       enddo!k
      enddo!i
!
! *********************************************************************
! ******************** end of subroutine backsub **********************
! *********************************************************************
!
      !print*,'backsub 2'
      return
      end
