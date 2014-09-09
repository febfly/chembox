      subroutine pderiv(ncsp)
!
      use mod_comode
      implicit none
      !input variables
      integer, intent(in) :: ncsp
      !local variables
      integer :: iarry,nondiag,nondiag1,npdl,npdh,nkn,ja,jb,jc,k,iar,n
      integer :: ial
      real(kind=DP) :: fracr1   
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!       PPPPPPP  DDDDDDD  EEEEEEE  RRRRRRR  IIIIIII  V       V 
!       P     P  D     D  E        R     R     I      V     V
!       PPPPPPP  D     D  EEEEEEE  RRRRRRR     I       V   V 
!       P        D     D  E        R  R        I        V V  
!       P        DDDDDDD  EEEEEEE  R    R   IIIIIII      V 
!
! *********************************************************************
! * this subroutine puts the partial derivatives of each ordinary     *  
! * differential equation into a matrix. the form of the matrix       *
! * equation is                                                       *
! *                       p = i - h x bo x j                          *
! *                                                                   *
! * where i = identity matrix, h = time-step, bo = coefficient        *
! * corresponding to the order of the method, and j is the jacobian   *
! * matrix of partial derivatives.                                    * 
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call pderiv.f from smvgear.f with                                * 
! *     ncs  = 1..ncsgas for gas chemistry                            *
! *     ncsp = ncs        for daytime   gas chem                      *  
! *     ncsp = ncs   +ics for nighttime gas chem                      *  
! *********************************************************************
!
! *********************************************************************
! *                        initialize matrix                          * 
! *********************************************************************
! cc2      = array of iarray units holding values of each maxtrix
!            position actually used.
!            cc2 = p = i - delt * aset(nqq,1) * partial derivatives.   
! urate    = term of jacobian (j) = partial derivative 
! iarray   = total number of matrix positions filled after mat. processes
! irma,b,c = species # of each reactant 
! ischan   = original order of matrix  
! ktloop   = number of grid-cells in a grid-block
! nondiag1 = 1 + # of final matrix positions, excluding diagonal terms,
!            filled after all matrix processes.
! npderiv  = counter of number of times this routine is called
! r1delt   = -aset(nqq,1) * time step = -coefficient of method * dt 
! rrate    = reaction rate coefficient
!
! example of how partial derivatives are placed in an array:
! ---------------------------------------------------------
!
! species:          a,   b,   c
! concentrations:  [a], [b], [c]       
!
! reactions:    1) a           --> b      j 
!               2) a  + b      --> c      k1 
!               3  a  + b + c  --> d      k2  
!
! first         d[a] / dt  =  -j[a] - k1[a][b] - k2[a][b][c]
! derivatives:  d[b] / dt  =  +j[a] - k1[a][b] - k2[a][b][c]
!               d[c] / dt  =        + k1[a][b] - k2[a][b][c]
!               d[d] / dt  =                   + k2[a][b][c]
!
! predictor matrix (p) = i - h * b * j:
!  j = jacobian matrix of partial derivates
!  i = identity matrix
!  h = time-step
!  b = coefficient of method
!  r = h * b = -r1delt  
!
!                a                 b                 c                 d 
!     ___________________________________________________________________
!   | 
! a | 1-r(-j-k1[b]-k2[b][c])  -r(-k1[a]-k2[a][c])   -r(-k2[a][b])       0    
!   | 
! b |  -r(+j-k1[b]-k2[b][c]) 1-r(-k1[a]-k2[a][c])   -r(-k2[a][b])       0  
!   | 
! c |  -r(  +k1[b]-k2[b][c])  -r(+k1[a]-k2[a][c])  1-r(-k2[a][b])       0 
!   |
! d |  -r(        +k2[b][c])  -r(      +k2[a][c])   -r(+k2[a][b])       1 
!    
!
! *********************************************************************
! *********         calculate partial derivatives            **********
! *********    and sum up partial derivative loss terms      ********** 
! *********************************************************************
!
      !print*,'check pderiv 1'
!      npderiv             = npderiv + 1
      iarry               = iarray(ncsp) 
      nondiag             = iarry - ischan  
      nondiag1            = nondiag + 1
      nfdh1               = nfdh2 + ioner(ncsp) 
      npdl                = npdlo(ncsp)
      npdh                = npdhi(ncsp)
!
! *********************************************************************
! *    partial derivatives for rates with three active loss terms     *
! *********************************************************************
!
      do     nkn       = 1, nfdh3
       ja              = irma(nkn)
       jb              = irmb(nkn)
       jc              = irmc(nkn)
       do     k        = 1, ktloop
        urate(k,nkn,1) = rrate(k,nkn) * cnew(k,jb) * cnew(k,jc) 
        urate(k,nkn,2) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jc)
        urate(k,nkn,3) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jb)
       enddo
      enddo
!
! *********************************************************************
! *    partial derivatives for rates with two active loss terms       *
! *********************************************************************
!
      do     nkn       = nfdl2, nfdh2
       ja              = irma(nkn)
       jb              = irmb(nkn)
       do     k        = 1, ktloop
        urate(k,nkn,1) = rrate(k,nkn) * cnew(k,jb) 
        urate(k,nkn,2) = rrate(k,nkn) * cnew(k,ja) 
       enddo
      enddo
!
! *********************************************************************
! *    partial derivatives for rates with one active loss term        *
! *********************************************************************
!
      do     nkn       = nfdl1, nfdh1 
       do     k        = 1, ktloop
        urate(k,nkn,1) = rrate(k,nkn) 
       enddo
      enddo
!
! *********************************************************************
! * put partial derivatives production and loss terms in matrix array * 
! *********************************************************************
! fracpl = -1. for all reactants 
!        = +1. or +fraction for all products
!   
      do     iar     = 1, nondiag
       do     k      = 1, ktloop
        cc2(k,iar)   = 0.
       enddo
      enddo
!
      do     iar     = nondiag1, iarry
       do     k      = 1, ktloop
        cc2(k,iar)   = 1.
       enddo
      enddo
!
      do     n       = npdl, npdh
       nkn           = nkpdterm(n) 
       iar           = ipospd(  n)
       ial           = iialpd(  n)
       fracr1        = fracpl(  n) * r1delt  
       do     k      = 1, ktloop 
        cc2(k,iar)   = cc2(k,iar) + fracr1 * urate(k,nkn,ial)   
       enddo
      enddo
!
! *********************************************************************
! ********************* end of subroutine pderiv **********************
! *********************************************************************
!
      !print*,'check pderiv 2'
      return 
      end
