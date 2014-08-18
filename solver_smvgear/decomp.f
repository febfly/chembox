      subroutine decomp(ncsp)
!
      use mod_comode
      implicit none
      !input variables
      integer, intent(in) :: ncsp
      
      !local variables
      integer :: j,ijt,ij,il5,ih5,il4,ih4,il3,ih3,il2,ih2,il1,ih1
      integer :: ic,ik0,ik1,ik2,ik3,ik4,kj0,kj1,kj2,kj3,kj4,k,iar
      integer :: jl,jh,jc,ija
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!         DDDDDDD  EEEEEEE  CCCCCCC  OOOOOOO  M     M  PPPPPPP 
!         D     D  E        C        O     O  MM   MM  P     P 
!         D     D  EEEEEEE  C        O     O  M M M M  PPPPPPP 
!         D     D  E        C        O     O  M  M  M  P  
!         DDDDDDD  EEEEEEE  CCCCCCC  OOOOOOO  M     M  P
! 
! *********************************************************************
! ************** decompose the sparse matrix **************************
! *********************************************************************
!
! *********************************************************************
! * this subroutine decomposes the matrix "p" into the matrix "a" in  *
! * order to solve the linear set of equations ax = b for x, which is *
! * a correction vector. ax = b is solved in subroutine backsub.f     *
! * above, the original matrix "p" is                                 *
! *                                                                   *
! *                       p = i - h x bo x j,                         * 
! *                                                                   *
! * where i = identity matrix, h = time-step, bo = a coefficient that * 
! * depends on the order of the integration method, and j is the      * 
! * matrix of partial derivatives. see press et al. (1992) numerical  *
! * recipes cambridge university press, for a better description of   *
! * the l-u decompostion process                                      *
! *                                                                   *
! * this l-u decompostion process uses sparse-matrix techniques,      *
! * vectorizes around the grid-cell dimension, and uses no partial    *
! * pivoting. tests by sherman & hindmarsh (1980) lawrence livermore  * 
! * rep. ucrl-84102 and by us have confirmed that the removal of      *
! * partial pivoting has little effect on results.                    *
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call decomp.f from smvgear.f with                                * 
! *     ncs  = 1..ncsgas for gas chemistry                            *
! *     ncsp = ncs        for daytime   gas chem                      *  
! *     ncsp = ncs   +ics for nighttime gas chem                      *  
! *********************************************************************
!
! ktloop   = number of grid-cells in a grid-block
! ischan   = original order of matrix  
! cc2      = array of iarray units holding values of each matrix
!            position actually used. originally,
!            cc2 = p = i - delt * aset(nqq,1) * partial derivatives.   
!            however, cc2 is decomposed here
!
! *********************************************************************
! ***                  first loop of l-u decompostion               *** 
! *********************************************************************
! sum 1,2,3,4, or 5 terms at a time to improve vectorization 
!
      do j       = 1, ischan
       do ijt    = ijtlo(j,ncsp), ijthi(j,ncsp)
        ij           = ijval(ijt)
        il5          = idl5( ijt)
        ih5          = idh5( ijt)
        il4          = idl4( ijt)  
        ih4          = idh4( ijt)
        il3          = idl3( ijt)  
        ih3          = idh3( ijt)
        il2          = idl2( ijt)  
        ih2          = idh2( ijt)
        il1          = idl1( ijt)  
        ih1          = idh1( ijt)

! ********************* sum 5 terms at a time *************************
!
        do ic    = il5, ih5
         ik0         = ikdeca(ic)
         ik1         = ikdecb(ic)
         ik2         = ikdecc(ic)
         ik3         = ikdecd(ic)
         ik4         = ikdece(ic)
         kj0         = kjdeca(ic)
         kj1         = kjdecb(ic)
         kj2         = kjdecc(ic)
         kj3         = kjdecd(ic)
         kj4         = kjdece(ic)
         do k    = 1, ktloop 
          cc2(k,ij)  = cc2(k,ij) - cc2(k,ik0) * cc2(k,kj0) 
     1                           - cc2(k,ik1) * cc2(k,kj1) 
     2                           - cc2(k,ik2) * cc2(k,kj2) 
     3                           - cc2(k,ik3) * cc2(k,kj3) 
     4                           - cc2(k,ik4) * cc2(k,kj4) 
         enddo !k
       enddo!ic
!
! ********************* sum 4 terms at a time *************************
!
        do ic    = il4, ih4
         ik0         = ikdeca(ic)
         ik1         = ikdecb(ic)
         ik2         = ikdecc(ic)
         ik3         = ikdecd(ic)
         kj0         = kjdeca(ic)
         kj1         = kjdecb(ic)
         kj2         = kjdecc(ic)
         kj3         = kjdecd(ic)
         do k    = 1, ktloop 
          cc2(k,ij)  = cc2(k,ij) - cc2(k,ik0) * cc2(k,kj0) 
     1                           - cc2(k,ik1) * cc2(k,kj1) 
     2                           - cc2(k,ik2) * cc2(k,kj2) 
     3                           - cc2(k,ik3) * cc2(k,kj3) 
         enddo!k
        enddo!ic 
!
! ********************* sum 3 terms at a time *************************
!
        do ic    = il3, ih3
         ik0         = ikdeca(ic)
         ik1         = ikdecb(ic)
         ik2         = ikdecc(ic)
         kj0         = kjdeca(ic)
         kj1         = kjdecb(ic)
         kj2         = kjdecc(ic)
         do k    = 1, ktloop 
          cc2(k,ij)  = cc2(k,ij) - cc2(k,ik0) * cc2(k,kj0) 
     1                           - cc2(k,ik1) * cc2(k,kj1) 
     2                           - cc2(k,ik2) * cc2(k,kj2) 
         enddo!k
        enddo!ic
!
! ********************* sum 2 terms at a time *************************
!
        do ic    = il2, ih2
         ik0         = ikdeca(ic)
         ik1         = ikdecb(ic)
         kj0         = kjdeca(ic)
         kj1         = kjdecb(ic)
         do k    = 1, ktloop 
          cc2(k,ij)  = cc2(k,ij) - cc2(k,ik0) * cc2(k,kj0) 
     1                           - cc2(k,ik1) * cc2(k,kj1) 
         enddo!k 
        enddo!ic
! 
! ********************* sum 1 term  at a time *************************
!
        do ic    = il1, ih1
         ik0         = ikdeca(ic)
         kj0         = kjdeca(ic)
         do k    = 1, ktloop 
          cc2(k,ij)  = cc2(k,ij) - cc2(k,ik0) * cc2(k,kj0) 
         enddo!k 
        enddo!ic
!
       enddo!ijt
!
! *********************************************************************
! *    vdiag = 1 / current diagonal term of the decomposed matrix     *
! *********************************************************************
!
       iar          = jarrdiag(j,ncsp)
       do k     = 1, ktloop 
        vdiag(k,j)  = 1.0 / cc2(k,iar) 
       enddo!k
!
! *********************************************************************
! ***               second loop of decompostion                     *** 
! *********************************************************************
! jzeroa  = identifies the array position of each jloz1..jhiz1 term
!
       jl           = jloz1(j,ncsp)
       jh           = jhiz1(j,ncsp)
       do jc    = jl, jh
        ija         = jzeroa(jc)
        do k    = 1, ktloop 
         cc2(k,ija) = cc2(k,ija) * vdiag(k,j)  
        enddo!k
       enddo!jc
!
      enddo!j
!
! *********************************************************************
! ********************* end of subroutine decomp  ********************* 
! *********************************************************************
!
      return
      end
