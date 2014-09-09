      subroutine update (ncs, ncsp)
!
      use mod_comode
      implicit none
      ! Input variables
      integer, intent(in) :: ncs, ncsp
      ! Local variables
      integer :: j,nkn,kloop,i,nk,nh,ispc1,ispc2,ispc3
      real(kind=DP) :: tofday,hourang,sinfunc      
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***                        (650) 723-6836                         *** 
! *********************************************************************
!
!        U     U  PPPPPPP  DDDDDD      A    TTTTTTT  EEEEEEE
!        U     U  P     P  D     D    A A      T     E  
!        U     U  PPPPPPP  D     D   A   A     T     EEEEEEE 
!        U     U  P        D     D  AAAAAAA    T     E  
!        UUUUUUU  P        DDDDDD  A       A   T     EEEEEE
!
! *********************************************************************
! * this subroutine updates photorates and arb emissions rates for    *
! * each time-step. photorates are included in first and partial      *
! * derivative equations while emissions rates are included in first  *
! * derivate equations only. since the emissions rates are constant   * 
! * for a given time step and location (although they change each     *
! * time step and location, they are put into the first derivative    *
! * term of subfun.f only (not into partial derivative terms. every   *
! * integration time-step, emissions are recalculated.                * 
! *********************************************************************
!
! *********************************************************************
! * update photo-rates and other parmeters because the time changed.  *
! * note that a time change could correspond to either a successful   *
! * or failed step                                                    * 
! *********************************************************************
! prate    = photorate (sec-1) from most recent call to radtran.f
! prate1   = photorate (sec-1) from previous call to radtran.f
! pratkd   = prate - prate1 / tinterval   
! rrate    = prate1 + xelaps * (prate - prate1) 
! xelaps   = elapsed time during interval  
! ifprat   = 1: use scaled photorates from photrate.dat (itestgear.eq.0) 
!          = 0: use photorates from globchem.dat (itestgear > 0)  
!
! *********************************************************************
! **************          update photorates             *************** 
! ****************** interpolate between two values *******************
! *********************************************************************
!
!      print*,'check update 1'
      xelaplast             = xelaps!??GEOS-Chem do not have this!!
!
! *********************************************************************
! *      set rates where photoreaction has no active loss term        *
! *********************************************************************
! jold = mappl(jold) for inactive species 
!
      do    i            = 1, nolosp(ncsp) 
       nk                = nknlosp(i,ncs)
       nkn               = newfold(nk,ncs) 
       nh                = nkn + nallr 
       do    kloop       = 1, ktloop
        trate(kloop,nkn) =  rrate(kloop,nkn)
        trate(kloop,nh)  = -rrate(kloop,nkn)
       enddo
      enddo
!
! *********************************************************************
! *              print out chemical rates and stop                    *
! *********************************************************************
!
      if (iprates.eq.1) then
       do    i           = 1, nallrat(ncs)
        nk               = ncequat(i,ncs)
        nkn              = newfold(nk,ncs) 
        ispc1            = irm(1,nk,ncs) 
        ispc2            = irm(2,nk,ncs) 
        ispc3            = irm(3,nk,ncs)
        if (ispc3.lt.0)          ispc3 = 0 
        if (ispc1.gt.nspec(ncs)) ispc1 = 0
        if (ispc2.gt.nspec(ncs)) ispc2 = 0
        if (ispc3.gt.nspec(ncs)) ispc3 = 0
        write(iout,11)i,nk,namencs(ispc1,ncs), namencs(ispc2,ncs), 
     1                namencs(ispc3,ncs), rrate(1,nkn)
       enddo!i
       stop
      endif
 11   format(i4,1x,i4,1x,3a15,1x,1pe13.6)
!
! *********************************************************************
! ******************** end of subroutine update.f *********************
! *********************************************************************
!      print*,'check update 2'
!
      return
      end
!
