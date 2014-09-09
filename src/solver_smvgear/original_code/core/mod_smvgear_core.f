!===========================================================
! Module: mod_smvgear_core
! Core functions that realize the smvgear algorithmn
! Original code written by Mark Jacobson
! Organized into a module by Yuzhong Zhang, 9/8/2014 
!===========================================================
      module mod_smvgear_core
      use mod_smvgear_comode
      implicit none
      !public functions
      public :: smvgear
      public :: jsparse

      !private functions
      private :: ksparse
      private :: backsub
      private :: decomp
      private :: pderiv
      private :: subfun
      private :: update


!==========================================================
      subroutine smvgear(ncs, ifsun,flag)
      !input variables
      integer, intent(in) :: ncs
      integer, intent(in) :: ifsun !1. daytime; 2. nighttime
      integer, intent(out):: flag
      
      !local variables
      integer :: ncsp
      integer :: jfail,ischan1,iabove,kloop,idoub,jrestar,jnew,ifsuccess
      integer :: k,jspc,k1,k2,k3,k4,k5,nqqold,jeval,js1,nqqisc
      integer :: lloopa,lloopb,jloop,mloop,m1,m2,jold,i1,j,i,j1,j2,j3,j4
      integer :: j5,l3,jb,jg1,kstepisc,nqisc,i2,nslp,kstep
      !integer :: ifail, lfail,nfail

      real(kind=DP) :: nylowdec,order,hrmax,yfac,errinit,reltol1,reltol2,reltol3
      real(kind=DP) :: abtoler1,abtoler2,hratio,asn1,rdelmax,cnw,cnewylow,errymax
      real(kind=DP) :: rmstop,delt1,enqq,eup,edwn,conp3,conp2,conp1,hmtim,rdelta
      real(kind=DP) :: conc3j3,conc4j4,conc10j5,conc5j5,drate,rmserrp,der2max
      real(kind=DP) :: rmsrat,dcon,rdeltup,asnqqj,der3max,rdeltsm,der1max,rdeltdn
      real(kind=DP) :: consmult,tinterval 
      
      ! Variables from "comode.h" which are only ever used in "smvgear.f"
      ! Remove them from "comode.h" and the THREADPRIVATE declarations
      ! (bmy, 7/28/03)
      integer   :: nsteps
      integer, dimension(KBLOOP,5) :: kgrp
      integer, dimension(KBLOOP)   :: iabovk
      real(kind=DP), dimension(KBLOOP) :: dely, errhold, yabst
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
! *********************************************************************
! *********************************************************************
!
!  SSSSSSS   M     M  V       V  GGGGGGG   EEEEEEE      A      RRRRRRR
!  S         MM   MM   V     V   G         E           A A     R     R  
!  SSSSSSS   M M M M    V   V    G  GGGG   EEEEEEE    A   A    RRRRRRR
!        S   M  M  M     V V     G     G   E         AAAAAAA   R  R
!  SSSSSSS   M     M      V      GGGGGGG   EEEEEEE  A       A  R    R
!
! *********************************************************************
!                    version:      smvgear ii
!                    last update:  august, 1997
! *********************************************************************
!
! *********************************************************************
! * smvgear is a gear-type integrator that solves first-order ordin-  *
! * ary differential equations with initial value boundary conditions.*
! * smvgear differs from an original gear code in that it uses sparse *
! * matrix and vectorization techniques to improve speed. much        * 
! * of the speed up in this program is due to sparse matrix           *
! * techniques and vectorization.                                     *
! *                                                                   *
! * this version includes re-ordering of grid-cells prior to each     *
! * time-interval. the purpose of the reordering is to group cells    *
! * with stiff equations together and those with non-stiff equations  *
! * this reordering can save signifcant computer time                 *
! * (e.g. speed the code by a factor of two or more), depending on    *
! * the variation in stiffness throughout the grid-domain. when the   *
! * stiffness is the same throughout the grid-domain (e.g. if all     *
! * concentrations and rates are the same), then re-ordering is       *
! * unnecessary and will not speed solutions.                         *
! *                                                                   *
! * this version includes a variable absolute error tolerance.        *
! * the absolute tolerance is recalculated every few gear time steps. *
! *                                                                   *
! * this version contains different sets of chemistry for             *
! * different regions of the atmosphere. thus, urban, free trop-      *
! * ospheric, and stratospheric chemistry can be solved during the    *
! * same model run.                                                   * 
! *                                                                   *
! * references:                                                       *
! * -----------                                                       *
! *                                                                   * 
! * jacobson m. z. (1999) fundamentals of atmospheric modeling.       *
! *  cambridge university press, new york, 656 p.                     *
! *                                                                   * 
! * jacobson m. z. (1998) improvement of smvgear ii on vector and     *
! *  scalar machines through absolute error tolerance control.        *    
! *  atmos. environ. 32, 791 - 796                                    *
! *                                                                   * 
! * jacobson m. z. (1995) computation of global photochemistry        *
! *  with smvgear ii. atmos. environ., 29a, 2541 - 2546               *
! *                                                                   *
! * jacobson m. z. (1994) developing, coupling, and applying a gas,   *
! *  aerosol, transport, and radiation model to studying urban        *
! *  and regional air pollution. ph. d. thesis, university of         *
! *  california, los angeles.                                         *
! *                                                                   *
! * jacobson m. z. and turco r. p. (1994) smvgear: a sparse-          * 
! *  matrix, vectorized gear code for atmospheric models.             *
! *  atmos. environ. 28a, 273 - 284.                                  * 
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call smvgear from physproc for gas chem w/ ncs = 1..ncsgas       *
! *                                                                   *
! *********************************************************************
! *                                                                   *
! * the origins of the gear integrator used in smvgear are found in   *
! *                                                                   *    
! * gear c. w. (1971) numerical initial value problems in ordinary    *  
! *  differential equations. prentice-hall, nj, pp. 158-166.          * 
! *                                                                   *    
! *********************************************************************
! *                                                                   *      
! * finally, in subroutine smvgear.f, the following ideas originated  *
! *   from lsodes, the livermore solver for ordinary differential     *
! *   with sparse matrices (hindmarsh a. c. and sherman a. h.):       *
! *                                                                   *      
! *  (a) predicting the first time-step;                              *
! *  (b) determining corrector convergence differently from in        *
! *      gear's original code (goc)                                   *
! *  (c) determining error differently from in goc                    *
! *  (d) summing up the pascal matrix differently from in goc         *      
! *                                                                   *      
! * references for the 1987 lsodes version include:                   *
! *                                                                   *      
! * sherman a. h. and hindmarsh a. c. (1980) gears: a package for     *
! *  the solution of sparse, stiff ordinary differential equations.   *
! *  lawrence livermore laboratory report ucrl-84102.                 *   
! *                                                                   *      
! * hindmarsh a. c. (1983) odepack, a systematized collection of      *
! *  ode solvers. in scientific computing, r.s. stepleman et al.,     *
! *  eds., north-holland, amsterdam, pp. 55 - 74.                     *
! *                                                                   *      
! *********************************************************************
!
! *********************************************************************
! *************** here are some parameter definitions *****************
! *********************************************************************
!                                                                          
! abst2     = 1. / timeinterval**2   (sec-2) (set in reader.f) 
! asn1      = the value of aset(nqq,1)
! cest      = stores value of dtlos when idoub = 1
! chold     = 1 / (reltol * cnew + abtol). multiply
!             chold by local errors in different error tests.
! cnew      = stores concentration (y [estimated])
! conc      = an array of length ischan * (maxord+1) that carries the
!             derivatives of cnew, scaled by delt**j/factorial(j),
!             where j is the j-th derivative. j varies from 1 to nqq,
!             which is the current order of the method.
!             e.g. conc(jspc,2) stores delt * y' (estimated)                   
! delt      = current time-step (s) length during a time-interval 
! drate     = parameter which used to determine whether convergence 
!             has occurred
! dtlos     = an array of length ischan, used for the accumulated
!             corrections.  on a successful return, dtlos(kloop,i) contains
!             the estimated one-step local error in cnew.
! edwn      = pertst**2 * order for one order lower than current order
! enqq      = pertst**2 * order for current order
! errmax    = relative error tolerance (see chold). set in m.dat.
!             eps should be < 1.0. for speedy and reliable results, 
!             10**-3 is reasonable. for many decimal places of accuracy,
!             decrease eps. 
! eup       = pertst**2 * order for one order higher than current order 
! fracdec   = fraction the time-step is decreased if convergence test fails
! gloss     = value of first derivatives on output from subfun.
!           = right-side of equation on input to backsub.f 
!           = error term (solution from backsub.f) on output from backsub
! hmax      = the maximum allowable value of delt
! hmin      = the minimum allowable value of delt
! hrmax     = maximum relative change in delt*aset(1) before pderiv is called.
! hratio    = relative change in delt * aset(1) each change in step or order
!             when abs(hratio-1) > hrmax, reset jeval = 1 to call pderiv
! iabovk    = number of species whose concentrations are larger than yabst
! idoub     = records the number of steps since the last change in step size   
!             or order.  it must be at least kstep = nqq+1 before doubling is 
!             allowed. 
! ifail     = number of times the corrector failed to converge while the
!             jacobian was old (pderiv not called during the last test)
! ifsuccess = identifies whether step is successful (=1) or not (=0)
! ifsun     = identifies whether sun is up (=1) or down (=2)
! ischan    = the number of first-order equations to solve = # of species = 
!             order of original matrix. ischan has a different value
!             for day and night and for gas- chemistry.
! isreord   = 1: calc initial stiffness before running code to reorder cells
!                in this case, use photorates for end of time-interval
!           = 0: do normal calculations
! jeval     = 1  --> call pderiv the next time through the corrector steps.
!           = 0  --> last step was successful and do not need to call pderiv
!           = -1 --> pderiv just called, and do not need to call again
!             until jeval switched to 1. 
! jrestar   = counts number of times smvgear starts over at order 1
!             because of excessive failures.
! lfail     = number of times the accumulated error test failed
! kstep     = nqq + 1
! ktloop    = number of grid-cells in a grid-block
! maxord    = the maximum allowable order of the integration method
! mbetween  = the maximum allowable number of steps between calls to pderiv
! mstep     = the maximum allowable number of corrector iterations
! ncs       = 1..ncsgas for gas chemistry                            
! ncsp      = ncs       for daytime   gas chem            
!           = ncs + ics for nighttime gas chem           
! nfail     = number of times correcter fails to converge after pderiv
!             was just called
! npderiv   = total number of times that matrix is evaluated (pderiv)
! npdtot    = number of calls to pderiv routine, over all time
! nsftot    = number of calls to subfun routine, over all time
! nslp      = the last time-step number during which pderiv was called 
! nsttot    = total number of successful time-steps, over all time 
! nsubfun   = total number of times subfun is called
! nsteps    = total number of successful time-steps taken
! nqq       = order of the integration method. it varies between 1 and maxord. 
! nqqisc    = nqq * ischan
! nqqold    = value of nqq during last time-step
! order     = floating point value of ischan, the order of number of odes.
! pderiv    = name of routine to evaluate the jacobian matrix (j) 
!             and p = i - delt * aset(1) * j 
! perts2    = coefficients used in selecting the step and order (see
!             jsparse.f) note that perts2 = original pertst**2     
! rdelmax   = the maximum factor by which delt can be increased in a single 
!             step.  as in lsodes, set it to 1e4 initially to compensate 
!             for the small initial delt, but then set it to 10 after 
!             successful steps and 2 after unsuccessful steps
! rdelt     = factor (time-step ratio) by which we increase or decrease delt
! rdeltdn   = time-step ratio at one order lower than current order 
! rdeltsm   = time-step ratio at current order 
! rdeltup   = time-step ratio at one order higher than current order 
! rmsrat    = ratio of current to previous rms scaled error. if this
!             ratio decreases, then convergence is occuring.
! subfun    = name of routine to solve first derivatives.
!           = evaluates derivatives in the special form f = y'(est)
!           = f(x,y,estimated), where f is the right hand side of the
!             differential equation. 
! tinterval = total seconds in a time-interval
! timremain = remaining time in an interval 
! told      = stores the last value of xelaps in case the current step fails
! xelaps    = elapsed time in an interval (s)
! abtol     = absolute error tolerance 
!             if abtol is too small, then integration will take too long.
!             if abtol too large, convergence will be too easy and errors
!             will accumulate, the time-step may be cut too small, and
!             the integration may stop (delt < hmin or floating point
!             exception in decomp.f).
!             typical gas-phase values of abstol are 10**3 cm-3
!             typical aq -phase values of abstol are 10**-13 to 10**-15 m l-1
! yfac      = 1.0 originially, but is decreased if excessive failures occur
!             in order to reduce absolute error tolerance 
! *********************************************************************
!
      ncsp = ncs + (ifsun-1) * ICS
      
!      nsubfun   = 0 
!      npderiv   = 0 
      nsteps    = 0
!      ifail     = 0 
      jfail     = 0 
!      lfail     = 0 
!      nfail     = 0 
      nylowdec  = 0
      tinterval = timeintv(ncs)
      ischan    = ischang( ncs)
      ischan1   = ischan - 1
!
      order     = float(ischan)
!
      iabove    = int(order * 0.4)       
!
      do     kloop   = 1, ktloop
       iabovk(kloop) = iabove
      enddo
!
      hrmax     = 0.3
      hmax      = hmaxuse( ncsp)
      yfac      = 1.0 
      errinit   = min(errmax(ncs),1.0d-03)
!
! *********************************************************************
!            start time interval or re-enter after total failure
! *********************************************************************
!
 120  idoub     = 2
      nslp      = mbetween
      jrestar   = 0 
      delt      = 0.
      xelaps    = 0. 
      xelaplast = -1. 
      told      = 0.
      timremain = tinterval
      reltol1   = yfac   / errinit
      reltol2   = yfac   / errmax(ncs)
      reltol3   = 1.     / errmax(ncs)
      abtoler1  = abtol(6,ncs) * reltol1
      abtoler2  = abtol(6,ncs) * reltol2
!
! *********************************************************************
!                  initialize concentration array 
! *********************************************************************
! corig = original concentrations, which do not change in smvgear
! cnew  = final concentrations, calculated in smvgear
!
      do     jnew         = 1, ischan
       do     kloop       = 1, ktloop 
        cnew( kloop,jnew) = corig(kloop,jnew)  
       enddo
      enddo
!
! *********************************************************************
!  re-enter here if total failure or if restarting with new cell block 
! *********************************************************************
!
 140  hratio    = 0. 
      asn1      = 1.
      ifsuccess = 1 
      rdelmax   = 1.0e+04
! *********************************************************************
!                         initialize photrates 
! *********************************************************************
!
      if (ifsun.eq.1)                            call update(ncs, ncsp)
!
! *********************************************************************
!               initialize first derivative for chemistry 
! *********************************************************************
!
                                                 call subfun(ncsp)
!
! *********************************************************************
!                determine initial absolute error tolerance 
! *********************************************************************
! iabovk  = number of species whose concentrations are larger than yabst
! isreord = 1: calc initial stiffness before running code to reorder cells
!              in this case, use photorates for end of time-interval
!         = 2: do normal calculations
! kgrp    = counts number of concentrations above abtol(i), i = 1..  
! yabst   = absolute error tolerance (molec. cm-3 for gases) 
! abtol   = pre-defined absolute error tolerances 
!

      do     kloop    = 1, ktloop
       errhold(kloop) = 0.
      enddo
!
      !!if (isreord.ne.1) then
!
       do     k              = 1, 5        
        do     kloop         = 1, ktloop
         kgrp(kloop,k)       = 0
        enddo
       enddo
!
       do     jspc           = 1, ischan
        do     kloop         = 1, ktloop
         cnw                 = cnew(kloop,jspc)
         if (cnw.gt.abtol(1,ncs)) then
          kgrp(kloop,1)      = kgrp(kloop,1) + 1
         elseif (cnw.gt.abtol(2,ncs)) then
          kgrp(kloop,2)      = kgrp(kloop,2) + 1
         elseif (cnw.gt.abtol(3,ncs)) then
          kgrp(kloop,3)      = kgrp(kloop,3) + 1
         elseif (cnw.gt.abtol(4,ncs)) then
          kgrp(kloop,4)      = kgrp(kloop,4) + 1
         elseif (cnw.gt.abtol(5,ncs)) then
          kgrp(kloop,5)      = kgrp(kloop,5) + 1
         endif
        enddo!jspc
       enddo!kloop
!
       do     kloop         = 1, ktloop
        k1                  = kgrp(kloop,1)
        k2                  = kgrp(kloop,2) + k1
        k3                  = kgrp(kloop,3) + k2 
        k4                  = kgrp(kloop,4) + k3 
        k5                  = kgrp(kloop,5) + k4 
        if (k1.gt.iabovk(kloop)) then
         yabst(kloop)       = abtol(1,ncs)
        elseif (k2.gt.iabovk(kloop)) then
         yabst(kloop)       = abtol(2,ncs)
        elseif (k3.gt.iabovk(kloop)) then
         yabst(kloop)       = abtol(3,ncs)
        elseif (k4.gt.iabovk(kloop)) then
         yabst(kloop)       = abtol(4,ncs)
        elseif (k5.gt.iabovk(kloop)) then
         yabst(kloop)       = abtol(5,ncs)
        else
         yabst(kloop)       = abtol(6,ncs)
        endif
       enddo!kloop
!
       do     jspc      = 1, ischan     
        do     kloop    = 1, ktloop 
         cnewylow       = cnew(kloop,jspc) + yabst(kloop) * reltol1
         errymax        = gloss(kloop,jspc) / cnewylow
         errhold(kloop) = errhold(kloop) + errymax * errymax
        enddo
       enddo
!
      !!else
!
! *********************************************************************
!          use lowest absolute error tolerance when reordering 
!          if reordering, set errmx2 then return to physproc.f 
! *********************************************************************
! abtoler1 = yfac * abtol(6,ncs) / min(errmax,1.0e-03) 
! 
       !!do 144 jspc      = 1, ischan     
        !!do 143 kloop    = 1, ktloop 
         !!errymax        = gloss(kloop,jspc)/(cnew(kloop,jspc)+abtoler1)
         !!errhold(kloop) = errhold(kloop) + errymax * errymax
!! 143    continue
!! 144   !!continue
!
       !!if (isreord.eq.1) then 
        !!do 150 kloop           = 1, ktloop 
         !!errmx2(jlooplo+kloop) = errhold(kloop)
!! 150    !!continue
!
!!        return
!!       endif
!!      endif
!
! *********************************************************************
!               calculate initial time step size (s) 
! *********************************************************************
! sqrt(errhold / [errinit * order]) = rmsnorm of error scaled to errinit 
!                                * cnew + abtol/reltol
!
      rmstop         = 0.
!
      do     kloop   = 1, ktloop
       if (errhold(kloop).gt.rmstop) rmstop = errhold(kloop)   
      enddo
!
      delt1        = sqrt(errinit / (abst2(ncs) + rmstop / order)) 
      delt         = max(min(delt1,timremain,hmax),hmin) 
!
! *********************************************************************
!                      set initial order to 1
! *********************************************************************
!
      nqqold       = 0
      nqq          = 1 
      jeval        = 1
      rdelt        = 1.0
!
! *********************************************************************
! *   store initial concentration and first derivatives x time-step   * 
! *********************************************************************
!
      do    jspc        = 1, ischan 
       js1               = ischan + jspc
       do     kloop      = 1, ktloop
        conc(kloop,jspc) = cnew(kloop,jspc)
        conc(kloop,js1)  = delt * gloss(kloop,jspc) 
       enddo
      enddo
!
! *********************************************************************
! ** update coefficients of the order. nqq is the order. aset and    **
! ** perts2 are defined in subroutine ksparse. note that perts2      **
! ** is the original pertst**2                                       **
! *********************************************************************
!
 170  if (nqq.ne.nqqold) then
       nqqold                = nqq 
       kstep                 = nqq + 1
       hratio                = hratio * aset(nqq,1) / asn1
       asn1                  = aset(nqq,1)  
       enqq                  = perts2(nqq,1) * order
       eup                   = perts2(nqq,2) * order 
       edwn                  = perts2(nqq,3) * order  
       conp3                 = 1.4 / ( eup**enqq3(nqq))
       conp2                 = 1.2 / (enqq**enqq2(nqq))
       conp1                 = 1.3 / (edwn**enqq1(nqq))
       nqqisc                = nqq * ischan 
      endif
!
! *********************************************************************
!   limit size of rdelt, then recalculate new time step and update
! hratio. use hratio to determine whether pderiv should be called again 
! *********************************************************************
!
      hmtim         = min(hmax,timremain) 
      rdelt         = min(rdelt,rdelmax,hmtim/delt)
      delt          = delt   * rdelt 
      hratio        = hratio * rdelt 
      xelaps        = xelaps + delt  
!
      if (abs(hratio-1.0).gt.hrmax.or.nsteps.ge.nslp) jeval = 1
!
! *********************************************************************
!      if time step < hmin, tighten absoloute error tolerance and 
!          restart integration at beginning of time interval
! *********************************************************************
!
      if (delt.lt.hmin) then
       write(iout,12) delt,ktloop,ncs,time,
     1                 timremain,yfac,errmax(ncs)
       nylowdec  = nylowdec + 1
       yfac      = yfac  * 0.01
!
       if (nylowdec.eq.10) then
!        lloopa      = 1
!        lloopb      = ktloop 
!        write(iout,13)
!
!        do     kloop = 1, ktloop 
!         jloop       = jreorder(jlooplo+kloop)
!         k           = (jloop - 1) / nloop + 1
!         mloop       = jloop - (k - 1) * nloop
!         m1          = (mloop - 1) / nlong + 1 
!         m2          = mloop - (m1 - 1) * nlong
!         write(iout,16) m1, m2, k, errhold(kloop)
!        enddo
!
!        do     jnew = 1, ischan 
!         jold       = inewold(jnew,ncs) 
!         write(iout,17) jnew,ncs,namencs(jold,ncs),corig(lloopa,jnew), 
!     1                   corig(lloopb,jnew)  
!        enddo
        flag = 2
        return
       endif  
!
       goto 120
      endif
! 
! *********************************************************************
! * if the delt is different than during the last step (if rdelt ne   *
! * 1), then scale the derivatives                                    *
! *********************************************************************
!
      if (rdelt.ne.1.0) then
       rdelta            = 1.0
       i1                = 1
       do     j          = 2, kstep 
        rdelta           = rdelta * rdelt 
        i1               = i1 + ischan 
        do     i         = i1, i1 + ischan1
         do     kloop    = 1, ktloop
          conc(kloop,i)  = conc(kloop,i) * rdelta  
         enddo
        enddo
       enddo
      endif
!
! *********************************************************************
! * update photo rates because the time changed.                      *
! * note that a time change could correspond to either a successful   *
! * or failed step                                                    * 
! *********************************************************************
!
      if (ifsun.eq.1.and.xelaps.ne.xelaplast)  call update(ncs,ncsp)
!
! *********************************************************************
! * if the last step was successful, reset rdelmax = 10 and update    *
! * the chold array with current values of cnew.                      * 
! *********************************************************************
!
      if (ifsuccess.eq.1) then
       rdelmax             = 10.
!
! *********************************************************************
!                determine new absolute error tolerance 
! *********************************************************************
! kgrp    = counts number of concentrations above abtol(i), i = 1..  
! yabst   = absolute error tolerance (molec. cm-3 for gases) 
! abtol   = pre-defined absolute error tolerances 
!
       if (mod(nsteps,3).eq.2) then
        do     k              = 1, 5        
         do     kloop         = 1, ktloop
          kgrp(kloop,k)       = 0
         enddo
        enddo
!
        do     jspc           = 1, ischan
         do     kloop         = 1, ktloop
          cnw                 = cnew(kloop,jspc)
          if (cnw.gt.abtol(1,ncs)) then
           kgrp(kloop,1)      = kgrp(kloop,1) + 1
          elseif (cnw.gt.abtol(2,ncs)) then
           kgrp(kloop,2)      = kgrp(kloop,2) + 1
          elseif (cnw.gt.abtol(3,ncs)) then
           kgrp(kloop,3)      = kgrp(kloop,3) + 1
          elseif (cnw.gt.abtol(4,ncs)) then
           kgrp(kloop,4)      = kgrp(kloop,4) + 1
          elseif (cnw.gt.abtol(5,ncs)) then
           kgrp(kloop,5)      = kgrp(kloop,5) + 1
          endif
         enddo
        enddo
!
        do     kloop         = 1, ktloop
         k1                  = kgrp(kloop,1)
         k2                  = kgrp(kloop,2) + k1
         k3                  = kgrp(kloop,3) + k2 
         k4                  = kgrp(kloop,4) + k3 
         k5                  = kgrp(kloop,5) + k4 
         if (k1.gt.iabovk(kloop)) then
          yabst(kloop)       = abtol(1,ncs)
         elseif (k2.gt.iabovk(kloop)) then
          yabst(kloop)       = abtol(2,ncs)
         elseif (k3.gt.iabovk(kloop)) then
          yabst(kloop)       = abtol(3,ncs)
         elseif (k4.gt.iabovk(kloop)) then
          yabst(kloop)       = abtol(4,ncs)
         elseif (k5.gt.iabovk(kloop)) then
          yabst(kloop)       = abtol(5,ncs)
         else
          yabst(kloop)       = abtol(6,ncs)
         endif
        enddo
       endif
!
       do     jspc         = 1, ischan
        do     kloop       = 1, ktloop
         chold(kloop,jspc) = reltol3 / (max(cnew(kloop,jspc),0d0) 
     1                     + yabst(kloop) * reltol2) 
        enddo
       enddo 
!
      endif 
!     endif ifsuccess.eq.1
!
! *********************************************************************
! * compute the predicted concentration and derivatives by multiply-  *
! * ing previous values by the pascal triangle matrix.                * 
! *********************************************************************
! this set of operations is equivalent to the reverse of loop 419.
! the expansion of the pascal triangle matrix was calculated by b. schwartz.
! the first derivative multiplied by the time step is the sum 
! of terms added to conc(kloop,i)
!
      if (nqq.eq.1) then
       do     i          = 1, ischan
        j                = i + ischan
        do     kloop     = 1, ktloop
         conc(kloop,i)   = conc(kloop,i) + conc(kloop,j)
        enddo
       enddo
! 
      elseif (nqq.eq.2) then
!
       do     i         = 1, ischan
        j1              = i  + ischan
        j2              = j1 + ischan 
        do     kloop    = 1, ktloop
         conc(kloop, i) = conc(kloop, i) + conc(kloop,j1)
     1                  +                  conc(kloop,j2) 
         conc(kloop,j1) = conc(kloop,j1) + conc(kloop,j2) * 2. 
        enddo
       enddo
! 
      elseif (nqq.eq.3) then
!
       do     i         = 1,   ischan
        j1              = i  + ischan
        j2              = j1 + ischan
        j3              = j2 + ischan
        do     kloop    = 1, ktloop
         conc3j3        = conc(kloop,j3) * 3. 
         conc(kloop, i) = conc(kloop, i) + conc(kloop,j1) 
     1                  + conc(kloop,j2) + conc(kloop,j3)
         conc(kloop,j1) = conc(kloop,j1) + conc(kloop,j2) * 2. + conc3j3 
         conc(kloop,j2) = conc(kloop,j2) + conc3j3
        enddo
       enddo
!
      elseif (nqq.eq.4) then
!
       do     i         = 1, ischan
        j1              = i  + ischan
        j2              = j1 + ischan
        j3              = j2 + ischan
        j4              = j3 + ischan
        do     kloop    = 1, ktloop
         conc3j3        = conc(kloop,j3) * 3. 
         conc4j4        = conc(kloop,j4) * 4. 
         conc(kloop, i) = conc(kloop, i) + conc(kloop,j1) 
     1                  + conc(kloop,j2) + conc(kloop,j3) 
     2                  + conc(kloop,j4)
         conc(kloop,j1) = conc(kloop,j1) + conc(kloop,j2) * 2. + conc3j3 
     1                                   + conc4j4
         conc(kloop,j2) = conc(kloop,j2) + conc3j3 + conc(kloop,j4) * 6. 
         conc(kloop,j3) = conc(kloop,j3) + conc4j4
        enddo
       enddo
!
      elseif (nqq.eq.5) then
!
       do     i         = 1, ischan
        j1              = i  + ischan
        j2              = j1 + ischan
        j3              = j2 + ischan
        j4              = j3 + ischan
        j5              = j4 + ischan
        do     kloop    = 1, ktloop
         conc3j3        = conc(kloop,j3) * 3. 
         conc4j4        = conc(kloop,j4) * 4. 
         conc5j5        = conc(kloop,j5) * 5. 
         conc10j5       = conc5j5 + conc5j5 
         conc(kloop, i) = conc(kloop, i) + conc(kloop,j1) 
     1                  + conc(kloop,j2) + conc(kloop,j3) 
     2                  + conc(kloop,j4) + conc(kloop,j5)
         conc(kloop,j1) = conc(kloop,j1) + conc(kloop,j2) * 2. + conc3j3 
     1                                   + conc4j4 + conc5j5
         conc(kloop,j2) = conc(kloop,j2) + conc3j3 + conc(kloop,j4) * 6. 
     1                                   + conc10j5
         conc(kloop,j3) = conc(kloop,j3) + conc4j4 + conc10j5
         conc(kloop,j4) = conc(kloop,j4) + conc5j5
        enddo
       enddo
      endif
!
! *********************************************************************
! ************************** correction loop **************************
! * take up to 3 corrector iterations. test convergence by requiring  *
! * that changes be less than the rms norm weighted by chold.         * 
! * accumulate the correction in the array dtlos(). it equals the     *
! * the j-th derivative of conc() multiplied by delt**kstep /         *
! * (factorial(kstep-1)*aset(kstep)); thus, it is proportional to the *
! * actual errors to the lowest power of delt present (delt**kstep)   *
! *********************************************************************
!
 220  l3                  = 0
      do     jspc         = 1, ischan
       do     kloop       = 1, ktloop
        cnew(kloop,jspc)  = conc(kloop,jspc)
        dtlos(kloop,jspc) = 0.
       enddo
      enddo
!
! *********************************************************************
! * if jeval = 1, re-evaluate predictor matrix p = i - h * aset(1) *j *
! * before starting the corrector iteration. after calling pderiv,    * 
! * set jeval = -1 to prevent recalling pderiv unless necessary later.*
! * call decomp to decompose the matrix                               *
! *********************************************************************
!
      if (jeval.eq.1) then
       r1delt   = -asn1 * delt
!
                                                call pderiv(ncsp)
!
                                                call decomp(ncsp) 
       jeval    = -1 
       hratio   = 1.0
       nslp     = nsteps + mbetween
       drate    = 0.7
      endif
!
! *********************************************************************
! *   evaluate the first derivative using corrected values of cnew    * 
! *********************************************************************
!
 270                                            call subfun(ncsp)
!
! *********************************************************************
! * in the case of the chord method, compute error (gloss) from the   * 
! * corrected calculation of the first derivative                     * 
! *********************************************************************
!
      do     jspc         = 1,  ischan
       j                  = jspc + ischan
       do     kloop       = 1, ktloop
        gloss(kloop,jspc) = delt * gloss(kloop,jspc)  
     1                    - (conc(kloop,j) + dtlos(kloop,jspc))
       enddo
      enddo
!
! *********************************************************************
! * solve the linear system of equations with the corrector error.    *
! * backsub.f solves backsubstitution over matrix of partial derivs.  * 
! *********************************************************************
!
                                                call backsub(ncsp)
!
! *********************************************************************
! * sum-up the accumulated error, correct the concentration with the  * 
! * error, and begin to calculate the rmsnorm of the error relative   *  
! * to chold.                                                         *
! *********************************************************************
!
      do    kloop     = 1, ktloop
       dely(kloop)     = 0.
      enddo
!
      if (asn1.eq.1.0) then 
       do     i         = 1, ischan
        do     kloop    = 1, ktloop
         dtlos(kloop,i) = dtlos(kloop,i)  + gloss(kloop,i)
         cnew(kloop,i)  = conc(kloop,i)   + dtlos(kloop,i)
         errymax        = gloss(kloop,i)  * chold(kloop,i)
         dely(kloop)    = dely(kloop)     + errymax * errymax
        enddo
       enddo
      else
       do     i         = 1, ischan
        do     kloop    = 1, ktloop
         dtlos(kloop,i) = dtlos(kloop,i)  + gloss(kloop,i)
         cnew(kloop,i)  = conc(kloop,i)   + asn1  *  dtlos(kloop,i)
         errymax        = gloss(kloop,i)  * chold(kloop,i)
         dely(kloop)    = dely(kloop)     + errymax * errymax
        enddo
       enddo
      endif
!
! *********************************************************************
! * set the previous rms error and calculate the new rms error.       *
! * if dcon < 1, then sufficient convergence has occurred. otherwise, *
! * if the ratio of the current to previous rmserr is decreasing,     *
! * iterate more. if it is not, then the convergence test failed      *
! *********************************************************************
!
      rmserrp           = rmserr
      der2max           = 0.
!
      do     kloop      = 1, ktloop
       if (dely(kloop).gt.der2max) der2max = dely(kloop)   
      enddo
!
      rmserr             = sqrt(der2max / order)
!
      l3                 = l3 + 1
!
      if (l3.gt.1) then
       rmsrat            = rmserr / rmserrp
       drate             = max(0.2 * drate, rmsrat)
      endif
!
      dcon               = rmserr * min(conpst(nqq),conp15(nqq)*drate)
!
! *********************************************************************
!       if convergence occurs, go on to check accumulated error
! *********************************************************************
!
      if (dcon .le. 1.0) then
       goto 390 
!
! *********************************************************************
!  if nonconvergence after one step, re-evaluate first derivative with
!                           new values of cnew
! *********************************************************************
!
!     elseif (l3.lt.mstep.and.(l3.eq.1.or.rmsrat.le.0.9)) then
      elseif (l3.eq.1) then
       goto 270
!
! *********************************************************************
! *             the corrector iteration failed to converge            *
! * if the jacobian matrix is more than one step old, update the      *
! * jacobian and try convergence again. if the jacobian is current,   * 
! * then reduce the time-step, re-set the accumulated derivatives to  *
! * their values before the failed step, and retry with the smaller   *
! * step.                                                             *
! *********************************************************************
!
      elseif (jeval .eq. 0) then
!       ifail           = ifail + 1
       jeval           = 1
       goto 220
      endif
!
!      nfail            = nfail + 1
      rdelmax          = 2.0
      jeval            = 1
      ifsuccess        = 0
      xelaps           = told
      rdelt            = fracdec
!
! *********************************************************************
!               subtract off derivatives previously added 
! *********************************************************************
! this set of operations is equivalent to loop 419.
!
      if (nqq.eq.1) then
       do     i        = 1, ischan
        j              = i + ischan
        do     kloop   = 1, ktloop
         conc(kloop,i) = conc(kloop,i) - conc(kloop,j)
        enddo
       enddo
!
      elseif (nqq.eq.2) then
!
       do     i         = 1, ischan
        j1              = i  + ischan
        j2              = j1 + ischan 
        do     kloop    = 1, ktloop
         conc(kloop, i) = conc(kloop, i) - conc(kloop,j1)-conc(kloop,j2) 
         conc(kloop,j1) = conc(kloop,j1) - conc(kloop,j2) * 2. 
        enddo
       enddo
!
      elseif (nqq.eq.3) then
       do     i         = 1,   ischan
        j1              = i  + ischan
        j2              = j1 + ischan
        j3              = j2 + ischan
        do     kloop    = 1, ktloop
         conc3j3        = conc(kloop,j3) * 3. 
         conc(kloop, i) = conc(kloop, i) - conc(kloop,j1) 
     1                  - conc(kloop,j2) - conc(kloop,j3)
         conc(kloop,j1) = conc(kloop,j1) - conc(kloop,j2) * 2. - conc3j3 
         conc(kloop,j2) = conc(kloop,j2) - conc3j3
        enddo
       enddo
!
      elseif (nqq.eq.4) then
!
       do     i         = 1, ischan
        j1              = i  + ischan
        j2              = j1 + ischan
        j3              = j2 + ischan
        j4              = j3 + ischan
        do     kloop    = 1, ktloop
         conc3j3        = conc(kloop,j3) * 3. 
         conc4j4        = conc(kloop,j4) * 4. 
         conc(kloop, i) = conc(kloop, i) - conc(kloop,j1) 
     1                  - conc(kloop,j2) - conc(kloop,j3)
     2                  - conc(kloop,j4)
         conc(kloop,j1) = conc(kloop,j1) - conc(kloop,j2) * 2. - conc3j3 
     1                                   - conc4j4
         conc(kloop,j2) = conc(kloop,j2) - conc3j3 - conc(kloop,j4) * 6. 
         conc(kloop,j3) = conc(kloop,j3) - conc4j4
        enddo
       enddo
!
      elseif (nqq.eq.5) then
!
       do     i         = 1, ischan
        j1              = i  + ischan
        j2              = j1 + ischan
        j3              = j2 + ischan
        j4              = j3 + ischan
        j5              = j4 + ischan
        do     kloop    = 1, ktloop
         conc3j3        = conc(kloop,j3) * 3. 
         conc4j4        = conc(kloop,j4) * 4. 
         conc5j5        = conc(kloop,j5) * 5. 
         conc10j5       = conc5j5 + conc5j5 
         conc(kloop, i) = conc(kloop, i) - conc(kloop,j1) - 
     1                    conc(kloop,j2) - conc(kloop,j3) - 
     2                    conc(kloop,j4) - conc(kloop,j5)
         conc(kloop,j1) = conc(kloop,j1) - conc(kloop,j2) * 2. - conc3j3 
     1                                   - conc4j4 - conc5j5
         conc(kloop,j2) = conc(kloop,j2) - conc3j3 - conc(kloop,j4) * 6. 
     2                                   - conc10j5
         conc(kloop,j3) = conc(kloop,j3) - conc4j4 - conc10j5
         conc(kloop,j4) = conc(kloop,j4) - conc5j5
        enddo
       enddo
      endif
!
      goto 170
!
! *********************************************************************
! *               the corrector iteration converged                   *  
! * set jeval = 0 so that it does not need to be called the next step *
! * if all else goes well. next, test the accumulated error from the  * 
! *                    convergence process, above                     * 
! *********************************************************************
!
 390  jeval          = 0
!
      if (l3.gt.1) then
       do     kloop  = 1, ktloop
        dely(kloop)  = 0.
       enddo
!
       do     jspc   = 1, ischan     
        do     kloop = 1, ktloop 
         errymax     = dtlos(kloop,jspc) * chold(kloop,jspc)
         dely(kloop) = dely(kloop) + errymax * errymax
        enddo
       enddo
!
       der2max       = 0.
!
       do     kloop  = 1, ktloop
        if (dely(kloop).gt.der2max) der2max = dely(kloop)   
       enddo
      endif
!
! *********************************************************************
! *                 the accumulated error test failed                 *  
! * in all cases, re-set the derivatives to their values before the   *
! * last time-step. next                                              *    
! *                                                                   * 
! * (a)   re-estimate a time-step at the same or one lower order and  *   
! *        retry the step.                                            *  
! * (b)   if the first attempts fail, retry the step at fracdec x     * 
! *        the prior step                                             *  
! * (c)   if this fails, re-set the order to 1 and go back to the     * 
! *        beginning, at order = 1, because errors of the wrong order *
! *        have accumulated                                           *
! *********************************************************************
! 
      if (der2max.gt.enqq) then
       xelaps           = told
!       lfail            = lfail + 1
       jfail            = jfail  + 1
!
       i1               = nqqisc + 1
       do     jb        = 1, nqq
        i1              = i1 - ischan
        do     i        = i1, nqqisc
         j              = i + ischan
         do     kloop   = 1, ktloop
          conc(kloop,i) = conc(kloop,i) - conc(kloop,j)
         enddo
        enddo
       enddo
!
       rdelmax          = 2.0
       if (jfail.le.6) then
        ifsuccess       = 0 
        rdeltup         = 0.0
        goto 540
       elseif (jfail.le.20) then
        ifsuccess       = 0 
        rdelt           = fracdec
        goto 170 
       else
        delt            = delt * 0.1
        rdelt           = 1.
        jfail           = 0 
        jrestar         = jrestar + 1
        idoub           = 5
!
        do    jspc        = 1, ischan
         do    kloop      = 1, ktloop
          cnew(kloop,jspc) = conc(kloop,jspc)
         enddo
        enddo
!
        write(iout,14) delt, xelaps 
        if (jrestar.eq.100) then
         write(iout,15)
         flag = 3
         return
        endif
!
        goto 140
       endif
!
      else
!
! *********************************************************************
! *             all successful steps come through here                * 
! *                                                                   * 
! * after a successful step, update the concentration and all deriv-  *
! * atives, reset told, set ifsuccess = 1, increment nsteps, and      *
! * reset jfail = 0.                                                  *
! *********************************************************************
!
       jfail            = 0
       ifsuccess        = 1
       nsteps           = nsteps + 1
       told             = xelaps
!
! *********************************************************************
!
       i1               = 1 
       do     j         = 2, kstep
        i1              = i1 + ischan
        asnqqj          = aset(nqq,j) 
        do     jspc     = 1, ischan
         i              = jspc + i1 - 1 
         do     kloop   = 1, ktloop
          conc(kloop,i) = conc(kloop,i) + asnqqj * dtlos(kloop,jspc)
         enddo
        enddo
       enddo
!
       if (asn1.eq.1.0) then 
        do     jspc         = 1, ischan
         do     kloop       = 1, ktloop
          conc( kloop,jspc) = conc( kloop,jspc) + dtlos(kloop,jspc)
         enddo
        enddo
       else
        do     jspc         = 1, ischan
         do     kloop       = 1, ktloop
          conc( kloop,jspc) = conc( kloop,jspc) + asn1*dtlos(kloop,jspc)
         enddo
        enddo
       endif 
!
! *********************************************************************
!          exit smvgear if a time interval has been completed
! *********************************************************************
!
       timremain        = tinterval - xelaps
       if (timremain.le.1.0e-06) goto 650
!
! *********************************************************************
! * idoub counts the number of successful steps before re-testing the *
! * step-size and order                                               *
! *                                                                   * 
! * if idoub > 1, decrease idoub and go on to the next time-step with *
! *               the current step-size and order.                    * 
! * if idoub = 1, store the value of the error (dtlos) for the time-  *
! *               step prediction, which will occur when idoub = 0,   *
! *               but go on to the next step with the current step-   *
! *               size and order.                                     * 
! * if idoub = 0, test the time-step and order for a change.          *  
! *********************************************************************
!
       if (idoub.gt.1) then
        idoub               = idoub - 1
        if (idoub.eq.1) then
         do     jspc        = 1, ischan, 2 
          jg1               = jspc + 1
          do     kloop      = 1, ktloop
           cest(kloop,jspc) = dtlos(kloop,jspc)
           cest(kloop,jg1)  = dtlos(kloop,jg1)
          enddo
         enddo
        endif
        rdelt               = 1.0
        goto 170
       endif
!
      endif 
!     endif der2max.gt.enqq
!
! *********************************************************************
! *         test whether to change the step-size and order            * 
! * determine the time-step at (a) one order lower than, (b) the same *
! * order as, and (c) one order higher than the current order. in the *
! * case of multiple grid-cells in a grid-block, find the minimum     *
! * step-size among all the cells for each of the orders. then, in    *
! * all cases, choose the longest time-step among the three steps     *
! * paired with orders, and choose the order allowing this longest    *
! * step.                                                             *
! *********************************************************************
!
! *********************************************************************
! * estimate the time-step ratio (rdeltup) at one order higher than   *
! * the current order. if nqq >= maxord, then we do not allow the     * 
! *                        order to increase.                         *
! *********************************************************************
!
      if (nqq.lt.maxord) then  
       do     kloop  = 1, ktloop
        dely(kloop)  = 0.
       enddo
!
       do     jspc   = 1, ischan     
        do     kloop = 1, ktloop 
         errymax     = (dtlos(kloop,jspc) - cest(kloop,jspc)) * 
     1                  chold(kloop,jspc)
         dely(kloop) = dely(kloop) + errymax * errymax
        enddo
       enddo
!
       der3max       = 0.
!
       do     kloop  = 1, ktloop
        if (dely(kloop).gt.der3max) der3max = dely(kloop)   
       enddo
!       
       rdeltup       = 1.0 / (conp3*der3max**enqq3(nqq)+1.4e-6)
      else
       rdeltup       = 0.0
      endif
!
! *********************************************************************
! *    estimate the time-step ratio (rdeltsm) at the current order    *
! *      we calculated der2max during the error tests earlier         *  
! *********************************************************************
!
 540  rdeltsm        = 1.0 / (conp2*der2max**enqq2(nqq)+1.2e-6)
!
! *********************************************************************
! * estimate the time-step ratio (rdeltdn) at one order lower than    *
! * the current order. if nqq = 1, then we cannot test a lower order. * 
! *********************************************************************
!
      if (nqq.gt.1) then
       do     kloop  = 1, ktloop
        dely(kloop)  = 0.
       enddo
!
       kstepisc      = (kstep - 1) * ischan
       do     jspc   = 1, ischan     
        i            = jspc + kstepisc
        do     kloop = 1, ktloop 
         errymax     = conc(kloop,i) * chold(kloop,jspc)
         dely(kloop) = dely(kloop)   + errymax * errymax
        enddo
       enddo
!
       der1max       = 0.
!
       do     kloop = 1, ktloop
        if (dely(kloop).gt.der1max) der1max = dely(kloop)   
       enddo
        
       rdeltdn       = 1.0 / (conp1*der1max**enqq1(nqq)+1.3e-6)
!
      else
       rdeltdn       = 0.
      endif
!
! *********************************************************************
! * find the largest of the predicted time-steps ratios of each order * 
! *********************************************************************
!
      rdelt        = max(rdeltup,rdeltsm,rdeltdn)
!
! *********************************************************************
! * if the last step was successful and rdelt is small, keep the      *
! * current step and order and allow three successful steps before    *
! * re-checking the time-step and order.                              *
! *********************************************************************
!
      if (rdelt.lt.1.1.and.ifsuccess.eq.1) then
       idoub         = 3
       goto 170
!
! *********************************************************************
! * if the maximum time-step ratio is that of one order lower than    *  
! * the current order, decrease the order. do not minimize rdelt      *
! * to =< 1 when ifsuccess = 0 since this is less efficient.          *
! *********************************************************************
!
      elseif (rdelt.eq.rdeltdn) then
       nqq              = nqq - 1
!
! *********************************************************************
! * if the maximum time-step ratio is that of one order higher than   *  
! * the current order, increase the order and add a derivative term   *
! * for the higher order.                                             *  
! *********************************************************************
!
      elseif (rdelt.eq.rdeltup) then
       consmult         = aset(nqq,kstep) / float(kstep)
       nqq              = kstep
       nqisc            = nqq * ischan
       do     jspc      = 1, ischan, 2     
        jg1             = jspc + 1
        i1              = jspc + nqisc 
        i2              = jg1  + nqisc 
        do     kloop    = 1, ktloop 
         conc(kloop,i1) = dtlos(kloop,jspc) * consmult
         conc(kloop,i2) = dtlos(kloop,jg1)  * consmult
        enddo
       enddo
      endif
!
! *********************************************************************
! * if the last two steps have failed, re-set idoub to the current    *
! * order + 1. do not minimize rdelt if jfail.ge.2 since tests show   *
! * that this merely leads to additional computations.                *
! *********************************************************************
!
      idoub             = nqq + 1
!
      goto 170  
!
! *********************************************************************
! *                      update counters                              *
! *********************************************************************
!
 650   continue 
       flag = 0
!      nsftot             = nsftot   + nsubfun 
!      npdtot             = npdtot   + npderiv
!      nsttot             = nsttot   + nsteps 
!      ifailtot           = ifailtot + ifail
!      nfailtot           = nfailtot + nfail
!      lfailtot           = lfailtot + lfail
!
! *********************************************************************
!                             formats
! *********************************************************************
!
 12   format('smvgear: delt= ',1pe8.2,' too low dec yfac. kblk, ', 
     1       'ktloop, ncs, time, timremain, yfac, ', 
     2       'eps = ',/3(1x,i4),2x,4(1pe9.3,1x))
 13   format('smvgear: too many decreases of yfac ')
 14   format('delt dec to =',e13.5,'; time ',e13.5,' because ',
     1       'excessive errors')                         
 15   format('smvgear: stop because of excessive errors.')   
 16   format('m1,m2,k,err = ',3(i4),2x,1pe10.4)
 17   format('conc when stop = ',2(i4,1x),a14,2(1x,1pe10.2))
!
! *********************************************************************
! ***************     end of subroutine smvgear     *******************
! *********************************************************************
!
      return
      endsubroutine smvgear


!==============================================================

      subroutine jsparse(ncs)
      !Input variables
      integer, intent(in) :: ncs
      !Local variables
      integer :: ncsp, ifsun
      integer :: nrept,i,j,nar,nk,k,ireact,l,ipo,nochang,jold,jnew
      integer :: minvalu,iminold,iminnew,inew,iold,nklast,ial,ire
      integer :: nmo,nol,isdiff,ib,jspcl,ispc1,ispc2,ispc3,iap,iprod
      integer :: ipr,lfrac,ngn,kprods,kdif,npl,ic,nk1,ntwo,icb,icd
      integer :: nkn,igr,isp,nsp,ngr,ngtsum,nltsum,ngsum,nlsum,ngfsum
      integer :: n,jgas,na,ihireac,jal,jre,jpr
      integer :: knumporl,nccount,nremain,nfive,nfour,nthree,none,mc
      integer :: ir,jr,iar,jp,jspc

      real(kind=DP) :: rfrac,alfrac,diff,tnumgna,tnumgn
      real(kind=DP) :: tnumls,sumgn,tsumgna,tnumlsa

      !integer, save :: npltot,nplfun,nfrcoun,npdcoun
     
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!                                                                         
!        JJ  SSSSSSS  PPPPPPP     A      RRRRRRR  SSSSSSS  EEEEEEE
!         J  S        P     P    A A     R     R  S        E 
!         J  SSSSSSS  PPPPPPP   A   A    RRRRRRR  SSSSSSS  EEEEEEE
!   J     J        S  P        AAAAAAA   R  R           S  E
!   JJJJJJJ  SSSSSSS  P       A       A  R    R   SSSSSSS  EEEEEEE 
!
! *********************************************************************
! * this routine sets up sparse-matrix and other arrays for smvgear   *
! * (sparse-matrix vectorized gear-code. it sets arrays for gas-      *
! * -phase, aqueous-phase, and any other type of chemistry. it also   *
! * sets arrays for both day and night chemistry of each type.        *
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call jsparse.f from readchem.f with                              * 
! *     ncs  = 1..ncsgas for gas chemistry                            *
! *********************************************************************
!
! *********************************************************************
! ******* sets up arrays for gas- and aqueous-phase chemistry  ******** 
! * includes arrays for calculating first derivatives, partial deriv- *
! * atives, matrix decompostion, and matrix back-substitution. first, *
! * jsparse re-orders the ordinary differential equations to maximize *
! * the number of zeros in the matrix of partial derivatives. it      *
! * later sets arrays to eliminate all calculations involving a zero. * 
! *********************************************************************
* 
! ntspec    = total number of active + inactive species.
! nspec     = total number of active species.
! nmreac    = 3 = maximum number of active reactants in a reaction 
! nallreac  = 4 = total reactant positions in a reaction 
! nmprod    = 5 = maximun number of active products in a reaction 
! nprodlo   = nallreac  + 1 = lowest product position number. 
! nprodhi   = nallreac + nmprod = highest product position number. 
!
! *********************************************************************
! * determine how many partial deriv terms are needed for each species*
! *********************************************************************
! ifrepro   = 1 then species is lost and reproduced in reaction nk 
! irm       = species # of each react or product in each nk reaction
! isaporl   = counts partial derivative terms for each species
! fkoef     = 1, 2, fraction, or more = # of a given reactant or products
!             e.g. reaction      a + b  --> 2c + 0.34d  + d 
!             value of fkoef     1   1      2    0.34     1     
! ncs       = 1..ncsgas for gas chemistry                            
! ncsp      = ncs        for daytime   gas chem 
!           = ncs   +ics for nighttime gas chem            
! nk        = reaction # of each reaction 
! nrates    = number of kinetic (non-photo) rate coefficients
! ntrates   = number of kinetic plus photo  rate coefficients
! nallrat   = number of kinetic plus photo reaction rates  
!
       
       ncsp                        = ncs + ics
       nrept                       = 0
!
       do    i                     = 1, mxgsaer 
        isaporl( i)                = 0
       enddo
       
!
       do    i                     = 1, maxgl 
        newnk(i)                   = 0
       enddo
!
       do    i                     = 1, mxgsaer
        do    j                    = 1, mxgsaer
         isparder(i,j)             = 0
        enddo
       enddo
!
       
       do    nar                 = 1, nallrat(ncs)
        nk                        = ncequat(nar,ncs) 
        
        if (nk.le.nrates(ncs))      nallrat(ncsp) = nar
        do    k                   = 1, nmreac  
         ireact                   = irm(k,nk,ncs)
         if (ireact.gt.0.and.ireact.le.nspec(ncs)) then
          do    l                 = 1, nprodhi  
           ipo                    = irm(l,nk,ncs)
           if ((l.le.nmreac.or.l.ge.nprodlo).and.ipo.gt.0.and.
     1          ipo.le.nspec(ncs)) isparder(ipo,ireact) = 1 
          enddo!l
         endif 
        enddo!k
       enddo!nar
!      continue nar = 1, nallrat
!
       do    ireact                = 1, ntspec(ncs)
        do    ipo                  = 1, ntspec(ncs)
         if (isparder(ipo,ireact).eq.1) isaporl(ipo)=isaporl(ipo)+1
        enddo
       enddo
!
! *********************************************************************
! *  re-arrage species array so that all species with at least one    * 
! *  partial derivative term appear first, and those with zero        *
! *  appear last.                                                     * 
! *********************************************************************
! ischang = number of original nspec species with at least one pd term. 
! inewold = original species number of each new jnew species 
! mappl   = new species number for chemistry of each original jold species 
!
       nochang                 = nspec(ncs) 
       do     jold             = 1, ntspec(ncs)
        if (jold.gt.nspec(ncs)) then 
         mappl(jold,ncs)       = jold
         inewold(jold,ncs)     = jold   
        elseif (isaporl(jold).gt.0) then
         ischang(ncs)          = ischang(ncs) + 1
         jnew                  = ischang(ncs) 
         inewold(jnew,ncs)     = jold   
         mappl(jold,ncs)       = jnew  
        else
         inewold(nochang,ncs)  = jold   
         mappl(jold,ncs)       = nochang  
         nochang               = nochang - 1
        endif
       enddo!jold
!
! *********************************************************************
! *  re-arrage species in ischang array so that species with the      *
! *  fewest partial derivative terms combined are placed first,       *
! *  and those with the most appear last. however, species with zero  *
! *  partial derivative terms still appear after all ischang species  *
! *********************************************************************
! 
       do     jnew             = 1, ischang(ncs)
        jold                   = inewold(jnew,ncs)
        minvalu                = isaporl(jold)
        iminold                = jold 
        iminnew                = jnew
        do     inew            = jnew+1, ischang(ncs)
         iold                  = inewold(inew,ncs)
         if (isaporl(iold).lt.minvalu) then
          minvalu              = isaporl(iold)
          iminold              = iold  
          iminnew              = inew
         endif
        enddo!inew
        inewold(iminnew,ncs)   = jold  
        inewold(jnew,ncs)      = iminold  
        mappl(jold,ncs)        = iminnew    
        mappl(iminold,ncs)     = jnew    
       enddo!jnew
!
! *********************************************************************
! *                    count gross and net loss                       * 
! *********************************************************************
! ioner    = number of reactions with one active reactant
! itwor    = number of reactions with two active reactants 
! ithrr    = number of reactions with three active reactants 
! nkoner   = reaction number of each ioner reaction 
! nktwor   = reaction number of each itwor reaction 
! nkthrr   = reaction number of each ithrr reaction 
! numlost  = every occurence of a loss (active & inactive spec) 
! numloss  = every net occurence of a loss where the species is not 
!            reproduced in the same reaction. (active & inactive species)
! jloss    = reaction number of each net loss occurrence
! irm2     = identifies each new active species number in each reaction
! numkial  = number of reactions with either 1, 2, or 3 active reactants  
! nksdt    = reaction number of each numkial reaction 
! nruse    = 1,2,3 if reaction has 1, 2, or 3 active reactants, respectively.
! nrrep    = 0 for each of two reactions where the reactants are identical.
!            if more than two reactions have the same reactants, nrrep = 0
!            for the first two reactions only.
!          = 1,2,3 if reaction has 1, 2, or 3 reactants, respectively.
! nmoth    = # of occurrences where inactive spec appears in rate equation
!            excludes third bodies in array nm3bod (e.g., o2, n2, m, when
!            these species do not lose concentration in the reaction)
! nreacoth = reaction number of each nmoth occurrence
! lgasbino = old species number of each inactive species  
!
       nolosp(ncsp)              = 0
       nklast                    = 0 
!
       do     nar                = 1, nallrat(ncs)
        nk                       = ncequat(nar,ncs)
!
! *********************************************************************
! ***  determine occurrences of inactive species in rate equations  ***
! *             set array to identify active loss species             *
! *********************************************************************
!
        ial                      = 0
!
        do     jspc              = 1, mxgsaer
         aporl(jspc)             = 0.
        enddo
!
        do     j                 = 1, nmreac  
         ireact                  = irm(j,nk,ncs)
         if (ireact.gt.0) then
          ire                    = mappl(ireact,ncs)
!
          aporl(ire)             = aporl(ire) - 1.
          numlost(ire,ncs)       = numlost(ire,ncs) + 1
!
          if (ire.le.nspec(ncs)) then
!
           ial                   = ial + 1
           irm2(ial,nk,ncs)      = ire
!
          elseif (ire.gt.nspec(ncs)) then
!
           if (nk.le.nrates(ncs)) then
            nmoth(ncs)           = nmoth(ncs) + 1
            nmo                  = nmoth(ncs)
            nreacoth(nmo,ncs)    = nk
            lgasbino(nmo,ncs)    = ireact
           else
            nolosp(ncs)          = nolosp(ncs) + 1
            nol                  = nolosp(ncs)
            nknlosp(nol,ncs)     = nk
            losinacp(nol,ncs)    = ireact
           endif
!
          endif
         endif        
!
        enddo
!
! *********************************************************************
! *  set arrays to identify reactions with at least one active loss   *
! *********************************************************************
!
        if (ial.gt.0) then
         nruse(nk,ncs)      = ial 
         nrrep(nk,ncs)      = ial
!
         if (ial.eq.1) then
          ioner(ncs)             = ioner(ncs) + 1
          nkoner(ioner(ncs),ncs) = nk 
         elseif (ial.eq.2) then
          itwor(ncs)             = itwor(ncs) + 1
          nktwor(itwor(ncs),ncs) = nk 
         elseif (ial.eq.3) then 
          ithrr(ncs)             = ithrr(ncs) + 1
          nkthrr(ithrr(ncs),ncs) = nk 
         endif
!
! *********************************************************************
! * compare two consecutive reactions. if the species (but not rates) * 
! * are the same, then save multiplications in subfun.f               *
! *********************************************************************
!
         if (nklast.gt.0) then 
          if (nruse(nklast,ncs).eq.ial) then  
           isdiff           = 0 
           do     ib        = 1, ial  
            jspcl           = irm2(ib,nklast,ncs) 
            jspc            = irm2(ib,nk    ,ncs) 
            if (jspcl.ne.jspc) isdiff = 1 
           enddo
           if (isdiff.eq.0.and.nrrep(nklast,ncs).ne.0) then 
            nrrep(nk,ncs)     = 0
            nrrep(nklast,ncs) = 0
            nrept           = nrept + 1
            ispc1           = irm2(1,nk,ncs) 
            ispc2           = irm2(2,nk,ncs) 
            ispc3           = irm2(3,nk,ncs)
            if (ispc1.gt.0) ispc1 = inewold(ispc1,ncs)
            if (ispc2.gt.0) ispc2 = inewold(ispc2,ncs)
            if (ispc3.gt.0) ispc3 = inewold(ispc3,ncs)
            !print*,iout
            write(iout,1) nrept, nk,namencs(ispc1,ncs), 
     +                      namencs(ispc2,ncs), namencs(ispc3,ncs)
  1         format('repeat reactants: ',i5,i5,3(1x,a14))
           endif 
          endif 
         endif 
!
! *********************************************************************
! *   determine the number of reactions with zero active loss terms   * 
! *********************************************************************
! nolosrat = number of active reactions with no loss terms 
! nolosrn  = reaction number of each reaction with no loss terms
!
        elseif (ial.eq.0) then
         nolosrat(ncs)         = nolosrat(ncs) + 1
         nol                   = nolosrat(ncs)
         nolosrn(nol,ncs)      = nk
        endif
!       endif ial.gt.0 
!
! *********************************************************************
! * count gross and net production and set a partial derivative array * 
! *********************************************************************
! numgaint = every occurence of a production (active & inactive spec) 
! numgain  = every net occurence of a production where the species is 
!            not lost in the same reaction. (active & inactive spec)
! iaprod   = number of active products in each nk reaction. used
!            to calculate partial derivatives in pderiv.f. 
! irm2     = new species # of each active product in each nk reaction
! ngnfrac  = number of fractionated production terms in entire reaction set
! ignfrac  = identifies product of each fractionaed production term
! nkgnfrac = identifies reaction number of each fractionated production term
! fracp    = actual fraction of each fractionated production term
!
        iap                       = nprodlo - 1
        do     k                  = nprodlo, nprodhi  
         iprod                    = irm(k,nk,ncs)
         if (iprod.gt.0) then
          ipr                     = mappl(iprod,ncs)
          rfrac                   = fkoef(k,nk,ncs)
          lfrac                   = int(rfrac + smal1) 
          alfrac                  = float(lfrac)
          diff                    = abs(rfrac-alfrac)
!
! ******************** production term is a fraction ******************
!
          if (diff.gt.smal1) then 
           if (ipr.le.nspec(ncs)) then 
            ngnfrac(ncs)          = ngnfrac(ncs) + 1 
            ngn                   = ngnfrac(ncs) 
            ignfrac( ngn,ncs)     = ipr 
            nkgnfrac(ngn,ncs)     = nk 
            fracp(   ngn,ncs)     = rfrac  
           endif 
           kprods                 = 1
           numgfrt( ipr,ncs)      = numgfrt( ipr,ncs) + 1
           fracgain(ipr,ncs)      = fracgain(ipr,ncs) + rfrac 
!
! ******************* production term is non-fraction *****************
!
          else
           aporl(ipr)             = aporl(ipr) + rfrac
           kprods                 = lfrac
           numgaint(ipr,ncs)      = numgaint(ipr,ncs) + lfrac
           fkoef(k,nk,ncs)        = 1.
          endif  
!
! ******************* identify all production terms *******************
!
          if (ipr.le.nspec(ncs)) then
           do     l               = 1, kprods
            iap                   = iap + 1
            iaprod(nk,ncs)        = iap
            irm2(iap,nk,ncs)      = ipr
            fk2( iap,nk,ncs)      = fkoef(k,nk,ncs)
           enddo!l
          endif
!
         endif
!
        enddo!k
!
! *********************************************************************
! *  find net prod and loss terms for all but fractionated products   * 
! *********************************************************************
!
         do     jspc              = 1, ntspec(ncs)
          if (abs(aporl(jspc)).lt.smal1) then
           kdif                   = 0 
!
          elseif (aporl(jspc).gt.0.) then 
           kdif                   = int(aporl(jspc) + 0.00001)
           do     l               = 1, kdif 
            numgain(jspc,ncs)     = numgain(jspc,ncs) + 1
            numporl(jspc,ncs)     = numporl(jspc,ncs) + 1
            npl                   = numporl(jspc,ncs)
            jporl(jspc,npl,ncs)   = nk + ntrates(ncs)  
           enddo
          else 
           kdif                   = -int(aporl(jspc) - 0.00001)
           do     l               = 1, kdif  
            numloss(jspc,ncs)     = numloss(jspc,ncs) + 1
            numporl(jspc,ncs)     = numporl(jspc,ncs) + 1
            npl                   = numporl(jspc,ncs)
            jporl(jspc,npl,ncs)   = nk 
           enddo
          endif 
!
          if (nk.le.nrates(ncs)) then 
           numloss(jspc,ncsp)     = numloss(jspc,ncs)
           numgain(jspc,ncsp)     = numgain(jspc,ncs)
           numporl(jspc,ncsp)     = numporl(jspc,ncs)
          endif 
!
         enddo!jspc
!
         if (nk.le.nrates(ncs)) then 
          nolosrat(ncsp)          = nolosrat(ncs) 
          ngnfrac( ncsp)          = ngnfrac( ncs)
          ioner(   ncsp)          = ioner(   ncs)
         endif 
!
         nklast                   = nk 
!
       enddo!n
!      continue n = 1, ntrates
!
! *********************************************************************
! * set array for reordering rates from 3..2..1..0 body reactions     *
! *********************************************************************
! inorep   = last reordered reaction number prior to sets of two
!            reactions with two reactants  
! noldfnew = old reaction rate # corresp. to each reordered reaction
! newfold  = new reaction rate # corresp. to each original rate number
!
       ic                 = 0
       do     i           = 1, ithrr(ncs)   
        ic                = ic + 1
        nk                = nkthrr(i,ncs)
        nk1               = nk + ntrates(ncs) 
        noldfnew(ic, ncs) = nk
        newfold( nk, ncs) = ic
        newfold( nk1,ncs) = ic + nallrat(ncs) 
       enddo!i
!
       ntwo                = ithrr(ncs) + itwor(ncs) 
       icb                 = ntwo + 1 
       do     i            = 1, itwor(ncs)   
        nk                 = nktwor(i,ncs)
        nk1                = nk + ntrates(ncs) 
        if (nrrep(nk,ncs).gt.0) then  
         ic                = ic + 1
         icd               = ic
        else 
         icb               = icb - 1
         icd               = icb
        endif 
        noldfnew(icd, ncs) = nk
        newfold( nk,  ncs) = icd 
        newfold( nk1, ncs) = icd + nallrat(ncs) 
       enddo!i
!
       inorep(ncs)         = ic 
       ic                  = ntwo 
       do     i            = 1, ioner(ncs)   
        ic                 = ic + 1
        nk                 = nkoner(i,ncs)
        nk1                = nk + ntrates(ncs) 
        noldfnew(ic, ncs)  = nk
        newfold( nk, ncs)  = ic
        newfold( nk1,ncs)  = ic + nallrat(ncs) 
       enddo!i
!
       do     i            = 1, nolosrat(ncs)   
        ic                 = ic + 1
        nk                 = nolosrn(i,ncs)
        nk1                = nk + ntrates(ncs) 
        noldfnew(ic, ncs)  = nk
        newfold( nk, ncs)  = ic
        newfold( nk1,ncs)  = ic + nallrat(ncs) 
       enddo!i
!
       if (ic.ne.nallrat(ncs)) then 
        write(iout,2) ic, nallrat(ncs)
        stop 
       endif 
!
! *********************************************************************
!                set a slightly more efficient photo array 
! *********************************************************************
!
       do     j          = 1, jphotrat(ncs)
        nk               = nkphotrat(j,ncs)  
        nkn              = newfold(nk,ncs)
        nknphotrt(j,ncs) = nkn 
       enddo!j
!
  2    format('jsparse: ic ne nallrat =',2(i5))
!
! *********************************************************************
! ****** determine number of species with gross/net losses/gains ******
! *********************************************************************
! nspcsolv = # of active species with at least one gross loss
! isolvspc = species number of each nspcsolv species
! isgainr  = # of active species with at least one net chem gain 
! igainr   = species number of each isgainr species
! isgaine  = # of active species with at least 1 net chem gain 
! igainr   = species number of each isgainr species
! nogaine  = # of active species with zero net chem or gains 
! ngaine   = species number of each nogaine species
! isporl   = # of active species with at least one net production
!            or loss term for smvgear.
! iporl    = species number of each isporl species
!
       do     jold             = 1, nspec(ncs) 
        jnew                   = mappl(jold,ncs)
!
        if (numgain(jnew,ncs).gt.0) then
         isgainr(ncs)          = isgainr(ncs) + 1
         igr                   = isgainr(ncs)
         igainr(igr,ncs)       = jnew  
        endif 
!
        if (numporl(jnew,ncs).gt.0) then
         isporl(ncs)           = isporl(ncs) + 1 
         isp                   = isporl(ncs)
         iporl(isp,ncs)        = jnew 
        endif
!
        if (numlost(jnew,ncs).gt.0) then
         nspcsolv(ncs)         = nspcsolv(ncs) + 1
         nsp                   = nspcsolv(ncs)
         isolvspc(nsp,ncs)     = jnew  
        endif
!
        if (numgain(jnew,ncs).gt.0.or.fracgain(jnew,ncs).gt.0) then
         isgaine(ncs)         = isgaine(ncs) + 1
         igr                  = isgaine(ncs)
         igaine(igr,ncs)      = jnew  
        elseif (numloss(jnew,ncs).gt.0) then 
         nogaine(ncs)         = nogaine(ncs) + 1
         ngr                  = nogaine(ncs)
         ngaine(ngr,ncs)      = jnew    
        endif
!
       enddo!jold
!
! *********************************************************************
! ********  check dimensions resulting from gains and losses  *********
! *********************************************************************
!
       ngtsum   = 0
       nltsum   = 0
       ngsum    = 0
       nlsum    = 0
       ngfsum   = 0
       do     k = 1, ntspec(ncs)
        j       = inewold(k,ncs) 
        ngtsum  = ngtsum + numgaint(k,ncs) 
        nltsum  = nltsum + numlost( k,ncs) 
        ngsum   = ngsum  + numgain( k,ncs) 
        nlsum   = nlsum  + numloss( k,ncs) 
        ngfsum  = ngfsum + numgfrt( k,ncs) 
        if (numgaint(k,ncs)   .gt.   maxgl .or.
     1      numlost( k,ncs)   .gt.   maxgl) then 
         write(iout,6) namencs(j,ncs),numgaint(k,ncs),maxgl,
     1                                  numlost( k,ncs),maxgl
         stop
        endif
       enddo!k
!
       if (ioreac.eq.1) then
        write(iout,*)
        write(iout,4)
        do     k = 1, ntspec(ncs) 
         j       = inewold(k,ncs) 
         write(iout,5) 
     1         namencs( j,ncs),numgaint(k,ncs),numgain( k,ncs),
     1         numlost( k,ncs),numloss( k,ncs),numgaint(k,ncs) 
     1        -numlost( k,ncs)-numgain( k,ncs)+numloss( k,ncs),
     1         fracgain(k,ncs),numgfrt( k,ncs)
        enddo!k
        write(iout,5) 'overall       ',ngtsum, ngsum, nltsum, nlsum,
     1                   ngtsum - nltsum - ngsum + nlsum, 0., ngfsum
       endif
!
       if (nmoth(  ncs).gt.maxgl2.or.nolosp(ncs).gt.maxgl3.or.
     1     ngnfrac(ncs).gt.maxgl5) then
        write(iout,7) maxgl2, nmoth(  ncs), maxgl3, nolosp(ncs),
     1                  maxgl5, ngnfrac(ncs)
        stop
       endif
!
! *********************************************************************
! *       check whether chemical system is atom-conservative          *
! *********************************************************************
! jmbcomp = species number for each species in a mass bal. group
! mbcomp  = counts the number of mass balance species in each m.b group
! nmasbal = number of mass balance groups (e.g. s, n, c are groups)
! wtmb(1) = number of atoms of a given mass balance species per molecule 
!
       write(iout,3) chemtyp(ncs)
!
!       if (ncs.le.ncsgas) then
! 
! ----------------------------   gas-phase   -------------------------- 
!
!        do 385 n     = 1, nmasbal 
!         if (mbcomp(n,mb1).gt.0) then 
!          tnumgn     = 0
!          tnumls     = 0
!          write(iout,325) namemb(n)
!          do 380 j   = 1, mbcomp(n,mb1) 
!           jgas      = jmbcomp(n,j,mb1)
!           jnew      = mappl(jgas,ncs)
!           sumgn     = numgain(jnew,ncs) + fracgain(jnew,ncs) 
!           tnumgna   = sumgn             * wtmb(n,jgas,mb1)  
!           tnumlsa   = numloss(jnew,ncs) * wtmb(n,jgas,mb1) 
!          tnumgn    = tnumgn + tnumgna
!           tnumls    = tnumls + tnumlsa
!           write(iout,320) namencs(jgas,ncs), tnumgna, tnumlsa, 0 
! 380      continue
!          write(iout,370) tnumgn, tnumls, tnumgn - tnumls 
!         endif 
! 385    continue
!       endif
!
!
  3    format(/'change in moles due to ',a14,' chemistry')
!325   format('mass balance group              = ',a14)
!320   format('gains/losses for ',a14,' = ',2(f8.3),i5)     
!370   format('total gains - losses            = ',3(f8.3)) 
!375   format(/'# kinetic reactions: ',i5,' photorates: ',i5,
!    1        ' total: ',i5) 
  4    format('spec            numgt   numg   numlt  numl    ngt-nlt- ',
     1        'ng+nl fracgn numgft')
  5    format(a14,4(2x,i5),7x,i5,3x,f10.3,i5)
  6    format('jsparse: spec ',a6,' dimens exceeded. either numgaint ', 
     1        'numloss,numgain, or numlost > maxgl ',
     2        4(i5,1x)) 
  7    format('one of the dimensions below is too small:',/,
     1        'dimension: maxgl2   =  ',i4,' variable: nmoth    = ',i4/  
     2        'dimension: maxgl3   =  ',i4,' variable: nolosp   = ',i4/
     3        'dimension: maxgl5   =  ',i4,' variable: ngnfrac  = ',i4)  
!
! *********************************************************************
! *********************************************************************
! **        set arrays to take advantage of sparse matrices          ** 
! *********************************************************************
! *********************************************************************
!
! ifsun  = 1 then day-chemistry;  = 2 then night chemistry
! ncsp   = ncs       for daytime   trop-gas, strat-gas chem  
! ncsp   = ncs + icp for nighttime trop-gas, strat-gas chem  
!
! lzero    = 1 if an array spot is filled with a non-zero value. lzero
!            is updated as we simulate the order of calculations during
!            a practice l-u decomposition
! mxgsaer  = larger of igas, iaerty
!
!
      if (ifnone.eq.0) then
       ifnone                 = 1
       nplfun                 = 0 
       nfrcoun                = 0 
       npdcoun                = 0 
       npltot                 = 0 
      endif
!
      do     ifsun            = 1, 2 
       ncsp                   = (ifsun - 1) * ics + ncs
!
       do     i               = 1, mxgsaer
        do     j              = 1, mxgsaer
         lzero(j,i)           = 0
        enddo!j
        lzero(i,i)            = 1
       enddo!i
!
       do     na              = 1, nallrat(ncsp)
        nk                    = ncequat(na,ncs)
        ihireac               = nruse(  nk,ncs)
        do     ial            = 1, ihireac
         ire                  = irm2(ial,nk,ncs)
         do    jal            = 1, ihireac 
          jre                 = irm2(jal,nk,ncs)
          lzero(jre,ire)      = 1
         enddo!jal
         do     iap           = nprodlo, iaprod(nk,ncs)
          jpr                 = irm2(iap,nk,ncs)
          lzero(jpr,ire)      = 1
         enddo!iap
        enddo!ial
       enddo!na
!
! *********************************************************************
! *   set decomposition and back-substitution sparse-matrix arrays    *
! *********************************************************************
!
       call ksparse(ncs, ncsp)
!
! *********************************************************************
! *    set arrays to improve efficiency of first-derivative calcs     * 
! *********************************************************************
! *********************************************************************
! **   set arrays for kinetic and photo production and loss rates    **
! *********************************************************************
!
       npllo(ncsp)         = npltot + 1
       do     i            = 1, isporl(ncs)
        jspc               = iporl(i,ncs)
        knumporl           = numporl(jspc,ncsp) 
        nccount            = 0 
        npltot             = npltot + 1
        nremain            = knumporl
        nfive              = (nremain + 0.0001) / 5 
        nremain            =  nremain - nfive   * 5 
        nfour              = (nremain + 0.0001) / 4 
        nremain            =  nremain - nfour   * 4
        nthree             = (nremain + 0.0001) / 3  
        nremain            =  nremain - nthree  * 3 
        ntwo               = (nremain + 0.0001) / 2   
        nremain            =  nremain - ntwo    * 2  
        none               = (nremain + 0.0001)  
        nremain            =  nremain - none
!
        jspnpl(npltot)     = jspc 
        npl5(  npltot)     = nplfun       + 1
        nph5(  npltot)     = nplfun       + nfive  
        npl4(  npltot)     = nph5(npltot) + 1
        nph4(  npltot)     = nph5(npltot) + nfour   
        npl3(  npltot)     = nph4(npltot) + 1
        nph3(  npltot)     = nph4(npltot) + nthree
        npl2(  npltot)     = nph3(npltot) + 1
        nph2(  npltot)     = nph3(npltot) + ntwo
        npl1(  npltot)     = nph2(npltot) + 1
        nph1(  npltot)     = nph2(npltot) + none
        nplfun             = nph1(npltot)
!
        do     n           = 1, knumporl 
         nk                = jporl(jspc,n,ncs) 
         newnk(n)          = newfold(nk,ncs)  
        enddo!n
!
        do     mc          = npl5(npltot), nph5(npltot)
         lossra(mc)        = newnk(nccount+1) 
         lossrb(mc)        = newnk(nccount+2) 
         lossrc(mc)        = newnk(nccount+3) 
         lossrd(mc)        = newnk(nccount+4) 
         lossre(mc)        = newnk(nccount+5) 
         nccount           = nccount + 5
        enddo
!
        do mc          = npl4(npltot), nph4(npltot)
         lossra(mc)        = newnk(nccount+1) 
         lossrb(mc)        = newnk(nccount+2) 
         lossrc(mc)        = newnk(nccount+3) 
         lossrd(mc)        = newnk(nccount+4) 
         nccount           = nccount + 4  
        enddo!mc
!
        do     mc          = npl3(npltot), nph3(npltot)
         lossra(mc)        = newnk(nccount+1) 
         lossrb(mc)        = newnk(nccount+2) 
         lossrc(mc)        = newnk(nccount+3) 
         nccount           = nccount + 3  
        enddo!mc
!
        do     mc          = npl2(npltot), nph2(npltot)
         lossra(mc)        = newnk(nccount+1) 
         lossrb(mc)        = newnk(nccount+2) 
         nccount           = nccount + 2    
        enddo!mc
!
        do     mc          = npl1(npltot), nph1(npltot)
         lossra(mc)        = newnk(nccount+1) 
         nccount           = nccount + 1     
        enddo!mc
!
       enddo!i
       nplhi(ncsp)         = npltot
!
! *********************************************************************
! *              set array for fractionated products                  *  
! *********************************************************************
!
       nfrlo(ncsp)          = nfrcoun + 1 
       do     i             = 1, ngnfrac(ncsp)
        jspc                = ignfrac(i,ncs)
        nfrcoun             = nfrcoun + 1 
        jspcnfr(nfrcoun)    = jspc 
        nk                  = nkgnfrac(i,ncs)  
        nknfr(  nfrcoun)    = newfold(nk,ncs)
        fracnfr(nfrcoun)    = fracp(i,ncs)
       enddo!i 
       nfrhi(ncsp)          = nfrcoun
!
! *********************************************************************
! * set arrays to improve efficiency of partial derivative calcs      * 
! *********************************************************************
!
       npdlo(ncsp)           = npdcoun + 1
!
       do     na             = 1, nallrat(ncsp) 
        nk                   = ncequat(na,ncs) 
        ihireac              = nruse(  nk,ncs) 
!
        do     ial           = 1, ihireac
         ir                  = irm2(ial,nk,ncs)
         do     jal          = 1, ihireac 
          jr                 = irm2(jal,nk,ncs)
          iar                = jarraypt(jr,ir)
          npdcoun            = npdcoun + 1 
          nkpdterm(npdcoun)  = newfold(nk,ncs)  
          ipospd(  npdcoun)  = iar 
          iialpd(  npdcoun)  = ial  
          fracpl(  npdcoun)  = -1.
         enddo!jal
!
         do     iap          = nprodlo, iaprod(nk,ncs)
          jp                 = irm2(iap,nk,ncs)
          iar                = jarraypt(jp,ir)
          npdcoun            = npdcoun + 1 
          nkpdterm(npdcoun)  = newfold(nk,ncs)  
          ipospd(  npdcoun)  = iar 
          iialpd(  npdcoun)  = ial  
          fracpl(  npdcoun)  = fk2(iap,nk,ncs)  
         enddo!iap
        enddo!ial
       enddo!na
!
       npdhi(ncsp)          = npdcoun
!
! *********************************************************************
! **        check dimensions and print out array savings             ** 
! *********************************************************************
!
       if (npltot  .gt. mxcount4  .or. nplfun   .gt. mxcount4 .or.
     1     nfrcoun .gt. mxcount4 .or.  npdcoun  .gt. mxcount2) then
        write(iout,8) mxcount4, npltot,    mxcount4, nplfun,
     1                  mxcount4, nfrcoun,   mxcount2, npdcoun
        stop
       endif
!
      enddo!ifsun
!     continue ifsun = 1, 2
!
        
  8   format('one of the dimensions below is too small:',/,
     1       'dimension: mxcount4 =  ',i5,' variable: npltot   = ',i5,/,
     2       'dimension: mxcount4 =  ',i5,' variable: nplfun   = ',i5,/,
     3       'dimension: mxcount4 =  ',i5,' variable: nfrcoun  = ',i5,/,
     4       'dimension: mxcount2 =  ',i5,' variable: npdcoun  = ',i5)
!
! *********************************************************************
! ********************** end of subroutine jsparse ********************
! *********************************************************************
!
      return                                                             
      endsubroutine jsparse

!============================================================

      subroutine ksparse(ncs, ncsp)
      !input variables
      integer, intent(in) :: ncs, ncsp
      
      !local variables
      integer :: kount0a,kount0,icnta,icntb
      integer :: kcnta,kcntb,mcnta,mcntv,iarray2,j,k,j1,i,i1,i2,kntarray
      integer :: izil,nremain,nfive,nfour,nthree,ntwo,none,ic,ka,kb,kc,kd
      integer :: ia,kzil,mc,jcnta,jcntb,mcntb,ke,mzil
      real(kind=DP), dimension(MORDER,3) :: pertdat

      !integer, save :: mcnt,kcnt,icnt,jcnt,mztot,ijtot,kztot,idecomp
      !integer, save :: mccount,iccount,jccount,kccount,kbsub,mbsub
      
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!                                                                         
!   K    K  SSSSSSS  PPPPPPP     A      RRRRRRR  SSSSSSS  EEEEEEE
!   K  K    S        P     P    A A     R     R  S        E 
!   KK      SSSSSSS  PPPPPPP   A   A    RRRRRRR  SSSSSSS  EEEEEEE
!   K  K          S  P        AAAAAAA   R  R           S  E
!   K    K  SSSSSSS  P       A       A  R    R   SSSSSSS  EEEEEEE 
!
! *********************************************************************
! * this routine sets up sparse-matrix and other arrays. it also      *
! * sets arrays for gas-phase, aqueous-phase, or any other type       *
! * of chemistry. finally, it sets arrays for both day and night      *
! * chemistry of each type.                                           *
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call ksparse.f from jsparse.f with                               * 
! *     ncs  = 1..ncsgas for gas chemistry                            *
! *********************************************************************
!
! *********************************************************************
! * sets up arrays for decomposition / back-substitution of sparse    *
! * matrices by removing all calculations involving a zero.           *
! *********************************************************************
!
      
!
! *********************************************************************
! *********************************************************************
! **        set arrays to take advantage of sparse matrices          ** 
! *********************************************************************
! *********************************************************************
!
! ifsun    = 1 then day-chemistry;  = 2 then night chemistry
! ncsp     = ncs        for daytime   gas chem  
! ncsp     = ncs   +ics for nighttime gas chem  
!
! kount0a  = # initial matrix spots filled w/o  sparse-matrix reductions 
! kount0   = # initial matrix spots filled with sparse-matrix reductions 
! kntarray = # final matrix spots filled w/o  sparse-matrix reductions 
! iarray2  = # final matrix spots filled with sparse-matrix reductions 
! icnta    = # operations in decomp loop 1 w/o  sparse-matrix reductions 
! icntb    = # operations in decomp loop 1 with sparse-matrix reductions 
! jcnta    = # operations in decomp loop 2 w/o  sparse-matrix reductions 
! jcntb    = # operations in decomp loop 2 with sparse-matrix reductions 
! kcnta    = # operations in back-sup loop 1 w/o  sparse-matrix reductions 
! kcntb    = # operations in back-sub loop 1 with sparse-matrix reductions 
! mcnta    = # operations in back-sup loop 2 w/o  sparse-matrix reductions 
! mcntb    = # operations in back-sub loop 2 with sparse-matrix reductions 
!
! lzero    = 1 if an array spot is filled with a non-zero value. lzero
!            is updated as we simulate the order of calculations during
!            a practice l-u decomposition
!
       if (ifnever.eq.0) then
        ifnever               = 1
        icnt                  = 0 
        jcnt                  = 0 
        kcnt                  = 0 
        mcnt                  = 0 
        iccount               = 0
        jccount               = 0
        kccount               = 0 
        mccount               = 0 
        idecomp               = 0
        kbsub                 = 0
        mbsub                 = 0
        ijtot                 = 0 
        kztot                 = 0 
        mztot                 = 0 
       endif
!
       kount0a                = 0
       kount0                 = 0
       icnta                  = 0
       icntb                  = 0
       jcnta                  = 0
       jcntb                  = 0
       kcnta                  = 0
       kcntb                  = 0
       mcnta                  = 0
       mcntb                  = 0
       iarray2                = 0 
!
       do     j               = 1, ischang(ncs)
        do     k              = 1, ischang(ncs)
         kount0a              = kount0a + 1
         if (lzero(k,j).eq.1) kount0 = kount0 + 1
         jarraypt(k,j)        = 0
        enddo
       enddo
!
! *********************************************************************
! **                arrays for decomposition (ludcmp)                ** 
! *********************************************************************
! izilch = # of calculations with non-zero values during matrix decomp
! izero  = each occurrence of each izilch calculation
!
       do     j                  = 1, ischang(ncs)
        jzilch(j)                = 0 
        j1                       = j - 1
!
! ------------------- first loop of decompostion ----------------------
!
        do     i            = 2, ischang(ncs) 
         izilch(j,i)        = 0 
         i1                 = j1 
         if (i.le.j1) i1    = i - 1
         do     k           = 1, i1
          icnta             = icnta + 1
          if (lzero(i,k).eq.1.and.lzero(k,j).eq.1) then
           izilch(j,i)      = izilch(j,i) + 1
           icnt             = icnt  + 1
           icntb            = icntb + 1
           izerok(icnt)     = k   
           lzero(i,j)       = 1 
          endif
         enddo!k
        enddo!i
!
! ------------------- second loop of decompostion ---------------------
!
! jzilch  = # of calculations with non-zero values to fill lower
!           part of decomposed matrix   
!
        do     i            = j+1, ischang(ncs) 
         jcnta              = jcnta + 1
         if (lzero(i,j).eq.1) then
          jzilch(j)         = jzilch(j) + 1
          jcnt              = jcnt  + 1
          jcntb             = jcntb + 1
          jzero(jcnt)       = i  
         endif 
        enddo!i
       enddo!j
!
! *********************************************************************
! **              arrays for back-substitution (lubksb)              ** 
! *********************************************************************
! jzilch and kzilch have same number of total elements
! both contain non-zeros in lower trianglular matrix 
!
!
! ------------------ first loop of back-substitution ------------------
!
       do     i             = 2, ischang(ncs)
        kzilch(i)           = 0 
        i1                  = i - 1
        do     j            = 1, i1    
         kcnta              = kcnta + 1
         if (lzero(i,j).eq.1) then 
          kzilch(i)         = kzilch(i) + 1
          kcntb             = kcntb     + 1
          kcnt              = kcnt      + 1 
          iarray2           = iarray2   + 1
          kzero(kcnt)       = j
          jarraypt(i,j)     = iarray2 
         endif
        enddo!j
       enddo!i
!
! ----------------- second loop of back-substitution ------------------
!
! mzilch contains non-zeros for upper triangular matrix, where back-
! substitution occurs. 
!
       do     i             = ischang(ncs), 1, -1
        mzilch(i)           = 0 
        i2                  = i + 1
        do     j            = i+1, ischang(ncs)
         mcnta              = mcnta + 1
         if (lzero(i,j).eq.1) then 
          mzilch(i)         = mzilch(i) + 1
          mcntb             = mcntb     + 1
          mcnt              = mcnt      + 1
          iarray2           = iarray2   + 1
          mzero(mcnt)       = j
          jarraypt(i,j)     = iarray2 
         endif
        enddo!j
       enddo!i
!
! *********************************************************************
! * fill jarraypt with remaining array points (along diagonal)        *
! *********************************************************************
!
       do     i             = 1, ischang(ncs) 
        iarray2             = iarray2 + 1
        jarraypt(i,i)       = iarray2 
       enddo
! 
       iarray(ncsp)         = iarray2 
       kntarray             = kcnta + mcnta + ischang(ncs)
!
! *********************************************************************
! *** change izero and jzero arrays so their values point to new    ***
! ***              array positions defined in jarraypt              *** 
! *********************************************************************
!
! jarraypt = identifies the one-dimensional array point for each two-
!            dimensional point i,j
! iarray   = the length of the one-dimensional array holding all
!            sparse matrix points = sparse-matrix dimension
! izer2    = used to identify the 1-d array point for each k,j value
!            found in the first major loop of matrix decomposition
! izero    = used to find the 1-d array point for each i,k value
!            found in the same loop.
!
       do     j                  = 1, ischang(ncs)
!
! ------------------- first loop of decompostion ----------------------
!
        ijtlo(j,ncsp)             = ijtot + 1
        do     i                  = 2, ischang(ncs)
         izil                     = izilch(j,i) 
         if (izil.gt.0) then 
          ijtot                   = ijtot + 1
          nremain                 = izil
          nfive                   = (nremain + 0.0001) / 5 
          nremain                 =  nremain - nfive   * 5 
          nfour                   = (nremain + 0.0001) / 4 
          nremain                 =  nremain - nfour   * 4
          nthree                  = (nremain + 0.0001) / 3  
          nremain                 =  nremain - nthree  * 3 
          ntwo                    = (nremain + 0.0001) / 2   
          nremain                 =  nremain - ntwo    * 2  
          none                    = (nremain + 0.0001)  
          nremain                 =  nremain - none
!
          ijval(ijtot)            = jarraypt(i,j) 
          idl5( ijtot)            = idecomp     + 1
          idh5( ijtot)            = idecomp     + nfive  
          idl4( ijtot)            = idh5(ijtot) + 1
          idh4( ijtot)            = idh5(ijtot) + nfour   
          idl3( ijtot)            = idh4(ijtot) + 1
          idh3( ijtot)            = idh4(ijtot) + nthree
          idl2( ijtot)            = idh3(ijtot) + 1
          idh2( ijtot)            = idh3(ijtot) + ntwo
          idl1( ijtot)            = idh2(ijtot) + 1
          idh1( ijtot)            = idh2(ijtot) + none
          idecomp                 = idh1(ijtot)
!
          do     ic          = idl5(ijtot), idh5(ijtot)
           ka                = izerok(iccount+1) 
           kb                = izerok(iccount+2) 
           kc                = izerok(iccount+3) 
           kd                = izerok(iccount+4) 
           ke                = izerok(iccount+5) 
           iccount           = iccount + 5
           ikdeca(ic)        = jarraypt(i,ka)
           ikdecb(ic)        = jarraypt(i,kb)
           ikdecc(ic)        = jarraypt(i,kc)
           ikdecd(ic)        = jarraypt(i,kd)
           ikdece(ic)        = jarraypt(i,ke)
           kjdeca(ic)        = jarraypt(ka,j)
           kjdecb(ic)        = jarraypt(kb,j)
           kjdecc(ic)        = jarraypt(kc,j)
           kjdecd(ic)        = jarraypt(kd,j)
           kjdece(ic)        = jarraypt(ke,j)
          enddo!ic
!
          do     ic          = idh5(ijtot) + 1, idh4(ijtot)
           ka                = izerok(iccount+1) 
           kb                = izerok(iccount+2) 
           kc                = izerok(iccount+3) 
           kd                = izerok(iccount+4) 
           iccount           = iccount + 4 
           ikdeca(ic)        = jarraypt(i,ka)
           ikdecb(ic)        = jarraypt(i,kb)
           ikdecc(ic)        = jarraypt(i,kc)
           ikdecd(ic)        = jarraypt(i,kd)
           kjdeca(ic)        = jarraypt(ka,j)
           kjdecb(ic)        = jarraypt(kb,j)
           kjdecc(ic)        = jarraypt(kc,j)
           kjdecd(ic)        = jarraypt(kd,j)
          enddo!ic 
!
          do     ic          = idh4(ijtot) + 1, idh3(ijtot)
           ka                = izerok(iccount+1) 
           kb                = izerok(iccount+2) 
           kc                = izerok(iccount+3) 
           iccount           = iccount + 3 
           ikdeca(ic)        = jarraypt(i,ka)
           ikdecb(ic)        = jarraypt(i,kb)
           ikdecc(ic)        = jarraypt(i,kc)
           kjdeca(ic)        = jarraypt(ka,j)
           kjdecb(ic)        = jarraypt(kb,j)
           kjdecc(ic)        = jarraypt(kc,j)
          enddo!ic
!
          do     ic          = idh3(ijtot) + 1, idh2(ijtot)
           ka                = izerok(iccount+1) 
           kb                = izerok(iccount+2) 
           iccount           = iccount + 2 
           ikdeca(ic)        = jarraypt(i,ka)
           ikdecb(ic)        = jarraypt(i,kb)
           kjdeca(ic)        = jarraypt(ka,j)
           kjdecb(ic)        = jarraypt(kb,j)
          enddo!ic
!
          do     ic          = idh2(ijtot) + 1, idh1(ijtot)
           ka                = izerok(iccount+1) 
           iccount           = iccount + 1 
           ikdeca(ic)        = jarraypt(i,ka)
           kjdeca(ic)        = jarraypt(ka,j)
          enddo

         endif 
        enddo!i
!
        ijthi(j,ncsp)           = ijtot 
!
! ------------------ diagonal term of decompostion --------------------
!
        jarrdiag(j,ncsp)        = jarraypt(j,j)
!
! ------------------- second loop of decompostion ---------------------
!
        jloz1(j,ncsp)           = jccount + 1
        do     i                = 1, jzilch(j) 
         jccount                = jccount + 1
         ia                     = jzero(jccount)
         jzeroa(jccount)        = jarraypt(ia,j)
        enddo!i
        jhiz1(j,ncsp)           = jccount 
!
       enddo!j
!
! *********************************************************************
! ** create more back-substitution arrays to increase efficiency     ** 
! *********************************************************************
!
! ------------------ first loop of back-substitution ------------------
!
       kztlo(ncsp)          = kztot + 1 
       do     i             = 2, ischang(ncs)
        kzil                = kzilch(i) 
        if (kzil.gt.0) then 
         kztot              = kztot + 1
         nremain            = kzil
         nfive              = (nremain + 0.0001) / 5 
         nremain            =  nremain - nfive   * 5 
         nfour              = (nremain + 0.0001) / 4 
         nremain            =  nremain - nfour   * 4
         nthree             = (nremain + 0.0001) / 3  
         nremain            =  nremain - nthree  * 3 
         ntwo               = (nremain + 0.0001) / 2   
         nremain            =  nremain - ntwo    * 2  
         none               = (nremain + 0.0001)  
         nremain            =  nremain - none
!
         ikztot(kztot)      = i 
         kbl5( kztot)       = kbsub       + 1
         kbh5( kztot)       = kbsub       + nfive  
         kbl4( kztot)       = kbh5(kztot) + 1
         kbh4( kztot)       = kbh5(kztot) + nfour   
         kbl3( kztot)       = kbh4(kztot) + 1
         kbh3( kztot)       = kbh4(kztot) + nthree
         kbl2( kztot)       = kbh3(kztot) + 1
         kbh2( kztot)       = kbh3(kztot) + ntwo
         kbl1( kztot)       = kbh2(kztot) + 1
         kbh1( kztot)       = kbh2(kztot) + none
         kbsub              = kbh1(kztot)
!
         do     kc          = kbl5(kztot), kbh5(kztot)
          kzeroa(kc)        = kzero(kccount+1) 
          kzerob(kc)        = kzero(kccount+2) 
          kzeroc(kc)        = kzero(kccount+3) 
          kzerod(kc)        = kzero(kccount+4) 
          kzeroe(kc)        = kzero(kccount+5) 
          kccount           = kccount + 5
         enddo!kc 
!
         do     kc          = kbl4(kztot), kbh4(kztot)
          kzeroa(kc)        = kzero(kccount+1) 
          kzerob(kc)        = kzero(kccount+2) 
          kzeroc(kc)        = kzero(kccount+3) 
          kzerod(kc)        = kzero(kccount+4) 
          kccount           = kccount + 4 
         enddo!kc
!
         do     kc          = kbl3(kztot), kbh3(kztot)
          kzeroa(kc)        = kzero(kccount+1) 
          kzerob(kc)        = kzero(kccount+2) 
          kzeroc(kc)        = kzero(kccount+3) 
          kccount           = kccount + 3  
         enddo!kc
!
         do     kc          = kbl2(kztot), kbh2(kztot)
          kzeroa(kc)        = kzero(kccount+1) 
          kzerob(kc)        = kzero(kccount+2) 
          kccount           = kccount + 2  
         enddo!kc
!
         do     kc          = kbl1(kztot), kbh1(kztot)
          kzeroa(kc)        = kzero(kccount+1) 
          kccount           = kccount + 1  
         enddo!kc
        endif 
       enddo!i
       kzthi(ncsp)          = kztot 
!
! ----------------- second loop of back-substitution ------------------
!
       do     i             = ischang(ncs), 1, -1
        mzil                = mzilch(i) 
        if (mzil.gt.0) then 
         mztot              = mztot + 1
         nremain            = mzil
         nfive              = (nremain + 0.0001) / 5 
         nremain            =  nremain - nfive   * 5 
         nfour              = (nremain + 0.0001) / 4 
         nremain            =  nremain - nfour   * 4
         nthree             = (nremain + 0.0001) / 3  
         nremain            =  nremain - nthree  * 3 
         ntwo               = (nremain + 0.0001) / 2   
         nremain            =  nremain - ntwo    * 2  
         none               = (nremain + 0.0001)  
         nremain            =  nremain - none
!
         imztot(i,ncsp)     = mztot  
         mbl5(  mztot)      = mbsub       + 1
         mbh5(  mztot)      = mbsub       + nfive  
         mbl4(  mztot)      = mbh5(mztot) + 1
         mbh4(  mztot)      = mbh5(mztot) + nfour   
         mbl3(  mztot)      = mbh4(mztot) + 1
         mbh3(  mztot)      = mbh4(mztot) + nthree
         mbl2(  mztot)      = mbh3(mztot) + 1
         mbh2(  mztot)      = mbh3(mztot) + ntwo
         mbl1(  mztot)      = mbh2(mztot) + 1
         mbh1(  mztot)      = mbh2(mztot) + none
         mbsub              = mbh1(mztot)
!
         do     mc          = mbl5(mztot), mbh5(mztot)
          mzeroa(mc)        = mzero(mccount+1) 
          mzerob(mc)        = mzero(mccount+2) 
          mzeroc(mc)        = mzero(mccount+3) 
          mzerod(mc)        = mzero(mccount+4) 
          mzeroe(mc)        = mzero(mccount+5) 
          mccount           = mccount + 5
         enddo!mc
!
         do     mc          = mbl4(mztot), mbh4(mztot)
          mzeroa(mc)        = mzero(mccount+1) 
          mzerob(mc)        = mzero(mccount+2) 
          mzeroc(mc)        = mzero(mccount+3) 
          mzerod(mc)        = mzero(mccount+4) 
          mccount           = mccount + 4 
         enddo!mc
!
         do     mc          = mbl3(mztot), mbh3(mztot)
          mzeroa(mc)        = mzero(mccount+1) 
          mzerob(mc)        = mzero(mccount+2) 
          mzeroc(mc)        = mzero(mccount+3) 
          mccount           = mccount + 3  
         enddo!mc
!
         do     mc          = mbl2(mztot), mbh2(mztot)
          mzeroa(mc)        = mzero(mccount+1) 
          mzerob(mc)        = mzero(mccount+2) 
          mccount           = mccount + 2  
         enddo!mc
!
         do     mc          = mbl1(mztot), mbh1(mztot)
          mzeroa(mc)        = mzero(mccount+1) 
          mccount           = mccount + 1  
         enddo!mc
        endif 
       enddo!i
!
! *********************************************************************
! **        check dimensions and print out array savings             ** 
! *********************************************************************
!
       if (icnt    .gt. mxcount2  .or. jcnt     .gt. mxcount3   .or. 
     1     kcnt    .gt. mxcount3  .or. mcnt     .gt. mxcount3   .or. 
     2     iccount .gt. mxcount2  .or. jccount  .gt. mxcount3   .or. 
     3     kccount .gt. mxcount3  .or. mccount  .gt. mxcount3   .or. 
     4     ijtot   .gt. mxcount3  .or. idecomp  .gt. mxcount3   .or.  
     5     kztot   .gt. mxcount4  .or. kbsub    .gt. mxcount4   .or.  
     6     mztot   .gt. mxcount4  .or. mbsub    .gt. mxcount4   .or.  
     7     iarray2 .gt. mxarray) then
!
        write(iout,9)
     1  mxcount2, icnt,        mxcount3, jcnt,
     2  mxcount3, kcnt,        mxcount3, mcnt,
     3  mxcount2, iccount,     mxcount3, jccount,
     4  mxcount3, kccount,     mxcount3, mccount,
     5  mxcount3, ijtot,       mxcount3, idecomp,
     6  mxcount4, kztot,       mxcount4, kbsub,
     7  mxcount4, mztot,       mxcount4, mbsub,
     8  mxarray,  iarray2     
        stop
       endif
!
  9    format('one of the dimensions below is too small:',/,
     1        'dimension: mxcount2 = ',i6,' variable: icnt     = ',i6,/,
     2        'dimension: mxcount3 = ',i6,' variable: jcnt     = ',i6,/,
     3        'dimension: mxcount3 = ',i6,' variable: kcnt     = ',i6,/,
     4        'dimension: mxcount3 = ',i6,' variable: mcnt     = ',i6,/,
     5        'dimension: mxcount2 = ',i6,' variable: iccount  = ',i6,/,
     6        'dimension: mxcount3 = ',i6,' variable: jccount  = ',i6,/,
     7        'dimension: mxcount3 = ',i6,' variable: kccount  = ',i6,/,
     8        'dimension: mxcount3 = ',i6,' variable: mccount  = ',i6,/,
     9        'dimension: mxcount3 = ',i6,' variable: ijtot    = ',i6,/,
     1        'dimension: mxcount3 = ',i6,' variable: idecomp  = ',i6,/,
     2        'dimension: mxcount4 = ',i6,' variable: kztot    = ',i6,/,
     3        'dimension: mxcount4 = ',i6,' variable: kbsub    = ',i6,/,
     4        'dimension: mxcount4 = ',i6,' variable: mztot    = ',i6,/,
     5        'dimension: mxcount4 = ',i6,' variable: mbsub    = ',i6,/,
     6        'dimension: mxarray  = ',i6,' variable: iarray2  = ',i6)
!
        write(iout,10)ncsp,kount0a,kount0,kntarray,iarray2,icnta,icntb,
     1                 jcnta,jcntb,kcnta,kcntb,mcnta,mcntb 
!
 10     format(/'param    poss matrix points -- nonzeros -- ncsp=',i4/
     1          'initmat  ',4x,i10,9x,i10/
     2          'finmat   ',4x,i10,9x,i10/
     3          'decomp1  ',4x,i10,9x,i10/
     4          'decomp2  ',4x,i10,9x,i10/
     5          'backsb1  ',4x,i10,9x,i10/
     6          'backsb2  ',4x,i10,9x,i10/)
!
! *********************************************************************
! *           set coefficients of the integration method              * 
! *********************************************************************
!
! parameters used in smvgear
! --------------------------
! pertst   = coefficients used to select the step-size and order. thus,
!            only about one-percent accuracy needed. see gear(1971)
!            or hindmarsh '73 ucid-30059.  
! pertdat  = pertst
! perts2   = pertst**2
! aset     = parameters for determining the order of the integration method
!            and for calculation the matrix p.
! mstep    = maximum number of corrector iterations allowed
! hmin     = minimum time-step allowed (sec)
! maxord   = maximum order of the method used
! mbetween = maximum number of steps between calls to pderiv
! nqq      = order of the integration method
!
      if (ifdid.eq.0) then 
       ifdid = 1 
!
       pertdat(:,1)= (/ 
     +       2.0,  4.5,  7.333, 10.42,   13.7,     17.15,     1.0/)
       pertdat(:,2) =(/
     +       3.0,  6.0,  9.167, 12.5,    15.98,     1.0,      1.0/)
       pertdat(:,3) =(/
     +       1.0,  1.0,  0.5,    0.1667,  0.04133,  0.008267, 1.0/ )
!
! adams-moulton coefficients
!
!    2       2.0, 12.0, 24.0,   37.89,   53.33,    70.08,    87.97,
!    4      12.0, 24.0, 37.89,  53.33,   70.08,    87.97,     1.0,
!    6       1.0,  1.0,  2.0,    1.0,     0.3157,   0.07407,  0.0139 / 
!
       mstep          = 3
       hmin           = 1.0e-15 
       maxord         = 5
       mbetween       = 50
!
       do     i        = 1, 3
        do     nqq     = 1, morder
         pertst(nqq,i) = pertdat(nqq,i)
         perts2(nqq,i) = pertst( nqq,i) * pertst(nqq,i) 
        enddo
       enddo
!
       do     nqq     = 1, 7
        enqq1(nqq)    = 0.5 / float(nqq    )
        enqq2(nqq)    = 0.5 / float(nqq + 1) 
        enqq3(nqq)    = 0.5 / float(nqq + 2)
        conpst(nqq)   = 1.0 / (pertst(nqq,1) * enqq3(nqq)) 
        conp15(nqq)   = 1.5 * conpst(nqq)
       enddo
!
       do     i2   = 1, 6 
        aset(i2,2) = 1.0
        aset(i2,8) = 0.
       enddo
!
       aset(1,1)   = 1.0                                                        
!
       aset(2,1)   = 2.0    /    3.0  
       aset(2,3)   = 1.0    /    3.0  
!
       aset(3,1)   = 6.0    /   11.0
       aset(3,3)   = 6.0    /   11.0   
       aset(3,4)   = 1.0    /   11.0                                            
!
       aset(4,1)   = 12.0   /   25.0                                            
       aset(4,3)   =   .70    
       aset(4,4)   =   .20   
       aset(4,5)   =   .020  
!
       aset(5,1)   =   60.0 /  137.0
       aset(5,3)   =  225.0 /  274.0  
       aset(5,4)   =   85.0 /  274.0  
       aset(5,5)   =   15.0 /  274.0  
       aset(5,6)   =    1.0 /  274.0                                            
!
       aset(6,1)   =  180.0 /  441.0 
       aset(6,3)   =  406.0 /  441.0   
       aset(6,4)   =  735.0 / 1764.0  
       aset(6,5)   =  175.0 / 1764.0
       aset(6,6)   =   21.0 / 1764.0  
       aset(6,7)   =    1.0 / 1764.0  
!
      endif
!     endif ifdid.eq.0
!
! *********************************************************************
! ********************** end of subroutine ksparse ********************
! *********************************************************************
!
      return                                                             
      endsubroutine ksparse

!============================================================

      subroutine backsub(ncsp)
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
      endsubroutine backsub

!============================================================

      subroutine decomp(ncsp)
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
      endsubroutine decomp

!============================================================

      subroutine pderiv(ncsp)
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

!============================================================

      subroutine subfun(ncsp)
      !input variables
      integer, intent(in) :: ncsp
      !local variables
      integer :: nkn,ja,jb,jc,nh,k,nk2,nh2,jspc,npl,nl5,nh5,nl4,nh4,nl3,nh3
      integer :: nl2,nl1,nh1,nc,nk0,nk1,nk3,nk4,n
      integer :: nk,i,jnew,kloop
      real(kind=DP) :: concmult,fracn
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!      SSSSSSS  U     U  BBBBBBB  FFFFFFF  U     U  N     N  
!      S        U     U  B     B  F        U     U  NN    N 
!      SSSSSSS  U     U  BBBBBBB  FFF      U     U  N  N  N 
!            S  U     U  B     B  F        U     U  N    NN
!      SSSSSSS  UUUUUUU  BBBBBBB  F        UUUUUUU  N     N
!
! *********************************************************************
! *  this subroutine evaluates the first derivative of each ordinary  *  
! *                  differential equation (ode)                      * 
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call subfun.f from smvgear.f with                                * 
! *     ncs  = 1..ncsgas for gas chemistry                            *
! *     ncsp = ncs        for daytime   gas chem                      *  
! *     ncsp = ncs   +ics for nighttime gas chem                      *  
! *********************************************************************
!
! example
! -------
!
! species:         a,   b,   c
! concentrations: [a], [b], [c]
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
! *********************************************************************
!
! concmult  = product of concentrations in a rate. if two  
!             consecutive reactions have the same species reacting
!             (eg a + b --> c and a + b --> d + e) then use the 
!             same value of concmult for both reactions.
! cnew      = init (and final) species conc (# cm-3-air or moles l-1-h2o)
! gloss     = first derivative = sum of prod. minus loss rates for a species
! irma,b,c  = locates reordered active species numbers  
! ischan    = number of odes.
! lossra..  = reaordered reaction rate numbers for each loss (and prod) term
! ktloop    = number of grid-cells in a grid-block
! nsubfun   = counts the number of times this routine is called
! rrate     = forward rate coefficient 
!           = s-1                                 for rates with 1 reactant
!           = l-h2o mole-1 s-1  or cm**3 #-1 s-1  for rates with 2 reactants 
!           = l**2-h2o m-2 s-1  or cm**6 #-2 s-1  for rates with 3 reactants 
! trate     = reaction rate  moles l-1 -h2o s-1 or # cm-3 s-1 
! 
! *********************************************************************
! *                      set rates of reaction                        *
! *********************************************************************
!
!
!      print*,'check subfun1'
!      nsubfun        = nsubfun + 1
      nfdh1          = nfdh2 + ioner(ncsp) 
!
! *********************************************************************
! *     first derivatives for rates with three active loss terms      *
! *********************************************************************
!
      do     nkn     = 1, nfdh3  
       ja            = irma(nkn)
       jb            = irmb(nkn)
       jc            = irmc(nkn)
       nh            = nkn + nallr
       do     k      = 1, ktloop
        trate(k,nkn) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jb)*cnew(k,jc) 
        trate(k,nh) = -trate(k,nkn) 
       enddo
      enddo
!
! *********************************************************************
! *     first derivatives for rates with two active loss terms        *
! *********************************************************************
!
      do     nkn     = nfdl2, nfdrep
       ja            = irma(nkn)
       jb            = irmb(nkn)
       nh            = nkn + nallr
       do     k      = 1, ktloop
        trate(k,nkn) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jb) 
        trate(k,nh) = -trate(k,nkn) 
       enddo
      enddo
!
! *********************************************************************
! *     first derivatives for rates with two active loss terms and    *
! *     where the subsequent reaction has the same reactants but a    *
! *     different rate.                                               *
! *********************************************************************
!
      do     nkn     = nfdrep1, nfdh2, 2
       ja            = irma(nkn)
       jb            = irmb(nkn)
       nk2           = nkn + 1
       nh            = nkn + nallr
       nh2           = nk2 + nallr
       do     k      = 1, ktloop
        concmult     = cnew(k,ja)   * cnew(k,jb) 
        trate(k,nkn) = rrate(k,nkn) * concmult
        trate(k,nk2) = rrate(k,nk2) * concmult
        trate(k,nh)  = -trate(k,nkn) 
        trate(k,nh2) = -trate(k,nk2) 
       enddo
      enddo
!
! *********************************************************************
! *     first derivatives for rates with one active loss term         *
! *********************************************************************
!
      do     nkn     = nfdl1, nfdh1
       ja            = irma(nkn)
       nh            = nkn + nallr
       do     k      = 1, ktloop
        trate(k,nkn) = rrate(k,nkn) * cnew(k,ja) 
        trate(k,nh) = -trate(k,nkn) 
       enddo
      enddo
!
! *********************************************************************
! *                  initialize first derivative = 0                  *
! *********************************************************************
!
      do     jspc      = 1, ischan
       do     k        = 1, ktloop
        gloss(k,jspc)  = 0. 
       enddo
      enddo
!
! *********************************************************************
! * sum net (not reproduced) kinetic and photo gains and losses for   *
! * each species.                                                     *
! *********************************************************************
! sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!
      do     npl       = npllo(ncsp), nplhi(ncsp)
       jspc            = jspnpl(npl)
       nl5             = npl5(  npl)
       nh5             = nph5(  npl)
       nl4             = npl4(  npl)
       nh4             = nph4(  npl)
       nl3             = npl3(  npl)
       nh3             = nph3(  npl)
       nl2             = npl2(  npl)
       nh2             = nph2(  npl)
       nl1             = npl1(  npl)
       nh1             = nph1(  npl)
!
! ***********************  sum 5 terms at a time  ********************* 
!
       do     nc       = nl5, nh5
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        nk2            = lossrc(nc)
        nk3            = lossrd(nc)
        nk4            = lossre(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
     1                 - trate(k,nk1)  - trate(k,nk2)
     2                 - trate(k,nk3)  - trate(k,nk4)  
        enddo
       enddo
!
! ***********************  sum 4 terms at a time  ********************* 
!
       do     nc       = nl4, nh4 
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        nk2            = lossrc(nc)
        nk3            = lossrd(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
     1                 - trate(k,nk1)  - trate(k,nk2)
     2                 - trate(k,nk3)  
        enddo
       enddo
!
! ***********************  sum 3 terms at a time  ********************* 
!
       do     nc       = nl3, nh3  
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        nk2            = lossrc(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
     1                 - trate(k,nk1)  - trate(k,nk2)
        enddo
       enddo
!
! ***********************  sum 2 terms at a time  ********************* 
!
       do     nc       = nl2, nh2   
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0) 
     1                 - trate(k,nk1)      
        enddo
       enddo
!
! ***********************  sum 1 term at a time  ********************** 
!
       do     nc       = nl1, nh1    
        nk0            = lossra(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
        enddo
       enddo
      enddo!npl
!
! *********************************************************************
! *  sum production term for reactions where products fractionated    *
! *********************************************************************
!
      do     n         = nfrlo(ncsp), nfrhi(ncsp)
       jspc            = jspcnfr(n)
       nkn             = nknfr(  n)
       fracn           = fracnfr(n)
       do     k        = 1, ktloop
        gloss(k,jspc)  = gloss(k,jspc) + fracn * trate(k,nkn)       
       enddo
      enddo
!
! *********************************************************************
! **********************  end of subroutine subfun  *******************
! *********************************************************************
!
!      print*,'check subfun2'
      return
      endsubroutine subfun

!============================================================

      subroutine update (ncs, ncsp)
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
      endsubroutine update
!
      endmodule mod_smvgear_core
