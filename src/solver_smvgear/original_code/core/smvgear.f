      subroutine smvgear(ncs, ifsun,flag)
!
      use mod_comode
      implicit none
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
      end
