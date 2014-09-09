      subroutine jsparse(ncs)
!
      use mod_comode
      implicit none
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
      end
