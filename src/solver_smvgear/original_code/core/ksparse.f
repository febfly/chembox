      subroutine ksparse(ncs, ncsp)
!
      use mod_comode
      implicit none
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
      end
