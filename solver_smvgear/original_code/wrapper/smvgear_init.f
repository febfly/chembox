C
C File:   smvgear_init.f
C Author: Yuzhong
C
C Created on July 24, 2014, 11:43 PM
C
      subroutine smvgear_init
      use mod_comode
      implicit none

      !=========================================================!
      !            Initialize variables                         !
      !=========================================================!     
      ioreac   = 1
      iprates  = 0
      ifdid    = 0
      ifnever  = 0
      ifnone   = 0
!      nsftot   = 0
!      npdtot   = 0
!      nsttot   = 0
!      ifailtot = 0
!      lfailtot = 0
!      nfailtot = 0
      
      rmserr   = 0.
            
      namencs(0:MXGSAER,1:ICS) =" "
      ntspec   (1:ICS) = 0
      nmoth    (1:ICS) = 0
      jphotrat (1:ICS) = 0
      isgainr  (1:ICS) = 0
      isporl   (1:ICS) = 0
      nogaine  (1:ICS) = 0 
      !nouse(    1:ICS) = 0
      nspec(    1:ICS) = 0
      ntrates(  1:ICS) = 0  
      isgaine(  1:ICS) = 0 
      nspcsolv( 1:ICS) = 0 
      ischang(  1:ICS) = 0 
      nrates(   1:ICS) = 0
      itwor(    1:ICS) = 0
      ithrr(    1:ICS) = 0 
      inorep(   1:ICS) = 0
      
      nolosp(   1:ICP) = 0  
      ngnfrac(  1:ICP) = 0
      nolosrat( 1:ICP) = 0
      iarray(   1:ICP) = 0
      nallrat(  1:ICP) = 0
      kztlo(    1:ICP) = 0
      kzthi(    1:ICP) = 0
      ioner(    1:ICP) = 0
      npllo(    1:ICP) = 0
      nplhi(    1:ICP) = 0
      nfrlo(    1:ICP) = 0
      nfrhi(    1:ICP) = 0
      npdlo(    1:ICP) = 0
      npdhi(    1:ICP) = 0      
       
      fracp(   1:MAXGL,1:ICS) = 0  
      ignfrac( 1:MAXGL,1:ICS) = 0
      nkgnfrac(1:MAXGL,1:ICS) = 0     
        
      nreacoth(1:MAXGL2,1:ICS) = 0
      lgasbino(1:MAXGL2,1:ICS) = 0
        
      losinacp(1:MAXGL3,1:ICS) = 0
      
      fracgain(1:MXGSAER,1:ICS) = 0.
      numlost( 1:MXGSAER,1:ICS) = 0        
      numgfrt( 1:MXGSAER,1:ICS) = 0 
      numgaint(1:MXGSAER,1:ICS) = 0
      ngaine(  1:MXGSAER,1:ICS) = 0
      igainr(  1:MXGSAER,1:ICS) = 0
      iporl(   1:MXGSAER,1:ICS) = 0
      igaine(  1:MXGSAER,1:ICS) = 0
      isolvspc(1:MXGSAER,1:ICS) = 0
      inewold( 1:MXGSAER,1:ICS) = 0
      mappl(   1:MXGSAER,1:ICS) = 0    
      
      numloss( 1:MXGSAER,1:ICP) = 0        
      numgain( 1:MXGSAER,1:ICP) = 0
      numporl( 1:MXGSAER,1:ICP) = 0 
      
      nolosrn( 1:NMTRATE,1:ICS) = 0
      nruse(   1:NMTRATE,1:ICS) = 0
      nrrep(   1:NMTRATE,1:ICS) = 0
      ncequat( 1:NMTRATE,1:ICS) = 0 
      noldfnew(1:NMTRATE,1:ICS) = 0
      newfold( 1:NMTRATE*2,1:ICS) = 0
      nkoner(  1:NMTRATE,1:ICS) = 0
      nktwor(  1:NMTRATE,1:ICS) = 0
      nkthrr(  1:NMTRATE,1:ICS) = 0
      
      nkphotrat(1:IPHOT,1:ICS) = 0 
      nknphotrt(1:IPHOT,1:ICS) = 0   
      
      jarrdiag(1:MXGSAER,1:ICP)  = 0
      jloz1(   1:MXGSAER,1:ICP)  = 0
      jhiz1(   1:MXGSAER,1:ICP)  = 0
      ijtlo(   1:MXGSAER,1:ICP)  = 0
      ijthi(   1:MXGSAER,1:ICP)  = 0
      imztot(  1:MXGSAER,1:ICP)  = 0
      
      irm(  1:NMRPROD,1:NMTRATE,1:ICS) = 0
      irm2( 1:NMRPROD,1:NMTRATE,1:ICS) = 0
      fkoef(1:NMRPROD,1:NMTRATE,1:ICS) = 0. 
      fk2(  1:NMRPROD,1:NMTRATE,1:ICS) = 0. 
      
      jporl(1:MXGSAER, 1:MAXGL, 1:ICS) = 0
      
      cnew( 1:KBLOOP,1:MXGSAER)       = 0.
      cest( 1:KBLOOP,1:MXGSAER)       = 0.
      gloss(1:KBLOOP,1:MXGSAER)       = 0.
      chold(1:KBLOOP,1:MXGSAER)       = 0.
      vdiag(1:KBLOOP,1:MXGSAER)       = 0.
      dtlos(1:KBLOOP,1:MXGSAER)       = 0.
      corig(1:KBLOOP,1:MXGSAER)       = 0.      
      
      conc( 1:KBLOOP, 1:MXGSAER*7)    = 0.
      rrate(1:KBLOOP, 1:NMTRATE)      = 0.
      urate(1:KBLOOP, 1:NMTRATE, 1:3) = 0.
      
      trate(1:KBLOOP, 1:NMTRATE*2)    = 0.
      
      !pratk1(1:KBLOOP, 1:IPHOT)       = 0.
      !pratkd(1:KBLOOP, 1:IPHOT)       = 0.
      
      cc2(   1:KBLOOP, 0:MXARRAY)     = 0.
      !kgrp(  1:KBLOOP, 1:5)           = 0
      

      endsubroutine smvgear_init

