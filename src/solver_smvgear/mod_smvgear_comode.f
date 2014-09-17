!=====================================================================
! Module: mod_smvgear_comode
! Common variables for smvgear solver
! Original code by Mark Jacobson
! Organized into a module by Yuzhong Zhang, 9/8/2014
!=====================================================================
      module mod_smvgear_comode
      use module_model_parameter,only:
     +    MDP=>DP, MAX_NSPEC, MAX_STR1, MAX_NREAC, MAX_NPROD,
     +    MAX_NRXN,MAX_NPHOTO
      implicit none
      !integer, parameter :: DP = 8       
      integer, parameter :: DP = MDP

! *********************************************************************
!                   COMMON VARIABLES: PARAMETERS
! *********************************************************************        
! ****************** coordinate-system parameters *********************
! KBLOOP   = maximum number of grid points in a vectorized block 
!            should range from 512 (below which vectorization decreases)
!            to 1024 (above which, array space is limited) 
      integer, parameter :: KBLOOP  = 15 

! ************************* gas-phase parameters **********************
! IGAS    = maximum number of gases, active + inactive 
! NMRATE  = maximum number of rates constants (max # reactions)
! NMTRATE = maximum number of kinetic + photo reactions 
! IPHOT   = maximum number of photo rates
! MAXGL   = maximum number of gains/losses for given array 
! MAXGL2  = 1/2 maxgl
! MAXGL3  = 1/5 maxgl
! MAXGL5  = maxgl unless huge reaction set, then about 2.5 x maxgl
! NMRPORD = maximum number of species in a reaction rate 
! ICS     = number of smvgear equation sets: 3 gaschem + 1 aqchem + 1 growth 
! ICP     = number of smvgear chemistry sets x 2 + 1 (for growth)
! MORDER  = maximum order for gear parameters for dimension purposes
      !integer, parameter ::  IGAS    =   139 
      integer, parameter ::  IGAS    =   MAX_NSPEC
      !integer, parameter ::  IPHOT   =    52
      integer, parameter ::  IPHOT   =   MAX_NPHOTO
      !integer, parameter ::  NMRATE  =   300
      integer, parameter ::  NMRATE  =   MAX_NRXN - MAX_NPHOTO
      integer, parameter ::  ICS     =     1
      integer, parameter ::  ICP     = ICS*2
      integer, parameter ::  NCSGAS  = ICS

      integer, parameter ::  NMTRATE = NMRATE + IPHOT
      integer, parameter ::  MORDER  = 7            
      integer, parameter ::  MAXGL   = 140
      integer, parameter ::  MAXGL2  = 70          
      integer, parameter ::  MAXGL3  = 35
      integer, parameter ::  MAXGL4  = 10          
      integer, parameter ::  MAXGL5  = MAXGL
      !integer, parameter ::  NMRPROD = 16              
      integer, parameter ::  NMRPROD = MAX_NREAC + MAX_NPROD

! ****************** parameters to minimize array space ***************
! MXCOUNT2,3 = arrays sizes used to minimize matrix space
! MXARRAY    = maximum length of matrix put in one-dimensional array
      integer, parameter ::  MXGSAER  = IGAS  
      integer, parameter ::  MXARRAY  = 2100         
!
      integer, parameter ::  MXCOUNT2 = MXGSAER * 80  
      integer, parameter ::  MXCOUNT3 = MXGSAER * 25
      integer, parameter ::  MXCOUNT4 = MXGSAER * 8 
      
      integer, parameter ::  STRLEN   = 14     
      !integer, parameter ::  NMREAC   = 3
      integer, parameter ::  NMREAC   = MAX_NREAC - 1
      !integer, parameter ::  NALLREAC = NMREAC + 1
      integer, parameter ::  NALLREAC = MAX_NREAC
      integer, parameter ::  NMPROD   = NMRPROD - NALLREAC
      integer, parameter ::  NPRODLO  = NALLREAC + 1
      integer, parameter ::  NPRODHI  = NPRODLO + NMPROD - 1
!
! *********************************************************************
!                   COMMON VARIABLES: STRINGS
! *********************************************************************
      character(len=STRLEN), dimension(0:MXGSAER, ICS) :: namencs
      character(len=STRLEN), dimension(ICS)            :: chemtyp     
!
! *********************************************************************
!                   COMMON VARIABLES: NON-ARRAY                 
! *********************************************************************
!      integer :: ktloop ,nlat, nlong, nlayer, nvert, nloop, ntloop
      integer :: jlooplo, ktloop      
      integer :: ioreac, iprates    
      real(kind=DP) :: fracdec,   time
      real(kind=DP) :: smal1,     smal2,     smal3
      integer ::   ifdid,   ifnever,  ifnone,   icnt,     jcnt,
     1             kcnt,    mcnt,    iccount, jccount, kccount, mccount,  
     2             idecomp, kbsub,   mbsub,   ijtot,   kztot,   mztot, 
     3             nplfun,  nfrcoun, npdcoun, npltot
      integer  ::  ischan,  nfdh3,   nfdl2,   nfdh2,  nfdl1,  nfdh1,
     1             nfdrep,  nfdrep1, nfdl0,   nallr   
      real(kind=DP) ::  hmax,  hmin,  r1delt,  delt,  timremain, xelaps,
     1                  told,  rdelt, xelaplast, rmserr     
      integer ::  nqq,     mstep,  maxord, mbetween!,nsubfun,  npderiv,
!     1            nsftot,  npdtot,    nsttot,   nfail,  ifail,   lfail
!     2            ifailtot, lfailtot, nfailtot 
      integer ::iout      
!
! *********************************************************************
!                     COMMON VARIABLES:  ARRAY                 
! *********************************************************************
      !real(kind=DP), dimension(ITLOOP) :: errmx2
      !integer, dimension(ITLOOP) :: jreorder
      integer, dimension(ICS) ::
     1  nmoth,  ntspec,   jphotrat, isgainr,   isporl,  nogaine,    
     2  nspec,  ntrates,  isgaine,  nspcsolv,  ischang, nrates,  
     3  itwor,  ithrr,    inorep     
      integer, dimension(ICP) ::
     1  nolosp,  ngnfrac, nolosrat, iarray,  nallrat, kztlo, kzthi,
     2  ioner,   npllo,   nplhi,    nfrlo,   nfrhi,   npdlo, npdhi 
      real(kind=DP), dimension(ICS)    :: abst2, errmax, timeintv
      real(kind=DP), dimension(ICP)    :: hmaxuse
      real(kind=DP), dimension(6, ICS) :: abtol
      real(kind=DP), dimension(MXGSAER) ::  aporl
      integer, dimension(NMTRATE, ICS) ::     
     1  iaprod,   nolosrn,  nruse,   nrrep,  ncequat, 
     2  noldfnew, nkoner,   nktwor,  nkthrr
      integer, dimension(NMTRATE) :: irma,  irmb,  irmc
      integer, dimension(NMTRATE*2, ICS) :: newfold
      integer, dimension(MAXGL) :: newnk
      real(kind=DP), dimension(MAXGL5, ICS) :: fracp
      integer, dimension(MAXGL2, ICS) ::  nreacoth, lgasbino
      integer, dimension(MAXGL3, ICS) ::  nknlosp,  losinacp
      integer, dimension(MAXGL5, ICS) ::  ignfrac,  nkgnfrac
      real(kind=DP), dimension(KBLOOP, MXGSAER) ::
     1  cnew,   cest,    gloss,   chold,  vdiag, dtlos,  corig
      real(kind=DP), dimension(KBLOOP, MXGSAER*7) ::   conc
      real(kind=DP), dimension(KBLOOP, NMTRATE) ::     rrate
      real(kind=DP), dimension(KBLOOP, NMTRATE,3) ::   urate 
      real(kind=DP), dimension(KBLOOP, NMTRATE*2) ::   trate
      real(kind=DP), dimension(KBLOOP, 0:MXARRAY) ::  cc2
      real(kind=DP), dimension(IPHOT, ICS) :: 
     1      nkphotrat, npphotrat, nknphotrt
      real(kind=DP), dimension(MXGSAER, ICS) :: fracgain
      integer, dimension(MXGSAER, ICS) ::
     1  numlost, numgfrt,  numgaint, ngaine,  igainr, 
     2  iporl,   igaine,   isolvspc, inewold,  mappl
      integer, dimension(MXGSAER) :: isaporl
      integer, dimension(MXGSAER, ICP) :: numporl, numloss, numgain
      integer, dimension(MXGSAER, MAXGL, ICS) :: jporl
      integer, dimension(MXGSAER, MXGSAER) :: isparder
      integer, dimension(MXGSAER) :: jzilch,  kzilch,  mzilch 
      integer, dimension(MXGSAER, MXGSAER) :: lzero, jarraypt, izilch
      integer, dimension(MXGSAER, ICP) ::
     1  jarrdiag, jloz1, jhiz1,  ijtlo, ijthi, imztot 
      real(kind=DP), dimension(MXCOUNT4) ::  fracnfr
      real(kind=DP), dimension(MXCOUNT2) ::  fracpl
      integer, dimension(MXCOUNT3) ::
     1  jzero,   kzero,   mzero,   jzeroa, ikdeca,  kjdeca,   
     2  ikdecc,  kjdecc,  ikdecd,  kjdecd, ikdecb,  kjdecb,
     3  ikdece,  kjdece,  ijval      
      integer, dimension(MXCOUNT2) ::
     1  izerok,  ipospd,  iialpd ,  nkpdterm
      integer, dimension(MXCOUNT4) ::
     1   lossra,  lossrb, lossrc,  lossrd, lossre,
     2   kzeroa,  mzeroa, kzerob,  mzerob, kzeroc,  mzeroc,  
     3   kzerod,  mzerod, kzeroe,  mzeroe, ikztot,  jspnpl,   
     4   nknfr ,  jspcnfr 
      integer, dimension(MXCOUNT3) ::
     1  idh5,  idh4,  idh3, idh2,  idh1,  idl5,
     2  idl4,  idl3,  idl2, idl1
      integer, dimension(MXCOUNT4) ::
     1  kbh5,  kbh4,  kbh3, kbh2,  kbh1,  kbl5,
     2  kbl4,  kbl3,  kbl2, kbl1,  
     3  mbh5,  mbh4,  mbh3, mbh2,  mbh1,  mbl5,
     4  mbl4,  mbl3,  mbl2, mbl1,  
     5  nph5,  nph4,  nph3, nph2,  nph1,  npl5,
     6  npl4,  npl3,  npl2, npl1  
      real(kind=DP), dimension(NMRPROD, NMTRATE, ICS) :: fkoef, fk2
      integer, dimension(NMRPROD, NMTRATE, ICS) :: irm, irm2
      real(kind=DP), dimension(10, 8) :: aset
      real(kind=DP), dimension(MORDER) ::
     1  enqq2,  enqq3,  conpst,  enqq1,  conp15
      real(kind=DP), dimension(MORDER, 3) :: perts2, pertst
!      integer, dimension(ICS)         :: nblockuse
!      integer, dimension(KBLOOP, ICS) :: jlowvar, ktlpvar
!$OMP THREADPRIVATE(ischan,nfdh3, nfdl2,nfdh2, nfdl1, nfdh1)
!$OMP THREADPRIVATE(nfdrep, nfdrep1, nfdl0, nallr,time)
!$OMP THREADPRIVATE(ktloop, jlooplo)
!$OMP THREADPRIVATE(hmax, r1delt, delt, timremain,xelaps) 
!$OMP THREADPRIVATE(told, rdelt, xelaplast, rmserr) 
!$OMP THREADPRIVATE(cc2, cnew, cest, gloss, chold, vdiag)
!$OMP THREADPRIVATE(dtlos,conc, rrate,urate,trate,corig)
!$OMP THREADPRIVATE(nqq)
!$OMP THREADPRIVATE(irma, irmb, irmc)
 
      endmodule mod_smvgear_comode

