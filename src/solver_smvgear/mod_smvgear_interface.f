!=========================================================================
! Module mod_smvgear_interface
! Written by Yuzhong Zhang, 9/8/2014
!=========================================================================
      module mod_smvgear_interface
      implicit none

      !public functions
      public :: smvgear_setup
      public :: smvgear_solve

      !private functions
      private :: smvgear_init
      private :: smvgear_setpara
      private :: smvgear_setchem
      contains

!=========================================================================
      subroutine smvgear_setup(ncs,chemintv,n_spec, n_active,specname_len,
     +           specname,  n_rxn, n_photo, n_mxreac, n_mxprod, reac,
     +           prod, reac_coef, prod_coef)
      use mod_smvgear_comode, only : DP, STRLEN, IGAS, IPHOT, NMRATE, 
     +           IPHOT, NMREAC, NMPROD, ICS, iout
      use mod_smvgear_core,only : jsparse
      implicit none

      real(kind=DP),intent(in) :: chemintv
      integer, intent(in) :: ncs, n_spec, n_active, specname_len
      integer, intent(in) :: n_rxn, n_photo, n_mxreac, n_mxprod
      character(len=specname_len),dimension(n_spec) :: specname
      integer, dimension(n_mxreac,n_rxn),intent(in) :: reac
      integer, dimension(n_mxprod,n_rxn),intent(in) :: prod
      real(kind=DP),dimension(n_mxreac, n_rxn), intent(in) :: reac_coef
      real(kind=DP),dimension(n_mxprod, n_rxn), intent(in) :: prod_coef

      character(len=STRLEN) :: typename

      iout = 6
      typename='CHEM'

      !Check dimensions
      if (ncs.gt.ICS) then
         write(iout,*) 'Error: ncs is greater then ICS',
     +                 ncs, ICS
         stop
      endif

      if (n_spec.gt.IGAS) then
         write(iout,*) 'Error: n_spec is greater than IGAS',
     +                  n_spec,IGAS
         stop
      endif

      if (specname_len.ne.STRLEN) then
         write(iout,*) 'Error: specname_len is not equal to STRLEN',
     +                 specname_len, STRLEN
         stop
      endif

      if (n_rxn.gt.NMRATE+IPHOT) then
         write(iout,*) 'Error: n_rxn is greater than NMRATE+IPHOT',
     +                 n_rxn, NMRATE+IPHOT
         stop
      endif
      if (n_photo.gt.IPHOT) then
         write(iout,*) 'Error: n_photo is greater than IPHOT',
     +                 n_photo, IPHOT
         stop
      endif

      if (n_mxprod.gt.NMPROD) then
         write(iout,*) 'Error: n_mxprod is greater than NMPROD',
     +                 n_mxprod, NMPROD
         stop
      endif

      if (n_mxreac.gt.NMREAC) then
         write(iout,*) 'Error: n_mxprod is greater than NMREAC',
     +                 n_mxprod, NMPROD
         stop
      endif

      !initialize arrays
      call smvgear_init

      !set parameters, e.g., time step, error tolerance, etc.
      call smvgear_setpara(ncs,chemintv)

      !set chemistry structure
      call smvgear_setchem(ncs, n_spec, n_active, specname_len,specname,
     +     n_rxn, n_photo, n_mxreac, n_mxprod, reac, prod, reac_coef,
     +     prod_coef, typename)

      !setup sparse matrix
      call jsparse(ncs)

      endsubroutine smvgear_setup

!====================================================================
      subroutine smvgear_solve (ncs, ifsun, n_grid, conc0, rate_const, 
     +                          conc1, flag)
!===================================================================
! INPUT
!     ncs      : index of chemistry set
!     ifsun= 1 : daytime, solve photolysis
!          = 2 : nighttime, solve no photolysis
!     n_grid   : # of grids
!     n_spec   : # of species
!     n_rxn    : # of reactions
!     conc0    : inital concentration of species, [n_grid, n_spec]
!     rate_const: rate constants of reactions, [n_grid, n_rxn]
! OUTPUT
!     conc1    : final concentration at the end of the step, [n_grid,
!     n_spec]  
!     flag = 0 : solve succesfully
!          = 2 : solver fails. nylowdec in smvgear.f is larger than
!          threshold
!          = 3 : solver fails. jrestar  in smvgear.f is larger than
!          threshold
!==============================================================================
          use mod_smvgear_comode
          use mod_smvgear_core, only : smvgear
          implicit none
          !Input variables
          integer, intent(in) :: ncs, ifsun
          integer, intent(in) :: n_spec, n_rxn, n_grid
          real(kind=DP), dimension(KBLOOP, IGAS), intent(in) :: conc0
          real(kind=DP), dimension(KBLOOP, NMRATE), intent(in) :: rate_const

          !Output variables
          real(kind=DP), dimension(KBLOOP, IGAS), intent(out) :: conc1
          integer, intent(out)                          :: flag
          !character(len=MSGLEN), intent(out)            :: message

          !Local variables
          real(kind=DP):: chemintv
          integer      :: nkn, nk, nh, kloop
          integer      :: jnew, jold

          chemintv = timeintv(ncs)
          time = chemintv

          !Set concentration
          do jold = 1, ntspec(ncs)
             jnew = mappl(jold, ncs)
             corig(1:n_grid, jnew) = conc0(1:n_grid, jold)
          enddo!jold          

          !Set rate constant
          nfdh3   = ithrr(ncs)
          nfdl2   = nfdh3  + 1
          nfdrep  = inorep(ncs)
          nfdrep1 = nfdrep + 1
          nfdh2   = nfdh3  + itwor(ncs)
          nfdl1   = nfdh2  + 1
          nfdh1   = nfdh2  + ioner(ncs)
          nfdl0   = nfdh1  + 1
          nallr   = nallrat(ncs)

          do nkn = 1, nallr
              nk = noldfnew(nkn,ncs)
              irma(nkn) = irm2(1, nk, ncs)
              irmb(nkn) = irm2(2, nk, ncs)
              irmc(nkn) = irm2(3, nk, ncs)
          enddo!nkn

          do nkn = 1, nallr
              nk = noldfnew(nkn, ncs)
              rrate(1:n_grid, nkn) = rate_const(1:n_grid, nk)
          enddo!nkn

          do nkn = nfdl0, nallr
              nh = nkn + nallr
              trate(1:n_grid, nkn) = rrate(1:n_grid, nkn)
              trate(1:n_grid, nh ) =-rrate(1:n_grid, nkn)
          enddo!nkn

          !Solve
          ktloop = n_grid
          jlooplo = 1
          call smvgear(ncs, ifsun,flag)

          !Update concentration
          conc1(:,:) = 0d0
          if (flag.eq.0) then !if solver succeeds, update concentration
             do jnew = 1, ischang(ncs)
                jold = inewold(jnew, ncs)
                do kloop = 1, n_grid
                   conc1(kloop, jold) = max(cnew(kloop, jnew), smal2)
                   ! conc1(kloop, jold) = cnew(kloop, jnew)
                enddo!kloop
             enddo!jnew
          else !if solver fails, reset 
             conc1(:,:) = conc0(:,:)
          endif
      endsubroutine smvgear_solve


!================================================================
      subroutine smvgear_init
      use mod_smvgear_comode
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

!===================================================================
      subroutine smvgear_setpara(ncs, chemintv)
      use mod_smvgear_comode
      implicit none
      !Input variables
      integer, intent(in) :: ncs
      real(kind=DP), intent(in) :: chemintv
      !Local variables
      !real(kind=DP) :: chemintv
      real(kind=DP) :: errmaxu, ylowu, yhiu, hmaxdayu, hmaxnit
      real(kind=DP) :: abhi, ablo
      integer       :: i

      !=========================================================!
      !            Set parameters for smvgear                   !
      !=========================================================!
      !chemintv = 900.
      !ifreord  = 1
      fracdec  = 0.25
      errmaxu  = 1e-3
      ylowu    = 1e3
      yhiu     = 1e7
      hmaxdayu = 9e2
      hmaxnit  = 2e3
      rmserr   = 0.
      iout     = 6

      smal1    = 1.0d-06
      smal2    = 1.0d-100
      smal3    = 1.0d-50

      errmax(ncs)       = errmaxu
      timeintv(ncs)     = chemintv
      abst2 (ncs)       = 1./(chemintv*chemintv)
      hmaxuse (ncs)     = hmaxdayu
      hmaxuse (ICS+ncs) = hmaxnit

      !Set abtol
      abtol(1, ncs) = yhiu
      abtol(6, ncs) = ylowu
      abhi = log10(abtol(1,ncs))
      ablo = log10(abtol(6,ncs))
      do i = 2, 5
         abtol(i, ncs) = 10.**(ablo + (abhi - ablo) * float(6 - i) / 5.)
      enddo
      endsubroutine smvgear_setpara

!======================================================================

      subroutine smvgear_setchem(ncs, n_spec, n_active, specname_len,specname, 
     +   n_rxn,  n_photo, n_mxreac, n_mxprod, reac, prod, reac_coef,
     +   prod_coef, chemtypename)
      use mod_smvgear_comode
      implicit none
      !Input variables
      integer, intent(in) :: ncs
      integer, intent(in) :: n_spec, n_active, n_rxn, n_photo
      integer, intent(in) :: specname_len, n_mxreac, n_mxprod

      character(len=STRLEN), dimension(n_spec), intent(in) :: specname
      !integer,                     dimension(n_spec), intent(in) ::
      !spectype

      integer, dimension(n_mxreac, n_rxn), intent(in) :: reac
      integer, dimension(n_mxprod, n_rxn), intent(in) :: prod
      real(kind=DP), dimension(n_mxreac, n_rxn), intent(in) :: reac_coef
      real(kind=DP), dimension(n_mxprod, n_rxn), intent(in) :: prod_coef
      character(len=STRLEN), intent(in) :: chemtypename

      !Output variable
      !character(len=255) :: message     

      !Local variables
      integer :: i, j, j1

      !=========================================================!
      !            Set chemistry mechanism information          !
      !=========================================================!
      !Set comode variable
      ntspec (ncs) = n_spec
      nspec  (ncs) = n_active
      namencs(:,ncs) = ""
      namencs(1:n_spec,ncs) = specname(1:n_spec)
      chemtyp(ncs) = chemtypename

      nallrat(ncs) = n_rxn
      ntrates(ncs) = n_rxn
      nallrat(ncs+ICS) = n_rxn - n_photo
      nrates(ncs)  = n_rxn - n_photo
      jphotrat(ncs) = n_photo


      do i = 1, n_rxn
         ncequat(i,ncs) = i
      enddo !i

      do i = 1, n_photo
          nkphotrat(i,ncs) = nrates(ncs) + i
      enddo!i

      irm(:,:,ncs) = 0
      fkoef(:,:,ncs) = 0.
      do i = 1, n_rxn
         do j = 1, n_mxreac
            irm  (j, i, ncs) = reac (j, i)
            fkoef (j, i, ncs) = reac_coef(j, i)
         enddo
         do j = 1, n_mxprod
            j1 = j + nprodlo -1
            irm (j1, i, ncs) = prod (j, i)
            fkoef(j1, i, ncs) = prod_coef(j, i)
         enddo
      enddo
      endsubroutine smvgear_setchem

      endmodule mod_smvgear_interface
