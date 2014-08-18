      module mod_smvgear 
      implicit none
      public :: smvgear_setup
      public :: smvgear_solve

      private :: backsub
      private :: decomp
      private :: jsparse
      private :: ksparse
      private :: pderiv
      private :: smvgear
      private :: subfun
      private :: update
      private :: smvgear_init
      private :: smvgear_setchem
      private :: smvgear_setpara
      private :: smvgear_solveeq

      contains

!PUBLIC SUBROUTINES
!======================================================================
! smvgear_setup
!=====================================================================
      subroutine smvgear_setup(chemintv,n_spec, n_active, specname_len,
     +           specname,  n_rxn, n_photo, n_mxreac, n_mxprod, reac,
     +           prod, reac_coef, prod_coef)
      use mod_comode, only : DP, STRLEN, IGAS, IPHOT, NMRATE, IPHOT,
     +                       NMREAC, NMPROD, iout
      real(kind=DP),intent(in) :: chemintv
      integer, intent(in) :: n_spec, n_active, specname_len
      integer, intent(in) :: n_rxn, n_photo, n_mxreac, n_mxprod
      character(len=specname_len),dimension(n_spec) :: specname
      integer, dimension(n_mxreac,n_rxn),intent(in) :: reac
      integer, dimension(n_mxprod,n_rxn),intent(in) :: prod
      real(kind=DP),dimension(n_mxreac, n_rxn), intent(in) :: reac_coef
      real(kind=DP),dimension(n_mxprod, n_rxn), intent(in) :: prod_coef

      character(len=STRLEN) :: typename
      integer :: ncs
        
    
      ncs = 1  
      iout = 6
      typename='CHEM'

      !Check dimensions
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

!======================================================================
! smvgear_solve
!=====================================================================
      subroutine smvgear_solve(ifsun, n_grid, n_spec, n_rxn, conc0,
     +           rate_const, conc1, flag, ifcheck)
      use mod_comode, only: KBLOOP, IGAS, NMRATE, IPHOT, iout
      !Input variables
      integer, intent(in) :: ifsun
      integer, intent(in) :: n_spec, n_rxn, n_grid
      real(kind=DP), dimension(n_grid, n_spec), intent(in) :: conc0
      real(kind=DP), dimension(n_grid, n_rxn), intent(in) :: rate_const
      logical             :: ifcheck

      !Output variables
      real(kind=DP), dimension(n_grid, n_spec), intent(out) :: conc1
      integer, intent(out) :: flag

      !Local variables
      integer :: ncs

      ncs = 1

      !Check dimensions and value of ifsun
      if (ifcheck.eq..true.) then
         if (n_grid.gt.KBLOOP) then
            write(iout,*) 'n_grid is greater than KBLOOP',
     +                    n_grid, KBLOOP
            stop
         endif
         if (ifsun.ne.1.and.ifsun.ne.2) then
            write(iout,*) 'Error: ifsun can be 1 or 2, 1:day, 2:night',
                           ifsun
            stop
         endif
      endif
 
      !Call solver
      call smvgear_solveeq(ncs, ifsun, n_grid, n_spec, n_rxn, conc0,
     +                     rate_const, conc1, flag)
      

      endsubroutine smvgear_solve

!PRIVATE SUBROUTINES
!======================================================================
! smvgear_init
!=====================================================================
      subroutine smvgear_init
      use mod_comode
      implicit none

      ioreac   = 1
      iprates  = 0
      ifdid    = 0
      ifnever  = 0
      ifnone   = 0
      rmserr   = 0.

      namencs(0:MXGSAER,1:ICS) =" "
      ntspec   (1:ICS) = 0
      nmoth    (1:ICS) = 0
      jphotrat (1:ICS) = 0
      isgainr  (1:ICS) = 0
      isporl   (1:ICS) = 0
      nogaine  (1:ICS) = 0 
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
      cc2(   1:KBLOOP, 0:MXARRAY)     = 0.

      endsubroutine smvgear_init

!======================================================================
! smvgear_setpara
!=====================================================================
      subroutine smvgear_setpara(ncs, chemintv)
      use mod_comode
      implicit none
      !Input variables
      integer, intent(in) :: ncs
      real(kind=DP), intent(in) :: chemintv
      !Local variables
      real(kind=DP) :: errmaxu, ylowu, yhiu, hmaxdayu, hmaxnit
      real(kind=DP) :: abhi, ablo
      integer       :: i

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
! smvgear_setchem
!=====================================================================
      subroutine smvgear_setchem(ncs, n_spec, n_active, specname_len, specname, 
     +   n_rxn,  n_photo, n_mxreac, n_mxprod, reac, prod, reac_coef, prod_coef,
     +   chemtypename)
      use mod_comode
      implicit none
      !Input variables
      integer, intent(in) :: ncs
      integer, intent(in) :: n_spec, n_active, n_rxn, n_photo
      integer, intent(in) :: specname_len, n_mxreac, n_mxprod

      character(len=STRLEN), dimension(n_spec), intent(in) :: specname

      integer, dimension(n_mxreac, n_rxn), intent(in) :: reac
      integer, dimension(n_mxprod, n_rxn), intent(in) :: prod
      real(kind=DP), dimension(n_mxreac, n_rxn), intent(in) :: reac_coef
      real(kind=DP), dimension(n_mxprod, n_rxn), intent(in) :: prod_coef
      character(len=STRLEN), intent(in) :: chemtypename

      !Local variables
      integer :: i, j, j1

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

!======================================================================
! backsub
!=====================================================================


      subroutine backsub(ncsp)
!
      use mod_comode
      implicit none
      !input variables
      integer, intent(in) :: ncsp

      !local variables
      integer :: ij,i,kzt,kl5,kh5,kl4,kh4,kl3,kh3,kl2,kh2,kl1,kh1,kc
      integer :: j0,ij0,ij1,ij2,ij3,ij4,j1,j2,j3,j4,k,mzt,ml5,mh5,ml4,mh4
      integer :: ml3,mh3,ml2,mh2,ml1,mh1,mc

      endmodule mod_smvgear
