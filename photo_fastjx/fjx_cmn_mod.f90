!------------------------------------------------------------------------------
!    'FJX_CMN_MOD.f90'  for fast-JX code v 7.0+ (prather 9/12)                  !
!------------------------------------------------------------------------------
!
! NB - ALL of these common variables are set paramters,
!    They are NOT to be used as variables for a local solution
!    Thus this entire set is 'in' only after it is initialized
!-----------------------------------------------------------------------
!
! !DESCRIPTION: FJX_CMN contains fast-JX variables
!
!
! !INTERFACE:
!
      MODULE FJX_CMN_MOD

      IMPLICIT NONE
      PUBLIC

!-----------------------------------------------------------------------

      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
      INTEGER, PARAMETER ::  JXL_=100, JXL1_=JXL_+1
      ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid (mid-level)
      INTEGER, PARAMETER ::  JXL2_=2*JXL_+2
      ! WX_  = dim = no. of wavelengths in input file
      INTEGER, PARAMETER ::  WX_=18
      ! X_   = dim = max no. of X-section data sets (input data)
      INTEGER, PARAMETER ::  X_=72
      ! A_   = dim = no. of Aerosol/cloud Mie sets (input data)
      INTEGER, PARAMETER ::  A_=40
      ! C_   = dim = no. of cld-data sets (input data)
      INTEGER, PARAMETER ::  C_=16
      ! W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
      INTEGER, PARAMETER ::  W_=18    ! W_= 08, 12 or 18
      ! N_  = no. of levels in Mie scattering arrays
      !     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
      INTEGER, PARAMETER ::  N_=601
      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      INTEGER, PARAMETER ::  M_=4
      ! M2_ = 2*M_ = 8, replaces MFIT
      INTEGER, PARAMETER ::  M2_=2*M_
!-----------------------------------------------------------------------
      ! 4 Gauss pts = 8-stream
      REAL*8, DIMENSION(M_), PARAMETER  ::                            &
                        EMU = [.06943184420297d0, .33000947820757d0,  &
                               .66999052179243d0, .93056815579703d0]
      REAL*8, DIMENSION(M_), PARAMETER  ::                            &
                        WT  = [.17392742256873d0, .32607257743127d0,  &
                               .32607257743127d0, .17392742256873d0]
!-----------------------------------------------------------------------

      ! ZZHT: scale height (cm)
      REAL*8, PARAMETER   :: ZZHT = 5.d5
      ! RAD: Radius of Earth (cm)
      REAL*8, PARAMETER   :: RAD = 6375.d5
      ! ATAU: heating rate (factor increase from one layer to the next)
      REAL*8, PARAMETER   :: ATAU = 1.120d0
      ! ATAU0: minimum heating rate
      REAL*8, PARAMETER   :: ATAU0 = 0.010d0
      ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
      INTEGER, PARAMETER  :: JTAUMX = (N_ - 4*JXL_)/2

!---- Variables in file 'FJX_spec.dat' (RD_XXX)

      ! WBIN: Boundaries of wavelength bins
      REAL*8  WBIN(WX_+1)
      ! WL: Centres of wavelength bins - 'effective wavelength'
      REAL*8  WL(WX_)
      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      REAL*8  FL(WX_)

      REAL*8  QO2(WX_,3)   ! QO2: O2 cross-sections
      REAL*8  QO3(WX_,3)   ! QO3: O3 cross-sections
      REAL*8  Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      REAL*8  QQQ(WX_,3,X_)
      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      REAL*8  QRAYL(WX_+1)
      ! TQQ: Temperature for supplied cross sections
      REAL*8  TQQ(3,X_)
      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      INTEGER LQQ(X_)

      ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*6  TITLEJX(X_)
      ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1  SQQ(X_)

!---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)

      ! QAA: Aerosol scattering phase functions
      REAL*8  QAA(5,A_)
      ! WAA: 5 Wavelengths for the supplied phase functions
      REAL*8  WAA(5,A_)
      ! PAA: Phase function: first 8 terms of expansion
      REAL*8  PAA(8,5,A_)
      ! RAA: Effective radius associated with aerosol type
      REAL*8  RAA(A_)
      ! SAA: Single scattering albedo
      REAL*8  SAA(5,A_)
      ! DAA: density (g/cm^3)
      REAL*8  DAA(A_)
      ! NAA: Number of categories for scattering phase functions
      INTEGER NAA

!---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)

      ! QCC: Cloud scattering phase functions
      REAL*8  QCC(5,C_)
      ! WCC: 5 Wavelengths for supplied phase functions
      REAL*8  WCC(5,C_)
      ! PCC: Phase function: first 8 terms of expansion
      REAL*8  PCC(8,5,C_)
      ! RCC: Effective radius associated with cloud type
      REAL*8  RCC(C_)
      ! SCC: Single scattering albedo
      REAL*8  SCC(5,C_)
      ! DCC: density (g/cm^3)
      REAL*8  DCC(C_)
      ! NCC: Number of categories for cloud scattering phase functions
      INTEGER NCC

!---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)

      ! WMM: U Michigan aerosol wavelengths
      REAL*8  WMM(6)
      ! UMAER: U Michigan aerosol data sets
      REAL*8  UMAER(3,6,21,33)

!---- Variables in file 'atmos_std.dat' (RD_PROF)

      ! T and O3 reference profiles
      REAL*8, DIMENSION(51,18,12) :: TREF, OREF


      INTEGER NJX,NW1,NW2

      END MODULE FJX_CMN_MOD
