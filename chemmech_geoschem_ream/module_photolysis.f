      module module_photolysis
      implicit none
      integer,parameter     :: mxnjrxn=80
      integer,parameter     :: mxnbin=7
      integer,parameter     :: mxnjval=50
      integer,parameter     :: mxnpdep=7

      integer               :: nbin
      integer               :: njval
      integer               :: npdep
      integer               :: njrxn

      real*8,dimension(mxnjrxn)                :: jfacta
      integer,dimension(mxnjrxn)               :: branch
      character(len=7),dimension(mxnjrxn)      :: jlabel 
      character(len=7),dimension(mxnjrxn)      :: rnames 

      real*8,dimension(3,mxnjval)              :: tqq
      real*8,dimension(mxnbin,2,mxnjval-3)     :: qqq
      real*8,dimension(mxnbin,3)               :: qo2,qo3,q1d
      integer                                  :: nw1, nw2
      real*8,dimension(mxnbin+1)               :: wbin
      real*8,dimension(mxnbin)                 :: wl
      real*8,dimension(mxnbin)                 :: fl
      real*8,dimension(mxnbin)                 :: qrayl
      real*8,dimension(mxnbin)                 :: qbc
      character(len=7),dimension(mxnjval)      :: titlej
      integer,dimension(mxnpdep)               :: pdepf
      real*8,dimension(mxnbin,mxnpdep)         :: zpdep
      real*8,dimension(mxnbin,3)               :: mglypdep

      integer,dimension(mxnjrxn)               :: jind
      integer,dimension(mxnjval)               :: jpdep

      real*8,dimension(mxnjval)                :: valj
      real*8,dimension(mxnjrxn)                :: jrate

      public :: photolysis_read
      public :: photolysis_jrate
      public :: photolysis_index
 
      private :: rd_js
      private :: rd_tjpl
      private :: flint
      private :: xseco3
      private :: xsec1d
      private :: xseco2
      private :: tfaca_f
      private :: tfac_f
      private :: tfac0_f
      private :: qq1_f
      private :: qq2_f
      contains

      subroutine photolysis_read(path_radj,path_jvspec)
      character(len=*)   :: path_radj,path_jvspec
      integer,parameter  :: u=78
      call rd_js  (u,path_radj)
      call rd_tjpl(u,path_jvspec)
      endsubroutine photolysis_read

C-----------------------------------------------------------------------
c  Reread the ratj.d file to map photolysis rate to reaction
c  Read in quantum yield 'jfacta' and fastj label 'jlabel'
C-----------------------------------------------------------------------
c     jfacta    Quantum yield (or multiplication factor) for photolysis
c     jlabel    Reference label identifying appropriate J-value to use
c     ipr       Photolysis reaction counter - should total 'jppj'
C-----------------------------------------------------------------------
      subroutine rd_js(nj1,namfil)
      integer,intent(in)             :: nj1
      integer                        :: ipr, i
      character(len=*),intent(in)    :: namfil
      character(len=120)             :: cline
     
      !Reread the ratj.d file to map photolysis rate to reaction
      ipr=0
      open(nj1,file=trim(namfil),status='old',form='formatted')

      cline="         "
      do while(cline(2:5).ne.'9999')
         read(nj1,'(a)') cline
         if (cline(1:1).eq.'#'.or.cline(5:5).eq.'$') cycle
         if (cline(2:5).eq.'9999') exit
         ipr=ipr+1

         !Read in quantum yield jfacta and fastj label jlabel
         read(cline(79:83),'(f5.1)') jfacta(ipr)
         read(cline(86:92),'(a7)')   jlabel(ipr)
         jfacta(ipr)=jfacta(ipr)/100.d0

         !read reaction names and branch numbers
         read (cline(7:13),"(a7)") rnames(ipr)
         rnames(ipr) = trim(rnames(ipr))
         branch(ipr) = 1
         do i=1,ipr-1
            if (rnames(ipr) == rnames(i)) branch(ipr) = branch(i) + 1
         enddo
      enddo
      close(nj1)
      njrxn = ipr

      !Print details to standard output
c      write(6,1100) ipr
c      write(6,1200) (i, jlabel(i), jfacta(i),i=1,ipr)

 1100 format(' Fast-J Photolysis Scheme: considering ',i2,' reactions')
 1200 format(3x,10(3(i2,': ',a7,' (Q.Y. ',f5.3,') '),/,3x))
      endsubroutine rd_js
    
c-----------------------------------------------------------------------
c  Read in wavelength bins, solar fluxes, Rayleigh parameters, temperature-
c  dependent cross sections and Rayleigh/aerosol scattering phase functions
c  with temperature dependences. Current data originates from JPL'97
C-----------------------------------------------------------------------
c     NAMFIL   Name of spectral data file (jv_spec.dat)
c     NJ1      Channel number for reading data file
c     NJVAL    Number of species to calculate J-values for
c     NWWW     Number of wavelength bins, from NW1:NW2
c     WBIN     Boundaries of wavelength bins
c     WL       Centres of wavelength bins - 'effective wavelength'
c     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
c     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
c     QBC      Black Carbon absorption extinct. (specific cross-sect.) (m2/g)
c     QO2      O2 cross-sections
c     QO3      O3 cross-sections
c     Q1D      O3 => O(1D) quantum yield
c     TQQ      Temperature for supplied cross sections
c     QQQ      Supplied cross sections in each wavelength bin (cm2)
c     NAA      Number of categories for scattering phase functions
c     QAA      Aerosol scattering phase functions
c     NK       Number of wavelengths at which functions supplied (set as 4)
c     WAA      Wavelengths for the NK supplied phase functions
c     PAA      Phase function: first 8 terms of expansion
c     RAA      Effective radius associated with aerosol type
c     SSA      Single scattering albedo
c
c     npdep    Number of pressure dependencies
c     zpdep    Pressure dependencies by wavelength bin
c     jpdep    Index of cross sections requiring pressure dependence
c     lpdep    Label for pressure dependence
c  NOTES:
c  (1 ) Updated to include new pressure-dependancy function for GLYX and MGLY. 
c        (tmf, 1/7/09)
c  (2 ) Added a pressure-dependancy function selector 'pdepf'. 
c        (tmf, ccc, 1/7/09)
C-----------------------------------------------------------------------

      subroutine rd_tjpl(nj1,namfil)
      character(len=*)   :: namfil
      integer            :: nj1

      integer i, j, k, iw, nk, nqqq
      character(len=7)   :: lpdep(7)
      character(len=7)   :: title

      tqq=0d0
      qqq=0d0
      open(nj1, file=trim(namfil))
      read(nj1,'(a)') title
      read(nj1,'(10x,14i5)') njval,nbin,nw1,nw2
      !check
      if (njval.gt.mxnjval) then
         print*,'njval',njval,'should be less than mxnjval',mxnjval
         stop
      endif
      if (nbin.gt.mxnbin) then
         print*,'nbin',nbin,'should be less than mxnbin',mxnbin
         stop
      endif
C------------NQQQ = no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)
C- NQQQ is changed to NJVAL-1 because there are 2 dummy species at the end
C used for acetone pressure dependency only. (ccc, 4/20/09)
C- prior to 4/20/09    
      read(nj1,102) (wbin(iw),iw=1,nbin)
      read(nj1,102) (wbin(iw+1),iw=1,nbin)
      read(nj1,102) (wl(iw),iw=1,nbin)
      read(nj1,102) (fl(iw),iw=1,nbin)
      read(nj1,102) (qrayl(iw),iw=1,nbin)
      read(nj1,102) (qbc(iw),iw=1,nbin)   !  from loiusse et al. [jgr, 1996]
C---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      do k=1,3
        read(nj1,103) title,tqq(k,1), (qo2(iw,k),iw=1,nbin)
        if (k.eq.1) titlej(1)=title
      enddo
      do k=1,3
        read(nj1,103) title,tqq(k,2), (qo3(iw,k),iw=1,nbin)
        if (k.eq.1) titlej(2)=title
      enddo
      do k=1,3
        read(nj1,103) title,tqq(k,3), (q1d(iw,k),iw=1,nbin)
        if (k.eq.1) titlej(3)=title
      enddo

C---Read remaining species:  X-sections at 2 T's
      nqqq=njval-3+2
      do j=1,nqqq
        read(nj1,103) title,tqq(1,j+3),(qqq(iw,1,j),iw=1,nbin)
        titlej(j+3)=title
        read(nj1,103) title,tqq(2,j+3),(qqq(iw,2,j),iw=1,nbin)
      enddo
      read(nj1,'(a)') title

c---Pressure dependencies
      read(nj1,104) npdep
      do k=1,npdep
         read(NJ1,105) lpdep(k), pdepf(k), (zpdep(iw,k),iw=1,nbin)
         !--------------------------------------
         ! Special treatment for MGLY pressure dependency
         ! (tmf, 11/16/06)
         !--------------------------------------
         if ( pdepf(k) .eq. 4 ) then
            ! pass zpdep to mglypdep
            mglypdep(:,1) = zpdep(:,k)
            read(nj1,105) lpdep(k), pdepf(k), (mglypdep(iw,2),iw=1,nbin)
            read(nj1,105) lpdep(k), pdepf(k), (mglypdep(iw,3),iw=1,nbin)
         endif
      enddo
      read(nj1,'(a)') title
      close(nj1)

      jind=0
      jpdep=0
C---Set mapping index
      print*,'njval',njval,'njrxn',njrxn
      do j=1,njval
         do k=1,njrxn
            if (jlabel(k).eq.titlej(j)) jind(k)=j
         enddo
         do k=1,npdep
            if (lpdep(k).eq.titlej(j)) jpdep(j)=k
         enddo
      enddo
      do k=1,njrxn
         if (jfacta(k).eq.0.d0)
     &      write(6,*) 'Not using photolysis reaction ',k
         if (jind(k).eq.0) then
            if (jfacta(k).eq.0.d0) then
               jind(k)=1
            else
               write(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
            stop
          endif
        endif
      enddo

  102 FORMAT(10X,7E10.3)
  103 FORMAT(A7,F3.0,7E10.3)
  104 FORMAT(13x,i2)
  105 FORMAT(A7,2x,I1,7E10.3)

      endsubroutine rd_tjpl


c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      subroutine photolysis_jrate(fff,Temp,Pres,iday)
      real*8,dimension(mxnbin),intent(in):: fff !Actinic flux #/cm2/s
      real*8,intent(in) :: Temp, Pres !K,hPa
      integer,optional  :: iday! julian day

C
      real*8 qo2tot, qo3tot, qo31d, qo33p, qqqt
      real*8 solf, tfact

C     Parameters for Solar distance compensation
      real*8,parameter ::  pi=3.14159265358979324d0

C     Physical constants
      real*8,parameter ::   Na=6.02217d23
      real*8,parameter ::   r=8.3143d0

C     Add Pressure dependancy function selector PF. (tmf, 1/7/09) 
      integer i, j, k, l, pf
      real*8 qptemp

C     For new pressure-dependency algorithm: (tmf, 1/7/09) 
      real*8 xp, xa, xb, xc

C     For new pressure dependency algo. for acetone
      real*8 tfaca,tfac0,tfac1, tfac2
      real*8 qqqa , qq1a , qq1b, qq2

C     Scale actinic flux (FFF) by Solar distance factor (SOLF)
      if (present(iday))then
         solf=1d0-(0.034d0*cos(dble(iday-172)*2d0*pi/365d0))
      else 
         solf=1d0
      endif

      valj=0d0
      jrate=0d0

      do k=nw1,nw2                       ! using model 't's here
         qo2tot= xseco2(k,Temp)
         valj(1) = valj(1) + qo2tot*fff(k)
         qo3tot= xseco3(k,Temp)
         qo31d = xsec1d(k,Temp)*qo3tot
         qo33p = qo3tot - qo31d
         valj(2) = valj(2) + qo33p*fff(k)
         valj(3) = valj(3) + qo31d*fff(k)
       enddo

C------Calculate remaining J-values with T-dep X-sections 
       do j=4,njval
         valj(j) = 0.d0
         tfact = 0.d0
         l = jpdep(j)

C        To choose different forms of pres. dependancy. (ccc, 4/20/09)
         if ( l.ne.0 ) pf = pdepf(l)

         if(tqq(2,j).gt.tqq(1,j)) tfact = max(0.d0,min(1.d0,
     $        (Temp-tqq(1,j))/(tqq(2,j)-tqq(1,j)) ))

C        FAST_JX introduces a new pres. dependancy for acetone (ccc, 4/20/09)
C        Special calculations for the temperature interpolation factors
         if ( pf.eq.2 ) then
            tfaca=tfaca_f(Temp, j      )
            tfac0=tfac0_f(Temp, j+1    )
            tfac1=tfac_f (Temp, njval+1)
            tfac2=tfac_f (Temp, njval+2)
         else if ( pf.eq.3 ) then
            tfaca=tfaca_f(Temp, j-1    )
            tfac0=tfac0_f(Temp, j      )
         endif

         do k=nw1,nw2
           qqqt = qqq(k,1,j-3) + (qqq(k,2,j-3) - qqq(k,1,j-3))*tfact
           if(l.eq.0) then
              valj(j) = valj(j) + qqqt*fff(k)
           else

              ! Select pressure dependancy function (tmf, 1/31/06)
              if (pf .eq. 1) then
C----------------------------------------------------------------------
C Prior to 9/17/99
C Original form for acetaldehyde P-dep -- believed to be incorrect (pjc)
C             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*
C     $                   (1.d0+zpdep(K,L)*(pj(i)+pj(i+1))*0.5d0)
C----------------------------------------------------------------------
C Essentially the change is the replacement of the factor
C
C   (1 + a P)     with               1
C                           ---------------------
C                             (1 + b density)
C
C where a and b are constants, P is pressure, and density is the 
C density of air in molec-cm(-3)   (pjc, 9/17/99)
C----------------------------------------------------------------------
              valj(j)=valj(j)+qqqt*fff(k)/(1 +
     1                 (zpdep(k,l)*Na*1d-6 /(R*Temp)) *
     2                 Pres*1d2)
             else if ( pf .eq. 4 ) then
C-----------------------------------------------------------------------
C For MGLY
C       y = a + ( b * exp(-p/c) )
C    where y is the ratio between Omega(p) / Omega(p=0);
C          x is the atmospheric pressure [Pa]
C          a,b,c are MGLYPDEP(:,1), MGLYPDEP(:,2), MGLYPDEP(:,3)
C-----------------------------------------------------------------------
                 xp = Pres*1d2   ! pressure [Pa]
                 xa = mglypdep( k, 1 )
                 xb = mglypdep( k, 2 )
                 xc = mglypdep( k, 3 )
                 qptemp = 1d0

                 if ( abs( xc ) .ge. 1d-10 ) then
                    qptemp = xa + ( xb * exp(-xp/xc) )
                 endif

                 valj(j) = valj(j) + qqqt*fff(k)*qptemp

              else if ( pf.eq.2 ) then
C             Acetone pressure dependency from FAST-JX (ccc, 4/20/09)
C             J1(acetone-a) ==> CH3CO + CH3
C             Special values for Xsect
                 qqqa = qq1_f (tfaca, j      , k            )
                 qq2  = qq2_f (tfac0, j+1    , k, Temp      )
                 qq1a = qq1_f (tfac1, njval+1, k            )
                 qq1b = qq1_f (tfac2, njval+2, k            ) * 4.d-20

                 valj(j) = valj(j) + fff(k)*qqqa *
     1            (1.d0-qq2)/(qq1a + (qq1b*na*1d-6 /(r*Temp)) *
     2            Pres*1d2)
              else if ( PF.eq.3 ) then
C             Second acetone pressure dependency from FAST-JX (ccc, 4/20/09)
C             J2(acetone-b) ==> CH3 + CO + CH3
C             Special values for Xsect
                 qqqa = qq1_f (tfaca, j-1    , k            )
                 qq2  = qq2_f (tfac0, j      , k, temp      )
                 valj(j) = valj(j) + fff(k)*qqqa*qq2
              else
                 print*,'Formula for presurre dependent photolysis'
     +                  //'(pdepf=',pf,') not specified'
                 stop
              endif
           endif
         enddo
       enddo
       valj(1:njval)=valj(1:njval)*solf
       do j=1,njrxn
          jrate(j)=valj(jind(j))*jfacta(j)
          !print*,j,jrate(j)
       enddo
      endsubroutine photolysis_jrate

c-----------------------------------------------------------------------
c  three-point linear interpolation function
c-----------------------------------------------------------------------
      function flint (tint,t1,t2,t3,f1,f2,f3)
      real*8 tint,t1,t2,t3,f1,f2,f3
      real*8 flint
      if (tint .le. t2)  then
        if (tint .le. t1)  then
          flint  = f1
        else
          flint = f1 + (f2 - f1)*(tint -t1)/(t2 -t1)
        endif
      else
        if (tint .ge. t3)  then
          flint  = f3
        else
          flint = f2 + (f3 - f2)*(tint -t2)/(t3 -t2)
        endif
      endif
      endfunction flint

C-----------------------------------------------------------------------
c  Quantum yields for O3 --> O2 + O(1D) interpolated across 3 temps
C-----------------------------------------------------------------------
      function xsec1d(k,ttt)
      integer k
      real*8 ttt, xsec1d
      xsec1d =
     f  flint(ttt,tqq(1,3),tqq(2,3),tqq(3,3),q1d(k,1),q1d(k,2),q1d(k,3))
      endfunction xsec1d

c-----------------------------------------------------------------------
c  Cross-sections for O2 interpolated across 3 temps; No S_R Bands yet!
C-----------------------------------------------------------------------
      function xseco2(k,ttt)
      integer k
      real*8 ttt, xseco2
      xseco2 =
     f  flint(ttt,tqq(1,1),tqq(2,1),tqq(3,1),qo2(k,1),qo2(k,2),qo2(k,3))
      endfunction xseco2


c-----------------------------------------------------------------------
c  Cross-sections for O3 for all processes interpolated across 3 temps
C-----------------------------------------------------------------------
      function xseco3(k,ttt)
      integer k
      real*8 ttt, xseco3
      xseco3  =
     f  flint(ttt,tqq(1,2),tqq(2,2),tqq(3,2),qo3(k,1),qo3(k,2),qo3(k,3))
      endfunction xseco3

c-----------------------------------------------------------------------
! !DESCRIPTION: Calculates temperature interpolation factors for acetone
c-----------------------------------------------------------------------
      function tfaca_f(ttt, iv)
      ! Index of the specie in jv_spec.dat (should be between 4 and NJVAL)
      integer :: iv
      ! Temperature in 1 grid box
      real*8  :: ttt
      ! Temperature interpolation factor
      real*8  :: tfaca_f

      tfaca_f = (ttt-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))
      tfaca_f = max(0.d0, min(1.d0, tfaca_f))
      endfunction tfaca_f

      function tfac0_f(ttt, iv)
      integer :: iv
      real*8  :: ttt
      real*8  :: tfac0_f
      tfac0_f = ( (ttt-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv)) )**2
      if (ttt .lt. tqq(1,iv)) then
         tfac0_f = (ttt - 210.d0)/(tqq(1,iv)-210.d0)
      endif
      tfac0_f = max(0.d0, min(1.d0, tfac0_f))
      end function tfac0_f

      
      function tfac_f(ttt, iv)
      integer :: iv
      real*8  :: ttt
      real*8  :: tfac_f
      real*8  :: tt200
      tt200 = min(300.d0, max(200.d0, ttt))
      tfac_f = (tt200-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))
      endfunction tfac_f

! !DESCRIPTION: This routine computes the cross-section for acetone.
      function qq2_f(tfac0, iv, k, ttt)
      integer :: iv
      integer :: k
      real*8  :: ttt,tfac0
      ! xsect (total abs) for acetone
      real*8  :: qq2_f
      qq2_f  = qqq(k,1,iv-3) + (qqq(k,2,iv-3)-qqq(k,1,iv-3))*tfac0
      if (ttt .lt. tqq(1,iv)) then
         qq2_f = qqq(k,1,iv-3)*tfac0
      endif
      end function qq2_f

      function qq1_f(tfac, iv, k)
      integer :: iv
      integer :: k
      real*8  :: tfac
      real*8  :: qq1_f
      qq1_f = qqq(k,1,iv-3) + (qqq(k,2,iv-3)-qqq(k,1,iv-3))*tfac
      end function qq1_f

c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      function photolysis_index(specname,p)result(ind)
      character(len=*) :: specname
      real*8           :: p
      integer          :: ind
      integer          :: ibrch, n

      ind=0
      ibrch=floor((p-floor(p))*10d0+0.5d0)
      do n=1,njrxn
         if (specname.eq.rnames(n).and.ibrch.eq.branch(n))then
            ind=n
            exit
         endif
      enddo
         !print*,trim(specname),ind
      
      
      endfunction
      endmodule module_photolysis

