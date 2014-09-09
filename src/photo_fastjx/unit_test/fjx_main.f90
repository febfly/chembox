!>>>>>>>>>>>>>>>>current code revised to JX ver 7.1 (10/13)<<<<<<<<<<<<
      program standalone

      USE FJX_CMN_MOD

      USE FJX_SUB_MOD

      USE FJX_INIT_MOD

      include 'FJX_notes.f90'

      implicit none

      real*8, parameter::  PI180 = 3.141592653589793d0/180.d0
      real*8, parameter::  MASFAC=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

      real*8, dimension(L1_) :: ETAA,ETAB, RI,TI,CLDP,AER1,AER2
      real*8 GMTAU,ALBEDO,  XLNG,YLAT,XGRD,YGRD
      real*8 PSURF, SCALEH

      integer, dimension(L1_):: NCLD,NAA1,NAA2
      integer  MONTH
      integer J,K,L, IDAY
!-----------------------------------------------------------------------

!--------------------key params sent to fast JX-------------------------
!---SZA = solar zenith angle, U0 = cos(SZA)
!---SOLF = solar flux factor for sun-earth distance
!---REFLB = Lambertian reflectivity at the Lower Boundary
      real*8                    :: U0,SZA,REFLB,SOLF
!---ZPJ(layer, JV#) = 2D array of J-values, indexed to CTM chemistry code
      real*8, dimension(JVL_,JVN_)::   ZPJ
!---turn internal PHOTO_JX print (unit=6) on/off
      logical                     :: LPRTJ
!---Independent Column Atmosphere data passed to PHOTO_JX:
!--- P = edge press (hPa), Z = edge alt (cm), T = layer temp (K)
!--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
!---
!--- R = layer rel.hum.(fraction)
!--- CLDWP = cloud water path (g/m2), AERSP = aerosol path (g/m2)
!--- NDXCLD = cld index type (liq & size, ice & scatt-phase)
!--- NDXAER = aerosol index type
!--- NB. clouds have only a single type within an ICA layer, pick liquid or ice
!---     aerosols are dimensioned with up to AN_ different types in an ICA layer
      real*8,  dimension(L1_+1)  :: PPP,ZZZ
      real*8,  dimension(L1_  )  :: TTT,DDD,RRR,OOO,LWP,IWP,REFFL,REFFI
      real*8,  dimension(L1_,AN_):: AERSP
      integer, dimension(L1_  )  :: NDXCLD
      integer, dimension(L1_,AN_):: NDXAER

      real*8        PCLD,IWC,FACTOR

!---these are the key results coming out of fast-JX core
!---they need to be dimensioned here and passed to fast-JX to fill.
      integer, parameter             ::  NJX_ = 100
      character*6, dimension(NJX_)   ::  TITLJXX
      real*8, dimension(L_,NJX_)     ::  VALJXX
      integer :: NJXX

!-----------------------------------------------------------------------
!---fast-JX:  INIT_JX is called only once to read in and store all fast-JX data:
      call INIT_FJX (TITLJXX,NJX_,NJXX)
!-----------------------------------------------------------------------


!-------READ IN test atmosphere:
      open(1,file='FJX_test.dat',status='old')
      read(1,*)
      read(1,*) IDAY
      read(1,*) MONTH
      read(1,*) GMTAU
      read(1,*) YLAT
      read(1,*) XLNG
         YGRD = YLAT*PI180
         XGRD = XLNG*PI180
      read(1,*) PSURF
      read(1,*) ALBEDO
      read(1,*)
!---for L_=57-layer EC - top layer, L_+1, does not use cloud and aerosol data
      do L=1,L_+1
        read(1,'(i3,1x,2f11.6,2x,2f5.1,f7.3,i4,f7.3,i4,f7.3,i4)')     &
          J,ETAA(L),ETAB(L),TI(L),RI(L),CLDP(L),NCLD(L),              &
                 AER1(L),NAA1(L),AER2(L),NAA2(L)
      enddo
      close(1)

!-----------------------------------------------------------------------
!---fast-JX:  SOLAR_JX is called only once per grid-square to set U0, etc.

       call SOLAR_JX(GMTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)
!-----------------------------------------------------------------------

        write(6,'(a,f8.3,f8.5)')'solar zenith angle, solar-f',SZA,SOLF
        write(6,'(a,f8.3,f8.3)') 'lat/lng',YLAT,XLNG

      do L = 1,L1_
        PPP(L) = ETAA(L) + ETAB(L)*PSURF
      enddo

!-----------------------------------------------------------------------
!---fast-JX:  ACLIM_FJX sets up climatologies for O3, T, Density and Z.

      call ACLIM_FJX (YLAT, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1_)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
!-----------------------------------------------------------------------

!---load CTM-based ICA atmospheric column data on top of climatology (1:L_)
!--e.g.
!---to convert kg (STT) in grid cell w/AREA (m2) to # molecules/cm^2
!      D_O3(I,J,L) = 6.023d26*STT(I,J,L,1)/48.d0  *1.d-4 /AREAXY(I,J)
!---to convert kg (STT) in grid cell w/AREA (m2) to PATH (g/m^2)
!      P_AERSL(I,J,L) = STT(I,J,L,1)*1.d-3/AREAXY(I,J)


        LWP(:)  =  0.d0
        IWP(:)  =  0.d0
        REFFL(:) = 0.d0
        REFFI(:) = 0.d0
        AERSP(:,:)  = 0.d0
        NDXAER(:,:) = 0
      do L = 1,L1_
        PPP(L) = ETAA(L) + ETAB(L)*PSURF
      enddo
      do L = 1,L_
        TTT(L) = TI(L)
        RRR(L) = RI(L)
        if (TTT(L) .gt. 253.d0) then
          LWP(L) = CLDP(L)
        else
          IWP(L) = CLDP(L)
        endif
        NDXAER(L,1) = NAA1(L)
        AERSP(L,1)  = AER1(L)
        NDXAER(L,2) = NAA2(L)
        AERSP(L,2)  = AER2(L)
      enddo
        ZZZ(1) = 0.d0
      do L = 1,L_
        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
        ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
      enddo
        ZZZ(L1_+1) = ZZZ(L1_) + 5.d5
        REFLB = ALBEDO
        LPRTJ = .true.

!>>>R-effective of clouds determined by main code, not FJX
!   REFF determined by user - some recommendations below
!       REFFI is a simple function of ice water content IWC (g/m3, 0.0001 to 0.1)
!          IWC = IWP / delta-Z (of layer in m, approx OK)
!   Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)
!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
      do L = 1,L_
       if (IWP(L) .gt. 1.d-5) then
         IWC = IWP(L) *100.d0 / (ZZZ(L+1)-ZZZ(L))
         IWC = max(0.001d0, IWC)
         REFFI(L) = 164.d0 * IWC**0.23d0
!       write(6,'(a,i3,3f10.4)') 'ICE:',L,IWP(L),IWC,REFFI(L)
       endif
       if (LWP(L) .gt. 1.d-5) then
         PCLD = 0.5d0*(PPP(L)+PPP(L+1))
         FACTOR = min(1.d0, max(0.d0, (PCLD-610.d0)/200.d0))
         REFFL(L) = 9.60d0*FACTOR + 12.68d0*(1.-FACTOR)
       endif
      enddo

!        call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!---fast-JX:  PHOTO_JX is called once per ICA, calculates all J-values

        call PHOTO_JX                                                 &
          (U0,SZA,REFLB,SOLF, LPRTJ, PPP,ZZZ,TTT,DDD,RRR,OOO,         &
          LWP,IWP,REFFL,REFFI,AERSP,NDXAER,L1_,AN_,       VALJXX,NJX_)
!-----------------------------------------------------------------------

!---map the J-values from fast-JX onto CTM (ZPJ) using JIND & JFACTA
       do L = 1,L_
         do J = 1,NRATJ
          if (JIND(J).gt.0) then
            ZPJ(L,J) = VALJXX(L,JIND(J))*JFACTA(J)
          else
            ZPJ(L,J) = 0.d0
          endif
         enddo
       enddo

!---Printout J's:
        write(6,'(a)')
        write(6,'(a)') ' >>>>CTM J-values taken from fast-JX----'
        write(6,'(a,i3,a,f4.1,a,f7.2,a,2f7.2,a,f8.5)')                &
        '  Day=',IDAY,' UT(hr)=',GMTAU,'  SZA=',SZA,                  &
       '  LAT x LONG=',YLAT,XLNG,' SolarFlx=',SOLF
        write(6,'(1x,a,72(i6,3x))') ' ',(K, K=1,NRATJ)
        write(6,'(1x,a,72(a6,3x))') 'L=  ',(JLABEL(K), K=1,NRATJ)
       do L=JVL_,1,-1
        write(6,'(i3,1p, 72e9.2)') L,(ZPJ(L,K),K=1,NRATJ)
       enddo

      stop
      end
