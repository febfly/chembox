!=========================================================================
! Module       : module_geoschem_rxntype
! Functionality: Define reaction type and the way to calculate reaction rate 
!                using GEOS-Chem chemical mechnism
! Written by Yuzhong Zhang, 9/8/2014
!=========================================================================
      module module_geoschem_rxntype
      use module_model_parameter,only: DP
      implicit none
      integer, parameter :: MAX_NTYPE = 30
      integer, parameter :: MAX_NPARA = 21
      integer, parameter :: STRLEN    = 50

      integer :: ntype
      
      character(len=STRLEN),dimension(MAX_NTYPE) :: type_name
      character(len=1),dimension(MAX_NTYPE)      :: symbol
      character(len=STRLEN),dimension(MAX_NTYPE):: comment
      integer,dimension(MAX_NTYPE)               :: nparameter
      integer,dimension(MAX_NPARA,MAX_NTYPE)     :: parameter_pos
      integer,dimension(MAX_NTYPE)               :: preceeding_type
      integer,dimension(MAX_NTYPE)               :: succeeding_type
      integer,dimension(MAX_NTYPE)               :: if_photo_type
      integer,dimension(MAX_NTYPE)               :: nrxnline

      !pulic functions
      public  :: rxn_rate
      public  :: rxntype_id
      public  :: rxntype_init

      !private functions (helper of rate calcuation for some special reactions)
      private :: r1
      private :: fyhoro
      private :: fyrno3
      private :: arsl1k
      contains

!=========================================================================
! subroutine: rxntype_init
!             Defines reaction type
!             The definition includes
!               1. name of the reaction type 
!                  GEOS-Chem chemical mechnism input file 
!               2. symbol used in the GEOS-Chem chemical mechnism input file
!                  to denote the special reaction type
!               3. comment of the reaction type
!               4. the reaction type that appears prior to this reaction
!                  type in the GEOSChem chemical mechnism input file; 0
!                  means no requirement for the preceeding reaction
!                  type
!               5. the reaction type that appears after this reaction
!                  type in the GEOS-Chem chemical mechnism input file; 0
!                  means no requirement for the succeeding reaction
!                  type
!               6. if this is a photolysis reaction type; 0 means normal
!                  reaction type; 1 means photolysis.
!               7. 
!              ATTENTION: when modifying, be sure to modify ntype to a
!                  proper value
!=========================================================================
      subroutine rxntype_init
      parameter_pos(:,:)     = 0
      ntype                  = 16

      type_name(1)           = 'Normal'
      symbol(1)              = " "
      comment(1)             = "A*(300/T)**B*EXP(C/T)"
      preceeding_type(1)     = 0
      succeeding_type(1)     = 0
      if_photo_type(1)       = 0
      nrxnline(1)            = 1

      type_name(2)           = 'Pressure dependent 3 body'
      symbol(2)              = "P"
      comment(2)             = "Kh*G**(1/(1+(log(Kh/Kl))^2))"
      preceeding_type(2)     = 0
      succeeding_type(2)     = 0      
      if_photo_type(2)       = 0
      nrxnline(2)            = 2

      type_name(3)           = 'Addition branch of RO2+NO'
      symbol(3)              = "B"
      comment(3)             = "k*(1-fyrno3)"
      preceeding_type(3)     = 0
      succeeding_type(3)     = 4
      if_photo_type(3)       = 0
      nrxnline(3)            = 2

      type_name(4)           = 'Abstaction branch of RO2+NO'
      symbol(4)              = "A"
      comment(4)             = "k*fyrno3"
      preceeding_type(4)     = 3
      succeeding_type(4)     = 0
      if_photo_type(4)       = 0
      nrxnline(4)            = 2

      type_name(5)           = 'Equilibrim reaction (reverse)'
      symbol(5)              = "E"
      comment(5)             = ""
      preceeding_type(5)     = 2
      succeeding_type(5)     = 0
      if_photo_type(5)       = 0
      nrxnline(5)            = 1

      type_name(6)           = 'OH+HNO3'
      symbol(6)              = "X"
      comment(6)             = "K1 + K3[M] / (1 + K3[M]/K2)"
      preceeding_type(6)     = 0
      succeeding_type(6)     = 0
      if_photo_type(6)       = 0
      nrxnline(6)            = 3

      type_name(7)           = 'OH+CO'
      symbol(7)              = "Y"
      comment(7)             = ""
      preceeding_type(7)     = 0
      succeeding_type(7)     = 0
      if_photo_type(7)       = 0
      nrxnline(7)            = 1

      type_name(8)           = 'HO2/NO3+HO2'
      symbol(8)              = "Z"
      comment(8)             = "(K1 + K2)*(1+1.4E-21*[H2O]*EXP(2200/T))"
      preceeding_type(8)     = 0
      succeeding_type(8)     = 0
      if_photo_type(8)       = 0
      nrxnline(8)            = 2

      type_name(9)           = 'GLYX+OH/NO3->HO2+CO'
      symbol(9)              = "C"
      comment(9)             = "k1*([o2]+3.5d18)/(2*[o2]+3.5d18)"
      preceeding_type(9)     = 0
      succeeding_type(9)     = 0
      if_photo_type(9)       = 0
      nrxnline(9)            = 1

      type_name(10)          = 'GLYX+OH/NO3->GLCO3'
      symbol(10)             = "D"
      comment(10)            = "k1*[o2]/(2*[o2]+3.5d18)"
      preceeding_type(10)    = 0
      succeeding_type(10)    = 0
      if_photo_type(10)      = 0
      nrxnline(10)           = 1

      type_name(11)          = 'Aerosol surface loss'
      symbol(11)             = "K"
      comment(11)            = ""
      preceeding_type(11)    = 0
      succeeding_type(11)    = 0
      if_photo_type(11)      = 0
      nrxnline(11)           = 1

      type_name(12)          = 'EO2->HO2+GLYC'
      symbol(12)             = "H"
      comment(12)            = "k1*FYHORO"
      preceeding_type(12)    = 0
      succeeding_type(12)    = 0
      if_photo_type(12)      = 0
      nrxnline(12)           = 1

      type_name(13)          = 'EO2->HO2+2CH2O'
      symbol(13)             = "F"
      comment(13)            = "k1*(1-FYHORO)"
      preceeding_type(13)    = 0
      succeeding_type(13)    = 0
      if_photo_type(13)      = 0
      nrxnline(13)           = 1

      type_name(14)          = 'Temperature dependent branching'
      symbol(14)             = "V"
      comment(14)            = "k1/(1+k2)"
      preceeding_type(14)    = 0
      succeeding_type(14)    = 0
      if_photo_type(14)      = 0
      nrxnline(14)           = 2

      type_name(15)          = 'Photolysis'
      symbol(15)             = " "
      comment(15)            = ""
      preceeding_type(15)    = 0
      succeeding_type(15)    = 0
      if_photo_type(15)      = 1
      nrxnline(15)           = 1

      type_name(16)          = 'O3 photolysis'
      symbol(16)             = "A"
      comment(16)            = ""
      preceeding_type(16)    = 0
      succeeding_type(16)    = 0
      if_photo_type(16)      = 1
      nrxnline(16)           = 1
      endsubroutine rxntype_init

!=========================================================================
! function  : rxntype_id
!             return the index of reaction type defined in rxntype_init,
!             given the symbol and if it is a photolysis reaction.
!=========================================================================
      function rxntype_id(sym,ifphoto) result(id)
       character(len=1) :: sym
       integer          :: id
       integer          :: i,st,ed
       integer          :: ifphoto
       st=1
       ed=ntype
       do i = st, ed
          if (symbol(i).eq.sym.and.if_photo_type(i).eq.ifphoto) exit
       enddo
       id = i
       if (id.eq.ed+1) then
          print*,'Error: cannot find a reaction type defined in'
     +          // ' rxntype_init that uses a symbol',sym
          stop
       endif
       !if (ifphoto.eq.1) id =15!overwrite type 16...
      endfunction rxntype_id

!=========================================================================
! subroutine: rxn_rate
!             defines the way to calcuate reaction rate for different
!             reaction type
! input     :
!             1. nr        : number of reactions
!             2. tp        : reaction type [nr]
!             3. para      : parameter list [MAX_NPARA, nr]
!             4. Temp      : temperature (K)
!             5. Pres      : pressure (Pa)
!             6. O2        : O2 concentration (molecules/cm3)
!             7. N2        : N2 concentration (molecules/cm3)
!             8. H2O       : H2O concentration (molecules/cm3)
!             9. aer_area  : aerosol surface area (cm2/cm3)
!            10. aer_radius: aerosol radius (cm)
!            11. denair    : density of air (molecules/cm3)
! output    :
!             1. rate_c    : rate constant [nr]
!=========================================================================
      subroutine rxn_rate(nr, tp, para, Temp, Pres, 
     +                  O2, N2, denair, H2O, aer_area, aer_radius, 
     +                  rate_c) 
      !in out variables
      integer, intent(in)                              :: nr
      integer, dimension(nr), intent(in)               :: tp
      real(kind=DP),dimension(MAX_NPARA,nr),intent(in) :: para
      real(kind=DP),intent(in) :: Temp, Pres, O2, N2, H2O
      real(kind=DP),intent(in) :: aer_area, aer_radius, denair
      real(kind=DP),dimension(nr),intent(out)      :: rate_c

      !local variables
      integer                            :: ir, itp
      real(kind=DP),dimension(MAX_NPARA) :: p
      real(kind=DP)                      :: rate_const
      real(kind=DP)            :: M
      real(kind=DP)            :: kh, kl, xyrat, blog, fexp
      real(kind=DP)            :: k1, k2, k3, kf
      real(kind=DP)            :: klo1,khi1,xyrat1,fexp1,kco1
      real(kind=DP)            :: klo2,khi2,xyrat2,fexp2,kco2
      real(kind=DP)            :: blog1,blog2     
      real(kind=DP)            :: stk, stkcf, sqm
     
      !loop through reactions
      M = denair
      rate_c(:) = 0d0 
      do ir = 1, nr
         itp = tp(ir)
         p(:)= para(:,ir)
         select case (itp)
           !Arrhenius formulation
           case (1) 
                rate_const = r1(p(1:3),Temp)
           !P: pressure dependent
           case (2)
                kh=r1(p(1:3),Temp)*M
                kl=r1(p(4:6),Temp)
                xyrat=kh/kl
                blog=log10(xyrat)
                fexp=1d0/(1d0+blog*blog)
                rate_const=kh*p(7)**fexp/(1d0+xyrat)/M
           !B
           case (3)
                k1=r1(p(1:3),Temp)
                rate_const=k1*(1e0-fyrno3(p(4),denair,Temp))           
           !A
           case (4)
                k1=r1(p(1:3),Temp)
                rate_const=k1*fyrno3(p(4),denair,Temp)        
           !E
           case (5)
                k1=r1(p(1:3),Temp)
                kf=rate_c(ir-1)
                rate_const=kf*M/k1
           !X
           case (6)
                k1=r1(p(1:3),Temp)
                k2=r1(p(4:6),Temp)
                k3=r1(p(7:9),Temp)*M
                rate_const=k1+k3/(1d0+k3/k2)
           !V
           case (7)
                klo1=5.9d-33*(300/Temp)**(1.4d0)
                khi1=1.1d-12*(300/Temp)**(-1.3d0)
                xyrat1=klo1*M/khi1
                blog1=log10(xyrat1)
                fexp1=1.d0/(1.d0+blog1*blog1)
                kco1=klo1*M*0.6**fexp1/(1.d0+xyrat1)
                klo2=1.5d-13*(300/Temp)**(-0.6d0)
                khi2=2.1d09 *(300/Temp)**(-6.1d0)
                xyrat2=klo2*M/khi2
                blog2=log10(xyrat2)
                fexp2=1.d0/(1.d0+blog2*blog2)
                kco2=klo2*0.6**fexp2/(1.d0+xyrat2)
                rate_const=kco1+kco2                            
           !Z
           case (8)
                k1=r1(p(1:3),Temp)
                k2=r1(p(4:6),Temp)
                rate_const=(k1+k2*M)*(1d0+1.4d-21*H2O*exp(2200d0/Temp))
           !C
           case (9)
                k1=r1(p(1:3),Temp)
                rate_const=k1*(O2+3.5d18)/(2d0*O2+3.5d18)               
           !D
           case (10)
                k1=r1(p(1:3),Temp)
                rate_const=k1*O2/(2d0*O2+3.5d18)
           !K
           case (11)
                sqm=sqrt(p(1))
                stkcf=p(2)
                stk=sqrt(Temp)
                rate_const=arsl1k(aer_area,aer_radius,denair,stkcf,stk,sqm) 
           !H
           case (12)
                k1=r1(p(1:3),Temp)
                rate_const=k1*fyhoro(denair,Temp)    
           !F
           case (13)
                k1=r1(p(1:3),Temp)
                rate_const=k1*(1d0-fyhoro(denair,Temp))
           !V
           case (14)
                k1=r1(p(1:3),Temp)
                k2=r1(p(1:3),Temp)
                rate_const=k1/(1d0+k2)
           !Photolysis
           case (15:16)
                rate_const= 0d0
         endselect
         rate_c(ir) = rate_const
      enddo !ir

      endsubroutine rxn_rate

!===========================================================================
! function : r1
!            calculate rate using Arrhenius formulation
!===========================================================================
      function r1(para, Temp) result(k)
      real(kind=DP),dimension(3),intent(in) :: para
      real(kind=DP),             intent(in) :: Temp
      real(kind=DP)                         :: k
      k = para(1)
      if (para(2).ne.0e0) k = k*(3d2/Temp)**para(2)
      if (para(3).ne.0e0) k = k*exp(para(3)/Temp)
      endfunction r1

!===========================================================================
! Helper functions for special rates
!===========================================================================

!===========================================================================
! $Id: fyrno3.f,v 4.18 2001/08/31 15:17:53 bmy v4.18 $
!===========================================================================
      function fyrno3(xcarbn,zdnum,tt)result(rfyrno3)
!***************************************************************************
!          Returns organic nitrate yields YN = RKA/(RKA+RKB)
!          from RO2+NO reactions as a function
!          of the number N of carbon atoms.
!          Updated following Atkinson 1990.
!***************************************************************************
      real(kind=DP) :: xcarbn
      real(kind=DP) :: zdnum,tt
      real(kind=DP) :: yyyn,xxyn,aaa,rarb,zzyn,xf,alpha,y300,beta,xminf,xm0
      real(kind=DP) :: rfyrno3

      data y300,alpha,beta,xm0,xminf,xf/.826,1.94e-22,.97,0.,8.1,.411/
c
      xxyn = alpha*exp(beta*xcarbn)*zdnum*((300./tt)**xm0)
      yyyn = y300*((300./tt)**xminf)
      aaa=log10(xxyn/yyyn)
      zzyn = 1./(1.+ aaa*aaa )
      rarb = (xxyn/(1.+ (xxyn/yyyn)))*(xf**zzyn)
      rfyrno3 = rarb/(1. + rarb)

      ! Change ISN2 yield from RIO2+NO to 4.4% after Chen et al, 1998
      ! This reaction is tagged by a 99 for the number of carbons in chem.dat
      ! (amf, 4/20/00)
      if (xcarbn == 9.90e+01) then
         rfyrno3 = 0.044d0
      endif
      endfunction fyrno3

!============================================================================
! $Id: fyhoro.f,v 1.1.1.1 2005/02/06 17:26:09 tmf Exp $
!============================================================================
      function fyhoro( zdnum, tt )result(rfyhoro)
!
!******************************************************************************
!  Function FYHORO returns returns the branching ratio between 
!  HOC2H4O oxidation and dissociation:
!  (1) HOC2H4 --O2--> HO2 + GLYC
!  (2) HOC2H4 ------> HO2 + 2CH2O
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ZDNUM  (REAL*8) : Air density   [molec/cm3 ]
!  (2 ) TT     (REAL*8) : Temperature   [K         ]
!    
!  NOTES: 
!  (1 ) Branching ratio calculation (tmf, 2/6/05).
!
!  REFERENCES:
!  (1 ) Orlando et al., 1998. Laboratory and theoretical study of the oxy 
!        radicals in the OH- and Cl-initiated oxidation of ethene, 
!        J. Phys. Chem. A, 102, 8116-8123.
!  (2 ) Orlando et al., 2003. The atmospheric chemistry of alkoxy radicals,
!        Chem. Rev., 103, 4657-4689.
!******************************************************************************
      ! arguments
      real(kind=DP), intent(in) :: zdnum, tt
      real(kind=DP)             :: rfyhoro
      ! local variables
      real(kind=DP)             :: k1, k2, o2dnum
      o2dnum = zdnum * 0.21d0
      k1     = 6.0d-14 * exp(-550.d0/tt) * o2dnum
      k2     = 9.5d+13 * exp(-5988.d0/tt)
      rfyhoro = k1 / (k1 + k2)
      endfunction fyhoro

!============================================================================
      function arsl1k(area, radius, denair, stkcf, stk, sqm)result(rarsl1k)
C***************************************************************************C
C*      This function calculates the 1st-order loss rate of species        *C
C*         on wet aerosol surface, with formula from Dentener's            *C
C*         thesis, pp14.               jyl, July 1, 1994                   *C
C*         arsl1k (/s) = area / [ radius/dfkg + 4./(stkcf * xmms) ]        *C
C***************************************************************************C
      real(kind=DP) ::  stkcf
      real(kind=DP) ::  area, radius, stk, sqm, denair
      real(kind=DP) ::  dfkg
      real(kind=DP) ::  rarsl1k
C* area: sfc area of wet aerosols per volume of air (cm2/cc);
C* radius: radius of wet aerosol (cm), order of 0.01-10 um;
c          note that radius here is rd, not ro ;
C* dfkg: gas phase diffusion coefficient (cm2/s), order of 0.1;
C* denair: density of air in #/cc.
C* stkcf: sticking coefficient (no unit), order of 0.1;
C* stk: sqrt of temperature (K);
C* sqm: sqrt of molecular weight (g/mol);
C* xmms: mean molecular speed (cm/s), = sqrt(8R*TK/pi/M) for Maxwell 
C*       distribution.

      if(area .lt. 0. .or. radius .lt. 0.) then
c default 
         rarsl1k = 1.d-3
      else
         dfkg  = 9.45d17/denair * stk * sqrt(3.472d-2 + 1.d0/(sqm*sqm))
         rarsl1k=area/(radius/dfkg + 2.749064e-4*sqm/(stkcf*stk))
      endif
      endfunction arsl1k

      endmodule module_geoschem_rxntype
