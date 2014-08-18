      module module_reaction
      implicit none
!       integer, parameter             :: funcname_max = 20
       integer, parameter              :: type_max = 30
       integer, parameter              :: string_max = 200

       integer, parameter              :: rxn_max  = 500
!       integer, parameter              :: photo_max = 200
       integer, parameter              :: reac_max = 4
       integer, parameter              :: prod_max = 20
       integer, parameter              :: para_max = 16

       type reaction
         integer :: nreac
         integer :: nprod
         character(len=1) :: status
         integer :: sn
         integer,dimension(reac_max) :: reacs
         integer,dimension(prod_max) :: prods
         real*8 ,dimension(prod_max) :: prod_coefs
         integer :: r_type
         real*8, dimension(para_max)   :: paras
         integer :: rindex          !for photolysis reactions
       endtype reaction

       type(reaction),dimension(rxn_max) :: rxn
       integer                           :: nr
       integer                           :: nphoto
       real*8,dimension(rxn_max)         :: rate_cst

       type reaction_type
         integer :: np
         character(len=50):: type_name
         character(len=1) :: sym1
         character(len=1) :: sym2
         character(len=50):: expres
       endtype reaction_type

       type(reaction_type), dimension (type_max)  :: rxn_type
       integer                                    :: n_tp

       public :: rxn_findtype
       public :: rxn_init
       public :: rxn_add
       public :: rxn_update_rates
       public :: rxn_rate
 
       private :: rxn_readtype       
       private :: r1
       private :: r2       
      contains
!=========================================================================
       subroutine rxn_readtype
       rxn_type(1)%type_name = 'Normal'
       rxn_type(1)%np        = 3
       rxn_type(1)%sym1      = " "
       rxn_type(1)%sym2      = " "
       rxn_type(1)%expres    = "A*(300/T)**B*EXP(C/T)"

       rxn_type(2)%type_name = 'Pressure dependent 3 body'
       rxn_type(2)%np        = 7
       rxn_type(2)%sym1      = "P"
       rxn_type(2)%sym2      = "P"
       rxn_type(2)%expres    = "Kh*G**(1/(1+(log(Kh/Kl))^2))"
      
       rxn_type(3)%type_name = 'Addition branch of RO2+NO'
       rxn_type(3)%np        = 4
       rxn_type(3)%sym1      = "B"
       rxn_type(3)%sym2      = "B"
       rxn_type(3)%expres    ='k*(1-fyrno3)'
       !k1*(1-FYRNO3)

       rxn_type(4)%type_name = 'Abstaction branch of RO2+NO'
       rxn_type(4)%np        = 4
       rxn_type(4)%sym1      = "A"
       rxn_type(4)%sym2      = "A"
       rxn_type(4)%expres    ='k*fyrno3'
       !k1*FYRNO3

       rxn_type(5)%type_name = 'Equilibrim reaction (reverse)'
       rxn_type(5)%np        = 3
       rxn_type(5)%sym1      = "E"
       rxn_type(5)%sym2      = "E"
       rxn_type(5)%expres    =' '

       rxn_type(6)%type_name = 'OH+HNO3'
       rxn_type(6)%np        = 9
       rxn_type(6)%sym1      = "X"
       rxn_type(6)%sym2      = "X"
       rxn_type(6)%expres    =' '
       !K1 + K3[M] / (1 + K3[M]/K2)       

       rxn_type(7)%type_name = 'OH+CO'
       rxn_type(7)%np        = 3
       rxn_type(7)%sym1      = "Y"
       rxn_type(7)%sym2      = "Y"
       rxn_type(7)%expres    =' '
       !Complicated

       rxn_type(8)%type_name = 'HO2/NO3+HO2'
       rxn_type(8)%np        = 6
       rxn_type(8)%sym1      = "Z"
       rxn_type(8)%sym2      = "Z"
       rxn_type(8)%expres    =' '
       !(K1 + K2)*(1+1.4E-21*[H2O]*EXP(2200/T))
       
       rxn_type(9)%type_name = 'GLYX+OH/NO3->HO2+CO'
       rxn_type(9)%np        = 3 
       rxn_type(9)%sym1      = "C"
       rxn_type(9)%sym2      = "C"
       rxn_type(9)%expres    =' '
       !k1*([o2]+3.5d18)/(2*[o2]+3.5d18)

       rxn_type(10)%type_name = 'GLYX+OH/NO3->GLCO3'
       rxn_type(10)%np        = 3 
       rxn_type(10)%sym1      = "D"
       rxn_type(10)%sym2      = "D"
       rxn_type(10)%expres    =' '
       !k1*[o2]/(2*[o2]+3.5d18)

       rxn_type(11)%type_name = 'Aerosol surface'
       rxn_type(11)%np        = 3 
       rxn_type(11)%sym1      = "K"
       rxn_type(11)%sym2      = "K"
       rxn_type(11)%expres    =' '
       
       rxn_type(12)%type_name = 'EO2->HO2+GLYC'
       rxn_type(12)%np        = 3 
       rxn_type(12)%sym1      = "H"
       rxn_type(12)%sym2      = "H"
       rxn_type(12)%expres    =' '
       !k1*FYHORO

       rxn_type(13)%type_name = 'EO2->HO2+2CH2O'
       rxn_type(13)%np        = 3 
       rxn_type(13)%sym1      = "F"
       rxn_type(13)%sym2      = "F"
       rxn_type(13)%expres    =' '
       !k1*(1-FYHORO)

       rxn_type(14)%type_name = 'Temperature dependent branching'
       rxn_type(14)%np        = 6
       rxn_type(14)%sym1      = "V"
       rxn_type(14)%sym2      = "V"
       rxn_type(14)%expres    =' '
       !k1/(1+k2)

       rxn_type(15)%type_name = 'DMS+OH+O2'
       rxn_type(15)%np        = 3
       rxn_type(15)%sym1      = ""
       rxn_type(15)%sym2      = "G"
       rxn_type(15)%expres    =' '
       !K1/(1+K2*O2)
       
       rxn_type(16)%type_name = 'photolysis'
       rxn_type(16)%np        = 3
       rxn_type(16)%sym1      = " "
       rxn_type(16)%sym2      = " "
       rxn_type(16)%expres    =' '

       rxn_type(17)%type_name = 'O3 photolysis'
       rxn_type(17)%np        = 3
       rxn_type(17)%sym1      = "A"
       rxn_type(17)%sym2      = "Q"
       rxn_type(17)%expres    =' '

       n_tp = 17
       endsubroutine rxn_readtype

!=========================================================================
       function rxn_findtype(sym,modelid,ifphoto) result(id)
       character(len=1) :: sym
       integer          :: modelid
       integer          :: id
       integer          :: i,st,ed
       logical          :: ifphoto
       st=1
       ed=15
       if (ifphoto) then
          st=16
          ed=17
       endif
       do i = st, ed
          if (modelid.eq.1) then
             if (rxn_type(i)%sym1.eq.sym) exit
          else 
             if (rxn_type(i)%sym2.eq.sym) exit
          endif
       enddo
       id = i
       if (id.eq.ed+1) id = 1
       if (id.eq.17) id =16
       endfunction 

!=========================================================================
       subroutine rxn_init
       nr = 0
       nphoto = 0
       call rxn_readtype
       endsubroutine rxn_init

!=========================================================================
       subroutine rxn_add(rtemp,ifphoto)
       type(reaction) :: rtemp
       logical        :: ifphoto
       if (ifphoto) nphoto = nphoto + 1
       nr = nr + 1
       rxn(nr)%nreac = rtemp%nreac
       rxn(nr)%nprod = rtemp%nprod
       rxn(nr)%status = rtemp%status
       rxn(nr)%sn     = rtemp%sn
       rxn(nr)%reacs = rtemp%reacs
       rxn(nr)%prods = rtemp%prods
       rxn(nr)%prod_coefs = rtemp%prod_coefs
       rxn(nr)%r_type = rtemp%r_type
       rxn(nr)%paras = rtemp%paras
       rxn(nr)%rindex =rtemp%rindex
       endsubroutine rxn_add

!=========================================================================
       subroutine rxn_update_rates(Temp,Pres,M,O2,H2O,aer_area,aer_radius)
       integer :: i,tp,np
       real*8  :: paraF1,paraF2
       real*8  :: Temp,Pres,M,O2,H2O,aer_area,aer_radius

       paraF1=0d0
       paraF2=0d0
       do i =1, nr
          tp=rxn(i)%r_type
          if (tp.eq.0) then 
             print*,i,rxn(i)%reacs(1:rxn(i)%nreac)
             stop
          endif
          np=rxn_type(tp)%np
          if (tp.eq.5) paraF1=rate_cst(i-1)
          if (tp.eq.11) then
             paraF1=aer_area
             paraF2=aer_radius
          endif

          rate_cst(i)=rxn_rate(tp,np,rxn(i)%paras(1:np),Temp,
     +                       Pres=Pres,M=M,H2O=H2O,O2=O2,
     +                       paraI1=rxn(i)%rindex,
     +                       paraF1=paraF1,paraF2=paraF2)
       enddo

       endsubroutine rxn_update_rates

!=========================================================================
       function rxn_rate (tp, np, p, Temp, Pres, O2, N2, M, H2O,
     +                            paraI1, paraI2, paraL1, paraL2,
     +                            paraF1, paraF2, paraF3, paraF4,
     +                            nparaF5,paraF5)
     + result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       integer, intent(in)             :: tp
C       character(len=string_max),intent(out) :: expres

       real*8                          :: rate_const
 
       real*8                          :: Temp
       real*8,optional                 :: Pres, O2, N2, M,H2O
       integer, optional               :: paraI1, paraI2
       logical, optional               :: paraL1, paraL2
       real*8,optional                 :: paraF1, paraF2, paraF3,paraF4
       integer,optional                :: nparaF5
       real*8,dimension(nparaF5),optional:: paraF5
       
       select case (tp)
         case (1) 
           rate_const = r1(np, p,Temp)
         case (2)
           rate_const = r2(np, p,Temp, M)
         case (3)
           rate_const = r3(np, p,Temp, M)
         case (4)
           rate_const = r4(np, p,Temp, M)
         case (5)
           rate_const = r5(np, p,Temp,M,paraF1) !paraF1,forward rate
         case (6)
           rate_const = r6(np, p,Temp,M)
         case (7)
           rate_const = r7(np, p, Temp,M)
         case (8)
           rate_const = r8(np,p, Temp,M,H2O)
         case (9)
           rate_const = r9(np,p,Temp,O2)
         case (10)
           rate_const = r10(np,p,Temp,O2)
         case (11)
           rate_const = r11(np,p,Temp,paraF1,paraF2,M)
           !area,radius,denair
         case (12)
           rate_const = r12(np,p,Temp,M)
         case (13)
           rate_const = r13(np,p,Temp,M)
         case (14)
           rate_const = r14(np,p,Temp)
         case (15)
           rate_const = 0d0
         case (16)
              rate_const = r16(paraI1)
         case default
           rate_const = r1(np, p,Temp)           
       endselect
       
       endfunction rxn_rate

!=========================================================================
!=========================================================================
!=========================================================================
!Functions to calculate rate constants
!      Arrhenius formulation
       function r1 (np, p,Temp) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp
       real*8                            :: rate_const
       rate_const = p(1)
       if (p(2).ne.0e0) rate_const = rate_const*(3e2/Temp)**p(2)
       if (p(3).ne.0e0) rate_const = rate_const*exp(p(3)/Temp)   
       endfunction r1 

!      2.P: Pressure dependent
       function r2 (np, p, Temp,M) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,M
       real*8                            :: rate_const
       real*8                            :: kh,kl,xyrat,blog,fexp
       kh=r1(3,p(1:3),Temp)*M
       kl=r1(3,p(4:6),Temp)
       xyrat=kh/kl
       blog=log10(xyrat)
       fexp=1d0/(1d0+blog*blog)
       rate_const=kh*p(7)**fexp/(1d0+xyrat)/M
       endfunction

!      3.B:
       function r3(np,p,Temp,denair) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,denair
       real*8                            :: rate_const
       real*8                            :: k1
       k1=r1(3,p(1:3),Temp)
       rate_const=k1*(1e0-fyrno3(p(4),denair,Temp))
       endfunction

!      4.A
       function r4(np,p,Temp,denair) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,denair
       real*8                            :: rate_const
       real*8                            :: k1
       k1=r1(3,p(1:3),Temp)
       rate_const=k1*fyrno3(p(4),denair,Temp)
       endfunction

!      5.E
       function r5(np,p,Temp,M,kf) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,M,kf
       real*8                            :: rate_const
       real*8                            :: k1
       k1=r1(3,p(1:3),Temp)
       rate_const=kf*M/k1
       endfunction

!      6.X,??
       function r6(np,p,Temp,M) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,M
       real*8                            :: rate_const
       real*8                            :: k1, k2, k3
       k1=r1(3,p(1:3),Temp)
       k2=r1(3,p(4:6),Temp)
       k3=r1(3,p(7:9),Temp)*M
       rate_const=k1+k3/(1d0+k3/k2)
       endfunction

!      7.Y,??
       function r7(np,p,Temp,M) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,M!K,molec/cc
       real*8                            :: rate_const
       real*8                            :: klo1,khi1,xyrat1,fexp1,kco1
       real*8                            :: klo2,khi2,xyrat2,fexp2,kco2
       real*8                            :: blog1,blog2
!       real*8,parameter                  :: cst=0.6d0*9.871d-7
!       k1=r1(3,p(1:3),Temp)
       !rate_const=k1*(1d0+cst*Pres) Old formulation
       !update following JPL2006
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
       endfunction

!      8.Z
       function r8(np,p,Temp,M,H2O) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,H2O,M!K,molec/cc
       real*8                            :: rate_const
       real*8                            :: k1,k2
       k1=r1(3,p(1:3),Temp)
       k2=r1(3,p(4:6),Temp)
       rate_const=(k1+k2*M)*(1d0+1.4d-21*H2O*exp(2200d0/Temp))
       endfunction 

!      9.C
       function r9(np,p,Temp,O2) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,O2
       real*8                            :: rate_const
       real*8                            :: k1
       k1=r1(3,p(1:3),Temp)
       rate_const=k1*(O2+3.5d18)/(2d0*O2+3.5d18)
       endfunction

!      10.D
       function r10(np,p,Temp,O2) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,O2
       real*8                            :: rate_const
       real*8                            :: k1
       k1=r1(3,p(1:3),Temp)
       rate_const=k1*O2/(2d0*O2+3.5d18)
       endfunction

!      11.K
       function r11(np,p,Temp,area,radius,denair)result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,area,radius,denair
       real*8                            :: rate_const
       real*8                            :: stk,stkcf,sqm
       sqm=sqrt(p(1))
       stkcf=p(2)
       stk=sqrt(Temp)
       rate_const=arsl1k(area,radius,denair,stkcf,stk,sqm)
       endfunction

!      12.H 
       function r12(np,p,Temp,denair) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,denair
       real*8                            :: rate_const
       real*8                            :: k1
       k1=r1(3,p(1:3),Temp)
       rate_const=k1*fyhoro(denair,Temp)
       endfunction

!      13.F
       function r13(np,p,Temp,denair) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp,denair
       real*8                            :: rate_const
       real*8                            :: k1
       k1=r1(3,p(1:3),Temp) 
       rate_const=k1*(1d0-fyhoro(denair,Temp))
       endfunction

!      14.V
       function r14(np,p,Temp) result(rate_const)
       integer, intent(in)             :: np
       real*8, dimension(np), intent(in) :: p
       real*8, intent(in)                :: Temp
       real*8                            :: rate_const
       real*8                            :: k1,k2
       k1=r1(3,p(1:3),Temp)
       k2=r1(3,p(4:6),Temp)
       rate_const=k1/(1d0+k2)
       endfunction

!      16.photolysis
       function r16(ind) result(rate_const)
       use module_photolysis,only: jrate
       integer, intent(in)             :: ind
c       real*8, dimension(np), intent(in) :: p
c       real*8, intent(in)                :: Temp,pres
c       integer,intent(in)                :: nbin
c       real*8,dimension(nbin),intent(in) :: fff
c       integer,optional                  :: iday
       real*8                            :: rate_const
       rate_const=jrate(ind)
       endfunction
!============================================================================
!============================================================================
!============================================================================
!Helper functions for special rates
C $Id: fyrno3.f,v 4.18 2001/08/31 15:17:53 bmy v4.18 $
      function fyrno3(xcarbn,zdnum,tt)result(rfyrno3)
C***************************************************************************C
C          Returns organic nitrate yields YN = RKA/(RKA+RKB)
C          from RO2+NO reactions as a function
C          of the number N of carbon atoms.
C          Updated following Atkinson 1990.
C***************************************************************************C
C
      real*8 xcarbn
      real*8 zdnum,tt
      real*8 yyyn,xxyn,aaa,rarb,zzyn,xf,alpha,y300,beta,xminf,xm0
      real*8 rfyrno3

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

      endfunction

! $Id: fyhoro.f,v 1.1.1.1 2005/02/06 17:26:09 tmf Exp $
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
      real*8, intent(in) :: zdnum, tt
      real*8             :: rfyhoro
      ! local variables
      real*8             :: k1, k2, o2dnum
      o2dnum = zdnum * 0.21d0
      k1     = 6.0d-14 * exp(-550.d0/tt) * o2dnum
      k2     = 9.5d+13 * exp(-5988.d0/tt)
      rfyhoro = k1 / (k1 + k2)
      endfunction

      function arsl1k(area, radius, denair, stkcf, stk, sqm)result(rarsl1k)
C***************************************************************************C
C*      This function calculates the 1st-order loss rate of species        *C
C*         on wet aerosol surface, with formula from Dentener's            *C
C*         thesis, pp14.               jyl, July 1, 1994                   *C
C*         arsl1k (/s) = area / [ radius/dfkg + 4./(stkcf * xmms) ]        *C
C***************************************************************************C
      real*8 stkcf
      real*8 area, radius, stk, sqm, denair
      real*8 dfkg
      real*8 rarsl1k
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
      endfunction

      endmodule module_reaction
