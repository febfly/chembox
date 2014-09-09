      subroutine subfun(ncsp)
!
      use mod_comode
      implicit none
      !input variables
      integer, intent(in) :: ncsp
      !local variables
      integer :: nkn,ja,jb,jc,nh,k,nk2,nh2,jspc,npl,nl5,nh5,nl4,nh4,nl3,nh3
      integer :: nl2,nl1,nh1,nc,nk0,nk1,nk3,nk4,n
      integer :: nk,i,jnew,kloop
      real(kind=DP) :: concmult,fracn
!
! *********************************************************************
! ************        written by mark jacobson (1993)      ************
! ***             (c) copyright, 1993 by mark z. jacobson           *** 
! ***       u.s. copyright office registration no. txu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!      SSSSSSS  U     U  BBBBBBB  FFFFFFF  U     U  N     N  
!      S        U     U  B     B  F        U     U  NN    N 
!      SSSSSSS  U     U  BBBBBBB  FFF      U     U  N  N  N 
!            S  U     U  B     B  F        U     U  N    NN
!      SSSSSSS  UUUUUUU  BBBBBBB  F        UUUUUUU  N     N
!
! *********************************************************************
! *  this subroutine evaluates the first derivative of each ordinary  *  
! *                  differential equation (ode)                      * 
! *                                                                   *
! * how to call subroutine:                                           *
! * ----------------------                                            *
! *  call subfun.f from smvgear.f with                                * 
! *     ncs  = 1..ncsgas for gas chemistry                            *
! *     ncsp = ncs        for daytime   gas chem                      *  
! *     ncsp = ncs   +ics for nighttime gas chem                      *  
! *********************************************************************
!
! example
! -------
!
! species:         a,   b,   c
! concentrations: [a], [b], [c]
!
! reactions:    1) a           --> b      j 
!               2) a  + b      --> c      k1 
!               3  a  + b + c  --> d      k2  
!
! first         d[a] / dt  =  -j[a] - k1[a][b] - k2[a][b][c]
! derivatives:  d[b] / dt  =  +j[a] - k1[a][b] - k2[a][b][c]
!               d[c] / dt  =        + k1[a][b] - k2[a][b][c]
!               d[d] / dt  =                   + k2[a][b][c]
!
! *********************************************************************
!
! concmult  = product of concentrations in a rate. if two  
!             consecutive reactions have the same species reacting
!             (eg a + b --> c and a + b --> d + e) then use the 
!             same value of concmult for both reactions.
! cnew      = init (and final) species conc (# cm-3-air or moles l-1-h2o)
! gloss     = first derivative = sum of prod. minus loss rates for a species
! irma,b,c  = locates reordered active species numbers  
! ischan    = number of odes.
! lossra..  = reaordered reaction rate numbers for each loss (and prod) term
! ktloop    = number of grid-cells in a grid-block
! nsubfun   = counts the number of times this routine is called
! rrate     = forward rate coefficient 
!           = s-1                                 for rates with 1 reactant
!           = l-h2o mole-1 s-1  or cm**3 #-1 s-1  for rates with 2 reactants 
!           = l**2-h2o m-2 s-1  or cm**6 #-2 s-1  for rates with 3 reactants 
! trate     = reaction rate  moles l-1 -h2o s-1 or # cm-3 s-1 
! 
! *********************************************************************
! *                      set rates of reaction                        *
! *********************************************************************
!
!
!      print*,'check subfun1'
!      nsubfun        = nsubfun + 1
      nfdh1          = nfdh2 + ioner(ncsp) 
!
! *********************************************************************
! *     first derivatives for rates with three active loss terms      *
! *********************************************************************
!
      do     nkn     = 1, nfdh3  
       ja            = irma(nkn)
       jb            = irmb(nkn)
       jc            = irmc(nkn)
       nh            = nkn + nallr
       do     k      = 1, ktloop
        trate(k,nkn) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jb)*cnew(k,jc) 
        trate(k,nh) = -trate(k,nkn) 
       enddo
      enddo
!
! *********************************************************************
! *     first derivatives for rates with two active loss terms        *
! *********************************************************************
!
      do     nkn     = nfdl2, nfdrep
       ja            = irma(nkn)
       jb            = irmb(nkn)
       nh            = nkn + nallr
       do     k      = 1, ktloop
        trate(k,nkn) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jb) 
        trate(k,nh) = -trate(k,nkn) 
       enddo
      enddo
!
! *********************************************************************
! *     first derivatives for rates with two active loss terms and    *
! *     where the subsequent reaction has the same reactants but a    *
! *     different rate.                                               *
! *********************************************************************
!
      do     nkn     = nfdrep1, nfdh2, 2
       ja            = irma(nkn)
       jb            = irmb(nkn)
       nk2           = nkn + 1
       nh            = nkn + nallr
       nh2           = nk2 + nallr
       do     k      = 1, ktloop
        concmult     = cnew(k,ja)   * cnew(k,jb) 
        trate(k,nkn) = rrate(k,nkn) * concmult
        trate(k,nk2) = rrate(k,nk2) * concmult
        trate(k,nh)  = -trate(k,nkn) 
        trate(k,nh2) = -trate(k,nk2) 
       enddo
      enddo
!
! *********************************************************************
! *     first derivatives for rates with one active loss term         *
! *********************************************************************
!
      do     nkn     = nfdl1, nfdh1
       ja            = irma(nkn)
       nh            = nkn + nallr
       do     k      = 1, ktloop
        trate(k,nkn) = rrate(k,nkn) * cnew(k,ja) 
        trate(k,nh) = -trate(k,nkn) 
       enddo
      enddo
!
! *********************************************************************
! *                  initialize first derivative = 0                  *
! *********************************************************************
!
      do     jspc      = 1, ischan
       do     k        = 1, ktloop
        gloss(k,jspc)  = 0. 
       enddo
      enddo
!
! *********************************************************************
! * sum net (not reproduced) kinetic and photo gains and losses for   *
! * each species.                                                     *
! *********************************************************************
! sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!
      do     npl       = npllo(ncsp), nplhi(ncsp)
       jspc            = jspnpl(npl)
       nl5             = npl5(  npl)
       nh5             = nph5(  npl)
       nl4             = npl4(  npl)
       nh4             = nph4(  npl)
       nl3             = npl3(  npl)
       nh3             = nph3(  npl)
       nl2             = npl2(  npl)
       nh2             = nph2(  npl)
       nl1             = npl1(  npl)
       nh1             = nph1(  npl)
!
! ***********************  sum 5 terms at a time  ********************* 
!
       do     nc       = nl5, nh5
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        nk2            = lossrc(nc)
        nk3            = lossrd(nc)
        nk4            = lossre(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
     1                 - trate(k,nk1)  - trate(k,nk2)
     2                 - trate(k,nk3)  - trate(k,nk4)  
        enddo
       enddo
!
! ***********************  sum 4 terms at a time  ********************* 
!
       do     nc       = nl4, nh4 
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        nk2            = lossrc(nc)
        nk3            = lossrd(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
     1                 - trate(k,nk1)  - trate(k,nk2)
     2                 - trate(k,nk3)  
        enddo
       enddo
!
! ***********************  sum 3 terms at a time  ********************* 
!
       do     nc       = nl3, nh3  
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        nk2            = lossrc(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
     1                 - trate(k,nk1)  - trate(k,nk2)
        enddo
       enddo
!
! ***********************  sum 2 terms at a time  ********************* 
!
       do     nc       = nl2, nh2   
        nk0            = lossra(nc)
        nk1            = lossrb(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0) 
     1                 - trate(k,nk1)      
        enddo
       enddo
!
! ***********************  sum 1 term at a time  ********************** 
!
       do     nc       = nl1, nh1    
        nk0            = lossra(nc)
        do     k       = 1, ktloop
         gloss(k,jspc) = gloss(k,jspc) - trate(k,nk0)       
        enddo
       enddo
      enddo!npl
!
! *********************************************************************
! *  sum production term for reactions where products fractionated    *
! *********************************************************************
!
      do     n         = nfrlo(ncsp), nfrhi(ncsp)
       jspc            = jspcnfr(n)
       nkn             = nknfr(  n)
       fracn           = fracnfr(n)
       do     k        = 1, ktloop
        gloss(k,jspc)  = gloss(k,jspc) + fracn * trate(k,nkn)       
       enddo
      enddo
!
! *********************************************************************
! **********************  end of subroutine subfun  *******************
! *********************************************************************
!
!      print*,'check subfun2'
      return
      end
