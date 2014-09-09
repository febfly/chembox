!
! File:   itf_smvgear_solve.f
! Author: Yuzhong
!
! Created on July 26, 2014, 3:49 PM
!
      subroutine itf_smvgear_solve (ncs, ifsun, n_grid, n_spec, n_rxn, conc0, rate_const, conc1, flag)

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
!     conc1    : final concentration at the end of the step, [n_grid, n_spec]  
!     flag = 0 : solve succesfully
!          = 2 : solver fails. nylowdec in smvgear.f is larger than threshold
!          = 3 : solver fails. jrestar  in smvgear.f is larger than threshold
!==============================================================================
          use mod_comode
          implicit none
          !Input variables
          integer, intent(in) :: ncs, ifsun
          integer, intent(in) :: n_spec, n_rxn, n_grid
          real(kind=DP), dimension(n_grid, n_spec), intent(in) :: conc0
          real(kind=DP), dimension(n_grid, n_rxn), intent(in) :: rate_const
          
          !Output variables
          real(kind=DP), dimension(n_grid, n_spec), intent(out) :: conc1
          integer, intent(out)                          :: flag
          !character(len=MSGLEN), intent(out)            :: message
                 
          !Local variables
          real(kind=DP):: chemintv
          integer      :: nkn, nk, nh, kloop
          integer      :: jnew, jold
         
!         print*,'in solve' 
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
                         
      endsubroutine itf_smvgear_solve
