      module module_chem
      implicit none
      public :: do_chem
      contains

      subroutine do_chem
      use module_model_parameter,only:DP, MAX_NRXN, MAX_NSPEC, MAX_BLK
      use module_model_option,only:option_chemmech, option_solver
      use module_domain_common,only:
     +           NIJK, grid_1_3_i, grid_1_3_j, grid_1_3_k
      use module_met_common,only:
     +           temperature, pressure, water_vapor, airdensity,
     +           aer_area, aer_radius, ifsun
      use module_chemmech_common,only:ninactrxn, inactrxn
      use module_conc_common,only:gas_conc
      use module_geoschem_rxntype,only: geos_rxnrate=>rxn_rate
      use module_ream_rxntype,only: ream_rxnrate=>rxn_rate
      use mod_smvgear_interface,only:smvgear_solve

      use module_chemmech_common,only:spec_getid
      logical, save :: firsttime=.true.
      integer,parameter :: MAX_NBLK = NIJK/MAX_BLK + 1
      integer,save      :: nblk
      integer,dimension(MAX_NBLK),save :: grid_st, blksize
      integer       :: iblk, igrid, irxn
      integer       :: gind, gi, gj, gk, g0, g1

      real(kind=DP),dimension(MAX_NRXN) :: rrate_tmp
      real(kind=DP),dimension(MAX_BLK, MAX_NRXN) :: rrate
      real(kind=DP)     :: O2, N2
      real(kind=DP),dimension(MAX_BLK, MAX_NSPEC) :: conc0, conc1
      integer       :: flag, ifsunflag

      !Initialize the loop block
      if (firsttime) then
         nblk = ceiling(float(NIJK)/MAX_BLK)
         do iblk = 1, nblk
            grid_st(iblk) = 1+ MAX_BLK * (iblk-1)
            if (iblk.ne.nblk) then
               blksize(iblk) = MAX_BLK
            else
               blksize(iblk) = NIJK - MAX_BLK * (nblk-1)
            endif
         enddo
         firsttime=.false.
      endif

!$OMP PARALLEL DO
!$OMP+DEFAULT(SHARED)
!$OMP+PRIVATE(iblk, igrid, gind, gi, gj, gk, O2, N2, irxn,
!$OMP+        rrate, rrate_tmp, conc0, conc1, g0, g1, flag,
!$OMP+        ifsunflag)
      !Loop through all grids and solve
      do iblk = 1, nblk

         !update rate constant
         ifsunflag = 2 !nighttime
         do igrid = 1, blksize(iblk)
            gind = grid_st(iblk) + igrid - 1
            gi = grid_1_3_i(gind)
            gj = grid_1_3_j(gind)
            gk = grid_1_3_k(gind)
            O2 = 0.21 * airdensity(gi,gj,gk)
            N2 = 0.78 * airdensity(gi,gj,gk)
            if (ifsun(gi,gj).eq.1) ifsunflag = 1 !daytime chem
            if (option_chemmech.eq.1) then !geoschem
               call geos_rxnrate(temperature(gi,gj,gk),
     +                           pressure(gi,gj,gk),
     +                           O2, N2, 
     +                           airdensity(gi,gj,gk),
     +                           water_vapor(gi,gj,gk),
     +                           aer_area(gi,gj,gk),
     +                           aer_radius(gi,gj,gk),
     +                           rrate_tmp)
            elseif (option_chemmech.eq.2) then !ream
               call ream_rxnrate(temperature(gi,gj,gk),
     +                           pressure(gi,gj,gk),
     +                           O2, N2, 
     +                           airdensity(gi,gj,gk),
     +                           water_vapor(gi,gj,gk),
     +                           aer_area(gi,gj,gk),
     +                           aer_radius(gi,gj,gk),
     +                           rrate_tmp)
            endif

            !update emission rate

            !update dry depotion rate

            !Multiply by concentration of inactive species for smvgear solver
            if (option_solver.eq.1) then !smvgear
               do irxn = 1, ninactrxn
                  rrate_tmp(inactrxn(1,irxn))= rrate_tmp(inactrxn(1,irxn))*
     +                                    gas_conc(gind, inactrxn(2,irxn))
               enddo
            endif

            rrate(igrid,:)=rrate_tmp(:)
            conc0(igrid,:)=gas_conc(gind,:)
         enddo!igrid

         !solve chemistry OPD
!         if (option_solver.eq.1) then !smvgear
            call smvgear_solve(1,ifsunflag, blksize(iblk), conc0, rrate,
     +                         conc1, flag)
            if (flag.ne.0) then 
               print '(A,3(X,I3))','Warning: smvgear failed near',gi,gj,gk
            endif
            g0 = grid_st(iblk)
            g1 = g0 + blksize(iblk) - 1
            gas_conc(g0:g1,:) = conc1(1:blksize(iblk),:)
!         endif
      enddo!iblk
!$OMP END PARALLEL DO
      !print*,'chemistry finished'
      print*,gas_conc(1,spec_getid('O3')),gas_conc(1,spec_getid('NO2'))
      endsubroutine do_chem

      endmodule module_chem
