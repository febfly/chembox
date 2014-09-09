C
C File:   smvgear_setchem.f
C Author: Yuzhong
C
C Created on July 25, 2014, 2:04 PM
C
      subroutine smvgear_setchem(ncs, n_spec, n_active, specname_len, specname, 
     +   n_rxn,  n_photo, n_mxreac, n_mxprod, reac, prod, reac_coef, prod_coef, chemtypename)
      use mod_comode
      implicit none
      !Input variables
      integer, intent(in) :: ncs
      integer, intent(in) :: n_spec, n_active, n_rxn, n_photo
      integer, intent(in) :: specname_len, n_mxreac, n_mxprod

      character(len=STRLEN), dimension(n_spec), intent(in) :: specname
      !integer,                     dimension(n_spec), intent(in) :: spectype

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

