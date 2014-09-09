      subroutine itf_smvgear_setup(ncs,chemintv,n_spec, n_active, specname_len,
     +           specname,  n_rxn, n_photo, n_mxreac, n_mxprod, reac,
     +           prod, reac_coef, prod_coef)
      use mod_comode, only : DP, STRLEN, IGAS, IPHOT, NMRATE, IPHOT,
     +                       NMREAC, NMPROD, iout
      real(kind=DP),intent(in) :: chemintv
      integer, intent(in) :: ncs, n_spec, n_active, specname_len
      integer, intent(in) :: n_rxn, n_photo, n_mxreac, n_mxprod
      character(len=specname_len),dimension(n_spec) :: specname
      integer, dimension(n_mxreac,n_rxn),intent(in) :: reac
      integer, dimension(n_mxprod,n_rxn),intent(in) :: prod
      real(kind=DP),dimension(n_mxreac, n_rxn), intent(in) :: reac_coef
      real(kind=DP),dimension(n_mxprod, n_rxn), intent(in) :: prod_coef

      character(len=STRLEN) :: typename
          
      iout = 6
      typename='CHEM'

      !Check dimensions
      if (ncs.gt.ICS) then 
         write(iout,*) 'Error: ncs is greater then ICS',
     +                 ncs, ICS
         stop
      endif

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

      endsubroutine itf_smvgear_setup
