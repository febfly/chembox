      subroutine itf_chemmech_setup
      !global modules
      use module_model_parameter, only : DP
      use module_model_cheminfo
      !use module_model_path

      !local modules
      use module_spec
      use module_reaction
      use module_ream, only:ream_read
      use module_photolysis, only:photolysis_read
      implicit none
      integer :: i
      character(len=255) :: pream,fream,fratj,fjspe

      pream='/data12/yzhang425/REAM/REAMSource/v2010_julyo3/Run/'
      fream='chem.dat'
      fratj='ratj.d'
      fjspe='jv_spec.dat'

      !Read input files
      call photolysis_read(trim(pream)//trim(fratj),trim(pream)//trim(fjspe))
      call ream_read(trim(pream)//trim(fream))

      !Initialize variables in the global module "module_model_cheminfo"
      spec_name(:)         = ''
      spec_flag(:)         = ''
      spec_default_conc(:) = 0.
      rxn_nreac(:)         = 0
      rxn_nprod(:)         = 0
      rxn_reac (:,:)       = 0
      rxn_prod (:,:)       = 0
      rxn_reac_coef(:,:)   = 0.
      rxn_prod_coef(:,:)   = 0.
      rxn_flag (:)         = ''

      !Set data in global module
      nspec_all      = ns
      nspec_active   = nactive
      nspec_inactive = ninactive
      spec_name(1:nspec_all) = spec(1:ns)%name
      spec_flag(1:nspec_all) = spec(1:ns)%status
      spec_default_conc(1(1:nspec_all) = spec(1:ns)%default_conc
    
      nrxn_all   = nr
      nrxn_photo = nphoto
      rxn_nreac(1:nrxn_all)  = rxn(1:nr)%nreac
      rxn_nprod(1:nrxn_all)  = rxn(1:nr)%nprod
      rxn_flag (1:nrxn_all)  = rxn(1:nr)%status      

      do i = 1, nrxn_all
         rxn_reac(:,i) = rxn(i)%reacs(:)
         rxn_prod(:,i) = rxn(i)%prods(:)
         rxn_reac_coef(1:rxn(i)%nreac,i) = 1.
         rxn_prod_coef(:,i)              = rxn(i)%prod_coefs(:)
      enddo
      
      endsubroutine itf_chemmech_setup
