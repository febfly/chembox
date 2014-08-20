      module module_reaction
      use module_model_parameter, only :: DP, MAX_NRXN, MAX_NRXNTYPE, 
     +                                    MAX_NREAC, MAX_NPROD
      implicit none
       integer, parameter              :: rxn_max  = MAX_NRXN
       integer, parameter              :: reac_max = MAX_NREAC
       integer, parameter              :: prod_max = MAX_NPROD

       type reaction
         integer :: nreac
         integer :: nprod
         character(len=1) :: status
         integer :: sn
         integer,dimension(reac_max) :: reacs
         integer,dimension(prod_max) :: prods
         real(kind=DP) ,dimension(prod_max) :: prod_coefs
         integer :: r_type
         real(kind=DP), dimension(para_max)   :: paras !paras for rate calc
         integer :: rindex          !for photolysis reactions
       endtype reaction

       type(reaction),dimension(rxn_max) :: rxn
       integer                           :: nr
       integer                           :: nphoto
       real(kind=DP),dimension(rxn_max)         :: rate_cst

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

      endmodule module_reaction
