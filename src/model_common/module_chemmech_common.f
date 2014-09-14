      module module_chemmech_common
      use module_model_parameter,only : MDP => DP,MMAX_NSPEC=>MAX_NSPEC,
     +                         MMAX_NRXN=>MAX_NRXN, MMAX_STR1=>MAX_STR1,
     +                      MMAX_NREAC=>MAX_NREAC,MMAX_NPROD=>MAX_NPROD      
      implicit none
      !parameters
      integer,parameter :: DP        = MDP
      integer,parameter :: MAX_NSPEC = MMAX_NSPEC
      integer,parameter :: MAX_NRXN  = MMAX_NRXN
      integer,parameter :: MAX_STR1  = MMAX_STR1
      integer,parameter :: MAX_NREAC = MMAX_NREAC
      integer,parameter :: MAX_NPROD = MMAX_NPROD

      integer,parameter :: MAX_NINACTRXN = 100
      integer,parameter :: MAX_NEMISRXN  = 30
      integer,parameter :: MAX_NDEPRXN   = 30
      integer,parameter :: MAX_NPHOTORXN = 100

      !species information
      integer :: nspec, ninactive, nactive
      character(len=MAX_STR1),dimension(MAX_NSPEC)   :: spec_name
      character(len=1),dimension(MAX_NSPEC)          :: spec_status
      real(kind=DP),dimension(MAX_NSPEC)             :: spec_defconc

      !reaction information
      integer :: nrxn
      integer, dimension(MAX_NRXN)                 :: nreac, nprod
      integer, dimension(MAX_NREAC,MAX_NRXN)       :: reac_id
      integer, dimension(MAX_NPROD,MAX_NRXN)       :: prod_id
      real(kind=DP),dimension(MAX_NPROD,MAX_NRXN)  :: prod_coefs      

      !More reaction information
      integer :: ninactrxn
      integer,dimension(2, MAX_NINACTRXN) :: inactrxn
      integer :: nemisrxn
      integer,dimension(MAX_NEMISRXN)     :: emisrxn
      integer :: ndeprxn
      integer,dimension(MAX_NDEPRXN)      :: deprxn
      integer :: nphotorxn
      integer,dimension(MAX_NPHOTORXN)    :: photorxn

      !utitility function
      public :: spec_getid
      
      contains

      !get the ID of a species given its name
       function spec_getid (s) result (id)
       character(len=*),intent(in)        :: s
       integer                            :: id
       integer                            :: n
  
       do n = nspec, 1, -1
          if (adjustl(trim(s)) .eq. spec_name(n)) exit
       enddo
       id = n 
       endfunction spec_getid

      endmodule
