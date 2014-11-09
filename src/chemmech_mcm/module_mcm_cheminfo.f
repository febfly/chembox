      module module_mcm_cheminfo
      use module_model_parameter,only:
     +    MAX_NSPEC,MAX_NRXN,MAX_STR1,
     +    MAX_NREAC,MAX_NPROD
      implicit none
      integer, parameter :: SL=15
      integer, parameter :: MNS=3487
      integer, parameter :: MNR=20000
      integer, parameter :: MNRT=3
      integer, parameter :: MNPD=6

      public :: mcm_read

      contains

      subroutine mcm_read
      use module_chemmech_common,only:
     +  nspec, ninactive, nactive, spec_name, spec_defconc,spec_status,
     +  nreac, nprod, reac_id, prod_id, prod_coefs, ninactrxn,
     +  nemisrxn, ndeprxn, nphotorxn, nrxn

      character (len=70) :: text
      integer            :: mcm_ns, mcm_nr
      character (len=SL),dimension(MNS) :: mcm_specname 
      integer,dimension(MNRT, MNR)      :: mcm_reac
      integer,dimension(MNR)            :: mcm_nreac
      integer,dimension(MNPD, MNR)      :: mcm_prod
      integer,dimension(MNR)            :: mcm_nprod      

      character (len=SL) :: nametmp
      integer            :: id, id2, nstmp, ir
      integer            :: ifok

      mcm_specname=''
      mcm_reac=0
      mcm_nreac=0
      mcm_prod=0
      mcm_nprod=0
      mcm_ns=0

      open(unit=20,file='mcm_chem.txt')
      text=''
      do while (text.ne.'Species present in subset')
         read(20,'(a)') text
      enddo

      !read species list
      read(20,*) id,nametmp     
      ifok=0
      do while (ifok.eq.0)
         mcm_specname(id)=nametmp
         mcm_ns = mcm_ns + 1
         read(20,*,iostat=ifok) id,nametmp         
      enddo

      !skip reaction rate calc
      text=''
      do while (text.ne.
     +   ' Subset reactants file - reaction number, species number')
         read(20,'(a)') text
      enddo

      !read reactants
      read(20,*) nstmp, mcm_nr
      if (nstmp.ne.mcm_ns) then 
         print*,'Error: in reading mcm file',nstmp,mcm_ns
         stop
      endif

      read(20,*) id, id2
      ifok=0
      do while (ifok.eq.0)
         mcm_nreac(id)=mcm_nreac(id)+1
         if (mcm_nreac(id).gt.MNRT) then
            print*,'Error: MNRT is too small,',MNRT
            stop
         endif
         mcm_reac(mcm_nreac(id),id)=id2
         read(20,*,iostat=ifok) id, id2
      enddo


      !read products
      text=''
      do while (text.ne.
     +   ' Subset products file - reaction number, species number')
         read(20,'(a)') text
      enddo
      read(20,*) id, id2
      ifok=0
      do while (ifok.eq.0)
         mcm_nprod(id)=mcm_nprod(id)+1
         if (mcm_nprod(id).gt.MNPD) then
            print*,'Error: MNPD is too small,',MNPD,id
            stop
         endif
         mcm_prod(mcm_nprod(id),id)=id2
         read(20,*,iostat=ifok) id, id2
      enddo
      close(20)


      !set variables in chemmech_common
      nspec       =mcm_ns
      ninactive   =0
      nactive     =mcm_ns-ninactive
      nrxn        =mcm_nr

      spec_name   =mcm_specname
      spec_status(1:nspec)='A'     !all species are active
      spec_defconc(1:nspec)=1d-18  !mixing ratio;default value

      nreac       = mcm_nreac
      nprod       = mcm_nprod
      reac_id     = mcm_reac
      prod_id     = mcm_prod

      do ir = 1, nrxn
         prod_coefs(1:nprod(ir),ir)=1. !all reactions are atomic
      enddo

      ninactrxn = 0
      nemisrxn  = 0
      ndeprxn   = 0
      nphotorxn = 0 !no need for special treatment for photolysis

      endsubroutine mcm_read

!================================================================
      endmodule module_mcm_cheminfo 
