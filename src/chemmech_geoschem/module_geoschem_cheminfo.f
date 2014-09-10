!=========================================================================
! Module   : module_geoschem_cheminfo
!            1. stores the information about chemical species and reactions
!            2. defines subroutines that initialize and add the data
! Written by Yuzhong Zhang, 9/8/2014
!=========================================================================
      module module_geoschem_cheminfo
      use module_model_parameter, only : DP, MAX_NSPEC, MAX_STR1,
     +    MAX_NRXN, MAX_NREAC, MAX_NPROD
      use module_geoschem_rxntype, only : MAX_NPARA
      implicit none

      !species information
      integer :: ns, ninactive, nactive
      character(len=MAX_STR1),dimension(MAX_NSPEC) :: specname
      character(len=1),dimension(MAX_NSPEC)        :: status
      real(kind=DP),dimension(MAX_NSPEC)           :: def_conc

      !reaction information
      integer :: nr, nphoto
      integer, dimension(MAX_NRXN)     :: nreac, nprod
      character(len=1),dimension(MAX_NSPEC)        :: rxn_symbol
      integer, dimension(MAX_NREAC,MAX_NRXN)       :: reacs
      integer, dimension(MAX_NPROD,MAX_NRXN)       :: prods
      real(kind=DP),dimension(MAX_NPROD,MAX_NRXN)  :: prod_coefs
      integer, dimension(MAX_NRXN)                 :: r_type
      real(kind=DP),dimension(MAX_NPARA,MAX_NRXN)  :: paras

      !public functions
      public :: geoschem_cheminfo_init
      public :: spec_add
      public :: spec_finish_add
      public :: spec_getid
      public :: rxn_add

      !private functions
      private :: spec_ifduplicate
      private :: treat_str

      contains

!=========================================================================
      subroutine geoschem_cheminfo_init
      ns = 0
      nactive = 0
      ninactive = 0
      specname(:) = ''
      status(:) = ''
      def_conc(:) = 0.

      nr = 0
      nphoto = 0
      nreac(:) = 0
      nprod(:) = 0
      rxn_symbol(:) = ''
      reacs(:,:) = 0
      prods(:,:) = 0
      prod_coefs(:,:) = 0.
      r_type(:) = 0
      paras(:,:) = 0.
      endsubroutine geoschem_cheminfo_init

!=========================================================================

      subroutine spec_add(s,stat,conc)
      character(len=*),intent(in)        :: s
      character(len=1)                   :: stat
      real(kind=DP)                      :: conc
      
      integer                            :: id
      logical                            :: ifdup
      character(len=MAX_STR1)            :: s1

      s1=treat_str(s)
      ifdup=spec_ifduplicate(s1)

      if (ifdup) then
         print*,'Error:duplicate definition of species'
         stop
      endif

      ns = ns + 1
      if (stat.eq.'A') then
         nactive = nactive + 1
         specname(nactive) = s1
         status(nactive)   = stat
         def_conc(nactive) = conc
      elseif (stat.eq.'I'.or.stat.eq.'D') then
         ninactive = ninactive + 1
         id = MAX_NSPEC-ninactive+1
         specname(id) = s1
         status(id)   = stat
         def_conc(id) = conc
      else
         print*,'Error:cannot recognize stat',stat
         stop
      endif
      endsubroutine spec_add

!=========================================================================
       subroutine spec_finish_add
       if (ninactive.gt.0) then 
          specname(nactive+1:ns) = specname(MAX_NSPEC-ninactive+1:MAX_NSPEC)
          status(nactive+1:ns)   = status(MAX_NSPEC-ninactive+1:MAX_NSPEC)
          def_conc(nactive+1:ns) = def_conc(MAX_NSPEC-ninactive+1:MAX_NSPEC)
       endif
       endsubroutine spec_finish_add
       
!=========================================================================

       function spec_getid (s) result (id)
       character(len=*),intent(in)        :: s
       integer                            :: id
       integer                            :: n
       character(len=MAX_STR1)             :: s1

       s1 = treat_str(s)
       do n = ns, 1, -1
          if (s1 .eq. specname(n)) exit
       enddo
       id = n 
       endfunction spec_getid

!=========================================================================
!      private function, only useful in the spec_add stage
       function spec_ifduplicate (s) result(ifd)
       character(len=MAX_STR1),intent(in)        :: s
       logical                                   :: ifd
       integer                                   :: n 

       ifd=.false.
       do n = nactive,1,-1
          if (s.eq.specname(n)) then
              ifd=.true.
              exit
          endif
       enddo
       if (ifd) return
       do n = MAX_NSPEC-ninactive+1, MAX_NSPEC
          if (s.eq.specname(n)) then
              ifd=.true. 
              exit 
          endif 
       enddo
       endfunction spec_ifduplicate
!=========================================================================
       subroutine rxn_add(ifp,n1, n2, id1, id2, cf2,  typeid, pr)
       integer, intent(in) :: ifp
       integer, intent(in) :: n1, n2
       integer, dimension(MAX_NREAC),intent(in) :: id1
       integer, dimension(MAX_NPROD),intent(in) :: id2
       real(kind=DP), dimension(MAX_NPROD),intent(in) :: cf2
!       character(len=1),intent(in)              :: flag
       integer,intent(in)                       :: typeid
       real(kind=DP), dimension(MAX_NPARA),intent(in) :: pr       

!       if (flag.eq.'D') return !dead reactions
       nr = nr + 1
       if (ifp.eq.1) nphoto = nphoto + 1
       nreac(nr) = n1
       nprod(nr) = n2
       reacs(:,nr) = id1(:)
       prods(:,nr) = id2(:)
       prod_coefs(:,nr) = cf2(:)
!       rxn_symbol(nr) = flag
       r_type (nr) = typeid
       paras(:,nr) = pr(:)

       endsubroutine rxn_add
!=========================================================================
       function treat_str (s) result (s1)
     
       character(len=*), intent(in)   :: s
       character(len=MAX_STR1)        :: s1
       integer                        :: ll
       ll = len(s)
       s1 = ""
       if (ll.le.MAX_STR1) then
          s1(1:ll) = s(1:ll)
       else 
          s1(1:MAX_STR1) = s(1:MAX_STR1)
       endif
       endfunction treat_str
      endmodule module_geoschem_cheminfo
