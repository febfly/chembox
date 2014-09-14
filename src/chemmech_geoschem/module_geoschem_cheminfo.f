!=========================================================================
! Module   : module_geoschem_cheminfo
!            1. stores the information about chemical species and reactions
!            2. defines subroutines that initialize and add the data
! Written by Yuzhong Zhang, 9/8/2014
!=========================================================================
      module module_geoschem_cheminfo
      !parameters
      use module_chemmech_common, only : DP, MAX_NSPEC, MAX_STR1,
     +    MAX_NRXN, MAX_NREAC, MAX_NPROD, MAX_NINACTRXN, MAX_NEMISRXN,
     +    MAX_NDEPRXN

      !species information
      use module_chemmech_common, only : 
     +           ns=>nspec, ninactive, nactive, specname=>spec_name,
     +           status=>spec_status, def_conc=>spec_defconc, 
     +           nphoto=>nphotorxn

      !reaction information
      use module_chemmech_common, only : 
     +           nr=>nrxn, nreac, nprod, reacs=>reac_id, prods=>prod_id,
     +           prod_coefs

      !special reaction index (Emission, drydep, photolysis, etc)
      use module_chemmech_common, only :
     +           ninactrxn, nemisrxn, ndeprxn,nphotorxn,
     +           inactrxn,emisrxn,deprxn,photorxn

      !utility function
      use module_chemmech_common, only : spec_getid

      !reaction rate type and parameters
      use module_geoschem_common, only : MAX_NPARA, r_type, paras

      implicit none

      !public functions
      public :: cheminfo_init
      public :: spec_add
      public :: spec_finish_add
      public :: rxn_add
      public :: rxn_finish_add

      !private functions
      private :: spec_ifduplicate
      private :: treat_str

      contains

!=========================================================================
      subroutine cheminfo_init
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
      reacs(:,:) = 0
      prods(:,:) = 0
      prod_coefs(:,:) = 0.
      r_type(:) = 0
      paras(:,:) = 0.
      endsubroutine cheminfo_init

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
      if (ns.gt.MAX_NSPEC) then
         print*,'Error:parameter MAX_NSPEC is too small'
         stop
      endif

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
       !Move inactive and dead species to the end of the list
       if (ninactive.gt.0) then 
          specname(nactive+1:ns) = specname(MAX_NSPEC-ninactive+1:MAX_NSPEC)
          status(nactive+1:ns)   = status(MAX_NSPEC-ninactive+1:MAX_NSPEC)
          def_conc(nactive+1:ns) = def_conc(MAX_NSPEC-ninactive+1:MAX_NSPEC)
       endif
       endsubroutine spec_finish_add
       
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
       if (nr.gt.MAX_NRXN) then
          print*,'Error: parameter MAX_NRXN is too small'
          stop
       endif

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
       subroutine rxn_finish_add
       integer :: ir, is, sid
       ninactrxn=0
       nemisrxn =0
       ndeprxn  =0
       do ir=1,nr
          do is=1,nreac(ir)
             !reactant is 'EMISSION',i.e. a emission rxn
             sid = reacs(is,ir)
             if (specname(sid).eq.'EMISSION') then
                nemisrxn = nemisrxn + 1
                if (nemisrxn.gt.MAX_NEMISRXN) then
                   print*,'Error: MAX_NEMISRXN is too small'
                   stop
                endif
                emisrxn(nemisrxn) = ir
             !reactant is an inactive species
             elseif (status(sid).eq.'I') then
                ninactrxn = ninactrxn + 1
                if (ninactrxn.gt.MAX_NINACTRXN) then
                   print*,'Error: MAX_NINACTRXN is too small'
                   stop
                endif
                inactrxn(1,ninactrxn) = ir
                inactrxn(2,ninactrxn) = sid
             endif
          enddo

          !a product is 'DRYDEP',i.e. a drydep rxn
          do is=1,nprod(ir)
             sid = prods(is,ir)
             if (specname(sid).eq.'DRYDEP') then
                ndeprxn  = ndeprxn  + 1
                if (ndeprxn.gt.MAX_NDEPRXN) then
                   print*,'Error: MAX_NDEPRXN is too small'
                   stop
                endif
                deprxn(ndeprxn)   = ir
             endif
          enddo
       enddo
       endsubroutine rxn_finish_add
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
