      module module_spec
       implicit none
       !variables
       integer, parameter, private              :: ns_max = 200
       integer, parameter, private              :: str_len = 15
       type species
            character(len=str_len)              :: name
            character(len=1)                    :: status
            real                                :: default_conc
       endtype species

       type(species),dimension(ns_max) :: spec
       integer                                  :: ns      
       integer                                  :: ninactive,nactive 
 
       !functions
       public                             :: spec_init
       public                             :: spec_add
       public                             :: spec_finish_add
       public                             :: spec_getid
       private                            :: treat_str
       private                            :: spec_ifduplicate

      contains
!=========================================================================
       subroutine spec_init
       ns = 0
       nactive = 0
       ninactive = 0
       endsubroutine

!=========================================================================
       function spec_add(s,stat,conc) result (id)
       character(len=*),intent(in)        :: s
       character(len=1), optional         :: stat
       real   , optional                  :: conc

       logical                            :: ifdup
       integer                            :: id
       character(len=1)                   :: stat1
       real                               :: conc1
       character(len=str_len)             :: s1

       if (present(stat)) then 
           stat1 = stat 
       else
           stat1 = 'A'
       endif
        
       if (present(conc)) then
          conc1 = conc 
       else 
          conc1 = 1e-20
       endif

       s1 = treat_str(s)
c       id = spec_getid(s1)
       ifdup = spec_ifduplicate(s)
       if ( ifdup.eq..false. ) then
          ns = ns + 1
          if (ns.gt.ns_max) then
             print*,'increase ns_max',ns_max
             stop
          endif
          if (stat1.eq.'A') then !active species
             nactive = nactive + 1
             spec(nactive)%name = s1
             spec(nactive)%status = stat1
             spec(nactive)%default_conc = conc1
             id = nactive
          else if (stat1.eq.'I'.or.stat1.eq.'D') then !inactive species
             ninactive = ninactive + 1
             spec(ns_max-ninactive+1)%name = s1
             spec(ns_max-ninactive+1)%status = stat1
             spec(ns_max-ninactive+1)%default_conc = conc1
             id = ns_max-ninactive+1
          else
             print*,'cannot recognize stat1',stat1
             stop
          endif
       endif
       
       endfunction spec_add

!=========================================================================
       subroutine spec_finish_add
       if (ninactive.gt.0) then 
          spec(nactive+1:ns) = spec(ns_max-ninactive+1:ns_max)
       endif
       endsubroutine spec_finish_add
       
!=========================================================================
!      private function, only useful in the spec_add stage
       function spec_ifduplicate (s) result(ifd)
       character(len=*),intent(in)        :: s
       logical                            :: ifd
       integer                            :: n 
       character(len=str_len)             :: s1

       s1 = treat_str(s)
       ifd=.false.
       do n = nactive,1,-1
          if (s1.eq.spec(n)%name) then
              ifd=.true.
              exit
          endif
       enddo
       if (ifd) return
       do n = ns_max-ninactive+1, ns_max 
          if (s1.eq.spec(n)%name) then
              ifd=.true. 
              exit 
          endif 
       enddo
       endfunction

!=========================================================================
       function spec_getid (s) result (id)
       character(len=*),intent(in)        :: s
       integer                            :: id
       integer                            :: n
       character(len=str_len)             :: s1

       s1 = treat_str(s)
       do n = ns, 1, -1
          if (s1 .eq. spec(n)%name) exit
       enddo
       id = n 

       endfunction spec_getid

!=========================================================================
       function treat_str (s) result (s1)
       character(len=*), intent(in)   :: s
       character(len=str_len)         :: s1
       integer                        :: ll
       ll = len(s)
       s1 = ""
       if (ll.le.str_len) then
          s1(1:ll) = s(1:ll)
       else 
          s1(1:str_len) = s(1:str_len)
       endif
       endfunction treat_str

      endmodule module_spec
