      module module_analyze
      implicit none
      integer,parameter :: nmax=500
      integer,parameter :: nr_max=100
      integer,parameter :: ns_max=200

c     analyze_do1
      integer,dimension(ns_max,nr_max) :: table_s2r
      integer,dimension(ns_max)        :: ct_s2r
c--------------------------------------------------
c      analyze_do2
      integer,dimension(ns_max*nr_max,2) :: r2p
      integer,dimension(ns_max)          :: pos_r2p
      integer,dimension(ns_max*nr_max,2) :: p2r
      integer,dimension(ns_max)          :: pos_p2r
   
c--------------------------------------------------

      integer,parameter :: loopmax=110
     
      contains
      subroutine analyze_do
c      use module_spec,only:spec,ns
      use module_reaction,only:rxn,nr,nphoto
      integer :: ir, is, inds
      ct_s2r(:) = 0
      table_s2r(:,:)=0
      do ir=1,nr
         do is=1,rxn(ir)%nreac
c            print*,ir,is,rxn(ir)%reacs(is)
            inds = rxn(ir)%reacs(is)
            ct_s2r(inds) = ct_s2r(inds) + 1
            table_s2r(inds,ct_s2r(inds)) = ir
         enddo
      enddo
      
      endsubroutine analyze_do

c      subroutine analyze_do2

c      endsubroutine analyze_do2

c      function analyze_genpath(a,b,p,np) result(s)
c      integer :: a, b, s,np
c      integer,dimension(nmax) :: p
c      integer :: is, ie, i

c      is=pos_r2p(a)
c      ie=pos_r2p(a+1)-1      
c      do i = is, ie
c         if (r2p(i,1).eq.b) then
c            np = np+1
c            p(np) = r2p(i,2)
c            s = 1
c            return
c         endif
c      enddo

c      endfunction analyze_genpath

      subroutine analyze_cascade(tc,ec,ne,oc,no)
      use module_spec,only: spec_getid,spec,ns
      
      character(len=15) :: tc
      integer           :: t
      integer           :: ne
      character(len=15),dimension(ne) :: ec
      integer,dimension(nmax)::e
      integer           :: no
      character(len=15),dimension(nmax) :: oc
      integer,dimension(nmax)::o

      integer :: i
      t=spec_getid(tc)
      do i=1,ne
         e(i)=spec_getid(ec(i))
      enddo
      call analyze_cascade_a(t,e,ne,o,no)
      do i=1,no
         oc(i)=spec(o(i))%name
      enddo

      endsubroutine analyze_cascade

      subroutine analyze_cascade_a(t,e,ne,o,no)
      use module_spec,only:spec,ns
      use module_reaction,only:rxn,nr,nphoto
      integer,intent(in) :: t
      integer,intent(in) :: ne
      integer,dimension(nmax),intent(in) :: e
      integer,intent(out) :: no
      integer,dimension(nmax),intent(out) :: o
      
      integer :: iloop,ie,ip,is,ir,it,nos,ife
      iloop=0

      o(1)=t
      no=1
      nos=1

      do while (iloop.le.loopmax.and.nos.le.no)
         is = o(nos)
         do it = 1, ct_s2r(is)
            ir = table_s2r(is,it)
            do ip = 1, rxn(ir)%nprod
               ife=0
               do ie = 1, ne
                  if (e(ie).eq.rxn(ir)%prods(ip)) then
                     ife=1
                     exit
                  endif
               enddo
              if (ife.eq. 0) call analyze_addlist(o,no,rxn(ir)%prods(ip))
            enddo
         enddo
         nos = nos + 1
         iloop = iloop + 1
      enddo

      endsubroutine analyze_cascade_a

      subroutine analyze_addlist(list,nl,e)
      integer :: nl
      integer,dimension(nmax) :: list
      integer :: e
      integer :: i,ifadd

      ifadd=1
      do i=nl,1,-1
         if (list(i).eq.e) then
             ifadd=0
             exit
         endif
      enddo
      if (ifadd.eq.1) then
         nl = nl + 1
         list(nl) = e
      endif
 
      endsubroutine analyze_addlist

      endmodule module_analyze
