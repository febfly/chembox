      module module_forpython
      implicit none
      contains
      subroutine fp_write(tc,oc,ec,ne)
      use module_spec,only:spec_getid
      character(len=15) :: tc,oc
      integer           :: ne
      character(len=15) :: ec(ne)
      integer           :: t, o ,i
      integer           :: e(ne)
      t=spec_getid(tc)
      o=spec_getid(oc)
      do i=1,ne
         e(i)=spec_getid(ec(i))
      enddo
      call fp_write_1(t,o,e,ne)
      endsubroutine

      subroutine fp_write_1(t, o, e, ne)
      use module_reaction,only:rxn,nr
      integer ir, i1, i2,ie,s1,s2, t, o
      integer ne
      integer e(ne)

      open(unit=40,file='chem.csv')
      write(40,'(i3)') t
      write(40,'(i3)') o
      do ir=1,nr
         do i1=1,rxn(ir)%nreac
            s1=rxn(ir)%reacs(i1)
            if (ifexcld(e,ne,s1)) cycle
            do i2=1,rxn(ir)%nprod
               s2=rxn(ir)%prods(i2)
               if (ifexcld(e,ne,s2)) cycle
               write(40,'3(i3,x)')s1,s2,ir
            enddo   
         enddo
      enddo
      close(40)
      endsubroutine fp_write_1

      function ifexcld(e,ne,t) result(ifok)
      integer ne, t, i
      integer e(ne)
      logical ifok
      ifok=.false.
      do i=1,ne
         if (t.eq.e(i)) then
            ifok=.true.
            exit
         endif
      enddo
      endfunction ifexcld
      endmodule
