      module module_kpp
      implicit none
      public  :: kpp_write
      private :: kpp_def
      private :: kpp_spc
      private :: kpp_rxn
      contains

      subroutine kpp_write(modelname)
      character(len=*) :: modelname
      call kpp_def(modelname)
      call kpp_spc(modelname)
      call kpp_rxn(modelname)
      endsubroutine 

      subroutine kpp_def(modelname)
      character(len=*) :: modelname
      character(len=1000) :: f
      integer             :: u

      u=78
      f=trim(modelname)//'.def'
      open(unit=u,file=trim(f))
      write(u,*) '#include '//trim(modelname)//'.spc'
      write(u,*) '#include '//trim(modelname)//'.eqn'
      write(u,*) ''
      write(u,*) '#LANGUAGE Fortran90'
      write(u,*) '#DOUBLE ON'
      write(u,*) '#JACOBIAN SPARSE_LU_ROW'
      write(u,*) '#HESSIAN ON'
      write(u,*) '#STOICMAT ON'
      write(u,*) ''
      write(u,*) '#INTEGRATOR rosenbrock'
      write(u,*) '#DRIVER general'
      write(u,*) ''

      !inline global
      write(u,*) '#INLINE F90_GLOBAL'
      write(u,*) ''
      write(u,*) '#ENDINLINE'
      write(u,*) ''

      !inline rates
      write(u,*) '#INLINE F90_RATES'
      write(u,*) ' function get_rate(i) result(r)'
      write(u,*) ' use module_reaction, only:rate_cst'
      write(u,*) ' real(kind=dp) :: r'
      write(u,*) ' integer  :: i'
      write(u,*) ' r = rate_cst(i)'
      write(u,*) ' endfunction get_rate'
      write(u,*) '#ENDINLINE'
      write(u,*) ''

      !inline util
      write(u,*) '#INLINE F90_UTIL'
      write(u,*) ' subroutine set_fix(F_USER,F)'
      write(u,*) ' real(kind=dp) :: F_USER(NFIX), F(NFIX)'
      write(u,*) ' F(:) = F_USER(:)'
      write(u,*) ' endsubroutine set_fix'
      write(u,*) ' subroutine get_fix(F,F_USER)'
      write(u,*) ' real(kind=dp) :: F_USER(NFIX), F(NFIX)'
      write(u,*) ' F_USER(:) = F(:)'
      write(u,*) ' endsubroutine get_fix' 
      write(u,*) '#ENDINLINE'
      close(u)
      endsubroutine

      subroutine kpp_spc(modelname)
      use module_spec,only:ns,spec,nactive,ninactive

      character(len=*)    :: modelname
      character(len=1000) :: f
      integer             :: u, i 
      f=trim(modelname)//'.spc'
      u=78
      open(unit=u,file=trim(f))
      write(u,*) "#include atoms"
      write(u,*) ""
      write(u,*) "#DEFVAR"
      do i = 1, ns
         if (spec(i)%status.eq.'A') then
            write(u,'(A15,A1,X,A7)') spec(i)%name,"=",'IGNORE;'    
         endif
      enddo 
      write(u,*) "#DEFFIX"

      do i = 1, ns
         if (spec(i)%status.ne.'A') then
            write(u,'(A15,A1,X,A7)') spec(i)%name,"=",'IGNORE;'    
         endif
      enddo
      close(u)
      endsubroutine kpp_spc


      subroutine kpp_rxn(modelname)
      use module_spec,only:ns,spec
      use module_reaction,only:nr,rxn,nphoto
 
      character(len=*)    :: modelname
      character(len=1000) :: f,str
      character(len=5)    :: str_coef
      character(len=4)    :: str_rid
      integer             :: u, i, j
      u=78
      f=trim(modelname)//'.eqn'
      open(unit=u,file=trim(f))
      write(u,*) "#EQUATIONS"
      write(u,*) ""

      do i = 1, nr
         if (rxn(i)%status.ne.'A') cycle
         write(str,'(a2,i3.3,a1)') '<R',i,'>'
         do j = 1, rxn(i)%nreac
            str=trim(str)//' '//trim(spec(rxn(i)%reacs(j))%name)
            if (j.ne.rxn(i)%nreac) str=trim(str)//' +'
         enddo
         str=trim(str)//' = '
         do j = 1, rxn(i)%nprod
            write(str_coef,'(f5.3)')rxn(i)%prod_coefs(j)
            str=trim(str)//' '//str_coef//' '//trim(spec(rxn(i)%prods(j))%name)
            if (j.ne.rxn(i)%nprod) str=trim(str)//' +'
         enddo
         write(str_rid,'(i4)') i
         str=trim(str)//'  : get_rate('// str_rid //');'
         write(u,*) trim(str)
      enddo 
      close(u)
      endsubroutine kpp_rxn
      endmodule module_kpp

