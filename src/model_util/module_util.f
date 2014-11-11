      module module_util
      integer,parameter: DP=8
      real(kind=DP),parameter :: EPS = 1d-8
      public:: isint
      contains
      
!============================================================
! This function checks if a double number is a
      function isint(x) result(r)
      real(kind=DP),intent(in) :: x
      logical                  :: r
      integer                  :: xint

      r=.false.
      xint=nint(x)
      if (x.le.xint+EPS.and.x.ge.xint-EPS) r=.true.
          
      endfunction isint

      endmodule module_util
