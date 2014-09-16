      module util_time_mod
      implicit none
      integer, parameter :: DP=8
      public :: jd2ymdhms

      contains
      subroutine jd2ymdhms(jd, y,m,d,h,mi,s)
      use julday_mod, only:caldate
      real(kind=DP) :: jd, s, rem
      integer :: y,m,d,h,mi,ymd,hms
      call caldate(jd,ymd,hms)
      y    = int(dble(ymd/1d4))
      rem  = dble(ymd) - (dble(y)*1d4)
      m    = int(rem/1d2)
      rem  = rem-(dble(m)*1d2)
      d    = int(rem)
      h    = int(dble(hms/1d4))
      rem  = dble(hms) - (dble(h)*1d4)
      mi   = int(rem/1d2)
      rem  = rem-(dble(mi)*1d2)
      s    = rem
      endsubroutine jd2ymdhms
      endmodule util_time_mod
