      module module_time_common
      use module_model_parameter,only:DP,TSCHEM=>TS_CHEM_MIN,
     +    TSM3D=>TS_MODEL_MIN
      implicit none
!=======================================================================
! t_total_min   : total simulation time in minute
! t_elapse_min  : time at the end of current step in minute
! t_remain_min  : time remaining at the end of current step in minute
! t_laststep_min: time at the end of last step and at the beginning of 
!                 the current one, in minute
! nstep         : total number of time steps
! istep         : steps finished at the end of current step
! t_st_min      : astronomical Julian day * 1440 at the beginning of 
!                 simulation
! t_ed_min      :astronomical Julian day * 1440 at the end of simulation
! 
!=======================================================================

      integer,parameter:: FLAG_END=2, FLAG_LAST=1, FLAG_NORMAL=0

      real(kind=DP) :: t_total_min
      real(kind=DP) :: t_remain_min, t_elapse_min
      integer       :: ntstep, itstep
      real(kind=DP) :: t_st_min, t_ed_min
      integer       :: year, jday, hour, minute, month, day
      real(kind=DP) :: second
      !integer       :: year_ls, jday_ls, 

      real(kind=DP) :: ts_model_min
      real(kind=DP) :: ts_chem_min

      public :: time_step_init
      public :: time_step_next

      contains

!=======================================================================
! time_step_init
!     1. Initialize time step
!     2. Calculate simulation time in minute
!     3. Initialize other time variables such as year, month, day, 
!        hour, etc.
!=======================================================================
      subroutine time_step_init(y0,m0,d0,h0,mi0,s0,
     +                          y1,m1,d1,h1,mi1,s1)
      use julday_mod, only : julday
      integer, intent(in) :: y0, m0, d0, h0, mi0
      integer, intent(in) :: y1, m1, d1, h1, mi1
      real(kind=DP),intent(in) :: s0, s1
      real(kind=DP) :: frac_day0, frac_day1

      ts_model_min = TSM3D
      ts_chem_min  = TSCHEM
      frac_day0 = d0+float(h0)/24d0+float(mi0)/1440d0+s0/86400d0
      frac_day1 = d1+float(h1)/24d0+float(mi1)/1440d0+s1/86400d0
      t_st_min = julday(y0,m0,frac_day0)*1440d0
      t_ed_min = julday(y1,m1,frac_day1)*1440d0
      if (t_st_min.ge.t_ed_min) then
         print*,'Error: starting time is later than ending time'
         stop
      else 
         t_total_min    = t_ed_min - t_st_min
         t_remain_min   = t_total_min
         t_elapse_min   = 0d0
         ntstep  = ceiling(t_total_min/ts_model_min)
         itstep = 0
         year   = y0
         month  = m0
         day    = d0
         minute = mi0
         second = s0
         jday   = int(t_st_min/1440d0 - julday(y0,1,1d0) + 1)
      endif
      
      endsubroutine time_step_init

!======================================================================
! time_step_next
!     move to the next time step of ts_model_min
!     flag = 0, normal
!          = 1, last step
!          = 2, simulation finished   
!======================================================================
      function time_step_next() result(flag)
      use julday_mod,only:julday
      use util_time_mod,only:jd2ymdhms
      integer :: ymd, hms, flag
      real(kind=DP) :: jd, rem
      itstep = itstep + 1
      if (itstep.lt.ntstep) then 
         t_remain_min  = t_remain_min - ts_model_min
         t_elapse_min  = t_elapse_min + ts_model_min
         jd = (t_st_min+t_elapse_min)/1440d0
         call jd2ymdhms(jd,year,month,day,hour,minute,second)
         jday = int(t_elapse_min/1440d0 - julday(year,1,1d0) + 1)
         flag = FLAG_NORMAL
      elseif (itstep.eq.ntstep) then
         ts_model_min = t_total_min - t_elapse_min
         ts_chem_min  = ts_model_min
         t_remain_min = 0
         t_elapse_min = t_total_min
         jd = (t_st_min+t_elapse_min)/1440d0
         call jd2ymdhms(jd,year,month,day,hour,minute,second)
         jday = int(t_elapse_min/1440d0 - julday(year,1,1d0) + 1)
         flag = FLAG_LAST
      elseif (itstep.gt.ntstep) then 
         print*,'Simulation finished'
         flag = FLAG_END
      endif
      if (flag.ne.FLAG_END) print "(A,i4,'-',i2,'-',i2,x,i2,':',i2)",
     +      'Current time:',year,month,day,hour,minute
      endfunction time_step_next

      endmodule module_time_common
