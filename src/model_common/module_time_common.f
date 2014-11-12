      module module_time_common
      use module_model_parameter,only: DP,TS_MODEL, TS_CHEM, 
     +    TS_OUTPUT
      implicit none

!=======================================================================
! timestep(NPROCESS)   : timestep of each process in second
! timestep_n(NPROCESS) : ratio between a timestep and the base timestep
!                        must be an integer >= 1
! ifdoit(NPROCESS)     : if the model needs to compute a particular
!                        process in the following time step
!-----------------------------------------------------------------------
!
! modle_status  : the status of the model iteration
!   = FLAG_START  : the first loop
!   = FLAG_NORMAL : a normal loop
!   = FLAG_END    : itstep has reached ntstep, the simulation has
!                   finished
!
!-----------------------------------------------------------------------
! t_total_s     : total simulation time in second
! t_elapse_s    : time at the end of current step in second
! t_remain_s    : time remaining at the end of current step in second
! ntstep        : total number of time steps
! itstep        : count of time steps, from 0-ntstep-1
!                 when itstep >= ntstep, model_status become FLAG_END
! t_st_s        : astronomical Julian day * 86500 at the beginning of 
!                 simulation
! t_ed_s        :astronomical Julian day * 86500 at the end of simulation
! 
!=======================================================================

      !time step control
      integer, parameter :: NPROCESS = 3
      integer, parameter :: IP_BASE = 1, IP_CHEM  = 2, IP_OUTPUT = 3
      real(kind=DP), dimension(NPROCESS) :: timestep
      integer, dimension(NPROCESS)       :: timestep_n
      logical, dimension(NPROCESS)       :: ifdoit
      integer                            :: ntstep, itstep

      !model status control
      integer           :: model_status
      integer, parameter:: FLAG_START=0, FLAG_END=2, FLAG_NORMAL=1

      !time variables
      real(kind=DP) :: t_total_s
      real(kind=DP) :: t_remain_s, t_elapse_s
      real(kind=DP) :: t_st_s, t_ed_s
      integer       :: year_st, hour_st, minute_st, month_st,
     +                 day_st
      integer       :: year_ed, hour_ed, minute_ed, month_ed,
     +                 day_ed
      integer       :: year, month, day, hour, minute
      real(kind=DP) :: jday, second_st, second_ed, second


      public :: time_init
      public :: time_update

      private:: time_settimestep
      private:: time_settime
      
      contains

!======================================================================
! time_init
!     Initialize the module
!=====================================================================
      subroutine time_init
      real(kind=DP) :: nts_chem_tmp, nts_output_tmp

      model_status = FLAG_START
      !setup time step for varied process
      timestep(IP_BASE)   = TS_MODEL
      timestep(IP_CHEM)   = TS_CHEM
      timestep(IP_OUTPUT) = TS_OUTPUT
      call time_settimestep

      !initialize time variables
      call time_settime

      endsubroutine time_init

!======================================================================      
! time_update
!     Update the time variables at the end of an iteration and control 
!     what process to compute in the next time step
!======================================================================
      subroutine time_update
      use julday_mod,   only: julday
      use util_time_mod,only: jd2ymdhms
   
      integer :: ymd, hms
      real(kind=DP) :: jd, rem
      integer :: i
 
      itstep = itstep + 1

      !update time variables for next time step
      if (itstep.lt.ntstep) then
         t_remain_s = t_remain_s - timestep(IP_BASE)
         t_elapse_s = t_elapse_s + timestep(IP_BASE)
         jd = (t_st_s + t_elapse_s)/86400d0
         call jd2ymdhms(jd, year, month, day, hour, minute, second)
         jday = jd - julday(year,1, 1d0)
         model_status = FLAG_NORMAL
      elseif (itstep.ge.ntstep) then
         print*, 'Simulation finished'
         model_status = FLAG_END 
      endif
 
      !determine whether process i will be executed in next time step
      do i=1,NPROCESS
         if (mod(itstep,timestep_n(i)).eq.0) then
            ifdoit(i)=.true.
         else 
            ifdoit(i)=.false.
         endif
      enddo

      endsubroutine time_update

!======================================================================      
! time_settimestep
!     Set up the timestep control when the module is initialized
!======================================================================
      subroutine time_settimestep
      use module_util, only: isint
      integer       :: i
      real(kind=DP) :: tmp

      !base time step
      if (timestep(IP_BASE).le.0) then
         print*,'Error: base time step must be greater than 0'
         stop
      endif

      !calculate the relationship with base time step for each process
      !timestep_n is t
      do i=1,NPROCESS
         if (timestep(i).le.0) timestep(i)=timestep(IP_BASE)
         tmp = timestep(i)/timestep(IP_BASE)
         
         if (.not.isint(tmp)) then
            print*,'Error: the time step for each process must be a'
     +        //'multiple of the base time step', 
     +        timestep(i), timestep(IP_BASE)
            stop
         endif

         timestep_n(i)=nint(tmp)

         if (timestep_n(i).lt.1) then
            print*,'Error: the time step for each process must be '
     +        // 'equal to or greater than the base time step',
     +        timestep(i), timestep(IP_BASE)
            stop
         endif
      enddo
      endsubroutine time_settimestep

!======================================================================      
! time_settime
!     Set the time variables when the module is initialized
!======================================================================
      subroutine time_settime
      use julday_mod,   only: julday
      real(kind=DP) :: frac_day_st, frac_day_ed
      
      frac_day_st = day_st+float(hour_st)/24d0+float(minute_st)/1440d0
     +            + second_st/86400d0
      frac_day_ed = day_ed+float(hour_ed)/24d0+float(minute_ed)/1440d0
     +            + second_ed/86400d0
      t_st_s      = julday(year_st, month_st, frac_day_st)*86400d0
      t_ed_s      = julday(year_ed, month_ed, frac_day_ed)*86400d0

      if (t_st_s.ge.t_ed_s) then
         print*,'Error: starting time is later than ending time'
         stop
      endif

      t_total_s  = t_ed_s - t_st_s
      t_remain_s = t_total_s
      t_elapse_s = 0d0

      ntstep     = ceiling(t_total_s / timestep(IP_BASE))
      itstep     = 0
      year       = year_st
      month      = month_st
      day        = day_st
      minute     = minute_st
      second     = second_st
      jday       = julday(year_st, month_st, frac_day_st) 
     +           - julday(year_st, 1, 1d0) 

      ifdoit     = .true.
      endsubroutine time_settime      
!======================================================================      

      endmodule module_time_common
