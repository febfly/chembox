      module module_time_common
      use module_model_parameter,only: DP,TS_MODEL, TS_CHEM, 
     +    TS_OUTPUT
      implicit none

      integer, parameter :: NPROCESS = 3
      integer, parameter :: IP_BASE = 1, IP_CHEM  = 2, IP_OUTPUT = 3
      real(kind=DP), dimension(NPROCESS) :: timestep
      integer, dimension(NPROCESS)       :: timestep_n
      logical, dimension(NPROCESS)       :: ifdoit

!=======================================================================
! t_total_s     : total simulation time in second
! t_elapse_s    : time at the end of current step in second
! t_remain_s    : time remaining at the end of current step in second
! nstep         : total number of time steps
! istep         : steps finished at the end of current step
! t_st_s        : astronomical Julian day * 86500 at the beginning of 
!                 simulation
! t_ed_s        :astronomical Julian day * 86500 at the end of simulation
! 
!=======================================================================
!      integer, parameter:: FLAG_ST=0, FLAG_END=2, FLAG_NORMAL=1

      real(kind=DP) :: t_total_s
      real(kind=DP) :: t_remain_s, t_elapse_s
      real(kind=DP) :: t_st_s, t_ed_s
      integer       :: year_st, jday_st, hour_st, minute_st, month_st,
     +                 day_st, second_st
      integer       :: year, jday, month, day, hour, minute, second
      integer       :: n_timestep, i_timestep

      public :: time_init
      public :: time_update
      
      contains

      subroutine time_init
      use julday_mod, only: julday
      use ???, only: isint

      real(kind=DP) :: nts_chem_tmp, nts_output_tmp

      !setup time step for varied process
      ts_model_s  = TS_MODEL
      ts_chem_s   = TS_CHEM
      ts_output_s = TS_OUTPUT

      nts_chem_tmp= ts_chem_s/ts_model_s
      if (isint(nts_chem_tmp)) 


      endsubroutine time_init

!======================================================================      
      subroutine time_settimestep
      use ???, only: isint
      integer       :: i
      real(kind=DP) :: tmp

      !base time step
      if (timestep(IP_BASE).le.0) then
         print*,'Error: base time step must be greater than 0'
         stop
      endif

      !calculate the relationship with base time step for each process
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
      enddo
      endsubroutine time_settimestep

      endmodule module_time_common
