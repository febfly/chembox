      program chembox
      use module_init, only:do_initialize
      use module_time_common, only:time_step_next
      implicit none 
      integer :: time_flag
!Initialization
      call do_initialize

      time_flag = time_step_next()
      do while (time_flag.ne.2) 
         call do_met !??
         call do_chem_constraint !??
         call do_chem
         time_flag = time_step_next()
      enddo
      endprogram chembox
