      program chembox
      use module_init, only:do_initialize
      use module_chem, only:do_chem
      use module_time_common, only:time_step_next
      use module_met_common, only:update_met
      implicit none 
      integer :: time_flag

      !Initialization
      call do_initialize

      time_flag = time_step_next()
      do while (time_flag.ne.2) 
         call update_met
!         call do_chem_constraint !??
         call do_chem
         time_flag = time_step_next()
      enddo


      endprogram chembox
