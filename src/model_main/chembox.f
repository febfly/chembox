      program chembox
      use module_init, only:do_initialize
      use module_chem, only:do_chem
      use module_output,only:do_output
      use module_time_common, only:ifdoit, model_status, time_update, 
     +                             IP_CHEM, IP_OUTPUT
      use module_met_common, only:update_met
      use module_model_option,only:path_input,path_output
 
      implicit none 
      integer :: time_flag

      call getatt(1,path_input)
      call getatt(2,path_output)

      !Initialization
      call do_initialize

      do while (time_flag.ne.FLAG_END) 
         if (if_time_output) call do_output
         call update_met
         if (ifdoit(IP_CHEM))   call do_chem
         if (ifdoit(IP_OUTPUT)) call do_output
         call time_update()
      enddo


      endprogram chembox
