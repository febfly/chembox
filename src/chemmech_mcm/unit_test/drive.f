      program drive
      use module_mcm_cheminfo,only:mcm_read
      use module_mcm_rate,only:constants
      implicit none
      real(kind=8),dimension(20000) :: rconst
      real(kind=8),dimension(3487)  :: conc
      real(kind=8)                  :: time
      conc(:)=1d-12
      time=0
      call mcm_read
      call constants(rconst,time,conc)
      endprogram
