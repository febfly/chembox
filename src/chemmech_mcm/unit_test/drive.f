      program drive
      use module_mcm_cheminfo,only:mcm_read
      use module_mcm_rate,only:constants,zenith
      implicit none
      real(kind=8),dimension(20000) :: rconst
      real(kind=8),dimension(3487)  :: conc
      real(kind=8)                  :: time, z,c,s
      integer                       :: i
      conc(:)=1d-12
      time=0
      do i=0,48
         call zenith(z,c,s,i*3600d0)
         print*,i,z
      enddo

      call mcm_read
      call constants(rconst,time,conc)
      endprogram
