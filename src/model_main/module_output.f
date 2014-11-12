      module module_output
      use module_model_parameter,only: DP, NIJK, MAX_STR1
      implicit none
      contains
      subroutine do_output
      use module_chemmech_common,only : spec_getid
      use module_met_common,only : airdensity
      use module_time_common,only : t_elapse_s
      use module_conc_common,only : gas_conc
      use module_domain_common,only : grid_1_3_i,grid_1_3_j,grid_1_3_k
      use module_model_option,only: path_output

      integer,parameter :: NOUT = 9, NPLOUT=1
      character(len=MAX_STR1),dimension(NOUT) :: names
      character(len=MAX_STR1),dimension(NPLOUT) :: plnames
      integer,dimension(NOUT),save :: ind
      integer,dimension(NPLOUT),save :: ind_pl
      logical,save :: first=.true.
      integer      :: is, ig, i, j, k
      character(len=255) ::  filename, outstring
      integer,save :: seqno = 1

      names=(/'O3  ','NO  ','NO2 ','PAN ','CO  ','CH2O','OH  ','HO2 ',
     +        'HNO3'/)
!      plnames=(/'PHOX1 '/)
      if (first) then
         do is=1,NOUT 
            ind(is)=spec_getid(names(is))
            if (ind(is).eq.0) then 
               print*,'Error:Cannot find output species!'
               stop
            endif
         enddo
 
         !==temporary
!         do is=1,NPLOUT
!            ind_pl(is)=spec_getid(plnames(is))
!         enddo
         !==temporary

         write(filename,"(a7,i2.2,a4)")'output.txt'
         open(unit=12,file=trim(path_output)//trim(filename))
         outstring=trim(names(1))

         do is=2,NOUT
            write(outstring,'(a,a1,a)')trim(outstring),',',trim(names(is))
         enddo
         write(12,*) trim(outstring)

         first=.false.
      endif

       seqno = seqno + 1


      !==temporary
!      do is=1,NPLOUT
!         write(outstring,'(a,a1,a)')trim(outstring),',',trim(plnames(is))
!      enddo
      !==temporary

 
      do ig=1,NIJK
         i=grid_1_3_i(ig)
         j=grid_1_3_j(ig)
         k=grid_1_3_k(ig)
         write(outstring,'(f5.1)')
     +        t_elapse_s/3600.
         do is=1,NOUT
            write(outstring,'(a,a1,e10.3)') trim(outstring),',',
     +        gas_conc(ig,ind(is))/airdensity(i,j,k)*1d9
         enddo
!         do is=1,NPLOUT
!            write(outstring,'(a,a1,e10.3)') trim(outstring),',',
!     +        gas_conc(ig,ind_pl(is))/airdensity(i,j,k)*1d9/ts_output_min !ppbv/min
!            gas_conc(ig,ind_pl(is))=0d0
!         enddo
         write(12,*) trim(outstring)
      enddo
!      close(12)
      endsubroutine do_output
      endmodule module_output
