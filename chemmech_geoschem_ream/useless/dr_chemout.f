      program dr_chemout
      use module_reaction
      use module_spec
      use module_ream,only:ream_read
      implicit none
      integer :: specset(64),specset_f(37),rxnset(151),ir,rr
      character(len=1024) :: pream
      pream="/data12/yzhang425/REAM/REAMSource/v2011_daq_12/Run/"
     +    //"chem.dat"
      call ream_read(trim(pream))

      open(unit=10,file='chemout.txt')
      read(10,'(64(i3,x))') specset
      read(10,'(37(i3,x))') specset_f
      read(10,'(151(i3,x))') rxnset
      close(10)

      print*,spec(specset)%name
      print*,spec(specset_f)%name

      do ir=1,151
         rr=rxnset(ir)
         print*,spec(rxn(rr)%reacs(1:2))%name,'=',spec(rxn(rr)%prods(1:rxn(rr)%nprod))%name
      enddo

      endprogram
