      program dr_geos2kpp
      use module_geos,only:geos_read
      use module_kpp,only:kpp_write
      character(len=1024) :: pgeos,pream,fgeos,fream

      pgeos='/data13/yzhang425/GEOS-CHEM/runtest/'
      fgeos='globchem.dat'

      call geos_read(trim(pgeos)//trim(fgeos))
      call kpp_write('geos')
      endprogram
