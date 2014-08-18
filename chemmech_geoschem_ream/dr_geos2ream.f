      program dr_geos2ream
      use module_geos,only:geos_read
      use module_ream,only:ream_write
      character(len=1024) :: pgeos,pream,fgeos,fream

      pgeos='/data13/yzhang425/GEOS-CHEM/standard/'
      fgeos='globchem.dat'

      call geos_read(trim(pgeos)//trim(fgeos))
      call ream_write
      endprogram
