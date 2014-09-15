       subroutine sub_chemmech_read
       use module_model_option,only:option_chemmech
       use module_geoschem_io,only:read_geoschem
       use module_ream_io,only:read_ream
       character(len=255) :: filename

       if (option_chemmech.eq.1) then
          filename='globchem.dat'
          call read_geoschem(trim(filename))
       elseif (option_chemmech.eq.2) then
          filename='chem.dat'
          call read_ream(trim(filename))
       endif

       endsubroutine sub_chemmech_read
