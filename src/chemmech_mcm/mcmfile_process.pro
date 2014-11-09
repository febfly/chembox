pro mcmfile_process,file

  openr,11,file
  openw,12,'mcm_chem.txt'
  openw,13,'mcm_rate.f'
  openw,14,'mcm_photolysis.txt'
  text=''
  while text ne ' Rate expressions fortran routine' do begin
    printf,12,text
    readf,11,text
  endwhile

  readf,11,text
  while strmid(text,1,14) ne 'photolysis.txt' do begin
    printf,13,text
    readf,11,text
  endwhile

  printf,14,text
  readf,11,text
  while strmid(text,1,56) ne 'Fortran routine for calculation of rhs of rate equations' do begin
    printf,14,text
    readf,11,text
  endwhile

  while strmid(text,1,21) ne 'Subset reactants file' do begin
    readf,11,text
  endwhile

  printf,12,text
  readf,11,text
  while ~eof(11) do begin
    printf,12,text
    readf,11,text
  endwhile
  printf,12,text
    
  close,11
  close,12
  close,13
  close,14

end
