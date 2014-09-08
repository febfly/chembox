pro joinfunction
    func=['smvgear.f','jsparse.f','ksparse.f','backsub.f','decomp.f','pderiv.f','subfun.f','update.f']
    openw,11,'mod_smvgear_core.f'
    str=''
    for i=0,n_elements(func)-1 do begin
        openr,12,func(i)
        while ~eof(12) do begin
          readf,12,str
          printf,11,str
        endwhile
        close,12
        printf,11,''
        printf,11,'!============================================================'
        printf,11,''
    endfor
    close,11

end
