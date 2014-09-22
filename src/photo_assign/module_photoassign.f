      module module_photoassign
      use module_model_parameter,only: DP,KX
      use module_model_parameter,only: MAX_NPHOTO
      implicit none

      integer,parameter                         :: NVERT = KX

      real(kind=DP)                             :: zenfixed
      real(kind=DP),dimension(MAX_NPHOTO,NVERT) :: prate
      
      contains

!============================================================

      subroutine photoassign_read(filename)
      use module_chemmech_common,only:nrxn, nphotorxn,spec_name, prod_id
      character(len=*) :: filename  
      integer          :: j, k, jst, jn, kn, id
      real(kind=DP),dimension(NVERT) :: phold
      open(unit=12,file=filename)
      read(12,370) zenfixed
      do k=1,NVERT
         read(12,380) jst,phold(k),jn,kn
         if (kn.ne.k) then 
            print*,'Error: photorate layer mismatch!'
            stop
         endif
      enddo
      close(12)

      do j=1,nphotorxn
         id = nrxn - nphotorxn + j
         if (jst.eq.spec_name(prod_id(1,id)).and.prate(j).eq.0) then
            do k=1,NVERT
               prate(k,j)=phold(k)
            enddo
         endif
      enddo

 370  format(25x,0pf6.4/)
 380  format(a14,1x,1pe10.4,i5,i7)

      endsubroutine photoassign_read
      endmodule module_photoassign
