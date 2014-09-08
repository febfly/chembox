C
C File:   smvgear_setpara.f
C Author: Yuzhong
C
C Created on July 25, 2014, 2:32 PM
C
      subroutine smvgear_setpara(ncs, chemintv)
      use mod_comode
      implicit none
      !Input variables
      integer, intent(in) :: ncs
      real(kind=DP), intent(in) :: chemintv
      !Local variables
      !real(kind=DP) :: chemintv
      real(kind=DP) :: errmaxu, ylowu, yhiu, hmaxdayu, hmaxnit
      real(kind=DP) :: abhi, ablo
      integer       :: i
            
      !=========================================================!
      !            Set parameters for smvgear                   !
      !=========================================================!
      !chemintv = 900.
      !ifreord  = 1
      fracdec  = 0.25
      errmaxu  = 1e-3
      ylowu    = 1e3
      yhiu     = 1e7
      hmaxdayu = 9e2
      hmaxnit  = 2e3
      rmserr   = 0.
      iout     = 6
      
      smal1    = 1.0d-06 
      smal2    = 1.0d-100
      smal3    = 1.0d-50
      
      errmax(ncs)       = errmaxu
      timeintv(ncs)     = chemintv
      abst2 (ncs)       = 1./(chemintv*chemintv)
      hmaxuse (ncs)     = hmaxdayu
      hmaxuse (ICS+ncs) = hmaxnit
      
      !Set abtol
      abtol(1, ncs) = yhiu
      abtol(6, ncs) = ylowu         
      abhi = log10(abtol(1,ncs))
      ablo = log10(abtol(6,ncs))         
      do i = 2, 5
          abtol(i, ncs) = 10.**(ablo + (abhi - ablo) * float(6 - i) / 5.)
      enddo
      
      endsubroutine smvgear_setpara
