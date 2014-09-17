      module module_met_common
      use module_model_parameter,only:DP, IX, JX, KX
      implicit none
      real(kind=DP),dimension(IX,JX,KX) :: temperature!K
      real(kind=DP),dimension(IX,JX,KX) :: pressure   !Pa
      real(kind=DP),dimension(IX,JX,KX) :: water_vapor!molec/cm3
      real(kind=DP),dimension(IX,JX,KX) :: aer_area   !cm2/cm3
      real(kind=DP),dimension(IX,JX,KX) :: aer_radius !cm
      real(kind=DP),dimension(IX,JX,KX) :: airdensity !molec/cm3
      integer, dimension(IX,JX)         :: ifsun
    
      public :: update_met

      contains
      subroutine update_met
      real(kind=DP),parameter :: Na=6.02d23, R=8.314
      temperature(:,:,:) = 300.
      pressure(:,:,:)    = 1d5
      airdensity = pressure*Na/R/temperature/1d6
      aer_area = 12e-7
      aer_radius = 1e-4
      water_vapor = 0.001*airdensity
      ifsun = 1

      endsubroutine update_met
      endmodule module_met_common
