      module module_domain_common
      use module_model_parameter,only: DP, IX,JX,KX,NIJ,NIJK
      implicit none
!      real(kind=DP),dimension(IX,JX) :: lat, lon
      integer,dimension(IX,JX,KX) :: grid_3_1
      integer,dimension(IX,JX)    :: grid_2_1
      integer,dimension(NIJK)     :: grid_1_3_i
      integer,dimension(NIJK)     :: grid_1_3_j
      integer,dimension(NIJK)     :: grid_1_3_k
      integer,dimension(NIJ)      :: grid_1_2_i
      integer,dimension(NIJ)      :: grid_1_2_j

      real(kind=DP)  :: height, width, length

      public :: calc_domain_index
      contains

      subroutine calc_domain_index
      integer :: i, j, k, n
      do k=1,KX
        do j=1,JX
          do i=1,IX
             n=(k-1)*IX*JX + (j-1)*IX +i
             grid_3_1(i,j,k)=n
             grid_1_3_i(n)=i
             grid_1_3_j(n)=j
             grid_1_3_k(n)=k
          enddo
        enddo
      enddo

      do j=1,JX
        do i=1,IX
           n=(j-1)*IX + i
           grid_2_1(i,j) = n
           grid_1_2_i(n) = i
           grid_1_2_j(n) = j
        enddo 
      enddo

      endsubroutine calc_domain_index
      endmodule module_domain_common
