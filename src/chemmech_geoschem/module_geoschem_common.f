      module module_geoschem_common
      use module_chemmech_common,only : DP,MAX_NRXN
      implicit none
      !parameters
      integer, parameter :: MAX_NTYPE = 30
      integer, parameter :: MAX_NPARA = 21
      integer, parameter :: STRLEN    = 50
      integer, parameter :: SYMLEN    = 2

      !reaction type for rate calculation 
      integer :: ntype
      character(len=STRLEN),dimension(MAX_NTYPE) :: type_name
      character(len=SYMLEN),dimension(MAX_NTYPE) :: symbol
      character(len=STRLEN),dimension(MAX_NTYPE) :: comment
      integer,dimension(MAX_NTYPE)               :: preceeding_type
      integer,dimension(MAX_NTYPE)               :: succeeding_type
      integer,dimension(MAX_NTYPE)               :: if_photo_type      

      !reaction type and parameters for each reaction
      integer, dimension(MAX_NRXN)                 :: r_type
      real(kind=DP),dimension(MAX_NPARA,MAX_NRXN)  :: paras
      endmodule module_geoschem_common
