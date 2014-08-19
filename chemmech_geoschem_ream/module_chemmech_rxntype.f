      module module_ream_rxntype
      use module_model_parameter,only: DP
      implicit none
      integer, parameter :: MAX_NTYPE = 30
      integer, parameter :: MAX_NPARA = 15
      integer, parameter :: STRLEN    = 50

      integer :: ntype
      
      character(len=STRLEN),dimension(MAX_NTYPE) :: type_name
      character(len=1),dimension(MAX_NTYPE)      :: symbol
      character(len=STRLEN),,dimension(MAX_NTYPE):: comment
      integer,dimension(MAX_NTYPE)               :: nparameter
      integer,dimension(MAX_NPARA,MAX_NTYPE)     :: parameter_pos
      integer,dimension(MAX_NTYPE)               :: preceeding_type
      integer,dimension(MAX_NTYPE)               :: succeeding_type
      integer,dimension(MAX_NTYPE)               :: if_photo_type

      private :: rxntype_define
      contains

      subroutine rxntype_define
      parameter_pos(:,:)     = 0
      ntype                  = 16

      type_name(1)           = 'Normal'
      nparameter(1)          = 3
      parameter_pos(1:3,1)   = (/1,2,3/)
      symbol(1)              = " "
      comment(1)             = "A*(300/T)**B*EXP(C/T)"
      preceeding_type(1)     = 0
      succeeding_type(1)     = 0
      if_photo_type(1)       = 0

      type_name(2)           = 'Pressure dependent 3 body'
      nparameter(2)          = 7
      parameter_pos(1:3,2)   = (/1,2,3/)
      symbol(2)              = "P"
      comment(2)             = "Kh*G**(1/(1+(log(Kh/Kl))^2))"
      preceeding_type(2)     = 0
      succeeding_type(2)     = 0      
      if_photo_type(2)       = 0

      type_name(3)           = 'Addition branch of RO2+NO'
      nparameter(3)          = 4
      parameter_pos(1:3,3)   = (/1,2,3/)
      symbol(3)              = "B"
      comment(3)             = "k*(1-fyrno3)"
      preceeding_type(3)     = 0
      succeeding_type(3)     = 4
      if_photo_type(3)       = 0

      type_name(4)           = 'Abstaction branch of RO2+NO'
      nparameter(4)          = 4
      parameter_pos(1:3,4)   = (/1,2,3/)
      symbol(4)              = "A"
      comment(4)             = "k*fyrno3"
      preceeding_type(4)     = 3
      succeeding_type(4)     = 0
      if_photo_type(4)       = 0

      type_name(5)           = 'Equilibrim reaction (reverse)'
      nparameter(5)          = 3
      parameter_pos(1:3,5)   = (/1,2,3/)
      symbol(5)              = "E"
      comment(5)             = ""
      preceeding_type(5)     = 2
      succeeding_type(5)     = 0
      if_photo_type(5)       = 0

      type_name(6)           = 'OH+HNO3'
      nparameter(6)          = 9
      parameter_pos(1:3,6)   = (/1,2,3/)
      symbol(6)              = "X"
      comment(6)             = "K1 + K3[M] / (1 + K3[M]/K2)"
      preceeding_type(6)     = 0
      succeeding_type(6)     = 0
      if_photo_type(6)       = 0

      type_name(7)           = 'OH+CO'
      nparameter(7)          = 3
      parameter_pos(1:3,7)   = (/1,2,3/)
      symbol(7)              = "Y"
      comment(7)             = ""
      preceeding_type(7)     = 0
      succeeding_type(7)     = 0
      if_photo_type(7)       = 0

      type_name(8)           = 'HO2/NO3+HO2'
      nparameter(8)          = 6
      parameter_pos(1:3,8)   = (/1,2,3/)
      symbol(8)              = "Z"
      comment(8)             = "(K1 + K2)*(1+1.4E-21*[H2O]*EXP(2200/T))"
      preceeding_type(8)     = 0
      succeeding_type(8)     = 0
      if_photo_type(8)       = 0

      type_name(9)           = 'GLYX+OH/NO3->HO2+CO'
      nparameter(9)          = 3
      parameter_pos(1:3,9)   = (/1,2,3/)
      symbol(9)              = "C"
      comment(9)             = "k1*([o2]+3.5d18)/(2*[o2]+3.5d18)"
      preceeding_type(9)     = 0
      succeeding_type(9)     = 0
      if_photo_type(9)       = 0

      type_name(10)          = 'GLYX+OH/NO3->GLCO3'
      nparameter(10)         = 3
      parameter_pos(1:3,10)  = (/1,2,3/)
      symbol(10)             = "D"
      comment(10)            = "k1*[o2]/(2*[o2]+3.5d18)"
      preceeding_type(10)    = 0
      succeeding_type(10)    = 0
      if_photo_type(10)      = 0

      type_name(11)          = 'Aerosol surface loss'
      nparameter(11)         = 3
      parameter_pos(1:3,11)  = (/1,2,3/)
      symbol(11)             = "K"
      comment(11)            = ""
      preceeding_type(11)    = 0
      succeeding_type(11)    = 0
      if_photo_type(11)      = 0

      type_name(12)          = 'EO2->HO2+GLYC'
      nparameter(12)         = 3
      parameter_pos(1:3,12)  = (/1,2,3/)
      symbol(12)             = "H"
      comment(12)            = "k1*FYHORO"
      preceeding_type(12)    = 0
      succeeding_type(12)    = 0
      if_photo_type(12)      = 0

      type_name(13)          = 'EO2->HO2+2CH2O'
      nparameter(13)         = 3
      parameter_pos(1:3,13)  = (/1,2,3/)
      symbol(13)             = "F"
      comment(13)            = "k1*(1-FYHORO)"
      preceeding_type(13)    = 0
      succeeding_type(13)    = 0
      if_photo_type(13)      = 0

      type_name(14)          = 'Temperature dependent branching'
      nparameter(14)         = 3
      parameter_pos(1:3,14)  = (/1,2,3/)
      symbol(14)             = "V"
      comment(14)            = "k1/(1+k2)"
      preceeding_type(14)    = 0
      succeeding_type(14)    = 0
      if_photo_type(14)      = 0

      type_name(15)          = 'Photolysis'
      nparameter(15)         = 3
      parameter_pos(1:3,15)  = (/1,2,3/)
      symbol(15)             = " "
      comment(15)            = ""
      preceeding_type(15)    = 0
      succeeding_type(15)    = 0
      if_photo_type(15)      = 1

      type_name(16)          = 'O3 photolysis'
      nparameter(16)         = 3
      parameter_pos(1:3,16)  = (/1,2,3/)
      symbol(16)             = "A"
      comment(16)            = ""
      preceeding_type(16)    = 0
      succeeding_type(16)    = 0
      if_photo_type(16)      = 1
      endsubroutine rxntype_define

!=========================================================================
      function rxntype_id(sym,ifphoto) result(id)
       character(len=1) :: sym
       integer          :: id
       integer          :: i,st,ed
       logical          :: ifphoto
       st=1
       ed=ntype
       do i = st, ed
          if (sumbol(i).eq.sym.and.if_photo_type(i).eq.ifphoto) exit
       enddo
       id = i
       if (id.eq.ed+1)   id = 1
       !if (ifphoto.eq.1) id =15!overwrite type 16...
      endfunction rxntype_id
!=========================================================================

      function rxn_rate(tp, np, p, Temp, Pres, O2, N2, M, H2O,
     +                  aer_area, aer_radius, paraI1, paraF1)
     +result(rate_const) 
      integer, intent(in)                    :: tp, np
      real(kind=DP),dimension(np),intent(in) :: p
      real(kind=DP) :: Temp, Pres, O2, N2, M, H2O
      real(kind=DP) :: aer_area, aer_radius, paraF1
      real(kind=DP) :: 
      endfunction rxn_rate
     
      endmodule module_ream_rxntype
