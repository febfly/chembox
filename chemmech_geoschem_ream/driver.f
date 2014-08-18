      program driver
      use module_spec
      use module_reaction
      use module_geos,only:geos_read
      use module_ream,only:ream_read
      use module_photolysis,only:photolysis_read,photolysis_jrate
      implicit none
      real*8              :: Temp, Pres, H2O, M, O2,aer_area,aer_radius
      integer             :: iday
      real*8,dimension(7) :: fff
 
      character(len=1024) :: pgeos,pream,fgeos,fream,fratj,fjspe
      character(len=20)   :: fmt1,fmt2
      integer             :: i,n,l
      Temp=300
      Pres=1000
      H2O =1e13
      M   =2.5e19
      O2  =0.2*M
      aer_area=5e-11!cm2/cm2
      aer_radius=2e-4!cm
    
      fff(1)=1.6e12
      fff(2)=1.2e14
      fff(3)=2.5e14
      fff(4)=1e15
      fff(5)=4.9e15
      fff(6)=1.97e16
      fff(7)=2.1e17


      pgeos='/data13/yzhang425/GEOS-CHEM/standard/'
      fgeos='globchem.dat'
      fratj='ratj.d'
      fjspe='jv_spec.dat'

      pream='/data12/yzhang425/REAM/REAMSource/v2010_julyo3/Run/'
      fream='chem.dat'

      call photolysis_read(trim(pgeos)//trim(fratj),
     +                     trim(pgeos)//trim(fjspe))

      call geos_read(trim(pgeos)//trim(fgeos))
      do i=1,ns
         print*,i,spec(i)%name
      enddo

      call photolysis_jrate(fff,Temp,Pres)
      call rxn_update_rates(Temp,Pres,M,O2,H2O,aer_area,aer_radius)

      do i=1,nr
         n=rxn(i)%nreac
         print*,i,spec(rxn(i)%reacs(1:n))%name,rate_cst(i)
      enddo
      print*,nphoto
      endprogram
