      module module_geos
      implicit none
 
      contains
       subroutine geos_read(filename)
       use module_spec,only:spec_init,spec_add,spec_finish_add,spec
       use module_reaction
       use module_photolysis,only:photolysis_index
       character(len=*)     :: filename
       integer,parameter  :: u=78
       character*512 :: READIN
       integer       :: nreac, nprod,tpid,iord
       character(len=1) :: stat, tp,tp0 
       integer,parameter :: nread=24
       real        :: pinp(nread)
       character*14::xinp(nread)
c       character*3 ::yinp
       integer     :: i,j,id
       type(reaction) :: rtemp 
       real, dimension(28):: para
!       real, dimension(3) :: paratemp
       character(len=14),dimension(4):: reac_name
       integer, dimension(4)::          reac_id
       character(len=14),dimension(16):: prod_name
       integer, dimension(16)::          prod_id
       real,    dimension(16)::          coef

       
       call spec_init
       call rxn_init

       open(unit=u, file=trim(filename), status='old')

       readin=""
       do while (trim(readin).ne.'BEGIN')
          read(u,*) readin
c          print*,trim(readin)
       enddo
       read(u,*) readin
       read(u,*) readin

       !specie
       read(u,50) stat,xinp(1),(pinp(j),j=1,5)
       read(u,*)  readin
       do while(trim(xinp(1)).ne.'END')
          id=spec_add(xinp(1),stat,pinp(4))   
c          print*,id,xinp(1),stat,pinp(4)
          read(u,50) stat,xinp(1),(pinp(j),j=1,5)
          read(u,*)  readin
       enddo
       call spec_finish_add
      
       !read reaction
       readin=""
       do while (trim(readin).ne.'BEGIN')
          read(u,*) readin
c          print*,trim(readin)
       enddo

       iord=0
       do while (iord.ne.9999)
          call geos_readline(u,reac_name,reac_id,nreac,prod_name,prod_id,
     +         nprod,coef,stat,iord,para,tp,tpid,.false.)
          if (iord.eq.9999) exit
          if (stat.eq.'D') cycle
          rtemp%nreac = nreac
          rtemp%reacs(1:nreac)=reac_id(1:nreac)
          rtemp%nprod = nprod
          rtemp%prods(1:nprod)= prod_id(1:nprod)
          rtemp%prod_coefs(1:nprod) = coef(1:nprod)
          rtemp%status = stat
          rtemp%sn     = iord
          rtemp%r_type = tpid
          rtemp%rindex = 0

      if (tp.eq.'P') then           
          rtemp%paras(1:3) = para(1:3)
          rtemp%paras(4:6) = para(8:10)
!          paratemp(1:3)    = para(8:10)
          rtemp%paras(7)   = para(5)          
      elseif (tp.eq.'E') then
          if (tp0.ne.'P') then
             print*,'last reaction shuld have type P'
             stop
          endif
          rtemp%paras(1:3) = para(1:3)
      elseif (tp.eq.'B') then
          rtemp%paras(1:3) = para(1:3)
          rtemp%paras(4:6) = para(8:10)
      elseif (tp.eq.'A'.and.stat.ne.'D') then
          if (tp0.ne.'B') then
              print*,'last reaction should have type B'
              print*,spec(rtemp%reacs(1:rtemp%nreac))%name
              print*,spec(rtemp%prods(1:rtemp%nprod))%name
              stop
          endif
          rtemp%paras(1:3) = para(1:3)
          rtemp%paras(4) = para(8)
      elseif (tp.eq.'Z') then
          rtemp%paras(1:3) = para(1:3)
          rtemp%paras(4:6) = para(8:10)
      elseif (tp.eq.'V') then
          rtemp%paras(1:3) = para(1:3)
          rtemp%paras(4:6) = para(8:10)
      elseif (tp.eq.'X') then
          rtemp%paras(1:3) = para(1:3)
          rtemp%paras(4:6) = para(8:10)
          rtemp%paras(7:9) = para(15:17)
      elseif (tp.eq.'G') then
          rtemp%paras(1:3) = para(1:3)
          rtemp%paras(4:6) = para(8:10)
      else
          rtemp%paras(1:3) = para(1:3)
      endif
       call rxn_add(rtemp,.false.)
       tp0=tp
       enddo       

       !read photolysis
       readin=""
       do while (trim(readin).ne.'BEGIN')
          read(u,*) readin
       enddo

       iord=0
       do while (iord.ne.9999)
          call geos_readline(u,reac_name,reac_id,nreac,prod_name,prod_id,
     +         nprod,coef,stat,iord,para,tp,tpid,.true.)
          if (iord.eq.9999) exit
          rtemp%nreac = nreac
          rtemp%reacs(1:nreac)=reac_id(1:nreac)
          rtemp%nprod = nprod
          rtemp%prods(1:nprod)= prod_id(1:nprod)
          rtemp%prod_coefs(1:nprod) = coef(1:nprod)
          rtemp%status = stat
          rtemp%sn     = iord
          rtemp%r_type = tpid
          rtemp%paras(1) = para(1)
          rtemp%rindex = photolysis_index(reac_name(1),dble(para(1)))
          call rxn_add(rtemp,.true.)
       enddo
       close(u)
  50   format(A1,1X,A14,3X,0PF6.2,4(1PE10.3))

       endsubroutine geos_read
!=========================================================================
       subroutine geos_readline(u,reac_name,reac_id,nreac,prod_name,prod_id,
     +            nprod,coef,stat,iord,para,tp,tpid,ifphoto)
       use module_spec,only:spec_getid
       use module_reaction
       integer :: u, iord,iordx
       logical :: ifphoto
       character(len=1) :: stat, statx,tp, tpx
       !character(len=2) :: tptmp
       real, dimension(28):: para
       real               :: a, b, e, f, g 
       integer            :: c, d
       character(len=14),dimension(4):: reac_name
       integer, dimension(4)::          reac_id
       character(len=14),dimension(16):: prod_name
       integer, dimension(16)::          prod_id
       real,    dimension(16)::          coef
       integer               ::          tpid
       character(len=20)     :: comment
       integer               :: i, j, nreac,nprod,id
       character(len=1),dimension(16) :: dummychar
       real,            dimension(4) :: dummyreal

       read(u,51) stat,iord,a,b,c,d,tp,e,f,g,comment
       para(1)=a
       para(2)=b
       para(3)=float(c)
       para(4)=float(d)
       para(5)=e
       para(6)=f
       para(7)=g
       tpid = rxn_findtype(tp,2,ifphoto)
       if (tp.eq.'Z'.or.tp.eq.'V'.or.tp.eq.'P'.or.tp.eq.'V'.or.tp.eq.'X'
     & .or.tp.eq.'B'.or.tp.eq.'A'.or.tp.eq.'G') then
           read(u,51) statx,iordx,a,b,c,d,tpx,e,f,g,comment
           para(8)=a
           para(9)=b
           para(10)=float(c)
           para(11)=float(d)
           para(12)=e
           para(13)=f
           para(14)=g
           if(tp.eq.'X') then
              read(u,51) statx,iordx,a,b,c,d,tpx,e,f,g,comment
              para(15)=a
              para(16)=b
              para(17)=float(c)
           endif
       endif

       read(u,52)(dummychar(i),dummyreal(i),reac_name(i),i=1,4)
       if (iord.eq.9999) return

       nreac=0
       do j = 1,4
          if (reac_name(j).ne." ")then
              nreac = nreac + 1
              id = spec_getid(reac_name(j))
              if (id.eq.0) then
                print*,'cannot find reactions in species list:geos-chem'
                print*,reac_name(j)
                 stop
              endif
              reac_id(nreac) = id
          endif
       enddo

       read(u,53)(dummychar(i),coef(i),prod_name(i),i=1,16)
       nprod=0
       do j=  1,16
          if (prod_name(j).ne." ".and.coef(j).ne.0.) then
             nprod = nprod + 1
             id = spec_getid(prod_name(j))
             if (id.eq.0) then
                print*,'cannot find product in species list:geos'
                stop
             endif
             prod_id(nprod) = id
          endif
       enddo

  51   format(A1,1X,I4,1X,ES8.2,1X,ES8.1,1X,I6,1X,I1,1X,A1,1X,F6.2,1X,
     1       2(F6.0,1X),A20)
  52   format(4(A1,0PF5.3,A14))
  53   format(4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14)/
     1        4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14))
       endsubroutine geos_readline
     
      endmodule module_geos
