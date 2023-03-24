! *** ******************************************************************
      subroutine readAndStoreRslFFR(MPprocNum,fileName,hpcprefix,
     &                              procIDVals,vrtxIDVals,
     &                              boundaryYp,bndYpWallfn,
     &                              animeStat,rsOpenStat,rsCloseStat,
     &                              wallZero,RelativeVelocity,rtrnStat)
!      subroutine readAndStoreRslFFR(MPprocNum,fileName,hpcprefix,
!     &                              procIDVals,boundaryYp,bndYpWallfn,
!     &                              animeStat,rsOpenStat,rsCloseStat,
!     &                              wallZero,RelativeVelocity,rtrnStat)
! *** ******************************************************************
      use FVdata, only : timeFV
      use FFRdata,only : nvrtx,ncell,numVars,nameVars,
     &                   numBVars,nameBVars,bnddata,
     &                   iterFFR,ncvFFR,ncvfacFFR,ncompFFR,ncompallFFR,
     &                   nrnsxFFR,ieul2ph,ical_sufFFR,timeFFR,
     &                   icon_cvx,icon_cvg,
     &                   ical_sldFFR,NMATFFR,ishaft,rotati,rot_ang,
     &                   end,begin,mat_no,nofld,
     &                   lvcell,lacell,cord,totdata,node_id,
     &                   cord_bak,
     &                   NBOUND,SFBOUN,NFBOUN,IBFACE,IFFACE,NBFS,
     &                   LBC_INDEX,LBC_ICV,
     &                   SDOT_suf,WDOT_gas,depoeth_spd,molefr_suf,
     &                   nsdot,ndeps,num_ratout,nmolfr,blk_thick,nthick,
     &                   ical_MHD,npotnxFFR,NFLIDxFFR,radflag,
     &                   idrdpFFR,compFFR,ical_WDRFFR,N_INJ
      use FFRreaddata, only : iuvw_ave_rms_rex,iuvwt_rex,
     &                        ip_avex,it_avex,ip_rmsx,it_rmsx,imu_avex,
     &                        icomp_avex,irans_avex,imu_rmsx,
     &                        icomp_rmsx,irans_rmsx,iuvw_rans_rex,
      ! Tagami
     &                        iuvwc_rex,
      ! imagaT
     &                        iwallvar,imaxminx
      use FFRreaddata, only : ianim,ianim_uvw,ianim_p,ianim_r,ianim_t,
     &                        ianim_GAS_WDOT,ianim_SDOT_suf,
     &                        ianim_molefr_suf
!
      implicit none
!     
      integer,     intent(in)  :: MPprocNum
      character*80,intent(in)  :: fileName
      character*80,intent(in)  :: hpcprefix
      logical,     intent(in)  :: procIDVals
      logical,     intent(in)  :: vrtxIDVals
!
      logical,     intent(in)  :: boundaryYp
      character*80,intent(in)  :: bndYpWallfn
      logical,     intent(in)  :: rsOpenStat
      logical,     intent(in)  :: rsCloseStat
      logical,     intent(out) :: animeStat
! onishi
      logical,     intent(in)  :: wallZero
      logical,     intent(in)  :: RelativeVelocity
      integer,     intent(out) :: rtrnStat
!
      integer :: CPU
      integer,save :: ifl,ifla_allo=0
      integer :: nnode,ncvTotal,icv,jcv
      integer,allocatable :: ivp(:)
      integer :: i,j,k,l,num,knum,tempI1, ierr
      integer :: IIMAT,IMAT,iv,ic,N_MHD,IMAT_U
      real(8) :: up,yp,utau,rho,rmu,rnu,xxx,yyy
      real(8),parameter :: SML=1.0d-25
      real(8),allocatable :: tempReal1(:),tempReal2(:,:),tempReal3(:)
      real(8),allocatable :: tempReal4(:)
! sliding
      real(8) :: alpha
      integer,allocatable :: iiflag(:)
      real*8  :: unit(3),tha,bb(3,3),rbb(3,3),radi(3),vr(3),dr,v0(3)
! onishi
      logical :: tmpflag
      integer :: nb,nbcnd,bnum,IBFL,IBFS,IBFE,ICFL,iii,IMHD
      integer :: nnCV=0
      character*80,allocatable :: boundName(:)
!
      character*80  :: newVarName
      character*255 :: filePath
      character*255 :: hpcPath    ! extern function
      integer :: icord=0
      ! Tagami
      character(80) :: char_var
      ! imagaT
      
      integer :: vertex_cen,cell_cen,ibot,u_func(1:200),uf_mm,
     &     ista_momeryx,ical_prtx
!
! --- -------------------------------------------------------
!     read 'reslt.frontflow' for every CPU, respectively
! --- -------------------------------------------------------
      ncvTotal=0
      
      do 800 CPU=1,MPprocNum
! ---    if rsOpenStat==.true., open result file
!        (at first on animation steps).
         ifl=30+CPU
         if(rsOpenStat) then
!           create file name/fp
            if(MPprocNum>1) then
!              Multi CPU
               filePath=hpcPath(CPU-1,fileName,hpcprefix)
            else
!              Serial CPU
               filePath=fileName
            end if
            write(*,*) ' #### Open result files...'
            write(*,*) ' #### file=',trim(filePath)
            OPEN(ifl,FILE=filePath,FORM='UNFORMATTED',
     &           action='read',iostat=ierr)
            if(ierr/=0) then
               write(*,*)'ABEND:Cannot open ',filePath
               stop
            endif
         end if

! --- -------------------------------------------------------
! ---    read header of FFR ---------------------------------
! --- -------------------------------------------------------
!
         num_ratout=0
         WDOT_gas=0
         if(.not.animeStat) then ! anime==.false.
            read(ifl)iterFFR,timeFFR,NMATFFR,ncvFFR,ncvfacFFR,ncompFFR,
     &           nrnsxFFR,npotnxFFR,numVars,
     &           NFLIDxFFR,idrdpFFR,compFFR,
     &           numBVars,nbcnd,ncompallFFR,
     &           num_ratout,WDOT_gas,radflag,
     &           icon_cvx,vertex_cen,cell_cen,ibot,uf_mm,u_func(1:uf_mm)
     &           ,ical_prtx
            if(uf_mm/=200) stop 'uf_mm/=200'
            if(icon_cvx/=icon_cvg) then
              if(icon_cvx==cell_cen) then
                write(*,*) 'ERR-MSG: result file is [cell-center]'
                write(*,*) 'ERR-MSG: grid file is [vertex-center]'
                write(*,*) 
     &           'ERR-MSG: re-running [prefflow] by "CV=cell"'
                write(*,*) 'ERR-MSG: then running ffr2viz'
              else
                write(*,*) 'ERR: result file is [vertex-center]'
                write(*,*) 'ERR-MSG: grid file is [cell-center]'
                write(*,*) 
     &       'ERR-MSG: re-running [prefflow] by "CV=vertex"'
                write(*,*) 'ERR-MSG: then running ffr2viz'
              endif
              stop 'ERR at: readAndStoreRsiFFR'
            endif
            if(icon_cvx==icon_cvg) then
              if(ifla_allo==0) then
                allocate(tempReal3(0:nvrtx),stat=ierr)
                if(ierr.ne.0) then
                  write(*,*) 'stop-1 at allocating tempReal3',
     &              'in readAndStoreRslFFR.'
                  stop
                endif
                allocate(tempReal4(0:nvrtx),stat=ierr)
                if(ierr.ne.0) then
                  write(*,*) 'stop-1 at allocating tempReal4',
     &                'in readAndStoreRslFFR.'
                 stop
                endif
              endif
              ifla_allo=1
            endif
            read(ifl) 
     &           ical_sldFFR,ical_sufFFR,ieul2ph,ical_MHD,ical_WDRFFR,
     &           N_INJ
            if(ierr/=0) then
               rtrnStat=1
               return
            end if

         else                    ! anime==.true.
!----------------
! --- animation
!----------------
            read(ifl, iostat=ierr) iterFFR,timeFFR,NMATFFR,ncvFFR,
     &           ncvfacFFR,ncompFFR,nrnsxFFR,
     &           numVars ,ncompallFFR     !,ieul2ph
            if(ierr/=0) then
               rtrnStat=-1
               animeStat=.false.
               close(ifl)
               return
            end if
!
            read(ifl)ical_sldFFR,ical_sufFFR
            if(ierr/=0) then
               rtrnStat=1
               return
            end if

         end if
! ---
! ---
! ---    sliding mesh  ---------------------------------
! ---
! DEBUG
         if(ical_sldFFR==1.or.ical_sldFFR==2.or.ical_sldFFR==3
     &    .or.ical_sldFFR==4) then
!            if(CPU==1)then
              allocate(ishaft(1:NMATFFR),rotati(1:NMATFFR))
              allocate(end(1:3,1:NMATFFR),begin(1:3,1:NMATFFR))
              allocate(mat_no(1:NMATFFR))
              allocate(nofld(-NMATFFR:NMATFFR))
              allocate(rot_ang(1:NMATFFR))
!            endif

            nofld(:)=0
            do I=1,NMATFFR
               read(ifl) IIMAT,IMAT,IMAT_U,ishaft(I),rotati(I),
     &              end(1:3,I),begin(1:3,I),rot_ang(I)
               mat_no(I)=IMAT
           write(*,*)'>>>',IMAT,IMAT_U,NMATFFR
!               nofld(IMAT)=IMAT_U
               nofld(I)=IMAT_U
            enddo
         endif


!
         SDOT_suf=0
         molefr_suf=0
         depoeth_spd=0
         nsdot=0
         ndeps=0
         nmolfr=0
         blk_thick=0
         nthick=0
         if(animeStat) then
           if(ical_sufFFR==1) then
             allocate(ianim_GAS_WDOT(ncompFFR),
     &                ianim_SDOT_suf(1:ncompallFFR),
     &                ianim_molefr_suf(1:ncompallFFR))
             read(ifl) ianim_GAS_WDOT(1:ncompFFR),
     &                 ianim_SDOT_suf(1:ncompallFFR),
     &                 ianim_molefr_suf(1:ncompallFFR),
     &                 nsdot,ndeps,nmolfr
           endif
         else
           if(ical_sufFFR==1) then
             read(ifl) SDOT_suf,WDOT_gas,depoeth_spd,molefr_suf,
     &           nsdot,ndeps,nmolfr,num_ratout,blk_thick,nthick
           endif
         endif
         if(ierr/=0) then
           rtrnStat=1
           return
         endif
!DEBUG
         if(ical_sldFFR==1.or.ical_sldFFR==2.or.ical_sldFFR==3) then
            if(ierr/=0) stop 8877
            allocate(iiflag(0:nvrtx),stat=ierr)
            if(ierr/=0) stop 8877
! --- initialize
            cord(:,:)=cord_bak(:,:)
!
            do 300 IIMAT=1,NMATFFR
               IMAT=mat_no(IIMAT)
               IMAT_U=nofld(IIMAT)
               if(ishaft(IIMAT)==3) goto 820 !cycle
               iiflag=0
               alpha=rot_ang(IIMAT)
         write(*,*)
         write(*,'(1x,a,I4)') 'Material NO.= ',IMAT
         write(*,'(1x,a,E10.4)')
     &   ' #### SLINDING ANGLE = ',alpha/3.1415926d0*180.d0
         write(*,'(1x,a,3(E10.4,1x))') 
     &              ' #### SLINDING Origin= ',begin(:,IMAT)
               do 310 ic=1,ncell
                  if(lacell(ic)==IMAT_U) then
                     do i=1,8
                        iv=lvcell(i,ic)
! onishi
                        if(iv.gt.0) then
                           iiflag(iv)=1
                        end if
!                        iiflag(iv)=1
                     enddo
                  endif
 310           continue
c! ***
c               cord(:,:)=cord_bak(:,:)
c! ***
               do 320 iv=1,nvrtx
                  if(iiflag(iv)==1) then
                     xxx=(cord(1,iv)-begin(1,IIMAT))*cos(alpha)
     &                  -(cord(2,iv)-begin(2,IIMAT))*sin(alpha)
                     yyy=(cord(1,iv)-begin(1,IIMAT))*sin(alpha)
     &                  +(cord(2,iv)-begin(2,IIMAT))*cos(alpha)
                     cord(1,iv)=xxx+begin(1,IIMAT)
                     cord(2,iv)=yyy+begin(2,IIMAT)
                    !only z-axis: cord(3,iv)=
                  endif
 320           continue
 300        continue
            deallocate(iiflag)
            icord=1
         endif
!
!         if(ical_sldFFR==1) then
!            deallocate(ishaft,rotati,rot_ang)
!            deallocate(end,begin)
!            deallocate(mat_no)
!         endif
! ---
         allocate(icomp_avex(1:ncompFFR),icomp_rmsx(1:ncompFFR))
         ! Tagami
         allocate(iuvwc_rex(1:ncompFFR))
         iuvwc_rex(:) = 0
         ! imagaT
         
!         allocate(icomp_avex(1:ncompFFR),icomp_rmsx(1:ncompFFR),
!     &            nameVars(1:numVars))
         icomp_avex=0
         icomp_rmsx=0
!
         if(nrnsxFFR.gt.0) then
            allocate(irans_avex(1:nrnsxFFR),irans_rmsx(1:nrnsxFFR),
     &               iuvw_rans_rex(1:nrnsxFFR))
            irans_avex=0
            irans_rmsx=0
            iuvw_rans_rex=0
         endif
!
!
!
! ---    for each wall boundary
         if(numBVars>0) then
            allocate(boundName(1:nbcnd))
            allocate(iwallvar(1:nbcnd))
            allocate(LBC_INDEX(0:nbcnd))
            iwallvar=0
            LBC_INDEX=0

            do i = 1 ,nbcnd
               read(ifl) boundName(i)
            end do
            read(ifl) (iwallvar(nb), nb=1,nbcnd)
            read(ifl) (LBC_INDEX(nb),nb=1,nbcnd)
!
            allocate(LBC_ICV(1:LBC_INDEX(nbcnd)))
            LBC_ICV=0
!           boundary CVs
            do nb=1,nbcnd
               IBFS=LBC_INDEX(nb-1)+1
               IBFE=LBC_INDEX(nb)
               if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
                  read(ifl) (LBC_ICV(IBFL),IBFL=IBFS,IBFE)
               end if
            end do
         else
            numBVars=0
         end if
!
         
!
! --     variable names
         allocate(nameVars(1:(numVars+numBVars)))
!         allocate(nameVars(1:numVars))
         do i = 1 , numVars
            read(ifl) nameVars(i)
            ! Tagami
            char_var = trim(nameVars(i))
            if(char_var(1:6) == "AVE_u_") then
               write(nameVars(i),'(2A)') trim(char_var),";Flux"
            endif
            !imagaT
!           ... if old type name, convert new one.
            select case(trim(adjustl(nameVars(i))))
!            case('pressu')
!               nameVars(i)='Pressure'
            ! Tagami
            case('AVE_u_C2H4')
               nameVars(i)="AVE_u_C2H4;Flux"
            ! imagaT
            case('K')
               nameVars(i)='RANS_K'
            case('eps')
               nameVars(i)='RANS_eps'
            end select
         end do
!
         if(numBVars>0) then
            allocate(nameBVars(1:numBVars))
            do i = 1 , numBVars
               read(ifl) nameBVars(i)
            end do
         end if
!
! ---
!
         timeFV=real(timeFFR)
         write(*,*) ' #### #### Read DATA at Time =',timeFV
         write(*,*) ' #### #### Read DATA at Tstep=',iterFFR
!      
!        initialize
         num = 0
         nnode=0
!        allocate array for strage variables
         if(CPU==1) then
!           *** totdata
            knum=numVars
!           --- for option: -----------------------------
!           *** wall variables
            if(numBVars>0) knum=knum+numBVars
!           *** vars, for extra infomation.
            if(procIDVals) knum=knum+1
! DEBUG
            if(vrtxIDVals) knum=knum+1
!
            if(boundaryYp) knum=knum+1
            if(RelativeVelocity) knum=knum+3
!
            allocate(totdata(1:knum,0:nvrtx),stat=ierr)
            if(ierr.ne.0) then
               write(*,*) 'stop at allocating totdata(',
     &              knum,nvrtx,') in readAndStoreRslFFR.'
               stop
            end if
            totdata(:,:) = 0.0
         end if  ! if(CPU==1)
!        *** ivp, vertex id for every CPU vertexs (1:ncvFFR)
!            ivp(1:nnode)= vertex id
!         allocate(ivp(1:ncvFFR),stat=ierr) !zhang
         allocate(ivp(1:nvrtx),stat=ierr)
         if(ierr.ne.0) then
            stop 'stop at allocating in readAndStoreRslFFR.ivp'
         end if
         ivp = 0
         icv = 0
!
         if(MPprocNum>1) then
!           for HPC
            do i=1,nvrtx  !????
               if(icv>ncvFFR) then
                  write(*,*) ' Warning: does not match CV num',
     &             ' for each CPU, in xxx.frontflow and test.m.part.x.',
     &             ' ncvFFR:',ncvFFR,' icv:',icv
                  stop
               end if
               if(node_id(i)==CPU) then
                  icv=icv+1
                  ivp(icv)=i
               end if
            end do
         else
!           for Single
            do i=1,nvrtx
              icv=icv+1
              ivp(icv)=i
            end do
         end if
         ncvTotal=ncvTotal+icv
!        --- check array size : -----------------------------
         if(ncvTotal>=nvrtx) then
            write(*,*)' Warning: Array size is different.',
     &           nvrtx,'->',ncvTotal
            write(*,*) '*** may be, you are using a buffle boundary.'
            write(*,*) '    you should use FrontFlow native',
     &           'grid file(-gf GF).'
            if(MPprocNum==1) then
               nnode=nvrtx
            else
               nnode=ncvFFR-(ncvTotal-nvrtx)
            endif
         else
            nnode=ncvFFR
         end if
         if(nnode<1) stop 'Error: reading result file.'
!
!
!
! --- ---------------------------------------------------------------
!        Read and store physical data of FFR
! ---    Default
         if( .not. animeStat) then ! anime==.false. #################
!
         allocate(tempReal1(0:ncvFFR),stat=ierr)
         tempReal1(0)=0.d0
         if(ierr.ne.0) stop 'stop at allocating in read result.1'
         !     1.prs
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal1,tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         !     2.pp0
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal1,tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         !     3.rho
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal1,tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         !     4.tmp
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal1,tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         !     5.rmut
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal1,tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         
         !     6.7.8.vel
         allocate(tempReal2(1:3,0:ncvFFR),stat=ierr)
         if(ierr.ne.0) stop 'stop at allocating in read result.vel'
         read(ifl) ((tempReal2(jcv,icv),jcv=1,3),icv=1,ncvFFR)
         num=num+1
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal2(1,1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal2(1,:),tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         num=num+1
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal2(2,1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal2(2,:),tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         num=num+1
         if(icon_cvx==vertex_cen) then
           totdata(num,ivp(1:nnode))=tempReal2(3,1:nnode)
         else
           call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal2(3,:),tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
         endif
         deallocate(tempReal2)

         !     9.yys
         allocate(tempReal2(1:ncompFFR,0:ncvFFR),stat=ierr)
         if(ierr.ne.0) stop 'stop at allocating in read result.yys'
         read(ifl) ((tempReal2(jcv,icv),jcv=1,ncompFFR),icv=1,ncvFFR)
         do i=1, ncompFFR
            num=num+1
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal2(i,:),tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
         end do
         deallocate(tempReal2)
!
         if(ibot/=0) then
           tempReal1(0)=0.d0
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &     tempReal3,tempReal1,tempReal4,lvcell)
           totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
         endif
!
         if(ical_MHD==1.or.ical_MHD==2) then
           N_MHD=32
           allocate(tempReal2(N_MHD,0:ncvFFR),stat=ierr)
           if(ierr.ne.0) stop 'stop at allocating in read result.MHD'
           iii=0
           DO IMHD=1,2
           do i=1,3   
             iii=iii+1
             read(ifl) (tempReal2(iii,icv),icv=1,ncvFFR)
           enddo
           do i=1,3
             iii=iii+1
             read(ifl) (tempReal2(iii,icv),icv=1,ncvFFR)
           enddo
           do i=1,3
             iii=iii+1
             read(ifl) (tempReal2(iii,icv),icv=1,ncvFFR)
           enddo
           do i=1,3
             iii=iii+1
             read(ifl) (tempReal2(iii,icv),icv=1,ncvFFR)
           enddo
           do i=1,3
             iii=iii+1
             read(ifl) (tempReal2(iii,icv),icv=1,ncvFFR)
           enddo
           iii=iii+1
           read(ifl) (tempReal2(iii,icv),icv=1,ncvFFR)
           enddo
           
           IF(iii/=N_MHD) stop 'iii/=N_MHD'
           do i=1,N_MHD
             num=num+1
             if(icon_cvx==vertex_cen) then
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
             else
               call ave_node(nvrtx,ncvFFR,ncell,
     &         tempReal3,tempReal2(i,:),tempReal4,lvcell)
               totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
             endif
           enddo
           deallocate(tempReal2)
         endif
!
         if(ical_sufFFR==1) then
           if(SDOT_suf==1) then
             allocate(tempReal2(1:nsdot,0:ncvFFR),stat=ierr)
             if(ierr.ne.0) stop 'stop at allocating in read surface'
             read(ifl) ((tempReal2(jcv,icv),jcv=1,nsdot),
     &                                                 icv=1,ncvFFR)
             do i=1,nsdot
             num=num+1
             if(icon_cvx==vertex_cen) then
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
             else
               call ave_node(nvrtx,ncvFFR,ncell,
     &         tempReal3,tempReal2(i,:),tempReal4,lvcell)
               totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
             endif
             enddo
             deallocate(tempReal2)
           endif
           if(molefr_suf==1) then
             allocate(tempReal2(1:nmolfr,0:ncvFFR),stat=ierr)
             if(ierr.ne.0) stop 'stop at allocating in read surface'
             read(ifl) ((tempReal2(jcv,icv),jcv=1,nmolfr),
     &                                                 icv=1,ncvFFR)
             do i=1,nmolfr
             num=num+1
             if(icon_cvx==vertex_cen) then
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
             else
               call ave_node(nvrtx,ncvFFR,ncell,
     &         tempReal3,tempReal2(i,:),tempReal4,lvcell)
               totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
             endif
             enddo
             deallocate(tempReal2)
           endif
           if(depoeth_spd==1) then
             allocate(tempReal2(1:ndeps,0:ncvFFR),stat=ierr)
             if(ierr.ne.0) stop 'stop at allocating in read surface'
             read(ifl) ((tempReal2(jcv,icv),jcv=1,ndeps),
     &                                                 icv=1,ncvFFR)
             do i=1,ndeps
             num=num+1
             if(icon_cvx==vertex_cen) then
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
             else
               call ave_node(nvrtx,ncvFFR,ncell,
     &         tempReal3,tempReal2(i,:),tempReal4,lvcell)
               totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
             endif
             enddo
             deallocate(tempReal2)
           endif
           if(blk_thick==1) then
             allocate(tempReal2(1:nthick,0:ncvFFR),stat=ierr)
             if(ierr.ne.0) stop 'stop at allocating in read surface'
             read(ifl) ((tempReal2(jcv,icv),jcv=1,nthick),
     &                                                 icv=1,ncvFFR)
             do i=1,nthick
             num=num+1
             if(icon_cvx==vertex_cen) then
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
             else
               call ave_node(nvrtx,ncvFFR,ncell,
     &         tempReal3,tempReal2(i,:),tempReal4,lvcell)
               totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
             endif
             enddo
             deallocate(tempReal2)
           endif
         endif
!
         if(ical_WDRFFR==1) then
!           allocate(tempReal2(0:ncvFFR,N_INJ+1),stat=ierr)
           do i=1,N_INJ+1
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
           enddo
           
         endif
!
         if(WDOT_gas==1) then
             allocate(tempReal2(1:ncompFFR,0:ncvFFR),stat=ierr)
             if(ierr.ne.0) stop 'stop at allocating in read surface'
             read(ifl) ((tempReal2(jcv,icv),jcv=1,ncompFFR),
     &                                                 icv=1,ncvFFR)
             do i=1, ncompFFR
             num=num+1
             if(icon_cvx==vertex_cen) then
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
             else
               call ave_node(nvrtx,ncvFFR,ncell,
     &         tempReal3,tempReal2(i,:),tempReal4,lvcell)
               totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
             endif
             enddo
             deallocate(tempReal2)
         endif
         if(num_ratout>0) then
             allocate(tempReal2(1:num_ratout,0:ncvFFR),stat=ierr)
             if(ierr.ne.0) stop 'stop at allocating in read surface'
             read(ifl) ((tempReal2(jcv,icv),jcv=1,num_ratout),
     &                                                 icv=1,ncvFFR)
             do i=1,num_ratout
             num=num+1
             if(icon_cvx==vertex_cen) then
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
             else
               call ave_node(nvrtx,ncvFFR,ncell,
     &         tempReal3,tempReal2(i,:),tempReal4,lvcell)
               totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
             endif
             enddo
             deallocate(tempReal2)
         endif
! --- -----------------------------------------------------------
! ---    ieul2ph
         if (ieul2ph == 1) then
            !   pressure2
            num=num+1
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal1,tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
            !   rho2
            num=num+1
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal1,tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
            !   tmp2
            num=num+1
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal1,tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
            !   rmut2
            num=num+1
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal1,tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
            !   vel2
            
            allocate(tempReal2(1:3,0:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.vel2'
            read(ifl) ((tempReal2(jcv,icv),jcv=1,3),icv=1,ncvFFR)
            num=num+1
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal2(1,1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal2(1,:),tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
            num=num+1
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal2(2,1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal2(2,:),tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
            num=num+1
            if(icon_cvx==vertex_cen) then
              totdata(num,ivp(1:nnode))=tempReal2(3,1:nnode)
            else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal2(3,:),tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
            endif
            deallocate(tempReal2)

            !     yys2
            allocate(tempReal2(1:ncompFFR,0:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.yys2'
            read(ifl) ((tempReal2(jcv,icv),jcv=1,ncompFFR),icv=1,ncvFFR)
            do i=1,ncompFFR
               num=num+1
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal2(i,:),tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            enddo
            deallocate(tempReal2)
         end if
!
! ---    aks
         if (nrnsxFFR>0) then
            allocate(tempReal2(1:nrnsxFFR,0:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.aks'
            do i=1 , nrnsxFFR
            read(ifl) (tempReal2(i,icv),icv=1,ncvFFR)
            enddo
!            read(ifl) ((tempReal2(jcv,icv),jcv=1,nrnsxFFR),icv=1,ncvFFR)
            do i=1,nrnsxFFR
               num=num+1
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal2(i,:),tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end do
            deallocate(tempReal2)
         end if
! ---  potential
         if(npotnxFFR>0) then
           allocate(tempReal2(1:npotnxFFR,0:ncvFFR),stat=ierr)
           if(ierr.ne.0) stop 'stop at allocating in read npotnxFFR'
           do i=1,npotnxFFR
!!!!!           read(ifl) ((tempReal2(jcv,icv),jcv=1,npotnxFFR),icv=1,ncvFFR)
            read(ifl) (tempReal2(i,icv),icv=1,ncvFFR)
           enddo
           do i=1,npotnxFFR
           num=num+1
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
           else
              call ave_node(nvrtx,ncvFFR,ncell,
     &        tempReal3,tempReal2(i,:),tempReal4,lvcell)
              totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
           enddo
           deallocate(tempReal2)
         endif
!
         if(radflag.EQ.1) then
!	   allocate(tempReal1(0:ncvFFR),stat=ierr)
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
!           deallocate(tempReal1)
         endif
!

         if(NFLIDxFFR>0) then
           allocate(tempReal2(1:NFLIDxFFR,0:ncvFFR),stat=ierr)
           if(ierr.ne.0) stop 'stop at allocating in read NFLIDxFFR'
           read(ifl) ((tempReal2(jcv,icv),jcv=1,NFLIDxFFR),icv=1,ncvFFR)
           do i=1,NFLIDxFFR
           num=num+1
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal2(i,:),tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
           enddo
           deallocate(tempReal2)
         endif
!
         if(idrdpFFR==compFFR) then
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
         endif
!
         if(MPprocNum>1) then
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
         endif
!
         if(ical_prtx==1) then
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
         endif
!
         if(u_func(2)>=1.or.u_func(3)>=1) then
           num=num+1
           read(ifl) (tempReal1(icv),icv=1,ncvFFR)
           if(icon_cvx==vertex_cen) then
             totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
           else
             call ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
             totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
           endif
         endif
!
         deallocate(tempReal1)
!
! --- -----------------------------------------------------------
! ---    Ave,RMS
         if (num .lt. numVars) then
            read(ifl) iuvw_ave_rms_rex,ip_avex,imu_avex,
     &                it_avex,it_rmsx,ip_rmsx,imu_rmsx,
     &                iuvwt_rex,imaxminx,ista_momeryx,ical_prtx
            if(ncompFFR>0) then
               read(ifl) (icomp_avex(i),i=1,ncompFFR),
     &                   (icomp_rmsx(i),i=1,ncompFFR)
               ! Tagami
               read(ifl) (iuvwc_rex(i),i=1,ncompFFR)
               ! imagaT
            end if
            if(nrnsxFFR>0) then
               read(ifl) (irans_avex(i),i=1,nrnsxFFR),
     &                   (irans_rmsx(i),i=1,nrnsxFFR),
     &                   (iuvw_rans_rex(i),i=1,nrnsxFFR)
            end if
!            read(ifl) iuvw_ave_rms_rex,
!     &           ip_avex,it_avex,
!     &           (icomp_avex(i)    ,i=1,ncompFFR),
!     &           (irans_avex(i)    ,i=1,nrnsxFFR),
!     &           it_rmsx,ip_rmsx,
!     &           (icomp_rmsx(i)    ,i=1,ncompFFR),
!     &           (irans_rmsx(i)    ,i=1,nrnsxFFR),
!     &           iuvwt_rex,
!     &           (iuvw_rans_rex(i) ,i=1,nrnsxFFR)
!
!
            allocate(tempReal1(0:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.rms'
            tempReal1(0)=0.d0
            if (iuvw_ave_rms_rex==1) then
               do i=1,9
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               end do
            end if
!
            if (ip_avex==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if
!
            if(ip_rmsx==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if
!Mu
            if (imu_avex==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
              else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            endif
!
            if(imu_rmsx==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            endif
!
            if(imaxminx==1) then
              do i=1,8
              num=num+1
              read(ifl) (tempReal1(icv),icv=1,ncvFFR)
              enddo
            endif
!
            if (it_avex==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if
!
            if (it_rmsx==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if
!
            if (iuvwt_rex==1) then
               do i=1,3
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               enddo
            end if
!
            do i=1,ncompFFR
               if(icomp_avex(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               endif
               if(icomp_rmsx(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               endif
            end do
            ! Tagami
            do j = 1, ncompFFR
               if(iuvwc_rex(j) == 1) then
                  do i = 1, 3
                     num=num+1
                     read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                     if(icon_cvx==vertex_cen) then
                       totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                     else
                       call ave_node(nvrtx,ncvFFR,ncell,
     &                 tempReal3,tempReal1,tempReal4,lvcell)
                       totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                     endif
                  enddo
               endif
            enddo
            ! imagaT
!
            do i=1,nrnsxFFR
               if(irans_avex(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               endif
               if(irans_rmsx(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               endif
            end do

            if (ical_prtx==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if            

!            if(MPprocNum>1) then
!               read(hmcpu,*) tempI1   ! total CPU number
!            end if
!
            deallocate(tempReal1)

         end if
         end if ! ###############################################
! --- -----------------------------------------------------------
! ---    Animation
         if(animeStat) then ! anime==.true.
! DEBUG
            read(ifl) ianim,ianim_uvw,ianim_p,ianim_r,ianim_t
            if(nrnsxFFR>0) then
               read(ifl) (irans_avex(i),i=1,nrnsxFFR)
            end if
            if(ncompFFR>0) then
               read(ifl) (icomp_avex(i),i=1,ncompFFR)
            end if
!            read(ifl) ianim,ianim_uvw,ianim_p,
!     &           ianim_r,ianim_t,(irans_avex(i),i=1,nrnsxFFR),
!     &           (icomp_avex(i),i=1,ncompFFR)

            allocate(tempReal1(0:ncvFFR),stat=ierr)
            tempReal1(0)=0.d0
            if(ierr.ne.0) stop 'stop at allocating in read result.anim'

            !     anime.vel
            if(ianim_uvw==1) then
               allocate(tempReal2(1:3,0:ncvFFR),stat=ierr)
               if(ierr.ne.0) stop 'stop at allocating in result.anm.vel'
               read(ifl) ((tempReal2(jcv,icv),jcv=1,3),icv=1,ncvFFR)
               num=num+1
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal2(1,1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal2(1,:),tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
               num=num+1
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal2(2,1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal2(2,:),tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
               num=num+1
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal2(3,1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal2(3,1:nnode),tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
               deallocate(tempReal2)
            end if

            !     anime.prs
            if (ianim_p==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if
            !     anime.rho
            if (ianim_r==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if
            !     anime.tmp
            if (ianim_t==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               if(icon_cvx==vertex_cen) then
                 totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               else
                 call ave_node(nvrtx,ncvFFR,ncell,
     &           tempReal3,tempReal1,tempReal4,lvcell)
                 totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
               endif
            end if

            !     anime.aks
            do i=1,nrnsxFFR
               if (irans_avex(i)==1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               end if
            end do

            !     anime.yys
            do i=1,ncompFFR
               if (icomp_avex(i)==1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                  else
                    call ave_node(nvrtx,ncvFFR,ncell,
     &              tempReal3,tempReal1,tempReal4,lvcell)
                    totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                  endif
               end if
            end do
!
            if(ical_sufFFR==1) then
              if(nsdot>0) then
                do i=1,nsdot
                num=num+1
                read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                if(icon_cvx==vertex_cen) then
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                else
                  call ave_node(nvrtx,ncvFFR,ncell,
     &            tempReal3,tempReal1,tempReal4,lvcell)
                  totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                endif
                enddo
              endif

              if(nmolfr>0) then
                do i=1,nmolfr
                num=num+1
                read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                if(icon_cvx==vertex_cen) then
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                else
                  call ave_node(nvrtx,ncvFFR,ncell,
     &            tempReal3,tempReal1,tempReal4,lvcell)
                  totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                endif
                enddo
              endif

              if(ndeps>0) then
                do i=1,ndeps
                num=num+1
                read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                if(icon_cvx==vertex_cen) then
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                else
                  call ave_node(nvrtx,ncvFFR,ncell,
     &            tempReal3,tempReal1,tempReal4,lvcell)
                  totdata(num,ivp(1:nnode))=tempReal3(1:nnode) 
                endif
                enddo
              endif
!
              if(nthick>0) then
                do i=1,nthick
                num=num+1
                read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                if(icon_cvx==vertex_cen) then
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                else
                  call ave_node(nvrtx,ncvFFR,ncell,
     &            tempReal3,tempReal1,tempReal4,lvcell)
                  totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                endif
                enddo
              endif
!
              do i=1,ncompFFR
                if(ianim_GAS_WDOT(i)==1) then
                   num=num+1
                   read(ifl,ERR=2000) (tempReal1(icv),icv=1,ncvFFR)
                   if(icon_cvx==vertex_cen) then
                     totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
                   else
                     call ave_node(nvrtx,ncvFFR,ncell,
     &               tempReal3,tempReal1,tempReal4,lvcell)
                     totdata(num,ivp(1:nnode))=tempReal3(1:nnode)
                   endif
                endif
              enddo
              goto 2001
 2000         stop 777
 2001         continue
              deallocate(ianim_GAS_WDOT,ianim_SDOT_suf,ianim_molefr_suf)
            endif

!
            deallocate(tempReal1)
!
         end if
!
! --
!
         if(numBVars>0) then
            allocate(tempReal1(0:LBC_INDEX(nbcnd)),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.bnd'
            tempReal1(0)=0.d0
            bnum=0
!            allocate(bnddata(1:numBVars,1:nvrtx))
!            bnddata=0.d0
!
! ---       Wall normal
!
            do nb=1,nbcnd
               IBFS=LBC_INDEX(nb-1)+1
               IBFE=LBC_INDEX(nb)
               if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
               do i=1,3
                  read(ifl) (tempReal1(IBFL),IBFL=IBFS,IBFE)
                  if(icon_cvx==vertex_cen) then
                    totdata(num+i,ivp(LBC_ICV(IBFS:IBFE)))=
     &                 tempReal1(IBFS:IBFE)
                   else
                     stop 'stop at BC'
                   endif
               enddo
               end if
            end do
            num=num+3
!
!
! ---       yplus
!
            num=num+1
            do nb=1,nbcnd
               IBFS=LBC_INDEX(nb-1)+1
               IBFE=LBC_INDEX(nb)
               if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
                  read(ifl) (tempReal1(IBFL),IBFL=IBFS,IBFE)
                  if(icon_cvx==vertex_cen) then
                    totdata(num,ivp(LBC_ICV(IBFS:IBFE)))=
     &                 tempReal1(IBFS:IBFE)
                  else
                  endif
               end if
            end do
!
! ---       wall shear force
! 
            do nb=1,nbcnd
               IBFS=LBC_INDEX(nb-1)+1
               IBFE=LBC_INDEX(nb)
               if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
               do i=1,3
                  read(ifl) (tempReal1(IBFL),IBFL=IBFS,IBFE)
                  if(icon_cvx==vertex_cen) then
                    totdata(num+i,ivp(LBC_ICV(IBFS:IBFE)))=
     &                 tempReal1(IBFS:IBFE)
                  else
                  endif
               enddo
               end if
            end do
            num=num+3
!
! ---      wall pres force
!
            do nb=1,nbcnd
               IBFS=LBC_INDEX(nb-1)+1
               IBFE=LBC_INDEX(nb)
               if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
               do i=1,3
                  read(ifl) (tempReal1(IBFL),IBFL=IBFS,IBFE)
                  if(icon_cvx==vertex_cen) then
                    totdata(num+i,ivp(LBC_ICV(IBFS:IBFE)))=
     &                 tempReal1(IBFS:IBFE)
                  else
                  endif
               enddo
               end if
            end do
            num=num+3
c!
c! ---      wall total force
c!
c            do nb=1,nbcnd
c               IBFS=LBC_INDEX(nb-1)+1
c               IBFE=LBC_INDEX(nb)
c               if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
c               do i=1,3
c                  read(ifl) (tempReal1(IBFL),IBFL=IBFS,IBFE)
c                  totdata(num+i,ivp(LBC_ICV(IBFS:IBFE)))=
c     &                 tempReal1(IBFS:IBFE)
c               enddo
c               end if
c            end do
c            num=num+3
!
!           conbine wall variables
            do i = 1 , numBVars
               nameVars(numVars+i)=nameBVars(i)
            end do
            numVars=numVars+numBVars
!            nameBVars=''
!            numBVars=0
!
            deallocate(nameBVars)
            deallocate(tempReal1)
!
         end if
! ---    if rsCloseStat==.true., close result file
!        (at last on animation steps).
         if(rsCloseStat) then         
            close(ifl)
         end if


! --- -- option:procID --------------------------------------
         if(procIDVals) then
            num=num+1
            newVarName='CPU_ID'
            call addVal(newVarName)
            totdata(num,ivp(1:nnode))=real(CPU)
         end if
! --- -- option:vrtxID --------------------------------------
         if(vrtxIDVals) then
            num=num+1
            newVarName='vertexID'
            call addVal(newVarName)
            totdata(num,ivp(1:nnode))=real(ivp(1:nnode))
         end if
! --- -- Set wall velocity to 0 -----------------------------
         if(wallZero) then
            do l=1,numVars
               tmpflag=.false.
               select case(trim(adjustl(nameVars(l))))
               case('velo_u;Velocity')
                  tmpflag=.true.
               case('velo_v')
                  tmpflag=.true.
               case('velo_w')
                  tmpflag=.true.
               case('aver_u;Average_velocity')
                  tmpflag=.true.
               case('aver_v')
                  tmpflag=.true.
               case('aver_w')
                  tmpflag=.true.
c               case('RMS_u')
c                  tmpflag=.true.
c               case('RMS_v')
c                  tmpflag=.true.
c               case('RMS_w')
c                  tmpflag=.true.
c               case('Re_uv')
c                  tmpflag=.true.
c               case('Re_vw')
c                  tmpflag=.true.
c               case('Re_uw')
c                  tmpflag=.true.
               end select

               if(tmpflag .and. numBVars>0) then
               do nb=1,nbcnd
                  IBFS=LBC_INDEX(nb-1)+1
                  IBFE=LBC_INDEX(nb)
                  if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
                     write(*,*) '*** set ',trim(adjustl(nameVars(l))),
     &                      ' to zero on:',trim(adjustl(boundName(nb)))
!
                     totdata(l  ,ivp(LBC_ICV(IBFS:IBFE)))=0.0d0
                  end if
               end do
               elseif(tmpflag .and. numBVars==0) then
                  do k=1,NBOUND
                  write(*,*) '*** set ',trim(adjustl(nameVars(l))),
     &                   ' to zero on:',trim(adjustl(SFBOUN(k)))
                  do j=1,NBFS
                  if(IBFACE(j)==k) then
                    do i=1,4
                     if(IFFACE(i,j)>0) then
                     totdata(l  ,ivp(IFFACE(i,j)))=0.0d0
                     end if
                    end do
                  end if
                  end do
                  end do
               end if
!
            end do
         end if
! --- -- Create reference vector(for relative velocity) -------------
! DEBUG
         if( (ical_sldFFR==1.or.ical_sldFFR==2.or.ical_sldFFR==3) 
     &                    .and. RelativeVelocity ) then
!         if(ical_sldFFR==1.and.RelativeVelocity) then
!            num=num+1
            newVarName='rpmx;rpm_velocity'
            call addVal(newVarName)
            newVarName='rpmy'
            call addVal(newVarName)
            newVarName='rpmz'
            call addVal(newVarName)
!
            allocate(iiflag(0:nvrtx),stat=ierr)
            if(ierr/=0) stop 8877
            do 200 IIMAT=1,NMATFFR
               IMAT=mat_no(IIMAT)
               IMAT_U=nofld(IIMAT)
               if(ishaft(IIMAT)/=1) cycle
               iiflag=0
               if(MPprocNum>1) then ! HPC
                  do 210 ic=1,ncell
                  if(lacell(ic)==IMAT_U) then
                     do i=1,8
                        iv=lvcell(i,ic)
                        if(node_id(iv)==CPU .and. iv.gt.0) then
                           iiflag(iv)=1
                        end if
                     enddo
                  endif
 210              continue
               else ! Single CPU
                  do 211 ic=1,ncell
                  if(lacell(ic)==IMAT_U) then
                     do i=1,8
                        iv=lvcell(i,ic)
                        if(iv.gt.0) then
                           iiflag(iv)=1
                        end if
                     enddo
                  endif
 211              continue
               endif
               do 220 iv=1,nvrtx
                  if(iiflag(iv)==1) then
                     unit(1)=end(1,IIMAT)-begin(1,IIMAT)
                     unit(2)=end(2,IIMAT)-begin(2,IIMAT)
                     unit(3)=end(3,IIMAT)-begin(3,IIMAT)
                     radi(1)=cord(1,iv)-begin(1,IIMAT)
                     radi(2)=cord(2,iv)-begin(2,IIMAT)
                     radi(3)=cord(3,iv)-begin(3,IIMAT)
                     call AXB_UNIT_C(unit,radi,vr)
                     dr=radi(1)*unit(1)+radi(2)*unit(2)+radi(3)*unit(3)
                     radi(1)=radi(1)-dr*unit(1)
                     radi(2)=radi(2)-dr*unit(2)
                     radi(3)=radi(3)-dr*unit(3)
                     dr=dsqrt( radi(1)*radi(1)
     &                        +radi(2)*radi(2)
     &                        +radi(3)*radi(3) )
                     v0(1)=dr*rotati(IIMAT)*vr(1)
                     v0(2)=dr*rotati(IIMAT)*vr(2)
                     v0(3)=dr*rotati(IIMAT)*vr(3)
!
                     totdata(num+1,iv)=v0(1)
                     totdata(num+2,iv)=v0(2)
                     totdata(num+3,iv)=v0(3)
!                     totdata(num+1,ivp(iv))=v0(1)
!                     totdata(num+2,ivp(iv))=v0(2)
!                     totdata(num+3,ivp(iv))=v0(3)
               endif
 220           continue
 200        continue
            deallocate(iiflag)
            num=num+3
! DEBUG
         else if( (ical_sldFFR/=1.and.ical_sldFFR/=2.and.
     &             ical_sldFFR/=3)
     &                           .and. RelativeVelocity ) then
!         else if(ical_sldFFR/=1.and.RelativeVelocity) then
            write(*,*) 'WRN: you can not create rpm reference velocity',
     &           'on NO SLIDING case.'
         endif
! --- -----------------------------------------------------------
!        after-processing
         rtrnStat=0

         if (num .ne. numVars) then
            write(*,*) ' #### 3D ERR: num:',num,'.ne. numVars:',numVars
            rtrnStat=1
            return
         end if
!
!        if last CPU
         if (CPU == MPprocNum) then
            if(MPprocNum>1 .and. .not.animeStat) then ! anime==.false.
               deallocate(node_id)
            end if
         else
            deallocate(icomp_avex,icomp_rmsx,nameVars)
            ! Tagami
            deallocate(iuvwc_rex)
            ! imagaT
         end if
         if(nrnsxFFR.gt.0) then
            deallocate(irans_avex,irans_rmsx,iuvw_rans_rex)
         endif

!        for next CPU
         deallocate(ivp)
         if(numBVars>0) then
            deallocate(boundName,iwallvar,LBC_INDEX,LBC_ICV)
            numBVars=0
         end if
!
 820     continue
         if(ical_sldFFR==1.or.ical_sldFFR==2.or.ical_sldFFR==3
     &    .or.ical_sldFFR==4) then
!            if(CPU==MPprocNum) then
              deallocate(ishaft,rotati,rot_ang)
              deallocate(nofld,end,begin)
              deallocate(mat_no)
!            endif
         endif

 800  end do                    ! do CPU=1,MPprocNum
!
! --- -------------------------------------------------------

      rtrnStat=0
      return
      end subroutine readAndStoreRslFFR
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE AXB_UNIT_C(A,B,C)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------
!     A X B = C 
!---------------------------------------
      real*8 ,intent(in)   :: A(3),B(3)
      real*8 ,intent(OUT)  :: C(3)
      real*8 :: D
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
!
      C(3)=(A(1)*B(2)-A(2)*B(1))
      C(2)=(A(3)*B(1)-A(1)*B(3))
      C(1)=(A(2)*B(3)-A(3)*B(2))
      D=DSQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))+SML
      C(3)=C(3)/D
      C(2)=C(2)/D
      C(1)=C(1)/D
!
      RETURN
      END subroutine AXB_UNIT_C
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ave_node(nvrtx,ncvFFR,ncell,
     &       tempReal3,tempReal1,tempReal4,lvcell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------
!     A X B = C 
!---------------------------------------
      integer,intent(in)   :: nvrtx,ncvFFR
      real*8 ,intent(in)   :: tempReal1(0:ncvFFR)
      real*8 ,intent(inOUT):: tempReal3(0:nvrtx)
      real*8 ,intent(inOUT):: tempReal4(0:nvrtx)
      integer,intent(in)   :: lvcell(1:8,1:ncell)
      integer              :: IC,IV,I
!
      tempReal3(:)=0.d0
      tempReal4(:)=0.d0
      do IC=1,ncvFFR
         do I=1,8
         IV=lvcell(I,IC)
         tempReal3(IV)=tempReal3(IV)+tempReal1(IC)
         tempReal4(IV)=tempReal4(IV)+1.d0*dble(IV)/dble(IV+1.d-20)
         enddo
      enddo
      tempReal3(1:nvrtx)=tempReal3(1:nvrtx)/tempReal4(1:nvrtx)
!
      RETURN
      END subroutine ave_node
!
