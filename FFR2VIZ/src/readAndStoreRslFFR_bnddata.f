! *** ******************************************************************
      subroutine readAndStoreRslFFR(MPprocNum,fileName,hpcprefix,
     &                              procIDVals,boundaryYp,bndYpWallfn,
     &                              animeStat,rsOpenStat,rsCloseStat,
     &                              rtrnStat)
! *** ******************************************************************
      use FVdata, only : timeFV
      use FFRdata,only : nvrtx,ncell,numVars,nameVars,
     &                   numBVars,nameBVars,bnddata,
     &                   iterFFR,ncvFFR,ncvfacFFR,ncompFFR,
     &                   nrnsxFFR,ieul2ph,timeFFR,
     &                   ical_sldFFR,NMATFFR,ishaft,rotati,rot_ang,
     &                   end,begin,mat_no,
     &                   lvcell,lacell,cord,totdata,node_id,
     &                   LBC_INDEX,LBC_ICV,npotnxFFR
      use FFRreaddata, only : iuvw_ave_rms_rex,iuvwt_rex,
     &                        ip_avex,it_avex,ip_rmsx,it_rmsx,
     &                        icomp_avex,irans_avex,
     &                        icomp_rmsx,irans_rmsx,iuvw_rans_rex,
     &                        ivar_bnd
      use FFRreaddata, only : ianim,ianim_uvw,ianim_p,ianim_r,ianim_t
!
      implicit none
!     
      integer,     intent(in)  :: MPprocNum
      character*80,intent(in)  :: fileName
      character*80,intent(in)  :: hpcprefix
      logical,     intent(in)  :: procIDVals
      logical,     intent(in)  :: boundaryYp
      character*80,intent(in)  :: bndYpWallfn
      logical,     intent(in)  :: rsOpenStat
      logical,     intent(in)  :: rsCloseStat
      logical,     intent(out) :: animeStat
      integer,     intent(out) :: rtrnStat
!
      integer :: CPU
      integer,save :: ifl
      integer :: nnode,ncvTotal,icv,jcv
      integer,allocatable :: ivp(:)
      integer :: i,j,l,num,knum,tempI1, ierr
      integer :: IIMAT,IMAT,iv,ic
      real(8) :: up,yp,utau,rho,rmu,rnu,xxx,yyy
      real(8),parameter :: SML=1.0d-25
      real(8),allocatable :: tempReal1(:),tempReal2(:,:)
! sliding
      integer :: IIMAT,IMAT,iv,ic
      real(8) :: alpha
      integer,allocatable :: iiflag(:)
! boundary variables
      integer :: nb,nbcnd,bnum,IBFL,IBFS,IBFE
      character*80,allocatable :: boundName(:)
!
      character*80  :: newVarName
      character*255 :: filePath
      character*255 :: hpcPath    ! extern function
      integer :: icord=0
!
! --- -------------------------------------------------------
!     read 'reslt.frontflow' for every CPU, respectively
! --- -------------------------------------------------------
      ncvTotal=0
      if(MPprocNum+30>231) then
         write(*,*)'ABEND:Too many CPUs data exists.'
         stop
      end if
      do CPU=1,MPprocNum
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
         if(.not.animeStat) then ! anime==.false.
            read(ifl)iterFFR,timeFFR,NMATFFR,ncvFFR,ncvfacFFR,ncompFFR,
     &           nrnsxFFR,numVars,numBVars,ieul2ph
!            read(ifl)iterFFR,timeFFR,NMATFFR,ncvFFR,ncvfacFFR,ncompFFR,
!     &           nrnsxFFR,numVars,ieul2ph
            read(ifl)ical_sldFFR
            if(ical_sldFFR==1.or.ical_sldFFR==2.or.ical_sldFFR==4) then
               allocate(ishaft(1:NMATFFR),rotati(1:NMATFFR))
               allocate(end(1:3,1:NMATFFR),begin(1:3,1:NMATFFR))
               allocate(mat_no(1:NMATFFR))
               allocate(rot_ang(1:NMATFFR))
               do I=1,NMATFFR
                  read(ifl) IIMAT,IMAT,ishaft(IMAT),rotati(IMAT),
     &                 end(1:3,IMAT),begin(1:3,IMAT),rot_ang(IMAT)
                  mat_no(I)=IMAT
               enddo
            endif
            if(ierr/=0) then
               rtrnStat=1
               return
            end if
! ---
! ---    sliding mesh  ---------------------------------
! ---
            if((ical_sldFFR==1.or.ical_sldFFR==2.or.ical_sldFFR==4)
     &        .and.icord==0) then
               if(ierr/=0) stop 8877
               allocate(iiflag(0:nvrtx),stat=ierr)
               if(ierr/=0) stop 8877
               do 100 IIMAT=1,NMATFFR
                  IMAT=mat_no(IIMAT)
                  if(ishaft(IMAT)/=1) cycle
                  iiflag=0
                  alpha=rot_ang(IMAT)
         write(*,*) ' #### SLINDING ANGLE= ',alpha/3.1415926d0*180.d0

                  do 110 ic=1,ncell
                     if(lacell(ic)==IMAT) then
                        do i=1,8
                           iv=lvcell(i,ic)
! onishi
                           if(iv.gt.0) then
                              iiflag(iv)=1
                           end if
!                           iiflag(iv)=1
                        enddo
                     endif
 110              continue
                  do 120 iv=1,nvrtx
                     if(iiflag(iv)==1) then
                        xxx=cord(1,iv)*cos(alpha)
     &                     -cord(2,iv)*sin(alpha)
                        yyy=cord(1,iv)*sin(alpha)
     &                     +cord(2,iv)*cos(alpha)
                        cord(1,iv)=xxx
                        cord(2,iv)=yyy
!      print*,'D',alpha,cos(alpha),-sin(alpha),sin(alpha),cos(alpha)

                       !only z-axis: cord(3,iv)=
                     endif
 120              continue
 100           continue
               deallocate(iiflag)
               icord=1
            endif
            if(ical_sldFFR==1.or.ical_sldFFR==2.or.ical_sldFFR==4) then
               deallocate(ishaft,rotati,rot_ang)
               deallocate(end,begin)
               deallocate(mat_no)
            endif


         else                    ! anime==.true.
! animation
            read(ifl, iostat=ierr) iterFFR,timeFFR,ncvFFR,
     &           ncvfacFFR,ncompFFR,nrnsxFFR,
     &           numVars !,ieul2ph
            if(ierr/=0) then
               rtrnStat=-1
               animeStat=.false.
               close(ifl)
               return
            end if
         end if





! ---
         allocate(icomp_avex(1:ncompFFR),icomp_rmsx(1:ncompFFR),
     &            nameVars(1:numVars))
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
         do i = 1 , numVars
            read(ifl) nameVars(i)
! onishi
!           ... if old type name, convert new one.
            select case(trim(adjustl(nameVars(i))))
!            case('pressu')
!               nameVars(i)='Pressure'
            case('K')
               nameVars(i)='RANS_K'
            case('eps')
               nameVars(i)='RANS_eps'
            end select
         end do
! --
! boundary variables
         if(numBVars>0) then
            read(ifl) nbcnd
            allocate(boundName(1:nbcnd))
            allocate(nameBVars(1:numBVars))
            do i = 1 ,nbcnd
            read(ifl) boundName(i)
            end do
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
!           *** vars, for extra infomation.
            if(procIDVals) knum=knum+1
            if(boundaryYp) knum=knum+1
            allocate(totdata(1:knum,1:nvrtx),stat=ierr)
            if(ierr.ne.0) then
               write(*,*) 'stop at allocating totdata(',
     &              knum,nvrtx,') in readAndStoreRslFFR.'
               stop
            end if
            totdata = 0.0
         end if  ! if(CPU==1)

!        *** ivp, vertex id for every CPU vertexs (1:ncvFFR)
!            ivp(1:nnode)= vertex id
         allocate(ivp(1:ncvFFR),stat=ierr)
         if(ierr.ne.0) then
            stop 'stop at allocating in readAndStoreRslFFR.ivp'
         end if
         ivp = 1
!         ivp = 0
         icv=0
!
         if(MPprocNum>1) then
!           for HPC
            do i=1,nvrtx
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
         if(ncvTotal>nvrtx) then
            write(*,*)' Warning: Array size is different.',
     &           nvrtx,'->',ncvTotal
            write(*,*) '*** may be, you are using a buffle boundary.'
            write(*,*) '    you should use FrontFlow native',
     &           'grid file(-gf GF).'
            if(MPprocNum==1) then
               nnode=nvrtx
            else
               nnode=ncvFFR-(ncvTotal-nvrtx)
            end if
         else
            nnode=ncvFFR
         end if
         if(nnode<1) stop 'Error: reading result file.'


! --- ---------------------------------------------------------------
!        Read and store physical data of FFR
! ---    Default
         if( .not. animeStat) then ! anime==.false. #################
!
         allocate(tempReal1(1:ncvFFR),stat=ierr)
         if(ierr.ne.0) stop 'stop at allocating in read result.1'

         !     1.prs
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         !     2.pp0
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         !     3.rho
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         !     4.tmp
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
         !     5.rmut
         num=num+1
         read(ifl) (tempReal1(icv),icv=1,ncvFFR)
         totdata(num,ivp(1:nnode))=tempReal1(1:nnode)

         !     6.7.8.vel
         allocate(tempReal2(1:3,1:ncvFFR),stat=ierr)
         if(ierr.ne.0) stop 'stop at allocating in read result.vel'
         read(ifl) ((tempReal2(jcv,icv),jcv=1,3),icv=1,ncvFFR)

!         do icv=1,ncvFFR  !ZHANG8
!            if(iiflag(icv)==1) then
!            xxx=tempReal2(1,Icv)*cos(alpha)
!     &         -tempReal2(2,Icv)*sin(alpha)
!            yyy=tempReal2(1,Icv)*sin(alpha)
!     &         +tempReal2(2,Icv)*cos(alpha)
!            tempReal2(1,Icv)=xxx
!            tempReal2(2,Icv)=yyy
!         ENDIF
!         enddo
!         deallocate(iiflag)



         num=num+1
         totdata(num,ivp(1:nnode))=tempReal2(1,1:nnode)
         num=num+1
         totdata(num,ivp(1:nnode))=tempReal2(2,1:nnode)
         num=num+1
         totdata(num,ivp(1:nnode))=tempReal2(3,1:nnode)
         deallocate(tempReal2)

         !     9.yys
         allocate(tempReal2(1:ncompFFR,1:ncvFFR),stat=ierr)
         if(ierr.ne.0) stop 'stop at allocating in read result.yys'
         read(ifl) ((tempReal2(jcv,icv),jcv=1,ncompFFR),icv=1,ncvFFR)
         do i=1, ncompFFR
            num=num+1
            totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
         end do
         deallocate(tempReal2)

! --- -----------------------------------------------------------
! ---    ieul2ph
         if (ieul2ph == 1) then
            !   pressure2
            num=num+1
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            !   rho2
            num=num+1
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            !   tmp2
            num=num+1
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            !   rmut2
            read(ifl) (tempReal1(icv),icv=1,ncvFFR)
            totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            !   vel2
            allocate(tempReal2(1:3,1:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.vel2'
            read(ifl) ((tempReal2(jcv,icv),jcv=1,3),icv=1,ncvFFR)
            num=num+1
            totdata(num,ivp(1:nnode))=tempReal2(1,1:nnode)
            num=num+1
            totdata(num,ivp(1:nnode))=tempReal2(2,1:nnode)
            num=num+1
            totdata(num,ivp(1:nnode))=tempReal2(3,1:nnode)
            deallocate(tempReal2)

            !     yys2
            allocate(tempReal2(1:ncompFFR,1:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.yys2'
            read(ifl) ((tempReal2(jcv,icv),jcv=1,ncompFFR),icv=1,ncvFFR)
            do i=1,ncompFFR
               num=num+1
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
            enddo
            deallocate(tempReal2)
         end if
!
! ---    aks
         if (nrnsxFFR>0) then
            allocate(tempReal2(1:nrnsxFFR,1:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.aks'
            read(ifl) ((tempReal2(jcv,icv),jcv=1,nrnsxFFR),icv=1,ncvFFR)
            do i=1 , nrnsxFFR
               num=num+1
               totdata(num,ivp(1:nnode))=tempReal2(i,1:nnode)
            end do
            deallocate(tempReal2)
         end if
!
         deallocate(tempReal1)
!

! --- -----------------------------------------------------------
! ---    Ave,RMS
         if (num .lt. numVars) then

            read(ifl) iuvw_ave_rms_rex,
     &           ip_avex,it_avex,
     &           (icomp_avex(i)    ,i=1,ncompFFR),
     &           (irans_avex(i)    ,i=1,nrnsxFFR),
     &           it_rmsx,ip_rmsx,
     &           (icomp_rmsx(i)    ,i=1,ncompFFR),
     &           (irans_rmsx(i)    ,i=1,nrnsxFFR),
     &           iuvwt_rex,
     &           (iuvw_rans_rex(i) ,i=1,nrnsxFFR)
!
!
            allocate(tempReal1(1:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.rms'
            if (iuvw_ave_rms_rex==1) then
               do i=1,9
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               end do
            end if
!
            if (ip_avex==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
!
            if (ip_rmsx==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
!
            if (ip_avex==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
!
            if (ip_rmsx==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
!
            if (it_avex==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
!
            if (it_rmsx==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
!
            if (iuvwt_rex==1) then
               do i=1,3
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               enddo
            end if
!
!
            do i=1,ncompFFR
               if(icomp_avex(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               endif
               if(icomp_rmsx(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               endif
            end do
!
            do i=1,nrnsxFFR
               if(irans_avex(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               endif
               if(irans_avex(i).eq.1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               endif
            end do
!
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
            read(ifl) ianim,ianim_uvw,ianim_p,
     &           ianim_r,ianim_t,(irans_avex(i),i=1,nrnsxFFR),
     &           (icomp_avex(i),i=1,ncompFFR)

            allocate(tempReal1(1:ncvFFR),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.anim'

            !     anime.vel
            if (ianim_uvw==1) then
               allocate(tempReal2(1:3,1:ncvFFR),stat=ierr)
               if(ierr.ne.0) stop 'stop at allocating in result.anm.vel'
               read(ifl) ((tempReal2(jcv,icv),jcv=1,3),icv=1,ncvFFR)
               num=num+1
               totdata(num,ivp(1:nnode))=tempReal2(1,1:nnode)
               num=num+1
               totdata(num,ivp(1:nnode))=tempReal2(2,1:nnode)
               num=num+1
               totdata(num,ivp(1:nnode))=tempReal2(3,1:nnode)
               deallocate(tempReal2)
            end if

            !     anime.prs
            if (ianim_p==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
            !     anime.rho
            if (ianim_r==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if
            !     anime.tmp
            if (ianim_t==1) then
               num=num+1
               read(ifl) (tempReal1(icv),icv=1,ncvFFR)
               totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
            end if

            !     anime.aks
            do i=1,nrnsxFFR
               if (irans_avex(i)==1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               end if
            end do

            !     anime.yys
            do i=1,ncompFFR
               if (icomp_avex(i)==1) then
                  num=num+1
                  read(ifl) (tempReal1(icv),icv=1,ncvFFR)
                  totdata(num,ivp(1:nnode))=tempReal1(1:nnode)
               end if
            end do

            deallocate(tempReal1)

         end if
!
! --
! DEBUG
         if(numBVars>0) then
            bnum=0
            allocate(ivar_bnd(1:nbcnd))
            allocate(bnddata(1:numBVars,1:nvrtx))
            allocate(LBC_INDEX(0:nbcnd))
            ivar_bnd=0
            bnddata=0.d0
            LBC_INDEX=0
!
! ---       flags for each boundary
!
            read(ifl) (ivar_bnd(nb), nb=1,nbcnd)
            read(ifl) (LBC_INDEX(nb),nb=1,nbcnd)
!
            allocate(LBC_ICV(1:LBC_INDEX(nbcnd)))
            allocate(tempReal1(1:LBC_INDEX(nbcnd)),stat=ierr)
            if(ierr.ne.0) stop 'stop at allocating in read result.bnd'
            LBC_ICV=0

! ---       boundary CVs
            print *, 'AAAAAAA 00',LBC_INDEX
            do nb=1,nbcnd
               IBFS=LBC_INDEX(nb-1)+1
               IBFE=LBC_INDEX(nb)
               if(ivar_bnd(nb).eq.1) then
               read(ifl) (LBC_ICV(IBFL),IBFL=IBFS,IBFE)
               end if
            end do
!
! ---      Wall normal
!
            do nb=1,nbcnd
               IBFS=LBC_INDEX(nb-1)+1
               IBFE=LBC_INDEX(nb)
               if(ivar_bnd(nb).eq.1) then
               do i=1,3
                  bnum=bnum+1
                  read(ifl) (tempReal1(IBFL),IBFL=IBFS,IBFE)
                  bnddata(bnum,ivp(LBC_ICV(IBFS:IBFE)))=
     &                 tempReal1(IBFS:IBFE)
               enddo
               end if
            end do
!
! DEBUG
            print *, 'AAAAAAA 02'
            deallocate(ivar_bnd,LBC_INDEX,LBC_ICV)
            deallocate(tempReal1)

         end if

! ---    if rsCloseStat==.true., close result file
!        (at last on animation steps).
         if(rsCloseStat) then         
            close(ifl)
         end if

! onishi
! --- -- option:procID --------------------------------------
         if(procIDVals) then
            num=num+1
            newVarName='CPU_ID'
            call addVal(newVarName)
            totdata(num,ivp(1:nnode))=real(CPU)
         end if

! --- -----------------------------------------------------------
!        after-processing
         rtrnStat=0

         if (num .ne. numVars) then
            write(*,*) ' #### BC ERR: num:',num,'.ne. numVars:',numVars
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
         end if
         if(nrnsxFFR.gt.0) then
            deallocate(irans_avex,irans_rmsx,iuvw_rans_rex)
         endif

!        for next CPU
         deallocate(ivp)

      end do      ! do CPU=1,MPprocNum
!
! --- -------------------------------------------------------

      rtrnStat=0
      return
      end subroutine readAndStoreRslFFR
