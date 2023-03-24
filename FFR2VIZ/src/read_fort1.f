!     read_fort1_boundary
!     read_fort1_hpc
!     read_fort1_fluid
! *** *********************************************************
      subroutine read_fort1_boundary
! *** *********************************************************
      use FFRreaddata,only : cntlnam
      use SCRYUpre, only   : NREGN,NKND
      use FLUENTdata ,only : F_NREGN,F_NKND,F_UVW
      use FFRdata, only    : NBOUND,SFBOUN
!
      implicit none
!
      integer :: i,j,k,ios,ros,rbcnd
      integer :: ifl=20
      integer,parameter :: lc=10
      integer,parameter :: mbctim=10
      integer,parameter :: mcomp=200, mrans=50
      character(lc) :: kind
      character(lc) :: vel,temp,conc,rans
      integer :: no,side,ntime,profile,open_air
      real    :: time(mbctim)
      real    :: p(mbctim),u(mbctim),v(mbctim),w(mbctim)
      real    :: t(mbctim),ys(mcomp*mbctim)
      real    :: u2(mbctim),v2(mbctim),w2(mbctim),htc2,mtc2,vsgm2(2)
      real    :: t2(mbctim),ys2(mcomp*mbctim)
      real    :: aks(mrans*mbctim)
      real    :: cycle,htc,mtc,vsgm(2),rot(4)
      character*80 :: name,name2
      character*1  :: axis
      real    :: xoffset,yoffset,zoffset
      namelist /boundary/ no,side,kind,ntime,
     &                    vel,temp,conc,rans,
     &                    time,profile,p,u,v,w,t,ys,aks,
     &                    cycle,htc,mtc,vsgm,rot,
     &                    name,name2,axis,xoffset,yoffset,zoffset,
     &                    u2,v2,w2,t2,ys2,htc2,mtc2,vsgm2,
     &                    open_air
!
! --- -----------------------------------------------------
      inquire(file=cntlnam,number=ios)
      if(ios>=0) then
         write(*,*)
     &        cntlnam,' has been already opened as number ',ios
      end if
!
      open(ifl,file=cntlnam,form='formatted',iostat=ios)
! --- ----------------------------------------------------------
! onishi
!     check NBOUND in FFRdata, and NREGN in module SCRYUpre
      NREGN=max(NREGN,NBOUND)
      if(NREGN<=0) then
         write(*,*) 'Warning: there is no boundary.'
      end if
!
      allocate(  NKND(1:NREGN)) ;     NKND=' '  ! SCRYU
      allocate(F_NKND(1:NREGN)) ;   F_NKND=' '  ! FLUENT
      allocate(F_UVW(1:3,1:NREGN)); F_UVW=0.0   ! FLUENT

!
      ros=1
      rbcnd=0
      rewind ifl
      do while(ros>=0)
         read(ifl,boundary,iostat=ros)
         if( ros.lt.0 ) exit
         rbcnd=rbcnd+1
!
         do i = 1 , NREGN
            if (trim(adjustl(SFBOUN(i))) == trim(adjustl(name))) then
                 NKND(i) = trim(adjustl(kind))
               F_NKND(i) = trim(adjustl(kind))
               if(NKND(i)(1:4)=='inle') then
                  F_UVW(1,i) = DBLE(u(1))
                  F_UVW(2,i) = DBLE(v(1))
                  F_UVW(3,i) = DBLE(w(1))
               end if
            end if
         end do
      end do
!
! --- ----------------------------------------------------------
      close(ifl)

!     check boundary num.
      if(rbcnd/=NREGN) then
         write(*,*) 'Warning: boundary num=',NREGN,'is not matched ',
     &        'in grid-file and ',cntlnam,':',rbcnd
      end if

!
      return
      end subroutine read_fort1_boundary



! *** *********************************************************
      subroutine read_fort1_hpc(PartfileName)
! *** *********************************************************
      use FFRreaddata,only : cntlnam
!      use FFRdata
!      use AVSdata

      implicit none
c
      character(80),intent(out) :: PartfileName
c
      character(80)     :: mtsfil,boundary,vertex,wall,subdir
      character(80)     :: hpc_boundary,hpc_comm,hpc_vertex,ucdfile
      character(80)     :: hpc_anim
      character(80)     :: hpc_source,hpc_result
      character(80)     :: hpc_initial,hpc_restart,hpc_force
c
      integer           :: NPE,walflg,ucdflg
      integer :: i,ii,ios
      integer :: ifl=20
c
      namelist /hpc_cntl/
     &        NPE,mtsfil,boundary,vertex,wall,subdir,
     &        hpc_vertex,hpc_boundary,hpc_comm,hpc_force,
     &        hpc_source,hpc_result,hpc_initial,hpc_restart,
     &        hpc_anim,
     &        walflg,ucdflg,ucdfile
!
      npe = 1 ; walflg = 0 ; ucdflg = 0
      hpc_boundary = ' '
      hpc_comm     = ' '
      hpc_vertex   = ' '
      mtsfil       = 'test.m'  ! default
      wall         = ' '
      subdir       = ' '
      ucdfile      = ' '
c
      write(*,*) ' #### Reading hpc_cntl in ',cntlnam
      open(ifl,file=cntlnam,status='old',iostat=ios)
      if(ios/=0) then
        write(*,*)'Cannot open file :',cntlnam
      end if
      read(ifl,hpc_cntl,iostat=ios)
      if(ios/=0) then
         write(*,*) ' Error:reading ',cntlnam,'(hpc_cntl):',mtsfil
         mtsfil='test.m'
      end if
      close(ifl)
c
      PartfileName=mtsfil
!      write(*,*) ' #### Node number obtained from ',cntlnam,'='
!     &           ,nvrtx
c
      return
      end subroutine read_fort1_hpc




! *** *********************************************************
      subroutine read_fort1_fluid(rmu)
! *** *********************************************************
      use FFRreaddata,only : cntlnam
!
      implicit none
!
      real*8,intent(out) :: rmu

      integer       :: IMAT_U,Poisson_start
      character(10) :: muopt,wall_distance
      real*8        :: Prandtl,Schmidt,mu,t0,c
      logical       :: visc,diff_ty,diff_ke
      real*8        :: unrelax_P=1.0,unrelax_T=1.0
      character(10) :: cnv_vv,cnv_ty,cnv_ke
      character(10) :: limitr_vv,limitr_ty,limitr_ke
      real*8        :: Prandtl2,Schmidt2,mu2,t02,c2
!
      integer,parameter :: iundef=-huge(1)
      real*8 ,parameter :: undef =-huge(1.)
      real*8  :: rnck=0.3d0,bldfct=0.8D0
!
      integer :: ios
      integer :: ifl=20
!
      namelist /fluid/ IMAT_U,muopt,
     &                 wall_distance,Poisson_start,
     &                 cnv_vv,cnv_ty,cnv_ke,
     &                 limitr_vv,limitr_ty,limitr_ke,
     &                 rnck,bldfct,
     &                 visc,diff_ty,diff_ke,
     &                 unrelax_P,unrelax_T,
     &                 Prandtl ,Schmidt ,mu ,t0 ,c,
     &                 Prandtl2,Schmidt2,mu2,t02,c2

      IMAT_U = iundef
      mu     = undef

      open(ifl,file=cntlnam,status='old',iostat=ios)
      if(ios/=0) then
        write(*,*)'Cannot open file :',cntlnam
      end if
!     ----- fluid ----
      do while(.true.)
         read(ifl,fluid,iostat=ios)
         if(ios/=0) then
            write(*,*) ' Error:reading ',cntlnam,'(fluid)'
            exit
         end if
         if(IMAT_U==1) exit
      end do
      close(ifl)

!     multi material is not supported.
      if(IMAT_U/=1) then
         write(*,*) '### error : data error'
         write(*,*) 'lack of data : IMAT_U'
      end if

      rmu=mu

      return
      end subroutine read_fort1_fluid
