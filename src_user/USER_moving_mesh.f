!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE user_moving_mesh(ical_mvmsh,
     &        iter,time,
     &        NVRTX,MXVRTX,MXFACE,MXCELL,
     &        movflg,cord,lbface,lacell,lvface,lvcell,lcface,
     &        NBOUND,NBFS,SFBOUN,NFBOUN,IFFACE,icall)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_usersub,only : usr_dum1,usr_dum2,usr_dum3,
     &                          usr_dum4,usr_dum5,usr_dum6
      use module_io,  only : ifle,ifll,gdScale
!
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      integer,intent(in)    :: ical_mvmsh,iter,icall
      real*8 ,intent(in)    :: time
      integer,intent(in)    :: NVRTX,MXVRTX,MXFACE,MXCELL
      integer,intent(inout) :: movflg(MXVRTX)
      real*8 ,intent(inout) :: cord(3,MXVRTX)
      integer,intent(in)    :: NBOUND,NBFS
      integer,intent(in)    :: NFBOUN(NBOUND)
      integer,intent(in)    :: IFFACE(4,NBFS,NBOUND)
      integer,intent(in)    :: lbface(2,MXFACE)
      integer,intent(in)    :: lvcell(8,MXCELL)
      integer,intent(in)    :: lvface(4,MXFACE)    
      integer,intent(in)    :: lcface(2,MXFACE)
      integer,intent(in)    :: lacell(  MXCELL)
      character*80,intent(in) :: SFBOUN(0:NBOUND)
!_____________________________________________________________________
!
!     ical_mvmsh :
!               =1  Running [Moving] solver 
!               =2  Checking CV-GEOM of [Moving] mesh no-running solver
!               =3  Viewing [Moving] & [Removing/Adding] mesh
!               =4  Checking CV-GEOM of [Moving] & [Removing/Adding] mesh 
!               =5  Running [Moving] & [Removing/Adding] mesh solver
!               =0  No [Moving] & [Removing/Adding] Mesh solver
!_____________________________________________________________________
!     
!     movflg: Grid type or layer-number
!             movflg(IV)=0 : no moving vertex
!             movflg(IV)>0 : layer-number (moving grid)
!             movflg(IV)<0 : removing/add 
!_____________________________________________________________________
!
! --- [LOCAL ENTITIES]
!
      integer :: mb,iv,I,J,layer,layera,iloc,nb,ic1,ic2,iv3

      integer :: i_umat1=1,i_umat2=2
!_____________________________________________________________________
!
! --- [user array]
!
      real*8 ,parameter :: SML=1.d-10
      real*8 ,parameter :: Ztopini=0.359d-3
      real*8 ,parameter :: RRR=0.04d0,rho=1.d0/3.d0,rpm=3600.d0
      real*8 ,parameter :: LINER=0.095d0,pi=3.1415926d0
      integer :: iter_ctl,ierr
      real*8,save,allocatable :: z_ini(:)
      real*8            :: Ztop,the,rps=rpm/60.d0
      real*8            :: the0=pi+pi/6.d0,dum1,dum2
!
!---------------------------------------------------------------------
! --- (Step-1: icall=1) First call for layer defination 
!---------------------------------------------------------------------
!
      if(icall==1) then
!
        allocate(z_ini(MXVRTX),stat=ierr)
        if(ierr/=0) then
           write(*,*) ' ERR: allocate in user_moving_mesh'
           stop
        endif
!
        movflg(:)=0
        z_ini(:)=0.d0
!
        do iv=1,NVRTX
        z_ini(iv)=cord(3,iv)
        enddo
      endif
!
!-----------------------
! --- moving mesh
!-----------------------
      if(icall==2) then 
! 
       the=the0+rps*2.d0*pi*time 
       dum1=RRR*((1.d0-cos(the))-1.d0/rho
     &     *(1.d0-sqrt(1.d0-rho**2*sin(the)**2))) 
       Ztop=0.08d0-dum1
       do iv=1,NVRTX
       dum2=LINER-(LINER-Ztop)*(LINER-z_ini(iv))/(LINER-Ztopini)
       cord(3,iv)=dum2
       enddo
!
       print*,'the=,Ztop=',the*360/2/pi,Ztop,dum1
      endif
!_____________________________________________________________________
      return
      end SUBROUTINE user_moving_mesh
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

