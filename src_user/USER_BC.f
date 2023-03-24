!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_BC(mph,iter,time,IBF,no,prt_S_iter,
     &           nrans,NCOMP,
     &           area_inl,area_cent,iters,itere,
     &           ICVL,t_cell,p_cell,v_cell,aks_cell,rho_cell,yys_cell,
     &           ydis,velinl,p,t,keps,yys,mv_surf)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$1$$$$$$$
      use module_usersub,only : usr_dum1,usr_dum2,usr_dum3,
     &                          usr_dum4,usr_dum5,usr_dum6      
!
!     User subroutine for [inlet],[wall],[outlet] boundary condition
!
!     IBF    : Globle BC Index (in)
!     no     : BC no. defined in fflow.ctl
!     area_inl(1:3) : unit vector of IBF (outward)
!     area_inl(4) :  area of IBF
!     area_cent(3): coordinate of IBF center
!     ydis   : Minimum distance to wall from IBF
!
! 
!     velinl : Velocity of u,v,w at IBF 
!     t      : if t='D' => then t=temperature[K] 
!            : if t='transfer' => then t=temperature[K] 
!            : if t='N' => then t=[W/m^2]
!     p      : Pressure and temperature at IBF
!     keps   : K & e 
!     yys    : Mass fraction of components 
!     prt_S_iter : particle injection start step 
!=====================================================================
!
      implicit none  
!
! --- [dummy arguments]
!
      integer,intent(in)  :: iter,mph,iters,itere,prt_S_iter
      real*8 ,intent(in)  :: time
      integer,intent(in)  :: IBF,no,nrans,NCOMP
      real*8 ,intent(in)  :: ydis
      real*8 ,intent(in)  :: area_inl(4),area_cent(3)
      integer,intent(in)  :: ICVL
      real*8 ,intent(in)  :: t_cell
      real*8 ,intent(in)  :: p_cell
      real*8 ,intent(in)  :: v_cell(3)
      real*8 ,intent(in)  :: aks_cell(nrans)
      real*8 ,intent(in)  :: rho_cell
      real*8 ,intent(in)  :: yys_cell(NCOMP)
      real*8 ,intent(out) :: velinl(3),p,t,keps(nrans),yys(NCOMP),
     &                       mv_surf
!
! --- [local entities]
!
      real(8),parameter :: um=5.0d0  !1.d0,5.d0           ! common in all direction
      real(8),parameter :: zh=74.60,zmin=0.d0   ! common in all direction
!      real(8),parameter :: angle = 337.5d0   ! wind direction NNW
      real(8),parameter :: angle = 235.0d0   ! wind direction SW
!      real(8),parameter :: angle = 157.5d0   ! wind direction SSE
!      real(8),parameter :: angle =  45.0d0   ! wind direction NE
!
!      real(8),parameter :: angle =  67.5d0   ! wind direction ENE  !?
!      real(8),parameter :: angle = 315.0d0   ! wind direction NW    !?
!-----------------------------------------------------------------
!      real(8),parameter :: angle = 82.d0   ! wind direction W
!      real(8),parameter :: angle = 90.d0   ! wind direction E
!      real(8),parameter :: angle = 360.d0   ! wind direction N
!      real(8),parameter :: angle = 180.d0   ! wind direction S
!      real(8),parameter :: angle = 45.d0   ! wind direction NE
!      real(8),parameter :: angle = 135.d0   ! wind direction SE
!      real(8),parameter :: angle = 225.d0   ! wind direction SW
!      real(8),parameter :: angle = 315.d0   ! wind direction NW
!-----------------------------------------------------------------
      real(8) :: dum1,dum2,z
      integer :: I,J,K
      integer,save :: iflag=0
      real(8) :: ang_x,ang_y
      real(8) :: theta
      real(8),parameter :: pi = 3.1415926535d0
!
! NNE-ENE
      if(no > 4) return

      theta = pi*angle/180.d0
      ang_x = -sin(theta)        
      ang_y = -cos(theta)        
      !
      z=area_cent(3)
      if(z-zmin<=0.d0) then
        dum1=0.d0
      else
         dum1=((z-zmin)/zh)**0.2
      endif
!
      if(dum1>1.d0)  dum1=1.d0
      dum2=dum1*um
!
      velinl(1)= dum2*ang_x
      velinl(2)= dum2*ang_y
      velinl(3)=0.d0
!
      return
      end subroutine user_BC
!

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_BC_VECTOR
     &  (mph,iter,time,deltt,iters,itere,timee,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,DISINL,SFAREA,SFCENT,locmsh,
     &  tmp,rho,vel,yys,aks,prs,XTA,
     &  prsbnd,aksbnd,
     &  velbnd,tmpbnd,yysbnd,htcbnd,mtcbnd,angbnd,SHUTFL)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments] 
!
      use module_dimension
      use module_constant
      use module_boundary,only : chkncomp,chknrans,nobcnd,
     &                           wdbcnd,adbcnd,iscomp,israns,distrb,
     &                           LBC_INDEX,MAT_BCIDX,nbcnd,kdbcnd,
     &                           kdintr,kdilet,kdolet,
     &                           kxnone,kdbuff,kdintr,kdfire,rotwall,
     &                           kdcvd,sizein,kdshutr,
     &                           endw,beginw,omega,
     &                           kdpres,openout,kdstag,
     &                           masflg,masflx,imasflg,sumare,
     &                           vofang,shut_strt,ktEWF,kyEWF
      use module_dimnsn,only   : xredc,yredc,zredc
      use module_usersub,only  : usrno,usryes,inlusr,oulusr,walusr,
     &                           stcusr
      use module_scalar,  only : rns_scl,ike
      use module_model,   only : icaltb,ke,ke_low,RNG,CHEN,ical_vect,
     &                           ical_prt
      use module_metrix,  only : yys_cell=>rcomp,aks_temp
      use module_species, only  : r_wm,gascns
      use module_dimension
      use module_constant
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: iter,mph,iters,itere
      real*8 ,intent(in)  :: time,deltt,timee
      integer,intent(in)  :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)  :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)  :: LCYCSF    (   MXSSFBC)
      logical,INTENT(IN)  :: mat_cal   (   0:MXMAT)
      real*8 ,intent(out) :: prsbnd    (   MXSSFBC)
      real*8 ,intent(out) :: velbnd    (   MXSSFBC,3)
      real*8 ,intent(out) :: tmpbnd    (   MXSSFBC)
      real*8 ,intent(out) :: yysbnd    (   MXSSFBC,MXcomp)
      real*8 ,intent(out) :: aksbnd    (   MXSSFBCR,MXrans)
      real*8 ,intent(out) :: htcbnd    (   MXSSFBC)
      real*8 ,intent(out) :: mtcbnd    (   MXSSFBC)
      real*8 ,intent(in)  :: DISINL    (   MXSSFBC)
      real*8 ,intent(in)  :: SFAREA    ( 4,MXCVFAC)
      real*8 ,intent(in)  :: SFCENT    ( 3,MXCVFAC)
      real*8 ,intent(in)  :: vel   (  MXALLCV,3)
      real*8 ,intent(in)  :: rho   (  MXALLCV)
      real*8 ,intent(in)  :: tmp   (  MXALLCV)
      real*8 ,intent(in)  :: yys   (  MXALLCV,MXCOMP)
      real*8 ,intent(in)  :: aks   (  MXALLCVR,mxrans)
      real*8 ,intent(in)  :: prs   (  MXALLCV)
      integer,intent(in)  :: locmsh(  MXSSFBC)
      real*8 ,intent(out) :: XTA   (  MXCVFAC)
      real*8 ,intent(out) :: angbnd(  MXSSFBC_VOF)
      integer,intent(inout)  :: SHUTFL (MXSSFBC_SHT)
!
! --- [local entities] 
!
      real*8  :: dum1,dum2,dum3,dum4
      integer :: idum1,idum2,idum3,idum4,nb,IIMAT,IBFS,IBFE,kd,no
      integer :: ibfl,icfl
      logical :: lkdwall
      real(8),parameter :: t1 = 6.d-3
      real(8),parameter :: u1 = 1200.5584d0
      
!
! --- ---------------------------------------------------------------
!
      do 4000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      kd=kdbcnd(0,nb)
      no=nobcnd(nb)
      if(no==1) then
        if(time < t1 ) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          velbnd(IBFL,1)= 0.d0
          velbnd(IBFL,2)= 0.d0
          velbnd(IBFL,3)=-SFAREA(3,ICFL)*u1/t1*time
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          velbnd(IBFL,1)= 0.d0
          velbnd(IBFL,2)= 0.d0
          velbnd(IBFL,3)=-SFAREA(3,ICFL)*u1
          enddo
        endif
      endif

 4000 enddo

      
      return
      end subroutine user_BC_VECTOR
