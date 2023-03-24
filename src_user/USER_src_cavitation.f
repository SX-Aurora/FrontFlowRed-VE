!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine  user_src_cavitation
     &   (FrontSTR,Axis_DIR,BC_no,CAVI_B,
     &    BC_name,deltt,iter,time,IMAT_U,
     &    POWER,vel,p,dens,dens2,T_WALL,vol,area,area_vect,
     &    T,void,vap_mass_fraction,eva_mass,
     &    IBFS,IBFE,MXSSFBC,IBFL1,LBC_SSF,
     &    MXCVFAC,SFAREA,SFCENT,mv_surf,ampl)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     FrontSTR   : FLAG 
!     Axis_DIR   : Force direction
!     CAVI_B     : BC number on which Cavitation occur 
!     POWER      : % of Cavitation
!     BC_no      : BC number
!     BC_name    : Boundary name
!     deltt      : dt
!     iter       : time step
!     IMAT_U     : Material No. of BC_name
!     u,v,w      : Velocities of BC_name   [m/s]
!     u2,v2,w2   : Velocities of BC_name   [m/s]
!     p          : pressure  [N/m^3]
!     dens       : Density at BC_name [kg/m^3]
!     dens2      : Density at BC_name [kg/m^3]
!     T_WALL     : wall Temperature at BC_name [K]
!     T          : Temperature at ICV [K]
!     T2         : Temperature at ICV [K]
!     eva_mass   : [kg/s/m^3]
!---------------------------------------------------------------------
!
!--------------------------------------------------------------------- 
      use module_initial ,only : rho0,rho02
      use module_species,only  : r_wm,gascns
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)      :: FrontSTR
      integer,intent(in)      :: BC_no
      integer,intent(in)      :: CAVI_B
      integer,intent(in)      :: Axis_DIR
      integer,intent(in)      :: IBFS,IBFE,MXSSFBC,MXCVFAC,IBFL1
      real*8 ,intent(in)      :: POWER
      character(*),intent(in) :: BC_name
      real*8 ,intent(in)      :: deltt,time,vol,area,area_vect(3)
      integer,intent(in)      :: IMAT_U,iter
      real*8 ,intent(inout)   :: vel(3),p,dens,dens2
      real*8 ,intent(in)      :: T,T_WALL
      real*8 ,intent(in)      :: void
      real*8 ,intent(inout)   :: vap_mass_fraction
      real*8 ,intent(out)     :: eva_mass
!
      real*8 ,intent(in)      :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)      :: SFCENT(3, MXCVFAC)
      INTEGER,INTENT(IN)      :: LBC_SSF(  MXSSFBC)
      REAL*8 ,INTENT(OUT)     :: mv_surf,ampl
      real*8,external :: Pvap
!---------------------------------------------------------------------
! --- [local entities]
!---------------------------------------------------------------------      
      real*8 :: dum2,dum3,dum4,fact=1.D-3,dis,big,fact1=1.D-3
      real*8,save :: dum1(3)
      integer :: idum,INI_FLAG=0,IBFL,NP=1
      integer :: ierr=0,ios=0,ifl=30,NNN,II,ICFL,IIN
      real*8,save,allocatable :: xyz(:,:),amp(:)
      integer,save,allocatable :: ISF(:)
      real*8 :: fff=47.d3,vvv,pi=3.1415!*2.d0
     &      ,rho_v=0.01708d0
     &          ,omega,
     &          uuu_n
      
!---------------------------------------------------------------------      
      IF(INI_FLAG==0.and.FrontSTR==1) then 
        INI_FLAG=1
        open(ifl,file='amp.dat',form='formatted',status='old',
     &       action='read',iostat=ios)
        read(ifl,*) NNN
        ALLOCATE(xyz(NNN,3),amp(NNN),ISF(IBFS:IBFE))
        do 1000 II=1,NNN
        read(ifl,*) IDUM,dum2,dum1(1),dum1(2),dum1(3)
        xyz(II,1)=dum1(1)*fact1
        xyz(II,2)=dum1(2)*fact1
        xyz(II,3)=dum1(3)*fact1
        amp(II)=dum2*fact
 1000   enddo
        close(ifl)
!
        do IBFL=IBFS,IBFE
        big=1.d10
        ICFL=LBC_SSF(IBFL) 
        dum1(:)=SFCENT(:,ICFL) 

        do II=1,NNN 
        dis=sqrt((dum1(1)-xyz(II,1))**2
     &          +(dum1(2)-xyz(II,2))**2
     &          +(dum1(3)-xyz(II,3))**2)
        if(dis<big) then
          big=dis
          IIN=II
        endif
        enddo
        ISF(IBFL)=IIN 
        enddo 
        if(Axis_DIR==1) then
          dum1(:)=0.d0
          dum1(1)=1.d0
        elseif(Axis_DIR==2) then
          dum1(:)=0.d0
          dum1(2)=1.d0
        elseif(Axis_DIR==3) then
          dum1(:)=0.d0
          dum1(3)=1.d0
        else
          call FFRABORT(1,'ERR: MUST: Axis = 1,2 or 3')
        endif 
      endif 
!
      if(BC_no==CAVI_B) then 
        if(FrontSTR==0) then 
          uuu_n=1.d-4*fff/pi 
          mv_surf=-POWER**NP*uuu_n*area*abs(area_vect(Axis_DIR))
          vap_mass_fraction=1.d0 
          eva_mass=0.d0 
        else
!
          uuu_n=abs(amp(ISF(IBFL1)))*fff/pi/2.d0   !pi 
          ampl=amp(ISF(IBFL1)) 
!--------------------------------------------(1) 
!          mv_surf=-POWER**NP*uuu_n*area*abs(area_vect(Axis_DIR)) 
!          vap_mass_fraction=1.d0 
!          eva_mass=0.d0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(2) 
!          vap_mass_fraction=1.d0 
!          vel(1)=0.d0
!          vel(2)=0.d0
!          vel(3)=0.d0
!          vel(Axis_DIR)=-uuu_n*POWER**NP*area_vect(Axis_DIR) 
!          mv_surf=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (3)
          vap_mass_fraction=1.d0 
          dum3=abs( 
     &       dum1(1)*area_vect(1)
     &      +dum1(2)*area_vect(2)
     &      +dum1(3)*area_vect(3)) 
          mv_surf=-uuu_n*area        !*dum3
          dum2=Pvap(T)/(gascns*r_wm(1)*T) 
          eva_mass=POWER**NP*dum2*abs(mv_surf)/vol*dum3
          vel(1)=uuu_n*area_vect(1)
          vel(2)=uuu_n*area_vect(2)
          vel(3)=uuu_n*area_vect(3)
        endif
      endif
!
      return
   
!--------------------------------------------------------------------
      end subroutine user_src_cavitation





!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine  user_src_cavitation_Olympus1
     &   (FrontSTR,Axis_DIR,BC_no,CAVI_B,
     &    BC_name,deltt,iter,time,IMAT_U,
     &    POWER,u,v,w,p,dens,dens2,T_WALL,vol,area,area_vect,
     &    T,void,vap_mass_fraction,eva_mass,
     &    IBFS,IBFE,MXSSFBC,IBFL1,LBC_SSF,
     &    MXCVFAC,SFAREA,SFCENT,mv_surf,ampl)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     FrontSTR   : FLAG 
!     Axis_DIR   : Force direction
!     CAVI_B     : BC number on which Cavitation occur
!     POWER      : % of Cavitation
!     BC_no      : BC number
!     BC_name    : Boundary name
!     deltt      : dt
!     iter       : time step
!     IMAT_U     : Material No. of BC_name
!     u,v,w      : Velocities of BC_name   [m/s]
!     u2,v2,w2   : Velocities of BC_name   [m/s]
!     p          : pressure  [N/m^3]
!     dens       : Density at BC_name [kg/m^3]
!     dens2      : Density at BC_name [kg/m^3]
!     T_WALL     : wall Temperature at BC_name [K]
!     T          : Temperature at ICV [K]
!     T2         : Temperature at ICV [K]
!     voit1,voit2: 
!     eva_mass   : [kg/s/m^3]
!---------------------------------------------------------------------
!
!--------------------------------------------------------------------- 
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)      :: FrontSTR
      integer,intent(in)      :: BC_no
      integer,intent(in)      :: CAVI_B
      integer,intent(in)      :: Axis_DIR
      integer,intent(in)      :: IBFS,IBFE,MXSSFBC,MXCVFAC,IBFL1
      real*8 ,intent(in)      :: POWER
      character(*),intent(in) :: BC_name
      real*8 ,intent(in)      :: deltt,time,vol,area,area_vect(3)
      integer,intent(in)      :: IMAT_U,iter
      real*8 ,intent(in)      :: u,v,w,p,dens,dens2
      real*8 ,intent(in)      :: T,T_WALL
      real*8 ,intent(in)      :: void
      real*8 ,intent(inout)   :: vap_mass_fraction
      real*8 ,intent(out)     :: eva_mass
!
      real*8 ,intent(in)      :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)      :: SFCENT(3, MXCVFAC)
      INTEGER,INTENT(IN)      :: LBC_SSF(  MXSSFBC)
      REAL*8 ,INTENT(OUT)     :: mv_surf,ampl
!---------------------------------------------------------------------
! --- [local entities]
!---------------------------------------------------------------------      
      real*8 :: dum1(3),dum2,dum3,dum4,fact=1.D-3,dis,big,fact1=1.D-3
      integer :: idum,INI_FLAG=0,IBFL
      integer :: ierr=0,ios=0,ifl=30,NNN,II,ICFL,IIN
      real*8,save,allocatable :: xyz(:,:),amp(:)
      integer,save,allocatable :: ISF(:)
      real*8 :: fff=47.d3,vvv,pi=2.d0*3.1415,omega
!---------------------------------------------------------------------      
      IF(INI_FLAG==0.and.BC_no==CAVI_B) then
        INI_FLAG=1
        open(ifl,file='amp.dat',form='formatted',status='old',
     &       action='read',iostat=ios)
        read(ifl,*) NNN
        ALLOCATE(xyz(NNN,3),amp(NNN),ISF(IBFS:IBFE))
        do 1000 II=1,NNN
        read(ifl,*) IDUM,dum2,dum1(1),dum1(2),dum1(3)
        xyz(II,1)=dum1(1)*fact1
        xyz(II,2)=dum1(2)*fact1
        xyz(II,3)=dum1(3)*fact1
        amp(II)=dum2*fact
 1000   enddo
        close(ifl)
!
        do IBFL=IBFS,IBFE
        big=1.d10
        ICFL=LBC_SSF(IBFL) 
        dum1(:)=SFCENT(:,ICFL)
        do II=1,NNN
        dis=sqrt((dum1(1)-xyz(II,1))**2
     &          +(dum1(2)-xyz(II,2))**2
     &          +(dum1(3)-xyz(II,3))**2)
        if(dis<big) then
          big=dis
          IIN=II
        endif
        enddo
        ISF(IBFL)=IIN
        enddo
      endif

      if(BC_no==CAVI_B) then
        if(FrontSTR==0) then
          eva_mass=abs(area_vect(Axis_DIR))*
     &     POWER*0.01708/3.14*0.00005d0/vol*area*47.D3
        else
          eva_mass=
     &     0.01708*POWER*amp(ISF(IBFL1))/3.14/vol*area*fff
!!!!!     &     POWER*0.01708/3.14*amp(ISF(IBFL1))/vol*area*47.D3  !OK
          eva_mass=0.d0
!
! --------------------------------
!
          omega=pi*fff
          vvv=omega*amp(ISF(IBFL1))*area*0.01708/3.14
!*cos(pi*iter/360.d0)
!          mv_surf=-vvv*abs(area_vect(Axis_DIR))*POWER
          mv_surf=-vvv*POWER
          ampl=amp(ISF(IBFL1))
          vap_mass_fraction=1.d0
        endif
      endif
!
      return
   
!--------------------------------------------------------------------
      end subroutine user_src_cavitation_Olympus1
