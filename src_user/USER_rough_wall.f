!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_rough_wall
     &	(icall,no,iter,deltt,time,
     &   u_cell,v_cell,w_cell,
     &   x_unit,y_unit,z_unit,
     &   area,x_cent,y_cent,z_cent,dis_cell,
     &   t_cell,t_wall,rho_cell,rmu_cell,rmut_cell,roughness)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!------------------------
! --- [module arguments]
!------------------------
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_model,only   : idrdp,mach0,comp
!
      implicit none
!
!------------------------
! --- [dummy arguments]
!------------------------
!
      integer,intent(inout) :: icall
      integer,intent(in)    :: iter,no
      real*8 ,intent(in)    :: deltt,time
      real*8 ,intent(in)    :: u_cell,v_cell,w_cell
      real*8 ,intent(in)    :: x_unit,y_unit,z_unit
      real*8 ,intent(in)    :: area
      real*8 ,intent(in)    :: x_cent,y_cent,z_cent
      real*8 ,intent(in)    :: dis_cell
      real*8 ,intent(in)    :: t_cell,t_wall,rho_cell
      real*8 ,intent(in)    :: rmu_cell,rmut_cell
      real*8 ,intent(inout) :: roughness

!----------------------
! --- [local entities]
!----------------------
      real(8) :: dum1,dum2=0.5d0
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL
      real*8 ,parameter	:: r1pn=1.d0/3.d0
      integer :: ICFL,ICV1,ICV2,II,I,ICFE,ICFS,ICVLA,ICVLB
!---------------------------------
! --- Explain using of subroutine
!---------------------------------
!     Attention : This subroutine should be called ONLY 
!                 by setting vel='sata-rough'
! 
!     icall :       Flag of calling this user subroutine by user
!                   =1 : this subroutine called Only one time 
!                   =0 : this subroutine called every time step
!     no    : BC no. defined in fflow.ctl
!     iter  : time step number 
!     deltt : time interval (s)
!     time  : time (s)
!     u_cell: Velocity of u,v,w at first cell center (m/s)
!     v_cell:
!     w_cell:
!     x_unit: unit vector of wall-face (outward)
!     y_unit: 
!     z_unit:
!     area  : area (m^2)
!     x_cent: coordinate of wall-face
!     y_cent:
!     z_cent: 
!     dis_cell : Minimum distance to wall from wall
!     t_cell: temp. of first cell (K)
!     t_wall: wall temp. of first cell (K)
!     rho_cell: density  of first cell (kg/m^3)
!     rmu_cell,rmut_cell (Pa.s)
!     roughness: roughness (m)
!---------------------------------
! --- 
!---------------------------------
      if(no==7) then  !
        icall=1
        roughness=1.3D-5*sin(x_cent*3.1415/0.5d0)
      endif
      if(no==8) icall=1

      return
      end subroutine user_rough_wall
