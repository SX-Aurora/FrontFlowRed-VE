!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_src_e2p
     &   (BC_no,BC_name,deltt,iter,time,IMAT_U,
     &   u,v,w,u2,v2,w2,p,dens,dens2,T_WALL,vol,
     &   T,T2,void1,void2,eva_mass)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
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
      integer,intent(in)      :: BC_no
      character(*),intent(in) :: BC_name
      real*8 ,intent(in)      :: deltt,time,vol
      integer,intent(in)      :: IMAT_U,iter
      real*8 ,intent(in)      :: u,v,w,p,dens
      real*8 ,intent(in)      :: u2,v2,w2,dens2
      real*8 ,intent(in)      :: T,T2,T_WALL
      real*8 ,intent(in)      :: void1,void2
      real*8 ,intent(out)     :: eva_mass
!---------------------------------------------------------------------
! --- [local entities]
!---------------------------------------------------------------------      
      real*8 :: dum1,dum2,dum3
      real(8),parameter :: void=0.0001
!      integer :: 
!---------------------------------------------------------------------      
      if(BC_no==2) then
        dum1=void1*vol*dens
        eva_mass=void*dum1/deltt/vol
      endif
!--------------------------------------------------------------------
      end subroutine user_src_e2p
!
