!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_src_uvw(deltt,iter,time,ICV,IMAT_U,iphs,
     &                      u,v,w,t,dens,
     &                      suu,suv,suw,utauu)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     iter        : time step
!     deltt       : dt
!     time        : time
!
!     ICV         : Index of Control Volume
!     IMATU       : Material No. of ICV
!     iphs        : No. of Euler two phase (=1 or 2)
!     u,v,w       : Velocities of ICV   [m/s]
!     dens        : Density of ICV    [kg/m^3]
!     utauu       : Friction velocity of wall surface nearest to ICV
!                   [m/s]
!
!     suu,suv,suw : Transport equations sorce terms of u,v,w [N/m^3]
!     
!     Example     : suu=1 for channel flow
!
!---------------------------------------------------------------------
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter
      integer,intent(in)    :: ICV,IMAT_U,iphs
      real*8 ,intent(in)    :: utauu
      real*8 ,intent(in)    :: u,v,w,t,dens
      real*8 ,intent(out)   :: suu,suv,suw
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! --- [local entities]
      suu=1.d0
      suv=0.d0
      suw=0.d0
      
      return
!
      end subroutine user_src_uvw
