!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_src_fire
     &        (no,BC_name,deltt,iter,time,IMAT_U,iphs,ncomp,
     &         u,v,w,p,dens,comp,T_WALL,
     &         T,eva_T,eva_comp,eva_mass,eva_heat)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     no         : Boundary number by user
!     BC_name    : Boundary name
!     deltt      : dt
!     iter       : time step
!     IMAT_U     : Material No. of BC_name
!     iphs       : No. of Euler two phase (=1 or 2)
!     ncomp      : Total Species number
!     u,v,w      : Velocities of BC_name   [m/s]
!     p          : pressure  [N/m^3]
!     dens       : Density at BC_name [kg/m^3]
!     comp       : Mass friction at BC_name  [100%]
!     T_WALL     : wall Temperature at BC_name [K]
!     T          : Temperature at BC_name [K]
!     eva_comp   : mass fraction of evaporation [100%]
!     eva_mass   : evaporation mass rate [kg/m^2/s]
!     eva_T      : evaporation gas Temperature [K]
!     eva_heat   : latent heat [W/m^3=J/s/m^3] 
!                : [+] heat release
!                : [-] heat absorption
!--------------------------------------------------------------------- 
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)      :: deltt,time
      integer,intent(in)      :: no
      integer,intent(in)      :: IMAT_U,iphs,ncomp,iter
      character(*),intent(in) :: BC_name
      real*8 ,intent(in)      :: u,v,w,p,dens,comp(ncomp),T_WALL
      real*8 ,intent(in)      :: T
      real*8 ,intent(out)     :: eva_comp(ncomp),eva_mass
      real*8 ,intent(out)     :: eva_T,eva_heat
!-------------------------------------------------------------------- 
!--------------------------------------------------------------------
! --- [local entities]
!
      real*8 :: dum=0
      integer:: idum
      real(8),parameter :: pi=3.14
      real*8,parameter :: totms=8.16d0/408000.d0/0.81d0
!
!      
!--------------------------------------------------------------------
! --- ***(CAL_6)*** square 900[mm] with gas ejector
!      if (time <= 30.d0) then
!        eva_mass = totms*40.d0*time
!      else if (time <= 60.d0) then
!        eva_mass = totms*(40.d0/3.0d0*time + 800.0d0)
!      else if (time <= 160.d0) then
!        eva_mass = totms*(6.d0*time + 1240.d0)
!      else if (time <= 190.d0) then
!        eva_mass = totms*(2200.d0)
!      else if (time <= 220.d0) then
!        eva_mass = totms*(-40.d0/3.0d0*time + 4733.d0)
!      else if (time <= 250.d0) then
!        eva_mass = totms*(-160.d0/3.0d0*time + 13533.d0)
!      else
!        eva_mass = totms*(-15.d0*time + 3950.d0)
!      end if
! hex_fire:
!      if(no==5) then ! !4.5MW
!        eva_comp(:)=0.d0
!        eva_comp(1)=1.d0 
!        eva_mass=0.08942d0/6.4d0 ![kg/m^2/s]
!        eva_T=400.d0  !K
!        eva_heat=0.d0
!      endif
!      return
!
      
      if(no==8) then !smoke-1 !4.5MW
        eva_comp(:)=0.d0
        eva_comp(1)=1.d0
        eva_mass=0.08942d0/6.4d0 ![kg/m^2/s]
        eva_T=400.d0  !K
        eva_heat=0.d0
      endif
!      return

      if(no==9)then !smoke-9 !4.5MW
        eva_comp(:)=0.d0
        eva_comp(1)=1.d0
        eva_mass=0.08942d0/6.4d0 ![kg/m^2/s]
        eva_T=400.d0  !K
        eva_heat=0.d0
      elseif(no==10)then !smoke-3 !2.5MW
        eva_comp(:)=0.d0
        eva_comp(1)=1.d0
        eva_mass=0.04967d0/6.4d0 ![kg/m^2/s]
        eva_T=400.d0  !K
        eva_heat=0.d0
      elseif(no==11)then !smoke-4 !4.5MW
        eva_comp(:)=0.d0
        eva_comp(1)=1.d0
        eva_mass=0.08942d0/6.4d0 ![kg/m^2/s]
        eva_T=400.d0  !K
        eva_heat=0.d0 
      elseif(no==12)then !smoke-5 !2.5MW
        eva_comp(:)=0.d0
        eva_comp(1)=1.d0
        eva_mass=0.04967d0/6.4d0 ![kg/m^2/s]
        eva_T=400.d0  !K
        eva_heat=0.d0
      endif
      return
!
        eva_mass = totms*(2200.d0)
        eva_comp(2:ncomp)=0.0d0
        eva_comp(1)=1.0              !  [100%]
        eva_T=400.d0
! --- ****************      
!
!!      if(T_wall.gt.350.d0) then
!        eva_mass=0.2d0               !  [kg/m^2/s]
!        eva_comp(1)=1.0              !  [100%]
!        eva_T=400.d0
!!      endif
!
!
      return
!     
      end subroutine user_src_fire
