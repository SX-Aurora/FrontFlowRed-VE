!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine USER_ini_particle
     &      (IP,INI,iter,time,dT,ncomp,
     &       DIM_P,RHO_P,XYZ_P,Vel_P,Ys_P,T_U)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     User subroutine for particle initial condition
!     (For Lagrangian two phase flow)
!
!     IP     : particle number
!     INI    : Injector number
!     iter   : time step
!     time   : time
!     dTtime : Dt
!     ncomp  : Species Number
!     DIM_P  : Particle Diameter
!     RHO_P  : Particle density
!     XYZ_P  : Particle coordinate
!     Vel_P  : Particle velcity
!     Ys_P   : Particle componate
!     T_U    : Particle tempture
!
!
!=====================================================================
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IP,INI,iter,ncomp
      real*8 ,intent(in)    :: time,dT
      real*8 ,intent(inout) :: DIM_P,RHO_P,XYZ_P(3),Vel_P(3)
      real*8 ,intent(inout) :: Ys_P(ncomp),T_U
!=====================================================================
! --- [local entities] 
!     
!      if(iter<105) then   !CHO-TEST4
!        T_U=800.d0
!      endif
!      return
!
      if(iter<500) then 
        T_U=1000.d0
      endif
!
      return
      end subroutine USER_ini_particle 
!
