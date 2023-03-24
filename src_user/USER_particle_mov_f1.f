!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine USER_particle_mov_f1
     &    (INI,IP,DT_P,U_slip,U_blow,Re_P,RE_B,TAU_P_r,P_mass,U_G,
     &     U_P,rho_G,rho_P,F1)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     User subroutine for particle drag coefficience F1
!
!=====================================================================
!
      IMPLICIT NONE
!      
      integer,intent(in)    :: INI,IP
      real*8 ,intent(in)    :: DT_P
      real*8 ,intent(in)    :: U_slip,U_blow,Re_P,RE_B,TAU_P_r,P_mass
      real*8 ,intent(in)    :: U_G(3),U_P(3),rho_G,rho_P
      real*8 ,intent(out)   :: F1
!
! --- [local entities]
!
      REAL*8  :: dum1,dum2


!
      
!
      return
      end subroutine
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
