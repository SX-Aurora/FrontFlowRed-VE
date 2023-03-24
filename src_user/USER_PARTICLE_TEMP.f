!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE USER_PARTICLE_TEMP(MXPRT,IP,N_injctr,DELT_M,
     &           DP_M,TP_M,RHOP_M,MASSP_M,VELP_M,
     &           TMPP,DIAP,
     &           G_T_M,G_VEL_M,rmu_G_M,G_RHO_M,G_RAMD_M,G_CP_M
     &      )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      IMPLICIT REAL*8 (A-H,O-Z)
!!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)    :: MXPRT,IP,N_injctr
      REAL*8,INTENT(IN)     :: DELT_M
      REAL*8,INTENT(INOUT)  :: DP_M,TP_M,RHOP_M,MASSP_M,VELP_M(3)
      REAL*8,INTENT(IN)     :: TMPP(MXPRT),DIAP(MXPRT)
!
      REAL*8,INTENT(IN)     :: G_T_M,G_VEL_M(3),rmu_G_M,G_RHO_M
      REAL*8,INTENT(IN)     :: G_RAMD_M,G_CP_M
!      
!
! --- [local entities]
!
      PARAMETER (ND=10)
      PARAMETER (N=10)
      REAL*8,save,allocatable  ::  T(:,:),XP(:,:)
      DIMENSION XPNEW(ND),TN(ND),XPP(N),TT(N)
      DIMENSION RHO(ND),CP(ND),RAMDA(ND)
      DIMENSION VELP(3),VELGAS(3)
      
      INTEGER,save  ::  iflagini=0
      INTEGER :: NSUB=1,NNSUB
!------------------------------------------------------------
!
      return
!
      end SUBROUTINE USER_PARTICLE_TEMP
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
