CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE USER_particle_HeatExch(time,DT_IN,INI,IPS,IPE,
     &  INSIDE,JOUT,PARCEL,TMPP,
     &  PMASS,PYS,
     &  DYDTP,DHDTP,EVAP_Y,EVAP_H,CVVOLM,MASSYS,
     &  yys,prs,tmp,rho,rhop,rmu,vel,uvwp,rds,DIA,rmd,
     &  ierror)
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension 
      USE module_hpcutil     
      USE module_constant,only : SML
      use module_io,only       : ifle,ifll
      use module_particle,only : latent
      use module_particle,only : PP_sat,TT_vap,TT_boil,PP_0,TT_0
      use module_particle,only : EVAP_COMP
      use module_species,only  : gascns,r_wm,wm
      USE module_particle,only : Evap_CHO,Avap_CHO,ical_la_coal
      use module_particle,only : IFLAG_COAL,VS_CHO,V_CHO,V_H2O,VS_H2O,
     &                           VS_HCN,V_HCN
      use module_particle,only : cp_pa
      use module_species,only  : gascns
!
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)    :: INI,IPS,IPE
      REAL*8,INTENT(IN)     :: time,DT_IN
      INTEGER,INTENT(IN)    :: INSIDE   (      MXPRT)
      INTEGER,INTENT(INout) :: JOUT     (      MXPRT)
      REAL*8,INTENT(IN)     :: PARCEL   (      MXPRT)
      REAL*8,INTENT(IN)     :: TMPP     (      MXPRT)
      REAL*8,INTENT(IN)     :: CVVOLM   (      MXALLCV)
      REAL*8,INTENT(INOUT)  :: PYS      (      MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: PMASS    (      MXPRT)
!
      REAL*8,INTENT(INOUT)  :: DYDTP    (MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: DHDTP    (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: EVAP_Y   (      MXALLCV_P,MXCOMP)
      REAL*8,INTENT(INOUT)  :: EVAP_H   (      MXALLCV_P,2)
      REAL*8,INTENT(INOUT)  :: MASSYS   (      MXALLCV,MXcomp)

      INTEGER,INTENT(INOUT) :: ierror
!
      REAL*8,INTENT(IN)     :: RHO      (      MXALLCV)
      REAL*8,INTENT(IN)     :: RMU      (      MXALLCV)
      REAL*8,INTENT(IN)     :: VEL      (      MXALLCV,3)
      REAL*8,INTENT(IN)     :: YYS      (      MXALLCV,MXCOMP)
      REAL*8,INTENT(IN)     :: TMP      (      MXALLCV)
      REAL*8,INTENT(IN)     :: RHOP     (      MXPRT)
      REAL*8,INTENT(IN)     :: PRS      (      MXALLCV)
      REAL*8,INTENT(IN)     :: UVWP     (      MXPRT,3)
      REAL*8,INTENT(IN)     :: RDS      (      MXALLCV,MXcomp)
      REAL*8,INTENT(IN)     :: RMD      (      MXALLCV)
      REAL*8, intent(in)    :: DIA      (MXPRT)
!
! --- [local entities]
!
      integer :: IP,ICFL
      INTEGER :: IPCV,ICOM
      REAL*8  :: dum1,dum2,dum3,pdum,DT0,Kv,vol,dum4,pi
      REAL*8,parameter :: f_latent_pt=1.d0
            
      ierror=0
      pi=4.d0*datan(1.d0)
!
      return
      end subroutine USER_particle_HeatExch
