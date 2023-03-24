!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine reset_vars(imode,LVEDGE,LCYCSF,
     &           MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,mat_cal,MAT_DCIDX,
     &           pp0,aks,
     &           rva,yys,vel,hhh,tmp,prs,
     &           rva2,yys2,vel2,hhh2,tmp2,prs2
     &           )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_dimnsn,only   : xredc,yredc,zredc
      use module_model,   only : idrdp,incomp,mach0,
     &                           icaltb,sles,dles,lles
      use module_scalar,only   : rns_scl
      use module_Euler2ph,only : ieul2ph
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: imode
      integer,intent(in)   :: LVEDGE   (2, MXCVFAC)
      integer,intent(in)   :: LCYCSF   (   MXSSFBC)
      INTEGER,INTENT(IN)   :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)   :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)   :: MAT_NO   (   0:MXMAT)
      integer,intent(in)   :: MAT_CFIDX(   0:MXMAT)
      logical,INTENT(IN)   :: mat_cal   (  0:MXMAT)
      INTEGER,INTENT(IN)   :: MAT_DCIDX(     0:MXMAT)
      real*8,intent(inout) :: prs    (     MXALLCV,2)
      real*8,intent(inout) :: pp0    (     MXALLCV,2)
      real*8,intent(inout) :: aks(         MXALLCVR,MXrans,2)
!
      real*8,intent(inout) :: yys(         MXALLCV,mxcomp,2)
      real*8,intent(inout) :: vel    (     MXALLCV,3,2)
      real*8,intent(inout) :: rva    (     MXCVFAC,2)
      REAL*8,INTENT(INOUT) :: hhh  (       MXALLCV,2)
      REAL*8,INTENT(INOUT) :: tmp  (       MXALLCV,2)
      REAL*8,INTENT(INOUT) :: VEL2 (       MXALLCV2,3,2)
      REAL*8,INTENT(INOUT) :: YYS2 (       MXALLCV2,MXCOMP,2)
      REAL*8,INTENT(INOUT) :: RVA2 (       MXCVFAC2,2)   
      REAL*8,INTENT(INOUT) :: hhh2 (       MXALLCV2,2)
      REAL*8,INTENT(INOUT) :: tmp2 (       MXALLCV2,2)
      real*8,intent(inout) :: prs2   (     MXALLCV2,2)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,ICOM,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
!
!
      if(imode==1) goto 1000
      if(imode==2) goto 2000
      if(ical_vect) then
        if( xredc ) vel(:,1,:)=0.d0
        if( yredc ) vel(:,2,:)=0.d0
        if( zredc ) vel(:,3,:)=0.d0
        return
      endif
!------------------------------------------
!-< 1. Reset variables in solid part >-
!------------------------------------------
!--< 1.1 mass flux located at cell face >--
!------------------------------------------
      do 100 IIMAT=1,NMAT  !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) goto 100
      IMAT=MAT_NO(IIMAT)
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(IMAT.lt.0 ) then
        rva(ICFS:ICFE,1)=0.d0
        rva(ICFS:ICFE,2)=0.d0
      endif
  100 continue
!
      if(ieul2ph>0) then
        do 110 IIMAT=1,NMAT  !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 110
        IMAT=MAT_NO(IIMAT)
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        if(IMAT.lt.0) then
          rva2(ICFS:ICFE,1)=0.d0
          rva2(ICFS:ICFE,2)=0.d0
        endif
  110   continue
      endif
!--------------------------------------------
!--< 1.2 variables located at cell center >--
!--------------------------------------------
      do 200 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT.lt.0) then
        vel(ICVS:ICVE,:,:)=0.d0
        yys(ICVS:ICVE,1,:)=0.d0
        do ICOM=2,ncomp
        yys(ICVS:ICVE,ICOM,:)=0.d0
        enddo
        prs(ICVS:ICVE,:)=0.d0
        pp0(ICVS:ICVE,:)=0.d0
        if(rns_scl) then
          aks(ICVS:ICVE,:,:)=0.d0
        endif
      endif
 200  continue
!
      if(ieul2ph>0) then
        do 210 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT.lt.0) then
          vel2(ICVS:ICVE,:,:)=0.d0
          yys2(ICVS:ICVE,1,:)=0.d0
          do ICOM=2,ncomp
          yys2(ICVS:ICVE,ICOM,:)=0.d0
          enddo
          prs2(ICVS:ICVE,:)=0.d0
        endif
 210    enddo
!
      endif
!
!-< 2. Reset velocity in reduced direction >-
!
      if( xredc ) vel(:,1,:)=0.d0
      if( yredc ) vel(:,2,:)=0.d0
      if( zredc ) vel(:,3,:)=0.d0
      if(ieul2ph>0) then
        if( xredc ) vel2(:,1,:)=0.d0
        if( yredc ) vel2(:,2,:)=0.d0
        if( zredc ) vel2(:,3,:)=0.d0
      endif
      return
!-----------------
! --- imode==1
!-----------------
 1000 continue
      do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)

        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(IMAT.lt.0) then
          hhh(ICVS:ICVE,:)=hhh2(ICVS:ICVE,:)
          hhh(IDCS:IDCE,:)=hhh2(IDCS:IDCE,:)
          tmp(ICVS:ICVE,:)=tmp2(ICVS:ICVE,:)
          tmp(IDCS:IDCE,:)=tmp2(IDCS:IDCE,:)
        endif
        if(ieul2ph==1) then
          prs2(ICVS:ICVE,:)=prs(ICVS:ICVE,:) !5555
        endif
      enddo

      return
!
 2000 continue
!-----------------
! --- imode==2
!-----------------
      do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          hhh2(ICVS:ICVE,:)=hhh(ICVS:ICVE,:)
          hhh2(IDCS:IDCE,:)=hhh(IDCS:IDCE,:)
          tmp2(ICVS:ICVE,:)=tmp(ICVS:ICVE,:)
          tmp2(IDCS:IDCE,:)=tmp(IDCS:IDCE,:)
        endif
      enddo

      return
      end subroutine reset_vars
