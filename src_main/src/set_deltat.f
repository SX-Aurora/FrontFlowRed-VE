!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine set_deltat(mph,iter,
     &  LVEDGE,LBC_SSF,LCYCSF,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  CVVOLM,CVCENT,SFAREA,
     &  vscc,sgm,wifsld,LCYCOLD,vctr,ccc,wiface,
     &  rva,rho,rmx,rmut,deltt)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       vscc(:) - temp memory
!       sgm(:)  - temp memory !zhang-cvd
!       rmx(:)  - diffusivity (max(mu,rmd/cp))
! --- [module arguments]
!
      use module_dimension
      use module_material,only : ivscvv,ivscty,ivscke,
     &                           relaxC,relaxC_step,relaxci
      use module_hpcutil,only  : my_rank,NPE
      use module_deltat, only  : ideltt,const,dtmax,dtsafe
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_metrix,only   : tmpfac=>d2vect
      use module_Euler2ph,only : ieul2ph
      use module_deltat  ,only : maxcou,maxcou2,icvmxc
      use module_time    ,only : i_steady
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mph,iter
      integer,intent(in)  :: LVEDGE(2, MXCVFAC)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      integer,intent(in)  :: LCYCSF(   MXSSFBC)
      INTEGER,INTENT(IN)  :: MAT_NO(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_CV(   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_DCIDX(0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (0:MXMAT)
      real*8 ,intent(in)  :: CVVOLM(   MXALLCV)
      real*8 ,intent(in)  :: CVCENT(3, MXALLCV)
      real*8 ,intent(in)  :: SFAREA(4, MXCVFAC)
      REAL*8 ,INTENT(IN)  :: WIFACE(   MXCVFAC)
      real*8 ,intent(in)  :: rva   (   MXCVFAC,2)
      real*8 ,intent(in)  :: rho   (   MXALLCV)
      real*8 ,intent(in)  :: ccc   (   MXALLCV)
      real*8 ,intent(in)  :: rmx   (   MXALLCV)
      real*8 ,intent(in)  :: rmut  (   MXALLCV)
      real*8 ,intent(inout)  :: vscc(  MXALLCV)
      real*8 ,intent(inout)  :: sgm (  MXALLCV)
      real*8 ,intent(out)   :: deltt
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!
!
! --- [local entities]
!
      real*8,parameter :: vdfct=2.d0
      real*8  :: vscf,dx,dy,dz,vdt,dum1
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,IE
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE,
     &           IMATX
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP,myid
      real*8 :: delt_t,wi1,wi2
!
!
!-< 1. Constant time step >-
      do IIMAT=1,NMAT
      relaxC(IIMAT)=min(1.d0,(1.d0-relaxCi(IIMAT))*dble(iter)
     &             /dble(relaxC_step(IIMAT))+relaxCi(IIMAT))
      enddo
!

      deltt=dtmax
      delt_t=dtmax
      if(ideltt.eq.const.and.(i_steady==1.or.i_steady==2)) then
        maxcou=0.d0
        maxcou2=0.d0
        return
      endif
      if(i_steady==3) then
        deltt=dtmax
        return
      endif
      if(iter==1) then
        return
      endif
      if(mph==2) then
        delt_t=deltt
!        call FFRABORT(1,'ERR: FFR_2P NOT support [auto] option')
      endif
!---------------------
! --- VECTOR Computor
!---------------------
      if(ical_vect) then !NNOTVECT
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V  !index_c(myid)+1,index_c(myid+1)
        vscc(ICVL)=(rmx(ICVL)+vdfct*rmut(ICVL))
        sgm(ICVL)=0.d0
        enddo
! --- 
        call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &                 mat_cal,vscc)
! --- 
!CIDR NODEP
        do ICFL=ICFS_V,ICFE_V!index_f(myid)+1,index_f(myid+1)
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        tmpfac(ICFL,1)=max(vscc(ICVLA),vscc(ICVLB))*SFAREA(4,ICFL)
     &    *((CVCENT(1,ICVLB)-CVCENT(1,ICVLA))*SFAREA(1,ICFL)
     &     +(CVCENT(2,ICVLB)-CVCENT(2,ICVLA))*SFAREA(2,ICFL)
     &     +(CVCENT(3,ICVLB)-CVCENT(3,ICVLA))*SFAREA(3,ICFL))
     &    /((CVCENT(1,ICVLB)-CVCENT(1,ICVLA))**2
     &     +(CVCENT(2,ICVLB)-CVCENT(2,ICVLA))**2
     &     +(CVCENT(3,ICVLB)-CVCENT(3,ICVLA))**2)
        enddo
!CIDR NODEP
        do ICVL=ICVS_V,ICVE_V!index_c(myid)+1,index_c(myid+1)
          do IE=1,vctr(ICVL,0)
          sgm(ICVL)=sgm(ICVL)
     &      +max(0.d0,dble(vctr(ICVL,IE)/ABS(vctr(ICVL,IE)))
     &      *rva(abs(vctr(ICVL,IE)),1))
     &      +abs(tmpfac(abs(vctr(ICVL,IE)),1))
          enddo
        enddo
!!$omp end parallel do
        vdt=dtsafe/dtmax
        vdt=-1.d0
!CIDR NODEP
        do ICVL=ICVS_V,ICVE_V
        vdt=max(vdt,sgm(ICVL)/(rho(ICVL)*CVVOLM(ICVL)))
        enddo
!
        deltt=min(dtsafe/vdt,delt_t) !deltt=dtsafe/vdt
        maxcou=deltt*vdt
        if(NPE.gt.1) call hpcrmin(deltt)
        return
      endif
!
!
!-< 2. Auto time step >-
!
! --- < 2.1 set diffusivity >--
!
      do 400 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 400
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
!
      IMATX=IMAT
      if(IMAT<0)  IMATX=-IMAT
      if(ivscvv(IMATX).gt.0.or.ivscty(IMATX).
     &    gt.0.or.ivscke(IMATX).gt.0) then
        n=1
      else
        n=0
      endif
      dum1=dble(n)
      if(IMAT.gt.0) then
        DO 430 ICVL=ICVS,ICVE
!  1)
!        vscc(ICVL)=vdfct*dum1*(rmx(ICVL)+rmut(ICVL)) !???suginaka
!  2)
        vscc(ICVL)=dum1*(rmx(ICVL)+vdfct*rmut(ICVL))  ! default
!  3)
!        vscc(ICVL)=vdfct*dum1*(rmut(ICVL))             !zhang-cvd

 430    continue
      else
        DO 440 ICVL=ICVS,ICVE
        vscc(ICVL)=rmx(ICVL)
 440    continue
      endif
 400  continue
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,vscc)
!
!--< 2.2 estimate CFL condition >--
!
      sgm(:)=0.d0
      if(ieul2ph>0) then
        do 210 IIMAT=1,NMAT 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
        dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
        dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
        vscf=max(vscc(ICVLA),vscc(ICVLB))*SFAREA(4,ICFL)
     &      *(dx*SFAREA(1,ICFL)+dy*SFAREA(2,ICFL)+dz*SFAREA(3,ICFL))
     &      /(dx*dx+dy*dy+dz*dz)
        dum1=rva(ICFL,1)!*rho(ICVLA)
        sgm(ICVLA)=sgm(ICVLA)-min(0.d0,dum1)+abs(vscf)
        sgm(ICVLB)=sgm(ICVLB)+max(0.d0,dum1)+abs(vscf)
        enddo
  210   enddo
      else
        do IIMAT=1,NMAT 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
        dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
        dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
        vscf=max(vscc(ICVLA),vscc(ICVLB))*SFAREA(4,ICFL)
     &      *(dx*SFAREA(1,ICFL)+dy*SFAREA(2,ICFL)+dz*SFAREA(3,ICFL))
     &      /(dx*dx+dy*dy+dz*dz)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dum1=wi1*ccc(ICVLA)*rho(ICVLA)**2+wi2*ccc(ICVLB)*rho(ICVLB)**2
!        dum1=rva(ICFL,1)
        dum1=rva(ICFL,1)+sqrt(dum1)
        sgm(ICVLA)=sgm(ICVLA)-min(0.d0,dum1)+abs(vscf)
        sgm(ICVLB)=sgm(ICVLB)+max(0.d0,dum1)+abs(vscf)
        enddo
        enddo
      endif
!
      vdt=dtsafe/dtmax
      vdt=-1.d0
      do 220 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 220 
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      vdt=max(vdt,sgm(ICVL)/(rho(ICVL)*CVVOLM(ICVL)))
      enddo
  220 continue
!
      if(mph==2) then
         deltt=min(dtsafe/vdt,delt_t)
         maxcou2=deltt*vdt
      else
        deltt=min(dtsafe/vdt,delt_t) !dtsafe/vdt
        maxcou=deltt*vdt
      endif
!
      if(NPE.gt.1) then
        call hpcrmin(deltt)
      endif
!
      end subroutine set_deltat
!
