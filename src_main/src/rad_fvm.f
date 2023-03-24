!
!      subroutine radfvm_admin()
!      subroutine getproperty()
!      subroutine calbounemi()
!      subroutine radesrterm()
!      subroutine makeradequ()
!      subroutine calradqfvm()
!      subroutine userdefradpro()
!      subroutine getgabsorb()
!      subroutine getptlepro()
!      subroutine getdivangle()
!      subroutine getdivweight()
!      subroutine scaphasefun()
!      subroutine wallreflectemi()
!      subroutine wallincienergy()
!      subroutine fvmscavalue()
!      subroutine wallanispower()
!      subroutine reflectdirect()
!-----------------------------------------------------------------------
!      for radiation heat transfer by FVM/DOM methods 
!-----------------------------------------------------------------------
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine radfvm_admin
     &  (ismpl,iphs,deltt,iter,p_grdc,time,ictl,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rva,rho,vel,prs,pp0,
     &  kdbt,vctr,cps,cp,S_inten,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LCYCOLD,wifsld,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,PAbsorb,PScatter,GWAbsorb,
     &  radinten,
     &  tmpbnd,yysbnd,htcbnd,mtcbnd,radbnd,tmp,yys,
     &  iterR,repsR,aepsR,errR,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only      : ifle,ifll
      use module_model   ,only: idrdp,mach0,comp,ical_t,ical_reac
      use module_flags   ,only: intgty,eulere
      use module_material,only: icnvty,lmtrty,ivscty,cal,
     &                          nflud,nofld,lclsd,nsold,nosld,
     &                          relaxh,rg_mark
      use module_param   ,only: yslw
      USE module_usersub, ONLY: src_t,src_r,usrno,usryes,src_fire
      use module_species, only: sw,acpk,c_p,enthalpy_ref,hform,spcnam,
     &                          wm,r_wm,Ri,spcno
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd,
     &                          phs_idx,phs_com,heat
      use module_time,    only: steady
      use module_metrix,only  : ys,rcomp,eva_comp,sinj_comp,erryys,ahk
      use module_cgsolver,only : aepsbcg,repsbcg,iterbcg
      use module_time,   only :  iters,itere
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IQ =>IW2K1
      use module_metrix,only   : bb
      use module_metrix,only   : aalw
      use module_material,only : ical_sld
      use module_model,only    : ical_vect
      use module_rad  ,only    : NGauss
!
! ---- Appended for Radiation Heat Transfer
!
!Ver.1 By Monte-Carlo | Zone Method,  
!In both the methods, heat exchange is calculated via exchange-area, which is
!appoximated in PREFFLOW. Ver.1 can solve radiations with constant properties.
!Developed in Apr.-Dec. 2005.
!
!Ver.2 By Finite-Volume-Method (or say Discrete Ordinate Method, Discrete Transfer Method)
!In principle, any kind of radiation problem could be solved.
!To be developed in 2006.
!
      use module_rad,    only  : gasmodel,ndiv
      use module_radsxf, only  : DivW,DivA
      use module_material,only : radmat,radfludflag
!
! 1. Update temperature & mass fraction
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      real*8 ,intent(inout) :: p_grdc
      integer,intent(in)    :: iter,iphs,ismpl
      integer,intent(inout) :: ictl
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: rva       (   MXCVFAC)
      real*8 ,intent(in)    :: rho       (   MXALLCV,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3  )
      real*8 ,intent(in)    :: prs       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(in)    :: yysbnd(       MXssfbc,MXcomp)
      real*8 ,intent(in)    :: htcbnd    (   MXssfbc)
      real*8 ,intent(in)    :: mtcbnd    (   MXssfbc)
      real*8 ,intent(in)    :: radbnd    (   MXssfbc)
      real*8 ,intent(in)    :: tmp       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: yys       (   MXALLCV,MXcomp,2)
      INTEGER,INTENT(INOUT) :: kdbt      (   MXCVFAC)
      integer,intent(out)   :: ierror,iterR(NDIV)
      real*8 ,intent(out)   :: repsR(NDIV),aepsR(NDIV)
      real*8 ,intent(out)   :: errR(NDIV)
      integer,intent(in)    :: LCYCOLD   (   MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld    (   MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      REAL*8 ,INTENT(INOUT)   :: cps     (        MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT)   :: cp      (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: S_inten (        MXALLCV)
      REAL*8 ,INTENT(INout)   :: P_rad      (MXALLCV_RAD,2)
      REAL*8 ,INTENT(INOUT)   :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: RadHeat    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: radinten  (NDNG,MXALLCV_RAD)
!
! --- [local entities]
!
      real*8  :: dum1,dum2,dum3
      integer :: kdv,kdt,kdy,kdk,kdp,ndf,nq,IMAX,no
      integer :: ICOM,IMD,ICH,IFLD,ICTP,ierr1=0
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: IMODE,IMAT_U,IBFS,IBFE,IBFL,IQMXCV
      integer :: NB,KD,IDC,ICV
      integer :: ICVA,ICVB,IBFP,ICFP,ICVP,IDCP
      integer :: iniusr,IDIV,iterq,ierr
      integer :: mark_not,ii0,idiv0,JSEG
      real*8  :: qvof,evapmdl,divqmax
      real*8  :: aepsq,repsq
      character*80 :: BC_name
!
      IQMXCV=MXALLCV
      allocate(IQ(IQMXCV,2),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in radfvm_admin'
        call FFRABORT(1,'radfvm_admin')
      endif
!

! --- cal. properties
!---------------------------------------------------------
!
      if(gasmodel(1:3)=='SLW') then
        CALL GetProSLW(deltt,iter,time,
     &		LVEDGE,LCYCSF,LBC_SSF,
     &		SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &		rho,vel,prs,pp0,yys,
     &		MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &          GWAbsorb,PAbsorb,PScatter,sumcoef,
     &		tmpbnd,tmp,ierror)
      elseif(gasmodel(1:4)=='WSGG') then
        CALL GetProWSGG(deltt,iter,time,
     &		LVEDGE,LCYCSF,LBC_SSF,
     &		SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &		rho,vel,prs,pp0,cps,cp,yys,
     &          MAT_NO,MAT_CV,MAT_INDEX,
     &		MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &          GWAbsorb,PAbsorb,PScatter,sumcoef,
     &		tmpbnd,tmp,ierror)
      else
        CALL GetProperty(deltt,iter,time,
     &	LVEDGE,LCYCSF,LBC_SSF,
     &	SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &	rho,vel,prs,pp0,cps,cp,yys,
     &  MAT_NO,MAT_CV,MAT_INDEX,
     &	MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,PAbsorb,
     &  PScatter,GWAbsorb,radinten,
     &	tmpbnd,tmp,ierror)
      ENDIF
!

!----------------------------------------------------------------
! --- cal. boundary emissions 
!----------------------------------------------------------------
      CALL CalBounEmi(ngauss,deltt,iter,time,
     &	LVEDGE,LCYCSF,LBC_SSF,
     &	SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &	rho,vel,prs,pp0,cps,cp,MAT_NO,MAT_CV,MAT_INDEX,
     &	MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,
     &  PAbsorb,PScatter,GWAbsorb,radinten,
     &	tmpbnd,tmp,ierror)
!
      DO 1100 IDIV=1,NDIV
      DO 1000 JSEG=1,ngauss !
      ii0=(JSEG-1)*NALLCV
      idiv0=(JSEG-1)*NDIV
!ii0,idiv0
!---------------------------------------------------------
! --- cal. emission and scattering source terms
!---------------------------------------------------------
      S_inten(:)=0.d0
      BB(:)=0.d0
      CALL RadESRTerm(JSEG,idiv,deltt,iter,time,
     &	LVEDGE,LCYCSF,LBC_SSF,
     &	SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &	MAT_NO,MAT_CV,MAT_INDEX,
     &	MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  BB,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,PAbsorb,
     &  PScatter,GWAbsorb,radinten,
     &	tmpbnd,tmp,ierror)
!
!---------------------------------------------------------
! --- make solve equation, AX=B 
!---------------------------------------------------------
      CALL MakeRadEqu(JSEG,idiv,deltt,iter,time,
     &	LVEDGE,LCYCSF,LBC_SSF,
     &	SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &	rho,vel,prs,pp0,cps,cp,kdbt,MAT_NO,MAT_CV,MAT_INDEX,
     &	MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  S_inten,BB,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,
     &  PAbsorb,PScatter,GWAbsorb,radinten,
     &	tmpbnd,tmp,ierror)
!---------------------------------------------------------
! --- solver AX=B by ICCG or Bi-CGSTAB 
!---------------------------------------------------------
!
      aepsq=aepsbcg
      repsq=repsbcg
      iterq=iterbcg
!
!
!
      if(ical_sld==0) then
	call utl_bcgstb('r',IDIV,
     &  MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &	NALLCV,NCV,NCVIN,MAXIE, IQ,IAL,aalw,
     &	MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &	bb,aepsq,repsq,iterq,ierr1,IMODE)
      else
	call utl_bcgstb_MAT1('r',
     &	MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &	NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &	MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &	bb,aepsq,repsq,iterq,ierr1,IMODE)
      endif
!
      iterR(IDIV)=iterq
      aepsR(IDIV)=aepsq
      repsR(IDIV)=repsq
!
      do ICV=1,NCV
	  radinten(idiv,ICV)=bb(ICV)
      enddo
      
 1000 enddo      !Iseg
 1100 CONTINUE   !idiv = 1,Ndiv
!---------------------------------------------------------
! --- cal. divq and qw
!---------------------------------------------------------
      CALL CalRadQFVM(ngauss,deltt,iter,time,
     &	LVEDGE,LCYCSF,LBC_SSF,
     &	SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &	rho,vel,prs,pp0,cps,cp,kdbt,MAT_NO,MAT_CV,MAT_INDEX,
     &	MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,
     &  PAbsorb,PScatter,GWAbsorb,radinten,
     &	tmpbnd,tmp,ierror)
!
! --- 
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,NALLCV,NCV,RadHeat)
	  CALL SOLVER_SEND_RECV(1,NALLCV,NCV,RadHeatFlux)
      ENDIF
!
      deallocate(IQ)
!
      return
!
 9999 continue
      if(my_rank.eq.ROOT) write(ifle,*) '(radfvm_admin)'
      ierror=1
!
      end subroutine radfvm_admin
!
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE GetProperty(deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rho,vel,prs,pp0,cps,cp,yys,
     &  MAT_NO,MAT_CV,MAT_INDEX,
     &  MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,
     &  PAbsorb,PScatter,GWAbsorb,radinten,
     &  tmpbnd,tmp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- [module arguments]
      use module_dimension
      use module_constant
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd
      use module_species,only  : spcnam,wm,Ri,spcno
      use module_rad
      use module_material,  only: radmat,radfludflag
      use module_boundary, only : radprop, radwalltype
      use module_radsxf,  only  : DivW,DivA
      use module_usersub,  only  : src_rad,usryes
      use module_metrix,only  : aayy => rcomp !=>aayy
!
      IMPLICIT NONE   !REAL*8(A-H,O-Z)
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: rho       (   MXALLCV,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3  )
      real*8 ,intent(in)    :: prs       (   MXALLCV,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV,2)
      REAL*8 ,INTENT(INOUT) :: cps       (   MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp        (   MXALLCV)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(inout) :: tmp       (   MXALLCV,2)
      real*8 ,intent(inout) :: yys       (   MXALLCV,MXcomp,2)
      REAL*8 ,INTENT(IN)       :: P_rad      (MXALLCV_RAD,2)
      REAL*8 ,INTENT(INOUT)   :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: RadHeat    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: radinten  (NDNG,MXALLCV_RAD)
      integer,intent(out)   :: ierror
!
!
!
      INTEGER :: iimat,imat,icvs,icve,icvl,icv,
     &           icomp,nb,kd,kdt,kdy,ibfs,ibfe,ibfl,icfl,idc
!
! --- [local entities]
!
!      REAL*8  :: aayy(ncomp)
!-----------------------------------------------------------------------
!for media
c	radmat(:,:)
c	!1: gabsorb, 2: pabsorb, 3: pscatter, 4:scabeta, 5: ptcdiamter
c	radfludtype(:,:)
c       !1: gastype, 2: scatype, 3: ptccal, 4: fgpflag
!gastype
!	!1: gray_gas,  2: real_gas
!scatype	
!   1,	Linear-anisotropic scatter / isotropic scatter, OK
!   2,  large-diffuse particle scatter, OK
!   3,	Rayleigh scatter, OK
!   4,	Mie scatter
!   5,  user defined scatter
! for wall 
! radprop(1,nb)	: rademi
! radprop(2,nb)	: radabsorb
! radwalltype(nb) : radtype
!-----------------------------------------------------------------------
      ierror=0
      GWAbsorb = 0.d0
      PAbsorb  = 0.d0
      PScatter = 0.d0
!
      IF(src_rad.eq.usryes)	 then
	do IIMAT=1,NMAT    !ICV=1,NCV
	if(.not.mat_cal(IIMAT)) cycle
	IMAT=MAT_NO(IIMAT)
	if(IMAT.LE.0) cycle		!not fluid
	ICVS=MAT_CVEXT(IIMAT-1)+1
	ICVE=MAT_CVEXT(IIMAT)
	do ICVL=ICVS,ICVE
	  ICV=MAT_CV(ICVL)
	  do Icomp=1,ncomp
            aayy(Icomp)=yys(ICVL,Icomp,1)
	  end do
	  call user_rad_pro(ICV,ncomp,spcno,CVCENT(:,ICVL),
     &	  tmp(ICVL,1),prs(ICVL,1),rho(ICVL,1),aayy,
     &	  GWAbsorb(ICVL),PAbsorb(ICVL),PScatter(ICVL))
	  do Icomp=1,ncomp
             yys(ICVL,Icomp,1)=aayy(Icomp)
	  end do
	end do
	end do
	GOTO 7777
      end if
!inner cv
      do IIMAT=1,NMAT    !ICV=1,NCV
	if(.not.mat_cal(IIMAT)) cycle
	IMAT=MAT_NO(IIMAT)
	if(IMAT.LE.0) cycle		!not fluid
	ICVS=MAT_CVEXT(IIMAT-1)+1
	ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
	GWAbsorb(ICVL) = RadMat(1,IIMAT)
	end do
	IF(radfludflag(3,IMAT).eq.0) THEN	!specify
          do ICVL=ICVS,ICVE
	  PAbsorb(ICVL) = RadMat(2,IIMAT)
	  PScatter(ICVL) = RadMat(3,IIMAT)
	  end do
	ELSE	!calculate
c	do ICVL=ICVS,ICVE
c	CALL GetPtcPro(PAbsorb(ICVL),PScatter(ICVL),
c     &	RadMat(5,IIMAT),tmp(ICVL,1),pp0(ICVL,1),yys(ICVL,1))
c	end do
        END IF
      end do
!boundary dummy cv
 7777   continue

      do nb=1,nbcnd
	IIMAT=MAT_BCIDX(nb,1)
	if(.not.mat_cal(IIMAT)) cycle
	kd=kdbcnd(0,nb)
        kdt=kdbcnd(2,nb)
	kdy=kdbcnd(3,nb)
	IBFS=LBC_INDEX(nb-1)+1
	IBFE=LBC_INDEX(nb)
	do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  GWAbsorb(IDC)=radprop(1,nb)
	enddo
      end do
8888	RETURN
!-------> end of GetProperty()
	END SUBROUTINE GetProperty
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE CalBounEmi(nseg,deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rho,vel,prs,pp0,cps,cp,MAT_NO,MAT_CV,MAT_INDEX,
     &  MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,
     &  PAbsorb,PScatter,GWAbsorb,radinten,
     &  tmpbnd,tmp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd,ktneum,kttrns,kdprdc,kdtchi,
     &                          kdsld,kdovst
      use module_rad    ,only:  ndiv
      use module_material,  only: radmat,radfludflag
      use module_radsxf,   only : frad,wrad
      use module_boundary, only : radprop, radwalltype
      use module_radsxf,  only  : DivW,DivA
      IMPLICIT NONE   !REAL*8(A-H,O-Z)
	
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: nseg,iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: rho       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3  )
      real*8 ,intent(in)    :: prs       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV  ,2)
      REAL*8 ,INTENT(INOUT) :: cps   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp    (       MXALLCV)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(inout) :: tmp       (   MXALLCV  ,2)
      REAL*8 ,INTENT(INOUT)   :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: RadHeat    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: radinten  (NDNG,MXALLCV_RAD)
      REAL*8 ,INTENT(IN)      :: P_rad     (MXALLCV_RAD,2)
      integer,intent(out)     :: ierror
!
! --- [local entities]
!
      real*8,PARAMETER :: SHIGMA=5.67d-8,
     &                    PAI=3.14159265d0,FUSHE=1.8048D-8
      integer :: ii0,idiv0,nb
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: kdv,kdt,kdy,kdk,kdp,ndf,nq,IMAX,no,iimat2
      integer :: imat1,imat2,kd,ibfs,ibfe,ibfl,icv,idc,icfp,icvp,idcp
!
      INTEGER ::  Iseg
      REAL*8  ::  WN(3)
!--------------------------------------------------------------------
!radwalltype:		1,diffuse-wall	,OK 
!			2,mirror-wall	,OK
!			3,directional-diffuse wall,OK
!			4,user defined
!--------------------------------------------------------------------
!boundary cv
      ierror=0
      DO 1200 Iseg=1,nseg
        ii0 = (Iseg-1)*NALLCV
	idiv0=(Iseg-1)*NDIV
        DO 1000 nb=1,nbcnd
	IIMAT=MAT_BCIDX(nb,1)
	IIMAT2=MAT_BCIDX(nb,2)
	IMAT1=MAT_NO(IIMAT)
	IMAT2=MAT_NO(IIMAT2)
	if(.not.mat_cal(IIMAT)) cycle
	kd=kdbcnd(0,nb)
	kdt=kdbcnd(2,nb)
	IBFS=LBC_INDEX(nb-1)+1
	IBFE=LBC_INDEX(nb)
	IF(kd.EQ.kdprdc) CYCLE
	IF(kd.EQ.kdilet.OR.kd.EQ.kdolet.OR.kd.EQ.kdtchi) THEN
!
! --- WRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),IDC)
!
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  WN(1) = SFAREA(1,ICFL)
	  WN(2) = SFAREA(2,ICFL)
	  WN(3) = SFAREA(3,ICFL)
	  FRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),ICV)
	  CALL WallReflectEmi(1,NDIV,
     &	    FRAD,WN,tmp(IDC,1),WRAD,
     &	    radprop(1,nb),sumcoef(ii0+IDC)) 
          radinten((idiv0+1):(idiv0+ndiv),IDC)=WRAD(1:NDIV) 
	  end do
	ELSE IF(kd.EQ.kdsymm) THEN
!--------------------------------------------
! --- treat as mirror surface (gray_degree=0)
!--------------------------------------------
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  WN(:) = SFAREA(1:3,ICFL)	!directing IDC
	  FRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),ICV)
!zhang??/	  CALL WallReflectEmi(2,NDIV,
	  CALL WallReflectEmi(2,NDIV,
     &	    FRAD,WN,0.d0,WRAD,
     &	    0.d0,sumcoef(ii0+IDC))
	  radinten((idiv0+1):(idiv0+ndiv),IDC)=WRAD(1:NDIV)
	  end do
	ELSE IF(kd.EQ.kdintr.or.kd.EQ.kdsld.or.kd==kdovst) THEN
	  if(min(IMAT1,IMAT2).GT.0.or.max(IMAT1,IMAT2).LT.0) then
            cycle
	  else
!            IF(IMAT1.GT.IMAT2) THEN		!imat1 is fluid
              do IBFL=IBFS,IBFE
	      ICFL=LBC_SSF(IBFL)
	      ICV=LVEDGE(1,ICFL)
	      IDC=LVEDGE(2,ICFL)
              WN(1) = SFAREA(1,ICFL) !directing IDC (outward)
	      WN(2) = SFAREA(2,ICFL)
	      WN(3) = SFAREA(3,ICFL)
	      FRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),ICV)
              CALL WallReflectEmi(radwalltype(nb),NDIV,
     &	      FRAD,WN,tmp(IDC,1),WRAD,
     &	      radprop(1,nb),sumcoef(ii0+IDC))
              radinten((idiv0+1):(idiv0+ndiv),IDC)=WRAD(1:NDIV)
	      end do
!	    ELSE
!	      do IBFL=IBFS,IBFE
!              ICFP=LCYCSF(IBFL)
!	      ICVP=LVEDGE(1,ICFP)
!	      IDCP=LVEDGE(2,ICFP)
!	      WN(1) = SFAREA(1,ICFP)
!	      WN(2) = SFAREA(2,ICFP)
!	      WN(3) = SFAREA(3,ICFP)
!	      FRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),ICVP)
!              CALL WallReflectEmi(radwalltype(nb),NDIV,
!     &	      FRAD,WN,tmp(IDCP,1),WRAD,
!     &	      radprop(1,nb),sumcoef(ii0+IDCP))
!              radinten((idiv0+1):(idiv0+ndiv),IDCP)=WRAD(1:NDIV)
!	      enddo
!	    ENDIF
          endif
        ELSE	!wall or other real wall-surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!revised by Jiang, 2006/11/21	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(kdt.eq.kttrns) then
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  WN(1) = SFAREA(1,ICFL)		!directing IDC
	  WN(2) = SFAREA(2,ICFL)
	  WN(3) = SFAREA(3,ICFL)
	  FRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),ICV)
     	  CALL WallReflectEmi(radwalltype(nb),NDIV,
     &       FRAD,WN,tmp(ICV,1),WRAD,
     &	     radprop(1,nb),sumcoef(ii0+IDC))
	  radinten((idiv0+1):(idiv0+ndiv),IDC)=WRAD(1:NDIV)
          end do
        elseif(kdt.NE.ktneum) then  !jiang, debud 2008/03/14
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  WN(1) = SFAREA(1,ICFL)	!directing IDC
	  WN(2) = SFAREA(2,ICFL)
	  WN(3) = SFAREA(3,ICFL)
	  FRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),ICV)
     	  CALL WallReflectEmi(radwalltype(nb),NDIV,
     &	                      FRAD,WN,tmp(IDC,1),WRAD,
     &	                      radprop(1,nb),sumcoef(ii0+IDC))
	  radinten((idiv0+1):(idiv0+ndiv),IDC)=WRAD(1:NDIV) 
          end do
        else
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  WN(1) = SFAREA(1,ICFL)	!directing IDC
	  WN(2) = SFAREA(2,ICFL)
	  WN(3) = SFAREA(3,ICFL)
	  FRAD(1:NDIV)=radinten((idiv0+1):(idiv0+ndiv),ICV)
     	  CALL WallReflectEmi(radwalltype(nb),NDIV,
     &	                      FRAD,WN,tmp(ICV,1),WRAD,
     &	                      radprop(1,nb),sumcoef(ii0+IDC))
	  radinten((idiv0+1):(idiv0+ndiv),IDC)=WRAD(1:NDIV) 
          end do
        endif
        endif
1000	CONTINUE
1200	CONTINUE


	RETURN
!-------> end of CalBounEmi()
!
	END SUBROUTINE CalBounEmi
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE RadESRTerm(Iseg,idiv,deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  MAT_NO,MAT_CV,MAT_INDEX,
     &  MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  BB,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,PAbsorb,
     &  PScatter,GWAbsorb,radinten,
     &  tmpbnd,tmp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd

!      use module_metrix,only   : bb
      use module_rad    ,only   : NDIV
      use module_material, only:  radmat,radfludflag
!      use module_radsxf,   only : GWAbsorb, PAbsorb, PScatter
!      use module_radsxf,   only : radinten,sumcoef!,P_rad
      use module_boundary, only : radprop, radwalltype
      use module_radsxf,  only  : DivW,DivA
      use module_metrix,only   : aalw						
!
      IMPLICIT NONE 
!
      real*8,PARAMETER :: SHIGMA=5.67d-8,
     &                    PAI=3.14159265d0,FUSHE=1.8048D-8
      real*8  :: cosv,beta,wfai
      integer :: ii0,idiv0,nb,jdiv
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: kdv,kdt,kdy,kdk,kdp,ndf,nq,IMAX,no,iimat2
      integer :: imat1,imat2,kd,ibfs,ibfe,ibfl,icv,idc,icfp,icvp,idcp
      
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: Iseg,idiv,iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(inout) :: tmp       (   MXALLCV  ,2)
      REAL*8 ,INTENT(out)   :: BB     (      MXCV)
      REAL*8 ,INTENT(INOUT)   :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: RadHeat    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: radinten  (NDNG,MXALLCV_RAD)
      REAL*8 ,INTENT(IN)      :: P_rad      (MXALLCV_RAD,2)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      real*8  :: dum1,dum2
!
      aalw(:,:)=0.d0
      ii0 = (iseg-1)*NALLCV
      idiv0=(iseg-1)*NDIV
      do 1000 IIMAT=1,NMAT    !ICV=1,NCV
	if(.not.mat_cal(IIMAT)) cycle
	IMAT=MAT_NO(IIMAT)
	ICVS=MAT_CVEXT(IIMAT-1)+1
	ICVE=MAT_INDEX(IIMAT)
	IF(IMAT.LT.0) CYCLE
	do ICVL=ICVS,ICVE
        dum1=tmp(ICVL,1)
        dum2=(GWAbsorb(ii0+ICVL)+PAbsorb(ICVL))*sumcoef(ii0+ICVL)
	BB(ICVL)=FUSHE*dum1**4*dum2
	enddo
!
        wfai=0.d0
        DO JDIV=1,NDIV
        cosv=DivA(1,IDIV)*DivA(1,JDIV)
     &      +DivA(2,IDIV)*DivA(2,JDIV)
     &	    +DivA(3,IDIV)*DivA(3,JDIV)
        CALL FVMScaValue(cosv,radfludflag(2,IMAT),beta,wfai)
        wfai=wfai*DivW(JDIV)
!
        do ICVL=ICVS,ICVE
        IF(PScatter(ICVL).GT.0.d0) THEN
          dum1=PScatter(ICVL)/(4.0*PAI)*wfai*sumcoef(ii0+ICVL)
          BB(ICVL)=BB(ICVL)+dum1*radinten(idiv0+JDIV,ICVL)
        else
          dum1=PScatter(ICVL)/(4.0*PAI)*wfai*sumcoef(ii0+ICVL)  
          BB(ICVL)=BB(ICVL)+dum1*radinten(idiv0+JDIV,ICVL)    !zhangrad???
        ENDIF
        enddo
        ENDDO  !zhangrad
!
        do ICVL=ICVS,ICVE
          dum1=P_rad(ICVL,1)
          dum2=P_rad(ICVL,2) 
          BB(ICVL)=BB(ICVL)*CVVOLM(ICVL)
     %  +dum1!-dum2*radinten(idiv0+IDIV,ICVL) !zhangrad
        enddo
 1000	ENDDO
!
	RETURN
!
!-------> end of RadESRTerm()
!
	END SUBROUTINE RadESRTerm
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE MakeRadEqu(ISEG,idiv,deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rho,vel,prs,pp0,cps,cp,kdbf,MAT_NO,MAT_CV,MAT_INDEX,
     &  MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  S_inten,BB,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,
     &  PAbsorb,PScatter,GWAbsorb,radinten,
     &  tmpbnd,tmp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IQ =>IW2K1
      use module_metrix,only   : aalw
      use module_metrix,only   : aax,ipx
      use module_material,only : ical_sld
      use module_model,only    : ical_vect,nthrds,ical_MHD
      use module_rad
      use module_material, only:  radmat,radfludflag
      use module_boundary, only : radprop,radwalltype
      use module_radsxf,  only  : DivW,DivA

      IMPLICIT NONE
!
      real*8,PARAMETER :: SHIGMA=5.67d-8,
     &                    PAI=3.14159265d0,FUSHE=1.8048D-8
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iseg,idiv,iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: rho       (   MXALLCV ,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3  )
      real*8 ,intent(in)    :: prs       (   MXALLCV ,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV  ,2)
      REAL*8 ,INTENT(INOUT) :: cps       (   MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp        (   MXALLCV)
      INTEGER,INTENT(INOUT) :: kdbf      (   MXCVFAC)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(in)    :: tmp       (   MXALLCV  ,2)
      real*8 ,intent(inout) :: S_inten   (   MXALLCV)
      integer,intent(out)   :: ierror
      REAL*8,INTENT(IN)     :: P_rad     (MXALLCV_RAD,2)
      REAL*8 ,INTENT(INOUT)   :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: RadHeat    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: radinten  (NDNG,MXALLCV_RAD)
      REAL*8 ,intent(inout)   :: BB    (MXCV)
!
!------private
!
      
      REAL*8 :: DV(3),AngleDiv(3),dum1,coefva,dum2,coefvb,cosva,iaflg
      integer :: ii0,idiv0,nb,icva,icvb,icx,k,m,ibflg
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: kdv,kdt,kdy,kdk,kdp,ndf,nq,IMAX,no,iimat2
      integer :: imat1,imat2,kd,ibfs,ibfe,ibfl,icv,idc,icfp,icvp,idcp
!
      
! --- Reset
!
      ii0=(Iseg-1)*NALLCV
      idiv0=(Iseg-1)*NDIV
!	
      IQ=0
      aalw(NCV+1:NALLCV,0)=0.d0
      DO IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      IMAT=MAT_NO(IIMAT)
      IF(IMAT.GT.0) THEN           !fluid 
        do ICVL=ICVS,ICVE
        dum2=P_rad(ICVL,2) 
        dum1=(GWAbsorb(ii0+ICVL)+PAbsorb(ICVL)
     &	    +PScatter(ICVL))*CVVOLM(ICVL)
	aalw(ICVL,0)=aalw(ICVL,0)+dum1+dum2   !zhangrad
        !BB(ICVL)=BB(ICVL)-dum1*radinten(idiv0+IDIV,ICVL)   !zhangrad
        enddo
      ELSE				! solid
	do ICVL=ICVS,ICVE
	aalw(ICVL,0)=1.0
	bb(ICVL)=0.d0
        end do	
      END IF	
      END DO
!
      CALL GetDivAngle(IDIV,NDIV,AngleDiv)
!
! --- calculation of intensity of radiation
!
      DO IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(IMAT.LT.0) CYCLE
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
!      DV(1)=CVCENT(1,ICVB)-CVCENT(1,ICVA)
!      DV(2)=CVCENT(2,ICVB)-CVCENT(2,ICVA)
!      DV(3)=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!      dum1=dsqrt(DV(1)**2+DV(2)**2+DV(3)**2)
      dum1=dsqrt(AngleDiv(1)**2+AngleDiv(2)**2+AngleDiv(3)**2)
      coefva=(AngleDiv(1)*SFAREA(1,ICFL)
     &	     +AngleDiv(2)*SFAREA(2,ICFL)
     &	     +AngleDiv(3)*SFAREA(3,ICFL))*SFAREA(4,ICFL)/dum1
      dum1=max(coefva,0.d0)*radinten(idiv0+IDIV,ICVA)
     &    +min(coefva,0.d0)*radinten(idiv0+IDIV,ICVB)
      S_inten(ICVA)=S_inten(ICVA)+dum1
      S_inten(ICVB)=S_inten(ICVB)-dum1
      enddo
      enddo

!
! --- 
!
      call bc_kdbt_rad(LBC_SSF,LCYCSF,mat_cal,kdbf)
!
      DO IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
	ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
	if(IMAT.LT.0) CYCLE
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).ge.99) cycle
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        DV(1)=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        DV(2)=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        DV(3)=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!--------------------------zhang-jiang
        dum1=dsqrt(AngleDiv(1)**2+AngleDiv(2)**2+AngleDiv(3)**2)
        coefva=(AngleDiv(1)*SFAREA(1,ICFL)
     &	       +AngleDiv(2)*SFAREA(2,ICFL)
     &	       +AngleDiv(3)*SFAREA(3,ICFL))*SFAREA(4,ICFL)/dum1
        dum1=max(coefva,0.d0)-min(coefva,0.d0)
        coefvb=-coefva
!	DV(1) = SFCENT(1,ICFL)-CVCENT(1,ICVA)
!	DV(2) = SFCENT(2,ICFL)-CVCENT(2,ICVA)
!	DV(3) = SFCENT(3,ICFL)-CVCENT(3,ICVA)
!	cosva =	 DV(1)*SFAREA(1,ICFL)
!     &		+DV(2)*SFAREA(2,ICFL)
!     &		+DV(3)*SFAREA(3,ICFL)
!	if(cosva.ge.0) then
! 	  iaflg = 1
!	else
!	  iaflg = -1
!	end if
!	ibflg=-iaflg
!	coefva = (AngleDiv(1)*SFAREA(1,ICFL)
!     &	   +AngleDiv(2)*SFAREA(2,ICFL)
!     &	   +AngleDiv(3)*SFAREA(3,ICFL))*iaflg*SFAREA(4,ICFL)
!	coefvb = -coefva
!--------------------------
!        if(kdbf(ICFL)==2) dum1=0.d0
        if(kdbf(ICFL).eq.0) then
          ICVA=msk(ICVA)
          ICVB=msk(ICVB)
          IQ(ICVA,1)=IQ(ICVA,1)+1
          IAL(ICVA,IQ(ICVA,1))=ICVB
          aalw(ICVA,0)=aalw(ICVA,0)+max(coefva,0.d0)
          !aalw(ICVA,0)=aalw(ICVA,0)-min(coefva,0.d0)
          aalw(ICVA,IQ(ICVA,1))=min(coefva,0.d0)
          !aalw(ICVA,IQ(ICVA,1))=-max(coefva,0.d0)

	  IQ(ICVB,1)=IQ(ICVB,1)+1
	  IAL(ICVB,IQ(ICVB,1))=ICVA

          aalw(ICVB,0)=aalw(ICVB,0)-min(coefva,0.d0)
          !aalw(ICVB,0)=aalw(ICVB,0)+max(coefva,0.d0)
          aalw(ICVB,IQ(ICVB,1))=-max(coefva,0.d0)
          !aalw(ICVB,IQ(ICVB,1))=min(coefva,0.d0)
          
        else	                            
          aalw(ICVA,0)=aalw(ICVA,0)+max(coefva,0.d0)
          ! aalw(ICVA,0)=aalw(ICVA,0)-min(coefva,0.d0)
!          aalw(ICVA,0)=aalw(ICVA,0)+dum1   
          bb(ICVA)=bb(ICVA)+max(-coefva,0.d0)      !zhangrad
     &		*radinten(idiv0+idiv,ICVB)
        endif
        enddo
      ENDDO
!
!-------------------------------------------------------
!--< 2.4 rearrange order to split lower & upper part >--
!-------------------------------------------------------
      do 240 IIMAT=1,NMAT   !ICV=1,NALLCV
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 255 ICVL=ICVS,ICVE
        icx=IQ(ICVL,1)
        k=0
        do 241 m=1,icx
        if(IAL(ICVL,m).lt.ICVL.and.IAL(ICVL,m)/=0) then
! --- lower part: IQ(ICVL,1)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  241   continue
        IQ(ICVL,1)=k
        do 242 m=1,icx
        if(IAL(ICVL,m).gt.ICVL.and.IAL(ICVL,m)/=0) then
! --- upper part: IQ(ICVL,2)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  242   continue
        IQ(ICVL,2)=k
        aalw(ICVL,1:)=0.d0
        IAL(ICVL,:)=0
        do 243 m=1,k
          aalw(ICVL,m)=aax(m)
          IAL(ICVL,m)=ipx(m)
 243    continue
 255    continue
 240    enddo
!
! --- 
!
        !BB(1:NCV)=BB(1:NCV)-S_inten(1:NCV) !zhangrad
!
        IF(NPE.GT.1) THEN
  	  CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,AALW(:,0))
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
        ENDIF
!
        return
	END SUBROUTINE MakeRadEqu
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE CalRadQFVM(nseg,deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rho,vel,prs,pp0,cps,cp,kdbf,MAT_NO,MAT_CV,MAT_INDEX,
     &  MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,
     &  PAbsorb,PScatter,GWAbsorb,radinten,
     &  tmpbnd,tmp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd,kttrns,kdtchi,kdsld,kdprdc
     &                          ,kdovst
      use module_rad     , only:  NDIV
      use module_material, only:  radmat,radfludflag
      use module_boundary, only : radprop, radwalltype
      use module_radsxf,  only  : DivW,DivA
      use module_radsxf,   only : frad,wrad
!
      IMPLICIT NONE    !REAL*8(A-H,O-Z)
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: nseg,iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: rho       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3  )
      real*8 ,intent(in)    :: prs       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV  ,2)
      REAL*8 ,INTENT(INOUT) :: cps   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp    (       MXALLCV)
      INTEGER,INTENT(INOUT) :: kdbf      (   MXCVFAC)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(inout) :: tmp       (   MXALLCV  ,2)
      REAL*8 ,INTENT(INOUT)   :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: RadHeat    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: radinten  (NDNG,MXALLCV_RAD)
      REAL*8,INTENT(IN)     :: P_rad     (MXALLCV_RAD,2)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      real*8,PARAMETER :: SHIGMA=5.67d-8,
     &                    PAI=3.14159265d0,FUSHE=1.8048D-8
      REAL*8	 ::	WN(3)
      REAL*8	 ::	twheat(nbcnd,3),ven,vinc
      INTEGER  :: idiv0,ii0,iseg,idiv
      integer :: nb,icva,icvb,icx,k,m,idca,idcb
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: kdv,kdt,kdy,kdk,kdp,ndf,nq,IMAX,no,iimat2
      integer :: imat1,imat2,kd,ibfs,ibfe,ibfl,icv,idc,icfp,icvp,idcp
      
!
!-----------------------------------------------------------------
!radwalltype:		1,  diffuse-wall	,OK 
!			2,  mirror-wall		,OK 
!			3,  directional-diffuse wall	,OK 
!			4,  user defined
!-----------------------------------------------------------------
!
      RadHeat=0.d0
      RadHeatFlux=0.d0
      do  IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      ven=4.0*shigma*tmp(ICVL,1)**4
      RadHeatFlux(ICVL)=0.d0
      do iseg=1,nseg   !nseg=1 & ngauss
      vinc=0.d0
      ii0=(Iseg-1)*NALLCV
      idiv0=(Iseg-1)*NDIV
      DO idiv=1,ndiv
      vinc=vinc+radinten(idiv0+IDIV,ICVL)*DivW(idiv)
      ENDDO
      vinc=(ven*sumcoef(ii0+ICVL)-vinc)
     &	*(GWAbsorb(ii0+ICVL)+PAbsorb(ICVL))
      RadHeatFlux(ICVL)=RadHeatFlux(ICVL)+vinc
      enddo
      RadHeat(ICVL)=RadHeatFlux(ICVL)*CVVOLM(ICVL)
      enddo
      enddo


	DO 1000 nb=1,nbcnd
	IIMAT=MAT_BCIDX(nb,1)
	IIMAT2=MAT_BCIDX(nb,2)
	IMAT1=MAT_NO(IIMAT)
	IMAT2=MAT_NO(IIMAT2)
	if(.not.mat_cal(IIMAT)) cycle
	kd=kdbcnd(0,nb)
        kdt=kdbcnd(2,nb)
	IBFS=LBC_INDEX(nb-1)+1
	IBFE=LBC_INDEX(nb)
	IF(kd.EQ.kdsymm) CYCLE
	IF(kd.EQ.kdilet.OR.kd.EQ.kdolet.OR.kd.EQ.kdtchi) THEN
!treat as black surface (diffuse surface)
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  WN(1:3) = SFAREA(1:3,ICFL)		!directing IDC
	  ven=0.d0
	  do Iseg=1,nseg
	  ii0=(Iseg-1)*NALLCV
	  idiv0=(Iseg-1)*NDIV
	  FRAD(1:NDIV) = radinten(idiv0+1:idiv0+ndiv,ICV)
          CALL WallInciEnergy(radwalltype(nb),NDIV,
     &	  FRAD,WN,vinc)
	  ven=ven+vinc
	  end do
	  RadHeatFlux(IDC)=shigma*(tmp(IDC,1)**4-ven)*radprop(2,nb)
	  RadHeat(IDC)=RadHeatFlux(IDC)*SFAREA(4,ICFL)
	  end do
       ELSE IF(kd.EQ.kdintr.or.kd.EQ.kdsld.or.kd==kdovst) THEN	!interface
	  if(min(IMAT1,IMAT2).LT.0) then	
            IF(IMAT1.GT.IMAT2) THEN
	      do IBFL=IBFS,IBFE
	      ICFL=LBC_SSF(IBFL)
	      ICV=LVEDGE(1,ICFL)
	      IDC=LVEDGE(2,ICFL)
              ICFP=LCYCSF(IBFL)
	      ICVP=LVEDGE(1,ICFP)
	      IDCP=LVEDGE(2,ICFP)
              RadHeat(IDC)=shigma*tmp(IDC,1)**4
              WN(1)=SFAREA(1,ICFL)
	      WN(2)=SFAREA(2,ICFL)
	      WN(3)=SFAREA(3,ICFL)
	      ven=0.d0
              do Iseg=1,nseg
		ii0=(Iseg-1)*NALLCV
		idiv0=(Iseg-1)*NDIV
          	FRAD(1:NDIV)=radinten(idiv0+1:idiv0+ndiv,ICV)
		CALL WallInciEnergy(radwalltype(nb),NDIV,
     &	        FRAD,WN,vinc)
		ven=ven+vinc
              enddo
              RadHeatFlux(IDC)=(RadHeat(IDC)-ven)*radprop(2,nb)
	      RadHeat(IDC)=RadHeatFlux(IDC)*SFAREA(4,ICFL)
!	      RadHeatFlux(IDCP)=RadHeatFlux(IDC)
!	      RadHeat(IDCP)=RadHeat(IDC)
              enddo
            ELSE
	      do IBFL=IBFS,IBFE
              ICFP=LCYCSF(IBFL)
	      ICVP=LVEDGE(1,ICFP)
              IDCP=LVEDGE(2,ICFP)
              RadHeat(IDCP)=shigma*tmp(IDCP,1)**4
	      WN(1) = SFAREA(1,ICFP)
	      WN(2) = SFAREA(2,ICFP)
              WN(3) = SFAREA(3,ICFP)
	      ven=0.d0
	      do Iseg=1,nseg
	      ii0=(Iseg-1)*NALLCV
              idiv0=(Iseg-1)*NDIV
	      FRAD(1:NDIV) = radinten(idiv0+1:idiv0+ndiv,ICVP)
	      CALL WallInciEnergy(radwalltype(nb),NDIV,
     &		  FRAD,WN,vinc)
              ven=ven+vinc
	      end do
              RadHeatFlux(IDCP)=(RadHeat(IDCP)-ven)*radprop(2,nb)
	      RadHeat(IDCP)=RadHeatFlux(IDCP)*SFAREA(4,ICFP)
	      ICFL=LBC_SSF(IBFL)
	      ICV=LVEDGE(1,ICFL)
	      IDC=LVEDGE(2,ICFL)
	      RadHeatFlux(IDC)=RadHeatFlux(IDCP)
	      RadHeat(IDC)=RadHeat(IDCP)
	      enddo
	    ENDIF
	  end if
	ELSE IF(kd.eq.kdprdc) THEN
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
	  IDCA=LVEDGE(2,ICFL)
	  ICFP=LCYCSF(IBFL)
	  ICVB=LVEDGE(1,ICFP)
	  IDCB=LVEDGE(2,ICFP)
	  RadHeat(IDCB)=0.d0
	  WN(1:3)=SFAREA(1:3,ICFL)		!directing IDC
	  RadHeat(IDCA)=0.d0
	  do Iseg=1,nseg
	  ii0=(Iseg-1)*NALLCV
	  idiv0=(Iseg-1)*NDIV
	  FRAD(1:NDIV) = radinten(idiv0+1:idiv0+ndiv,ICVA)
	  CALL WallInciEnergy(1,NDIV,
     &	  FRAD,WN,vinc)
          RadHeat(IDCA)=RadHeat(IDCA)+vinc
	  enddo
	  WN(1:3) = SFAREA(1:3,ICFP)		!directing IDC
	  RadHeat(IDCB)=0.d0
	  do Iseg=1,nseg
          ii0=(Iseg-1)*NALLCV
	  idiv0=(Iseg-1)*NDIV
	  FRAD(1:NDIV) = radinten(idiv0+1:idiv0+ndiv,ICVB)
	  CALL WallInciEnergy(1,NDIV,
     &	  FRAD,WN,vinc)
          RadHeat(IDCB)=RadHeat(IDCB)+vinc
	  enddo
	  RadHeatFlux(IDCA)=(RadHeat(IDCA)-RadHeat(IDCB))
	  RadHeat(IDCA)=RadHeatFlux(IDCA)*SFAREA(4,ICFL)
	  RadHeatFlux(IDCB)=-RadHeatFlux(IDCA)
	  RadHeat(IDCB)=-RadHeat(IDCA)
	  end do
	ELSE !wall, or other kind wall		
	  IF(kdt.eq.kttrns) then
            do IBFL=IBFS,IBFE
	    ICFL=LBC_SSF(IBFL)
	    ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            RadHeat(IDC)=shigma*tmp(ICV,1)**4
            WN(1:3) = SFAREA(1:3,ICFL)
	    ven=0.d0
	    do Iseg=1,nseg
            ii0=(Iseg-1)*NALLCV
	    idiv0=(Iseg-1)*NDIV
	    FRAD(1:NDIV) = radinten(idiv0+1:idiv0+ndiv,ICV)
	    CALL WallInciEnergy(radwalltype(nb),NDIV,
     &	    FRAD,WN,vinc)
	    ven=ven+vinc
            end do
            RadHeatFlux(IDC)=(RadHeat(IDC)-ven)*radprop(2,nb)
	    RadHeat(IDC)=RadHeatFlux(IDC)*SFAREA(4,ICFL)
	    end do
	  else
	    do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
	    ICV=LVEDGE(1,ICFL)
	    IDC=LVEDGE(2,ICFL)
            RadHeat(IDC)=shigma*tmp(IDC,1)**4
            WN(1:3) = SFAREA(1:3,ICFL)		!directing IDC
            ven=0.d0
	    do Iseg=1,nseg
            ii0=(Iseg-1)*NALLCV
	    idiv0=(Iseg-1)*NDIV
	    FRAD(1:NDIV) = radinten(idiv0+1:idiv0+ndiv,ICV)
	    CALL WallInciEnergy(radwalltype(nb),NDIV,
     &	    FRAD,WN,vinc)
	    ven=ven+vinc
	    enddo
            RadHeatFlux(IDC)=(RadHeat(IDC)-ven)*radprop(2,nb)
	    RadHeat(IDC)=RadHeatFlux(IDC)*SFAREA(4,ICFL)
            enddo
          endif
	ENDIF

1000	CONTINUE
!
! only for benchmark
!	
! --- !
      twheat=0.d0
      DO nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      if(.not.mat_cal(IIMAT)) cycle
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      do IBFL=IBFS,IBFE
      ICFL=LBC_SSF(IBFL)
      ICV=LVEDGE(1,ICFL)
      IDC=LVEDGE(2,ICFL)
      twheat(nb,3)= twheat(nb,3)+SFAREA(4,ICFL)
      twheat(nb,1)= twheat(nb,1)+RadHeat(IDC)
      enddo
      if(twheat(nb,3).gt.0.d0) then
        twheat(nb,1)=twheat(nb,1)/twheat(nb,3)
      endif
      ENDDO
!      if(my_rank.eq.0) write(ifll,*) 'Radiation Calculation by FVM'
!
      RETURN
!-----------------------------
!-------> end of CalRadQFVM()
!-----------------------------
      END SUBROUTINE CalRadQFVM

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE 	GetDivAngle(IDIV,NDIV,NV)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	IMPLICIT NONE     !REAL*8(A-H,O-Z)

        INTEGER,INTENT(INOUT) :: IDIV,NDIV

	REAL*8,INTENT(INOUT) :: NV(3)
	REAL*8	:: NV8(3),NV24(3,3),NV48(3,6),NV80(3,10)
        INTEGER :: nd,imod,kv,k1,jv,j1,iv
	


	DATA  NV8  /0.5773503, 0.5773503,  0.5773503 /

	DATA  NV24 /0.2958759, 0.2958759,  0.9082483, 
     &			0.9082483, 0.2958759,  0.2958759, 
     &			0.2958759, 0.9082483,  0.2958759 /

	DATA  NV48 /0.1838670, 0.1838670,  0.9656013,			!1
     &			0.6950514, 0.1838670,  0.6950514,			!2
     &			0.9656013, 0.1838670,  0.1838670,			!3
     &			0.1838670, 0.6950514,  0.6950514,			!4
     &			0.6950514, 0.6950514,  0.1838670,			!5
     &			0.1838670, 0.9656013,  0.1838670 /			!6

	DATA  NV80 /0.1422555, 0.1422555,  0.9795543,			!1
     &			0.5773503, 0.1422555,  0.8040087,			!2
     &			0.8040087, 0.1422555,  0.5773503,			!3
     &			0.9795543, 0.1422555,  0.1422555,			!4
     &			0.1422555, 0.5773503,  0.8040087,			!5
     &			0.5773503, 0.5773503,  0.5773503,			!6
     &			0.8040087, 0.5773503,  0.1422555,			!7
     &			0.1422555, 0.8040087,  0.5773503,			!8
     &			0.5773503, 0.8040087,  0.1422555,			!9
     &			0.1422555, 0.9795543,  0.1422555 /			!10
	CHECK_NDIV: SELECT CASE (NDIV)
	CASE (8)
	NV(1)=NV8(1)
	NV(2)=NV8(2)
	NV(3)=NV8(3)
	ND = 1
	CASE (24)
	IMOD = MOD(IDIV-1,3)+1
	NV(1) = NV24(1,IMOD)
	NV(2) = NV24(2,IMOD)
	NV(3) = NV24(3,IMOD)
	ND = 3
	CASE (48)
	IMOD = MOD(IDIV-1,6)+1
	NV(1) = NV48(1,IMOD)
	NV(2) = NV48(2,IMOD)
	NV(3) = NV48(3,IMOD)
	ND = 6
        CASE (80)
	IMOD = MOD(IDIV-1,10)+1
	NV(1) = NV80(1,IMOD)
	NV(2) = NV80(2,IMOD)
	NV(3) = NV80(3,IMOD)
	ND = 10
	END SELECT CHECK_NDIV
!
	KV=(IDIV-1)/(ND*4)
	K1=IDIV-KV*ND*4
	JV=(K1-1)/(ND*2)
	J1=K1-JV*ND*2
	IV = (J1-1)/ND
	NV(1)=NV(1)*(1.0-2.0*IV)
	NV(2)=NV(2)*(1.0-2.0*JV)
	NV(3)=NV(3)*(1.0-2.0*KV)
!
	RETURN
	END
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE 	GetDivWeight(IDIV,NDIV,weight)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	IMPLICIT NONE    !REAL*8(A-H,O-Z)
        INTEGER,INTENT(INOUT) :: IDIV,NDIV
	REAL*8,INTENT(INOUT) :: weight
	REAL*8	:: W8,W24(3),W48(6),W80(10)
        INTEGER :: imod


	DATA  W8  /1.5707963d0/

	DATA  W24 / 0.5235987d0, 0.5235987d0, 0.5235987d0 /
     
	DATA  W48 / 0.1609517d0, 0.3626469d0, 0.1609517d0,
     &			0.3626469d0, 0.3626469d0, 0.1609517d0 /


	DATA  W80 /0.1712359d0, 0.0992284d0, 0.0992284d0, 0.1712359d0,
     &		0.0992284d0, 0.4617179d0, 0.0992284d0, 0.0992284d0,
     &			0.0992284d0, 0.1712359d0 / 

	

	CHECK_NDIV: SELECT CASE (NDIV)
		
		CASE (8)
			weight = W8

		CASE (24)
			IMOD = MOD(IDIV-1,3)+1
			weight = W24(IMOD)

		CASE (48)
			IMOD = MOD(IDIV-1,6)+1
			weight = W48(IMOD)

		CASE (80)
			IMOD = MOD(IDIV-1,10)+1
			weight = W80(IMOD)

	END SELECT CHECK_NDIV

	RETURN

!------> end of GetDivWeight()

	END
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ScaPhaseFun(IDIV,JDIV,NDIV,IScaKind,beta,value)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT NONE     !REAL*8(A-H,O-Z)
!      IMPLICIT NONE
      integer,intent(in)  :: NDIV,JDIV,IDIV,IScaKind
      
      REAL*8,intent(inout)  :: beta,value


      REAL*8,PARAMETER :: PAI=3.14169265d0
	
      REAL*8 ::	NV1(3),NV2(3) ,Qsca=0.d0
      REAL*8 ::	weight,cosv
      
	
      CALL GetDivAngle(IDIV,NDIV,NV1)
      CALL GetDivAngle(JDIV,NDIV,NV2)
      CALL GetDivWeight(JDIV,NDIV,weight)
	
      cosv = NV1(1)*NV2(1)+NV1(2)*NV2(2)+NV1(3)*NV2(3)
      CALL FVMScaValue(cosv,IScaKind,beta,Qsca)
	
      value = weight*Qsca

      RETURN

!------> end of ScaPhaseFun()

      END SUBROUTINE ScaPhaseFun

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE WallReflectEmi(WType,NDIV,FRAD,WN,TT,WRAD,graydeg,
     &						   sumcoef)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	use module_radsxf, only  : DivW,DivA,ID,WK
	
        IMPLICIT NONE
        
	INTEGER ,intent(in) :: WType
        integer,intent(in)  :: NDIV
        real*8 ,intent(in)  :: WN(3),FRAD(NDIV)
        real*8 ,intent(inout)  :: WRAD(NDIV),sumcoef
        real*8 ,intent(in)  :: graydeg,TT
!
! --- 
!
	INTEGER	:: IDw,IDIV,jdiv
	REAL*8 :: Wmax,NV1(3),NV2(3) ,cosv1,wref,cosv2,qref,scsumcoef
	REAL*8 :: value,weight,cosv,T4
!
	REAL*8,PARAMETER :: PAI=3.14169265d0
     &                  ,SHIGMA=5.67d-8, FUSHE=1.8048d-8
!
!------------------------------------
!WN is outward normal direction 
!------------------------------------
	T4=TT**4
!---------------
!diffuse wall 
!---------------
	IF(WType.EQ.1) THEN
	  DO IDIV=1,NDIV
          WRAD(IDIV)=0.d0
          cosv1=-(DivA(1,IDIV)*WN(1)+DivA(2,IDIV)*WN(2)
     &         +DivA(3,IDIV)*WN(3))
          if(cosv1.le.0.d0) CYCLE
          
          WRAD(IDIV)=graydeg*FUSHE*T4*sumcoef
          if(graydeg.eq.1.d0) cycle
          wref=0.d0
	  DO JDIV=1,NDIV
          cosv=DivA(1,JDIV)*WN(1)
     &        +DivA(2,JDIV)*WN(2)
     &	      +DivA(3,JDIV)*WN(3)
          cosv=max(cosv,0.d0)
          wref=wref+DivW(JDIV)*FRAD(JDIV)*cosv
          ENDDO
          WRAD(IDIV)=WRAD(IDIV)+wref*(1.d0-graydeg)/PAI

	  ENDDO
	  RETURN
        ENDIF
!
!mirror wall
!
	ID=0
	IDw=0
	Wmax=0.d0
	WK=0.d0
	IF(WType.EQ.2) THEN
          DO IDIV=1,NDIV
          WRAD(IDIV)=0.d0
          cosv=(DivA(1,IDIV)*WN(1)+DivA(2,IDIV)*WN(2)
     &	      +DivA(3,IDIV)*WN(3))
          cosv1=-cosv
          IF(Wmax.LT.cosv1) THEN
	    Wmax=cosv1
	    IDw=IDIV
          ENDIF
          if(cosv.le.0.d0) CYCLE
          NV1(:)=DivA(:,IDIV)
          CALL ReflectDirect(NV1,WN)
          DO JDIV= 1,NDIV
          cosv2=NV1(1)*DivA(1,JDIV)+NV1(2)
     &	  *DivA(2,JDIV)+NV1(3)*DivA(3,JDIV)
	  IF(cosv2.gt.WK(IDIV)) then
            WK(IDIV)=cosv2
            ID(IDIV)=JDIV
          ENDIF
          END DO
          IDw=IDIV   !zhang????
	  WRAD(IDw)=graydeg*shigma*T4*sumcoef
          END DO 
	  if(graydeg.eq.1.d0) return
	  DO JDIV=1,NDIV
          IDIV=ID(JDIV)
          IF(IDIV.GT.0) WRAD(IDIV)=WRAD(IDIV)
     &		 +FRAD(JDIV)*(1.0-graydeg) !NOTE: no "*weight/pai" !!
          ENDDO
        ENDIF
!
	IF(WType.EQ.3) THEN	
	  DO IDIV=1,NDIV
          WRAD(IDIV)=0.d0
          cosv1=-(DivA(1,IDIV)*WN(1)+DivA(2,IDIV)*WN(2)
     &    +DivA(3,IDIV)*WN(3))
          if(cosv1.le.0.d0) CYCLE
!------------------------------------------------------
!cosv1 is correction to emission
!------------------------------------------------------
          WRAD(IDIV) = graydeg*FUSHE*T4*cosv1*sumcoef
          if(graydeg.eq.1.d0) cycle
          wref=0.d0
          DO JDIV= 1,NDIV
          cosv=DivA(1,JDIV)*WN(1)+DivA(2,JDIV)*WN(2)
     &		  +DivA(3,JDIV)*WN(3)
          IF(cosv.gt.0) THEN
            Qref=1.d0	!/cosv*cosv
!------------------------------------------------------
!/cosv is correction to reflectivity
!*cosv is correction to surface-area
!------------------------------------------------------
            wref=wref+DivW(JDIV)*FRAD(JDIV)*Qref
          ENDIF
          ENDDO
!------------------------------------------------------
!cosv1 is correction to reflected energy intensity
!------------------------------------------------------
          WRAD(IDIV)=WRAD(IDIV)+wref*(1.d0-graydeg)*cosv1/PAI
	  ENDDO
        ENDIF
        
	RETURN
!-----------------------------------
!------> end of WallReflectEmi()
!-----------------------------------
	END

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE WallInciEnergy(WTYPE,NDIV,FRAD,WN,vinc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	use module_radsxf, only  : DivW,DivA
!					
	IMPLICIT NONE     !REAL*8(A-H,O-Z)
!
        INTEGER,INTENT(INOUT) :: WType,NDIV
        real*8,INTENT(INOUT) :: FRAD(NDIV),WN(3)
        real*8,INTENT(INOUT) :: vinc

        real*8,PARAMETER :: SHIGMA=5.67d-8,
     &                    PAI=3.14159265d0,FUSHE=1.8048D-8
	REAL*8 :: weight,cosv,graydeg
	REAL*8 :: NV1(3),NV2(3)
        INTEGER :: jdiv
!--------------------
	vinc=0 
!--------------------
! --- diffuse wall 
!--------------------
	IF(WType.EQ.1.OR.WType.EQ.2) THEN	
          DO JDIV=1,NDIV
          cosv=DivA(1,JDIV)*WN(1)+DivA(2,JDIV)*WN(2)
     &			   +DivA(3,JDIV)*WN(3)
          cosv=max(cosv,0.d0)
          vinc=vinc+DivW(JDIV)*FRAD(JDIV)*cosv
          ENDDO
	ENDIF
!
	IF(WType.EQ.3) THEN	
          DO JDIV=1,NDIV
          cosv=DivA(1,JDIV)*WN(1)+DivA(2,JDIV)*WN(2)
     &			   +DivA(3,JDIV)*WN(3)
          cosv=max(cosv,0.d0)
          vinc=vinc+DivW(JDIV)*FRAD(JDIV)*cosv*cosv
          ENDDO
        ENDIF
!
	RETURN
!
	END SUBROUTINE WallInciEnergy
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE FVMScaValue(cosv,ptype,beta,fani)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	IMPLICIT NONE  !REAL*8(A-H,O-Z)
!
        real*8 ,INTENT(INOUT) :: cosv,beta,fani
	INTEGER,INTENT(IN) :: PTYPE
!
	real*8,PARAMETER :: SHIGMA=5.67d-8,
     &                      PAI=3.14159265d0,FUSHE=1.8048D-8
	REAL*8 :: f,Q0,Qout,sinv,eta

!---------------------------------------------------------
!PType = 1 :  Isotropic / linear-anisotropic scattering		
!      = 2 :  Large-Diffuse-Particle Scattering
!      = 3 :  Rayleigh Scattering
!---------------------------------------------------------
        beta=0.d0
        if(PTYPE==1) then
          fani=1.d0+beta*cosv
        elseif(PTYPE==2) then
          cosv=min(cosv,1.d0)
	  cosv=max(cosv,-1.d0)
	  sinv=DSQRT(1.d0-cosv*cosv)
	  eta=Dacos(cosv)
          fani=8.d0/(3.d0*PAI)*(sinv-eta*cosv)
        elseif(PTYPE==3) then
          fani=0.75d0*(1.d0+cosv*cosv)
        endif

!	CHECK_PTYPE1: SELECT CASE (PTYPE)
!	CASE (1)
!	fani=1.d0+beta*cosv
!	CASE (2)
!	cosv=min(cosv,1.d0)
!	cosv=max(cosv,-1.d0)
!	sinv=DSQRT(1.d0-cosv*cosv)
!	eta=Dacos(cosv)
!	fani=8.d0/(3.d0*PAI)*(sinv-eta*cosv)
!	CASE (3)
!	fani=0.75d0*(1.d0+cosv*cosv)
!	END SELECT CHECK_PTYPE1
!
	RETURN
	END SUBROUTINE FVMScaValue
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE WallAnisPower(cosv,WType,Qout)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	IMPLICIT none

        INTEGER,INTENT(INOUT) ::  WType
        real*8,INTENT(INOUT) :: cosv,Qout
        REAL*8, PARAMETER :: PAI=3.14159265d0

!	INTEGER WType
!	REAL*8	cosv,Qout

	
	CHECK_WTYPE1: SELECT CASE (WTYPE)
		CASE (1)
			Qout = 1.d0/PAI

		CASE (2)
			!mirror surface, if vv!=0, the anisotropic reflection is zero
			Qout=1.d0
	
		CASE (3)
			
			!directional-gray surface, alfa = alfa0*cos(seta)
			!							  alfa0
			!						alfa	|	  alfa
			!					     `		|  	  /
			!					      `		|    /
			!					       `	|   /
			!					        `	seta
			!						     `	| /
			!					__________`_|/__________
			!					///////////////////////////		wall
	
			Qout = Qout*cosv			!????
		
		
	END SELECT CHECK_WTYPE1


	RETURN

!------>
	END
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE ReflectDirect(Vin,NV)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	IMPLICIT none  !IMPLICIT REAL*8(A-H,O-Z)

	
	real*8,INTENT(INOUT) :: Vin(3),NV(3)
        real*8  :: f1,f2
	
	!NV is the outward normal direction
	!Vin is the incident direction, and at return, it is the relfected direction
	
	f1 = Vin(1)*NV(1)+Vin(2)*NV(2)+Vin(3)*NV(3)
	f2 = NV(1)*NV(1)+NV(2)*NV(2)+NV(3)*NV(3)
	
	Vin(:)=Vin(:)-2.d0*f1/f2*NV(:)
	
		!							    
		!						Vin		|	  Vin`
		!					      `		|  	  /
		!					       `	|    /
		!					        `	|   /
		!					         `	!  /
		!						      `	| /
		!					____________|/_________
		!								|
		!								|
		!								|
		!								NV


	RETURN

!------>
	END
!
!
!-----------------------------------------------------------------------
!	Developed for radiation heat transfer by MC/ZONE methods
!-----------------------------------------------------------------------
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GetHeatFlux(iter,NPE,NBOUN,NMAT,NALLCV,tmp,
     &    MAT_CV,RadHeatFlux,RadHeat)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	
      use module_dimension, only  : MXALLCV,MXALLCV_RAD
      use module_radsxf, only  : RDValue,RDIndex,RDN1,RDN2,RDId,
     &                           NRD,MRD,NRD0,MRD0,MID,StackOpt
      use module_radsxf, only  : EB4,QSUM
      use module_radsxf, only  : VJOBMD,NJOBMD
      use module_radsxf, ONLY  : IndexToRD,CoefToRD,MatBounEnd
      use module_radsxf, ONLY  : IndexToExpt,RaToExpt,ExptID,EBcomm
      use module_radsxf, only  : NPALL,NDatExpt,NumTran,NumTranIn,
     &                           NMatBoun
      use module_boundary,only : ktneum,kttrns,kyneum,kytrns,kdintr,
     &                           kdcvd,kydirc,kdilet,nbcnd,kdbcnd
!------- for benchmark only 
      use module_radsxf, only  : NWDAT,WDAT,WPNUM
      use module_time,   only :  iters,itere
!array VJOBMD ::
!(:,1)	:	Area / Volume
!(:,2)  :	Property
!(:,3)  :	QW
!(:,4)  :	EB (SHIGMA*T^4)
!array NJobMD ::
!(:,1)  :   ird <-- icv
!(:,2)	:   ird --> icv	
!(:,3)  :   mark RAD export communication vertex
!
      IMPLICIT NONE
!
      integer,intent(in)    :: iter,NPE,NBOUN,NMAT,NALLCV,
     &                         MAT_CV(MXALLCV)
      real*8 ,intent(out)   :: RadHeatFlux(MXALLCV_RAD)
      real*8 ,intent(out)   :: RadHeat    (MXALLCV_RAD)
      real*8 ,intent(in)    :: tmp(MXALLCV)




      real*8,PARAMETER :: PAI=3.141592653589793d0,SHIGMA=5.67d-8
      integer :: ICV,I,IV,IRD,ID,J,ID1,J1,IP,idt
      real*8  :: revdt,sumqw,sumtg,qw1,tg1,qw2,tg2
!
!      REAL*8  tmp(MXALLCV)
!      INTEGER MAT_CV(MXALLCV)
      CHARACTER*80 FileNam
C-----------------------------------------------------------------------------------------C
C	Energy conservation equation for unknown element:
C			PiiXi = sum(Dij*Xj) + Bi	-----
C		  			 |
C						 |
C			Aii*Xi + sum(Aij*Xj) = Bi	<----
C			Aij = -Dij		(for i!=j)
C			Aii = Pii-Dii
C			Bi = HSource			
C
C	RadHeat for known element:
C	RadHeat = Aii*shigma*Ti^4 + sum(Aij*shigma*Tj^4) - Bi		
C	dimension:	RadHeat/AreaVol	[W/m^3] for fluid
C			[W/m^2] for boundary
C			
C	Positive RadHeat means Qgo > Qcome
C-----------------------------------------------------------------
C	1. Value (Emissive power) transfer from fluid to radiation
C-----------------------------------------------------------------
	IF(iter.EQ.iters) THEN
	  DO ICV=1,NALLCV
	  EBcomm(ICV)=shigma*tmp(icv)**4
	  END DO
	ENDIF
!
	DO I=1,MatBounEnd(NMAT+NBOUN)
	  ICV=NJobMD(I,1)
	  VJobMD(I,4)=shigma*tmp(icv)**4
	END DO
!
	DO I=MatBounEnd(NBOUN+NMAT)+1,NPALL
	  ICV=NJobMD(I,1)
          VJobMD(I,4) = EBcomm(ICV)
	END DO
!-----------------------------------------------------
! --- VJobMD(NALLCV+1:NPALL,4) is from communication
!-----------------------------------------------------
	EB4=0.d0
	DO 100 I=1,NumTran
	IV = IndexToRD(1,I)
	IRD = IndexToRD(2,I)
	EB4(IRD)=EB4(IRD)+CoefToRD(I)*VJobMD(IV,4)
100	CONTINUE
C --------------------------------------
C	2. radiation heat transfer 
C --------------------------------------
	VJOBMD(:,3)=0.d0
	IF(StackOpt(1:4).EQ.'LTRI') THEN
	  DO 200 I = 1,NRD
	  ID = RDIndex(I)
	  VJOBMD(I,3)=RDValue(ID+I)*EB4(I)
          DO J=1,I-1
   	  VJOBMD(I,3)=VJOBMD(I,3)+RDValue(ID+J)*EB4(J)
		!Aij < 0
	  END DO
!				
	  DO J=I+1,NRD
	  ID1 = RDIndex(J)
	  VJOBMD(I,3) = VJOBMD(I,3)+RDValue(ID1+I)*EB4(J)
		!Aji = Aij < 0
          END DO
200	  CONTINUE

	ELSE
	  DO 220 I = 1,NRD
	  ID = RDIndex(I)
	  VJOBMD(I,3) = RDValue(ID+1)*EB4(I)
          DO J=ID+2,RDIndex(I+1)
	  J1 = RDId(J)
	  VJOBMD(I,3) = VJOBMD(I,3)+RDValue(J)*EB4(J1)
		!Aij < 0
	  END DO
220	  CONTINUE
	END IF


C-----------------------------------------------------
C	3. radiation heat generation for fluid vertex
C-----------------------------------------------------

	VJobMD(:,4)=0.d0
	DO 300 I=1,NumTran
	IV = IndexToRD(1,I)
	IRD = IndexToRD(2,I)
	VJobMD(IV,4)=VJobMD(IV,4)+CoefToRD(I)*VJOBMD(IRD,3)
300	CONTINUE
	
	RadHeat = 0.d0
	DO 310 ICV=1,NALLCV
	IV = NJobMD(ICV,2)
	IF(IV.GT.0) RadHeat(ICV) = VJobMD(IV,4)
310	CONTINUE
	


C---------------------------------------------------------------
C   4	. imaginary emissive power for export vertex 
!         / communication in parallel computation
C---------------------------------------------------------------
	IF(NDatExpt.LE.0) GOTO 5678
	
	DO IRD=1,NRD
	IF(ExptID(IRD).GT.0) EB4(IRD)=0.d0			!=1, or 2
	END DO

	VJobMD(:,3)=0.d0
	IF(StackOpt(1:4).EQ.'LTRI') THEN
	DO 400 I = 1,NRD
	IF(ExptID(I).EQ.0) CYCLE
	ID = RDIndex(I)
	VJOBMD(I,3) = RDValue(ID+I)*EB4(I)
	DO J=1,I-1
	VJOBMD(I,3) = VJOBMD(I,3)+RDValue(ID+J)*EB4(J)
!Aij < 0
	END DO
!
	DO J=I+1,NRD
	ID1 = RDIndex(J)
	VJOBMD(I,3) = VJOBMD(I,3)+RDValue(ID1+I)*EB4(J)
                       !Aji = Aij < 0
	END DO
			
400	CONTINUE

	ELSE

	DO 420 I = 1,NRD
	IF(ExptID(I).EQ.0) CYCLE
	ID = RDIndex(I)
	VJOBMD(I,3) = RDValue(ID+1)*EB4(I)
	DO J=ID+2,RDIndex(I+1)
	J1 = RDId(J)
	VJOBMD(I,3) = VJOBMD(I,3)+RDValue(J)*EB4(J1)
!Aij < 0
	END DO
	
420	CONTINUE

	END IF

	VJobMD(:,4)=0.d0
	DO 430 I=1,NDatExpt
        IP = IndexToExpt(1,I)
	IRD = IndexToExpt(2,I)
	VJobMD(IP,4) = VJobMD(IP,4) + RAToExpt(I)*VJOBMD(IRD,3)
430	CONTINUE
	

	EBComm=0.d0
	DO ICV=1,NALLCV
	EBcomm(ICV)=shigma*tmp(icv)**4
	END DO
	DO 440 ICV=1,NALLCV
	IP = MAT_CV(ICV)
	IF(NJobMD(IP,3).EQ.1) EBcomm(ICV)=-VJobMD(IP,4)
440	CONTINUE
!
C---------------------------------------------------------------------------
C   5	. custom treatment / heat flux divergence and surface heat flux (qw)
C---------------------------------------------------------------------------
!
5678	RadHeatFlux(:)=0.d0		!use as heat flux Qw
	DO ICV=1,NALLCV
	IV = NJobMD(ICV,2)
	IF(IV.GT.0) RadHeatFlux(ICV) = RadHeat(ICV)/VJobMD(IV,1)
	END DO

	GOTO 9999
!for benchmark computation
	if(itere-iters.gt.101) then
	  idt=100
	else
	  idt=itere-iters-1
	end if

	if(iter.LE.itere-idt) RETURN
		
	revdt = 1.d0/idt

	DO I=1,NPALL
	J = WDAT(I,1)			!qw, ave_y
	IF(J.GT.0) THEN
	ICV =  NJobMD(I,1)
	QSUM(J,1)=QSUM(J,1)+VJOBMD(ICV,3)/dble(WPNUM(J,1))
        END IF
		
	J = WDAT(I,2)
	IF(J.GT.0) THEN			!t, ave_y
	ICV =  NJobMD(I,1)
	QSUM(J,2)=QSUM(J,2)+tmp(ICV)/dble(WPNUM(J,2))
	END IF

	J = WDAT(I,3)			!qw, ave_z
	IF(J.GT.0) THEN
	ICV =  NJobMD(I,1)
	QSUM(J,3)=QSUM(J,3)+VJOBMD(ICV,3)/dble(WPNUM(J,3))
	END IF
		
	J = WDAT(I,4)
	IF(J.GT.0) THEN			!t,ave_z
	ICV =  NJobMD(I,1)
	QSUM(J,4)=QSUM(J,4)+tmp(ICV)/dble(WPNUM(J,4))
	END IF

	END DO

	
	if(iter.EQ.itere) then
	sumqw = 0.d0
	sumtg = 0.d0
	FileNam='HotWall.dat'
	open(10,file=FileNam)
	write(10,'(I4)') NWDAT
	DO I=1,NWDat
		qw1 = QSUM(I,1)*revdt
		tg1 = QSUM(I,2)*revdt
		qw2 = QSUM(I,3)*revdt
		tg2 = QSUM(I,4)*revdt

		WRITE(10,'(I5,4F12.4)') I,qw1,tg1,qw2,tg2
		sumqw = sumqw + qw2		!+qw1
		sumtg = sumtg + tg2		!+tg1
		END DO
		write(10,*) 'Global average'
		WRITE(10,'(4F12.4)') sumqw/NWDat,sumtg/NWDat
		close(10)
		
	end if





9999	RETURN

!---------> end of GetHeatFlux()

	END SUBROUTINE GetHeatFlux

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE GetProWSGG(deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rho,vel,prs,pp0,cps,cp,yys,
     &  MAT_NO,MAT_CV,MAT_INDEX,
     &  MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  GWAbsorb,PAbsorb,PScatter,sumcoef,
     &  tmpbnd,tmp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd,ktneum,ktdirc,kttrns
 
      use module_species,only  : spcnam,wm,r_wm,Ri,spcno
      use module_boundary, only : radprop, radwalltype
      use module_material,  only: radmat,radfludflag
      use module_rad
!      use module_radsxf,   only : GWAbsorb, PAbsorb, PScatter
      use module_radsxf,   only : DivW,DivA
      use module_radsxf,   only : raddum   !,sumcoef,raddum
      use module_usersub,  only : src_rad,usryes	
      use module_metrix,only    : rcomp
      use module_model,only     : mach0,idrdp,incomp
      
!
      implicit none	

!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(inout) :: rho       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3  )
      real*8 ,intent(inout) :: prs       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV  ,2)
      REAL*8 ,INTENT(INOUT) :: cps   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp    (       MXALLCV)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(inout) :: tmp       (   MXALLCV  ,2)
      real*8 ,intent(inout) :: yys       (   MXALLCV,MXcomp,2)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      integer,intent(out)   :: ierror
!------------------------
! --- [local entities]
!------------------------
      REAL*8  :: yref(ncomp),vall,ymol,pfenbu,kosesb,ratio,
     &		 Cabs,yref_all,pref,Tref,tmin,tmax,tw
      REAL*8  :: CabsSeg(0:100),ffenbu(0:100),CabsLog(0:100),
     &		 CabsAve(100,0:Ncomp),fave(100,0:Ncomp),
     &		 ymol_all(0:5),peff(0:5)
      REAL*8  :: aa1,aa2,aii,bjj,cmin,cmax,lenth
      INTEGER :: Icomp,I,J,K,ii0,Iseg,Nseg,Igauss,IWSGG=5
      REAL*8  :: ymol_t
      
      INTEGER :: MarkNo(2),iimat,imat,icvs,icve,icvl,icv,ii,nb,
     &           kd,kdt,kdy
      INTEGER :: ibfs,ibfe,ibfl,icfl,idc
!----------------------------------------------------------------------------
!for media
c	radmat(:,:)
c	!1: gabsorb, 2: pabsorb, 3: pscatter, 4:scabeta, 5: ptcdiamter
c	radfludtype(:,:)
c	!1: gastype, 2: scatype, 3: ptccal, 4: fgpflag
!----------------------------------------------------------------------------
!gastype
!1: gray_gas,  2: real_gas
!scatype	
!		    1,	Linear-anisotropic scatter / isotropic scatter, OK
!			2,  large-diffuse particle scatter, OK
!			3,	Rayleigh scatter, OK
!			4,	Mie scatter
!			5,  user defined scatter
!----------------------------------------------------------------------------
!for wall
!radprop(1,nb)	:	 rademi
!radprop(2,nb)	:	 radabsorb
!radwalltype(nb) :	 radtype
!----------------------------------------------------------------------------
!
      call FFRABORT(1,'MSG: NOT support WSGG model')
      if(incomp==idrdp) then
        call FFRABORT(1,'ERR: WSGG NOT support incompressible')
      endif
!        

      GWAbsorb(:)=0.d0
      PAbsorb(:)=0.d0
      PScatter(:)=0.d0
      sumcoef(:)=0.d0
      TotCoef(:)=0.d0
!----------------------------------------------------------------------------
c	write(*,*) 'I am here, Ngauss=',Ngauss
c	write(*,*) 'awsgg=',awsgg(:,1:2)
c	write(*,*) 'bwsgg=',bwsgg(:,:,1)
!----------------------------------------------------------------------------
!
!------get user define properties and/or parameters
!----------------------------------------------------------------------------
      IF(src_rad.eq.usryes) then
        do IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.LE.0) cycle		!not fluid
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        ICV=MAT_CV(ICVL)
        do Icomp=1,ncomp
        rcomp(Icomp)=yys(ICVL,Icomp,1)
        enddo
        call user_rad_pro(ICV,ncomp,spcno,CVCENT(:,ICVL),
     &  tmp(ICVL,1),pp0(ICVL,1),rho(ICVL,1),rcomp,
     &	GWAbsorb(ICVL),PAbsorb(ICVL),PScatter(ICVL))
        do Icomp=1,ncomp
          yys(ICVL,Icomp,1)=rcomp(Icomp)
        enddo
        enddo
        enddo
      endif
!
! --- 
!
      do IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
	IMAT=MAT_NO(IIMAT)
	if(IMAT.LE.0) cycle		!not fluid
	ICVS=MAT_CVEXT(IIMAT-1)+1
	ICVE=MAT_CVEXT(IIMAT)
	tref=0.d0
	pref=0.d0
	vall=0.0
	tmin=1.d10
	tmax=0.d0
	do ICVL=ICVS,ICVE
        tref=tref+tmp(ICVL,1)*CVVOLM(ICVL)
        pref=pref+pp0(ICVL,1)*CVVOLM(ICVL)
        vall=vall+CVVOLM(ICVL)
        tmin=min(tmin,tmp(ICVL,1))
        tmax=max(tmax,tmp(ICVL,1))
	enddo
	tref=tref/vall
	pref=pref/vall
	yref_all=0.0

	DO icomp=1,ncomp
        yref(icomp)=0.d0
	do ICVL=ICVS,ICVE
        yref(Icomp)=yref(Icomp)+yys(ICVL,Icomp,1)*rho(ICVL,1)
     &	*r_wm(Icomp)*CVVOLM(ICVL)
        enddo
        yref(Icomp)=yref(Icomp)/vall	!molar concentration
        yref_all=yref_all+yref(Icomp)   !total molar concentration
	ENDDO
!
	if(yref_all.le.1.0d-10) cycle	!radiation is negligable
!
	DO ICVL=ICVS,ICVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	 NOTE:
! The following code is according to the Ref.[5] ASME-JHT, 117, 359-365,1995
! Any uncaustious changement may reduce the accuracy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ymol_all=0.d0
	ymol_t=0.d0
        do Icomp=1,ncomp
	ymol=yys(ICVL,Icomp,1)*r_wm(Icomp)*rho(ICVL,1)
	ymol_t=ymol_t+ymol
	II=ParaFlag(spcno(Icomp))
	ymol_all(II)=ymol_all(II)+ymol
        enddo
        if(sum(ymol_all(1:5)).lt.1.0d-10) CYCLE
        do II=0,5
        peff(II)=pp0(ICVL,1)*(ymol_all(II)/ymol_t)/1.013d5
        enddo
!
        bjj=0.d0
        DO Igauss=1,Ngauss-1
        aii=bwsgg(1,IGauss,1)
	DO I=2,4
	aii=aii+bwsgg(I,IGauss,1)*tmp(ICVL,1)**(I-1)
	ENDDO
	ffenbu(Igauss)=aii
	bjj=bjj+aii
	ENDDO
!
	if(bjj.GT.1.d0+1.0d-2) GOTO 9999
	ffenbu(Ngauss)=max(0.d0,1.d0-bjj)
        DO Igauss=1,Ngauss-1  !zhang-1
	ii0=NALLCV*(Igauss-1)
	GWAbsorb(ii0+ICVL)=awsgg(Igauss,1)*peff(1)
     &		          +awsgg(Igauss,2)*peff(2)
	sumcoef(ii0+ICVL)=ffenbu(Igauss)
        ENDDO
	ii0=NALLCV*(Ngauss-1)
	GWAbsorb(ii0+ICVL)=0.0
	sumcoef(ii0+ICVL)=ffenbu(Ngauss)
	ENDDO	!ICVL
!gray particle
	if(src_rad.ne.usryes) then
	  DO ICVL=ICVS,ICVE
          PAbsorb(ICVL)=radmat(2,IIMAT)
          PScatter(ICVL)=radmat(3,IIMAT)
	  END DO
	endif
	enddo
!	
! --- boundary dummy cv
!
 7777 do nb=1,nbcnd
	IIMAT=MAT_BCIDX(nb,1)
	if(.not.mat_cal(IIMAT)) cycle
	kd=kdbcnd(0,nb)
	kdt=kdbcnd(2,nb)
	kdy=kdbcnd(3,nb)
	IBFS=LBC_INDEX(nb-1)+1
	IBFE=LBC_INDEX(nb)
	if(kdt.eq.ktneum.or.kdt.eq.kttrns) then
          do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  aa1=0.d0
	  DO Iseg=1,ngauss
          ii0=NALLCV*(Iseg-1)
          GWAbsorb(ii0+IDC) = radprop(1,nb)
          sumcoef(ii0+IDC) = sumcoef(ii0+ICV)		!??? OK
          aa1=aa1+sumcoef(ii0+IDC)
          ENDDO
	  ii0=NALLCV*(Ngauss-1)
	  sumcoef(ii0+IDC) = 1.0-aa1
          enddo
	else
!		 
!at this version,
!we suppose only one fluid block, the tref,yref... are not changed
!new version for multi-fluids is under developing.
	  ICFL=LBC_SSF(IBFS)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  tw=tmp(IDC,1)
	  ffenbu = 0.d0
	  if(yref_all.le.1.d-10) then
            raddum(:)=1.0/Ngauss
	  else
            ymol_all=0.0
            do Icomp=1,ncomp
            II=ParaFlag(spcno(Icomp))
            ymol_all(II)=ymol_all(II)+yref(Icomp)
            end do
            bjj=0.d0
            DO Igauss=1,Ngauss-1
            aii=bwsgg(1,IGauss,1)
            DO I=2,4
              aii=aii+bwsgg(I,IGauss,1)*tw**(I-1)
            ENDDO
            raddum(Igauss)=aii
            bjj=bjj+aii
            ENDDO
! --- write(*,*) 'absorb=',Igauss,bjj,raddum(Ngauss)
            if(bjj.GT.1.d0+1.0d-2) GOTO 9999
            raddum(Ngauss)=max(0.d0,1.d0-bjj)
            end if
!	
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            aa1=0.0
            DO Igauss=1,ngauss
            ii0=NALLCV*(Igauss-1)
            GWAbsorb(ii0+IDC)=radprop(1,nb)
            sumcoef(ii0+IDC)=raddum(Igauss)
            ENDDO
            enddo
          endif
	enddo
!		
	do Igauss=1,Ngauss
	ii0=NALLCV*(Igauss-1)
	TotCoef(Igauss)=sum(sumcoef((ii0+1):(ii0+NALLCV)))
	end do
	aa1=sum(TotCoef(1:Ngauss))
	TotCoef=TotCoef/aa1
!
8888	RETURN
!
9999	WRITE(ifle,*) 'ERR in GetProWSGG()'
	WRITE(ifle,*) 'Sum of weight > 1'
!
	call FFRABORT(1,'GetProWSGG')
!
	END SUBROUTINE GetProWSGG

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE GetProSLW(deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rho,vel,prs,pp0,yys,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  GWAbsorb,PAbsorb,PScatter,sumcoef,
     &  tmpbnd,tmp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd,ktneum,ktdirc,kttrns
!
      use module_species,only : spcnam,wm,r_wm,Ri,spcno
      use module_boundary,only : radprop, radwalltype
      use module_material,  only: radmat,radfludflag
      use module_rad
      use module_radsxf,   only : DivW,DivA
      use module_radsxf,   only : raddum   !sumcoef,raddum
      use module_usersub,  only : src_rad,usryes
      use module_metrix,only    : rcomp,yref=>eva_comp
!
      implicit none	
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: rho    (   MXALLCV,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3)
      real*8 ,intent(in)    :: prs       (   MXALLCV,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV,2)
      real*8 ,intent(in)    :: tmpbnd    (   MXssfbc)
      real*8 ,intent(in)    :: tmp       (   MXALLCV,2)
      real*8 ,intent(inout) :: yys       (   MXALLCV,MXcomp,2)
      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
      integer,intent(out)   :: ierror

!
! --- [local entities]
!
!
      REAL*8 ::	vall,ymol,pfenbu,kosesb,ratio,
     &			Cabs,yref_all,ymol_all,pref,Tref,tmin,tmax
      REAL*8 ::	CabsSeg(0:100),ffenbu(0:100),CabsLog(0:100),
     &			CabsAve(100,0:Ncomp),fave(100,0:Ncomp)
      REAL*8 ::	aa1,aa2,aii,bjj,cmin,cmax,lenth,tw,aa3,ymol_t,
     &  cabs1,cabs2,f1
      INTEGER::   Icomp,I,J,K,ii0,Iseg,Nseg,Igauss
      INTEGER :: iimat,imat,icvs,icve,icvl,icv,ii,nb,kd,kdt,kdy,istep
      INTEGER :: ibfs,ibfe,ibfl,icfl,idc
!	
!for media
c	radmat(:,:)
c	!1: gabsorb, 2: pabsorb, 3: pscatter, 4:scabeta, 5: ptcdiamter
c	radfludtype(:,:)
c	!1: gastype, 2: scatype, 3: ptccal, 4: fgpflag
	!gastype
	!1: gray_gas,  2: real_gas
	!scatype	
!1, Linear-anisotropic scatter / isotropic scatter, OK
!2,  large-diffuse particle scatter, OK
!3,	Rayleigh scatter, OK
!4,	Mie scatter
!5,  user defined scatter
				

!for wall
!radprop(1,nb)	:	 rademi
!radprop(2,nb)	:	 radabsorb
!radwalltype(nb) :	 radtype
!
! --- 
!
	GWAbsorb=0.d0
	PAbsorb=0.d0
	PScatter=0.d0
	sumcoef=0.d0

	cmin=1.d10
	cmax=-1.d0
	do I=1,ncomp
        J=spcno(I)
        if(J.EQ.0.OR.ParaFlag(J).NE.1) cycle 
        if(cminmax(2,J).gt.cminmax(1,J)) then
	cmin = min(cmin,cminmax(1,J))
	cmax = max(cmax,cminmax(2,J))
        end if
	end do
!
	Nseg = max(50,min(100,ngauss*16))
	Nseg = max(1,(Nseg/ngauss))*ngauss
	lenth=(log(cmax)-log(cmin))/nseg
!	
	do I=0,nseg
	CabsSeg(I)=exp(log(cmin)+lenth*I)
	end do
!

	IF(src_rad.eq.usryes)	 then
	do IIMAT=1,NMAT    !ICV=1,NCV
	if(.not.mat_cal(IIMAT)) cycle
	IMAT=MAT_NO(IIMAT)
	if(IMAT.LE.0) cycle		!not fluid
	ICVS=MAT_CVEXT(IIMAT-1)+1
	ICVE=MAT_CVEXT(IIMAT)
	do ICVL=ICVS,ICVE
	  ICV=MAT_CV(ICVL)
	  do Icomp=1,ncomp
	  rcomp(Icomp)=yys(ICVL,Icomp,1)
	  end do
	  call user_rad_pro(ICV,ncomp,spcno,CVCENT(:,ICVL),
     &			tmp(ICVL,1),prs(ICVL,1),rho(ICVL,1),rcomp,
     &			GWAbsorb(ICVL),PAbsorb(ICVL),PScatter(ICVL))
!zhang	  do Icomp=1,ncomp
!	  yys(ICVL,Icomp,1)=rcomp(Icomp)
!	  end do
	enddo
	enddo
	endif
	
	do IIMAT=1,NMAT    !ICV=1,NCV
	if(.not.mat_cal(IIMAT)) cycle
	IMAT=MAT_NO(IIMAT)
	if(IMAT.LE.0) cycle		!not fluid
	ICVS=MAT_CVEXT(IIMAT-1)+1
	ICVE=MAT_CVEXT(IIMAT)
1010	tref=0.d0
	pref=0.d0
	vall=0.0
	tmin=1.d10
	tmax=0.d0
	do ICVL=ICVS,ICVE
	  tref=tref + tmp(ICVL,1)*CVVOLM(ICVL)
	  pref=pref + pp0(ICVL,1)*CVVOLM(ICVL)
	  vall=vall+CVVOLM(ICVL)
	  tmin = min(tmin,tmp(ICVL,1))
	  tmax = max(tmax,tmp(ICVL,1))
	end do
	tref=tref/vall
	pref=pref/vall
!------------------------------------------------------
! --- 
!------------------------------------------------------
	yref_all=0.d0
	DO Icomp=1,ncomp
	  yref(Icomp)=0.d0
	  do ICVL=ICVS,ICVE
	  yref(Icomp)=yref(Icomp)+ yys(ICVL,Icomp,1)*rho(ICVL,1)
     &			*r_wm(Icomp)*CVVOLM(ICVL)
	  enddo
	  yref(Icomp)=yref(Icomp)/vall	!molar concentration
	  yref_all=yref_all+yref(Icomp)!total molar concentration
	ENDDO
!
	if(yref_all.le.1.0d-10) cycle	!radiation is negligable
!------------------------------------------
! --- optimize minimum and maximum
!------------------------------------------
	raddum=0.d0
	Iseg=0
	aa3=0.d0
	do while (aa3.LT.5.d-3.and.Iseg.LT.Nseg)
	  aa3=0.d0
	  ymol_all=0
	  do Icomp=1,Ncomp
          if(spcno(Icomp).le.0.or.ParaFlag(spcno(Icomp)).EQ.0) CYCLE
	  if(CabsSeg(Iseg).lt.cminmax(1,spcno(Icomp))) then
	    aa1=0.d0
	    aa2=0.d0
   	  else
	    CALL  GetDistriSLW(Tmin,Tref,pref,yref(Icomp),
     &		yref_all,CabsSeg(Iseg),spcno(Icomp),aa1,
     &		almn(:,:,:,spcno(Icomp)),blmn(:,:,:,spcno(Icomp)))
	    CALL  GetDistriSLW(Tmax,Tref,pref,yref(Icomp),
     &		yref_all,CabsSeg(Iseg),spcno(Icomp),aa2,
     &		almn(:,:,:,spcno(Icomp)),blmn(:,:,:,spcno(Icomp)))
	  end if
	  aa3=aa3+max(aa1,aa2)*yref(Icomp)
	  ymol_all=ymol_all+yref(Icomp)
	  end do
	  if(ymol_all.gt.0.0) aa3=aa3/ymol_all
	  Iseg=Iseg+1
	end do
	Iseg = max(Iseg-2,0)
	cmin=CabsSeg(Iseg)
!
	aa3=1.d0
	Iseg=Nseg
	do while (aa3.GT.0.995.and.Iseg.GT.0)
	  aa3=0.d0
	  ymol_all=0
	  do Icomp=1,Ncomp
	  if(spcno(Icomp).le.0.or.ParaFlag(spcno(Icomp)).EQ.0) CYCLE
	  if(CabsSeg(Iseg).lt.cminmax(1,spcno(Icomp))) then
            aa1=0.d0
	    aa2=0.d0
  	  else
	    CALL  GetDistriSLW(Tmin,Tref,pref,yref(Icomp),
     &		yref_all,CabsSeg(Iseg),spcno(Icomp),aa1,
     &		almn(:,:,:,spcno(Icomp)),blmn(:,:,:,spcno(Icomp)))
	    CALL  GetDistriSLW(Tmax,Tref,pref,yref(Icomp),
     &		yref_all,CabsSeg(Iseg),spcno(Icomp),aa2,
     &		almn(:,:,:,spcno(Icomp)),blmn(:,:,:,spcno(Icomp)))
	  end if
	  aa3=aa3+min(aa1,aa2)*yref(Icomp)
	  ymol_all=ymol_all+yref(Icomp)
	  end do
	  if(ymol_all.gt.0.0) aa3=aa3/ymol_all
	  Iseg=Iseg-1
	end do
	Iseg = min(Iseg+2,Nseg)
	cmax=CabsSeg(Iseg)
	Nseg = max(50,min(100,ngauss*16))
	Nseg = max(1,Nseg/ngauss)*ngauss
	lenth=(log(cmax)-log(cmin))/nseg
	do Iseg=0,nseg
	CabsSeg(Iseg)=exp(log(cmin)+lenth*Iseg)
	end do
!-----------------------------------------------------
! ---optimize the CabsSeg segment for small NGauss
!-----------------------------------------------------
	IF(Ngauss.LT.Nseg) THEN
	  CALL  OptSegSLW(Tref,pref,yref,yref_all,CabsSeg,CabsAve,
     &				spcno,ncomp,NSEG,Ngauss,charlen)
	ELSE
	  DO Icomp=0,Ncomp
	  DO Igauss=1,Ngauss
	    CabsAve(Iseg,Icomp) = exp(0.5*( log(CabsSeg(Igauss))
     &		+log(CabsSeg(Igauss-1)) ) )
	  END DO
	  END DO
	END IF
	fave=0.d0
	do Igauss=1,Ngauss
	  aa1=0.0
	  ymol_all=0.0
	  do Icomp=1,Ncomp
	  if(spcno(Icomp).le.0.or.ParaFlag(spcno(Icomp)).EQ.0) CYCLE
	  if(yref(Icomp).lt.1.0d-10) CYCLE
	  CALL  GetDistriSLW(Tref,Tref,pref,yref(Icomp),
     &	        yref_all,CabsAve(Igauss,Icomp),spcno(Icomp),
     &		fave(Igauss,Icomp),almn(:,:,:,spcno(Icomp)),
     &		blmn(:,:,:,spcno(Icomp)))
	  aa1 = aa1+fave(Igauss,Icomp)*yref(Icomp)
	  ymol_all=ymol_all+yref(Icomp)
	  end do
	  fave(Igauss,0)=aa1/ymol_all
	end do
!
	if(CabsSeg(0).gt.0.d0) CabsLog(0)=log(CabsSeg(0))
	do Igauss=1,Ngauss
	  CabsLog(Igauss) = Log(CabsSeg(Igauss))
	  DO Icomp=0,Ncomp
	  IF(CabsAve(Igauss,Icomp).GT.0.0) 
     &		CabsAve(Igauss,Icomp) = Log(CabsAve(Igauss,Icomp))
	  END DO
	end do
!
	DO ICVL=ICVS,ICVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	 NOTE:
! The following code is according to the Ref.[5]	ASME-JHT, 117, 359-365,1995
! Any uncaustious changement may reduce the accuracy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ymol_t=0.d0
	  raddum(:)=0.d0
	  do 100 Icomp=1,ncomp
	  ymol=yys(ICVL,Icomp,1)*r_wm(Icomp)*rho(ICVL,1)
	  ymol_t=ymol_t+ymol
	  ffenbu(Icomp)=ymol
 100      end do
	  if(ymol_t.lt.1.0d-10) CYCLE
  !only for CO2/H2O mixture
         
	  CALL GetMixFunSLW(Tref,pref,yref,yref_all,tmp(ICVL,1),
     &		raddum,CabsSeg,ffenbu(1:ncomp),spcno,ncomp,
     &		NSEG,Ngauss)
	  do 102 Igauss=1,ngauss
	    ii0=NALLCV*(Igauss-1)
	    sumcoef(ii0+ICVL)=raddum(Igauss)
 102      end do
!
	  ymol_all = 0.d0
	  do Icomp=1,ncomp
	  ymol=yys(ICVL,Icomp,1)*r_wm(Icomp)*rho(ICVL,1)
	  if(spcno(Icomp).le.0.or.ParaFlag(spcno(Icomp)).EQ.0) CYCLE
	  if(ymol.lt.1.0d-10) CYCLE
	  ymol_all=ymol_all+ymol
	  GOTO 1050
	  ffenbu=0.d0	
	  do Igauss= 1,Ngauss			!0,Nseg
	  if(cminmax(1,spcno(Icomp)).gt.CabsSeg(Igauss)) cycle
	  if(cminmax(2,spcno(Icomp)).lt.CabsSeg(Igauss)) cycle
	  kosesb=0.0
!-------------------------
	  DO K=0,3		!l
          bjj=0.0
	  DO J=0,3		!m
	  aii=0.0	
	  DO I=0,2	!n
	  aii=aii+blmn(i,j,k,spcno(icomp))
     &   	  *((tmp(ICVL,1)*4.d-4)**I)
	  END DO
	  bjj=bjj+aii*(CabsLog(Igauss)**J)
	  END DO
	  kosesb=kosesb+bjj*(yref(Icomp)/yref_all)**(K+1)
	  END DO
!-------------------------
	  pfenbu=0.0
	  DO K=0,3		!l
	  bjj=0.d0
	  DO J=0,3		!m
	  aii=0.0
	  DO I=0,3	!n
	  aii=aii+almn(i,j,k,spcno(icomp))*((tref*4.d-4)**I)
	  END DO
	  bjj=bjj+aii*((tmp(ICVL,1)*4.d-4)**J)
	  END DO
	  pfenbu=pfenbu+bjj*(CabsLog(Igauss)-kosesb)**K
	  END DO
	  ffenbu(Igauss)=0.5*tanh(pfenbu)+0.5
	  end do
!-------------------------
	  ffenbu(Ngauss)=1.d0
1050	  do 1233 Igauss=1,ngauss
	  Cabs = CabsAve(Igauss,Icomp)		!log value
	  Cabs1= CabsLog(Igauss-1)-10.0
          Cabs2= CabsLog(Igauss)+10.0
          istep=0
1111	  istep=istep+1
          if(istep>10) goto 1234
	  kosesb=0.0					!a Bi-Section Loop
!-------------------------
          DO K=0,3		!l
	  bjj=0.0
          DO J=0,3		!m
	  aii=0.0	
	  DO I=0,2	!n
	  aii=aii+blmn(i,j,k,spcno(icomp))*((Tref*4.d-4)**I)
	  END DO
	  bjj=bjj+aii*(Cabs**J)
	  END DO
	  kosesb=kosesb+bjj*(ymol/ymol_t)**(K+1)		!NOTE:
	  END DO
	  pfenbu=0.0
!-------------------------
	  DO K=0,3		!l
	  bjj=0.d0
	  DO J=0,3		!m
	  aii=0.0
	  DO I=0,3	!n
	  aii=aii+almn(i,j,k,spcno(icomp))*
     &		 ((tmp(ICVL,1)*4.d-4)**I)
	  END DO
	  bjj=bjj+aii*((Tref*4.d-4)**J)
	  END DO
	  pfenbu=pfenbu+bjj*(Cabs-kosesb)**K
	  END DO
!-------------------------
	  f1=0.5*tanh(pfenbu)+0.5
	  aa1 = (f1-fave(Igauss,Icomp))/fave(Igauss,Icomp)
	  IF(aa1.GT.5.0d-3) THEN			!Cabs1---OK----Cabs------Cabs2
	    Cabs2= Cabs
	    Cabs=(Cabs+Cabs1)*0.5
	    GOTO 1111
	  ELSE IF(aa1.LT.-5.0d-3) THEN	!Cabs1----Cabs----OK-----Cabs2
	    Cabs1=Cabs
	    Cabs=(Cabs+Cabs2)*0.5
	    GOTO 1111
	  END IF
1234	  ii0 = NALLCV*(Igauss-1)
	  GWAbsorb(ii0+ICVL)=GWAbsorb(ii0+ICVL)+exp(Cabs)*ymol
 1233     end do
	  end do  !Icomp
	END DO		!ICVL
!----------------------------
! --- gray particle
!----------------------------  
	if(src_rad.ne.usryes) then
	  DO ICVL=ICVS,ICVE
	  PAbsorb(ICVL)=radmat(2,IIMAT)
	  PScatter(ICVL)=radmat(3,IIMAT)
	  END DO
	end if

	end do
!----------------------------
! --- boundary dummy cv
!----------------------------
 7777 do nb=1,nbcnd
	IIMAT=MAT_BCIDX(nb,1)
	if(.not.mat_cal(IIMAT)) cycle
	kd=kdbcnd(0,nb)
	kdt=kdbcnd(2,nb)
	kdy=kdbcnd(3,nb)
	IBFS=LBC_INDEX(nb-1)+1
	IBFE=LBC_INDEX(nb)
		
	if(kdt.eq.ktneum) then		!.or.kdt.eq.kttrns
          do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  aa1= 0.d0
	  DO Iseg=1,ngauss
		ii0=NALLCV*(Iseg-1)
		GWAbsorb(ii0+IDC) = radprop(1,nb)
		sumcoef(ii0+IDC) = sumcoef(ii0+ICV)		!??? OK
	        aa1 = aa1 + sumcoef(ii0+IDC)
          END DO
		  ii0=NALLCV*(Ngauss-1)
		  sumcoef(ii0+IDC) = sumcoef(ii0+IDC)+(1.0-aa1)
          end do
	else
!----------------------------------------------------------------		 
		 !at this version,
		 !we suppose only one fluid block, the tref,yref... are not changed
		 !new version for multi-fluids is under developing.
!----------------------------------------------------------------		 
	  ICFL=LBC_SSF(IBFS)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  tw=tmp(IDC,1)
	  ffenbu = 0.d0
	  if(yref_all.le.1.d-10) then
	    raddum(:)=1.d0/dble(Ngauss)
	  else
            ymol_all=0.d0
	    raddum(:)=0.d0
	    do Icomp=1,ncomp
            if(spcno(Icomp).le.0.or.ParaFlag(spcno(Icomp)).EQ.0) CYCLE
            ymol_all = ymol_all + yref(IComp)		!!!!!!!!!
            ffenbu = 0.d0
            do Igauss=1,ngauss
            if(cminmax(1,spcno(Icomp)).gt.CabsSeg(Igauss)) cycle
            if(cminmax(2,spcno(Icomp)).lt.CabsSeg(Igauss)) cycle
            kosesb=0.d0
            DO K=0,3
	    bjj=0.0
	    DO J=0,3		!m
	    aii=0.0
	    DO I=0,2	!n
	    aii=aii+blmn(i,j,k,spcno(icomp))*((tw*4.d-4)**I)
	    END DO
	    bjj=bjj+aii*(CabsLog(Igauss)**J)
            END DO
            kosesb=kosesb+bjj*(yref(Icomp)/yref_all)**(K+1)
	    END DO
	    pfenbu=0.0
	    DO K=0,3		!l
	    bjj=0.d0
	    DO J=0,3		!m
	    aii=0.0
	    DO I=0,3	!n
	    aii=aii+almn(i,j,k,spcno(icomp))*((tref*4.d-4)**I)
	    END DO
	    bjj=bjj+aii*((tw*4.d-4)**J)
            END DO
	    pfenbu=pfenbu+bjj*(CabsLog(Igauss)-kosesb)**K
	    END DO
	    ffenbu(Igauss)=0.5*tanh(pfenbu)+0.5
            end do		!Igauss
            ffenbu(Ngauss)=1.d0
	    do Igauss=1,ngauss
              raddum(Igauss)=raddum(Igauss)
     &	     +(ffenbu(Igauss)-ffenbu(Igauss-1))*yref(icomp)
	    end do
            enddo  !Icomp
!	
           aa1=0
           do Igauss=1,ngauss
	   raddum(igauss)=raddum(Igauss)/ymol_all		!NOTE HERE !!!
	   aa1=aa1+raddum(igauss)
	   end do
         end if
!
	 do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  aa1=0.0
	  DO Igauss=1,ngauss
	  ii0=NALLCV*(Igauss-1)
	  GWAbsorb(ii0+IDC) = radprop(1,nb)
	  sumcoef(ii0+IDC) = raddum(Igauss)			!??? OK,
	  END DO
         enddo
!
	end if
	end do
!
	RETURN


!-------> end of GetProSLW()

	END SUBROUTINE GetProSLW

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE  GetDistriSLW(Tb,Tg,P,Y,Yall,CabsSeg,spcno,ffenbu,
     &						 almn,blmn)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	implicit none
!
! --- [dummy arguments]
!	
	INTEGER ,intent(in)  ::	spcno
	REAL*8  ,intent(inout) :: almn(0:3,0:3,0:3), blmn(0:3,0:3,0:3)
        REAL*8  ,intent(inout) :: Tb,Tg,CabsSeg,Y,P,Yall,ffenbu
!
! --- [local entities]
!
        INTEGER :: i,j,k
	REAL*8	kosesb,sss,aii,bjj,pfenbu,CabsLog
	
	
	sss = 0.d0
	ffenbu = 0.d0
	CabsLog=log(CabsSeg)
			
c	if(cminmax(1,spcno(Icomp)).gt.CabsSeg(Iseg)) cycle
c	if(cminmax(2,spcno(Icomp)).lt.CabsSeg(Iseg)) cycle
	kosesb=0.d0
	DO 1120 K=0,3		!l
					   
		bjj=0.d0
		DO 1110 J=0,3		!m
						
			aii=0.d0	
			DO 1100 I=0,2	!n
				aii=aii+blmn(i,j,k)*((Tb*4.d-4)**I)
1100			CONTINUE

			bjj=bjj+aii*(CabsLog**J)
1110		CONTINUE

		kosesb=kosesb+bjj*((Y/Yall)**(K+1))
1120	CONTINUE
	
c	write(*,*) 'hh',spcno,kosesb,blmn(1,1,1),blmn(1,2,3)
c	pause

	pfenbu=0.d0
	DO 1220 K=0,3		!l
				  
		bjj=0.d0
		DO 1210 J=0,3		!m
		  
			aii=0.d0
			DO 1200 I=0,3	!n
				aii=aii+almn(i,j,k)*((Tg*4.d-4)**I)
1200			CONTINUE
					
			bjj=bjj+aii*((Tb*4.d-4)**J)
1210		CONTINUE
				  
		pfenbu=pfenbu+bjj*(CabsLog-kosesb)**K
1220	CONTINUE

	ffenbu=0.5*tanh(pfenbu)+0.d5
	

!-------> end of GetDistriSLW()

	END SUBROUTINE GetDistriSLW

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE  OptSegSLW(Tref,pref,yref,yref_all,CabsSeg,CabsAve,
     &				spcno,ncomp,NSEG,Ngauss,L)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	use module_rad, only : almn,blmn,ParaFlag

!	IMPLICIT REAL*8(A-H,O-Z)
	implicit none
        INTEGER ,intent(in)  ::	spcno(ncomp),ncomp,nseg,ngauss
        REAL*8  ,intent(inout) :: CabsSeg(0:NSeg),yref(ncomp),
     &           CabsAve(100,0:Ncomp)
        REAL*8  ,intent(inout) :: L,Tref,pref,yref_all
        

!	INTEGER ::	spcno(ncomp),ncomp,nseg,ngauss
!	REAL*8	::  CabsSeg(0:NSeg),yref(ncomp),CabsAve(100,0:Ncomp)
	REAL*8	::  Cabs,segf,ymol_all,kosesb,lenth,dff
	REAL*8	::  sumcoef(0:Nseg,0:Ncomp),ffenbu(0:Nseg),
     &			Absorb(0:Nseg),Emi(0:Nseg)
	INTEGER ::  ISEG,Igauss,I,J,K,Icomp,ii,iseg1
        REAL*8	::  bjj, aii,pfenbu,atotal,v1,v2,v3,vmin,dv,bemi
	
	if(yref_all.LE.1.d-10) return

	sumcoef = 0.d0
	Absorb= 0.d0
	ffenbu = 0.d0
	ymol_all = 0.d0

	do Icomp=1,ncomp
	  
	  if(spcno(Icomp).le.0.or.ParaFlag(spcno(Icomp)).EQ.0) cycle
	  if(yref(Icomp).lt.1.0d-10) cycle

	  ymol_all=ymol_all + yref(Icomp)

	  do Iseg=0, Nseg-1			!0,Nseg
	
		kosesb=0.d0
		
		DO 120 K=0,3		!l
		   bjj=0.d0
		   DO 110 J=0,3		!m
			aii=0.d0	
			DO 100 I=0,2	!n
			  aii=aii+blmn(i,j,k,spcno(icomp))*((Tref*4.0d-4)**I)
100			CONTINUE
		    bjj=bjj+aii*(log(CabsSeg(Iseg))**J)
110		   CONTINUE
		   kosesb=kosesb+bjj*(yref(Icomp)/yref_all)**(K+1)
120		CONTINUE
		
				
		pfenbu=0.d0
		DO 220 K=0,3		!l
		  bjj=0.d0
		  DO 210 J=0,3		!m
			aii=0.d0
			DO 200 I=0,3	!n
			 aii=aii+almn(i,j,k,spcno(icomp))*((Tref*4.d-4)**I)
200			CONTINUE
		    bjj=bjj+aii*((Tref*4.0d-4)**J)
210		  CONTINUE
		  pfenbu=pfenbu+bjj*(log(CabsSeg(Iseg))-kosesb)**K
220		CONTINUE

		ffenbu(Iseg)=0.5d0*tanh(pfenbu)+0.5d0

	!write(*,*) Iseg,CabsSeg(Iseg),kosesb,pfenbu,ffenbu(Iseg)
	
	  end do		!Iseg

	  ffenbu(Nseg)=1.d0
	  
	  sumcoef(0,Icomp)=ffenbu(0)
	  sumcoef(0,0)=sumcoef(0,0)+ffenbu(0)*yref(Icomp)

	  do Iseg=1,Nseg
	    sumcoef(Iseg,Icomp)=ffenbu(Iseg)-ffenbu(Iseg-1)
          sumcoef(Iseg,0)=sumcoef(Iseg,0)+sumcoef(Iseg,Icomp)
     &		*yref(Icomp)

	    Cabs = exp(0.5d0*(log(CabsSeg(Iseg))+log(CabsSeg(Iseg-1))))
	    Absorb(Iseg)=Absorb(Iseg)+yref(icomp)*Cabs
	  end do

	end do
	
		
	Emi = 0.d0
	atotal = 0
	sumcoef(0,0)=sumcoef(0,0)/ymol_all
	do Iseg=1,Nseg
		sumcoef(Iseg,0)=sumcoef(Iseg,0)/ymol_all
		atotal=atotal+sumcoef(Iseg,0)
!Emi(Iseg)=Emi(Iseg-1)+sumcoef(Iseg)*(1.d0-Exp(-Absorb(Iseg)*L))

c	write(*,*) 'old', Iseg,Emi(Iseg),Absorb(Iseg),sumcoef(Iseg),
c     &	L,Exp(-Absorb(Iseg)*L)

	end do
	
	
	!make segment
	dff=(1.d0-sumcoef(0,0))/Ngauss
	ffenbu = 0.d0
	Iseg1=0
	do Igauss=1,Ngauss
		v1=dff*Igauss+sumcoef(0,0)
		v2=0.d0
		vmin=100.d0
		do Iseg=0,Nseg
			v2 = v2+sumcoef(Iseg,0)
			dv = abs(v1-v2)
			if(dv.lt.vmin) then
				vmin=dv
				ffenbu(Igauss)=CabsSeg(Iseg)
				II=Iseg
				v3=v2
			end if
	
		end do

		v1=0.0
		v2=0.0
		bemi=0.0
		IF(IGauss.EQ.Ngauss) II=Nseg
		do Iseg=Iseg1+1,II
			v1 = v1 + Absorb(Iseg)*sumcoef(Iseg,0)
			v2 = v2 + sumcoef(Iseg,0)
		bemi=bemi+sumcoef(Iseg,0)*(1.d0-Exp(-Absorb(Iseg)*L))
		end do
		v3 = v1/v2
		Cabs=1.d0-(bemi/v2)
		CabsAve(Igauss,0)=-log(Cabs)/L/ymol_all

		DO Icomp=1,Ncomp
!------------
! --- ZD			
!------------
        if(spcno(Icomp).le.0.or.ParaFlag(spcno(Icomp)).EQ.0) cycle
			if(yref(Icomp).lt.1.0d-10) cycle

			v1=0.d0
			v2=0.d0
			bemi=0.d0
			IF(IGauss.EQ.Ngauss) II=Nseg
			do Iseg=Iseg1+1,II
			 v1 = v1 + Absorb(Iseg)*sumcoef(Iseg,Icomp)
			 v2 = v2 + sumcoef(Iseg,Icomp)
	 bemi=bemi+sumcoef(Iseg,Icomp)*(1.d0-Exp(-Absorb(Iseg)*L))
			end do
			v3 = v1/v2
			!vemi = v2*(1.d0-Exp(-v3*L))
			Cabs=1.d0-(bemi/v2)
			CabsAve(Igauss,Icomp)=-log(Cabs)/L/ymol_all

!write(*,*) 'opt',Igauss,II,CabsAve(Igauss,Icomp),ffenbu(Igauss)

		END DO

		Iseg1=II
	
	end do
	

	CabsSeg(1:Ngauss-1)=ffenbu(1:Ngauss-1)
	CabsSeg(Ngauss)=CabsSeg(Nseg)

	return

	!optimize
c	demi= Emi(Nseg)/Ngauss
c	Iseg1=0
c	emi0=0.0
c	ffenbu = 0.d0
c	ffenbu(0)=CabsSeg(0)
c	ffenbu(NGauss)=CabsSeg(Nseg)
c	do Igauss=1,Ngauss-1
c		bemi = demi*Igauss
c		vmin = 1.0e10*Emi(Nseg)
c		v1 = 0.0
c		v2 = 0.0
c		do Iseg=Iseg1+1,Nseg
c			v1 = v1 + Absorb(Iseg)*sumcoef(Iseg)
c			v2 = v2 + sumcoef(Iseg)
c			v3 = v1/v2
c			vemi = v2*(1.0-Exp(-v3*L))+emi0
c			dv = vemi-bemi
c			if(dv.gt.0) then
c			ffenbu(Igauss)=exp(0.5*(log(CabsSeg(Iseg))
c     &				+log(CabsSeg(Iseg-1))))
c				ii= iseg
c				emi0=vemi
c				Iseg1=ii
c				CabsAve(Igauss)= v3/ymol_all
c				exit
c			end if
c	
c		end do
	
c	write(*,*) 'Igauss=',Igauss,II,CabsAve(Igauss),v3,ymol_all
c
c	
c	end do

c	v1 = 0.d0
c	v2 = 0.d0	
c	do Iseg=Iseg1+1,Nseg
c		v1 = v1 + Absorb(Iseg)*sumcoef(Iseg)
c		v2 = v2 + sumcoef(Iseg)
c	end do
c	v3 = v1/v2
c	CabsAve(Ngauss)= V3/ymol_all
	
c	write(*,*) 'Igauss=',Ngauss,CabsAve(Igauss),v3,ymol_all

c	do Igauss=1,Ngauss
c		CabsSeg(Igauss)=ffenbu(Igauss)
c	end do
c	
c	
c	!pause
c
!-------> end of OptSegSLW()

	END SUBROUTINE OptSegSLW
		
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE  GetMixFunSLW(Tref,pref,yref,yref_all,tmp,sumcoef,
     &				CabsSeg,ymol,spcno,ncomp,NSEG,
     &				Ngauss)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	use module_rad, only : almn,blmn,ParaFlag
	implicit none

        INTEGER ,intent(in)    :: spcno(ncomp),ncomp,nseg,ngauss
        REAL*8  ,intent(inout) :: CabsSeg(0:Ngauss),yref(ncomp)
        REAL*8  ,intent(inout) :: ymol(ncomp)
        REAL*8  ,intent(inout) :: sumcoef(0:Ngauss)
        REAL*8  ,intent(inout) :: Tref,pref,yref_all,tmp

!	INTEGER ::	spcno(ncomp),ncomp,nseg,ngauss
!	REAL*8	::  CabsSeg(0:Ngauss),yref(ncomp),ymol(ncomp)
	REAL*8	::  Cabs,segf,ymol_all,kosesb,lenth
	REAL*8	::  ffenbu(0:Nseg),
     &			CabsLog(0:Nseg),CabsAve(Nseg)
	INTEGER ::  ISEG,Igauss,I,J,K,Icomp,NS,no(2),iseg1,iseg2
        REAL*8	::  absorb,rr,ftotal,pfenbu,bjj,aii,ff0,cmix,cmax
        REAL*8	::  ff1,cco2,ch2o,ffc,ff2,ff
!
	if(yref_all.LE.1.d-10) return

	sumcoef = 0.d0
	Absorb= 0.d0
	ffenbu = 0.d0
	ymol_all = 0.d0

	do Icomp=1,ncomp
	if(spcno(Icomp).eq.1) no(1)=icomp
	if(spcno(Icomp).eq.2) no(2)=icomp
	end do
	NS=max(1,Nseg/Ngauss)
	CabsLog(0)=log(CabsSeg(0))
!	rr=ymol(no(2))/(ymol(no(1))+ymol(no(2)))	!CO2
        rr=0.5d0
        
	ftotal=0.0
	Icomp=2		!CO2
	kosesb=0.d0
	pfenbu=0.d0
!----------------------------------------
	DO K=0,3		!l
        bjj=0.d0
	DO J=0,3		!m
	aii=0.d0
	DO I=0,3	!n
	aii=aii+almn(i,j,k,spcno(icomp))*((tref*4.d-4)**I)
	END DO
	bjj=bjj+aii*((tmp*4.d-4)**J)
	END DO
	pfenbu=pfenbu+bjj*(CabsLog(0)-kosesb)**K
	END DO
	ff0=0.5*tanh(pfenbu)+0.5d0
!
	do Igauss= 1,Ngauss
	Iseg1= (Igauss-1)*NS
	Iseg2= Igauss*NS
	lenth=(log(CabsSeg(Igauss))-log(CabsSeg(Igauss-1)))/NS
	do Iseg=Iseg1+1,Iseg2
	CabsLog(Iseg)=CabsLog(Iseg-1)+lenth
	CabsAve(Iseg)=Exp(CabsLog(Iseg-1)+lenth*0.5d0)
        end do
	end do

	do Igauss= 1,Ngauss
	Cmix=CabsSeg(Igauss)
	Cmax=(Cmix-(1.d0-rr)*CabsSeg(0))/rr
	!Cmax=(CabsSeg(Ngauss)-(1.d0-rr)*CabsSeg(0))/rr
	lenth=(log(Cmax)-log(CabsSeg(0)))/NSeg
	ftotal = 0.d0
	ff1=ff0
	do Iseg=1,Nseg
	Cco2=log(CabsSeg(0))+lenth*Iseg
	Ch2o=(Cmix-rr*exp(Cco2))/(1.d0-rr)
	Ch2o=log(Ch2o)
	icomp=1   !no(1)		!H2O
	kosesb=0.d0
	DO K=0,3		!l
	bjj=0.d0
	DO J=0,3		!m
	aii=0.0	
	DO I=0,2	!n
	  aii=aii+blmn(i,j,k,spcno(icomp))*((tmp*4.d-4)**I)
	END DO
	bjj=bjj+aii*(Ch2o**J)
        END DO
        kosesb=kosesb+bjj*(yref(Icomp)/yref_all)**(K+1)
	END DO
!		    
	pfenbu=0.d0
	DO K=0,3		!l
	  bjj=0.d0
	  DO J=0,3		!m
	  aii=0.d0
	  DO I=0,3	!n
	  aii=aii+almn(i,j,k,spcno(icomp))*((tref*4.d-4)**I)
	  END DO
	  bjj=bjj+aii*((tmp*4.d-4)**J)
	  END DO
	  pfenbu=pfenbu+bjj*(Ch2o-kosesb)**K
	END DO
        

	ffc=0.5d0*tanh(pfenbu)+0.5d0

	Icomp=no(2)		!CO2
	kosesb=0.0
	pfenbu=0.0
	DO K=0,3		!l
	  bjj=0.d0
	  DO J=0,3		!m
	  aii=0.d0
	  DO I=0,3	!n
	  aii=aii+almn(i,j,k,spcno(icomp))*((tref*4.d-4)**I)
	  END DO
	  bjj=bjj+aii*((tmp*4.d-4)**J)
	  END DO
	  pfenbu=pfenbu+bjj*(Cco2-kosesb)**K
	END DO
	ff2=0.5d0*tanh(pfenbu)+0.5d0
	ff=(ff2-ff1)*ffc
	ftotal=ftotal+ff
	ff1=ff2
        end do
	sumcoef(Igauss)=ftotal
        end do
	do Igauss=Ngauss,1,-1
	sumcoef(Igauss)=sumcoef(Igauss)-sumcoef(Igauss-1)
	end do

	END SUBROUTINE GetMixFunSLW
!


