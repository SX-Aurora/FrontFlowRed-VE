!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine vel_admin(
     &  iphs,iter,ismpl,deltt,nsmpl,rsdc,time,ictl,volctr,
     &  LVEDGE,LCYCSF,LFUTAU,LBC_SSF,LCYCOLD,wifsld,OPPANG,
     &  SFAREA,SFCENT,CVCENT,CVVOLM,CVVOL0,wiface,FRSTCV,DISALL,
     &  velopp,rhoopp,prsopp,aks,fric_c,dvddopp,MASS_I,SDOT,xta,pp0,
     &  rmu,rmut,tmp,rva,vel,rho,prs,yys,utau,dvdd,
     &  grdc,diag,kdbv,rvd,velf,
     &  rmue,cofd,dsclv,dvdt,bdyfrc,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  FIELD_U,MOMCELL,
     &  angbnd,t_dir,gggg,
     &  iterv,repsv,aepsv,errv,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-------------------------------------------------
!     rmue  <=   rmx 
!     cofd  <=   cp 
!     dsclv <=   cr 
!-------------------------------------------------
!
! --- [module arguments] 
!
      use module_dimension
      use module_constant
      use module_hpcutil
!
      use module_io,only       : ifle
      use module_model,only    : idrdp,incomp,mach0,comp,ical_vect,
     &                           nthrds,iLBF_P,
     &                           PEFC,PAFC,MCFC,SOFC,AFC,MDFC,
     &                           ical_prt,u_func,drag_coef
      use module_flags,only    : 
     &                           intgvv,euleri,eulere,Adams_Bashforth,
     &                           Adams_Moulton,Crank_Nicolson,
     &                           Runge_Kutta,
     &                           icalrv,face
      use module_material,only : nflud,nofld,lclsd,nsold,nosld,
     &                           ivscvv,cal,icnvvv,lmtrvv,
     &                         ical_sld,rotati,ishaft,end,begin,
     &                           rot_ang,
     &                           porous,N_pors,C0,C1,ical_porous,
     &                           C2p,alpha_co,ino,ipower,idarcy,
     &                           idarcy2,
     &                           porosty,relaxv,domegadt
      use module_model,only    : ICVREF,PREFM,ical_MHD,ical_dens
      use module_usersub,ONLY  : src_uvw,src_r,usrno,usryes,src_fire
      use module_boundary,only : kdbcnd,kdilet,kvlglw,kvnslp,kxnone,
     &                           kdintr,boundName,kdfire,kdbuff,kdsld,
     &                           kdolet,kdprdc,kdsymm,kdtchi,kdcvd,
     &                           distrb,pkdwall,nobcnd,
     &                           phs_idx,phs_com,
     &                           LBC_INDEX,nbcnd,MAT_BCIDX,openout,
     &                           imasflg,masbc,dvel,surfreac,
     &                           ical_poroBC,poro,alpha_coBC,C2pBC,
     &                           thkposBC,prosty
      use module_dimnsn,only   : xredc,yredc,zredc
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_gravity ,only : ggg,ramb,tamb,beta,beta2,buoyancy,
     &                           ramb2,tamb2,iave_rho_T
      use module_vof     ,only : ical_vof,intervof,LG,LS,change_phase,
     &                           grvty,byncy,modul
      use module_initial ,only : rho0,rho02
      use module_chemreac,ONLY : ical_suf
      use module_scalar,  only : ivof,iaph,ical_FC,ical_s,ical_cavi
      use module_metrix,only   : vof=>W1K11
      use module_metrix,only   : rcomp,sinj_comp,eva_comp
      use module_species, only : wm
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_metrix  ,only : vsm,rsm,tsm,rsm_prim
!      use module_material,only : KMAT_S,MAT_S
      use module_time,    only : i_steady
      use module_particle,only : latent
      use module_metrix,only  : bdyf,vctr,multiR
      use module_scalar,only  : Axis_DIR,FSTR_F,CAVI_B,PWR
      USE module_usersub,ONLY : src_BC_CAVI
!
! 1.  Update velocity field
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,ismpl,nsmpl,iphs
      real*8 ,intent(inout) :: deltt,rsdc,time,volctr
      integer,intent(inout) :: ictl
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO(      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV(      MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)      
      integer,intent(in)    :: LFUTAU    (     MXCV)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (  MXALLCV)
      real*8 ,intent(in)    :: wiface    (  MXCVFAC)
      real*8 ,intent(in)    :: FRSTCV    (  MXSSFBC)
      real*8 ,intent(in)    :: DISALL    (  MXCV)
      real*8 ,intent(in)    :: xta       (  MXCVFAC)
      real*8 ,intent(inout) :: rmu       (  MXALLCV)
      real*8 ,intent(inout) :: rmut      (  MXALLCV)
      real*8 ,intent(in)    :: tmp       (  MXALLCV)
      real*8 ,intent(inout) :: rva       (  MXCVFAC,2)
      real*8 ,intent(inout) :: vel       (  MXALLCV,3,2)
      real*8 ,intent(inout) :: velf      (  MXCVFAC_B,3,2)
      real*8 ,intent(in)    :: yys       (  MXALLCV,MXcomp)
      real*8 ,intent(inout) :: rho       (  MXALLCV,2)
      real*8 ,intent(inout) :: prs       (  MXALLCV,2)
      real*8 ,intent(inout) :: pp0       (  MXALLCV)
      real*8 ,intent(inout) :: grdc      (  MXALLCV,3,3)
      real*8 ,intent(inout) :: diag      (  MXALLCV)
      real*8 ,intent(inout) :: dvdd      (     MXCV,3,2)
      real*8 ,intent(inout) :: dvdt      (  MXALLCV,3)
      REAL*8 ,INTENT(INOUT) :: rvd       (  MXCVFAC,3)
      REAL*8 ,INTENT(INOUT) :: cofd      (  MXALLCV)
      REAL*8 ,INTENT(INOUT) :: dsclv     (  MXALLCV)
      REAL*8 ,INTENT(INOUT) :: rmue      (  MXALLCV)
      real*8 ,intent(in)    :: utau      (0:MXSSFBC)
      real*8 ,intent(inout) :: aks       (  MXALLCVR,mxrans,2)
      REAL*8 ,INTENT(IN)    :: SDOT      (MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,intent(in)    :: velopp
     &        ((MXALLCV*(iphs-1)+MXALLCV2*(2-iphs)),3,2)
      REAL*8 ,intent(in)    :: rhoopp
     &        ((MXALLCV*(iphs-1)+MXALLCV2*(2-iphs)),2)
      REAL*8 ,intent(in)    :: prsopp
     &        ((MXALLCV*(iphs-1)+MXALLCV2*(2-iphs)),2)
      real*8 ,intent(in)    :: fric_c    ( MXALLCV2)
!
      REAL*8 ,intent(in)    :: dvddopp
     &        ((MXCV*(iphs-1)+MXCV2*(2-iphs)),3,2)
      INTEGER,INTENT(IN)    :: kdbv       (MXCVFAC)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
      
      integer,intent(out)   :: ierror,iterv(3)
      real*8 ,intent(out)   :: repsv(3),aepsv(3),errv(3)
!      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(inout) :: bdyfrc(MXCV_B,3,2)
      real*8 ,intent(in)    :: MASS_I     (  MXALLCV2)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
      REAL*8 ,INTENT(inout)    :: MOMCELL(MXALLCV_P,4)
!
      real*8 ,intent(inout) :: angbnd(MXSSFBC_VOF)
      real*8,intent(in)     :: t_dir(MXALLCV_VOF,3)
      real*8 ,intent(in)    :: gggg
!
! --- [local entities]
!
      real*8 ,parameter :: r6=1.d0/6.d0,r2=0.5D0,SML1=1.d-15
!
      real*8  :: duma(1),rr,vol,ad_comp
      real*8  :: sinj_r,sinj_u,sinj_v,sinj_w,sinj_t
      real*8  :: dxx,dyx,dzx,dll
      real*8  :: u,v,w,dens,suu,suv,suw,utauu,p,T
      real*8  :: wi1,wi2,grx,gry,grz,gf0,pref,rmax,rmax1,rmax2
      real*8  :: dx,dy,dz,dpAB,dl,dl2,dvu,dvv,dvw,UCONV,grdsf
      real*8  :: ru,rv,rw,errabs(3),errabp
      real*8  :: dum1,dum2,dum3,dum4,dum5,dum6,
     &           rdeltt,timtrm,ER,ERP,ALPHP,delt
      real*8  :: rvaicf,rvaicfp,calph
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,nx,ierr1,ndf,nq
      integer :: ICOM,IMD,ICH,IFLD,ICTP
      integer :: IMAT,IIMAT,IMAT_U,ICVS,ICVE,ICVL,IWL
      integer :: ICFS,ICFE,ICFL
      integer :: IO,IN,IP,IS,IW
      integer :: IC,IV,IE,ICF,ICV,IBF,NB,KD,IDC,ICVFIX
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ICVLA,ICVLB,IBFS,IBFE,IMODE,IDCS,IDCE,IBFL
      integer :: ipcom,icoms,icome,iph
      real*8  :: gf1,dlvect,tha,r(3),vr(3),betax,tambx
      real*8  :: velref(3),shaftn(3),radius(3),bb(3,3),rbb(3,3)
      integer :: ICODE,myid,no
      logical :: explicit,viscous,noviscs
      character*80 :: BC_name
      real*8  :: unitmass,T_WALL,eva_heat,eva_T,eva_mass,FMASS,dum_typ
      real*8  :: dum_grv,rambx,rambx_op,grvi,porodum
      integer :: nb2,IBFS1,IBFE1,iphs_1,iphs_2,KMAT
      integer :: ndiag,idum
      real(8) :: rotup_inertial(3)
!
!------------
! ---
!------------
!
!      ALLOCATE (dvdt(MXALLCV,3) ,stat=ierr1)
!
! --- < 1.1 Preliminary set >-
      rdeltt=1.d0/deltt
      if(i_steady==3) then
        rdeltt=1.d0      ! vel_admin=0.d0   !1212
      endif

      iphs_1=iphs
      iphs_2=3-iphs
!
      explicit=intgvv.eq.eulere.OR.intgvv.eq.Adams_Bashforth.or.
     &         intgvv.eq.Runge_Kutta
!
      if(intgvv.eq.Adams_Bashforth) then
        calph=0.5d0
        if( iter.eq.1 ) calph=0.d0
      elseif(intgvv.eq.Adams_Moulton) then
        calph=1.d0/12.d0
        if( iter.eq.1 ) calph=0.d0
      elseif(intgvv.eq.Crank_Nicolson) then
        calph=0.5d0
        if( iter.eq.1 ) calph=0.d0
      else
        calph=0.d0
      endif
!
      if(ical_MHD==0) dvdt(:,1:3)=0.0D0   !dvdt=>[N/m^3]
      diag=0.D0
      rmue=0.D0
      dsclv=-1.d0
!
      ierror=0
      ICODE=0
!
!--------------------------
! --- Gravity term: 
!--------------------------
!
      if(buoyancy.eq.1.and..not.
     &  (ical_vof==1.or.ieul2ph>0.or.ical_FC==PEFC)) then
        call gravity(dvdt(:,1:3))
        if(iLBF_P==5) then
          do iv=1,3
          bdyfrc(:,iv,1)=dvdt(:,iv)
          enddo
          dvdt(:,:)=0.d0
        endif
      endif
!------------------------
! --- large body-force --
!------------------------
!-----------------------------------------------------
! --- User Source term: Not including density source
!-----------------------------------------------------
      if(src_uvw.eq.usryes) then
        do 150 IIMAT=1,NMAT   !ICN=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT.gt.0) then
          IMAT_U=nofld(IMAT)
          
          DO 151 ICVL=ICVS,ICVE
          ICV=MAT_CV(ICVL)
          if(NPE.gt.1) ICV=NOD_IDg(ICV)
          IBFL=LFUTAU(ICVL)
          utauu=utau(IBFL)
          u=vel(ICVL,1,1)
          v=vel(ICVL,2,1)
          w=vel(ICVL,3,1)
          t=tmp(ICVL)
          dens=rho(ICVL,1)
          suu=0.d0
          suv=0.d0
          suw=0.d0
          call user_src_uvw
     &       (deltt,iter,time,ICV,IMAT_U,iphs,u,v,w,t,dens,
     &                    suu,suv,suw,utauu)
          dvdt(ICVL,1)=dvdt(ICVL,1)+suu
          dvdt(ICVL,2)=dvdt(ICVL,2)+suv
          dvdt(ICVL,3)=dvdt(ICVL,3)+suw
 151      continue
        elseif(IMAT.lt.0) then
          IMAT_U=nosld(-IMAT)
        endif
 150    continue
      endif
!--------------------------
! --- User rho soure term
!--------------------------
      if(src_r.eq.usryes) then
        do 155 IIMAT=1,NMAT   !ICN=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT.gt.0) then
          IMAT_U=nofld(IMAT)
          DO 156 ICVL=ICVS,ICVE
          ICV=MAT_CV(ICVL)
          if(NPE.gt.1) ICV=NOD_IDg(ICV)
          u=vel(ICVL,1,1)
          v=vel(ICVL,2,1)
          w=vel(ICVL,3,1)
          t=tmp(ICVL)
          p=prs(ICVL,1) 
          dens=rho(ICVL,1)
          sinj_u=0.d0
          sinj_v=0.d0
          sinj_w=0.d0
          sinj_r=0.d0
          sinj_t=tmp(ICVL)
          sinj_comp(:)=0.d0
          vol=CVVOLM(ICVL)
          do icom=1,ncomp
            rcomp(icom)=yys(ICVL,icom)
          enddo
          dxx=CVCENT(1,ICVL)
          dyx=CVCENT(2,ICVL)
          dzx=CVCENT(3,ICVL)
          call user_src_r(3,ictl,
     &    deltt,iter,time,ICV,IMAT_U,iphs,ncomp,dxx,dyx,dzx,
     &    u,v,w,T,p,dens,rcomp,vol,
     &         sinj_r,sinj_u,sinj_v,sinj_w,sinj_t,sinj_comp)
          suu=0.d0
          suv=0.d0
          suw=0.d0
          if(sinj_r.gt.SML) then 
            suu=sinj_r*(sinj_u-0.d0*u)  ! sinj_r=[kg/(s*m^3)]
            suv=sinj_r*(sinj_v-0.d0*v)  ! zhang-cvd
            suw=sinj_r*(sinj_w-0.d0*w)
            diag(ICVL)=diag(ICVL)+sinj_r
!
!            dvdt(ICVL,1)=dvdt(ICVL,1)+suu
!            dvdt(ICVL,2)=dvdt(ICVL,2)+suv
!            dvdt(ICVL,3)=dvdt(ICVL,3)+suw
! ---- removed gravity
            dvdt(ICVL,1)=suu
            dvdt(ICVL,2)=suv
            dvdt(ICVL,3)=suw
          elseif(sinj_r.lt.-SML) then !NOT finished:
            suu=sinj_r*vel(ICVL,1,1)
            suv=sinj_r*vel(ICVL,2,1)
            suw=sinj_r*vel(ICVL,3,1)
            diag(ICVL)=diag(ICVL)-sinj_r
            dvdt(ICVL,1)=dvdt(ICVL,1)+suu
            dvdt(ICVL,2)=dvdt(ICVL,2)+suv
            dvdt(ICVL,3)=dvdt(ICVL,3)+suw
          endif
 156      continue
        elseif(IMAT.lt.0) then
          IMAT_U=nosld(-IMAT)
        endif
 155  continue
      endif
!-------------------------
! --- canopy model
!-------------------------
      if(u_func(3)==1) then
        do IIMAT=1,NMAT   !ICN=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT.gt.0) then
          IMAT_U=nofld(IMAT)
          DO ICVL=ICVS,ICVE
            idum=multiR(ICVL)
            if(idum>=21.and.idum<=50) then

            u=vel(ICVL,1,1)
            v=vel(ICVL,2,1)
            w=vel(ICVL,3,1)
            dens=rho(ICVL,1)
            dum1=dens*drag_coef(idum)*sqrt(u**2+v**2+w**2)
            suu=dum1*u
            suv=dum1*v
            suw=dum1*w
            dvdt(ICVL,1)=dvdt(ICVL,1)-suu
            dvdt(ICVL,2)=dvdt(ICVL,2)-suv
            dvdt(ICVL,3)=dvdt(ICVL,3)-suw
            diag(ICVL)=dum1
            endif
            enddo

        elseif(IMAT.lt.0) then
          IMAT_U=nosld(-IMAT)
        endif
        enddo

      endif
!---------------------------------------
! --- Fire wall evaporation source term
!---------------------------------------
      if(src_fire.eq.usryes) then
        do nb=1,nbcnd
        BC_name=boundName(nb,1)
        IIMAT=MAT_BCIDX(nb,1)
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT.lt.0) cycle
        IMAT_U=nofld(IMAT)
        kd=kdbcnd(0,nb)
        no=nobcnd(nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd.eq.kdfire) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL) 
          ICVL=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          u=vel(ICVL,1,1)
          v=vel(ICVL,2,1)
          w=vel(ICVL,3,1)
          T=tmp(ICVL)
          T_WALL=tmp(IDC)
          if(idrdp.eq.mach0.or.idrdp.eq.comp) then
            p=prs(ICVL,1)
          else
            p=prs(ICVL,1)+pp0(ICVL)
          endif
          dens=rho(ICVL,1)
          do icom=1,ncomp
            rcomp(icom)=yys(ICVL,icom)
          enddo
          eva_mass=0.d0
          eva_comp(:)=0.d0
          eva_heat=0.d0
          eva_T=tmp(ICVL)
!
          call user_src_fire
     &   (no,BC_name,deltt,iter,time,IMAT_U,iphs,ncomp,
     &    u,v,w,p,dens,rcomp,T_WALL,
     &    T,eva_T,eva_comp,eva_mass,eva_heat)
!
          if(eva_mass.gt.0.d0) then  ! unitmass=[kg/(s*m^3)]
            unitmass=eva_mass*SFAREA(4,ICFL)/CVVOLM(ICVL)
            suu=unitmass*vel(ICVL,1,1)
            suv=unitmass*vel(ICVL,2,1)
            suw=unitmass*vel(ICVL,3,1)
            dvdt(ICVL,1)=dvdt(ICVL,1)+suu
            dvdt(ICVL,2)=dvdt(ICVL,2)+suv
            dvdt(ICVL,3)=dvdt(ICVL,3)+suw
          else
          endif
          enddo
        endif
        enddo
      endif
!------------------------------------------
! --- surface reaction (sticking function)
!------------------------------------------
      if(ical_suf==1) then
	do nb=1,nbcnd
        BC_name=boundName(nb,1)
        IIMAT=MAT_BCIDX(nb,1)
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle
        IMAT_U=nofld(IMAT)
        kd=kdbcnd(0,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(surfreac(nb)>0) then
          iph=1
          icoms=phs_idx(iph-1)+1
          icome=phs_idx(iph)
          do ipcom=icoms,icome
            icom=phs_com(ipcom)
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICVL=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            vol=CVVOLM(ICVL)
	    ad_comp=SDOT(IBFL,icom)*SFAREA(4,ICFL)*wm(icom)/vol
            suu=ad_comp*vel(ICVL,1,1)
            suv=ad_comp*vel(ICVL,2,1)
            suw=ad_comp*vel(ICVL,3,1)
            dvdt(ICVL,1)=dvdt(ICVL,1)+suu
            dvdt(ICVL,2)=dvdt(ICVL,2)+suv
            dvdt(ICVL,3)=dvdt(ICVL,3)+suw
            if(ad_comp>SML) then
	    else
              diag(ICVL)=diag(ICVL)-ad_comp
            endif
            enddo
          enddo
        endif
        enddo
      endif
!-------------------------------------------------------------------
! --- VOF (buoyancy force, BUT not consider phase change) zhang123 
!-------------------------------------------------------------------
      if(ical_vof==1) then  
      endif
!
!-----------------------------
! --- Sliding Frame (rotation)
!---------------------  zhang8
      if(ical_sld/=0) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.le.0) cycle
        if(ishaft(IMAT)==1) then 
          ICVS=MAT_CVEXT(IIMAT-1)+1 
          ICVE=MAT_CVEXT(IIMAT)
          shaftn(1)=end(1,IMAT)-begin(1,IMAT) 
          shaftn(2)=end(2,IMAT)-begin(2,IMAT) 
          shaftn(3)=end(3,IMAT)-begin(3,IMAT) 
          dum1=dsqrt(shaftn(1)**2+shaftn(2)**2+shaftn(3)**2) 
          shaftn(1)=shaftn(1)/dum1
          shaftn(2)=shaftn(2)/dum1
          shaftn(3)=shaftn(3)/dum1
          DO ICVL=ICVS,ICVE 
          radius(1)=CVCENT(1,ICVL)-begin(1,IMAT)   !IMAT??? 
          radius(2)=CVCENT(2,ICVL)-begin(2,IMAT)
          radius(3)=CVCENT(3,ICVL)-begin(3,IMAT)
! ----------------------------------------
          rotup_inertial(3)=  -rho(ICVL,1)*domegadt(IMAT)*
     &                (shaftn(1)*radius(2)-shaftn(2)*radius(1))
          rotup_inertial(2)=  -rho(ICVL,1)*domegadt(IMAT)*
     &                (shaftn(3)*radius(1)-shaftn(1)*radius(3))
          rotup_inertial(1)=  -rho(ICVL,1)*domegadt(IMAT)*
     &                (shaftn(2)*radius(3)-shaftn(3)*radius(2))
!----------------------------------------
          dum3=radius(1)*shaftn(1)
     &        +radius(2)*shaftn(2)
     &        +radius(3)*shaftn(3)
          radius(1)=radius(1)-dum3*shaftn(1)
          radius(2)=radius(2)-dum3*shaftn(2)
          radius(3)=radius(3)-dum3*shaftn(3)
          vr(1)=vel(ICVL,1,1)
          vr(2)=vel(ICVL,2,1)
          vr(3)=vel(ICVL,3,1)
!!!!          call cal_AXBEQC(shaftn,vr,velref,dum2) 
!
          velref(3)=(shaftn(1)*vr(2)-shaftn(2)*vr(1))
          velref(2)=(shaftn(3)*vr(1)-shaftn(1)*vr(3))
          velref(1)=(shaftn(2)*vr(3)-shaftn(3)*vr(2))
!
          dum1=-2.d0*rho(ICVL,1)*velref(1)*rotati(IMAT)
     &              +rho(ICVL,1)*radius(1)*rotati(IMAT)**2
     &              +rotup_inertial(1)
          dum2=-2.d0*rho(ICVL,1)*velref(2)*rotati(IMAT)
     &              +rho(ICVL,1)*radius(2)*rotati(IMAT)**2
     &              +rotup_inertial(2)
          dum3=-2.d0*rho(ICVL,1)*velref(3)*rotati(IMAT)
     &              +rho(ICVL,1)*radius(3)*rotati(IMAT)**2
     &              +rotup_inertial(3)
          dvdt(ICVL,1)=dvdt(ICVL,1)+dum1
          dvdt(ICVL,2)=dvdt(ICVL,2)+dum2
          dvdt(ICVL,3)=dvdt(ICVL,3)+dum3
          enddo
        endif
        enddo
      endif
!--------------------
! --- porous media
!--------------------
      if(ical_porous==1) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.le.0) cycle
        IMAT_U=nofld(IMAT)
        if(IMAT_U>N_pors) then
          if(porous(IMAT)==ipower) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            DO ICVL=ICVS,ICVE
            dum4=dsqrt(vel(ICVL,1,1)**2
     &                +vel(ICVL,2,1)**2
     &                +vel(ICVL,3,1)**2)**(C1(IMAT)-1.D0)
            dum1=-C0(IMAT)*vel(ICVL,1,1)*dum4
            dum2=-C0(IMAT)*vel(ICVL,2,1)*dum4
            dum3=-C0(IMAT)*vel(ICVL,3,1)*dum4
            DUM5=-C0(IMAT)*dum4
            dvdt(ICVL,1)=dvdt(ICVL,1)+dum1
            dvdt(ICVL,2)=dvdt(ICVL,2)+dum2
            dvdt(ICVL,3)=dvdt(ICVL,3)+dum3
            diag(ICVL)=diag(ICVL)-DUM5
            enddo
          elseif(porous(IMAT)==idarcy) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            DO ICVL=ICVS,ICVE
            dum4=dsqrt(vel(ICVL,1,1)**2
     &                +vel(ICVL,2,1)**2
     &                +vel(ICVL,3,1)**2)
            DUM5=rmu(ICVL)/alpha_co(IMAT)+rho(ICVL,1)*C2p(IMAT)*dum4
            dum1=vel(ICVL,1,1)*DUM5
            dum2=vel(ICVL,2,1)*DUM5
            dum3=vel(ICVL,3,1)*DUM5
!
            dvdt(ICVL,1)=dvdt(ICVL,1)-dum1
            dvdt(ICVL,2)=dvdt(ICVL,2)-dum2
            dvdt(ICVL,3)=dvdt(ICVL,3)-dum3
            diag(ICVL)=diag(ICVL)+DUM5
            enddo
          elseif(porous(IMAT)==idarcy2) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            DO ICVL=ICVS,ICVE
            DUM5=rmu(ICVL)/alpha_co(IMAT)*porosty(IMAT)
            dum1=vel(ICVL,1,1)*DUM5
            dum2=vel(ICVL,2,1)*DUM5
            dum3=vel(ICVL,3,1)*DUM5
!
            dvdt(ICVL,1)=dvdt(ICVL,1)-dum1
            dvdt(ICVL,2)=dvdt(ICVL,2)-dum2
            dvdt(ICVL,3)=dvdt(ICVL,3)-dum3
            diag(ICVL)=diag(ICVL)+DUM5
            enddo
          endif
        endif
        enddo
      endif
!
! --- 
!
      if(ical_poroBC==1) then 
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(poro(nb)==1) then
          
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL) 
          ICFP=LCYCSF(IBFL)
          ICVL=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDC=LVEDGE(2,ICFL)
          dum4=(vel(ICVL,1,1)*SFAREA(1,ICFL)
     &         +vel(ICVL,2,1)*SFAREA(2,ICFL)
     &         +vel(ICVL,3,1)*SFAREA(3,ICFL))
          dum3=(vel(ICVP,1,1)*SFAREA(1,ICFP)
     &         +vel(ICVP,2,1)*SFAREA(2,ICFP)
     &         +vel(ICVP,3,1)*SFAREA(3,ICFP))
          if(dum4>dum3) then
            ICF=ICFL
            ICV=ICVL
            IC=IDC
          else
            ICF=ICFP
            ICV=ICVP
            IC=IDCP
          endif
!          dum1=CVVOLM(ICVL)+CVVOLM(ICVP)
          dum1=CVVOLM(ICV)
          gf1=SFAREA(4,ICF)*thkposBC(nb)/dum1
          dum4=sqrt(vel(ICV,1,1)**2    !*SFAREA(1,ICF)
     &             +vel(ICV,2,1)**2    !*SFAREA(2,ICF)
     &             +vel(ICV,3,1)**2)    !*SFAREA(3,ICF))
!          u=dum4*SFAREA(1,ICF)
!          v=dum4*SFAREA(2,ICF)
!          w=dum4*SFAREA(3,ICF)
!          dum6=dsqrt(u**2
!     &              +v**2
!     &              +w**2)
          DUM5=(rmu(ICV)/alpha_coBC(nb)+0.5d0*rho(ICV,1)*C2pBC(nb) 
     &          *dum4)*gf1
          dum1=DUM5*vel(ICV,1,1)
          dum2=DUM5*vel(ICV,2,1)
          dum3=DUM5*vel(ICV,3,1)
          dvdt(ICV,1)=dvdt(ICV,1)-dum1
          dvdt(ICV,2)=dvdt(ICV,2)-dum2
          dvdt(ICV,3)=dvdt(ICV,3)-dum3
          diag(ICV)=diag(ICV)+DUM5
!
!          dum1=DUM5*vel(ICVP,1,1)
!          dum2=DUM5*vel(ICVP,2,1)
!          dum3=DUM5*vel(ICVP,3,1)
!          dvdt(ICVP,1)=dvdt(ICVP,1)-dum1
!          dvdt(ICVP,2)=dvdt(ICVP,2)-dum2
!          dvdt(ICVP,3)=dvdt(ICVP,3)-dum3
!          diag(ICVP)=diag(ICVP)+DUM5
!
          enddo
        endif
        enddo
      endif
!---------------------------------------
! --- 
!---------------------------------------
      if(ical_cavi==1) then 
      endif
!---------------------------
! --- body force save 
!---------------------------3333
      if(iLBF_P==1.or.iLBF_P==2) then
        if(ieul2ph>0) then
        else
          do iv=1,3
          bdyfrc(:,iv,1)=dvdt(:,iv)
          enddo
        endif
      endif
!
!      do iv=1,3
!      bdyf(:,iv)=dvdt(:,iv)
!      enddo
!
!---------------------------------------------------- 
!     Moment source of PARTICLE (TWO-WAY METHOD) 
!---------------------------------------------------- 
      IF(ismpl.EQ.1.and.ical_prt>=1) then !.and.iter>=ical_prt) THEN
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICVL=ICVS,ICVE
          FMASS=CVVOLM(ICVL)!*rho(ICVL,2)
          dvdt(ICVL,1)=dvdt(ICVL,1)+MOMCELL(ICVL,1)/FMASS
          dvdt(ICVL,2)=dvdt(ICVL,2)+MOMCELL(ICVL,2)/FMASS
          dvdt(ICVL,3)=dvdt(ICVL,3)+MOMCELL(ICVL,3)/FMASS
          diag(ICVL)=diag(ICVL)+MOMCELL(ICVL,4)/FMASS
        ENDDO
        ENDDO
      ENDIF
!
!------------------------------------------------------------------
! --- < 2.1 set time difference, gravity, pressure &source term >--
!------------------------------------------------------------------
!
!-------------------
! --- source term: 
!-------------------
!
      if(ical_vect) then  !chang1
        DO ICVL=ICVS_V,ICVE_V    !index_c(myid)+1,index_c(myid+1)
          diag(ICVL)=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
     &                          +CVVOLM(ICVL)*diag(ICVL)
          dvdt(ICVL,1)=CVVOLM(ICVL)*(dvdt(ICVL,1))
          dvdt(ICVL,2)=CVVOLM(ICVL)*(dvdt(ICVL,2))
          dvdt(ICVL,3)=CVVOLM(ICVL)*(dvdt(ICVL,3))
        enddo

      else
        if(ieul2ph>0) then 

        elseif(ical_vof==1.or.ical_cavi==1.or.ical_dens==4) then   ! zhang-vof
        else
          do 236 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          DO 237 ICVL=ICVS,ICVE
          timtrm=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          diag(ICVL)=timtrm+CVVOLM(ICVL)*diag(ICVL)
 237      enddo
          do 238 IV=1,3
          DO ICVL=ICVS,ICVE
          dvdt(ICVL,IV)=CVVOLM(ICVL)*(dvdt(ICVL,IV))
          enddo
 238      enddo
 236      enddo
        endif
      endif
!
! --- 
!
!-------------------------------------------
! --- < 2.2 viscous term >-- diffusion term
!-------------------------------------------
!
      if(intgvv.eq.Runge_Kutta) then
        if(ical_vect)  call FFRABORT
     &  (1,'MSG: Vector-Ver. NOT support Runge_Kutta')
        do 160 IIMAT=1,NMAT     !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do IV=1,3
        dvdd(ICVS:ICVE,IV,1)=dvdt(ICVS:ICVE,IV)
        enddo
 160    enddo
        dvdt(:,1:3)=0.d0
        dvdd(:,:,2)=0.d0
      endif
!
!      if(viscous) then
        if(ieul2ph>0) then
        else
          ICODE=1
          call vel_diff(iphs,ismpl,iter,rdeltt,
     &    LVEDGE,LBC_SSF,LCYCSF,kdbv,LFUTAU,
     &    MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,FRSTCV,
     &    grdc,rvd,aks,
     &    LCYCOLD,wifsld,OPPANG,vctr,
     &    rmu,rmut,rho(:,1),vel,utau,dsclv,dvdt(:,1:3),rmue,ICODE)
        endif
!      endif
!
!
       do 163 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).OR.IMAT.lt.0) cycle
        if(ivscvv(IMAT).ne.cal) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          rmue(ICVS:ICVE)=0.d0
          dsclv(ICVS:ICVE)=0.d0
          rmue(IDCS:IDCE)=0.d0
          dsclv(IDCS:IDCE)=0.d0
        endif
 163    continue
!------------------------------
! --- < 2.3 convection term >--
!------------------------------
      if(ieul2ph>0) then
      else
         call conv_term_vel
     &   (iphs,icnvvv,lmtrvv,deltt,
     &    LVEDGE,LBC_SSF,LCYCSF,vctr,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &    rvd(:,1:2),dvdd(:,:,1),diag,
     &    LCYCOLD,wifsld,
     &    grdc,rva(:,1),vel,dvdt(:,1:3))
      endif
!--------------------------------------------------
! --- 4th Order Runge_Kutta explicit time integral
!--------------------------------------------------
!
      if(intgvv.eq.Runge_Kutta) then
        ICODE=0
        call grad_cell(1,15,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,prs(:,1),grdc(:,:,1))
!---------------------------------------------------------------
! NOTE:
! --- dvdt(:,:)  : convection+diffusion term
! --- dvdd(:,:,1): source term(=dp/dx+F)
! --- dvdd(:,:,2): sum of 4-step's convection+diffusion term
!---------------------------------------------------------------
        viscous=.true.
!----------------------------
! --- 1st Step: FAI(*,n+1/2) 
!----------------------------
!
        if(ieul2ph>0) then
        else
          do 811 IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdd(ICVL,IV,1)=dvdd(ICVL,IV,1)-CVVOLM(ICVL)*grdc(ICVL,IV,1)
          dvdd(ICVL,IV,2)=dvdd(ICVL,IV,2)+r6*dvdt(ICVL,IV)
          dvdt(ICVL,IV)=r2*(dvdt(ICVL,IV)+dvdd(ICVL,IV,1))
     &           +vel(ICVL,IV,2)*rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          vel(ICVL,IV,1)=dvdt(ICVL,IV)/diag(ICVL)
          enddo
          enddo
 811      enddo
        endif
!----------------------------------
! --- 2nd Step: FAI(**,n+1/2)
!----------------------------------
        dvdt(:,1:3)=0.d0
        if(viscous) then
          call vel_diff(iphs,ismpl,iter,rdeltt,
     &    LVEDGE,LBC_SSF,LCYCSF,kdbv,LFUTAU,
     &    MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,FRSTCV,
     &    grdc,rvd,aks,
     &    LCYCOLD,wifsld,OPPANG,vctr,
     &    rmu,rmut,rho(:,1),vel,utau,dsclv,dvdt(:,1:3),rmue,ICODE)
        endif
        call conv_term_vel
     &   (iphs,icnvvv,lmtrvv,deltt,
     &    LVEDGE,LBC_SSF,LCYCSF,vctr,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &    rvd(:,1:2),dvdd(:,:,1),diag,
     &    LCYCOLD,wifsld,
     &    grdc,rva(:,1),vel,dvdt(:,1:3)
     &    )
        if(ieul2ph>0) then
        else
          do 821 IIMAT=1,NMAT !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdd(ICVL,IV,2)=dvdd(ICVL,IV,2)+2.d0*r6*dvdt(ICVL,IV)
          dvdt(ICVL,IV)=r2*(dvdt(ICVL,IV)+dvdd(ICVL,IV,1))
     &               +vel(ICVL,IV,2)*rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          vel(ICVL,IV,1)=dvdt(ICVL,IV)/diag(ICVL)
          enddo
          enddo
 821      continue
        endif
!----------------------------------
! --- 3rd Step: FAI(*,n+1)
!----------------------------------
        dvdt(:,1:3)=0.d0
        if(viscous) then
          call vel_diff(iphs,ismpl,iter,rdeltt,
     &    LVEDGE,LBC_SSF,LCYCSF,kdbv,LFUTAU,
     &    MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,FRSTCV,
     &    grdc,rvd,aks,
     &    LCYCOLD,wifsld,OPPANG,vctr,
     &    rmu,rmut,rho(:,1),vel,utau,dsclv,dvdt(:,1:3),rmue,ICODE)
        endif
          call conv_term_vel
     &   (iphs,icnvvv,lmtrvv,deltt,
     &    LVEDGE,LBC_SSF,LCYCSF,vctr,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &    rvd(:,1:2),dvdd(:,:,1),diag,
     &    LCYCOLD,wifsld,
     &    grdc,rva(:,1),vel,dvdt(:,1:3)
     &    )
        if(ieul2ph>0) then
        else
          do 831 IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdd(ICVL,IV,2)=dvdd(ICVL,IV,2)+2.d0*r6*dvdt(ICVL,IV)
          dvdt(ICVL,IV)=dvdt(ICVL,IV)+dvdd(ICVL,IV,1)
     &              +vel(ICVL,IV,2)*rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          vel(ICVL,IV,1)=dvdt(ICVL,IV)/diag(ICVL)
          enddo
          enddo
 831      continue
        endif
!----------------------------------
! --- 4th Step:
!----------------------------------
        dvdt(:,1:3)=0.d0
        if(viscous) then
          call vel_diff(iphs,ismpl,iter,rdeltt,
     &    LVEDGE,LBC_SSF,LCYCSF,kdbv,LFUTAU,
     &    MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,FRSTCV,
     &    grdc,rvd,aks,
     &    LCYCOLD,wifsld,OPPANG,vctr,
     &    rmu,rmut,rho(:,1),vel,utau,dsclv,dvdt(:,1:3),rmue,ICODE)
        endif
        call conv_term_vel
     &   (iphs,icnvvv,lmtrvv,deltt,
     &    LVEDGE,LBC_SSF,LCYCSF,vctr,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &    rvd(:,1:2),dvdd(:,:,1),diag,
     &    LCYCOLD,wifsld,
     &    grdc,rva(:,1),vel,dvdt(:,1:3)
     &    )
        if(ieul2ph>0) then
        else
          do 841 IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) goto 841
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdd(ICVL,IV,2)=dvdd(ICVL,IV,2)+r6*dvdt(ICVL,IV)
          dvdt(ICVL,IV)=dvdd(ICVL,IV,2)+dvdd(ICVL,IV,1)
! --- vel(:,:,1): old time-term
          vel(ICVL,IV,1)=vel(ICVL,IV,2)
     &                  *rho(ICVL,2)*CVVOL0(ICVL)*rdeltt/diag(ICVL)
          enddo
          enddo
 841      enddo
        endif
      endif
!----------------------------------
! --- < 2.4 pressure gradient >-- 
!---------------------------------- 3333
!----------------------3333
! --- Body force 
!----------------------
      if(.false.) then
        do 324 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        kd=kdbcnd(0,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd==kxnone)then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dum1=bdyfrc(ICVA,1,iphs)*(SFCENT(1,ICFL)-CVCENT(1,ICVA))
     &        +bdyfrc(ICVA,2,iphs)*(SFCENT(2,ICFL)-CVCENT(2,ICVA))
     &        +bdyfrc(ICVA,3,iphs)*(SFCENT(3,ICFL)-CVCENT(3,ICVA))
          prs(ICVB,1)=prs(ICVA,1)+dum1
          enddo
        endif
 324    enddo
      endif
!     
      if(iLBF_P==1) then     !3333
        call grad_cell_body(1,16,
     &   MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &   LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,
     &   prs(:,1),bdyfrc(:,:,iphs),grdc(:,:,1))
      elseif(iLBF_P==3) then   !6666
        call grad_cell_least(1,18,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,rva(:,1),
     &  prs(:,1),grdc(:,:,1))
      elseif(iLBF_P==0.or.iLBF_P==2.or.iLBF_P==4.or.iLBF_P==5) then      !6666
        ICVS=115
        if(iter==12181) ICVS=116
        call grad_cell(1,ICVS,
     &   MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &   LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,prs(:,1),grdc(:,:,1))
      endif
!
      call corre_grdc(LBC_SSF,LVEDGE,SFAREA,mat_cal,LCYCSF,grdc)
!--------------------------------------------
! --- Adams-Bashforth 
!--------------------------------------------
!
      IF(ical_vect) THEN  !chang2
        if(ismpl.eq.1) then
!CIDR NODEP
          do ICVL=ICVS_V,ICVE_V    !index_c(myid)+1,index_c(myid+1)
          dvdd(ICVL,1,1)=dvdt(ICVL,1)
          dvdd(ICVL,2,1)=dvdt(ICVL,2)
          dvdd(ICVL,3,1)=dvdt(ICVL,3)
          enddo
        endif
        if(intgvv.eq.Adams_Bashforth) then
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V!index_c(myid)+1,index_c(myid+1)
          dvdt(ICVL,1)=(1.d0+calph)*dvdt(ICVL,1)
     &                     -calph*dvdd(ICVL,1,2)
          dvdt(ICVL,2)=(1.d0+calph)*dvdt(ICVL,2)
     &                     -calph*dvdd(ICVL,2,2)
          dvdt(ICVL,3)=(1.d0+calph)*dvdt(ICVL,3)
     &                     -calph*dvdd(ICVL,3,2)
          enddo
        elseif(intgvv.eq.Adams_Moulton) then
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V!index_c(myid)+1,index_c(myid+1)
          dvdt(ICVL,1)=calph*(5.d0*dvdt(ICVL,1)
     &                  +8.d0*dvdd(ICVL,1,1)
     &              -dvdd(ICVL,1,2))
          dvdt(ICVL,2)=calph*(5.d0*dvdt(ICVL,2)
     &                  +8.d0*dvdd(ICVL,2,1)
     &              -dvdd(ICVL,2,2))
          dvdt(ICVL,3)=calph*(5.d0*dvdt(ICVL,3)
     &                  +8.d0*dvdd(ICVL,3,1)
     &              -dvdd(ICVL,3,2))
          enddo
!          enddo
        elseif(intgvv.eq.Crank_Nicolson) then
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V !index_c(myid)+1,index_c(myid+1)
          dvdt(ICVL,1)=calph*dvdt(ICVL,1)+calph*dvdd(ICVL,1,1)
          dvdt(ICVL,2)=calph*dvdt(ICVL,2)+calph*dvdd(ICVL,2,1)
          dvdt(ICVL,3)=calph*dvdt(ICVL,3)+calph*dvdd(ICVL,3,1)
          enddo
        endif
        if(intgvv.ne.Runge_Kutta) then
!CIDR NODEP
          do ICVL=ICVS_V,ICVE_V   !index_c(myid)+1,index_c(myid+1)
          dvdt(ICVL,1)=dvdt(ICVL,1)
     &                +rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
     &                *(vel(ICVL,1,2)-vel(ICVL,1,1))
     &            -CVVOLM(ICVL)*grdc(ICVL,1,1)
          dvdt(ICVL,2)=dvdt(ICVL,2)
     &                +rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
     &                *(vel(ICVL,2,2)-vel(ICVL,2,1))
     &            -CVVOLM(ICVL)*grdc(ICVL,2,1)
          dvdt(ICVL,3)=dvdt(ICVL,3)
     &                +rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
     &                *(vel(ICVL,3,2)-vel(ICVL,3,1))
     &            -CVVOLM(ICVL)*grdc(ICVL,3,1)
          enddo
        endif
      else
        if(intgvv.eq.Adams_Bashforth) then
          if(ismpl.eq.1) then
            DO 124 IIMAT=1,NMAT   !ICV=1,NCV
            if(.not.mat_cal(IIMAT)) cycle
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do IV=1,3
            do ICVL=ICVS,ICVE
            dvdd(ICVL,IV,1)=dvdt(ICVL,IV)
            enddo
            enddo
 124        enddo
          endif
          DO 134 IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do 135 IV=1,3
          do ICVL=ICVS,ICVE
          dvdt(ICVL,IV)=(1.d0+calph)*dvdt(ICVL,IV)
     &                 -calph*dvdd(ICVL,IV,2)
          enddo
 135      continue
 134      continue
!--------------------------------------------
! --- Adams_Moulton
!--------------------------------------------
        elseif(intgvv.eq.Adams_Moulton) then
          if(ismpl.eq.1) then
          DO 136 IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdd(ICVL,IV,1)=dvdt(ICVL,IV)
          enddo
          enddo
 136      continue
          endif
          DO 139 IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdt(ICVL,IV)=calph*(5.d0*dvdt(ICVL,IV)+8.d0*dvdd(ICVL,IV,1)
     &              -dvdd(ICVL,IV,2))
          enddo
          enddo
 139      continue
!--------------------------------------------
! --- Crank_Nicolson
!--------------------------------------------
        elseif(intgvv.eq.Crank_Nicolson) then
          if(ismpl.eq.1) then
          DO 137 IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdd(ICVL,IV,1)=dvdt(ICVL,IV)
          enddo
          enddo
 137      continue
          endif
          DO 138 IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do IV=1,3
          do ICVL=ICVS,ICVE
          dvdt(ICVL,IV)=calph*dvdt(ICVL,IV)+calph*dvdd(ICVL,IV,2) 
          enddo
          enddo
 138      enddo
        endif
!--------------------------------------------
! --- time difference term and pressure term:
!--------------------------------------------

        if(intgvv.ne.Runge_Kutta) then
          if(ieul2ph>0) then
          elseif(ical_vof==1) then
          elseif(ical_dens==4) then
          else
            dum1=1.d0
!            if(iLBF_P==2.or.iLBF_P==3) dum1=0.d0
            if(iLBF_P==2) dum1=0.d0
            do 140 IIMAT=1,NMAT
              if(.not.mat_cal(IIMAT)) cycle
              IMAT=MAT_NO(IIMAT)
              if(IMAT<0) cycle
              porodum=1.d0
              if(porous(IMAT)==idarcy2) porodum=porosty(IMAT)
              ICVS=MAT_CVEXT(IIMAT-1)+1
              ICVE=MAT_CVEXT(IIMAT)
              ICVE=MAT_INDEX(IIMAT)
              do 141 IV=1,3
              do ICVL=ICVS,ICVE
              timtrm=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
              dvdt(ICVL,IV)=dvdt(ICVL,IV)
     &           +timtrm*(vel(ICVL,IV,2)-vel(ICVL,IV,1))
     &           -porodum*dum1*CVVOLM(ICVL)*grdc(ICVL,IV,1)
              enddo
 141          enddo
 140        enddo
          endif
        endif
      endif

!-------------------------------
! --- < 2.2 clear solid part >--
!-------------------------------
!      if(.not.ical_vect) then
        do 120 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT.lt.0) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do IV=1,3
            DO 125 ICVL=ICVS,ICVE
            dvdt(ICVL,IV)=0.d0
 125        enddo
            enddo
          endif
 120    enddo
!      endif
!5555-1
      if(ieul2ph>0) then
      endif
!------------------------
! --- porous madia
!------------------------
      if(ical_porous==1) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.le.0) cycle
        IMAT_U=nofld(IMAT)
        if(IMAT_U>N_pors) then
          if(porous(IMAT)==idarcy2) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            DO ICVL=ICVS,ICVE
            dvdt(ICVL,1)=dvdt(ICVL,1)*porosty(IMAT)
            dvdt(ICVL,2)=dvdt(ICVL,2)*porosty(IMAT)
            dvdt(ICVL,3)=dvdt(ICVL,3)*porosty(IMAT)
            enddo
          endif
        endif
        enddo
      endif
      if(ical_poroBC==1) then 
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(poro(nb)==1) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL) 
          ICFP=LCYCSF(IBFL)
          ICVL=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)
          dvdt(ICVL,1)=dvdt(ICVL,1)*prosty(nb)
          dvdt(ICVL,2)=dvdt(ICVL,2)*prosty(nb)
          dvdt(ICVL,3)=dvdt(ICVL,3)*prosty(nb)
!          dvdt(ICVP,1)=dvdt(ICVP,1)*prosty(nb)
!          dvdt(ICVP,2)=dvdt(ICVP,2)*prosty(nb)
!          dvdt(ICVP,3)=dvdt(ICVP,3)*prosty(nb)
          enddo
        endif
        enddo
      endif
!-------------------------
! --- 
!-------------------------
!      if(ical_FC==PEFC) then
!        do IIMAT=1,NMAT
!        IMAT=MAT_NO(IIMAT)
!        if(IMAT.le.0) cycle
!        IMAT_U=nofld(IMAT)
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_CVEXT(IIMAT)
!        DO ICVL=ICVS,ICVE
!          dum1=1.d0-aks(ICVL,ical_s,1)
!          dvdt(ICVL,1)=dvdt(ICVL,1)*dum1
!          dvdt(ICVL,2)=dvdt(ICVL,2)*dum1
!          dvdt(ICVL,3)=dvdt(ICVL,3)*dum1
!        enddo
!        enddo
!      endif
!-----------------------------
! --- < 3. Time integration >-
!--------------------------------------------
! --- Euler implicit method in time-integral 
!--------------------------------------------
      if(intgvv.eq.euleri.or.
     &   intgvv.eq.Adams_Moulton.or.
     &   intgvv.eq.Crank_Nicolson) then
!
        call bc_velcofd(LVEDGE,LBC_SSF,mat_cal,cofd)
!
        duma=0.d0
	ndf=1
	nq=3
        IMODE=5
        if(ieul2ph>0) then
        else
          rvd=0.d0
          if(ical_vect) then  !NOVECT
            ndiag=1
            IMODE=5
            call solve_cnvdif_vect
     & (.true.,1,ndf,nq,'v',time,ndiag,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdbv,rva,rvd(:,1:2),rmue,diag,
     &  LCYCOLD,wifsld,OPPANG,
     &  cofd,dsclv,duma,dvdt(:,1:3),deltt,vctr,
     &  iterv,repsv,aepsv,iter,ierr1,IMODE)
          else
            call solve_cnvdif
     &  (.true.,1,ndf,nq,'v',time,relaxv,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,SFCENT,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdbv,rva,rvd(:,1:2),rmue,diag,
     &  LCYCOLD,wifsld,OPPANG,CVVOLM,
     &  cofd,dsclv,duma,dvdt(:,1:3),deltt,
     &  iterv,repsv,aepsv,iter,ierr1,IMODE)
          endif
        endif
        if( ierr1.ne.0 ) goto 9999
      elseif(intgvv.eq.eulere
     &       .or.intgvv.eq.Adams_Bashforth
     &       .or.intgvv.eq.Runge_Kutta 
     &       )then
!----------------------------------------------------------------
! --- Euler and Adams_Bashforth explicit method in time-integral
!----------------------------------------------------------------
        if(ical_vect) then
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V!index_c(myid)+1,index_c(myid+1)
          dvdt(ICVL,1)=dvdt(ICVL,1)/diag(ICVL)
          dvdt(ICVL,2)=dvdt(ICVL,2)/diag(ICVL)
          dvdt(ICVL,3)=dvdt(ICVL,3)/diag(ICVL)
          enddo
        else
          do IIMAT=1,NMAT   !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do 206 IV=1,3
          do ICVL=ICVS,ICVE
          dum1=1.d0/diag(ICVL)
          dvdt(ICVL,IV)=dum1*dvdt(ICVL,IV)
          enddo
 206      enddo
 205      enddo
        endif
      endif
!-------------
! --- E2P
!-------------5555-????
      if(ieul2ph>0.and..false.) then
      endif
!
!----------------------------------
! --- Convergence Residual error:
!----------------------------------
      if(ical_vect) then
        errv=ZERO
        errabs=ZERO
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V!index_c(myid)+1,index_c(myid+1)
        errv  (1)=errv  (1)+(dvdt(ICVL,1)* dvdt(ICVL,1))
        errabs(1)=errabs(1)+( vel(ICVL,1,1)*vel(ICVL,1,1))
        errv  (2)=errv  (2)+(dvdt(ICVL,2)* dvdt(ICVL,2))
        errabs(2)=errabs(2)+( vel(ICVL,2,1)*vel(ICVL,2,1))
        errv  (3)=errv  (3)+(dvdt(ICVL,3)* dvdt(ICVL,3))
        errabs(3)=errabs(3)+( vel(ICVL,3,1)*vel(ICVL,3,1))
        enddo
      else
        errv=ZERO
        errabs=ZERO
        Do 220 IIMAT=1,NMAT   !ICV=1,NCVIN
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        do 222 iv=1,3
        do ICVL=ICVS,ICVE
        errv(iv)=errv(iv)+(dvdt(ICVL,iv)*dvdt(ICVL,iv))
        errabs(iv)=errabs(iv)+(vel(ICVL,iv,1)*vel(ICVL,iv,1))
        enddo
 222    continue
 220    continue
      endif
!
      IF(NPE.GT.1) THEN
        do i=1,3
        CALL hpcrsum(errv(i))
        CALL hpcrsum(errabs(i))
        enddo
      ENDIF
!
      do 225 i=1,3
      errv(i)=dsqrt(errv(i))/(dsqrt(errabs(i))+SML1)
  225 continue
!
      if( xredc ) errv(1)=0.d0
      if( yredc ) errv(2)=0.d0
      if( zredc ) errv(3)=0.d0
!----------------------------------------------------
! --- correcting velocity: (not have pressure term)  
!----------------------------------------------------
!
      if(ical_vect) then  !chang3
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V    !index_c(myid)+1,index_c(myid+1)
        vel(ICVL,1,1)=vel(ICVL,1,1)+dvdt(ICVL,1)*relaxv(1)
        vel(ICVL,2,1)=vel(ICVL,2,1)+dvdt(ICVL,2)*relaxv(1)
        vel(ICVL,3,1)=vel(ICVL,3,1)+dvdt(ICVL,3)*relaxv(1)
        enddo
      else
        do 212 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) goto 212
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 210 iv=1,3
        do ICVL=ICVS,ICVE
        vel(ICVL,iv,1)=vel(ICVL,iv,1)+dvdt(ICVL,iv)*relaxv(IIMAT)
        enddo
 210    continue
 212    continue
      endif
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (3,MXALLCV,NCV,vel(:,:,1))
      ENDIF
!
!--------------------
! --- Solid material
!--------------------
!
      if(ical_vect) then
      else
        do 126 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          DO IV=1,3
          vel(ICVS:ICVE,IV,1)=0.d0
          enddo         
        endif
 126    continue
      endif
!
      return
!
 9999 continue
!
      if(my_rank.eq.ROOT) write(ifle,*) '(vel_admin)'
      ierror=1
!
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!
!< #1. calculate gravity term >
!
!=================================================
      subroutine gravity(gtrm)
!=================================================
!
      use module_gravity ,only : ggg,ramb,tamb,beta,iave_rho_T,
     &                           ramb2,tamb2,beta2
      use module_material,only : lclsd,iclosd!,KMAT_S,MAT_S,iclosd
      use module_metrix  ,only : vsm,rsm,tsm
!
      implicit none
!
      real*8,intent(inout) :: gtrm(MXALLCV,1:3)
      real*8  :: gg1(3),gg2(3),dr1,dr2,betax,rambx,tambx
      integer :: i,l,n,KMAT,ierr1=0
!
! --- 
!
      vsm=0.d0
      rsm=0.d0
      tsm=0.d0
      if(iave_rho_T==1.or.iave_rho_T==3) then 
        do 100 IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle 
          IMAT=MAT_NO(IIMAT) 
          if(IMAT.gt.0) then 
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do 110 ICVL=ICVS,ICVE
            vsm(IIMAT)=vsm(IIMAT)+CVVOLM(ICVL)
            rsm(IIMAT)=rsm(IIMAT)+CVVOLM(ICVL)*rho(ICVL,1)
            tsm(IIMAT)=tsm(IIMAT)+CVVOLM(ICVL)*tmp(ICVL)
 110        enddo
          endif
 100    enddo
!
        if(NPE.gt.1) then
!          DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          call hpcrsum(vsm(IIMAT))
!          call hpcrsum(rsm(IIMAT))
!          call hpcrsum(tsm(IIMAT))
!          vsm(0)=0.d0
!          rsm(0)=0.d0
!          tsm(0)=0.d0
!          enddo
          DO IIMAT=1,NMAT
          call hpcrsum(vsm(IIMAT))
          call hpcrsum(rsm(IIMAT))
          call hpcrsum(tsm(IIMAT))
          enddo
        endif
!
        do 200 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.le.0) cycle
        if(vsm(IIMAT).gt.SML) then
          rsm(IIMAT)=rsm(IIMAT)/vsm(IIMAT)
          tsm(IIMAT)=tsm(IIMAT)/vsm(IIMAT)
        endif
 200    enddo
!
      endif
!
                          gg1=ggg
      if(idrdp.eq.incomp) gg1=0.d0
!
      gg2=ggg-gg1
!
      do 120 IIMAT=1,NMAT !ICV=1,NCV
      if(.not.mat_cal(IIMAT))  cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.gt.0) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(iave_rho_T==1) then
          if(ICAL_VECT) then
            do i=1,3
            do ICVL=ICVS,ICVE
            dr1=rho(ICVL,1)-rsm(IIMAT)
            dr2=rho(ICVL,1)*beta*(tmp(ICVL)-tsm(IIMAT))
            gtrm(ICVL,i)=gtrm(ICVL,i)+gg1(i)*dr1-gg2(i)*dr2
            enddo
            enddo
          else
            do 125 ICVL=ICVS,ICVE
            dr1=rho(ICVL,1)-rsm(IIMAT)
            dr2=rho(ICVL,1)*beta*(tmp(ICVL)-tsm(IIMAT))
            do 121 i=1,3 
            gtrm(ICVL,i)=gtrm(ICVL,i)+gg1(i)*dr1-gg2(i)*dr2 
 121        enddo
 125        enddo
          endif
        elseif(iave_rho_T==2) then
          do i=1,3
          do ICVL=ICVS,ICVE
          gtrm(ICVL,i)=gtrm(ICVL,i)+ggg(i)*rho(ICVL,1)
          enddo
          enddo
        elseif(iave_rho_T==3) then 
          do i=1,3 
          do ICVL=ICVS,ICVE
          gtrm(ICVL,i)=gtrm(ICVL,i)+ggg(i)*(rho(ICVL,1)-rsm(IIMAT))
          enddo
          enddo
        elseif(iave_rho_T==4) then 
          rambx=ramb
          do i=1,3
          do ICVL=ICVS,ICVE
          dr1=rho(ICVL,1)-rambx
          gtrm(ICVL,i)=gtrm(ICVL,i)+ggg(i)*dr1
          enddo
          enddo
        else
          betax=beta
          tambx=tamb
          rambx=ramb
          if(ical_vect) then
            do i=1,3
            do ICVL=ICVS,ICVE
            dr1=rho(ICVL,1)-rambx
            dr2=rho(ICVL,1)*betax*(tmp(ICVL)-tambx)
            gtrm(ICVL,i)=gtrm(ICVL,i)+gg1(i)*dr1-gg2(i)*dr2
            enddo
            enddo
          else
            do 127 ICVL=ICVS,ICVE
            dr1=rho(ICVL,1)-rambx
            dr2=rho(ICVL,1)*betax*(tmp(ICVL)-tambx)

            do 128 i=1,3
            gtrm(ICVL,i)=gtrm(ICVL,i)+gg1(i)*dr1-gg2(i)*dr2
 128        enddo
 127        enddo
          endif
!
!          if(iLBF_P==1.or.iLBF_P==2.or.iLBF_P==3) then !3333
          if(iLBF_P==1.or.iLBF_P==2) then !3333
            do ICVL=IDCS,IDCE
            dr1=rho(ICVL,1)-rambx
            dr2=rho(ICVL,1)*betax*(tmp(ICVL)-tambx)
            do i=1,3
            gtrm(ICVL,i)=gtrm(ICVL,i)+gg1(i)*dr1-gg2(i)*dr2
            enddo
            enddo
          endif
        endif
      endif
 120  continue
!
      end subroutine gravity
!
      end subroutine vel_admin
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine rva_admin(iphs,iter,ismpl,time,deltt,nsmpl,volctr,
     &  LVEDGE,LCYCSF,LFUTAU,LBC_SSF,
     &  SFAREA,SFCENT,CVCENT,CVVOLM,CVVOL0,wiface,FRSTCV,DISALL,
     &  rva,vel,rho,prs,xta,tmp,
     &  grdc,rvx,aks,bdyfrc,velf,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  prsbnd,
     &  LCYCOLD,wifsld,OPPANG,locmsh,FIELD_U,
     &  ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
!
      use module_io,only       : ifle,ifll
      use module_model,only    : idrdp,incomp,mach0,comp,ical_vect
     &                           ,nthrds,iLBF_P,PEFC,sles,SDES,
     &                            dles,icaltb,expfact,ical_dens
      use module_flags,only    : icalrv,face,intgvv,euleri,eulere,
     &                           Adams_Bashforth,
     &                           Adams_Moulton,Crank_Nicolson,
     &                           Runge_Kutta
      use module_material,only : nflud,nofld,lclsd,nsold,nosld,
     &                     ical_sld,rotati,ishaft,end,begin,rot_ang,
     &                     porosty,relaxv,relaxp,relaxrc
      use module_boundary,only : kdbcnd,kdilet,kvlglw,kvnslp,kxnone,
     &                           kdintr,kdfire,kdstag,kdpres,kdbuff,
     &                           kdolet,kdprdc,kdsymm,kdtchi,kdcvd,
     &                           distrb,rotsld,idis,LBC_pair,kdshutr,
     &                           LBC_INDEX,nbcnd,MAT_BCIDX,openout,
     &                           stagvel,kdsld,ical_buff,pkdwall,
     &                           masflg,sumflx,masflx,sumare,masbc,dvel
     &                           ,kdovst
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_vof     ,only : ical_vof
      use module_scalar,  only : ivof,iaph,ical_cavi,icavi,ivold,
     &                           ical_FC,ical_s
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_gravity ,only : ggg
      use module_initial ,only : rho0,rho02
      use module_FUEL    ,only : No_Mem,No_AGDL,
     &                           No_CGDL,vap_no
      use module_time,    only : i_steady
      use module_metrix,only   : SHUTFL,vctr
      use module_metrix,only   : msk
      use module_metrix,only  : bdyf      
!
! 1.  Update velocity field
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,ismpl,nsmpl,iphs
      real*8 ,intent(in)    :: deltt,time,volctr
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO(      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV(      MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)      
      integer,intent(in)    :: LFUTAU    (     MXCV)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (  MXALLCV)
      real*8 ,intent(in)    :: wiface    (  MXCVFAC)
      real*8 ,intent(in)    :: FRSTCV    (  MXSSFBC)
      real*8 ,intent(in)    :: DISALL    (  MXCV)
      real*8 ,intent(inout)    :: xta       (  MXCVFAC)
      real*8 ,intent(inout) :: rva       (  MXCVFAC,2)
      real*8 ,intent(inout) :: vel       (  MXALLCV,3,2)
      real*8 ,intent(inout) :: velf      (  MXCVFAC_B,3,2)
      real*8 ,intent(inout) :: rho       (  MXALLCV,2)
      real*8 ,intent(inout) :: prs       (  MXALLCV,2) 
      real*8 ,intent(inout) :: grdc      (  MXALLCV,3,3)
      REAL*8 ,INTENT(INOUT) :: rvx       (  MXALLCV,3)
      real*8 ,intent(inout) :: aks       (  MXALLCVR,mxrans)
      real*8 ,intent(in)    :: tmp       (  MXALLCV)
      real*8 ,intent(in)    :: PRSBND    (  MXSSFBC)
!
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
      integer,intent(out)   :: ierror
!      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      integer,intent(in)    :: locmsh (  MXSSFBC)
      real*8 ,intent(inout) :: bdyfrc(MXCV_B,3,2)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
!
! --- [local entities]
!
      real*8  :: dx2,dy2,dz2,dll,vel_av(3),vell(3,3),usum,rk,
     &           dxx,dyx,dzx
      real*8  :: dxo,dyo,dzo,avsfn,avsfo,dln,dlo !,expfact=1.d0
      real*8  :: wi1,wi2,grx,gry,grz,gf0,gf1,gf1o
      real*8  :: dx,dy,dz,dpAB,dl,dlvect!,dlvecto
      real*8  :: ru,rv,rw,ru1,rv1,rw1,ru2,rv2,rw2
      real*8  :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,
     &           costh1,sinth1,costh2,sinth2,th(2),avsf
      real*8  :: ru3,rv3,rw3,pors
      real*8  :: sf11,sf12,sf13,dl1
      real*8  :: ds1,ds2,ds3
      real*8  :: rvaicf,rvaicfp,calph,dump,rdeltt
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,nx,ierr1,ndf,nq
      integer :: ICOM,IMD,ICH,IFLD,ICTP,icfo,icvb1,icvbo
      integer :: IMAT,IIMAT,IMAT_U,ICVS,ICVE,ICVL,IWL,IMAT2,IIMAT2
      integer :: ICFS,ICFE,ICFL
!      integer :: IO,IN,IP,IS,IW
      integer :: IC,IV,IE,ICF,ICV,IBF,NB,KD,IDC,ICVFIX,IDCA,IDCB,IDCBO
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ICVLA,ICVLB,IBFS,IBFE,IMODE,IDCS,IDCE,IBFL
!      integer :: ICVBO,ICFO,ICVB1
      real*8  :: wii1,wii2,wii1O,wii2O
      real*8  :: grdf(3),gf2,gf3
      real*8  :: radius(3),velref(3),unit(3,2),
     &           rbb(3,3,2),fbb(3,3,2),grdff(3,3),grdfA(3),grdfB(3),
     &           dxA,dyA,dzA,dxB,dyB,dzB
      integer :: IIMATS(2),IMATS(2),ISLD,ISLD2,myid,imasflg=1
      real*8  :: org_x1,org_y1,org_x2,org_y2,duf,dum_mas,dum_are,alpha,
     &           dum_2,dum_1
      integer :: IBFS1,IBFE1,nb2
!     
!---------------------------------
! --- < 4. Mass Flux for CV face>-
!---------------------------------
!     
      rdeltt=deltt
      if(i_steady==3) then
        rdeltt=1.d0*relaxrc(1)        !1.d0          !rva_admin
      endif
!
      expfact=0.d0
      if(intgvv.eq.euleri.or.
     &   intgvv.eq.Adams_Moulton.or.
     &   intgvv.eq.Crank_Nicolson) then
        expfact=1.d0
      endif
!
      if(ical_sld>0) expfact=0.d0
      if(ical_buff==1) expfact=0.d0
      if(iLBF_P==2.and.ieul2ph>0) expfact=0.d0
      if(iLBF_P>0.and.ieul2ph>0) expfact=0.d0
      if(ical_cavi==1.or.ical_dens==4) expfact=0.d0
      if(ical_vof==1) expfact=0.d0
      if(ical_FC==PEFC) expfact=0.d0
      if(intgvv.eq.Adams_Moulton.and.
     &  (icaltb==sles.or.icaltb==SDES.OR.icaltb==dles)) then
        expfact=0.d0
      endif
      if(ical_vect) expfact=0.d0  !chang4
!
!
      if(iLBF_P==1) then  !   3333
        call grad_cell_body(1,16,
     &   MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &   LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,
     &   prs(:,1),bdyfrc(:,:,iphs),grdc(:,:,1))
      elseif(iLBF_P==3) then  !6666
        call grad_cell_least(1,18,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,rva(:,1),
     &  prs(:,1),grdc(:,:,1))
      elseif(iLBF_P==0.or.iLBF_P==2.or.iLBF_P==4.or.iLBF_P==5) then      !6666
        call grad_cell(1,17,
     &   MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &   LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,prs(:,1),grdc(:,:,1))
      endif
!
!      if(ieul2ph==1) then
!        call grad_cell(1,17,
!     &   MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
!     &   LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,aks(:,iaph(iphs)),grdc(:,:,2))
!      endif
!
!----------------------------------------------
! --- Boundary condition setting [Dummy cell] -
!----------------------------------------------!
      call dc_symprv
     & (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  LCYCOLD,wifsld,OPPANG,
     &  SFAREA,SFCENT,grdc(:,:,1),1,0)  !zh000
!
      call dc_symvel
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    LCYCOLD,wifsld,OPPANG,
     &    SFAREA,SFCENT,vel(:,:,1),0,1)
!
      call corre_grdc(LBC_SSF,LVEDGE,SFAREA,mat_cal,LCYCSF,grdc)
!
!
!
      if(.false.) then
        do nb=1,nbcnd
        if(pkdwall(nb)) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVL=LVEDGE(1,ICFL)
          dum1=SFAREA(1,ICFL)*bdyf(ICVL,1)
     &        +SFAREA(2,ICFL)*bdyf(ICVL,2)
     &        +SFAREA(3,ICFL)*bdyf(ICVL,3)
          do IV=1,3
          grdc(ICVL,IV,1)=grdc(ICVL,IV,1)+SFAREA(IV,ICFL)*dum1
          enddo
          enddo
        endif
        enddo
      endif
!
!----------------------------------------------------------
! --- < 4.1 temporary velocity at cell center > ---
! --- Remove pressure term from moment equation.(Rhie-Chow)
!----------------------------------------------------------
      if(ical_vect) then  !chang4
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V
          rvx(ICVL,1)=rho(ICVL,1)*vel(ICVL,1,1)
     &               -rho(ICVL,2)*vel(ICVL,1,2)*volctr
     &                   +rdeltt*grdc(ICVL,1,1)

          rvx(ICVL,2)=rho(ICVL,1)*vel(ICVL,2,1)
     &               -rho(ICVL,2)*vel(ICVL,2,2)*volctr
     &                   +rdeltt*grdc(ICVL,2,1)

          rvx(ICVL,3)=rho(ICVL,1)*vel(ICVL,3,1)
     &               -rho(ICVL,2)*vel(ICVL,3,2)*volctr
     &                   +rdeltt*grdc(ICVL,3,1)
        enddo
        call dc_symprv
     &    (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &     LCYCOLD,wifsld,OPPANG,
     &     SFAREA,SFCENT,rvx,2,2)
      else
        if(ieul2ph>0) then 
        elseif(ical_vof==1) then
        elseif(ical_cavi==1.or.ical_FC==PEFC) then
          do IIMAT=1,NMAT
            IMAT=MAT_NO(IIMAT)
            if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            pors=porosty(IMAT) 
            do IV=1,3
            do ICVL=ICVS,ICVE
            rvx(ICVL,IV)=
     &                   rho(ICVL,1)*vel(ICVL,IV,1)
     &                  +rdeltt*grdc(ICVL,IV,1)*pors
            enddo
            enddo
          enddo
          call dc_symprv
     &    (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &     LCYCOLD,wifsld,OPPANG,
     &     SFAREA,SFCENT,rvx,2,2)
        elseif(ical_dens==4) then
        elseif(idrdp==comp) then
          dum_2=1.d0
          if(iLBF_P==2) dum_2=0.d0 
          do IIMAT=1,NMAT
            if(.not.mat_cal(IIMAT)) cycle
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do IV=1,3
            do ICVL=ICVS,ICVE
            rvx(ICVL,IV)=rho(ICVL,1)*vel(ICVL,IV,1)
     &                  -rho(ICVL,2)*vel(ICVL,IV,2)*volctr
     &          +dum_2*rdeltt*grdc(ICVL,IV,1)*relaxv(IIMAT)
            enddo
            enddo
          enddo
          call dc_symprv
     &    (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &     LCYCOLD,wifsld,OPPANG,
     &     SFAREA,SFCENT,rvx,2,2)
          
        else
          dum_2=1.d0
!          if(iLBF_P==2.or.iLBF_P==3) dum_2=0.d0   !7777
          if(iLBF_P==2) dum_2=0.d0   !7777
          do 300 IIMAT=1,NMAT  !ICV=1,NCV
            if(.not.mat_cal(IIMAT)) cycle
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do 302 IV=1,3
            do ICVL=ICVS,ICVE
            rvx(ICVL,IV)=rho(ICVL,1)*vel(ICVL,IV,1)
     &                  -rho(ICVL,2)*vel(ICVL,IV,2)*volctr
     &          +dum_2*rdeltt*grdc(ICVL,IV,1)*relaxv(IIMAT)
            enddo
 302        enddo
 300      enddo
          call dc_symprv
     &    (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &     LCYCOLD,wifsld,OPPANG,
     &     SFAREA,SFCENT,rvx,2,2)
        endif
      endif
!
!-----------------------------------------------------------
! --- < 4.2 temporary velocity at cell face >--(Rhie-Chow)  
!-----------------------------------------------------------
!
      if(ical_vect) then   !NOVECT
!CIDR NODEP
        IIMAT=1
        DO ICFL=ICFS_V,ICFE_V
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
          dx=dx/dl
          dy=dy/dl
          dz=dz/dl
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          gf1=(prs(ICVB,1)-prs(ICVA,1))/(dl*dlvect)
     &       *(relaxv(IIMAT))
          ru=(wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1))
          rv=(wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2))
          rw=(wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3))
          grx=wi1*grdc(ICVA,1,1)+wi2*grdc(ICVB,1,1)
          gry=wi1*grdc(ICVA,2,1)+wi2*grdc(ICVB,2,1)
          grz=wi1*grdc(ICVA,3,1)+wi2*grdc(ICVB,3,1)
          dump=(SFAREA(1,ICFL)-dx)*grx
     &        +(SFAREA(2,ICFL)-dy)*gry
     &        +(SFAREA(3,ICFL)-dz)*grz
!
          gf2=SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw
          rva(ICFL,1)=rva(ICFL,2)*volctr+     !????zh
     &                rho(ICVA,1)*max(0.d0,xta(ICFL))
     &               +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &               +SFAREA(4,ICFL)*
     &                (gf2-rdeltt*(gf1+expfact*dump))
        enddo
!-------------------------------------------------------------------
      else
!-----------------------------
! --- Distance weigh factor
!-----------------------------
        if(ieul2ph>0) then       !  (1: 2P)
        elseif(iLBF_P==0.and.ieul2ph==0.and.ical_vof==1) then   !7777
        elseif(iLBF_P==0.and.ieul2ph==0.and.ical_dens==4) then   !7777
        elseif(iLBF_P==0.and.ieul2ph==0.and.idrdp==comp) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          pors=porosty(IMAT)
          IMAT_U=nofld(IMAT)
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
          dx=dx/dl
          dy=dy/dl
          dz=dz/dl
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          gf1=(prs(ICVB,1)-prs(ICVA,1))/(dl*dlvect)*(relaxv(IIMAT))
!zhang5  (7)
!          ru=(wi1*vel(ICVA,1,1)+wi2*vel(ICVB,1,1))
!          rv=(wi1*vel(ICVA,2,1)+wi2*vel(ICVB,2,1))
!          rw=(wi1*vel(ICVA,3,1)+wi2*vel(ICVB,3,1))
!          gf2=SFAREA(4,ICFL)*
!     &   (SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw)
!zhang5   (8)
          ru=(wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1))
          rv=(wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2))
          rw=(wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3))
          gf2=(SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw)
! --- ----------
          grx=wi1*grdc(ICVA,1,1)+wi2*grdc(ICVB,1,1)
          gry=wi1*grdc(ICVA,2,1)+wi2*grdc(ICVB,2,1)
          grz=wi1*grdc(ICVA,3,1)+wi2*grdc(ICVB,3,1)
!
          dump=(SFAREA(1,ICFL)-dx)*grx
     &        +(SFAREA(2,ICFL)-dy)*gry
     &        +(SFAREA(3,ICFL)-dz)*grz
!
! zhang5 (7) 
!          rva(ICFL,1)=
!     &                rho(ICVA,1)*max(0.d0,gf2+xta(ICFL))
!     &               +rho(ICVB,1)*min(0.d0,gf2+xta(ICFL))
!     &               +SFAREA(4,ICFL)*
!     &               (-rdeltt*pors*(gf1+expfact*dump)) 
! zhang5  (8)
          rva(ICFL,1)=
     &                rho(ICVA,1)*max(0.d0,xta(ICFL))
     &               +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &               +SFAREA(4,ICFL)*
     &               (gf2-rdeltt*pors*(gf1+expfact*dump)) 
!---------
          enddo
          enddo
         
        elseif((iLBF_P==0.or.iLBF_P==5.or.iLBF_P==3)
     &          .and.ieul2ph==0) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          pors=porosty(IMAT)
          IMAT_U=nofld(IMAT)
!          if(ical_FC==PEFC) then
!           expfact=0.d0
!           if(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem) then
!             expfact=1.d0
!           endif
!          endif
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!
!-------------------------------------------------------------------
!
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
          dx=dx/dl
          dy=dy/dl
          dz=dz/dl
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          gf1=(prs(ICVB,1)-prs(ICVA,1))/(dl*dlvect)*
     &           (relaxv(IIMAT))

          ru=(wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1))
          rv=(wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2))
          rw=(wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3))
          grx=wi1*grdc(ICVA,1,1)+wi2*grdc(ICVB,1,1)
          gry=wi1*grdc(ICVA,2,1)+wi2*grdc(ICVB,2,1)
          grz=wi1*grdc(ICVA,3,1)+wi2*grdc(ICVB,3,1)
!
          dump=(SFAREA(1,ICFL)-dx)*grx
     &        +(SFAREA(2,ICFL)-dy)*gry
     &        +(SFAREA(3,ICFL)-dz)*grz
!
          gf2=SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw
!
          rva(ICFL,1)=rva(ICFL,2)*volctr+     !????zh
     &                rho(ICVA,1)*max(0.d0,xta(ICFL))
     &               +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &               +SFAREA(4,ICFL)*
     &               (gf2-rdeltt*pors*(gf1+expfact*dump)) 

          enddo
          enddo
        elseif(iLBF_P==1.and.ieul2ph==0) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          pors=porosty(IMAT)
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          gf1=(prs(ICVB,1)-prs(ICVA,1))/dlvect
!
          dx=-(SFCENT(1,ICFL)-CVCENT(1,ICVA))     !fA
          dy=-(SFCENT(2,ICFL)-CVCENT(2,ICVA))
          dz=-(SFCENT(3,ICFL)-CVCENT(3,ICVA))
          dx2=-(SFCENT(1,ICFL)-CVCENT(1,ICVB))    !fB
          dy2=-(SFCENT(2,ICFL)-CVCENT(2,ICVB))
          dz2=-(SFCENT(3,ICFL)-CVCENT(3,ICVB))
          dum1=(wi1*dx+wi2*dx2)
          dum2=(wi1*dy+wi2*dy2)
          dum3=(wi1*dz+wi2*dz2)
!
          dum4=(bdyfrc(ICVB,1,iphs)-bdyfrc(ICVA,1,iphs))/dl!vect
          ds1=dum4*SFAREA(1,ICFL)
          ds2=dum4*SFAREA(2,ICFL)
          ds3=dum4*SFAREA(3,ICFL)
          dum5=dum1*ds1+dum2*ds2+dum3*ds3
          vel_av(1)=-dum5
!
          dum4=(bdyfrc(ICVB,2,iphs)-bdyfrc(ICVA,2,iphs))/dl!vect
          ds1=dum4*SFAREA(1,ICFL)
          ds2=dum4*SFAREA(2,ICFL)
          ds3=dum4*SFAREA(3,ICFL)
          dum5=dum1*ds1+dum2*ds2+dum3*ds3
          vel_av(2)=-dum5
!
          dum4=(bdyfrc(ICVB,3,iphs)-bdyfrc(ICVA,3,iphs))/dl!vect
          ds1=dum4*SFAREA(1,ICFL)
          ds2=dum4*SFAREA(2,ICFL)
          ds3=dum4*SFAREA(3,ICFL)
          dum5=dum1*ds1+dum2*ds2+dum3*ds3
          vel_av(3)=-dum5
!
          ru=(wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1))
          rv=(wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2))
          rw=(wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3))
          gf2=SFAREA(1,ICFL)*ru
     &       +SFAREA(2,ICFL)*rv
     &       +SFAREA(3,ICFL)*rw
          dump=SFAREA(1,ICFL)*vel_av(1)
     &        +SFAREA(2,ICFL)*vel_av(2)
     &        +SFAREA(3,ICFL)*vel_av(3)
!---------------------------------------------------------------
          rva(ICFL,1)=rva(ICFL,2)*volctr+     !????zh
     &                rho(ICVA,1)*max(0.d0,xta(ICFL))
     &               +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &               +SFAREA(4,ICFL)
     &               *(gf2-rdeltt*pors*(gf1-0.d0*expfact*dump)) !3333
!---------------------------------------------------------------
          enddo
          enddo
!        elseif((iLBF_P==2.or.iLBF_P==3).and.ieul2ph==0) then   !7777
        elseif((iLBF_P==2).and.ieul2ph==0) then   !7777
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          pors=porosty(IMAT) 
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          ru=(wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1))
          rv=(wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2))
          rw=(wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3))
!
          gf2=SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw
          rva(ICFL,1)=rva(ICFL,2)*volctr+ 
     &                rho(ICVA,1)*max(0.d0,xta(ICFL))
     &               +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &          +SFAREA(4,ICFL)*gf2 
          enddo
          enddo
        elseif(iLBF_P==4.and.ieul2ph==0) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dum1=wi1*rho(ICVA,1)+wi2*rho(ICVB,1)
          ru=velf(ICFL,1,1)*dum1
          rv=velf(ICFL,2,1)*dum1
          rw=velf(ICFL,3,1)*dum1

          gf2=SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw
          rva(ICFL,1)=
     &                rho(ICVA,1)*max(0.d0,xta(ICFL))
     &               +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &          +SFAREA(4,ICFL)*gf2 
          enddo
          enddo
        endif
      endif
!     
!----------------------------------
! --- <4.3 BC treatment for rva >--
!----------------------------------
!
!--------------------
! --- Euler 2 phase 
!--------------------
      if(ieul2ph>0) then
      else
!
!------------------
! --- Single phase 
!------------------
!
        do 340 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        IIMAT2=MAT_BCIDX(nb,2)
        if(.not.mat_cal(IIMAT)) cycle
        kd=kdbcnd(0,nb)
        kdv=kdbcnd(1,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd==kdilet.or.kd==kdstag) then
          if(ical_dens==4) then
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICVA=LVEDGE(1,ICFL)
            ICVB=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)*
     &          (SFAREA(1,ICFL)*vel(ICVB,1,1)
     &          +SFAREA(2,ICFL)*vel(ICVB,2,1)
     &          +SFAREA(3,ICFL)*vel(ICVB,3,1))
     &          +xta(ICFL)
            rva(ICFL,1)=
     &            +rho(ICVA,1)*max(0.d0,dum1)
     &            +rho(ICVB,1)*min(0.d0,dum1)
            enddo
          endif
        elseif(kd.eq.kdpres) then 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(ICV,1,1)
     &        +SFAREA(2,ICFL)*vel(ICV,2,1)
     &        +SFAREA(3,ICFL)*vel(ICV,3,1))+xta(ICFL)
          rva(ICFL,1)=rho(ICV,1)*max(0.d0,dum1)
     &               +rho(IDC,1)*min(0.d0,dum1)
          enddo
        elseif(kd==kdtchi) then
           do IBFL=IBFS,IBFE
           ICFL=LBC_SSF(IBFL)
           ICFP=LCYCSF(IBFL)
           ICVA=LVEDGE(1,ICFL)
           ICVB=LVEDGE(2,ICFL)
           dum1=SFAREA(4,ICFL)*
     &         (SFAREA(1,ICFL)*vel(ICVB,1,1)
     &         +SFAREA(2,ICFL)*vel(ICVB,2,1)
     &         +SFAREA(3,ICFL)*vel(ICVB,3,1))
     &         +xta(ICFL)
           rva(ICFL,1)=
     &             +rho(ICVA,1)*max(0.d0,dum1)
     &             +rho(ICVB,1)*min(0.d0,dum1)
           enddo
        elseif(kd.eq.kdolet.and.(openout(nb)==0)) then 
!          dum2=1.d0/1.4d0
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          ICV=ICVA
!          ICV=ICVB
!          dum3=(rho(ICVA,1)**1.4d0*prs(ICVB,1)/prs(ICVA,1))**dum2
          dum1=SFAREA(4,ICFL)*
     &        (SFAREA(1,ICFL)*vel(ICV,1,1)
     &        +SFAREA(2,ICFL)*vel(ICV,2,1)
     &        +SFAREA(3,ICFL)*vel(ICV,3,1))
     &        +xta(ICFL)
          rva(ICFL,1)=
     &            +rho(ICVA,1)*max(0.d0,dum1)
     &            +rho(ICVB,1)*min(0.d0,dum1)
!     &            +dum3*min(0.d0,dum1)
          enddo
        elseif(kd.eq.kdolet.and.
     &    (openout(nb)==8.or.openout(nb)==10)) then 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          ICV=ICVB
          dum1=SFAREA(4,ICFL)*
     &        (SFAREA(1,ICFL)*vel(ICV,1,1)
     &        +SFAREA(2,ICFL)*vel(ICV,2,1)
     &        +SFAREA(3,ICFL)*vel(ICV,3,1))
     &        +xta(ICFL)
          rva(ICFL,1)=
     &            +rho(ICVA,1)*max(0.d0,dum1)
     &            +rho(ICVB,1)*min(0.d0,dum1)
          enddo
        elseif(kd.eq.kdolet.and.(openout(nb)==7.or.openout(nb)==9)) then 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          ICV=ICVB
          dum1=SFAREA(4,ICFL)*
     &        (SFAREA(1,ICFL)*vel(ICV,1,1)
     &        +SFAREA(2,ICFL)*vel(ICV,2,1)
     &        +SFAREA(3,ICFL)*vel(ICV,3,1))
     &        +xta(ICFL)
          rva(ICFL,1)=
     &            +rho(ICVA,1)*max(0.d0,dum1)
     &            +rho(ICVB,1)*min(0.d0,dum1)
!          rva(ICFL)=
!     &             +rho(IDC)*min(0.d0,dum1)
          enddo
        elseif(kd.eq.kdolet.and.openout(nb).eq.1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          ICV=ICVB
          dum1=SFAREA(4,ICFL)*
     &        (SFAREA(1,ICFL)*vel(ICV,1,1)
     &        +SFAREA(2,ICFL)*vel(ICV,2,1)
     &        +SFAREA(3,ICFL)*vel(ICV,3,1))
     &        +xta(ICFL)
          rva(ICFL,1)=
     &            +rho(ICVA,1)*max(0.d0,dum1)
     &            +rho(ICVB,1)*min(0.d0,dum1)
          enddo
        elseif(kd.eq.kdolet.and.openout(nb).eq.2) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          ICV=ICVB
          dum1=SFAREA(4,ICFL)*
     &        (SFAREA(1,ICFL)*vel(ICV,1,1)
     &        +SFAREA(2,ICFL)*vel(ICV,2,1)
     &        +SFAREA(3,ICFL)*vel(ICV,3,1))
     &        +xta(ICFL)
          if(dum1<0.d0) then
            rva(ICFL,1)=0.d0
          else
            if(prs(ICV,1)>prsbnd(IBFL)) then
              rva(ICFL,1)=rho(ICV,1)*dum1
            else
              rva(ICFL,1)=0.d0
            endif
          endif
          enddo
        elseif(kd.eq.kdolet.and.(openout(nb)==4.or.openout(nb)==3)) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(ICV,1,1)
     &        +SFAREA(2,ICFL)*vel(ICV,2,1)
     &        +SFAREA(3,ICFL)*vel(ICV,3,1))+xta(ICFL)
          rva(ICFL,1)=
     &            +rho(ICV,1)*max(0.d0,dum1)
     &            +rho(IDC,1)*min(0.d0,dum1)
          enddo
        elseif(kd.eq.kdolet.and.openout(nb)==6) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
!          ICV=IDC
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
          rva(ICFL,1)=
     &            +rho(ICV,1)*max(0.d0,dum1)
     &            +rho(IDC,1)*min(0.d0,dum1)
          enddo
        elseif(kd.eq.kdsymm) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          rva(ICFL,1)=0.d0
          enddo
        elseif(kd.eq.kxnone.or.
     &         kd.eq.kdfire.or.
     &         kd.eq.kdcvd) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)!????2
!          dum1=xta(ICFL)!????2  cavitation
          rva(ICFL,1)=
     &            +rho(ICV,1)*max(0.d0,dum1)
     &            +rho(IDC,1)*min(0.d0,dum1)
!          rva(ICFL,1)=0.d0
          enddo
        elseif(kd.eq.kdprdc) then
          if(idis(nb)==0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            rvaicf=0.5d0*(rva(ICFL,1)-rva(ICFP,1))
            rva(ICFL,1)=  rvaicf
            rva(ICFP,1)= -rvaicf
            enddo
          elseif(idis(nb)==1) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICVA=LVEDGE(1,ICFL)
            ICVB=LVEDGE(2,ICFL)
!            
            dum1=SFAREA(4,ICFL)*0.5d0*
     &     (SFAREA(1,ICFL)*(vel(ICVA,1,1)+vel(ICVB,1,1))
     &     +SFAREA(2,ICFL)*(vel(ICVA,2,1)+vel(ICVB,2,1))
     &     +SFAREA(3,ICFL)*(vel(ICVA,3,1)+vel(ICVB,3,1)))
     &     +xta(ICFL)
            rva(ICFL,1)=
     &              +rho(ICVA,1)*max(0.d0,dum1)
     &              +rho(ICVB,1)*min(0.d0,dum1)
            enddo
          endif
        elseif(kd.eq.kdbuff) then
          do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            rvaicf=0.5d0*(rva(ICFL,1)-rva(ICFP,1))
            rva(ICFL,1)=  rvaicf
            rva(ICFP,1)= -rvaicf
          enddo
        elseif(kd.eq.kdshutr) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          if(SHUTFL(IBFL)==0) then
            rvaicf=0.5d0*(rva(ICFL,1)-rva(ICFP,1))
            rva(ICFL,1)=  rvaicf
            rva(ICFP,1)= -rvaicf
          else
            rva(ICFL,1)=0.D0
            rva(ICFP,1)=0.D0
          endif
          enddo

        elseif(kd.eq.kdintr) then 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          rva(ICFL,1)=0.d0
          rva(ICFP,1)=0.d0
          enddo
        elseif(kd==kdovst) then
          do IBFL=IBFS,IBFE 
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LCYCSF(IBFL)
!
          dum1=SFAREA(4,ICFL)*
     &        (SFAREA(1,ICFL)*(vel(IDC,1,1))
     &        +SFAREA(2,ICFL)*(vel(IDC,2,1))
     &        +SFAREA(3,ICFL)*(vel(IDC,3,1)))
     &        +xta(ICFL)
!
!         dum1=SFAREA(4,ICFL)*0.5d0*
!     &        (SFAREA(1,ICFL)*(vel(IDC,1,1)+vel(ICV,1,1))
!     &        +SFAREA(2,ICFL)*(vel(IDC,2,1)+vel(ICV,2,1))
!     &        +SFAREA(3,ICFL)*(vel(IDC,3,1)+vel(ICV,3,1))
!     &        )
!     &       +xta(ICFL)
!
          rva(ICFL,1)=
     &             +rho(ICV,1)*max(0.d0,dum1)
     &            +rho(ICVP,1)*min(0.d0,dum1)
          enddo
        elseif(kd==kdsld.and.idis(nb)==1) then
          dum_2=1.d0 
          if(iLBF_P==2) dum_2=0.d0
          do 2200 ISLD=1,2
          IF(ISLD==1) THEN
            ISLD2=2
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=IBFS+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
            elseif(idis(nb)>=1) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_pair(nb)
            endif
          ELSE
            ISLD2=1
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2
              IBFE=LBC_INDEX(nb)
            elseif(idis(nb)>=1) then
              IBFS=LBC_pair(nb)+1
              IBFE=LBC_INDEX(nb)
            endif
          ENDIF
!
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          th(ISLD)=rot_ang(IMATS(ISLD))
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2
     &              +unit(2,ISLD2)**2
     &              +unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
          th(ISLD2)=rot_ang(IMATS(ISLD2))
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
          org_x1=begin(1,IMATS(ISLD))*costh1
     &          -begin(2,IMATS(ISLD))*sinth1
          org_y1=begin(1,IMATS(ISLD))*sinth1
     &          +begin(2,IMATS(ISLD))*costh1
          org_x2=begin(1,IMATS(ISLD2))*costh2
     &          -begin(2,IMATS(ISLD2))*sinth2
          org_y2=begin(1,IMATS(ISLD2))*sinth2
     &          +begin(2,IMATS(ISLD2))*costh2
!                                
          do 2210 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(1,ICFP)
          ICVB1=LVEDGE(2,ICFL)  !??????????????????????????
          ICVBO=LVEDGE(1,ICFO)
!          
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
!
          dum1=CVCENT(1,ICVA)
          dum2=CVCENT(2,ICVA)
          dum3=CVCENT(3,ICVA)
          avsf=SFAREA(4,ICFL)
!----------------------------------------------
          dxo=CVCENT(1,ICVB1)
          dyo=CVCENT(2,ICVB1)
          dzo=CVCENT(3,ICVB1)
          dxo=dxo-dum1
          dyo=dyo-dum2
          dzo=dzo-dum3
          dll=dsqrt(dxo*dxo+dyo*dyo+dzo*dzo+SML)
          dxo=dxo/dll
          dyo=dyo/dll
          dzo=dzo/dll          
          dlvect=abs(dxo*SFAREA(1,ICFL)
     &              +dyo*SFAREA(2,ICFL)
     &              +dzo*SFAREA(3,ICFL))
!----------------------------------------------
          dum1=prs(ICVB,1)*wi1+prs(ICVBO,1)*wi2
          gf1=avsf*rdeltt*(dum1-prs(ICVA,1))/(dll*dlvect)
!
          sf11=rbb(1,1,ISLD)*SFAREA(1,ICFL)
     &        +rbb(1,2,ISLD)*SFAREA(2,ICFL)
     &        +rbb(1,3,ISLD)*SFAREA(3,ICFL)
          sf12=rbb(2,1,ISLD)*SFAREA(1,ICFL)
     &        +rbb(2,2,ISLD)*SFAREA(2,ICFL)
     &        +rbb(2,3,ISLD)*SFAREA(3,ICFL)
          sf13=rbb(3,1,ISLD)*SFAREA(1,ICFL)
     &        +rbb(3,2,ISLD)*SFAREA(2,ICFL)
     &        +rbb(3,3,ISLD)*SFAREA(3,ICFL)
!----------------------------------
! --- ICVA
          vell(1,1)=rbb(1,1,ISLD)*vel(ICVA,1,1)
     &             +rbb(1,2,ISLD)*vel(ICVA,2,1)
     &             +rbb(1,3,ISLD)*vel(ICVA,3,1)
          vell(1,2)=rbb(2,1,ISLD)*vel(ICVA,1,1)
     &             +rbb(2,2,ISLD)*vel(ICVA,2,1)
     &             +rbb(2,3,ISLD)*vel(ICVA,3,1)
          vell(1,3)=rbb(3,1,ISLD)*vel(ICVA,1,1)
     &             +rbb(3,2,ISLD)*vel(ICVA,2,1)
     &             +rbb(3,3,ISLD)*vel(ICVA,3,1)
! --- ICVB
          vell(2,1)=rbb(1,1,ISLD2)*vel(ICVB,1,1)
     &             +rbb(1,2,ISLD2)*vel(ICVB,2,1)
     &             +rbb(1,3,ISLD2)*vel(ICVB,3,1)
          vell(2,2)=rbb(2,1,ISLD2)*vel(ICVB,1,1)
     &             +rbb(2,2,ISLD2)*vel(ICVB,2,1)
     &             +rbb(2,3,ISLD2)*vel(ICVB,3,1)
          vell(2,3)=rbb(3,1,ISLD2)*vel(ICVB,1,1)
     &             +rbb(3,2,ISLD2)*vel(ICVB,2,1)
     &             +rbb(3,3,ISLD2)*vel(ICVB,3,1)
! --- ICVBO
          vell(3,1)=rbb(1,1,ISLD2)*vel(ICVBO,1,1)
     &             +rbb(1,2,ISLD2)*vel(ICVBO,2,1)
     &             +rbb(1,3,ISLD2)*vel(ICVBO,3,1)
          vell(3,2)=rbb(2,1,ISLD2)*vel(ICVBO,1,1)
     &             +rbb(2,2,ISLD2)*vel(ICVBO,2,1)
     &             +rbb(2,3,ISLD2)*vel(ICVBO,3,1)
          vell(3,3)=rbb(3,1,ISLD2)*vel(ICVBO,1,1)
     &             +rbb(3,2,ISLD2)*vel(ICVBO,2,1)
     &             +rbb(3,3,ISLD2)*vel(ICVBO,3,1)
!
          grdff(1,1)=rbb(1,1,ISLD)*grdc(ICVA,1,1)
     &              +rbb(1,2,ISLD)*grdc(ICVA,2,1)
     &              +rbb(1,3,ISLD)*grdc(ICVA,3,1)
          grdff(1,2)=rbb(2,1,ISLD)*grdc(ICVA,1,1)
     &              +rbb(2,2,ISLD)*grdc(ICVA,2,1)
     &              +rbb(2,3,ISLD)*grdc(ICVA,3,1)
          grdff(1,3)=rbb(3,1,ISLD)*grdc(ICVA,1,1)
     &              +rbb(3,2,ISLD)*grdc(ICVA,2,1)
     &              +rbb(3,3,ISLD)*grdc(ICVA,3,1)
!
          grdff(2,1)=rbb(1,1,ISLD2)*grdc(ICVB,1,1)
     &              +rbb(1,2,ISLD2)*grdc(ICVB,2,1)
     &              +rbb(1,3,ISLD2)*grdc(ICVB,3,1)
          grdff(2,2)=rbb(2,1,ISLD2)*grdc(ICVB,1,1)
     &              +rbb(2,2,ISLD2)*grdc(ICVB,2,1)
     &              +rbb(2,3,ISLD2)*grdc(ICVB,3,1)
          grdff(2,3)=rbb(3,1,ISLD2)*grdc(ICVB,1,1)
     &              +rbb(3,2,ISLD2)*grdc(ICVB,2,1)
     &              +rbb(3,3,ISLD2)*grdc(ICVB,3,1)
!
          grdff(3,1)=rbb(1,1,ISLD2)*grdc(ICVBO,1,1)
     &              +rbb(1,2,ISLD2)*grdc(ICVBO,2,1)
     &              +rbb(1,3,ISLD2)*grdc(ICVBO,3,1)
          grdff(3,2)=rbb(2,1,ISLD2)*grdc(ICVBO,1,1)
     &              +rbb(2,2,ISLD2)*grdc(ICVBO,2,1)
     &              +rbb(2,3,ISLD2)*grdc(ICVBO,3,1)
          grdff(3,3)=rbb(3,1,ISLD2)*grdc(ICVBO,1,1)
     &              +rbb(3,2,ISLD2)*grdc(ICVBO,2,1)
     &              +rbb(3,3,ISLD2)*grdc(ICVBO,3,1)
!
          ru1=(rho(ICVA,1) *vell(1,1)+rdeltt*grdff(1,1))
          ru2=(rho(ICVB,1) *vell(2,1)+rdeltt*grdff(2,1))
          ru3=(rho(ICVBO,1)*vell(3,1)+rdeltt*grdff(3,1))
          ru=0.5d0*sf11*(ru1+ru2*wi1+ru3*wi2)
!
          rv1=(rho(ICVA,1) *vell(1,2)+rdeltt*grdff(1,2))
          rv2=(rho(ICVB,1) *vell(2,2)+rdeltt*grdff(2,2))
          rv3=(rho(ICVBO,1)*vell(3,2)+rdeltt*grdff(3,2))
          rv=0.5d0*sf12*(rv1+rv2*wi1+rv3*wi2)
!
          rw1=(rho(ICVA,1) *vell(1,3)+rdeltt*grdff(1,3))
          rw2=(rho(ICVB,1) *vell(2,3)+rdeltt*grdff(2,3))
          rw3=(rho(ICVBO,1)*vell(3,3)+rdeltt*grdff(3,3))
          rw=0.5d0*sf13*(rw1+rw2*wi1+rw3*wi2)
          gf2=avsf*(ru+rv+rw)
!-------------------------------------------------------------------------
          rvaicf=rho(ICVA,1)*max(0.d0,xta(ICFL))
     &          +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &          +gf2-dum_2*(gf1)
          rva(ICFL,1)=rvaicf
!
 2210     enddo
 2200     enddo

        elseif(kd==kdsld.and.idis(nb)==2) then
          dum_2=0.d0 
          if(iLBF_P==2) dum_2=0.d0
          do ISLD=1,2
          IF(ISLD==1) THEN
            ISLD2=2
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=IBFS+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
            elseif(idis(nb)>=1) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_pair(nb)
            endif
          ELSE
            ISLD2=1
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2
              IBFE=LBC_INDEX(nb)
            elseif(idis(nb)>=1) then
              IBFS=LBC_pair(nb)+1
              IBFE=LBC_INDEX(nb)
            endif
          ENDIF
!
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          th(ISLD)=rot_ang(IMATS(ISLD))
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2
     &              +unit(2,ISLD2)**2
     &              +unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
          org_x1=begin(1,IMATS(ISLD))*costh1
     &          -begin(2,IMATS(ISLD))*sinth1
          org_y1=begin(1,IMATS(ISLD))*sinth1
     &          +begin(2,IMATS(ISLD))*costh1
          org_x2=begin(1,IMATS(ISLD2))*costh2
     &          -begin(2,IMATS(ISLD2))*sinth2
          org_y2=begin(1,IMATS(ISLD2))*sinth2
     &          +begin(2,IMATS(ISLD2))*costh2
!                                
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(1,ICFP)
          ICVB1=LVEDGE(2,ICFL)  !??????????????????????????
          ICVBO=LVEDGE(1,ICFO)
!          
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
!
          th(ISLD2)=OPPANG(IBFL)
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          dum1=CVCENT(1,ICVA)
          dum2=CVCENT(2,ICVA)
          dum3=CVCENT(3,ICVA)
          avsf=SFAREA(4,ICFL)
!----------------------------------------------
          dxo=CVCENT(1,ICVB1)
          dyo=CVCENT(2,ICVB1)
          dzo=CVCENT(3,ICVB1)
          dxo=dxo-dum1
          dyo=dyo-dum2
          dzo=dzo-dum3
          dll=dsqrt(dxo*dxo+dyo*dyo+dzo*dzo+SML)
          dxo=dxo/dll
          dyo=dyo/dll
          dzo=dzo/dll          
          dlvect=abs(dxo*SFAREA(1,ICFL)
     &              +dyo*SFAREA(2,ICFL)
     &              +dzo*SFAREA(3,ICFL))
!----------------------------------------------
          dum1=prs(ICVB,1)*wi1+prs(ICVBO,1)*wi2
          gf1=avsf*rdeltt*(dum1-prs(ICVA,1))/(dll*dlvect)
!
          sf11=rbb(1,1,ISLD)*SFAREA(1,ICFL)
     &        +rbb(1,2,ISLD)*SFAREA(2,ICFL)
     &        +rbb(1,3,ISLD)*SFAREA(3,ICFL)
          sf12=rbb(2,1,ISLD)*SFAREA(1,ICFL)
     &        +rbb(2,2,ISLD)*SFAREA(2,ICFL)
     &        +rbb(2,3,ISLD)*SFAREA(3,ICFL)
          sf13=rbb(3,1,ISLD)*SFAREA(1,ICFL)
     &        +rbb(3,2,ISLD)*SFAREA(2,ICFL)
     &        +rbb(3,3,ISLD)*SFAREA(3,ICFL)
!----------------------------------
! --- ICVA
          vell(1,1)=rbb(1,1,ISLD)*vel(ICVA,1,1)
     &             +rbb(1,2,ISLD)*vel(ICVA,2,1)
     &             +rbb(1,3,ISLD)*vel(ICVA,3,1)
          vell(1,2)=rbb(2,1,ISLD)*vel(ICVA,1,1)
     &             +rbb(2,2,ISLD)*vel(ICVA,2,1)
     &             +rbb(2,3,ISLD)*vel(ICVA,3,1)
          vell(1,3)=rbb(3,1,ISLD)*vel(ICVA,1,1)
     &             +rbb(3,2,ISLD)*vel(ICVA,2,1)
     &             +rbb(3,3,ISLD)*vel(ICVA,3,1)
! --- ICVB
          vell(2,1)=rbb(1,1,ISLD2)*vel(ICVB,1,1)
     &             +rbb(1,2,ISLD2)*vel(ICVB,2,1)
     &             +rbb(1,3,ISLD2)*vel(ICVB,3,1)
          vell(2,2)=rbb(2,1,ISLD2)*vel(ICVB,1,1)
     &             +rbb(2,2,ISLD2)*vel(ICVB,2,1)
     &             +rbb(2,3,ISLD2)*vel(ICVB,3,1)
          vell(2,3)=rbb(3,1,ISLD2)*vel(ICVB,1,1)
     &             +rbb(3,2,ISLD2)*vel(ICVB,2,1)
     &             +rbb(3,3,ISLD2)*vel(ICVB,3,1)
! --- ICVBO
          vell(3,1)=rbb(1,1,ISLD2)*vel(ICVBO,1,1)
     &             +rbb(1,2,ISLD2)*vel(ICVBO,2,1)
     &             +rbb(1,3,ISLD2)*vel(ICVBO,3,1)
          vell(3,2)=rbb(2,1,ISLD2)*vel(ICVBO,1,1)
     &             +rbb(2,2,ISLD2)*vel(ICVBO,2,1)
     &             +rbb(2,3,ISLD2)*vel(ICVBO,3,1)
          vell(3,3)=rbb(3,1,ISLD2)*vel(ICVBO,1,1)
     &             +rbb(3,2,ISLD2)*vel(ICVBO,2,1)
     &             +rbb(3,3,ISLD2)*vel(ICVBO,3,1)
!
          grdff(1,1)=rbb(1,1,ISLD)*grdc(ICVA,1,1)
     &              +rbb(1,2,ISLD)*grdc(ICVA,2,1)
     &              +rbb(1,3,ISLD)*grdc(ICVA,3,1)
          grdff(1,2)=rbb(2,1,ISLD)*grdc(ICVA,1,1)
     &              +rbb(2,2,ISLD)*grdc(ICVA,2,1)
     &              +rbb(2,3,ISLD)*grdc(ICVA,3,1)
          grdff(1,3)=rbb(3,1,ISLD)*grdc(ICVA,1,1)
     &              +rbb(3,2,ISLD)*grdc(ICVA,2,1)
     &              +rbb(3,3,ISLD)*grdc(ICVA,3,1)
!
          grdff(2,1)=rbb(1,1,ISLD2)*grdc(ICVB,1,1)
     &              +rbb(1,2,ISLD2)*grdc(ICVB,2,1)
     &              +rbb(1,3,ISLD2)*grdc(ICVB,3,1)
          grdff(2,2)=rbb(2,1,ISLD2)*grdc(ICVB,1,1)
     &              +rbb(2,2,ISLD2)*grdc(ICVB,2,1)
     &              +rbb(2,3,ISLD2)*grdc(ICVB,3,1)
          grdff(2,3)=rbb(3,1,ISLD2)*grdc(ICVB,1,1)
     &              +rbb(3,2,ISLD2)*grdc(ICVB,2,1)
     &              +rbb(3,3,ISLD2)*grdc(ICVB,3,1)
!
          grdff(3,1)=rbb(1,1,ISLD2)*grdc(ICVBO,1,1)
     &              +rbb(1,2,ISLD2)*grdc(ICVBO,2,1)
     &              +rbb(1,3,ISLD2)*grdc(ICVBO,3,1)
          grdff(3,2)=rbb(2,1,ISLD2)*grdc(ICVBO,1,1)
     &              +rbb(2,2,ISLD2)*grdc(ICVBO,2,1)
     &              +rbb(2,3,ISLD2)*grdc(ICVBO,3,1)
          grdff(3,3)=rbb(3,1,ISLD2)*grdc(ICVBO,1,1)
     &              +rbb(3,2,ISLD2)*grdc(ICVBO,2,1)
     &              +rbb(3,3,ISLD2)*grdc(ICVBO,3,1)
!
          ru1=(rho(ICVA,1) *vell(1,1)+rdeltt*grdff(1,1))
          ru2=(rho(ICVB,1) *vell(2,1)+rdeltt*grdff(2,1))
          ru3=(rho(ICVBO,1)*vell(3,1)+rdeltt*grdff(3,1))
          ru=0.5d0*sf11*(ru1+ru2*wi1+ru3*wi2)
!
          rv1=(rho(ICVA,1) *vell(1,2)+rdeltt*grdff(1,2))
          rv2=(rho(ICVB,1) *vell(2,2)+rdeltt*grdff(2,2))
          rv3=(rho(ICVBO,1)*vell(3,2)+rdeltt*grdff(3,2))
          rv=0.5d0*sf12*(rv1+rv2*wi1+rv3*wi2)
!
          rw1=(rho(ICVA,1) *vell(1,3)+rdeltt*grdff(1,3))
          rw2=(rho(ICVB,1) *vell(2,3)+rdeltt*grdff(2,3))
          rw3=(rho(ICVBO,1)*vell(3,3)+rdeltt*grdff(3,3))
          rw=0.5d0*sf13*(rw1+rw2*wi1+rw3*wi2)
          gf2=avsf*(ru+rv+rw)
!-------------------------------------------------------------------------
          rvaicf=rho(ICVA,1)*max(0.d0,xta(ICFL))
     &          +rho(ICVB,1)*min(0.d0,xta(ICFL))
     &          +gf2-dum_2*(gf1)
          rva(ICFL,1)=rvaicf
          enddo
          enddo
        endif  
 340    enddo
      endif
!

      return
!
 9999 continue
!
      if(my_rank.eq.ROOT) write(ifle,*) '(rva_dmin)'
      ierror=1
!
      end subroutine rva_admin
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hys_admin(ismpl,iphs,deltt,iter,p_grdc,time,ictl,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rva,rmd,rds,rmut,rho,vel,prs,pp0,wdot,rmu,
     &  utau,frstcv,
     &  EVAP_Y,EVAP_H,
     &  hhh,hhs,cps,cp,rmdc,dsclt,dscly,cofd,
     &  grdc,diagy,diagh,kdbt,kdby,rvd,
!     &  aks,vctr,mdot,ccc,
     &  aks,mdot,ccc,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  SUFBND,SDOT,LCYCOLD,wifsld,OPPANG,FIELD_U,
     &  tmpbnd,yysbnd,htcbnd,mtcbnd,radbnd,HTFLUX,
     &  tmp,yys,MASS_I,
     &  RadHeat,RadHeatFlux,pabsorb,pscatter,
     &  itery,repsy,aepsy,erry,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       rmdc    <= cr
!       dsclt   <= rvx(1,:)
!       dscly   <= rvx(2,:)
!       cofd    <= rvx(3,:)
!       diagy   <= diag
!       diagh   <= rmx
!ALPN=aks(:,iaph(1),1)
!ALPO=aks(:,iaph(1),2)
!vof=aks(:,ivof,1)
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only      : ifle,ifll
      use module_model   ,only: idrdp,mach0,comp,ical_t,ical_reac,
     &                          PEFC,PAFC,MCFC,SOFC,AFC,MDFC,iheat,
     &                          ical_vect,ical_dens,u_func,drag_coef
      use module_flags   ,only: intgty,eulere
      use module_material,only: icnvty,lmtrty,ivscty,cal,
     &                          nflud,nofld,lclsd,nsold,nosld,
     &                          relaxh,tsat,hgsat,hlsat,BGCOMP,
     &                          idarcy2,porosty,porous,ical_porous,
     &                          N_pors,
     &                          relaxys,
     &                         ical_sld,rotati,ishaft,end,begin
      use module_param   ,only: yslw
      USE module_usersub, ONLY: src_t,src_r,usrno,usryes,src_fire
      use module_species, only: sw,acpk,c_p,enthalpy_ref,hform,spcnam,
     &                          wm,r_wm,Ri,gascns,h_ref,enthalpy_r,
     &                          acpk2,I_H2O,I_PROTON,spcno,act
      use module_chemreac,only: nneq,ical_suf
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd,surfreac,
     &                          phs_idx,phs_com,heat
      use module_time,    only: steady,i_steady
      use module_simple,  only: simeer
      use module_scalar  ,only: ivof,ike,iaph,ical_cavi,
     &                          ical_FC,ical_s,ical_prp
      use module_vof     ,only: ical_vof,intervof,LG,LS,evap_hh,
     &                          change_phase
      use module_Euler2ph,only: ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_metrix,only  : rdsc =>d2work1,vctr
      use module_metrix,only  : dydt =>d2work2
      use module_metrix,only  : rvds =>d2work3
      use module_metrix,only  : dhdt =>d1work1
      use module_metrix,only  : yys_tmp=>W1K11
      use module_metrix,only  : ys,rcomp,eva_comp,sinj_comp,erryys,ahk
      use module_FUEL  ,only  : No_Mem,vap_no
      use module_rad
      use module_radsxf, only : QSUM,EBcomm,NRD,NRD0
      use module_radsxf, only : GWAbsorb    !,PAbsorb,PScatter
      use module_usersub,only : src_rad,usryes
      use module_model,   only : ical_prt
      use module_initial ,only : p0
      use module_metrix,only  : multiR
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
      real*8 ,intent(inout) :: rmd       (   MXALLCV)
      real*8 ,intent(inout) :: rds       (   MXALLCV,MXcomp)
      real*8 ,intent(in)    :: rmut      (   MXALLCV)
      REAL*8 ,INTENT(IN)    :: RMU   (        MXALLCV)
      real*8 ,intent(inout)    :: rho       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3  )
      real*8 ,intent(in)    :: prs       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: pp0       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: wdot  (       MXALLCV,MXcomp)

      real*8 ,intent(in)    :: EVAP_Y(       MXALLCV_P,MXCOMP)
      real*8 ,intent(in)    :: EVAP_H(       MXALLCV_P,2)
      REAL*8 ,INTENT(INOUT) :: hhh   (       MXALLCV  ,2)
      REAL*8 ,INTENT(INOUT) :: hhs   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cps   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp    (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: rmdc  (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: dsclt (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: dscly (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cofd  (       MXALLCV)
      REAL*8 ,INTENT(IN)    :: SUFBND(       MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(INOUT) :: SDOT  (       MXSSFBC_SUF,MXCOMPALL)
      real*8 ,intent(in)    :: tmpbnd(       MXssfbc)
      real*8 ,intent(in)    :: yysbnd(       MXssfbc,MXcomp)
      real*8 ,intent(in)    :: htcbnd    (   MXssfbc)
      real*8 ,intent(inout) :: HTFLUX    (   MXssfbc)
      real*8 ,intent(in)    :: mtcbnd    (   MXssfbc)
      real*8 ,intent(in)    :: radbnd    (   MXssfbc)
      real*8 ,intent(inout) :: tmp       (   MXALLCV  ,2)
      real*8 ,intent(inout) :: yys       (   MXALLCV,MXcomp,2)
      real*8 ,intent(inout) :: grdc      (   MXALLCV,3,3)
      real*8 ,intent(inout) :: aks       (   MXALLCVR,mxrans)
      real*8 ,intent(inout) :: diagy     (   MXALLCV)
      real*8 ,intent(inout) :: diagh     (   MXALLCV)
      real*8 ,intent(inout) :: rvd       (   MXCVFAC,3)
      INTEGER,INTENT(INOUT) :: kdbt      (   MXCVFAC)
      INTEGER,INTENT(INOUT) :: kdby      (   MXCVFAC)
      integer,intent(out)   :: ierror,itery(0:MXcomp)
      real*8 ,intent(out)   :: repsy(0:MXcomp),aepsy(0:MXcomp)
      real*8 ,intent(out)   :: erry(0:MXcomp)
      integer,intent(in)    :: LCYCOLD   (   MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld    (   MXSSFBC_SLD)!OPPANG
      real*8 ,intent(in)    :: OPPANG    (   MXSSFBC_SLD)
!      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)      
      real*8 ,intent(in)    :: MASS_I   (    MXALLCV2)
      REAL*8 ,INTENT(IN)    :: mdot(MXCV,3)
      REAL*8 ,INTENT(IN)   :: CCC   (        MXALLCV)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
      real*8 ,intent(inout) :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT) :: RadHeat    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT) :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT) :: PScatter   (MXALLCV_RAD)
!
      REAL*8 ,INTENT(IN)    :: FRSTCV(         MXSSFBC)
      REAL*8 ,INTENT(IN)    :: UTAU  (       0:MXSSFBC)
! 
! --- [local entities] 
! 
      real*8  :: unstdy_h=1.d0,dumm=1.d0
      real*8  :: duma(1),cpp,dhh,vol
      real*8  :: sinj_r,sinj_u,sinj_v,sinj_w,sinj_t,ad_comp
      real*8  :: eva_mass,unitmass,T_WALL
      real*8  :: sumcop,heat1,eva_heat,eva_T
      real*8  :: u,v,w,T,dens,suu,suv,suw,p
      real*8  :: vdlt,dum1,dum2,dum3,wi1,wi2
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,ndf,nq,IMAX,no
      integer :: ICOM,IMD,ICH,IFLD,ICTP,ierr1=0
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: IMODE,IMAT_U,IBFS,IBFE,IBFL
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP,ISEC
      integer :: ipcom,icoms,icome,iph
      logical :: viscous,noviscs,inoheat=.false.
      real*8  :: grx,gry,grz,qvof,evapmdl,dumx,dum,dum_typ
      character*80 :: BC_name
      real*8  :: dhdtmax
      integer :: Mark_rad_done,IIMAT2,Icomp,iters
      real*8  :: twheat(nbcnd,3),trad,rr,ppp,shigma=5.67e-8,vv1,vv2,vv3
      real*8  :: aayy(ncomp),rdeltt,vr(3),radius(3),shaftn(3),velref(3)
      integer :: ndiag,idum
!
      ALLOCATE(rdsc(MXALLCV,MXcomp),stat=ierr1)
      ALLOCATE(dhdt(MXALLCV),stat=ierr1)
      ALLOCATE(dydt(MXALLCV,MXcomp),stat=ierr1)
      ALLOCATE(rvds(MXCVFAC,MXcomp),stat=ierr1)
      if(ierr1/=0) call FFRABORT(1,'allocating error in hys_admin')
!
!----------------------------------------------------------------
! --- 
!----------------------
!
      diagy=0.d0
      diagh=0.d0
      kdbt=0.d0
      kdby=0.d0
      do ICOM=1,ncomp
      hhs(:,ICOM)=0.d0
      cps(:,ICOM)=0.d0
      dydt(:,ICOM)=0.d0
      rdsc(:,ICOM)=0.d0
      enddo
      cp=0.d0
      dsclt=0.d0
      dscly=0.d0
      do I=1,3
      rvd(:,I)=0.d0
      enddo
      dhdt=0.d0
      vdlt=1.d0/deltt
      rdeltt=1.d0/deltt
      if(i_steady==3) then
        rdeltt=1.d0       !1212
      endif
!2222
      if(.not.sw) then
       if(iphs==1) then
        do ICOM=1,ncomp
        do n=1,5
          ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk(0,ICOM)
        enddo
       else
        do ICOM=1,ncomp
        do n=1,5
          ahk(n,ICOM)=acpk2(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk2(0,ICOM)
        enddo
       endif
      endif
!-------------------------------
! --- ZONE or MC
!-------------------------------
      Mark_Rad_Done=0
!
      IF(ical_t.AND.radflag==1.and.radmodel(1:3).NE.'FVM') THEN
        Mark_Rad_Done=1
	IF(MAT_CVEXT(NMAT).NE.NRD0) THEN
	  WRITE(ifle,*) 'Error in Radiation, NCV != NRD0'
          WRITE(ifle,*) MAT_CVEXT(NMAT), '!=', NRD0
	  call FFRABORT(1,'ERR: Error in Radiation')
	ENDIF
!
! --- RadHeat :	J/S
!
        CALL GetHeatFlux(iter,NPE,nbcnd,NMAT,nallcv,
     &    tmp(:,1),MAT_CV,RadHeatFlux,RadHeat)
        IF(NPE.GT.1) THEN
          CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,RadHeat)
	  CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,RadHeatFlux)
	  CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,EBcomm)
        ENDIF
!
	twheat = 0.d0
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
	  twheat(nb,3) =twheat(nb,3)+SFAREA(4,ICFL)
	  twheat(nb,1) = twheat(nb,1)+RadHeat(IDC) 
	enddo
!
	if(twheat(nb,3).gt.0.d0) then
          twheat(nb,1)=twheat(nb,1)/twheat(nb,3)
	endif
	ENDDO
!	
!	if(my_rank.eq.0) write(ifll,*) 'Radiation Calculation ...'
!       
	if(src_rad.eq.usryes) then
          do IIMAT=1,NMAT   !ICN=1,NCV
            if(.not.mat_cal(IIMAT)) cycle
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
	    do ICVL=ICVS,ICVE
	      ICV=MAT_CV(ICVL)
	      do Icomp=1,ncomp
              aayy(Icomp)=yys(ICVL,Icomp,1)
	      end do
              call user_rad_pro(ICV,ncomp,spcno,CVCENT(:,ICVL),
     &	      tmp(ICVL,1),prs(ICVL,1),rho(ICVL,1),aayy,
     &	      GWAbsorb(ICVL),PAbsorb(ICVL),PScatter(ICVL))
              do Icomp=1,ncomp
		yys(ICVL,Icomp,1)=aayy(Icomp)
              enddo
	      enddo
	    enddo
	endif
!
      ENDIF
!-----------------------------------------------------------------
! --- < 1.1 calculate enthalpy & specific heat from temperature>--
!-----------------------------------------------------------------
!
      if(ical_t.and.ical_prp) then 
        if(sw) then
          call cal_t2cp(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &          tmp(:,1),yys(:,:,1),prs(:,1),rho(:,1),
     &          rmu,rmd,cps,cp,rmdc)
          if(ical_dens==4) then
          else
            call cal_t2h(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &          tmp(:,1),yys(:,:,1),hhs,cofd,rho(:,1),prs(:,1))
!     &          tmp(:,1),yys(:,:,1),hhs,hhh(:,1),rho(:,1),prs(:,1))
          endif
        else
          if(ical_vof==1) then
          elseif(ical_cavi==1) then
            call set_trnsprp_t2h_cavi(0,
     &      LVEDGE,LBC_SSF,LCYCSF,
     &      MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,wifsld,LCYCOLD,
     &      diagy,hhs,cps,cp,cofd,
     &      tmp(:,1),yys(:,:,1),aks,rho(:,1),dhdt,rmd,rds,rmdc)
!     &      tmp(:,1),yys(:,:,1),aks,rho(:,1),hhh(:,1),rmd,rds,rmdc)
            dhdt=0.d0
            diagy=0.d0
          else
            call cal_t2hcp(iphs,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &          tmp(:,1),yys(:,:,1),hhs,cps,cofd,cp,rmdc,prs)
!     &          tmp(:,1),yys(:,:,1),hhs,cps,hhh(:,1),cp,rmdc,prs)
          endif
        endif
      endif
!2222
!-------------------------------------------------
! --- < 1.2 time difference term & source term >--
!-------------------------------------------------
!      call src_setxxx(1,'tmp',rho(:,1),tmp(:,1),dhdt,diagh)
!
!      call src_setxxx(ncomp,'yys',rho(:,1),yys(:,:,1),dydt,diagy)
!--------------------------------------------
! --- User rho source term subroutine 
!--------------------------------------------
!
      if(src_r.eq.usryes) then
         do 155 IIMAT=1,NMAT   !ICN=1,NCV
         if(.not.mat_cal(IIMAT)) goto 155
         IMAT=MAT_NO(IIMAT)
         ICVS=MAT_CVEXT(IIMAT-1)+1
         ICVE=MAT_CVEXT(IIMAT)
         if(IMAT.gt.0) then
           IMAT_U=nofld(IMAT)
           DO 156 ICVL=ICVS,ICVE
           ICV=MAT_CV(ICVL)
           if(NPE.gt.1) ICV=NOD_IDg(ICV)
           u=vel(ICVL,1)
           v=vel(ICVL,2)
           w=vel(ICVL,3)
           T=tmp(ICVL,1)
           p=prs(ICVL,1) 
           do icom=1,ncomp
           rcomp(icom)=yys(ICVL,icom,1)
           sinj_comp(icom)=yys(ICVL,icom,1)
           enddo
           dens=rho(ICVL,1)
           sinj_u=0.d0
           sinj_v=0.d0
           sinj_w=0.d0
           sinj_r=0.d0
           sinj_t=tmp(ICVL,1)
           vol=CVVOLM(ICVL)
           grx=CVCENT(1,ICVL)
           gry=CVCENT(2,ICVL)
           grz=CVCENT(3,ICVL)
           call user_src_r(2,ictl,
     &     deltt,iter,time,ICV,IMAT_U,iphs,ncomp,grx,gry,grz,
     &     u,v,w,T,p,dens,rcomp,vol,
     &     sinj_r,sinj_u,sinj_v,sinj_w,sinj_t,sinj_comp)
!
           sumcop=SML
           do icom=1,ncomp
             sumcop=sumcop+sinj_comp(icom)
           enddo
           sinj_comp(:)=sinj_comp(:)/sumcop
           if(sinj_r.gt.SML) then  !zhang-cvd
             do icom=1,ncomp   !sinj_r=[kg/(s*m^3)] [m^3/s]
             dydt(ICVL,ICOM)=dydt(ICVL,ICOM)+
     &                       sinj_r*(sinj_comp(icom)-0.d0*rcomp(icom))
             enddo             !dhdt=[J/s/m^3];cp=J/K/kg
!
             if(.not.sw) then
               dhh=0.d0
               do icom=1,ncomp
               dhh=dhh+sinj_comp(ICOM)*(((((
     &         ahk(5,ICOM) *sinj_t
     &        +ahk(4,ICOM))*sinj_t
     &        +ahk(3,ICOM))*sinj_t
     &        +ahk(2,ICOM))*sinj_t
     &        +ahk(1,ICOM))*sinj_t
     &        +ahk(0,ICOM))
              enddo
              dhh=sinj_r*dhh
            else  !zhang-cvd:
              dhh=sinj_r*enthalpy_r(ncomp,sinj_t,T,sinj_comp)
!             dhh=sinj_r*(enthalpy_ref(ncomp,sinj_t,sinj_comp)
!     &             -0.d0*enthalpy_ref(ncomp,T,rcomp))
!              hhh(ICVL,1)=enthalpy_ref(ncomp,sinj_t,sinj_comp)
            endif
            dhdt(ICVL)=dhdt(ICVL)+dhh  !zhang-cvd 
            diagh(ICVL)=diagh(ICVL)+sinj_r
            diagy(ICVL)=diagy(ICVL)+sinj_r
          elseif(sinj_r.lt.-SML) then
             do icom=1,ncomp
             dydt(ICVL,ICOM)=dydt(ICVL,ICOM)
     &                      +sinj_r*yys(ICVL,icom,1)
             enddo
             if(.not.sw) then
               dhh=0.d0
               do icom=1,ncomp
               dhh=dhh+rcomp(ICOM)*(((((
     &         ahk(5,ICOM) *T
     &        +ahk(4,ICOM))*T
     &        +ahk(3,ICOM))*T
     &        +ahk(2,ICOM))*T
     &        +ahk(1,ICOM))*T
     &        +ahk(0,ICOM))
               enddo
               dhh=sinj_r*dhh
             else
               dhh=sinj_r*enthalpy_ref(ncomp,T,rcomp)
             endif
             dhdt(ICVL)=dhdt(ICVL)+dhh
             diagh(ICVL)=diagh(ICVL)-sinj_r
             diagy(ICVL)=diagy(ICVL)-sinj_r
           endif
 156       continue
         elseif(IMAT.lt.0) then
           IMAT_U=nosld(-IMAT)
         endif
 155    continue
      endif
!------------------------------
! --- User enthalp source term
!------------------------------
       if(src_t.eq.usryes.and.ical_t) then
         do 170 IIMAT=1,NMAT   !ICN=1,NCV
         if(.not.mat_cal(IIMAT)) goto 170
         IMAT=MAT_NO(IIMAT)
         ICVS=MAT_CVEXT(IIMAT-1)+1
         ICVE=MAT_CVEXT(IIMAT)
         if(IMAT.gt.0) then
           IMAT_U=nofld(IMAT)
           DO 175 ICVL=ICVS,ICVE
           ICV=MAT_CV(ICVL)
           if(NPE.gt.1) ICV=NOD_IDg(ICV)
           u=vel(ICVL,1)
           v=vel(ICVL,2)
           w=vel(ICVL,3)
           T=tmp(ICVL,2)
           do icom=1,ncomp
           rcomp(icom)=yys(ICVL,icom,1)
           enddo
           dens=rho(ICVL,1)
           heat1=0.d0
           dhh=0.d0
           grx=CVCENT(1,ICVL)
           gry=CVCENT(2,ICVL)
           grz=CVCENT(3,ICVL)
           vol=CVVOLM(ICVL)
           call user_src_T(deltt,iter,time,
     &          ICV,IMAT_U,iphs,ncomp,grx,gry,grz,
     &          u,v,w,dens,rcomp,vol,T,
     &          heat1)         !heat1=[J/s/m^3]
           dhh=(T-tmp(ICVL,2))*cp(ICVL)*rho(ICVL,1)*rdeltt
           dhdt(ICVL)=dhdt(ICVL)+heat1+dhh
 175       continue
         elseif(IMAT.lt.0) then
           IMAT_U=nosld(-IMAT)
         endif
 170     continue
       endif
!---------------------------------------
! --- Fire wall evaporation source term
!---------------------------------------
      if(src_fire.eq.usryes) then
        do nb=1,nbcnd
        BC_name=boundName(nb,1)
        IIMAT=MAT_BCIDX(nb,1)
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle
        IMAT_U=nofld(IMAT)
        kd=kdbcnd(0,nb)
        no=nobcnd(nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd.eq.kdfire) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVL=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          u=vel(ICVL,1)
          v=vel(ICVL,2)
          w=vel(ICVL,3)
          if(idrdp.eq.mach0.or.idrdp.eq.comp) then
            p=prs(ICVL,1)
          else
            p=prs(ICVL,1)+pp0(ICVL,1)
          endif
          T=tmp(ICVL,1)
          T_WALL=tmp(IDC,1)
          dens=rho(ICVL,1)
          do icom=1,ncomp
            rcomp(icom)=yys(ICVL,icom,1)
            eva_comp(icom)=yys(ICVL,icom,1)
          enddo
          eva_mass=0.d0
          eva_comp(:)=0.d0
          eva_heat=0.d0
          eva_T=tmp(ICVL,1)
!
          call user_src_fire
     &   (no,BC_name,deltt,iter,time,IMAT_U,iphs,ncomp,
     &    u,v,w,p,dens,rcomp,T_WALL,
     &    T,eva_T,eva_comp,eva_mass,eva_heat)
!
                                     ! eva_mass=[kg/(s*m^2)]
          if(eva_mass.gt.SML) then   ! unitmass=[kg/(s/m^3)]
            unitmass=eva_mass*SFAREA(4,ICFL)/CVVOLM(ICVL)
            sumcop=SML
            do icom=1,ncomp
              sumcop=sumcop+eva_comp(icom)
            enddo
            eva_comp(:)=eva_comp(:)/sumcop
            do icom=1,ncomp   !dydt=[kg/(s*m^3)]
              dydt(ICVL,ICOM)=dydt(ICVL,ICOM)+
     &             unitmass*(eva_comp(icom)-0.d0*rcomp(icom)) !cvd
            enddo             
            if(.not.sw) then
              dhh=0.d0
              do icom=1,ncomp
              dhh=dhh+eva_comp(ICOM)*(((((
     &        ahk(5,ICOM) *eva_T
     &       +ahk(4,ICOM))*eva_T
     &       +ahk(3,ICOM))*eva_T
     &       +ahk(2,ICOM))*eva_T
     &       +ahk(1,ICOM))*eva_T
     &       +ahk(0,ICOM))
              enddo
              dhh=dhh*unitmass
            else
              dhh=unitmass*(enthalpy_ref(ncomp,eva_T,eva_comp)
     &           -enthalpy_ref(ncomp,T,eva_comp)) !cvd
            endif
! --- dhdt=[J/s/m^3];cp=J/K/kg; dhh=[J/s/m^3]
            dhdt(ICVL)=dhdt(ICVL)+dhh+eva_heat
            diagh(ICVL)=diagh(ICVL)+0.d0*sinj_r
            diagy(ICVL)=diagy(ICVL)+0.d0*sinj_r
          else
          endif
          enddo
        endif
        enddo
      endif
!------------------------------------------
! --- production term from surface reaction
!------------------------------------------
      if(ical_suf==1) then
        if(ical_FC>0) then
!  	  do nb=1,nbcnd
!          IIMAT=MAT_BCIDX(nb,1)
!          IMAT=MAT_NO(IIMAT)
!          IMAT_U=nofld(abs(IMAT))
!          ISEC=0
!          if(IMAT_U==No_Mem) ISEC=1
!          kd=kdbcnd(0,nb)
!          IBFS=LBC_INDEX(nb-1)+1
!          IBFE=LBC_INDEX(nb)
!          
!          if(surfreac(nb)>0) then
!            if(heat(nb)==0) then
!              dum1=0.d0
!            else
!              dum1=1.d0
!            endif
!            iph=1
!            icoms=phs_idx(iph-1)+1
!            icome=phs_idx(iph)
!            do ipcom=icoms,icome
!              icom=phs_com(ipcom)
!              do IBFL=IBFS,IBFE
!              ICFL=LBC_SSF(IBFL)
!              ICFP=LCYCSF(IBFL)
!              ICV=LVEDGE(1,ICFL)
!              IDC=LVEDGE(2,ICFL)
!              ICVP=LVEDGE(1,ICFP)
!              ICVL=ICVP*ISEC+ICV*(1-ISEC)
!              T=tmp(ICVL,1)
!              vol=CVVOLM(ICVL)
!              ad_comp=SDOT(IBFL,icom)
!! --- ad_comp=mol/m^2/s * m^2 * kg/mol / m^3 =[kg/(s*m^3)]
!              if(.not.sw) then
!                 dhh=ad_comp*(((((
!     &           ahk(5,ICOM) *T
!     &          +ahk(4,ICOM))*T
!     &          +ahk(3,ICOM))*T
!     &          +ahk(2,ICOM))*T
!     &          +ahk(1,ICOM))*T
!     &          +ahk(0,ICOM))
!              else
!                dhh=ad_comp*h_ref(icom,t)*Ri(icom)
!              endif
!              dhdt(ICVL)=dhdt(ICVL)+dhh*dum1
!              dydt(ICVL,ICOM)=dydt(ICVL,ICOM)+ad_comp
!              if(ad_comp>SML) then	    
!              else
!                diagh(ICVL)=diagh(ICVL)-ad_comp*dum1
!                diagy(ICVL)=diagy(ICVL)-ad_comp
!              endif
!              if(ICOM==I_PROTON) then
!                dydt(ICVL,ICOM)=0.d0
!              endif
!              enddo
!            enddo
!          endif
!          enddo
        else
  	  do nb=1,nbcnd
          BC_name=boundName(nb,1)
          IIMAT=MAT_BCIDX(nb,1)
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle
          IMAT_U=nofld(IMAT)
          kd=kdbcnd(0,nb)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(surfreac(nb)>0) then
            if(heat(nb)==0) then
              dum1=0.d0
            else
              dum1=1.d0
            endif
            iph=1
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            do ipcom=icoms,icome
              icom=phs_com(ipcom)
              do 300 IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              T=tmp(ICVL,1)
              vol=CVVOLM(ICVL)
              ad_comp=SDOT(IBFL,icom)*SFAREA(4,ICFL)*wm(icom)/vol
! --- ad_comp=mol/m^2/s * m^2 * kg/mol / m^3 =[kg/(s*m^3)]
              if(.not.sw) then
                 dhh=ad_comp*(((((
     &           ahk(5,ICOM) *T
     &          +ahk(4,ICOM))*T
     &          +ahk(3,ICOM))*T
     &          +ahk(2,ICOM))*T
     &          +ahk(1,ICOM))*T
     &          +ahk(0,ICOM))
              else
                dhh=ad_comp*h_ref(icom,t)*Ri(icom)
              endif
              dhdt(ICVL)=dhdt(ICVL)+dhh*dum1
              dydt(ICVL,ICOM)=dydt(ICVL,ICOM)+ad_comp
              if(ad_comp>SML) then	    
              else
                diagh(ICVL)=diagh(ICVL)-ad_comp*dum1
                diagy(ICVL)=diagy(ICVL)-ad_comp
              endif
 300          enddo
            enddo
          endif
          enddo
        endif
      endif
!-----------------------------------------------
! --- VOF tow phase flow
! --- dydt and dhdt source term, when using VOF 
! --- phase change : dydt
! --- latent heat by phase change : dhdt
!----------------------------------------------
      if(ical_vof==1.and.change_phase.eq.1) then
      endif

!2222
!------------
! --- E2P 
!------------
      if(ieul2ph>0.and.ical_t) then 
      endif
!
! --- 
!
      if(ical_FC==PEFC) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICVL=ICVS,ICVE
        eva_T=tmp(ICVL,1)
        dhh=-(hgsat(IMAT)-hlsat(IMAT))*mdot(ICVL,1)
        IF(ical_t) dhdt(ICVL)=dhdt(ICVL)+dhh
        dydt(ICVL,vap_no)=dydt(ICVL,vap_no)+mdot(ICVL,1) 
                         !dydt=[kg/(s*m^3)]
        enddo
        enddo
      endif
!---------------------------------------------------------
! --- Source term, wdot(:,:)=[kg/(s*m^3)], dydt)=[kg/s] 
!---------------------------------------------------------
      unstdy_h=1.d0
      if(iheat==0) unstdy_h=0.d0
      if(ieul2ph>0) then
      else
        if(ical_vect)  then
          rdsc(:,1)=0.d0
          yys_tmp(:)=0.d0
          do IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
!          ICVE=MAT_INDEX(IIMAT)
          ICVE=MAT_CVEXT(IIMAT)
          do ICOM=1,ncomp
          do ICVL=ICVS,ICVE
          rdsc(ICVL,1)=rdsc(ICVL,1)+
     &      dble(ACT(ICOM))*wdot(ICVL,ICOM)*hform(ICOM)
          dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          dydt(ICVL,ICOM)=dum1*(yys(ICVL,ICOM,2)-yys(ICVL,ICOM,1))
     &           +CVVOLM(ICVL)*(
     &            dydt(ICVL,ICOM)
     &           +wdot(ICVL,ICOM))
          yys_tmp(ICVL)=
     &     min(yys_tmp(ICVL),wdot(ICVL,ICOM)/max(yslw,yys(ICVL,ICOM,2)))
          enddo
          enddo
!
          do ICVL=ICVS,ICVE
          dhdt(ICVL)=dhdt(ICVL)-rdsc(ICVL,1)*unstdy_h
          dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          diagy(ICVL)=dum1+CVVOLM(ICVL)*(diagy(ICVL)-yys_tmp(ICVL))
          enddo
          enddo
        else   !elseif(ical_vect) then
          do 120 IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) goto 120
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          dum2=0.d0
          dhh=0.d0
          do 121 ICOM=1,ncomp
          dhh=dhh+dble(ACT(ICOM))*wdot(ICVL,ICOM)*hform(ICOM)
          dydt(ICVL,ICOM)=dum1*(yys(ICVL,ICOM,2)-yys(ICVL,ICOM,1))
     &           +CVVOLM(ICVL)*(
     &            dydt(ICVL,ICOM)
     &           +wdot(ICVL,ICOM))
          dum2=min(dum2,dble(ACT(ICOM))*wdot(ICVL,ICOM)
     &        /max(yslw,yys(ICVL,ICOM,2)))
  121     enddo
          diagy(ICVL)=dum1+CVVOLM(ICVL)*(diagy(ICVL)-dum2)
          dhdt(ICVL)=dhdt(ICVL)-dhh*unstdy_h
          enddo
  120     enddo
        endif
!
        if(ical_prt>=1) then
          do IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle 
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          do ICOM=1,ncomp
          dydt(ICVL,ICOM)=dydt(ICVL,ICOM)+CVVOLM(ICVL)*EVAP_Y(ICVL,ICOM)
          enddo
          dhh=EVAP_H(ICVL,1)
          dhdt(ICVL)=dhdt(ICVL)+dhh*unstdy_h
          diagh(ICVL)=diagh(ICVL)+EVAP_H(ICVL,2)/cp(ICVL)
!          diagy(ICVL)=diagy(ICVL)+CVVOLM(ICVL)*(-dum2)      !zhangrad
          enddo
          enddo
        endif
      endif
!
      if(i_steady==3.and..false.) then
        do IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
        do ICOM=1,ncomp
        dydt(ICVL,ICOM)=dydt(ICVL,ICOM)-dum1*yys(ICVL,ICOM,2)
        enddo
        dhdt(ICVL)=dhdt(ICVL)-dum1*hhh(ICVL,2)
        enddo
        enddo        
      endif
!
      if(ical_FC>0.and.I_PROTON/=0) then
        dydt(:,I_PROTON)=0.d0
      endif
!
      if(ical_t.and.ical_prp) then
        if(ieul2ph>0) then

        else
          do IIMAT=1,NMAT               !ICV=1,NCV 
          if(.not.mat_cal(IIMAT)) cycle 
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
!
          dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          diagh(ICVL)=dum1+CVVOLM(ICVL)*diagh(ICVL)
!
          dhdt(ICVL)=dum1*(hhh(ICVL,2)-hhh(ICVL,1))       !dhdt=[J/s]
     &            +CVVOLM(ICVL)*dhdt(ICVL)
!
          enddo
          enddo
        endif
      endif
!
!2222
!-------------------------
!--< 1.5 pressure term >--
!-------------------------
! --- / zero Mach no. approx. /
!
      if(idrdp.eq.mach0) then 
        do 150 IIMAT=1,NMAT   !ICV=1,NCV 
        if(.not.mat_cal(IIMAT)) goto 150 
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dhdt(ICVL)=dhdt(ICVL)
     &            +CVVOLM(ICVL)*(pp0(ICVL,1)-pp0(ICVL,2))*rdeltt
        enddo
  150   continue
!----------------------
! --- / compressible / 
!----------------------
      elseif(idrdp.eq.comp) then
        call grad_cell(1,13,
     &       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,prs,grdc(:,:,1))
!        call grad_cell(1,13,
!     &       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
!     &       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,rho,grdc(:,:,2))
!
!      do nb=1,nbcnd
!      IIMAT=MAT_BCIDX(nb,1)
!      if(.not.mat_cal(IIMAT)) goto 325
!      kd=kdbcnd(0,nb)
!      IBFS=LBC_INDEX(nb-1)+1
!      IBFE=LBC_INDEX(nb)
!      if(kd.eq.kdolet) then
!        do IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        ICV=LVEDGE(1,ICFL)
!        IDC=LVEDGE(2,ICFL)
!        grdc(ICV,1,1)=0.d0
!        grdc(ICV,2,1)=0.d0
!        grdc(ICV,3,1)=0.d0
!        enddo
!      endif
!      enddo
!-----------------------------------------
! --- pressure term for enthalpy equation
!-----------------------------------------
!
      if(.true.) then    !=>true
      do 151 IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 151
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_CVEXT(IIMAT)
        ICVE=MAT_INDEX(IIMAT)
        do ICVL=ICVS,ICVE
        dhdt(ICVL)=dhdt(ICVL)+
     &     p_grdc* 
     &     CVVOLM(ICVL)*(
     &    +(prs(ICVL,1)-prs(ICVL,2))*rdeltt
     &    +(vel(ICVL,1)*grdc(ICVL,1,1)
     &     +vel(ICVL,2)*grdc(ICVL,2,1)
     &     +vel(ICVL,3)*grdc(ICVL,3,1))
     &     )
        enddo
 151  enddo
      endif
!
      if(.true.) then !true false okabe  !=>true
        if(.false.) then !=>false
          grdc(:,1,1)=0.d0
          do IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle 
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          grdc(ICVA,1,1)=grdc(ICVA,1,1)-rva(ICFL)
          grdc(ICVB,1,1)=grdc(ICVB,1,1)+rva(ICFL)
          enddo
          enddo
          
          do IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          if(IMAT<0) cycle 
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
          dhdt(ICVL)=dhdt(ICVL)+grdc(ICVL,1,1)/rho(ICVL,1)*prs(ICVL,1)
          enddo
          enddo
        endif

        if(.true.) then  !=>true 
          call grad_cell(3,13,
     &       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vel,grdc)
!
          do IIMAT=1,NMAT 
          if(.not.mat_cal(IIMAT)) cycle 
          IMAT=MAT_NO(IIMAT)
          if(IMAT<0) cycle 
          ICVS=MAT_CVEXT(IIMAT-1)+1 
!          ICVE=MAT_CVEXT(IIMAT)
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
        
          dum1=
     &     p_grdc* 
     &     CVVOLM(ICVL)*( 
!
!     &     +(grdc(ICVL,1,1)
!     &      +grdc(ICVL,2,2)
!     &      +grdc(ICVL,3,3))*prs(ICVL,1)
!new 
     &    +(2.d0*grdc(ICVL,1,1)*grdc(ICVL,1,1)
     &     +2.d0*grdc(ICVL,2,2)*grdc(ICVL,2,2)
     &     +2.d0*grdc(ICVL,3,3)*grdc(ICVL,3,3)
     &     +(grdc(ICVL,2,1)+grdc(ICVL,1,2))**2
     &     +(grdc(ICVL,3,1)+grdc(ICVL,1,3))**2
     &     +(grdc(ICVL,3,2)+grdc(ICVL,2,3))**2)
     &     *(rmut(ICVL)+rmu(ICVL))
!     &     *(rmu(ICVL))
!old
!     &    +(grdc(ICVL,1,1)*grdc(ICVL,1,1)
!     &     +grdc(ICVL,2,2)*grdc(ICVL,2,2)
!     &     +grdc(ICVL,3,3)*grdc(ICVL,3,3)
!     &     +2.d0*grdc(ICVL,2,1)*grdc(ICVL,2,1)
!     &     +2.d0*grdc(ICVL,3,1)*grdc(ICVL,3,1)
!     &     +2.d0*grdc(ICVL,3,2)*grdc(ICVL,3,2))
!!     &     *(rmut(ICVL)+rmu(ICVL))
!     &     *(rmu(ICVL))

     &      )
            dhdt(ICVL)=dhdt(ICVL)+dum1
          enddo
          enddo
          endif

        endif
      endif
!
!
      if(ical_sld/=0.and..false.) then
!      if(ical_sld/=0) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.le.0) cycle
        if(ishaft(IMAT)==1) then 
          ICVS=MAT_CVEXT(IIMAT-1)+1 
          ICVE=MAT_CVEXT(IIMAT)
          shaftn(1)=end(1,IMAT)-begin(1,IMAT) 
          shaftn(2)=end(2,IMAT)-begin(2,IMAT) 
          shaftn(3)=end(3,IMAT)-begin(3,IMAT) 
          dum1=dsqrt(shaftn(1)**2+shaftn(2)**2+shaftn(3)**2) 
          shaftn(1)=shaftn(1)/dum1
          shaftn(2)=shaftn(2)/dum1
          shaftn(3)=shaftn(3)/dum1
          DO ICVL=ICVS,ICVE 
          radius(1)=CVCENT(1,ICVL)-begin(1,IMAT)   !IMAT??? 
          radius(2)=CVCENT(2,ICVL)-begin(2,IMAT)
          radius(3)=CVCENT(3,ICVL)-begin(3,IMAT)
!
          dum3=radius(1)*shaftn(1)
     &        +radius(2)*shaftn(2)
     &        +radius(3)*shaftn(3)
          radius(1)=radius(1)-dum3*shaftn(1)
          radius(2)=radius(2)-dum3*shaftn(2)
          radius(3)=radius(3)-dum3*shaftn(3)
          vr(:)=vel(ICVL,:)
          call cal_AXBEQC(shaftn,radius,velref,dum2)
!          
          dum1=(rotati(IMAT)*dum2)**2*rho(ICVL,1)
          dum2=(vr(1)**2+vr(2)**2+vr(3)**2)*rho(ICVL,1)
          dhdt(ICVL)=dhdt(ICVL)
     &              +0.5d0*(dum1-dum2)*rdeltt*CVVOLM(ICVL)
          ENDDO
        endif
        enddo
      endif
!--------------------------
!--< 1.4 convection term >-
!--------------------------,dydt(100)
      if(ical_t.and.ical_prp) then
        call conv_term   !NNOTVECT
     & (icnvty,lmtrty,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  rvd(:,1:2),vel,
     &  grdc(:,:,1),rva,hhh,dhdt,
     &  -1)
      endif
!
      do 140 ICOM=1,ncomp
      IMODE=ICOM
      yys_tmp(:)=yys(:,ICOM,1)
      call conv_term
     & (icnvty,lmtrty,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  rvd(:,1:2),vel,
     &  grdc(:,:,1),rva,yys_tmp,dydt(:,ICOM),
     &  IMODE)
  140 continue
!2222
!-------------------------------
! --- 
!-------------------------------
      if(ical_porous==1) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.le.0) cycle
        IMAT_U=nofld(IMAT)
        if(IMAT_U>N_pors) then
          if(porous(IMAT)==idarcy2) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            dum1=porosty(IMAT)
            do ICOM=1,ncomp
            DO ICVL=ICVS,ICVE
            dhdt(ICVL)=dhdt(ICVL)*dum1
            dydt(ICVL,ICOM)=dum1*dydt(ICVL,ICOM)
            enddo
            enddo
          endif
        endif
        enddo
      endif
!
!-------------------------------
! --- 
!-------------------------------
      if(ical_FC==PEFC) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.le.0) cycle
        IMAT_U=nofld(IMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICVL=ICVS,ICVE
        dum1=1.d0-aks(ICVL,ical_s)
        dhdt(ICVL)=dhdt(ICVL)*dum1
        do ICOM=1,ncomp
        dydt(ICVL,ICOM)=dum1*dydt(ICVL,ICOM)
        enddo
        enddo
        enddo
      endif
!
!--------------------------
!--< 1.3 diffusion term >--
!--------------------------
!
!      viscous=.false.
!      noviscs=.false.
!      do 162 IIMAT=1,NMAT
!      IMAT=MAT_NO(IIMAT)
!      if(IMAT.lt.0) goto 162
!      if(ivscty(IMAT).eq.cal) then
!        viscous=.true.
!      else
!        noviscs=.true.
!      endif
! 162  continue
!      if(viscous) then
!
      call hys_diff(iphs,   !NNOTVECT
     &  LVEDGE,LCYCSF,LBC_SSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  LCYCOLD,wifsld,OPPANG,rho(:,1),
     &  rmd,rds,rmu,rmut,tmp,yys,hhs,
     &  utau,frstcv,vel,
     &  tmpbnd,yysbnd,htcbnd,mtcbnd,radbnd,HTFLUX,
     &  grdc(:,:,1),rvd(:,1),rvd(:,2),cp,vctr,aks,
     &  dhdt,dydt,rmdc,rdsc,rvds,dsclt,dscly
     & )
      do 163 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        if(IMAT>0) then
          IMAT_U=nofld(abs(IMAT))
          if(ivscty(IMAT)/=cal) then
            rmdc(ICVS:ICVE)=0.d0
            rdsc(ICVS:ICVE,:)=0.d0
            rvds(ICFS:ICFE,:)=0.d0
            dsclt(ICVS:ICVE)=0.d0
            dscly(ICVS:ICVE)=0.d0

            rmdc(IDCS:IDCE)=0.d0
            rdsc(IDCS:IDCE,:)=0.d0
            rvds(IDCS:IDCE,:)=0.d0
            dsclt(IDCS:IDCE)=0.d0
            dscly(IDCS:IDCE)=0.d0
          endif
        else
          IMAT_U=nosld(abs(IMAT))
        endif
 163  enddo
!----------------------------
!--< 1.6 Clear solid part >--
!----------------------------
      do 160 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 160
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICOM=1,NCOMP
        dydt(ICVS:ICVE,ICOM)=0.d0
        enddo
      endif
 160  continue
!
      if(ical_FC>0) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT>0) then
          IMAT_U=nofld(IMAT)
          if(IMAT_U==No_Mem) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICOM=1,NCOMP
            if(ICOM/=I_H2O) then
              dydt(ICVS:ICVE,ICOM)=0.d0
            endif
            enddo
          endif
        endif
        enddo
      endif
!--------------------------
! --- 
!--------------------------yys(2,:,1)
      if(ical_t.AND.radflag==1) then
        IF(iter-iters.lt.20) then
	  do IIMAT=1,NMAT    !ICV=1,NCV
	  IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
	  do ICVL=ICVS,ICVE
	  dhdtmax=1.d0
	  if(abs(RadHeat(ICVL)).gt.1.d-2) then
	    dhdtmax=abs(dhdt(ICVL)/RadHeat(ICVL))
	    if(dhdtmax.gt.0.1) then
              dhdtmax=1.d0
            endif
          end if
          dhdt(ICVL)=dhdt(ICVL)-RadHeat(ICVL)*dhdtmax
          enddo
          enddo
	ELSE
	  do IIMAT=1,NMAT    !ICV=1,NCV
            IMAT=MAT_NO(IIMAT)
	    ICVS=MAT_CVEXT(IIMAT-1)+1
	    ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
            dhdt(ICVL)=dhdt(ICVL)-RadHeat(ICVL) 
	    enddo
	  end do
	ENDIF
      endif
!-----------------------------------------------
!-< 2. Update enthalpy & mass fracion >-
! --- 
!-----------------------------------------------
!      if(u_func(3)==1) then
!        do IIMAT=1,NMAT   !ICN=1,NCV
!        if(.not.mat_cal(IIMAT)) cycle
!        IMAT=MAT_NO(IIMAT)
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_CVEXT(IIMAT)
!        if(IMAT.gt.0) then
!          IMAT_U=nofld(IMAT)
!          DO ICVL=ICVS,ICVE
!            idum=multiR(ICVL)
!            if(idum>=21.and.idum<=50) then
!              dhdt(ICVL)=0.d0
!              dydt(ICVL,:)=0.d0
!            endif
!          enddo
!        elseif(IMAT.lt.0) then
!          IMAT_U=nosld(-IMAT)
!        endif
!        enddo        
!      endif
!--------------------
!--< 2.1 explicit >--
!--------------------
      if(intgty.eq.eulere ) then
        do IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICOM=1,ncomp
        do ICVL=ICVS,ICVE
        dum1=1.d0/diagy(ICVL)
        dydt(ICVL,ICOM)=dum1*dydt(ICVL,ICOM)
        enddo
        enddo
        enddo
!
        if(ical_t.and.ical_prp) then
          do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
            dhdt(ICVL)=dhdt(ICVL)/diagh(ICVL)
          enddo
          enddo
        endif
!--------------------
!--< 2.2 implicit >--
!--------------------
      else
        cofd=1.d0
!--------------------------------
! / set boundary condition flag /
!--------------------------------
        call bc_kdbty(LBC_SSF,LCYCSF,mat_cal,kdbt,kdby)
!---------------
! / enthalpy /
!---------------
        if(ical_t.and.ical_prp) then 
          do 215 IIMAT=1,NMAT     !ICV=1,NALLCV 
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT)) goto 215
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
            rmdc(ICVL)=rmdc(ICVL)/cp(ICVL)
          enddo
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          do ICVL=IDCS,IDCE
            rmdc(ICVL)=rmdc(ICVL)/cp(ICVL)
          enddo
 215      continue
!
          do 220 IIMAT=1,NMAT     !ICF=1,NCVFAC
          if(.not.mat_cal(IIMAT)) goto 220
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVB=LVEDGE(2,ICFL)
          ICVA=LVEDGE(1,ICFL)
          dsclt(ICVB)=dsclt(ICVB)/cp(ICVA)
          enddo
  220     continue
!
          rvd=0.d0   ! temperature unboundness, if rvd/=0.d0
	  ndf=1      ! it is NOT necessary in delt-h eq. 
	  nq=1
          IMODE=1
!
          if(ical_vect) then
            ndiag=1   !integer :: ndiag
            call solve_cnvdif_vect
     &      (.false.,mxssfbc,ndf,nq,'h',time,ndiag,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdbt,rva,rvd(:,1:2),rmdc,diagh,
     &  LCYCOLD,wifsld,OPPANG,
     &  cofd,dsclt,htcbnd,dhdt,deltt,vctr,
     &  itery(0:),repsy(0:),aepsy(0:),iter,ierror,IMODE)
          else
            call solve_cnvdif
     &      (.false.,mxssfbc,ndf,nq,'h',time,relaxh,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,SFCENT,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdbt,rva,rvd(:,1:2),rmdc,diagh,
     &  LCYCOLD,wifsld,OPPANG,CVVOLM,
     &  cofd,dsclt,htcbnd,dhdt,deltt,
     &  itery(0:),repsy(0:),aepsy(0:),iter,ierror,IMODE)
          if(ierror.ne.0) goto 9999
          endif
        endif
!
!------------------------
! --- / mass fraction /
!------------------------
        ndf=ncomp
        nq=ncomp
        rvd=0.d0
        duma=0.d0
        IMODE=2
        if(ical_vect) then
           ndiag=1
           call solve_cnvdif_vect
     &    (.true.,1,ndf,nq,'y',time,ndiag,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdby,rva,rvd(:,1:2),rdsc,diagy,
     &  LCYCOLD,wifsld,OPPANG,
     &  cofd,dscly,duma,dydt,deltt,vctr,
     &  itery(1:),repsy(1:),aepsy(1:),iter,ierror,IMODE)
        else
          call solve_cnvdif
     &    (.true.,1,ndf,nq,'y',time,relaxys,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,SFCENT,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdby,rva,rvd(:,1:2),rdsc,diagy,
     &  LCYCOLD,wifsld,OPPANG,CVVOLM,
     &  cofd,dscly,duma,dydt,deltt,
     &  itery(1:),repsy(1:),aepsy(1:),iter,ierror,IMODE)
        endif
!
        if(ierror/=0) call FFRABORT(1,'ERR: at solve_cnvdif =y')
!
!        do icom=1,ncomp
!        if(itery(icom)==0) then
!          do IIMAT=1,NMAT    !ICV=1,NCV
!          if(.not.mat_cal(IIMAT)) cycle
!          IMAT=MAT_NO(IIMAT)
!          ICVS=MAT_CVEXT(IIMAT-1)+1
!          ICVE=MAT_CVEXT(IIMAT)
!          do ICVL=ICVS,ICVE
!          dydt(ICVL,ICOM)=dydt(ICVL,ICOM)/diagy(ICVL)
!          enddo
!          enddo
!        endif
!        enddo
!
      endif
!
!-------------
! --- E2P 
!-------------
!
      if(ieul2ph>0) then
      endif
!-------------------------------------------------------------------
!--< 2.3 reset temp. & mass frac. at (n+1) time level >-- hhh[J/kg]
!-------------------------------------------------------------------
      if(.false.) then   !OK  .fasle.
        do IIMAT=1,NMAT    !ICV=1,NCV
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        hhh(ICVL,1)=hhh(ICVL,1)+dhdt(ICVL)*relaxh(IIMAT)    !hhh=[J/kg]
!
        dum2=0.d0
        do ICOM=1,ncomp
        dum1=max(0.d0,yys(ICVL,ICOM,1)+relaxys(IIMAT)*dydt(ICVL,ICOM))
        dum2=dum2+dum1*dble(ACT(ICOM))
        ys(ICOM)=dum1
        enddo
        dum2=1.d0/dum2
        do ICOM=1,ncomp
        if(ACT(ICOM)==0) then
          yys(ICVL,ICOM,1)=ys(ICOM)
        else
          yys(ICVL,ICOM,1)=ys(ICOM)*dum2
        endif
        enddo
        enddo
        enddo
      endif
!----------------------
! --- BackGrand gas
!----------------------
!
      if(ical_vect) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT>0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          hhh(ICVL,1)=hhh(ICVL,1)+dhdt(ICVL)*relaxh(IIMAT)
          rdsc(ICVL,BGCOMP(IMAT))=max(0.d0,yys(ICVL,BGCOMP(IMAT),1)
     &      +relaxys(IIMAT)*dydt(ICVL,BGCOMP(IMAT)))
          enddo
          yys_tmp(:)=0.d0
          do ICOM=1,ncomp
          if(ICOM==BGCOMP(IMAT)) cycle
          do ICVL=ICVS,ICVE
          dum1=max(0.d0,yys(ICVL,ICOM,1)+relaxys(IIMAT)*dydt(ICVL,ICOM))
          rdsc(ICVL,ICOM)=dum1
          yys_tmp(ICVL)=yys_tmp(ICVL)+dum1
          enddo
          enddo
          do ICVL=ICVS,ICVE
          if(yys_tmp(ICVL)-1.d0>SML)then
            yys_tmp(ICVL)=1.d0/(yys_tmp(ICVL)+rdsc(ICVL,BGCOMP(IMAT)))
          else
            rdsc(ICVL,BGCOMP(IMAT))=max(0.d0,(1.d0-yys_tmp(ICVL)))
            yys_tmp(ICVL)=1.d0
          endif
          enddo
          do ICOM=1,ncomp
          do ICVL=ICVS,ICVE
          yys(ICVL,ICOM,1)=yys_tmp(ICVL)*rdsc(ICVL,ICOM)
          enddo
          enddo
        else
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          yys(ICVS:ICVE,:,1)=0.d0
         hhh(ICVS:ICVE,1)=hhh(ICVS:ICVE,1)+dhdt(ICVS:ICVE)*relaxh(IIMAT) 
        endif
        enddo
      else
        if(.false.) then    !.false.
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT)) cycle
          if(IMAT>0) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
            hhh(ICVL,1)=hhh(ICVL,1)+dhdt(ICVL)*relaxh(IIMAT)
            dum3=max(0.d0,yys(ICVL,BGCOMP(IMAT),1)
     &        +relaxys(IIMAT)*dydt(ICVL,BGCOMP(IMAT)))
            ys(BGCOMP(IMAT))=dum3
            dum2=0.d0
            dum3=0.d0
            do 231 ICOM=1,ncomp
            if(ICOM==BGCOMP(IMAT)) cycle
            dum1=max(0.d0,yys(ICVL,ICOM,1)
     &          +relaxys(IIMAT)*dydt(ICVL,ICOM))
            dum2=dum2+dum1*dble(ACT(ICOM))
            dum3=dum3+dum1
            ys(ICOM)=dum1
  231       continue
            if((dum2-1.d0)>SML) then
              dum2=1.d0/(dum2+ys(BGCOMP(IMAT)))
!              dum3=1.d0/(dum3+ys(BGCOMP(IMAT)))
              do 232 ICOM=1,ncomp 
              if(ACT(ICOM)==0) then 
!                yys(ICVL,ICOM,1)=ys(ICOM)*dum3 
              else 
                yys(ICVL,ICOM,1)=ys(ICOM)*dum2 
              endif 
  232         enddo 
              ys(BGCOMP(IMAT))=0.d0  
            else
              yys(ICVL,:,1)=ys(:) 
              yys(ICVL,BGCOMP(IMAT),1)=max(0.d0,(1.d0-dum2))
            endif
            dum3=1.d0/(dum3+ys(BGCOMP(IMAT)))
            do ICOM=1,ncomp
            if(ACT(ICOM)==1) cycle
            yys(ICVL,ICOM,1)=ys(ICOM)*dum3
            enddo
            enddo
          else
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            yys(ICVS:ICVE,:,1)=0.d0
            hhh(ICVS:ICVE,1)=
     &      hhh(ICVS:ICVE,1)+dhdt(ICVS:ICVE)*relaxh(IIMAT) 
          endif
          enddo
        elseif(.true.) then    !  .true.
          do IIMAT=1,NMAT  
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT)) cycle 
          if(IMAT>0) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
            hhh(ICVL,1)=hhh(ICVL,1)+dhdt(ICVL)*relaxh(IIMAT)    !hhh=[J/kg]
            dum2=0.d0
            dum3=0.d0
            do ICOM=1,ncomp
            dum1=max(0.d0,yys(ICVL,ICOM,1)+relaxys(IIMAT)*dydt(ICVL,ICOM))
            dum2=dum2+dum1*dble(ACT(ICOM))
            dum3=dum3+dum1
            ys(ICOM)=dum1
            enddo
            dum2=1.d0/dum2
            dum3=1.d0/dum3
            do ICOM=1,ncomp
            if(ACT(ICOM)==0) then
              yys(ICVL,ICOM,1)=ys(ICOM)*dum3
            else
              yys(ICVL,ICOM,1)=ys(ICOM)*dum2
            endif
            enddo
            enddo
          else   !bate
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            yys(ICVS:ICVE,:,1)=0.d0
            hhh(ICVS:ICVE,1)=
     &      hhh(ICVS:ICVE,1)+dhdt(ICVS:ICVE)*relaxh(IIMAT) 
          endif
          enddo
        endif
      endif
! check mass 
      if(.false.) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        dum1=0.d0
        do ICOM=1,ncomp
        if(ICOM/=1) cycle
        do ICVL=ICVS,ICVE
        dum1=dum1+yys(ICVL,ICOM,1)*rho(ICVL,1)*CVVOLM(ICVL)
        enddo
        enddo
        enddo
        print*,'mass in gas',dum1,(dum1-dumm)/dumm
        dumm=dum1
      endif
!
!---------------------------
! --- Convergence error: 
!---------------------------
      if(ical_vect) then
      erry(:)=ZERO
      erryys=ZERO
      Do IIMAT=1,NMAT    !ICV=1,NCVIN
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)

      do ICVL=ICVS,ICVE
      erry(0)=erry(0)+dhdt(ICVL)*dhdt(ICVL)
      erryys(0)=erryys(0)+hhh(ICVL,1)*hhh(ICVL,1)
      enddo

      do ICOM=1,ncomp
      do ICVL=ICVS,ICVE
      erry(ICOM)=erry(ICOM)+dydt(ICVL,ICOM)*dydt(ICVL,ICOM)
      erryys(ICOM)=erryys(ICOM)+yys(ICVL,ICOM,1)*yys(ICVL,ICOM,1)
      enddo
      enddo
!
      do ICOM=0,ncomp
      erry(ICOM)=dsqrt(erry(ICOM))/(dsqrt(erryys(ICOM))+1.d-8)
      enddo
      enddo

      else
      erry(:)=ZERO
      erryys=ZERO
      Do 320 IIMAT=1,NMAT    !ICV=1,NCVIN
      if(.not.mat_cal(IIMAT)) goto 320
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      erry(0)=erry(0)+dhdt(ICVL)*dhdt(ICVL)
      erryys(0)=erryys(0)+hhh(ICVL,1)*hhh(ICVL,1)
      do 325 ICOM=1,ncomp
      erry(ICOM)=erry(ICOM)+dydt(ICVL,ICOM)*dydt(ICVL,ICOM)
      erryys(ICOM)=erryys(ICOM)+yys(ICVL,ICOM,1)*yys(ICVL,ICOM,1)
 325  continue
      enddo
  320 continue
!
      do 328 ICOM=0,ncomp
      erry(ICOM)=dsqrt(erry(ICOM))/(dsqrt(erryys(ICOM))+1.d-8)
 328  continue
      endif
!
!------------------------ -----------------------------
! --- calculate temperature from enthalpy & specific
!------------------------ -----------------------------
!
!      if(ical_t.and..NOT.inoheat) then
      if(ical_t.and.ical_prp) then
        if(ical_vof==1) then
        elseif(ical_cavi==1) then
        else
          call cal_h2t(iphs,NCV,
     &    MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &    hhh(:,1),yys(:,:,1),tmp(:,1),prs(:,1),rho(:,1))
        endif
      endif
!
!
      DEALLOCATE (rdsc,dhdt,dydt,rvds)  
!
      return
!
 9999 continue
!
      if(my_rank.eq.ROOT) write(ifle,*) '(hys_admin)' 
      ierror=1
!
      end subroutine hys_admin
!
      
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$LCYCSF$$$$$$$$$$$$
      subroutine corre_grdc(LBC_SSF,LVEDGE,SFAREA,mat_cal,LCYCSF,grdc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_io,only      : ifle,ifll
      use module_boundary,only: nbcnd,kdbcnd,boundName,
     &                          MAT_BCIDX,LBC_INDEX,kdolet,kdilet,
     &                          kxnone,kdfire,kdintr,kdsymm,kdcvd,
     &                          nobcnd,surfreac,openout,kdprdc,
     &                          idis,kdsld,kdtchi
      use module_model,only   : ical_dens,idrdp,comp
!
      implicit none
!    
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      real*8 ,intent(inout) :: grdc      (   MXALLCV,3,3)
!
! --- [local entities]
!
      integer :: NB,IIMAT,IBFS,IBFE,IBFL,kd,ICFL,ICV,IDC,ICFP,IDCP,ICVP
      real*8  :: grx,gry,grz,dum1
      
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      do 324 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) goto 324
        kd=kdbcnd(0,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd.eq.kdolet.and.openout(nb)==10)then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &         grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &         grdc(ICV,3,1)*SFAREA(3,ICFL)
          grx=grdc(ICV,1,1)-dum1*SFAREA(1,ICFL)
          gry=grdc(ICV,2,1)-dum1*SFAREA(2,ICFL)
          grz=grdc(ICV,3,1)-dum1*SFAREA(3,ICFL)
          grdc(IDC,1,1)=grx
          grdc(IDC,2,1)=gry
          grdc(IDC,3,1)=grz

!          grdc(ICV,1,1)=grx
!          grdc(ICV,2,1)=gry 
!          grdc(ICV,3,1)=grz
          enddo
        ELSEif(kd.eq.kdolet)then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &         grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &         grdc(ICV,3,1)*SFAREA(3,ICFL)
          grx=grdc(ICV,1,1)-dum1*SFAREA(1,ICFL)
          gry=grdc(ICV,2,1)-dum1*SFAREA(2,ICFL)
          grz=grdc(ICV,3,1)-dum1*SFAREA(3,ICFL)
          grdc(IDC,1,1)=grx
          grdc(IDC,2,1)=gry
          grdc(IDC,3,1)=grz
          enddo
        ELSEif(ical_dens==4.and.(kd.eq.kdilet.or.kd.eq.kdtchi)) then

        ELSEif(idrdp.eq.comp.and.(kd.eq.kdilet.or.kd.eq.kdtchi)) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &         grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &         grdc(ICV,3,1)*SFAREA(3,ICFL)
          grx=grdc(ICV,1,1)-dum1*SFAREA(1,ICFL)
          gry=grdc(ICV,2,1)-dum1*SFAREA(2,ICFL)
          grz=grdc(ICV,3,1)-dum1*SFAREA(3,ICFL)
          grdc(ICV,1,1)=grx!*0.d0
          grdc(ICV,2,1)=gry!*0.d0
          grdc(ICV,3,1)=grz!*0.d0
          grdc(IDC,1,1)=grx!*0.d0
          grdc(IDC,2,1)=gry!*0.d0
          grdc(IDC,3,1)=grz!*0.d0
          enddo
        ELSEif(kd.eq.kdilet.or.kd.eq.kdtchi) then
          cycle
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &         grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &         grdc(ICV,3,1)*SFAREA(3,ICFL)
          grdc(ICV,1,1)=dum1*SFAREA(1,ICFL)
          grdc(ICV,2,1)=dum1*SFAREA(2,ICFL)
          grdc(ICV,3,1)=dum1*SFAREA(3,ICFL)
          grdc(IDC,1,1)=dum1*SFAREA(1,ICFL)
          grdc(IDC,2,1)=dum1*SFAREA(2,ICFL)
          grdc(IDC,3,1)=dum1*SFAREA(3,ICFL)
          enddo

        elseif(kd.eq.kxnone.or.kd.eq.kdfire.or.
     &         kd.eq.kdintr.or.kd.eq.kdcvd) then
! --- wall (???zhang)
!           CYCLE
!          if(iLBF_P==2.or.iLBF_P==1) cycle  !5555
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &         grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &         grdc(ICV,3,1)*SFAREA(3,ICFL)
          grdc(IDC,1,1)=dum1*SFAREA(1,ICFL)
          grdc(IDC,2,1)=dum1*SFAREA(2,ICFL)
          grdc(IDC,3,1)=dum1*SFAREA(3,ICFL)
          enddo 
        elseif(kd==kdprdc.and.idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &         grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &         grdc(ICV,3,1)*SFAREA(3,ICFL)
          grdc(IDC,1,1)=dum1*SFAREA(1,ICFL)
          grdc(IDC,2,1)=dum1*SFAREA(2,ICFL)
          grdc(IDC,3,1)=dum1*SFAREA(3,ICFL)
          enddo
        elseif(kd==kdsld.and.idis(nb)>=1) then !
        endif
 324  continue
!
      
      return
      end subroutine corre_grdc
