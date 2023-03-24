!	subroutine hys_diff()
!       subroutine vel_diff()
!       subroutine vel_diff_e2p
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hys_diff(iphs,
     &  LVEDGE,LCYCSF,LBC_SSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  LCYCOLD,wifsld,OPPANG,rho,
     &  rmd,rds,rmu,rmut,tmp,yys,hhs,
     &  utau,frstcv,vel,
     &  tmpbnd,yysbnd,htcbnd,mtcbnd,radbnd,HTFLUX,
     &  grdc,dflxt,dflxh,cp,vctr,aks,
     &  dhdt,dydt,rmdc,rdsc,dflxy,dsclt,dscly
     &  )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     dflxt(:) <= rvd(1,:)
!     dflxh(:) <= rvd(1,2)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_Euler2ph,only : ieul2ph
      use module_turbparm,only : prturb,scturb
!
      use module_metrix,only   : dyc  =>d2work4       !dyc  =>W2K1
      use module_metrix,only   : dhc  =>W1K8          !W1K8
      use module_metrix,only   : tmpfac=>d2vect
      use module_metrix,only   : iwork1

      use module_model,only    : ical_vect,PEFC,ical_thmoDIFF
!
      use module_FUEL  ,only   : vap_no,No_Mem
      use module_material,only : nflud,nofld,ivscty,cal,BGCOMP
      use module_scalar  ,ONLY : ical_FC
      use module_boundary,only : kdintr,kdbcnd,MAT_BCIDX,LBC_INDEX,nbcnd
      use module_species,only  : act
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V
      use module_model,only    : expfact
      use module_metrix,only   : TH_DIFF
!
! 1.  Calculate thermal conductivity & diffusion term 
!
      implicit none
!
! --- [dummy arguments]
!
      INTEGER,INTENT(IN)    :: iphs
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)    
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(inout) :: rmd       (   MXALLCV)
      real*8 ,intent(in) :: rds       (   MXALLCV,mxcomp)
      real*8 ,intent(inout) :: rmut      (   MXALLCV)
      real*8 ,intent(in) :: rmu      (   MXALLCV)
      real*8 ,intent(in)    :: tmp       (   MXALLCV)
      real*8 ,intent(inout)    :: yys       (   MXALLCV,mxcomp)
      real*8 ,intent(in)    :: hhs       (   MXALLCV,mxcomp)
      real*8 ,intent(in)    :: tmpbnd    (   mxssfbc)
      real*8 ,intent(in)    :: yysbnd    (   mxssfbc,mxcomp)
      real*8 ,intent(inout) :: htcbnd    (   mxssfbc)
      real*8 ,intent(inout) :: HTFLUX    (   mxssfbc)
      real*8 ,intent(in)    :: mtcbnd    (   mxssfbc)
      real*8 ,intent(in)    :: radbnd    (   mxssfbc)
      real*8 ,intent(inout) :: dhdt      (   MXALLCV)
      real*8 ,intent(inout) :: dydt      (   MXALLCV,mxcomp)
      real*8 ,intent(inout) :: grdc      (   MXALLCV,3)
      real*8 ,intent(inout) :: dflxt     (   MXCVFAC)
      real*8 ,intent(inout) :: dflxh     (   MXCVFAC)
      real*8 ,intent(out)   :: rmdc      (   MXALLCV)          !=>solve
      real*8 ,intent(out)   :: rdsc      (   MXALLCV,mxcomp)   !=>solve
      real*8 ,intent(out)   :: dflxy     (   MXCVFAC,mxcomp)
      real*8 ,intent(out)   :: dsclt     (   MXALLCV)
      real*8 ,intent(out)   :: dscly     (   MXALLCV)
      real*8 ,intent(in)    :: cp        (   MXALLCV)
      integer,intent(in)    :: LCYCOLD   (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG    (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(in)    :: rho       (   MXALLCV)
      real*8 ,intent(in)    :: aks       (   MXALLCVR,mxrans)

      REAL*8 ,INTENT(IN)   :: FRSTCV(        MXSSFBC)
      REAL*8 ,INTENT(IN)   :: UTAU  (        0:MXSSFBC)
      real*8 ,intent(in)   :: vel       (   MXALLCV,3  )
!
! --- [local entities]
!
      real*8  :: dd1,dd2,dd0,dum1,dll
      real*8,parameter :: eps=tiny(1.d0)*1.d2,SML=1.D-25,dirr=0.d0
      integer :: i,j,k,l,m,n,ierr=0
      integer :: IMAT,IMAT2,IIMAT,IIMAT2,IMAT_U,IMAT_U2,
     &   ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE,ICFP
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD,ICM
      integer :: nb,kd,kdt,kdy,kdk,kdp,IE
      real*8  :: wi1,wi2,grx,gry,grz,dx,dy,dz,dl,dlvect,Fs,TC,TP
      real*8  :: grdf(3),dxx,dyx,dzx,gf0,gf1,gf2,gf3,dum2,dum3,dum4,dum5,dum6
!
!-< 2. Calculate diffusion flux >-
!
      allocate(dyc(MXALLCV,mxcomp),stat=ierr)
      if(ierr.ne.0)  call FFRABORT(1,'hys_diff/allocation [dyc]')
!
!
      
!--< 2.1 inner region >--
!----------------------------------------------
! --- prturb : Turbulence Prandtl number 
! --- scturb : Turbulence Schmidt number rmd
! --- enthalpy diffusion term: 
!----------------------------------------------
      dflxy=0.d0
!
      dum1=1.d0/prturb
      do 200 IIMAT=1,NMAT     !ICV=1,NALLCV
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 200
      if(IMAT.gt.0) then        
        if(ivscty(IMAT).eq.cal) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          rmdc(ICVL)=rmd(ICVL)+rmut(ICVL)*dum1*cp(ICVL)
          enddo
!
          ICVS=MAT_DCIDX(IIMAT-1)+1
          ICVE=MAT_DCIDX(IIMAT)
!
          do ICVL=ICVS,ICVE
          rmdc(ICVL)=rmd(ICVL)+rmut(ICVL)*dum1*cp(ICVL)
          enddo
        else
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          rmdc(ICVS:ICVE)=0.d0
          ICVS=MAT_DCIDX(IIMAT-1)+1
          ICVE=MAT_DCIDX(IIMAT)
          rmdc(ICVS:ICVE)=0.d0
        endif
      elseif(IMAT.lt.0) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        rmdc(ICVL)=rmd(ICVL)
        enddo
        ICVS=MAT_DCIDX(IIMAT-1)+1
        ICVE=MAT_DCIDX(IIMAT)
        do ICVL=ICVS,ICVE
        rmdc(ICVL)=rmd(ICVL)
        enddo
      endif
  200 enddo
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmdc)
!
      dum1=1.d0/scturb
!
      do 210 IIMAT=1,NMAT     !   ICV=1,NALLCV
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 210

      if(IMAT>0) then
      if(ivscty(IMAT).eq.cal) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 211 ICOM=1,ncomp
        do ICVL=ICVS,ICVE
        rdsc(ICVL,ICOM)=rds(ICVL,ICOM)+rmut(ICVL)*dum1
        enddo
  211   enddo
        ICVS=MAT_DCIDX(IIMAT-1)+1
        ICVE=MAT_DCIDX(IIMAT)
        do 212 ICOM=1,ncomp
        do ICVL=ICVS,ICVE
        rdsc(ICVL,ICOM)=rds(ICVL,ICOM)+rmut(ICVL)*dum1
        enddo
 212    enddo
      else
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        rdsc(ICVS:ICVE,:)=0.d0
        ICVS=MAT_DCIDX(IIMAT-1)+1
        ICVE=MAT_DCIDX(IIMAT)
        rdsc(ICVS:ICVE,:)=0.d0
      endif
      else
        ICVS=MAT_DCIDX(IIMAT-1)+1
        ICVE=MAT_DCIDX(IIMAT)
        do ICOM=1,ncomp
        do ICVL=ICVS,ICVE
        rdsc(ICVL,ICOM)=0.d0
        enddo
        enddo
      endif
 210  enddo
!
      
!

      do ICOM=1,ncomp
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &     mat_cal,rdsc(:,ICOM))
      enddo
!
!      call dc_symprs(ncomp,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
!     &     mat_cal,rdsc)
!--------------------------------------
! --- Soret effect
!--------------------------------------
      dyc(1:MXALLCV,1:mxcomp)=0.d0
      if(ical_thmoDIFF==1) then
        dsclt(:)=tmp(:)
        call grad_cell(1,10,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,dsclt,grdc)
        call dc_symprv
     &  (3,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  LCYCOLD,wifsld,OPPANG,
     &  SFAREA,SFCENT,grdc,1,0)

!
        do 500 ICOM=1,ncomp
        do IIMAT=1,NMAT 
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        if(IMAT<0) then
          cycle
        endif

        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        dd1=TH_DIFF(ICVA,ICOM)
        dd2=TH_DIFF(ICVB,ICOM)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dum2=wi1*tmp(ICVA)+wi2*tmp(ICVB)
        dd0=dd1*dd2/(dd1*wi1+dd2*wi2+SML)
!        dd0=wi1*dd1+wi2*dd2
        dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)  
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl
        dlvect=abs(dx*SFAREA(1,ICFL)
     &            +dy*SFAREA(2,ICFL)
     &            +dz*SFAREA(3,ICFL))
        grx=wi1*grdc(ICVA,1)+wi2*grdc(ICVB,1)
        gry=wi1*grdc(ICVA,2)+wi2*grdc(ICVB,2)
        grz=wi1*grdc(ICVA,3)+wi2*grdc(ICVB,3)
        gf1=(dsclt(ICVB)-dsclt(ICVA))/dl/dum2
        gf2=grx*(SFAREA(1,ICFL)-dx)
     &     +gry*(SFAREA(2,ICFL)-dy)
     &     +grz*(SFAREA(3,ICFL)-dz)
        dflxy(ICFL,ICOM)=dd0*gf1!+expfact*gf2)
!------------------------------------------
        enddo
        enddo
 500    enddo
!---------------------------
! --- Soret effect
!---------------------------
        if(.false.) then
        do 540 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle 
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        dd0=0.d0
        dd1=-1.d0
        do 541 ICOM=1,ncomp
        dd0=dd0+dflxy(ICFL,ICOM)*dble(act(icom))
        dd2=abs(dflxy(ICFL,ICOM))*dble(act(icom))
        if(dd2.gt.dd1) then
          dd1=dd2
          ICM=ICOM
        endif
  541   enddo
        dflxy(ICFL,ICM)=dflxy(ICFL,ICM)-dd0
        enddo
  540   enddo
! --- Soret effect
        elseif(.true.) then  !BackGrand gas method :BGCOMP(IMAT)
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle 
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        dd0=0.d0
        dd1=-1.d0
        do ICOM=1,ncomp
        dd0=dd0+dflxy(ICFL,ICOM)*dble(act(icom))
        enddo
        dflxy(ICFL,BGCOMP(IMAT))=-dd0
!        dflxy(ICFL,BGCOMP(IMAT))=dflxy(ICFL,BGCOMP(IMAT))-dd0
        enddo
        enddo
        endif
!---------------------------
! --- 
!---------------------------
        call bc_thmodif(mxcomp,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  dflxy)
!
! --- 
!
        do 550 IIMAT=1,NMAT   !ICF=1,NCVFAC
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        do 551 ICOM=1,ncomp
        dum1=(dflxy(ICFL,ICOM))*SFAREA(4,ICFL)
        dyc(ICVA,ICOM)=dyc(ICVA,ICOM)+dum1
        dyc(ICVB,ICOM)=dyc(ICVB,ICOM)-dum1
  551   enddo
        enddo
  550   enddo
      endif  !endif(ical_thmoDIFF==1)
!---------------------------
! --- temperature : dt/dx
!---------------------------
      call grad_cell(1,9,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,tmp,grdc)
!
      call dc_symprv
     &(1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & LCYCOLD,wifsld,OPPANG,
     & SFAREA,SFCENT,grdc,1,0)
!
      do 220 IIMAT=1,NMAT    !ICF=1,NCVFAC
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 220
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(IMAT>0) then
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        dd1=rmdc(ICVA)
        dd2=rmdc(ICVB)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dd0=dd1*dd2/(dd1*wi1+dd2*wi2+SML)
!
        dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)  
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl
        dlvect=abs(dx*SFAREA(1,ICFL)
     &            +dy*SFAREA(2,ICFL)
     &            +dz*SFAREA(3,ICFL))
        grx=wi1*grdc(ICVA,1)+wi2*grdc(ICVB,1)
        gry=wi1*grdc(ICVA,2)+wi2*grdc(ICVB,2)
        grz=wi1*grdc(ICVA,3)+wi2*grdc(ICVB,3)
        gf1=(tmp(ICVB)-tmp(ICVA))/dl   !(dlvect*dl+SML)
!-------------------------------------(1)
!        grdf(1)=SFAREA(1,ICFL)*gf1+grx*(SFAREA(1,ICFL)-dx/dlvect)
!        grdf(2)=SFAREA(2,ICFL)*gf1+gry*(SFAREA(2,ICFL)-dy/dlvect)
!        grdf(3)=SFAREA(3,ICFL)*gf1+grz*(SFAREA(3,ICFL)-dz/dlvect)
!        dflxt(ICFL)=dd0*
!     &        (SFAREA(1,ICFL)*grdf(1)
!     &        +SFAREA(2,ICFL)*grdf(2)
!     &        +SFAREA(3,ICFL)*grdf(3))
!-------------------------------------(2)
        gf2=grx*(SFAREA(1,ICFL)-dx)
     &     +gry*(SFAREA(2,ICFL)-dy)
     &     +grz*(SFAREA(3,ICFL)-dz)
!        dflxt(ICFL)=dd0*(gf1+expfact*gf2)
        dflxt(ICFL)=dd0*(gf1+gf2*dirr)
!
!-------------------------------------
!
        enddo
      elseif(IMAT.lt.0.and.iphs==1) then
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        dd1=rmdc(ICVA)
        dd2=rmdc(ICVB)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dd0=dd1*dd2/(dd1*wi1+dd2*wi2+SML)
        dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl
        dlvect=abs(dx*SFAREA(1,ICFL)
     &            +dy*SFAREA(2,ICFL)
     &            +dz*SFAREA(3,ICFL))
        grx=wi1*grdc(ICVA,1)+wi2*grdc(ICVB,1)
        gry=wi1*grdc(ICVA,2)+wi2*grdc(ICVB,2)
        grz=wi1*grdc(ICVA,3)+wi2*grdc(ICVB,3)
        gf1=(tmp(ICVB)-tmp(ICVA))/(dl+SML)   !(dlvect*dl+SML)
        gf2=grx*(SFAREA(1,ICFL)-dx)
     &     +gry*(SFAREA(2,ICFL)-dy)
     &     +grz*(SFAREA(3,ICFL)-dz)
        dflxt(ICFL)=dd0*(gf1)
        enddo
      endif
  220 enddo
!
      dflxh(:)=0.d0
!
! --- Mass fraction
!
!
      do 231 ICOM=1,ncomp
!
      call grad_cell(1,10,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,yys(:,ICOM),grdc)

!
      call dc_symprv
     &(1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  LCYCOLD,wifsld,OPPANG,
     &  SFAREA,SFCENT,grdc,1,0)
!
      do 230 IIMAT=1,NMAT    !ICF=1,NCVFAC
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(IMAT<0) then
        dflxy(ICFS:ICFE,:)=0.d0
        cycle
      endif
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
      dd1=rdsc(ICVA,ICOM)
      dd2=rdsc(ICVB,ICOM)
      wi1=wiface(ICFL)
      wi2=1.d0-wiface(ICFL)
      dd0=dd1*dd2/(dd1*wi1+dd2*wi2+SML)
!zhang-cvd
      dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
      dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
      dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
      dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)  
      dx=dx/dl
      dy=dy/dl
      dz=dz/dl
      dlvect=abs(dx*SFAREA(1,ICFL)
     &          +dy*SFAREA(2,ICFL)
     &          +dz*SFAREA(3,ICFL))
      grx=wi1*grdc(ICVA,1)+wi2*grdc(ICVB,1)
      gry=wi1*grdc(ICVA,2)+wi2*grdc(ICVB,2)
      grz=wi1*grdc(ICVA,3)+wi2*grdc(ICVB,3)
      gf1=(yys(ICVB,ICOM)-yys(ICVA,ICOM))/dl   !(dlvect*dl+SML)
!------------------------------------------(1)
!      grdf(1)=SFAREA(1,ICFL)*gf1+grx*(SFAREA(1,ICFL)-dx/dlvect)
!      grdf(2)=SFAREA(2,ICFL)*gf1+gry*(SFAREA(2,ICFL)-dy/dlvect)
!      grdf(3)=SFAREA(3,ICFL)*gf1+grz*(SFAREA(3,ICFL)-dz/dlvect)
!      dflxy(ICFL,ICOM)=dd0*
!     &      (SFAREA(1,ICFL)*grdf(1)
!     &      +SFAREA(2,ICFL)*grdf(2)
!     &      +SFAREA(3,ICFL)*grdf(3))
!------------------------------------------(2)
      gf2=grx*(SFAREA(1,ICFL)-dx)
     &   +gry*(SFAREA(2,ICFL)-dy)
     &   +grz*(SFAREA(3,ICFL)-dz)
!!!!!      dflxy(ICFL,ICOM)=dd0*(gf1+expfact*gf2)
      dflxy(ICFL,ICOM)=dd0*(gf1+gf2*dirr)   !coal
!!!!      dflxy(ICFL,ICOM)=dd0*(gf1)
!------------------------------------------
      enddo
  230 enddo
!
  231 enddo
!
!
!----------------------------------
! --- < 2.2 boundary conditions >--
!----------------------------------
      dsclt=0.d0
      dscly=0.d0
!
      call bc_hysdif(
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,mat_no,
     &  sfarea,cp,utau,frstcv,
     &  vel,tmp,yys,tmpbnd,yysbnd,htcbnd,mtcbnd,HTFLUX,
     &  rho,rmu,aks,rmut,
     &  dflxt,dflxy,dsclt,dscly)
!
!-------------------------------------
! --- < 2.3 modify diffusion flux >---
!-------------------------------------
!
      dflxh(:)=-1.d0       !temp-memory
      if(ical_VECT) then  !NOVECT
        tmpfac(:,1)=0.d0
        iwork1(:)=0
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICOM=1,ncomp
        do ICFL=ICFS,ICFE
        tmpfac(ICFL,1)=tmpfac(ICFL,1)+dflxy(ICFL,ICOM)
        if(abs(dflxy(ICFL,ICOM))>dflxh(ICFL)) then
          dflxh(ICFL)=abs(dflxy(ICFL,ICOM))
          iwork1(ICFL)=ICOM
        endif
        enddo
        enddo
!
        do ICFL=ICFS,ICFE
        dflxy(ICFL,iwork1(ICFL))=dflxy(ICFL,iwork1(ICFL))-tmpfac(ICFL,1)
        enddo
        enddo
      else
        if(.true.) then  !true for coal  !true YYYY1
        do 240 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        dd0=0.d0
        dd1=-1.d0
        do 241 ICOM=1,ncomp
        dd0=dd0+dflxy(ICFL,ICOM)*dble(act(icom))
        dd2=abs(dflxy(ICFL,ICOM))*dble(act(icom))
        if(dd2.gt.dd1) then
          dd1=dd2
          ICM=ICOM
        endif
  241   enddo
        dflxy(ICFL,ICM)=dflxy(ICFL,ICM)-dd0 
        enddo
  240   enddo
        elseif(.false.) then  !.false.  !false YYYY2
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle 
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        dd0=0.d0
        dd1=-1.d0
        do ICOM=1,ncomp
        dd0=dd0+dflxy(ICFL,ICOM)*dble(act(icom))
        enddo
        dflxy(ICFL,BGCOMP(IMAT))=dflxy(ICFL,BGCOMP(IMAT))-dd0
        enddo
        enddo
        endif
      endif
!
      if(ical_FC==PEFC) then
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdintr) then
        IIMAT=MAT_BCIDX(nb,1)
        IIMAT2=MAT_BCIDX(nb,2)
        IMAT=MAT_NO(IIMAT)
        IMAT2=MAT_NO(IIMAT2)
        IMAT_U=nofld(abs(IMAT))
        IMAT_U2=nofld(abs(IMAT2))
        if(IMAT_U==No_Mem.or.IMAT_U2==No_Mem) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          dd1=rdsc(ICV,vap_no)
          dd2=rdsc(ICVP,vap_no)
          dd0=dd1*dd2/(dd1+dd2+SML)
          Fs=0.5d0*(CVVOLM(ICVP)**(1/3.d0)+CVVOLM(ICV)**(1/3.d0))
          TC=yys(ICV,vap_no)
          TP=yys(ICVP,vap_no)
          dum1=dd0/Fs*(TP-TC)
          dflxy(ICFL,vap_no)= dum1
          dflxy(ICFP,vap_no)=-dum1
          enddo
        endif
        endif
        enddo
      endif
!
!
!
!---------------------------------
!--< 2.4 diffusion of enthalpy >--
!---------------------------------
!
      dflxh=0.d0
      if(ical_vect) then   !NOVECT

      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICOM=1,ncomp
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
      dflxh(ICFL)=dflxh(ICFL)
     &   +max(0.d0,dflxy(ICFL,ICOM))*hhs(ICVB,ICOM)
     &   +min(0.d0,dflxy(ICFL,ICOM))*hhs(ICVA,ICOM)
      
      enddo
      enddo
      enddo

      else

      if(.true.) then   !true  YYYY3
      do 250 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)

      do ICFL=ICFS,ICFE     !vector
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
!      dd0=0.d0
      do 251 ICOM=1,ncomp
!      dd0=dd0+
!     &   +max(0.d0,dflxy(ICFL,ICOM))*hhs(ICVA,ICOM)*dble(act(icom))
!     &   +min(0.d0,dflxy(ICFL,ICOM))*hhs(ICVB,ICOM)*dble(act(icom))
      dflxh(ICFL)=dflxh(ICFL)
     &   +max(0.d0,dflxy(ICFL,ICOM))*hhs(ICVB,ICOM)*dble(act(icom))
     &   +min(0.d0,dflxy(ICFL,ICOM))*hhs(ICVA,ICOM)*dble(act(icom))
  251 enddo
!      dd0=dd0+dflxt(ICFL)
!      ICM=BGCOMP(IMAT)
!      dflxh(ICFL)=dd0-
!     &   +max(0.d0,dflxy(ICFL,ICM))*hhs(ICVA,ICM)*dble(act(icm))
!     &   +min(0.d0,dflxy(ICFL,ICM))*hhs(ICVB,ICM)*dble(act(icm))
      enddo
  250 enddo
      endif
!
      endif
!--------------------------------------------------------------
!      call bc_thmodif(1,
!     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,                        
!     &  dflxh)
!--------------------------------------------------------------
!------------------------------
!--< 2.5 interface boundary >--
!------------------------------
!
      call grad_cell(1,11,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,tmp,grdc)
!
      call bc_intrfc(iphs,
     &  MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  SFAREA,CVCENT,CVVOLM,SFCENT,
     &  rho,tmp,aks,grdc,rmdc,htcbnd,radbnd,tmpbnd,HTFLUX,rmu,dflxh,
     &  dflxt,dsclt) 
!--------------------------------------------------------------
!-< 3. Take thermal conductivity & diffusion term into count >-
!--------------------------------------------------------------
      dhc(:)=0.d0
      if(ical_vect) then   !NOVECT
        do IE=1,MAXIE
        DO ICVL=ICVS_V,ICVE_V
        IF(ABS(vctr(ICVL,IE))>0) then
        dhc(ICVL)=dhc(ICVL)
     &       -sign(1,vctr(ICVL,IE))
     &       *(dflxt(abs(vctr(ICVL,IE)))+dflxh(abs(vctr(ICVL,IE))))
     &       *SFAREA(4,abs(vctr(ICVL,IE)))
        endif
        enddo
        enddo
!
        do ICOM=1,ncomp
        do IE=1,MAXIE
      
        DO ICVL=ICVS_V,ICVE_V
        IF(ABS(vctr(ICVL,IE))>0) then
          dyc(ICVL,ICOM)=dyc(ICVL,ICOM)
     &        -sign(1,vctr(ICVL,IE))
     &       *dflxy(abs(vctr(ICVL,IE)),ICOM)
     &       *SFAREA(4,abs(vctr(ICVL,IE)))
        endif
        enddo
        enddo
        enddo
      else
        do 300 IIMAT=1,NMAT   !ICF=1,NCVFAC 
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle 
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE 
        ICVA=LVEDGE(1,ICFL) 
        ICVB=LVEDGE(2,ICFL) 
        dum1=(dflxt(ICFL)+dflxh(ICFL))*SFAREA(4,ICFL) 
        
        dhc(ICVA)=dhc(ICVA)+dum1
        dhc(ICVB)=dhc(ICVB)-dum1
        if(IMAT<0) cycle
        do 301 ICOM=1,ncomp
        dum1=(dflxy(ICFL,ICOM))*SFAREA(4,ICFL)
        dyc(ICVA,ICOM)=dyc(ICVA,ICOM)+dum1
        dyc(ICVB,ICOM)=dyc(ICVB,ICOM)-dum1
  301   enddo
        enddo
  300   enddo
!
        if(.false.) then    !true   !false   !YYYY4
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle 
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dd0=0.d0
        do ICOM=1,ncomp
!        if(ICOM==BGCOMP(IMAT)) cycle
!!!!!!!!!!!        dum1=dyc(ICVL,ICOM)*hhs(ICVL,ICOM)
        dd0=dd0+dyc(ICVL,ICOM)
!!!!!!!!!!!        dhc(ICVL)=dhc(ICVL)+dum1
        enddo
        dyc(ICVL,BGCOMP(IMAT))=dyc(ICVL,BGCOMP(IMAT))-dd0
        enddo
        enddo
        endif
      endif
!
      if(ical_vect) then   !NOVECT
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
           dhdt(ICVL)=dhdt(ICVL)+dhc(ICVL)
        enddo
!
        do ICOM=1,ncomp
        do ICVL=ICVS,ICVE
        dydt(ICVL,ICOM)=dydt(ICVL,ICOM)+dyc(ICVL,ICOM)
        enddo
        enddo
        enddo
      else
        do 310 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) goto 310
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        dhdt(ICVL)=dhdt(ICVL)+dhc(ICVL)
        do 311 ICOM=1,ncomp
        dydt(ICVL,ICOM)=dydt(ICVL,ICOM)+dyc(ICVL,ICOM)
 311    enddo
 330    enddo
 310    enddo
      endif
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (ncomp,MXALLCV,NCV,dydt)
      ENDIF
!
      deallocate(dyc)
!
      return
!
      end subroutine hys_diff
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine vel_diff(iphs,ismpl,iter,rdeltt,
     &  LVEDGE,LBC_SSF,LCYCSF,kdbv,LFUTAU,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,FRSTCV,
     &  grdc,dflx,aks,
     &  LCYCOLD,wifsld,OPPANG,vctr,
     &  rmu,rmut,rho,vel,utau,dsclv,dvdt,rmue,ICODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       dflx  <=  rvd 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. Calculate viscous term
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only         : dflxTmpFFlag,lenfnm
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX,
     &                             kdprdc,kdsymm,kdilet,kdolet,kvnslp,
     6                             kdtchi
      use module_Euler2ph,only   : ieul2ph
      use module_scalar,  only   : iaph,ical_cavi,ivold
      use module_model,only      : ical_vect,nthrds
      use module_vector,only     : ICVS_V,ICVE_V,
     &                             ICFS_V,ICFE_V,
     &                             ICVSIN_V,ICVEIN_V,
     &                             IDCS_V,IDCE_V,index_c,index_f
      use module_model,only      : expfact
      use module_fluidforce,only : F_fric,Nforce,ForceFlag,F_NUM,
     &                             NO_F_WALL
      use module_metrix,only     : DW2K_VECT,
     &                             tmpfac=>d2vect,
     &                             difv=>d2work4
!
      implicit none
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: rdeltt
      integer,intent(in)    :: ismpl,ICODE,iter,iphs
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF( MXSSFBC)
      integer,intent(in)    :: LCYCSF(  MXSSFBC)
      integer,intent(in)    :: LFUTAU(  MXCV   )
      integer,intent(in)    :: kdbv  (  MXCVFAC)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3,MXCVFAC)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(  MXALLCV)
      real*8 ,intent(in)    :: FRSTCV(  MXSSFBC)
      real*8 ,intent(inout) :: grdc  (  MXALLCV,3,3)
      real*8 ,intent(inout) :: dflx  (  MXCVFAC,3)
      real*8 ,intent(inout) :: rmu   (  MXALLCV)
      real*8 ,intent(inout) :: rmut  (  MXALLCV)
      real*8 ,intent(in)    :: rho   (  MXALLCV)
      real*8 ,intent(inout) :: utau  (  0:MXSSFBC)
      real*8 ,intent(inout) :: vel   (  MXALLCV,3,2)
      real*8 ,intent(inout) :: dvdt  (  MXALLCV,3)
      real*8 ,intent(out)   :: rmue  (  MXALLCV)
      integer,intent(in)    :: LCYCOLD (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld  (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG  (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(in)    :: aks   (MXALLCVR,mxrans)
      real*8 ,intent(inout)    :: dsclv     (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: dd1,dd2,dum1,rmuf,UCONV,dxx,dyx,dzx,grdf(3)
      integer :: i,j,k,l,m,n,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD,IV,IE
      integer :: IMAT_U,II,IS
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,nbs
      real*8  :: wi1,wi2,grx,gry,grz,dx,dy,dz,dl,dlvect,dll
      real*8  :: dum,dum11,dum22,dum0,gf1,gf2,gf3,grdf1(3),grdf2(3,3)
      integer :: IDIM,ierr=0
      character(lenfnm) :: fnam1,fnam2
!
! --- 
!
!      if(ieul2ph.eq.1) then
!        rmu2(:)=rmu(:)
!        rmut2(:)=rmut(:)
!        do ICV=1,NALLCV
!        rmu(ICV)=rmu2(ICV)/rho(ICV)   !ALP(ICV)
!        rmut(ICV)=rmut2(ICV)/rho(ICV) !*ALP(ICV)
!        enddo
!      endif
!
!-< 1. Calculate gradient of u,v,w >-
!
!--< 1.1 gradient at cell center >--
!
      call grad_cell(3,12,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vel(:,:,1),grdc)
!
!--< 1.2 set dummy cell >--
!
      call dc_vgrad(MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,
     &              LCYCOLD,wifsld,OPPANG,
     &              grdc)
!-------------------------------------------
!-< 2. Calculate viscous flux >-
!-------------------------------------------
      if(ical_vect) then   !NOVECT
!CIDR NODEP
        do ICVL=ICVS_V,ICVE_V   !index_c(myid)+1,index_c(myid+1)
        rmue(ICVL)=rmu(ICVL)+rmut(ICVL)
        enddo
!CIDR NODEP
        DO ICVL=IDCS_V,IDCE_V   !index_c(myid)+1,index_c(myid+1)
        rmue(ICVL)=rmu(ICVL)+rmut(ICVL)
        enddo
      else
        do 200 IIMAT=1,NMAT    !
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        rmue(ICVL)=rmu(ICVL)+rmut(ICVL)

 330    continue
        ICVS=MAT_DCIDX(IIMAT-1)+1
        ICVE=MAT_DCIDX(IIMAT)
        do 332 ICVL=ICVS,ICVE
        rmue(ICVL)=rmu(ICVL)+rmut(ICVL)
 332    enddo
 200    enddo
      endif
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &                mat_cal,rmue)
!--------------------------
! --- BC Mue 
!--------------------------
      call bc_veldif(
     &  LVEDGE,LBC_SSF,LCYCSF,LFUTAU,mat_cal,MAT_NO,
     &  SFAREA,CVCENT,FRSTCV,vel,rho,rmu,rmut,utau,dsclv,
     &  dflx,rmue,1)
!-------------------------------
!--< 2.2 boundary conditions >--
!-------------------------------
! --- diffusion flux 
!-------------------------------
      if(ical_vect) then   !NOVECT
        do IIMAT=1,NMAT 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl
        tmpfac(ICFL,1)=dl*abs(dx*SFAREA(1,ICFL)  !dlvect
     &                       +dy*SFAREA(2,ICFL)
     &                       +dz*SFAREA(3,ICFL))
        enddo
        enddo
!
        do IIMAT=1,NMAT 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do 415 IDIM=1,3 
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl


        rmuf=rmue(ICVA)*rmue(ICVB)/
     6      (rmue(ICVA)*wiface(ICFL)
     &      +rmue(ICVB)*(1.d0-wiface(ICFL))+SML)

        grx=      wiface(ICFL) *grdc(ICVA,1,IDIM)
     &     +(1.d0-wiface(ICFL))*grdc(ICVB,1,IDIM)
        gry=      wiface(ICFL) *grdc(ICVA,2,IDIM)
     &     +(1.d0-wiface(ICFL))*grdc(ICVB,2,IDIM)
        grz=      wiface(ICFL) *grdc(ICVA,3,IDIM)
     &     +(1.d0-wiface(ICFL))*grdc(ICVB,3,IDIM)
        gf1=(vel(ICVB,IDIM,1)-vel(ICVA,IDIM,1))/tmpfac(ICFL,1)
        gf2=grx*(SFAREA(1,ICFL)-dx)
     &     +gry*(SFAREA(2,ICFL)-dy)
     &     +grz*(SFAREA(3,ICFL)-dz)
        dflx(ICFL,IDIM)=rmuf*(gf1+gf2)
        enddo
 415    enddo
        enddo

        do IIMAT=1,NMAT 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        rmuf=rmue(ICVA)*rmue(ICVB)/
     6      (rmue(ICVA)*wiface(ICFL)
     &      +rmue(ICVB)*(1.d0-wiface(ICFL))+SML)
        grdf2(1,1)=wiface(ICFL) *grdc(ICVA,1,1)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,1,1)
        grdf2(2,1)=wiface(ICFL) *grdc(ICVA,2,1)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,2,1)
        grdf2(3,1)=wiface(ICFL) *grdc(ICVA,3,1)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,3,1)

        grdf2(1,2)=wiface(ICFL) *grdc(ICVA,1,2)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,1,2)
        grdf2(2,2)=wiface(ICFL) *grdc(ICVA,2,2)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,2,2)
        grdf2(3,2)=wiface(ICFL) *grdc(ICVA,3,2)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,3,2)

        grdf2(1,3)=wiface(ICFL) *grdc(ICVA,1,3)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,1,3)
        grdf2(2,3)=wiface(ICFL) *grdc(ICVA,2,3)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,2,3)
        grdf2(3,3)=wiface(ICFL) *grdc(ICVA,3,3)
     &      +(1.d0-wiface(ICFL))*grdc(ICVB,3,3)


        dflx(ICFL,1)= dflx(ICFL,1)+
     &   rmuf*(SFAREA(1,ICFL)*grdf2(1,1)
     &        +SFAREA(2,ICFL)*grdf2(1,2)
     &        +SFAREA(3,ICFL)*grdf2(1,3))
        dflx(ICFL,2)= dflx(ICFL,2)+
     &   rmuf*(SFAREA(1,ICFL)*grdf2(2,1)
     &        +SFAREA(2,ICFL)*grdf2(2,2)
     &        +SFAREA(3,ICFL)*grdf2(2,3))
        dflx(ICFL,3)= dflx(ICFL,3)+
     &   rmuf*(SFAREA(1,ICFL)*grdf2(3,1)
     &        +SFAREA(2,ICFL)*grdf2(3,2)
     &        +SFAREA(3,ICFL)*grdf2(3,3))
        enddo
        enddo
!------------------
! --- Sclar
!------------------
      else
      do 221 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) goto 221
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
      dd1=rmue(ICVA)
      dd2=rmue(ICVB)
      wi1=wiface(ICFL)
      wi2=1.d0-wiface(ICFL)
      rmuf=dd1*dd2/(dd1*wi1+dd2*wi2+SML)
!      rmuf=dd1*wi1+dd2*wi2
      dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
      dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
      dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
      dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
      dx=dx/dl
      dy=dy/dl
      dz=dz/dl
      dlvect=abs(dx*SFAREA(1,ICFL)
     &          +dy*SFAREA(2,ICFL)
     &          +dz*SFAREA(3,ICFL))
      do 405 IDIM=1,3  !u v w 
      grx=wi1*grdc(ICVA,1,IDIM)+wi2*grdc(ICVB,1,IDIM)
      gry=wi1*grdc(ICVA,2,IDIM)+wi2*grdc(ICVB,2,IDIM)
      grz=wi1*grdc(ICVA,3,IDIM)+wi2*grdc(ICVB,3,IDIM)
      gf1=(vel(ICVB,IDIM,1)-vel(ICVA,IDIM,1))/(dlvect*dl)
!------------------------------------- (1)
!      grdf(1)=grx*(SFAREA(1,ICFL)-dx)
!      grdf(2)=gry*(SFAREA(2,ICFL)-dy)
!      grdf(3)=grz*(SFAREA(3,ICFL)-dz)
!      gf2=
!     &      (SFAREA(1,ICFL)*grdf(1)
!     &      +SFAREA(2,ICFL)*grdf(2)
!     &      +SFAREA(3,ICFL)*grdf(3))
!------------------------------------- (2)
      gf2=grx*(SFAREA(1,ICFL)-dx)
     &   +gry*(SFAREA(2,ICFL)-dy)
     &   +grz*(SFAREA(3,ICFL)-dz)
!-------------------------------------
      grdf1(IDIM)=gf1+gf2!*expfact
      grdf2(1,IDIM)=grx
      grdf2(2,IDIM)=gry
      grdf2(3,IDIM)=grz
 405  continue
      do IDIM=1,3
      dflx(ICFL,IDIM)=rmuf*grdf1(IDIM)+
     &   rmuf*(SFAREA(1,ICFL)*grdf2(IDIM,1)
     &        +SFAREA(2,ICFL)*grdf2(IDIM,2)
     &        +SFAREA(3,ICFL)*grdf2(IDIM,3))
      ENDDO
      enddo
 221  continue
      endif
!------------------
! --- BC: dflx
!------------------
      call bc_veldif(
     &  LVEDGE,LBC_SSF,LCYCSF,LFUTAU,mat_cal,MAT_NO,
     &  SFAREA,CVCENT,FRSTCV,vel,rho,rmu,rmut,utau,dsclv,
     &  dflx,rmue,2)
!-------------------------------------
!-< 3. Take viscous term into count >-
!-------------------------------------new
      if(ical_vect) then   !NOVECT
        do IE=1,MAXIE
        DO ICVL=ICVS_V,ICVE_V
        IF(ABS(vctr(ICVL,IE))>0) then
        dvdt(ICVL,1)=dvdt(ICVL,1)
     &              -sign(1,vctr(ICVL,IE))
     &              *(dflx(abs(vctr(ICVL,IE)),1)
     &              *SFAREA(4,abs(vctr(ICVL,IE))))
        dvdt(ICVL,2)=dvdt(ICVL,2)
     &              -sign(1,vctr(ICVL,IE))
     &              *(dflx(abs(vctr(ICVL,IE)),2)
     &              *SFAREA(4,abs(vctr(ICVL,IE))))
        dvdt(ICVL,3)=dvdt(ICVL,3)
     &              -sign(1,vctr(ICVL,IE))
     &              *(dflx(abs(vctr(ICVL,IE)),3)
     &              *SFAREA(4,abs(vctr(ICVL,IE))))
        endif
        enddo
        enddo
!!$omp end parallel do
      else
        do 300 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 300
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        do 301 IV=1,3
        dum1=dflx(ICFL,IV)*SFAREA(4,ICFL)
        dvdt(ICVA,IV)=dvdt(ICVA,IV)+dum1
        dvdt(ICVB,IV)=dvdt(ICVB,IV)-dum1
 301    continue
        enddo
 300    continue
      endif
!
      if(ForceFlag) then
        F_fric(:,:)=0.d0
        do I=1,Nforce
        IS=F_NUM(I-1)+1
        IE=F_NUM(I)
        do II=IS,IE
        nb=NO_F_WALL(II)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        F_fric(:,I)=F_fric(:,I)-dflx(ICFL,:)*SFAREA(4,ICFL)
        enddo
        enddo
        enddo
      endif
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (3,MXALLCV,NCV,dvdt(:,1:3))
      ENDIF
!
!      if(ieul2ph.eq.1) then
!        rmu(:)=rmu2(:)
!        rmut(:)=rmut2(:)
!      endif
!

      return
      end subroutine vel_diff
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine vel_diff_e2p(iphs,ismpl,iter,rdeltt,
     &  LVEDGE,LBC_SSF,LCYCSF,kdbv,LFUTAU,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,FRSTCV,
     &  grdc,dflx,aks,
     &  LCYCOLD,wifsld,OPPANG,vctr,
     &  rmu,rmut,rho,vel,utau,dsclv,dvdt,rmue,ICODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       dflx  <=  rvd 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. Calculate viscous term
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only : dflxTmpFFlag,lenfnm
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX,
     &                             kdprdc,kdsymm,kdilet,kdolet,kvnslp,
     6                             kdtchi
      use module_Euler2ph,only   : ieul2ph
      use module_scalar,  only   : iaph,ical_cavi,ivold
      use module_model,only  : ical_vect,nthrds
      use module_vector,only : ICVS_V,ICVE_V,
     &                         ICFS_V,ICFE_V,
     &                         ICVSIN_V,ICVEIN_V,
     &                         IDCS_V,IDCE_V,index_c,index_f
      use module_model,only  : expfact
!
      implicit none
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: rdeltt
      integer,intent(in)    :: ismpl,ICODE,iter,iphs
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF( MXSSFBC)
      integer,intent(in)    :: LCYCSF(  MXSSFBC)
      integer,intent(in)    :: LFUTAU(  MXCV   )
      integer,intent(in)    :: kdbv  (  MXCVFAC)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3,MXCVFAC)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(  MXALLCV)
      real*8 ,intent(in)    :: FRSTCV(  MXSSFBC)
      real*8 ,intent(inout) :: grdc  (  MXALLCV,3,3)
      real*8 ,intent(inout) :: dflx  (  MXCVFAC,3)
      real*8 ,intent(inout) :: rmu   (  MXALLCV)
      real*8 ,intent(inout) :: rmut  (  MXALLCV)
      real*8 ,intent(in)    :: rho   (  MXALLCV)
      real*8 ,intent(inout) :: utau  (  0:MXSSFBC)
      real*8 ,intent(inout) :: vel   (  MXALLCV,3,2)
      real*8 ,intent(inout) :: dvdt  (  MXALLCV,3)
      real*8 ,intent(out)   :: rmue  (  MXALLCV)
      integer,intent(in)    :: LCYCOLD (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld  (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG  (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(in)    :: aks   (MXALLCVR,mxrans)
      real*8 ,intent(inout)    :: dsclv     (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: dd1,dd2,dum1,rmuf,UCONV,dxx,dyx,dzx,grdf(3)
      integer :: i,j,k,l,m,n,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD,IV,IE
      integer :: IMAT_U
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,nbs
      real*8  :: wi1,wi2,grx,gry,grz,dx,dy,dz,dl,dlvect,dll
      real*8  :: dum,dum11,dum22,dum0,gf1,gf2,gf3,grdf1(3),grdf2(3,3)
      integer :: IDIM,ierr=0
      character(lenfnm) :: fnam1,fnam2
!
! --- 
!
!-< 1. Calculate gradient of u,v,w >-
!
!--< 1.1 gradient at cell center >--
!
      call grad_cell(3,12,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vel(:,:,1),grdc)
!
!--< 1.2 set dummy cell >--
!
      call dc_vgrad(MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,
     &              LCYCOLD,wifsld,OPPANG,
     &              grdc)
!-------------------------------------------
!-< 2. Calculate viscous flux >-
!-------------------------------------------
      do 200 IIMAT=1,NMAT    !
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        rmue(ICVL)=(rmu(ICVL)+rmut(ICVL))
 330    continue
        ICVS=MAT_DCIDX(IIMAT-1)+1
        ICVE=MAT_DCIDX(IIMAT)
        do 332 ICVL=ICVS,ICVE
        rmue(ICVL)=rmu(ICVL)+rmut(ICVL)
 332    enddo
 200  enddo
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &                mat_cal,rmue)
!--------------------------
! --- BC Mue
!--------------------------
      call bc_veldif(
     &  LVEDGE,LBC_SSF,LCYCSF,LFUTAU,mat_cal,MAT_NO,
     &  SFAREA,CVCENT,FRSTCV,vel,rho,rmu,rmut,utau,dsclv,
     &  dflx,rmue,1)
!-------------------------------
!--< 2.2 boundary conditions >--
!-------------------------------
! --- diffusion flux
!-------------------------------
      do 221 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) goto 221
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
      dd1=rmue(ICVA)
      dd2=rmue(ICVB)
      wi1=wiface(ICFL)
      wi2=1.d0-wiface(ICFL)
      rmuf=dd1*dd2/(dd1*wi1+dd2*wi2+SML)
      dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
      dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
      dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
      dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
      dx=dx/dl
      dy=dy/dl
      dz=dz/dl
      dlvect=abs(dx*SFAREA(1,ICFL)
     &          +dy*SFAREA(2,ICFL)
     &          +dz*SFAREA(3,ICFL))
      do 405 IDIM=1,3  !u v w 
      grx=wi1*grdc(ICVA,1,IDIM)+wi2*grdc(ICVB,1,IDIM)
      gry=wi1*grdc(ICVA,2,IDIM)+wi2*grdc(ICVB,2,IDIM)
      grz=wi1*grdc(ICVA,3,IDIM)+wi2*grdc(ICVB,3,IDIM)
      gf1=(vel(ICVB,IDIM,1)-vel(ICVA,IDIM,1))/(dlvect*dl)
!-------------------------------------(2)
      gf2=grx*(SFAREA(1,ICFL)-dx)
     &   +gry*(SFAREA(2,ICFL)-dy)
     &   +grz*(SFAREA(3,ICFL)-dz)
!-------------------------------------
      grdf1(IDIM)=gf1+gf2!*expfact
      grdf2(1,IDIM)=grx
      grdf2(2,IDIM)=gry
      grdf2(3,IDIM)=grz
 405  continue
      do IDIM=1,3
      dflx(ICFL,IDIM)=rmuf*grdf1(IDIM)+
     &   rmuf*(SFAREA(1,ICFL)*grdf2(IDIM,1)
     &        +SFAREA(2,ICFL)*grdf2(IDIM,2)
     &        +SFAREA(3,ICFL)*grdf2(IDIM,3))
      ENDDO
      enddo
 221  continue
!------------------
! --- BC: dflx
!------------------
      call bc_veldif(
     &  LVEDGE,LBC_SSF,LCYCSF,LFUTAU,mat_cal,MAT_NO,
     &  SFAREA,CVCENT,FRSTCV,vel,rho,rmu,rmut,utau,dsclv,
     &  dflx,rmue,2)
!-------------------------------------
!-< 3. Take viscous term into count >-
!-------------------------------------
      do 300 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 300
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dum1=aks(ICVA,iaph(iphs))*wi1+aks(ICVB,iaph(iphs))*wi2
        do 301 IV=1,3
        dum1=dum1*dflx(ICFL,IV)*SFAREA(4,ICFL)
        dvdt(ICVA,IV)=dvdt(ICVA,IV)+dum1
        dvdt(ICVB,IV)=dvdt(ICVB,IV)-dum1
 301    continue
        enddo
 300  continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (3,MXALLCV,NCV,dvdt(:,1:3))
      ENDIF
!
      return
      end subroutine vel_diff_e2p
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
