!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!      subroutine METR_FFLOW
!      subroutine metric_admin()
!      subroutine metric_intrpw()
!      subroutine list_fluid
!      subroutine list_bcin
!      subroutine wall_distance
!      subroutine read_geom
!      subroutine list_sliding
!      subroutine allocat_array
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine METR_FFLOW(
     &  LVEDGE,LVRTCV,LFUTAU,
     &  LBC_SSF,LCYCSF,FRSTCV,DISINL,
     &  SFAREA,SFCENT,CVVOLM,CVVOL0,CVCENT,wiface,DISALL,xta,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
!     &  LCV_CV,vctr,vctrp,NEIGHBOR,DNORM,FIELD_U,
     &  LCV_CV,vctrp,NEIGHBOR,DNORM,FIELD_U,
     &  Pstart,LCYCOLD,wifsld,OPPANG,locmsh,
     &  KDBF,grdc,
     &  ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension
      USE module_dimension,only : IVDIM
      use module_hpcutil
      use module_io,       only : ifli,ifll,ifle,gdScale
      use module_nowtime,  only : time,iter
      use module_vector   ,ONLY : input_vect    =>inputdata 
! 
      use module_vector   ,ONLY : MXUP,MXLW,N2
      use module_metrix,only    : msk,W2K1,IW2K2,
     &                            ipx,aax,W2K2,d1work2,W1K6,
     &                            iwork4,iwork5,iwork6,iwork7,d2vect,
     &                            IW2K_VECT,DW2K_VECT,J_SRC,iwork1,
     &                            TH_DIFF,tmpsld!,mask
!      use module_metrix,only    : IQ =>IW2K1
      use module_metrix,only    : movflg,cord,lacell,lvcell,
     &                            lvface,lbface,lcface,lfcell,LEFACE,
     &                            area,volume,gface,gcell,listbc
      use module_metrix,only    : bdyfrc,bdyf,vctr
      use module_metrix,only    : vflc,yplusf,rva_s
      use module_boundary,only  : SFBOUN,NBOUND,kdshutr,
     &                            IFFACE,NFBOUN,NBFS,nbcnd,kdbcnd,
     &                            pkdwall,kxnone,kdintr,kdfire,kdcvd,
     &                            kdsymm,kdovst,kdsld,kdbuff,mat_bcidx
      use module_model, only    : ical_mvmsh,iLBF_P,ical_thmoDIFF
      use module_model, only    : PEFC,dles,icaltb,monitor_stp
      use module_model,only     : ical_vect,nthrds,ical_prt,u_func,
     &                            drag_coef
      use module_scalar,ONLY    : ical_FC,x1max,x1min,
     &                                    x2max,x2min,x3max,x3min
      use module_scalar,ONLY    : idifflm,ORG,comb_mdl,
     &                                    ixi,ixi_2,ixi_x,ixi_c
      use module_rans,only    : akslw,akshg
      use module_chemreac,only  : iBV,nneq
!      use module_particle,only  : nneqP
      use module_flags,    only : euleri,intgvv,intgty,intgke
      use module_Dynamic_SGS,only : lij,mij
      use module_vof,only : ical_vof
      use module_metrix,only    : t_dir,ANGBND,STA_SAV,lreac
      use module_les  ,only : n_ave,iuvwc_re,
     &                        iuvw_ave_rms_re,
     &                        ip_ave,it_ave,imu_ave,
     &                        icomp_ave,irans_ave,
     &                        it_rms,ip_rms,imu_rms,
     &                        icomp_rms,irans_rms,
     &                        iuvwt_re,imaxmin,
     &                        iuvw_rans_re,ista_momery
      use module_radsxf,only : RadHeat,RadHeatFlux,radinten,
     &                         P_rad,PAbsorb,PScatter
      use module_material,only : radfludflag
      use module_rad   , only  : radflag,ndiv
      use module_metrix,only   : SHUTFL
      use module_material,only : ical_sld,rg_mark,lclsd,iclosd
      use module_metrix,only   : matsld
      use module_movegrid,only : ipiston
      use module_metrix,only   : WDRBND,multiR
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(out)      :: ierror
!
      integer,intent(out)      :: LVEDGE(2,  MXCVFAC)
      integer,intent(out)      :: LVRTCV(    MXALLCV)
      integer,intent(out)      :: LFUTAU(    MXCV)
      integer,intent(out)      :: MAT_NO(0:MXMAT)
      integer,intent(out)      :: LCV_CV(MXALLCV)
      integer,intent(out)      :: MAT_CV(MXALLCV)
      integer,intent(out)      :: MAT_INDEX( 0:MXMAT)
      integer,intent(out)      :: MAT_CVEXT( 0:MXMAT)
      integer,intent(out)      :: MAT_DCIDX( 0:MXMAT)
      integer,intent(out)      :: MAT_CFIDX( 0:MXMAT)
      logical,INTENT(out)      :: mat_cal  ( 0:MXMAT)
      integer,intent(out)      :: Pstart(      MXMAT)
      integer,intent(out)      :: LBC_SSF(   mxssfbc)
      integer,intent(inout)    :: LCYCSF(    mxssfbc)
      integer,intent(inout)    :: locmsh(    mxssfbc)
      integer,intent(inout)    :: LCYCOLD(MXSSFBC_SLD)
!
!      integer,intent(inout)    :: vctr(MXCV_V,0:MXBND_V)
      integer,intent(inout)    :: vctrp(MXCV_VP,0:MXBND_VP)
      integer,intent(inout)    :: NEIGHBOR(MEP,MXALLCV_P)
      real*8 ,intent(inout)    :: DNORM   (MXCVFAC_P)
!
      real*8 ,intent(out)      :: FIELD_U(MXCV_D,NFLID)
      real*8 ,intent(out)      :: grdc   (   MXALLCV,3)
!
      integer,intent(out)      :: KDBF  (    MXCVFAC)
      real*8 ,intent(out)      :: FRSTCV(    mxssfbc)
      real*8 ,intent(out)      :: DISINL(    MXSSFBC)
!
      real*8 ,intent(out)      :: SFAREA(4,  MXCVFAC)
      real*8 ,intent(out)      :: SFCENT(3,  MXCVFAC)
      real*8 ,intent(out)      :: CVVOLM(    MXALLCV)
      real*8 ,intent(out)      :: CVVOL0(    MXALLCV)
      real*8 ,intent(out)      :: CVCENT(3,  MXALLCV)
      real*8 ,intent(out)      :: wiface(    MXCVFAC)
      real*8 ,intent(out)      :: wifsld(MXSSFBC_SLD)
      real*8 ,intent(out)      :: OPPANG(MXSSFBC_SLD)
      real*8 ,intent(out)      :: DISALL(    MXCV)
      real*8 ,intent(out)      :: xta   (    MXCVFAC)
!
      real*8  :: delt_tmp=-1.d0
      integer :: icom,i,nrnsx,shutt=0,kd,nb,NSLD=0,KMAT,iimat
      integer :: IIMAT1,IIMAT2,MATSAM1,MATSAM2,idum1,idum2
      real*8  :: xyz(3),dum1
      integer :: IMAT,ICVL,ICVS,ICVE

!      integer,save,allocatable   :: sldbcmat (:,:)
! --- ZERO Initialization
!
!      allocate(sldbcmat(nbcnd,2))
!      sldbcmat(:,:)=0
      LVEDGE=0
      LVRTCV=0
      LFUTAU=0      
!
      LCYCSF=0
      LCYCOLD=0
      LBC_SSF=0
      FRSTCV=0.d0
      DISINL=0.d0
      wifsld=1.d0
      OPPANG=0.d0
!
      SFAREA=0.d0
      SFCENT=0.d0
      CVVOLM=0.d0
      CVVOL0=0.d0
      CVCENT=0.d0
      wiface=0.d0
      DISALL=0.d0
      xta   =0.d0
      locmsh=0
!
!
!
       MXALLCV_RAD=1
!       Ngauss_R=1
!       NDIV_R=1
       NDNG=1
!       MXALLNG=1
       IF(radflag.EQ.1) THEN
         MXALLCV_RAD=MXALLCV
         IF(rg_mark.eq.1) THEN
           NDNG=Ngauss_R*NDIV_R
         else
           NDNG=NDIV
         ENDIF
       endif
       ALLOCATE(P_rad(MXALLCV_RAD,2),stat=ierror)
       ALLOCATE(RadHeat(MXALLCV_RAD),stat=ierror)
       ALLOCATE(RadHeatFlux(MXALLCV_RAD),stat=ierror)
       ALLOCATE(PAbsorb(MXALLCV_RAD),stat=ierror)
       ALLOCATE(PScatter(MXALLCV_RAD),stat=ierror)
!       ALLOCATE(radinten(NDNG,MXALLCV_RAD),stat=ierror)

       P_rad(:,:)=0.d0
       RadHeat=0.d0
       RadHeatFlux=0.d0
       PAbsorb=0.d0
       PScatter=0.d0
!       radinten=0.d0
!
!-< I. Pre-process >-
!
!
!-< 1. Set up topological data >-
!
!--< 1.1 input data >--
!
      if(NPE==1) then
!-----------------------------
! --- Input Serial data
!-----------------------------
        call read_geom (
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,
     &  LBC_SSF,LCYCSF,LCV_CV,locmsh,
     &  LVEDGE,LVRTCV,SFAREA,SFCENT,CVVOLM,CVCENT)
!
      elseif(NPE>1) then
!-----------------------------
! --- Input HPC data
!-----------------------------
        call read_geom(
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,
     &  LBC_SSF,LCYCSF,LCV_CV,locmsh,
     &  LVEDGE,LVRTCV,SFAREA,SFCENT,CVVOLM,CVCENT)
      else
        write(ifle,*) '    ###    ERR : CPU number = ',NPE
        call FFRABORT(1,'METR_FFLOW')
      endif
!
      allocate(bdyf(MXALLCV,3),stat=ierror)
        if(ierror/=0) call 
     &  FFRABORT(1,'ERR:allocate bdyf in METR_FFLOW')
        bdyf(:,:)=0.d0
!-----------------------------
! --- 
!-----------------------------
!      call corr_area(1,
!     &  MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,
!     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,grdc)
!-----------------------------
! --- 
!-----------------------------
!-----------------------------
! --- 
!-----------------------------
      vctr  =0
      vctrp =0
      call input_vect(1,
     &     ifli,ifll,ifle,ical_vect,nthrds,iLBF_P,ical_prt,
     &     NCV,NCVIN,NALLCV,NCVFAC,NMAT,NSSFBC,MXMAT,MAXIE,
     &     MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,
     &     vctr,vctrp,NEIGHBOR,DNORM,SFAREA,SFCENT,
     &     MXCV_V,MXBND_V,MXCV_VP,MXBND_VP,
     &     msk,MXALLCV,LVEDGE,MXCVFAC,1,
     &     MEP,MXALLCV_P,MXCVFAC_P,
     &     ierror)
!----------------------------
! --- 
!----------------------------
      if(ical_vect) then
! 2D
        if(intgvv==euleri.or.intgty==euleri.or.intgke==euleri) then
          IVDIM=8
        else
          IVDIM=4
        endif
!
        deallocate(W1K6 )
        allocate(IW2K2(MXCV,MAXIE)    ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate IW2K2 error in METR_FFLOW')
        allocate(W2K1(MXCV,0:MAXIE)   ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate W2K1 error in METR_FFLOW')
        allocate(W2K2(MXCV+N2,IVDIM)  ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate W2K2 error in METR_FFLOW')
        allocate(iwork7(0:MXCV)       ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate iwork7 error in METR_FFLOW')
        allocate(iwork5(MXCV)         ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate iwork5 error in METR_FFLOW')
        allocate(iwork6(0:MXCV)       ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate iwork6 error in METR_FFLOW')
        allocate(iwork4(MXCV)         ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate iwork4 error in METR_FFLOW')
        allocate(d1work2(MXCV+N2)     ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate d1work2 error in METR_FFLOW')
        allocate(W1K6(0:MXCV)         ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate W1K6 error in METR_FFLOW')
        if(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
          allocate(d2vect(MXCVFAC,1)    ,stat=ierror)
          if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate d2vect error in METR_FFLOW')
          allocate(tmpsld(MXSSFBC_SLD,9),stat=ierror)
          if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate tmpsld error in METR_FFLOW')
        else
          allocate(d2vect(MXCVFAC,1)    ,stat=ierror)
          if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate d2vect error in METR_FFLOW')
        endif
        allocate(IW2K_VECT(MXCV,MAXIE),stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate IW2K_VECT error in METR_FFLOW')
        allocate(DW2K_VECT(MXCV,MAXIE),stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate DW2K_VECT error in METR_FFLOW')
        allocate(iwork1(MXCVFAC)      ,stat=ierror)
        if(ierror/=0) 
     &  call FFRABORT(1,'ERR: allocate iwork1 error in METR_FFLOW')
!
        IW2K2(:,:)=0
        W2K1(:,:)=0.d0
        W2K2(:,:)=0.d0
        iwork7(:)=0
        iwork5(:)=0
        iwork6(:)=0
        iwork4(:)=0
        iwork1(:)=0
        d1work2(:)=0
        d2vect(:,:)=0.d0
        IW2K_VECT(:,:)=0
        DW2K_VECT(:,:)=0.d0
        W1K6(:)=0
!
        if(ierror/=0) call 
     &  FFRABORT(1,'ERR:allocate ical_vect in METR_FFLOW')
      endif
!
!
      IF(iLBF_P==1.or.iLBF_P==2.or.iLBF_P==5) THEN
        allocate(bdyfrc(MXCV_B,3,2),stat=ierror)
        if(ierror/=0) call 
     &  FFRABORT(1,'ERR:allocate bdyfrc in METR_FFLOW')
        bdyfrc(:,:,:)=0.d0
      else
        allocate(bdyfrc(1,3,2),stat=ierror)
      ENDIF
!
      allocate(WDRBND(MXCV_WDR,N_inj),stat=ierror)
      if(ierror/=0) call 
     &  FFRABORT(1,'ERR:allocate rst,vflc,yplusf in METR_FFLOW')
!-----------------------------
! --- flamelet model
!-----------------------------
      allocate(vflc(MXCV_F),yplusf(MXCV_F),stat=ierror)
      allocate(rva_s(MXCVFAC_F),stat=ierror)
      if(ierror/=0) call 
     &  FFRABORT(1,'ERR:allocate rst,vflc,yplusf in METR_FFLOW')
!
! --- 
!
      if(idifflm/=ORG) then
          if(comb_mdl==1) then
            akshg(ixi)=x1max
            akslw(ixi)=x1min
            akshg(ixi_2)=x2max
            akslw(ixi_2)=x2min
            akshg(ixi_X)=x3max
            akslw(ixi_X)=x3min
          elseif(comb_mdl==2) then
            akshg(ixi)=x1max
            akslw(ixi)=x1min
            akshg(ixi_2)=x2max
            akslw(ixi_2)=x2min
            akshg(ixi_C)=x3max
            akslw(ixi_C)=x3min
          endif
      endif
!
!
      if(ista_momery==1) then
        n_ave=0
        if(iuvw_ave_rms_re==1) then
          n_ave=9
        endif
        if(ip_ave==1) then
          n_ave=n_ave+1
        endif
        if(imu_ave==1) then
          n_ave=n_ave+1
        endif
        if(ip_rms.eq.1) then
          n_ave=n_ave+1
        endif
        if(imu_rms.eq.1) then
          n_ave=n_ave+1
        endif
        if(imaxmin==1) then
          n_ave=n_ave+8
        endif
        if(it_ave==1) then
          n_ave=n_ave+1
        endif
        if(it_rms.eq.1) then
          n_ave=n_ave+1
        endif
        ! Debug Tagami
        if(iuvwt_re == 1 ) then
           n_ave=n_ave+3
        endif
        ! Debug end
        do icom=1,ncomp
           if(icomp_ave(icom).eq.1) then
             n_ave=n_ave+1
           endif
           if(icomp_rms(icom).eq.1) then
             n_ave=n_ave+1
           endif
        enddo
        do i = 1, ncomp
           if( iuvwc_re(i) == 1 ) then
              n_ave=n_ave+3
           endif
        enddo
        do i=1,nrans !nrnsx
          if(irans_ave(i)==1) then
            n_ave=n_ave+1
          endif
          if(irans_rms(i).eq.1) then
            n_ave=n_ave+1
          endif
          if(iuvw_rans_re(i).eq.1) then
            n_ave=n_ave+1
          endif
        enddo
        if(ical_prt==1) then
          n_ave=n_ave+1
        endif
        allocate(STA_SAV(NCV,n_ave),stat=ierror)
        if(ierror/=0) call FFRABORT(1,'ERR: allocating STA_SAV error')
      endif

!
      if(ical_FC==PEFC.and.iBV==1) then
        allocate(J_SRC(MXALLCVP,2),stat=ierror)
      else
        allocate(J_SRC(1,2),stat=ierror)
      endif
!
      if(ical_thmoDIFF==1) then
        allocate(TH_DIFF(MXALLCV,ncomp),stat=ierror)
      else
        allocate(TH_DIFF(1,ncomp),stat=ierror)
      endif
!-------------------
! --- Reaction flag
!-------------------
      allocate(lreac(nneq),stat=ierror)
!      allocate(lreacP(nneqP),stat=ierror)
!
      if(u_func(2)==1.or.u_func(3)==1) then
        allocate (multiR(MXCV),stat=ierror)
        multiR(:)=0
        drag_coef(:)=0.d0
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT.lt.0) cycle
        do ICVL=ICVS,ICVE
        dum1=0.d0
        xyz(1)=CVCENT(1,ICVL)
        xyz(2)=CVCENT(2,ICVL)
        xyz(3)=CVCENT(3,ICVL)
        idum1=1
        call USER_UserDefine_region(IMAT,ICVL,xyz,idum1,dum1)
        if(u_func(3)==1) then
          if(idum1>=21.and.idum1<=50) then
            drag_coef(idum1)=dum1
          endif
        endif
        if (idum1==3) u_func(2)=2
        multiR(ICVL)=int(idum1)
        enddo
        enddo
      endif
!
      if(icaltb==dles.or.u_func(2)==2) then
        allocate(lij(MXALLCV,6)   ,stat=ierror)
        allocate(mij(MXALLCV,6)   ,stat=ierror)
      endif


!-----------------------------
! --- 
!-----------------------------
      IF(ical_mvmsh/=0) THEN
        if(NPE>1) then
!          call read_move_ini_grid_HPC()
           call FFRABORT(1,'ERR: NOT support HPC for moving mesh')
        else
          call read_move_ini_grid(cord,lacell,lvcell,
     &         lvface,lbface,lcface,lfcell,LEFACE,listbc)
!
          if(ipiston==1) then
            call piston(ical_mvmsh,iter,time,
     &        NVRTX,MXVRTX_m,MXFACE_m,MXCELL_m,
     &        movflg,cord,lbface,lacell,lvface,lvcell,lcface,
     &        NBOUND,NBFS,SFBOUN,NFBOUN,IFFACE,1)
          else
            call user_moving_mesh(ical_mvmsh,iter,time,
     &        NVRTX,MXVRTX_m,MXFACE_m,MXCELL_m,
     &        movflg,cord,lbface,lacell,lvface,lvcell,lcface,
     &        NBOUND,NBFS,SFBOUN,NFBOUN,IFFACE,1)
          endif
        endif
!
      ENDIF
!
! --- vof
!
      allocate(t_dir(MXALLCV_VOF,3),stat=ierror)
      allocate(ANGBND(MXSSFBC_VOF),stat=ierror)
      if(ierror/=0) call 
     &  FFRABORT(1,'ERR:allocate ical_vof in METR_FFLOW')
!------------------
! --- shuuter
!------------------
      shutt=0
      do nb=1,nbcnd
      kd=kdbcnd(0,nb)
      pkdwall(nb)=kd==kxnone.or.
     &            kd==kdintr.or.
     &            kd==kdfire.or.
     &            kd==kdsymm.or.
     &            kd==kdcvd

      if(kd==kdshutr) then
        shutt=1
      endif
      enddo
!
      if(shutt==1) then
         MXSSFBC_SHT=mxssfbc
      else
         MXSSFBC_SHT=mxssfbc   !1
      endif
      allocate(SHUTFL(MXSSFBC_SHT),stat=ierror)

      if(ierror/=0) call 
     &  FFRABORT(1,'ERR:allocate SHUTFL in METR_FFLOW')
      SHUTFL(:)=0
!------------------------------
! --- radition
!------------------------------
      call list_fluid
     &  (LVEDGE,MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,LBC_SSF,
     &   Pstart,KDBF)
!------------------------------
! --- 
!------------------------------  
!     
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcimin(lclsd(IIMAT))
!        CALL hpcimin(iclosd(IIMAT))
!        lclsd(0)=10
!        iclosd(0)=10
!        enddo
        DO IIMAT=1,NMAT 
        CALL hpcimin(lclsd(IIMAT))
        CALL hpcimin(iclosd(IIMAT))
        enddo
      endif

!
!
!
      matsld(:,:)=0
!      sldbcmat(:,:)=0

      do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdovst.or.kd==kdsld.or.
     &     kd==kdbuff.OR.kd==kdshutr) then
          IIMAT1=MAT_BCIDX(nb,1)
          IIMAT2=MAT_BCIDX(nb,2)
          if(IIMAT1==0) cycle
          if(IIMAT2==0) cycle
          matsld(IIMAT1,1)=1
          matsld(IIMAT2,1)=1
        endif
      enddo
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcimax(matsld(IIMAT))
!        matsld(0)=0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcimax(matsld(IIMAT,1))
        enddo
      endif
!
      idum1=100000
      idum2=100000
      DO IIMAT=1,NMAT
      if(matsld(IIMAT,1)==1) then
        idum1=min(idum1,lclsd(IIMAT))
      endif
      if(matsld(IIMAT,1)==1) then
        idum2=min(idum2,iclosd(IIMAT))
      endif
      enddo
!
      DO IIMAT=1,NMAT
      if(matsld(IIMAT,1)==1) then
        lclsd(IIMAT)=idum1
      endif
      if(matsld(IIMAT,1)==1) then
        iclosd(IIMAT)=idum2
      endif
      enddo
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcimin(lclsd(IIMAT))
!        CALL hpcimin(iclosd(IIMAT))
!        enddo
        DO IIMAT=1,NMAT 
        CALL hpcimin(lclsd(IIMAT))
        CALL hpcimin(iclosd(IIMAT))
        enddo
      endif
!
      call mat_flag(0,Pstart,mat_cal)
!-----------------------------
!-< 3. Calculate metrics >-
!-----------------------------
      call metric_admin(1,iter,time,
     &  LVEDGE,LVRTCV,LFUTAU,
     &  LBC_SSF,LCYCSF,FRSTCV,DISINL,
     &  SFAREA,SFCENT,CVVOLM,CVVOL0,CVCENT,wiface,DISALL,xta,
     &  MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LCYCOLD,wifsld,OPPANG,
     &  movflg,cord,lacell,lvcell,
     &  lvface,lbface,lcface,lfcell,LEFACE,listbc,locmsh,
     6  area,volume,gface,gcell,
     &  delt_tmp,ierror)
      if(ierror.ne.0) goto 9999
!
!      deallocate(sldmat)
      return
!
 9999 continue
      write(ifle,*) '(METR_FFLOW)'
      ierror=1
      return
      end subroutine METR_FFLOW
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_admin(icall,iter,time,
     &  LVEDGE,LVRTCV,LFUTAU,
     &  LBC_SSF,LCYCSF,FRSTCV,DISINL,
     &  SFAREA,SFCENT,CVVOLM,CVVOL0,CVCENT,wiface,DISALL,xta,
     &  MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LCYCOLD,wifsld,OPPANG,
     &  movflg,cord,lacell,lvcell,
     &  lvface,lbface,lcface,lfcell,LEFACE,listbc,locmsh,
     6  area,volume,gface,gcell,
     &  deltt,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension
      use module_hpcutil
      use module_io,      only : ifle,ifll
      use module_movegrid,only : ngrid,dtmin,grdtim=>time
      use module_boundary,only : kdbcnd,kxnone,lovstbc,lsldbc
      use module_boundary,only : SFBOUN,NBOUND,
     &                           IFFACE,NFBOUN,NBFS,nbcnd,LBC_INDEX
      use module_model,   only : icaltb,sles,dles,lles,ke,RSM
      use module_material,only : cvdist,ical_sld
      use module_model,   only : ical_mvmsh
      use module_metrix,only   : SFCENT_OLD =>d2work2
      use module_model,only    : ical_vect
      use module_movegrid,only : ipiston
!
      implicit none 
!
! --- [dummy arguments]
!
      integer,intent(in)      :: iter,icall
      real*8 ,intent(in)      :: time
      integer,intent(inout)   :: LVEDGE(2,  MXCVFAC)
      integer,intent(in)      :: LVRTCV(    MXALLCV)
      integer,intent(in)      :: LFUTAU(    MXCV)
      logical,INTENT(IN)      :: mat_cal  ( 0:MXMAT)
      integer,intent(in)      :: MAT_NO(    0:MXMAT)
      integer,intent(in)      :: MAT_CV(    MXALLCV)
      integer,intent(in)      :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)      :: MAT_DCIDX( 0:MXMAT)
      integer,intent(in)      :: MAT_CFIDX( 0:MXMAT)
      integer,intent(in)      :: LBC_SSF(   mxssfbc)
      integer,intent(inout)   :: LCYCSF(    mxssfbc)
!
      real*8 ,intent(out)     :: FRSTCV(    mxssfbc)
      real*8 ,intent(out)     :: DISINL(    MXSSFBC)

      real*8 ,intent(inout)   :: SFAREA(4,  MXCVFAC)
      real*8 ,intent(inout)   :: SFCENT(3,  MXCVFAC)
      real*8 ,intent(inout)   :: CVVOLM(    MXALLCV)
      real*8 ,intent(inout)   :: CVVOL0(    MXALLCV)
      real*8 ,intent(inout)   :: CVCENT(3,  MXALLCV)
      real*8 ,intent(out)     :: wiface(    MXCVFAC)
      real*8 ,intent(out)     :: DISALL(    MXCV)
      real*8 ,intent(INout)   :: xta   (    MXCVFAC)
      real*8 ,intent(in)      :: deltt
      integer,intent(out)     :: ierror
!
      integer,intent(inout)   :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(inout)   :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(inout)   :: OPPANG     (MXSSFBC_SLD)
!
      integer,intent(inout)   :: movflg     (  MXVRTX_m)
      real*8 ,intent(inout)   :: cord       (3,MXVRTX_m)
      integer,intent(inout)   :: lacell     (  MXCELL_m)
      integer,intent(inout)   :: lvcell     (8,MXCELL_m)
      integer,intent(inout)   :: lvface     (4,MXFACE_m)      
      integer,intent(inout)   :: lbface     (2,MXFACE_m)
      integer,intent(inout)   :: lcface     (2,MXFACE_m)
      integer,intent(inout)   :: lfcell     (7,MXCELL_m)
      integer,intent(inout)   :: LEFACE     (5,MXFACE_m)
      integer,intent(inout)   :: listbc     (4,MXSSFBC_old_m)
      real*8 ,intent(inout)   :: area       (4,MXFACE_m)
      real*8 ,intent(inout)   :: volume     (  MXCELL_m)
      real*8 ,intent(inout)   :: gface      (3,MXFACE_m)
      real*8 ,intent(inout)   :: gcell      (3,MXCELL_m)
      integer,intent(inout)   :: locmsh     (   MXSSFBC)
!
! --- [local entities] 
!
      integer :: nb,kd
      integer :: IBFS,IBFE,IBFL,ICV
      integer :: IMAT,IIMAT,ICFS,ICFE,ICFL,ICVS,ICVE,ICVL
      integer :: i,n,icyc,lgd,lgdx,ierr1=0,ierr2
      integer :: lgrid=-1
      real*8  :: dtgrd=0.d0
      real*8  :: tim1,tcyc,tgd,dt1,dt2,dum1,yy,vol
      real*8 ,parameter :: SML=1.D-25,ZERO=0.D0,GREAT=1.D25 
      logical,save :: static=.false.
!
      ierror=0
!
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
        if(lsldbc) then
          if(ical_vect) then
            call list_sliding_V(icall,iter,time,deltt,
     &      MAT_NO,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &      LBC_SSF,SFCENT,LCYCOLD,wifsld,OPPANG,
     &      LCYCSF,LVEDGE,CVVOLM,CVCENT,SFAREA)
          else
            call list_sliding(icall,iter,time,deltt,
     &      MAT_NO,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &      LBC_SSF,SFCENT,LCYCOLD,wifsld,OPPANG,
     &      LCYCSF,LVEDGE,CVVOLM,CVCENT,SFAREA)
          endif
        endif
        if(lovstbc) then 
          call list_overset(icall,iter,time,deltt,
     &      MAT_NO,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &      LBC_SSF,SFCENT,
!     &      LCYCOLD,wifsld,OPPANG,
     &      LCYCSF,LVEDGE,CVVOLM,CVCENT,SFAREA)
        endif
      endif
!
      if(static) then
        return
      endif
!
! --- 
!
      do nb=1,nbcnd
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      do IBFL=IBFS,IBFE
      ICFL=LBC_SSF(IBFL)
      ICV=LVEDGE(1,ICFL)
      if(ICV<=NCVIN) then
        locmsh(IBFL)=1
      endif
      enddo
      enddo
!
!
      if(.not.(ngrid>=2.or.ical_mvmsh/=0)) static=.true.
!
      if(ical_mvmsh/=0) then
!--------------------------------------------------------------
! --- removing & add coordinate for next time mesh moving 
!--------------------------------------------------------------
!-------------------------------------------
! --- coordinate for next time mesh moving  
!-------------------------------------------
        if(ipiston==1) then
          call piston(ical_mvmsh,iter,time,
     &        NVRTX,MXVRTX_m,MXFACE_m,MXCELL_m,
     &        movflg,cord,lbface,lacell,lvface,lvcell,lcface,
     &        NBOUND,NBFS,SFBOUN,NFBOUN,IFFACE,2)
          else
            call user_moving_mesh(ical_mvmsh,iter,time,
     &        NVRTX,MXVRTX_m,MXFACE_m,MXCELL_m,
     &        movflg,cord,lbface,lacell,lvface,lvcell,lcface,
     &        NBOUND,NBFS,SFBOUN,NFBOUN,IFFACE,2)
          endif
!
        if(ical_mvmsh==1.or.
     &     ical_mvmsh==2.or.
     &     ical_mvmsh==4.or.
     &     ical_mvmsh==5) then
!--------------------------------------------------------------
! --- calculate new time-step's mesh area, volume & center 
!--------------------------------------------------------------
          call metric_3ds
     &   (icall,cord,lacell,lvcell,
     &    lvface,lbface,lcface,lfcell,
     &    area,volume,gface,gcell,ierr1)
!-------------------------------------------
! --- CV area, volume & center 
!-------------------------------------------
          allocate(SFCENT_OLD(3,MXCVFAC),stat=ierr1)
          if(ierr1/=0) call FFRABORT(1,'ERR: allocate SFCENT_OLD')
          do i=1,3
          SFCENT_OLD(i,:)=SFCENT(i,:)
          enddo
          call metric_CV 
     &   (LBC_SSF,lcface,lbface,lacell,
     &    cord,area,volume,gface,gcell,
     &    LVEDGE,LEFACE,LVRTCV,listbc,MAT_CV,
     &    SFAREA,SFCENT,CVVOLM,CVCENT,
     &    ierr2)
!-------------------------------------------
! --- 
!-------------------------------------------
          call dc_metric
     &   (LVEDGE,LBC_SSF,LCYCSF,SFCENT,SFAREA,CVCENT,CVVOLM,2)
!-------------------------------------------
! --- 
!-------------------------------------------
          dum1=deltt
          call metric_speed(iter,MAT_CFIDX,SFCENT_OLD,SFCENT,SFAREA,xta,dum1)
          deallocate(SFCENT_OLD)
!
        endif
!
      else
!
!      if(deltt.gt.dtmin) goto 9001
!
!      if(lgrid.lt.0) then
!        xta=0.d0
!      endif
!
!-< 1. Search present position >-
!
!                                         icyc=0
!      if( grdtim(ngrid+1).gt.grdtim(1) ) icyc=1
!
!      tgd=time
!      if(icyc.gt.0) then
!        tim1=grdtim(1)
!        tcyc=grdtim(ngrid+1)-tim1
!        tgd=tgd-tim1
!        dum1=tgd/tcyc
!        tgd=tgd-dble(int(dum1))*tcyc
!        if( tgd.lt.0.d0 ) tgd=tgd+tcyc
!        tgd=tgd+tim1
!      endif
!      lgdx=icyc
!      do 100 n=icyc+1,ngrid
!      if(grdtim(n).lt.tgd) lgdx=n
!  100 continue
!      lgd=max(1,min(ngrid+icyc-1,lgdx))
!
!      dt1=0.d0
!      dt2=1.d0
!      if( lgdx.eq.lgrid ) goto 1001
!
!
!-< 2. Calculate metrics >-
!
!--< 2.3 calculate metrics of grid speed >--
!
!      if(lgd.eq.lgdx) then
!        dum1=grdtim(lgd+1)-grdtim(lgd)
!      else
!      endif
!
!      dt1=max(0.d0,min(1.d0,dtgrd/deltt))
!      dt2=1.d0-dt1
!
! 1001 continue
!      lgrid=lgdx
!      dtgrd=grdtim(lgdx+1)-tgd
!
!-< 5. Set metrics at present time >-
!
!--< 5.2 interpolation in time >--
!
!      dt2=max(0.d0,min(1.d0,
!     &    (tgd-grdtim(lgd))/(grdtim(lgd+1)-grdtim(lgd))))
!      dt1=1.d0-dt2
      endif
!
! --- HPC for exter vertex
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,CVVOLM)
        ierr1=0
        do i=1,3
!        CVVOL0(:)=CVCENT(i,:)
        CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,CVCENT(i,:))
!        CVCENT(i,:)=CVVOL0(:)
        enddo
!        CVVOL0=CVVOLM
      ENDIF
!
!--< 5.4 set interpolation factor >--
!
       call metric_intrpw(1,MAT_CFIDX,CVVOLM,
     &  LVEDGE,SFCENT,CVCENT,SFAREA,wiface)
!
!-< 6. Deallocate arrays in case of static grid >-
!
!-< 9. wall distance
!
      
      call wall_distance
     &    (LVEDGE,LBC_SSF,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,
     &     SFAREA,SFCENT,CVCENT,CVVOLM,
     &     LFUTAU,FRSTCV,DISALL)

!
! --- list Inlet BC
!
      call list_bcin(LBC_SSF,SFAREA,SFCENT,DISINL)
!
      return
!
 9001 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'time increment is bigger than minimum ',
     &           'time interval of grid tables'
      write(ifle,*) 'time increment = ',deltt
      write(ifle,*) 'minimum time interval = ',dtmin
      goto 9999
 9999 continue
      write(ifle,*) '(metric_admin)'
      ierror=1
!
      end subroutine metric_admin
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_intrpw(IMODE,MAT_CFIDX,CVVOLM,
     &  LVEDGE,SFCENT,CVCENT,SFAREA,wiface)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. Calculate weight rate for calculation of interpolate face value
!
! --- [module arguments]
!
      USE module_dimension
      use module_hpcutil,  only : NPE,my_rank
      use module_model,only     : ical_vect,nthrds
      use module_vector,only    : ICVS_V,ICVE_V,
     &                            ICFS_V,ICFE_V,
     &                            ICVSIN_V,ICVEIN_V,
     &                            IDCS_V,IDCE_V,index_c,index_f
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: IMODE
      integer,intent(in)  :: MAT_CFIDX(0:MXMAT)
      integer,intent(in)  :: LVEDGE(2, MXCVFAC)
      real*8 ,intent(in)  :: SFCENT(3, MXCVFAC)
      real*8 ,intent(in)  :: CVCENT(3, MXALLCV)
      real*8 ,intent(in)  :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)  :: CVVOLM(   MXALLCV)
      real*8 ,intent(out) :: wiface(   MXCVFAC)
!
! --- [local entities]
!
      real*8  :: dr1,dr2
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      integer :: IBFS,IBFE,IBFL
      integer :: IMAT,IIMAT,ICFS,ICFE,ICFL,ICVS,ICVE,ICVL
      integer :: ICVLA,ICVLB
!
! --- vertex center
!
      if(ical_vect) then
!!$omp do private(ICFL,ICVLA,ICVLB)
!CIDR NODEP        
        DO ICFL=1,NCVFAC
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        wiface(ICFL) =dsqrt((CVCENT(1,ICVLB)-SFCENT(1,ICFL))**2
     &                     +(CVCENT(2,ICVLB)-SFCENT(2,ICFL))**2
     &                     +(CVCENT(3,ICVLB)-SFCENT(3,ICFL))**2)
     &              /(dsqrt((CVCENT(1,ICVLA)-SFCENT(1,ICFL))**2
     &                     +(CVCENT(2,ICVLA)-SFCENT(2,ICFL))**2
     &                     +(CVCENT(3,ICVLA)-SFCENT(3,ICFL))**2)
     &               +dsqrt((CVCENT(1,ICVLB)-SFCENT(1,ICFL))**2
     &                     +(CVCENT(2,ICVLB)-SFCENT(2,ICFL))**2
     &                     +(CVCENT(3,ICVLB)-SFCENT(3,ICFL))**2)+SML)
        enddo
!!$omp end do
      else!IF(IMODE==1) then
        do 100 IIMAT=1,NMAT
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dr1=dsqrt( (CVCENT(1,ICVLA)-SFCENT(1,ICFL))**2
     &            +(CVCENT(2,ICVLA)-SFCENT(2,ICFL))**2
     &            +(CVCENT(3,ICVLA)-SFCENT(3,ICFL))**2)
        dr2=dsqrt( (CVCENT(1,ICVLB)-SFCENT(1,ICFL))**2
     &            +(CVCENT(2,ICVLB)-SFCENT(2,ICFL))**2
     &            +(CVCENT(3,ICVLB)-SFCENT(3,ICFL))**2)
        wiface(ICFL)=dr2/(dr1+dr2+SML)
        enddo
  100   enddo
!      elseIF(IMODE==2) then
!        do 200 IIMAT=1,NMAT
!        ICFS=MAT_CFIDX(IIMAT-1)+1
!        ICFE=MAT_CFIDX(IIMAT)
!        do ICFL=ICFS,ICFE
!        ICVLA=LVEDGE(1,ICFL)
!        ICVLB=LVEDGE(2,ICFL)
!        dr1=CVVOLM(ICVLA)
!        dr2=CVVOLM(ICVLB)
!        wiface(ICFL)=dr2/(dr1+dr2+SML)
!        enddo
!  200   enddo
      endif
!
      end subroutine metric_intrpw
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_fluid
     &  (LVEDGE,MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,LBC_SSF,
     &   Pstart,KDBF)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil,only : NPE,my_rank
      use module_io,only      : ifle,ifll
      use module_model,only   : comp,idrdp,incomp,mach0,iHPC_close
      use module_material,only: deficld,PFstart,PSstart
      use module_metrix ,only : iclos=>msk
      use module_chemreac,only: ireq,nneq,ovalfr
      use module_time,only    : steady,MHD_steady
      use module_vof,only     : ical_vof
      use module_model,only   : ical_MHD
!
! 2. Set flag for closed domain
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)    :: LBC_SSF(   mxssfbc)
      integer,intent(in)    :: MAT_NO(    0:MXMAT)
      integer,intent(in)    :: MAT_CV(    MXALLCV)
      integer,intent(in)    :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)    :: MAT_DCIDX( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX( 0:MXMAT)
      integer,intent(in)    :: LVEDGE(2,  MXCVFAC)
      integer,intent(out)   :: KDBF     ( MXCVFAC)
      integer,intent(out)   :: Pstart(      MXMAT)
!
! --- [local entities]
!
      integer :: jclos (  0:100 )
!      integer :: MATC_SLD(100),MATC_FLD(100)
      integer :: nb,kd,ierror
      integer :: IBFS,IBFE,IBFL
      integer :: IMAT,IIMAT,ICFS,ICFE,ICFL,ICVS,ICVE,ICVL
      integer :: ICVLA,ICVLB,ioval
      integer :: ierr1,ierr,nfludx,NSLIDX,j
!
      ierror=0
      ierr=0
      ierr1=0
!
!-< 3. Set flag for closed domain >-
!
      call bc_kdbp_cell(LBC_SSF,KDBF)
!
      iclos(:)=0
!
!      NFLUDX=0
!      NSLIDX=0
!      MATC_FLD=0
!      MATC_SLD=0
!      do 300 IIMAT=1,NMAT 
!      IMAT=MAT_NO(IIMAT)
!      if(IMAT.gt.0) then
!        NFLUDX=NFLUDX+1
!        MATC_FLD(NFLUDX)=IIMAT
!!        Pstart(IIMAT)=PFstart(IMAT)
!      elseif(IMAT.lt.0) then
!        NSLIDX=NSLIDX+1
!        MATC_SLD(NSLIDX)=IIMAT
!!        Pstart(IIMAT)=PSstart(-IMAT)
!      endif
! 300  continue
!
      jclos(:)=0
!
!      if(NFLUDX+NSLIDX.ne.NMAT) then
!        write(*,*) '### error: NFLUDX-NSLIDX is not',NFLUDX,NSLIDX,NMAT
!        CALL FFRABORT(1,'list_fluid')
!      endif
!
!      IF(ical_vect) then
        
!      else
      do 310 IIMAT=1,NMAT
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(KDBF(ICFL).eq.1) then
! --- Outlet BC
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          iclos(ICVLA)=1
          iclos(ICVLB)=1
        elseif(KDBF(ICFL).eq.2) then
! --- pressure or stag-pressure  !NOT finished 
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          iclos(ICVLA)=2
          iclos(ICVLB)=2
        elseif(KDBF(ICFL).eq.4) then
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          iclos(ICVLA)=4
          iclos(ICVLB)=4
        endif
        enddo
  310 enddo
!      endif
!
      do 320 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      if(IMAT.gt.0) then
        if(iclos(ICVL).eq.1) then
          jclos(IMAT)=1
        elseif(iclos(ICVL).eq.2) then
          jclos(IMAT)=2
        elseif(iclos(ICVL).eq.3) then
          jclos(IMAT)=3
        elseif(iclos(ICVL).eq.4) then
          jclos(IMAT)=4
        endif
      endif
      enddo
  320 continue
!
      ioval=0
      if (nneq.gt.0) then
        do j=1,nneq
        if(ireq(j).eq.ovalfr) then
          ioval=1       !lclsd(IIMAT)=0
        endif
        enddo
      endif
!
      call deficld
     &  (NMAT,my_rank,ifle,NFLUDX,MXMAT,MAT_NO,
     &   idrdp,incomp,mach0,comp,jclos,Pstart,
     &   ioval,ical_vof,steady,iHPC_close,MHD_steady,
     &   ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      return
!
 9999 continue
      write(ifle,*) '(list_fluid)'
      ierror=1
!
      end subroutine list_fluid
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_bcin(LBC_SSF,SFAREA,SFCENT,DISINL)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension
      use module_hpcutil,only  : NPE,my_rank
      use module_io,only       : ifle,ifll
      use module_boundary,only : kdbcnd,kdilet,kvlglw,kvnslp,kxnone,
     &                           kdolet,kdprdc,kdsymm,kdtchi,distrb,
     &                           LBC_INDEX,nbcnd,kvrogf,kvmodl,kvEWF,
     &                           stagare,kdstag,kdpres,kvsataW
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)    :: LBC_SSF(   MXSSFBC)
      real*8, intent(out)   :: DISINL(    MXSSFBC)
      real*8 ,intent(in)    :: SFAREA(4,  MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3,  MXCVFAC)
!
! --- [local entities]
!
      real*8  :: DISMIN
      real*8 ,parameter :: SML=1.D-25,ZERO=0.D0,GREAT=1.D25
!      integer :: NBCINLo,NBCCYCo,NBCOUTo,NBCSYMo,NBCTCIo
      integer :: IWL,ICFmin
      real*8  :: dx,dy,dz,dis,rdum1
      integer :: IBFS,IBFE,IBFL,IBFS1,IBFE1,IBFL1,ICFL,ICFL1
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,nb1
!
!
      DISINL=0.D0
      stagare=0.D0
!
!
! --- Distance from INLET 2D cell to wall's edge
!
!      if(distrb.eq.0) goto 3000
!
      do 2000 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if(kd.eq.kdilet) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do 2300 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        DISMIN=GREAT
        ICFmin=-1
        do 2100  nb1=1,nbcnd
        kdv=kdbcnd(1,nb1)
        if(kdv==kvlglw.or.kdv.eq.kvnslp.or.
     &     kdv==kvrogf.or.kdv==kvmodl.or.kdv==kvEWF.or.kdv==kvsataW)
     &  then
          IBFS1=LBC_INDEX(nb1-1)+1
          IBFE1=LBC_INDEX(nb1)
          do 2200 IBFL1=IBFS1,IBFE1
          ICFL1=LBC_SSF(IBFL1)
          dx=SFCENT(1,ICFL1)-SFCENT(1,ICFL)
          dy=SFCENT(2,ICFL1)-SFCENT(2,ICFL)
          dz=SFCENT(3,ICFL1)-SFCENT(3,ICFL)
          dis=dsqrt(dx*dx+dy*dy+dz*dz)
          if(dis.lt.dismin) then
            dismin=dis
            ICFmin=ICFL1
          end if
 2200     continue
        endif
 2100   continue
        if(ICFmin>0) then
        DISINL(IBFL)=
     &     abs(SFAREA(1,ICFmin)*(SFCENT(1,ICFL)-SFCENT(1,ICFmin))
     &       + SFAREA(2,ICFmin)*(SFCENT(2,ICFL)-SFCENT(2,ICFmin))
     &       + SFAREA(3,ICFmin)*(SFCENT(3,ICFL)-SFCENT(3,ICFmin)))
        end if
 2300   continue
      endif
!
!      IBFS=LBC_INDEX(nb-1)+1
!      IBFE=LBC_INDEX(nb)
!      do IBFL=IBFS,IBFE
!      ICFL=LBC_SSF(IBFL)
!      stagare(nb)=stagare(nb)+SFAREA(4,ICFL)
!      enddo
 2000 continue
!
 3000 continue
!
!
!
      do nb=1,nbcnd
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      do IBFL=IBFS,IBFE
      ICFL=LBC_SSF(IBFL)
      stagare(nb)=stagare(nb)+SFAREA(4,ICFL)
      enddo
      enddo
!      
      return
!
      end subroutine list_bcin
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine wall_distance
     &    (LVEDGE,LBC_SSF,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,
     &     SFAREA,SFCENT,CVCENT,CVVOLM,
     &     LFUTAU,FRSTCV,DISALL)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
!      use module_hpcutil,only  : NPE,my_rank
      use module_io,only       : ifle,ifll,lenfnm,getfil
      use module_boundary,only : kdbcnd,kvnslp,kvfslp,kxnone,nbcnd,
     &                           LBC_INDEX,kvlglw,kvrogf,kvmodl,kvEWF,
     &                           kdshutr,kvsataW
      use module_model,   only : icaltb,ke,sles,dles,lles,RSM,SDES,
     &                           ical_mvmsh
      use module_material,only : cvdist
      
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: LVEDGE(2,  MXCVFAC)
      integer,intent(in)  :: LBC_SSF   (MXSSFBC)
      integer,intent(in)  :: MAT_CV(    MXALLCV)
      integer,intent(in)  :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)  :: MAT_DCIDX( 0:MXMAT)
      integer,intent(in)  :: MAT_NO(    0:MXMAT)
      real*8 ,intent(in)  :: SFAREA(4,  MXCVFAC)
      real*8 ,intent(in)  :: SFCENT(3,  MXCVFAC)
      real*8 ,intent(in)  :: CVCENT(3,  MXALLCV)
      real*8 ,intent(in)  :: CVVOLM(    MXALLCV)
      integer,intent(out) :: LFUTAU(    MXCV)
      real*8 ,intent(out) :: FRSTCV(    mxssfbc)
      real*8 ,intent(out) :: DISALL(    MXCV)
!
! --- [local entities]
!
      integer :: IBFS,IBFE,IBFL
      integer :: IMAT,IIMAT,ICFS,ICFE,ICFL,ICVS,ICVE,ICVL,NCVX
      integer :: ICVLA,ICVLB
      integer :: NB,KD,kdv
      integer :: NTAU,ICFMIN
      integer :: ios=0,ifl
      logical :: wallfile=.false.
      real*8  :: dx,dy,dz,dismin,dis1,dum1,yy,vol,yy1,yy2,yy3
      real*8 ,parameter :: SML=1.D-25,ZERO=0.D0,GREAT=1.D25
      character(len=1) :: recal
      character(lenfnm),save :: fnam
!
!--< 1. Measure Distance for first CV .
!
!
! --- Distance for First CV
!
!      FRSTCV=1.D0
!      do 600 nb=1,nbcnd
!      IBFS=LBC_INDEX(nb-1)+1
!      IBFE=LBC_INDEX(nb)
!      do IBFL=IBFS,IBFE
!      ICFL=LBC_SSF(IBFL)
!      ICVL=LVEDGE(1,ICFL)
!!1)
!      yy2=abs(SFAREA(1,ICFL)*(CVCENT(1,ICVL)-SFCENT(1,ICFL))
!     &       +SFAREA(2,ICFL)*(CVCENT(2,ICVL)-SFCENT(2,ICFL))
!     &       +SFAREA(3,ICFL)*(CVCENT(3,ICVL)-SFCENT(3,ICFL)))
!!2)
!      yy=CVVOLM(ICVL)**(1.d0/3.d0)
!!3)    
!!      yy1=dsqrt((CVCENT(1,ICVL)-SFCENT(1,ICFL))**2
!!     &         +(CVCENT(2,ICVL)-SFCENT(2,ICFL))**2
!!     &         +(CVCENT(3,ICVL)-SFCENT(3,ICFL))**2)
!      FRSTCV(IBFL)=max(1.d-6,min(yy,yy2))
!      enddo
!  600 continue
!--------------------
! --- List Wall Face
!--------------------
      NBCWAL=0
      do 750 nb=1,nbcnd
      kdv=kdbcnd(1,nb)
      kd=kdbcnd(0,nb)
      if(kd==kdshutr) cycle
      if(
     &   kdv==kvnslp.or.
     &   kdv==kvlglw.or.
     &   kdv==kvrogf.or.
     &   kdv==kvsataW.or.
     &   kdv==kvmodl.or.
     &   kdv==kvEWF.or.
     
!     &   kdv==kvfslp.or.
     &   kd==kxnone) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        NBCWAL=NBCWAL+1
        enddo
      endif
 750  enddo
!
      if(NBCWAL.NE.NBCWAL_L) then
        write(ifle,*) 'error'
        write(ifle,*) 'NBCWAL is NOT equal to NBCWAL_L',NBCWAL,NBCWAL_L
        WRITE(ifle,*) ' ####    ','Please run [prefflow] program again'
        write(ifle,*) 'at subroutine :wall_distance'
        CALL FFRABORT(1,'in wall_distance')
      endif
!
! --- < 1. Measure Distance for ALL CV .
!
      LFUTAU(:)=0
      DISALL(:)=great
      if(NPE.eq.1) then
        call getfil(ifl,fnam,'distance')
        fnam='distance'
        if(NBCWAL.gt.0) then
          inquire(file=fnam,exist=wallfile)
          if(wallfile.and..not.(ical_mvmsh/=0)) then
 1200       write(ifll,*) ' ####    Wall Distance File existed, ',
     &      ' Will you calculate again and over write it ? (n/y)'
            read(*,*) recal
            if(recal=='y'.or.recal=='Y') then
              goto 1000
            elseif(recal=='n'.or.recal=='N') then
              open (ifl,FILE=fnam,FORM='UNFORMATTED',
     &            action='read',iostat=ios)
              read(ifl) NCVX
              if(NCVX.ne.NCV) then
                write(ifle,*) ' #### ERR: reading wall file error' 
                call FFRABORT(2,'wall_distance')
              endif
              read (ifl,ERR=1100) (DISALL(ICVL),ICVL=1,NCV)
              read (ifl,ERR=1100) (LFUTAU(ICVL),ICVL=1,NCV)
              close(ifl)
              goto 2000
 1100         call FFRABORT(1,'File [distance] reading error')
            else
              write(ifll,*) ' ####    Input [n] or [y]'
              goto 1200
            endif
          endif
 1000     continue
!
! --- Measure distance to no-slip wall
!
          if(.not.(ical_mvmsh/=0)) write(ifll,3010) my_rank
          do 700 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT.lt.0) goto 700
!          if(icaltb.eq.sles.or.cvdist(IMAT).or.icaltb.eq.SDES) then
          if(icaltb.eq.sles.or.cvdist(IMAT).or.icaltb.eq.SDES
     &      .or. icaltb==dles
     &      ) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do 740 ICVL=ICVS,ICVE
            dismin=GREAT
            dis1=GREAT
            NTAU=0
            do 720 nb=1,nbcnd
            kdv=kdbcnd(1,nb)
            if(kdv.eq.kvnslp.or.kdv.eq.kvlglw.or.kdv==kvrogf.or.
     &         kdv==kvmodl.or.kdv==kvEWF.or.kdv==kvsataW) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_INDEX(nb)
              do 710 IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              dx=SFCENT(1,ICFL)-CVCENT(1,ICVL)
              dy=SFCENT(2,ICFL)-CVCENT(2,ICVL)
              dz=SFCENT(3,ICFL)-CVCENT(3,ICVL)
              dis1=
     &        sqrt((CVCENT(1,ICVL)-SFCENT(1,ICFL))**2
     &            +(CVCENT(2,ICVL)-SFCENT(2,ICFL))**2
     &            +(CVCENT(3,ICVL)-SFCENT(3,ICFL))**2)
! 
              if(dis1.lt.dismin) then
                dismin=dis1
                NTAU=IBFL
                ICFMIN=ICFL
              endif
 710          continue
            endif
 720        continue
!
            if(NTAU.eq.0) goto 700
            DISALL(ICVL)=
     &      abs(SFAREA(1,ICFMIN)*(CVCENT(1,ICVL)-SFCENT(1,ICFMIN))
     &         +SFAREA(2,ICFMIN)*(CVCENT(2,ICVL)-SFCENT(2,ICFMIN))
     &         +SFAREA(3,ICFMIN)*(CVCENT(3,ICVL)-SFCENT(3,ICFMIN)))
            LFUTAU(ICVL)=NTAU
            IF(MOD(ICVL,2000).EQ.0) then
              if(.not.(ical_mvmsh/=0))  
     &        write(ifll,3000) ICVL,dismin,NTAU,my_rank
            endif
 740      continue
          endif
 700      continue
!
          open (ifl,FILE='distance',FORM='UNFORMATTED',action='write',
     &          iostat=ios)
          write(ifl) NCV
          write(ifl) (DISALL(ICVL),ICVL=1,MXCV)
          write(ifl) (LFUTAU(ICVL),ICVL=1,MXCV)
          if(.not.(ical_mvmsh/=0)) write(ifll,3020)  my_rank
          close(ifl)
        endif
      elseif(NPE.gt.1) then
!fortran90  ifl=30+my_rank
        call getfil(ifl,fnam,'distance')
        fnam=adjustl(WALLin)
        open (ifl,FILE=fnam,FORM='UNFORMATTED',action='read',iostat=ios)
        read (ifl,ERR=2100) NCVX
        if(NCVX/=NCV) stop 8877
        read (ifl,ERR=2100) (DISALL(ICVL),ICVL=1,NCV)
        read (ifl,ERR=2100) (LFUTAU(ICVL),ICVL=1,NCV)
        close(ifl)
      endif
 2000 continue
!------------------------
! --- Distance for First CV
!------------------------

      FRSTCV=1.D0
      do 600 nb=1,nbcnd
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      do IBFL=IBFS,IBFE
      ICFL=LBC_SSF(IBFL)
      ICVL=LVEDGE(1,ICFL)
!1)
      yy3=abs(SFAREA(1,ICFL)*(CVCENT(1,ICVL)-SFCENT(1,ICFL))
     &       +SFAREA(2,ICFL)*(CVCENT(2,ICVL)-SFCENT(2,ICFL))
     &       +SFAREA(3,ICFL)*(CVCENT(3,ICVL)-SFCENT(3,ICFL)))
!      yy2=3.d0*CVVOLM(ICVL)/SFAREA(ICFL,4)
!
      yy=CVVOLM(ICVL)**(1.d0/3.d0)
!
      dum1=max(1.d-16,min(yy,yy3))
!      dum1=0.5d0*CVVOLM(ICVL)/SFAREA(4,ICFL)

      FRSTCV(IBFL)=dum1
!      DISALL(ICVL)=dum1             !FRSTCV(IBFL)
      
      enddo 
  600 enddo
!
 3000 format(12X,'Processing ... NO. ',I8,' FINISHED,   Dis. & BC ',
     &4X,E10.4,I8,2X,' at my_rank',I4)
 3010 format(2X,28('*'),' Start measuring all Distance at my_rank= ',
     &I4,' !!! ',28('*'))
 3020 format(2X,28('*'),' End   measuring all Distance at my_rank= ',
     &I4,' !!! ',28('*'))
!
      return 
!
 2100 call FFRABORT(1,'Reading wall distance file error')
      end subroutine wall_distance
!   
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_geom(
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,
     &  LBC_SSF,LCYCSF,LCV_CV,locmsh,
     &  LVEDGE,LVRTCV,SFAREA,SFCENT,CVVOLM,CVCENT)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_boundary,only : LBC_INDEX,MAT_BCIDX,
     &                           nbcnd,kdbcnd,kdprdc,LBC_pair,kdsld
      use module_material,only : !KMAT_S,MAT_ALL,MAT_S,
     &                           nosld,nofld,
     &                           rg_mark
      use module_io,only       : ifle,ifll,lenfnm,getfil
!
!
      use module_rad   , only  : radflag,radmodel,gasmodel,ndiv,
     &                           paraflag,almn,blmn,awsgg,bwsgg,cwsgg,
     &                           ngauss
      use module_radsxf, only  : RDValue,RDIndex,RDN1,RDN2,RDId,
     &			         NRD,MRD,NRD0,MRD0,MID,StackOpt,ID,WK
      use module_radsxf, only  : RadHeat,RadHeatFlux,EB4,QSUM,NPALL,
     &                           WRad,FRAD,raddum
      use module_radsxf, only  : VJOBMD,NJOBMD
      use module_radsxf, ONLY  : IndexToRD,CoefToRD,MatBounEnd
      use module_radsxf, ONLY  : IndexToExpt,RaToExpt,ExptID,EBcomm
      use module_radsxf, only  : NDatExpt,NumTran,NumTranIn,NMatBoun
      use module_radsxf, only  : GWAbsorb, PAbsorb, PScatter,
     &				 radinten,sumcoef,P_rad
      use module_radgroup, only : IflagGroup
      use module_radsxf, only  : DivW,DivA
      use MODULE_ARRAY_FFLOW ,only : iterR,repsR,aepsR,errR
      use module_model,   only : vertex_cen,cell_cen,icon_cv
!
!
!
      use module_material,only : radmat,radfludflag
      use module_species,only  : spcnam,spcno
      use module_usersub, only : src_rad,para_rad,usryes
      use module_model,  only  : ical_prt
!
      implicit none
!
! --- [module arguments]
!
      integer,intent(out)      :: LVEDGE(2,  MXCVFAC)
      integer,intent(out)      :: LVRTCV(    MXALLCV)
      integer,intent(out)      :: MAT_NO(    0:MXMAT)
      integer,intent(out)      :: MAT_CV(    MXALLCV)
      integer,intent(out)      :: LCV_CV(    MXALLCV)
      integer,intent(out)      :: MAT_INDEX( 0:MXMAT)
      integer,intent(out)      :: MAT_CVEXT( 0:MXMAT)
      integer,intent(out)      :: MAT_DCIDX( 0:MXMAT)
      integer,intent(out)      :: MAT_CFIDX( 0:MXMAT)
      integer,intent(out)      :: LBC_SSF (MXSSFBC)
      integer,intent(out)      :: LCYCSF  (MXSSFBC)
      integer,intent(out)      :: locmsh  (MXSSFBC)

      real*8 ,intent(out)      :: SFAREA(4,MXCVFAC)
      real*8 ,intent(out)      :: SFCENT(3,MXCVFAC)
      real*8 ,intent(out)      :: CVVOLM(  MXALLCV)
      real*8 ,intent(out)      :: CVCENT(3,MXALLCV)
!
! --- [local entities]
!  
      integer :: nb,nbcndx,NMATX,IIMATX,ICVSX,ICVEX,ICFSX,ICFEX,
     &           IBFSX,IBFEX,IDCSX,IDCEX,ICFP
      integer :: ICVL,ICVS,ICVE,ICFL,ICFS,ICFE,IDCL,IDCS,IDCE
      integer :: IBFL,IBFS,IBFE,ios,I,ICV,neib
      integer :: NOUT,inum,istart,k,JCVS,JCVE
      integer :: IC,KMAT,IMAT,IIMAT,KIMAT,LINDEX,ICL,nbx,kd
      integer :: ifl,IIMAT_U,IMAT_U,ICVA,ICVB,radflagx,icon_cvX
      character(lenfnm),save :: fnam
      real*8  :: dx,dy,dz,dum1
!
      INTEGER :: IERR=0,IDIV
!
      MAT_NO(0)=0
      call getfil(ifl,fnam,'geom')
      if(NPE.gt.1) then
        fnam=adjustl(GEOMin)
      else
        fnam='geom'
      endif
!
      NMATX=0
      open (ifl,FILE=fnam,FORM='unformatted',status='unknown',
     &            iostat=ios)
      rewind(ifl)
      read(ifl)   icon_cvX
      if(icon_cvX/=icon_cv) then
        write(ifle,'(1X,2a)') 
     &   'ERR: Control Volume algorithm defined is ',
     &  'different from defined when running prefflow'
        write(ifle,'(1X,a)') 'MSG: re-running prefflow or :'
        if(icon_cvX==vertex_cen) then
          write(ifle,'(1X,a)') 'MSG: defining as: [&model/CV="vertex"]'
        else
          write(ifle,'(1X,a)') 'MSG: efining as: [&model/CV="cell"]'
        endif
        call FFRABORT
     &   (1,'STOP AT:Control Volume algorithm define error ')
      endif
      read(ifl)   NMATX
      read(ifl)  (MAT_NO(IIMAT),IIMAT=1,NMATX)
      read(ifl)  (MAT_INDEX(IIMAT),IIMAT=0,NMATX)
      read(ifl)  (MAT_CVEXT(IIMAT),IIMAT=0,NMATX)
      read(ifl)  (MAT_DCIDX(IIMAT),IIMAT=0,NMATX)
      read(ifl)  (MAT_CFIDX(IIMAT),IIMAT=0,NMATX)
      do 650 IIMAT=1,NMATX
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      read(ifl)  IIMATX,ICVSX,ICVEX
      if(ICVSX<=ICVEX) then
        read(ifl) (MAT_CV(ICVL),ICVL=ICVS,ICVE)
        read(ifl) (LVRTCV(ICVL),ICVL=ICVS,ICVE)
        read(ifl) (CVVOLM(ICVL),ICVL=ICVS,ICVE)
        read(ifl)((CVCENT(I,ICVL),I=1,3),ICVL=ICVS,ICVE)
      endif
      IDCS=MAT_DCIDX(IIMAT-1)+1
      IDCE=MAT_DCIDX(IIMAT)
      read(ifl)  IIMATX,IDCSX,IDCEX
      if(IDCSX<=IDCEX) then
        read(ifl) (MAT_CV(ICVL),ICVL=IDCS,IDCE)
        read(ifl) (LVRTCV(ICVL),ICVL=IDCS,IDCE)
        read(ifl) (CVVOLM(ICVL),ICVL=IDCS,IDCE)
        read(ifl)((CVCENT(I,ICVL),I=1,3),ICVL=IDCS,IDCE)
      endif
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      read(ifl)  IIMATX,ICFSX,ICFEX
      if(ICFSX<ICFEX) then
        read(ifl)((LVEDGE(I,ICFL),I=1,2),ICFL=ICFS,ICFE)
        read(ifl)((SFAREA(I,ICFL),I=1,4),ICFL=ICFS,ICFE)
        read(ifl)((SFCENT(I,ICFL),I=1,3),ICFL=ICFS,ICFE)
      endif
 650  continue
!
      read(ifl) nbcndx
      read(ifl)(LBC_INDEX(nb),nb=0,nbcndx)
      read(ifl)(LBC_pair(nb),nb=0,nbcndx)
      read(ifl)((MAT_BCIDX(nb,i),nb=1,nbcndx),i=1,2)
      DO nb=1,nbcndx
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)

      read(ifl) nbx,IBFSX,IBFEX
      if(IBFSX<=IBFEX) then
        read(ifl)(LBC_SSF(IBFL),IBFL=IBFS,IBFE)
        read(ifl)(LCYCSF(IBFL),IBFL=IBFS,IBFE)
      endif
      enddo
      read(ifl)(LCV_CV(ICV),ICV=1,NALLCV)
!
      READ(ifl) radflagx
      if(radflagx/=radflag) then
        call FFRABORT(1,'ERR: re-run prefflow for radiation')
      endif

!
!
      IF(ierr.NE.0) call FFRABORT(1,'read_geom--1')
!
! --- 
!
      IF(radflag.EQ.1) THEN
        ALLOCATE(WRad(NDIV), FRAD(NDIV),raddum(Ngauss))
        IF(radmodel(1:3).EQ.'FVM') THEN
          ALLOCATE(QSUM(200,4),
     &	  DivA(3,NDIV),DivW(NDIV),ID(NDIV),WK(NDIV),stat=ierr)
          IF(ierr.NE.0) call FFRABORT(1,'read_geom-0')
!          rg_mark=0
!          do IIMAT=1,NMAT    !ICV=1,NCV
!          IMAT=MAT_NO(IIMAT)
!          if(IMAT.LE.0) CYCLE
!          if(radfludflag(1,IMAT).eq.2) then 
!            rg_mark=1
!            exit
!          endif
!          enddo
!
          IF(rg_mark==1.and.para_rad.eq.usryes) THEN
            CALL  User_rad_GasPara(gasmodel,ParaFlag,almn,
     &			blmn,awsgg,bwsgg,cwsgg)   !  !zhangrad
          ENDIF
!
       ALLOCATE(iterR(NDIV))
       ALLOCATE(repsR(NDIV))
       ALLOCATE(aepsR(NDIV))
       ALLOCATE(errR(NDIV))
       IF(rg_mark.eq.1) THEN
         Ngauss_R=Ngauss
         NDIV_R=NDIV
         NDNG=Ngauss_R*NDIV_R
         MXALLNG=MXALLCV*Ngauss_R
       else
         Ngauss_R=Ngauss
         NDIV_R=NDIV
         NDNG=NDIV
         MXALLNG=MXALLCV
       ENDIF
!
       IF(rg_mark.eq.1) THEN
	    IF(gasmodel(1:3).eq.'SLW'.OR.gasmodel(1:4).eq.'WSGG') THEN
	      ALLOCATE(
     &                 radinten(NDNG,MXALLCV_RAD),
     &	               GWAbsorb(MXALLNG),
     &	               sumcoef(MXALLNG),stat=ierr)
              IF(ierr.NE.0) call FFRABORT(1,'read_geom-1')
            ELSE
              ALLOCATE(
     &                 radinten(NDNG,MXALLCV_RAD),
     &                 GWAbsorb(MXALLNG),
     &                 sumcoef(MXALLNG),stat=ierr)
              IF(ierr.NE.0) call FFRABORT(1,'read_geom-2')
            ENDIF
       ELSE
            ALLOCATE(
     &                 radinten(NDNG,MXALLCV_RAD),
     &               GWAbsorb(MXALLNG),
     &	             sumcoef(MXALLNG),stat=ierr)
            IF(ierr.NE.0) call FFRABORT(1,'read_geom-3')
       ENDIF
          
          IF(ierr/=0) THEN
            WRITE(ifle,*) 'Error in read_geom(): '
	    WRITE(ifle,*) '-1- can not allocate radiation arrays'
	    call FFRABORT(1,'read_geom')
          END IF
!		
          QSUM = 0.d0
	  radinten = 0.d0
	  RadHeat=0.d0
	  RadHeatFlux=0.d0
	  sumcoef=1.0
	  do idiv=1,ndiv
          CALL GetDivAngle(IDIV,NDIV,DivA(:,IDIV))   !  !zhangrad
          CALL GetDivWeight(IDIV,NDIV,DivW(IDIV))
	  end do
      ELSE
          READ(ifl,end=9999) StackOpt
	  READ(ifl) IflagGroup
	  READ(ifl) NRD,NRD0,MRD,MID
	  READ(ifl) NPALL,NDatExpt,NumTran,NMatBoun
	  WRITE(ifll,'(a,I6,a,I6)')  
     &	  'Vertex: Rad=',NRD0, ', fluid=',MAT_CVEXT(NMATX)
          IF(NMatBoun.NE.NMATX+nbcndx+1) THEN
     	    WRITE(ifle,*) 'Num of fluids and bounds different'
	    call FFRABORT(1,'read_geom')
      ENDIF
!
!Exchange area
!
      ALLOCATE(RDValue(MRD),RDIndex(NRD+1),stat=ierr)
      IF(ierr.NE.0) THEN
        WRITE(ifle,*) 'Error in read_geom(): '
        WRITE(ifle,*) '-2- can not allocate radiation arrays'
        call FFRABORT(1,'read_geom')
      END IF
!		
      READ(ifl) (RDIndex(I),I=1,NRD+1)
      READ(ifl) (RDValue(I),I=1,MRD)
!
	  IF(StackOpt(1:4).EQ.'BAND') THEN
            ALLOCATE(RDN1(NRD),RDN2(NRD),RDId(MRD),stat=ierr)
            IF(ierr.NE.0) THEN
              WRITE(ifle,*) 'Error in read_geom(): '
	      WRITE(ifle,*) '-3- can not allocate radiation arrays'
	      call FFRABORT(1,'read_geom')
	    ENDIF
!
            READ(ifl) (RDN1(I),I=1,NRD)
            READ(ifl) (RDN2(I),I=1,NRD)
            READ(ifl) (RDId(I),I=1,MRD)
          ENDIF
!property	
!VJobMD(:,5:6) is only for benchmark computation
          ALLOCATE(VJobMD(MAX(NPALL,MXALLCV_RAD),4),
     &	  NJOBMD(MAX(NPALL,MXALLCV_RAD),3),stat=ierr)
	  IF(ierr.NE.0) THEN
	    WRITE(ifle,*) 'Error in read_geom(): '
	    WRITE(ifle,*) '-4- can not allocate radiation arrays'
	    call FFRABORT(1,'read_geom')
	  END IF

	  READ(ifl) (NJobMD(I,1),I=1,NPALL)	!Vindex, ICV->IRD
	  READ(ifl) (NJobMD(I,2),I=1,NALLCV)  !Index, IRD->ICV
	  READ(ifl) (VJobMD(I,1),I=1,NPALL)	!AreaVol
	  READ(ifl) (VJobMD(I,2),I=1,NPALL)	!Prop
!
!data transfer between fluid and radiation
!
          ALLOCATE(MatBounEnd(0:NMatBoun),IndexToRD(2,NumTran),
     &	  CoefToRD(NumTran), stat=ierr)
	  IF(ierr.NE.0) THEN
	    WRITE(ifle,*) 'Error in read_geom(): '
	    WRITE(ifle,*) '-5- can not allocate radiation arrays'
	    call FFRABORT(1,'read_geom')
	  END IF
!
	  READ(ifl) (MatBounEnd(I),I=0,NMATX+nbcndx+1)
	  READ(ifl) (CoefToRD(I),I=1,NumTran)
	  READ(ifl) (IndexToRD(1,I),I=1,NumTran) !from ICV / IP
	  READ(ifl) (IndexToRD(2,I),I=1,NumTran) !to IRD
!
!data for communication
!
          IF(NDatExpt.GT.0) THEN
	    ALLOCATE(IndexToExpt(2,NDatExpt),RAToExpt(NDatExpt), 
     &	    ExptID(NRD),EBcomm(MXALLCV_RAD),stat=ierr)
	    IF(ierr.NE.0) THEN
	      WRITE(ifle,*) 'Error in read_geom(): '
	      WRITE(ifle,*) '-6- can not allocate radiation arrays'
              call FFRABORT(1,'read_geom')
	    END IF
!
            READ(ifl) (RaToExpt(I),I=1,NDatExpt)
	    READ(ifl) (IndexToExpt(1,I),I=1,NDatExpt)
            READ(ifl) (IndexToExpt(2,I),I=1,NDatExpt)			!to IRD
	    READ(ifl) (ExptID(I),I=1,NRD)						!mark IRD 
!		
            NJobMD(:,3)=0
            DO I=1,NDatExpt
	      NJobMD(IndexToExpt(1,I),3)=1
            END DO
            EBcomm=0.d0
            END IF
!		
            ALLOCATE(QSUM(200,4),EB4(NRD))
	    IF(ierr.NE.0) THEN
	      WRITE(ifle,*) 'Error in read_geom(): '
	      WRITE(ifle,*) '-7- can not allocate radiation arrays'
	      call FFRABORT(1,'read_geom')
	    ENDIF
            QSUM=0.d0
            RadHeat=0.d0
            RadHeatFlux=0.d0
	  END IF
      ENDIF


      close(ifl) 
!
      if(NPE.gt.1) then 
!----------------------
! --- read comm file
!----------------------
        call comm
!----------------------
! --- check CV Number
!----------------------
        do IIMAT=1,NMAT 
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        JCVE=MAT_INDEX(IIMAT)
        do neib= 1, NEIBPETOT
        istart=EXPORT_index(neib-1)
        inum  =EXPORT_index(neib)-istart
        do k=istart+1,istart+inum
        ICVL=EXPORT_item(k)
        if(ICVL.ge.ICVS.and.ICVL.le.ICVE) then
          if(ICVL.gt.JCVE) then
            write(ifle,*) ' ### ERR : CV Number sequence error' 
            write(ifle,*) ' Please contact FrontFlow developer' 
            call FFRABORT(1,'read_geom')
          endif
        endif
        enddo
        enddo
        enddo
!
        NOUT=0
        do IIMAT=1,NMAT 
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        JCVE=MAT_INDEX(IIMAT)
        NOUT=NOUT+ICVE-JCVE
        do neib= 1, NEIBPETOT
        istart=IMPORT_index(neib-1)          !EXPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart     !EXPORT_index(neib)-istart
        do k=istart+1,istart+inum
        ICVL=IMPORT_item(k)
        if(ICVL.ge.ICVS.and.ICVL.le.ICVE) then
          if(ICVL.lt.JCVE) then
            write(ifle,*) ' ### ERR : CV Number sequence error' 
            write(ifle,*) ' Please contact FrontFlow developer' 
            call FFRABORT(2,'read_geom')
          endif
          if(ICVL.gt.NCV) then
            write(ifle,*) ' ### ERR : CV Number sequence error' 
            write(ifle,*) ' Please contact FrontFlow developer' 
            call FFRABORT(3,'read_geom')
          endif
        endif
        enddo
        enddo
        enddo
        if(IMPORT_index(NEIBPETOT).ne.NOUT) then
          write(ifle,*) ' ### ERR : CV Number sequence error' 
          write(ifle,*) ' Please contact FrontFlow developer' 
          call FFRABORT(1,'read_geom')
        endif
!
      endif
!
! --- check material number
!
!      if(KMAT_S.ne.MXMAT) then
!        write(ifle,*) ' ### ERR : Material list wrong',KMAT_S,MXMAT
!        write(ifle,'(1x,2a,2I4)') 
!     &  ' ### ERR : Material Number NOT matched between ',
!     &  'fflow.ctl and user-grid fiel :',KMAT_S,MXMAT
!        call FFRABORT(2,'read_geom')
!      endif
!      MAT_S(:)=0
!      DO KMAT=1,KMAT_S
!      IMAT_U=MAT_ALL(KMAT)
!      do IIMAT=1,NMAT
!      IIMAT_U=MAT_NO(IIMAT)
!      if(IMAT_U.eq.IIMAT_U) then
!        MAT_S(KMAT)=IIMAT
!      endif
!      enddo
!      ENDDO
!

      return
 9910 WRITE(ifle,*) 'Error in read_geom(): '
      WRITE(ifle,*) '	 -1- can not allocate radiation arrays'
	
 9999 WRITE(ifle,*) 'Error in read_geom: read radiation data'
      WRITE(ifle,*) 'data not available, please cal. at PREFFLOW'
!
      return
!
      end subroutine read_geom
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine allocat_array
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      USE module_dimension
      use module_io,      only : ifle,ifll
      use module_metrix,only   : MHD_A,MHD_FAI,MHD_RVA,SIGMA,MHDbnd,
     &                           tempMHD,MHD_CRNT,MHD_CRT0,mat_coil,
     &                           MHD_FAI0
      use module_metrix,only   : movflg,cord,lacell,lvcell,
     &                           lvface,lbface,lcface,lfcell,LEFACE,
     &                           area,volume,gface,gcell,listbc
      use module_model,only    : ical_MHD,ical_mvmsh
!
      implicit none
!
! --- [local entities]
!
      integer :: ierr1=0
!
      allocate(MHD_A   (MXCV_MHD,3,2,2) ,stat=ierr1)
      allocate(MHD_FAI (MXCV_MHD,2,2)   ,stat=ierr1)
      allocate(MHD_RVA (MXFACE_MHD)     ,stat=ierr1)
      allocate(SIGMA   (MXCV_MHD)       ,stat=ierr1)
      allocate(MHDbnd  (MXBND_MHD,4,2)  ,stat=ierr1)
      allocate(tempMHD (MXCV_MHD,3,2)   ,stat=ierr1)
      allocate(MHD_CRNT(MXCV_MHD,3,2)   ,stat=ierr1)
      allocate(MHD_CRT0(MXCV_MHD,3,2)   ,stat=ierr1)
      allocate(MHD_FAI0 (MXCV_MHD,2,2)   ,stat=ierr1)
      allocate(mat_coil(MXMAT)   ,stat=ierr1)
!      allocate(GaugeC0(MXMAT)   ,stat=ierr1)
!
      if(ierr1/=0) call 
     &    FFRABORT(1,'ERR: allocate MHD arraies in allocat_array')
!
! --- Moving grid array
!
      allocate(movflg(  MXVRTX_m) ,stat=ierr1)
      allocate(cord  (3,MXVRTX_m) ,stat=ierr1)
      allocate(lacell(  MXCELL_m) ,stat=ierr1)
      allocate(lvcell(8,MXCELL_m) ,stat=ierr1)
      allocate(lvface(4,MXFACE_m) ,stat=ierr1)
      allocate(lbface(2,MXFACE_m) ,stat=ierr1)
      allocate(lcface(2,MXFACE_m) ,stat=ierr1)
      allocate(lfcell(7,MXCELL_m) ,stat=ierr1)
      allocate(LEFACE(5,MXFACE_m) ,stat=ierr1)
      allocate(area  (4,MXFACE_m) ,stat=ierr1)
      allocate(volume(  MXCELL_m) ,stat=ierr1)
      allocate(gface (3,MXFACE_m) ,stat=ierr1)
      allocate(gcell (3,MXCELL_m) ,stat=ierr1)
      allocate(listbc(4,MXSSFBC_old_m) ,stat=ierr1)
      if(ierr1/=0) call 
     &    FFRABORT(1,'ERR: allocate MOVE arraies in allocat_array')
!
      return
!
      end subroutine allocat_array
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_sliding(icall,iter,time,deltt,
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LBC_SSF,SFCENT,LCYCOLD,wifsld,OPPANG,
     &  LCYCSF,LVEDGE,CVVOLM,CVCENT,SFAREA)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,       only: ifll,ifle,lenfnm
      use module_material,only : rotati,ishaft,end,begin,nsplit,
     &                           ical_sld,rot_ang,sav_ang
      use module_boundary,only : kdsld,boundName,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           idis,LBC_pair,nblade
      use module_model, only   : monitor_stp
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,icall
      real*8 ,intent(in)    :: time,deltt
      integer,intent(in)    :: MAT_NO   ( 0:MXMAT)
      integer,intent(in)    :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)    :: MAT_DCIDX( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  ( 0:MXMAT)
      real*8 ,intent(in)    :: SFCENT(3,  MXCVFAC)
      integer,intent(in)    :: LBC_SSF  ( MXSSFBC)
      integer,intent(inout) :: LCYCSF   ( MXSSFBC)
      integer,intent(in)    :: LVEDGE  (2,MXCVFAC)
      real*8 ,intent(inout) :: CVVOLM  (  MXALLCV)
      real*8 ,intent(inout) :: CVCENT  (3,MXALLCV)
      real*8 ,intent(inout) :: SFAREA(4,  MXCVFAC)      
!
      integer,intent(inout) :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(inout) :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(inout) :: OPPANG     (MXSSFBC_SLD)
!
! --- [local entities]
!
      character(lenfnm),save :: fnam
      integer :: IIMAT1,IIMAT2,IMAT1,IMAT2
      integer :: isdB,isdA,NSLDA,NSLDB,IBFL,IBFS,IBFE,ICFL,ICFP
      integer :: IBFLS,ICFLS,ICFPS,ICFPSS,ICFL1A,ICFL2A,ICFL2B
      integer :: nb,kd,NSLD,IBF,i,ICV,IDC,ICVP,IDCP,ICFOLD
      integer :: ierr=0,ifl=-1,ios
      real*8  :: unit(3),shft(3),SFT,rps,sgn,alpha(2),th(3),dum1,dum2
      real*8  :: dmin,dmin1,rrr,XA,YA,ZA,XB,YB,ZB,
     &           costh1,sinth1,costh2,sinth2,
     &           xo,yo,zo,dro,drn
      real*8  :: org_x1,org_y1,org_x2,org_y2,r1(3),
     &           tha_pich,r_ang
      logical :: search1,search2
      integer :: IBFS1,IBFE1,IBFS2,IBFE2,ICFL1,ICFL2,IBFL1,IBFL2,
     &           IBFL2A,IBFL1A,ICVBO,IPICH
      integer,save :: I_SLD2=0
!
! --- 
!
      if((ical_sld==1.or.ical_sld==2.or.ical_sld==2)
     &  .and.MXSSFBC_SLD==1) then
        call FFRABORT
     &   (1,'ERR: Sliding Mesh dimen. error, re-run [prefflow]')
      endif
      if(deltt<0.d0) then
        LCYCOLD=LCYCSF
        I_SLD2=I_SLD2+1
      endif
      if(ical_sld==2.and.icall==2) then
        I_SLD2=I_SLD2+1
      endif
      if(I_SLD2>=3.and.ical_sld==2) return
!
      LCYCOLD=LCYCSF
      if(ierr/=0) call FFRABORT(1,'allocating error in sliding')
!
      alpha=0.d0
      
      if(mod(iter,monitor_stp)==0.and.my_rank==root) then
        write(ifll,'(2x,a,I4)') 
     &    '|    Sliding-Mesh calculating ...at my_rank=' ,my_rank
      endif
!
! --- 
!
      do 100 nb=1,nbcnd
      IIMAT1=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT1)) cycle
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT1=MAT_NO(IIMAT1)
      IMAT2=MAT_NO(IIMAT2)
      kd=kdbcnd(0,nb)
      if(kd==kdsld) then
        if(idis(nb)==0) then
          IBFS1=LBC_INDEX(nb-1)+1
          IBFE1=IBFS1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
          IBFS2=LBC_INDEX(nb-1)+1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2
          IBFE2=LBC_INDEX(nb)
        elseif(idis(nb)>=1) then
          IBFS1=LBC_INDEX(nb-1)+1
          IBFE1=LBC_pair(nb)
          IBFS2=LBC_pair(nb)+1
          IBFE2=LBC_INDEX(nb)
        endif
        search1=.false.
        search2=.false.
        alpha(1)=rot_ang(IMAT1)
        alpha(2)=rot_ang(IMAT2)
        if(ishaft(IMAT1)==1.and.IMAT1>0) then
          search1=.true.
        endif
!
        if(ishaft(IMAT2)==1.and.IMAT2>0) then
          search2=.true.
        endif
        search1=.true.
        search2=.true.
!
        costh1=dcos(alpha(1))
        sinth1=dsin(alpha(1))
        costh2=dcos(alpha(2))
        sinth2=dsin(alpha(2))

        org_x1=begin(1,IMAT1)*costh1-begin(2,IMAT1)*sinth1
        org_y1=begin(1,IMAT1)*sinth1+begin(2,IMAT1)*costh1
        org_x2=begin(1,IMAT2)*costh2-begin(2,IMAT2)*sinth2
        org_y2=begin(1,IMAT2)*sinth2+begin(2,IMAT2)*costh2
!
        if(search1.or.search2) then
          do IBFL1=IBFS1,IBFE1
            ICFL1=LBC_SSF(IBFL1)
            ICFOLD=LCYCOLD(IBFL1)
!
            xa=SFCENT(1,ICFL1)*costh1
     &        -SFCENT(2,ICFL1)*sinth1-org_x1
            ya=SFCENT(1,ICFL1)*sinth1
     &        +SFCENT(2,ICFL1)*costh1-org_y1
            za=SFCENT(3,ICFL1)          !only z-axis
! --- old
            xo=SFCENT(1,ICFOLD)*costh2
     &        -SFCENT(2,ICFOLD)*sinth2-org_x2
            yo=SFCENT(1,ICFOLD)*sinth2
     &        +SFCENT(2,ICFOLD)*costh2-org_y2
            zo=SFCENT(3,ICFOLD)         !only z-axis
            dmin=1.d15
            dmin1=1.d15
            ICFL2A=0
            IBFL2A=0
!
            do IPICH=1,nblade(nb,2)  !one-pitch
               
            tha_pich=2.d0*pi/dble(nblade(nb,2))*dble(IPICH-1)
            dum2=1.d0
            if(nblade(nb,2)>1) dum2=0.d0
            do IBFL2=IBFS2,IBFE2
              ICFL2=LBC_SSF(IBFL2)
              xb=SFCENT(1,ICFL2)*dcos(alpha(2)+tha_pich)        !costh2
     &          -SFCENT(2,ICFL2)*dsin(alpha(2)+tha_pich)-org_x2 !sinth2-org_x2
              yb=SFCENT(1,ICFL2)*dsin(alpha(2)+tha_pich)        !sinth2
     &          +SFCENT(2,ICFL2)*dcos(alpha(2)+tha_pich)-org_y2 !costh2-org_y2
              zb=SFCENT(3,ICFL2)        !only z-axis
              dum1=dum2*((xb-xa)*(xo-xa)+(yb-ya)*(yo-ya))   
              SFT =((xa-xb)**2+(ya-yb)**2+(za-zb)**2)  
              if( SFT.lt.dmin
     &           .and.dum1<=SML   !okabe
!    &           .and.ICFL2/=ICFOLD
     &          ) then
                dmin=SFT
                ICFL2A=ICFL2
                r_ang=alpha(2)+tha_pich
              endif

              if(SFT.lt.dmin1
!     &          .and.ICFL2/=ICFOLD
     &         ) then
                dmin1=SFT
                IBFL2A=ICFL2
                r_ang=alpha(2)+tha_pich
              endif
            enddo
            enddo
            LCYCSF(IBFL1)=ICFL2A
            if(ICFL2A==0) then
              LCYCSF(IBFL1)=IBFL2A
            endif
            if(idis(nb)==2) OPPANG(IBFL1)=r_ang
          enddo
!----------------------
! --- 2nd SLD
!----------------------
          NSLDB=0
          do IBFL2=IBFS2,IBFE2
            ICFL2=LBC_SSF(IBFL2)
            ICFOLD=LCYCOLD(IBFL2)
            NSLDB=NSLDB+1
            xa=SFCENT(1,ICFL2)*costh2
     &        -SFCENT(2,ICFL2)*sinth2-org_x2
            ya=SFCENT(1,ICFL2)*sinth2
     &        +SFCENT(2,ICFL2)*costh2-org_y2
            za=SFCENT(3,ICFL2) !only z-axis
!
            xo=SFCENT(1,ICFOLD)*costh1
     &        -SFCENT(2,ICFOLD)*sinth1-org_x1
            yo=SFCENT(1,ICFOLD)*sinth1
     &        +SFCENT(2,ICFOLD)*costh1-org_y1
            zo=SFCENT(3,ICFOLD)        !only z-axis
            dmin=1.d15
            dmin1=1.d15
            ICFL1A=0
            IBFL1A=0
!
            do IPICH=1,nblade(nb,1)  !one-pitch
            tha_pich=2.d0*pi/dble(nblade(nb,1))*dble(IPICH-1)
            dum2=1.d0
            if(nblade(nb,1)>1) dum2=0.d0
            do IBFL1=IBFS1,IBFE1
              ICFL1=LBC_SSF(IBFL1)
              xb=SFCENT(1,ICFL1)*cos(alpha(1)+tha_pich)        !costh1
     &          -SFCENT(2,ICFL1)*sin(alpha(1)+tha_pich)-org_x1 !sinth1-org_x1
              yb=SFCENT(1,ICFL1)*sin(alpha(1)+tha_pich)        !sinth1
     &          +SFCENT(2,ICFL1)*COS(alpha(1)+tha_pich)-org_y1 !costh1-org_y1 
              zb=SFCENT(3,ICFL1) !only z-axis
              dum1=dum2*((xb-xa)*(xo-xa)+(yb-ya)*(yo-ya))
              SFT =((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
              if(SFT.lt.dmin
!     &           .and.ICFL1/=ICFOLD
     &           .and.dum1<=SML
     &          ) then
                dmin=SFT
                ICFL1A=ICFL1
                r_ang=alpha(1)+tha_pich
              ENDIF
              if(SFT.lt.dmin1
!     &          .and.ICFL1/=ICFOLD
     &         ) then
                dmin1=SFT
                IBFL1A=ICFL1
                r_ang=alpha(1)+tha_pich
              endif
            enddo
            ENDDO
            LCYCSF(IBFL2)=ICFL1A
            if(ICFL1A==0)  LCYCSF(IBFL2)=IBFL1A
            if(idis(nb)==2) OPPANG(IBFL2)=r_ang
          enddo
!
          do IBFL1=IBFS1,IBFE1
            ICFL2A=LCYCSF(IBFL1)
            ICFOLD=LCYCOLD(IBFL1)
            if(ICFL2A==0) then
              LCYCSF(IBFL1)=ICFOLD
            endif
          enddo
!
          do IBFL2=IBFS2,IBFE2
            ICFL1A=LCYCSF(IBFL2)
            ICFOLD=LCYCOLD(IBFL2)
            if(ICFL1A==0) then
              LCYCSF(IBFL2)=ICFOLD
            endif
          enddo
!
!--------------------------------------------------------
! --- Nearest CV to LCYCSF; Second Near CV to LCYCOLD
!--------------------------------------------------------
!
          if(nblade(nb,1)==1) then
            do IBFL1=IBFS1,IBFE1
            ICFL=LBC_SSF(IBFL1)
            ICFOLD=LCYCOLD(IBFL1)
            ICFP=LCYCSF(IBFL1)
              if(ICFP==0.or.ICFOLD==0) 
     &          call FFRABORT(1,'ICFOLD=0 or ICFP=0')
            xa=SFCENT(1,ICFL)*costh1
     &        -SFCENT(2,ICFL)*sinth1-org_x1
            ya=SFCENT(1,ICFL)*sinth1
     &        +SFCENT(2,ICFL)*costh1-org_y1
            za=SFCENT(3,ICFL)          !only z-axis
            xo=SFCENT(1,ICFOLD)*costh2
     &        -SFCENT(2,ICFOLD)*sinth2-org_x2
            yo=SFCENT(1,ICFOLD)*sinth2
     &        +SFCENT(2,ICFOLD)*costh2-org_y2
            zo=SFCENT(3,ICFOLD)          !only z-axis
            xb=SFCENT(1,ICFP)*costh2
     &        -SFCENT(2,ICFP)*sinth2-org_x2
            yb=SFCENT(1,ICFP)*sinth2
     &        +SFCENT(2,ICFP)*costh2-org_y2
            zb=SFCENT(3,ICFP)         !only z-axis
            dum1=((xb-xa)**2+(yb-ya)**2+(zb-za)**2)
            dum2=((xo-xa)**2+(yo-ya)**2+(zo-za)**2)
            if(dum1<dum2) then
              LCYCSF(IBFL1)=ICFP
              LCYCOLD(IBFL1)=ICFOLD
            else
              LCYCSF(IBFL1)=ICFOLD
              LCYCOLD(IBFL1)=ICFP
            endif
            enddo

          else
            LCYCOLD(IBFS1:IBFE1)=LCYCSF(IBFS1:IBFE1)
          endif
!
          if(nblade(nb,2)==1) then
            do IBFL2=IBFS2,IBFE2
            ICFL=LBC_SSF(IBFL2)
            ICFOLD=LCYCOLD(IBFL2)
            ICFP=LCYCSF(IBFL2)
            if(ICFP==0.or.ICFOLD==0) 
     &      call FFRABORT(1,'ICFOLD=0 or ICFP=0')
            xa=SFCENT(1,ICFL)*costh2
     &        -SFCENT(2,ICFL)*sinth2-org_x2
            ya=SFCENT(1,ICFL)*sinth2
     &        +SFCENT(2,ICFL)*costh2-org_y2
            za=SFCENT(3,ICFL)          !only z-axis
            xo=SFCENT(1,ICFOLD)*costh1
     &        -SFCENT(2,ICFOLD)*sinth1-org_x1
            yo=SFCENT(1,ICFOLD)*sinth1
     &        +SFCENT(2,ICFOLD)*costh1-org_y1
            zo=SFCENT(3,ICFOLD)          !only z-axis
            xb=SFCENT(1,ICFP)*costh1
     &        -SFCENT(2,ICFP)*sinth1-org_x1
            yb=SFCENT(1,ICFP)*sinth1
     &        +SFCENT(2,ICFP)*costh1-org_y1
            zb=SFCENT(3,ICFP)        !only z-axis
            dum1=((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
            dum2=((xo-xa)**2+(yo-ya)**2+(zo-za)**2)
            if(dum1<dum2) then
              LCYCSF(IBFL2)=ICFP
              LCYCOLD(IBFL2)=ICFOLD
            else
              LCYCSF(IBFL2)=ICFOLD
              LCYCOLD(IBFL2)=ICFP
            endif
            enddo
          else
            LCYCOLD(IBFS2:IBFE2)=LCYCSF(IBFS2:IBFE2)
          endif
!----------------------
! --- Weighted coeff.          
!----------------------
          if(nblade(nb,1)==1) then
            do IBFL1=IBFS1,IBFE1
            ICFL=LBC_SSF(IBFL1)
            ICFP=LCYCSF(IBFL1)
            ICFOLD=LCYCOLD(IBFL1)
            xa=SFCENT(1,ICFL)*costh1
     &        -SFCENT(2,ICFL)*sinth1-org_x1
            ya=SFCENT(1,ICFL)*sinth1
     &        +SFCENT(2,ICFL)*costh1-org_y1
            za=SFCENT(3,ICFL)         !only z-axis
!
            xb=SFCENT(1,ICFP)*costh2
     &        -SFCENT(2,ICFP)*sinth2-org_x2
            yb=SFCENT(1,ICFP)*sinth2
     &        +SFCENT(2,ICFP)*costh2-org_y2
            zb=SFCENT(3,ICFP)         !only z-axis
!
            xo=SFCENT(1,ICFOLD)*costh2
     &        -SFCENT(2,ICFOLD)*sinth2-org_x2
            yo=SFCENT(1,ICFOLD)*sinth2
     &        +SFCENT(2,ICFOLD)*costh2-org_y2
            zo=SFCENT(3,ICFOLD)          !only z-axis
!
            dro=dsqrt((xo-xa)**2+(yo-ya)**2+(zo-za)**2)   !now
            drn=dsqrt((xb-xa)**2+(yb-ya)**2+(zb-za)**2)   !old
!
            wifsld(IBFL1)=dro/(drn+dro+SML)
!
            enddo
          else
            wifsld(IBFS1:IBFE1)=0.5d0
          endif
!
          if(nblade(nb,2)==1) then
            do IBFL2=IBFS2,IBFE2
            ICFL=LBC_SSF(IBFL2)
            ICFP=LCYCSF(IBFL2)
            ICFOLD=LCYCOLD(IBFL2)
            xa=SFCENT(1,ICFL)*costh2
     &        -SFCENT(2,ICFL)*sinth2-org_x2
            ya=SFCENT(1,ICFL)*sinth2
     &        +SFCENT(2,ICFL)*costh2-org_y2
            za=SFCENT(3,ICFL)        !only z-axis
!
            xb=SFCENT(1,ICFP)*costh1
     &        -SFCENT(2,ICFP)*sinth1-org_x1
            yb=SFCENT(1,ICFP)*sinth1
     &        +SFCENT(2,ICFP)*costh1-org_y1
            zb=SFCENT(3,ICFP)        !only z-axis
!
            xo=SFCENT(1,ICFOLD)*costh1
     &        -SFCENT(2,ICFOLD)*sinth1-org_x1
            yo=SFCENT(1,ICFOLD)*sinth1
     &        +SFCENT(2,ICFOLD)*costh1-org_y1
            zo=SFCENT(3,ICFOLD)      !only z-axis
!
            dro=dsqrt((xo-xa)**2+(yo-ya)**2+(zo-za)**2)
            drn=dsqrt((xb-xa)**2+(yb-ya)**2+(zb-za)**2)
            wifsld(IBFL2)=dro/(drn+dro+SML)
            enddo
          else
            wifsld(IBFS2:IBFE2)=0.5d0
          endif
        endif
      endif
 100  continue
!      
      if(mod(iter,monitor_stp)==0.and.my_rank==root) then
        write(ifll,'(2x,a,I4)') 
     &  '|    Sliding-Mesh calculating end !!!',my_rank
      endif
!
      return
!
      end subroutine list_sliding
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_sliding_V(icall,iter,time,deltt,
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LBC_SSF,SFCENT,LCYCOLD,wifsld,OPPANG,
     &  LCYCSF,LVEDGE,CVVOLM,CVCENT,SFAREA)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,       only: ifll,ifle,lenfnm
      use module_material,only : rotati,ishaft,end,begin,nsplit,
     &                           ical_sld,rot_ang,sav_ang
      use module_boundary,only : kdsld,boundName,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           idis,LBC_pair,nblade
      use module_metrix,only   : tmpfac=>d2vect,tmpsld!,mask
      use module_model,only    : ical_vect,monitor_stp
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,icall
      real*8 ,intent(in)    :: time,deltt
      integer,intent(in)    :: MAT_NO   ( 0:MXMAT)
      integer,intent(in)    :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)    :: MAT_DCIDX( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  ( 0:MXMAT)
      real*8 ,intent(in)    :: SFCENT(3,  MXCVFAC)
      integer,intent(in)    :: LBC_SSF  ( MXSSFBC)
      integer,intent(inout) :: LCYCSF   ( MXSSFBC)
      integer,intent(in)    :: LVEDGE  (2,MXCVFAC)
      real*8 ,intent(inout) :: CVVOLM  (  MXALLCV)
      real*8 ,intent(inout) :: CVCENT  (3,MXALLCV)
      real*8 ,intent(inout) :: SFAREA(4,  MXCVFAC)      
!
      integer,intent(inout) :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(inout) :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(inout) :: OPPANG     (MXSSFBC_SLD)
!
! --- [local entities]
!
      integer :: IIMAT1,IIMAT2,IMAT1,IMAT2
      integer :: isdB,isdA,NSLDA,NSLDB,IBFL,IBFS,IBFE,ICFL,ICFP
      integer :: IBFLS,ICFLS,ICFPS,ICFPSS,ICFL1A,ICFL2A,ICFL2B
      integer :: nb,kd,NSLD,IBF,i,ICV,IDC,ICVP,IDCP,ICFOLD
      integer :: ierr=0,ifl=-1,ios
      real*8  :: unit(3),shft(3),SFT,rps,sgn,alpha(2),th(3),dum1,dum2
      real*8  :: dmin,dmin1,rrr,XA,YA,ZA,XB,YB,ZB,
     &           costh1,sinth1,costh2,sinth2,
     &           xo,yo,zo,dro,drn
      real*8  :: org_x1,org_y1,org_x2,org_y2,r1(3),
     &           tha_pich,r_ang
      integer :: IBFS1,IBFE1,IBFS2,IBFE2,ICFL1,ICFL2,IBFL1,IBFL2,
     &           IBFL2A,IBFL1A,ICVBO,IPICH,NO
      integer,save :: I_SLD2=0
      integer :: IMLC(1)
      logical,save :: MSK
!
! --- 
!
      if((ical_sld==1.or.ical_sld==2.or.ical_sld==2)
     &  .and.MXSSFBC_SLD==1) then
        call FFRABORT
     &   (1,'ERR: Sliding Mesh dimen. error, re-run [prefflow]')
      endif
      if(deltt<0.d0) then
        LCYCOLD=LCYCSF
        I_SLD2=I_SLD2+1
      endif
      if(ical_sld==2.and.icall==2) then
        I_SLD2=I_SLD2+1
      endif
      if(I_SLD2>=3.and.ical_sld==2) return
!
      LCYCOLD=LCYCSF
      if(ierr/=0) call FFRABORT(1,'allocating error in sliding')
!
      alpha=0.d0
      
      if(mod(iter,monitor_stp)==0.and.my_rank==root) then
        write(ifll,'(2x,a,I4)') 
     &    '|    Sliding-Mesh calculating ...at my_rank=' ,my_rank
      endif
!
! --- 
!
      do 100 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if(kd==kdsld) then
        IIMAT1=MAT_BCIDX(nb,1)
        IIMAT2=MAT_BCIDX(nb,2)
        IMAT1=MAT_NO(IIMAT1)
        IMAT2=MAT_NO(IIMAT2)
!
        IBFS1=LBC_INDEX(nb-1)+1
        IBFE1=LBC_pair(nb)
        IBFS2=LBC_pair(nb)+1
        IBFE2=LBC_INDEX(nb)
        if(IBFS1>IBFE1)  cycle
        if(IBFS2>IBFE2)  cycle
        alpha(1)=rot_ang(IMAT1)
        alpha(2)=rot_ang(IMAT2)
!
        costh1=dcos(alpha(1))
        sinth1=dsin(alpha(1))
        costh2=dcos(alpha(2))
        sinth2=dsin(alpha(2))
!
        org_x1=begin(1,IMAT1)*costh1-begin(2,IMAT1)*sinth1
        org_y1=begin(1,IMAT1)*sinth1+begin(2,IMAT1)*costh1
        org_x2=begin(1,IMAT2)*costh2-begin(2,IMAT2)*sinth2
        org_y2=begin(1,IMAT2)*sinth2+begin(2,IMAT2)*costh2
!------------------------
! --- 1st SLD
!------------------------
        do IBFL2=IBFS2,IBFE2
           ICFL2=LBC_SSF(IBFL2)
           tmpsld(IBFL2,7)=     ! xb=
     &           SFCENT(1,ICFL2)*costh2
     &          -SFCENT(2,ICFL2)*sinth2-org_x2
           tmpsld(IBFL2,8)=! yb=
     &           SFCENT(1,ICFL2)*sinth2
     &          +SFCENT(2,ICFL2)*costh2-org_y2
           tmpsld(IBFL2,9)=! zb=
     &           SFCENT(3,ICFL2)      !only z-axis
        enddo

!
        do 1000 IBFL1=IBFS1,IBFE1
        ICFL1=LBC_SSF(IBFL1)
        ICFOLD=LCYCOLD(IBFL1)
        xa=SFCENT(1,ICFL1)*costh1
     &    -SFCENT(2,ICFL1)*sinth1-org_x1
        ya=SFCENT(1,ICFL1)*sinth1
     &    +SFCENT(2,ICFL1)*costh1-org_y1
        za=SFCENT(3,ICFL1)          !only z-axis

        xo=SFCENT(1,ICFOLD)*costh2
     &    -SFCENT(2,ICFOLD)*sinth2-org_x2
        yo=SFCENT(1,ICFOLD)*sinth2
     &    +SFCENT(2,ICFOLD)*costh2-org_y2
        zo=SFCENT(3,ICFOLD)         !only z-axis
!
        IMLC(1)=0
!
        do 1100 IBFL2=IBFS2,IBFE2
           tmpsld(IBFL2,1)=  !dum1 
     &                    ((tmpsld(IBFL2,7)-xa)
     &                    *(xo-xa)
     &                    +(tmpsld(IBFL2,8)-ya)
     &                    *(yo-ya))
           tmpsld(IBFL2,2)=   !SFT
     &                    ((xa-tmpsld(IBFL2,7))**2
     &                    +(ya-tmpsld(IBFL2,8))**2
     &                    +(za-tmpsld(IBFL2,9))**2) 
 1100   enddo
!
        do IBFL2=IBFS2,IBFE2   !zhang
           if(tmpsld(IBFL2,1)>SML) then
           tmpsld(IBFL2,2)=1.d10
           endif
        enddo
!
        NO=IBFE2-IBFS2+1
        OPPANG=1.d10
        OPPANG(1:NO)=tmpsld(IBFS2:IBFE2,2)
        IMLC(1)=MinLOC(OPPANG,dim=1)
        if(IMLC(1)>0) then
          LCYCSF(IBFL1)=LBC_SSF(IBFS2-1+IMLC(1))
        endif
 1000 enddo
!----------------------
! --- 2nd SLD
!----------------------
!
      do IBFL1=IBFS1,IBFE1
      ICFL1=LBC_SSF(IBFL1)
      tmpsld(IBFL1,7)=! xb=
     &         SFCENT(1,ICFL1)*costh1
     &        -SFCENT(2,ICFL1)*sinth1-org_x1
      tmpsld(IBFL1,8)=! yb=
     &         SFCENT(1,ICFL1)*sinth1
     &        +SFCENT(2,ICFL1)*costh1-org_y1 
      tmpsld(IBFL1,9)=! zb=
     &         SFCENT(3,ICFL1) !only z-axis
      enddo
!
      do 2000 IBFL2=IBFS2,IBFE2
      ICFL2=LBC_SSF(IBFL2)
      ICFOLD=LCYCOLD(IBFL2)
      xa=SFCENT(1,ICFL2)*costh2
     &        -SFCENT(2,ICFL2)*sinth2-org_x2
      ya=SFCENT(1,ICFL2)*sinth2
     &        +SFCENT(2,ICFL2)*costh2-org_y2
      za=SFCENT(3,ICFL2) !only z-axis
!
      xo=SFCENT(1,ICFOLD)*costh1
     &        -SFCENT(2,ICFOLD)*sinth1-org_x1
      yo=SFCENT(1,ICFOLD)*sinth1
     &        +SFCENT(2,ICFOLD)*costh1-org_y1
      zo=SFCENT(3,ICFOLD)        !only z-axis
!
      IMLC(1)=0
!
      DO 2100 IBFL1=IBFS1,IBFE1
      ICFL1=LBC_SSF(IBFL1)
       tmpsld(IBFL1,1)=  !dum1
     &                ((tmpsld(IBFL1,7)-xa)
     &                *(xo-xa)
     &                +(tmpsld(IBFL1,8)-ya)
     &                *(yo-ya))
       tmpsld(IBFL1,2)=  !SFT =
     &                ((xa-tmpsld(IBFL1,7))**2
     &                +(ya-tmpsld(IBFL1,8))**2
     &                +(za-tmpsld(IBFL1,9))**2)
 2100 enddo
!
      do IBFL1=IBFS1,IBFE1   !zhang
        if(tmpsld(IBFL1,1)>SML) then
          tmpsld(IBFL1,2)=1.d10
        endif
      enddo
!
      NO=IBFE1-IBFS1+1
      OPPANG=1.d10
!
      OPPANG(1:NO)=tmpsld(IBFS1:IBFE1,2)
      IMLC(1)=MINLOC(OPPANG,dim=1)
!
      if(IMLC(1)>0) then
        LCYCSF(IBFL2)=LBC_SSF(IBFS1-1+IMLC(1))
      endif
 2000 enddo
!
      do IBFL1=IBFS1,IBFE1
      ICFL2A=LCYCSF(IBFL1)
      ICFOLD=LCYCOLD(IBFL1)
      if(ICFL2A==0) then
        LCYCSF(IBFL1)=ICFOLD
      endif
      enddo

!
!
      do IBFL2=IBFS2,IBFE2
      ICFL1A=LCYCSF(IBFL2)
      ICFOLD=LCYCOLD(IBFL2)
      if(ICFL1A==0) then
        LCYCSF(IBFL2)=ICFOLD
      endif
      enddo
!
!--------------------------------------------------------
! --- Nearest CV to LCYCSF; Second Near CV to LCYCOLD
!--------------------------------------------------------
!
      if(nblade(nb,1)==1) then
        do IBFL1=IBFS1,IBFE1
        ICFL=LBC_SSF(IBFL1)
        ICFOLD=LCYCOLD(IBFL1)
        ICFP=LCYCSF(IBFL1)
        if(ICFP==0.or.ICFOLD==0) 
     &      call FFRABORT(1,'ICFOLD=0 or ICFP=0')
        tmpsld(IBFL1,7)=  !xa=
     %                  SFCENT(1,ICFL)*costh1
     &                 -SFCENT(2,ICFL)*sinth1-org_x1
        tmpsld(IBFL1,8)=  !ya=
     %                  SFCENT(1,ICFL)*sinth1
     &                 +SFCENT(2,ICFL)*costh1-org_y1
        tmpsld(IBFL1,9)=  !za=
     %                  SFCENT(3,ICFL)          !only z-axis
!
        tmpsld(IBFL1,1)=  !xo=
     %                  SFCENT(1,ICFOLD)*costh2
     &                 -SFCENT(2,ICFOLD)*sinth2-org_x2
        tmpsld(IBFL1,2)=  !yo=
     %                  SFCENT(1,ICFOLD)*sinth2
     &                 +SFCENT(2,ICFOLD)*costh2-org_y2
        tmpsld(IBFL1,3)=  !zo=
     %                  SFCENT(3,ICFOLD)          !only z-axis
!
        tmpsld(IBFL1,4)=  !xb=
     %                  SFCENT(1,ICFP)*costh2
     &                 -SFCENT(2,ICFP)*sinth2-org_x2
        tmpsld(IBFL1,5)=  !yb=
     %                  SFCENT(1,ICFP)*sinth2
     &                 +SFCENT(2,ICFP)*costh2-org_y2
        tmpsld(IBFL1,6)=  !zb=
     %                  SFCENT(3,ICFP)         !only z-axis
        enddo
!
        do IBFL1=IBFS1,IBFE1
        ICFL=LBC_SSF(IBFL1)
        ICFOLD=LCYCOLD(IBFL1)
        ICFP=LCYCSF(IBFL1)
        dum2=((tmpsld(IBFL1,1)-tmpsld(IBFL1,7))**2
     &       +(tmpsld(IBFL1,2)-tmpsld(IBFL1,8))**2
     &       +(tmpsld(IBFL1,3)-tmpsld(IBFL1,9))**2)

        dum1=((tmpsld(IBFL1,4)-tmpsld(IBFL1,7))**2
     &       +(tmpsld(IBFL1,5)-tmpsld(IBFL1,8))**2
     &       +(tmpsld(IBFL1,6)-tmpsld(IBFL1,9))**2)

        if(dum1<dum2) then
          LCYCSF(IBFL1)=ICFP
          LCYCOLD(IBFL1)=ICFOLD
        else
          LCYCSF(IBFL1)=ICFOLD
          LCYCOLD(IBFL1)=ICFP
        endif
        enddo
!
        do IBFL1=IBFS1,IBFE1
        ICFL=LBC_SSF(IBFL1)
        ICFP=LCYCSF(IBFL1)
        ICFOLD=LCYCOLD(IBFL1)
        dro=dsqrt((tmpsld(IBFL1,1)-tmpsld(IBFL1,7))**2
     &           +(tmpsld(IBFL1,2)-tmpsld(IBFL1,8))**2
     &           +(tmpsld(IBFL1,3)-tmpsld(IBFL1,9))**2)   !now

        drn=dsqrt((tmpsld(IBFL1,4)-tmpsld(IBFL1,7))**2
     &           +(tmpsld(IBFL1,5)-tmpsld(IBFL1,8))**2
     &           +(tmpsld(IBFL1,6)-tmpsld(IBFL1,9))**2)   !old
        wifsld(IBFL1)=dro/(drn+dro+SML)
        enddo
      else
        LCYCOLD(IBFS1:IBFE1)=LCYCSF(IBFS1:IBFE1)
        wifsld(IBFS1:IBFE1)=0.5d0
      endif
!---------------------
! --- 
!---------------------
      if(nblade(nb,2)==1) then
         do IBFL2=IBFS2,IBFE2
         ICFL=LBC_SSF(IBFL2)
         ICFOLD=LCYCOLD(IBFL2)
         ICFP=LCYCSF(IBFL2)
         if(ICFP==0.or.ICFOLD==0) 
     &      call FFRABORT(1,'ICFOLD=0 or ICFP=0')
!
         tmpsld(IBFL2,7)=  !xa=
     &                 SFCENT(1,ICFL)*costh2
     &                -SFCENT(2,ICFL)*sinth2-org_x2
         tmpsld(IBFL2,8)=  !ya=
     &                 SFCENT(1,ICFL)*sinth2
     &                +SFCENT(2,ICFL)*costh2-org_y2
         tmpsld(IBFL2,9)=  !za=
     &                 SFCENT(3,ICFL)          !only z-axis

         tmpsld(IBFL2,1)=  !xo=
     &                  SFCENT(1,ICFOLD)*costh1
     &                 -SFCENT(2,ICFOLD)*sinth1-org_x1
         tmpsld(IBFL2,2)=  !yo=
     &                  SFCENT(1,ICFOLD)*sinth1
     &                 +SFCENT(2,ICFOLD)*costh1-org_y1
         tmpsld(IBFL2,3)=  !zo=
     &                 SFCENT(3,ICFOLD)          !only z-axis
!
         tmpsld(IBFL2,4)=  !xb=
     &                 SFCENT(1,ICFP)*costh1
     &                 -SFCENT(2,ICFP)*sinth1-org_x1
         tmpsld(IBFL2,5)=  !yb=
     &                 SFCENT(1,ICFP)*sinth1
     &                 +SFCENT(2,ICFP)*costh1-org_y1
         tmpsld(IBFL2,6)=  !zb=
     &                 SFCENT(3,ICFP)        !only z-axis
!
         enddo
!
         do IBFL2=IBFS2,IBFE2
         ICFL=LBC_SSF(IBFL2)
         ICFOLD=LCYCOLD(IBFL2)
         ICFP=LCYCSF(IBFL2)
         dum2=((tmpsld(IBFL2,1)-tmpsld(IBFL2,7))**2
     &        +(tmpsld(IBFL2,2)-tmpsld(IBFL2,8))**2
     &        +(tmpsld(IBFL2,3)-tmpsld(IBFL2,9))**2)

         dum1=((tmpsld(IBFL2,4)-tmpsld(IBFL2,7))**2
     &        +(tmpsld(IBFL2,5)-tmpsld(IBFL2,8))**2
     &        +(tmpsld(IBFL2,6)-tmpsld(IBFL2,9))**2)

         if(dum1<dum2) then
           LCYCSF(IBFL2)=ICFP
           LCYCOLD(IBFL2)=ICFOLD
         else
           LCYCSF(IBFL2)=ICFOLD
           LCYCOLD(IBFL2)=ICFP
         endif
         enddo
!
         do IBFL2=IBFS2,IBFE2
         ICFL=LBC_SSF(IBFL2)
         ICFP=LCYCSF(IBFL2)
         ICFOLD=LCYCOLD(IBFL2)
         dro=dsqrt((tmpsld(IBFL2,1)-tmpsld(IBFL2,7))**2
     &            +(tmpsld(IBFL2,2)-tmpsld(IBFL2,8))**2
     &            +(tmpsld(IBFL2,3)-tmpsld(IBFL2,9))**2)   !now

         drn=dsqrt((tmpsld(IBFL2,4)-tmpsld(IBFL2,7))**2
     &            +(tmpsld(IBFL2,5)-tmpsld(IBFL2,8))**2
     &            +(tmpsld(IBFL2,6)-tmpsld(IBFL2,9))**2)   !old
         wifsld(IBFL2)=dro/(drn+dro+SML)
         enddo
       else
         LCYCOLD(IBFS2:IBFE2)=LCYCSF(IBFS2:IBFE2)
         wifsld(IBFS2:IBFE2)=0.5d0
       endif
!----------------------
! --- Weighted coeff.          
!----------------------
      endif
 100  continue
!
      if(mod(iter,monitor_stp)==0.and.my_rank==root) then
        write(ifll,'(2x,a,I4)') 
     &  '|    Sliding-Mesh calculating end !!!',my_rank
      endif
!
      return
!
      end subroutine list_sliding_V
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_overset(
     &  icall,iter,time,deltt,
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LBC_SSF,SFCENT,
!     &  LCYCOLD,wifsld,OPPANG,
     &  LCYCSF,LVEDGE,CVVOLM,CVCENT,SFAREA)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,       only: ifll,ifle,lenfnm
      use module_material,only : rotati,ishaft,end,begin,nsplit,
     &                           ical_sld,rot_ang,sav_ang
      use module_boundary,only : kdsld,boundName,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           idis,LBC_pair,nblade,kdovst
      use module_model, only   : monitor_stp,ical_mvmsh
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,icall
      real*8 ,intent(in)    :: time,deltt
      integer,intent(in)    :: MAT_NO   ( 0:MXMAT)
      integer,intent(in)    :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)    :: MAT_DCIDX( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  ( 0:MXMAT)
      real*8 ,intent(in)    :: SFCENT(3,  MXCVFAC)
      integer,intent(in)    :: LBC_SSF  ( MXSSFBC)
      integer,intent(inout) :: LCYCSF   ( MXSSFBC)
      integer,intent(in)    :: LVEDGE  (2,MXCVFAC)
      real*8 ,intent(inout) :: CVVOLM  (  MXALLCV)
      real*8 ,intent(inout) :: CVCENT  (3,MXALLCV)
      real*8 ,intent(inout) :: SFAREA(4,  MXCVFAC)      
!
!      integer,intent(inout) :: LCYCOLD    (MXSSFBC_SLD)
!      real*8 ,intent(inout) :: wifsld     (MXSSFBC_SLD)
!      real*8 ,intent(inout) :: OPPANG     (MXSSFBC_SLD)
!
! --- [local entities]
!
      integer :: IIMAT1,IIMAT2,IMAT1,IMAT2
      integer :: isdB,isdA,NSLDA,NSLDB,IBFL,IBFS,IBFE,ICFL,ICFP
      integer :: IBFLS,ICFLS,ICFPS,ICFPSS,ICFL1A,ICFL2A,ICFL2B
      integer :: nb,kd,NSLD,IBF,i,ICV,IDC,ICVP,IDCP,ICFOLD
      integer :: ierr=0,ifl=-1,ios
      integer :: ICVS,ICVE,ICVL,ICVL_M
      real*8  :: unit(3),shft(3),SFT,rps,sgn,alpha(2),th(3),dum1,dum2
      real*8  :: dmin,dmin1,rrr,XA,YA,ZA,XB,YB,ZB,xyz(3),
     &           costh1,sinth1,costh2,sinth2,
     &           xo,yo,zo,dro,drn
      real*8  :: org_x1,org_y1,org_x2,org_y2,r1(3),
     &           tha_pich,r_ang
      logical :: search1,search2
      integer :: IBFS1,IBFE1,IBFS2,IBFE2,ICFL1,ICFL2,IBFL1,IBFL2,
     &           IBFL2A,IBFL1A,ICVBO,IPICH
      integer,save :: I_SLD2=0
!
! --- 
!
      alpha=0.d0
      if(deltt<0.d0) then
      else
        if(ical_sld==2.and.ical_mvmsh==0) return
      endif
!      
      if(mod(iter,monitor_stp)==0.and.my_rank==root) then
        write(ifll,'(2x,a,I4)') 
     &    '|    Overset-Mesh calculating ...at my_rank=' ,my_rank
      endif
!--------------------
! --- rotation case
!--------------------
      do 100 nb=1,nbcnd
      IIMAT1=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT1)) cycle
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT1=MAT_NO(IIMAT1)
      IMAT2=MAT_NO(IIMAT2)
      kd=kdbcnd(0,nb)
!
      if(kd==kdovst) then
        if(IMAT1==IMAT2) then
          write(ifle,*) 
     &     'ERR: Same Material No. on two side of Overset BC'
           write(ifll,*) 'MSG: Overset BC No= ',nb
           write(ifll,*) 
     &     'MSG: re-run [prefflow] or Call your FFR supportor'
           call FFRABORT(1,'ERR: in list_overset') 
        endif
!
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
!
        search1=.true.
        search2=.true.
!
        alpha(1)=rot_ang(IMAT1)
        alpha(2)=rot_ang(IMAT2)
!
        costh1=cos(alpha(1))
        sinth1=sin(alpha(1))
        costh2=cos(alpha(2))
        sinth2=sin(alpha(2))
!
        org_x1=begin(1,IMAT1)*costh1-begin(2,IMAT1)*sinth1
        org_y1=begin(1,IMAT1)*sinth1+begin(2,IMAT1)*costh1
        org_x2=begin(1,IMAT2)*costh2-begin(2,IMAT2)*sinth2
        org_y2=begin(1,IMAT2)*sinth2+begin(2,IMAT2)*costh2
!
        if(search1.or.search2) then 
          do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            do i=1,3
            shft(i)=(CVCENT(i,ICV)-SFCENT(i,ICFL))
            enddo
            dum1=shft(1)*SFAREA(1,ICFL)
     &          +shft(2)*SFAREA(2,ICFL)
     &          +shft(3)*SFAREA(3,ICFL)
!         
            do i=1,3
!              unit(i)=SFCENT(i,ICFL)-dum1*SFAREA(i,ICFL) 
              unit(i)=SFCENT(i,ICFL) 
            enddo 
            xa=unit(1)*costh1
     &        -unit(2)*sinth1
            ya=unit(1)*sinth1
     &        +unit(2)*costh1
            za=unit(3)           !only z-axis rotation 
!
            xyz(1)=CVCENT(1,IDC)
            xyz(2)=CVCENT(2,IDC)
            xyz(3)=CVCENT(3,IDC)

            ICVS=MAT_CVEXT(IIMAT2-1)+1
            ICVE=MAT_CVEXT(IIMAT2)
            ICVL_M=0
            dmin=1.d10
            do ICVL=ICVS,ICVE
            xo=CVCENT(1,ICVL)*costh2
     &        -CVCENT(2,ICVL)*sinth2-org_x2
            yo=CVCENT(1,ICVL)*sinth2
     &        +CVCENT(2,ICVL)*costh2-org_y2
            zo=CVCENT(3,ICVL)        !only z-axis
            dum2=(xa-xo)**2+(ya-yo)**2+(za-zo)**2
            if(dum2<dmin) then
              dmin=dum2
              ICVL_M=ICVL
              xyz(1)=xo
              xyz(2)=yo
              xyz(3)=zo
            endif
            enddo            
            if(ICVL_M/=0) then
              LCYCSF(IBFL)=ICVL_M
!              CVCENT(1,IDC)=xyz(1)
!              CVCENT(2,IDC)=xyz(2)
!              CVCENT(3,IDC)=xyz(3)
!              CVCENT(1,IDC)=SFCENT(1,ICFL)
!              CVCENT(2,IDC)=SFCENT(2,ICFL)
!              CVCENT(3,IDC)=SFCENT(3,ICFL)
            endif
          enddo
        endif
      endif
 100  enddo
!
      write(ifll,'(2x,a,I4)') 
     &  '|    Overset-Mesh calculating end !!!',my_rank
!
      return
!
      end subroutine list_overset
!
