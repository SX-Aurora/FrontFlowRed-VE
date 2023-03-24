!
!      subroutine write_result
!      subroutine write_restart
!      subroutine write_force_history
!      subroutine write_anim
!      subroutine write_cdcl
!      subroutine write_probe
!      subroutine write_fluid_force
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_result(
     & iter,time,delt,ffexit,outputfile,vlasum,
     & pp0,prs,rho,tmp,yys,vel,aks,rmut,rva,WDOT,ccc,
     & MHD_A,MHD_FAI,MHD_CRNT,MHD_CRT0,BBB,
     & LORENZ,
     & rho2,tmp2,rmut2,vel2,yys2,prs2,
     & SUFBND,SDOT,MOLFRC,BLKTHK,
     & mat_cal,lcycold,wifsld,oppang,wiface,
     & LCV_CV,LVEDGE,LBC_SSF,LCYCSF,
     & MAT_CV,MAT_NO,MAT_CVEXT,MAT_INDEX,CVCENT,
     & rmu,utau,SFAREA,SFCENT,FRSTCV,CVVOLM,
     & cord,lacell,lvcell,
     & POTNAL,POFLUX,POTNBC,FIELD_U,DISALL,WDRBND,
     & ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    ICV=MAT_CV(ICVL);  ICVG=NOD_IDg(ICV)  ;LCV_CV(ICV)=ICVL
!
!    ICV : Local CV number in HPC
!    ICVG: Globle CV number (1 CPU)
!    ICVL: material-resorted Local CV number in HPC
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only       : getfil,lenfnm,ifle,ifll
      use module_output,only   : outchk,fnmlst,none,yes,mlt_rslt,
     &                           BC_output
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,kdcvd,
     &                           kxnone,kdbuff,boundName,kdolet,
     &                           phs_idx,phs_dns,IDX_SUFRAC,
     &                           phs_nbc,phs_typ,phs_snm,phs_nam,
     &                           gasphase,surphase,blkphase,phs_com,
     &                           phs_inifrac,num_site,blk_dns,surfreac,
     &				 kdintr,kdsld
      use module_les   ,only   : ISTART,n_ave,ista_momery,
     &                           iuvw_ave_rms_re,
     &                           ip_ave,imu_ave,it_ave,imaxmin,
     &                           icomp_ave,irans_ave,
     &                           it_rms,ip_rms,imu_rms,
     &                           icomp_rms,irans_rms,
     &                           iuvwt_re,
     &                           iuvw_rans_re,
     &                           SDOT_suf,WDOT_gas,depoeth_spd,
     &                           rat_E_S,num_ratout,molefr_suf,
     &                           mol_rslt,blk_thick
     &                           ,iuvwc_re
      use module_species,only  : spcnam,wm,r_wm
      use module_scalar,only   : sclname,rns_scl,iaph,ivold,icavi,
     &                           ical_cavi,ivof
      use module_model, only   : idrdp,incomp,mach0,icaltb,KE2S,
     &                           ke,sles,dles,lles,dns,RSM,SDES,
     &                           ical_field,idrdp,comp,u_func,uf_mm,
     &                           ical_prt
      use module_Euler2ph,only : ieul2ph
      USE module_movegrid,only : ical_mov
      use module_material,only : ical_sld,rotati,ishaft,end,begin,
     &                           rot_ang,nofld,nosld
      use module_chemreac,ONLY : ical_suf,nneq,vreq
      use module_metrix,only   : vlr =>d2work1,multiR
      use module_metrix,only   : temp=>W1K1
      use module_fluidforce,only : ForceFlag,iwallvar
      use module_model,   only   : ical_MHD,ical_mvmsh
      use module_scalar,only   : pot_scl,ip_pefc,potname
      use module_scalar,only   : calgeq,calxi,ihflmlt,igeq,ixi,
     &                           ORG,idifflm
      use module_VOF     ,only : ical_vof
      use module_rad, only     : radflag,radmodel
      use module_radsxf, only  : RadHeat,RadHeatFlux
      use module_time,    only : steady
      use module_metrix,only   : STA_SAV
      use module_model,   only : vertex_cen,cell_cen,icon_cv
      use module_model,   only : ical_WDR
      use module_boundary,only : partkd,kpart_s
      USE module_particle,only : starttime,Rh_d,func_d
      use module_trace,only    : PRT_VF
      
!
!  1. Write result file
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter
      real*8 ,intent(in)    :: time,delt
      logical,intent(in)    :: ffexit,outputfile
      integer,intent(in)    :: LCV_CV (   MXALLCV)
      integer,intent(in)    :: LVEDGE (2, MXCVFAC)
      INTEGER,INTENT(IN)    :: LBC_SSF(   MXSSFBC)
      integer,intent(in)    :: LCYCSF (   MXSSFBC)
      integer,intent(in)    :: MAT_CV (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_NO (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_INDEX( 0:MXMAT)
      real*8 ,intent(in)    :: pp0    (   MXALLCV)
      real*8 ,intent(in)    :: prs    (   MXALLCV)
      real*8 ,intent(in)    :: rho    (   MXALLCV)
      real*8 ,intent(in)    :: tmp    (   MXALLCV)
      real*8 ,intent(in)    :: yys    (   MXALLCV,MXcomp)
      real*8 ,intent(inout) :: vel    (   MXALLCV,3)
      real*8 ,intent(in)    :: aks    (   MXALLCVR,MXrans)
      real*8 ,intent(in)    :: rmut   (   MXALLCV)
      real*8 ,intent(in)    :: rva    (   MXCVFAC)
      real*8 ,intent(in)    :: WDOT   (   MXALLCV,MXCOMP)
      real*8 ,intent(in)    :: CCC    (        MXALLCV)
      real*8 ,intent(in)    :: RHO2   (   MXALLCVC)
      real*8 ,intent(in)    :: tmp2   (   MXALLCV2)
      real*8 ,intent(in)    :: RMUT2  (   MXALLCV2)
      real*8 ,intent(in)    :: VEL2   (   MXALLCV2,3)
      real*8 ,intent(in)    :: YYS2   (   MXALLCV2,MXCOMP)
      real*8 ,intent(in)    :: prs2   (   MXALLCV2)
      real*8 ,intent(inout) :: vlasum (   MXCV)
      real*8 ,intent(in)    :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: SUFBND (MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(IN)    :: SDOT   (MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(IN)    :: MOLFRC (MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(IN)    :: BLKTHK (MXSSFBC_SUF,MXPHASE)
      real*8 ,intent(in)    :: rmu    (   MXALLCV)
      real*8 ,intent(in)    :: utau   ( 0:MXSSFBC)
      real*8 ,intent(in)    :: SFAREA (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT (3, MXCVFAC)
      real*8 ,intent(in)    :: FRSTCV (   MXSSFBC)
      real*8 ,intent(in)    :: MHD_A   (  MXCV_MHD,3,2)
      real*8 ,intent(in)    :: LORENZ  (  MXCV_MHD,3,2)
      real*8 ,intent(in)    :: MHD_FAI (  MXCV_MHD,2)
      real*8 ,intent(in)    :: MHD_CRNT(  MXCV_MHD,3,2)
      real*8 ,intent(in)    :: MHD_CRT0(  MXCV_MHD,3,2)
      real*8 ,intent(in)    :: BBB     (  MXCV_MHD,3,2)
      real*8 ,intent(in)    :: cord    (3,MXVRTX_m)
      integer,intent(in)    :: lacell  (  MXCELL_m)
      integer,intent(in)    :: lvcell  (8,MXCELL_m)
      REAL*8,INTENT(IN)  :: POTNAL  (  MXALLCVP,MXPOTN)
      REAL*8,INTENT(IN)  :: POFLUX  (  MXCVFACP,MXPOTN)
      REAL*8,INTENT(IN)  :: POTNBC  (  MXSSFBCP,MXPOTN)
      real*8 ,intent(in) :: FIELD_U   (    MXCV_D,NFLID)
      REAL*8 ,INTENT(IN) :: CVVOLM(        MXALLCV)
!
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)    
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
      REAL*8 ,INTENT(IN)    :: WIFACE(        MXCVFAC)
      REAL*8 ,INTENT(IN)    :: DISALL(        MXCV)
      REAL*8 ,INTENT(IN)    :: WDRBND(        MXCV_WDR,N_inj)
!      
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      character(lenfnm),save :: fname,fname1,fnam,fnam1,fnamsuf,fnamgrid
      integer,save :: ifl=-1,iflmvg,ierr=0,ifl1,nsdot=0,ndeps=0,
     &                neq=0,iflsuf,nmolfr=0,nthick=0
      integer,parameter :: noutvar=500
      character(20) :: NM_HDN_B(3,2),NM_HDN_C(3,2),NM_HDN_C0(3,2),
     &                 NM_HDN_A(3,2),NM_HDN_L(3,2)
      character(20) :: NM_MHD_F(2)
      character(80) :: nameVar(noutvar)=' '
      character(8) :: text,NMHDN(3,2)
!
      logical :: fexist=.false.,unstdy=.false.
      integer :: i,j,l,n,out,nrnsx,numVar,icom,iicom,ierr1=0,IMD
      integer :: npotnx,NFLIDx
      integer :: IIMAT,IMAT,ICV,ICF,ios,ipcom,icoms,icome,iph,ii
      integer :: ntp,nbx,idum,ieq,idum2,IRC
      integer :: ICVL,ICVS,ICVE,nb,kd,IBFS,IBFE,IBFL,ICFL,IDC
!
      integer :: iuvw_ave_rms_rex,IMHD,kdp
      integer :: ip_avex,it_avex,imu_avex,imaxminx,ical_prtx
      integer :: ip_rmsx,it_rmsx,imu_rmsx,ista_momeryx
      integer :: iuvwt_rex,IMAT_U
      real*8  :: unit(3),tha,bb(3,3),rbb(3,3),radi(3),vr(3),dr,v0(3)
      real*8  :: sinth,costh,dum1,dum2,dum3
!
      integer,allocatable :: icomp_avex(:),irans_avex(:)      
      integer,allocatable :: icomp_rmsx(:),irans_rmsx(:)
      integer,allocatable :: iuvw_rans_rex(:)
      integer,allocatable :: iuvwc_rex(:)
      integer,allocatable :: amii(:),ami(:)
      character*80        :: nameWallVar(11)=' '
      integer             :: numWallVar,isw
      real(8)             :: uwall(1:3),uwallnrm,ydis
      real*8              :: yp,rnu
      real*8,parameter    :: SML1=1.d-8
      INTEGER :: ICVP,ICFP,IDCP,in_ave
      integer :: iin_ave
      integer,save :: ibot=0
!
      do nb=1,nbcnd
      if(trim(boundName(nb,1))=='bot') then
        ibot=nb
      endif
      enddo
!
!
!
      ALLOCATE(ami(mxcomp),amii(mxcompall),stat=ierr1)
!
      ierror=0
!
      unstdy=icaltb.eq.sles.or.icaltb.eq.dles.or.
     &       icaltb.eq.lles.or.icaltb.eq.dns.or.icaltb.eq.SDES.or.
     &       .not.steady
!
      call outchk(fnmlst(2),iter,time,delt,out)
      if(out.eq.none) then
        DEALLOCATE(ami,amii)
        return
      endif
!
      if(ifl.lt.0) then
        if(NPE.gt.1) then
          call getfil(ifl,fnam,'result')
          fnam=RESLTout
          if(fnam.eq.' ') call FFRABORT(1,'write_result')
!
          call getfil(iflmvg,fnam1,'move_grid')
          fnam1=RESLTgrid
          if(fnam1.eq.' ') call FFRABORT(1,'write_move_grid')
        else
          call getfil(ifl,fnam,'result')
          if(fnam.eq.' ') 
     &        call FFRABORT(1,'call getfil at write_result')
!     
          call getfil(iflmvg,fnam1,'move_grid')
          if(fnam1.eq.' ') 
     &        call FFRABORT(1,'call getfil at write_move_grid')
        endif
        fname=adjustl(fnam)
        fnam=TRIM(TRIM(fname)//'.'//'frontflow')
        fname=trim(fnam)
!
        fname1=adjustl(fnam1)
        fnamgrid=TRIM(TRIM(fname1)//'.'//'frontflow')
        fname1=trim(fnamgrid)
        
      endif
      
      if(.not.ffexit.and.out.ne.yes.and..not.outputfile) then
         DEALLOCATE(ami,amii)
         return
      endif
!
!
!
!
      IF(idrdp/=incomp.and.(calgeq.or.calxi).and.idifflm==ORG) THEN
      call flamelet_ys_org
     &   (MAT_CVEXT,MAT_CAL,MAT_NO,aks,rho,tmp,yys)        
      ENDIF
      

!
! --- 
!
      if(ical_sld>0.or.mlt_rslt==yes.or.(ical_mvmsh/=0)) then
        write(text,'(i8)') iter
        text=adjustl(text)
        fnam=TRIM(TRIM(fname)//'_'//TRIM(text))
        fnamgrid=TRIM(TRIM(fname1)//'_'//TRIM(text))
      endif
!
! --- print out moving grid coordinates
!
      if(ical_mvmsh/=0) then
        if(NPE>1) then
          call FFRABORT(1,'ERR: NOT support HPC move grid')
!        call output_move_grid_hpc(iflmvg,fnamgrid,cord,lacell,lvcell)
        else
          call output_move_grid(iflmvg,fnamgrid,cord,lacell,lvcell)
        endif
      endif
      if(ical_mvmsh==2.or.ical_mvmsh==3.or.ical_mvmsh==4) then
        DEALLOCATE(ami,amii)
        return
      endif
!
!-< 1. Print out result >-
!
!
!-< 1.1. Output result into external file >-
!
                       nrnsx=nrans
      if(.not.rns_scl) nrnsx=0
                       npotnx=npotn
      if(.not.pot_scl) npotnx=0
                        NFLIDx=NFLID
      if(ical_field==0) NFLIDx=0
      
!
!------------------------------------
! --- Specified numVar & nameVar(:) 
!------------------------------------
!
      nameVar(1)='pressure'
      nameVar(2)='back_prs'
      nameVar(3)='densty'
      nameVar(4)='temperature'
      nameVar(5)='tur_mu'
      nameVar(6)='velo_u;Velocity'
      nameVar(7)='velo_v'
      nameVar(8)='velo_w'
      numVar=5+3
      if(mol_rslt==1) then
        do i=1,ncomp
        numVar=numVar+1
        nameVar(numVar)='GasMolFr_'//TRIM(adjustl(spcnam(i)))
        enddo
      else
        do i=1,ncomp
        numVar=numVar+1
        nameVar(numVar)='GasMasFr_'//TRIM(adjustl(spcnam(i)))
        enddo
      endif
      if(ibot/=0) then
        numVar=numVar+1
        nameVar(numVar)='Z_bot'
      endif
!-----------------------------------------
! --- MHD 
!-----------------------------------------
      if(ical_MHD==1.or.ical_MHD==2) then
!
        NM_HDN_B(1,1)='MHD_Brx;MHD_Br'
        NM_HDN_B(2,1)='MHD_Bry'
        NM_HDN_B(3,1)='MHD_Brz'
        NM_HDN_B(1,2)='MHD_Bix;MHD_Bi'
        NM_HDN_B(2,2)='MHD_Biy'
        NM_HDN_B(3,2)='MHD_Biz'
!
        NM_HDN_C0(1,1)='MHD_J0rx;MHD_J0r'
        NM_HDN_C0(2,1)='MHD_J0ry'
        NM_HDN_C0(3,1)='MHD_J0rz'
        NM_HDN_C0(1,2)='MHD_J0ix;MHD_J0i'
        NM_HDN_C0(2,2)='MHD_J0iy'
        NM_HDN_C0(3,2)='MHD_J0iz'        
!
        NM_HDN_C(1,1)='MHD_Jrx;MHD_Jr'
        NM_HDN_C(2,1)='MHD_Jry'
        NM_HDN_C(3,1)='MHD_Jrz'
        NM_HDN_C(1,2)='MHD_Jix;MHD_Ji'
        NM_HDN_C(2,2)='MHD_Jiy'
        NM_HDN_C(3,2)='MHD_Jiz'
!
        
!
        NM_HDN_A(1,1)='MHD_Arx;MHD_Ar'
        NM_HDN_A(2,1)='MHD_Ary'
        NM_HDN_A(3,1)='MHD_Arz'
        NM_HDN_A(1,2)='MHD_Aix;MHD_Ai'
        NM_HDN_A(2,2)='MHD_Aiy'
        NM_HDN_A(3,2)='MHD_Aiz'
!
        NM_HDN_L(1,1)='MHD_Lrx;LORENZr'
        NM_HDN_L(2,1)='MHD_Lry'
        NM_HDN_L(3,1)='MHD_Lrz'
!
        NM_HDN_L(1,2)='MHD_Lix;LORENZi'
        NM_HDN_L(2,2)='MHD_Liy'
        NM_HDN_L(3,2)='MHD_Liz'
!
        NM_MHD_F(1)  ='MHD_FAIr'
        NM_MHD_F(2)  ='MHD_FAIi'
!
        DO IMHD=1,2
        do i=1,3
        numVar=numVar+1
        nameVar(numVar)=NM_HDN_B(i,IMHD)
        enddo
        do i=1,3
        numVar=numVar+1
        nameVar(numVar)=NM_HDN_C0(i,IMHD)
        enddo
        do i=1,3
        numVar=numVar+1
        nameVar(numVar)=NM_HDN_C(i,IMHD)
        enddo

        do i=1,3
        numVar=numVar+1
        nameVar(numVar)=NM_HDN_A(i,IMHD)
        enddo

        do i=1,3
        numVar=numVar+1
        nameVar(numVar)=NM_HDN_L(i,IMHD)
        enddo

        numVar=numVar+1
        nameVar(numVar)=NM_MHD_F(IMHD)
        enddo
      endif
!-----------------------------------------
! --- surface module output
!-----------------------------------------
      if(ical_suf==1) then
        if(SDOT_suf==1) then
          iicom=0
          do nb=1,nbcnd
!
          ami(:)=0
          do IRC=1,nneq
          if(IDX_SUFRAC(nb,IRC)) then
            iph=1
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            do  ipcom=icoms,icome
            icom=phs_com(ipcom)
            if(vreq(ICOM,2,IRC)>SML.or.vreq(ICOM,1,IRC)>SML) then
              ami(icom)=1
            endif
            enddo
          endif
          enddo
!    
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1.or.surfreac(nb)==2) then
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==gasphase) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ami(icom)==1)then
                iicom=iicom+1
                numVar=numVar+1
                nameVar(numVar)=
     & 'SurMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/m^2/s]_'//
     &    trim(boundName(nb,1))
              endif
              enddo
            elseif(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              numVar=numVar+1
              nameVar(numVar)=
     & 'SurMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/m^2/s]_'//
     &    trim(phs_nam(iph))
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              numVar=numVar+1
              nameVar(numVar)=
     & 'SurMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/m^2/s]_'//
     &    trim(phs_nam(iph))
              enddo
            endif
            enddo
          endif
          enddo
          nsdot=iicom
        endif
!
        if(molefr_suf==1) then
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              numVar=numVar+1
              nameVar(numVar)=
     &   'SiteMolFr_'//TRIM(adjustl(spcnam(icom)))//
     &    '_'//trim(phs_nam(iph))
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              numVar=numVar+1
              nameVar(numVar)=
     &    'BulkMolFr_'//TRIM(adjustl(spcnam(icom)))//
     &    '_'//trim(phs_nam(iph))
              enddo
            endif
            enddo
          endif
          enddo
          nmolfr=iicom
!          if(nmolfr/=NCOMP_SUF) then
!            call FFRABORT(1,'ERR: nmolfr/=NCOMP_SUF')
!          endif
        endif
!
        if(depoeth_spd==1.or.blk_thick==1) then
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase.and.nbx==nb) then
!              do ipcom=icoms,icome
!              icom=phs_com(ipcom)
!              iicom=iicom+1
!              numVar=numVar+1
!              nameVar(numVar)=
!     &        TRIM(adjustl(spcnam(icom)))//'_[m/s]'
!              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              numVar=numVar+1
              nameVar(numVar)=
     &   'BulkEH/DP_'//TRIM(adjustl(spcnam(icom)))//'_[m/s]_'//
     &    trim(phs_nam(iph))
              enddo
            endif
            enddo
          endif
          enddo
          ndeps=iicom
        endif
!
!        if(ndeps/=NCOMP_SUF) then
!          call FFRABORT(1,'ERR: ndeps/=NCOMP_SUF')
!        endif
!
!
!
        if(blk_thick==1) then
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase.and.nbx==nb) then
            elseif(ntp==blkphase.and.nbx==nb) then
              iicom=iicom+1
              numVar=numVar+1
              nameVar(numVar)=
     &   'Bulk_Thickness_'//trim(phs_nam(iph))//'_[m]_'//
     &    trim(phs_nam(iph))
            endif
            enddo
          endif
          enddo
          nthick=iicom
        endif
        
      endif
!
      if(ical_WDR==1) then
        do ICOM=1,N_inj+1
          write(text,'(i4)') ICOM
          numVar=numVar+1
          if(ICOM<=N_inj) then
            nameVar(numVar)='WDR_'//trim(text)//'[%]'
          else
            nameVar(numVar)='WDR_total_[%]'
          endif
        enddo
      endif
!
      if(WDOT_gas==1) then
          do icom=1,ncomp
          numVar=numVar+1
          nameVar(numVar)=
     & 'GasMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/(s*m^3)]'
          enddo
      endif
!
      if(num_ratout>0) then
          iicom=0
          amii=0
          do ieq=1,nneq
          do icom=1,ncompall
            do i=1,num_ratout
            idum=rat_E_S(i,icom)
            if(idum==ieq) then
              iicom=iicom+1
              amii(iicom)=i
              write(text,'(i3)') ieq
              text=adjustl(text)
              numVar=numVar+1
              nameVar(numVar)=
     &    'EqMolRate_'//TRIM(adjustl(spcnam(icom)))
     &         //'_eq_'//TRIM(text)//'_[mole/(s*m^3)]'
            endif
            enddo
          enddo
          enddo
          neq=iicom
          if(neq/=num_ratout) then
            call FFRABORT(1,'ERR: neq/=num_ratout')
          endif
      endif
!
!-----------------------------------------
! --- Two-Phase flow 
!-----------------------------------------
!
      if(ieul2ph>0) then 
        numVar=numVar+1
        if(ieul2ph==1) then
          nameVar(numVar)='averaged_density'
        else
          nameVar(numVar)='pressure2'
        endif
        numVar=numVar+1
        nameVar(numVar)='densty2'
        numVar=numVar+1
        nameVar(numVar)='temperature2'
        numVar=numVar+1
        nameVar(numVar)='tur_mu2'
        numVar=numVar+1
        nameVar(numVar)='velo_u2;Velocity2'
        numVar=numVar+1
        nameVar(numVar)='velo_v2'
        numVar=numVar+1
        nameVar(numVar)='velo_w2'
        do i=1,ncomp
        numVar=numVar+1
        nameVar(numVar)=TRIM(adjustl(spcnam(i)))//'2'
        enddo
      endif
!
      
!
      do i=1,nrnsx
        numVar=numVar+1
        nameVar(numVar)=TRIM(adjustl(sclname(i)))
      enddo
!
      do i=1,npotnx
        numVar=numVar+1
        nameVar(numVar)=TRIM(adjustl(potname(i)))
      enddo
!
      if(radflag==1) then
	numVar=numVar+1
        nameVar(numVar)='Raddivq'
	numVar=numVar+1
        nameVar(numVar)='Radqw'
      endif
!
      if(ical_field>0) then
        do i=1,NFLIDx
        write(text,'(i8)') i
        numVar=numVar+1
        nameVar(numVar)='FIELD'//'_'//TRIM(adjustl(text))
        enddo
      endif
!
      IF(idrdp==comp) then
        numVar=numVar+1
        nameVar(numVar)='Mach_Number'
        numVar=numVar+1
        nameVar(numVar)='Sound_speed'
      endif
!
      if(NPE.gt.1) then
        numVar=numVar+1
        nameVar(numVar)='HPC_region'
      endif
!
      if(ical_prt==1) then
        numVar=numVar+1
        nameVar(numVar)='particle_VF'
      endif
!
      if(u_func(2)>=1.or.u_func(3)==1) then
        numVar=numVar+1
        nameVar(numVar)='Multi_region'
      endif
!
      if(NPE.gt.1) then
        fnam1=TRIM(statis)
      else
        fnam1='statis'
      endif
      inquire(file=fnam1,exist=fexist)
      if(iter.gt.ISTART.and.(fexist.or.ista_momery==1).and.unstdy) then 
        if(iuvw_ave_rms_re.eq.1) then 
          numVar=numVar+1
          nameVar(numVar)='aver_u;Average_velocity'
          numVar=numVar+1
          nameVar(numVar)='aver_v'
          numVar=numVar+1
          nameVar(numVar)='aver_w'
          numVar=numVar+1
          nameVar(numVar)='RMS_u '
          numVar=numVar+1
          nameVar(numVar)='RMS_v '
          numVar=numVar+1
          nameVar(numVar)='RMS_w '
          numVar=numVar+1
          nameVar(numVar)='Re_uv '
          numVar=numVar+1
          nameVar(numVar)='Re_vw '
          numVar=numVar+1
          nameVar(numVar)='Re_uw '
        endif
!
        if(ip_ave==1) then
          numVar=numVar+1
          nameVar(numVar)='aver_p'
        endif
!
        if(ip_rms.eq.1) then
          numVar=numVar+1
          nameVar(numVar)='RMS_p '
        endif
!
        if(imu_ave==1) then
          numVar=numVar+1
          nameVar(numVar)='aver_Mu'
        endif
!
        if(imu_rms.eq.1) then
          numVar=numVar+1
          nameVar(numVar)='RMS_Mu'
        endif
! 
        if(imaxmin==1) then
          numVar=numVar+1
          nameVar(numVar)='MAX_u;MAX_velocity'
          numVar=numVar+1
          nameVar(numVar)='MAX_v'
          numVar=numVar+1
          nameVar(numVar)='MAX_w'
          numVar=numVar+1
          nameVar(numVar)='MIN_u;MIN_velocity'
          numVar=numVar+1
          nameVar(numVar)='MIN_v'
          numVar=numVar+1
          nameVar(numVar)='MIN_w'
          numVar=numVar+1
          nameVar(numVar)='MAX_p'
          numVar=numVar+1
          nameVar(numVar)='MIN_p'
        endif
!
        if(it_ave==1) then
          numVar=numVar+1
          nameVar(numVar)='aver_T'
        endif
!
        if(it_rms.eq.1) then
          numVar=numVar+1
          nameVar(numVar)='RMS_T '
        endif
!
        if(iuvwt_re.eq.1) then
          numVar=numVar+1
          nameVar(numVar)='Re_uT '
          numVar=numVar+1
          nameVar(numVar)='Re_vT '
          numVar=numVar+1
          nameVar(numVar)='Re_wT '
        endif
!
        do i=1,ncomp
        if(icomp_ave(i).eq.1) then
          numVar=numVar+1
          nameVar(numVar)='AVE_'//TRIM(adjustl(spcnam(i)))
        endif
        if(icomp_rms(i).eq.1)then
          numVar=numVar+1
          nameVar(numVar)='RMS_'//TRIM(adjustl(spcnam(i)))
        endif
        enddo
!
        ! Tagami
        do i=1,ncomp
           if(iuvwc_re(i).eq.1) then
              numVar=numVar+1
              nameVar(numVar)='AVE_u_'//TRIM(adjustl(spcnam(i)))
              numVar=numVar+1
              nameVar(numVar)='AVE_v_'//TRIM(adjustl(spcnam(i)))
              numVar=numVar+1
              nameVar(numVar)='AVE_w_'//TRIM(adjustl(spcnam(i)))
           endif
        enddo
        ! imagaT
!
        do i=1,nrnsx
        if(irans_ave(i).eq.1) then
          numVar=numVar+1
          nameVar(numVar)='AVE_'//TRIM(adjustl(sclname(i)))
        endif
!
        if(irans_rms(i).eq.1) then
          numVar=numVar+1
          nameVar(numVar)='RMS_'//TRIM(adjustl(sclname(i)))
        endif
        if(iuvw_rans_re(i).eq.1) then
          numVar=numVar+1
          nameVar(numVar)='Re_u_'//TRIM(adjustl(sclname(i)))
          numVar=numVar+1
          nameVar(numVar)='Re_v_'//TRIM(adjustl(sclname(i)))
          numVar=numVar+1
          nameVar(numVar)='Re_w_'//TRIM(adjustl(sclname(i)))
        endif
        enddo
!
        if(ical_prt==1) then
          numVar=numVar+1
          nameVar(numVar)='ave_particle_VF'
        endif
      endif
!
      if(BC_output==yes) then
         numWallVar=10
      else
         numWallVar=0
      endif
      if(numWallVar.gt.0) then
         nameWallVar(1) ='nx;wall_normal'
         nameWallVar(2) ='ny'
         nameWallVar(3) ='nz'
         nameWallVar(4) ='yplus'
         nameWallVar(5) ='Fsx;wall_shear_force'
         nameWallVar(6) ='Fsy'
         nameWallVar(7) ='Fsz'
         nameWallVar(8) ='Fpx;wall_press_force'
         nameWallVar(9) ='Fpy'
         nameWallVar(10)='Fpz'
!
!         write(ifll,'(1x,a)') ' MSG: Write Boundary Variables'
!         do nb=1,numofFForceWall
!         write(ifll,'(1x,2a)') ' MSG: Wall=',trim(nameofFForceWall(nb))
!         enddo
!
      endif
!-----------------------------------------
! --- Output result file [***.frontflow]
!-----------------------------------------
      open(ifl,file=fnam,form='unformatted',status='unknown',iostat=ios)
      write(ifl) iter,time,NMAT,NCV,NCVFAC,ncomp,
     &           nrnsx,npotnx,numVar,
     &           NFLIDx,idrdp,comp,
     &           numWallVar,nbcnd,ncompall,
     &           num_ratout,WDOT_gas,radflag,
     &           icon_cv,vertex_cen,cell_cen,ibot,uf_mm,u_func(1:uf_mm),
     &           ical_prt


      write(ifl) ical_sld,ical_suf,ieul2ph,ical_MHD,ical_WDR,N_INJ
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==3.or.ical_sld==4) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT>0) then
          IMAT_U=nofld(IMAT)
        else
          IMAT_U=nofld(IMAT)
        endif
        write(ifl) IIMAT,IMAT,IMAT_U,ishaft(IMAT),rotati(IMAT),
     &             end(:,IMAT),begin(:,IMAT),rot_ang(IMAT)
        enddo
      endif
      if(ical_suf==1) then
        write(ifl) SDOT_suf,WDOT_gas,depoeth_spd,molefr_suf,
     &  nsdot,ndeps,nmolfr,num_ratout,blk_thick,nthick
      endif
!
! --- 
!
      if(BC_output==yes) then
        do nb=1,nbcnd
          write(ifl) boundName(nb,1)
        end do
        write(ifl) (iwallvar(nb), nb=1,nbcnd)
        write(ifl) (LBC_INDEX(nb),nb=1,nbcnd)
! --- boundary CVs
        do nb=1,nbcnd
           IBFS=LBC_INDEX(nb-1)+1
           IBFE=LBC_INDEX(nb)
           if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
           write(ifl) (MAT_CV(LVEDGE(1,LBC_SSF(IBFL))),
     &                                 IBFL=IBFS,IBFE)
           endif
        enddo
      endif
!      do nb=1,nbcnd
!        IBFS=LBC_INDEX(nb-1)+1
!        IBFE=LBC_INDEX(nb)
!        do IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        IDC=LVEDGE(2,ICFL)
!        ICV=LVEDGE(1,ICFL)
!        enddo
!      enddo
!-----------------------
! --- variable names
!-----------------------
      do i=1,numVar
      write(ifl) nameVar(i)
      enddo
!
      if(numWallVar.gt.0) then
         do i=1,numWallVar
         write(ifl) nameWallVar(i)
         enddo
      end if
!

      write(ifl) (prs(LCV_CV(ICV)),ICV=1,NCV)
      
      write(ifl) (pp0(LCV_CV(ICV)),ICV=1,NCV)
      write(ifl) (rho(LCV_CV(ICV)),ICV=1,NCV)
      write(ifl) (tmp(LCV_CV(ICV)),ICV=1,NCV)
      write(ifl) (rmut(LCV_CV(ICV)),ICV=1,NCV)
!--------------------------------------
! --- using temp memory for saving
!--------------------------------------
      if(ieul2ph>0.or.
     &  ical_sld==1.or.
     &  ical_sld==2.or.
!!     &  ical_sld==3.or.
     &  ical_sld==4
     &  ) then
        allocate(vlr(MXALLCV,3),stat=ierr)
        if(ierr.ne.0) then
          write(ifle,*) 'allocating array error in write_result=>sld'
          call FFRABORT(1,'allocating array at write_result')
        endif
        if(ieul2ph>0) then
           do i=1,3
           vlr(:,i)=vel(:,i)*aks(:,iaph(1))/(aks(:,iaph(1))+SML1)
           enddo
        endif
        if(ical_sld==1.or.ical_sld==2.or.ical_sld==3.or.
     &     ical_sld==4) then
          do  IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          tha=rot_ang(IMAT)
          unit(1)=end(1,IMAT)-begin(1,IMAT)
          unit(2)=end(2,IMAT)-begin(2,IMAT)
          unit(3)=end(3,IMAT)-begin(3,IMAT)
          dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
          unit(:)=unit(:)/dum1
          CALL rotth(unit,tha,bb)

          do i=1,3
          do j=1,3
          rbb(i,j)=bb(j,i)
          enddo
          enddo
          costh=cos(tha)
          sinth=sin(tha)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
          radi(1)=CVCENT(1,ICVL)-begin(1,IMAT)
          radi(2)=CVCENT(2,ICVL)-begin(2,IMAT)
          radi(3)=CVCENT(3,ICVL)-begin(3,IMAT)
          call AXB_UNIT_C(unit,radi,vr)
          dr=radi(1)*unit(1)+radi(2)*unit(2)+radi(3)*unit(3)
          radi(1)=radi(1)-dr*unit(1)
          radi(2)=radi(2)-dr*unit(2)
          radi(3)=radi(3)-dr*unit(3)
          dr=dsqrt(radi(1)*radi(1)+radi(2)*radi(2)+radi(3)*radi(3))
!
          v0(1)=dr*rotati(IMAT)*vr(1)
          v0(2)=dr*rotati(IMAT)*vr(2)
          v0(3)=dr*rotati(IMAT)*vr(3)
!          v0(:)=0.d0  !zhang
!
          vlr(ICVL,1)=rbb(1,1)*(vel(ICVL,1)+v0(1))
     &               +rbb(1,2)*(vel(ICVL,2)+v0(2))
     &               +rbb(1,3)*(vel(ICVL,3)+v0(3))
          vlr(ICVL,2)=rbb(2,1)*(vel(ICVL,1)+v0(1))
     &               +rbb(2,2)*(vel(ICVL,2)+v0(2))
     &               +rbb(2,3)*(vel(ICVL,3)+v0(3))
          vlr(ICVL,3)=rbb(3,1)*(vel(ICVL,1)+v0(1))
     &               +rbb(3,2)*(vel(ICVL,2)+v0(2))
     &               +rbb(3,3)*(vel(ICVL,3)+v0(3))
          enddo
          enddo
          write(ifl) ((vlr(LCV_CV(ICV),i),i=1,3),ICV=1,NCV)
        else
          write(ifl) ((vlr(LCV_CV(ICV),i),i=1,3),ICV=1,NCV)
        endif
        deallocate(vlr)
      else
	write(ifl) ((vel(LCV_CV(ICV),i),i=1,3),ICV=1,NCV)
      endif 
!----------------------------
! --- species 
!----------------------------
      if(mol_rslt==1) then
        allocate(vlr(MXALLCV,ncomp),stat=ierr)
        if(ierr.ne.0) then
          write(ifle,*) 'allocating array error in write_result=>yys'
          call FFRABORT(1,'allocating array at write_result')
        endif
        do 200 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        do 210 ICVL=ICVS,ICVE
        dum1=0.d0
        do ICOM=1,ncomp
          dum1=dum1+yys(ICVL,ICOM)*r_wm(ICOM)
        enddo
        do ICOM=1,ncomp
        dum2=wm(ICOM)*dum1
        vlr(ICVL,ICOM)=yys(ICVL,ICOM)/dum2
        enddo
 210    enddo
 200    enddo
!
        write(ifl) ((vlr(LCV_CV(ICV),i),i=1,ncomp),ICV=1,NCV)
        deallocate(vlr)
      else
        write(ifl) ((yys(LCV_CV(ICV),i),i=1,ncomp),ICV=1,NCV)
!????
!        dum1=0.d0
!        do IIMAT=1,NMAT
!        IMAT=MAT_NO(IIMAT)
!        if(IMAT<0) cycle
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        do ICVL=ICVS,ICVE
!        dum1=dum1+yys(ICVL,4)*rho(ICVL)*CVVOLM(ICVL)
!        enddo
!        enddo


      endif
!
! --- 
!
      if(ibot/=0) then
        allocate(vlr(MXALLCV,1),stat=ierr)
        vlr(:,1)=0.d0
        IBFS=LBC_INDEX(ibot-1)+1
        IBFE=LBC_INDEX(ibot)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(1,ICFL)
        vlr(IDC,1)=SFCENT(3,ICFL)
        enddo
        write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
        deallocate(vlr)
      endif
!-----------------------------------------
! --- MHD 
!-----------------------------------------
      if(ical_MHD==1.or.ical_MHD==2) then
        
!        call dc_symprvMHD
!     &(2,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
!     & LCYCOLD,wifsld,OPPANG,
!     & SFAREA,SFCENT,BBB,2)
        DO IMHD=1,2
        do i=1,3   !1)
        write(ifl) (BBB(LCV_CV(ICV),i,IMHD),ICV=1,NCV)
        enddo
        do i=1,3   !2)
        write(ifl) (MHD_CRT0(LCV_CV(ICV),i,IMHD),ICV=1,NCV)
        enddo
        do i=1,3   !3)
        write(ifl) (MHD_CRNT(LCV_CV(ICV),i,IMHD),ICV=1,NCV)
        enddo
        do i=1,3   !4)
        write(ifl) (MHD_A(LCV_CV(ICV),i,IMHD),ICV=1,NCV)
        enddo
        do i=1,3   !5)
        write(ifl) (LORENZ(LCV_CV(ICV),i,IMHD),ICV=1,NCV)
        enddo
        write(ifl) (MHD_FAI(LCV_CV(ICV),IMHD),ICV=1,NCV)
        enddo
      endif
!---------------------------
! --- Save surface result
!---------------------------
      if(ical_suf==1) then
        if(SDOT_suf==1) then
          allocate(vlr(MXCV,nsdot),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'allocating array error in write_result=>suf'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          iicom=0
          do nb=1,nbcnd
!
          ami(:)=0
          do IRC=1,nneq
          if(IDX_SUFRAC(nb,IRC)) then
            iph=1
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            do  ipcom=icoms,icome
            icom=phs_com(ipcom)
            if(vreq(ICOM,2,IRC)>SML.or.vreq(ICOM,1,IRC)>SML) then
              ami(icom)=1
            endif
            enddo
          endif
          enddo
!
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1.or.surfreac(nb)==2) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==gasphase) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ami(icom)==1)then
                iicom=iicom+1
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICVL=LVEDGE(1,ICFL)
                dum2=SDOT(IBFL,icom)
                vlr(ICVL,iicom)=dum2         !dum2=mol/m^2/s
                enddo
              endif
              enddo
            elseif(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              dum2=SDOT(IBFL,icom)         !dum2=mol/m^2/s
              vlr(ICVL,iicom)=dum2         !dum2=mol/m^2/s
              enddo
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              dum2=SDOT(IBFL,icom)         !dum2=mol/m^2/s
              vlr(ICVL,iicom)=dum2         !dum2=mol/m^2/s
              enddo
              enddo            
            endif
            enddo
          endif
          enddo

          if(nsdot/=iicom) then
            call FFRABORT(1,'ERR: nsdot/=iicom')
          endif
          write(ifl) ((vlr(LCV_CV(ICV),i),i=1,iicom),ICV=1,NCV)
          deallocate(vlr)
        endif
!
        if(molefr_suf==1) then
          allocate(vlr(MXCV,nmolfr),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'allocating array error in write_result=>suf'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              vlr(ICVL,iicom)=MOLFRC(IBFL,icom)
              enddo
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              vlr(ICVL,iicom)=MOLFRC(IBFL,icom)
              enddo
              enddo       
            endif
            enddo
          endif
          enddo
          if(nmolfr/=iicom) then
            call FFRABORT(1,'ERR: nmolfr/=iicom')
          endif
          write(ifl) ((vlr(LCV_CV(ICV),i),i=1,iicom),ICV=1,NCV)
          deallocate(vlr)
        endif
!
        if(depoeth_spd==1) then
          allocate(vlr(MXCV,ndeps),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'allocating array error in write_result=>suf'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(NB)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase) then
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              iicom=iicom+1
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              dum2=SUFBND(IBFL,icom)
              vlr(ICVL,iicom)=dum2
              enddo
              enddo
            endif
            enddo
          endif
          enddo
          if(ndeps/=iicom) then
            call FFRABORT(1,'ERR: ndeps/=iicom')
          endif
          write(ifl) ((vlr(LCV_CV(ICV),i),i=1,iicom),ICV=1,NCV)
          deallocate(vlr)
        endif
!
        if(blk_thick==1) then
          allocate(vlr(MXCV,ndeps),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'allocating array error in write_result=>suf'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase) then
            elseif(ntp==blkphase.and.nbx==nb) then
              iicom=iicom+1
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              dum2=BLKTHK(IBFL,iph)
              vlr(ICVL,iicom)=dum2
              enddo
            endif
            enddo
          endif
          enddo
          if(nthick/=iicom) then
            call FFRABORT(1,'ERR: nthick/=iicom')
          endif
          write(ifl) ((vlr(LCV_CV(ICV),i),i=1,iicom),ICV=1,NCV)
          deallocate(vlr)
        endif
      endif
!
      if(ical_WDR==1) then
        allocate(vlr(MXCV,2),stat=ierr)
        if(ierr.ne.0) then
          write(ifle,*) 
     &    'allocating array error in write_result=>ical_WDR'
          call FFRABORT(1,'allocating array at write_result')
        endif
        vlr(:,:)=0.d0

        do I=1,N_inj
        do ICV=1,NCV
        dum2=Rh_d(I)*max(0.d0,time-starttime(I))
        dum1=WDRBND(ICV,I)/(dum2+SML)
        vlr(ICV,1)=dum1
        vlr(ICV,2)=vlr(ICV,2)+dum1*func_d(I)
        enddo
        write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
        enddo
        write(ifl) (vlr(LCV_CV(ICV),2),ICV=1,NCV)
        deallocate(vlr)
      endif
!
      if(WDOT_gas==1) then
          write(ifl) 
     &    ((r_wm(i)*WDOT(LCV_CV(ICV),i),i=1,ncomp),ICV=1,NCV)
      endif
!
      if(num_ratout>0) then
          allocate(vlr(MXCV,neq),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'ERR: allocate in num_ratout'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          call getfil(iflsuf,fnamsuf,'eqrate_wk')
!fortarn90          iflsuf=ifl+NPE
          if(NPE>1) then
            fnamsuf=trim(suf_wk)
!fortarn90            iflsuf=iflsuf+my_rank
          else
            fnamsuf='eqrate_wk'
          endif
          open(iflsuf,file=fnamsuf,form='unformatted',
     &       status='unknown',iostat=ios)
          if(ios/=0) then
            call FFRABORT(1,'ERR: open eqrate_wk error in result sub')
          endif
          read(iflsuf) vlr(1:NCV,1:num_ratout)
          close(iflsuf) 
          write(ifl) ((vlr(LCV_CV(ICV),amii(iicom)),iicom=1,neq),
     &                                                ICV=1,NCV)
          deallocate(vlr)
      endif
!
      if(ieul2ph>0) then 
        if(ieul2ph==1) then
          allocate(vlr(MXALLCV,1),stat=ierr)
          vlr(:,1)=rho(:)*aks(:,iaph(1))+rho2(:)*aks(:,iaph(2))
          write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
          deallocate(vlr)
        else
          write(ifl) (prs2(LCV_CV(ICV)),ICV=1,NCV)
        endif
        write(ifl) (rho2(LCV_CV(ICV)),ICV=1,NCV)
        write(ifl) (tmp2(LCV_CV(ICV)),ICV=1,NCV)
        write(ifl) (rmut2(LCV_CV(ICV)),ICV=1,NCV)
        allocate(vlr(MXALLCV,3),stat=ierr)
        if(ierr.ne.0) then
          write(ifle,*) 'allocating array error in write_result=>sld'
          call FFRABORT(1,'allocating array at write_result')
        endif
        do i=1,3
        vlr(:,i)=vel2(:,i)*aks(:,iaph(2))/(aks(:,iaph(2))+SML1)
        enddo
        write(ifl) ((vlr(LCV_CV(ICV),i),i=1,3),ICV=1,NCV)

        write(ifl) ((yys2(LCV_CV(ICV),i),i=1,ncomp),ICV=1,NCV)
        deallocate(vlr)
      endif
!
      if( nrnsx.gt.0 ) then
        if(ical_cavi==1.or.ical_vof==1) then
          allocate(vlr(NCV,1),stat=ierr)
        endif
        do IMD=1,nrnsx
        if(ical_cavi==1.and.(IMD==ivold.or.IMD==icavi)) then
          if(IMD==icavi) dum1=1.d-6
          if(IMD==ivold) dum1=1.d-4
          do ICV=1,NCV
            if(aks(ICV,IMD)<dum1) then
              vlr(ICV,1)=0.d0
            else
              vlr(ICV,1)=aks(ICV,IMD)
            endif
          enddo
          write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
!        elseif(ical_vof==1.and.IMD==ivof) then
!          do ICV=1,NCV
!          dum1=aks(ICV,IMD)
!          if(dum1-0.5d0>0.d0) then
!            vlr(ICV,1)=1.d0
!          else
!            vlr(ICV,1)=0.d0
!          endif
!          enddo
!          write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
        else
          write(ifl) (aks(LCV_CV(ICV),IMD),ICV=1,NCV)
        endif
        enddo
        if(ical_cavi==1.or.ical_vof==1) then
          deallocate(vlr)
        endif
      endif
      if(npotnx.gt.0) then
        allocate(vlr(MXALLCV,npotnx),stat=ierr)
        vlr(:,:)=POTNAL(:,:)
        do i=1,npotnx
        
!        call dc_out_potn(i,npotnx,LVEDGE,LBC_SSF,LCYCSF,POTNBC,
!     &                   wiface,mat_cal,MAT_NO,vlr)
        write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
        enddo
        deallocate(vlr)
      endif
!
      if(radflag.eq.1) then
        write(ifl) (RadHeatFlux(LCV_CV(ICV)),ICV=1,NCV)
	ALLOCATE(vlr(NALLCV,1))
	vlr=0.d0
	do nb=1,nbcnd	
	  IBFS=LBC_INDEX(nb-1)+1
	  IBFE=LBC_INDEX(nb)
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  vlr(ICV,1)=RadHeatFlux(IDC)
!RadHeatFlux(IDC) is QW  ::  radiation generated heat flux
!NOTE: positive wall heat flux means radiation_heat_loss of wall
          enddo
          if(kd.eq.kdintr.or.kd.eq.kdsld)  then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
	    ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            ICFP=LCYCSF(IBFL)
            ICVP=LVEDGE(1,ICFP)
            IDCP=LVEDGE(2,ICFP)
            vlr(ICVP,1)=RadHeatFlux(IDC)
            enddo
          endif
        enddo
        write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
        deallocate(vlr)
      endif
      
      if(NFLIDx>0) then
        
        write(ifl) ((FIELD_U(LCV_CV(ICV),i),i=1,NFLIDx),ICV=1,NCV)
      endif
!
!
      if(idrdp==comp) then
        ALLOCATE(vlr(NALLCV,2))
        do  IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        do ICVL=ICVS,ICVE
        dum1=vel(ICVL,1)**2+vel(ICVL,2)**2+vel(ICVL,3)**2
        vlr(ICVL,1)=dsqrt(dum1/CCC(ICVL))
        vlr(ICVL,2)=dsqrt(CCC(ICVL))
        enddo
        enddo
        write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
        write(ifl) (vlr(LCV_CV(ICV),2),ICV=1,NCV)
        deallocate(vlr)
      endif
!
      IF(NPE.gt.1) then
        ALLOCATE(vlr(NALLCV,1))
        vlr(:,1)=dble(my_rank)
        write(ifl) (vlr(LCV_CV(ICV),1),ICV=1,NCV)
        deallocate(vlr)
      endif
!
      if(ical_prt==1) then
        write(ifl) (PRT_VF(LCV_CV(ICV)),ICV=1,NCV)
      endif
!
      if(u_func(2)>=1.or.u_func(3)==1) then
        write(ifl) (dble(multiR(LCV_CV(ICV))),ICV=1,NCV)
      endif

!      write(ifl) (rva(ICF),ICF=1,NCVFAC)
!      close(ifl)
!----------------------
!-< 1.2. statistics  >-
!----------------------
      call getfil(ifl1,fnam1,'statis')
!fortarn90      ifl1=ifl+NPE
      if(NPE.gt.1) then
        fnam1=TRIM(statis)
!fortarn90        ifl1=ifl1+my_rank
      else
        fnam1='statis'
      endif
!
      if(.not.unstdy) then
        goto 2000
      endif
!
      if(iter.lt.ISTART) then
        inquire(file=fnam1,exist=fexist)
        if(fexist) call unlink(fnam1)
        goto 2000
      elseif(iter.eq.ISTART) then
        goto 2000
      elseif(iter.gt.ISTART) then
         if( ista_momery == 1 ) then
            inquire(file=fnam1,exist=fexist)
            if(fexist) call unlink(fnam1)
         elseif( ista_momery /= 1 ) then
            inquire(file=fnam1,exist=fexist)
            if(.not.fexist) then
               goto 2000
            endif
         endif

!
        ALLOCATE(icomp_avex(ncomp),irans_avex(nrans),stat=ierr1)
        ALLOCATE(icomp_rmsx(ncomp),irans_rmsx(nrans),stat=ierr1)
        ALLOCATE(iuvw_rans_rex(nrans),stat=ierr1)
!
        allocate(iuvwc_rex(ncomp),stat=ierr1)
!
        allocate(vlr(MXCV,3),stat=ierr)
        if(ierr.ne.0) then
          write(ifle,*) 'allocating array error in write_result'
          call FFRABORT(1,'allocating array at write_result')
        endif
!
        if( ista_momery /= 1 ) then
           open(ifl1,file=fnam1,form='unformatted',
     &          status='unknown',iostat=ios)
           iuvw_ave_rms_rex=0
           ip_avex=0;imu_avex=0;it_avex=0;imaxminx=0
           icomp_avex(:)=0;irans_avex(:)=0
           it_rmsx=0;ip_rmsx=0;imu_rmsx=0
           icomp_rmsx(:)=0;irans_rmsx(:)=0
           iuvwt_rex=0
           iuvw_rans_rex(:)=0
           ! Tagami
           iuvwc_rex(:)=0
           ! imagaT
           read(ifl1) iuvw_ave_rms_rex,ip_avex,imu_avex,
     &                it_avex,it_rmsx,ip_rmsx,imu_rmsx,
     &          iuvwt_rex,imaxminx,ista_momeryx,ical_prtx
           if(ncomp>0) then
              read(ifl1) (icomp_avex(i),i=1,ncomp),
     &             (icomp_rmsx(i),i=1,ncomp)
           endif
           ! Tagami
           if(ncomp > 0 ) then
              read(ifl1) (iuvwc_rex(i),i=1,ncomp)
           endif
           ! imagaT
           if(nrnsx>0) then
              read(ifl1) (irans_avex(i),i=1,nrnsx),
     &             (irans_rmsx(i),i=1,nrnsx),
     &             (iuvw_rans_rex(i),i=1,nrnsx)
           endif

        else
           ista_momeryx = ista_momery
           iuvw_ave_rms_rex=iuvw_ave_rms_re
           ip_avex=ip_ave;it_avex=it_ave;imaxminx=imaxmin
           imu_avex=imu_ave
           icomp_avex(:)=icomp_ave(:);irans_avex(:)=irans_ave(:)
           it_rmsx=it_rms;ip_rmsx=ip_rms;imu_rmsx=imu_rms
           icomp_rmsx(:)=icomp_rms(:);irans_rmsx(:)=irans_rms(:)
           iuvwt_rex=iuvwt_re
           iuvw_rans_rex(:)=iuvw_rans_re(:)
           ! Tagami
           iuvwc_rex(:)=iuvwc_re(:)
           ! imagaT
           ical_prtx=ical_prt
        endif

        write(ifl) iuvw_ave_rms_re,ip_ave,imu_ave,
     &             it_ave,it_rms,ip_rms,imu_rms,
     &             iuvwt_re,imaxmin,ista_momery,ical_prtx
        if(ncomp>0) then
           write(ifl) (icomp_ave(i),i=1,ncomp),
     &          (icomp_rms(i),i=1,ncomp)
        endif
        ! Tagami
        if( ncomp>0 ) then
           write(ifl) (iuvwc_re(i),i=1,ncomp)
        endif
        ! imagaT
        if(nrnsx>0) then
           write(ifl) (irans_ave(i),i=1,nrnsx),
     &          (irans_rms(i),i=1,nrnsx),
     &          (iuvw_rans_re(i),i=1,nrnsx)
        endif        

!
! --- ------------------------------------------------
!
        in_ave=0
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momeryx==1) then
            do i=1,3
            in_ave=in_ave+1
            vlr(1:NCV,i)=STA_SAV(1:NCV,in_ave)
            enddo
          else
            do i=1,3
            read(ifl1) (vlr(ICV,i),ICV=1,NCV)
            enddo
          endif
          if(ical_sld==1.or.ical_sld==2.or.
     &       ical_sld==3.or.ical_sld==4) then
            if(ierr.ne.0) then
              write(ifle,*) 
     &       'allocating array error in write_result=>sld'
              call FFRABORT(1,'allocating array at cal_statis')
            endif
            do IIMAT=1,NMAT
              IMAT=MAT_NO(IIMAT)
              tha=rot_ang(IMAT)
              unit(1)=end(1,IMAT)-begin(1,IMAT)
              unit(2)=end(2,IMAT)-begin(2,IMAT)
              unit(3)=end(3,IMAT)-begin(3,IMAT)
              dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
              unit(:)=unit(:)/dum1
              CALL rotth(unit,tha,bb)
              do ii=1,3
              do j=1,3
              rbb(ii,j)=bb(j,ii)
              enddo
              enddo
              costh=cos(tha)
              sinth=sin(tha)
              ICVS=MAT_CVEXT(IIMAT-1)+1
              ICVE=MAT_INDEX(IIMAT)
              do ICVL=ICVS,ICVE
              radi(1)=CVCENT(1,ICVL)-begin(1,IMAT)
              radi(2)=CVCENT(2,ICVL)-begin(2,IMAT)
              radi(3)=CVCENT(3,ICVL)-begin(3,IMAT)
              call AXB_UNIT_C(unit,radi,vr)
              dr=radi(1)*unit(1)+radi(2)*unit(2)+radi(3)*unit(3)
              radi(1)=radi(1)-dr*unit(1)
              radi(2)=radi(2)-dr*unit(2)
              radi(3)=radi(3)-dr*unit(3)
              dr=dsqrt(radi(1)*radi(1)+radi(2)*radi(2)+radi(3)*radi(3))
              v0(1)=dr*rotati(IMAT)*vr(1)
              v0(2)=dr*rotati(IMAT)*vr(2)
              v0(3)=dr*rotati(IMAT)*vr(3)
!              v0=0.d0
!    
              radi(1)=rbb(1,1)*(vlr(ICVL,1)+v0(1))
     &               +rbb(1,2)*(vlr(ICVL,2)+v0(2))
     &               +rbb(1,3)*(vlr(ICVL,3)+v0(3))
              radi(2)=rbb(2,1)*(vlr(ICVL,1)+v0(1))
     &               +rbb(2,2)*(vlr(ICVL,2)+v0(2))
     &               +rbb(2,3)*(vlr(ICVL,3)+v0(3))
              radi(3)=rbb(3,1)*(vlr(ICVL,1)+v0(1))
     &               +rbb(3,2)*(vlr(ICVL,2)+v0(2))
     &               +rbb(3,3)*(vlr(ICVL,3)+v0(3))
              vlr(ICVL,:)=radi(:)
              
              enddo
            enddo
          endif
        else
          do I=1,3
          vlr(:,i)=vel(:,i)
          enddo
        endif
!
        if(iuvw_ave_rms_re.eq.1) then
          do i=1,3
          write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
          enddo
        endif
!
        if(iuvw_ave_rms_rex.eq.1) then
          if(ical_sld==1.or.ical_sld==2.or.
     &       ical_sld==3.or.ical_sld==4) then
            do IIMAT=1,NMAT
              IMAT=MAT_NO(IIMAT)
              tha=rot_ang(IMAT)
              unit(1)=end(1,IMAT)-begin(1,IMAT)
              unit(2)=end(2,IMAT)-begin(2,IMAT)
              unit(3)=end(3,IMAT)-begin(3,IMAT)
              dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
              unit(:)=unit(:)/dum1
              CALL rotth(unit,tha,bb)
              do ii=1,3
              do j=1,3
              rbb(ii,j)=bb(j,ii)
              enddo
              enddo
              costh=cos(tha)
              sinth=sin(tha)
              ICVS=MAT_CVEXT(IIMAT-1)+1
              ICVE=MAT_INDEX(IIMAT)
              do ICVL=ICVS,ICVE
              radi(1)=CVCENT(1,ICVL)-begin(1,IMAT)
              radi(2)=CVCENT(2,ICVL)-begin(2,IMAT)
              radi(3)=CVCENT(3,ICVL)-begin(3,IMAT)
              call AXB_UNIT_C(unit,radi,vr)
              dr=radi(1)*unit(1)+radi(2)*unit(2)+radi(3)*unit(3)
              radi(1)=radi(1)-dr*unit(1)
              radi(2)=radi(2)-dr*unit(2)
              radi(3)=radi(3)-dr*unit(3)
              dr=dsqrt(radi(1)*radi(1)+radi(2)*radi(2)+radi(3)*radi(3))
              v0(1)=dr*rotati(IMAT)*vr(1)
              v0(2)=dr*rotati(IMAT)*vr(2)
              v0(3)=dr*rotati(IMAT)*vr(3)
              radi(1)=rbb(1,1)*(vlr(ICVL,1))
     &                   +rbb(2,1)*(vlr(ICVL,2))
     &                   +rbb(3,1)*(vlr(ICVL,3))
     &                   -v0(1)
              radi(2)=rbb(1,2)*(vlr(ICVL,1))
     &                   +rbb(2,2)*(vlr(ICVL,2))
     &                   +rbb(3,2)*(vlr(ICVL,3))
     &                   -v0(2)
              radi(3)=rbb(1,3)*(vlr(ICVL,1))
     &                   +rbb(2,3)*(vlr(ICVL,2))
     &                   +rbb(3,3)*(vlr(ICVL,3))
     &                   -v0(3)
              vlr(ICVL,:)=radi(:)
              enddo
            enddo
          endif
        endif
        do 110 i=1,3
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=SQRT(ABS(vlasum(ICV)-vlr(ICV,i)*vlr(ICV,i)))
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
 110    continue
!
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=vlasum(ICV)-vlr(ICV,1)*vlr(ICV,2)
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=vlasum(ICV)-vlr(ICV,2)*vlr(ICV,3)
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=vlasum(ICV)-vlr(ICV,1)*vlr(ICV,3)
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
! --- ------------------------------------------------
!pressure
        if(ip_avex==1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(ip_ave==1) then
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
        if(ip_rmsx.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            temp(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (temp(ICV),ICV=1,NCV)
          endif
        else
          temp(:)=0.d0
        endif
!
        if(ip_rms.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=dSQRT(abs(temp(ICV)-vlasum(ICV)*vlasum(ICV)))
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!-----------------------------
!Mu
!-----------------------------
        if(imu_avex==1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
        if(imu_ave==1) then
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
        if(imu_rmsx.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            temp(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (temp(ICV),ICV=1,NCV)
          endif
        else
          temp(:)=0.d0
        endif
!
        if(imu_rms.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=dSQRT(abs(temp(ICV)-vlasum(ICV)*vlasum(ICV)))
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
! --- ------------------------------------------------
!
        if(imaxminx.eq.1) then
          if(ista_momeryx==1) then
            do i=1,3
            in_ave=in_ave+1
            vlr(1:NCV,i)=STA_SAV(1:NCV,in_ave)
            enddo
          else
            do i=1,3
            read(ifl1) (vlr(ICV,i),ICV=1,NCV)
            enddo
          endif
        else
          do I=1,3
          vlr(:,i)=vel(:,i)
          enddo
        endif
        if(imaxmin.eq.1) then
          do i=1,3
          write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
          enddo
        endif
!
        if(imaxminx.eq.1) then
          if(ista_momeryx==1) then
            do i=1,3
            in_ave=in_ave+1
            vlr(1:NCV,i)=STA_SAV(1:NCV,in_ave)
            enddo
          else
            do i=1,3
            read(ifl1) (vlr(ICV,i),ICV=1,NCV)
            enddo
          endif
        else
          do I=1,3
          vlr(:,i)=vel(:,i)
          enddo
        endif
        if(imaxmin.eq.1) then
          do i=1,3
          write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
          enddo
        endif
!
        if(imaxminx.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(imaxmin.eq.1) then
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif

        if(imaxminx.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
        if(imaxmin.eq.1) then
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
! --- ------------------------------------------------
!
        if(it_avex==1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
        if(it_ave==1) then
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif      
!
        if(it_rmsx.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            temp(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (temp(ICV),ICV=1,NCV)
          endif
        else
          temp(:)=0.d0
        endif
        if(it_rms.eq.1) then
          do ICV=1,NCV
          temp(ICV)=
     &        dSQRT(abs(temp(ICV)-vlasum(ICV)*vlasum(ICV)))
          enddo
          write(ifl) (temp(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
        do i=1,3
        if(iuvwt_rex.eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            temp(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (temp(ICV),ICV=1,NCV)
          endif
        else
          temp(:)=0.d0
        endif
        if(iuvwt_re.eq.1) then
          do ICV=1,NCV
!          vlasum(ICV)=abs(temp(ICV)-vlasum(ICV)*vlr(ICV,i))
            vlasum(ICV)=temp(ICV)-vlasum(ICV)*vlr(ICV,i)
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
        enddo
!
! ---------------------------------------------------
        iin_ave = in_ave
        do i=1,ncomp
           if(icomp_avex(i).eq.1) then
              if(ista_momeryx==1) then
                 in_ave=in_ave+1
                 vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
              else
                 read(ifl1) (vlasum(ICV),ICV=1,NCV)
              endif
           else
              vlasum(:)=0.d0
           endif
           if(icomp_ave(i).eq.1) then
              write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
           endif
!     
           if(icomp_rmsx(i).eq.1) then
              if(ista_momeryx==1) then
                 in_ave=in_ave+1
                 temp(1:NCV)=STA_SAV(1:NCV,in_ave)
              else
                 read(ifl1) (temp(ICV),ICV=1,NCV)
              endif
           else
              temp(:)=0.d0
           endif
           if(icomp_rms(i).eq.1) then
              do ICV=1,NCV
                 temp(ICV)=dSQRT(abs(temp(ICV)-vlasum(ICV)*vlasum(ICV)))
              enddo
              write(ifl) (temp(LCV_CV(ICV)),ICV=1,NCV)
           endif
        enddo
        !---------------------------------
        ! Tagami
        do icom = 1, ncomp
           if(icomp_avex(icom)/=1) then
              write(6,*) "i can't calculate mole flux."
           else
              iin_ave = iin_ave+1
           endif
           do i = 1, 3
              if( iuvwc_rex(icom) == 1 ) then
                 if(ista_momery==1) then
                    in_ave=in_ave+1
                    temp(1:ncv)=STA_SAV(1:NCV,in_ave)
                    vlasum(1:ncv)=STA_SAV(1:ncv,iin_ave)  ! <-mole comp(icom)
                 else
                    write(6,*) "i can't calculate mole flux."
                    read(ifl1) (temp(ICV),ICV=1,NCV)
                 endif
              else
                 temp(:) = 0.d0
              endif
              if( iuvwc_re(icom) == 1 ) then
                 ! bar(u*c)-bar(u)*bar(c)
!                 do icv = 1, ncv
!                    temp(icv) = temp(icv)
!                    temp(icv)=abs(temp(icv)-vlasum(icv)*vlr(icv,i))
!                    temp(icv)=temp(icv)-vlasum(icv)*vlr(icv,i)
!                 enddo
                 write(ifl) (temp(lcv_cv(icv)),icv=1,ncv)
              endif
           enddo
        enddo
        ! imagaT
        !---------------------------------
!
! ---------------------------------------------------
!
        do i=1,nrnsx
        if(irans_avex(i).eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
        if(irans_ave(i).eq.1) then
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
        if(irans_rmsx(i).eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            temp(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (temp(ICV),ICV=1,NCV)
          endif
        else
          temp(:)=0.d0
        endif
        if(irans_rms(i).eq.1) then
          do ICV=1,NCV
          temp(ICV)=dSQRT(abs(temp(ICV)-vlasum(ICV)*vlasum(ICV)))
          enddo
          write(ifl) (temp(LCV_CV(ICV)),ICV=1,NCV)
        endif
!
        do l=1,3
        if(iuvw_rans_rex(i).eq.1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            temp(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (temp(ICV),ICV=1,NCV)
          endif
        else
          temp(:)=0.d0
        endif
        if(iuvw_rans_re(i).eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=abs(temp(ICV)-vlasum(ICV)*vlr(ICV,l))
          enddo
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif
        enddo
        enddo
!
        if(ical_prtx==1) then
          if(ista_momeryx==1) then
            in_ave=in_ave+1
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(ical_prt==1) then
          write(ifl) (vlasum(LCV_CV(ICV)),ICV=1,NCV)
        endif        
! --- ------------------------------------------------
        if(ista_momery /= 1) close(ifl1)
        deallocate(vlr)
        DEALLOCATE(icomp_avex,irans_avex,icomp_rmsx,irans_rmsx)
        DEALLOCATE(iuvw_rans_rex)
        deallocate(iuvwc_rex)
      endif
! --- ------------------------------------------------

      DEALLOCATE(ami,amii)
!
! --- 
! 
 2000 continue
!
!----------------------------
! --- BC face variable
!----------------------------
      if(BC_output==yes) then
        allocate(vlr(1:3,1:MXCVFAC))
        vlr=0.d0
!--------------------------------------
! ---   Face normal vector ------------1
!--------------------------------------
        do nb=1,nbcnd
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
            do i=1,3
              write(ifl) (-1.d0*SFAREA(i,LBC_SSF(IBFL)),
     &                                            IBFL=IBFS,IBFE)
            end do
          end if
        end do
!----------------------------------------
! ---   Wall function y+ at 1st grid ----
!----------------------------------------
        do nb=1,nbcnd
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              yp=FRSTCV(IBFL)
              rnu=rmu(ICV)/rho(ICV)
              vlr(1,ICFL)=utau(IBFL)*yp/rnu ! yplus
            end do
            write(ifl) (vlr(1,LBC_SSF(IBFL)),IBFL=IBFS,IBFE)
          end if
        end do
!----------------------------------------
! ---   Friction force at wall ----------
!----------------------------------------
        do nb=1,nbcnd
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              uwall(1:3)=vel(ICV,1:3)-vel(IDC,1:3)
              uwallnrm=uwall(1)*SFAREA(1,ICFL)
     &                +uwall(2)*SFAREA(2,ICFL)
     &                +uwall(3)*SFAREA(3,ICFL)
              uwall(1:3)=uwall(1:3)-uwallnrm*SFAREA(1:3,ICFL)
              uwallnrm=dsqrt(uwall(1)**2+uwall(2)**2+uwall(3)**2)+SML
              ydis=FRSTCV(IBFL)
              vlr(1:3,ICFL)=rmut(ICV)*(uwall(1:3)/ydis)
            end do
            do i=1,3
               write(ifl) (vlr(i,LBC_SSF(IBFL)),IBFL=IBFS,IBFE)
            end do
          end if
        end do
!----------------------------------------
! ---    press force at wall ------------
!----------------------------------------
        do nb=1,nbcnd
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              vlr(1:3,ICFL)=prs(IDC)*SFAREA(1:3,ICFL)
            end do
            do i=1,3
               write(ifl) (vlr(i,LBC_SSF(IBFL)),IBFL=IBFS,IBFE)
            end do
          end if
        end do
        deallocate(vlr)
      endif
!
      close(ifl)      
!
      return
!
 9991 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'lack of result file'
      write(ifle,*) '(write_result)'
      ierror=1
      end subroutine write_result
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      subroutine calcWallPress(MXCVFAC,MXSSFBC,MXALLCV,nbcnd,
!     &     iwallvar,LVEDGE,LBC_SSF,LBC_INDEX,prs,SFAREA,
!     &     ftmp,fpress,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!      implicit none
!
!      integer,intent(in)    :: MXCVFAC,MXSSFBC,MXALLCV
!      integer,intent(in)    :: nbcnd
!      integer,intent(in)    :: iwallvar(nbcnd)
!      integer,intent(in)    :: LVEDGE (2, MXCVFAC)
!      INTEGER,INTENT(IN)    :: LBC_SSF(   MXSSFBC)
!      INTEGER,INTENT(IN)    :: LBC_INDEX( nbcnd)
!      real*8 ,intent(in)    :: prs    (   MXALLCV)
!      real*8 ,intent(in)    :: SFAREA (4, MXCVFAC)
!      real*8 ,intent(inout) :: ftmp   (3, MXSSFBC)
!      real*8 ,intent(inout) :: fpress(1:3,nbcnd+1)
!      integer,intent(inout) :: ierror
!
!      integer :: ICVL,ICVS,ICVE,nb,kd,IBFS,IBFE,IBFL,ICFL,IDC
!      real*8  :: dum1
!
!      ierror=0
!
! --- Calculate surface integral of pressure force at walls.
! --- Data on dummy cells are used as pressure at wall.
!
!      do nb=1,nbcnd
!         IBFS=LBC_INDEX(nb-1)+1
!         IBFE=LBC_INDEX(nb)
!         if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
!         do IBFL=IBFS,IBFE
!            ICFL=LBC_SSF(IBFL)
!            IDC=LVEDGE(2,ICFL)
!            dum1=prs(IDC)
!!            dum1=prs(IDC)*SFAREA(4,ICFL)
!            ftmp(1,IBFL)=dum1*SFAREA(1,ICFL)
!            ftmp(2,IBFL)=dum1*SFAREA(2,ICFL)
!            ftmp(3,IBFL)=dum1*SFAREA(3,ICFL)
!            fpress(1,nb)=fpress(1,nb)+ftmp(1,IBFL)*SFAREA(4,ICFL)
!            fpress(2,nb)=fpress(2,nb)+ftmp(2,IBFL)*SFAREA(4,ICFL)
!            fpress(3,nb)=fpress(3,nb)+ftmp(3,IBFL)*SFAREA(4,ICFL)
!!            fpress(1,nb)=fpress(1,nb)+ftmp(1,IBFL)
!!            fpress(2,nb)=fpress(2,nb)+ftmp(2,IBFL)
!!            fpress(3,nb)=fpress(3,nb)+ftmp(3,IBFL)
!         end do
!         end if
!      end do
!
!      return
!
!      end subroutine calcWallPress
!
!!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      subroutine calcWallShear(MXCVFAC,MXSSFBC,MXALLCV,nbcnd,
!     &     iwallvar,LVEDGE,LBC_SSF,LBC_INDEX,
!     &     vel,rmut,utau,SFAREA,SFCENT,FRSTCV,
!     &     ftmp,fshear,ierror)
!!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!      use module_constant
!      implicit none
!
!      integer,intent(in)    :: MXCVFAC,MXSSFBC,MXALLCV
!      integer,intent(in)    :: nbcnd
!      integer,intent(in)    :: iwallvar(nbcnd)
!      integer,intent(in)    :: LVEDGE (2, MXCVFAC)
!      INTEGER,INTENT(IN)    :: LBC_SSF(   MXSSFBC)
!      INTEGER,INTENT(IN)    :: LBC_INDEX( nbcnd)
!      real*8 ,intent(in)    :: vel    ( MXALLCV,3)
!      real*8 ,intent(in)    :: rmut   (   MXALLCV)
!      real*8 ,intent(in)    :: utau   ( 0:MXSSFBC)
!      real*8 ,intent(in)    :: SFAREA (4, MXCVFAC)
!      real*8 ,intent(in)    :: SFCENT (3, MXCVFAC)
!      real*8 ,intent(in)    :: FRSTCV (   MXSSFBC)
!      real*8 ,intent(inout) :: ftmp   (3, MXSSFBC)
!      real*8 ,intent(inout) :: fshear(1:3,nbcnd+1)
!      integer,intent(inout) :: ierror
!
!      integer :: ICV,nb,IBFS,IBFE,IBFL,ICFL,IDC
!      real*8 ::  ux,uy,uz,up,yp,rnu,tauw
!      real*8  :: dum1
!
!      ierror=0
!
! --- Calculate surface integral of shear force at walls.
! --- Using Wall fraction velocity: utau at non-slip wall
!
!      do nb=1,nbcnd
!         IBFS=LBC_INDEX(nb-1)+1
!         IBFE=LBC_INDEX(nb)
!         if(iwallvar(nb).eq.1 .and. IBFS<=IBFE) then
!         do IBFL=IBFS,IBFE
!            ICFL=LBC_SSF(IBFL)
!            ICV=LVEDGE(1,ICFL)
!            IDC=LVEDGE(2,ICFL)
!            ux=vel(ICV,1)-vel(IDC,1)
!            uy=vel(ICV,2)-vel(IDC,2)
!            uz=vel(ICV,3)-vel(IDC,3)
!            up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
!            ux=ux-up*SFAREA(1,ICFL)
!            uy=uy-up*SFAREA(2,ICFL)
!            uz=uz-up*SFAREA(3,ICFL)
!            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
!            yp=FRSTCV(IBFL)
!            tauw=rmut(ICV)*up/yp
!            tauw=rmu(ICV)*up/yp
!
!            dum1=tauw
!            ftmp(1,IBFL)=dum1*ux/up
!            ftmp(2,IBFL)=dum1*uy/up
!            ftmp(3,IBFL)=dum1*uz/up
!            fshear(1,nb)=fshear(1,nb)+ftmp(1,IBFL)*SFAREA(4,ICFL)
!            fshear(2,nb)=fshear(2,nb)+ftmp(2,IBFL)*SFAREA(4,ICFL)
!            fshear(3,nb)=fshear(3,nb)+ftmp(3,IBFL)*SFAREA(4,ICFL)
!         end do
!         end if
!      end do
!
!      return
!
!      end subroutine calcWallShear
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_restart(
     &  iter,time,delt,ffexit,outputfile,dvdd,
     &  pp0,prs,rho,hhh,yys,vel,
     &  aks,rva,
     &  MHD_A,MHD_FAI,MHD_CRNT,
     &  MOLFRC,SITDEN,BLKTHK,
     &  rho2,tmp2,yys2,vel2,prs2,hhh2,rva2,dvdd2,
     &  cord,lacell,lvcell,
     &  POTNAL,WDRBND,
     &  ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- Moving Mesh has not been used 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_io,    only  : getfil,lenfnm,ifll,ifle
      use module_output,only  : outchk,fnmlst,none,yes,kopt
      use module_model, only  : idrdp,mach0,comp,icaltb,
     &                          sles,dles,lles,ical_MHD
      use module_flags,only   : intgvv,euleri,eulere,Adams_Bashforth
      use module_scalar,only  : rns_scl
      use module_Euler2ph,only: ieul2ph
      use module_chemreac,ONLY: ical_suf
      use module_material,only : ical_sld,rotati,ishaft,end,begin,
     &                           rot_ang,nflud
      use module_model,   only : ical_mvmsh
      use module_scalar,only   : pot_scl,ip_pefc,potname
      use module_metrix,only   : STA_SAV
      use module_les   ,only   : ISTART,n_ave,ista_momery
      use module_model,   only : ical_WDR
!
!  1. Write restart file
!
      implicit none
!
!--- [dummy arguments]
!
      integer,intent(in)  :: iter
      real*8 ,intent(in)  :: time,delt
      logical,intent(in)  :: ffexit,outputfile
      real*8 ,intent(in)  :: pp0    (   MXALLCV)
      real*8 ,intent(in)  :: prs    (   MXALLCV)
      real*8 ,intent(in)  :: rho    (   MXALLCV)
      real*8 ,intent(in)  :: hhh    (   MXALLCV)
      real*8 ,intent(in)  :: yys(       MXALLCV,MXcomp)
      real*8 ,intent(in)  :: vel    (   MXALLCV,3)
      real*8 ,intent(in)  :: aks(       MXALLCVR,MXrans)
      real*8 ,intent(in)  :: rva    (   MXCVFAC)
      real*8 ,intent(in)  :: MOLFRC (MXSSFBC_SUF,MXCOMPALL)
      real*8 ,intent(in)  :: SITDEN (MXSSFBC_SUF,MXPHASE)
      real*8 ,intent(in)  :: BLKTHK (MXSSFBC_SUF,MXPHASE)
      real*8 ,intent(in)  :: dvdd   (   MXCV,   3,2)
      real*8 ,intent(in)  :: rho2   (   MXALLCVC)
      real*8 ,intent(in)  :: tmp2   (   MXALLCV2)
      real*8 ,intent(in)  :: yys2   (   MXALLCV2,MXCOMP)
      real*8 ,intent(in)  :: vel2   (   MXALLCV2,3)
      real*8 ,intent(in)  :: rva2   (   MXCVFAC2)
      real*8 ,intent(in)  :: prs2   (   MXALLCV2)
      real*8 ,intent(in)  :: hhh2   (   MXALLCV2)
      real*8 ,intent(in)  :: dvdd2  (   MXCV2   ,3,2)
      real*8 ,intent(in)  :: MHD_A   (  MXCV_MHD,3,2)
      real*8 ,intent(in)  :: MHD_FAI (  MXCV_MHD,2)
      real*8 ,intent(in)  :: MHD_CRNT(  MXCV_MHD,3,2)
      real*8 ,intent(in)  :: cord    (3,MXVRTX_m)
      integer,intent(in)  :: lacell  (  MXCELL_m)
      integer,intent(in)  :: lvcell  (8,MXCELL_m)
      REAL*8,INTENT(IN)   :: POTNAL  (  MXALLCVP,MXPOTN)
      REAL*8,INTENT(IN)   :: WDRBND  (  MXCV_WDR,N_inj)
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(lenfnm),save :: fnam
      character(lenfnm),save :: fname
      character(8)           :: text
      integer :: out,nrnsx,npotnx
      integer :: ifl=-1,NSS
      integer :: i,j,k,l,m,n,ios,IMHD,IMVMSH=0
      integer :: ICOM,IMD,ICM,ICV,ICF,IBFL,IMAT,iph
!
! --- 
!
!
      ierror=0
      if(ical_mvmsh/=0) then
        IMVMSH=ical_mvmsh
      else
        IMVMSH=0
      endif
!
      call outchk(fnmlst(1),iter,time,delt,out)
      if(out.eq.none) return
!
      if(ifl.lt.0) then
        if(NPE.gt.1) then
          call getfil(ifl,fnam,'restart')
          fnam=RESTRTout
          if(fnam.eq.' ') call FFRABORT(1,'write_restart')
        else
          call getfil(ifl,fnam,'restart')
          if(fnam.eq.' ') call FFRABORT(1,'write_restart')
        endif
        fname=fnam
        fname=adjustl(fname)
        fnam=adjustl(fnam)
      endif
!
      if(.not.ffexit .and. out.ne.yes.and..not.outputfile) return
!
! --- 
!
      if(kopt(1).eq.1.or.kopt(1).eq.3) then
        write(text,'(i8)') iter
        text=adjustl(text)
        fnam=TRIM(TRIM(fname)//'_'//TRIM(text))
      endif
!
                  nrnsx=0
      if(rns_scl) nrnsx=nrans
!
                       npotnx=npotn
      if(.not.pot_scl) npotnx=0
      
!
      open(ifl,file=fnam,form='unformatted',status='unknown',iostat=ios)
      write(ifl) NALLCV,NCV,ncomp,nrnsx,npotnx,NCVFAC,NMAT,ncompall
      write(ifl) iter,ieul2ph,ical_suf,ical_sld,ical_MHD,time,IMVMSH,
     &           ista_momery,ical_WDR,N_inj
      if(ical_sld/=0) write(ifl) (rot_ang(IMAT),IMAT=1,nflud)
      write(ifl) (pp0(ICV),ICV=1,NALLCV)
      write(ifl) (prs(ICV),ICV=1,NALLCV)
      write(ifl) (hhh(ICV),ICV=1,NALLCV)
      write(ifl) 
     &           ((yys(ICV,ICOM),ICOM=1,ncomp),ICV=1,NALLCV)
      write(ifl) ((vel(ICV,i),i=1,3),ICV=1,NALLCV)
      if(nrnsx.gt.0) then
        write(ifl) ((aks(ICV,IMD),IMD=1,nrnsx),ICV=1,NALLCV)
      endif
      if(npotnx.gt.0) then
        write(ifl) ((POTNAL(ICV,IMD),IMD=1,npotnx),ICV=1,NALLCV)
      endif
      write(ifl) (rva(ICF),ICF=1,NCVFAC)
!
      write(ifl) (((dvdd(ICV,i,k),i=1,3),ICV=1,NCV),k=1,2)
!
      if(ieul2ph>0) then
        write(ifl) (tmp2(ICV),ICV=1,NALLCV)
        write(ifl) 
     &            ((yys2(ICV,ICOM),ICOM=1,ncomp),ICV=1,NALLCV)
        write(ifl) ((vel2(ICV,i),i=1,3),ICV=1,NALLCV)
        write(ifl) (prs2(ICV),ICV=1,NALLCV)
        write(ifl) (rva2(ICF),ICF=1,NCVFAC)
        write(ifl) (((dvdd2(ICV,i,k),i=1,3),ICV=1,NCV),k=1,2)
        write(ifl) (hhh2(ICV),ICV=1,NALLCV)
      endif
!
      if(ical_suf==1) then
        NSS=MIN(NSSFBC,MXSSFBC_SUF)
        write(ifl) ((MOLFRC(IBFL,icom),icom=1,ncompall),IBFL=1,NSS)
        write(ifl) ((SITDEN(IBFL,iph),iph=1,NPHASE),IBFL=1,NSS)
        write(ifl) ((BLKTHK(IBFL,iph),iph=1,NPHASE),IBFL=1,NSS)
      endif
!
      if(ical_WDR==1) then
        write(ifl) ((WDRBND(ICV,I),ICV=1,NCV),I=1,N_inj)
      endif
!
      if(ical_MHD==1.or.ical_MHD==2) then
        DO IMHD=1,2
        write(ifl) ((MHD_A(ICV,I,IMHD),ICV=1,NALLCV),I=1,3)
        write(ifl)  (MHD_FAI(ICV,IMHD),ICV=1,NALLCV)
        write(ifl) ((MHD_CRNT(ICV,I,IMHD),ICV=1,NALLCV),I=1,3)
        enddo
      endif
!
      if(ical_mvmsh/=0) then
        write(ifl) ((cord(I,ICV),ICV=1,NVRTX),I=1,3)
        write(ifl) (lacell(ICV),ICV=1,NCELL)
        write(ifl) ((lvcell(I,ICV),ICV=1,NCELL),I=1,8)
      endif
!
      if(ista_momery==1.and.iter>=ISTART) then
        write(ifl) ((STA_SAV(ICV,i),i=1,n_ave),ICV=1,NCV)
      endif
!
      if(my_rank.eq.root) then
        write(ifll,*)
        write(ifll,6010)
        write(ifll,6000) iter,time
        write(ifll,6010)
      endif
 6000 format(2X,'|',10X,' output restart file, iter=(',i10,
     & '), time=(',1pd12.5,'Sec.)',34X,'|')

 6010 format(2X,'|',108('*'),'|')
!
      close(ifl)
!
      return
!
 9991 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'lack of restart file'
      write(ifle,*) '(write_restart)'
      ierror=1
!
      end subroutine write_restart
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_force_history1
     &    (iter,time,delt,ffexit,
     &     LVEDGE,LBC_SSF,SFAREA,SFCENT,
     &     vel,dflx,prs,sundspd,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_io,only         : getfil,lenfnm,ifle,ifll,cntlnam
      use module_output,only     : outchk,fnmlst,none,yes,kopt,
     &                             start_force
      use module_flow_sound1,only : lsound,FX,FY,FZ,nsond,iname,isound,
     &                             names
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX
      use module_les     ,only   : ISTART
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: iter
      real*8 ,intent(in)  :: time,delt
      real*8 ,intent(in)  :: sundspd             ! Speed of sound.
      logical,intent(in)  :: ffexit
      integer,intent(in)  :: LVEDGE (2, MXCVFAC)
      integer,intent(in)  :: LBC_SSF(   MXSSFBC)
      real*8 ,intent(in)  :: SFAREA (4, MXCVFAC)
      real(8),intent(in)  :: SFCENT(3,MXCVFAC)
      real(8),intent(in) :: vel(MXALLCV,3)
      real(8),intent(in) :: dflx(MXCVFAC,3)
      real*8 ,intent(in)  :: prs    (   MXALLCV)
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,nb,kd,IBFS,IBFE,IBFL
      integer :: ICV,ICF,IBF,nbs,IIMAT,ICFL,IDC,ios
      integer :: ifl=-1,out
      logical :: firstopen=.true.,fexist=.false.,title=.false.
      CHARACTER*80 :: tempStrings
      character(lenfnm),save :: fnam
      character(lenfnm),save :: fname
      character(6)           :: text
      real*8  :: dum1
!
! --- 
!
      FX=0.d0
      FY=0.d0
      FZ=0.d0
      ierror=0
      if(.not.lsound) return
!
      if(my_rank.eq.root) then
        if(ifl.lt.0) then
          call getfil(ifl,fnam,'force')
          if(fnam.eq.' ') then
            write(ifle,*) 'Force file name is not defined in ',cntlnam
            call FFRABORT(1,'write_force_history1')
          endif
          fname=fnam
          fname=adjustl(fname)
          fnam=adjustl(fnam)
          inquire(file=fnam,exist=fexist)
          if(fexist) then
            if(kopt(4).eq.1) then
              if(iter.le.1+nint(start_force)) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            elseif(kopt(4).eq.3) then
              if(time.le.delt+start_force) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            endif
          endif
        endif
      endif
!
      call outchk(fnmlst(4),iter,time,delt,out)
      if(out.eq.none) return
      if(.not.ffexit .and. out.ne.yes) return
!
! --- 
!
      if(my_rank.eq.root) then
        if(kopt(4).eq.2.or.kopt(4).eq.4) then
          write(text,'(i6)') iter
          text=adjustl(text)
          fnam=TRIM(TRIM(fname)//'_'//TRIM(text))
        endif
      endif
!
! --- Flow 
!
      do 1000 nbs=1,nsond
      nb=iname(nbs)
      IIMAT=MAT_BCIDX(nb,1)
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      do 1100 IBFL=IBFS,IBFE
      ICFL=LBC_SSF(IBFL)
      IDC=LVEDGE(2,ICFL)
      dum1=prs(IDC)*SFAREA(4,ICFL)
      FX(nbs)=FX(nbs)+dum1*SFAREA(1,ICFL)
      FY(nbs)=FY(nbs)+dum1*SFAREA(2,ICFL)
      FZ(nbs)=FZ(nbs)+dum1*SFAREA(3,ICFL)
 1100 continue
 1000 continue
! --- 
      do 1200 nbs=1,nsond
      FX(nsond+1)=FX(nsond+1)+FX(nbs)
      FY(nsond+1)=FY(nsond+1)+FY(nbs)
      FZ(nsond+1)=FZ(nsond+1)+FZ(nbs)
 1200 continue
!
      if(NPE.gt.1) then
        do 1300 nbs=1,nsond+1
          call hpcrsum(FX(nbs))
          call hpcrsum(FY(nbs))
          call hpcrsum(FZ(nbs))
 1300   continue
      endif
!
      if(my_rank.eq.root) then
        open(ifl,file=fnam,form='formatted',
     &  status='unknown',iostat=ios,position='append')
        if(title) then
          write(ifl,5000) nsond
 5000     format(1X,'time, timestep_iter, num_wall= ',I4)
          do nbs=1,nsond
          write(ifl,*) names(nbs)(:len_trim(names(nbs))),
     &                 '   FX,   FY,   FZ'
          enddo
          write(ifl,*) 'total','   FX,   FY,   FZ'
        endif
        write(ifl,'(1pe12.5,1x,I8)') time,iter
        do 100 nbs=1,nsond+1
        write(ifl,2100) FX(nbs),FY(nbs),FZ(nbs)
 100    continue
        close(ifl)
        title=.false.
      endif
!
 2000 FORMAT(A80)
 2100 format(3(e12.5,1x))
      return
      end subroutine write_force_history1
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_force_history
     &    (iter,time,delt,ffexit,
     &     MAT_NO,LVEDGE,LBC_SSF,SFAREA,SFCENT,
     &     vel,prs,soundspd,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only         : getfil,lenfnm,ifle,ifll
      use module_output,only     : outchk,fnmlst,none,yes,kopt,
     &                             start_force
      use module_flow_sound
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX
      use module_les  ,only      : ISTART
      use module_io,only         : dflxTmpFFlag
      use module_material,only   : ical_sld,rotati,ishaft,end,begin,
     &                             rot_ang,rot_angold
!
      implicit none
!--- [dummy arguments]
!
      integer,intent(in) :: iter
      real*8 ,intent(in) :: time,delt
      logical,intent(in) :: ffexit
!     
      INTEGER,INTENT(IN) :: MAT_NO(  0:MXMAT)
      integer,intent(in) :: LVEDGE(2,MXCVFAC)
      integer,intent(in) :: LBC_SSF(MXSSFBC)
      real(8),intent(in) :: SFAREA(4,MXCVFAC)
      real(8),intent(in) :: SFCENT(3,MXCVFAC)
      real(8),intent(in) :: vel(MXALLCV,3,2)
      real(8),intent(in) :: prs(MXALLCV,2)
      real(8),intent(in) :: soundspd
      integer,intent(out) :: ierror

!
! --- [local entities]
!
      logical,save :: headerOut=.true.
      logical :: fexist=.false.
      integer,save :: ifl=-1
      integer :: returnStat
      character(lenfnm),save :: fnam
      integer :: srOut,ierr,ios
      integer :: i,j,isw,iso,nb,IBFS,IBFE,IBFL,ICFL,IDC,IMAT,IIMAT
      real*8,allocatable :: fx(:),fy(:),fz(:)
!      real*8 ::
!     &          fx(1:numofSoundWall+1),
!     &          fy(1:numofSoundWall+1),
!     &          fz(1:numofSoundWall+1)
      real*8  :: unit(3),tha,tha2,rbb(3,3),bb(3,3),walcod(3),bb2(3,3),
     &           rbb2(3,3),costh2,sinth2
      logical,save :: calSWallCenter=.true.
      logical,save :: title=.false.
      character*6  :: text
      logical      :: fileExistFlag
      real(8)      :: org_x,org_y      
!
 2502 format(3(i8,1x),3(e16.9,1x))
 2503 format(e16.9,1x,i16,1x,e16.9)
!
      ierror=0
      if(.not.soundObservFlag) return
! --- Calculation of the center of sound wall FROM HERE ----------------
      if(calSWallCenter) then
        totalSFaceArea=0.d0
        soundWallCenterPosit=0.d0
        do isw=1,numofSoundWall
          nb=soundWall_bcNo(isw)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            IDC=LVEDGE(2,ICFL)
            totalSFaceArea(isw)=totalSFaceArea(isw)+SFAREA(4,ICFL)
            soundWallCenterPosit(isw,1:3)=
     &                          soundWallCenterPosit(isw,1:3)+
     &                          SFCENT(1:3,ICFL)*SFAREA(4,ICFL)
          end do
          totalSFaceArea(numofSoundWall+1)=
     &            totalSFaceArea(numofSoundWall+1)+totalSFaceArea(isw)
          soundWallCenterPosit(numofSoundWall+1,1:3)=
     &            soundWallCenterPosit(numofSoundWall+1,1:3)+
     &                                    soundWallCenterPosit(isw,1:3)
        end do
        if(NPE.gt.1) then
          call hpcrasum(totalSFaceArea,numofSoundWall+1)
          call hpcrasum(soundWallCenterPosit(:,1),numofSoundWall+1)
          call hpcrasum(soundWallCenterPosit(:,2),numofSoundWall+1)
          call hpcrasum(soundWallCenterPosit(:,3),numofSoundWall+1)
        end if
        soundWallCenterPosit(:,1)=
     &            soundWallCenterPosit(:,1)/totalSFaceArea(:)
        soundWallCenterPosit(:,2)=
     &            soundWallCenterPosit(:,2)/totalSFaceArea(:)
        soundWallCenterPosit(:,3)=
     &            soundWallCenterPosit(:,3)/totalSFaceArea(:)
        calSWallCenter=.false.
      end if
! --- Calculation of the center of sound wall TO HERE -----------------
!
! --- Preparation of output file FROM HERE ----------------------------
!     Decide output file name and logical number. Erase nonused files.
      if(my_rank.eq.root) then
        if(ifl.lt.0) then
          call getfil(ifl,fnam,'force')
          if(fnam.eq.' ') then
            write(ifle,*) 'Force file name is not defined in fflow.ctl'
            call FFRABORT(1,'write_force_history')
          endif
          fnam=adjustl(fnam)
          inquire(file=fnam,exist=fexist)
          if(fexist) then
            if(kopt(4).eq.1) then
              if(iter.le.1+nint(start_force)) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            elseif(kopt(4).eq.3) then
              if(time.le.delt+start_force) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            endif
          endif
        endif
      end if
!     Check whether output data or not.
      call outchk(fnmlst(4),iter,time,delt,srOut)
      if(srOut.eq.none) return
      if(.not.ffexit .and. srOut.ne.yes) return

!     Add step number to file names (Not used here)
      if(my_rank.eq.root) then
        if(kopt(4).eq.2.or.kopt(4).eq.4) then
          write(text,'(i6)') iter
          text=adjustl(text)
          fnam=TRIM(TRIM(fnam)//'_'//TRIM(text))
        endif
      endif
!     Output header of fluid force data.
      if(my_rank.eq.root) then
        if(headerOut) then
          inquire(file=fnam,exist=fileExistFlag)
          if(.not.(fileExistFlag)) then
            call ffr_fh_headerout(returnStat)
            if(returnStat/=0) then
              ierror=1
              return
            end if
          end if
          headerOut=.false.
        end if
      end if
! --- Preparation of output file TO HERE -------------------------------
! --- Calculation and output of the fluid force FROM HERE --------------
      allocate(fx(1:numofSoundWall+1))
      allocate(fy(1:numofSoundWall+1))
      allocate(fz(1:numofSoundWall+1))
      fx(:)=0.d0;fy(:)=0.d0;fz(:)=0.d0
!
      if(my_rank.eq.root) then
        open(ifl,file=fnam,form='formatted',
     &  status='unknown',iostat=ios,position='append')
      end if
!     Output time, iterations and the sound speed.
      if(my_rank.eq.root) then
        write(ifl,*)'#CUSTOM_DATA'
        write(ifl,*)'FFR_FORCE_FIELD'
        write(ifl,*)'FFR Force History data field'
        write(ifl,*)0
        write(ifl,*)'#FFR_FORCE_DATA'
        write(ifl,2503)time,iter,soundspd
      end if

      do iso=1,numofSoundObserver
! --- Calculate fluid force.
        if(soundModelCode(iso)==1) then
! --- Using Curle's equation.
          call calcForcebyCurle(iso)
        elseif(soundModelCode(iso)==2) then
! --- Using FWH equation.
          call calcForcebyFWH(iso)
        else
          write(ifle,*) 
     &      'Warning:Sound model is not defined for observer',iso
          write(ifle,*) 
     &      'Fluid force is not calculated for this observer'
          fx(:)=0.d0;fy(:)=0.d0;fz(:)=0.d0
        end if

        fx(numofSoundWall+1)=0.d0
        fy(numofSoundWall+1)=0.d0
        fz(numofSoundWall+1)=0.d0
        
        do isw=1,numofSoundWall
          fx(numOfSoundWall+1)=fx(numOfSoundWall+1)+fx(isw)
          fy(numOfSoundWall+1)=fy(numOfSoundWall+1)+fy(isw)
          fz(numOfSoundWall+1)=fz(numOfSoundWall+1)+fz(isw)
        end do
        if(NPE.gt.1) then
          call hpcrasum_0(fx,numofSoundWall+1,0)
          call hpcrasum_0(fy,numofSoundWall+1,0)
          call hpcrasum_0(fz,numofSoundWall+1,0)
        endif
!       Output fluid force.
        if(my_rank.eq.root) then
          write(ifl,2502)0,iso,soundModelCode(iso),fx(numOfSoundWall+1),
     &                        fy(numOfSoundWall+1),fz(numOfSoundWall+1)
          do isw=1,numofSoundWall
            write(ifl,2502)isw,iso,soundModelCode(iso),
     &                             fx(isw),fy(isw),fz(isw)
          end do
        end if
!        end do
      end do
!     Close file.
      if(my_rank.eq.root) then
        write(ifl,*)'#FFR_FORCE_DATA_END'
        close(ifl)
      end if
! --- Calculation and output of the fluid force TO HERE ----------------
      deallocate(fx,fy,fz)
      ierror=0
      return
      
      contains
      
!=======================================================================
      subroutine calcForcebyCurle(obId)
!=======================================================================
! --- Calculating fluid force by Curle's equation FROM HERE ------------
!=======================================================================
      integer,intent(in) :: obId
      real(8) :: dum1
      character(lenfnm) :: fnam1,fnam2
!
      fx(1:numofSoundWall)=0.d0
      fy(1:numofSoundWall)=0.d0
      fz(1:numofSoundWall)=0.d0
!
! --- Read friction force data
!      allocate(dflxNow(MXCVFAC,3))
!      if(NPE.gt.1) then
!        fnam1=adjustl(dflxTmpFname1)
!        fnam2=adjustl(dflxTmpFname2)
!      else
!        fnam1='dflx_tmp1.dat'
!        fnam2='dflx_tmp2.dat'
!      endif
!      if(dflxTmpFFlag) then
!        open(71,file=fnam2,form='unformatted',
!     &            action='read',status='UNKNOWN')
!      else
!        open(71,file=fnam1,form='unformatted',
!     &            action='read',status='UNKNOWN')
!      endif
!      read(71)dflxNow
!      close(71)
! 
      do isw=1,numofSoundWall
        nb=soundWall_bcNo(isw)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
! 
! --- Calculate surface integral of pressure at sound source walls.
! --- Data on dummy cells are used as pressure at wall.
! 
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          dum1=prs(IDC,1)*SFAREA(4,ICFL)
          fx(isw)=fx(isw)+dum1*SFAREA(1,ICFL)
          fy(isw)=fy(isw)+dum1*SFAREA(2,ICFL)
          fz(isw)=fz(isw)+dum1*SFAREA(3,ICFL)
        enddo
      enddo
! 
      end subroutine calcForcebyCurle
! 
!=======================================================================
      subroutine calcForcebyFWH(obId) 
!=======================================================================
! --- Calculating fluid force by FWH equation FROM HERE ----------------
!=======================================================================
      integer,intent(in) :: obId
      real(8) :: dum1,dum2,M_r,M_rprv,r_i(1:3),
     &           rabs,dumMr,Mrinv,costh,sinth,costh2,sinth2,
     &         v0(3,2),r(3),dr,vr(3),r_i2(1:3),r2(3),walcod2(3),
     &         duma(3),n(3),n2(3)
!      character(lenfnm) :: fnam1,fnam2
!
      fx(1:numofSoundWall)=0.d0
      fy(1:numofSoundWall)=0.d0
      fz(1:numofSoundWall)=0.d0
!
      do isw=1,numofSoundWall
        nb=soundWall_bcNo(isw)
        IIMAT=MAT_BCIDX(nb,1)
        IMAT=MAT_NO(IIMAT)
        tha=rot_ang(IMAT)
        unit(1)=end(1,IMAT)-begin(1,IMAT)
        unit(2)=end(2,IMAT)-begin(2,IMAT)
        unit(3)=end(3,IMAT)-begin(3,IMAT)
        dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
        unit(:)=unit(:)/dum1
        CALL rotth(unit,tha,bb)
        tha2=rot_angold(IMAT)
        CALL rotth(unit,tha2,bb2)
        do i=1,3
        do j=1,3
        rbb(i,j)=bb(j,i)
        enddo
        enddo
        do i=1,3
        do j=1,3
        rbb2(i,j)=bb2(j,i)
        enddo
        enddo
        costh=cos(tha)
        sinth=sin(tha)
        costh2=cos(tha2)
        sinth2=sin(tha2)
        org_x=begin(1,IMAT)*costh-begin(2,IMAT)*sinth
        org_y=begin(1,IMAT)*sinth+begin(2,IMAT)*costh
!
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
!------------------------------------------------------------
! --- Calculate direction vector from a sound source wall
!     to an observer point.
!------------------------------------------------------------
!
          n(1)=SFCENT(1,ICFL)-begin(1,IMAT)
          n(2)=SFCENT(2,ICFL)-begin(2,IMAT)
          n(3)=SFCENT(3,ICFL)-begin(3,IMAT)
!!!!!!          call AXB_UNIT_C(unit,n,vr)
!          
          vr(3)=(unit(1)*n(2)-unit(2)*n(1))
          vr(2)=(unit(3)*n(1)-unit(1)*n(3))
          vr(1)=(unit(2)*n(3)-unit(3)*n(2))
          dr=dsqrt(vr(1)*vr(1)+vr(2)*vr(2)+vr(3)*vr(3))
          vr(:)=vr(:)/dr
!
          dr=n(1)*unit(1)
     &      +n(2)*unit(2)
     &      +n(3)*unit(3)
          n(1)=n(1)-dr*unit(1)
          n(2)=n(2)-dr*unit(2)
          n(3)=n(3)-dr*unit(3)
!
          dr=dsqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
!
          v0(1,1)=dr*rotati(IMAT)*vr(1)+vel(IDC,1,1)
          v0(2,1)=dr*rotati(IMAT)*vr(2)+vel(IDC,2,1)
          v0(3,1)=dr*rotati(IMAT)*vr(3)+vel(IDC,3,1)
!
          duma(1)=rbb(1,1)*v0(1,1)
     &           +rbb(1,2)*v0(2,1)
     &           +rbb(1,3)*v0(3,1)
!
          duma(2)=rbb(2,1)*v0(1,1)
     &           +rbb(2,2)*v0(2,1)
     &           +rbb(2,3)*v0(3,1)
!
          duma(3)=rbb(3,1)*v0(1,1)
     &           +rbb(3,2)*v0(2,1)
     &           +rbb(3,3)*v0(3,1)
!
          v0(1,1)=duma(1)
          v0(2,1)=duma(2)
          v0(3,1)=duma(3)
!
          v0(1,2)=dr*rotati(IMAT)*vr(1)+vel(IDC,1,2)
          v0(2,2)=dr*rotati(IMAT)*vr(2)+vel(IDC,2,2)
          v0(3,2)=dr*rotati(IMAT)*vr(3)+vel(IDC,3,2)
!
          duma(1)=rbb2(1,1)*v0(1,2)
     &           +rbb2(1,2)*v0(2,2)
     &           +rbb2(1,3)*v0(3,2)
!
          duma(2)=rbb2(2,1)*v0(1,2)
     &           +rbb2(2,2)*v0(2,2)
     &           +rbb2(2,3)*v0(3,2)
!
          duma(3)=rbb2(3,1)*v0(1,2)
     &           +rbb2(3,2)*v0(2,2)
     &           +rbb2(3,3)*v0(3,2)
!
          v0(1,2)=duma(1)
          v0(2,2)=duma(2)
          v0(3,2)=duma(3)
!
          n(1)=rbb(1,1)*(SFAREA(1,ICFL))
     &        +rbb(1,2)*(SFAREA(2,ICFL))
     &        +rbb(1,3)*(SFAREA(3,ICFL))
!
          n(2)=rbb(2,1)*(SFAREA(1,ICFL))
     &        +rbb(2,2)*(SFAREA(2,ICFL))
     &        +rbb(2,3)*(SFAREA(3,ICFL))
!
          n(3)=rbb(3,1)*(SFAREA(1,ICFL))
     &        +rbb(3,2)*(SFAREA(2,ICFL))
     &        +rbb(3,3)*(SFAREA(3,ICFL))
!
          n2(1)=rbb2(1,1)*(SFAREA(1,ICFL))
     &         +rbb2(1,2)*(SFAREA(2,ICFL))
     &         +rbb2(1,3)*(SFAREA(3,ICFL))
!
          n2(2)=rbb2(2,1)*(SFAREA(1,ICFL))
     &         +rbb2(2,2)*(SFAREA(2,ICFL))
     &         +rbb2(2,3)*(SFAREA(3,ICFL))
!
          n2(3)=rbb2(3,1)*(SFAREA(1,ICFL))
     &         +rbb2(3,2)*(SFAREA(2,ICFL))
     &         +rbb2(3,3)*(SFAREA(3,ICFL))
!
          walcod(1)=SFCENT(1,ICFL)*costh
     &             -SFCENT(2,ICFL)*sinth-org_x
          walcod(2)=SFCENT(1,ICFL)*sinth
     &             +SFCENT(2,ICFL)*costh-org_y
          walcod(3)=SFCENT(3,ICFL)
!
          walcod2(1)=SFCENT(1,ICFL)*costh2
     &              -SFCENT(2,ICFL)*sinth2-org_x
          walcod2(2)=SFCENT(1,ICFL)*sinth2
     &              +SFCENT(2,ICFL)*costh2-org_y
          walcod2(3)=SFCENT(3,ICFL)
!
          r_i(1:3)=soundObservPosit(obId,1:3)-walcod(1:3)
          r_i2(1:3)=soundObservPosit(obId,1:3)-walcod2(1:3)
          rabs=dsqrt(dot_product(r_i,r_i))+SML
         
!------------------------------------------------------------
! --- Calculate the relative Mach number of sound source.
!------------------------------------------------------------
          M_r=dot_product(r_i(1:3),v0(1:3,1))      !vel(IDC,1:3,1))
     &                                      /(SML+rabs*soundspd)
          M_rprv=dot_product(r_i2(1:3),v0(1:3,2))   !vel(IDC,1:3,2))
     &                                      /(SML+rabs*soundspd)
!------------------------------------------------------------
! --- Calculate fluid force.
! --- calculated because the previous time step data is not 
! --- available. Output data is setted to ZERO.
!------------------------------------------------------------
          dum2=1.d0/(SML+(rabs**2)*((1.d0-M_r)**2))
          Mrinv=(1.d0/(1.d0-M_r))
!          dumMr=(1.d0/((1.d0-M_r)**2))*(M_r-M_rprv)/delt
          dumMr=(M_r-M_rprv)/delt
!------------------------------------------------------------
! --- Calculate surface integral of FWH equation. 
!------------------------------------------------------------
!
          dum1=(r_i(1)*dum2)*(
     &    ((n(1)*prs(IDC,1)-n2(1)*prs(IDC,2)))/delt
     &     +n(1)*Mrinv*dumMr*(prs(IDC,1)))
          fx(isw)=fx(isw)+dum1*SFAREA(4,ICFL)
!
          dum1=(r_i(2)*dum2)*(
     &    ((n(2)*prs(IDC,1)-n2(2)*prs(IDC,2)))/delt
     &     +n(2)*Mrinv*dumMr*(prs(IDC,1)))
          fy(isw)=fy(isw)+dum1*SFAREA(4,ICFL)
!
          dum1=(r_i(3)*dum2)*(
     &    ((n(3)*prs(IDC,1)-n2(3)*prs(IDC,2)))/delt
     &     +n(3)*Mrinv*dumMr*(prs(IDC,1)))
          fz(isw)=fz(isw)+dum1*SFAREA(4,ICFL)
!
        enddo
      enddo

      end subroutine calcForcebyFWH
!
!=======================================================================
      subroutine ffr_fh_headerout(statRtrn)      
!=======================================================================
!
      integer,intent(out) :: statRtrn

 2501 format(i8,1x,3(e16.9,1x))
!     Output format for the position of an observer point.
!
      open(ifl,file=fnam,form='formatted',
     &  status='unknown',iostat=ios)
      
      write(ifl,'(a)')'#A_GF_V2'
      write(ifl,'(a)')'#NEW_SET'
      write(ifl,'(a)')'FrontFlow Red Force History data'
      write(ifl,'(a)')'#CUSTOM_DATA'
      write(ifl,'(a)')'FFR_FORCE_FIELD'
      write(ifl,'(a)')'FFR Force History Header field'
      write(ifl,*)0
      write(ifl,'(a)')'#FFR_FORCE_HEADER'
      write(ifl,*)numofSoundWall+1
      write(ifl,*)numofSoundObserver
      write(ifl,*)2
      write(ifl,'(a)')'#FFR_FORCE_HEADER_END'
      write(ifl,'(a)')'#CUSTOM_DATA'
      write(ifl,'(a)')'FFR_FORCE_FIELD'
      write(ifl,'(a)')'FFR Force History Sound Model field'
      write(ifl,*)0
      write(ifl,'(a)')'#FFR_FORCE_MODEL'
      write(ifl,'(i8,1x,a)')1,'Curle'
      write(ifl,'(i8,1x,a)')2,'FWH'
      write(ifl,'(a)')'#FFR_FORCE_MODEL_END'
      write(ifl,'(a)')'#CUSTOM_DATA'
      write(ifl,'(a)')'FFR_FORCE_FIELD'
      write(ifl,'(a)')'FFR Force History Wall Label field'
      write(ifl,*)0
      write(ifl,'(a)')'#FFR_FORCE_WALLS'
      write(ifl,'(i8,1x,a)')0,'total'
      do iso=1,numofSoundWall
        write(ifl,'(i8,1x,a)')iso,nameofSoundWall(iso)
      end do
      write(ifl,2501)0,soundWallCenterPosit(numofSoundWall+1,1:3)
      do isw=1,numofSoundWall
        write(ifl,2501) isw,soundWallCenterPosit(isw,1:3)
      end do
      write(ifl,'(a)')'#FFR_FORCE_WALLS_END'
      write(ifl,'(a)')'#CUSTOM_DATA'
      write(ifl,'(a)')'FFR_FORCE_FIELD'
      write(ifl,'(a)')'FFR Force History Observer field'
      write(ifl,*)0
      write(ifl,'(a)')'#FFR_FORCE_OBSERVER'
      do iso=1,numofSoundObserver
        write(ifl,2501) iso,soundObservPosit(iso,:)
      end do
      write(ifl,'(a)')'#FFR_FORCE_OBSERVER_END'

      close(ifl)
      statRtrn=0
      end subroutine ffr_fh_headerout
      end subroutine write_force_history
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_anim
     &           (iter,time,deltt,ffexit,LCV_CV,LVEDGE,LBC_SSF,
     &                pp0,prs,rho,tmp,yys,vel,aks,
     &                MAT_NO,MAT_CVEXT,MAT_INDEX,CVCENT,
     &                SUFBND,SDOT,MOLFRC,WDOT,
     &                ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_model,only   : icaltb
      use module_io,only      : getfil,lenfnm,ifle,ifll,cntlnam
      use module_output,only  : outchk,fnmlst,none,yes,kopt,tps,
     &                          start_force
      use module_anim  ,only  : ianim,ianim_uvw,ianim_p,
     &                          ianim_r,ianim_t,ianim_rans,ianim_comp,
!
     6                          ianim_GAS_WDOT,
     &                          ianim_SDOT_suf,
     &                          ianim_molefr_suf,
     &                          ianim_depoeth_spd
!
      use module_scalar,only  : rns_scl
      use module_material,only: ical_sld,rotati,ishaft,end,begin,
     &                          rot_ang
      use module_metrix,only  : vlr =>d2work1
      use module_les   ,only  : mol_rslt
      use module_species,only : spcnam,wm,r_wm
      use module_boundary,only: nbcnd,kdbcnd,LBC_INDEX,kdcvd,
     &                          kxnone,kdbuff,boundName,
     &                          phs_idx,phs_dns,IDX_SUFRAC,
     &                          phs_nbc,phs_typ,phs_snm,phs_nam,
     &                          gasphase,surphase,blkphase,phs_com,
     &                          phs_inifrac,num_site,blk_dns,surfreac
      use module_chemreac,ONLY : ical_suf,nneq,vreq
      use module_scalar,only   : pot_scl,ip_pefc,potname
!
      implicit none
!
!
! --- [dummy arguments]
!
      integer,intent(in)  :: iter
      real*8 ,intent(in)  :: time,deltt
      logical,intent(in)  :: ffexit
      integer,intent(in)  :: LCV_CV (   MXALLCV)
      INTEGER,INTENT(IN)  :: LBC_SSF(   MXSSFBC)
      integer,intent(in)  :: LVEDGE (2, MXCVFAC)
      real*8 ,intent(in)  :: pp0    (   MXALLCV)
      real*8 ,intent(in)  :: prs    (   MXALLCV)
      real*8 ,intent(in)  :: rho    (   MXALLCV)
      real*8 ,intent(in)  :: tmp    (   MXALLCV)
      real*8 ,intent(in)  :: yys(       MXALLCV,MXcomp)
      real*8 ,intent(in)  :: vel    (   MXALLCV,3)
      real*8 ,intent(in)  :: aks(       MXALLCVR,MXrans)
      INTEGER,INTENT(IN)  :: MAT_NO (   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_INDEX( 0:MXMAT)
      real*8 ,intent(in)  :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: SUFBND (MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(IN)    :: SDOT   (MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(IN)    :: MOLFRC (MXSSFBC_SUF,MXCOMPALL)
      real*8 ,intent(inout)    :: WDOT   (   MXALLCV,MXCOMP)
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer,allocatable :: ami(:)
      character(lenfnm),save :: fnam
      character(lenfnm),save :: fname
      integer,save :: nsdot=0,ndeps=0,
     &                nmolfr=0
      character(3)           :: text
      character(80) :: nameVar(200)=' '
      integer :: out,nrnsx,numVar,npotnx
      logical :: fexist=.false.
      integer :: ifl=-1,ierr=0,ierr1=0
      integer :: i,j,k,l,m,n,ios,nb,iicom,IRC,iph,ntp,nbx,ipcom,
     &           icoms,icome,kd,IBFS,IBFE
      integer :: ICOM,IMD,ICM,ICV,ICF,IIMAT,IMAT,ICVS,ICVE,
     &           ICVL,IBFL,ICFL
      real*8  :: tha,unit(3),rbb(3,3),costh,sinth,radi(3),dr,
     &           v0(3),vr(3),bb(3,3),dum1,dum2
!
      if(ianim.eq.0) return
!
      if(ifl.lt.0) then
        if(NPE.gt.1) then
          call getfil(ifl,fnam,'anim')
          fnam=animfile
!!!          ifl=ifl+my_rank
        else
          call getfil(ifl,fnam,'anim')
        endif
        if(fnam.eq.' ') call FFRABORT(1,
     &         'Set &file/anim ')
        fname=fnam
        fname=adjustl(fname)
        fnam=adjustl(fnam)
        inquire(file=fnam,exist=fexist)
        if(fexist) then
          if(kopt(5).eq.1) then
            if(iter.le.1+nint(start_force)) then
              call unlink(fnam)
            endif
          elseif(kopt(5).eq.3) then
            if(time.le.deltt+start_force) then
              call unlink(fnam)
            endif
          endif
        endif
        if(kopt(5).eq.2.or.kopt(5).eq.4) then
          write(*,*) ' ### ERR: type=[inter_i] must be satisfied for ',
     &               'anim &output in ',cntlnam
          call FFRABORT(1,'write_anim') 
        endif
      endif
!
      ierror=0
!

!
      call outchk(fnmlst(5),iter,time,deltt,out)
      if(out.eq.none) return
!
      if(.not.ffexit .and. out.ne.yes) return
!
                       nrnsx=nrans
      if(.not.rns_scl) nrnsx=0
                       npotnx=npotn
      if(.not.pot_scl) npotnx=0
!
! --- 
!
      numVar=0
      if(ianim_uvw.eq.1) then
        numVar=numVar+3
        nameVar(1)='velo_u;Velocity'
        nameVar(2)='velo_v'
        nameVar(3)='velo_w'
      endif
!
      if(ianim_p.eq.1) then
        numVar=numVar+1
        nameVar(numVar)='pressure'
      endif
      if(ianim_r.eq.1) then
        numVar=numVar+1
        nameVar(numVar)='density'
      endif
!
      if(ianim_t.eq.1) then
        numVar=numVar+1
        nameVar(numVar)='Temperature'
      endif
!    
      do i=1,nrnsx
      if(ianim_rans(i).eq.1) then
        numVar=numVar+1
        write(text,'(i2.2)') i
        nameVar(numVar)='scl_'//TRIM(adjustl(text))
      endif
      enddo
!
!      do i=1,npotnx
!        numVar=numVar+1
!        nameVar(numVar)=TRIM(adjustl(potname(i)))
!      enddo      
!
      if(mol_rslt==1) then
        do i=1,ncomp
        if(ianim_comp(i).eq.1) then
          numVar=numVar+1
          write(text,'(i3.3)') i
          nameVar(numVar)='GasMolFr_'//TRIM(adjustl(spcnam(i)))
        endif
        enddo
      else
        do i=1,ncomp
        if(ianim_comp(i).eq.1) then
          numVar=numVar+1
          write(text,'(i3.3)') i
          nameVar(numVar)='GasMasFr_'//TRIM(adjustl(spcnam(i)))
        endif
        enddo
      endif
!
!-----------------------------------------
! --- surface module output
!-----------------------------------------
      if(ical_suf==1) then
        ALLOCATE(ami(ncomp),stat=ierr1)
!        if(ianim_SDOT_suf==1) then
!----------------------------------------
          iicom=0
          do nb=1,nbcnd
!
          ami(:)=0
          do IRC=1,nneq
          if(IDX_SUFRAC(nb,IRC)) then
            iph=1
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            do  ipcom=icoms,icome
            icom=phs_com(ipcom)
            if(vreq(ICOM,2,IRC)>SML.or.vreq(ICOM,1,IRC)>SML) then
              ami(icom)=1
            endif
            enddo
          endif
          enddo
!    
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1.or.surfreac(nb)==2) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==gasphase) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ami(icom)==1.and.ianim_SDOT_suf(icom)==1)then
                iicom=iicom+1
                numVar=numVar+1
                nameVar(numVar)=
     & 'SurMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/m^2/s]_'//
     &    trim(boundName(nb,1))
              elseif(ianim_SDOT_suf(icom)==1) then
                write(ifll,'(2a)') 
     &          'WRN: animation NOT output for',trim(spcnam(icom))
              endif
              enddo
            elseif(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_SDOT_suf(icom)==1) then
                iicom=iicom+1
                numVar=numVar+1
                nameVar(numVar)=
     & 'SurMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/m^2/s]_'//
     &    trim(phs_nam(iph))
              endif
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_SDOT_suf(icom)==1) then
                iicom=iicom+1
                numVar=numVar+1
                nameVar(numVar)=
     & 'SurMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/m^2/s]_'//
     &    trim(phs_nam(iph))
              endif
              enddo
            endif
            enddo
          endif
          enddo
          nsdot=iicom
!----------------------------------------
!        endif
!
!        if(ianim_molefr_suf==1) then
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_molefr_suf(icom)==1) then
                iicom=iicom+1
                numVar=numVar+1
                nameVar(numVar)=
     &   'SiteMolFr_'//TRIM(adjustl(spcnam(icom)))//'_'//
     &    trim(phs_nam(iph))
              endif
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_molefr_suf(icom)==1) then
                iicom=iicom+1
                numVar=numVar+1
                nameVar(numVar)=
     &    'BulkMolFr_'//TRIM(adjustl(spcnam(icom)))//'_'//
     &    trim(phs_nam(iph))
              endif
              enddo
            endif
            enddo
          endif
          enddo
          nmolfr=iicom
!        endif
!
!        if(ianim_depoeth_spd==1) then
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_depoeth_spd(icom)==1) then
                write(ifll,'(a)') 
     &        'WRN: Surface Phase species NOT depos_etch_speed'
              endif
              enddo
!              iicom=iicom+1
!              numVar=numVar+1
!              nameVar(numVar)=
!     &        TRIM(adjustl(spcnam(icom)))//'_[m/s]'
!              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_depoeth_spd(icom)==1) then
              iicom=iicom+1
                numVar=numVar+1
                nameVar(numVar)=
     &   'BulkEH/DP_'//TRIM(adjustl(spcnam(icom)))//'_[m/s]_'//
     &    trim(phs_nam(iph))
              endif
              enddo
            endif
            enddo
          endif
          enddo
!        endif
        ndeps=iicom
!
        do icom=1,ncomp
        if(ianim_GAS_WDOT(icom)==1) then
          numVar=numVar+1
          nameVar(numVar)=
     & 'GasMolRate_'//TRIM(adjustl(spcnam(icom)))//'_[mole/(s*m^3)]'
        endif
        enddo
!
      endif
!
      open(ifl,file=fnam,form='unformatted',status='unknown',
     &     iostat=ios,position='append')
      write(ifl) iter,time,NMAT,NCV,NCVFAC,ncomp,nrnsx,numVar,
     &           ncompall
      write(ifl) ical_sld,ical_suf
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==3.or.ical_sld==4) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        write(ifl) IIMAT,IMAT,ishaft(IMAT),rotati(IMAT),
     &             end(:,IMAT),begin(:,IMAT),rot_ang(IMAT)
        enddo
      endif
      if(ical_suf==1) then
        write(ifl) ianim_GAS_WDOT(1:ncomp),
     &             ianim_SDOT_suf(1:ncomp+ncomp_suf),
     &             ianim_molefr_suf(1:ncomp+ncomp_suf),
     &             nsdot,ndeps,nmolfr
      endif
!
! --- 
!
      do i=1,numVar
      write(ifl) nameVar(i)
      enddo
      write(ifl) ianim,ianim_uvw,ianim_p,ianim_r,ianim_t
      if(nrnsx>0) then
         write(ifl) (ianim_rans(i),i=1,nrnsx)
      end if
      if(ncomp>0) then
         write(ifl) (ianim_comp(i),i=1,ncomp)
      end if

!      write(ifl) ianim,ianim_uvw,ianim_p,
!     &           ianim_r,ianim_t,(ianim_rans(i),i=1,nrnsx),
!     &           (ianim_comp(i),i=1,ncomp)
      if(ianim_uvw==1) then
!        if(ical_sld==1.or.ical_sld==2.or.ical_sld==3) then
        if(ical_sld==1.or.ical_sld==2) then
           allocate(vlr(MXALLCV,3),stat=ierr)
           do  IIMAT=1,NMAT
           IMAT=MAT_NO(IIMAT)
!           if(ishaft(IMAT)==1) then
             tha=rot_ang(IMAT)
             unit(1)=end(1,IMAT)-begin(1,IMAT)
             unit(2)=end(2,IMAT)-begin(2,IMAT)
             unit(3)=end(3,IMAT)-begin(3,IMAT)
             dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
             unit(:)=unit(:)/dum1
             CALL rotth(unit,tha,bb)
             do i=1,3
             do j=1,3
             rbb(i,j)=bb(j,i)
             enddo
             enddo
             costh=cos(tha)
             sinth=sin(tha)
             ICVS=MAT_CVEXT(IIMAT-1)+1
             ICVE=MAT_INDEX(IIMAT)
             do ICVL=ICVS,ICVE
             radi(1)=CVCENT(1,ICVL)-begin(1,IMAT)
             radi(2)=CVCENT(2,ICVL)-begin(2,IMAT)
             radi(3)=CVCENT(3,ICVL)-begin(3,IMAT)
             call AXB_UNIT_C(unit,radi,vr)
             dr=radi(1)*unit(1)+radi(2)*unit(2)+radi(3)*unit(3)
             radi(1)=radi(1)-dr*unit(1)
             radi(2)=radi(2)-dr*unit(2)
             radi(3)=radi(3)-dr*unit(3)
             dr=dsqrt(radi(1)*radi(1)+radi(2)*radi(2)+radi(3)*radi(3))
             v0(1)=dr*rotati(IMAT)*vr(1)
             v0(2)=dr*rotati(IMAT)*vr(2)
             v0(3)=dr*rotati(IMAT)*vr(3)
             vlr(ICVL,1)=rbb(1,1)*(vel(ICVL,1)+v0(1))
     &                  +rbb(1,2)*(vel(ICVL,2)+v0(2))
     &                  +rbb(1,3)*(vel(ICVL,3)+v0(3))
             vlr(ICVL,2)=rbb(2,1)*(vel(ICVL,1)+v0(1))
     &                  +rbb(2,2)*(vel(ICVL,2)+v0(2))
     &                  +rbb(2,3)*(vel(ICVL,3)+v0(3))
             vlr(ICVL,3)=rbb(3,1)*(vel(ICVL,1)+v0(1))
     &                  +rbb(3,2)*(vel(ICVL,2)+v0(2))
     &                  +rbb(3,3)*(vel(ICVL,3)+v0(3))
             enddo
!          else
!             ICVS=MAT_CVEXT(IIMAT-1)+1
!             ICVE=MAT_INDEX(IIMAT)
!             do i=1,3
!             vlr(ICVS:ICVE,i)=vel(ICVS:ICVE,i)
!             enddo
!          endif
          enddo
          write(ifl) ((vlr(LCV_CV(ICV),i),i=1,3),ICV=1,NCV)
          deallocate(vlr)
        else
          write(ifl) ((vel(LCV_CV(ICV),i),i=1,3),ICV=1,NCV)
        endif
      endif
      if(ianim_p==1) write(ifl) (prs(LCV_CV(ICV)),ICV=1,NCV)
      if(ianim_r==1) write(ifl) (rho(LCV_CV(ICV)),ICV=1,NCV)
      if(ianim_t==1) write(ifl) (tmp(LCV_CV(ICV)),ICV=1,NCV)
!
      if(nrnsx>0) then
      do i=1,nrnsx
      if(ianim_rans(i)==1) then
        write(ifl) (aks(LCV_CV(ICV),i),ICV=1,NCV)
      endif
      enddo
      endif
!
      if(mol_rslt==1) then
        allocate(vlr(MXALLCV,ncomp),stat=ierr)
        if(ierr.ne.0) then
          write(ifle,*) 'allocating array error in write_result=>yys'
          call FFRABORT(1,'allocating array at write_animation')
        endif
        do 200 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        do 210 ICVL=ICVS,ICVE
        dum1=0.d0
        do ICOM=1,ncomp
          dum1=dum1+yys(ICVL,ICOM)*r_wm(ICOM)
        enddo
        do ICOM=1,ncomp
          dum2=wm(ICOM)*dum1
          vlr(ICVL,ICOM)=yys(ICVL,ICOM)/dum2
        enddo
 210    enddo
 200    enddo
!
        do i=1,ncomp
        if(ianim_comp(i)==1)then
          write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
        endif
        enddo
        deallocate(vlr)
      else
        do i=1,ncomp
        if(ianim_comp(i)==1)then
          write(ifl) (yys(LCV_CV(ICV),i),ICV=1,NCV)
        endif
        enddo
      endif
!
!---------------------------------------------------------
!
      if(ical_suf==1) then
        if(nsdot>0) then
          allocate(vlr(MXCV,nsdot),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'allocating array error in write_result=>suf'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          iicom=0
          do nb=1,nbcnd
!
          ami(:)=0
          do IRC=1,nneq
          if(IDX_SUFRAC(nb,IRC)) then
            iph=1
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            do  ipcom=icoms,icome
            icom=phs_com(ipcom)
            if(vreq(ICOM,2,IRC)>SML.or.vreq(ICOM,1,IRC)>SML) then
              ami(icom)=1
            endif
            enddo
          endif
          enddo
!
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==gasphase) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ami(icom)==1.and.ianim_SDOT_suf(icom)==1)then
                iicom=iicom+1
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICVL=LVEDGE(1,ICFL)
                dum2=SDOT(IBFL,icom)
                vlr(ICVL,iicom)=dum2         !dum2=mol/m^2/s
                enddo
              endif
              enddo
            elseif(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_SDOT_suf(icom)==1) then
                iicom=iicom+1
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICVL=LVEDGE(1,ICFL)
                dum2=SDOT(IBFL,icom)         !dum2=mol/m^2/s
                vlr(ICVL,iicom)=dum2         !dum2=mol/m^2/s
                enddo
              endif
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_SDOT_suf(icom)==1) then
                iicom=iicom+1
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICVL=LVEDGE(1,ICFL)
                dum2=SDOT(IBFL,icom)         !dum2=mol/m^2/s
                vlr(ICVL,iicom)=dum2         !dum2=mol/m^2/s
                enddo
              endif
              enddo            
            endif
            enddo
          endif
          enddo
          if(nsdot/=iicom) then
            call FFRABORT(1,'ERR: nsdot/=iicom')
          endif
          do i=1,nsdot
          write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
          enddo
          deallocate(vlr)
        endif
!
        if(nmolfr>0) then
          allocate(vlr(MXCV,nmolfr),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'allocating array error in write_result=>suf'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_molefr_suf(icom)==1) then
                iicom=iicom+1
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICVL=LVEDGE(1,ICFL)
                vlr(ICVL,iicom)=MOLFRC(IBFL,icom)
                enddo
              endif
              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_molefr_suf(icom)==1) then
                iicom=iicom+1
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICVL=LVEDGE(1,ICFL)
                vlr(ICVL,iicom)=MOLFRC(IBFL,icom)
                enddo
              endif
              enddo            
            endif
            enddo
          endif
          enddo
          if(nmolfr/=iicom) then
            call FFRABORT(1,'ERR: nmolfr/=iicom')
          endif
          do i=1,nmolfr
          write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
          enddo
          deallocate(vlr)
        endif
!
        if(ndeps>0) then
          allocate(vlr(MXCV,ndeps),stat=ierr)
          if(ierr.ne.0) then
            write(ifle,*) 'allocating array error in write_result=>suf'
            call FFRABORT(1,'allocating array at write_result')
          endif
          vlr(:,:)=0.d0
          iicom=0
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
!          if(kd==kdcvd) then
          if(surfreac(nb)==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            do iph=1,nphase
            icoms=phs_idx(iph-1)+1
            icome=phs_idx(iph)
            nbx=phs_nbc(iph)
            ntp=phs_typ(iph)
            if(ntp==surphase) then
!              do ipcom=icoms,icome
!              icom=phs_com(ipcom)
!              iicom=iicom+1
!              do IBFL=IBFS,IBFE
!              ICFL=LBC_SSF(IBFL)
!              ICVL=LVEDGE(1,ICFL)
!              dum2=SDOT(IBFL,icom)        !dum2=mol/m^2/s
!              vlr(ICVL,iicom)=0.d0        !dum2=mol/m^2/s
!              enddo
!              enddo
            elseif(ntp==blkphase.and.nbx==nb) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              if(ianim_depoeth_spd(icom)==1) then
                iicom=iicom+1
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICVL=LVEDGE(1,ICFL)
                dum2=SUFBND(IBFL,icom)
                vlr(ICVL,iicom)=dum2
                enddo
              endif
              enddo
            endif
            enddo
          endif
          enddo
!
          if(ndeps/=iicom) then
            call FFRABORT(1,'ERR: ndeps/=iicom')
          endif
          do i=1,ndeps
            write(ifl) (vlr(LCV_CV(ICV),i),ICV=1,NCV)
          enddo
          deallocate(vlr)
        endif
!
        do icom=1,ncomp
        if(ianim_GAS_WDOT(icom)==1) then
          write(ifl) 
     &    (r_wm(icom)*WDOT(LCV_CV(ICV),icom),ICV=1,NCV)
        endif
        enddo
!
        DEALLOCATE(ami)
      endif
!

!
!---------------------------------------------------------
      close(ifl)
!
      return
!
 9991 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'lack of restart file name'
      write(ifle,*) '(write_anim)'
      ierror=1
!
      return
!
      end subroutine write_anim
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_cdcl(iter,time,delt,ffexit,
     &     LVEDGE,LBC_SSF,SFAREA,SFCENT,
     &     prs,rho,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only         : getfil,lenfnm,ifle,ifll
      use module_output,only     : outchk,fnmlst,none,yes,kopt,
     &                             start_cdcl
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX
      use module_cdcl
      use module_FFRGFwords

      implicit none
!
!--- [dummy arguments]
!
      integer,intent(in) :: iter
      real*8 ,intent(in) :: time,delt
      logical,intent(in) :: ffexit
!     
      integer,intent(in) :: LVEDGE(2,MXCVFAC)
      integer,intent(in) :: LBC_SSF(MXSSFBC)
      real(8),intent(in) :: SFAREA(4,MXCVFAC)
      real(8),intent(in) :: SFCENT(3,MXCVFAC)
      real(8),intent(in) :: prs(MXALLCV,2)
      real(8),intent(in) :: rho(MXALLCV,2)
      integer,intent(out):: ierror
!
! --- [local entities]
!
      real(8)      :: mflowVel(1:3)
      logical,save :: headerOut=.true.
      logical      :: fexist=.false.
      integer,save :: ifl=-1
      integer      :: returnStat
      character(lenfnm),save :: fnam
      integer      :: srOut,ierr,ios
      integer      :: ic,jc,itmp,nb
      real(8)      :: cdValue(1:numofCdClWall+1),
     &                clValue(1:numofCdClWall+1)
      real(8)      :: dragForce(1:numofCdClWall+1)
      real(8)      :: liftForce(1:numofCdClWall+1,1:3)
      real(8)      :: rtmp1,rtmp2
      real(8),save :: mvnorm(1:3),ndotu
      integer      :: IBFS,IBFE,IBFL,ICFL,IDC
      logical,save :: calcAreas=.true.
      character*80 :: formatKey,formatKeytmp
      character*6  :: text
      logical      :: fileExistFlag
      logical,save :: title=.false.
!
! --- temporally, inflow velocity  is assigned by fort.1 file.
!
      if(.not.cdcloutputFlag) return
!
      mflowVel(1:3)=inflowVelCdCl(1:3)
      
! --- calculate projected area of object FROM HERE ---------------------
      if(calcAreas) then
        allocate(projAreaS(1:numofCdClWall+1))
        projAreaS(:)=0.d0
        mvnorm(1:3)=mflowVel(1:3)/
     &  (SML+dsqrt(dot_product(mflowVel,mflowVel)))
        do ic=1,numofCdClWall
          nb=CdClWall_bcNo(ic)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            IDC=LVEDGE(2,ICFL)
            ndotu=dot_product(SFAREA(1:3,ICFL),mvnorm(1:3))
            projAreaS(ic)=projAreaS(ic)+
     &            0.5d0*dabs(SFAREA(4,ICFL)*ndotu)
          end do
        end do
        projAreaS(numofCdClWall+1)=sum(projAreaS(1:numofCdClWall))
!
        if(NPE.gt.1) then
          call hpcrasum(projAreaS,numofCdClWall+1)
        end if
!
        calcAreas=.false.
      end if
! --- calculate projected area of object TO HERE ---------------------
          
! --- Calculate Drag and Lift force FROM HERE --------------------------
      dragForce(:)=0.d0
      liftForce(:,:)=0.d0
      mvnorm(1:3)=mflowVel(1:3)/
     &  (SML+dsqrt(dot_product(mflowVel,mflowVel)))
      do ic=1,numofCdClWall
        nb=CdClWall_bcNo(ic)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ndotu=dot_product(SFAREA(1:3,ICFL),mvnorm(1:3))
          dragForce(ic)=dragForce(ic)-
     &            SFAREA(4,ICFL)*prs(IDC,1)*ndotu
     &            /rho(IDC,1)

          liftForce(ic,1:3)=liftForce(ic,1:3)+
     &            SFAREA(4,ICFL)*prs(IDC,1)*
     &                  (SFAREA(1:3,ICFL)-ndotu*mvnorm(1:3))
     &            /rho(IDC,1)
        end do
      end do
      dragForce(numofCdClWall+1)=sum(dragForce(1:numofCdClWall))
      do ic=1,3
        liftForce(numofCdClWall+1,ic)=sum(liftForce(1:numofCdClWall,ic))
      end do
      if(NPE.gt.1) then
        call hpcrasum(dragForce,numofCdClWall+1)
        call hpcrasum(liftForce(:,1),numofCdClWall+1)
        call hpcrasum(liftForce(:,2),numofCdClWall+1)
        call hpcrasum(liftForce(:,3),numofCdClWall+1)
      end if
! --- Calculate Drag and Lift force TO HERE --------------------------

! --- Calculate CD and CL force FROM HERE ------------------------------
      cdValue(:)=0.d0
      clValue(:)=0.d0
      cdValue(1:numofCdClWall+1)=dragForce(1:numofCdClWall+1)/
     &      (dot_product(mflowVel,mflowVel)*projAreaS*0.5d0)
      do ic=1,numofCdClWall+1
        clValue(ic)=
     &      dsqrt(dot_product(liftForce(ic,1:3),liftForce(ic,1:3)))/
     &      (dot_product(mflowVel,mflowVel)*projAreaS(ic)*0.5d0)
      end do
! --- Calculate CD and CL force TO HERE ------------------------------
!
! --- output ----------------------------------------------------------
! --- Preparation of output file FROM HERE ----------------------------
!     Decide output file name and logical number. Erase nonused files.
!
      if(my_rank.eq.root) then
        if(ifl.lt.0) then
          call getfil(ifl,fnam,'cdcl')
          if(fnam.eq.' ') then
            write(ifle,*) 'CdCl file name is not defined in ffloe.ctl'
            call FFRABORT(1,'write_cdcl')
          endif
          fnam=adjustl(fnam)
          inquire(file=fnam,exist=fexist)
          if(fexist) then
            if(kopt(6).eq.1) then
              if(iter.le.1+nint(start_cdcl)) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            elseif(kopt(6).eq.3) then
              if(time.le.delt+start_cdcl) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            endif
          endif
        endif
      end if

!     Check whether output data or not.
      call outchk(fnmlst(6),iter,time,delt,srOut)
      if(srOut.eq.none) return
      if(.not.ffexit .and. srOut.ne.yes) return

!     Add step number to file names (Not used here)
      if(my_rank.eq.root) then
        if(kopt(6).eq.2.or.kopt(6).eq.4) then
          write(text,'(i6)') iter
          text=adjustl(text)
          fnam=TRIM(TRIM(fnam)//'_'//TRIM(text))
        endif
      endif

!     Output header of fluid force data.
      if(my_rank.eq.root) then
        if(headerOut) then
          inquire(file=fnam,exist=fileExistFlag)
          if(.not.(fileExistFlag)) then
            call ffr_cdcl_headerout(returnStat)
            if(returnStat/=0) then
              ierror=1
              return
            end if
          end if
          headerOut=.false.
        end if
      end if
      
!      output format
      write(formatKeytmp,*)numofCdClWall+1
      write(formatKey,*)'(i8,1x,',trim(adjustl(formatKeytmp)),
     &                                                '(1x,e16.9))'
! --- Preparation of output file TO HERE -------------------------------

      if(my_rank.eq.root) then
        open(ifl,FILE=fnam,form='formatted',position='append')
        write(ifl,'(a)')customData_FFGF
        write(ifl,'(a)')ffrTseries_FFGF
        write(ifl,'(a)')'FFR Time Series data field'
        write(ifl,*)0
        write(ifl,'(a)')tseriesData_FFGF
        write(ifl,*)1
        write(ifl,'(e16.9,1x,i16)') time,iter
        write(ifl,formatKey)1,(cdValue(ic),ic=1,numofCdClWall+1)
        write(ifl,formatKey)2,(clValue(ic),ic=1,numofCdClWall+1)
!
        write(ifl,'(a)')tseriesDatae_FFGF
        close(ifl)
      end if
      return
      
      contains
!=======================================================================
      subroutine ffr_cdcl_headerout(statRtrn)
      integer,intent(out) :: statRtrn

!     Output format for the position of an observer point.
      
      open(ifl,file=fnam,form='formatted',
     &  status='unknown',iostat=ios)
      write(ifl,'(a)')asciiv2_FFGF
      write(ifl,'(a)')newSet_FFGF
      write(ifl,'(a)')'FrontFlow Red Time Series data'
      write(ifl,'(a)')customData_FFGF
      write(ifl,'(a)')ffrTseries_FFGF
      write(ifl,'(a)')'FFR Time Series Header field'
      write(ifl,*)0
      write(ifl,'(a)')tseriesHead_FFGF
      write(ifl,*)1
      write(ifl,*)2
      write(ifl,*)1,'f',numofCdClWall+1,'CD'
      write(ifl,*)2,'f',numofCdClWall+1,'CL'
      write(ifl,'(a)')tseriesHeade_FFGF
      write(ifl,'(a)')customData_FFGF
      write(ifl,'(a)')ffrTseries_FFGF
      write(ifl,'(a)')'FFR Time Series observer field'
      write(ifl,*)0
      write(ifl,'(a)')tseriesObserv_FFGF
      write(ifl,*)1
      write(ifl,*)2
      write(ifl,'(i4,a,1x,3(e16.9,1x))')1,' I',0.0,0.0,0.0
      write(ifl,'(i4,a,1x,3(e16.9,1x))')2,' I',0.0,0.0,0.0
      write(ifl,'(a)')tseriesObserve_FFGF
      close(ifl)
      statRtrn=0
      end subroutine ffr_cdcl_headerout
!=======================================================================
      
      end subroutine write_cdcl

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_probe(iter,time,delt,ffexit,
     &  MAT_CVEXT,MAT_INDEX,CVCENT,
     &  prs,pp0,rho,tmp,yys,vel,aks,
     &  ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_FFRGFwords
      use module_io,only         : getfil,lenfnm,ifle,ifll
      use module_output,only     : outchk,fnmlst,none,yes,kopt,
     &                             start_prob,ProbeFiletype
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX
      use module_probe
      use module_les   ,only     : mol_rslt
      use module_species,only  : spcnam,wm,r_wm

      implicit none
!
!--- [dummy arguments] 
!
      integer,intent(in)    :: iter
      real*8 ,intent(in)    :: time,delt
      logical,intent(in)    :: ffexit
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_INDEX( 0:MXMAT)
      real*8 ,intent(in)    :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: prs    (   MXALLCV)
      real*8 ,intent(in)    :: pp0    (   MXALLCV)
      real*8 ,intent(in)    :: rho    (   MXALLCV)
      real*8 ,intent(in)    :: tmp    (   MXALLCV)
      real*8 ,intent(in)    :: yys    (   MXALLCV,MXcomp)
      real(8),intent(in)    :: vel    (   MXALLCV,3)
      real*8 ,intent(in)    :: aks    (   MXALLCVR,MXrans)
      integer,intent(out)   :: ierror
!
! --- [local entities]
! 
      logical,save :: headerOut=.true.
      logical,save :: searchCVflag=.true.
      logical      :: fexist=.false.
      integer,save :: ifl=-1
      integer      :: returnStat,IIMAT
      integer      :: icvl,icve,icvs
      integer      :: srOut,ierr,ios
      integer      :: ic,jc,itmp,nb,probCV,ICOM
      real(8)      :: probDist,rtmp1,rtmp2,dum1
      real(8),save :: mvnorm(1:3),ndotu
      logical,save :: calcAreas=.true.
      character*80 :: formatKey,formatKeytmp
      character*6  :: text
      logical      :: fileExistFlag
      logical,save :: title=.false.
      character(lenfnm),save :: fnam,fnam1
      integer,allocatable :: probCV_MPtmp(:)
      real(8),allocatable :: probDist_MPtmp(:)
      integer      :: iii,jjj
      real(8)      :: ave,ave_time
      
!-------------------------------------------
! --- search CV index at probe position.
!-------------------------------------------
      if(.not.probeFlag) return
!
      if(searchCVflag) then
        if(NPE.GT.1) then
          allocate(probCV_MPtmp(1:NPE))
          allocate(probDist_MPtmp(1:NPE))
          probCV_MPtmp=0
          probDist_MPtmp=0.d0
        end if
!
        do IC=1,numofProbes
          do IIMAT=1,NMAT
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_INDEX(IIMAT)
            if(IIMAT==1) then
              probDist=
     &            dsqrt((CVCENT(1,icvs)-probePosit(ic,1))**2+
     &                  (CVCENT(2,icvs)-probePosit(ic,2))**2+
     &                  (CVCENT(3,icvs)-probePosit(ic,3))**2)
              probCV=ICVS
            endif
            do ICVL=ICVS,ICVE
              rtmp1=
     &            dsqrt((CVCENT(1,ICVL)-probePosit(ic,1))**2+
     &                  (CVCENT(2,ICVL)-probePosit(ic,2))**2+
     &                  (CVCENT(3,ICVL)-probePosit(ic,3))**2)
              if(probDist>rtmp1) then
                probDist=rtmp1;     probCV=icvl
              end if
            end do
          end do
          if(NPE==1) then
            probePositIn(ic,1:3)=cvcent(1:3,probCV)
            probeCVIndex(ic,1)=probCV
            probeCVIndex(ic,2)=root
          else if(NPE>1) then
            call hpc_igather(NPE,probCV,probCV_MPtmp)
            call hpc_dgather(NPE,probDist,probDist_MPtmp)
            probeCVIndex(ic,2:2)=minloc(probDist_MPtmp)-1
            probeCVIndex(ic,1)=probCV_MPtmp(probeCVIndex(ic,2)+1)
            if(probeCVIndex(ic,2)/=root) then
              if(my_rank==root) then
                call hpc_drecv(3,probePositIn(ic,1:3),
     &                                    probeCVIndex(ic,2),ic)
              else if(my_rank==probeCVIndex(ic,2)) then
                call hpc_dsend(3,cvcent(1:3,probeCVIndex(ic,1)),root,ic)
              end if
            else
              probePositIn(ic,1:3)=cvcent(1:3,probeCVIndex(ic,1))
            end if
          endif
          if(MY_RANK==ROOT) then
            write(*,*)'Output Position'
            write(*,*)ic,probeCVIndex(ic,2),probeCVIndex(ic,1)
            write(*,*)probePositIn(ic,1:3)
            write(*,*)probePosit(ic,1:3)
          end if
        end do
        if(NPE>1) then
          deallocate(probCV_MPtmp);   deallocate(probDist_MPtmp)
        end if
        searchCVflag=.false.
      end if
!      
      do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then
          if(probeDataType(ic,1)==1) then 
!---------------
! --- pressure
!---------------
            probeData(ic,1)=prs(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==2) then      
!---------------
! --- Velocity
!---------------
            probeData(ic,1:3)=vel(probeCVIndex(ic,1),1:3)
          else if(probeDataType(ic,1)==3) then      
!---------------
! --- temparature
!---------------
            probeData(ic,1)=tmp(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==4) then      
!---------------
! --- density
!---------------
            probeData(ic,1)=rho(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==5) then      
!---------------
! --- thermal pressure
!---------------
            probeData(ic,1)=pp0(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==6) then      
!--------------------------
! --- Mass/Mole Fraction
!--------------------------
            if(mol_rslt==0) then
              probeData(ic,1)=yys(probeCVIndex(ic,1),NO_YS(ic))
            else
              dum1=0.d0  !yys(probeCVIndex(ic,1),NO_YS(ic))
              do ICOM=1,ncomp
                dum1=dum1+yys(probeCVIndex(ic,1),ICOM)*r_wm(ICOM)
              enddo
              dum1=wm(NO_YS(ic))*dum1
              probeData(ic,1)=yys(probeCVIndex(ic,1),NO_YS(ic))/dum1
            endif
          else if(probeDataType(ic,1)==7) then 
!---------------
! --- averaged pressure
!---------------
            iii=MOD(iter,ave_step(ic))
            if(iii.eq.0)iii=ave_step(ic)
            probeData(ic,iii)=prs(probeCVIndex(ic,1))
          end if
       endif
!
        if(probeCVIndex(ic,2)/=root) then
          if(my_rank==root) then
            call hpc_drecv(probeDataType(ic,2),
     &            probeData(ic,1:probeDataType(ic,2)),
     &                        probeCVIndex(ic,2),ic)
          else if(MY_RANK==probeCVIndex(ic,2)) then
            call hpc_dsend(probeDataType(ic,2),
     &            probeData(ic,1:probeDataType(ic,2)),
     &                        root,ic)
          end if
        end if
      end do
!
! --- --------------------------------------------
! --- Preparation of output file FROM HERE
!     Decide output file name and logical number. 
!     Erase non-used files.
! --- --------------------------------------------
!
      if(MY_RANK.EQ.ROOT) then
        if(ifl.lt.0) then
          call getfil(ifl,fnam,'probe')
          if(fnam.eq.' ') then
            write(ifle,*) 'Probe file name is not defined in fort.1'
            call FFRABORT(1,'write_probe')
          endif
          fnam=adjustl(fnam)
          inquire(file=fnam,exist=fexist)
          if(fexist) then
            if(kopt(7).eq.1) then
              if(iter.le.1+nint(start_prob)) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            elseif(kopt(7).eq.3) then
              if(time.le.delt+start_prob) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            endif
          endif
        endif
      end if
!----------------------------------------
! --- Check whether output data or not.
!----------------------------------------
      call outchk(fnmlst(7),iter,time,delt,srOut)
      if(srOut.eq.none) return
      if(.not.ffexit .and. srOut.ne.yes) return
!---------------------------------------------------
! --- Add step number to file names (Not used here)
!---------------------------------------------------
      if(MY_RANK.EQ.ROOT) then
        if(kopt(7).eq.2.or.kopt(7).eq.4) then
          write(text,'(i6)') iter
          text=adjustl(text)
          fnam=TRIM(TRIM(fnam)//'_'//TRIM(text))
        endif
      endif
! ----------------------------------
! --- Output header of probe data.
! ----------------------------------
      if(MY_RANK.EQ.ROOT) then
        if(headerOut) then
          inquire(file=fnam,exist=fileExistFlag)
          if(.not.(fileExistFlag)) then
            if(ProbeFiletype==2) then
! --- Excel
              do ic=1,numofProbes
              fnam1=trim(fnam)//'_'//trim(probeDataLabel(ic))
              if(probeDataLabel(ic)(1:1)=='M') then
               fnam1=trim(fnam1)//'_'//TRIM(adjustl(spcnam(NO_YS(ic))))
              endif
              open(ifl,file=fnam1,form='formatted',
     &            status='unknown',iostat=ios)
              write(ifl,'(a)')'#probe positions: '
              write(ifl,'(i16,a9,3(e16.9,1x))') ic,' (x,y,z):',
     &               probePositIn(ic,1:3)
              close(ifl)
              enddo
            else
! --- GF
              call ffr_probe_headerout(returnStat)
              if(returnStat/=0) then
                ierror=1
                return
              end if
            end if
          end if
          headerOut=.false.
        endif
      endif
! ---------------------------------------
! --- Preparation of output file TO HERE 
! ---------------------------------------
      if( probeDataType(1,1).ne.7 ) then 
        if(MY_RANK.EQ.ROOT) then
! --- Excel
          if(ProbeFiletype==2) then
            do ic=1,numofProbes
              fnam1=trim(fnam)//'_'//trim(probeDataLabel(ic))
              open(ifl,FILE=fnam1,form='formatted',position='append')
              write(formatKey,'(a,I4,a)')'(e16.9,1x,i16,1x,',
     &          probeDataType(ic,2),'(e16.9,1x))'

              write(ifl,formatKey) time,iter,
     &          (probeData(ic,jc),jc=1,probeDataType(ic,2))
              close(ifl)
            enddo
          else
            open(ifl,FILE=fnam,form='formatted',position='append')
            write(ifl,'(a)')trim(customData_FFGF)
            write(ifl,'(a)')trim(ffrTseries_FFGF)
            write(ifl,'(a)')'FFR Time Series data field'
            write(ifl,*)0
            write(ifl,'(a)')trim(tseriesData_FFGF)
            write(ifl,*) 1
            write(ifl,'(e16.9,1x,i16)')time,iter
            do ic=1,numofProbes
              write(formatKeytmp,*) probeDataType(ic,2)
              write(formatKey,*)'(i8,',trim(adjustl(formatKeytmp)),
     &           '(1x,e16.9))'
              write(ifl,formatKey)ic,
     &           (probeData(ic,jc),jc=1,probeDataType(ic,2))
            end do
            write(ifl,'(a)')trim(tseriesDatae_FFGF)
            close(ifl)
          endif
        endif

      else

        iii=MOD(iter,ave_step(1))
        if(iii.eq.0)then

          if(MY_RANK.EQ.ROOT) then
! --- Excel
            if(ProbeFiletype==2) then
              do ic=1,numofProbes
                fnam1=trim(fnam)//'_'//trim(probeDataLabel(ic))
                open(ifl,FILE=fnam1,form='formatted',position='append')
                write(formatKey,'(a,I4,a)')'(e16.9,1x,i16,1x,',
     &            probeDataType(ic,2),'(e16.9,1x))'
                ave=0.0
                do jjj=1,ave_step(ic)
                  ave=ave+probeData(ic,jjj)
                enddo
                ave=ave/ave_step(ic)
                ave_time=time-(delt*ave_step(ic))/2.0
                write(ifl,formatKey) ave_time,iter,ave
                close(ifl)
              enddo
            else
              open(ifl,FILE=fnam,form='formatted',position='append')
              write(ifl,'(a)')trim(customData_FFGF)
              write(ifl,'(a)')trim(ffrTseries_FFGF)
              write(ifl,'(a)')'FFR Time Series data field'
              write(ifl,*)0
              write(ifl,'(a)')trim(tseriesData_FFGF)
              write(ifl,*) 1
              ave_time=time-(delt*ave_step(1))/2.0
              write(ifl,'(e16.9,1x,i16)')ave_time,iter
              do ic=1,numofProbes
                write(formatKeytmp,*) probeDataType(ic,2)
                write(formatKey,*)'(i8,',trim(adjustl(formatKeytmp)),
     &            '(1x,e16.9))'
                ave=0.0
                do jjj=1,ave_step(ic)
                  ave=ave+probeData(ic,jjj)
                enddo
                ave=ave/ave_step(ic)
                write(ifl,formatKey)ic,ave
              end do
              write(ifl,'(a)')trim(tseriesDatae_FFGF)
              close(ifl)
            endif
          endif

        endif

      endif

      return
      
      contains
!=======================================================================
      subroutine ffr_probe_headerout(statRtrn)
!=======================================================================
      integer,intent(out) :: statRtrn

!     Output format for the position of an observer point.
      
      open(ifl,file=fnam,form='formatted',
     &  status='unknown',iostat=ios)
      write(ifl,'(a)')asciiv2_FFGF
      write(ifl,'(a)')newSet_FFGF
      write(ifl,'(a)')'FrontFlow Red Time Series data'
      write(ifl,'(a)')customData_FFGF
      write(ifl,'(a)')ffrTseries_FFGF
      write(ifl,'(a)')'FFR Time Series Header field'
      write(ifl,*)0
      write(ifl,'(a)')tseriesHead_FFGF
      write(ifl,*)1
      write(ifl,*)numofProbes
      do ic=1,numofProbes
        if(probeDataType(ic,1).ne.7) then
          write(ifl,'(i5,a,i4,1x,a)')ic,' f',probeDataType(ic,2),
     &                                    probeDataLabel(ic)
        else
          write(ifl,'(i5,a,i4,1x,a)')ic,' f',1, probeDataLabel(ic)
        endif
      end do
      write(ifl,'(a)')tseriesHeade_FFGF
      write(ifl,'(a)')customData_FFGF
      write(ifl,'(a)')ffrTseries_FFGF
      write(ifl,'(a)')'FFR Time Series observer field'
      write(ifl,*)0
      write(ifl,'(a)')tseriesObserv_FFGF
      write(ifl,*)1
      write(ifl,*)numofProbes
      do ic=1,numofProbes
        write(ifl,'(i5,a,1x,3(e16.9,1x))')ic,' P',probePositIn(ic,1:3)
      end do
      write(ifl,'(a)')tseriesObserve_FFGF
      close(ifl)
      statRtrn=0
      end subroutine ffr_probe_headerout
!=======================================================================
      
      end subroutine write_probe

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_probe_para(iter,time,delt,ffexit,
     &  MAT_CVEXT,MAT_INDEX,CVCENT,
     &  prs,pp0,rho,tmp,yys,vel,aks,
     &  ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_FFRGFwords
      use module_io,only         : getfil,lenfnm,ifle,ifll
      use module_output,only     : outchk,fnmlst,none,yes,kopt,
     &                             start_prob,ProbeFiletype
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX
      use module_probe
      use module_les   ,only     : mol_rslt
      use module_species,only  : spcnam,wm,r_wm

      implicit none
!
!--- [dummy arguments] 
!
      integer,intent(in)    :: iter
      real*8 ,intent(in)    :: time,delt
      logical,intent(in)    :: ffexit
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_INDEX( 0:MXMAT)
      real*8 ,intent(in)    :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: prs    (   MXALLCV)
      real*8 ,intent(in)    :: pp0    (   MXALLCV)
      real*8 ,intent(in)    :: rho    (   MXALLCV)
      real*8 ,intent(in)    :: tmp    (   MXALLCV)
      real*8 ,intent(in)    :: yys    (   MXALLCV,MXcomp)
      real(8),intent(in)    :: vel    (   MXALLCV,3)
      real*8 ,intent(in)    :: aks    (   MXALLCVR,MXrans)
      integer,intent(out)   :: ierror
!
! --- [local entities]
! 
      logical,save :: headerOut=.true.
      logical,save :: searchCVflag=.true.
      logical      :: fexist=.false.
      integer,save :: ifl=-1
      integer      :: returnStat,IIMAT
      integer      :: icvl,icve,icvs
      integer      :: srOut,ierr,ios
      integer      :: ic,jc,itmp,nb,probCV,ICOM
      real(8)      :: probDist,rtmp1,rtmp2,dum1
      real(8),save :: mvnorm(1:3),ndotu
      logical,save :: calcAreas=.true.
      character*80 :: formatKey,formatKeytmp
      character*6  :: text
      logical      :: fileExistFlag
      logical,save :: title=.false.
      character(lenfnm),save :: fnam,fnam1
      integer,allocatable :: probCV_MPtmp(:)
      real(8),allocatable :: probDist_MPtmp(:)
      integer      :: iii,jjj
      real(8)      :: ave,ave_time
      
!-------------------------------------------
! --- search CV index at probe position.
!-------------------------------------------
      if(.not.probeFlag) return
!
      if(searchCVflag) then
        if(NPE.GT.1) then
          allocate(probCV_MPtmp(1:NPE))
          allocate(probDist_MPtmp(1:NPE))
          probCV_MPtmp=0
          probDist_MPtmp=0.d0
        end if
!
        do IC=1,numofProbes
          do IIMAT=1,NMAT
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_INDEX(IIMAT)
            if(IIMAT==1) then
              probDist=
     &            dsqrt((CVCENT(1,icvs)-probePosit(ic,1))**2+
     &                  (CVCENT(2,icvs)-probePosit(ic,2))**2+
     &                  (CVCENT(3,icvs)-probePosit(ic,3))**2)
              probCV=ICVS
            endif
            do ICVL=ICVS,ICVE
              rtmp1=
     &            dsqrt((CVCENT(1,ICVL)-probePosit(ic,1))**2+
     &                  (CVCENT(2,ICVL)-probePosit(ic,2))**2+
     &                  (CVCENT(3,ICVL)-probePosit(ic,3))**2)
              if(probDist>rtmp1) then
                probDist=rtmp1;     probCV=icvl
              end if
            end do
          end do
          if(NPE==1) then
            probePositIn(ic,1:3)=cvcent(1:3,probCV)
            probeCVIndex(ic,1)=probCV
            probeCVIndex(ic,2)=root
          else if(NPE>1) then
            call hpc_igather(NPE,probCV,probCV_MPtmp)
            call hpc_dgather(NPE,probDist,probDist_MPtmp)
            probeCVIndex(ic,2:2)=minloc(probDist_MPtmp)-1
            probeCVIndex(ic,1)=probCV_MPtmp(probeCVIndex(ic,2)+1)
            if(probeCVIndex(ic,2)/=root) then
              if(my_rank==root) then
                call hpc_drecv(3,probePositIn(ic,1:3),
     &                                    probeCVIndex(ic,2),ic)
              else if(my_rank==probeCVIndex(ic,2)) then
                call hpc_dsend(3,cvcent(1:3,probeCVIndex(ic,1)),root,ic)
              end if
            else
              probePositIn(ic,1:3)=cvcent(1:3,probeCVIndex(ic,1))
            end if
          endif
          call hpcrbcast_stream(root,probePositIn(ic,1),1)   !!!modify
          call hpcrbcast_stream(root,probePositIn(ic,2),1)   !!!modify
          call hpcrbcast_stream(root,probePositIn(ic,3),1)   !!!modify
          if(MY_RANK==ROOT) then
            write(*,*)'Output Position'
            write(*,*)ic,probeCVIndex(ic,2),probeCVIndex(ic,1)
            write(*,*)probePositIn(ic,1:3)
            write(*,*)probePosit(ic,1:3)
          end if
        end do
        if(NPE>1) then
          deallocate(probCV_MPtmp);   deallocate(probDist_MPtmp)
        end if
        searchCVflag=.false.
      end if
!      
      do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then
          if(probeDataType(ic,1)==1) then 
!---------------
! --- pressure
!---------------
            probeData(ic,1)=prs(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==2) then      
!---------------
! --- Velocity
!---------------
            probeData(ic,1:3)=vel(probeCVIndex(ic,1),1:3)
          else if(probeDataType(ic,1)==3) then      
!---------------
! --- temparature
!---------------
            probeData(ic,1)=tmp(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==4) then      
!---------------
! --- density
!---------------
            probeData(ic,1)=rho(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==5) then      
!---------------
! --- thermal pressure
!---------------
            probeData(ic,1)=pp0(probeCVIndex(ic,1))
          else if(probeDataType(ic,1)==6) then      
!--------------------------
! --- Mass/Mole Fraction
!--------------------------
            if(mol_rslt==0) then
              probeData(ic,1)=yys(probeCVIndex(ic,1),NO_YS(ic))
            else
              dum1=0.d0  !yys(probeCVIndex(ic,1),NO_YS(ic))
              do ICOM=1,ncomp
                dum1=dum1+yys(probeCVIndex(ic,1),ICOM)*r_wm(ICOM)
              enddo
              dum1=wm(NO_YS(ic))*dum1
              probeData(ic,1)=yys(probeCVIndex(ic,1),NO_YS(ic))/dum1
            endif
          else if(probeDataType(ic,1)==7) then 
!---------------
! --- averaged pressure
!---------------
            iii=MOD(iter,ave_step(ic))
            if(iii.eq.0)iii=ave_step(ic)
            probeData(ic,iii)=prs(probeCVIndex(ic,1))
          end if
       endif
!
!!!        if(probeCVIndex(ic,2)/=root) then
!!!          if(my_rank==root) then
!!!            call hpc_drecv(probeDataType(ic,2),
!!!     &            probeData(ic,1:probeDataType(ic,2)),
!!!     &                        probeCVIndex(ic,2),ic)
!!!          else if(MY_RANK==probeCVIndex(ic,2)) then
!!!            call hpc_dsend(probeDataType(ic,2),
!!!     &            probeData(ic,1:probeDataType(ic,2)),
!!!     &                        root,ic)
!!!          end if
!!!        end if
      end do
!
! --- --------------------------------------------
! --- Preparation of output file FROM HERE
!     Decide output file name and logical number. 
!     Erase non-used files.
! --- --------------------------------------------
!
!!!      if(MY_RANK.EQ.ROOT) then
        if(ifl.lt.0) then
          call getfil(ifl,fnam,'probe')
          fnam=PROBEout
          if(fnam.eq.' ') then
            write(ifle,*) 'Probe file name is not defined in fort.1'
            call FFRABORT(1,'write_probe')
          endif
          fnam=adjustl(fnam)
          inquire(file=fnam,exist=fexist)
          if(fexist) then
            if(kopt(7).eq.1) then
              if(iter.le.1+nint(start_prob)) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            elseif(kopt(7).eq.3) then
              if(time.le.delt+start_prob) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            endif
          endif
        endif
!!!      end if
!----------------------------------------
! --- Check whether output data or not.
!----------------------------------------
      call outchk(fnmlst(7),iter,time,delt,srOut)
      if(srOut.eq.none) return
      if(.not.ffexit .and. srOut.ne.yes) return
!---------------------------------------------------
! --- Add step number to file names (Not used here)
!---------------------------------------------------
!      if(MY_RANK.EQ.ROOT) then
!        if(kopt(7).eq.2.or.kopt(7).eq.4) then
!          write(text,'(i6)') iter
!          text=adjustl(text)
!          fnam=TRIM(TRIM(fnam)//'_'//TRIM(text))
!        endif
!      endif
! ----------------------------------
! --- Output header of probe data.
! ----------------------------------
!      if(MY_RANK.EQ.ROOT) then
        if(headerOut) then
          inquire(file=fnam,exist=fileExistFlag)
          if(.not.(fileExistFlag)) then
            if(ProbeFiletype==2) then
! --- Excel
              do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then
              fnam1=trim(fnam)//'_'//trim(probeDataLabel(ic))
              if(probeDataLabel(ic)(1:1)=='M') then
               fnam1=trim(fnam1)//'_'//TRIM(adjustl(spcnam(NO_YS(ic))))
              endif
              open(ifl,file=fnam1,form='formatted',
     &            status='unknown',iostat=ios)
              write(ifl,'(a)')'#probe positions: '
              write(ifl,'(i16,a9,3(e16.9,1x))') ic,' (x,y,z):',
     &               probePositIn(ic,1:3)
              close(ifl)
         endif
              enddo
            else
! --- GF
              call ffr_probe_para_headerout(returnStat)
              if(returnStat/=0) then
                ierror=1
                return
              end if
            end if
          end if
          headerOut=.false.
        endif
!      endif
! ---------------------------------------
! --- Preparation of output file TO HERE 
! ---------------------------------------
      if( probeDataType(1,1).ne.7 ) then 
 !       if(MY_RANK.EQ.ROOT) then
! --- Excel
          if(ProbeFiletype==2) then
            do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then
              fnam1=trim(fnam)//'_'//trim(probeDataLabel(ic))
              open(ifl,FILE=fnam1,form='formatted',position='append')
              write(formatKey,'(a,I4,a)')'(e16.9,1x,i16,1x,',
     &          probeDataType(ic,2),'(e16.9,1x))'

              write(ifl,formatKey) time,iter,
     &          (probeData(ic,jc),jc=1,probeDataType(ic,2))
              close(ifl)
        endif
            enddo
          else
            open(ifl,FILE=fnam,form='formatted',position='append')
            write(ifl,'(a)')trim(customData_FFGF)
            write(ifl,'(a)')trim(ffrTseries_FFGF)
            write(ifl,'(a)')'FFR Time Series data field'
            write(ifl,*)0
            write(ifl,'(a)')trim(tseriesData_FFGF)
            write(ifl,*) 1
            write(ifl,'(e16.9,1x,i16)')time,iter
            do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then
              write(formatKeytmp,*) probeDataType(ic,2)
              write(formatKey,*)'(i8,',trim(adjustl(formatKeytmp)),
     &           '(1x,e16.9))'
              write(ifl,formatKey)ic,
     &           (probeData(ic,jc),jc=1,probeDataType(ic,2))
        endif
            end do
            write(ifl,'(a)')trim(tseriesDatae_FFGF)
            close(ifl)
          endif
 !       endif

      else

        iii=MOD(iter,ave_step(1))
        if(iii.eq.0)then

!          if(MY_RANK.EQ.ROOT) then
! --- Excel
            if(ProbeFiletype==2) then
              do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then
                fnam1=trim(fnam)//'_'//trim(probeDataLabel(ic))
                open(ifl,FILE=fnam1,form='formatted',position='append')
                write(formatKey,'(a,I4,a)')'(e16.9,1x,i16,1x,',
     &            probeDataType(ic,2),'(e16.9,1x))'
                ave=0.0
                do jjj=1,ave_step(ic)
                  ave=ave+probeData(ic,jjj)
                enddo
                ave=ave/ave_step(ic)
                ave_time=time-(delt*ave_step(ic))/2.0
                write(ifl,formatKey) ave_time,iter,ave
                close(ifl)
        endif
              enddo
            else
              open(ifl,FILE=fnam,form='formatted',position='append')
              write(ifl,'(a)')trim(customData_FFGF)
              write(ifl,'(a)')trim(ffrTseries_FFGF)
              write(ifl,'(a)')'FFR Time Series data field'
              write(ifl,*)0
              write(ifl,'(a)')trim(tseriesData_FFGF)
              write(ifl,*) 1
              ave_time=time-(delt*ave_step(1))/2.0
              write(ifl,'(e16.9,1x,i16)')ave_time,iter
              do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then
                write(formatKeytmp,*) probeDataType(ic,2)
                write(formatKey,*)'(i8,',trim(adjustl(formatKeytmp)),
     &            '(1x,e16.9))'
                ave=0.0
                do jjj=1,ave_step(ic)
                  ave=ave+probeData(ic,jjj)
                enddo
                ave=ave/ave_step(ic)
                write(ifl,formatKey)ic,ave
        endif
              end do
              write(ifl,'(a)')trim(tseriesDatae_FFGF)
              close(ifl)
            endif
!          endif

        endif

      endif

      return
      
      contains
!=======================================================================
      subroutine ffr_probe_para_headerout(statRtrn)
!=======================================================================
      integer,intent(out) :: statRtrn
      integer :: num_local       !!!modify

!     Output format for the position of an observer point.
      

      open(ifl,file=fnam,form='formatted',
     &  status='unknown',iostat=ios)
      write(ifl,'(a)')asciiv2_FFGF
      write(ifl,'(a)')newSet_FFGF
      write(ifl,'(a)')'FrontFlow Red Time Series data'
      write(ifl,'(a)')customData_FFGF
      write(ifl,'(a)')ffrTseries_FFGF
      write(ifl,'(a)')'FFR Time Series Header field'
      write(ifl,*)0
      write(ifl,'(a)')tseriesHead_FFGF
      write(ifl,*)1

      num_local=0       !!!modify
      do ic=1,numofProbes                       !!!modify
        if(my_rank==probeCVIndex(ic,2)) then
          num_local=num_local+1
        endif
      end do
!!!      write(ifl,*)numofProbes
      write(ifl,*)num_local

      do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then   !!!modify
          if(probeDataType(ic,1).ne.7) then
            write(ifl,'(i5,a,i4,1x,a)')ic,' f',probeDataType(ic,2),
     &                                         probeDataLabel(ic)
          else
            write(ifl,'(i5,a,i4,1x,a)')ic,' f',1, probeDataLabel(ic)
          endif
        endif
      end do
      write(ifl,'(a)')tseriesHeade_FFGF
      write(ifl,'(a)')customData_FFGF
      write(ifl,'(a)')ffrTseries_FFGF
      write(ifl,'(a)')'FFR Time Series observer field'
      write(ifl,*)0
      write(ifl,'(a)')tseriesObserv_FFGF
      write(ifl,*)1
!!!      write(ifl,*)numofProbes
      write(ifl,*)num_local   !!!modify

      do ic=1,numofProbes
        if(my_rank==probeCVIndex(ic,2)) then   !!!modify
          write(ifl,'(i5,a,1x,3(e16.9,1x))')ic,' P',probePositIn(ic,1:3)
        endif
      end do

      write(ifl,'(a)')tseriesObserve_FFGF
      close(ifl)
      statRtrn=0
      end subroutine ffr_probe_para_headerout
!=======================================================================
      
      end subroutine write_probe_para


!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine write_fluid_force
     &    (iter,time,delt,ffexit,
     &     MAT_NO,LVEDGE,LBC_SSF,SFAREA,SFCENT,
     &     FRSTCV,vel,prs,rmut,rho,utau,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_FFRGFwords
      use module_io,only         : getfil,lenfnm,ifle,ifll
      use module_output,only     : outchk,fnmlst,none,yes,kopt,
     &                             start_fforce,fforceFiletype
      use module_fluidforce,only  : ForceFlag,Nforce,UNIT_M,
     &                             NM_F,NM_F_WALL,NO_F_WALL,F_NUM,
     &                             end_mnt,begin_mnt,
     &                             F_fric,F_pres,F_SUM,M_SUM
      use module_boundary,only   : kdbcnd,LBC_INDEX,nbcnd,MAT_BCIDX
      use module_les  ,only      : ISTART
      use module_io,only         : dflxTmpFFlag
      use module_material,only   : ical_sld,rotati,ishaft,end,begin,
     &                             rot_ang
!
      implicit none
!
!--- [dummy arguments]
!
      integer,intent(in) :: iter
      real*8 ,intent(in) :: time,delt
      logical,intent(in) :: ffexit
!     
      INTEGER,INTENT(IN) :: MAT_NO(  0:MXMAT)
      integer,intent(in) :: LVEDGE(2,MXCVFAC)
      integer,intent(in) :: LBC_SSF(MXSSFBC)
      real(8),intent(in) :: SFAREA(4,MXCVFAC)
      real(8),intent(in) :: SFCENT(3,MXCVFAC)
      real(8),intent(in) :: FRSTCV(MXSSFBC)
      real(8),intent(in) :: vel(MXALLCV,3,2)
      real(8),intent(in) :: prs(MXALLCV,2)
      real(8),intent(in) :: rmut(MXALLCV)
      real(8),intent(in) :: utau(0:MXSSFBC)
      real*8 ,intent(in) :: rho(MXALLCV ,2)
      integer,intent(out) :: ierror
!
! --- [local entities]
!
!----------------------------------------------------
      integer :: isw,nb,IBFS,IBFE,IBFL,ICFL,ICV,IDC,srOut,ios,ic,jc
      integer,save :: ifl=-1
      character(lenfnm),save :: fnam
      character(lenfnm) :: fnam1
      logical,save :: headerOut=.true.,allo=.true.
      character*6 :: text
      logical :: fileExistFlag
      logical :: fexist=.false.
      integer :: returnStat,I,II,IS,IE,K
      logical,save,allocatable :: title1(:)
      logical,save :: title=.true.
      character(80):: formatKey,formatKeytmp
      real(8) :: dum1(1:3),tauw(3),uwall(1:3),uwallnrm,dum2
!
      tauw=0.d0
!
! 
 2502 format(3(i8,1x),3(e16.9,1x))
 2503 format(e16.9,1x,i16,1x,e16.9)
!
      ierror=0
      if(.NOT.ForceFlag) return
      if(allo) then
        allo=.false.
        allocate(title1(Nforce))
        title1(:)=.false.
        if(iter<=1+nint(start_fforce)) title1(:)=.true.
      endif
!----------------------
! --- force file
!----------------------
      if(my_rank.eq.root) then
        if(ifl.lt.0) then
          call getfil(ifl,fnam,'fforce')
          if(fnam.eq.' ') then
            write(ifle,*) 'Force file name is not defined in fort.1'
            call FFRABORT(1,'write_fluid_force')
          endif
          fnam=adjustl(fnam)
          inquire(file=fnam,exist=fexist)
          if(fexist) then
            if(kopt(8).eq.1) then
              if(iter.le.1+nint(start_fforce)) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            elseif(kopt(8).eq.3) then
              if(time.le.delt+start_fforce) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            else
              if(iter.le.1+nint(start_fforce)) then
                call unlink(fnam(:len_trim(fnam)))
                title=.true.
              endif
            endif
          endif
        endif
      endif
!
      if(my_rank==root.and.fforceFiletype==2) then
        do ic=1,Nforce
          fnam1=trim(fnam)//'_'//trim(NM_F(ic))
          fexist=.false.
          inquire(file=fnam1,exist=fexist)
          if(iter.le.1+nint(start_fforce)) then
            call unlink(fnam1(:len_trim(fnam1)))
            title1(ic)=.true.
          endif
        enddo
      endif
      
!------------------------------------------------------
!     Check whether output data or not.
!------------------------------------------------------
      call outchk(fnmlst(8),iter,time,delt,srOut)
      if(srOut.eq.none) return
      if(.not.ffexit .and. srOut.ne.yes) return
!------------------------------------------------------
! --- Add step number to file names (Not used here)
!------------------------------------------------------
      if(my_rank.eq.root) then
        if(kopt(8).eq.2.or.kopt(8).eq.4) then
          write(text,'(i6)') iter
          text=adjustl(text)
          fnam=TRIM(TRIM(fnam)//'_'//TRIM(text))
        endif
      endif
!------------------------------------------------------
! --- Output header of fluid force data.
!------------------------------------------------------
      if(my_rank.eq.root) then
        if(headerOut) then
          inquire(file=fnam,exist=fileExistFlag)
          if(.not.(fileExistFlag)) then
            call ffr_ff_headerout(returnStat)
            if(returnStat/=0) then
              ierror=1
              return
            end if
          end if
          headerOut=.false.
        endif
      endif
!----------------------------------------------------------------------
!---------------
! --- force
!---------------
      if(ForceFlag) then
!        F_fric(:,:)=0.d0
        F_pres(:,:)=0.d0
        F_SUM(:,:)=0.d0
        M_SUM=0.d0
        do I=1,Nforce
        IS=F_NUM(I-1)+1
        IE=F_NUM(I)
        do II=IS,IE
        nb=NO_F_WALL(II)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
! --- pressure
! --- fraction
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        dum1(1:3)=prs(IDC,1)*SFAREA(1:3,ICFL)*SFAREA(4,ICFL)        !pressure
        F_pres(1:3,I)=F_pres(1:3,I)+dum1(1:3)
!
!        uwall(1:3)=vel(ICV,1:3,1)-vel(IDC,1:3,1)
!        uwallnrm=uwall(1)*SFAREA(1,ICFL)
!     &          +uwall(2)*SFAREA(2,ICFL)
!     &          +uwall(3)*SFAREA(3,ICFL)
!        uwall(1:3)=uwall(1:3)-uwallnrm*SFAREA(1:3,ICFL)     
!        uwallnrm=dsqrt(uwall(1)*uwall(1)
!     &                +uwall(2)*uwall(2)
!     &                +uwall(3)*uwall(3))
!        uwall(1:3)=uwall(1:3)/(uwallnrm+SML)
!        tauw(1:3)=utau(IBFL)**2*rho(ICV,1)*uwall(1:3)*SFAREA(4,ICFL)  ! friction
!        F_fric(:,I)=F_fric(:,I)+tauw(:)
! --- Moment (r x F)
!        dum1(1:3)=tauw(1:3)+dum1(1:3)
        uwall(1:3)=SFCENT(1:3,ICFL)-begin_mnt(1:3,I)
        call cal_AXBEQC(uwall,dum1,tauw,dum2)


        do k=1,3
        M_SUM(k,I)=M_SUM(k,I)+tauw(k)
        enddo
        M_SUM(4,I)=M_SUM(4,I)
     &            +tauw(1)*UNIT_M(1,I)
     &            +tauw(2)*UNIT_M(2,I)
     &            +tauw(3)*UNIT_M(3,I)

        enddo
! --- sum
        enddo
        F_SUM(:,I)=F_pres(:,I)+F_fric(:,I)
        enddo
      endif
!---------------------------------------------------------------
!---------------------------------------------------------------
! --- Preparation of output file TO HERE -----------------------
!---------------------------------------------------------------
! --- Calculation and output of the fluid force FROM HERE ------
!     Calculate fluid force.
!---------------------------------------------------------------
      if(NPE.gt.1) then
        do I=1,Nforce
        call hpcrasum_0(F_SUM(:,I),3,0)
        call hpcrasum_0(F_fric(:,I),3,0)
        call hpcrasum_0(F_pres(:,I),3,0)
        call hpcrasum_0(M_SUM(:,I),4,0)
        enddo
      endif
!      dum2=0.062d0*0.5d0*1.2d0*50.d0**2
!      write(*,'(1X,a,4E14.5)') 'CD_P,CD_F,CD=',F_pres(1,1)/dum2,
!     & F_fric(1,1)/dum2,F_SUM(1,1)/dum2
!      write(*,'(1X,a,4E14.5)') 'CL_P,CL_F,CL=',F_pres(2,1)/dum2,
!     & F_fric(2,1)/dum2,F_SUM(2,1)/dum2
!      write(*,'(1X,a,4E14.5)') 'CS_P,CS_F,CS=',F_pres(3,1)/dum2,
!     & F_fric(3,1)/dum2,F_SUM(3,1)/dum2
!---------------------------------------------------------------
!     Output fluid force.
!---------------------------------------------------------------
      if(my_rank.eq.root) then
        if(fforceFiletype==1) then
          open(ifl,file=fnam,form='formatted',
     &         status='unknown',iostat=ios,position='append')
!-------------------------------
!     Output time, iterations.
!-------------------------------
          write(ifl,'(a)')customData_FFGF
          write(ifl,'(a)')ffrTseries_FFGF
          write(ifl,'(a)')'FFR Time Series data field'
          write(ifl,*)0
          write(ifl,'(a)')tseriesData_FFGF
          write(ifl,*)1
          write(ifl,'(e16.9,1x,i16)')time,iter
!-------------------------------
!     Output data 
!-------------------------------
          do ic=1,Nforce
            write(ifl,'(i8,1x,13(e16.9,1x))')ic,F_pres(1:3,ic),
     &            F_fric(1:3,ic),F_SUM(1:3,ic),M_SUM(1:4,ic)
          end do
          write(ifl,'(a)')tseriesDatae_FFGF
          close(ifl)
        else if(fforceFiletype==2) then
! --- 'Excel' format
           do ic=1,Nforce
           fnam1=trim(fnam)//'_'//trim(NM_F(ic))
           open(ifl,file=fnam1,form='formatted',
     &             status='unknown',position='append')

           if(title1(ic)) then
             title1(ic)=.false.
             write(ifl,'(1X,6a)') 'time, iter, no, ',
     &            'F_pres_x, F_pres_y, F_pres_z, ',

     &            'F_fric_x, F_fric_y, F_fric_z, ',
     &            'F_SUM_x,  F_SUM_y,   F_SUM_z, ',
     &            'F_SUM_x,  F_SUM_y,   F_SUM_z, ',
     &            'M_SUM_x,  M_SUM_y,   M_SUM_z, '
           endif
           write(ifl,'(e16.9,1x,i16,i8,1x,13(e16.9,1x))') 
     &            time,iter,ic,F_pres(1:3,ic),
     &            F_fric(1:3,ic),F_SUM(1:3,ic),M_SUM(1:4,ic)
           close(ifl)
           end do
        end if
      end if
!--------------------------------------------------------------
! --- Calculation and output of the fluid force TO HERE -------
!--------------------------------------------------------------
      ierror=0
      return
      
      contains
!=======================================================================
      subroutine ffr_ff_headerout(statRtrn)
!=======================================================================
      integer,intent(out) :: statRtrn
 2501 format(i8,1x,3(e16.9,1x))
!     Output format for the position of an observer point.

      if(fforceFiletype==1) then
! --- FFR-GF ASCII Ver 2.0 --------------
        open(ifl,file=fnam,form='formatted',
     &             status='unknown',iostat=ios)
        write(ifl,'(a)')asciiv2_FFGF
        write(ifl,'(a)')newSet_FFGF
        write(ifl,'(a)')'FrontFlow Red Time Series data'
        write(ifl,'(a)')customData_FFGF
        write(ifl,'(a)')ffrTseries_FFGF
        write(ifl,'(a)')'FFR Time Series Header field'
        write(ifl,*)0
        write(ifl,'(a)')tseriesHead_FFGF
        write(ifl,*)1
!        write(ifl,*) Nforce+1   !  numofFForceWall+1
        write(ifl,*) Nforce   !  numofFForceWall+1
        do ic=1,Nforce
          write(ifl,'(i4,a,i4,1x,a)')ic,' f ',13,NM_F(ic)
        end do
        write(ifl,'(a)')tseriesHeade_FFGF

        write(ifl,'(a)')customData_FFGF
        write(ifl,'(a)')ffrAxis_FFGF
        write(ifl,'(a)')'FFR Axis data field'
        write(ifl,*)0
        write(ifl,'(a)')axisData_FFGF
        write(ifl,*)1
        write(ifl,*)Nforce
        do ic=1,Nforce
          write(ifl,'(i4,a,6(e16.9,1x))')ic,' V ',
     &                    begin_mnt(1:3,ic),end_mnt(1:3,ic)
        end do
        write(ifl,'(a)')axisDatae_FFGF
        
        write(ifl,'(a)')customData_FFGF
        write(ifl,'(a)')ffrText_FFGF
        write(ifl,'(a)')'FFR Text data field'
        write(ifl,*)0
        write(ifl,'(a)')textData_FFGF
        write(ifl,*)1
        write(ifl,'(i4,1x,a)')1,'Columm Data'
        write(ifl,'(i4,1x,a)')1,'   1   Time'
        write(ifl,'(i4,1x,a)')1,'   2   Press_x'
        write(ifl,'(i4,1x,a)')1,'   3   Press_z'
        write(ifl,'(i4,1x,a)')1,'   4   Press_z'
        write(ifl,'(i4,1x,a)')1,'   5   Friction_x'
        write(ifl,'(i4,1x,a)')1,'   6   Friction_y'
        write(ifl,'(i4,1x,a)')1,'   7   Friction_z'
        write(ifl,'(i4,1x,a)')1,'   8   Total_x'
        write(ifl,'(i4,1x,a)')1,'   9   Total_y'
        write(ifl,'(i4,1x,a)')1,'  10   Total_z'
        write(ifl,'(i4,1x,a)')1,'  11   Moment_x'
        write(ifl,'(i4,1x,a)')1,'  12   Moment_y'
        write(ifl,'(i4,1x,a)')1,'  13   Moment_z'
        write(ifl,'(i4,1x,a)')1,'  14   Moment_axis'
        write(ifl,'(i4,1x,a)')-1,'Text data end.'
        write(ifl,'(a)')textDatae_FFGF
        close(ifl)
      else if(fforceFiletype==2) then
! ! --- Space separated text data. ------------------
!  5050 format(1X,'# num_wall= ',I4)
!         write(formatKeytmp,*)numofFForceWall+numofFForceBody
!         fnam1=trim(fnam) // '.shear'
!         open(ifl,file=fnam1,form='formatted',
!      &             status='unknown',iostat=ios)
!         write(ifl,*) '##### WALL SHEAR FORCE HISTORY ####'
!         write(ifl,5050) numofFForceWall+numofFForceBody
!         write(formatKey,*)'(a3,i4,a6,4x,i4,a6,7x,',
!      &          trim(adjustl(formatKeytmp)),'(i4,a1,a45,1x))'
!         write(ifl,formatKey) '# ',1,':time',2,':step',
!      &       ((ic-1)*3+3,':',nameofFForceWall(ic),ic=1,numofFForceWall),
!      &       ((numofFForceWall+ic-1)*3+3,':',nameofFForceBody(ic),
!      &                                          ic=1,numofFForceBody)
!         close(ifl)
!         fnam1=trim(fnam) // '.press'
!         open(ifl,file=fnam1,form='formatted',
!      &             status='unknown',iostat=ios)
!         write(ifl,*) '##### WALL PRESS FORCE HISTORY ####'
!         write(ifl,5050) numofFForceWall+numofFForceBody
!         write(formatKey,*)'(a3,i4,a6,4x,i4,a6,7x,',
!      &          trim(adjustl(formatKeytmp)),'(i4,a1,a45,1x))'
!         write(ifl,formatKey) '# ',1,':time',2,':step',
!      &       ((ic-1)*3+3,':',nameofFForceWall(ic),ic=1,numofFForceWall),
!      &       ((numofFForceWall+ic-1)*3+3,':',nameofFForceBody(ic),
!      &                                          ic=1,numofFForceBody)
!         close(ifl)
!         fnam1=trim(fnam) // '.total'
!         open(ifl,file=fnam1,form='formatted',
!      &             status='unknown',iostat=ios)
!         write(ifl,*) '##### WALL TOTAL FORCE HISTORY ####'
!         write(ifl,5050) numofFForceWall+numofFForceBody
!         write(formatKey,*)'(a3,i4,a6,4x,i4,a6,7x,',
!      &          trim(adjustl(formatKeytmp)),'(i4,a1,a45,1x))'
!         write(ifl,formatKey) '# ',1,':time',2,':step',
!      &       ((ic-1)*3+3,':',nameofFForceWall(ic),ic=1,numofFForceWall),
!      &       ((numofFForceWall+ic-1)*3+3,':',nameofFForceBody(ic),
!      &                                          ic=1,numofFForceBody)
!         close(ifl)
!         fnam1=trim(fnam) // '.moment'
!         open(ifl,file=fnam1,form='formatted',
!      &             status='unknown',iostat=ios)
!         write(ifl,*) '##### MOMENT HISTORY ####'
!         write(ifl,5050) numofFForceWall+numofFForceBody
!         write(formatKey,*)'(a3,i4,a6,4x,i4,a6,7x,',
!      &          trim(adjustl(formatKeytmp)),'(i4,a1,a45,1x))'
!         write(ifl,formatKey) '# ',1,':time',2,':step',
!      &       ((ic-1)*3+3,':',nameofFForceWall(ic),ic=1,numofFForceWall),
!      &       ((numofFForceWall+ic-1)*3+3,':',nameofFForceBody(ic),
!      &                                          ic=1,numofFForceBody)
!         close(ifl)
! !---------------------------------------------------
      end if
      statRtrn=0
      end subroutine ffr_ff_headerout
      
      end subroutine write_fluid_force
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
