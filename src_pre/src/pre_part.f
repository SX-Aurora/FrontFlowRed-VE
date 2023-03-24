!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine part
     &(mcell,mvrtx,mface,medge,mssfbc,
     & ncell,ncelb,nvrtx,nface,nedge,nssfbc,ncomp,ncomp_suf,npotn,
     & nphase,nrans,NMAT,MXMAT,NBCWAL,
     & NCV,NCVFAC,NALLCV,IBCCYC,IBCTCI,IBCSLD,IBCINT,IEMAX,
     & lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,lbface,listpr,
     & LVEDGE,LEFACE,lacell,
     & LBCSSF,LCYCSF,LMAT,locmsh,listbc,
     & cord,area,volume,gface,gcell,
     & SFAREA,SFCENT,CVVOLM,CVCENT,
     & time,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_hpc_input,only : hpc_control => inputdata
      use module_hpc,only : NPETOT
      use module_io,only  : ifli,ifll,ifle
      use module_rad
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mcell,mvrtx,mface,medge,mssfbc
      integer,intent(in)     :: ncell,ncelb,nvrtx,nface,nedge,nssfbc
      integer,intent(in)     :: ncomp,nrans,NMAT,MXMAT,ncomp_suf,nphase
      integer,intent(in)     :: NCV,NCVFAC,NALLCV,IBCCYC,IEMAX,IBCTCI,
     &                          IBCSLD,IBCINT,NBCWAL,npotn
      integer,intent(out)    :: ierror
!
      integer,intent(inout)  :: lvcell    (8,mcell)
      integer,intent(inout)  :: lfcell    (7,mcell)
      integer,intent(inout)  :: lvface    (4,mface)
      integer,intent(inout)  :: lcface    (2,mface)
      integer,intent(inout)  :: LVRTCV    (  mvrtx)
      integer,intent(inout)  :: LCVFAC    (4,mface)
!
      integer,intent(inout)  :: lbface    (2,mface)
      integer,intent(inout)  :: listpr  (0:3,mvrtx)
      integer,intent(inout)  :: LVEDGE    (2,medge)
      integer,intent(inout)  :: LEFACE    (5,mface)
      integer,intent(inout)  :: lacell    (  mcell)
      integer,intent(inout)  :: LBCSSF    ( mssfbc)
      integer,intent(inout)  :: LCYCSF    ( mssfbc)
      integer,intent(inout)  :: locmsh    ( mssfbc)
      integer,intent(inout)  :: LMAT      (  mvrtx)
      integer,intent(inout)  :: listbc    (4,mssfbc)
!
      real*8 ,intent(inout)  :: cord  (3,mvrtx)
      real*8 ,intent(inout)  :: area  (4,mface)
      real*8 ,intent(inout)  :: volume(  mcell)
      real*8 ,intent(inout)  :: gface (3,mface)
      real*8 ,intent(inout)  :: gcell (3,mcell)
      real*8 ,intent(inout)  :: SFAREA(4,medge)
      real*8 ,intent(inout)  :: SFCENT(3,medge)
      real*8 ,intent(inout)  :: CVVOLM(  mvrtx)
      real*8 ,intent(inout)  :: CVCENT(3,mvrtx)
      real*8 ,intent(in)     :: time
!
! --- [local entities]
!
      integer :: icpu,ierr
      integer :: ncomp_h,nrans_h,ncomp_suf_h,nphase_h,npotn_h
      integer :: ncell_h,ncelb_h,nvrtx_h,nface_h,nedge_h,NCV_h
      integer :: NCVFAC_h,NALLCV_h,NCVIN_H,IEMAX_h,nssfbc_h,nssfbc_oldh
      integer :: NBCINL_h,NBCWAL_h,NBCOUT_h,NBCSYM_h,MAXWAL
      integer :: NBCTCI_h,IBCTCI_h
      integer :: NBCCYC_h,IBCCYC_h
      integer :: NBCSLD_h,IBCSLD_h
      integer :: NBCINT_h,IBCINT_h
      integer :: NMAT_h,IVBCTOT_h
      
      integer :: ncelb_local
      integer :: rewall=0
      integer :: NBCpair=0,imode=0
!
      ierr=0
      IF(NPETOT.eq.1) return
      MAXWAL=NBCWAL
!--------------------------------------------
! --- Read HPC control file: CONT.HPC
!--------------------------------------------
      
      write(ifll,5000) 
      imode=0
      IF(radflag.EQ.1.AND.radmodel(1:3).NE.'FVM') imode=1
      write (*,'(/,a)') '    ###    READ HPC CONTROL FILE ...'
      call hpc_control(ifli,ifll,ifle,NPETOT,ierr)
      if(ierr.ne.0) goto 9999
!--------------------------------------------
! --- Read Serial grid and BC data for HPC
!--------------------------------------------
      write(ifll,5000) 
      write (*,'(/,a)') '    ###    SERIAL CPU GRID/BC FILE ...'
!
      call DATA
     & (ncell,nvrtx,nface,nedge,nssfbc,NCV,NCVFAC,NALLCV,
     &  mvrtx,mcell,medge,mssfbc,
     &  cord,lvcell,lacell) 
!--------------------------------------------
! --- 
!--------------------------------------------
      write(ifll,5000) 
      write (ifll,'(/,a)') '    ###    PARTITIONING DOMAIN ...'
      call PARASET(medge,nedge,mssfbc,IEMAX,
     &           LBCSSF,LCYCSF,LVEDGE)
!--------------------------------------------
! --- 
!--------------------------------------------
      write(ifll,5000) 
      write(ifll,'(/,a)')
     &'    ###    MAKE CYCLIC BOUNDARY CONNECTIVITY ...'
      if(BC_CYCL_tot>0.or.BC_TCIN_tot>0.or.BC_SLID_tot>0.or.
     &   BC_INTR_tot>0)
     & call CYCLE
     & (ncell,nvrtx,nface,nedge,nssfbc,NCV,NCVFAC,NALLCV,IBCCYC,IBCTCI,
     & IBCSLD,IBCINT,
     & medge,mssfbc,mface,mcell,
     & lbface,lcface,lvface)
!--------------------------------------------
! --- 
!--------------------------------------------
      write (ifll,5000) 
      write (ifll,'(/,a)') 
     &'    ###    MAKE SHORT DISTANCE CELL CONNECTIVITY ...'
      write(ifll,*)
!???      call WALL(mvrtx,cord)
!--------------------------------------------
! --- 
!--------------------------------------------
      write(ifll,5000) 
      write(ifll,'(/,a)')
     & '    ###    MAKE COMMUNICATION AND WRITE DOWN ...'
      call COMM(mcell,mface,nface,lvcell,lvface,lcface,lbface)
!--------------------------------------------
! --- 
!--------------------------------------------
      write(ifll,5000) 
      write (ifll,'(/,a)') 
     & '    ###    OUTPUT HPC MESH/BC FILE ...'
      call OUTPUT_hpc_grid
     &  (nface,ncelb,NCV,nedge,nssfbc,NCVFAC,NALLCV,
     &   ncell,mcell,mface,medge,mvrtx,
     &   lfcell,lcface,LEFACE,lvface,cord,lvcell,lacell,
     &   lbface)
!
      write (ifll,'(/,a)') '    ###    PARTITION  SUCCESSFULLY'
      write(ifll,5000)
!--------------------------------------------
! --- deallocate BC arraies for globle mesh
!--------------------------------------------
      write (ifll,'(/,a)') '    ###    TESTING HPC DOMAINS'
      if(BC_INLT_tot.gt.0) then
        deallocate(BC_INLT)
        deallocate(BC_IV3D)
      endif
      if(BC_WALL_tot.gt.0) deallocate(BC_WALL)
      if(BC_SYMT_tot.gt.0) deallocate(BC_SYMT)
      if(BC_CYCL_tot.gt.0) deallocate(BC_CYCL)
      if(BC_INTR_tot.gt.0) deallocate(BC_INTR)
      if(BC_BODY_tot.gt.0) deallocate(BC_BODY)
      if(BC_FREE_tot.gt.0) deallocate(BC_FREE)
      if(BC_TCIN_tot.gt.0) deallocate(BC_TCIN)
      if(BC_MWAL_tot.gt.0) then
        deallocate(BC_MWAL)
        deallocate(BC_MV3D)
      endif
!
! --- 
!
      deallocate(ICELTYP)
      deallocate(NODE_ID)
!      deallocate(WALLnode)
!      deallocate(WALLDIST)
      deallocate(NEIBPE)
      deallocate(NEIBPEinv)
      deallocate(NEIBPETOT)
!
      if(isufCYCLtot.gt.0) then
        deallocate(SUFtoCELL)
        deallocate(iSUF_NODE)
        deallocate(ISUF_NODE_P)
      endif
!
      deallocate(ELMstack)
      deallocate(ELMitem)
      deallocate(NODstack)
      deallocate(NODitem)
      deallocate(ELMtotL)
!
!--------------------------------------------
! --- Calculate each cpu's dimension size
!--------------------------------------------
!
      do 200 icpu=1,PETOT
!
      BC_INLT_tot=0
      BC_WALL_tot=0
      BC_SYMT_tot=0
      BC_CYCL_tot=0
      BC_INTR_tot=0
      BC_BODY_tot=0
      BC_FREE_tot=0
      BC_MWAL_tot=0
      BC_TCIN_tot=0
!
! --- 
!
      nvrtx_h=nvrtx_hpc(icpu)
      ncell_h=ncell_hpc(icpu)
      nface_h=nface_hpc(icpu)
      ncelb_h=ncelb_hpc(icpu)
      NCV_h=NCV_hpc(icpu)
      nedge_h=nedge_hpc(icpu)
      NCVIN_h=NODtotL(icpu)
!
! --- 
!
      write (ifll,5000)
      write (ifll,5010) icpu,icpu-1
      write (ifll,5000)
      write(ifle,*) '    ###    NCVIN_H= ',NODTOTL(ICPU)
      WRITE(IFLE,*) '    ###    NCV_H  = ',NCV_H
      WRITE(IFLE,*) '    ###    NVRTX_H= ',NVRTX_H
      WRITE(IFLE,*) '    ###    NCELL_H= ',NCELL_H
!
      NALLCV_h=0
      NCVFAC_h=0
      nssfbc_h=0
      nssfbc_oldh=0
      IEMAX_h=0
      IBCCYC_h=0
      IBCTCI_h=0
      NMAT_h=0
!--------------------------------
! --- Read HPC mesh and BC file
!--------------------------------
      lvcell=0
      lfcell=0
      lvface=0
      lcface=0
      LVRTCV=0
      LCVFAC=0
      lbface=0
      listpr=0
      LVEDGE=0
      LEFACE=0
      lacell=0
      LBCSSF=0
      LCYCSF=0
      LMAT=0
!
      cord=0.D0
      area=0.D0
      volume=0.D0
      gface=0.D0
      gcell=0.D0
      SFAREA=0.D0
      SFCENT=0.D0
      CVVOLM=0.D0
      CVCENT=0.D0
!--------------------------------
! --- Read HPC mesh (grid) data
!--------------------------------
      write(ifll,*) '    ###    read_hpc_grid ...'
      call read_hpc_grid !jiang
     &  (imode,icpu,mcell,mvrtx,mface,nvrtx_h,ncell_h,NBCpair,
     &   cord,lacell,lvcell,LEFACE,ierror)
!
      if( ierror.ne.0 ) goto 9999
!--------------------------------------------
! --- nface_h;ncelb_h,NCV_h 
!--------------------------------------------
      write(ifll,*) '    ###    LIST_FACE_HPC ...'
      call list_face_hpc
     & (mcell,mface,mvrtx,nvrtx_h,ncell_h,nface_h,ncelb_h,NCV_h,
     &  ncelb_local,
     &  lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,ierr)
        ncelb_hpc(icpu)=ncelb_local
        ncelb_h=ncelb_local
      if(ierr.ne.0) goto 9999
!--------------------------------------------
! --- IBCCYC_h;
!--------------------------------------------
      write(ifll,*) '    ###    LIST_FACEPR_HPC ...'
      call list_facepr_hpc
     &  (mcell,mface,nface_h,mvrtx,nvrtx_h,ncell_h,mssfbc,
     &   IBCCYC_h,NBCpair,
     &   lvface,lcface,lbface,lacell,listpr,LEFACE,ierr)
      LEFACE=0
      if( ierr.ne.0 ) goto 9999
!--------------------------------------------
! --- 
!--------------------------------------------
      write(ifll,*) '    ###    LIST_FACEBC_HPC ...'
      call list_facebc_hpc
     & (mface,nface_h,mvrtx,nvrtx_h,ncell_h,
     &  lvface,lbface,lcface)
!--------------------------------------------
! --- NCV_h
!--------------------------------------------
!      call list_facein_hpc
!     & (mcell,mface,medge,mvrtx,ncell_h,nface_h,ncelb_h,nvrtx_h,NCV_h,
!     &  lacell,lfcell,lvface,lcface,lbface,
!     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierr)
!      if(ierr.ne.0) goto 9999
!--------------------------------------------
! --- NCV_h
!--------------------------------------------
!      call list_facenn_hpc
!     & (mcell,mface,medge,mvrtx,ncell_h,nface_h,ncelb_h,nvrtx_h,NCV_h,
!     &  lfcell,lvface,lcface,lbface,
!     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierr)
!      if( ierr.ne.0 ) goto 9999
!
!--------------------------------------------
! --- nedge_h
!--------------------------------------------
      write(ifll,*) '    ###    LIST_EDGE_HPC ...'
      call list_edge_hpc
     & (mcell,mface,medge,mvrtx,
     &  ncell_h,nvrtx_h,nface_h,nedge_h,ncelb_h,NCV_h,
     &  lvcell,lfcell,lvface,lcface,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierr)
      if(ierr.ne.0) goto 9999
!--------------------------------------------
! --- lacell(mcell)
!--------------------------------------------
      write(ifll,*) '    ###    LIST_FACECHK_HPC ...'
      call list_facechk_hpc
     & (mface,mcell,nface_h,ncell_h,mvrtx,nvrtx_h,
     &  lvface,lcface,lbface,lacell,ierr)
      if(ierr.ne.0) goto 9999
!--------------------------------------------
! --- 
!--------------------------------------------
      write(ifll,*) '    ###    BC_CHECK_HPC ...'
      call bc_check_hpc
     & (mface,mcell,nface_h,mvrtx,nvrtx_h,
     &  lbface,lcface,lvface,lacell,ierr)
      if(ierr.ne.0) goto 9999
!--------------------------------------------
! --- lacell(mcell)
!--------------------------------------------
      write(ifll,*) '    ###    LIST_SOLID_HPC ...'
      call list_solid_hpc(mcell,ncell_h,ncelb_h,lacell,ierr)
      if(ierr.ne.0) goto 9999
!--------------------------------------------
! --- 
!--------------------------------------------
      write(ifll,*) '    ###    LIST_FLUID_HPC ...'
      call list_fluid_hpc
     & (mcell,mface,nface_h,ncell_h,ncelb_h,mvrtx,nvrtx_h,
     &  lbface,lcface,lacell,ierr)
      if(ierr.ne.0) goto 9999
!
!--------------------------------------------
! --- LMAT(mvrtx);LBCSSF(mssfbc);LCYCSF(mssfbc)
! --- nssfbc,NCVFAC,NALLCV,IEMAX,NMAT
!--------------------------------------------
!
      write(ifll,*) '    ###    METRIC_ADMIN ...'
      call metric_admin
     & (mcell,mvrtx,mface,medge,
     &  ncell_h,nvrtx_h,nface_h,nedge_h,nssfbc_h,nssfbc_oldh,NCV_h,
     &  mssfbc,NCVFAC_h,NALLCV_h,IEMAX_h,NMAT_h,ncelb_h,
     &  lcface,lvface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,LBCSSF,LCYCSF,listbc,listpr,LMAT,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,lvcell,lfcell,
     &  time,-1.d0,ierr,1)
        if(ierr.ne.0) goto 9999
!
      NALLCV_hpc(icpu)=NALLCV_h
      NCVFAC_hpc(icpu)=NCVFAC_h
      nssfbc_hpc(icpu)=nssfbc_h
      IEMAX_hpc(icpu)=IEMAX_h
!
!--------------------------------------------
! --- NBCINL,NBCWAL,NBCOUT,NBCCYC,NBCSYM
!--------------------------------------------
      NBCINL_h=0
      NBCWAL_h=0
      NBCOUT_h=0
      NBCCYC_h=0
      NBCSYM_h=0
      NBCTCI_h=0
      NBCSLD_h=0
      NBCINT_h=0
!
      call list_BC_hpc
     & (mssfbc,nssfbc_h,NBCINL_h,NBCWAL_h,NBCOUT_h,NBCCYC_h,NBCSYM_h,
     &  NBCTCI_h,NBCSLD_h,NBCINT_h,LBCSSF)
      NBCINL_hpc(icpu)=NBCINL_h
      NBCWAL_hpc(icpu)=NBCWAL_h
      NBCCYC_hpc(icpu)=NBCCYC_h
      NBCOUT_hpc(icpu)=NBCOUT_h
      NBCSYM_hpc(icpu)=NBCSYM_h
      NBCTCI_hpc(icpu)=NBCTCI_h
      NBCSLD_hpc(icpu)=NBCSLD_h
      NBCINT_hpc(icpu)=NBCINT_h
!
!-------------------------------------------------------------
! --- wall distance
!-------------------------------------------------------------
! 
      write(ifll,'(a)') '###    list_wall_distance...'
      call list_wall_distance
!     &      (icpu,rewall,NCV_h,NBCWAL_h,nssfbc_h,NEDGE_h,
     &      (icpu,rewall,NCV_h,MAXWAL,nssfbc_h,NEDGE_h,
     &       mvrtx,mssfbc,medge,
     &       LMAT,LBCSSF,SFCENT,CVCENT,SFAREA)
!
!
! jiang
      IF(radflag.EQ.1.AND.radmodel(1:3).NE.'FVM') THEN
	write(ifll,'(a)') '###    list_radiation_data...'
	CALL advance_rad_hpc(icpu,mvrtx,mcell,nvrtx_h,ncell_h,
     &	cord,lacell,lvcell)
      ENDIF
!-------------------------------------------------------------
! --- Resort multi-material and Write down geom file for HPC 
!-------------------------------------------------------------
      write (ifll,'(/,a)') '    ###    list_output_geom ...' 
      call list_output_geom
     & (icpu,rewall,mcell,mvrtx,mface,medge,mssfbc,NBCWAL_h,
     &  NCV_h,nvrtx_h,ncell_h,NMAT_h,NCVFAC_h,
     &  NALLCV_h,nedge_h,nssfbc_h,nssfbc_oldh,
     &  NCVIN_h,
     &  LMAT,LVEDGE,LVRTCV,LBCSSF,LCYCSF,locmsh,listbc,listpr,
     &  SFAREA,SFCENT,CVVOLM,CVCENT
     &  )
      NMAT_hpc(icpu)=NMAT_h
!--------------------------------------------
! --- 
!--------------------------------------------
      ncomp_hpc(icpu)=ncomp
      ncomp_suf_hpc(icpu)=ncomp_suf
      nrans_hpc(icpu)=nrans
      nphase_hpc(icpu)=nphase
      npotn_hpc(icpu)=npotn
!
      write (ifll,5000)
!
 200  continue
!
!-------------------------------------------------------
! --- output dimen.parm file in hpc_0000=>hpc_0007 etc.
!-------------------------------------------------------
! --- max. dimension size
!-------------------------------------------------------
!
        ncomp_h=0
        ncomp_suf_h=0
        nphase_h=0
        npotn_h=0
        nrans_h=0
        IVBCTOT_h=0
        call dimsn_hpc(nvrtx_h,ncell_h,ncelb_h,nface_h,nedge_h,
     &  NCVIN_h,NCV_h,NALLCV_h,NCVFAC_h,
     &  nssfbc_h,ncomp_h,ncomp_suf_h,nphase_h,npotn_h,
     &  nrans_h,NMAT_h,IEMAX_h,IVBCTOT_h,NBCWAL_h,
     &  NBCSLD_h,NBCINT_h)
        write(ifll,5000)
!
      do 300 icpu=1,PETOT
        call wrtdmn
     & (nvrtx_h,ncell_h,ncelb_h,nface_h,nedge_h,
     &  NCVIN_h,NCV_h,NALLCV_h,NCVFAC_h,
     &  nssfbc_h,nssfbc_oldh,ncomp_h,ncomp_suf_h,nphase_h,npotn_h,
     &  nrans_h,MXMAT,
     &  IEMAX_h,IVBCTOT_h,NBCWAL_h,
     &  NBCSLD_h,NBCINT_h,NPART,icpu)
      write (ifll,5000)
!
 300  continue
!-------------------------------------------------------
! --- Reordering comm file because of multi-material
!-------------------------------------------------------
      write(ifll,*) '    ###    REORDERING COMM FILE ...'
!
      call re_comm(mssfbc)
!
!      write(ifll,*) '    ###    HPC: list_output_sld ......'
!      call list_output_sld
!     & (NPETOT,nssfbc,nedge,mssfbc,medge,LBCSSF,LCYCSF)
!-------------------------------------------------------
! --- deallocate 
!-------------------------------------------------------
!
      deallocate(GRIDout)
      deallocate(BCout)
      deallocate(COMMout)
      deallocate(WORKFIL)
      deallocate(GEOMout)
      deallocate(RODRICV)
      deallocate(hpcdmn)
      deallocate(NODtotL)
      deallocate(SLIDNG)
!---------------------------------
! --- dimension deallocable array
!---------------------------------
      deallocate(nvrtx_hpc)
      deallocate(ncell_hpc)
      deallocate(ncelb_hpc)
      deallocate(nface_hpc)
      deallocate(nedge_hpc)
      deallocate(NCV_hpc)
      deallocate(NALLCV_hpc)
      deallocate(NCVFAC_hpc)
      deallocate(nssfbc_hpc)
      deallocate(ncomp_hpc)
      deallocate(ncomp_suf_hpc)
      deallocate(npotn_hpc)
      deallocate(nphase_hpc)
      deallocate(nrans_hpc)
      deallocate(NMAT_hpc)
      deallocate(IEMAX_hpc)
      deallocate(NBCINL_hpc)
      deallocate(NBCWAL_hpc)
      deallocate(NBCCYC_hpc)
      deallocate(NBCOUT_hpc)
      deallocate(NBCSYM_hpc)
      deallocate(NBCTCI_hpc)
      deallocate(NBCSLD_hpc)
      deallocate(NBCINT_hpc)
      deallocate(IVBCTOT_HPC) 
!
      return
 5000 format(2X,'|',108('='),'|')
 5010 format(2X,'|',40(' '),'Matrix for cpu = ',I4,7X,
     &  'my_rank= ',I4,27(' '),'|')
 9999 continue
      write(ifle,*) '(part)'
!
      ierror=1
!--------------------------------------------
      end subroutine part
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dimsn_hpc(nvrtx,ncell,ncelb,nface,nedge,
     &  NCVIN,NCV,NALLCV,NCVFAC,
     &  nssfbc,ncomp,ncomp_suf,nphase,npotn,nrans,
     &  NMAT,IEMAX,IVBCTOT,NBCWAL,NBCSLD,NBCINT)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_hpc,only : NPETOT
      use module_io, only : ifli,ifll,ifle
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout)     :: 
     &  nvrtx,ncell,ncelb,nface,nedge,ncomp_suf,nphase,
     &  NCVIN,NCV,NALLCV,NCVFAC,nssfbc,ncomp,nrans,NMAT,IEMAX,IVBCTOT,
     &  NBCWAL,NBCSLD,NBCINT,npotn
!
! --- [local entities]
!
      integer :: icpu
!
      nvrtx=0
      ncell=0
      ncelb=0
      nface=0
      nedge=0
      NCVIN=0
      NCV=0
      NALLCV=0
      NCVFAC=0
      nssfbc=0
      ncomp=0
      ncomp_suf=0
      npotn=0
      nphase=0
      nrans=0
      NMAT=0
      IEMAX=0
      IVBCTOT=0
      do 100 icpu=1,PETOT
      nvrtx=MAX(nvrtx,nvrtx_hpc(icpu))
      ncell=MAX(ncell,ncell_hpc(icpu))
      ncelb=MAX(ncelb,ncelb_hpc(icpu))
      nface=MAX(nface,nface_hpc(icpu))
      nedge=MAX(nedge,nedge_hpc(icpu))
      NCVIN=MAX(NCVIN,NODtotL(icpu))
      NCV=MAX(NCV,NCV_hpc(icpu))
      NALLCV=MAX(NALLCV,NALLCV_hpc(icpu))
      NCVFAC=MAX(NCVFAC,NCVFAC_hpc(icpu))
      nssfbc=MAX(nssfbc,nssfbc_hpc(icpu))
      ncomp=MAX(ncomp,ncomp_hpc(icpu))
      ncomp_suf=MAX(ncomp_suf,ncomp_suf_hpc(icpu))
      npotn=MAX(npotn,npotn_hpc(icpu))
      nphase=MAX(nphase,nphase_hpc(icpu))
      nrans=MAX(nrans,nrans_hpc(icpu))
      NMAT=MAX(NMAT,NMAT_hpc(icpu))
      IEMAX=MAX(IEMAX,IEMAX_hpc(icpu))
      IVBCTOT=MAX(IVBCTOT,IVBCTOT_HPC(icpu))
      NBCWAL=MAX(NBCWAL,NBCWAL_hpc(icpu))
      NBCSLD=MAX(NBCSLD,NBCSLD_hpc(icpu))
      NBCINT=MAX(NBCINT,NBCINT_hpc(icpu))
 100  continue
!
      return
!
      end subroutine dimsn_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
     
