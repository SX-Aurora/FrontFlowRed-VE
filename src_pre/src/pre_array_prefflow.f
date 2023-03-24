!======================================================================!
!                                                                      !
! Software Name :                                                      !
!     FrontFlow/red   (Ver. 1.0) (Cell-Center,FACE-BASE)               !
!                                                                      !
!     Main Program : um0_FFRMAIN                                       !
!                                                                      !
!                          Written by Masayuki KAKEI,                  !
!                                     Huilai ZHANG, Osamu KITAMURA     !
!                                     Yoshinobu YAMADE, Masato IDA,    !
!                                     Takeshi UNEMURA, Eisuke YAMADA,  !
!                                     Nobuyuki TANIGUCHI               !
!                                     2003/03/25                       !
!                                                                      !
!     FontFlow/red    (Ver. 1.2) (EDGE-BASE,HPC)                       !
!                          Modified by                                 !
!                                     Huilai ZHANG                     !
!                                     Takeshi UNEMURA                  !
!                                     NAKAJIMA???                      !
!                                                                      !
!     Contact address : The University of Tokyo, FSIS Project          !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 2.0)                           !
!                                                                      !
! Main Program  : FRONTFLOW    (Vertex-Center,EDGE-BASE,HPC,           !
!                                  Multi-Domain)                       !
!                                                                      !
!                          Written by Huilai ZHANG                     !
!                                     Takeshi UNEMURA                  !
!                                     Nobuyuki TANIGUCHI               !
!                                     2003/03/25                       !
!                                                                      !
!     Contact address : The University of Tokyo, FSIS Project          !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 3.0)                           !
!                                                                      !
! Main Program  : FRONTFLOW    (Vertex-Center,EDGE-BASE,HPC,           !
!                                  Multi-Domain)                       !
!                                                                      !
!                          Written by Huilai ZHANG                     !
!                                     Takeshi UNEMURA                  !
!                                     Nobuyuki TANIGUCHI               !
!                                     2007/12/03                       !
!                                                                      !
!     Contact address : The University of Tokyo, FSIS Project          !
!                                                                      !
!======================================================================!

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      MODULE MODULE_ARRAY_PreFFLOW
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- Cell connectivity
!
      integer,allocatable :: lvcell    (:,:)
      integer,allocatable :: lacell    (  :)
      integer,allocatable :: lfcell    (:,:)
      integer,allocatable :: lvface    (:,:)
      integer,allocatable :: lcface    (:,:)
      integer,allocatable :: lbface    (:,:)
!
! --- Cell Matrix 
!
      real*8,allocatable  :: cord      (:,:)
!
! --- CV connectivity
!
      integer,allocatable :: LEFACE    (:,:)
      integer,allocatable :: LCVFAC    (:,:)
      integer,allocatable :: LVEDGE    (:,:)
      integer,allocatable :: LVRTCV    (  :)
      integer,allocatable :: LBCSSF    (  :)
      integer,allocatable :: LCYCSF    (  :)
      integer,allocatable :: locmsh    (  :)
      integer,allocatable :: listpr    (:,:)
      integer,allocatable :: listbc    (:,:)
      integer,allocatable :: LMAT      (  :)
!
      END MODULE MODULE_ARRAY_PreFFLOW
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_partitioner
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
        use module_hpc_input,
     &                 only : NPART,UCDFLAG,WALLFLAG,
     &                        GRIDFIL,BCFIL,INPFIL,UCDFIL,
     &                        WALFIL,DIRHED,
     &                        HEADERg, HEADERb, HEADERc, HEADERw
        implicit REAL*8 (A-H,O-Z)
!
! --- GLOBAL
!
        integer             :: INODTOT,ICELTOT,IEDGTOT
        integer             :: PETOT
        integer             :: isufCYCLtot,isufPAIRtot
        
! --- FILE
        character (len=80)  :: METISFIL
        character (len=80)  :: HEADER1,HEADER2
        character*80        :: LINE
! --- File name:
        character (len=80), allocatable :: GRIDout(:)
        character (len=80), allocatable :: BCout(:)
        character (len=80), allocatable :: COMMout(:)
        character (len=80), allocatable :: WORKFIL(:)
        character (len=80), allocatable :: hpcdmn(:)
        character (len=80), allocatable :: GEOMout(:)
        character (len=80), allocatable :: RODRICV(:)
        character (len=80), allocatable :: walldis(:)
        character (len=80), allocatable :: SLIDNG(:)
!
        integer,allocatable :: ICELTYP(:)
!
        character (len=80), allocatable ::   RADout(:) !Jiang
!
! --- PARALLEL
!
        integer,allocatable :: ISUF_NODE(:,:)
        integer,allocatable :: SUFtoCELL(:,:)
        integer,allocatable :: ISUF_NODE_P(:,:,:)
!
        integer,allocatable :: ISUF_NODE_pair(:,:)
        integer,allocatable :: SUFtoCELL_pair(:,:)
        integer,allocatable :: ISUF_NODE_O(:,:,:)
!
        integer,allocatable :: NODE_ID(:,:)
!
        integer,allocatable :: WALLnode(:)
        real*8 ,allocatable :: WALLDIST(:)
!
        integer,allocatable :: NEIBPE(:,:)
        integer,allocatable :: NEIBPEinv(:,:)
        integer,allocatable :: NEIBPETOT(:)
!
        integer,allocatable :: IMPORT_index(:)
        integer,allocatable :: IMPORT_item(:)
!
        integer,allocatable :: EXPORT_index(:)
        integer,allocatable :: EXPORT_item(:)
!
        integer,allocatable :: ELMstack(:)
        integer,allocatable :: ELMitem(:)
!
        integer,allocatable :: NODstack(:)
        integer,allocatable :: NODitem(:)
!
        integer,allocatable :: NODstackN(:)
        integer,allocatable :: NODitemN(:)
!
        integer,allocatable :: NODtotL(:)
        integer,allocatable :: NODtotX(:)
        integer,allocatable :: NODtotW(:)
!
        integer,allocatable :: ELMtotL(:)
!-----------
! --- WORK
!-----------
        integer,allocatable :: WK1(:), WK2(:), WK3(:), WK4(:)
        integer,allocatable :: WK5(:), WK6(:),IW1(:,:),IW2(:,:)
        integer,allocatable :: givbcnd(:),glvbcnd(:),pMap(:)
        logical,allocatable :: tempFlag(:)
!
! --- BOUNDARYs
!
        integer,allocatable :: BC_INLT(:)
        integer,allocatable :: BC_WALL(:)
        integer,allocatable :: BC_SYMT(:)
        integer,allocatable :: BC_CYCL(:,:)
        integer,allocatable :: BC_INTR(:,:)
        integer,allocatable :: BC_BODY(:)
        integer,allocatable :: BC_FREE(:)
        integer,allocatable :: BC_MWAL(:)
        integer,allocatable :: BC_TCIN(:)
!
        real*8,allocatable  :: BC_IV3D(:,:), BC_MV3D(:,:)
!
        integer             :: BC_INLT_tot, BC_WALL_tot, BC_FREE_tot
        integer             :: BC_SYMT_tot, BC_CYCL_tot, BC_BODY_tot
        integer             :: BC_MWAL_tot, BC_TCIN_tot, BC_SLID_tot
        integer             :: BC_INTR_tot
        integer             :: NSUF_NODE_P,NSUF_NODE_O
!
! --- dimension
!
        integer,allocatable :: nvrtx_hpc(:)
        integer,allocatable :: ncell_hpc(:)
        integer,allocatable :: ncelb_hpc(:)
        integer,allocatable :: nface_hpc(:)
        integer,allocatable :: nedge_hpc(:)
        integer,allocatable :: NCV_hpc(:)
        integer,allocatable :: NALLCV_hpc(:)
        integer,allocatable :: NCVFAC_hpc(:)
        integer,allocatable :: nssfbc_hpc(:)
        integer,allocatable :: ncomp_hpc(:)
        integer,allocatable :: ncomp_suf_hpc(:)
        integer,allocatable :: npotn_hpc(:)
        integer,allocatable :: nphase_hpc(:)
        integer,allocatable :: nrans_hpc(:)
        integer,allocatable :: NMAT_hpc(:)
        integer,allocatable :: IEMAX_hpc(:)
        integer,allocatable :: NBCINL_hpc(:)
        integer,allocatable :: NBCWAL_hpc(:)
        integer,allocatable :: NBCCYC_hpc(:)
        integer,allocatable :: NBCOUT_hpc(:)
        integer,allocatable :: NBCSYM_hpc(:)
        integer,allocatable :: NBCTCI_hpc(:)
        integer,allocatable :: NBCSLD_hpc(:)
        integer,allocatable :: NBCINT_hpc(:)
        integer,allocatable :: IVBCTOT_HPC(:)
!
      end module module_partitioner
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_cood
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      real*8,save,allocatable :: cordsq  (:,:,:)
      real*8,save,allocatable :: areasq  (:,:,:)
      real*8,save,allocatable :: volumesq(:,:)
      real*8,save,allocatable :: gfacesq (:,:,:)
      real*8,save,allocatable :: gcellsq (:,:,:)
!
! --- Cell Matrix 
!
      real*8,allocatable  :: volume    (  :)
      real*8,allocatable  :: gcell     (:,:)
      real*8,allocatable  :: area      (:,:)
      real*8,allocatable  :: gface     (:,:)
!
! --- CV Matrix
!
      real*8,allocatable  :: SFAREA    (:,:)
      real*8,allocatable  :: SFCENT    (:,:)
      real*8,allocatable  :: CVVOLM    (  :)
      real*8,allocatable  :: CVCENT    (:,:)
!
      end module module_cood
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine pre_namelist_admin
     & (mcell,mvrtx,mface,medge,mssfbc,ncell,nvrtx,ncomp,nrans,
     &  ncomp_suf,nphase,npotn,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io       ,only : ifli,ifll,ifle,openfiles,cntlnam,
     &                            gdformat
      use module_io       ,only : input_io       =>inputdata
      use module_model    ,only : input_model    =>inputdata
      use module_flags    ,only : input_flags    =>inputdata
      use module_time     ,only : input_time     =>inputdata
      use module_deltat   ,only : input_deltat   =>inputdata
      use module_simple   ,only : input_simple   =>inputdata
      use module_gravity  ,only : input_gravity  =>inputdata
      use module_cgsolver ,only : input_cgsolver =>inputdata
      use module_rans     ,only : input_ransmodel  =>inputke
      use module_les      ,only : input_les      =>inputdata
      use module_species  ,only : input_species_gas  =>inputdata_gas
      use module_species  ,only : input_species_suf  =>inputdata_suf
      use module_chemcntl ,only : input_chemcntl =>inputdata
      use module_chemreac ,only : input_chemreac =>inputdata
      use module_boundary ,only : input_boundary =>inputdata
      use module_source   ,only : input_source   =>inputdata
      USE module_material ,ONLY : INPUT_MAT      =>INPUTDATA
      use module_movegrid ,only : input_movegrid =>inputdata
      use module_output   ,only : input_output   =>inputdata
      use module_dimnsn   ,only : input_dimnsn   =>inputdata
      use module_debug    ,only : input_debug    =>inputdata
      use module_usersub  ,only : input_user     =>inputdata
      use module_hpc      ,only : input_hpc      =>inputdata
      use module_flow_sound,only: INPUT_sound    =>inputdata
      use module_anim     ,only : input_anim     =>inputdata
      USE MODULE_BOUNDARY ,ONLY : nbcnd,boundName,kdbcnd,kxnone,
     &                            kdfire,kdbuff,kdcvd,kdshutr,kdpors
      use module_Euler2ph ,ONLY : input_eul2ph   =>inputdata
      use module_vof      ,ONLY : input_vof      =>inputdata
      use module_scalar   ,ONLY : input_scalar   =>inputdata
      use module_param    ,ONLY : smallNrm
      use module_Euler2ph ,ONLY : NPHS
      use module_rad   ,ONLY : input_rad   =>inputdata  !jiang
      USE MODULE_BOUNDARY ,ONLY : lsldbc,lovstbc
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: mcell,mvrtx,mface,medge,mssfbc
      integer,intent(inout) :: ncell,nvrtx,ncomp,nrans,ncomp_suf,
     &                         nphase,npotn
      integer,intent(inout) :: ierror
!
! --- [local entities]
!
      logical :: lcomp,lrans,E2P=.false.,kemdl,RSMmdl,
     &           lvof=.false.,
     &           lsld=.false.,
     &           lowke=.false.,
     &           lsuf=.false.,
     &           lRNG=.false.,
     &           lCHEN=.false.,
     &           l2SKE=.false.,
     &           lSDES=.false.,
     &           lpotn=.false.,
     &           RANSMDL=.false.,
     &           lKLES=.false.
      integer :: nbcndx=0,idens=0
      integer :: lFC=0

!
!-< file names & open files >-
!
      call input_io(smallNrm,ierror)
      if(ierror.ne.0) goto 9999
      call openfiles(ierror)
      if(ierror.ne.0) goto 9999
!
!-< model flags : defining [nrans]>-
!
      call input_model
     & (ifli,ifll,ifle,lcomp,lrans,kemdl,lowke,lRNG,lCHEN,l2SKE,
     &  RANSMDL,
     &  lSDES,lFC,lKLES,
     &  RSMmdl,nrans,cntlnam,idens,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< array sizes >-
!
      call nml_sizes(ifli,ifll,ifle,cntlnam,
     & mcell,mvrtx,mface,medge,mssfbc,
     & ncell,nvrtx,ncomp,ncomp_suf,lFC,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< Euler two Phase flow 
!
      call input_eul2ph(ifli,ifll,ifle,cntlnam,E2P,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< VOF two Phase flow
!
      call input_vof(ifli,ifll,ifle,cntlnam,lvof,ierror)
!
!-< Scalar equation : defining [nrans]>
!
      call input_scalar
     &   (ifli,ifll,ifle,cntlnam,
     &    lcomp,lrans,E2P,kemdl,lowke,lRNG,lCHEN,l2SKE,lSDES,lKLES,
     &    RSMmdl,RANSMDL,lvof,lFC,npotn,
     &    nrans,NPHS,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< data for every material >-
!
      call input_MAT(ifli,ifll,ifle,cntlnam,lsld,lrans,lovstbc,
     &               ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< some flags >-
!
      call input_flags(ifli,ifll,ifle,cntlnam,lrans,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< start/end of time integrarion >-
!
      call input_time(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< data to estimate time increment >-
!
      call input_deltat(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< control data for SIMPLE >-
!
      call input_simple(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< data for gravity >-
!
      call input_gravity(ifli,ifll,ifle,cntlnam,lcomp,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< parameter for cg solver >-
!
      call input_cgsolver(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< parameter of k-e model >-
!
      call input_ransmodel(ifli,ifll,ifle,cntlnam,nrans,kemdl,
     & lowke,lRNG,lCHEN,l2SKE,lKLES,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< species data : ncomp >- 
!
      call input_species_gas
     &        (ifli,ifll,ifle,cntlnam,ncomp,E2P,lvof,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< species data : ncomp_suf >- 
!
      call input_species_suf
     &        (ifli,ifll,ifle,cntlnam,ncomp_suf,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< parameter of LES model >-
!
      call input_les(ifli,ifll,ifle,cntlnam,nrans,ncomp,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< chemical reaction >-
!
      call input_chemcntl(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
      call input_chemreac(ifli,ifll,ifle,cntlnam,lsuf,
     &    ncomp,ncomp_suf,lcomp,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< boundary condition >-
!
      call input_boundary(ifli,ifll,ifle,cntlnam,lsuf,idens,lcomp,
     &      ncomp,nrans,nbcndx,nphase,lrans,E2P,lsld,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< source term >-
!
      call input_source(ifli,ifll,ifle,cntlnam,
     &         ncomp,nrans,lrans,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< control data for moving grid >-
!
      call input_movegrid(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< control data to output >-
!
      call input_output(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< dimension of computational domain >-
!
      call input_dimnsn(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< debug flags >-
!
      call input_debug(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< User subroutine >-
!
      call input_user(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< HPC in cntlnam=fort.1>-
!
      call input_hpc(ifli,ifll,ifle,cntlnam,nbcndx,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< flow sound>-
!
      call INPUT_sound
     & (ifli,ifll,ifle,cntlnam,
     &  nbcnd,boundName,kdbcnd,kxnone,kdbuff,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< animation >-
!
      call input_anim(ifli,ifll,ifle,cntlnam,nrans,ncomp,ierror)
      if( ierror.ne.0 ) goto 9999
!
!-< radiation >- Jiang
!
      call input_rad(ifli,ifll,ifle,cntlnam,ierror)
      return
!
 9999 continue
!
      write(ifle,'(a)') '(pre_namelist_admin)'
      ierror=1
      stop 'at pre_namelist_admin'
      end subroutine pre_namelist_admin
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine pre_main
     & (mcell,mvrtx,mface,medge,mssfbc,
     &  ncell,nvrtx,nvrtxnw,ncomp,ncomp_suf,nphase,nrans,npotn,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,      only : ifli,ifll,ifle,gdScale,
     &                           ffrgridOutFlag,ffgFormatKey,
     &                           gdformat
      use module_nowtime, only : time
      use module_time,    only : iters
      use module_hpc,     only : NPETOT
      use module_cood,    only : volume,gcell,area,gface
      use module_cood,    only : SFAREA,SFCENT,CVVOLM,CVCENT
      use module_model,only    : ical_mvmsh,vertex_cen,cell_cen,icon_cv
      use MODULE_ARRAY_PreFFLOW
      use module_rad,	  only:  radmodel,radflag  !Jiang
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: mcell,mvrtx,mface,mssfbc,medge
      integer,intent(in)   :: ncell,ncomp,nrans,ncomp_suf,nphase,npotn
      integer,intent(inout):: nvrtx,nvrtxnw
      integer,intent(out)  :: ierror
!
! --- [local entities]
!--------------------------------------
! --- CV's Physical Arraies:
!--------------------------------------
!
      integer :: ierr=0
      integer :: nface,ncelb,nedge,nssfbc,NCVIN,NCV,NCVFAC,NALLCV,NMAT
      integer :: nssfbc_old
      integer :: IEMAX,MXMAT
      integer :: NBCINL,NBCWAL,NBCCYC,NBCOUT,NBCSYM,NBCTCI,NBCSLD,
     &           NBCINT
      integer :: IVBCTOT
      integer :: IBCCYC,IBCTCI,IBCSLD,IBCINT,IC,IS,I,IIS
      integer :: rewall=0
!
!
!
      IBCINT=0
      IBCCYC=0
      IBCTCI=0
      IBCSLD=0
      
!--------------------------------------
! --- Allocating gloable arries
!--------------------------------------
! --- Cell connectivity
!--------------------------------------
      ALLOCATE (lvcell(8,mcell),stat=ierr)
      if(ierr/=0) stop 'stop at allocating lvcell'
      ALLOCATE (lacell(  mcell),stat=ierr)
      if(ierr/=0) stop 'stop at allocating lacell'
      ALLOCATE (lfcell(7,mcell),stat=ierr)
      if(ierr/=0) stop 'stop at allocating lfcell'
      ALLOCATE (lvface(4,mface),stat=ierr)
      if(ierr/=0) stop 'stop at allocating lvface'
      ALLOCATE (lcface(2,mface),stat=ierr)
      if(ierr/=0) stop 'stop at allocating lcface'
      ALLOCATE (lbface(2,mface),stat=ierr)
      if(ierr/=0) stop 'stop at allocating lbface'
!--------------------------------------
! --- CV connectivity
!--------------------------------------
      ALLOCATE (LEFACE(5,mface),stat=ierr)
      if(ierr/=0) stop 'stop at allocating LEFACE'
      ALLOCATE (LCVFAC(4,mface),stat=ierr)
      if(ierr/=0) stop 'stop at allocating LCVFAC'
      ALLOCATE (LVEDGE(2,medge),stat=ierr)
      if(ierr/=0) stop 'stop at allocating LVEDGE'
      ALLOCATE (LVRTCV(  mvrtx),stat=ierr)
      if(ierr/=0) stop 'stop at allocating LVRTCV'
      ALLOCATE (LBCSSF( mssfbc),stat=ierr)
      if(ierr/=0) stop 'stop at allocating LBCSSF'
      ALLOCATE (LCYCSF( mssfbc),stat=ierr)
      if(ierr/=0) stop 'stop at allocating LCYCSF'
      ALLOCATE (locmsh( mssfbc),stat=ierr)
      if(ierr/=0) stop 'stop at allocating locmsh'
      ALLOCATE (listpr(0:3,mvrtx),stat=ierr)
      if(ierr/=0) stop 'stop at allocating listpr'
      ALLOCATE (listbc(4,mssfbc),stat=ierr)
      if(ierr/=0) stop 'stop at allocating listbc'
!-------------------
! --- Cell Matrix
!-------------------
      ALLOCATE  (volume(  mcell),stat=ierr)
      if(ierr/=0) stop 'stop at allocating volume'
      ALLOCATE  (gcell (3,mcell),stat=ierr)
      if(ierr/=0) stop 'stop at allocating stop'
      ALLOCATE  (area  (4,mface),stat=ierr)
      if(ierr/=0) stop 'stop at allocating area'
      ALLOCATE  (gface (3,mface),stat=ierr)
      if(ierr/=0) stop 'stop at allocating gface'
      ALLOCATE  (cord  (3,mvrtx),stat=ierr)
      if(ierr/=0) stop 'stop at allocating cord'
!-------------------
! --- CV Matrix
!-------------------
      ALLOCATE  (SFAREA(4,medge),stat=ierr)
      if(ierr/=0) stop 'stop at allocating SFAREA'
      ALLOCATE  (SFCENT(3,medge),stat=ierr)
      if(ierr/=0) stop 'stop at allocating SFCENT'
      ALLOCATE  (CVVOLM(mvrtx),stat=ierr)
      if(ierr/=0) stop 'stop at allocating CVVOLM'
      ALLOCATE  (CVCENT(3,mvrtx),stat=ierr)
      if(ierr/=0) stop 'stop at allocating CVCENT'
!-------------------
! --- CV-Face's 
!-------------------
!-------------------
! --- CV Material
!-------------------
      ALLOCATE  (LMAT(mvrtx),stat=ierr)
      if(ierr/=0) stop 'stop at allocating LMAT'
!
      write(ifll,*)
      write(ifll,*)
      write(ifll,5000)
      write(ifll,5010)
      write(ifll,5000)
!--------------------------------------
! --- ZERO Initialization
!--------------------------------------
      write(ifll,'(1x,a)') '###    INITIALIZING, PLEASE WAITING...'
!
      lvcell=0
      lacell=0
      lfcell=0
      lvface=0
      lcface=0
      lbface=0
!
      LEFACE=0
      LCVFAC=0
      LVEDGE=0  
      LVRTCV=0
      LBCSSF=0
      LCYCSF=0
      locmsh=0
      listpr=0
      listbc=0
!
      volume=0.d0
      gcell =0.d0
      area  =0.d0
      gface =0.d0
      cord  =0.d0
!
      SFAREA=0.d0
      SFCENT=0.d0
      CVVOLM=0.d0
      CVCENT=0.d0
!
      LMAT=0
!
      IVBCTOT=0
      IEMAX=0
      NBCINL=0
      NBCWAL=0
      NBCCYC=0
      NBCOUT=0
      NBCSYM=0
      nface=0
      ncelb=0
      nedge=0
      nssfbc=0
      nssfbc_old=0
      NCV=0
      NCVFAC=0
      NALLCV=0
      NMAT=0
      ierror=0
!------------------------------------
! --- < I. Pre-process >- 
!------------------------------------
! --- < 1. Set up topological data >-
!-----------------------------------------------
! --- < 1.1 input data :Read Serial grid data >-
!-----------------------------------------------
!
      call read_grid_data(mcell,mface,mvrtx,ncell,nvrtx,cord,lacell,
     & lvcell,lfcell,lcface,lvface,ierror)
      if(ierror.ne.0) goto 9999
!
!------------------------------------------------------------------------
! --- Creat new vertex between diff-material (interface and sliding BC)
! --- define: nvrtxnw and (set : nvrtx=nvrtxnw)
!------------------------------------------------------------------------
!
      if(gdformat/='NST'.and.gdformat/='GF') then
        if(icon_cv==vertex_cen) then
          write(ifll,'(a)') '###    list_nwvrtx... '
          call list_nwvrtx
     &    (mcell,mface,mvrtx,nvrtx,ncell,nface,ncelb,nvrtxnw,
     &    lvcell,lacell,cord,CVCENT,
     &    ierror)
          if(ierror.ne.0) goto 9999
        endif
      endif
!------------------
! --- Scaling geom
!------------------
      cord(:,:)=gdScale*cord(:,:)
!--------------------
! --- scale factor
!-------------------
      call list_gridfact(mvrtx,mcell,ncell,nvrtx,cord,lacell,lvcell)
!------------------
! --- 
!------------------
      write(ifll,'(a)') '###    storeGridData... '
      call storeGridData(nvrtx,mvrtx,mcell,ncell,
     &   cord,lacell,ierror)
      if(ierror.ne.0) goto 9999
!
!----------------------------------------
! --- < 1.2 make up face, edge lists >---
!----------------------------------------
! --- Geometry face listing 
!----------------------------------------
      write(ifll,'(a)') '###    LIST FACE ... '
!
      call list_face
     & (mcell,mface,mvrtx,nvrtx,ncell,nface,ncelb,NCV,
     &  lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,lbface,
     &  cord,ierror)
      if(ierror.ne.0) goto 9999
!
!-------------------------------------------------------
! --- Periodic BC face listing ("Periodic","Sliding") --
!-------------------------------------------------------
!
      write(ifll,'(a)') '###    LIST FACEPR ...'
      call list_facepr
     & (mface,nface,mvrtx,nvrtx,ncell,mcell,mssfbc,
     &  IBCCYC,IBCSLD,IBCINT,
     &  lvface,lcface,lbface,lacell,listpr,cord,ierror)
      if( ierror.ne.0 ) goto 9999
!------------------------------------
! --- Other Physical BC face listing
!------------------------------------
!      
      write(ifll,'(a)') '###    LIST FACEBC ...'
      call list_facebc
     & (mface,nface,mvrtx,nvrtx,ncell,mcell,
     &  lvface,lbface,lcface,lacell)
!--------------------------------------------------------
! --- Different material inter face listing [interface] 
!--------------------------------------------------------
      write(ifll,'(a)') '###    LIST FACEIN ...'
      call list_facein
     & (mcell,mface,medge,mvrtx,ncell,nface,ncelb,nvrtx,NCV,
     &  lacell,lfcell,lvface,lcface,lbface,lvcell,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierror)
      if(ierror.ne.0) goto 9999
!--------------------------------------------------------
! --- Buffle BC face listing
!--------------------------------------------------------
!      write(ifll,'(a)') '###    LIST FACEBF ...'
!      call list_facebf
!     & (mface,nface,mvrtx,nvrtx,ncell,mcell,ncelb,
!     &  lvcell,lfcell,lvface,lbface,lcface,lacell,cord,
!     &  NCV,LVRTCV,LCVFAC,ierror)
      if( ierror.ne.0 ) goto 9999
!--------------------------------------------------------
! --- Inner BC face listing
!--------------------------------------------------------
!      write(ifll,'(a)') '###    LIST FACENN ...'
!      call list_facenn
!     & (mcell,mface,medge,mvrtx,ncell,nface,ncelb,nvrtx,NCV,
!     &  lfcell,lvface,lcface,lbface,
!     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierror)
!      if(ierror.ne.0) goto 9999
!-----------------------------------------------
! --- Geometry edge listing, make list 'nedge'
!-----------------------------------------------
      write(ifll,'(a)') '###    LIST EDGE ...'
      call list_edge
     & (mcell,mface,medge,mvrtx,ncell,nvrtx,nface,nedge,ncelb,NCV,
     &  lvcell,lfcell,lvface,lcface,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierror)
      if(ierror.ne.0) goto 9999
!------------------
! --- Check face 
!------------------
      write(ifll,'(a)') '###    LIST FACECHK ...'
      call list_facechk
     & (mface,mcell,nface,ncell,mvrtx,nvrtx,
     &  lvface,lcface,lbface,lacell,ierror)
      if(ierror.ne.0) goto 9999
!------------------
! --- Check BC 
!------------------
      write(ifll,'(a)') '###    BC CHECK ...'
      call bc_check
     & (mface,mcell,nface,mvrtx,nvrtx,
     &  lbface,lcface,lvface,lacell,ierror)
      if(ierror.ne.0) goto 9999
!------------------------------------
! --- < 2. Reset attribute data >----
!------------------------------------
      write(ifll,'(a)') '###    LIST SOLID ...'
      call list_solid(mcell,ncell,ncelb,lacell,ierror)
      if(ierror.ne.0) goto 9999
!
      write(ifll,'(a)') '###    LIST FLUID ...'
      call list_fluid
     & (mcell,mface,nface,ncell,ncelb,mvrtx,nvrtx,
     &  lbface,lcface,lacell,ierror)
      if(ierror.ne.0) goto 9999
!------------------------------
! --- < 3. Calculate metrics >-
!---------------------------------------------------------
! --- Also make SSF list and Counter it into egde-list
!---------------------------------------------------------
      write(ifll,'(a)') '###    MAKE METRIC ......'
      call metric_admin
     & (mcell,mvrtx,mface,medge,
     &  ncell,nvrtx,nface,nedge,nssfbc,nssfbc_old,NCV,
     &  mssfbc,NCVFAC,NALLCV,IEMAX,NMAT,ncelb,
     &  lcface,lvface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,LBCSSF,LCYCSF,listbc,listpr,LMAT,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,lvcell,lfcell,
     &  time,-1.d0,ierror,0)
      if(ierror.ne.0) goto 9999
!-------------------------------------------------------------------
! --- List "Touch-inlet" and "Interface" for HPC : "lbface(:,:)"
! --- ("Interface" for HPC NOT finished)
!-------------------------------------------------------------------
      if(NPETOT.gt.1) then
        write(ifll,'(a)') '###    LIST PAIR ......'
        call list_pair
     &       (ncell,nface,mcell,mface,IBCTCI,lbface,gface)
      endif
!---------------------------------------------------------------------
! --- Match "touch-inlet","interface" and "Sliding" BC: "LCYCSF(:,:)"
!---------------------------------------------------------------------
!      write(ifll,'(a)') '###    LIST TOUCH_INLET ...'
!      write(ifll,'(a)') '###    LIST FLUID/SOLIDINTERFACE ...'
!      call list_touch_inlet(mssfbc,NSSFBC,medge,nedge,mvrtx,
!     &                      LVEDGE,LMAT,LCYCSF,LBCSSF,SFCENT)
!-------------------------
! --- Every BC dimension
!-------------------------
      call list_BCDIM
     & (mssfbc,nssfbc,NBCINL,NBCWAL,NBCOUT,NBCCYC,NBCSYM,NBCTCI,NBCSLD,
     &  NBCINT,LBCSSF)
!-------------------
!-< 4. Input data >-
!--------------------------------
! --- < 4.1 initial condition >--
!--------------------------------
      write(ifll,'(a)') '###    READ_INITIAL...'
      call read_initial
     & (NCV,NALLCV,NCVFAC,ncomp,nrans,
     & iters,ierror)
      if( ierror.ne.0 ) goto 9999
!
      write(ifll,'(a)') '###    READ_SOURCE...'
      call read_source(mcell,ncell,ierror)
      if( ierror.ne.0 ) goto 9999
!-------------------
! --- wall distance
!-------------------
      if(NPETOT==1) then
        write(ifll,'(a)') '###    list_wall_distance...'
        call list_wall_distance
     &      (0,rewall,NCV,NBCWAL,nssfbc,NEDGE,
     &       mvrtx,mssfbc,medge,
     &       LMAT,LBCSSF,SFCENT,CVCENT,SFAREA)
      endif
!-------------------------------------
! --- NASTRAN BC grid connectivity
!-------------------------------------
      if(gdformat=='NST') then
        write(ifll,'(a)') '###    list_BC_wall_connectivity...'
        call list_BC_connectivity(nface,mface,lvface,lbface)
      endif
!
!---------------------------------------------------------------
! --- Calculate the Exchange-Area for radiation heat transfer
! --- Jiang 
!---------------------------------------------------------------
      IF(NPETOT==1.AND.radflag==1.and.radmodel(1:3)/='FVM') THEN
	 write(ifll,'(a)') '###   list_radiation_data...'
	 CALL advance_rad
     &   (mvrtx,mcell,mface,nvrtx,ncell,nface,cord,lacell,
     &	 lvcell,lbface,lvface)
		
	END IF
!------------------------------------------------------------------
! --- List material and OUTPUT GEOM FILE for Serial (wall distance)
!------------------------------------------------------------------
!
      NCVIN=NCV  !  if(icon_cv==cell_cen) NCVIN=ncell
      write(ifll,'(a)') '###    list_output_geom...'
      call list_output_geom
     & (0,rewall,mcell,mvrtx,mface,medge,mssfbc,NBCWAL,
     &  NCV,nvrtx,ncell,NMAT,NCVFAC,NALLCV,nedge,nssfbc,nssfbc_old,
     &  NCVIN,
     &  LMAT,LVEDGE,LVRTCV,LBCSSF,LCYCSF,locmsh,listbc,listpr,
     &  SFAREA,SFCENT,CVVOLM,CVCENT)
!
      if(ical_mvmsh>0)
     &   then
        write(ifll,'(a)') '###    list_output_ini_movegrid...'
        call list_output_ini_movegrid
     & (mvrtx,mcell,mface,mssfbc,nvrtx,ncell,nface,nssfbc,nssfbc_old,
     &  cord,
     &  lacell,lvcell,lvface,lbface,
     &  lcface,lfcell,LEFACE,listbc)
      endif
!
!--------------------------------------
! --- OUTPUT GRID FILE (GF)
!--------------------------------------
!
      write(ifll,'(a)') '###    list_output_grid...'
       call list_output_grid
     & (mvrtx,mcell,mface,mssfbc,nvrtx,ncell,nface,nssfbc,nssfbc_old,
     &  cord,
     &  lacell,lvcell,lvface,lbface,
     &  lcface,lfcell)
!
      write(ifll,'(a)') '###    list_output_grid FFR-GF...'
      if(ffrgridOutFlag) then
        if(ffgFormatKey=='A') then
          call list_output_grid_A
     &       (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell)
        else if(ffgFormatKey=='U') then
          call list_output_grid_U
     &       (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell)
        else
          write(ifll,'(a)') 'Not yet support '
          ierror=1;return
        end if
      end if
!---------------------------------------------------------
! --- Serial cpu : dimen.h_serial (mxcell,mxvrtx ...)
!---------------------------------------------------------
      write(ifll,5000)
      call wrtdmn
     & (nvrtx,ncell,ncelb,nface,nedge,NCVIN,
     & NCV,NALLCV,NCVFAC,nssfbc,nssfbc_old,
     % ncomp,ncomp_suf,nphase,npotn,
     & nrans,NMAT,IEMAX,
     & 0,NBCWAL,NBCSLD,NBCINT,NPETOT,0)
      MXMAT=NMAT
!--------------------------------------
! --- HPC Domain Partitioning
!--------------------------------------
      IF(NPETOT.gt.1) then
        write(ifll,5000)
        write(ifll,5020)
!
        if(icon_cv==vertex_cen) then
          call part
     &  (mcell,mvrtx,mface,medge,mssfbc,
     &   ncell,ncelb,nvrtx,nface,nedge,nssfbc,ncomp,ncomp_suf,npotn,
     &   nphase,nrans,NMAT,MXMAT,NBCWAL,
     &   NCV,NCVFAC,NALLCV,IBCCYC,IBCTCI,IBCSLD,IBCINT,IEMAX,
     &   lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,lbface,listpr,
     &   LVEDGE,LEFACE,lacell,
     &   LBCSSF,LCYCSF,LMAT,locmsh,listbc,
     &   cord,area,volume,gface,gcell,
     &   SFAREA,SFCENT,CVVOLM,CVCENT,
     &   time,ierror)
         if(ierror.ne.0) goto 9999
        else
          write(ifll,*) 'ERR: HPC ONLY support vertex-center '
          stop 
        endif
      endif
!
!----------------------
! --- List Sliding BC
!----------------------
!
!      write(ifll,'(a)') '###    Serial: list_output_sld...'
!      call list_output_sld
!     & (1,nssfbc,nedge,mssfbc,medge,LBCSSF,LCYCSF,SFCENT,SFAREA)
!
!-------------------
! --- END
!-------------------
      return
!
 9999 continue
!
 5000 format('|',108('='),'|')
 5010 format('|',46(' '),' LIST SERIAL CPU MESH ',40(' '),'|')
 5020 format('|',46(' '),' LIST HPC    CPU MESH ',40(' '),'|')
!
      write(ifle,*) '(pre_main)'
      ierror=1
      end subroutine pre_main
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine wrtdmn
     & (nvrtx,ncell,ncelb,nface,nedge,NCVIN,
     &  NCV,NALLCV,NCVFAC,nssfbc,nssfbc_old,
     &  ncomp,ncomp_suf,nphase,npotn,
     &  nrans,NMAT,IEMAX,
     &  IVBCTOT,NBCWAL,NBCSLD,NBCINT,NPETOT,icpu)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     NPETOT=1: serial computation
!     NPETOT>1: HPC computation
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only          : ifli,ifll,ifle
!      use module_partitioner,only : hpcdmn
      use module_Euler2ph,only    : KE2P,ieul2ph,NPHS
      use module_scalar,only      : KRANS
      use module_boundary,only    : KSUF
      use module_partitioner
      use module_scalar,only      : ical_POTN,Kpotn
      use module_model,only       : KMHD
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: nvrtx,ncell,ncelb,nface,nedge
      integer,intent(in)  :: icpu,NPETOT
      integer,intent(in)  :: NCVIN,NCV,NALLCV,NCVFAC
      integer,intent(in)  :: nssfbc,nssfbc_old,
     &                       ncomp,nrans,NMAT,ncomp_suf,nphase,npotn
      integer,intent(in)  :: IEMAX,IVBCTOT,NBCWAL,NBCSLD,NBCINT
!
!
! --- [local entities]
!
      character(80)        :: dimh
      integer :: ios
      integer :: LH=10,LENTH
      integer :: KSLD=0
!
! --- 
!
      IF(NPETOT.gt.1.and.icpu.gt.0) THEN
! --- HPC dimension in sub-directory : [dimen.parm]
        dimh=hpcdmn(icpu)
        dimh=adjustL(dimh)
        LENTH=len_trim(dimh)
        open(LH,file=dimh,form='formatted',status='UNKNOWN',iostat=ios)
        if(ios.ne.0) then
          write(ifll,*) 'ERR: file: dimen.parm, open error'
          write(ifll,*) 'iostat =',ios
          stop
        endif
        write(LH,100) 
        write(LH,1000) max(nvrtx  ,1),  nvrtx_hpc(icpu)
        write(LH,1001) max(ncell  ,1),  ncell_hpc(icpu)
        write(LH,1002) max(ncelb  ,1),  ncelb_hpc(icpu)
        write(LH,1003) max(nface  ,1),  nface_hpc(icpu)
        write(LH,1004) max(nedge  ,1),  nedge_hpc(icpu)
        write(LH,1018) max(NCVIN  ,1),  NODtotL(icpu)
        write(LH,1005) max(NCV    ,1),    NCV_hpc(icpu)
        write(LH,1006) max(NALLCV ,1), NALLCV_hpc(icpu)
        write(LH,1007) max(NCVFAC ,1), NCVFAC_hpc(icpu)
        write(LH,1008) max(nssfbc ,1), nssfbc_hpc(icpu)
        write(LH,1009) max(ncomp  ,1),  ncomp_hpc(icpu)
        write(LH,1013) max(ncomp_suf,1),ncomp_suf_hpc(icpu)
        write(LH,1017) max(npotn,1),npotn_hpc(icpu)
        write(LH,1016) max(nphase ,1),  nphase_hpc(icpu)
        write(LH,1010) max(nrans  ,1),  nrans_hpc(icpu)
        write(LH,1011) max(NMAT   ,1),   NMAT_hpc(icpu)
        write(LH,2170) max(nssfbc_OLD ,1), nssfbc_hpc(icpu)
        if(NBCSLD>0) then
          write(LH,1012) max(IEMAX,1)+3,  IEMAX_hpc(icpu)+3
        else
          write(LH,1012) max(IEMAX  ,1),  IEMAX_hpc(icpu)
        endif
        write(LH,1019) max(IVBCTOT,1),IVBCTOT_hpc(icpu)
        write(LH,1014) max(NBCWAL,1) , NBCWAL_hpc(icpu)
        write(LH,1015) max(NBCSLD,1) , NBCSLD_hpc(icpu)
      ELSEIF(icpu.eq.0) THEN
! --- Serial cpu dimension file [dimen.h_serial] in work-directory
        dimh='dimen.h_serial'
        dimh=adjustL(dimh)
        LENTH=len_trim(dimh)
        open(LH,file=dimh,form='formatted',status='UNKNOWN',iostat=ios)
        if(ios.ne.0) then
          write(ifll,'(1x,a)') 'ERR: file: dimen.h_serial, open error'
          write(ifll,'(1x,a,I4)') 'iostat =',ios
          stop
        endif
        write(LH,100) 
        write(LH,1000) max(nvrtx  ,1),  nvrtx
        write(LH,1001) max(ncell  ,1),  ncell
        write(LH,1002) max(ncelb  ,1),  ncelb
        write(LH,1003) max(nface  ,1),  nface
        write(LH,1004) max(nedge  ,1),  nedge
        write(LH,1018) max(NCVIN  ,1),  NCVIN
        write(LH,1005) max(NCV    ,1),    NCV
        write(LH,1006) max(NALLCV ,1), NALLCV
        write(LH,1007) max(NCVFAC ,1), NCVFAC
        write(LH,1008) max(nssfbc ,1), nssfbc
        write(LH,1009) max(ncomp  ,1),  ncomp
        write(LH,1013) max(ncomp_suf,1),ncomp_suf
        write(LH,1017) max(npotn,1),npotn
        write(LH,1016) max(nphase ,1),  nphase
        write(LH,1010) max(nrans  ,1),  nrans
        write(LH,1011) max(NMAT   ,1),   NMAT
        write(LH,2170) max(nssfbc_old ,1),nssfbc_old
        if(NBCSLD>0) then
          write(LH,1012) max(IEMAX,1)+3,  IEMAX+3
        else
          write(LH,1012) max(IEMAX  ,1),  IEMAX
        endif
        write(LH,1019) max(IVBCTOT,1),IVBCTOT
        write(LH,1014) max(NBCWAL,1) , NBCWAL
        write(LH,1015) max(NBCSLD,1) , NBCSLD
      ENDIF
!     
      if(NBCSLD>0) then
        KSLD=1
      endif
      write(LH,2100) 
      write(LH,2107) KSLD*nssfbc+1-KSLD
!
      write(LH,2100) 
      write(LH,2108) KSUF
      write(LH,2109) KSUF*nssfbc+1-KSUF
!
      write(LH,2100) 
      write(LH,2101) KRANS
      write(LH,2102) KRANS*NCVIN+1-KRANS 
      write(LH,2103) KRANS*NCV+1-KRANS   
      write(LH,2104) KRANS*NALLCV+1-KRANS
      write(LH,2105) KRANS*NCVFAC+1-KRANS
      write(LH,2106) KRANS*nssfbc+1-KRANS
!
      write(LH,1999) 
      write(LH,2001) KE2P
      write(LH,2007) NPHS
      write(LH,2002) KE2P*NCVIN+1-KE2P 
      write(LH,2003) KE2P*NCV+1-KE2P 
      write(LH,2004) KE2P*NALLCV+1-KE2P
      write(LH,2005) KE2P*NCVFAC+1-KE2P
      write(LH,2006) KE2P*nssfbc+1-KE2P
!
      write(LH,2100) 
      write(LH,2201) Kpotn
      write(LH,2202) Kpotn*NCVIN+1-Kpotn 
      write(LH,2203) Kpotn*NCV+1-Kpotn   
      write(LH,2204) Kpotn*NALLCV+1-Kpotn
      write(LH,2205) Kpotn*NCVFAC+1-Kpotn
      write(LH,2206) Kpotn*nssfbc+1-Kpotn
      write(LH,200) 
      close(LH)
!

!
 100  FORMAT('&DIMENSION_SIZE')
 200  FORMAT('/')
 1000 FORMAT(6X,' MXVRTX  =',I10,', NNVRTX =',I10)
 1001 FORMAT(6X,' MXCELL  =',I10,', NNCELL =',I10)
 1002 FORMAT(6X,' MXCELB  =',I10,', NNCELB =',I10)
 1003 FORMAT(6X,' MXFACE  =',I10,', NNFACE =',I10)
 1004 FORMAT(6X,' MXEDGE  =',I10,', NNEDGE =',I10)
 1005 FORMAT(6X,' MXCV    =',I10,', NNCV   =',I10)
 1006 FORMAT(6X,' MXALLCV =',I10,', NNALLCV=',I10)
 1007 FORMAT(6X,' MXCVFAC =',I10,', NNCVFAC=',I10)
 1008 FORMAT(6X,' MXSSFBC =',I10,', NNSSFBC=',I10)
 1009 FORMAT(6X,' MXCOMP  =',I10,', NNCOMP =',I10)
 1010 FORMAT(6X,' MXRANS  =',I10,', NNRANS =',I10)
 1011 FORMAT(6X,' MXMAT   =',I10,', NNMAT  =',I10)
 1012 FORMAT(6X,' MAXIE   =',I10,', NNIE   =',I10)
 1013 FORMAT(6X,' MXCOMP_SUF  =',I10,', NNCOMP_SUF =',I10)
 1017 FORMAT(6X,' MXpotn  =',I10,', NNpotn =',I10)
 1016 FORMAT(6X,' MXphase =',I10,', NNphase=',I10)
 1014 FORMAT(6X,' MXBCWAL =',I10,', NNBCWAL=',I10)

 1015 FORMAT(6X,' MXBCSLD =',I10,', NNBCSLD=',I10)
! 1017 FORMAT(6X,' MXBCSYM =',I10,', NNBCSYM=',I10)
 1018 FORMAT(6X,' MXCVIN  =',I10,', NNCVIN =',I10)
 1019 FORMAT(6X,' MXBCTOT =',I10,', NNBCTOT=',I10)
 1020 FORMAT(6X,' MXBCTCI =',I10,', NNBCTCI=',I10)
!
 1999 FORMAT(6X)
 2001 FORMAT(6X,' KE2P    =',I10)
 2007 FORMAT(6X,' NPHS    =',I10)
 2002 FORMAT(6X,' MXCVIN2 =',I10)
 2003 FORMAT(6X,' MXCV2   =',I10)
 2004 FORMAT(6X,' MXALLCV2=',I10)
 2005 FORMAT(6X,' MXCVFAC2=',I10)
 2006 FORMAT(6X,' MXSSFBC2=',I10)
!
 2100 FORMAT(6X)
 2101 FORMAT(6X,' KRANS   =',I10)
 2102 FORMAT(6X,' MXCVINR =',I10)
 2103 FORMAT(6X,' MXCVR   =',I10)
 2104 FORMAT(6X,' MXALLCVR=',I10)
 2105 FORMAT(6X,' MXCVFACR=',I10)
 2106 FORMAT(6X,' MXSSFBCR=',I10)
!
 2107 FORMAT(6X,' MXSSFBC_SLD=',I10)
!
 2201 FORMAT(6X,' KPOTN   =',I10)
 2202 FORMAT(6X,' MXCVINP =',I10)
 2203 FORMAT(6X,' MXCVP   =',I10)
 2204 FORMAT(6X,' MXALLCVP=',I10)
 2205 FORMAT(6X,' MXCVFACP=',I10)
 2206 FORMAT(6X,' MXSSFBCP=',I10)
!
!
 2108 FORMAT(6X,' KSUF       =',I10)
 2109 FORMAT(6X,' MXSSFBC_SUF=',I10)
!
 2170 FORMAT(6X,' MXSSFBC_OLD =',I10,', NNSSFBC_OLD=',I10)
!
 2000 FORMAT(19I10)
!
      write(ifll,5010) dimh(1:LENTH)
 5010 format(5X,'###',4X,a,' file has been created successfully')
!
      return
      end subroutine wrtdmn

