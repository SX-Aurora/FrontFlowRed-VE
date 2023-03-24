!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_FFRGFwords
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character*80,parameter :: asciiv2_FFGF="#A_GF_V2"
      character*80,parameter :: unformv2_FFGF="#U_GF_V2"
      character*80,parameter :: binaryv2_FFGF="#B_GF_V2"
      character*80,parameter :: newSet_FFGF="#NEW_SET"
      character*80,parameter :: fileend_FFGF="#ENDFILE"
      character*80,parameter :: customData_FFGF="#CUSTOM_DATA"
      integer,parameter :: unkownDatasize_FFGF=0
!     For grids.
      character*80,parameter :: ffrGrid_FFGF="FFR_GRID_FIELD"
      character*80,parameter :: gridHead_FFGF="#FFR_GRID_HEADER"
      character*80,parameter :: gridHeade_FFGF="#FFR_GRID_HEADER_END"
      character*80,parameter :: elemType_FFGF="#ELEM_TYPE_LIST"
      character*80,parameter :: elemTypee_FFGF="#ELEM_TYPE_LIST_END"
      character*80,parameter :: bndryType_FFGF="#BNDRY_TYPE_LIST"
      character*80,parameter :: bndryTypee_FFGF="#BNDRY_TYPE_LIST_END"
      character*80,parameter :: vrtxCord_FFGF="#VRTX_CORD"
      character*80,parameter :: vrtxCorde_FFGF="#VRTX_CORD_END"
      character*80,parameter :: bndryFace_FFGF="#BNDRY_FACE"
      character*80,parameter :: bndryFacee_FFGF="#BNDRY_FACE_END"
      character*80,parameter :: elemVrtx_FFGF="#ELEM_VRTX"
      character*80,parameter :: elemVrtxe_FFGF="#ELEM_VRTX_END"
      character*80,parameter :: undefMatName="MAT"
      integer,parameter :: undefBndryType=0
      integer,parameter :: defaultCellTypeID=1
!     For fulid force.
      character*80,parameter :: ffrForce_FFGF="FFR_FORCE_FIELD"
      character*80,parameter :: forceHead_FFGF="#FFR_FORCE_HEADER"
      character*80,parameter :: forceHeade_FFGF="#FFR_FORCE_HEADER_END"
      character*80,parameter :: forceWall_FFGF="#FFR_FORCE_WALLS"
      character*80,parameter :: forceWalle_FFGF="#FFR_FORCE_WALLS_END"
      character*80,parameter :: forceModel_FFGF="#FFR_FORCE_MODEL"
      character*80,parameter :: forceModele_FFGF="#FFR_FORCE_MODEL_END"
      character*80,parameter :: forceObserv_FFGF="#FFR_FORCE_OBSERVER"
      character*80,parameter :: forceObserve_FFGF=
     &                                    "#FFR_FORCE_OBSERVER_END"
      character*80,parameter :: forceData_FFGF="#FFR_FORCE_DATA"
      character*80,parameter :: forceDatae_FFGF="#FFR_FORCE_DATA_END"


      end module module_FFRGFwords
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- -------------------------------------------------------------------
      module FFRreaddata
! --- -------------------------------------------------------------------
      character(LEN=9),save :: cntlnam='fflow.ctl'

      integer :: ianim,ianim_uvw,ianim_p,ianim_r,ianim_t
      integer,allocatable :: ianim_GAS_WDOT(:),
     &                       ianim_SDOT_suf(:),
     &                       ianim_molefr_suf(:),
     &                       ianim_depoeth_spd(:)
!
      integer :: iuvw_ave_rms_rex
      integer :: ip_avex,it_avex,imu_avex
      integer :: ip_rmsx,it_rmsx,imu_rmsx
      integer :: iuvwt_rex,imaxminx
!      integer :: iflag
      integer,allocatable :: icomp_avex(:),irans_avex(:)
      integer,allocatable :: icomp_rmsx(:),irans_rmsx(:)
      integer,allocatable :: iuvw_rans_rex(:)
      integer,allocatable :: iuvwc_rex(:)
! wallvar
      integer,allocatable :: iwallvar(:)
      end module FFRreaddata
! --- -------------------------------------------------------------------
      module SCRYUheaderTag
! --- -------------------------------------------------------------------
      character*8, parameter :: HEADER = 'CRDL-FLD'
      character*8, parameter :: SOLV   = 'SCRYUTET'
      character*80,parameter :: COMENT = ' '
      integer,parameter      :: IVSN   = 2
      integer,parameter      :: ITYP   = 1
      integer :: L2D3D=3,LFORT=0
      end module SCRYUheaderTag
! --- -------------------------------------------------------------------
      module SCRYUdata
! --- -------------------------------------------------------------------
      integer :: IDAT,NCYC,NDATA,LVCT
      integer :: IRECL,LNX,LCORD,NGFAX,NBNN,NTTS
      integer :: LKAI,NTTE,NRTY
      real(4) :: TIME,AMNUM
      character*32 :: TITLE,OverlapStart_n,OverlapEnd
      character*32 :: LNAM,IRTY
      character*12 :: LRGN,MRGN
      integer,allocatable :: IPTYP(:),NDFA(:),IPMAT(:),sumnod(:)
      end module SCRYUdata
! --- -------------------------------------------------------------------
      module SCRYUpre
! --- -------------------------------------------------------------------
      character*12 :: LREGN
      integer :: NNODS,NELEM,NREGN,NREGI,IDVER,NGR
     &          ,NDE,NE
      integer, allocatable :: LGR(:),MAT(:),NDNO(:),IE(:),IFA(:),
     &                        IEtot(:,:),IFAtot(:,:),
     &                        fnode(:),FTYP(:,:)
      real(4), allocatable :: xcrd(:),ycrd(:),zcrd(:)
      character*1, allocatable :: IETYP(:)
      character*32,allocatable :: NKND(:)
      end module SCRYUpre
! --- -------------------------------------------------------------------
      module FLUENTdata
! --- -------------------------------------------------------------------
      character(10),parameter  :: F_VERSION='fluent6.x'
      integer                  :: F_NREGN   ! same as SCRYUpre
      character*32,allocatable :: F_NKND(:)
      real(8),allocatable      :: F_UVW(:,:)
!
      integer,allocatable :: F_NFBOUN(:)
      end module FLUENTdata
!
! --- -------------------------------------------------------------------
      module AVSheaderTag
! --- -------------------------------------------------------------------
!  AVS Header tags.
      character*20,parameter :: AVS_NAME = '# UCD (AVS/Express)'
      character*20,parameter :: AVS_STEP = '1'
      character*20,parameter :: AVS_CYCL = 'data'
      character*20,parameter :: AVS_STPN = 'step1 (fixed time)'
      end module AVSheaderTag
!
!-----------------------------------------------------------------------
      module FVheaderTag
!-----------------------------------------------------------------------
!  FV-UNS Header tags.
      integer,parameter :: FV_MAGIC =66051
      integer,parameter :: FV_NODES =1001
      integer,parameter :: FV_FACES =1002
      integer,parameter :: FV_ELEMENTS    =1003
      integer,parameter :: FV_VARIABLES   =1004
      integer,parameter :: FV_BNDRY_VARS  =1006
      integer,parameter :: FV_TET_ELEM_ID       =1
      integer,parameter :: FV_HEX_ELEM_ID       =2
      integer,parameter :: FV_PRISM_ELEM_ID     =3
      integer,parameter :: FV_PYRA_ELEM_ID      =4
      integer,parameter :: A_WALL         =7
      integer,parameter :: NOT_A_WALL     =0
      character*80,parameter :: FV_NAME = 'FIELDVIEW'
      end module FVheaderTag
!
!-----------------------------------------------------------------------
      module FVdata
!-----------------------------------------------------------------------
!  Grid parameters
      integer :: numGGrp      ! Number of grid group
      integer,allocatable :: sfRsltFlg(:)! Boundary variables flag
      integer,allocatable :: sfClkFlg(:) ! Boundary surface clockness flag
!  Simulation parameters
      real(4) :: timeFV,fsmachFV,alphaFV,reFV
      integer,allocatable :: hcell(:)
!
      end module FVdata

!-----------------------------------------------------------------------
      module ENSIGHTdata
!-----------------------------------------------------------------------
!  Simulation parameters
      real(4),allocatable :: timeES(:)
!
      end module ENSIGHTdata
!
!-----------------------------------------------------------------------
      module FFRdata
!-----------------------------------------------------------------------
!  Grid parameter
      integer :: nvrtx
      integer :: ncell
      integer :: NBOUND
      CHARACTER*80,allocatable :: SFBOUN(:)
      real(8) :: gdScale=1.d0
      real(8),allocatable :: cord(:,:)
! DEBUG
      real(8),allocatable :: cord_bak(:,:)
! onishi
      integer :: NFCE   ! Total number of faces
      integer :: NBFS
      integer,allocatable :: NFBOUN(:),IBFACE(:),IFFACE(:,:)
      integer,allocatable :: lacell(:),lvcell(:,:)
! onishi
      integer,allocatable :: kmesh(:)
      real(8),allocatable :: totdata(:,:)
      real(8),allocatable :: bnddata(:,:)
!
      real(8),allocatable :: wallDist(:)
! following are used only in FLUENT output
      integer :: nface
      integer,allocatable :: lfcell(:,:)
      integer,allocatable :: lvface(:,:)
      integer,allocatable :: lcface(:,:)
      integer,allocatable :: lbface(:,:)
      integer,allocatable :: LVRTCV(  :)
      integer,allocatable :: LCVFAC(:,:)
      integer :: nbcnd
      integer,parameter :: numbname=2
      integer,allocatable :: boundIDMap(:,:)
      integer,allocatable :: ivbcnd(:)
      integer,allocatable :: lvbcnd(:)

!  Variable paramegers
      integer :: numVars      ! Number of variables
      character*80,allocatable :: nameVars(:)   ! Names of variables
      integer :: numBVars     ! Number of boundary variables
      character*80,allocatable :: nameBVars(:)  ! Names of boundary variables

!  FFR header
      integer :: iterFFR,ncvFFR,ncvfacFFR,ncompFFR,nrnsxFFR,ieul2ph,
     &           ncompallFFR,icon_cvx,icon_cvg
      integer :: npotnxFFR,NFLIDxFFR
      integer :: idrdpFFR,compFFR
      integer :: NMATFFR,ical_sldFFR=0,ical_sufFFR=0,
     &           ical_MHD,ical_mvmsh,ical_WDRFFR=0,N_INJ
      integer :: SDOT_suf,WDOT_gas,depoeth_spd,molefr_suf,
     &           nsdot,ndeps,nmolfr,num_ratout,blk_thick,nthick
      real(8) :: timeFFR
!
      integer,allocatable :: ishaft(:),mat_no(:),nofld(:)
      real(8),allocatable :: rotati(:),end(:,:),begin(:,:),rot_ang(:)
!
      integer :: MPprocNum    ! Number of CPUs
      integer :: outGGrpNum   ! Output grid number
      logical,allocatable :: outGGFlg(:)  ! Output grid flag
      integer,allocatable :: node_id(:)   ! CPU No. of every vrtx.
! wallvar
      integer,allocatable :: LBC_INDEX(:) ! Boundary Faces pointer
      integer,allocatable :: LBC_ICV(:)   ! Boundary Face -> CV pointer
      integer:: radflag

      end module FFRdata
!
!-----------------------------------------------------------------------
      module ProgData
!-----------------------------------------------------------------------
!      integer :: MPprocNum    ! Number of CPUs
!      integer :: outGGrpNum   ! Output grid number
!      logical,allocatable :: outGGFlg(:)  ! Output grid flag
      !   Y+/u+
      logical :: ypAverage    ! Yplus average output flag
      logical :: ypLogDiv
      integer :: divnumYPAV
      integer :: wdisCalcType ! 0: read file 1:calc channel 2:calc Cylinder
      real(4) :: axis1,axis2,Radius,scaleFactor
      character*1 :: axisTag
      character*80 :: wdisYPAVfn,ypAvupfn
      real(4),allocatable :: ypVal(:),upVal(:,:)
      integer,allocatable :: ypAvCnt(:)
      real(4) :: ReTau,uTau,wdMax,wdMin,dyp
      type gridDataMP
        integer :: procID
        integer :: numNodes    ! Number of nodes
        real(4),pointer :: XYZcord(:,:) ! xyz coordinates
        integer :: numCell      ! Number of cell
        integer,pointer :: ndListCel(:,:)      ! List of nodes in cell
        integer,pointer :: mtTypeCel(:)      ! Material type
        integer :: numNodeBReg
        integer,pointer :: nodeBRegPtr(:) ! Number of nodes in boundary region
        integer,pointer :: nodeBReg(:) ! List of nodes in boundary region
        real(8),pointer :: wallDistMP(:)
      end type
!
      end module ProgData
!-----------------------------------------------------------------------
