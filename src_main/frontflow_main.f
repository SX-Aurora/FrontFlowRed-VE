!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 3.1)                           !
!                                                                      !
!                 This code is re-developed by NuFD CO.LTD             !
!     Main Program : FRONTFLOW                                         !
!                          Both Cell-Center and Vertex-Center          !
!                                                                      !
!                          Developer: Huilai   ZHANG                   !
!                          Adviser  :                                  !
!                                     Nobuyuki OSHIMA                  !
!                                     Makoto   TSUBOKURA               !
!                                     Takuji   TOMINAGA                !
!                                     Toshio   KOBAYASHI               !
!                                     2008/06/15                       !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 3.0beta)                       !
!                                                                      !
!     Main Program : FRONTFLOW   (Vector Routine (for vector computer) !
!                                 Lagrangian 2P[Particle Tracking]     !
!                                 LES Flamelet conbustion model        ! 
!                                 Radiation [FVM,MCM & ZONE]           !  
!                                 Dynamical SGS modle                  !
!                                                                      !
!                          Written by Huilai   ZHANG                   !
!                                     Yuyan    Jiang                   !
!                                     Takuji   TOMINAGA                !
!                                     Takeshi  UNEMURA                 !
!                                     Nobuyuki OSHIMA                  !
!                                     Makoto   TSUBOKURA               !
!                                     2007/06/01                       !
!                                                                      !
!     Contact address : The University of Tokyo, RSS21 Project         !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 2.8a)                          !
!                                                                      !
!     Main Program : FRONTFLOW    (Vector Routine, DES,                !
!                                      Serval RANS Models)             !
!                                                                      !
!                          Written by Huilai   ZHANG                   !
!                                     Takuji   TOMINAGA                !
!                                     Keiji    ONISHI                  !
!                                     Eisuke   YAMADA,                 !
!                                     Takafumi SUGINAKA                !
!                                     Yuyan    Jiang                   !
!                                     Takeshi  UNEMURA                 !
!                                     Nobuyuki OSHIMA                  !
!                                     Makoto   TSUBOKURA               !
!                                     2006/06/15                       !
!                                                                      !
!     Contact address : The University of Tokyo, RSS21 Project         !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 2.0)                           !
!                                                                      !
!     Main Program : FRONTFLOW    (Vertex-Center,EDGE-BASE,HPC,        !
!                                  Multi-Domain[thermal solid coupled])!
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

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_dimension
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- array dimension
!
      INTEGER :: NNVRTX ,NVRTX_L  ,MXVRTX, NVRTX
      INTEGER :: NNCELL ,NCELL_L  ,MXCELL, NCELL
      INTEGER :: NNCELB ,NCELB_L  ,MXCELB, NCELB
      INTEGER :: NNFACE ,NFACE_L  ,MXFACE, NFACE
      INTEGER :: NNEDGE ,NEDGE_L  ,MXEDGE, NEDGE
      INTEGER :: NNCVIN ,NCVIN_L  ,MXCVIN, NCVIN
      INTEGER :: NNCV   ,NCV_L    ,MXCV   ,NCV
      INTEGER :: NNALLCV,NALLCV_L ,MXALLCV,NALLCV
      INTEGER :: NNCVFAC,NCVFAC_L ,MXCVFAC,NCVFAC
      INTEGER :: NNSSFBC,NSSFBC_L ,MXSSFBC,NSSFBC
      INTEGER :: NNSSFBC_OLD,NSSFBC_OLD_L,MXSSFBC_OLD,NSSFBC_OLD
      INTEGER :: NNCOMP ,NCOMP_L  ,MXCOMP, NCOMP
      INTEGER :: NNCOMP_SUF ,NCOMP_SUF_L,MXCOMP_SUF, NCOMP_SUF
      INTEGER :: NNRANS ,NRANS_L  ,MXRANS, NRANS
      INTEGER :: NNPOTN ,NPOTN_L  ,MXPOTN, NPOTN
      INTEGER :: NNMAT  ,NMAT_L   ,MXMAT,  NMAT
      INTEGER :: NNIE   ,IEMAX_L  ,MAXIE,  IEMAX
      INTEGER :: NNBCTOT,NBCTOT_L ,MXBCTOT,NBCTOT
      INTEGER :: NNBCWAL,NBCWAL_L ,MXBCWAL,NBCWAL
      INTEGER :: NNBCSLD,NBCSLD_L ,MXBCSLD,NBCSLD
      
      INTEGER :: MXCV_V,MXBND_V,MXCV_VP,MXBND_VP,IVDIM
      INTEGER :: MXCV_D,NFLID
      INTEGER :: MXCV_MHD,MXBND_MHD,MXFACE_MHD
      INTEGER :: MXVRTX_m,MXCELL_m,MXCELB_m,MXFACE_m,MXEDGE_m,MXSSFBC_m,
     &           MXSSFBC_OLD_m
      INTEGER :: MXCV_B
      INTEGER :: MXCV_F,MXCVFAC_F,MXCVFAC_B
!
      INTEGER :: MXSSFBC_VOF
      INTEGER :: MXALLCV_VOF
!
      INTEGER :: NG=10
!
      INTEGER :: MXALLCVC
!
      INTEGER :: MXSSFBC_SLD
!
      INTEGER :: MXSSFBC_SHT
!      INTEGER :: MXSSFBC_OLD,NBCSLD_OLD
!
      INTEGER :: MXCV_WDR,N_inj !MXSSFBC_WDR
!
      INTEGER :: KSUF
      INTEGER :: MXSSFBC_SUF
      INTEGER :: MXCOMPALL,NCOMPALL
      INTEGER :: MXPHASE,NPHASE,NNPHASE,NPHASE_L
!
      INTEGER :: KRANS
      INTEGER :: MXCVINR
      INTEGER :: MXCVR
      INTEGER :: MXALLCVR
      INTEGER :: MXCVFACR
      INTEGER :: MXSSFBCR
!
      INTEGER :: KE2P
      INTEGER :: NPHS
      INTEGER :: MXCVIN2
      INTEGER :: MXCV2
      INTEGER :: MXALLCV2
      INTEGER :: MXCVFAC2
      INTEGER :: MXSSFBC2
!
      INTEGER :: KPOTN
      INTEGER :: MXCVINP
      INTEGER :: MXCVP
      INTEGER :: MXALLCVP
      INTEGER :: MXCVFACP
      INTEGER :: MXSSFBCP
!
      INTEGER :: MXPRT,MEP,MXALLCV_P,MXCVFAC_P,NPRT,MEVENT
!
      INTEGER :: Ngauss_R,NDIV_R
      INTEGER :: MXALLCV_RAD,NDNG,MXALLNG
!
      INTEGER :: HPC_dead,A_dead,BC_dead,P_Reced,MASS_dead
!
      namelist /DIMENSION_SIZE/ 
     &     NNVRTX ,MXVRTX,
     &     NNCELL ,MXCELL,
     &     NNCELB ,MXCELB,
     &     NNFACE ,MXFACE,
     &     NNEDGE ,MXEDGE,
     &     NNCVIN ,MXCVIN,
     &     NNCV   ,MXCV,
     &     NNALLCV,MXALLCV,
     &     NNCVFAC,MXCVFAC,
     &     NNSSFBC,MXSSFBC,
     &     NNSSFBC_OLD,MXSSFBC_OLD,
     &     NNCOMP ,MXCOMP,
     &     NNCOMP_SUF,MXCOMP_SUF,
     &     NNPOTN ,MXPOTN,
     &     NNPHASE,MXPHASE,
     &     NNRANS ,MXRANS,
     &     NNMAT  ,MXMAT,
     &     NNIE   ,MAXIE,
     &     NNBCTOT,MXBCTOT,
     &     NNBCWAL,MXBCWAL,
     &     NNBCSLD,MXBCSLD,
!
     &     MXSSFBC_SLD,
!
     &     KSUF,
     &     MXSSFBC_SUF,
     &     MXPHASE,
!
     &     KRANS,
     &     MXCVINR,
     &     MXCVR,
     &     MXALLCVR,
     &     MXCVFACR,
     &     MXSSFBCR,
!
     &     KE2P,
     &     NPHS,
     &     MXCVIN2,
     &     MXCV2,
     &     MXALLCV2,
     &     MXCVFAC2,
     &     MXSSFBC2,
!
     &     KPOTN,
     &     MXCVINP,
     &     MXCVP,
     &     MXALLCVP,
     &     MXCVFACP,
     &     MXSSFBCP    
!
!/////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine n_size(NPE,Dmns,
     &   ivector,iMHD,iMvmsh,ivof,
     &   iBODY,FLMLT,iCAVIT,iPEFC,iprtcle,iinj,
     &   ifld,MXPRTX,ICALL)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: NPE,ICALL,ivector,iMHD,
     &                      iMvmsh,iBODY,FLMLT,MXPRTX,
     &                      iCAVIT,ifld,iPEFC,iprtcle,iinj,ivof
      character(len=80),intent(inout) :: Dmns
!
! --- [local entities]
!
      integer             :: ifl=11,ios=0
!
      if(ICALL.eq.0) then
        Dmns='dimen.h_serial'
        open(ifl,file=Dmns,form='formatted',status='UNKNOWN',iostat=ios)
        if(ios.ne.0) call FFRABORT(1,'ERR: [dimen.h_serial] open error')
        rewind ifl
        read(ifl,DIMENSION_SIZE,iostat=ios)
        if(ios.ne.0) call FFRABORT(1,'ERR: [dimen.h_serial] read error')
        close(ifl)
        NRANS=NNRANS
        NCOMP=NNCOMP
        NCOMP_SUF=NNCOMP_SUF
        NPHASE=NNPHASE
        NCOMPALL=NNCOMP_SUF+NNCOMP
        MXCOMPALL=NNCOMP_SUF+NNCOMP
        NPOTN=NNPOTN
      elseif(ICALL.eq.1) then
        Dmns='dimen.h_serial'
        open(ifl,file=Dmns,form='formatted',status='UNKNOWN',iostat=ios)
        if(ios.ne.0) call FFRABORT(1,'ERR: [dimen.h_serial] open error')
        rewind ifl
        read(ifl,DIMENSION_SIZE,iostat=ios)
        if(ios.ne.0) call FFRABORT(1,'ERR: [dimen.h_serial] read error')
        close(ifl)
! --- Only for serial cpu
          NVRTX_L =NNVRTX
          NCELL_L =NNCELL
          NCELB_L =NNCELB
          NFACE_L =NNFACE
          NEDGE_L =NNEDGE
          NCVIN_L =NNCVIN
          NCV_L   =NNCV
          NALLCV_L=NNALLCV
          NCVFAC_L=NNCVFAC

          NSSFBC_L=NNSSFBC
          NSSFBC_OLD_L=NNSSFBC_OLD

          NCOMP_L =NNCOMP
          NCOMP_SUF_L =NNCOMP_SUF
          NRANS_L =NNRANS
          NPOTN_L =NNPOTN

          NMAT_L  =NNMAT
          IEMAX_L =NNIE

          NBCTOT_L=NNBCTOT
          NBCWAL_L=NNBCWAL
          NBCSLD_L=NNBCSLD
          NPHASE_L =NNPHASE
!
          NVRTX =NVRTX_L
          NCELL =NCELL_L
          NCELB =NCELB_L
          NFACE =NFACE_L
          NEDGE =NEDGE_L
          NCVIN =NCVIN_L
          NCV   =NCV_L
          NALLCV=NALLCV_L

          NCVFAC=NCVFAC_L
          NSSFBC=NSSFBC_L
          NSSFBC_OLD=NSSFBC_OLD_L
          NCOMP =NCOMP_L
          NCOMP_SUF =NCOMP_SUF_L

          NRANS =NRANS_L
          NPOTN =NPOTN_L
          NMAT  =NMAT_L
          IEMAX =IEMAX_L

          NBCWAL=NBCWAL_L
          NBCSLD=NBCSLD_L
          
          NBCTOT=NBCTOT_L
          NCOMPALL=NNCOMP_SUF+NNCOMP
          MXCOMPALL=NNCOMP_SUF+NNCOMP
          NPHASE=NPHASE_L
!
          if(ivector==1.or.iBODY==3) then  !.or.iprtcle==1) then
            MXCV_V=MXALLCV
            MXBND_V=MAXIE
          else
            MXCV_V=1
            MXBND_V=1
          endif
!
!
          if(iMHD==1.or.iMHD==2) then
            MXCV_MHD=MXALLCV
            MXBND_MHD=MXSSFBC
            MXFACE_MHD=MXCVFAC
          else
            MXCV_MHD=1
            MXBND_MHD=1
            MXFACE_MHD=1
          endif
!
          if(iMvmsh/=0) then
            MXVRTX_m=MXVRTX
            MXCELL_m=MXCELB!MXCELL
            MXFACE_m=MXFACE
            MXEDGE_m=MXEDGE
            MXSSFBC_m=MXSSFBC
            MXSSFBC_old_m=MXSSFBC_old
          else
            MXVRTX_m=1
            MXCELL_m=1
            MXFACE_m=1
            MXEDGE_m=1
            MXSSFBC_m=1
            MXSSFBC_old_m=1
          endif
!
          if(iBODY==1.or.iBODY==2.or.iBODY==5) then
            MXCV_B=MXALLCV
          else
            MXCV_B=1
          endif
!
          if(iBODY==4) then
            MXCVFAC_B=MXCVFAC
          else
            MXCVFAC_B=1
          endif

          if(FLMLT==1) then
            MXCV_F=MXCV
            MXCVFAC_F=MXCVFAC
          else
            MXCV_F=1
            MXCVFAC_F=1
          endif
!
          if(KE2P>0.or.iCAVIT==1.or.iPEFC>0) then
            MXALLCVC=MXALLCV
          else
            MXALLCVC=1
          endif
!
          if(ifld>0) then
            NFLID=ifld
            MXCV_D=MXALLCV
          else
            NFLID=1
            MXCV_D=1
          endif
!
          if(iprtcle==1) then
            MXPRT=MXPRTX
            MEVENT=500
            MEP=MAXIE
            MXALLCV_P=MXALLCV
            MXCVFAC_P=MXCVFAC
            MXCV_WDR=1
            N_inj=max(1,iinj)
          elseif(iprtcle==2) then
            MXPRT=MXPRTX
            MEVENT=500
            MEP=MAXIE
            MXALLCV_P=MXALLCV
            MXCVFAC_P=MXCVFAC
            MXCV_WDR=MXCV
            N_inj=max(1,iinj)
          else
            MXPRT=1
            MEVENT=1
            MEP=1
            MXALLCV_P=1
            MXCVFAC_P=1
            MXCV_WDR=1
            N_inj=max(1,iinj)
          endif
!
          if(iprtcle==1) then
            MXCV_VP=MXALLCV
            MXBND_VP=MAXIE
          elseif(iprtcle==2) then
            MXCV_VP=MXALLCV
            MXBND_VP=MAXIE
          else
            MXCV_VP=1
            MXBND_VP=1
          endif

!
          if(ivof==1)then
            MXSSFBC_VOF=MXSSFBC
            MXALLCV_VOF=MXALLCV
          else
            MXSSFBC_VOF=1
            MXALLCV_VOF=1
          endif
!
      endif

!
      return
!
      end subroutine n_size
!
!=================================================
      subroutine n_size_hpc
     & (NPE,my_rank,Dmns,ivector,iMHD,iMvmsh,ivof,
     &  iBODY,FLMLT,iCAVIT,iPEFC,iprtcle,iinj,ifld,MXPRTX,ICALL)
!=================================================
      implicit none
!
! --- [dummy argument]
!
      integer,intent(in)  :: NPE,my_rank,ICALL,ivector,MXPRTX,
     6                       iMHD,iMvmsh,iBODY,ivof
      integer,intent(in)  :: FLMLT,iCAVIT,ifld,iPEFC,iprtcle,iinj
      character(len=80),intent(in) :: Dmns
!
! --- [local entities]
!
      integer             :: ifl,ios=0
! --- Only for HPC!
      ifl=11!+my_rank
      open(ifl,file=Dmns,form='formatted',status='UNKNOWN',iostat=ios)
      if(ios.ne.0) call FFRABORT(1,'ERR: [dimen.parm] open error')
      rewind ifl
      read(ifl,DIMENSION_SIZE,iostat=ios)
      if(ios.ne.0) call FFRABORT(1,'ERR: [dimen.parm] read error')
      close(ifl)
!
      NVRTX_L =NNVRTX
      NCELL_L =NNCELL
      NCELB_L =NNCELB
      NFACE_L =NNFACE
      NEDGE_L =NNEDGE
      NCVIN_L =NNCVIN
      NCV_L   =NNCV
      NALLCV_L=NNALLCV
      NCVFAC_L=NNCVFAC

      NSSFBC_L=NNSSFBC
      NSSFBC_OLD_L=NNSSFBC_OLD

      NCOMP_L =NNCOMP
      NCOMP_SUF_L =NNCOMP_SUF
      NRANS_L =NNRANS
      NPOTN_L =NNPOTN

      NMAT_L  =NNMAT
      IEMAX_L =NNIE

      NBCTOT_L=NNBCTOT
      NBCWAL_L=NNBCWAL
      NBCSLD_L=NNBCSLD
      NPHASE_L =NNPHASE
!
      NVRTX =NVRTX_L
      NCELL =NCELL_L
      NCELB =NCELB_L
      NFACE =NFACE_L
      NEDGE =NEDGE_L
      NCVIN =NCVIN_L
      NCV   =NCV_L
      NALLCV=NALLCV_L

      NCVFAC=NCVFAC_L
      NSSFBC=NSSFBC_L
      NSSFBC_OLD=NSSFBC_OLD_L
      NCOMP =NCOMP_L
      NCOMP_SUF =NCOMP_SUF_L

      NRANS =NRANS_L
      NPOTN =NPOTN_L
      NMAT  =NMAT_L
      IEMAX =IEMAX_L

      NBCWAL=NBCWAL_L
      NBCSLD=NBCSLD_L

      NBCTOT=NBCTOT_L
      NCOMPALL=NNCOMP_SUF+NNCOMP
      MXCOMPALL=NNCOMP_SUF+NNCOMP
      NPHASE=NPHASE_L
!
      if(ivector==1.or.iBODY==3) then
        MXCV_V=MXALLCV
        MXBND_V=MAXIE
      else
        MXCV_V=1
        MXBND_V=1
      endif
!
!
      if(iMHD==1.or.iMHD==2) then
        MXCV_MHD=MXALLCV
        MXBND_MHD=MXSSFBC
        MXFACE_MHD=MXCVFAC
      else
        MXCV_MHD=1
        MXBND_MHD=1
        MXFACE_MHD=1
      endif
!
      if(iMvmsh/=0) then
        MXVRTX_m=MXVRTX
        MXCELL_m=MXCELB   !MXCELL
        MXFACE_m=MXFACE
        MXEDGE_m=MXEDGE
        MXSSFBC_m=MXSSFBC
        MXSSFBC_old_m=MXSSFBC_old
      else
        MXVRTX_m=1
        MXCELL_m=1
        MXFACE_m=1
        MXEDGE_m=1
        MXSSFBC_m=1
        MXSSFBC_old_m=1
      endif
!
      if(iBODY==1.or.iBODY==2.or.iBODY==5) then
        MXCV_B=MXALLCV
      else
        MXCV_B=1
      endif
!
      if(iBODY==4) then
        MXCVFAC_B=MXCVFAC
      else
        MXCVFAC_B=1
      endif
!
      if(FLMLT==1) then
        MXCV_F=MXCV
        MXCVFAC_F=MXCVFAC
      else
        MXCV_F=1
        MXCVFAC_F=1
      endif
!
      if(KE2P>0.or.iCAVIT==1.or.iPEFC>0) then
        MXALLCVC=MXALLCV
      else
        MXALLCVC=1
      endif
!
      if(ifld>0) then
        NFLID=ifld
        MXCV_D=MXALLCV
      else
        NFLID=1
        MXCV_D=1
      endif
!
      if(iprtcle==1) then
        MXPRT=MXPRTX
        MEVENT=500
        MEP=MAXIE
        MXALLCV_P=MXALLCV
        MXCVFAC_P=MXCVFAC
        MXCV_WDR=1
        N_inj=max(1,iinj)
      elseif(iprtcle==2) then
        MXPRT=MXPRTX
        MEVENT=500
        MEP=MAXIE
        MXALLCV_P=MXALLCV
        MXCVFAC_P=MXCVFAC
        MXCV_WDR=MXCV
        N_inj=max(1,iinj)
      else
        MXPRT=1
        MEP=1
        MEVENT=1
        MXALLCV_P=1
        MXCVFAC_P=1
        MXCV_WDR=1
        N_inj=max(1,iinj)
      endif
!
      if(iprtcle==1) then
        MXCV_VP=MXALLCV
        MXBND_VP=MAXIE
      elseif(iprtcle==2) then
        MXCV_VP=MXALLCV
        MXBND_VP=MAXIE
      else
        MXCV_VP=1
        MXBND_VP=1
      endif
!
      if(ivof==1)then
        MXSSFBC_VOF=MXSSFBC
        MXALLCV_VOF=MXALLCV
      else
        MXSSFBC_VOF=1
        MXALLCV_VOF=1
      endif
!
      return
!
      end subroutine n_size_hpc
!
      end module module_dimension
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      MODULE MODULE_ARRAY_FFLOW
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!---------------------
! --- CV CONNECTIVITY
!---------------------
      INTEGER,allocatable :: LVEDGE    (:,:)
      INTEGER,allocatable :: LVRTCV    (:)
      INTEGER,allocatable :: LFUTAU    (:)
!---------------------
! --- MATERIAL ARRAY
!---------------------
      INTEGER,allocatable :: MAT_NO    (:)
      INTEGER,allocatable :: MAT_CV    (:)
      INTEGER,allocatable :: LCV_CV    (:)
      INTEGER,allocatable :: MAT_INDEX (:)
      INTEGER,allocatable :: MAT_CVEXT (:)
      INTEGER,allocatable :: MAT_DCIDX (:)
      INTEGER,allocatable :: MAT_CFIDX (:)
      INTEGER,allocatable :: Pstart    (:)
      logical,allocatable :: mat_cal   (:)
!-------------------------------------------------
! --- CV'S BC Physical and CONNECTIVITY parameters
!--------------------------------------------------
      INTEGER,allocatable :: LBC_SSF   (:)
      INTEGER,allocatable :: LCYCSF    (:)
      INTEGER,allocatable :: locmsh    (:)
!
      INTEGER,allocatable :: LCYCOLD   (:)
      REAL(8),allocatable :: wifsld    (:)
      REAL(8),allocatable :: OPPANG    (:)
!
      REAL(8),allocatable :: FRSTCV   (:)
      REAL(8),allocatable :: DISINL   (:)
!
! --- BC LIST
!
!-------------------------------
! --- CV MATRIX
!-------------------------------
      REAL(8),allocatable :: SFAREA   (:,:)
      REAL(8),allocatable :: SFCENT   (:,:)
      REAL(8),allocatable :: CVVOLM   (:)
      REAL(8),allocatable :: CVVOL0   (:)
      REAL(8),allocatable :: CVCENT   (:,:)
      REAL(8),allocatable :: WIFACE   (:)
      REAL(8),allocatable :: DISALL   (:)
!-------------------------------
! --- MOVING MESH PARAMETERS
!-------------------------------
      REAL(8),allocatable  :: XTA      (:)
!-------------------------------
! --- CV'S PHYSICAL ARRAIES:
!-------------------------------
      REAL(8),allocatable :: VEL      (:,:,:)
      REAL(8),allocatable :: VELF     (:,:,:)
      REAL(8),allocatable :: RVA      (:,:)
      REAL(8),allocatable :: RHO      (:,:)
      REAL(8),allocatable :: PRS      (:,:)
      REAL(8),allocatable :: PP0      (:,:)
      REAL(8),allocatable :: DVDD     (:,:,:)
      REAL(8),allocatable :: RMUT     (:)
      REAL(8),allocatable :: RMU      (:)
      real(8),allocatable :: hhh      (:,:)
      real(8),allocatable :: cp_w     (:)
      real(8),allocatable :: cr_w     (:)
      REAL(8),allocatable :: RMX_w    (:)
      REAL(8),allocatable :: GRDC_W   (:,:,:)
      REAL(8),allocatable :: DIAG_W   (:)
      REAL(8),allocatable :: rvx_W    (:,:)
      REAL(8),allocatable :: rvd_W    (:,:)
      INTEGER,allocatable :: KDBT_W    (:)
      INTEGER,allocatable :: KDBY_W    (:)
!-------------------------------
! --- RANS model:
!-------------------------------
      REAL(8),allocatable  :: AKS      (:,:,:)
!-------------------------------
! --- Euler 2 Phase (ieul2ph=1)
!-------------------------------
!      REAL(8),allocatable  :: VEL2     (:,:,:,:)
!      REAL(8),allocatable  :: tmp2     (:,:,:)
!      REAL(8),allocatable  :: RVA2     (:,:,:)

!      REAL(8),allocatable  :: RHO2     (:,:,:)
!      REAL(8),allocatable  :: RMUT2    (:,:)
!      REAL(8),allocatable  :: RMU2     (:,:)
!      REAL(8),allocatable  :: YYS2     (:,:,:,:)
!      REAL(8),allocatable  :: RDS2     (:,:,:)
!      REAL(8),allocatable  :: WDOT2    (:,:,:)
!      REAL(8),allocatable  :: RMD2     (:,:)
!      REAL(8),allocatable  :: ccc2     (:,:)
!      REAL(8),allocatable  :: dvdd2    (:,:,:,:)
!
!      REAL(8),allocatable  :: UTAU2    (:,:)
!      REAL(8),allocatable  :: VELBND2  (:,:,:)
!      REAL(8),allocatable  :: TMPBND2  (:,:)
!      REAL(8),allocatable  :: YYSBND2  (:,:,:)
!      REAL(8),allocatable  :: HTCBND2  (:,:)
!      REAL(8),allocatable  :: MTCBND2  (:,:)
!      REAL(8),allocatable  :: RADBND2  (:,:)
!
!
      REAL(8),allocatable  :: VEL2     (:,:,:)
      REAL(8),allocatable  :: tmp2     (:,:)
      REAL(8),allocatable  :: RVA2     (:,:)
      REAL(8),allocatable  :: PRS2     (:,:)

      REAL(8),allocatable  :: RHO2     (:,:)
      REAL(8),allocatable  :: RMUT2    (:)
      REAL(8),allocatable  :: RMU2     (:)
      REAL(8),allocatable  :: YYS2     (:,:,:)
      REAL(8),allocatable  :: RDS2     (:,:)
      REAL(8),allocatable  :: WDOT2    (:,:)
      REAL(8),allocatable  :: RMD2     (:)
      REAL(8),allocatable  :: ccc2     (:)
      REAL(8),allocatable  :: dvdd2    (:,:,:)
      real(8),allocatable ::  hhh2     (:,:)
      
!
      REAL(8),allocatable  :: UTAU2    (:)
      REAL(8),allocatable  :: VELBND2  (:,:)
      REAL(8),allocatable  :: TMPBND2  (:)
      REAL(8),allocatable  :: YYSBND2  (:,:)
      REAL(8),allocatable  :: HTCBND2  (:)
      REAL(8),allocatable  :: MTCBND2  (:)
      REAL(8),allocatable  :: RADBND2  (:)
!
!
      REAL(8),allocatable  :: HEATL_C  (:)
      REAL(8),allocatable  :: HEATG_C  (:)
      REAL(8),allocatable  :: HLS      (:)
      REAL(8),allocatable  :: HGS      (:)
      REAL(8),allocatable  :: TMP_SAT  (:)
      REAL(8),allocatable  :: FRIC_C   (:)
      REAL(8),allocatable  :: MASS_I   (:)
!MXCVFAC2
!      REAL*8  :: RVAA1     (  :) !! rho*v*area*alpha1(phase 1)
!      REAL*8  :: RVAA2     (  :) !! rho*v*area*alpha2(phase 2)
!---------------------------
! --- Surface reaction
!---------------------------
      REAL(8),allocatable  :: SUFBND   (:,:)
      REAL(8),allocatable  :: SDOT     (:,:)
      REAL(8),allocatable  :: MOLFRC   (:,:,:)
      REAL(8),allocatable  :: MOLCEN   (:,:)
      REAL(8),allocatable  :: SITDEN   (:,:,:)
      REAL(8),allocatable  :: BLKTHK   (:,:,:)
!

!
!-------------------------------
! --- speciese
!-------------------------------
!
      REAL(8),allocatable  :: TMP      (:,:)
      REAL(8),allocatable  :: YYS      (:,:,:)
      REAL(8),allocatable  :: RDS      (:,:)
      REAL(8),allocatable  :: WDOT     (:,:)
      REAL(8),allocatable  :: RMD      (:)
!------------------------------
! --- speciese (working array)
!------------------------------
      real(8),allocatable  :: hhs_w    (:,:)
      real(8),allocatable  :: cps_w    (:,:)
!-------------------------------
! --- Sound speed (compressible)
!-------------------------------
      REAL(8),allocatable  :: CCC      (:)
!-------------------------------
! --- wall friction velocity
!-------------------------------
      REAL(8),allocatable  :: UTAU     (:)
!-------------------------------
! --- POTENTIAL
!-------------------------------
      REAL(8),allocatable  :: POTNAL   (:,:,:)
      REAL(8),allocatable  :: POFLUX   (:,:)
      REAL(8),allocatable  :: POTNBC   (:,:)
!      REAL(8),allocatable  :: CURENT   (:,:)
!
      REAL(8),allocatable  :: VELBND   (:,:)
      REAL(8),allocatable  :: PRSBND   (:)
      REAL(8),allocatable  :: TMPBND   (:)
      REAL(8),allocatable  :: YYSBND   (:,:)
!
      REAL(8),allocatable  :: AKSBND   (:,:)
      REAL(8),allocatable  :: HTCBND   (:)
      REAL(8),allocatable  :: HTFLUX   (:)
      REAL(8),allocatable  :: MTCBND   (:)
      REAL(8),allocatable  :: RADBND   (:)
!
!-------------------------------
! --- STATISTICS AVERAGE (LES)
!-------------------------------
      REAL(8),allocatable  :: VLASUM   (:)
!--------------------------
! --- Control parameter
!--------------------------
      integer,allocatable :: iptfix(:)
      integer,allocatable :: ipfix(:)
      integer,allocatable :: itery(:)  ,
     &                       iterk(:),
     &                       iterR(:)
      real(8),allocatable :: repsy(:)  ,
     &                       repsk(:),
     &                       repsR(:)
      real(8),allocatable :: aepsy(:)  ,
     &                       aepsk(:),
     &                       aepsR(:)
      real(8),allocatable :: erry (:)  , 
     &                       errk (:),
     &                       errR (:)
      real(8),allocatable :: errsdy(:) ,
     &                       errsdk(:)
!
! --- Euler two phase
!
!      integer,allocatable  :: iterv2(:,:)
!      real(8),allocatable  :: repsv2(:,:),
!     &                        aepsv2(:,:),
!     &                        errv2(:,:),
!     &                        errsdv2(:,:)
!      integer,allocatable  :: itery2(:,:)
!      real(8),allocatable  :: repsy2(:,:)
!      real(8),allocatable  :: aepsy2(:,:)
!      real(8),allocatable  :: erry2 (:,:)
!      real(8),allocatable  :: errsdy2(:,:)
!
      integer,allocatable  :: iterv2(:)
      real(8),allocatable  :: repsv2(:),
     &                        aepsv2(:),
     &                        errv2(:),
     &                        errsdv2(:)
      integer,allocatable  :: itery2(:)
      real(8),allocatable  :: repsy2(:)
      real(8),allocatable  :: aepsy2(:)
      real(8),allocatable  :: erry2 (:)
      real(8),allocatable  :: errsdy2(:)
!
      real(8),allocatable  :: reps_POTN(:),aeps_POTN(:),err_POTN(:)
      integer,allocatable  :: iter_POTN(:)
!
      
!
      END MODULE MODULE_ARRAY_FFLOW
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_metrix
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- vectorozation 
! 
      integer,save,allocatable  :: vctr(:,:)
      integer,save,allocatable  :: vctrp(:,:)
      real(8),save,allocatable  :: FIELD_U(:,:)
      real(8),save,allocatable  :: J_SRC(:,:)
      real(8),save,allocatable  :: TH_DIFF(:,:)
! --- static(for CG solver)
      integer,save,allocatable  :: msk (:)
      integer,save,allocatable  :: IAL (:,:)
      integer,save,allocatable  :: INL (:)
      integer,save,allocatable  :: INU (:)
      real(8),save,allocatable  :: aalw(:,:)
      real(8),save,allocatable  :: bb  (:)
! --- static (CG solver, BiCGstab)
      real(8),allocatable       :: W1K1(:)
      real(8),allocatable       :: W1K2(:)
      real(8),allocatable       :: W1K3(:)
      real(8),allocatable       :: W1K4(:)
      real(8),allocatable       :: W1K5(:)
      real(8),allocatable       :: W1K6(:)
      real(8),allocatable       :: W1K7(:)
      real(8),allocatable       :: W1K8(:)
      real(8),allocatable       :: W1K9(:)
      real(8),allocatable       :: W1K10(:)
      real(8),allocatable       :: W1K11(:)
!
! --- Dynamic (for CG solver) 
!
      real(8),allocatable       :: W2K1(:,:)
      real(8),allocatable       :: W2K2(:,:)
!
! --- Dynamic (only in hpcmw.f)
!
      integer,allocatable       :: iwork1(:),iwork2(:),iwork3(:)
      integer,allocatable       :: iwork4(:),iwork5(:),
     &                             iwork6(:),iwork7(:)
      integer,allocatable       :: IW2K1(:,:),IW2K2(:,:),
     &                             IW2K_VECT(:,:)
! --- Dynamic (for ***_admin subroutine)
      real(8),allocatable       :: d2work1(:,:),
     &                             d2work2(:,:),
     &                             d2work3(:,:),
     &                             d2work4(:,:),
     &                             d2work5(:,:),
     &                             d2vect(:,:),DW2K_VECT(:,:)
      real(8),allocatable       :: d1work1(:),d1work2(:)
      real(8),allocatable       :: tmpsld(:,:)
!      logical,save,allocatable  :: mask(:)
!
! --- dummy
!
      real(8),allocatable       :: ys(:),
     &                             rcomp(:),mu_i(:),k_i(:),
     &                             dum_c1(:),
     &                             difu(:),thmdif(:),
     &                             eva_comp(:),
     &                             sinj_comp(:),
     &                             aks_temp(:),
     &                             erryys(:),
     &                             ahk(:,:),
     &                             ahk2(:,:),
     &                             acpk_w(:,:),
     &                             errrans(:),
     &                             dvdc(:),
     &                             dpdt(:),
     &                             rsdc(:),
     &                             vsm(:),
     &                             rsm_prim(:),
     &                             rsm(:),
     &                             tsm(:),
     %  QINL(:),QOUT(:),OFFSET(:),SINL(:),SOUT(:),
     &  c1(:),c2(:),c3(:),rmaxs(:),rmaxs0(:),
     6  alpha(:),beta(:),t(:),s(:),u(:),v(:),omega(:)
      integer,allocatable       :: matsld(:,:)
      integer,allocatable       :: ipx(:),NCOL1(:),NCOL2(:)
      real(8),allocatable       :: aax(:)
! --- HPC send_recv
      real*8 ,save, allocatable :: WS(:),WR(:)
      integer,save, allocatable :: sta1(:,:)
      integer,save, allocatable :: sta2(:,:)
      integer,save, allocatable :: req1(:  )
      integer,save, allocatable :: req2(:  )
! --- MHD
      real*8 ,save,allocatable  :: MHD_A   (:,:,:,:)
      real*8 ,save,allocatable  :: MHD_FAI (:,:,:)
      real*8 ,save,allocatable  :: MHD_RVA (:)
      real*8 ,save,allocatable  :: MHD_CRT0(:,:,:)
      real*8 ,save,allocatable  :: MHD_FAI0 (:,:,:)
      real*8 ,save,allocatable  :: MHD_CRNT(:,:,:)
      real*8 ,save,allocatable  :: SIGMA   (:)
      real*8 ,save,allocatable  :: MHDbnd  (:,:,:)
      real*8 ,save,allocatable  :: tempMHD (:,:,:)
      integer,save, allocatable :: mat_coil(:)
!      integer,save, allocatable :: GaugeC0(:)
! --- MOVE MESH
      integer,save,allocatable  :: movflg(:)
      real*8 ,save,allocatable  :: cord  (:,:)!
      integer,save,allocatable  :: lacell(:)
      integer,save,allocatable  :: lvcell(:,:)

      integer,save,allocatable  :: lvface(:,:)!lvface(4,mface)
      integer,save,allocatable  :: lbface(:,:)!lbface(2,mface)
      integer,save,allocatable  :: lcface(:,:)!lcface(2,mface)
      integer,save,allocatable  :: lfcell(:,:) !lfcell(7,mcell)
      integer,save,allocatable  :: LEFACE(:,:)!LEFACE(5,MXFACE_m)
      integer,save,allocatable  :: listbc(:,:)!listbc(4,mssfbc)
      real*8 ,save,allocatable  :: area  (:,:)!(4,mface)
      real*8 ,save,allocatable  :: volume(:)  !(  mcell)
      real*8 ,save,allocatable  :: gface (:,:)!(3,mface)
      real*8 ,save,allocatable  :: gcell (:,:)!(3,mcell)
!
! --- body force 
!
      real*8 ,save,allocatable  :: bdyfrc(:,:,:) !(:,3)
      real*8 ,save,allocatable  :: bdyf(:,:)
! --- flamelet
      real*8 ,save,allocatable  :: vflc(:),yplusf(:)
      real*8 ,save,allocatable  :: rva_s(:)
!
      logical,save,allocatable  :: lreac(:)
      logical,save,allocatable  :: lreacP(:)
!
! --- VOF      
!
      real*8 ,save,allocatable  :: t_dir(:,:)
      REAL*8 ,save,ALLOCATABLE  :: ANGBND(:) 
!
      REAL*8 ,save,ALLOCATABLE  :: STA_SAV(:,:)
!
      integer,save,allocatable  :: SHUTFL(:)
!
      integer,save ::         iter_gag(3,2)
!
      REAL(8),allocatable  :: WDRBND   (:,:)
!
      REAL(8),allocatable  :: rghnesx   (:)
!
      integer,allocatable  :: multiR(:)
!
      end module module_metrix
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_Dynamic_SGS
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      REAL(8),allocatable :: lij(:,:),mij(:,:)
!
      end module module_Dynamic_SGS
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      MODULE module_trace 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      USE MODULE_DIMENSION
!      integer, parameter :: MP= 100000, MEP = 8
!--------------------------------------------------------------
!   <<<<<<< ELEMENT'S CONNECTIVITY SPECIFYING LISTS >>>>>>>
!--------------------------------------------------------------
      integer,save,allocatable :: NEIGHBOR(:,:)		
!      integer,save,allocatable :: NVCRSIGN(:,:) 
!      integer,save,allocatable :: LFCV(:,:)
!-------------------------------------------------------------
! --- <<<<<<< ARRAYS FOR PARTICLES' PROPERTY >>>>>>>
!-------------------------------------------------------------
      real(8) ,save,allocatable :: CC(:)
      real(8) ,save,allocatable :: HEIGHT0(:)
      real(8) ,save,allocatable :: HEIGHT1(:)
      real(8) ,save,allocatable :: VELN(:)
      
      real(8) ,save,allocatable :: XYZP(:,:)
      real(8) ,save,allocatable :: UVWP(:,:)
      real(8) ,save,allocatable :: PARCEL(:)
      real(8) ,save,allocatable :: DIA(:)
      real(8) ,save,allocatable :: PMASS(:,:)
      real(8) ,save,allocatable :: HIS_P(:)
      real(8) ,save,allocatable :: PYS(:,:,:)
      real(8) ,save,allocatable :: HEIGHT(:,:)
      real(8) ,save,allocatable :: FU(:,:)
      real(8) ,save,allocatable :: FV(:,:)
      real(8) ,save,allocatable :: FW(:,:)
      real(8) ,save,allocatable :: RHOP(:)
      real(8) ,save,allocatable :: tmpp(:,:)
      integer,save,allocatable  :: INSIDE(:)
      integer,save,allocatable  :: JOUT(:)
      real(8) ,save,allocatable :: DNORM(:)
      real(8) ,save,allocatable :: MOMCELL(:,:)
      real(8) ,save,allocatable :: EVAP_Y(:,:)
      real(8) ,save,allocatable :: EVAP_H(:,:)
      integer,save,allocatable  :: IPCELL(:),NEIBCELL(:)
      real(8) ,save,allocatable :: DYDTP(:,:)
      real(8) ,save,allocatable :: DHDTP(:,:)
!
! --- dum array
!
      real(8) ,save,allocatable :: WK1_P(:,:)
      real(8) ,save,allocatable :: PRT_VF(:)
      integer,save,allocatable  :: IWK_P(:)
!
!
!      INTEGER, PARAMETER :: MEVEVT = 500
!
      end MODULE module_trace 
      
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_initial 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      real(8) ,save,allocatable  :: t0(:),y0(:,:)
      real(8) ,save,allocatable  :: p0(:)
      real(8) ,save,allocatable  :: u0(:),v0(:),w0(:)
      real(8) ,save,allocatable  :: vsgm0(:,:),rho0(:)
      real(8) ,save,allocatable  :: aks0(:,:)
      
!
      real(8) ,save,allocatable  :: t02(:),y02(:,:)
      real(8) ,save,allocatable  :: u02(:),v02(:),w02(:)
      real(8) ,save,allocatable  :: rho02(:)
!
!      integer ,save,allocatable  :: vap_no(:) 
      real(8) ,save,allocatable  :: Ke_FC(:),Kc_FC(:),Pv_FC(:),
     &                              Kap_FC(:),sigma_FC(:),contang_FC(:),
     &                              elecCondFC(:),ionCondFC(:)
!      integer ,save,allocatable  :: osmoticFC(:),surf_tenFC(:)
      integer ,save,allocatable  :: surf_tenFC(:)
!
      integer ,save,allocatable  :: icvprs(:) 
!
      real(8) ,save,allocatable  :: fai_iMHD(:,:),A_iMHD(:,:,:)
      real(8) ,save,allocatable  :: eps_iMHD(:),eps0_iMHD(:)
      real(8) ,save,allocatable  :: mu_iMHD(:),mu0_iMHD(:)
      real(8) ,save,allocatable  :: sigma_iMHD(:),sigma0_iMHD(:)
      real(8) ,save,allocatable  :: Current_iMHD(:,:,:)
      real(8) ,save,allocatable  :: OMEGA_iMHD(:)
!
      real(8) ,save,allocatable  :: potn0(:,:),potn_c(:,:)
!
      end module module_initial
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_radsxf
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      INTEGER :: NRD,MRD,NRD0,MRD0,MID
      INTEGER :: NumTran,NumTranIn,NDatExpt,NPALL,NMatBoun
      CHARACTER*20 :: StackOpt
      REAL*8 :: AveVs(10)
      REAL*8, ALLOCATABLE :: VJOBMD(:,:)
!--------------------------------------
!(:,1)	:	Area / Volume
!(:,2)  :	Property
!(:,3)  :	QW
!(:,4)  :	Timag (for communication)
!--------------------------------------
      INTEGER,ALLOCATABLE:: NJOBMD(:,:),ID(:)
      REAL*8, ALLOCATABLE:: DivW(:),DivA(:,:)
      
      REAL*8, ALLOCATABLE:: CoefToRD(:),RaToExpt(:),EBcomm(:)
      INTEGER,ALLOCATABLE:: IndexToRD(:,:),IndexToExpt(:,:),
     &                      MatBounEnd(:),ExptID(:),WK(:)
!--------------------------------------
!----- only for benchmark
!--------------------------------------
      INTEGER,ALLOCATABLE :: WDAT(:,:),WPNUM(:,:)
      INTEGER::	NWDat
!--------------------------------------
!temperature / heat flux
!--------------------------------------
      REAL*8, ALLOCATABLE::  !Timag(:),RadBB(:),TPtc(:),
     &       QSUM(:,:),EB4(:),RadHeat(:),RadHeatFlux(:)
      REAL*8, ALLOCATABLE::  radinten(:,:),sumcoef(:)
!--------------------------------------
!------ Exchange-Area (Configure-factor)
!--------------------------------------
      INTEGER, ALLOCATABLE:: RDId(:), RDN1(:),RDN2(:),RDIndex(:)
      REAL*8, ALLOCATABLE::  RDValue(:)!,RDV2(:)
!--------------------------------------
!RDN1	 :   off-diagonal-line elements in Lower-part
!RDN2	 :   all off-diagonal-line elements
!--------------------------------------
	
!      REAL*8,ALLOCATABLE::  ExtinCoef(:),Albedo(:)
      REAL*8,ALLOCATABLE::  GWAbsorb(:),PAbsorb(:),PScatter(:)
      REAL*8,ALLOCATABLE::  P_rad(:,:)
!
      REAL*8,ALLOCATABLE::  WRad(:), FRAD(:),raddum(:)
!
!------ DUMMY ARRAYS
!
      end module module_radsxf
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_radgroup
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	
      INTEGER:: IFlagGroup,NumGPTran
!      real*8,allocatable:: GpProp(:)
!      INTEGER,allocatable:: GpCNum(:),WSID(:,:)
!      REAL*8, ALLOCATABLE::  COEF(:)
!      INTEGER,ALLOCATABLE::  NIndexToGp(:,:)
      end module module_radgroup
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      module module_radwork
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	
!      real*8,allocatable:: VWK1 (:),VWK2 (:),VWK3 (:),VWK4 (:)
!      real*8,allocatable:: VWKM1(:,:),VWKM2(:,:)
!      real*8,allocatable:: VWK3D(:,:,:)
!      INTEGER,allocatable:: WK1(:),WK2(:),WK3(:),WK4(:),WK5(:),WK6(:)
!      INTEGER,allocatable:: WKM1(:,:),WKM2(:,:),WKM3(:,:)
!      end module module_radwork
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program FRONTFLOW
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. compresible/incompresible fluid flow solver
!     on unstructured grid
!
!  2. Vertex Center Finite Volume Method With "Edge-Base" Structure
!
!  3. HPC by MPI
!
!=====================================================================
!
! --- [module arguments]
!
      use module_dimension
      use module_ARRAY_FFLOW
      use module_metrix
      use module_initial
      USE module_trace 
      USE module_radsxf
!
      implicit none
!
! --- [local entities]
!
      integer             :: ICALL,license=0
      integer             :: ierror=0,NPE=1,my_rank=0
      integer             :: ierr=0
      character (len=80)  :: Dmns
      
!
      integer             :: iterv(3)
      real(8)             :: repsv(3),aepsv(3),errv(3),errsdv(3)
!
      integer             :: iterp
      real(8)             :: repsp,aepsp,errp
!
      real(8)             :: reps_FAI(2) ,aeps_FAI(2) ,err_FAI(2),
     &                       reps_A(3,2),aeps_A(3,2),err_A(3,2)
      integer             :: iter_FAI(2),iter_A(3,2)
      real(8)             :: errsda(3,2),errsdfai(2)
!
!      real(8)             :: reps_POTN(MXPOTN),aeps_POTN(MXPOTN)
!
      integer,save        :: ivector=0,iMHD=0,iMvmsh=0,iBODY=0,FLMLT=0,
     &                       iCAVIT=0,iPEFC=0,ifld=0,iprtcle=0,iinj=0,
     &                       ivof=0
!
      integer,save        :: MXPRTX
!=====================================================================
! --- Start
!=====================================================================
      call unlink('STOP')
      call unlink('stop')
      call unlink('output')
!
!-----------------------------------------
! --- Intializing MPI for HPC computation 
!-----------------------------------------
      call HPC_INIT(NPE,my_rank)
      call fftime(0)
      if(my_rank.eq.0) call title
!
      ICALL=0
      call n_size(NPE,Dmns,
     &            ivector,iMHD,iMvmsh,ivof,
     &            iBODY,FLMLT,iCAVIT,iPEFC,iprtcle,iinj,
     &            ifld,MXPRTX,ICALL) 
!----------------------------------------------------------------
! --- Input initial, Boundary Condition and check dimension size
!----------------------------------------------------------------
      call INPUT_FFLOW(Dmns,ivector,iMHD,iMvmsh,ivof,iBODY,FLMLT,
     &                 iCAVIT,iPEFC,iprtcle,iinj,ifld,MXPRTX,ierror)
      if( ierror.ne.0 ) call FFRABORT(1,'ERR: in INPUT_FFLOW')
!
!
!
!---------------------------------
! --- set local N-size (dimension)
!---------------------------------
      if(NPE.eq.1) then
        ICALL=1
        call n_size(NPE,Dmns,
     &     ivector,iMHD,iMvmsh,ivof,
     &     iBODY,FLMLT,iCAVIT,iPEFC,iprtcle,iinj,
     &     ifld,MXPRTX,ICALL)
      elseif(NPE.gt.1) then
        ICALL=1
        call n_size_hpc
     &  (NPE,my_rank,Dmns,ivector,iMHD,iMvmsh,ivof,
     &   iBODY,FLMLT,iCAVIT,iPEFC,iprtcle,iinj,ifld,MXPRTX,ICALL)
      else
        call FFRABORT(1,'ERR: IN MPI INITIALIZATION FRONTFLOW')
      endif
!
!=====================================================================
! --- Allocating all array
!
      ALLOCATE (ipfix(MXMAT)                  ,stat=ierr)
      ALLOCATE (iptfix(MXMAT)                 ,stat=ierr)
      ALLOCATE (itery(0:ncomp) ,iterk(mxrans) ,stat=ierr)
      ALLOCATE (repsy(0:ncomp) ,repsk(mxrans) ,stat=ierr)
      ALLOCATE (aepsy(0:ncomp) ,aepsk(mxrans) ,stat=ierr)
      ALLOCATE (erry(0:ncomp)  ,errk(mxrans)  ,stat=ierr)
      ALLOCATE (errsdy(0:ncomp),errsdk(mxrans),stat=ierr)
!
! --- Euler two phase
!
!      ALLOCATE (itery2(0:ncomp,NPHS)        ,stat=ierr)
!      ALLOCATE (repsy2(0:ncomp,NPHS)        ,stat=ierr)
!      ALLOCATE (aepsy2(0:ncomp,NPHS)        ,stat=ierr)
!      ALLOCATE (erry2(0:ncomp,NPHS)         ,stat=ierr)
!      ALLOCATE (errsdy2(0:ncomp,NPHS)       ,stat=ierr)
!      ALLOCATE (iterv2(3,NPHS),stat=ierr)
!      ALLOCATE (repsv2(3,NPHS),stat=ierr)
!      ALLOCATE (aepsv2(3,NPHS),stat=ierr)
!      ALLOCATE (errv2(3,NPHS),stat=ierr)
!      ALLOCATE (errsdv2(3,NPHS),stat=ierr)
!
      ALLOCATE (itery2(0:ncomp)        ,stat=ierr)
      ALLOCATE (repsy2(0:ncomp)        ,stat=ierr)
      ALLOCATE (aepsy2(0:ncomp)        ,stat=ierr)
      ALLOCATE (erry2(0:ncomp)         ,stat=ierr)
      ALLOCATE (errsdy2(0:ncomp)       ,stat=ierr)
      ALLOCATE (iterv2(3)              ,stat=ierr)
      ALLOCATE (repsv2(3)              ,stat=ierr)
      ALLOCATE (aepsv2(3)              ,stat=ierr)
      ALLOCATE (errv2(3)               ,stat=ierr)
      ALLOCATE (errsdv2(3)             ,stat=ierr)
!
      if(ierr/=0) call FFRABORT(1,'stop at allocating array')
!
      ALLOCATE (reps_POTN(MXPOTN)            ,stat=ierr)
      ALLOCATE (aeps_POTN(MXPOTN)            ,stat=ierr)
      ALLOCATE (err_POTN(MXPOTN)             ,stat=ierr)
      ALLOCATE (iter_POTN(MXPOTN)            ,stat=ierr)
!
      ALLOCATE (LVEDGE    (2,MXCVFAC)  ,stat=ierr)   !!!!
      ALLOCATE (LVRTCV    (MXALLCV)    ,stat=ierr)
      ALLOCATE (LFUTAU    (MXCV)       ,stat=ierr)
      ALLOCATE (MAT_NO    (0:MXMAT)    ,stat=ierr)
      ALLOCATE (MAT_CV    (MXALLCV)    ,stat=ierr)
      ALLOCATE (LCV_CV    (MXALLCV)    ,stat=ierr)
      ALLOCATE (MAT_INDEX (0:MXMAT)    ,stat=ierr)
      ALLOCATE (MAT_CVEXT (0:MXMAT)    ,stat=ierr)
      ALLOCATE (MAT_DCIDX (0:MXMAT)    ,stat=ierr)
      ALLOCATE (MAT_CFIDX (0:MXMAT)    ,stat=ierr)
      ALLOCATE (Pstart    (MXMAT)      ,stat=ierr)
      ALLOCATE (mat_cal   (0:MXMAT)    ,stat=ierr)
      ALLOCATE (LBC_SSF   (MXSSFBC)    ,stat=ierr)
      ALLOCATE (LCYCSF    (MXSSFBC)    ,stat=ierr)
      ALLOCATE (LCYCOLD   (MXSSFBC_SLD),stat=ierr)
      ALLOCATE (locmsh    (MXSSFBC)    ,stat=ierr)
      ALLOCATE (wifsld    (MXSSFBC_SLD),stat=ierr)
      ALLOCATE (OPPANG    (MXSSFBC_SLD),stat=ierr)
      ALLOCATE (FRSTCV    (MXSSFBC)    ,stat=ierr)
      ALLOCATE (DISINL    (MXSSFBC)    ,stat=ierr)
      ALLOCATE (SFAREA    (4,MXCVFAC)  ,stat=ierr)  !!!!
      ALLOCATE (SFCENT    (3,MXCVFAC)  ,stat=ierr)  !!!!
      ALLOCATE (CVVOLM    (MXALLCV)    ,stat=ierr)
      ALLOCATE (CVVOL0    (MXALLCV)    ,stat=ierr)
      ALLOCATE (CVCENT    (3,MXALLCV)  ,stat=ierr)  !!!!
      ALLOCATE (WIFACE    (MXCVFAC)    ,stat=ierr)
      ALLOCATE (DISALL    (MXCV)       ,stat=ierr)
      ALLOCATE (XTA       (MXCVFAC)    ,stat=ierr)
      ALLOCATE (VEL       (MXALLCV,3,2),stat=ierr)
      ALLOCATE (VELF      (MXCVFAC_B,3,2),stat=ierr)
      ALLOCATE (RVA       (MXCVFAC,2)  ,stat=ierr)
      ALLOCATE (RHO       (MXALLCV,2)  ,stat=ierr)
      ALLOCATE (PRS       (MXALLCV,2)  ,stat=ierr)
      ALLOCATE (PP0       (MXALLCV,2)  ,stat=ierr)
      ALLOCATE (DVDD      (MXCV,3,2)   ,stat=ierr)
      ALLOCATE (RMUT      (MXALLCV  )  ,stat=ierr)
      ALLOCATE (RMU       (MXALLCV  )  ,stat=ierr)
      ALLOCATE (hhh       (MXALLCV,2)  ,stat=ierr)
      ALLOCATE (cp_w      (MXALLCV)    ,stat=ierr)
      ALLOCATE (cr_w      (MXALLCV)    ,stat=ierr)
      ALLOCATE (RMX_w     (MXALLCV)    ,stat=ierr)
      ALLOCATE (GRDC_W    (MXALLCV,3,3),stat=ierr)
      ALLOCATE (DIAG_W    (MXALLCV)    ,stat=ierr)
      ALLOCATE (rvx_W     (MXALLCV,3)  ,stat=ierr)
      ALLOCATE (rvd_W     (MXCVFAC,3)  ,stat=ierr)
      ALLOCATE (KDBT_W    (MXCVFAC)    ,stat=ierr)
      ALLOCATE (KDBY_W    (MXCVFAC)    ,stat=ierr)
      ALLOCATE (AKS       (MXALLCVR,MXRANS,2),stat=ierr)
!
      ALLOCATE (POTNAL    (MXALLCVP,MXPOTN,2),stat=ierr)
      ALLOCATE (POFLUX    (MXCVFACP,MXPOTN),stat=ierr)
      ALLOCATE (POTNBC    (MXSSFBCP,MXPOTN),stat=ierr)
!
!      ALLOCATE (VEL2      (MXALLCV2,3,2,NPHS),stat=ierr)
!      ALLOCATE (PRS2      (MXALLCV2,2,,NPHS) ,stat=ierr)
!      ALLOCATE (tmp2      (MXALLCV2,2,NPHS)  ,stat=ierr)
!      ALLOCATE (RVA2      (MXCVFAC2,2,NPHS)  ,stat=ierr)
!      ALLOCATE (RHO2      (MXALLCVC,2,NPHS)  ,stat=ierr)
!      ALLOCATE (RMUT2     (MXALLCV2  ,NPHS)  ,stat=ierr)
!      ALLOCATE (RMU2      (MXALLCV2  ,NPHS)  ,stat=ierr)
!      ALLOCATE (YYS2      (MXALLCV2,MXCOMP,2,NPHS),stat=ierr)
!      ALLOCATE (RDS2      (MXALLCV2,MXCOMP,NPHS)  ,stat=ierr)
!      ALLOCATE (WDOT2     (MXALLCV2,MXCOMP,NPHS)  ,stat=ierr)
!      ALLOCATE (RMD2      (MXALLCV2,NPHS)         ,stat=ierr)
!      ALLOCATE (ccc2      (MXALLCV2,NPHS)         ,stat=ierr)
!      ALLOCATE (dvdd2     (MXCV2,3,2,NPHS)        ,stat=ierr)
!      ALLOCATE (UTAU2     (0:MXSSFBC2,NPHS)       ,stat=ierr)
!      ALLOCATE (VELBND2   (MXSSFBC2,3,NPHS)       ,stat=ierr)
!      ALLOCATE (TMPBND2   (MXSSFBC2,NPHS)         ,stat=ierr)
!      ALLOCATE (YYSBND2   (MXSSFBC2,MXCOMP,NPHS)  ,stat=ierr)
!      ALLOCATE (HTCBND2   (MXSSFBC2,NPHS)         ,stat=ierr)
!      ALLOCATE (MTCBND2   (MXSSFBC2,NPHS)         ,stat=ierr)
!      ALLOCATE (RADBND2   (MXSSFBC2,NPHS)         ,stat=ierr)
!
      ALLOCATE (VEL2      (MXALLCV2,3,2),stat=ierr)
      ALLOCATE (PRS2      (MXALLCV2,2),stat=ierr)
      ALLOCATE (tmp2      (MXALLCV2,2)  ,stat=ierr)
      ALLOCATE (RVA2      (MXCVFAC2,2)  ,stat=ierr)
      ALLOCATE (RHO2      (MXALLCVC,2)  ,stat=ierr)
      ALLOCATE (RMUT2     (MXALLCV2  )  ,stat=ierr)
      ALLOCATE (RMU2      (MXALLCV2  )  ,stat=ierr)
      ALLOCATE (YYS2      (MXALLCV2,MXCOMP,2),stat=ierr)
      ALLOCATE (RDS2      (MXALLCV2,MXCOMP)  ,stat=ierr)
      ALLOCATE (WDOT2     (MXALLCV2,MXCOMP)  ,stat=ierr)
      ALLOCATE (RMD2      (MXALLCV2)         ,stat=ierr)
      ALLOCATE (ccc2      (MXALLCV2)         ,stat=ierr)
      ALLOCATE (hhh2      (MXALLCV2,2)         ,stat=ierr)
      ALLOCATE (dvdd2     (MXCV2,3,2)        ,stat=ierr)
      ALLOCATE (UTAU2     (0:MXSSFBC2)       ,stat=ierr)
      ALLOCATE (VELBND2   (MXSSFBC2,3)       ,stat=ierr)
      ALLOCATE (TMPBND2   (MXSSFBC2)         ,stat=ierr)
      ALLOCATE (YYSBND2   (MXSSFBC2,MXCOMP)  ,stat=ierr)
      ALLOCATE (HTCBND2   (MXSSFBC2)         ,stat=ierr)
      ALLOCATE (MTCBND2   (MXSSFBC2)         ,stat=ierr)
      ALLOCATE (RADBND2   (MXSSFBC2)         ,stat=ierr)
!
      ALLOCATE (HEATL_C   (MXALLCV2)         ,stat=ierr)
      ALLOCATE (HEATG_C   (MXALLCV2)         ,stat=ierr)
      ALLOCATE (HLS       (MXALLCV2)         ,stat=ierr)
      ALLOCATE (HGS       (MXALLCV2)         ,stat=ierr)
      ALLOCATE (TMP_SAT   (MXALLCV2)         ,stat=ierr)
      ALLOCATE (FRIC_C    (MXALLCV2)         ,stat=ierr)
      ALLOCATE (MASS_I    (MXALLCV2)         ,stat=ierr)
!
      ALLOCATE (TMP      ( MXALLCV,2)        ,stat=ierr)
      ALLOCATE (YYS      ( MXALLCV,MXCOMP,2) ,stat=ierr)
      ALLOCATE (RDS      ( MXALLCV,MXCOMP)   ,stat=ierr)
      ALLOCATE (WDOT     ( MXALLCV,MXCOMP)   ,stat=ierr)
      ALLOCATE (RMD      ( MXALLCV)          ,stat=ierr)
      ALLOCATE (hhs_w    ( MXALLCV,MXcomp)   ,stat=ierr)
      ALLOCATE (cps_w    ( MXALLCV,MXcomp)   ,stat=ierr)
      ALLOCATE (CCC      ( MXALLCV)          ,stat=ierr)
      ALLOCATE (UTAU     ( 0:MXSSFBC)        ,stat=ierr)
!
      ALLOCATE (VELBND   ( MXSSFBC,3)        ,stat=ierr)
      ALLOCATE (PRSBND   ( MXSSFBC)          ,stat=ierr)
      ALLOCATE (TMPBND   ( MXSSFBC)          ,stat=ierr)
      ALLOCATE (YYSBND   ( MXSSFBC,MXCOMP)   ,stat=ierr)
      ALLOCATE (AKSBND   ( MXSSFBCR,MXRANS)  ,stat=ierr)
      ALLOCATE (HTCBND   ( MXSSFBC)          ,stat=ierr)
      ALLOCATE (HTFLUX   ( MXSSFBC)          ,stat=ierr)
      ALLOCATE (MTCBND   ( MXSSFBC)          ,stat=ierr)
      ALLOCATE (RADBND   ( MXSSFBC)          ,stat=ierr)
      ALLOCATE (VLASUM   ( MXCV)             ,stat=ierr)
      ALLOCATE (SUFBND   ( MXSSFBC_SUF,MXCOMPALL)  ,stat=ierr)
      ALLOCATE (SDOT     ( MXSSFBC_SUF,MXCOMPALL)  ,stat=ierr)
      ALLOCATE (MOLFRC   ( MXSSFBC_SUF,MXCOMPALL,2),stat=ierr)
      ALLOCATE (MOLCEN   ( MXSSFBC_SUF,MXCOMPALL)  ,stat=ierr)
      ALLOCATE (SITDEN   ( MXSSFBC_SUF,MXPHASE,2)  ,stat=ierr)
      ALLOCATE (BLKTHK   ( MXSSFBC_SUF,MXPHASE,2)  ,stat=ierr)
!
!      ALLOCATE (WDRBND   ( MXCV_WDR)      ,stat=ierr)
!
      ALLOCATE (msk      (MXALLCV)           ,stat=ierr)
      ALLOCATE (IAL      (MXALLCV,MAXIE)     ,stat=ierr)
      ALLOCATE (INL      (MXCV)              ,stat=ierr)
      ALLOCATE (INU      (MXCV)              ,stat=ierr)
      ALLOCATE (aalw     (MXALLCV,0:MAXIE)   ,stat=ierr)
      ALLOCATE (bb       (MXCV)              ,stat=ierr)
      ALLOCATE (W1K1     (MXCV)              ,stat=ierr)
      ALLOCATE (W1K2     (MXCV)              ,stat=ierr)
      ALLOCATE (W1K3     (MXCV)              ,stat=ierr)
      ALLOCATE (W1K4     (MXCV)              ,stat=ierr)
      ALLOCATE (W1K5     (MXCV)              ,stat=ierr)
      ALLOCATE (W1K6     (MXCV)              ,stat=ierr)
      ALLOCATE (W1K7     (MXALLCV)           ,stat=ierr)
      ALLOCATE (W1K8     (MXALLCV)           ,stat=ierr)
      ALLOCATE (W1K9     (MXALLCV)           ,stat=ierr)
      ALLOCATE (W1K10    (MXALLCV)           ,stat=ierr)
      ALLOCATE (W1K11    (MXALLCV)           ,stat=ierr)
!
      ALLOCATE (t0       (MXMAT)             ,stat=ierr)
      ALLOCATE (y0       (MXMAT,MXcomp)      ,stat=ierr)
      ALLOCATE (p0       (MXMAT)             ,stat=ierr)
      ALLOCATE (u0       (MXMAT),v0(MXMAT),w0(MXMAT),stat=ierr)
      ALLOCATE (vsgm0    (MXMAT,2),rho0(MXMAT)      ,stat=ierr)
      ALLOCATE (aks0     (MXMAT,MXrans)             ,stat=ierr)
!
      ALLOCATE (potn0    (MXMAT,MXpotn)             ,stat=ierr)
      ALLOCATE (potn_c   (MXMAT,MXpotn)             ,stat=ierr)
!
!      ALLOCATE (t02(MXMAT,NPHS),
!     &          y02(MXMAT,MXcomp,NPHS),stat=ierr)
!      ALLOCATE (u02(MXMAT,NPHS),
!     &          v02(MXMAT,NPHS),
!     &          w02(MXMAT,NPHS),stat=ierr)
!      ALLOCATE (rho02(MXMAT,NPHS),stat=ierr)
!
      ALLOCATE (t02(MXMAT),
     &          y02(MXMAT,MXcomp),stat=ierr)
      ALLOCATE (u02(MXMAT),
     &          v02(MXMAT),
     &          w02(MXMAT),stat=ierr)
      ALLOCATE (rho02(MXMAT)       ,stat=ierr)
!
      ALLOCATE (
     &          Ke_FC(MXMAT),
     &          Kc_FC(MXMAT),
     &          Pv_FC(MXMAT),
     &          Kap_FC(MXMAT),
     &          sigma_FC(MXMAT),contang_FC(MXMAT),
     &          surf_tenFC(MXMAT),
     &          elecCondFC(MXMAT),ionCondFC(MXMAT)
     &                       ,stat=ierr)
!
      ALLOCATE (icvprs   (MXMAT)       ,stat=ierr)
!
      ALLOCATE (fai_iMHD (MXMAT,2)                  ,stat=ierr)
      ALLOCATE (A_iMHD   (MXMAT,3,2)                ,stat=ierr)
      ALLOCATE (eps_iMHD (MXMAT)                    ,stat=ierr)
      ALLOCATE (eps0_iMHD(MXMAT)                    ,stat=ierr)
      ALLOCATE (mu_iMHD  (MXMAT)                    ,stat=ierr)
      ALLOCATE (mu0_iMHD (MXMAT)                    ,stat=ierr)
      ALLOCATE (sigma_iMHD(MXMAT)                   ,stat=ierr)
      ALLOCATE (sigma0_iMHD(MXMAT)                  ,stat=ierr)
      ALLOCATE (Current_iMHD(MXMAT,3,2)             ,stat=ierr)
      ALLOCATE (OMEGA_iMHD(MXMAT)                   ,stat=ierr)
!
      ALLOCATE (vctr     (MXCV_V,0:MXBND_V)         ,stat=ierr)
      ALLOCATE (vctrp    (MXCV_VP,0:MXBND_VP)       ,stat=ierr)
      ALLOCATE (FIELD_U    (MXCV_D,NFLID)         ,stat=ierr)
!
!
!
      ALLOCATE (NEIGHBOR (MEP,  MXALLCV_P),stat=ierr) 
      ALLOCATE (CC       (    MEP)        ,stat=ierr)
      ALLOCATE (HEIGHT0  (    MEP)        ,stat=ierr)
      ALLOCATE (HEIGHT1  (    MEP)        ,stat=ierr)
      ALLOCATE (VELN     (    MEP)        ,stat=ierr)
!
      ALLOCATE (XYZP     (    MXPRT,3)    ,stat=ierr)
      ALLOCATE (UVWP     (    MXPRT,3)    ,stat=ierr)
      ALLOCATE (PARCEL   (MXPRT)          ,stat=ierr)
      ALLOCATE (DIA      (MXPRT)          ,stat=ierr)
      ALLOCATE (PMASS    (MXPRT,2)        ,stat=ierr)
      ALLOCATE (HIS_P    (MXPRT)        ,stat=ierr)
      ALLOCATE (PYS      (MXPRT,MXCOMP,2) ,stat=ierr)
      ALLOCATE (HEIGHT   (MEP,MXPRT)      ,stat=ierr)
      ALLOCATE (FU       (MXPRT,2)      ,stat=ierr)
      ALLOCATE (FV       (MXPRT,2)      ,stat=ierr)
      ALLOCATE (FW       (MXPRT,2)      ,stat=ierr)
      ALLOCATE (RHOP     (MXPRT)          ,stat=ierr)
      ALLOCATE (TMPP     (MXPRT,2)        ,stat=ierr)
      ALLOCATE (INSIDE   (MXPRT)          ,stat=ierr)
      ALLOCATE (JOUT     (MXPRT)          ,stat=ierr)
!
      ALLOCATE (DYDTP    (MXPRT,MXCOMP)   ,stat=ierr)
      ALLOCATE (DHDTP    (MXPRT,2)          ,stat=ierr)
!
      ALLOCATE (DNORM    (MXCVFAC_P)      ,stat=ierr)
      ALLOCATE (MOMCELL  (MXALLCV_P,4)    ,stat=ierr)
      ALLOCATE (EVAP_Y   (MXALLCV_P,MXCOMP),stat=ierr)
      ALLOCATE (EVAP_H   (MXALLCV_P,2),stat=ierr)
      ALLOCATE (IPCELL   (MXALLCV_P)      ,stat=ierr)
      ALLOCATE (NEIBCELL (MXALLCV_P)      ,stat=ierr)
      
!
! --- dummy
!
      ALLOCATE (ys       (mxcomp)      ,stat=ierr)
      ALLOCATE (rcomp    (mxcomp)      ,stat=ierr)
      ALLOCATE (dum_c1   (mxcomp)      ,stat=ierr)
      ALLOCATE (mu_i     (mxcomp)      ,stat=ierr)
      ALLOCATE (k_i      (mxcomp)      ,stat=ierr)
      ALLOCATE (difu     (mxcomp)     ,stat=ierr)
      ALLOCATE (thmdif   (mxcomp)     ,stat=ierr)

      ALLOCATE (eva_comp (mxcomp)      ,stat=ierr)
      ALLOCATE (sinj_comp(mxcomp)      ,stat=ierr)
      ALLOCATE (aks_temp (mxrans)      ,stat=ierr)
      ALLOCATE (erryys   (0:mxcomp)    ,stat=ierr)
      ALLOCATE (ahk      (0:5,mxcomp)  ,stat=ierr)
      ALLOCATE (ahk2     (0:5,mxcomp)  ,stat=ierr)
      ALLOCATE (acpk_w   (0:5,mxcomp)  ,stat=ierr)
      ALLOCATE (errrans  (0:mxrans)    ,stat=ierr)
      ALLOCATE (dvdc     (0:MXMAT)     ,stat=ierr)
      ALLOCATE (rsdc     (0:MXMAT)     ,stat=ierr)
      ALLOCATE (dpdt     (MXMAT)       ,stat=ierr)
      ALLOCATE (matsld   (0:MXMAT,NG)     ,stat=ierr)
      ALLOCATE (vsm      (0:MXMAT)     ,stat=ierr)
      ALLOCATE (rsm      (0:MXMAT)     ,stat=ierr)
      ALLOCATE (rsm_prim (0:MXMAT)     ,stat=ierr)
      ALLOCATE (tsm      (0:MXMAT)     ,stat=ierr)
      ALLOCATE (QINL     (0:MXMAT),
     &          QOUT     (0:MXMAT),
     &          OFFSET   (0:MXMAT),
     &          SINL     (0:MXMAT),
     &          SOUT     (0:MXMAT)     ,stat=ierr)


      ALLOCATE (c1       (0:MXMAT),
     &          c2       (0:MXMAT),
     &          c3       (0:MXMAT),
     &          rmaxs    (0:MXMAT),
     &          rmaxs0   (0:MXMAT),
     6          alpha    (  MXMAT),
     &          beta     (  MXMAT),
     &          t        (0:MXMAT),
     &          s        (0:MXMAT),
     &          u        (0:MXMAT),
     &          v        (0:MXMAT),
     &          omega    (  MXMAT)     ,stat=ierr)
      ALLOCATE (ipx(MAXIE),aax(MAXIE)  ,stat=ierr)
      ALLOCATE (NCOL1(MAXIE),NCOL2(MAXIE)  ,stat=ierr)
      if(ierr/=0) call FFRABORT(1,'stop at allocating ')
!
!--------------------------------------------
! --- Creat list, Control Volume and metrix
!--------------------------------------------
!
      call allocat_array 
!
      call METR_FFLOW(
     &  LVEDGE,LVRTCV,LFUTAU,
     &  LBC_SSF,LCYCSF,FRSTCV,DISINL,
     &  SFAREA,SFCENT,CVVOLM,CVVOL0,CVCENT,wiface,DISALL,xta,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
!     &  LCV_CV,vctr,vctrp,NEIGHBOR,DNORM,FIELD_U,
     &  LCV_CV,vctrp,NEIGHBOR,DNORM,FIELD_U,
     &  Pstart,LCYCOLD,wifsld,OPPANG,locmsh,
     &  KDBY_W,GRDC_W,
     &  ierror)
!
      if(ierror.ne.0) goto 9999
!
      call mat_flag(0,Pstart,mat_cal)
!
      call main_FFLOW(
     &  LVEDGE,LVRTCV,LFUTAU,
     &  LBC_SSF,LCYCSF,FRSTCV,UTAU,DISINL,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,LCV_CV,
     &  Pstart,mat_cal,
     &  SFAREA,SFCENT,CVVOLM,CVVOL0,CVCENT,wiface,DISALL,xta,
     &  rva,rho,prs,pp0,tmp,yys,vel,velf,aks,rmut,rmu,rmd,rds,
     &  wdot,ccc,
     &  rmx_w,hhh,hhs_w,cps_w,cp_w,cr_w,diag_w,kdbt_w,kdby_w,
     &  rvx_w,rvd_W,grdc_w,
!     &  vlasum,dvdd,vctr,vctrp,FIELD_U,
     &  vlasum,dvdd,vctrp,FIELD_U,
     &  prsbnd,tmpbnd,yysbnd,velbnd,aksbnd,
     Y  htcbnd,mtcbnd,radbnd,HTFLUX,
     &  SUFBND,SDOT,MOLFRC,MOLCEN,SITDEN,BLKTHK,!WDRBND,
     &  LCYCOLD,wifsld,OPPANG,locmsh,
     &  ipfix,itery,iterk,iterv,iterp,repsy,repsk,repsv,repsp,
     &  aepsy,aepsk,aepsv,aepsp,erry,errk,errv,errp,
     &  errsdy,errsdk,errsdv,
!MHD
     &  reps_FAI,aeps_FAI,err_FAI,reps_A,
     &  aeps_A,err_A,iter_FAI,iter_A,errsdA,errsdFAI,
!Particle
     1  DNORM,NEIGHBOR,
     2  INSIDE,JOUT,PMASS,DIA,XYZP,UVWP,PARCEL,RHOP,PYS,TMPP,
     3  HEIGHT,FU,FV,FW,
     4  MOMCELL,EVAP_Y,EVAP_H,IPCELL,NEIBCELL,DYDTP,DHDTP,HIS_P,
! Euler two phase
     &  VEL2,PRS2,tmp2,RVA2,RHO2,RMUT2,RMU2,YYS2,RDS2,
     &  WDOT2,RMD2,ccc2,dvdd2,
     &  HEATL_C,HEATG_C,HLS,HGS,TMP_SAT,FRIC_C,MASS_I,UTAU2,hhh2,
     &  VELBND2,TMPBND2,YYSBND2,
     &  HTCBND2,MTCBND2,RADBND2,
     &  iptfix,POTNAL,POFLUX,POTNBC,
     &  iter_POTN,reps_POTN,aeps_POTN,err_POTN,
     &  itery2,iterv2,repsy2,repsv2,
     &  aepsy2,aepsv2,erry2,errv2,
     &  errsdy2,errsdv2,
!
     &  RadHeatFlux,P_rad,PAbsorb,PScatter,
     &  RadHeat,
     &  ierror)

      if(ierror.ne.0) goto 9999
!-------------
! --- End time
!-------------
      call fftime(999)
!
      CALL FFRABORT(0,'FRONTFLOW')
!
 9999 continue
!
      CALL FFRABORT(1,'Abnormaly Stop in FRONTFLOW')
!
      end program FRONTFLOW
