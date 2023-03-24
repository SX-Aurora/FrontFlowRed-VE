!
!      module module_include
!      module module_hpcutil
!      module module_hpc_input
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      module module_include
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- array dimension
!
!      include 'dimen.h_serial'
!
!      end module module_include
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_hpc_input
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_hpc_input)'
      integer             :: NPART,UCDFLAG,WALLFLAG
      character (len=80)  :: GRIDFIL,BCFIL,INPFIL,UCDFIL
      character (len=80)  :: WALFIL,DIRHED
      character (len=80)  :: HEADERg, HEADERb, HEADERc, HEADERw
      character (len=80)  :: HEADERs,HEADERi,HEADERr,HEADERt
      character (len=80)  :: HEADER_Pr,HEADER_Pi
      character (len=80)  :: INITin_P,RESTRTout_P
      character (len=80)  :: HEADERa,HEADERm
      character (len=80)  :: GRIDin,BCin,COMMin
      character (len=80)  :: SRCin,INITin,RESLTout,RESTRTout,RESLTgrid
      character (len=80)  :: FORCEout,GEOMin,WALLin
      character (len=80)  :: statis,animfile,sld_wk,suf_wk
      character (len=80) :: WLDISout
      character (len=80) :: dflxTmpFname1,dflxTmpFname2
      character (len=80)  :: HEADER_pro,PROBEout     !!!modify
!
!
!/////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata(my_rank,NPE,ifli,ifll,ifle,
     &                     cntlnam,Dmns,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy qrgument]
!
      integer,intent(in)  :: ifli,ifll,ifle,my_rank,NPE
      character,intent(in) :: cntlnam
      integer,intent(out) :: ierror
      character (len=80),intent(out)  :: Dmns
!
! --- [local entities]
!
      integer :: ierr1=0
      character(80),parameter :: blank=' '
!
! --- [namelist]
!
      NPART=1
      UCDFLAG=0
      WALLFLAG=0
      GRIDFIL=blank
      BCFIL=blank
      INPFIL=blank
      UCDFIL=blank
      WALFIL=blank
      DIRHED=blank
      HEADERg=blank
      HEADERb=blank
      HEADERc=blank
      HEADERw=blank
      GRIDin=blank
      BCin=blank
      COMMin=blank
      Dmns=blank
      HEADERa=blank
!
      return
 9999 continue
      ierror=1
      write(ifle,*) modnam,subnam
      end subroutine inputdata
!
      end module module_hpc_input
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_hpcutil
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This module is dummy hpc module
!
        use module_hpc_input,
     &                 only : NPART,UCDFLAG,WALLFLAG,
     &                        GRIDFIL,BCFIL,INPFIL,UCDFIL,
     &                        WALFIL,DIRHED,
     &                        HEADERg, HEADERb, HEADERc, HEADERw,
     &                        HEADERs,HEADERi,HEADERr,HEADERt,
     &                        HEADERa,HEADERm,
     &                        HEADER_Pr,HEADER_Pi,
     &                        INITin_P,RESTRTout_P,
     &                        GRIDin,BCin,COMMin,
     &                        SRCin,INITin,RESLTout,RESTRTout,RESLTgrid,
     &                        FORCEout,GEOMin,WALLin,
     &                        statis,animfile,sld_wk,suf_wk
     &                        ,WLDISout
     &                        ,dflxTmpFname1,dflxTmpFname2
     &                        ,HEADER_pro,PROBEout !!!modify
!
      integer,public :: PETOT
      integer,public :: my_rank
      integer,public :: SOLVER_COMM
      integer,public :: errno
      integer,public :: NEIBPETOT=1
      integer,public,allocatable:: NEIBPE(:)
      integer,public :: NPE,ROOT
!
      integer ,allocatable :: NOD_IDg(:)
      integer ,allocatable :: IMPORT_index(:), IMPORT_item(:),
     &                        EXPORT_index(:), EXPORT_item(:)
      integer ,allocatable :: IMPORT_indexIP(:),
     &                        EXPORT_indexIP(:)
      integer ,allocatable :: IBC_P_NO(:,:,:,:)
      integer ,allocatable :: IBC_P_TOL(:,:,:)
!
      end module module_hpcutil
!
