! *** *************************************************************
      subroutine append_resultfluent(fname_root,ifl,gridOnly,status)
! *** *************************************************************
      use FLUENTdata, only : F_NREGN,F_NKND,F_UVW,F_NFBOUN
      use FFRdata, only : nvrtx,ncell,nface,cord,NBFS,NFCE,
     &                    NBOUND,NFBOUN,SFBOUN,
     &                    numVars,numBVars,nameVars,
     &                    lvcell,lacell,kmesh,totdata,
     &                    lfcell,lvface,lcface,lbface,LVRTCV,LCVFAC
!
      implicit none
!
      character(80),intent(in) :: fname_root
      integer,intent(in) :: ifl
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
!      
      integer :: i,j,k,l,ios
      integer :: IC,IC1,IC2,IS,IV,numIV,nb
!
      real(8) :: cell_center_dbl(3),tmp_dbl
      integer :: isub,numsub,isub_size,izone,ivar(3)
       integer :: numZone,itype,face_cnt,buffle_flag
      integer,allocatable :: NSZONE(:),NEZONE(:)
      integer :: element_type,element_faces,face_nodes
      character(20),allocatable :: nameZone(:)
      character(15) :: NOTE_TYPE
      character(10) :: NOTE_SUBSECTION
      character(80) :: NOTE_ZONE
      character(80) :: fname
      character(10) :: tmpchar10
      integer,parameter :: ssize=256,linew=80
      character(len=ssize) :: sln
!
!
      end subroutine append_resultfluent
!
