! **********************************************************************
      subroutine append_resultensight(fname_root,ifl,CPUnum,gridOnly,
     &                                animeStat,animeExt,animeStep,
     &                                status)
! **********************************************************************
      use FFRdata, only : nvrtx,ncell,cord,NBFS,NBOUND,NFBOUN,SFBOUN,
     &                    numVars,nameVars,numBVars,lvcell,
     &                    IBFACE,IFFACE,kmesh,totdata
!
      implicit none
!
      character*80,intent(in) :: fname_root
      integer,intent(in) :: ifl,CPUnum
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
      logical,intent(in) :: animeStat
      character*32,intent(in) :: animeExt
      integer,intent(in) :: animeStep
!
      integer :: grd,iMAT,iv
      integer :: i,j,k,l,ios,numSfBRG,nnum
      integer :: type,tempElems,tempFaces
      character*80 :: elem_type,face_type
      integer,allocatable :: ip(:)
      integer,allocatable :: ib(:),it(:)
      character*80 :: fname
      integer :: int_4
      character*80 :: char_80
      integer :: sclnum,vctnum,clen,strptr
      integer,allocatable :: iscl(:),ivct(:)
      character*80,allocatable :: char_vct(:)
!
! --- 
      end subroutine append_resultensight
