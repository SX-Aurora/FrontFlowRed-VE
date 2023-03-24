      subroutine writeAVSHeader(fname,ifl,status)

!
      use AVSheaderTag, only : AVS_NAME,AVS_STEP,AVS_CYCL,AVS_STPN
      use FFRdata, only : nvrtx,ncell
      implicit none

      character*80,intent(in) :: fname
      integer,intent(in) :: ifl
      integer,intent(out) :: status
      
      integer :: ios
      integer :: i

!     open AVS format file and write header
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeAVSHeader)File ',trim(fname),
     &            ' is already opened as number ',ios
        status=1; return
      end if
      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeAVSHeader)Cannot create ',trim(fname)
        status=1; return
      end if
      write(ifl,115) trim(AVS_NAME)
      write(ifl,115) trim(AVS_STEP)
      write(ifl,115) trim(AVS_CYCL)
      write(ifl,115) trim(AVS_STPN)
      write(ifl,'(2i8)') nvrtx, ncell
  115 format(a,' ')
      
      close(ifl)
      status=0
      end subroutine writeAVSHeader
