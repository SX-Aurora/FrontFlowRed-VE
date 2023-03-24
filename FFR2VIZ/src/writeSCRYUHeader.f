      subroutine writeSCRYUHeader(fname,ifl,status)
!
      use SCRYUheaderTag, only : HEADER,SOLV,COMENT,IVSN,ITYP,
     &                           L2D3D,LFORT
      use SCRYUdata, only : IDAT,NCYC,TIME,AMNUM
      use FFRdata,   only : iterFFR,timeFFR
!
      implicit none
!
      character*80,intent(in) :: fname
      integer,intent(in) :: ifl
      integer,intent(out) :: status
      
      character*8 :: cdate
      integer :: i,ios
!
! --- Today's date ---
      call date_and_time(cdate)
      read(cdate,*) IDAT
      write(*,*) ' #### Todays date is',IDAT
!
! --- NCYC: cycle, TIME: time
!
      NCYC  = iterFFR
      TIME  = timeFFR
!
      AMNUM = 1.0e+20    !(kessoku-chi?)
!
!     open SCRYU format file and write header
! --- -----------------------------------------------------
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeSCRYUHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
      open(ifl,file=fname,form='unformatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeSCRYUHeader)Cannot create ',trim(fname)
        status=1; return
      end if
! --- -----------------------------------------------------
!
      write(ifl) HEADER
      write(ifl) IVSN
      write(ifl) SOLV
      write(ifl) ITYP,L2D3D,LFORT
      write(ifl) IDAT
      write(ifl) COMENT
      write(ifl) NCYC,TIME
      write(ifl) AMNUM
!
      close(ifl)
      status=0
!
      return
      end subroutine writeSCRYUHeader
