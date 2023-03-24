!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine writeFVHeader(fname,ifl,status)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use FVheaderTag, only : FV_MAGIC,FV_NAME
      use FVdata, only : sfRsltFlg,sfClkFlg,
     &                   timeFV,fsmachFV,alphaFV,reFV
      use FFRdata, only : NBOUND,SFBOUN,
     &                    numVars,nameVars,numBVars,nameBVars
      implicit none
!
      character*80,intent(in) :: fname
      integer,intent(in) :: ifl
      integer,intent(out) :: status
      logical :: opn=.false.
!      
      integer :: ios
      integer :: i
! 
! --- open FieldView format file and write header
!     
!
      inquire(file=fname,OPENED=opn)
!      inquire(file=fname,number=ios)
!      if(ios>=0) then
      if(opn) then
        write(*,*)'(writeFVHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
! --- Open FV result file
!
!      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
      open(ifl,file=fname,form='unformatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeFVHeader)Cannot create ',trim(fname)
        status=1; return
      end if
!
! --- OUTPUT TO FIELDVIEW DATA FILE
!
!
! --- OUTPUTING HEADDERS
!
      write(ifl)FV_MAGIC
      write(ifl)FV_NAME
      write(ifl)2,5
      write(ifl)timeFV,fsmachFV,alphaFV,reFV
      write(ifl)1
      write(ifl)NBOUND
!
! --- OUTPUTING BOUNDARY-FACE NAMES
!
      do 3010 i=1,NBOUND
        write(ifl)sfRsltFlg(i),sfClkFlg(i),SFBOUN(i)
        write(*,*) ' #### Boundary number and name :',i
     &               ,trim(SFBOUN(i))
 3010 continue
!
! --- OUTPUTING VARIABLE NAMES
!
      write(ifl) numVars
      do 3020 i=1,numVars
        write(ifl) nameVars(i)
        write(*,*) ' #### Variables number and name :',i
     &               ,trim(nameVars(i))
 3020 continue
!
! --- 
      write(ifl) numBVars
      do 3030 i=1,numBVars
        write(ifl)nameBVars(i)
 3030 continue
      
      close(ifl)
      status=0
      end subroutine writeFVHeader

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine writeFVHeader_ASCII(fname,ifl,status)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use FVheaderTag, only : FV_MAGIC,FV_NAME
      use FVdata, only : sfRsltFlg,sfClkFlg,
     &                   timeFV,fsmachFV,alphaFV,reFV
      use FFRdata, only : NBOUND,SFBOUN,
     &                    numVars,nameVars,numBVars,nameBVars
      implicit none
!
      character*80,intent(in) :: fname
      integer,intent(in) :: ifl
      integer,intent(out) :: status
!      
      integer :: ios
      integer :: i
! 
! --- open FieldView format file and write header
!     
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeFVHeader_ASCII)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
! --- Open FV result file
!
      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
!      open(ifl,file=fname,form='unformatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeFVHeader_ASCII)Cannot create ',trim(fname)
        status=1; return
      end if
!
! --- OUTPUT TO FIELDVIEW DATA FILE
!
!
! --- OUTPUTING HEADDERS
!
C      write(ifl,*)FV_MAGIC
      write(ifl,'(A,X,I2,X,I2)')trim(FV_NAME),2,4
!      write(ifl,'(A)'      )trim(FV_NAME)
!      write(ifl,'(I2,X,I2)')2,5
      write(ifl,'(A)'      )'CONSTANTS'
      write(ifl,*)timeFV
      write(ifl,*)fsmachFV
      write(ifl,*)alphaFV
      write(ifl,*)reFV
!      write(ifl,*          )timeFV,fsmachFV,alphaFV,reFV
      write(ifl,'(A)'      )'GRIDS'
      write(ifl,'(I5)'     )1
!
! --- OUTPUTING BOUNDARY-FACE NAMES
!
      write(ifl,'(A)'      )'Boundary Table'
      write(ifl,'(I8)'     )NBOUND
      do 3010 i=1,NBOUND
        write(ifl,'(I1,X,A)') 0, trim(SFBOUN(i))
!        write(ifl,'(3(I2,X),A)') 0, sfRsltFlg(i),sfClkFlg(i),
!     &                           trim(SFBOUN(i))
        write(*,*) ' #### Boundary number and name :',i
     &               ,trim(SFBOUN(i))
 3010 continue
!
! --- OUTPUTING VARIABLE NAMES
!
      write(ifl,'(A)'      )'Variable Names'
      write(ifl,'(I8)'     ) numVars
      do 3020 i=1,numVars
        write(ifl,'(A)') trim(nameVars(i))
        write(*,*) ' #### Variables number and name :',i
     &               ,trim(nameVars(i))
 3020 continue
!
! --- 
!      write(ifl,'(A)'      )'Boundary Variable Names'
!      write(ifl,'(I8)'     ) numBVars
!      do 3030 i=1,numBVars
!        write(ifl,'(A)')trim(nameBVars(i))
! 3030 continue
      
      close(ifl)
      status=0
      end subroutine writeFVHeader_ASCII

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine writeFVHeader_SPLIT(fname_root,ifl,status)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use FVheaderTag, only : FV_MAGIC,FV_NAME
      use FVdata, only : sfRsltFlg,sfClkFlg,
     &                   timeFV,fsmachFV,alphaFV,reFV
      use FFRdata, only : NBOUND,SFBOUN,
     &                    numVars,nameVars,numBVars,nameBVars
      implicit none
!
      character*80,intent(in) :: fname_root
      integer,intent(in) :: ifl
      integer,intent(out) :: status
!      
      integer :: ios
      integer :: i
      character*80 :: fname
!
!     ******** grid file *********
! 
! --- open FieldView format file and write header
!     
      fname=trim(adjustl(fname_root)) // '.gfv'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeFVHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
! --- Open FV file
!
!      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
      open(ifl,file=fname,form='unformatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeFVHeader)Cannot create ',trim(fname)
        status=1; return
      end if
!
! --- OUTPUT TO FIELDVIEW DATA FILE
!
!
! --- OUTPUTING HEADDERS
!
      write(ifl)FV_MAGIC
      write(ifl)FV_NAME
      write(ifl)2,7,1,0 ! FV_GRIDS_FILE
      write(ifl)1
      write(ifl)NBOUND
!
! --- OUTPUTING BOUNDARY-FACE NAMES
!
      do 3010 i=1,NBOUND
        write(ifl)sfRsltFlg(i),sfClkFlg(i),SFBOUN(i)
        write(*,*) ' #### Boundary number and name :',i
     &               ,trim(SFBOUN(i))
 3010 continue
!
      close(ifl)
!
!     ******** result file *********
!
      fname=trim(adjustl(fname_root)) // '.rfv'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeFVHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
! --- Open FV file
!
!      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
      open(ifl,file=fname,form='unformatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeFVHeader)Cannot create ',trim(fname)
        status=1; return
      end if
!
! --- OUTPUT TO FIELDVIEW DATA FILE
!
!
! --- OUTPUTING HEADDERS
!
      write(ifl)FV_MAGIC
      write(ifl)FV_NAME
      write(ifl)2,7,2,0 ! FV_RESULTS_FILE
      write(ifl)timeFV,fsmachFV,alphaFV,reFV
      write(ifl)1
!
! --- OUTPUTING VARIABLE NAMES
!
      write(ifl) numVars
      do 3020 i=1,numVars
        write(ifl) nameVars(i)
        write(*,*) ' #### Variables number and name :',i
     &               ,trim(nameVars(i))
 3020 continue
!
! --- 
      write(ifl) numBVars
      do 3030 i=1,numBVars
        write(ifl)nameBVars(i)
 3030 continue
      
      close(ifl)

      status=0
      end subroutine writeFVHeader_SPLIT

