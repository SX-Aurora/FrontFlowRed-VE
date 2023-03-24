!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine writeENSIGHTHeader(fname_root,ifl,gridOnly,
     &                              animeStat,animeExt,animeStep,status)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use FVdata, only : timeFV
      use ENSIGHTdata, only : timeES
      use FFRdata, only : numVars,nameVars
      implicit none
!
      character*80,intent(in) :: fname_root
      integer,intent(in) :: ifl
      logical,intent(in) :: gridOnly
      logical,intent(in) :: animeStat
      character*32,intent(in) :: animeExt
      integer,intent(in) :: animeStep

      integer,intent(out) :: status
!      
      integer :: ios
      integer :: i,j,l
      integer :: today(8)
      character*80 :: fname
      character*80 :: char_80
      character*1024 :: char_1024
      integer :: sclnum,vctnum,clen,strptr
      integer,allocatable :: iscl(:),ivct(:)
      character*80,allocatable :: char_vct(:)

!     ******** grid file *********
! 
! --- open EnSight6 format file and write header
!     
      fname=trim(adjustl(fname_root)) // '.geo'
      if(animeStat .or. animeStep>0) then
         fname=trim(adjustl(fname)) // animeExt
      end if
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeENSIGHTHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
! --- Open EnSight file
!
!      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
      open(ifl,file=fname,form='unformatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeENSIGHTHeader)Cannot create ',trim(fname)
        status=1; return
      end if
!
! --- OUTPUT TO ENSIGHT GEOM FILE
!
!
! --- OUTPUTING HEADDERS
!
      char_80='Fortran Binary'
      write(ifl)char_80
      char_80='FrontFlow/Red ver.2 result to EnSight6 files'
      write(ifl)char_80
      char_80='converted from ffr2viz -rf ES option'
      write(ifl)char_80
      char_80='node id assign'
      write(ifl)char_80
      char_80='element id assign'
      write(ifl)char_80
!
      close(ifl)


!
!     ******** EnSight6 Case File *********
!
      if(animeStat) then
         if(animeStep==0) then  ! animation first step
            allocate(timeES(0:100000),stat=ios) ! MAX files
            if(ios/=0) then
               stop 'ERR: allocation in writeENSIGHTHeader'
            end if
            timeES=0.0
            timeES(animeStep)=timeFV
         elseif(animeStep>0) then
            timeES(animeStep)=timeFV
            return  ! *** note!
         end if
      end if

      fname=trim(adjustl(fname_root)) // '.case'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeENSIGHTHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
! --- Open EnSight Case file
!
      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
!      open(ifl,file=fname,form='unformatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeENSIGHTHeader)Cannot create ',trim(fname)
        status=1; return
      end if
!
! --- OUTPUT TO ENSIGHT CASE FILE
!
!
! --- OUTPUTING HEADDERS
!
      call date_and_time(VALUES=today)
      write(ifl,1001) '# '
      write(ifl,1002) '# ',today(1:3),' ',today(5:7)
      write(ifl,1001) '# FrontFlow/Red v2 to EnSight6'
      write(ifl,1003) '# Case File:',trim(fname)
      write(ifl,1001) ''
!
!
! --- OUTPUTING GEOMETRY NAMES
!
      write(ifl,1001) 'FORMAT'
      write(ifl,1001) ''
      write(ifl,1001) 'type:  ensight'
      write(ifl,1001) ''
      write(ifl,1001) 'GEOMETRY'
      write(ifl,1001) ''
      char_80=trim(adjustl(fname_root)) // '.geo'
      if(animeStat .or. animeStep>0) then
         char_80=trim(adjustl(char_80)) // '****'
         write(ifl,1003) 'model: 1 ',trim(char_80)
      else
         write(ifl,1003) 'model:   ',trim(char_80)
      end if
      
      write(ifl,1001) ''
!
! --- OUTPUTING VARIABLE NAMES
!
      if(.not. gridOnly) then
      allocate(iscl(1:numVars),ivct(1:numVars))
      allocate(char_vct(1:numVars))
      iscl=0
      ivct=0
      sclnum=0
      vctnum=0
      char_vct=''
      i=1
      do while(i<=numVars)
!     search string ';' as 'vector' in variable names.
      char_80=nameVars(i)
      clen = len_trim(adjustl(char_80))
      do j=1,clen
         strptr=scan(char_80(j:clen),';')
         if(strptr>=1.and.strptr<=clen) then
            vctnum=vctnum+1
            ivct(vctnum)=i
            char_vct(vctnum)=char_80(strptr+1:clen)
            i=i+3
            exit
         else
            sclnum=sclnum+1
            iscl(sclnum)=i
            i=i+1
            exit
         end if
      end do
      end do
      if(sclnum+vctnum*3/=numVars) then
         write(*,*) 'ERR: count vector names'
         stop
      end if


      write(ifl,1001) 'VARIABLE'
      write(ifl,1001) ''
      do l=1,sclnum
      i=iscl(l)
      char_1024=trim(adjustl(fname_root)) // '_' //
     &          trim(adjustl(nameVars(i))) // '.scl'
      if(animeStat .or. animeStep>0) then
         char_1024=trim(adjustl(char_1024)) // '****'
         write(ifl,1004) 'constant per case: 1 ','PI ','3.141592653'
         write(ifl,1004) 'scalar per node: 1 ',trim(nameVars(i)),'\t',
     &                   trim(char_1024)
      else
         write(ifl,1004) 'constant per case:   ','PI ','3.141592653'
         write(ifl,1004) 'scalar per node:   ',trim(nameVars(i)),'\t',
     &                   trim(char_1024)
      end if
      end do
      do l=1,vctnum
      i=ivct(l)
      char_1024=trim(adjustl(fname_root)) // '_' //
     &          trim(adjustl(char_vct(l))) // '.vct'
      if(animeStat .or. animeStep>0) then
         char_1024=trim(adjustl(char_1024)) // '****'
         write(ifl,1004) 'vector per node: 1 ',trim(char_vct(l)),'\t',
     &                   trim(char_1024)
      else
         write(ifl,1004) 'vector per node:   ',trim(char_vct(l)),'\t',
     &                   trim(char_1024)
      end if
      end do

      deallocate(iscl,ivct,char_vct)
      end if
!
!
! ---(OUTPUTING TIME ARRAYS)
!
      close(ifl)

      status=0

 1001 format(A)
 1002 format(A,I4,'/',I2,'/',I2,A,I2,':',I2,':',I2)
 1003 format(A,A)
 1004 format(A,A,A,A)
 1005 format(A,I10)
      end subroutine writeENSIGHTHeader

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine writeENSIGHTHeader_TIME(fname_root,ifl,
     &                                   animeStat,animeStep,status)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use ENSIGHTdata, only : timeES
      implicit none
!
      character*80,intent(in) :: fname_root
      integer,intent(in) :: ifl
      logical,intent(in) :: animeStat
      integer,intent(in) :: animeStep
      integer,intent(out) :: status
!      
      integer :: ios
      integer :: i,j,l
      character*80 :: fname

!
!     ******** EnSight6 Case File *********
!
      fname=trim(adjustl(fname_root)) // '.case'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeENSIGHTHeader_TIME)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
! --- Open EnSight Case file
!
      open(ifl,file=fname,form='formatted',action='write',
     &                    position='append',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeENSIGHTHeader_TIME)Cannot create ',trim(fname)
        status=1; return
      end if
!
! --- OUTPUT TO ENSIGHT CASE FILE
!
!
! --- OUTPUTING TIME ARRAYS
!
      if(.not.animeStat .and. animeStep>0) then  ! last animation loop
         write(ifl,1001) ''
         write(ifl,1001) 'TIME'
         write(ifl,1001) ''
         write(ifl,1005) 'time set: ',1
         write(ifl,1005) 'number of steps: ',animeStep
         write(ifl,1005) 'filename start number: ',0
         write(ifl,1005) 'filename increment: ',1
         write(ifl,*) 'time values: ',timeES(0:animeStep-1)
      end if
      close(ifl)

      status=0

 1001 format(A)
 1005 format(A,I10)
      end subroutine writeENSIGHTHeader_TIME
