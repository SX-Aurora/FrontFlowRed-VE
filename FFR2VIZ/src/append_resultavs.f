! *** *************************************************************
      subroutine append_resultavs(fname,ifl,CPUnum,procTag,
     &                    gridOnly,gridFormat,status)
! *** *************************************************************
      use FFRdata, only : nvrtx,ncell,ncompFFR,numVars,nameVars,
     &                    lvcell,cord,kmesh,totdata
      implicit none
!
      character*80,intent(in) :: fname
      character*2, intent(in) :: gridFormat
      integer,intent(in) :: ifl,CPUnum
      logical,intent(in) :: procTag(0:CPUnum-1)
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
!
      integer :: i,j,k,l,ios,numSfBRG,nnum
      real*8 :: tmpdbl
! DEBUG
      character*80 :: formatKey,formatKeytmp
!
      nnum = 0
!
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(append_resultavs)File ',trim(fname),
     &            ' has been opened as number ',ios
        status=1; return
      end if
      open(ifl,file=fname,form='formatted',action='write',
     &                              position='append',iostat=ios)
      if(ios/=0) then
        write(*,*)'(append_resultavs) Cannot open ',trim(fname)
        status=1; return
      end if
!      
      do nnum = 1 , nvrtx
         write(ifl,fmt='(1i8,3f13.7)') nnum,cord(1,nnum),
     &           cord(2,nnum),cord(3,nnum)
      end do
!
! --- --------------------------------------------------------------
       do nnum = 1 , ncell
         select case(kmesh(nnum))
         case(1)
            write(ifl,fmt='(1i8,1x,i1,a4,4i8)') nnum,0,' tet'
     &           ,lvcell(3,nnum),lvcell(4,nnum)
     &           ,lvcell(1,nnum),lvcell(2,nnum)

         case(2)
            write(ifl,fmt='(1i8,1x,i1,a4,8i8)') nnum,0,' hex'
     &           ,lvcell(5,nnum),lvcell(6,nnum)
     &           ,lvcell(7,nnum),lvcell(8,nnum)
     &           ,lvcell(1,nnum),lvcell(2,nnum)
     &           ,lvcell(3,nnum),lvcell(4,nnum)

         case(3)
            write(ifl,fmt='(1i8,1x,i1,a6,6i8)') nnum,0,' prism'
     &           ,lvcell(4,nnum),lvcell(5,nnum)
     &           ,lvcell(6,nnum),lvcell(1,nnum)
     &           ,lvcell(2,nnum),lvcell(3,nnum)

         case(4)
            write(ifl,fmt='(1i8,1x,i1,a4,5i8)') nnum,0,' pyr'
     &           ,lvcell(5,nnum),lvcell(1,nnum)
     &           ,lvcell(2,nnum),lvcell(3,nnum)
     &           ,lvcell(4,nnum)

         case default
            stop 'Incomplete element description'
         end select

       end do
!
! --- --------------------------------------------------------------
      if (.not.(gridOnly)) then
!        write out variables, added a 'velovity magnitude'
! DEBUG
         write(ifl,'(i8,1x,i2)') numVars+1,0
         write(formatKeytmp,*) numVars+1
         write(formatKey,*)'(i8,1x,',
     &        trim(adjustl(formatKeytmp)),'(i8,1x))'
         write(ifl,formatKey) numVars+1,(1,i=1,numVars+1)
c         write(ifl,*) numVars+1,' ',0
c         write(ifl,*) numVars+1,' ',(1,' ',i=1,numVars+1)
!         write(ifl,*) numVars,' ',0
!         write(ifl,*) numVars,' ',(1,' ',i=1,numVars)
!
         write(ifl,*) trim(nameVars(1))//', Pa'
         write(ifl,*) trim(nameVars(2))//', Pa'
         write(ifl,*) trim(nameVars(3))//', kg/m^3'
         write(ifl,*) trim(nameVars(4))//', K'
         write(ifl,*) trim(nameVars(5))//', N'
         write(ifl,*) trim(nameVars(6))//', m/s'
         write(ifl,*) trim(nameVars(7))//', m/s'
         write(ifl,*) trim(nameVars(8))//', m/s'
!
         do i = 9 , numVars
            write(ifl,*) trim(nameVars(i))//', Non-D'
         end do
         write(ifl,*) 'Mag_vel'//', m/s'
!
         do nnum = 1 , nvrtx
            tmpdbl=sqrt(totdata(6,nnum)*totdata(6,nnum)+
     &                  totdata(7,nnum)*totdata(7,nnum)+
     &                  totdata(8,nnum)*totdata(8,nnum))
! DEBUG
         write(formatKeytmp,*) numVars
         write(formatKey,*)'(i8,1x,',
     &        trim(adjustl(formatKeytmp)),'(E16.9,1x),E16.9)'
         write(ifl,formatKey) nnum,
     &        (totdata(i,nnum),i=1,numVars),tmpdbl
c            write(ifl,*) nnum,' ',(totdata(i,nnum),' ',i=1,numVars),
c     &                   tmpdbl
!            write(ifl,*) nnum,' ',(totdata(i,nnum),' ',i=1,numVars)
         end do

      else  ! if (.not.(gridOnly)) then
         if(numVars<=0) then
            write(*,*) 'Warning:there is no dummy vars for AVS.'
            status=1
            return
         end if
         write(ifl,fmt='(a2,2x,a1)') '1','0'
         write(ifl,fmt='(10(a2,1x))') '1','1'
         write(ifl,*) trim(nameVars(1))
         do nnum = 1 , nvrtx
            write(ifl,fmt='(1i8,f14.5)') nnum,totdata(1,nnum)
         end do
      end if
!
! --- --------------------------------------------------------------
      close(ifl)
      status=0
      end subroutine append_resultavs
