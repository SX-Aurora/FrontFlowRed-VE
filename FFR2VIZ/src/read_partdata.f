! *** ******************************************************************
      subroutine read_partdata(MPprocNum,rtrnStat)
! *** ******************************************************************
      use FFRdata, only : nvrtx,node_id
!
      implicit none
!     
      integer,     intent(in)  :: MPprocNum
      integer,     intent(out) :: rtrnStat
!
      integer :: ifl, CPU
      integer :: i,j,ierr
!
      character*32  :: hmcpu
      character*80  :: PartfileName
      character*255 :: filePath
      character*255 :: hpcPath    ! extern function
!
!     *** re-writed for all of following 20040604 -- by onishi
!     read 'CONT.HPC' and 'test.m.part.x'
!       for deciding vertex-id(node_id) on HPC
      if(MPprocNum>1) then
         allocate(node_id(nvrtx),stat=ierr)
         if(ierr.ne.0) then
            stop 'stop at allocating in readAndStoreRslFFR.nid'
            rtrnStat=1
         end if
         call read_fort1_hpc(PartfileName)

!        creating metis part file name 'test.m.part.x',
!          'hmcpu', a character of total CPU num.
         write(hmcpu,*) MPprocNum
         PartfileName = adjustl(PartfileName)
         filePath = trim(PartfileName)//'.part'//
     &        '.'//trim(adjustl(hmcpu))
         write(*,*) ' #### Reading part file:',trim(filePath)
!
         ifl=30
         open(ifl,file=filePath,iostat=ierr)
         if(ierr/=0) then
            write(*,*)'Cannot open file :',filePath
            rtrnStat=1
            return
         end if
         do i = 1 , nvrtx
            read(ifl,*) CPU
            node_id(i) = CPU + 1
         end do
         close(ifl)
      end if

      rtrnStat=0
      return
      end subroutine read_partdata
