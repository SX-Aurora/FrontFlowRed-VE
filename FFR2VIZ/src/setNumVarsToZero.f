      subroutine setNumVarsToZero(MPprocNum,rtrnStat)
      use FFRdata, only : nvrtx,numVars,nameVars,totdata,node_id
!
      implicit none
!      
      integer,     intent(in)  :: MPprocNum
      integer,     intent(out) :: rtrnStat
!
      integer :: i,ierr
!      numVars=0
      numVars=1
      allocate(nameVars(1:numVars))
      allocate(totdata(1:numVars,1:nvrtx),stat=ierr)
      if(ierr.ne.0) then
         write(*,*) 'stop at allocating totdata(',
     &        numVars,nvrtx,') in setNumVarsToZero.'
         rtrnStat=1
         return
      end if

      if(MPprocNum>1) then
         nameVars(1)='CPU_ID'
         do i=1, nvrtx
            totdata(1,i)=node_id(i)
         end do
      else
         nameVars(1)='vertexID'
         do i=1, nvrtx
            totdata(1,i)=i
         end do
      end if

      rtrnStat=0
      return
      end subroutine setNumVarsToZero
