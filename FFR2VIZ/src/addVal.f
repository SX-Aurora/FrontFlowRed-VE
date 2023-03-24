!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine addVal(name)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use FFRdata, only : numVars,nameVars
!
      implicit none
!      
      character*80,intent(in)  :: name
      character*80,allocatable :: nameVarTmp(:)
!      
      if(numVars>0) then
        allocate(nameVarTmp(1:numVars))
        nameVarTmp(1:numVars)=nameVars(1:numVars)
        deallocate(nameVars)
      end if
      numVars=numVars+1
      allocate(nameVars(1:numVars))
      nameVars=' '
      if(numVars>1) then
        nameVars(1:numVars-1)=nameVarTmp(1:numVars-1)
      end if
      nameVars(numVars)=name
!      
      deallocate(nameVarTmp)
      end subroutine addVal
