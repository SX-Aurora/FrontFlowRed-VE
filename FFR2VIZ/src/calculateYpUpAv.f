      subroutine calculateYpUpAv(ifl,rtrnStat)
!      use FVdata
      use FFRdata, only  : nvrtx,totdata,wallDist
      use ProgData, only : ypval,upval,ypAvCnt,divnumYPAV,
     &                     ypAvupfn,ReTau,uTau,scaleFactor
!
      implicit none
!      
      integer,intent(in)  :: ifl
      integer,intent(out) :: rtrnStat
      integer :: i,j,k,ios
      integer,parameter :: lowLim=10
      !upVal :: 1:y+ 2:ui 3:vi 4:wi 5:um 6:vm 7:wm 8:uu 9:vv 10:ww
      !vars  ::      3:ui 4:vi 5:wi 6:um 7:vm 8:wm 9:uu 10:vv 11:ww
      allocate(upVal(0:divnumYPAV,1:10))
      allocate(ypAvCnt(1:divnumYPAV))
      
      upVal=0.d0; ypAvCnt=0

      do i=1,nvrtx
         do j=1,divnumYPAV
            if((ypVal(j-1)-wallDist(i))*
     &           (ypVal(j)-wallDist(i))<=0.d0) then
               do k=2,10
                  upVal(j,k)=upVal(j,k)+totdata(k+1,i)
               end do
               upVal(j,1)=upVal(j,1)+wallDist(i)
               ypAvCnt(j)=ypAvCnt(j)+1
               exit
            end if
         end do
      end do
      
      do i=1,divnumYPAV
         if(ypAvCnt(i)>0) then
            upVal(i,1:10)=upVal(i,1:10)/real(ypAvCnt(i))
         end if
      end do
      
      upVal(:,1)=upVal(:,1)*ReTau*scaleFactor
      upVal(:,2:7)=upVal(:,2:7)/uTau
      
      open(ifl,file=ypAvupfn,form='formatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'Cannot open file',trim(ypAvupfn)
        rtrnStat=1;     return
      end if
 1000 format(10E15.6)
      write(ifl,1000)upVal(0,:)
      do i=1,divnumYPAV
        if(ypAvCnt(i)>=lowLim) then
        write(ifl,1000)upVal(i,:)
        end if
      end do
      close(ifl)
      rtrnStat=0
      end subroutine calculateYpUpAv
