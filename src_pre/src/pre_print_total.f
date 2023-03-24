!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine printi(ctitl,dbg,ns,nko)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_nowtime,only : iter,time
      implicit none
      character(*),intent(in) :: ctitl
      integer     ,intent(in) :: ns,nko
      integer     ,intent(in) :: dbg(ns,nko)
!
      integer  i,n,i1,i2
!
      write(*,6001) ctitl,iter,time
      do 100 n=1,nko
      i2=0
  101 continue
      i1=i2+1
      i2=min(i2+10,ns)
      if( i1.gt.1 ) then
        write(*,'(10x,a,10i10)') ':',(dbg(i,n),i=i1,i2)
      else
        write(*,'(i10,a,10i10)') n,':',(dbg(i,n),i=i1,i2)
      endif
      if( i2.lt.ns ) goto 101
 100  continue
      write(*,*)
!
      return
 6001 format('array ',a,' : iter =',i10,', time =',1pd12.5)
      end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine printj(ctitl,dbg,ns,nko)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_nowtime,only : iter,time
      implicit none
      character(*),intent(in) :: ctitl
      integer     ,intent(in) :: ns,nko
      integer     ,intent(in) :: dbg(ns,nko)
!
      integer  i,j,n,i1,i2
!
      do 100 n=1,ns
      write(*,6001) ctitl,iter,time
      write(*,6002) n
      write(*,'(5x,10(i10,1x))') (i,i=1,10)
      write(*,'(a,10(a))') '-----',('-----------',i=1,10)
      i2=0
      j=0
  103 continue
      i1=i2+1
      i2=min(i2+10,nko)
      write(*,'(i4,a,10(i10,1x))') j,':',(dbg(n,i),i=i1,i2)
      j=j+10
      if( i2.lt.nko ) goto 103
  100 continue
      write(*,*)
!
      return
 6001 format('array ',a,' : iter =',i10,', time =',1pd12.5)
 6002 format('N = ',i4)
      end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine printw(ctitl,dbg,ns,nko)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_nowtime,only : iter,time
      implicit none
      character(*),intent(in) :: ctitl
      integer     ,intent(in) :: ns,nko
      real*8      ,intent(in) :: dbg(ns,nko)
!
      integer  i,n,i1,i2
!
      write(*,6001) ctitl,iter,time
      do 100 n=1,nko
      i2=0
  101 continue
      i1=i2+1
      i2=min(i2+10,ns)
      if( i1.gt.1 ) then
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
!        write(*,'(10x,a,10(1pe12.5,x))') ':',(dbg(i,n),i=i1,i2)
        write(*,'(10x,a,10(1pe12.5,1x))') ':',(dbg(i,n),i=i1,i2)
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
      else
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
!        write(*,'(i10,a,10(1pe12.5,x))') n,':',(dbg(i,n),i=i1,i2)
        write(*,'(i10,a,10(1pe12.5,1x))') n,':',(dbg(i,n),i=i1,i2)
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
      endif
      if( i2.lt.ns ) goto 101
 100  continue
      write(*,*)
!
      return
 6001 format('array ',a,' : iter =',i10,', time =',1pd12.5)
      end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine printx(ctitl,dbg,ns,nko)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_nowtime,only : iter,time
      implicit none
      character(*),intent(in) :: ctitl
      integer     ,intent(in) :: ns,nko
      real*8      ,intent(in) :: dbg(ns,nko)
!
      integer  i,j,n,i1,i2
!
      do 100 n=1,ns
      write(*,6001) ctitl,iter,time
      write(*,6002) n
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
!      write(*,'(10(i12,x))') (i,i=1,10)
      write(*,'(10(i12,1x))') (i,i=1,10)
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
      write(*,'(a,10(a))') '-----',('-------------',i=1,10)
      i2=0
      j=0
  103 continue
      i1=i2+1
      i2=min(i2+10,nko)
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
!      write(*,'(i4,a,10(1pe12.5,x))') j,':',(dbg(n,i),i=i1,i2)
      write(*,'(i4,a,10(1pe12.5,1x))') j,':',(dbg(n,i),i=i1,i2)
!++++++Modified by T.Unemura 040520+++++++++++++++++++++++++
      j=j+10
      if( i2.lt.nko ) goto 103
  100 continue
      write(*,*)
!
      return
 6001 format('array ',a,' : iter =',i10,', time =',1pd12.5)
 6002 format('N = ',i4)
      end
