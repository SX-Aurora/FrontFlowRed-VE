!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_ini_E2P
     &    (ICV,IMAT_U,ncomp,nrans,NPHS,iters,times,x,y,z,waldis,vol,
     &     u2,v2,w2,p2,dens2,ys2,t2,alpha)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     User subroutine for inlet boundary condition
!     (For Euler two phase flow)
!
!=====================================================================
!
! --- [dummy arguments]
!
      integer,intent(in)    :: NPHS,ICV,IMAT_U
      integer,intent(in)    :: iters,ncomp,nrans
      real*8 ,intent(in)    :: times,x,y,z,waldis,vol
      real*8 ,intent(inout) :: u2,v2,w2,dens2,t2,p2 
      real*8 ,intent(inout) :: ys2(ncomp)
      real*8 ,intent(inout) :: alpha(NPHS)
!=====================================================================
! --- [local entities] 
!     
!      if(y>1.d0) then
!      if(z>0.15d0) then
      if(x<0.d0.and.y<2.d0) then
!
        alpha(1)=1.0d0
        alpha(2)=0.0d0
      else
        alpha(1)=0.0d0
        alpha(2)=1.0d0
      endif
      return
      end subroutine user_ini_E2P
