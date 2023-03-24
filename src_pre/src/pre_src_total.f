!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine src_rho(NALLCV,rsrc)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_source,only : nscnd,kdscnd,kdmass,
     & icscnd,lcscnd,wdscnd
!
      implicit none
!
! 1. Set source term for density
!
! --- [dummy arguments]
!
      integer,intent(in)  :: NALLCV
      real*8 ,intent(out) :: rsrc(NALLCV)
!
! --- [local entities]
!
      real*8  :: rhox
      integer :: l,m
!
!
!
      rsrc=0.d0
      if( nscnd.lt.1 ) return
!
      do 1000 m=1,nscnd
      if( kdscnd(m).eq.kdmass ) then
        rhox=wdscnd(1,m)
        do 100 l=icscnd(m-1)+1,icscnd(m)
        rsrc(lcscnd(l))=rhox
  100   continue
      endif
 1000 continue
!
      end subroutine src_rho
!
