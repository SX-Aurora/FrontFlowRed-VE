
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_ini
     &          (ICV,IMAT_U,ncomp,nrans,npotn,
     &           iters,times,x,y,z,waldis,vol,
     &           u,v,w,p,dens,ys,t,aks,potn)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     User subroutine for inlet boundary condition
!
!=====================================================================
!
      use module_io,  only     : ifle,ifll,gdScale
      use module_movegrid,only : ipiston,omeini,mv_dir
      use module_model,only    : encomp
      
! --- [dummy arguments]
!
      IMPLICIT NONE
!
      integer,intent(in)    :: ICV,IMAT_U
      integer,intent(in)    :: iters,ncomp,nrans,npotn
      real*8 ,intent(in)    :: times,x,y,z,waldis,vol
      real*8 ,intent(inout) :: u,v,w,p,dens,t
      real*8 ,intent(inout) :: ys(ncomp),aks(nrans),potn(npotn)
!
!=====================================================================
! --- [local entities] 
      real(8),parameter :: turbMagFactor=7.d0
      integer :: mainDrivID=1
      real(8),parameter :: mainDrivRe=1102.210813d0
      real(8),parameter :: mainDrivUmax=0.23d0
      real(8),parameter :: mainDrivRadi=0.075d0
      real(8),parameter :: mainDrivTemp=321.15d0
      character*1,parameter :: mainDrivDirecTag='z'

      integer :: subDrivID=2
      real(8),parameter :: subDrivRe=1525.006563d0
      real(8),parameter :: subDrivUmax=1.0d0
      real(8),parameter :: subDrivRadi=0.025d0
      real(8),parameter :: subDrivTemp=306.15d0
      character*1,parameter :: subDrivDirecTag='y'

      real*8            :: yplus,r1,r2
      integer           :: n,i,j,idm
      integer,save      :: iflag=-1,nmax,ifl2,ifl3
      real*8 ,save      :: ds(7,128)
      real(8) :: velProf,eqivVelPeak
!
      real(8),parameter :: prs1=1.d0
      real(8),parameter :: prs5=1.0d-1
      real(8),parameter :: dens1=1.d0
      real(8),parameter :: dens5=1.25d-1
      real(8),parameter :: vel1=0.d0
      real(8),parameter :: vel5=0.d0
      real(8),parameter :: x0=0.5d0
!--------------------------------------------------------------------- 
! 
      real*8,parameter  :: PI=3.141592653589793238,z0=20.d0
      real*8 :: dum1,dum2,dum(3),unit(3),vvv(3),D,rrr
!---------------------------------------------------------------------

! --- [local entities] 
!=====================================================================
!      real*8            :: yplus,r1,r2
!      integer           :: n,i,j,idm
!      integer,save      :: iflag=-1,nmax,ifl2,ifl3
!      real*8 ,save      :: ds(7,128)
!      real(8) :: velProf,eqivVelPeak,Ret=395.d0,turbMagFactor=4.d0
      real(8) :: dist,Ret=395.d0
!=====================================================================
!  -- read dns.dat

      
! moving piston 
!      if(ipiston==1) then
        dum(1)=x
        dum(2)=y
        dum(3)=z
        unit(:)=0.d0
        unit(mv_dir)=omeini
        vvv(1:3)=0.d0
        rrr=sqrt(x**2+y**2)
        call cal_AXBEQC(unit,dum,vvv,D)
        u=vvv(1)
        v=vvv(2)
        w=0.d0
        p=0.d0
        aks(:)=0.d0
        return
!      endif


      return

      end subroutine user_ini
