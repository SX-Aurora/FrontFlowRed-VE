!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_src_r(icall,ictl,deltt,iter,time,
     &           ICV,IMAT_U,iphs,ncomp,x_cv,y_cv,z_cv,
     &           u,v,w,T,p,dens,comp,vol,
     &           sinj_r,sinj_u,sinj_v,sinj_w,sinj_t,sinj_comp)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     icall       : No. of calling this subroutine in one time-step
!                   =1 : 1st Calling for defining injector CV number
!                        in every time-step(ictl=0) 
!                        or  just one time (ictl=1)
!                   =2 : 2nd Calling from enthalpy and Ys transport eq.
!                        in every time-step
!                   =3 : 3rd Calling from u,v,w transport eq.
!                        in every time-step
!                   =4 : 4th Calling from continuty eq.
!                        in every time-step
!                   =5 : 5th Calling for 2PHYS : NOT Finished
!                        in every time-step
!                   (Total number in one time-step: Three)
!
!     ictl        : controlling flag
!                   =1 : '[icall=1] calling' just only once
!                   =0 : '[icall=1] calling' will be carried in 
!                        every time-step
!     ICV         : Index of Control Volume
!     deltt       : dt
!     time        : time
!     IMATU       : Material No. of ICV
!     iphs        : No. of Euler two phase (=1 or 2)
!     ncomp       : Total Species number at ICV
!     vol         : volume of ICV
!     x_cv,y_cv,z_cv : coordinate of CV
!     u,v,w       : Velocities at ICV [m/s]
!     dens        : Density at ICV [kg/m^3]
!     T           : Temperature at ICV  [K]
!     p           : pressure  [N/m^3]
!     comp        : Mass friction at ICV [100%]
!     utauu       : Friction velocity of wall surface nearest to ICV
!     suu,suv,suw : Transport equations sorce terms of u,v,w [N/m^3]
!     sinj_r      : Source term of continuity Eq. : [kg/s/m^3]
!     sinj_u,sinj_v,sinj_w : Injection velocity of sinj_r : [m/s]
!     sinj_t      : Injection temperature of sinj_r : [K]
!     sinj_comp   : Mass Friction of Injection Chemical species [100%]
!                   SUM of sinj_comp(1:ncomp) Must be 1.0
!
!---------------------------------------------------------------------
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter,icall
      integer,intent(inout) :: ictl
      integer,intent(in)    :: ICV,IMAT_U,iphs,ncomp
      real*8 ,intent(in)    :: x_cv,y_cv,z_cv
      real*8 ,intent(in)    :: u,v,w,dens,T,comp(ncomp),vol
      real*8 ,intent(in)    :: p
      real*8 ,intent(out)   :: sinj_r
      real*8 ,intent(out)   :: sinj_u,sinj_v,sinj_w
      real*8 ,intent(out)   :: sinj_t
      real*8 ,intent(out)   :: sinj_comp(ncomp)
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! --- [local entities]
!
      integer,save      :: ICV_inj=-1
      real*8            :: dum1,dum2,dum3,dum4
!
!      real*8,parameter  :: x_inj=0.2d0,y_inj=1.d0,z_inj=3.1d0  !BackStep !zhang
!      real*8,parameter  :: x_inj=1.1d0,y_inj=0.d0,z_inj=0.5d0
!      real*8,parameter  :: x_inj=1.d0,y_inj=0.d0,z_inj=0.01d0
      real*8,parameter  :: x_inj=0.1d0,y_inj=0.d0,z_inj=0.2d0
!-3.0, 2.0, 3.1
!
!
      integer           :: idum
      real*8,save       :: dismin=1.d10
      real(8),parameter :: pi=3.14d0
!
      ictl=1
!
      if(icall==1) then
        if(ICV==1) then
          dismin=1.d10
        endif
        dum1=(x_cv-x_inj)**2+(y_cv-y_inj)**2+(z_cv-z_inj)**2
        if(dum1<dismin) then
          ICV_inj=ICV
          dismin=dum1
        endif
      endif
!
      if(icall/=1) then
        sinj_comp(:)=0.d0

        if(ICV==ICV_inj.and.iter==5) then
          sinj_comp(:)=0.d0
          sinj_r = dens!*vol !zhang
          sinj_u=0.d0
          sinj_v=0.d0
          sinj_w=0.d0
          sinj_t=300.d0
          sinj_comp(1)=1.d0
          print*,'SSSSSSSSSS',sinj_r*deltt*vol,dens*vol
        endif
      endif

      return
!
      end subroutine user_src_r
