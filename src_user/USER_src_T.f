!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_src_T(deltt,iter,time,
     &           ICV,IMAT_U,iphs,ncomp,x_cv,y_cv,z_cv,
     &           u,v,w,dens,comp,vol,
     &           T,heat)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     iter        : time step
!     deltt       : dt
!     time        : time
!
!     ICV         : Index of Control Volume
!     IMATU       : Material No. of ICV
!     iphs        : No. of Euler two phase (=1 or 2)
!     ncomp       : Species number at ICV
!     u,v,w       : Velocities at ICV   [m/s]
!     dens        : Density at ICV    [kg/m^3]
!     comp        : Mass friction at ICV   [100%]
!     T           : Temperature at ICV [K]
!     heat        : Heater [W/m^3=J/s/m^3]
!     x_cv,y_cv,z_cv : coordinate of CV
!     vol         : volume of ICV
!
!---------------------------------------------------------------------
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time,vol
      integer,intent(in)    :: iter
      integer,intent(in)    :: ICV,IMAT_U,iphs,ncomp
      real*8 ,intent(in)    :: u,v,w,dens,comp(ncomp)
      real*8 ,intent(in)    :: x_cv,y_cv,z_cv
      real*8 ,intent(inout) :: T
      real*8 ,intent(out)   :: heat
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! --- [local entities]
!
      real*8 :: dmin=1.D10,dum
      real*8 :: x_ing=.0d0,y_ing=.0d0,z_ing=1.d0
      integer:: idum
      real(8),parameter :: pi=3.14
      integer,save :: first_step=0,ICV_ING
!
! --- (x_ing,y_ing,z_ing)=> ingnation point
!
      if(first_step==0.or.first_step==iter) then
        first_step=iter
      else
        dum=(z_ing-z_cv)
        if(abs(dum)>0.01d0.or.iter>650) return
        if(T<1800) then
          T=1800.0
        end if
      endif

      heat=0.d0
      return
!
      end subroutine user_src_T
