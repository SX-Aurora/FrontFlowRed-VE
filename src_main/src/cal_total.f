!
!      subroutine cal_dspfnc()
!      subroutine cal_eddy_viscosity
!      subroutine cal_rhoccc
!      subroutine cal_rva()
!      subroutine cal_h2t()   : a5 + a7
!      subroutine cal_vofh2t  : a5 for VOF 
!      subroutine cal_t2hcp() : a5
!      subroutine cal_t2cp    : a7
!      subroutine cal_real
!      subroutine cal_t2h     : a7
!      subroutine cal_temp()  : 
!      subroutine cal_usetar()
!      subroutine cal_statis   
!      subroutine cal_t2hcp_dc 
!      subroutine cal_t2h_dc
!      subroutine cal_vofdens
!      subroutine cal_csfvof
!      subroutine cal_AXBEQC
!      subroutine cal_molefraction
!      subroutine cal_stickcoef
!      subroutine cal_stickrate
!      subroutine cal_injcoord
!      subroutine CONTROL_RPM
!      subroutine cal_flamelet   !vector=>OK
!      subroutine cal_dmpbv      !vector=>OK
!      subroutine cal_trbbv      !vector=>OK
!      subroutine cal_bv         !vector=>OK
!      subroutine cal_prpgf      !vector=>OK
!      subroutine cal_cavi_h2t
!      subroutine cal_Dynamic_SGS
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
      subroutine cal_dspfnc(
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  LCYCOLD,wifsld,vctr,
     &  grdc,dtau,rmut,rk,vel,dspf)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     grdc(:,:,1)=>dtau(:,:)
!     grdc(:,1,2)=>dutau(:)
! 1. Calculate dissipation function: 2*Sij*Sij
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_metrix,only    : dutau  => W1K9
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: LVEDGE(2,MXCVFAC)
      integer,intent(in)   :: LBC_SSF( MXSSFBC)
      integer,intent(in)   :: LCYCSF(  MXSSFBC)
      INTEGER,INTENT(IN)   :: MAT_CV(  MXALLCV)
      INTEGER,INTENT(IN)   :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)   :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)   :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)   :: mat_cal(  0:MXMAT)
      integer,intent(in)   :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)   :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)   :: wiface(  MXCVFAC)
      real*8 ,intent(in)   :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)   :: CVVOLM(  MXALLCV)
      real*8 ,intent(in)   :: SFCENT(3,MXCVFAC)
      real*8 ,intent(in)   :: rmut  (  MXALLCV)
      real*8 ,intent(in)   :: rk    (  MXALLCV)
      real*8 ,intent(in)   :: vel   (  MXALLCV,3)
      real*8 ,intent(out)  :: dspf  (  MXALLCV)
      real*8 ,intent(inout):: dtau  (  MXALLCV,3)
      real*8 ,intent(inout):: grdc  (  MXALLCV,3,3)
      integer,intent(in)   :: LCYCOLD (MXSSFBC_SLD)
      real*8 ,intent(in)   :: wifsld  (MXSSFBC_SLD)
      integer,intent(in)   :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      real*8,parameter :: r2p3=2.d0/3.d0
      real*8  :: tau(3,3),vv(3),dxx,dyx,dzx,dll
      integer :: IMAT,IIMAT
      integer :: i,j,k,l,m,n,IV
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: ICVS,ICVE,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICVL,ICV,IDC,ICVP,IDCP
      real*8  :: dd1,dd2,rmuf,rkf,dum1,wi1,wi2
      real*8  :: grx,gry,grz,dx,dy,dz,gf0,dl
      real*8  :: grdf(3,3)
      integer :: IDIM
!
!-< 1. Calculate velocity gradient at cell center >- 
!
      call grad_cell(3,1,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vel,grdc)
!
!-< 2. Set dummy cell >-
!
      call dc_vgrad(MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,
     &              LCYCOLD,wifsld,     
     &              grdc)
!
!-< 3. Calculate velocity gradient at cell face >-
!
!-< 4. Calculate dissipation function >-
!
! --- grdc(:,:,1) => dtau; grdc(1,:,1) => dutau
!
      dtau(1:MXALLCV,1:3)=0.d0
      dutau(:)=0.d0
      do 400 IIMAT=1,NMAT    !ICF=1,NCVFAC
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 400
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do 404 ICFL=ICFS,ICFE
      ICVLA=LVEDGE(1,ICFL)
      ICVLB=LVEDGE(2,ICFL)
      dd1=rmut(ICVLA)
      dd2=rmut(ICVLB)
      wi1=wiface(ICFL)
      wi2=1.d0-wiface(ICFL)
      rmuf=dd1*dd2/(dd1*wi1+dd2*wi2)
!
      dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
      dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
      dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
      dl=dx*dx+dy*dy+dz*dz+SML
!
      do 405 IDIM=1,3
      grx=wi1*grdc(ICVLA,1,IDIM)+wi2*grdc(ICVLB,1,IDIM)
      gry=wi1*grdc(ICVLA,2,IDIM)+wi2*grdc(ICVLB,2,IDIM)
      grz=wi1*grdc(ICVLA,3,IDIM)+wi2*grdc(ICVLB,3,IDIM)
      gf0=(
     &    (vel(ICVLB,IDIM)-vel(ICVLA,IDIM))
     &   -(grx*dx+gry*dy+grz*dz)
     &    )/dl
      grdf(1,IDIM)=grx+dx*gf0
      grdf(2,IDIM)=gry+dy*gf0
      grdf(3,IDIM)=grz+dz*gf0
 405  continue
!
      do 401 IV=1,3
      tau(iv,1)=rmuf*(grdf(iv,1)+grdf(1,iv))
      tau(iv,2)=rmuf*(grdf(iv,2)+grdf(2,iv))
      tau(iv,3)=rmuf*(grdf(iv,3)+grdf(3,iv))
      vv(iv)=wi1*vel(ICVLA,iv)+wi2*vel(ICVLB,iv)
  401 continue
      rkf=r2p3*(wi1*rk(ICVLA)+wi2*rk(ICVLB))
      do 402 i=1,3
      tau(i,i)=tau(i,i)-rkf
  402 continue
      do 403 i=1,3
      dum1=(SFAREA(1,ICFL)*tau(i,1)
     &     +SFAREA(2,ICFL)*tau(i,2)
     &     +SFAREA(3,ICFL)*tau(i,3))*SFAREA(4,ICFL)
      dtau(ICVLA,i)=dtau(ICVLA,i)+dum1
      dtau(ICVLB,i)=dtau(ICVLB,i)-dum1
  403 continue
      dum1=(SFAREA(1,ICFL)*
     &     (vv(1)*tau(1,1)
     &     +vv(2)*tau(2,1)
     &     +vv(3)*tau(3,1))
     &     +SFAREA(2,ICFL)*
     &     (vv(1)*tau(1,2)
     &     +vv(2)*tau(2,2)
     &     +vv(3)*tau(3,2))
     &     +SFAREA(3,ICFL)*
     &     (vv(1)*tau(1,3)
     &     +vv(2)*tau(2,3)
     &     +vv(3)*tau(3,3))
     &    )*SFAREA(4,ICFL)
      dutau(ICVLA)=dutau(ICVLA)+dum1
      dutau(ICVLB)=dutau(ICVLB)-dum1
 404  continue
 400  continue
!
      do 410 IIMAT=1,NMAT   !ICV=1,NCV
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 410
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      dspf(ICVL)=max(0.d0,dutau(ICVL)-
     &         (vel(ICVL,1)*dtau(ICVL,1)
     &         +vel(ICVL,2)*dtau(ICVL,2)
     &         +vel(ICVL,3)*dtau(ICVL,3)))
     &  /CVVOLM(ICVL)
      enddo
 410  continue
!
!
      return
!
      end subroutine cal_dspfnc
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_h2t(mph,NCVX,   !15768
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,hhh,yys,tmp,
     &  prs,rho)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_species,only  : sw
      use module_material,only : cpsld,rmdsld,nofld
      use module_io,      only : ifll,ifle
      use module_metrix  ,only : yi=>dum_c1,ahk,ahk2,acpk_w
      use module_metrix,only   : W1K9,W1K8!,W1K7
      use module_model,only    : ical_vect,ical_dens
      use module_species,only  : gascns
!
! 1. Calculate temperature from enthalpy & mass fraction
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: NCVX,mph
      integer,intent(in)    :: MAT_NO(    0:MXMAT)
      integer,intent(in)    :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)    :: MAT_DCIDX( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(   0:MXMAT)
      real*8 ,intent(in)    :: hhh    (   MXALLCV)
      real*8 ,intent(in)    :: yys(       MXALLCV,mxcomp)
      real*8 ,intent(inout) :: tmp    (   MXALLCV)
      real*8 ,intent(in)    :: prs     (  MXALLCV)
      real*8 ,intent(in)    :: rho     (  MXALLCV)
!
! --- [local entities]
!
      integer,parameter :: lpmax=50
      real*8 ,parameter :: eps=1.d-10 !eps=1.d-10
      integer :: i,j,k,l,m,n,ll
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP,IMAT_U
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD,ierr1=0
      real*8  :: h0,tt,dt,hh,cp,y1,t1,hho
      real*8  :: uu,hmin,hmino,dTeps=1.d-5
!
!-------------------------------------------------
! --- 
!------------------------------------------------
!
      do 200 IIMAT=1,NMAT                 !ICV=1,NCVX
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT.gt.0) then
        IMAT_U=nofld(IMAT)
        if(sw) then
          if(ical_vect) then  !NOVECT
            call a7_V(ICVS,ICVE,IMAT,IMAT_U)
            if(NCVX>NCV) then        ! --- Dumy CV
              IDCS = MAT_DCIDX(IIMAT-1)+1
              IDCE = MAT_DCIDX(IIMAT)
              call a7_V(IDCS,IDCE,IMAT,IMAT_U)
            endif
          else
            if(ical_dens==4) then
            else
              call a7(ICVS,ICVE,IMAT,IMAT_U)
              if(NCVX>NCV) then        ! --- Dumy CV
                IDCS = MAT_DCIDX(IIMAT-1)+1
                IDCE = MAT_DCIDX(IIMAT)
                call a7(IDCS,IDCE,IMAT,IMAT_U)
              endif
            endif
          endif
        else
          if(ical_vect) then  !NOVECT
            call a5_V(mph,ICVS,ICVE,IMAT)
            if(NCVX>NCV) then        ! --- Dumy CV
              IDCS=MAT_DCIDX(IIMAT-1)+1
              IDCE=MAT_DCIDX(IIMAT)
              call a5_V(mph,IDCS,IDCE,IMAT)
            endif
          else
            call a5(mph,ICVS,ICVE,IMAT)
            if(NCVX>NCV) then        ! --- Dumy CV
              IDCS=MAT_DCIDX(IIMAT-1)+1
              IDCE=MAT_DCIDX(IIMAT)
              call a5(mph,IDCS,IDCE,IMAT)
            endif
          endif
        endif
      else                    ! solid
        call cal_solid(ICVS,ICVE)
        if(NCVX>NCV) then          ! --- Dumy CV
          IDCS = MAT_DCIDX(IIMAT-1)+1
          IDCE = MAT_DCIDX(IIMAT)
          call cal_solid(IDCS,IDCE)
        endif
      endif
 200  continue
!
      return
!//////////////////////////////////////////////////////////////////////
      contains
!
!=============================================================
      subroutine a5(mph,do_start,do_end,IMAT)
!=============================================================
      use module_species,only  : acpk,acpk2
      integer,intent(in) :: mph,do_start,do_end,IMAT
      real*8  :: dum1,dum2,dum3
!
      if(mph==1) then
        do ICOM=1,ncomp
        do n=1,5
        ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk(0,ICOM)
        enddo
        acpk_w=acpk
      else
        do ICOM=1,ncomp
        do n=1,5
        ahk(n,ICOM)=acpk2(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk2(0,ICOM)
        enddo
        acpk_w=acpk2
      endif
!
      dum1=1.d10
!
      do ICVL=do_start,do_end
        do ICOM=1,ncomp
        yi(ICOM)=yys(ICVL,ICOM)
        enddo
        tt=tmp(ICVL)
        do k=1,lpmax
          hh = 0.d0
          cp = 0.d0
          do ICOM=1,ncomp
!
            hh = hh+yi(ICOM)*(((((
     &        ahk(5,ICOM) *tt
     &       +ahk(4,ICOM))*tt
     &       +ahk(3,ICOM))*tt
     &       +ahk(2,ICOM))*tt
     &       +ahk(1,ICOM))*tt
     &       +ahk(0,ICOM))
!
            cp = cp+yi(ICOM)*((((
     &        acpk_w(5,ICOM )*tt
     &       +acpk_w(4,ICOM))*tt
     &       +acpk_w(3,ICOM))*tt
     &       +acpk_w(2,ICOM))*tt
     &       +acpk_w(1,ICOM))
          enddo
!
          dt = (hhh(ICVL)-hh)/cp
          tt = max(0.5d0*tt,tt+dt)
!
          if(abs(dt)<eps*tt) goto 254 
!
        enddo
        write(ifle,'(a)') '### program error -1- (cal_h2t: a5)'
        write(ifle,'(1X,a)') 
     & 'MSG : mph,IMAT,ICVL,hhh(ICVL),hh,cp,tt,tmp(ICVL)='
        write(ifle,'(1X,3I10,5E16.4)') 
     &   mph,IMAT,ICVL,hhh(ICVL),hh,cp,tt,tmp(ICVL)
        CALL FFRABORT(1,'dt too large for cal_h2t')
  254   continue
!        if(tt<dum1) then
!          dum1=tt 
!          dum2=hhh(ICVL) 
!          dum3=hh 
!        endif
        tmp(ICVL)=tt
      enddo
      end subroutine a5 

!=============================================================
      subroutine a5_V(mph,do_start,do_end,IMAT)
!=============================================================
      use module_species,only  : acpk,acpk2
      integer,intent(in) :: mph,do_start,do_end,IMAT
      integer,parameter :: nstep=5
!
      if(mph==1) then
        do ICOM=1,ncomp
        do n=1,5
        ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk(0,ICOM)
        enddo
        acpk_w=acpk
      else
        do ICOM=1,ncomp
        do n=1,5
        ahk(n,ICOM)=acpk2(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk2(0,ICOM)
        enddo
        acpk_w=acpk2
      endif
!
      do k=1,nstep
      W1K9(do_start:do_end)=0.d0   !hh
      W1K8(do_start:do_end)=0.d0   !cp
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      tt=tmp(ICVL)
      W1K9(ICVL)=W1K9(ICVL)+yys(ICVL,ICOM)*(((((
     &        ahk(5,ICOM) *tt
     &       +ahk(4,ICOM))*tt
     &       +ahk(3,ICOM))*tt
     &       +ahk(2,ICOM))*tt
     &       +ahk(1,ICOM))*tt
     &       +ahk(0,ICOM))
      W1K8(ICVL)=W1K8(ICVL)+yys(ICVL,ICOM)*((((
     &        acpk_w(5,ICOM )*tt
     &       +acpk_w(4,ICOM))*tt
     &       +acpk_w(3,ICOM))*tt
     &       +acpk_w(2,ICOM))*tt
     &       +acpk_w(1,ICOM))
      enddo
      enddo
      do ICVL=do_start,do_end
      tt=tmp(ICVL)
      dt=(hhh(ICVL)-W1K9(ICVL))/W1K8(ICVL)
      tt=max(0.5d0*tt,tt+dt)
      tmp(ICVL)=tt
      enddo
      enddo
!
      return
      end subroutine a5_V
!
!=================================================
      subroutine a7(do_start,do_end,IMAT,IMAT_U)
!=================================================
      use module_species,only : c_p,enthalpy,
     &                          href,enthalpy_ref,
     &                          act
      integer,intent(in) :: do_start,do_end,IMAT,
     &                      IMAT_U
      real*8 :: tt1,dum1
!
      yi(:)=0.d0
      do ICVL=do_start,do_end
        do ICOM=1,ncomp
        yi(ICOM)=yys(ICVL,ICOM)*dble(ACT(ICOM))  ! mass fraction
        enddo
        tt=tmp(ICVL)                      ! temperature [K]
        hho=hhh(ICVL)
!
        do k=1,lpmax
          hh=enthalpy_ref(ncomp,tt,yi)
          cp=c_p(ncomp,tt,yi)
          dT = (hho-hh)/cp
          tt=max(0.5*tt,tt+dt)
          if(abs(dt)<dTeps*tt ) exit 
        enddo
!
        if(k==lpmax+1) then
          call error_message(IMAT,ICVL,ncomp,yi,dt,hh,hho,cp,tt,k)
          cycle
        endif
        tmp(ICVL) = tt
      enddo
      end subroutine a7
!

!=================================================
      subroutine a7_real(do_start,do_end,IMAT,IMAT_U)
!=================================================
      use module_species,only : Tc,Pc,KIJ,omga,r_wm,Ri
      use module_species,only : c_p,href,rho_real_PR,
     &                          enthalpy_ref,
     &                          act
      integer,intent(in) :: do_start,do_end,IMAT,IMAT_U
      real*8 :: tt1,dfunc,func
      real*8           :: TTT,PPP,VVV,
     &                    dum1,dum2,dum3,hho1
      integer :: j,NR,itr,mx_ycomp
      real*8  :: tol=1.0d-8,DVDP_R
      real*8  :: dtt,mx_y,delt_cp,delt_H,delt_C
!
! --- ----------
!
! ---
!
      
      end subroutine a7_real

!===================================================
      subroutine a7_V(do_start,do_end,IMAT,IMAT_U)
!===================================================
      use module_species,only : href,Ri,t3,a7
      integer,parameter :: nstep=5
      integer,intent(in) :: do_start,do_end,IMAT,IMAT_U
      real*8 :: tt1,ttt,hh,cp,dT,dum1,h_reff,c_pi
      integer :: rj
!
      do k=1,nstep
      W1K9(do_start:do_end)=0.d0   !enthalpy_ref
      W1K8(do_start:do_end)=0.d0
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      ttt=tmp(ICVL)
      rj=(3-sign(1,int(ttt-t3(2,ICOM))))/2
!
      h_reff=a7(6,rj,icom)
     &   +ttt*(a7(1,rj,icom)
     &   +ttt*(a7(2,rj,icom)*0.5d0
     &   +ttt*(a7(3,rj,icom)*1.d0/3.d0
     &   +ttt*(a7(4,rj,icom)*0.25d0
     &   +ttt*(a7(5,rj,icom)*0.2d0)))))-href(icom)
      W1K9(ICVL)=W1K9(ICVL)+yys(ICVL,ICOM)*h_reff*Ri(icom)
      c_pi=a7(1,rj,ICOM)
     &   +ttt*(a7(2,rj,ICOM)
     &   +ttt*(a7(3,rj,ICOM)
     &   +ttt*(a7(4,rj,ICOM)
     &   +ttt* a7(5,rj,ICOM))))
      W1K8(ICVL)=W1K8(ICVL)+c_pi*yys(ICVL,ICOM)*Ri(ICOM)
      enddo
      enddo
!
      do ICVL=do_start,do_end
      ttt=tmp(ICVL)
      dt=(hhh(ICVL)-W1K9(ICVL))/W1K8(ICVL)
      ttt=max(0.5d0*ttt,ttt+dt)
      tmp(ICVL)=ttt  !tmp(ICVL)+!(1.d0-W1K7(ICVL))*dt
      enddo
!
      enddo
!
      end subroutine a7_V
!      
!
!==============================================
      subroutine cal_solid(ICVS,ICVE)
!==============================================
      use module_material,only : cpsld
      integer,intent(in) :: ICVS,ICVE
!
      tmp(ICVS:ICVE) = hhh(ICVS:ICVE)/cpsld(-IMAT)
      end subroutine cal_solid
!
!==============================================
      subroutine error_message(IMAT,ICVL,ncomp,yy,dt,h,h0,cp,tt,k)
!==============================================
      use module_model,   only : ical_prt
      integer,intent(in) :: ICVL,k,IMAT,ncomp
      real*8, intent(in) :: dt,h,h0,cp,tt,yy(ncomp)
      CHARACTER*6  :: for1
      character*50:: formt
!
      write(for1,'(I3)') ncomp
      write(formt,'(a,a,a,a)') "'(a,","1x,",trim(for1),"(E15.8))'"
!
      write(ifle,'(1x,a)') "ERR: program error -2- subroutine a7_real"
      write
     & (ifle,
     & "(1x,a,i3,a,i8,a,e10.2,a,e10.2,a,e10.2,
     &   a,e10.2,a,e10.2,a,e10.2,i4)")
     &     "IMAT_U= ",IMAT_U," ICVL =",ICVL,"  dT =",dt,
     &     " ttt=", tt,  " h =", h,
     &     "  h0 =",h0,"  cp =",cp,"  tmp0 =",tmp(ICVL),k
      do icom=1,ncomp
      write(ifle,'(1x,a,I4,E15.8)') 'Ys(i)=',icom,yy(icom) 
      enddo

!      if(ical_prt<1) call FFRABORT(1,'a7_real MSG: dt too large')
      end subroutine error_message
!
      end subroutine cal_h2t
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_h2t1(mph,NCVX,
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,hhh,yys,tmp)   !15769
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_species,only  : sw
      use module_material,only : cpsld,rmdsld
      use module_io,      only : ifll,ifle
      use module_metrix  ,only : yy=>ys,ahk,ahk2,acpk_w
!
! 1. Calculate temperature from enthalpy & mass fraction
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: NCVX,mph
      integer,intent(in)    :: MAT_NO(    0:MXMAT)
      integer,intent(in)    :: MAT_CVEXT( 0:MXMAT)
      integer,intent(in)    :: MAT_DCIDX( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(   0:MXMAT)
      real*8 ,intent(in)    :: hhh    (   MXALLCV)
      real*8 ,intent(in)    :: yys(       MXALLCV,mxcomp)
      real*8 ,intent(inout) :: tmp    (   MXALLCV)
!
! --- [local entities]
!
      integer,parameter :: lpmax=100
      real*8 ,parameter :: eps=1.d-10 !eps=1.d-10
      integer :: i,j,k,l,m,n,ll
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD,ierr1=0
      real*8  :: h0,tt,dt,hh,cp,y1,t1,hho
      real*8  :: uu,hmin,hmino,dTeps=1.d-5 !jiang
!
!-------------------------------------------------
! --- 
!-------------------------------------------------
!
      do 200 IIMAT=1,NMAT                 !ICV=1,NCVX
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT.gt.0) then
        if(sw) then
          call a7(ICVS,ICVE,IMAT)
          if(NCVX>NCV) then        ! --- Dumy CV
            IDCS = MAT_DCIDX(IIMAT-1)+1
            IDCE = MAT_DCIDX(IIMAT)
            call a7(IDCS,IDCE,IMAT)
          endif
        else
          call a5(mph,ICVS,ICVE,IMAT)
          if(NCVX>NCV) then        ! --- Dumy CV
            IDCS=MAT_DCIDX(IIMAT-1)+1
            IDCE=MAT_DCIDX(IIMAT)
            call a5(mph,IDCS,IDCE,IMAT)
          endif
        endif
      else                    ! solid
        call cal_solid(ICVS,ICVE)
        if(NCVX>NCV) then          ! --- Dumy CV
          IDCS = MAT_DCIDX(IIMAT-1)+1
          IDCE = MAT_DCIDX(IIMAT)
          call cal_solid(IDCS,IDCE)
        endif
      endif
 200  continue
!
      return
!//////////////////////////////////////////////////////////////////////
      contains
!
!=============================================================
      subroutine a5(mph,do_start,do_end,IMAT)
!=============================================================
      use module_species,only  : acpk,acpk2
      integer,intent(in) :: mph,do_start,do_end,IMAT
!
      if(mph==1) then
        do ICOM=1,ncomp
        do n=1,5
        ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk(0,ICOM)
        enddo
        acpk_w=acpk
      else
        do ICOM=1,ncomp
        do n=1,5
        ahk(n,ICOM)=acpk2(n,ICOM)/dble(n)
        enddo
        ahk(0,ICOM)=acpk2(0,ICOM)
        enddo
        acpk_w=acpk2
      endif
!
      do ICVL=do_start,do_end
        do ICOM=1,ncomp
        yy(ICOM)=yys(ICVL,ICOM)
        enddo
        tt=tmp(ICVL)
        do k=1,lpmax
          hh = 0.d0
          cp = 0.d0
          do ICOM=1,ncomp
!
            hh = hh+yy(ICOM)*(((((
     &        ahk(5,ICOM) *tt
     &       +ahk(4,ICOM))*tt
     &       +ahk(3,ICOM))*tt
     &       +ahk(2,ICOM))*tt
     &       +ahk(1,ICOM))*tt
     &       +ahk(0,ICOM))
!
            cp = cp+yy(ICOM)*((((
     &        acpk_w(5,ICOM )*tt
     &       +acpk_w(4,ICOM))*tt
     &       +acpk_w(3,ICOM))*tt
     &       +acpk_w(2,ICOM))*tt
     &       +acpk_w(1,ICOM))
          enddo
!
          dt = (hhh(ICVL)-hh)/cp
          tt = max(0.5d0*tt,tt+dt)


          if(abs(dt)<eps*tt) goto 254
!
        enddo
        write(ifle,'(a)') '### program error -1- (cal_h2t: a5)'
        write(ifle,'(1X,a)') 
     & 'MSG : mph,IMAT,ICVL,hhh(ICVL),hh,cp,tt,tmp(ICVL)='
        write(ifle,'(1X,3I10,5E16.4)') 
     &   mph,IMAT,ICVL,hhh(ICVL),hh,cp,tt,tmp(ICVL)
        CALL FFRABORT(1,'dt too large for cal_h2t')
  254   continue
        tmp(ICVL)=tt
      enddo
      end subroutine a5


!revised by jiang yuyan 
!==============================================
      subroutine a7(do_start,do_end,IMAT)
!==============================================
      use module_species,only : c_p,enthalpy,href,enthalpy_ref
	use module_species,only : t_bi1,t_bi2	!jiang
      integer,intent(in) :: do_start,do_end,IMAT
      real*8 :: tt1,dum1, tt2,hh1,hh2
!
      do ICVL=do_start,do_end

	  tt1=	t_bi1	!+ 300
	  tt2=	t_bi2	!t_bi2
	  	 
	  do ICOM=1,ncomp
          yy(ICOM)=yys(ICVL,ICOM)           ! mass fraction
	  enddo
        tt=tmp(ICVL)                        ! temperature [K]
        hho=hhh(ICVL)

!revised by Jiang yuyan
	  
	  hh1=enthalpy_ref(ncomp,tt1,yy)
	  hh2=enthalpy_ref(ncomp,tt2,yy)

	  if(hho.lt.hh1.and.hho.gt.hh2) then
	         
	  do k=1,lpmax   ! lpmax=100
          hh=enthalpy_ref(ncomp,tt,yy)
          cp=c_p(ncomp,tt,yy)
          dT = (hho-hh)/cp
          tt=max(0.5*tt,tt+dt)
          if(abs(dt)<dTeps*tt ) exit
c          dt=(hho-hh)/cp
c          tt1=tt+dt
c          if(k>50)  tt1=tt+abs(dt)
c          tt=max(0.5*tt,tt1)
c          if(abs(dt)<eps*tt ) exit
        enddo
 
        if(k==lpmax+1) then
          call error_message(IMAT,ICVL,dt,hh,hho,cp,tt,k)
        endif

	  tmp(ICVL) = tt

	  end if
      enddo
      end subroutine a7
!


!
!==============================================
      subroutine a7_zhang(do_start,do_end,IMAT)
!==============================================
      use module_species,only : c_p,enthalpy,href,enthalpy_ref
      integer,intent(in) :: do_start,do_end,IMAT
      real*8 :: tt1,dum1
!
      do ICVL=do_start,do_end
        do ICOM=1,ncomp
          yy(ICOM)=yys(ICVL,ICOM)           ! mass fraction
        enddo
        tt=tmp(ICVL)                        ! temperature [K]
        hho=hhh(ICVL)
        do k=1,lpmax   ! lpmax=100
          hh=enthalpy_ref(ncomp,tt,yy)
          cp=c_p(ncomp,tt,yy)
          dT = (hho-hh)/cp
          tt=max(0.5*tt,tt+dt)
          if(abs(dt)<dTeps*tt ) exit
!          dt=(hho-hh)/cp
!          tt1=tt+dt
!          if(k>50)  tt1=tt+abs(dt)
!          tt=max(0.5*tt,tt1)
!          if(abs(dt)<eps*tt ) exit
        enddo
        if(k==lpmax+1) then
          call error_message(IMAT,ICVL,dt,hh,hho,cp,tt,k)
        endif
        tmp(ICVL) = tt
      enddo
      end subroutine a7_zhang
!
!==============================================
      subroutine cal_solid(ICVS,ICVE)
!==============================================
      use module_material,only : cpsld
      integer,intent(in) :: ICVS,ICVE
!
      tmp(ICVS:ICVE) = hhh(ICVS:ICVE)/cpsld(-IMAT)
      end subroutine cal_solid
!
!==============================================
      subroutine error_message(IMAT,ICVL,dt,h,h0,cp,tt,k)
!==============================================
      integer,intent(in) :: ICVL,k,IMAT
      real*8, intent(in) :: dt,h,h0,cp,tt
!
      write(ifle,'(a)') "ERR: program error -2- MAT (cal_h2t: a7)"
      write(6,"(a,i5,a,i5,a,e10.2,a,e10.2,a,e10.2,a,e10.2,a,e10.2,i4)")
     &     "IMAT =",IMAT," ICVL =",ICVL,"  dT =",dt,"  h =",h,
     &     "  h0 =",h0,"  cp =",cp,"  tmp =",tt,k
      call FFRABORT(1,'cal_h2t MSG: dt too large')
      end subroutine error_message
!
      end subroutine cal_h2t1


!
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_eddy_viscosity(iter,deltt,time,
     &  LVEDGE,LBC_SSF,LCYCSF,LFUTAU,
     &  MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,MAT_INDEX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,DISALL,FRSTCV,GRDC,OPPANG,
     &  rho,vel,rmu,rmut,rmut2,aks,tmp,utau,wifsld,LCYCOLD,vctr,FIELD_U,
     &  yplusf,vflc)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!15800
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,      only : ifll,ifle
      use module_les,only      : csles,apls
      use module_turbparm,only : akappa,awfnc,yplsm,
     &                           Rey_s,Amu,dRey
      use module_boundary,only : kdbcnd,kvfslp,kvlglw,LBC_INDEX,nbcnd,
     &                           kvnslp,MAT_BCIDX,kvrogf,rghnes,kvmodl,
     &                           kvEWF,kvsataW,nobcnd
      use module_model,   only : idrdp,incomp,mach0,icaltb,KLES,
     &                           ke,sles,dles,lles,noturb,RSM,ke_low,
     &                           RNG,CHEN,SDES,cnst,ical_EWT,KE2S,u_func
      use module_rans    ,only : cmu,Ck
      use module_rans,only     : akdes,cb1,cb2,sgmDES,cw1,cw2,cw3,cv1,
     &                           cDES,L2s
      use module_scalar,  only : ike,ides,calgeq,calxi,POLYCDP,mCD,
     &                           fitvalue,ikles
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_material,only : turbmu,turbmu2
      use module_Euler2ph,only : ieul2ph
      use module_metrix,only   : rghnesx,multiR
      use module_usersub,only  : usrno,usryes,sata_ruf_wall
!
! 1. Calculate eddy viscosity of Smagorinsky model
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)  :: deltt,time
      integer,intent(in)  :: iter
      integer,intent(in)  :: LVEDGE (2,MXCVFAC)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      INTEGER,INTENT(IN)  :: MAT_NO(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_CV(   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_DCIDX(0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_INDEX(0:MXMAT)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      integer,intent(in)  :: LFUTAU (    MXCV )
!
      real*8 ,intent(in)  :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)  :: SFCENT(3, MXCVFAC)
      real*8 ,intent(in)  :: wiface(   MXCVFAC)
      real*8 ,intent(in)  :: CVCENT(3, MXALLCV)
      real*8 ,intent(in)  :: CVVOLM(   MXALLCV)
      real*8 ,intent(in)  :: DISALL(   MXCV   )
      real*8 ,intent(in)  :: FRSTCV(   MXSSFBC)
      real*8 ,intent(inout) :: grdc(   MXALLCV,3,3)
      real*8 ,intent(in)  :: rho   (   MXALLCV)
      real*8 ,intent(in)  :: vel   (   MXALLCV,3)
      real*8 ,intent(in)  :: rmu   (   MXALLCV)
      real*8 ,intent(out) :: rmut  (   MXALLCV)
      real*8 ,intent(out) :: rmut2 (   MXALLCV2)
      real*8 ,intent(in)  :: aks   (  MXALLCVR,MXRANS)
      real*8 ,intent(inout):: utau (0: MXSSFBC)
      real*8 ,intent(in)  :: tmp   (MXALLCV)
      real*8 ,intent(in)  :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)  :: LCYCOLD( MXSSFBC_SLD)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(inout) :: FIELD_U  (MXCV_D,NFLID)
      real*8 ,intent(inout):: yplusf(MXCV_F)
      real*8 ,intent(inout):: vflc(MXCV_F)
      real*8 ,intent(in)  :: OPPANG(   MXSSFBC_SLD)
!
! --- [local entities]
!
      real*8 ,parameter :: r1pn=1.d0/3.d0,R26=1.D0/26.D0
      real*8  :: dspd,utaul,dum1,dum2

      real*8  :: Rey,rmd_E,lmu,cl,AAA,CD1,CD2,lll,CDx

      real*8  :: ux,uy,uz,up,yp,rnu,uplu,ux1,uy1,uz1
      real*8  :: SIJ11,SIJ22,SIJ33,SIJ12,SIJ13,SIJ23,fmu,
     &           Ry,rans_k,rans_e,Rt,SML1=1.d-10
      integer :: no,nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      real*8  :: del,uutau,rkkk,xxx,fv1des
      real*8  :: ustar,duplus,vk,yre,uplus,yplus,awfnc1,eps
      real*8  :: cfflc=0.15d0,dum3
      integer,save :: icall=0,flag=0,ierr=0,user_rfn=0
      real*8  :: t_cell,t_wall,rhoo,rmu_m,rmu_t,area,cntr(3)
      integer,save,allocatable  :: flagsata(:)
!
!------------------------------------------------
! --- Calculating Wall fraction velocity: utau() 
!------------------------------------------------
!
      utau=0.d0
      utau(0)=0.d0
!      
      uutau=0.d0
      do 100 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      kdv=kdbcnd(1,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kdv.eq.kvlglw.or.kdv==kvmodl) then 
        IF(ical_vect) then  !NOVECT
!          if(icaltb==ke.or.icaltb==RNG.or.icaltb==CHEN) then
!            do IBFL=IBFS,IBFE
!            ICFL=LBC_SSF(IBFL)
!            ICV=LVEDGE(1,ICFL)
!            IDC=LVEDGE(2,ICFL)
!            ux=vel(ICV,1)-vel(IDC,1)
!            uy=vel(ICV,2)-vel(IDC,2)
!            uz=vel(ICV,3)-vel(IDC,3)
!            up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
!            ux=ux-up*SFAREA(1,ICFL)
!            uy=uy-up*SFAREA(2,ICFL)
!            uz=uz-up*SFAREA(3,ICFL)
!            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
!            yp=FRSTCV(IBFL)
!            rnu=rmu(ICV)/rho(ICV)
!            utaul=cmu**0.25d0*dsqrt(aks(ICV,ike(1)))
!            utau(IBFL)=utaul
!            enddo
!          else !else=>if(icaltb==ke
!CIDR NODEP
            vk=1.d0/akappa
            uplus=one
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            ux=vel(ICV,1)-vel(IDC,1)
            uy=vel(ICV,2)-vel(IDC,2)
            uz=vel(ICV,3)-vel(IDC,3)
            up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
            ux=ux-up*SFAREA(1,ICFL)
            uy=uy-up*SFAREA(2,ICFL)
            uz=uz-up*SFAREA(3,ICFL)
            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
            yp=FRSTCV(IBFL)
            rnu=rmu(ICV)/rho(ICV)
            utaul=utau(IBFL)
!---------------
! --- log-law
!---------------
            yre=abs(up*yp/rnu)
            if(yre<1.d0) then
              utaul=dsqrt(rnu*up/yp)
            else
              ustar=vk*log(yre)+awfnc
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
              uplus=uplus+duplus
              utaul=up/(uplus+SML)
            endif
            yplus=yp*utaul/(rnu+SML)
            utau(IBFL)=utaul
            if(yplus<yplsm) then
              utau(IBFL)=dsqrt(rnu*up/yp)
            endif
            enddo
!          endif
        else   !IF(ical_vect) then
          if(icaltb==KE2S.or.icaltb==ke.or.icaltb==RNG.or.icaltb==CHEN) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            ux=vel(ICV,1)-vel(IDC,1)
            uy=vel(ICV,2)-vel(IDC,2)
            uz=vel(ICV,3)-vel(IDC,3)
            up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
            ux=ux-up*SFAREA(1,ICFL)
            uy=uy-up*SFAREA(2,ICFL)
            uz=uz-up*SFAREA(3,ICFL)
            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
            yp=FRSTCV(IBFL)
            rnu=rmu(ICV)/rho(ICV)
            utaul=cmu**0.25d0*dsqrt(aks(ICV,ike(1)))
            utaul=utau(IBFL)  ! nikki
            call cal_ustar(yp,rnu,akappa,awfnc,up,utaul,2,1) !
            utau(IBFL)=utaul
            enddo
          else
            do 1100 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            ux=vel(ICV,1)-vel(IDC,1)
            uy=vel(ICV,2)-vel(IDC,2)
            uz=vel(ICV,3)-vel(IDC,3)
            up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
            ux=ux-up*SFAREA(1,ICFL)
            uy=uy-up*SFAREA(2,ICFL)
            uz=uz-up*SFAREA(3,ICFL)
            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
            yp=FRSTCV(IBFL)
            rnu=rmu(ICV)/rho(ICV)
            utaul=utau(IBFL)
            call cal_ustar(yp,rnu,akappa,awfnc,up,utaul,2,2) !1
            utau(IBFL)=utaul
 1100       enddo
          endif
        endif
      elseif(kdv==kvsataW) then
        if(flag==0.and.sata_ruf_wall==usryes) then
          flag=1
          ALLOCATE(rghnesx(MXSSFBC),stat=ierr)
          if(ierr/=0) call FFRABORT(1,'ALLOCATE rghnesx error')
          ALLOCATE(flagsata(nbcnd),stat=ierr)
          if(ierr/=0) call FFRABORT(1,'ALLOCATE flagsata error')
          rghnesx(:)=0.d0
          flagsata(:)=0
        endif
        if(sata_ruf_wall==usryes) then
          no=nobcnd(nb)
          if(flagsata(nb)==0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            rghnesx(IBFL)=rghnes(nb)
            ux=vel(ICV,1)
            uy=vel(ICV,2)
            uz=vel(ICV,3)
            ux1=-SFAREA(1,ICFL)
            uy1=-SFAREA(2,ICFL)
            uz1=-SFAREA(3,ICFL)
            area=SFAREA(4,ICFL)
            cntr(1)=SFCENT(1,ICFL)
            cntr(2)=SFCENT(2,ICFL)
            cntr(3)=SFCENT(3,ICFL)
            t_cell=tmp(ICV)
            t_wall=tmp(IDC)
            rhoo=rho(ICV)
            rmu_m=rmu(ICV)
            rmu_t=rmut(ICV)
            
            del=FRSTCV(IBFL)
            dum1=rghnes(nb)
            call user_rough_wall
     &      (flagsata(nb),no,iter,deltt,time,
     &      ux,uy,uz,
     &      ux1,uy1,uz1,area,cntr(1),cntr(2),cntr(3),del,
     &      t_cell,t_wall,rhoo,rmu_m,rmu_t,dum1)
            rghnesx(IBFL)=dum1
          enddo
          endif
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ux=vel(ICV,1)-vel(IDC,1)
          uy=vel(ICV,2)-vel(IDC,2)
          uz=vel(ICV,3)-vel(IDC,3)
          up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
          ux=ux-up*SFAREA(1,ICFL)
          uy=uy-up*SFAREA(2,ICFL)
          uz=uz-up*SFAREA(3,ICFL)
          up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
          yp=FRSTCV(IBFL)
          utau(IBFL)=up*0.42d0/log(yp/rghnesx(IBFL)+1.d0)
          enddo
!          
        else
          dum1=rghnes(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ux=vel(ICV,1)-vel(IDC,1)
          uy=vel(ICV,2)-vel(IDC,2)
          uz=vel(ICV,3)-vel(IDC,3)
          up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
          ux=ux-up*SFAREA(1,ICFL)
          uy=uy-up*SFAREA(2,ICFL)
          uz=uz-up*SFAREA(3,ICFL)
          up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
          yp=FRSTCV(IBFL)
          utau(IBFL)=up*0.42d0/log(yp/dum1+1.d0)
          enddo
        endif
      elseif(kdv==kvrogf) then
!CIDR NODEP
        vk=1.d0/0.42d0
        awfnc1=8.5d0
        dum1=rghnes(nb)          
        uplus=one
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ux=vel(ICV,1)-vel(IDC,1)
          uy=vel(ICV,2)-vel(IDC,2)
          uz=vel(ICV,3)-vel(IDC,3)
          up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
          ux=ux-up*SFAREA(1,ICFL)
          uy=uy-up*SFAREA(2,ICFL)
          uz=uz-up*SFAREA(3,ICFL)
          up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
          yp=FRSTCV(IBFL)
          rnu=dum1*rmu(ICV)/rho(ICV)
          utaul=utau(IBFL)
!---------------
! --- log-law
!---------------
          yre=abs(up*yp/rnu)
          if(yre<1.d0) then
            utaul=dsqrt(rnu*up/yp)
          else
            ustar=vk*log(yre)+awfnc1
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
            uplus=uplus+duplus
            utaul=up/(uplus+SML)
          endif
          utau(IBFL)=utaul
          yplus=yp*utaul/(rnu+SML)
          if(yplus<yplsm) then
            utau(IBFL)=dsqrt(rnu*up/yp)
          endif
        enddo 
      elseif(kdv.eq.kvnslp) then
!        if(icaltb==sles.or.icaltb==dles.or.
!     %     icaltb==lles.and.icaltb==SDES) then    !huilai
!          call FFRABORT(1,'ERR: Set wall in [log-law] for LES')
!        endif
        do 2100 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ux=vel(ICV,1)-vel(IDC,1)
        uy=vel(ICV,2)-vel(IDC,2)
        uz=vel(ICV,3)-vel(IDC,3)
        up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
        ux=ux-up*SFAREA(1,ICFL)
        uy=uy-up*SFAREA(2,ICFL)
        uz=uz-up*SFAREA(3,ICFL)
        up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
        yp=FRSTCV(IBFL)
        rnu=rmu(ICV)/rho(ICV)
        utaul=dsqrt(rnu*up/yp)+SML
        utau(IBFL)=utaul
 2100   enddo
      endif
 100  enddo
!
! --- Wall Fraction Monitor
!
      if(icaltb.eq.noturb) return
!
!
!
      if(icaltb==dles.or.u_func(2)==2) then 
        call cal_Dynamic_SGS
     &   (LVEDGE,LBC_SSF,LCYCSF,LFUTAU,
     &  MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,DISALL,FRSTCV,GRDC,OPPANG,
     &  rho,vel,rmu,rmut,aks,tmp,utau,wifsld,LCYCOLD,vctr,FIELD_U)
!
        call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmut)
        if(u_func(3)==1) 
     &    call  cal_canopy(iter,deltt,time,
     &    MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,MAT_INDEX,
     &    rmu,rmut)
        if(u_func(2)/=2) return 
      endif
!
! --- < 1. Calculate dissipation function >- 2*SijSij
!
      call grad_cell(3,21,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vel,grdc)
!
!--- < 2. Calculate eddy viscosity >-
!
      if(icaltb.eq.sles.and.ical_vect) then  !-1
!CIDR NODEP
!        DO ICVL=ICVS_V,ICVE_V
        do IIMAT=1,NMAT 
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT<0) cycle
        DO ICVL=ICVS,ICVE
        IBFL=LFUTAU(ICVL)
        rnu=rmu(ICVL)/rho(ICVL)
        del=min(CVVOLM(ICVL)**r1pn,DISALL(ICVL))
        yplus=utau(IBFL)*DISALL(ICVL)/rnu
        dum1=1.d0
        if(yplus<10.3d0) dum1=0.d0
        rmut(ICVL)=dum1*rho(ICVL)*(csles*del    !CVVOLM(ICVL)**r1pn
     &       *(ONE-EXP(-(utau(LFUTAU(ICVL))*DISALL(ICVL)
     &             /(rmu(ICVL)/rho(ICVL)))*R26)))**2
     &  *dsqrt((2.d0*GRDC(ICVL,1,1)**2
     &         +2.d0*GRDC(ICVL,2,2)**2
     &         +2.d0*GRDC(ICVL,3,3)**2
     &         +((GRDC(ICVL,1,2)+GRDC(ICVL,2,1))**2
     &         + (GRDC(ICVL,1,3)+GRDC(ICVL,3,1))**2
     &         + (GRDC(ICVL,2,3)+GRDC(ICVL,3,2))**2)))
        enddo
        enddo
        call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmut)
        if(calgeq) then
          yplusf=0.0d0
          vflc  =0.0d0
          DO 432 ICVL=ICVS,ICVE
            SIJ11=grdc(ICVL,1,1)
            SIJ22=GRDC(ICVL,2,2)
            SIJ33=GRDC(ICVL,3,3)
            SIJ12=GRDC(ICVL,1,2)+GRDC(ICVL,2,1)
            SIJ13=GRDC(ICVL,1,3)+GRDC(ICVL,3,1)
            SIJ23=GRDC(ICVL,2,3)+GRDC(ICVL,3,2)
            dspd=(2.d0*SIJ11*SIJ11
     &           +2.d0*SIJ22*SIJ22
     &           +2.d0*SIJ33*SIJ33
     &           +(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
            IBFL=LFUTAU(ICVL)
            rnu=rmu(ICVL)/rho(ICVL)
            del=CVVOLM(ICVL)**r1pn
            yplus=utau(IBFL)*DISALL(ICVL)/rnu
            yplusf(ICVL)=yplus
            vflc(ICVL) =cfflc*del*dsqrt(dspd)
 432      enddo
        endif
!
        if(u_func(3)==1) call  cal_canopy(iter,deltt,time,
     &    MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,MAT_INDEX,
     &    rmu,rmut)
        return
      else
        do 410 IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 410
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT<0) cycle
        if(icaltb.eq.sles) then
! --- LES Smagorinsky model: dspd=(2*Sij*Sij); Sij=0.5(dui/dxj+duj/dxi)
          DO 430 ICVL=ICVS,ICVE
          SIJ11=grdc(ICVL,1,1)
          SIJ22=GRDC(ICVL,2,2)
          SIJ33=GRDC(ICVL,3,3)
          SIJ12=GRDC(ICVL,1,2)+GRDC(ICVL,2,1)
          SIJ13=GRDC(ICVL,1,3)+GRDC(ICVL,3,1)
          SIJ23=GRDC(ICVL,2,3)+GRDC(ICVL,3,2)
          dspd=(2.d0*SIJ11*SIJ11
     &         +2.d0*SIJ22*SIJ22
     &         +2.d0*SIJ33*SIJ33
     &         +(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
          IBFL=LFUTAU(ICVL)
          rnu=rmu(ICVL)/(rho(ICVL)+SML)
          del=min(CVVOLM(ICVL)**r1pn,DISALL(ICVL))
!          del=CVVOLM(ICVL)**r1pn
          yplus=utau(IBFL)*DISALL(ICVL)/(rnu+SML)
          
          fmu=(ONE-EXP(-yplus*R26))
          if(yplus<10.3d0) fmu=0.d0
!          if(yplus<8.d0) fmu=0.d0
          rmut(ICVL)=rho(ICVL)*(csles*del*fmu)**2*dsqrt(dspd)
 430      enddo
          if(calgeq) then
            yplusf=0.0d0
            vflc  =0.0d0
            DO ICVL=ICVS,ICVE
            SIJ11=grdc(ICVL,1,1)
            SIJ22=GRDC(ICVL,2,2)
            SIJ33=GRDC(ICVL,3,3)
            SIJ12=GRDC(ICVL,1,2)+GRDC(ICVL,2,1)
            SIJ13=GRDC(ICVL,1,3)+GRDC(ICVL,3,1)
            SIJ23=GRDC(ICVL,2,3)+GRDC(ICVL,3,2)
            dspd=(2.d0*SIJ11*SIJ11
     &           +2.d0*SIJ22*SIJ22
     &           +2.d0*SIJ33*SIJ33
     &           +(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
            IBFL=LFUTAU(ICVL)
            rnu=rmu(ICVL)/rho(ICVL)
            del=CVVOLM(ICVL)**r1pn
            yplus=utau(IBFL)*DISALL(ICVL)/rnu
            yplusf(ICVL)=yplus
            vflc(ICVL) =cfflc*del*dsqrt(dspd)
            enddo
          endif
        elseif(icaltb==KE2S.or.icaltb==ke.or.
     &         icaltb==RNG.or.icaltb==CHEN) then
! --- K-E two-layer model
          if(ical_EWT==1) then 
            DO ICVL=ICVS,ICVE 
            rans_k=aks(ICVL,ike(1)) 
            AAA=abs((0.2d0*Rey_s)/tanh(0.98d0))
            rnu=rmu(ICVL)/rho(ICVL) 
!Norris & Reynolds S------------------------
!            Rey=DISALL(ICVL)*sqrt(rans_k)/rnu  !ReY  
!            cl=akappa/cmu**0.75*DISALL(ICVL) 
!            fmu=(1.d0-exp(-Rey/Amu)) 
!            rmut(ICVL)=fmu*rho(ICVL)*cmu*rans_k**2 
!     &                /(aks(ICVL,ike(2))+SML) 
!Hassid etc S -----------------------------
!            lll=akappa*DISALL(ICVL)
!            Rey=lll*sqrt(rans_k)/rnu  !Rel 
!            fmu=(1.d0-exp(-Rey/Amu))  
!!            dum1=fmu*cmu**0.25*rho(ICVL)*lll*dsqrt(rans_k) 
!!            dum1=fmu*rho(ICVL)*cmu*rans_k**2 
!!     &                /(aks(ICVL,ike(2))+SML)
!            rmut(ICVL)=fmu*rho(ICVL)*cmu*rans_k**2 
!     &                /(aks(ICVL,ike(2))+SML)
!Wolfstein F----- ------------------------
            Rey=DISALL(ICVL)*sqrt(rans_k)/rnu  !ReY 
            fmu=(1.d0-exp(-Rey/Amu))  
            cl=akappa/cmu**0.75*DISALL(ICVL) 
            Lmu=Cl*(1.d0-exp(-Rey/Amu)) 
            dum1=rho(ICVL)*cmu*Lmu*dsqrt(rans_k) 
            dum2=min(1.d0,rho(ICVL)*cmu*rans_k**2   !6666 fmu 
     &                /(aks(ICVL,ike(2))+SML))
            rmd_E=0.5d0*(1.d0+TANH((Rey-Rey_s)/AAA)) 
            rmut(ICVL)=rmd_E*dum2+(1.d0-rmd_E)*dum1 
!-----------    ------------------------
            enddo
          else
            if(icaltb==KE2S) then 
              DO ICVL=ICVS,ICVE
              dum1=aks(ICVL,ike(1))
              dum2=aks(ICVL,ike(2))
              rnu=rmu(ICVL)/rho(ICVL) 
              Rey=sqrt(dum1)*L2s/rnu

              if(Rey>=89.d0) then
                dum3=-0.00261*Rey+0.5185d0
                if(dum3<0.d0) dum3=0.01d0
              elseif(Rey<4.d0) then
                dum3=1.23d0  !CDp
              else
                call fitVALUE(Rey,POLYCDP(:),mCD+1,dum3)
              endif

              CDx=dum3**(4.d0/3.d0)
              CDx=max(cmu,CDx)
              rmut(ICVL)=min(1.d0,rho(ICVL)*CDx*dum1**2
     &                /(dum2+SML))
              enddo
            else
              DO 440 ICVL=ICVS,ICVE
              rmut(ICVL)=min(1.d0,rho(ICVL)*cmu*aks(ICVL,ike(1))**2
     &                /(aks(ICVL,ike(2))+SML))
 440          enddo
            endif
          endif
        elseif(icaltb==ke_low) then
          DO ICVL=ICVS,ICVE
          rans_k=aks(ICVL,ike(1))
          rans_e=aks(ICVL,ike(2))
          rnu=rmu(ICVL)/rho(ICVL)
          Ry=dsqrt(rans_k)*DISALL(ICVL)/rnu 
          Rt=rans_k*rans_k*rho(ICVL)/(rmu(ICVL)*(rans_e)+SML1)
          fmu=(ONE-exp(-0.0198d0*Ry))*(1.d0+5.29d0/Ry)
          rmut(ICVL)=fmu*rho(ICVL)*cmu*rans_k**2/(rans_e+SML)
          enddo
        elseif(icaltb==SDES) then
          DO ICVL=ICVS,ICVE
            rans_k=aks(ICVL,ides)
            xxx=rans_k/rmu(ICVL)*rho(ICVL)
            fv1des=xxx**3/(xxx**3+cv1**3)
            rmut(ICVL)=rho(ICVL)*fv1des*rans_k
          enddo
        elseif(icaltb==KLES) then 
          if(u_func(2)==0) then
            DO ICVL=ICVS,ICVE
            IBFL=LFUTAU(ICVL)
            rans_k=dsqrt(aks(ICVL,ikles))
            rnu=rmu(ICVL)/(rho(ICVL)+SML)
            del=min(CVVOLM(ICVL)**r1pn,DISALL(ICVL))
            yplus=utau(IBFL)*DISALL(ICVL)/(rnu+SML)
            fmu=(ONE-EXP(-yplus*R26))
            rmut(ICVL)=rho(ICVL)*rans_k*del
            enddo
          else  !u_func(2)=1 or 2
            DO ICVL=ICVS,ICVE
            if(multiR(ICVL)==1) then  !
              IBFL=LFUTAU(ICVL)
              rans_k=dsqrt(aks(ICVL,ikles))
              rnu=rmu(ICVL)/(rho(ICVL)+SML)
              del=min(CVVOLM(ICVL)**r1pn,DISALL(ICVL))
              yplus=utau(IBFL)*DISALL(ICVL)/(rnu+SML)
              fmu=(ONE-EXP(-yplus*R26))
              rmut(ICVL)=rho(ICVL)*rans_k*del
            endif
            enddo
          endif
          if(u_func(2)==1) then
            DO ICVL=ICVS,ICVE
            if(multiR(ICVL)==2) then  !
              SIJ11=grdc(ICVL,1,1)
              SIJ22=GRDC(ICVL,2,2)
              SIJ33=GRDC(ICVL,3,3)
              SIJ12=GRDC(ICVL,1,2)+GRDC(ICVL,2,1)
              SIJ13=GRDC(ICVL,1,3)+GRDC(ICVL,3,1)
              SIJ23=GRDC(ICVL,2,3)+GRDC(ICVL,3,2)
              dspd=(2.d0*SIJ11*SIJ11
     &             +2.d0*SIJ22*SIJ22
     &             +2.d0*SIJ33*SIJ33
     &             +(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
              IBFL=LFUTAU(ICVL)
              rnu=rmu(ICVL)/(rho(ICVL)+SML)
              del=min(CVVOLM(ICVL)**r1pn,DISALL(ICVL))
              yplus=utau(IBFL)*DISALL(ICVL)/(rnu+SML)
              fmu=(ONE-EXP(-yplus*R26))
              if(yplus<10.3d0) fmu=0.d0
              rmut(ICVL)=rho(ICVL)*(csles*del*fmu)**2*dsqrt(dspd)
            endif
            enddo
          endif
        elseif(icaltb==cnst) then
          DO ICVL=ICVS,ICVE
            rmut(ICVL)=turbmu(IMAT)
          enddo
          if(ieul2ph>0) then
            rmut2(ICVL)=turbmu2(IMAT)
          endif
        elseif(icaltb==RSM) then
        endif
 410    enddo
      endif

!
! --- Dummy cell for BC:
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmut)
!
      if(u_func(3)==1) 
     &    call  cal_canopy(iter,deltt,time,
     &    MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,MAT_INDEX,
     &    rmu,rmut)
      return
      end subroutine cal_eddy_viscosity
!
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_rhoccc
     &           (mph,NCVX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &            hhh,hhs,cps,cp,cr,
     &            aks,tmp,yys,pp0,prs,rho,rmu,rmd,RHO2,VEL,ccc,
     &            yysbnd,prsbnd,tmpbnd,LBC_SSF,lvedge,imode
     &            )
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     ccc: (sound-speed)^2
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! zhang-cvd prs=>pp0
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_initial 
      use module_material,only : rsld,ical_sld,relaxC
      use module_species,only  : sw,gascns,wm,r_wm,I_H2OL,act
      use module_model,   only : ical_topt,ical_t,mach0,idrdp,comp,
     &                           ical_week,ical_dens,incomp,PEFC,
     &                           ical_dens
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_gravity ,only : beta,tamb,ramb,pamb,ysamb,iave_rho_T
      use module_scalar,only   : ical_cavi,icavi,Pres_c,
     &                           P_str_v,Tref,T0_ref,KL,Rv,ivold,
     &                           ical_FC,ical_s
      use module_initial ,only : rho0,rho02
      use module_metrix  ,only : W1K9
      use module_boundary,only : chkncomp,chknrans,nobcnd,
     &                           wdbcnd,adbcnd,iscomp,israns,distrb,
     &                           LBC_INDEX,MAT_BCIDX,nbcnd,kdbcnd,
     &                           kdintr,kdilet,kdolet,
     &                           kxnone,kdbuff,kdintr,kdfire,rotwall,
     &                           kdcvd,sizein,kdshutr,
     &                           endw,beginw,omega,
     &                           kdpres,openout,kdstag,
     &                           masflg,masflx,imasflg,sumare,
     &                           vofang,shut_strt,ktEWF,kyEWF
!
! 1. Calculate density & speed of sound
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: NCVX,mph,imode
      real*8 ,intent(in)    :: pp0     (  MXALLCV)
      real*8 ,intent(in)    :: prs     (  MXALLCV,2)
      real*8 ,intent(out)   :: rho     (  MXALLCV,2)
      real*8 ,intent(out)   :: rho2    (  MXALLCVC)
      real*8 ,intent(out)   :: ccc     (  MXALLCV)
      real*8 ,intent(out)   :: tmp     (  MXALLCV)
      real*8 ,intent(in)    :: yys     (  MXALLCV,mxcomp)
      real*8 ,intent(inout) :: aks     (  MXALLCVR,MXRANS)
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   ( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  ( 0:MXMAT)
! --- temporary array 
      REAL*8 ,INTENT(INOUT) :: hhh   (    MXALLCV)
      REAL*8 ,INTENT(INOUT) :: hhs   (    MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cps   (    MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp    (    MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cr    (    MXALLCV)
      real*8 ,intent(in)    :: vel   (    MXALLCV,3)
      real*8 ,intent(in)    :: RMD   (        MXALLCV)
      real*8 ,intent(in)    :: RMU   (        MXALLCV)
      REAL*8 ,INTENT(IN) :: PRSBND(        MXSSFBC)
      REAL*8 ,INTENT(IN) :: TMPBND(        MXSSFBC)
      REAL*8 ,INTENT(IN) :: YYSBND(        MXSSFBC,MXCOMP)
      integer,intent(in)   :: LBC_SSF( MXSSFBC)
      integer,intent(in)   :: LVEDGE(2,MXCVFAC)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,myid
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ICOM
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: ibfs,ibfe,ibfl
      REAL*8  :: pressure,crl,dum1,dum2,dum3,dum4
      REAL*8,parameter :: gamm=6.722d-5,const=0.615d0
      integer,SAVE :: ICAL=0
!
! --- 
!
      if(idrdp==incomp) then 
        do IIMAT=1,NMAT
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        IDCS = MAT_DCIDX(IIMAT-1)+1
        IDCE = MAT_DCIDX(IIMAT)
        if(IMAT>0) then 
          if(ical_dens==1) then
            call cal_fluid_1(IIMAT,ICVS,ICVE)
            if(NCVX>NCV) call cal_fluid_1(IIMAT,IDCS,IDCE)
          elseif(ical_dens==2) then
            call cal_fluid_2(IIMAT,ICVS,ICVE)
            if(NCVX>NCV) call cal_fluid_2(IIMAT,IDCS,IDCE)
          elseif(ical_dens==3) then
            call cal_fluid_3(IIMAT,ICVS,ICVE)
            if(NCVX>NCV) call cal_fluid_3(IIMAT,IDCS,IDCE)
          elseif(ical_dens==5) then
            call cal_fluid_5(IIMAT,IDCS,IDCE)
            if(NCVX>NCV) call cal_fluid_5(IIMAT,IDCS,IDCE)
          endif
          if(ical_cavi==1) then
            call cal_fluid_cavi(IIMAT,ICVS,ICVE,1)
            if(NCVX>NCV) call cal_fluid_cavi(IIMAT,IDCS,IDCE,2)
          endif
        endif
        enddo
        if(.not.ical_week) return
      endif
!
      hhh=0.d0
      hhs(1:MXALLCV,1:MXcomp)=0.d0
      cps(1:MXALLCV,1:MXcomp)=0.d0
      cp=0.d0
      cr=0.d0
!
      if(sw) then
        if(ical_dens==4) then !zhang5
        else
          call cal_t2cp(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
!     &                tmp,yys,prs,rho,rmu,rmd,cps,cp,cr)
     &                tmp,yys,pp0,rho,rmu,rmd,cps,cp,cr)
        endif
      else
        call cal_t2hcp(mph,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
!     &                 tmp,yys,hhs,cps,hhh,cp,cr,prs)
     &                 tmp,yys,hhs,cps,hhh,cp,cr,pp0)
      endif
!
      if(idrdp==comp.and.ical_dens==4) then 
      endif
!
!
      if(ical_week) then
        if(ICAL_VECT) then   !NOVECT
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V   !index_c(myid)+1,index_c(myid+1)
          ccc(ICVL)=relaxC(1)**2*cp(ICVL)
     &             /(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
          enddo
!CIDR NODEP
          DO ICVL=IDCS_V,IDCE_V!index_c(myid)+1,index_c(myid+1)
          ccc(ICVL)=relaxC(1)**2*cp(ICVL)
     &             /(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
          enddo
        else 
          do IIMAT=1,NMAT
          IMAT = MAT_NO(IIMAT)
          ICVS = MAT_CVEXT(IIMAT-1)+1
          ICVE = MAT_CVEXT(IIMAT)
          IDCS = MAT_DCIDX(IIMAT-1)+1
          IDCE = MAT_DCIDX(IIMAT)
          do ICVL=ICVS,ICVE
          ccc(ICVL)=relaxC(IIMAT)**2*
     &             cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
          enddo
          do ICVL=ICVS,ICVE
          ccc(ICVL)=relaxC(IIMAT)**2*
     &             cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
          enddo
          
          do ICVL=IDCS,IDCE
          ccc(ICVL)=relaxC(IIMAT)**2*
     &             cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
          enddo
          enddo
        endif
        return
      endif
! --- heat insulation 
      if(ical_topt==2) then
        callFFRABORT(1,'ERR: 8888')
        do IIMAT=1,NMAT
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        IDCS = MAT_DCIDX(IIMAT-1)+1
        IDCE = MAT_DCIDX(IIMAT)
        do ICVL=ICVS,ICVE
        pressure=max(pp0(ICVL),1.d0)
        rho(ICVL,1)=gamm*pressure**const
        ccc(ICVL)=relaxC(IIMAT)**2*
     &             cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
        enddo
        do ICVL=IDCS,IDCE
        pressure=max(pp0(ICVL),1.d0)
        rho(ICVL,1)=gamm*pressure**const
        ccc(ICVL)=relaxC(IIMAT)**2*
     &             cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
        enddo
        enddo
        return
      endif
!
!------------------------------------------------
! --- Calculate specific heat from temperature
!------------------------------------------------
!
      do IIMAT=1,NMAT   !ICV=1,NCVX ???
      if( .not.mat_cal(IIMAT) ) cycle
      IMAT = MAT_NO(IIMAT)
      ICVS = MAT_CVEXT(IIMAT-1)+1
      ICVE = MAT_CVEXT(IIMAT)
      IDCS = MAT_DCIDX(IIMAT-1)+1
      IDCE = MAT_DCIDX(IIMAT)
      if(IMAT>0) then ! fluid
        if(ical_cavi==1) then
          call cal_fluid_cavi(IIMAT,ICVS,ICVE,1)
          call cal_fluid_cavi(IIMAT,IDCS,IDCE,2)
        elseif(ical_dens==5) then
          call cal_fluid_5(IIMAT,ICVS,ICVE)
          if(NCVX>NCV) call cal_fluid_5(IIMAT,IDCS,IDCE)
        else
          if(ICAL_VECT) then  !NOVECT
            call cal_fluid_VR(ICVS,ICVE)
            if(NCVX>NCV) call cal_fluid_VR(IDCS,IDCE)
          else
            call cal_fluid(ICVS,ICVE,1)
            if(NCVX>NCV) call cal_fluid(IDCS,IDCE,2)
!            call cal_fluid(IDCS,IDCE,2)

            if(imasflg) then
              do nb=1,nbcnd
              kd=kdbcnd(0,nb)
              if(kd==kdolet.and.(openout(nb)==8)) then
                IBFS=LBC_INDEX(nb-1)+1
                IBFE=LBC_INDEX(nb)

                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                IDC=LVEDGE(2,ICFL)
                dum2=0.d0
!
                do ICOM=1,ncomp 
                dum2=dum2+yysbnd(IBFL,ICOM)*r_wm(ICOM)
                enddo
!
                dum3=prsbnd(IBFL)/(gascns*dum2*tmpbnd(IBFL)) 
                rho(IDC,1)=dum3
                enddo

              endif
              enddo
            endif

          endif
        endif
      else            ! solid
        do ICVL=ICVS,ICVE
          rho(ICVL,1) = rsld(-IMAT)
          ccc(ICVL) = 1.d6
        enddo
      endif
      enddo
!
      return
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!=================================================
      subroutine cal_fluid(do_start,do_end,imode)
!=================================================
      integer,intent(in) :: do_start,do_end,imode
      REAL*8  :: dum1
!
      do ICVL=do_start,do_end 
        crl=0.d0 
        dum2=(cp(ICVL)-cr(ICVL)) 
        do ICOM=1,ncomp 
        crl=crl+yys(ICVL,ICOM)*r_wm(ICOM)*act(ICOM)
        enddo
! [m^2/s^2]
! [kg/m^3]
! 
! ---------------------- gascns*crl=cr(ICVL)
        rho(ICVL,1)=pp0(ICVL)/(gascns*crl*tmp(ICVL))
        ccc(ICVL)=relaxC(IIMAT)**2*
     &             cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL) 
!        ccc(ICVL)=gascns*crl*tmp(ICVL) 
! ---------------------- 
!        rho(ICVL,1)=prs(ICVL)/(gascns*crl*tmp(ICVL)) 
!        ccc(ICVL)=cp(ICVL)/(gascns*gascns*crl*tmp(ICVL)) 
! ---------------------- 
!        ccc(ICVL) = gascns*tmp(ICVL)*1.4d0 
! !       ccc(ICVL) = cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL) 
                      ! [m^2/s^2] 
!        ccc(ICVL) = cp(ICVL)/(gascns*cr(ICVL)*tmp(ICVL)) 
!        rho(ICVL,1)=1.4d0*prs(ICVL)/tmp(ICVL)+1.d0/tmp(ICVL) 
!
!!        rho(ICVL,1)=prs(ICVL)/(gascns*crl*tmp(ICVL))
!        rho(ICVL,1)=1.4d0*ccc(ICVL)*prs(ICVL)/tmp(ICVL) 
!        ccc(ICVL)=gascns*tmp(ICVL)*1.4d0 
!        ccc(ICVL)=prs(ICVL)*1.4d0/rho(ICVL,1) 
!        
!prs(ICVL)*1.4d0/rho(ICVL)      !gascns*tmp(ICVL)*1.4d0 
      enddo 
!
!
!

      IF(ICAL==0.AND.iave_rho_T==4) THEN 
        crl=0.d0
        do ICOM=1,ncomp
        crl=crl+ysamb(ICOM)*r_wm(ICOM)
        enddo
        ramb=pamb/(gascns*crl*tamb)
        ICAL=1
      ENDIF
!
      end subroutine cal_fluid
!
!=================================================
      subroutine cal_fluid_VR(do_start,do_end)
!=================================================
      integer,intent(in) :: do_start,do_end
!
!
!      W1K9(do_start:do_end)=0.d0
!      do ICOM=1,ncomp
!      do ICVL=do_start,do_end
!      W1K9(ICVL)=W1K9(ICVL)+yys(ICVL,ICOM)*r_wm(ICOM)
!      enddo
!      enddo
!
      do ICVL=do_start,do_end
!      rho(ICVL,1)=pp0(ICVL)/(gascns*W1K9(ICVL)*tmp(ICVL))
      rho(ICVL,1)=pp0(ICVL)/(cr(ICVL)*tmp(ICVL))
      ccc(ICVL) = relaxC(IIMAT)**2*
     &             cp(ICVL)/(cp(ICVL)-cr(ICVL))*cr(ICVL)*tmp(ICVL)
      enddo
! 
      IF(ICAL==0.AND.iave_rho_T==4) THEN
        crl=0.d0
        do ICOM=1,ncomp
        crl=crl+ysamb(ICOM)*r_wm(ICOM)
        enddo
        ramb=pamb/(gascns*crl*tamb)
        ICAL=1
      ENDIF
!
      end subroutine cal_fluid_VR

!=====================================================
      subroutine cal_fluid_1(IIMAT,do_start,do_end)
!=====================================================
      integer,intent(in) :: IIMAT,do_start,do_end
!
      if(ncomp==1) return
      do ICVL=do_start,do_end
        crl=0.d0
        do ICOM=1,ncomp 
        dum1=p0(IIMAT)/(gascns*r_wm(ICOM)*tmp(ICVL))
!       dum1=pp0(ICVL)/(gascns*r_wm(ICOM)*tmp(ICVL))

        crl=crl+yys(ICVL,ICOM)/dum1
        enddo
        rho(ICVL,1)=1.d0/crl
                      ! [kg/m^3]
      enddo
      end subroutine cal_fluid_1
!
!===================================================
      subroutine cal_fluid_2(IIMAT,do_start,do_end)
!===================================================
      integer,intent(in) :: IIMAT,do_start,do_end
!
      if(abs(rho0(IIMAT)-ramb)>1.d-10)
     &     call FFRABORT(1,'ERR: Initial [dens] /= rho')
      do ICVL=do_start,do_end
        crl=1.d0+beta*(tmp(ICVL)-tamb)
        rho(ICVL,1)=rho0(IIMAT)/crl
                      !  [kg/m^3] 
      enddo
      end subroutine cal_fluid_2
!
!===================================================
      subroutine cal_fluid_3(IIMAT,do_start,do_end)
!===================================================
      integer,intent(in) :: IIMAT,do_start,do_end
!
      if(ncomp==1) return
      do ICVL=do_start,do_end
        crl=0.d0
        do ICOM=1,ncomp 
          dum1=yys(ICVL,ICOM)/(wm(ICOM)-1.41d0*tmp(ICVL))
          crl=crl+dum1
        enddo
        rho(ICVL,1)=1.d0/crl
                      ! [kg/m^3] 
      enddo
      end subroutine cal_fluid_3
!
!=============================================================
      subroutine cal_fluid_cavi(IIMAT,do_start,do_end,IMODE)
!=============================================================
      integer,intent(in) :: IIMAT,do_start,do_end,IMODE
      real*8 :: Pv,TTT
      real*8,external :: Pvap
!
      do ICVL=do_start,do_end
        TTT=tmp(ICVL)
        Pv=Pvap(TTT)
!--------------------------------------(1)
        dum1=aks(ICVL,icavi)
        crl=dum1*rho0(IIMAT)/
     &     (dum1*rho0(IIMAT)+(1.d0-dum1)*rho02(IIMAT))
        dum1=dmin1(1.d0,dmax1(0.d0,crl))
!        aks(ICVL,ivold)=dum1
!--------------------------------------(2)
        dum3=max(prs(ICVL,1),Pv)
        dum1=KL*(1.d0-aks(ICVL,icavi))*dum3*(tmp(ICVL)+T0_ref)
        dum2=Rv*aks(ICVL,icavi)*(dum3+Pres_c)*tmp(ICVL)
        crl=1.d0/(dum1+dum2)
        rho(ICVL,1)=dum3*(dum3+Pres_c)*crl  !4444
                      ! [kg/m^3] 
        crl=(rho(ICVL,1)-(prs(ICVL,1)+Pres_c)/KL/(tmp(ICVL)+T0_ref))
     &  /(prs(ICVL,1)/Rv/tmp(ICVL) 
     &   -(prs(ICVL,1)+Pres_c)/KL/(tmp(ICVL)+T0_ref)) 
!        aks(ICVL,ivold)=dmin1(1.d0,dmax1(0.d0,crl))   !4444
!---------------- 
        crl=dum2/(dum1+dum2)
        aks(ICVL,ivold)=dmin1(1.d0,dmax1(0.d0,crl))   !4444
!-------------------------------------- 
! --- 
!
        dum2=0.d0
        do ICOM=1,ncomp
        dum2=dum2+yys(ICVL,ICOM)*r_wm(ICOM)
        enddo
!        dum2=max(Pv,pp0(ICVL))/(gascns*dum2*TTT)
        dum2=Pv/(gascns*dum2*TTT)
!        dum3=rho0(IIMAT)*(1.d0-dum1)+dum2*dum1

!        rho(ICVL,1)=dum3
        rho2(ICVL)=dum2

!--------------------------------------
        dum3=1.222d-3*Pres_c             !water: c^2
        dum4=4.d0/3.d0*461.52d0*TTT      !steam: c^2
        dum1=aks(ICVL,ivold)
        ccc(ICVL)=relaxC(IIMAT)**2*
     &             (dum3*(1.d0-dum1)+dum4*dum1)
!--------------------------------------
      enddo
      end subroutine cal_fluid_cavi
!
!===================================================
      subroutine cal_fluid_5(IIMAT,do_start,do_end)
!===================================================
      integer,intent(in) :: IIMAT,do_start,do_end
      real*8,save :: water_s=1234.d0,gas_s=340.d0
      real*8 :: yys_H2OL,dum5
      
!
      if(ical_FC/=PEFC) 
     &  call FFRABORT(1,'ERR: density =5=>Fuel_cell=/PEFC') 
      do ICVL=do_start,do_end 
      dum1=aks(ICVL,ical_s)
      dum2=1.d0-aks(ICVL,ical_s)
      crl=0.d0
      yys_H2OL=yys(ICVL,I_H2OL)
      do ICOM=1,ncomp 
      if(ICOM==I_H2OL) cycle
      crl=crl+r_wm(ICOM)*yys(ICVL,ICOM)
      enddo
      dum5=yys(ICVL,I_H2OL)
      crl=crl*gascns*tmp(ICVL)/prs(ICVL,1)+yys_H2OL/rho02(IIMAT)  !???
      dum4=1.d0/crl
      crl=dum1*rho02(IIMAT)+dum2*dum4
      rho(ICVL,1)=crl
      rho2(ICVL)=dum4
      ccc(ICVL)=relaxC(IIMAT)**2*
     &             (1.d0/((1.d0-yys_H2OL)/gas_s+yys_H2OL/water_s)*dum2
     &         +dum1*water_s)
      enddo
      end subroutine cal_fluid_5
!
      end subroutine cal_rhoccc
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_rva(
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,wiface,SFCENT,
     &  LCYCOLD,wifsld,OPPANG,
     &  xta,rho,vel,rvx,rva)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. Calculate mass flux at cell face
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_Euler2ph,only : ieul2ph
!
      implicit none
! --- [dummy arguments]
!
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT (3,MXCVFAC)
      real*8 ,intent(in)    :: xta    (  MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(inout) :: rho    (  MXALLCV)
      real*8 ,intent(in)    :: vel    (  MXALLCV,3)
      real*8 ,intent(out)   :: rva    (  MXCVFAC)
      REAL*8 ,INTENT(INOUT) :: rvx    (  MXALLCV,3)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
!
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,IV,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD
      integer :: IMAT_U
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      real*8  :: ru,rv,rw,wi1,wi2
!
! --- 
!
      call dc_symprv
     &(1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & LCYCOLD,wifsld,OPPANG,
     & SFAREA,SFCENT,vel,0,1)
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rho)
!
! --- 
!
!
! --- Face velocity is averaged from cell center.
!
      if(ical_vect) then  !NOVECT
        do iv=1,3
        rvx(ICVS_V:ICVE_V,iv)=
     &  rho(ICVS_V:ICVE_V)*vel(ICVS_V:ICVE_V,iv)
        rvx(IDCS_V:IDCE_V,iv)=
     &  rho(IDCS_V:IDCE_V)*vel(IDCS_V:IDCE_V,iv)
        enddo
      else
!        if(ieul2ph==0) then
        do 100 IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 100
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        do iv=1,3
        rvx(ICVS:ICVE,iv)=rho(ICVS:ICVE)*vel(ICVS:ICVE,iv)
        rvx(IDCS:IDCE,iv)=rho(IDCS:IDCE)*vel(IDCS:IDCE,iv)
        enddo
  100   enddo
!        endif
      endif
!
      if(ical_vect) then  !NOVECT
!CIDR NODEP
        do ICFL=ICFS_V,ICFE_V
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        ru=wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1)
        rv=wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2)
        rw=wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3)
        rva(ICFL)=(rho(ICVA)*max(0.d0,xta(ICFL))
     &            +rho(ICVB)*min(0.d0,xta(ICFL)))
     &         +SFAREA(4,ICFL)*
     &   (ru*SFAREA(1,ICFL)+rv*SFAREA(2,ICFL)+rw*SFAREA(3,ICFL))
        enddo
      else
        if(ieul2ph>0) then
          do 200 IIMAT=1,NMAT     !ICF=1,NCVFAC
          if(.not.mat_cal(IIMAT)) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          ru=wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1)
          rv=wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2)
          rw=wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3)
          rva(ICFL)=(rho(ICVA)*max(0.d0,xta(ICFL))
     &              +rho(ICVB)*min(0.d0,xta(ICFL)))
     &         +SFAREA(4,ICFL)*
     &     (ru*SFAREA(1,ICFL)+rv*SFAREA(2,ICFL)+rw*SFAREA(3,ICFL))
          enddo
  200     enddo
        else
          do IIMAT=1,NMAT     !ICF=1,NCVFAC
          if(.not.mat_cal(IIMAT)) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          ru=wi1*rvx(ICVA,1)+wi2*rvx(ICVB,1)
          rv=wi1*rvx(ICVA,2)+wi2*rvx(ICVB,2)
          rw=wi1*rvx(ICVA,3)+wi2*rvx(ICVB,3)
          rva(ICFL)=(rho(ICVA)*max(0.d0,xta(ICFL))
     &              +rho(ICVB)*min(0.d0,xta(ICFL)))
     &         +SFAREA(4,ICFL)*
     &     (ru*SFAREA(1,ICFL)+rv*SFAREA(2,ICFL)+rw*SFAREA(3,ICFL))
          enddo
          enddo
        endif
      endif
!
      call debug
      return
!
!/////////////////////////////////////////////////////////////////////
      contains
!
!=================================================
      subroutine debug
!=================================================
      use module_debug,only : idebug
      if( idebug(11).eq.0 ) return
      end subroutine debug
!
      end subroutine cal_rva
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_t2hcp(mph,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,hhs,cps,hhh,cp,cr,prs)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_species,only  : gascns,acpk,wm,Ri,acpk2,act
      use module_material,only : cpsld
      use module_model,only    : ical_vect,nthrds,ical_dens
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_metrix  ,only : ahk,acpk_w
      use module_scalar,only   : Rv
!
! 1. Calculate enthalpy & specific heat from temperature
!
      implicit none
! --- [dummy arguments]
!
      INTEGER,INTENT(IN)  :: mph
      real*8 ,intent(in)  :: tmp    (   MXALLCV)
      real*8 ,intent(in)  :: yys(       MXALLCV,mxcomp)
!
      INTEGER,INTENT(IN)  :: MAT_CVEXT (0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_DCIDX (0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO  (  0:MXMAT)
      logical,INTENT(IN)  :: mat_cal (  0:MXMAT)
      real*8 ,intent(out) :: hhs(       MXALLCV,mxcomp)
      real*8 ,intent(out) :: cps(       MXALLCV,mxcomp)
      real*8 ,intent(out) :: hhh    (   MXALLCV)
      real*8 ,intent(out) :: cp     (   MXALLCV)
      real*8 ,intent(out) :: cr     (   MXALLCV)
      real*8 ,intent(in)    :: prs     (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: tt,hhx,cpx,crx
      integer :: i,j,k,l,m,n,ierr1=0,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD
      integer :: IMAT_U
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
!
! --- 
!
      if(mph==1) then
        do 100 ICOM=1,ncomp
        ahk(0,ICOM)=acpk(0,ICOM)
        do 105 n=1,5
        ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
  105   enddo
  100   enddo
        acpk_w=acpk
      else
        do ICOM=1,ncomp
        ahk(0,ICOM)=acpk2(0,ICOM)
        do n=1,5
        ahk(n,ICOM)=acpk2(n,ICOM)/dble(n)
        enddo
        enddo
        acpk_w=acpk2
      endif
!
      if(ical_vect) then  !NOVECT
        do ICOM=1,ncomp
!CIDR NODEP
          DO ICVL=1,NALLCV
          hhs(ICVL,ICOM)=((((
     &      ahk(5,ICOM) *tmp(ICVL)
     &     +ahk(4,ICOM))*tmp(ICVL)
     &     +ahk(3,ICOM))*tmp(ICVL)
     &     +ahk(2,ICOM))*tmp(ICVL)
     &     +ahk(1,ICOM))*tmp(ICVL)
     &     +ahk(0,ICOM)
          cps(ICVL,ICOM)=(((
     &      acpk_w(5,ICOM) *tmp(ICVL)
     &     +acpk_w(4,ICOM))*tmp(ICVL)
     &     +acpk_w(3,ICOM))*tmp(ICVL)
     &     +acpk_w(2,ICOM))*tmp(ICVL)
     &     +acpk_w(1,ICOM)
          ENDDO
        enddo
!
!CIDR NODEP
        DO ICVL=1,NALLCV   !index_c(myid)+1,index_c(myid+1)
          hhh(ICVL)=0.d0
          cp(ICVL)=0.d0
          cr(ICVL)=0.d0
        enddo
!
!CIDR NODEP
        do ICOM=1,ncomp
          DO ICVL=1,NALLCV    !index_c(myid)+1,index_c(myid+1)
          hhh(ICVL)=hhh(ICVL)+yys(ICVL,ICOM)*hhs(ICVL,ICOM)
          cp(ICVL)=cp(ICVL)+yys(ICVL,ICOM)*cps(ICVL,ICOM)
          cr(ICVL)=cr(ICVL)+yys(ICVL,ICOM)*Ri(ICOM)
          enddo
        enddo
      else
        do IIMAT=1,NMAT    !ICV=1,NCVX
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
!
        if(IMAT>0) then
          if(ical_dens==4)then
            call FFRABORT(1,'ical_dens==4: Real Gas /=cal_t2hcp')
          else
            call cal_fluid(ICVS,ICVE)
            call cal_fluid(IDCS,IDCE)
          endif
        else
! --- solid:
          call cal_solid(ICVS,ICVE)
! --- solid dummy CV
!        if(NCVX.gt.NCV) then
          call cal_solid(IDCS,IDCE)
        endif
        enddo
      endif
!
!      DEALLOCATE(ahk)
      return
!//////////////////////////////////////////////////////////////////
      contains
!
!===============================================
      subroutine cal_fluid(do_start,do_end)
!===============================================
      integer,intent(in) :: do_start,do_end
      do ICVL=do_start,do_end
        tt=tmp(ICVL)
        hhx=0.d0
        cpx=0.d0
        crx=0.d0
        do ICOM=1,ncomp
            hhs(ICVL,ICOM)=((((
     &      ahk(5,ICOM) *tt
     &     +ahk(4,ICOM))*tt
     &     +ahk(3,ICOM))*tt
     &     +ahk(2,ICOM))*tt
     &     +ahk(1,ICOM))*tt
     &     +ahk(0,ICOM)
          cps(ICVL,ICOM)=(((
     &      acpk_w(5,ICOM) *tt
     &     +acpk_w(4,ICOM))*tt
     &     +acpk_w(3,ICOM))*tt
     &     +acpk_w(2,ICOM))*tt
     &     +acpk_w(1,ICOM)
          hhx=hhx+yys(ICVL,ICOM)*hhs(ICVL,ICOM)*act(icom)
          cpx=cpx+yys(ICVL,ICOM)*cps(ICVL,ICOM)*act(icom)
          crx=crx+yys(ICVL,ICOM)*Ri(ICOM)*act(icom)
        enddo
        hhh(ICVL)=hhx
        cp(ICVL)=cpx
        cr(ICVL)=crx
      enddo
      end subroutine cal_fluid
!
!===============================================
      subroutine cal_solid(do_start,do_end)
!===============================================
      use module_material,only : cpsld
      integer,intent(in) :: do_start,do_end
!
      do ICVL=do_start,do_end
        cr(ICVL)  = 1.d0
        cp(ICVL)  = cpsld(-IMAT)
        hhh(ICVL) = cpsld(-IMAT)*tmp(ICVL)
        do ICOM=1,ncomp
          cps(ICVL,ICOM) = cp(ICVL)
          hhs(ICVL,ICOM) = hhh(ICVL)
        enddo
      enddo
      end subroutine cal_solid
!
      end subroutine cal_t2hcp
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_t2hcp_dc(mph,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,hhh)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_species,only  : gascns,acpk,wm,Ri,acpk2
      use module_material,only : cpsld
      use module_model,only    : ical_vect,nthrds
      use module_metrix  ,only : ahk,acpk_w
!
! 1. Calculate enthalpy & specific heat from temperature
!
      implicit none
! --- [dummy arguments]
!
      INTEGER,INTENT(IN)  :: mph
      real*8 ,intent(in)  :: tmp    (   MXALLCV)
      real*8 ,intent(in)  :: yys(       MXALLCV,mxcomp)
!
      INTEGER,INTENT(IN)  :: MAT_DCIDX (0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO  (  0:MXMAT)
      logical,INTENT(IN)  :: mat_cal (  0:MXMAT)
      real*8 ,intent(out) :: hhh    (   MXALLCV)
!
! --- [local entities]
!
      real*8  :: tt,hhx,cpx,crx,hhs
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL
      integer :: IMODE,ICOM,IMD,ierr1=0
      integer :: IMAT_U
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
!
! --- 
!
      if(mph==1) then
        do 100 ICOM=1,ncomp
        ahk(0,ICOM)=acpk(0,ICOM)
        do 105 n=1,5
        ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
  105   enddo
  100   enddo
        acpk_w=acpk
      else
        do ICOM=1,ncomp
        ahk(0,ICOM)=acpk2(0,ICOM)
        do n=1,5
        ahk(n,ICOM)=acpk2(n,ICOM)/dble(n)
        enddo
        enddo
        acpk_w=acpk
      endif
!
      do IIMAT=1,NMAT    !ICV=1,NCVX
        if( .not.mat_cal(IIMAT) ) cycle
        IMAT=MAT_NO(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
!
        if(IMAT>0) then
          if(ical_vect) then  !NOVECT
            call cal_fluid_V(IDCS,IDCE)
          else
            call cal_fluid(IDCS,IDCE)
          endif
        else
! --- solid dummy CV
          call cal_solid(IDCS,IDCE)
        endif
      enddo
!
      return
!//////////////////////////////////////////////////////////////////
      contains
!
!===============================================
      subroutine cal_fluid(do_start,do_end)
!===============================================
      integer,intent(in) :: do_start,do_end
!
!      IF(ical_vect) call FFRABORT(1,'ERR:cal_t2hcp_dc/cal_fluid')
      do ICVL=do_start,do_end
        tt=tmp(ICVL)
        hhx=0.d0
        do ICOM=1,ncomp
          hhs=((((ahk(5,ICOM)*tt+ahk(4,ICOM))*tt+ahk(3,ICOM))*tt
     &     +ahk(2,ICOM))*tt+ahk(1,ICOM))*tt+ahk(0,ICOM)
          hhx=hhx+yys(ICVL,ICOM)*hhs
        enddo
        hhh(ICVL)=hhx
      enddo
      end subroutine cal_fluid
!
!===============================================
      subroutine cal_fluid_V(do_start,do_end)
!===============================================
      integer,intent(in) :: do_start,do_end
!
      hhh(do_start:do_end)=0.d0
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      tt=tmp(ICVL)
      hhs=((((ahk(5,ICOM)*tt+ahk(4,ICOM))*tt+ahk(3,ICOM))*tt
     &     +ahk(2,ICOM))*tt+ahk(1,ICOM))*tt+ahk(0,ICOM)
      hhh(ICVL)=hhh(ICVL)+yys(ICVL,ICOM)*hhs
      enddo
      enddo
      end subroutine cal_fluid_V
!===============================================
      subroutine cal_solid(do_start,do_end)
!===============================================
      use module_material,only : cpsld
      integer,intent(in) :: do_start,do_end
!
      do ICVL=do_start,do_end
        hhh(ICVL) = cpsld(-IMAT)*tmp(ICVL)
      enddo
      end subroutine cal_solid
!
      end subroutine cal_t2hcp_dc
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_t2cp(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,prs,rho,rmu,rmd,cpi,cp,cr)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- 1. Calculate specific heat from temperature
      use module_dimension
      use module_constant
      use module_metrix  ,only : yi=>dum_c1
      use module_model,only   : ical_vect,ical_dens
      use module_species,only  : gascns
      use module_metrix,only   : W1K8
!
      implicit none 
!
! --- [dummy arguments]
!
      integer,intent(in)  :: MAT_CVEXT (0:MXMAT)
      integer,intent(in)  :: MAT_DCIDX (0:MXMAT)
      integer,intent(in)  :: MAT_NO    (0:MXMAT)
      logical,intent(in)  :: mat_cal   (0:MXMAT)
      real*8 ,intent(in)  :: tmp       (MXALLCV)
                                ! temperature [K]
      real*8 ,intent(in)  :: yys (       MXALLCV,mxcomp)
                                ! mass fraction of species [-]
      real*8 ,intent(out) :: cpi(       MXALLCV,mxcomp)
                                ! specific heat of species [J/(kg K)]
      real*8 ,intent(out) :: cp        (MXALLCV)
                                ! specific heat of mixture gas [J/(kg K)]
      real*8 ,intent(out) :: cr        (MXALLCV)
!                               ! [J/(kg K)]
      real*8 ,intent(in)    :: rho   (        MXALLCV)
      real*8 ,intent(in)    :: prs   (        MXALLCV)
      real*8 ,intent(in)    :: RMD   (        MXALLCV)
      real*8 ,intent(in)    :: RMU   (        MXALLCV)
! --- [local entities]
!
      integer :: n,s,IMAT,IIMAT,ICVS,ICVE,IDCS,IDCE,ICVL,icom,ierr1=0
!
!
      if(ical_vect) then  !vector  !NOVECT
        do IIMAT=1,NMAT    !ICV=1,NCVX
        if( .not.mat_cal(IIMAT) ) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(IMAT>0) then
          call cal_fluid_V1(ICVS,ICVE,1)
          call cal_fluid_V1(IDCS,IDCE,2)
        else                            ! solid
          call cal_solid(ICVS,ICVE)
          call cal_solid(IDCS,IDCE)
        endif
        enddo
      else
        do IIMAT=1,NMAT    !ICV=1,NCVX
        if( .not.mat_cal(IIMAT) ) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
!
        if(IMAT>0) then
          if(ical_dens==4) then
          else
            call cal_fluid(ICVS,ICVE,1)
            call cal_fluid(IDCS,IDCE,2)
          endif
        else                            ! solid
          call cal_solid(ICVS,ICVE)
          call cal_solid(IDCS,IDCE)
        endif
        enddo
      endif
!
      RETURN
!////////////////////////////////////////////////////////////////
      contains
!
!=================================================t2cp
      subroutine cal_fluid(do_start,do_end,imode)
!=================================================
      use module_species,only : c_pi,c_p,sigma_yr,ACT!,sigma_yr
      integer,intent(in) :: do_start,do_end,imode
!
!      do ICVL=do_start,do_end
!      do ICOM=1,ncomp
!      cpi(ICVL,ICOM)=1579.d0
!      enddo
!      cp(ICVL)=1579.d0
!      cr(ICVL)=sigma_yr(ncomp,yi(1:ncomp))  ! [J/(kg K)]
!      enddo
!      return
!
      do ICVL=do_start,do_end
      do ICOM=1,ncomp
! --- specific heat of species [J/(kg K)]
      cpi(ICVL,ICOM)=c_pi(tmp(ICVL),ICOM)*dble(ACT(ICOM))
      yi(ICOM)=yys(ICVL,ICOM)*dble(ACT(ICOM))
! --- specific heat of mixture gas [J/(kg K)]
      enddo
      cp(ICVL)=c_p(ncomp,tmp(ICVL),yi(1:ncomp))
      cr(ICVL)=sigma_yr(ncomp,yi(1:ncomp))  ! [J/(kg K)]
      enddo
      end subroutine cal_fluid
!
!

!=====================================================
      subroutine cal_fluid_V1(do_start,do_end,imode) 
!=====================================================
!      use module_species,only : c_pi,c_p,sigma_yr
      use module_species,only : a7,Ri,t3
      integer,intent(in) :: do_start,do_end,imode
      real*8             :: ttt,c_pi
      integer   :: rj
!
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      ttt=tmp(ICVL)
      rj=(3-sign(1,int(ttt-t3(2,ICOM))))/2
      c_pi=a7(1,rj,ICOM)
     &        +ttt*(a7(2,rj,ICOM)
     &        +ttt*(a7(3,rj,ICOM)
     &        +ttt*(a7(4,rj,ICOM)
     &        +ttt* a7(5,rj,ICOM))))
      cpi(ICVL,ICOM)=c_pi*Ri(ICOM)
      enddo
      enddo
!
      cp(do_start:do_end)=0.d0
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      cp(ICVL)=cp(ICVL)+yys(ICVL,ICOM)*cpi(ICVL,ICOM)
      enddo
      enddo
!
      cr(do_start:do_end)=0.d0
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      cr(ICVL)=cr(ICVL)+yys(ICVL,ICOM)*Ri(ICOM)
      enddo
      enddo
!
      end subroutine cal_fluid_V1
!
!============================================
      subroutine cal_solid(do_start,do_end)
!============================================
      use module_material,only : cpsld
      integer,intent(in) :: do_start,do_end
      do ICVL=do_start,do_end
        cr(ICVL) = 1.d0
        cp(ICVL) = cpsld(-IMAT)
        cpi(ICVL,1:ncomp) = cp(ICVL)
      enddo
      end subroutine cal_solid
!
      end subroutine cal_t2cp
!
!------------------------------------------
! --- Calculate enthalpy from temperature
!------------------------------------------
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_t2h(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &                   tmp,yys,hi,h,rho,prs)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_constant
      use module_metrix,only  : yi =>dum_c1
      use module_model,only   : ical_vect,ical_dens
      use module_species,only : r_wm
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: MAT_CVEXT (0:MXMAT)
      integer,intent(in)  :: MAT_DCIDX (0:MXMAT)
      integer,intent(in)  :: MAT_NO    (0:MXMAT)
      logical,intent(in)  :: mat_cal   (0:MXMAT)
      real*8 ,intent(inout)  :: tmp    (MXALLCV)
                             ! temperature [K]
      real*8 ,intent(in)  :: yys (       MXALLCV,mxcomp)
                             ! mass fraction of species [-]
      real*8 ,intent(out) :: hi (       MXALLCV,mxcomp)
                             ! enthalpy of species [J/kg]
      real*8 ,intent(out) :: h         (MXALLCV)
                             ! enthalpy [J/kg]
      real*8 ,intent(in)    :: rho   (       MXALLCV)
      real*8 ,intent(in)    :: prs   (       MXALLCV)

!
! --- [local entities]
!
      integer :: n,s,IMAT,IIMAT,ICVS,ICVE,IDCS,IDCE,ICVL,icom
!
      if(ical_vect) then  !NOVECT
        do IIMAT=1,NMAT    !ICV=1,NCVX
        if(.not.mat_cal(IIMAT)) cycle
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        IDCS = MAT_DCIDX(IIMAT-1)+1
        IDCE = MAT_DCIDX(IIMAT)
        if(IMAT.gt.0) then
          call cal_fluid_V2(ICVS,ICVE)
          call cal_fluid_V2(IDCS,IDCE)
        else                            ! solid
          call cal_solid(ICVS,ICVE)
          call cal_solid(IDCS,IDCE)
        endif
        enddo
      else
        do IIMAT=1,NMAT    !ICV=1,NCVX
        if(.not.mat_cal(IIMAT)) cycle
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        IDCS = MAT_DCIDX(IIMAT-1)+1
        IDCE = MAT_DCIDX(IIMAT)
        if(IMAT.gt.0) then
          if(ical_dens==4) then
          else
            call cal_fluid(ICVS,ICVE)
            call cal_fluid(IDCS,IDCE)
          endif
        else                            ! solid
          call cal_solid(ICVS,ICVE)
          call cal_solid(IDCS,IDCE)
        endif
        enddo
      endif
!
!///////////////////////////////////////////////////////////////////
      contains
!
!================================================
      subroutine cal_fluid(do_start,do_end)
!================================================
      use module_species,only : h_t,Ri,
     &                          enthalpy_ref,
     &                          href,h_ref,
     &                          enthalpy,gascns
!
      implicit none
      integer,intent(in) :: do_start,do_end
!      real*8, intent(inout) :: yyi(mxcomp)
      real*8  :: tt
!
      do ICVL=do_start,do_end
      tt=tmp(ICVL)
      do icom=1,ncomp
!-------------------------------
! --- enthalpy of species [J/kg]
!-------------------------------
        yi(ICOM)=yys(ICVL,ICOM)
        hi(ICVL,icom)=h_ref(icom,tt)*Ri(icom)
      enddo
      h(ICVL)=enthalpy_ref(ncomp,tt,yi(1:ncomp))
                                       ! enthalpy [J/kg]
      enddo
!
      end subroutine cal_fluid
!
!
!================================================
      subroutine cal_fluid_V2(do_start,do_end)
!================================================
!      use module_species,only : h_t,Ri,enthalpy_ref,href,h_ref
      use module_species,only : href,Ri,t3,a7,act
!
      implicit none      
      integer,intent(in) :: do_start,do_end
      real*8  :: ttt,h_reff
      integer :: rj
!
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      ttt=tmp(ICVL)
      rj=(3-sign(1,int(ttt-t3(2,ICOM))))/2
      h_reff   =a7(6,rj,icom)
     &   +ttt*(a7(1,rj,icom)
     &   +ttt*(a7(2,rj,icom)*0.5d0
     &   +ttt*(a7(3,rj,icom)*1.d0/3.d0
     &   +ttt*(a7(4,rj,icom)*0.25d0
     &   +ttt*(a7(5,rj,icom)*0.2d0)))))-href(icom)
      hi(ICVL,icom)=h_reff*Ri(icom)
      enddo
      enddo
!
      h(do_start:do_end)=0.d0
      do icom=1,ncomp
      do ICVL=do_start,do_end
      h(ICVL)=h(ICVL)+dble(act(icom))*yys(ICVL,ICOM)*hi(ICVL,icom)!*Ri(icom)
      enddo
      enddo
!
      end subroutine cal_fluid_V2

!================================================
      subroutine cal_solid(do_start,do_end)
!================================================
      use module_material,only : cpsld
      integer,intent(in) :: do_start,do_end
      do ICVL=do_start,do_end
      h(ICVL) = cpsld(-IMAT)*tmp(ICVL) ! enthalpy [J/kg]
      enddo
      end subroutine cal_solid
!
      end subroutine cal_t2h
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_t2h_dc(mph,MAT_DCIDX,MAT_NO,mat_cal,
     &                   tmp,yys,h,prs,rho)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_constant
      use module_metrix  ,only : yi =>dum_c1 !yyi=>ys      
      use module_model,only    : ical_vect,ical_dens
      use module_species,only : r_wm
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mph
      integer,intent(in)  :: MAT_DCIDX (0:MXMAT)
      integer,intent(in)  :: MAT_NO    (0:MXMAT)
      logical,intent(in)  :: mat_cal   (0:MXMAT)
      real*8 ,intent(in)  :: tmp       (MXALLCV)
                             ! temperature [K]
      real*8 ,intent(in)  :: yys (       MXALLCV,mxcomp)
                             ! mass fraction of species [-]
      real*8 ,intent(out) :: h         (MXALLCV)
                             ! enthalpy [J/kg]
      real*8 ,intent(in)    :: rho   (       MXALLCV)
      real*8 ,intent(in)    :: prs   (       MXALLCV)

!
! --- [local entities]
!
      integer :: n,s,IMAT,IIMAT,ICVS,ICVE,IDCS,IDCE,ICVL,icom,ierr1=0
!
!
      if(ical_vect)then
        do IIMAT=1,NMAT    !ICV=1,NCVX
        if(.not.mat_cal(IIMAT)) cycle
        IMAT = MAT_NO(IIMAT)
        IDCS = MAT_DCIDX(IIMAT-1)+1
        IDCE = MAT_DCIDX(IIMAT)
        if(IMAT.gt.0) then
          call cal_fluid_VDC(IDCS,IDCE)
        else                            ! solid
          call cal_solid(IDCS,IDCE)
        endif
        enddo
      else
        do IIMAT=1,NMAT    !ICV=1,NCVX
        if(.not.mat_cal(IIMAT)) cycle
        IMAT = MAT_NO(IIMAT)
        IDCS = MAT_DCIDX(IIMAT-1)+1
        IDCE = MAT_DCIDX(IIMAT)
!
        if(IMAT.gt.0) then
           
           if(ical_dens==4)then
           else
              call cal_fluid(IDCS,IDCE)
           endif
        else                    ! solid
          call cal_solid(IDCS,IDCE)
        endif
        enddo
      endif
!
      RETURN
!
!///////////////////////////////////////////////////////////////////
      contains
!
!================================================
      subroutine cal_fluid(do_start,do_end)
!================================================
      use module_species,only : h_t,Ri,enthalpy_ref,href,h_ref
      integer,intent(in) :: do_start,do_end
      real*8  :: tt
!
      do ICVL=do_start,do_end
      tt=tmp(ICVL)
      yi(1:ncomp)=yys(ICVL,1:ncomp)
      h(ICVL)=enthalpy_ref(ncomp,tt,yi(1:ncomp))
      enddo
!
      end subroutine cal_fluid



!================================================
      subroutine cal_fluid_VDC(do_start,do_end)
!================================================
      use module_species,only : href,Ri,t3,a7
!
      implicit none      
      integer,intent(in) :: do_start,do_end
      real*8  :: ttt,h_reff
      integer :: rj
!
      h(do_start:do_end)=0.d0
      do ICOM=1,ncomp
      do ICVL=do_start,do_end
      ttt=tmp(ICVL)
      rj=(3-sign(1,int(ttt-t3(2,ICOM))))/2
      h_reff   =a7(6,rj,icom)
     &   +ttt*(a7(1,rj,icom)
     &   +ttt*(a7(2,rj,icom)*0.5d0
     &   +ttt*(a7(3,rj,icom)*1.d0/3.d0
     &   +ttt*(a7(4,rj,icom)*0.25d0
     &   +ttt*(a7(5,rj,icom)*0.2d0)))))-href(icom)
      h(ICVL)=h(ICVL)+yys(ICVL,ICOM)*h_reff*Ri(icom)
      enddo
      enddo
!
      end subroutine cal_fluid_VDC

!================================================
      subroutine cal_solid(do_start,do_end)
!================================================
      use module_material,only : cpsld
      integer,intent(in) :: do_start,do_end
!
      do ICVL=do_start,do_end
      h(ICVL) = cpsld(-IMAT)*tmp(ICVL) ! enthalpy [J/kg]
      enddo
      end subroutine cal_solid
!
      end subroutine cal_t2h_dc
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_temp
     &(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,yys,pp0,prs,rho,tmp,p_grdc)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!zhang-cvd prs=>pp0
! --- [module arguments]  zhang-cvd
!
      use module_dimension
      use module_constant
      use module_initial 
      use module_species,only : gascns,wm,sigma_yr,sw,r_wm,act,
     &                          rho_real_PR
      use module_model,  only : ical_topt,ical_t,mach0,idrdp,comp,
     &                          ical_dens
      use module_chemreac,ONLY : ical_suf
      use module_metrix  ,only : yyi=>ys,W1K9
      use module_model,only    : ical_vect
!      use module_particle,only : IFLAG_COAL,icoal
!
! 1. Calculate temperature from equation of state
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(inout) :: p_grdc
      INTEGER,INTENT(INOUT) :: MAT_NO   (0:MXMAT)
      INTEGER,INTENT(INOUT) :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX( 0:MXMAT)
      logical,INTENT(INOUT) :: mat_cal  (  0:MXMAT)
      real*8 ,intent(in)    :: yys(       MXALLCV,mxcomp)
      real*8 ,intent(inout) :: pp0      ( MXALLCV)
      real*8 ,intent(inout) :: prs      ( MXALLCV)
      real*8 ,intent(in)    :: rho      ( MXALLCV)
      real*8 ,intent(out)   :: tmp      ( MXALLCV)
!
! --- [local entities]
!
      real*8  :: cr,cr1,dum1,dum2,dum3,dum4,dum5,TTT,PPP,VVV
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IBFS,IBFE,IBFL,ierr1=0
      integer :: IMODE,ICOM,IMD
      integer :: IMAT_U
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
!
!      real*8,allocatable  :: yyi(:)
!
!      ALLOCATE(yyi(mxcomp),stat=ierr1)
!
! --- 
!
!zhang-cvd
      if(idrdp==comp.and.ical_suf==1) then
        p_grdc=1.d0
        do IIMAT=1,NMAT
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        IDCS = MAT_DCIDX(IIMAT-1)+1
        IDCE = MAT_DCIDX(IIMAT)
        do ICVL=ICVS,ICVE
        if(prs(ICVL)<1.d0) then
          pp0(ICVL)=1.d0
          prs(ICVL)=1.d0
          p_grdc=0.d0
        endif
        enddo
        do ICVL=IDCS,IDCE
          pp0(ICVL)=max(pp0(ICVL),1.d0)
          prs(ICVL)=max(prs(ICVL),1.d0)
        enddo
        enddo
      endif
!
!zhang-cvd
!      if(idrdp==mach0.and.IFLAG_COAL/=icoal) then
!      if(idrdp==mach0.or.idrdp==comp) then
      if(idrdp==mach0.or.idrdp==comp) then   !OK 
        if(ical_vect) then  !NOVECT
          do IIMAT=1,NMAT    !ICV=1,NCV
          IMAT=MAT_NO(IIMAT)
          if(IMAT.lt.0.or..not.mat_cal(IIMAT)) cycle 
          ICVS=MAT_CVEXT(IIMAT-1)+1 
          ICVE=MAT_CVEXT(IIMAT)
          W1K9(ICVS:ICVE)=0.d0
          do ICOM=1,ncomp
          do ICVL=ICVS,ICVE
          W1K9(ICVL)=W1K9(ICVL)+yys(ICVL,ICOM)*r_wm(ICOM)*act(icom)
          enddo
          enddo
          do ICVL=ICVS,ICVE
          tmp(ICVL)=pp0(ICVL)/(gascns*W1K9(ICVL)*rho(ICVL))
          enddo
          enddo
        else
          if(ical_dens==4) then
          else
            do 100 IIMAT=1,NMAT    !ICV=1,NCV
            IMAT=MAT_NO(IIMAT)
            if(IMAT.lt.0) goto 100 
            if(.not.mat_cal(IIMAT)) goto 100 
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
              cr=0.d0
              do 101 ICOM=1,ncomp
              cr=cr+yys(ICVL,ICOM)*r_wm(ICOM)*act(icom)
  101         continue
              tmp(ICVL)=pp0(ICVL)/(gascns*cr*rho(ICVL))
            enddo
 100        enddo
          endif
        endif
      endif
!
!      DEALLOCATE(yyi)
      return
      end subroutine cal_temp
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_ustar(yp,rnu,akappa,awfnc,u,utaul,icall,ical)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_constant
      use module_turbparm,only : yplsm
! --- uplus=uplus; u=u (in); u=utau (out)
!
!  1. Solve friction velocity at wall making use of log-law
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: icall,ical
      real*8, intent(in)    :: yp,rnu
      real*8, intent(in)    :: akappa,awfnc
      real*8, intent(inout) :: u
      real*8, intent(out)   :: utaul
!
!
! --- [local entities]
!
      real*8,parameter :: r1p2=1.d0/2.d0,r1p6=1.d0/6.d0,eps=1.d-8
      real*8,parameter :: r16=0.16d0,r32=0.032d0,r1a=1.d0
      real*8  :: yplus,duplus,yre,C0,ak1,ak2,ak3,vk,utau,uplus,ustar
      real*8  :: upl,gutau,ggutau,spald,spaldd
      integer :: lop
!
      uplus=one
      vk=1.d0/akappa
!
      if(icall==1) then
!^^^^^^^^^^^^^^
! --- log-law: 
!^^^^^^^^^^^^^^
        yre=abs(u*yp/rnu)
        if(yre.le.ONE) then
          utaul=dsqrt(rnu*u/yp)
          return
        endif
        ustar=vk*log(yre)+awfnc
        do 100 lop=1,100
        duplus=(ustar-uplus-vk*log(uplus))/(1.d0+vk/uplus)
!
        uplus=uplus+duplus
        if(duplus.lt.eps*uplus) goto 101
  100   enddo
        
        write(*,*) '### program error -1- (cal_ustar)',icall,ical
        CALL FFRABORT(1,'cal_ustar/wall-function not conv.')
 101    continue
        utau=u/(uplus+SML)
        yplus=yp*utau/(rnu+SML)
        if(yplus.lt.yplsm) then
          utau=dsqrt(rnu*u/yp)
        endif
        utaul=utau
!
        return
      else if(icall==2) then
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-< Spalding-law >
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- yRe=abs(u*yp/rnu); yplus=yp*u/rnu=yRe/uplus
        yRe=abs(u*yp/rnu+SML)
        if(yre.le.ONE) then
          utaul=dsqrt(rnu*u/yp)
          return
        endif
!
        ak1=akappa
        ak2=akappa**2
        ak3=akappa**3
        c0=exp(-ak1*awfnc)
!
        do 200 lop=1,100
        spald=uplus
     &       +c0*(exp(ak1*uplus)-1.d0
     &       -ak1*uplus-r1p2*ak2*uplus**2
     &       -r1p6*ak3*uplus**3)-yre/uplus
        spaldd=1.d0+c0*(ak1*exp(ak1*uplus)-ak1
     &       -ak2*uplus-r1p2*ak3*uplus**2)+yre/uplus**2
        duplus=-spald/(spaldd+SML)
        uplus=uplus+duplus
        
        if(abs(duplus).lt.eps*uplus) goto 201
  200   continue

        write(*,*) '### program error -2- (cal_ustar) icall=2',ical
        CALL FFRABORT(1,'cal_ustar')
  201   continue
        utaul=u/(uplus+SML)
!        u=utaul
        return
!-------------------------------------------------------------------
      elseif(icall==3) then
        utau=dsqrt(rnu*u/yp)+SML
        DO 300 lop=1,10
        upl=U/UTAU
        gutau=upl-UTAU*yp/(rnu+SML)
     &       +EXP(-akappa*awfnc)*(EXP(akappa*upl)
     &       -r1a-akappa*upl-r1p2*akappa*upl*akappa*upl
     &       -(akappa*upl)*(akappa*upl)*(akappa*upl)*r1p6)
        ggutau=-upl/(UTAU+SML)-yp/(rnu+SML)-EXP(-akappa*awfnc)
     &      *(akappa*upl/(UTAU+SML)*EXP(akappa*upl)
     &      -akappa*upl/(UTAU+SML)-r16*upl*upl/(UTAU+SML)
     &      -r32*upl*upl*upl/(UTAU+SML))
        UTAU=UTAU-gutau/(ggutau+SML)
 300    CONTINUE
        utaul=UTAU
        return
      endif
!
      end subroutine cal_ustar
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_statis(iter,vel,vlasum,
     &           CVCENT,MAT_NO,MAT_CVEXT,MAT_INDEX,
     &           pp0,prs,rho,tmp,yys,aks,rmut)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only    : lenfnm,ifle,ifll,getfil
      use module_model,only : idrdp,incomp,mach0,icaltb,noturb,
     &                         ke,sles,dles,lles,dns,RSM,SDES,ical_prt
     
      use module_les  ,only : ISTART,n_ave,
     &                        iuvw_ave_rms_re,
     &                        ip_ave,it_ave,imu_ave,
     &                        icomp_ave,irans_ave,
     &                        it_rms,ip_rms,imu_rms,
     &                        icomp_rms,irans_rms,
     &                        iuvwt_re,imaxmin,
     &                        iuvw_rans_re,ista_momery,
     &                        iuvwc_re,mol_rslt
!
      use module_scalar,only : sclname,rns_scl
      use module_material,only : ical_sld,rotati,ishaft,end,begin,
     &                           rot_ang
      use module_time,    only : steady
      use module_metrix,only   : STA_SAV, vlr=>d2work1
      use module_species,only  : spcnam,wm,r_wm
      use module_trace,only    : PRT_VF
!
!  1. statistics
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: iter
      real*8 ,intent(in)     :: vel    (   MXALLCV,3)
      real*8 ,intent(inout)  :: vlasum (   MXCV)
      real*8 ,intent(in)     :: pp0    (   MXALLCV)
      real*8 ,intent(in)     :: prs    (   MXALLCV)
      real*8 ,intent(in)     :: rho    (   MXALLCV)
      real*8 ,intent(in)     :: tmp    (   MXALLCV)
      real*8 ,intent(in)     :: yys(       MXALLCV,MXcomp)
      real*8 ,intent(in)     :: aks(       MXALLCVR,MXrans)
      real*8 ,intent(in)     :: rmut   (   MXALLCV)
      real*8 ,intent(in)     :: CVCENT (3, MXALLCV)
      INTEGER,INTENT(IN)     :: MAT_NO (   0:MXMAT)
      INTEGER,INTENT(IN)     :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)     :: MAT_INDEX( 0:MXMAT)
      
!
! --- [local entities]
!
      integer :: i,j,l,ICV,ifl1,ifl2,nrnsx,ios,icom,ierr1=0,ierr=0
      integer :: rename,status,ICVS,ICVE,ICVL,IIMAT,IMAT
      character(lenfnm),save :: fnam1,fnam2
      logical :: fexist=.false.,unstdy=.false.
      real*8  :: titer,titer0
      real*8  :: aveles(3),rmsles(6)
      real*8  :: avet,avep
      real*8  :: rmst,rmsp
      real*8  :: reuvwt(3)
      integer :: in_ave,iin_ave
!
      real*8  :: unit(3),tha,bb(3,3),rbb(3,3),radi(3),vr(3),dr,v0(3)
      real*8  :: sinth,costh,dum1,dum2,dum3
!
      integer :: iuvw_ave_rms_rex
      integer :: ip_avex,it_avex,imu_avex
      integer :: ip_rmsx,it_rmsx,imu_rmsx,ical_prtx
      integer :: iuvwt_rex,imaxminx,ista_momeryx
!
      real*8,allocatable  :: rmscomp(:),rmsrans(:)
      real*8,allocatable  :: avecomp(:),averans(:)
      real*8,allocatable  :: rerans(:)
      integer,allocatable :: icomp_avex(:),irans_avex(:)
      integer,allocatable :: icomp_rmsx(:),irans_rmsx(:)
      integer,allocatable :: iuvw_rans_rex(:)
      integer,allocatable :: iuvwc_rex(:)
!
      ALLOCATE(rmscomp(ncomp),rmsrans(nrans),
     &         avecomp(ncomp),averans(nrans),rerans(nrans),
     &         icomp_avex(ncomp),irans_avex(nrans),
     &         icomp_rmsx(ncomp),irans_rmsx(nrans),
     &         iuvw_rans_rex(nrans),stat=ierr1)
      allocate(iuvwc_rex(ncomp),stat=ierr1)
!
! --- 
!
      fexist=.false.
      call getfil(ifl1,fnam1,'statis')
      call getfil(ifl2,fnam2,'statis.tmp')
      if(NPE.gt.1) then
        fnam1=TRIM(statis)
        fnam2=TRIM(TRIM(statis)//'.'//'tmp')
      else
        fnam1='statis'
        fnam2='statis'//'.'//'tmp'
      endif
!
                       nrnsx=nrans
      if(.not.rns_scl) nrnsx=0
!
      unstdy=icaltb.eq.sles.or.icaltb.eq.dles.or.
     &       icaltb.eq.lles.or.icaltb.eq.dns.or.icaltb.eq.SDES.or.
     &       (.not.steady)
!      inquire(file=fnam1,exist=fexist)
      if( ista_momery /= 1 ) then
         inquire(file=fnam1,exist=fexist)

         if(.not.unstdy) then
            if(fexist) call unlink(fnam1)
            DEALLOCATE(rmscomp,rmsrans,avecomp,averans,rerans,
     &           icomp_avex,irans_avex,icomp_rmsx,irans_rmsx,
     &           iuvw_rans_rex)
            return
         endif
      endif


!      if(.not.unstdy) then
!        if(fexist) call unlink(fnam1)
!        DEALLOCATE(rmscomp,rmsrans,avecomp,averans,rerans,
!     &           icomp_avex,irans_avex,icomp_rmsx,irans_rmsx,
!     &           iuvw_rans_rex)
!        return
!      endif
!
      if( ista_momery == 1 .and. iter < ISTART ) return

      if( ista_momery /= 1 ) then
         if(iter.lt.ISTART) then
            inquire(file=fnam1,exist=fexist)
            if(fexist) call unlink(fnam1)
            DEALLOCATE(rmscomp,rmsrans,avecomp,averans,rerans,
     &           icomp_avex,irans_avex,icomp_rmsx,irans_rmsx,
     &           iuvw_rans_rex)
            return
         endif
      endif
      
      allocate(vlr(mxallcv,ncomp),stat=ierr1)
      if(mol_rslt==1) then
         if(ierr /= 0 ) then
        write(ifle,*) 'allocating array error in write_result=>yys'
        call FFRABORT(1,'allocating array at write_result')
         endif
         do 200 IIMAT = 1, NMAT
            IMAT=MAT_NO(IIMAT)
            if(IMAT<0) cycle
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_INDEX(IIMAT)
            do 210 ICVL=ICVS,ICVE
               dum1=0.d0
               do ICOM=1,ncomp
                  dum1=dum1+yys(ICVL,ICOM)*r_wm(ICOM)
               enddo
               do ICOM=1,ncomp
                  dum2=wm(ICOM)*dum1
                  vlr(ICVL,ICOM)=yys(ICVL,ICOM)/dum2
               enddo
 210        enddo
 200     enddo
      else
         vlr(1:mxallcv,1:ncomp) = yys(1:mxallcv,1:ncomp)
      endif



!      if(iter.lt.ISTART) then
!        inquire(file=fnam1,exist=fexist)
!        if(fexist) call unlink(fnam1)
!        DEALLOCATE(rmscomp,rmsrans,avecomp,averans,rerans,
!     &           icomp_avex,irans_avex,icomp_rmsx,irans_rmsx,
!     &           iuvw_rans_rex)
!        return
!      endif
!
! --- 
!
!      allocate(vlr(MXALLCV,3),stat=ierr)
!      if(ical_sld==1.or.ical_sld==2.or.ical_sld==3) then
!        if(ierr.ne.0) then
!          write(ifle,*) 'allocating array error in write_result=>sld'
!          call FFRABORT(1,'allocating array at cal_statis')
!        endif
!        do  IIMAT=1,NMAT
!          IMAT=MAT_NO(IIMAT)
!          tha=rot_ang(IMAT)
!          unit(1)=end(1,IMAT)-begin(1,IMAT)
!          unit(2)=end(2,IMAT)-begin(2,IMAT)
!          unit(3)=end(3,IMAT)-begin(3,IMAT)
!          CALL rotth(unit,tha,bb)
!          do i=1,3
!          do j=1,3
!          rbb(i,j)=bb(j,i)
!          enddo
!          enddo
!          costh=cos(tha)
!          sinth=sin(tha)
!          ICVS=MAT_CVEXT(IIMAT-1)+1
!          ICVE=MAT_INDEX(IIMAT)
!          do ICVL=ICVS,ICVE
!          radi(1)=CVCENT(1,ICVL)-begin(1,IMAT)
!          radi(2)=CVCENT(2,ICVL)-begin(2,IMAT)
!          radi(3)=CVCENT(3,ICVL)-begin(3,IMAT)
!          call AXB_UNIT_C(unit,radi,vr)
!          dr=radi(1)*unit(1)+radi(2)*unit(2)+radi(3)*unit(3)
!          radi(1)=radi(1)-dr*unit(1)
!          radi(2)=radi(2)-dr*unit(2)
!          radi(3)=radi(3)-dr*unit(3)
!          dr=dsqrt(radi(1)*radi(1)+radi(2)*radi(2)+radi(3)*radi(3))
!          v0(1)=dr*rotati(IMAT)*vr(1)
!          v0(2)=dr*rotati(IMAT)*vr(2)
!          v0(3)=dr*rotati(IMAT)*vr(3)
!
!          vlr(ICVL,1)=rbb(1,1)*(vel(ICVL,1)+v0(1))
!     &               +rbb(1,2)*(vel(ICVL,2)+v0(2))
!     &               +rbb(1,3)*(vel(ICVL,3)+v0(3))
!          vlr(ICVL,2)=rbb(2,1)*(vel(ICVL,1)+v0(1))
!     &               +rbb(2,2)*(vel(ICVL,2)+v0(2))
!     &               +rbb(2,3)*(vel(ICVL,3)+v0(3))
!          vlr(ICVL,3)=rbb(3,1)*(vel(ICVL,1)+v0(1))
!     &               +rbb(3,2)*(vel(ICVL,2)+v0(2))
!     &               +rbb(3,3)*(vel(ICVL,3)+v0(3))
!          enddo
!        enddo
!      else
!        do i=1,3
!        vlr(1:NCV,i)=vel(1:NCV,i)
!        enddo
!      endif
!
! --- 
!
      if(ista_momery==1) then
         
      endif
!
      if(iter.eq.ISTART) then

        if( ista_momery /= 1 ) then
            inquire(file=fnam1,exist=fexist)
            if(fexist) call unlink(fnam1)
            inquire(file=fnam2,exist=fexist)
            if(fexist) call unlink(fnam2)
!
            open(ifl1,file=fnam1,form='unformatted',
     &           status='unknown',iostat=ios)
            write(ifl1) iuvw_ave_rms_re,ip_ave,imu_ave,
     &                  it_ave,it_rms,ip_rms,imu_rms,
     &                  iuvwt_re,imaxmin,ista_momery,ical_prt
            if(ncomp>0) then
               write(ifl1) (icomp_ave(i),i=1,ncomp),
     &                     (icomp_rms(i),i=1,ncomp)
            end if
!
            if(ncomp>0) then
               write(ifl1) (iuvwc_re(i),i=1,ncomp)
            endif
!
            if(nrnsx>0) then
               write(ifl1) (irans_ave(i),i=1,nrnsx),
     &              (irans_rms(i),i=1,nrnsx),
     &              (iuvw_rans_re(i),i=1,nrnsx)
            end if

        endif

!        inquire(file=fnam1,exist=fexist)
!        if(fexist) call unlink(fnam1)
!        inquire(file=fnam2,exist=fexist)
!        if(fexist) call unlink(fnam2)
!
!        open(ifl1,file=fnam1,form='unformatted',
!     &       status='unknown',iostat=ios)
!
!        write(ifl1) iuvw_ave_rms_re,ip_ave,it_ave,it_rms,ip_rms,
!     &              iuvwt_re,imaxmin,ista_momery
!        if(ncomp>0) then
!           write(ifl1) (icomp_ave(i),i=1,ncomp),
!     &                 (icomp_rms(i),i=1,ncomp)
!        end if
!        if(nrnsx>0) then
!           write(ifl1) (irans_ave(i),i=1,nrnsx),
!     &                 (irans_rms(i),i=1,nrnsx),
!     &                 (iuvw_rans_re(i),i=1,nrnsx)
!        end if
!        write(ifl1) iuvw_ave_rms_re,
!     &             ip_ave,it_ave,
!     &             (icomp_ave(i),i=1,ncomp),(irans_ave(i),i=1,nrnsx),
!     &             it_rms,ip_rms,
!     &             (icomp_rms(i),i=1,ncomp),(irans_rms(i),i=1,nrnsx),
!     &             iuvwt_re,
!     &             (iuvw_rans_re(i),i=1,nrnsx)
!
        in_ave=0
        if(iuvw_ave_rms_re==1) then
          if(ista_momery==1) then
            if(n_ave<9) 
     &      call FFRABORT(1,'ERR:n_ave<9, call your supportor ')
            do i=1,3
            STA_SAV(1:NCV,i)=vel(1:NCV,i)
            enddo
            do i=1,6
            STA_SAV(1:NCV,i+3)=0.d0   !zhangcor
            enddo
            in_ave=9
          else
            do i=1,3
            write(ifl1) (vel(ICV,i),ICV=1,NCV)
            enddo
            vlasum=0.d0
            do i=1,6
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
            enddo
          endif
        endif
!
        if(ip_ave==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=prs(1:NCV)
          else
            write(ifl1) (prs(ICV),ICV=1,NCV)
          endif
        endif
!
        if(ip_rms.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=0.d0
          else
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!Mu
        if(imu_ave==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=rmut(1:NCV)
          else
            write(ifl1) (rmut(ICV),ICV=1,NCV)
          endif
        endif
!
!Mu
        if(imu_rms.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=0.d0
          else
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(imaxmin==1) then
          if(ista_momery==1) then
            do i=1,3
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=vel(1:NCV,i) 
            enddo
            do i=1,3
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=vel(1:NCV,i) 
            enddo
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=prs(1:NCV)
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=prs(1:NCV)
          else
            do i=1,3
            write(ifl1) (vel(ICV,i),ICV=1,NCV)
            enddo
            do i=1,3
            write(ifl1) (vel(ICV,i),ICV=1,NCV)
            enddo
            write(ifl1) (prs(ICV),ICV=1,NCV)
            write(ifl1) (prs(ICV),ICV=1,NCV)
          endif
        endif
!
        if(it_ave==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=tmp(1:NCV)
          else
            write(ifl1) (tmp(ICV),ICV=1,NCV)
          endif
        endif
!
        if(it_rms.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=0.d0
          else
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(iuvwt_re.eq.1) then
          if(ista_momery==1) then
            do i=1,3
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=0.d0       !!write(ifl1) (vlasum(ICV),ICV=1,NCV)
            enddo
          else
            do i=1,3
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
            enddo
          endif
        endif
!
        do icom=1,ncomp
        if(icomp_ave(icom).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=vlr(1:NCV,icom)
          else
            write(ifl1) (vlr(ICV,icom),ICV=1,NCV)
          endif
        endif
!
        if(icomp_rms(icom).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=0.d0
          else
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
        enddo
!
        do icom = 1, ncomp
           if(iuvwc_re(icom)==1) then
              if(ista_momery==1) then
                 do i = 1, 3
                    in_ave=in_ave+1
                    if(in_ave>n_ave) 
     &       call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
                    STA_SAV(1:NCV,in_ave)=0.d0
                 enddo
              else
                 do i = 1, 3
                    write(ifl1) (vlasum(ICV),ICV=1,NCV)
                 enddo
              endif
           endif
        enddo
!
        do i=1,nrnsx
        if(irans_ave(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=aks(1:NCV,i)
          else
            write(ifl1) (aks(ICV,i),ICV=1,NCV)
          endif
        endif
!
        if(irans_rms(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=0.d0
          else
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(iuvw_rans_re(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=0.d0
          else
            write(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
        enddo
!
        if(ical_prt==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,in_ave)=PRT_VF(1:NCV)
          else
            write(ifl1) (PRT_VF(ICV),ICV=1,NCV)
          endif
        endif        
!
        if(in_ave/=n_ave) then
          write(ifle,*) 'ERR:FFR bug: in_ave=,n_ave=  ',in_ave,n_ave
          call FFRABORT(1,'ERR:in_ave/=n_ave, call your supportor ')
        endif
!
        if( ista_momery /= 1 ) then
           close(ifl1)
        endif
      else
!
        titer=dble(iter-ISTART)
        titer0=titer-1.d0
!
        if( ista_momery == 1 ) then
           iuvw_ave_rms_rex=iuvw_ave_rms_re
           ip_avex=ip_ave
           imu_avex=imu_ave
           imu_rmsx=imu_rms
           it_avex=it_ave
           icomp_avex(:)=icomp_ave(:)
           irans_avex(:)=irans_ave(:)
           it_rmsx=it_rms;ip_rmsx=ip_rms
           icomp_rmsx(:)=icomp_rms(:)
           irans_rmsx(:)=irans_rms(:)
           iuvwt_rex=iuvwt_re
           imaxminx=imaxmin
           iuvw_rans_rex(:)=iuvw_rans_re(:)
           ista_momeryx=ista_momery
!           ! Tagami
           iuvwc_rex(:) = iuvwc_re(:)
!           ! imagaT
           ical_prtx=ical_prt
        else
           iuvw_ave_rms_rex=0
           ip_avex=0
           imu_avex=0
           imu_rmsx=0
           it_avex=0
           icomp_avex(:)=0
           irans_avex(:)=0
           it_rmsx=0;ip_rmsx=0
           icomp_rmsx(:)=0
           irans_rmsx(:)=0
           iuvwt_rex=0
           imaxminx=0
           iuvw_rans_rex(:)=0
           ista_momeryx=0
!           ! Tagami
           iuvwc_rex(:) = 0
!           ! imagaT
           ical_prtx=0
        endif

        
        if(ista_momery /= 1 ) then
           inquire(file=fnam1,exist=fexist)
           if(.not.fexist) then
              write(ifle,*) '   ### ERR: Statise file does NOT exist'
              write(ifle,*) '   ### Set [NSTART] > ',iter
              call FFRABORT(1,'cal_statis')
           endif
!
           open(ifl1,file=fnam1,
     &          form='unformatted',status='unknown',iostat=ios)
           open(ifl2,file=fnam2,
     &          form='unformatted',status='unknown',iostat=ios)
           rewind ifl1
           read(ifl1) iuvw_ave_rms_rex,ip_avex,imu_avex,it_avex,
     &                it_rmsx,ip_rmsx,imu_rmsx,
     &                iuvwt_rex,imaxminx,ista_momeryx,ical_prtx
           if(ista_momeryx/=ista_momery) then
              call FFRABORT(1,'ERR: Can NOT change ')
           endif
           if(ncomp>0) then
              read(ifl1) (icomp_avex(i),i=1,ncomp),
     &                   (icomp_rmsx(i),i=1,ncomp)
           end if
!
           if(ncomp>0) then
               read(ifl1) (iuvwc_rex(i),i=1,ncomp)
           endif
!
           if(nrnsx>0) then
              read(ifl1) (irans_avex(i),i=1,nrnsx),
     &             (irans_rmsx(i),i=1,nrnsx),
     &             (iuvw_rans_rex(i),i=1,nrnsx)
           end if

           write(ifl2)  iuvw_ave_rms_re,ip_ave,imu_ave,
     &                  it_ave,it_rms,ip_rms,imu_rms,
     &                  iuvwt_re,imaxmin,ista_momery,ical_prt
           if(ncomp>0) then
              write(ifl2) (icomp_ave(i),i=1,ncomp),
     &                    (icomp_rms(i),i=1,ncomp)
           end if
!
           if(ncomp>0) then
               write(ifl2) (iuvwc_re(i),i=1,ncomp)
           endif
           if(nrnsx>0) then
              write(ifl2) (irans_ave(i),i=1,nrnsx),
     &             (irans_rms(i),i=1,nrnsx),
     &             (iuvw_rans_re(i),i=1,nrnsx)
           end if
        endif
!
!        inquire(file=fnam1,exist=fexist)
!        if(.not.fexist) then
!          write(ifle,*) '   ### ERR: Statise file does NOT exist'
!          write(ifle,*) '   ### Set [NSTART] > ',iter
!          call FFRABORT(1,'cal_statis')
!        endif
!
!        open(ifl1,file=fnam1,
!     &       form='unformatted',status='unknown',iostat=ios)
!        open(ifl2,file=fnam2,
!     &       form='unformatted',status='unknown',iostat=ios)
!        rewind ifl1
!        read(ifl1) iuvw_ave_rms_rex,ip_avex,it_avex,it_rmsx,ip_rmsx,
!     &             iuvwt_rex,imaxminx,ista_momeryx
!        if(ista_momeryx/=ista_momery) then
!          call FFRABORT(1,'ERR: Can NOT change ')
!        endif
!        if(ncomp>0) then
!           read(ifl1) (icomp_avex(i),i=1,ncomp),
!     &                (icomp_rmsx(i),i=1,ncomp)
!        end if
!        if(nrnsx>0) then
!           read(ifl1) (irans_avex(i),i=1,nrnsx),
!     &                (irans_rmsx(i),i=1,nrnsx),
!     &                (iuvw_rans_rex(i),i=1,nrnsx)
!        end if
!
!        write(ifl2)  iuvw_ave_rms_re,ip_ave,it_ave,it_rms,ip_rms,
!     &               iuvwt_re,imaxmin,ista_momery
!        if(ncomp>0) then
!           write(ifl2) (icomp_ave(i),i=1,ncomp),
!     &                 (icomp_rms(i),i=1,ncomp)
!        end if
!        if(nrnsx>0) then
!           write(ifl2) (irans_ave(i),i=1,nrnsx),
!     &                 (irans_rms(i),i=1,nrnsx),
!     &                 (iuvw_rans_re(i),i=1,nrnsx)
!        end if
!1 --- ////
        in_ave=0
        iin_ave=0
        do 100 i=1,3
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave)
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          do ICV=1,NCV
          vlasum(ICV)=vel(ICV,i)
          enddo
        endif
!2
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vel(ICV,i))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
 100    continue
!
        
        do 110 i=1,3
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vel(ICV,i)*vel(ICV,i))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
 110    continue
!
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vel(ICV,1)*vel(ICV,2))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vel(ICV,2)*vel(ICV,3))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(iuvw_ave_rms_rex.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(iuvw_ave_rms_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vel(ICV,3)*vel(ICV,1))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
! --- ////
        if(ip_avex==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(ip_ave==1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+prs(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(ip_rmsx.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(ip_rms.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+prs(ICV)*prs(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
! Mu ---
!
        if(imu_avex==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(imu_ave==1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+rmut(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!!
        if(imu_rmsx.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(imu_rms.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+rmut(ICV)*rmut(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
! --- ////
        do i=1,3
          if(imaxminx==1) then
            if(ista_momery==1) then
              in_ave=in_ave+1
              if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
              vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
            else
              read(ifl1) (vlasum(ICV),ICV=1,NCV)
            endif
          else
            vlasum(:)=vel(:,i)
          endif
!
          if(imaxmin==1) then
            do ICV=1,NCV
            vlasum(ICV)=max(vlasum(ICV),vel(ICV,i))
            enddo
            if(ista_momery==1) then
              iin_ave=iin_ave+1
              if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
              STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
            else
              write(ifl2) (vlasum(ICV),ICV=1,NCV)
            endif
          endif
        enddo
!
        do i=1,3
          if(imaxminx==1) then
            if(ista_momery==1) then
              in_ave=in_ave+1
              if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
              vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
            else
              read(ifl1) (vlasum(ICV),ICV=1,NCV)
            endif
          else
            vlasum(:)=vel(:,i)
          endif
!
          if(imaxmin==1) then
            do ICV=1,NCV
            vlasum(ICV)=min(vlasum(ICV),vel(ICV,i))
            enddo
            if(ista_momery==1) then
              iin_ave=iin_ave+1
              if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
              STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
            else
              write(ifl2) (vlasum(ICV),ICV=1,NCV)
            endif
          endif
        enddo
!
        if(imaxminx==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=prs(:)
        endif
!
        if(imaxmin==1) then
          do ICV=1,NCV
            vlasum(ICV)=max(vlasum(ICV),prs(ICV))
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(imaxminx==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=prs(:)
        endif
!
        if(imaxmin==1) then
          do ICV=1,NCV
            vlasum(ICV)=min(vlasum(ICV),prs(ICV))
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
! --- ////
        if(it_avex==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(it_ave==1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+tmp(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(it_rmsx.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(it_rms.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+tmp(ICV)*tmp(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        do i=1,3
        if(iuvwt_rex.eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(iuvwt_re.eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vel(ICV,i)*tmp(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
        enddo
! --- ////
        do i=1,ncomp
        if(icomp_avex(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(icomp_ave(i).eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vlr(ICV,i))/titer
!vlasum(ICV)=(vlasum(ICV)*titer0+yys(ICV,i))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!

        if(icomp_rmsx(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(icomp_rms(i).eq.1) then
           do ICV=1,NCV
           vlasum(ICV)=(vlasum(ICV)*titer0
     &           +vlr(ICV,i)*vlr(ICV,i))/titer
           enddo
           if(ista_momery==1) then
             iin_ave=iin_ave+1
             if(iin_ave>n_ave) 
     &       call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
                 STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
           else
             write(ifl2) (vlasum(ICV),ICV=1,NCV)
           endif
        endif
        enddo
!///
        ! Tagami
        do icom = 1, ncomp
           do i = 1, 3
              if( iuvwc_rex(icom) == 1 ) then
                 if(ista_momery==1) then
                    in_ave=in_ave+1
                    if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
                    vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
                 else
                    read(ifl1) (vlasum(ICV),ICV=1,NCV)
                 endif
              else
                 vlasum(:)=0.d0
              endif
!
              if( iuvwc_re(icom) == 1 ) then
                 do ICV = 1, NCV
                    vlasum(ICV)=(vlasum(ICV)*titer0
     &                   +vel(ICV,i)*vlr(ICV,icom))/titer
                 enddo
                 if(ista_momery==1) then
                    iin_ave=iin_ave+1
                    if( in_ave /= iin_ave ) then
                       call FFRABORT(1,"ERR:iin_ave/=in_ave")
                    endif
                    if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
                    STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
!                    STA_SAV(1:NCV,in_ave)=vlasum(1:NCV)
                 else
                    write(ifl2) (vlasum(ICV),ICV=1,NCV)
                 endif
              endif
           enddo
        enddo




! --- ////
        do i=1,nrnsx
        if(irans_avex(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(irans_ave(i).eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+aks(ICV,i))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        if(irans_rmsx(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(irans_rms(i).eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+aks(ICV,i)*aks(ICV,i))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
!
        do l=1,3
        if(iuvw_rans_rex(i).eq.1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(iuvw_rans_re(i).eq.1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+vel(ICV,l)*aks(ICV,i))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif
        enddo
        enddo
!
        if(ical_prtx==1) then
          if(ista_momery==1) then
            in_ave=in_ave+1
            if(in_ave>n_ave) 
     &      call FFRABORT(1,'ERR:in_ave>n_ave, call your supportor ')
            vlasum(1:NCV)=STA_SAV(1:NCV,in_ave)
          else
            read(ifl1) (vlasum(ICV),ICV=1,NCV)
          endif
        else
          vlasum(:)=0.d0
        endif
!
        if(ical_prt==1) then
          do ICV=1,NCV
          vlasum(ICV)=(vlasum(ICV)*titer0+PRT_VF(ICV))/titer
          enddo
          if(ista_momery==1) then
            iin_ave=iin_ave+1
            if(iin_ave>n_ave) 
     &      call FFRABORT(1,'ERR:iin_ave>n_ave, call your supportor ')
            STA_SAV(1:NCV,iin_ave)=vlasum(1:NCV)
          else
            write(ifl2) (vlasum(ICV),ICV=1,NCV)
          endif
        endif        
! --- ////
        if(ista_momery /= 1 ) then
           close(ifl1)
           close(ifl2)
!
           call unlink(fnam1)
           status=rename(fnam2,fnam1)
           if(status.ne.0) then
              write(*,'(1X,a,I8)') 
     &             'ERR: rename functio error, status= ',status
           call FFRABORT(1,'  ### ERR: rename function in cal_statis')
           endif
        endif

!

!        titer=dble(iter-ISTART)
!        titer0=titer-1.d0
!        do 100 ICV=1,NCV
!        do 101 l=1,3
!        vlasum(l,ICV)=(vlasum(l,ICV)*titer0+vel(ICV,l))/titer
! 101    continue
!        do 102 l=1,3
!        vlrsum(l,ICV)=(vlrsum(l,ICV)*titer0+vel(ICV,l)*vel(ICV,l))/titer
! 102    continue
!        vlrsum(4,ICV)=(vlrsum(4,ICV)*titer0+vel(ICV,1)*vel(ICV,2))/titer
!        vlrsum(5,ICV)=(vlrsum(5,ICV)*titer0+vel(ICV,2)*vel(ICV,3))/titer
!        vlrsum(6,ICV)=(vlrsum(6,ICV)*titer0+vel(ICV,3)*vel(ICV,1))/titer
! 100    continue
!
      endif
!
      
      DEALLOCATE (rmscomp,rmsrans,avecomp,averans,rerans,
     &           icomp_avex,irans_avex,icomp_rmsx,irans_rmsx,
     &           iuvw_rans_rex)
!
      deallocate(vlr)
      deallocate(iuvwc_rex)


      return
!
      end subroutine cal_statis
!
!

!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE cal_AXBEQC(A,B,C,D)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------
!     A X B = C 
!---------------------------------------
      real*8 ,intent(in)   :: A(3),B(3)
      real*8 ,intent(OUT)  :: C(3),D
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
      C(3)=(A(1)*B(2)-A(2)*B(1))
      C(2)=(A(3)*B(1)-A(1)*B(3))
      C(1)=(A(2)*B(3)-A(3)*B(2))
      D=DSQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))
!
      RETURN
      END subroutine cal_AXBEQC
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE cal_injcoord(deltt,time,iter,ictl,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  CVVOLM,CVCENT)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      USE module_usersub, ONLY : src_r,usrno,usryes
      use module_material,only : nofld
      use module_metrix  ,only : rcomp,ocomp=>ys
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter
      real*8 ,intent(in)    :: deltt,time
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      integer,intent(inout) :: ictl
!
! --- [local entities]
!
      integer             :: IMAT,IIMAT,ICVS,ICVE,ICVL,IMAT_U,ICV
!      real*8,allocatable  :: rcomp(:),ocomp(:)
      real*8              :: dum0,dum1,dum2,dum3,dum4,vol,xxx,yyy,zzz
!
      if(src_r==usrno) return
!
!      allocate(rcomp(MXCOMP),ocomp(MXCOMP))
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT.gt.0) then
         IMAT_U=nofld(IMAT)
         DO ICVL=ICVS,ICVE
         ICV=MAT_CV(ICVL)
         if(NPE.gt.1) ICV=NOD_IDg(ICV)
           xxx=CVCENT(1,ICVL)
           yyy=CVCENT(2,ICVL)
           zzz=CVCENT(3,ICVL)
           vol=CVVOLM(ICVL)
           rcomp(:)=0.d0
           ocomp(:)=0.d0
           call  user_src_r(1,ictl,
     &     deltt,iter,time,ICV,IMAT_U,1,ncomp,xxx,yyy,zzz,
     &     0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,rcomp,vol,
     &     dum0,dum1,dum2,dum3,dum4,ocomp)
           enddo
       endif
      enddo
!      
!      deallocate(rcomp,ocomp)
!
      return
      end SUBROUTINE cal_injcoord
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine CONTROL_RPM(time,deltt,iter,MAT_NO)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_boundary,only : idis,nbcnd,kdbcnd,MAT_BCIDX,
     &                           LBC_INDEX,kdsld,rotsld
      use module_material,only : strtim,rotup,rotati,rotN0,rot_ang,
     &                           ishaft,nflud,sav_ang,rot_angold,
     &                           domegadt
      use module_constant
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter
      real(8),intent(in)    :: time,deltt
      INTEGER,INTENT(IN)    :: MAT_NO(0:MXMAT)
!
! --- [local entities]
!
      real(8)  :: appr,drpm
      integer  :: imat,IIMAT1,IIMAT2,IMAT1,IMAT2,nb,kd
!
      do imat=1,nflud
        appr=(dble(iter)-strtim(imat))/(rotup(imat)+1.d-10)
        rotati(imat)=rotN0(imat)*min(max(0.d0,appr),1.d0)
        sav_ang(imat)=sav_ang(imat)+deltt*rotati(imat)
        rot_ang(imat)=rot_ang(imat)+deltt*rotati(imat)
! only for linear 
        if(iter>=strtim(imat) .and.
     &      iter<=strtim(imat)+rotup(imat)) then
          domegadt(imat)=rotN0(imat)/(dble(rotup(imat))*deltt)
        else
          domegadt(imat)=0.d0
        end if
!        IF(MOD(ITER,20000)==0) THEN
!          rot_ang(imat)=sav_ang(imat)
!        ENDIF
        ishaft(imat)=0
        if(abs(rotati(imat))>SML) ishaft(imat)=1
        rot_angold(imat)=rot_ang(imat)-deltt*rotati(imat)
      enddo
!
      do nb=1,nbcnd
      kd=kdbcnd(0,nb)
      rotsld(nb)=0
      if(kd==kdsld) then
        IIMAT1=MAT_BCIDX(nb,1)
        IIMAT2=MAT_BCIDX(nb,2)
        IMAT1=MAT_NO(IIMAT1)
        IMAT2=MAT_NO(IIMAT2)
        drpm=abs(rotati(imat1)-rotati(imat2))
        if(drpm>SML) then
          rotsld(nb)=1
        endif
      endif
      enddo
!
      end subroutine CONTROL_RPM
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      Subroutines for flamelet calculation
!       by T.Tominaga
!
!      1. "cal_bv"
!            : subroutine to set laminar burning velocity
!      2. "cal_trvbv"
!            : subroutine to evaluate turbulent effect on burning velocity
!      3. "cal_dmpbv"
!            : subroutine to evaluate wall damping of burning velocity
!      4. "cal_prpgf"
!            : subroutine to evaluate propagation flux
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_prpgf
     &          (MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,mat_cal
     &          ,LVEDGE,LBC_SSF,LCYCSF,LCYCOLD,OPPANG
     &          ,SFAREA,SFCENT,wiface,wifsld,CVVOLM
     &          ,rva,akstmp,flmspd,grdc,vctr,rva_s)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_scalar,only  : calgeq,calxi,caltemp,calh,
     &                          igeq,ixi,itemp,ihflmlt
!
      implicit none
!
      integer,intent(in)    :: MAT_CV    (   MXALLCV)
      integer,intent(in)    :: MAT_CVEXT (   0:MXMAT)
      integer,intent(in)    :: MAT_NO    (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      integer,intent(in)    :: LCYCOLD   (MXSSFBC_SLD)

      logical,intent(in)    :: mat_cal   (   0:MXMAT)

      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: wifsld    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)

      real*8 ,intent(in)    :: akstmp    (  MXALLCVR)
      real*8 ,intent(in)    :: flmspd    (  MXALLCV)
      real*8 ,intent(in)    :: rva       (   MXCVFAC)
      real*8 ,intent(inout) :: grdc      (   MXALLCV,3,3)
      real*8 ,intent(out)   :: rva_s     (   MXCVFAC_F)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      integer :: IIMAT,IMAT
      integer :: ICVL,ICVS,ICVE
      integer :: ICFL,ICFS,ICFE
      integer :: ICVA,ICVB
      real*8  :: grx,gry,grz,dum1
      real*8  :: wi1,wi2
      real*8  :: grdmg
!---------------------------
! --- Cal. gradient of G
!---------------------------
      call grad_cell(1,51,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,akstmp,grdc(:,:,1))
!
      call dc_symprv
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    LCYCOLD,wifsld,OPPANG,
     &    SFAREA,SFCENT,grdc(:,:,1),1,0)
!
      do IIMAT=1,NMAT           !ICV=1,NCV
       if(.not.mat_cal(IIMAT)) CYCLE
       IMAT=MAT_NO(IIMAT)
       ICVS=MAT_CVEXT(IIMAT-1)+1
       ICVE=MAT_CVEXT(IIMAT)
       do ICVL=ICVS,ICVE
!-----------------------------------
! --- Magnitude of normal vector
!-----------------------------------
        grdmg = dsqrt(grdc(ICVL,1,1)*grdc(ICVL,1,1)
     &               +grdc(ICVL,2,1)*grdc(ICVL,2,1)
     &               +grdc(ICVL,3,1)*grdc(ICVL,3,1))
! --- For ZERO-gradient
        grdmg = max(grdmg,1.0d-10)
! --- Unit normal vector
        grdc(ICVL,1,1) = grdc(ICVL,1,1) / grdmg
        grdc(ICVL,2,1) = grdc(ICVL,2,1) / grdmg
        grdc(ICVL,3,1) = grdc(ICVL,3,1) / grdmg
       enddo
      enddo
!------------------------------------------------------
! --- Cal. Curvature of iso-surface
! --- Cal. gradient of unit normal vector of iso-surface
!------------------------------------------------------
      call grad_cell(1,52,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,grdc(:,1,1),grdc(:,:,2))
!
      call dc_symprv
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    LCYCOLD,wifsld,OPPANG,
     6    SFAREA,SFCENT,grdc(:,:,2),1,0)
!
      grdc(:,1,3) = grdc(:,1,2)
!
!
      call grad_cell(1,53,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,grdc(:,2,1),grdc(:,:,2))
!
      call dc_symprv
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    LCYCOLD,wifsld,OPPANG,
     &    SFAREA,SFCENT,grdc(:,:,2),1,0)
!
      grdc(:,2,3) = grdc(:,2,2)
!
!
      call grad_cell(1,53,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,grdc(:,3,1),grdc(:,:,2))
!
      call dc_symprv
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    LCYCOLD,wifsld,OPPANG,
     &    SFAREA,SFCENT,grdc(:,:,2),1,0)
!
      grdc(:,3,3) = grdc(:,3,2)
!---------------------------------------------
! --- Cal. Curvature of iso-surface
!---------------------------------------------
      do IIMAT=1,NMAT           !ICV=1,NCV
       if(.not.mat_cal(IIMAT)) CYCLE
       IMAT=MAT_NO(IIMAT)
       ICVS=MAT_CVEXT(IIMAT-1)+1
       ICVE=MAT_CVEXT(IIMAT)
       do ICVL=ICVS,ICVE
        dum1 = grdc(ICVL,1,3)
     &        +grdc(ICVL,2,3)
     &        +grdc(ICVL,3,3)
        grdc(ICVL,3,3) =dum1
       enddo
      enddo
!---------------------------------------------
! --- Cal. Propagation direction
!---------------------------------------------
      do IIMAT=1,NMAT     !ICF=1,NCVFAC
       if(.not.mat_cal(IIMAT)) CYCLE
       ICFS=MAT_CFIDX(IIMAT-1)+1
       ICFE=MAT_CFIDX(IIMAT)
       do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
!
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
!---------------------------------------------
! --- Propagation velocity @ surface-center
!---------------------------------------------
        grx  = wi1*flmspd(ICVA)*grdc(ICVA,1,1)
     &        +wi2*flmspd(ICVB)*grdc(ICVB,1,1)
        gry  = wi1*flmspd(ICVA)*grdc(ICVA,2,1)
     &        +wi2*flmspd(ICVB)*grdc(ICVB,2,1)
        grz  = wi1*flmspd(ICVA)*grdc(ICVA,3,1)
     &        +wi2*flmspd(ICVB)*grdc(ICVB,3,1)
!---------------------------------------------
! --- Cal. (Convection + Propagation) velocity
!---------------------------------------------
        rva_s(ICFL)=rva(ICFL)
     &             -(SFAREA(1,ICFL)*grx
     &              +SFAREA(2,ICFL)*gry
     &              +SFAREA(3,ICFL)*grz)*SFAREA(4,ICFL)
!
       enddo
      enddo
!
      return
!
      end subroutine cal_prpgf
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      Subroutine to set laminar burning velocity
!
!       Chemical conditions
!           1. Equivalence ratio "\phi"
!               a) defined in module_scalar (Mixture "A" or "A and B")
!                  : phiA(,phiB)
!               b) calculate from mixture fraction
!                  : aks(:,ixi,:)
!           2. Pressure "P_0" : pp0
!           3. Unburned mixture temperature "T_u"
!               a) defined in module_scalar : tmpu
!               b) calculate from h_flamelet : aks(:,ihfl,:)
!
!          ==>  Laminar burning velocity : flmspd(ICVL)
!               !!ATTENTION!!
!                 "flmspd" is an array for turbulent B.V..
!                 Here, however, it store laminar B.V. temporarily.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_bv(MAT_CVEXT,MAT_CAL,MAT_INDEX,aks,pp0,flmspd)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
      use module_dimension
      use module_constant
      use module_scalar,only  :  calgeq,calxi,caltemp,calh ! -- logical flags
     &                          ,igeq,ixi,itemp,ihflmlt    ! -- Indexes
!
     &                          ,fd_lbv1,fd_lbv2,nlbv1,nlbv2, ! -- Poly. Coef.
     &                           lbv_XXi,XXi0
      use module_model,only   : ical_vect
      use module_metrix,only  : W1K9,W1K8
!
      implicit none
!
      real*8 ,intent(in)  :: aks (MXALLCVR,MXRANS,2)
      real*8 ,intent(in)  :: pp0 (MXALLCV)
      real*8 ,intent(out) :: flmspd(MXALLCV)
!
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      LOGICAL,INTENT(IN)    :: MAT_CAL   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
!
      real*8 :: xi
      real*8 :: bvsl,bvsl1
!
      integer :: IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,IDCL,ilbv
!
! --- [Initialize]
      flmspd= 0.0d0
      bvsl= 0.0d0
!
      if(ical_vect) then 
        if (calxi) then   !partial-pre-mixed
!-------------------------------------
! --- [Set laminar burning velocity]
! ----- All Cell
!-------------------------------------
          W1K8(:)=0.D0
          W1K9(:)=0.D0
          do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_INDEX(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ilbv=nlbv1,0,-1
          DO ICVL=ICVS,ICVE
            xi=max(0.d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            !call phi_poly_nopdf(xi,fd_lbv1(0:nlbv1),nlbv1,bvsl1)
            W1K8(ICVL)=W1K8(ICVL)*xi+fd_lbv1(ilbv)
          ENDDO
          enddo
          do ilbv=nlbv2,0,-1
          DO ICVL=ICVS,ICVE
            xi=max(0.d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            !call phi_poly_nopdf(xi,fd_lbv2(0:nlbv2),nlbv2,bvsl1)
            W1K9(ICVL)=W1K9(ICVL)*xi+fd_lbv2(ilbv)
          ENDDO
          enddo
          do ICVL=ICVS,ICVE
          xi=max(0.d0,aks(ICVL,ixi,1))
          xi=min(1.0d0,xi)
          if(xi<=lbv_XXi(2).and.xi>=lbv_XXi(1)) then
            flmspd(ICVL)=W1K8(ICVL)
          elseif(xi>=lbv_XXi(2).and.
     &      xi<=lbv_XXi(3).and.lbv_XXi(3)>0.d0) then
            flmspd(ICVL)=W1K9(ICVL)
          else
            flmspd(ICVL)=0.d0
          endif
          enddo

          enddo
        else !pre-mixed 
!-------------------------------------
! --- [Set laminar burning velocity]
! ----- All Cell
!-------------------------------------
          do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_INDEX(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
            xi = XXi0
            call phi_poly_nopdf(xi,fd_lbv1(0:nlbv1),nlbv1,bvsl)
            flmspd(ICVL) = bvsl
          enddo
          enddo
        endif
      else  !.NOT.ical_vect
!
! --- 
!
        if (calxi) then   !partial-pre-mixed
!-------------------------------------
! --- [Set laminar burning velocity]
! ----- All Cell
!-------------------------------------
          do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_INDEX(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
!-------------------------------------
! --- [Set local equivalence ratio] 
!-------------------------------------
            xi = max(0.d0,aks(ICVL,ixi,1))
            xi = min(1.0d0,xi)
            if(xi<=lbv_XXi(2).and.xi>=lbv_XXi(1)) then
              call phi_poly_nopdf(xi,fd_lbv1(0:nlbv1),nlbv1,bvsl1)
               bvsl=bvsl1
            elseif(xi>=lbv_XXi(2).and.
     &           xi<=lbv_XXi(3).and.lbv_XXi(3)>0.d0) then
              call phi_poly_nopdf(xi,fd_lbv2(0:nlbv2),nlbv2,bvsl1)
              bvsl=bvsl1
            else
              bvsl=0.d0
            endif
            flmspd(ICVL) = bvsl
          enddo
          enddo
        else !pre-mixed 
!-------------------------------------
! --- [Set laminar burning velocity]
! ----- All Cell
!-------------------------------------
          do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_INDEX(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
            xi = XXi0
            call phi_poly_nopdf(xi,fd_lbv1(0:nlbv1),nlbv1,bvsl)
            flmspd(ICVL) = bvsl
          enddo
          enddo
        endif

      endif
!
      contains
!=====================================================================
      subroutine phi_poly_nopdf(xi,fldat,nfl,phi)
!=====================================================================
      integer,intent(in)   :: nfl
      real*8,intent(in)    :: xi
      real*8,intent(inout) :: phi
      real*8,intent(in)    :: fldat(0:nfl)
      integer :: ifl
!
      phi=0.0d0
      do ifl=nfl,0,-1
        phi=fldat(ifl)+phi*xi
      enddo
!
      end subroutine phi_poly_nopdf
!
      end subroutine cal_bv
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      Subroutine to evaluate turbulent effect on burning velocity
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_trbbv(MAT_CVEXT,MAT_CAL,vflc,flmspd)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
      use module_dimension
      use module_constant
      use module_scalar,only  :  calgeq,calxi,caltemp,calh ! -- logical flags
     &                          ,igeq,ixi,itemp,ihflmlt,    ! -- Indexes
     &                          istmdl
      use module_model,only   : ical_vect
      use module_metrix,only  : W1K9,W1K8,W1K10
!
      implicit none
!
      real*8 ,intent(in)    :: vflc(MXALLCV)
      real*8 ,intent(inout) :: flmspd(MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      LOGICAL,INTENT(IN)    :: MAT_CAL   (   0:MXMAT)
!
      real*8  :: cfdamk,cfgmma
      real*8  :: vel_rat         ! Velocity Ratio   ( u'/s_L)
      real*8  :: trbfct          ! TuRBulent FaCTor (s_T/s_L)
!
! --- for s_T NIBUNGI method
      real*8  :: STres
      real*8  :: STSLt,STSLb,dum1
!
      integer :: IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,IDCL
!      integer :: istmdl
      integer :: nstsl
!
! --- choose model
      istmdl = 2   ! 1:Damk\"oler  2:Yakhot
!
! --- model coefficients
      cfgmma = 1.0d0
      cfdamk = 0.7d0
!
      
      if(istmdl.eq.1) then
!***********************************************************************
!      Based on Damk\"ohler Assumption
!***********************************************************************
! ----- Internal Cell
        do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
            trbfct = 1.0d0 + cfdamk*vflc(ICVL)/(flmspd(ICVL)+SML)
            trbfct = trbfct**cfgmma
            flmspd(ICVL)=flmspd(ICVL)*trbfct
          enddo
        enddo
!
      elseif(istmdl.eq.2) then
!***********************************************************************
!      Based on Yakhot model
!***********************************************************************
! ----- Internal Cell
        if(ical_vect) then
          do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          W1K8(ICVL)=10.0d-10                        !=>STSLb
          vel_rat=vflc(ICVL)/(flmspd(ICVL)+SML)
          W1K9(ICVL)=dmax1(dabs(vel_rat),exp(1.0d0)) !=>STSLt
          enddo

          do nstsl=1,20
          do ICVL=ICVS,ICVE
          vel_rat=vflc(ICVL)/(flmspd(ICVL)+SML)  !vel_rat
          trbfct=(W1K9(ICVL)+W1K8(ICVL))*0.5d0+SML
          STres=trbfct-dexp(vel_rat*vel_rat/trbfct/trbfct)
          if(STres.gt.0.0d0) W1K9(ICVL) = trbfct
          if(STres.lt.0.0d0) W1K8(ICVL) = trbfct
!          W1K10(ICVL)=min(trbfct,15.0d0)
          enddo
          enddo

          do ICVL=ICVS,ICVE
          dum1=(W1K9(ICVL)+W1K8(ICVL))*0.5d0
          trbfct = min(dum1,15.0d0)
          flmspd(ICVL)=flmspd(ICVL)*trbfct
          enddo

          enddo

        else  !(.not.ical_vect)
          do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
!
! ----- ST nibungi method START
!
            vel_rat =  vflc(ICVL)/(flmspd(ICVL)+SML)
            STSLt = dmax1(dabs(vel_rat),exp(1.0d0))
            STSLb = 10.0d-10
            do nstsl=1,50
               trbfct = (STSLt+STSLb)*0.5d0
               STres  = trbfct
     &                  - dexp(vel_rat*vel_rat/trbfct/trbfct)
               if(STres.gt.0.0d0) STSLt = trbfct
               if(STres.lt.0.0d0) STSLb = trbfct
               if(dabs(STres).le.1.0d-3) goto 110
            enddo
 110        continue
!
! ----- ST nibungi method END
!
            trbfct = min(trbfct,15.0d0)
            flmspd(ICVL)=flmspd(ICVL)*trbfct
          enddo
          enddo
          
        endif
!
      endif
!
      end subroutine cal_trbbv
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      Subroutine to evaluate damping function for burning velocity
!
!      Damping Function "f = 1-exp[-(y^+)/A]"
!
!      yplus : Wall distance y^+
!      fdmp  : Damping factor "f"
!      cfdmp : Coefficient "A" in f = 1-exp[-(y^+)/A]
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_dmpbv(MAT_CVEXT,MAT_CAL,yplus,flmspd)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
      use module_dimension
      use module_constant
!
      implicit none
!
      real*8 ,intent(in)    :: yplus (MXALLCV )
      real*8 ,intent(inout) :: flmspd(MXALLCV)
!
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      LOGICAL,INTENT(IN)    :: MAT_CAL   (   0:MXMAT)
!
      real*8  :: fdmp,cfdmp
!
      integer :: IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,IDCL
!
! --- [Set parameter]
      cfdmp = 25.0d0    !- Coefficient "A"
!
! --- [Set laminar burning velocity]
! ----- Internal Cell
        do IIMAT=1,NMAT    !ICV=1,NCV
          if(.not.mat_cal(IIMAT)) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
            fdmp = 1.0d0 - EXP(-yplus(ICVL)/cfdmp)
            flmspd(ICVL)=fdmp*flmspd(ICVL)
          enddo
        enddo
!
      end subroutine cal_dmpbv
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Perfectly Changed by T.Tominaga 2005/02/22 for 2-scalar
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      Subroutine to set Temperature "T" and Density "\rho"
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_flamelet
     &     (iter,MAT_CVEXT,MAT_DCIDX,MAT_CAL,MAT_NO,aks,rho,tmp,YYS)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_scalar,only  :  calgeq,calxi,caltemp,calh ! -- logical flags
     &                          ,igeq,ixi,itemp,ihflmlt    ! -- Indexes
     &                          ,fd_rho,fd_tmp
     &                          ,nrho,ntmp
      use module_scalar,only  :  XXi0,
     &                           Xfuel_rho,
     &                           Xfuel_temp,
     &                           XOxidant_rho,
     &                           XOxidant_temp,
     &                           XXi_max,igT_iter_flm
      use module_model,only   : ical_vect
      use module_metrix,only  : W1K9,W1K8
! --- Add T.Tominaga END
!
      implicit none
!

      real*8 ,intent(in)    :: aks  (MXALLCVR,MXRANS,2)
      real*8 ,intent(out)   :: rho  (MXALLCV ,2)
      real*8 ,intent(out)   :: tmp  (MXALLCV ,2)
      real*8 ,intent(inout) :: yys  (MXALLCV,mxcomp)
!

      INTEGER,INTENT(IN)    :: iter,MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      LOGICAL,INTENT(IN)    :: MAT_CAL   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(       0:MXMAT)
!
      real*8 :: fdmp,cfdmp
      real*8 :: gggeq
      real*8 :: xi
      real*8 :: dum
!
      real*8 :: rho_u,rho_b,rho_f,rho_o,sum_rho,dumXXi
      real*8 :: tmp_u,tmp_b,tmp_f,tmp_o,sum_tmp
!

      integer :: IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,IDCL
      integer :: i,itmp,irho
!
! --- Calculate Density & Temperature of Pure Fuel and Oxidant
!
!
      if(XXi_max<0.d0) then
        tmp_o = fd_tmp(0)
        rho_o = fd_rho(0)
        sum_tmp=0.0d0
        sum_rho=0.0d0
        do i=ntmp,0,-1
        sum_tmp=fd_tmp(i)+sum_tmp
        enddo
        do i=nrho,0,-1
        sum_rho=fd_rho(i)+sum_rho
        enddo
        tmp_f = sum_tmp
        rho_f = sum_rho
      elseif(calgeq.and..not.calxi) then
        tmp_o = fd_tmp(0)
        rho_o = fd_rho(0)
        sum_tmp=0.0d0
        sum_rho=0.0d0
        do i=ntmp,0,-1
        sum_tmp=fd_tmp(i)+sum_tmp
        enddo
        do i=nrho,0,-1
        sum_rho=fd_rho(i)+sum_rho
        enddo
        tmp_f = sum_tmp
        rho_f = sum_rho
      else
        tmp_o = XOxidant_temp
        rho_o = XOxidant_rho
        tmp_f = Xfuel_temp
        rho_f = Xfuel_rho
      endif
!
      if(ical_vect) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
!        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
!
        if(calgeq.and.calxi) then
          W1K8(:)=0.D0
          W1K9(:)=0.D0
          dumXXi=XXi_max
          if(XXi_max<0.d0) dumXXi=1.d0
!
          do itmp=ntmp,0,-1
          DO ICVL=ICVS,ICVE
            xi=min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
            !call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
            W1K8(ICVL)=W1K8(ICVL)*xi+fd_tmp(itmp)
          ENDDO
          enddo
!          
          do irho=nrho,0,-1
          DO ICVL=ICVS,ICVE
            xi=min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
            !call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
            W1K9(ICVL)=W1K9(ICVL)*xi+fd_rho(irho)
          ENDDO
          enddo
!
          DO ICVL=ICVS,ICVE
            xi=min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
            gggeq = max(0.0d0,aks(ICVL,igeq,1))
            gggeq = min(1.0d0,gggeq)
            !call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)
            tmp_u=tmp_o+(tmp_f-tmp_o)*xi
            !call phi_arithmetic(rho_o,rho_f,xi,rho_u)
            rho_u=rho_o+(rho_f-rho_o)*xi
            tmp_b=W1K8(ICVL)
            rho_b=W1K9(ICVL)
            !call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
            dum=tmp_u+(tmp_b-tmp_u)*gggeq
            tmp(ICVL,1)=dum
            !call phi_arithmetic(rho_u,rho_b,gggeq,dum)
            dum=rho_u+(rho_b-rho_u)*gggeq
            rho(ICVL,1)=1.0d0/dum
          ENDDO
! --- BC
          do itmp=ntmp,0,-1
          DO ICVL=IDCS,IDCE
            xi=min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
            !call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
            W1K8(ICVL)=W1K8(ICVL)*xi+fd_tmp(itmp)
          ENDDO
          enddo
!          
          do irho=nrho,0,-1
          DO ICVL=IDCS,IDCE
            xi=min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
            !call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
            W1K9(ICVL)=W1K9(ICVL)*xi+fd_rho(irho)
          ENDDO
          enddo

          do ICVL=IDCS,IDCE
            xi=min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
            gggeq = max(0.0d0,aks(ICVL,igeq,1))
            gggeq = min(1.0d0,gggeq)
            !call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)
            tmp_u=tmp_o+(tmp_f-tmp_o)*xi
            !call phi_arithmetic(rho_o,rho_f,xi,rho_u)
            rho_u=rho_o+(rho_f-rho_o)*xi
            tmp_b=W1K8(ICVL)
            rho_b=W1K9(ICVL)
            !call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
            dum=tmp_u+(tmp_b-tmp_u)*gggeq
            tmp(ICVL,1)=dum
            !call phi_arithmetic(rho_u,rho_b,gggeq,dum)
            dum=rho_u+(rho_b-rho_u)*gggeq
            rho(ICVL,1)=1.0d0/dum
          enddo
!
        elseif(calxi) then
!
          W1K8(:)=0.D0
          W1K9(:)=0.D0
          do itmp=ntmp,0,-1
          DO ICVL=ICVS,ICVE
            xi=max(0.0d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            
            !call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
            W1K8(ICVL)=W1K8(ICVL)*xi+fd_tmp(itmp)
          ENDDO
          enddo
          do irho=nrho,0,-1
          DO ICVL=ICVS,ICVE
            xi=max(0.0d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            !call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
            W1K9(ICVL)=W1K9(ICVL)*xi+fd_rho(irho)
          ENDDO
          enddo
          do ICVL=ICVS,ICVE
            xi=max(0.0d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            gggeq=1.0d0
            tmp_u=tmp_o+(tmp_f-tmp_o)*xi
            rho_u=rho_o+(rho_f-rho_o)*xi
            tmp_b=W1K8(ICVL)
            rho_b=W1K9(ICVL)
            !call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
            dum=tmp_u+(tmp_b-tmp_u)*gggeq
            tmp(ICVL,1)=dum
            !call phi_arithmetic(rho_u,rho_b,gggeq,dum)
            dum=rho_u+(rho_b-rho_u)*gggeq
            rho(ICVL,1)=1.0d0/dum
          enddo
! --- BC
          do itmp=ntmp,0,-1
          DO ICVL=IDCS,IDCE
            xi=max(0.0d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            
            !call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
            W1K8(ICVL)=W1K8(ICVL)*xi+fd_tmp(itmp)
          ENDDO
          enddo
          do irho=nrho,0,-1
          DO ICVL=IDCS,IDCE
            xi=max(0.0d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            !call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
            W1K9(ICVL)=W1K9(ICVL)*xi+fd_rho(irho)
          ENDDO
          enddo

          do ICVL=IDCS,IDCE
            xi=max(0.0d0,aks(ICVL,ixi,1))
            xi=min(1.0d0,xi)
            gggeq=1.0d0
            tmp_u=tmp_o+(tmp_f-tmp_o)*xi
            rho_u=rho_o+(rho_f-rho_o)*xi
            tmp_b=W1K8(ICVL)
            rho_b=W1K9(ICVL)
            !call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
            dum=tmp_u+(tmp_b-tmp_u)*gggeq
            tmp(ICVL,1)=dum
            !call phi_arithmetic(rho_u,rho_b,gggeq,dum)
            dum=rho_u+(rho_b-rho_u)*gggeq
            rho(ICVL,1)=1.0d0/dum
          enddo
!
        elseif(calgeq) then
!
          W1K8(:)=0.D0
          W1K9(:)=0.D0
          do itmp=ntmp,0,-1
          DO ICVL=IDCS,IDCE
            xi=XXi0
            !call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
            W1K8(ICVL)=W1K8(ICVL)*xi+fd_tmp(itmp)
          ENDDO
          enddo
          do irho=nrho,0,-1
          DO ICVL=IDCS,IDCE
            xi=XXi0
            !call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
            W1K9(ICVL)=W1K9(ICVL)*xi+fd_rho(irho)
          ENDDO
          enddo
!
          do ICVL=IDCS,IDCE
          xi=XXi0
          gggeq=max(0.0d0,aks(ICVL,igeq,1))
          gggeq=min(1.0d0,gggeq)
          tmp_u=tmp_o+(tmp_f-tmp_o)*xi
          rho_u=rho_o+(rho_f-rho_o)*xi
          tmp_b=W1K8(ICVL)
          rho_b=W1K9(ICVL)
          !call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
          dum=tmp_u+(tmp_b-tmp_u)*gggeq
          tmp(ICVL,1)=dum
          !call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          dum=rho_u+(rho_b-rho_u)*gggeq
          rho(ICVL,1)=1.0d0/dum
          enddo
! --- BC
          do itmp=ntmp,0,-1
          DO ICVL=ICVS,ICVE
            xi=XXi0
            !call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
            W1K8(ICVL)=W1K8(ICVL)*xi+fd_tmp(itmp)
          ENDDO
          enddo
          do irho=nrho,0,-1
          DO ICVL=ICVS,ICVE
            xi=XXi0
            !call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
            W1K9(ICVL)=W1K9(ICVL)*xi+fd_rho(irho)
          ENDDO
          enddo
!
          do ICVL=ICVS,ICVE
          xi=XXi0
          gggeq=max(0.0d0,aks(ICVL,igeq,1))
          gggeq=min(1.0d0,gggeq)
          tmp_u=tmp_o+(tmp_f-tmp_o)*xi
          rho_u=rho_o+(rho_f-rho_o)*xi
          tmp_b=W1K8(ICVL)
          rho_b=W1K9(ICVL)
          !call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
          dum=tmp_u+(tmp_b-tmp_u)*gggeq
          tmp(ICVL,1)=dum
          !call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          dum=rho_u+(rho_b-rho_u)*gggeq
          rho(ICVL,1)=1.0d0/dum
          enddo
!
        endif
!
        enddo
      else !(.NOT.ical_vect)
!-----------------------------
! --- scalar routine
!-----------------------------
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
!        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
!
        if(calgeq.and.calxi) then
          dumXXi=XXi_max
          if(XXi_max<0.d0) dumXXi=1.d0
          do ICVL=ICVS,ICVE
          xi    = min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
          gggeq = max(0.0d0,aks(ICVL,igeq,1))
          gggeq = min(1.0d0,gggeq)
!--- calculate density & temperature in unburnt side ------------------
          call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)
          call phi_arithmetic(rho_o,rho_f,xi,rho_u)
!
!--- calculate density & temperature in burnt side --------------------
          call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
          call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
!
!--- calculate density & temperature ----------------------------------
          call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
          tmp(ICVL,1)=dum
          call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          rho(ICVL,1)=1.0d0 / dum
          enddo
!
          do ICVL=IDCS,IDCE
          xi    = min(max(0.0d0,aks(ICVL,ixi,1)),dumXXi)
          gggeq = max(0.0d0,aks(ICVL,igeq,1))
          gggeq = min(1.0d0,gggeq)

          call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)

          call phi_arithmetic(rho_o,rho_f,xi,rho_u)

          call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)

          call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)

          call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)

          tmp(ICVL,1)=dum
          call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          rho(ICVL,1)=1.0d0 / dum
          enddo
!
        elseif(calxi) then

!
          do ICVL=ICVS,ICVE
          xi    = max(0.0d0,aks(ICVL,ixi,1))
          xi    = min(1.0d0,xi)
          gggeq = 1.0d0
!--- calculate density & temperature in unburnt side ------------------
          call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)
          call phi_arithmetic(rho_o,rho_f,xi,rho_u)
!
!--- calculate density & temperature in burnt side --------------------
          call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
          call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
!
!--- calculate density & temperature ----------------------------------
          call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
          tmp(ICVL,1)=dum
          call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          rho(ICVL,1)=1.0d0 / dum
          enddo

!
          do ICVL=IDCS,IDCE
          xi    = max(0.0d0,aks(ICVL,ixi,1))
          xi    = min(1.0d0,xi)
          gggeq = 1.0d0 
          call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)
          call phi_arithmetic(rho_o,rho_f,xi,rho_u)
          call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
          call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
          call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
          tmp(ICVL,1)=dum
          call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          rho(ICVL,1)=1.0d0 / dum
          enddo
!
        elseif(calgeq) then
!
          do ICVL=IDCS,IDCE
          xi = XXi0
          gggeq = max(0.0d0,aks(ICVL,igeq,1))
          gggeq = min(1.0d0,gggeq)
!--- calculate density & temperature in unburnt side ------------------
          call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)
          call phi_arithmetic(rho_o,rho_f,xi,rho_u)
!
!--- calculate density & temperature in burnt side --------------------
          call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
          call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
!
!--- calculate density & temperature ----------------------------------
          call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
          tmp(ICVL,1)=dum
          call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          rho(ICVL,1)=1.0d0 / dum
          enddo
!
          do ICVL=ICVS,ICVE
          xi = XXi0  !0.0d0
          gggeq = max(0.0d0,aks(ICVL,igeq,1))
          gggeq = min(1.0d0,gggeq)
          call phi_arithmetic(tmp_o,tmp_f,xi,tmp_u)
          call phi_arithmetic(rho_o,rho_f,xi,rho_u)
          call phi_poly_nopdf(xi,fd_tmp(0:ntmp),ntmp,tmp_b)
          call phi_poly_nopdf(xi,fd_rho(0:nrho),nrho,rho_b)
          call phi_arithmetic(tmp_u,tmp_b,gggeq,dum)
          tmp(ICVL,1)=dum
          call phi_arithmetic(rho_u,rho_b,gggeq,dum)
          rho(ICVL,1)=1.0d0 / dum
          enddo
!
        endif
!
        enddo

      endif
!
      call flamelet_ys_org
     &   (MAT_CVEXT,MAT_CAL,MAT_NO,aks,rho,tmp,yys)
!
      return
!/////////////////////////////////////////////////////////////////////
      contains
!=====================================================================
      subroutine phi_poly_nopdf(xi,fldat,nfl,phi)
!=====================================================================
      integer,intent(in)   :: nfl
      real*8,intent(in)    :: xi
      real*8,intent(inout) :: phi
      real*8,intent(in)    :: fldat(0:nfl)
      integer :: ifl
!
      phi=0.0d0
      do ifl=nfl,0,-1
        phi=fldat(ifl)+phi*xi
      enddo
!
      end subroutine phi_poly_nopdf
!=====================================================================
      subroutine phi_poly_pdf(gamma1,gamma2,fldat,nfl,phi)
!=====================================================================
      integer,intent(in)   :: nfl
      real*8,intent(in)    :: gamma1,gamma2
      real*8,intent(inout) :: phi
      real*8,intent(in)    :: fldat(0:nfl)
      integer :: ifl
!
      phi=0.0d0
      do ifl=nfl,0,-1
        phi=fldat(ifl)+phi*(gamma1+ifl)/(gamma1+gamma2+ifl)
      enddo
!
      end subroutine phi_poly_pdf
!
!======================================================================
      subroutine phi_arithmetic(phi0,phi1,prm,phi)
!======================================================================
      real*8,intent(in)    :: phi0,phi1,prm
      real*8,intent(inout) :: phi
!
      phi=0.0d0
      phi=phi0+(phi1-phi0)*prm
!
      end subroutine phi_arithmetic
!
!======================================================================
      subroutine phi_harmonic(phi0,phi1,prm,phi)
!======================================================================
      real*8,intent(in)    :: phi0,phi1,prm
      real*8,intent(inout) :: phi
!
      phi=0.0d0
      phi=phi0*phi1/(phi1+(phi0-phi1)*prm)
!
      end subroutine phi_harmonic
!
      end subroutine cal_flamelet
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      real(8) function Pvap(TTT)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      real*8,intent(in) :: TTT
!
      real*8 :: dum1,dum2,dum3,dum4
!
      dum1=(TTT-483.16d0)**2
      dum2=7.21379d0+(1.152d-5-4.787d-9*TTT)*dum1
      dum3=1.d0-647.31d0/TTT
      dum4=exp(dum2*dum3)
      Pvap=22.13D6*dum4
!
      return
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      end function Pvap
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_Dynamic_SGS
     &	(LVEDGE,LBC_SSF,LCYCSF,LFUTAU,
     &	MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &	SFAREA,SFCENT,wiface,CVCENT,CVVOLM,DISALL,FRSTCV,GRDC,OPPANG,
     &	rho,vel,rmu,rmut,aks,tmp,utau,wifsld,LCYCOLD,vctr,FIELD_U)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!------------------------
! --- [module arguments]
!------------------------
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_metrix,only  :	mvel1  =>d2work1
      use module_metrix,only  :	mvel2  =>d2work2
      use module_metrix,only  :	abs_dspd => W1K1
      use module_metrix,only  :	abs_sfil => W1K2
      use module_Dynamic_SGS,only : lij,mij
      use module_les ,only : dlalpha
      use module_metrix,only  :	mvel3 =>d2work3
      use module_metrix,only  :	mvel4 =>d2work4
      use module_model,only   : idrdp,mach0,comp
      use module_metrix,only  : R_deltV=>W1K8 

!
      implicit none
!
!------------------------
! --- [dummy arguments]
!------------------------
!
      integer,intent(in)  :: LVEDGE (2,MXCVFAC)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      INTEGER,INTENT(IN)  :: MAT_NO (  0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_DCIDX(0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      integer,intent(in)  :: LFUTAU (	  MXCV)
!
      real*8 ,intent(in)  :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)  :: SFCENT(3, MXCVFAC)
      real*8 ,intent(in)  :: wiface(   MXCVFAC)
      real*8 ,intent(in)  :: CVCENT(3, MXALLCV)
      real*8 ,intent(in)  :: CVVOLM(   MXALLCV)
      real*8 ,intent(in)  :: DISALL(   MXCV   )
      real*8 ,intent(in)  :: FRSTCV(   MXSSFBC)
      real*8 ,intent(inout) :: grdc(   MXALLCV,3,3)
      real*8 ,intent(in)  :: rho   (   MXALLCV)
      real*8 ,intent(in)  :: vel   (   MXALLCV,3)
      real*8 ,intent(in)  :: rmu   (   MXALLCV)
      real*8 ,intent(out) :: rmut  (   MXALLCV)
      real*8 ,intent(in)  :: aks   (   MXALLCVR,MXRANS)
      real*8 ,intent(inout):: utau (0: MXSSFBC)
      real*8 ,intent(in)  :: tmp   (   MXALLCV)
      real*8 ,intent(in)  :: wifsld(   MXSSFBC_SLD)
      integer,intent(in)  :: LCYCOLD(  MXSSFBC_SLD)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(in)  :: OPPANG(   MXSSFBC_SLD)
      real*8 ,intent(inout) :: FIELD_U  (MXCV_D,NFLID)
!----------------------
! --- [local entities]
!----------------------
      real(8) :: sij(7)
      real(8) :: delh,a1,b1
      real(8) :: csdles
      real(8) :: cs_max,cs_min,cs_ave,dum1,dum2,dum3,dum4
      integer :: nd,iv,ierr1=0,ierr2=0
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL
      real*8 ,parameter	:: r1pn=1.d0/3.d0
      real(8) :: Lkk
      integer :: ICFL,ICV1,ICV2,II,I,ICFE,ICFS,ICVLA,ICVLB
      integer :: fltr=1
!
!------------------------------------------------
! ---
      allocate(mvel1(MXALLCV,3),mvel2(MXALLCV,3),stat=ierr1)
      if(ierr1/=0) call
     &	 FFRABORT(1,'ERR:allocate error	in Dynamic_SGS')
      allocate(mvel3(MXALLCV,3),mvel4(MXALLCV,3),stat=ierr1)
      if(ierr1/=0) call
     &	 FFRABORT(1,'ERR:allocate error-1 in Dynamic_SGS')

!
!------------------------------------
! --- hat(deltV)/deltV 
!------------------------------------
!
      grdc=0.d0
      do 100 IIMAT=1,NMAT            !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=CVVOLM(ICVLA)+CVVOLM(ICVLB) 
        grdc(ICVLA,1,1)=grdc(ICVLA,1,1)+dum1
        grdc(ICVLB,1,1)=grdc(ICVLB,1,1)+dum1
        grdc(ICVLA,3,1)=grdc(ICVLA,3,1)+1.d0
        grdc(ICVLB,3,1)=grdc(ICVLB,3,1)+1.d0
        enddo
 100  enddo
!
      do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT<0) cycle
        do ICVL=ICVS,ICVE
        dum1=grdc(ICVL,1,1)-(grdc(ICVL,3,1)-1.d0)*CVVOLM(ICVL)
        if(dum1<0.d0) 
     &    call FFRABORT(1,'ERR: hat_bar error, call supportor')
        R_deltV(ICVL)=(dum1/CVVOLM(ICVL))**r1pn
        enddo
      enddo
!
!------------------
! --- bar(ui*uj)
!------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT<0) cycle
      do ICVL=ICVS,ICVE
      Lij(ICVL,1)=vel(ICVL,1)**2
      Lij(ICVL,2)=vel(ICVL,2)**2
      Lij(ICVL,3)=vel(ICVL,3)**2
      Lij(ICVL,4)=vel(ICVL,1)*vel(ICVL,2)
      Lij(ICVL,5)=vel(ICVL,1)*vel(ICVL,3)
      Lij(ICVL,6)=vel(ICVL,2)*vel(ICVL,3)
      enddo
      enddo
!-----------------------
! ---- Lij=hat(ui*uj)
!-----------------------
      do II=1,3
      mvel2(:,II)=Lij(:,II)
      enddo
      call dc_symprs(3,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &		     mat_cal,mvel2)
!
      if(fltr==1) then
        call hat_bar(3,2,dlalpha,R_deltV,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      else
        call top_hat(3,2,dlalpha,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      endif
!
      do II=1,3
      Lij(:,II)=mvel1(:,II)
      enddo
!
      mvel2(:,1:3)=Lij(:,4:6)
      call dc_symprs(3,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &		     mat_cal,mvel2)
      if(fltr==1) then
        call hat_bar(3,2,dlalpha,R_deltV,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      else
        call top_hat(3,2,dlalpha,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      endif
      Lij(:,4:6)=mvel1(:,1:3)
!----------------
! --- hat(ui)
!----------------!
      if(fltr==1) then
        call hat_bar(3,2,dlalpha,R_deltV,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       vel,grdc,mvel1)
      else
        call top_hat(3,2,dlalpha,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      endif
!------------------------------------------------
! --- Mij=-dlalpha*abs(hat(Sij))*hat(Sij)
!------------------------------------------------
      call grad_cell(3,2,
     &	MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,mvel1,grdc)
!
      
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      Sij(1)=grdc(ICVL,1,1)
      Sij(2)=grdc(ICVL,2,2)
      Sij(3)=grdc(ICVL,3,3)
      Sij(4)=(grdc(ICVL,1,2)+grdc(ICVL,2,1))
      Sij(5)=(grdc(ICVL,1,3)+grdc(ICVL,3,1))
      Sij(6)=(grdc(ICVL,2,3)+grdc(ICVL,3,2))
      Sij(7)=R_deltV(ICVL)**2*dsqrt(2.d0*(Sij(1)**2			!S
     &			+Sij(2)**2
     &			+Sij(3)**2)
     &		 +(Sij(4)**2+Sij(5)**2+Sij(6)**2))
      Mij(ICVL,1)=-Sij(7)*Sij(1)		!S*Sij
      Mij(ICVL,2)=-Sij(7)*Sij(2)
      Mij(ICVL,3)=-Sij(7)*Sij(3)
      Mij(ICVL,4)=-Sij(7)*Sij(4)
      Mij(ICVL,5)=-Sij(7)*Sij(5)
      Mij(ICVL,6)=-Sij(7)*Sij(6)
      enddo
      enddo      
!------------------------------------------------
! --- Lij= Lij-hat(ui)hat(uj)
!------------------------------------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      Lij(ICVL,1)=Lij(ICVL,1)-mvel1(ICVL,1)**2
      Lij(ICVL,2)=Lij(ICVL,2)-mvel1(ICVL,2)**2
      Lij(ICVL,3)=Lij(ICVL,3)-mvel1(ICVL,3)**2
      Lij(ICVL,4)=Lij(ICVL,4)-mvel1(ICVL,1)*mvel1(ICVL,2)
      Lij(ICVL,5)=Lij(ICVL,5)-mvel1(ICVL,1)*mvel1(ICVL,3)
      Lij(ICVL,6)=Lij(ICVL,6)-mvel1(ICVL,2)*mvel1(ICVL,3)
      enddo
      enddo

!----------------
! --- Sij
!----------------
      call grad_cell(3,2,
     &	MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vel,grdc)
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      Sij(1)=grdc(ICVL,1,1)
      Sij(2)=grdc(ICVL,2,2)
      Sij(3)=grdc(ICVL,3,3)
      Sij(4)=(grdc(ICVL,1,2)+grdc(ICVL,2,1))
      Sij(5)=(grdc(ICVL,1,3)+grdc(ICVL,3,1))
      Sij(6)=(grdc(ICVL,2,3)+grdc(ICVL,3,2))
      Sij(7)=dsqrt(2.d0*(Sij(1)**2			!S
     &			+Sij(2)**2
     &			+Sij(3)**2)
     &		 +(Sij(4)**2+Sij(5)**2+Sij(6)**2))
      abs_dspd(ICVL)=Sij(7)
      mvel3(ICVL,1)=Sij(7)*Sij(1)		!S*Sij
      mvel3(ICVL,2)=Sij(7)*Sij(2)
      mvel3(ICVL,3)=Sij(7)*Sij(3)
      mvel4(ICVL,1)=Sij(7)*Sij(4)
      mvel4(ICVL,2)=Sij(7)*Sij(5)
      mvel4(ICVL,3)=Sij(7)*Sij(6)
      enddo
      enddo
!-------------------------
! --- Mij=hat(|Sij|*Sij)
!-------------------------
      mvel2(:,1:3)=mvel3(:,1:3)
      call dc_symprs(3,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &		     mat_cal,mvel2)
      if(fltr==1) then
      call hat_bar(3,2,dlalpha,R_deltV,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      else
        call top_hat(3,2,dlalpha,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      endif
      mvel3(:,1:3)=mvel1(:,1:3)
!
      mvel2(:,1:3)=mvel4(:,1:3)
      call dc_symprs(3,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &		     mat_cal,mvel2)
      if(fltr==1) then
      call hat_bar(3,2,dlalpha,R_deltV,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &        mvel2,grdc,mvel1)
      else
        call top_hat(3,2,dlalpha,
     &	       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &	       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &         LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &	       mvel2,grdc,mvel1)
      endif
      mvel4(:,1:3)=mvel1(:,1:3)
!---------------------
! --- Mij=hat(|Sij|*Sij)+Mij
!---------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      Mij(ICVL,1)=Mij(ICVL,1)+mvel3(ICVL,1) 
      Mij(ICVL,2)=Mij(ICVL,2)+mvel3(ICVL,2) 
      Mij(ICVL,3)=Mij(ICVL,3)+mvel3(ICVL,3) 
      Mij(ICVL,4)=Mij(ICVL,4)+mvel4(ICVL,1) 
      Mij(ICVL,5)=Mij(ICVL,5)+mvel4(ICVL,2) 
      Mij(ICVL,6)=Mij(ICVL,6)+mvel4(ICVL,3) 
      enddo
      enddo      
!
      mvel1(:,1)=0.d0
      mvel2(:,1)=0.d0
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      a1=   Lij(ICVL,1)*Mij(ICVL,1)
     &     +Lij(ICVL,2)*Mij(ICVL,2)
     &	   +Lij(ICVL,3)*Mij(ICVL,3)
     &     +2.d0*Lij(ICVL,4)*Mij(ICVL,4)
     &	   +2.d0*Lij(ICVL,5)*Mij(ICVL,5)
     &	   +2.d0*Lij(ICVL,6)*Mij(ICVL,6)
      mvel1(ICVL,1)=a1
      b1=   Mij(ICVL,1)*Mij(ICVL,1)
     &	   +Mij(ICVL,2)*Mij(ICVL,2)
     &	   +Mij(ICVL,3)*Mij(ICVL,3)
     &	   +2.d0*Mij(ICVL,4)*Mij(ICVL,4)
     &	   +2.d0*Mij(ICVL,5)*Mij(ICVL,5)
     &	   +2.d0*Mij(ICVL,6)*Mij(ICVL,6)
      mvel2(ICVL,1)=b1
      enddo
      enddo
!----------------------
! --- space averaging
!----------------------
      mvel3(:,:)=0.d0
      mvel4(:,:)=0.d0
      do IIMAT=1,NMAT 
      if(.not.mat_cal(IIMAT)) cycle 
      ICFS=MAT_CFIDX(IIMAT-1)+1 
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=CVVOLM(ICVLA)*mvel1(ICVLA,1)
     &      +CVVOLM(ICVLB)*mvel1(ICVLB,1)
        dum2=CVVOLM(ICVLA)*mvel2(ICVLA,1)
     &      +CVVOLM(ICVLB)*mvel2(ICVLB,1)
        dum3=CVVOLM(ICVLA)+CVVOLM(ICVLB)
        mvel3(ICVLA,1)=mvel3(ICVLA,1)+dum1
        mvel3(ICVLB,1)=mvel3(ICVLB,1)+dum1
        mvel3(ICVLA,2)=mvel3(ICVLA,2)+dum2
        mvel3(ICVLB,2)=mvel3(ICVLB,2)+dum2
        mvel3(ICVLA,3)=mvel3(ICVLA,3)+1.d0
        mvel3(ICVLB,3)=mvel3(ICVLB,3)+1.d0
        mvel4(ICVLA,1)=mvel4(ICVLA,1)+dum3
        mvel4(ICVLB,1)=mvel4(ICVLB,1)+dum3
      enddo
      enddo
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT<0) cycle
      do ICVL=ICVS,ICVE
      dum4=(mvel3(ICVL,3)-1.D0)
      DUM1=mvel3(ICVL,1)-mvel1(ICVL,1)*CVVOLM(ICVL)*dum4
      DUM2=mvel3(ICVL,2)-mvel2(ICVL,1)*CVVOLM(ICVL)*dum4
      DUM3=mvel4(ICVL,1)-CVVOLM(ICVL)*dum4

      mvel1(ICVL,1)=DUM1/DUM3
!     &  mvel3(ICVL,1)-mvel1(ICVL,1)*CVVOLM(ICVL)*(mvel3(ICVL,3)-1.D0)
      mvel2(ICVL,1)=DUM2/DUM3
!     &  mvel3(ICVL,2)-mvel2(ICVL,1)*CVVOLM(ICVL)*(mvel3(ICVL,3)-1.D0)
      
      enddo
      enddo

!--------------------
! --- (Cs*deltV)**2
!--------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      dum1=mvel1(ICVL,1)        !MijLij
      dum2=mvel2(ICVL,1)        !MijLij
      csdles=dum1/(dum2+SML)
      csdles=max(0.d0,csdles)
      dum1=rho(ICVL)*csdles*abs_dspd(ICVL)
      rmut(ICVL) = min(1.d0,dum1)

      FIELD_U(ICVL,1)=csdles/(CVVOLM(ICVL))**r1pn !Cs
      enddo
      enddo
      deallocate(mvel1,mvel2)
      deallocate(mvel3,mvel4)
!
      return
!
      end subroutine cal_Dynamic_SGS


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE cal_wall_dir(LVEDGE,CVCENT,MAT_NO,LBC_SSF,
     &  SFAREA,MAT_CV,MAT_CVEXT,MAT_CFIDX,mat_cal,t_dir)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_scalar,only  : ivof
      use module_vof     ,only: ical_vof,intervof,LG,LS,change_phase
      use module_metrix,only  : akstmp=>d1work1
      use module_boundary,only : kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdsld,kdfire,kdbuff,kxnone,
     &                           nbcnd,LBC_INDEX,MAT_BCIDX,kdintr
     &                           ,kdpres,kdstag,boundName
     &                           ,kdcvd,idis,LBC_pair,rotsld
      use module_gravity ,only : ggg
      implicit none
!
! --- [dummy arguments]
! 
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3,MXALLCV)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      INTEGER,INTENT(IN)    :: MAT_NO(       0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV(         MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      real*8 ,intent(inout) :: t_dir     (MXALLCV,3)
      logical :: lkdwall
!
! --- [local entities]
!
      integer :: IMAT,ICFS,ICFE,ICFL,ICV,ICVS,ICVE,nb,kd
      integer :: IBFS,IBFE,IBFL
      integer :: IIMAT,ICVA,ICVB,icv_dum( MXALLCV,6)
      integer :: i_min,i_max,j_min,j_max,k_min,k_max,i,l
      real*8  :: n1,n2,n3,g1,g2,g3,cos_ang,gggg
      real*8  :: a11,a12,a13,a21,a22,a23,a31,a32,a33
      real*8  :: det_A
!
      t_dir=0.d0
      gggg=(ggg(1)*ggg(1)+ggg(2)*ggg(2)+ggg(3)*ggg(3))**0.5
      do 300 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 300
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 300
      kd=kdbcnd(0,nb)
      lkdwall=kd==kxnone.or.
     &        kd==kdbuff.or.
     &        kd==kdintr.or.
     &        kd==kdfire.or.
     &        kd==kdcvd
      if(lkdwall) then
      do IBFL=IBFS,IBFE                      
      ICFL=LBC_SSF(IBFL)
      ICVA=LVEDGE(1,ICFL)
      n1=-SFAREA(1,ICFL)
      n2=-SFAREA(2,ICFL)
      n3=-SFAREA(3,ICFL)
      g1=ggg(1)/gggg
      g2=ggg(2)/gggg
      g3=ggg(3)/gggg
      cos_ang=n1*g1+n2*g2+n3*g3
      cos_ang=(1.d0-cos_ang*cos_ang)**0.5
      a11=n2*g3-n3*g2
      a12=n3*g1-n1*g3
      a13=n1*g2-n2*g1
      a21=n1
      a22=n2
      a23=n3
      a31=g1
      a32=g2
      a33=g3
      det_A=a11*a22*a33+a21*a32*a13+a31*a12*a23
     &     -a11*a32*a23-a31*a22*a13-a21*a12*a33
      if(abs(det_A).lt.1.d-10) then
      t_dir(ICVA,1)=0.d0
      t_dir(ICVA,2)=0.d0
      t_dir(ICVA,3)=0.d0
      else
      t_dir(ICVA,1)=1.d0/det_A*(a12*a23-a13*a22)*cos_ang     
      t_dir(ICVA,2)=1.d0/det_A*(a13*a21-a11*a23)*cos_ang  
      t_dir(ICVA,3)=1.d0/det_A*(a11*a22-a12*a21)*cos_ang  
      end if
      enddo
      end if
  300 continue
      return     
      end subroutine cal_wall_dir
!


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine flamelet_ys_ORG
     &  (MAT_CVEXT,MAT_CAL,MAT_NO,aks,rho,tmp,yys)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_scalar,only  :  calgeq,calxi,caltemp,calh 
     &                          ,igeq,ixi,itemp,ihflmlt    
     &                          ,fd_rho,fd_tmp
     &                          ,nrho,ntmp
      use module_scalar,only  :  XXi0,
     &                           Xfuel_rho,
     &                           Xfuel_temp,
     &                           XOxidant_rho,
     &                           XOxidant_temp,
     &                           XXi_max
      use module_scalar,only  :  POLYC_YS,fitVALUE,mm,flag_flm_ys
      use module_model,only   :  ical_vect
      use module_material,only:  BGCOMP
      use module_constant
      use module_species,only : spcnam
      use module_scalar,only  : tags,NVAR_ORG,ORGtable,NVAR_COMP,IZ
      use module_io,      only : ifll,ifle
!
!
!
      implicit none
!
      real*8 ,intent(in)    :: aks  (MXALLCVR,MXRANS,2)
      real*8 ,intent(out)   :: rho  (MXALLCV ,2)
      real*8 ,intent(out)   :: tmp  (MXALLCV ,2)
      real*8 ,intent(inout) :: yys  (MXALLCV,mxcomp)
!
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
!      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      LOGICAL,INTENT(IN)    :: MAT_CAL   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(       0:MXMAT)
!
      real*8 :: fdmp,cfdmp
      real*8 :: gggeq
      real*8 :: xi
      real*8 :: dum
!
      real*8 :: rho_u,rho_b,rho_f,rho_o,sum_rho,dumXXi,
     &          dum3,dum2,dum1
      real*8 :: tmp_u,tmp_b,tmp_f,tmp_o,sum_tmp,dum_y1,dum_y2
!

      integer :: IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,IDCL
      integer :: i,j,itmp,irho,IMAT,ICOM
      integer,save :: flag=0
!
      if(.NOT.flag_flm_ys) return 
!
      if(flag==0) then
        do ICOM=NVAR_ORG,NCOMP,-1 
        if(TRIM(adjustl(tags(ICOM)))/=spcnam(ICOM-NCOMP+1))
     &  then
          write(ifle,'(1X,a,I4,2a)') 
     &     'species no=',ICOM,TRIM(adjustl(tags(ICOM))),
     &      trim(spcnam(ICOM-NCOMP+1))
          write(ifle,'(1X,2a)') 
     %     'species no is different between fflow.ctl and',
     &      '&flamelet/flamelet_file'
          call FFRABORT(1,'MSG: reset species no in fflow.ctl')
        endif
        enddo
        flag=1
      endif
!-----------------------------
! --- BackGrand gas
!-----------------------------
      do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          yys(ICVS:ICVE,1:NCOMP)=0.d0
          cycle
        endif
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
!
        do ICOM=1,NCOMP
        do ICVL=ICVS,ICVE
        xi=aks(ICVL,ixi,1)
        dum2=0.d0
        do J=2,IZ
        dum1=ORGtable(J,1)
        
        if(dum1>=xi) then
          dum2=ORGtable(J-1,1)
          dum_y1=ORGtable(J,NVAR_COMP+icom)
          dum_y2=ORGtable(J-1,NVAR_COMP+icom)
          dum3=(dum_y1-dum_y2)/(dum1-dum2)*(xi-dum2)+dum_y2
        exit
        endif
        enddo
        yys(ICVL,ICOM)=dum3 
        enddo
!
!        dum3=0.d0
!        do 100 ICOM=1,ncomp
!        dum3=dum3+yys(ICVL,ICOM)
! 100    enddo
!        yys(ICVL,:)=yys(ICVL,:)/dum3
!
        enddo
!
      enddo
!
      return
!
      end subroutine flamelet_ys_ORG

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cal_canopy(iter,deltt,time,
     &    MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,MAT_INDEX,
     &    rmu,rmut)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_metrix,only  : multiR
      USE module_dimension
      use module_hpcutil
      use module_io,     only : ifli,ifll,ifle,gdScale
!
!
      implicit none
!
!
! --- [dummy arguments]
!
      real*8 ,intent(in)  :: deltt,time
      integer,intent(in)  :: iter
      INTEGER,INTENT(IN)  :: MAT_NO(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_CV(   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_DCIDX(0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_INDEX(0:MXMAT)
!
      real*8 ,intent(in)  :: rmu   (   MXALLCV)
      real*8 ,intent(out) :: rmut  (   MXALLCV)      
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: idum
      

      do 100 IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT<0) cycle 
        DO ICVL=ICVS,ICVE
        idum=multiR(ICVL)
        if(idum>=21.and.idum<=50) then 
          rmut(ICVL)=0.d0
        endif
        enddo
 100  enddo

      
      return
      end subroutine cal_canopy
