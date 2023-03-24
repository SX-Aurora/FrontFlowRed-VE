!
!      subroutine bc_hysdif() 
!      subroutine bc_ransdif()
!      subroutine bc_veldif() 
!      subroutine bc_intrfc() 
!      subroutine bc_intrfcimp() 
!      subroutine bc_kdbpv() 
!      subroutine bc_kdbp_cell 
!      subroutine bc_kdbty() 
!      subroutine bc_pp0()
!      subroutine bc_prs()
!      subroutine bc_prdmsk()
!      subroutine bc_setbnd()
!      subroutine bc_setbnd2()
!      subroutine utl_boxmlr(x1,x2) 
!      subroutine bc_rans() 
!      subroutine bc_alpha()
!      subroutine bc_rva()
!      subroutine bc_tys()
!      subroutine bc_vel()
!      subroutine bc_velcofd()
!      subroutine bc_wallke()
!
!--------------------------------------------------------------------- 
!
!      INLET :       rav>0  : outflow : Dirichlet > vel,Ys,T,aks
!                    rva<0  : inflow  : Dirichlet > vel,Ys,T,aks
!
!      OUTLET:       rav>0  : outflow : Neumann   > vel,Ys,T,aks
!                    rva<0  : inflow  : Dirichlet > vel,Ys,T,aks
!
!      OUTLET-OPEN : rav>0  : outflow : Neumann   > vel,Ys,T,aks
!                    rva<0  : inflow  : Dirichlet > Ys,T,aks
!                                     : Neumann   > vel
!---------------------------------------------------------------------
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_thmodif(ndim,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  dflxy)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : ktneum,kttrns,ktEWF,kyneum,kyEWF,
     &                           kytrns,kdintr,
     &                           kdcvd,kydirc,kdilet,kdolet,openout,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           kdshutr,kdbuff
      use module_species ,only : wm
!
! 1. Set boundary condition for diffusion flux of energy & mass
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim
      integer,intent(in)  :: LVEDGE    (2,  MXCVFAC)
      integer,intent(in)  :: LBC_SSF   (    MXSSFBC)
      logical,INTENT(IN)  :: mat_cal   (    0:MXMAT)
      real*8 ,intent(inout) :: dflxy   (    MXCVFAC,ndim)
      integer,intent(in)  :: LCYCSF    (    MXSSFBC)
!
! --  [local entities]
!
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,IBFS,IBFE,IBFL,ICFL,ICV,IDC,ICVP,IDCP
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ICFP,ICOM
      real*8  :: dum3,dum2,dum1

!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      kdt=kdbcnd(2,nb)
      kdy=kdbcnd(3,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      do 1100 IBFL=IBFS,IBFE
      ICFL=LBC_SSF(IBFL)
      dflxy(ICFL,:)=0.d0
 1100 enddo
 1000 enddo
!
      end subroutine bc_thmodif
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_hysdif(LVEDGE,LBC_SSF,LCYCSF,mat_cal,mat_no,
     &  sfarea,cp,utau,frstcv,
     &  vel,tmp,yys,tmpbnd,yysbnd,htcbnd,mtcbnd,HTFLUX,
     &  rho,rmu,aks,rmut,
     &  dflxt,dflxy,dsclt,dscly)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : ktneum,kttrns,ktEWF,kyneum,kyEWF,
     &                           kytrns,kdintr,
     &                           kdcvd,kydirc,kdilet,kdolet,openout,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           kdshutr,kdbuff
      use module_species ,only : wm
      use module_rad, only     : radflag,radmodel
      use module_radsxf, only  : RadHeat,RadHeatFlux
      use module_metrix,only   : SHUTFL
      use module_turbparm,only : yplsm,prturb,scturb,akappa,E_vel
      use module_rans    ,only : cmu
      use module_material ,only : prlmnr,sclmnr,nflud
      use module_scalar,only  : ike
      use module_model,only   : icaltb,ke,RNG,CHEN,KE2S,idrdp,comp,ke_low
      
!
! 1. Set boundary condition for diffusion flux of energy & mass
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: LVEDGE    (2,  MXCVFAC)
      integer,intent(in)  :: LBC_SSF   (    MXSSFBC)
      real*8 ,intent(in)  :: tmp       (    MXALLCV)
      real*8 ,intent(in)  :: yys       (    MXALLCV,MXCOMP)
      real*8 ,intent(in)  :: tmpbnd    (    MXSSFBC)
      real*8 ,intent(in)  :: yysbnd    (    MXSSFBC,MXCOMP)
      real*8 ,intent(inout)  :: htcbnd    (    MXSSFBC)
      real*8 ,intent(inout)  :: HTFLUX    (    MXSSFBC)
      real*8 ,intent(in)  :: mtcbnd    (    MXSSFBC)
      logical,INTENT(IN)  :: mat_cal   (    0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO    (   0:MXMAT)      
      real*8 ,intent(inout) :: dflxt   (    MXCVFAC)
      real*8 ,intent(inout) :: dflxy   (    MXCVFAC,MXCOMP)
      real*8 ,intent(out) :: dsclt     (    MXALLCV)
      real*8 ,intent(out) :: dscly     (    MXALLCV)
      integer,intent(in)  :: LCYCSF    (    MXSSFBC)
!
      real*8 ,intent(in) :: vel       (  MXALLCV,3)
      real*8 ,intent(in)  :: rho       (    MXALLCV)
      real*8 ,intent(in)  :: rmu       (    MXALLCV)
      real*8 ,intent(in)  :: aks       (    MXALLCVR,MXRANS)   
!
      real*8 ,intent(in) :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in) :: cp        (  MXALLCV)
      real*8 ,intent(in) :: utau     (0:MXSSFBC)
      real*8 ,intent(in) :: FRSTCV    (  MXSSFBC)
      real*8 ,intent(in) :: rmut      (   MXALLCV)
      
!
! --  [local entities]
!
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,IBFS,IBFE,IBFL,ICFL,ICV,IDC,ICVP,IDCP
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ICFP,ICOM
      real*8  :: ratio,dradt,dum3,dum2,dum1
      real*8  :: yplus,AAA,BBB,CCC,PPP,ux,uy,uz,up,Pr,yp,Uc

!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      kdt=kdbcnd(2,nb)
      kdy=kdbcnd(3,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdshutr.or.kd==kdbuff) goto 2000
      if(kdt.eq.ktneum) then
        do 1100 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        dflxt(ICFL)=tmpbnd(IBFL)
 1100   enddo
      elseif(kdt.eq.kttrns) then 
        IF(radflag.EQ.1) THEN
          do 1200 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
	  dradt=RadHeatFlux(IDC)/htcbnd(IBFL)
          dflxt(ICFL)=htcbnd(IBFL)*(tmp(IDC)-dradt-tmp(ICV))
          dsclt(IDC)=htcbnd(IBFL)
          HTFLUX(IBFL)=dflxt(ICFL)
 1200     enddo
        else
          do 1201 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dflxt(ICFL)=htcbnd(IBFL)*(tmp(IDC)-tmp(ICV))
          HTFLUX(IBFL)=dflxt(ICFL)
          dsclt(IDC)=htcbnd(IBFL)
 1201     enddo
        endif
      elseif(kdt==ktEWF) then
        if(icaltb/=ke.and.
     &     icaltb/=RNG.and.
     &     icaltb/=KE2S.and.
     &     icaltb/=ke_low.and.
     &     icaltb/=CHEN
     &   ) then
          call FFRABORT(1,'ERR: temp=ktEWF ONLY for KE,RNG and CHEN')
        endif 
!
        do IBFL=IBFS,IBFE 
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
!
        yp=FRSTCV(IBFL)
        ux=vel(ICV,1)-vel(IDC,1)
        uy=vel(ICV,2)-vel(IDC,2)
        uz=vel(ICV,3)-vel(IDC,3)
        up=ux*SFAREA(1,ICFL)
     &    +uy*SFAREA(2,ICFL)
     &    +uz*SFAREA(3,ICFL)
        ux=ux-up*SFAREA(1,ICFL)
        uy=uy-up*SFAREA(2,ICFL)
        uz=uz-up*SFAREA(3,ICFL)
        up=(ux*ux+uy*uy+uz*uz)+SML 
        dum1=(aks(ICV,ike(1)))
        dum2=cmu**(0.25d0)*dsqrt(aks(ICV,ike(1)))
        Pr=prlmnr(IMAT) 
        yplus=yp*rho(ICV)*dum2/(rmu(ICV))+SML
        if(yplus<yplsm*Pr) then
          dum1=min(up,dum1)
          !dum1=up**2 
          BBB=0.5d0*rho(ICV)*Pr*dum2*dum1
          if(idrdp/=comp) BBB=0.d0
          AAA=rho(ICV)*CP(ICV)*dum2
          dum1=((tmp(IDC)-tmp(ICV))*AAA-BBB)/(Pr*yplus)
          dflxt(ICFL)=dum1
          dsclt(IDC)=AAA/(Pr*yplus)
          htcbnd(IBFL)=AAA/(Pr*yplus)
        else
          dum1=Pr/prturb
          PPP=9.24d0*((dum1)**0.75d0-1.d0)
     &              *(1.d0+0.28d0*exp(-0.007d0*dum1))
          AAA=rho(ICV)*CP(ICV)*dum2
          BBB=prturb*(log(E_vel*yplus)/akappa+PPP)
          !!!!!          Uc=yplsm*Pr*utau(IBFL)
          Uc=log(E_vel*yplsm)/akappa*utau(IBFL)
          CCC=0.5d0*rho(ICV)*dum2*(prturb*up+(Pr-prturb)*Uc**2)
          if(idrdp/=comp) CCC=0.d0
          dflxt(ICFL)=((tmp(IDC)-tmp(ICV))*AAA-CCC)/BBB
          dsclt(IDC)=AAA/BBB
          htcbnd(IBFL)=AAA/BBB
        endif
        HTFLUX(IBFL)=dflxt(ICFL)
        enddo
!
        if(kd==kdintr) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICFP=LCYCSF(IBFL)
          if(rmu(ICV)<SML) then
            HTFLUX(IBFL)=-dflxt(ICFP)
          endif
          enddo
        endif
      else
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        HTFLUX(IBFL)=dflxt(ICFL)
        enddo
      endif
!
      if(kdy==kyneum) then
        do ICOM=1,ncomp
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        dflxy(ICFL,ICOM)=yysbnd(IBFL,ICOM)
        enddo
        enddo
      elseif(kdy.eq.kytrns) then
        do 102 ICOM=1,ncomp
        do 1400 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        dflxy(ICFL,ICOM)=mtcbnd(IBFL)*(yys(IDC,ICOM)-yys(ICV,ICOM))
        dscly(IDC)=mtcbnd(IBFL)
 1400   enddo
  102   enddo
      endif
!---------------------------------------
! --- Only 'kd.eq.kdintr' for inlet BC 
!---------------------------------------
      if(kd==kdilet) then
        do ICOM=1,ncomp
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        dflxy(ICFL,ICOM)=0.d0
        enddo
        enddo
      endif
!
      if(kd==kdolet.and.openout(nb)/=5) then   !7777
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          dflxt(ICFL)=0.d0
          dflxy(ICFL,:)=0.d0
          enddo
      endif

 2000 continue

      if(kd==kdintr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        dflxy(ICFL,:)=0.d0
        dflxy(ICFP,:)=0.d0
        enddo
      endif
      if(kd==kdshutr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        if(SHUTFL(IBFL)==0) then
          
        else
          dflxy(ICFL,:)=0.d0
          dflxy(ICFP,:)=0.d0
          dflxt(ICFL)=0.d0
          dflxt(ICFP)=0.d0
          dsclt(ICFL)=0.d0
          dsclt(ICFP)=0.d0
          dscly(ICFL)=0.d0
          dscly(ICFP)=0.d0
        endif
        enddo
      endif
!
      IF(radflag.NE.1) GOTO 1000
!
      if(kdt==ktneum) then !.or.kdt==ktEWF.or.kdt==kttrns) then
        ratio=-1.d0 
!RadHeatFlux(IDC) is QW  ::  radiation generated heat flux
!NOTE: Qr>0, means heat emit from wall 
!      qr is emitted to the whole field,
!      for isolate boundary, qf+qr = 0 ==> qf=-qr
!	   qf is the local conduction heat flux (wall -> fluid)
!      that means the radiation loss is conpensated by conduction from fluid to wall
!NOTE: Qr>0, means heat emit from wall (loss of wall)
!      qt=h*(Tw-Tf),   qt=qc-qr,   qc=h*(Tw`-Tf)
!      the absorbed radiation energy change the Tw -> Tw`, i.e., Tw`=Tw + qr/h
         do IBFL=IBFS,IBFE
         ICFL=LBC_SSF(IBFL)
         ICV=LVEDGE(1,ICFL)
         IDC=LVEDGE(2,ICFL)
         dflxt(ICFL)=dflxt(ICFL)+RadHeatFlux(IDC)*ratio
!RadHeatFlux(IDC) is QW  ::  radiation generated heat flux
!NOTE: here is minus not plus, because positive wall heat flux 
!      means radiation_heat_loss of wall, this should be compensated
!	   by heat transfer from fluid to wall
         enddo
      endif
 1000 enddo
!
      end subroutine bc_hysdif
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_intrfc(iphs,
     &  MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  SFAREA,CVCENT,CVVOLM,SFCENT,
     &  rho,tmp,aks,grdt,rmdc,htcbnd,radbnd,tmpbnd,HTFLUX,rmu,
     &  dflxh,dflxt,dsclt) 
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     This subroutine is for interface BC 
!
! --- [module arguments] 
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdintr,ktneum,kttrns,ktEWF,kdbuff,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           kdshutr
      use module_io,only       : ifll,ifle
      use module_Euler2ph,only : ieul2ph
      use module_scalar,  only : iaph
      use module_rad, only     : radflag,radmodel
      use module_radsxf, only  : RadHeat,RadHeatFlux
      use module_metrix,only  : SHUTFL      
!
! 1.  Set boundary condition for  interface boundary
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iphs
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3,MXALLCV)
      integer,intent(in)    :: MAT_NO    (  0:MXMAT)
      real*8 ,intent(inout)    :: tmp       (  MXALLCV)
      real*8 ,intent(in)    :: grdt      (  MXALLCV,3)
      real*8 ,intent(in)    :: rmdc      (  MXALLCV)
      real*8 ,intent(in)    :: htcbnd    (  MXSSFBC)
      real*8 ,intent(in)    :: tmpbnd    (   mxssfbc)
      real*8 ,intent(in)    :: radbnd    (  MXSSFBC)
      real*8 ,intent(inout) :: dflxh     (  MXCVFAC)
      real*8 ,intent(inout) :: dflxt     (  MXCVFAC)
      real*8 ,intent(out)   :: dsclt     (  MXALLCV)
      real*8 ,intent(in)    :: rho       (  MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)    :: aks       (  MXALLCVR,mxrans)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: HTFLUX    (   mxssfbc)
      real*8 ,intent(in)    :: rmu   (  MXALLCV)
!
! --- [local entities]
!
!      logical :: fldsld
      integer :: i,j,k,l,m,n
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IMAT,IIMAT,mcvf,IMAT1,IMAT2,mcv,ndc
      integer :: IIMAT2,lfac
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      real*8  :: dx,dy,dz,drs,htr,flx,ga(2),gb(2),
     &           dum1,dum2,dum3,dum4,dum5,dum6,dl,Tw,TC,TP,wi1,wi2,
     &           FS,FS1,FS2,Tmax,Tmin,heatfct
      integer,save :: intr_type=1
      real*8	   :: ratio
!
! --- interface BC: +nb=>IMAT1 //// -nb=>IMAT2 //// (IMAT1>IMAT2) 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT1=MAT_NO(IIMAT)
      IMAT2=MAT_NO(IIMAT2)
      kd=kdbcnd(0,nb)
      if(kd.eq.kdintr) then
        if(IMAT1.gt.IMAT2) then 
          heatfct=1.d0
        else
          heatfct=-1.d0
        endif
!
!        if(IMAT1.gt.IMAT2) then
        if((IMAT1>0.and.IMAT2<0).or.(IMAT1<0.and.IMAT2>0)) then
            intr_type=1
        elseif(IMAT1>0.and.IMAT2>0) then
            intr_type=2
        elseif(IMAT1<0.and.IMAT2<0) then
            intr_type=3
        endif
!        endif
!
        if(intr_type==1.or.intr_type==2) then
        else
          call FFRABORT(1,'MSG:NOT Support')
        endif
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do 1100 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        tmp(IDCP)=tmp(ICV)
        tmp(IDC)=tmp(ICVP)
!
        Fs=0.5d0*(CVVOLM(ICVP)**(1/3.d0)+CVVOLM(ICV)**(1/3.d0))
!-----------------------
! ICV  - fluid cell
! ICVP - solid cell
!-----------------------
        kdt=kdbcnd(2,nb)
        if(kdt==ktEWF) then
          flx=HTFLUX(IBFL)
        elseif(kdt.eq.ktneum) then 
!if(Neumann) heat-transfer fluid=>solid then +value
!
          dsclt(ICV)=0.d0
          dsclt(ICVP)=0.d0
          TC=tmp(ICV)
          TP=tmp(ICVP)
!
! --- 
!
          dum1=1.d0/
     &        (1.d0/(htcbnd(IBFL)+SML)
     &          +Fs/rmdc(ICVP)
     &          +Fs/rmdc(ICV)
     &          )*(TP-TC)
          flx=-tmpbnd(IBFL)*heatfct
!
        elseif(kdt.eq.kttrns) then 
         if(.true.) then
          dum5=0.d0
          if(rmu(ICV)>0) dum5=1.d0
          dsclt(ICV)=0.d0
          dsclt(ICVP)=0.d0
          TC=tmp(ICV)
          TP=tmp(ICVP)
          dum1=1.d0/htcbnd(IBFL)
          dum2=Fs/rmdc(ICV)
          dum3=Fs/rmdc(ICVP)
          dum4=(TP-TC)*(dum1+dum5*dum2+(1.d0-dum5)*dum3)
     &        /(dum1*dum2+dum2*dum3+dum1*dum3)
          dum6=(dum1+dum5*dum2+(1.d0-dum5)*dum3)
     &        /(dum1*dum2+dum2*dum3+dum1*dum3)
          flx=dum4
         else
          call FFRABORT
     &    (1,'ERR: [interface temp="Neumann" NOT supported ')
!          dsclt(ICV)=0.d0
!          dsclt(ICVP)=0.d0
          TC=tmp(ICV)
          TP=tmp(ICVP)
          htr=htcbnd(IBFL)
          dum1=rmdc(ICVP)*TP/Fs+htr*TC
          dum2=htr+rmdc(ICVP)/Fs
          Tw=dum1/dum2
          Tmax=max(TC,TP)
          Tmin=min(TC,TP)
          Tw=max(min(Tw,Tmax),Tmin)
          dum1=1.d0/
     &        (1.d0/(htcbnd(IBFL)+SML)
     &          +Fs/rmdc(ICVP)
     &          +Fs/rmdc(ICV)
     &          )*(TP-TC)
          flx=dum1
         endif
        else
          write(ifle,'(1X,2a)')
     &      'ERR: [interface] BC only for [temp="Neumann"] or',
     &      '[temp="transfer"]'
          CALL FFRABORT(1,'stop at bc_intrfc')
        endif
!       
        if(ieul2ph>0) then
          flx=flx*aks(ICV,iaph(iphs))
        endif
!
        dflxt(ICFL)= flx
!        dflxt(ICFP)=-flx
!        dsclt(ICV)=dum6
        dflxh(ICFL)=0.d0
        dflxh(ICFP)=0.d0
 1100   enddo
        IF(radflag.NE.1) cycle
!
!      IF(IMAT1.GT.IMAT2) THEN
!          do IBFL=IBFS,IBFE
!	  ICFP=LCYCSF(IBFL)
!	  ICVP=LVEDGE(1,ICFP)
!	  IDCP=LVEDGE(2,ICFP)
!	  ratio = 1.d0
!          dflxt(ICFP)=dflxt(ICFP)-RadHeatFlux(IDCP)*ratio
!          enddo
!      ELSE
!          ratio = 1.d0
!          IF(IMAT1.GT.IMAT2) ratio = 1.d0
!	  do IBFL=IBFS,IBFE
!	  ICFL=LBC_SSF(IBFL)
!	  ICV=LVEDGE(1,ICFL)
!	  IDC=LVEDGE(2,ICFL)
!          ratio = 0.d0
!          if(rmu(ICV)<0) ratio = 1.d0

!          ratio = 0.d0
!          if(rmu(ICV)>0) ratio = 1.d0

!          ratio = 0.d0
!          if(rmu(ICV)>0) ratio = -1.d0

!          ratio = 0.d0
!          if(rmu(ICV)<0) ratio = -1.d0

!          ratio = 0.d0
!          if(rmu(ICV)<0) ratio = -1.d0

!	  dflxt(ICFL)=dflxt(ICFL)-RadHeatFlux(IDC)*ratio
!	  enddo
!      ENDIF

      elseif(kd==kdbuff) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        dflxh(ICFL)=0.d0
        dflxh(ICFP)=0.d0
        dum1=0.5d0*(dflxt(ICFL)-dflxt(ICFP))
        dflxt(ICFL)=dum1
        dflxt(ICFP)=-dum1
        enddo
!------appended by jiang yuyan for radiation heat transfer
!interface or sliding interface, only meaningful for fluid-solid one
!the above "if(.not.mat_cal(IIMAT)) cycle" ensure that only fluid-side
!is considered
!if(min(IMAT1,IMAT2).GT.0.or.max(IMAT1,IMAT2).LT.0) cycle
        IF(radflag.NE.1) cycle
!
        IF(IMAT1.GT.IMAT2) THEN
          do IBFL=IBFS,IBFE
	  ICFP=LCYCSF(IBFL)
	  ICVP=LVEDGE(1,ICFP)
	  IDCP=LVEDGE(2,ICFP)
	  ratio = 1.d0
          dflxt(ICFP)=dflxt(ICFP)-RadHeatFlux(IDCP)*ratio
          enddo
	ELSE
	  do IBFL=IBFS,IBFE
	  ICFL=LBC_SSF(IBFL)
	  ICV=LVEDGE(1,ICFL)
	  IDC=LVEDGE(2,ICFL)
	  ratio = 1.d0
	  dflxt(ICFL)=dflxt(ICFL)-RadHeatFlux(IDC)*ratio
	  enddo
        ENDIF
      elseif(kd==kdshutr) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        if(SHUTFL(IBFL)==0) then
          dflxh(ICFL)=0.d0
          dflxh(ICFP)=0.d0
          dum1=0.5d0*(dflxt(ICFL)-dflxt(ICFP))
          dflxt(ICFL)=dum1
          dflxt(ICFP)=-dum1
        else
          dflxh(ICFL)=0.d0
          dflxh(ICFP)=0.d0
          dflxt(ICFL)=0.d0
          dflxt(ICFP)=0.d0
        endif
        enddo
!
      endif
 1000 enddo
!     
      return
!
      end subroutine bc_intrfc
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_intrfcimp
     & (NSSF,IQMXCV,
     &  MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,alx,htcbnd,IQ,ip,aa)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     This subroutine is for interface BC
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only    : kdintr,ktneum,kttrns,ktEWF,
     &                              kdbuff,kdshutr,
     &                              nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX
!
! 1.  Set boundary condition for interface boundary
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: NSSF,IQMXCV
      integer,intent(in)  :: LVEDGE (2,MXCVFAC)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      integer,intent(in)  :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      real*8 ,intent(in)  :: alx    (  MXALLCV)
      real*8 ,intent(in)  :: htcbnd (  NSSF)
      integer,intent(out) :: IQ     (  IQMXCV)
      integer,intent(out) :: ip     (  MXALLCV,MAXIE)
      real*8 ,intent(out) :: aa     (  MXALLCV,0:MAXIE)
!
! --- [local entities]
!
      logical :: fldsld
      real*8  :: htr,ga(2)
      integer :: iptrn
      integer :: IMAT,IIMAT,IMAT1,IMAT2,IIMAT2
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!-------------------------------------------
! --- 
!-------------------------------------------
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT1=MAT_NO(IIMAT)
      IMAT2=MAT_NO(IIMAT2)
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdintr) then
        do 1100 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        if(IMAT1.gt.IMAT2) then
          if(IMAT1.gt.0) then
            fldsld=.true.
          else
            fldsld=.false.
          endif
!
          ga(1)=alx(IDC)
          ga(2)=alx(IDCP)
!
          iptrn=1
          kdt=kdbcnd(2,nb)
          if(kdt.eq.ktEWF) then
          elseif(kdt.eq.ktneum) then
          elseif(kdt.eq.kttrns) then
            if(fldsld) then
              ga(1)=htcbnd(IBFL)
            else
              htr=htcbnd(IBFL)
              iptrn=2
            endif
          else
            write(*,*)
     &      '### program error -1- (bc_intrfcimp)'
            CALL FFRABORT(1,'bc_intrfcimp')
          endif
!
          if(iptrn.eq.1) then
            aa(IDC,0)= ga(1)+ga(2)
            aa(IDC,1)=-ga(1)
            aa(IDC,2)=-ga(2)
            ip(IDC,1)= ICV
            ip(IDC,2)= ICVP
            IQ(IDC)  = 2
            aa(IDCP,0)= ONE
            aa(IDCP,1)=-ONE
            ip(IDCP,1)= IDC
            IQ(IDCP)  = 1
          else
            aa(IDC,0)= ga(1)+htr
            aa(IDC,1)=-ga(1)
            aa(IDC,2)=-htr
            ip(IDC,1)= ICV
            ip(IDC,2)= IDCP
            IQ(IDC)  = 2
            aa(IDCP,0)= ga(2)+htr
            aa(IDCP,1)=-ga(2)
            aa(IDCP,2)=-htr
            ip(IDCP,1)= ICVP
            ip(IDCP,2)= IDC
            IQ(IDCP)  = 2
          endif
        endif
 1100   continue
      endif
 1000 continue
!
      end subroutine bc_intrfcimp
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_kdbpv(LBC_SSF,LCYCSF,mat_cal,MAT_NO,kdbp,kdbv)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------------------------------------
!        Inlet    Outlet   Symm    Cyc      Wall  sliding CV_face
!    (touch-inlet)
! kdbp     3(1)     1       2(1)    0        2(1)   -1       0
! kdbv     1        2       3       0        3      -1       0
!        pressure staganation
! kdbp     1        3
!---------------------------------------------------------------------
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil,only  : my_rank
      use module_boundary,only : kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdsld,kdfire,kdbuff,kxnone,
     &                           nbcnd,LBC_INDEX,MAT_BCIDX,kdintr
     &                           ,kdpres,kdstag,boundName,kdshutr
     &                           ,kdcvd,idis,LBC_pair,rotsld,openout
     &                           ,kdovst,masflg
      use module_material,only : lclsd,ical_sld,nofld
      use module_model,only    : idrdp,incomp,mach0,comp,ical_vect,PEFC,
     &                            ical_dens
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_scalar  ,only: ical_FC
      use module_FUEL  ,only  : No_Mem
      use module_metrix,only  : SHUTFL
!
! 1. Set boundary condition flag for pressure & velocity
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT) 
      INTEGER,INTENT(IN)    :: MAT_NO    (  0:MXMAT)
      integer,intent(out)   :: kdbp      (  MXCVFAC)
      integer,intent(out)   :: kdbv      (  MXCVFAC)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMAT2,IIMAT2,IMAT_U,IMAT_U2
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,IPR
      logical :: lkdwall
!
      kdbp=0
      kdbv=0
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      lkdwall=kd==kxnone.or.
!     &        (kd==kdbuff).or.
     &        kd==kdintr.or.
     &        kd==kdfire.or.
     &        kd==kdcvd
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdsymm) then
        if(lclsd(IIMAT)==0.or.lclsd(IIMAT)==2.or.lclsd(IIMAT)==3) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2   !zhang8(=2) 
          kdbv(ICFL)=3   !3
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2
          kdbv(ICFL)=3     !3 
          enddo
        endif
      elseif(kd.ne.kdprdc.and.kd/=kdbuff.and.kd/=kdshutr) then
        if(kd.eq.kdolet) then
          if(lclsd(IIMAT)==1) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=1   !???zhang1
            enddo
          elseif(lclsd(IIMAT)==0.or.lclsd(IIMAT)==2) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1 
            enddo
          elseif(lclsd(IIMAT)==3) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1
            enddo
          endif
          if(openout(nb)==5) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1 !1  E2P  !=1:outlet ;=2:wall
            enddo
          endif
          if(openout(nb)==6) then 
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1
            enddo
          endif
          if(openout(nb)==10) then 
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1
            enddo
          endif
          if(openout(nb)==7) then   !4444
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1   
            enddo
          endif
          if(openout(nb)==8) then   !4444
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=3
            enddo
          endif
          if(openout(nb)==9) then   !4444
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1   !1
            enddo
          endif
          
        elseif(kd.eq.kdpres) then
          do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=1
          enddo
        elseif(kd.eq.kdstag) then
          if(idrdp.eq.incomp.or.idrdp.eq.mach0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=1
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=1  !???=1
            enddo
          endif
! --- 
        elseif(kd==kdilet.or.kd==kdtchi) then
          if(lclsd(IIMAT)==0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3  !zhang8  =3  ???4
            enddo
          elseif(lclsd(IIMAT)==2) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3  !zhang8  =3
            enddo
          elseif(lclsd(IIMAT)==3) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3
            ENDDO
          endif
        elseif(kd==kdsld) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=0    !=2  zhang!!!  4444
          enddo
        elseif(kd==kdovst) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=999  !999
          enddo
!
!        elseif(kd==kdsld.and.rotsld(nb)==0) then  !ical_sld==2  
!          if(idis(nb)==0) then
!            do IBFL=IBFS,IBFE
!            ICFL=LBC_SSF(IBFL)
!            kdbp(ICFL)=0
!            enddo
!          elseif(idis(nb)==1) then
!            do IBFL=IBFS,IBFE
!            ICFL=LBC_SSF(IBFL)
!            kdbp(ICFL)=2
!            enddo
!          endif
!
        elseif(lkdwall) then
          if(lclsd(IIMAT)==0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2  !2
            enddo
          elseif(lclsd(IIMAT)==2) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2
            enddo
          elseif(lclsd(IIMAT)==3) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2
            enddo
          endif
          if(kd==kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            kdbp(ICFP)=2
            enddo
          endif
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2
          enddo
        endif
!--------------
! --- velocity
!--------------
        if(kd.eq.kdilet.or.kd.eq.kdtchi) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbv(ICFL)=1
          enddo
        elseif(kd==kdsld) then  !ical_sld==1
          if(ical_vect) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=0  !   0 !3   4444
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=0  !
            enddo
          endif
        elseif(kd==kdovst) then
           do IBFL=IBFS,IBFE
           ICFL=LBC_SSF(IBFL)
           kdbv(ICFL)=2
           enddo
!        elseif(kd==kdsld.and.rotsld(nb)==0) then !ical_sld==2
!          do IBFL=IBFS,IBFE
!          ICFL=LBC_SSF(IBFL)
!          kdbv(ICFL)=0
!          enddo
        elseif(kd.eq.kdpres) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbv(ICFL)=2
          enddo
        elseif(kd.eq.kdstag) then
          if(idrdp.eq.incomp.or.idrdp.eq.mach0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=2
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=1  !   2!   1
            enddo
          endif
        elseif(kd.eq.kdolet) then
          do 1600 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbv(ICFL)=2    !2
 1600     continue

          if(openout(nb)==6) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=2    !1
            enddo
          endif
          if(openout(nb)==8.and.masflg(nb)==1) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=1
            enddo
          endif
          if(openout(nb)==10) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=1
            enddo
          endif
!
!          if(openout(nb)==5) then
!            if(kdphs(mph)==kdphs_l) then
!              do IBFL=IBFS,IBFE
!              ICFL=LBC_SSF(IBFL)
!              kdbv(ICFL)=3   !wall
!              enddo
!            else
!              do IBFL=IBFS,IBFE
!              ICFL=LBC_SSF(IBFL)
!              kdbv(ICFL)=2   !outlet
!              enddo
!            endif
!          endif
!
        elseif(kd.eq.kdintr) then
          if(idis(nb)==0) then
            do 1950 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            kdbv(ICFL)=3
            kdbv(ICFP)=3
 1950       continue
          elseif(idis(nb)==1) then
            DO IPR=1,2
            if(IPR==1) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_pair(nb)
            else
              IBFS=LBC_pair(nb)+1
              IBFE=LBC_INDEX(nb)
            endif
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbv(ICFL)=3
            enddo
            ENDDO
          endif
        else
          do 1700 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbv(ICFL)=3
 1700     continue
        endif
      elseif(kd.eq.kdprdc) then
        IF(ical_vect) then
!CIDR VECTOR
          do IBFL=IBFS,IBFE
          kdbv(LCYCSF(IBFL))=999
          kdbp(LCYCSF(IBFL))=999
          enddo
        ELSE
        if(idis(nb)==0) then
          do 1800 IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          kdbv(ICFP)=999   !okabe
          kdbp(ICFP)=999
 1800     continue
        elseif(idis(nb)==1) then   !zhang5
          DO IPR=1,2
          if(IPR==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_pair(nb)
          else
            IBFS=LBC_pair(nb)+1
            IBFE=LBC_INDEX(nb)
          endif
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=1!999  !zhang???
          kdbv(ICFL)=1!999
          enddo
          ENDDO
        endif
        ENDIF
      elseif(kd.eq.kdbuff) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          kdbv(ICFP)=999
          kdbp(ICFP)=999
          kdbp(ICFL)=0
        enddo
!        if(ical_FC==PEFC) then
!          IIMAT2=MAT_BCIDX(nb,2)
!          IMAT=MAT_NO(IIMAT)
!          IMAT2=MAT_NO(IIMAT2)
!          IMAT_U=nofld(abs(IMAT))
!          IMAT_U2=nofld(abs(IMAT2))
!          if(IMAT_U==No_Mem.or.IMAT_U2==No_Mem) then 
!            do IBFL=IBFS,IBFE
!            ICFP=LCYCSF(IBFL)
!            ICFL=LBC_SSF(IBFL)
!            kdbp(ICFP)=2
!            kdbp(ICFL)=2
!            enddo
!          endif
!        endif
      elseif(kd.eq.kdshutr) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          if(SHUTFL(IBFL)==0) then
            kdbv(ICFP)=999
            kdbp(ICFP)=999
          else
            kdbv(ICFP)=3
            kdbp(ICFP)=2
          endif
        enddo
      endif
 1000 continue
!
      return
      end subroutine bc_kdbpv
!

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_kdbpv_2p(mph,LBC_SSF,LCYCSF,mat_cal,kdbv)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------------------------------------
!        Inlet    Outlet   Symm    Cyc      Wall  sliding CV_face
!    (touch-inlet)
! kdbp     3(1)     1       2(1)    0        2(1)   -1       0
! kdbv     1        2       3       0        3      -1       0
!        pressure staganation
! kdbp     1        3
!---------------------------------------------------------------------
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdsld,kdfire,kdbuff,kxnone,
     &                           nbcnd,LBC_INDEX,MAT_BCIDX,kdintr
     &                           ,kdpres,kdstag,boundName,kdshutr
     &                           ,kdcvd,idis,LBC_pair,rotsld,openout
     &                           ,kdovst
      use module_material,only : lclsd,ical_sld
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
!
! 1. Set boundary condition flag for pressure & velocity
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mph
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT) 
      integer,intent(out)   :: kdbv      (  MXCVFAC)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,IPR
      logical :: lkdwall
!
      do nb=1,nbcnd
      if(openout(nb)/=5) cycle
      kd=kdbcnd(0,nb)
      if(kd/=kdolet) cycle
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kdphs(mph)==kdphs_l) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbv(ICFL)=3  !wall     3
        enddo
      else
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbv(ICFL)=2   !outlet     2
        enddo
      endif
      enddo

      return
      end subroutine bc_kdbpv_2p



!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_kdbp(iphs,LBC_SSF,LCYCSF,mat_cal,kdbp)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------------------------------------
!        Inlet    Outlet   Symm    Cyc      Wall  sliding CV_face
!    (touch-inlet)
! kdbp     3(1)     1       2(1)    0        2(1)   -1       0
! kdbv     1        2       3       0        3      -1       0
!        pressure staganation
! kdbp     1        3
!---------------------------------------------------------------------
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil,only  : my_rank
      use module_boundary,only : kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdsld,kdfire,kdbuff,kxnone,
     &                           nbcnd,LBC_INDEX,MAT_BCIDX,kdintr
     &                           ,kdpres,kdstag,boundName,kdshutr
     &                           ,kdcvd,idis,LBC_pair,rotsld,openout
     &                           ,kdovst
      use module_material,only : lclsd,ical_sld
      use module_model,only    : idrdp,incomp,mach0,comp,ical_vect
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_metrix,only  : SHUTFL
!
! 1. Set boundary condition flag for pressure & velocity
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iphs
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT) 
      integer,intent(out)   :: kdbp      (  MXCVFAC)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,IPR
      logical :: lkdwall
!
      kdbp=0
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      lkdwall=kd==kxnone.or.
!     &        kd==kdbuff.or.
     &        kd==kdintr.or.
     &        kd==kdfire.or.
     &        kd==kdcvd
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdsymm) then
        if(lclsd(IIMAT)==0.or.lclsd(IIMAT)==2.or.lclsd(IIMAT)==3) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2   !zhang8(=2)
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2
          enddo
        endif
      elseif(kd.ne.kdprdc.and.kd/=kdbuff.and.kd/=kdshutr) then
        if(kd.eq.kdolet) then
          if(lclsd(IIMAT)==1) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=1   !???zhang1
            enddo
          elseif(lclsd(IIMAT)==0.or.lclsd(IIMAT)==2) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1 
            enddo
          elseif(lclsd(IIMAT)==3) then
            do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=1
            enddo
          endif
          if(openout(nb)==5) then
            if(kdphs(iphs)==kdphs_l) then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              kdbp(ICFL)=2  !1  E2P: 1 for dp
              enddo
            endif
          endif
        elseif(kd.eq.kdpres) then
          do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=1
          enddo
        elseif(kd.eq.kdstag) then
          if(idrdp.eq.incomp.or.idrdp.eq.mach0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=1 
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3
            enddo
          endif
        elseif(kd==kdilet.or.kd==kdtchi) then
          if(lclsd(IIMAT)==0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3  !zhang8  =3  ???4
            enddo
          elseif(lclsd(IIMAT)==2) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3  !zhang8  =3
            enddo
          elseif(lclsd(IIMAT)==3) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=3
            ENDDO
          endif
        elseif(kd==kdovst) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2
          enddo
        elseif(kd==kdsld) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2    !=2!4  zhang!!!
          enddo
        elseif(lkdwall) then
          if(lclsd(IIMAT)==0) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2  !2
            enddo
          elseif(lclsd(IIMAT)==2) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2
            enddo
          elseif(lclsd(IIMAT)==3) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbp(ICFL)=2
            enddo
          endif
          if(kd==kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            kdbp(ICFP)=2
            enddo
          endif
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=2
          enddo
        endif
!
      elseif(kd.eq.kdprdc) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFL)=1!999  !zhang???
          enddo
      elseif(kd.eq.kdbuff) then
          stop 8877
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          kdbp(ICFP)=999  !zhang???
          enddo
      elseif(kd.eq.kdshutr) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          if(SHUTFL(IBFL)==0) then
            kdbp(ICFP)=999  !zhang???
          else
            kdbp(ICFP)=999 
          endif
          enddo
!      endif
      endif
 1000 continue
!
      return
      end subroutine bc_kdbp

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_kdbp_cell(LBC_SSF,KDBF)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdbcnd,nbcnd,kdolet,LBC_INDEX,kdilet
     &                          ,kdpres,kdstag,openout,kdsld
     &                          ,kdovst
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: LBC_SSF(MXSSFBC)
      integer,intent(out) :: KDBF   (MXCVFAC)
!
! --- [local entities]
!
      integer :: nb,kd
      integer :: IBFS,IBFE,IBFL,ICFL,IIMAT,IMAT
!
! --- pressure BC
!
      KDBF=0
      do 1000 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if(kd.eq.kdolet.and.openout(nb)/=2) then  !zhang-cvd
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        KDBF(ICFL)=1
        enddo
      elseif(kd.eq.kdolet.and.openout(nb)==2) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        KDBF(ICFL)=0
        enddo
      elseif(kd.eq.kdilet) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        KDBF(ICFL)=3
        enddo
      elseif(kd.eq.kdpres.or.kd.eq.kdstag) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        KDBF(ICFL)=2
        enddo
!      elseif(kd.eq.kdsld) then
!        IBFS=LBC_INDEX(nb-1)+1
!        IBFE=LBC_INDEX(nb)
!        do IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        KDBF(ICFL)=4
!        enddo
      endif
 1000 continue
!
      return
      end subroutine bc_kdbp_cell
!

!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_kdbty(LBC_SSF,LCYCSF,mat_cal,kdbt,kdby)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------------------------------------
!        Direc    neum   Symm    Cyc      Wall  sliding CV_face
!
! kdbt     1        2      2
! kdby     1        2      2
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_boundary,only : kdprdc,kdsymm,kdintr,ktneum,
     &                             ktEWF,kyEWF,
     &                           kyneum,
     &                           kdsld,idis,LBC_pair,rotsld,kdbuff,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           kdshutr
     &                           ,kdovst,kdolet,openout
      use module_metrix,only   : SHUTFL
      use module_model,only    : ical_vect
!
! 1. Set boundary condition flag for temp. & mass frac.
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      integer,intent(out) :: kdbt   (  MXCVFAC)
      integer,intent(out) :: kdby   (  MXCVFAC)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
! --- 
!
      kdbt=0
      kdby=0
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdsymm) then
        do 1100 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbt(ICFL)=2   !0
        kdby(ICFL)=2   !0
 1100   continue
      elseif(kd.eq.kdintr) then
          do 1300 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          kdbt(ICFP)=0   !2  !0
          kdbt(ICFL)=0   !2  !0
 1300     continue
      elseif(kd.ne.kdprdc.and.kd.ne.kdsld.and.kd/=kdovst
     &  .and.(kd/=kdbuff.and.kd/=kdshutr)) then
! --- kdbt
         if(kdbcnd(2,nb).eq.ktneum) then
            do 1200 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbt(ICFL)=2
 1200       continue
         else
            do 1400 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbt(ICFL)=1   !1 zhang
 1400       continue
         endif
! --- kdby
        if(kdbcnd(3,nb).eq.kyneum) then
          do 1500 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdby(ICFL)=2
 1500     continue
          if (kd.eq.kdintr) then
            do 1550 IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            kdby(ICFP)=2
 1550       continue
          endif
        else
          do 1600 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdby(ICFL)=1
 1600     continue
          if (kd.eq.kdintr) then
            do 1560 IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            kdby(ICFP)=1
 1560       continue
          endif
        endif
      elseif(kd==kdprdc) then
        if(idis(nb)==0) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          kdbt(ICFP)=999
          kdby(ICFP)=999
          enddo
        elseif(idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbt(ICFL)=999 !!!!!hui
          kdby(ICFL)=999 !!!!!hui
          enddo
        endif
      elseif(kd==kdbuff) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          kdbt(ICFP)=0   !999 !!!!!hui
          kdby(ICFP)=0   !999 !!!!!hui
        enddo
      elseif(kd==kdshutr) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          if(SHUTFL(IBFL)==0) then
            kdbt(ICFP)=0   !999 !!!!!hui
            kdby(ICFP)=0   !999 !!!!!hui
          else
            kdbt(ICFL)=999 !!!!!hui
            kdby(ICFL)=999 !!!!!hui
            kdbt(ICFP)=999 !!!!!hui
            kdby(ICFP)=999 !!!!!hui
          endif
        enddo
      elseif(kd==kdsld) then
        if(ical_vect) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          kdbt(ICFL)=0 !999!-1  !
          kdby(ICFL)=0
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          kdbt(ICFL)=0 !999!-1
          kdby(ICFL)=0
          enddo
        endif
      elseif(kd==kdovst) then
         do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbt(ICFL)=999
          kdby(ICFL)=999
          enddo
      endif
      if(kd==kdolet.and.(openout(nb)==10.or.
     &   openout(nb)==6.or.openout(nb)==8)) then
           do IBFL=IBFS,IBFE
           ICFL=LBC_SSF(IBFL)
           kdbt(ICFL)=1
           kdby(ICFL)=1
           enddo
      endif
 1000 continue
!
      return
      end subroutine bc_kdbty
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_kdbt_rad(LBC_SSF,LCYCSF,mat_cal,kdbt)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------------------------------------
!        Direc    neum   Symm    Cyc      Wall  sliding CV_face
!
! kdbt     1        2      2
! kdby     1        2      2
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdprdc,kdsymm,kdintr,ktneum,
     &                            ktEWF,kyEWF,
     &                           kyneum,
     &                           kdsld,idis,LBC_pair,rotsld,kdbuff,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           kdshutr
     &                           ,kdovst
      use module_metrix,only  : SHUTFL
!
! 1. Set boundary condition flag for temp. & mass frac.
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      integer,intent(out) :: kdbt   (  MXCVFAC)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
! --- 
!
      kdbt=0
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdsymm) then
        do 1100 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbt(ICFL)=2   !0
 1100   continue
      elseif(kd.eq.kdintr) then
          do 1300 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          kdbt(ICFP)=0
          kdbt(ICFL)=0
 1300     continue
      elseif(kd.ne.kdprdc.and.kd.ne.kdsld.and.kd/=kdovst
     &  .and.kd/=kdbuff.and.kd/=kdshutr) then !???
! --- kdbt
         if(kdbcnd(2,nb).eq.ktneum) then
            do 1200 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbt(ICFL)=2   !2
 1200       continue
         else
            do 1400 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbt(ICFL)=1
 1400       continue
         endif
      elseif(kd==kdprdc) then
        if(idis(nb)==0) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          kdbt(ICFP)=999
          enddo
        elseif(idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbt(ICFL)=999 !!!!!hui
          enddo
        endif
      elseif(kd==kdbuff) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          kdbt(ICFP)=0   !999 !!!!!hui
        enddo
      elseif(kd==kdshutr) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          if(SHUTFL(IBFL)==0) then
            kdbt(ICFP)=0  !999 !!!!!hui
          else
            kdbt(ICFP)=999
            kdbt(ICFL)=999
          endif
        enddo
      elseif(kd==kdsld) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        kdbt(ICFP)=0 !999!-1
        enddo
      elseif(kd==kdovst) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbt(ICFL)=999 !999!-1
        enddo
      endif
 1000 continue
!
      return
      end subroutine bc_kdbt_rad

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_pp0
     &  (LVEDGE,mat_cal,
     &   MAT_CV,MAT_CVEXT,MAT_NO,LBC_SSF,
     &   prsbnd,pp0)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!! BUG in this subroutine
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdolet,openout,wdbcnd,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX
!
! 1. Set thermodynamic pressure in case of zero Mach no. approx.
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      INTEGER,INTENT(IN)    :: MAT_CV(   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      real*8 ,intent(in)    :: prsbnd (  MXSSFBC)
      real*8 ,intent(inout) :: pp0    (  MXALLCV)
!
! --- [local entities]
!
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ierr1=0
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFS,ICFE
      real*8,allocatable  :: pp0x(:)
!
      ALLOCATE(pp0x(MXMAT),stat=ierr1)
!
      do 100 IIMAT=1,NMAT  !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) goto 100
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT.gt.0) then
        do ICVL=ICVS,ICVE
        pp0x(IIMAT)=pp0(ICVL)
        enddo
      endif
 100  continue
!
      do 200 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 200
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdolet.and.openout(nb)/=2) then !zhang-cvd
        do IBFL=IBFS,IBFE  
        if(IMAT.gt.0) then
          pp0x(IIMAT)=prsbnd(IBFL)
        endif
        enddo
      elseif(kd.eq.kdolet.and.openout(nb)==2) then
        if(IMAT.gt.0) then
          if(pp0x(IIMAT)>wdbcnd(1,nb)) then
            pp0x(IIMAT)=wdbcnd(1,nb)
          else
            do IBFL=IBFS,IBFE  
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            pp0x(IIMAT)=pp0(ICV)  !not finished
            enddo
          endif
        endif
      endif
 200  continue
!
      do 300 IIMAT=1,NMAT   !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) goto 300
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(IMAT.gt.0) then
        do ICVL=ICVS,ICVE
        pp0(ICVL)=pp0x(IIMAT)
      enddo
      endif
 300  continue
!
      deALLOCATE(pp0x)
!
      return
!
      end subroutine bc_pp0
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_prdmsk
     &          (MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &           LVEDGE,LBC_SSF,LCYCSF,msk)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdprdc,kdsld,kdintr,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           idis,LBC_pair,rotsld,kdbuff,kdshutr
     &                           ,kdovst
      use module_hpcutil,only  : my_rank
      use module_metrix,only   : SHUTFL
      use module_metrix,only   : tmpsld
      use module_model ,only   : ical_vect
!
! 1. Index of real cell for dummy cell on periodic boundary
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: LVEDGE (2,MXCVFAC)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)  :: MAT_NO(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_DCIDX(0:MXMAT)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      integer,intent(out) :: msk   (  MXALLCV)
!
! --- [local entities]
!
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,IPR
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFS,ICFE
!
      msk(:)=0
!
      do 100 IIMAT=1,NMAT    !ICV=1,NALLCV
      if(.not.mat_cal(IIMAT)) goto 100 
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      msk(ICVL)=ICVL
      enddo
      IDCS=MAT_DCIDX(IIMAT-1)+1
      IDCE=MAT_DCIDX(IIMAT)
      do ICVL=IDCS,IDCE
      msk(ICVL)=ICVL
      enddo
  100 continue
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      if(kd==kdprdc.or.kd==kdintr.or.kd==kdbuff) then !2222
        if(idis(nb)==0) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do 1200 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          msk(IDC)=ICVP
          msk(IDCP)=ICV
 1200     enddo
        elseif(idis(nb)==1) then
          DO IPR=1,2
          if(IPR==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_pair(nb)
          else
            IBFS=LBC_pair(nb)+1
            IBFE=LBC_INDEX(nb)
          endif
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          ICV=LVEDGE(1,ICFL)
          msk(IDC)=ICVP  !ICVP
          enddo
          ENDDO
        endif
      
      elseif(kd==kdsld) then

        if(ical_vect) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          tmpsld(IBFL,1)=ICVP
          enddo

          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          msk(IDC)=INT(tmpsld(IBFL,1))
          enddo
        else
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          msk(IDC)=ICVP
          enddo
        endif
      elseif(kd==kdovst) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
          ICVP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          msk(IDC)=ICVP
        enddo
      elseif(kd==kdshutr) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        if(SHUTFL(IBFL)==0) then
          msk(IDC)=ICVP
          msk(IDCP)=ICV
        else
          msk(IDC)=ICV
          msk(IDCP)=ICVP
        endif
        enddo
      endif
 1000 ENDDO
!
      return
      end subroutine bc_prdmsk
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_prs
     & (nflud,ipfix,deltt,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,FRSTCV,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,LBC_SSF,mat_cal,
     &  LVEDGE,LCYCSF,vel,prsbnd,tmpbnd,rho,ccc,yys,
     &  wifsld,LCYCOLD,
     &  prs,pp0,tmp,rva)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_model,only    : idrdp,incomp,mach0,comp,ical_vect,
     &                           ical_dens
      use module_boundary,only : kdprdc,kdolet,kdintr,kdsld,rotsld,
     &                           LBC_pair,idis,kdbuff,kdshutr,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           openout,kdpres,kdstag,stagvel,stagare,
     &                           wdbcnd,kvnslp,kvlglw,kvmodl,
     &                           kvrogf,kvEWF,kvsataW
     &                           ,kdilet
     &                           ,kdovst,masflg
      use module_species ,only : acpk,gascns,r_wm,sw,r_wm,c_pi
      use module_initial ,only : p0
      use module_io,only       : ifll,ifle
      use module_material,only : lclsd,ical_sld,dpmxv,incmp_prs,rotati,
     &                           ishaft,end,begin,nsplit
      use module_metrix,only   : tmpfac=>d2vect,tmpsld
      use module_nowtime, only : iter,time
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_metrix,only  : SHUTFL
!
! 1. Set boundary condition for pressure
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nflud
      real*8 ,intent(in)    :: deltt
      INTEGER,INTENT(IN)    :: MAT_CV(   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(  0:MXMAT)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF( MXSSFBC)
      integer,intent(in)    :: LCYCSF(  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal( 0:MXMAT)
      real*8 ,intent(in)    :: SFCENT(3,MXCVFAC)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(  MXALLCV)
      integer,intent(in)    :: ipfix (  MXMAT)
      real*8 ,intent(inout) :: vel   (  MXALLCV,3,2)
      real*8 ,intent(in)    :: rho   (  MXALLCV)
      real*8 ,intent(in)    :: tmp   (  MXALLCV)
      real*8 ,intent(in)    :: tmpbnd(  MXSSFBC)
      real*8 ,intent(in)    :: prsbnd(  MXSSFBC)
      real*8 ,intent(in)    :: ccc   (  MXALLCV)
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(in)    :: rva   (  MXCVFAC)
      real*8 ,intent(in)    :: FRSTCV(  MXSSFBC)
      real*8 ,intent(inout) :: prs   (  MXALLCV)
      real*8 ,intent(inout) :: pp0   (  MXALLCV)
      real*8 ,intent(in)    :: yys   (  MXALLCV,MXcomp)
!
! --- [local entities]
!
      real*8  :: psum,dum1,up,rmach_2,usum,dum2,dum3,coeff
      real*8  :: wi1,wi2,gamm=1.4d0,r(3),v0(3),vr(3),unit(3),dr
      real*8  :: unitP(3),dum4
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICFO,ICVA,ICVB,ICVBO,IIMAT2,IMAT2,icom
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFS,ICFE,iflag=0
      logical :: lkdwall
      real*8  :: rk=1.4d0,ttt
!
      IF(idrdp==comp) then
        do 2000 IIMAT=1,NMAT !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 2000
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(idrdp.eq.comp) then
          pp0(ICVS:ICVE)=prs(ICVS:ICVE)
          pp0(IDCS:IDCE)=prs(IDCS:IDCE)
        endif
 2000   continue
      endif
!
! ----------------------------------------------------------------
!
      iflag=iflag+1
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT=MAT_NO(IIMAT)
      IMAT2=MAT_NO(IIMAT2)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc) then
! --- periodic BC
        if(idis(nb)==0) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          prs(IDC)=prs(ICVP)
          prs(IDCP)=prs(ICV)
          pp0(IDC)=pp0(ICVP)
          pp0(IDCP)=pp0(ICV)
          enddo
        elseif(idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          prs(IDC)=0.d0!prs(ICV)
          pp0(IDC)=0.d0!pp0(ICV)
          enddo
        endif
      elseif(kd==kdbuff) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
!          dum1=0.5d0*(prs(ICVP)+prs(ICV))
!          prs(ICV)=dum1
!          prs(ICVP)=dum1
!          prs(IDC)=dum1
!          prs(IDCP)=dum1
          prs(IDC)=prs(ICVP)
          prs(IDCP)=prs(ICV)
          pp0(IDC)=pp0(ICVP)
          pp0(IDCP)=pp0(ICV)
        enddo
      elseif(kd==kdshutr) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          if(SHUTFL(IBFL)==0) then
            prs(IDC)=prs(ICVP)
            prs(IDCP)=prs(ICV)
            pp0(IDC)=pp0(ICVP)
            pp0(IDCP)=pp0(ICV)
          else
            prs(IDC)=prs(ICV)    !prs(ICVP)
            prs(IDCP)=prs(ICVP)  !prs(ICV)
            pp0(IDC)=pp0(ICV)    !pp0(ICVP)
            pp0(IDCP)=pp0(ICVP)  !pp0(ICV)
          endif
        enddo
      elseif(kd==kdovst) then 
        do IBFL=IBFS,IBFE 
        ICFL=LBC_SSF(IBFL)
        ICVP=LCYCSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
!
!        prs(IDC)=prs(ICVP)
!        pp0(IDC)=pp0(ICVP)
!
        if(rva(ICFL)>0.d0) then
!          prs(IDC)=prs(ICV)       !prs(ICVP)
!          pp0(IDC)=pp0(ICV)       !pp0(ICVP)
!          prs(IDC)=prs(ICVP)
!          pp0(IDC)=pp0(ICVP)
        else 
!          prs(IDC)=prs(ICVP)
!          pp0(IDC)=pp0(ICVP)
!          prs(IDC)=prs(ICV)       !prs(ICVP)
!          pp0(IDC)=pp0(ICV)       !pp0(ICVP)
        endif
        prs(IDC)=prs(ICV)       !prs(ICVP)
        pp0(IDC)=pp0(ICV)       !pp0(ICVP)
!        prs(IDC)=prs(ICVP)
!        pp0(IDC)=pp0(ICVP)
        enddo
      elseif(kd==kdsld) then
       IF(ical_vect) then
         do IBFL=IBFS,IBFE
         ICFP=LCYCSF(IBFL)
         ICFO=LCYCOLD(IBFL)
         ICVB=LVEDGE(1,ICFP)
         ICVBO=LVEDGE(1,ICFO)
         tmpsld(IBFL,1)=wifsld(IBFL)*prs(ICVB)
     &                 +(1.d0-wifsld(IBFL))*prs(ICVBO)
         tmpsld(IBFL,2)=pp0(ICVB)
         enddo

         do IBFL=IBFS,IBFE
         ICFL=LBC_SSF(IBFL)
         IDC=LVEDGE(2,ICFL)
         prs(IDC)=tmpsld(IBFL,1)
         pp0(IDC)=tmpsld(IBFL,2)
         enddo
       else
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        ICFO=LCYCOLD(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(1,ICFP)
        ICVBO=LVEDGE(1,ICFO)
        wi1=wifsld(IBFL)
        wi2=1.d0-wi1
        prs(IDC)=wi1*prs(ICVB)+wi2*prs(ICVBO)
        pp0(IDC)=pp0(ICVB)
        enddo
       endif
      elseif(kd==kdintr) then
!NEC$ ivdep
        do 1300 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          prs(IDC)=prs(ICV)
          prs(IDCP)=prs(ICVP)
          pp0(IDC)=pp0(ICV)
          pp0(IDCP)=pp0(ICVP)
 1300   continue
!      elseif(kd==kdtchi) then
! --- touch-inlet BC
!        do 1300 IBFL=IBFS,IBFE
!          ICFL=LBC_SSF(IBFL)
!          IDC=LVEDGE(2,ICFL)
!          ICFP=LCYCSF(IBFL)
!          ICVP=LVEDGE(1,ICFP)
!          prs(IDC)=prs(ICVP)
!          pp0(IDC)=pp0(ICVP)
! 1300   continue
      else             !zhang-cvd
        if(idrdp==comp)then
          if(kd==kdolet) then
            if(openout(nb)==1) then
              DO IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              dum1=
     &            (SFAREA(1,ICFL)*vel(IDC,1,1)
     &            +SFAREA(2,ICFL)*vel(IDC,2,1)
     &            +SFAREA(3,ICFL)*vel(IDC,3,1))
              if(dum1.gt.0.d0) then  ! outgoing
                psum=prsbnd(IBFL)!-dum1*dum1*0.5d0*rho(IDC)
              else                   ! incoming
                psum=prsbnd(IBFL)-dum1*dum1*0.5d0*rho(IDC)
              endif
              prs(IDC)=psum
              pp0(IDC)=psum
              ENDDO
            elseif(openout(nb)==2) then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              prs(IDC)=min(prsbnd(IBFL),prs(ICV))
              pp0(IDC)=prs(IDC)
              enddo
            elseif(openout(nb)==7) then !4444 
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              
              up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &          +vel(ICV,2,1)*SFAREA(2,ICFL)
     &          +vel(ICV,3,1)*SFAREA(3,ICFL)
              dum1=up*up/ccc(ICV)
!              if(prs(ICV)>prsbnd(IBFL)) then
              if(dum1<1.d0) then
                prs(IDC)=prsbnd(IBFL)  
                pp0(IDC)=pp0(ICV) 
              else
                prs(IDC)=prs(ICV) 
                pp0(IDC)=prs(ICV) 
              endif
              prs(IDC)=min(prsbnd(IBFL),prs(ICV))
              pp0(IDC)=prs(ICV) 
              enddo
            elseif(openout(nb)==9) then !4444 
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              
              up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &          +vel(ICV,2,1)*SFAREA(2,ICFL)
     &          +vel(ICV,3,1)*SFAREA(3,ICFL)
              dum1=up*up/ccc(ICV)
              if(dum1<1.d0) then
                prs(IDC)=prsbnd(IBFL)  
                pp0(IDC)=pp0(ICV) 
              else
                prs(IDC)=prs(ICV) 
                pp0(IDC)=prs(ICV) 
              endif

              prs(IDC)=prs(ICV) 
              pp0(IDC)=prs(ICV) 

              enddo

            elseif(openout(nb)==6.or.openout(nb)==10) then     !4444 
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              up=vel(IDC,1,1)*SFAREA(1,ICFL)
     &          +vel(IDC,2,1)*SFAREA(2,ICFL)
     &          +vel(IDC,3,1)*SFAREA(3,ICFL)
              ttt=tmp(IDC)
              if(sw) then
                rk=0.d0
                do ICOM=1,ncomp
                dum1=c_pi(ttt,ICOM)
                rk=rk+yys(IDC,ICOM)*(dum1/(dum1-gascns*r_wm(ICOM)))
                enddo
              else
                rk=0.d0
                do ICOM=1,ncomp
                dum1=(((
     &            acpk(5,ICOM) *ttt
     &           +acpk(4,ICOM))*ttt
     &           +acpk(3,ICOM))*ttt
     &           +acpk(2,ICOM))*ttt
     &           +acpk(1,ICOM)
                rk=rk+yys(IDC,ICOM)*(dum1/(dum1-gascns*r_wm(ICOM)))
                enddo
              endif
              dum2=rk/(1.d0-rk)
              dum3=abs(up)/sqrt(ccc(IDC))
              dum1=prsbnd(IBFL)*
     &            (dum3*(rk-1.d0)/2.d0+1.d0)**dum2
              if(dum3>1.d0) then
                prs(IDC)=max(dum1,prs(ICV))
                pp0(IDC)=dum1  !max(dum1,prs(ICV))   !  pp0(ICV)
              else
                prs(IDC)=prs(ICV)  !max(dum1,prs(ICV))
                pp0(IDC)=pp0(ICV)  !max(dum1,prs(ICV))   !  pp0(ICV)
              endif
              enddo
            elseif(openout(nb)==8) then     !4444 
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              prs(IDC)=prs(ICV)
              pp0(IDC)=prsbnd(IBFL) 
              enddo
            else
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              prs(IDC)=prsbnd(IBFL) 
              pp0(IDC)=pp0(ICV)
              enddo
            endif
          elseif(kd==kdpres) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            prs(IDC)=prsbnd(IBFL)
            pp0(IDC)=pp0(ICV)
            enddo
          elseif(kd==kdstag) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            prs(IDC)=prs(ICV)
            pp0(IDC)=pp0(ICV)
            enddo
!          elseif(kd==kdilet.and.(ical_dens==4.or.masflg(nb)==1) then
          elseif(kd==kdilet.and.ical_dens==4) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
!            prs(IDC)=prs(ICV)     !max(prs(ICV),prsbnd(IBFL))  (2)
!            pp0(IDC)=prsbnd(IBFL) !max(prs(ICV),prsbnd(IBFL))

!            prs(IDC)=prs(ICV)   !max(prs(ICV),prsbnd(IBFL))  (1)
!            pp0(IDC)=max(prs(ICV),prsbnd(IBFL)) 

            prs(IDC)=prs(ICV)   !,prsbnd(IBFL)) ! (3)
            pp0(IDC)=prsbnd(IBFL) !max(prs(ICV),prsbnd(IBFL))

!             prs(IDC)=prs(ICV) ! (4)
!             pp0(IDC)=prs(ICV)
            enddo
          elseif(kd==kdilet.and.masflg(nb)==1) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            prs(IDC)=prs(ICV)     !max(prs(ICV),prsbnd(IBFL))
            pp0(IDC)=prsbnd(IBFL) !max(prs(ICV),prsbnd(IBFL))  !prsbnd(IBFL)
            enddo
          else
!NEC$ ivdep
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            prs(IDC)=prs(ICV)
            pp0(IDC)=pp0(ICV)
            enddo
          endif
        elseif(idrdp==incomp.or.idrdp==mach0) then
          if(kd==kdolet) then
            if(openout(nb)==3.and.lclsd(IIMAT)==1) then
              call FFRABORT(1,'open_air can NOT be 3')
            endif
!
            if(openout(nb)==1) then
              DO IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              dum1=
     &            (SFAREA(1,ICFL)*vel(IDC,1,1)
     &            +SFAREA(2,ICFL)*vel(IDC,2,1)
     &            +SFAREA(3,ICFL)*vel(IDC,3,1))
              if(dum1.gt.0.d0) then  ! outgoing
                psum=0.d0 
              else                   ! incoming
                psum=-dum1*dum1*0.5d0*rho(IDC)
              endif
              prs(IDC)=psum  !prs(ICV)
              pp0(IDC)=pp0(ICV)
              ENDDO
            elseif(openout(nb)==2) then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              prs(IDC)=prs(ICV)
              pp0(IDC)=pp0(ICV)
              enddo
            elseif(openout(nb)==4.and.lclsd(IIMAT)==0) then  !sliding fan
              unit(1)=end(1,IMAT)-begin(1,IMAT)
              unit(2)=end(2,IMAT)-begin(2,IMAT)
              unit(3)=end(3,IMAT)-begin(3,IMAT)
              dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
              unit(1)=unit(1)/dum1
              unit(2)=unit(2)/dum1
              unit(3)=unit(3)/dum1
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              dum3=
     &          (SFAREA(1,ICFL)*vel(IDC,1,1)
     &          +SFAREA(2,ICFL)*vel(IDC,2,1)
     &          +SFAREA(3,ICFL)*vel(IDC,3,1))
              r(1)=SFCENT(1,ICFL)-begin(1,IMAT)
              r(2)=SFCENT(2,ICFL)-begin(2,IMAT)
              r(3)=SFCENT(3,ICFL)-begin(3,IMAT)
              call AXB_UNIT_C(unit,r,vr)
              dr=r(1)*unit(1)
     &          +r(2)*unit(2)
     &          +r(3)*unit(3)
              r(1)=r(1)-dr*unit(1)
              r(2)=r(2)-dr*unit(2)
              r(3)=r(3)-dr*unit(3)
              dr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
              v0(1)=dr*rotati(IMAT)*vr(1)
              v0(2)=dr*rotati(IMAT)*vr(2)
              v0(3)=dr*rotati(IMAT)*vr(3)
              dum1=-0.5d0*rho(ICV)*(v0(1)**2+v0(2)**2+v0(3)**2)
              dum2=-0.5d0*rho(ICV)!*dum3*dum3
     &            *(vel(ICV,1,1)**2
     &             +vel(ICV,2,1)**2
     &             +vel(ICV,3,1)**2)
              if(dum3>0.d0) then
                prs(IDC)=prs(ICV)
                pp0(IDC)=pp0(ICV)
              else
                prs(IDC)=prs(ICV)+dum2
                pp0(IDC)=prs(ICV)+dum2
              endif
              enddo
            elseif(openout(nb)==3.and.lclsd(IIMAT)==0) then
              do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICV=LVEDGE(1,ICFL)
                IDC=LVEDGE(2,ICFL)
                prs(IDC)=0.d0
                pp0(IDC)=pp0(ICV)
              enddo
            elseif(openout(nb)==5) then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              prs(IDC)=prs(ICV)
              pp0(IDC)=pp0(ICV)
              enddo
            else
!NEC$ ivdep
              IF(ical_vect) then
                do IBFL=IBFS,IBFE
                tmpfac(IBFL,1)=prs(LVEDGE(1,LBC_SSF(IBFL)))
                enddo
                do IBFL=IBFS,IBFE
                prs(LVEDGE(2,LBC_SSF(IBFL)))=tmpfac(IBFL,1)
                enddo
                do IBFL=IBFS,IBFE
                tmpfac(IBFL,1)=pp0(LVEDGE(1,LBC_SSF(IBFL)))
                enddo
                do IBFL=IBFS,IBFE
                pp0(LVEDGE(2,LBC_SSF(IBFL)))=tmpfac(IBFL,1)
                enddo
              else
                do IBFL=IBFS,IBFE
                ICFL=LBC_SSF(IBFL)
                ICV=LVEDGE(1,ICFL)
                IDC=LVEDGE(2,ICFL)
                prs(IDC)=prs(ICV)
                pp0(IDC)=pp0(ICV)
                enddo
              endif
            endif
          elseif(kd==kdpres) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            prs(IDC)=prsbnd(IBFL)-incmp_prs(IIMAT) 
            pp0(IDC)=pp0(ICV) 
            enddo
          elseif(kd==kdstag) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=0.5d0*rho(ICV)*(vel(ICV,1,1)**2+
     &                           vel(ICV,2,1)**2+
     &                           vel(ICV,3,1)**2)
            prs(IDC)=prsbnd(IBFL)-incmp_prs(IIMAT)-dum1
            pp0(IDC)=pp0(ICV)
            enddo
          else
            kdv=kdbcnd(1,nb)
            lkdwall=kdv==kvnslp.or.kdv==kvlglw.or.
     &              kdv==kvrogf.or.kdv==kvmodl.or.
     &              kdv==kvEWF.or.kdv==kvsataW
            if(lkdwall) then
!NEC$ ivdep
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              dum1=SFAREA(1,ICFL)*(vel(ICV,1,1)-vel(IDC,1,1))
     &            +SFAREA(2,ICFL)*(vel(ICV,2,1)-vel(IDC,2,1))
     &            +SFAREA(3,ICFL)*(vel(ICV,3,1)-vel(IDC,3,1))
              dum2=rho(ICV)*dum1*dum1
              dum1=dum1/(abs(dum1)+SML)
              prs(IDC)=prs(ICV)!+dum1*dum2
              pp0(IDC)=pp0(ICV)!+dum1*dum2
              enddo
            else
!NEC$ ivdep
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              prs(IDC)=prs(ICV)
              pp0(IDC)=pp0(ICV)
              enddo
            endif
          endif
        endif
      endif
 1000 continue

!
!      IF(idrdp==comp.and.ical_dens/=4) then
!        do 2000 IIMAT=1,NMAT !ICV=1,NCV
!        if(.not.mat_cal(IIMAT)) goto 2000
!        IMAT=MAT_NO(IIMAT)
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_CVEXT(IIMAT)
!        IDCS=MAT_DCIDX(IIMAT-1)+1
!        IDCE=MAT_DCIDX(IIMAT)
!        if(idrdp.eq.comp) then
!          pp0(ICVS:ICVE)=prs(ICVS:ICVE)
!          pp0(IDCS:IDCE)=prs(IDCS:IDCE)
!        endif
! 2000   continue
!      endif
!
      end subroutine bc_prs
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_rans(LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &                   rva,rva2,aksbnd,aks,kdbk)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdilet,kdolet,kdtchi,kdintr,
     &                           kdsymm,kdprdc,kkdirc,kkneum,
     &                           kdstag,kdpres,kdsld,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &     rotsld,idis,LBC_pair,kdbuff,kdshutr
     &                           ,kdovst
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_scalar,only   : iaph
      use module_model,only    : ical_VECT
      use module_metrix,only   : tmpfac=>d2vect,tmpsld
      use module_metrix,only   : SHUTFL
!
! 1. Set boundary condition for k & epsilon
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: rva       (  MXCVFAC)
      real*8 ,intent(in)    :: rva2      (  MXCVFAC2)
      real*8 ,intent(in)    :: aksbnd    (  MXSSFBCR,MXrans)
      real*8 ,intent(inout) :: aks       (  MXALLCVR,MXrans)
      integer,intent(out)   :: kdbk      (  MXCVFAC)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMD
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: IPR
      real*8  :: dum1
!
! --- 
!
      kdbk(:)=0
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdsymm) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbk(ICFL)=2    !2
        enddo
      elseif(kd.ne.kdprdc.and.kd/=kdbuff.and.kd/=kdshutr) then
        if(kdbcnd(4,nb).eq.kkneum) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          kdbk(ICFL)=2
          enddo
          if(kd.eq.kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            kdbk(ICFP)=2
            enddo
          endif
        else
          if(kd==kdsld) then    ! zhang1111  4444
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbk(ICFL)=2
            enddo
          elseif(kd==kdovst) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbk(ICFL)=999
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            kdbk(ICFL)=1   !1
            enddo
          endif
!
          if (kd.eq.kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            kdbk(ICFP)=1
            enddo
          endif
        endif
      elseif(kd==kdprdc.or.kd==kdbuff) then
        if(idis(nb)==0) then
          do 1800 IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          kdbk(ICFP)=999
 1800     continue
        elseif(idis(nb)==1) then   !zhang5
          DO IPR=1,2
          if(IPR==1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_pair(nb)
          else
            IBFS=LBC_pair(nb)+1
            IBFE=LBC_INDEX(nb)
          endif
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          kdbk(ICFL)=999  !!!!!hui
          enddo
          ENDDO
        endif
      elseif(kd==kdshutr) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          if(SHUTFL(IBFL)==0) then
            kdbk(ICFP)=999
          else
            kdbk(ICFP)=999
            kdbk(ICFL)=999
          endif
        enddo
      endif
 1000 enddo
!
!-< 2. Set value >-
!
!
!--< 2.1 set dummy cells >--
!     
      do 2000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdtchi) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(1,ICFP)
        IDC=LVEDGE(2,ICFL)
        do 103 IMD=1,nrans
        aks(IDC,IMD)=aks(ICVP,IMD)
 103    continue
        enddo
      elseif(kd==kdprdc.or.kd==kdbuff) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 102 IMD=1,nrans
        aks(IDC,IMD)=aks(ICVP,IMD)
        aks(IDCP,IMD)=aks(ICV,IMD)
 102    continue
        enddo
      elseif(kd==kdshutr) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        if(SHUTFL(IBFL)==0) then
          do IMD=1,nrans
          aks(IDC,IMD)=aks(ICVP,IMD)
          aks(IDCP,IMD)=aks(ICV,IMD)
          enddo
        else
          do IMD=1,nrans
          aks(IDC,IMD)=aks(ICV,IMD)!aks(ICVP,IMD)
          aks(IDCP,IMD)=aks(ICVP,IMD)!aks(ICV,IMD)
          enddo
        endif
        enddo
      elseif(kd==kdsld) then
        if(ical_vect) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IMD=1,nrans
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          tmpsld(IBFL,1)=aks(ICVP,IMD)
          enddo

          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          aks(IDC,IMD)=tmpsld(IBFL,1)
          enddo
          enddo
        else
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        DO IMD=1,nrans
        aks(IDC,IMD)=aks(ICVP,IMD)
        ENDDO
        enddo
        endif
      elseif(kd.eq.kdovst) then 
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICVP=LCYCSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        DO IMD=1,nrans
        aks(IDC,IMD)=aks(ICVP,IMD)
        ENDDO
        enddo
      elseif(kd.eq.kdintr) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 105 IMD=1,nrans
        aks(IDC,IMD)=aks(ICV,IMD)
        aks(IDCP,IMD)=aks(ICVP,IMD)
 105    enddo
        enddo
      else
        if(ical_vect) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IMD=1,nrans
        do IBFL=IBFS,IBFE
        tmpfac(IBFL,1)=aks(LVEDGE(1,LBC_SSF(IBFL)),IMD)
        enddo
        do IBFL=IBFS,IBFE
        aks(LVEDGE(2,LBC_SSF(IBFL)),IMD)=tmpfac(IBFL,1)
        enddo
        enddo

        else

        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        do 101 IMD=1,nrans
        aks(IDC,IMD)=aks(ICV,IMD)
  101   enddo
        enddo
        endif
      endif
!
!--< 2.2 inlet,outlet >--
!
      if(kd.eq.kdilet.OR.kd.eq.kdstag.or.kd.eq.kdpres)
     &  then 
        if(ical_VECT) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IMD=1,nrans
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
        if(rva(ICFL).lt.0.d0) then   !7777
          aks(IDC,IMD)=aksbnd(IBFL,IMD)
        endif
        enddo
        enddo
        else
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
        if(rva(ICFL).lt.0.d0) then   !7777
          do 110 IMD=1,nrans
          aks(IDC,IMD)=aksbnd(IBFL,IMD)
  110     enddo
        endif
        enddo
        endif
      elseif(kd.eq.kdolet) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IMD=1,nrans
        if(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) then
          if(IMD==iaph(1)) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=rva(ICFL)!+rva2(ICFL)
            if(dum1.lt.0.d0) then
              aks(IDC,IMD)=aksbnd(IBFL,IMD) !aks(ICV,IMD)
            endif
            enddo
          elseif(IMD==iaph(2)) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=rva2(ICFL)
            if(dum1.lt.0.d0) then
              aks(IDC,IMD)=aksbnd(IBFL,IMD) !aks(ICV,IMD)
            endif
            enddo
          endif
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          if(rva(ICFL).lt.0.d0) then
            aks(IDC,IMD)=aksbnd(IBFL,IMD)
          else
            aks(IDC,IMD)=aks(ICV,IMD)
          endif
          enddo
        endif
        enddo
      endif
!
!------------------
!--< 2.3 others >--
!------------------
!
      kdk=kdbcnd(4,nb)
      if(kdk.eq.kkdirc) then
        do IMD=1,NRANS
!        if(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) cycle
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        aks(IDC,IMD)=aksbnd(IBFL,IMD)
        enddo
        enddo
      endif
 2000 enddo
!
      return
      end subroutine bc_rans
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_alpha(LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &                   aksbnd,aks)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdilet,kdolet,kdtchi,kdintr,
     &                           kdsymm,kdprdc,kkdirc,kkneum,
     &                           kdstag,kdpres,kdsld,kdbuff,kdshutr,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &     rotsld,idis,LBC_pair,kdovst
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_scalar,only  : iaph
!
! 1. Set boundary condition for k & epsilon
!
      implicit none
!
! --- [dummy arguments]
!
      INTEGER,INTENT(IN)    :: MAT_NO(       0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF    (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      real*8 ,intent(in)    :: aksbnd(       MXSSFBCR,MXrans)
      real*8 ,intent(inout) :: aks   (       MXALLCVR,MXrans)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMD
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: IPR
!
! --- 
!
!
      do IMD=1,NRANS
      if(IMD==iaph(1).or.IMD==iaph(2)) then
        do 1000 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) goto 1000
        kd=kdbcnd(0,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
!
!-< 1. Set value >-
!
        if(kd.eq.kdtchi) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(1,ICFP)
          IDC=LVEDGE(2,ICFL)
          aks(IDC,IMD)=aks(ICVP,IMD)
          enddo
        elseif(kd==kdprdc.or.kd==kdbuff) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          aks(IDC,IMD)=aks(ICVP,IMD)
          aks(IDCP,IMD)=aks(ICV,IMD)
          enddo
        elseif(kd==kdovst) then
          call FFRABORT(1,'ERR: kdovst for e2p BC')
        elseif(kd==kdsld) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          aks(IDC,IMD)=aks(ICVP,IMD)
          enddo
        elseif(kd.eq.kdintr) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          aks(IDC,IMD)=aks(ICV,IMD)
          aks(IDCP,IMD)=aks(ICV,IMD)
          aks(ICVP,IMD)=aks(ICV,IMD)
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          aks(IDC,IMD)=min(1.d0,aks(ICV,IMD))
          enddo
        endif
!
!--< 2.2 inlet,outlet >--
!
        if(kd.eq.kdilet) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          aks(IDC,IMD)=aksbnd(IBFL,IMD)
          enddo
        endif
!
        if(kd.eq.kdolet) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          aks(IDC,IMD)=aksbnd(IBFL,IMD)
          enddo
        endif
 1000   enddo
!
      endif
      enddo
!
      return
      end subroutine bc_alpha
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_ransdif
     &  (LBC_SSF,LCYCSF,mat_cal,aksbnd,rva,rva2,dflxk)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kkneum,kdintr,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           kdilet,kdolet,kdstag,kdpres,kdsld,
     &                           kdshutr,kdbuff,kklglw
     &                           ,kdovst
      use module_VOF     ,only : ical_vof
      use module_scalar  ,only : ivof,iaph
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_scalar,only   : iaph,ike
      use module_metrix,only   : SHUTFL
!
! 1. Set boundary condition for diffusion flux of
!    dependent variables of RANS model
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)    ::  LBC_SSF   (  MXSSFBC)
      integer,intent(in)    ::  LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    ::  mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    ::  aksbnd    (  MXSSFBCR,mxrans)
      real*8 ,intent(in)    ::  rva       (  MXCVFAC)
      real*8 ,intent(in)    ::  rva2      (  MXCVFAC2)
      real*8 ,intent(inout) ::  dflxk     (  MXCVFACR,mxrans)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMD
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer,parameter :: phs_1=1,phs_2=2
!
!
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      kdk=kdbcnd(4,nb)
      kd=kdbcnd(0,nb)
      if(kd==kdshutr.or.kd==kdbuff) goto 2000
      if(kdk==kkneum) then
        do 100 IMD=1,nrans
        if(ical_vof.eq.1.and.IMD.eq.ivof) then
           cycle
        elseif(ieul2ph>0.and.
     &        (IMD==iaph(phs_1).or.IMD==iaph(phs_2))) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          dflxk(ICFL,IMD)=0.d0
          enddo
          if (kd.eq.kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            dflxk(ICFP,IMD)=0.d0
            ENDDO
          endif  
        else
          
          do 1100 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          dflxk(ICFL,IMD)=aksbnd(IBFL,IMD)
 1100     enddo
!          if() 
          if (kd==kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            dflxk(ICFP,IMD)=aksbnd(IBFL,IMD)
            ENDDO
          endif
        endif
  100   enddo
      elseif(kdbcnd(4,nb)==kklglw) then
         do IBFL=IBFS,IBFE
           ICFL=LBC_SSF(IBFL)
           dflxk(ICFL,ike(1))=0.d0
         ENDDO
      endif
!
      if(kd==kdsld) then
        do IMD=1,nrans
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        dflxk(ICFL,IMD)=0.d0
        enddo
        enddo
      endif
      if(kd==kdovst) then
        do IMD=1,nrans
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        dflxk(ICFL,IMD)=0.d0
        enddo
        enddo
      endif
!
 2000 continue
      if(kd==kdshutr) then
        do IBFL=IBFS,IBFE
        ICFP=LCYCSF(IBFL)
        ICFL=LBC_SSF(IBFL)
        if(SHUTFL(IBFL)==0) then
        else
          dflxk(ICFP,IMD)=0.d0
          dflxk(ICFL,IMD)=0.d0
        endif
        ENDDO
      endif
!
      if(kd.eq.kdolet) then
        do IMD=1,nrans
        if(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) then
          if(IMD==iaph(1)) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            if(rva(ICFL).lt.0.d0) then
              dflxk(ICFL,IMD)=0.0d0
            endif
            enddo
          elseif(IMD==iaph(2)) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            if(rva2(ICFL).lt.0.d0) then
              dflxk(ICFL,IMD)=0.d0
            endif
            enddo
          endif
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          if(rva(ICFL).lt.0.d0) then
            dflxk(ICFL,IMD)=0.d0
          endif
          enddo
        endif
        enddo
      endif
 1000 continue
!
      return
!
      end subroutine bc_ransdif
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_rva
     & (mph,time,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &  LCYCOLD,wifsld,OPPANG,
     &  CVCENT,SFCENT,wiface,locmsh,aks,
     &  SFAREA,rho,vel,xta,velbnd,prsbnd,aksbnd,prs,rva,imode)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil ,only : NPE
      use module_boundary,only : kdprdc,kdilet,kdolet,kdtchi,kdintr,
     &                           kdbuff,kdsymm,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           openout,kdstag,kdpres,kdsld,kdshutr,
     &                           rotsld,idis,LBC_pair,masflg,sumflx,
     &                           sumare,dvel
     &                           ,kdovst
!
      use module_material,only : ical_sld,rotati,ishaft,end,begin,
     &                           rot_ang
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_scalar,  only : iaph,ical_cavi
      use module_model,only    : ical_vect,ical_dens
      use module_metrix,only   : tmpfac=>d2vect
      use module_metrix,only  : SHUTFL
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)   :: time
      INTEGER,INTENT(IN)   :: imode,mph
      INTEGER,INTENT(IN)   :: MAT_NO(  0:MXMAT)
      integer,intent(in)   :: LVEDGE(2,MXCVFAC)
      integer,intent(in)   :: LBC_SSF( MXSSFBC)
      integer,intent(in)   :: LCYCSF(  MXSSFBC)
      logical,INTENT(IN)   :: mat_cal( 0:MXMAT)
      real*8 ,intent(in)   :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)   :: rho   (  MXALLCV)
      real*8 ,intent(inout):: vel   (  MXALLCV,3,2)
      real*8 ,intent(in)   :: xta   (  MXCVFAC)
      real*8 ,intent(in)   :: velbnd(  MXSSFBC,3)
      real*8 ,intent(in)   :: prsbnd ( MXSSFBC)
      real*8 ,intent(in)   :: prs    ( MXALLCV)
      real*8 ,intent(inout):: rva    ( MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3,MXALLCV)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      integer,intent(in)    :: locmsh     (MXSSFBC)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wiface     (MXCVFAC)
      REAL*8 ,INTENT(IN)    :: AKSBND     (MXSSFBCR,MXRANS)
      real*8 ,intent(inout) :: aks        (MXALLCVR,mxrans)
!
! --- [local entities]
!
      real*8  :: rvaicf,rvaicfp
      real*8  :: costh1,sinth1,costh2,sinth2,th(2),avsf
      real*8  :: ru,rv,rw,ru1,rv1,rw1,ru2,rv2,rw2,ru3,rv3,rw3,gf2,
     &           sfn1,sfn2,sfn3,sfo1,sfo2,sfo3
      integer :: IMAT,IIMAT,i,j,k,IVA,IVB,IDCA,IDCB,IDCBO,IMAT2,IIMAT2
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP,ICVB1
      real*8  :: radius(3),velref(3),unit(3,2),
     &           rbb(3,3,2),fbb(3,3,2),grdff(3,3),vell(3,3)
      integer :: IIMATS(2),IMATS(2),
     &           ISLD,ISLD2,ICFO,ICVA,ICVB,ICVBO,ICVAO
      real*8  :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,wi1,wi2,
     &           ds1,ds2,ds3,ds11,ds12,ds13,dx,dy,dz,ds21,ds22,ds23,
     &           dxo,dyo,dzo,sf1,sf2,sf3,
     &           ds31,ds32,ds33,wii1,wii2,wii1O,wii2O,
     &           avsfn,avsfo,dl,dlvect,gf1,grdf(3),vel_av(3)
      real*8  :: org_x1,org_y1,org_x2,org_y2,grx,gry,grz
      real*8  :: dl1,dl2,dl3,sf11,sf12,sf13
!
!----------------------------------------------------------
! --- 
!----------------------------------------------------------
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdprdc) then   !zhang???
        if(ieul2ph>0) call FFRABORT(1,'ERR: E2P NOT Support')
        if(idis(nb)==0) then
          if(ical_vect) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            tmpfac(IBFL,1)=0.5d0*(rva(ICFL)-rva(ICFP))
            enddo
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            rva(ICFL) = tmpfac(IBFL,1)
            rva(ICFP) =-tmpfac(IBFL,1)
            enddo
          else
            do 1100 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            rvaicf=0.5d0*(rva(ICFL)-rva(ICFP))
            rva(ICFL) = rvaicf
            rva(ICFP) =-rvaicf
 1100       enddo
          endif
        elseif(idis(nb)==1) then 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)*0.5d0*
     &        (SFAREA(1,ICFL)*(vel(IDC,1,1)+vel(ICV,1,1))
     &        +SFAREA(2,ICFL)*(vel(IDC,2,1)+vel(ICV,2,1))
     &        +SFAREA(3,ICFL)*(vel(IDC,3,1)+vel(ICV,3,1)))
     &        +xta(ICFL)
          rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
          enddo
        endif
      elseif(kd==kdbuff) then 
        do IBFL=IBFS,IBFE 
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          rvaicf=0.5d0*(rva(ICFL)-rva(ICFP)) 
          rva(ICFL) =  rvaicf
          rva(ICFP) = -rvaicf
        enddo
      elseif(kd==kdshutr) then 
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          if(SHUTFL(IBFL)==0) then
            rvaicf=0.5d0*(rva(ICFL)-rva(ICFP))
            rva(ICFL) =  rvaicf
            rva(ICFP) = -rvaicf
          else
            rva(ICFL) =0.d0
            rva(ICFP) =0.d0
          endif
        enddo
      elseif(kd==kdovst) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)  !LCYCSF(IBFL)  !attention
          ICVP=LCYCSF(IBFL)
!zhang5
!          dum1=SFAREA(4,ICFL)*0.5d0*
!     &        (SFAREA(1,ICFL)*(vel(IDC,1,1)+vel(ICV,1,1))
!     &        +SFAREA(2,ICFL)*(vel(IDC,2,1)+vel(ICV,2,1))
!     &        +SFAREA(3,ICFL)*(vel(IDC,3,1)+vel(ICV,3,1)))
!     &        +xta(ICFL)
!
          dum1=SFAREA(4,ICFL)*
     &        (SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))
     &        +xta(ICFL)
!
          rva(ICFL)=rho(ICV)*max(0.d0,dum1)
     &             +rho(ICVP)*min(0.d0,dum1)
!
        enddo
      elseif(kd==kdsld) then
          if(ieul2ph>0) call FFRABORT(1,'ERR: E2P NOT Support')
          do 2300 ISLD=1,2
          IF(ISLD==1) THEN
            ISLD2=2
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=IBFS+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
            elseif(idis(nb)>=1) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_pair(nb)
            endif
          ELSE
            ISLD2=1
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2
              IBFE=LBC_INDEX(nb)
            elseif(idis(nb)>=1) then
              IBFS=LBC_pair(nb)+1
              IBFE=LBC_INDEX(nb)
            endif
          ENDIF
!
          do 2310 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
!----------------------------------
          avsf=SFAREA(4,ICFL)
          sf11=SFAREA(1,ICFL)
          sf12=SFAREA(2,ICFL)
          sf13=SFAREA(3,ICFL)
          vell(1,1)=vel(ICVA,1,1)
          vell(1,2)=vel(ICVA,2,1)
          vell(1,3)=vel(ICVA,3,1)
          vell(2,1)=vel(ICVB,1,1)
          vell(2,2)=vel(ICVB,2,1)
          vell(2,3)=vel(ICVB,3,1)
          dum1=0.25d0*(rho(ICVA)+rho(ICVB))
          ru=(vell(1,1)+vell(2,1))*sf11*dum1
          rv=(vell(1,2)+vell(2,2))*sf12*dum1
          rw=(vell(1,3)+vell(2,3))*sf13*dum1
          gf2=avsf*(ru+rv+rw)
!-------------------------------------------------------------------------
          rvaicf=rho(ICVA)*max(0.d0,xta(ICFL))
     &          +rho(ICVB)*min(0.d0,xta(ICFL))
     &          +gf2
          rva(ICFL)=rvaicf
!
 2310     enddo
 2300     enddo
!
!          IMAT=MAT_NO(IIMAT)
!          IMAT2=MAT_NO(IIMAT2)
!          if(ishaft(IMAT)==0.and.ishaft(IMAT2)==0) then
!            IBFS=LBC_INDEX(nb-1)+1
!            IBFE=LBC_INDEX(nb)
!            do IBFL=IBFS,IBFE
!            ICFL=LBC_SSF(IBFL)
!            ICFP=LCYCSF(IBFL)
!            rvaicf=0.5d0*(rva(ICFL)-rva(ICFP))
!            rva(ICFL)= rvaicf
!            rva(ICFP)=-rvaicf
!            enddo
!          endif
!
      elseif(kd.eq.kdtchi) then
        if(ieul2ph>0) call FFRABORT(1,'ERR: E2P NOT Support')
          do 1400 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
          rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
 1400     enddo
          
      elseif(kd.eq.kdilet) then 
        if(ieul2ph>0) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum2=AKSBND(IBFL,iaph(mph))
          dum1=SFAREA(4,ICFL)*dum2
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
          rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
          enddo
        elseif(ical_dens==4) then
          
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
          rva(ICFL)=
     &             +rho(IDC)*min(0.d0,dum1)
          enddo
        else
!          dl1=0.d0
!          dum2=0.d0
          do 1200 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
          rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
!          if(nb==8) then
!            dl1=dl1+rva(ICFL)
!            dum2=dum2+SFAREA(4,ICFL)
!          endif
 1200     enddo
!          if(nb==8) then
!             print*,'BBBBBBBBb',nb,dum2,dl1,
!     &         rho(LVEDGE(1,LBC_SSF(IBFS+10))),LVEDGE(1,LBC_SSF(IBFS+10)),
!     &       LVEDGE(2,LBC_SSF(IBFS+10))
 !            print*,'BBBBBBBBb-1',nb,rho(LVEDGE(2,LBC_SSF(IBFS+10)))
!          endif
        endif
      elseif(kd.eq.kdolet) then
        if(ieul2ph==0) then      
          if(openout(nb)==6.or.openout(nb)==8) then 
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
            rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
            enddo
          elseif(openout(nb)==7.or.openout(nb)==10) then 
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
!            rva(ICFL)=rho(ICV)*max(0.d0,dum1)
            rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
            enddo
          elseif(openout(nb)==9) then 
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
!            rva(ICFL)=rho(ICV)*max(0.d0,dum1)
            rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
            enddo
          elseif(openout(nb)==1) then 
            do IBFL=IBFS,IBFE 
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL) 
            rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
            enddo
          elseif(openout(nb)==2) then 
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)
     &         *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &          +SFAREA(2,ICFL)*vel(IDC,2,1)
     &          +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
            if(dum1<0.d0) then
              rva(ICFL)=0.d0
            else
              if(prs(ICV)>prsbnd(IBFL)) then
                rva(ICFL)=rho(ICV)*dum1
              else
                rva(ICFL)=0.d0
              endif
            endif
            enddo
          elseif(openout(nb)==4.or.
     &           openout(nb)==3) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)
     &         *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &          +SFAREA(2,ICFL)*vel(IDC,2,1)
     &          +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
            rva(ICFL)=
     &            +rho(ICV)*max(0.d0,dum1)
     &            +rho(IDC)*min(0.d0,dum1)
            enddo
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
!            rva(ICFL)=rho(ICV)*max(0.d0,dum1)
!     &               +rho(IDC)*min(0.d0,dum1)
            dum2=rho(ICV)*max(0.d0,dum1)
     &               +rho(IDC)*min(0.d0,dum1)
            enddo
          endif
        else
          if(openout(nb)==5) then
            if(kdphs(mph)==kdphs_l) then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              rva(ICFL)=0.d0
              enddo
            else
              do IBFL=IBFS,IBFE 
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              dum2=aks(ICV,iaph(mph))
              dum1=SFAREA(4,ICFL)*dum2
     &           *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &            +SFAREA(2,ICFL)*vel(IDC,2,1)
     &            +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
              rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
              enddo
            endif
          else
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            dum2=aks(ICV,iaph(mph))
            dum3=aks(IDC,iaph(mph))
            dum1=SFAREA(4,ICFL)*dum2
     &         *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &          +SFAREA(2,ICFL)*vel(IDC,2,1)
     &          +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)
             rva(ICFL)=dum2*rho(ICV)*max(0.d0,dum1)
     &                +dum3*rho(IDC)*min(0.d0,dum1)
!            rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
            enddo
          endif
        endif
      elseif(kd.eq.kdstag.or.kd.eq.kdpres) then
        if(ieul2ph>0) call FFRABORT(1,'ERR: E2P NOT Support')
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,2)
     &        +SFAREA(2,ICFL)*vel(IDC,2,2)
     &        +SFAREA(3,ICFL)*vel(IDC,3,2))+xta(ICFL)
          rva(ICFL)=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
        enddo 
      elseif(kd.eq.kdintr) then 
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        rva(ICFL)=0.d0
        rva(ICFP)=0.d0
        enddo
      else
        if(kd.ne.kdolet) then
          if(kd==kdsymm) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            rva(ICFL)=0.D0
            enddo
          else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          dum1=SFAREA(4,ICFL)
     &       *(SFAREA(1,ICFL)*vel(IDC,1,1)
     &        +SFAREA(2,ICFL)*vel(IDC,2,1)
     &        +SFAREA(3,ICFL)*vel(IDC,3,1))+xta(ICFL)!????2
!          dum1=xta(ICFL)!????2  cAVITATION
          rva(ICFL)=
     &            +rho(ICV)*max(0.d0,dum1)
     &            +rho(IDC)*min(0.d0,dum1)
!          rva(ICFL)=0.d0
          enddo
          endif
        endif
      endif
!
!      if(ical_cavi==1) then
!        IF(NB==1) THEN
!        do IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        rva(ICFL)=-0.01708*1.59D-5
!        enddo
!        ENDIF
!      endif
!    
 1000 continue
!
      end subroutine bc_rva
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
      subroutine bc_setbnd(mph,iter,time,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,DISINL,SFAREA,SFCENT,locmsh,
     &  tmp,rho,vel,yys,aks,prs,XTA,ccc,
     &  prsbnd,aksbnd,
     &  velbnd,tmpbnd,yysbnd,htcbnd,mtcbnd,angbnd,SHUTFL)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
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
     &                           ,machno,machflg,ifixmach
      use module_dimnsn,only   : xredc,yredc,zredc
      use module_usersub,only  : usrno,usryes,inlusr,oulusr,walusr,
     &                           stcusr
      use module_scalar,  only : rns_scl,ike
      use module_model,   only : icaltb,ke,ke_low,RNG,CHEN,ical_vect,
     &                           ical_prt,KE2S
      use module_metrix,  only : yys_cell=>rcomp,aks_temp
      use module_time,    only : iters,itere,timee
      use module_species, only  : r_wm,gascns
!
! 1. Set boundary condition onto arrays
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: iter,mph
      real*8 ,intent(in)  :: time,deltt
      integer,intent(in)  :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)  :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)  :: LCYCSF    (   MXSSFBC)
      logical,INTENT(IN)  :: mat_cal   (   0:MXMAT)
      real*8 ,intent(out) :: prsbnd    (   MXSSFBC)
      real*8 ,intent(out) :: velbnd    (   MXSSFBC,3)
      real*8 ,intent(out) :: tmpbnd    (   MXSSFBC)
      real*8 ,intent(out) :: yysbnd    (   MXSSFBC,MXcomp)
      real*8 ,intent(out) :: aksbnd    (   MXSSFBCR,MXrans)
      real*8 ,intent(out) :: htcbnd    (   MXSSFBC)
      real*8 ,intent(out) :: mtcbnd    (   MXSSFBC)
      real*8 ,intent(in)  :: DISINL    (   MXSSFBC)
      real*8 ,intent(in)  :: SFAREA    ( 4,MXCVFAC)
      real*8 ,intent(in)  :: SFCENT    ( 3,MXCVFAC)
      real*8 ,intent(in)  :: vel   (  MXALLCV,3)
      real*8 ,intent(in)  :: rho   (  MXALLCV)
      real*8 ,intent(in)  :: tmp   (  MXALLCV)
      real*8 ,intent(in)  :: yys   (  MXALLCV,MXCOMP)
      real*8 ,intent(in)  :: aks   (  MXALLCVR,mxrans)
      real*8 ,intent(in)  :: prs   (  MXALLCV)
      integer,intent(in)  :: locmsh(  MXSSFBC)
      real*8 ,intent(out) :: XTA   (  MXCVFAC)
      real*8 ,intent(out) :: angbnd(  MXSSFBC_VOF)
      integer,intent(inout)  :: SHUTFL (MXSSFBC_SHT)
      REAL*8 ,INTENT(IN)   :: CCC  (  MXALLCV)
!
! --- [local entities] 
!
      logical :: userBC,uinl,uout,uwal,ustc,mflg
      logical :: lkdwall
      real*8  :: vmx,sgm,r1,r2,velinl(3),t,p,ydis,pmin
      real*8  :: dx,dy,dz,dx1,dy1,dz1,dis,xxx(3) 
      real*8  :: area_inl(4),radius(3),shaftn(3),dum1,dum2,dum3
      real*8  :: t_cell,p_cell,v_cell(3),rho_cell,dum_V,mv_surf
      integer,save  :: isize=0
      real*8,save   :: tur_int=0.05d0 
      integer :: i,j,k,l,m,n,no,iv,ierr1=0,ICV,IDC
      integer :: IMAT,IIMAT,ICOM,IMD
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ICF 
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICVL,IDCL,ICVP,IDCP 
      integer :: mxcvfacx,mxssfbc1,idum1,idum2 
!
      real*8,allocatable  :: keps(:),yys1(:)
!
      ALLOCATE(keps(mxrans),yys1(ncomp),stat=ierr1)
!
      call chkncomp(ncomp,'(bc_setbnd)')
      call chknrans(nrans,'(bc_setbnd)')
!
      if(isize==0.and.rns_scl.and.
     & (icaltb==KE2S.or.icaltb==ke.or.
     &  icaltb==ke_low.or.icaltb==RNG.or.icaltb==CHEN)) 
     &  then
        isize=1
        sizein(:)=0.d0
        do 2000 nb=1,nbcnd
        dum1=0.d0
        IIMAT=MAT_BCIDX(nb,1)
        kd=kdbcnd(0,nb)
        if(kd/=kdilet.and.kd/=kdstag.and.kd/=kdpres) cycle
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(IBFS>IBFE) cycle
        ICFL=LBC_SSF(IBFS)
        do 2200 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        dum1=dum1+SFAREA(4,ICFL)
 2200   enddo
        sizein(nb)=2.d0*dsqrt(dum1/3.1415d0)
 2000   enddo
!
      endif
!
      do 1000 nb=1,nbcnd   !IBF=1,NSSFBC
      IIMAT=MAT_BCIDX(nb,1)
      kd=kdbcnd(0,nb)
      if(.not.mat_cal(IIMAT)) cycle
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
!
      IF(ICAL_VECT) then
        do IBFL=IBFS,IBFE
        prsbnd(IBFL  )=wdbcnd(1,nb)
        velbnd(IBFL,1)=wdbcnd(2,nb)
        velbnd(IBFL,2)=wdbcnd(3,nb)
        velbnd(IBFL,3)=wdbcnd(4,nb)
        tmpbnd(  IBFL)=wdbcnd(5,nb)
        enddo
!
        do ICOM=1,ncomp
        do IBFL=IBFS,IBFE
        yysbnd(IBFL,ICOM)=wdbcnd(iscomp+ICOM,nb)
        enddo
        enddo
!
        if(rns_scl) then
          do IMD=1,nrans
          do IBFL=IBFS,IBFE
          aksbnd(IBFL,IMD)=wdbcnd(israns+IMD,nb)
          enddo
          enddo
!
          dum1=wdbcnd(israns+ike(1),nb)+wdbcnd(israns+ike(2),nb)
!
          if((icaltb==KE2S.or.icaltb==ke.or.icaltb==ke_low.or.
     &     icaltb==RNG.or.icaltb==CHEN).and.
     &     abs(dum1)<SML) then
            if(kd==kdilet) then
              do IBFL=IBFS,IBFE
               dum1=tur_int*
     &        (wdbcnd(2,nb)*wdbcnd(2,nb)
     &        +wdbcnd(3,nb)*wdbcnd(3,nb)
     &        +wdbcnd(4,nb)*wdbcnd(4,nb))
              aksbnd(IBFL,ike(1))=dum1
              aksbnd(IBFL,ike(2))=
     &                  0.09**(0.75D0)*dum1**(1.5d0)/(SML+sizein(nb))
              enddo
            elseif(kd==kdstag.or.kd==kdpres)then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFS)
              IDC=LVEDGE(2,ICFL)
               dum1=tur_int*
     &        (vel(IDC,1)**2+vel(IDC,2)**2+vel(IDC,3)**2)
              aksbnd(IBFL,ike(1))=dum1
              aksbnd(IBFL,ike(2))=
     &                  0.09**(0.75D0)*dum1**(1.5d0)/(SML+sizein(nb))
              enddo
            endif
          endif
        endif
!
        do IBFL=IBFS,IBFE
          htcbnd(IBFL)=adbcnd(1,nb)
          mtcbnd(IBFL)=adbcnd(2,nb)
          if(MXSSFBC_VOF==MXSSFBC) angbnd(IBFL)=vofang(nb)
        enddo
!
        if(adbcnd(3,nb).gt.0.d0) then
          if( adbcnd(4,nb).gt.0.d0 ) then
            do IBFL=IBFS,IBFE
            vmx=sqrt(velbnd(IBFL,1)**2
     &            +velbnd(IBFL,2)**2
     &            +velbnd(IBFL,3)**2)
            sgm=adbcnd(3,nb)*vmx
            vmx=adbcnd(4,nb)*vmx
            do iv=1,3
            call utl_boxmlr(r1,r2)
            velbnd(IBFL,iv)=velbnd(IBFL,iv)+max(-vmx,min(vmx,sgm*r1))
            enddo
            enddo
          endif
        endif
      else
!---------------------
! --- Scalar computer
!---------------------
        do 1200 IBFL=IBFS,IBFE
        prsbnd(IBFL  )=wdbcnd(1,nb)
        velbnd(IBFL,1)=wdbcnd(2,nb)
        velbnd(IBFL,2)=wdbcnd(3,nb)
        velbnd(IBFL,3)=wdbcnd(4,nb)
        tmpbnd(  IBFL)=wdbcnd(5,nb)

        do 100 ICOM=1,ncomp
        yysbnd(IBFL,ICOM)=wdbcnd(iscomp+ICOM,nb)
  100   enddo
!
        if(rns_scl) then
          do 101 IMD=1,nrans
          aksbnd(IBFL,IMD)=wdbcnd(israns+IMD,nb)
  101     continue
          dum1=wdbcnd(israns+ike(1),nb)+wdbcnd(israns+ike(2),nb)
          if((icaltb==KE2S.or.icaltb==ke.or.icaltb==ke_low.or.
     &     icaltb==RNG.or.icaltb==CHEN).and.
     &     abs(dum1)<SML) then
            if(kd==kdilet) then
!              do IBFL=IBFS,IBFE
               dum1=tur_int*
     &        (wdbcnd(2,nb)*wdbcnd(2,nb)
     &        +wdbcnd(3,nb)*wdbcnd(3,nb)
     &        +wdbcnd(4,nb)*wdbcnd(4,nb))
              aksbnd(IBFL,ike(1))=dum1
              aksbnd(IBFL,ike(2))=
     &                  0.09**(0.75D0)*dum1**(1.5d0)/(SML+sizein(nb))
!              enddo
            elseif(kd==kdstag.or.kd==kdpres) then
!              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFS)
              IDC=LVEDGE(2,ICFL)
               dum1=tur_int*
     &        (vel(IDC,1)**2+vel(IDC,2)**2+vel(IDC,3)**2)
              aksbnd(IBFL,ike(1))=dum1
              aksbnd(IBFL,ike(2))=
     &                  0.09**(0.75D0)*dum1**(1.5d0)/(SML+sizein(nb))
!              enddo
            endif
          endif
        endif
!
!        do IBFL=IBFS,IBFE
        htcbnd(IBFL)=adbcnd(1,nb)
        mtcbnd(IBFL)=adbcnd(2,nb)
        if(MXSSFBC_VOF==MXSSFBC) angbnd(IBFL)=vofang(nb)
!
        if(adbcnd(3,nb).gt.0.d0) then
          if(adbcnd(4,nb).gt.0.d0 ) then
            vmx=sqrt(velbnd(IBFL,1)**2
     &              +velbnd(IBFL,2)**2
     &              +velbnd(IBFL,3)**2)
            sgm=adbcnd(3,nb)*vmx
            vmx=adbcnd(4,nb)*vmx
            do 110 iv=1,3
            call utl_boxmlr(r1,r2)
            velbnd(IBFL,iv)=velbnd(IBFL,iv)+max(-vmx,min(vmx,sgm*r1))
  110       enddo
          endif
        endif
!        enddo
 1200   enddo
      endif
!
      if(kd.eq.kdintr) then
        do 1300 IBFL=IBFS,IBFE
        prsbnd(IBFL  )=wdbcnd(1,nb)
        velbnd(IBFL,1)=wdbcnd(2,nb)
        velbnd(IBFL,2)=wdbcnd(3,nb)
        velbnd(IBFL,3)=wdbcnd(4,nb)
        tmpbnd(IBFL  )=wdbcnd(5,nb)
!!
        do 120 ICOM=1,ncomp
        yysbnd(IBFL,ICOM)=wdbcnd(iscomp+ICOM,nb)
  120   continue
!!
        if(rns_scl) then
          do 121 IMD=1,nrans
          aksbnd(IBFL,IMD)=wdbcnd(israns+IMD,nb)
  121     enddo
        endif
!
        htcbnd(  IBFL)=adbcnd(1,nb)
        mtcbnd(  IBFL)=adbcnd(2,nb)
!
        if(adbcnd(3,nb).gt.0.d0) then
          if(adbcnd(4,nb).gt.0.d0 ) then
            vmx=sqrt(velbnd(IBFL,1)**2
     &              +velbnd(IBFL,2)**2
     &              +velbnd(IBFL,3)**2)
            sgm=adbcnd(3,nb)*vmx
            vmx=adbcnd(4,nb)*vmx
            do 112 iv=1,3
            call utl_boxmlr(r1,r2)
            velbnd(IBFL,iv)=velbnd(IBFL,iv)+max(-vmx,min(vmx,sgm*r1))
 112        enddo
          endif
        endif        
 1300   enddo
      endif
!
 1000 enddo
!
! --- rot wall
!
      DO 3000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      lkdwall=kd==kxnone.or.
!     &        kd==kdbuff.or.
     &        kd==kdintr.or.
     &        kd==kdfire.or.
     &        kd==kdcvd.or.
     &        kd==kdolet
      if(lkdwall.and.rotwall(nb)==1) then 
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        shaftn(1)=endw(1,nb)-beginw(1,nb)
        shaftn(2)=endw(2,nb)-beginw(2,nb)
        shaftn(3)=endw(3,nb)-beginw(3,nb)
        dum1=dsqrt(shaftn(1)**2+shaftn(2)**2+shaftn(3)**2)
        shaftn(1)=shaftn(1)/dum1
        shaftn(2)=shaftn(2)/dum1
        shaftn(3)=shaftn(3)/dum1
        if(ical_vect) then
          do IBFL=IBFS,IBFE    !vector
          ICFL=LBC_SSF(IBFL)
          dum1=(SFCENT(1,ICFL)-beginw(1,nb))*shaftn(1)
     &        +(SFCENT(2,ICFL)-beginw(2,nb))*shaftn(2)
     &        +(SFCENT(3,ICFL)-beginw(3,nb))*shaftn(3)
          radius(1)=(SFCENT(1,ICFL)-beginw(1,nb))-dum1*shaftn(1)
          radius(2)=(SFCENT(2,ICFL)-beginw(2,nb))-dum1*shaftn(2)
          radius(3)=(SFCENT(3,ICFL)-beginw(3,nb))-dum1*shaftn(3)
          velinl(3)=(shaftn(1)*radius(2)-shaftn(2)*radius(1))
          velinl(2)=(shaftn(3)*radius(1)-shaftn(1)*radius(3))
          velinl(1)=(shaftn(2)*radius(3)-shaftn(3)*radius(2))
          velbnd(IBFL,1)=velbnd(IBFL,1)+velinl(1)*omega(nb)
          velbnd(IBFL,2)=velbnd(IBFL,2)+velinl(2)*omega(nb)
          velbnd(IBFL,3)=velbnd(IBFL,3)+velinl(3)*omega(nb)
          enddo 
        else
          do 3100 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          radius(1)=SFCENT(1,ICFL)-beginw(1,nb)
          radius(2)=SFCENT(2,ICFL)-beginw(2,nb)
          radius(3)=SFCENT(3,ICFL)-beginw(3,nb)
          dum1=radius(1)*shaftn(1)
     &        +radius(2)*shaftn(2)
     &        +radius(3)*shaftn(3)
          radius(1)=radius(1)-dum1*shaftn(1)
          radius(2)=radius(2)-dum1*shaftn(2)
          radius(3)=radius(3)-dum1*shaftn(3)
          call cal_AXBEQC(shaftn,radius,velinl,dum1)
          velbnd(IBFL,1)=velbnd(IBFL,1)+velinl(1)*omega(nb)
          velbnd(IBFL,2)=velbnd(IBFL,2)+velinl(2)*omega(nb)
          velbnd(IBFL,3)=velbnd(IBFL,3)+velinl(3)*omega(nb)
 3100     enddo
        endif
      endif
 3000 ENDDO
!
!
!
!
      if(imasflg) then
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        mflg=(kd==kdilet.or.kd==kdstag.or.
!     &  (kd==kdolet.and.(openout(nb)==7)).or.
     &  (kd==kdolet.and.(openout(nb)==8)))
     &  .and.masflg(nb)>=1
        if(mflg) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(openout(nb)==8) then
            dum1=0.d0
            do IBFL=IBFS,IBFE
            if(locmsh(IBFL)==1) then
              ICFL=LBC_SSF(IBFL)
              IDC=LVEDGE(2,ICFL)
              dum2=0.d0
              do ICOM=1,ncomp 
              dum2=dum2+yysbnd(IBFL,ICOM)*r_wm(ICOM)
              enddo
              dum3=prsbnd(IBFL)/(gascns*dum2*tmpbnd(IBFL)) 
              dum1=dum1+SFAREA(4,ICFL)*dum3  !dum3=rho(IDC)
            endif
            enddo
          else
            dum1=0.d0
            if(masflg(nb)==1) then
              do IBFL=IBFS,IBFE
              if(locmsh(IBFL)==1) then 
                ICFL=LBC_SSF(IBFL)
                ICV=LVEDGE(1,ICFL)
                IDC=LVEDGE(2,ICFL)
                dum1=dum1+SFAREA(4,ICFL)*rho(IDC) 
              endif
              enddo
            elseif(masflg(nb)==2) then
              do IBFL=IBFS,IBFE
              if(locmsh(IBFL)==1) then 
                ICFL=LBC_SSF(IBFL)
                ICV=LVEDGE(1,ICFL)
                IDC=LVEDGE(2,ICFL)
                dum1=dum1+SFAREA(4,ICFL)
              endif
              enddo
            endif
          endif
          sumare(nb)=dum1+SML
!
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          dum1=masflx(nb)/sumare(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          velbnd(IBFL,1)=dum1*SFAREA(1,ICFL)
          velbnd(IBFL,2)=dum1*SFAREA(2,ICFL)
          velbnd(IBFL,3)=dum1*SFAREA(3,ICFL)
          enddo
        endif
        enddo
      endif
!
      if(ifixmach) then
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        mflg=(kd==kdilet.and.machflg(nb)==1).or.
     &       (kd==kdolet.and.machflg(nb)==1.and.openout(nb)==10)
        
        if(mflg) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          dum1=-machno(nb)*sqrt(CCC(IDC))
          velbnd(IBFL,1)=dum1*SFAREA(1,ICFL)
          velbnd(IBFL,2)=dum1*SFAREA(2,ICFL)
          velbnd(IBFL,3)=dum1*SFAREA(3,ICFL)

          enddo
        endif
        enddo
      endif
!

!
! --- define velocity profile in INLET
!
      userBC=.false.
      uinl=inlusr.eq.usryes
      uout=oulusr.eq.usryes
      uwal=walusr.eq.usryes
      ustc=stcusr.eq.usryes
      userBC=uinl.or.uout.or.uwal.or.ustc
!
!------------------------------------------------------------
!------------ for Compressor Labyrinth leakage --------------
!------------------------------------------------------------
!      if(.false.) then
!        call user_BC_leak(
!     &     mph,iter,time,mat_cal,
!     &     LBC_SSF,LVEDGE,DISINL,prsbnd,velbnd,tmpbnd,
!     &     aksbnd,aks,yysbnd,yys,SFAREA,SFCENT,
!     &     tmp,prs,vel,rho)
!      endif
!------------------------------------------------------------
!------------------------------------------------------------


      if(userBC) then
!        call FFRABORT(1,'MSG:NOT Support E2P')
! --- User defined Vel. Profile
       if(ical_vect) then
         call user_BC_VECTOR
     &  (mph,iter,time,deltt,iters,itere,timee,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,DISINL,SFAREA,SFCENT,locmsh,
     &  tmp,rho,vel,yys,aks,prs,XTA,
     &  prsbnd,aksbnd,
     &  velbnd,tmpbnd,yysbnd,htcbnd,mtcbnd,angbnd,SHUTFL
     &  )
       else
        idum1=iters
        idum2=itere
        do 4000 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        kd=kdbcnd(0,nb)
        lkdwall=kd==kxnone.or.
!     &          kd==kdbuff.or.
     &          kd==kdintr.or.
     &          kd==kdfire.or.
     &          kd==kdcvd
        if((kd.eq.kdilet.and.uinl).or.
     &     (kd.eq.kdolet.and.uout).or.
     &     (lkdwall.and.uwal).or.
     &     (kd.eq.kdpres.and.ustc)
     &) then
          pmin=1.d25
          do 2420 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVL=LVEDGE(1,ICFL)
          IDCL=LVEDGE(2,ICFL)
!
          no=nobcnd(nb)
          ydis=DISINL(IBFL)
          p=prsbnd(IBFL  )   !wdbcnd(1,nb)
          velinl(1)=velbnd(IBFL,1)   !wdbcnd(2,nb)
          velinl(2)=velbnd(IBFL,2)
          velinl(3)=velbnd(IBFL,3)
          t=tmpbnd(  IBFL)   !wdbcnd(5,nb)
          if(rns_scl) then
            do 2400 IMD=1,nrans
            keps(IMD)=aksbnd(IBFL,IMD)
            aks_temp(IMD)=aks(ICVL,IMD)
 2400       enddo
          endif
          do 2450 ICOM=1,NCOMP
          yys1(ICOM)=yysbnd(IBFL,ICOM)
          yys_cell(ICOM)=yys(ICVL,ICOM)
 2450     enddo
          area_inl(1)=SFAREA(1,ICFL)
          area_inl(2)=SFAREA(2,ICFL)
          area_inl(3)=SFAREA(3,ICFL)
          area_inl(4)=SFAREA(4,ICFL)
          xxx(1)=SFCENT(1,ICFL)
          xxx(2)=SFCENT(2,ICFL)
          xxx(3)=SFCENT(3,ICFL)
          t_cell=tmp(ICVL)
          p_cell=prs(ICVL)
          v_cell(1:3)=vel(ICVL,1:3)
          rho_cell=rho(ICVL)
          p=prsbnd(IBFL)
          mv_surf=0.d0
          t=tmpbnd(IBFL)
          call user_BC(mph,iter,time,ICFL,no,ical_prt,
     &      mxrans,MXcomp,area_inl,xxx,idum1,idum2,
     &      ICVL,t_cell,p_cell,v_cell,aks_temp,rho_cell,yys_cell,
     &      ydis,velinl,p,t,keps,yys1,mv_surf)
          prsbnd(  IBFL)=p
          pmin=min(p,pmin)
          velbnd(IBFL,1)=velinl(1)
          velbnd(IBFL,2)=velinl(2)
          velbnd(IBFL,3)=velinl(3)
          tmpbnd(  IBFL)=t 
          XTA(ICFL)=mv_surf 
          if(rns_scl) then
            do 2410 IMD=1,nrans
            aksbnd(IBFL,IMD)=keps(IMD)
 2410       enddo
          endif
 2420     enddo
          wdbcnd(1,nb)=pmin
        endif
 4000   enddo
       endif
      endif
!
      do nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      kd=kdbcnd(0,nb)
      if(.not.mat_cal(IIMAT)) cycle
      if(kd==kdshutr) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        no=nobcnd(nb)
!        dum1=shut_strt(nb)
        dum2=deltt
        dum3=time
        idum1=iter
        idum2=MXSSFBC_SHT
        MXCVFACX=MXCVFAC
        MXSSFBC1=MXSSFBC
        call USER_SHUTTER(idum1,dum3,dum2,no,idum2,MXSSFBC1,MXCVFACX,
     &       IBFS,IBFE,SHUTFL,SFCENT,LBC_SSF)
      endif
      enddo
!
      if( xredc ) velbnd(:,1)=0.d0
      if( yredc ) velbnd(:,2)=0.d0
      if( zredc ) velbnd(:,3)=0.d0
!
      DEALLOCATE(keps,yys1)
!
      end subroutine bc_setbnd
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_setbnd2(mph,iter,time,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,DISINL,SFAREA,SFCENT,!locmsh,
     &  prsbnd,velbnd2,tmpbnd2,yysbnd2,htcbnd2,mtcbnd2)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : chkncomp,chknrans,nobcnd,
     &                           wdbcnd2,adbcnd2,iscomp,israns,distrb,
     &                           LBC_INDEX,MAT_BCIDX,nbcnd,kdbcnd,
     &                           kdintr
      use module_dimnsn,only   : xredc,yredc,zredc
      use module_usersub,only  : usrno,usryes,inlusr
      use module_scalar,  only : rns_scl
!
! 1. Set boundary condition onto arrays
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mph,iter
      real*8 ,intent(in)  :: time
      integer,intent(in)  :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)  :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)  :: LCYCSF    (   MXSSFBC)
      logical,INTENT(IN)  :: mat_cal   (   0:MXMAT)
      real*8 ,intent(in)  :: DISINL    (   MXSSFBC)
      real*8 ,intent(in)  :: SFAREA    ( 4,MXCVFAC)
      real*8 ,intent(in)  :: SFCENT    ( 3,MXCVFAC)
!
      real*8 ,intent(out) :: velbnd2   (   MXSSFBC2,3)
      real*8 ,intent(out) :: tmpbnd2   (   MXSSFBC2)
      real*8 ,intent(out) :: yysbnd2(      MXSSFBC2,MXcomp)
      real*8 ,intent(out) :: htcbnd2   (   MXSSFBC2)
      real*8 ,intent(out) :: mtcbnd2   (   MXSSFBC2)
      real*8 ,intent(in) :: prsbnd    (   MXSSFBC)
!      integer,intent(in)    :: locmsh (  MXSSFBC)
!
! --- [local entities]
!
      logical :: uinl
      real*8  :: vmx,sgm,r1,r2,velinl(3),t,p,ydis,keps(nrans),yys(ncomp)
      integer :: i,j,k,l,m,n,no,iv
      integer :: IMAT,IIMAT,ICOM,IMD
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
!
!
      call chkncomp(ncomp,'(bc_setbnd2)')
      call chknrans(nrans,'(bc_setbnd2)')
!
      do 1000 nb=1,nbcnd   !IBF=1,NSSFBC
      IIMAT=MAT_BCIDX(nb,1)
      kd=kdbcnd(0,nb)
      if(.not.mat_cal(IIMAT)) goto 1000
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      do 1200 IBFL=IBFS,IBFE
      velbnd2(IBFL,1)=wdbcnd2(2,nb)
      velbnd2(IBFL,2)=wdbcnd2(3,nb)
      velbnd2(IBFL,3)=wdbcnd2(4,nb)
      tmpbnd2(  IBFL)=wdbcnd2(5,nb)
!
      do 100 ICOM=1,ncomp
      yysbnd2(IBFL,ICOM)=wdbcnd2(iscomp+ICOM,nb)
  100 continue
!
      htcbnd2(  IBFL)=adbcnd2(1,nb)
      mtcbnd2(  IBFL)=adbcnd2(2,nb)
!
      if(adbcnd2(3,nb).gt.0.d0) then
        if( adbcnd2(4,nb).gt.0.d0 ) then
          vmx=sqrt(velbnd2(IBFL,1)**2
     &            +velbnd2(IBFL,2)**2
     &            +velbnd2(IBFL,3)**2)
          sgm=adbcnd2(3,nb)*vmx
          vmx=adbcnd2(4,nb)*vmx
          do 110 iv=1,3
          call utl_boxmlr(r1,r2)
          velbnd2(IBFL,iv)=velbnd2(IBFL,iv)+max(-vmx,min(vmx,sgm*r1))
  110     continue
        endif
      endif
 1200 continue
!
      if(kd.eq.kdintr) then
        do 1300 IBFL=IBFS,IBFE
        ICFP=LCYCSF(IBFL)
        velbnd2(ICFP,1)=wdbcnd2(2,nb)
        velbnd2(ICFP,2)=wdbcnd2(3,nb)
        velbnd2(ICFP,3)=wdbcnd2(4,nb)
        tmpbnd2(  ICFP)=wdbcnd2(5,nb)
!
        do 120 ICOM=1,ncomp
        yysbnd2(ICFP,ICOM)=wdbcnd2(iscomp+ICOM,nb)
  120   enddo
!
        htcbnd2(  ICFP)=adbcnd2(1,nb)
        mtcbnd2(  ICFP)=adbcnd2(2,nb)
!
        if(adbcnd2(3,nb).gt.0.d0) then
          if( adbcnd2(4,nb).gt.0.d0 ) then
            vmx=sqrt(velbnd2(ICFP,1)**2
     &              +velbnd2(ICFP,2)**2
     &              +velbnd2(ICFP,3)**2)
            sgm=adbcnd2(3,nb)*vmx
            vmx=adbcnd2(4,nb)*vmx
            do 112 iv=1,3
            call utl_boxmlr(r1,r2)
            velbnd2(ICFP,iv)=velbnd2(ICFP,iv)+max(-vmx,min(vmx,sgm*r1))
 112     continue
          endif
        endif
 1300   continue
      endif
 1000 continue
!
      if( xredc ) velbnd2(:,1)=0.d0
      if( yredc ) velbnd2(:,2)=0.d0
      if( zredc ) velbnd2(:,3)=0.d0
!
      end subroutine bc_setbnd2
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_boxmlr(x1,x2)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. normal random number generator (Box-Mullar method)
!
! [dummy arguments]
      real*8,intent(out) :: x1,x2
! [local entities]
      real*8,parameter :: pi2=2.d0*3.14159265358979d0
      real*8 :: r,t
!
!
!
      call utl_random(x1)
      call utl_random(x2)
      r=sqrt(-2.d0*log(x1))
      t=pi2*x2
      x1=r*sin(t)
      x2=r*cos(t)
      end subroutine utl_boxmlr
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
      subroutine bc_tys(mph,
     &   LVEDGE,LCYCSF,LBC_SSF,mat_cal,MAT_DCIDX,MAT_NO,
     &   LCYCOLD,wifsld,
     &   rva,tmpbnd,yysbnd,vel,rho,ccc,tmp,yys,hhh,prs)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil,only  : my_rank 
      use module_boundary,only : kdprdc,kdilet,kdolet,kdtchi,kdintr,
     &                           ktdirc,kttrns,kydirc,kytrns,kdfire,
     &                           kdsld,kdcvd,kdbuff,kyneum,kdshutr,
     &                           ktEWF,kyEWF,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                        openout,kdpres,kdstag,stagvel,stagare,
     &                           rotsld,idis,LBC_pair
     &                           ,kdovst
!
      use module_species,only  : sw
      use module_species ,only : acpk,gascns,wm,r_wm
      use module_model,   only : ical_t,ical_vect
      use module_metrix,only   : tmpfac=>d2vect,tmpsld
      use module_metrix,only  : SHUTFL
      use module_metrix,only   : tmpsld
!
! 1. Set boundary condition for temperature & mass fraction
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mph
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      real*8 ,intent(in)    :: rva       (   MXCVFAC)
      real*8 ,intent(in)    :: tmpbnd    (   MXSSFBC)
      real*8 ,intent(in)    :: yysbnd    (   MXSSFBC,MXcomp)
      real*8 ,intent(in)    :: vel       (   MXALLCV,3,2)
      real*8 ,intent(in)    :: rho       (   MXALLCV)
      real*8 ,intent(inout) :: tmp       (   MXALLCV)
      real*8 ,intent(inout) :: yys       (   MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: hhh       (   MXALLCV)
      real*8 ,intent(in)    :: ccc       (   MXALLCV)
      integer,intent(in)    :: LCYCOLD   (   MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      real*8 ,intent(in)    :: prs   (       MXALLCV)
! 
! --- [local entities]
!
      integer :: IMAT,IIMAT,ICOM,ierr1=0,i
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP,
     &           ICFO,ICVBO
      REAL*8  :: g=9.8d0,t_stat,dt_v,rk,uvel,wi1,wi2,
     &           dum1,dum2,dum3
      REAL*8  :: grx,gry,grz
      logical :: tflg
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      kdt=kdbcnd(2,nb)
      kdy=kdbcnd(3,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdprdc.or.kd==kdbuff)THEN 
!--< 1.1 periodic BC >--
        do 1200 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        tmp(IDC)=tmp(ICVP)
        tmp(IDCP)=tmp(ICV)
        hhh(IDC)=hhh(ICVP)
        hhh(IDCP)=hhh(ICV)
        do 102 ICOM=1,ncomp
        yys(IDC,ICOM)=yys(ICVP,ICOM)
        yys(IDCP,ICOM)=yys(ICV,ICOM)
  102   continue
 1200   continue
      elseif(kd==kdshutr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        if(SHUTFL(IBFL)==0) then
          tmp(IDC)=tmp(ICVP)
          tmp(IDCP)=tmp(ICV)
          hhh(IDC)=hhh(ICVP)
          hhh(IDCP)=hhh(ICV)
          do ICOM=1,ncomp
          yys(IDC,ICOM)=yys(ICVP,ICOM)
          yys(IDCP,ICOM)=yys(ICV,ICOM)
          enddo
        else
          tmp(IDC)=tmp(ICV)!tmp(ICVP)
          tmp(IDCP)=tmp(ICVP)!tmp(ICV)
          hhh(IDC)=hhh(ICV)!hhh(ICVP)
          hhh(IDCP)=hhh(ICVP)!hhh(ICV)
          do ICOM=1,ncomp
          yys(IDC,ICOM)=yys(ICV,ICOM)!yys(ICVP,ICOM)
          yys(IDCP,ICOM)=yys(ICVP,ICOM)!yys(ICV,ICOM)
          enddo
        endif
        enddo
      elseif(kd==kdovst) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVP=LCYCSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          tmp(IDC)=tmp(ICVP)
          hhh(IDC)=hhh(ICVP)
          do ICOM=1,ncomp
          yys(IDC,ICOM)=yys(ICVP,ICOM)
          ENDDO
        ENDDO
      elseif(kd==kdsld) then
        if(ical_vect) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          tmpsld(IBFL,1)=tmp(ICVP)
          tmpsld(IBFL,2)=hhh(ICVP)
          enddo
          
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          tmp(IDC)=tmpsld(IBFL,1)
          hhh(IDC)=tmpsld(IBFL,2)
          enddo

          
          do ICOM=1,ncomp
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          tmpsld(IBFL,1)=yys(ICVP,ICOM)
          ENDDO
          
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          yys(IDC,ICOM)=tmpsld(IBFL,1)
          enddo
          ENDDO
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)
          ICVBO=LVEDGE(1,ICFO)
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
          tmp(IDC)=tmp(ICVP)
          hhh(IDC)=hhh(ICVP)
          do ICOM=1,ncomp
          yys(IDC,ICOM)=yys(ICVP,ICOM)
          ENDDO
          ENDDO
        endif
      elseif(kd.eq.kdtchi) then
!--< 1.2 touch-inlet BC >--
        do 1400 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        hhh(IDC)=hhh(ICVP)
        tmp(IDC)=tmp(ICVP)
        do 104 ICOM=1,ncomp
        yys(IDC,ICOM)=yys(ICVP,ICOM)
 104    enddo
 1400   enddo
      elseif(kd.eq.kdintr) then
        do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
!            ICFP=LCYCSF(IBFL)
!            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            ICVP=LVEDGE(1,ICFP)
!            IDCP=LVEDGE(2,ICFP) 
            tmp(IDC)=tmp(ICVP)
            tmp(IDCP)=tmp(ICV)
            do ICOM=1,ncomp
            yys(IDC,ICOM)=yys(ICVP,ICOM)
            enddo
        enddo
      else
!--< 1.3 Other BC >--
        if(kd==kdilet.or.
     &    (kd==kdolet.and.(
     &     openout(nb)==6.or.
     &     openout(nb)==7.or.
     &     openout(nb)==10.or.
     &     openout(nb)==1))) then
          do 1300 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          if(rva(ICFL).lt.0.d0.or.abs(rva(ICFL))<SML) then
            tmp(IDC)=tmpbnd(IBFL)
            do 110 ICOM=1,ncomp
            yys(IDC,ICOM)=yysbnd(IBFL,ICOM)
  110       enddo
          else
            tmp(IDC)=tmp(ICV)
            hhh(IDC)=hhh(ICV)
            do 103 ICOM=1,ncomp
            yys(IDC,ICOM)=yys(ICV,ICOM)
  103       enddo
          endif
 1300     enddo
        elseif(kd==kdolet.and.openout(nb)==8) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          tmp(IDC)=tmpbnd(IBFL)
          do ICOM=1,ncomp
          yys(IDC,ICOM)=yysbnd(IBFL,ICOM)
          enddo
          enddo
        elseif(kd==kdolet.and.openout(nb)==2) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
! --- include open-air BC (openout(nb).eq.1)
          tmp(IDC)=tmp(ICV)
          hhh(IDC)=hhh(ICV)
          do ICOM=1,ncomp
          yys(IDC,ICOM)=yys(ICV,ICOM)
          enddo
          enddo
        else
          IF(ical_vect) then
            do IBFL=IBFS,IBFE
              tmpfac(IBFL,1)=tmp(LVEDGE(1,LBC_SSF(IBFL)))
!              tmp(LVEDGE(2,LBC_SSF(IBFL)))=tmp(LVEDGE(1,LBC_SSF(IBFL)))
            enddo
            
            do IBFL=IBFS,IBFE
            tmp(LVEDGE(2,LBC_SSF(IBFL)))=tmpfac(IBFL,1)
            enddo
            
            do IBFL=IBFS,IBFE
            tmpfac(IBFL,1)=hhh(LVEDGE(1,LBC_SSF(IBFL)))
!              hhh(LVEDGE(2,LBC_SSF(IBFL)))=hhh(LVEDGE(1,LBC_SSF(IBFL)))
            enddo
            
            do IBFL=IBFS,IBFE
            hhh(LVEDGE(2,LBC_SSF(IBFL)))=tmpfac(IBFL,1)
            enddo
            
            DO ICOM=1,ncomp
              do IBFL=IBFS,IBFE
              tmpfac(IBFL,1)=yys(LVEDGE(1,LBC_SSF(IBFL)),ICOM)
              enddo
              do IBFL=IBFS,IBFE
              yys(LVEDGE(2,LBC_SSF(IBFL)),ICOM)=tmpfac(IBFL,1)
              enddo
            enddo
          else
            do 1100 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            tmp(IDC)=tmp(ICV)
            hhh(IDC)=hhh(ICV)
            do 101 ICOM=1,ncomp
            yys(IDC,ICOM)=yys(ICV,ICOM)
  101       enddo
 1100       enddo
          endif
        endif
      endif
!
! --- 
!
      tflg=kd==kdintr.or.
     &     kd==kdbuff.or.
     &     kd==kdshutr.or.
     &     kd==kdsld.or.
     &     kd==kdovst
!
      if((kdt==ktdirc.or.kdt==kttrns.or.kdt==ktEWF)
     &        .and.(.not.tflg)) then
        if(kd.eq.kdintr.and.kdt.eq.ktdirc) 
     &  call FFRABORT(1,'interface NOT Dir or Tra')
        do 1500 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
        tmp(IDC)=tmpbnd(IBFL)
 1500   enddo
      endif
!
! --- 
!
      if(kdy.eq.kydirc.or.kdy.eq.kytrns.or.kdy==kyEWF) then
        IF(ical_vect) then
          do ICOM=1,ncomp
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          yys(IDC,ICOM)=yysbnd(IBFL,ICOM)
          ENDDO
          enddo
        else
          if(kd/=kdintr) then
          do 1600 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          yys(IDC,:)=yysbnd(IBFL,:)
 1600     enddo
          endif
        endif

!        if(kd.eq.kdintr) then
!           do IBFL=IBFS,IBFE
!           ICFP=LCYCSF(IBFL)
!           IDC=LVEDGE(2,ICFP)
!           yys(IDC,:)=yysbnd(IBFL,:)
!           enddo
!        endif
      endif
!
      if(kd.eq.kdpres.or.kd.eq.kdstag) then
        if(ical_vect) then
          do ICOM=1,ncomp
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          dum1=sign(1.d0,rva(ICFL))
          dum2=0.5d0*(1.d0+dum1)
          dum3=0.5d0*(1.d0-dum1)
          yys(IDC,ICOM)=dum2*yys(ICV,ICOM)
     &                 +dum3*yysbnd(IBFL,ICOM)
          enddo
          enddo
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          dum1=sign(1.d0,rva(ICFL))
          dum2=0.5d0*(1.d0+dum1)
          dum3=0.5d0*(1.d0-dum1)
          tmp(IDC)=dum2*tmp(ICV)+dum3*tmpbnd(IBFL)
          hhh(IDC)=dum2*hhh(ICV)
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          if(rva(ICFL).gt.0.d0) then
            do ICOM=1,ncomp
            yys(IDC,ICOM)=yys(ICV,ICOM)
            enddo
            tmp(IDC)=tmp(ICV)
            hhh(IDC)=hhh(ICV)
          else
            do ICOM=1,ncomp
            yys(IDC,ICOM)=yysbnd(IBFL,ICOM)
            enddo
            tmp(IDC)=tmpbnd(IBFL)
          endif
          enddo
        endif
      endif
 1000 continue
!
      if(ical_t) then
        if(sw) then
          call cal_t2h_dc(mph,MAT_DCIDX,MAT_NO,mat_cal,
     &       tmp,yys,hhh,prs,rho)
        else
          call cal_t2hcp_dc(mph,MAT_DCIDX,MAT_NO,mat_cal,
     &        tmp,yys,hhh)
        endif     
      endif
!
      return
!
      end subroutine bc_tys
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_vel
     & (mph,deltt,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,SFCENT,
     &  LCYCOLD,wifsld,OPPANG,locmsh,
     &  xta,rva,rho,velbnd,prsbnd,prs,ccc,aksbnd,
     &  vel,tmp,yys,imode)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil ,only : NPE,my_rank
      use module_model,only    : idrdp,incomp,mach0,comp,ical_vect
      use module_boundary,only : kdilet,kdolet,kvnslp,kvlglw,kvfslp,
     &                           kdintr,kdtchi,kvmodl,kvEWF,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           openout,kdpres,kdstag,stagvel,rotwall,
     &                           wdbcnd,kdprdc,kvrogf,set_rotmtrx,
     &                           masflg,sumflx,masflx,sumare,dvel,
     &                           imasflg,kdbuff,kdsld
     &                           ,kdovst,kvsataW,
     &                           machno,machflg
      use module_material,only : dpmxv,incmp_prs,rotati,ishaft,end,
     &                           begin,nsplit
      use module_species ,only : gascns,r_wm,acpk,sw,c_pi
      use module_time,    only : steady
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_scalar,  only : iaph
      use module_metrix,only   : tmpsld
!
! 1. Set boundary condition for velocity
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MAT_NO (  0:MXMAT),imode,mph
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      real*8 ,intent(in)    :: deltt
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT (3,MXCVFAC)
      real*8 ,intent(in)    :: xta    (  MXCVFAC)
      real*8 ,intent(inout) :: rva    (  MXCVFAC)
      real*8 ,intent(in)    :: rho    (  MXALLCV)
      real*8 ,intent(in)    :: velbnd (  MXSSFBC,3)
      real*8 ,intent(in)    :: prsbnd (  MXSSFBC)
      real*8 ,intent(in)    :: prs    (  MXALLCV)
      real*8 ,intent(inout) :: vel    (  MXALLCV,3,2)
      real*8 ,intent(in)    :: ccc    (  MXALLCV)
      REAL*8 ,INTENT(IN)    :: AKSBND (  MXSSFBCR,MXRANS)
      integer,intent(in)    :: LCYCOLD(  MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (  MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG (  MXSSFBC_SLD)
      integer,intent(in)    :: locmsh (  MXSSFBC)
      real*8 ,intent(in)    :: tmp    (  MXALLCV)
      real*8 ,intent(in)    :: yys    (  MXALLCV,MXcomp)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ICOM,nb2
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP,I
      real*8  :: prd,vrd,ux,uy,uz,up,dum1,dum2,dum3,dum4,usum,rk=1.4d0,
     &           ttt,r(3),v0(3),vr(3),unit(3),dr,rot(3,3),
     &           dum_mas,dum_are
      integer :: ICVA,ICVB,ICVLA,ICVLB
!
      call dc_symprv
     &(1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & LCYCOLD,wifsld,OPPANG,
     & SFAREA,SFCENT,vel(:,:,1),0,1)

! --- 
!
      if(imasflg) then
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd.eq.kdtchi) then
          dvel(nb)=0.d0
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          dum_mas=0.d0
          dum_are=0.d0
          call set_rotmtrx(nbcnd,kd,nb,rot)
          if(masflg(nb)>=1) then
            dum3=0.d0
            do IBFL=IBFS,IBFE
            if(locmsh(IBFL)==1) then
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            ICVP=LVEDGE(1,ICFP)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            ux= rot(1,1)*vel(ICVP,1,1)
     &         +rot(1,2)*vel(ICVP,2,1)
     &         +rot(1,3)*vel(ICVP,3,1)
            uy= rot(2,1)*vel(ICVP,1,1)
     &         +rot(2,2)*vel(ICVP,2,1)
     &         +rot(2,3)*vel(ICVP,3,1)
            uz= rot(3,1)*vel(ICVP,1,1)
     &         +rot(3,2)*vel(ICVP,2,1)
     &         +rot(3,3)*vel(ICVP,3,1)
            dum1=SFAREA(4,ICFL)
     &         *(SFAREA(1,ICFL)*ux
     &          +SFAREA(2,ICFL)*uy
     &          +SFAREA(3,ICFL)*uz)+xta(ICFL)
            if(masflg(nb)==1) then
              dum2=rho(ICV)*max(0.d0,dum1)+rho(IDC)*min(0.d0,dum1)
              dum_mas=dum_mas+dum2
!              dum_are=dum_are+SFAREA(4,ICFL)*rho(ICV)
              dum_are=dum_are+SFAREA(4,ICFL)*rho(IDC)
              dum3=dum3+SFAREA(4,ICFL)
            else
              dum2=max(0.d0,dum1)+min(0.d0,dum1)
              dum_mas=dum_mas+dum2
              dum_are=dum_are+SFAREA(4,ICFL)
              dum3=dum3+SFAREA(4,ICFL)
            endif
            endif
            enddo
            if(NPE.gt.1) then
              CALL hpcrsum(dum_mas)
              CALL hpcrsum(dum_are)
            endif
            sumflx(nb)=dum_mas
            sumare(nb)=dum_are
            dvel(nb)=(sumflx(nb)-masflx(nb))/(sumare(nb)+SML)
          endif
!
          dum_mas=0.d0
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          dum1=rot(1,1)*vel(ICVP,1,1)
     &        +rot(1,2)*vel(ICVP,2,1)
     &        +rot(1,3)*vel(ICVP,3,1)
     &        -SFAREA(1,ICFL)*dvel(nb)
          dum2=rot(2,1)*vel(ICVP,1,1)
     &        +rot(2,2)*vel(ICVP,2,1)
     &        +rot(2,3)*vel(ICVP,3,1)
     &        -SFAREA(2,ICFL)*dvel(nb)
          dum3=rot(3,1)*vel(ICVP,1,1)
     &        +rot(3,2)*vel(ICVP,2,1)
     &        +rot(3,3)*vel(ICVP,3,1)
     &        -SFAREA(3,ICFL)*dvel(nb)
          vel(IDC,1,1) =dum1
          vel(IDC,2,1) =dum2
          vel(IDC,3,1) =dum3
          enddo
        endif
        enddo
      endif
!----------------
! --- BC
!----------------
      do 1000 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) goto 1000
        kd=kdbcnd(0,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd.eq.kdilet) then
          IMAT=MAT_NO(IIMAT)
          if(ishaft(IMAT)==0) then
            if(ieul2ph>0) then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              dum2=AKSBND(IBFL,iaph(mph))
              vel(IDC,1,1)=velbnd(IBFL,1)
              vel(IDC,2,1)=velbnd(IBFL,2)
              vel(IDC,3,1)=velbnd(IBFL,3)
              enddo
            else
              do 1100 IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              vel(IDC,1,1)=velbnd(IBFL,1)
              vel(IDC,2,1)=velbnd(IBFL,2)
              vel(IDC,3,1)=velbnd(IBFL,3)
 1100         enddo
            endif
          elseif(ishaft(IMAT)==1) then
            unit(1)=end(1,IMAT)-begin(1,IMAT)
            unit(2)=end(2,IMAT)-begin(2,IMAT)
            unit(3)=end(3,IMAT)-begin(3,IMAT)
            dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
            unit(1)=unit(1)/dum1
            unit(2)=unit(2)/dum1
            unit(3)=unit(3)/dum1
            
            do 1150 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            r(1)=SFCENT(1,ICFL)-begin(1,IMAT)
            r(2)=SFCENT(2,ICFL)-begin(2,IMAT)
            r(3)=SFCENT(3,ICFL)-begin(3,IMAT)
            call AXB_UNIT_C(unit,r,vr)
            dr=r(1)*unit(1)
     &        +r(2)*unit(2)
     &        +r(3)*unit(3)
            r(1)=r(1)-dr*unit(1)
            r(2)=r(2)-dr*unit(2)
            r(3)=r(3)-dr*unit(3)
            dr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
            v0(1)=dr*rotati(IMAT)*vr(1)
            v0(2)=dr*rotati(IMAT)*vr(2)
            v0(3)=dr*rotati(IMAT)*vr(3)
            vel(IDC,1,1)=velbnd(IBFL,1)-v0(1)
            vel(IDC,2,1)=velbnd(IBFL,2)-v0(2)
            vel(IDC,3,1)=velbnd(IBFL,3)-v0(3)
 1150       enddo
          endif
        endif
!
      if(kd==kdolet) then
        if(openout(nb)==1) then
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          vel(IDC,1,1)=up*SFAREA(1,ICFL)
          vel(IDC,2,1)=up*SFAREA(2,ICFL)
          vel(IDC,3,1)=up*SFAREA(3,ICFL)
 1250     continue
        elseif(openout(nb)==6) then 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          up=velbnd(IBFL,1)*SFAREA(1,ICFL)
     &      +velbnd(IBFL,2)*SFAREA(2,ICFL)
     &      +velbnd(IBFL,3)*SFAREA(3,ICFL)
!          vel(IDC,1,1)=velbnd(IBFL,1) 
!          vel(IDC,2,1)=velbnd(IBFL,2) 
!          vel(IDC,3,1)=velbnd(IBFL,3) 
          vel(IDC,1,1)=vel(ICV,1,1)
          vel(IDC,2,1)=vel(ICV,2,1)
          vel(IDC,3,1)=vel(ICV,3,1)
          enddo
        elseif(openout(nb)==8.or.openout(nb)==10) then 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          up=velbnd(IBFL,1)*SFAREA(1,ICFL)
     &      +velbnd(IBFL,2)*SFAREA(2,ICFL)
     &      +velbnd(IBFL,3)*SFAREA(3,ICFL)
          vel(IDC,1,1)=velbnd(IBFL,1) 
          vel(IDC,2,1)=velbnd(IBFL,2) 
          vel(IDC,3,1)=velbnd(IBFL,3) 
          enddo
        elseif(openout(nb)==7)then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
!
          vrd=SFAREA(1,ICFL)*velbnd(IBFL,1)
     &       +SFAREA(2,ICFL)*velbnd(IBFL,2)
     &       +SFAREA(3,ICFL)*velbnd(IBFL,3)
!
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          if(rva(ICFL).lt.0.d0) then
            prd=rva(ICFL)/(rho(IDC)*SFAREA(4,ICFL))-vrd
            vel(IDC,1,1)=velbnd(IBFL,1)+prd*SFAREA(1,ICFL)
            vel(IDC,2,1)=velbnd(IBFL,2)+prd*SFAREA(2,ICFL)
            vel(IDC,3,1)=velbnd(IBFL,3)+prd*SFAREA(3,ICFL)
          endif
          enddo
        elseif(openout(nb)==0) then 
          do 1200 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
!
          vrd=SFAREA(1,ICFL)*velbnd(IBFL,1)
     &       +SFAREA(2,ICFL)*velbnd(IBFL,2)
     &       +SFAREA(3,ICFL)*velbnd(IBFL,3)
!
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          if(rva(ICFL).lt.0.d0) then
            prd=rva(ICFL)/(rho(IDC)*SFAREA(4,ICFL))-vrd
! --- Dirichlet condition
            vel(IDC,1,1)=velbnd(IBFL,1)+prd*SFAREA(1,ICFL)
            vel(IDC,2,1)=velbnd(IBFL,2)+prd*SFAREA(2,ICFL)
            vel(IDC,3,1)=velbnd(IBFL,3)+prd*SFAREA(3,ICFL)
! --- mass flux=0 
!            vel(IDC,1,1)=vel(ICV,1,1)-up*SFAREA(1,ICFL)
!            vel(IDC,2,1)=vel(ICV,2,1)-up*SFAREA(2,ICFL)
!            vel(IDC,3,1)=vel(ICV,3,1)-up*SFAREA(3,ICFL)
          endif
 1200     enddo
        elseif(openout(nb)==2) then  !zhang-cvd
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          if(up<0.d0) then
            vel(IDC,1,1)=0.d0
            vel(IDC,2,1)=0.d0
            vel(IDC,3,1)=0.d0
          else
            if(prs(ICV)>prsbnd(IBFL)) then
              vel(IDC,1,1)=up*SFAREA(1,ICFL)
              vel(IDC,2,1)=up*SFAREA(2,ICFL)
              vel(IDC,3,1)=up*SFAREA(3,ICFL)
            else
              vel(IDC,1,1)=0.d0
              vel(IDC,2,1)=0.d0
              vel(IDC,3,1)=0.d0
            endif
          endif
          enddo
        elseif(openout(nb)==4) then  !sliding fan
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          if(up<0.d0) then
            vel(IDC,1,1)=vel(ICV,1,1)-up*SFAREA(1,ICFL)
            vel(IDC,2,1)=vel(ICV,2,1)-up*SFAREA(2,ICFL)
            vel(IDC,3,1)=vel(ICV,3,1)-up*SFAREA(3,ICFL)
          else
            vel(IDC,1,1)=vel(ICV,1,1)
            vel(IDC,2,1)=vel(ICV,2,1)
            vel(IDC,3,1)=vel(ICV,3,1)
          endif
          enddo
        elseif(openout(nb)==3) then !vector
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          vel(IDC,1,1)=up*SFAREA(1,ICFL)
          vel(IDC,2,1)=up*SFAREA(2,ICFL)
          vel(IDC,3,1)=up*SFAREA(3,ICFL)
          enddo
        elseif(openout(nb)==9) then !vector
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          vel(IDC,1,1)=up*SFAREA(1,ICFL)
          vel(IDC,2,1)=up*SFAREA(2,ICFL)
          vel(IDC,3,1)=up*SFAREA(3,ICFL)
          enddo
        elseif(openout(nb)==5) then
          if(ieul2ph>0) then
            if(kdphs(mph)==kdphs_l) then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &          +vel(ICV,2,1)*SFAREA(2,ICFL)
     &          +vel(ICV,3,1)*SFAREA(3,ICFL)
              prd=up+xta(ICFL)/SFAREA(4,ICFL)
              vel(IDC,1,1)=vel(ICV,1,1)-prd*SFAREA(1,ICFL)
              vel(IDC,2,1)=vel(ICV,2,1)-prd*SFAREA(2,ICFL)
              vel(IDC,3,1)=vel(ICV,3,1)-prd*SFAREA(3,ICFL)
              enddo
            else
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICV=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              vel(IDC,1,1)=vel(ICV,1,1)
              vel(IDC,2,1)=vel(ICV,2,1)
              vel(IDC,3,1)=vel(ICV,3,1)
              enddo
            endif
          else
            call FFRABORT(1,'ERR: open_air=5 is only for E2P')
          endif
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
!
          vrd=SFAREA(1,ICFL)*velbnd(IBFL,1)
     &       +SFAREA(2,ICFL)*velbnd(IBFL,2)
     &       +SFAREA(3,ICFL)*velbnd(IBFL,3)
!
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
          if(rva(ICFL).lt.0.d0) then
            prd=rva(ICFL)/(rho(IDC)*SFAREA(4,ICFL))-vrd
            vel(IDC,1,1)=velbnd(IBFL,1)+prd*SFAREA(1,ICFL)
            vel(IDC,2,1)=velbnd(IBFL,2)+prd*SFAREA(2,ICFL)
            vel(IDC,3,1)=velbnd(IBFL,3)+prd*SFAREA(3,ICFL)
          endif
          enddo
        endif
      elseif(idrdp.eq.comp.and.kd.eq.kdpres) then
        if(steady) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          vel(IDC,1,1)=vel(ICV,1,1)
          vel(IDC,2,1)=vel(ICV,2,1)
          vel(IDC,3,1)=vel(ICV,3,1)
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          vel(IDC,1,1)=vel(ICV,1,1)
          vel(IDC,2,1)=vel(ICV,2,1)
          vel(IDC,3,1)=vel(ICV,3,1)
          enddo
        endif
      elseif(idrdp.eq.comp.and.kd.eq.kdstag) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
!        ttt=tmp(ICV)
        ttt=tmp(IDC)
        if(sw) then
          rk=0.d0
          do ICOM=1,ncomp
          dum1=c_pi(ttt,ICOM)
          rk=rk+yys(IDC,ICOM)*(dum1/(dum1-gascns*r_wm(ICOM)))
          enddo
        else
          rk=0.d0
          do ICOM=1,ncomp
          dum1=(((
     &      acpk(5,ICOM) *ttt
     &     +acpk(4,ICOM))*ttt
     &     +acpk(3,ICOM))*ttt
     &     +acpk(2,ICOM))*ttt
     &     +acpk(1,ICOM)
          rk=rk+yys(IDC,ICOM)*(dum1/(dum1-gascns*r_wm(ICOM)))
          enddo
        endif
        dum1=(rk-1.d0)/rk
        dum3=max(prsbnd(IBFL)/prs(ICV),1.d0)
!!!!!!!!!!!!!!
!        dum3=prsbnd(IBFL)/prs(ICV)
!        dum2=(dum3**dum1-1.d0)*2.d0/(rk-1.d0)*ttt*gascns
!!!!!!!!!!!!!!
        dum2=(2.d0/(rk-1.d0)*((dum3**dum1-1.d0))*ccc(ICV))
        dum3=(dum3-1.d0)/abs(1.d0-dum3+SML)
        
        usum=dsqrt(abs(dum2))  !*ccc(ICV)
        vel(IDC,1,1)=-SFAREA(1,ICFL)*usum!*dum3
        vel(IDC,2,1)=-SFAREA(2,ICFL)*usum!*dum3
        vel(IDC,3,1)=-SFAREA(3,ICFL)*usum!*dum3
        enddo
!        do IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        ICV=LVEDGE(1,ICFL)
!        IDC=LVEDGE(2,ICFL)
!        vel(IDC,1,1)=vel(ICV,1,1)
!        vel(IDC,2,1)=vel(ICV,2,1)
!        vel(IDC,3,1)=vel(ICV,3,1)
!        enddo
      elseif((idrdp.eq.incomp.or.idrdp.eq.mach0)
     &      .and.kd.eq.kdpres) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        vel(IDC,1,1)=vel(ICV,1,1)
        vel(IDC,2,1)=vel(ICV,2,1)
        vel(IDC,3,1)=vel(ICV,3,1)
        enddo
      elseif((idrdp.eq.incomp.or.idrdp.eq.mach0)
     &        .and.kd.eq.kdstag) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        vel(IDC,1,1)=vel(ICV,1,1)
        vel(IDC,2,1)=vel(ICV,2,1)
        vel(IDC,3,1)=vel(ICV,3,1)
        enddo
      elseif(kd.eq.kdtchi.and.masflg(nb)==0) then
        call set_rotmtrx(nbcnd,kd,nb,rot)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do 1450 IBFL=IBFS,IBFE
        if(locmsh(IBFL)==0) cycle
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        vel(IDC,1,1)= rot(1,1)*vel(ICVP,1,1)
     &               +rot(1,2)*vel(ICVP,2,1)
     &               +rot(1,3)*vel(ICVP,3,1)
        vel(IDC,2,1)= rot(2,1)*vel(ICVP,1,1)
     &               +rot(2,2)*vel(ICVP,2,1)
     &               +rot(2,3)*vel(ICVP,3,1)
        vel(IDC,3,1)= rot(3,1)*vel(ICVP,1,1)
     &               +rot(3,2)*vel(ICVP,2,1)
     &               +rot(3,3)*vel(ICVP,3,1)
 1450   continue
      elseif(kd==kdovst) then
        cycle
      elseif(kd==kdsld) then  !sldz
        cycle
        IF(ical_vect) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          tmpsld(IBFL,1)=vel(ICVP,1,1)
          tmpsld(IBFL,2)=vel(ICVP,2,1)
          tmpsld(IBFL,3)=vel(ICVP,3,1)
          enddo
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          vel(IDC,1,1)=tmpsld(IBFL,1)
          vel(IDC,2,1)=tmpsld(IBFL,2)
          vel(IDC,3,1)=tmpsld(IBFL,3)
          enddo
        else
          do IBFL=IBFS,IBFE    !zhangsld
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          vel(IDC,1,1)=vel(ICVP,1,1)
          vel(IDC,2,1)=vel(ICVP,2,1)
          vel(IDC,3,1)=vel(ICVP,3,1)
          enddo
        endif
      endif
!
      kdv=kdbcnd(1,nb)
      if(kdv==kvnslp.or.kdv==kvlglw.or.kdv==kvrogf.or.
     &   kdv==kvmodl.or.kdv==kvEWF.or.kdv==kvsataW) then
        do 1300 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        vrd=SFAREA(1,ICFL)*velbnd(IBFL,1)
     &     +SFAREA(2,ICFL)*velbnd(IBFL,2)
     &     +SFAREA(3,ICFL)*velbnd(IBFL,3)
        prd=vrd+xta(ICFL)/SFAREA(4,ICFL)  !????2 xta
        vel(IDC,1,1)=velbnd(IBFL,1)-prd*SFAREA(1,ICFL)
        vel(IDC,2,1)=velbnd(IBFL,2)-prd*SFAREA(2,ICFL)
        vel(IDC,3,1)=velbnd(IBFL,3)-prd*SFAREA(3,ICFL)
 1300   continue
        if (kd.eq.kdintr) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          IDCP=LVEDGE(2,ICFP)
          IDC=LVEDGE(2,ICFL)
          vrd=SFAREA(1,ICFP)*velbnd(IBFL,1)
     &       +SFAREA(2,ICFP)*velbnd(IBFL,2)
     &       +SFAREA(3,ICFP)*velbnd(IBFL,3)
          prd=vrd+xta(ICFP)/SFAREA(4,ICFP)
          vel(IDCP,1,1)=velbnd(IBFL,1)-prd*SFAREA(1,ICFP)
          vel(IDCP,2,1)=velbnd(IBFL,2)-prd*SFAREA(2,ICFP)
          vel(IDCP,3,1)=velbnd(IBFL,3)-prd*SFAREA(3,ICFP)
          vel(IDC,1,1)=velbnd(IBFL,1)+prd*SFAREA(1,ICFL)
          vel(IDC,2,1)=velbnd(IBFL,2)+prd*SFAREA(2,ICFL)
          vel(IDC,3,1)=velbnd(IBFL,3)+prd*SFAREA(3,ICFL)
          enddo
        endif
      elseif(kdv==kvfslp) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        vrd=SFAREA(1,ICFL)*vel(ICV,1,1)
     &     +SFAREA(2,ICFL)*vel(ICV,2,1)
     &     +SFAREA(3,ICFL)*vel(ICV,3,1)
        vel(IDC,1,1)=vel(ICV,1,1)-vrd*SFAREA(1,ICFL)
        vel(IDC,2,1)=vel(ICV,2,1)-vrd*SFAREA(2,ICFL)
        vel(IDC,3,1)=vel(ICV,3,1)-vrd*SFAREA(3,ICFL)
        enddo
      endif
!
      ! zhang 070830 buffle for Tohyo Seikan
      if (kd==kdbuff) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do I=1,3
        dum1=0.5d0*(vel(ICVP,I,1)+vel(ICV,I,1))
        vel(ICVP,I,1)=dum1
        vel(ICV,I,1)=dum1
        vel(IDC,I,1)=dum1
        vel(IDCP,I,1)=dum1
        enddo
        enddo
      endif
      ! End
!
      if(kd.eq.kdolet) then
        IMAT=MAT_NO(IIMAT)
        if(ishaft(IMAT)==1) then
          unit(1)=end(1,IMAT)-begin(1,IMAT)
          unit(2)=end(2,IMAT)-begin(2,IMAT)
          unit(3)=end(3,IMAT)-begin(3,IMAT)
          dum1=dsqrt(unit(1)**2+unit(2)**2+unit(3)**2)
          unit(1)=unit(1)/dum1
          unit(2)=unit(2)/dum1
          unit(3)=unit(3)/dum1
          do 1660 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          r(1)=SFCENT(1,ICFL)-begin(1,IMAT)
          r(2)=SFCENT(2,ICFL)-begin(2,IMAT)
          r(3)=SFCENT(3,ICFL)-begin(3,IMAT)
          call AXB_UNIT_C(unit,r,vr)
          dr=r(1)*unit(1)
     &      +r(2)*unit(2)
     &      +r(3)*unit(3)
          r(1)=r(1)-dr*unit(1)
          r(2)=r(2)-dr*unit(2)
          r(3)=r(3)-dr*unit(3)
          dr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
          v0(1)=dr*rotati(IMAT)*vr(1)
          v0(2)=dr*rotati(IMAT)*vr(2)
          v0(3)=dr*rotati(IMAT)*vr(3)
          up=vel(ICV,1,1)*SFAREA(1,ICFL)
     &      +vel(ICV,2,1)*SFAREA(2,ICFL)
     &      +vel(ICV,3,1)*SFAREA(3,ICFL)
!          vel(IDC,1,1)=SFAREA(1,ICFL)*up-v0(1)
!          vel(IDC,2,1)=SFAREA(2,ICFL)*up-v0(2)
!          vel(IDC,3,1)=SFAREA(3,ICFL)*up-v0(3)
!Unemura
          vel(IDC,1,1)=vel(ICV,1,1)
          vel(IDC,2,1)=vel(ICV,2,1)
          vel(IDC,3,1)=vel(ICV,3,1)
 1660     enddo
        endif
      endif
 1000 continue
!
      end subroutine bc_vel
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_velcofd(LVEDGE,LBC_SSF,mat_cal,cofd)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : set_rotmtrx,kdprdc,kdsld,
     &                           rotsld,idis,LBC_pair,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX
     &                           ,kdovst
!
! 1. Set diagonal element of rotational matrix
!    for periodic boundary
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)  :: LVEDGE (2,MXCVFAC)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      logical,INTENT(IN)  :: mat_cal(  0:MXMAT)
      real*8 ,intent(out) :: cofd   (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: aa(3,3)
      integer :: IMAT,IIMAT,ierr1=0
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
      cofd=1.d0
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      if(kd.eq.kdprdc) then
        call set_rotmtrx(nbcnd,kd,nb,aa)
        if(idis(nb)==0) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          cofd(IDC)=max(0.d0,min(aa(1,1),aa(2,2),aa(3,3)))
          enddo
        elseif(idis(nb)==1) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          cofd(IDC)=max(0.d0,min(aa(1,1),aa(2,2),aa(3,3)))
          enddo
        endif
      elseif(kd==kdsld) then  
        !????????????????   zhang????
        call set_rotmtrx(nbcnd,kd,nb,aa)
        if(idis(nb)==0.or.idis(nb)>=1) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          cofd(IDC)=max(0.d0,min(aa(1,1),aa(2,2),aa(3,3)))
          enddo
        endif
      endif
 1000 continue
!
      return
!
      end subroutine bc_velcofd
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_veldif(
     &  LVEDGE,LBC_SSF,LCYCSF,LFUTAU,mat_cal,MAT_NO,
     &  SFAREA,CVCENT,FRSTCV,vel,rho,rmu,rmut,utau,dsclv,
     &  dflx,rmue,IMODE)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdbcnd,kvfslp,kvlglw,kvnslp,kdintr,
     &                           kdolet,kdilet,kdtchi,kdsymm,
     &                           kvrogf,kvEWF,kvsataW,
     &                           LBC_INDEX,kdsld,rotsld,idis,kdcvd,
     &                           kdbuff,LBC_pair,kvmodl,kdshutr,kdbuff,
     &                           nbcnd,MAT_BCIDX,LBC_INDEX,rotwall
     &                           ,kdovst,kdstag
!
      use module_nowtime, only : iter,time
      use module_turbparm,only : akappa,awfnc
      use module_material,only : ical_sld,rot_ang
      use module_model,only    : icaltb,ke,ke_low,RNG,CHEN,SDES,SLES,
     &                           dles,ical_tetra,KE2S
      use module_metrix,only   : SHUTFL
      use module_model,only    : idrdp,incomp,mach0,comp
!
! 1. Set boundary condition for viscous flux
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)    :: IMODE
      INTEGER,INTENT(IN)    :: MAT_NO(0:MXMAT)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF( MXSSFBC)
      integer,intent(in)    :: LCYCSF ( MXSSFBC)      
      integer,intent(in)    :: LFUTAU (    MXCV)
      logical,INTENT(IN)    :: mat_cal( 0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: FRSTCV(  MXSSFBC)
      real*8 ,intent(in)    :: vel   (  MXALLCV,3)
      real*8 ,intent(in)    :: rho   (  MXALLCV)
      real*8 ,intent(in)    :: rmu   (  MXALLCV)
      real*8 ,intent(inout) :: utau  (0:MXSSFBC)
      real*8 ,intent(inout) :: dflx  (  MXCVFAC,3)
      real*8 ,intent(inout) :: rmue  (  MXALLCV)
      real*8 ,intent(in)    :: rmut  (  MXALLCV)
      real*8 ,intent(inout)    :: dsclv     (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: ux,uy,uz,up,up1,yp,rnu,tau,utaul,ux1,uy1,uz1
      real*8  :: rvaicf,rvaicfp,dum1,yplus,yplus_s=60.d0,gamm
      
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,iflag
      integer :: IMAT,IIMAT,IDIM
      integer :: nb,kd,ISLD
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!------------------------------------------------
! --- dflx 
!------------------------------------------------
      if(IMODE.eq.2) then 
!--------------------------------------
! --- Calculation of effec. mu
!--------------------------------------
        do 1000 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
!        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle 
        kd=kdbcnd(0,nb)
        kdv=kdbcnd(1,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd==kdshutr.or.kd==kdbuff) cycle
        if(kdv.eq.kvfslp) then
!--------------------------------------
! --- / free slip /
!--------------------------------------
          do 1100 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          dflx(ICFL,1)=0.d0
          dflx(ICFL,2)=0.d0
          dflx(ICFL,3)=0.d0
 1100     continue
          if(kd.eq.kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            dflx(ICFP,1)=0.d0
            dflx(ICFP,2)=0.d0
            dflx(ICFP,3)=0.d0
            enddo
          endif
!        elseif(kdv==kvEWF) then
!          do IBFL=IBFS,IBFE
!          utaul=utau(IBFL)
!          yp=FRSTCV(IBFL)
!          yplus=utaul*yp*rho(ICV)/rmu(ICV)
!          gamm=-0.01d0*yplus**4/(1.d0+5.d0*yplus)
!          uplus_lam=
!          uplus_trb=
!          enddo
        elseif(kdv==kvlglw.or.kdv==kvrogf.or.kdv==kvsataW) then
!--------------------------------------
! --- / log-law wall / 
!--------------------------------------
          do 1200 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            ux=vel(ICV,1)-vel(IDC,1)
            uy=vel(ICV,2)-vel(IDC,2)
            uz=vel(ICV,3)-vel(IDC,3)
            up=ux*SFAREA(1,ICFL)
     &        +uy*SFAREA(2,ICFL)
     &        +uz*SFAREA(3,ICFL)
            ux=ux-up*SFAREA(1,ICFL)
            uy=uy-up*SFAREA(2,ICFL)
            uz=uz-up*SFAREA(3,ICFL)
            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
            ux=ux/up
            uy=uy/up
            uz=uz/up
            utaul=utau(IBFL)
            tau=rho(ICV)*utaul*utaul
            yp=FRSTCV(IBFL)
            dsclv(IDC)=rmue(IDC)/yp
            dflx(ICFL,1)=-tau*ux
            dflx(ICFL,2)=-tau*uy
            dflx(ICFL,3)=-tau*uz
 1200     enddo
        
! --- / log-law wall for interface BC/
          if(kd==kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            ICV=LVEDGE(1,ICFP)
            IDC=LVEDGE(2,ICFP)
            ux=vel(ICV,1)-vel(IDC,1)
            uy=vel(ICV,2)-vel(IDC,2)
            uz=vel(ICV,3)-vel(IDC,3)
            up=ux*SFAREA(1,ICFP)+uy*SFAREA(2,ICFP)+uz*SFAREA(3,ICFP)
            ux=ux-up*SFAREA(1,ICFP)
            uy=uy-up*SFAREA(2,ICFP)
            uz=uz-up*SFAREA(3,ICFP)
            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
            ux=ux/up
            uy=uy/up
            uz=uz/up
            utaul=utau(IBFL)
            yp=FRSTCV(IBFL)
            dsclv(IDC)=rmue(IDC)/yp
            tau=rho(ICV)*utaul*utaul
            dflx(ICFP,1)=-tau*ux
            dflx(ICFP,2)=-tau*uy
            dflx(ICFP,3)=-tau*uz
            enddo
          endif
        elseif(kdv==kvnslp) then
!-------------------
! --- / NO-slip /
!-------------------
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ux=vel(ICV,1)-vel(IDC,1)
          uy=vel(ICV,2)-vel(IDC,2)
          uz=vel(ICV,3)-vel(IDC,3)
!
          up=dsqrt(ux*ux+uy*uy+uz*uz+SML)
          ux=ux/up
          uy=uy/up
          uz=uz/up
!
!          ux=vel(ICV,1)-vel(IDC,1)
!          uy=vel(ICV,2)-vel(IDC,2)
!          uz=vel(ICV,3)-vel(IDC,3)
!          up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
!          ux=ux-up*SFAREA(1,ICFL)
!          uy=uy-up*SFAREA(2,ICFL)
!          uz=uz-up*SFAREA(3,ICFL)
!          up=dsqrt(ux*ux+uy*uy+uz*uz+SML)
!          ux=ux/up
!          uy=uy/up
!          uz=uz/up
!
          yp=FRSTCV(IBFL)
!!!!!          tau=rmue(IDC)*up/yp  !"dens==4", "comp"
          tau=rmue(ICV)*up/yp
          dsclv(IDC)=rmue(IDC)/yp
          dflx(ICFL,1)=-tau*ux
          dflx(ICFL,2)=-tau*uy
          dflx(ICFL,3)=-tau*uz
          enddo
! --- No-slip for interface BC
          if (kd.eq.kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            ICV=LVEDGE(1,ICFP)
            IDC=LVEDGE(2,ICFP)
            ux=vel(ICV,1)-vel(IDC,1)
            uy=vel(ICV,2)-vel(IDC,2)
            uz=vel(ICV,3)-vel(IDC,3)
            up=ux*SFAREA(1,ICFP)+uy*SFAREA(2,ICFP)+uz*SFAREA(3,ICFP)
            ux=ux-up*SFAREA(1,ICFP)
            uy=uy-up*SFAREA(2,ICFP)
            uz=uz-up*SFAREA(3,ICFP)
            up=dsqrt(ux*ux+uy*uy+uz*uz)+SML
            ux=ux/up
            uy=uy/up
            uz=uz/up
            yp=FRSTCV(IBFL)
            dsclv(IDC)=rmue(IDC)/yp
            tau=rmue(IDC)*up/yp
            dflx(ICFP,1)=-tau*ux
            dflx(ICFP,2)=-tau*uy
            dflx(ICFP,3)=-tau*uz
            enddo
          endif
        elseif(kdv==kvmodl) then
          if(icaltb==SDES
     &   .or.icaltb==SLES
     &   .or.icaltb==dles
     &   .or.icaltb==ke
     &   .or.icaltb==RNG
     &   .or.icaltb==KE2S
     &   .or.icaltb==CHEN) then
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
            ux=ux/up
            uy=uy/up
            uz=uz/up
            yp=FRSTCV(IBFL)
            dsclv(IDC)=rmue(IDC)/yp
            tau=rmue(ICV)*up/yp  !
            dflx(ICFL,1)=-tau*ux
            dflx(ICFL,2)=-tau*uy
            dflx(ICFL,3)=-tau*uz
            enddo
          else
            call FFRABORT
     &    (1,'ERR: wall BC [model] is NOT supported for DNS')
          endif
        endif
        IF(kd==kdsld) then !zhangsld
          cycle
          do 1400 ISLD=1,2
          IF(ISLD==1) THEN
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_pair(nb)
          ELSE
            IBFS=LBC_pair(nb)+1
            IBFE=LBC_INDEX(nb)
          ENDIF
          do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            do IDIM=1,3
            rvaicf=dflx(ICFL,IDIM)
            rvaicfp=dflx(ICFP,IDIM)
            dum1=0.5d0*(rvaicf-rvaicfp)
            dflx(ICFL,IDIM)=dum1
            enddo
          enddo
 1400     enddo
        ENDIF
        IF(kd==kdovst) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          dflx(ICFL,:)=0.d0
          enddo
        ENDIF
!        if(kd==kdsld.and.ical_sld==1) then
!          do 1200 ISLD=1,2
!          IF(ISLD==1) THEN
!            if(idis(nb)==0) then
!              IBFS=LBC_INDEX(nb-1)+1
!              IBFE=IBFS+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
!            elseif(idis(nb)==1) then
!              IBFS=LBC_INDEX(nb-1)+1
!              IBFE=LBC_pair(nb)
!            endif
!          ELSE
!            if(idis(nb)==0) then
!              IBFS=LBC_INDEX(nb-1)+1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2
!              IBFE=LBC_INDEX(nb)
!            elseif(idis(nb)==1) then
!              IBFS=LBC_pair(nb)+1
!              IBFE=LBC_INDEX(nb)
!            endif
!          ENDIF
!----------------------------------------------------------
!          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
!          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
!          if(ishaft(IMATS(ISLD))==1) then
!            th(ISLD)=rot_ang(IMATS(ISLD))
!            unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
!            unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
!            unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
!            CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
!            do i=1,3
!            do j=1,3
!            rbb(i,j,ISLD)=bb(j,i,ISLD)
!            enddo
!            enddo
          
!            DO IBFL=IBFS,IBFE
!            ICFL=LBC_SSF(IBFL)
!            ICV=LVEDGE(1,ICFL)
!            IDC=LVEDGE(2,ICFL)
!            radius(1)=CVCENT(1,IDC)-begin(1,IMATS(ISLD))
!            radius(2)=CVCENT(2,IDC)-begin(2,IMATS(ISLD))
!            radius(3)=CVCENT(3,IDC)-begin(3,IMATS(ISLD))
!            dum3=radius(1)*unit(1)
!     &          +radius(2)*unit(2)
!     &          +radius(3)*unit(3)
!            radius(1)=radius(1)-dum3*unit(1)
!            radius(2)=radius(2)-dum3*unit(2)
!            radius(3)=radius(3)-dum3*unit(3)
!            velref(:)=vel(IDC,:)
!            call cal_AXBEQC(unit,velref,velref,dum2)
!! --- force=[N/m^3]
!            force(1)=-2.d0*rho(IDC)*velref(1)*rotati(IMATS(ISLD))
!     &                    +rho(IDC)*radius(1)*rotati(IMATS(ISLD))**2
!            force(2)=-2.d0*rho(IDC)*velref(2)*rotati(IMATS(ISLD))
!     &                    +rho(IDC)*radius(2)*rotati(IMATS(ISLD))**2
!            force(3)=-2.d0*rho(IDC)*velref(3)*rotati(IMATS(ISLD))
!     &                    +rho(IDC)*radius(3)*rotati(IMATS(ISLD))**2
!            enddo
!          endif
! 1200     enddo
!        endif
 1000   enddo
!
        do nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        kd=kdbcnd(0,nb)
        kdv=kdbcnd(1,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd.eq.kdshutr.or.kd==kdbuff) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            ICFL=LBC_SSF(IBFL)
!            if(SHUTFL(IBFL)==0) then
!            else
              if(kdv==kvlglw.or.kdv==kvrogf.or.kdv==kvsataW) then
              elseif(kdv.eq.kvfslp) then


                dflx(ICFL,1)=0.d0
                dflx(ICFL,2)=0.d0
                dflx(ICFL,3)=0.d0
                dflx(ICFP,1)=0.d0
                dflx(ICFP,2)=0.d0
                dflx(ICFP,3)=0.d0
              elseif(kdv.eq.kvnslp) then
              elseif(kdv==kvmodl) then
              else

              endif
!            endif
            enddo
        endif
        enddo
!------------------------------------------------
! --- rmue
!------------------------------------------------
      elseif(IMODE.eq.1) then 
        do 2000 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        kd=kdbcnd(0,nb)
        kdv=kdbcnd(1,nb)
        if(kd==kdshutr.or.kd==kdbuff) cycle
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kdv.eq.kvfslp) then
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- / free slip /
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          rmue(IDC)=0.d0
          enddo
        elseif(kdv.eq.kvlglw.or.kdv==kvrogf.or.kdv==kvsataW) then
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- / log-law wall /
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          do 2200 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          yp=FRSTCV(IBFL)
          ux=vel(ICV,1)-vel(IDC,1)
          uy=vel(ICV,2)-vel(IDC,2)
          uz=vel(ICV,3)-vel(IDC,3)
          up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
          ux=ux-up*SFAREA(1,ICFL)
          uy=uy-up*SFAREA(2,ICFL)
          uz=uz-up*SFAREA(3,ICFL)
          up=dsqrt(ux*ux+uy*uy+uz*uz)+SML 
          utaul=utau(IBFL)
          yplus=yp*utaul*rho(ICV)/(rmu(ICV))
          if(yplus.lt.11.6d0) then 
            rmue(IDC)=rmu(IDC)
          else
            rmue(IDC)=utaul*utaul*rho(IDC)*yp/(up)
          endif
 2200     continue
          if (kd.eq.kdintr) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            ICV=LVEDGE(1,ICFP)
            IDC=LVEDGE(2,ICFP)
            yp=FRSTCV(IBFL)
            ux=vel(ICV,1)-vel(IDC,1)
            uy=vel(ICV,2)-vel(IDC,2)
            uz=vel(ICV,3)-vel(IDC,3)
            up=ux*SFAREA(1,ICFP)+uy*SFAREA(2,ICFP)+uz*SFAREA(3,ICFP)
            ux=ux-up*SFAREA(1,ICFP)
            uy=uy-up*SFAREA(2,ICFP)
            uz=uz-up*SFAREA(3,ICFP)
            up=sqrt(ux*ux+uy*uy+uz*uz)+SML
            utaul=utau(IBFL)
            yplus=yp*utaul*rho(ICV)/(rmu(ICV))
            if(yplus.lt.11.6d0) then 
              rmue(IDC)=rmu(IDC)
            else
              rmue(IDC)=utaul*utaul*rho(IDC)*yp/(up)
            endif
            enddo
          endif
        elseif(kdv.eq.kvnslp) then
!^^^^^^^^^^^^^^^^
! --- No-slip 
!^^^^^^^^^^^^^^^^
          do 2300 IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICV=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)
            rmue(IDC)=rmu(IDC)
 2300     continue
          if (kd.eq.kdintr) then  !.and..false.) then
            do IBFL=IBFS,IBFE
            ICFP=LCYCSF(IBFL)
            ICV=LVEDGE(1,ICFP)
            IDC=LVEDGE(2,ICFP)
            rmue(IDC)=rmu(IDC)
            enddo
          endif
        endif
 2000   enddo
      endif
!--------------------------------
! --- diffusion flux  zero 
!--------------------------------
      if(IMODE.eq.2.or.IMODE.eq.1) then
        if(ical_tetra==1) then
          do nb=1,nbcnd
          IIMAT=MAT_BCIDX(nb,1)
          if(.not.mat_cal(IIMAT)) cycle
          kd=kdbcnd(0,nb)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(
!     &       kd.eq.kdolet.or.
!!!!     &     kd.eq.kdilet.or.   !nikki
     &       kd.eq.kdsymm.or.  
!     &       kd.eq.kdstag.or.  
     &       kd.eq.kdtchi
     &    ) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            dflx(ICFL,1)=0.d0
            dflx(ICFL,2)=0.d0
            dflx(ICFL,3)=0.d0
            enddo
            endif
          enddo
        else
          do 325 nb=1,nbcnd
          IIMAT=MAT_BCIDX(nb,1)
          if(.not.mat_cal(IIMAT)) cycle
          kd=kdbcnd(0,nb)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          if(
     &       kd.eq.kdolet.or.
!!!     &     (kd.eq.kdilet.and.idrdp.eq.comp).or.   !nikki
     &       kd.eq.kdsymm.or.  
!     &       kd.eq.kdstag.or.  
     &       kd.eq.kdtchi) then
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            dflx(ICFL,1)=0.d0
            dflx(ICFL,2)=0.d0
            dflx(ICFL,3)=0.d0
            enddo
            endif
 325      enddo
        endif
!
      endif
!
      return
!
      end subroutine bc_veldif
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_wallke_VECT(big,LVEDGE,LBC_SSF,LCYCSF,
     &   mat_cal,MAT_NO,
     &   SFAREA,SFCENT,CVCENT,CVVOLM,rmu,utau1,rho,
     &   vel,aks,aksrc,diagk)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_constant
      use module_boundary,only : kklglw,kdintr,kvfslp,kdilet,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX
      use module_rans,only     : cmu
      use module_turbparm,only : akappa,awfnc
      use module_scalar,only   : ike
      use module_model,only    : icaltb,ke,ke_low,RNG,CHEN,KE2S
      use module_metrix,only   : tmpfac=>d2vect
!
! 1. Set k & epsilon of cell adjucent to log-law wall
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)  :: big
      integer,intent(in)  :: LVEDGE (2, MXCVFAC)
      integer,intent(in)  :: LBC_SSF (  MXSSFBC)
      integer,intent(in)  :: LCYCSF  (  MXSSFBC)      
      logical,intent(in)  :: mat_cal(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO(      0:MXMAT)
      real*8 ,intent(in)  :: SFAREA (4, MXCVFAC)
      real*8 ,intent(in)  :: SFCENT (3, MXCVFAC)
      real*8 ,intent(in)  :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)  :: CVVOLM (   MXALLCV)
      real*8 ,intent(in)  :: rmu    (   MXALLCV)
      real*8 ,intent(in)  :: rho    (   MXALLCV)
      real*8 ,intent(in)  :: vel    (   MXALLCV,3)
      real*8 ,intent(in)  :: aks(       MXALLCVR,MXrans)
      real*8 ,intent(out) :: aksrc(     MXALLCVR,MXrans)
      real*8 ,intent(inout) :: diagk  ( MXALLCVR,MXrans)
      real*8 ,intent(in)    :: utau1 (0:MXSSFBC)
!
! --- [local entities]
!
      real*8  :: ux,uy,uz,yp,yp2,
     &           rnu,uplus,up,utau,dumy,rans_k,rans_e
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!--------------------------------
! --- High Reynolds model
!--------------------------------
      do 1000 nb=1,nbcnd   !IBF=1,nssfbc 
      IIMAT=MAT_BCIDX(nb,1)
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT).or.IMAT<0) goto 1000
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      kdv=kdbcnd(1,nb) 
      kd=kdbcnd(0,nb) 
      if((icaltb==ke.or.icaltb==RNG.or.icaltb==CHEN.or.icaltb==KE2S)
     &   .and.kdbcnd(4,nb).eq.kklglw) then 
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
!
        ux=vel(ICV,1)-vel(IDC,1)
        uy=vel(ICV,2)-vel(IDC,2)
        uz=vel(ICV,3)-vel(IDC,3)
        up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL)
        ux=ux-up*SFAREA(1,ICFL)
        uy=uy-up*SFAREA(2,ICFL)
        uz=uz-up*SFAREA(3,ICFL)
        up=dsqrt(ux*ux+uy*uy+uz*uz)
        yp2=dsqrt((SFCENT(1,ICFL)-CVCENT(1,ICV))**2
     &           +(SFCENT(2,ICFL)-CVCENT(2,ICV))**2
     &           +(SFCENT(3,ICFL)-CVCENT(3,ICV))**2)
        yp=CVVOLM(ICV)**(1.d0/3.d0)
!        yp=max(1.d-6,min(yp,yp2))       !(1)
        utau=cmu**0.25d0*dsqrt(aks(ICV,ike(1)))
!        utau=utau1(IBFL)
        aksrc(ICV,ike(1))=(utau*utau/sqrt(cmu)-aks(ICV,ike(1)))*big
        aksrc(ICV,ike(2))=(utau**3/(akappa*yp+SML)-aks(ICV,ike(2)))*big


!        aksrc(ICV,ike(1))=(utau*utau/sqrt(cmu))*big
!        aksrc(ICV,ike(2))=(utau**3/(akappa*yp+SML))*big

!
        
!        aksrc(ICV,ike(1))=0.d0
!        aksrc(ICV,ike(2))=(utau**3/(akappa*yp+SML)-aks(ICV,ike(2)))*big
!
        diagk(ICV,ike(1))=big/rho(ICV)
        diagk(ICV,ike(2))=big/rho(ICV)
        enddo
      endif
      if((icaltb==ke.or.icaltb==RNG.or.icaltb==CHEN.or.icaltb==KE2S)
     &    .and.(kd==kdilet))
     &   then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        aksrc(ICV,ike(1))=(aks(IDC,ike(1))-aks(ICV,ike(1)))*big
        aksrc(ICV,ike(2))=(aks(IDC,ike(2))-aks(ICV,ike(2)))*big
        diagk(ICV,ike(1))=big/rho(ICV)
        diagk(ICV,ike(2))=big/rho(ICV)
        enddo
      endif
!----------------------------
! --- low Reynolds model
!----------------------------
      if(icaltb.eq.ke_low.and.kdbcnd(4,nb).eq.kklglw) then
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
        up=dsqrt(ux*ux+uy*uy+uz*uz)
        yp2=dsqrt((SFCENT(1,ICFL)-CVCENT(1,ICV))**2
     &           +(SFCENT(2,ICFL)-CVCENT(2,ICV))**2
     &           +(SFCENT(3,ICFL)-CVCENT(3,ICV))**2)
        yp=CVVOLM(ICV)**(1.d0/3.d0)
        yp=max(1.d-6,min(yp,yp2))
        rnu=rmu(ICV)/rho(ICV)
        aksrc(ICV,ike(2))=(2.d0*rnu*aks(ICV,ike(1))/yp**2
     &    -aks(ICV,ike(2)))*big
        diagk(ICV,ike(2))=big/rho(ICV)
        enddo
      endif
!----------------------------
 1000 continue
      return
      end subroutine bc_wallke_VECT   !15767
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_wallke(big,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &   SFAREA,SFCENT,CVCENT,CVVOLM,FRSTCV,
     &   rmu,rho,vel,aks,aksrc,diagk)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kklglw,kdintr,kvfslp,kdilet,kdsld,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                           kdolet,kdpres
     &                           ,kdovst
      use module_rans,only     : cmu
      use module_turbparm,only : akappa,awfnc,yplsm,E_vel
      use module_scalar,only   : ike
      use module_model,only    : icaltb,ke,ke_low,RNG,CHEN,ical_EWT,KE2S
!
! 1. Set k & epsilon of cell adjucent to log-law wall
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)  :: big
      integer,intent(in)  :: LVEDGE (2, MXCVFAC)
      integer,intent(in)  :: LBC_SSF (  MXSSFBC)
      integer,intent(in)  :: LCYCSF  (  MXSSFBC)      
      logical,intent(in)  :: mat_cal(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO(      0:MXMAT)
      real*8 ,intent(in)  :: SFAREA (4, MXCVFAC)
      real*8 ,intent(in)  :: SFCENT (3, MXCVFAC)
      real*8 ,intent(in)  :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)  :: CVVOLM (   MXALLCV)
      real*8 ,intent(in)  :: rmu    (   MXALLCV)
      real*8 ,intent(in)  :: rho    (   MXALLCV)
      real*8 ,intent(in)  :: vel    (   MXALLCV,3)
      real*8 ,intent(in)  :: aks(       MXALLCVR,MXrans)
      real*8 ,intent(out) :: aksrc(     MXALLCVR,MXrans)
      real*8 ,intent(inout) :: diagk  ( MXALLCVR,MXrans)
      REAL*8 ,INTENT(IN)  :: FRSTCV(        MXSSFBC)
!
! --- [local entities] 
!
      real*8  :: ux,uy,uz,yp,yp2,dum1,dum2,dum3,dum4,Yv,lll,
     &           rnu,uplus,up,utau,dumy,rans_k,rans_e
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
      real*8,parameter :: Rey_s=200.d0,Amu=70.d0
      real*8  :: Rey,rmd_E,AAA,dRey=0.2d0*Rey_s,lmu,cl,
     &           ddiag_k,src_k
!--------------------------------
! --- High Reynolds model 
!--------------------------------
      do 1000 nb=1,nbcnd   !IBF=1,nssfbc 
      IIMAT=MAT_BCIDX(nb,1)
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT).or.IMAT<0) goto 1000 
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      kdv=kdbcnd(1,nb) 
      kd=kdbcnd(0,nb) 
      if((icaltb==KE2S.or.icaltb==ke.or.icaltb==RNG.or.icaltb==CHEN) 
     &   .and.kdbcnd(4,nb).eq.kklglw) then 
        AAA=abs(dRey/tanh(0.98d0))
        
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
! 
        ux=vel(ICV,1)-vel(IDC,1)
        uy=vel(ICV,2)-vel(IDC,2)
        uz=vel(ICV,3)-vel(IDC,3)
        up=ux*SFAREA(1,ICFL)+uy*SFAREA(2,ICFL)+uz*SFAREA(3,ICFL) 
        ux=ux-up*SFAREA(1,ICFL)
        uy=uy-up*SFAREA(2,ICFL)
        uz=uz-up*SFAREA(3,ICFL)
        
        up=dsqrt(ux*ux+uy*uy+uz*uz)
        yp=FRSTCV(IBFL)
        rnu=rmu(ICV)/rho(ICV) 
!
!        yp=CVVOLM(ICV)**(1.d0/3.d0) 
!        if(kdv==kvfslp) then 
!          utau=cmu**0.25d0*dsqrt(aks(ICV,ike(1)))
!        else 
!          call cal_ustar(yp,rnu,akappa,awfnc,up,utau,2,3)
!        endif 
!-------------------------------------------------------------- 
!0) 
!        rans_k=aks(ICV,ike(1)) 
!        rans_e=aks(ICV,ike(2)) 
!        utau=cmu**0.25d0*dsqrt(rans_k)
!        aksrc(ICV,ike(2))=((cmu**0.75d0*rans_k**1.5d0)
!     &     /(akappa*yp+SML)-rans_e)*big
!        diagk(ICV,ike(2))=big/rho(ICV)
!1-1) 
! 
        if(ical_EWT==1.and..false.) then 
!        if(ical_EWT==1) then 
          IF(.false.) THEN 

          rans_k=aks(ICV,ike(1)) 
!
          utau=cmu**0.25d0*dsqrt(rans_k)
          dum1=utau 
          call cal_ustar(yp,rnu,akappa,awfnc,up,dum1,2,4) 
          utau=dum1   !min(utau,dum1) 
          
          lll=akappa*yp+SML

          rans_e=cmu**0.75d0*rans_k**1.5d0/(lll) 
          ddiag_k=cmu**0.75d0*rans_k**0.5d0/(lll) 
!
!!Tw      
          dum2=utau*yp/rnu*E_vel !Y+*E
          dum1=rho(ICV)*cmu**0.25d0*akappa  !Tw
     &        *dsqrt(rans_k)*up/log(dum2)
          dum2=cmu**0.25d0*dsqrt(rans_k)/(lll) !dVdn 
!dVdn
          dum1=utau**2*rho(ICV)   !Tw
!
          src_k=dum1*dum2
!          src_k=utau**3/lll*rho(ICV)
!
!          if(yp<yv) then 

          aksrc(ICV,ike(1))=src_k-rans_e*rho(ICV)
          diagk(ICV,ike(1))=ddiag_k*rho(ICV) 

!          else
!          dum2=utau**3/(akappa*yp+SML)
!          aksrc(ICV,ike(2))=(dum2-aks(ICV,ike(2)))*big
!          diagk(ICV,ike(2))=big/rho(ICV)
!          dum1=utau*utau/sqrt(cmu)   
!          aksrc(ICV,ike(1))=(dum1-aks(ICV,ike(1)))*big
!          diagk(ICV,ike(1))=big/rho(ICV)
!          endif

!
! --- --- 
!     
          ENDIF 
        else 
          IF(.true.) THEN

          utau=cmu**0.25d0*dsqrt(aks(ICV,ike(1))) 

          dum1=utau
          call cal_ustar(yp,rnu,akappa,awfnc,up,dum1,2,5) 
          utau=dum1  !min(utau,dum1)


          dum1=utau*utau/sqrt(cmu)   
          dum2=utau**3/(akappa*yp+SML)
!
          aksrc(ICV,ike(1))=(dum1-aks(ICV,ike(1)))*big
          aksrc(ICV,ike(2))=(dum2-aks(ICV,ike(2)))*big
          diagk(ICV,ike(1))=big/rho(ICV)
          diagk(ICV,ike(2))=big/rho(ICV)
          
          endif
        endif
!
!1) --- OK 
!!!!!!!!!!!!!!!!!!!! 
!        utau=cmu**0.25d0*dsqrt(aks(ICV,ike(1))) 
!        call cal_ustar(yp,rnu,akappa,awfnc,up,utau,1,6) ! 
!        dum2=(utau**3/(akappa*yp+SML)-aks(ICV,ike(2)))*big 
!        aksrc(ICV,ike(1))=(utau*utau/sqrt(cmu)-aks(ICV,ike(1)))*big 
!        aksrc(ICV,ike(2))=dum2
!        diagk(ICV,ike(1))=big/rho(ICV) 
!        diagk(ICV,ike(2))=big/rho(ICV) 
!2) --- 
!        utau=cmu**0.25d0*dsqrt(aks(ICV,ike(1)))
!        call cal_ustar(yp,rnu,akappa,awfnc,up,utau,2,7)
!        aksrc(ICV,ike(2))=((cmu**0.75d0*aks(ICV,ike(1))**1.5d0)
!     &     /(akappa*yp+SML)-aks(ICV,ike(2)))*big
!        aksrc(ICV,ike(1))=(utau*utau/sqrt(cmu)-aks(ICV,ike(1)))*big
!        diagk(ICV,ike(1))=big/rho(ICV)
!        diagk(ICV,ike(2))=big/rho(ICV)
!3) --- 0) 
!        utau=cmu**0.25d0*dsqrt(aks(ICV,ike(1)))
!        call cal_ustar(yp,rnu,akappa,awfnc,up,utau,2,8)
!        aksrc(ICV,ike(2))=((cmu**0.75d0*aks(ICV,ike(1))**1.5d0)
!     &     /(akappa*yp+SML)-aks(ICV,ike(2)))*big
!        diagk(ICV,ike(2))=big/rho(ICV)
!4) --- 
!        utau=cmu**0.25d0*dsqrt(aks(ICV,ike(1)))
!        call cal_ustar(yp,rnu,akappa,awfnc,up,utau,2,9)
!        aksrc(ICV,ike(2))=(utau**3/(akappa*yp+SML)-aks(ICV,ike(2)))*big
!        diagk(ICV,ike(2))=big/rho(ICV)
! -------------
        enddo
      endif
!
      if((icaltb==KE2S.or.icaltb==ke.or.icaltb==RNG.or.icaltb==CHEN)
     & .and.(kd==kdsld.OR.kd==kdovst))
     &  then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          dum1=0.d0
          dum2=0.1d0
          aksrc(ICV,ike(1))=(dum1-aks(ICV,ike(1)))*big 
          aksrc(ICV,ike(2))=(dum2-aks(ICV,ike(2)))*big 
          diagk(ICV,ike(1))=big/rho(ICV) 
          diagk(ICV,ike(2))=big/rho(ICV) 
          enddo
      endif

      if(icaltb==KE2S.or.icaltb==ke.or.icaltb==RNG.or.icaltb==CHEN) THEN 
        IF(kd==kdilet) THEN 
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          aksrc(ICV,ike(1))=(aks(IDC,ike(1))-aks(ICV,ike(1)))*big
          aksrc(ICV,ike(2))=(aks(IDC,ike(2))-aks(ICV,ike(2)))*big
          diagk(ICV,ike(1))=big/rho(ICV)
          diagk(ICV,ike(2))=big/rho(ICV)
          enddo
        ELSEIF(kd==kdolet.or.kd==kdpres) THEN
          do IBFL=IBFS,IBFE !watanabe
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=0.d0 
          dum2=0.d0 
          aksrc(ICV,ike(1))=(dum1-aks(ICV,ike(1)))*big 
          aksrc(ICV,ike(2))=(dum2-aks(ICV,ike(2)))*big 
          diagk(ICV,ike(1))=big/rho(ICV)
          diagk(ICV,ike(2))=big/rho(ICV)
          enddo
        ENDIF
      endif
!----------------------------
! --- low Reynolds model
!----------------------------
      if(icaltb.eq.ke_low.and.kdbcnd(4,nb).eq.kklglw) then
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
        up=dsqrt(ux*ux+uy*uy+uz*uz)
        yp2=abs(SFAREA(1,ICFL)*(SFCENT(1,ICFL)-CVCENT(1,ICV))
     &         +SFAREA(2,ICFL)*(SFCENT(2,ICFL)-CVCENT(2,ICV))
     &         +SFAREA(3,ICFL)*(SFCENT(3,ICFL)-CVCENT(3,ICV)))
!        yp=dsqrt((SFCENT(1,ICFL)-CVCENT(1,ICV))**2
!     &          +(SFCENT(2,ICFL)-CVCENT(2,ICV))**2
!     &          +(SFCENT(3,ICFL)-CVCENT(3,ICV))**2)
        yp=CVVOLM(ICV)**(1.d0/3.d0)
        yp=max(1.d-6,min(yp,yp2))
        rnu=rmu(ICV)/rho(ICV)
        aksrc(ICV,ike(2))=(2.d0*rnu*aks(ICV,ike(1))/yp**2
     &    -aks(ICV,ike(2)))*big
        diagk(ICV,ike(2))=big/rho(ICV)
        enddo
      endif
!----------------------------
 1000 continue
!
      return
      end subroutine bc_wallke
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine USER_SHUTTER(iter,time,dt,BCNO,MXSSFBCX,
     &            MXSSFBC1,MXCVFACX,IBFS,IBFE,
     &            SHUTFL,SFCENT,LBC_SSF)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none  
!
! --- [dummy arguments]
!
      real*8 ,intent(in)   :: time,dt
      integer,intent(in)   :: iter
      integer,intent(in)   :: BCNO,MXSSFBCX,MXSSFBC1,MXCVFACX
      integer,intent(in)   :: IBFS,IBFE
      real*8 ,intent(in)   :: SFCENT(3,MXCVFACX)
      integer,intent(in)   :: LBC_SSF(MXSSFBC1)
      integer,intent(out)  :: SHUTFL(MXSSFBCX)
!
!--------------------------------------------------------------------- 
!
! --- [local entities] 
!
      real*8,parameter :: shudown=0.02d0,door_H=2.d0
      real*8,parameter :: shu_spd=door_H/shudown
      real*8 :: shu_start=0.05d0
      real*8  :: dum1,dum2,dum3,x(3)
      integer :: idum1,IBFL,ICFL

      SHUTFL(IBFS:IBFE)=0   !open
      return
      if(time>shu_start) then
        dum2=time-shu_start
        dum1=dum2*shu_spd
        dum3=door_H-dum1
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        x(1)=SFCENT(1,ICFL)
        x(2)=SFCENT(2,ICFL)
        x(3)=SFCENT(3,ICFL)
!        if(x(3)>dum3) then
        if(x(2)>dum3) then
          SHUTFL(IBFL)=1
        endif
        enddo
      endif

      return
      end subroutine USER_SHUTTER



!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_SLD(nko,LVEDGE,LBC_SSF,LCYCSF,LCYCOLD,
     &                     mat_cal,aa)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments] 
!
      use module_dimension
      use module_boundary,only : kdbcnd,kdprdc,LBC_INDEX,nbcnd,kdtchi,
     &                           MAT_BCIDX,kdintr,kdsld,kdcvd,idis,
     &                           kdbuff,kdshutr,lbc_pair
     &                           ,kdovst
!      use module_model,only    : ical_vect
!      use module_metrix,only   : tmpfac=>d2vect
!      use module_metrix,only   : SHUTFL
!
! 1. Set scalar at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(inout) :: aa     ( MXALLCV,nko,nko)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,I
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC
      integer :: ISLD,ISLD2
!
! --- 
!
      do nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc) then
      elseif(kd==kdbuff) then
      elseif(kd==kdshutr) then
      elseif(kd==kdsld) then
        do ISLD=1,2
        IF(ISLD==1) THEN
          ISLD2=2
          if(idis(nb)==0) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=IBFS+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
          elseif(idis(nb)>=1) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_pair(nb)
          endif
        ELSE
          ISLD2=1
          if(idis(nb)==0) then
            IBFS=LBC_INDEX(nb-1)+1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2
            IBFE=LBC_INDEX(nb)
          elseif(idis(nb)>=1) then
            IBFS=LBC_pair(nb)+1
            IBFE=LBC_INDEX(nb)
          endif
        ENDIF
!
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
!        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
        aa(ICV,1:nko,1:nko)=0.d0
        enddo
        enddo
      else
      endif
      enddo
!
      return
!
      end subroutine bc_SLD


