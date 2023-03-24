!
!    subroutine src_rho
!    subroutine src_setxxx
!    subroutine keps_src
!    subroutine cavit_src
!    subroutine C_dash_src
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine src_rho(MXALLCV,rsrc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_source,only : nscnd,kdscnd,kdmass,
     &                        icscnd,lcscnd,wdscnd
!
! 1. Set source term for density
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: MXALLCV
      real*8 ,intent(out) :: rsrc(MXALLCV)
!
! --- [local entities]
!
      real*8  :: rhox
      integer :: l,m
!
!
!
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
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine src_setxxx(nval,cval,rho,dval,src,tau)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_source,only : nscnd,kdscnd,kdpvol,kdpmas,kdmass,
     & icscnd,lcscnd,iscomp,israns,wdscnd,chkncomp,chknrans
!
!
!
      implicit none
! --- [dummy arguments]
!
      character(*),intent(in)  :: cval
      integer     ,intent(in)  :: nval
      real*8      ,intent(in)  :: rho    (  MXALLCV)
      real*8      ,intent(in)  :: dval(     MXALLCV,nval)
      real*8      ,intent(out) :: src (     MXALLCV,nval)
      real*8      ,intent(out) :: tau    (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: wdt(nval)
      real*8  :: rhox,taux
      integer :: i,l,m,n,kd,ls,le,isval
!
!
!
      src=0.d0
      tau=0.d0
      if( nscnd.lt.1 ) return
!
      if( cval.eq.'vel' ) then
        isval=1
      elseif( cval.eq.'tmp' ) then
        isval=4
      elseif( cval.eq.'yys' ) then
        call chkncomp(nval,'(src_setxxx)')
        isval=iscomp
      elseif( cval.eq.'aks' ) then
        call chknrans(nval,'(src_setxxx)')
        isval=israns
      else
        write(*,*) '### program error -1- (src_setxxx)'
        CALL FFRABORT(1,'src_setxxx')
      endif
!
! --- 
!
      do 1000 m=1,nscnd
      ls=icscnd(m-1)+1
      le=icscnd(m)
      kd=kdscnd(m)
!
      do 100 i=1,nval
      wdt(i)=wdscnd(isval+i,m)
  100 continue
!
!-< 1. source per unit volume >-
!
      if( kd.eq.kdpvol ) then
        do 200 l=ls,le
          n=lcscnd(l)
          do 201 i=1,nval
            src(n,i)=wdt(i)
  201     continue
  200   continue
      endif
!
!-< 2. source per unit mass >-
!
      if( kd.eq.kdpmas ) then
        do 210 l=ls,le
          n=lcscnd(l)
          do 211 i=1,nval
            src(n,i)=rho(n)*wdt(i)
  211     continue
  210   continue
      endif
!
!-< 3. mass source >-
!
      if( kd.eq.kdmass ) then
      rhox=wdscnd(1,m)
      if( rhox.ge.0.d0 ) then
        do 300 i=1,nval
          wdt(i)=rhox*wdt(i)
  300   continue
        do 310 l=ls,le
          n=lcscnd(l)
          do 311 i=1,nval
            src(n,i)=wdt(i)
  311     continue
  310   continue
      else
        taux=-rhox
        do 320 l=ls,le
          n=lcscnd(l)
          tau(n)=taux
          do 321 i=1,nval
            src(n,i)=rhox*dval(n,i)
  321     continue
  320   continue
      endif
      endif
!
 1000 continue
!
      end subroutine src_setxxx
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine keps_src(
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  DISALL,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  grdc,dtau,rmut,rmu,rho,vel,tmp,aks,aksrc,
     &  LCYCOLD,wifsld,utau,LFUTAU,
     &  diagk,vctr)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_gravity,only   : ggg
      use module_rans,only      : cep1,cep2,cep3,cep4,cep5,eta0,beta,
     &                            cmu,Ck,Ce
      use module_rans,only      : akdes,cb1,cb2,sgmDES,cw1,cw2,cw3,cv1,
     &                            cDES,L2s
      use module_turbparm,only  : prturb,akappa,awfnc,yplsm,
     &                            Rey_s,Amu,dRey
      use module_scalar,only    : ike,ides,polycdp,mcd,fitVALUE,ikles
      use module_model,only     : icaltb,ke,ke_low,RNG,CHEN,SDES,
     &                            ical_EWT,KE2S,KLES
      use module_metrix,only    : ROOT_S => W1K7
      use module_metrix,only    : dspf   => W1K8
      use module_metrix,only    : dudx => W1K9
      use module_material,only  : ical_sld
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LVEDGE (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF(   MXSSFBC)
      integer,intent(in)    :: LCYCSF (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV(    MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal ( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4, MXCVFAC)
      real*8 ,intent(in)    :: wiface (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (   MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3, MXCVFAC)
      real*8 ,intent(in)    :: DISALL (   MXCV)
      real*8 ,intent(inout) :: grdc   (   MXALLCV,3,3)
      real*8 ,intent(inout) :: dtau   (   MXALLCV,3)
      real*8 ,intent(in)    :: rmut   (   MXALLCV)
      real*8 ,intent(in)    :: rmu       (   MXALLCV)
      real*8 ,intent(in)    :: rho    (   MXALLCV)
      real*8 ,intent(in)    :: vel    (   MXALLCV,3)
      real*8 ,intent(in)    :: tmp    (   MXALLCV)
      real*8 ,intent(in)    :: aks    (   MXALLCVR,MXrans)
      real*8 ,intent(inout) :: aksrc  (   MXALLCVR,MXrans)
      real*8 ,intent(inout) :: diagk  (   MXALLCVR,MXrans)
      integer,intent(in)    :: LCYCOLD(MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      REAL*8 ,INTENT(IN)    :: utau (0:MXSSFBC)
      integer,intent(in)    :: LFUTAU (    MXCV )
!
! --- [local entities]
!
      real*8,parameter :: r2p3=2.d0/3.d0,SML1=1.d-10   !
      real*8,parameter :: ct1=1.d0,ct2=2.d0,ct3=1.1d0,ct4=2.d0
      real*8  :: vv,hh,gtx,rep,wep,prod,prod1,rkk,dum1,dum2,
     &           dum3,yp,CDX,yplus,del
      real*8 ,parameter :: r1pn=1.d0/3.d0,R26=1.D0/26.D0
      real*8  :: Rey,rmd_E,lmu,cl,rnu,AAA,CD1,CD2,lll,fmu

      real*8  :: SIJ11,SIJ22,SIJ33,SIJ12,SIJ13,SIJ23,dspd,rans_k,rans_e
      real*8  :: f2,Rt,rng_trm,eta,chen_trm,R_fact,C_fact,Yprod1,yv
      real*8  :: fv1des,fv2des,Scal,Ddis,Stud,rdes,gdes,fwdes,dens,Pprod
     &          ,Yprod,dume1,rans_nu,xxx,ft1,ft2,gt,fv3des,cv2,nu
      integer :: i,j,k,l,m,n,IBFL
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
!
!-< 1. Calculate dissipation function >-
!
!---------------------------------------------------------------------
! --- production term (grad form):P=2SijSij=2dui/dxj*Sij
!---------------------------------------------------------------------
!
      call grad_cell(3,19,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vel,grdc)
!
      if(ical_sld/=0) then
        call bc_SLD(3,LVEDGE,LBC_SSF,LCYCSF,LCYCOLD,
     &                     mat_cal,grdc)
      endif
!
      if(icaltb==SDES.or.icaltb==KLES) then
        do 600 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO 630 ICVL=ICVS,ICVE
        SIJ12=GRDC(ICVL,1,2)-GRDC(ICVL,2,1)
        SIJ13=GRDC(ICVL,1,3)-GRDC(ICVL,3,1)
        SIJ23=GRDC(ICVL,2,3)-GRDC(ICVL,3,2)
        dspd=dsqrt(2.d0*(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
        ROOT_S(ICVL)=dspd
        SIJ11=grdc(ICVL,1,1)
        SIJ22=GRDC(ICVL,2,2)
        SIJ33=GRDC(ICVL,3,3)
        SIJ12=GRDC(ICVL,1,2)+GRDC(ICVL,2,1)
        SIJ13=GRDC(ICVL,1,3)+GRDC(ICVL,3,1)
        SIJ23=GRDC(ICVL,2,3)+GRDC(ICVL,3,2)
        dspd=(2.d0*SIJ11*SIJ11
     &       +2.d0*SIJ22*SIJ22
     &       +2.d0*SIJ33*SIJ33
     &       +(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
        dspf(ICVL)=dspd  !aks(ICVL,ides)
 630    enddo
 600    enddo
      else
        do 400 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle 
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO 430 ICVL=ICVS,ICVE
        SIJ11=GRDC(ICVL,1,1)
        SIJ22=GRDC(ICVL,2,2)
        SIJ33=GRDC(ICVL,3,3)
        SIJ12=GRDC(ICVL,1,2)+GRDC(ICVL,2,1)
        SIJ13=GRDC(ICVL,1,3)+GRDC(ICVL,3,1)
        SIJ23=GRDC(ICVL,2,3)+GRDC(ICVL,3,2)
        dspd=(2.d0*SIJ11*SIJ11                !dspd=>P
     &       +2.d0*SIJ22*SIJ22
     &       +2.d0*SIJ33*SIJ33
     &       +(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
        dum1=SIJ11+SIJ22+SIJ33
        rans_e=aks(ICVL,ike(2))
		rans_k=aks(ICVL,ike(1))         ! add debag  OSHIMA 2021.10.25 by NuFD 
        rkk=r2p3*(rho(ICVL)*rans_k+rmut(ICVL)*dum1)*dum1
        dum2=cep4*dum1*rho(ICVL)*rans_e !/cep1
        dspf(ICVL)=dspd*rmut(ICVL)-rkk   !+dum2
        ROOT_S(ICVL)=dspd  !dum2                !dspd
        dudx(ICVL)=dum2
 430    enddo
 400    enddo
      endif
!
!
!---------------------------------------------------------------------
! --- production term (intergration form) 
!---------------------------------------------------------------------
!      rk=0.d0
!      call cal_dspfnc(
!     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
!     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
!     &  LCYCOLD,wifsld,vctr,
!     &  grdc,dtau,rmut,rk,vel,dspf)
!
!----------------------------------------------
!-< 2. Calculate temperature gradient >-
!----------------------------------------------
      call grad_cell(1,20,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,tmp,grdc)
!
!-< 3. Calculate source term >-
!
      if(icaltb==KE2S.or.icaltb==ke.or.
     &   icaltb==RNG.or.icaltb==CHEN) then
        if(icaltb==ke) then
           R_fact=0.d0
           C_fact=0.d0
        elseif(icaltb==RNG) then
           R_fact=1.d0
           C_fact=0.d0
        elseif(icaltb==CHEN) then
           R_fact=0.d0
           C_fact=1.d0
        elseif(icaltb==KE2S) then
           R_fact=0.d0
           C_fact=0.d0
        endif
        wep=cep2+cep2
        do 300 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT<0) cycle 
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        rans_k=aks(ICVL,ike(1))
        rans_e=aks(ICVL,ike(2))
        prod=dspf(ICVL)
        hh=rans_e/(rans_k+SML1)
        gtx=(ggg(1)*grdc(ICVL,1,1)
     &      +ggg(2)*grdc(ICVL,2,1)
     &      +ggg(3)*grdc(ICVL,3,1))
     &      *rmut(ICVL)/(prturb*tmp(ICVL)) 
        rep=rho(ICVL)*rans_e 
! --- RNG 
        eta=dsqrt(ROOT_S(ICVL))*rans_k/(rans_e+SML) 
        rng_trm=
     &    R_fact*cmu*eta**3*(1.d0-eta/eta0)/(1.d0+beta*eta**3) 
! --- CHEN 
!        chen_trm=
!     &  C_fact*ROOT_S(ICVL)**2*rmut(ICVL)**2/rans_k/rho(ICVL)
        chen_trm=
     &  C_fact*cmu*ROOT_S(ICVL)**2*rmut(ICVL)*rans_k/(rans_e+SML)
! --- epson source term 
        if(icaltb==KE2S) then
          rnu=rmu(ICVL)/rho(ICVL) 
          Rey=sqrt(rans_k)*L2s/rnu
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
          cep1=1.44d0*Rey**0.5/(15.d0*CDx**(3.d0/8.d0))
          prod1=cep1*prod*(rans_e/rnu**3)**0.25d0
!          prod1=1.44d0*prod*hh*rans_k/(15.d0*(rnu*rans_e)**0.75+SML)
          dum1=prod1
     &       +(cep3*max(0.d0,gtx)
     &       -(cep2+rng_trm)*rep
     &        )*hh
     &       +dudx(ICVL)
        else
          prod1=cep1*prod
          dum1=(prod1
     &       +cep3*max(0.d0,gtx)
     &       -(cep2+rng_trm)*rep
     &        )*hh
     &       +dudx(ICVL)
     &       +cep5*chen_trm!*hh !zhang5 
        endif
        aksrc(ICVL,ike(2))=aksrc(ICVL,ike(2))+dum1 
! --- k source term 
        aksrc(ICVL,ike(1))=aksrc(ICVL,ike(1))+(prod+gtx-rep)
! --- diag term  
        diagk(ICVL,ike(1))=diagk(ICVL,ike(1))+hh*rho(ICVL) 
        dum2=hh*wep*rho(ICVL) 
        diagk(ICVL,ike(2))=diagk(ICVL,ike(2))+dum2
        enddo
  300   continue
!
        if(ical_EWT==1) then
          AAA=abs(dRey/tanh(0.98d0))
          Cl=akappa/cmu**0.75d0     !Cl 
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          IBFL=LFUTAU(ICVL)
          rans_e=aks(ICVL,ike(2))
          yp=DISALL(ICVL)
          rans_k=aks(ICVL,ike(1)) 
          rnu=rmu(ICVL)/rho(ICVL) 
!Norris & Reynolds S------------------------
!          Rey=yp*sqrt(rans_k)/rnu  !Rey 
!          lll=akappa*yp/cmu**0.75d0
!          fmu=(1.d0-exp(-Rey/Amu))  
!          cl=(1.d0+5.3d0/Rey)
!          dum2=rans_k**1.5d0/(lll+SML)*cl 
!Hassid etc S------------------------------
!          lll=akappa*yp
!          Rey=lll*sqrt(rans_k)/rnu  !Rel 
!          fmu=(1.d0-exp(-Rey/Amu))  
!          CD1=cmu**0.75d0
!          CD2=0.336d0
!          AAA=abs((0.2d0*Rey_s)/tanh(0.98d0))
!          cl=(CD1*fmu+CD2/Rey)
!          dum2=cl*rans_k**1.5d0/(lll+SML)
!Wolfstein F---------------------------- 
          Rey=yp*sqrt(rans_k)/rnu  !Rey 
          fmu=(1.d0-exp(-Rey/Amu))  
          cl=akappa/cmu**0.75 
          lll=yp*cl*(1.d0-exp(-Rey/(cl*2.d0))) 
          dum2=rans_k**1.5d0/(lll+SML) 
!----------------------------- 
          yv=rnu*yplsm/utau(IBFL) 
!          IF(fmu<0.99D0) THEN 
          IF(Rey<200.d0) THEN 
            aksrc(ICVL,ike(2))=(dum2-aks(ICVL,ike(2)))*GREAT 
            diagk(ICVL,ike(2))=GREAT/rho(ICVL)
          ENDIF
          enddo
          enddo
        endif
      elseif(icaltb==ke_low) then
!---------------------------
! --- Low Re number model 
!---------------------------
        wep=cep2+cep2
        do 500 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        rans_k=aks(ICVL,ike(1))
        rans_e=aks(ICVL,ike(2))
        Rt=rans_k*rans_k*rho(ICVL)/(rmu(ICVL)*(rans_e)+SML1)
        nu=rmu(ICVL)/rho(ICVL)
        Rey=DISALL(ICVL)*dsqrt(rans_k)/nu
        f2=1.d0-0.3d0*exp(-Rt*Rt)
        prod=dspf(ICVL)
        prod1=dspf(ICVL)+rmut(ICVL)*1.33*(1.d0-0.3d0*exp(-Rt*Rt))
     &    *(dspf(ICVL)*2.D0*rmu(ICVL)/rmut(ICVL)*rans_k/DISALL(ICVL)**2)
     &    *exp(-0.00375*Rey**2)
        hh=rans_e/(rans_k+SML1)
        gtx=(ggg(1)*grdc(ICVL,1,1)
     &      +ggg(2)*grdc(ICVL,2,1)
     &      +ggg(3)*grdc(ICVL,3,1))
     &     *rmut(ICVL)/(prturb*tmp(ICVL))
        rep=rho(ICVL)*rans_e
!
        aksrc(ICVL,ike(2))=aksrc(ICVL,ike(2))+
     &      (cep1*prod1+cep3*max(0.d0,gtx)-f2*cep2*rep)*hh
        aksrc(ICVL,ike(1))=aksrc(ICVL,ike(1))+prod+gtx-rep
!
        diagk(ICVL,ike(2))=diagk(ICVL,ike(2))+hh*wep*f2
        diagk(ICVL,ike(1))=diagk(ICVL,ike(1))+hh
        enddo
  500   enddo
!----------------------------------------
! Near wall treatment for Low Re model
!----------------------------------------
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          IBFL=LFUTAU(ICVL)
          rans_e=aks(ICVL,ike(2))
          yp=DISALL(ICVL)
          rans_k=aks(ICVL,ike(1)) 
          rnu=rmu(ICVL)/rho(ICVL) 
          Rey=yp*sqrt(rans_k)/rnu  !Rey 
          fmu=(1.d0-exp(-Rey/Amu))  
          dum2=2.d0*rnu*rans_k/yp**2
          IF(fmu<0.99D0) THEN 
!          IF(Rey<200.d0) THEN 
            aksrc(ICVL,ike(2))=(dum2-aks(ICVL,ike(2)))*GREAT 
            diagk(ICVL,ike(2))=GREAT/rho(ICVL)
          ENDIF
          enddo
        enddo
      elseif(icaltb==SDES) then
        do 700 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        dume1=1.d0/akdes**2
        do ICVL=ICVS,ICVE
        dens=rho(ICVL)
        rans_nu=aks(ICVL,ides)
        xxx=rans_nu/rmu(ICVL)*dens
        fv1des=xxx**3/(xxx**3+cv1**3)
!
        cv2=5.d0
        fv2des=1.d0-(xxx/(1.d0+xxx*fv1des)) !old
!        fv2des=1.d0/(1.d0+xxx/cv2)**3      !modified
        fv3des=1.d0                                 !old
!        fv3des=(1.d0+xxx*fv1des)*(1.d0-fv2des)/xxx !modified
!
        Scal=CVVOLM(ICVL)**(1.d0/3.d0)
        Ddis=min(DISALL(ICVL),Scal*cDES) 
        Stud=fv3des*ROOT_S(ICVL)
     &      +rans_nu*fv2des*dume1/Ddis**2
        rdes=rans_nu/Stud/Ddis**2*dume1
        gdes=rdes+cw2*(rdes**6-rdes)
        fwdes=gdes*((1.d0+cw3**6)/(gdes**6+cw3**6))**(1.d0/6.d0)

        ft2=0.d0                   !old
!        ft2=ct3*exp(-ct4*xxx**2)  !modified
        Pprod=cb1*dens*Stud*rans_nu*(1.d0-ft2)
        Yprod=dens*(rans_nu**2/Ddis**2)*(cw1*fwdes-cb1*dume1*ft2)
        Yprod1=dens*(rans_nu/Ddis**2)*cw1*fwdes
!        
        aksrc(ICVL,ides)=aksrc(ICVL,ides)+Pprod-Yprod
        diagk(ICVL,ides)=diagk(ICVL,ides)+Yprod1
!
        enddo
  700   enddo
      elseif(icaltb==KLES) then
        do 800 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)

        do ICVL=ICVS,ICVE
        dens=rho(ICVL)
        rnu=rmu(ICVL)/(dens+SML)
        IBFL=LFUTAU(ICVL)
        yplus=utau(IBFL)*DISALL(ICVL)/(rnu+SML)
        fmu=(ONE-EXP(-yplus*R26))
        del=min(CVVOLM(ICVL)**r1pn,DISALL(ICVL))

        Scal=del !CVVOLM(ICVL)**(1.d0/3.d0)



        dum1=Ce*aks(ICVL,ikles)**(0.5d0)/Scal*dens
        dspd=Ce*aks(ICVL,ikles)**(3.d0/2.d0)/Scal*dens

!        Pprod=dspf(ICVL)*rmut(ICVL)
        Pprod=dspf(ICVL)*Ck*Scal*aks(ICVL,ikles)**(0.5d0)*fmu

        aksrc(ICVL,ikles)=aksrc(ICVL,ikles)+Pprod-dspd
        diagk(ICVL,ikles)=diagk(ICVL,ikles)+dum1

        enddo
  800   enddo
        
      endif
!
      
!
      return
!
      end subroutine keps_src
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine cavit_src(iter,deltt,time,imode,ns,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,rho2,tmp,aks,diagk,src,prs,
     &  LCYCOLD,wifsld)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_constant
      use module_scalar,only   : ical_cavi,icavi,Pres_c,
     &                           P_str_v,Tref,T0_ref,KL,Rv,ivold,
     &                           iterCVAI,CeA_PP,CeA_MM
      use module_initial ,only : rho0,rho02
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kxnone
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: LVEDGE (2, MXCVFAC),imode,ns,iter
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal ( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4, MXCVFAC)
      real*8 ,intent(in)    :: wiface (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (   MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3, MXCVFAC)
      real*8 ,intent(in)    :: rho    (   MXALLCV)
      real*8 ,intent(in)    :: rho2   (   MXALLCVC)
      real*8 ,intent(in)    :: tmp    (   MXALLCV)
      real*8 ,intent(inout) :: aks    (   MXALLCVR,MXrans)
      real*8 ,intent(inout) :: src    (   MXALLCVR,ns)
      real*8 ,intent(inout) :: diagk(     MXALLCVR,ns)
      integer,intent(in)    :: LCYCOLD(   MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (   MXSSFBC_SLD)
      real*8 ,intent(in)    :: prs(MXALLCV)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,nb,kd
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
     &           ICV,IBFL,IBFS,IBFE
      real*8  :: dum1,dum2,dum3,dum4,rhol_rhov,Pv,ppp,TTT,rhodum,
     &           void
      real*8,external :: Pvap
!
      if(iter<iterCVAI) return
      if(imode==0.and..true.) then
        do 600 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO 630 ICVL=ICVS,ICVE
        ppp=prs(ICVL)                !max(0.d0,prs(ICVL))
        TTT=tmp(ICVL)
        void=aks(ICVL,ivold)
        Pv=Pvap(TTT)
        dum4=0.d0
        if(ppp<Pv) then   
          dum3=Pv*KL*(tmp(ICVL)+T0_ref)
!!!!!!!          rhol_rhov=(ppp+Pres_c)*Rv*TTT/dum3
          rhol_rhov=rho0(IIMAT)/rho2(ICVL)
          dum2=void*(1.d0-void)*void*(1.d0-void)
          dum1=CeA_PP*dum2*rhol_rhov*(Pv-ppp)
     &      /sqrt(2.d0*3.14d0*Rv*TTT)
!CeA_PP=Ce*Ca
          dum4=dum1
        else
          dum3=rho0(IIMAT)/
     &     (aks(ICVL,icavi)*rho0(IIMAT)
     &    +(1.d0-aks(ICVL,icavi))*rho2(ICVL))

          dum2=void*(1.d0-void)*void*(1.d0-void)

          dum1=CeA_MM*dum2*(Pv-ppp)/sqrt(2.d0*3.14d0*Rv*TTT)
          dum4=dum1
          
          dum2=CeA_MM*void*(1.d0-void)**2*(Pv-ppp)  !Pv<ppp
     &      /sqrt(2.d0*3.14d0*Rv*TTT)*dum2

          diagk(ICVL,icavi)=diagk(ICVL,icavi)-dum2

        endif
        src(ICVL,icavi)=src(ICVL,icavi)+dum4
 630    enddo
 600    enddo
      elseif(imode==1) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle 
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICVL=ICVS,ICVE
        dum1=0.d0
        src(ICVL,icavi)=src(ICVL,icavi)+dum1
        enddo
        enddo
      endif
!------------------
! --- BC 
!------------------
!
      if(.false.) then
        do nb=1,nbcnd
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        IIMAT=MAT_BCIDX(nb,1)
        kd=kdbcnd(0,nb)
        if(kd==kxnone) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVL=LVEDGE(1,ICFL)
          src(ICVL,icavi)=0.d0   !src(ICVL,icavi)+dum1
          enddo
        endif
        enddo
      endif
!
!
      end subroutine cavit_src
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine C_dash_src(
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  DISALL,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  grdc,dtau,rmut,rmu,rho,vel,tmp,aks,yys,aksrc,
     &  LCYCOLD,wifsld,utau,LFUTAU,
     &  diagk,vctr)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_gravity,only   : ggg
      use module_rans,only      : cep1,cep2,cep3,cep4,cep5,eta0,beta,
     &                            cmu
      use module_rans,only      : sgmaks
      use module_turbparm,only  : prturb,akappa,awfnc,yplsm,
     &                            Rey_s,Amu,dRey
      use module_scalar,only    : ike,ike_c
      use module_model,only     : icaltb,ke,ke_low,RNG,CHEN,SDES,
     &                            ical_EWT,KE2S
      use module_material,only  : ical_sld
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LVEDGE (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF(   MXSSFBC)
      integer,intent(in)    :: LCYCSF (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV(    MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal ( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4, MXCVFAC)
      real*8 ,intent(in)    :: wiface (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (   MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3, MXCVFAC)
      real*8 ,intent(in)    :: DISALL (   MXCV)
      real*8 ,intent(inout) :: grdc   (   MXALLCV,3,3)
      real*8 ,intent(inout) :: dtau   (   MXALLCV,3)
      real*8 ,intent(in)    :: rmut   (   MXALLCV)
      real*8 ,intent(in)    :: rmu       (   MXALLCV)
      real*8 ,intent(in)    :: rho    (   MXALLCV)
      real*8 ,intent(in)    :: vel    (   MXALLCV,3)
      real*8 ,intent(in)    :: tmp    (   MXALLCV)
      real*8 ,intent(in)    :: aks    (   MXALLCVR,MXrans)
      real*8 ,intent(inout) :: aksrc  (   MXALLCVR,MXrans)
      real*8 ,intent(inout) :: diagk  (   MXALLCVR,MXrans)


      real*8 ,intent(in)    :: yys    (   MXALLCV,MXCOMP,2)
      integer,intent(in)    :: LCYCOLD(MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      REAL*8 ,INTENT(IN)    :: utau (0:MXSSFBC)
      integer,intent(in)    :: LFUTAU (    MXCV )
!
! --- [local entities]
!
      real*8,parameter :: r2p3=2.d0/3.d0,SML1=1.d-10   !
      real*8  :: dum1,dum2,dum3

      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
!
!-< 1. Calculate >-
!
      if(MXCOMP/=2) then
        call FFRABORT
     &  (1,'ERR: [atmospheric diffusion model] NOT support [ncomp/=2]')
      endif
!
      call grad_cell(1,19,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,yys(:,2,1),grdc(:,:,1))
!
      do 600 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      DO 630 ICVL=ICVS,ICVE
      dum1=(grdc(ICVL,1,1)**2
     &     +grdc(ICVL,2,1)**2
     &     +grdc(ICVL,3,1)**2)
     &     *rmut(ICVL)/sgmaks(ike_c)
      dum2=aks(ICVL,ike(2))/(0.8d0*aks(ICVL,ike(1)))
      dum3=-dum2*aks(ICVL,ike_c)  ! Rodi model

      aksrc(ICVL,ike_c)=aksrc(ICVL,ike_c)+dum1+dum3
      diagk(ICVL,ike_c)=diagk(ICVL,ike_c)+dum2
 630  enddo
 600  enddo

      return
!
      end subroutine C_dash_src


