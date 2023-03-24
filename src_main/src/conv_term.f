!
!     subroutine conv_term
!     subroutine conv_term_vel
!     subroutine conv_term_vof
!     subroutine conv_term_vel_e2p
!     subroutine conv_term_e2p
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine conv_term
     & (icnv,lmtr,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  qfv,vel,
     &  grdc,rva,q,dqdt,
     &  IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     qfv(1:2,:)  <=  rvd(1:2,:)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine if for Calculating convective term
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_boundary,only: nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                          kdolet,kdilet,kdsld
      use module_metrix,only  : phi =>W1K8
      use module_metrix,only  : cnvv=>W1K9
      use module_metrix,only  : cnvq=>W1K7
      USE module_dimension
      use module_material,only : nflud,venkat,rnck,MUSCL_k,slope,
     &                           up1st,up2nd,up3rd,cnt2nd,bldf,
     &                           usi2nd,engcsv,c3bld,mscl
      use module_Euler2ph,only : ieul2ph 
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_model,only    : ical_vect
!
      implicit none
!
! 1. Calculate convection term
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt
      integer,intent(in)    :: icnv(nflud),lmtr(nflud),IMODE
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3,MXCVFAC)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(  MXALLCV)
      real*8 ,intent(inout) :: qfv   (  MXCVFAC,2)
      real*8 ,intent(in)    :: rva   (  MXCVFAC)
      real*8 ,intent(in)    :: q     (  MXALLCV)
      real*8 ,intent(inout) :: dqdt  (  MXALLCV)
      real*8 ,intent(inout) :: grdc  (  MXALLCV,3)
      real*8 ,intent(in)    :: vel   (  MXALLCV,3)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      real*8  :: dq  (2)
      real*8 ,parameter :: ak=0.5d0
      real*8  :: dqp,dqm,dqpp,dqmm,dqpm,wi1,wi2
      real*8  :: dx,dy,dz,dab,da,db
      real*8  :: epsw,sgn,dum1,qqq,alhpa,dum2,dum3
      integer :: i,j,k,l,m,n
      integer :: nb,kd,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC
      integer :: IBFS,IBFE,IBFL,IE,ICFP,ICVP
      integer :: ierr=0
      logical :: icnTVD
      real*8 ,save :: vol_cov,dxx,dyx,dzx,dl
!
! --- START 
!
      
      icnTVD=.false.
!
      do 1000 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 1000
      IF(icnv(IMAT)==up2nd.OR.
     &   icnv(IMAT)==up3rd.OR.
     &   icnv(IMAT)==mscl.OR.
     &   icnv(IMAT)==usi2nd
     &  ) THEN
!        if(ical_vect) then
!          call FFRABORT(1,'ERR: VectorVer ONLY support c2d,1st')
!        else
          icnTVD=.true.
!        endif
      endif
      if(icnv(IMAT)==c3bld) then
         call FFRABORT(1,'ERR: [c3d] NOT support  [ys,t,k]')
      endif
 1000 continue
      if(NPE.gt.1) call hpclor(icnTVD)
      if(icnTVD) then
        call grad_cell(1,6,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,q,grdc)
      endif
!
!-< 1. 1st order upwind >-
!
      do 100 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 100
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 100
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(icnv(IMAT).eq.up1st) then
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        qfv(ICFL,1)=q(ICVLA)
        qfv(ICFL,2)=q(ICVLB)
        enddo
      endif
  100 enddo
!-----------------------------------------------------------------
!-< 1. 2nd order centeral >-: bldf=0.D0: 1stUD; bldf=1.D0: 2ndCD;
!-----------------------------------------------------------------
      do 150 IIMAT=1,NMAT  
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(icnv(IMAT).eq.cnt2nd) then
        dum1=bldf(IMAT)
        dum2=1.D0-dum1
        if(IMODE.eq.-1) then
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=q(ICVLA)*wi1+q(ICVLB)*wi2
          qfv(ICFL,1)=q(ICVLA)
          qfv(ICFL,2)=q(ICVLB)
          enddo
        else
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=q(ICVLA)*wi1+q(ICVLB)*wi2
          qfv(ICFL,1)=dum1*qqq+dum2*q(ICVLA)
          qfv(ICFL,2)=dum1*qqq+dum2*q(ICVLB)
          enddo
        endif
      endif
  150 enddo
!
!-< 3. Higher order with Venkatakrishnan limiter >-
!
!
!--< 3.1 calculate gradient in cell >--
!
!
!--< 3.2 calculate value on face with limited gradient >--
!        qmax(:)<=cnvq(:)
!        qmin(:)<=cnvv(:)
!
      phi(:)=1.d0
      do 250 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) goto 250
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 250
      if(icnv(IMAT)==up2nd.or.icnv(IMAT)==up3rd.or.
     &   icnv(IMAT)==mscl) then
        if(lmtr(IMAT).eq.venkat.or.lmtr(IMAT).eq.slope) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          cnvq(ICVS:ICVE)=q(ICVS:ICVE)
          cnvv(ICVS:ICVE)=q(ICVS:ICVE)
          cnvq(IDCS:IDCE)=q(IDCS:IDCE)
          cnvv(IDCS:IDCE)=q(IDCS:IDCE)
          do 255 ICFL=ICFS,ICFE
          do 260 I=1,2
          ICVL=LVEDGE(I,ICFL)
          dq(I)=grdc(ICVL,1)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &         +grdc(ICVL,2)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &         +grdc(ICVL,3)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
 260      continue
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          do 257 I=1,2
          ICVL=LVEDGE(I,ICFL)
          if(abs(dq(i)).gt.SML) then
            cnvq(ICVL)=max(cnvq(ICVL),q(ICVLA),q(ICVLB))
            cnvv(ICVL)=min(cnvv(ICVL),q(ICVLA),q(ICVLB))
            dqm=dq(i)
            epsw=rnck*rnck*rnck*CVVOLM(ICVL)
            sgn=sign(1.d0,dqm)
            dqp=max(sgn*(cnvq(ICVL)-q(ICVL)),
     &              sgn*(cnvv(ICVL)-q(ICVL)))
            dqpp=dqp*dqp
            dqmm=dqm*dqm
            dqpm=abs(dqp*dqm)
            if(lmtr(IMAT).eq.venkat) then
              phi(ICVL)=min(phi(ICVL),
     &         (dqpp+dqpm+dqpm+epsw)/(dqpp+dqmm+dqmm+dqpm+epsw+SML))
            elseif(lmtr(IMAT).eq.slope) then
              epsw=1.d-3*rnck
              phi(ICVL)=min(phi(ICVL),(dqpm+dqpm+epsw)/(dqpp+dqmm+epsw))
            endif
	  else
 	    phi(ICVL)=1.0d0
          endif
 257      continue
 255      continue
        endif
      elseif(icnv(IMAT)==usi2nd) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          phi(ICVS:ICVE)=1.d0
      endif
 250  continue
!
      do 200 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 200
      if(icnv(IMAT).eq.up2nd) then
        cnvq(ICVS:ICVE)=q(ICVS:ICVE)
        cnvv(ICVS:ICVE)=q(ICVS:ICVE)
        do 210 ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        do 201 I=1,2
        ICVL=LVEDGE(I,ICFL)
        cnvq(ICVL)=max(cnvq(ICVL),q(ICVLA),q(ICVLB))
        cnvv(ICVL)=min(cnvv(ICVL),q(ICVLA),q(ICVLB))
        dq(I)=grdc(ICVL,1)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &       +grdc(ICVL,2)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &       +grdc(ICVL,3)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
  201   enddo
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        qfv(ICFL,1)=q(ICVLA)+phi(ICVLA)*dq(1)
        qfv(ICFL,2)=q(ICVLB)+phi(ICVLB)*dq(2)
 210    continue
      endif
 200  continue
!
!-< 4. 3rd TVD Scheme >
!
      do 420 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) goto 420
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 420
      if(icnv(IMAT)==up3rd.or.icnv(IMAT)==mscl) then
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dab=q(ICVLB)-q(ICVLA)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
        dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
        dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
        da=grdc(ICVLA,1)*dx+grdc(ICVLA,2)*dy+grdc(ICVLA,3)*dz
        db=grdc(ICVLB,1)*dx+grdc(ICVLB,2)*dy+grdc(ICVLB,3)*dz
        qfv(ICFL,1)=q(ICVLA)+0.5d0*phi(ICVLA)*((1.d0
     &      -MUSCL_k(IMAT)*phi(ICVLA))*da+MUSCL_k(IMAT)*phi(ICVLA)*dab)
        qfv(ICFL,2)=q(ICVLB)-0.5d0*phi(ICVLB)*((1.d0
     &      -MUSCL_k(IMAT)*phi(ICVLB))*db+MUSCL_k(IMAT)*phi(ICVLB)*dab)
        enddo
      endif
  420 continue
!----------------------------------------------------------
!-< 4. Upstream Shift Interpolation (usi) 2nd Scheme >
!----------------------------------------------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      if(icnv(IMAT)==usi2nd) then
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=SFAREA(1,ICFL)
        dum2=SFAREA(2,ICFL)
        dum3=SFAREA(3,ICFL)
        dab=vel(ICVLA,1)*dum1+vel(ICVLA,2)*dum2+vel(ICVLA,3)*dum3
        dx=dum1*dab*deltt
        dy=dum2*dab*deltt
        dz=dum3*dab*deltt
        dq(1)=
     &        (grdc(ICVLA,1)*dx
     &        +grdc(ICVLA,2)*dy
     &        +grdc(ICVLA,3)*dz)
        dab=vel(ICVLB,1)*dum1+vel(ICVLB,2)*dum2+vel(ICVLB,3)*dum3
        dx=dum1*dab*deltt
        dy=dum2*dab*deltt
        dz=dum3*dab*deltt
        dq(2)=
     &        (grdc(ICVLB,1)*dx
     &        +grdc(ICVLB,2)*dy
     &        +grdc(ICVLB,3)*dz)
        qfv(ICFL,1)=q(ICVLA)+dq(1)
        qfv(ICFL,2)=q(ICVLB)-dq(2)
        enddo
      endif
      enddo
!
!-< 5. Calculate convection term >-
!
      do 300 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      kd=kdbcnd(0,nb)
      if(icnv(IMAT).eq.cnt2nd) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        qfv(ICFL,1)=q(ICV)
        qfv(ICFL,2)=q(IDC)
        enddo
      elseif(icnv(IMAT).ne.up1st) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        qfv(ICFL,2)=q(IDC)
        enddo
      endif
        if(kd.eq.kdsld) then  !sldz
          cycle
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)

          qfv(ICFL,2)=q(ICVP)
          qfv(ICFL,1)=q(ICV)
          enddo
        endif
!      if(kd==kdolet) then !15763
!        do IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        ICV=LVEDGE(1,ICFL)
!        qfv(ICFL,2)=q(ICV)
!        enddo
!      endif
  300 continue
!
      if(ical_vect) then
        cnvq=0.d0 
        cnvv=0.d0 
        do IE=1,MAXIE 
        DO ICVL=ICVS_V,ICVE_V 
        IF(ABS(vctr(ICVL,IE))>0) then 
          dum1=max(0.d0,rva(ABS(vctr(ICVL,IE))))
     &        *qfv(ABS(vctr(ICVL,IE)),1)
     &        +min(0.d0,rva(ABS(vctr(ICVL,IE))))
     &        *qfv(ABS(vctr(ICVL,IE)),2)
          cnvq(ICVL)=cnvq(ICVL)
     &        -sign(1,vctr(ICVL,IE))*dum1
          cnvv(ICVL)=cnvv(ICVL)
     &         -sign(1,vctr(ICVL,IE))*rva(ABS(vctr(ICVL,IE)))
        endif
        enddo
        enddo
!
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dqdt(ICVL)=dqdt(ICVL)+(q(ICVL)*cnvv(ICVL)-cnvq(ICVL))
        enddo 
        enddo
      else
       if(IMODE==29.or.IMODE==30) then
        cnvq=0.d0
        cnvv=0.d0
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &      +min(0.d0,rva(ICFL))*qfv(ICFL,2)
        cnvq(ICVLA)=cnvq(ICVLA)+dum1
        cnvq(ICVLB)=cnvq(ICVLB)-dum1
        cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)
        cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)
        enddo
        enddo
!
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dqdt(ICVL)=dqdt(ICVL)-cnvv(ICVL)
        enddo
        enddo
       elseif(IMODE==72) then
        cnvq=0.d0
        cnvv=0.d0
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &      +min(0.d0,rva(ICFL))*qfv(ICFL,2)
        cnvq(ICVLA)=cnvq(ICVLA)+dum1
        cnvq(ICVLB)=cnvq(ICVLB)-dum1
        cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)
        cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)
        enddo
        enddo
!
        do IIMAT=1,NMAT   !IMODE==72
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dqdt(ICVL)=dqdt(ICVL)+(q(ICVL)*cnvv(ICVL)-cnvq(ICVL))
!+(q(ICVL)*cnvv(ICVL)-cnvq(ICVL))
!-cnvq(ICVL)
!+(q(ICVL)*cnvv(ICVL)-cnvq(ICVL))
!-q(ICVL)*cnvv(ICVL)   
        enddo
        enddo

       else
        cnvq=0.d0 
        cnvv=0.d0 
        do 310 IIMAT=1,NMAT    !ICF=1,NCVFAC
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &      +min(0.d0,rva(ICFL))*qfv(ICFL,2)
        cnvq(ICVLA)=cnvq(ICVLA)+dum1
        cnvq(ICVLB)=cnvq(ICVLB)-dum1
        cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)
        cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)
        enddo
  310   enddo

!
        do 320 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        dqdt(ICVL)=dqdt(ICVL)+(q(ICVL)*cnvv(ICVL)-cnvq(ICVL))
 330    enddo 
 320    enddo
       endif
      endif
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,dqdt)
      ENDIF
!
      return
      end subroutine conv_term
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine conv_term_vel
     & (iphs,icnv,lmtr,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  qfv,dvdd,diag,
     &  LCYCOLD,wifsld,
     &  grdc,rva,q,dqdt
     &  )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------
!     qfv(1:2,:)  <=  rvd(1:2,:)
!-----------------------------------------
! --- This subroutine if for Calculating convective term 
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
!      use module_metrix,only : q2  =>W2K1
      use module_metrix,only : phi =>W1K8
      use module_metrix,only : cnvv=>W1K9
      use module_metrix,only : cnvq=>W1K7 !cnvq=>W1K10
      use module_scalar,only : iaph
      use module_model,only  : ical_vect,nthrds
      use module_vector,only : ICVS_V,ICVE_V,
     &                         ICFS_V,ICFE_V,
     &                         ICVSIN_V,ICVEIN_V,
     &                         IDCS_V,IDCE_V,index_c,index_f
      use module_material,only : ical_sld,porosty,idarcy2,porous
      use module_metrix,only   : msk
!
! --- [module arguments]
!
      use module_material,only   : nflud,venkat,rnck,MUSCL_k,slope,
     &                             up1st,up2nd,up3rd,cnt2nd,bldf,
     &                             usi2nd,engcsv,c3bld,mscl
      use module_boundary,only   : nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                             kdolet,kdsld,kdcvd,kdilet,kdprdc,
     &                             kdbuff,kdovst
!      use module_Euler2ph,only   : ieul2ph
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt
      integer,intent(inout) :: icnv(nflud)
      integer,intent(in)    :: lmtr(nflud),iphs
      integer,intent(in)    :: LVEDGE(2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF(   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3, MXCVFAC)
      real*8 ,intent(in)    :: wiface(   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(   MXALLCV)
      real*8 ,intent(inout) :: qfv   (   MXCVFAC,2)
      real*8 ,intent(inOUT) :: rva   (   MXCVFAC)
      real*8 ,intent(inout) :: q     (   MXALLCV,3)
      real*8 ,intent(inout) :: diag  (   MXALLCV)
      real*8 ,intent(inout) :: dqdt  (   MXALLCV,3)
      real*8 ,intent(inout) :: grdc  (   MXALLCV,3,3)
      real*8 ,intent(inout) :: dvdd  (   MXCV,3)
      integer,intent(in)    :: LCYCOLD  (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld   (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      real*8  :: dq(2),rrva,rrvb
      real*8 ,parameter :: ak=0.5d0   !(0.3333,0.5,1.0)
      real*8  :: dqp,dqm,dqpp,dqmm,dqpm,wi1,wi2
      real*8  :: dx,dy,dz,dab,da,db
      real*8  :: epsw,sgn,dum1,qqq,alhpa,dum2,dum3,dum4
      integer :: i,j,k,l,m,n,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP,ICVBO
      integer :: IBFS,IBFE,IBFL,ICFP,ICFO
      integer :: IMODE,ICOM,IMD,IE
      integer :: nb,kd,kdt,kdy,kdk,kdp
      integer :: ierr=0
      logical :: icnTVD
      integer,save :: MAT_SLD
!
!
!      IF(ical_sld==1.or.ical_sld==2.or.ical_sld==4)then
!        MAT_SLD=icnv(1)
!        do IIMAT=1,NMAT
!          IMAT=MAT_NO(IIMAT)
!          if(IMAT>0) then
!            icnv(IMAT)=icnv(1)
!          else
!            call FFRABORT(1,'SLIDING CON. ERR')
!          endif
!        enddo
!      endif
!
! --- 
!
      icnTVD=.false.
      do 500 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IF(icnv(IMAT)==up2nd.OR.
     &   icnv(IMAT)==up3rd.OR.
     &   icnv(IMAT)==mscl.OR.
     &   icnv(IMAT)==usi2nd.OR.
     &   icnv(IMAT)==c3bld
     &  ) THEN
!        if(ical_vect) then
!          call FFRABORT(1,'ERR: VectorVer ONLY support c2d,1st')
!        else
          icnTVD=.true.
!        endif
      endif
 500  continue
      if(NPE.gt.1) call hpclor(icnTVD)
!
      if(icnTVD) then
        call grad_cell(3,7,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,q,grdc)
      ENDIF
!----------------------------------------------------
!
!
!
!----------------------------------------------------
!----------------------------------------------
! --- Scalar CPU
!----------------------------------------------
!      else
!---------------------------------
        do 1000 l=1,3
!---------------------------------
!-< 1. 1st order upwind >-
!---------------------------------
        do 100 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        if(icnv(IMAT).eq.up1st) then
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          qfv(ICFL,1)=q(ICVLA,l)
          qfv(ICFL,2)=q(ICVLB,l)
          enddo
        endif
 100    enddo
!-----------------------------------------------------------------
!-< 1. 2nd order centeral >-: bldf=0.D0: 1stUD; bldf=1.D0: 2ndCD;
!-----------------------------------------------------------------
        do 150 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        if(icnv(IMAT).eq.cnt2nd) then
          dum1=bldf(IMAT)
          dum2=1.D0-dum1
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=q(ICVLA,l)*wi1+q(ICVLB,l)*wi2
          qfv(ICFL,1)=dum1*qqq+dum2*q(ICVLA,l)
          qfv(ICFL,2)=dum1*qqq+dum2*q(ICVLB,l)
          enddo
        endif
  150   enddo
!--------------------------
! --- 
!--------------------------
        do 160 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        if(icnv(IMAT).eq.engcsv) then
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qfv(ICFL,1)=q(ICVLB,l)*wi2
          qfv(ICFL,2)=q(ICVLA,l)*wi1
          enddo
        endif
 160    enddo
!
!
!-< 3. Higher order with Venkatakrishnan limiter >-
!
! --- / calc. limiter function : phi(ICVL)/
!
        phi(:)=1.d0
        do 250 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 250
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 250
        if(icnv(IMAT)==up2nd.or.icnv(IMAT)==up3rd.or.
     &  icnv(IMAT)==mscl.or.icnv(IMAT)==c3bld) 
     &  then
          if(lmtr(IMAT).eq.venkat.or.lmtr(IMAT).eq.slope) then
            ICFS=MAT_CFIDX(IIMAT-1)+1
            ICFE=MAT_CFIDX(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            IDCS=MAT_DCIDX(IIMAT-1)+1
            IDCE=MAT_DCIDX(IIMAT)
            cnvq(ICVS:ICVE)=q(ICVS:ICVE,l)
            cnvv(ICVS:ICVE)=q(ICVS:ICVE,l)
            cnvq(IDCS:IDCE)=q(IDCS:IDCE,l)
            cnvv(IDCS:IDCE)=q(IDCS:IDCE,l)
            do 255 ICFL=ICFS,ICFE
            do 260 I=1,2
            ICVL=LVEDGE(I,ICFL)
            dq(I)=grdc(ICVL,1,l)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &           +grdc(ICVL,2,l)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &           +grdc(ICVL,3,l)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
 260        continue
            ICVLA=LVEDGE(1,ICFL)
            ICVLB=LVEDGE(2,ICFL)
            do 257 I=1,2
            ICVL=LVEDGE(I,ICFL)
            if(abs(dq(i)).gt.SML) then
              cnvq(ICVL)=max(cnvq(ICVL),q(ICVLA,l),q(ICVLB,l))
              cnvv(ICVL)=min(cnvv(ICVL),q(ICVLA,l),q(ICVLB,l))
              dqm=dq(i)
              epsw=rnck*rnck*rnck*CVVOLM(ICVL)
              sgn=sign(1.d0,dqm)
              dqp=max(sgn*(cnvq(ICVL)-q(ICVL,l)),
     &                sgn*(cnvv(ICVL)-q(ICVL,l)))
              dqpp=dqp*dqp
              dqmm=dqm*dqm
              dqpm=abs(dqp*dqm)
              if(lmtr(IMAT).eq.venkat) then
                phi(ICVL)=min(phi(ICVL),
     &           (dqpp+dqpm+dqpm+epsw)/(dqpp+dqmm+dqmm+dqpm+epsw+SML))
              elseif(lmtr(IMAT).eq.slope) then
                epsw=1.d-3*rnck
                phi(ICVL)=min(phi(ICVL)
     &                   ,(dqpm+dqpm+epsw)/(dqpp+dqmm+epsw))
              endif
	    else
 	      phi(ICVL)=1.0d0
            endif
 257        continue
 255        continue
          endif
        elseif(icnv(IMAT)==usi2nd) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            phi(ICVS:ICVE)=1.d0

        endif
 250    continue
!
!--< 3.2 calculate value on face with limited gradient >--
!
! --- 2nd TVD:
!
! --- phi(ICVL)=0.0d0=> 1st Upwind; phi(ICVL)=1: 2nd Upwind
!
        do 200 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 200
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 200
        if(icnv(IMAT).eq.up2nd) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do 210 ICFL=ICFS,ICFE
          do 201 I=1,2
          ICVL=LVEDGE(I,ICFL)
          dq(I)=grdc(ICVL,1,l)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &         +grdc(ICVL,2,l)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &         +grdc(ICVL,3,l)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
  201     enddo
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          qfv(ICFL,1)=q(ICVLA,l)+phi(ICVLA)*dq(1)
          qfv(ICFL,2)=q(ICVLB,l)+phi(ICVLB)*dq(2)
 210      enddo
        endif
 200    enddo
!
!-< 4. 3rd TVD Scheme >
!
! --- ak=1/3: 3rd Upwind; ak=-1: 2nd Upwind;
!
        do 420 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 420
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 420
        if(icnv(IMAT)==up3rd.or.icnv(IMAT)==mscl) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do 430 ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dab=q(ICVLB,l)-q(ICVLA,l)
          dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
          dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
          dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
          da=grdc(ICVLA,1,l)*dx+grdc(ICVLA,2,l)*dy+grdc(ICVLA,3,l)*dz
          db=grdc(ICVLB,1,l)*dx+grdc(ICVLB,2,l)*dy+grdc(ICVLB,3,l)*dz
          qfv(ICFL,1)=q(ICVLA,l)+0.5d0*phi(ICVLA)*((1.d0
     &       -MUSCL_k(IMAT)*phi(ICVLA))*da+MUSCL_k(IMAT)*phi(ICVLA)*dab)
          qfv(ICFL,2)=q(ICVLB,l)-0.5d0*phi(ICVLB)*((1.d0
     &       -MUSCL_k(IMAT)*phi(ICVLB))*db+MUSCL_k(IMAT)*phi(ICVLB)*dab)
 430      enddo
        elseif(icnv(IMAT)==c3bld) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          dum1=bldf(IMAT)
          dum2=1.D0-dum1
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dab=q(ICVLB,l)-q(ICVLA,l)
          dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
          dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
          dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
          da=grdc(ICVLA,1,l)*dx+grdc(ICVLA,2,l)*dy+grdc(ICVLA,3,l)*dz
          db=grdc(ICVLB,1,l)*dx+grdc(ICVLB,2,l)*dy+grdc(ICVLB,3,l)*dz
          dum3=q(ICVLA,l)+0.5d0*phi(ICVLA)*((1.d0
     &      -MUSCL_k(IMAT)*phi(ICVLA))*da+MUSCL_k(IMAT)*phi(ICVLA)*dab)
          dum4=q(ICVLB,l)-0.5d0*phi(ICVLB)*((1.d0
     &      -MUSCL_k(IMAT)*phi(ICVLB))*db+MUSCL_k(IMAT)*phi(ICVLB)*dab)
          qqq=q(ICVLA,l)*wi1+q(ICVLB,l)*wi2
          qfv(ICFL,1)=qqq*dum1+dum2*dum3
          qfv(ICFL,2)=qqq*dum1+dum2*dum4
          enddo
        endif
 420    continue
!----------------------------------------------------------
!-< 4. Upstream Shift Interpolation (usi) 2nd Scheme >
!----------------------------------------------------------
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        if(icnv(IMAT)==usi2nd) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          dum1=SFAREA(1,ICFL)
          dum2=SFAREA(2,ICFL)
          dum3=SFAREA(3,ICFL)
          dab=q(ICVLA,1)*dum1+q(ICVLA,2)*dum2+q(ICVLA,3)*dum3
          dx=dum1*dab*deltt
          dy=dum2*dab*deltt
          dz=dum3*dab*deltt
!
          dq(1)=
     &          (grdc(ICVLA,1,l)*dx
     &          +grdc(ICVLA,2,l)*dy
     &          +grdc(ICVLA,3,l)*dz)
          dab=q(ICVLB,1)*dum1+q(ICVLB,2)*dum2+q(ICVLB,3)*dum3
          dx=dum1*dab*deltt
          dy=dum2*dab*deltt
          dz=dum3*dab*deltt
          dq(2)=
     &        (grdc(ICVLB,1,l)*dx
     &        +grdc(ICVLB,2,l)*dy
     &        +grdc(ICVLB,3,l)*dz)
!
          qfv(ICFL,1)=q(ICVLA,l)+dq(1)
          qfv(ICFL,2)=q(ICVLB,l)-dq(2)
          enddo
        endif
        enddo

!
!-< 5. Calculate convection term >-
!
        do 300 nb=1,nbcnd    !IBF=1,NSSFBC
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 300
        kd=kdbcnd(0,nb)
!        if(kd==kdprdc) cycle
        if(icnv(IMAT).eq.cnt2nd) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          qfv(ICFL,1)=q(ICV,l)
          qfv(ICFL,2)=q(IDC,l)
          enddo
        elseif(icnv(IMAT)/=up1st) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          qfv(ICFL,2)=q(IDC,l)
          enddo
        endif
        if(kd.eq.kdsld) then  !sldz
          cycle
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)

          qfv(ICFL,2)=q(ICVP,l)
          qfv(ICFL,1)=q(ICV,l)
          enddo
        elseif(kd.eq.kdovst) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LCYCSF(IBFL)
          qfv(ICFL,2)=q(IDC,l)
          enddo
        endif
!
        if(kd==kdprdc.or.kd==kdbuff) then  !okabe 
          IBFS=LBC_INDEX(nb-1)+1 
          IBFE=LBC_INDEX(nb) 
          do IBFL=IBFS,IBFE 
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          qfv(ICFL,2)=q(ICVP,l)
          qfv(ICFP,2)=q(ICV,l)
          enddo
        endif 
! 
!        if(kd==kdbuff) then  !okabe 15764  
!          IBFS=LBC_INDEX(nb-1)+1
!          IBFE=LBC_INDEX(nb)
!          do IBFL=IBFS,IBFE
!          ICFL=LBC_SSF(IBFL)
!          ICV=LVEDGE(1,ICFL)
!          qfv(ICFL,2)=q(,l)
!          enddo
!        endif
  300   enddo
! 
        cnvq(:)=0.d0
        cnvv(:)=0.d0
        if(ical_vect) then
            do IE=1,MAXIE
            DO ICVL=ICVS_V,ICVE_V
            IF(ABS(vctr(ICVL,IE))>0) then
              dum1=max(0.d0,rva(ABS(vctr(ICVL,IE))))
     &            *qfv(ABS(vctr(ICVL,IE)),1)
     &            +min(0.d0,rva(ABS(vctr(ICVL,IE))))
     &            *qfv(ABS(vctr(ICVL,IE)),2)
              cnvq(ICVL)=cnvq(ICVL)
     &                  -sign(1,vctr(ICVL,IE))
     &                  *dum1
              cnvv(ICVL)=cnvv(ICVL)
     &                  -sign(1,vctr(ICVL,IE))
     &                  *rva(ABS(vctr(ICVL,IE)))
            endif
            enddo
            enddo
        else
        do 310 IIMAT=1,NMAT    !ICF=1,NCVFAC 
        if(.not.mat_cal(IIMAT)) cycle 
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        if(icnv(IMAT).eq.engcsv) then
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          
          dum1=rva(ICFL)*qfv(ICFL,1)
          dum2=rva(ICFL)*qfv(ICFL,2)
          cnvq(ICVLA)=cnvq(ICVLA)+dum1
          cnvq(ICVLB)=cnvq(ICVLB)-dum2
          cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)*0.d0
          cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)*0.d0
          enddo
        else
          do ICFL=ICFS,ICFE
            ICVLA=LVEDGE(1,ICFL)
            ICVLB=LVEDGE(2,ICFL)
            dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &          +min(0.d0,rva(ICFL))*qfv(ICFL,2)
            cnvq(ICVLA)=cnvq(ICVLA)+dum1
            cnvq(ICVLB)=cnvq(ICVLB)-dum1
            cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)
            cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)
          enddo
        endif
  310   enddo
        endif
!
! --- d(r*ui*uj)=uj*d(r*ui)+r*ui*d(uj):
!
        do 320 IIMAT=1,NMAT         !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 320
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 320
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        dqdt(ICVL,l)=dqdt(ICVL,l)
     &     +(q(ICVL,l)*cnvv(ICVL)-cnvq(ICVL))
 330    continue
 320    continue
!
 1000   continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (3,MXALLCV,NCV,dqdt(:,1:3))
       ENDIF
!
      return
      end subroutine conv_term_vel
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine conv_term_vel1
     & (iphs,icnv,lmtr,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  qfv,dvdd,diag,
     &  LCYCOLD,wifsld,
     &  grdc,rva,q,dqdt
     &  )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------
!     qfv(1:2,:)  <=  rvd(1:2,:)
!-----------------------------------------
! --- This subroutine if for Calculating convective term 
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
!      use module_metrix,only : q2  =>W2K1
      use module_metrix,only : phi =>W1K8
      use module_metrix,only : cnvv=>W1K9
      use module_metrix,only : cnvq=>W1K7 !cnvq=>W1K10
      use module_scalar,only : iaph
      use module_model,only  : ical_vect,nthrds
      use module_vector,only : ICVS_V,ICVE_V,
     &                         ICFS_V,ICFE_V,
     &                         ICVSIN_V,ICVEIN_V,
     &                         IDCS_V,IDCE_V,index_c,index_f
      use module_material,only : ical_sld,porosty,idarcy2,porous
      use module_metrix,only   : msk
!
! --- [module arguments]
!
      use module_material,only   : nflud,venkat,rnck,MUSCL_k,slope,
     &                             up1st,up2nd,up3rd,cnt2nd,bldf,
     &                             usi2nd,engcsv,c3bld,mscl
      use module_boundary,only   : nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                             kdolet,kdsld,kdcvd,kdilet,kdprdc,
     &                             kdbuff
!      use module_Euler2ph,only   : ieul2ph
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt
      integer,intent(inout) :: icnv(nflud)
      integer,intent(in)    :: lmtr(nflud),iphs
      integer,intent(in)    :: LVEDGE(2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF(   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3, MXCVFAC)
      real*8 ,intent(in)    :: wiface(   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(   MXALLCV)
      real*8 ,intent(inout) :: qfv   (   MXCVFAC,2)
      real*8 ,intent(inOUT) :: rva   (   MXCVFAC)
      real*8 ,intent(inout) :: q     (   MXALLCV,3)
      real*8 ,intent(inout) :: diag  (   MXALLCV)
      real*8 ,intent(inout) :: dqdt  (   MXALLCV,3)
      real*8 ,intent(inout) :: grdc  (   MXALLCV,3,3)
      real*8 ,intent(inout) :: dvdd  (   MXCV,3)
      integer,intent(in)    :: LCYCOLD  (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld   (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      real*8  :: dq(2),rrva,rrvb
      real*8 ,parameter :: ak=0.5d0   !(0.3333,0.5,1.0)
      real*8  :: dqp,dqm,dqpp,dqmm,dqpm,wi1,wi2
      real*8  :: dx,dy,dz,dab,da,db
      real*8  :: epsw,sgn,dum1,qqq,alhpa,dum2,dum3,dum4
      integer :: i,j,k,l,m,n,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP,ICVBO
      integer :: IBFS,IBFE,IBFL,ICFP,ICFO
      integer :: IMODE,ICOM,IMD,IE
      integer :: nb,kd,kdt,kdy,kdk,kdp
      integer :: ierr=0
      logical :: icnTVD
      integer,save :: MAT_SLD
!
!
!      IF(ical_sld==1.or.ical_sld==2.or.ical_sld==4)then
!        MAT_SLD=icnv(1)
!        do IIMAT=1,NMAT
!          IMAT=MAT_NO(IIMAT)
!          if(IMAT>0) then
!            icnv(IMAT)=icnv(1)
!          else
!            call FFRABORT(1,'SLIDING CON. ERR')
!          endif
!        enddo
!      endif
!
! --- 
!
!---------------------------------
        do 1000 l=1,3
!-----------------------------------------------------------------
!-< 1. 2nd order centeral >-: bldf=0.D0: 1stUD; bldf=1.D0: 2ndCD;
!-----------------------------------------------------------------
        do 150 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        if(icnv(IMAT).eq.cnt2nd) then
          dum1=bldf(IMAT)
          dum2=1.D0-dum1
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=q(ICVLA,l)*wi1+q(ICVLB,l)*wi2
          qfv(ICFL,1)=dum1*qqq+dum2*q(ICVLA,l)
          qfv(ICFL,2)=dum1*qqq+dum2*q(ICVLB,l)
          enddo
        endif
  150   enddo
!--------------------------
! --- 
!--------------------------
!
!-< 5. Calculate convection term >-
!
        do 300 nb=1,nbcnd    !IBF=1,NSSFBC
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 300
        kd=kdbcnd(0,nb)
!        if(kd==kdprdc) cycle
        if(icnv(IMAT).eq.cnt2nd) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          qfv(ICFL,1)=q(ICV,l)
          qfv(ICFL,2)=q(IDC,l)
          enddo
        elseif(icnv(IMAT)/=up1st) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          qfv(ICFL,2)=q(IDC,l)
          enddo
        endif
        if(kd.eq.kdsld) then

          IIMAT=MAT_BCIDX(nb,1)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          qfv(ICFL,2)=q(IDC,l)
          qfv(ICFL,1)=q(ICV,l)
          enddo
        endif
!
  300   enddo
! 
        cnvq(:)=0.d0
        cnvv(:)=0.d0
        do IE=1,MAXIE
          DO ICVL=ICVS_V,ICVE_V
            IF(ABS(vctr(ICVL,IE))>0) then
              dum1=max(0.d0,rva(ABS(vctr(ICVL,IE))))
     &            *qfv(vctr(ICVL,IE),1)
     &            +min(0.d0,rva(ABS(vctr(ICVL,IE))))
     &            *qfv(vctr(ICVL,IE),2)
              cnvq(ICVL)=cnvq(ICVL)
     &                  -sign(1,vctr(ICVL,IE))
     &                  *dum1
              cnvv(ICVL)=cnvv(ICVL)
     &                  -sign(1,vctr(ICVL,IE))
     &                  *rva(ABS(vctr(ICVL,IE)))
            endif
          enddo
        enddo
!
! --- d(r*ui*uj)=uj*d(r*ui)+r*ui*d(uj):
!
        do 320 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) goto 320
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 320
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        dqdt(ICVL,l)=dqdt(ICVL,l)
     &     +(q(ICVL,l)*cnvv(ICVL)-cnvq(ICVL))
 330    continue
 320    continue
!
 1000   continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (3,MXALLCV,NCV,dqdt(:,1:3))
       ENDIF
!
      return
      end subroutine conv_term_vel1

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine conv_term_vel_e2p
     & (iphs,icnv,lmtr,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  qfv,dvdd,diag,
     &  LCYCOLD,wifsld,
     &  grdc,rva,q,dqdt,aks
     &  )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------
!     qfv(1:2,:)  <=  rvd(1:2,:)
!-----------------------------------------
! --- This subroutine if for Calculating convective term
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
!      use module_metrix,only : q2  =>W2K1
      use module_metrix,only : phi =>W1K8
      use module_metrix,only : cnvv=>W1K9
      use module_metrix,only : cnvq=>W1K7 !cnvq=>W1K10
      use module_scalar,only : iaph
      use module_model,only  : ical_vect,nthrds
      use module_vector,only : ICVS_V,ICVE_V,
     &                         ICFS_V,ICFE_V,
     &                         ICVSIN_V,ICVEIN_V,
     &                         IDCS_V,IDCE_V,index_c,index_f
      use module_material,only : ical_sld,porosty,idarcy2,porous,MUSCL_k
!
! --- [module arguments]
!
      use module_material,only   : nflud,venkat,rnck,slope,
     &                             up1st,up2nd,up3rd,cnt2nd,bldf,
     &                             usi2nd,engcsv,mscl
      use module_boundary,only   : nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                             kdolet,kdsld,kdcvd,kdilet,kdprdc
      use module_Euler2ph,only   : ieul2ph
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt
      integer,intent(inout) :: icnv(nflud)
      integer,intent(in)    :: lmtr(nflud),iphs
      integer,intent(in)    :: LVEDGE(2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF(   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3, MXCVFAC)
      real*8 ,intent(in)    :: wiface(   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(   MXALLCV)
      real*8 ,intent(inout) :: qfv   (   MXCVFAC,2)
      real*8 ,intent(inOUT) :: rva   (   MXCVFAC)
      real*8 ,intent(inout) :: q     (   MXALLCV,3)
      real*8 ,intent(inout) :: diag  (   MXALLCV)
      real*8 ,intent(inout) :: dqdt  (   MXALLCV,3)
      real*8 ,intent(inout) :: grdc  (   MXALLCV,3,3)
      real*8 ,intent(inout) :: dvdd  (   MXCV,3)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(in)    :: aks   (MXALLCVR,mxrans)
!
! --- [local entities]
!
      real*8  :: dq(2),rrva,rrvb
      real*8 ,parameter :: ak=0.3333d0   !(0.3333,0.5,1.0)
      real*8  :: dqp,dqm,dqpp,dqmm,dqpm,wi1,wi2
      real*8  :: dx,dy,dz,dab,da,db
      real*8  :: epsw,sgn,dum1,qqq,alhpa,dum2,dum3
      integer :: i,j,k,l,m,n,myid
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP,ICVBO
      integer :: IBFS,IBFE,IBFL,ICFP,ICFO
      integer :: IMODE,ICOM,IMD,IE
      integer :: nb,kd,kdt,kdy,kdk,kdp
      integer :: ierr=0
      logical :: icnTVD
      integer,save :: sldcnv=0,MAT_SLD
!
!
      IF(ical_sld==1.or.ical_sld==2.or.ical_sld==4)then
        MAT_SLD=icnv(1)
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT>0) then
            icnv(IMAT)=icnv(1)
          else
            call FFRABORT(1,'SLIDING CON. ERR')
          endif
        enddo
        sldcnv=1
      endif
!
! --- 
!
      icnTVD=.false.
      do 500 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IF(icnv(IMAT)==up2nd.OR.
     &   icnv(IMAT)==up3rd.OR.
     &   icnv(IMAT)==mscl.OR.
     &   icnv(IMAT)==usi2nd
     &  ) THEN
        icnTVD=.true.
      endif
 500  continue
      if(NPE.gt.1) call hpclor(icnTVD)
!      
      if(icnTVD) then
        call grad_cell(3,7,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,q,grdc)
      ENDIF
!----------------------------------------------------
!
!
!
!----------------------------------------------------
      if(ical_vect) then
        call FFRABORT(1,'ERR: Vector NOT support E2P')
      else
        do 1000 l=1,3
!
        do 100 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        if(icnv(IMAT).eq.up1st) then
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          qfv(ICFL,1)=q(ICVLA,l)
          qfv(ICFL,2)=q(ICVLB,l)
          enddo
        endif
 100    enddo
!
!-< 1. 2nd order centeral >-: bldf=0.D0: 1stUD; bldf=1.D0: 2ndCD;
!
        do 150 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        if(icnv(IMAT).eq.cnt2nd ) then
          dum1=bldf(IMAT)
          dum2=1.D0-dum1
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=q(ICVLA,l)*wi1+q(ICVLB,l)*wi2
          qfv(ICFL,1)=dum1*qqq+dum2*q(ICVLA,l)
          qfv(ICFL,2)=dum1*qqq+dum2*q(ICVLB,l)
          enddo
        endif
  150   enddo
!
!
!
        do 160 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        if(icnv(IMAT).eq.engcsv) then
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qfv(ICFL,1)=q(ICVLB,l)*wi2
          qfv(ICFL,2)=q(ICVLA,l)*wi1
          enddo
        endif
 160    enddo
!
!
!-< 3. Higher order with Venkatakrishnan limiter >-
!
! --- / calc. limiter function : phi(ICVL)/
!
        phi(:)=1.d0
        do 250 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 250
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 250
        if(icnv(IMAT)==up2nd.or.icnv(IMAT)==up3rd.or.
     &     icnv(IMAT)==mscl) then
          if(lmtr(IMAT).eq.venkat.or.lmtr(IMAT).eq.slope) then
            ICFS=MAT_CFIDX(IIMAT-1)+1
            ICFE=MAT_CFIDX(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            IDCS=MAT_DCIDX(IIMAT-1)+1
            IDCE=MAT_DCIDX(IIMAT)
            cnvq(ICVS:ICVE)=q(ICVS:ICVE,l)
            cnvv(ICVS:ICVE)=q(ICVS:ICVE,l)
            cnvq(IDCS:IDCE)=q(IDCS:IDCE,l)
            cnvv(IDCS:IDCE)=q(IDCS:IDCE,l)
            do 255 ICFL=ICFS,ICFE
            do 260 I=1,2
            ICVL=LVEDGE(I,ICFL)
            dq(I)=grdc(ICVL,1,l)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &           +grdc(ICVL,2,l)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &           +grdc(ICVL,3,l)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
 260        continue
            ICVLA=LVEDGE(1,ICFL)
            ICVLB=LVEDGE(2,ICFL)
            do 257 I=1,2
            ICVL=LVEDGE(I,ICFL)
            if(abs(dq(i)).gt.SML) then
              cnvq(ICVL)=max(cnvq(ICVL),q(ICVLA,l),q(ICVLB,l))
              cnvv(ICVL)=min(cnvv(ICVL),q(ICVLA,l),q(ICVLB,l))
              dqm=dq(i)
              epsw=rnck*rnck*rnck*CVVOLM(ICVL)
              sgn=sign(1.d0,dqm)
              dqp=max(sgn*(cnvq(ICVL)-q(ICVL,l)),
     &                sgn*(cnvv(ICVL)-q(ICVL,l)))
              dqpp=dqp*dqp
              dqmm=dqm*dqm
              dqpm=abs(dqp*dqm)
              if(lmtr(IMAT).eq.venkat) then
                phi(ICVL)=min(phi(ICVL),
     &           (dqpp+dqpm+dqpm+epsw)/(dqpp+dqmm+dqmm+dqpm+epsw+SML))
              elseif(lmtr(IMAT).eq.slope) then
                epsw=1.d-3*rnck
                phi(ICVL)=min(phi(ICVL)
     &                   ,(dqpm+dqpm+epsw)/(dqpp+dqmm+epsw))
              endif
	    else
 	      phi(ICVL)=1.0d0
            endif
 257        continue
 255        continue
          endif
        elseif(icnv(IMAT)==usi2nd) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            phi(ICVS:ICVE)=1.d0

        endif
 250    continue
!
!--< 3.2 calculate value on face with limited gradient >--
!
! --- 2nd TVD:
!
! --- phi(ICVL)=0.0d0=> 1st Upwind; phi(ICVL)=1: 2nd Upwind
!
        do 200 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 200
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 200
        if(icnv(IMAT).eq.up2nd) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do 210 ICFL=ICFS,ICFE
          do 201 I=1,2
          ICVL=LVEDGE(I,ICFL)
          dq(I)=grdc(ICVL,1,l)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &         +grdc(ICVL,2,l)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &         +grdc(ICVL,3,l)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
  201     enddo
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          qfv(ICFL,1)=q(ICVLA,l)+phi(ICVLA)*dq(1)
          qfv(ICFL,2)=q(ICVLB,l)+phi(ICVLB)*dq(2)
 210      enddo
        endif
 200    enddo
!
!-< 4. 3rd TVD Scheme >
!
! --- ak=1/3: 3rd Upwind; ak=-1: 2nd Upwind;
!
        do 420 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) goto 420
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 420
        if(icnv(IMAT)==up3rd.or.icnv(IMAT)==mscl) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do 430 ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dab=q(ICVLB,l)-q(ICVLA,l)
          dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
          dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
          dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
          da=grdc(ICVLA,1,l)*dx+grdc(ICVLA,2,l)*dy+grdc(ICVLA,3,l)*dz
          db=grdc(ICVLB,1,l)*dx+grdc(ICVLB,2,l)*dy+grdc(ICVLB,3,l)*dz
          qfv(ICFL,1)=q(ICVLA,l)+0.5d0*phi(ICVLA)*((1.d0
     &       -MUSCL_k(IMAT)*phi(ICVLA))*da+MUSCL_k(IMAT)*phi(ICVLA)*dab)
          qfv(ICFL,2)=q(ICVLB,l)-0.5d0*phi(ICVLB)*((1.d0
     &      -MUSCL_k(IMAT)*phi(ICVLB))*db+MUSCL_k(IMAT)*phi(ICVLB)*dab)
 430      continue
        endif
 420    continue
!
!----------------------------------------------------------
!-< 4. Upstream Shift Interpolation (usi) 2nd Scheme >
!----------------------------------------------------------
!
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        if(icnv(IMAT)==usi2nd) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          dum1=SFAREA(1,ICFL)
          dum2=SFAREA(2,ICFL)
          dum3=SFAREA(3,ICFL)
          dab=q(ICVLA,1)*dum1+q(ICVLA,2)*dum2+q(ICVLA,3)*dum3
          dx=dum1*dab*deltt
          dy=dum2*dab*deltt
          dz=dum3*dab*deltt
!
          dq(1)=
     &          (grdc(ICVLA,1,l)*dx
     &          +grdc(ICVLA,2,l)*dy
     &          +grdc(ICVLA,3,l)*dz)
          dab=q(ICVLB,1)*dum1+q(ICVLB,2)*dum2+q(ICVLB,3)*dum3
          dx=dum1*dab*deltt
          dy=dum2*dab*deltt
          dz=dum3*dab*deltt
          dq(2)=
     &        (grdc(ICVLB,1,l)*dx
     &        +grdc(ICVLB,2,l)*dy
     &        +grdc(ICVLB,3,l)*dz)
!
          qfv(ICFL,1)=q(ICVLA,l)+dq(1)
          qfv(ICFL,2)=q(ICVLB,l)-dq(2)
          enddo
        endif
        enddo

!
!-< 5. Calculate convection term >-
!
        do 300 nb=1,nbcnd    !IBF=1,NSSFBC
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 300
        kd=kdbcnd(0,nb)
        if(kd==kdprdc) cycle
        if(icnv(IMAT).eq.cnt2nd) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          qfv(ICFL,1)=q(ICV,l)
          qfv(ICFL,2)=q(IDC,l)
          enddo
        elseif(icnv(IMAT)/=up1st) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          qfv(ICFL,2)=q(IDC,l)
          enddo
        endif
        if(kd.eq.kdsld) then
          IIMAT=MAT_BCIDX(nb,1)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          qfv(ICFL,2)=q(IDC,l)
          qfv(ICFL,1)=q(ICV,l)
          enddo
        endif
  300   enddo
!
        cnvq(:)=0.d0
        cnvv(:)=0.d0
        do 310 IIMAT=1,NMAT    !ICF=1,NCVFAC 
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        if(icnv(IMAT).eq.engcsv) then
!          do ICFL=ICFS,ICFE
!          ICVLA=LVEDGE(1,ICFL)
!          ICVLB=LVEDGE(2,ICFL)
!          dum1=rva(ICFL)*qfv(ICFL,1)
!          dum2=rva(ICFL)*qfv(ICFL,2)
!          cnvq(ICVLA)=cnvq(ICVLA)+dum1
!          cnvq(ICVLB)=cnvq(ICVLB)-dum2
!          cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)
 !         cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)
!          enddo
        else
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dum2=aks(ICVA,iaph(iphs))*wi1+aks(ICVB,iaph(iphs))*wi2
          dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &        +min(0.d0,rva(ICFL))*qfv(ICFL,2)
          cnvq(ICVA)=cnvq(ICVA)+dum1!*dum2
          cnvq(ICVB)=cnvq(ICVB)-dum1!*dum2
          cnvv(ICVA)=cnvv(ICVA)+rva(ICFL)!*dum2
          cnvv(ICVB)=cnvv(ICVB)-rva(ICFL)!*dum2
          enddo
        endif
  310   enddo
!
! --- d(r*ui*uj)=uj*d(r*ui)+r*ui*d(uj):
!
        do 320 IIMAT=1,NMAT         !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 320
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 320
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        dqdt(ICVL,l)=dqdt(ICVL,l)+(q(ICVL,l)*cnvv(ICVL)-cnvq(ICVL))
 330    continue
 320    continue
!
!        if(ieul2ph>0) then
!          do IIMAT=1,NMAT
!          IMAT=MAT_NO(IIMAT)
!          if(IMAT.lt.0) cycle
!          ICVS=MAT_CVEXT(IIMAT-1)+1
!          ICVE=MAT_CVEXT(IIMAT)
!          do ICVL=ICVS,ICVE
!          dvdd(ICVL,l)=q(ICVL,l)*cnvv(ICVL)
!          enddo
!          enddo
!        endif
!
 1000   continue
      endif
!
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (3,MXALLCV,NCV,dqdt(:,1:3))
       ENDIF
!
      return
      end subroutine conv_term_vel_e2p

!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine conv_term_vof
     & (deltt,LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  qfv,rho,
     &  grdc,rva,vof,dqdt,vel,
     &  IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     qfv(1:2,:)  <=  rvd(1:2,:) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine if for Calculating convective term
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_boundary,only: nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                          kdolet
      use module_metrix,only  : phi =>W1K8
      use module_metrix,only  : cnvv=>W1K9
      use module_metrix,only  : cnvq=>W1K7 !cnvq=>W1K10
      use module_vof,   only  : cnvvof,up1st,up2nd,cnt2nd,cicsam,
     &                          den1,den2,bldf_vof,hric,
     &                          donr,intervof
      use module_initial ,only : rho0,rho02
      
!
      implicit none
!
! 1. Calculate convection term
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt
      integer,intent(in)    :: IMODE
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT (3,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      real*8 ,intent(inout) :: qfv    (  MXCVFAC,2)
      real*8 ,intent(in)    :: rva    (  MXCVFAC)
      real*8 ,intent(in)    :: vof    (  MXALLCV)
      real*8 ,intent(in)    :: rho    (  MXALLCV)
      real*8 ,intent(inout) :: dqdt   (  MXALLCV)
      real*8 ,intent(inout) :: grdc   (  MXALLCV,3)
      real*8 ,intent(in)    :: vel    (  MXALLCV,3)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      real*8  :: dq(2),vof_f,ru,rv,rw,gf2,vof1A,vof2A,vof1B,vof2B
      real*8 ,parameter :: ak=0.3333d0
      real*8  :: dqp,dqm,dqpp,dqmm,dqpm,wi1,wi2
      real*8  :: dx,dy,dz,dab,da,db,donrsig
      real*8  :: epsw,sgn,dum1,qqq,alhpa,
     &        dum2,dum3,dum4,dum5,dumC,dumflux
      real*8  :: dd,dl,dum_U,the,CD,CfUQ,Gammf,Betaf,CfCBC,Cf
      integer :: i,j,k,l,m,n
      integer :: nb,kd,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD,ICVDN,ICVAC
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVA,ICVB,ICV,IDC
      integer :: IBFS,IBFE,IBFL
      integer :: ierr=0
      logical :: icnTVD
      real*8  :: vn,cu,alp_d,alp_cbc,alp_uq,alp_f,gamma_f,beta_f 
      real*8  :: dum0,dum_cos,dum_cos2,root_cos,vof_U               
      integer :: ICVLA,ICVLB,ICVLA2,ICVLB2,ICV_U,ICV_D,ICV_A,imode_c
!
! --- START
!
      imode_c=1  !=0
      if(.NOT.(cnvvof==cicsam.or.cnvvof==hric)) then
!
        if(cnvvof==up1st) then
          do IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          if(IMAT.lt.0) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
            ICVA=LVEDGE(1,ICFL)
            ICVB=LVEDGE(2,ICFL)
            qfv(ICFL,1)=vof(ICVA)
            qfv(ICFL,2)=vof(ICVB)
          enddo
          enddo
        elseif(cnvvof==cnt2nd) then
          do IIMAT=1,NMAT    !ICF=1,NCVFAC
          if(.not.mat_cal(IIMAT)) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          IMAT=MAT_NO(IIMAT)
          if(IMAT.lt.0) cycle
          dum1=bldf_vof
          dum2=1.D0-dum1
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=vof(ICVA)*wi1+vof(ICVB)*wi2
          qfv(ICFL,1)=dum1*qqq+dum2*vof(ICVA)
          qfv(ICFL,2)=dum1*qqq+dum2*vof(ICVB)
          enddo
          enddo
        else!if(cnvvof==up2nd) then
          call FFRABORT(1,'ERR:VOF NOT support 2nd upwind shceme')
        endif
!
        do nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        kd=kdbcnd(0,nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        qfv(ICFL,2)=vof(IDC)
        enddo
        enddo
!
        cnvq=0.d0
        cnvv=0.d0
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &      +min(0.d0,rva(ICFL))*qfv(ICFL,2)
        dumC=wi1*rho(ICVA)+wi2*rho(ICVB)

        cnvq(ICVA)=cnvq(ICVA)+dum1/dumC
        cnvq(ICVB)=cnvq(ICVB)-dum1/dumC
        cnvv(ICVA)=cnvv(ICVA)+rva(ICFL)/dumC
        cnvv(ICVB)=cnvv(ICVB)-rva(ICFL)/dumC

!        dum2=max(0.d0,rva(ICFL))*qfv(ICFL,2)
!     &      +min(0.d0,rva(ICFL))*qfv(ICFL,1)
!        cnvq(ICVA)=cnvq(ICVA)+(dum1+dum2)/dumC
!        cnvq(ICVB)=cnvq(ICVB)-(dum1+dum2)/dumC
!        cnvv(ICVA)=cnvv(ICVA)+rva(ICFL)/dumC
!        cnvv(ICVB)=cnvv(ICVB)-rva(ICFL)/dumC
        enddo
        enddo
!
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dqdt(ICVL)=dqdt(ICVL)+(vof(ICVL)*cnvv(ICVL)-cnvq(ICVL))
!        dqdt(ICVL)=dqdt(ICVL)-cnvq(ICVL)
        enddo
        enddo
!
      elseif(cnvvof==cicsam.and.imode_c==0) then
!
! --- 
!
        call grad_cell(1,8,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vof,grdc)
!
!-< 1. 1st order upwind >-
!
        do 100 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do 110 ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        if(rva(ICFL)>0.d0) then
          ICVDN=ICVA
          ICVAC=ICVB
        else
          ICVDN=ICVB
          ICVAC=ICVA
        endif
        dum1=vof(ICVAC)
        dx=(CVCENT(1,ICVAC)-CVCENT(1,ICVDN))
        dy=(CVCENT(2,ICVAC)-CVCENT(2,ICVDN))
        dz=(CVCENT(3,ICVAC)-CVCENT(3,ICVDN))
        dl=dsqrt(dx**2+dy**2+dz**2)
        dd=dsqrt(grdc(ICVDN,1)**2+grdc(ICVDN,2)**2+grdc(ICVDN,3)**2)
        dum2=(dx*grdc(ICVDN,1)
     &       +dy*grdc(ICVDN,2)
     &       +dz*grdc(ICVDN,3))
        dum_U=max(0.d0,min((dum1-2.d0*dum2),1.d0))
        CD=(vof(ICVDN)-dum_U)/(vof(ICVAC)-dum_U+SML)

        Cu=abs(rva(ICFL)*deltt/(CVVOLM(ICVA)*rho(ICVDN)))
        if(CD>0.d0.and.CD<1.d0) then
          CfCBC=min(1.d0,CD/(Cu+SML))
          CfUQ=min(8.d0*Cu+(1.d0-Cu)*(6.d0*CD+3.d0),CfCBC)
        else
          CfCBC=CD
          CfUQ=CD
        endif
        the=acos(abs(dum2)/(dl*dd+SML))
        Gammf=min((cos(2.d0*the)+1.d0)*0.5d0,1.d0)
        Cf=Gammf*CfCBC+(1.d0-Gammf)*CfUQ
        Betaf=(Cf-CD)/(1.d0-CD+SML)
        qfv(ICFL,1)=(1.d0-Betaf)*vof(ICVDN)+Betaf*vof(ICVAC)
!        qfv(ICFL,2)=Betaf*vof(ICVAC)
 110    enddo
 100    enddo
!------------------------------------
! --- 
!
      elseif(cnvvof==cicsam.and.imode_c==1) then !suginaka
        call grad_cell(1,8,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vof,grdc)
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=SFAREA(1,ICFL)
        dum2=SFAREA(2,ICFL)
        dum3=SFAREA(3,ICFL)
        vn=vel(ICVLA,1)*dum1+vel(ICVLA,2)*dum2+vel(ICVLA,3)*dum3
        vn=abs(vn)   
        do i=1,2
        if(i.eq.1) then  
          ICV_D=ICVLA
          ICV_A=ICVLB
        else
          ICV_D=ICVLB
          ICV_A=ICVLA
        end if
        dx=CVCENT(1,ICV_A)-CVCENT(1,ICV_D)
        dy=CVCENT(2,ICV_A)-CVCENT(2,ICV_D)
        dz=CVCENT(3,ICV_A)-CVCENT(3,ICV_D)
        dum1=(dx*dx+dy*dy+dz*dz)**0.5
!        vof_U=vof(ICV_A)-2.d0*grdc(ICV_D)*dum1
        vof_U=vof(ICV_A)-2.d0*vof(ICV_D)*dum1
        vof_U=max(min(vof_U,1.d0),0.d0)
!        dum0=abs(vof(ICV_A)-vof(ICV_U))
        dum0=abs(vof(ICV_A)-vof_U)
        if(dum0 > 1.d-6) then
!          alp_d=(vof(ICV_D)-vof(ICV_U))/(vof(ICV_A)-vof(ICV_U))  
          alp_d=(vof(ICV_D)-vof_U)/(vof(ICV_A)-vof_U)  
          if(CVVOLM(ICV_D) < 1.d-6) then
            cu=1.d-6
          else
            cu=max(vn*SFAREA(4,ICFL)*deltt/CVVOLM(ICV_D),1.d-6)
          end if
! Hyper-C
          if(alp_d >= 0.d0 .and. alp_d <=  1.d0) then
            alp_cbc=min(alp_d/cu,1.d0)
          else 
            alp_cbc=alp_d   
          end if
! Ultimate-Quickest
          if(alp_d >= 0.d0 .and. alp_d <=  1.d0) then
            dum1=(8.d0*cu*alp_d+
     &          (1.d0-cu)*(6.d0*alp_d+3.d0))/8.d0
            alp_uq=min(dum1,alp_cbc)
          else 
            alp_uq=alp_d   
          end if
          dum1=-(grdc(ICV_D,1)*SFAREA(1,ICFL)
     &          +grdc(ICV_D,2)*SFAREA(2,ICFL)
     &          +grdc(ICV_D,3)*SFAREA(3,ICFL))
          dum2= (grdc(ICV_D,1)*grdc(ICV_D,1)
     &          +grdc(ICV_D,2)*grdc(ICV_D,2)
     &          +grdc(ICV_D,3)*grdc(ICV_D,3))**0.5 + 1.d-6
          dum3=dum1/dum2
          dum_cos=abs(dum3)
          dum_cos2=2.d0*dum_cos*dum_cos-1.d0
          gamma_f=min((dum_cos2+1.d0)*0.5d0,1.d0)
          alp_f=gamma_f*alp_cbc+(1.d0-gamma_f)*alp_uq
! local boundedness criteria
          if(alp_d >= 0.d0 .and. alp_d <= 1.d0 ) then
            alp_f=min(alp_f,1.d0)
            alp_f=max(alp_f,alp_d)
          else
            alp_f=alp_d
          end if
          if(alp_d > 0.99999) then
            beta_f=0.d0
          else
            beta_f=(alp_f-alp_d)/(1.d0-alp_d) 
          end if
          dq(i)=(1.d0-beta_f)*vof(ICV_D)+beta_f*vof(ICV_A)
        else
!          dq(i)=vof(ICV_D)  
          dq(i)=0.d0
        end if
        enddo
        qfv(ICFL,1)=dq(1)
        qfv(ICFL,2)=dq(2)
        enddo
        enddo
      elseif(cnvvof==hric) then !suginaka
!
        call grad_cell(1,8,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,vof,grdc)
!
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=SFAREA(1,ICFL)
        dum2=SFAREA(2,ICFL)
        dum3=SFAREA(3,ICFL)
        vn=vel(ICVLA,1)*dum1+vel(ICVLA,2)*dum2+vel(ICVLA,3)*dum3
        vn=abs(vn)   
        do i=1,2
        if(i.eq.1) then  
!          ICV_U=ICVLA2
          ICV_D=ICVLA
          ICV_A=ICVLB
        else
!          ICV_U=ICVLB2
          ICV_D=ICVLB
          ICV_A=ICVLA
        end if
        dx=CVCENT(1,ICV_A)-CVCENT(1,ICV_D)
        dy=CVCENT(2,ICV_A)-CVCENT(2,ICV_D)
        dz=CVCENT(3,ICV_A)-CVCENT(3,ICV_D)
        dum1=(dx*dx+dy*dy+dz*dz)**0.5
        vof_U=vof(ICV_A)-2.d0*vof(ICV_D)*dum1
        vof_U=max(min(vof_U,1.d0),0.d0)
!        dum0=abs(vof(ICV_A)-vof(ICV_U))
        dum0=abs(vof(ICV_A)-vof_U)
        if(dum0 > 1.d-6) then
!          alp_d=(vof(ICV_D)-vof(ICV_U))/(vof(ICV_A)-vof(ICV_U))  
          alp_d=(vof(ICV_D)-vof_U)/(vof(ICV_A)-vof_U)  
          if(alp_d < 0.d0) then
            alp_f=alp_d
          else if(alp_d < 0.5d0) then
            alp_f=2.d0*alp_d
          else if(alp_d < 1.d0 ) then
            alp_f=1.d0
          else
            alp_f=alp_d
          end if
          if(CVVOLM(ICV_D) < 1.d-6) then
            cu=0.d0
          else
            cu=vn*SFAREA(4,ICFL)*deltt/CVVOLM(ICV_D)
          end if
          if(cu >= 0.3d0 .and. cu < 0.7d0) then
            alp_f=alp_d+(alp_f-alp_d)*(0.7d0-cu)/0.4d0
          else if(cu >= 0.7d0) then
            alp_f=alp_d
          end if
          dum1=-(grdc(ICV_D,1)*SFAREA(1,ICFL)
     &          +grdc(ICV_D,2)*SFAREA(2,ICFL)
     &          +grdc(ICV_D,3)*SFAREA(3,ICFL))
!
          dum2= (grdc(ICV_D,1)*grdc(ICV_D,1)
     &          +grdc(ICV_D,2)*grdc(ICV_D,2)
     &          +grdc(ICV_D,3)*grdc(ICV_D,3))**0.5+1.d-5
          dum3=dum1/dum2
          dum3=abs(dum3)
          root_cos=dum3**0.5
          alp_f=alp_f*root_cos+alp_d*(1.d0-root_cos)
!          dq(i)=alp_f*(vof(ICV_A)-vof(ICV_U))+vof(ICV_U)
          dq(i)=alp_f*(vof(ICV_A)-vof_U)+vof_U
        else
          dq(i)=0.d0
        end if
        end do
        qfv(ICFL,1)=dq(1)
        qfv(ICFL,2)=dq(2)
        enddo
        enddo
      endif
!
      if(cnvvof==cicsam.or.cnvvof==hric) then
        do 300 nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) goto 300
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 300
        kd=kdbcnd(0,nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        qfv(ICFL,2)=vof(IDC)
        enddo
  300   continue
!
        cnvq=0.d0
        cnvv=0.d0
        do 310 IIMAT=1,NMAT    !ICF=1,NCVFA
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) goto 310 
        if(IMAT.lt.0) goto 310
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
!---------------------------------------------------
        vof1A=vof(ICVA)
        vof2A=1.d0-vof1A
        vof1B=vof(ICVB)
        vof2B=1.d0-vof1B
!---------------------------------------------------
        dumC=wi1*rho(ICVA)+wi2*rho(ICVB)
        dumflux=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &         +min(0.d0,rva(ICFL))*qfv(ICFL,2)
!---------------------------------------------------
       cnvq(ICVA)=cnvq(ICVA)+dumflux/dumC
       cnvq(ICVB)=cnvq(ICVB)-dumflux/dumC
       cnvv(ICVA)=cnvv(ICVA)+rva(ICFL)/dumC
       cnvv(ICVB)=cnvv(ICVB)-rva(ICFL)/dumC
!---------------------------------------------------
        enddo
  310   continue
!
! --- d(r*ui*uj)=uj*d(r*ui)+r*ui*d(uj):
!
        do 320 IIMAT=1,NMAT   !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 320
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) goto 320
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 330 ICVL=ICVS,ICVE
        dqdt(ICVL)=dqdt(ICVL)+(VOF(ICVL)*cnvv(ICVL)-cnvq(ICVL))
!        dqdt(ICVL)=dqdt(ICVL)-cnvq(ICVL)
 330    continue
 320    continue
      endif
!

!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,dqdt)
      ENDIF
!
      return
      end subroutine conv_term_vof
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine conv_term_e2p
     & (IMD,icnv,lmtr,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  qfv,vel,diagq,
     &  grdc,rva,rho,q,dqdt,rvatmp,
!     #  qp,
     &  IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     qfv(1:2,:)  <=  rvd(1:2,:)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine if for Calculating convective term
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_boundary,only: nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &                          kdolet,kdilet
      use module_metrix,only  : phi =>W1K8
      use module_metrix,only  : cnvv=>W1K9
      use module_metrix,only  : cnvq=>W1K7 !cnvq=>W1K10
      USE module_dimension
      use module_material,only : nflud,venkat,rnck,MUSCL_k,slope,
     &                           up1st,up2nd,up3rd,cnt2nd,bldf,
     &                           usi2nd,engcsv ,mscl
      use module_Euler2ph,only : ieul2ph
      use module_scalar,  only : iaph
!
      implicit none
!
! 1. Calculate convection term 
!
! --- [dummy arguments] 
!
      real*8 ,intent(in)    :: deltt
      integer,intent(in)    :: IMD,icnv(nflud),lmtr(nflud),IMODE
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV   (MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT(3,MXCVFAC)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(  MXALLCV)
      real*8 ,intent(inout) :: qfv   (  MXCVFAC,2)
      real*8 ,intent(inout) :: rva   (  MXCVFAC)
      real*8 ,intent(inout) :: rvatmp(  MXCVFAC)
      real*8 ,intent(in)    :: q     (  MXALLCV)
      real*8 ,intent(in)    :: rho   (  MXALLCV)
      real*8 ,intent(inout) :: dqdt  (  MXALLCV)
      real*8 ,intent(inout) :: diagq (  MXALLCV)
      real*8 ,intent(inout) :: grdc  (  MXALLCV,3)
      real*8 ,intent(in)    :: vel   (  MXALLCV,3)
      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!      real*8 ,intent(in)    :: qp    (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: dq  (2)
      real*8 ,parameter :: ak=0.3333d0
      real*8  :: dqp,dqm,dqpp,dqmm,dqpm,wi1,wi2
      real*8  :: dx,dy,dz,dab,da,db
      real*8  :: epsw,sgn,dum1,qqq,alhpa,dum2,dum3,ru,rv,rw
      integer :: i,j,k,l,m,n
      integer :: nb,kd,kdt,kdy,kdk,kdp
      integer :: ICOM
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC
      integer :: IBFS,IBFE,IBFL
      integer :: ierr=0
      logical :: icnTVD
      real*8 ,save :: vol_cov,dxx,dyx,dzx,dl
!
! --- START
!
      icnTVD=.false.

      do 1000 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 1000
      IF(icnv(IMAT)==up2nd.OR.
     &   icnv(IMAT)==up3rd.OR.
     &   icnv(IMAT)==mscl.OR.
     &   icnv(IMAT)==usi2nd
!     &   IMODE==29.or.IMODE==30
     &  ) THEN
        icnTVD=.true.
      endif
 1000 continue
      if(NPE.gt.1) call hpclor(icnTVD)
      if(icnTVD) then
        call grad_cell(1,6,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,q,grdc)
      endif
!
!-< 1. 1st order upwind >-
!
      do 100 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 100
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 100
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(icnv(IMAT).eq.up1st) then
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        qfv(ICFL,1)=q(ICVLA)
        qfv(ICFL,2)=q(ICVLB)
        enddo
      endif
  100 enddo
!
!-< 1. 2nd order centeral >-: bldf=0.D0: 1stUD; bldf=1.D0: 2ndCD;
!
      do 150 IIMAT=1,NMAT  
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      if(icnv(IMAT).eq.cnt2nd) then
        dum1=bldf(IMAT)
        dum2=1.D0-dum1
        if(IMODE.eq.1) then
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)

          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=q(ICVLA)*wi1+q(ICVLB)*wi2
          qfv(ICFL,1)=q(ICVLA)
          qfv(ICFL,2)=q(ICVLB)
          enddo
        else
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          qqq=q(ICVLA)*wi1+q(ICVLB)*wi2
          qfv(ICFL,1)=dum1*qqq+dum2*q(ICVLA)
          qfv(ICFL,2)=dum1*qqq+dum2*q(ICVLB)
          enddo
        endif
      endif
  150 enddo
!
!-< 3. Higher order with Venkatakrishnan limiter >-
!
!
!--< 3.1 calculate gradient in cell >--
!
!
!--< 3.2 calculate value on face with limited gradient >--
!        qmax(:)<=cnvq(:)
!        qmin(:)<=cnvv(:)
!
      phi(:)=1.d0
      do 250 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) goto 250
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 250
      if(icnv(IMAT)==up2nd.or.icnv(IMAT)==up3rd.or.
     &   icnv(IMAT)==mscl) then
        if(lmtr(IMAT).eq.venkat.or.lmtr(IMAT).eq.slope) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          cnvq(ICVS:ICVE)=q(ICVS:ICVE)
          cnvv(ICVS:ICVE)=q(ICVS:ICVE)
          cnvq(IDCS:IDCE)=q(IDCS:IDCE)
          cnvv(IDCS:IDCE)=q(IDCS:IDCE)
          do 255 ICFL=ICFS,ICFE
          do 260 I=1,2
          ICVL=LVEDGE(I,ICFL)
          dq(I)=grdc(ICVL,1)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &         +grdc(ICVL,2)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &         +grdc(ICVL,3)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
 260      continue
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          do 257 I=1,2
          ICVL=LVEDGE(I,ICFL)
          if(abs(dq(i)).gt.SML) then
            cnvq(ICVL)=max(cnvq(ICVL),q(ICVLA),q(ICVLB))
            cnvv(ICVL)=min(cnvv(ICVL),q(ICVLA),q(ICVLB))
            dqm=dq(i)
            epsw=rnck*rnck*rnck*CVVOLM(ICVL)
            sgn=sign(1.d0,dqm)
            dqp=max(sgn*(cnvq(ICVL)-q(ICVL)),
     &              sgn*(cnvv(ICVL)-q(ICVL)))
            dqpp=dqp*dqp
            dqmm=dqm*dqm
            dqpm=abs(dqp*dqm)
            if(lmtr(IMAT).eq.venkat) then
              phi(ICVL)=min(phi(ICVL),
     &         (dqpp+dqpm+dqpm+epsw)/(dqpp+dqmm+dqmm+dqpm+epsw+SML))
            elseif(lmtr(IMAT).eq.slope) then
              epsw=1.d-3*rnck
              phi(ICVL)=min(phi(ICVL),(dqpm+dqpm+epsw)/(dqpp+dqmm+epsw))
            endif
	  else
 	    phi(ICVL)=1.0d0
          endif
 257      continue
 255      continue
        endif
      elseif(icnv(IMAT)==usi2nd) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          phi(ICVS:ICVE)=1.d0
      endif
 250  continue
!
      do 200 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 200
      if(icnv(IMAT).eq.up2nd) then
        cnvq(ICVS:ICVE)=q(ICVS:ICVE)
        cnvv(ICVS:ICVE)=q(ICVS:ICVE)
        do 210 ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        do 201 I=1,2
        ICVL=LVEDGE(I,ICFL)
        cnvq(ICVL)=max(cnvq(ICVL),q(ICVLA),q(ICVLB))
        cnvv(ICVL)=min(cnvv(ICVL),q(ICVLA),q(ICVLB))
        dq(I)=grdc(ICVL,1)*(SFCENT(1,ICFL)-CVCENT(1,ICVL))
     &       +grdc(ICVL,2)*(SFCENT(2,ICFL)-CVCENT(2,ICVL))
     &       +grdc(ICVL,3)*(SFCENT(3,ICFL)-CVCENT(3,ICVL))
  201   enddo
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        qfv(ICFL,1)=q(ICVLA)+phi(ICVLA)*dq(1)
        qfv(ICFL,2)=q(ICVLB)+phi(ICVLB)*dq(2)
 210    continue
      endif
 200  continue
!
!-< 4. 3rd TVD Scheme >
!
      do 420 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) goto 420
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 420
      if(icnv(IMAT)==up3rd.or.icnv(IMAT)==mscl) then
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dab=q(ICVLB)-q(ICVLA)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
        dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
        dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
        da=grdc(ICVLA,1)*dx+grdc(ICVLA,2)*dy+grdc(ICVLA,3)*dz
        db=grdc(ICVLB,1)*dx+grdc(ICVLB,2)*dy+grdc(ICVLB,3)*dz
        qfv(ICFL,1)=q(ICVLA)+0.5d0*phi(ICVLA)*((1.d0
     &     -MUSCL_k(IMAT)*phi(ICVLA))*da+MUSCL_k(IMAT)*phi(ICVLA)*dab)
        qfv(ICFL,2)=q(ICVLB)-0.5d0*phi(ICVLB)*((1.d0
     &     -MUSCL_k(IMAT)*phi(ICVLB))*db+MUSCL_k(IMAT)*phi(ICVLB)*dab)
        enddo
      endif
  420 continue
!----------------------------------------------------------
!-< 4. Upstream Shift Interpolation (usi) 2nd Scheme >
!----------------------------------------------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      if(icnv(IMAT)==usi2nd) then
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=SFAREA(1,ICFL)
        dum2=SFAREA(2,ICFL)
        dum3=SFAREA(3,ICFL)
        dab=vel(ICVLA,1)*dum1+vel(ICVLA,2)*dum2+vel(ICVLA,3)*dum3
        dx=dum1*dab*deltt
        dy=dum2*dab*deltt
        dz=dum3*dab*deltt
        dq(1)=
     &        (grdc(ICVLA,1)*dx
     &        +grdc(ICVLA,2)*dy
     &        +grdc(ICVLA,3)*dz)
        dab=vel(ICVLB,1)*dum1+vel(ICVLB,2)*dum2+vel(ICVLB,3)*dum3
        dx=dum1*dab*deltt
        dy=dum2*dab*deltt
        dz=dum3*dab*deltt
        dq(2)=
     &        (grdc(ICVLB,1)*dx
     &        +grdc(ICVLB,2)*dy
     &        +grdc(ICVLB,3)*dz)
        qfv(ICFL,1)=q(ICVLA)+dq(1)
        qfv(ICFL,2)=q(ICVLB)-dq(2)
        enddo
      endif
      enddo
!
!-< 5. Calculate convection term >-
!
      do 300 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      kd=kdbcnd(0,nb)
      if(icnv(IMAT).eq.cnt2nd) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        qfv(ICFL,1)=q(ICV)
        qfv(ICFL,2)=q(IDC)
        enddo
      elseif(icnv(IMAT).ne.up1st) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        qfv(ICFL,2)=q(IDC)
        enddo
      else
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        qfv(ICFL,2)=q(IDC)
        enddo
      endif
  300 continue
!
!
      if(IMODE==29.or.IMODE==30) then
        cnvq=0.d0
        cnvv=0.d0
        phi=0.d0
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT.lt.0) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
!
        dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &      +min(0.d0,rva(ICFL))*qfv(ICFL,2)
!
        dum2=wi1/rho(ICVLA)+wi2/rho(ICVLB)
        cnvq(ICVLA)=cnvq(ICVLA)+dum1 !*dum2
        cnvq(ICVLB)=cnvq(ICVLB)-dum1 !*dum2
        cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL) !+rvatmp(ICFL) !*dum2
        cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL) !-rvatmp(ICFL) !*dum2
        enddo
        enddo
!
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dqdt(ICVL)=dqdt(ICVL)-cnvq(ICVL)!+cnvv(ICVL)*q(ICVL)
!        diagq(ICVL)=diagq(ICVL)+cnvv(ICVL) 
        enddo
        enddo
!
      else
      endif
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,dqdt)
      ENDIF
!
      return
      end subroutine conv_term_e2p
