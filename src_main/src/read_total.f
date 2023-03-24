!
!      subroutine read_cell
!      subroutine read_initial
!      subroutine read_novertex
!      subroutine read_source
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_cell(lacell,lvcell,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [dummy arguments]
!
      use module_dimension
      use module_io,only : lenfnm,ifle,getfil
!
! 1. Read cell file
!
! --- [dummy arguments]
!
      integer,intent(out) :: lacell(  mxcelb)
      integer,intent(out) :: lvcell(8,mxcell)
      integer,intent(out) :: ierror
!
!
! --- [local entities]
!
      integer,allocatable :: mncv(:)
!
      character(lenfnm) :: fnam
      integer :: idat(3)
      integer :: i,k,n,IC
      integer :: nrec,ncelx,ifl,ios,jerr,idum1
!
      allocate(mncv(mxcell))
!
!
!-< 1. Initial set >-
!
      ierror=0
!
      call getfil(ifl,fnam,'cell')
!
      do 100 IC=1,ncell
      lacell(IC)=-999999
      do 101 i=1,8
      lvcell(i,IC)=0
  101 continue
  100 continue
!
!
!-< 2. Input data >-
!
      nrec =0
      ncelx=0
  200 continue
      read(ifl,iostat=ios) (idat(i),i=1,3)
      if(ios.lt.0) goto 201
      call errmsg0
      if(ierror.ne.0) goto 9999
      ncelx=ncelx+1
      n=min(ncelx,ncell)
      k=min(idat(2),8)
      lacell(n)=idat(3)
      mncv   (n)=idat(2)
      if(k.lt.1) goto 9001
      read(ifl,iostat=ios) (lvcell(i,n),i=1,k)
      call errmsg0
      if(ierror.ne.0) goto 9999
      goto 200
  201 continue
!
      if(ncelx.ne.ncell) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of cells is not matched'
        write(ifle,*) 'in cell file    :',ncelx
        write(ifle,*) 'in control data :',ncell
        write(ifle,*) 'name of cell file =',fnam(:len_trim(fnam))
        goto 9999
      endif
!
!
!-< 3. Check data >-
!
      do 300 n=1,ncell
!
! / attribute of cell /
      if( lacell(n).eq.0 .or. lacell(n).gt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'attribute of cell =',lacell(n)
        write(ifle,*) 'it must be < 0 or = 1'
        write(ifle,*) 'cell no. =',n
        write(ifle,*) 'name of cell file =',fnam(:len_trim(fnam))
        goto 9999
      endif
!
! / no. of vertices /
      if( mncv(n).gt.8 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of vertices =',mncv(n)
        write(ifle,*) 'it must be < 9'
        write(ifle,*) 'cell no. =',n
        write(ifle,*) 'name of cell file =',fnam(:len_trim(fnam))
        goto 9999
      endif
!
! / range of vertex no. /
      do 301 i=1,mncv(n)
      if( lvcell(i,n).lt.1 .or. lvcell(i,n).gt.nvrtx ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'vertex no. =',lvcell(i,n)
        write(ifle,*) 'it must be > 0 and <=',nvrtx
        write(ifle,*) 'cell no. =',n
        write(ifle,*) 'name of cell file =',fnam(:len_trim(fnam))
        goto 9999
      endif
  301 continue
!
! / duplicated vertex no. /
      do 302 i=1,mncv(n)-1
      do 303 k=i+1,mncv(n)
      if( lvcell(i,n).eq.lvcell(k,n) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'the cell has duplicated vertex no.'
        write(ifle,*) 'vertex no. =',lvcell(i,n)
        write(ifle,*) 'cell no. =',n
        write(ifle,*) 'name of cell file =',fnam(:len_trim(fnam))
        goto 9999
      endif
  303 continue
  302 continue
!
  300 continue
!
      deallocate(mncv)
      call debug
      return
!
 9001 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'no. of vertices =',k
      write(ifle,*) 'it must be > 0'
      write(ifle,*) 'cell no. =',n
      write(ifle,*) 'name of cell file =',fnam(:len_trim(fnam))
      goto 9999
 9999 continue
      write(ifle,*) '(read_cell)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine errmsg0
!=================================================
      nrec=nrec+1
      if( ios.eq.0 ) return
      if( ierror.ne.0 ) return
      if( ios.gt.0 ) then
        write(ifle,*) '### error-1 : file read error'
      else
        write(ifle,*) '### error : end of file detected unexpectedly'
      endif
      write(ifle,*) 'file name = ',fnam(:len_trim(fnam))
      write(ifle,*) 'record no.2 =',nrec
      write(ifle,*) 'iostat =',ios
      ierror=1
      end subroutine errmsg0
!=================================================
      subroutine debug
!=================================================
      use module_debug,only : idebug
      if( idebug(1).eq.0 ) return
!
      end subroutine debug
!
      end subroutine read_cell
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_initial(
     & iters,times,
     & CVVOLM,CVCENT,DISALL,
     & MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,mat_cfidx,
     & LBC_SSF,LCYCSF,mat_cal,LVEDGE,wiface,
     & dvdd,pp0,prs,tmp,yys,vel,
     & velf,rva,
     & aks,rho,hhh,
     & MHD_A,MHD_FAI,SIGMA,MHD_CRT0,MHD_CRNT,
     & dvdd2,rho2,tmp2,yys2,
     & vel2,PRS2,rva2,hhh2,
     & MOLFRC,SITDEN,BLKTHK,wifsld,LCYCOLD,OPPANG,
     & cord,lacell,lvcell,
     & POTNAL,WDRBND,
     & ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     reading restart file
!
! --- [dummy arguments]
!
      use module_hpcutil
      use module_dimension
!      use module_constant
      use module_io,only      : lenfnm,ifle,getfil,ifll
      use module_model,only   : idrdp,mach0,comp,incomp,icaltb,
     &                          sles,dles,lles,SDES,KE2S
      use module_rans,only    : akslw,sgmaks
      use module_dimnsn,only  : xredc,yredc,zredc
      USE module_usersub,ONLY : iniusr,usrno,usryes,ini_MHD_user,
     &                          iniusr_E2P,usryesyes
      use module_material,only: nflud,nofld
      use module_material,only: sthrmu
      use module_initial 
      use module_scalar  ,only: rns_scl,iaph,ike,ivof,ides,ivold,icavi,
     &                          ical_cavi,
     &                          ical_FC,ical_s
      use module_model   ,only: icaltb,ke,ke_low,RNG,CHEN,ical_MHD,
     &                          ical_mvmsh,iLBF_P
      use module_Euler2ph,only: ieul2ph
      use module_VOF,     only: ical_vof
      use module_chemreac,ONLY: ical_suf
      use module_material,only: nsold,nosld,rsld,ical_sld,rotati,
     &                          ishaft,end,begin,
     &                          rot_ang,OMGA
      use module_boundary,only: nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,kdcvd,
     &                          phs_idx,phs_dns,
     &                          phs_nbc,phs_typ,phs_snm,phs_nam,
     &                          gasphase,surphase,blkphase,phs_com,
     &                          phs_inifrac,blk_dns,phs_thk,surfreac
      use module_species,only : r_wm,spcnam,ACT
      use module_scalar,only  : pot_scl,ip_pefc,potname
      use module_metrix,only  : STA_SAV
      use module_les   ,only   : ISTART,n_ave,ista_momery
      use module_model,   only : ical_WDR
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: iters
      real*8 ,intent(out)   :: times
      INTEGER,INTENT(IN)    :: MAT_NO(    0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV(    MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX( 0:MXMAT)
      integer,intent(in)    :: LVEDGE (2, MXCVFAC)
      logical,INTENT(IN)    :: mat_cal(   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF(   MXSSFBC)
      INTEGER,INTENT(IN)    :: LCYCSF (   MXSSFBC)
      real*8 ,intent(in)    :: CVVOLM (   MXALLCV)
      real*8 ,intent(in)    :: CVCENT (3, MXALLCV)
      real*8 ,intent(in)    :: DISALL (      MXCV)
      real*8 ,intent(out)   :: dvdd   (      MXCV,3,2)
      real*8 ,intent(out)   :: pp0    (   MXALLCV)
      real*8 ,intent(out)   :: prs    (   MXALLCV)
      real*8 ,intent(out)   :: tmp    (   MXALLCV)
      real*8 ,intent(out)   :: hhh    (   MXALLCV)
      real*8 ,intent(out)   :: yys    (   MXALLCV,mxcomp)
      real*8 ,intent(out)   :: vel    (   MXALLCV,3)
      real*8 ,intent(out)   :: velf   (   MXCVFAC_B,3)
      real*8 ,intent(out)   :: rva    (   MXCVFAC)
      real*8 ,intent(out)   :: rho    (   MXALLCV)
      real*8 ,intent(out)   :: aks    (   MXALLCVR,mxrans)
      REAL*8, INTENT(out)   :: POTNAL  (  MXALLCVP,MXPOTN)
      real*8 ,intent(out)   :: rho2   (   MXALLCVC)
      real*8 ,intent(out)   :: tmp2   (   MXALLCV2)
      real*8 ,intent(out)   :: hhh2   (   MXALLCV2)
      real*8 ,intent(out)   :: yys2   (   MXALLCV2,MXCOMP)
      real*8 ,intent(out)   :: vel2   (   MXALLCV2,3)
      real*8 ,intent(out)   :: PRS2   (   MXALLCV2)
      real*8 ,intent(out)   :: rva2   (   MXCVFAC2)
      real*8 ,intent(out)   :: dvdd2  (   MXCV2   ,3,2)
      real*8 ,intent(out)   :: MOLFRC ( MXSSFBC_SUF,MXCOMPALL,2)
      real*8 ,intent(inout) :: SITDEN ( MXSSFBC_SUF,MXPHASE)
      real*8 ,intent(inout) :: BLKTHK ( MXSSFBC_SUF,MXPHASE)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld ( MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG ( MXSSFBC_SLD)
      real*8 ,intent(inout) :: MHD_A  ( MXCV_MHD,3,2,2)
      real*8 ,intent(inout) :: MHD_FAI( MXCV_MHD,2,2)
      real*8 ,intent(inout) :: SIGMA  ( MXCV_MHD)
      real*8 ,intent(inout) :: MHD_CRT0(MXCV_MHD,3,2)
      real*8 ,intent(inout) :: MHD_CRNT(MXCV_MHD,3,2)
      real*8 ,intent(inout) :: cord    (3,MXVRTX_m)
      integer,intent(inout) :: lacell  (  MXCELL_m)
      integer,intent(inout) :: lvcell  (8,MXCELL_m)
      REAL*8 ,INTENT(IN)    :: WIFACE(        MXCVFAC)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      integer,intent(inout) :: WDRBND    (   MXCV_WDR,N_inj)
!     
      integer,intent(out)   :: ierror
!
! -- [local entities]
!
      character(lenfnm) :: fnam
      integer :: ifl=-1,ios,nrec,iter,itera,ieul2phx=0,itersld=-1,
     &           ical_sufx,linike=0,ical_MHDx,IMVMSHx,ista_momeryx,
     &           ical_WDRx,N_injx
      integer :: ierr1,i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: NCVX,ncmpx,nrnsx,NCVFAX,IIMAT,IMAT,ICVS,ICVE,NALLCVX,
     &           NMATX,ical_sldx,ncompallx,npotnx
      integer :: ICV,ICOM,IMD,ICF,IMAT_U,ICFL,ICVL,IDCS,IDCE
      real*8  :: t1,y1(mxcomp),p1,u1,v1,w1,rho1,aks1(mxrans)
      real*8  :: potn1(mxpotn)
      real*8  :: vr(3),unit(3),radi(3),rbb(3,3),bb(3,3),dr
      real*8  :: t2,y2(mxcomp),u2,v2,w2,dens2,p2
      real*8  :: x,y,z,waldis,vol,rk,re
      real*8  :: time,sum,vmx,sgm,r1,r2,dum1,dum2
      real*8,parameter  :: SML=1.d-25
      logical :: inifil,ldum
      integer :: nb,kd,nbx,ntp,icoms,icome,iph,ipcom
      integer :: IBFS,IBFE,IBFL
      real(8) :: org_x,org_y,tha,costh,sinth,v00(3),v_temp(3)
      real*8  :: MHD_Ai(3,2),MHD_FAIi(2),MHD_SIGMAi,MHD_CURRENTi(3,2)
      integer :: IMHD,icfs,icfe,icva,icvb
      real*8  :: alpha(2)
!
!-< 1. Uniform flow field >-
!
      ierror=0
!
!--------------------------------------------------------------------
! --- iters < 0: uniform initial flow condition
! --- iters = 0: User defined initial flow condition by initial file
! --- iters = 0 & iniusr.eq.usryes : User by user subroutine
! --- iters > 0: restart from former running
!--------------------------------------------------------------------
!
      t0(:)=300.d0
      y0(:,1)=1.d0
      y0(:,2:ncomp)=0.d0
      p0(:)=0.d0
      u0(:)=0.d0
      v0(:)=0.d0
      w0(:)=0.d0
      rho0(:)=1.d0
      vsgm0(:,:)=0.d0
      aks0(:,:)=0.d0
      if(icaltb==SDES) then
        aks0(:,ides)=sthrmu(:)
      endif
      potn0(:,:)=0.d0
!
      t02(:)=300.d0
      y02(:,1)=1.d0
      y02(:,2:ncomp)=0.d0
      u02(:)=0.d0
      v02(:)=0.d0
      w02(:)=0.d0
      rho02(:)=1.d0
!
      fai_iMHD(:,:)=0.d0
      A_iMHD(:,:,:)=0.d0
      eps_iMHD(:)=0.d0
      eps0_iMHD(:)=0.d0
      mu_iMHD(:)=0.d0
      mu0_iMHD(:)=0.d0
      sigma_iMHD(:)=0.d0
      sigma0_iMHD(:)=0.d0
      Current_iMHD(:,:,:)=0.d0
      OMEGA_iMHD(:)=0.d0
!
      call nml_init
     &  (MAT_NO,t0,y0,p0,u0,v0,w0,rho0,vsgm0,aks0,potn0,
     &   t02,y02,u02,v02,w02,rho02,
     &   fai_iMHD,A_iMHD,eps_iMHD,eps0_iMHD,mu_iMHD,
     &   mu0_iMHD,sigma_iMHD,sigma0_iMHD,Current_iMHD,OMEGA_iMHD,
     &   ierr1)
!
      if(ierr1/=0) goto 9999
!
      if(iters.eq.-1.and.iniusr.eq.usryes) then
        write(ifle,'(1X,a)') 
     &   'ERR: Reset start=0 when using initial user subroutine'
        call FFRABORT(1,'initial')
      endif
!
!
      itersld=iters
      if(iters.gt.0) goto 3000   ! restart
!------------------------------------------------
! --- <1> namelist initialization (in fflow.ctl)
!------------------------------------------------
      iters=0
      times=0.d0
!
      do 110 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)   !  if(IMAT.lt.0) goto 110
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      IDCS=MAT_DCIDX(IIMAT-1)+1
      IDCE=MAT_DCIDX(IIMAT)
      prs(ICVS:ICVE)=0.d0
      pp0(ICVS:ICVE)=0.d0
      prs(IDCS:IDCE)=0.d0
      pp0(IDCS:IDCE)=0.d0
!
      if(IMAT.gt.0)  then
        rho(ICVS:ICVE)=rho0(IIMAT) 
        rho(IDCS:IDCE)=rho0(IIMAT) 
      elseif(IMAT.lt.0) then
        rho(ICVS:ICVE)=rsld(-IMAT)
        rho(IDCS:IDCE)=rsld(-IMAT)
      endif
      if(idrdp.eq.comp) then
        prs(ICVS:ICVE)=p0(IIMAT)
        prs(IDCS:IDCE)=p0(IIMAT)
      elseif(idrdp.eq.mach0.or.idrdp==incomp) then
        pp0(ICVS:ICVE)=p0(IIMAT)
        pp0(IDCS:IDCE)=p0(IIMAT)
      endif
      tmp(ICVS:ICVE)=t0(IIMAT)
      tmp(IDCS:IDCE)=t0(IIMAT)
      do 100 ICOM=1,ncomp
      if(IMAT.gt.0) then
        yys(ICVS:ICVE,ICOM)=y0(IIMAT,ICOM)
        yys(IDCS:IDCE,ICOM)=y0(IIMAT,ICOM)
      else
        yys(ICVS:ICVE,ICOM)=0.d0
        yys(IDCS:IDCE,ICOM)=0.d0
      endif
  100 continue
      vel(ICVS:ICVE,1)=u0(IIMAT)
      vel(ICVS:ICVE,2)=v0(IIMAT)
      vel(ICVS:ICVE,3)=w0(IIMAT)
      vel(IDCS:IDCE,1)=u0(IIMAT)
      vel(IDCS:IDCE,2)=v0(IIMAT)
      vel(IDCS:IDCE,3)=w0(IIMAT)
      if(min(vsgm0(IIMAT,1),vsgm0(IIMAT,2)).gt.0.d0) then
        sgm=vsgm0(IIMAT,1)
        vmx=vsgm0(IIMAT,2)
        do 102 l=1,3
        call utl_boxmlr(r1,r2)
        vel(ICVS:ICVE,l)=vel(ICVS:ICVE,l)+max(-vmx,min(vmx,sgm*r1))
        vel(IDCS:IDCE,l)=vel(IDCS:IDCE,l)+max(-vmx,min(vmx,sgm*r1))
  102   continue
      endif
      if(rns_scl) then
        if(IMAT<0) then
          do IMD=1,nrans
          aks(ICVS:ICVE,IMD)=0.d0
          aks(IDCS:IDCE,IMD)=0.d0
          enddo
        else
          do 103 IMD=1,nrans
          aks(ICVS:ICVE,IMD)=aks0(IIMAT,IMD)
          aks(IDCS:IDCE,IMD)=aks0(IIMAT,IMD)
!
          if((icaltb==ke.or.
     &        icaltb==ke_low.or.
     &        icaltb==RNG.or.icaltb==CHEN.or.icaltb==KE2S).and.
     &       (IMD.eq.ike(1).or.IMD.eq.ike(2))) then 
            rk=aks0(IIMAT,ike(1))
            re=aks0(IIMAT,ike(2))
            linike=0
            if(abs(rk+re)<SML) linike=1
            if(linike==1)then
!             u1=(u0(IIMAT)**2+v0(IIMAT)**2+w0(IIMAT)**2)*0.5d0
              u1=(u0(IIMAT)**2+v0(IIMAT)**2+w0(IIMAT)**2)**0.5d0
              if(abs(u1)<SML) u1=0.1d0

!              rk=0.1d0*2**1.5d0*u1**2   !I=0.1
!              re=0.09D0**0.75d0*rk**1.5d0/2.0d0 !l=2,
!              rk=0.1d0*1.5d0*u1*u1
              rk=(0.1d0*u1)**2*1.5d0
              re=rk**(1.5d0)   !0.09D0**0.75*kk**1.5/2.0
              aks(ICVS:ICVE,ike(1))=rk
              aks(ICVS:ICVE,ike(2))=re
              aks(IDCS:IDCE,ike(1))=rk
              aks(IDCS:IDCE,ike(2))=re
            endif
          elseif(icaltb==SDES.and.IMD==ides) then
            aks(ICVS:ICVE,IMD)=sthrmu(IMAT)
            aks(IDCS:IDCE,IMD)=sthrmu(IMAT)
          
          endif
  103     enddo
          if(ical_cavi==1) then
            aks(ICVS:ICVE,ivold)=aks(ICVS:ICVE,icavi)
     &                          *rho02(IIMAT)/rho0(IIMAT)
            aks(IDCS:IDCE,ivold)=aks(IDCS:IDCE,icavi)
     &                          *rho02(IIMAT)/rho0(IIMAT)
          endif
          if(ical_FC>0) then
            aks(ICVS:ICVE,ical_s)=aks0(IIMAT,ical_s)
          endif
        endif
      endif
!
      if(pot_scl) then
        do IMD=1,npotn
        POTNAL(ICVS:ICVE,IMD)=potn0(IIMAT,IMD)
        POTNAL(IDCS:IDCE,IMD)=potn0(IIMAT,IMD)
        enddo
      endif
 110  enddo
!----------------------
! --- Euler two phase 
!----------------------
      if(ieul2ph>0) then
        do 120 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)   !  if(IMAT.lt.0) goto 120
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        prs2(ICVS:ICVE)=0.d0
        prs2(IDCS:IDCE)=0.d0
        rho2(ICVS:ICVE)=rho02(IIMAT) 
        rho2(IDCS:IDCE)=rho02(IIMAT)
        if(IMAT.lt.0) then
          rho2(ICVS:ICVE)=rsld(-IMAT)
          rho2(IDCS:IDCE)=rsld(-IMAT)
        endif
        tmp2(ICVS:ICVE)=t02(IIMAT)
        tmp2(IDCS:IDCE)=t02(IIMAT)
        if(IMAT.lt.0) then
          tmp2(ICVS:ICVE)=t0(IIMAT)
          tmp2(IDCS:IDCE)=t0(IIMAT)
        endif
        if(idrdp.eq.comp) then
          prs2(ICVS:ICVE)=p0(IIMAT)
          prs2(IDCS:IDCE)=p0(IIMAT)
        endif
        do ICOM=1,ncomp
        yys2(ICVS:ICVE,ICOM)=y02(IIMAT,ICOM)
        yys2(IDCS:IDCE,ICOM)=y02(IIMAT,ICOM)
        enddo
        vel2(ICVS:ICVE,1)=u02(IIMAT)
        vel2(ICVS:ICVE,2)=v02(IIMAT)
        vel2(ICVS:ICVE,3)=w02(IIMAT)
        vel2(IDCS:IDCE,1)=u02(IIMAT)
        vel2(IDCS:IDCE,2)=v02(IIMAT)
        vel2(IDCS:IDCE,3)=w02(IIMAT)
!
        do IMD=1,nrans
        if(iaph(1)==IMD.or.iaph(2)==IMD) then
          aks(ICVS:ICVE,IMD)=aks0(IIMAT,IMD)
          aks(IDCS:IDCE,IMD)=aks0(IIMAT,IMD)
        endif
        enddo
 120    enddo
      endif
!
      if(ical_vof.eq.1) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(IMAT.gt.0) then
          rho(ICVS:ICVE)=aks0(IIMAT,ivof)*rho0(IIMAT)
     &                  +(1.d0-aks0(IIMAT,ivof))*rho02(IIMAT)
          rho(IDCS:IDCE)=aks0(IIMAT,ivof)*rho0(IIMAT)
     &                  +(1.d0-aks0(IIMAT,ivof))*rho02(IIMAT)
        endif
        enddo
      endif
!
!
      if(ical_MHD.eq.1.or.ical_MHD.eq.2) then
      endif
!
      if(ical_suf==1) then
        SITDEN=0.d0
        MOLFRC=0.d0
        BLKTHK=0.d0
        do nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        kd=kdbcnd(0,nb)
!        if(kd==kdcvd) then
        if(surfreac(nb)==1) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do 500 iph=1,nphase
          icoms=phs_idx(iph-1)+1
          icome=phs_idx(iph)    
          nbx=phs_nbc(iph)
          ntp=phs_typ(iph)
!
!          SITDEN(IBFS:IBFE,iph)=0.d0
!          do ipcom=icoms,icome
!          icom=phs_com(ipcom)
!          MOLFRC(IBFS:IBFE,icom,1)=0.d0
!          enddo
!
          if(ntp==surphase.and.nbx==nb) then
            dum1=0.d0
            do ipcom=icoms,icome
            icom=phs_com(ipcom)
            dum1=dum1+phs_inifrac(ipcom)
            enddo
            dum2=phs_dns(iph)
!
            if(dum1>SML.and.dum2>SML) then
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              do IBFL=IBFS,IBFE
              MOLFRC(IBFL,icom,1)=phs_inifrac(ipcom)/dum1
              enddo
              enddo
              SITDEN(IBFS:IBFE,iph)=phs_dns(iph)   !mol/m^2
            else
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              MOLFRC(IBFS:IBFE,icom,1)=0.d0
              enddo
              SITDEN(IBFS:IBFE,iph)=0.d0
            endif
          elseif(ntp==blkphase.and.nbx==nb) then
            dum1=0.d0
            do ipcom=icoms,icome
            icom=phs_com(ipcom)
            dum1=dum1+phs_inifrac(ipcom)
            enddo
            dum2=phs_thk(iph)
            BLKTHK(IBFS:IBFE,iph)=dum2
!
            if(dum1>SML.and.dum2>SML) then 
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              do IBFL=IBFS,IBFE
              MOLFRC(IBFL,icom,1)=phs_inifrac(ipcom)/dum1
              enddo
              enddo
              do IBFL=IBFS,IBFE
              dum1=0.d0
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              dum1=dum1+
     &             blk_dns(ipcom)*MOLFRC(IBFL,icom,1)*r_wm(icom)
              enddo
              enddo
              SITDEN(IBFS:IBFE,iph)=dum1*phs_thk(iph)   !mol/m^2 ,initial
            else
              SITDEN(IBFS:IBFE,iph)=0.d0
              do ipcom=icoms,icome
              icom=phs_com(ipcom)
              MOLFRC(IBFS:IBFE,icom,1)=0.d0
              enddo
            endif
          endif
 500      enddo
        endif
        enddo
      endif
!
      if(iters.lt.0) goto 2000
!
 1000 continue
!--------------------------------------------------
! --- <2> User subroutine for initial condition
!--------------------------------------------------
      if(iters==0) then
        write(ifll,'(1X,a)') 
     &  'MSG: User initial subroutine [user_ini.f] is used' 
        iters=0    
        times=0.d0    
!
        call nml_init
     &   (MAT_NO,t0,y0,p0,u0,v0,w0,rho0,vsgm0,aks0,potn0,
     &    t02,y02,u02,v02,w02,rho02,
     &    fai_iMHD,A_iMHD,eps_iMHD,eps0_iMHD,mu_iMHD,
     &    mu0_iMHD,sigma_iMHD,sigma0_iMHD,Current_iMHD,OMEGA_iMHD,
     &    ierr1) 
        if( ierr1.ne.0 ) goto 9999 
!
        do 510 IIMAT=1,NMAT 
        IMAT=MAT_NO(IIMAT)     ! if(IMAT.lt.0) goto 510 
        ICVS=MAT_CVEXT(IIMAT-1)+1 
        ICVE=MAT_CVEXT(IIMAT) 
        IDCS=MAT_DCIDX(IIMAT-1)+1 
        IDCE=MAT_DCIDX(IIMAT) 
        prs(ICVS:ICVE)=0.d0 
        pp0(ICVS:ICVE)=0.d0 
        prs(IDCS:IDCE)=0.d0 
        pp0(IDCS:IDCE)=0.d0 
        rho(ICVS:ICVE)=rho0(IIMAT) 
        rho(IDCS:IDCE)=rho0(IIMAT) 
        if(IMAT.lt.0) then
          rho(ICVS:ICVE)=rsld(-IMAT) 
          rho(IDCS:IDCE)=rsld(-IMAT) 
        endif 
        vel(ICVS:ICVE,1)=u0(IIMAT) 
        vel(ICVS:ICVE,2)=v0(IIMAT) 
        vel(ICVS:ICVE,3)=w0(IIMAT) 
        vel(IDCS:IDCE,1)=u0(IIMAT) 
        vel(IDCS:IDCE,2)=v0(IIMAT) 
        vel(IDCS:IDCE,3)=w0(IIMAT) 
        if(idrdp.eq.comp) then
          prs(ICVS:ICVE)=p0(IIMAT) 
          prs(IDCS:IDCE)=p0(IIMAT) 
        elseif(idrdp.eq.mach0.or.idrdp==incomp) then
          pp0(ICVS:ICVE)=p0(IIMAT) 
          pp0(IDCS:IDCE)=p0(IIMAT) 
        endif
        if(min(vsgm0(IIMAT,1),vsgm0(IIMAT,2)).gt.0.d0 ) then
          sgm=vsgm0(IIMAT,1)
          vmx=vsgm0(IIMAT,2)
          do 402 l=1,3
          call utl_boxmlr(r1,r2)
          vel(ICVS:ICVE,l)=vel(ICVS:ICVE,l)+max(-vmx,min(vmx,sgm*r1))
          vel(IDCS:IDCE,l)=vel(ICVS:ICVE,l)+max(-vmx,min(vmx,sgm*r1))
  402     continue
        endif
        tmp(ICVS:ICVE)=t0(IIMAT)
        tmp(IDCS:IDCE)=t0(IIMAT)
        do 400 ICOM=1,ncomp
        yys(ICVS:ICVE,ICOM)=y0(IIMAT,ICOM)
        yys(IDCS:IDCE,ICOM)=y0(IIMAT,ICOM)
  400   continue
        if(rns_scl) then
          do 403 IMD=1,nrans
          aks(ICVS:ICVE,IMD)=aks0(IIMAT,IMD)
          aks(IDCS:IDCE,IMD)=aks0(IIMAT,IMD)
  403     continue
        endif
 510    continue
!
!----------------------
! --- Euler two phase 
!----------------------
!
        if(ieul2ph>0) then
          do 520 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)     ! if(IMAT.lt.0) goto 510
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          prs2(ICVS:ICVE)=0.d0
          prs2(IDCS:IDCE)=0.d0
          rho2(ICVS:ICVE)=rho02(IIMAT)
          rho2(IDCS:IDCE)=rho02(IIMAT)
          if(IMAT.lt.0) then
            rho2(ICVS:ICVE)=rsld(-IMAT)
            rho2(IDCS:IDCE)=rsld(-IMAT)
          endif
          vel2(ICVS:ICVE,1)=u02(IIMAT) 
          vel2(ICVS:ICVE,2)=v02(IIMAT) 
          vel2(ICVS:ICVE,3)=w02(IIMAT) 
          vel2(IDCS:IDCE,1)=u02(IIMAT) 
          vel2(IDCS:IDCE,2)=v02(IIMAT) 
          vel2(IDCS:IDCE,3)=w02(IIMAT) 
          if(idrdp.eq.comp) then
            prs2(ICVS:ICVE)=p0(IIMAT) 
            prs2(IDCS:IDCE)=p0(IIMAT) 
          endif
          tmp2(ICVS:ICVE)=t02(IIMAT)
          tmp2(IDCS:IDCE)=t02(IIMAT)
          if(IMAT.lt.0) then
            tmp2(ICVS:ICVE)=t0(IIMAT)
            tmp2(IDCS:IDCE)=t0(IIMAT)
          endif
          do ICOM=1,ncomp
          yys2(ICVS:ICVE,ICOM)=y02(IIMAT,ICOM)
          yys2(IDCS:IDCE,ICOM)=y02(IIMAT,ICOM)
          enddo
 520      continue
        endif
!---------
! --- VOF
!---------
        if(ical_vof.eq.1) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          if(IMAT.gt.0)  then
            rho(ICVS:ICVE)=aks0(IIMAT,ivof)*rho0(IIMAT)
     &                  +(1.d0-aks0(IIMAT,ivof))*rho02(IIMAT)
            rho(IDCS:IDCE)=aks0(IIMAT,ivof)*rho0(IIMAT)
     &                  +(1.d0-aks0(IIMAT,ivof))*rho02(IIMAT)
          endif
          enddo
        endif
!---------
!--- CAVI
!---------
!        if(ical_cavi==1) then
!          do IIMAT=1,NMAT
!          IMAT=MAT_NO(IIMAT)
!          ICVS=MAT_CVEXT(IIMAT-1)+1
!          ICVE=MAT_CVEXT(IIMAT)
!          IDCS=MAT_DCIDX(IIMAT-1)+1
!          IDCE=MAT_DCIDX(IIMAT)
!          aks(ICVS:ICVE,ivold)=aks(ICVS:ICVE,icavi)
!     &                          *rho02(IIMAT)/rho0(IIMAT)
!          aks(IDCS:IDCE,ivold)=aks(IDCS:IDCE,icavi)
!     &                          *rho02(IIMAT)/rho0(IIMAT)
!          enddo
!        endif
!------------
! --- MHD
!------------
        IF(ical_MHD.eq.1.or.ical_MHD.eq.2) then
        endif
!
!-----------------------------------
! --- call User initial subroutine
!-----------------------------------
!
        if(iniusr==usryes) then 
          do 122 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT.gt.0) then
            IMAT_U=nofld(IMAT)
          elseif(IMAT.lt.0) then
            IMAT_U=nosld(-IMAT)
          endif
          t1=t0(IIMAT)
          y1(1:ncomp)=y0(IIMAT,1:ncomp)
          if(idrdp.eq.comp) then
            p1=p0(IIMAT)
          elseif(idrdp.eq.mach0) then
          else
            p1=0.d0
          endif
          rho1=rho0(IIMAT)
          if(IMAT.lt.0) rho1=rsld(-IMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          do ICVL=ICVS,ICVE
          ICV=MAT_CV(ICVL)
          if(NPE.gt.1) ICV=NOD_IDg(ICV)
          waldis=DISALL(ICVL)
          x=CVCENT(1,ICVL)
          y=CVCENT(2,ICVL)
          z=CVCENT(3,ICVL)
          vol=CVVOLM(ICVL)
          aks1(1:nrans)=aks0(IIMAT,1:nrans)
          u1=u0(IIMAT)
          v1=v0(IIMAT)
          w1=w0(IIMAT)
          p1=p0(IIMAT)
          rho1=rho0(IIMAT)
          y1(:)=y0(IIMAT,:)
          t1=t0(IIMAT)
          aks1(:)=aks0(IIMAT,:)
          potn1(:)=potn0(IIMAT,:)
          call user_ini
     &       (ICV,IMAT_U,ncomp,nrans,npotn,
     &        iters,times,x,y,z,waldis,vol,
     &        u1,v1,w1,p1,rho1,y1,t1,aks1,potn1)
          vel(ICVL,1)=u1
          vel(ICVL,2)=v1
          vel(ICVL,3)=w1
!
          if(idrdp.eq.incomp.or.idrdp.eq.mach0) then
            prs(ICVL)=0.d0
          elseif(idrdp.eq.comp) then
            prs(ICVL)=p1
          endif
          rho(ICVL)=rho1
!
          tmp(ICVL)=t1
          do icom=1,ncomp
          yys(ICVL,icom)=y1(icom)
          enddo
!
          if(rns_scl) then
            do IMD=1,nrans
            aks(ICVL,IMD)=aks1(IMD)
            enddo
          endif
          if(pot_scl) then
            do IMD=1,npotn
            POTNAL(ICVL,IMD)=potn1(IMD)
            enddo
          endif
          enddo
 122      continue
!-------------------------------------
! --- DC initial for user subroutine
!-------------------------------------
         call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rho)
         call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,prs)
         call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,pp0)
         call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,tmp)
         if(rns_scl) then
           call dc_symprs(nrans,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,aks)
         endif
         if(pot_scl) then
           call dc_symprs(npotn,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,POTNAL)
         endif
         call dc_symprs(ncomp,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,yys)
        endif
!------------------
!--- Euler two phase
!------------------
        if(iniusr_E2P==usryes) then
          if(ieul2ph>0) then
            do 130 IIMAT=1,NMAT
            IMAT=MAT_NO(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do 135 ICVL=ICVS,ICVE
            ICV=MAT_CV(ICVL)
            if(IMAT.gt.0) then
              IMAT_U=nofld(IMAT)
            elseif(IMAT.lt.0) then
              IMAT_U=nosld(-IMAT)
            endif
            t2=t02(IIMAT)
            y2(1:ncomp)=y02(IIMAT,1:ncomp)
            dens2=rho02(IIMAT)
            if(IMAT.lt.0) dens2=rsld(-IMAT)
            if(NPE.gt.1) ICV=NOD_IDg(ICV)
            waldis=DISALL(ICVL)
            x=CVCENT(1,ICVL)
            y=CVCENT(2,ICVL)
            z=CVCENT(3,ICVL)
            vol=CVVOLM(ICVL)
            u2=u02(IIMAT)
            v2=v02(IIMAT)
            w2=w02(IIMAT)
            p2=p0(IIMAT)
            alpha(1)=aks(ICVL,iaph(1))
            alpha(2)=aks(ICVL,iaph(2))
            call user_ini_E2P
     &      (ICV,IMAT_U,ncomp,nrans,2,iters,times,x,y,z,waldis,vol,
     &       u2,v2,w2,p2,dens2,y2,t2,alpha)
            vel2(ICVL,1)=u2
            vel2(ICVL,2)=v2
            vel2(ICVL,3)=w2
            if(idrdp.eq.incomp.or.idrdp.eq.mach0) then
              prs2(ICVL)=0.d0
            elseif(idrdp.eq.comp) then
              prs2(ICVL)=p2
            endif
            rho2(ICVL)=dens2
            tmp2(ICVL)=t2
            if(IMAT<0) then
              tmp2(ICVL)=t0(IIMAT)
            endif
            do icom=1,ncomp
            yys2(ICVL,icom)=y2(icom)
            enddo
!
            if(rns_scl) then
              aks(ICVL,iaph(1))=alpha(1)
              aks(ICVL,iaph(2))=alpha(2)
            endif

 135        enddo
 130        enddo

!-------------------------------------
! --- DC initial for user subroutine
!-------------------------------------
            call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rho2)
            call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,tmp2)
            call dc_symprs(ncomp,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,yys2)
            call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,prs2)
          endif
        endif
      endif   ! end user subroutine
!-----------------------------------
! --- MHD user initial Condition
!-----------------------------------
      if((ical_MHD.eq.1.or.ical_MHD.eq.2).and.ini_MHD_user==usryes) then
         do 230 IIMAT=1,NMAT
         IMAT=MAT_NO(IIMAT)
         ICVS=MAT_CVEXT(IIMAT-1)+1
         ICVE=MAT_CVEXT(IIMAT)
         do 235 ICVL=ICVS,ICVE
         ICV=MAT_CV(ICVL)
         if(IMAT.gt.0) then
           IMAT_U=nofld(IMAT)
         elseif(IMAT.lt.0) then
           IMAT_U=nosld(-IMAT)
         endif
         DO I=1,3
         MHD_Ai(I,1)=A_iMHD(IIMAT,I,1)
         MHD_Ai(I,2)=A_iMHD(IIMAT,I,2)
         MHD_CURRENTi(I,1)=Current_iMHD(IIMAT,I,1)
         MHD_CURRENTi(I,2)=Current_iMHD(IIMAT,I,2)
         ENDDO
         MHD_FAIi(1)=fai_iMHD(IIMAT,1)
         MHD_FAIi(2)=fai_iMHD(IIMAT,2)
         MHD_SIGMAi=sigma_iMHD(IIMAT)
         if(NPE.gt.1) ICV=NOD_IDg(ICV)
         waldis=DISALL(ICVL)
         x=CVCENT(1,ICVL)
         y=CVCENT(2,ICVL)
         z=CVCENT(3,ICVL)
         vol=CVVOLM(ICVL)
         aks1(1:nrans)=aks0(IIMAT,1:nrans)
         u1=vel(ICVL,1)
         v1=vel(ICVL,2)
         w1=vel(ICVL,3)
         p1=prs(ICVL)
         rho1=rho(ICVL)
         y1(:)=y0(IIMAT,:)
         t1=t0(IIMAT)
!
!         call usr_ini_MHD
!     &    (ICV,IMAT_U,ncomp,nrans,iters,times,x,y,z,waldis,vol,
!     &     u1,v1,w1,p1,rho1,y1,t1,aks1,
!     &     MHD_Ai,MHD_CURRENTi,MHD_FAIi,MHD_SIGMAi)
         DO I=1,3
         MHD_A(ICVL,I,1,1)=MHD_Ai(I,1)
         MHD_A(ICVL,I,2,1)=MHD_Ai(I,2)
         MHD_CRT0(ICVL,I,1)=MHD_CURRENTi(I,1)
         MHD_CRT0(ICVL,I,2)=MHD_CURRENTi(I,2)
         enddo
         MHD_FAI(ICVL,1,1)=MHD_FAIi(1)
         MHD_FAI(ICVL,2,1)=MHD_FAIi(2)
         SIGMA(ICVL)=MHD_SIGMAi
 235     enddo
 230     enddo
      endif
!
      goto 2000
!-------------------------------------------------------
! --- <3> User defined initial field file or restart file
!-------------------------------------------------------
 3000 continue
!
      if(ifl.lt.0) then
        if(NPE.gt.1) then
          ! get file channel number : ifl
          call getfil(ifl,fnam,'restart')
          fnam=INITin
!fortarn90          ifl=ifl+my_rank
        else
          call getfil(ifl,fnam,'initial')
        endif
      endif
!
      inifil=.false.
      if(fnam.eq.' ') then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of initial file name in control file'
        goto 9999
      elseif(fnam.ne.' ') then
        inquire(file=fnam,exist=inifil)
        if(.not.inifil) then
          write(ifle,*) '### error : file error'
          write(ifle,*) 'lack of initial file ',TRIM(adjustl(fnam))
          goto 9999
        endif
      endif
!
      iter=-999
      nrec=0
  200 continue
!
      do 220 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      IDCS=MAT_DCIDX(IIMAT-1)+1
      IDCE=MAT_DCIDX(IIMAT)
      rho(ICVS:ICVE)=rho0(IIMAT)
      rho(IDCS:IDCE)=rho0(IIMAT)
      tmp(ICVS:ICVE)=t0(IIMAT)
      tmp(IDCS:IDCE)=t0(IIMAT)
 220  continue
!
      if(ical_MHD.eq.1.or.ical_MHD.eq.2) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(IMAT.gt.0.or.IMAT.lt.0)  then
          DO I=1,3
          MHD_CRT0(ICVS:ICVE,I,1)=Current_iMHD(IIMAT,I,1)
          MHD_CRT0(ICVS:ICVE,I,2)=Current_iMHD(IIMAT,I,2)
          ENDDO
          SIGMA(ICVS:ICVE)=sigma_iMHD(IIMAT)
          OMGA(IIMAT)=OMEGA_iMHD(IIMAT)*2.d0*3.14d0
        endif
        enddo
      endif
!
!
      open(ifl,file=fnam,form='unformatted',status='old',iostat=ios)
      read(ifl,iostat=ios) NALLCVX,NCVX,ncmpx,nrnsx,npotnx,NCVFAX,NMATX,
     &                     ncompallx
      call errmsg0
      if( ierr1.ne.0 ) goto 9998
      itera=iter
      read(ifl,iostat=ios) 
     &     iter,ieul2phx,ical_sufx,ical_sldx,ical_MHDx,times,IMVMSHx,
     &     ista_momeryx,ical_WDRx,N_injx
!------------------------------------------------------------
!      read(ifl,iostat=ios) iter,ieul2phx,times
!      ical_sufx=0
!      print*,'Check read_total.f !!!'
!------------------------------------------------------------
!
      if(ical_WDR/=ical_WDRx) then
        write(ifle,'(1x,a)') ' MSG: [ical_WDR] has been changed '
        write(ifle,'(1x,a)') ' MSG: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
!        call FFRABORT(1,'ERR: reading restart file')
      endif
!
      if(N_inj/=N_injx) then
        write(ifle,'(1x,a)') ' MSG: [N_injx] has been changed '
        write(ifle,'(1x,a)') ' MSG: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
!        call FFRABORT(1,'ERR: reading restart file')
      endif
!
      if(ical_sld/=ical_sldx) then
        write(ifle,'(1x,a)') ' ### ERR: [ical_sld] has been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(ical_MHDx/=ical_MHD) then
        write(ifle,'(1x,a)') ' ### ERR: [ical_MHD] has been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if((IMVMSHx/=0.and.ical_mvmsh==0).or.
     &   (IMVMSHx==0.and.ical_mvmsh/=0)) then
        write(ifle,'(1x,a)') ' ### ERR: [Cal_mvmsh] has been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(ieul2phx/=ieul2ph) then
        write(ifle,'(1x,a)') ' ### ERR: [Euler2ph] has been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(ical_suf/=ical_sufx) then
        write(ifle,'(1x,a)') ' ### ERR: Surface reaction been changed'
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(ncmpx/=ncomp) then
        write(ifle,'(1x,a)') ' ### ERR: Species number been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(ncompallx/=ncompall) then
        write(ifle,'(1x,a)') ' ### ERR: Species number been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(nrnsx/=nrans) then
        write(ifle,'(1x,a)') ' ### ERR: Scalar number been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(npotnx/=npotn) then
        write(ifle,'(1x,a)') 
     &  ' ### ERR: Scalar Potential number been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if(NALLCVX/=NALLCV) then
        write(ifle,'(1x,a)') ' ### ERR: Node number been changed '
        write(ifle,'(1x,a)') ' ### ERR: Can NOT restart.'
        write(ifle,'(1x,2a)') 
     &    'MSG: Really read restart file? :',trim(fnam)
        call FFRABORT(1,'reading restart file')
      endif
      if( iter.ne.iters ) then
        write(ifle,'(1x,a)') '### error : data error'
        write(ifle,'(1x,a)') 'specified iteration counts is not found'
        write(ifle,'(1X,a,I10)') 'in initial file step no.= ',iter
        write(ifle,'(1X,a,I10)') 'step no. in control file =:',iters
        write(ifle,'(1x,a)') 'Please reset your control file'
        goto 9998
!
      elseif( iter.eq.iters ) then
      endif
!
      call errmsg0
      call sizchk(ierr1)
!
      if(ical_sldx/=0) read(ifl,iostat=ios) (rot_ang(IMAT),IMAT=1,nflud)
!
      read(ifl,iostat=ios) (pp0(ICVL),ICVL=1,NALLCVX)
      call errmsg0
      read(ifl,iostat=ios) (prs(ICVL),ICVL=1,NALLCVX)
      call errmsg0
      read(ifl,iostat=ios) (hhh(ICVL),ICVL=1,NALLCVX)
      call errmsg0
      read(ifl,iostat=ios) 
     &                ((yys(ICVL,ICOM),ICOM=1,ncmpx),ICVL=1,NALLCVX)
      call errmsg0
      read(ifl,iostat=ios) ((vel(ICVL,l),l=1,3),ICVL=1,NALLCVX)
      call errmsg0
      if(nrnsx.gt.0 ) then
        read(ifl,iostat=ios)  
     &  ((aks(ICVL,IMD),IMD=1,nrnsx),ICVL=1,NALLCVX)
        call errmsg0
      endif
      if(npotnx.gt.0 ) then
        read(ifl,iostat=ios) 
     &  ((POTNAL(ICVL,IMD),IMD=1,npotnx),ICVL=1,NALLCVX)
        call errmsg0
      endif
      if(NCVFAX.gt.0) then
        read(ifl,iostat=ios) (rva(ICFL),ICFL=1,NCVFAX)
        call errmsg0
      endif
!
      read(ifl,iostat=ios) (((dvdd(ICVL,i,k),i=1,3),ICVL=1,NCVX),k=1,2)
      call errmsg0
!
      if(ieul2phx>0) then
        read(ifl,iostat=ios) (tmp2(ICV),ICV=1,NALLCVX)
        call errmsg0
        read(ifl,iostat=ios)  
     &           ((yys2(ICV,ICOM),ICOM=1,ncmpx),ICV=1,NALLCVX)
        call errmsg0
        read(ifl,iostat=ios) ((vel2(ICV,i),i=1,3),ICV=1,NALLCVX)
        call errmsg0
        read(ifl,iostat=ios) (prs2(ICV),ICV=1,NALLCVX)
        call errmsg0
        read(ifl,iostat=ios) (rva2(ICF),ICF=1,NCVFAX)
        call errmsg0
        read(ifl,iostat=ios) (((dvdd2(ICV,i,k),i=1,3),ICV=1,NCVX),k=1,2)
        call errmsg0
        read(ifl,iostat=ios) (hhh2(ICV),ICV=1,NALLCVX)
        call errmsg0
      endif
!
      if(ical_sufx==1) then
        read(ifl,iostat=ios)
     &        ((MOLFRC(IBFL,icom,1),icom=1,ncompallx),IBFL=1,NSSFBC)
        read(ifl,iostat=ios)
     &        ((SITDEN(IBFL,iph),iph=1,nphase),IBFL=1,NSSFBC)
        read(ifl,iostat=ios) 
     &        ((BLKTHK(IBFL,iph),iph=1,nphase),IBFL=1,NSSFBC)
        
        call errmsg0
      endif
!
      if(ical_WDRx==1) then
        read(ifl,iostat=ios) 
     &         ((WDRBND(ICV,I),ICV=1,NCVX),I=1,N_injx)
      endif
!
      if(ical_MHDx==1.or.ical_MHDx==2) then
        DO IMHD=1,2
        read(ifl,iostat=ios) ((MHD_A(ICV,I,IMHD,1),ICV=1,NALLCV),I=1,3)
        read(ifl,iostat=ios)  (MHD_FAI(ICV,IMHD,1),ICV=1,NALLCV)
        read(ifl,iostat=ios) ((MHD_CRNT(ICV,I,IMHD),ICV=1,NALLCV),I=1,3)
        enddo
      endif
!
      if(IMVMSHx/=0) then
        read(ifl,iostat=ios) ((cord(I,ICV),ICV=1,NVRTX),I=1,3)
        read(ifl,iostat=ios) (lacell(ICV),ICV=1,NCELL)
        read(ifl,iostat=ios) ((lvcell(I,ICV),ICV=1,NCELL),I=1,8)
      endif
!
      if(ista_momeryx/=ista_momery.and.iter>=ISTART) then
        call FFRABORT(1,'ERR: ista_momeryx')
      endif
      if(ista_momeryx==1.and.iter>=ISTART) then
        read(ifl,iostat=ios) ((STA_SAV(ICV,i),i=1,n_ave),ICV=1,NCV)
      endif
!
      close(ifl)
!
! --- if(iniusr==usryes.and.iters.gt.0) then
!
      if(iniusr==usryesyes) then
        write(ifll,'(1X,a)') 
     &   'MSG: restart file is used for initial condition'
        write(ifll,'(1X,2a)') 
     &        'MSG: and then, user-subroutine is ',
     &        'used for correcting initial condition'
!
          do 420 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT.gt.0) then
            IMAT_U=nofld(IMAT)
          elseif(IMAT.lt.0) then
            IMAT_U=nosld(-IMAT)
          endif
          t1=t0(IIMAT)
          y1(1:ncomp)=y0(IIMAT,1:ncomp)
          if(idrdp.eq.comp) then
            p1=p0(IIMAT)
          elseif(idrdp.eq.mach0) then
          else
            p1=0.d0
          endif
          rho1=rho0(IIMAT)
          if(IMAT.lt.0) rho1=rsld(-IMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          do ICVL=ICVS,ICVE
          ICV=MAT_CV(ICVL)
          if(NPE.gt.1) ICV=NOD_IDg(ICV)
          waldis=DISALL(ICVL)
          x=CVCENT(1,ICVL)
          y=CVCENT(2,ICVL)
          z=CVCENT(3,ICVL)
          vol=CVVOLM(ICVL)
          if(rns_scl) then
            aks1(1:nrans)=aks(ICVL,1:nrans)
          else
            aks1(1:nrans)=aks0(IIMAT,1:nrans)
          endif
          u1=vel(ICVL,1)
          v1=vel(ICVL,2)
          w1=vel(ICVL,3)
          p1=prs(ICVL)
          rho1=rho(ICVL)
          y1(:)=yys(ICVL,:)
          t1=tmp(ICVL)
          if(pot_scl) then
            potn1(:)=POTNAL(ICVL,:)
          else
            potn1(:)=potn0(IIMAT,:)
          endif
          call user_ini
     &       (ICV,IMAT_U,ncomp,nrans,npotn,
     &        iters,times,x,y,z,waldis,vol,
     &        u1,v1,w1,p1,rho1,y1,t1,aks1,potn1)
          vel(ICVL,1)=u1
          vel(ICVL,2)=v1
          vel(ICVL,3)=w1
!
          if(idrdp.eq.incomp.or.idrdp.eq.mach0) then
            prs(ICVL)=0.d0
          elseif(idrdp.eq.comp) then
            prs(ICVL)=p1
          endif
          rho(ICVL)=rho1
!
          tmp(ICVL)=t1
          do icom=1,ncomp
          yys(ICVL,icom)=y1(icom)
          enddo
!
          if(rns_scl) then
            do IMD=1,nrans
            aks(ICVL,IMD)=aks1(IMD)
            enddo
          endif
          if(pot_scl) then
            do IMD=1,npotn
            POTNAL(ICVL,IMD)=potn1(IMD)
            enddo
          endif
          enddo
 420      continue
        
      endif
!
!
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(3,MXALLCV,NCV,vel)
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,prs)
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,pp0)
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,hhh)
        CALL SOLVER_SEND_RECV(ncmpx,MXALLCV,NCV,yys)
        if( nrnsx.gt.0 ) then
          CALL SOLVER_SEND_RECV(nrnsx,MXALLCV,NCV,aks)
        endif
        if( npotnx.gt.0 ) then
          CALL SOLVER_SEND_RECV(npotnx,MXALLCV,NCV,POTNAL)
        endif
        if(ieul2phx>0) then
          CALL SOLVER_SEND_RECV(3,MXALLCV,NCV,vel2)
          CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,tmp2)
          CALL SOLVER_SEND_RECV(ncmpx,MXALLCV,NCV,yys2)
          CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,hhh2)
          CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,prs2)
        endif
      ENDIF
!
      if( ierror.ne.0 ) goto 9999
!
      if(rns_scl .and. nrnsx.lt.1 ) then
        do 225 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        do 215 IMD=1,nrans
        aks(ICVS:ICVE,IMD)=aks0(IIMAT,IMD)
        aks(IDCS:IDCE,IMD)=aks0(IIMAT,IMD)
 215    continue
 225    continue
      endif
!
!-------------------------
!-< 3. Common procedure >-
!-------------------------
 2000 continue
!      iters=0
!
      if((ical_sld==1.or.ical_sld==2)!.or.ical_sld==4)
     &  .and.itersld.eq.-1) then

        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
!        if(ishaft(IMAT)==1) then
          unit(1)=end(1,IMAT)-begin(1,IMAT)
          unit(2)=end(2,IMAT)-begin(2,IMAT)
          unit(3)=end(3,IMAT)-begin(3,IMAT)
          tha=rot_ang(IMAT)
          CALL rotth(unit,tha,bb)
          do i=1,3
          do j=1,3
          rbb(i,j)=bb(j,i)
          enddo
          enddo
          costh=cos(tha)
          sinth=sin(tha)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          org_x=begin(1,IMAT)*costh-begin(2,IMAT)*sinth
          org_y=begin(1,IMAT)*sinth+begin(2,IMAT)*costh

          do ICVL=ICVS,ICVE
          radi(1)=CVCENT(1,ICVL)-begin(1,IMAT)
          radi(2)=CVCENT(2,ICVL)-begin(2,IMAT)
          radi(3)=CVCENT(3,ICVL)-begin(3,IMAT)
          call AXB_UNIT_C(unit,radi,vr)
          dr=radi(1)*unit(1)+radi(2)*unit(2)+radi(3)*unit(3)
          radi(1)=radi(1)-dr*unit(1)
          radi(2)=radi(2)-dr*unit(2)
          radi(3)=radi(3)-dr*unit(3)
          dr=dsqrt(radi(1)*radi(1)+radi(2)*radi(2)+radi(3)*radi(3))
          v00(1)=dr*rotati(IMAT)*vr(1)
          v00(2)=dr*rotati(IMAT)*vr(2)
          v00(3)=dr*rotati(IMAT)*vr(3)
          v_temp(1)=vel(ICVL,1)
          v_temp(2)=vel(ICVL,2)
          v_temp(3)=vel(ICVL,3)
          vel(ICVL,1)=rbb(1,1)*(v_temp(1))!+v00(1))
     &               +rbb(2,1)*(v_temp(2))!+v00(2))
     &               +rbb(3,1)*(v_temp(3))!+v00(3))
     &               -v00(1)
          vel(ICVL,2)=rbb(1,2)*(v_temp(1))!+v00(1))
     &               +rbb(2,2)*(v_temp(2))!+v00(2))
     &               +rbb(3,2)*(v_temp(3))!+v00(3))
     &               -v00(2)
          vel(ICVL,3)=rbb(1,3)*(v_temp(1))!+v00(1))
     &               +rbb(2,3)*(v_temp(2))!+v00(2))
     &               +rbb(3,3)*(v_temp(3))!+v00(3))
     &               -v00(3)
          enddo
!!        endif
        enddo
      endif
!-------------------------------------
!--< 3.1 mass fraction & pressure >--
!-------------------------------------
      if(iters.lt.1) then
        do 300 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        do ICVL=ICVS,ICVE
        sum=0.d0
        do 301 ICOM=1,ncomp
        yys(ICVL,ICOM)=max(0.d0,yys(ICVL,ICOM))
        sum=sum+yys(ICVL,ICOM)*dble(ACT(ICOM))
  301   continue
        sum=1.d0/(sum+SML)
        do 302 ICOM=1,ncomp
        if(ACT(ICOM)==0) cycle
        yys(ICVL,ICOM)=yys(ICVL,ICOM)*sum
  302   continue
        enddo
!
        do ICVL=IDCS,IDCE
        sum=0.d0
        do ICOM=1,ncomp
        yys(ICVL,ICOM)=max(0.d0,yys(ICVL,ICOM))
        sum=sum+yys(ICVL,ICOM)*dble(ACT(ICOM))
        enddo
        sum=1.d0/(sum+SML)
        do ICOM=1,ncomp
        if(ACT(ICOM)==0) cycle
        yys(ICVL,ICOM)=yys(ICVL,ICOM)*sum
        enddo
        enddo
  300   continue
!
!
!
        if(ieul2ph>0) then
          do IIMAT=1,NMAT
            IMAT=MAT_NO(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            IDCS=MAT_DCIDX(IIMAT-1)+1
            IDCE=MAT_DCIDX(IIMAT)
            do ICVL=ICVS,ICVE
              sum=0.d0
              do ICOM=1,ncomp
                yys(ICVL,ICOM)=max(0.d0,yys(ICVL,ICOM))
                sum=sum+yys(ICVL,ICOM)
              enddo
              sum=1.d0/(sum+SML)
              do ICOM=1,ncomp
                yys(ICVL,ICOM)=yys(ICVL,ICOM)*sum
              enddo
            enddo
            do ICVL=IDCS,IDCE
              sum=0.d0
              do ICOM=1,ncomp
                yys(ICVL,ICOM)=max(0.d0,yys(ICVL,ICOM))
                sum=sum+yys(ICVL,ICOM)
              enddo
              sum=1.d0/(sum+SML)
              do ICOM=1,ncomp
                yys(ICVL,ICOM)=yys(ICVL,ICOM)*sum
              enddo
            enddo
          enddo
        endif
      endif
!
      if(ical_vof.eq.1) then
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          do ICVL=ICVS,ICVE
          aks(ICVL,ivof)=min(max(aks(ICVL,ivof),0.d0),1.d0)
          enddo
          do ICVL=IDCS,IDCE
          aks(ICVL,ivof)=min(max(aks(ICVL,ivof),0.d0),1.d0)
          enddo
        enddo
      endif
!
      if(ical_cavi==1) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        aks(ICVS:ICVE,ivold)=aks(ICVS:ICVE,icavi)
     &                      *rho02(IIMAT)/rho0(IIMAT)
        aks(IDCS:IDCE,ivold)=aks(IDCS:IDCE,icavi)
     &                       *rho02(IIMAT)/rho0(IIMAT)
        enddo
      endif
!
      if( idrdp.eq.comp ) pp0=prs
!
!--< 3.2 velocity >--
!
      if( xredc ) vel(:,1)=0.d0
      if( yredc ) vel(:,2)=0.d0
      if( zredc ) vel(:,3)=0.d0
      if(ieul2ph>0) then
        if( xredc ) vel2(:,1)=0.d0
        if( yredc ) vel2(:,2)=0.d0
        if( zredc ) vel2(:,3)=0.d0
      endif
!
! --- < 3.3 k & epsilon >--
!
      if(rns_scl.and.
     &   (icaltb==KE2S.or.icaltb==ke.or.
     &    icaltb==ke_low.or.icaltb==RNG.or.icaltb==CHEN)
     &  ) then
        if( iters.lt.1 ) then
          do 340 IIMAT=1,NMAT
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          IDCS=MAT_DCIDX(IIMAT-1)+1
          IDCE=MAT_DCIDX(IIMAT)
          do ICVL=ICVS,ICVE
            do 320 IMD=1,nrans
              aks(ICVL,IMD)=max(akslw(IMD),aks(ICVL,IMD))
  320       continue
          enddo
          do ICVL=IDCS,IDCE
            do IMD=1,nrans
              aks(ICVL,IMD)=max(akslw(IMD),aks(ICVL,IMD))
            enddo
          enddo
 340      continue
        endif
      else
      endif
!
      if(iLBF_P==4) then
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          dum1=wiface(ICFL)
          dum2=1.d0-wiface(ICFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          velf(ICFL,1)=dum1*vel(ICVA,1)+dum2*vel(ICVB,1)
          velf(ICFL,2)=dum1*vel(ICVA,2)+dum2*vel(ICVB,2)
          velf(ICFL,3)=dum1*vel(ICVA,3)+dum2*vel(ICVB,3)
          enddo
        enddo
      endif
!
      call debug
      return
!
 9998 continue
      write(ifle,*) 'name of initial file =',fnam(:len_trim(fnam))
 9999 continue
      write(ifle,*) '(read_initial)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine errmsg0
!=================================================
      nrec=nrec+1
      if( ios.eq.0 ) return
      if( ierror.ne.0 ) return
      if( ios.gt.0 ) then
        write(ifle,*) '### error : restart file read error'
      else
        write(ifle,*) '### error : end of file detected unexpectedly'
      endif
      write(ifle,*) 'file name = ',fnam
      write(ifle,*) 'record no.3 =',nrec
      write(ifle,*) 'iostat =',ios
      ierror=1
      end subroutine errmsg0
!=================================================
      subroutine sizchk(ierr)
!=================================================
      integer,intent(out) :: ierr
      ierr=0
      if( NALLCVX.ne.NALLCV) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of cells is not matched'
        write(ifle,*) 'This restart file version NOT used '
        ierr=1
      endif
      if( ncmpx.ne.ncomp ) then
        write(ifle,*) '### error : initial file'
        write(ifle,*) 'no. of chemical species is not matched'
        write(ifle,*) 'in initial file :',ncmpx
        write(ifle,*) 'in control data :',ncomp
        ierr=1
      endif
      if( NCVFAX.gt.0 ) then
        if( NCVFAX.ne.NCVFAC ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of cell faces is not matched'
        write(ifle,*) 'in initial file :',NCVFAX
        write(ifle,*) 'in this case    :',NCVFAC
        ierr=1
        endif
      else
        if( iters.gt.0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'mass flux is not stored in initial file'
        ierr=1
        endif
      endif
      end subroutine sizchk
!=================================================
      subroutine debug
!=================================================
      use module_debug,only : idebug
      if( idebug(1).eq.0 ) return
      end subroutine debug
!
      end subroutine read_initial
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_novertex(nvrtx,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only : lenfnm,ifle,getfil
!
! 1. Read vertex file to check no. of vertices
!
! --- [dummy arguments]
!
      integer,intent(in)  :: nvrtx
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(lenfnm) :: fnam
      integer :: ifl,nrec,ngrd,ios,nvrtz
!
!
!
      ierror=0
!
      call getfil(ifl,fnam,'vertex')
!
      nrec=0
      ngrd=0
  100 continue
      read(ifl,iostat=ios) nvrtz
      if( ios.lt.0 ) goto 101
      call errmsg0
      read(ifl,iostat=ios)
      call errmsg0
      if( ierror.ne.0 ) goto 9999
      ngrd=ngrd+1
      if( nvrtz.ne.nvrtx ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of vertices is not matched'
        write(ifle,*) 'in vertex file  :',nvrtz
        write(ifle,*) 'in control data :',nvrtx
        write(ifle,*) 'grid no. =',ngrd
        write(ifle,*) 'name of vertex file =',fnam(:len_trim(fnam))
        goto 9999
      endif
      goto 100
  101 continue
      if( nrec.lt.1 ) call errmsg0
      if( ierror.ne.0 ) goto 9999
!
      return
!
 9999 continue
      write(ifle,*) '(read_novertex)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine errmsg0
!=================================================
      nrec=nrec+1
      if( ios.eq.0 ) return
      if( ierror.ne.0 ) return
      if( ios.gt.0 ) then
      write(ifle,*) '### error-2 : file read error'
      else
      write(ifle,*) '### error : end of file detected unexpectedly'
      endif
      write(ifle,*) 'file name = ',fnam(:len_trim(fnam))
      write(ifle,*) 'record no.4 =',nrec
      write(ifle,*) 'iostat =',ios
      ierror=1
      end subroutine errmsg0
!
      end subroutine read_novertex
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_source(NCV,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only      : lenfnm,ifle,getfil,ifll
      use module_source,only  : nscnd,strdomain
      USE module_usersub,ONLY : src_uvw,usrno,usryes
!
! 1. Input source file
!
! --- [dummy arguments]
!
      integer,intent(in)  :: NCV
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(lenfnm)   :: fnam
      integer,allocatable :: nosdmn(:)
      integer,allocatable :: icsdmn(:)
      integer,allocatable :: lcsdmn(:)
      integer :: k,n,ifl,ios,ierr1,ierr
      integer :: no,kc,kcs,kce
      integer :: nrec,nsdmn,ncsdmn
      logical :: srcfil
!
      return
!
      ierror=0
!
      if(src_uvw.eq.usryes) then
        write(ifll,*) ' ####    User Subroutine [user/usrsrc.f] ',
     &  'is used for source terms by priority'
        return
      endif
!
      ierr=0
      allocate(nosdmn(NCV),icsdmn(0:NCV),lcsdmn(NCV),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in read_source'
        call FFRABORT(1,'read_source')
      endif
!
! --- 
!
      srcfil=.false.
      call getfil(ifl,fnam,'source')
      if( nscnd.gt.0 .and. fnam.eq.' ' ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of source file name in control file'
        goto 9999
      elseif( nscnd.lt.1 .and. fnam.ne.' ' ) then
        write(ifle,*) '### warning : ',
     & 'The source file ', trim(fnam),' is not used'
      elseif(nscnd.gt.0 .and. fnam.ne.' ' ) then
        inquire(file=fnam,exist=srcfil)
        if(.not.srcfil) then
          write(ifle,*) '### error : file error'
          write(ifle,*) 'lack of source file',TRIM(adjustl(fnam))
          goto 9999
        endif
      elseif( fnam.eq.' ' ) then
        return
      endif
!
!-< 1. Input data >-
!
      nrec=0
      n   =0
      kce =0
      icsdmn(0)=0
      open(ifl,file=fnam,form='unformatted',status='old',iostat=ios)
      rewind ifl
  100 continue
      read(ifl,iostat=ios) no,kc
      if( ios.lt.0 ) goto 101
      call errmsg0
      if( ierror.ne.0 ) goto 9999
      n=n+1
      kcs=kce+1
      kce=kce+kc
      if( n.gt.NCV ) goto 9001
      if( kce.gt.NCV ) goto 9002
      nosdmn(n)=no
      icsdmn(n)=kce
      read(ifl,iostat=ios) (lcsdmn(k),k=kcs,kce)
      call errmsg0
      if( ierror.ne.0 ) goto 9999
      goto 100
  101 continue
      nsdmn =n
      ncsdmn=kce
      close(ifl)
!
!
!-< 2. Check data >-
!
      do 300 n=1,nsdmn
!
! / domain no. /
      do 301 k=n+1,nsdmn
      if( nosdmn(n).eq.nosdmn(k) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'duplicated source domain no.'
        write(ifle,*) 'source domain no. =',n
        write(ifle,*) 'name of source file =',fnam(:len_trim(fnam))
        goto 9999
      endif
  301 continue
!
! / range of cell no. /
      do 302 k=icsdmn(n-1)+1,icsdmn(n)
      if( lcsdmn(k).lt.1 .or. lcsdmn(k).gt.NCV ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'cell no. =',lcsdmn(k)
        write(ifle,*) 'it must be > 0 and <=',NCV
        write(ifle,*) 'source domain no. =',n
        write(ifle,*) 'name of source file =',fnam(:len_trim(fnam))
        goto 9999
      endif
  302 continue
!
  300 continue
!
!
!-< 4. Store data into module >-
!
      call strdomain(ifle,NCV,nsdmn,ncsdmn,
     & nosdmn,icsdmn,lcsdmn,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      call debug
!
      deallocate(nosdmn,icsdmn,lcsdmn)
!
      return
!
 9001 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'no. of source domains exceeds NCV'
      write(ifle,*) 'NCV =',NCV
      write(ifle,*) 'source domain no. =',n
      write(ifle,*) 'name of source file =',fnam(:len_trim(fnam))
      goto 9999
 9002 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'no. of cells in source domains exceeds NCV'
      write(ifle,*) 'NCV =',NCV
      write(ifle,*) 'source domain no. =',n
      write(ifle,*) 'name of source file =',fnam(:len_trim(fnam))
      goto 9999
 9999 continue
      write(ifle,*) '(read_source)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine errmsg0
!=================================================
      nrec=nrec+1
      if( ios.eq.0 ) return
      if( ierror.ne.0 ) return
      if( ios.gt.0 ) then
      write(ifle,*) '### error-3 : file read error'
      else
      write(ifle,*) '### error : end of file detected unexpectedly'
      endif
      write(ifle,*) 'file name = ',fnam(:len_trim(fnam))
      write(ifle,*) 'record no.5 =',nrec
      write(ifle,*) 'iostat =',ios
      ierror=1
      end subroutine errmsg0
!=================================================
      subroutine debug
!=================================================
      use module_debug,only : idebug
      if( idebug(1).eq.0 ) return

      end subroutine debug
!
      end subroutine read_source
!
!
