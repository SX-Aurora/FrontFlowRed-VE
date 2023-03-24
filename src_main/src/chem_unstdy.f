!
!      subroutine chem_admin_unstdy
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine chem_admin_unstdy
     & (iter,deltt,
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,MAT_CFIDX,MAT_CV,
     &  LBC_SSF,LCYCSF,LVEDGE,wiface,wifsld,LCYCOLD,OPPANG,vctr,
     &  SFAREA,SFCENT,CVCENT,CVVOLM,
     &  vel,prs,rho,pp0,aks,tmp,yys,wdot,
     &  SDOT,MOLFRC,MOLCEN,SITDEN,
     &  zzs,cp,cr,FIELD_U,minr,grdc,
     &  ierror)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     zzs    <= hhs  (zzs=mol concentration [mole/m^3], or [mole/m^2])
!       
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle,lenfnm,getfil
      use module_param,only    : yslw
      use module_species,only  : gascns,wm,acpk,acpk2,
     &                           r_wm,sw,spcnam,gascns_1,
     &                           chk_spec => chkncomp
      use module_chemcntl,only : itintg,msizfc,divint,
     &                           abserr,relerr
      use module_chemreac,only : nneq,ireq,vreq,preq,UNIT_SI,
     &                           chk_chem => chkncomp,
     &                           ovrall,elemnt,ebrkup,ovalfr,userreac,
     &                           stick_sf,LH_ER,Bohm,elemarbi,
     &                           BohmE,BohmY,const_R,Butler_Volmer,
     &                           Zeldovich,oval_ebrkup,SSFRRM,
     &                           kind_chem,igas,isurface,ig_iter,
     &                           stp_iter,
     &                           ical_suf
      use module_Euler2ph,only : ieul2ph
      use module_model   ,only : icaltb,ke,ical_reac,ical_surf
      use module_material,only : nofld,nflud,relaxys
      use module_boundary,only : kdbcnd,kdcvd,nbcnd,MAT_BCIDX,
     &                           IDX_SUFRAC,
     &                           phs_idx,phs_dns,
     &                           phs_nbc,phs_typ,phs_snm,phs_nam,
     &                           gasphase,surphase,blkphase,
     &                           LBC_INDEX,phs_com,num_site,surfreac
      use module_les   ,only   : rat_E_S,num_ratout
      use module_metrix,only   : rate_eq=>d2work4,lreac
      use module_particle,only : icoal,IFLAG_COAL
!--------------------------------------
! 1. Integrate chemical reaction term
!--------------------------------------
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter
      real*8 ,intent(in)    :: deltt
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   ( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  ( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CFIDX( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV(    MXALLCV)
      INTEGER,INTENT(IN)    :: LBC_SSF(   MXSSFBC)
      integer,intent(in)    :: LCYCSF (   MXSSFBC)

      REAL*8 ,INTENT(IN)   :: WIFACE(        MXCVFAC)
      real*8 ,intent(in)   :: wifsld   (     MXSSFBC_SLD)
      integer,intent(in)   :: LCYCOLD  (     MXSSFBC_SLD)
      real*8 ,intent(in)   :: OPPANG   (     MXSSFBC_SLD)
      integer,intent(in)   :: vctr     (     MXCV_V,0:MXBND_V)

      INTEGER,INTENT(IN)    :: LVEDGE (2, MXCVFAC)



      REAL*8 ,INTENT(IN)    :: CVCENT (3, MXALLCV)
      REAL*8 ,INTENT(IN)    :: CVVOLM(    MXALLCV)
      REAL*8 ,INTENT(IN)    :: SFAREA (4, MXCVFAC)
      REAL*8 ,INTENT(IN)    :: SFCENT (3, MXCVFAC)
      real*8 ,intent(in)    :: rho     (  MXALLCV)
      real*8 ,intent(in)    :: pp0     (  MXALLCV)
      real*8 ,intent(in)    :: aks     (  MXALLCVR,mxrans)
      real*8 ,intent(in)    :: tmp     (  MXALLCV)
      real*8 ,intent(in)    :: vel     (  MXALLCV,3)
      real*8 ,intent(in)    :: prs     (  MXALLCV)
      real*8 ,intent(in)    :: yys     (  MXALLCV,MXcomp)
      real*8 ,intent(inout) :: wdot    (  MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: zzs     (  MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp      (  MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cr      (  MXALLCV)
      REAL*8 ,INTENT(INOUT) :: SDOT    (  MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(IN)    :: MOLFRC(  MXSSFBC_SUF,MXCOMPALL,2)
      REAL*8 ,INTENT(INOUT) :: MOLCEN(    MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(IN)    :: SITDEN  (  MXSSFBC_SUF,MXPHASE)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
      real*8 ,intent(inout) :: minr    (  MXALLCV)
      real*8 ,intent(inout) :: grdc      (   MXALLCV,3,3)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer :: IIMAT,IMAT,IRC,ICVS,ICVE,ICVL,ICOM,ierr1=0
      integer :: IBFS,IBFE,IBFL,ICFL,IDC,i,idum2,idum
      integer :: nb,kd,nbx,ntp,icoms,icome,iph,ipcom,icomm
      integer :: iphg
      integer,save :: Neqrat,ifl=20,ios=0
      real*8,save  :: dum1
      character(lenfnm),save :: fnam
!
      integer,allocatable  :: ami(:)
!      logical,allocatable  :: lreac(:)
!
!---------------------------------------------------------------------
! --- W2K1
!---------------------------------------------------------------------
!
      neqrat=max(1,num_ratout)
!
      ALLOCATE(ami(0:mxcompall),stat=ierr1)
      ALLOCATE(rate_eq(MXCV,Neqrat),stat=ierr1)
      if(ierr1/=0) then
        call FFRABORT(1,'ERR: ALLOCATE error for rate_eq')
      endif
!---------------------------------------------------------------------
      if(MXSSFBC_SUF>MXALLCV) then
        call FFRABORT(1,'ERR: Dimension error in chem_unstdy.')
      endif
!
      ierror=0
      lreac(:)=.false.
      do IRC=1,nneq
        lreac(IRC)=iter>ig_iter(IRC).and.iter<stp_iter(IRC)
      enddo
!
!--------------------------------------
! --- zzs: mol concentration [mol/m^3]
!--------------------------------------
!
      do 100 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 100
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICOM=1,ncomp
      do ICVL=ICVS,ICVE
      zzs(ICVL,ICOM)=rho(ICVL)*yys(ICVL,ICOM)*r_wm(ICOM)
      enddo
      enddo
  100 enddo
!------------------------------------------------
! --- Surface phase mole concentration [mol/m^2] 
!------------------------------------------------
      MOLCEN(:,:)=0.d0 
      if(ical_surf) then
        do 110 nb=1,nbcnd
!
        ami(:)=0
        do IRC=1,nneq
        if(IDX_SUFRAC(nb,IRC)) then
          iph=1
          icoms=phs_idx(iph-1)+1
          icome=phs_idx(iph)
          do  ipcom=icoms,icome
          icom=phs_com(ipcom)
          if(vreq(ICOM,2,IRC)>SML.or.vreq(ICOM,1,IRC)>SML) then
            ami(icom)=1
          endif
          enddo
        endif
        enddo
!
        IIMAT=MAT_BCIDX(nb,1)
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle
        kd=kdbcnd(0,nb)
        if(surfreac(nb)>0) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do 120 iph=1,nphase
          nbx=phs_nbc(iph)
          ntp=phs_typ(iph)
          icoms=phs_idx(iph-1)+1
          icome=phs_idx(iph)
          if(ntp==gasphase) then
            do ipcom=icoms,icome
            icom=phs_com(ipcom)
            if(ami(icom)==1)then
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              ICVL=LVEDGE(1,ICFL)
              IDC=LVEDGE(2,ICFL)
              MOLCEN(IBFL,icom)=zzs(ICVL,ICOM)
              enddo
            endif
            enddo
          elseif(ntp==surphase.and.nbx==nb) then ! SITE SPECIES
            do ipcom=icoms,icome
            icom=phs_com(ipcom)
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
            ICVL=LVEDGE(1,ICFL)
            IDC=LVEDGE(2,ICFL)                   !MOLCEN=[mol/m^2]
            MOLCEN(IBFL,icom)=MOLFRC(IBFL,icom,1)*SITDEN(IBFL,iph)
     &                     /(num_site(ipcom)+SML)
            enddo
            enddo
          elseif(ntp==blkphase.and.nbx==nb)then
            icoms=phs_idx(iph-1)+1 !icom=phs_idx(iph-1)+1,phs_idx(iph)
            icome=phs_idx(iph)     
            do ipcom=icoms,icome
            icom=phs_com(ipcom)
            do IBFL=IBFS,IBFE
            ICFL=LBC_SSF(IBFL)
!            area=SFAREA(4,ICFL)    
! Bluk density is activity and defined in mol fraction
            MOLCEN(IBFL,icom)=MOLFRC(IBFL,icom,1)   !*SITDEN(IBFL,iph)
            ENDDO
            enddo
          elseif(nbx==nb) then
            call FFRABORT(1,'ERR: NOT clear ERR at chem_admin_unstdy')
          endif
 120      enddo
        endif
 110    enddo
      endif
!-------------------------------------
! --- gas phase reaction
!-------------------------------------
      if(ical_reac) then
        do 200 IRC=1,nneq 
        ami(:)=0
        do icom=1,ncompall
        do i=1,num_ratout
        idum=rat_E_S(i,icom)
        if(idum==IRC.and.kind_chem(IRC)==igas) then
          ami(icom)=i
        endif
        enddo
        enddo
!-------------------------------------
! --- 
!-------------------------------------
!      const_R=gascns
!-------------------------------------
!
        if(ireq(IRC)==oval_ebrkup) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
            endif
          endif
        elseif(ireq(IRC)==ovrall) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              call overall_un_rate(IRC,Neqrat,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,zzs,yys,tmp,pp0,cr,wdot,rate_eq)
            endif
          endif
        elseif(ireq(IRC)==elemnt.or.ireq(IRC)==elemarbi) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              call elementary_un_rate(IRC,Neqrat,nneq,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,zzs,yys,tmp,cr,cp,wdot,minr,rate_eq)
            endif
          endif
        elseif(ireq(IRC)==ebrkup) then 
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              call eddybreakup_un_rate2(IRC,Neqrat,deltt,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,yys,rho,aks,cr,cp,wdot,rate_eq)  !  R
            endif
          endif
        elseif(ireq(IRC)==ovalfr) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              call overall_un_rate(IRC,Neqrat,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,zzs,yys,tmp,pp0,cr,wdot,rate_eq)
            endif
          endif
        elseif(ireq(IRC)==userreac) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              call USER_CHEM_RATE(
     &                       MXMAT,MXCOMP,MXALLCV,NMAT,NCOMP,nflud,
     &                       mxcomp_suf,ncomp_suf,
     &                       nofld,MAT_NO,MAT_CVEXT,MAT_DCIDX,CVCENT,
     &                       IRC,nneq,vreq,kind_chem,spcnam,wm,
     &                       gascns,
     &                       zzs,yys,tmp,cr,cp,wdot)
            endif
          endif

        elseif(ireq(IRC)==SSFRRM) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
            endif
          endif
        elseif(ireq(IRC)==stick_sf) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              write(ifle,'(a,I4)') 'ERR: Reaction no= ',IRC
              call FFRABORT
     &        (1,'ERR: [stick] model is surface reaction model')
            endif
          endif
        elseif(ireq(IRC)==LH_ER) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              write(ifle,'(a,I4)') 'ERR: Reaction no= ',IRC
              call FFRABORT
     &        (1,'ERR: [LH_ER] model is surface reaction model')
            endif
          endif
        elseif(ireq(IRC)==Bohm.or.ireq(IRC)==BohmE.or.ireq(IRC)==BohmY)
     &  then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              write(ifle,'(a,I4)') 'ERR: Reaction no= ',IRC
              call FFRABORT
     &        (1,'ERR: [Bohm] model is surface reaction model')
            endif
          endif
        elseif(ireq(IRC)==Butler_Volmer) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              write(ifle,'(a,I4)') 'ERR: Reaction no= ',IRC
              call FFRABORT
     &        (1,'ERR: [Butler_Volmer] model is surface reaction model')
            endif
          endif
        elseif(ireq(IRC)==Zeldovich) then
          if(lreac(IRC)) then
            if(kind_chem(IRC)==igas) then
              call FFRABORT(1,'Zeldovich_NO is NOT supported now')
            endif
          endif
        endif
 200    enddo
      endif
!------------------------------
! --- surface reaction rate 
!------------------------------
      if(ical_surf) then
      do 300 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle
      kd=kdbcnd(0,nb)
      if(surfreac(nb)>0) then
!      if(kd==kdcvd) then
        do 360 IRC=1,nneq
!
        ami(:)=0
        do icom=1,ncompall   
        do i=1,num_ratout    
        idum=rat_E_S(i,icom) 
        if(idum==IRC.and.
     &     kind_chem(IRC)== isurface 
     &     .and.IDX_SUFRAC(nb,IRC)
     &   ) then
          ami(icom)=i
        endif
        enddo
        enddo
!
        if(IDX_SUFRAC(nb,IRC)) then
          if(ireq(IRC)==ovrall) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
                write(ifle,'(a,I4)') 'ERR: Reaction no= ',IRC
                call FFRABORT(1,'NOT support for overall surface')
              endif
            endif
          elseif(ireq(IRC)==elemnt.or.ireq(IRC)==elemarbi) then
            if(lreac(IRC)) then 
              if(kind_chem(IRC)==isurface) then
              endif
            endif
          elseif(ireq(IRC)==ebrkup) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
                write(ifle,'(a,I4)') 'ERR: Reaction no= ',IRC
                call FFRABORT(1,'NOT support for eddy-breakup surface')
              endif
            endif
          elseif(ireq(IRC)==ovalfr) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
                write(ifle,'(a,I4)') 'ERR: Reaction no= ',IRC
                call FFRABORT(1,'Can NOT use overall-fire for surface')
              endif
            endif
          elseif(ireq(IRC)==userreac) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
              endif
            endif
          elseif(ireq(IRC)==stick_sf) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
              endif
            endif
          elseif(ireq(IRC)==LH_ER) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
              endif
            endif
          elseif(ireq(IRC)==Bohm.or.ireq(IRC)==BohmE.or.
     &                              ireq(IRC)==BohmY) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
              endif
            endif
          elseif(ireq(IRC)==Butler_Volmer) then
            if(lreac(IRC)) then
              if(kind_chem(IRC)==isurface) then
              endif
            endif
          endif
        endif
 360    enddo
      endif
 300  enddo
      endif
!
!
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do icom=1,ncomp
        wdot(ICVS:ICVE,icom)=wdot(ICVS:ICVE,icom)*relaxys(IIMAT)
      enddo
      enddo
      SDOT=relaxys(1)*SDOT  
!----------------------------
! --- save rate intemp file :
!----------------------------
      if(num_ratout>0) then
        call getfil(ifl,fnam,'eqrate_wk')
        if(NPE>1) then
          fnam=trim(suf_wk)
!!!!!!          ifl=20
!!frotran          ifl=ifl+my_rank
        else
          fnam='eqrate_wk'
        endif
        open(ifl,file=fnam,form='unformatted',
     &       status='unknown',iostat=ios)
        if(ios/=0) then
          call FFRABORT(1,'ERR: open eqrate_wk error')
        endif
        write(ifl) rate_eq(1:NCV,1:num_ratout)
        close(ifl)
      endif
!
      deALLOCATE(ami,rate_eq)
!
      return
      end subroutine chem_admin_unstdy
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine overall_un_rate(IRC,Neqrat,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,zzs,yys,tmp,pp0,aak,w,rate_eq)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     zzs : mole  concentration [mol/m^3]
!     yys : mass fraction [%]
!     w(ICV,ICOM)=[kg/m^3/s]
!
      use module_species,only  : gascns,wm,gascns_1,spcnam 
      use module_chemreac,only : vreq,preq,creq,UNIT_SI,const_R 
      use module_dimension 
      use module_hpcutil 
      use module_constant 
      use module_les   ,only   : rat_E_S 
!
! 1. Calculate overall reaction rate
!
! 2. If isw>0, calculate Jacobian matrix of reaction rate
!
      implicit none
!
! --- [dummy arguments] 
!
      integer,intent(in)    :: IRC,Neqrat
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      INTEGER,intent(in)    :: ami(0:ncompall)
      real*8 ,intent(in)    :: zzs(  MXALLCV,MXcomp)
      real*8 ,intent(in)    :: yys(  MXALLCV,MXcomp)
      real*8 ,intent(in)    :: tmp(  MXALLCV)
      real*8 ,intent(in)    :: pp0(  MXALLCV)
      real*8 ,intent(inout) :: aak(  MXALLCV)
      real*8 ,intent(inout) :: w  (  MXALLCV,MXcomp)
      real*8 ,intent(inout) :: rate_eq(MXCV,Neqrat)
!
! --- [local entities]
!
      real*8  :: dum1,dum2,dum3
      integer :: IMODE,ICOM,IMD,ICM,ICV,n,s
      integer :: IIMAT,IMAT,ICVS,ICVE
      integer :: icomm
!
      
!------------------------------------------------------------------
! --- forward rate constant: aak(:)=k=A*T^alpha*p^beta*exp(-E/R/T)
!------------------------------------------------------------------
!
!
!
      do 100 IIMAT=1,NMAT   !ICV=1,nko
      if(.not.mat_cal(IIMAT)) goto 100
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 100
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICV=ICVS,ICVE
      aak(ICV)=preq(1,1,IRC)
     &      *(tmp(ICV)**preq(2,1,IRC))
     &      *(pp0(ICV)**preq(3,1,IRC))
     &      *exp(-preq(4,1,IRC)/(const_R(IRC)*tmp(ICV)))
      enddo
  100 continue
!---------------------------------------------------------------
! --- forward reaction rate: aak(:)=rate=k*[fuel]^Csf*[O2]^Cso2 
!---------------------------------------------------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      do 200 ICOM=1,ncomp
      if(vreq(ICOM,1,IRC).gt.SML) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(creq(ICOM,IRC)<0.d0) then
          do ICV=ICVS,ICVE
          if(zzs(ICV,ICOM)<1.d-9) then
            aak(ICV)=0.d0
          else
! --- [mol/(m3 s)]
            aak(ICV)=aak(ICV)*(zzs(ICV,ICOM)**creq(ICOM,IRC))
          endif
          enddo
        else
          do ICV=ICVS,ICVE
          aak(ICV)=aak(ICV)*(zzs(ICV,ICOM)**creq(ICOM,IRC))   ! [mol/(m3 s)]
          enddo
        endif
      endif
  200 enddo
      enddo
!------------------------- -------------------------------
! --- net value of reaction rate :w(ICV,ICOM)=[kg/m^3/s] 
!------------------------- -------------------------------
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
        do ICOM=1,ncomp
        if(abs(vreq(ICOM,1,IRC))<SML.and.
     &     abs(vreq(ICOM,2,IRC))<SML) cycle
        dum1 = 0.d0
        if(abs(vreq(ICOM,1,IRC))>SML)    ! species in l.h.s?
     &    dum1 = -wm(ICOM)*vreq(ICOM,1,IRC)           ! [kg/mol]
        if(abs(vreq(ICOM,2,IRC))>SML)    ! species in r.h.s?
     &    dum1 = dum1+wm(ICOM)*vreq(ICOM,2,IRC)       ! [kg/mol]
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICV=ICVS,ICVE
! --- [kg/(m3 s)] <- [mol/(m3 s)]
          w(ICV,ICOM)=w(ICV,ICOM)+dum1*aak(ICV) 
        enddo
! --- 
        if(ami(icom)/=0) then
          icomm=ami(icom)
          do ICV=ICVS,ICVE
          rate_eq(ICV,icomm)=dum1*aak(ICV)
          enddo
        endif
        enddo
      enddo
      return
      end subroutine overall_un_rate
!
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine elementary_un_rate(IRC,Neqrat,nneq,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,mc_i,yi,tmp,aak1,aak2,w,minr,rate_eq)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_species,only  : wm,judg_t2,spcnam
      use module_chemreac,only : revised_mc,vreq,preq,M_3rd,const_R,
     &                           P_3rd,UNIT_SI
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_les   ,only   : rat_E_S
! 1. Calculate elementary reaction rate
!
! 2. If isw>0, calculate Jacobian matrix of reaction rate
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IRC,Neqrat,nneq
      INTEGER,INTENT(IN)    :: MAT_CVEXT (0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO    (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (0:MXMAT)
      INTEGER,intent(in)    :: ami       (0:ncompall)
      real*8 ,intent(in)    :: mc_i(MXALLCV,MXcomp) ! mol concentr. [mol/m3]
      real*8 ,intent(in)    :: yi(MXALLCV,MXcomp)   ! mass fraction
      real*8 ,intent(in)    :: tmp      (MXALLCV)
      real*8 ,intent(inout) :: aak1     (MXALLCV)
      real*8 ,intent(inout) :: aak2     (MXALLCV)
      real*8 ,intent(inout) :: w(MXALLCV,MXcomp)    ! net value of reaction
                                                    ! rate [kg/(m3 s)]
      real*8 ,intent(inout) :: rate_eq(MXCV,Neqrat)
                                                    ! rate_eq [mole/(m^3 s)]
      real*8 ,intent(inout) :: minr    (  MXALLCV)
!
! --- [local entities]
!
! aak1(:) ! reaction rate (constant) [mol/m3/s]
! aak2(:) ! rever reaction rate
! 
      real*8  :: dum1,dum2,dum3,third
      real*8 ,parameter :: undef = -huge(1.d0),SML1=1.d-20
      integer :: i,j,k,l,m,ICVL,ICV,ICOM,ICM,s,icomm
      integer :: IIMAT,IMAT,ICVS,ICVE
!
      aak1(:)=0.d0
      aak2(:)=0.d0
      third=1.d0
!
!--------------------------------------------------------
! --- if every branch point of temperature is same value
!--------------------------------------------------------
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        call rate_constant(ICVS,ICVE,IRC)
        call reaction_rate(ICVS,ICVE,IRC)
        call net_value(ICVS,ICVE,IRC)
      enddo
!
      contains
!//////////////////////////////////////////////////////////////////
!
!--------------------------------------------------
! --- forward/backward rate constant
!===================================================
      subroutine rate_constant(ICVS,ICVE,IRC)
!===================================================
      use module_species,only : gascns
      use module_chemreac,only : sigma_vreq,dg_t,lreq,dg_RT,
     &                           kind_prs
      integer,intent(in)    :: ICVS,ICVE,IRC
      real*8 :: RT,dum_Linde,GASR_T,TTT
!
      if(kind_prs(IRC)==0) then
        do ICVL=ICVS,ICVE
        if(UNIT_SI(IRC)==2) then
          RT=tmp(ICVL)   !Oefelein model
        else
          RT=const_R(IRC)*tmp(ICVL)
        endif
        GASR_T=tmp(ICVL)*gascns
        TTT=tmp(ICVL)
!--------------------------------------
! --- calculate forward rate constant
!--------------------------------------
        aak1(ICVL)=preq(1,1,IRC)
     &     *(tmp(ICVL)**preq(2,1,IRC))
     &     *exp(-preq(3,1,IRC)/RT)
        if(abs(preq(1,2,IRC)).lt.SML1) then
!--------------------------------------
! --- NOT consider backward reaction
!--------------------------------------
          aak2(ICVL)=0.d0
        elseif(preq(1,2,IRC)==undef) then
!---------------------------------------------------
! --- calculate backward rate constant from forward
!---------------------------------------------------
          if(judg_t2(1)) then
            dum2=exp(dg_RT(TTT,IRC))
          else
            dum2=exp(dg_t(ncomp,TTT,IRC)/TTT)
          endif
          dum1=(GASR_T/101325.0d0)**sigma_vreq(IRC)
     &       *dum2
          aak2(ICVL)=aak1(ICVL)*dum1
        else
!---------------------------------------------------
! --- if constant values are given for backward
!---------------------------------------------------
          aak2(ICVL)=preq(1,2,IRC)*(TTT**preq(2,2,IRC))
     &       *exp(-preq(3,2,IRC)/RT)
        endif
        enddo
      else
        do ICVL=ICVS,ICVE
        RT=const_R(IRC)*tmp(ICVL)
        GASR_T=tmp(ICVL)*gascns
        TTT=tmp(ICVL)
!--------------------------------------
! --- calculate forward rate constant
!--------------------------------------
        aak1(ICVL)=preq(1,1,IRC)
     &     *(tmp(ICVL)**preq(2,1,IRC))
     &     *exp(-preq(3,1,IRC)/RT)
!--------------------------------------
! --- Lindemann form
!--------------------------------------
        dum_Linde=Lindemann(ICVL,aak1(ICVL),RT,IRC,1)
        aak1(ICVL)=dum_Linde
!
        if(abs(preq(1,2,IRC)).lt.SML1 ) then
!--------------------------------------
! --- NOT consider backward reaction
!--------------------------------------
          aak2(ICVL)=0
        elseif(preq(1,2,IRC)==undef) then
!---------------------------------------------------
! --- calculate backward rate constant from forward 
!---------------------------------------------------
          dum1=(GASR_T/101325.0d0)**sigma_vreq(IRC)
     &       *exp(dg_t(ncomp,TTT,IRC)/TTT)
          aak2(ICVL)=dum_Linde*dum1
        else
!---------------------------------------------------
! --- if constant values are given for backward
!---------------------------------------------------
          aak2(ICVL)=preq(1,2,IRC)*(tmp(ICVL)**preq(2,2,IRC))
     &       *exp(-preq(3,2,IRC)/RT)
        endif
!---------------------------------------------------
! --- Lindemann form
!---------------------------------------------------
!        if(P_3rd(IRC,2)==1) then
!          aak2(ICVL)=Lindemann(ICVL,aak2(ICVL),RT,IRC,2)
!        endif
        enddo
      endif
      end subroutine rate_constant
!
!--------------------------------------------------------
! --- forward/backward reaction rate [mol/m3/s]
!========================================================
      subroutine reaction_rate(ICVS,ICVE,IRC)
!========================================================
      integer,intent(in)    :: ICVS,ICVE,IRC
!
      if(.false.)  then
      minr(:)=1.d10
      do icom=1,ncomp
      if(abs(vreq(icom,1,IRC)).lt.SML1) cycle
      do ICVL=ICVS,ICVE
      dum1=aak1(ICVL)*(mc_i(ICVL,icom)**vreq(icom,1,IRC))
      aak1(ICVL)=dum1
      dum2=mc_i(ICVL,icom)/abs(vreq(icom,1,IRC))
      minr(ICVL)=min(dum1,dum2,minr(ICVL))
      enddo
      enddo

      do ICVL=ICVS,ICVE
      aak1(ICVL)=minr(ICVL)
      enddo

!
      minr(:)=1.d10
      do icom=1,ncomp
      if(abs(vreq(icom,2,IRC)).lt.SML1) cycle
      do ICVL=ICVS,ICVE
      dum1=aak2(ICVL)*(mc_i(ICVL,icom)**vreq(icom,2,IRC))
      aak2(ICVL)=dum1
      dum2=mc_i(ICVL,icom)/abs(vreq(icom,2,IRC))
      minr(ICVL)=min(dum1,dum2,minr(ICVL))
      enddo
      enddo
      do ICVL=ICVS,ICVE
      aak2(ICVL)=minr(ICVL)
      enddo
!
      else
! --- 
      if(UNIT_SI(IRC)==2) then
        do icom=1,ncomp
        if(abs(vreq(icom,1,IRC)).lt.SML1) cycle
        do ICVL=ICVS,ICVE
        aak1(ICVL)=aak1(ICVL)*(1.d-3*mc_i(ICVL,icom))**vreq(icom,1,IRC)
        enddo
        enddo
!
        do icom=1,ncomp
        if(abs(vreq(icom,2,IRC)).lt.SML1) cycle
        do ICVL=ICVS,ICVE
        aak2(ICVL)=aak2(ICVL)*(1.d-3*mc_i(ICVL,icom))**vreq(icom,2,IRC)
        enddo
        enddo
      else
        do icom=1,ncomp
        if(abs(vreq(icom,1,IRC)).lt.SML1) cycle
        do ICVL=ICVS,ICVE
        aak1(ICVL)=aak1(ICVL)*(mc_i(ICVL,icom)**vreq(icom,1,IRC))
        enddo
        enddo
!
        do icom=1,ncomp
        if(abs(vreq(icom,2,IRC)).lt.SML1) cycle
        do ICVL=ICVS,ICVE
        aak2(ICVL)=aak2(ICVL)*(mc_i(ICVL,icom)**vreq(icom,2,IRC))
        enddo
        enddo
      endif
!
      endif
      end subroutine reaction_rate
!
!--------------------------------------------------------
! --- calculate Jacobian matrix of reaction rate
!
!--------------------------------------------------------
! --- net value of reaction rate of species [kg/(m3 s)]
!========================================================
      subroutine net_value(ICVS,ICVE,IRC)
!========================================================
      use module_chemreac,only : sub_vreq,mreq
      integer,intent(in)      :: ICVS,ICVE,IRC
      real*8  :: unit_ef=1.d0
!
      if(UNIT_SI(IRC)==2) then
        unit_ef=1.d3
      else
        unit_ef=1.d0
      endif
!
      if(M_3rd(IRC)==1.and.P_3rd(IRC,1)/=1) then
        if(UNIT_SI(IRC)==2) then
          do ICVL=ICVS,ICVE
          dum2=revised_mc(ncomp,IRC,1.d-3*mc_i(ICVL,1:ncomp))
          aak1(ICVL)=(aak1(ICVL)-aak2(ICVL))*dum2
          enddo
        else
          do ICVL=ICVS,ICVE
          dum2=revised_mc(ncomp,IRC,mc_i(ICVL,1:ncomp))
          aak1(ICVL)=(aak1(ICVL)-aak2(ICVL))*dum2
          enddo
        endif
      else
        do ICVL=ICVS,ICVE
        aak1(ICVL)=(aak1(ICVL)-aak2(ICVL))
        enddo
      endif
!

      do icom=1,ncomp
        if(abs(sub_vreq(icom,IRC)).lt.SML1) cycle
        dum3=wm(icom)*sub_vreq(icom,IRC)*unit_ef
        do ICVL=ICVS,ICVE
!---------------------------------
! --- net value of each equation  
!---------------------------------
          w(ICVL,icom)=w(ICVL,icom)+dum3*aak1(ICVL)
        enddo
!
        if(ami(icom)/=0) then
          dum3=sub_vreq(icom,IRC)
          icomm=ami(icom)
          do ICV=ICVS,ICVE
          rate_eq(ICV,icomm)=dum3*aak1(ICV)
          enddo
        endif
      enddo
      end subroutine net_value
!
! --- forward reaction rate derived from Lindemann form
!=======================================================
      real*8 function Lindemann(ICVL,k_high,RT,q,rev)
!=======================================================
      use module_chemreac,only :
     &    lreq,		! parameters for low pressure
     &    kind_prs,	! flag of Lindemann form
     &    revised_mc,   ! revised mol concentration by 'mreq' [mol/m3]
     &    troe,         ! parameters for Troe form
     &    falloff,activated
      integer,intent(in) ::
     &    ICVL,         ! counter
     &    q,            ! classification number of equation
     &    rev		! switch of reaction direction
      real*8,intent(in) ::
     &    k_high,	! forward reaction rate of high pressure
     &    RT		! (universal gas constant)*(temperature) [J/mol]
      real*8 :: k_low	! forward reaction rate of low pressure
      real*8 :: Pr	! reduced pressure [mol/m3]
      real*8 :: SUM_M
!
!
!
      k_low = lreq(1,rev,q)
     &   *(tmp(ICVL)**lreq(2,rev,q))
     &   *exp(-lreq(3,rev,q)/RT)
      Pr=k_low*revised_mc(ncomp,q,mc_i(ICVL,1:ncomp))/k_high   !+SML
      if(kind_prs(q)==falloff) then	! unimolecular(fall-off)
        Lindemann=k_high*Pr/(1.d0+Pr)
      else			        ! bimolecular(activated)
        Lindemann=k_low/(1.d0+Pr)
      endif
      if(troe(1,rev,q)/=undef) then
        Lindemann=Lindemann*f_Troe(tmp(ICVL),Pr,q,rev)
      endif
      end function Lindemann
!=======================================================
      real*8 function f_Troe(T,Pr,q,rev)
!=======================================================
      use module_chemreac,only :
     &    troe		! parameters for Troe form
      real*8,intent(in) ::
     &    T,		! temperature [k]
     &    Pr		! reduced pressure [mol/m3]
      real*8 :: dumf

      integer,intent(in) ::
     &    q,		! classification number of equation
     &    rev		! switch of reaction direction
      real*8 :: c,d,n
      f_Troe=(1.d0-troe(1,rev,q))*exp(-T/troe(2,rev,q))
     &      +troe(1,rev,q)*exp(-T/troe(3,rev,q))
     &      +exp(-troe(4,rev,q)/T)
      dumf=log10(f_Troe)
      c = -0.4d0-0.67d0*dumf
      d = 0.14d0
      n = 0.75d0-1.27d0*dumf
      f_Troe=dumf/(1.d0+((log10(Pr)+c)/(n-d*(log10(Pr)+c)))**2)
      f_Troe=10.d0**f_Troe
      end function f_Troe
!
      end subroutine elementary_un_rate
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$15762
      subroutine eddybreakup_un_rate(IRC,Neqrat,deltt,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,yys,rho,aks,aak,yi,w,rate_eq)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_species,only  : gascns,wm
      use module_chemreac,only : vreq,preq,const_R
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_scalar,  only : ike
      use module_les   ,only   : rat_E_S
      use module_model,only   : icaltb,ke,RSM,ke_low,RNG,CHEN,KE2S
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IRC,Neqrat
      real*8 ,intent(in)    :: deltt
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      INTEGER,intent(in)    :: ami      (0:ncompall)
      real*8 ,intent(in)    :: yys      (MXALLCV,MXcomp)
      real*8 ,intent(in)    :: rho      (MXALLCV)
      real*8 ,intent(in)    :: aks      (MXALLCVR,mxrans)
      real*8 ,intent(inout) :: aak      (MXALLCV)
      real*8 ,intent(inout) :: yi       (MXALLCV)
      real*8 ,intent(inout) :: w        (MXALLCV,MXcomp)
      real*8 ,intent(inout) :: rate_eq(MXCV,Neqrat)
                                                    ! rate_eq [mole/(m^3 s)]
!
! --- [local entities]
!
      real*8  :: cr1,cr2,r,dd2,fc1,fc2,dum1,rr,ratio
      integer :: ifc,ioc,ICV,ICOM,icomm
      integer :: IIMAT,IMAT,ICVS,ICVE,ICVL
!
      if(icaltb==ke.or.
     &   icaltb==ke_low.or.
     &   icaltb==RNG.or.
     &   icaltb==KE2S.or.
     &   icaltb==CHEN)
     &  then
      else
        call FFRABORT(1,'ERR: eddybreakup model need KE Model')
      endif
!
      ifc=nint(preq(1,1,IRC))
      ioc=nint(preq(2,1,IRC))
      cr1=preq(3,1,IRC)           
      cr2=preq(4,1,IRC)           
      r  =preq(5,1,IRC)           
      r=wm(ioc)*vreq(ioc,1,IRC)/(wm(ifc)*vreq(ifc,1,IRC))
!
      dd2=preq(6,1,IRC)
      fc1=1.d0/r
      fc2=cr2/(1.d0+r)
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      yi(ICVS:ICVE)=1.d20
      do 110 ICOM=1,ncomp
      if(vreq(ICOM,2,IRC).gt.SML) then
        do 111 ICV=ICVS,ICVE
        yi(ICV)=min(yi(ICV),yys(ICV,ICOM))
  111   continue
      endif
  110 continue
!
      do 200 ICV=ICVS,ICVE
      ratio= aks(ICV,ike(2))/(aks(ICV,ike(1))+SML)
      if(ratio.gt.1.d4.or.aks(ICV,ike(1)).gt.1.d5) then
        aak(icv)=0.0
      else if(yi(ICV).gt.SML.and.fc2>0.d0) then
        aak(ICV)=cr1*rho(ICV)*ratio/(wm(ifc)*vreq(ifc,1,IRC))
     &      *min(yys(ICV,ifc),yys(ICV,ioc)*fc1, yi(ICV)*fc2+dd2)
      else
        aak(ICV)=cr1*rho(ICV)*ratio/(wm(ifc)*vreq(ifc,1,IRC))
     &      *min(yys(ICV,ifc),yys(ICV,ioc)*fc1)
      endif
  200 enddo
!
      do 400 ICOM=1,ncomp
      if(vreq(ICOM,2,IRC).ne.vreq(ICOM,1,IRC)) then
        dum1=wm(ICOM)*(vreq(ICOM,2,IRC)-vreq(ICOM,1,IRC))
        w(ICVS:ICVE,ICOM)=w(ICVS:ICVE,ICOM)+dum1*aak(ICVS:ICVE)

        if(ami(icom)/=0) then
          dum1=(vreq(ICOM,2,IRC)-vreq(ICOM,1,IRC))
          icomm=ami(icom)
          rate_eq(ICVS:ICVE,icomm)=dum1*aak(ICVS:ICVE)
        endif
      endif
  400 continue

!
      enddo
!
      return
      end subroutine EddyBreakup_un_Rate

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine eddybreakup_un_rate1(IRC,Neqrat,deltt,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,yys,rho,aks,aak,yi,w,rate_eq)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_species,only  : gascns,wm
      use module_chemreac,only : vreq,preq,const_R
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_scalar,  only : ike
      use module_les   ,only   : rat_E_S
      use module_model,only   : icaltb,ke,RSM,ke_low,RNG,CHEN,KE2S
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IRC,Neqrat
      real*8 ,intent(in)    :: deltt
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      INTEGER,intent(in)    :: ami      (0:ncompall)
      real*8 ,intent(in)    :: yys      (MXALLCV,MXcomp)
      real*8 ,intent(in)    :: rho      (MXALLCV)
      real*8 ,intent(in)    :: aks      (MXALLCVR,mxrans)
      real*8 ,intent(inout) :: aak      (MXALLCV)
      real*8 ,intent(inout) :: yi       (MXALLCV)
      real*8 ,intent(inout) :: w        (MXALLCV,MXcomp)
      real*8 ,intent(inout) :: rate_eq(MXCV,Neqrat)
                                                    ! rate_eq [mole/(m^3 s)]
!
! --- [local entities]
!
      real*8  :: cr1,cr2,r,dd2,fc1,fc2,dum1,rr,ratio
      integer :: ifc,ioc,ICV,ICOM,icomm
      integer :: IIMAT,IMAT,ICVS,ICVE,ICVL
!
      if(icaltb==ke.or.
     &   icaltb==ke_low.or.
     &   icaltb==RNG.or.
     &   icaltb==KE2S.or.
     &   icaltb==CHEN)
     &  then
      else
        call FFRABORT(1,'ERR: eddybreakup model need KE Model')
      endif
!
!
      ifc=nint(preq(1,1,IRC))
      ioc=nint(preq(2,1,IRC))     
      cr1=preq(3,1,IRC)           
      cr2=preq(4,1,IRC)           
      r  =preq(5,1,IRC)           
      r=wm(ioc)*vreq(ioc,1,IRC)/(wm(ifc)*vreq(ifc,1,IRC))
      dd2=preq(6,1,IRC)
      fc1=1.d0/r
      fc2=cr2/(1.d0+r)
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      yi(ICVS:ICVE)=1.d20
      do 110 ICOM=1,ncomp
      if(vreq(ICOM,2,IRC).gt.SML) then
        do 111 ICV=ICVS,ICVE
        yi(ICV)=min(yi(ICV),yys(ICV,ICOM))
  111   continue
      endif
  110 continue
!

      do 200 ICV=ICVS,ICVE
      if(abs(yi(ICV))<SML) then
        aak(ICV)=cr1*rho(ICV)*aks(ICV,ike(2))
     &       /(aks(ICV,ike(1))*wm(ifc)*vreq(ifc,1,IRC))
     &      *min(yys(ICV,ifc),
     &           yys(ICV,ioc)*fc1)
      else
        aak(ICV)=cr1*rho(ICV)*aks(ICV,ike(2))
     &       /(aks(ICV,ike(1))*wm(ifc)*vreq(ifc,1,IRC))
     &      *min(yys(ICV,ifc),
     &           yys(ICV,ioc)*fc1,
     &           yi(ICV)*fc2+dd2 ) 
      endif
  200 continue
!!
!
      do 400 ICOM=1,ncomp
      if(vreq(ICOM,2,IRC).ne.vreq(ICOM,1,IRC)) then
        dum1=wm(ICOM)*(vreq(ICOM,2,IRC)-vreq(ICOM,1,IRC))
        w(ICVS:ICVE,ICOM)=w(ICVS:ICVE,ICOM)+dum1*aak(ICVS:ICVE)
!
        if(ami(icom)/=0) then
          dum1=(vreq(ICOM,2,IRC)-vreq(ICOM,1,IRC))
          icomm=ami(icom)
          rate_eq(ICVS:ICVE,icomm)=dum1*aak(ICVS:ICVE)
        endif
!
      endif
  400 continue
!
      enddo
!
      return
      end subroutine EddyBreakup_un_Rate1


!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine eddybreakup_un_rate2(IRC,Neqrat,deltt,
     &                       MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &                       ami,yys,rho,aks,aak,yi,w,rate_eq)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_species,only  : gascns,wm
      use module_chemreac,only : vreq,preq,const_R
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_scalar,  only : ike
      use module_les   ,only   : rat_E_S
      use module_model,only   : icaltb,ke,RSM,ke_low,RNG,CHEN,KE2S
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IRC,Neqrat
      real*8 ,intent(in)    :: deltt
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  (0:MXMAT)
      INTEGER,intent(in)    :: ami      (0:ncompall)
      real*8 ,intent(in)    :: yys      (MXALLCV,MXcomp)
      real*8 ,intent(in)    :: rho      (MXALLCV)
      real*8 ,intent(in)    :: aks      (MXALLCVR,mxrans)
      real*8 ,intent(inout) :: aak      (MXALLCV)
      real*8 ,intent(inout) :: yi       (MXALLCV)
      real*8 ,intent(inout) :: w        (MXALLCV,MXcomp)
      real*8 ,intent(inout) :: rate_eq(MXCV,Neqrat)
                                                    ! rate_eq [mole/(m^3 s)]
!
! --- [local entities]
!
      real*8  :: cr1,cr2,r,dd2,fc1,fc2,dum1,rr,ratio
      integer :: ifc,ioc,ICV,ICOM,icomm
      integer :: IIMAT,IMAT,ICVS,ICVE,ICVL
!
      if(icaltb==ke.or.
     &   icaltb==ke_low.or.
     &   icaltb==RNG.or.
     &   icaltb==KE2S.or.
     &   icaltb==CHEN)
     &  then
      else
        call FFRABORT(1,'ERR: eddybreakup model need KE Model')
      endif
!
!
      ifc=nint(preq(1,1,IRC))
      ioc=nint(preq(2,1,IRC))     
      cr1=preq(3,1,IRC)           
      cr2=preq(4,1,IRC)           
      r  =preq(5,1,IRC)           
      r=wm(ioc)*vreq(ioc,1,IRC)/(wm(ifc)*vreq(ifc,1,IRC))
      dd2=preq(6,1,IRC)
      fc1=1.d0/r
      fc2=cr2/(1.d0+r)
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      yi(ICVS:ICVE)=1.d20
      do 110 ICOM=1,ncomp
      if(vreq(ICOM,2,IRC).gt.SML) then
        do 111 ICV=ICVS,ICVE
        yi(ICV)=min(yi(ICV),yys(ICV,ICOM))
  111   continue
      endif
  110 continue
!

      do 200 ICV=ICVS,ICVE
      ratio=min(aks(ICV,ike(2))/
     &    (aks(ICV,ike(1))+SML),1.d0/deltt)*rho(ICV)
      if(abs(yi(ICV))<SML) then
        aak(ICV)=cr1*ratio/(wm(ifc)*vreq(ifc,1,IRC))
     &      *min(yys(ICV,ifc),yys(ICV,ioc)*fc1)
      else
        aak(ICV)=cr1*ratio/(wm(ifc)*vreq(ifc,1,IRC))
     &      *min(yys(ICV,ifc),yys(ICV,ioc)*fc1, yi(ICV)*fc2+dd2)
      endif
  200 continue
!
      do 400 ICOM=1,ncomp
      if(vreq(ICOM,2,IRC).ne.vreq(ICOM,1,IRC)) then
        dum1=(wm(ICOM)*(vreq(ICOM,2,IRC)-vreq(ICOM,1,IRC)))
        w(ICVS:ICVE,ICOM)=w(ICVS:ICVE,ICOM)+dum1*aak(ICVS:ICVE)
!
        if(ami(icom)/=0) then
          dum1=(vreq(ICOM,2,IRC)-vreq(ICOM,1,IRC))
          icomm=ami(icom)
          rate_eq(ICVS:ICVE,icomm)=dum1*aak(ICVS:ICVE)
        endif
!
      endif
  400 continue
!
      enddo
!
      return
      end subroutine EddyBreakup_un_Rate2


!
!
