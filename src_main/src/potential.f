!
!     subroutine potential_admin 
!     subroutine bc_potential
!     subroutine bc_poten_flux
!     subroutine potential_src
!     subroutine dc_sym_potn
!     subroutine bc_setbnd_poten
!     subroutine dc_symprv_potn
!     subroutine J_SRC_PEFC_part
!     subroutine ramd
!     subroutine conv_term_FC
!    subroutine FC_src
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine potential_admin(deltt,iter,iterPOTEN,time,ismpl,nsmpl,
     &   ipfix,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,FRSTCV,disall,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  wifsld,LCYCOLD,FIELD_U,
     &  POTNAL,POFLUX,POTNBC,
     &  kdbp,dp,diag,sigma,rho,rho2,iptfix,aks,
     &  grdc,vctr,rva,
     &  yys,tmp,prs,OPPANG,
     &  iter_POTN,reps_POTN,aeps_POTN,err_POTN
     &  )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-------------------------------------------------
!     dp    <=   RMX(:) 
!-------------------------------------------------
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle
      use module_material,only : nflud,nofld,nosld,porosty
      use module_scalar  ,ONLY : ip_pefc,potname,
     &                           ical_FC,ical_s
      use module_model,   only : PEFC,PAFC,MCFC,SOFC,AFC,MDFC
      use module_initial,only  : elecCondFC,ionCondFC,potn_c
      use module_initial ,only : Ke_FC,Kc_FC,Pv_FC,
     &                           sigma_FC,contang_FC
      use module_boundary,only : kpdirc,kpneum,kpintr,kdpotn,vapotn,
     &                           kpindr
      use module_boundary,only : nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX
      use module_metrix  ,only : diagp => d2work2
      use module_metrix  ,only : srcp  => d2work1
      use module_metrix,only   : rcomp,eva_comp
      use module_species,only  : wm,r_wm,gascns
      use module_metrix,only   : J_SRC
      use module_FUEL  ,only   : vap_no,No_Mem,No_AGDL,No_CGDL,
     &                           Tau_diff,OPEN_V
      
!
      implicit none
!
! --- [dummy arguments] 
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter,ismpl,nsmpl,iterPOTEN
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: FRSTCV    (   MXssfbc)
      real*8 ,intent(in)    :: DISALL    (   MXCV)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      REAL*8,INTENT(INOUT)  :: POTNAL    (MXALLCVP,MXPOTN,2)
      REAL*8,INTENT(INOUT)  :: POFLUX    (MXCVFACP,MXPOTN)
      REAL*8,INTENT(INOUT)  :: POTNBC    (MXSSFBCP,MXPOTN)
      INTEGER,INTENT(INOUT) :: kdbp      (   MXCVFAC)
      real*8 ,intent(inout) :: dp        (   MXALLCV)
      real*8 ,intent(inout) :: diag      (   MXALLCV)
      real*8 ,intent(inout) :: sigma     (   MXALLCV)
      real*8 ,intent(in)    :: rho       (   MXALLCV)
      real*8 ,intent(in)    :: rho2      (   MXALLCVC)
      integer,intent(inout) :: iptfix     (   MXMAT)
      integer,intent(in)    :: LCYCOLD   (   MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld    (   MXSSFBC_SLD)
      real*8 ,intent(in)    :: aks       (   MXALLCVR,mxrans)
      real*8 ,intent(inout) :: grdc      (   MXALLCV,3)
!
      real*8 ,intent(out)   :: err_POTN(MXPOTN)
      real*8 ,intent(out)   :: reps_POTN(MXPOTN),aeps_POTN(MXPOTN)
      integer,intent(out)   :: iter_POTN(MXPOTN)
      real*8 ,intent(inout) :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
      real*8 ,intent(in)    :: yys(MXALLCV,MXcomp)
      real*8 ,intent(in)    :: tmp(MXALLCV)
      real*8 ,intent(in)    :: prs(MXALLCV)
      real*8 ,intent(in)    :: rva(MXCVFAC)
      real*8 ,intent(in) :: OPPANG   (       MXSSFBC_SLD)
      integer,intent(in) :: ipfix     (    MXMAT)
!
! --- [local entities]
!
      real*8  :: dum1,dum2,dum3,dum4,dum5,errab_POTN(MXPOTN),dumA,dumB
      real*8  :: dx,dy,dz,grx,gry,grz,wi1,wi2,err_loc
      integer :: ierr1=0,ierror=0,IREAC=1
      real*8  :: ru,rv,rw,dl,gf1,gf2,relax_pt
      integer :: ICVA,ICVB,IIMAT,IMAT,ICFL,ICFS,ICFE,IMAT_U,ICOM
      integer :: ICVS,ICVE,ICVL,IDCS,IDCE
      integer :: IMODE,IPOT,iter_p,IBFS,IBFE,IBFL,nb,IDC,kdp,ISOL,ISING
      character :: sn
      real*8,parameter :: eps=1.d-10
      integer,parameter:: nter_potn=100
      integer,save :: iclose=0
!
      real*8  :: Xw,aaa,ppp,TTT
!
      ALLOCATE (srcp (MXALLCVP,MXPOTN),stat=ierr1)
      ALLOCATE (diagp(MXALLCVP,MXPOTN),stat=ierr1)
!
      if(ierr1/=0) 
     &    call FFRABORT(1,'allocating error in potential_admin')
!
      relax_pt=deltt
      srcp(:,:) =0.d0
      diagp(:,:)=0.d0
      dp(:)     =0.d0
      sigma(:)  =0.d0
      diag(:)   =0.d0
!
!-------------------
! --- set BC 
!-------------------
!
      call bc_setbnd_poten(iter,time,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &   POTNBC
     &       )
! 
      if(ical_FC==PEFC) then
        IREAC=1+iterPOTEN
        do IPOT=1,Npotn
        if(IPOT==ip_pefc(1).or.IPOT==ip_pefc(2)) then
          if(.true.) then
            call J_SRC_PEFC_part(iter,IREAC,deltt,time,IPOT,
     &        LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &        MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &        yys,tmp,prs,rho,POTNAL,aks,srcp,diagp,J_SRC,FIELD_U)
          endif
        endif
        enddo
      endif
!
!-------------------
! --- BC condition
!-------------------
!
      do 1000 IPOT=1,Npotn
!
      call bc_potential
     &  (IPOT,iclose,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &   SFCENT,CVCENT,SFAREA,
     &   POFLUX,POTNBC,POTNAL(:,:,1),tmp)
!
!----------------------
! --- source term
!----------------------
!
      call potential_src(iter,deltt,time,IPOT,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  srcp,diagp)
!
!--------------------------
! --- 
!--------------------------
!
      dp(:)=POTNAL(:,IPOT,1)
!
      call grad_cell(1,97,
     &   MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &   LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,dp,
     &   grdc)
!
      call dc_symprv_potn
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    SFAREA,SFCENT,grdc,0)
!
!--------------------------
! --- 
!--------------------------
!
      call PEFC_sigma(IPOT,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,aks,prs,
     &  yys,tmp,sigma,FIELD_U
     &   )
!
       call dc_sigma_potn(IPOT,1,LVEDGE,LBC_SSF,LCYCSF,
     &                 mat_cal,MAT_NO,sigma)
!
      POFLUX(:,IPOT)=0.d0
      do 100 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dumA=sigma(ICVA)
        dumB=sigma(ICVB)
!        dum2=1.d0/(wi1/dumA+wi2/dumB)
        dum2=(wi1*dumA+wi2*dumB)
        dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl
        gf1=(dp(ICVB)-dp(ICVA))/dl
        ru=(wi1*grdc(ICVA,1)+wi2*grdc(ICVB,1))
        rv=(wi1*grdc(ICVA,2)+wi2*grdc(ICVB,2))
        rw=(wi1*grdc(ICVA,3)+wi2*grdc(ICVB,3))
        gf2=ru*(SFAREA(1,ICFL)-dx)
     &     +rv*(SFAREA(2,ICFL)-dy)
     &     +rw*(SFAREA(3,ICFL)-dz)
! --- 
        POFLUX(ICFL,IPOT)=SFAREA(4,ICFL)*(gf1-gf2)*dum2
      enddo
 100  enddo
!

      call bc_poten_flux(IPOT,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &   SFCENT,SFAREA,POTNBC,POFLUX,sigma,dp,CVCENT,kdbp)
!
      dp(:)=0.d0
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle 
      IMAT=MAT_NO(IIMAT)
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
      dum1=POFLUX(ICFL,IPOT)
      dp(ICVA)=dp(ICVA)+dum1
      dp(ICVB)=dp(ICVB)-dum1
      enddo
      enddo
!      
!------------------------
! --- 
!------------------------
!
      dum1=1.d0
      if(ical_FC==PEFC) dum1=1.d0  !important
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle 
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      diag(ICVL)=diagp(ICVL,IPOT)*CVVOLM(ICVL)
      dp(ICVL)=dp(ICVL)+dum1*srcp(ICVL,IPOT)*CVVOLM(ICVL)   !important
      enddo
      enddo
!
      if(.false.) then
      if(ical_FC==PEFC.and.IPOT==ip_pefc(2).and..false.) then
!      if(ical_FC==PEFC.and.IPOT==ip_pefc(2)) then
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)   !????
          if(IMAT.gt.0) then 
            IMAT_U=nofld(IMAT)
          elseif(IMAT.lt.0) then
            IMAT_U=nosld(-IMAT)
          endif
!          iptfix(IIMAT)=0
          if(IMAT_U==No_Mem) then 
            dum1=0.d0
            dum2=0.d0
            do ICVL=ICVS,ICVE
            dum1=dum1+dp(ICVL)
            dum2=dum2+CVVOLM(ICVL)
            enddo
            dum1=dum1/dum2
            do ICVL=ICVS,ICVE
            dp(ICVL)=dp(ICVL)-dum1*CVVOLM(ICVL)
            enddo
          else
            dp(ICVS:ICVE)=0.d0
          endif
        enddo
      elseif(ical_FC==PEFC.and.IPOT==ip_pefc(1)) then 
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)   !????
          if(IMAT.gt.0) then
            IMAT_U=nofld(IMAT)
          elseif(IMAT.lt.0) then
            IMAT_U=nosld(-IMAT)
          endif
          if(IMAT_U==No_Mem.or.
     &       IMAT_U==No_AGDL.or.
     &       IMAT_U==No_CGDL.or.
     &       IMAT_U<0) then
          else
            dp(ICVS:ICVE)=0.d0
          endif
        enddo
      endif
      endif
!
      IMODE=7
      iter_POTN(IPOT)=0
      ISOL=1     !=1
      if(ISOL==0) then
        call solve_poisson_unsymm_e2p
!     &    (iter,nflud,deltt,time,
     &    (iter,nflud,1.D0,time,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    iptfix,kdbp,diag,dp,sigma,rho2,aks,rva,tmp,
     &    LCYCOLD,wifsld,
     &    iter_POTN(IPOT),reps_POTN(IPOT),aeps_POTN(IPOT),IMODE,ierror)
      else
        ISING=0
        if(IPOT==ip_pefc(2).and.ical_FC==PEFC) ISING=2
        if(IPOT==ip_pefc(1).and.ical_FC==PEFC) ISING=1  !   1
        call solve_poisson_e2p
!     &    (iter,nflud,deltt,time,ISING,
     &    (iter,nflud,1.D0,time,ISING,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    iptfix,kdbp,diag,dp,sigma,rho2,aks,
     &    LCYCOLD,wifsld,
     &    iter_POTN(IPOT),reps_POTN(IPOT),aeps_POTN(IPOT),IMODE,ierror)

      endif
!
      POTNAL(1:NCV,IPOT,1)=POTNAL(1:NCV,IPOT,1)
     &                    +dp(1:NCV)
!
      err_POTN(IPOT)=0.d0
      errab_POTN(IPOT)=0.d0
      Do 200 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      err_POTN(IPOT)=err_POTN(IPOT)+(dp(ICVL)*dp(ICVL))
      errab_POTN(IPOT)=errab_POTN(IPOT)+(POTNAL(ICVL,IPOT,1)**2)
      enddo
  200 enddo
!
      if(NPE.gt.1) then
        call hpcrsum(err_POTN(IPOT))
        call hpcrsum(errab_POTN(IPOT))
      endif
!
 1000 enddo
!
      deALLOCATE(srcp,diagp,stat=ierr1)
      return
!--------------------------------------------------
      end subroutine potential_admin
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_potential
     &  (IPOT,iclose,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &   SFCENT,CVCENT,SFAREA,
     &   POFLUX,POTNBC,POTNAL,tmp)
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
     &     rotsld,idis,LBC_pair,kdbuff,mem_cat
      use module_boundary,only : kpdirc,kpneum,kpintr,kdpotn,vapotn,
     &                           kpindr
      use module_scalar  ,ONLY : ip_pefc,ical_FC
      use module_model,   only : PEFC,PAFC,MCFC,SOFC,AFC,MDFC
      use module_material,only : nflud,nofld
      use module_FUEL    ,only : No_Mem
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: iclose
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC),IPOT
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: MAT_NO    (   0:MXMAT)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      REAL*8,INTENT(INOUT)  :: POTNAL    (MXALLCVP,MXPOTN)
      REAL*8,INTENT(IN)     :: POFLUX    (MXCVFACP,MXPOTN)
      REAL*8,INTENT(IN)     :: POTNBC    (MXSSFBCP,MXPOTN)
      real*8 ,intent(in)    :: tmp(MXALLCV)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMD,IMAT2,IIMAT2,IMAT_U,IMAT_U2
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ICFM,ICFS,ISEC
      integer :: IIPOT
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      REAL*8  :: xxx,yyy,zzz,dum1,dum2,TTT,OPEN_V1
!
!
!
      iclose=0
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT=MAT_NO(IIMAT)
      IMAT2=MAT_NO(IIMAT2)
      IMAT_U=nofld(abs(IMAT))
      IMAT_U2=nofld(abs(IMAT2))
      kdp=kdpotn(IPOT,nb)
      if(.not.mat_cal(IIMAT)) cycle
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
!
      if(kdp.eq.kpdirc) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        POTNAL(IDC,IPOT)=POTNBC(IBFL,IPOT) 
        enddo
        if(IIMAT/=IIMAT2) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          POTNAL(IDCP,IPOT)=POTNBC(IBFL,IPOT) 
          enddo
        endif
      elseif(kdp.eq.kpneum) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        dum1=POTNBC(IBFL,IPOT)
        xxx=SFCENT(1,ICFL)-CVCENT(1,ICV)
        yyy=SFCENT(2,ICFL)-CVCENT(2,ICV)
        zzz=SFCENT(3,ICFL)-CVCENT(3,ICV)
        dum2=xxx*SFAREA(1,ICFL)
     &      +yyy*SFAREA(2,ICFL)
     &      +zzz*SFAREA(3,ICFL)
        POTNAL(IDC,IPOT)=POTNAL(ICV,IPOT)+dum1*dum2
        enddo
        if(IIMAT/=IIMAT2) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          dum1=POTNBC(IBFL,IPOT)
          xxx=SFCENT(1,ICFP)-CVCENT(1,ICVP)
          yyy=SFCENT(2,ICFP)-CVCENT(2,ICVP)
          zzz=SFCENT(3,ICFP)-CVCENT(3,ICVP)
          dum2=xxx*SFAREA(1,ICFP)
     &        +yyy*SFAREA(2,ICFP)
     &        +zzz*SFAREA(3,ICFP)
          POTNAL(IDCP,IPOT)=POTNAL(ICVP,IPOT)+dum1*dum2
          enddo
        endif
      elseif(kdp.eq.kpintr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        POTNAL(IDC,IPOT)=POTNAL(ICVP,IPOT)
        POTNAL(IDCP,IPOT)=POTNAL(ICV,IPOT)
        enddo
      elseif(kdp.eq.kpindr) then 
        if(INT(vapotn(IPOT,nb))>0) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        IIPOT=INT(abs(vapotn(IPOT,nb)))
        POTNAL(IDC,IPOT)=POTNAL(ICVP,IIPOT)
        POTNAL(IDCP,IPOT)=POTNAL(ICV,IIPOT)
        enddo
        else
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        POTNAL(IDC,IPOT)=POTNAL(ICV,IPOT)
        POTNAL(IDCP,IPOT)=POTNAL(ICVP,IPOT)
        enddo
        endif
      endif
 1000 enddo
!
      return
!
      end subroutine bc_potential
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_poten_flux
     & (IPOT,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &   SFCENT,SFAREA,POTNBC,POFLUX,sigma,dp,CVCENT,kdbp)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_boundary,only : kdilet,kdolet,kdtchi,kdintr,
     &                           kdsymm,kdprdc,kkdirc,kkneum,
     &                           kdstag,kdpres,kdsld,mem_cat,
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,
     &     rotsld,idis,LBC_pair,kdbuff
      use module_boundary,only : kpdirc,kpneum,kpintr,kdpotn,vapotn,
     &                           kpindr
      use module_model,   only : PEFC,PAFC,MCFC,SOFC,AFC,MDFC
      use module_scalar  ,ONLY : ip_pefc,ical_FC
      use module_FUEL    ,only : No_Mem
      use module_material,only : nflud,nofld
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MAT_NO    (   0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC),IPOT
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      REAL*8,INTENT(IN)     :: POTNBC    (MXSSFBCP,MXPOTN)
      REAL*8,INTENT(INOUT)  :: POFLUX    (MXCVFACP,MXPOTN)
      real*8 ,intent(in) :: sigma (MXALLCV)
      integer,intent(out)   :: kdbp      (MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in) :: dp (MXALLCV)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMD,IMAT2,IIMAT2,IMAT_U,IMAT_U2
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,ICFM,ICFS,ISEC,ICVA,ICVB
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      REAL*8  :: xxx,yyy,zzz,dum1,dum2,dumA,dumB,dx,dy,dz,dl
      integer :: IIPOT
!
      kdbp(:)=0
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT=MAT_NO(IIMAT)
      IMAT2=MAT_NO(IIMAT2)
      IMAT_U=nofld(abs(IMAT))
      IMAT_U2=nofld(abs(IMAT2))
      ISEC=0
      if(IMAT_U==No_Mem) ISEC=1
      if(.not.mat_cal(IIMAT)) cycle
      kdp=kdpotn(IPOT,nb)
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
!
!
!
      if(kdp.eq.kpdirc) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbp(ICFL)=1
        enddo
        if(IIMAT/=IIMAT2) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          kdbp(ICFP)=1
          enddo
        endif
      elseif(kdp.eq.kpneum) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        kdbp(ICFL)=2     ! zhang
        enddo
        if(IIMAT/=IIMAT2) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          kdbp(ICFP)=2
          enddo
        endif
      elseif(kdp.eq.kpintr) then !???? 
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        kdbp(ICFL)=0
        kdbp(ICFP)=99   !0 important
        enddo
!      elseif(kdp.eq.kpindr) then
!        if(INT(vapotn(IPOT,nb))>0) then
!          do IBFL=IBFS,IBFE
!          ICFL=LBC_SSF(IBFL)
!          ICFP=LCYCSF(IBFL)
!          kdbp(ICFL)=1
!          kdbp(ICFP)=1
!          enddo
!        else
!          do IBFL=IBFS,IBFE
!          ICFL=LBC_SSF(IBFL)
!          ICFP=LCYCSF(IBFL)
!          kdbp(ICFL)=1
!          kdbp(ICFP)=1
!          enddo
!        endif
      endif  
!
      if(kdp.eq.kpdirc) then
!        do IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        ICFP=LCYCSF(IBFL)
!        dum1=0.5d0*(POFLUX(ICFL,IPOT)-POFLUX(ICFP,IPOT))
!        POFLUX(ICFL,IPOT)=0.d0*dum1
!        POFLUX(ICFP,IPOT)=-0.d0*dum1
!        enddo
      elseif(kdp.eq.kpneum) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        POFLUX(ICFL,IPOT)=SFAREA(4,ICFL)*vapotn(IPOT,nb) 
        enddo
        if(IIMAT/=IIMAT2) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          POFLUX(ICFP,IPOT)=SFAREA(4,ICFP)*vapotn(IPOT,nb) 
          enddo
        endif
      elseif(kdp.eq.kpintr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        dum1=0.5d0*(POFLUX(ICFL,IPOT)-POFLUX(ICFP,IPOT))
        POFLUX(ICFL,IPOT)= dum1
        POFLUX(ICFP,IPOT)=-dum1
        enddo
      elseif(kdp.eq.kpindr) then 
        if(INT(vapotn(IPOT,nb))>0) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          POFLUX(ICFL,IPOT)=0.d0
          POFLUX(ICFP,IPOT)=0.d0
          enddo
        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          POFLUX(ICFL,IPOT)=0.d0
          POFLUX(ICFP,IPOT)=0.d0
          enddo

        endif
      endif

 1000 enddo
!
      end subroutine bc_poten_flux
!
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine potential_src(iter,deltt,time,IPOT,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  srcp,diagp
     &   )
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter,IPOT
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3,MXCVFAC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal ( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(inout) :: srcp (MXALLCVP,MXPOTN)
      real*8 ,intent(inout) :: diagp(MXALLCVP,MXPOTN)
!
!
!


!
      return
      end subroutine potential_src
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_setbnd_poten
     &   (iter,time,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &   POTNBC)
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
     &     rotsld,idis,LBC_pair,kdbuff
      use module_boundary,only : kpdirc,kpneum,kpintr,kdpotn,vapotn,
     &                           kpindr
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter
      REAL*8,INTENT(IN)     :: time
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF    (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      REAL*8,INTENT(OUT)    :: POTNBC    (MXSSFBCP,MXPOTN)
!
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IPOT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
!
!
      do IPOT=1,Npotn
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kdp=kdpotn(IPOT,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
!
      if(kdp.eq.kpdirc) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        POTNBC(IBFL,IPOT)=vapotn(IPOT,nb)
        enddo
      elseif(kdp.eq.kpneum) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        POTNBC(IBFL,IPOT)=vapotn(IPOT,nb)
        enddo
      endif
 1000 enddo
      enddo
!
!
      return
      end subroutine bc_setbnd_poten
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_sigma_potn(IPOT,nko,LVEDGE,LBC_SSF,LCYCSF,
     &                     mat_cal,MAT_NO,aa)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_boundary,only : kdbcnd,kdprdc,LBC_INDEX,nbcnd,kdtchi,
     &                           MAT_BCIDX,kdintr,kdsld,kdcvd,idis,
     &                           kdbuff,kdpotn,kpdirc,kpneum,kpintr,
     &                           kpindr
      use module_model,only    : ical_vect
!      use module_metrix,only   : tmpfac=>d2vect
      use module_material,only : nflud,nofld
!
! 1. Set scalar at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko,IPOT
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      integer,intent(in)    :: MAT_NO    (   0:MXMAT)
      real*8 ,intent(inout) :: aa        (  MXALLCV,nko)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMD,IMAT2,IIMAT2,IMAT_U,IMAT_U2,I
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,ICVBO,ICFO
      real*8  :: wi1,wi2
!
! --- 
!
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
! --- Periodic & interface BC:
        if(idis(nb)==0) then
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do 201 k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          aa(IDCP,k)=aa(ICV,k)
  201     continue
 1250     continue
        elseif(idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
!
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          do k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          enddo
          enddo
        endif
      elseif(kd==kdbuff.or.kd==kdintr) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do k=1,nko
          aa(IDC,k)=aa(ICVP,k)   !5555
          aa(IDCP,k)=aa(ICV,k)
          enddo
        enddo
      elseif(kd==kdsld) then
      else
        do 1150 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 101 k=1,nko
          aa(IDC,k)=aa(ICV,k)
  101     continue
 1150     continue
      endif
!
      kdp=kdpotn(IPOT,nb)
      if(kdp.eq.kpdirc.or.kdp.eq.kpneum) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do k=1,nko
          aa(IDC,k)=aa(ICV,k)
          enddo
        enddo
        if(IIMAT/=IIMAT2) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFP)
          IDC=LVEDGE(2,ICFP)
          do k=1,nko
          aa(IDC,k)=aa(ICV,k)
          enddo
          enddo
        endif
      elseif(kdp==kpintr.or.kdp==kpindr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do k=1,nko
        aa(IDC,k)=aa(ICV,k)   !5555
        aa(IDCP,k)=aa(ICVP,k)
        enddo
        enddo
      endif
!
 1000 continue
!
      return
!
      end subroutine dc_sigma_potn
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_out_potn(iko,nko,LVEDGE,LBC_SSF,LCYCSF,POTNBC,
     &                     wiface,mat_cal,MAT_NO,aa)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_boundary,only : kdbcnd,kdprdc,LBC_INDEX,nbcnd,kdtchi,
     &                           MAT_BCIDX,kdintr,kdsld,kdcvd,idis,
     &                           kdbuff,kdpotn,kpdirc,kpneum
      use module_boundary,only : kpdirc,kpneum,kpintr,kdpotn,vapotn,
     &                           kpindr
      use module_material,only : nflud,nofld
!
! 1. Set scalar at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko,iko
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      integer,intent(in)    :: MAT_NO    (   0:MXMAT)
      REAL*8 ,INTENT(IN)    :: WIFACE    (   MXCVFAC)
      real*8 ,intent(inout) :: aa        (   MXALLCV,nko)
      REAL*8,INTENT(IN)  :: POTNBC  (  MXSSFBCP,MXPOTN)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,IMD,IMAT2,IIMAT2,IMAT_U,IMAT_U2,I
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,ICVBO,ICFO
      real*8  :: wi1,wi2,dum1
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IIMAT2=MAT_BCIDX(nb,2)
      IMAT=MAT_NO(IIMAT)
      IMAT2=MAT_NO(IIMAT2)
      if(.not.mat_cal(IIMAT)) cycle 
      kd=kdbcnd(0,nb)
      kdp=kdpotn(iko,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kdp==kpdirc) then
! --- Periodic & interface BC:
        do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          aa(IDC,iko)=POTNBC(IBFL,iko)
 1250   continue
      elseif(kdp==kpintr) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          aa(IDC,iko)=aa(ICVP,iko)
          aa(IDCP,iko)=aa(ICV,iko)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dum1=wi1*aa(IDC,iko)+wi2*aa(IDCP,iko)
          aa(ICVP,iko)=dum1
          aa(ICV,iko)=dum1
        enddo
      elseif(kdp==kpneum) then
      else
        do 1150 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          aa(ICV,iko)=aa(IDC,iko)
 1150     continue
      endif
!
 1000 continue
!
      return
!
      end subroutine dc_out_potn

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_symprv_potn
     &(nko,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & SFAREA,SFCENT,aa,imode)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     imode==0 : velocity
!     imode==1 : grdc
!     imode==2 : rvx
!
! --- [module arguments]
!
      USE module_dimension
      use module_hpcutil,only  : my_rank
      use module_boundary,only : set_rotmtrx,nbcnd,LBC_INDEX,
     &                           kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdintr,kdsld,MAT_BCIDX,rotsld,
     &                           idis,LBC_pair,kdbuff
      use module_material,only : rotati,ishaft,end,begin,nsplit,
     &                           rot_ang
      use module_nowtime ,only : iter,time
      use module_material,only : ical_sld
      use module_model,only    : ical_vect
!      use module_metrix,only   : tmpfac=>d2vect
!
! 1. Set vector at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko,imode
      integer,intent(in)    :: MAT_NO    (  0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(inout) :: aa        (  MXALLCV,3,nko)
!
! --- [local entities]
!
      real*8  :: prd
      real*8  :: rot(3,3)
      real*8  :: shft(3),ufix(3,2)
      real*8  :: unit(3,2),th(2),bb(3,3,2)
      real*8  :: rbb(3,3,2)
      real*8  :: urot(3,2)
      real*8  :: vr(3,2),r(3,2)
      real*8  :: v0(3,2),dr(2),vell(3),wi1,wi2,dum1,dum2,dum3
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2)
      integer :: i,j,ICVBO,ICFO
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k,ISLD,ISLD2
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc)then  !zhang???
         rot=0.d0
         call set_rotmtrx(nbcnd,kd,nb,rot)
!--< 1.1 periodic boundary >--
        do 1150 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 101 k=1,nko
!
        aa(IDC,1,k)=rot(1,1)*aa(ICVP,1,k)
     &             +rot(1,2)*aa(ICVP,2,k)
     &             +rot(1,3)*aa(ICVP,3,k)
        aa(IDC,2,k)=rot(2,1)*aa(ICVP,1,k)
     &             +rot(2,2)*aa(ICVP,2,k)
     &             +rot(2,3)*aa(ICVP,3,k)
        aa(IDC,3,k)=rot(3,1)*aa(ICVP,1,k)
     &             +rot(3,2)*aa(ICVP,2,k)
     &             +rot(3,3)*aa(ICVP,3,k)
!
        if(idis(nb)==0) then   !zhang???
          aa(IDCP,1,k)=rot(1,1)*aa(ICV,1,k)
     &                +rot(2,1)*aa(ICV,2,k)
     &                +rot(3,1)*aa(ICV,3,k)
          aa(IDCP,2,k)=rot(1,2)*aa(ICV,1,k)
     &                +rot(2,2)*aa(ICV,2,k)
     &                +rot(3,2)*aa(ICV,3,k)
          aa(IDCP,3,k)=rot(1,3)*aa(ICV,1,k)
     &                +rot(2,3)*aa(ICV,2,k)
     &                +rot(3,3)*aa(ICV,3,k)
        endif
 101    continue
 1150   continue
      elseif(kd==kdbuff.or.kd==kdintr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do k=1,nko
!
        aa(IDC,1,k)=aa(ICVP,1,k)   !5555
        aa(IDC,2,k)=aa(ICVP,2,k)
        aa(IDC,3,k)=aa(ICVP,3,k)
        aa(IDCP,1,k)=aa(ICV,1,k)
        aa(IDCP,2,k)=aa(ICV,2,k)
        aa(IDCP,3,k)=aa(ICV,3,k)

        enddo
        enddo
      elseif(kd==kdsld.and.idis(nb)==1) then
      elseif(kd==kdsld.and.idis(nb)==2) then
      elseif(kd.eq.kdsymm) then
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 201 k=1,nko
          prd=2.d0*(SFAREA(1,ICFL)*aa(ICV,1,k)
     &             +SFAREA(2,ICFL)*aa(ICV,2,k)
     &             +SFAREA(3,ICFL)*aa(ICV,3,k))
          aa(IDC,1,k)=aa(ICV,1,k)-prd*SFAREA(1,ICFL)
          aa(IDC,2,k)=aa(ICV,2,k)-prd*SFAREA(2,ICFL)
          aa(IDC,3,k)=aa(ICV,3,k)-prd*SFAREA(3,ICFL)
  201     continue
 1250     continue
!      elseif(kd.eq.kdintr) then
!        do 1550 IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        ICV=LVEDGE(1,ICFL)
!        IDC=LVEDGE(2,ICFL)
!        ICFP=LCYCSF(IBFL)
!        ICVP=LVEDGE(1,ICFP)
!        IDCP=LVEDGE(2,ICFP)
!        do 501 k=1,nko
!        aa(IDC,1,k)=aa(ICV,1,k)
!        aa(IDC,2,k)=aa(ICV,2,k)
!        aa(IDC,3,k)=aa(ICV,3,k)
!
!        aa(IDCP,1,k)=aa(ICVP,1,k)
!        aa(IDCP,2,k)=aa(ICVP,2,k)
!        aa(IDCP,3,k)=aa(ICVP,3,k)
! 501    continue
! 1550   continue
      elseif(kd/=kdtchi) then
        do 1350 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 301 k=1,nko
          aa(IDC,1,k)=aa(ICV,1,k)
          aa(IDC,2,k)=aa(ICV,2,k)
          aa(IDC,3,k)=aa(ICV,3,k)
  301     continue
 1350   continue
      endif
 1000 continue
!
      return
      end subroutine dc_symprv_potn
!

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine J_SRC_PEFC_part(iter,IREAC,deltt,time,IPOT,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  yys,tmp,prs,rho,POTNAL,aks,srcp,diagp,J_SRC,FIELD_U)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle
      use module_scalar  ,ONLY : ip_pefc,potname,
     &                           ical_FC,ical_s,iterPOTEN
      use module_chemreac,only : Butler_Volmer,iBV,
     &           BV_ref_I,BV_a,BV_para,BV_ref_molfrc,vreq,
     &                           Faraday,nneq,ig_iter,
     &                           kind_chem,isurface,ireq
      use module_species,only  : r_wm,gascns,I_H2,I_O2,I_H2O
      use module_boundary,only : kdbcnd,kdcvd,nbcnd,MAT_BCIDX,
     &                           IDX_SUFRAC,
     &                           phs_idx,phs_dns,
     &                           phs_nbc,phs_typ,phs_snm,phs_nam,
     &                          gasphase,surphase,blkphase,FC_CELL_VOLT,
     &                           LBC_INDEX,phs_com,num_site,surfreac,
     &                           kdbuff,kdfire,kdcvd,kdintr,mem_cat
      use module_FUEL    ,only : No_Mem,No_AGDL,No_CGDL,OPEN_V,Tau_diff,
     &                           A_cat_th,C_cat_th
      use module_material,only : nflud,nofld
      use module_metrix,only   : lreac
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter,IPOT,IREAC
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3,MXCVFAC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal ( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: yys(MXALLCV,MXcomp)
      real*8 ,intent(in)    :: tmp(MXALLCV)
      real*8 ,intent(in)    :: prs(MXALLCV)
      real*8 ,intent(in)    :: rho(MXALLCV)
      REAL*8 ,INTENT(IN)    :: POTNAL (MXALLCVP,MXPOTN,2)
      real*8 ,intent(in)    :: aks    (MXALLCVR,mxrans)
      real*8 ,intent(inout) :: srcp (MXALLCVP,MXPOTN)
      real*8 ,intent(inout) :: diagp(MXALLCVP,MXPOTN)
      real*8 ,intent(inout) :: J_SRC (MXALLCVP,2)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
!
! --- [local entities]
!
      real*8  :: dum1,dum2,dum3,dum4,dum5,over,TTT,over_CA
      real*8  :: Cref_H2=1.d6,Cref_O2=1.d6
      real*8  :: factor=1.d0,dumm,SBYV,OPEN_V1
      integer :: ICVS,ICVE,ICVL,IDCS,IDCE,IMAT_U,ICOM,ICVL1,IMAT_U2
      integer :: IIMAT,IIMAT2,IMAT,IMAT2,ISEC,nb,kd,IIPOT,ICVLM,ICVLP
      integer :: IBFS,IBFE,IBFL,ICFL,IDC,IRC,ICFP,ICVP,IDCP,ICV,ICVLS
!      logical,allocatable  :: lreac(:)
      integer :: ierr1=0
      logical :: surf_reac=.false.
!
! --- ------------------------------------------------------
!
      if(iBV==0) return
      IIPOT=0
      if(IPOT==ip_pefc(1)) then
        IIPOT=1
        factor=-1.d0
      elseif(IPOT==ip_pefc(2)) then
        IIPOT=2
        factor=1.d0
      endif
!
      J_SRC(:,IPOT)=0.d0
      if(iBV==0.or.iter<=iterPOTEN+1) return
!      
      lreac(:)=.false.
      do IRC=1,nneq
        lreac(IRC)=iter>ig_iter(IRC)
      enddo
!
!-----------------------------------------
! --- dum1 => mol concentration [mol/m^3]
!-----------------------------------------
!
      do 100 nb=1,nbcnd
        kd=kdbcnd(0,nb)
        surf_reac=(kd==kdbuff).or.
     &            (kd==kdfire).or.
     &            (kd==kdcvd).or.
     &            (kd==kdintr)
        if(.not.surf_reac) cycle
        if(surfreac(nb)==2.or.surfreac(nb)==1) then
          do 200 IRC=1,nneq
            if(IDX_SUFRAC(nb,IRC)) then
              if(ireq(IRC)==Butler_Volmer) then
                if(lreac(IRC)) then
                  if(kind_chem(IRC)==isurface) then
        IIMAT=MAT_BCIDX(nb,1)
        IIMAT2=MAT_BCIDX(nb,2)
        IMAT=MAT_NO(IIMAT)
        IMAT2=MAT_NO(IIMAT2)
        IMAT_U=nofld(abs(IMAT))
        IMAT_U2=nofld(abs(IMAT2))
        ISEC=0  !1=>0
        if(IMAT_U==No_Mem) ISEC=1   !0=>1
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(IRC==1) then
!          if(IIPOT==1) factor=-1.d0    !(2)
!          if(IIPOT==2) factor=1.d0     !(4)
          do 300 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          ICVLM=ICV*ISEC+ICVP*(1-ISEC)
          ICVLP=ICV*(1-ISEC)+ICVP*ISEC !IDC*ISEC+IDCP*(1-ISEC)
          ICVL=ICVP*ISEC+ICV*(1-ISEC)
          ICVLS=ICVLM
          if(ip_pefc(1)==IPOT) ICVLS=ICVLP
          TTT=tmp(ICVLP)
          dum1=min(0.d0,POTNAL(ICVLP,ip_pefc(1),2))
! --- 
          dum2=min(POTNAL(ICVLM,ip_pefc(2),2),
     &         POTNAL(ICVLP,ip_pefc(1),2))
          over=min(max(dum1-dum2,0.d0),-dum2)

          dum1=POTNAL(ICVLP,ip_pefc(1),2)
          dum2=POTNAL(ICVLM,ip_pefc(2),2)
          
          if(iter>IREAC) then
            over=max(0.d0,dum1-dum2)
            over=dum1-dum2
          else
            over=-dum1
          endif
!          if(IBFL==IBFS) then
!          endif
! ------------------------------- 
! --- anode 
! --- --------------------------- 
! --- BV Formule 
! --- --------------------------- 
!          dum1=rho(ICVLP)*yys(ICVLP,I_H2)*r_wm(I_H2)
!          dum2=(dum1/BV_ref_molfrc(I_H2,1,IRC))**BV_para(I_H2,1,IRC)
!          dum2=BV_ref_I(IRC)*dum2*(1.d0-aks(ICVLP,ical_s))
!          dum3= exp(BV_a(1,IRC)*Faraday/(gascns*TTT)*over)
!!          dum4=-exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA)
!          dum4=-exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over)
!          dum5=SFAREA(4,ICFL)*C_cat_th/CVVOLM(ICVLS)
!          dumm=dum2*(dum3+dum4)
! ------------------------------
! --- Linearzed BV 
! ------------------------------
          dum1=rho(ICVLP)*yys(ICVLP,I_H2)*r_wm(I_H2)
          dum3=(dum1/BV_ref_molfrc(I_H2,1,IRC))**BV_para(I_H2,1,IRC)
          dum2=BV_ref_I(IRC)*dum3*(1.d0-aks(ICVLP,ical_s)) 
          dum3=BV_a(1,IRC)+BV_a(2,IRC)
          dum5=SFAREA(4,ICFL)*C_cat_th/CVVOLM(ICVLS)
          dumm=dum3*over*dum2*Faraday/(gascns*TTT) 
!   ------------------------------
!          if(IBFL==IBFS) then
!          endif
          J_SRC(ICVL,IIPOT)=dumm                     !(A/m^3) 
          srcp(ICVLS,IPOT)=factor*dumm*dum5
          FIELD_U(ICVLS,2)=over
          FIELD_U(ICVLP,2)=over
 300      enddo
        elseif(IRC==2) then
          do 400 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          ICVLM=ICV*ISEC+ICVP*(1-ISEC)
          ICVLP=ICV*(1-ISEC)+ICVP*ISEC  !IDC*ISEC+IDCP*(1-ISEC)
          ICVL=ICVP*ISEC+ICV*(1-ISEC)
          ICVLS=ICVLM
          if(ip_pefc(1)==IPOT) ICVLS=ICVLP
          TTT=tmp(ICVLP)
!          OPEN_V1=0.0025*TTT+0.2329d0
          OPEN_V1=1.23d0-0.9d-3*(TTT-298.d0)
!------------------------------
! --- cathode 
! --- -------------------------
! --- BV Formule
! --- -------------------------
!          dum1=POTNAL(ICVLP,ip_pefc(1),2)
!          over_CA=dum1-OPEN_V1
!          dum1=rho(ICVLP)*yys(ICVLP,I_O2)*r_wm(I_O2)
!          dum2=(dum1/BV_ref_molfrc(I_O2,1,IRC))**BV_para(I_O2,1,IRC)
!          dum1=rho(ICVLP)*yys(ICVLP,I_H2O)*r_wm(I_H2O)
!          dum2=dum2
!     &        *(dum1/BV_ref_molfrc(I_H2O,2,IRC))**BV_para(I_H2O,2,IRC)
!          dum2=BV_ref_I(IRC)*dum2*(1.d0-aks(ICVLP,ical_s))
!          dum3= exp(BV_a(1,IRC)*Faraday/(gascns*TTT)*over_CA)
!          dum4=-exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA)
!          SBYV=SFAREA(4,ICFL)*C_cat_th/CVVOLM(ICVLS) 
!          dumm=dum2*(dum3+dum4)
!-------------------------------
! --- Linearzed BV -1
!-------------------------------
!          over_CA=-OPEN_V1!*0.5
!          dum1=rho(ICVLP)*yys(ICVLP,I_O2)*r_wm(I_O2)
!          dum2=(dum1/BV_ref_molfrc(I_O2,1,IRC))
!          dum2=-BV_ref_I(IRC)*dum2*(1.d0-aks(ICVLP,ical_s))
!          dum3=!exp(-16456.d0*(1.d0/TTT+1.d0/353.15))*
!     &        exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA)
!          SBYV=SFAREA(4,ICFL)*C_cat_th/CVVOLM(ICVLS) 
!          dumm=dum2*dum3
!-------------- ----------------
! --- Linearzed BV -2 
!-------------- ----------------
          dum1=min(OPEN_V1,max(FC_CELL_VOLT,POTNAL(ICVLP,ip_pefc(1),2)))
          dum3=dum1-OPEN_V1
          dum2=max(min(0.d0,POTNAL(ICVLM,ip_pefc(2),2)),dum3)
          over_CA=max(min(0.d0,dum1-dum2-OPEN_V1),-OPEN_V1)
!
          dum1=POTNAL(ICVLP,ip_pefc(1),2)
          dum2=POTNAL(ICVLM,ip_pefc(2),2)
!
          if(iter>IREAC) then
            over_CA=dum1-dum2-OPEN_V1
          else
            over_CA=FC_CELL_VOLT-OPEN_V1
          endif
!
          dum1=rho(ICVLP)*yys(ICVLP,I_O2)*r_wm(I_O2) 
          dum2=(dum1/BV_ref_molfrc(I_O2,1,IRC))**BV_para(I_O2,1,IRC) 
          dum2=BV_ref_I(IRC)*dum2*(1.d0-aks(ICVLP,ical_s)) 
          dum4=exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA) 
          SBYV=SFAREA(4,ICFL)*C_cat_th/CVVOLM(ICVLS) 
          dumm=-dum2*dum4 
!          if(IBFL==IBFS) then
!            print*,'Cathode overpotential',IPOT,over_CA,factor*dumm
!          endif
          J_SRC(ICVL,IIPOT)=dumm
          srcp(ICVLS,IPOT)=factor*dumm*SBYV
          FIELD_U(ICVLS,2)=over_CA
          FIELD_U(ICVLP,2)=over_CA
!--------------------------------
 400      enddo
        endif
                  endif
                endif
              endif
            endif
 200      enddo
        endif
 100  enddo
!
      return
      end subroutine J_SRC_PEFC_part


!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine PEFC_sigma(IPOT,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,aks,prs,
     &  yys,tmp,sigma,FIELD_U
     &   )
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle
      use module_initial,only  : elecCondFC,ionCondFC,potn_c
      use module_initial ,only : Ke_FC,Kc_FC,Pv_FC,
     &                           sigma_FC,contang_FC
      use module_scalar  ,ONLY : ip_pefc,potname,
     &                           ical_FC,ical_s
      use module_FUEL  ,only   : vap_no,No_Mem,Tau_diff,OPEN_V
      use module_metrix,only   : rcomp,eva_comp
      use module_species,only  : wm,r_wm,gascns
      use module_model,   only : PEFC,PAFC,MCFC,SOFC,AFC,MDFC
      use module_material,only : nflud,nofld,nosld,porosty
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IPOT
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3,MXCVFAC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal ( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: aks(MXALLCVR,mxrans)
      real*8 ,intent(in)    :: prs(MXALLCV)
      real*8 ,intent(in)    :: yys(MXALLCV,MXcomp)
      real*8 ,intent(in)    :: tmp(MXALLCV)
      real*8 ,intent(inout) :: sigma (MXALLCV)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
!
!
! --- [local entities]
!
      real*8  :: dum1,dum2,dum3,dum4,dum5,errab_POTN(MXPOTN),ppp,TTT
      integer :: ICVA,ICVB,IIMAT,IMAT,ICFL,ICFS,ICFE,IMAT_U,ICOM
      integer :: k,ICVS,ICVE,IDCS,IDCE,ICVL
!
!
!
      sigma(:)=0.d0
      if(IPOT==ip_pefc(1).and.ical_FC==PEFC) then      !'FAI_S'
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        IMAT_U=nofld(abs(IMAT))
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        sigma(ICVS:ICVE)=MAX(elecCondFC(IIMAT),SML)
        FIELD_U(ICVS:ICVE,1)=sigma(ICVS:ICVE)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        sigma(IDCS:IDCE)=MAX(elecCondFC(IIMAT),SML)
        if(elecCondFC(IIMAT)<0.d0.and.IMAT_U==No_Mem) then
          write(ifle,'(1X,a)') 
     &    'ERR: Negative electron conductivity in fflow.ctl'
          write(ifle,'(1X,a,I4)') 'MSG: in IMAT_U = ',IMAT_U
          call FFRABORT
     &  (1,'ERR: Electron conductivity negative value')
        endif
        enddo
      elseif(IPOT==ip_pefc(2).and.ical_FC==PEFC) then  !'FAI_M'
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          IMAT_U=nofld(IMAT)
        elseif(IMAT.lt.0) then
          IMAT_U=nosld(-IMAT)
        endif
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(ionCondFC(IIMAT)<0.d0.and.IMAT_U==No_Mem) then
          do ICVL=ICVS,ICVE
          ppp=prs(ICVL)
          TTT=tmp(ICVL)
          rcomp(:)=yys(ICVL,:)
          dum1=aks(ICVL,ical_s)
          dum2=Pv_FC(IIMAT)
          dum3=aks(ICVL,ical_s)
!
          call ramd(rcomp,ppp,TTT,dum2,dum3,
     &              dum4,r_wm,ncomp,vap_no,dum5)
          if(dum5<0.326d0/0.514d0) dum5=1.d0
! --- 
          dum2=(0.514d0*dum5-0.326d0)
     &        *exp(1268.d0*(1.d0/303.d0-1.d0/TTT))
!          dum2=porosty(IMAT)**Tau_diff*dum2
          sigma(ICVL)=max(-ionCondFC(IIMAT),dum2)
!          sigma(ICVL)=dum2
          FIELD_U(ICVL,1)=sigma(ICVL)
          enddo
        elseif(IMAT_U==No_Mem) then
          sigma(ICVS:ICVE)=max(ionCondFC(IIMAT),SML)
          sigma(IDCS:IDCE)=MAX(ionCondFC(IIMAT),SML)
          FIELD_U(ICVS:ICVE,1)=sigma(ICVS:ICVE)

        else
          sigma(ICVS:ICVE)=max(ionCondFC(IIMAT),SML)!SML!GREAT
          sigma(IDCS:IDCE)=max(ionCondFC(IIMAT),SML)!SML!GREAT
        endif
        enddo
      else
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        sigma(ICVS:ICVE)=potn_c(IIMAT,IPOT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        sigma(IDCS:IDCE)=potn_c(IIMAT,IPOT)
        enddo
      endif
!
      end subroutine PEFC_sigma
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine ramd(rcomp,ppp,TTT,Pv,sss,aaa,r_wm,ncomp,vap_no,ramd1) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      integer,intent(in) :: ncomp,vap_no
      real*8,intent(in)  :: rcomp(ncomp),ppp,TTT,Pv,r_wm(ncomp),sss
      real*8,intent(out) :: aaa,ramd1
!
      real*8 :: dum1,dum2,dum3,dum4,dum5,dum6,Xw
      integer:: ICOM,IIMAT
!
      Xw=0.d0
      do ICOM=1,ncomp
        Xw=Xw+rcomp(ICOM)*r_wm(ICOM) 
      enddo
      dum2=Xw/r_wm(vap_no)
      Xw=rcomp(vap_no)/dum2
      dum1=ppp*Xw
      if(Pv<0.d0) then
        dum5=-2.1794d0
     &        +0.02953d0*(TTT-273.15d0)
     &        -9.1837D-5*(TTT-273.15d0)**2
     &        +1.4454D-7*(TTT-273.15d0)**3
        dum5=1.d5*10.d0**dum5
      else
        dum5=Pv
      endif
      aaa=dum1/dum5+2.d0*sss
!
      if(aaa<1.d0) then
        ramd1=0.043d0+17.81d0*aaa-39.85d0*aaa**2+36.d0*aaa**3 
      elseif(aaa>=1.d0.and.aaa<3.d0) then
        ramd1=14.d0+1.4d0*(aaa-1.d0) 
      else
        ramd1=16.8d0
      endif
!
      return
      end subroutine ramd
      
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine conv_term_FC
     & (icnv,lmtr,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  qfv,vel,dggg,POTNAL,
     &  grdc,rva,q,dqdt,rho,rho2,
     &  tmp,prs,yys,aks,FIELD_U,
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
      use module_metrix,only  : sigma =>W1K8
      use module_metrix,only  : cnvv=>W1K9
      use module_metrix,only  : rcomp,cnvq=>W1K7 !cnvq=>W1K10,rcomp
      USE module_dimension
      use module_material,only : nflud,venkat,rnck,slope,
     &                           up1st,up2nd,up3rd,cnt2nd,bldf,
     &                           usi2nd,sthrmu2,sthrmu,porosty,
     &                           nofld
      use module_Euler2ph,only : ieul2ph
      use module_initial ,only : rho0,rho02,Ke_FC,Kc_FC,Pv_FC,
     &                           sigma_FC,contang_FC,Kap_FC,
     &                           surf_tenFC,ionCondFC
      use module_gravity ,only : ggg
      use module_species,only  : wm,r_wm
      use module_FUEL    ,only : No_Mem,No_AGDL,No_CGDL,osmoticFC,
     &                           vap_no,osmo_Nd,Tau_diff
      use module_chemreac,only : Faraday
      use module_scalar  ,only : ip_pefc,ical_s
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
      real*8 ,intent(inout) :: dggg  (  MXALLCV,3)
      REAL*8 ,INTENT(IN)    :: POTNAL(MXALLCVP,MXPOTN)
      REAL*8 ,INTENT(IN)    :: tmp ( MXALLCV)
      REAL*8 ,INTENT(IN)    :: YYS ( MXALLCV,MXCOMP)
      REAL*8 ,INTENT(IN)    :: prs ( MXALLCV)
      REAL*8 ,INTENT(IN)    :: aks ( MXALLCVR,MXRANS)
      REAL*8,INTENT(INOUT)  :: RHO2( MXALLCVC)
      REAL*8,INTENT(INOUT)  :: RHO ( MXALLCV)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
!
! --- [local entities]
!
      real*8  :: gf1,gf2,dlvect,grx,gry,grz
      real*8 ,parameter :: ak=1.d0/0.3333d0,drop_d=0.001
      real*8  :: dqp,dqm,dqpp,dqmm,dqpm,wi1,wi2
      real*8  :: dx,dy,dz,dab,da,db
      real*8  :: epsw,sgn,dum1,qqq,alhpa,dum2,dum3,dum4
      integer :: i,j,k,l,m,n
      integer :: nb,kd,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,IMAT_U
      integer :: IBFS,IBFE,IBFL
      integer :: ierr=0,IPOT
      logical :: icnTVD
      real*8 ,save :: dxx,dyx,dzx,dl,TTT,PPP,dumA,dumB
      real*8  :: Ramda_l,Ramda_g,ru,rv,rw
!
! --- START
!
!
!-< 1. 1st order upwind only for PEFC>- 
!
      
      do 100 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 100
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) goto 100
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        qfv(ICFL,1)=q(ICVLA)
        qfv(ICFL,2)=q(ICVLB)
      enddo
!
  100 enddo
!-----------------------------------
!-< 5. Calculate convection term >-
!-----------------------------------
      cnvq=0.d0
      cnvv=0.d0
      do 310 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle

      IMAT_U=nofld(IMAT)
!      if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL)) cycle
      if(IMAT_U==No_Mem) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVLA=LVEDGE(1,ICFL)
      ICVLB=LVEDGE(2,ICFL)
      dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)
     &    +min(0.d0,rva(ICFL))*qfv(ICFL,2)
      cnvq(ICVLA)=cnvq(ICVLA)+dum1
      cnvq(ICVLB)=cnvq(ICVLB)-dum1
      cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)
      cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)
      enddo
  310 enddo
!
      do 320 IIMAT=1,NMAT    !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IMAT_U=nofld(IMAT)
!      if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL)) cycle
      if(IMAT_U==No_Mem) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do 330 ICVL=ICVS,ICVE
      dqdt(ICVL)=dqdt(ICVL)+2.d0*(q(ICVL)*cnvv(ICVL)-cnvq(ICVL))
 330  enddo
 320  enddo
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,dqdt)
      ENDIF
!-------------------------------
! --- 
!-------------------------------
      cnvq=0.d0
      cnvv=0.d0
      do 610 IIMAT=1,NMAT    !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IMAT_U=nofld(IMAT)      
!      if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL)) cycle
      if(IMAT_U==No_Mem) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVLA=LVEDGE(1,ICFL)
      ICVLB=LVEDGE(2,ICFL)
      dum1=max(0.d0,rva(ICFL))*qfv(ICFL,1)**2
     &      +min(0.d0,rva(ICFL))*qfv(ICFL,2)**2
      cnvq(ICVLA)=cnvq(ICVLA)+dum1
      cnvq(ICVLB)=cnvq(ICVLB)-dum1
      cnvv(ICVLA)=cnvv(ICVLA)+rva(ICFL)
      cnvv(ICVLB)=cnvv(ICVLB)-rva(ICFL)
      enddo
  610 enddo
!
      do 620 IIMAT=1,NMAT    !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IMAT_U=nofld(IMAT)
!      if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL)) cycle
      if(IMAT_U==No_Mem) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do 630 ICVL=ICVS,ICVE
      dqdt(ICVL)=dqdt(ICVL)-(q(ICVL)**2*cnvv(ICVL)-cnvq(ICVL))
 630  enddo
 620  enddo
!-------------------------------
! --- porosity 
!-------------------------------
      do IIMAT=1,NMAT    !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT>0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          dum1=porosty(IMAT)
          do ICVL=ICVS,ICVE
          dqdt(ICVL)=dqdt(ICVL)*dum1
          enddo
        else
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          dqdt(ICVS:ICVE)=0.d0
        endif
      enddo
!-------------------------------
! ---- g term
!-------------------------------
      dggg(:,:)=0.d0
      do IIMAT=1,NMAT    !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IMAT_U=nofld(IMAT)      
      if(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem) then
!          if(IMAT_U==No_Mem) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        dxx=rho02(IIMAT)
        dx=sthrmu2(IMAT)
        do ICVL=ICVS,ICVE
        dyx=rho2(ICVL)   !rho0(IIMAT)
        dy=sthrmu(IMAT)
        dum2=(rho02(IIMAT)-rho2(ICVL))*Kap_FC(IIMAT)
        dum1=q(ICVL) 
        Ramda_l=dum1*(2.d0-dum1) 
        Ramda_g=1.d0-Ramda_l
        dum3=dxx*dum1/dx+dyx*(1.d0-dum1)/dy
        dggg(ICVL,:)=Ramda_l*Ramda_g*dum2*ggg(:)*dum3
!        dggg(ICVL,:)=dum2*ggg(:)*dum3
        enddo
      else
!        cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dum1=q(ICVL)
        dum2=rho(ICVL)
        dum3=dsqrt((rho02(IIMAT)-dum2)/dum2*ak*drop_d)
        dggg(ICVL,:)=dum3*ggg(:)*dum1*dum2
        enddo
      endif
      enddo
!
      call dc_symprv_potn
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    SFAREA,SFCENT,dggg,0)
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IMAT_U=nofld(IMAT)
!      if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem)) 
!     &  cycle
!      if(IMAT_U==No_Mem) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
      wi1=wiface(ICFL)
      wi2=1.d0-wiface(ICFL)
      ru=(wi1*dggg(ICVA,1)+wi2*dggg(ICVB,1))
      rv=(wi1*dggg(ICVA,2)+wi2*dggg(ICVB,2))
      rw=(wi1*dggg(ICVA,3)+wi2*dggg(ICVB,3))
      dum1=SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw
      qfv(ICFL,1)=SFAREA(4,ICFL)*dum1
      enddo
      enddo
!
      cnvv=0.d0
      do IIMAT=1,NMAT    !ICF=1,NCVFAC
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      if(IMAT.lt.0) cycle
      IMAT_U=nofld(IMAT)
!      if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem)) 
!     &  cycle
!      if(IMAT_U==No_Mem) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVLA=LVEDGE(1,ICFL)
      ICVLB=LVEDGE(2,ICFL)
      cnvv(ICVLA)=cnvv(ICVLA)+qfv(ICFL,1)
      cnvv(ICVLB)=cnvv(ICVLB)-qfv(ICFL,1)
      enddo
      enddo
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      IMAT_U=nofld(IMAT)
!      if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem)) 
!     &  cycle
!      if(IMAT_U==No_Mem) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      dqdt(ICVL)=dqdt(ICVL)-cnvv(ICVL)!*Ramda_l*Ramda_g
      enddo
      enddo
!-----------------------------------------------
! ---- i term (Electro-osmotic drag)   :No_Mem 
!-----------------------------------------------
      if(osmoticFC==1) then
        cnvv(:)=POTNAL(:,ip_pefc(2))   ! FAI_M
        call grad_cell(1,6,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cnvv,grdc)
        call dc_symprv_potn
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    SFAREA,SFCENT,grdc,0)
!
        IPOT=ip_pefc(2)
        call PEFC_sigma(IPOT,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,aks,prs,
     &  yys,tmp,sigma,FIELD_U
     &   )
!
        call  dc_sigma_potn(IPOT,1,LVEDGE,LBC_SSF,LCYCSF,
     &                 mat_cal,MAT_NO,sigma)
!
        do IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        FIELD_U(ICVL,3:5)=-grdc(ICVL,1:3)*sigma(ICVL)
        enddo
        enddo
!
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        IMAT_U=nofld(IMAT)
        if(IMAT_U==No_Mem) then
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dum1=osmo_Nd*wm(vap_no)/Faraday
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dumA=sigma(ICVA)
          dumB=sigma(ICVB)
          dum2=dumA*dumB/(wi1*dumB+wi2*dumA+SML)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
          dx=dx/dl
          dy=dy/dl
          dz=dz/dl
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          grx=wi1*grdc(ICVA,1)+wi2*grdc(ICVB,1)
          gry=wi1*grdc(ICVA,2)+wi2*grdc(ICVB,2)
          grz=wi1*grdc(ICVA,3)+wi2*grdc(ICVB,3)
          gf2=grx*(SFAREA(1,ICFL)-dx)
     &       +gry*(SFAREA(2,ICFL)-dy)
     &       +grz*(SFAREA(3,ICFL)-dz)
          gf1=(POTNAL(ICVB,ip_pefc(2))-POTNAL(ICVA,ip_pefc(2)))/dl
          dum1=-dum1*(gf1+gf2)*dum2*SFAREA(4,ICFL)
          dqdt(ICVA)=dqdt(ICVA)-dum1*q(ICVA)
          dqdt(ICVB)=dqdt(ICVB)+dum1*q(ICVB)
          enddo
        endif
        enddo
      endif
!-------------------------------
      return
      end subroutine conv_term_FC

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine FC_src(imode,ns,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,rho2,tmp,aks,src,diagk,prs,yys,mDOT,
     &  LCYCOLD,wifsld)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_dimension
      use module_constant
      use module_initial ,only : rho0,rho02,Ke_FC,Kc_FC,Pv_FC,
     &                           sigma_FC,contang_FC
      use module_initial ,only : Kap_FC
      use module_scalar,  only : ical_FC,ical_s
      use module_species,only  : wm,r_wm,gascns
      use module_metrix,only   : rcomp,eva_comp
      use module_material,only : porosty,nofld
      use module_FUEL    ,only : No_Mem,No_AGDL,No_CGDL,vap_no
      
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LVEDGE (2, MXCVFAC),imode,ns
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
      real*8 ,intent(in)    :: rho    (   MXALLCV)
      real*8 ,intent(in)    :: rho2   (   MXALLCVC)
      real*8 ,intent(in)    :: tmp    (   MXALLCV)
      real*8 ,intent(inout) :: aks    (   MXALLCVR,MXrans)
      real*8 ,intent(inout) :: diagk  (   MXALLCVR,MXrans)

      real*8 ,intent(inout) :: src    (   MXALLCVR,ns)
      integer,intent(in)    :: LCYCOLD(   MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (   MXSSFBC_SLD)
      real*8 ,intent(in)    :: prs(MXALLCV)
      real*8 ,intent(in)    :: yys    (   MXALLCV,MXcomp)
      real*8 ,intent(inout) :: mDOT   (   MXALLCV,3)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,ICOM,IMAT_U
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      real*8  :: dum1,dum2,dum3,dum4,rhol_rhov,Pv,ppp,TTT,Xw,dum5,dum6

!
      if(imode==0) then  !rans_admin
!--------------
! --- mdot
!--------------
        do 600 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle
        IMAT_U=nofld(IMAT)
!        if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem))
!     &  cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO 630 ICVL=ICVS,ICVE
        ppp=prs(ICVL)
        TTT=tmp(ICVL)
        Xw=0.d0
        rcomp(:)=yys(ICVL,:)
        do ICOM=1,ncomp
          Xw=Xw+rcomp(ICOM)*r_wm(ICOM)
        enddo
        dum2=wm(vap_no)*Xw
        Xw=rcomp(vap_no)/dum2
        dum1=ppp*Xw
        dum2=(1.d0-aks(ICVL,ical_s))*porosty(IMAT)
        if(Pv_FC(IIMAT)<0.d0) then
          dum5=
     &        -2.1794d0
     &        +0.02953d0*(TTT-273.15d0)
     &        -9.1837D-5*(TTT-273.15d0)**2
     &        +1.4454D-7*(TTT-273.15d0)**3
          dum5=1.d5*10.d0**dum5
        else
          dum5=Pv_FC(IIMAT)
        endif
        if(dum1>dum5) then  ! conde
          dum3=wm(vap_no)*Kc_FC(IIMAT)*porosty(IMAT)*Xw*
     &        (dum1-dum5)/gascns/TTT
          src(ICVL,ical_s)=src(ICVL,ical_s)+dum3   ![kg/m^3/s]
        else                ! evapor
          dum2=aks(ICVL,ical_s)
          if(dum2<SML15) dum2=0.d0 
          dum3=Ke_FC(IIMAT)*porosty(IMAT)*dum2*rho02(IIMAT)*
     &        (dum1-dum5)
          dum4=Ke_FC(IIMAT)*porosty(IMAT)*rho02(IIMAT)*
     &        (dum1-dum5)/rho02(IIMAT)
          src(ICVL,ical_s)=src(ICVL,ical_s)+dum3   ![kg/m^3/s]
          diagk(ICVL,ical_s)=diagk(ICVL,ical_s)+dum4
        endif   !dum3=mdot=> dum3>0:conde  ; dum3<0:evapor 
        mDOT(ICVL,1)=-dum3
 630    enddo
 600    enddo
      endif
!
!      do 700 IIMAT=1,NMAT
!        IMAT=MAT_NO(IIMAT)
!        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle
!        IMAT_U=IMAT  !nofld(IMAT)
!        if(.NOT.(IMAT_U==No_Mem)) cycle
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_CVEXT(IIMAT)
!        DO 730 ICVL=ICVS,ICVE
!        dum3=
!     & yys(ICVL,vap_no)*rho2(ICVL)*(1.d0-aks(ICVL,ical_s))*porosty(IMAT)
!        src(ICVL,ical_s)=src(ICVL,ical_s)+dum3
!        mDOT(ICVL,1)=-dum3
! 730    enddo
! 700  enddo
!---------------------
! --- 
!---------------------
      end subroutine FC_src
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine J_SRC_PEFC_part1(iter,deltt,time,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM, 
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  yys,tmp,prs,rho,POTNAL,aks,srcp,diagp,J_SRC,FIELD_U)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle
      use module_scalar  ,ONLY : ip_pefc,potname,
     &                           ical_FC,ical_s,iterPOTEN
      use module_chemreac,only : Butler_Volmer,iBV,
     &           BV_ref_I,BV_a,BV_para,BV_ref_molfrc,vreq,
     &                           Faraday,nneq,ig_iter,
     &                           kind_chem,isurface,ireq
      use module_species,only  : r_wm,gascns,I_H2,I_O2,I_H2O
      use module_boundary,only : kdbcnd,kdcvd,nbcnd,MAT_BCIDX,
     &                           IDX_SUFRAC,
     &                           phs_idx,phs_dns,
     &                           phs_nbc,phs_typ,phs_snm,phs_nam,
     &                          gasphase,surphase,blkphase,FC_CELL_VOLT,
     &                           LBC_INDEX,phs_com,num_site,surfreac,
     &                           kdbuff,kdfire,kdcvd,kdintr,mem_cat
      use module_FUEL    ,only : No_Mem,No_AGDL,No_CGDL,OPEN_V,Tau_diff,
     &                           A_cat_th,C_cat_th
      use module_material,only : nflud,nofld
      use module_metrix,only   : lreac
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3,MXCVFAC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal ( 0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: yys(MXALLCV,MXcomp)
      real*8 ,intent(in)    :: tmp(MXALLCV)
      real*8 ,intent(in)    :: prs(MXALLCV)
      real*8 ,intent(in)    :: rho(MXALLCV)
      REAL*8 ,INTENT(IN)    :: POTNAL (MXALLCVP,MXPOTN,2)
      real*8 ,intent(in)    :: aks    (MXALLCVR,mxrans)
      real*8 ,intent(inout) :: srcp (MXALLCVP,MXPOTN)
      real*8 ,intent(inout) :: diagp(MXALLCVP,MXPOTN)
      real*8 ,intent(inout) :: J_SRC (MXALLCVP,2)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
!
! --- [local entities]
!
      real*8  :: dum1,dum2,dum3,dum4,dum5,over,TTT,over_CA
      real*8  :: Cref_H2=1.d6,Cref_O2=1.d6
      real*8  :: factor=1.d0,dumm,SBYV,OPEN_V1
      integer :: ICVS,ICVE,ICVL,IDCS,IDCE,IMAT_U,ICOM,ICVL1,IMAT_U2
      integer :: IIMAT,IIMAT2,IMAT,IMAT2,ISEC,nb,kd,IIPOT,ICVLM,ICVLP
      integer :: IBFS,IBFE,IBFL,ICFL,IDC,IRC,ICFP,ICVP,IDCP,ICV,ICVLS
!      logical,allocatable  :: lreac(:)
      integer :: ierr1=0
      logical :: surf_reac=.false.
!
! --- ------------------------------------------------------
      return

!      if(IPOT==2) then
!        print*,'Ion Potential' 
!      else
!        print*,'Electron Potential'
!      endif
!
      if(iBV==0) return
!      IIPOT=0
!      if(IPOT==ip_pefc(1)) then
!        IIPOT=1
!        factor=-1.d0
!      elseif(IPOT==ip_pefc(2)) then
!        IIPOT=2
!        factor=1.d0
!      endif
!
      J_SRC(:,1:2)=0.d0
      if(iBV==0.or.iter<=iterPOTEN+1) return
      
      lreac(:)=.false.
      do IRC=1,nneq
        lreac(IRC)=iter>ig_iter(IRC)
      enddo
!
!-----------------------------------------
! --- dum1 => mol concentration [mol/m^3] 
!-----------------------------------------
!
      do 100 nb=1,nbcnd 
        kd=kdbcnd(0,nb)
        surf_reac=(kd==kdbuff).or.
     &            (kd==kdfire).or.
     &            (kd==kdcvd).or.
     &            (kd==kdintr)
        if(.not.surf_reac) cycle
        if(surfreac(nb)==2.or.surfreac(nb)==1) then
          do 200 IRC=1,nneq
            if(IDX_SUFRAC(nb,IRC)) then
              if(ireq(IRC)==Butler_Volmer) then
                if(lreac(IRC)) then
                  if(kind_chem(IRC)==isurface) then
        IIMAT=MAT_BCIDX(nb,1)
        IIMAT2=MAT_BCIDX(nb,2)
        IMAT=MAT_NO(IIMAT)
        IMAT2=MAT_NO(IIMAT2)
        IMAT_U=nofld(abs(IMAT))
        IMAT_U2=IMAT2  !nofld(abs(IMAT2))
        ISEC=0  !1=>0 
        if(IMAT_U==No_Mem) ISEC=1   !0=>1
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(IRC==1) then
          do 300 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)
!
          TTT=tmp(ICVLP)
          OPEN_V1=0.0025*TTT+0.2329d0
          dum1=POTNAL(ICVLP,ip_pefc(1),2)
          dum2=POTNAL(ICVLM,ip_pefc(2),2)
          over=dum1-dum2
!          over=0.2d0  !0.05
! ------------------------------- 
! --- anode 
! --- --------------------------- 
! --- BV Formule 
! --- --------------------------- 
!          dum1=rho(ICVL)*yys(ICVL,I_H2)*r_wm(I_H2)
!          dum2=(dum1/BV_ref_molfrc(I_H2,1,IRC))**BV_para(I_H2,1,IRC
!          dum2=BV_ref_I(IRC)*dum2*(1.d0-aks(ICVL,ical_s))
!          dum3= exp(BV_a(1,IRC)*Faraday/(gascns*TTT)*over)
!          dum4=-exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA)
!          dum4=0.d0
!          dumm=dum2*(dum3+dum4)
!          J_SRC(ICVL,IIPOT)=dumm
!          srcp(ICVL,IPOT)=-factor*dumm
! ------------------------------
! --- Linearzed BV 
! ------------------------------
          dum1=rho(ICVLP)*yys(ICVLP,I_H2)*r_wm(I_H2)
          dum3=(dum1/BV_ref_molfrc(I_H2,1,IRC))**BV_para(I_H2,1,IRC)
          dum2=BV_ref_I(IRC)*dum3*(1.d0-aks(ICVLP,ical_s)) 
          dum3=BV_a(1,IRC)+BV_a(2,IRC)
          dum5=SFAREA(4,ICFL)*A_cat_th/CVVOLM(ICVLS)
          dumm=dum3*over*dum2*Faraday/(gascns*TTT)   !*SBYV 

!          dumm=1.d8
!          J_SRC(ICVL,IIPOT)=dumm             !(A/m^3)
!          srcp(ICVLS,IPOT)=factor*dumm*dum5
          if(IBFL==IBFS) then
            print*,'============================================'
            print*,'Anode side, BC no=',nb
            print*,'Reaction no=',IRC
          endif
 300      enddo
        elseif(IRC==2) then
          do 400 IBFL=IBFS,IBFE 
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)

          TTT=tmp(ICV)
          OPEN_V1=0.0025*TTT+0.2329d0
!          OPEN_V1=1.23d0-0.9d-3*(TTT-298.d0)   !OPEN_V1=0.0025*TTT+0.2329d0
!------------------------------
! --- cathode 
! --- -------------------------
! --- BV Formule
! --- -------------------------
!          dum1=rho(ICVL)*yys(ICVL,I_O2)*r_wm(I_O2)
!          dum2=(dum1/BV_ref_molfrc(I_O2,1,IRC))**BV_para(I_O2,1,IRC)
!          dum1=rho(ICVL)*yys(ICVL,I_H2O)*r_wm(I_H2O)
!          dum2=dum2  !??? H2O-G 
!     &        *(dum1/BV_ref_molfrc(I_H2O,2,IRC))**BV_para(I_H2O,2,IRC)
!          dum2=-BV_ref_I(IRC)*dum2*(1.d0-aks(ICVL,ical_s))
!          dum3= exp(BV_a(1,IRC)*Faraday/(gascns*TTT)*over)
!          dum4=-exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA)
!          dum3=0.d0
!          dumm=dum2*(dum3+dum4)
!          J_SRC(ICVL,IIPOT)=dumm
!          srcp(ICVL,IPOT)=-factor*dumm
!------------------------------
! --- Linearzed BV -1
!------------------------------
!          dum1=rho(ICVL)*yys(ICVL,I_O2)*r_wm(I_O2)
!          dum2=(dum1/BV_ref_molfrc(I_O2,1,IRC))
!          dum2=-BV_ref_I(IRC)*dum2*(1.d0-aks(ICVL,ical_s))
!          dum3=exp(-16456.d0*(1.d0/TTT+1.d0/353.15))
!     &        *exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA)
!          dumm=dum2*dum3!*SBYV
!          J_SRC(ICVL,IIPOT)=dumm
!          srcp(ICVL,IPOT)=factor*dumm
!-------------- ----------------
! --- Linearzed BV -1 
!-------------- ----------------
          dum1=POTNAL(ICVLP,ip_pefc(1),2)
          dum2=POTNAL(ICVLM,ip_pefc(2),2)
          over_CA=dum1-dum2-OPEN_V1
          if(IBFL==IBFS) then
            print*,'============================================'
            print*,'Cathode side, BC no=',nb
            print*,'Reaction no=',IRC
          endif
!          over_CA=-0.2d0  !-0.35
          dum1=rho(ICVLP)*yys(ICVLP,I_O2)*r_wm(I_O2)
          dum2=(dum1/BV_ref_molfrc(I_O2,1,IRC))**BV_para(I_O2,1,IRC)
          dum2=BV_ref_I(IRC)*dum2*(1.d0-aks(ICVLP,ical_s))
          dum4=exp(-BV_a(2,IRC)*Faraday/(gascns*TTT)*over_CA)
          dum5=SFAREA(4,ICFL)*C_cat_th/CVVOLM(ICVLS)
          dumm=dum2*dum4

!          dumm=-1.d8
!          J_SRC(ICVL,IIPOT)=dumm
!          srcp(ICVLS,IPOT)=factor*dumm*dum5
!------------------------------

 400      enddo

        endif
                  endif
                endif
              endif
            endif
 200      enddo
        endif
 100  enddo
!
      return
      end subroutine J_SRC_PEFC_part1

