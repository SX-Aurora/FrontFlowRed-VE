!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     subroutine solve_cnvdif
!     subroutine solve_poisson
!     subroutine solve_poisson_unsymm
!     subroutine solve_poisson_unsymm_vect
!     subroutine solve_poisson_unsymm_e2p
!     subroutine solve_poisson_e2p
!     subroutine solve_cnvdif_e2p1
!     subroutine solve_cnvdif_vect
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_cnvdif
     & (clrsld,NSSF,ndf,nq,cn,time,relax,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,SFCENT,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdbf,rva,rvd,dif,diag,
     &  LCYCOLD,wifsld,OPPANG,CVVOLM,
     &  cofd,dscl,htcbnd,rhs,deltt,
     &  iterph,reps,aeps,iter,ierror,IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine is for calculating NS equations
      use module_dimension
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IQ =>IW2K1 
      use module_metrix,only   : bb
      use module_metrix,only   : aalw
      use module_metrix,only   : dmax=>W1K7
      use module_metrix,only   : dmin=>W1K8
      use module_metrix,only   : alx =>W1K9
      use module_material,only : ical_sld
      use module_model,only    : ical_vect,nthrds,ical_MHD
      use module_material,only : relax_MHD
      use module_time,    only : i_steady
!
! --- [module arguments]
!
      use module_hpcutil
      use module_io,only       : ifll,ifle,cntlnam
      use module_cgsolver,only : aepsbcg,repsbcg,iterbcg,
     &                           aepsbcg_MHD,repsbcg_MHD,
     &                           iterbcg_MHD
      use module_flags,only    : 
     &                           Adams_Moulton,intgvv,Crank_Nicolson
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,ical_buff,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang
      use module_species,only  : spcnam
      use module_scalar ,only  : sclname
      use module_metrix,only   : aax,ipx
      use module_vector,ONLY   : MXUP,MXLW
      use module_time    ,only : MHD_steady
!
      implicit none
!
! --- [dummy arguments]
!
      logical,intent(in)    :: clrsld
      integer,intent(in)    :: NSSF,ndf,nq,IMODE,iter
      real*8 ,intent(in)    :: deltt,time,relax(0:MXMAT)
      character,intent(in)  :: cn
      integer,intent(in)    :: kdbf  (  MXCVFAC)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF(  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV(  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: rva   (  MXCVFAC)
      real*8 ,intent(inout) :: rvd   (  MXCVFAC,2)
      real*8 ,intent(in)    :: dif   (  MXALLCV,ndf)
      real*8 ,intent(in)    :: diag  (  MXALLCV)
      real*8 ,intent(in)    :: cofd  (  MXALLCV)
      real*8 ,intent(in)    :: dscl  (  MXALLCV)
      real*8 ,intent(in)    :: htcbnd(  NSSF)
      real*8 ,intent(inout) :: rhs   (  MXALLCV,nq)
      real*8 ,intent(in)    :: rho   (  MXALLCV,2)
      integer,intent(out)   :: ierror,iterph(nq)
      real*8 ,intent(out)   :: reps(nq),aeps(nq)
      integer,intent(in)    :: LCYCOLD(MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG (MXSSFBC_SLD)
      real*8 ,intent(in)    :: CVVOLM (    MXALLCV)
      real*8 ,intent(in)    :: SFCENT (3,  MXCVFAC)
!
! --- [local entities]
!
      character(20) :: cnx
!
      real*8 ,parameter   :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8 ,parameter   :: third=1.d0/3.d0
      real*8  :: dumC,wi1,wi2,dlvect
      real*8  :: dum1,dum2,dum3,dum4,th(2),costh1,sinth1,costh2,sinth2
      integer :: iti,icx,iterq,ierr,ierr1,kl,ku
      real*8  :: dx,dy,dz,alf,alfd,alfo,vlp,vlm,vdd,cof,alfn,dl
      real*8  :: wifa,wifb,vaa,vab,vbb,vba,dl1,dxo,dyo,dzo
      real*8  :: aepsq,repsq,calph,RMAXS
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf
      real*8  :: avsfo,dln,dl2,avsfn,alfdn,alfdo,alfe,alfen,
     &           alfeo,vlpN,vlmN,vlpO,vlmO
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD,ICH,IFLD,IBFL,IBFS,IBFE,IE,IV,ICFO,ICVBO
      integer :: NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ICFS,ICFE,ICFL
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICV,INLO,
     &           IIMATS(2),IMATS(2),ISLD,ISLD2,nl,nu,INLB
      integer :: IALLOCA=0,IQMXCV,ICV_D,ICV_A
      logical :: lg_IE
      real*8  :: org_x1,org_y1,org_x2,org_y2,
     &           cofA,cofB,alfA,alfB
!
      ipx(:)=0
      aax(:)=0.d0
!
      ierr=0
      IALLOCA=0
      ierror=0
      ierr1=0
      calph=1.0d0
      if(cn.eq.'v') then
        if(intgvv.eq.Adams_Moulton) then
          calph=5.0d0/12.d0
          if( iter.eq.1 ) calph=1.0d0
        elseif(intgvv.eq.Crank_Nicolson) then
          calph=0.5d0
          if( iter.eq.1 ) calph=1.0d0
        endif
      endif
!
      IQMXCV=MXCV
      if(.not.clrsld) IQMXCV=MXALLCV
      allocate(IQ(IQMXCV,2),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in solve_cnvdif'
        call FFRABORT(1,'solve_cnvdif')
      endif
!
!--------------------------
! --- re-allocate array  --
!--------------------------
      IALLOCA=0
 1000 continue
      lg_IE=.false.
      if(IALLOCA==1) then
        allocate (aax(MAXIE),
     &            ipx(MAXIE),
     &            IAL(MXALLCV,MAXIE),
     &            aalw(MXALLCV,0:MAXIE),
     &            stat=ierr1)
        if(ierr1.ne.0) then
          write(ifle,*) 
     &   'ERR: re-allocating error in solve_cnvdif'
          call FFRABORT(1,'solve_cnvdif')
      endif
      endif
!---------------------------------------------------------
!-< 1. Set maximul & minimul diffusivity for simplicity >-
!---------------------------------------------------------
      if(IMODE/=3) then
        do 100 IIMAT=1,NMAT   !ICV=1,NALLCV
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
          dx=dif(ICVL,1)
          dy=dx
          do m=2,ndf
            dx=max(dx,dif(ICVL,m))
            dy=min(dy,dif(ICVL,m))
          enddo
          dmax(ICVL)=dx
          dmin(ICVL)=dy
        enddo
!
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        do ICVL=IDCS,IDCE
        dx=dif(ICVL,1)
        dy=dx
        do 102 m=2,ndf
        dx=max(dx,dif(ICVL,m))
        dy=min(dy,dif(ICVL,m))
  102   continue
        dmax(ICVL)=dx
        dmin(ICVL)=dy
        enddo
  100   continue
      endif
!-----------------------------------
!-< 2. Make up coefficient matrix >-
!-----------------------------------
!--< 2.1 all over the domain >--    
!-----------------------------------

      if(ical_sld/=0) then
        call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!
      do 200 IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        IQ(ICVL,:)=0
        aalw(ICVL,0)=diag(ICVL)
        enddo
        do ICVL=ICVS,ICVE
        do 203 m=1,MAXIE
          aalw(ICVL,m)=0.d0
          IAL(ICVL,m)=0
 203    enddo
        enddo
!
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        aalw(IDCS:IDCE,1:MAXIE)=0.d0
        aalw(IDCS:IDCE,0)=diag(IDCS:IDCE)
        IAL(IDCS:IDCE,1:MAXIE)=0
  200 enddo
!
      if(IMODE==3) then  !
! --- VOF solver
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.99) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dumC=(wi1*rho(ICVA,1)+wi2*rho(ICVB,1))
          if(cn.eq.'a') then !dumC=1.d0
            vlp=max(0.d0,rva(ICFL))
            vlm=min(0.d0,rva(ICFL))
!            vlp=rva(ICFL)
!            vlm=rva(ICFL)
          else
            vlp=rva(ICFL)/dumC
            vlm=rva(ICFL)/dumC
          endif
!
          if(kdbf(ICFL).eq.0) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)-vlm
            aalw(ICVB,0)=aalw(ICVB,0)+vlp
            aalw(ICVA,IQ(ICVA,1))= vlm
            aalw(ICVB,IQ(ICVB,1))=-vlp
!
!            aalw(ICVA,0)=aalw(ICVA,0)+vlp 
!            aalw(ICVB,0)=aalw(ICVB,0)-vlm 
!            aalw(ICVA,IQ(ICVA,1))= +vlm   
!            aalw(ICVB,IQ(ICVB,1))= -vlp   
!
          else
           aalw(ICVA,0)=aalw(ICVA,0)-vlm
!           aalw(ICVA,0)=aalw(ICVA,0)+vlp
          endif
        endif
        enddo
        enddo
      elseif(IMODE==6) then
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.99) then
!
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
!
          unit(1,1)=CVCENT(1,ICVB)-SFCENT(1,ICFL)
          unit(2,1)=CVCENT(2,ICVB)-SFCENT(2,ICFL)
          unit(3,1)=CVCENT(3,ICVB)-SFCENT(3,ICFL)
          unit(1,2)=SFCENT(1,ICFL)-CVCENT(1,ICVA)
          unit(2,2)=SFCENT(2,ICFL)-CVCENT(2,ICVA)
          unit(3,2)=SFCENT(3,ICFL)-CVCENT(3,ICVA)
!
          th(1)=(unit(1,1)**2+unit(2,1)**2+unit(3,1)**2) 
          th(2)=(unit(1,2)**2+unit(2,2)**2+unit(3,2)**2) 
!
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA) 
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA) 
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA) 
!          dl=dx**2+dy**2+dz**2 
          dl=th(1)+th(2) 
          
!
          dlvect=(dx*SFAREA(1,ICFL)
     &           +dy*SFAREA(2,ICFL)
     &           +dz*SFAREA(3,ICFL))
!
!          alf=SFAREA(4,ICFL)/(dlvect+SML) 
          alf=dlvect*SFAREA(4,ICFL)/(dl+SML) 
! --- outlet: 
!
          if(kdbf(ICFL).eq.2) alf=0.d0
          alfA=alf*dif(ICVA,1)
          alfB=alf*dif(ICVB,1)
!
          vlp=max(0.d0,rva(ICFL))
          vlm=min(0.d0,rva(ICFL))
          cofA=cofd(ICVA)
          cofB=cofd(ICVB)
          if(kdbf(ICFL)==0) then
! --- regular face and periodic BC face 
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)

!
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
!
            if(IQ(ICVA,1)>MAXIE.or.IQ(ICVB,1)>MAXIE) then
              lg_IE=.true.
              cycle
            endif
!
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+(alfA-vlm)*calph
            aalw(ICVB,0)=aalw(ICVB,0)+(alfB+vlp)*calph
            
!
            aalw(ICVA,IQ(ICVA,1))=
     &                 (-alfA+vlm)*cofA*calph
!
            aalw(ICVB,IQ(ICVB,1))=
     &                 (-alfB-vlp)*cofB*calph
          else
!--- 1
            aalw(ICVA,0)=aalw(ICVA,0)+(alfA-vlm)*calph
          endif
        ENDIF
        ENDDO
        ENDDO
      elseif(IMODE==2) then
        do IIMAT=1,NMAT 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        
        if(kdbf(ICFL).lt.99) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dlvect=ABS(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          alf=SFAREA(4,ICFL)/(dlvect+SML)
!--------------
! --- outlet:
!--------------
          if(kdbf(ICFL).eq.2) alf=0.d0
!
          alfd=max(alf*max(dmax(ICVA),dmax(ICVB)), !15801 
     &                 dscl(ICVB)*SFAREA(4,ICFL))
          alfo=alf*min(dmin(ICVA),dmin(ICVB))
         
!
          vlp=max(0.d0,rva(ICFL))
          vlm=min(0.d0,rva(ICFL))
          vdd=max(rvd(ICFL,1),-rvd(ICFL,2))
          cof=cofd(ICVB)
          alx(ICVB)=alf

          if(kdbf(ICFL)==0) then
! --- regular face and periodic BC face
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            if(ICVB>IQMXCV) then
            endif
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
!
            if(IQ(ICVA,1)>MAXIE.or.IQ(ICVB,1)>MAXIE) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!1) Particle-MASS=OK; large-DT=DAME  (TANNO)
!            aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm+vdd)*calph  !OK
!            aalw(ICVB,0)=aalw(ICVB,0)+(alfd+vlp+vdd)*calph  !OK
!            aalw(ICVA,IQ(ICVA,1))=
!     &                 (-alfo-vlp+rvd(ICFL,2))*cof*calph
!            aalw(ICVB,IQ(ICVB,1))=
!     &                 (-alfo+vlm-rvd(ICFL,1))*cof*calph
!A) large-DT=OK
            aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm+vdd)*calph 
            aalw(ICVB,0)=aalw(ICVB,0)+(alfd+vlp+vdd)*calph 
            aalw(ICVA,IQ(ICVA,1))=
     &                 (-alfo+vlm+rvd(ICFL,2))*cof*calph
            aalw(ICVB,IQ(ICVB,1))=
     &                 (-alfo-vlp-rvd(ICFL,1))*cof*calph
!3) large-DT=OK (velocity is DAME)
!            aalw(ICVA,0)=aalw(ICVA,0)+(alfd+vlp+vdd)*calph 
!            aalw(ICVB,0)=aalw(ICVB,0)+(alfd-vlm+vdd)*calph 
!            aalw(ICVA,IQ(ICVA,1))=
!     &                 (-alfo+vlm+rvd(ICFL,2))*cof*calph
!            aalw(ICVB,IQ(ICVB,1))=
!     &                 (-alfo-vlp-rvd(ICFL,1))*cof*calph
!4) large-DT=DAME
!            aalw(ICVA,0)=aalw(ICVA,0)+(alfd+vlp+vdd)*calph 
!            aalw(ICVB,0)=aalw(ICVB,0)+(alfd-vlm+vdd)*calph 
!            aalw(ICVA,IQ(ICVA,1))=
!     &                 (-alfo-vlp+rvd(ICFL,2))*cof*calph
!            aalw(ICVB,IQ(ICVB,1))=
!     &                 (-alfo+vlm-rvd(ICFL,1))*cof*calph
!5) NEW
!            dum1=vlp-vlm
!            aalw(ICVA,0)=aalw(ICVA,0)+(alfd+dum1+vdd)*calph 
!            aalw(ICVB,0)=aalw(ICVB,0)+(alfd+dum1+vdd)*calph 
!            aalw(ICVA,IQ(ICVA,1))=
!     &                 (-alfo-dum1+rvd(ICFL,2))*cof*calph
!            aalw(ICVB,IQ(ICVB,1))=
!     &                 (-alfo-dum1-rvd(ICFL,1))*cof*calph
          else
!--- 1
!1)
             aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm+vdd)*calph  !OK
!A)
!             aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm+vdd)*calph  !OK      
!3)
!             aalw(ICVA,0)=aalw(ICVA,0)+(alfd+vlp+vdd)*calph
!4)
!             aalw(ICVA,0)=aalw(ICVA,0)+(alfd+vlp+vdd)*calph
          endif
        endif
        enddo
        enddo
      
      else
        do IIMAT=1,NMAT 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        
        if(kdbf(ICFL).lt.99) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dlvect=ABS(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          alf=SFAREA(4,ICFL)/(dlvect+SML)
!--------------
! --- outlet:
!--------------
          if(kdbf(ICFL).eq.2) alf=0.d0
!
          alfd=max(alf*max(dmax(ICVA),dmax(ICVB)), !15801 
     &                 dscl(ICVB)*SFAREA(4,ICFL))
          alfo=alf*min(dmin(ICVA),dmin(ICVB))
         
!
          vlp=max(0.d0,rva(ICFL))
          vlm=min(0.d0,rva(ICFL))
          vdd=max(rvd(ICFL,1),-rvd(ICFL,2))
          cof=cofd(ICVB)
          alx(ICVB)=alf

          if(kdbf(ICFL)==0) then
! --- regular face and periodic BC face
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            if(ICVB>IQMXCV) then
            endif
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
!
            if(IQ(ICVA,1)>MAXIE.or.IQ(ICVB,1)>MAXIE) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm+vdd)*calph  !OK
            aalw(ICVB,0)=aalw(ICVB,0)+(alfd+vlp+vdd)*calph  !OK
!
            aalw(ICVA,IQ(ICVA,1))=
!     &                 (-alfo-vlp+rvd(ICFL,2))*cof*calph
     &                 (-alfo+vlm+rvd(ICFL,2))*cof*calph
!
            aalw(ICVB,IQ(ICVB,1))=
!     &                 (-alfo+vlm-rvd(ICFL,1))*cof*calph
     &                 (-alfo-vlp-rvd(ICFL,1))*cof*calph

          else
!--- 1
            aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm+vdd)*calph  !OK
          endif
        endif
        enddo
        enddo
      endif
!-------------------
! --- check MAXIE 
!-------------------
      if(NPE>1) then
        call hpclor(lg_IE)
      endif
      if(lg_IE) then
          INLO=0
          do ICV=1,NCV
            INLB=IQ(ICV,1)
            if(INLB>INLO) then
              INLO=INLB
            endif
          enddo
          if(NPE>1) then
            call hpcimax(INLO)
          endif
          if(my_rank==root) then
            write(ifle,'(1x,a,I12,a,I12)')
     &         'MSG:MAX. IQ(ICV,1)= ',INLO,' MAXIE= ',MAXIE
            write(ifle,'(1x,a,a,I12,2a)')'MSG: ',
     &  'Increasing-1 [MAXIE] and [MAXIE] to :',INLO,' at ',trim(cn)
          endif
          MAXIE=INLO
          deallocate(aax,ipx,IAL,aalw,stat=ierr1)
          if(ierr1.ne.0) then
            write(ifle,*) 
     & 'ERR: deallocating error in solve_cnvdif'
            call FFRABORT(1,'MSG: Call your FFR supportor')
          endif
          IALLOCA=1
          goto 1000
      endif
!-----------------------------------------
!--< 2.2 Clear dummy cell & solid part >--
!-----------------------------------------
      IQ(NCV:,1)=0
      do 220 IDC=NCV+1,NALLCV
!      IQ(IDC,1)=0
      aalw(IDC,0)=1.d0
  220 continue
!
      if(clrsld.and.IMODE/=6) then
        do 225 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          DO 230 ICVL=ICVS,ICVE
          IQ(ICVL,1)=0
          aalw(ICVL,0)=1.d0
 230      continue          
        endif
 225    continue
      endif
!
!------------------------------
!--< 2.3 interface boundary >--
!------------------------------
!     
      if(.not.clrsld) then
        call bc_intrfcimp
     &   (NSSF,IQMXCV,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    alx,htcbnd,IQ,IAL,aalw)
      endif
!
!-------------------------
!--< 2.4 sliding mesh  >--
!-------------------------
!
!      if(cn/='R') then
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
        IIMATS=0
        IMATS=0
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdsld.and.idis(nb)==1) then
          do 400 ISLD=1,2
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
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          th(ISLD)=rot_ang(IMATS(ISLD))
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          th(ISLD2)=rot_ang(IMATS(ISLD2))
!
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2
     &              +unit(2,ISLD2)**2
     &              +unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
          org_x1=begin(1,IMATS(ISLD))*costh1
     &          -begin(2,IMATS(ISLD))*sinth1
          org_y1=begin(1,IMATS(ISLD))*sinth1
     &          +begin(2,IMATS(ISLD))*costh1
          org_x2=begin(1,IMATS(ISLD2))*costh2
     &          -begin(2,IMATS(ISLD2))*sinth2
          org_y2=begin(1,IMATS(ISLD2))*sinth2
     &          +begin(2,IMATS(ISLD2))*costh2
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(1,ICFP)
!
          ICVBO=LVEDGE(1,ICFO)
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          dum1=CVCENT(1,ICVA)*cos(th(ISLD))   !costh1
     &        -CVCENT(2,ICVA)*sin(th(ISLD))   !sinth1-org_x1
          dum2=CVCENT(1,ICVA)*sin(th(ISLD))   !sinth1
     &        +CVCENT(2,ICVA)*cos(th(ISLD))   !costh1-org_y1
          dum3=CVCENT(3,ICVA)
!
          dx  =CVCENT(1,ICVB)*cos(th(ISLD2))  !costh2
     &        -CVCENT(2,ICVB)*sin(th(ISLD2))  !sinth2-org_x2
          dy  =CVCENT(1,ICVB)*sin(th(ISLD2))  !sinth2
     &        +CVCENT(2,ICVB)*cos(th(ISLD2))  !costh2-org_y2
          dz  =CVCENT(3,ICVB)
!
          dx=dx-dum1
          dy=dy-dum2
          dz=dz-dum3
          dl1=dsqrt(dx*dx+dy*dy+dz*dz+SML)
!
          dxo =CVCENT(1,ICVBO)*cos(th(ISLD2))    !costh2
     &        -CVCENT(2,ICVBO)*sin(th(ISLD2))    !sinth2-org_x2
          dyo =CVCENT(1,ICVBO)*sin(th(ISLD2))    !sinth2
     &        +CVCENT(2,ICVBO)*cos(th(ISLD2))    !costh2-org_y2
          dzo =CVCENT(3,ICVBO)
!
          dxo=dxo-dum1
          dyo=dyo-dum2
          dzo=dzo-dum3
          dl2=dsqrt(dxo*dxo+dyo*dyo+dzo*dzo+SML)
!
          sf1=rbb(1,1,ISLD)*SFAREA(1,ICFL)
     &       +rbb(1,2,ISLD)*SFAREA(2,ICFL)
     &       +rbb(1,3,ISLD)*SFAREA(3,ICFL)
          sf2=rbb(2,1,ISLD)*SFAREA(1,ICFL)
     &       +rbb(2,2,ISLD)*SFAREA(2,ICFL)
     &       +rbb(2,3,ISLD)*SFAREA(3,ICFL)
          sf3=rbb(3,1,ISLD)*SFAREA(1,ICFL)
     &       +rbb(3,2,ISLD)*SFAREA(2,ICFL)
     &       +rbb(3,3,ISLD)*SFAREA(3,ICFL)
!
          avsf=SFAREA(4,ICFL)
          avsfn=wi1*avsf
          avsfo=wi2*avsf
!
          dx=dx/dl1
          dy=dy/dl1
          dz=dz/dl1
          dxo=dxo/dl2
          dyo=dyo/dl2
          dzo=dzo/dl2
!
          dln=abs(dx*sf1+dy*sf2+dz*sf3)+SML
!
          alf=1.d0/dl1  !(dln*dl1)
          alfn=alf*avsfn
          alfo=alf*avsfo
          alf =alfn+alfo
!
          alfd =(alf*max(dmax(ICVA),dmax(ICVB),dmax(ICVBO)))
          alfdn=max(alfn*max(dmax(ICVA),dmax(ICVB)) ,dscl(ICVB)*avsfn)
          alfdo=max(alfo*max(dmax(ICVA),dmax(ICVBO)),dscl(ICVB)*avsfo)
!
          alfe= alf*min(dmin(ICVA),dmin(ICVB),dmin(ICVBO))
          alfen=alf*min(dmin(ICVA),dmin(ICVB))
          alfeo=alf*min(dmin(ICVA),dmin(ICVBO))
!
          vlp=max(0.d0,rva(ICFL))
          vlm=min(0.d0,rva(ICFL))
!
          vlpN=wi1*max(0.d0,rva(ICFP))
          vlmN=wi1*min(0.d0,rva(ICFP))
          vlpO=wi2*max(0.d0,rva(ICFO))
          vlmO=wi2*min(0.d0,rva(ICFO))
!
          aalw(ICVA,0)=aalw(ICVA,0)+(-vlm)*calph       !1)
          aalw(ICVB,0)=aalw(ICVB,0)+(+vlpn)*calph     !2)
          aalw(ICVBO,0)=aalw(ICVBO,0)+(+vlpo)*calph   !3)
!
          if(IQ(ICVA,1)+1>MAXIE) cycle
          IQ(ICVA,1)=min(IQ(ICVA,1)+1,MAXIE)
          IAL(ICVA,IQ(ICVA,1))=ICVB
          aalw(ICVA,IQ(ICVA,1))=vlp  !    (+vlmn)*calph
!5
          if(IQ(ICVA,1)+1>MAXIE) cycle
          IQ(ICVA,1)=min(IQ(ICVA,1)+1,MAXIE)
          IAL(ICVA,IQ(ICVA,1))=ICVBO
          aalw(ICVA,IQ(ICVA,1))=vlm   !(+vlmo)*calph
!
          enddo
 400      enddo
          !alf
        elseif(kd==kdsld.and.idis(nb)==2) then
          do 500 ISLD=1,2
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
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          th(ISLD)=rot_ang(IMATS(ISLD))
!new
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
!          th(ISLD2)=rot_ang(IMATS(ISLD2))
!
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2+unit(2,ISLD2)**2+unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
!
!          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
!          do i=1,3
!          do j=1,3
!          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
!          enddo
!          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
          org_x1=begin(1,IMATS(ISLD))*costh1
     &          -begin(2,IMATS(ISLD))*sinth1
          org_y1=begin(1,IMATS(ISLD))*sinth1
     &          +begin(2,IMATS(ISLD))*costh1
          org_x2=begin(1,IMATS(ISLD2))*costh2
     &          -begin(2,IMATS(ISLD2))*sinth2
          org_y2=begin(1,IMATS(ISLD2))*sinth2
     &          +begin(2,IMATS(ISLD2))*costh2
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(1,ICFP)   
          
          ICVBO=LVEDGE(1,ICFO)
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          th(ISLD2)=OPPANG(IBFL)
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          dum1=CVCENT(1,ICVA)*cos(th(ISLD))   !costh1
     &        -CVCENT(2,ICVA)*sin(th(ISLD))   !sinth1-org_x1
          dum2=CVCENT(1,ICVA)*sin(th(ISLD))   !sinth1
     &        +CVCENT(2,ICVA)*cos(th(ISLD))   !costh1-org_y1
          dum3=CVCENT(3,ICVA)
!
          dx  =CVCENT(1,ICVB)*cos(th(ISLD2))  !costh2
     &        -CVCENT(2,ICVB)*sin(th(ISLD2))  !sinth2-org_x2
          dy  =CVCENT(1,ICVB)*sin(th(ISLD2))  !sinth2
     &        +CVCENT(2,ICVB)*cos(th(ISLD2))  !costh2-org_y2
          dz  =CVCENT(3,ICVB)
!
          dx=dx-dum1
          dy=dy-dum2
          dz=dz-dum3
          dl1=dsqrt(dx*dx+dy*dy+dz*dz+SML)
!
          dxo =CVCENT(1,ICVBO)*cos(th(ISLD2))    !costh2
     &        -CVCENT(2,ICVBO)*sin(th(ISLD2))    !sinth2-org_x2
          dyo =CVCENT(1,ICVBO)*sin(th(ISLD2))    !sinth2
     &        +CVCENT(2,ICVBO)*cos(th(ISLD2))    !costh2-org_y2
          dzo =CVCENT(3,ICVBO)
!
          dxo=dxo-dum1
          dyo=dyo-dum2
          dzo=dzo-dum3
          dl2=dsqrt(dxo*dxo+dyo*dyo+dzo*dzo+SML)
!
          sf1=rbb(1,1,ISLD)*SFAREA(1,ICFL)
     &       +rbb(1,2,ISLD)*SFAREA(2,ICFL)
     &       +rbb(1,3,ISLD)*SFAREA(3,ICFL)
          sf2=rbb(2,1,ISLD)*SFAREA(1,ICFL)
     &       +rbb(2,2,ISLD)*SFAREA(2,ICFL)
     &       +rbb(2,3,ISLD)*SFAREA(3,ICFL)
          sf3=rbb(3,1,ISLD)*SFAREA(1,ICFL)
     &       +rbb(3,2,ISLD)*SFAREA(2,ICFL)
     &       +rbb(3,3,ISLD)*SFAREA(3,ICFL)
!
          avsf=SFAREA(4,ICFL)
          avsfn=wi1*avsf
          avsfo=wi2*avsf
!
          dx=dx/dl1
          dy=dy/dl1
          dz=dz/dl1
          dxo=dxo/dl2
          dyo=dyo/dl2
          dzo=dzo/dl2
!
          dln=abs(dx*sf1+dy*sf2+dz*sf3)+SML
!
          alf=1.d0/(dln*dl1)
          alfn=alf*avsfn
          alfo=alf*avsfo
          alf =alfn+alfo
!
          alfd =(alf*max(dmax(ICVA),dmax(ICVB),dmax(ICVBO)))
          alfdn=max(alfn*max(dmax(ICVA),dmax(ICVB)) ,dscl(ICVB)*avsfn)
          alfdo=max(alfo*max(dmax(ICVA),dmax(ICVBO)),dscl(ICVB)*avsfo)
!
          alfe= alf*min(dmin(ICVA),dmin(ICVB),dmin(ICVBO))
          alfen=alf*min(dmin(ICVA),dmin(ICVB))
          alfeo=alf*min(dmin(ICVA),dmin(ICVBO))
!
          vlp=max(0.d0,rva(ICFL))
          vlm=min(0.d0,rva(ICFL))
!
          vlpN=wi1*max(0.d0,rva(ICFP))
          vlmN=wi1*min(0.d0,rva(ICFP))
          vlpO=wi2*max(0.d0,rva(ICFO))
          vlmO=wi2*min(0.d0,rva(ICFO))
!
          aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm)*calph       !1)
          aalw(ICVB,0)=aalw(ICVB,0)+(alfdn+vlpn)*calph     !2)
          aalw(ICVBO,0)=aalw(ICVBO,0)+(alfdo+vlpo)*calph   !3)
!
          if(IQ(ICVA,1)+1>MAXIE) cycle
          IQ(ICVA,1)=min(IQ(ICVA,1)+1,MAXIE)
          IAL(ICVA,IQ(ICVA,1))=ICVB
          aalw(ICVA,IQ(ICVA,1))=(-alfen+vlmn)*calph
!5
          if(IQ(ICVA,1)+1>MAXIE) cycle
          IQ(ICVA,1)=min(IQ(ICVA,1)+1,MAXIE)
          IAL(ICVA,IQ(ICVA,1))=ICVBO
          aalw(ICVA,IQ(ICVA,1))=(-alfeo+vlmo)*calph
!
          enddo
 500      enddo
          
        elseif(kd==kdovst) then
          cycle
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVP=LCYCSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
!
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
!
          alf=SFAREA(4,ICFL)/dum1
!
          alfd =(alf*max(dmax(ICVA),dmax(ICVB),dmax(ICVBO)))
          alfo=alf*min(dmin(ICVA),dmin(ICVB))
!
          vlp=max(0.d0,rva(ICFL))
          vlm=min(0.d0,rva(ICFL))
!
          ICVA=msk(ICVA)
          ICVB=msk(ICVB)
!
          aalw(ICVA,0)=aalw(ICVA,0)+(alfd-vlm)*calph       !1)
          aalw(ICVB,0)=aalw(ICVB,0)+(alfd+vlp)*calph     !2)
!
          if(IQ(ICVA,1)+1>MAXIE) cycle
          IQ(ICVA,1)=min(IQ(ICVA,1)+1,MAXIE)
          IAL(ICVA,IQ(ICVA,1))=ICVB
          aalw(ICVA,IQ(ICVA,1))=-alfo-vlp
!
          if(IQ(ICVB,1)+1>MAXIE) cycle
          IQ(ICVB,1)=min(IQ(ICVB,1)+1,MAXIE)
          IAL(ICVB,IQ(ICVB,1))=ICVA
          aalw(ICVB,IQ(ICVB,1))=-alfo+vlm
!
          enddo
          
        ENDIF
        enddo
      ENDIF
!      endif
!
!-------------------------------------------------------
!--< 2.4 rearrange order to split lower & upper part >--
!-------------------------------------------------------
        do 240 IIMAT=1,NMAT   !ICV=1,NALLCV  !8888
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 255 ICVL=ICVS,ICVE
        icx=IQ(ICVL,1)
        k=0
        do 241 m=1,icx
        if(IAL(ICVL,m).lt.ICVL.and.IAL(ICVL,m)/=0) then
! --- lower part: IQ(ICVL,1)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  241   continue
        IQ(ICVL,1)=k
        do 242 m=1,icx
        if(IAL(ICVL,m).gt.ICVL.and.IAL(ICVL,m)/=0) then
!        if(IAL(ICVL,m).gt.ICVL) then
! --- upper part: IQ(ICVL,2)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  242   continue
        IQ(ICVL,2)=k
        aalw(ICVL,1:)=0.d0
        IAL(ICVL,:)=0
        do 243 m=1,k
          aalw(ICVL,m)=aax(m)
          IAL(ICVL,m)=ipx(m)
  243   continue
!        if(k.ne.icx) then
!         write(*,'(1x,a,3I8)') 'ICVL= k=, icx= ',ICVL,k,icx
!         write(*,'(1x,a,4x,a)') 
!     &   '### program error -1- (solve_cnvdif)',cn 
!         CALL FFRABORT(1,'solve_cnvdif') 
!        endif
 255    continue
 240    enddo
!----------------------------
!-< 3. Solve linear system >-
!----------------------------
      bb(:)=0.d0
!
      if(i_steady==3) then
        do IIMAT=1,NMAT
         IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICVL=ICVS,ICVE
        aalw(ICVL,0)=aalw(ICVL,0)/relax(IIMAT)
        enddo
        enddo
      endif
!--------------------------------------------
! --- 
!--------------------------------------------
!
      do 300 m=1,nq
!
      do 310 IIMAT=1,NMAT    !ICV=1,NCV
      if(mat_cal(IIMAT)) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        bb(ICVL)=rhs(ICVL,m)
        enddo
      endif
  310 continue
!----------
! --- MHD
!----------
      if(IMODE==6.and.ical_MHD/=0.and.MHD_steady==1) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICVL=ICVS,ICVE
        aalw(ICVL,0)=aalw(ICVL,0)/relax_MHD(IIMAT)
        enddo
        enddo
      endif

! --- 
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,AALW(:,0))
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      if(IMODE==6) then
        aepsq=aepsbcg_MHD
        repsq=repsbcg_MHD
        iterq=iterbcg_MHD
      else
        aepsq=aepsbcg
        repsq=repsbcg
        iterq=iterbcg
      endif

!
!      if(IMODE==6) then
!        call utl_bcgstb(cn,m,
!     &     MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
!     &     NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
!     &     MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
!     &     bb,aepsq,repsq,iterq,ierr1,IMODE)
!      else
!
        if(ical_sld==0.and.ical_buff==0.and.ical_MHD==0) then 
          call utl_bcgstb(cn,m,
     &     MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &     NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &     MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &     bb,aepsq,repsq,iterq,ierr1,IMODE)
!
        else  !multi-material fluid
!
          call utl_bcgstb_MAT1(cn,
     &   MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aepsq,repsq,iterq,ierr1,IMODE)
        endif
!
!      endif
!
      IF(NPE.gt.1) THEN
        CALL hpcimax(iterq)
      ENDIF
!
      iterph(m)=iterq
      aeps(m)=aepsq
      reps(m)=repsq
!
      if(my_rank==root.and.mod(iter,100)==0) then
        if(iterq==0) then
          if(cn=='h') then
            write(ifle,'(2a)') ' ### WRN: BCGSTB solver NOT run at ',
     &      'enthalpy equation'
          elseif(cn=='k') then
            write(ifle,'(2a,I4,2x,a)') 
     &     ' ### WRN: BCGSTB solver NOT run at ',
     &      'scalar equation no= ',m,trim(sclname(m))
          elseif(cn=='y') then
            write(ifle,'(2a,I4,2x,a)') 
     &    ' ### WRN: BCGSTB solver NOT run at ',
     &      'species equation no= ',m,trim(spcnam(m))
          elseif(cn=='v') then
            write(ifle,'(2a,I4)') 
     &    ' ### WRN: BCGSTB solver NOT run at ',
     &      'velocity equation no= ',m
          endif
        write(ifle,'(2a)')
     &     ' ### MSG: Reduce [aepsbcg] and [repsbcg] in ',
     &  trim(cntlnam)
        endif
      endif
!--------------------------------------------
! --- 
!--------------------------------------------
      if(iterq/=0) then
        do 320 IIMAT=1,NMAT    !ICV=1,NCV
        if(mat_cal(IIMAT)) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          rhs(ICVL,m)=bb(ICVL)
          enddo
        endif
  320   enddo
      endif
!
      if(ierr1.ne.0) then
        call mknam(m,cn,iti,cnx)
        write(ifll,6000) cnx(:iti),iterq,repsq,aepsq,ierr1
        write(ifle,*) '### error : bicgstab not converged'
        goto 9999
      endif
!
  300 continue
!
!      if(.not.ical_vect) then
         deallocate(IQ)
!      endif
!
      return
 6000 format('*** bicg/',a,' : iterations=',i10,
     * ' / maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     * ' / return code=',i5)
!
 9999 continue
      write(ifle,*) '(solve_cnvdif)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine mknam(m,cn,j,cnx)
!=================================================
      integer     ,intent(in)  :: m
      character(*),intent(in)  :: cn
      integer     ,intent(out) :: j
      character(*),intent(out) :: cnx
      integer :: i
!
      cnx=cn
      if(nq.gt.1) write(cnx(11:),'(i10)') m
      j=0
      do 100 i=1,len(cnx)
      if(cnx(i:i).ne.' ' ) then
        j=j+1
        cnx(j:j)=cnx(i:i)
      endif
  100 continue
      end subroutine mknam
!
      end subroutine solve_cnvdif
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_cnvdif_e2p
     & (clrsld,NSSF,ndf,nq,cn,time,iph,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdbf,rva,dif,diag,
     &  LCYCOLD,wifsld,OPPANG,aks,
     &  cofd,dscl,htcbnd,rhs,deltt,
     &  iterph,reps,aeps,iter,ierror,IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine is for calculating NS equations
      use module_dimension
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IQ =>IW2K1
      use module_metrix,only   : bb
      use module_metrix,only   : aalw
      use module_metrix,only   : dmax=>W1K7
      use module_metrix,only   : dmin=>W1K8
      use module_metrix,only   : alx =>W1K9
      use module_material,only : ical_sld
      use module_model,only    : ical_vect,nthrds,ical_MHD
      use module_material,only : relax_MHD
!
! --- [module arguments]
!
      use module_hpcutil
      use module_io,only       : ifll,ifle,cntlnam
      use module_cgsolver,only : aepsbcg,repsbcg,iterbcg
      use module_flags,only    : 
     &                           Adams_Moulton,intgvv,Crank_Nicolson
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,ical_buff,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang
      use module_species,only  : spcnam
      use module_scalar ,only  : sclname
      use module_metrix,only   : aax,ipx
      use module_vector,ONLY   : MXUP,MXLW
      use module_time    ,only : MHD_steady
!
      implicit none
!
! --- [dummy arguments]
!
      logical,intent(in)    :: clrsld
      integer,intent(in)    :: NSSF,ndf,nq,IMODE,iter,iph
      real*8 ,intent(in)    :: deltt,time
      character,intent(in)  :: cn
      integer,intent(in)    :: kdbf  (  MXCVFAC)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF(  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV(  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: rva   (  MXCVFAC)
      real*8 ,intent(in)    :: dif   (  MXALLCV,ndf)
      real*8 ,intent(in)    :: diag  (  MXALLCV)
      real*8 ,intent(in)    :: cofd  (  MXALLCV)
      real*8 ,intent(in)    :: dscl  (  MXALLCV)
      real*8 ,intent(in)    :: htcbnd(  NSSF)
      real*8 ,intent(inout) :: rhs   (  MXALLCV,nq)
      real*8 ,intent(in)    :: rho   (  MXALLCV,2)
      integer,intent(out)   :: ierror,iterph(nq)
      real*8 ,intent(out)   :: reps(nq),aeps(nq)
      integer,intent(in)    :: LCYCOLD(MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG (MXSSFBC_SLD)
      real*8 ,intent(in)    :: aks    (MXALLCVR,mxrans) 
!
! --- [local entities]
!
      character(20) :: cnx
!
      real*8 ,parameter   :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0,
     &                       SML1=1.D-15
      real*8  :: dumC,wi1,wi2,dlvect
      real*8  :: dum1,dum2,dum3,dum4,th(2),costh1,sinth1,costh2,sinth2
      integer :: it,icx,iterq,ierr,ierr1,kl,ku
      real*8  :: dx,dy,dz,alf,alfd,alfo,vlp,vlm,vdd,cof,alfn,dl
      real*8  :: wifa,wifb,vaa,vab,vbb,vba,dl1,dxo,dyo,dzo
      real*8  :: aepsq,repsq,calph,RMAXS
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf
      real*8  :: avsfo,dln,dl2,avsfn,alfdn,alfdo,alfe,alfen,
     &           alfeo,vlpN,vlmN,vlpO,vlmO
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD,ICH,IFLD,IBFL,IBFS,IBFE,IE,IV,ICFO,ICVBO
      integer :: NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ICFS,ICFE,ICFL
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICV,INLO,
     &           IIMATS(2),IMATS(2),ISLD,ISLD2,nl,nu,INLB
      integer :: IALLOCA=0,IQMXCV
      logical :: lg_IE
      real*8  :: org_x1,org_y1,org_x2,org_y2,
     &           cofA,cofB,alfA,alfB
!
      ipx(:)=0
      aax(:)=0.d0
!
      ierr=0
      IALLOCA=0
      ierror=0
      ierr1=0
      calph=1.0d0
!
      IQMXCV=MXCV
      if(.not.clrsld) IQMXCV=MXALLCV
      allocate(IQ(IQMXCV,2),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in solve_cnvdif_e2p'
        call FFRABORT(1,'solve_cnvdif_e2p')
      endif
!
!---------------------------------------------------------
!-< 1. Set maximul & minimul diffusivity for simplicity >-
!---------------------------------------------------------
      do 100 IIMAT=1,NMAT   !ICV=1,NALLCV
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
          dx=dif(ICVL,1)
          dy=dx
          do m=2,ndf
            dx=max(dx,dif(ICVL,m))
            dy=min(dy,dif(ICVL,m))
          enddo
          dmax(ICVL)=dx
          dmin(ICVL)=dy
        enddo
!
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        do ICVL=IDCS,IDCE
        dx=dif(ICVL,1)
        dy=dx
        do 102 m=2,ndf
        dx=max(dx,dif(ICVL,m))
        dy=min(dy,dif(ICVL,m))
  102   continue
        dmax(ICVL)=dx
        dmin(ICVL)=dy
        enddo
  100 continue
!-----------------------------------
!-< 2. Make up coefficient matrix >-
!-----------------------------------
!--< 2.1 all over the domain >--
!-----------------------------------
      if(ical_sld/=0) then
      call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!
      do 200 IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        IQ(ICVL,:)=0
        aalw(ICVL,0)=diag(ICVL)
        enddo
        do ICVL=ICVS,ICVE
        do 203 m=1,MAXIE
          aalw(ICVL,m)=0.d0
          IAL(ICVL,m)=0
 203    enddo
        enddo
!
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        aalw(IDCS:IDCE,1:MAXIE)=0.d0
        aalw(IDCS:IDCE,0)=diag(IDCS:IDCE)
        IAL(IDCS:IDCE,1:MAXIE)=0
  200 enddo
!
      do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        
        if(kdbf(ICFL).lt.99) then
          wi1=wiface(ICFL)     !wiface(ICFL)
          wi2=1.d0-wi1
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dlvect=ABS(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          dum2=aks(ICVA,iph)
          dum3=aks(ICVB,iph)
          dum4=wi1*dum2+wi2*dum3
          alf=SFAREA(4,ICFL)/(dlvect+SML)
          if(kdbf(ICFL).eq.2) alf=0.d0
          alfd=dum4*max(alf*max(dmax(ICVA),dmax(ICVB)),
     &                 dscl(ICVB)*SFAREA(4,ICFL))
          alfo=dum4*alf*min(dmin(ICVA),dmin(ICVB))

          vlp=max(0.d0,rva(ICFL))*dum4
          vlm=min(0.d0,rva(ICFL))*dum4

!          alfd=max(alf*max(dmax(ICVA),dmax(ICVB)),
!     &                 dscl(ICVB)*SFAREA(4,ICFL))
!          alfo=alf*min(dmin(ICVA),dmin(ICVB))

!          vlp=max(0.d0,rva(ICFL))
!          vlm=min(0.d0,rva(ICFL))

          
          cof=cofd(ICVB)
          alx(ICVB)=alf
          if(kdbf(ICFL)==0) then
! --- regular face and periodic BC face
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
!
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+(alfd+vlp)!*dum2
            aalw(ICVB,0)=aalw(ICVB,0)+(alfd-vlm)!*dum3
            aalw(ICVA,IQ(ICVA,1))=(-alfo+vlm)*cof!*dum2
            aalw(ICVB,IQ(ICVB,1))=(-alfo-vlp)*cof!*dum3
          else
!--- 1
            aalw(ICVA,0)=aalw(ICVA,0)+(alfd+vlp)!*dum2
          endif
        endif
        enddo
      enddo
!-------------------
! --- check MAXIE
!-------------------
      if(NPE>1) then
        call hpclor(lg_IE)
      endif
!-----------------------------------------
!--< 2.2 Clear dummy cell & solid part >--
!-----------------------------------------
      IQ(NCV:,1)=0
      do 220 IDC=NCV+1,NALLCV
!      IQ(IDC,1)=0
      aalw(IDC,0)=1.d0
  220 continue
!
      do 225 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          DO 230 ICVL=ICVS,ICVE
          IQ(ICVL,1)=0
          aalw(ICVL,0)=1.d0
 230      continue          
        endif
 225  continue
!------------------------------
!--< 2.3 interface boundary >--
!------------------------------
!     
      if(.not.clrsld) then
        call bc_intrfcimp
     &   (NSSF,IQMXCV,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    alx,htcbnd,IQ,IAL,aalw)
      endif
!--------------------------------------------------
! --- rearrange order to split lower & upper part 
!--------------------------------------------------
      do 240 IIMAT=1,NMAT   !ICV=1,NALLCV
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 255 ICVL=ICVS,ICVE
        icx=IQ(ICVL,1)
        k=0
        do 241 m=1,icx
        if(IAL(ICVL,m).lt.ICVL.and.IAL(ICVL,m)/=0) then
! --- lower part: IQ(ICVL,1)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  241   enddo
        IQ(ICVL,1)=k
        do 242 m=1,icx
        if(IAL(ICVL,m).gt.ICVL.and.IAL(ICVL,m)/=0) then
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  242   enddo
        IQ(ICVL,2)=k
        aalw(ICVL,1:)=0.d0
        IAL(ICVL,:)=0
        do 243 m=1,k
          aalw(ICVL,m)=aax(m)
          IAL(ICVL,m)=ipx(m)
  243   enddo
 255    enddo
 240  enddo
!----------------------------
!-< 3. Solve linear system
!----------------------------
      bb(:)=0.d0
      do 300 m=1,nq
      do 310 IIMAT=1,NMAT    !ICV=1,NCV
      if(mat_cal(IIMAT)) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        bb(ICVL)=rhs(ICVL,m)
        enddo
      endif
  310 continue
!----------
! --- MHD
!----------
      do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICVL=ICVS,ICVE
        if(aalw(ICVL,0)<SML1) then
           aalw(ICVL,0)=1.D0
           bb(ICVL)=0.D0
        endif
        enddo
      enddo
!
! --- 
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,AALW(:,0))
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      aepsq=aepsbcg
      repsq=repsbcg
      iterq=iterbcg
!
      if(ical_sld==0.and.ical_buff==0) then  
!
        call utl_bcgstb(cn,m,
     &     MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &     NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &     MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &     bb,aepsq,repsq,iterq,ierr1,IMODE)
!
      else  !multi-material fluid

        call utl_bcgstb_MAT1(cn,
     &   MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aepsq,repsq,iterq,ierr1,IMODE)

      endif
!
      IF(NPE.gt.1) THEN
        CALL hpcimax(iterq)
      ENDIF
!
      iterph(m)=iterq
      aeps(m)=aepsq
      reps(m)=repsq
!
      if(my_rank==root.and.mod(iter,100)==0) then
        if(iterq==0) then
          if(cn=='h') then
            write(ifle,'(2a)') ' ### WRN: BCGSTB solver NOT run at ',
     &      'enthalpy equation'
          elseif(cn=='k') then
            write(ifle,'(2a,I4,2x,a)') 
     &     ' ### WRN: BCGSTB solver NOT run at ',
     &      'scalar equation no= ',m,trim(sclname(m))
          elseif(cn=='y') then
            write(ifle,'(2a,I4,2x,a)') 
     &    ' ### WRN: BCGSTB solver NOT run at ',
     &      'species equation no= ',m,trim(spcnam(m))
          elseif(cn=='v') then
            write(ifle,'(2a,I4)') 
     &    ' ### WRN: BCGSTB solver NOT run at ',
     &      'velocity equation no= ',m
          endif
        write(ifle,'(2a)')
     &     ' ### MSG: Reduce [aepsbcg] and [repsbcg] in ',
     &  trim(cntlnam)
        endif
      endif
!
! --- 
!
      if(iterq/=0) then
        do 320 IIMAT=1,NMAT    !ICV=1,NCV
        if(mat_cal(IIMAT)) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          rhs(ICVL,m)=bb(ICVL)
          enddo
        endif
  320   continue
      endif
!
      if(ierr1.ne.0) then
        write(ifll,6000) cnx(:it),iterq,repsq,aepsq,ierr1
        write(ifle,*) '### error : bicgstab not converged'
        goto 9999
      endif
!
  300 continue
!
       deallocate(IQ)
!
      return
 6000 format('*** bicg/',a,' : iterations=',i10,
     * ' / maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     * ' / return code=',i5)
!
 9999 continue
      write(ifle,*) '(solve_cnvdif_e2p)'
      ierror=1
!
      end subroutine solve_cnvdif_e2p
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_poisson_unsymm_e2p
     & (iter,nflud,deltt,time,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  ipfix,kdbf,diag,rhs,rho,rho2,aks,rva,ccc,
     &  LCYCOLD,wifsld,
     &  iterp,repsp,aepsp,IMODE,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine is for calculating NS equations
      use module_dimension
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IQ =>IW2K1
      use module_metrix,only   : bb
      use module_metrix,only   : aalw
      use module_material,only : ical_sld
      
!
! --- [module arguments]
!
      use module_hpcutil
      use module_io,only       : ifll,ifle,cntlnam
      use module_initial ,only : rho0,rho02
      use module_cgsolver,only : aepscg,repscg,itercg
      use module_flags,only    :
     &                           Adams_Moulton,intgvv,Crank_Nicolson
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,kdprdc,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang
      use module_species,only  : spcnam,gascns
      use module_scalar ,only  : sclname,ivold,ical_s,iaph
      use module_metrix,only   : aax,ipx
      use module_model,only    : ical_vect
      use module_scalar,  only : iaph,ivof
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nflud,iter,IMODE
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: ipfix  (    MXMAT) 
      integer,intent(in)    :: kdbf   (  MXCVFAC)
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: diag   (  MXALLCV)
      real*8 ,intent(inout) :: rhs    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      integer,intent(out)   :: ierror,iterp
      real*8 ,intent(out)   :: repsp,aepsp
      real*8 ,intent(inout) :: rho    (  MXALLCV)
      real*8 ,intent(in)    :: ccc    (  MXALLCV)
!
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: rho2   (  MXALLCVC)
      real*8 ,intent(in)    :: aks    (  MXALLCVR,mxrans)
      real*8 ,intent(inout) :: rva    (  MXCVFAC)
!
! --- [local entities]
!
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8  :: dumC,wi1,wi2,dlvect
      real*8  :: dum1,dum2,dum3,dum4,th(2),costh1,sinth1,costh2,sinth2
      integer :: it,icx,iterq,ierr,ierr1
      real*8  :: dx,dy,dz,alf,alfd,alfo,vlp,vlm,vdd,cof,alfn,dl
      real*8  :: wifa,wifb,vaa,vab,vbb,vba,dl1,dxo,dyo,dzo
      real*8  :: aepsq,repsq,calph,RMAXS
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf
      real*8  :: avsfo,dln,dl2,dl3,avsfn,alfdn,alfdo,alfe,alfen,
     &           alfA,alfB,vlp1,vlm1,
     &           avsfn1,avsfo1,aph1A,aph1B,aph2A,aph2B,dumA,dumB
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,IPR
      integer :: IBFL,IBFS,IBFE,ICFO,ICVBO
      integer :: NB,KD,IDC,IV
      integer :: ICVA,ICVB,IVA,IVB,IBFP,ICFP,ICVP,IDCP,ICVAO
      integer :: ICFS,ICFE,ICFL
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICV,INLO,
     &           IIMATS(2),IMATS(2),ISLD,ISLD2,nl,nu,INLB,IA,IAA,IAO,IAB
      logical :: lg_IE
      integer :: IALLOCA=0
      character(20) :: cn
      integer :: IQMXCV,ICV_D,ICV_A
!
      ierr=0
      ierror=0
      ierr1=0
      IALLOCA=0
!
      IQMXCV=MXCV
!      
      allocate(IQ(IQMXCV,2),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error solve_poisson_unsymm_e2p'
        call FFRABORT(1,'solve_poisson_unsymm_e2p')
      endif
!---------------------------------------------------------
!-< 1. Set maximul & minimul diffusivity for simplicity >-
!---------------------------------------------------------
      IALLOCA=0
 1000 continue
      lg_IE=.false.

      if(IALLOCA==1) then
        allocate (aax(MAXIE),
     &            ipx(MAXIE),
     &            IAL(MXALLCV,MAXIE),
     &            aalw(MXALLCV,0:MAXIE),
     &            stat=ierr1)
        if(ierr1.ne.0) then
          write(ifle,*) 
     &         'ERR: re-allocating error in solve_poisson_unsymm_e2p'
          call FFRABORT(1,'solve_poisson_unsymm_e2p')
      endif
      endif
!-----------------------------------
!-< 2. Make up coefficient matrix >-
!-----------------------------------
!--< 2.1 all over the domain >--
!-----------------------------------
      if(ical_sld/=0) then
        call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!
      IQ(:,:)=0
      aalw(:,0)=0.d0
      do 200 IIMAT=1,NMAT   !ICV=1,NALLCV
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      aalw(ICVL,0)=diag(ICVL)
      enddo
  200 enddo
!
      if(IMODE==6) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
!          aph1A=1.d0-aks(ICVA,ivold)
!          aph1B=1.d0-aks(ICVB,ivold)
!          aph2A=aks(ICVA,ivold)
!          aph2B=aks(ICVB,ivold)
!          dumA=aph1A/rho0(IIMAT)+aph2A/rho2(ICVA)
!          dumB=aph1B/rho0(IIMAT)+aph2B/rho2(ICVB)
          dumA=rho(ICVA)
          dumB=rho(ICVB)
          dum4=wi1*dumA+wi2*dumB

          dum1=rho(ICVA)*ccc(ICVA)*rho(ICVA)
          dum2=rho(ICVB)*ccc(ICVB)*rho(ICVB)
          dum3=wi1*dum1+wi2*dum2


          vlp1=rva(ICFL)/dum3  !cavit
          vlm1=rva(ICFL)/dum3
          vlp=max(0.d0,vlp1)
          vlm=min(0.d0,vlm1)

          alf=SFAREA(4,ICFL)/(dlvect)*deltt/dum4  !cavit
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE .or. IQ(ICVB,1)>MAXIE ) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf-vlm!+vlp
            aalw(ICVB,0)=aalw(ICVB,0)+alf+vlp!-vlm
!
            aalw(ICVA,IQ(ICVA,1))=-alf+vlm
            aalw(ICVB,IQ(ICVB,1))=-alf-vlp
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf-vlm!+vlp
          endif
        endif
        enddo
        enddo
      elseif(IMODE==8) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          aph1A=aks(ICVA,ical_s)
          aph1B=aks(ICVB,ical_s)
          aph2A=1.d0-aks(ICVA,ical_s)
          aph2B=1.d0-aks(ICVB,ical_s)
          dumA=aph1A/rho02(IIMAT)+aph2A/rho2(ICVA)
          dumB=aph1B/rho02(IIMAT)+aph2B/rho2(ICVB)
          dum4=wi1*dumA+wi2*dumB
          alf=SFAREA(4,ICFL)/(dlvect)*deltt*dum4
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE .or. IQ(ICVB,1)>MAXIE ) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            aalw(ICVA,IQ(ICVA,1))=-alf
            aalw(ICVB,IQ(ICVB,1))=-alf
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
        enddo
        enddo
      elseif(IMODE==9) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          aph1A=aks(ICVA,ivof)
          aph1B=aks(ICVB,ivof)
          aph2A=1.d0-aks(ICVA,ivof)
          aph2B=1.d0-aks(ICVB,ivof)
          dumA=aph1A/rho0(IIMAT)+aph2A/rho02(IIMAT)
          dumB=aph1B/rho0(IIMAT)+aph2B/rho02(IIMAT)
          dum4=wi1*dumA+wi2*dumB
!
!          if(rva(ICFL)>0.d0) then
!            ICV_D=ICVA
!            ICV_A=ICVB
!          else
!            ICV_D=ICVB
!            ICV_A=ICVA
!          endif
!          dum1=aks(ICV_D,ivof)
!          if(dum1<SML) then
!            dum4=rho02(IIMAT)
!          elseif((1.d0-dum1)<SML) then
!            dum4=rho0(IIMAT)
!          else
!            dum4=rho0(IIMAT)*aks(ICV_A,ivof)
!     &         +rho02(IIMAT)*(1.d0-aks(ICV_A,ivof))
!          endif
!          
          alf=SFAREA(4,ICFL)/(dlvect)*deltt*dum4
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE.or.IQ(ICVB,1)>MAXIE) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            aalw(ICVA,IQ(ICVA,1))=-alf
            aalw(ICVB,IQ(ICVB,1))=-alf
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
        enddo
        enddo
      elseif(IMODE==7) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          dumA=rho(ICVA)
          dumB=rho(ICVB)
          dum4=wi1*dumA+wi2*dumB
!
          alf=SFAREA(4,ICFL)/dlvect*dum4
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE.or.IQ(ICVB,1)>MAXIE) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            aalw(ICVA,IQ(ICVA,1))=-alf
            aalw(ICVB,IQ(ICVB,1))=-alf
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
        enddo
        enddo
      elseif(IMODE==5) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          aph1A=aks(ICVA,iaph(1))
          aph1B=aks(ICVB,iaph(1))
          aph2A=aks(ICVA,iaph(2))
          aph2B=aks(ICVB,iaph(2))
          dumA=aph1A/rho(ICVA)+aph2A/rho2(ICVA)
          dumB=aph1B/rho(ICVB)+aph2B/rho2(ICVB)

!          dumA=aph1A+aph2A
!          dumB=aph1B+aph2B
          dum4=wi1*dumA+wi2*dumB
          alf=deltt*SFAREA(4,ICFL)/(dlvect)*dum4!*deltt
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE.or.IQ(ICVB,1)>MAXIE ) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf!*dumA
            aalw(ICVB,0)=aalw(ICVB,0)+alf!*dumB
!!
            aalw(ICVA,IQ(ICVA,1))=-alf!*dumB
            aalw(ICVB,IQ(ICVB,1))=-alf!*dumA
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf!*dumA
          endif
        endif
        enddo
        enddo
!
      endif


      if(.false.) then
      if(IMODE==10) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          alf=SFAREA(4,ICFL)/(dlvect)*deltt
          dum1=rho(ICVA)*ccc(ICVA)
          dum2=rho(ICVB)*ccc(ICVB)
          dum3=wi1/dum1+wi2/dum2
          
          vlp1=rva(ICFL)*dum3
          vlm1=rva(ICFL)*dum3
!
          vlp=max(0.d0,vlp1)
          vlm=min(0.d0,vlm1)

          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE.or.IQ(ICVB,1)>MAXIE ) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf+vlp
            aalw(ICVB,0)=aalw(ICVB,0)+alf-vlm
!
            aalw(ICVA,IQ(ICVA,1))=-alf+vlm 
            aalw(ICVB,IQ(ICVB,1))=-alf-vlp 
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf+vlp
          endif
        endif
        enddo
        enddo

      else
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          aph1A=aks(ICVA,iaph(1))
          aph1B=aks(ICVB,iaph(1))
          aph2A=aks(ICVA,iaph(2))
          aph2B=aks(ICVB,iaph(2))
          dumA=aph1A/rho(ICVA)+aph2A/rho2(ICVA)
          dumB=aph1B/rho(ICVB)+aph2B/rho2(ICVB)
          dum4=wi1*dumA+wi2*dumB
          alf=dum4*SFAREA(4,ICFL)/(dlvect)*deltt
!          alfA=dumA*SFAREA(4,ICFL)/(dlvect)*deltt
!          alfB=dumB*SFAREA(4,ICFL)/(dlvect)*deltt
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE .or. IQ(ICVB,1)>MAXIE ) then
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            aalw(ICVA,IQ(ICVA,1))=-alf
            aalw(ICVB,IQ(ICVB,1))=-alf
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
        enddo
        enddo
      endif
      endif
!-------------------
! --- check MAXIE
!-------------------
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
        IIMATS=0
        IMATS=0
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdsld) then
          do 400 ISLD=1,2
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
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          th(ISLD)=rot_ang(IMATS(ISLD))
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          th(ISLD2)=rot_ang(IMATS(ISLD2))
!
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2+unit(2,ISLD2)**2+unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          sf1=SFAREA(1,ICFL)
          sf2=SFAREA(2,ICFL)
          sf3=SFAREA(3,ICFL)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl1=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl1=deltt/dl1
!
          IVA=LVEDGE(1,ICFP)
          IVB=LVEDGE(2,ICFP)
          dx=CVCENT(1,IVB)-CVCENT(1,IVA)
          dy=CVCENT(2,IVB)-CVCENT(2,IVA)
          dz=CVCENT(3,IVB)-CVCENT(3,IVA)
          sf1=SFAREA(1,ICFP)
          sf2=SFAREA(2,ICFP)
          sf3=SFAREA(3,ICFP)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl2=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
!          dl2=deltt/dl2
!
          ICVAO=LVEDGE(1,ICFO)
          ICVBO=LVEDGE(2,ICFO)
          dx=CVCENT(1,ICVBO)-CVCENT(1,ICVAO)
          dy=CVCENT(2,ICVBO)-CVCENT(2,ICVAO)
          dz=CVCENT(3,ICVBO)-CVCENT(3,ICVAO)
          sf1=SFAREA(1,ICFO)
          sf2=SFAREA(2,ICFO)
          sf3=SFAREA(3,ICFO)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl3=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
!          dl3=deltt/dl3
!
          avsf =dl1*SFAREA(4,ICFL)
          avsfn=dl1*SFAREA(4,ICFL)*wi1
          avsfo=dl1*SFAREA(4,ICFL)*wi2
!
          IQ(ICVA,1)=IQ(ICVA,1)+1      !1)
          IAA=IQ(ICVA,1)
!
          IQ(IVA,1)=IQ(IVA,1)+1        !2)
          IA=IQ(IVA,1)
!
          IQ(ICVA,1)=IQ(ICVA,1)+1
          IAB=IQ(ICVA,1)
!
          IQ(ICVAO,1)=IQ(ICVAO,1)+1    !3)
          IAO=IQ(ICVAO,1)
!
!          aalw(ICVA,0)=aalw(ICVA,0)+avsf
          if(IAA<=MAXIE.and.IA<=MAXIE) then
            IAL(ICVA,IAA)=IVA
            aalw(ICVA,IAA)=-avsfn
            IAL(IVA,IA)=ICVA
            aalw(IVA,IA)=-avsfn
            aalw(ICVA,0)=aalw(ICVA,0)+avsf
            aalw(IVA,0)=aalw(IVA,0)+avsf
          else
            lg_IE=.true.
          endif
!
          if(IAB<=MAXIE.and.IAO<=MAXIE) then
            IAL(ICVA,IAB)=ICVAO
            aalw(ICVA,IAB)=-avsfo
            IAL(ICVAO,IAO)=ICVA
            aalw(ICVAO,IAO)=-avsfo
            aalw(ICVA,0)=aalw(ICVA,0)+avsf
            aalw(ICVAO,0)=aalw(ICVAO,0)+avsf
          else
            lg_IE=.true.
          endif
!
         enddo
 400     enddo
        ENDIF
        enddo
!
!---------------------------
! --- Coef. of Ap>=SIGMA(An)
!---------------------------
        do ICV=1,NCV
        RMAXS=0.d0
        do IV=1,IQ(ICV,1)
          RMAXS=RMAXS-aalw(ICV,IV)
        enddo
        if(aalw(ICV,0)<RMAXS) then
          aalw(ICV,0)=RMAXS
        endif
        enddo
      endif
      if(NPE>1) then
        call hpclor(lg_IE)
      endif
      if(lg_IE) then
          INLO=0
          do ICV=1,NCV
            INLB=IQ(ICV,1)
            if(INLB>INLO) then
              INLO=INLB
            endif
          enddo
          if(NPE>1) then
            call hpcimax(INLO)
          endif
          if(my_rank==root) then
            write(ifle,'(1x,a,I12,a,I12)')
     &         'MSG:MAX. IQ(ICV,1)= ',INLO,' MAXIE= ',MAXIE
            write(ifle,'(1x,a,a,I12,2a)')'MSG: ',
     &    'Increasing-2 [MAXIE] and [MAXIE] to :',INLO,' at ',trim(cn)
          endif
          IEMAX=INLO
          MAXIE=INLO
          deallocate(aax,ipx,IAL,aalw,stat=ierr1)
          if(ierr1.ne.0) then
            write(ifle,*) 
     &           'ERR: deallocating error in solve_poisson_unsymm_e2p'
            call FFRABORT(1,'MSG: Call your FFR supportor')
          endif
          IALLOCA=1
          goto 1000
      endif
!-----------------------------------------
!
!
!-------------------------
!--< 2.4 sliding mesh  >--
!-------------------------
!
!-------------------------------------------------------
!--< 2.4 rearrange order to split lower & upper part >--
!-------------------------------------------------------
      do 240 IIMAT=1,NMAT   !ICV=1,NALLCV
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 255 ICVL=ICVS,ICVE
        icx=IQ(ICVL,1)
        k=0
        do 241 m=1,icx
        if(IAL(ICVL,m).lt.ICVL) then
! --- lower part: IQ(ICVL,1)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  241   enddo
        IQ(ICVL,1)=k
        do 242 m=1,icx
        if(IAL(ICVL,m).gt.ICVL) then
! --- upper part: IQ(ICVL,2)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  242   continue
        IQ(ICVL,2)=k
        do 243 m=1,k
          aalw(ICVL,m)=aax(m)
          IAL(ICVL,m)=ipx(m)
  243   continue
        if( k.ne.icx ) then
          write(*,*) '### program error -1- (solve_poisson_unsymm_e2p)'
          CALL FFRABORT(1,'solve_poisson_unsymm_e2p')
        endif
 255    continue
 240  enddo
!----------------------------
!-< 3. Solve linear system >-
!----------------------------
      bb(:)=0.d0
!
      do 310 IIMAT=1,NMAT    !ICV=1,NCV
      if(mat_cal(IIMAT)) then
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        bb(ICVL)=rhs(ICVL)
        enddo
      endif
  310 enddo
!
!-----------------------------------------
!--< 2.2 Clear dummy cell & solid part >--
!-----------------------------------------
!
      IQ(NCV+1:IQMXCV,1)=0
      do 220 IDC=NCV+1,IQMXCV
      aalw(IDC,0)=1.d0
  220 continue
!
      do 225 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0.and.IMODE/=7) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          DO 230 ICVL=ICVS,ICVE
          IQ(ICVL,1)=0
          aalw(ICVL,0)=1.d0
          bb(ICVL)=0.d0
 230      continue          
        elseif(IMODE==7) then
!          ICVS=MAT_CVEXT(IIMAT-1)+1
!          ICVE=MAT_CVEXT(IIMAT)
!          DO ICVL=ICVS,ICVE
!          dum1=0.d0
!          do I=1,MAXIE
!          dum1=dum1+aalw(ICVL,I)
!          enddo
!          aalw(ICVL,0)=abs(dum1)
!          enddo
       endif
 225  continue
      

! --- 
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,AALW(:,0))
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      aepsq=aepscg
      repsq=repscg
      iterq=itercg
!
      m=1
      cn='p'
!      cn='bicg_poisson'
      if(NMAT==1) then
        call utl_bcgstb(cn,m,
     &   MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aepsq,repsq,iterq,ierr1,IMODE)
      else
        call utl_bcgstb_MAT1(cn,
     &   MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aepsq,repsq,iterq,ierr1,IMODE)
      endif
!
      IF(NPE.gt.1) THEN
        CALL hpcimax(iterq)
      ENDIF
!
      aepsp=aepsq
      repsp=repsq
      iterp=iterq
!
      if(my_rank==root.and.mod(iter,100)==0) then
        if(iterq==0) then
          write(ifle,'(2a)') ' ### WRN: BCGSTB solver NOT run at ',
     &   'unsymm poisson'
        endif
      endif
!
! --- 
!
      if(iterq/=0) then
        do 320 IIMAT=1,NMAT    !ICV=1,NCV
        if(mat_cal(IIMAT)) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          rhs(ICVL)=bb(ICVL)
          enddo
        endif
  320   continue
      endif
!
      deallocate(IQ)
!
      return
!
      write(ifll,6000) cn,iterq,repsq,aepsq,ierr1
 6000 format('ERR: bicg_poisson/',a,' : iterations=',i10,
     * ' / maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     * ' / return code=',i5)
!
 9999 continue
      write(ifle,*) '(solve_poisson_unsymm_e2p)'
      ierror=1
!
      end subroutine solve_poisson_unsymm_e2p
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_poisson 
     & (iter,nflud,deltt,time,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,SFCENT,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  ipfix,kdbf,diag,rhs,rho,
     &  LCYCOLD,wifsld,OPPANG,
     &  iterp,repsp,aepsp,IMODE,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- CG SOLVER
!
      use module_metrix,only : msk
      use module_metrix,only : IAL
      use module_metrix,only : IAU=>IW2K2
      use module_metrix,only : INL
      use module_metrix,only : INU
      use module_metrix,only : aalw
      use module_metrix,only : bb
      use module_metrix,only : aax,ipx
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_io,only       : ifll,ifle
      use module_cgsolver,only : aepscg,repscg,itercg
      use module_material,only : lclsd,ical_sld
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,kdprdc,ical_buff,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang,relaxp
      use module_model,only    : ical_week
      use module_model,only    : ical_vect,ical_MHD
!
      implicit none
!
! 1. Implicit operation for Poisson equation
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nflud,iter,IMODE
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: ipfix  (    MXMAT) 
      integer,intent(in)    :: kdbf   (  MXCVFAC)
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: diag   (  MXALLCV)
      real*8 ,intent(inout) :: rhs    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      integer,intent(out)   :: ierror,iterp
      real*8 ,intent(out)   :: repsp,aepsp
      real*8 ,intent(inout) :: rho    (  MXALLCV)
!
      integer,intent(in)    :: LCYCOLD(MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG (MXSSFBC_SLD)
      real*8 ,intent(in)    :: SFCENT (3, MXCVFAC)
!
! --- [local entities]
!
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8  :: dx,dy,dz,alf,aepsq,repsq,dl,ER,pref,anb,RMAXS
      real*8  :: dxo,dyo,dzo,dlvect,dl1,dl2,dl3
      real*8  :: dum1,dum2,dum3,dum4,th(2),costh1,sinth1,costh2,sinth2
      real*8  :: sc1,sc2,sc3
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf
      real*8  :: wi1,wi2,dln,dlo,alfn,alfo,avsfn,avsfo,avsfn1,avsfo1
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,ipref,INLB,INLA,INLO
      integer :: ICOM,IMD,ICH,IFLD,ICTP,II,IFLAG
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2),ISLD,ISLD2
      integer :: ICVS,ICVE,ICVL,IDCS,IDCE
      integer :: ICFS,ICFE,ICFL,IBFL,IBFS,IBFE
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP,ICFO,
     &           ICVAO,ICVBO,ICVB1,IPR
      integer :: nl,nu,IL,IU,iterq,ierr,ierr1
      logical :: lg_IE
      integer :: IALLOCA=0
      integer,save :: IICAL=0
!      integer,allocatable :: ipx(:)
      character(20)       :: cn
!
!--------------------------------------------------------------
!
      ierr1=0
      ierror=0
      IALLOCA=0
!
      allocate (IAU(MXCV,MAXIE),stat=ierr1)
      if(ierr1.ne.0) then
        write(ifle,*) 'allocating IAU(:,:) error in solve_poisson'
        call FFRABORT(1,'solve_poisson')
      endif
!
!--------------------------
! --- re-allocate array  --
!--------------------------
      IALLOCA=0
 1000 continue
      lg_IE=.false.
      if(IALLOCA==1) then
        allocate (IAU(MXCV,MAXIE),
     &            ipx(MAXIE),
     &            IAL(MXALLCV,MAXIE),
     &            aalw(MXALLCV,0:MAXIE),
     &            aax(MAXIE),
     &            stat=ierr1)
        if(ierr1.ne.0) then
          write(ifle,*) 'ERR: re-allocating error in solve_poisson' 
          call FFRABORT(1,'solve_poisson')
      endif
      endif
!-----------------------------------
!-< 1. Make up coefficient matrix >-
!-----------------------------------
!--< 1.1 all over the domain >--
!-----------------------------------
      if(ical_sld/=0)  then
        call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!-----------------------------------
      bb(:)=0.0d0
      INL(:)=0
      INU(:)=0
      IAL(:,:)=0
      IAU(:,:)=0
      aalw(:,0:MAXIE)=0.d0
      do 100 IIMAT=1,NMAT     !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) goto 100
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      aalw(ICVL,0)=diag(ICVL)
      enddo
  100 continue

!----------------------------------------------
      if(IMODE==6) then
        DO IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
!          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
!          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
!          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!          dl=abs(dx*SFAREA(1,ICFL)
!     &          +dy*SFAREA(2,ICFL)
!     &          +dz*SFAREA(3,ICFL))+SML
!
          unit(1,1)=CVCENT(1,ICVB)-SFCENT(1,ICFL)
          unit(2,1)=CVCENT(2,ICVB)-SFCENT(2,ICFL)
          unit(3,1)=CVCENT(3,ICVB)-SFCENT(3,ICFL)
          unit(1,2)=SFCENT(1,ICFL)-CVCENT(1,ICVA)
          unit(2,2)=SFCENT(2,ICFL)-CVCENT(2,ICVA)
          unit(3,2)=SFCENT(3,ICFL)-CVCENT(3,ICVA)
          th(1)=dsqrt(unit(1,1)**2+unit(2,1)**2+unit(3,1)**2)
          th(2)=dsqrt(unit(1,2)**2+unit(2,2)**2+unit(3,2)**2)
          dl=th(1)+th(2)
          alf=deltt*SFAREA(4,ICFL)/dl
          if(kdbf(ICFL).ne.1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            nl=min(ICVA,ICVB)
            nu=max(ICVA,ICVB)
!
            INL(nl)=INL(nl)+1
            IAL(nl,INL(nl))=nu
!
            INU(nu)=INU(nu)+1
            IAU(nu,INU(nu))=nl
!
            aalw(nl,INL(nl))=-alf
          elseif(kdbf(ICFL).eq.1) then
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
        enddo
        enddo
      else
        DO 110 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do 126 ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)+SML
          dl=abs(dx*SFAREA(1,ICFL)
     &          +dy*SFAREA(2,ICFL)
     &          +dz*SFAREA(3,ICFL))+SML
          alf=deltt*SFAREA(4,ICFL)/dl
!------------------
          if(kdbf(ICFL).ne.1) then
! --- regular (& periodic) face
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            nl=min(ICVA,ICVB)
            nu=max(ICVA,ICVB)
!
            INL(nl)=INL(nl)+1
            IAL(nl,INL(nl))=nu
!
            INU(nu)=INU(nu)+1
            IAU(nu,INU(nu))=nl
!
            aalw(nl,INL(nl))=-alf
          elseif(kdbf(ICFL).eq.1) then
! --- outlet BC:kdbf(ICFL)=1
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
 126    continue
 110    continue

      endif
!----------------------------------
! --- 
!----------------------------------
      do nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if(kd==kdprdc.and.idis(nb)==1) then  !zhang???
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
          ICFP=LCYCSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          rhs(ICVA)=0.d0
          enddo
        enddo
      endif
      enddo
!
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
        IIMATS=0
        IMATS=0
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdsld.and.idis(nb)==1) then
          do 400 ISLD=1,2
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
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          th(ISLD)=rot_ang(IMATS(ISLD))
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          th(ISLD2)=rot_ang(IMATS(ISLD2))
!
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2+unit(2,ISLD2)**2+unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          sf1=SFAREA(1,ICFL)
          sf2=SFAREA(2,ICFL)
          sf3=SFAREA(3,ICFL)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
!          dl1=deltt/dum1
          dl1=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl1=deltt/dl1
!
          IVA=LVEDGE(1,ICFP)
          IVB=LVEDGE(2,ICFP)
          dx=CVCENT(1,IVB)-CVCENT(1,IVA)
          dy=CVCENT(2,IVB)-CVCENT(2,IVA)
          dz=CVCENT(3,IVB)-CVCENT(3,IVA)
          sf1=SFAREA(1,ICFP)
          sf2=SFAREA(2,ICFP)
          sf3=SFAREA(3,ICFP)
!
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl2=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl2=SFAREA(4,ICFP)/dl2*deltt
!
          ICVAO=LVEDGE(1,ICFO)
          ICVBO=LVEDGE(2,ICFO)
          dx=CVCENT(1,ICVBO)-CVCENT(1,ICVAO)
          dy=CVCENT(2,ICVBO)-CVCENT(2,ICVAO)
          dz=CVCENT(3,ICVBO)-CVCENT(3,ICVAO)
          sf1=SFAREA(1,ICFO)
          sf2=SFAREA(2,ICFO)
          sf3=SFAREA(3,ICFO)
!
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl3=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl3=SFAREA(4,ICFO)/dl3*deltt
!
          avsf =dl1*SFAREA(4,ICFL)
          avsfn=dl1*SFAREA(4,ICFL)*wi1
          avsfo=dl1*SFAREA(4,ICFL)*wi2
          avsfn1=dl2  !*wi1
          avsfo1=dl3  !*wi2
!
          nl=min(ICVA,IVA)
          nu=max(ICVA,IVA)
          INLB=INL(nl)+1
          if(INLB<=MAXIE) then
            INL(nl)=INLB
            IAL(nl,INLB)=nu
            if(nl==ICVA) then
              aalw(nl,INLB)=-avsf
              aalw(nu,0)=aalw(nu,0)+avsf
              aalw(nl,0)=aalw(nl,0)+avsf
            else
              aalw(nl,INLB)=-avsfn
              aalw(nu,0)=aalw(nu,0)+avsfn
              aalw(nl,0)=aalw(nl,0)+avsfn
            endif
          else
             INL(nl)=INLB
             lg_IE=.true.
          endif
!
          nl=min(ICVA,ICVAO)
          nu=max(ICVA,ICVAO)
          INLO=INL(nl)+1
!
          if(INLO<=MAXIE) then
            INL(nl)=INLO
            IAL(nl,INLO)=nu
            if(nl==ICVA) then
              aalw(nl,INLB)=-avsf
              aalw(nu,0)=aalw(nu,0)+avsf
              aalw(nl,0)=aalw(nl,0)+avsf
            else
              aalw(nl,INLO)=-avsfo
              aalw(nu,0)=aalw(nu,0)+avsfo
              aalw(nl,0)=aalw(nl,0)+avsfo
            endif
          else
            INL(nl)=INLO
            lg_IE=.true.
          endif
          enddo
 400      enddo
        elseif(kd==kdsld.and.idis(nb)==2) then
          do 500 ISLD=1,2
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
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          th(ISLD)=rot_ang(IMATS(ISLD))
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
!          th(ISLD2)=rot_ang(IMATS(ISLD2))
!
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2+unit(2,ISLD2)**2+unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
!          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
!          do i=1,3
!          do j=1,3
!          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
!          enddo
!          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          th(ISLD2)=OPPANG(IBFL)
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          sf1=SFAREA(1,ICFL)
          sf2=SFAREA(2,ICFL)
          sf3=SFAREA(3,ICFL)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl1=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl1=deltt/dl1
!
          IVA=LVEDGE(1,ICFP)
          IVB=LVEDGE(2,ICFP)
          dx=CVCENT(1,IVB)-CVCENT(1,IVA)
          dy=CVCENT(2,IVB)-CVCENT(2,IVA)
          dz=CVCENT(3,IVB)-CVCENT(3,IVA)
          sf1=SFAREA(1,ICFP)
          sf2=SFAREA(2,ICFP)
          sf3=SFAREA(3,ICFP)
!
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl2=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl2=SFAREA(4,ICFP)/dl2*deltt
!
          ICVAO=LVEDGE(1,ICFO)
          ICVBO=LVEDGE(2,ICFO)
          dx=CVCENT(1,ICVBO)-CVCENT(1,ICVAO)
          dy=CVCENT(2,ICVBO)-CVCENT(2,ICVAO)
          dz=CVCENT(3,ICVBO)-CVCENT(3,ICVAO)
          sf1=SFAREA(1,ICFO)
          sf2=SFAREA(2,ICFO)
          sf3=SFAREA(3,ICFO)
!
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl3=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl3=SFAREA(4,ICFO)/dl3*deltt
!
          avsf =dl1*SFAREA(4,ICFL)
          avsfn=dl1*SFAREA(4,ICFL)*wi1
          avsfo=dl1*SFAREA(4,ICFL)*wi2
          avsfn1=dl2  !*wi1
          avsfo1=dl3  !*wi2
!
          nl=min(ICVA,IVA)
          nu=max(ICVA,IVA)
          INLB=INL(nl)+1
          if(INLB<=MAXIE) then
            INL(nl)=INLB
            IAL(nl,INLB)=nu
            if(nl==ICVA) then
              aalw(nl,INLB)=-avsf
              aalw(nu,0)=aalw(nu,0)+avsf
              aalw(nl,0)=aalw(nl,0)+avsf
            else
              aalw(nl,INLB)=-avsfn
              aalw(nu,0)=aalw(nu,0)+avsfn
              aalw(nl,0)=aalw(nl,0)+avsfn
            endif
          else
             INL(nl)=INLB
             lg_IE=.true.
          endif
!
          nl=min(ICVA,ICVAO)
          nu=max(ICVA,ICVAO)
          INLO=INL(nl)+1
!
          if(INLO<=MAXIE) then
            INL(nl)=INLO
            IAL(nl,INLO)=nu
            if(nl==ICVA) then
              aalw(nl,INLB)=-avsf
              aalw(nu,0)=aalw(nu,0)+avsf
              aalw(nl,0)=aalw(nl,0)+avsf
            else
              aalw(nl,INLO)=-avsfo
              aalw(nu,0)=aalw(nu,0)+avsfo
              aalw(nl,0)=aalw(nl,0)+avsfo
            endif
          else
            INL(nl)=INLO
            lg_IE=.true.
          endif
          enddo
 500      enddo
!
        elseif(kd==kdovst) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVP=LCYCSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
!
          avsf =deltt*SFAREA(4,ICFL)/dum1
!
          ICVA=msk(ICVA)
          ICVB=msk(ICVB)
          nl=min(ICVA,ICVB)
          nu=max(ICVA,ICVB)
!          INLA=INL(ICVA)+1
!          INLB=INL(ICVB)+1
          if(INLA<=MAXIE) then
!            INL(ICVA)=INLA
!            IAL(ICVA,INLA)=ICVB
!            INL(ICVB)=INLB
!            IAL(ICVB,INLB)=ICVA
!            aalw(ICVA,INLA)=-avsf
!            aalw(ICVB,INLB)=-avsf
            aalw(ICVA,0)=aalw(ICVA,0)+avsf
            aalw(ICVB,0)=aalw(ICVB,0)+avsf

!            INL(nl)=INLB
!            IAL(nl,INLB)=nu
!            if(nl==ICVA) then
!              aalw(nl,INLB)=-avsf
!              aalw(nu,0)=aalw(nu,0)+avsf
!              aalw(nl,0)=aalw(nl,0)+avsf
!            else
!              aalw(nl,INLB)=-avsf
!              aalw(nu,0)=aalw(nu,0)+avsf
!              aalw(nl,0)=aalw(nl,0)+avsf
!            endif
          endif
!
          enddo
          
        ENDIF
        enddo
      endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!-------------------
! --- check MAXIE
!-------------------
      if(NPE>1) then
        call hpclor(lg_IE)
      endif
      if(lg_IE) then
          INLO=0
          do ICV=1,NCV
            INLB=INL(ICV)
            if(INLB>INLO) then
              INLO=INLB
            endif
          enddo
          if(NPE>1) then
            call hpcimax(INLO)
          endif
          if(my_rank==root) then
            write(ifle,*)' ### MSG:MAX. INL(NL)= ',INLO,' MAXIE= ',MAXIE
            write(ifle,*)' ### MSG: ',
     &         'Increasing-3 [MAXIE] and [MAXIE] to :',INLO
          endif
          IEMAX=INLO
          MAXIE=INLO
          deallocate(IAU,ipx,IAL,aalw,aax)
          if(ierr1.ne.0) then
            write(ifle,*) 'ERR: allocating error in solve_poisson'
            call FFRABORT(1,'MSG: Call your FFR supportor')
          endif
          IALLOCA=1
          goto 1000
      endif
!---------------------------
! --- Coef. of Ap>=SIGMA(An)
!---------------------------
      do 200 ICV=1,NCV
        RMAXS=0.d0
        do 215 IV=1,INL(ICV)
          RMAXS=RMAXS-aalw(ICV,IV)
 215    enddo
        ipx(:)=0
        do IV=1,INU(ICV)-1
          IE=IAU(ICV,IV)
          do i=IV+1,INU(ICV)
            if(IAU(ICV,i)==IE) then
              ipx(i)=1
            endif
          enddo
        enddo
        do 220 IV=1,INU(ICV)
          if(ipx(IV)==1) cycle
          IE=IAU(ICV,IV)
          do 225 j=1,INL(IE)
            if(IAL(IE,j)==ICV) then
              RMAXS=RMAXS-aalw(IE,j)
            endif
 225      enddo
 220    enddo
        if(aalw(ICV,0)<RMAXS) then
          aalw(ICV,0)=RMAXS
        endif
 200  enddo
!      endif
!
!----------------------------
!--< 1.2 Clear solid part >--
!----------------------------
!
      IF(IMODE==6) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT) 
        bb(ICVS:ICVE)=rhs(ICVS:ICVE)
!        DO ICVL=ICVS,ICVE
!        if(ABS(aalw(ICVL,0))<SML) then
!          bb(ICVL)=0.d0
!          aalw(ICVL,0)=1.d0
!        endif
!        enddo
        enddo
      else
        do 120 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT.gt.0) then
          bb(ICVS:ICVE)=rhs(ICVS:ICVE)
        else
          bb(ICVS:ICVE)=0.d0
          INL(ICVS:ICVE)=0
          aalw(ICVS:ICVE,0)=1.d0
        endif
 120    enddo
      endif
!
!--< 1.3 Pressure fix point >-
!

      do 130 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        if(ipfix(IIMAT).gt.0) then
! --- no outlet BC:
          n=ipfix(IIMAT)
          bb(n)=0.d0
          INL(n)=0
          aalw(n,0)=1.d0
          do 131 i=1,INU(n)
            m=IAU(n,i)
            do 132 j=1,INL(m)
            if(IAL(m,j).eq.n) aalw(m,j)=0.d0
  132       enddo
  131     enddo
        endif
  130 enddo
!
      DEALLOCATE(IAU)
!
! --- 
!
      IF(NPE.GT.1) THEN 
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,AALW(:,0))
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
!-< 2. Solve linear system >-
!
      aepsq=aepscg
      repsq=repscg
      iterq=itercg
!
! --- ICCG solver (2D array) 
!
      if(IMODE==6) then
        if(NPE.gt.1) then
          call utl_iccg_hpc(iter,
     &      MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &      MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &      INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        else
          call utl_iccg
     &      (iter,MXCV,NCV,NMAT,MAXIE,MXMAT,MXALLCV,
     &      MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &      INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        endif
      else
!!!!!      if(ical_sld==0.and.ical_buff==0.and.ical_MHD>0) then 
      if((ical_sld==0.and.ical_buff==0).or.ical_MHD>0) then !important
        if(NPE.gt.1) then
          call utl_iccg_hpc(iter,
     &      MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &      MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &      INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        else
          call utl_iccg
     &      (iter,MXCV,NCV,NMAT,MAXIE,MXMAT,MXALLCV,
     &      MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &      INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        endif
      else   !strong coupling
        if(NPE.gt.1) then
          call utl_iccg_hpc_MAT(iter,ical_sld,
     &    MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &    MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &    INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        else
! --- < GOOD for serial and HPC (ical_sld==1 or 2)  > -------------------
          call utl_iccg_hpc_MAT(iter,ical_sld,
     &    MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &    MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &    INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        endif
      endif
      endif
!
      aepsp=aepsq
      repsp=repsq
      iterp=iterq
!
!=======================================================================
!
      do 300 IIMAT=1,NMAT    !ICV=1,NCV
      if(mat_cal(IIMAT)) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        rhs(ICVL)=bb(ICVL)
        enddo
      endif
  300 continue

!
      if( ierr1.ne.0 ) then
        write(ifll,6000) my_rank,iterq,repsq,aepsq
        write(ifle,*) ' ### error : iccg not converged'
        call FFRABORT(1,'iccg not converged')
      endif
!=======================================================================
      return
!
 6000 format(/,2X,'|',108('-'),'|',
     &  /,2X,' ### MSG: my_rank: ',i10,' iccg_iter= ',i10,
     &  /,2X,' ### MSG: maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     &  /,2X,'|',108('-'),'|')
!
 9999 continue
      write(ifle,*) '(solve_poisson)'
      ierror=1
      end subroutine solve_poisson
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_poisson_e2p
     & (iter,nflud,deltt,time,ISING,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  ipfix,kdbf,diag,rhs,rho,rho2,aks,
     &  LCYCOLD,wifsld,
     &  iterp,repsp,aepsp,IMODE,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- CG SOLVER
!
      
      use module_metrix,only : msk
      use module_metrix,only : IAL
      use module_metrix,only : IAU=>IW2K2
      use module_metrix,only : INL
      use module_metrix,only : INU
      use module_metrix,only : aalw
      use module_metrix,only : bb,ipx
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_io,only       : ifll,ifle
      use module_cgsolver,only : aepscg,repscg,itercg,
     &                           aepscg_p,repscg_p,itercg_p
      use module_material,only : lclsd,ical_sld,nofld,nosld
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,kdprdc,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang
      use module_model,only    : ical_week
      use module_model,only    : ical_vect,ical_MHD
      use module_scalar,  only : iaph
      use module_FUEL  ,only   : vap_no,No_Mem,Tau_diff,OPEN_V
!
      implicit none 
!
! 1. Implicit operation for Poisson equation
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nflud,iter,IMODE,ISING
      real*8 ,intent(in)    :: deltt,time
      integer,intent(inout)    :: ipfix  (    MXMAT) 
      integer,intent(in)    :: kdbf   (  MXCVFAC)
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: diag   (  MXALLCV)
      real*8 ,intent(inout) :: rhs    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      integer,intent(out)   :: ierror,iterp
      real*8 ,intent(inout)   :: repsp,aepsp
      real*8 ,intent(in)    :: rho    (  MXALLCV)
      real*8 ,intent(in)    :: rho2   (  MXALLCVC)
      real*8 ,intent(in)    :: aks    (  MXALLCVR,mxrans)
!
      integer,intent(in)    :: LCYCOLD(MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (MXSSFBC_SLD)
!
! --- [local entities]
!
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8  :: dx,dy,dz,alf,aepsq,repsq,dl,ER,pref,anb,RMAXS
      real*8  :: dxo,dyo,dzo,dlvect,dl1,dl2,dl3,dumA,dumB
      real*8  :: dum1,dum2,dum3,dum4,th(2),costh1,sinth1,costh2,sinth2
      real*8  :: sc1,sc2,sc3
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf
      real*8  :: wi1,wi2,dln,dlo,alfn,alfo,avsfn,avsfo,avsfn1,avsfo1
      real*8  :: aph1A,aph2A,aph1B,aph2B
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,ipref,INLB,INLO
      integer :: ICOM,IMD,ICH,IFLD,ICTP,II,IFLAG
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2),ISLD,ISLD2
      integer :: ICVS,ICVE,ICVL,IDCS,IDCE
      integer :: ICFS,ICFE,ICFL,IBFL,IBFS,IBFE
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP,ICFO,
     &           ICVAO,ICVBO,ICVB1,IPR
      integer :: nl,nu,IL,IU,iterq,ierr,ierr1
      logical :: lg_IE
      integer :: IALLOCA=0,IMAT_U
      integer,save :: IICAL=0
!      integer,allocatable :: ipx(:)
      character(20)       :: cn
!
!--------------------------------------------------------------
!
      ierr1=0
      ierror=0
      IALLOCA=0
!
      allocate (IAU(MXCV,MAXIE),stat=ierr1)
      if(ierr1.ne.0) then
        write(ifle,*) 'allocating IAU(:,:) error in solve_poisson'
        call FFRABORT(1,'solve_poisson')
      endif
!
!--------------------------
! --- re-allocate array  --
!--------------------------
!
      IALLOCA=0
 1000 continue
      lg_IE=.false.
      if(IALLOCA==1) then
        allocate (IAU(MXCV,MAXIE),
     &            ipx(MAXIE),
     &            IAL(MXALLCV,MAXIE),
     &            aalw(MXALLCV,0:MAXIE),
     &            
     &            stat=ierr1)
        if(ierr1.ne.0) then
          write(ifle,*) 'ERR: re-allocating error in solve_poisson'
          call FFRABORT(1,'solve_poisson')
      endif
      endif
!-----------------------------------
!-< 1. Make up coefficient matrix >-
!-----------------------------------
!--< 1.1 all over the domain >--
!-----------------------------------
      if(ical_sld/=0) then 
        call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!-----------------------------------
      bb(:)=0.0d0
      INL(:)=0
      INU(:)=0
      IAL(:,:)=0
      IAU(:,:)=0
      aalw(:,0:MAXIE)=0.d0
      do 100 IIMAT=1,NMAT     !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) goto 100
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      aalw(ICVL,0)=diag(ICVL)
      enddo
  100 continue
!----------------------------------------------
      if(IMODE==5.OR.IMODE==6) then
        DO 110 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do 126 ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)+SML
          dl=abs(dx*SFAREA(1,ICFL)
     &          +dy*SFAREA(2,ICFL)
     &          +dz*SFAREA(3,ICFL))+SML
          aph1A=aks(ICVA,iaph(1))
          aph1B=aks(ICVB,iaph(1))
          aph2A=aks(ICVA,iaph(2))
          aph2B=aks(ICVB,iaph(2))
          dumA=aph1A/rho(ICVA)+aph2A/rho2(ICVA)
          dumB=aph1B/rho(ICVB)+aph2B/rho2(ICVB)
          dum4=wi1*dumA+wi2*dumB
          alf=dum4*deltt*SFAREA(4,ICFL)/dl
!------------------
          if(kdbf(ICFL).ne.1) then
! --- regular (& periodic) face
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            nl=min(ICVA,ICVB)
            nu=max(ICVA,ICVB)
!
            INL(nl)=INL(nl)+1
            IAL(nl,INL(nl))=nu
!
            INU(nu)=INU(nu)+1
            IAU(nu,INU(nu))=nl
!
            aalw(nl,INL(nl))=-alf
          elseif(kdbf(ICFL).eq.1) then
! --- outlet BC:kdbf(ICFL)=1
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
 126    enddo
 110    enddo
      ELSEIF(IMODE==7) THEN
        DO IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)+SML
          dl=abs(dx*SFAREA(1,ICFL)
     &          +dy*SFAREA(2,ICFL)
     &          +dz*SFAREA(3,ICFL))+SML
          
          dumA=rho(ICVA)
          dumB=rho(ICVB)
          dum4=wi1*dumA+wi2*dumB
!          dum4=1.d0/(SML+wi1/dumA+wi2/dumB)
          alf=dum4*SFAREA(4,ICFL)/dl!*deltt
!------------------
          if(kdbf(ICFL).ne.1) then
! --- regular (& periodic) face
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            nl=min(ICVA,ICVB)
            nu=max(ICVA,ICVB)
!
            INL(nl)=INL(nl)+1
            INU(nu)=INU(nu)+1
            if(INL(nl)>MAXIE.or.INU(nu)>MAXIE) then
              lg_IE=.true.
              cycle
            endif
!
            IAL(nl,INL(nl))=nu
            IAU(nu,INU(nu))=nl
!
            aalw(nl,INL(nl))=-alf
          elseif(kdbf(ICFL).eq.1) then
! --- outlet BC:kdbf(ICFL)=1
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
        enddo
        enddo
      ENDIF
!
!----------------------------------
! --- 
!----------------------------------
!
      IF(IMODE/=7) THEN
        do 120 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT.gt.0) then
          bb(ICVS:ICVE)=rhs(ICVS:ICVE)
        else
          bb(ICVS:ICVE)=0.d0
          INL(ICVS:ICVE)=0
          aalw(ICVS:ICVE,0)=1.d0
        endif
 120    enddo
      elseif(IMODE==7)then
         
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        bb(ICVS:ICVE)=rhs(ICVS:ICVE)
        enddo
        if(ISING==2) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT.gt.0) then
            IMAT_U=nofld(IMAT)
          else
            IMAT_U=nosld(-IMAT)
          endif
          if(IMAT_U==No_Mem) cycle 
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          bb(ICVS:ICVE)=0.d0
          aalw(ICVS:ICVE,0)=1.d0
          INL(ICVS:ICVE)=0
          enddo
        elseif(ISING==1) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT.gt.0) then
            IMAT_U=nofld(IMAT)
          else
            IMAT_U=nosld(-IMAT)
          endif
          if(IMAT_U==No_Mem) then 
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            bb(ICVS:ICVE)=0.d0 
            aalw(ICVS:ICVE,0)=1.d0 
            INL(ICVS:ICVE)=0 
          endif
          enddo 
        endif
      ENDIF
!
! --- 
!
      if(NPE>1) then
        call hpclor(lg_IE)
      endif
!
      if(lg_IE) then
          INLO=0
          do ICV=1,NCV
            INLB=INL(ICV)
            if(INLB>INLO) then
              INLO=INLB
            endif
          enddo
          if(NPE>1) then
            call hpcimax(INLO)
          endif
          if(my_rank==root) then
            write(ifle,*)' ### MSG:MAX. INL(NL)= ',INLO,' MAXIE= ',MAXIE
            write(ifle,*)' ### MSG: ',
     &         'Increasing-4 [MAXIE] and [MAXIE] to :',INLO
          endif
          IEMAX=INLO
          MAXIE=INLO
          deallocate(IAU,ipx,IAL,aalw)
          if(ierr1.ne.0) then
            write(ifle,*) 'ERR: allocating error in solve_poisson'
            call FFRABORT(1,'MSG: Call your FFR supportor')
          endif
          IALLOCA=1
          goto 1000
      endif
!--< 1.3 Pressure fix point >-
!
      do 130 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        if(ipfix(IIMAT).gt.0) then
          n=ipfix(IIMAT)
          bb(n)=0.d0
          INL(n)=0
          aalw(n,0)=1.d0
          do 131 i=1,INU(n)
            m=IAU(n,i)
            do 132 j=1,INL(m)
            if(IAL(m,j).eq.n) aalw(m,j)=0.d0
  132       enddo
  131     enddo
        endif
  130 enddo
!
      DEALLOCATE(IAU)
!
! --- 
!
      IF(NPE.GT.1) THEN 
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,AALW(:,0))
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
!-< 2. Solve linear system >-
!
      if(IMODE==7) then
        aepsq=aepscg_p
        repsq=repscg_p
        iterq=itercg_p
      else
        aepsq=aepscg
        repsq=repscg
        iterq=itercg
      endif
!
! --- ICCG solver (2D array) 
!
      if(NMAT==1) then
        if(NPE.gt.1) then
          call utl_iccg_hpc(iter,
     &      MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &      MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &      INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        else
          call utl_iccg
     &      (iter,MXCV,NCV,NMAT,MAXIE,MXMAT,MXALLCV,
     &      MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &      INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        endif
      else
        if(NPE.gt.1) then
          call utl_iccg_hpc_MAT(iter,ical_sld,
     &    MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &    MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &    INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)
        else
! --- < GOOD for serial and HPC (ical_sld==1 or 2)  > -------------------
          call utl_iccg_hpc_MAT(iter,ical_sld,
     &    MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &    MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &    INL,IAL,aalw,bb,aepsq,repsq,iterq,IMODE,ierr1)

        endif
      endif
!
      aepsp=aepsq
      repsp=repsq
      iterp=iterq
!
!=======================================================================
!
      do 300 IIMAT=1,NMAT    !ICV=1,NCV
      if(mat_cal(IIMAT)) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        rhs(ICVL)=bb(ICVL)
        enddo
      endif
  300 continue

!
      if( ierr1.ne.0 ) then
        write(ifll,6000) my_rank,iterq,repsq,aepsq
        write(ifle,*) ' ### error : iccg not converged'
        call FFRABORT(1,'iccg not converged')
      endif
!=======================================================================
      return
!
 6000 format(/,2X,'|',108('-'),'|',
     &  /,2X,' ### MSG: my_rank: ',i10,' iccg_iter= ',i10,
     &  /,2X,' ### MSG: maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     &  /,2X,'|',108('-'),'|')
!
 9999 continue
      write(ifle,*) '(solve_poisson)'
      ierror=1
      end subroutine solve_poisson_e2p
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_poisson_unsymm
     & (iter,nflud,deltt,time,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,SFCENT,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  ipfix,kdbf,diag,rhs,rho,rva,ccc,tmp,vel,
     &  LCYCOLD,wifsld,
     &  iterp,repsp,aepsp,IMODE,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine is for calculating NS equations
      use module_dimension
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IQ =>IW2K1
      use module_metrix,only   : bb
      use module_metrix,only   : aalw
      use module_material,only : ical_sld
!
! --- [module arguments]
!
      use module_hpcutil
      use module_io,only       : ifll,ifle,cntlnam
      use module_cgsolver,only : aepscg,repscg,itercg
      use module_flags,only    :
     &                           Adams_Moulton,intgvv,Crank_Nicolson
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,kdprdc,ical_buff,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang
      use module_species,only  : spcnam,gascns
      use module_scalar ,only  : sclname
      use module_metrix,only   : aax,ipx
      use module_model,only    : ical_vect
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nflud,iter,IMODE
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: ipfix  (    MXMAT) 
      integer,intent(in)    :: kdbf   (  MXCVFAC)
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: diag   (  MXALLCV)
      real*8 ,intent(inout) :: rhs    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      integer,intent(out)   :: ierror,iterp
      real*8 ,intent(out)   :: repsp,aepsp
      real*8 ,intent(inout) :: rho    (  MXALLCV)
      real*8 ,intent(in)    :: ccc    (  MXALLCV)
      real*8 ,intent(inout) :: rva    (  MXCVFAC)
      real*8 ,intent(inout) :: vel    (  MXALLCV ,3)
!
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: tmp    (  MXALLCV)

      real*8 ,intent(in)    :: SFCENT (3, MXCVFAC)
!
! --- [local entities]
!
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8  :: dumC,wi1,wi2,dlvect
      real*8  :: th(2),costh1,sinth1,costh2,sinth2
      integer :: it,icx,iterq,ierr,ierr1
      real*8  :: dx,dy,dz,alf,alfd,alfo,vlp,vlm,vdd,cof,alfn,dl
      real*8  :: wifa,wifb,vaa,vab,vbb,vba,dl1,dxo,dyo,dzo
      real*8  :: aepsq,repsq,calph,RMAXS
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf
      real*8  :: ru,rv,rw,dum1,dum2,dum3,dum4,dum5
      real*8  :: avsfo,dln,dl2,dl3,avsfn,alfdn,alfdo,alfe,alfen,
     &           alfeo,cofB,cofO,vlpn,vlmn,vddn,vlpo,vlmo,vddo,
     &           avsfn1,avsfo1,vlp1,vlm1
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,IPR
      integer :: IBFL,IBFS,IBFE,ICFO,ICVBO
      integer :: NB,KD,IDC,IV
      integer :: ICVA,ICVB,IVA,IVB,IBFP,ICFP,ICVP,IDCP,ICVAO
      integer :: ICFS,ICFE,ICFL
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICV,INLO,
     &           IIMATS(2),IMATS(2),ISLD,ISLD2,nl,nu,INLB,IA,IAA,IAO,IAB
      logical :: lg_IE
      integer :: IALLOCA=0
      character(20) :: cn
      integer :: IQMXCV
!
      ierr=0
      ierror=0
      ierr1=0
      IALLOCA=0
!
      IQMXCV=MXCV
      allocate(IQ(IQMXCV,2),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in solve_poisson_unsymm'
        call FFRABORT(1,'solve_poisson_unsymm')
      endif
!---------------------------------------------------------
!-< 1. Set maximul & minimul diffusivity for simplicity >-
!---------------------------------------------------------
      IALLOCA=0
 1000 continue
      lg_IE=.false.
      if(IALLOCA==1) then
        allocate (
     &            ipx(MAXIE),
     &            IAL(MXALLCV,MAXIE),
     &            aalw(MXALLCV,0:MAXIE),
     &            aax(MAXIE),
     &            stat=ierr1)
        if(ierr1.ne.0) then
          write(ifle,*) 
     &    'ERR: re-allocating error in unsymm_poisson'
          call FFRABORT(1,'solve_poisson')
        ENDIF
      endif
!-----------------------------------
!-< 2. Make up coefficient matrix >-
!-----------------------------------
!--< 2.1 all over the domain >--
!-----------------------------------
      if(ical_sld/=0) then
        call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!
      IQ(:,:)=0
      aalw(:,0)=0.d0
      do 200 IIMAT=1,NMAT   !ICV=1,NALLCV
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      aalw(ICVL,0)=diag(ICVL)
      enddo
  200 enddo
!
      if(IMODE==10) then
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dlvect=ABS(dx*SFAREA(1,ICFL)   !?????????
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          alf=SFAREA(4,ICFL)/(dlvect)*deltt
          ru=(wi1*vel(ICVA,1)+wi2*vel(ICVB,1))
          rv=(wi1*vel(ICVA,2)+wi2*vel(ICVB,2))
          rw=(wi1*vel(ICVA,3)+wi2*vel(ICVB,3))

          dum4=SFAREA(4,ICFL)*
     &    (SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw)

          dum3=wi1/ccc(ICVA)+wi2/ccc(ICVB)
          dum5=dum4*dum3    !zhang5 (9) 
!
          vlp=max(0.d0,dum5)
          vlm=min(0.d0,dum5)
          

          dum1=ccc(ICVA)*rho(ICVA)
          dum2=ccc(ICVB)*rho(ICVB)
!          vlp1=rva(ICFL)/dum1
!          vlm1=rva(ICFL)/dum2
          dum3=wi1/ccc(ICVA)+wi2/ccc(ICVB)
          dum5=wi1/rho(ICVA)+wi2/rho(ICVB)
          dum4=rva(ICFL)*dum3*dum5   !*deltt
          vlp=max(0.d0,dum4)
          vlm=min(0.d0,dum4)


          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE .or. IQ(ICVB,1)>MAXIE ) then 
              write(ifle,*) '(skip solve_poisson_unsymm)', 
     &              kdbf(ICFL),ICVA,ICVB
!              deallocate(aax,ipx,IAL,aalw)
!              return
              lg_IE=.true.
              cycle
            endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf-vlm !+vlp
            aalw(ICVB,0)=aalw(ICVB,0)+alf+vlp !-vlm
!
!            aalw(ICVA,0)=aalw(ICVA,0)+alf+vlp
!            aalw(ICVB,0)=aalw(ICVB,0)+alf-vlm
!
            aalw(ICVA,IQ(ICVA,1))=-alf-vlp !+vlm
            aalw(ICVB,IQ(ICVB,1))=-alf+vlm !-vlp

!            aalw(ICVA,IQ(ICVA,1))=-alf+vlm
!            aalw(ICVB,IQ(ICVB,1))=-alf-vlp
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf-vlm !+vlp
!            aalw(ICVA,0)=aalw(ICVA,0)+alf+vlp

          endif
        endif
        enddo
        enddo
!
      elseif(IMODE==11) then
        do IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dlvect=ABS(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          dum1=rho(ICVA)
          dum2=rho(ICVB)
          dum3=wi1*dum1+wi2*dum2
          alf=SFAREA(4,ICFL)/(dlvect)*deltt/dum3
          dum1=rho(ICVA)*ccc(ICVA)
          dum2=rho(ICVB)*ccc(ICVB)

          dum3=wi1*dum1+wi2*dum2
          vlp1=rva(ICFL)/dum3
          vlm1=rva(ICFL)/dum3
          vlp=max(0.d0,vlp1)+alf*rho(ICVA)
          vlm=min(0.d0,vlm1)-alf*rho(ICVB)
          alf=0.d0
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE .or. IQ(ICVB,1)>MAXIE ) then 
              write(ifle,*) '(skip solve_poisson_unsymm)', 
     &              kdbf(ICFL),ICVA,ICVB
              lg_IE=.true.
              cycle
           endif
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf+vlp
            aalw(ICVB,0)=aalw(ICVB,0)+alf-vlm
!
            aalw(ICVA,IQ(ICVA,1))=-alf-vlp  !+vlm
            aalw(ICVB,IQ(ICVB,1))=-alf+vlm
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf+vlp
          endif
        endif
        enddo
        enddo


      ELSE
        do IIMAT=1,NMAT        !ICF=1,NCVFAC 
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        if(kdbf(ICFL).lt.2) then
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dlvect=ABS(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          alf=SFAREA(4,ICFL)/(dlvect)*deltt
          if(kdbf(ICFL)/=1) then
            ICVA=msk(ICVA)
            ICVB=msk(ICVB)
            IQ(ICVA,1)=IQ(ICVA,1)+1
            IQ(ICVB,1)=IQ(ICVB,1)+1
            if(IQ(ICVA,1)>MAXIE .or. IQ(ICVB,1)>MAXIE ) then 
              write(ifle,*) '(skip solve_poisson_unsymm)', 
     &              kdbf(ICFL),ICVA,ICVB
!              deallocate(aax,ipx,IAL,aalw)
!              return
              lg_IE=.true.
              cycle
            endif
!
            IAL(ICVA,IQ(ICVA,1))=ICVB
            IAL(ICVB,IQ(ICVB,1))=ICVA
!
            aalw(ICVA,0)=aalw(ICVA,0)+alf
            aalw(ICVB,0)=aalw(ICVB,0)+alf
!
            aalw(ICVA,IQ(ICVA,1))=-alf
            aalw(ICVB,IQ(ICVB,1))=-alf
          else
            aalw(ICVA,0)=aalw(ICVA,0)+alf
          endif
        endif
        enddo
        enddo
      ENDIF
!-----------------------------------------
! --- 
!-----------------------------------------
      do nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if(kd==kdprdc.and.idis(nb)==1) then
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
          ICFP=LCYCSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          IVA=LVEDGE(1,ICFP)
          IVB=LVEDGE(2,ICFP)
!
          dx=CVCENT(1,IVB)-CVCENT(1,IVA)
          dy=CVCENT(2,IVB)-CVCENT(2,IVA)
          dz=CVCENT(3,IVB)-CVCENT(3,IVA)
          dl=   abs(dx*SFAREA(1,ICFP)
     &             +dy*SFAREA(2,ICFP)
     &             +dz*SFAREA(3,ICFP))
!          
          dum2=SFAREA(4,ICFP)/dl*deltt
          alfn=dum1/(dum1+dum2)
          aalw(ICVA,0)=aalw(ICVA,0)+alfn
          enddo
        enddo
      endif
      enddo
!------------------------------------------------------------
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then 
        IIMATS=0
        IMATS=0
        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdsld) then
          do 400 ISLD=1,2
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
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          th(ISLD)=rot_ang(IMATS(ISLD))
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          th(ISLD2)=rot_ang(IMATS(ISLD2))
!
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
!
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2+unit(2,ISLD2)**2+unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          sf1=SFAREA(1,ICFL)
          sf2=SFAREA(2,ICFL)
          sf3=SFAREA(3,ICFL)
!          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
!          dl1=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
!          dl1=deltt/dl1

          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl1=deltt/dum1
!
          IVA=LVEDGE(1,ICFP)
          IVB=LVEDGE(2,ICFP)
          dx=CVCENT(1,IVB)-CVCENT(1,IVA)
          dy=CVCENT(2,IVB)-CVCENT(2,IVA)
          dz=CVCENT(3,IVB)-CVCENT(3,IVA)
          sf1=SFAREA(1,ICFP)
          sf2=SFAREA(2,ICFP)
          sf3=SFAREA(3,ICFP)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl2=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl2=SFAREA(4,ICFP)/dl2*deltt
!
          ICVAO=LVEDGE(1,ICFO)
          ICVBO=LVEDGE(2,ICFO)
          dx=CVCENT(1,ICVBO)-CVCENT(1,ICVAO)
          dy=CVCENT(2,ICVBO)-CVCENT(2,ICVAO)
          dz=CVCENT(3,ICVBO)-CVCENT(3,ICVAO)
          sf1=SFAREA(1,ICFO)
          sf2=SFAREA(2,ICFO)
          sf3=SFAREA(3,ICFO)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dl3=(dx*sf1+dy*sf2+dz*sf3)**2/dum1
          dl3=SFAREA(4,ICFO)/dl3*deltt
!
          avsf =dl1*SFAREA(4,ICFL) 
          avsfn=dl1*SFAREA(4,ICFL)*wi1 
          avsfo=dl1*SFAREA(4,ICFL)*wi2 
          avsfn1=dl2  !*wi1
          avsfo1=dl3  !*wi2
!
          IQ(ICVA,1)=IQ(ICVA,1)+1      !1)
          IAA=IQ(ICVA,1)
!
          IQ(IVA,1)=IQ(IVA,1)+1        !2)
          IA=IQ(IVA,1)
!
          IQ(ICVA,1)=IQ(ICVA,1)+1
          IAB=IQ(ICVA,1)
!
          IQ(ICVAO,1)=IQ(ICVAO,1)+1    !3)
          IAO=IQ(ICVAO,1)
!
          aalw(ICVA,0)=aalw(ICVA,0)+avsf
          if(IAA<=MAXIE.and.IA<=MAXIE) then
            IAL(ICVA,IAA)=IVA
            aalw(ICVA,IAA)=-avsfn
            IAL(IVA,IA)=ICVA
            aalw(IVA,IA)=-avsfn
            aalw(IVA,0)=aalw(IVA,0)+avsfn !avsf
          else
            lg_IE=.true.
          endif
!
          if(IAB<=MAXIE.and.IAO<=MAXIE) then
            IAL(ICVA,IAB)=ICVAO
            IAL(ICVAO,IAO)=ICVA
            aalw(ICVA,IAB)=-avsfo
            aalw(ICVAO,IAO)=-avsfo
            aalw(ICVAO,0)=aalw(ICVAO,0)+avsfo !avsf
          else
            lg_IE=.true.
          endif
!
         enddo
 400     enddo


        elseif(kd==kdovst) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          IAA=0
          IA=0
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVP=LCYCSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
!
!          ICVA=msk(ICVA)
!          ICVB=msk(ICVB)
!
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &         +dz*SFAREA(3,ICFL))
!
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
!
          dl1=deltt/dlvect
!
          avsf=dl1*SFAREA(4,ICFL)
!
          ICVA=msk(ICVA)
          ICVB=msk(ICVB)
!
          if(IQ(ICVA,1)+1<=MAXIE.and.IQ(ICVB,1)+1<=MAXIE) then
            IQ(ICVA,1)=IQ(ICVA,1)+1 
            IAA=IQ(ICVA,1) 
            aalw(ICVA,0)=aalw(ICVA,0)+avsf 
            IAL(ICVA,IAA)=ICVB 
            aalw(ICVA,IAA)=-avsf 
! 
!            IQ(ICVB,1)=IQ(ICVB,1)+1  
!            IA=IQ(ICVB,1) 
!            aalw(ICVB,0)=aalw(ICVB,0)+avsf 
!            IAL(ICVB,IA)=ICVA 
!            aalw(ICVB,IA)=-avsf 
! 
          else 
            lg_IE=.true. 
            IQ(ICVA,1)=IQ(ICVA,1)+1 
            IQ(ICVB,1)=IQ(ICVB,1)+1  
          endif
!
         enddo

        ENDIF
        enddo
      endif
!
      if(NPE>1) then
        call hpclor(lg_IE)
      endif

      if(lg_IE) then
          INLO=0
          do ICV=1,NCV
          INLB=IQ(ICV,1)
          if(INLB>INLO) then
            INLO=INLB
          endif
          enddo
          if(NPE>1) then
             call hpcimax(INLO)
          endif
          if(my_rank==root) then
            write(ifle,*)
     &      ' ### MSG:MAX. IQ(ICV,1)= ',INLO,' MAXIE= ',MAXIE
            write(ifle,*)' ### MSG: ',
     &      'Increasing-5 [MAXIE] and [MAXIE] to :',INLO
          endif
          IEMAX=INLO
          MAXIE=INLO
          deallocate(IAL,aalw,ipx,aax)
          if(ierr1.ne.0) then
            write(ifle,*) 'ERR: allocating error in solve_poisson'
            call FFRABORT(1,'MSG: Call your FFR supportor')
          endif
          IALLOCA=1
          goto 1000
      endif
!---------------------------
! --- Coef. of Ap>=SIGMA(An)
!---------------------------
      do ICV=1,NCV
        RMAXS=0.d0
        do IV=1,IQ(ICV,1)
          RMAXS=RMAXS-aalw(ICV,IV)
        enddo
        if(aalw(ICV,0)<RMAXS) then
          aalw(ICV,0)=RMAXS
        endif
      enddo
!      endif
!-----------------------------------------
!--< 2.2 Clear dummy cell & solid part >--
!-----------------------------------------
      IQ(NCV+1:IQMXCV,1)=0
      do 220 IDC=NCV+1,IQMXCV
      aalw(IDC,0)=1.d0    !1.d0
  220 continue
!
      do 225 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          DO 230 ICVL=ICVS,ICVE
          IQ(ICVL,1)=0
          aalw(ICVL,0)=1.d0
 230      continue          
        endif
 225  continue
!
!-------------------------
!--< 2.4 sliding mesh  >--
!-------------------------
!
!-------------------------------------------------------
!--< 2.4 rearrange order to split lower & upper part >--
!-------------------------------------------------------
!      if(.false.) then
      do 240 IIMAT=1,NMAT   !ICV=1,NALLCV
        if(.not.mat_cal(IIMAT)) cycle 
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 255 ICVL=ICVS,ICVE
        icx=IQ(ICVL,1)
        k=0
        do 241 m=1,icx
        if(IAL(ICVL,m).lt.ICVL) then 
! --- lower part: IQ(ICVL,1)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  241   enddo
        IQ(ICVL,1)=k
        do 242 m=1,icx
        if(IAL(ICVL,m).gt.ICVL) then
! --- upper part: IQ(ICVL,2)
          k=k+1
          aax(k)=aalw(ICVL,m)
          ipx(k)=IAL(ICVL,m)
        endif
  242   continue
        IQ(ICVL,2)=k
        do 243 m=1,k
          aalw(ICVL,m)=aax(m)
          IAL(ICVL,m)=ipx(m)
  243   continue
        if( k/=icx ) then
          write(*,*) 'MSG: k=, icx=',k,icx,ICVL
          write(*,*) 'MSG: program error -1- (solve_poisson_unsymm)'
          CALL FFRABORT(1,'solve_poisson_unsymm')
        endif
 255    continue
 240  enddo
!      endif
!----------------------------
!-< 3. Solve linear system >-
!----------------------------
      bb(:)=0.d0
!
      do 310 IIMAT=1,NMAT    !ICV=1,NCV
      if(mat_cal(IIMAT)) then
        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_CVEXT(IIMAT)
        ICVE=MAT_INDEX(IIMAT)
        do ICVL=ICVS,ICVE
        bb(ICVL)=rhs(ICVL)
        enddo
      endif
  310 continue
!
! --- 
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,AALW(:,0))
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      aepsq=aepscg
      repsq=repscg
      iterq=itercg
!
      m=1
      cn='bicg_poisson'
      if(ical_sld==0.and.ical_buff==0) then   
        call utl_bcgstb(cn,m,
     &   MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aepsq,repsq,iterq,ierr1,IMODE)
      else
        call utl_bcgstb_MAT1(cn,
     &   MXALLCV,MXCV,MAXIE,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,IAL,aalw,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aepsq,repsq,iterq,ierr1,IMODE)
      endif
!
      IF(NPE.gt.1) THEN
        CALL hpcimax(iterq)
      ENDIF
!
      aepsp=aepsq
      repsp=repsq
      iterp=iterq
!
      if(my_rank==root.and.mod(iter,100)==0) then
        if(iterq==0) then
          write(ifle,'(2a)') ' ### WRN: BCGSTB solver NOT run at ',
     &   'unsymm poisson'
        endif
      endif
!
! --- 
!
      if(iterq/=0) then
        do 320 IIMAT=1,NMAT    !ICV=1,NCV
        if(mat_cal(IIMAT)) then
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          rhs(ICVL)=bb(ICVL)
          enddo
        endif
  320   continue
      endif
!
      deallocate(IQ)
!
      return
!
      write(ifll,6000) cn,iterq,repsq,aepsq,ierr1
 6000 format('ERR: bicg_poisson/',a,' : iterations=',i10,
     * ' / maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     * ' / return code=',i5)
!
 9999 continue
      write(ifle,*) '(solve_poisson_unsymm)'
      ierror=1
!
      end subroutine solve_poisson_unsymm
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_poisson_unsymm_vect
     & (iter,nflud,deltt,time,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  ipfix,kdbf,diag,rhs,rho,ccc,rva,
     &  LCYCOLD,wifsld,SFCENT,!vctr,
     &  iterp,repsp,aepsp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine is for calculating NS equations
      use module_dimension
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IAU=>IW2K2
      use module_metrix,only   : INL
      use module_metrix,only   : INU
!
      use module_metrix,only   : bb
      use module_metrix,only   : aalw
      use module_metrix,only   : aaup=>W2K1
!
      use module_metrix,only   : D  =>W1K6    !W1K1 W1K6
      use module_metrix,only   : X  =>d1work2
      use module_metrix,only   : Z  =>W1K2
      use module_metrix,only   : W  =>W2K2
      use module_metrix,only   : TEMP=>W1K3
!
      use module_metrix,only   : IW=>iwork4  !used in move.f 

      use module_metrix,only   : IT=>iwork5
      use module_metrix,only   : IJT=>iwork6
      use module_metrix,only   : LN=>iwork7

      use module_metrix,only   : aax,ipx
      use module_metrix,only   : tmpfac=>d2vect
      use module_metrix,only   : IW2K_VECT,DW2K_VECT
      use module_metrix,only   : vctrp
      use module_trace ,only   : NEIGHBOR,DNORM
!
      use module_vector  ,ONLY : input_vect    =>inputdata
      use module_metrix,only   : vctr
!
! --- [module arguments]
!
      use module_hpcutil
      use module_io,only       : ifli,ifll,ifle,cntlnam
      use module_cgsolver,only : aepscg,repscg,itercg
      use module_flags,only    :
     &                           Adams_Moulton,intgvv,Crank_Nicolson
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,kdprdc,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang
      use module_species,only  : spcnam
      use module_scalar ,only  : sclname
      use module_material,only : ical_sld
      use module_model,only    : ical_vect,nthrds,ical_prt,iLBF_P,
     &                           idrdp,comp
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f,
     &                           NCOLOR,NO,
     &                           MXUP,MXLW,N2
!
      implicit none
!------------------------
! --- [dummy arguments] 
!------------------------
      integer,intent(in)    :: nflud,iter
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: ipfix  (    MXMAT) 
      integer,intent(inout) :: kdbf   (  MXCVFAC)
      integer,intent(in)    :: LVEDGE (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF(  MXSSFBC)
      integer,intent(in)    :: LCYCSF (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV (  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA (4,MXCVFAC)
      real*8 ,intent(in)    :: wiface (  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT (3,MXALLCV)
      real*8 ,intent(in)    :: diag   (  MXALLCV)
      real*8 ,intent(inout) :: rhs    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOLM (  MXALLCV)
      integer,intent(out)   :: ierror,iterp
      real*8 ,intent(out)   :: repsp,aepsp
      real*8 ,intent(in)    :: rho    (  MXALLCV)
      real*8 ,intent(in)    :: ccc    (  MXALLCV)
!
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      REAL*8 ,INTENT(IN)    :: SFCENT (3, MXCVFAC) 
      real*8 ,intent(in)    :: rva    (   MXCVFAC)
!      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      real*8  :: dlvect,vlp,vlm,comdum=0.d0
      real*8  :: dum1,dum2,dum3,dum4
      real*8  :: dx,dy,dz,alf
      real*8  :: aepsq,repsq,calph
      integer :: i,j,k,l,m,n
      integer :: IMODE,IDC
      integer :: ICVA,ICVB,IMODE1
      integer :: ICFS,ICFE,ICFL
!      
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IE,myid,IV,iti=0,
     &           idum,nu,nl
      integer :: icx,iterq,ierr,ierr1,kl,ku,nb,kd,ISLD,ISLD2
      integer :: IBFS, IBFE,IIMATS(2),IMATS(2),ICVBO,ICVAO,INLB,IVA,
     &           IVB,ICV,INLO,IBFL,ICFO,ICFP
      character(20) :: cnx
      real*8,parameter :: SML=1.d-15
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf,
     &           wi1,wi2,th(2),RMAXS
      real*8  :: avsfn,avsfo,avsfn1,dl2,dl1,dl3,avsfo1,
     &           sinth1,sinth2,costh1,costh2
      logical :: lg_IE
      integer :: IA,IAB,IAO,IAA
      integer :: IALLOCA=0,IBFS1,IBFS2,IBFE1,IBFE2
      logical :: calsld
!

!      if(NPE==1) calsld=.false.
      if(idrdp.eq.comp) comdum=1.d0
      comdum=0.d0
      ierr=0
      ierror=0
      ierr1=0
      IALLOCA=0
!--------------------------
! --- re-allocate array  --
!--------------------------
      IALLOCA=0
 1000 continue
      lg_IE=.false.
      if(IALLOCA==1) then
!        call FFRABORT(1,'1239317')
        MXBND_V=MAXIE
        allocate (IAU(MXCV,MAXIE),stat=ierr1)
        allocate (ipx(MAXIE),stat=ierr1)
        allocate (IAL(MXALLCV,MAXIE),stat=ierr1)
        allocate (aalw(MXALLCV,0:MAXIE),stat=ierr1)
        allocate (aaup(MXCV,0:MAXIE),stat=ierr1)
        allocate (vctr(MXCV_V,0:MXBND_V),stat=ierr1)
        allocate (aax(MAXIE),stat=ierr1)
!
        allocate (IW2K_VECT(MXCV,MAXIE),stat=ierr1)
        allocate (DW2K_VECT(MXCV,MAXIE),stat=ierr1)
!
        if(ierr1.ne.0) then
          write(ifle,*) 'ERR: re-allocating error in solve_poisson'
          call FFRABORT(1,'solve_poisson')
        endif
        vctr=0
        call input_vect(2,
     &     ifli,ifll,ifle,ical_vect,nthrds,iLBF_P,ical_prt,
     &     NCV,NCVIN,NALLCV,NCVFAC,NMAT,NSSFBC,MXMAT,MAXIE,
     &     MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,
     &     vctr,vctrp,NEIGHBOR,DNORM,SFAREA,SFCENT,
     &     MXCV_V,MXBND_V,MXCV_V,MXBND_VP,
     &     msk,MXALLCV,LVEDGE,MXCVFAC,1,
     &     MEP,MXALLCV_P,MXCVFAC_P,
     &     ierror)
      endif
!-----------------------------------
!-< 2. Make up coefficient matrix >-
!-----------------------------------
!--< 2.1 all over the domain >--
!-----------------------------------
!
!-----------------------------------
! --- 
!-----------------------------------
      if(ical_sld/=0) then 
        call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!
      aalw(1:MXALLCV,0:MAXIE)=0.d0
      aaup(1:MXCV,0:MAXIE)=0.d0
!
!      if(NEW_TSTP) then
        IAL(1:MXALLCV,1:MAXIE)=0
        IAU(1:MXCV,1:MAXIE)=0
        INL(:)=0
        INU(:)=0
        IT(:)=0
        IJT(:)=0
        LN(:)=0
!      endif
!
      D(0)=0.d0
      D(1:NCV)=diag(1:NCV)
      IW(:)=0
!
!-------------------------------------------------------
!CIDR VECTOR
      do ICFL=ICFS_V,ICFE_V
        IF(kdbf(ICFL)/=1.and.kdbf(ICFL)/=0) kdbf(ICFL)=2
      enddo
!CIDR VECTOR
      do ICFL=ICFS_V,ICFE_V 
!        if(kdbf(ICFL)>=2) cycle
        ICVA=msk(LVEDGE(1,ICFL))
        ICVB=msk(LVEDGE(2,ICFL))
        tmpfac(ICFL,1)=SFAREA(4,ICFL)
     &        *deltt/(abs((
     &            CVCENT(1,ICVB)
     &           -CVCENT(1,ICVA))*SFAREA(1,ICFL)
     &          +(CVCENT(2,ICVB)
     &           -CVCENT(2,ICVA))*SFAREA(2,ICFL)
     &          +(CVCENT(3,ICVB)
     &           -CVCENT(3,ICVA))*SFAREA(3,ICFL))+SML)
      enddo 
!---------------------------------------------------------
      do IE=1,MAXIE   !kdbf=0,1,2 
CIDR VECTOR 
        DO ICVL=ICVS_V,ICVE_V
        if(vctr(ICVL,IE)==0) cycle        
        if(kdbf(ABS(vctr(ICVL,IE)))>=2) cycle 

        dum1=rho(ICVL)*ccc(ICVL)
        vlp=max(0.d0,rva(ABS(vctr(ICVL,IE)))/dum1)
        vlm=min(0.d0,rva(ABS(vctr(ICVL,IE)))/dum1)
!
        IF(vctr(ICVL,IE)<0) then       !ICVB 
          D(ICVL)=D(ICVL)
     &     +(2.d0-dble(kdbf(ABS(vctr(ICVL,IE)))))/
     &      (2.d0-dble(kdbf(ABS(vctr(ICVL,IE))))+SML) 
     &     *(tmpfac(ABS(vctr(ICVL,IE)),1)-vlm*comdum)
        elseIF(vctr(ICVL,IE)>0) then   !ICVA 
          D(ICVL)=D(ICVL)
     &     +(2.d0-dble(kdbf(ABS(vctr(ICVL,IE)))))/
     &      (2.d0-dble(kdbf(ABS(vctr(ICVL,IE))))+SML)
     &     *(tmpfac(ABS(vctr(ICVL,IE)),1)+vlp*comdum)
        ENDIF
        ENDDO
!
!        do IIMAT=1,NMAT
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!!CDIR NODEP
!        do ICVL=ICVS,ICVE
!        IF(ABS(vctr(ICVL,IE))>0) then
!          D(ICVL)=D(ICVL) 
!     &     +(2.d0-dble(kdbf(ABS(vctr(ICVL,IE)))))/
!     &      (2.d0-dble(kdbf(ABS(vctr(ICVL,IE))))+SML)
!     &     *tmpfac(ABS(vctr(ICVL,IE)),1)
!        ENDIF
!        enddo
!        enddo
!
      ENDDO
!
      DO IE=1,MAXIE   !vctr(ICVL,0) 
!CIDR VECTOR
        DO ICVL=ICVS_V,ICVE_V
!        if(kdbf(ABS(vctr(ICVL,IE)))>=1) cycle 
        IF(ABS(vctr(ICVL,IE))>0) then 
          INL(ICVL)=INL(ICVL)+max(0,(1-kdbf(ABS(vctr(ICVL,IE)))))
        ENDIF
        ENDDO
      ENDDO
!
      DO IE=1,MAXIE 
        DO ICVL=ICVS_V,ICVE_V
        if(vctr(ICVL,IE)==0) cycle 
        if(kdbf(ABS(vctr(ICVL,IE)))>=1) cycle 
        dum1=rho(ICVL)*ccc(ICVL)
        vlp=max(0.d0,rva(ABS(vctr(ICVL,IE)))/dum1)
        vlm=min(0.d0,rva(ABS(vctr(ICVL,IE)))/dum1)
        IF(vctr(ICVL,IE)<0) then  !ICVB 
!CIDR VECTOR 
          IAL(ICVL,IE)=max(0,(1-kdbf(ABS(vctr(ICVL,IE)))))*
     &               msk(LVEDGE((3-(sign(1,vctr(ICVL,IE))))/2,
     &               ABS(vctr(ICVL,IE))))
          aalw(ICVL,IE)=
     &       -(1.d0-dble(kdbf(ABS(vctr(ICVL,IE)))))/
     &        (1.d0-dble(kdbf(ABS(vctr(ICVL,IE))))+SML)
     &        *(tmpfac(ABS(vctr(ICVL,IE)),1))
        elseIF(vctr(ICVL,IE)>0) then   !ICVA 
          IAL(ICVL,IE)=max(0,(1-kdbf(ABS(vctr(ICVL,IE)))))*
     &               msk(LVEDGE((3-(sign(1,vctr(ICVL,IE))))/2,
     &               ABS(vctr(ICVL,IE))))
          aalw(ICVL,IE)=
     &       -(1.d0-dble(kdbf(ABS(vctr(ICVL,IE)))))/
     &        (1.d0-dble(kdbf(ABS(vctr(ICVL,IE))))+SML)
     &        *(tmpfac(ABS(vctr(ICVL,IE)),1))
        ENDIF
        ENDDO
      ENDDO
!
!-------------------------
!--< 2.4 sliding mesh  >--
!-------------------------
!
      IF(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
        INLB=0
        
        IIMATS=0
        IMATS=0
        do 420 nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdsld) then
          do 400 ISLD=1,2
          ISLD2=3-ISLD
          IBFS1=(LBC_INDEX(nb-1)+1)*(2-ISLD)
          IBFE1=(LBC_pair(nb))*(2-ISLD)
          IBFS2=(LBC_pair(nb)+1)*(2-ISLD2)
          IBFE2=(LBC_INDEX(nb))*(2-ISLD2)
!
          do 410 IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)

          dum1=dsqrt(dx*dx+dy*dy+dz*dz)

!          dum1=
!     &        (abs((
!     &            CVCENT(1,ICVB)
!     &           -CVCENT(1,ICVA))*SFAREA(1,ICFL)
!     &          +(CVCENT(2,ICVB)
!     &           -CVCENT(2,ICVA))*SFAREA(2,ICFL)
!     &          +(CVCENT(3,ICVB)
!     &           -CVCENT(3,ICVA))*SFAREA(3,ICFL))+SML)

          dl1=deltt/dum1
!
          IVA=LVEDGE(1,ICFP)
          ICVAO=LVEDGE(1,ICFO)

          avsf =dl1*SFAREA(4,ICFL) 
          avsfn=dl1*SFAREA(4,ICFL)*wi1 
          avsfo=dl1*SFAREA(4,ICFL)*wi2 
!
          INL(ICVA)=INL(ICVA)+1      !1)
          IAA=INL(ICVA)
!
          INL(IVA)=INL(IVA)+1        !2)
          IA=INL(IVA)
!
          INL(ICVA)=INL(ICVA)+1
          IAB=INL(ICVA)
!
          INL(ICVAO)=INL(ICVAO)+1    !3)
          IAO=INL(ICVAO)
!
          D(ICVA)=D(ICVA)+avsf
          if(IAA<=MAXIE.and.IA<=MAXIE) then
            IAL(ICVA,IAA)=IVA
            aalw(ICVA,IAA)=-avsfn
            IAL(IVA,IA)=ICVA
            aalw(IVA,IA)=-avsfn
            D(IVA)=D(IVA)+avsfn !avsf
          else
            lg_IE=.true.
          endif
!
          if(IAB<=MAXIE.and.IAO<=MAXIE) then
            IAL(ICVA,IAB)=ICVAO
            aalw(ICVA,IAB)=-avsfo
            IAL(ICVAO,IAO)=ICVA
            aalw(ICVAO,IAO)=-avsfo
            D(ICVAO)=D(ICVAO)+avsfo
          else
            lg_IE=.true.
          endif          
 410      enddo
 400      enddo
        ENDIF
 420    enddo
      endif
!
!-------------------
! --- check MAXIE 
!-------------------
!
      if(NPE>1) then
        call hpclor(lg_IE)
      endif
!
      if(lg_IE) then
          INLO=0
          do ICV=1,NCV
            INLB=INL(ICV)
            if(INLB>INLO) then
              INLO=INLB
            endif
          enddo
          if(NPE>1) then
            call hpcimax(INLO)
          endif
          if(my_rank==root) then
            write(ifle,*)' ### MSG:MAX. INL(NL)= ',INLO,' MAXIE= ',MAXIE
            write(ifle,*)' ### MSG: ',
     &         'Increasing [MAXIE] and [MAXIE] to :',INLO
          endif
          IEMAX=INLO
          MAXIE=INLO
!
          deallocate(IAU,ipx,IAL,aalw,aaup,vctr,aax,IW2K_VECT,DW2K_VECT)
!
          if(ierr1.ne.0) then
            write(ifle,*) 'ERR: allocating error in solve_poisson'
            call FFRABORT(1,'MSG: Call your FFR supportor')
          endif
          IALLOCA=1
          goto 1000
      endif
!      endif
!
!-------------------------------------
!--< rearrange order to split lower & upper part >--
!-------------------------------------
!
!----------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      if(.false.) then
      IW2K_VECT(:,1:2)=0
      DO 2000 IE=1,MAXIE

      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      if(IAL(ICVL,IE)/=0) then
        IW2K_VECT(IAL(ICVL,IE),1)=1
      ENDIF
      enddo
      enddo
!
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1 
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      if(IAU(ICVL,IE)/=0) then
        IW2K_VECT(IAU(ICVL,IE),1)=1
      ENDIF
      enddo
      enddo

 2000 ENDDO
!
      do 3000 IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1 
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      if(IW2K_VECT(ICVL,1)==0) then
        D(ICVL)=1.d0
        INL(ICVL)=0
        INU(ICVL)=0
        BB(ICVL)=0.d0
        IAL(ICVL,:)=0
        IAU(ICVL,:)=0
        aalw(ICVL,:)=0.d0
        aaup(ICVL,:)=0.d0
      endif
      enddo
!
      do  IE=1,MAXIE
      ICVS=MAT_CVEXT(IIMAT-1)+1 
      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE
      if(IAL(ICVL,IE)==0) cycle
      if(IW2K_VECT(IAL(ICVL,IE),1)==0) then
        IAL(ICVL,IE)=0
        aalw(ICVL,IE)=0.d0
!        INL(ICVL)=INL(ICVL)-1
      endif
      if(IAU(ICVL,IE)==0) cycle
      if(IW2K_VECT(IAU(ICVL,IE),1)==0) then
        IAU(ICVL,IE)=0
        aaup(ICVL,IE)=0.d0
!        INU(ICVL)=INU(ICVL)-1
      endif
      enddo
      enddo
 3000 enddo
!
      endif


      IW(:)=0
      do  IE=1,MAXIE
      IW2K_VECT(:,IE)=0
      DW2K_VECT(:,IE)=0.d0
      enddo
!
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      if(IAL(ICVL,IE).GT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW(ICVL)=IW(ICVL)+1
      ENDIF
      ENDDO
      DO ICVL=1,NCV
      IF(IAL(ICVL,IE).GT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW2K_VECT(ICVL,IW(ICVL))=IAL(ICVL,IE)
        DW2K_VECT(ICVL,IW(ICVL))=aalw(ICVL,IE)
      ENDIF
      ENDDO
      ENDDO
      
!
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      aaup(ICVL,IE)=DW2K_VECT(ICVL,IE)
      IAU(ICVL,IE)=IW2K_VECT(ICVL,IE)
      ENDDO
      ENDDO
      INU(:)=IW(:)
!
      do  IE=1,MAXIE
      IW2K_VECT(:,IE)=0
      DW2K_VECT(:,IE)=0.d0
      enddo
      IW(:)=0
!
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      IF(IAL(ICVL,IE).LT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW(ICVL)=IW(ICVL)+1
      ENDIF
      ENDDO
!
      DO ICVL=1,NCV
      IF(IAL(ICVL,IE).LT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW2K_VECT(ICVL,IW(ICVL))=IAL(ICVL,IE)
        DW2K_VECT(ICVL,IW(ICVL))=aalw(ICVL,IE)
      ENDIF
      ENDDO
      ENDDO
!
      aalw(1:MXALLCV,0:MAXIE)=0.d0
      IAL(1:MXALLCV,1:MAXIE)=0
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      aalw(ICVL,IE)=DW2K_VECT(ICVL,IE)
      IAL(ICVL,IE)=IW2K_VECT(ICVL,IE)
      ENDDO
      ENDDO
      INL(:)=IW(:)
      IW(:)=0
!
!----------------------------
!-< 3. Solve linear system >-
!----------------------------
!
!      D(NCV+1:)=1.d0
!      INL(NCV+1:)=0
!      INU(NCV+1:)=0
!      bb(NCV+1:)=0.d0
!
      
!-----------------------------------------------------
!
!-----------------------------------------------------
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV0(1,0,MXCV,NCV,D)
      endif
! --- 
      do ICVL=1,NCV
        bb(ICVL)=rhs(ICVL)
      enddo
! --- 
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
!----------------------------------------------------
      calsld=.true.!.true.  
      call CMICCG(calsld,IAL,INL,MAXIE,IAU,INU,MAXIE,NCVIN,
     &       MAT_INDEX,MAT_CVEXT,MAT_NO,ical_sld,NMAT,MXMAT,
     &       MXALLCV,MXCV,NCV,N2,LN,NO,IT,IJT,IW,IW2K_VECT(:,1:2)
     &       )
!
!------------------------------------------------------------------
! -----------------------------------------------------------------
!------------------------------------------------------------------
!
      IMODE=1
!
      aepsq=aepscg
      repsq=repscg
      iterq=itercg
!
      cnx='VICCG'
      IMODE1=1   !1
!
      call VICCG(IMODE1,NCOLOR,calsld,IVDIM,
     &           MXALLCV,MXCV,NCV,NCVIN,MAXIE,MAXIE,N2,NO,
     &           MAT_INDEX,MAT_CVEXT,MAT_NO,ical_sld,NMAT,MXMAT,
     &                     aalw,IAL,INL,
     &                     aaup,IAU,INU,
     &                     BB,X,LN,D,Z,
     +                     IT,IJT,W,IW2K_VECT(:,1:2),TEMP,IW,
     &                     iterq,repsq,aepsq,ierr1)
!
!-----------------------------------------------------
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X(1:))
      endif
      
      DO ICVL=1,NCV
       rhs(ICVL)=X(ICVL)    !!!!!!*0.5d0 
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcimax(iterq)
      ENDIF
!
      aepsp=aepsq
      repsp=repsq
      iterp=iterq
!
      if(my_rank==root.and.mod(iter,100)==0) then
        if(iterq==0) then
          write(ifle,'(2a)') ' ### WRN: BCGSTB solver NOT run at ',
     &   'unsymm_vect poisson'
        endif
      endif
! --- --------------------------------------------------------------
      return
!
      write(ifll,6000) cnx,iterq,repsq,aepsq,ierr1
 6000 format('ERR: solve_poisson_unsymm_vect ',a,' : 
     &  iterations=',i10,
     * ' / maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     * ' / return code=',i5)
!
 9999 continue
      write(ifle,*) '(solve_poisson_unsymm_vect)'
      ierror=1
!
      end subroutine solve_poisson_unsymm_vect
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine solve_cnvdif_vect
     & (clrsld,NSSF,ndf,nq,cn,time,ndiag,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,kdbf,rva,tempf,dif,diag,
     &  LCYCOLD,wifsld,OPPANG,
     &  cofd,dscl,htcbnd,rhs,deltt,vctr,
     &  iterph,reps,aeps,iter,ierror,IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This subroutine is for calculating NS equations
      use module_dimension
      use module_metrix,only   : msk
      use module_metrix,only   : IAL
      use module_metrix,only   : IAU=>IW2K2
      use module_metrix,only   : INL
      use module_metrix,only   : INU
      use module_metrix,only   : bb
      use module_metrix,only   : aalw
      use module_metrix,only   : aaup=>W2K1
!
      use module_metrix,only   : dmax=>W1K7
      use module_metrix,only   : dmin=>W1K8
!
      use module_metrix,only   : D  =>W1K6    !W1K1 W1K6
      use module_metrix,only   : X  =>d1work2
      use module_metrix,only   : Z  =>W1K2
      use module_metrix,only   : W  =>W2K2
      use module_metrix,only   : TEMP=>W1K3
!
      use module_metrix,only   : IW=>iwork4  !used in move.f

      use module_metrix,only   : IT=>iwork5
      use module_metrix,only   : IJT=>iwork6
      use module_metrix,only   : LN=>iwork7

      use module_metrix,only   : aax,ipx
      use module_metrix,only   : tmpfac=>d2vect
      use module_metrix,only   : IW2K_VECT,DW2K_VECT
!
! --- [module arguments]
!
      use module_material,only : ical_sld
      use module_model,only    : ical_vect,nthrds,ical_MHD
      use module_material,only : relax_MHD
      use module_hpcutil
      use module_io,only       : ifll,ifle,cntlnam
      use module_cgsolver,only : aepsbcg,repsbcg,iterbcg
      use module_flags,only    : 
     &                           Adams_Moulton,intgvv,Crank_Nicolson
      use module_boundary,only : nbcnd,kdbcnd,LBC_INDEX,MAT_BCIDX,
     &                           kdsld,idis,LBC_pair,ical_buff,
     &                           kdovst
      use module_material,only : rotati,end,begin,rot_ang
      use module_species,only  : spcnam
      use module_scalar ,only  : sclname
      use module_metrix,only   : aax,ipx
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f,
     &                           NCOLOR,NO,
     &                           MXUP,MXLW,N2
      use module_time  ,only   : MHD_steady
      use module_material,only : ical_sld
!
      implicit none
!
! --- [dummy arguments]
!
      logical,intent(in)    :: clrsld
      integer,intent(in)    :: NSSF,ndf,nq,IMODE,iter,ndiag
      real*8 ,intent(in)    :: deltt,time
      character,intent(in)  :: cn
      integer,intent(inout) :: kdbf  (  MXCVFAC)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF( MXSSFBC)
      integer,intent(in)    :: LCYCSF(  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV(  MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: rva   (  MXCVFAC)
      real*8 ,intent(inout) :: tempf (  MXCVFAC,2)
      real*8 ,intent(inout) :: dif   (  MXALLCV,ndf)
      real*8 ,intent(in)    :: diag  (  MXALLCV,ndiag)
      real*8 ,intent(in)    :: cofd  (  MXALLCV)
      real*8 ,intent(in)    :: dscl  (  MXALLCV)
      real*8 ,intent(in)    :: htcbnd(  NSSF)
      real*8 ,intent(inout) :: rhs   (  MXALLCV,nq)
      real*8 ,intent(in)    :: rho   (  MXALLCV,2)
      integer,intent(out)   :: ierror,iterph(nq)
      real*8 ,intent(out)   :: reps(nq),aeps(nq)
      integer,intent(in)    :: LCYCOLD(MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG (MXSSFBC_SLD)
      integer,intent(in) :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      character(20) :: cnx
!
      real*8 ,parameter   :: XHALF=0.50D0,SML=1.D-15,ZERO=0.D0
      real*8  :: dumC,wi1,wi2,dlvect
      real*8  :: dum1,dum2,dum3,dum4,th(2),costh1,sinth1,costh2,sinth2
      integer :: icx,iterq,ierr,ierr1,kl,ku
      real*8  :: dx,dy,dz,alf,alfd,alfo,vlp,vlm,vdd,cof,alfn,dl
      real*8  :: wifa,wifb,vaa,vab,vbb,vba,dl1,dxo,dyo,dzo
      real*8  :: aepsq,repsq,calph,RMAXS
      real*8  :: unit(3,2),fbb(3,3,2),rbb(3,3,2),sf1,sf2,sf3,avsf
      real*8  :: avsfo,dln,dl2,avsfn,alfdn,alfdo,alfe,alfen,
     &           alfeo,vlpN,vlmN,vlpO,vlmO
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,iti
      integer :: ICOM,IMD,ICH,IFLD,IBFL,IBFS,IBFE,IE,IV,ICFO,ICVBO
      integer :: NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ICFS,ICFE,ICFL,nn
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICV,INLO,
     &           IIMATS(2),IMATS(2),ISLD,ISLD2,nl,nu,INLB
      integer :: IALLOCA=0,IQMXCV,IMODE1,icvao,IBFS1,IBFS2,IBFE1,IBFE2
      real*8  :: org_x1,org_y1,org_x2,org_y2,
     &           cofA,cofB,alfA,alfB
      logical :: lg_IE
      logical :: calsld
      integer :: IFLAGG=0
!
      
      calsld=.true.
      ipx(:)=0
      aax(:)=0.d0
!
      ierr=0
      ierror=0
      ierr1=0
      IALLOCA=0
!
      calph=1.0d0
      if(cn.eq.'v') then
        if(intgvv.eq.Adams_Moulton) then
          calph=5.0d0/12.d0
          if( iter.eq.1 ) calph=1.0d0
        elseif(intgvv.eq.Crank_Nicolson) then
          calph=0.5d0
          if( iter.eq.1 ) calph=1.0d0
        endif
      endif
      
!--------------------------
! --- re-allocate array  --
!--------------------------
!
      IALLOCA=0
      lg_IE=.false.
!---------------------------------------------------------
!-< 1. Set maximul & minimul diffusivity for simplicity >-
!---------------------------------------------------------
!      if(IMODE/=3.and.cn/='R') then
      if(IMODE/=3) then
        do IIMAT=1,NMAT   !ICV=1,NALLCV
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
          dx=dif(ICVL,1)
          dy=dx
          do m=2,ndf
            dx=max(dx,dif(ICVL,m))
            dy=min(dy,dif(ICVL,m))
          enddo
          dmax(ICVL)=dx
          dmin(ICVL)=dy
        enddo
!
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        do ICVL=IDCS,IDCE
        dx=dif(ICVL,1)
        dy=dx
        do 102 m=2,ndf
        dx=max(dx,dif(ICVL,m))
        dy=min(dy,dif(ICVL,m))
  102   continue
        dmax(ICVL)=dx
        dmin(ICVL)=dy
        enddo
        enddo
      endif
!2222
!-----------------------------------
!-< 2. Make up coefficient matrix >-
!-----------------------------------
!--< 2.1 all over the domain >--
!-----------------------------------
!
      if(ical_sld/=0) then
        call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
      endif
!
      aalw(1:MXALLCV,0:MAXIE)=0.d0
      aaup(1:MXCV,0:MAXIE)=0.d0
!
!      if(NEW_TSTP) then
        IAL(1:MXALLCV,1:MAXIE)=0
        IAU(1:MXCV,1:MAXIE)=0
        INL(:)=0
        INU(:)=0
        IT(:)=0
        IJT(:)=0
        LN(:)=0
!      endif
!
      D(0)=0.d0
      if(ndiag>1) THEN
        DO ICVL=ICVS_V,ICVE_V
        D(ICVL)=0.d0
        enddo
      ELSE
        D(1:NCV)=diag(1:NCV,1)
      ENDIF
!
!      IW(:)=0
      tmpfac(:,1)=0.d0
      tempf(:,1)=0.d0
      tempf(:,2)=0.d0
!
      do ICFL=ICFS_V,ICFE_V
        IF(kdbf(ICFL)==3) kdbf(ICFL)=1
      enddo
!
      do ICFL=ICFS_V,ICFE_V
        IF(kdbf(ICFL)/=2.and.kdbf(ICFL)/=0.and.kdbf(ICFL)/=1) 
     &  kdbf(ICFL)=3
      enddo
!
!
!
      do 100 ICFL=ICFS_V,ICFE_V
      ICVA=msk(LVEDGE(1,ICFL))
      ICVB=msk(LVEDGE(2,ICFL))
!                dum1=(abs(
!     &           (CVCENT(1,ICVB)-CVCENT(1,ICVA))*SFAREA(1,ICFL)
!     &          +(CVCENT(2,ICVB)-CVCENT(2,ICVA))*SFAREA(2,ICFL)
!     &          +(CVCENT(3,ICVB)-CVCENT(3,ICVA))*SFAREA(3,ICFL))+SML)

      tmpfac(ICFL,1)=SFAREA(4,ICFL)   !alf
     &        /(abs(
     &           (CVCENT(1,ICVB)-CVCENT(1,ICVA))*SFAREA(1,ICFL)
     &          +(CVCENT(2,ICVB)-CVCENT(2,ICVA))*SFAREA(2,ICFL)
     &          +(CVCENT(3,ICVB)-CVCENT(3,ICVA))*SFAREA(3,ICFL))+SML)



!     &        *dum1/(abs(
!     &           (CVCENT(1,ICVB)-CVCENT(1,ICVA))
!     &          +(CVCENT(2,ICVB)-CVCENT(2,ICVA))
!     &          +(CVCENT(3,ICVB)-CVCENT(3,ICVA)))+SML)**2
!     &
 100  enddo
!
      do 200 ICFL=ICFS_V,ICFE_V
      if(kdbf(ICFL)==2) then
        tmpfac(ICFL,1)=0.d0
      endif
 200  enddo
!
      do 300 ICFL=ICFS_V,ICFE_V 
        tempf(ICFL,1)=              !alfd  7777
     &  max(tmpfac(ICFL,1)*
     &  max(dmax(LVEDGE(1,ICFL)),dmax(LVEDGE(2,ICFL))),
     &  dscl(LVEDGE(2,ICFL))*SFAREA(4,ICFL))
 300  enddo
!
      if(IMODE==2) then  !2=>3
        do ICFL=ICFS_V,ICFE_V
        tempf(ICFL,2)=tempf(ICFL,1)              !alfo  7777
        enddo
      else
        do 340 ICFL=ICFS_V,ICFE_V
        tempf(ICFL,2)=              !alfo  7777
     &  tmpfac(ICFL,1)*min(dmin(LVEDGE(1,ICFL)),dmin(LVEDGE(2,ICFL)))
 340    enddo
      endif
!
      do IE=1,MAXIE 
!CIDR VECTOR
      DO ICVL=ICVS_V,ICVE_V
      if(vctr(ICVL,IE)==0) cycle
      IF(vctr(ICVL,IE)<0) then 
        ICFL=-vctr(ICVL,IE)
        KD=kdbf(ICFL)
        D(ICVL)=D(ICVL)
     &     +(3.d0-dble(KD))/(3.d0-dble(KD)+SML)   !ICVB
!     &     *(tempf(ICFL,1)+max(0.d0,rva(ICFL)))*calph
     &     *(tempf(ICFL,1)-min(0.d0,rva(ICFL)))*calph
      elseIF(vctr(ICVL,IE)>0) then
        ICFL=vctr(ICVL,IE)
        KD=kdbf(ICFL)
        D(ICVL)=D(ICVL)
     &     +(3.d0-dble(KD))/(3.d0-dble(KD)+SML)   !ICVA
     &     *(2.d0-dble(KD))/(2.d0-dble(KD)+SML)
     &     *(1.d0-dble(KD))/(1.d0-dble(KD)+SML)
!     &     *(tempf(ICFL,1)-min(0.d0,rva(ICFL)))*calph
     &     *(tempf(ICFL,1)+max(0.d0,rva(ICFL)))*calph
      ENDIF
      ENDDO
      ENDDO
!????
!
      DO IE=1,MAXIE
      do 500 ICVL=ICVS_V,ICVE_V
      if(vctr(ICVL,IE)==0) cycle
      IF(vctr(ICVL,IE)<0) then                    !ICVB
        ICFL=-vctr(ICVL,IE)
        KD=kdbf(ICFL)
        aalw(ICVL,IE)=
     &  max(0.d0,dble(1-KD))
!     &   *(-tempf(ICFL,2)-max(0.d0,rva(ICFL)))
     &   *(-tempf(ICFL,2)+min(0.d0,rva(ICFL)))    !OK
     &         *cofd(LVEDGE(2,ICFL))*calph
      elseIF(vctr(ICVL,IE)>0) then                !ICVA
        ICFL=vctr(ICVL,IE)
        KD=kdbf(ICFL)
        aalw(ICVL,IE)=
     &  max(0.d0,dble(1-KD))
!     &  *(-tempf(ICFL,2)+min(0.d0,rva(ICFL)))
     &  *(-tempf(ICFL,2)-max(0.d0,rva(ICFL)))     !OK
     &         *cofd(LVEDGE(1,ICFL))*calph 
      endif
 500  enddo
      enddo
!2222
!
!      if(NEW_TSTP) then
        DO IE=1,MAXIE   !           vctr(ICVL,0)
!CIDR VECTOR
        DO ICVL=ICVS_V,ICVE_V
        IF(ABS(vctr(ICVL,IE))>0) then
          INL(ICVL)=INL(ICVL)+max(0,(1-kdbf(ABS(vctr(ICVL,IE)))))
          IAL(ICVL,IE)=max(0,(1-kdbf(ABS(vctr(ICVL,IE)))))*
     &    msk(LVEDGE((3-(sign(1,vctr(ICVL,IE))))/2,
     &        ABS(vctr(ICVL,IE))))
        ENDIF
        ENDDO
        ENDDO
!      ENDIF
!
!-------------------------
!--< 2.4 sliding mesh  >--
!-------------------------
!-----------------------------------------
!--< 2.2 Clear dummy cell & solid part >--
!-----------------------------------------
      IF(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
        do 421 nb=1,nbcnd
        kd=kdbcnd(0,nb)
        if(kd==kdsld) then
          do 400 ISLD=1,2
          ISLD2=3-ISLD
          IBFS1=(LBC_INDEX(nb-1)+1)*(2-ISLD)
          IBFE1=(LBC_pair(nb))*(2-ISLD)
          IBFS2=(LBC_pair(nb)+1)*(2-ISLD2)
          IBFE2=(LBC_INDEX(nb))*(2-ISLD2)
          
          do 410 IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
!
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(1,ICFP)
!
          ICVAO=LVEDGE(1,ICFO)
!
          wi1=wifsld(IBFL)     !wiface(ICFL)
          wi2=1.d0-wi1
!
          vlp=max(0.d0,rva(ICFL))
          vlm=min(0.d0,rva(ICFL))
!
          vlpN=wi1*max(0.d0,rva(ICFP))
          vlmN=wi1*min(0.d0,rva(ICFP))
          vlpO=wi2*max(0.d0,rva(ICFO))
          vlmO=wi2*min(0.d0,rva(ICFO))
!
          INLB=INL(ICVB)+1
!          INLB=0
          INLO=INL(ICVA)+1
          if(INLB<=MAXIE.and.INLO<=MAXIE) then

            INL(ICVB)=INLB
            IAL(ICVB,INLB)=ICVA
            INL(ICVA)=INLO
            IAL(ICVA,INLO)=ICVB
!
            D(ICVA)=D(ICVA)+(-vlm       )*calph         !!!okabe
            D(ICVB)=D(ICVB)+(     +vlpn )*calph         !!!okabe
!!okabe            D(ICVA)=D(ICVA)+(-vlm +vlp)*wi1      !!!okabe
!!okabe            D(ICVB)=D(ICVB)+(-vlm +vlp)*wi1      !!!okabe
!!okabe            aalw(ICVB,INLB)=-vlp*wi1             !!!okabe
!!okabe            aalw(ICVA,INLO)=+vlm*wi1             !!!okabe
!            aalw(ICVB,INLB)=+vlm                       !!!okabe
            aalw(ICVA,INLO)=-vlp                        !!!okabe
          endif
!
          INLB=INL(ICVAO)+1
!          INLB=0
          INLO=INL(ICVA)+1
          if(INLB<=MAXIE.and.INLO<=MAXIE) then
            INL(ICVAO)=INLB
            IAL(ICVAO,INLB)=ICVA
            INL(ICVA)=INLO
            IAL(ICVA,INLO)=ICVAO
!     
!            D(ICVA)=D(ICVA)  +(vlp-vlm)*calph          !!!okabe
            D(ICVAO)=D(ICVAO)+(vlpo   )*calph           !!!okabe
!!okabe            D(ICVA)=D(ICVA)  +(vlp-vlm)*wi2      !!!okabe
!!okabe            D(ICVAO)=D(ICVAO)+(vlp-vlm)*wi2      !!!okabe
!!okabe            aalw(ICVAO,INLB)=-vlp*wi2            !!!okabe
!!okabe            aalw(ICVA,INLO) =+vlm*wi2            !!!okabe
!            aalw(ICVAO,INLB)=-vlp                      !!!okabe
            aalw(ICVA,INLO) =+vlm                       !!!okabe
          endif


 410      enddo
 400      enddo
        ENDIF
 421  enddo
!
!-------------------
! --- check MAXIE 
!-------------------
!
      endif
!-------------------------------------------------------
!--< rearrange order to split lower & upper part >--
!-------------------------------------------------------
      IW(:)=0
      do  IE=1,MAXIE
      IW2K_VECT(:,IE)=0
      DW2K_VECT(:,IE)=0.d0
      enddo
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      if(IAL(ICVL,IE).GT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW(ICVL)=IW(ICVL)+1
      ENDIF
      ENDDO
      DO ICVL=1,NCV
      IF(IAL(ICVL,IE).GT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW2K_VECT(ICVL,IW(ICVL))=IAL(ICVL,IE)
        DW2K_VECT(ICVL,IW(ICVL))=aalw(ICVL,IE)
      ENDIF
      ENDDO
      ENDDO
!
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      aaup(ICVL,IE)=DW2K_VECT(ICVL,IE)
      IAU(ICVL,IE)=IW2K_VECT(ICVL,IE)
      ENDDO
      ENDDO
      INU(:)=IW(:)
!
      do  IE=1,MAXIE
      IW2K_VECT(:,IE)=0
      DW2K_VECT(:,IE)=0.d0
      enddo
      IW(:)=0
!
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      IF(IAL(ICVL,IE).LT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW(ICVL)=IW(ICVL)+1
      ENDIF
      ENDDO
!
      DO ICVL=1,NCV
      IF(IAL(ICVL,IE).LT.ICVL.and.IAL(ICVL,IE)/=0) then
        IW2K_VECT(ICVL,IW(ICVL))=IAL(ICVL,IE)
        DW2K_VECT(ICVL,IW(ICVL))=aalw(ICVL,IE)
      ENDIF
      ENDDO
      ENDDO
!
      aalw(1:MXALLCV,0:MAXIE)=0.d0
      IAL(1:MXALLCV,1:MAXIE)=0
      DO IE=1,MAXIE
      DO ICVL=1,NCV
      aalw(ICVL,IE)=DW2K_VECT(ICVL,IE)
      IAL(ICVL,IE)=IW2K_VECT(ICVL,IE)
      ENDDO
      ENDDO
      INL(:)=IW(:)
      IW(:)=0
! 
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV0(1,0,MXCV,NCV,D)
      endif
!
!-------------------------------------------------------
!
!      IW2K_VECT(:,1:2)=0
!      DO 2000 IE=1,MAXIE

!      do IIMAT=1,NMAT
!      ICVS=MAT_CVEXT(IIMAT-1)+1
!      ICVE=MAT_INDEX(IIMAT)
!      do ICVL=ICVS,ICVE
!      if(IAL(ICVL,IE)/=0) then
!        IW2K_VECT(IAL(ICVL,IE),1)=1
!      ENDIF
!      enddo
!      enddo
!
!      do IIMAT=1,NMAT
!      ICVS=MAT_CVEXT(IIMAT-1)+1 
!      ICVE=MAT_INDEX(IIMAT)
!      do ICVL=ICVS,ICVE
!      if(IAU(ICVL,IE)/=0) then
!        IW2K_VECT(IAU(ICVL,IE),1)=1
!      ENDIF
!      enddo
!      enddo

! 2000 ENDDO
!
!      do 3000 IIMAT=1,NMAT
!      ICVS=MAT_CVEXT(IIMAT-1)+1 
!      ICVE=MAT_CVEXT(IIMAT)
!      do ICVL=ICVS,ICVE
!      if(IW2K_VECT(ICVL,1)==0) then
!        D(ICVL)=1.d0
!        INL(ICVL)=0
!        INU(ICVL)=0
!        BB(ICVL)=0.d0
!        IAL(ICVL,:)=0
!        IAU(ICVL,:)=0
!        aalw(ICVL,:)=0.d0
!        aaup(ICVL,:)=0.d0
!      endif
!      enddo
!
!      do  IE=1,MAXIE
!      ICVS=MAT_CVEXT(IIMAT-1)+1 
!      ICVE=MAT_CVEXT(IIMAT)
!      do ICVL=ICVS,ICVE
!      if(IAL(ICVL,IE)==0) cycle
!      if(IW2K_VECT(IAL(ICVL,IE),1)==0) then
!        IAL(ICVL,IE)=0
!        aalw(ICVL,IE)=0.d0
!      endif
!      if(IAU(ICVL,IE)==0) cycle
!      if(IW2K_VECT(IAU(ICVL,IE),1)==0) then
!        IAU(ICVL,IE)=0
!        aaup(ICVL,IE)=0.d0
!      endif
!      enddo
!      enddo
! 3000 enddo
!-------------------------------------------------------
      calsld=.true.
      call CMICCG(calsld,IAL,INL,MAXIE,IAU,INU,MAXIE,NCVIN,
     &       MAT_INDEX,MAT_CVEXT,MAT_NO,ical_sld,NMAT,MXMAT,
     &       MXALLCV,MXCV,NCV,N2,LN,NO,IT,IJT,IW,IW2K_VECT(:,1:2)
     &       )
!
      if(ndiag>1) THEN
        DW2K_VECT(1:NCV,1)=D(1:NCV)
      endif
!
!----------------------------
!-< 3. Solve linear system >-
!----------------------------
!
      do 600 m=1,nq
!
      X(:)=0.d0
      if(ndiag>1) THEN
        nn=min(m,ndiag)
        do ICVL=1,NCV
        D(ICVL)=DW2K_VECT(ICVL,1)+diag(ICVL,nn)
        enddo
        if(ndiag>1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,D)
        endif
      endif
!
      do ICVL=1,NCV
      bb(ICVL)=rhs(ICVL,m)
      enddo
!
! --- 
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
!------------------------------------------------------------------
! --- 
!------------------------------------------------------------------
!
      aepsq=aepsbcg
      repsq=repsbcg
      iterq=iterbcg
!
      cnx='VBiCG'
      IMODE1=2     !2
!
      calsld=.true.
      call VICCG_B(m,cn,nq,ndiag,IMODE1,NCOLOR,calsld,IVDIM,
     &           MXALLCV,MXCV,NCV,NCVIN,MAXIE,MAXIE,N2,NO,
     &           MAT_INDEX,MAT_CVEXT,MAT_NO,ical_sld,NMAT,MXMAT,
     &                     aalw,IAL,INL,
     &                     aaup,IAU,INU,
     &                     BB,X,LN,D,Z,
     +                     IT,IJT,W,IW2K_VECT(:,1:2),TEMP,IW,
     &                     iterq,repsq,ierr1)
!
!
      IF(NPE.gt.1) THEN
        CALL hpcimax(iterq)
      ENDIF
!
      iterph(m)=iterq
      aeps(m)=aepsq
      reps(m)=repsq
!
      if(my_rank==root.and.mod(iter,100)==0) then
        if(iterq==0) then
          if(cn=='h') then
            write(ifle,'(2a)') ' ### WRN: BCGSTB solver NOT run at ',
     &      'enthalpy equation'
          elseif(cn=='k') then
            write(ifle,'(2a,I4,2x,a)')
     &     ' ### WRN: BCGSTB solver NOT run at ',
     &      'scalar equation no= ',m,trim(sclname(m))
          elseif(cn=='y') then
            write(ifle,'(2a,I4,2x,a)') 
     &    ' ### WRN: BCGSTB solver NOT run at ',
     &      'species equation no= ',m,trim(spcnam(m))
          elseif(cn=='v') then
            write(ifle,'(2a,I4)') 
     &    ' ### WRN: BCGSTB solver NOT run at ',
     &      'velocity equation no= ',m
          endif
        write(ifle,'(2a)')
     &     ' ### MSG: Reduce [aepsbcg] and [repsbcg] in ',
     &  trim(cntlnam)
        endif
      endif
!
! --- 
!
      if(iterq/=0) then
        do ICVL=1,NCV
        rhs(ICVL,m)=X(ICVL)*0.5d0 
        enddo
      endif
!
      if(ierr1.ne.0) then
        call mknam(m,cn,iti,cnx)
        write(ifll,6000) cnx(:iti),iterq,repsq,aepsq,ierr1
        write(ifle,*) '### error : bicgstab not converged'
        goto 9999
      endif
!
  600 continue
!
      return
 6000 format('*** bicg/',a,' : iterations=',i10,
     * ' / maximul r.&a. error=',1pe15.7,' &',1pe15.7,
     * ' / return code=',i5)
!
 9999 continue
      write(ifle,*) '(solve_cnvdif_vect)'
      ierror=1
!!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine mknam(m,cn,j,cnx)
!=================================================
      integer     ,intent(in)  :: m
      character(*),intent(in)  :: cn
      integer     ,intent(out) :: j
      character(*),intent(out) :: cnx
      integer :: i
!
      cnx=cn
      if(nq.gt.1) write(cnx(11:),'(i10)') m
      j=0
      do 100 i=1,len(cnx)
      if(cnx(i:i).ne.' ' ) then
        j=j+1
        cnx(j:j)=cnx(i:i)
      endif
  100 continue
      end subroutine mknam
!
      end subroutine solve_cnvdif_vect
