!
!      subroutine metric_admin()
!      subroutine metric_3ds()
!      subroutine metric_CV
!      subroutine metric_CV_cell
!      subroutine Merge_SSF
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_admin
     & (mcell,mvrtx,mface,medge,
     &  ncell,nvrtx,nface,nedge,nssfbc,nssfbc_old,NCV,
     &  mssfbc,NCVFAC,NALLCV,IEMAX,NMAT,ncelb,
     &  lcface,lvface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,LBCSSF,LCYCSF,listbc,listpr,LMAT,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,lvcell,lfcell,
     &  time,deltt,ierror,icode)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only       : ifle,ifll,cntlnam
      use module_movegrid,only : ngrid,dtmin,grdtim=>time
      use module_model,only    : vertex_cen,cell_cen,icon_cv
!
! 1. Read vertex file
!
! 2. Calculate metrics
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)       :: mface,mcell,mvrtx,nvrtx,nface,ncell
      integer,intent(in)       :: ncelb
      integer,intent(in)       :: medge,nedge,NCV,mssfbc,icode
      integer,intent(out)      :: NMAT,nssfbc_old
      integer,intent(in)       :: lcface(2,mface)
      integer,intent(in)       :: lvface(4,mface)
      integer,intent(in)       :: lbface(2,mface)
      integer,intent(in)       :: lacell(  mcell)
      integer,intent(out)      :: LMAT  (  mvrtx)
      integer,intent(in)       :: LEFACE(5,mface)
      integer,intent(in)       :: LVRTCV(  mvrtx)
      integer,intent(inout)    :: LVEDGE(2,medge)
      integer,intent(out)      :: LBCSSF( mssfbc)
      integer,intent(in)       :: listpr(0:3,mvrtx)
      integer,intent(out)      :: listbc(4,mssfbc)
      integer,intent(out)      :: LCYCSF( mssfbc)
      integer,intent(in)       :: lvcell(8,mcell)
      integer,intent(in)       :: lfcell(7,mcell)
      real*8 ,intent(inout)    :: cord  (3,mvrtx)
!
      real*8 ,intent(out)      :: area  (4,mface)
      real*8 ,intent(out)      :: volume(  mcell)
      real*8 ,intent(out)      :: gface (3,mface)
      real*8 ,intent(out)      :: gcell (3,mcell)
      real*8 ,intent(out)      :: SFAREA(4,medge)
      real*8 ,intent(out)      :: SFCENT(3,medge)
      real*8 ,intent(out)      :: CVVOLM(  mvrtx)
      real*8 ,intent(out)      :: CVCENT(3,mvrtx)
      real*8 ,intent(in)       :: time,deltt
      integer,intent(out)      :: ierror,nssfbc,NCVFAC,NALLCV,IEMAX
!
! --- [local entities]
!
      logical :: static=.false.
      integer :: lgrid=-1,ihpc=1,lacell0,ic,IV,IS,IIS
      real*8  :: dtgrd=0.d0
      integer :: i,n,icyc,lgd,lgdx,ierr1,ierr2
      real*8  :: tim1,tcyc,tgd,dt1,dt2,dum1
!
! --- check lacell/=0
!
      lacell0=0
      do ic=1,ncell
      if(lacell(ic).eq.0) then
        lacell0=lacell0+1
        write(ifle,'(a,I8)')
     &  ' ### ERR: NOT defined material number: ',ic
      endif
      enddo
      if(lacell0.gt.0) then
        write(ifle,'(a)') 
     &  ' ### ERR: some cell NOT defined material Number'
        write(ifle,'(a,I8)')
     &  ' ### ERR: NOT defined material : ',lacell0
        stop 'metric_admin'
      endif
!
      ierror=0
      ierr1=0
      if(static.and.ihpc.eq.0) return
      if(ngrid.lt.2) static=.true.
      if(deltt.gt.dtmin) goto 9001
!
!      if(lgrid.lt.0) then
!        allocate( cordsq  (3,mvrtx,2),stat=ierr1 )
!        if( ierr1.ne.0 ) then
!          write(ifle,*) '### error : allocating cordsq(:,:,:)'
!          goto 9002
!        endif
!        allocate( areasq  (4,mface,2),stat=ierr1 )
!        if( ierr1.ne.0 ) then
!          write(ifle,*) '### error : allocating areasq(:,:,:)'
!          goto 9002
!        endif
!        allocate( volumesq(  mcell,2),stat=ierr1 )
!        if( ierr1.ne.0 ) then
!          write(ifle,*) '### error : allocating volumesq(:,:)'
!          goto 9002
!        endif
!        allocate( gfacesq (3,mface,2),stat=ierr1 )
!        if( ierr1.ne.0 ) then
!          write(ifle,*) '### error : allocating gfacesq(:,:,:)'
!          goto 9002
!        endif
!        allocate( gcellsq (3,mcell,2),stat=ierr1)
!        if( ierr1.ne.0 ) then
!          write(ifle,*) '### error : allocating gcellsq(:,:,:)'
!          goto 9002
!        endif
!
!        cordsq  =0.d0
!        areasq  =0.d0
!        volumesq=0.d0
!        gfacesq =0.d0
!        gcellsq =0.d0
!      endif
!
!-< 1. Search present position >-
!
                                         icyc=0
      if( grdtim(ngrid+1).gt.grdtim(1) ) icyc=1
!
      tgd=time
      if(icyc.gt.0) then
        tim1=grdtim(1)
        tcyc=grdtim(ngrid+1)-tim1
        tgd=tgd-tim1
        tgd=tgd-dble(int(tgd/tcyc))*tcyc
        if(tgd.lt.0.d0) tgd=tgd+tcyc
        tgd=tgd+tim1
      endif
      lgdx=icyc
      do 100 n=icyc+1,ngrid
      if( grdtim(n).lt.tgd ) lgdx=n
  100 continue
      lgd=max(1,min(ngrid+icyc-1,lgdx))
!
      dt1=0.d0
      dt2=1.d0
      if(lgdx==lgrid.and.ihpc==0) goto 1001
!
!
!-< 2. Calculate metrics >-
!
!--< 2.1 read vertex file >--
!
!      do 210 IV=1,nvrtx
!        do 211 i=1,3
!          cordsq(i,IV,1)=cord(i,IV)
!          cordsq(i,IV,2)=cord(i,IV)
!  211   continue
!  210 continue
!
!--< 2.2 calculate area, volume & center >--
!
!      do 200 n=1,2
      

      write(ifll,'(a,I8)') '###    LIST CELL ...',n
      call metric_3ds
     & (mcell,mvrtx,mface,ncell,nface,
     &  lcface,lvface,lbface,lacell,cord,
     &  area,volume,gface,gcell,
     &  lvcell,lfcell,ierr1)
      if( ierr1.ne.0 ) goto 9999
! 200  continue
!
!--< 2.3 calculate metrics of grid speed >--
!
      dt1=max(0.d0,min(1.d0,dtgrd/deltt))
      dt2=1.d0-dt1
!
 1001 continue
      lgrid=lgdx
      dtgrd=grdtim(lgdx+1)-tgd
!
!-< 5. Set metrics at present time >-
!
!
!--< 5.2 interpolation in time >--
!
      dt2=max(0.d0,min(1.d0,
     &    (tgd-grdtim(lgd))/(grdtim(lgd+1)-grdtim(lgd)) ))
      dt1=1.d0-dt2
!
!      do 510 IV=1,nvrtx
!      do 511 i=1,3
!      cord(i,IV)=dt1*cordsq(i,IV,1)+dt2*cordsq(i,IV,2)
!  511 continue
!  510 continue
!
!      do 520 IC=1,ncell
!      volume(IC)=dt1*volumesq(IC,1)+dt2*volumesq(IC,2)
!      do 521 i=1,3
!      gcell(i,IC)=dt1*gcellsq(i,IC,1)+dt2*gcellsq(i,IC,2)
!  521 continue
!  520 continue
!
!      do 530 IS=1,nface
!      do 531 i=1,4
!      area(i,IS)=dt1*areasq(i,IS,1)+dt2*areasq(i,IS,2)
!  531 continue
!      do 532 i=1,3
!      gface(i,IS)=dt1*gfacesq(i,IS,1)+dt2*gfacesq(i,IS,2)
!  532 continue
!  530 continue
!------------------------
! --- CV Metrix: 
!------------------------ 
      write(ifll,'(a)') '###    METRIC_CV ...'
      if(icon_cv==vertex_cen) then
        call metric_CV
     & (mcell,mvrtx,mface,medge,mssfbc,ncell,nface,nedge,
     &  NCV,nssfbc,NCVFAC,NALLCV,IEMAX,nssfbc_old,
     &  lcface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,LBCSSF,LCYCSF,listbc,listpr,LMAT,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,
     &  ierr2)
      else
        call metric_CV_cell
     & (mcell,mvrtx,mface,medge,mssfbc,ncell,nface,nedge,
     &  NCV,nssfbc,NCVFAC,NALLCV,nssfbc_old,ncelb,
     &  lcface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,LBCSSF,LCYCSF,listbc,listpr,LMAT,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,
     &  ierr2)
      endif
!------------------------
! --- Merge SSF on BC 
!------------------------
      write(ifll,'(a)') '###    MERGE_SSF ...' 
      if(icon_cv==vertex_cen) then
        call Merge_SSF(mssfbc,medge,nedge,mvrtx,
     &           NCV,nssfbc,nssfbc_old,NCVFAC,NALLCV,
     &           LMAT,LBCSSF,LCYCSF,LVEDGE,listbc,SFCENT,SFAREA)
      endif
!
      IEMAX=0
      call list_edgvtx(2,medge,NCVFAC,NALLCV,LVEDGE,IEMAX)
!---------------------------------------------------------------------
! --- Match "touch-inlet","interface" and "Sliding" BC: "LCYCSF(:,:)"
!---------------------------------------------------------------------
      write(ifll,'(a)') '###    LIST TOUCH_INLET ...' 
      write(ifll,'(a)') '###    LIST FLUID/SOLID INTERFACE ...' 
      call list_touch_inlet(mssfbc,NSSFBC,medge,nedge,mvrtx,
     &                      LVEDGE,LMAT,LCYCSF,LBCSSF,SFCENT)
!------------------------
! --- set dummy cells >--
!------------------------
      write(ifll,'(a)') '###    dc_metric (After merge SSF)...'
      call dc_metric
     & (NALLCV,NCVFAC,NSSFBC,nedge,medge,mssfbc,mvrtx,
     &  LVEDGE,LBCSSF,LCYCSF,SFCENT,SFAREA,CVCENT,CVVOLM,2,icode)
!-------------------------------------------------
!-< 6. Deallocate arrays in case of static grid >-
!-------------------------------------------------
!      if( static.and.ihpc.eq.0 )
!     & deallocate( cordsq,areasq,volumesq,gfacesq,gcellsq)
!
      return
!
 9001 continue
      write(ifle,'(a)') '### error : data error'
      write(ifle,'(2a)') 
     &         'time increment is bigger than minimum ',
     &         'time interval of grid tables'
      write(ifle,'(a,F14.4)') 'time increment = ',deltt
      write(ifle,'(a,F14.4)') 'minimum time interval = ',dtmin
      goto 9999
 9002 continue
      write(ifle,'(a)') '### error : allocation failed'
      goto 9999
 9999 continue
      write(ifle,'(a)') '(metric_admin)'
      ierror=1
!
      end subroutine metric_admin
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_CV
     & (mcell,mvrtx,mface,medge,mssfbc,ncell,nface,nedge,
     &  NCV,nssfbc,NCVFAC,NALLCV,IEMAX,nssfbc_old,
     &  lcface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,LBCSSF,LCYCSF,listbc,listpr,LMAT,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,
     &  ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only          : ifle,ifll,cntlnam
      use module_partitioner,only : lvrtx=>IW1
      use module_partitioner,only : ip   =>WK1
      use module_partitioner,only : nclv =>WK2
      use module_boundary,only    : kdbcnd,kdintr,kdprdc,kdsld,idis,
     &                              kdbuff,kdshutr,kdpors,kdovst
      use module_model,only    : ical_mvmsh
!
      implicit none
!
! 1.  Calculate metrics for each face & cell
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mcell,mvrtx,mface,medge,mssfbc
      integer,intent(in)     :: ncell,nface,nedge,NCV
      integer,intent(in)     :: lcface(2,mface)
      integer,intent(in)     :: lbface(2,mface)
      integer,intent(in)     :: lacell(  mcell)
      integer,intent(in)     :: LEFACE(5,mface)
      integer,intent(in)     :: listpr(0:3,mvrtx)
      integer,intent(out)    :: listbc(4,mssfbc)
      integer,intent(inout)  :: LVEDGE(2,medge)
      integer,intent(in)     :: LVRTCV(  mvrtx)
      integer,intent(out)    :: LBCSSF( mssfbc)
      integer,intent(out)    :: LCYCSF( mssfbc)
      integer,intent(out)    :: LMAT  (  mvrtx)
!      
      real*8 ,intent(in)     :: cord  (3,mvrtx)
      real*8 ,intent(in)     :: area  (4,mface)
      real*8 ,intent(in)     :: volume(  mcell)
      real*8 ,intent(in)     :: gface (3,mface)
      real*8 ,intent(in)     :: gcell (3,mcell)
      real*8 ,intent(out)    :: SFAREA(4,medge)
      real*8 ,intent(out)    :: SFCENT(3,medge)
      real*8 ,intent(out)    :: CVVOLM(  mvrtx)
      real*8 ,intent(out)    :: CVCENT(3,mvrtx)
      integer,intent(out)    :: ierror,nssfbc,NCVFAC,NALLCV,IEMAX,
     &                          nssfbc_old
!
! --- [local entities]
!
      integer :: idmax,LMAT0,SMLV
      integer :: KMAT(2),IDC(2)
      integer :: i,j,k,l,m,n,kd,kdv,kdt,kdy,kdk,kdp,nbp,ISP
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ILV,ILS,ILE,ILCV,IEE
!
      real*8  :: gfc(3),gc(3),gec(3)
      real*8  :: gfc1(3),gc1(3),gec1(3),gfc2(3),gc2(3),gec2(3)
      real*8  :: DIRC,SSFS,SSV
      real*8  :: SSF(3),GSSF(3),GSSA(3),GSSB(3)
      real*8  :: GSSA1(3),GSSB1(3)
      real*8  :: GSSA2(3),GSSB2(3)
      real*8  :: VEDGE1(3),VFE1(3),VCE1(3)
      real*8  :: VEDGE2(3),VFE2(3),VCE2(3)
      real*8  :: VCE(3),VEDGE(3),VFE(3),UEDGE(3)
      real*8  :: gvA(3),gvB(3),VCEA(3),VCEB(3)
      real*8  :: gvA2(3),gvB2(3),VCEA2(3),VCEB2(3)
      real*8  :: gvA1(3),gvB1(3),VCEA1(3),VCEB1(3)
      real*8  :: SSFA(3),SSFB(3),SSFSA,SSFSB
      real*8  :: SSFA1(3),SSFB1(3),SSFSA1,SSFSB1
      real*8  :: SSFA2(3),SSFB2(3),SSFSA2,SSFSB2
      real*8  :: SMSSFX,SMSSFY,SMSSFZ,SUMSSF,SMSSF,SUMSSV
      real*8  :: SUMSWX,SUMSWY,SUMSWZ,GGG1,GGG2,GGG3
      real*8  :: SSVXA,SSVYA,SSVZA,SSVXB,SSVYB,SSVZB,DIRC1,DIRC2
      real*8  :: dum1,dum2
      logical :: ICYC
      integer :: ierr1=0,IE2,ICVBP2,ICFA1,ICFB1,ICFA2,ICFB2,IDCA1,IDCB1
      integer :: ICFA,ICFB,IDCA,IDCB,ISSFA,ISSFB,ILE1,IE1,ICVA1,ICVB1,
     &          IVA1,IVB2,IVB1,ILE2,ICVA2,ICVB2,IVA2,IA,IB,IVAP2,IVBP2,
     &           ICVAP2,IDCA2,IDCB2,ISSFA1,ISSFB1,ISSFA2,ISSFB2
      integer :: NDUP=0,NBCSSF
      real*8 ,parameter :: r1p6=1.d0/6.d0
      real*8 ,parameter :: r1p3=1.d0/3.d0
      real*8 ,parameter :: r1p4=0.25d0
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      allocate(lvrtx(4,mface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating lvrtx(:,:) in metric_CV'
      allocate(ip(0:nedge),stat=ierr1) 
      if(ierr1.ne.0) stop'stop at allocating ip(:) in metric_CV'
      allocate(nclv(4*nface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating nclv(:) in metric_CV'
!
!------------------------------------
! --- initializing array of CVCENT
!------------------------------------
!
      CVCENT=0.D0
      CVVOLM=0.D0
      LMAT=0
      LCYCSF=0
      LBCSSF=0
      ICYC=.false.
!
!----------------------------------------------------
! --- CONTROL VOLUME (CV) METRIX CALCULATION
!----------------------------------------------------
! --- make-up egde-face connectivity
!----------------------------------------------------
!
      lvrtx(1:4,:)=LEFACE(1:4,:)
      call list_onvrtx(4,mface,nedge,nface,lvrtx,ip,nclv,idmax)
!
!=====================================================================
! --- <1> define inner Sub-Face (SF) vector (CV face)
!=====================================================================
      DO 400 IE=1,nedge
!
      ICVA=LVEDGE(1,IE)
      ICVB=LVEDGE(2,IE)
      IVA=LVRTCV(ICVA)
      IVB=LVRTCV(ICVB)
!
      SMSSFX=ZERO
      SMSSFY=ZERO
      SMSSFZ=ZERO
!
      SUMSWX=ZERO
      SUMSWY=ZERO
      SUMSWZ=ZERO
!
      SUMSSV=ZERO
      smssf=ZERO
!
      SSVXA=ZERO
      SSVYA=ZERO
      SSVZA=ZERO
!
      SSVXB=ZERO
      SSVYB=ZERO
      SSVZB=ZERO
!
      do 410 IEE=ip(IE-1)+1,ip(IE)
!
! --- ALL CELL-FACE AROUND EDGE(IE)
!
      IS=nclv(IEE)
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      IDC(1)=IC1
      IDC(2)=IC2
      KMAT(1)=lacell(IC1)
      KMAT(2)=lacell(IC2)
      IF(KMAT(1).ne.KMAT(2)) then
	write(ifle,'(a,I12)') 
     &  '### error : two Cells No. Different on face ',IS
        write(ifle,'(a,2I10,2I4)') 
     &  'IC1, IC2, IMAT1, IMAT2 ',IC1,IC2,KMAT(1),KMAT(2)
      else
        IMAT=lacell(IC1)
      endif
      DO 450 J=1,2
      IC=lcface(J,IS)
!
      IF(IC.GT.NCELL) goto 450
!
! --- THREE CENTERS: CELL, CELL-FACE, EDGE
!
      DO 420 I=1,3
      gfc(I)=gface(I,IS)
      gc(I)=gcell(I,IC)
      gec(I)=XHALF*(cord(I,IVA)+cord(I,IVB))
!
! --- HALF EDGE VECTOR: (DIRECTION: A->B)
      VEDGE(I)=XHALF*(cord(I,IVB)-cord(I,IVA))
!
! --- THREE VECTORS OF TWO SUB-SUB-FACE (SSF)
!
! --- edge-center TO cell-face-center
      VFE(I)=gfc(I)-gec(I)
!
! --- edge-center TO cell-center
      VCE(I)=gc(I)-gec(I)
!
! --- CENTER OF SSF
!
      GSSF(I)=r1p3*(gfc(I)+gc(I)+gec(I))
!
! --- CENTER OF SSVA AND SSVB
!
      GSSA(I)=r1p4*(gfc(I)+gc(I)+gec(I)+cord(I,IVA))
      GSSB(I)=r1p4*(gfc(I)+gc(I)+gec(I)+cord(I,IVB))
!
 420  CONTINUE
! ===
      DIRC=DSQRT(VEDGE(1)*VEDGE(1)+VEDGE(2)*VEDGE(2)+VEDGE(3)*VEDGE(3))
      UEDGE(1)=VEDGE(1)/DIRC
      UEDGE(2)=VEDGE(2)/DIRC
      UEDGE(3)=VEDGE(3)/DIRC
!
! --- SSF vector: SSF
!
      CALL AXBEQC(VFE,VCE,SSF,SSFS)
! === 
      DIRC=SSF(1)*UEDGE(1)
     &    +SSF(2)*UEDGE(2)
     &    +SSF(3)*UEDGE(3)
!
! === 
!      SSF(1)=DIRC*UEDGE(1)
!      SSF(2)=DIRC*UEDGE(2)
!      SSF(3)=DIRC*UEDGE(3)
!
      IF(DIRC.LT.ZERO) THEN
        SSF(1)=-SSF(1)
        SSF(2)=-SSF(2)
        SSF(3)=-SSF(3)
      ENDIF
! 
      SMSSFX=SMSSFX+SSF(1)
      SMSSFY=SMSSFY+SSF(2)
      SMSSFZ=SMSSFZ+SSF(3)
!
      SUMSWX=SUMSWX+GSSF(1)*SSFS
      SUMSWY=SUMSWY+GSSF(2)*SSFS
      SUMSWZ=SUMSWZ+GSSF(3)*SSFS
!
! --- SUB-SUB-VOLUME: SSV
!
      SSV=r1p3*(SSF(1)*VEDGE(1)
     &         +SSF(2)*VEDGE(2)
     &         +SSF(3)*VEDGE(3))
!
      SUMSSV=SUMSSV+SSV
      smssf=smssf+SSFS
!
      GGG1=GSSA(1)*SSV
      GGG2=GSSA(2)*SSV
      GGG3=GSSA(3)*SSV
      SSVXA=SSVXA+GGG1
      SSVYA=SSVYA+GGG2
      SSVZA=SSVZA+GGG3
!
      GGG1=GSSB(1)*SSV
      GGG2=GSSB(2)*SSV
      GGG3=GSSB(3)*SSV
      SSVXB=SSVXB+GGG1
      SSVYB=SSVYB+GGG2
      SSVZB=SSVZB+GGG3
!
 450  CONTINUE
      IF(IDC(1).gt.ncell.and.IDC(2).gt.ncell) then
	write(ifle,'(a,I10)') 
     &   'ERR: Two Cells are Dumy Cell on face ',IS
        stop
      endif
 410  CONTINUE
!
! --- CV'S SUB-FACE VECTOR AND AREA
!
      SUMSSF=SQRT(SMSSFX*SMSSFX
     &           +SMSSFY*SMSSFY
     &           +SMSSFZ*SMSSFZ)
!
      SFAREA(1,IE)=SMSSFX/(SUMSSF+SML)
      SFAREA(2,IE)=SMSSFY/(SUMSSF+SML)
      SFAREA(3,IE)=SMSSFZ/(SUMSSF+SML)
      SFAREA(4,IE)=SUMSSF
!
! --- SF CENTER
!
      SFCENT(1,IE)=SUMSWX/(SMSSF+SML)  !SMSSF
      SFCENT(2,IE)=SUMSWY/(SMSSF+SML)
      SFCENT(3,IE)=SUMSWZ/(SMSSF+SML)
!
! --- 
!
      if(LMAT(ICVA).ne.0.and.IMAT.ne.LMAT(ICVA)) then
        NDUP=NDUP+1
        write(ifle,'(a,2I10)') 
     &   'ERR: IMAT,LMAT(ICVA)=',IMAT,LMAT(ICVA)
        write(ifle,'(a,I10)') 'ERR: Duplicating Vertex-1 :',ICVA
        write(ifle,*)
      else
        LMAT(ICVA)=IMAT
      endif
!
      CVVOLM(ICVA)=CVVOLM(ICVA)+SUMSSV
      CVCENT(1,ICVA)=CVCENT(1,ICVA)+SSVXA
      CVCENT(2,ICVA)=CVCENT(2,ICVA)+SSVYA
      CVCENT(3,ICVA)=CVCENT(3,ICVA)+SSVZA
! --- 
      if(LMAT(ICVB).ne.0.and.IMAT.ne.LMAT(ICVB)) then
        NDUP=NDUP+1
        write(ifle,'(a,2I10)') 
     &   'ERR: IMAT,LMAT(ICVB)=',IMAT,LMAT(ICVB)
        write(ifle,'(a,I8)') 'ERR: Duplicating Vertex-2 :',ICVB
      else
        LMAT(ICVB)=IMAT        
      endif
      CVVOLM(ICVB)=CVVOLM(ICVB)+SUMSSV
      CVCENT(1,ICVB)=CVCENT(1,ICVB)+SSVXB
      CVCENT(2,ICVB)=CVCENT(2,ICVB)+SSVYB
      CVCENT(3,ICVB)=CVCENT(3,ICVB)+SSVZB
!
 400  CONTINUE
!
      if(NDUP.gt.0) then
        write(ifle,'(a)') 
     & 'MSG: STOP at 400 Loop in subroutine metric_CV'
        write(ifle,'(a,I4)') 
     &      'ERR: Duplicating Vertex number: ',NDUP
        write(ifle,'(a)') 'MSG: Call you FFR supportor'
        stop ': ERR-1: Duplicating Vertex'
      endif
!------------------------
! --- check LMAT
!------------------------
!
      LMAT0=0
      SMLV=0
      do ICV=1,NCV
        if(LMAT(ICV).eq.0) then
          write(ifle,'(a,I8)') 
     &    'ERR: NOT defined material numberICV :',
     &        ICV
          LMAT0=LMAT0+1
        endif
        if(CVVOLM(ICV).lt.SML) then
        SMLV=SMLV+1
        endif
      enddo
      if(LMAT0.gt.0) then
        write(ifle,'(a)') 
     &  'ERR: some cell NOT defined material Number'
        write(ifle,'(a,I4)') 
     &  'ERR: NOT defined material cell no.: ',LMAT0
        stop 'metric_CV-1-CELL'
      endif
      if(SMLV.gt.0) then
        write(ifle,'(a)') 
     &  'ERR: some vrtex NOT belong to any cell'
        write(ifle,'(a,E14.4)') 
     &  'ERR: free vertex no. : ',SMLV
        stop 'metric_CV-1-CELL'
      endif
!
!------------------------
! --- Weighted CV center:
!------------------------
!
      DO 480 IV=1,NCV

      CVCENT(1,IV)=CVCENT(1,IV)/(CVVOLM(IV)+SML)
      CVCENT(2,IV)=CVCENT(2,IV)/(CVVOLM(IV)+SML)
      CVCENT(3,IV)=CVCENT(3,IV)/(CVVOLM(IV)+SML)
 480  CONTINUE
!
!=====================================================================
! --- <2> Defining SSF on BC
!=====================================================================
!
      nssfbc=0
      NCVFAC=nedge
      NALLCV=NCV
!
      DO 500 IS=1,nface
!
      nb=lbface(1,IS)
! --- If not BC face goto 500:
      if(nb==0) goto 500
      kd=kdbcnd(0,abs(nb))
      if(nb.lt.0.and.kd==kdprdc.and.idis(abs(nb))==0) goto 500
      ISP=lbface(2,IS)
      IC2=lcface(2,IS)
      IMAT=lacell(IC2)
!
      if(ISP==0.or.(kd==kdintr).or.(kd==kdbuff).or.(kd==kdshutr).or.
     &             (kd==kdsld).or.(kd==kdpors)) then !kd==kdovst
!---------------------------------------------------------
! --- 1) Inlet, 2) Outlet, 3) symmetric BC, 4) touchinlet
! --- 5) interface, 6) Sliding
!---------------------------------------------------------
        do 510 ILE=1,LEFACE(5,IS)
        IE=LEFACE(ILE,IS)
        ICVA=LVEDGE(1,IE)
        ICVB=LVEDGE(2,IE)
        IVA=LVRTCV(ICVA)
        IVB=LVRTCV(ICVB)
!
        nssfbc=nssfbc+2
        NCVFAC=NCVFAC+2
        NALLCV=NALLCV+2
        if(NCVFAC>medge)  then
          write(*,*) 'ERR: '
          stop 'Increase [medge]  in cntlnam '
        endif
        if(NALLCV>mvrtx)  then
          write(*,*) 'ERR: '
          stop 'Increase [mvrtx]  in cntlnam '
        endif
        if(nssfbc>mssfbc) then
          write(*,*) 'ERR: '
          stop 'Increase [mssfbc] in cntlnam '
        endif
!
        ICFA=NCVFAC-1
        ICFB=NCVFAC
        IDCA=NALLCV-1
        IDCB=NALLCV
        ISSFA=ICFA-nedge
        ISSFB=ICFB-nedge
!
        DO 530 I=1,3
        gfc(I)=gface(I,IS)
        gvA(I)=cord(I,IVA)
        gvB(I)=cord(I,IVB)
        gec(I)=XHALF*(cord(I,IVA)+cord(I,IVB))
! --- Cell Face Vector:
        VEDGE(I)=area(I,IS)
! --- edge-center TO cell-face-center
        VFE(I)=gfc(I)-gec(I)
! --- edge-center TO A & B
        VCEA(I)=gvA(I)-gec(I)
        VCEB(I)=gvB(I)-gec(I)
! --- CENTER OF SSF A AND B
        GSSA(I)=r1p3*(gfc(I)+gvA(I)+gec(I))
	GSSB(I)=r1p3*(gfc(I)+gvB(I)+gec(I))
!
 530    CONTINUE
!
! --- SSF vector: SSF
!
        CALL AXBEQC(VFE,VCEA,SSFA,SSFSA)
        DIRC=SSFA(1)*VEDGE(1)
     &      +SSFA(2)*VEDGE(2)
     &      +SSFA(3)*VEDGE(3)
        IF(DIRC.LT.ZERO) THEN
          SSFA(1)=-SSFA(1)
          SSFA(2)=-SSFA(2)
          SSFA(3)=-SSFA(3)
        ENDIF
!
        CALL AXBEQC(VFE,VCEB,SSFB,SSFSB)
        DIRC=SSFB(1)*VEDGE(1)
     &      +SSFB(2)*VEDGE(2)
     &      +SSFB(3)*VEDGE(3)
        IF(DIRC.LT.ZERO) THEN
          SSFB(1)=-SSFB(1)
          SSFB(2)=-SSFB(2)
          SSFB(3)=-SSFB(3)
        ENDIF
!      
        SFAREA(4,ICFA)=SSFSA
        SFAREA(1,ICFA)=SSFA(1)/SSFSA
        SFAREA(2,ICFA)=SSFA(2)/SSFSA
        SFAREA(3,ICFA)=SSFA(3)/SSFSA
!
        SFAREA(4,ICFB)=SSFSB
        SFAREA(1,ICFB)=SSFB(1)/SSFSB
        SFAREA(2,ICFB)=SSFB(2)/SSFSB
        SFAREA(3,ICFB)=SSFB(3)/SSFSB
!
        SFCENT(1,ICFA)=GSSA(1)
        SFCENT(2,ICFA)=GSSA(2)
        SFCENT(3,ICFA)=GSSA(3)
!
        SFCENT(1,ICFB)=GSSB(1)
        SFCENT(2,ICFB)=GSSB(2)
        SFCENT(3,ICFB)=GSSB(3)
!
        LVEDGE(1,ICFA)=ICVA
        LVEDGE(2,ICFA)=IDCA
        LVEDGE(1,ICFB)=ICVB
        LVEDGE(2,ICFB)=IDCB
!
        LBCSSF(ISSFA)=nb
        listbc(1,ISSFA)=IE
        listbc(2,ISSFA)=IS
        listbc(3,ISSFA)=IVA
        LBCSSF(ISSFB)=nb
        listbc(1,ISSFB)=IE
        listbc(2,ISSFB)=IS
        listbc(3,ISSFB)=IVB
!
        if(LMAT(IDCA).ne.0.and.IMAT.ne.LMAT(IDCA)) then
          NDUP=NDUP+1
          write(ifle,'(a,I10)') 
     &     'ERR: BC:',nb
          write(ifle,'(a,I10)') 'ERR: Duplicating Vertex-3 :',IDCA
        else
          LMAT(IDCA)=IMAT
        endif
        if(LMAT(IDCB).ne.0.and.IMAT.ne.LMAT(IDCB)) then
          NDUP=NDUP+1
          write(ifle,'(a,I10)') 'ERR: BC:',nb
          write(ifle,'(a,I10)') 'ERR: Duplicating Vertex-4 :',IDCB
        else
          LMAT(IDCB)=IMAT
        endif
!      
  510   CONTINUE
!
      ELSEIF(ISP.GT.0.and.(kd==kdprdc.and.idis(nb)==0)) THEN
!---------------------------------------------------------------
! --- 1) periodic; 2) interface(S&L); 3) inner(buffer); etc  BC:
!---------------------------------------------------------------
!        if(kd==kdprdc.and.idis(nb)==0) then
          nbp=lbface(1,ISP)
          if(nbp+nb.ne.0) then
            write(ifle,*) 
     &      'ERR: IS & ISP is not periodic pair',IS,ISP
            write(ifle,*) 
     &      '### error : BC no.=',nbp,'',nb
            stop 'in metric_CV'
          endif
          do 610 ILE1=1,LEFACE(5,IS)
          IE1=LEFACE(ILE1,IS)
          ICVA1=LVEDGE(1,IE1)
          ICVB1=LVEDGE(2,IE1)
          IVA1=LVRTCV(ICVA1)
          IVB1=LVRTCV(ICVB1)
!---------------------
! --- PERIODIC FACE: 
!---------------------
          DO 670 ILE2=1,LEFACE(5,ISP)
          ICYC=.false.
          IE2=LEFACE(ILE2,ISP)
          ICVA2=LVEDGE(1,IE2)
          ICVB2=LVEDGE(2,IE2)
          IVA2=LVRTCV(ICVA2)
          IVB2=LVRTCV(ICVB2)
!
          DO 672 IA=1,listpr(0,IVA1)
          IF(listpr(IA,IVA1).EQ.IVA2) THEN
            DO 674 IB=1,listpr(0,IVB1)
            IF(listpr(IB,IVB1).EQ.IVB2) THEN
              IVAP2=IVA2
              IVBP2=IVB2
              ICVAP2=ICVA2
              ICVBP2=ICVB2
              ICYC=.TRUE.
              GOTO 680
            ENDIF
  674       CONTINUE
          ENDIF
  672     CONTINUE
!
          DO 676 IA=1,listpr(0,IVA1)
          IF(listpr(IA,IVA1).EQ.IVB2) THEN
            DO 678 IB=1,listpr(0,IVB1)
            IF(listpr(IB,IVB1).EQ.IVA2) THEN
              IVAP2=IVB2
              IVBP2=IVA2
              ICVAP2=ICVB2
              ICVBP2=ICVA2
              ICYC=.TRUE.
              GOTO 680
            ENDIF
  678       CONTINUE
          ENDIF
  676     CONTINUE
! --- IF NOT FOUND GOTO NEXT EDGE OF FACE ISP
          GOTO 670
! --- FOUND EDGE OF FACE ISP
 680      CONTINUE
!-----------------------------------------------------------
! --- Search couterpart edge, then deciding couterpart SSF
!-----------------------------------------------------------
          nssfbc=nssfbc+4
          NCVFAC=NCVFAC+4
          NALLCV=NALLCV+4
          if(NCVFAC>medge)  stop 'Increase [medge]  in cntlnam'
          if(NALLCV>mvrtx)  stop 'Increase [mvrtx]  in cntlnam'
          if(nssfbc>mssfbc) stop 'Increase [mssfbc] in cntlnam'
!
          ICFA1=NCVFAC
          ICFB1=NCVFAC-1
          ICFA2=NCVFAC-2
          ICFB2=NCVFAC-3
!
          IDCA1=NALLCV
          IDCB1=NALLCV-1
          IDCA2=NALLCV-2
          IDCB2=NALLCV-3
!
          ISSFA1=ICFA1-nedge
          ISSFB1=ICFB1-nedge
          ISSFA2=ICFA2-nedge
          ISSFB2=ICFB2-nedge
!
          LCYCSF(ISSFA1)=ISSFA2
          LCYCSF(ISSFB1)=ISSFB2
          LCYCSF(ISSFA2)=ISSFA1
          LCYCSF(ISSFB2)=ISSFB1
!
	  IF(ICYC) THEN
            DO 630 I=1,3
            gfc1(I)=gface(I,IS)
            gvA1(I)=cord(I,IVA1)
            gvB1(I)=cord(I,IVB1)
            gec1(I)=XHALF*(cord(I,IVA1)+cord(I,IVB1))
!
            gfc2(I)=gface(I,ISP)
            gvA2(I)=cord(I,IVAP2)
            gvB2(I)=cord(I,IVBP2)
            gec2(I)=XHALF*(cord(I,IVAP2)+cord(I,IVBP2))
!
! --- Cell Face Vector:
            VEDGE1(I)=area(I,IS)
            VEDGE2(I)=area(I,ISP)
!
! --- edge-center TO cell-face-center
            VFE1(I)=gfc1(I)-gec1(I)
            VFE2(I)=gfc2(I)-gec2(I)
!
! --- edge-center TO A & B
            VCEA1(I)=gvA1(I)-gec1(I)
            VCEB1(I)=gvB1(I)-gec1(I)
            VCEA2(I)=gvA2(I)-gec2(I)
            VCEB2(I)=gvB2(I)-gec2(I)
!
! --- CENTER OF SSF A AND B
!
            GSSA1(I)=r1p3*(gfc1(I)+gvA1(I)+gec1(I))
            GSSB1(I)=r1p3*(gfc1(I)+gvB1(I)+gec1(I))
            GSSA2(I)=r1p3*(gfc2(I)+gvA2(I)+gec2(I))
            GSSB2(I)=r1p3*(gfc2(I)+gvB2(I)+gec2(I))
!
 630  CONTINUE
!---------------------
! --- SSF vector: SSF
!---------------------
            CALL AXBEQC(VFE1,VCEA1,SSFA1,SSFSA1)
            CALL AXBEQC(VFE2,VCEA2,SSFA2,SSFSA2)
            DIRC1=SSFA1(1)*VEDGE1(1)
     &           +SSFA1(2)*VEDGE1(2)
     &           +SSFA1(3)*VEDGE1(3)
            DIRC2=SSFA2(1)*VEDGE2(1)
     &           +SSFA2(2)*VEDGE2(2)
     &           +SSFA2(3)*VEDGE2(3)
            IF(DIRC1.LT.ZERO) THEN
              SSFA1(1)=-SSFA1(1)
              SSFA1(2)=-SSFA1(2)
              SSFA1(3)=-SSFA1(3)
            ENDIF
            IF(DIRC2.LT.ZERO) THEN
              SSFA2(1)=-SSFA2(1)
              SSFA2(2)=-SSFA2(2)
              SSFA2(3)=-SSFA2(3)
            ENDIF
!
            CALL AXBEQC(VFE1,VCEB1,SSFB1,SSFSB1)
            CALL AXBEQC(VFE2,VCEB2,SSFB2,SSFSB2)
            DIRC1=SSFB1(1)*VEDGE1(1)
     &           +SSFB1(2)*VEDGE1(2)
     &           +SSFB1(3)*VEDGE1(3)
            DIRC2=SSFB2(1)*VEDGE2(1)
     &           +SSFB2(2)*VEDGE2(2)
     &           +SSFB2(3)*VEDGE2(3)
            IF(DIRC1.LT.ZERO) THEN
              SSFB1(1)=-SSFB1(1)
              SSFB1(2)=-SSFB1(2)
              SSFB1(3)=-SSFB1(3)
            ENDIF
            IF(DIRC2.LT.ZERO) THEN
              SSFB2(1)=-SSFB2(1)
              SSFB2(2)=-SSFB2(2)
              SSFB2(3)=-SSFB2(3)
            ENDIF
!      
            SFAREA(4,ICFA1)=SSFSA1
            SFAREA(1,ICFA1)=SSFA1(1)/SSFSA1
            SFAREA(2,ICFA1)=SSFA1(2)/SSFSA1
            SFAREA(3,ICFA1)=SSFA1(3)/SSFSA1
!
            SFAREA(4,ICFA2)=SSFSA2
            SFAREA(1,ICFA2)=SSFA2(1)/SSFSA2
            SFAREA(2,ICFA2)=SSFA2(2)/SSFSA2
            SFAREA(3,ICFA2)=SSFA2(3)/SSFSA2
!
            SFAREA(4,ICFB1)=SSFSB1
            SFAREA(1,ICFB1)=SSFB1(1)/SSFSB1
            SFAREA(2,ICFB1)=SSFB1(2)/SSFSB1
            SFAREA(3,ICFB1)=SSFB1(3)/SSFSB1
!
            SFAREA(4,ICFB2)=SSFSB2
            SFAREA(1,ICFB2)=SSFB2(1)/SSFSB2
            SFAREA(2,ICFB2)=SSFB2(2)/SSFSB2
            SFAREA(3,ICFB2)=SSFB2(3)/SSFSB2
!
            SFCENT(1,ICFA1)=GSSA1(1)
            SFCENT(2,ICFA1)=GSSA1(2)
            SFCENT(3,ICFA1)=GSSA1(3)
!
            SFCENT(1,ICFA2)=GSSA2(1)
            SFCENT(2,ICFA2)=GSSA2(2)
            SFCENT(3,ICFA2)=GSSA2(3)
!
            SFCENT(1,ICFB1)=GSSB1(1)
            SFCENT(2,ICFB1)=GSSB1(2)
            SFCENT(3,ICFB1)=GSSB1(3)
!
            SFCENT(1,ICFB2)=GSSB2(1)
            SFCENT(2,ICFB2)=GSSB2(2)
            SFCENT(3,ICFB2)=GSSB2(3)
!
            LVEDGE(1,ICFA1)=ICVA1
            LVEDGE(2,ICFA1)=IDCA1
            LVEDGE(1,ICFA2)=ICVAP2
            LVEDGE(2,ICFA2)=IDCA2
!
            LVEDGE(1,ICFB1)=ICVB1
            LVEDGE(2,ICFB1)=IDCB1
            LVEDGE(1,ICFB2)=ICVBP2
            LVEDGE(2,ICFB2)=IDCB2
!
            LBCSSF(ISSFA1)=nb
            listbc(1,ISSFA1)=IE
            listbc(2,ISSFA1)=IS
            listbc(3,ISSFA1)=IVA1
            LBCSSF(ISSFA2)=nbp
            listbc(1,ISSFA2)=IE
            listbc(2,ISSFA2)=ISP
            listbc(3,ISSFA2)=IVAP2
!
            LBCSSF(ISSFB1)=nb
            listbc(1,ISSFB1)=IE
            listbc(2,ISSFB1)=IS
            listbc(3,ISSFB1)=IVB1
            LBCSSF(ISSFB2)=nbp
            listbc(1,ISSFB2)=IE
            listbc(2,ISSFB2)=ISP
            listbc(3,ISSFB2)=IVBP2
!
            if(LMAT(IDCA1).ne.0.and.IMAT.ne.LMAT(IDCA1)) then
              NDUP=NDUP+1
              write(ifle,'(a,I10)') 'ERR: BC:',nb
              write(ifle,'(a,I10)') 'ERR: Duplicating Vertex-5 :',IDCA1
            else
              LMAT(IDCA1)=IMAT
            endif
            if(LMAT(IDCB1).ne.0.and.IMAT.ne.LMAT(IDCB1)) then
              NDUP=NDUP+1
              write(ifle,'(a,I10)') 'ERR: BC:',nb
              write(ifle,'(a,I10)') 'ERR: Duplicating Vertex-6 :',IDCB1
            else
              LMAT(IDCB1)=IMAT
            endif
            if(LMAT(IDCA2).ne.0.and.IMAT.ne.LMAT(IDCA2)) then
              NDUP=NDUP+1
              write(ifle,'(a,I10)') 'ERR: BC:',nb
              write(ifle,'(a,I10)') 'ERR: Duplicating Vertex-7 :',IDCA2
            else
              LMAT(IDCA2)=IMAT
            endif
            if(LMAT(IDCB2).ne.0.and.IMAT.ne.LMAT(IDCB2)) then
              NDUP=NDUP+1
              write(ifle,'(a,I10)') 'ERR: BC:',nb
              write(ifle,'(a,I10)') 'ERR: Duplicating Vertex-8 :',IDCB2
            else
              LMAT(IDCB2)=IMAT
            endif
!
          ENDIF
!      
  670     CONTINUE
  610     CONTINUE
!        elseif(kd.eq.kdintr) then
!        endif
!
      ENDIF
!
  500 CONTINUE
!
      LMAT0=0
      do ICV=NCV+1,NALLCV
        if(LMAT(ICV).eq.0) then
          LMAT0=LMAT0+1
        endif
      enddo
      if(LMAT0.gt.0) then
        write(ifle,'(a)') 'ERR: some cell NOT defined material Number'
        write(ifle,'(a,I4)') ' ### ERR: NOT defined material : ',LMAT0
        stop 'metric_CV-2-BC'
      endif
!
      if(NDUP.gt.0) then
         write(ifle,'(a,I10)') 'ERR: BC:',nb
        write(ifle,'(a,I10)') 'ERR: Duplicating Vertex number: ',NDUP
        stop ': ERR: Duplicating Vertex'
      endif
!

      if(NCVFAC.gt.medge) then
       write(ifle,*) 'ERR: medge= ',medge,' is not large enough'
       write(ifle,*) 'Please make [medge] > or = ',NCVFAC
       stop
      endif
      if(NALLCV.gt.mvrtx) then
       write(ifle,*) '### error : mvrtx= ',mvrtx,' is not large enough'
       write(ifle,*) 'Please make [mvrtx] > or = ',NALLCV
       stop
      endif
      if(nssfbc.gt.mssfbc) then
        write(ifle,'(a,I10,a)') 
     &   'ERR: mssfbc= ',mssfbc,' is not large enough'
        write(ifle,'(a,I10)') 'Please make [mssfbc] > or = ',nssfbc
        stop
      endif
!
! --- 
!
      if(ical_mvmsh>0) then
        NBCSSF=0
        DO 600 IS=1,nface
        nb=lbface(1,IS)
        if(nb==0) cycle 
        do 650 ILE=1,LEFACE(5,IS)
        NBCSSF=NBCSSF+2
  650   enddo
  600   enddo
        if(NBCSSF/=nssfbc) stop 'MOVE MESH ERROR'
      endif
      nssfbc_old=nssfbc
!
      deallocate(lvrtx,ip,nclv)
!
      ierror=0
!
      return
!
      end subroutine metric_CV
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_CV_cell
     & (mcell,mvrtx,mface,medge,mssfbc,ncell,nface,nedge,
     &  NCV,nssfbc,NCVFAC,NALLCV,nssfbc_old,ncelb,
     &  lcface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,LBCSSF,LCYCSF,listbc,listpr,LMAT,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,
     &  ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only          : ifle,ifll,cntlnam
      use module_partitioner,only : lvrtx=>IW1
      use module_partitioner,only : ip   =>WK1
      use module_partitioner,only : nclv =>WK2
      use module_boundary,only    : kdbcnd,kdintr,kdprdc,kdsld,idis,
     &                              kdbuff,kdshutr,kdpors,kdovst
      use module_parameter,  only : NIFACE
      use module_model,only    : ical_mvmsh
!
      implicit none
!
! 1.  Calculate metrics for each face & cell
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mcell,mvrtx,mface,medge,mssfbc,ncelb
      integer,intent(in)     :: ncell,nface,nedge,NCV
      integer,intent(in)     :: lcface(2,mface)
      integer,intent(in)     :: lbface(2,mface)
      integer,intent(in)     :: lacell(  mcell)
      integer,intent(in)     :: LEFACE(5,mface)
      integer,intent(in)     :: listpr(0:3,mvrtx)
      integer,intent(out)    :: listbc(4,mssfbc)
      integer,intent(inout)  :: LVEDGE(2,medge)
      integer,intent(in)     :: LVRTCV(  mvrtx)
      integer,intent(out)    :: LBCSSF( mssfbc)
      integer,intent(out)    :: LCYCSF( mssfbc)
      integer,intent(out)    :: LMAT  (  mvrtx)
!      
      real*8 ,intent(in)     :: cord  (3,mvrtx)
      real*8 ,intent(in)     :: area  (4,mface)
      real*8 ,intent(in)     :: volume(  mcell)
      real*8 ,intent(in)     :: gface (3,mface)
      real*8 ,intent(in)     :: gcell (3,mcell)
      real*8 ,intent(out)    :: SFAREA(4,medge)
      real*8 ,intent(out)    :: SFCENT(3,medge)
      real*8 ,intent(out)    :: CVVOLM(  mvrtx)
      real*8 ,intent(out)    :: CVCENT(3,mvrtx)
      integer,intent(out)    :: ierror,nssfbc,NCVFAC,NALLCV,
     &                          nssfbc_old
!
! --- [local entities]
!
      integer :: idmax,LMAT0,SMLV
      integer :: KMAT(2),IDC(2)
      integer :: i,j,k,l,m,n,kd,kdv,kdt,kdy,kdk,kdp,nbp,ISP
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ILV,ILS,ILE,ILCV,IEE
!
      logical :: ICYC

      real*8 ,parameter :: r1p6=1.d0/6.d0
      real*8 ,parameter :: r1p3=1.d0/3.d0
      real*8 ,parameter :: r1p4=0.25d0
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
      
      NCVFAC=nface
      NALLCV=ncelb
      nssfbc=nface-NIFACE
      
!
      do IS=1,NCVFAC
      SFAREA(1,IS)=area(1,IS)
      SFAREA(2,IS)=area(2,IS)
      SFAREA(3,IS)=area(3,IS)
      SFAREA(4,IS)=area(4,IS)
      enddo
!
      do IS=1,NCVFAC
      SFCENT(1,IS)=gface(1,IS)
      SFCENT(2,IS)=gface(2,IS)
      SFCENT(3,IS)=gface(3,IS)
      enddo
!
      do IC=1,NALLCV
      CVVOLM(IC)=volume(IC)
      LMAT(IC)=lacell(IC)
      CVCENT(1,IC)=gcell(1,IC)
      CVCENT(2,IC)=gcell(2,IC)
      CVCENT(3,IC)=gcell(3,IC)
      enddo
!
      do IS=1,nface
      LVEDGE(1,IS)=lcface(1,IS)
      LVEDGE(2,IS)=lcface(2,IS)
      enddo
!
      do IS=1,nssfbc
        LBCSSF(IS)=lbface(1,IS+nedge)
      enddo
!------------------------------------
! --- initializing array of CVCENT
!------------------------------------
!
!      CVCENT(:,:)=0.D0
!      CVVOLM(:)=0.D0

!      LMAT(:)=0
!      LCYCSF(:,:)=0
!      LBCSSF(:)=0
!      ICYC=.false.
      
!      SFAREA(:,:)=0.d0
!      SFCENT(:,:)=0.d0
!      LVEDGE(:,:)=0   ! BC only

!=====================================================================
! --- <2> Defining SSF on BC
!=====================================================================
!
!      nssfbc=0
!      NCVFAC=nedge
!      NALLCV=NCV
!
      ierror=0
!
      return
!
      end subroutine metric_CV_cell
!

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_3ds
     & (mcell,mvrtx,mface,ncell,nface,
     &  lcface,lvface,lbface,lacell,cord,
     &  area,volume,gface,gcell,lvcell,lfcell,
     &  ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only : ifle
      use module_partitioner,only : lvrtx =>IW1
!      use module_partitioner,only : lfcell=>IW2      
!
! 1. Calculate metrics for each face & cell
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mface,mcell,mvrtx,nface,ncell
      integer,intent(in)     :: lcface(2,mface)
      integer,intent(in)     :: lvface(4,mface)
      integer,intent(in)     :: lbface(2,mface)
      integer,intent(in)     :: lacell(  mcell)
      integer,intent(in)     :: lvcell(8,mcell)
      integer,intent(in)     :: lfcell(7,mcell)
      real*8 ,intent(in)     :: cord  (3,mvrtx)
      real*8 ,intent(out)    :: area  (4,mface)
      real*8 ,intent(out)    :: volume(  mcell)
      real*8 ,intent(out)    :: gface (3,mface)
      real*8 ,intent(out)    :: gcell (3,mcell)
      integer,intent(out)    :: ierror
!
! --- [local entities]
!
!      integer :: lvrtx(4,mface)
!
      integer :: i,j,k,m,n,k1,k2,nf,ie,iemax
      integer :: IC,IS,IV,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE
!
      
      integer,parameter :: itr(2,2)=
     &   reshape( source=(/1,3,2,4/), shape=(/2,2/) )
      real*8 ,parameter :: r1p6=1.d0/6.d0
      real*8 ,parameter :: r1p3=1.d0/3.d0
      real*8 ,parameter :: r1p4=0.25d0
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8  :: r(3,4),r1(3),r2(3),r3(3),r4(3),rg(4),aa(3,4)
      real*8  :: ss1,ss2,sx,sy,sz,vv
      real*8  :: sx1,sy1,sz1,sx2,sy2,sz2,vol,dum1
      integer,parameter :: lvfcel(4,6,4)=reshape( source=
     &  (/1,3,2,0, 2,3,4,0, 3,1,4,0, 4,1,2,0, 0,0,0,0, 0,0,0,0,
     &    1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0, 0,0,0,0,
     &    1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2,0, 4,5,6,0, 0,0,0,0,
     &    1,5,8,4, 2,3,7,6, 1,2,6,5, 3,4,8,7, 1,4,3,2, 5,6,7,8/),
     &  shape=(/4,6,4/) )
      integer :: lvf (4,6,4)
!
      ierror=0
!
! --- initializing array of CVCENT
!
      do 112 k=1,4
      do 112 j=1,6
      do 111 i=1,4
      lvf(i,j,k)=lvfcel(i,j,k)
 111  continue
      if(lvf(4,j,k).lt.1) lvf(4,j,k)=8
 112  continue
!----------------------------------------------------
! --- CELL METRIX CALCULATION
!----------------------------------------------------
!
      allocate(lvrtx(4,mface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating lvrtx(:,:) in metric_3ds'
!      allocate(lfcell(7,mcell),stat=ierr1) 
!      if(ierr1.ne.0)stop'stop at allocating lfcell(:,:) in metric_3ds'
!
!-< 1. Calculate area (face vector) & center of face >-
!
      do 100 IS=1,nface
!
      do 101 ILV=1,4
      lvrtx(ILV,IS)=lvface(ILV,IS)
      rg(ILV)=0.d0
  101 continue
! --- if/not 3p-face?
      if(lvrtx(4,IS).lt.1) then
        lvrtx(4,IS)=lvrtx(1,IS)
      endif
      do 102 ILV=1,4
      do 102 i=1,3
      r(i,ILV)=cord(i,lvrtx(ILV,IS))
  102 continue
!
!--< 1.1 cell-face area (multiplied by 2) >--
!
! --- assert: Outside-ward for IC1 cell
!
!
      do 110 i=1,3
      r1(i)=r(i,3)-r(i,1)
      r2(i)=r(i,4)-r(i,2)
  110 continue
      sx=r1(2)*r2(3)-r1(3)*r2(2)
      sy=r1(3)*r2(1)-r1(1)*r2(3)
      sz=r1(1)*r2(2)-r1(2)*r2(1)
      area(1,IS)=sx
      area(2,IS)=sy
      area(3,IS)=sz
      area(4,IS)=sqrt(sx*sx+sy*sy+sz*sz)
      if( area(4,IS).le.0.d0 ) goto 9001
!
!--< 1.2 cell-face center >--
!
      do 120 k1=1,2
      k2=3-k1
      do 121 i=1,3
      r1(i)=r(i,itr(1,k2))-r(i,itr(k1,k1))
      r2(i)=r(i,itr(2,k2))-r(i,itr(k1,k1))
      dum1 =r(i,itr(1,k2))+r(i,itr(2,k2))
      r3(i)=r(i,itr(k1,k1))+dum1
      r4(i)=r(i,itr(k2,k1))+dum1
  121 continue
      sx1=r1(2)*r2(3)-r1(3)*r2(2)
      sy1=r1(3)*r2(1)-r1(1)*r2(3)
      sz1=r1(1)*r2(2)-r1(2)*r2(1)
      ss1=sign(1.d0,sx*sx1+sy*sy1+sz*sz1)
     &   *sqrt(sx1*sx1+sy1*sy1+sz1*sz1)
      sx2=sx-sx1
      sy2=sy-sy1
      sz2=sz-sz1
      ss2=sign(1.d0,sx*sx2+sy*sy2+sz*sz2)
     &   *sqrt(sx2*sx2+sy2*sy2+sz2*sz2)
      do 122 i=1,3
      rg(i)=rg(i)+ss1*r3(i)+ss2*r4(i)
  122 continue
      rg(4)=rg(4)+ss1+ss2
  120 continue
      if( rg(4).le.0.d0 ) goto 9001
      dum1=r1p3/rg(4)
      do 123 i=1,3
! --- Face center's coordinates
      gface(i,IS)=rg(i)*dum1
  123 continue
!
  100 continue
!
!
!-< 2. Calculate volume/center of cell >-
! --- set local lfcell
!
!      lfcell=0
!      do 290 IS=1,nface
!      do 290 i=1,2
!      IC=lcface(i,IS)
!      ILS=lfcell(7,IC)+1
!      lfcell(ILS,IC)=IS
!      lfcell(7,IC)=ILS
!  290 continue
!
      do 200 IC=1,ncell
!
!--< 2.1 cell volume >--
!
      vol=0.d0
      do 201 j=1,4
      do 201 i=1,3
      aa(i,j)=0.d0
  201 continue
!
      do 210 ILS=1,lfcell(7,IC)
      IS=lfcell(ILS,IC)
      do 211 i=1,3
      r1(i)=cord(i,lvrtx(1,IS))
     &     +cord(i,lvrtx(2,IS))
     &     +cord(i,lvrtx(3,IS))
      if(lvrtx(4,IS).gt.0) then
         r1(i)=r1(i)+cord(i,lvrtx(4,IS))
      end if
  211 continue
      vv=r1(1)*area(1,IS)+r1(2)*area(2,IS)+r1(3)*area(3,IS)
!
! --- Distinguish IC1 and IC2 for not getting negative volume.
!
      if(IC.eq.lcface(1,IS)) then
! --- IC1 cell
        vol=vol+vv
        do 212 i=1,3
        do 213 j=1,3
        aa(i,j)=aa(i,j)+gface(i,IS)*area(j,IS)
  213   continue
        aa(i,4)=aa(i,4)+gface(i,IS)*vv
  212   continue
      else
! --- IC2 cell
        vol=vol-vv
        do 214 i=1,3
        do 215 j=1,3
        aa(i,j)=aa(i,j)-gface(i,IS)*area(j,IS)
  215   continue
        aa(i,4)=aa(i,4)-gface(i,IS)*vv
  214   continue
      endif
  210 continue
!
      vol=r1p4*vol
      volume(IC)=r1p6*vol
      if(vol.le.SML) then
         write(ifle,*) 'ERR: Nagitive volume'
         write(*,1100) IC,(lvcell(i,IC),i=1,8)
         do i=1,8
         if(lvcell(i,IC).gt.0) then
            write(*,1200) i,lvcell(i,IC),(cord(j,lvcell(i,IC)),j=1,3)
         endif
         enddo
         do ILS=1,lfcell(7,IC)
         IS=lfcell(ILS,IC)
         write(ifle,1000) ILS,(lvface(i,IS),i=1,4)
         write(ifle,1000) ILS,(lvcell(lvf(i,ILS,4),IC),i=1,4)
         enddo
 1000    format(4x,'MSG: ',I8,' face: ',4I8)
 1100    format(4x,'MSG: IC= ',I8,' IV= ',8I8)
 1200    format(4x,'MSG: cord: ',2I8,3F15.6)
         write(ifle,*) '### error : invalid coordinate'
         write(ifle,*) 'volume of cell is not positive'
         write(ifle,*) 'cell no. =',IC,volume(IC),lacell(IC)
         stop 9002
      endif
!
!      scale(IC)=volume(IC)**(1.d0/3.d0)
!
!--< 2.2 center >--
!
      do 220 i=1,3
      aa(i,i)=aa(i,i)+vol
      aa(i,4)=r1p4*aa(i,4)
  220 continue
!
      do 221 i=1,3
      if(aa(i,i).le.0.d0) then
       write(ifle,*) '### error : invalid coordinate'
       write(ifle,*) 'volume of cell is not positive'
       write(ifle,*) 'cell no. =',IC,volume(IC),lacell(IC)
       stop 90021
      endif
      do 222 j=i+1,4
      aa(i,j)=aa(i,j)/aa(i,i)
  222 continue
      do 223 k=i+1,3
      do 223 j=i+1,4
      aa(k,j)=aa(k,j)-aa(i,j)*aa(k,i)
  223 continue
  221 continue
!
      gcell(3,IC)=aa(3,4)
      gcell(2,IC)=aa(2,4)- aa(2,3)*gcell(3,IC)
      gcell(1,IC)=aa(1,4)-(aa(1,2)*gcell(2,IC)+aa(1,3)*gcell(3,IC))
!
  200 continue
!
!
!-< 3. Normalizing face vector to unit vector>-
!
      do 300 IS=1,nface
      dum1=1.d0/area(4,IS)
      do 301 i=1,3
      area(i,IS)=area(i,IS)*dum1
  301 continue
      area(4,IS)=XHALF*area(4,IS)
  300 continue
!
! --- CHECK VOLUME/=0.D0
!
!
      deallocate(lvrtx)   !,lfcell)
!
      call debug
!
!
      return
!
 9001 continue
      ie=3+min(1,abs(lvrtx(4,IC)-lvrtx(1,IC)))
      write(ifle,*) '### error : invalid coordinate'
      write(ifle,*) 'area of face is not positive'
      write(ifle,*) 'face no. =',IC
      write(ifle,*) 'vertices of face =',(lvrtx(i,IC),i=1,ie)
      goto 9999
!
 9002 continue
      write(ifle,*) '### error : invalid coordinate'
      write(ifle,*) 'volume of cell is not positive'
      write(ifle,*) 'cell no. =',IC,volume(IC),lacell(IC)
      goto 9999
!
 9999 continue
      write(ifle,*) '(metric_3ds)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(3).eq.0 ) return
      call printw('area/metric_3ds',area,4,nface)
      call printx('volume/metric_3ds',volume,1,ncell)
      call printw('gface/metric_3ds',gface,3,nface)
      call printw('gcell/metric_3ds',gcell,3,ncell)
      end subroutine debug
!
      end subroutine metric_3ds
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE AXBEQC(A,B,C,D)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------
!     C=0.5(A X B)
!---------------------------------------
      real*8 ,intent(in)   :: A(3),B(3)
      real*8 ,intent(OUT)  :: C(3),D
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
!
!      D=0.D0
!      C=0.D0
      C(3)=XHALF*(A(1)*B(2)-A(2)*B(1))
      C(2)=XHALF*(A(3)*B(1)-A(1)*B(3))
      C(1)=XHALF*(A(2)*B(3)-A(3)*B(2))
      D=DSQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))
!      C(3)=C(3)/D
!      C(2)=C(2)/D
!      C(1)=C(1)/D
!
      RETURN
      END subroutine AXBEQC
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine Merge_SSF(mssfbc,medge,nedge,mvrtx,
     &           NCV,nssfbc,nssfbc_old,NCVFAC,NALLCV,
     &           LMAT,LBCSSF,LCYCSF,LVEDGE,listbc,SFCENT,SFAREA)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Merge SSF on BC
!
! --- [module arguments]
!
      use module_io,only          : ifle,ifll,angle_ssf
      use module_boundary,only    : nbcnd,kdbcnd
      use module_partitioner,only : WK1,WK2,WK3, WK4,WK5, WK6,
     &                              SFAREA_L=>BC_IV3D,
     &                              SFCENT_L=>BC_MV3D
      
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mssfbc,medge,nedge,NCV,mvrtx
      integer,intent(inout) :: nssfbc,NCVFAC,NALLCV,nssfbc_old
      integer,intent(inout) :: LBCSSF(  mssfbc)
      integer,intent(inout) :: LCYCSF(  mssfbc)
      integer,intent(inout) :: LVEDGE(2,medge)
      integer,intent(inout) :: LMAT  (  mvrtx)
      integer,intent(inout) :: listbc(4,mssfbc)
      real*8 ,intent(inout) :: SFCENT(3,medge)
      real*8 ,intent(inout) :: SFAREA(4,medge)
!
! --- [local entities]
!
      integer :: MMAT(-100:100),MAT_NO(200),ILOG,i
      integer :: MAT_BCIDX(0:200)
      integer :: NSF=1,MXNSF=0,MXIDCL=0
      integer,allocatable :: nb(:),mark(:),IBFLO(:),ILMARK(:)
      integer,allocatable :: DCNSSF(:)
      integer :: ISSF,IBF,ICF,ICV,IDCSSF,NDCSSF,IDC,IDCS,IDCE
      integer :: J,NB1,NB2,NDCL,IDCLL,IDCL,NDCL2,IMARK,IBFNEW,IBFNEWP
      integer :: NBFNEW,JJ,NN1,NN2
!      integer :: IBFP,ierr1=0,nssf=24
!      real*8  :: SFFC(24,4),SFCT(24,3),VECT
      integer :: IBFP,ierr1=0,nssf=40       !onishi
      real*8  :: SFFC(40,4),SFCT(40,3),VECT

      real*8  :: SFFCT(4)
      real*8  :: SFCTT(3)
      real*8  :: SSF_agl,SSF
!-----------------------------------------
! --- 
!-----------------------------------------
      allocate(WK1(NCV),WK2(nssfbc),WK3(NCV),stat=ierr1)
      if(ierr1.ne.0) stop 'Error-1: allocation in Merge_SSF'
!
      if(angle_ssf/=20.d0) then
        write(ifll,'(a,E14.4)') 'MSG: angle_ssf(deg.): ',angle_ssf
      endif
      SSF_agl=cos(angle_ssf*3.1415926d0/180.d0)
!
! --- List SSF on BC
!
      WK1(:)=0
      DO 100 IBF=1,nssfbc
      ICF=NEDGE+IBF
      ICV=LVEDGE(1,ICF)
      WK1(ICV)=WK1(ICV)+1
 100  CONTINUE
!
      IDCSSF=0
      WK2(:)=0
      WK3(:)=0
      DO 200 ICV=1,NCV
      IF(WK1(ICV).GT.0) THEN
        IDCSSF=IDCSSF+1
        WK2(IDCSSF)=ICV
        WK3(ICV)=IDCSSF
      ENDIF
 200  CONTINUE
!
      NDCSSF=IDCSSF
      ALLOCATE(WK4(0:NDCSSF),stat=ierr1)
      if(ierr1.ne.0) stop 'Error-3: allocation in Merge_SSF'
      WK4(:)=0
      DO 400 IDCSSF=1,NDCSSF
      ICV=WK2(IDCSSF)
      WK4(IDCSSF)=WK1(ICV)
 400  CONTINUE
!
      do IDCSSF=1,NDCSSF
      nssf=max(WK4(IDCSSF),nssf)
      enddo
!
      DO IDCSSF=1,NDCSSF
      WK4(IDCSSF)=WK4(IDCSSF)+WK4(IDCSSF-1)
      ENDDO
! --- 
      ALLOCATE(WK5(WK4(NDCSSF)),stat=ierr1)
      if(ierr1.ne.0) stop 'Error-4: allocation in Merge_SSF'
!
      DO 410 IBF=1,nssfbc
      ICF=NEDGE+IBF
      ICV=LVEDGE(1,ICF)
      IDCSSF=WK3(ICV)
      WK4(IDCSSF-1)=WK4(IDCSSF-1)+1
      WK5(WK4(IDCSSF-1))=IBF
 410  CONTINUE
!
! --- WK3
!
      deallocate(WK3)
!
      WK4(:)=0
      DO 430 IDCSSF=1,NDCSSF
      ICV=WK2(IDCSSF)
      WK4(IDCSSF)=WK1(ICV)
 430  CONTINUE
!
! --- WK1
!
      deallocate(WK1)
      ALLOCATE(WK1(nssfbc),DCNSSF(NDCSSF),stat=ierr1)
      DCNSSF(:)=0
      if(ierr1.ne.0) stop 'Error-5: allocation in Merge_SSF'
!
      DO IDCSSF=1,NDCSSF
      WK4(IDCSSF)=WK4(IDCSSF)+WK4(IDCSSF-1)
      ENDDO
!
 445  continue
      allocate(nb(nssf),mark(nssf),IBFLO(nssf),stat=ierr1)
      WK1(:)=0
!---------------
! --- Merge ssf
!---------------
      NN1=0
      NN2=0
      do 647 IDCSSF=1,NDCSSF
      IDCS=WK4(IDCSSF-1)+1
      IDCE=WK4(IDCSSF)
      IDCL=0
      SFCT(:,:)=0.D0
      SFFC(:,:)=0.D0
      mark(:)=0
      nb(:)=0
      IBFLO(:)=0
! --- -------------------------------------
      do IDC=IDCS,IDCE
        IDCL=IDCL+1
        mark(IDCL)=IDCL
        IBF=WK5(IDC)

!      if(LBCSSF(IBF).eq.1 ) NN1=NN1+1
!      if(LBCSSF(IBF).eq.-1) NN2=NN2+1

        nb(IDCL)=LBCSSF(IBF)
        ICF=IBF+NEDGE
        IBFLO(IDCL)=IBF
        DO  J=1,4
          SFFC(IDCL,J)=SFAREA(J,ICF)
        enddo
        DO  J=1,3
          SFCT(IDCL,J)=SFCENT(J,ICF)
        enddo
      enddo
      MXIDCL=MAX(MXIDCL,IDCL)
! ------------- ------------- -------------
      NSF=1
      IBF=IBFLO(NSF)
      NDCL=1
 645  continue
      WK1(IBF)=-NSF
      mark(NDCL)=-NSF
      NB1=nb(NDCL)
!
      DO 640 IDCLL=NDCL+1,IDCL
      IBF=IBFLO(IDCLL)
      NB2=nb(IDCLL)
      IF(mark(IDCLL).NE.IDCLL) goto 640
      VECT=SFFC(NDCL,1)*SFFC(IDCLL,1)
     &    +SFFC(NDCL,2)*SFFC(IDCLL,2)
     &    +SFFC(NDCL,3)*SFFC(IDCLL,3)
      IF(VECT.GT.SSF_agl.AND.NB1.EQ.NB2) THEN
        mark(IDCLL)=-NSF
        WK1(IBF)=-NSF
      endIF
 640  CONTINUE
!
      DO 650 IDCLL=NDCL+1,IDCL
        IF(mark(IDCLL).eq.IDCLL) then
          NDCL2=mark(IDCLL)
          GOTO 666
        endif
 650  CONTINUE
!
      GOTO 648
 666  CONTINUE
      NDCL=NDCL2
      IBF=IBFLO(NDCL)
      NSF=NSF+1
      goto 645
 648  CONTINUE
      DCNSSF(IDCSSF)=NSF
      MXNSF=MAX(MXNSF,NSF)
      DO 316 IDCLL=1,IDCL
      ILOG=0
      do J=1,NSF
      IF(mark(IDCLL).EQ.-J)  then
        ILOG=1
      endif
      enddo
      if(ILOG.eq.0) then
        WRITE(ifle,'(a)') ' ###ERR: NOT SUPPORT THIS KIND MESH'
        STOP
      endif
 316  CONTINUE
 647  CONTINUE
!
      allocate(SFAREA_L(nssfbc,4),SFCENT_L(nssfbc,3),stat=ierr1)
      if(ierr1.ne.0) stop 'Error-6: allocation in Merge_SSF'
      allocate(WK6(nssfbc),stat=ierr1)
      if(ierr1.ne.0) stop 'Error-7: allocation in Merge_SSF'
      allocate(WK3(0:nssfbc),stat=ierr1)
      if(ierr1.ne.0) stop 'Error-8: allocation in Merge_SSF'
      allocate(ILMARK(MXNSF),stat=ierr1 )
      if(ierr1.ne.0) stop 'Error-9: allocation in Merge_SSF'
!
! --- 
!
      IBFNEW=0
      DO 500 IDCSSF=1,NDCSSF
      ICV=WK2(IDCSSF)
      IDCS=WK4(IDCSSF-1)+1
      IDCE=WK4(IDCSSF)
      ILMARK(:)=0
      NSF=DCNSSF(IDCSSF)
      do 555 JJ=1,NSF
      SFFCT=0.D0
      SFCTT=0.D0
      do 551 IDC=IDCS,IDCE
      IBF=WK5(IDC)
      IMARK=WK1(IBF)
      ICF=IBF+NEDGE
      IF(IMARK.EQ.-JJ) THEN
        if(ILMARK(JJ).eq.0) IBFNEW=IBFNEW+1
        ILMARK(JJ)=ILMARK(JJ)+1
        WK3(IBF)=IBFNEW
        DO J=1,3
          SFFCT(J)=SFFCT(J)+SFAREA(J,ICF)*SFAREA(4,ICF)
        ENDDO
        SFFCT(4)=SFFCT(4)+SFAREA(4,ICF)
        DO J=1,3
          SFCTT(J)=SFCTT(J)+SFCENT(J,ICF)*SFAREA(4,ICF)
        ENDDO
      endif
 551  continue
      IF(ILMARK(JJ).GT.0) THEN
        SSF=SFFCT(4)
        DO J=1,3
          SFCENT_L(IBFNEW,J)=SFCTT(J)/SSF
        ENDDO
        SSF=SFFCT(1)*SFFCT(1)
     &     +SFFCT(2)*SFFCT(2)
     &     +SFFCT(3)*SFFCT(3)
        SSF=SQRT(SSF)
        SFAREA_L(IBFNEW,4)=SSF
        DO J=1,3
          SFAREA_L(IBFNEW,J)=SFFCT(J)/SSF
        ENDDO
        WK6(IBFNEW)=ICV
      ENDIF
!
      do 554 IDC=IDCS,IDCE
      IBF=WK5(IDC)
      IMARK=WK1(IBF)
      if(IMARK.eq.0) then
        write(ifle,'(a)')
     &    'EER: Contact to your Supportor'
        stop 'in Merge_SSF -3- '
      ENDIF
 554  continue
 555  continue
!
 500  CONTINUE
      NBFNEW=IBFNEW
!
! --- WK2,WK4,WK5,WK1
!
      deallocate(WK2,WK4,WK5,WK1)
!
      allocate(WK2(NBFNEW),WK4(NBFNEW),stat=ierr1)
      if(ierr1.ne.0) stop 'Error-10: allocation in Merge_SSF'
!
      WK3(0)=0
      WK2=0
      WK4=0
      do 700 IBF=1,nssfbc
      IBFNEW=WK3(IBF)
      nb1=LBCSSF(IBF)
      IBFP=LCYCSF(IBF)
      IBFNEWP=WK3(IBFP)
      if(WK2(IBFNEW).eq.0) then
        WK2(IBFNEW)=nb1
      else
        if(nb1.ne.WK2(IBFNEW)) then
          write(*,*) ' ### EER: Same IBFNEW NO.,diff nb in Merge_SSF'
          stop 'in Merge_SSF -1- '
        endif
      endif
      if(WK4(IBFNEW).eq.0) then
        WK4(IBFNEW)=IBFNEWP
      else
        if(IBFNEWP.ne.WK4(IBFNEW)) then
          write(*,*) 
     &    ' ### EER: Contact to your Supportor'
          stop 'in Merge_SSF -2- '
        endif
      endif
 700  continue
!
      listbc(4,1:nssfbc_old)=WK3(1:nssfbc_old)
      deallocate(WK3)
!
      SFAREA(:,nedge+1:)=0.d0
      SFCENT(:,nedge+1:)=0.d0
      LVEDGE(:,nedge+1:)=0
      LMAT(nedge+1:)=0
      LBCSSF(:)=0
      LCYCSF(:)=0
!
      NCVFAC=nedge
      NALLCV=NCV
      do 660 IBFNEW=1,NBFNEW
      NCVFAC=NCVFAC+1
      NALLCV=NALLCV+1
      do J=1,3
      SFAREA(J,NCVFAC)=SFAREA_L(IBFNEW,J)
      SFCENT(J,NCVFAC)=SFCENT_L(IBFNEW,J)
      enddo
      
      SFAREA(4,NCVFAC)=SFAREA_L(IBFNEW,4)
      LVEDGE(1,NCVFAC)=WK6(IBFNEW)
      LVEDGE(2,NCVFAC)=NALLCV
      LMAT(NALLCV)=LMAT(WK6(IBFNEW))
      LBCSSF(IBFNEW)=WK2(IBFNEW)
      LCYCSF(IBFNEW)=WK4(IBFNEW)
 660  continue
      nssfbc=NBFNEW
!
      deallocate(WK2,WK4,WK6,SFAREA_L,SFCENT_L,DCNSSF,ILMARK)
      deallocate(nb,mark,IBFLO)
!
      return
!
      end subroutine Merge_SSF
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
