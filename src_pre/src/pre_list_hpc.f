!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!      subroutine list_face_hpc
!      subroutine list_edge_hpc
!      subroutine list_facebc_hpc
!      subroutine list_facechk_hpc
!      subroutine list_facein_hpc
!      subroutine list_facenn_hpc
!      subroutine list_facepr_hpc
!      subroutine list_fluid_hpc
!      subroutine list_solid_hpc
!      subroutine list_edgvtx_hpc
!      subroutine list_BC_hpc
!      subroutine bc_check_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_face_hpc
     & (mcell,mface,mvrtx,nvrtx_h,ncell_h,nface_h,ncelb_h,NCV_h,
     &  ncelb_local,
     & lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [ arguments]
!
      use module_io,only : ifle
      use module_partitioner,only : nclv=>WK1
      use module_partitioner,only : ityp=>WK2
      use module_partitioner,only : ip  =>WK3
      use module_partitioner,only : lvrtx=>IW1
!
! --- 1. Make up list of connectivity for cell & face
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mcell,mface,mvrtx,nvrtx_h,ncell_h
      integer,intent(in)  :: nface_h,ncelb_h,NCV_h
      integer,intent(out) :: ncelb_local
      integer,intent(in)  :: lvcell(8,mcell)
      integer,intent(out) :: lfcell(7,mcell)
      integer,intent(out) :: lvface(4,mface)
      integer,intent(out) :: lcface(2,mface)

      integer,intent(out) :: LVRTCV(  mvrtx)
      integer,intent(out) :: LCVFAC(4,mface)
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer           :: nface,ncelb,NCV
      integer,parameter :: nv2typ(8)=(/0,0,0,1,2,3,0,4/)
      integer,parameter :: nfacel(4)=(/4,5,5,6/)
! --- define local face vertex order
      integer,parameter :: lvfcel(4,6,4)=reshape( source=
     &  (/1,3,2,0, 2,3,4,0, 3,1,4,0, 4,1,2,0, 0,0,0,0, 0,0,0,0,
     &    1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0, 0,0,0,0,
     &    1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2,0, 4,5,6,0, 0,0,0,0,
     &    1,5,8,4, 2,3,7,6, 1,2,6,5, 3,4,8,7, 1,4,3,2, 5,6,7,8/),
     &  shape=(/4,6,4/))
!
!       nv2typ : type of cell according to face no.
!              : =1;tetrahedron, =2;pyramid, =3;prism, =4;hexahedron
!       nfacel : no. of faces in cell
!       lvfcel : dummy vertex no. in cell face
!
!      integer,allocatable :: lvrtx(:,:)
!      integer :: nclv(8*ncell_h)
!      integer :: ityp(  ncell_h)
!      integer :: ip  (0:nvrtx_h)
!
      integer :: lvf (4,6,4)
      integer :: lvrt(4),IS,IVV
      integer :: i,j,k,l,m,n
      integer :: j1,j2,n1,n2,m1,m2,lf1,ie,ierr1=0,ierr2
      integer :: kclv,nf,icmax,match,nch,nfacx,nclbx
      integer :: IC,IV,ILV,ILS,NF1,NF2,ILS1,ILS2,IC1,IC2,ICTP
      logical :: INNER
!
      allocate(nclv(8*ncell_h),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating nclv(:) in list_face_hpc'
      allocate(ityp(ncell_h),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ityp(:) in list_face_hpc'
      allocate(ip(0:nvrtx_h),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ip(:) in list_face_hpc'
!
! --- 
!
      NCV=NCV_h
      DO 50 IV=1,NCV
      LVRTCV(IV)=IV
  50  CONTINUE
!
! --- < 1. Set up some arrays >-
!
! --- < 1.1 make up "lvf" >--
!
      ierror=0
      ierr1=0
      ierr2=0
!
      do 100 k=1,4
      do 100 j=1,6
      do 101 i=1,4
      lvf(i,j,k)=lvfcel(i,j,k)
  101 continue
      if( lvf(4,j,k).lt.1 ) lvf(4,j,k)=8
  100 continue
!
! --- < 1.2 clear array >--
!
      do 110 IC=1,mcell
      do 110 i=1,7
      lfcell(i,IC)=0
  110 continue
!
      do 111 IS=1,mface
      lcface(1,IS)=0
      lcface(2,IS)=0
      do 111 i=1,4
      lvface(i,IS)=0
      LCVFAC(i,IS)=0
  111 continue
!
! --- < 1.3 set type & no. of faces >--
!
      do 120 IC=1,ncell_h
      kclv=0
      do 121 ILV=1,8
      if(lvcell(ILV,IC).gt.0) kclv=ILV
  121 continue
      if( kclv.lt.1 ) then
        goto 9001
      endif
      ityp(IC)=nv2typ(kclv)
      if( ityp(IC).lt.1 ) then
        goto 9001
      endif
      lfcell(7,IC)=nfacel(ityp(IC))
  120 continue
!
! --- < 2. Set up cell-face list >-
!
! --- < 2.1 make up list for cells connected with each vertex >--
!
      call list_onvrtx(8,mcell,nvrtx_h,ncell_h,
     & lvcell,ip,nclv,icmax)
!
      allocate(lvrtx(6,4*icmax),stat=ierr1)
      if(ierr1.ne.0) goto 9002
!
! --- < 2.2 Search internal faces >-
! --- icmax is maximun no. of cells which connect one vertex.
!
      nface=0
!
      do 200 IV=1,nvrtx_h
!
! --- / set vertex of faces in cell /
!
      nf=0
      do 210 IVV=ip(IV-1)+1,ip(IV)
      IC=nclv(IVV)
      ICTP=ityp(IC)
      do 211 ILS=1,lfcell(7,IC)
      if(lfcell(ILS,IC).gt.0) goto 211
      do 213 ILV=1,4
      lvrt(ILV)=lvcell(lvf(ILV,ILS,ICTP),IC)
  213 continue
!
      INNER=lvrt(1).EQ.IV
     &  .OR.lvrt(2).EQ.IV
     &  .OR.lvrt(3).EQ.IV
     &  .OR.lvrt(4).EQ.IV
      IF(INNER) THEN
        nf=nf+1
        do 212 ILV=1,4
        lvrtx(ILV,nf)=lvcell(lvf(ILV,ILS,ICTP),IC)
  212   continue
        lvrtx(5,nf)=ILS
        lvrtx(6,nf)=IC
      ENDIF
  211 continue
  210 continue
!
! --- matching procedure
!
      do 220 nf1=1,nf-1
      ILS1=lvrtx(5,nf1)
      IC1=lvrtx(6,nf1)
      if(lfcell(ILS1,IC1).gt.0) goto 220
      do 221 nf2=nf1+1,nf
      ILS2=lvrtx(5,nf2)
      IC2=lvrtx(6,nf2)
      if(lfcell(ILS2,IC2).gt.0) goto 221
      call list_fmatch(match,lvrtx(1,nf1),lvrtx(1,nf2),-1)
      if(match.lt.1) goto 221
      nface=nface+1
      nfacx=min(nface,mface)
      lfcell(ILS1,IC1)=nface
      lfcell(ILS2,IC2)=nface
      lcface(1,nfacx)=IC1
      lcface(2,nfacx)=IC2
! --- Set face-vertex-number-order in IC1's outside-ward order
      do 222 i=1,4
      lvface(i,nfacx)=lvrtx(i,nf1)
  222 continue
      goto 220
  221 continue
  220 continue
!
! --- / procedure only for check /
!
      do 230 nf1=1,nf-1
      ILS1=lvrtx(5,nf1)
      IC1=lvrtx(6,nf1)
      lf1=lfcell(ILS1,IC1)
      if(lf1.lt.1) lf1=-1
      do 231 nf2=nf1+1,nf
      ILS2=lvrtx(5,nf2)
      IC2=lvrtx(6,nf2)
      if(lfcell(ILS2,IC2).ne.lf1) then
        call list_fmatch(match,lvrtx(1,nf1),lvrtx(1,nf2),0)
        if(match.gt.0) then
           goto 9003
        ENDIF
      else
        nch=IC2
      endif
  231 continue
  230 continue
!
  200 continue
!
!-< 3. Search boundary faces >-
!
      ncelb=ncell_h
      do 300 IC=1,ncell_h
      k=ityp(IC)
      do 301 j=1,lfcell(7,IC)
      if(lfcell(j,IC).lt.1) then
        ncelb=ncelb+1
        nface=nface+1
        nclbx=min(ncelb,mcell)
        nfacx=min(nface,mface)
        lfcell(j,IC)=nface
! --- 'ncell~ncelb' is 2D BC cell no., CL1 is 3d-cell; CL2 is 2d-BC-cell
        lcface(1,nfacx)=IC
        lcface(2,nfacx)=ncelb
        do 302 i=1,4
        lvface(i,nfacx)=lvcell(lvf(i,j,k),IC)
  302   continue
      endif
  301 continue
  300 continue
!
      do 410 IS=1,nface
      do 410 i=1,4
      LCVFAC(i,IS)=lvface(i,IS)
  410 continue
!
      call list_facerr(ifle,mcell,mface,nface,ncelb,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if(nface_h.ne.nface) then
        write(ifle,*) '### error : nface not equal to nface_h'
        goto 9999
      else
        write(ifle,*) '    ###    NFACE_H= ',
     &               NFACE_H,' NFACE_LOCAL= ',NFACE
      endif
!
      write(ifle,*) '    ###    NCELB_H= ',
     &                NCELB_H,' NCELB_LOCAL= ',NCELB
      ncelb_local=ncelb
!
! --- 
!
      deallocate(lvrtx,nclv,ityp,ip)
!

      call debug
!
      return
!
 9001 continue
      write(ifle,*) '### error 1 : data error'
      write(ifle,*) 'no. of vertices =',kclv,ncell_h
      write(ifle,*) 'it must be 4,5,6 or 8'
      write(ifle,*) 'cell no. =',IC
      goto 9999
 9002 continue
      write(ifle,*) '### error 2 : allocation failed'
      goto 9999
 9003 continue
      write(ifle,*) '### error 3 : data error'
      call list_fmatch(match,lvrtx(1,nf1),lvrtx(1,nf2),-1)
      if( match.gt.0 ) then
        write(ifle,*) 'the three cells are connected in the same face'
        write(ifle,*) 'cell no. =',IC1,',',IC2,' and',nch
      else
        write(ifle,*) 'the two cells share their vetices ',
     &              'in incorrect order'
        write(ifle,*) 'cell no. =',IC1,' and',IC2
      endif
      ie=3+min(1,lvrtx(4,nf1))
      write(ifle,*) 'vertices of the face =',(lvrtx(i,nf1),i=1,ie)
      goto 9999
 9999 continue
      write(ifle,*) '(list_face)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine debug
!=================================================
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printi('lfcell/list_face',lfcell,7,mcell)
      call printi('lvface/list_face',lvface,4,mface)
      call printi('lcface/list_face',lcface,2,mface)
      end subroutine debug
!
      end subroutine list_face_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE LIST_EDGE_hpc
     & (MCELL,MFACE,MEDGE,mvrtx,
     &  NCELL_h,NVRTX_h,NFACE_h,NEDGE_h,NCELB_h,NCV_h,
     &  LVCELL,LFCELL,LVFACE,LCFACE,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,IERROR)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      USE MODULE_IO,ONLY : IFLE
      use module_partitioner,only : nclv=>WK1
      use module_partitioner,only : ityp=>WK2
      use module_partitioner,only : ip  =>WK3
      use module_partitioner,only : lvrtx=>IW1
!
!     1. MAKE UP LIST OF CONNECTIVITY EDGE-TO-FACE AND EDGE-TO-VERTEX
!
      implicit none
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)  :: MCELL,MFACE,MEDGE,mvrtx,NCV_h
      INTEGER,INTENT(IN)  :: NCELL_h,NVRTX_h,NFACE_h,NCELB_h
      INTEGER,INTENT(IN)  :: NEDGE_h
      INTEGER,INTENT(IN)  :: LVCELL(8,MCELL)
      INTEGER,INTENT(IN)  :: LFCELL(7,MCELL)
      INTEGER,INTENT(IN)  :: LVFACE(4,MFACE)
      INTEGER,INTENT(IN)  :: LCFACE(2,MFACE)
      integer,intent(in)  :: LVRTCV(  mvrtx)
      integer,intent(in)  :: LCVFAC(4,mface)
      INTEGER,INTENT(OUT) :: LEFACE(5,MFACE)
      INTEGER,INTENT(OUT) :: LVEDGE(2,MEDGE)
      INTEGER,INTENT(OUT) :: IERROR
!
!
! --- [LOCAL ENTITIES]
!
      INTEGER           :: NEDGE
      INTEGER,PARAMETER :: NF2TYP(4)=(/0,0,1,2/)
      INTEGER,PARAMETER :: NEDGEL(2)=(/3,4/)
! --- DEFINE LOCAL EDGE VERTEX ORDER
      INTEGER,PARAMETER :: LVEDGL(2,4,2)=RESHAPE( SOURCE=
     &            (/1,2, 2,3, 3,1, 0,0,
     &              1,2, 2,3, 3,4, 4,1/),SHAPE=(/2,4,2/))
!
!      INTEGER,ALLOCATABLE :: LVRTX(:,:)
!      INTEGER :: NCLV(4*NFACE_h)
!      INTEGER :: ITYP(NFACE_h)
!      INTEGER :: IP(0:NVRTX_h)
!
      INTEGER :: LVF(2,4,2)
      INTEGER :: LVRT(2)
      INTEGER :: LF1,IE,IERR1,IERR2
      INTEGER :: KCLV,IFMAX,MATCH,NCH,NEDGEX,NCLBX
      INTEGER :: I,J,K,L,M,N
      INTEGER :: IC,IS,IV,IVA,IVB,IC1,IC2
! delete already given attribute -- by onishi
!      INTEGER :: IC,IS,IV,IE,IVA,IVB,IC1,IC2
      INTEGER :: ILS,ILV,ILE,IVV,IEE,ILCV,ICV,ICTP
      INTEGER :: ILE1,ILE2,IS1,IS2,NE1,NE2,NE,IEGVT1,IEGVT2
      LOGICAL :: INNER,IEDGE
!
! ??? nedge has to be created in 'inner and interface', not finished, 
!
!
!-< 1. SET UP SOME ARRAYS >-
      allocate(nclv(4*NFACE_h),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating nclv(:) in LIST_EDGE_hpc'
      allocate(ityp(NFACE_h),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ityp(:) in LIST_EDGE_hpc'
      allocate(ip(0:NCV_h),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ip(:) in LIST_EDGE_hpc'
!
!--< 1.1 MAKE UP "LVF" >--
!
      IERROR=0
!
      DO 100 K=1,2
      DO 100 J=1,4
      DO 101 I=1,2
      LVF(I,J,K)=LVEDGL(I,J,K)
  101 CONTINUE
  100 CONTINUE
!
!--< 1.2 CLEAR ARRAY >--
!
      DO 110 IS=1,MFACE
      DO 110 ILE=1,5
      LEFACE(ILE,IS)=0
  110 CONTINUE
!
      DO 112 IE=1,MEDGE
      LVEDGE(1,IE)=0
      LVEDGE(2,IE)=0
  112 CONTINUE
      
!
!--<  1.3 SET-UP TYPE & LOCAL NO. OF EDGE >--
!
      DO 120 IS=1,NFACE_h
      KCLV=0
      DO 121 ILCV=1,4
      IF(LCVFAC(ILCV,IS).GT.0) KCLV=ILCV
  121 CONTINUE
      IF(KCLV.LT.1) GOTO 9001
      ITYP(IS)=NF2TYP(KCLV)
      IF(ITYP(IS).LT.1) GOTO 9001
      LEFACE(5,IS)=NEDGEL(ITYP(IS))
  120 CONTINUE
!
!
!-< 2. SET UP CELL-FACE LIST >-
!
!--< 2.1 MAKE UP LIST FOR FACES CONNECTED WITH EACH VERTEX >--
!NCV_h
      IFMAX=1
      CALL LIST_ONVRTX(4,MFACE,NVRTX_h,NFACE_h,LCVFAC,IP,NCLV,IFMAX)
!
      ALLOCATE(LVRTX(4,2*IFMAX),STAT=IERR1)
      IF(IERR1.NE.0) GOTO 9002
!
!
!--< 2.2 SEARCH INTERNAL FACES >-
!
      NEDGE=0
!
      DO 200 ICV=1,NVRTX_h
!
! --- SET VERTEX OF EDGES IN CELL 
!
      NE=0
      DO 210 IVV=IP(ICV-1)+1,IP(ICV)
      IS=NCLV(IVV)
      ICTP=ITYP(IS)
      DO 211 ILE=1,LEFACE(5,IS)
      IF(LEFACE(ILE,IS).GT.0) GOTO 211
      DO 213 ILCV=1,2
      LVRT(ILCV)=LCVFAC(LVF(ILCV,ILE,ICTP),IS)
  213 CONTINUE
      INNER=LVRT(1).EQ.ICV.OR.LVRT(2).EQ.ICV
      IF(INNER) THEN
        NE=NE+1
        DO 212 ILCV=1,2
        LVRTX(ILCV,NE)=LCVFAC(LVF(ILCV,ILE,ICTP),IS)
  212   CONTINUE
        LVRTX(3,NE)=ILE
        LVRTX(4,NE)=IS
      ENDIF
  211 CONTINUE
  210 CONTINUE
!
! --- MATCHING PROCEDURE
!
      DO 220 NE1=1,NE-1
      IEDGE=.TRUE.
      ILE1=LVRTX(3,NE1)
      IS1=LVRTX(4,NE1)
      IF(LEFACE(ILE1,IS1).GT.0) GOTO 220
      IF(ICV.EQ.LVRTX(1,NE1)) THEN
	IEGVT1=LVRTX(1,NE1)
	IEGVT2=LVRTX(2,NE1)
      ELSE
	IEGVT1=LVRTX(2,NE1)
	IEGVT2=LVRTX(1,NE1)
      ENDIF
      DO 221 NE2=NE1+1,NE
      ILE2=LVRTX(3,NE2)
      IS2=LVRTX(4,NE2)
      IF(LEFACE(ILE2,IS2).GT.0) GOTO 221
      CALL LIST_EMATCH(MATCH,LVRTX(1,NE1),LVRTX(1,NE2),0)
      IF(MATCH.LT.2) GOTO 221
      IF(IEDGE) THEN
        NEDGE=NEDGE+1
        NEDGEX=MIN(NEDGE,MEDGE)
        LEFACE(ILE1,IS1)=NEDGE
        LVEDGE(1,NEDGEX)=IEGVT1
        LVEDGE(2,NEDGEX)=IEGVT2
	IEDGE=.FALSE.
      ENDIF
      LEFACE(ILE2,IS2)=NEDGE
  221 CONTINUE
  220 CONTINUE
!
  200 CONTINUE
!
      deallocate(lvrtx,nclv,ityp,ip)
!
      CALL DEBUG
!
      write(ifle,*)'    ###    NEDGE_H= ',
     &                         NEDGE_H,' NEDGE_LOCAL= ',NEDGE
      if(NEDGE.ne.NEDGE_h) then
        write(ifle,*) '### error : NEDGE not equal to NEDGE_h'
        goto 9999
      endif
!
      RETURN
!
 9001 CONTINUE
      WRITE(IFLE,*) '### ERROR : DATA ERROR'
      WRITE(IFLE,*) 'NO. OF VERTICES =',KCLV
      WRITE(IFLE,*) 'IT MUST BE 3 OR 4'
      WRITE(IFLE,*) 'FACE NO. =',IS
      GOTO 9999
 9002 CONTINUE
      WRITE(IFLE,*) '### ERROR : ALLOCATION FAILED'
      GOTO 9999
 9999 CONTINUE
      WRITE(IFLE,*) '(LIST_EDGE_hpc)'
      IERROR=1
!
!///////////////////////////////////////////////////////////////////////
      CONTAINS
!=================================================
      SUBROUTINE DEBUG
!=================================================
      USE MODULE_DEBUG,ONLY : IDEBUG
      IF( IDEBUG(2).EQ.0 ) RETURN
      CALL PRINTI('LEFACE/LIST_EDGE_hpc',LEFACE,5,MFACE)
!      CALL PRINTI('LFEDGE/LIST_EDGE_hpc',LFEDGE,MXFACE,MEDGE)
      CALL PRINTI('LVEDGE/LIST_EDGE_hpc',LVEDGE,2,MEDGE)
      END SUBROUTINE DEBUG
!
      END SUBROUTINE LIST_EDGE_hpc
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facebc_hpc
     & (mface,nface_h,mvrtx,nvrtx_h,ncell_h,
     &  lvface,lbface,lcface)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [ arguments]
!
      use module_io,only : ifle
      use module_boundary,only : nbcnd,ivbcnd,lvbcnd,kdbcnd,kdprdc,
     &                           kdsld,kdintr,kdbuff,idis,ivpair
     &                           ,kdtchi,kdshutr,kdpors
      use module_partitioner,only : lbc=>WK1
!      use module_partitioner,only : lbfacx=>WK2
      use module_partitioner,only : lbc1=>WK3
!
! 1. Make up list to link boundary conditions
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mface,nface_h,nvrtx_h,mvrtx,ncell_h
      integer,intent(in)    :: lvface(4,mface)
      integer,intent(in)    :: lcface(2,mface)
      integer,intent(inout) :: lbface(2,mface)

! NOT finished for interface BC
!
! --- [local entities]
!
      integer :: i,m,n,nb,ierr1
      integer :: IC,IS,IV,IE,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC,kd,id
      logical :: logflg=.false.
!
      allocate(lbc(0:nvrtx_h),stat=ierr1) 
      if(ierr1.ne.0) 
     & stop 'stop at allocating lbc(:) in list_facebc_hpc'
      allocate(lbc1(0:nvrtx_h),stat=ierr1)
      if(ierr1/=0) stop 'stop at allocating lbc1(:) in list_facebc'
!      allocate(lbfacx(mface),stat=ierr1) 
!      if(ierr1.ne.0) 
!     & stop 'stop at allocating lbfacx(:) in list_facebc_hpc'
!
! --- 
!
      lbc(1:nvrtx_h)=0
      lbc(0)=1
      lbc1(1:nvrtx_h)=0
      lbc1(0)=-1
      
!
!      do 110 IS=1,nface_h
!      nb=lbface(1,IS)
!      lbfacx(IS)=nb
!  110 continue
!
      do 200 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      logflg=(kd==kdintr.and.idis(nb)==0).or.
     &       (kd==kdsld .and.idis(nb)==0).or.
     &       (kd==kdbuff.and.idis(nb)==0).or.
     &       (kd==kdpors.and.idis(nb)==0).or.
     &       (kd==kdprdc.and.idis(nb)==0).or.
     &       (kd==kdshutr.and.idis(nb)==0)
      if(.NOT.logflg) then
! --- non-periodic BC:
          do 201 i=ivbcnd(nb-1)+1,ivbcnd(nb)
          iv=lvbcnd(i)
          lbc(iv)=1
  201     continue
          do 210 IS=1,nface_h
          if(lcface(2,IS).le.ncell_h) goto 210
          do 211 i=1,4
          if(lbc(lvface(i,IS)).lt.1) goto 210
  211     continue
          lbface(1,IS)=nb
  210     continue
          do 202 i=ivbcnd(nb-1)+1,ivbcnd(nb)
          iv=lvbcnd(i)
          lbc(iv)=0
  202     continue
! --- 
          if(idis(nb)==1.and.kd/=kdtchi) then
! --- Paired Discontinuous kdintr,kdsld,kdprdc,BC
            do id=ivpair(nb)+1,ivbcnd(nb)
              iv=lvbcnd(id)
              lbc1(iv)=-1
            enddo
            do IS=1,nface_h
              if(lcface(2,IS).le.ncell_h) cycle
              do i=1,4
                if(lbc1(lvface(i,IS))/=-1) goto 341
              enddo
              lbface(1,IS)=-nb
 341          continue
            enddo 
            do id=ivpair(nb)+1,ivbcnd(nb)
              iv=lvbcnd(id)
              lbc1(iv)=0
            enddo
          endif
!      elseif(kdbcnd(0,nb)==kdprdc.and.idis(nb)==0) then
!          do 220 IS=1,nface_h
!          if(abs(lbfacx(IS)).eq.nb) lbface(1,IS)=lbfacx(IS)
!  220     continue
      endif
  200 continue
!
! --- define partitioning interface (domain interface): 
! --- whole domain BC face check has been carried out
!???
      do 310 IS=1,nface_h
      nb=lbface(1,IS)
      if(lcface(2,IS).gt.ncell_h.and.nb.eq.0) then
        lbface(1,IS)=0
      endif
  310 continue
!
      deallocate(lbc,lbc1)
!
      call debug
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printi('lbface/list_facebc',lbface,2,mface)
      end subroutine debug
!
      end subroutine list_facebc_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facechk_hpc
     & (mface,mcell,nface_h,ncell_h,mvrtx,nvrtx_h,
     &  lvface,lcface,lbface,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only : ifle
!
! 1.  Check validity of boundary face
!
! 2.  Set dummy cell of lacell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mface,mcell,mvrtx
      integer,intent(in)    :: nface_h,ncell_h,nvrtx_h
      integer,intent(in)    :: lvface(4,mface)
      integer,intent(in)    :: lcface(2,mface)
      integer,intent(inout)    :: lbface(2,mface)
      integer,intent(inout) :: lacell(  mcell)
      integer,intent(inout) :: ierror
!
!
! --- [local entities]
!
      integer :: i,IS,ie,IC1,IC2,nb
!
      do 100 IS=1,nface_h
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      nb=lbface(1,IS)
      if(IC2.gt.ncell_h) then
! --- Check undefined face:
!        if(nb.eq.0) goto 9001
! --- Dumy cell IC2 must be same material no. with IC1:
        lacell(IC2)=lacell(IC1)
      else
! onishi debug =>zhang???
        if(nb.ne.0) lbface(1,IS)=0
!        if(nb.ne.0) goto 9002
        if( lacell(IC1).gt.0 .and. lacell(IC2).lt.0 ) goto 9003
        if( lacell(IC1).lt.0 .and. lacell(IC2).gt.0 ) goto 9003
      endif
  100 continue
!
      return
!
 9001 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) '### error 1 : data error'
      write(ifle,*) 'no boundary condtion is specified ',
     &              'for boundary face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',IC1
      goto 9999
 9002 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) '### error 2 : data error'
      write(ifle,*) 'boundary condtion can not be specified ',
     &              'for internal face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',lcface(1,IS)
      write(ifle,*) 'boundary condition no. =',lbface(1,IS)
      goto 9999
 9003 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) '### error 3 : data error'
      write(ifle,*) 'boundary condtion must be specified ',
     &              'for interface of fluid & solid'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',IC1,IC2
      goto 9999
 9999 continue
      write(ifle,*) '(list_facechk_hpc)'
      ierror=1
      end subroutine list_facechk_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facein_hpc
     & (mcell,mface,medge,mvrtx,ncell,nface,ncelb,nvrtx,NCV,
     &  lacell,lfcell,lvface,lcface,lbface,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only : ifle
      use module_boundary,only : kdbcnd,kdintr,kdprdc,kdbuff,kdshutr,
     &                           kdpors
      use module_partitioner,only : lvrtx=>WK1
      use module_partitioner,only : LCVVRT=>WK2
      use module_partitioner,only : lbfacx=>WK3
!
! 1. Make up list of connectivity for interface boundary 
!    2D cell for:
!            1) Liquid&Solid;
!            2) Liquid$Liquid; 
!            3) Solid&Solid
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,mface,medge,mvrtx,ncell,nvrtx
      integer,intent(inout) :: nface,ncelb
      integer,intent(in)    :: lacell(  mcell)
      integer,intent(inout) :: lfcell(7,mcell)
      integer,intent(inout) :: lvface(4,mface)
      integer,intent(inout) :: lcface(2,mface)
      integer,intent(inout) :: lbface(2,mface)
      INTEGER,INTENT(INOUT) :: LEFACE(5,MFACE)
      INTEGER,INTENT(INOUT) :: LVEDGE(2,MEDGE)
      integer,intent(inout) :: LVRTCV(  mvrtx)
      integer,intent(inout) :: LCVFAC(4,mface)
      integer,intent(inout) :: NCV
      integer,intent(out)   :: ierror
!
!
! --- [local entities]
!
      integer :: i,j,ie,nfacx,nfacy,ierr1=0,ierr2,NnewCV
      integer :: nb,m1,m2,n1,n2,n1d,n2d,kd
      integer :: IC,IS,IV,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE,ILC,ICV,ILCV
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC,IIS,IDC1,IDC2
!
      ierror=0
      return
!
      allocate(lvrtx(mvrtx),stat=ierr1) 
      if(ierr1.ne.0) 
     &   stop 'stop at allocating lvrtx(:) in list_facein_hpc'
      allocate(LCVVRT(mvrtx),stat=ierr1) 
      if(ierr1.ne.0)  
     &   stop 'stop at allocating LCVVRT(:) in list_facein_hpc'
      allocate(lbfacx(mface),stat=ierr1) 
      if(ierr1.ne.0)  
     &   stop 'stop at allocating lbfacx(:) in list_facein_hpc'    
!
!
      ierror=0
!
      lbfacx=0
      lvrtx=0
      LCVVRT=0
!
! ??? => LCYCCV has not been finished

!
      nfacx=nface
!
      do 100 IS=1,nfacx
      nb=lbface(1,IS)
      if(nb.lt.1) goto 100
      kd=kdbcnd(0,nb)
! --- if not interface BC: 
      if(.not.(kd==kdintr.or.kd==kdbuff.or.kd==kdshutr)) goto 100
      IIS=lbface(2,IS)
      ie=3+min(1,lvface(4,IS))
! --- if periodic BC
      if(IIS.gt.0 ) goto 101
! --- 
! / internal face /
!
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
! --- if cell IC2 has been defined as BC:
      if(IC2.gt.ncell) goto 9001
! --- if cell IC1 and cell IC2 are same material:
      if(lacell(IC1).eq.lacell(IC2)) goto 9002
!
      nface=nface+1
      ncelb=ncelb+2
!
      IDC1=ncelb-1
      IDC2=ncelb
      IIS =min(nface,mface)
!
! --- Redifine lcface(1:2,IS)
!
      lcface(1,IS)=IC1
      lcface(2,IS)=IDC1
      lbface(1,IS)=nb
      lbface(2,IS)=IIS
!
      lcface(1,IIS)=IC2
      lcface(2,IIS)=IDC2
      lbface(1,IIS)=nb
      lbface(2,IIS)=IS
!
      lvface(4,IIS)=0
      do 110 i=1,ie
      lvface(i,IIS)=lvface(ie+1-i,IS)
  110 continue
      ILC=0
      do 120 i=1,lfcell(7,IC2)
      if( lfcell(i,IC2).eq.IS ) ILC=i
  120 continue
      if( ILC.lt.1 ) then
        write(*,*) '### program error -1- (list_facein_hpc)'
        stop 999
      endif
      lfcell(ILC,IC2)=IIS
      goto 100
!
! / face on periodic boundary /
!
  101 continue
      IC1=lcface(1,IS)
      IC2=lcface(1,IIS)
      lbfacx(IIS)=nb
! --- IC1 and IC2 are same material:
      if(lacell(IC1).eq.lacell(IC2) ) goto 9002
!
  100 continue
!
      call list_facerr(ifle,mcell,mface,nface,ncelb,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      do 200 IS=1,nfacx
      if(lbfacx(IS).gt.0) then
        nb=abs(lbface(1,IS))
	kd=kdbcnd(0,nb)
        if(kd.ne.kdprdc) goto 9003
        lbface(1,IS)=lbfacx(IS)
      endif
  200 continue
!
!
!---------------------------------------------------------------------
! --- Creat new edges and new CVs
!---------------------------------------------------------------------
!
      do 300 IIS=nfacx+1,nface
      do 310 ILV=1,4
      IV=lvface(ILV,IIS)
      if(IV.gt.1) then
        lvrtx(IV)=lvrtx(IV)+1
      endif
  310 CONTINUE
  300 CONTINUE
!
      NnewCV=0
!
      do 320 IV=1,nvrtx
      if(lvrtx(IV).gt.0) then
        NnewCV=NnewCV+1
        ICV=NnewCV+nvrtx
        LVRTCV(ICV)=IV
        LCVVRT(IV)=ICV
      endif
  320 CONTINUE
!
      NCV=NCV+NnewCV
!
      do 330 IIS=nfacx+1,nface
      do 330 ILCV=1,4
      IV=lvface(ILCV,IIS)
      ICV=LCVVRT(IV)
      if(ICV.gt.1) LCVFAC(ILCV,IIS)=ICV
  330 CONTINUE
!
      deallocate(lvrtx,LCVVRT,lbfacx)
!
      call debug
      return
!
 9001 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'the face in interface boundary is boundary ',
     &              'face, not internal face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',IC1
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9002 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'the two cells on interface boundary has ',
     &              'the same attribute with each other'
      write(ifle,*) 'cell no. =',IC1,' and',IC2
      write(ifle,*) 'attribute no. =',lacell(IC1)
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      if( IIS.gt.0 )
     &write(ifle,*) '                     &',(lvface(i,IIS),i=1,ie)
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9003 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) '### error : data error'
      write(ifle,*) 'counterpart face in interface/periodic boundary ',
     &              'has a condition other than periodic boundary'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',lcface(1,IS)
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9999 continue
      write(ifle,*) '(list_facein_hpc)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printi('lfcell/list_facein_hpc',lfcell,7,mcell)
      call printi('lvface/list_facein_hpc',lvface,4,mface)
      call printi('lcface/list_facein_hpc',lcface,2,mface)
      call printi('lbface/list_facein_hpc',lbface,2,mface)
      end subroutine debug
!
      end subroutine list_facein_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facenn_hpc
     & (mcell,mface,medge,mvrtx,ncell,nface,ncelb,nvrtx,NCV,
     &  lfcell,lvface,lcface,lbface,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only          : ifle
      use module_boundary,only    : kdbinr,icbinr,lcbinr
      use module_partitioner,only : lvrtx=>WK1
      use module_partitioner,only : LCVVRT=>WK2
!
! 1. Make up list of connectivity for inner boundary (buffer)
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,mface,ncell,nvrtx,mvrtx,medge
      integer,intent(inout) :: nface,ncelb
      integer,intent(inout) :: lfcell(7,mcell)
      integer,intent(inout) :: lvface(4,mface)
      integer,intent(inout) :: lcface(2,mface)
      integer,intent(inout) :: lbface(2,mface)
      INTEGER,INTENT(INOUT) :: LEFACE(5,MFACE)
      INTEGER,INTENT(INOUT) :: LVEDGE(2,MEDGE)
      integer,intent(inout) :: LVRTCV(  mvrtx)
      integer,intent(inout) :: LCVFAC(4,mface)
      integer,intent(inout) :: NCV
      integer,intent(out)   :: ierror
!
!
! --- [local entities]
!
      integer :: i,j,ie,nfacx,ierr1=0
      integer :: nb,nb1,nb2,m1,m2,n1,n2,n1d,n2d,kdn
      integer :: IC,IS,IV,IVA,IVB,IC1,IC2,IDC1,IDC2,IIS,NnewCV,ICV,ILCV
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC,IV1,IV2,i1,i2,i3,i4
!
      allocate(lvrtx(mvrtx),stat=ierr1) 
      if(ierr1.ne.0) 
     &   stop 'stop at allocating lvrtx(:) in list_facen_hpcn'
      allocate(LCVVRT(mvrtx),stat=ierr1) 
      if(ierr1.ne.0) 
     &   stop 'stop at allocating LCVVRT(:) in list_facenn_hpc'
!
!
      ierror=0
!
      lvrtx=0
      LCVVRT=0
! ??? => LCYCCV has not been finished
!
      nfacx=nface
!
!
      do 100 IS=1,nfacx
      nb=lbface(1,IS)
      if(nb.lt.1) goto 100
      kdn=kdbinr(1,nb)
      if(kdn.lt.1) goto 100
!
      ie=3+min(1,lvface(4,IS))
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
! --- IC2.gt.ncell: IC2 is BC 2D-cell
      if(IC2.gt.ncell) goto 9001
!
      nb1=0
      nb2=0
      do 101 i=icbinr(nb-1)+1,icbinr(nb)
      if(IC1.eq.lcbinr(i)) THEN
! --- If BC lcbinr(i) is belong to IC1
        nb1=1
      ELSEif(IC2.eq.lcbinr(i)) THEN
! --- If BC lcbinr(i) is belong to IC2
        nb2=1
      ENDIF
  101 continue
      if(nb1.eq.nb2) goto 9002
      if(kdbinr(2,nb).ne.nb1+nb2+nb2 ) then
        nb1=kdbinr(1,nb)
        nb2=nb
      else
        nb1=nb
        nb2=kdbinr(1,nb)
      endif
!
      nface=nface+1
      ncelb=ncelb+2
!
      IDC1=ncelb-1
      IDC2=ncelb
      IIS =min(nface,mface)
!
      lcface(1,IS)=IC1
      lcface(2,IS)=IDC1
      lbface(1,IS)=nb1
!
      lcface(1,IIS)=IC2
      lcface(2,IIS)=IDC2
      lbface(1,IIS)=nb2
!
      lvface(4,IIS)=0
      do 120 i=1,ie
      lvface(i,IIS)=lvface(ie+1-i,IS)
  120 continue
      j=0
      do 130 i=1,lfcell(7,IC2)
      if( lfcell(i,IC2).eq.IS ) j=i
  130 continue
      if( j.lt.1 ) then
        write(*,*) '### program error -1- (list_facenn_hpc)'
        stop 999
      endif
      lfcell(j,IC2)=IIS
!
  100 continue
!
      call list_facerr(ifle,mcell,mface,nface,ncelb,ierr1)
      if( ierr1.ne.0 ) goto 9999
!

!---------------------------------------------------------------------
! --- Creat new edges and new CVs
!---------------------------------------------------------------------
!
      do 300 IIS=nfacx+1,nface
      do 310 ILV=1,4
      IV=lvface(ILV,IIS)
      if(IV.gt.1) then
        lvrtx(IV)=lvrtx(IV)+1
      endif
  310 CONTINUE
  300 CONTINUE
      NnewCV=0
      do 320 IV=1,nvrtx
      if(lvrtx(IV).gt.0) then
        NnewCV=NnewCV+1
        ICV=NnewCV+nvrtx
        LVRTCV(ICV)=IV
        LCVVRT(IV)=ICV
      endif
  320 CONTINUE
!
      NCV=NCV+NnewCV
!
      do 330 IIS=nfacx+1,nface
      do 330 ILCV=1,4
      IV=lvface(ILCV,IIS)
      ICV=LCVVRT(IV)
      if(ICV.gt.1) LCVFAC(ILCV,IIS)=ICV
  330 CONTINUE
!
      deallocate(lvrtx,LCVVRT)
!
      call debug
      return
!
 9001 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'the face in inner boundary is boundary ',
     &              'face, not internal face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',IC1
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9002 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'the two cells on inner boundary is located ',
     &              'at the same side of the boundary'
      write(ifle,*) 'cell no. =',IC1,' and',IC2
      write(ifle,*) 'side no. =',2-nb1
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9999 continue
      write(ifle,*) '(list_facenn_hpc)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printi('lfcell/list_facenn_hpc',lfcell,7,mcell)
      call printi('lvface/list_facenn_hpc',lvface,4,mface)
      call printi('lcface/list_facenn_hpc',lcface,2,mface)
      call printi('lbface/list_facenn_hpc',lbface,2,mface)
      end subroutine debug
!
      end subroutine list_facenn_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facepr_hpc
     & (mcell,mface,nface_h,mvrtx,nvrtx_h,ncell_h,mssfbc,
     &  IBCCYC_h,NBCpair,
     &  lvface,lcface,lbface,lacell,listpr,LEFACE,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only : ifle
      use module_boundary,only : nbcnd,ivbcnd,lvbcnd,kdbcnd,kdprdc,
     &                           nobcnd,kdsld,kdintr,idis,kdbuff,
     &                           kdshutr,kdpors,
     &                           numbname,boundName,boundIDMap
      use module_partitioner,only : listva=>WK1
      use module_partitioner,only : listvb=>WK2
      use module_partitioner,only : listsa=>WK3
      use module_partitioner,only : listsb=>WK4
!
! 1. Make up list of connectivity for periodic boundary
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mface,nface_h,nvrtx_h,ncell_h,mvrtx,mssfbc
      integer,intent(in)  :: NBCpair
      integer,intent(in)  :: mcell
      integer,intent(in)  :: lvface(4,mface)
      integer,intent(in)  :: lcface(2,mface)
      integer,intent(in)  :: lacell(  mcell)
      integer,intent(out) :: lbface(2,mface)
      integer,intent(out) :: listpr(0:3,mvrtx)
      integer,intent(out) :: IBCCYC_h
      integer,intent(in)  :: LEFACE(5,mface)
!
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer :: list_d(4,91)
      integer :: i,j,k,m,n,nb,nob,kd,nbuf,nbufa,nbufb
      integer :: j1,j2,m1,m2,n1,n2,ierr1
      integer :: nvfc,icmax,match,ie,nf,mch,ip0,ipe
      integer :: IC,IS,IS1,IS2,IV,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL,mb,iis
      integer :: IBS,IPS,IDC,IV1,IV2,i1,i2,i3,i4,isum1,isum2
!
      IBCCYC_h=0
      ierror=0
      do 20 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if((kd==kdprdc.and.idis(nb)==0).or.
     &   (kd==kdsld.and.idis(nb)==0) .or.
     &   (kd==kdbuff.and.idis(nb)==0) .or.
     &   (kd==kdpors.and.idis(nb)==0) .or.
     &   (kd==kdintr.and.idis(nb)==0).or.
     &   (kd==kdshutr.and.idis(nb)==0))
     &    goto 21
   20 continue
      return
   21 continue
!
      allocate(listva(0:mvrtx),stat=ierr1) 
      if(ierr1.ne.0) 
     &   stop 'stop at allocating listva(:) in list_facepr_hpc'
      allocate(listvb(0:mvrtx),stat=ierr1) 
      if(ierr1.ne.0)  
     &   stop 'stop at allocating listvb(:) in list_facepr_hpc'
      allocate(listsa(mface),stat=ierr1) 
      if(ierr1.ne.0)  
     &   stop 'stop at allocating listsa(:) in list_facepr_hpc' 
      allocate(listsb(mface),stat=ierr1) 
      if(ierr1.ne.0)  
     &   stop 'stop at allocating listsb(:) in list_facepr_hpc' 
!
! --- 
!
      listva(0)=1
      listvb(0)=2
!
!
      do 10 IS=1,mface
      lbface(1,IS)=0
      lbface(2,IS)=0
   10 continue
!
!----------------------------------------------------
! --- Periodic Face listing (Cell Center)
!----------------------------------------------------
!
      listpr=0
      do 1000 nb=1,nbcnd
      nob=nobcnd(nb)
      kd=kdbcnd(0,nb)
      if(kd==kdprdc.and.idis(nb)==0) then
!------------------------------------------------
! --- 1> list peroidic BC 
!------------------------------------------------
        do 1010 IV=1,mvrtx
        listva(IV)=0
        listvb(IV)=0
 1010   continue
!
        do 1020 IS=1,mface
        listsa(IS)=0
        listsb(IS)=0
 1020   continue  
!..................................................
!-< 1. Make up vertex list :listva,listvb,listpr >-
!..................................................
        ipe=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)
        do 1100 i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listva(IV1)=1
        listvb(IV2)=2

        listpr(0,IV1)=listpr(0,IV1)+1
        if(listpr(0,IV1)>3) stop 'list-facepr-hpc-1'
        listpr(0,IV2)=listpr(0,IV2)+1
        if(listpr(0,IV2)>3) stop 'list-facepr-hpc-2'

        nbufa=listpr(0,IV1)
        nbufb=listpr(0,IV2)
        listpr(nbufa,IV1)=IV2
        listpr(nbufb,IV2)=IV1
 1100   continue
!..................................................
!-< 2. Make up surface list :listsa >-
!..................................................
        nbufa=0
        nbufb=0
        do 1200 IS=1,nface_h
        if(lcface(2,IS).le.ncell_h) goto 1200
        do 1201 i=1,4
        ILV=lvface(i,IS)
        if(listva(ILV).ne.1) goto 1210
 1201   continue
        nbufa=nbufa+1
        listsa(nbufa)=IS
        goto 1200
 1210   continue
!
        do 1211 i=1,4
        ILV=lvface(i,IS)
        if(listvb(ILV).ne.2) goto 1220
 1211   continue
        nbufb=nbufb+1
        listsb(nbufb)=IS
 1220   continue
!
 1200   continue
!
        if(nbufb.ne.nbufa) then
          write(ifle,'(1x,a,2x,2I8,2x,a,2x,I4)') 
     &    '### ERR-P1 : periodic BC face No. of A & B not equal= '
     &    ,nbufa,nbufb,'at BC No. ',nob
          stop
        else
          write(ifle,'(1x,a,2x,I8,2x,a,2x,I8,2x,a,2x,I4)') 
     &   '    ###    FACE-A NO.=',
     &    NBUFA,' FACE-B NO.=',NBUFB,'AT BC NO. ',NOB
        endif
!...........................
!-< 3. matching procedure >-
!...........................
        mch=0
        do 1300 n1=1,nbufa
        IS1=listsa(n1)
        do 1400 n2=1,nbufb
        IS2=listsb(n2)
        nbuf=0
        list_d(1,:)=0
        list_d(2,:)=0
        list_d(3,:)=0
        list_d(4,:)=0
        do 1301 i1=1,listpr(0,lvface(1,IS2))
        do 1301 i2=1,listpr(0,lvface(2,IS2))
        do 1301 i3=1,listpr(0,lvface(3,IS2))
        if(lvface(4,IS2).eq.0) then
          nbuf=nbuf+1
          list_d(1,nbuf)=listpr(i1,lvface(1,IS2))
          list_d(2,nbuf)=listpr(i2,lvface(2,IS2))
          list_d(3,nbuf)=listpr(i3,lvface(3,IS2))
          list_d(4,nbuf)=0
        else
          do 1304 i4=1,listpr(0,lvface(4,IS2))
          nbuf=nbuf+1
          list_d(1,nbuf)=listpr(i1,lvface(1,IS2))
          list_d(2,nbuf)=listpr(i2,lvface(2,IS2))
          list_d(3,nbuf)=listpr(i3,lvface(3,IS2))
          list_d(4,nbuf)=listpr(i4,lvface(4,IS2))
 1304     continue
        endif
 1301   continue
        do 1405 j=1,nbuf
        call list_fmatch(match,lvface(1,IS1),list_d(1,j),-1)
        if( match.ge.1 ) goto 1303
 1405   continue
 1400   continue
        goto 9001
 1303   continue
        mch=mch+1
        lbface(1,IS1)= nb
        lbface(1,IS2)=-nb
        lbface(2,IS1)=IS2
        lbface(2,IS2)=IS1
 1300   continue
!
        if(nbufa.eq.mch) then
          write(ifle,*) '    ###    NO. OF PERIDIC PAIRS ', MCH
        else
          write(ifle,*) '*** Eroor: No. of Peridic pairs must be ', nbufa
        endif
        IBCCYC_h=mch+IBCCYC_h
!
      elseif(kd==kdsld.and.idis(nb)==1) then
      elseif(kd==kdsld.and.idis(nb)==0) then
!------------------------------------------------
! --- 2> list sliding BC
!------------------------------------------------
        listva(1:mvrtx)=0
        listvb(1:mvrtx)=0
!
        listsa(1:mface)=0
        listsb(1:mface)=0
!..................................................
!-< 1. Make up vertex list :listva,listvb,listpr >-
!..................................................
        ipe=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)
        do 2100 i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listva(IV1)=1
        listvb(IV2)=2
        listpr(0,IV1)=listpr(0,IV1)+1
        if(listpr(0,IV1)>3) stop 'list-facepr-hpc-1'
        listpr(0,IV2)=listpr(0,IV2)+1
        if(listpr(0,IV2)>3) stop 'list-facepr-hpc-2'
        nbufa=listpr(0,IV1)
        nbufb=listpr(0,IV2)
        listpr(nbufa,IV1)=IV2
        listpr(nbufb,IV2)=IV1
 2100   continue
!..................................................
!-< 2. Make up surface list :listsa >-
!..................................................
        nbufa=0
        nbufb=0
        do 2200 IS=1,nface_h
        if(lcface(2,IS).le.ncell_h) goto 2200
        do 2201 i=1,4
        ILV=lvface(i,IS)
        if(listva(ILV).ne.1) goto 2210
 2201   continue
        nbufa=nbufa+1
        listsa(nbufa)=IS
        goto 2200
 2210   continue
!
        do 2211 i=1,4
        ILV=lvface(i,IS)
        if(listvb(ILV).ne.2) goto 2220
 2211   continue
        nbufb=nbufb+1
        listsb(nbufb)=IS
 2220   continue
!
 2200   continue
!
        if(nbufb.ne.nbufa) then
          write(ifle,'(a,2I8)')
     &  'ERR: Sliding BC face No. of A & B not equal= '
     &    ,nbufa,nbufb
          stop 'list_facepr=> Sliding BC'
        endif
!..................................................
!-< 3. matching procedure >-
!..................................................
        mch=0
        do 2300 n1=1,nbufa
        IS1=listsa(n1)
        do 2400 n2=1,nbufb
        IS2=listsb(n2)
        nbuf=0
        list_d(1,:)=0
        list_d(2,:)=0
        list_d(3,:)=0
        list_d(4,:)=0
        do 2301 i1=1,listpr(0,lvface(1,IS2))
        do 2301 i2=1,listpr(0,lvface(2,IS2))
        do 2301 i3=1,listpr(0,lvface(3,IS2))
        if(lvface(4,IS2).eq.0) then
          nbuf=nbuf+1
          list_d(1,nbuf)=listpr(i1,lvface(1,IS2))
          list_d(2,nbuf)=listpr(i2,lvface(2,IS2))
          list_d(3,nbuf)=listpr(i3,lvface(3,IS2))
          list_d(4,nbuf)=0
        else
          do 2304 i4=1,listpr(0,lvface(4,IS2))
          nbuf=nbuf+1
          list_d(1,nbuf)=listpr(i1,lvface(1,IS2))
          list_d(2,nbuf)=listpr(i2,lvface(2,IS2))
          list_d(3,nbuf)=listpr(i3,lvface(3,IS2))
          list_d(4,nbuf)=listpr(i4,lvface(4,IS2))
 2304     continue
        endif
 2301   continue
        do 2405 j=1,nbuf
        call list_fmatch(match,lvface(1,IS1),list_d(1,j),-1)
        if( match.ge.1 ) goto 2303
 2405   continue
 2400   continue
        goto 9001
 2303   continue
        mch=mch+1
        lbface(1,IS1)= nb
        lbface(1,IS2)=-nb
        lbface(2,IS1)=IS2
        lbface(2,IS2)=IS1
 2300   continue
!
        if(nbufa.eq.mch) then
          write(ifle,'(a,I8)')
     &   'MSG: NO. OF Sliding FACE PAIRS ', MCH
        else
          write(ifle,'(a,I8)')
     &   'ERR: NO. OF Sliding PAIRS MUST BE ', NBUFA
          stop
        endif
        DO i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listpr(:,IV1)=0
        listpr(:,IV2)=0
        ENDDO
!
      elseif((kd==kdintr.or.kd==kdbuff.or.kd==kdshutr.or.kd==kdpors)
     &   .and.idis(nb)==0) then 
!------------------------------------------------
! --- 3> list interface BC
!------------------------------------------------
        listva(1:mvrtx)=0
        listvb(1:mvrtx)=0
!
        listsa(1:mface)=0
        listsb(1:mface)=0
!..................................................
!-< 1. Make up vertex list :listva,listvb,listpr >-
!..................................................
        ipe=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)
        do 3100 i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listva(IV1)=1
        listvb(IV2)=2
        listpr(0,IV1)=listpr(0,IV1)+1
        if(listpr(0,IV1)>3) stop 'list-facepr-hpc-1'
        listpr(0,IV2)=listpr(0,IV2)+1
        if(listpr(0,IV2)>3) stop 'list-facepr-hpc-2'
        nbufa=listpr(0,IV1)
        nbufb=listpr(0,IV2)
        listpr(nbufa,IV1)=IV2
        listpr(nbufb,IV2)=IV1
 3100   enddo
!..................................................
!-< 2. Make up surface list :listsa >-
!..................................................
        nbufa=0
        nbufb=0
        do 3200 IS=1,nface_h
        if(lcface(2,IS).le.ncell_h) goto 3200
        if(lbface(2,IS)>0) goto 3200
!
!        do 3201 i=1,4
!        ILV=lvface(i,IS)
!        if(listva(ILV).ne.1) goto 3210
! 3201   enddo
!
!        do 3211 i=1,4
!        ILV=lvface(i,IS)
!        if(listvb(ILV).ne.2) goto 3220
! 3211   enddo
!
! 3220   continue
!
!        ISUM1=lvface(1,IS)+
!     &        lvface(2,IS)+
!     &        lvface(3,IS)+
!     &        lvface(4,IS)
!        do iis=1,NBCpair
!        mb=LEFACE(5,iis)
!        if(boundIDMap(nb,1)==mb) then
!          ISUM2=LEFACE(1,iis)+
!     &          LEFACE(2,iis)+
!     &          LEFACE(3,iis)+
!     &          LEFACE(4,iis)
!          if(ISUM1==ISUM2) then
!            nbufa=nbufa+1
 !           listsa(nbufa)=IS
!            goto 3250
!          endif
!        endif
! 3210   continue
!
!        if(boundIDMap(nb,2)==mb) then
!          ISUM2=LEFACE(1,iis)+
!     &          LEFACE(2,iis)+
!     &          LEFACE(3,iis)+
!     &          LEFACE(4,iis)
!          if(ISUM1==ISUM2) then
!            nbufb=nbufb+1
!            listsb(nbufb)=IS
!            goto 3250
!          endif
!        endif
!        enddo
! 3250   continue
!!
        do 3201 i=1,4
        ILV=lvface(i,IS)
        if(listva(ILV).ne.1) goto 3210
        
 3201   enddo
        nbufa=nbufa+1
        listsa(nbufa)=IS
        goto 3200
 3210   continue

        do 3211 i=1,4
        ILV=lvface(i,IS)
        if(listvb(ILV).ne.2) goto 3220
        if(lbface(2,IS)>0) goto 3200
 3211   enddo
        nbufb=nbufb+1
        listsb(nbufb)=IS
 3220   continue
!
 3200   enddo
!
!        if(nbufb/=nbufa) then
!          write(ifle,'(a,2I10,a,I4)')
!     &     "ERR: Interface BC face No. of A & B not equal= ",
!     &   nbufa,nbufb,' at BC No.= ',nb
!          write(ifle,'(a)') 'Contact your supportor'
!          stop 'list_facepr=> Interface BC'
!        endif
!..................................................
!-< 3. matching procedure >-
!..................................................
        mch=0
        do 3300 n1=1,nbufa
        IS1=listsa(n1)
        do 3400 n2=1,nbufb
        IS2=listsb(n2)
        if(lacell(lcface(1,IS1))==lacell(lcface(1,IS2))) cycle
        nbuf=0
        list_d(1,:)=0
        list_d(2,:)=0
        list_d(3,:)=0
        list_d(4,:)=0
        do 3301 i1=1,listpr(0,lvface(1,IS2))
        do 3301 i2=1,listpr(0,lvface(2,IS2))
        do 3301 i3=1,listpr(0,lvface(3,IS2))
        if(lvface(4,IS2).eq.0) then
          nbuf=nbuf+1
          list_d(1,nbuf)=listpr(i1,lvface(1,IS2))
          list_d(2,nbuf)=listpr(i2,lvface(2,IS2))
          list_d(3,nbuf)=listpr(i3,lvface(3,IS2))
          list_d(4,nbuf)=0
        else
          do 3304 i4=1,listpr(0,lvface(4,IS2))
          nbuf=nbuf+1
          list_d(1,nbuf)=listpr(i1,lvface(1,IS2))
          list_d(2,nbuf)=listpr(i2,lvface(2,IS2))
          list_d(3,nbuf)=listpr(i3,lvface(3,IS2))
          list_d(4,nbuf)=listpr(i4,lvface(4,IS2))
 3304     enddo
        endif
 3301   continue
        do 3405 j=1,nbuf
        call list_fmatch(match,lvface(1,IS1),list_d(1,j),-1)
        if( match.ge.1 ) goto 3303
 3405   enddo
 3400   enddo
!        goto 9001
        goto 3300
 3303   continue
        mch=mch+1
! --- Only for Solid/Fluid 
        if(kd==kdintr) then
        if(lacell(lcface(1,IS1))==lacell(lcface(1,IS2))) then
          write(ifle,'(a)') 'ERR: Interface BC error'
          stop 'Contact your supportor'
        endif
        if(lacell(lcface(1,IS1))>lacell(lcface(1,IS2))) then
          lbface(1,IS1)= nb
          lbface(1,IS2)=-nb
        else
          lbface(1,IS1)=-nb
          lbface(1,IS2)= nb
        endif
        lbface(2,IS1)=IS2
        lbface(2,IS2)=IS1
        else
          lbface(1,IS1)= nb
          lbface(1,IS2)=-nb
          lbface(2,IS1)=IS2
          lbface(2,IS2)=IS1
        endif
 3300   enddo
!
!        if(nbufa.eq.mch) then
!          write(ifle,'(a,I8,a,I4)') 
!     &   'MSG: NO. OF Interface FACE PAIRS ', MCH,
!     &   ' at BC No.=',nb
!        else
!          write(ifle,'(a,I8)')
!     &    'ERR: NO. OF Interface PAIRS MUST BE ', NBUFA
!          stop
!        endif
          write(ifle,'(a,I8,a,I4)') 
     &   'MSG: NO. OF Interface FACE PAIRS ', MCH,
     &   ' at BC No.=',nb
!
        DO i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listpr(:,IV1)=0
        listpr(:,IV2)=0
        ENDDO
!        
      endif
 1000 continue
!
      deallocate(listva,listvb,listsa,listsb)
!      deallocate(boundIDMap)
!
      call debug
!
      return
!
 9001 continue
      ie=3+min(1,lvface(4,IS1))
      write(ifle,*) '### error : data error'
      write(ifle,*) 'the face in periodic boundary has ',
     &              'no couterpart face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS1),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',lcface(1,IS1)
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9999 continue
      write(ifle,*) '(list_facepr_hpc)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printi('lbface/list_facepr_hpc',lbface,2,mface)
      end subroutine debug
!
      end subroutine list_facepr_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_fluid_hpc
     & (mcell,mface,nface_h,ncell_h,ncelb_h,mvrtx,nvrtx_h,
     &  lbface,lcface,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only          : ifle
      use module_material,only       : strdata,lclsd
      use module_partitioner,only : iclos=>WK1
      use module_partitioner,only : kdbp=>WK2
!
! 1.  Make up list "lacell" for multiple fluid domain
!     separated by solid part
!
! 2.  Set flag for closed domain
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,mface,nface_h,ncell_h,ncelb_h
      integer,intent(in)    :: nvrtx_h,mvrtx
      integer,intent(in)    :: lbface(2,mface)
      integer,intent(in)    :: lcface(2,mface)
      integer,intent(inout) :: lacell(  mcell)
      integer,intent(out)   :: ierror
!
!
! --- [local entities]
!
!      integer :: iclos (  mcell)
!      integer :: kdbp  (  mface)
!
      integer :: jclos (  0:100 )
      integer :: MMAT(-100:100),MATC_NO(200)
      integer :: MATC_SLD(100),MATC_FLD(100)
      integer :: i,j,m,n,la,ll,nf,nc
! delete already given attribute -- by onishi
!      integer :: i,j,m,n,la,ll,nf,nc,IC,IS
      integer :: nxx,ns,ne,nd,iflg
      integer :: ierr1,ierr,nflud,ihpc,nfludx,NSLIDX
      integer :: IC,IS,IV,IE,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL,ICFLD,NNMAT
      integer :: IBS,IPS,IDC,IIS,IIMAT
!
!
!
!-< 1. Preliminary set >-
!
      allocate(iclos(mcell),stat=ierr1) 
      if(ierr1.ne.0) 
     &   stop 'stop at allocating nclv(:) in list_fluid_hpc'
      allocate(kdbp  (  mface),stat=ierr1) 
      if(ierr1.ne.0)
     &   stop 'stop at allocating ityp(:) in list_fluid_hpc'
!
!--< 1.1 make up list local "lfcell" ; face no. surrounding cell >--
!
      ierror=0
      ierr=0
      ierr1=0
!
!-< 3. Set flag for closed domain >-
!
      call bc_kdbpv1(mface,nface_h,lbface,kdbp)
!
      iclos(:)=0
      MMAT=0
      do 100 IC=1,ncell_h
      IMAT=lacell(IC)
      MMAT(IMAT)=1
 100  continue
!
      if(MMAT(0).eq.1) then
        write(ifle,*) '### error: Some cells has NOT been defined'
        goto 9999
      endif
!
      NNMAT=0
      do 200 IMAT=-100,100
      if(MMAT(IMAT).eq.1) then
        NNMAT=NNMAT+1
        MATC_NO(NNMAT)=IMAT
      endif
 200  continue
!
      NFLUDX=0
      NSLIDX=0
      MATC_FLD=0
      MATC_SLD=0
      do 300 IIMAT=1,NNMAT
      IMAT=MATC_NO(IIMAT)
      if(IMAT.gt.0) then
        NFLUDX=NFLUDX+1
        MATC_FLD(NFLUDX)=IIMAT
      elseif(IMAT.lt.0) then
        NSLIDX=NSLIDX+1
        MATC_SLD(NSLIDX)=IIMAT
      endif
 300  continue
!
      jclos(:)=0
!
      do 310 IS=1,nface_h
      if( kdbp(IS).eq.1 ) then
! --- Outlet BC:
        IC1=lcface(1,IS)
        IC2=lcface(2,IS)
        iclos(IC1)=1
        iclos(IC2)=1
      endif
  310 continue
!
      do 320 IC=1,ncell_h
      IMAT=lacell(IC)
      if(IMAT.gt.0) then
        if(iclos(IC).eq.1) then
          jclos(IMAT)=1
        endif
      endif
  320 continue
!
! --- 
!
      ihpc=1
      deallocate(lclsd)
      call strdata
     &  (ifle,NFLUDX,MATC_FLD,MATC_NO,jclos,ihpc,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      deallocate(iclos,kdbp)  
!
      call debug
      return
!
 9999 continue
      write(ifle,*) '(list_fluid_hpc)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printj('lacell/list_fluid_hpc',lacell,1,mcell)
      end subroutine debug
!
      end subroutine list_fluid_hpc
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_solid_hpc(mcell,ncell_h,ncelb_h,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only      : ifle
      use module_material,only   : nsold,nosld
      use module_material,only   : nflud,nofld
!
!-< 1. Check solid no. & renumber lacell >-
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,ncell_h,ncelb_h
      integer,intent(inout) :: lacell(mcell)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer :: i,n,IC,isld,ifld
!
      ierror=0
!
      do 100 IC=1,ncelb_h
      if( lacell(IC).lt.0 ) then
        do 101 isld=1,nsold
        if( nosld(isld).eq.lacell(IC) ) then
          lacell(IC)=-isld
          goto 100
        endif
  101   continue
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of control data for solid part'
        write(ifle,*) 'solid no. =',lacell(IC)
        goto 9999
      elseif(lacell(IC).gt.0) then 
        do 102 ifld=1,nflud
        if( nofld(ifld).eq.lacell(IC) ) then
          lacell(IC)=ifld
          goto 100
        endif
  102   continue
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of control data for fluid part'
        write(ifle,*) 'fluid no. =',lacell(IC)
        goto 9999
      elseif(lacell(IC).eq.0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'Element IC= ',IC,' is lack of Material NO.'
        goto 9999
      endif
  100 continue
!
      return
!
 9999 continue
      write(ifle,*) '(list_solid_hpc)'
      ierror=1
!
      end subroutine list_solid_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_BC_hpc
     & (mssfbc,nssfbc_h,NBCINL_h,NBCWAL_h,NBCOUT_h,NBCCYC_h,NBCSYM_h,
     &  NBCTCI_h,NBCSLD_h,NBCINT_h,LBCSSF)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_boundary,only : kdbcnd,kdilet,kvlglw,kvnslp,kdtchi,
     &                           kxnone,kdolet,kdprdc,kdsymm,kdsld,
     &                           kdintr,idis,kvmodl,kvrogf,kdbuff,
     &                           kdshutr,kdpors
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mssfbc,nssfbc_h
      integer,intent(out)   :: NBCINL_h,NBCWAL_h,NBCOUT_h,NBCCYC_h
      integer,intent(out)   :: NBCSYM_h,NBCTCI_h,NBCSLD_h,NBCINT_h
      integer,intent(in)    :: LBCSSF(  mssfbc)
!
! --- [local entities]
!
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: NBCINL,NBCWAL,NBCOUT,NBCCYC,NBCSYM,NBCTCI,NBCSLD,
     &           NBCINT
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KDV,KD,IDC,IN,IBFI,ICFI
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
!
!
      NBCINL=0
      NBCWAL=0
      NBCOUT=0
      NBCCYC=0
      NBCSYM=0
      NBCTCI=0
      NBCSLD=0
      NBCINT=0
!
      do 1000 IBF=1,nssfbc_h
      nb=LBCSSF(IBF)
      if (nb.le.0) goto 1000
      kd=kdbcnd(0,nb)
      kdv=kdbcnd(1,nb)
      if(kd.eq.kdilet) then
        NBCINL=NBCINL+1
!???      elseif(kd.eq.kxnone) then
      elseif
     & (kdv==kvnslp.or.kdv==kvlglw.or.kd==kxnone.or.kdv==kvrogf.or.
     &   kdv==kvmodl) then
        if(LBCSSF(IBF)>0) then
          NBCWAL=NBCWAL+1
        endif
!        if((kd.eq.kdintr.or.kd==kdshutr).and.LBCSSF(IBF)>0) then
!          NBCWAL=NBCWAL+1
!        endif
      elseif(kd.eq.kdolet) then
        NBCOUT=NBCOUT+1
      elseif(kd.eq.kdprdc.and.idis(nb)==0) then
        NBCCYC=NBCCYC+1
      elseif(kd.eq.kdsymm) then
        NBCSYM=NBCSYM+1
      elseif(kd.eq.kdtchi.and.idis(nb)==0) then
        NBCTCI=NBCTCI+1
      elseif(kd.eq.kdsld) then
        NBCSLD=NBCSLD+1
      elseif((kd==kdintr.or.kd==kdbuff.or.kd==kdshutr.or.kd==kdpors)
     &     .and.idis(nb)==0) then
        NBCINT=NBCINT+1
      endif
 1000 continue
      NBCINL_h=NBCINL
      NBCWAL_h=NBCWAL
      NBCOUT_h=NBCOUT
      NBCCYC_h=NBCCYC
      NBCSYM_h=NBCSYM
      NBCTCI_h=NBCTCI
      NBCSLD_h=NBCSLD
      NBCINT_h=NBCINT
!
!
      return
      end subroutine list_BC_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_check_hpc
     & (mface,mcell,nface_h,mvrtx,nvrtx_h,
     & lbface,lcface,lvface,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only       : ifle
      use module_model,only    : icaltb,rns_scl
      use module_boundary,only : kdbcnd,kdprdc,kdsymm,kdintr,kdsld,
     &                           kdilet,kdolet,kdtchi,kxnone,kdpres,
     &                           kdstag,idis,kdbuff,kdshutr,kdpors,
     &                           kdovst
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mface,mcell,nface_h,nvrtx_h,mvrtx
      integer,intent(in)  :: lbface(2,mface)
      integer,intent(in)  :: lcface(2,mface)
      integer,intent(in)  :: lvface(4,mface)
      integer,intent(in)  :: lacell(  mcell)
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
!
!
!
      ierror=0
!
      do 1000 IS=1,nface_h
      nb=lbface(1,IS)
      if(nb.gt.0) then
        kd=kdbcnd(0,nb)
        if(kd.ne.kdsymm.and.kd.ne.kdprdc.and.kd.ne.kdtchi.and.
     &     kd.ne.kdsld.and.kd/=kdovst) then
!
          IC1=lcface(1,IS)
          kdv=kdbcnd(1,nb)
          kdt=kdbcnd(2,nb)
          kdy=kdbcnd(3,nb)
          kdk=kdbcnd(4,nb)
          kdp=kxnone
!
!-< 1. fluid face >-
!
          IMAT=lacell(IC1)
          if(IMAT.gt.0) then
            if(kd.ne.kdilet 
     &        .and.kd.ne.kdolet
     &        .and.kd.ne.kdpres
     &        .and.kd.ne.kdstag
     &        .and.kd/=kdbuff
     &        .and.kd/=kdshutr
     &        .and.kd/=kdpors) then
              call errmsg1(1,kdv,IC1,'velocity')
              call errmsg1(1,kdt,IC1,'temperature')
              call errmsg1(1,kdy,IC1,'mass fraction')
              if(rns_scl) call errmsg1(1,kdk,IC1,'RANS')
              if(ierror.ne.0 ) goto 9999
            endif
          endif
!
!-< 2. solid face >-
!
          if(IMAT.lt.0) then
            if(kd.eq.kdilet) kdv=kxnone+1
            if(kd.eq.kdolet) kdp=kxnone+1
            call errmsg1(2,kdt,IC1,'temperature')
            if(kd.ne.kdintr
     &    .and.kd/=kdbuff
     &    .and.kd/=kdshutr
     &    .and.kd/=kdpors) then
              call errmsg2(2,kdv,IC1,'velocity')
              call errmsg2(2,kdp,IC1,'pressure')
              call errmsg2(2,kdy,IC1,'mass fraction')
              if(rns_scl) call errmsg2(2,kdk,IC1,'RANS')
            endif
            if(ierror.ne.0) goto 9999
          endif
!
        endif
      endif
 1000 continue
!
      return
!
 9999 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) 'vertices of the face 1=',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'boundary condition no. =',nb
      write(ifle,*) '(bc_check_hpc)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine errmsg1(ifs,kd,n1,vnam)
!=================================================
      implicit none
      integer     ,intent(in) :: ifs,kd,n1
      character(*),intent(in) :: vnam
      character(5),parameter  :: cfs(2)=(/'fluid','solid'/)
      if( kd.ne.kxnone ) return
      write(ifle,*) '### error : data error'
      write(ifle,*) 'boundary condition for ',vnam,' is not ',
     &              'specified on ',cfs(ifs),' face'
      write(ifle,*) 'cell no. =',n1
      ierror=1
      end subroutine errmsg1
!=================================================
      subroutine errmsg2(ifs,kd,n1,vnam)
!=================================================
      implicit none
      integer     ,intent(in) :: ifs,kd,n1
      character(*),intent(in) :: vnam
      character(5),parameter  :: cfs(2)=(/'fluid','solid'/)
      if( kd.eq.kxnone ) return
      write(ifle,*) '### error : data error'
      write(ifle,*) 'boundary condition for ',vnam,' can not ',
     &              'be specified on ',cfs(ifs),' face'
      write(ifle,*) 'cell no. =',n1
      ierror=1
      end subroutine errmsg2
      end subroutine bc_check_hpc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

