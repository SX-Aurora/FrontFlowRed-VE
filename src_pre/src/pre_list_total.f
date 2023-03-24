!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!      subroutine list_face
!      subroutine list_edge
!      subroutine list_facebc
!      subroutine list_facechk
!      subroutine list_facein
!      subroutine list_facenn
!      subroutine list_facebf
!      subroutine list_facepr
!      subroutine list_facerr
!      subroutine list_fluid 
!      subroutine list_fmatch
!      subroutine list_ematch
!      subroutine list_onvrtx
!      subroutine list_solid 
!      subroutine list_edgvt
!      subroutine list_BCDIM 
!      subroutine list_output_geom
!      subroutine list_touch_inlet
!      subroutine list_pair
!      subroutine list_wall_distance
!      subroutine list_nwvrtx
!      subroutine list_output_sld
!      subroutine list_output_grid 
!      subroutine list_output_ini_movegrid 
!      subroutine list_output_grid_A
!      subroutine list_output_grid_U
!      subroutine list_gridfact
!      subroutine list_BC_connectivity
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_face
     & (mcell,mface,mvrtx,nvrtx,ncell,nface,ncelb,NCV,
     & lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,lbface,cord,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only          : ifle,ifll
      use module_partitioner,only : nclv =>WK1
      use module_partitioner,only : ityp =>WK2
      use module_partitioner,only : ip   =>WK3
      use module_partitioner,only : lvrtx=>IW1
      use module_parameter,  only : NIFACE
      use module_boundary,only    : nbcnd,ivbcnd,lvbcnd,NBOUND,thinw
      use module_model,only       : vertex_cen,cell_cen,icon_cv
!
! --- 1. Make up list of connectivity for cell & face
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mcell,mface,ncell,mvrtx
      integer,intent(out) :: nvrtx
      integer,intent(out) :: nface,ncelb
      integer,intent(inout)  :: lvcell(8,mcell)
      integer,intent(out) :: lfcell(7,mcell)
      integer,intent(out) :: lvface(4,mface)
      integer,intent(out) :: lcface(2,mface)
      integer,intent(out) :: LVRTCV(  mvrtx)
      integer,intent(out) :: LCVFAC(4,mface)
      integer,intent(out) :: lbface(2,mface)
      integer,intent(out) :: ierror,NCV
      real*8 ,intent(out) :: cord(  3,mvrtx)
!
! --- [local entities]
!
      integer,parameter :: nv2typ(8)=(/0,0,0,1,2,3,0,4/)
      integer,parameter :: nfacel(4)=(/4,5,5,6/)
! --- define local face vertex order
      integer,parameter :: lvfcel(4,6,4)=reshape( source=
     &  (/1,3,2,0, 2,3,4,0, 3,1,4,0, 4,1,2,0, 0,0,0,0, 0,0,0,0,
     &    1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0, 0,0,0,0,
     &    1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2,0, 4,5,6,0, 0,0,0,0,
     &    1,5,8,4, 2,3,7,6, 1,2,6,5, 3,4,8,7, 1,4,3,2, 5,6,7,8/),
     &  shape=(/4,6,4/) )
!
!       nv2typ : type of cell according to face no.
!              : =1;tetrahedron, =2;pyramid, =3;prism, =4;hexahedron
!       nfacel : no. of faces in cell
!       lvfcel : dummy vertex no. in cell face
!
      integer :: lvf(4,6,4),nb,nvrtxnw
      integer :: lvrt(4),LBC(nbcnd)
      integer :: i,j,k,l,m,n,IS,IVV,idum
      integer :: j1,j2,n1,n2,m1,m2,lf1,ie,ierr1=0,ierr2
      integer :: kclv,nf,icmax,match,nch,nfacx,nclbx
      integer :: IC,IV,ILV,ILS,NF1,NF2,ILS1,ILS2,IC1,IC2,ICTP
      integer :: ID,IFACE,IS2,IFACE2
      logical :: INNER
!
! --- 
!
      allocate(nclv(8*ncell),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating nclv(:) in list_face'
      allocate(ityp(ncell),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ityp(:) in list_face'
      allocate(ip(0:mvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ip(:) in list_face'

      if(icon_cv==vertex_cen) then
        NCV=NVRTX
      endif
!-------------------------------
! --- < 1. Set up some arrays >-
!-------------------------------
! --- < 1.1 make up "lvf"     >-
!-------------------------------
      ierror=0
      ierr1=0
      ierr2=0
!
      do 100 k=1,4
      do 100 j=1,6
      do 101 i=1,4
      lvf(i,j,k)=lvfcel(i,j,k)
  101 continue
      if(lvf(4,j,k).lt.1) lvf(4,j,k)=8
  100 continue
!
! --- < 1.2 clear array >--
!
      do IC=1,mcell
      do i=1,7
      lfcell(i,IC)=0
      enddo
      enddo
!
      do IS=1,mface
      lcface(1,IS)=0
      lcface(2,IS)=0
      do i=1,4
      lvface(i,IS)=0
      LCVFAC(i,IS)=0
      enddo
      enddo
!--------------------------------------
! --- < 1.3 set type & no. of faces >--
!--------------------------------------
      do 120 IC=1,ncell
      kclv=0
      do ILV=1,8
      if(lvcell(ILV,IC).gt.0) kclv=ILV
      enddo
      if(kclv.lt.1) then
        write(ifle,'(a,I10)') 'ERR: kclv= ',kclv
        goto 9001
      endif
      idum=nv2typ(kclv)
      ityp(IC)=idum
      if(idum.lt.1) then
         write(ifle,'(a,3I10)') 'ERR: ityp(IC)= ',ityp(IC),kclv,IC
         write(ifle,'(a,I10)') 'MSG: lvcell(1:8,IC)= ',lvcell(:,IC)
         goto 9001
      endif
      lfcell(7,IC)=nfacel(ityp(IC))
  120 enddo
!-----------------------------------
! --- < 2. Set up cell-face list >--
!----------------------------------------------------------------
! --- < 2.1 make up list for cells connected with each vertex >--
!----------------------------------------------------------------
      call list_onvrtx(8,mcell,nvrtx,ncell,
     & lvcell,ip,nclv,icmax)
!
      allocate(lvrtx(6,4*icmax),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating lvrtx in list_face'
!-----------------------------------
! --- < 2.2 Search internal faces >-
! --- icmax is maximun no. of cells which connect one vertex.
!-----------------------------------
      nface=0
!
      do 200 IV=1,nvrtx
!------------------------------------
! --- / set vertex of faces in cell /
!------------------------------------
      nf=0
      do 210 IVV=ip(IV-1)+1,ip(IV)
      IC=nclv(IVV)
      ICTP=ityp(IC)
      do 211 ILS=1,lfcell(7,IC)
      if(lfcell(ILS,IC).gt.0) goto 211
      do ILV=1,4
      lvrt(ILV)=lvcell(lvf(ILV,ILS,ICTP),IC)
      enddo
!
      INNER=lvrt(1).EQ.IV
     &  .OR.lvrt(2).EQ.IV
     &  .OR.lvrt(3).EQ.IV
     &  .OR.lvrt(4).EQ.IV
      IF(INNER) THEN
        nf=nf+1
        if(nf>4*icmax) stop 8866
        do ILV=1,4
        lvrtx(ILV,nf)=lvcell(lvf(ILV,ILS,ICTP),IC)
        enddo
        lvrtx(5,nf)=ILS
        lvrtx(6,nf)=IC
      ENDIF
  211 continue
  210 continue
!-----------------------------------
! --- / matching procedure /
!-----------------------------------
      do 220 nf1=1,nf-1
      ILS1=lvrtx(5,nf1)
      IC1=lvrtx(6,nf1)
      if(lfcell(ILS1,IC1).gt.0) goto 220
      do 221 nf2=nf1+1,nf
      ILS2=lvrtx(5,nf2)
      IC2=lvrtx(6,nf2)
!      if(IC1==IC2) goto 221
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
      do i=1,4
      lvface(i,nfacx)=lvrtx(i,nf1)
      enddo
      goto 220
  221 continue
  220 continue
!-----------------------------------
! --- / procedure only for check /
!-----------------------------------
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
        if(match.gt.0) goto 9003
      else
        nch=IC2
      endif
  231 continue
  230 continue
!
  200 continue
!
      write(ifll,'(a,I8)') 'MSG: INTER FACE => NIFACE= ',NFACE
      NIFACE=Nface
!-------------------------------
!-< 3. Search boundary faces >--
!-------------------------------
      ncelb=ncell
      do 300 IC=1,ncell
      k=ityp(IC)
      do 301 IS=1,lfcell(7,IC)
      if(lfcell(IS,IC).lt.1) then 
        ncelb=ncelb+1
        nface=nface+1
        nclbx=min(ncelb,mcell)
        nfacx=min(nface,mface)
!
        lfcell(IS,IC)=nface
        lcface(1,nfacx)=IC
        lcface(2,nfacx)=ncelb
        do 302 iv=1,4
        lvface(iv,nfacx)=lvcell(lvf(iv,IS,k),IC)
  302   continue
      endif
  301 continue
  300 continue
!---------------------------------------
! --- BC thin wall(same vertex number)
!---------------------------------------
      if(.false.) then
      nvrtxnw=nvrtx

!      ip(:)=0
!      do iv=1,nvrtx
!      ip(iv)=-iv 
!      enddo
!
!      do nb=1,nbcnd 
!      if(thinw(nb)==1) then
!        do id=ivbcnd(nb-1)+1,ivbcnd(nb)
!        iv=lvbcnd(id)
!        if(ip(iv)<0) then
!          ip(iv)=IV
!!        elseif(ip(iv)>0) then
!!          nvrtxnw=nvrtxnw+1
!!          ip(iv)=nvrtxnw
!        endif
!        enddo
!      endif
!      enddo
!--------------------------------
! --- 
!--------------------------------
      
      do nb=1,nbcnd 
      ip(0)=1
      if(thinw(nb)==1) then
        ip(1:nvrtx)=0 
        do id=ivbcnd(nb-1)+1,ivbcnd(nb)
        iv=lvbcnd(id)
        ip(iv)=IV
        ENDDO


        do 400 IC=1,ncell
        k=ityp(IC)
        do 401 IS=1,lfcell(7,IC)
        IFACE=lfcell(IS,IC)
        if(lcface(2,IFACE)>ncell.and.lbface(1,IFACE)>=0)  goto 401
        do i=1,4
        if(ip(lvface(i,IFACE)).lt.1) goto 401
        enddo
        if(lbface(1,IFACE)==0) then
          IC1=lcface(1,IFACE)
          IC2=lcface(2,IFACE)
        
          ncelb=ncelb+1
          lfcell(IS,IC)=IFACE
          lcface(1,IFACE)=IC1
          lcface(2,IFACE)=ncelb
          lbface(1,IFACE)=nb
          do iv=1,4
          lvface(iv,IFACE)=lvcell(lvf(iv,IS,k),IC)
          enddo
!
          ncelb=ncelb+1
          nface=nface+1
          nclbx=min(ncelb,mcell)
          nfacx=min(nface,mface)
          do IS2=1,lfcell(7,IC2)
          IFACE2=lfcell(IS2,IC2)
          if(IFACE2==IFACE) then
            lfcell(IS2,IC2)=nface
            exit
          endif
          enddo
          lcface(1,nfacx)=IC2
          lcface(2,nfacx)=ncelb
          lbface(1,nfacx)=-nb
        elseif(lbface(1,IFACE)<0.and.abs(lbface(1,IFACE))/=nb) then
!        elseif(lbface(1,IFACE)<0) then
          lbface(1,IFACE)=nb
!          do iv=1,4
!          if(lvf(iv,IS,k)==8) cycle
!          lvcell(lvf(iv,IS,k),IC)=ip(lvcell(lvf(iv,IS,k),IC))
!          enddo

          do iv=1,4
          lvface(iv,IFACE)=lvcell(lvf(iv,IS,k),IC)
          enddo
        endif
  401   continue
  400   continue
      endif
      enddo
      nvrtx=nvrtxnw
      endif

!
      if(icon_cv==vertex_cen) then
        NCV=NVRTX
      endif

      DO 50 IV=1,NCV
      LVRTCV(IV)=IV
  50  CONTINUE

!
!
!
      if(ncelb>mcell) then
        write(ifle,*) 'ERR: Increasing mcell'
        stop
      endif
!
      if(nface>mface) then
        write(ifle,*) 'ERR: Increasing mface'
        stop
      endif
!
      deallocate(lvrtx,nclv,ip,ityp)
      write(ifll,'(a,I10)') 'MSG: TOTAL FACE=> NFACE= ',NFACE
!
      do 410 IS=1,nface
      do 410 iv=1,4
      LCVFAC(iv,IS)=lvface(iv,IS)
  410 continue
!
      call list_facerr(ifle,mcell,mface,nface,ncelb,ierr1)
      if( ierr1.ne.0 ) goto 9999
!------------------
! --- undefined BC
!------------------
      allocate(ip(mface),nclv(0:mvrtx),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating ip(:)_2 in list_face'
!
      if(ivbcnd(nbcnd).eq.0) then
!!!!!!!        NBOUND=NBOUND+1
        ip(:)=0
        do nb=1,nbcnd-1
          nclv(:)=0
          nclv(0)=1
          do i=ivbcnd(nb-1)+1,ivbcnd(nb)
             iv=lvbcnd(i)
             nclv(iv)=1
          enddo
          do 500 IS=NIFACE,nface
            do ivv=1,4
              iv=lvface(ivv,IS)
              if(nclv(iv).eq.0) goto 500
            enddo
            ip(IS)=nb
 500      continue
        enddo
        allocate(ityp(ivbcnd(nbcnd-1)),stat=ierr1)
        do nb=1,nbcnd-1
          do i=ivbcnd(nb-1)+1,ivbcnd(nb)
          ityp(i)=lvbcnd(i)
          enddo
        enddo
        deallocate(lvbcnd)
!
        nclv(:)=0
        do IS=NIFACE,nface
          if(ip(IS).eq.0) then
            do ivv=1,4
            iv=lvface(ivv,IS)
            nclv(iv)=1
            enddo
          endif
        enddo
!
        ivv=ivbcnd(nbcnd-1)
        do iv=1,nvrtx
        if(nclv(iv).eq.1) then
          ivv=ivv+1
        endif
        enddo
        ivbcnd(nbcnd)=ivv
!
        allocate(lvbcnd(ivbcnd(nbcnd)),stat=ierr1)
        DO nb=1,nbcnd-1
          do i=ivbcnd(nb-1)+1,ivbcnd(nb)
            lvbcnd(i)=ityp(i)
          enddo
        ENDDO
        ivv=ivbcnd(nbcnd-1)
        do iv=1,nvrtx
          if(nclv(iv).eq.1) then
            ivv=ivv+1
            lvbcnd(ivv)=iv
          endif
        enddo
        deallocate(ityp)
      endif
!
      deallocate(ip,nclv)
!
      return
!
 9001 continue
      write(ifle,*) ' ### error : data error'
      write(ifle,*) 'no. of vertices =',kclv
      write(ifle,*) 'it must be 4,5,6 or 8'
      write(ifle,*) 'cell no. =',IC
      goto 9999
 9002 continue
      write(ifle,'(a)') ' ### error : allocation failed'
      goto 9999
 9003 continue
      write(ifle,*) ' ### error : data error'
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
      end subroutine list_face
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE LIST_EDGE
     & (MCELL,MFACE,MEDGE,mvrtx,NCELL,NVRTX,NFACE,NEDGE,NCELB,NCV,
     &  LVCELL,LFCELL,LVFACE,LCFACE,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,IERROR)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE MODULE_IO,ONLY : IFLE
      use module_partitioner,only : nclv=>WK1
      use module_partitioner,only : ityp=>WK2
      use module_partitioner,only : ip  =>WK3
      use module_partitioner,only : lvrtx=>IW1
      use module_model,only       : vertex_cen,cell_cen,icon_cv
      use module_parameter,  only : NIFACE      
!
!     1. MAKE UP LIST OF CONNECTIVITY EDGE-TO-FACE AND EDGE-TO-VERTEX
!
      implicit none
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)  :: MCELL,MFACE,MEDGE,mvrtx,NCV
      INTEGER,INTENT(IN)  :: NCELL,NVRTX,NFACE,NCELB
      INTEGER,INTENT(OUT) :: NEDGE
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
! --- [LOCAL ENTITIES]
!
      INTEGER,PARAMETER :: NF2TYP(4)=(/0,0,1,2/)
      INTEGER,PARAMETER :: NEDGEL(2)=(/3,4/)
! --- DEFINE LOCAL EDGE VERTEX ORDER
      INTEGER,PARAMETER :: LVEDGL(2,4,2)=RESHAPE( SOURCE=
     &            (/1,2, 2,3, 3,1, 0,0,
     &              1,2, 2,3, 3,4, 4,1/),SHAPE=(/2,4,2/))
!
      INTEGER :: LVF(2,4,2)
      INTEGER :: LVRT(2)
      INTEGER :: LF1,IE,IERR1=0,IERR2
      INTEGER :: KCLV,IFMAX,MATCH,NCH,NEDGEX,NCLBX
      INTEGER :: I,J,K,L,M,N
      INTEGER :: IC,IS,IV,IVA,IVB,IC1,IC2,ILCV,ICV,ICTP
      INTEGER :: ILS,ILV,ILE,IVV,IEE
      INTEGER :: ILE1,ILE2,IS1,IS2,NE1,NE2,NE,IEGVT1,IEGVT2
      LOGICAL :: INNER,IEDGE
      INTEGER :: NF
!
!-< 1. SET UP SOME ARRAYS >-
!
      if(icon_cv==cell_cen) then
      endif

      allocate(nclv(4*NFACE),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating nclv(:) in LIST_EDGE'
      allocate(ityp(NFACE),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ityp(:) in LIST_EDGE'
      allocate(ip(0:NCV),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ip(:) in LIST_EDGE'
!
!--< 1.1 MAKE UP "LVF" >--
!
      IERROR=0
!
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
      DO 120 IS=1,NFACE
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
!-< 2. SET UP CELL-FACE LIST >-
!
!--< 2.1 MAKE UP LIST FOR FACES CONNECTED WITH EACH VERTEX >--
!
      CALL LIST_ONVRTX(4,MFACE,NCV,NFACE,LCVFAC,IP,NCLV,IFMAX)
!
      ALLOCATE(LVRTX(4,2*IFMAX),STAT=IERR1)
      IF(IERR1.NE.0) GOTO 9002
!
!--< 2.2 SEARCH INTERNAL FACES >-
!
      NEDGE=0
!
      DO 200 ICV=1,NCV
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
      IS1=LVRTX(4,NE1)    !338264
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
        if(NEDGE>MEDGE) stop 'MSG: Increase MEDGE'
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
!
! --- check edge on face
!
      NCH=0
      do IS=1,nface
      do i=1,LEFACE(5,IS)
        if(LEFACE(i,IS).eq.0) then
          NCH=NCH+1
          write(IFLE,*) ' ERR-MSG: //local-edge no.=',i,
     &    ' //gl-face no=',IS,
     &    ' //No edge specified : ',LEFACE(i,IS),
     &    ' //Face vertex no.= ',LVFACE(:,IS),' //two cell no.=',
     &    LCFACE(1,IS),LCFACE(2,IS)
        endif
      enddo
      enddo
      if(NCH.gt.0) then
        stop 'ERR: some edge NOT specified on face at LIST_EDGE'
      endif
!
! --- check vertex on edgw
!
      NCH=0
      do IE=1,NEDGE
        if(LVEDGE(1,IE).eq.0.or.LVEDGE(2,IE).eq.0) then
          NCH=NCH+1
        endif
      enddo
      if(NCH.gt.0) then
        stop 'NO vertex on edge no'
      endif
!
      RETURN
!
 9001 CONTINUE
      WRITE(IFLE,*) ' ### ERROR : DATA ERROR'
      WRITE(IFLE,*) 'NO. OF VERTICES =',KCLV
      WRITE(IFLE,*) 'IT MUST BE 3 OR 4'
      WRITE(IFLE,*) 'FACE NO. =',IS
      GOTO 9999
 9002 CONTINUE
      WRITE(IFLE,*) ' ### ERROR : ALLOCATION FAILED'
      GOTO 9999
 9999 CONTINUE
      WRITE(IFLE,*) '(LIST_EDGE)'
      IERROR=1
!
      END SUBROUTINE LIST_EDGE
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facebc
     & (mface,nface,mvrtx,nvrtx,ncell,mcell,
     &  lvface,lbface,lcface,lacell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only          : ifle,ifll
      use module_boundary,only    : nbcnd,ivbcnd,lvbcnd,kdbcnd,kdprdc,
     &                              kdintr,kdbuff,kdsld,boundName,
     &                              kdtchi,idis,boundIDMap,ivpair,
     &                              kdshutr,kdpors,kdovst,
     &                              thinw
      use module_partitioner,only : lbc=>WK1
!???      use module_partitioner,only : lbfacx=>WK2
      use module_partitioner,only : lbc1=>WK3
      use module_parameter,  only : NIFACE
!
! 1. Make up list to link boundary conditions
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mface,nface,nvrtx,mvrtx,ncell,mcell
      integer,intent(in)    :: lvface(4,mface)
      integer,intent(in)    :: lcface(2,mface)
      integer,intent(inout) :: lbface(2,mface)
      integer,intent(in)    :: lacell(  mcell)
!
! --- [local entities]
!
      integer :: i,m,n,nb,mb,id,idum1,idum2
      integer :: IC,IS,IV,IE,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEEs
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC,ierr1=0,kd
      logical :: logflg=.false.
!
      allocate(lbc(0:nvrtx),stat=ierr1)
      if(ierr1/=0) stop 'stop at allocating lbc(:) in list_facebc'
      allocate(lbc1(0:nvrtx),stat=ierr1)
      if(ierr1/=0) stop 'stop at allocating lbc1(:) in list_facebc'
!
!      allocate(lbfacx(mface),stat=ierr1) 
!      if(ierr1/=0) 
!     &    stop 'stop at allocating lbfacx(:) in list_facebc' 
!
! --- 
!
      lbc(1:nvrtx)=0
      lbc(0)=1
      lbc1(1:nvrtx)=0
      lbc1(0)=-1
!
!      do 110 IS=1,nface
!      nb=lbface(1,IS)
!      lbfacx(IS)=nb
!  110 continue
!
      do 200 nb=1,nbcnd 
      kd=kdbcnd(0,nb)
!      if(thinw(nb)==1) cycle
      logflg=(kd==kdintr.and.idis(nb)==0).or. 
     &       (kd==kdsld .and.idis(nb)==0).or.
     &       (kd==kdbuff.and.idis(nb)==0).or.
     &       (kd==kdpors.and.idis(nb)==0).or.
     &       (kd==kdprdc.and.idis(nb)==0).or.
     &       (kd==kdshutr.and.idis(nb)==0)
      if(.NOT.logflg) then
! --- single BC 
          do 301 id=ivbcnd(nb-1)+1,ivbcnd(nb)
          iv=lvbcnd(id)
          lbc(iv)=1
  301     continue
          IFL=0 
          do 310 IS=1,nface 
!          if(lcface(2,IS).le.ncell.and.thinw(nb)/=1) goto 310 
          if(lcface(2,IS).le.ncell) goto 310 
          do 311 i=1,4
          if(lbc(lvface(i,IS)).lt.1) goto 310
  311     continue
          if(lbface(1,IS)/=0.and.lbface(1,IS)/=nb) then 
!            if(thinw(nb)==1) then
!!              lbface(1,IS)=nb
!              goto 310
!            else
              write(*,*)
              write(*,'(1X,a,1x,2I4)') 'MSG: lbface(1,IS),nb',
     &        lbface(1,IS),nb
              write(*,'(1X,a)') 'MSG: 2 Layer mesh is necessary, or: '
              write(*,'(1X,a)') 
     &    'MSG: Different BC name is necessary for Different material'
              write(*,'(1X,a,I8,4I8)') 'MSG:IS=, lvface(1:4,IS)',
     &        IS,lvface(:,IS)
              write(*,'(1X,a,I8)')'MSG: IC1=',lcface(1,IS)
!              stop 'STOP at list_facebc'
            endif
!          endif
!
          lbface(1,IS)=nb 
          IFL=IFL+1
  310     continue
          do 302 id=ivbcnd(nb-1)+1,ivbcnd(nb)
            iv=lvbcnd(id)
            lbc(iv)=0
  302     continue
          write(ifll,'(1X,a,1x,2I8)') 'MSG: nb=, BC-face=',nb,IFL
!          
          if(idis(nb)==1.and.kd/=kdtchi) then
! --- Paired Discontinuous kdintr,kdsld,kdprdc,BC
            do id=ivpair(nb)+1,ivbcnd(nb)
              iv=lvbcnd(id)
              lbc1(iv)=-1
            enddo
!
            do IS=1,nface
              if(lcface(2,IS).le.ncell) cycle
              do i=1,4
                if(lbc1(lvface(i,IS))/=-1) goto 341
              enddo
              lbface(1,IS)=-nb
 341          continue
            enddo 
!
            do id=ivpair(nb)+1,ivbcnd(nb)
              iv=lvbcnd(id)
              lbc1(iv)=0
            enddo
          endif
!      elseif(kd==kdprdc.and.idis(nb)==0) then
!        do 220 IS=1,nface
!        if(abs(lbfacx(IS))==nb) lbface(1,IS)=lbfacx(IS)
!  220   continue
      endif
  200 continue
!      do IS=1,nface
!        IC2=lcface(2,IS)
!        nb=lbface(1,IS)
!        if(IC2.gt.ncell) then
!          if(nb==0) then
!          endif
!        endif
!      enddo
!
      deallocate(lbc,lbc1)   !,lbfacx)
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
      end subroutine list_facebc
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facechk
     & (mface,mcell,nface,ncell,mvrtx,nvrtx,
     &  lvface,lcface,lbface,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only  : ifle
      use module_boundary,only  : nobcnd,thinw
!
! 1. Check validity of boundary face
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mface,mcell,nface,ncell,nvrtx,mvrtx
      integer,intent(in)    :: lvface(4,mface)
      integer,intent(in)    :: lcface(2,mface)
      integer,intent(in)    :: lbface(2,mface)
      integer,intent(inout) :: lacell(  mcell)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer :: i,IS,ie,IC1,IC2,nb,IBC
!
! --- 
!
      IBC=0
      do 100 IS=1,nface
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      nb=lbface(1,IS)
      if(IC2.le.IC1) then
        write(ifle,*) 'ERR : IC1 > IC2',IC1,IC2
      endif
      if(IC2.gt.ncell) then
! --- Check undefined face:
        if(nb.eq.0 ) then
          IBC=IBC+1
          ie=3+min(1,lvface(4,IS))
          write(ifle,*) 'ERR: data error'
          write(ifle,*) 'no boundary condition is specified ',
     &                  'for boundary face'
          write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
          write(ifle,*) 'cell no. who has the face =',IC1
          write(ifle,'(1X,a,I6)') 'MSG: nb= ',nb,IS
          goto 9999
        endif
! --- Dumy cell IC2 must be same material no. woth IC1:
        lacell(IC2)=lacell(IC1)
      else
!        if(nb.ne.0.and.thinw(nb)/=1) goto 9002
        if(nb.ne.0) goto 9002
!        if(nb.ne.0) lbface(1,IS)=0 !???
        if(lacell(IC1).gt.0.and.lacell(IC2).lt.0) goto 9003
        if(lacell(IC1).lt.0.and.lacell(IC2).gt.0) goto 9003
      endif
  100 continue
!      stop
!
      return
!
 9002 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) ' ### error : data error'
      write(ifle,*) 'boundary condition can not be specified ',
     &              'for internal face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',lcface(1,IS)
      write(ifle,*) 'boundary condition no. =',lbface(1,IS)
!onishi
      write(ifle,*) 'face no. =',IS,lbface(1,IS)
      goto 9999
 9003 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) ' ### error : data error'
      write(ifle,*) 'boundary condtion must be specified ',
     &              'for interface of fluid & solid'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',IC1,IC2
      goto 9999
 9999 continue
      write(ifle,'(a)') '(list_facechk)'
      ierror=1
      end subroutine list_facechk
!

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facein
     & (mcell,mface,medge,mvrtx,ncell,nface,ncelb,nvrtx,NCV,
     &  lacell,lfcell,lvface,lcface,lbface,lvcell,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only          : ifle,ifll,cntlnam
      use module_boundary,only    : kdbcnd,kdintr,kdprdc,nbcnd,kdbcnd,
     &                              idis
      use module_partitioner,only : lvrtx=>WK1
      use module_partitioner,only : LCVVRT=>WK2
      use module_partitioner,only : lbfacx=>WK3
      use module_partitioner,only : lc2=> WK4
!
! 1.  Make up list of connectivity for interface boundary
!     2D cell for:
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
      integer,intent(in)    :: lvcell(8,mcell)
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
! --- [local entities]
!
      integer :: i,j,ie,nfacx,nfacy,ierr1=0,ierr2,NnewCV,icmax
      integer :: nb,m1,m2,n1,n2,n1d,n2d,kd
      integer :: IC,IS,IV,IVA,IVB,IC1,IC2,ILC,ICV,ILCV
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL,IIIS
      integer :: IBS,IPS,IDC,IIS,IDC1,IDC2,lvtemp(4)=0
      integer :: MMAT(-100:100),MATC_NO(200),NNMAT,IIMAT,IIIMAT
      logical :: logitr=.false.
!
      ierror=0
!
      return
      do nb=1,nbcnd
      if(kdbcnd(0,nb)==kdintr.and.idis(nb)==0) then
        logitr=.true.
      endif
      enddo
      if(.not.logitr) return
!
      allocate(lvrtx(mvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating lvrtx(:) in list_facein'
      allocate(LCVVRT(mvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating LCVVRT(:) in list_facein'
      allocate(lbfacx(mface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating lbfacx(:) in list_facein'      
      allocate(lc2(nface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating lc2(:) in list_facein'  
!
      ierror=0
!
      lbfacx=0
      LCVVRT=0
!
! --- 
!
      nfacx=nface
!
      do 100 IS=1,nfacx
      nb=lbface(1,IS)
      if(nb.lt.1) goto 100
      kd=kdbcnd(0,nb)
      if(kd/=kdintr) goto 100
      IIS=lbface(2,IS)
      ie=3+min(1,lvface(4,IS))
! --- if periodic BC
      if(IIS.gt.0) goto 101
! --- 
! / internal face /
!
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      if(lacell(IC1).eq.lacell(IC2)) goto 9002
! --- if cell IC2 has been defined as BC:
      if(IC2.gt.ncell) goto 9001
      nface=nface+1
      ncelb=ncelb+2
!
      IDC1=ncelb-1
      IDC2=ncelb
      IIS =min(nface,mface)
!
! --- Redifine lcface(1:2,IS)
!
      if(lacell(IC1).gt.lacell(IC2)) then
        lcface(1,IS)=IC1
        lcface(2,IS)=IDC1
        lbface(1,IS)=nb
        lbface(2,IS)=IIS
!
        lcface(1,IIS)=IC2
        lcface(2,IIS)=IDC2
        lbface(1,IIS)=-nb
        lbface(2,IIS)=IS
!
        lvface(4,IIS)=0
        do 110 i=1,ie
        lvface(i,IIS)=lvface(ie+1-i,IS)
  110   continue
!
        ILC=0
        do 120 i=1,lfcell(7,IC2)
        if(lfcell(i,IC2).eq.IS) ILC=i
  120   continue
!
        if( ILC.lt.1 ) then
          write(ifle,'(a)') 'ERR: program error -1- (list_facein)'
          stop 999
        endif
        lfcell(ILC,IC2)=IIS
      elseif(lacell(IC1).lt.lacell(IC2)) then
        lcface(1,IS)=IC2
        lcface(2,IS)=IDC2
        lbface(1,IS)=nb
        lbface(2,IS)=IIS
!
        lcface(1,IIS)=IC1
        lcface(2,IIS)=IDC1
        lbface(1,IIS)=-nb
        lbface(2,IIS)=IS
!
        lvface(4,IIS)=0
        lvtemp(:)=0
        do 210 i=1,ie
        lvtemp(i)=lvface(ie+1-i,IS)
        lvface(i,IIS)=lvface(i,IS)
  210   continue
        do i=1,ie
        lvface(i,IS)=lvtemp(i)
        enddo
!
        ILC=0
        do 220 i=1,lfcell(7,IC1)
        if(lfcell(i,IC1).eq.IS) ILC=i
  220   continue
!
        if( ILC.lt.1 ) then
          write(ifle,'(a)') 'ERR: program error -2- (list_facein)'
          stop 998
        endif
        lfcell(ILC,IC1)=IIS
      else
        goto 9002
      endif
!
      goto 100
!
! / face on periodic boundary /
!
  101 continue
!
      IC1=lcface(1,IS)
      IC2=lcface(1,IIS)
      lbfacx(IIS)=nb
! --- IC1 and IC2 are same material:
      if(lacell(IC1).eq.lacell(IC2) ) goto 9002
!
  100 continue
!
      DO IS=1,NFACE
      DO  ILCV=1,4
      LCVFAC(ILCV,IS)=lvface(ILCV,IS)
      enddo
      enddo
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
!---------------------------------------------------------------------
! --- new CVs
!---------------------------------------------------------------------
!
      lvrtx=0
      do 300 IIS=nfacx+1,nface
      do 310 ILV=1,4
      IV=lvface(ILV,IIS)
      if(IV.gt.0) then
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
        ICV=NnewCV+NCV
        LVRTCV(ICV)=IV    ! important
        LCVVRT(IV)=ICV
      endif
  320 CONTINUE
!
      NCV=NCV+NnewCV
      if(NCV.gt.mvrtx) then
        write(ifle,'(a,I10)') 'ERR: [mvrtx] Must great than ',NCV
        write(ifle,'(3a,I12)') 
     &  'MSG: Reset [mvrtx] in ',cntlnam,' file to ',NCV+1
        stop
      endif 
      write(ifll,'(a,I12,a)') 'MSG: New CV Number:',NnewCV,
     &              ' Be Created in list_facein'
      write(ifll,'(a,I12)') 'MSG: Total CV Number:',NCV
      write(ifll,'(a,I12,a)') 'MSG: New Face Number:',nface-nfacx,
     &              ' Be Created in list_facein'
!
! --- 
!
      do IS=1,NFACE
        IC1=lcface(1,IS)
        IC2=lcface(2,IS)
        if(IC2.lt.ncell) then
          if(lacell(IC1).ne.lacell(IC2)) then
            write(ifle,'(a)') 
     &      'ERR: List inter-face of diff. mat. error'
            stop ' list_facein'
          endif
        endif
      enddo
!
      MMAT=0
      do IC=1,ncell
        IMAT=lacell(IC)
        MMAT(IMAT)=1
      enddo
      NNMAT=0
      do IMAT=-100,100
      if(MMAT(IMAT).eq.1) then
        NNMAT=NNMAT+1
        MATC_NO(NNMAT)=IMAT
      endif
      enddo
!
      do IIMAT=1,NNMAT
      lc2(:)=0
      IMAT=MATC_NO(IIMAT)
      do IIS=nfacx+1,nface
      IC1=lcface(1,IIS)
      IIIMAT=lacell(IC1)
      if(IIIMAT.eq.IMAT) then
        lc2(IC1)=1
      endif
      enddo
!
      do IC=1,NCELL
        if(lc2(IC).eq.1) then
          do I=1,lfcell(7,IC)
            IIS=lfcell(I,IC)
            IC1=lcface(1,IIS)
            IC2=lcface(2,IIS)
!
            do J=1,lfcell(7,IC1)
              IIIS=lfcell(J,IC1)
              do ILCV=1,4
              IV=lvface(ILCV,IIIS)
              if(IV.gt.0) then
                ICV=LCVVRT(IV)
                if(ICV.gt.0) then
                  LCVFAC(ILCV,IIIS)=ICV
                endif
              endif
              enddo
            enddo
!
            if(IC2.gt.ncell) cycle
!
            do J=1,lfcell(7,IC2)
              IIIS=lfcell(J,IC2)
              do ILCV=1,4
              IV=lvface(ILCV,IIIS)
              if(IV.gt.0) then
                ICV=LCVVRT(IV)
                if(ICV.gt.0) then
                  LCVFAC(ILCV,IIIS)=ICV
                endif
              endif
              enddo
            enddo
          enddo
        endif
      enddo
      enddo
!-----------------------------------------------------------------
      deallocate(lvrtx,LCVVRT,lbfacx,lc2)
      return
!
 9001 continue
      write(ifle,*) ' ### error-1 : data error'
      write(ifle,*) 'the face in interface boundary is other boundary',
     &              'face, not internal face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',IC1
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9002 continue
      write(ifle,*) ' ### error-2 : data error'
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
      write(ifle,*) ' ### error-3 : data error'
      write(ifle,*) 'counterpart face in interface/periodic boundary ',
     &              'has a condition other than periodic boundary'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',lcface(1,IS)
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9999 continue
      write(ifle,*) '(list_facein)'
      ierror=1
!
      end subroutine list_facein
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facenn
     & (mcell,mface,medge,mvrtx,ncell,nface,ncelb,nvrtx,NCV,
     &  lfcell,lvface,lcface,lbface,
     &  LVEDGE,LEFACE,LVRTCV,LCVFAC,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only          : ifle
      use module_boundary,only    : kdbinr,icbinr,lcbinr,kdbuff,kdpors,
     &                              kdbcnd,nbcnd,kdshutr
      use module_partitioner,only : lvrtx=>WK1
      use module_partitioner,only : LCVVRT=>WK2
      use module_parameter,  only : NIFACE
!
! 1. Make up list of connectivity for inner boundary (buffer)
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,mface,ncell,nvrtx,mvrtx
      integer,intent(inout) :: nface,ncelb,medge,NCV
      integer,intent(inout) :: lfcell(7,mcell)
      integer,intent(inout) :: lvface(4,mface)
      integer,intent(inout) :: lcface(2,mface)
      integer,intent(inout) :: lbface(2,mface)
      INTEGER,INTENT(INOUT) :: LEFACE(5,MFACE)
      INTEGER,INTENT(INOUT) :: LVEDGE(2,MEDGE)
      integer,intent(inout) :: LVRTCV(  mvrtx)
      integer,intent(inout) :: LCVFAC(4,mface)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer :: i,j,ie,nfacx,ierr1
      integer :: nb,nb1,nb2,m1,m2,n1,n2,n1d,n2d,kdn
      integer :: IC,IS,IV,IVA,IVB,IC1,IC2,IDC1,IDC2
      integer :: ILS,ILV,ILE,IVV,IEE,NnewCV,ICV,ILCV
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC,IIS
      logical :: logbfl=.false.
!
      ierror=0
      do nb=1,nbcnd
      if(kdbcnd(0,nb)==kdbuff.or.
     &   kdbcnd(0,nb)==kdshutr.or.
     &   kdbcnd(0,nb)==kdpors) then
        logbfl=.true.
      endif
      enddo
      if(.not.logbfl) return
!
      allocate(lvrtx(mvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating lvrtx(:) in list_facenn'
      allocate(LCVVRT(mvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating LCVVRT(:) in list_facenn'
!
      LCVVRT=0
!
! --- 
!
      lvrtx=0
      DO IS=1,nface
      nb=lbface(1,IS)
      if(nb.lt.1) cycle
      if(kdbcnd(0,nb).ne.kdbuff.and.kdbcnd(0,nb)/=kdshutr) cycle
      IC2=lcface(2,IS)
      if(IC2.le.ncell) then
        ie=3+min(1,lvface(4,IS))
        do i=1,ie
        iv=lvface(i,IS)
        lvrtx(iv)=1
        enddo
      endif
      ENDDO
!
      nfacx=nface
!
      do 100 IS=1,nfacx
      nb=lbface(1,IS)
      if(nb.lt.1) goto 100
      if(kdbcnd(0,nb)/=kdbuff.and.
     &   kdbcnd(0,nb)/=kdshutr.and.
     &   kdbcnd(0,nb)/=kdpors) goto 100
!
      ie=3+min(1,lvface(4,IS))
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
!
      if(IC2.le.ncell) then
        nface=nface+1
        ncelb=ncelb+2
!
        IDC1=ncelb-1
        IDC2=ncelb
        IIS =min(nface,mface)
!
        lcface(1,IS) =IC1
        lcface(2,IS) =IDC1
        lbface(1,IS) =nb
!
        lcface(1,IIS)=IC2
        lcface(2,IIS)=IDC2
        lbface(1,IIS)=nb
!
        lvface(4,IIS)=0
        do 120 i=1,ie
        lvface(i,IIS)=lvface(ie+1-i,IS)
  120   continue
!
        j=0
        do 130 i=1,lfcell(7,IC2)
        if(lfcell(i,IC2).eq.IS) j=i
  130   continue
!
        if(j.lt.1) then
          write(*,*) ' ### program error -1- (list_facenn)'
          stop 999
        endif
        lfcell(j,IC2)=IIS
      else
! --- NOT necessary for wall-buffle      
      endif
  100 continue
!
      call list_facerr(ifle,mcell,mface,nface,ncelb,ierr1)
      if( ierr1.ne.0 ) goto 9999
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
      write(ifle,*) ' ### error : data error'
      write(ifle,*) 'the face in inner boundary is boundary ',
     &              'face, not internal face'
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'cell no. who has the face =',IC1
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9002 continue
      write(ifle,*) ' ### error : data error'
      write(ifle,*) 'the two cells on inner boundary is located ',
     &              'at the same side of the boundary'
      write(ifle,*) 'cell no. =',IC1,' and',IC2
      write(ifle,*) 'vertices of the face =',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'boundary condition no. =',nb
      goto 9999
 9999 continue
      write(ifle,*) '(list_facenn)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printi('lfcell/list_facenn',lfcell,7,mcell)
      call printi('lvface/list_facenn',lvface,4,mface)
      call printi('lcface/list_facenn',lcface,2,mface)
      call printi('lbface/list_facenn',lbface,2,mface)
      end subroutine debug
!
      end subroutine list_facenn
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facepr
     & (mface,nface,mvrtx,nvrtx,ncell,mcell,mssfbc,
     &  IBCCYC,IBCSLD,IBCINT,
     &  lvface,lcface,lbface,lacell,listpr,cord,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only          : ifle
      use module_boundary,only    : nbcnd,ivbcnd,lvbcnd,kdbcnd,kdprdc,
     &                              kdsld,kdintr,idis,kdbuff,kdshutr,
     &                              kdpors
      use module_partitioner,only : listva=>WK1
      use module_partitioner,only : listvb=>WK2
      use module_partitioner,only : listsa=>WK3
      use module_partitioner,only : listsb=>WK4
!
! --- nbcnd: BC domain no.
!
! 1. Make up list of connectivity for periodic boundary
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mface,nface,nvrtx,ncell,mvrtx,mssfbc,
     &                       mcell
      integer,intent(in)  :: lvface(4,mface)
      integer,intent(in)  :: lcface(2,mface)
      integer,intent(in)  :: lacell(  mcell)
      integer,intent(out) :: lbface(2,mface)
      integer,intent(out) :: listpr(0:3,mvrtx)
      integer,intent(out) :: IBCCYC,IBCSLD,IBCINT
      real*8 ,intent(in)  :: cord(  3,mvrtx)
!
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer :: list_d(4,91)
      integer :: i,j,k,m,n,nb,kd,nbuf,nbufa,nbufb
      integer :: j1,j2,m1,m2,n1,n2,ierr1
      integer :: nvfc,icmax,match,ie,nf,mch,ip0,ipe
      integer :: IC,IS,IS1,IS2,IV,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE,IV1,IV2,i1,i2,i3,i4
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC
!
      IBCCYC=0
      IBCSLD=0
      ierror=0
      do 20 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if((kd==kdprdc.and.idis(nb)==0).or.
     &   (kd==kdsld .and.idis(nb)==0).or.
     &   (kd==kdbuff.and.idis(nb)==0).or.
     &   (kd==kdpors.and.idis(nb)==0).or.
     &   (kd==kdintr.and.idis(nb)==0).or.
     &   (kd==kdshutr.and.idis(nb)==0)
     &    ) 
     &    goto 21
   20 continue
      return
   21 continue
!
      allocate(listva(0:mvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating listva(:) in list_facepr'
      allocate(listvb(0:mvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating listvb(:) in list_facepr'
      allocate(listsa(mface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating listsa(:) in list_facepr' 
      allocate(listsb(mface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating listsb(:) in list_facepr' 
!
! --- 
!
      lbface(:,:)=0
!
      listva(0)=1
      listvb(0)=2
!
!----------------------------------------------------
! --- Periodic Face listing (Cell Center)
!----------------------------------------------------
!      
!
      listpr=0
      do 1000 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if(kd==kdprdc.and.idis(nb)==0) then
!------------------------
! --- list peroidic BC 
!------------------------
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
        if(listpr(0,IV1)>3) stop 'list-facepr-1'
        listpr(0,IV2)=listpr(0,IV2)+1
        if(listpr(0,IV2)>3) stop 'list-facepr-2'
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
        do 1200 IS=1,nface
        if(lcface(2,IS).le.ncell) goto 1200
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
        if(nbufb/=nbufa) then
          write(ifle,'(a,2I10)')
     &    'ERR: periodic BC face No. of A & B not equal= '
     &    ,nbufa,nbufb
          stop 'list_facepr'
        endif
!..................................................
!-< 3. matching procedure >- 
!..................................................
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
          write(ifle,'(a,I8)') 
     &    'MSG: NO. OF PERIDIC FACE PAIRS ', MCH
        else
          write(ifle,'(a,I8)') 
     &    'ERR : NO. OF PERIDIC PAIRS MUST BE ', NBUFA
          stop
        endif
        IBCCYC=IBCCYC+mch
!
      elseif(kd==kdprdc.and.idis(nb)==1) then
         
!
        
      elseif((kd==kdsld).and.idis(nb)==0) then
!------------------------
! --- list splinding BC
!------------------------
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
        if(listpr(0,IV1)>3) stop 'list-facepr-1'
        listpr(0,IV2)=listpr(0,IV2)+1
        if(listpr(0,IV2)>3) stop 'list-facepr-2'
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
        do 2200 IS=1,nface
        if(lcface(2,IS).le.ncell) goto 2200
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
     &    'ERR: Sliding BC face No. of A & B not equal= '
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
     &    'MSG: NO. OF Sliding FACE PAIRS ', MCH
        else
          write(ifle,'(a,I8)') 
     &    'ERR: NO. OF Sliding PAIRS MUST BE ', NBUFA
          stop
        endif
        IBCSLD=IBCSLD+mch
!
        DO i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listpr(:,IV1)=0
        listpr(:,IV2)=0
        ENDDO
!
      elseif(kd==kdsld.and.idis(nb)==1) then
!
      elseif((kd==kdintr.or.kd==kdbuff.or.kd==kdshutr.or.kd==kdpors)
     &      .and.idis(nb)==0) then 
!------------------------
! --- list interface BC
!------------------------
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
!
        do 3100 i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listva(IV1)=1
        listvb(IV2)=2
        listpr(0,IV1)=listpr(0,IV1)+1
        if(listpr(0,IV1)>3) stop 'list-facepr-1'
        listpr(0,IV2)=listpr(0,IV2)+1
        if(listpr(0,IV2)>3) stop 'list-facepr-2'
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
        do 3200 IS=1,nface
        if(lcface(2,IS).le.ncell) goto 3200
!        if(lbface(2,IS)>0.or.lbface(1,IS)/=0) goto 3210
        if(lbface(1,IS)/=0) goto 3200
        if(lbface(2,IS)>0) goto 3200
        do 3201 i=1,4
        ILV=lvface(i,IS)
        if(listva(ILV).ne.1) goto 3210 
 3201   continue 
        nbufa=nbufa+1
        listsa(nbufa)=IS
        goto 3200
 3210   continue
!
        do 3211 i=1,4
        ILV=lvface(i,IS)
        if(listvb(ILV).ne.2) goto 3220
 3211   continue
        nbufb=nbufb+1
        listsb(nbufb)=IS
 3220   continue
!
 3200   continue
!
        if(nbufb.ne.nbufa) then
          write(ifle,'(a,2I10,a,I4)')
     &    'ERR: Interface/Buffle BC face No. of A & B not equal= '
     &    ,nbufa,nbufb,' at BC No=',nb
          stop 'list_facepr=> Interface BC'
        endif
!..................................................
!-< 3. matching procedure >-
!..................................................
        mch=0
        do 3300 n1=1,nbufa
        IS1=listsa(n1)
        if(IS1==0) cycle
        do 3400 n2=1,nbufb
        IS2=listsb(n2)
        if(IS2==0) cycle
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
 3304     continue
        endif
 3301   continue
        do 3405 j=1,nbuf
        call list_fmatch(match,lvface(1,IS1),list_d(1,j),-1)
        if( match.ge.1 ) goto 3303
 3405   enddo
 3400   enddo
        goto 9001
 3303   continue
        mch=mch+1
        if(kd==kdintr) then
! --- Only for Solid/Fluid & Solid/Solid
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
          listsb(n2)=0
          listsa(n1)=0
       endif
 3300  continue
!
       if(kd==kdintr) then
       else
         do n1=1,nbufa
         if(listsa(n1)/=0) then
           write(ifle,'(a,I10)') 'ERR: NOT pair IS1 =',listsa(n1)
           IS1=listsa(n1)
           print*,'4 vertex= ',lvface(:,IS1)
           i1=lvface(1,IS1)
           i2=lvface(2,IS1)
           i3=lvface(3,IS1)
           i4=lvface(4,IS1)
           print*,'vertex NO= ',i1,i2,i3,i4
           print*,'vertex cord= ',cord(:,i1)
         endif
         enddo
!
         do n2=1,nbufb
         if(listsb(n2)/=0) then
           print*,'ERR: NOT pair IS2 =',listsb(n2)
           IS2=listsb(n2)
           print*,'4 vertex= ',lvface(:,IS2)
           i1=lvface(1,IS2)
           i2=lvface(2,IS2)
           i3=lvface(3,IS2)
           i4=lvface(4,IS2)
           print*,'vertex NO= ',i1,i2,i3,i4 
           print*,'vertex cord= ',cord(:,i1)
         endif
         enddo
       endif
!
        if(nbufa.eq.mch) then
          write(ifle,'(a,2I8)') 
     &   'MSG: NO. OF Interface FACE PAIRS and nb', MCH,nb
        else
          write(ifle,'(a,I8)') 
     &   'ERR: NO. OF Interface PAIRS MUST BE ', NBUFA
          stop
        endif
        IBCINT=IBCINT+mch
!
        DO i=ip0+1,ip0+ipe
        IV1=lvbcnd(i    )
        IV2=lvbcnd(i+ipe)
        listpr(:,IV1)=0
        listpr(:,IV2)=0
        ENDDO
!
      elseif(kd==kdintr.and.idis(nb)==1) then 
!
      endif
 1000 continue
!
      deallocate(listva,listvb,listsa,listsb)
!
      call debug
!
      return
!
 9001 continue
      ie=3+min(1,lvface(4,IS1))
      write(ifle,'(a)') 'ERR: data error'
      write(ifle,'(2a)') 'MSG: the face in periodic boundary has ',
     &              'no couterpart face'
      write(ifle,*) 'MSG: vertices of the face =',(lvface(i,IS1),i=1,ie)
      write(ifle,*) 'MSG: cell no. who has the face =',lcface(1,IS1)
      write(ifle,*) 'MSG: boundary condition no. =',nb
      goto 9999
 9999 continue
      write(ifle,*) '(list_facepr)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printi('lbface/list_facepr',lbface,2,mface)
      end subroutine debug
!
      end subroutine list_facepr
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facerr(ifle,mcell,mface,nface,ncelb,ierr)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      integer,intent(in)  :: ifle,mcell,mface,nface,ncelb
      integer,intent(out) :: ierr
!
      ierr=0
      if(nface.gt.mface) then
        write(ifle,'(1x,a)') 'MSG: array size is NOT unsufficient'
        write(ifle,'(1x,a,I12)') '[mface] must be >=',nface
        ierr=1
      endif
      if(ncelb.gt.mcell) then
        write(ifle,'(1x,a)') 'MSG: array size is NOT unsufficient'
        write(ifle,'(1x,a,I12)') ' [mcell] must be >=',ncelb
        ierr=1
      endif
      if( ierr.ne.0 ) write(ifle,'(1x,a)') '(list_facerr)'
      end subroutine list_facerr
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_fluid
     & (mcell,mface,nface,ncell,ncelb,mvrtx,nvrtx,
     &  lbface,lcface,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only          : ifle
      use module_material   ,only : strdata
      use module_partitioner,only : iclos=>WK1
      use module_partitioner,only : kdbp=>WK2
!
! 1. Make up list "lacell" for multiple fluid domain
!    separated by solid part
!
! 2. Set flag for closed domain
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,mface,nface,ncell,ncelb,nvrtx,mvrtx
      integer,intent(in)    :: lbface(2,mface)
      integer,intent(in)    :: lcface(2,mface)
      integer,intent(inout) :: lacell(  mcell)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
!      integer :: iclos (  mcell)
!      integer :: kdbp  (  mface)
!
      integer :: jclos (0:100)
      integer :: MMAT(-100:100),MATC_NO(200)
      integer :: MATC_SLD(100),MATC_FLD(100)
      integer :: i,j,m,n,la,ll,nf,nc
      integer :: nxx,ns,ne,nd,iflg
      integer :: ierr1=0,ierr,nflud,ihpc,nfludx,NSLIDX
      integer :: IC,IS,IV,IE,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL,ICFLD,NNMAT
      integer :: IBS,IPS,IDC,IIS,IIMAT
!
      allocate(iclos(mcell),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating iclos(:) in list_fluid'
      allocate(kdbp(mface),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating kdbp(:) in list_fluid'
!
!-< 1. Preliminary set >-
!
      ierror=0
      ierr=0
      ierr1=0
!
!-< 3. Set flag for closed domain >-
!
      call bc_kdbpv1(mface,nface,lbface,kdbp)
!
      iclos(:)=0
      MMAT=0
      do 100 IC=1,ncell
      IMAT=lacell(IC)
      MMAT(IMAT)=1
 100  continue
!
      if(MMAT(0).eq.1) then
        write(ifle,'(a)') 'ERR: Some cells has NOT been defined'
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
      do 310 IS=1,nface
      if( kdbp(IS).eq.1 ) then
! --- Outlet BC:
        IC1=lcface(1,IS)
        IC2=lcface(2,IS)
        iclos(IC1)=1
        iclos(IC2)=1
      endif
  310 continue
!
      do 320 IC=1,ncell
      IMAT=lacell(IC)
      if(IMAT.gt.0) then
        if(iclos(IC).eq.1) then
          jclos(IMAT)=1
        endif
      endif
  320 continue
!
      ihpc=0
      call strdata(ifle,NFLUDX,MATC_FLD,MATC_NO,jclos,ihpc,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      deallocate(iclos,kdbp)
!
      call debug
      return
!
 9999 continue
      write(ifle,*) '(list_fluid)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
      subroutine debug
      use module_debug,only : idebug
      if( idebug(2).eq.0 ) return
      call printj('lacell/list_fluid',lacell,1,mcell)
      end subroutine debug
!
      end subroutine list_fluid
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_fmatch(match,lvrtx1,lvrtx2,isw)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Two face is or not same one
!     isw=0: check 3p or 4p separatly 
!     isw<0: converse sequence check (A:1,2,3,4=>B:4,3,2,1)
!     isw>0: positive sequence check (A:1,2,3,4=>B:1,2,3,4)
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(out) :: match
      integer,intent(in)  :: lvrtx1(4),lvrtx2(4),isw
!
!
! --- [local entities]
!
      integer :: lvrtx(8)
      integer :: i,j,ie
!
!
!
      match=0
      ie=min(1,lvrtx1(4))
! --- compare 4p-face and 3p-face:
      if( ie.ne.min(1,lvrtx2(4)) ) return
! --- ie=0: for 3p-face /// ie=1: for 4p-face
!
      ie=3+ie
!
      if( isw.eq.0 ) then
        do 100 j=1,ie
        do 101 i=1,ie
        if( lvrtx1(i).eq.lvrtx2(j) ) match=match+1
  101   continue
  100   continue
        match=1/(abs(ie-match)+1)
        return
      endif
!
      do 200 i=1,ie
      lvrtx(i   )=lvrtx2(i)
      lvrtx(i+ie)=lvrtx2(i)
  200 continue
!
      if(isw.lt.0) then
        do 210 j=ie+ie+1,ie+2,-1
        do 211 i=1,ie
        if(lvrtx1(i).ne.lvrtx(j-i)) goto 210
  211   continue
        match=1
        return
  210   continue
      else
        do 220 j=0,ie-1
        do 221 i=1,ie
        if( lvrtx1(i).ne.lvrtx(j+i) ) goto 220
  221   continue
        match=1
        return
  220   continue
      endif
!
      return
      end subroutine list_fmatch
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_fmatch4(match,lvrtx1,lvrtx2,isw)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Two face is or not same one
!     isw=0: check 3p or 4p separatly 
!     isw<0: converse sequence check (A:1,2,3,4=>B:4,3,2,1)
!     isw>0: positive sequence check (A:1,2,3,4=>B:1,2,3,4)
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: match
      integer,intent(in)  :: lvrtx1(4),lvrtx2(4),isw
!
!
! --- [local entities]
!
      integer :: i,j,ie,imat
!
!
!
      match=0
      imat=0
!      ie=min(1,lvrtx1(4))
! --- compare 4p-face and 3p-face:
!      if( ie.ne.min(1,lvrtx2(4)) ) return
! --- ie=0: for 3p-face /// ie=1: for 4p-face
!
!      ie=3+ie
      ie=4
!
      do 100 j=1,ie
      do 101 i=1,ie
      if(lvrtx1(i).eq.lvrtx2(j)) imat=imat+1
  101 continue
  100 continue
      if(imat.ne.4) match=1
      return
!
      return
      end subroutine list_fmatch4
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_ematch(match,lvrtx1,lvrtx2,isw)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Two edge if/not match
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(out) :: match
      integer,intent(in)  :: lvrtx1(2),lvrtx2(2),isw
!
!
! --- [local entities]
!
      integer :: i,j,ie
!
      ie=2
      match=0
      if(isw.eq.0) then
        do 100 j=1,2
        do 101 i=1,2
        if(lvrtx1(i).eq.lvrtx2(j)) match=match+1
  101   continue
  100   continue
        return
      endif
!
      return
      end subroutine list_ematch
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_onvrtx(mv,mcell,nvrtx,ncell,
     &                       lvcell,ip,nclv,icmax)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_partitioner,only : kclv=>WK5
      use module_partitioner,only : iq  =>WK6
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mv,mcell,nvrtx,ncell
      integer,intent(in)  :: lvcell(mv,mcell)
      integer,intent(out) :: ip    (0:nvrtx)
      integer,intent(out) :: nclv  (mv*ncell)
      integer,intent(out) :: icmax
!
! --- [local entities]
!
      integer :: i,j,n,IC,IV,ILV,ierr1=0
!
      allocate(kclv(ncell) ,stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating kclv(:) in list_onvrtx'
      allocate(iq(nvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating iq(:) in list_onvrtx'
!
!
!-< 1. Count up no. of cells connected with each vertex >--
!
      do 100 IV=1,nvrtx
      iq(IV)=0
  100 continue
!
      do 110 IC=1,ncell
      j=0
      do 111 ILV=1,mv
      if(lvcell(ILV,IC).gt.0) j=ILV
  111 continue
      kclv(IC)=j
      do 112 ILV=1,j
      IV=lvcell(ILV,IC)
      iq(IV)=iq(IV)+1
  112 continue
  110 continue
!
!-< 2. Set pointer of array "nclv" >--
!
      ip(0)=0
      ip(1)=0
      icmax=iq(nvrtx)
      do 200 IV=1,nvrtx-1
      ip(IV+1)=ip(IV)+iq(IV)
      icmax=max(icmax,iq(IV))
  200 continue
!
!-< 3. Set cell no. connected with each vertex >--
!
      do 300 IC=1,ncell
      do 301 ILV=1,kclv(IC)
      IV=lvcell(ILV,IC)
      ip(IV)=ip(IV)+1
      nclv(ip(IV))=IC
  301 continue
  300 continue
!
      deallocate(kclv,iq)
!
      return
      end subroutine list_onvrtx
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_solid(mcell,ncell,ncelb,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only         : ifle
      use module_material,only   : nsold,nosld
      use module_material,only   : nflud,nofld
      use module_hpc,     only   : NPETOT
!      use module_partitioner,only : ichk  =>WK1
!
!-< 1. Check solid no. & renumber lacell >-
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,ncell,ncelb
      integer,intent(inout)    :: lacell(mcell)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer :: i,n,IC,isld,ifld,ierr1=0
!
!      allocate(ichk(mcell),stat=ierr1) 
!      if(ierr1.ne.0) stop 'stop at allocating ichk(:) in list_solid'
!      ichk=0
!
      ierror=0
!
      do 100 IC=1,ncelb

      if(lacell(IC).lt.0) then
        do 101 isld=1,nsold
        if( nosld(isld).eq.lacell(IC) ) then
          if(NPETOT.eq.1) lacell(IC)=-isld
          goto 100
        endif
  101   continue
        write(ifle,*) ' ### error : data error'
        write(ifle,*) 'lack of control data for solid part'
        write(ifle,*) 'solid no. =',lacell(IC)
        goto 9999
      elseif(lacell(IC).gt.0) then
        do 102 ifld=1,nflud
        if(nofld(ifld).eq.lacell(IC)) then
          if(NPETOT.eq.1) lacell(IC)=ifld
          goto 100
        endif
  102   continue
        write(ifle,*) ' ### error : data error'
        write(ifle,*) 'lack of control data for fluid part'
        write(ifle,*) 'fluid no. =',lacell(IC)
        goto 9999

      elseif(lacell(IC).eq.0) then
        write(ifle,*) ' ### error : data error'
        write(ifle,*) 'Element IC= ',IC,' is lack of Material NO.'
        goto 9999
      endif
  100 continue
!
!
!      do IC=1,ncell
!        if(ichk(IC)==0) then
!        write(ifle,'(1X,2a,I6,1X,a,I6)') 'ERR: Cell Material is NOT defined',
!     &  ' Cell Number=',IC,'Defined Material No =',lacell(IC)
!        stop 'ERR: STOP at list_solid'
!        endif
!      enddo
!
!      deallocate(ichk)
!
      return
!
 9999 continue
      write(ifle,*) '(list_solid)'
      ierror=1
!
      end subroutine list_solid
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_edgvtx
     &    (mv,medge,NCVFAC,NALLCV,lvcell,IEMAX)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_partitioner,only : iq  =>WK6
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mv,medge,NCVFAC,NALLCV
      integer,intent(in)  :: lvcell(mv,medge)
      integer,intent(out) :: iemax
!
!
! --- [local entities]
!
      integer :: i,j,n,IC,IV,ILV,ierr1=0,ICV,ICF
!
!
      allocate(iq(NALLCV),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating iq(:) in list_edgvtx'
!
!-< 1. Count up no. of cells connected with each vertex >--
!
      do 100 ICV=1,NALLCV
      iq(ICV)=0
  100 continue
!
      do 110 ICF=1,NCVFAC
      j=0
      do 111 ILV=1,mv
      if(lvcell(ILV,ICF).gt.0) j=ILV
  111 continue
      do 112 ILV=1,j
      ICV=lvcell(ILV,ICF)
      iq(ICV)=iq(ICV)+1
  112 continue
  110 continue
!
!-< 2. Maximum IEMAX>--
!
      IEMAX=0
      do 200 ICV=1,NALLCV
      iemax=max(iemax,iq(ICV))
  200 continue
!
      deallocate(iq)
!
      return
      end subroutine list_edgvtx
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_BCDIM
     & (mssfbc,nssfbc,NBCINL,NBCWAL,NBCOUT,NBCCYC,NBCSYM,NBCTCI,
     &  NBCSLD,NBCINT,LBCSSF)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_boundary,only : kdbcnd,kdilet,kvlglw,kvnslp,kvrogf,
     &                           kxnone,kdolet,kdprdc,kdsymm,kdtchi,
     &                           kdsld,kdintr,idis,kvmodl,kdbuff,
     &                           kdshutr,kdpors
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mssfbc,nssfbc
      integer,intent(out)   :: NBCINL,NBCWAL,NBCOUT,NBCCYC,NBCSYM
      integer,intent(out)   :: NBCTCI,NBCSLD,NBCINT
      integer,intent(in)    :: LBCSSF(  mssfbc)
!
! --- [local entities]
!
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
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
      do 1000 IBF=1,nssfbc
      nb=abs(LBCSSF(IBF))
!      if (nb.le.0) goto 1000
      kd=kdbcnd(0,nb)
      kdv=kdbcnd(1,nb)
      if(kd.eq.kdilet) then
        NBCINL=NBCINL+1
      elseif
     & (kdv==kvnslp.or.kdv==kvlglw.or.kd==kxnone.or.kdv==kvrogf.or.
     &   kdv==kvmodl) then
        if(LBCSSF(IBF)>0) then
          NBCWAL=NBCWAL+1
        endif
        if((kd.eq.kdintr.and.idis(nb)==1).and.LBCSSF(IBF)<0) then
          NBCWAL=NBCWAL+1
        endif
      elseif(kd.eq.kdolet) then
        NBCOUT=NBCOUT+1
      elseif(kd.eq.kdprdc) then
        NBCCYC=NBCCYC+1
      elseif(kd.eq.kdsymm) then
        NBCSYM=NBCSYM+1
      elseif(kd.eq.kdtchi) then
        NBCTCI=NBCTCI+1
      elseif(kd.eq.kdsld) then
        NBCSLD=NBCSLD+1
      endif
      if(kd.eq.kdintr.or.kd==kdbuff.or.kd==kdshutr.or.kd==kdpors) then
        NBCINT=NBCINT+1
      endif
 1000 continue
!
      return
!
      end subroutine list_BCDIM
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_output_geom
     & (icpu,rewall,mcell,mvrtx,mface,medge,mssfbc,NBCWAL,
     &  NCV,nvrtx,ncell,NMAT,NCVFAC,
     &  NALLCV,nedge,nssfbc,nssfbc_old,
     &  NCVIN,
     &  LMAT,LVEDGE,LVRTCV,LBCSSF,LCYCSF,locmsh,listbc,NWK,
     &  SFAREA,SFCENT,CVVOLM,CVCENT)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_io,only  : ifle,ifll
      use module_hpc,only : NPETOT,nsldT,nsldbc,nmsldL,nmsldG,nbsldG,
     &                      MAT_BCSLD
!      
      use module_partitioner,only : MAT_CV => WK1
      use module_partitioner,only : MAT_CF => WK2
      use module_partitioner,only : LCF_CF => WK3
      use module_partitioner,only : LCV_CV => WK4
      use module_partitioner,only : LBC_SSF=> WK5
      use module_partitioner,only : LFUTAU => wk6
      use module_partitioner,only : DISALL => WALLDIST
      use module_partitioner,only : LSLDSSF=> IW1
      use module_partitioner,only : LBC_pair=>pMap
!      use module_partitioner,only : NWK =>listpr
      use module_boundary,only    : nbcnd,LBC_INDEX,kdbcnd,kdprdc,
     &                              kdintr,boundName,kdsld,idis,
     &                              nobcnd,kdbuff,kdshutr,kdpors,
     &                              conntBC,kdovst
      use module_material,only    : ical_sld,nsold,nosld,nflud,nofld,
     &                              DONOR_M
      use module_model,only       : ical_mvmsh
      use module_material,only    : KMAT_S,MAT_ALL 
      use module_model,only       : vertex_cen,cell_cen,icon_cv
!
! --- 
!
      use module_rad   , only	  : radflag,radmodel
      use module_radsxf, only     : RDValue,RDIndex,RDN1,RDN2,RDId,
     &	                            NRD,MRD,NRD0,MRD0,MID,StackOpt
      use module_radsxf, only     : VIndex,AreaVol,PROP
      use module_radsxf, only     : IndexToRD,CoefToRD,MatBounEnd
      use module_radsxf, only     : IndexToExpt,RaToExpt,ExptID
      use module_radsxf, only     : NPALL,NDatExpt,NumTran,NumTranIN
!
!--- only for benchmark	 !jiang
      use module_radsxf , only  : NWDAT,WDAT,WPNUM
      use module_radgroup, only : IFlagGroup
!
      implicit none
!
! --- [dummy arguments]
!nssfbc_old
      integer,intent(in)  :: rewall,icpu,nssfbc_old
      integer,intent(in)  :: mcell,mvrtx,mface,medge,mssfbc
      integer,intent(in)  :: NCV,nvrtx,ncell,NCVFAC,NBCWAL
      integer,intent(in)  :: NALLCV,nedge,nssfbc,NCVIN
      integer,intent(out) :: NMAT
      integer,intent(in)  :: LMAT  (  mvrtx)
      integer,intent(inout)  :: LVEDGE(2,medge)
      integer,intent(in)     :: LVRTCV(  mvrtx)
      integer,intent(inout)  :: LBCSSF(mssfbc)
      integer,intent(inout)  :: LCYCSF (mssfbc)
      integer,intent(out)    :: NWK(0:3,mvrtx)
!      
      integer,intent(inout) :: listbc(4,mssfbc)
      integer,intent(inout) :: locmsh (mssfbc)
      real*8 ,intent(in)  :: SFAREA(4,medge)
      real*8 ,intent(in)  :: SFCENT(3,medge)
      real*8 ,intent(in)  :: CVVOLM(  mvrtx)
      real*8 ,intent(in)  :: CVCENT(3,mvrtx)
      
!
! --- [local entities]
!
      integer :: MMAT(-100:100),MAT_NO(200)
      integer :: MAT_CVIDX(0:200)
      integer :: MAT_CFIDX(0:200)
      integer :: MAT_CVEXT(0:200)
      integer :: MAT_BCIDX(0:200)
      integer :: MAT_DCIDX(0:200)
      integer,allocatable :: MAT_BC(:,:) !MAT_BC(nbcnd,2)
      integer,allocatable :: MAT_CHKBC(:,:) !MAT_CHKBC(nbcnd,2)
      integer :: ifl,LIN,LINDEX,LCVIN,kd
      character(len=80),save :: fnam,fnamwl
      integer :: ICV,NCVX,IMAT,IIMAT,KKMAT,ICVL,ierr1=0,IMAT_U
      integer :: ICV2,ICFP,KIMAT2,IIMAT2,ICVA,ICVB
      integer :: NEIBPETOT_S,nn,INTNODTOT,NODTOT,my_rank
      integer :: ios=0,LIDXBC,ICF,KIMAT,ICVLIN,ICFL,nb,IBFL
      integer :: IBF,IBFP,NBSSF,IBFS,IBFE,ICVS,ICVE,I,IDCS,IDCE
      integer :: ICFS,ICFE,ifd
      integer,allocatable :: nsld(:) !nsld(0:nbcnd)
      integer :: nbsld,j,isd
      real*8  :: dx,dy,dz,dum1
      integer :: ICV1,IRD,IDC,IDCP,IDC1
      integer :: NWK1=1,NWK2=2,NWK3=3
!
      allocate(MAT_CV(0:NALLCV),
     &         MAT_CF(0:NCVFAC),
     &         LCV_CV(0:NALLCV),
     &         LBC_SSF(0:nssfbc+1),
     &         LCF_CF(0:NCVFAC),
     &         LBC_pair(0:nbcnd),
     &         MAT_BC(nbcnd,2),
     &         MAT_CHKBC(nbcnd,2),
     &         nsld(0:nbcnd),
     &         stat=ierr1) 
      if(ierr1.ne.0) 
     & stop 'stop at allocating MAT_CV(:) in list_output_geom'

!
      MMAT(:)=2
      do 100 ICV=1,NCV
      IMAT_U=LMAT(ICV)
      MMAT(IMAT_U)=1
 100  continue
!
      NMAT=0
      MAT_NO(:)=0
      do 200 IMAT_U=-100,100
      if(MMAT(IMAT_U)==1) then
        NMAT=NMAT+1
        MAT_NO(NMAT)=IMAT_U
      endif
 200  continue 
!
      if(KMAT_S/=NMAT.and.icpu.eq.0) then
        stop 'ERR: KMAT_S/=NMAT in subroutine list_output_geom'
      endif
      NMAT=KMAT_S
!------------------------
! --- List CV 
!------------------------
      MAT_CVIDX=0
      MAT_DCIDX=0
      MAT_CVEXT=0
      do 310 IIMAT=1,NMAT
      LINDEX=0
      LIN=0 
      LCVIN=0
      IMAT_U=MAT_NO(IIMAT)
!      IMAT_U=MAT_ALL(IIMAT)
      do 300 ICV=1,NALLCV
      KKMAT=LMAT(ICV)
      if(KKMAT.eq.IMAT_U) then
        if(ICV.le.NCVIN) then
          LINDEX=LINDEX+1
        endif
        if(ICV.le.NCV) then
          LCVIN=LCVIN+1
        endif
        if(ICV.gt.NCV) then
          LIN=LIN+1
        endif
      endif
 300  continue
      MAT_CVIDX(IIMAT)=LINDEX
      MAT_CVEXT(IIMAT)=LCVIN
      MAT_DCIDX(IIMAT)=LIN
 310  continue
!
      do 325 IIMAT=1,NMAT
      MAT_CVEXT(IIMAT)=MAT_CVEXT(IIMAT)+MAT_CVEXT(IIMAT-1)
 325  continue
!
      do 320 IIMAT=1,NMAT
      MAT_CVIDX(IIMAT)=MAT_CVIDX(IIMAT)+MAT_CVEXT(IIMAT-1)
 320  continue
!
      MAT_DCIDX(0)=MAT_CVEXT(NMAT)
      do IIMAT=1,NMAT
      MAT_DCIDX(IIMAT)=MAT_DCIDX(IIMAT)+MAT_DCIDX(IIMAT-1)
      enddo
!
      if(MAT_CVEXT(NMAT).ne.NCV) then
        stop 'ERR: NCV'
      endif
!
      if(MAT_DCIDX(NMAT).ne.NALLCV) then
        write(ifle,*) 
     & 'ERR: MAT_CVIDX(NMAT).ne.NALLCV',MAT_CVIDX(NMAT),NALLCV
        stop
      else
        write(ifll,'(a,I4)') 'MSG: TOTAL MATERIAL NUMBER= ',NMAT
        do IIMAT=1,NMAT
        IMAT_U=MAT_NO(IIMAT)
!        IMAT_U=MAT_ALL(IIMAT)
        write(ifll,'(a,I4,a,I10)') 
     & 'MSG: MATERIAL CV NO.= (NCVIN)',IMAT_U,' => ',
     &  MAT_CVIDX(IIMAT)-MAT_CVEXT(IIMAT-1)
        write(ifll,'(a,I4,a,I10)') 
     & 'MSG: MATERIAL CV NO.= (NCV  )',IMAT_U,' => ',
     &  MAT_CVEXT(IIMAT)-MAT_CVEXT(IIMAT-1)
        enddo
      endif
!------------------------
! --- list CV FACE
!------------------------
      MAT_CFIDX=0
      MAT_BCIDX=0
      do 400 IIMAT=1,NMAT
      LINDEX=0
      LIDXBC=0
      IMAT_U=MAT_NO(IIMAT)
!      IMAT_U=MAT_ALL(IIMAT)
      do 410 ICF=1,NCVFAC
      ICV=LVEDGE(1,ICF)
      KIMAT=LMAT(ICV)
      if(KIMAT.eq.IMAT_U) then
        LINDEX=LINDEX+1
        if(ICF.gt.nedge) then
          LIDXBC=LIDXBC+1
        endif
      endif
 410  continue
      MAT_CFIDX(IIMAT)=LINDEX
      MAT_BCIDX(IIMAT)=LIDXBC
 400  continue
!
      do 420 IIMAT=1,NMAT
      MAT_CFIDX(IIMAT)=MAT_CFIDX(IIMAT)+MAT_CFIDX(IIMAT-1)
      MAT_BCIDX(IIMAT)=MAT_BCIDX(IIMAT)+MAT_BCIDX(IIMAT-1)
 420  continue
!
      if(MAT_CFIDX(NMAT).ne.NCVFAC) then
        write(ifle,'(a,2I10)') 
     & ' ### ERR: MAT_CFIDX(NMAT).ne.NCVFAC',
     &      MAT_CFIDX(NMAT),NCVFAC
        stop
      else
        write(ifll,'(a,I4)') 
     & 'MSG: TOTAL MATERIAL NUMBER= ',NMAT
        do IIMAT=1,NMAT
        IMAT_U=MAT_NO(IIMAT)
!        IMAT_U=MAT_ALL(IIMAT)
        write(ifll,'(a,I4,a,I10)') 
     & 'MSG: MATERIAL CV-FACE NO.= ',IMAT_U,' => ',
     &      MAT_CFIDX(IIMAT)
        enddo
      endif
!
      
      if(MAT_BCIDX(NMAT).ne.NSSFBC) then
        write(ifle,'(a,I10,I10)') 
     & ' ### ERR: MAT_BCIDX(NMAT).ne.NSSFBC',
     &      MAT_BCIDX(NMAT),NSSFBC
        stop
      endif
!----------------------------
! --- List LCV_CV and MAT_CV
!----------------------------
      LCV_CV=0
      MAT_CV=0
      do 550 IIMAT=1,NMAT
      ICVL=  MAT_CVIDX(IIMAT-1)
      ICVLIN=MAT_CVEXT(IIMAT-1)
      LIN=   MAT_DCIDX(IIMAT-1)
      IMAT_U=MAT_NO(IIMAT)
!      IMAT_U=MAT_ALL(IIMAT)
      do 560 ICV=1,NALLCV
      KKMAT=LMAT(ICV)
      if(KKMAT.eq.IMAT_U) then
        if(ICV.le.NCV) then
          ICVLIN=ICVLIN+1
          MAT_CV(ICVLIN)=ICV
          LCV_CV(ICV)=ICVLIN
        endif
        if(ICV.gt.NCV) then
          LIN=LIN+1
          MAT_CV(LIN)=ICV
          LCV_CV(ICV)=LIN
        endif
      endif
 560  continue
 550  continue
!----------------------------
! --- Test & check face
!----------------------------
      MAT_CF=0
      if(icon_cv==vertex_cen) then
        do IBF=1,nssfbc
         ICF=IBF+NEDGE
         NBSSF=LBCSSF(IBF)
         ICV=LVEDGE(1,ICF)
         ICV2=LVEDGE(2,ICF)
         MAT_CF(ICF)=1
         if(ICV2<NVRTX) stop 2221
         if(ICV>NVRTX) stop 2222
        enddo

        do ICF=1,NCVFAC
         if(MAT_CF(ICF)==0) then
           ICV=LVEDGE(1,ICF)
           ICV2=LVEDGE(2,ICF)
           if(ICV2>NVRTX) then
              stop 2223
           endif
           if(ICV>NVRTX) then
              stop 2224
           endif
         endif
        enddo
      else
        do IBF=1,nssfbc
         ICF=IBF+NEDGE
         NBSSF=LBCSSF(IBF)
         ICV=LVEDGE(1,ICF)
         ICV2=LVEDGE(2,ICF)
         MAT_CF(ICF)=1
         if(ICV2<NCV) stop 2221
         if(ICV>NCV) stop 2222
        enddo

        do ICF=1,NCVFAC
         if(MAT_CF(ICF)==0) then
           ICV=LVEDGE(1,ICF)
           ICV2=LVEDGE(2,ICF)
           if(ICV2>NCV) then
              stop 2223
           endif
           if(ICV>NCV) then
              stop 2224
           endif
         endif
        enddo
      endif
!----------------------------
! --- List MAT_CF and LCF_CF
!----------------------------
      MAT_CF=0
      LCF_CF=0
      do 570 IIMAT=1,NMAT
      IMAT_U=MAT_NO(IIMAT)
!      IMAT_U=MAT_ALL(IIMAT)
      ICFL=MAT_CFIDX(IIMAT-1)
      do 580 ICF=1,NCVFAC
      ICV=LVEDGE(1,ICF)
      KIMAT=LMAT(ICV)
      if(KIMAT.eq.IMAT_U) then
        ICFL=ICFL+1
        MAT_CF(ICFL)=ICF
        LCF_CF(ICF)=ICFL
      endif
 580  continue
 570  continue
!------------------------
! --- list BC FACE
!------------------------
      LBC_INDEX=0
      LBC_pair=0
      LBC_SSF=0
!
      do nb=1,nbcnd
      kd=kdbcnd(0,NB)
      LBC_INDEX(nb)=LBC_INDEX(nb-1)
      IBFL=0
      do IBF=1,nssfbc
      NBSSF=LBCSSF(IBF)
      if(nb.eq.NBSSF) then
        IBFL=IBFL+1
      ENDIF
      enddo
!
      if(kd==kdsld.or.
     &  (kd==kdprdc.and.idis(nb)==1).or.
     &  (kd==kdintr.and.idis(nb)==1)
     &   )then
        if(idis(nb)==1) LBC_pair(nb)=LBC_INDEX(nb)+IBFL
        do IBF=1,nssfbc
        NBSSF=LBCSSF(IBF)
        if(nb.eq.-NBSSF) then
          IBFL=IBFL+1
        endif
        enddo
      endif
      LBC_INDEX(nb)=LBC_INDEX(nb)+IBFL
      enddo
!
      do 120 nb=1,nbcnd
      if(LBC_pair(nb)==0) then
        LBC_pair(nb)=LBC_INDEX(nb)
      endif
 120  continue
!--------------------
! --- List LBC_SSF --
!--------------------
      DO 600 nb=1,nbcnd
      kd=kdbcnd(0,NB)
      IBFL=LBC_INDEX(nb-1)
      do IBF=1,nssfbc
      NBSSF=LBCSSF(IBF)
      if(nb==NBSSF) then
        IBFL=IBFL+1
        ICF=IBF+NEDGE
        ICFL=LCF_CF(ICF)
        LBC_SSF(IBFL)=ICFL
      ENDIF
      enddo
!
      if(kd==kdsld.or.
     &  (kd==kdprdc.and.idis(nb)==1).or.
     &  (kd==kdintr.and.idis(nb)==1)
     &  ) then
        do IBF=1,nssfbc
        NBSSF=LBCSSF(IBF)
        if(nb==-NBSSF) then
          IBFL=IBFL+1
          ICF=IBF+NEDGE
          ICFL=LCF_CF(ICF)
          LBC_SSF(IBFL)=ICFL
        endif
        enddo
      endif
 600  continue
!
!----------------------------------
! --- Check BC if croses materials
!----------------------------------
!
      MAT_CHKBC(:,:)=0
      DO 630 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdsld.and.idis(nb)==0) then
        IBFE=IBFS+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
      elseif(idis(nb)==1) then
        IBFE=LBC_pair(nb)
      endif
      do IBFL=IBFS,IBFE
      ICFL=LBC_SSF(IBFL)
      ICF=MAT_CF(ICFL)
      ICV=LVEDGE(1,ICF)
      KIMAT=LMAT(ICV)
      IBF=ICF-NEDGE
      ICFP=LCYCSF(IBF)+NEDGE
      ICV2=LVEDGE(1,ICFP)
      KIMAT2=LMAT(ICV2)
!
      if(KIMAT.ne.KIMAT2.and.kd==kdprdc) then
        write(ifle,*) 'ERR: NOT Supporting Case: ',
     &                'Peridoic BC belong to different material'
        stop
      endif
!
      if(kd==kdovst) then
         MAT_CHKBC(nb,2)=DONOR_M(KIMAT)
         MAT_CHKBC(nb,1)=KIMAT
         if(DONOR_M(KIMAT)==0) then
           write(*,*) 'MSG: over set BC used'
           write(*,*) 'ERR: NOT defined , IMAT_DONOR=',DONOR_M(KIMAT)
           write(*,*) 
     &        'MSG: Defining Donor Material no by  IMAT_DONOR in &fluid'
           stop 'list_output_geom'
         endif
      elseif(kd.ne.kdintr.and.kd.ne.kdsld.and.kd/=kdbuff
     &               .and.kd/=kdshutr.and.kd/=kdpors) then
        if(MAT_CHKBC(nb,1).eq.0) then
          MAT_CHKBC(nb,1)=KIMAT
          MAT_CHKBC(nb,2)=KIMAT
          write(ifll,'(a,I4,a,I4)') 'MSG: BC no: ',nobcnd(nb),
     &     ' IMAT_U=',KIMAT
        elseif(MAT_CHKBC(nb,1).ne.0.and.MAT_CHKBC(nb,1).ne.KIMAT) then
         write(ifle,'(2a,3I4,2a,1x,a,3I4)') 'ERR: One BC ',
     &  'cannot corss on different Materials ',MAT_CHKBC(nb,1),KIMAT,
     &   nobcnd(nb)
     &  ,'BC name: ',trim(boundName(nb,1)),trim(boundName(nb,2))
     &  ,KIMAT2,ICV2,ICV
         stop 'NOT Support this kind mesh, list_material-1'
        endif
      elseif(kd==kdintr.or.kd==kdbuff.or.kd==kdshutr) then
        if(KIMAT==KIMAT2) then
           write(ifll,'(a,2I4)') 
     &    'ERR: Same materials on two side of Interface BC',
     &     KIMAT,KIMAT2
           STOP 'Not supportor this case in list_geom'
        endif
        if(MAT_CHKBC(nb,1).eq.0) then
          MAT_CHKBC(nb,1)=KIMAT
          write(ifll,'(a,I4,a,I4)')
     &    'MSG: BC no-1: ',nobcnd(nb),' IMAT_U=',KIMAT
        elseif(MAT_CHKBC(nb,1).ne.0.and.MAT_CHKBC(nb,1).ne.KIMAT) then
          write(ifle,'(2a,3I6,2a,1x,a)') '  *** ERR ***  One BC ',
     &  'cannot corss on different Materials ',MAT_CHKBC(nb,1),KIMAT,
     &   nobcnd(nb)
     &  ,'BC name: ',trim(boundName(nb,1)),trim(boundName(nb,2))
          stop 'list_material-2'
        endif
        if(MAT_CHKBC(nb,2).eq.0) then
          MAT_CHKBC(nb,2)=KIMAT2
          write(ifll,'(a,I4,a,I4)')
     &    'MSG: BC no-2: ',nobcnd(nb),' IMAT_U=',KIMAT2
        elseif(MAT_CHKBC(nb,2).ne.0.and.MAT_CHKBC(nb,2).ne.KIMAT2) then
          write(ifle,'(a,3I8,a,a,2x,a)') 
     &    'ERR: Interface BC crosses two different materials',
     &    MAT_CHKBC(nb,2),KIMAT2,nobcnd(nb)
     &  ,' BC name: ',trim(boundName(nb,1)),trim(boundName(nb,2))
          stop 'NOT Support this kind mesh, list_material-3'
        endif
      elseif(kd==kdsld.or.kd==kdpors) then
        if(MAT_CHKBC(nb,1).eq.0) then
          MAT_CHKBC(nb,1)=KIMAT
          write(ifll,'(a,I4,a,I4)') 'MSG: BC no: ',nobcnd(nb),
     &    ' IMAT_U=',KIMAT
        elseif(MAT_CHKBC(nb,1).ne.0.and.MAT_CHKBC(nb,1).ne.KIMAT) then
          write(ifle,'(2a,3I4,1x,2a,1x,a)') '  *** ERR ***  One BC ',
     &  'cannot corss on different Materials ',MAT_CHKBC(nb,1),KIMAT,
     &   nobcnd(nb)
     &  ,'BC name: ',trim(boundName(nb,1)),trim(boundName(nb,2))
          stop 'list_material-2'
        endif
        if(MAT_CHKBC(nb,2).eq.0) then
          MAT_CHKBC(nb,2)=KIMAT2
          write(ifll,'(a,I4,a,I4)') 'MSG: BC no: ',nobcnd(nb),
     &    ' IMAT_U=',KIMAT2
        elseif(MAT_CHKBC(nb,2).ne.0.and.MAT_CHKBC(nb,2).ne.KIMAT2) then
           write(ifle,'(a,3I4,3a)') 
     &    'ERR: *** Sliding BC crosses two different materials',
     &     MAT_CHKBC(nb,2),KIMAT,nobcnd(nb),
     &    'BC name: ',trim(boundName(nb,1)),trim(boundName(nb,2))
           stop 'NOT Support this kind mesh, list_material-3'
        endif
      endif
      enddo
 630  continue
!
      MAT_BC(:,:)=0
      DO NB=1,nbcnd
      do i=1,2
      KIMAT=MAT_CHKBC(nb,i)
      DO IIMAT=1,NMAT
      IMAT_U=MAT_NO(IIMAT)
!      IMAT_U=MAT_ALL(IIMAT)
      IF(KIMAT.EQ.IMAT_U) THEN
        MAT_BC(NB,i)=IIMAT
      ENDIF
      ENDDO
      enddo
      ENDDO
!------------------------
!      DO NB=1,nbcnd
!      do i=1,2
!        IIMAT=MAT_BC(NB,i)
!        IMAT=MAT_NO(IIMAT)
!      enddo
!      enddo
!------------------------
! --- interface BC: 
!------------------------
!      do 640 IBF=1,NSSFBC
!        NB=LBCSSF(IBF)
!        if(NB.lt.0) goto 640
!        if(kdbcnd(0,NB).eq.kdintr) then
!          IBFP=LCYCSF(IBF)
!          if(LBCSSF(IBFP).lt.0) LBCSSF(IBFP)=-LBCSSF(IBFP)
!        endif
! 640  continue
!-----------------------------------------
! --- output geom file (Serial and HPC) 
!-----------------------------------------
      if(icpu.eq.0) then
        fnam='geom'
      else
        my_rank=icpu-1
        HEADER1='geom_hpc'
        call DEFINE_FILE_NAME(DIRHED,HEADER1,my_rank,GEOMout(icpu))
        fnam=GEOMout(icpu)
      endif
!
! --- 
!
      IF(icpu==0) then
        locmsh(:)=1
      else
        
      endif
!
!
!      NMAT=0
!      MAT_NO(:)=0
!      do IMAT=-100,100
!      if(IMAT>0) then
!        do ifd=1,nflud
!        if(nofld(ifd)==IMAT) then
!          NMAT=NMAT+1
!          MAT_NO(NMAT)=ifd
!        endif
!        enddo
!      elseif(IMAT<0) then
!        do ifd=1,nsold
!        if(nosld(ifd)==IMAT) then
!          NMAT=NMAT+1
!          MAT_NO(NMAT)=-ifd
!        endif
!        enddo
!      endif
!      enddo
!
      ifl=50
      open(ifl,file=fnam,FORM='unformatted',status='unknown',
     &            iostat=ios)
      write(ifl)   icon_cv
      write(ifl)   NMAT
      write(ifl)  (MAT_NO(IIMAT),IIMAT=1,NMAT)
!      write(ifl)  (MAT_ALL(IIMAT),IIMAT=1,NMAT)
      write(ifl)  (MAT_CVIDX(IIMAT),IIMAT=0,NMAT)
      write(ifl)  (MAT_CVEXT(IIMAT),IIMAT=0,NMAT)
      write(ifl)  (MAT_DCIDX(IIMAT),IIMAT=0,NMAT)
      write(ifl)  (MAT_CFIDX(IIMAT),IIMAT=0,NMAT)
!
!----------------------------------------------------------------
!

      do 650 IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      write(ifl)  IIMAT,ICVS,ICVE
      if(ICVS<=ICVE) then
        write(ifl)  (MAT_CV(ICVL),ICVL=ICVS,ICVE)
        write(ifl)  (LVRTCV(MAT_CV(ICVL)),ICVL=ICVS,ICVE)
        write(ifl)  (CVVOLM(MAT_CV(ICVL)),ICVL=ICVS,ICVE)
        write(ifl) ((CVCENT(I,MAT_CV(ICVL)),I=1,3),ICVL=ICVS,ICVE)
      endif
!
      IDCS=MAT_DCIDX(IIMAT-1)+1
      IDCE=MAT_DCIDX(IIMAT)
      write(ifl)  IIMAT,IDCS,IDCE
      if(IDCS<=IDCE) then
        write(ifl) (MAT_CV(ICVL),ICVL=IDCS,IDCE)
        write(ifl) (LVRTCV(MAT_CV(ICVL)),ICVL=IDCS,IDCE)
        write(ifl) (CVVOLM(MAT_CV(ICVL)),ICVL=IDCS,IDCE)
        write(ifl)((CVCENT(I,MAT_CV(ICVL)),I=1,3),ICVL=IDCS,IDCE)
      endif
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      write(ifl)  IIMAT,ICFS,ICFE
      if(ICFS<=ICFE) then
        write(ifl)((LCV_CV(LVEDGE(I,MAT_CF(ICFL))),
     &               I=1,2),ICFL=ICFS,ICFE)
        write(ifl)((SFAREA(I,MAT_CF(ICFL)),I=1,4),ICFL=ICFS,ICFE)
        write(ifl)((SFCENT(I,MAT_CF(ICFL)),I=1,3),ICFL=ICFS,ICFE)
      endif
 650  continue
!----------------------------------------------------------------
      write(ifl) nbcnd
      write(ifl) (LBC_INDEX(nb),nb=0,nbcnd)
      write(ifl) (LBC_pair(nb),nb=0,nbcnd)
      write(ifl)((MAT_BC(nb,i),nb=1,nbcnd),i=1,2)
      DO nb=1,nbcnd
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      write(ifl) nb,IBFS,IBFE 
      if(IBFS<=IBFE) then
        write(ifl)(LBC_SSF(IBFL),IBFL=IBFS,IBFE)
        write(ifl)
     &  (LCF_CF(LCYCSF(MAT_CF(LBC_SSF(IBFL))-NEDGE)+NEDGE),
     &             IBFL=IBFS,IBFE)
      endif
      enddo
      write(ifl)(LCV_CV(ICV),ICV=1,NALLCV)
      write(ifl) radflag
!
! ---------------------------------------------------
!
      IF(icpu.EQ.0.AND.NPETOT.GT.1) GOTO 1234
      IF(radflag==1.AND.radmodel(1:3).NE.'FVM') THEN
        IF(NALLCV.LT.MatBounEnd(NMAT+nbcnd)) THEN
          WRITE(ifle,*) 'vertex all for rad =',MatBounEnd(NMAT+nbcnd)
	  WRITE(ifle,*) 'vertex all for fluid =',NALLCV
	  STOP 'vertex number different'
        ENDIF
!
        allocate(LFUTAU(NALLCV),stat=ierr1)
        if(ierr1.ne.0)  stop 'allocating LFUTAU error-1'
        LFUTAU=0
!------------------------------------
! --- Rebuild the NPALL index to ICV 
!------------------------------------
        NWK(NWK2,:)=0
!
        DO IIMAT=1,NMAT
        ICV1 = MatBounEnd(IIMAT)-MatBounEnd(IIMAT-1)
	ICV2 = MAT_CVEXT(IIMAT)-MAT_CVEXT(IIMAT-1)
        IF(ICV1.GT.ICV2) THEN
          WRITE(ifle,*) 'vertex num is different for mat ',IIMAT
	  WRITE(ifle,*) ICV1,ICV2
	  STOP 
        ENDIF
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        NWK(NWK1,:)=0

	DO ICV=ICVS,ICVE
	I=MAT_CV(ICV)
	NWK(NWK1,I)=ICV
        ENDDO
        DO IRD=MatBounEnd(IIMAT-1)+1,MatBounEnd(IIMAT)
          I=VIndex(IRD)
          ICV=NWK(NWK1,I)
          IF(ICV.LE.0) STOP 'not mapped rad => fluid'
          VIndex(IRD)=ICV
          LFUTAU(I)=IRD
        ENDDO
        DO ICV=ICVS,ICVE
          I=MAT_CV(ICV)
          IRD=LFUTAU(I)
          IF(IRD.LE.0) STOP 'not mapped fluid => rad'
          NWK(NWK2,ICV)=IRD
        ENDDO
        ENDDO
!
!
!
        DO nb = 1,nbcnd
        ICV1 = MatBounEnd(NMAT+NB)-MatBounEnd(NMAT+NB-1)
        ICV2 = LBC_INDEX(nb)-LBC_INDEX(nb-1)
        IF(ICV1.GT.ICV2) THEN
          WRITE(ifle,*) 'vertex num is different for boun ',nb
          WRITE(ifle,*) ICV1,ICV2
          STOP
        END IF
        NWK(NWK1,:)=0  !NWK1=0
        NWK(NWK3,:)=0  !NWK3=0
        LFUTAU=0
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        DO IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV1=LVEDGE(1,MAT_CF(ICFL))			
!
        IDC1=LVEDGE(2,MAT_CF(ICFL))
        IDC =LCV_CV(IDC1)
        NWK(NWK1,ICV1)=IDC  !NWK1(ICV1)=IDC
!
! --- hpc -> cal. index		   ----|
!
        ENDDO
        kd=kdbcnd(0,nb)
        IF(kd.eq.kdprdc.or.kd.eq.kdintr.or.kd.eq.kdsld) THEN
          DO IBFL=IBFS,IBFE
          ICFP=LCYCSF(MAT_CF(LBC_SSF(IBFL))-NEDGE)+NEDGE
          ICV2=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          IDC=LCV_CV(IDCP)
          NWK(NWK1,ICV2)=IDC  !NWK1(ICV2)=IDC 
! --- hpc -> cal. index   ----|
          ICFL= LBC_SSF(IBFL)
          ICV1= LVEDGE(1,MAT_CF(ICFL))			
          NWK(NWK3,ICV1)=ICV2  !NWK3(ICV1) = ICV2
          NWK(NWK3,ICV2)=ICV1  !NWK3(ICV2) = ICV1
          ENDDO
        ENDIF
        DO IRD=MatBounEnd(NMAT+NB-1)+1, MatBounEnd(NMAT+NB)
        I=VIndex(IRD)
        IDC=NWK(NWK1,I)   !IDC=NWK1(I)
        IF(IDC.LE.0) THEN
          if(NWK(NWK3,I).GT.0) then
            I=NWK(NWK3,I)     !NWK3(I)
            IDC=NWK(NWK1,I)   !NWK1(I)
            if(IDC.LE.0) 
     &     STOP 'not mapped boundary,rad => both fluid'
          else
            STOP 'not mapped boundary,rad => fluid'
          endif
        ENDIF
        VIndex(IRD)=IDC ! --- rad -> fluid cal. index	<--|
        NWK(NWK2,IDC)=IRD   !NWK2(IDC)=IRD
        LFUTAU(I)=IRD     !VWK1(I)=IRD
!pair
        if(NWK(NWK3,I).GT.0) then
          I=NWK(NWK3,I)   !NWK3(I)
          IDC=NWK(NWK1,I)   !NWK1(I)
	  NWK(NWK2,IDC)=IRD  !NWK2(IDC)=IRD
	  LFUTAU(I)=IRD
        endif
        END DO
        DO IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV1=LVEDGE(1,MAT_CF(ICFL))
        ICV=LCV_CV(ICV1)
        IDC1=LVEDGE(2,MAT_CF(ICFL))
        IDC=LCV_CV(IDC1)
        IRD=LFUTAU(ICV1)
        if(IRD.LE.0) then
          STOP 'not mapped boundary, fluid => rad'
        endif
        NWK(NWK2,IDC)=IRD !NWK2(IDC)=IRD 
!
        if(NWK(NWK3,ICV1).GT.0) then
          ICV2=NWK(NWK3,ICV1)
	  IDC=NWK(NWK1,ICV2) 
	  NWK(NWK2,IDC)=IRD  
        endif
        END DO
        END DO
!node for communication in parallel computation (import nodes)
!the transfer data is imaginary temperature, for a node located in the interface of 
!two materials, we still consider they have one unique imaginary temperature.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!///////////////////////////////////////////////////////////////
! for communication,  export(ICV) --> import(ICV) --> (IP) --> RD
! the VIndex(NALLCV+1:NPALL) is kept unchanged 
!	  IndexToExpt(1,:) is kept unchanged 
!///////////////////////////////////////////////////////////////
        NWK(NWK1,:)=0
!
        DO ICV=1,MAT_CVEXT(NMAT)
        I=MAT_CV(ICV)
        NWK(NWK1,I)=ICV  !NWK1(I)=ICV
        END DO
!
        DO IRD=MatBounEnd(NMAT+nbcnd)+1, NPALL
        I=VIndex(IRD)
        IF(NWK(NWK1,I).LE.0) then
          STOP 'not mapped rad bound <=> fluid'
        endif
        VIndex(IRD)=NWK(NWK1,I) !VIndex(IRD)=NWK1(I)
!rad -> fluid cal. index	<--|
        ENDDO
! --- parameters
        WRITE(ifl) StackOpt
        WRITE(ifl) IFlagGroup
        WRITE(ifl) NRD,NRD0,MRD,MID
        WRITE(ifl) NPALL,NDatExpt,NumTran,NMAT+nbcnd+1
! --- Exchange area
        WRITE(ifl) (RDIndex(I),I=1,NRD+1)
        WRITE(ifl) (RDValue(I),I=1,MRD)
        IF(StackOpt(1:4).EQ.'BAND') THEN
	  WRITE(ifl) (RDN1(I),I=1,NRD)
	  WRITE(ifl) (RDN2(I),I=1,NRD)
	  WRITE(ifl) (RDId(I),I=1,MRD)
        ENDIF
! --- Value transfer from node to RD  (1)
        WRITE(ifl) (VINDEX(I),I=1,NPALL)
        WRITE(ifl) (NWK(NWK2,ICV),ICV=1,NALLCV)
        WRITE(ifl) (AreaVol(I),I=1,NPALL)
        WRITE(ifl) (PROP(I),I=1,NPALL)
! --- data transfer  (2)
        WRITE(ifl) (MatBounEnd(I),I=0,NMAT+nbcnd+1)
        WRITE(ifl) (CoefToRD(I),I=1,NumTran)
        WRITE(ifl) (IndexToRD(1,I),I=1,NumTran)	!from ICV / IP
        WRITE(ifl) (IndexToRD(2,I),I=1,NumTran)	!to IRD
! --- export nodes  (3)
        IF(NDatExpt.GT.0) THEN
	  WRITE(ifl) (RaToExpt(I),I=1,NDatExpt)
	  WRITE(ifl) (IndexToExpt(1,I),I=1,NDatExpt) !from IP 
	  WRITE(ifl) (IndexToExpt(2,I),I=1,NDatExpt) !to IRD
          WRITE(ifl) (ExptID(I),I=1,NRD)	     !mark IRD 
        ENDIF
! --- release arrays
        DEALLOCATE(VIndex,AreaVol,PROP)	
        DEALLOCATE(RDN1,RDN2,RDId,RDIndex,RDValue)
        DEALLOCATE(IndexToRD,CoefToRD,MatBounEnd)
        IF(NDatExpt.GT.0) DEALLOCATE(IndexToExpt,RaToExpt,ExptID)
        deallocate(LFUTAU)
      ENDIF
 1234 continue
      close(ifl)  !LFUTAU
!----------------------------------------------------------------
! --- Serial LCV_CV and Globle Sliding BC temp file
!----------------------------------------------------------------
      if(icpu==0) then
        open(11,file='MAT_CVG',FORM='unformatted',status='unknown',
     &            iostat=ios)
        write(11)  NCV
        write(11) (MAT_CV(ICVL),ICVL=1,NCV)
        close(11)
      endif
!-----------------------------
! --- HPC domain
!-----------------------------
      if(icpu.ne.0) then
        my_rank=icpu-1
        HEADER1='reod'
        call DEFINE_FILE_NAME(DIRHED,HEADER1,my_rank,RODRICV(icpu))
        fnam=RODRICV(icpu)
        ifl=50
        open(ifl,file=fnam,FORM='unformatted',status='unknown',
     &            iostat=ios)
        write(ifl) (LCV_CV(ICV),ICV=1,NALLCV)
        write(ifl) (MAT_CV(ICV),ICV=1,NALLCV)
!
        close(ifl)
!
      endif
!------------------------------------
! --- Serial wall distance reording
!------------------------------------
      if(rewall.eq.1) then
        if(icpu.eq.0) then
          fnamwl='distance'
        else
          my_rank=icpu-1
          HEADER1='wall_hpc'
          call DEFINE_FILE_NAME(DIRHED,HEADER1,my_rank,walldis(icpu))
          fnamwl=walldis(icpu)
        endif
!---------------------------------------------------------------------
        open (15,FILE=fnamwl,FORM='UNFORMATTED',status='unknown',
     &        iostat=ios)
        if(ios.ne.0) stop 'open distance file error'
        read (15,ERR=1200) NCVX
        allocate(DISALL(NCVX),LFUTAU(NCVX),stat=ierr1)
        DISALL=0.d0
        LFUTAU=0
        if(ierr1.ne.0)  stop 'allocating DISALL error'
        read (15,ERR=1200) (DISALL(ICV),ICV=1,NCVX)
        read (15,ERR=1200) (LFUTAU(ICV),ICV=1,NCVX)
        close(15,ERR=1200)
!---------------------------------------------------------------------
        open (15,FILE=fnamwl,FORM='UNFORMATTED',status='unknown',
     &        iostat=ios)
        write (15)  NCVX
        write (15) (DISALL(MAT_CV(ICVL)),ICVL=1,NCVX)
        write (15) (LFUTAU(MAT_CV(ICVL)),ICVL=1,NCVX)
        close(15)
!---------------------------------------------------------------------
        deallocate(DISALL,LFUTAU)
      endif
!
! --- 
!
      if(ical_mvmsh>0) then
        DO IBF=1,nssfbc_old
           ICF=NEDGE+listbc(4,IBF)
           ICFL=LCF_CF(ICF)
           listbc(4,IBF)=ICFL
           
        enddo
      endif
!
      if(icpu==0) then   !HPC speed-up
!      if(icpu==0.and..false.) then
        if(NPETOT>1) then
          LCF_CF(:)=0

          do IBF=1,nssfbc
          NBSSF=LBCSSF(IBF)
          if(LBC_SSF(IBF)>0) then
            LBCSSF(IBF)=MAT_CF(LBC_SSF(IBF))
          endif
          enddo
!
          do IBF=1,nssfbc
          if(LBC_SSF(IBF)>0) then
            LCF_CF(IBF)=LCYCSF(MAT_CF(LBC_SSF(IBF))-NEDGE)+NEDGE
          endif
          enddo
!HORI
          do IBF=1,nssfbc
          LCYCSF(IBF)=LCF_CF(IBF)
          enddo
! --- s
          DO nb=1,nbcnd
          kd=kdbcnd(0,NB)
          conntBC(nb)=
     &       (idis(nb)==0.and.
     &         (kd==kdintr.or.kd==kdbuff.or.
     &          kd==kdshutr.or.kd==kdpors)
     &       )
!     &        .or.
!     &       (kd==kdsld)
          if(conntBC(nb)) then
            IBFS=LBC_INDEX(nb-1)+1
            IBFE=LBC_INDEX(nb)
            DO IBFL=IBFS,IBFE
            ICFL=LBCSSF(IBFL)
            ICFP=LCYCSF(IBFL)
            ICV2=LVEDGE(1,ICFP)
            ICV1=LVEDGE(1,ICFL)
            LVEDGE(2,ICFP)=ICV1
            LVEDGE(2,ICFL)=ICV2
            enddo
          endif
          enddo
        endif
      endif
!
      deallocate(MAT_CV,MAT_CF,LCF_CF,LCV_CV,LBC_SSF,LBC_pair)
      deallocate(MAT_BC,MAT_CHKBC,nsld) 
!
      return
!
 1200 stop ' distance file read error in list_output_geom'
!
      end subroutine list_output_geom
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_touch_inlet(mssfbc,NSSFBC,medge,nedge,mvrtx,
     &           LVEDGE,LMAT,LCYCSF,LBCSSF,SFCENT)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     calculation of "LCYCSF"
!
! --- [module arguments]
!
      use module_io,      only : ifll,ifle
      use module_hpc,     only : NPETOT
      use module_boundary,only : kdtchi,kdprdc,kdintr,kdsld,kdbuff,
     &                           kdshutr,kdpors,
     &                           LBC_INDEX,nbcnd,kdbcnd,
     &                           boundIDMap,prdcAxis,prdcOffset,
     &                           numbname,boundName,nobcnd,idis
      use module_partitioner,only : LBC_SSF => WK6
      use module_partitioner,only : LBC_pair=> WK1
      
      use module_param ,only   : scalfctr,smallNrm
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mssfbc,medge,mvrtx,NSSFBC,nedge
      integer,intent(inout) :: LBCSSF( mssfbc)
      integer,intent(inout) :: LCYCSF( mssfbc)
      real*8 ,intent(in)    :: SFCENT(3,medge)
      integer,intent(in)    :: LMAT  (  mvrtx)
      integer,intent(in)    :: LVEDGE(2,medge)
!
! --- [local entities]
!
      character*80 :: name2,namemain
      REAL*8       :: x1,y1,z1,x2,y2,z2,dis,disoff
      REAL*8       :: dism
      integer      :: IBFS,IBFE,IBFL,IBF,ICF,IMAT1,IMAT2,IMAT3,IBFT1
      integer      :: ICVA,ICVB
      integer      :: IBFST,IBFET,IBFLT,IBFT,ICFT,IBFS1,
     &                IBFE1,IBFS2,IBFE2
      integer      :: nb,nbtch,kd,kd1,j,i,k,tchno=0,sldno=0
      integer      :: nbtouch=0,ierr1,IND,NBSSF,NBSSF1,IBFP,IBCTCI
      integer      :: IDUM,IMINBC
!----------------------------
! --- --------
!----------------------------
      allocate(LBC_SSF(nssfbc),stat=ierr1)
      if(ierr1.ne.0) 
     &     stop 'stop at allocating LBC_SSF(:) in list_touch_inlet'
      allocate(LBC_pair(0:nbcnd),stat=ierr1)
      if(ierr1.ne.0) 
     &     stop 'stop at allocating LBC_pair(:) in list_touch_inlet'
!
!----------------------------
! --- list "touch-inlet" BC 
!----------------------------
!
      LBC_SSF=0
      LBC_INDEX=0
!
      do 210 nb=1,nbcnd
      IND=0
      do 200 IBF=1,NSSFBC
      NBSSF=LBCSSF(IBF)
      if(nb==NBSSF) then
        IND=IND+1
      endif
 200  continue
      LBC_INDEX(nb)=IND
 210  continue
!
      do 220 nb=1,nbcnd
      LBC_INDEX(nb)=LBC_INDEX(nb)+LBC_INDEX(nb-1)
 220  continue
!
      DO 230 nb=1,nbcnd
      IBFL=LBC_INDEX(nb-1)
      do 240 IBF=1,NSSFBC
      NBSSF=LBCSSF(IBF)
      if(nb==NBSSF) then
        IBFL=IBFL+1
        LBC_SSF(IBFL)=IBF
      ENDIF
 240  continue
 230  continue
!------------------
! --- Touch inlet
!------------------
      do 100 nb=1,nbcnd
      nbtouch=0
      tchno=0
      kd=kdbcnd(0,nb)
      namemain=boundName(nb,2)
      if(kd==kdtchi.and.idis(nb)==0) then
        disoff=dsqrt(prdcOffset(nb,1)**2
     &              +prdcOffset(nb,2)**2
     &              +prdcOffset(nb,3)**2)
        do 110 j=1,numbname
        do 140 nbtch=1,nbcnd
        
!        if((IBFE-IBFS).ne.(IBFET-IBFST)) goto 140
        if(nbtch.ne.nb) then
          name2=boundName(nbtch,j)
          if(adjustl(name2).eq.adjustl(namemain)) then
            if(j.eq.1.or.j==2) then
                IBFS=LBC_INDEX(nb-1)+1
                IBFE=LBC_INDEX(nb)
                do 160 IBFL=IBFS,IBFE
                tchno=tchno+1
                IBF=LBC_SSF(IBFL)
                ICF=IBF+nedge
                x1=SFCENT(1,ICF)
                y1=SFCENT(2,ICF)
                z1=SFCENT(3,ICF)
!
                IBFST=LBC_INDEX(nbtch-1)+1
                IBFET=LBC_INDEX(nbtch)
                dism=1.d20
                IMINBC=0

                do 150 IBFLT=IBFST,IBFET
                IBFT=LBC_SSF(IBFLT)
                ICFT=IBFT+nedge
!               
                x2=SFCENT(1,ICFT)
                y2=SFCENT(2,ICFT)
                z2=SFCENT(3,ICFT)
                if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
                  dis=ABS(dsqrt((x1-x2)*(x1-x2)*0.d0
     &                       +(y1-y2)*(y1-y2)
     &                       +(z1-z2)*(z1-z2)))
                elseif((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
                  dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                       +(y1-y2)*(y1-y2)*0.d0
     &                       +(z1-z2)*(z1-z2)))
                elseif((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
                  dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                       +(y1-y2)*(y1-y2)
     &                       +(z1-z2)*(z1-z2)*0.d0))
                endif
                dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                       +(y1-y2)*(y1-y2)
     &                       +(z1-z2)*(z1-z2))-disoff)
!                if(dis.lt.smallNrm*scalfctr(nb)) then
!                  nbtouch=nbtouch+1
!                  if(LCYCSF(IBF).ne.0) stop 'list_touch_inlet'
!                  LCYCSF(IBF)=IBFT
!                endif
                if(dis.lt.dism) then
                  IMINBC=IBFT
                  dism=dis
                endif
 150            continue
                if(IMINBC/=0) then
                  LCYCSF(IBF)=IMINBC
                  nbtouch=nbtouch+1
                else
                  stop'ERR-1: Touch-Inlet BC NOT matched'
                endif
!                if(LCYCSF(IBF)==0) 
!     &            stop'ERR-1: Touch-Inlet BC NOT matched'
 160            continue
!            elseif(j.eq.2) then
!              write(ifle,'(2a)') 
!     &      ' ERR: In Touch-Inlet, Vice BC (name2) must be defined by',
!     &        ' other Main BC (name)'
!              stop 'at list_touch_inlet'
            endif
          endif
        endif
 140    continue
 110    continue
!
        if(nbtouch.ne.tchno) then
          write(ifle,'(a,2I8)') 
     &     'ERR: touch BC NOT matched ',tchno,nbtouch
          write(ifle,'(2a)') 
     &     'MSG: You should check .offset value in fflow.ctl.',
     &       ' You should make the value after scale change.'
          write(ifle,'(a)') 
     &    '    or check the order of cyclic/touch-inlet'
          stop
        else
          write(ifll,'(a,I8,a,I8)')
     &    'MSG: TOUCH-INLET PAIR-BC SSF FACE NUMBER : ',
     &    nbtouch,' at BC no=',nobcnd(nb)
        endif
        IBCTCI=IBCTCI+nbtouch
      endif
!--------------------------------
! --- Discontinuous touch-inlet BC
!--------------------------------
      if(kd==kdtchi.and.idis(nb)==1) then
        tchno=0
        do 111 j=1,numbname
        do 141 nbtch=1,nbcnd
        if(nbtch.ne.nb) then
          name2=boundName(nbtch,j)
          if(adjustl(name2).eq.adjustl(namemain)) then
            if(j==1.or.j==2) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_INDEX(nb)
              do 161 IBFL=IBFS,IBFE
              dism=1.d20
              tchno=tchno+1
              IBF=LBC_SSF(IBFL)
              ICF=IBF+nedge
              x1=SFCENT(1,ICF)
              y1=SFCENT(2,ICF)
              z1=SFCENT(3,ICF)
!
              IBFST=LBC_INDEX(nbtch-1)+1
              IBFET=LBC_INDEX(nbtch)
              do 151 IBFLT=IBFST,IBFET
              IBFT=LBC_SSF(IBFLT)
              ICFT=IBFT+nedge
!
              x2=SFCENT(1,ICFT)
              y2=SFCENT(2,ICFT)
              z2=SFCENT(3,ICFT)
              if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
                dis=ABS(dsqrt((x1-x2)*(x1-x2)*0.d0
     &                       +(y1-y2)*(y1-y2)
     &                       +(z1-z2)*(z1-z2)))
              elseif((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
                dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                       +(y1-y2)*(y1-y2)*0.d0
     &                       +(z1-z2)*(z1-z2)))
              elseif((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
                dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                       +(y1-y2)*(y1-y2)
     &                       +(z1-z2)*(z1-z2)*0.d0))
              endif
!              dis=ABS(dsqrt((x1-x2)*(x1-x2)
!     &                     +(y1-y2)*(y1-y2)
!     &                     +(z1-z2)*(z1-z2)))
              if(dis.lt.dism) then
                nbtouch=nbtouch+1
                dism=dis
                IDUM=IBFT    !LCYCSF(IBF)=IBFT
              endif
 151          continue
              LCYCSF(IBF)=IDUM
              if(LCYCSF(IBF)==0) 
     &          stop'ERR-2: Discontinuous Touch-Inlet BC NOT matched'
 161          continue
!            elseif(j==2) then
!              write(ifle,'(2a)') 
!     &' ERR: In Touch-Inlet, Vice BC (name2) must be defined by',
!     &        ' other Main BC (name)'
!              stop 'at list_touch_inlet'
            endif
          endif
        endif
 141    enddo
 111    enddo
        IBCTCI=IBCTCI+tchno
        if(tchno/=0) then
          write(ifll,'(a,I8,a,I8)')
     &    'MSG: Discontinuous TOUCH-INLET SSF FACE NUMBER : ',
     &    tchno,' at BC no=',nobcnd(nb)
        endif
        if(tchno/=0.and.nbtouch==0) then
          stop 'ERR: Discontinuous touch-inlet Sliding'
        endif
      endif
!
 100  continue
!

!-----------------------------------------------------------
! --- "Interface BC" 
! --- list fluid/solid; solid/solid; fluid/fluid; interface
!-----------------------------------------------------------
!
      LBC_SSF=0
      LBC_INDEX=0
!
      do nb=1,nbcnd
      IND=0
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdintr.and.(idis(nb)==0.or.idis(nb)==1)) then
        if(nb.eq.NBSSF) then
          IND=IND+1
        endif
      endif
      enddo
      LBC_INDEX(nb)=IND
      enddo
!
      do nb=1,nbcnd
      IBFST=0
      IBFET=0
      kd1=kdbcnd(0,nb)
      do IBF=1,NSSFBC
      NBSSF=LBCSSF(IBF)
      kd=kdbcnd(0,abs(NBSSF))
      if(kd==kdintr.and.idis(nb)==0) then
        if(nb==NBSSF) then
          IBFST=IBFST+1
        endif
        if(nb==-NBSSF) then
          IBFET=IBFET+1
        endif
      endif
      enddo
      if(kd1==kdintr.and.idis(nb)==0) then
        write(ifll,'(3I6)') nb,IBFST,IBFET
        IF(IBFST/=IBFET) THEN
          WRITE(ifle,'(a)') 'ERR: Interface BC Pair No. NOT equal'
          WRITE(ifle,'(a)') 'MSG: Increasing [angle_ssf] in fflow.ctl'
          stop
        ENDIF
      endif
      enddo
!
      do nb=1,nbcnd
      LBC_INDEX(nb)=LBC_INDEX(nb)+LBC_INDEX(nb-1)
      enddo
!
      DO nb=1,nbcnd
      IBFL=LBC_INDEX(nb-1)
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdintr.and.(idis(nb)==0.or.idis(nb)==1)) then
        if(nb.eq.NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
      endif
      enddo
      enddo
!
      do 300 nb=1,nbcnd
      tchno=0
      nbtouch=0
      kd=kdbcnd(0,nb)
      if(kd==kdintr.and.idis(nb)==0) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do 320 IBFL=IBFS,IBFE
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(NBSSF1.lt.0.or.NBSSF1/=nb) goto 320
        if(LCYCSF(IBF).ne.0) goto 320
        nbtouch=nbtouch+1
        ICF=IBF+nedge
        ICVA=LVEDGE(1,ICF)
        IMAT1=LMAT(ICVA)
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        do 330 IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          if(LCYCSF(IBF).ne.0) goto 330
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 330
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.smallNrm*scalfctr(nb)) then
            tchno=tchno+1
            LCYCSF(IBF)=IBFT
            dism=dis
          endif
 330    continue
        if(LCYCSF(IBF)==0) then
          disoff=1.d10
          do 335 IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 335
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.disoff) then
            disoff=dis
            IBFT1=IBFT
          endif
 335      continue
          tchno=tchno+1
          LCYCSF(IBF)=IBFT1
        endif
 320    continue
        if(tchno/=nbtouch) then
          write(ifle,'(a,2I8,a,I8)') 
     &   'ERR: Interface BC NOT matched ',tchno,nbtouch,
     &                ' at BC no.=',nb
          write(ifle,*) 
     &     'MSG: You should check [tolerance] in fflow.ctl'
          stop
        else
          write(ifll,'(a,I8,a,I8)') 
     &   'MSG: Interface PAIR-BC SF FACE NUMBER : ',
     &    nbtouch,' at BC No=',nb
        endif
!----------------------------------
! --- Discontinuous interface BC
!----------------------------------
      elseif(kd==kdintr.and.idis(nb)==1) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do 321 IBFL=IBFS,IBFE
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(abs(NBSSF1)/=nb) cycle
        if(LCYCSF(IBF)/=0) cycle
        nbtouch=nbtouch+1
        ICF=IBF+nedge
        ICVA=LVEDGE(1,ICF)
        IMAT1=LMAT(ICVA)
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        dism=1.d20
        do 331 IBFLT=IBFS,IBFE
          if(IBFL==IBFLT) cycle
          IBFT=LBC_SSF(IBFLT)
          if(LCYCSF(IBF).ne.0) cycle
          NBSSF=LBCSSF(IBFT)

          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          if(IMAT2==IMAT1) cycle
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.dism) then
            tchno=tchno+1
            IBFT1=IBFT
            dism=dis
            IMAT3=IMAT2
          endif
 331    enddo
        LCYCSF(IBF)=IBFT1
        if(LCYCSF(IBF)==0) then
          disoff=1.d10
          do 336 IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
!     if(NBSSF.ne.-nb) cycle
          if(IBFL==IBFLT) cycle
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.disoff) then
            disoff=dis
            IBFT1=IBFT
          endif
 336      enddo
          tchno=tchno+1
          LCYCSF(IBF)=IBFT1
        endif
 321    enddo
      endif
 300  continue
!
!-----------------------------------------------------------
! --- "Buffle BC" 
! --- list continue fluid/fluid interface
!-----------------------------------------------------------
!
      LBC_SSF=0
      LBC_INDEX=0
!
      do nb=1,nbcnd
      IND=0
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdbuff.and.idis(nb)==0) then
        if(nb.eq.NBSSF) then
          IND=IND+1
        endif
      endif
      enddo
      LBC_INDEX(nb)=IND
      enddo
!
      do nb=1,nbcnd
      IBFST=0
      IBFET=0
      kd1=kdbcnd(0,nb)
      do IBF=1,NSSFBC
      NBSSF=LBCSSF(IBF)
      kd=kdbcnd(0,abs(NBSSF))
      if(kd==kdbuff.and.idis(nb)==0) then
        if(nb==NBSSF) then
          IBFST=IBFST+1
        endif
        if(nb==-NBSSF) then
          IBFET=IBFET+1
        endif
      endif
      enddo
      if(kd1==kdbuff.and.idis(nb)==0) then
        write(ifll,'(1x,a,3I6)') 
     &  'MSG: Buffle nb= A= B=',nb,IBFST,IBFET
        IF(IBFST/=IBFET) THEN
          WRITE(ifle,'(a)') 'ERR: Buffle BC Pair No. NOT equal'
          WRITE(ifle,'(a)') 'MSG: Increasing [angle_ssf] in fflow.ctl'
          stop
        ENDIF
      endif
      enddo
!
      do nb=1,nbcnd
      LBC_INDEX(nb)=LBC_INDEX(nb)+LBC_INDEX(nb-1)
      enddo
!
      DO nb=1,nbcnd
      IBFL=LBC_INDEX(nb-1)
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdbuff.and.idis(nb)==0) then
        if(nb.eq.NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
      endif
      enddo
      enddo
!
      do 600 nb=1,nbcnd
      tchno=0
      nbtouch=0
      kd=kdbcnd(0,nb)
      if(kd==kdbuff.and.idis(nb)==0) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do 620 IBFL=IBFS,IBFE
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(NBSSF1.lt.0.or.NBSSF1/=nb) goto 620
        if(LCYCSF(IBF).ne.0) goto 620
        nbtouch=nbtouch+1
        ICF=IBF+nedge
        ICVA=LVEDGE(1,ICF)
        IMAT1=LMAT(ICVA)
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        do 630 IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          if(LCYCSF(IBF).ne.0) goto 630
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 630
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.smallNrm*scalfctr(nb)) then
            tchno=tchno+1
            LCYCSF(IBF)=IBFT
            dism=dis
          endif
 630    continue
        if(LCYCSF(IBF)==0) then
          disoff=1.d10
          do 635 IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 635
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.disoff) then
            disoff=dis
            IBFT1=IBFT
          endif
 635      continue
          tchno=tchno+1
          LCYCSF(IBF)=IBFT1
        endif
 620    continue
        if(tchno/=nbtouch) then
          write(ifle,'(a,2I8,a,I8)') 
     &   'ERR: Buffle BC NOT matched ',tchno,nbtouch,
     &                ' at BC no.=',nb
          write(ifle,*) 
     &     'MSG: You should check [tolerance] in fflow.ctl'
          stop
        else
          write(ifll,'(a,I8,a,I8)') 
     &   'MSG: Buffle PAIR-BC SF FACE NUMBER : ',
     &    nbtouch,' at BC No=',nb
        endif
      endif
!----------------------------------
! --- Discontinuous interface BC
!----------------------------------
      if(kd==kdbuff.and.idis(nb)==1) then
        stop 'ERR: Discontinuous Buffle BC is NOT supported'
      endif
 600  continue
!
!-----------------------------------------------------------
! --- "Shutter BC" 
! --- list fluid/solid; solid/solid; fluid/fluid; interface
!-----------------------------------------------------------
      LBC_SSF=0
      LBC_INDEX=0
!
      do nb=1,nbcnd
      IND=0
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdshutr.and.idis(nb)==0) then
        if(nb.eq.NBSSF) then
          IND=IND+1
        endif
      endif
      enddo
      LBC_INDEX(nb)=IND
      enddo
!
      do nb=1,nbcnd
      IBFST=0
      IBFET=0
      kd1=kdbcnd(0,nb)
      do IBF=1,NSSFBC
      NBSSF=LBCSSF(IBF)
      kd=kdbcnd(0,abs(NBSSF))
      if(kd==kdshutr.and.idis(nb)==0) then
        if(nb==NBSSF) then
          IBFST=IBFST+1
        endif
        if(nb==-NBSSF) then
          IBFET=IBFET+1
        endif
      endif
      enddo
      if(kd1==kdshutr.and.idis(nb)==0) then
        write(ifll,'(1x,a,3I6)') 
     &  'MSG: Buffle nb= A= B=',nb,IBFST,IBFET
        IF(IBFST/=IBFET) THEN
          WRITE(ifle,'(a)') 'ERR: Buffle BC Pair No. NOT equal'
          WRITE(ifle,'(a)') 'MSG: Increasing [angle_ssf] in fflow.ctl'
          stop
        ENDIF
      endif
      enddo
!
      do nb=1,nbcnd
      LBC_INDEX(nb)=LBC_INDEX(nb)+LBC_INDEX(nb-1)
      enddo
!
      DO nb=1,nbcnd
      IBFL=LBC_INDEX(nb-1)
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdshutr.and.idis(nb)==0) then
        if(nb.eq.NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
      endif
      enddo
      enddo
!
      do 700 nb=1,nbcnd
      tchno=0
      nbtouch=0
      kd=kdbcnd(0,nb)
      if(kd==kdshutr.and.idis(nb)==0) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do 720 IBFL=IBFS,IBFE
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(NBSSF1.lt.0.or.NBSSF1/=nb) goto 720
        if(LCYCSF(IBF).ne.0) goto 720
        nbtouch=nbtouch+1
        ICF=IBF+nedge
        ICVA=LVEDGE(1,ICF)
        IMAT1=LMAT(ICVA)
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        do 730 IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          if(LCYCSF(IBF).ne.0) goto 730
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 730
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.smallNrm*scalfctr(nb)) then
            tchno=tchno+1
            LCYCSF(IBF)=IBFT
            dism=dis
          endif
 730    continue
        if(LCYCSF(IBF)==0) then
          disoff=1.d10
          do 735 IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 735
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.disoff) then
            disoff=dis
            IBFT1=IBFT
          endif
 735      continue
          tchno=tchno+1
          LCYCSF(IBF)=IBFT1
        endif
 720    continue
        if(tchno/=nbtouch) then
          write(ifle,'(a,2I8,a,I8)') 
     &   'ERR: Buffle BC NOT matched ',tchno,nbtouch,
     &                ' at BC no.=',nb
          write(ifle,*) 
     &     'MSG: You should check [tolerance] in fflow.ctl'
          stop
        else
          write(ifll,'(a,I8,a,I8)') 
     &   'MSG: Porous PAIR-BC SF FACE NUMBER : ',
     &    nbtouch,' at BC No=',nb
        endif
      endif
!----------------------------------
! --- Discontinuous interface BC
!----------------------------------
      if(kd==kdshutr.and.idis(nb)==1) then
        stop 'ERR: Discontinuous Buffle BC is NOT supported'
      endif
 700  continue
!
!-----------------------------------------------------------
! --- "Porous BC" 
! --- list fluid/solid; solid/solid; fluid/fluid; interface
!-----------------------------------------------------------
      LBC_SSF=0
      LBC_INDEX=0
!
      do nb=1,nbcnd
      IND=0
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdpors.and.idis(nb)==0) then
        if(nb.eq.NBSSF) then
          IND=IND+1
        endif
      endif
      enddo
      LBC_INDEX(nb)=IND
      enddo
!
      do nb=1,nbcnd
      IBFST=0
      IBFET=0
      kd1=kdbcnd(0,nb)
      do IBF=1,NSSFBC
      NBSSF=LBCSSF(IBF)
      kd=kdbcnd(0,abs(NBSSF))
      if(kd==kdpors.and.idis(nb)==0) then
        if(nb==NBSSF) then
          IBFST=IBFST+1
        endif
        if(nb==-NBSSF) then
          IBFET=IBFET+1
        endif
      endif
      enddo
      if(kd1==kdpors.and.idis(nb)==0) then
        write(ifll,'(1x,a,3I6)') 
     &  'MSG: Porous nb= A= B=',nb,IBFST,IBFET
        IF(IBFST/=IBFET) THEN
          WRITE(ifle,'(a)') 'ERR: Porous BC Pair No. NOT equal'
          WRITE(ifle,'(a)') 'MSG: Increasing [angle_ssf] in fflow.ctl'
          stop
        ENDIF
      endif
      enddo
!
      do nb=1,nbcnd
      LBC_INDEX(nb)=LBC_INDEX(nb)+LBC_INDEX(nb-1)
      enddo
!
      DO nb=1,nbcnd
      IBFL=LBC_INDEX(nb-1)
      do IBF=1,NSSFBC
      NBSSF=ABS(LBCSSF(IBF))
      kd=kdbcnd(0,NBSSF)
      if(kd==kdpors.and.idis(nb)==0) then
        if(nb.eq.NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
      endif
      enddo
      enddo
!
      do nb=1,nbcnd
      tchno=0
      nbtouch=0
      kd=kdbcnd(0,nb)
      if(kd==kdpors.and.idis(nb)==0) then 
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(NBSSF1.lt.0.or.NBSSF1/=nb) cycle 
        if(LCYCSF(IBF).ne.0) cycle!goto 720
        nbtouch=nbtouch+1
        ICF=IBF+nedge
        ICVA=LVEDGE(1,ICF)
        IMAT1=LMAT(ICVA)
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
!        do IBFLT=IBFS,IBFE
!          IBFT=LBC_SSF(IBFLT)
!          if(LCYCSF(IBF).ne.0) cycle 
!          NBSSF=LBCSSF(IBFT)
!          if(NBSSF.ne.-nb) cycle
!          ICFT=IBFT+nedge
!          ICVB=LVEDGE(1,ICFT)
!          IMAT2=LMAT(ICVB)
!          x2=SFCENT(1,ICFT)
!          y2=SFCENT(2,ICFT)
!          z2=SFCENT(3,ICFT)
!          dis=ABS(dsqrt((x1-x2)*(x1-x2)
!     &                 +(y1-y2)*(y1-y2)
!     &                 +(z1-z2)*(z1-z2)))
!          if(dis.lt.smallNrm*scalfctr(nb)) then
!            tchno=tchno+1
!            LCYCSF(IBF)=IBFT
!            dism=dis
!          endif
!        enddo
        if(LCYCSF(IBF)==0) then
          disoff=1.d10
          do IBFLT=IBFS,IBFE
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) cycle
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                 +(y1-y2)*(y1-y2)
     &                 +(z1-z2)*(z1-z2)))
          if(dis.lt.disoff) then
            disoff=dis
            IBFT1=IBFT
          endif
          enddo
          tchno=tchno+1
          LCYCSF(IBF)=IBFT1
        endif
        enddo
        if(tchno/=nbtouch) then
          write(ifle,'(a,2I8,a,I8)') 
     &   'ERR: Porous BC NOT matched ',tchno,nbtouch,
     &                ' at BC no.=',nb
          write(ifle,*) 
     &     'MSG: You should check [tolerance] in fflow.ctl'
          stop
        else
          write(ifll,'(a,I8,a,I8)') 
     &   'MSG: Porous PAIR-BC SF FACE NUMBER : ',
     &    nbtouch,' at BC No=',nb
        endif
      endif
!----------------------------------
! --- Discontinuous interface BC 
!----------------------------------
      if(kd==kdpors.and.idis(nb)==1) then
        stop 'ERR: Discontinuous Porous BC is NOT supported'
      endif
      enddo
!
!
!------------------------
! --- sliding mesh BC: 
!------------------------
      LBC_SSF=0
      LBC_INDEX=0
      LBC_pair=0
!
      do nb=1,nbcnd
      kd=kdbcnd(0,NB)
      if(kd==kdsld) then
        LBC_INDEX(nb)=LBC_INDEX(nb-1)
        IND=0
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.NBSSF) then
          IND=IND+1
        endif
        enddo
!
        if(idis(nb)==1) LBC_pair(nb)=LBC_INDEX(nb)+IND
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.-NBSSF) then
          IND=IND+1
        endif
        enddo
        LBC_INDEX(nb)=LBC_INDEX(nb)+IND
      endif
      enddo
!
      DO nb=1,nbcnd
      kd=kdbcnd(0,NB)
      if(kd==kdsld) then
        IBFL=LBC_INDEX(nb-1)
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
        enddo
!
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.-NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
        enddo
      endif
      enddo
!
      DO 400 nb=1,nbcnd
      sldno=0
      tchno=0
      kd=kdbcnd(0,nb)
      if(kd==kdsld.and.idis(nb)==0) then
        IBFS1=LBC_INDEX(nb-1)+1
        IBFE1=IBFS1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
        IBFS2=IBFE1+1
        IBFE2=LBC_INDEX(nb)
        do 420 IBFL=IBFS1,IBFE1
        tchno=tchno+1
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(LCYCSF(IBF).ne.0) goto 420
        ICF=IBF+nedge
        ICVA=LVEDGE(1,ICF)
        IMAT1=LMAT(ICVA)
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        do 430 IBFLT=IBFS2,IBFE2
          IBFT=LBC_SSF(IBFLT)
          if(LCYCSF(IBFT).ne.0) goto 430
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 430
          ICFT=IBFT+nedge
          ICVB=LVEDGE(1,ICFT)
          IMAT2=LMAT(ICVB)
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=(dsqrt((x1-x2)*(x1-x2)
     &              +(y1-y2)*(y1-y2)
     &              +(z1-z2)*(z1-z2)))
          if(dis.lt.smallNrm*scalfctr(nb)) then
            sldno=sldno+1
            LCYCSF(IBF)=IBFT
            LCYCSF(IBFT)=IBF
          endif
 430    continue
        if(LCYCSF(IBF)==0) stop 'ERR: Sliding BC NOT matched'
        if(LCYCSF(IBFT)==0) stop 'ERR: Sliding BC NOT matched'
 420    continue
        if(tchno/=sldno) then
          write(ifle,'(a,2I8)') 
     &     'ERR: Sliding BC SF NOT matched : ',tchno,sldno
          write(ifle,'(a,I8)') 
     &     'MSG: BC no. in fflow.ctl : ',nobcnd(nb)
          stop
        else
         write(ifll,'(a,I8)') 
     &     'MSG: Sliding BC SSF face Number : ', sldno
        endif
      endif
!---------------------------------
! --- discontinuous sliding BC 
!---------------------------------
      if(kd==kdsld.and.idis(nb)==1) then
        IBFS1=LBC_INDEX(nb-1)+1
        IBFE1=LBC_pair(nb)
        IBFS2=LBC_pair(nb)+1
        IBFE2=LBC_INDEX(nb)
        tchno=0
        do 421 IBFL=IBFS1,IBFE1
        dism=1.d20
        tchno=tchno+1
        IBF=LBC_SSF(IBFL) 
        if(LCYCSF(IBF).ne.0) goto 421
        ICF=IBF+nedge
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        IDUM=0
        do 431 IBFLT=IBFS2,IBFE2
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 431
          ICFT=IBFT+nedge
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=(dsqrt((x1-x2)*(x1-x2)
     &              +(y1-y2)*(y1-y2)
     &              +(z1-z2)*(z1-z2)))
          if(dis.lt.dism) then
            dism=dis
            IDUM=IBFT
          endif
 431    continue
        LCYCSF(IBF)=IDUM
        if(LCYCSF(IBF)==0) 
     &    stop 'ERR: discontinuous Sliding BC-1 NOT matched'
 421    continue
        if(tchno/=0) then
          write(ifll,'(a,I8,a,I4,2a)')
     &    'MSG: Discontinuous Main sliding SSF FACE NUMBER : ',
     &    tchno,' BC no=',nobcnd(nb),
     &    ' BC name : ',trim(boundName(nb,1))
        elseif(NPETOT==1) then
          stop 'ERR-1: discontinuous Sliding'
        endif
!
        tchno=0
        do 422 IBFL=IBFS2,IBFE2
        dism=1.d20
        tchno=tchno+1
        IBF=LBC_SSF(IBFL)
        if(LCYCSF(IBF).ne.0) goto 422
        ICF=IBF+nedge
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        do 432 IBFLT=IBFS1,IBFE1
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.nb) goto 432
          ICFT=IBFT+nedge
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=(dsqrt((x1-x2)*(x1-x2)
     &              +(y1-y2)*(y1-y2)
     &              +(z1-z2)*(z1-z2)))
          if(dis.lt.dism) then
            dism=dis
            IDUM=IBFT
          endif
 432    continue
        LCYCSF(IBF)=IDUM
        if(LCYCSF(IBF)==0) 
     &    stop 'ERR: discontinuous Sliding BC-2 NOT matched'
 422    continue
        if(tchno/=0) then
          write(ifll,'(a,I8,a,I4,2a)')
     &    'MSG: Discontinuous Sub  sliding SSF FACE NUMBER : ',
     &     tchno,' BC no=',nobcnd(nb),
     &    ' BC name : ',trim(boundName(nb,2))
        elseif(NPETOT==1) then
          stop 'ERR-2: discontinuous Sliding'
        endif
      endif
 400  continue
!---------------------------------
! --- Discontinuous periodic BC
!---------------------------------
      LBC_SSF=0
      LBC_INDEX=0
      LBC_pair=0
!
      do nb=1,nbcnd
      kd=kdbcnd(0,NB)
      if(kd==kdprdc.and.idis(nb)==1) then
        LBC_INDEX(nb)=LBC_INDEX(nb-1)
        IND=0
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.NBSSF) then
          IND=IND+1
        endif
        enddo
!
        LBC_pair(nb)=LBC_INDEX(nb)+IND
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.-NBSSF) then
          IND=IND+1
        endif
        enddo
        LBC_INDEX(nb)=LBC_INDEX(nb)+IND
      endif
      enddo
!
      DO nb=1,nbcnd
      kd=kdbcnd(0,NB)
      if(kd==kdprdc.and.idis(nb)==1) then
        IBFL=LBC_INDEX(nb-1)
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
        enddo
        do IBF=1,NSSFBC
        NBSSF=LBCSSF(IBF)
        if(nb.eq.-NBSSF) then
          IBFL=IBFL+1
          LBC_SSF(IBFL)=IBF
        ENDIF
        enddo
      endif
      enddo
!
      DO 500 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if(kd==kdprdc.and.idis(nb)==1) then
        IBFS1=LBC_INDEX(nb-1)+1
        IBFE1=LBC_pair(nb)
        IBFS2=LBC_pair(nb)+1
        IBFE2=LBC_INDEX(nb)
        tchno=0
        do 521 IBFL=IBFS1,IBFE1
        dism=1.d20
        tchno=tchno+1
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(LCYCSF(IBF).ne.0) goto 521
        ICF=IBF+nedge
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        do 531 IBFLT=IBFS2,IBFE2
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.-nb) goto 531
          ICFT=IBFT+nedge
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=(dsqrt((x1-x2)*(x1-x2)
     &              +(y1-y2)*(y1-y2)
     &              +(z1-z2)*(z1-z2)))
          if(dis.lt.dism) then
            dism=dis
            IDUM=IBFT
          endif
 531    continue
        LCYCSF(IBF)=IDUM
        if(LCYCSF(IBF)==0) 
     &    stop 'ERR: discontinuous periodic BC-1 NOT matched'
 521    continue
        if(tchno/=0) then
          write(ifll,'(a,I8,a,I4,2a)')
     &    'MSG: Discontinuous Main periodic SSF FACE NUMBER : ',
     &    tchno,' BC no=',nobcnd(nb),
     &    ' BC name : ',trim(boundName(nb,1))
        elseif(NPETOT==1) then
          stop 'ERR-1: discontinuous periodic'
        endif
!
        tchno=0
        do 522 IBFL=IBFS2,IBFE2
        dism=1.d20
        tchno=tchno+1
        IBF=LBC_SSF(IBFL)
        NBSSF1=LBCSSF(IBF)
        if(LCYCSF(IBF).ne.0) goto 522
        ICF=IBF+nedge
        x1=SFCENT(1,ICF)
        y1=SFCENT(2,ICF)
        z1=SFCENT(3,ICF)
        do 532 IBFLT=IBFS1,IBFE1
          IBFT=LBC_SSF(IBFLT)
          NBSSF=LBCSSF(IBFT)
          if(NBSSF.ne.nb) goto 532
          ICFT=IBFT+nedge
          x2=SFCENT(1,ICFT)
          y2=SFCENT(2,ICFT)
          z2=SFCENT(3,ICFT)
          dis=(dsqrt((x1-x2)*(x1-x2)
     &              +(y1-y2)*(y1-y2)
     &              +(z1-z2)*(z1-z2)))
          if(dis.lt.dism) then
!!!!!          if(dis.lt.smallNrm*scalfctr(nb)) then
            dism=dis
            IDUM=IBFT
          endif
 532    continue
        LCYCSF(IBF)=IDUM
        if(LCYCSF(IBF)==0) 
     &    stop 'ERR: discontinuous periodic BC-2 NOT matched'
 522    continue
        if(tchno/=0) then
          write(ifll,'(a,I8,a,I4,2a)')
     &    'MSG: Discontinuous Sub  periodic SSF FACE NUMBER : ',
     &     tchno,' BC no=',nobcnd(nb),
     &    ' BC name : ',trim(boundName(nb,2))
        elseif(NPETOT==1) then
          stop 'ERR-2: discontinuous periodic'
        endif
      endif
 500  continue
!
      deallocate(LBC_SSF,LBC_pair)
!
      return
      end subroutine list_touch_inlet
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_pair
     &       (ncell,nface,mcell,mface,IBCTCI,lbface,gface)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     calculation of "lbface" for touch-inlet of HPC
!
! --- [module arguments]
!
      use module_io,      only  : ifll,ifle,cntlnam
      use module_boundary,only  : kdtchi,kdsld,
     &                            nbcnd,kdbcnd,
     &                            boundIDMap,prdcOffset,
     &                            numbname,boundName,kdintr,idis
      use module_param    ,only   : scalfctr,smallNrm
      use module_parameter,only   : NIFACE
      use module_partitioner,only : LBC_SSF=> WK6
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,mface,ncell,nface
      integer,intent(out)   :: IBCTCI
      integer,intent(inout) :: lbface(2,mface)
      real*8 ,intent(in)    :: gface (3,mface)
!
! --- [local entities]
!
      character*80 :: name2,namemain
      REAL*8       :: x1,y1,z1,x2,y2,z2,dis,disoff,dism
      integer      :: IBFS,IBFE,IBFL,IBF,ICF,ISB
      integer      :: LBC_IND(0:nbcnd)
      integer      :: IBFST,IBFET,IBFLT,IBFT,ICFT
      integer      :: nb,nbtch,j,i,k,tchno=0
      integer      :: nbtouch=0,ierr1=0
      integer      :: IS,NBSSF,kd,IST,ISP,NCHK,IMINBC
!
      LBC_IND=0
      do 210 nb=1,nbcnd
      ISB=0
      do 200 IS=NIFACE,nface
      NBSSF=lbface(1,IS)
      if(NBSSF.le.0) goto 200
      if(nb.eq.NBSSF) then
        ISB=ISB+1
      endif
 200  continue
      LBC_IND(nb)=ISB
 210  continue
      do 300 nb=1,nbcnd
      LBC_IND(nb)=LBC_IND(nb)+LBC_IND(nb-1)
 300  continue
!
      allocate(LBC_SSF(LBC_IND(nbcnd)),stat=ierr1) 
      if(ierr1.ne.0) 
     &     stop 'stop at allocating LBC_SSF(:) in list_pair'
      LBC_SSF=0
!
      DO 400 nb=1,nbcnd
      ISB=LBC_IND(nb-1)
      do 410 IS=NIFACE,nface
      NBSSF=lbface(1,IS)
      if(NBSSF.le.0) goto 410
      if(nb.eq.NBSSF) then
        ISB=ISB+1
        LBC_SSF(ISB)=IS
      ENDIF
 410  continue
 400  continue
!
      NCHK=0
      nbtouch=0
      tchno=0
      do 100 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      IBFS=LBC_IND(nb-1)+1
      IBFE=LBC_IND(nb)
      namemain=boundName(nb,2)
      if(kd==kdtchi.and.idis(nb)==0) then
        disoff=dsqrt(prdcOffset(nb,1)**2
     &              +prdcOffset(nb,2)**2
     &              +prdcOffset(nb,3)**2)
        do 110 j=1,numbname
        do 140 nbtch=1,nbcnd
        IBFST=LBC_IND(nbtch-1)+1
        IBFET=LBC_IND(nbtch)
        if((IBFE-IBFS).ne.(IBFET-IBFST)) goto 140
        if(nbtch.ne.nb) then
          name2=boundName(nbtch,j)
          if(adjustl(name2).eq.adjustl(namemain)) then
            if(j.eq.1) then
              IBFS=LBC_IND(nb-1)+1
              IBFE=LBC_IND(nb)
              do 225 IBFL=IBFS,IBFE
              tchno=tchno+1
              IS=LBC_SSF(IBFL)
              x1=gface(1,IS)
              y1=gface(2,IS)
              z1=gface(3,IS)
!
              dism=1.d20
              IMINBC=0
!
              IBFST=LBC_IND(nbtch-1)+1
              IBFET=LBC_IND(nbtch)
              do 235 IBFLT=IBFST,IBFET
              IST=LBC_SSF(IBFLT)
!
              x2=gface(1,IST)
              y2=gface(2,IST)
              z2=gface(3,IST)
              dis=ABS(dsqrt((x1-x2)*(x1-x2)
     &                     +(y1-y2)*(y1-y2)
     &                     +(z1-z2)*(z1-z2))-disoff)
              if(NCHK.eq.0) then
                NCHK=1
                write(ifll,*) 'MSG: The first CV OFFSET: ',
     &          ' X_offset= ',abs(x1-x2),
     &          ' Y_offset= ',abs(y1-y2),
     &          ' Z_offset= ',abs(z1-z2)
              endif
!              if(dis.lt.smallNrm*scalfctr(nb)) then
!                nbtouch=nbtouch+1
!                if(lbface(2,IS).ne.0) stop 'list_pair'
!                lbface(2,IS)=IST
!              endif
              if(dis.lt.dism) then
                dism=dis
                IMINBC=IST
              endif
 235          continue
              if(IMINBC/=0) then
                nbtouch=nbtouch+1
                lbface(2,IS)=IMINBC
              else
                stop 'list_pair'
              endif
 225          continue
            elseif(j.eq.2) then
              write(ifle,*) 
     &      ' ERR: In Touch-Inlet, Vice BC (name2) must be defined by',
     &        ' other Main BC (name)'
              stop
            endif
          endif
        endif
 140    continue
 110    continue
!      elseif(kd.eq.kdsld) then
!        disoff=1.D-15
!        do 310 j=1,numbname
!        do 340 nbsld=1,nbcnd
!        IBFST=LBC_IND(nbsld-1)+1
!        IBFET=LBC_IND(nbsld)
!        if((IBFE-IBFS).ne.(IBFET-IBFST)) goto 340
!        if(nbsld.ne.nb) then
!          name2=boundName(nbtch,j)
!          if(adjustl(name2).eq.adjustl(namemain)) then
!            if(j.eq.1) then
!              IBFS=LBC_IND(nb-1)+1
!              IBFE=LBC_IND(nb)
!              do IBFL=IBFS,IBFE
!              nosldng=nosldng+1
!              IS=LBC_SSF(IBFL)
!              x1=gface(1,IS)
!              y1=gface(2,IS)
!              z1=gface(3,IS)
!
!              IBFST=LBC_IND(nbsld-1)+1
!              IBFET=LBC_IND(nbsld)
!              do IBFLT=IBFST,IBFET
!              IST=LBC_SSF(IBFLT)
!
!              x2=gface(1,IST)
!              y2=gface(2,IST)
!              z2=gface(3,IST)
!              dis=ABS(dsqrt((x1-x2)*(x1-x2)
!     &                     +(y1-y2)*(y1-y2)
!     &                     +(z1-z2)*(z1-z2))-disoff)
!              if(dis.lt.1.D-10) then
!                nosldng1=nosldng1+1
!                if(lbface(2,IS).ne.0) stop 'list_pair'
!                lbface(2,IS)=IST
!              endif
!              enddo
!              enddo
!            endif
!          endif
!        endif
! 340    continue
! 310    continue
      endif
 100  continue
!
      if(nbtouch.ne.tchno) then
        write(ifle,*) ' ERR: touch BC NOT matched in list_pair',
     &                  tchno,nbtouch
       write(ifle,*)' Check offset value of touch-inlet BC in ',cntlnam
        stop
      else
        write(ifll,'(a,I10)')
     &   'MSG: TOUCH-INLET BC FACE NUMBER : ',tchno
      endif

      IBCTCI=nbtouch
      deallocate(LBC_SSF) 
!
      return
      end subroutine list_pair
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_wall_distance
     &          (icpu,rewall,NCV,NBCWAL,nssfbc,NEDGE,
     &           mvrtx,mssfbc,medge,
     &           LMAT,LBCSSF,SFCENT,CVCENT,SFAREA)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_io,         only : ifli,ifll,ifle
      use module_hpc,        only : NPETOT
      use module_boundary,   only : kdbcnd,kvnslp,kvfslp,kxnone,nbcnd,
     &                              LBC_INDEX,kvlglw,kvrogf,kvmodl
      use module_model,      only : icaltb,rns_scl,sles,dles,lles
      use module_material,   only : cvdist
      use module_partitioner,only : LFUTAU =>wk1
      use module_partitioner,only : DISALL =>WALLDIST
      use module_partitioner,only : LBC_SSF=>wk2
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: rewall
      integer,intent(in)  :: icpu,NCV,NBCWAL,nssfbc,NEDGE
      integer,intent(in)  :: mvrtx,mssfbc,medge
      integer,intent(in)  :: LMAT  (  mvrtx)
      integer,intent(in)  :: LBCSSF( mssfbc)
      real*8 ,intent(in)  :: SFCENT(3,medge)
      real*8 ,intent(in)  :: CVCENT(3,mvrtx)
      real*8 ,intent(in)  :: SFAREA(4,medge)
!
! --- [local entities]
!
      integer :: ierr=0
      integer :: ios=0
      integer :: NB,KD,kdv
      integer :: NTAU,ICFMIN
      logical :: wallfile=.false.
      real*8  :: dx,dy,dz,dismin,dis1,dum1,yy,vol
      real*8 ,parameter :: SML=1.D-25,ZERO=0.D0,GREAT=1.D25
      integer :: IBFS,IBFE,IBFL,ICFL,IBF,NBSSF,ICF,ICV,my_rank
      integer :: IMAT,IIMAT,ICFS,ICFE,ICVS,ICVE,ICVL
      character(len=1) :: recal
      character(len=80),save :: fnam
      integer,save :: calwall=1
!
! --- < 1. Measure Distance for ALL CV .
!
      allocate(LFUTAU(NCV),DISALL(NCV),LBC_SSF(nssfbc),stat=ierr)
      if(ierr/=0) stop 'stop at allocating in list_wall_distance'
      rewall=0
!
! --- list bc
!
      LBC_INDEX=0
      do 150 nb=1,nbcnd
      IBFL=0
      do IBF=1,nssfbc
      NBSSF=LBCSSF(IBF)
      if(nb.eq.NBSSF) then
        IBFL=IBFL+1
      ENDIF
      enddo
      LBC_INDEX(nb)=IBFL
 150  continue
!
      do 120 nb=1,nbcnd
      LBC_INDEX(nb)=LBC_INDEX(nb)+LBC_INDEX(nb-1)
 120  continue
!
      LBC_SSF(:)=0
      DO 600 nb=1,nbcnd
      IBFL=LBC_INDEX(nb-1)
      do IBF=1,nssfbc
      NBSSF=LBCSSF(IBF)
      if(nb.eq.NBSSF) then
        IBFL=IBFL+1
        ICF=IBF+NEDGE
        LBC_SSF(IBFL)=ICF
      ENDIF
      enddo
 600  continue
!-----------------------------------
! --- Measuring distance from wall
!-----------------------------------
      if(icpu.eq.0) then
        fnam='distance'
      else
        my_rank=icpu-1
        HEADER1='wall_hpc'
        call DEFINE_FILE_NAME(DIRHED,HEADER1,my_rank,walldis(icpu))
        fnam=walldis(icpu)
      endif
      LFUTAU(:)=0
      DISALL(:)=great
      if(icpu.eq.0) then
        if(NBCWAL.gt.0) then
          inquire(file=fnam,exist=wallfile)
          if(wallfile) then
 1200       write(ifll,'(2a)')
     &      ' ####    Wall Distance File existed, ',
     &      ' Will you calculate again and over write it ? (n/y)'
            read(*,*) recal
            if(recal=='y') then
              goto 1000
            elseif(recal=='n') then
              goto 2000
            else
              write(ifll,'(a)') ' ####    Input [n] or [y]'
              goto 1200
            endif
          endif
!--------------------------------------
! --- Measure distance to no-slip wall
!--------------------------------------
 1000     continue
          rewall=1
          write(ifll,3010)
          do 700 ICVL=1,NCV    !IIMAT=1,NMAT
          IMAT=LMAT(ICVL)
          if(IMAT.lt.0) goto 700
          if(icaltb.eq.sles.or.cvdist(IMAT)) then
            dismin=GREAT
            dis1=GREAT
            NTAU=0
            do 720 nb=1,nbcnd
            kdv=kdbcnd(1,nb)
            if(kdv.eq.kvnslp.or.kdv.eq.kvlglw.or.kdv==kvrogf.or.
     &         kdv==kvmodl) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_INDEX(nb)
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              dx=SFCENT(1,ICFL)-CVCENT(1,ICVL)
              dy=SFCENT(2,ICFL)-CVCENT(2,ICVL)
              dz=SFCENT(3,ICFL)-CVCENT(3,ICVL)
              dis1=
     &        sqrt((CVCENT(1,ICVL)-SFCENT(1,ICFL))**2
     &            +(CVCENT(2,ICVL)-SFCENT(2,ICFL))**2
     &            +(CVCENT(3,ICVL)-SFCENT(3,ICFL))**2)
! 
              if(dis1.lt.dismin) then
                dismin=dis1
                NTAU=IBFL
                ICFMIN=ICFL
              endif
              enddo
            endif
  720       continue
!
            if(NTAU.eq.0) goto 700
            DISALL(ICVL)=
     &      abs(SFAREA(1,ICFMIN)*(CVCENT(1,ICVL)-SFCENT(1,ICFMIN))
     &         +SFAREA(2,ICFMIN)*(CVCENT(2,ICVL)-SFCENT(2,ICFMIN))
     &         +SFAREA(3,ICFMIN)*(CVCENT(3,ICVL)-SFCENT(3,ICFMIN)))
            LFUTAU(ICVL)=NTAU
            IF(MOD(ICVL,2000).EQ.0) then
              write(ifll,3000) ICVL,dismin,NTAU
            endif
          endif
 700      continue
!
          open (11,FILE=fnam,FORM='UNFORMATTED',
     &          status='unknown',iostat=ios)
          write(11)  NCV
          write(11) (DISALL(ICV),ICV=1,NCV)
          write(11) (LFUTAU(ICV),ICV=1,NCV)
          close(11)
          write(ifll,3020)
        endif
      else
        if(calwall==2) then
!          goto 1100
        elseif(calwall==3) then
          goto 2000
        endif
        if(NBCWAL.gt.0) then
          if(calwall==2) goto 1100
          inquire(file=fnam,exist=wallfile)
          if(wallfile) then
 1210       write(ifll,'(2a)')
     &      ' ####    Wall Distance File existed, ',
     &      ' Will you calculate again and over write it ? (n/y)'
            read(*,*) recal
            if(recal=='y') then
              calwall=2
              goto 1100
            elseif(recal=='n') then
              calwall=3
              goto 2000
            else
              write(ifll,'(a)') ' ####    Input [n] or [y]'
              goto 1210
            endif
          endif
!--------------------------------------
! --- Measure distance to no-slip wall
!--------------------------------------
 1100     continue
          rewall=1
          write(ifll,3010)
          do 701 ICVL=1,NCV    !IIMAT=1,NMAT
          IMAT=LMAT(ICVL)
          if(IMAT.lt.0) goto 701
          if(icaltb.eq.sles.or.cvdist(IMAT)) then
            dismin=GREAT
            dis1=GREAT
            NTAU=0
            do 721 nb=1,nbcnd
            kdv=kdbcnd(1,nb)
            if(kdv.eq.kvnslp.or.kdv.eq.kvlglw.or.kdv==kvrogf.or.
     &         kdv==kvmodl) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_INDEX(nb)
              do IBFL=IBFS,IBFE
              ICFL=LBC_SSF(IBFL)
              dx=SFCENT(1,ICFL)-CVCENT(1,ICVL)
              dy=SFCENT(2,ICFL)-CVCENT(2,ICVL)
              dz=SFCENT(3,ICFL)-CVCENT(3,ICVL)
              dis1=
     &        sqrt((CVCENT(1,ICVL)-SFCENT(1,ICFL))**2
     &            +(CVCENT(2,ICVL)-SFCENT(2,ICFL))**2
     &            +(CVCENT(3,ICVL)-SFCENT(3,ICFL))**2)
! 
              if(dis1.lt.dismin) then
                dismin=dis1
                NTAU=IBFL
                ICFMIN=ICFL
              endif
              enddo
            endif
 721        enddo
!
            if(NTAU.eq.0) goto 701
            DISALL(ICVL)=
     &      abs(SFAREA(1,ICFMIN)*(CVCENT(1,ICVL)-SFCENT(1,ICFMIN))
     &         +SFAREA(2,ICFMIN)*(CVCENT(2,ICVL)-SFCENT(2,ICFMIN))
     &         +SFAREA(3,ICFMIN)*(CVCENT(3,ICVL)-SFCENT(3,ICFMIN)))
            LFUTAU(ICVL)=NTAU
            IF(MOD(ICVL,2000).EQ.0) then
              write(ifll,3000) ICVL,dismin,NTAU
            endif
          endif
 701      enddo
!
          open (11,FILE=fnam,FORM='UNFORMATTED',
     &          status='unknown',iostat=ios)
          write(11)  NCV
          write(11) (DISALL(ICV),ICV=1,NCV)
          write(11) (LFUTAU(ICV),ICV=1,NCV)
          close(11)
          write(ifll,3020)
        endif
      endif
!
 3000 format(12X,'Processing ... NO. ',I8,' FINISHED,   Dis. & BC ',
     &  4X,E10.4,I8,2X)
 3010 format(2X,28('*'),' Start measuring all Distance ',35('*'))
 3020 format(2X,28('*'),' End   measuring all Distance ',35('*'))
!
 2000 continue
!
      deallocate(LFUTAU,DISALL,LBC_SSF)
!
      return
!
      end subroutine list_wall_distance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_nwvrtx
     & (mcell,mface,mvrtx,nvrtx,ncell,nface,ncelb,nvrtxnw,
     &  lvcell,lacell,cord,CVCENT,
     &  ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only          : ifle,ifll,cntlnam
      use module_boundary,only    : nbcnd,ivbcnd,lvbcnd,
     &                              kxnone,kdbcnd,kdsld,kdprdc,
     &                              kdbuff,kdintr,NBOUND,numbname,
     &                              IFFACE,NFBOUN,SFBOUN,boundName,
     &                              idis,dum_bc,
     &                              kdshutr,kdpors,kdovst
      use module_partitioner,only : ino  =>WK1
      use module_partitioner,only : ivn  =>WK2
      use module_partitioner,only : ipv  =>WK3
      use module_partitioner,only : ibc  =>WK4
      use module_partitioner,only : IW1
      use module_partitioner,only : IBCF=>ISUF_NODE_O !ISUF_NODE_P
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mcell,mface,mvrtx
      integer,intent(inout)  :: nvrtx,nface,nvrtxnw
      integer,intent(in)     :: ncell,ncelb
      integer,intent(in)     :: lacell(  mcell)
      integer,intent(inout)  :: lvcell(8,mcell)
!      integer,intent(inout)  :: lfcell(7,mcell)
!      integer,intent(inout)  :: lvface(4,mface)
!      integer,intent(inout)  :: lcface(2,mface)
!      integer,intent(inout)  :: LVRTCV(  mvrtx)
!      integer,intent(inout)  :: LCVFAC(4,mface)
      real*8 ,intent(inout)  :: cord(3,mvrtx),CVCENT(3,mvrtx)
      integer,intent(out)    :: ierror
!
! --- [local entities]
!
      integer,parameter :: M1_MAT=-100,M2_MAT=150
      integer,parameter :: M12=M2_MAT+M1_MAT
!
      real*8,parameter :: SML=1.d-25,zero=0.d0
      integer          :: IC,IS,IV,ICVA,ICVB,ierr1=0
      integer          :: nb,mb,id,kd,ie,iv1,iv2,IVV,ILS,inv
      integer          :: i,j,jj,Int,ibcchk,NBCFAC,ibcnode
      integer :: IMAT,NMAT,IMAT1,IIMAT
      integer :: MMAT(M1_MAT:M2_MAT),MAT_NO(M12),nvtmat(0:M12)
      integer,allocatable :: MAT_BC(:),nbc(:),lflag(:)
      integer,allocatable :: namflg(:,:),llflag(:,:)
!
      logical :: logsld=.false.,logbuf=.false.,logdum,logshu=.false.
      
!--------------------------------
! --- 
!--------------------------------
      do nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if((kd==kdsld ).or.
     &   (kd==kdovst).or.
     &   (kd==kdintr).or.
     &   (kd==kdbuff).or.
!     &   (kd==kdpors).or.
     &   (kd==kdshutr)
!     &   (kd==kdprdc)
     &   ) then
         logsld=.true.
         if(kd==kdshutr) then
           logshu=.true.
         endif
!      elseif(kd==kdbuff) then
!         logbuf=.true.
      endif
      enddo
      if(.not.logsld.and..not.logbuf) return
      if(logbuf.and..not.logsld) goto 2000
!
      allocate(ino(mvrtx),stat=ierr1)
      if(ierr1.ne.0) stop 'allocate ino in list_nwvrtx'
      allocate(ipv(0:mvrtx),stat=ierr1)
      if(ierr1.ne.0) stop 'allocate ipv in list_nwvrtx'
      allocate(ivn(0:mvrtx),stat=ierr1)
      if(ierr1.ne.0) stop 'allocate ivn in list_nwvrtx'
!      allocate(IW1(8,mvrtx),stat=ierr1)
      allocate(IW1(8,mcell),stat=ierr1)
      if(ierr1.ne.0) stop 'allocate IW1 in list_nwvrtx'
      allocate(MAT_BC(NBOUND),nbc(0:NBOUND),lflag(NBOUND),stat=ierr1)
      if(ierr1.ne.0) stop 'allocate MAT_BC in list_nwvrtx'
      allocate(namflg(nbcnd,numbname),llflag(nbcnd,numbname),
     &                                                    stat=ierr1)
      if(ierr1.ne.0) stop 'allocate namflg in list_nwvrtx'
!-------------------------------------------------------------
! --- CHECK BC NAME between fflow.ctl and BC name in BC-file
!-------------------------------------------------------------
      namflg(:,:)=0
      lflag(:)=0
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j)/=' ') then
        do mb=1,NBOUND
        logdum=trim(adjustl(boundName(nb,j)))==trim(adjustl(SFBOUN(mb)))
        logdum=logdum.or.trim(adjustl(boundName(nb,j)))=='undefined'
        logdum=logdum.or.trim(adjustl(boundName(nb,j)))=='dummy'  !zhanghuilai
        if(logdum) then
          lflag(mb)=1
          namflg(nb,j)=1
        endif
        enddo
      endif
      enddo
      enddo
!.........................................
      do mb=1,NBOUND
      if(lflag(mb)==0) then
        write(ifle,'(2a)') 
     &  'ERR-1: NOT Defined BC name in fflow.ctl: ',
     &  trim(adjustl(SFBOUN(mb)))
        write(ifle,'(a)') 
     & 'ATT: Maybe your data transport is wrong between PC <=> Unix'
        write(ifle,'(a)') 
     & 'ATT: Be sure using: ASSIC or Binary form',
     & ' while transfer mesh file'
        stop
      endif
      enddo
!.........................................
      do j=1,numbname
      do nb=1,nbcnd
      if(namflg(nb,j)==0.and.boundName(nb,j)/=' ') then
        write(ifle,'(2a)') 
     &  'ERR: NOT Found BC name in user-BC grid file: ',
     &  trim(adjustl(boundName(nb,j)))
        stop
      endif
      enddo
      enddo
!.........................................
      namflg(:,:)=-1
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j).ne.' ') then
        do mb=1,NBOUND
        if(trim(adjustl(boundName(nb,j))).eq.trim(adjustl(SFBOUN(mb))))
     &  then
          namflg(nb,j)=mb
          exit
        endif
	enddo
      endif
      enddo
      enddo
!------------------------------------
! --- 
!------------------------------------
      NBCFAC=0
      do mb=1,NBOUND
        NBCFAC=NBCFAC+NFBOUN(mb)
      enddo
      allocate(IBCF(4,NBCFAC,NBOUND),stat=ierr1)
      if(ierr1.ne.0) stop 'allocate IBCF'
!
      IBCF(:,:,:)=0
      do mb=1,NBOUND
      do IS=1,NFBOUN(mb)
      do Int=1,4
        IBCF(Int,IS,mb)=IFFACE(Int,IS,mb)
      enddo
      enddo
      enddo
!
      do iv=1,nvrtx
      do i=1,3
      CVCENT(i,iv)=cord(i,iv)
      enddo
      enddo
!
!--------------------
! --- List materials
!--------------------
!
      MMAT(:)=0
      do 100 IC=1,ncell
      IMAT=lacell(IC)
      MMAT(IMAT)=1
 100  enddo
!
      NMAT=0
      MAT_NO(:)=0
      do 110 IMAT=M2_MAT,M1_MAT,-1
      if(MMAT(IMAT).eq.1) then
        NMAT=NMAT+1
        MAT_NO(NMAT)=IMAT
      endif
 110  enddo
!
      if(NMAT.lt.1) then
        write(ifle,*) 'ERR: At last one material NO. has been defined'
        stop
      ELSEif(NMAT.eq.1.and.logsld) then
        write(ifle,*) 'ERR: Interface BC , Sliding BC or Buffle BC',
     &       'Both of them are BCs between different materials'
        stop
      endif
!
!-----------------
! --- List BC face
!-----------------
!
      nbc(:)=0
      ipv(0)=1
      do 200 mb=1,NBOUND
!
      ipv(:)=0
      do IS=1,NFBOUN(mb)
      do Int=1,4
      IV=IFFACE(Int,IS,mb)
      ipv(IV)=1
      enddo
      enddo
!
      nbc(mb)=nbc(mb-1)
      do iv=1,nvrtx
      if(ipv(IV).gt.0) then
        nbc(mb)=nbc(mb)+1
      endif
      enddo
!
 200  enddo
!
      ALLOCATE(ibc(nbc(NBOUND)))
!
      nbc(:)=0
      do 210 mb=1,NBOUND
!
      ipv(:)=0
      do IS=1,NFBOUN(mb)
      do Int=1,4
      IV=IFFACE(Int,IS,mb)
      ipv(IV)=1
      enddo
      enddo
!
      nbc(mb)=nbc(mb-1)
      do iv=1,nvrtx
      if(ipv(IV).gt.0) then
        nbc(mb)=nbc(mb)+1
        ibc(nbc(mb))=iv
      endif
      enddo
!
 210  enddo
!-----------------------
! --- Creat new vertex
!-----------------------
      ino(:)=0
      nvtmat(:)=0
      IW1(:,:)=0
      llflag(:,:)=0
      lflag(:)=0
      ipv(:)=0
!
      do 300 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ipv(:)=0
      do 310 IC=1,ncell
      IMAT1=lacell(IC)
      if(IMAT1==IMAT) then
        do IVV=1,8
        IV=lvcell(IVV,IC)
        ipv(IV)=ipv(IV)+1
        enddo
      endif
 310  enddo
!
      nvtmat(IIMAT)=nvtmat(IIMAT-1)
!
      ivn(:)=0
      do 320 iv=1,nvrtx
      if(ipv(IV).gt.0) then
        nvtmat(IIMAT)=nvtmat(IIMAT)+1
        ivn(iv)=nvtmat(IIMAT)
        if(nvtmat(IIMAT).gt.mvrtx) then
          write(ifle,'(1X,2a)') 'ERR : Creat mvrtx in ',cntlnam
        endif
        ino(nvtmat(IIMAT))=iv
!        imt(nvtmat(IIMAT))=IMAT1
      endif
 320  enddo
!
      write(ifll,'(1x,2(a,I12))') 'MSG: Vertex NO=',
     &  nvtmat(IIMAT)-nvtmat(IIMAT-1),
     & ' in MAT=',IMAT
!
      ivn(0)=0
      do 330 IC=1,ncell
      IMAT1=lacell(IC)
      if(IMAT1==IMAT) then
        do IVV=1,8
        iv=lvcell(IVV,IC)
        IW1(IVV,IC)=ivn(iv)
        
        enddo
      endif
 330  continue
!
      do 340 mb=1,NBOUND
!
      ibcchk=0
      do id=nbc(mb-1)+1,nbc(mb)
      iv=ibc(id)
      if(ipv(iv).gt.0) then
        ibcchk=ibcchk+1
      endif
      enddo
!
      do IS=1,NFBOUN(mb)
        do Int=1,4
        iv=IBCF(Int,IS,mb)
        if(ivn(iv)>0) then
          ipv(IV)=2
        endif
        enddo
      enddo
      ibcnode=0
      do id=nbc(mb-1)+1,nbc(mb)
      iv=ibc(id)
      if(ipv(iv)==2) then
        ibcnode=ibcnode+1
      endif
      enddo
!
      if(ibcnode==ibcchk.and.
     &   ibcchk/=0.and.
     &   ibcnode==(nbc(mb)-nbc(mb-1)).and.
     &   lflag(mb)==0
     &   ) then
        if(dum_bc/=0.and.mb==NBOUND.and.lflag(NBOUND)/=IMAT) then
          do IS=1,NFBOUN(mb)
          do Int=1,4
          IV=IBCF(Int,IS,mb)
          if(ivn(IV)>0) then
            IFFACE(Int,IS,mb)=ivn(iv)
          endif
          enddo
          enddo
        else
          lflag(mb)=IMAT
          do IS=1,NFBOUN(mb)
          do Int=1,4
          IV=IBCF(Int,IS,mb)
          if(ivn(IV)>0) then
            IFFACE(Int,IS,mb)=ivn(iv)
          endif
          enddo
          enddo
        endif
      endif
!
 340  enddo
!
 300  continue
!
      do 400 IC=1,ncell
      do IVV=1,8
      lvcell(IVV,IC)=IW1(IVV,IC)
      enddo
 400  continue
!
      nvrtxnw=nvtmat(NMAT)
      write(ifll,*) 'MSG: Orginal Vertex number: ', nvrtx
      write(ifll,*) 'MSG: Current Vertex number: ', nvrtxnw
      write(ifll,*) 'MSG: Created Vertex number: ', nvrtxnw-nvrtx
      nvrtx=nvrtxnw
      if(nvrtxnw.gt.mvrtx) then
        write(ifle,'(1x,a,I12)') 'ERR : Creat mvrtx > '  ,nvrtxnw
        stop
      endif
!
      do iv=1,nvrtxnw
      do i=1,3
      cord(i,iv)=CVCENT(i,ino(iv))
      enddo
      enddo
!
      deallocate(ino)
      deallocate(ipv)
      deallocate(ivn)
      deallocate(ibc)
      deallocate(IBCF)
      deallocate(IW1)
!
      deallocate(MAT_BC,nbc,lflag,namflg,llflag)
      CVCENT(:,:)=0.d0
!
!-------------------------------------------------------------------
! --- Creat buffle layer new vertex
!-------------------------------------------------------------------
 2000 continue
!
!      NBCFAC=0
!      do mb=1,NBOUND
!        NBCFAC=NBCFAC+NFBOUN(mb)
!      enddo
!      allocate(IBCF(4,NBCFAC,NBOUND),stat=ierr1)
!      if(ierr1.ne.0) stop 'allocate IBCF in buffle layer new vertex'
!      IBCF(:,:,:)=IFFACE(:,:,:)
!      CVCENT(:,:)=cord(:,:)
!
!      do nb=1,nbcnd
!      kd=kdbcnd(0,nb)
!      if(kd==kdbuff) then
!        do mb=1,NBOUND
!        if(trim(adjustl(boundName(nb,1)))==trim(adjustl(SFBOUN(mb))))
!     &  then
!          do IS=1,NFBOUN(mb)
!            do Int=1,4
!            IV=IBCF(Int,IS,mb)
!            enddo
!          endif
!        endif
!        enddo
!      endif
!      enddo
!
!      deallocate(IBCF)
!
!-------------------------------------------------------------------
!
      ierror=0
      return
!
      end subroutine list_nwvrtx
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_output_sld
     & (NCPU,nssfbc,nedge,mssfbc,medge,LBCSSF,LCYCSF)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ////////////////  NOT USED    ///////////////
!
!     WK3(iisd)=isd
!     WK2(isdn)=isd
!     WK4(isdnn)=isd
!     INDXSLD(isd)=isdnn                       
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_io,only           : ifle,ifll
      use module_material,only     : ical_sld,rotati,ishaft,end,begin,
     &                               nsplit
      use module_hpc,only          : nsldT,nsldbc,NPETOT,nmsldL,nmsldG,
     &                               nbsldG,MAT_BCSLD
      use module_boundary,only     : nbcnd,LBC_INDEX,kdbcnd,
     &                               kdsld,idis
      use module_partitioner,only  : LSLD   => IW1
      use module_cood       ,only  : CSLD   => gfacesq
      use module_partitioner,only  : INDXSLD=> WK1
      use module_cood       ,only  : CSLDG  => volumesq
      use module_cood       ,only  : cylcod => area
      use module_partitioner,only  : WK2,WK3,WK4
      use module_param ,only   : scalfctr,smallNrm
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: nssfbc,nedge,mssfbc,medge,NCPU
      integer,intent(in) :: LBCSSF( mssfbc)
      integer,intent(in) :: LCYCSF( mssfbc)
!
! --- [local entities]
!
      integer :: my_rank,ifl,icpu,icpu1,ifchl,kmax1,kmax2,k,ios      
      integer :: nb,IBFL,IBF,NBSSF,kd,IBFS,IBFE,IBFP
      integer :: nbsld,i,j,isd,nt
      integer :: nbsldx,nsldx(0:nbcnd)
      integer :: nbsldGx,nsldG(0:nbcnd)
      integer :: nslid(0:NPETOT)
      integer :: nbx,ical_sldx,IBFSx,IBFEx,ierr1=0
      integer :: iisd,nbG,isdG,iisdG,nisd,isdn,NSLDTOT
      integer :: MAT_BCSLDL(nbcnd,2),IMAT1,IMAT2,nslp,ns
      integer :: NEWSLD(nbcnd),NMSLD(nbcnd),isdnn,nbx1
      real*8  :: dx,dy,dz,dx1,dy1,dz1,dl,dmin,x(3),x1(3)
      real*8  :: unit1(3),unit2(3),rpm,alpha,dalp,sgn
      real*8,parameter :: SML=1.d-15
      real*8  :: dt,ttt
!
      if(ical_sld/=1) then
        return
      endif
      deallocate(cylcod)
!-------------------------------------------------------------------
! --- 
!-------------------------------------------------------------------
      if(NCPU==1) then
        open (11,FILE='sliding',FORM='UNFORMATTED',
     &     status='unknown',iostat=ios)
        if(ios/=0) stop 'write sliding at list_output_sld'
! --- 
        nbsldGx=0
        nsldG(:)=0
        read(11) nbsldGx,(nsldG(nbx),nbx=1,nbsldGx)
        read(11) ((MAT_BCSLDL(nbx,i),i=1,2),nbx=1,nbsldGx)
!
        do nbx=1,nbsldGx
        nsldG(nbx)=nsldG(nbx-1)+nsldG(nbx)
        enddo
!
        IF(nbsldGx/=nbsldG) THEN
          WRITE(IFLE,*) 'ERR: Sliding Mesh reading error'
          stop 'list_output_sld'
        endif
        allocate(LSLD(nsldG(nbsldG),2),stat=ierr1)
        allocate(CSLD(nsldG(nbsldG),3,2),stat=ierr1)
        allocate(INDXSLD(nsldG(nbsldG)),stat=ierr1)
        allocate(CSLDG(nsldG(nbsldG),3),stat=ierr1)
        allocate(cylcod(nsldG(nbsldG),2),stat=ierr1)
        if(ierr1.ne.0) stop ': at allocating LSLD'
        LSLD(:,:)=0
        CSLD(:,:,:)=0.d0
        INDXSLD(:)=0
        CSLDG(:,:)=0.d0
!
        do nbx=1,nbsldG
        read(11)  (LSLD(isd,1),LSLD(isd,2),
     &             isd=nsldG(nbx-1)+1,nsldG(nbx))
        read(11) ((CSLD(isd,j,1),j=1,3),
     &             isd=nsldG(nbx-1)+1,nsldG(nbx))
        read(11) ((CSLD(isd,j,2),j=1,3),
     &             isd=nsldG(nbx-1)+1,nsldG(nbx))
        enddo
        close(11)
        CSLDG(:,:)=CSLD(:,:,1)
      endif
!--------------------------------------------------------------------
! --- 
!--------------------------------------------------------------------
      if(NCPU.gt.1) then
! --- 
        nslid(:)=0
        do icpu=1,NCPU
        do nbx=1,nsldbc(icpu)
        nslid(icpu)=nslid(icpu)+nsldT(nbx,icpu)
        enddo
        enddo
        do icpu=1,NCPU
        nslid(icpu)=nslid(icpu-1)+nslid(icpu)
        enddo
        NSLDTOT=nslid(NCPU)
        allocate(LSLD(NSLDTOT,2),stat=ierr1)
        allocate(CSLD(NSLDTOT,3,2),stat=ierr1)
        allocate(INDXSLD(NSLDTOT),stat=ierr1)
        allocate(CSLDG(NSLDTOT,3),stat=ierr1)
        allocate(cylcod(NSLDTOT,2),stat=ierr1)
        LSLD(:,:)=0
        CSLD(:,:,:)=0.d0
        INDXSLD(:)=0
        CSLDG(:,:)=0.d0
! --- 
        do 100 icpu=1,NCPU
        ifl=50
        open(ifl,file=SLIDNG(icpu),FORM='unformatted',
     &       status='unknown',iostat=ios)
        nbsldx=0
        nsldx(:)=0
        read(ifl) nbsldx,(nsldx(nbx),nbx=1,nbsldx)
        read(ifl) ((MAT_BCSLDL(nbx,i),i=1,2),nbx=1,nbsldx)
        do nbx=1,nbsldx
        nsldx(nbx)=nsldx(nbx-1)+nsldx(nbx)
        enddo
        do nbx=1,nbsldx
        read(ifl) (LSLD(isd,1),LSLD(isd,2),
     &   isd=nslid(icpu-1)+nsldx(nbx-1)+1,nsldx(nbx)+nslid(icpu-1))
        read(ifl) ((CSLD(isd,j,1),j=1,3),
     &   isd=nslid(icpu-1)+nsldx(nbx-1)+1,nsldx(nbx)+nslid(icpu-1))
        read(ifl) ((CSLD(isd,j,2),j=1,3),
     &   isd=nslid(icpu-1)+nsldx(nbx-1)+1,nsldx(nbx)+nslid(icpu-1))
        enddo
        close(ifl)
 100    continue
!
! --- 
!
        allocate(WK2(NSLDTOT),WK3(NSLDTOT),WK4(NSLDTOT),stat=ierr1)
!
        do 200 nbG=1,nbsldG
        IMAT1=MAT_BCSLD(nbG,1)
        IMAT2=MAT_BCSLD(nbG,2)
        INDXSLD(:)=0
        WK2(:)=0
        WK3(:)=0
        WK4(:)=0
        NMSLD(nbG)=0
        do 210 icpu=1,ncpu
        do nbx=1,nsldbc(icpu)
        if(trim(nmsldL(nbx,icpu)).eq.trim(nmsldG(nbG))) then
          do isd=nslid(icpu-1)+nsldT(nbx-1,icpu)+1,
     &           nslid(icpu-1)+nsldT(nbx,icpu)
          NMSLD(nbG)=NMSLD(nbG)+1
          WK2(NMSLD(nbG))=isd
          x(1)=CSLD(isd,1,1)
          x(2)=CSLD(isd,2,1)
          x(3)=CSLD(isd,3,1)
          DO 220 icpu1=icpu+1,ncpu
            do nbx1=1,nsldbc(icpu1)
            IF(trim(nmsldL(nbx1,icpu1)).eq.trim(nmsldG(nbG))) THEN
              do iisd=nslid(icpu1-1)+nsldT(nbx1-1,icpu1)+1,
     &                nslid(icpu1-1)+nsldT(nbx1,icpu1)
              x1(1)=CSLD(iisd,1,1)
              x1(2)=CSLD(iisd,2,1)
              x1(3)=CSLD(iisd,3,1)
              dx=x1(1)-x(1)
              dy=x1(2)-x(2)
              dz=x1(3)-x(3)
              dl=dx*dx+dy*dy+dz*dz
              if(dl.lt.smallNrm*scalfctr(nbG)) then
                WK3(isd)=iisd
              ENDIF
              enddo
            endif
            enddo
 220      continue
          enddo
        endif
        enddo
 210    continue
! --- 
        NEWSLD(nbG)=0
        DO isdn=1,NMSLD(nbG)
        isd=WK2(isdn)
        if(WK3(isd)==0) then
          NEWSLD(nbG)=NEWSLD(nbG)+1
          WK4(NEWSLD(nbG))=isd
          INDXSLD(isd)=NEWSLD(nbG)
        endif
        enddo
!
        DO isdn=1,NMSLD(nbG)
        isd=WK2(isdn)
        if(WK3(isd)/=0) then
          iisd=WK3(isd)
          isdnn=INDXSLD(iisd)
          INDXSLD(iisd)=isdnn
        endif
        enddo
!
        DO isdn=1,NMSLD(nbG)
        isd=WK2(isdn)
        iisd=WK3(isd)
        if(INDXSLD(iisd)==0) then
          write(ifle,*)
     &    'ERR: Sliding Bc Marging error, Connect your supportor'
          stop ': List_output_sld'
        endif
        enddo
!
        if(ishaft(IMAT1)==1) then
          do isdnn=1,NEWSLD(nbG)
          isd=WK4(isdnn)
          CSLDG(isdnn,:)=CSLD(isd,:,1)
          enddo
!
          do isdnn=1,NEWSLD(nbG)
          cylcod(isdnn,1)=CSLDG(isdnn,1)
          cylcod(isdnn,2)=CSLDG(isdnn,2)
          enddo
!
          unit1(1)=end(1,IMAT1)-begin(1,IMAT1)
          unit1(2)=end(2,IMAT1)-begin(2,IMAT1)
          unit1(3)=end(3,IMAT1)-begin(3,IMAT1)
          dx=unit1(1)*unit1(1)
     &      +unit1(2)*unit1(2)
     &      +unit1(3)*unit1(3)
          dx=dsqrt(dx)
          unit1(1)=unit1(1)/dx
          unit1(2)=unit1(2)/dx
          unit1(3)=unit1(3)/dx
          x(1)=0.d0    ! only support z-axis
          x(2)=0.d0
          x(3)=1.d0
          rpm=rotati(IMAT1)
          nslp=nsplit(IMAT1)
          sgn=
     &    sign(1.d0,(x(1)*unit1(1)+x(2)*unit1(2)+x(3)*unit1(3))*rpm)
          dalp=sgn*2.d0*3.1415926/dble(nslp)
          DO 230 ns=0,nslp
          do isdnn=1,NEWSLD(nbG)
            CSLDG(isdnn,1)=cylcod(isdnn,1)*cos(alpha)
     &                    -cylcod(isdnn,2)*sin(alpha)
            CSLDG(isdnn,2)=cylcod(isdnn,1)*sin(alpha)
     &                    +cylcod(isdnn,2)*cos(alpha)
          enddo
          alpha=alpha+dalp
! --- only z-axis
          do isdn=1,NMSLD(nbG)
          isd=WK2(isdn)
          isdnn=INDXSLD(isd)
          CSLD(isd,1,1)=CSLDG(isdnn,1)
          CSLD(isd,2,1)=CSLDG(isdnn,2)
          enddo
 230      continue
        endif
        if(ishaft(IMAT2)==1) then
          do isdnn=1,NEWSLD(nbG)
          isd=WK4(isdnn)
          CSLDG(isdnn,:)=CSLD(isd,:,2)
          enddo
!
          do isdnn=1,NEWSLD(nbG)
          cylcod(isdnn,1)=CSLDG(isdnn,1)
          cylcod(isdnn,2)=CSLDG(isdnn,2)
          enddo
!
          unit2(1)=end(1,IMAT2)-begin(1,IMAT2)
          unit2(2)=end(2,IMAT2)-begin(2,IMAT2)
          unit2(3)=end(3,IMAT2)-begin(3,IMAT2)
          dx1=unit2(1)*unit2(1)
     &       +unit2(2)*unit2(2)
     &       +unit2(3)*unit2(3)
          dx1=dsqrt(dx1)
          unit2(1)=unit2(1)/dx1
          unit2(2)=unit2(2)/dx1
          unit2(3)=unit2(3)/dx1
          x(1)=0.d0    ! only support z-axis
          x(2)=0.d0
          x(3)=1.d0
          rpm=rotati(IMAT2)
          nslp=nsplit(IMAT2)
          sgn=
     &    sign(1.d0,(x(1)*unit2(1)+x(2)*unit2(2)+x(3)*unit2(3))*rpm)
          dalp=sgn*2.d0*3.1415926/dble(nslp)
          alpha=0.d0
!
          DO 240 ns=0,nslp
          do isdnn=1,NEWSLD(nbG)
            CSLDG(isdnn,1)=cylcod(isdnn,1)*cos(alpha)
     &                    -cylcod(isdnn,2)*sin(alpha)
            CSLDG(isdnn,2)=cylcod(isdnn,1)*sin(alpha)
     &                    +cylcod(isdnn,2)*cos(alpha)
          enddo
          alpha=alpha+dalp
! --- only z-axis
          do isdn=1,NMSLD(nbG)
          isd=WK2(isdn)
          isdnn=INDXSLD(isd)
          CSLD(isd,1,1)=CSLDG(isdnn,1)
          CSLD(isd,2,1)=CSLDG(isdnn,2)
          enddo
 240      continue
        endif
!
 200    continue
!
        deallocate(WK2,WK3,WK4)
      endif
!
      deallocate(LSLD,CSLD,CSLDG,INDXSLD)
      return
!
      end subroutine list_output_sld
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_output_grid
     & (mvrtx,mcell,mface,mssfbc,nvrtx,ncell,nface,nssfbc,nssfbc_old,
     &  cord,
     &  lacell,lvcell,lvface,lbface,
     &  lcface,lfcell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_boundary,only : SFBOUN,NBOUND,boundIDMap,numbname,
     &                           IFFACE,NFBOUN,NAME_I
      use module_io,only       : gdScale
      use module_model,only    : ical_mvmsh
      use module_material,only : nflud,nofld,nsold,nosld
      use module_hpc,     only   : NPETOT
      use module_model,only       : vertex_cen,cell_cen,icon_cv
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: mvrtx,mcell,mface,mssfbc,nssfbc_old,nssfbc
      integer,intent(in) :: nvrtx,ncell,nface
      real*8 ,intent(inout) :: cord(3,mvrtx)
!
      integer,intent(inout) :: lacell(  mcell)
      integer,intent(in)    :: lvcell(8,mcell)
      integer,intent(in) :: lvface(4,mface)
      integer,intent(in) :: lbface(2,mface)
      integer,intent(in) :: lcface(2,mface)
      integer,intent(in) :: lfcell(7,mcell)
!      INTEGER,INTENT(in) :: LEFACE(5,MFACE)
!      INTEGER,INTENT(in) :: listbc(4,mssfbc)
!
! --- [local entities]
!
      integer :: i,j,IS,boundID,my_rank,iv,ic,mb
      integer :: ifl,ios=0
      integer :: NBFS,NBOUND_TEMP,ifld
      character(len=80),save :: fnam
      integer,allocatable :: IBFACE(:)
      integer,allocatable :: IBIDX(:)
!
! --- output grid file
!
      if(.true.) then

      allocate(IBIDX(0:NBOUND))
      NBFS=0
      IBIDX=0
      do mb=NAME_I,NBOUND 
      NBFS=NBFS+NFBOUN(mb) 
      IBIDX(mb)=NFBOUN(mb) 
      enddo 
!
      do mb=NAME_I,NBOUND 
        IBIDX(mb)=IBIDX(mb)+IBIDX(mb-1) 
      enddo

      else

      NBFS=0
      allocate(IBIDX(0:NBOUND))      
      IBIDX=0
      do mb=NAME_I,NBOUND 
      NBFS=NBFS+NFBOUN(mb) 
      IBIDX(mb)=NFBOUN(mb) 
      enddo 
!
      do mb=NAME_I,NBOUND 
        IBIDX(mb)=IBIDX(mb)+IBIDX(mb-1) 
      enddo
      endif      
!
      fnam='geom.frontflow'
      ifl=50
      open(ifl,file=fnam,FORM='unformatted',status='unknown',
     &            iostat=ios)
      if(ios/=0) then
        write(*,*)'*** Cannot create FrontFlowRed Grid File:',fnam
        return
      endif
!
      NBOUND_TEMP=NBOUND+(1-NAME_I)
!
      if(NPETOT.eq.1) then
        do ifld=1,nflud
        do IC=1,ncell
        if(lacell(IC)==ifld) then
        lacell(IC)=nofld(ifld)
        endif
        enddo
        enddo
      endif
!
      write(ifl)   nvrtx,ical_mvmsh,icon_cv
      write(ifl) ((cord(i,iv)/gdScale, i=1,3), iv=1,nvrtx)
      write(ifl)   ncell
      write(ifl)  (lacell(ic),ic=1,ncell)
      write(ifl) ((lvcell(i,ic), i=1,8), ic=1,ncell)
      write(ifl)   NBOUND_TEMP
      write(ifl)   NBFS
!
      do mb=NAME_I,NBOUND
        NBOUND_TEMP=mb+(1-NAME_I)
        write(ifl) NBOUND_TEMP,NFBOUN(mb),SFBOUN(mb)
        write(ifl) mb,IBIDX(mb)
        write(ifl) IFFACE(1:4,1:NFBOUN(mb),mb)
      enddo
!
      
!
      close(ifl)
!
! --- 
!
      return
      end subroutine list_output_grid


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_output_ini_movegrid
     & (mvrtx,mcell,mface,mssfbc,nvrtx,ncell,nface,nssfbc,nssfbc_old,
     &  cord,
     &  lacell,lvcell,lvface,lbface,
     &  lcface,lfcell,LEFACE,listbc
     &  )
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_boundary,only : SFBOUN,NBOUND,boundIDMap,numbname,
     &                           IFFACE,NFBOUN,NAME_I
      use module_io,only       : gdScale
      use module_model,only    : ical_mvmsh
      use module_model,   only : vertex_cen,cell_cen,icon_cv
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: mvrtx,mcell,mface,mssfbc
      integer,intent(in) :: nvrtx,ncell,nface,nssfbc,nssfbc_old
      real*8 ,intent(inout) :: cord(3,mvrtx)
!
      integer,intent(in) :: lacell(  mcell)
      integer,intent(in) :: lvcell(8,mcell)
      integer,intent(in) :: lvface(4,mface)
      integer,intent(in) :: lbface(2,mface)
      integer,intent(in) :: lcface(2,mface)
      integer,intent(in) :: lfcell(7,mcell)
      INTEGER,INTENT(in) :: LEFACE(5,MFACE)
      INTEGER,INTENT(in) :: listbc(4,mssfbc)
!
! --- [local entities]
!
      integer :: i,j,IS,boundID,my_rank,iv,ic,mb
      integer :: ifl,ios=0
      integer :: NBFS,NBOUND_TEMP
      character(len=80),save :: fnam
!
! --- 
!
      NBFS=0
      do mb=NAME_I,NBOUND
      NBFS=NBFS+NFBOUN(mb)
      enddo
!
! --- output grid file
!
      fnam='moveinigrid.frontflow'
      ifl=50
      open(ifl,file=fnam,FORM='unformatted',status='unknown',
     &            iostat=ios)
      if(ios/=0) then
        write(*,*)'*** Cannot create FrontFlowRed Grid File:',fnam
        return
      end if
!
      NBOUND_TEMP=NBOUND+(1-NAME_I)
!
      write(ifl)   nvrtx,ical_mvmsh,icon_cv
      write(ifl) ((cord(i,iv)/gdScale, i=1,3), iv=1,nvrtx)
      write(ifl)   ncell
      write(ifl)  (lacell(ic),ic=1,ncell)
      write(ifl) ((lvcell(i,ic), i=1,8), ic=1,ncell)
      write(ifl)   NBOUND_TEMP
      write(ifl)   NBFS
!
      do mb=NAME_I,NBOUND
        NBOUND_TEMP=mb+(1-NAME_I)
        write(ifl) NBOUND_TEMP,NFBOUN(mb),SFBOUN(mb)
        write(ifl) IFFACE(1:4,1:NFBOUN(mb),mb)
      end do
!
      write(ifl) ((lvface(i,iv),i=1,4),iv=1,nface)
      write(ifl) ((lbface(i,iv),i=1,2),iv=1,nface)
      write(ifl) ((lcface(i,iv),i=1,2),iv=1,nface)
      write(ifl) ((lfcell(i,iv),i=1,7),iv=1,ncell)
      write(ifl) ((LEFACE(i,iv),i=1,5),iv=1,nface)
      write(ifl) ((listbc(i,iv),i=1,4),iv=1,nssfbc_old)
!
      close(ifl)
!
! --- 
!
      return
      end subroutine list_output_ini_movegrid
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_output_grid_A
     & (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_FFRGFwords
      use module_io,only       : getfil,gdScale
      use module_boundary,only : SFBOUN,NBOUND,boundIDMap,numbname,
     &                           IFFACE,NFBOUN
!---------------------------------------------------------------------
!    SFBOUN : Boundary region name read from grid file.
!    NBOUND : Number of boundary regions read from grid file
!    boundIDMap : Indexmap from FFR-Internal boundary region ID to 
!                 grid file boundary ID
!    numbname : Maximum number of bounary region for one boundary 
!               condition (usually 2)
!    IFFACE : Vertices ID constructing boundary face.
!    NFBOUN : Number of faces in a boundary region read form grid file.
!---------------------------------------------------------------------
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: mvrtx,mcell
      integer,intent(in) :: nvrtx,ncell
      real*8 ,intent(in) :: cord(3,mvrtx)
      integer,intent(in) :: lacell(  mcell)
      integer,intent(in) :: lvcell(8,mcell)
! --- [local entities]
      character(len=80),save :: fnam
      integer :: ifl,ios=0
      integer :: numofMaterial
      integer,allocatable :: material_ID(:)
      character*80 :: tmpStr80_1,tmpStr80_2
      integer :: ic,jc

      call getfil(ifl,fnam,'ffrgrid')
      open(ifl,file=fnam,form='formatted',
     &  status='unknown',iostat=ios)
      if(ios/=0) then
        write(*,*)'*** Cannot create FrontFlowRed Grid File:',
     &                                                trim(fnam)
        return
      end if

!     preparation of grid data...
      call ffr_gr_prepData(ios)
      if(ios/=0) then
        write(*,*)'*** Error in preparation of grid data.'
        return
      end if
      
!     write FFR-GRID header field...
      call ffr_gr_headerout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid header.'
        return
      end if

!     write FFR-GRID element type field...
      call ffr_gr_elemtypeout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid element type.'
        return
      end if

!     write FFR-GRID boundary type field...
      call ffr_gr_bndryTypeout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid boundary type.'
        return
      end if

!     write FFR-GRID vertices coordinates field...
      call ffr_gr_vrtxCordout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid vertices coordinates.'
        return
      end if

!     write FFR-GRID boundary face list field...
      call ffr_gr_bndryFaceout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid boundary faces.'
        return
      end if

!     write FFR-GRID element vertices list field...
      call ffr_gr_elemVrtxout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid elements.'
        return
      end if
      

!     termination of FFR-GRID file output
      write(ifl)fileend_FFGF
      close(ifl)
!
      deallocate(IFFACE,NFBOUN)
!      
      return
      
      contains
!===============================================================      
      subroutine ffr_gr_prepData(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
!      count up material number
      numofMaterial=maxval(lacell(1:ncell))-minval(lacell(1:ncell))+1
      do ic=minval(lacell(1:ncell)),maxval(lacell(1:ncell))
        if(count(lacell(1:ncell)==ic)==0) numofMaterial=numofMaterial-1
      end do
!
! --- make material ID list
!
      allocate(material_ID(1:numofMaterial))
      jc=1
      do ic=minval(lacell(1:ncell)),maxval(lacell(1:ncell))
        if(count(lacell(1:ncell)==ic)/=0) then
          material_ID(jc)=ic; jc=jc+1
        end if
      end do
      
      statRtrn=0
      end subroutine ffr_gr_prepData
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_headerout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
 1010 format(a80)
 1020 format(1I8)
 1030 format(I8,' ',I8)
      
      write(ifl,1010)asciiv2_FFGF
      write(ifl,1010)newSet_FFGF
      tmpStr80_1='FrontFlow Red Grid data'
      write(ifl,1010)tmpStr80_1
      write(ifl,1010)customData_FFGF
      write(ifl,1010)ffrGrid_FFGF
      tmpStr80_1='FFR Grid header field'
      write(ifl,1010)tmpStr80_1
      write(ifl,1020)unkownDatasize_FFGF
      write(ifl,1010)gridHead_FFGF
      write(ifl,1020)fieldVersion
      write(ifl,1020)nvrtx
      write(ifl,1020)ncell
      write(ifl,1020)numofMaterial
      write(ifl,1020)NBOUND
      do ic=1,NBOUND
        write(ifl,1030)ic,NFBOUN(ic)
      end do
      write(ifl,1010)gridHeade_FFGF
      
      statRtrn=0
      end subroutine ffr_gr_headerout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_elemTypeout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
 1010 format(a80)
 1020 format(1I8)
 1040 format(I8,' ',a80)
      
      write(ifl,1010)customData_FFGF
      write(ifl,1010)ffrGrid_FFGF
      tmpStr80_1='FFR Grid element type field'
      write(ifl,1010)tmpStr80_1
      write(ifl,1020)unkownDatasize_FFGF
      write(ifl,1010)elemType_FFGF
      write(ifl,1020)fieldVersion
      write(ifl,1020)numofMaterial
      do ic=1,numofMaterial
!     In this version, name of material is not included from STAR-CD file,
!     default material name is used for every material.
        write(tmpStr80_1,*)material_ID(ic)
        write(tmpStr80_2,*)'MAT',trim(adjustl(tmpStr80_1))
        write(ifl,1040)material_ID(ic),tmpStr80_2
      end do
      write(ifl,1010)elemTypee_FFGF

      statRtrn=0
      end subroutine ffr_gr_elemTypeout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_bndryTypeout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
 1010 format(a80)
 1020 format(1I8)
 1050 format(I8,' ',I8,' ',a80,' ',a80)
      
      write(ifl,1010)customData_FFGF
      write(ifl,1010)ffrGrid_FFGF
      tmpStr80_1='FFR Grid boundary type field'
      write(ifl,1010)tmpStr80_1
      write(ifl,1020)unkownDatasize_FFGF
      write(ifl,1010)bndryType_FFGF
      write(ifl,1020)fieldVersion
      write(ifl,1020)NBOUND
      do ic=1,NBOUND
        tmpStr80_1=trim(adjustl(SFBOUN(ic)))
        tmpStr80_2=' '
        write(ifl,1050)ic,undefBndryType,tmpStr80_1,tmpStr80_2
      end do
      write(ifl,1010)bndryTypee_FFGF

      statRtrn=0
      end subroutine ffr_gr_bndryTypeout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_vrtxCordout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
 1010 format(a80)
 1020 format(1I8)
 1060 format(I8,' ',3E15.6)
      
      write(ifl,1010)customData_FFGF
      write(ifl,1010)ffrGrid_FFGF
      tmpStr80_1='FFR Grid vertices coordinate field'
      write(ifl,1010)tmpStr80_1
      write(ifl,1020)unkownDatasize_FFGF
      write(ifl,1010)vrtxCord_FFGF
      write(ifl,1020)fieldVersion
      write(ifl,1020)nvrtx
      do ic=1,nvrtx
        write(ifl,1060)ic,cord(1,ic)/gdScale,cord(2,ic)/gdScale,
     &                                    cord(3,ic)/gdScale
!       array cord is scaled by gdScale in parent program
!       divide by gdScale to rescale in FFR-GRID file.
      end do
      write(ifl,1010)vrtxCorde_FFGF

      statRtrn=0
      end subroutine ffr_gr_vrtxCordout
      
      
!---------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_bndryFaceout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
 1010 format(a80)
 1020 format(1I8)
 1030 format(I8,' ',I8)
 1070 format(5(I8,1X))
      
      do ic=1,NBOUND
        write(ifl,1010)customData_FFGF
        write(ifl,1010)ffrGrid_FFGF
        tmpStr80_1='FFR Grid Boundary Face field'
        write(ifl,1010)tmpStr80_1
        write(ifl,1020)unkownDatasize_FFGF
        write(ifl,1010)bndryFace_FFGF
        write(ifl,1020)fieldVersion
        write(ifl,1020)undefBndryType
        write(ifl,1020)ic
        write(ifl,1030)NFBOUN(ic),4
        do jc=1,NFBOUN(ic)
          write(ifl,1070)jc,IFFACE(1,jc,ic),IFFACE(2,jc,ic),
     &                        IFFACE(3,jc,ic),IFFACE(4,jc,ic)
        end do
        write(ifl,1010)bndryFacee_FFGF
      end do

      statRtrn=0
      end subroutine ffr_gr_bndryFaceout
      
      
!------------------------------------------------------------
      subroutine ffr_gr_elemVrtxout(statRtrn)
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
 1010 format(a80)
 1020 format(1I8)
 1080 format(11(I8,1X))

      write(ifl,1010)customData_FFGF
      write(ifl,1010)ffrGrid_FFGF
      tmpStr80_1='FFR Grid Element field'
      write(ifl,1010)tmpStr80_1
      write(ifl,1020)unkownDatasize_FFGF
      write(ifl,1010)elemVrtx_FFGF
      write(ifl,1020)fieldVersion
      write(ifl,1020)ncell
      do ic=1,ncell
        write(ifl,1080)ic,lvcell(1:8,ic),lacell(ic),defaultCellTypeID
      end do
      write(ifl,1010)elemVrtxe_FFGF

      statRtrn=0
      end subroutine ffr_gr_elemVrtxout

      end subroutine list_output_grid_A
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_output_grid_U
     & (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_FFRGFwords
      use module_io,only       : getfil,gdScale
      use module_boundary,only : SFBOUN,NBOUND,boundIDMap,numbname,
     &                           IFFACE,NFBOUN
!------------------------------------------------------------
!    SFBOUN     : Boundary region name read from grid file.
!    NBOUND     : Number of boundary regions read from grid file
!    boundIDMap : Indexmap from FFR-Internal boundary region ID 
!                 to grid file boundary ID
!    numbname   : Maximum number of bounary region for one boundary 
!                 condition (usually 2)
!    IFFACE     : Vertices ID constructing boundary face.
!    NFBOUN     : Number of faces in a boundary region read form grid file.
!------------------------------------------------------------
!
      implicit none

! --- [dummy arguments]
      integer,intent(in) :: mvrtx,mcell
      integer,intent(in) :: nvrtx,ncell
      real*8 ,intent(in) :: cord(3,mvrtx)
      integer,intent(in) :: lacell(  mcell)
      integer,intent(in) :: lvcell(8,mcell)
! --- [local entities]
      character(len=80),save :: fnam
      integer :: ifl,ios=0
      integer :: numofMaterial
      integer,allocatable :: material_ID(:)
      character*80 :: tmpStr80_1,tmpStr80_2
      integer :: ic,jc

      call getfil(ifl,fnam,'ffrgrid')
      open(ifl,file=fnam,form='unformatted',
     &  status='unknown',iostat=ios)
      if(ios/=0) then
        write(*,*)'*** Cannot create FrontFlowRed Grid File:',
     &                                                trim(fnam)
        return
      end if

!     preparation of grid data...
      call ffr_gr_prepData(ios)
      if(ios/=0) then
        write(*,*)'*** Error in preparation of grid data.'
        return
      end if

!     write FFR-GRID header field...
      call ffr_gr_headerout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid header.'
        return
      end if

!     write FFR-GRID element type field...
      call ffr_gr_elemtypeout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid element type.'
        return
      end if

!     write FFR-GRID boundary type field...
      call ffr_gr_bndryTypeout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid boundary type.'
        return
      end if

!     write FFR-GRID vertices coordinates field...
      call ffr_gr_vrtxCordout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid vertices coordinates.'
        return
      end if

!     write FFR-GRID boundary face list field...
      call ffr_gr_bndryFaceout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid boundary faces.'
        return
      end if

!     write FFR-GRID element vertices list field...
      call ffr_gr_elemVrtxout(ios)
      if(ios/=0) then
        write(*,*)'*** Cannot write FFR-Grid elements.'
        return
      end if
      

!     termination of FFR-GRID file output
      write(ifl)fileend_FFGF
      close(ifl)
!
      deallocate(IFFACE,NFBOUN)   
!
      return
!///////////////////////////////////////////////////////////////      
      contains
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_prepData(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
!      count up material number
      numofMaterial=maxval(lacell(1:ncell))-minval(lacell(1:ncell))+1
      do ic=minval(lacell(1:ncell)),maxval(lacell(1:ncell))
        if(count(lacell(1:ncell)==ic)==0) numofMaterial=numofMaterial-1
      end do
!      make material ID list
      allocate(material_ID(1:numofMaterial))
      jc=1
      do ic=minval(lacell(1:ncell)),maxval(lacell(1:ncell))
        if(count(lacell(1:ncell)==ic)/=0) then
          material_ID(jc)=ic; jc=jc+1
        end if
      end do
      
      statRtrn=0
      end subroutine ffr_gr_prepData
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_headerout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
      
      write(ifl)unformv2_FFGF
      write(ifl)newSet_FFGF
      tmpStr80_1='FrontFlow Red Grid data'
      write(ifl)tmpStr80_1
      write(ifl)customData_FFGF
      write(ifl)ffrGrid_FFGF
      tmpStr80_1='FFR Grid header field'
      write(ifl)tmpStr80_1
      write(ifl)unkownDatasize_FFGF
      write(ifl)gridHead_FFGF
      write(ifl)fieldVersion
      write(ifl)nvrtx
      write(ifl)ncell
      write(ifl)numofMaterial
      write(ifl)NBOUND
      do ic=1,NBOUND
        write(ifl)ic,NFBOUN(ic)
      end do
      write(ifl)gridHeade_FFGF
      
      statRtrn=0
      end subroutine ffr_gr_headerout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_elemTypeout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
      
      write(ifl)customData_FFGF
      write(ifl)ffrGrid_FFGF
      tmpStr80_1='FFR Grid element type field'
      write(ifl)tmpStr80_1
      write(ifl)unkownDatasize_FFGF
      write(ifl)elemType_FFGF
      write(ifl)fieldVersion
      write(ifl)numofMaterial
      do ic=1,numofMaterial
!     In this version, name of material is not included from STAR-CD file,
!     default material name is used for every material.
        write(tmpStr80_1,*)material_ID(ic)
        write(tmpStr80_2,*)'MAT',trim(adjustl(tmpStr80_1))
        write(ifl)material_ID(ic),tmpStr80_2
      end do
      write(ifl)elemTypee_FFGF

      statRtrn=0
      end subroutine ffr_gr_elemTypeout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_bndryTypeout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
      
      write(ifl)customData_FFGF
      write(ifl)ffrGrid_FFGF
      tmpStr80_1='FFR Grid boundary type field'
      write(ifl)tmpStr80_1
      write(ifl)unkownDatasize_FFGF
      write(ifl)bndryType_FFGF
      write(ifl)fieldVersion
      write(ifl)NBOUND
      do ic=1,NBOUND
        tmpStr80_1=trim(adjustl(SFBOUN(ic)))
        tmpStr80_2=' '
        write(ifl)ic,undefBndryType,tmpStr80_1,tmpStr80_2
      end do
      write(ifl)bndryTypee_FFGF

      statRtrn=0
      end subroutine ffr_gr_bndryTypeout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_vrtxCordout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
      
      write(ifl)customData_FFGF
      write(ifl)ffrGrid_FFGF
      tmpStr80_1='FFR Grid vertices coordinate field'
      write(ifl)tmpStr80_1
      write(ifl)unkownDatasize_FFGF
      write(ifl)vrtxCord_FFGF
      write(ifl)fieldVersion
      write(ifl)nvrtx
      do ic=1,nvrtx
        write(ifl)ic,cord(1,ic)/gdScale,cord(2,ic)/gdScale,
     &                                    cord(3,ic)/gdScale
!       array cord is scaled by gdScale in parent program
!       divide by gdScale to rescale in FFR-GRID file.
      end do
      write(ifl)vrtxCorde_FFGF

      statRtrn=0
      end subroutine ffr_gr_vrtxCordout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_bndryFaceout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1
      
      do ic=1,NBOUND
        write(ifl)customData_FFGF
        write(ifl)ffrGrid_FFGF
        tmpStr80_1='FFR Grid Boundary Face field'
        write(ifl)tmpStr80_1
        write(ifl)unkownDatasize_FFGF
        write(ifl)bndryFace_FFGF
        write(ifl)fieldVersion
        write(ifl)undefBndryType
        write(ifl)ic
        write(ifl)NFBOUN(ic),4
        do jc=1,NFBOUN(ic)
          write(ifl)jc,IFFACE(1,jc,ic),IFFACE(2,jc,ic),
     &                        IFFACE(3,jc,ic),IFFACE(4,jc,ic)
        end do
        write(ifl)bndryFacee_FFGF
      end do

      statRtrn=0
      end subroutine ffr_gr_bndryFaceout
      
      
!------------------------------------------------------------
!===============================================================      
      subroutine ffr_gr_elemVrtxout(statRtrn)
!===============================================================      
      integer,intent(out) :: statRtrn
      
      integer,parameter :: fieldVersion=1

      write(ifl)customData_FFGF
      write(ifl)ffrGrid_FFGF
      tmpStr80_1='FFR Grid Element field'
      write(ifl)tmpStr80_1
      write(ifl)unkownDatasize_FFGF
      write(ifl)elemVrtx_FFGF
      write(ifl)fieldVersion
      write(ifl)ncell
      do ic=1,ncell
        write(ifl)ic,lvcell(1:8,ic),lacell(ic),defaultCellTypeID
      end do
      write(ifl)elemVrtxe_FFGF

      statRtrn=0
      end subroutine ffr_gr_elemVrtxout

      end subroutine list_output_grid_U
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_gridfact(mvrtx,mcell,ncell,nvrtx,
     &   cord,lacell,lvcell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_boundary,only : NBOUND,IFFACE,NFBOUN
      use module_io,only       : ifle,ifll,gdScale
      use module_param ,only   : scalfctr,smallNrm
      use module_partitioner,only : ityp =>WK2
      use module_model,only       : Buffle_shft
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mvrtx,mcell,nvrtx,ncell
      real*8 ,intent(inout)  :: cord(3,mvrtx)
      integer,intent(inout)  :: lacell(mcell)
      integer,intent(in)     :: lvcell(8,mcell)
!
! --- [local entities]
!
      real*8 :: codmx1,codmx2,codmx3,codmn1,codmn2,codmn3
      real*8 :: facmn,dummn=0.d0
      real*8 :: dum1,dum2,dum3,bai
      integer:: nb,kd,mb,IS,Int,INOD,IV,ILV,IC,ierr1=0
      integer:: Nxt,NxtList(3)=(/2,3,1/),Nxt1
!

      allocate(scalfctr(0:NBOUND))
      scalfctr(:)=1.d0
!
      codmx1=cord(1,IFFACE(1,1,1))
      codmx2=cord(2,IFFACE(1,1,1))
      codmx3=cord(3,IFFACE(1,1,1))
      codmn1=cord(1,IFFACE(1,1,1))
      codmn2=cord(2,IFFACE(1,1,1))
      codmn3=cord(3,IFFACE(1,1,1))
!
!
      do 100 mb=1,NBOUND
!
      facmn=huge(dummn)
      do IS=1,NFBOUN(mb)-1
      do Int=1,3
        INOD=IFFACE(Int,IS,mb)
        IF(INOD==0) cycle
        dum1=cord(1,INOD)
        dum2=cord(2,INOD)
        dum3=cord(3,INOD)
        if(dum1>codmx1) then
          codmx1=dum1
        endif
        if(dum2>codmx2) then
          codmx2=dum2
        endif
        if(dum3>codmx3) then
          codmx3=dum3
        endif
        if(dum1<codmn1) then
          codmn1=dum1
        endif
        if(dum2<codmn2) then
          codmn2=dum2
        endif
        if(dum3<codmn3) then
          codmn3=dum3
        endif
!
        Nxt1=0
        Nxt=IFFACE(NxtList(Int),IS,mb) !   NxtList(Int)   !(/2,3,1/)
        if(Nxt==0) then
          Nxt=IFFACE(1,IS+1,mb)
          Nxt1=1
        endif
        dum1=cord(1,Nxt)-cord(1,INOD)
        dum2=cord(2,Nxt)-cord(2,INOD)
        dum3=cord(3,Nxt)-cord(3,INOD)
        dummn=dsqrt(dum1**2+dum2**2+dum3**2)
        if(dummn<facmn) then
           facmn=dummn
        endif
        IF(Nxt1==1) exit
      enddo
      enddo
      if(facmn<1.0d-16) facmn=1.0d-16
      scalfctr(mb)=facmn
 100  enddo
      scalfctr(0)=dsqrt((codmx1-codmn1)**2
     &                 +(codmx2-codmn2)**2
     &                 +(codmx3-codmn3)**2)
      
      do mb=0,NBOUND
      write(ifll,'(a,I8,E12.4)') 
     &  'MSG: BC. NO= and Scale Factor = ',mb,scalfctr(mb)
      enddo
!
!
!
!      if(.true.) then
      if(abs(Buffle_shft)>1.d-15) then
        bai=Buffle_shft
        dum1=gdScale*bai
        write(*,'(1X,a,E12.4)' ) 
     & 'MSG: in ffr_gr_elemVrtxout  0.2', dum1
        allocate(ityp(0:nvrtx),stat=ierr1) 
        
        if(ierr1.ne.0) 
     &  stop 'STOP at allocating ityp(:) in list_gridfact'
        ityp(:)=0
        do IC=1,ncell
        if(lacell(IC)>80) then
          do ILV=1,8
            IV=lvcell(ILV,IC)
            ityp(IV)=1
          enddo
        endif

        enddo
        ityp(0)=0
        do IV=1,nvrtx
        if(ityp(IV)==1) then
!          cord(1,IV)=
!          cord(2,IV)=
          cord(3,IV)=cord(3,IV)+dum1
        endif
        enddo
        deallocate(ityp)
      endif
!
      return
      end subroutine list_gridfact
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_BC_connectivity(nface,mface,lvface,lbface)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_boundary,only    : nbcnd,kdbcnd,SFBOUN,NBOUND,
     &                              boundName,IFFACE,NFBOUN,numbname,
     &                              boundIDMap
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mface,nface
      integer,intent(in)    :: lvface(4,mface)
      integer,intent(in)    :: lbface(2,mface)
!
! --- [local entities]
!
      integer :: kd,nb,mb,mb2,idoub,j,IS,IBS
!
        NFBOUN(:)=0
        IFFACE(:,:,:)=0
        do IS=1,nface
        nb=lbface(1,IS)
        if(nb==0) cycle
        do j=1,2
        mb=boundIDMap(abs(nb),j)
        if(mb/=-1) then
          NFBOUN(mb)=NFBOUN(mb)+1
          IBS=NFBOUN(mb)
          IFFACE(:,IBS,mb)=lvface(:,IS)
        endif
        enddo
        enddo
!
        return
!
      end subroutine list_BC_connectivity
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$







