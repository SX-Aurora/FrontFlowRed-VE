      subroutine setupGridGrobalInfo(gdFormat,rsFormat,rtrnStat)
      use FFRdata, only : nvrtx,ncell,nface,cord,NBFS,NFCE,
     &                    NBOUND,NFBOUN,SFBOUN,IBFACE,IFFACE,
     &                    numVars,numBVars,nameVars,
     &                    lvcell,lacell,kmesh,
     &                    lfcell,lvface,lcface,lbface,LVRTCV,LCVFAC,
     &                    nbcnd,ivbcnd,lvbcnd
      use FVdata, only : sfRsltFlg,sfClkFlg
      use FLUENTdata, only : F_NREGN,F_NFBOUN

      implicit none
!
      character*2,intent(in)  :: gdFormat,rsFormat
      integer,    intent(out) :: rtrnStat
!
      integer :: i,j,k,l
      integer :: iMAT,prevMAT
      integer :: ndCl(1:8)
!
      integer :: IC1,IC2,IS,nb
      integer :: mcell,mface,mvrtx
!      integer :: ncell,nvrtx,nface
      integer :: ncelb,NCV
      integer :: match,IVN0(4),IVN1(4),ierror

      
      allocate(sfRsltFlg(1:NBOUND))
      sfRsltFlg=0              ! Boundary variables flag
      allocate(sfClkFlg(1:NBOUND))
      sfClkFlg=0              ! Boundary surface clockness flag

!     make up kmesh, if not exist.
      if( .false. ) then
!      if( .not.(allocated(kmesh))) then
         write(*,*) ' Warning: kmesh is not created in read_griddata.'
         allocate(kmesh(1:ncell))
         do i=1,ncell
            k=0
            do j=1,8
               if(lvcell(j,i)>0) k=k+1
            end do
         end do
         select case(k)
         case(4)
            kmesh(i) = 1        ! tetra
         case(8)
            kmesh(i) = 2        ! hexa
         case(6)
            kmesh(i) = 3        ! prism
         case(5)
            kmesh(i) = 4        ! pyramid
         end select
      end if

!
!     make up lvcell order for FrontFlow format.
!     *** following should be use for tempolarily
!          to solve some kind of mesh order problem.
      if(.false.) then
!      if(gdFormat=='FF') then
      write(*,*) ' (setupGridGlobalInfo):convert lvcell'
      do i=1,ncell
         ndCl(1:8)=lvcell(1:8,i)
         select case(kmesh(i))
         case(1)   ! tetra
            lvcell(1,i)=ndCl(1)
            lvcell(2,i)=ndCl(2)
            lvcell(3,i)=ndCl(3)
            lvcell(4,i)=ndCl(4)
         case(2)   ! hexa
            lvcell(1,i)=ndCl(1)
            lvcell(2,i)=ndCl(2)
            lvcell(3,i)=ndCl(4)
            lvcell(4,i)=ndCl(3)
            lvcell(5,i)=ndCl(5)
            lvcell(6,i)=ndCl(6)
            lvcell(7,i)=ndCl(8)
            lvcell(8,i)=ndCl(7)
         case(3)   ! prism
            lvcell(1,i)=ndCl(1)
            lvcell(2,i)=ndCl(2)
            lvcell(3,i)=ndCl(3)
            lvcell(4,i)=ndCl(4)
            lvcell(5,i)=ndCl(5)
            lvcell(6,i)=ndCl(6)
         case(4)   ! pyramid
            lvcell(1,i)=ndCl(1)
            lvcell(2,i)=ndCl(2)
            lvcell(3,i)=ndCl(3)
            lvcell(4,i)=ndCl(4)
            lvcell(5,i)=ndCl(5)
         case default
            write(*,*) 'Error: unknown element type.'
            stop
            rtrnStat=1
         end select

      end do
      end if


!     ----------------------------------------------------
!     make up lfcell,lvface,lcface,lbface,LVRTCV,LCVFAC, 
!         for FLUENT output.
      if(rsFormat=='FL') then

!        set up dummy arguments
         mcell=ncell
         mvrtx=nvrtx
         mface=ncell*6
         allocate(lfcell(7,mcell),stat=ierror)
         allocate(lvface(4,mface),stat=ierror)
         allocate(lcface(2,mface),stat=ierror)
         allocate(lbface(2,mface),stat=ierror)
         allocate(LVRTCV(  1),stat=ierror)
         allocate(LCVFAC(4,1),stat=ierror)
!         allocate(LVRTCV(  mvrtx),stat=ierror)
!         allocate(LCVFAC(4,mface),stat=ierror)
         if(ierror.ne.0) stop 'stop at allocating dummy arguments,
     &        in append_resultfluent'
      
         lfcell=0
         lvface=0
         lcface=0
         lbface=0
         LVRTCV=0
         LCVFAC=0

!        check lvcell for FFR face rule( -1 -> 0 )
         do i=1,ncell
            do j=1,8
               if(lvcell(j,i)<0) then
                  lvcell(j,i)=0
               end if
            end do
         end do

!        ----------------------------------------
!        search faces of cell
!         [grid]  -> [result]
!          GM-x   ->  FL    : no need(need?)
!         (other) ->  FL    : need
!        ----------------------------------------
!        *** use lvcell in FFR face order
         write(*,*) '    ###    LIST FACE(FL) ... '
         call list_face_FLUENT
     &        (mcell,mface,mvrtx,nvrtx,ncell,nface,ncelb,NCV,
     &        lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,ierror)
         if( ierror.ne.0 ) goto 9999
 9999    continue

!        *** make ivbcnd,lvbcnd for list_facebc
!         ivbcnd=0
!         lvbcnd=0
         nbcnd   =NBOUND
         write(*,*) '    ###    storeGridData(FL) ... '
         call storeGridData_FLUENT(ierror)
         if(ierror/=0) then
            rtrnStat=1
            return
         end if

!        *** make lbface
         write(*,*) '    ###    LIST_FACEBC(FL) ... '
         call list_facebc_FLUENT
     &        (mface,nface,mvrtx,nvrtx,ncell,mcell,
     &        lvface,lbface,lcface,lacell)

!        *** set NFCE
         if(NFCE/=nface) then
            write(*,*) '    *** checked: set nface:',nface
            NFCE=nface
            if(nface<=0) then
               write(*,*) 'Error: there is no face.'
               stop
            end if
         end if

!        *** make up F_NFBOUN insted of NFBOUN for new lvface,lvcell...
         allocate(F_NFBOUN(0:NBOUND))
         F_NFBOUN=0
         do IS=1,nface
            IC1=lcface(1,IS)
            IC2=lcface(2,IS)
            nb=lbface(1,IS)
!           check
            if(nb>NBOUND .or. nb<0) then
               write(*,*) 'Error: IS=',IS,' has wrong boundary id.'
               stop
            end if
            if(nb==0 .and. IC2>ncell) then
!              *** unknown boundary
            else
               F_NFBOUN(nb)=F_NFBOUN(nb)+1
            end if
!            if(IC2>ncell) then
!               F_NFBOUN(nb)=F_NFBOUN(nb)+1
!            end if
         end do

!        *** check F_NREGN in module FLUENTdata
         if(F_NREGN<NBOUND) then
            F_NREGN=NBOUND
            write(*,*) '    *** checked: set F_NREGN:',F_NREGN
         end if
         if(F_NREGN<=0) then
            write(*,*) 'Error: there is no boundary.'
            stop
         end if
! ---    obtain BC information from (fflow.ctl)
         call read_fort1_boundary

      end if


      rtrnStat=0
      end subroutine setupGridGrobalInfo



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_face_FLUENT
     & (mcell,mface,mvrtx,nvrtx,ncell,nface,ncelb,NCV,
     & lvcell,lfcell,lvface,lcface,LVRTCV,LCVFAC,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! '!F' lines are comment outed from 'list_face'.
!      of src_022_vof_portable/src_pre/src/pre_list_total.f
!
!
! --- [module arguments]
!
!F      use module_io,only : ifle,ifll
!F      use module_partitioner,only : nclv =>WK1
!F      use module_partitioner,only : ityp =>WK2
!F      use module_partitioner,only : ip   =>WK3
!F      use module_partitioner,only : lvrtx=>IW1
!F      use module_parameter,  only : NIFACE
!F      use module_boundary,only    : nbcnd,ivbcnd,lvbcnd,NBOUND
!
! --- 1. Make up list of connectivity for cell & face
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mcell,mface,nvrtx,ncell,mvrtx
      integer,intent(out) :: nface,ncelb
      integer,intent(in)  :: lvcell(8,mcell)
      integer,intent(out) :: lfcell(7,mcell)
      integer,intent(out) :: lvface(4,mface)
      integer,intent(out) :: lcface(2,mface)
      integer,intent(out) :: LVRTCV(  1)
      integer,intent(out) :: LCVFAC(4,1)
!F      integer,intent(out) :: LVRTCV(  mvrtx)
!F      integer,intent(out) :: LCVFAC(4,mface)
      integer,intent(out) :: ierror,NCV
!F
!F    dummy arguments insted of using modules.
!F
      integer :: ifle=6,ifll=6
      integer,allocatable :: nclv(:)
      integer,allocatable :: ityp(:)
      integer,allocatable :: ip(:)
      integer,allocatable :: lvrtx(:,:)
      integer :: NIFACE
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
      integer :: lvf (4,6,4),nb
      integer :: lvrt(4)
!F      integer :: lvrt(4),LBC(nbcnd)
      integer :: i,j,k,l,m,n,IS,IVV
      integer :: j1,j2,n1,n2,m1,m2,lf1,ie,ierr1=0,ierr2
      integer :: kclv,nf,icmax,match,nch,nfacx,nclbx
      integer :: IC,IV,ILV,ILS,NF1,NF2,ILS1,ILS2,IC1,IC2,ICTP
      logical :: INNER
!
! --- 
!
      allocate(nclv(8*ncell),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating nclv(:) in list_face'
      allocate(ityp(ncell),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ityp(:) in list_face'
      allocate(ip(0:nvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating ip(:) in list_face'
!
!F      NCV=NVRTX
!F      DO 50 IV=1,NCV
!F      LVRTCV(IV)=IV
!F  50  CONTINUE
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
      if(lvf(4,j,k).lt.1) lvf(4,j,k)=8
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
!F      LCVFAC(i,IS)=0
  111 continue
!
! --- < 1.3 set type & no. of faces >--
!
      do 120 IC=1,ncell
      kclv=0
      do 121 ILV=1,8
      if(lvcell(ILV,IC).gt.0) kclv=ILV
  121 continue
      if(kclv.lt.1) then
        write(ifle,*) ' #### ERR: kclv= ',kclv
        goto 9001
      endif
      ityp(IC)=nv2typ(kclv)
      if(ityp(IC).lt.1) then
         write(ifle,*) ' #### ERR: ityp(IC)= ',ityp(IC),kclv
         goto 9001
      endif
      lfcell(7,IC)=nfacel(ityp(IC))
  120 continue
!
! --- < 2. Set up cell-face list >-
!
! --- < 2.1 make up list for cells connected with each vertex >--
!
      call list_onvrtx(8,mcell,nvrtx,ncell,
     & lvcell,ip,nclv,icmax)
!
      allocate(lvrtx(6,4*icmax),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating lvrtx in list_face'
!
! --- < 2.2 Search internal faces >-
! --- icmax is maximun no. of cells which connect one vertex.
!
      nface=0
!
      do 200 IV=1,nvrtx
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
! --- / matching procedure /
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
        if(match.gt.0) goto 9003
      else
        nch=IC2
      endif
  231 continue
  230 continue
!
  200 continue
!
      write(ifll,*) '    ###    INTER FACE : NIFACE= ',NFACE
      NIFACE=Nface
!
!-< 3. Search boundary faces >-
!
      ncelb=ncell
      do 300 IC=1,ncell
      k=ityp(IC)
      do 301 IS=1,lfcell(7,IC)
      if(lfcell(IS,IC).lt.1) then
        ncelb=ncelb+1
        nface=nface+1
        nclbx=min(ncelb,mcell)
        nfacx=min(nface,mface)
        lfcell(IS,IC)=nface
        lcface(1,nfacx)=IC
        lcface(2,nfacx)=ncelb
        do 302 iv=1,4
        lvface(iv,nfacx)=lvcell(lvf(iv,IS,k),IC)
  302   continue
      endif
  301 continue
  300 continue
!
      deallocate(lvrtx,nclv,ip,ityp)
      write(ifll,*) '    ###    TOTAL FACE: NFACE= ',NFACE
!
!F      do 410 IS=1,nface
!F      do 410 iv=1,4
!F      LCVFAC(iv,IS)=lvface(iv,IS)
!F  410 continue
!F!
!F      call list_facerr(ifle,mcell,mface,nface,ncelb,ierr1)
!F      if( ierr1.ne.0 ) goto 9999
!F!
!F      allocate(ip(mface),nclv(0:mvrtx),stat=ierr1)
!F      if(ierr1.ne.0) stop 'stop at allocating ip(:)_2 in list_face'
!F!
!F      if(ivbcnd(nbcnd).eq.0) then
!F        NBOUND=NBOUND+1
!F        ip(:)=0
!F        do nb=1,nbcnd-1
!F          nclv(:)=0
!F          nclv(0)=1
!F          do i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F             iv=lvbcnd(i)
!F             nclv(iv)=1
!F          enddo
!F          do 500 IS=NIFACE,nface
!F            do ivv=1,4
!F              iv=lvface(ivv,IS)
!F              if(nclv(iv).eq.0) goto 500
!F            enddo
!F            ip(IS)=nb
!F 500      continue
!F        enddo
!F        allocate(ityp(ivbcnd(nbcnd-1)),stat=ierr1)
!F        do nb=1,nbcnd-1
!F          do i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F          ityp(i)=lvbcnd(i)
!F          enddo
!F        enddo
!F        deallocate(lvbcnd)
!F!
!F        nclv(:)=0
!F        do IS=NIFACE,nface
!F          if(ip(IS).eq.0) then
!F            do ivv=1,4
!F            iv=lvface(ivv,IS)
!F            nclv(iv)=1
!F            enddo
!F          endif
!F        enddo
!F!
!F        ivv=ivbcnd(nbcnd-1)
!F        do iv=1,nvrtx
!F        if(nclv(iv).eq.1) then
!F          ivv=ivv+1
!F        endif
!F        enddo
!F        ivbcnd(nbcnd)=ivv
!F!
!F        allocate(lvbcnd(ivbcnd(nbcnd)),stat=ierr1)
!F        DO nb=1,nbcnd-1
!F          do i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F            lvbcnd(i)=ityp(i)
!F          enddo
!F        ENDDO
!F        ivv=ivbcnd(nbcnd-1)
!F        do iv=1,nvrtx
!F          if(nclv(iv).eq.1) then
!F            ivv=ivv+1
!F            lvbcnd(ivv)=iv
!F          endif
!F        enddo
!F        deallocate(ityp)
!F      endif
!F!
!F      deallocate(ip,nclv)
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
      write(ifle,*) ' ### error : allocation failed'
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
      end subroutine list_face_FLUENT
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_onvrtx(mv,mcell,nvrtx,ncell,
     & lvcell,ip,nclv,icmax)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! '!F' lines are comment outed from 'list_face'.
!      of src_022_vof_portable/src_pre/src/pre_list_total.f
!
!F      use module_partitioner,only : kclv=>WK5
!F      use module_partitioner,only : iq  =>WK6
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mv,mcell,nvrtx,ncell
      integer,intent(in)  :: lvcell(mv,mcell)
      integer,intent(out) :: ip    (0:nvrtx)
      integer,intent(out) :: nclv  (mv*ncell)
      integer,intent(inout) :: icmax
!F
!F    dummy arguments insted of using modules.
!F
      integer,allocatable :: kclv(:)
      integer,allocatable :: iq(:)
!
! --- [local entities]
!
!      integer :: kclv(ncell)   
!      integer :: iq(nvrtx)
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
!===============================================================
      subroutine storeGridData_FLUENT(ierrorFVG)
!===============================================================
! '!F' lines are comment outed from 'storeGridData'.
!      of src_022_vof_portable/src_pre/src/pre_read_griddata.f

!F    for FLUENT only
      use FFRreaddata,only : cntlnam
      use FFRdata,    only : nvrtx,NBOUND,SFBOUN,NFBOUN,
     &                       NBFS,IBFACE,IFFACE,
     &                       ivbcnd,lvbcnd,nbcnd,numbname,boundIDMap
      implicit none
!
! --- [module arguments]
!
      integer,intent(inout)     :: ierrorFVG
!
! --- [local entities]
!
      logical,allocatable :: tempFlag(:)
      integer,allocatable :: givbcnd(:),glvbcnd(:),pMap(:)
      real(8),parameter   :: smallNrm=1.0d-13
      real(8) :: max1cord,min1cord,max2cord,min2cord,bounChSize
      real(8)             :: dum1
      integer             :: IBIDa,IBIDb
      integer             :: tempInt,K,J,I,M,IV,Int,nb,mb,IS
      integer             :: nbcndde
!
!F    for FLUENT only
      character*80,allocatable :: boundName(:,:)  ! boundary names in fflow.ctl
      integer,parameter :: ifll=6,ifle=6
!
! --- 
!
!F    for FLUENT only
      allocate(boundName(nbcnd,numbname))
      boundName=' '
      do nb=1,nbcnd
         boundName(nb,1)=trim(adjustl(SFBOUN(nb)))
      end do
!F
!
      allocate(boundIDMap(nbcnd,numbname))
      boundIDMap=-1
!
      nbcndde=nbcnd
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j).eq.'undefined') then
        nbcndde=nbcnd-1
        boundIDMap(nbcnd,1)=0
        boundIDMap(nbcnd,2)=-1
!F        SFBOUN(0)='undefined'
      endif
      enddo
      enddo
!
      IF(nbcndde.LT.nbcnd-1) THEN
        WRITE(ifle,*)
     &   ' ### ERR: NAME [undefined] BC too many in ',cntlnam
        ierrorFVG=1
        return
      ENDIF
!
      if(NBOUND.lt.nbcndde) then
	WRITE(ifle,*) '(append_refultfluent)',
     &   'Number of boundary region is not matched',NBOUND,nbcnd
	ierrorFVG=1
        return
      endif
!
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j).ne.' ') then
        do mb=1,NBOUND
        if(trim(adjustl(boundName(nb,j))).eq.trim(adjustl(SFBOUN(mb))))
     & then
!        if(adjustl(boundName(nb,j)).eq.adjustl(SFBOUN(mb))) then
          boundIDMap(nb,j)=mb
          exit
        endif
	enddo
      end if
      end do
      end do
!
!F      allocate(tempFlag(0:NBOUND))
!F      tempFlag=.false.
!F!
!F      do j=1,numbname
!F      do nb=1,nbcnd
!F      mb=boundIDMap(nb,j)
!F      if(boundIDMap(nb,j).ne.-1) then
!F        if((boundIDMap(nb,j)<0).or.(boundIDMap(nb,j)>NBOUND)) then
!F	  WRITE(ifle,*) '(append_refultfluent)',
!F     &   'Wrong boundary ID.',boundIDMap(nb,j)
!F          ierrorFVG=1
!F	  return
!F        end if
!F        if(kdbcnd(0,nb).ne.kdtchi) then
!F          if(tempFlag(mb)) then
!F	    WRITE(ifle,*)  
!F     &      '(append_refultfluent) Boundary [',
!F     &      trim(SFBOUN(boundIDMap(nb,j))),'] is defined twice.'
!F            ierrorFVG=1
!F            return
!F          end if
!F        endif
!F        tempFlag(mb)=.true.
!F      end if
!F      end do
!F      end do
!F!
!F      do mb=1,NBOUND
!F      if(.not.(tempFlag(mb))) then
!F        WRITE(ifle,*)  
!F     & '(append_refultfluent)Boundary name: [',
!F     &   trim(SFBOUN(mb)),'] is NOT defined in ',cntlnam,' file.'
!F! onishi
!F        WRITE(ifle,*) '    *** please check "boundary" in ',cntlnam
!F        WRITE(ifle,*) '       1.miss type or missing in ',cntlnam
!F        WRITE(ifle,*) '       2.wrong charactor codes(CR<->CR+LF etc.',
!F     &       ' wrong FTP mode(ASCII/Binary) causes this error)'
!F        WRITE(ifle,*) '       3.wrong file, etc.'
!F	ierrorFVG=1
!F	return
!F      end if
!F      end do
!F      deallocate(tempFlag)
!
      allocate(tempFlag(0:nvrtx))
      allocate(givbcnd(0:NBOUND))
      givbcnd=0
!
      do 500 mb=1,NBOUND
      tempFlag(:)=.false.
      do 540 IS=1,NBFS
      if(IBFACE(IS)==mb) then
      do 510 Int=1,4
      IV=IFFACE(Int,IS)
      if(IV.gt.0) tempFlag(IV)=.true.
 510  continue
      endif
 540  continue
      do iv=1,nvrtx
      if(tempFlag(iv)) then
        givbcnd(mb)=givbcnd(mb)+1
      endif
      end do
 500  continue
!
      do 550 mb=1,NBOUND
!F      write(ifll,*) 'MSG: BC NODE ',mb,trim(sfboun(mb)),givbcnd(mb)
      givbcnd(mb)=givbcnd(mb)+givbcnd(mb-1)
 550  continue
!
      allocate(glvbcnd(givbcnd(NBOUND)))
      glvbcnd(:)=0
      tempInt=0
!
      do 600 mb=1,NBOUND
      tempFlag=.false.
      do 640 IS=1,NBFS
      if(IBFACE(IS)==mb) then
      do 610 Int=1,4
      IV=IFFACE(Int,IS)
      if(IV.gt.0) tempFlag(IV)=.true.
 610  continue
      end if
 640  continue
      do iv=1,nvrtx
      if(tempFlag(iv)) then
        tempInt=tempInt+1
        glvbcnd(tempInt)=iv
      endif
      enddo
      if(tempInt.ne.givbcnd(mb)) then
        WRITE(ifle,*)
     &  '(append_refultfluent)',
     &  'Number of boundary vertex is wrong.',mb,tempInt,givbcnd(mb)
	ierrorFVG=1
	return
      end if
 600  continue
!
!F      allocate(icbinr(0:nbcnd))
!F      icbinr=0
!F      allocate(lcbinr(icbinr(nbcnd)))
!F!
!F! --- define periodic BC and touch-inlet BC
!F!
!F      do 200 nb=1,nbcnd
!F      if(kdbcnd(0,nb).eq.kdprdc) then
!F!
!F! --- periodic BC:
!F!
!F	IBIDa=boundIDMap(nb,1)
!F 	IBIDb=boundIDMap(nb,2)
!F        max1cord=0.d0
!F        min1cord=0.d0
!F        max2cord=0.d0
!F        min2cord=0.d0
!F        bounChSize=1.d0
!F        if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
!F          do 700 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!F          max1cord=max(max1cord,cord(2,glvbcnd(j)))
!F          max2cord=max(max2cord,cord(3,glvbcnd(j)))
!F          min1cord=min(min1cord,cord(2,glvbcnd(j)))
!F          min2cord=min(min2cord,cord(3,glvbcnd(j)))
!F  700     continue
!F        else if((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
!F          do 710 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!F          max1cord=max(max1cord,cord(3,glvbcnd(j)))
!F          max2cord=max(max2cord,cord(1,glvbcnd(j)))
!F          min1cord=min(min1cord,cord(3,glvbcnd(j)))
!F          min2cord=min(min2cord,cord(1,glvbcnd(j)))
!F  710     continue
!F        else if((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
!F          do 720 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!F          max1cord=max(max1cord,cord(1,glvbcnd(j)))
!F          max2cord=max(max2cord,cord(2,glvbcnd(j)))
!F          min1cord=min(min1cord,cord(1,glvbcnd(j)))
!F          min2cord=min(min2cord,cord(2,glvbcnd(j)))
!F  720     continue
!F        end if
!F        bounChSize=abs(max1cord-min1cord)**2+abs(max2cord-min2cord)**2
!F	if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
!F     &     (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
!F	  WRITE(ifle,*)  
!F     &    '(append_refultfluent)',
!F     &    'Number of periodic boundary vertex is not matched.',
!F     &    (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
!F     &    (givbcnd(IBIDb)-givbcnd(IBIDb-1))
!F	  ierrorFVG=1
!F	  return
!F	end if
!F	allocate(pMap(givbcnd(IBIDb-1)+1: 
!F     &                   givbcnd(IBIDb)))
!F	pMap(:)=-1
!F	if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
!F	  do 210 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!F	  do 220 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!F          if(((cord(2,glvbcnd(k))-cord(2,glvbcnd(j))-
!F     &         prdcOffset(nb,2))**2+
!F     &        (cord(3,glvbcnd(k))-cord(3,glvbcnd(j))-
!F     &         prdcOffset(nb,3))**2)
!F     &         .lt.smallNrm*bounChSize) then
!F	    pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
!F     &      = glvbcnd(k)
!F	    exit
!F	  end if
!F 220      continue
!F 210      continue
!F	else if((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
!F	  do 230 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!F	  do 240 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!F          if(((cord(1,glvbcnd(k))-cord(1,glvbcnd(j))-
!F     &         prdcOffset(nb,1))**2+
!F     &         (cord(3,glvbcnd(k))-cord(3,glvbcnd(j))-
!F     &         prdcOffset(nb,3))**2)<smallNrm*bounChSize) then
!F	    pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
!F     &      = glvbcnd(k)
!F  	    exit
!F          end if
!F 240      continue
!F 230      continue
!F        else if((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
!F          do 250 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!F          do 260 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!F            if(((cord(2,glvbcnd(k))-cord(2,glvbcnd(j))-
!F     &           prdcOffset(nb,2))**2+
!F     &          (cord(1,glvbcnd(k))-cord(1,glvbcnd(j))-
!F     &          prdcOffset(nb,1))**2)<smallNrm*bounChSize) then
!F	    pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
!F     &      = glvbcnd(k)
!F	    exit
!F	  end if
!F 260      continue
!F 250      continue
!F	else
!F	  WRITE(ifle,*)  
!F     &    '(append_refultfluent)Periodic axis setting is wrong.',
!F     &    prdcAxis(nb)
!F	  ierrorFVG=1
!F	  return
!F	end if
!F        do 270 j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!F	if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
!F	  WRITE(ifle,*)  
!F     &    '(append_refultfluent)',
!F     &    'Periodic boundary mapping error (no vertex).'
!F	  ierrorFVG=1
!F	  return
!F	end if
!F 270    continue
!F        do 280 j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
!F        do 290 k=j+1,givbcnd(IBIDb)
!F	if(pMap(j)==pMap(k)) then
!F	  WRITE(ifle,*)  
!F     &    '(append_refultfluent)',
!F     &    'Periodic boundary mapping error (duplicated vertex).'
!F	  ierrorFVG=1
!F	  return
!F	end if
!F 290    continue
!F 280    continue
!F!
!F        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
!F     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!F!
!F	if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
!F	  do 300 j=1,givbcnd(IBIDa)-givbcnd(IBIDa-1)
!F	  if(((cord(2,glvbcnd(j+givbcnd(IBIDb-1)))-
!F     &         cord(2,glvbcnd(j+givbcnd(IBIDa-1)))-
!F     &         prdcOffset(nb,2))**2+
!F     &        (cord(3,glvbcnd(j+givbcnd(IBIDb-1)))-
!F     &         cord(3,glvbcnd(j+givbcnd(IBIDa-1)))-
!F     &         prdcOffset(nb,3))**2)>=smallNrm) then
!F            WRITE(ifle,*)  
!F     &      '(append_refultfluent)',
!F     &      'Periodic boundary mapping error (x-axis check).'
!F            ierrorFVG=1
!F            return
!F          end if
!F 300      continue
!F	else if((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
!F	  do 310 j=1,givbcnd(IBIDa)-givbcnd(IBIDa-1)
!F	  if(((cord(1,glvbcnd(j+givbcnd(IBIDb-1)))-
!F     &         cord(1,glvbcnd(j+givbcnd(IBIDa-1)))-
!F     &         prdcOffset(nb,1))**2+
!F     &        (cord(3,glvbcnd(j+givbcnd(IBIDb-1)))-
!F     &         cord(3,glvbcnd(j+givbcnd(IBIDa-1)))- 
!F     &         prdcOffset(nb,3))**2)>=smallNrm) then
!F            WRITE(ifle,*) '(append_refultfluent)',
!F     &      'Periodic boundary mapping error (y-axis check).'
!F            ierrorFVG=1
!F	    return
!F	  end if
!F 310      continue
!F        else if((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
!F	  do 320 j=1,givbcnd(IBIDa)-givbcnd(IBIDa-1)
!F	  if(((cord(2,glvbcnd(j+givbcnd(IBIDb-1)))-
!F     &         cord(2,glvbcnd(j+givbcnd(IBIDa-1)))-
!F     &         prdcOffset(nb,2))**2+
!F     &        (cord(1,glvbcnd(j+givbcnd(IBIDb-1)))-
!F     &         cord(1,glvbcnd(j+givbcnd(IBIDa-1)))-
!F     &         prdcOffset(nb,1))**2)>=smallNrm) then
!F	    WRITE(ifle,*)  
!F     &      '(append_refultfluent)',
!F     &      'Periodic boundary mapping error (z-axis check).'
!F	    ierrorFVG=1
!F	    return
!F	  end if
!F 320      continue
!F	end if
!F	deallocate(pMap)
!F!
!F      elseif(kdbcnd(0,nb).eq.kdsld) then
!F!
!F! --- slinding mesh:
!F!
!F       	IBIDa=boundIDMap(nb,1)
!F 	IBIDb=boundIDMap(nb,2)
!F        allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
!F	pMap(:)=-1
!F!
!F        if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
!F     &     (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
!F	  WRITE(ifle,*) '(append_refultfluent)',
!F     &    'Number of slinding boundary vertex is not matched.',
!F     &    (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
!F     &    (givbcnd(IBIDb)-givbcnd(IBIDb-1))
!F	  ierrorFVG=1
!F	  return
!F	end if
!F!
!F        do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!F	do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!F!
!F        dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
!F     &       (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
!F     &       (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
!F!
!F        if(dum1.lt.1.d-7) then
!F          pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))=glvbcnd(k)
!F  	exit
!F        endif
!F        enddo
!F        enddo
!F!
!F        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!F	if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
!F	  WRITE(ifle,*)
!F     &    '(append_refultfluent)',
!F     &    'Slinding boundary mapping error (no vertex).'
!F	  ierrorFVG=1
!F	  return
!F	end if
!F        enddo
!F!
!F        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
!F        do k=j+1,givbcnd(IBIDb)
!F	if(pMap(j)==pMap(k)) then
!F	  WRITE(ifle,*)  
!F     &    '(append_refultfluent)',
!F     &    'Slinding boundary mapping error (duplicated vertex).'
!F	  ierrorFVG=1
!F	  return
!F	endif
!F        enddo
!F        enddo
!F!
!F        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
!F     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!F!
!F        deallocate(pMap)
!F      endif
!F 200  continue
!F!
      allocate(ivbcnd(0:nbcnd))
      allocate(lvbcnd(givbcnd(NBOUND)))
      ivbcnd=0
      lvbcnd=0
      do 400 nb=1,nbcndde
      IBIDa=boundIDMap(nb,1)
      IBIDb=boundIDMap(nb,2)
      ivbcnd(nb)=ivbcnd(nb-1)
      ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDa)-givbcnd(IBIDa-1)
      if((ivbcnd(nb)-ivbcnd(nb-1))/=(givbcnd(IBIDa) 
     &  -givbcnd(IBIDa-1))) then
	WRITE(ifle,*)  
     &  '(append_refultfluent)',
     &  'Number of periodic boundary vertex is not matched.'
	ierrorFVG=1
	return
      end if
      lvbcnd(ivbcnd(nb-1)+1:ivbcnd(nb))=
     &  glvbcnd(givbcnd(IBIDa-1)+1:givbcnd(IBIDa))
!F      if(kdbcnd(0,nb)==kdprdc.or.
!F     &   kdbcnd(0,nb)==kdsld.or.
!F     &   kdbcnd(0,nb)==kdintr) then
!F	tempInt=ivbcnd(nb)
!F	ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDb)-givbcnd(IBIDb-1)
!F	lvbcnd(tempInt+1:ivbcnd(nb))=
!F     &  glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!F!      elseif(kdbcnd(0,nb)==kdintr) then
!F!	tempInt=ivbcnd(nb)
!F!	ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDb)-givbcnd(IBIDb-1)
!F!	lvbcnd(tempInt+1:ivbcnd(nb))= 
!F!     &  glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!F      end if
 400  continue
!
      deallocate(givbcnd)
      deallocate(glvbcnd)
      deallocate(tempFlag)
!F      deallocate(IFFACE)
!F      deallocate(NFBOUN)
      ierrorFVG=0
!
      end subroutine storeGridData_FLUENT
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine list_facebc_FLUENT
     & (mface,nface,mvrtx,nvrtx,ncell,mcell,
     &  lvface,lbface,lcface,lacell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! '!F' lines are comment outed from 'list_facebc'.
!      of src_022_vof_portable/src_pre/src/pre_list_total.f

!F    for FLUENT only
      use FFRdata,    only : SFBOUN,nbcnd,ivbcnd,lvbcnd

!F      use module_io,only          : ifle
!F      use module_boundary,only    : nbcnd,ivbcnd,lvbcnd,kdbcnd,kdprdc,
!F     &                              kdintr,kdbuff,kdsld,boundName
!F      use module_partitioner,only : lbc=>WK1
!F      use module_partitioner,only : lbfacx=>WK2
!F      use module_parameter,  only : NIFACE
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
      integer :: i,m,n,nb,mb
      integer :: IC,IS,IV,IE,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC,ierr1=0
!F    for FLUENT only
      integer,allocatable :: lbc(:)
!
      allocate(lbc(0:nvrtx),stat=ierr1) 
      if(ierr1.ne.0) stop 'stop at allocating lbc(:) in list_facebc'
!      allocate(lbfacx(mface),stat=ierr1) 
!      if(ierr1.ne.0) stop 'stop at allocating lbfacx(:) in list_facebc'
!
! --- 
!
      lbc(1:nvrtx)=0
      lbc(0)=1
!
!F      do 110 IS=1,nface
!F      nb=lbface(1,IS)
!F      lbfacx(IS)=nb
!F  110 continue
!F!
      do 200 nb=1,nbcnd
!F      if(kdbcnd(0,nb).ne.kdprdc) then
!F! --- non-periodic BC:
!F        if(kdbcnd(0,nb).ne.kdintr.and.
!F     &     kdbcnd(0,nb).ne.kdbuff.and.
!F     &     kdbcnd(0,nb).ne.kdsld
!F     &   ) then
          do 301 i=ivbcnd(nb-1)+1,ivbcnd(nb)
          iv=lvbcnd(i)
          lbc(iv)=1
  301     continue
          do 310 IS=1,nface
!F *** important for buffle
!F          if(lcface(2,IS).le.ncell) goto 310
          do 311 i=1,4
          if(lbc(lvface(i,IS)).lt.1) goto 310
  311     continue
          lbface(1,IS)=nb
  310     continue
          do 302 i=ivbcnd(nb-1)+1,ivbcnd(nb)
          iv=lvbcnd(i)
          lbc(iv)=0
  302     continue
!F!        elseif(kdbcnd(0,nb).eq.kdbuff) then
!F!          do i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F!            iv=lvbcnd(i)
!F!            lbc(iv)=1
!F!          enddo
!F!          do 322 IS=1,nface
!F!          do 321 i=1,4
!F!          if(lbc(lvface(i,IS)).lt.1) goto 322
!F! 321      continue
!F!          lbface(1,IS)=nb
!F! 322      continue
!F!          do 332 i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F!          iv=lvbcnd(i)
!F!          lbc(iv)=0
!F! 332      continue
!F        endif
!F      elseif(kdbcnd(0,nb).eq.kdprdc) then
!F        do 220 IS=1,nface
!F        if(abs(lbfacx(IS)).eq.nb) lbface(1,IS)=lbfacx(IS)
!F  220   continue
!F      endif
  200 continue
!
!F      do 300 nb=1,nbcnd
!F      if(kdbcnd(0,nb).ne.kdprdc) then
!F! --- non-periodic BC:
!F        if(kdbcnd(0,nb).eq.kdintr) then
!F          do 201 i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F          iv=lvbcnd(i)
!F          lbc(iv)=1
!F  201     continue
!F          do 210 IS=1,nface   !NIFACE
!F          IC1=lcface(1,IS)
!F          IC2=lcface(2,IS)
!F          if(lbface(1,IS).ne.0) goto 210
!F          if(lcface(2,IS).gt.ncell) goto 210
!F          if(lacell(IC1).eq.lacell(IC2)) goto 210
!F          do 211 i=1,4
!F          if(lbc(lvface(i,IS)).lt.1) goto 210
!F  211     continue
!F          lbface(1,IS)=nb
!F  210     continue
!F          do 202 i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F          iv=lvbcnd(i)
!F          lbc(iv)=0
!F  202     continue
!F        elseif(kdbcnd(0,nb).eq.kdsld) then
!F          do 401 i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F          iv=lvbcnd(i)
!F          lbc(iv)=1
!F  401     continue
!F          do 410 IS=1,nface
!F          IC1=lcface(1,IS)
!F          IC2=lcface(2,IS)
!F          if(lbface(1,IS).ne.0) goto 410
!F          if(IC2.le.ncell) goto 410
!F!          if(lacell(IC1).eq.lacell(IC2)) goto 410
!F          do 411 i=1,4
!F          if(lbc(lvface(i,IS)).lt.1) goto 410
!F  411     continue
!F          lbface(1,IS)=nb
!F  410     continue
!F          do 402 i=ivbcnd(nb-1)+1,ivbcnd(nb)
!F          iv=lvbcnd(i)
!F          lbc(iv)=0
!F  402     continue
!F        endif
!F      endif
!F  300 continue
!
      do mb=1,nbcnd
      m=0
      do IS=1,nface
        IC1=lcface(1,IS)
        IC2=lcface(2,IS)
        nb=lbface(1,IS)
        if(IC2.gt.ncell.and.nb.eq.0) then
!F        if(mb.eq.nb) then
          m=m+1
!          stop 9997
        endif
      enddo
!F    for FLUENT only
      print*,'           BC No.',mb,' ',trim(SFBOUN(mb)),m
!F      print*,'WWWWWWWWWWWWWWWWWW',mb,'',m,trim(boundName(mb,1)),'====',
!F     &       trim(boundName(mb,2))
      enddo

!
      deallocate(lbc)
!F      deallocate(lbc,lbfacx)
!F      call debug
!F!
!F!///////////////////////////////////////////////////////////////////////
!F      contains
!F!
!F      subroutine debug
!F      use module_debug,only : idebug
!F      if( idebug(2).eq.0 ) return
!F      call printi('lbface/list_facebc',lbface,2,mface)
!F      end subroutine debug
!F!
      end subroutine list_facebc_FLUENT
