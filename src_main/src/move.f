!
!     subroutine output_move_grid
!     subroutine read_move_ini_grid
!     subroutine metric_3ds
!     subroutine metric_CV
!     subroutine list_onvrtx
!     subroutine AXBEQC
!     subroutine dc_metric
!     subroutine metric_speed
!     subroutine piston
!     subroutine 
!     subroutine 
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine output_move_grid(iflmvg,fnamgrid,cord,lacell,lvcell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_io,only       : ifle,ifll,gdScale
      use module_boundary,only : SFBOUN,NBOUND,
     &                           IFFACE,NFBOUN,IBIDX,NBFS
      use module_model, only   : ical_mvmsh
      use module_model,only       : vertex_cen,cell_cen,icon_cv
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: iflmvg
      real*8 ,intent(in) :: cord  (3,MXVRTX_m)
      integer,intent(in) :: lacell(  MXCELL_m)
      integer,intent(in) :: lvcell(8,MXCELL_m)
      character(*),intent(in) :: fnamgrid
!
! --- [local entities]
!
      integer :: ios=0,NAME_I,NBOUND_TEMP
      integer :: iv,i,ic,mb,ical_mvmshx=0
!
!---------------------------------
! --- output grid file
! --- iflmvg='move_grid'
!---------------------------------
!
      if(ical_mvmsh/=0) ical_mvmshx=ical_mvmsh
      open(iflmvg,file=fnamgrid,FORM='unformatted',status='unknown',
     &            iostat=ios)
      if(ios/=0) then
        write(ifle,'(1X,2a)')
     &  '*** Cannot create FrontFlowRed Grid File:',fnamgrid
        call FFRABORT(1,'ERR: in output_move_grid')
      end if
!
      write(iflmvg)   nvrtx,ical_mvmshx,icon_cv
      write(iflmvg) ((cord(i,iv)/gdScale, i=1,3), iv=1,nvrtx)
      write(iflmvg)   ncell
      write(iflmvg)  (lacell(ic),ic=1,ncell)
      write(iflmvg) ((lvcell(i,ic), i=1,8), ic=1,ncell)
      write(iflmvg)   NBOUND
      write(iflmvg)   NBFS
!
!      allocate(IBIDX(0:NBOUND))
      IBIDX=0
      NAME_I=1
      do mb=NAME_I,NBOUND 
      IBIDX(mb)=NFBOUN(mb) 
      enddo 
      do mb=NAME_I,NBOUND 
        IBIDX(mb)=IBIDX(mb)+IBIDX(mb-1) 
      enddo
!
      do mb=NAME_I,NBOUND
        NBOUND_TEMP=mb+(1-NAME_I)
        write(iflmvg) NBOUND_TEMP,NFBOUN(mb),SFBOUN(mb)
        write(iflmvg) mb,IBIDX(mb)
        write(iflmvg) IFFACE(1:4,1:NFBOUN(mb),mb)
      enddo
!
!      do mb=1,NBOUND   !NAME_I,NBOUND
!        write(iflmvg) mb,NFBOUN(mb),SFBOUN(mb)
!        write(iflmvg) IFFACE(1:4,1:NFBOUN(mb),mb)
!      end do
!
      close(iflmvg)
!
!      deallocate(IBIDX)
! --- 
      return
      end subroutine output_move_grid
!     
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_move_ini_grid(cord,lacell,lvcell,
     &           lvface,lbface,lcface,lfcell,LEFACE,listbc)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_io,  only : ifle,ifll,gdScale
      use module_model, only : ical_mvmsh
      use module_metrix,only : idum =>IW2K1
      use module_boundary,only : SFBOUN,NBOUND,
     &                           IFFACE,NFBOUN,IBIDX,NBFS
      use module_model,   only : vertex_cen,cell_cen,icon_cv
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(inout) :: cord  (3,MXVRTX_m)
      integer,intent(inout) :: lacell(  MXCELL_m)
      integer,intent(inout) :: lvcell(8,MXCELL_m)
      integer,intent(inout) :: lvface(4,MXFACE_m)
      integer,intent(inout) :: lbface(2,MXFACE_m)
      integer,intent(inout) :: lcface(2,MXFACE_m)
      integer,intent(inout) :: lfcell(7,MXCELL_m)
      integer,intent(inout) :: LEFACE(5,MXFACE_m)
      integer,intent(inout) :: listbc(4,MXSSFBC_old_m)
!
! --- [local entities]
!
      integer :: i,j,IS,boundID,my_rank,iv,ic,mb,nvrtxx
      integer :: ifl,ios=0,ical_mvmshx,ierr1=0
      integer :: NBOUND_TEMP,NFBOUNx,mbx,icon_cvX
      character(len=80),save :: fnam,SFBOUNx
!
! --- output grid file
!
      fnam='moveinigrid.frontflow'
      ifl=50
      open(ifl,file=fnam,FORM='unformatted',status='unknown',
     &            iostat=ios)
      if(ios/=0) then
        write(ifle,'(1X,a)') 'ERR: read_move_ini_grid:',trim(fnam)
        call FFRABORT(1,'ERR: moveinigrid.frontflow')
      end if
!
!
      read(ifl)   nvrtxx,ical_mvmshx,icon_cvX
      if(icon_cvX/=icon_cv) then
        call FFRABORT(1,'ERR: icon_cvX/=icon_cv in movemesh read')
      endif
      read(ifl) ((cord(i,iv),i=1,3),iv=1,nvrtx)
      read(ifl)   ncell
      read(ifl)  (lacell(ic),ic=1,ncell)
      read(ifl) ((lvcell(i,ic), i=1,8), ic=1,ncell)
      read(ifl)   NBOUND_TEMP
      read(ifl)   NBFS
      cord=cord*gdScale
      NBOUND=NBOUND_TEMP
!
      allocate(NFBOUN(NBOUND),stat=ierr1)
      if(ierr1/=0) call FFRABORT(1,'ERR: allocate NFBOUN')
      allocate(IFFACE(4,NBFS,NBOUND),stat=ierr1)
      if(ierr1/=0) call FFRABORT(1,'ERR: allocate IFFACE')
      allocate(SFBOUN(0:NBOUND),stat=ierr1)
      if(ierr1/=0) call FFRABORT(1,'ERR: allocate IFFACE')
      allocate(IBIDX(0:NBOUND),stat=ierr1)
      IFFACE=0
      SFBOUN=' '
      IBIDX=0
!
      do mb=1,NBOUND
        read(ifl) mbx,NFBOUN(mb),SFBOUN(mb)
        SFBOUN(mb)=TRIM(adjustl(SFBOUN(mb)))  !SFBOUN(mb)
        read(ifl) IFFACE(1:4,1:NFBOUN(mb),mb)
      end do
!
      if(ical_mvmshx/=0) then
        if(.not.(ical_mvmsh/=0)) then
          call FFRABORT(1,'RER:')
        endif
      elseif(ical_mvmshx==0) then
        if(ical_mvmsh/=0) then
          call FFRABORT(1,'RER:')
        endif
      endif
      if(ical_mvmshx/=0) then
        read(ifl) ((lvface(i,iv),i=1,4),iv=1,nface)
        read(ifl) ((lbface(i,iv),i=1,2),iv=1,nface)
        read(ifl) ((lcface(i,iv),i=1,2),iv=1,nface)
        read(ifl) ((lfcell(i,iv),i=1,7),iv=1,ncell)
        read(ifl) ((LEFACE(i,iv),i=1,5),iv=1,nface)
        read(ifl) ((listbc(i,iv),i=1,4),iv=1,nssfbc_old)
      endif
      close(ifl)
!
      return
      end subroutine read_move_ini_grid 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_3ds
     & (icall,cord,lacell,lvcell,
     &  lvface,lbface,lcface,lfcell,
     &  area,volume,gface,gcell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_dimension
      use module_io,only       : ifle,ifll
      use module_metrix,only   : lvrtx =>IW2K1
      use module_model,   only : ical_mvmsh
      
!
! 1. Calculate metrics for each face & cell
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: icall
      real*8 ,intent(in)     :: cord  (3,MXVRTX_m)
      integer,intent(in)     :: lacell(  MXCELL_m)
      integer,intent(in)     :: lvcell(8,MXCELL_m) 
      integer,intent(in)     :: lvface(4,MXFACE_m)      
      integer,intent(in)     :: lbface(2,MXFACE_m)
      integer,intent(in)     :: lcface(2,MXFACE_m)
      integer,intent(in)     :: lfcell(7,MXCELL_m)
      real*8 ,intent(out)    :: area  (4,MXFACE_m)
      real*8 ,intent(out)    :: volume(  MXCELL_m)
      real*8 ,intent(out)    :: gface (3,MXFACE_m)
      real*8 ,intent(out)    :: gcell (3,MXCELL_m)
      integer,intent(out)    :: ierror
!
! --- [local entities] 
!
!
      integer :: i,j,k,m,n,k1,k2,nf,ie
      integer :: IC,IS,IV,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE,ierr1=0
!
      
      integer,parameter :: itr(2,2)=
     &   reshape(source=(/1,3,2,4/), shape=(/2,2/))
      real*8 ,parameter :: r1p6=1.d0/6.d0
      real*8 ,parameter :: r1p3=1.d0/3.d0
      real*8 ,parameter :: r1p4=0.25d0
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8  :: r(3,4),r1(3),r2(3),r3(3),r4(3),rg(4),aa(3,4)
      real*8  :: ss1,ss2,sx,sy,sz,vv
      real*8  :: sx1,sy1,sz1,sx2,sy2,sz2,vol,dum1,tolvol
      real*8            :: SML2
      integer,parameter :: lvfcel(4,6,4)=reshape( source=
     &  (/1,3,2,0, 2,3,4,0, 3,1,4,0, 4,1,2,0, 0,0,0,0, 0,0,0,0,
     &    1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0, 0,0,0,0,
     &    1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2,0, 4,5,6,0, 0,0,0,0,
     &    1,5,8,4, 2,3,7,6, 1,2,6,5, 3,4,8,7, 1,4,3,2, 5,6,7,8/),
     &  shape=(/4,6,4/) )
      integer :: lvf (4,6,4)
      
!
      ierror=0
      SML2=0.d0
      if(ical_mvmsh==4.or.ical_mvmsh==5) SML2=1.d-15
!
! --- initializing array of CVCENT
!
      do k=1,4
      do 112 j=1,6
      do 111 i=1,4
      lvf(i,j,k)=lvfcel(i,j,k)
 111  enddo
      if(lvf(4,j,k).lt.1) lvf(4,j,k)=8
 112  enddo
      enddo
!
!----------------------------------------------------
! --- CELL METRIX CALCULATION
!----------------------------------------------------
!
      allocate(lvrtx(4,MXFACE_m),stat=ierr1)
      if(ierr1/=0) call FFRABORT(1,'ERR:allocate lvrtx')
      lvrtx(:,:)=0
      if(ierr1.ne.0) call FFRABORT(1,'ERR: allocate lvrtx')
!
      if((ical_mvmsh==4.or.ical_mvmsh==5).and.icall==2) goto 5000
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
      do ILV=1,4
      do 102 i=1,3
      r(i,ILV)=cord(i,lvrtx(ILV,IS))
  102 enddo
      enddo
!
!--< 1.1 cell-face area (multiplied by 2) >--
!
! --- assert: Outside-ward for IC1 cell
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
!
      if( area(4,IS).le.0.d0 ) then
          write(ifle,*) '### error-1-1 : invalid coordinate'
          write(ifle,*) 'area of face is not positive'
          write(ifle,*) 'face no. =',IS
          write(ifle,*) 'vertices of face =',(lvrtx(i,IS),i=1,ie)
          call FFRABORT(1,'ERR: area<0')
      endif
!----------------------------
!--< 1.2 cell-face center >--
!----------------------------
      do 120 k1=1,2
      k2=3-k1
      do 121 i=1,3
      r1(i)=r(i,itr(1,k2))-r(i,itr(k1,k1))
      r2(i)=r(i,itr(2,k2))-r(i,itr(k1,k1))
      dum1 =r(i,itr(1,k2))+r(i,itr(2,k2))
      r3(i)=r(i,itr(k1,k1))+dum1
      r4(i)=r(i,itr(k2,k1))+dum1
  121 enddo
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
  122 enddo
      rg(4)=rg(4)+ss1+ss2
  120 enddo
!
      if(rg(4).le.0.d0 ) then
          write(ifle,*) '### error-1-2 : invalid coordinate'
          write(ifle,*) 'area of face is not positive'
          write(ifle,*) 'face no. =',IS
          write(ifle,*) 'vertices of face =',(lvrtx(i,IS),i=1,ie)
          call FFRABORT(1,'ERR: rg<0')
      endif
      dum1=r1p3/(rg(4)+SML2)
      do 123 i=1,3
! --- Face center's coordinates
      gface(i,IS)=rg(i)*dum1
  123 enddo
!
  100 enddo
!
      do 200 IC=1,ncell
!
!--< 2.1 cell volume >--
!
      vol=0.d0
      do j=1,4
      do 201 i=1,3
      aa(i,j)=0.d0
  201 continue
      enddo
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
  211 enddo
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
         write(ifle,'(1x,a)') 'ERR: Nagitive volume'
         write(ifle,1100) IC,(lvcell(i,IC),i=1,8)
         do i=1,8
         if(lvcell(i,IC).gt.0) then
            write(ifle,1200) i,lvcell(i,IC),(cord(j,lvcell(i,IC)),j=1,3)
         endif
         enddo
         do ILS=1,lfcell(7,IC)
         IS=lfcell(ILS,IC)
         write(ifle,1000) ILS,(lvface(i,IS),i=1,4)
         write(ifle,1000) ILS,(lvcell(lvf(i,ILS,4),IC),i=1,4)
         enddo
         write(ifle,*) '### error-1-3 : invalid coordinate'
         write(ifle,*) 'volume of cell is not positive'
         write(ifle,*) 'cell no. =',IC,volume(IC)
         call FFRABORT(1,'ERR: ')
      endif
!
!-----------------------
!--< 2.2 cell center >--
!-----------------------
!
      do 220 i=1,3
      aa(i,i)=aa(i,i)+vol
      aa(i,4)=r1p4*aa(i,4)
  220 enddo
!
      do 221 i=1,3
      if(aa(i,i).le.0.d0) then
        write(ifle,*) '### error-1-4 : invalid coordinate'
        write(ifle,*) 'volume of cell is not positive'
        write(ifle,*) 'cell no. =',IC,volume(IC)
        call FFRABORT(1,'ERR: ')
      endif
      do 222 j=i+1,4
      aa(i,j)=aa(i,j)/(aa(i,i)+SML2)
  222 enddo
      do k=i+1,3
      do 223 j=i+1,4
      aa(k,j)=aa(k,j)-aa(i,j)*aa(k,i)
  223 enddo
      enddo
  221 enddo
!
      gcell(3,IC)=aa(3,4)
      gcell(2,IC)=aa(2,4)- aa(2,3)*gcell(3,IC)
      gcell(1,IC)=aa(1,4)-(aa(1,2)*gcell(2,IC)+aa(1,3)*gcell(3,IC))
!
  200 enddo
!----------------------------------------------
!-< 3. Normalizing face vector to unit vector>-
!----------------------------------------------
      do 300 IS=1,nface
      dum1=1.d0/area(4,IS)
      do 301 i=1,3
      area(i,IS)=area(i,IS)*dum1
  301 continue
      area(4,IS)=XHALF*area(4,IS)
  300 continue
!----------------------------------------------------------------------
      deallocate(lvrtx)
      return
!----------------------------------------------------------------------
 5000 continue
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      do 500 IS=1,nface
!
      do 501 ILV=1,4
      lvrtx(ILV,IS)=lvface(ILV,IS)
      rg(ILV)=0.d0
  501 continue
! --- if/not 3p-face?
      if(lvrtx(4,IS).lt.1) then
        lvrtx(4,IS)=lvrtx(1,IS)
      endif
      do 502 ILV=1,4
      do i=1,3
      r(i,ILV)=cord(i,lvrtx(ILV,IS))
      enddo
  502 enddo
!
!--< 1.1 cell-face area (multiplied by 2) >--
!
! --- assert: Outside-ward for IC1 cell
!
      do 510 i=1,3
      r1(i)=r(i,3)-r(i,1)
      r2(i)=r(i,4)-r(i,2)
  510 continue
      sx=r1(2)*r2(3)-r1(3)*r2(2)
      sy=r1(3)*r2(1)-r1(1)*r2(3)
      sz=r1(1)*r2(2)-r1(2)*r2(1)
      area(4,IS)=sqrt(sx*sx+sy*sy+sz*sz)
      if(area(4,IS)<-SML2) 
     &  call FFRABORT(1,'ERR:area of face is not positive')
!----------------------------
!--< 1.2 cell-face center >--
!----------------------------
      do 520 k1=1,2
      k2=3-k1
      do 521 i=1,3
      r1(i)=r(i,itr(1,k2))-r(i,itr(k1,k1))
      r2(i)=r(i,itr(2,k2))-r(i,itr(k1,k1))
      dum1 =r(i,itr(1,k2))+r(i,itr(2,k2))
      r3(i)=r(i,itr(k1,k1))+dum1
      r4(i)=r(i,itr(k2,k1))+dum1
  521 enddo
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
      do 522 i=1,3
      rg(i)=rg(i)+ss1*r3(i)+ss2*r4(i)
  522 enddo
      rg(4)=rg(4)+ss1+ss2
  520 enddo
!
      if(rg(4).le.SML2) then
        if(rg(4)<-SML2) then
          call FFRABORT(1,'ERR:area of face is not positive')
        endif
        do i=1,3
        dum1=0.d0
        do ILV=1,4
        dum1=dum1+r(i,ILV)
        enddo
        gface(i,IS)=dum1*0.25d0
        enddo
      else
        dum1=r1p3/(rg(4)+SML2)
        do 523 i=1,3
! --- Face center's coordinates
        gface(i,IS)=rg(i)*dum1
  523   enddo
      endif
  500 enddo
!
      do 600 IC=1,ncell
!
!--< 2.1 cell volume >--
!
      vol=0.d0
      do j=1,4
      do 601 i=1,3
      aa(i,j)=0.d0
  601 continue
      enddo
!
      do 610 ILS=1,lfcell(7,IC)
      IS=lfcell(ILS,IC)
      do 611 i=1,3
      r1(i)=cord(i,lvrtx(1,IS))
     &     +cord(i,lvrtx(2,IS))
     &     +cord(i,lvrtx(3,IS))
      if(lvrtx(4,IS).gt.0) then
        r1(i)=r1(i)+cord(i,lvrtx(4,IS))
      end if
  611 enddo
      vv=(r1(1)*area(1,IS)+r1(2)*area(2,IS)+r1(3)*area(3,IS))*area(4,IS)
!-------------------------------------------------------------
! --- Distinguish IC1 and IC2 for not getting negative volume.
!-------------------------------------------------------------
      if(IC.eq.lcface(1,IS)) then
! --- IC1 cell
        vol=vol+vv
        do 612 i=1,3
        do 613 j=1,3
        aa(i,j)=aa(i,j)+gface(i,IS)*area(j,IS)*area(4,IS)
  613   continue
        aa(i,4)=aa(i,4)+gface(i,IS)*vv
  612   continue
      else
! --- IC2 cell
        vol=vol-vv
        do 614 i=1,3
        do 615 j=1,3
        aa(i,j)=aa(i,j)-gface(i,IS)*area(j,IS)*area(4,IS)
  615   continue
        aa(i,4)=aa(i,4)-gface(i,IS)*vv
  614   continue
      endif
  610 continue
!
      vol=r1p4*vol
      volume(IC)=r1p6*vol
      tolvol=volume(IC)
      if(tolvol.le.-SML2) call FFRABORT(1,'ERR:nagative volume') 
      if(tolvol.le.SML2) volume(IC)=0.d0
!-----------------------
!--< 2.2 cell center >--
!-----------------------
!
      do 620 i=1,3
      aa(i,i)=aa(i,i)+vol
      aa(i,4)=r1p4*aa(i,4)
  620 enddo
!
      if(abs(tolvol).le.SML2) then
        r1(:)=0.d0
        dum1=0.d0
        do ILS=1,lfcell(7,IC)
        IS=lfcell(ILS,IC)
        dum1=dum1+1.d0
        do i=1,3
        r1(i)=r1(i)+gface(i,IS)
        enddo
        enddo
        gcell(:,IC)=r1(:)/dum1
      else
        do 621 i=1,3
        do 622 j=i+1,4
        aa(i,j)=aa(i,j)/(aa(i,i)+SML2)
  622   enddo
        do k=i+1,3
        do 623 j=i+1,4
        aa(k,j)=aa(k,j)-aa(i,j)*aa(k,i)
  623   enddo
        enddo
  621   enddo
        gcell(3,IC)=aa(3,4)
        gcell(2,IC)=aa(2,4)- aa(2,3)*gcell(3,IC)
        gcell(1,IC)=aa(1,4)-(aa(1,2)*gcell(2,IC)+aa(1,3)*gcell(3,IC))
      endif
!
  600 enddo
!----------------------------------------------
!-< 3. Normalizing face vector to unit vector>-
!----------------------------------------------
      do 700 IS=1,nface
      area(4,IS)=XHALF*area(4,IS)
  700 continue
      
!      endif
!
      deallocate(lvrtx)
!
      return
!
 1000    format(4x,'MSG: ',I8,' face: ',4I8)
 1100    format(4x,'MSG: IC= ',I8,' IV= ',8I8)
 1200    format(4x,'MSG: cord: ',2I8,3F15.6)

!
      end subroutine metric_3ds
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_CV
     & (LBC_SSF,lcface,lbface,lacell,
     &  cord,area,volume,gface,gcell,
     &  LVEDGE,LEFACE,LVRTCV,listbc,MAT_CV,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,
     &  ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_io,only          : ifle,ifll,cntlnam
      use module_metrix,only      : lvrtx=>IW2K1
      use module_metrix,only      : ip   =>iwork1
      use module_metrix,only      : nclv =>iwork2
      use module_boundary,only    : kdbcnd,kdintr,kdprdc,kdsld,idis
!
      implicit none
! 1. Calculate metrics for each face & cell
!
! -- [dummy arguments]
!
      integer,intent(in)     :: lcface(2,mxface_m)
      integer,intent(in)     :: lbface(2,mxface_m)
      integer,intent(in)     :: lacell(  mxcell_m)
      integer,intent(in)     :: LEFACE(5,mxface_m)    !
      integer,intent(in)     :: LVEDGE(2,MXCVFAC)     !
      integer,intent(in)     :: LVRTCV(  MXALLCV)     !
      integer,intent(in)     :: listbc(4,MXSSFBC_old_m)   !
      integer,intent(in)     :: LBC_SSF( mxssfbc)
      integer,intent(in)     :: MAT_CV  (MXALLCV)
!
      real*8 ,intent(in)     :: cord  (3,mxvrtx_m)
      real*8 ,intent(in)     :: area  (4,mxface_m)
      real*8 ,intent(in)     :: volume(  mxcell_m)
      real*8 ,intent(in)     :: gface (3,mxface_m)
      real*8 ,intent(in)     :: gcell (3,mxcell_m)
      real*8 ,intent(out)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(out)    :: SFCENT(3,MXCVFAC)
      real*8 ,intent(out)    :: CVVOLM(  MXALLCV)
      real*8 ,intent(out)    :: CVCENT(3,MXALLCV)
      integer,intent(out)    :: ierror
!
! --- [local entities]
!
      integer :: idmax,LMAT0,SMLV
      integer :: KMAT(2),IDC(2)
      integer :: i,j,k,l,m,n,kd,kdv,kdt,kdy,kdk,kdp,nbp,ISP
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB
      integer :: ICVA,ICVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
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
      integer :: ierr1=0,IE2,ICVBP2,ICFA1,ICFB1,ICFA2,ICFB2,IDCA1,IDCB1
      integer :: ICFA,ICFB,IDCA,IDCB,ISSFA,ISSFB,ILE1,IE1,ICVA1,ICVB1,
     &          IVA1,IVB2,IVB1,ILE2,ICVA2,ICVB2,IVA2,IA,IB,IVAP2,IVBP2,
     &           ICVAP2,IDCA2,IDCB2,ISSFA1,ISSFB1,ISSFA2,ISSFB2
      integer :: ISSF,IBFL,ICFL
      integer :: IVA,IVB
      integer :: NDUP=0
      real*8 ,parameter :: r1p6=1.d0/6.d0
      real*8 ,parameter :: r1p3=1.d0/3.d0
      real*8 ,parameter :: r1p4=0.25d0
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
      allocate(lvrtx(4,mxface_m),stat=ierr1)
      if(ierr1.ne.0) 
     &  call FFRABORT(1,'stop at allocating lvrtx(:,:) in metric_CV')
      allocate(ip(0:nedge),stat=ierr1) 
      if(ierr1.ne.0) 
     &  call FFRABORT(1,'stop at allocating ip(:) in metric_CV')
      allocate(nclv(4*nface),stat=ierr1) 
      if(ierr1.ne.0) 
     &  call FFRABORT(1,'stop at allocating nclv(:) in metric_CV')
!
!------------------------------------
! --- initializing array of CVCENT
!------------------------------------
!
      CVCENT=0.D0
      CVVOLM=0.D0
      SFAREA=0.D0
      SFCENT=0.D0
      
!----------------------------------------------------
! --- CONTROL VOLUME (CV) METRIX CALCULATION
!----------------------------------------------------
! --- make-up egde-face connectivity
!----------------------------------------------------
      lvrtx(1:4,:)=LEFACE(1:4,:)
      call list_onvrtx(4,mxface_m,nedge,nface,lvrtx,ip,nclv,idmax)
!=====================================================================
! --- <1> define inner Sub-Face (SF) vector (CV face)
!=====================================================================
      DO 400 IE=1,nedge
      ICVA=LVEDGE(1,IE)
      ICVB=LVEDGE(2,IE)
      IVA=MAT_CV(ICVA)
      IVB=MAT_CV(ICVB)
      if(IVA>NCVIN.or.IVB>NCVIN) cycle
      IVA=LVRTCV(IVA)
      IVB=LVRTCV(IVB)
      SMSSFX=ZERO
      SMSSFY=ZERO
      SMSSFZ=ZERO
      SUMSWX=ZERO
      SUMSWY=ZERO
      SUMSWZ=ZERO
      SUMSSV=ZERO
      smssf=ZERO
      SSVXA=ZERO
      SSVYA=ZERO
      SSVZA=ZERO
      SSVXB=ZERO
      SSVYB=ZERO
      SSVZB=ZERO
      do 410 IEE=ip(IE-1)+1,ip(IE)
      IS=nclv(IEE)
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      IDC(1)=IC1
      IDC(2)=IC2
      DO 450 J=1,2
      IC=lcface(J,IS)
      IF(IC.GT.NCELL) goto 450
      DO 420 I=1,3
      gfc(I)=gface(I,IS)
      gc(I)=gcell(I,IC)
      gec(I)=XHALF*(cord(I,IVA)+cord(I,IVB))
      VEDGE(I)=XHALF*(cord(I,IVB)-cord(I,IVA))


      VFE(I)=gfc(I)-gec(I)
      VCE(I)=gc(I)-gec(I)
      GSSF(I)=r1p3*(gfc(I)+gc(I)+gec(I))
      GSSA(I)=r1p4*(gfc(I)+gc(I)+gec(I)+cord(I,IVA))
      GSSB(I)=r1p4*(gfc(I)+gc(I)+gec(I)+cord(I,IVB))
 420  CONTINUE
      DIRC=DSQRT(VEDGE(1)*VEDGE(1)+VEDGE(2)*VEDGE(2)+VEDGE(3)*VEDGE(3))
      UEDGE(1)=VEDGE(1)/DIRC
      UEDGE(2)=VEDGE(2)/DIRC
      UEDGE(3)=VEDGE(3)/DIRC
      CALL AXBEQC(VFE,VCE,SSF,SSFS)
      DIRC=SSF(1)*UEDGE(1)
     &    +SSF(2)*UEDGE(2)
     &    +SSF(3)*UEDGE(3)
      IF(DIRC.LT.ZERO) THEN
        SSF(1)=-SSF(1)
        SSF(2)=-SSF(2)
        SSF(3)=-SSF(3)
      ENDIF
      SMSSFX=SMSSFX+SSF(1)
      SMSSFY=SMSSFY+SSF(2)
      SMSSFZ=SMSSFZ+SSF(3)
      SUMSWX=SUMSWX+GSSF(1)*SSFS
      SUMSWY=SUMSWY+GSSF(2)*SSFS
      SUMSWZ=SUMSWZ+GSSF(3)*SSFS
      SSV=r1p3*(SSF(1)*VEDGE(1)
     &         +SSF(2)*VEDGE(2)
     &         +SSF(3)*VEDGE(3))
      SUMSSV=SUMSSV+SSV
      smssf=smssf+SSFS
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
 450  ENDDO
      IF(IDC(1).gt.ncell.and.IDC(2).gt.ncell) then
	write(ifle,'(1x,a,I10)') 
     &   'ERR: Two Cells are Dumy Cell on face ',IS
        call FFRABORT(1,'ERR:')
      endif
 410  ENDDO
      SUMSSF=SQRT(SMSSFX*SMSSFX
     &           +SMSSFY*SMSSFY
     &           +SMSSFZ*SMSSFZ)
      SFAREA(1,IE)=SMSSFX/(SUMSSF+SML)
      SFAREA(2,IE)=SMSSFY/(SUMSSF+SML)
      SFAREA(3,IE)=SMSSFZ/(SUMSSF+SML)
      SFAREA(4,IE)=SUMSSF
      SFCENT(1,IE)=SUMSWX/(SMSSF+SML)
      SFCENT(2,IE)=SUMSWY/(SMSSF+SML)
      SFCENT(3,IE)=SUMSWZ/(SMSSF+SML)
      
      CVVOLM(ICVA)=CVVOLM(ICVA)+SUMSSV
      CVCENT(1,ICVA)=CVCENT(1,ICVA)+SSVXA
      CVCENT(2,ICVA)=CVCENT(2,ICVA)+SSVYA
      CVCENT(3,ICVA)=CVCENT(3,ICVA)+SSVZA
      CVVOLM(ICVB)=CVVOLM(ICVB)+SUMSSV
      CVCENT(1,ICVB)=CVCENT(1,ICVB)+SSVXB
      CVCENT(2,ICVB)=CVCENT(2,ICVB)+SSVYB
      CVCENT(3,ICVB)=CVCENT(3,ICVB)+SSVZB
 400  ENDDO
!=====================================================================
! --- Volume of CV 
!=====================================================================
      DO 480 IV=1,NCV
      CVCENT(1,IV)=CVCENT(1,IV)/(CVVOLM(IV)+SML)
      CVCENT(2,IV)=CVCENT(2,IV)/(CVVOLM(IV)+SML)
      CVCENT(3,IV)=CVCENT(3,IV)/(CVVOLM(IV)+SML)
 480  ENDDO
!
!=====================================================================
! --- <2> Defining SSF on BC 
!=====================================================================
!
      do ISSF=1,NSSFBC_OLD
      IE=listbc(1,ISSF)
      IS=listbc(2,ISSF)
      IVA=listbc(3,ISSF)
!      IBFL=listbc(4,ISSF)
!      ICFL=LBC_SSF(IBFL)
      ICFL=listbc(4,ISSF)
      ICVA=LVEDGE(1,IE)
      ICVB=LVEDGE(2,IE)
      IVA1=LVRTCV(ICVA)
      IVB1=LVRTCV(ICVB)
!      if(IVB1==0.or.IVA1==0) cycle
      DO 630 I=1,3
      gfc(I)=gface(I,IS)
      gvA(I)=cord(I,IVA)
      gec(I)=XHALF*(cord(I,IVA1)+cord(I,IVB1))
      VEDGE(I)=area(I,IS)
      VFE(I)=gfc(I)-gec(I)
      VCEA(I)=gvA(I)-gec(I)
      GSSA(I)=r1p3*(gfc(I)+gvA(I)+gec(I))
 630  enddo

      CALL AXBEQC(VFE,VCEA,SSFA,SSFSA)
      DIRC=SSFA(1)*VEDGE(1)
     &    +SSFA(2)*VEDGE(2)
     &    +SSFA(3)*VEDGE(3)
      IF(DIRC.LT.ZERO) THEN
          SSFA(1)=-SSFA(1)
          SSFA(2)=-SSFA(2)
          SSFA(3)=-SSFA(3)
      ENDIF
      SFAREA(4,ICFL)=SFAREA(4,ICFL)+SSFSA
      SFAREA(1,ICFL)=SFAREA(1,ICFL)+SSFA(1)
      SFAREA(2,ICFL)=SFAREA(2,ICFL)+SSFA(2)
      SFAREA(3,ICFL)=SFAREA(3,ICFL)+SSFA(3)
      SFCENT(1,ICFL)=SFCENT(1,ICFL)+GSSA(1)*SSFSA
      SFCENT(2,ICFL)=SFCENT(2,ICFL)+GSSA(2)*SSFSA
      SFCENT(3,ICFL)=SFCENT(3,ICFL)+GSSA(3)*SSFSA
      enddo
!
!
      do IBFL=1,NSSFBC
      ICFL=LBC_SSF(IBFL)
      SSFS=SFAREA(1,ICFL)*SFAREA(1,ICFL)
     &    +SFAREA(2,ICFL)*SFAREA(2,ICFL)
     &    +SFAREA(3,ICFL)*SFAREA(3,ICFL)
      SSFS=DSQRT(SSFS)
      SSV=SFAREA(4,ICFL)
!      SFAREA(4,ICFL)=SSFS
      SFAREA(1,ICFL)=SFAREA(1,ICFL)/SSFS
      SFAREA(2,ICFL)=SFAREA(2,ICFL)/SSFS
      SFAREA(3,ICFL)=SFAREA(3,ICFL)/SSFS
      SFCENT(1,ICFL)=SFCENT(1,ICFL)/SSV
      SFCENT(2,ICFL)=SFCENT(2,ICFL)/SSV
      SFCENT(3,ICFL)=SFCENT(3,ICFL)/SSV
      enddo
!

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
      subroutine list_onvrtx(mv,mcell,nvrtx,ncell,
     &                       lvcell,ip,nclv,icmax)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_metrix,only      : kclv=>iwork3
      use module_metrix,only      : iq  =>iwork4
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
      if(ierr1.ne.0)  
     &  call FFRABORT(1,'stop at allocating kclv(:) in list_onvrtx')
      allocate(iq(nvrtx),stat=ierr1) 
      if(ierr1.ne.0)   
     &  call FFRABORT(1,'stop at allocating iq(:) in list_onvrtx')
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
      C(3)=XHALF*(A(1)*B(2)-A(2)*B(1))
      C(2)=XHALF*(A(3)*B(1)-A(1)*B(3))
      C(1)=XHALF*(A(2)*B(3)-A(3)*B(2))
      D=DSQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))
!
      RETURN
      END subroutine AXBEQC
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_metric
     & (LVEDGE,LBC_SSF,LCYCSF,SFCENT,SFAREA,CVCENT,CVVOLM,ICALL)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      USE module_dimension
      use module_hpcutil
      use module_boundary,only : set_rotmtrx,nbcnd,
     &                           kdbcnd,kdprdc,kdsymm,kdolet,kdilet,
     &                           kdsld,idis,kdintr,LBC_INDEX,MAT_BCIDX
     &                           ,kdbuff,kdovst
      use module_io,only       : ifle,ifll
!
! 1. Set metric component at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: ICALL
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  mxssfbc)
      integer,intent(in)    :: LCYCSF    (  mxssfbc)
      real*8 ,intent(inout) :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(inout) :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(inout) :: CVCENT    (3,MXALLCV)
      real*8 ,intent(inout) :: CVVOLM    (  MXALLCV)
!
! --- [local entities]
!
      real*8  :: aa(3,3)
      real*8  :: r1(3),r2(3),prd,dx,dy,dz,dlvect
      real*8  :: alf,dl,dum1
      real*8,parameter :: SML=1.d-15
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: ICV,IDC,IBFL,nb,IIMAT,kd,ICFL,IBFS,IBFE
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
!---------------------------------
!-< 1. Set values at dummy cell >-
!---------------------------------
      DO 1000 nb=1,nbcnd    !IBF=1,nssfbc
      IIMAT=MAT_BCIDX(nb,1)
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
!------------------------------
!--< 1.1 periodic boundary  >--
!------------------------------
      if(kd==kdprdc) then
        if(idis(nb)==0) then
          call set_rotmtrx(nbcnd,kd,nb,aa)
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          CVVOLM(IDC)=CVVOLM(ICVP)
          CVVOLM(IDCP)=CVVOLM(ICV)
          do 101 i=1,3
          r1(i)=(CVCENT(i,ICV)-SFCENT(i,ICFL))
          r2(i)=(CVCENT(i,ICVP)-SFCENT(i,ICFP))
  101     continue
          do 102 i=1,3
          CVCENT(i,IDC)=SFCENT(i,ICFL)
     &                +(aa(i,1)*r2(1)
     &                 +aa(i,2)*r2(2)
     &                 +aa(i,3)*r2(3))
          CVCENT(i,IDCP)=SFCENT(i,ICFP)
     &                 +(aa(1,i)*r1(1)
     &                  +aa(2,i)*r1(2)
     &                  +aa(3,i)*r1(3))
  102     continue
          enddo
        elseif(idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          CVVOLM(IDC)=CVVOLM(ICV)
          do i=1,3
          r1(i)=(CVCENT(i,ICV)-SFCENT(i,ICFL))
          enddo
          prd=r1(1)*SFAREA(1,ICFL)
     &       +r1(2)*SFAREA(2,ICFL)
     &       +r1(3)*SFAREA(3,ICFL)
          prd=prd+prd

          do i=1,3
          CVCENT(i,IDC)=SFCENT(i,ICFL)+(r1(i)-prd*SFAREA(i,ICFL))
          enddo
          enddo
        endif
      elseif(kd==kdsld.or.kd==kdovst) then
        do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          CVVOLM(IDC) =CVVOLM(ICV)
          r1(:)=SFCENT(:,ICFP)-CVCENT(:,ICVP)
          prd=dsqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
          do i=1,3
          CVCENT(i,IDC)=SFCENT(i,ICFL)+prd*SFAREA(i,ICFL)
          enddo
        enddo
!------------------------------
!--< 1.2 symmetric boundary >--
!------------------------------
      elseif(kd.eq.kdsymm) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        CVVOLM(IDC)=CVVOLM(ICV)
        do 201 i=1,3
        r1(i)=CVCENT(i,ICV)-SFCENT(i,ICFL)
  201   continue
        prd=r1(1)*SFAREA(1,ICFL)
     &     +r1(2)*SFAREA(2,ICFL)
     &     +r1(3)*SFAREA(3,ICFL)
        prd=prd+prd
        do 202 i=1,3
        CVCENT(i,IDC)=SFCENT(i,ICFL)+(r1(i)-prd*SFAREA(i,ICFL))
  202   continue
        enddo
      elseif(kd.eq.kdolet) then
!---------------------------
!--< 1.3 outlet boundary >--
!---------------------------
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        CVVOLM(IDC)=0.d0
        do 302 i=1,3
        CVCENT(i,IDC)=SFCENT(i,ICFL)
  302   continue
        enddo
      elseif(kd.eq.kdintr.or.kd.eq.kdbuff) then !???
        do IBFL=IBFS,IBFE
        ICFP=LCYCSF(IBFL)
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        CVVOLM(IDC)=CVVOLM(ICVP)
        CVVOLM(IDCP)=CVVOLM(ICV)
        do i=1,3
        CVCENT(i,IDC)=CVCENT(i,ICVP)
        CVCENT(i,IDCP)=CVCENT(i,ICV)
        enddo
        enddo
!       
!--< 1.3 other boundary (incluing "interface")>--
!
      else
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        CVVOLM(IDC)=0.d0
        do 301 i=1,3
        CVCENT(i,IDC)=SFCENT(i,ICFL)
  301   continue
        enddo
      endif
!
 1000 continue
!
! --- 
!      do 2000 IBF=1,nssfbc   !ICF=1,NCVFAC
!      ICFL=LBC_SSF(IBFL)
!      ICVA=LVEDGE(1,ICFL)
!      ICVB=LVEDGE(2,ICFL)
!      if(CVCENT(1,ICVB).eq.0.d0.and.
!     &   CVCENT(2,ICVB).eq.0.d0.and.
!     &   CVCENT(3,ICVB).eq.0.d0) then
!        write(ifle,'(1X,a)') 
!     &   ' ### ERR: Dummy Cell coordinates are NOT defined'
!        write(ifle,'(1x,a)') ' Please Contact FrontFlow developer'
!        call FFRABORT(1,'ERR: ')
!      endif
! 2000 continue
!---------------------
! --- check all CV
!---------------------
!      ICFP=0
!      do 3000 ICF=1,NCVFAC
!      if(ICF.gt.NEDGE) goto 3000
!      if(SFCENT(1,ICF).eq.0.d0.and.
!     &   SFCENT(2,ICF).eq.0.d0.and.
!     &   SFCENT(3,ICF).eq.0.d0) then
!        ICFP=ICFP+1
!      endif
! 3000 continue
!      if(ICFP.ge.2) then
!        write(ifle,*) 
!     &  ' ### ERR: CV Sub Face coordinates are NOT defined'
!        write(ifle,*) ' Please Contact FrontFlow developer'
!        stop 'stop at dc_metric'
!      endif
!
!      ICFP=0
!      do 3010 ICF=1,NCVFAC
!      if(abs(SFAREA(4,ICF)).le.SML) then
!        ICFP=ICFP+1
!      endif
! 3010 continue
!      if(ICFP.gt.0) then
!        write(ifle,*) 
!     &   ' ### ERR: CV Sub Face Area is ZERO N=',ICFP
!        write(ifle,*) ' Please Contact FrontFlow developer'
!        stop 'stop at dc_metric'
!      endif
!
!      do 4000 ICFL=1,NCVFAC
!      ICVA=LVEDGE(1,ICFL)
!      ICVB=LVEDGE(2,ICFL)
!      if(ICVB.lt.ICVA)then
!        call FFRABORT(1,'ERR: ICVB < ICVA')
!      endif
!      dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
!      dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
!      dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!      alf=(dx*SFAREA(1,ICFL)
!     &    +dy*SFAREA(2,ICFL)
!     &    +dz*SFAREA(3,ICFL))
!      if(alf.le.0.d0) then
!        write(ifle,'(a)') 
!     &  'WRN: Opposite Dir. between Edge vector and SF vector'
!        write(ifle,'(5x,a,I4,a,I10,E10.4)')
!     &   'ICALL=',ICALL,' At SF=',ICFL,alf
!        write(ifle,'(5x,a,I10,a,I10)') 
!     &   'SF between ICVA=',ICVA,'and ICVB=',ICVB
!      endif 
! 4000 continue
!
      return
!
      end subroutine dc_metric
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine metric_speed(iter,MAT_CFIDX,
     &           SFCENT_OLD,SFCENT,SFAREA,xta,dt)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_dimension
      use module_io    ,  only : ifle,ifll,gdScale
      use module_constant
      use module_Euler2ph,only : ieul2ph
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: iter,MAT_CFIDX( 0:MXMAT)
      real*8 ,intent(in)  :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)  :: SFCENT_OLD(3,MXCVFAC)
      real*8 ,intent(in)  :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(out) :: xta       (  MXCVFAC)
      real*8 ,intent(in)  :: dt
!
! --- [local entities]
!
      real*8  :: dum1,dum2,fct,dx(3)
      integer :: lvrtx(4)
      integer :: i,j,n,IIMAT,ICFS,ICFE,ICFL
!
!-------------------------------------
! 1. Calculate metrics of grid speed
!      integral( dX/dt*area )dt
!-------------------------------------
!
      fct=1.d0/(24.d0*dt)
      if(dt<0.d0) return
!
! --- 
!
      DO IIMAT=1,NMAT
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      DO ICFL=ICFS,ICFE
        dx(1)=SFCENT(1,ICFL)-SFCENT_OLD(1,ICFL)
        dx(2)=SFCENT(2,ICFL)-SFCENT_OLD(2,ICFL)
        dx(3)=SFCENT(3,ICFL)-SFCENT_OLD(3,ICFL)
        dum1=dx(1)*SFAREA(1,ICFL)
     &      +dx(2)*SFAREA(2,ICFL)
     &      +dx(3)*SFAREA(3,ICFL)
!
!        dum2=dsqrt(dx(1)**2+dx(2)**2+dx(3)**2)
!        dx(1)=dx(1)/(dum2+SML)
!        dx(2)=dx(2)/(dum2+SML)
!        dx(3)=dx(3)/(dum2+SML)
!
        xta(ICFL)=-dum1*SFAREA(4,ICFL)/dt
      ENDDO
      ENDDO
!      if(ieul2ph>0) then
!       xta=0.d0
!      endif
!
      end subroutine metric_speed

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE piston(ical_mvmsh,
     &        iter,time,
     &        NVRTX,MXVRTX,MXFACE,MXCELL,
     &        movflg,cord,lbface,lacell,lvface,lvcell,lcface,
     &        NBOUND,NBFS,SFBOUN,NFBOUN,IFFACE,icall)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,        only : ifle,ifll,gdScale
      use module_movegrid,  only : ratio,rps,LINER,the0,pisBC,mv_dir,RRR
      
!
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      integer,intent(in)    :: ical_mvmsh,iter,icall
      real*8 ,intent(in)    :: time
      integer,intent(in)    :: NVRTX,MXVRTX,MXFACE,MXCELL
      integer,intent(inout) :: movflg(MXVRTX)
      real*8 ,intent(inout) :: cord(3,MXVRTX)
      integer,intent(in)    :: NBOUND,NBFS
      integer,intent(in)    :: NFBOUN(NBOUND)
      integer,intent(in)    :: IFFACE(4,NBFS,NBOUND)
      integer,intent(in)    :: lbface(2,MXFACE)
      integer,intent(in)    :: lvcell(8,MXCELL)
      integer,intent(in)    :: lvface(4,MXFACE)    
      integer,intent(in)    :: lcface(2,MXFACE)
      integer,intent(in)    :: lacell(  MXCELL)
      character*80,intent(in) :: SFBOUN(0:NBOUND)
!_____________________________________________________________________
!
!     ical_mvmsh :
!               =1  Running [Moving] solver 
!               =2  Checking CV-GEOM of [Moving] mesh no-running solver
!               =3  Viewing [Moving] & [Removing/Adding] mesh
!               =4  Checking CV-GEOM of [Moving] & [Removing/Adding] mesh 
!               =5  Running [Moving] & [Removing/Adding] mesh solver
!               =0  No [Moving] & [Removing/Adding] Mesh solver
!_____________________________________________________________________
!     
!     movflg: Grid type or layer-number
!             movflg(IV)=0 : no moving vertex
!             movflg(IV)>0 : layer-number (moving grid)
!             movflg(IV)<0 : removing/add 
!_____________________________________________________________________
!
! --- [LOCAL ENTITIES]
!
      integer :: mb,iv,I,J,layer,layera,iloc,nb,ic1,ic2,iv3

      integer :: i_umat1=1,i_umat2=2
!_____________________________________________________________________
!
! --- [user array]
!
      real*8 ,parameter :: SML=1.d-10
!      real*8 ,parameter :: Ztopini=0.359d-3
      real*8 ,parameter :: pi=3.1415926d0
      integer :: iter_ctl,ierr=0,IB,NBC
      integer,save,allocatable :: BCIV(:)
      real*8,save,allocatable :: z_ini(:),Ztopini(:)
      real*8            :: Ztop,the
      real*8            :: dum1,dum2,dum(3),xyz(3)
!
!---------------------------------------------------------------------
! --- (Step-1: icall=1) First call for layer defination 
!---------------------------------------------------------------------
!
      if(icall==1) then
!
        allocate(z_ini(MXVRTX),Ztopini(MXVRTX),stat=ierr)
        if(ierr/=0) then
           write(*,*) ' ERR: allocate in z_ini'
           call FFRABORT(1,'STOP AT: subroutine piston')
        endif
!
        movflg(:)=0
        z_ini(:)=0.d0
        Ztopini(:)=0.d0
!
!
        do mb=1,NBOUND
        if(mb==pisBC) then
          do IB=1,NFBOUN(mb)
          do i=1,4
          iv=IFFACE(i,IB,mb)
          if(iv==0) cycle
          z_ini(iv)=1.d0
          enddo
          enddo
        endif
        enddo
!
        NBC=0
        do iv=1,NVRTX
          if(z_ini(iv)>0.5d0) then
            NBC=NBC+1
          endif
        enddo
!
        if(NBC/=0) then
          allocate(BCIV(NBC),stat=ierr)
          if(ierr/=0) call FFRABORT(1,'ERR: allocate error BCIV')
        else
          call FFRABORT(1,'ERR: NBC=0')
        endif
!
        NBC=0
        do iv=1,NVRTX
          if(z_ini(iv)>0.5d0) then
            NBC=NBC+1
            BCIV(NBC)=iv
          endif
        enddo
!
        do iv=1,NVRTX
        do i=1,NBC
        ib=BCIV(i)
        dum(:)=cord(1:3,ib)
        dum1=(dum(1)-cord(1,iv))**2*(1.d0-dble(mv_dir))/(1.d-20+1.d0-dble(mv_dir))
     &      +(dum(2)-cord(2,iv))**2*(2.d0-dble(mv_dir))/(1.d-20+2.d0-dble(mv_dir))
     &      +(dum(3)-cord(3,iv))**2*(3.d0-dble(mv_dir))/(1.d-20+3.d0-dble(mv_dir))
        if(dum1<1.d-15) then
          Ztopini(iv)=cord(mv_dir,ib)
        endif
        enddo
        enddo
!
        do iv=1,NVRTX
        z_ini(iv)=cord(mv_dir,iv)
        enddo

      endif
!
!-----------------------
! --- moving mesh
!-----------------------
      if(icall==2) then 
! 
       the=the0+rps*2.d0*pi*time 
       dum1=RRR*((1.d0-cos(the))-1.d0/ratio
     &     *(1.d0-sqrt(1.d0-ratio**2*sin(the)**2))) 
       Ztop=2.d0*RRR-dum1
       do iv=1,NVRTX
       dum2=LINER-(LINER-Ztop)*(LINER-z_ini(iv))/(LINER-Ztopini(iv))
       cord(mv_dir,iv)=dum2
       enddo
!
       write(ifll,'(1X,a,2E15.4)') 'the=,piston-top=',the*360/2/pi,Ztop

      endif
!_____________________________________________________________________
      return
      end SUBROUTINE piston
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
