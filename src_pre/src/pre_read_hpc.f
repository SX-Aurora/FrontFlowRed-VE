!
!      subroutine read_hpc_grid
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_hpc_grid
     &(iraddo,icpu,mcell,mvrtx,mface,nvrtx,ncell,NBCpair,
     & cord,lacell,lvcell,LEFACE,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_io      ,only : ifle,ifll
      use module_boundary,only : nbcnd,boundName,numbname,ivbcnd,
     &                           lvbcnd,icbinr,lcbinr,
     &                           kdbcnd,kdprdc,kdsld,kdintr,
     &                           prdcAxis,prdcOffset,ivpair,
     &                           boundIDMap,SFBOUN,NBOUND,idis,kdbuff,
     &                           set_rotmtrx1,kdshutr,kdpors
      use module_param ,only   : scalfctr,smallNrm
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)      :: icpu,mcell,mvrtx,nvrtx,ncell,mface,
     &                           iraddo !jiang
      integer,intent(out)     :: NBCpair
      integer,intent(out)     :: lacell(  mcell)
      integer,intent(out)     :: lvcell(8,mcell)
      real*8 ,intent(out)     :: cord  (3,mvrtx)
      integer,intent(out)     :: LEFACE(5,mface)
      integer,intent(inout)   :: ierror
!
! --- [local entities]
!
      integer :: nn0,nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8
      integer :: NNX,i,j,k,MAXBCV,nbs
      integer :: ic,iv
      integer :: ierr1,ihpc,iface,ii,idum(4)
      real(8) :: dum1,dum2,dum3,dumy(3), rot(3,3)
      integer             :: IBIDa,IBIDb,icL,ivv,kd
      integer             :: tempInt,M,Int,nb,mb,IS
      character*80 :: tmpStr80
!
!
      deallocate(ivbcnd,lvbcnd)
!-------------------
! --- Read Grid file
!-------------------
      open (11,file=GRIDout(icpu), status= 'unknown')
      read (11,*) HEADER1
      read (11,*) HEADER2
      read (11,*) nn3, NNX
      read (11,'(5e14.6)') ((cord(k,i),k=1,nn3),i=1,NNX)
      cord=DBLE(cord)
! ---  Local CONNECTIVITY
      read (11,*) HEADER2
      read (11,*) nn8, NNX
      read (11,'(6i12)') ((lvcell(k,i),k=1,nn8),i=1,NNX)
! --- Material Type
      read (11,*) HEADER2
      read (11,*) nn1,NNX
      read (11,'(6i12)') (lacell(icL),icL=1,NNX)
      read (11,*) HEADER2
!
      close (11)
!-------------------
! --- Read BC file
!-------------------
      open (11,file=BCout(icpu),status='unknown')
      read (11,*) HEADER1
      read (11,*) HEADER2
      read (11,*) NBOUND
!
      allocate(SFBOUN(NBOUND))
      allocate(boundIDMap(nbcnd,numbname))
      allocate( givbcnd(0:NBOUND), glvbcnd(IVBCTOT_HPC(icpu)+1), !onishi
     &          stat=ierr1 )
      allocate(ivbcnd(0:nbcnd), lvbcnd(IVBCTOT_HPC(icpu)))
      givbcnd=0
      glvbcnd=0
      ivbcnd=0
      lvbcnd=0
      ivpair(:)=0
!
      NBCpair=0
      do 100 mb=1,NBOUND
      iface=0;nbs=0
      read (11,'(a80)') tmpStr80
      SFBOUN(mb) = trim(adjustl(tmpStr80))
      read (11,'(3i12)') nn1,ivv,nbs
      givbcnd(mb)=givbcnd(mb-1)+ivv
      read (11,'(6i12)') (glvbcnd(i),i=givbcnd(mb-1)+1,givbcnd(mb))
      if(nbs>0) then
        do i=1,nbs
        read(11,'(4i12)') idum(1:4)
        NBCpair=NBCpair+1
        do j=1,4
        LEFACE(j,NBCpair)=idum(j)
        enddo
        LEFACE(5,NBCpair)=mb
        enddo
      endif
 100  continue
      close(11)
!
      boundIDMap=-1
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j).ne.' ') then
        do mb=1,NBOUND
        if(adjustl(boundName(nb,j)).eq.adjustl(SFBOUN(mb))) then
          boundIDMap(nb,j)=mb
          exit
        endif
	end do
      end if
      end do
      end do
!------------------------
! --- Reording glvbcnd
!------------------------
      do 200 nb=1,nbcnd
      if(kdbcnd(0,nb).eq.kdprdc.and.idis(nb)==0) then
!-------------------
! --- periodic BC: -
!-------------------
	IBIDa=boundIDMap(nb,1)
 	IBIDb=boundIDMap(nb,2)
	allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	pMap=-1
	if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
	  do 210 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	  do 220 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
          if(((cord(2,glvbcnd(k))-cord(2,glvbcnd(j))
     &        -prdcOffset(nb,2))**2
     &       +(cord(3,glvbcnd(k))-cord(3,glvbcnd(j))
     &        -prdcOffset(nb,3))**2
     &        ).lt.smallNrm*scalfctr(nb)) then
	    pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &      = glvbcnd(k)
	    exit
	  end if
 220      continue
 210      continue
	else if((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
	  do 230 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	  do 240 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	  if(((cord(1,glvbcnd(k))-cord(1,glvbcnd(j))- 
     &         prdcOffset(nb,1))**2+
     &        (cord(3,glvbcnd(k))-cord(3,glvbcnd(j))- 
     &         prdcOffset(nb,3))**2)<smallNrm*scalfctr(nb)) then
	    pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &      = glvbcnd(k)
  	    exit
          end if
 240      continue
 230      continue
        else if((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
          do 250 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
          do 260 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	  if(((cord(2,glvbcnd(k))-cord(2,glvbcnd(j))- 
     &         prdcOffset(nb,2))**2+ 
     &        (cord(1,glvbcnd(k))-cord(1,glvbcnd(j))- 
     &         prdcOffset(nb,1))**2)<smallNrm*scalfctr(nb)) then
	    pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &      = glvbcnd(k)
	    exit
	  end if
 260      continue
 250      continue
        elseif((prdcAxis(nb)=='r').or.(prdcAxis(nb)=='R')) then
            rot(:,:)=0.d0
            kd=kdbcnd(0,nb)
            call set_rotmtrx1(nbcnd,kd,nb,rot)
            do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)

            dum1=cord(1,glvbcnd(j))
            dum2=cord(2,glvbcnd(j))
            dum3=cord(3,glvbcnd(j))
            dumy(1)=rot(1,1)*dum1!cord(1,glvbcnd(j))
     &             +rot(2,1)*dum2!cord(2,glvbcnd(j))
     &             +rot(3,1)*dum3!cord(3,glvbcnd(j))
            dumy(2)=rot(1,2)*dum1!cord(1,glvbcnd(j))
     &             +rot(2,2)*dum2!cord(2,glvbcnd(j))
     &             +rot(3,2)*dum3!cord(3,glvbcnd(j))
            dumy(3)=rot(1,3)*dum1!cord(1,glvbcnd(j))
     &             +rot(2,3)*dum2!cord(2,glvbcnd(j))
     &             +rot(3,3)*dum3!cord(3,glvbcnd(j))
            do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
            dum1=(dumy(1)-cord(1,glvbcnd(k)))**2
     &          +(dumy(2)-cord(2,glvbcnd(k)))**2
     &          +(dumy(3)-cord(3,glvbcnd(k)))**2
            IF(dum1<smallNrm*scalfctr(nb)) then
              pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &        = glvbcnd(k)
              exit
            endif
            enddo
            enddo
	else
	  WRITE(ifle,*)
     &    '(read_hpc_grid)Periodic axis setting is wrong.',
     &    prdcAxis(nb)
	  stop
	end if
!
        do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
          ii=pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
          if(ii==-1) then
            stop 'peridic BC error-1'
          endif
        enddo
!
        do 280 j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
          do 290 k=j+1,givbcnd(IBIDb)
	  if(pMap(j)==pMap(k)) then
	    WRITE(ifle,'(1x,2a)')  
     &      'MSG: ',
     &      'Periodic boundary mapping error (duplicated vertex).'
	    stop
	  end if
 290      continue
 280    continue
!
!
        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!
	deallocate(pMap)
      elseif(kdbcnd(0,nb)==kdintr.and.idis(nb)==0) then
!------------------------
! --- interface BC mesh:
!------------------------
       	IBIDa=boundIDMap(nb,1)
 	IBIDb=boundIDMap(nb,2)
        allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	pMap(:)=-1
!
        if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &     (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	  WRITE(ifle,*) '(read_hpc)',
     &    'Number of Interface boundary vertex is not matched.',
     &    (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &    (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	  ierror=1
	  stop 'read_hpc'
        else
          WRITE(ifll,'(a,I8,a,I4)') 
     &   'MSG: Number of interface boundary NODE pair number =',
     &    givbcnd(IBIDa)-givbcnd(IBIDa-1),' at BC No.=',nb
	endif
!
        do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
        dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &       (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &       (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
!
        if(dum1.lt.smallNrm*scalfctr(nb)) then
          pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))=glvbcnd(k)
  	  exit
        endif
        enddo
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
	  WRITE(ifle,*)
     &    '(read_hpc_grid)',
     &    'Interface boundary mapping error (no vertex).'
	  ierror=1
	  stop 'read_hpc'
	end if
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
        do k=j+1,givbcnd(IBIDb)
	if(pMap(j)==pMap(k)) then
	  WRITE(ifle,*)  
     &     '(read_hpc_grid)',
     &     ' Interface boundary mapping error',
     &     ' (duplicated vertex on identical BC).'
	  ierror=1
	  stop 'read_hpc'
	endif
        enddo
        enddo
!
        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))

        deallocate(pMap)
      elseif(kdbcnd(0,nb)==kdbuff.and.idis(nb)==0) then
!------------------------
! --- Buffle BC mesh:
!------------------------
       	IBIDa=boundIDMap(nb,1)
 	IBIDb=boundIDMap(nb,2)
        allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	pMap(:)=-1
!
        if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &     (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	  WRITE(ifle,*) '(read_hpc)',
     &    'Number of Buffle boundary vertex is not matched.',
     &    (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &    (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	  ierror=1
	  stop 'read_hpc'
        else
          WRITE(ifll,'(a,I8,a,I4)') 
     &   'MSG: Number of Buffle boundary NODE pair number =',
     &    givbcnd(IBIDa)-givbcnd(IBIDa-1),' at BC No.=',nb
	endif
!
        do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
        dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &       (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &       (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
!
        if(dum1.lt.smallNrm*scalfctr(nb)) then
          pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))=glvbcnd(k)
  	  exit
        endif
        enddo
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
	  WRITE(ifle,*)
     &    '(read_hpc_grid)',
     &    'Buffle boundary mapping error (no vertex).'
	  ierror=1
	  stop 'read_hpc'
	end if
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
        do k=j+1,givbcnd(IBIDb)
	if(pMap(j)==pMap(k)) then
	  WRITE(ifle,*)  
     &     '(read_hpc_grid)',
     &     ' Buffle boundary mapping error',
     &     ' (duplicated vertex on identical BC).'
	  ierror=1
	  stop 'read_hpc'
	endif
        enddo
        enddo
!
        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
        deallocate(pMap)
      elseif(kdbcnd(0,nb)==kdshutr.and.idis(nb)==0) then
!------------------------
! --- Shutter BC mesh:
!------------------------
       	IBIDa=boundIDMap(nb,1)
 	IBIDb=boundIDMap(nb,2)
        allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	pMap(:)=-1
!
        if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &     (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	  WRITE(ifle,*) '(read_hpc)',
     &    'Number of Buffle boundary vertex is not matched.',
     &    (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &    (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	  ierror=1
	  stop 'read_hpc'
        else
          WRITE(ifll,'(a,I8,a,I4)') 
     &   'MSG: Number of Buffle boundary NODE pair number =',
     &    givbcnd(IBIDa)-givbcnd(IBIDa-1),' at BC No.=',nb
	endif
!
        do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
        dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &       (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &       (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
!
        if(dum1.lt.smallNrm*scalfctr(nb)) then
          pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))=glvbcnd(k)
  	  exit
        endif
        enddo
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
	  WRITE(ifle,*)
     &    '(read_hpc_grid)',
     &    'Buffle boundary mapping error (no vertex).'
	  ierror=1
	  stop 'read_hpc'
	end if
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
        do k=j+1,givbcnd(IBIDb)
	if(pMap(j)==pMap(k)) then
	  WRITE(ifle,*)  
     &     '(read_hpc_grid)',
     &     ' Buffle boundary mapping error',
     &     ' (duplicated vertex on identical BC).'
	  ierror=1
	  stop 'read_hpc'
	endif
        enddo
        enddo
!
        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
        deallocate(pMap)



      elseif(kdbcnd(0,nb)==kdpors.and.idis(nb)==0) then
!------------------------
! --- Porous BC mesh:
!------------------------
       	IBIDa=boundIDMap(nb,1)
 	IBIDb=boundIDMap(nb,2)
        allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	pMap(:)=-1
!
        if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &     (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	  WRITE(ifle,*) '(read_hpc)',
     &    'Number of Buffle boundary vertex is not matched.',
     &    (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &    (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	  ierror=1
	  stop 'read_hpc'
        else
          WRITE(ifll,'(a,I8,a,I4)') 
     &   'MSG: Number of Buffle boundary NODE pair number =',
     &    givbcnd(IBIDa)-givbcnd(IBIDa-1),' at BC No.=',nb
	endif
!
        do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
        dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &       (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &       (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
!
        if(dum1.lt.smallNrm*scalfctr(nb)) then
          pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))=glvbcnd(k)
  	  exit
        endif
        enddo
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
	  WRITE(ifle,*)
     &    '(read_hpc_grid)',
     &    'Buffle boundary mapping error (no vertex).'
	  ierror=1
	  stop 'read_hpc'
	end if
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
        do k=j+1,givbcnd(IBIDb)
	if(pMap(j)==pMap(k)) then
	  WRITE(ifle,*)  
     &     '(read_hpc_grid)',
     &     ' Buffle boundary mapping error',
     &     ' (duplicated vertex on identical BC).'
	  ierror=1
	  stop 'read_hpc'
	endif
        enddo
        enddo
!
        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
        deallocate(pMap)



      elseif(kdbcnd(0,nb)==kdsld.and.idis(nb)==0) then
!-----------------------
! --- sliding BC mesh:
!-----------------------
       	IBIDa=boundIDMap(nb,1)
 	IBIDb=boundIDMap(nb,2)
        allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	pMap(:)=-1
!
        if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &     (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	  WRITE(ifle,*) '(read_hpc)',
     &    'Number of sliding boundary vertex is not matched.',
     &    (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &    (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	  ierror=1
	  return
	endif
!
        do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
        dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &       (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &       (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
!
        if(dum1.lt.smallNrm*scalfctr(nb)) then
          pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))=glvbcnd(k)
  	  exit
        endif
        enddo
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
	  WRITE(ifle,*)
     &    '(read_hpc_grid)',
     &    'Sliding boundary mapping error (no vertex).'
	  ierror=1
	  return
	end if
        enddo
!
        do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
        do k=j+1,givbcnd(IBIDb)
	if(pMap(j)==pMap(k)) then
	  WRITE(ifle,*)  
     &     '(read_hpc_grid)',
     &     ' Sliding boundary mapping error',
     &     ' (duplicated vertex on identical BC).'
	  ierror=1
	  return
	endif
        enddo
        enddo
!
        glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &  pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))

        deallocate(pMap)
      end if
 200  continue
!
!--------------------------------------------
! --- 
!--------------------------------------------
      do 400 nb=1,nbcnd
      IBIDa=boundIDMap(nb,1)
      IBIDb=boundIDMap(nb,2)
      ivbcnd(nb)=ivbcnd(nb-1)
      ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDa)-givbcnd(IBIDa-1)
      if((ivbcnd(nb)-ivbcnd(nb-1))/=(givbcnd(IBIDa) 
     &  -givbcnd(IBIDa-1))) then
	WRITE(ifle,*)  
     &  '(read_hpc_grid)',
     &  'Number of periodic boundary vertex is not matched.'
        stop
      end if
      lvbcnd(ivbcnd(nb-1)+1:ivbcnd(nb))=
     &  glvbcnd(givbcnd(IBIDa-1)+1:givbcnd(IBIDa))
!
      if((kdbcnd(0,nb)==kdprdc.and.idis(nb)==0).or.
     &   (kdbcnd(0,nb)==kdsld.and.idis(nb)==0) .or.
     &   (kdbcnd(0,nb)==kdbuff.and.idis(nb)==0) .or.
     &   (kdbcnd(0,nb)==kdpors.and.idis(nb)==0) .or.
     &   (kdbcnd(0,nb)==kdshutr.and.idis(nb)==0) .or.
     &   (kdbcnd(0,nb)==kdintr.and.idis(nb)==0)) then
	tempInt=ivbcnd(nb)
	ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDb)-givbcnd(IBIDb-1)
	lvbcnd(tempInt+1:ivbcnd(nb))= 
     &  glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
        ivpair(nb)=tempInt
      elseif((kdbcnd(0,nb)==kdprdc.and.idis(nb)==1).or.
     &       (kdbcnd(0,nb)==kdsld .and.idis(nb)==1).or.
     &       (kdbcnd(0,nb)==kdbuff.and.idis(nb)==1).or.
     &       (kdbcnd(0,nb)==kdpors.and.idis(nb)==1).or.
     &       (kdbcnd(0,nb)==kdshutr.and.idis(nb)==1).or.
     &       (kdbcnd(0,nb)==kdintr.and.idis(nb)==1)) then
        tempInt=ivbcnd(nb)
	ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDb)-givbcnd(IBIDb-1)
	lvbcnd(tempInt+1:ivbcnd(nb))=
     &  glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
        ivpair(nb)=tempInt
      end if
 400  continue
!
      deallocate(givbcnd,glvbcnd)
!      deallocate(SFBOUN,boundIDMap)
      if(iraddo.EQ.0) deallocate(SFBOUN,boundIDMap)
      return
!
      end subroutine read_hpc_grid

