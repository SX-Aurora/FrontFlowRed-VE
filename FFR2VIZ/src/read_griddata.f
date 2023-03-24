!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_griddata(gridData,gdformat,ffgFormatKey,rsformat,
     &                         ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use FFRreaddata,only : cntlnam
      use FFRdata, only    : nvrtx,ncell,NBOUND,SFBOUN,NBFS,NFCE,NFBOUN,
     &                       IBFACE,IFFACE,lacell,lvcell,cord,kmesh,
     &                       ical_mvmsh,icon_cvg
      use SCRYUpre, only   : NDE
!
      implicit none

! 1. Input boundary file

! --- [dummy arguments]
      character*80,intent(in) :: gridData
      character*4,intent(in)  :: gdformat
      character*4,intent(in)  :: ffgFormatKey
      character*2,intent(in)  :: rsformat
      integer,intent(out)     :: ierror

! --- [local entities]
      character*84 :: fnam
      integer :: i,j,k,n,ios
      CHARACTER*60        :: FV_NAME,strID
      integer             :: major_version,minor_version,ifvreg
!
!
      ierror=0
      if(gdformat=='FV') then
!
!-[FieldView]-----------------------------------------------------------
!
! --- Read FieldView grid file created with Gridgen
!
      fnam=trim(adjustl(gridData))
      OPEN(20,FILE=fnam,FORM='FORMATTED',action='read',iostat=ios)
      if(ios/=0) then
        write(*,*)'Cannot open file :',trim(fnam),'.'
        ierror=1
        stop
      end if

!     check FV/UNS version
      READ(20,'(A60)') FV_NAME
      IF(FV_NAME(1:9) .ne. 'FIELDVIEW') THEN
	WRITE(*,'(1x,a)') 
     &    'MSG: This is NOT a FIELDVIEW ASCII file'
        ierror=1
	return
      END IF
      rewind(20)
      read(FV_NAME,*) strID,major_version,minor_version
      if(major_version==3) then
         fnam=trim(fnam) // '.fvreg'
         ifvreg=21
         OPEN(ifvreg,FILE=fnam,FORM='FORMATTED',
     &               action='read',iostat=ios)
         if(ios/=0) then
           write(*,'(1x,3a)')
     &     'WRG:  Region file cannot open:',trim(fnam)
           ifvreg=-1             
         end if

         call readFVGrids_30(20,ifvreg,ios)
         if(ifvreg/=-1) then
            close(ifvreg)
         end if
      elseif(major_version==2 .and. minor_version==4) then
         call readFVGrids(20,ios)
      end if
!
      close(20)
      if(ios/=0) then
        ierror=1
        stop 'Reading FV Grid'
      endif
!
!#ifdef GAMBIT
      else if(gdformat=='GB') then
!
!-[Gambit]-------------------------------------------------------
!
!
! ---  read GAMBIT(FLUENT) grid file (*.cas or *.msh)
      fnam=trim(adjustl(gridData))
      open(20,file=fnam,form='formatted',action='read',iostat=ios)
      if (ios/=0) then
        write(*,*) 'Cannot open file:',trim(fnam),'.'
        ierror=1
        stop
      end if
!
      call readGAMBITgrids(20,ios)
      if (ios/=0) then
        ierror=1
        return
      end if
      close(20)
!
      else if(gdformat=='FF') then
c
c  ---- read FrontFlowRed grid file (*.frontflow)
c
      fnam=trim(adjustl(gridData))
      open(20,file=fnam,form='unformatted',action='read',iostat=ios)
      if (ios/=0) then
        write(*,*) 'Cannot open file:',trim(fnam),'.'
        ierror=1
        stop
      end if
c
      call readFFRold(20,ios)
      if (ios/=0) then
        ierror=1
        return
      end if
      close(20)
c
      else if(gdformat=='GF') then
!
!-[General format]-----------------------------------------------------
!
! --- read GF grid file greated by Gridgen
!
        fnam=trim(adjustl(gridData))

        if(ffgFormatKey=='A') then
c warning: caracter size is differ 84/=80
          call readFFRGrid_A(20,fnam(1:80),ios)
!          call readFFRGrid_A
!     &      (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell,ios)
        else if(ffgFormatKey=='U') then
           call readFFRGrid(20,fnam(1:80),ios)
!           call readFFRGrid
!     &         (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell,ios)
        else
          ierror=1;return
        end if
        if(ios/=0) then
          ierror=1;return
        endif
!
! -[]----
!
      else
         write(*,*) 'Unknown -gf format.'
         stop
      end if

!     check reading data
      if(nvrtx<1 .or. ncell<1) then
         write(*,*) 'Error: no vrtx/cell found in read_griddata.'
         ierror=1
         return
      end if
!

!
      return

!//////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine readFVGrids_30(IUTFV,IUTRG,ierrorFVG)
!=================================================
!     Get following parameter
!     :  nvrtx,ncell,NCVIN,NCV
!
!     This subroutine read 'FieldView' grid format 
!     which is created with 'Gridgen'
!=================================================
!
      implicit none
!
! --- [module arguments]
!
      integer,intent(in)     :: IUTFV,IUTRG
      integer,intent(inout)  :: ierrorFVG
!
! --- [local entities]
!
      integer        :: NGRIDS,FNVAR,NME,NG,J,I,elemType
      CHARACTER*1024 :: tempStrings
      CHARACTER*60   :: FV_NAME,strID1,strID2
      real(4)        :: tempReal,tempReal2,tempReal3
      integer        :: tempInt,boundID,boundNodes
      integer      :: nvrtx_tmp,ncell_tmp,readStat
      integer      :: iNV
      integer      :: iFS,NBFS_tmp=0
      integer      :: iNE
      INTEGER,allocatable :: IBFACE_tmp(:)
      INTEGER,allocatable :: IFFACE_tmp(:,:)
!
      integer, parameter  :: maxdmn=100
      integer             :: idum,idum1
      character*10, allocatable :: dmnnm(:)
      integer, allocatable      :: dmn1(:),dmn2(:),dmn0(:)
! --- 
!
      iNV=0
      iFS=0
      iNE=0
      nvrtx_tmp=0
      ncell_tmp=0
      allocate(dmn1(maxdmn),dmn2(maxdmn),dmn0(maxdmn))
      allocate(dmnnm(maxdmn))
      dmn0=0
      dmn1=0
      dmn2=0

      WRITE(*,*) ' ####    Reading FV File...'
!
      READ(IUTFV,'(A60)') FV_NAME
      WRITE(*,*) ' ####    Reading field: ',trim(FV_NAME)
      IF(FV_NAME(1:9) .ne. 'FIELDVIEW') THEN
      WRITE(*,*) 
     &    ' ####    ',
     &    'This is NOT a FIELDVIEW ASCII file'
      ierrorFVG=1
      stop 'readFVGrids'
      END IF
!
      READ(IUTFV,2000) tempStrings ! comment
      READ(IUTFV,2000) tempStrings ! comment
      READ(IUTFV,2000) tempStrings
      read(tempStrings,*) strID1, NGRIDS
      WRITE(*,*) ' ####    ',
     &              'Reading field: ',trim(strID1)
      IF (trim(strID1).ne.'Grids') THEN
        WRITE(*,*) ' ####    ',
     & 'Section header exception happened at GRIDS'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      IF (NGRIDS.LT.1) THEN
        WRITE(*,*) ' ####    ',
     &                'Multi grid file is not supported'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      READ(IUTFV,2000) tempStrings
      read(tempStrings,*) strID1, strID2, NBOUND
      WRITE(*,*) ' ####    ',
     &              'Reading field: ',trim(strID1),' ',trim(strID2)
      IF(trim(strID1).ne.'Boundary'.or.trim(strID2).ne.'Table') THEN
        WRITE(*,*) ' ####    ', 
     &  'Section header exception happened at Boundary Table'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
 2010 FORMAT(I2,I2,I2,A80)
      allocate(SFBOUN(1:NBOUND))
      allocate(NFBOUN(NBOUND))
      NFBOUN(:)=0
      DO J = 1 , NBOUND
      READ(IUTFV,2010) tempInt,tempInt,tempInt, SFBOUN(J)
      WRITE(*,*) ' ####    ',
     &         'Reading boundary: ',trim(SFBOUN(j))
      end do
!
!     count node/element number
      nvrtx=0
      ncell=0
      do while(.true.)
        read(IUTFV,2000,iostat=readStat) tempStrings
        IF (readStat/=0) exit
        tempStrings=trim(tempStrings) // ' 0'
        read(tempStrings,*) strID1, tempInt
        IF (trim(strID1).eq.'Nodes') then
          nvrtx=nvrtx+tempInt
        ELSEIF (trim(strID1) .eq.'Elements') then
          do while(.true.)
            read(IUTFV,2000,iostat=readStat) tempStrings
            IF (readStat/=0) exit
            tempStrings=trim(tempStrings) // ' 0'
            read(tempStrings,*) strID1, tempInt
            IF (trim(strID1).eq.'Nodes') then
              backspace(IUTFV)
              exit
            end if
            ncell=ncell+1
          end do
        end if
      end do
      rewind(IUTFV)
      read(IUTFV,2000) tempStrings ! FV_NAME
      read(IUTFV,2000) tempStrings ! comment
      read(IUTFV,2000) tempStrings ! comment
      do while(.true.)
        read(IUTFV,2000) tempStrings
        tempStrings=trim(tempStrings) // ' 0'
        read(tempStrings,*) strID1, tempInt
        IF (trim(strID1).eq.'Nodes') then
          backspace(IUTFV)
          exit
        end if
      end do
      write(*,*)' #### Total NODE #: ',nvrtx
      write(*,*)' #### Detect ',ncell,'elements'
      allocate(cord(1:3,1:nvrtx))
      allocate(lacell(    1:ncell))
      allocate(lvcell(1:8,1:ncell))
      allocate(kmesh(     1:ncell))

!
! --- loop with NGRIDS
!
      do NG=1,NGRIDS
!
! --- 'Nodes'
!
      READ(IUTFV,2000) tempStrings
      read(tempStrings,*) strID1, tempInt
      WRITE(*,*) ' ####    ',
     &               'Reading field: ',trim(strID1)
      IF (trim(strID1).ne.'Nodes') THEN
	WRITE(*,*) ' ####    ',
     &                ' Section header exception happened at Nodes'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      iNV=nvrtx_tmp
      nvrtx_tmp=nvrtx_tmp+tempInt
!
      DO j = iNV+1, nvrtx_tmp
      READ(IUTFV,*) tempReal,tempReal2,tempReal3
      cord(1,j)=DBLE(tempReal)
      cord(2,j)=DBLE(tempReal2)
      cord(3,j)=DBLE(tempReal3)
      end do
!
! --- 'Boundary Faces'
!
      READ(IUTFV,2000) tempStrings
      read(tempStrings,*) strID1, strID2, NBFS_tmp
      WRITE(*,*) ' ####    ',
     &              'Reading field: ',trim(strID1),' ',trim(strID2)
      IF (trim(strID1).ne.'Boundary' .or. trim(strID2).ne.'Faces') THEN
	WRITE(*,*) ' ####    ',
     &           'Section header exception happened at Boundary Faces'
	ierrorFVG=1
	stop'readFVGrids'
      ENDIF
!
      if(NG==1) then
         iFS=0
         NBFS=NBFS_tmp
         allocate(IBFACE(NBFS))
         allocate(IFFACE(4,NBFS))
         IBFACE=0
         IFFACE=-1
      else
         iFS=NBFS
         allocate(IBFACE_tmp(NBFS))
         allocate(IFFACE_tmp(4,NBFS))
         IBFACE_tmp(1:NBFS)=IBFACE(1:NBFS)
         IFFACE_tmp(1:4,1:NBFS)=IFFACE(1:4,1:NBFS)
         NBFS=NBFS+NBFS_tmp
         deallocate(IBFACE)
         deallocate(IFFACE)
         allocate(IBFACE(NBFS))
         allocate(IFFACE(4,NBFS))
         IBFACE=0
         IFFACE=-1
         IBFACE(1:iFS)=IBFACE_tmp(1:iFS)
         IFFACE(1:4,1:iFS)=IFFACE_tmp(1:4,1:iFS)
         deallocate(IBFACE_tmp)
         deallocate(IFFACE_tmp)
      end if
      DO j=iFS+1,NBFS
	READ(IUTFV,2000)  tempStrings
	READ(tempStrings, *) boundID,boundNodes
	NFBOUN(boundID) = NFBOUN(boundID) + 1
        IBFACE(j)=boundID
	READ(tempStrings,*) boundID,boundNodes,
     & (IFFACE(i,j),i=1,boundNodes)
!       shift vertex No.
        do i=1,4
        if(IFFACE(i,j)>0) then
           IFFACE(i,j)=
     &     IFFACE(i,j)+iNV
        end if
        end do
      end do
!
      do j=iFS+1,NBFS
      do i=1,4
        if(IFFACE(i,j)==0) then
          write(*,*)'Node list contains zero.'
          stop'readFVGrids'
        end if
      end do
      end do
!
      do j=iFS+1,NBFS
      do i=1,4
        if(IFFACE(i,j)==-1) IFFACE(i,j)=0
      end do
      end do
!
! --- 'Elements'
!
      READ(IUTFV,2000) tempStrings
      WRITE(*,*) ' ####    ',
     &               'Reading field: ',trim(tempStrings)
      IF (trim(tempStrings).ne.'Elements') THEN
	WRITE(*,*) ' ####    ',
     &  'Section header exception happened at Elements'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      iNE=ncell_tmp
      dmn0(NG)=NG
      dmn1(NG)=iNE+1
!
      do
        READ(IUTFV,2000,iostat=readStat) tempStrings
        IF (readStat/=0) then
            ncell_tmp=iNE
            dmn2(NG)=iNE
            exit
        end if
        tempStrings=trim(tempStrings) // ' 0'
        read(tempStrings,*) strID1, tempInt
        IF (trim(strID1).eq.'Nodes') then
          ncell_tmp=iNE
          dmn2(NG)=iNE
          backspace(IUTFV)
          exit
        endif
	iNE=iNE+1

        READ(tempStrings, *)elemType
        select case(elemType)
        case(1)   ! --- Tet mesh
          kmesh(iNE) = 1
          READ(tempStrings, *) i, tempInt,
     &     lvcell(4,iNE), lvcell(1,iNE), lvcell(2,iNE), lvcell(3,iNE)
!          lacell(iNE)=NG
!          shift vertex No.
           lvcell(1:4,iNE)=lvcell(1:4,iNE)+iNV
        case(2)   ! --- Hex mesh
          kmesh(iNE) = 2
          READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(2,iNE), lvcell(4,iNE), lvcell(3,iNE),
     &     lvcell(5,iNE), lvcell(6,iNE), lvcell(8,iNE), lvcell(7,iNE)
!          lacell(iNE)=NG
!          shift vertex No.
           lvcell(1:8,iNE)=lvcell(1:8,iNE)+iNV
        case(3)   ! --- Prism mesh
          kmesh(iNE) = 3
          READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(4,iNE), lvcell(5,iNE), lvcell(2,iNE),
     &     lvcell(6,iNE), lvcell(3,iNE)
!          lacell(iNE)=NG
!          shift vertex No.
           lvcell(1:6,iNE)=lvcell(1:6,iNE)+iNV
        case(4)   ! --- Pyramid mesh
          kmesh(iNE) = 4
          READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(2,iNE), lvcell(3,iNE), lvcell(4,iNE),
     &     lvcell(5,iNE)
!          lacell(iNE)=NG
!          shift vertex No.
           lvcell(1:5,iNE)=lvcell(1:5,iNE)+iNV
        case default
          ierrorFVG=1
          stop 'Incompatible element detected'
        end select
      end do
!
! --- end loop with NGRIDS
!
      end do !do ng=1,NGRIDS

!
! --- define material number
!
      if(IUTRG/=-1) then
      READ(IUTRG,2000) tempStrings ! comment
      READ(IUTRG,2000) tempStrings ! comment
      READ(IUTRG,'(A60)') FV_NAME
      WRITE(*,'(1x,2a)') 'MSG: Reading field: ',trim(FV_NAME)
      IF(FV_NAME(1:5) .ne. 'FVREG') THEN
	WRITE(*,'(1x,a)') 
     &    'MSG: This is NOT a FIELDVIEW REGION file'
        dmnnm='fluid'
!	ierrorFVG=1
!	return
      ELSE
        READ(IUTRG,2000) tempStrings ! DATASET_COORD_TYPE
        do ng=1,NGRIDS
           READ(IUTRG,2000) tempStrings
           IF (trim(tempStrings).ne.'REGION') THEN
              WRITE(*,'(1x,2a)') 'MSG: ',
     &         ' Section header exception happened at REGION'
              dmnnm='fluid'
              exit
!              ierrorFVG=1
!              return
           END IF
           READ(IUTRG,2000) tempStrings ! NUM_GRIDS 1
           READ(tempStrings,'(A10)') dmnnm(ng)
           READ(IUTRG,2000) tempStrings ! NUM_GRIDS 1
           READ(IUTRG,2000) tempStrings ! (IMAT==NG)
        end do
      END IF
      else
        dmnnm='fluid'
      end if ! if(IUTRG/=-1) then

      do i=1,maxdmn
        if(dmn0(i)/=0) then
          do k=dmn1(i),dmn2(i)
            lacell(k)=dmn0(i)
          end do
        write(*,*) 'IMAT_U(=lacell)=',dmn0(i),' Region Name:'
     &      ,trim(dmnnm(i)),' cell# from',dmn1(i),' to',dmn2(i)
        endif
      enddo

      WRITE(*,*) ' #### !!!FINISH READING!!!'
2000  FORMAT(A1024)
      ierrorFVG=0
!
      end subroutine readFVGrids_30

!=================================================
      subroutine readFVGrids(IUTFV,ierrorFVG)
!=================================================
!     Get following parameter
!     :  nvrtx,ncell,NCVIN,NCV
!
!     This subroutine read 'FieldView' grid format 
!     which is created with 'Gridgen'
!=================================================
!
      implicit none
!
! --- [module arguments]
!
      integer,intent(in)     :: IUTFV
      integer,intent(inout)  :: ierrorFVG
!
! --- [local entities]
!
      integer        :: NGRIDS,FNVAR,NME,J,I,elemType
      CHARACTER*1024 :: tempStrings
      CHARACTER*60   :: FV_NAME
      real(4)        :: tempReal,tempReal2,tempReal3
      integer        :: tempInt,boundID,boundNodes
!
! --- 
!
      WRITE(*,*) ' ####    Reading FV File...'
!
      READ(IUTFV,'(A60)') FV_NAME
      WRITE(*,*) ' ####    Reading field: ',trim(FV_NAME)
      IF(FV_NAME(1:9) .ne. 'FIELDVIEW') THEN
      WRITE(*,*) 
     &    ' ####    ',
     &    'This is NOT a FIELDVIEW ASCII file'
      ierrorFVG=1
      stop 'readFVGrids'
      END IF
!
      READ(IUTFV,2000) tempStrings
      WRITE(*,*) 
     &   ' ####    Reading field: ',trim(tempStrings)
      IF (trim(tempStrings).ne. 'CONSTANTS') THEN
        WRITE(*,*) ' ####    ', 
     &       'Section header exception happened at CONSTANTS'
      ierrorFVG=1
      stop 'readFVGrids'
      END IF
      READ(IUTFV,*) tempReal	!	TIME
      READ(IUTFV,*) tempReal	!	FSMACH
      READ(IUTFV,*) tempReal	!	ALPHA
      READ(IUTFV,*) tempReal	!	RE


      READ(IUTFV,2000) tempStrings
      WRITE(*,*) ' ####    ',
     &              'Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'GRIDS') THEN
        WRITE(*,*) ' ####    ',
     & 'Section header exception happened at GRIDS'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      READ(IUTFV,*) NGRIDS
      IF (NGRIDS.GT.1) THEN
        WRITE(*,*) ' ####    ',
     &                'Multi grid file is not supported'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      READ(IUTFV,2000) tempStrings
      WRITE(*,*) ' ####    ',
     &              'Reading field: ',trim(tempStrings)
      IF(tempStrings.ne.'Boundary Table') THEN
        WRITE(*,*) ' ####    ', 
     &  'Section header exception happened at Boundary Table'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      READ(IUTFV,*) NBOUND
 2010 FORMAT(I1,A80)
      allocate(SFBOUN(1:NBOUND))
      DO J = 1 , NBOUND
      READ(IUTFV,2010) tempInt, SFBOUN(J)
      WRITE(*,*) ' ####    ',
     &         'Reading boundary: ',tempInt,trim(SFBOUN(j))
      end do
!
      READ(IUTFV,2000) tempStrings
      WRITE(*,*) ' ####    ',
     &               'Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'Variable Names') THEN
	WRITE(*,*)  ' ####    ', 
     &  'Section header exception happened at Variable Names'
	ierrorFVG=1
	stop'readFVGrids'
      END IF
!
      READ(IUTFV,*) FNVAR
      IF (FNVAR.GE.1)
     &  WRITE(*,*) ' ####    ',
     &                'Sorry, this version does not support variables'
      DO j = 1 ,FNVAR	
	READ(IUTFV,2000) tempStrings	
	WRITE(*,*) ' ####    ',
     &                'WARRNING! IGNORE :', trim(tempStrings)
      end do
!
! --- 'Nodes'
!
      READ(IUTFV,2000) tempStrings
      WRITE(*,*) ' ####    ',
     &               'Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'Nodes') THEN
	WRITE(*,*) ' ####    ',
     &                ' Section header exception happened at Nodes'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
      READ(IUTFV,*) nvrtx
!
      allocate(cord(1:3,1:nvrtx))
      DO j = 1, nvrtx
      READ(IUTFV,*) tempReal,tempReal2,tempReal3
      cord(1,j)=DBLE(tempReal)
      cord(2,j)=DBLE(tempReal2)
      cord(3,j)=DBLE(tempReal3)
      end do
!
! --- 'Boundary Faces'
!
      READ(IUTFV,2000) tempStrings
      WRITE(*,*) ' ####    ',
     &              'Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'Boundary Faces') THEN
	WRITE(*,*) ' ####    ',
     &           'Section header exception happened at Boundary Faces'
	ierrorFVG=1
	stop'readFVGrids'
      ENDIF
!
      READ(IUTFV,*) NBFS
      allocate(NFBOUN(NBOUND))
      allocate(IBFACE(NBFS))
      allocate(IFFACE(4,NBFS))
      NFBOUN=0
      IBFACE=0
      IFFACE=-1
      DO j=1,NBFS
	READ(IUTFV,2000)  tempStrings
	READ(tempStrings, *) boundID,boundNodes
	NFBOUN(boundID) = NFBOUN(boundID) + 1
        IBFACE(j)=boundID
	READ(tempStrings,*) boundID,boundNodes,
     & (IFFACE(i,j),i=1,boundNodes)
      end do
!
      do j=1,NBFS
      do i=1,4
        if(IFFACE(i,j)==0) then
          write(*,*)'Node list contains zero.'
          stop'readFVGrids'
        end if
      end do
      end do
!
      do j=1,NBFS
      do i=1,4
        if(IFFACE(i,j)==-1) IFFACE(i,j)=0
      end do
      end do
!
!
! --- 'Elements'
!
      READ(IUTFV,2000) tempStrings
      WRITE(*,*) ' ####    ',
     &               'Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'Elements') THEN
	WRITE(*,*) ' ####    ',
     &  'Section header exception happened at Elements'
	ierrorFVG=1
	stop 'readFVGrids'
      END IF
!
!     count element number
      NME=0
      do while(.true.)
        read(IUTFV,2000) tempStrings
        IF (tempStrings .eq.'Variables') exit
        NME=NME+1
      end do
      rewind(IUTFV)
      do while(.true.)
        read(IUTFV,2000) tempStrings
        IF (tempStrings .eq.'Elements') exit
      end do
      write(*,*)' #### Detect ',NME,'elements'
      ncell=NME
      
      allocate(lacell(1:ncell)); lacell=-1
      allocate(lvcell(1:8,1:ncell));   lvcell=-1 !***
      allocate(kmesh(1:ncell))
      do j=1,ncell
        READ(IUTFV,2000) tempStrings
        READ(tempStrings, *)elemType
        select case(elemType)
        case(1)   ! --- Tet mesh
          kmesh(j) = 1
          READ(tempStrings, *) i, tempInt,
     &     lvcell(4,j), lvcell(1,j), lvcell(2,j), lvcell(3,j)
          lacell(j)=1
        case(2)   ! --- Hex mesh
          kmesh(j) = 2
          READ(tempStrings, *) i, tempInt,
     &     lvcell(1,j), lvcell(2,j), lvcell(4,j), lvcell(3,j),
     &     lvcell(5,j), lvcell(6,j), lvcell(8,j), lvcell(7,j)
          lacell(j)=1
        case(3)   ! --- Prism mesh
          kmesh(j) = 3
          READ(tempStrings, *) i, tempInt,
     &     lvcell(1,j), lvcell(4,j), lvcell(5,j), lvcell(2,j),
     &     lvcell(6,j), lvcell(3,j)
          lacell(j)=1
        case(4)   ! --- Pyramid mesh
          kmesh(j) = 4
          READ(tempStrings, *) i, tempInt,
     &     lvcell(1,j), lvcell(2,j), lvcell(3,j), lvcell(4,j),
     &     lvcell(5,j)
          lacell(j)=1
        case default
          ierrorFVG=1
          stop 'Incompatible element detected'
        end select
      end do
      read(IUTFV,2000) tempStrings
      if (tempStrings/='Variables') then
        ierrorFVG=1
        stop'readFVGrids'
      end if
!
      write(*,*) ' #### Total NODE #: ',NDE
      WRITE(*,*) ' #### !!!FINISH READING!!!'
2000  FORMAT(A1024)
      ierrorFVG=0
!
      end subroutine readFVGrids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine readGAMBITgrids(IUTFV,ierrorGB)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
! --- [dummy arguments]
!
      integer, intent(in)       :: IUTFV
      integer, intent(inout)    :: ierrorGB
!
! --- [local entities]
!
      integer, allocatable      :: mfboun(:),nfpz(:)
      character*80, allocatable :: nmboun(:)
      character*10, allocatable :: dmnnm(:)
      integer, allocatable      :: dmn1(:),dmn2(:),dmn0(:)
! *** only in FFR2VIZ
!      integer, allocatable      :: kmesh(:)
!
      integer                   :: i,j,ii,jj,boundID,izn,izone,ic
      integer                   :: clen,isum,nn1,ier
      character*60              :: FV_NAME
      character*80              :: tempStrings,nextStrings
      character*1               :: tmpc7,tmpc8,elem_ty
      character*2               :: tmpc1,tmpc2
      character*40              :: tmpc3,tmpc4,face_id,fnum_id,zone_id
      character*9 :: tmpc5,tmpc6,tmpc9,tmpc10,tmpc11,tmpc12,tmpc13
     &               ,firs_id,last_id
! DEBUG
      integer  :: IMAT,NFCE,IS1,IS2,IS3
!      integer  :: NBFS,IMAT,NFCE,IS1,IS2,IS3
      integer  :: IV1,IV2,IV3,IV4,IV5,IV6,IV,IIS,IERR
      integer  :: pnface,psface,nf,nc,nbnd,bndsfCnt
      integer  :: tmpI1,tmpI2,tmpI3,tmpI4,tmpI5,tmpI6,tmpI7
      integer  :: nnf4(2),nnf6(6),nnf(6),IVN0(4)
      integer  :: ICNUM(2),IVNUM(4),IVN1(4),IVN2(4),IVN3(4),IVN4(4)
      real(4)  :: tempReal,tempReal2,tempReal3,aa
      integer :: IC1,IC2,IS
! *** only in FFR2VIZ
      integer, allocatable :: lvface(:,:),lcface(:,:),lfcell(:,:)
      integer  :: ifll,ifle
!
      integer   :: strptr
      character :: strID
      logical   :: gridgenflag
      integer, parameter :: maxboun=100, maxdmn=100
!      integer, parameter :: maxboun=40, maxdmn=30
      integer   :: nvrtx_start,nvrtx_end,nvrtxx
      nvrtx_start=0
      nvrtx_end=0
!
      nbnd=0
! *** only in FFR2VIZ
      ifll=6
      ifle=6
!      lacell=0
!      lvcell=0
!
! --- ----------------------------------------------------
      write(ifll,*) 'Reading GAMBIT grid file ...'
! --- ----------------------------------------------------
      boundID=0
      allocate(nmboun(1:maxboun)) ; nmboun=' '
      allocate(mfboun(maxboun))   ; mfboun=0
      allocate(nfpz(maxboun))     ; nfpz=0
!
      allocate(dmn1(maxdmn),dmn2(maxdmn),dmn0(maxdmn))
      allocate(dmnnm(maxdmn))
!
      dmn0=0
      dmn1=0
      dmn2=0
      dmnnm=' '
      IMAT=0
!
!     ###############################################
!     ##### read data regarding BC information ######
!     ###############################################
!
      read(IUTFV,'(A60)') FV_NAME
      write(ifll,*) 'First LINE in the FILE :',trim(FV_NAME)
      gridgenflag=.false.
      clen = len_trim(adjustl(FV_NAME))
!     search string 'Gridgen' in grid file.
      do i=1,clen
         strptr=scan(FV_NAME(i:clen),'G')
         if(strptr>=1.and.strptr<=clen-7) then
            if(FV_NAME(strptr:strptr+7)=='Gridgen') then
               gridgenflag=.true.
            end if
         end if
      end do
!
! --- 
!
! *** only in FFR2VIZ
!      lvface=0
!      lcface=0
!      lfcell=0
!
! --- -------------------------
!
      nf=0
      nvrtxx=0
      do
        read(IUTFV,'(A60)',iostat=ier) tempStrings
        if(ier/=0) exit
        read(IUTFV,'(A60)',iostat=ier) nextStrings
        if(ier/=0) nextStrings=''
        nextStrings=trim(adjustl(nextStrings))
        if(nextStrings(1:1)=='(' .and. len_trim(nextStrings)==1) then
!          *** skip the line which begin with' ('
        else
           backspace(IUTFV)
        end if
!
        if(tempStrings(1:1)== '(') then
          read(tempStrings(2:3),*,iostat=ier) tmpI1
          if(ier/=0) cycle
          if(tmpI1==13) then                      ! FACES
            write(*,*) 
            write(*,*) 'start read 13-zone'
            read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc10
            if(len_trim(tmpc5)>3) then
              face_id=tmpc5(5:5)
              fnum_id=tmpc9
            else
              face_id=tmpc6(2:3)
              fnum_id=tmpc10
            endif
            if(trim(face_id)=='0') then
              read(fnum_id,'(z8)') tmpI1
              NFCE=tmpI1
! *** only in FFR2VIZ
              allocate(lvface(1:4,1:NFCE))
              allocate(lcface(1:2,1:NFCE))
              lvface=0
              lcface=0
              write(ifll,*) 'FACE number: NFCE =',NFCE
            else
              if(len_trim(tmpc5)>3) then
                face_id=tmpc5(5:5)
                fnum_id=tmpc9
                read(tempStrings,*) 
     &                      tmpc5,tmpc6,tmpc9,tmpc10,tmpc11
                elem_ty=tmpc11(1:1)
                zone_id=tmpc5(5:6)
                firs_id=tmpc6
                last_id=tmpc9
              else
                face_id=tmpc6(2:3)
                fnum_id=tmpc10
                read(tempStrings,*) 
     &                       tmpc5,tmpc6,tmpc9,tmpc10,tmpc11,tmpc12
                elem_ty=tmpc12(1:1)
                zone_id=tmpc6(2:3)
                firs_id=tmpc9
                last_id=tmpc10
              endif
              read(zone_id,'(z8)') tmpI4   ! Zone-id
              read(firs_id,'(z8)') tmpI1   ! first-index
              read(last_id,'(z8)') tmpI2   ! last-index
              boundID=tmpI4                ! temporal zone-id
              psface=tmpI1
              pnface=tmpI2
!
! --- 0=mixed; 2=linear; 3=triangular; 4=quadrilateral -------
!
              write(*,*) 'start read lvface; lcface'
              write(*,*) 'face_ty= ',elem_ty
              write(*,*) 'face_ty= 0:mixed; 2:linear; 3:triangular; 4:quadrilateral'
              write(*,*) 'face no_start=; no_end=',psface,pnface
              if(elem_ty=='0') then
                do 100 i=psface,pnface
                nfpz(boundID)=nfpz(boundID)+1
                nf=nf+1
                read(IUTFV,'(A60)',iostat=ier) tempStrings
                read(tempStrings,*) tmpc5,tmpc6,tmpc9,
     &                              tmpc10,tmpc11,tmpc12
                read(tmpc5, '(z10)') tmpI1
                if(tmpI1.eq.3) then
                   read(tmpc6, '(z10)') tmpI2
                   read(tmpc9, '(z10)') tmpI3
                   read(tmpc10,'(z10)') tmpI4
                   read(tmpc11,'(z10)') tmpI5
                   read(tmpc12,'(z10)') tmpI6
                   mfboun(boundID)=mfboun(boundID)+1
                   lvface(1,nf)= tmpI2 ! node number
                   lvface(2,nf)= tmpI3
                   lvface(3,nf)= tmpI4
                   lcface(1,nf)= tmpI5 ! cell number
                   lcface(2,nf)= tmpI6
                elseif(tmpI1.eq.4) then
                   read(tempStrings,*) tmpc5,tmpc6,tmpc9,
     &                                 tmpc10,tmpc11,tmpc12,tmpc13
                   read(tmpc6, '(z10)') tmpI2
                   read(tmpc9, '(z10)') tmpI3
                   read(tmpc10,'(z10)') tmpI4
                   read(tmpc11,'(z10)') tmpI5
                   read(tmpc12,'(z10)') tmpI6
                   read(tmpc13,'(z10)') tmpI7
                   mfboun(boundID)=mfboun(boundID)+1
                   lvface(1,nf)= tmpI2 ! node number
                   lvface(2,nf)= tmpI3
                   lvface(3,nf)= tmpI4
                   lvface(4,nf)= tmpI5
                   lcface(1,nf)= tmpI6 ! cell number
                   lcface(2,nf)= tmpI7
                else
                   write(ifll,*) 
     &                   '## WARNING: face is not tri./quad. ##'
                end if
 100            enddo
              elseif(elem_ty=='2') then
                write(ifll,*) '<Read element-type: linear>'
                write(ifll,*) 'NOT YET SUPPORTED'
              elseif(elem_ty=='3') then
                do 110 i=psface,pnface
                nfpz(boundID)=nfpz(boundID)+1
                nf=nf+1
                read(IUTFV,'(A60)',iostat=ier) tempStrings
                read(tempStrings,*) tmpc6,tmpc9,tmpc10,tmpc11,tmpc12
                read(tmpc6, '(z10)') tmpI2
                read(tmpc9, '(z10)') tmpI3
                read(tmpc10,'(z10)') tmpI4
                read(tmpc11,'(z10)') tmpI5
                read(tmpc12,'(z10)') tmpI6
                mfboun(boundID)=mfboun(boundID)+1
                lvface(1,nf)=tmpI2    ! node number
                lvface(2,nf)=tmpI3
                lvface(3,nf)=tmpI4
                lcface(1,nf)=tmpI5    ! cell number
                lcface(2,nf)=tmpI6
 110            end do
                write(ifll,*) 
     &                '<Read element-type: triangular>: boundID=',
     &          boundID
              elseif(elem_ty=='4') then   ! modified01
                do 300 i=psface,pnface
                nf=nf+1
                nfpz(boundID)=nfpz(boundID)+1
                read(IUTFV,'(A60)') tempStrings
                read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc10,
     &                              tmpc11,tmpc12
                read(tmpc5, '(z8)') tmpI1
                read(tmpc6, '(z8)') tmpI2
                read(tmpc9, '(z8)') tmpI3
                read(tmpc10,'(z8)') tmpI4
                read(tmpc11,'(z8)') tmpI5
                read(tmpc12,'(z8)') tmpI6
                lvface(1,nf)=tmpI1    ! node number
                lvface(2,nf)=tmpI2
                lvface(3,nf)=tmpI3
                lvface(4,nf)=tmpI4
                lcface(1,nf)=tmpI5    ! cell number
                lcface(2,nf)=tmpI6
 300            enddo
              end if
              write(*,*) 'end   read lvface; lcface'
            end if
            write(*,*) 'end   read 13-zone'
          elseif(tmpI1==12) then      ! CELLS block
            write(*,*) 
            write(*,*) 'start read 12-zone'
            read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc10
            if(tmpc6(2:2)=='0') then
              read(tmpc10,'(z8)') tmpI1
              ncell=tmpI1
              write(ifll,*) 'CELL number: ncell =',tmpI1
              allocate(kmesh(1:tmpI1),stat=ier)
              if(ier/=0) stop 'stop at allocating in read_griddata(GB)'
              kmesh=0
            endif
!
            if(tmpc6(2:2)/='0') then
! --- 2 more string in 12) section, in "regular cell section"
              IMAT=IMAT+1
              read(tempStrings,*) 
     &                       tmpc5,tmpc6,tmpc9,tmpc10,tmpc11,tmpc12
              read(tmpc6(2:2),'(z8)') tmpI1      ! Zone-ID
              read(tmpc9,'(z8)')  tmpI2
              read(tmpc10,'(z8)') tmpI3
              clen=len_trim(adjustl(tmpc12))
              do i=1,clen
                 strptr=scan(tmpc12,'()')
                 if(strptr>=1.and.strptr<=clen) then
                    tmpc12=tmpc12(1:strptr-1)//tmpc12(strptr+1:clen)
                 end if
              enddo
              write(*,*) 'IMAT=',IMAT
              write(*,*) 'cell no_start=, no_end=',tmpI2,tmpI3
              read(tmpc12,'(z8)') tmpI4
              write(*,*) 'cell type= ',tmpI4
              write(*,*) 'cell type= 0:mixed; =2:tetra; =4:hexa; =5:pyramid; =6:prism'
! --- cell type
              if(tmpI4==2) then
!                tetra
                 kmesh(tmpI2:tmpI3)=2
              else if(tmpI4==4) then
!                hexa
                 kmesh(tmpI2:tmpI3)=4
              else if(tmpI4==5) then
!                pyramid
                 kmesh(tmpI2:tmpI3)=5
              else if(tmpI4==6) then
!                prism(wedge)
                 kmesh(tmpI2:tmpI3)=6
              else if(tmpI4==0) then
!                mixed
                 i=tmpI2
!ncell=13027644
!no_start=, no_end=1,13016964
                 do ! loop
                    if(i>tmpI3) exit
                    read(IUTFV,'(a)',ADVANCE='NO',EOR=101) tmpc7
 101                read(IUTFV,'(a)',ADVANCE='NO') tmpc7
                    do while(tmpc7==' ' .or. tmpc7=='(')
                       read(IUTFV,'(a)',ADVANCE='NO') tmpc7
!                       if(i>13016950) then
!                       endif
                    end do
                    if(tmpc7(1:1)==')') then
                       if(i-1/=tmpI3) then
                          write(*,*) 'Error: cell type num is not match'
                          ierrorGB=1
                          return
                       end if
                       exit
                    end if
                    read(tmpc7,'(z8)') tmpI7
                    if(tmpI7==2 .or. tmpI7==4 .or.
     &                 tmpI7==5 .or. tmpI7==6) then
                       kmesh(i)=tmpI7
                       i=i+1
                    else
                       write(*,*) 'Error: unknown type of cell(2)'
                       ierrorGB=1
                       return
                    end if
                 end do
              else
                 write(*,*) 'Error: unknown type of cell(1)'
                 ierrorGB=1
                 return
              end if
c! --- this subroutine allows only on hexahedron.
c              if(tmpI4.ne.4) then
c                 write(ifll,*) 
c     &                    '## WARNING: cells are not hexahedron ##'
c              end if
              dmn0(IMAT)=IMAT
              dmn1(IMAT)=tmpI2
              dmn2(IMAT)=tmpI3
            endif
            write(*,*) 'end   read 12-zone'
          elseif(tmpI1==10) then
            write(*,*) 
            write(*,*) 'start read 10-zone'
            read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc10
! *** only in FFR2VIZ
!            if(tmpc6(2:2)/='0') then
            if(tmpc6(2:2)=='0') then
              read(tmpc10,'(z8)') tmpI1
              write(*,*) 'NODE number: nvrtx =',tmpI1
              nvrtx=tmpI1
              allocate(cord(1:3,1:nvrtx))
            else
              read(tmpc9, '(z8)') tmpI1
              read(tmpc10,'(z8)') tmpI2
              write(*,*) 'start read cord'
              write(*,*) 'vertex no_start=; no_end=',tmpI1,tmpI2
              do i=tmpI1,tmpI2
              nvrtxx=nvrtxx+1
              read(IUTFV,*) tempReal,tempReal2,tempReal3
              cord(1,i)=dble(tempReal )
              cord(2,i)=dble(tempReal2)
              cord(3,i)=dble(tempReal3)
              enddo
              write(*,*) 'end   read cord'
            endif
            write(*,*)   'end   read 10-zone'
          elseif(tmpI1==39) then         !(maybe) boundary condition
            write(ifll,*) 'NOT YET SUPPORTED for Index 39'
          elseif(tmpI1==45) then         !boundary condition
            write(*,*) 
            write(*,*) 'start read 45-zone'
            read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc3
            read(tmpc6(2:3),*) tmpI1     ! Zone_ID
            if(trim(tmpc9)/='interior'.and.trim(tmpc9)/='fluid'.and.
     &         trim(tmpc9)/='solid') then
              nbnd=nbnd+1
              clen=len_trim(adjustl(tmpc3))
              do i=1,clen
                 strptr=scan(tmpc3,'()')
                 if(strptr>=1.and.strptr<=clen) then
                    tmpc3=tmpc3(1:strptr-1)//tmpc3(strptr+1:clen)
                 endif
              enddo
              if(gridgenflag) then
!                delete last ID(ex. '-3') for Gridgen .cas file
                 clen  = len_trim(adjustl(tmpc3))
                 strptr=0
                 do i=1,clen
                    strID=tmpc3(i:i)
                    if(strID=="-" .and. i+1<=clen) then
                       strptr=i
                    end if
                 end do
!                if strptr is near the end of 'tmpc3', almost 4 chars width
                 if(strptr/=0.and.strptr>=clen-4.and.strptr<=clen) then
!                   now, strptr is a char size of '-3'
                    strptr=int(log10(real(tmpI1)))+1
                    strptr=strptr+1
                    if(strptr>=1.and.strptr<=clen) then
                       tmpc3=tmpc3(1:clen-strptr)
                    end if
                 end if
              end if ! if(gridgenflag) then
              tmpc4=trim(tmpc3)
              nmboun(tmpI1)=trim(tmpc4)
            elseif(trim(tmpc9)=='fluid'.or.
     &             trim(tmpc9)=='solid') then
              clen=len_trim(adjustl(tmpc3))
              do i=1,clen
                strptr=scan(tmpc3,'()')
                if(strptr>=1.and.strptr<=clen) then
                  tmpc3=tmpc3(1:strptr-1)//tmpc3(strptr+1:clen)
                end if
              enddo
              if(gridgenflag) then
!               delete last ID(ex. '-3') for Gridgen .cas file
                clen  = len_trim(adjustl(tmpc3))
                strptr=0
                do i=1,clen
                   strID=tmpc3(i:i)
                   if(strID=="-" .and. i+1<=clen) then
                      strptr=i
                   end if
                end do
!               if strptr is near the end of 'tmpc3', almost 4 chars width
                if(strptr/=0.and.strptr>=clen-4.and.strptr<=clen) then
!                  now, strptr is a char size of '-3'
                   strptr=int(log10(real(tmpI1)))+1
                   strptr=strptr+1
                   if(strptr>=1.and.strptr<=clen) then
                      tmpc3=tmpc3(1:clen-strptr)
                   end if
                end if
              endif  ! if(gridgenflag) then
              tmpc4=trim(tmpc3)
              dmnnm(tmpI1)=trim(tmpc4)
            endif
            write(*,*) 'end   read 45-zone'
          endif
        endif
      enddo
! --- --------------------------------------------------------------
! *** only in FFR2VIZ
      if(nvrtx/=nvrtxx) then
         write(ifll,*) 
     &        'WRN: NVRTX should be :',nvrtxx,' in ',cntlnam
         nvrtx=nvrtxx
      endif
!!
!      if(nvrtx>mvrtx) then
!         write(ifll,*) 
!         write(ifle,*) 'ERR: mvrtx > ',nvrtx,' in ',cntlnam
!         stop
!      endif
!
!      if(nf.ne.NFCE) then
!        write(ifle,*) 'ERR: reading error in face number',nf,NFCE
!        stop
!      endif
!!
!      if(nf>mface) then
!        write(ifle,*) 'ERR: mface > ',nf,' in ',cntlnam
!        stop ' read_grid_data '
!      endif
!!
      write(ifll,*) 'Total face number=',nf
      write(ifll,*) 'BC number without INTERIOR and FLUID=',nbnd
!
      isum=0
      do izn=1,maxboun
        if(nmboun(izn)/=' ') then
          isum=isum+nfpz(izn)
        endif
      enddo
      bndsfCnt=isum
      NBFS=isum
      write(ifll,*) 'BC FACES number= ',bndsfCnt
!
! --- --------------------------------------------------------------
!
      allocate(SFBOUN(0:nbnd))
      allocate(NFBOUN(nbnd))
      allocate(IBFACE(bndsfCnt))
      allocate(IFFACE(4,bndsfCnt))
      NFBOUN=0
      IFFACE=0
      IBFACE=0
      SFBOUN=' '
      boundID=0
      do izn=1,maxboun
        if(nmboun(izn)/=' ') then
          boundID=boundID+1
          SFBOUN(boundID)=nmboun(izn)   ! name of BC
          NFBOUN(boundID)=nfpz(izn)
        endif
      enddo
!
! --- --------------------------------------------------------------
!
      write(ifll,*) 'BC NAME :'
      do i=1,nbnd
      write(ifll,*) 'BC FACE',i,trim(SFBOUN(i)),NFBOUN(i)
      end do
!
      deallocate(nfpz)
!
!    ##################################################
!    ###### read data again for node coordinates ######
!    ##################################################
!
      boundID=0
      nf=0
!
! --- Second
!
      rewind(IUTFV)
! --- --------------
      do
      read(IUTFV,'(A60)',iostat=ier) tempStrings
      if(ier/=0) exit
      if(tempStrings(1:1)== '(') then
        read(tempStrings(2:3),*,iostat=ier) tmpI1
        if(ier/=0) cycle
        if(tmpI1==13) then      ! FACES
          read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc10
          if(len_trim(tmpc5)>3) then
            face_id=tmpc5(5:5)
            fnum_id=tmpc9
          else
            face_id=tmpc6(2:3)
            fnum_id=tmpc10
          end if
          if(trim(face_id)=='0') then
          else
            if(len_trim(tmpc5)>3) then
              face_id= tmpc5(5:5)
              fnum_id= tmpc9
              read(tempStrings,*) 
     &                 tmpc5,tmpc6,tmpc9,tmpc10,tmpc11
              elem_ty= tmpc11(1:1)
              zone_id= tmpc5(5:6)
              firs_id= tmpc6
              last_id= tmpc9
            else
              face_id= tmpc6(2:3)
              fnum_id= tmpc10
              read(tempStrings,*) 
     &                  tmpc5,tmpc6,tmpc9,tmpc10,tmpc11,tmpc12
              elem_ty= tmpc12(1:1)
              zone_id= tmpc6(2:3)
              firs_id= tmpc9
              last_id= tmpc10
            end if
            read(zone_id,'(z8)') tmpI4      ! Zone-id
            read(firs_id,'(z8)') tmpI1      ! first-index
            read(last_id,'(z8)') tmpI2      ! last-index
            izone=tmpI4
            if(nmboun(tmpI4)/=' ') then
              boundID=boundID+1             !number of BC
              psface=tmpI1
              pnface=tmpI2
              if(elem_ty=='0') then
                write(ifll,*) '<Read element-type: mixed>'
                do i=psface,pnface
                nf=nf+1
                read(IUTFV,'(A60)',iostat=ier) tempStrings
                read(tempStrings,*) tmpc5,tmpc6,tmpc9,
     &                              tmpc10,tmpc11,tmpc12
                read(tmpc5, '(z10)') tmpI1
! 
                if(tmpI1.eq.3) then
                   read(tmpc6, '(z10)') tmpI2
                   read(tmpc9, '(z10)') tmpI3
                   read(tmpc10,'(z10)') tmpI4
                   read(tmpc11,'(z10)') tmpI5
                   read(tmpc12,'(z10)') tmpI6
                   IBFACE(nf)= boundID
                   IFFACE(1,nf)= tmpI2
                   IFFACE(2,nf)= tmpI3
                   IFFACE(3,nf)= tmpI4
                   IFFACE(4,nf)= 0
                elseif(tmpI1.eq.4) then
!                  if face has 4 nodes, one more char exists.
                   read(tempStrings,*) tmpc5,tmpc6,tmpc9,
     &                                 tmpc10,tmpc11,tmpc12,tmpc13
                   read(tmpc6, '(z10)') tmpI2
                   read(tmpc9, '(z10)') tmpI3
                   read(tmpc10,'(z10)') tmpI4
                   read(tmpc11,'(z10)') tmpI5
                   read(tmpc12,'(z10)') tmpI6
                   read(tmpc13,'(z10)') tmpI7
                   IBFACE(nf)= boundID
                   IFFACE(1,nf)= tmpI2
                   IFFACE(2,nf)= tmpI3
                   IFFACE(3,nf)= tmpI4
                   IFFACE(4,nf)= tmpI5
                end if
                end do
              elseif(elem_ty=='2') then
                write(ifll,*) '<Read element-type: linear>'
                write(ifll,*) 'NOT YET SUPPORTED'
              elseif(elem_ty=='3') then
                do i=psface,pnface
                nf=nf+1
                read(IUTFV,'(A60)',iostat=ier) tempStrings
                read(tempStrings,*) tmpc6,tmpc9,tmpc10,tmpc11,tmpc12
                read(tmpc6, '(z10)') tmpI2
                read(tmpc9, '(z10)') tmpI3
                read(tmpc10,'(z10)') tmpI4
                read(tmpc11,'(z10)') tmpI5
                read(tmpc12,'(z10)') tmpI6
                IBFACE(nf)= boundID
                IFFACE(1,nf)= tmpI2
                IFFACE(2,nf)= tmpI3
                IFFACE(3,nf)= tmpI4
                IFFACE(4,nf)= 0
                enddo
                write(ifll,*) 
     &              '<Read element-type: triangular> for IFFACE'
              elseif(elem_ty=='4') then
                do i=psface,pnface
                nf=nf+1
                read(IUTFV,'(A60)') tempStrings
                read(tempStrings,*) tmpc5,tmpc6,tmpc9,
     &                              tmpc10,tmpc11,tmpc12
                read(tmpc5, '(z8)') tmpI1
                read(tmpc6, '(z8)') tmpI2
                read(tmpc9, '(z8)') tmpI3
                read(tmpc10,'(z8)') tmpI4
                read(tmpc11,'(z8)') tmpI5
                read(tmpc12,'(z8)') tmpI6
                IBFACE(nf)= boundID
                IFFACE(1,nf)=tmpI1
                IFFACE(2,nf)=tmpI2
                IFFACE(3,nf)=tmpI3
                IFFACE(4,nf)=tmpI4
                enddo
              endif
            end if
          end if
        end if
      end if
      end do
!
! --- -------------------------
!
      NBOUND=boundID
      write(ifll,*) 'Face Number in BC:',nf,NFCE
      write(ifll,*) 'complete reading all the faces'
!
! *** only in FFR2VIZ
      allocate(lacell(1:ncell))
      allocate(lvcell(1:8,1:ncell))
      allocate(lfcell(1:7,1:ncell))
      lacell=0
      lvcell=0                 !????
      lfcell=0

! --- -------------------------
!     face -> cell construction for FFR order
!     IC1 < IC2 in Fluent
!
      do IS=1,NFCE
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      ICNUM(1)=IC1
      ICNUM(2)=IC2
      do II=1,2
      IC=ICNUM(II)
      if(IC.eq.0) cycle
!     Gambit IC2 face has a left-hand coordinates.
!      so, z-axis is inverted -> right-hand.
      if(II==2) then
! ---    reverse a face vrtx order
         if(lvface(4,IS).eq.0) then
            IVNUM(1)=lvface(1,IS)
            IVNUM(2)=lvface(3,IS)
            IVNUM(3)=lvface(2,IS)
         else
            IVNUM(1)=lvface(1,IS)
            IVNUM(2)=lvface(4,IS)
            IVNUM(3)=lvface(3,IS)
            IVNUM(4)=lvface(2,IS)
         end if
      else
         IVNUM(1:4)=lvface(1:4,IS)
      end if
!
      if(kmesh(IC)==2) then ! tetra ----------------------------
         if(lfcell(7,IC).eq.0) then
            lvcell(1,IC)=IVNUM(1)
            lvcell(2,IC)=IVNUM(2)
            lvcell(3,IC)=IVNUM(3)
            lfcell(1,IC)=IS
            lfcell(7,IC)=lfcell(7,IC)+1
         elseif(lfcell(7,IC).eq.1) then
            do IV=1,3
               if(lvcell(1,IC)/=IVNUM(IV).and.
     &            lvcell(2,IC)/=IVNUM(IV).and.
     &            lvcell(3,IC)/=IVNUM(IV) ) then
                  lvcell(4,IC)=IVNUM(IV)
                  lfcell(2,IC)=IS
                  lfcell(7,IC)=lfcell(7,IC)+1
                  exit ! important
               endif
            end do
         else
            lfcell(7,IC)=lfcell(7,IC)+1
            lfcell(lfcell(7,IC),IC)=IS
         endif

      elseif(kmesh(IC)==5) then ! pyramid -----------------------
         do IV=1,4
            if(lvface(IV,IS)/=0) then
               JJ=IV
            end if
         end do
         if(JJ==4) then ! base face(4 node)
            lvcell(1,IC)=IVNUM(1)
            lvcell(2,IC)=IVNUM(2)
            lvcell(3,IC)=IVNUM(3)
            lvcell(4,IC)=IVNUM(4)
         end if
         lfcell(7,IC)=lfcell(7,IC)+1
         lfcell(lfcell(7,IC),IC)=IS

      elseif(kmesh(IC)==6) then ! prism(wedge) ------------------
         do IV=1,4
            if(lvface(IV,IS)/=0) then
               JJ=IV
            end if
         end do
         if(JJ==3) then ! base face(4 node)
         if(lfcell(7,IC).eq.0) then
            lvcell(1,IC)=IVNUM(1)
            lvcell(2,IC)=IVNUM(2)
            lvcell(3,IC)=IVNUM(3)
            lfcell(1,IC)=IS
            lfcell(7,IC)=lfcell(7,IC)+1
         elseif(lfcell(7,IC).eq.1) then
            if(lvcell(1,IC)/=IVNUM(1).and.
     &         lvcell(1,IC)/=IVNUM(2).and.
     &         lvcell(1,IC)/=IVNUM(3).and.
     &         lvcell(2,IC)/=IVNUM(1).and.
     &         lvcell(2,IC)/=IVNUM(2).and.
     &         lvcell(2,IC)/=IVNUM(3).and.
     &         lvcell(3,IC)/=IVNUM(1).and.
     &         lvcell(3,IC)/=IVNUM(2).and.
     &         lvcell(3,IC)/=IVNUM(3)
     &        ) then
               lvcell(4,IC)=IVNUM(1)
               lvcell(6,IC)=IVNUM(2)
               lvcell(5,IC)=IVNUM(3)
               lfcell(2,IC)=IS
               lfcell(7,IC)=lfcell(7,IC)+1
            endif
         endif
         endif

      elseif(kmesh(IC)==4) then ! hexa --------------------------
         if(lfcell(7,IC).eq.0) then
            lvcell(1,IC)=IVNUM(1)
            lvcell(2,IC)=IVNUM(2)
            lvcell(3,IC)=IVNUM(3)
            lvcell(4,IC)=IVNUM(4)
            lfcell(1,IC)=IS
            lfcell(7,IC)=lfcell(7,IC)+1
         elseif(lfcell(7,IC).eq.1) then
            if(lvcell(1,IC)/=IVNUM(1).and.
     &         lvcell(1,IC)/=IVNUM(2).and.
     &         lvcell(1,IC)/=IVNUM(3).and.
     &         lvcell(1,IC)/=IVNUM(4).and.
     &         lvcell(3,IC)/=IVNUM(1).and.
     &         lvcell(3,IC)/=IVNUM(2).and.
     &         lvcell(3,IC)/=IVNUM(3).and.
     &         lvcell(3,IC)/=IVNUM(4)
     &        ) then
               lvcell(5,IC)=IVNUM(1)
               lvcell(8,IC)=IVNUM(2)
               lvcell(7,IC)=IVNUM(3)
               lvcell(6,IC)=IVNUM(4)
               lfcell(2,IC)=IS
               lfcell(7,IC)=lfcell(7,IC)+1
               lfcell(6,IC)=II ! temporary flag
            endif
         endif
      endif ! kmesh()
!
      enddo
      enddo
!
! --- Pyramid construction check
!
      do IC=1,ncell
      if(kmesh(IC)/=5) cycle
      do II=1,5
         IS=lfcell(II,IC)
         do IV=1,4
            if(lvface(IV,IS)/=0) then
               JJ=IV
            end if
         end do
         if(JJ==3) then ! side face(3 node)
            do IV=1,3
               if(lvcell(1,IC)/=lvface(IV,IS).and.
     &            lvcell(2,IC)/=lvface(IV,IS).and.
     &            lvcell(3,IC)/=lvface(IV,IS).and.
     &            lvcell(4,IC)/=lvface(IV,IS) ) then
                  lvcell(5,IC)=lvface(IV,IS)
                  exit
               endif
            end do
            exit ! important
         elseif(JJ/=4) then
            write(ifll,*) 'ERR: wrong pyramid face',IC,II,lvface(II,IC)
         end if
      end do
      end do
!
! --- Prism construction check
!
      do IS=1,NFCE
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      ICNUM(1)=IC1
      ICNUM(2)=IC2
      do 200 II=1,2
      IC=ICNUM(II)
      if(IC.eq.0) cycle
      if(kmesh(IC).ne.6) cycle ! important
      IS1=lfcell(1,IC)
      IS2=lfcell(2,IC)
      if(IS/=IS1.and.IS/=IS2) then
        lfcell(7,IC)=lfcell(7,IC)+1
        lfcell(lfcell(7,IC),IC)=IS
!
        if(lfcell(7,IC).eq.3) then
          IVNUM(1)=lvcell(4,IC)
          IVNUM(2)=lvcell(6,IC)
          IVNUM(3)=lvcell(5,IC)
          IV=lvface(1,IS)+lvface(2,IS)+lvface(3,IS)+lvface(4,IS)
          IVN0(1)=lvface(1,IS)
          IVN0(2)=lvface(2,IS)
          IVN0(3)=lvface(3,IS)
          IVN0(4)=lvface(4,IS)
!
! --- 1
          IV1=lvcell(1,IC)+lvcell(4,IC)+lvcell(5,IC)+lvcell(2,IC)
          IVN1(1)=lvcell(1,IC)
          IVN1(2)=lvcell(4,IC)
          IVN1(3)=lvcell(5,IC)
          IVN1(4)=lvcell(2,IC)
          IV2=lvcell(3,IC)+lvcell(6,IC)+lvcell(4,IC)+lvcell(1,IC)
          IVN2(1)=lvcell(3,IC)
          IVN2(2)=lvcell(6,IC)
          IVN2(3)=lvcell(4,IC)
          IVN2(4)=lvcell(1,IC)
          IV3=lvcell(2,IC)+lvcell(5,IC)+lvcell(6,IC)+lvcell(3,IC)
          IVN3(1)=lvcell(2,IC)
          IVN3(2)=lvcell(5,IC)
          IVN3(3)=lvcell(6,IC)
          IVN3(4)=lvcell(3,IC)
          IERR=0
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3)) then
            IERR=1
          else
            call list_fmatch4(IERR,IVN1,IVN0,0)
            if(IERR==0)  goto 200
            call list_fmatch4(IERR,IVN2,IVN0,0)
            if(IERR==0)  goto 200
            call list_fmatch4(IERR,IVN3,IVN0,0)
            if(IERR==0)  goto 200
          endif
          if(IERR.gt.0) then
             lvcell(6,IC)=IVNUM(1)
             lvcell(5,IC)=IVNUM(2)
             lvcell(4,IC)=IVNUM(3)
          ELSE
            goto 200
          endif
! --- 2
          IV1=lvcell(1,IC)+lvcell(4,IC)+lvcell(5,IC)+lvcell(2,IC)
          IVN1(1)=lvcell(1,IC)
          IVN1(2)=lvcell(4,IC)
          IVN1(3)=lvcell(5,IC)
          IVN1(4)=lvcell(2,IC)
          IV2=lvcell(3,IC)+lvcell(6,IC)+lvcell(4,IC)+lvcell(1,IC)
          IVN2(1)=lvcell(3,IC)
          IVN2(2)=lvcell(6,IC)
          IVN2(3)=lvcell(4,IC)
          IVN2(4)=lvcell(1,IC)
          IV3=lvcell(2,IC)+lvcell(5,IC)+lvcell(6,IC)+lvcell(3,IC)
          IVN3(1)=lvcell(2,IC)
          IVN3(2)=lvcell(5,IC)
          IVN3(3)=lvcell(6,IC)
          IVN3(4)=lvcell(3,IC)
          IERR=0
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3)) then
            IERR=2
          else
            call list_fmatch4(IERR,IVN1,IVN0,0)
            if(IERR==0)  goto 200
            call list_fmatch4(IERR,IVN2,IVN0,0)
            if(IERR==0)  goto 200
            call list_fmatch4(IERR,IVN3,IVN0,0)
            if(IERR==0)  goto 200
          endif
          if(IERR.gt.0) then
             lvcell(5,IC)=IVNUM(1)
             lvcell(4,IC)=IVNUM(2)
             lvcell(6,IC)=IVNUM(3)
          ELSE
            goto 200
          endif
!
! --- finial check
!
          IERR=0
          IV1=lvcell(1,IC)+lvcell(4,IC)+lvcell(5,IC)+lvcell(2,IC)
          IVN1(1)=lvcell(1,IC)
          IVN1(2)=lvcell(4,IC)
          IVN1(3)=lvcell(5,IC)
          IVN1(4)=lvcell(2,IC)
          IV2=lvcell(3,IC)+lvcell(6,IC)+lvcell(4,IC)+lvcell(1,IC)
          IVN2(1)=lvcell(3,IC)
          IVN2(2)=lvcell(6,IC)
          IVN2(3)=lvcell(4,IC)
          IVN2(4)=lvcell(1,IC)
          IV3=lvcell(2,IC)+lvcell(5,IC)+lvcell(6,IC)+lvcell(3,IC)
          IVN3(1)=lvcell(2,IC)
          IVN3(2)=lvcell(5,IC)
          IVN3(3)=lvcell(6,IC)
          IVN3(4)=lvcell(3,IC)
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3)) then
            IERR=4
          else
            call list_fmatch4(IERR,IVN1,IVN0,0)
            if(IERR==0)  goto 200
            call list_fmatch4(IERR,IVN2,IVN0,0)
            if(IERR==0)  goto 200
            call list_fmatch4(IERR,IVN3,IVN0,0)
            if(IERR==0)  goto 200
          endif
          if(IERR.gt.0) then
            write(ifll,*) 'IERR:',IERR
            write(*,1100) (lvcell(i,IC),i=1,6)
            write(*,1000) (lvface(i,IS),i=1,4)
            stop 'Contact your consultor: IC'
          else
            goto 200
          endif
        endif
      endif
 200  continue
      enddo
!
! --- Hexa construction check
!
      do IS=1,NFCE
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      ICNUM(1)=IC1
      ICNUM(2)=IC2
      do 210 II=1,2
      IC=ICNUM(II)
      if(IC.eq.0) cycle
      if(kmesh(IC).ne.4) cycle ! important
      IS1=lfcell(1,IC)
      IS2=lfcell(2,IC)
      if(IS/=IS1.and.IS/=IS2) then
        lfcell(7,IC)=lfcell(7,IC)+1
        lfcell(lfcell(7,IC),IC)=IS
!
        if(lfcell(7,IC).eq.3) then
          IVNUM(1)=lvcell(5,IC)
          IVNUM(2)=lvcell(8,IC)
          IVNUM(3)=lvcell(7,IC)
          IVNUM(4)=lvcell(6,IC)
          IV=lvface(1,IS)+lvface(2,IS)+lvface(3,IS)+lvface(4,IS)
          IVN0(1)=lvface(1,IS)
          IVN0(2)=lvface(2,IS)
          IVN0(3)=lvface(3,IS)
          IVN0(4)=lvface(4,IS)
!
! --- 1
          IV1=lvcell(1,IC)+lvcell(5,IC)+lvcell(6,IC)+lvcell(2,IC)
          IVN1(1)=lvcell(1,IC)
          IVN1(2)=lvcell(5,IC)
          IVN1(3)=lvcell(6,IC)
          IVN1(4)=lvcell(2,IC)
          IV2=lvcell(4,IC)+lvcell(8,IC)+lvcell(5,IC)+lvcell(1,IC)
          IVN2(1)=lvcell(4,IC)
          IVN2(2)=lvcell(8,IC)
          IVN2(3)=lvcell(5,IC)
          IVN2(4)=lvcell(1,IC)
          IV3=lvcell(3,IC)+lvcell(7,IC)+lvcell(8,IC)+lvcell(4,IC)
          IVN3(1)=lvcell(3,IC)
          IVN3(2)=lvcell(7,IC)
          IVN3(3)=lvcell(8,IC)
          IVN3(4)=lvcell(4,IC)
          IV4=lvcell(2,IC)+lvcell(6,IC)+lvcell(7,IC)+lvcell(3,IC)
          IVN4(1)=lvcell(2,IC)
          IVN4(2)=lvcell(6,IC)
          IVN4(3)=lvcell(7,IC)
          IVN4(4)=lvcell(3,IC)
          IERR=0
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3.or.IV==IV4)) then
            IERR=1
          else
            call list_fmatch4(IERR,IVN1,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN2,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN3,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN4,IVN0,0)
            if(IERR==0)  goto 210
          endif
          if(IERR.gt.0) then
             lvcell(8,IC)=IVNUM(1)
             lvcell(7,IC)=IVNUM(2)
             lvcell(6,IC)=IVNUM(3)
             lvcell(5,IC)=IVNUM(4)
          ELSE
            goto 210
          endif
! --- 2
          IV1=lvcell(1,IC)+lvcell(5,IC)+lvcell(6,IC)+lvcell(2,IC)
          IVN1(1)=lvcell(1,IC)
          IVN1(2)=lvcell(5,IC)
          IVN1(3)=lvcell(6,IC)
          IVN1(4)=lvcell(2,IC)
          IV2=lvcell(4,IC)+lvcell(8,IC)+lvcell(5,IC)+lvcell(1,IC)
          IVN2(1)=lvcell(4,IC)
          IVN2(2)=lvcell(8,IC)
          IVN2(3)=lvcell(5,IC)
          IVN2(4)=lvcell(1,IC)
          IV3=lvcell(3,IC)+lvcell(7,IC)+lvcell(8,IC)+lvcell(4,IC)
          IVN3(1)=lvcell(3,IC)
          IVN3(2)=lvcell(7,IC)
          IVN3(3)=lvcell(8,IC)
          IVN3(4)=lvcell(4,IC)
          IV4=lvcell(2,IC)+lvcell(6,IC)+lvcell(7,IC)+lvcell(3,IC)
          IVN4(1)=lvcell(2,IC)
          IVN4(2)=lvcell(6,IC)
          IVN4(3)=lvcell(7,IC)
          IVN4(4)=lvcell(3,IC)
          IERR=0
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3.or.IV==IV4)) then
            IERR=2
          else
            call list_fmatch4(IERR,IVN1,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN2,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN3,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN4,IVN0,0)
            if(IERR==0)  goto 210
          endif
          if(IERR.gt.0) then
             lvcell(7,IC)=IVNUM(1)
             lvcell(6,IC)=IVNUM(2)
             lvcell(5,IC)=IVNUM(3)
             lvcell(8,IC)=IVNUM(4)
          ELSE
            goto 210
          endif
!
! --- 3
!
          IV1=lvcell(1,IC)+lvcell(5,IC)+lvcell(6,IC)+lvcell(2,IC)
          IVN1(1)=lvcell(1,IC)
          IVN1(2)=lvcell(5,IC)
          IVN1(3)=lvcell(6,IC)
          IVN1(4)=lvcell(2,IC)
          IV2=lvcell(4,IC)+lvcell(8,IC)+lvcell(5,IC)+lvcell(1,IC)
          IVN2(1)=lvcell(4,IC)
          IVN2(2)=lvcell(8,IC)
          IVN2(3)=lvcell(5,IC)
          IVN2(4)=lvcell(1,IC)
          IV3=lvcell(3,IC)+lvcell(7,IC)+lvcell(8,IC)+lvcell(4,IC)
          IVN3(1)=lvcell(3,IC)
          IVN3(2)=lvcell(7,IC)
          IVN3(3)=lvcell(8,IC)
          IVN3(4)=lvcell(4,IC)
          IV4=lvcell(2,IC)+lvcell(6,IC)+lvcell(7,IC)+lvcell(3,IC)
          IVN4(1)=lvcell(2,IC)
          IVN4(2)=lvcell(6,IC)
          IVN4(3)=lvcell(7,IC)
          IVN4(4)=lvcell(3,IC)
          IERR=0
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3.or.IV==IV4)) then
            IERR=3
          else
            call list_fmatch4(IERR,IVN1,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN2,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN3,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN4,IVN0,0)
            if(IERR==0)  goto 210
          endif
          if(IERR.gt.0) then
             lvcell(6,IC)=IVNUM(1)
             lvcell(5,IC)=IVNUM(2)
             lvcell(8,IC)=IVNUM(3)
             lvcell(7,IC)=IVNUM(4)
          ELSE
            goto 210
          endif
!
! --- finial check
!
          IERR=0
          IV1=lvcell(1,IC)+lvcell(5,IC)+lvcell(6,IC)+lvcell(2,IC)
          IVN1(1)=lvcell(1,IC)
          IVN1(2)=lvcell(5,IC)
          IVN1(3)=lvcell(6,IC)
          IVN1(4)=lvcell(2,IC)
          IV2=lvcell(4,IC)+lvcell(8,IC)+lvcell(5,IC)+lvcell(1,IC)
          IVN2(1)=lvcell(4,IC)
          IVN2(2)=lvcell(8,IC)
          IVN2(3)=lvcell(5,IC)
          IVN2(4)=lvcell(1,IC)
          IV3=lvcell(3,IC)+lvcell(7,IC)+lvcell(8,IC)+lvcell(4,IC)
          IVN3(1)=lvcell(3,IC)
          IVN3(2)=lvcell(7,IC)
          IVN3(3)=lvcell(8,IC)
          IVN3(4)=lvcell(4,IC)
          IV4=lvcell(2,IC)+lvcell(6,IC)+lvcell(7,IC)+lvcell(3,IC)
          IVN4(1)=lvcell(2,IC)
          IVN4(2)=lvcell(6,IC)
          IVN4(3)=lvcell(7,IC)
          IVN4(4)=lvcell(3,IC)
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3.or.IV==IV4)) then
            IERR=4
          else
            call list_fmatch4(IERR,IVN1,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN2,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN3,IVN0,0)
            if(IERR==0)  goto 210
            call list_fmatch4(IERR,IVN4,IVN0,0)
            if(IERR==0)  goto 210
          endif
          if(IERR.gt.0) then
            write(ifll,*) 'IERR:',IERR
            write(*,1100) (lvcell(i,IC),i=1,8)
            write(*,1000) (lvface(i,IS),i=1,4)
            stop 'Contact your consultor: IC'
          else
            goto 210
          endif
        endif
      endif
 210  continue
      enddo
 1000 format(4x,'---',4I8)
 1100 format(4x,'===',8I8)
!
! --- ??? NOT necessary
!
      do IS=1,NFCE
      IC1=lcface(1,IS)
      IC2=lcface(2,IS)
      ICNUM(1)=IC1
      ICNUM(2)=IC2
      do II=1,2
      IC=ICNUM(II)
      if(IC.eq.0) cycle
      if(kmesh(IC).ne.4) cycle ! important
      IS1=lfcell(1,IC)
      IS2=lfcell(2,IC)
      IS3=lfcell(3,IC)
      if(IS/=IS1.and.IS/=IS2.and.IS/=IS3) then
        if(lfcell(7,IC).eq.3) then
          lfcell(4,IC)=IS
          lfcell(7,IC)=lfcell(7,IC)+1
          IV=lvface(1,IS)+lvface(2,IS)+lvface(3,IS)+lvface(4,IS)
          IV1=lvcell(1,IC)+lvcell(2,IC)+lvcell(5,IC)+lvcell(6,IC)
          IV2=lvcell(1,IC)+lvcell(4,IC)+lvcell(8,IC)+lvcell(5,IC)
          IV3=lvcell(3,IC)+lvcell(4,IC)+lvcell(7,IC)+lvcell(8,IC)
          IV4=lvcell(2,IC)+lvcell(3,IC)+lvcell(7,IC)+lvcell(6,IC)
          IERR=0
          if(.not.(IV==IV1.or.IV==IV2.or.IV==IV3.or.IV==IV4)) then
            IERR=1
            write(*,1000) (lvface(i,IS),i=1,4)
            write(*,1100) (lvcell(i,IC),i=1,8)
            STOP 'CHECK CELL'
          endif
        endif
      endif
      enddo
      enddo
!
! --- debug
!
      DO IC=1,ncell
        if(kmesh(IC).eq.2) then ! tetra
           II=4 ! face number
           JJ=4 ! vrtx number
        elseif(kmesh(IC).eq.4) then ! hexa
           II=6
           JJ=8
        elseif(kmesh(IC).eq.5) then ! pyramid
           II=5
           JJ=5
        elseif(kmesh(IC).eq.6) then ! prism
           II=5
           JJ=6
        else
           write(ifle,*) 'ERR: unknown cell type: IC, kmesh(IC): ',
     &          IC,kmesh(IC)
           stop
        end if
!
        if(lfcell(7,IC).ne.II) then
           write(ifle,*) ' ERR:',IC,' face num',lfcell(7,IC),' /=',II
           stop
        end if
        do iv=1,JJ
           if(lvcell(iv,IC).eq.0) then
              write(ifle,*) 'ERR: lvcell(iv,IC)=0: ',IC,iv
              stop 
           endif
        enddo
!
      enddo
      write(ifll,*) 'Completed connectity '
! +++++++++++++++++++++++++
!
      deallocate(nmboun)
      deallocate(mfboun)
!
!
! --- define material number
!
      do i=1,maxdmn
        if(dmn0(i)/=0) then
          do k=dmn1(i),dmn2(i)
            lacell(k)=dmn0(i)
          end do
        write(ifll,*) 'IMAT_U(=lacell)=',dmn0(i),' Region Name:'
     &      ,trim(dmnnm(i)),' cell# from',dmn1(i),' to',dmn2(i)
        endif
      enddo
!
! *** only in FFR2VIZ
!     convert kmesh value from FLUENT to FIELDVIEW
      do IC=1,ncell
        select case(kmesh(IC))
        case(2) ! tetra
           kmesh(IC)=1
        case(4) ! hexa
           kmesh(IC)=2
        case(5) ! pyramid
           kmesh(IC)=4
        case(6) ! prism
           kmesh(IC)=3
        case default
           write(ifle,*) 'ERR: unknown cell type: IC, kmesh(IC): ',
     &          IC,kmesh(IC)
           stop
        end select
      enddo

      deallocate(lfcell,lcface,lvface)
!      lfcell(:,:)=0
!      lcface(:,:)=0
!      lvface(:,:)=0
      deallocate(dmnnm)
!
      return
      end subroutine readGAMBITgrids

!#endif
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
!#ifdef SCRYU
!=================================================
!=================================================
      subroutine readFFRold(IUTFV,ierrorFVG)
!=================================================
!     Get following parameter
!     :  nvrtx,ncell,NCVIN,NCV
!
!     This subroutine read 'FrontFlowRed' grid format 
!     which is created with 'prefflow'
!=================================================
!
      implicit none
! --- [module arguments]
!
      integer,intent(in)     :: IUTFV
      integer,intent(inout)  :: ierrorFVG
!
! --- [local entities]
!
      integer           :: i,j,ILV,IS,IC,boundID
      CHARACTER*80      :: tempStrings
      integer           :: kclv,ical_mvmshx,mb
      integer,parameter :: nv2typ(8)=(/0,0,0,1,4,3,0,2/)
      integer,allocatable :: IBIDX(:),IFFACE_temp(:,:,:)
      
!
! --- 
!
      write(*,*) ' ####    Reading FFR Grid File...'
!
      read(IUTFV) nvrtx,ical_mvmshx,icon_cvg
      if(nvrtx<1) then
	write(*,*) ' ####    ',
     &             ' vertex data exception happened'
        write(*,*) '** may be wrong FFR Grid File'
	ierrorFVG=1
        stop 'ERR: readFFRGrids'
      end if
      allocate(cord(1:3,1:nvrtx))
      read(IUTFV) ((cord(i,j), i=1,3), j=1,nvrtx)

!     ---
      read(IUTFV) ncell
      if(ncell<1) then
	write(*,*) ' ####    ',
     &             ' cell data exception happened'
        write(*,*) '** may be wrong FFR Grid File'
	ierrorFVG=1
!	call abort(1,'readFFRGrids')
        stop 'ERR: readFFRGrids'
      end if
      allocate(lacell(1:ncell)); lacell=-1
      allocate(lvcell(1:8,1:ncell));   lvcell=-1 !***
      allocate(kmesh(1:ncell))
      read(IUTFV) (lacell(i),i=1,ncell)
      read(IUTFV) ((lvcell(i,j), i=1,8), j=1,ncell)
!----------------------
! --- set up cell type
!----------------------
      do IC=1,ncell
         kclv=0
         do ILV=1,8
            if(lvcell(ILV,IC).gt.0) kclv=ILV
         end do
         if( kclv.lt.1 ) goto 9001
         kmesh(IC)=nv2typ(kclv)
         if( (IC).lt.1 ) goto 9001
      end do
!     ---
      read(IUTFV) NBOUND
      read(IUTFV) NBFS

      if(NBOUND<1.or.NBFS<1) then
	write(*,*) ' ####    ',
     &             ' boundary data exception happened',
     &             ' NBOUND:',NBOUND,'NBFS:',NBFS
        write(*,*) '** may be wrong FFR Grid File'
	ierrorFVG=1
!	call abort(1,'readFFRGrids')
        stop 'ERR: readFFRGrids'
      end if

      allocate(SFBOUN(NBOUND))
      allocate(NFBOUN(NBOUND))
      allocate(IBFACE(NBFS))
      allocate(IFFACE(4,NBFS))
      allocate(IBIDX(0:NBOUND))
      IBIDX=0
      IS=0
      do i=1,NBOUND
      read(IUTFV) boundID,NFBOUN(boundID),SFBOUN(boundID)
      read(IUTFV) mb,IBIDX(boundID)
      read(IUTFV) IFFACE(1:4,IBIDX(boundID-1)+1:IBIDX(boundID))
      IBFACE(IBIDX(boundID-1)+1:IBIDX(boundID))=mb
      end do
      ierrorFVG=0
      return

 9001 continue
      write(*,*) ' ### error : data error'
      write(*,*) 'no. of vertices =',kclv
      write(*,*) 'it must be 4,5,6 or 8'
      write(*,*) 'cell no. =',IC
      goto 9999
 9999 continue
      write(*,*) '(readFFRGrids)'
      ierrorFVG=1


      end subroutine readFFRold
!
      end subroutine read_griddata

!
!+++++Created by T.Unemura 040521++++++++++++++++++++++++++++++++++++++
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine readFFRGrid_A(ifl,fnam,rtrnStat)
!     & (icpu,mvrtx,mcell,mface,nvrtx,ncell,nface,cord,
!     &  lacell,lvcell,rtrnStat)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_FFRGFwords
      use FFRdata, only : nvrtx,ncell,SFBOUN,NBOUND,IBFACE,IFFACE,
     &                    NFBOUN,NBFS,lvcell,lacell,cord,kmesh
!------------------------------------------------------------------------------------------
!    SFBOUN : Boundary region name read from grid file.
!    NBOUND : Number of boundary regions read from grid file
!    boundIDMap : Indexmap from FFR-Internal boundary region ID to grid file boundary ID
!    numbname : Maximum number of bounary region for one boundary condition (usually 2)
!    IFFACE : Vertices ID constructing boundary face.
!    NFBOUN : Number of faces in a boundary region read form grid file.
!------------------------------------------------------------------------------------------
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifl
      character*80,intent(in) :: fnam
      integer,intent(out)     :: rtrnStat
!
! --- [local entities]
!
c      character(len=80),save :: fnam
      integer :: ios=0
c      integer :: ifl,ios=0
      logical :: flag
      integer :: ic,jc
c      integer :: ic,jc,rtrnStat
      character*80 :: fileFormatID
      character*80 :: datasetID,commentStr,datatypeID,keywordStr,fieldID
      integer :: datasetSize
!
      integer :: intTmp1,intTmp2,intTmp3,intTmp4,intTmp5
      real(8) :: dblTmp1,dblTmp2,dblTmp3
      character*80 :: tmpStr80_1,tmpStr80_2
      character*255 :: tmpStr255
!------------------------------------------------------------------------
!     material ID and material name arrrays are alloated as local
!     variable in this version. These variable should be allocated
!     as module variable and shared in program.
!------------------------------------------------------------------------
      integer :: numofMaterial
      integer,allocatable :: material_ID(:)
      character*80,allocatable :: materialName(:)
!------------------------------------------------------------------------
!     number of boundary region. In this version, boundary definition
!     is not considered boundary condition so periodic or touchinlet 
!     boundary is defined with two boundary regions. In this case
!     numofBoundary is same with NBOUND
!------------------------------------------------------------------------
      integer :: numofBoundary
!------------------------------------------------------------------------
!     Maximum number of faces in a boundary region. This value is used
!     to allocate NFBOUN array. But this is not used when NFBOUN is refered.
!------------------------------------------------------------------------
      integer :: NBCE
!------------------------------------------------------------------------
!     These logical arrays are used to check for substituting values to 
!     all array elements.
!------------------------------------------------------------------------
      logical,allocatable :: chkflg_cord(:)
      logical,allocatable :: chkflg_IFFACE(:,:)
      logical,allocatable :: chkflg_lacell(:),chkflg_lvcell(:)
!
! --- 
!
! --- confirm logical file id ifl is not used.
      inquire(unit=ifl,opened=flag)
      if(flag) then
        write(*,*)'*** Error at (readFFRGrid)'
        write(*,*)'*** Unit number ',ifl,'is alread used.'
        rtrnStat=1
        return
      end if
!
! --- confirm the ffr-grid file exists.
!
      inquire(file=fnam,exist=flag)
      if(.not.(flag)) then
        write(*,*)'*** Error at (readFFRGrid)'
        write(*,*)'*** File ',trim(fnam),' does not exists.'
        rtrnStat=1
        return
      end if
!      
      open(ifl,file=fnam,form='formatted',
     &  status='unknown',action='read',iostat=ios)
      if(ios/=0) then
        write(*,*)'*** Cannot open FrontFlowRed Grid File:',
     &                                                trim(fnam)
        return
      end if
!
! --- read file format ID
!
      read(ifl,*) fileFormatID
      if(trim(adjustl(fileFormatID))==
     &            trim(adjustl(asciiv2_FFGF))) then
        call readFFRGridAGFv2(ios)
        if(ios/=0) then
          write(*,*)'*** Error in reading FFR-GRID'
          rtrnStat=1;return
        end if
      end if
      
      
      close(ifl)
      rtrnStat=0
      
      return
      
      contains
      
!////////////////////////////////////////////////////////////////
!========================================================
      subroutine readFFRGridAGFv2(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      integer :: ios2
      

!     read header field.
      print *, 'AAAA 01'
      read(ifl,*)datasetID
      if(trim(adjustl(datasetID))/=trim(adjustl(newSet_FFGF))) then
        write(*,*)'No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      read(ifl,'(a)')commentStr
      read(ifl,*)datatypeID
      read(ifl,*)keywordStr
      read(ifl,'(a)')commentStr
      read(ifl,*,iostat=ios2)datasetSize
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID header datasize'
          rtrnStat2=1;return
        end if
      read(ifl,*)fieldID
      if(trim(adjustl(fieldID))/=trim(adjustl(gridHead_FFGF))) then
        write(*,*)'***No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      call readFFRGridHeaderField(ios2)
      if(ios2/=0) then
        write(*,*)'***Error in reading header field.'
        rtrnStat2=1;return
      end if

!     allocate memory region
      call allocFFRGridMemory(ios2)
      if(ios2/=0) then
        write(*,*)'***Error in allocating memory.'
        rtrnStat2=1;return
      end if

!     Initialize grid data
      call initFFRGridMemory(ios2)
      if(ios2/=0) then
        write(*,*)'***Error in initializing memory.'
        rtrnStat2=1;return
      end if

!     read data field.
      do while(.true.)
        read(ifl,*)datatypeID
        if(trim(adjustl(datatypeID))==trim(adjustl(fileend_FFGF))) then
          rtrnStat2=0;  exit
        else if(trim(adjustl(datatypeID))==
     &                        trim(adjustl(newSet_FFGF))) then
          write(*,*)'***Second dataset is detected.'
!         In this version, FFR-GRID is construced by only one dataset
          rtrnStat2=1;  exit
        else if(trim(adjustl(datatypeID))==
     &                      trim(adjustl(customData_FFGF))) then
          read(ifl,*)keywordStr
          if(trim(adjustl(keywordStr))/=
     &                        trim(adjustl(ffrGrid_FFGF))) then
            write(*,*)'***Non grid datafield is deteced.'
!          In this version, field skipping function is not implemented.
            rtrnStat2=1;  exit
          end if
          read(ifl,'(a)')commentStr
          write(*,*)'    ###    Reading ',
     &            trim(adjustl(commentStr)),'...'
          read(ifl,*,iostat=ios2)datasetSize
            if(ios2/=0) then
              write(*,*)'*** Reading Error : FFR-GRID datasize'
              rtrnStat2=1;return
            end if
          read(ifl,*)fieldID
          if(trim(adjustl(fieldID))==
     &                  trim(adjustl(elemType_FFGF))) then
            call readFFRGridElemTypeField(ios2)
          else if(trim(adjustl(fieldID))==
     &                  trim(adjustl(bndryType_FFGF))) then
            call readFFRGridBndryTypeField(ios2)
          else if(trim(adjustl(fieldID))==
     &                  trim(adjustl(vrtxCord_FFGF))) then
            call readFFRGridVrtxCordField(ios2)
          else if(trim(adjustl(fieldID))==
     &                  trim(adjustl(bndryFace_FFGF))) then
            call readFFRGridBndryFaceField(ios2)
          else if(trim(adjustl(fieldID))==
     &                  trim(adjustl(elemVrtx_FFGF))) then
            call readFFRGridElemVrtxField(ios2)
          else
            write(*,*)'***Undefined datafield detected.',fieldID
            rtrnStat2=1;exit
          end if
          
        else
          write(*,*)'*** Undefined datatype is detected.',datatypeID
          rtrnStat2=1;  exit
        end if
      end do
      
!     check grid data
      if(rtrnStat2==0) then
        call checkFFRGridData(ios2)
        if(ios2/=0) then
          write(*,*)'*** Error in checking grid data.'
          rtrnStat2=1;return
        end if
      end if
      
      end subroutine readFFRGridAGFv2


!========================================================
      subroutine readFFRGridHeaderField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      character*80 :: fieldEndID
      
      read(ifl,*,iostat=ios2)fieldVersion
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID Header-1'
          rtrnStat2=1;return
        end if
      if(fieldVersion==1) then
!  ---  Read version 1 header field. FROM HERE -------------------------      
!       number of vertices
        read(ifl,*,iostat=ios2)intTmp1
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID header-2'
          rtrnStat2=1;return
        end if
c *** only in FFR2VIZ
        nvrtx=intTmp1
        if(nvrtx<1) then
          write(*,*)
     &      '*** Number of vertices are wrong. Check configuration.'
          write(*,*)'*** Control:',nvrtx,'  Grid:',intTmp1
          rtrnStat2=1;return
        end if
!       number of elements
        read(ifl,*,iostat=ios2)intTmp1
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID header-3'
          rtrnStat2=1;return
        end if
c *** only in FFR2VIZ
        ncell=intTmp1
        if(ncell<1) then
          write(*,*)
     &      '*** Number of elements are wrong. Check configuration.'
          write(*,*)'*** Control:',ncell,'  Grid:',intTmp1
          rtrnStat2=1;return
        end if
!       number of materials
        read(ifl,*,iostat=ios2)numofMaterial
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID header-4'
          rtrnStat2=1;return
        end if
!       number of boundary regions
        read(ifl,*,iostat=ios2)numofBoundary
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID header-5'
          rtrnStat2=1;return
        end if
        NBOUND=numofBoundary
!       number of boundary faces
        allocate(NFBOUN(1:NBOUND))
        NFBOUN=-1
        NBFS=0
        NBCE=0
        do ic2=1,NBOUND
          read(ifl,*,iostat=ios2)intTmp1,intTmp2
          if(ios2/=0) then
            write(*,*)'*** number of boundary region is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NBOUND)) then
            write(*,*)'*** Boundary region ID,',intTmp1,' is wrong.'
            rtrnStat2=1;return
          end if
          NFBOUN(intTmp1)=intTmp2
          NBFS=NBFS+NFBOUN(intTmp1)
        end do
        if(minval(NFBOUN)<0) then
          write(*,*)'*** Face number of boundary',minloc(NFBOUN),
     &                              ' is wrong:',minval(NFBOUN)
          rtrnStat2=1;return
        end if
!        NBFS=maxval(NFBOUN)
!  ---  Read version 1 header field. TO HERE -------------------------      
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(gridHeade_FFGF)))) then
        write(*,*)'*** Abnormal termination of header field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0
      
      end subroutine readFFRGridHeaderField

!-----------------------------------------------------------------------
!========================================================
      subroutine allocFFRGridMemory(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: ios2,ic2
      character*80 :: fieldEndID
      
!     material ID
      allocate(material_ID(1:numofMaterial))
!     material name
      allocate(materialName(1:numofMaterial))
      
!     *In this version material ID number and material name are 
!      included from FFR-GRID file but not sheared (only as local variable).
!      Only lacell data is effective as global data.

!     boundary type
!     boundary name1
!     boundary name2
!     *In this version boundary type and name2 are not
!      included from FFR-GRID file. Only name1 is loaded to SFBOUN array.
      allocate(SFBOUN(0:NBOUND))
!      SFBOUN number 0 element is setted to ' '

!     vertices coordinate
!     *vertices array is allocated external subroutine and passed as 
!      dummy argument.

!     face number on a bounary region.
!     *face number list is allocated in subroutine readFFRGridHeaderField.
!     vertices ID constructs boudary face.
      allocate(IBFACE(1:NBFS))
      allocate(IFFACE(1:4,1:NBFS))
!     vertices list of every element
!     material ID of elements
!     *vertices list lvcell and material ID lacell are allocated in external 
!      subroutine and passd as dummy argument.
!     cell type ID of elements
!     *In this version, cell type data is not included from FFR-GRID.

!     These logical arrays are used to check for substituting values to 
!     all array elements.
      allocate(chkflg_cord(1:nvrtx))
      allocate(chkflg_IFFACE(1:NBFS,1:NBOUND))
      allocate(chkflg_lacell(1:ncell),chkflg_lvcell(1:ncell))

c *** only in FFR2VIZ
c     allocate main arrays
      allocate(cord(1:3,1:nvrtx))
      allocate(lacell(1:ncell))
      allocate(lvcell(1:8,1:ncell))
      allocate(kmesh(1:ncell))

      rtrnStat2=0
      end subroutine allocFFRGridMemory

!-----------------------------------------------------------------------
!========================================================
      subroutine initFFRGridMemory(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: ios2,ic2
      character*80 :: fieldEndID
      
!     material ID
      material_ID=0
!     material name
      materialName=' '

!     boundary type
!     boundary name1
!     boundary name2
      SFBOUN(:)=' '

!     vertices coordinate
      cord(:,:)=0.d0

!     face number on a bounary region.
      IBFACE(:)=0
      IFFACE(:,:)=0

!     vertices list of every element
      lvcell(:,:)=0
!     material ID of elements
      lacell(:)=0
      
!     cell type ID of elements


!     check flag
      chkflg_cord=.false.
      chkflg_IFFACE=.false.
      chkflg_lacell=.false.
      chkflg_lvcell=.false.

c *** only in FFR2VIZ
      kmesh=0

      rtrnStat2=0
      end subroutine initFFRGridMemory

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridElemTypeField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2,jc2
      integer :: nmat_l

      character*80 :: fieldEndID
      
      read(ifl,*,iostat=ios2)fieldVersion
      if(ios2/=0) then
        write(*,*)'*** Reading Error : FFR-GRID Elemtype-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 Element type field FROM HERE -------------------
        read(ifl,*,iostat=ios2)nmat_l
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID Elemtype-2'
          rtrnStat2=1;return
        end if
        if((nmat_l<0).or.(nmat_l>numofMaterial)) then
          write(*,*)
     &      '*** Material number in element type field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nmat_l
          read(ifl,*,iostat=ios2)intTmp1,tmpStr80_1
          if(ios2/=0) then
            write(*,*)'*** line number of element type is wrong.'
            rtrnStat2=1;return
          end if
          if(intTmp1==0) then
            write(*,*)'*** Material ID number is ZERO.'
            rtrnStat2=1;return
          end if
!         search same ID definition or not defined position
          do jc2=1,numofMaterial
            if((material_ID(jc2)==intTmp1).or.(material_ID(jc2)==0)) 
     &                                                            exit
          end do
          material_ID(jc2)=intTmp2; materialName(jc2)=tmpStr80_1
        end do
!  ---  Read version 1 Element type field TO HERE -------------------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                  trim(adjustl(elemTypee_FFGF)))) then
        write(*,*)'*** Abnormal termination of element type field.'
        rtrnStat2=1;return
      end if  
      
      rtrnStat2=0

      end subroutine readFFRGridElemTypeField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridBndryTypeField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,jos2,ic2
      integer :: nbound_l
      
      character*80 :: fieldEndID
      
      read(ifl,*,iostat=ios2)fieldVersion
      if(ios2/=0) then
        write(*,*)'*** Reading Error : FFR-GRID Boundarytype-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 boundary type field. FROM HERE ------------------
        read(ifl,*,iostat=ios2)nbound_l
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID Boundarytype-2'
          rtrnStat2=1;return
        end if
        if((nbound_l<0).or.(nbound_l>numofBoundary)) then
          write(*,*)
     &      '*** Boundary number in boundary type field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nbound_l
          read(ifl,'(a)')tmpStr255
          read(tmpStr255,*,iostat=ios2)intTmp1,intTmp2,
     &                                    tmpStr80_1,tmpStr80_2
          if(ios2/=0) then
            read(tmpStr255,*,iostat=jos2)intTmp1,intTmp2,tmpStr80_1
            if(jos2/=0) then
              write(*,*)'*** line number of boundary type is wrong.'
              rtrnStat2=1;return
            end if
            tmpStr80_2=' '
          end if
          if((intTmp1<1).or.(intTmp1>numofBoundary)) then
            write(*,*)'*** Boundary ID number is wrong.'
            rtrnStat2=1;return
          end if
!         In this version, boundary condition ID (intTmp2) and 
!         boundary name2 (tmpStr80_2) is not used.
          SFBOUN(intTmp1)=tmpStr80_1
        end do
!  ---  Read version 1 boundary type field. TO HERE ------------------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(bndryTypee_FFGF)))) then
        write(*,*)'*** Abnormal termination of boundary type field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridBndryTypeField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridVrtxCordField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nvertex_l
      character*80 :: fieldEndID
      
      read(ifl,*,iostat=ios2)fieldVersion
      if(ios2/=0) then
        write(*,*)'*** Reading Error : FFR-GRID VertexCoordinate-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 vertices coordinate field. FROM HERE ---------
        read(ifl,*,iostat=ios2)nvertex_l
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID VertexCoordinate-2'
          rtrnStat2=1;return
        end if
        if((nvertex_l<0).or.(nvertex_l>nvrtx)) then
          write(*,*)'*** ',
     &        'Vertices number in vertex coordinate field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nvertex_l
          read(ifl,*,iostat=ios2)intTmp1,dblTmp1,dblTmp2,dblTmp3
          if(ios2/=0) then
            write(*,*)'*** line number of vertex coordinate is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>nvrtx)) then
            write(*,*)'*** Vertex ID number is wrong.'
            rtrnStat2=1;return
          end if
          cord(1,intTmp1)=dblTmp1
          cord(2,intTmp1)=dblTmp2
          cord(3,intTmp1)=dblTmp3
          chkflg_cord(intTmp1)=.true.
        end do
!  ---  Read version 1 vertices coordinate field. TO HERE -----------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(vrtxCorde_FFGF)))) then
        write(*,*)
     &      '*** Abnormal termination of vertices coordinates field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridVrtxCordField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridBndryFaceField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nface_l,nvrtxf_l,bndryTyp,bndryID
      character*80 :: fieldEndID
      
      read(ifl,*,iostat=ios2)fieldVersion
      if(ios2/=0) then
        write(*,*)'*** Reading Error : FFR-GRID BndryFace-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
! ---  Read version 1 boundary face field. FROM HERE -----------
        read(ifl,*,iostat=ios2)bndryTyp
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID BndryFace-2'
          rtrnStat2=1;return
        end if
        read(ifl,*,iostat=ios2)bndryID
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID BndryFace-3'
          rtrnStat2=1;return
        end if
        read(ifl,*,iostat=ios2)nface_l,nvrtxf_l
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID BndryFace-4'
          rtrnStat2=1;return
        end if
        if((bndryID<0).or.(bndryID>numofBoundary)) then
          write(*,*)'*** Boundary ID in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if((nface_l<0).or.(nface_l>NFBOUN(bndryID))) then
          write(*,*)'*** Face number in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if(nvrtxf_l/=4) then
          write(*,*)'*** Vertices number in a face is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nface_l
          read(ifl,*,iostat=ios2)
     &            intTmp1,intTmp2,intTmp3,intTmp4,intTmp5
          if(ios2/=0) then
            write(*,*)'*** line number of boundary face is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NFBOUN(bndryID))) then
            write(*,*)'*** Face ID number is wrong.'
            rtrnStat2=1;return
          end if
          NBCE=NBCE+1
          IBFACE(NBCE)= bndryID
          IFFACE(1,NBCE)=intTmp2
          IFFACE(2,NBCE)=intTmp3
          IFFACE(3,NBCE)=intTmp4
          IFFACE(4,NBCE)=intTmp5
          chkflg_IFFACE(intTmp1,bndryID)=.true.
        end do
! ---  Read version 1 boundary face field. TO HERE -----------
      else
        write(*,*)'****Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(bndryFacee_FFGF)))) then
        write(*,*)'*** Abnormal termination of boundary face field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridBndryFaceField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridElemVrtxField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nelem_l,vertxTmp(1:8)
      character*80 :: fieldEndID
c *** only in FFR2VIZ
      integer :: j,k

      read(ifl,*,iostat=ios2)fieldVersion
      if(ios2/=0) then
        write(*,*)'*** Reading Error : FFR-GRID Element-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 element field. FROM HERE -----------------
        read(ifl,*,iostat=ios2)nelem_l
        if(ios2/=0) then
          write(*,*)'*** Reading Error : FFR-GRID Element-2'
          rtrnStat2=1;return
        end if
        if((nelem_l<0).or.(nelem_l>ncell)) then
          write(*,*)'*** Element ID in element field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nelem_l
          read(ifl,*,iostat=ios2)intTmp1,vertxTmp(1:8),intTmp2,intTmp3
          if(ios2/=0) then
            write(*,*)'*** line number of element is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>ncell)) then
            write(*,*)'*** Element ID number is wrong.'
            rtrnStat2=1;return
          end if
c *** only in FFR2VIZ
          k=0
          do j=1,8
             if(vertxTmp(j)>0) then
                k=k+1
                lvcell(k,intTmp1)=vertxTmp(j)
             end if
          end do
          select case(k)
          case(4)
             kmesh(intTmp1) = 1 ! tetra
          case(8)
             kmesh(intTmp1) = 2 ! hexa
          case(6)
             kmesh(intTmp1) = 3 ! prism
          case(5)
             kmesh(intTmp1) = 4 ! pyramid
          end select
          lacell(intTmp1)=intTmp2
!          In this version cell type  data is not included from FFR-GRID.
          chkflg_lacell(intTmp1)=.true.
          chkflg_lvcell(intTmp1)=.true.
        end do
!  ---  Read version 1 element field. FROM TO -----------------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(elemVrtxe_FFGF)))) then
        write(*,*)'*** Abnormal termination of element field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridElemVrtxField


!-----------------------------------------------------------------------
!========================================================
      subroutine checkFFRGridData(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: ios2,ic2,jc2
      
!     check all elements of arrays are defined.
      do ic2=1,nvrtx
        if(.not.(chkflg_cord(ic2))) then
          write(*,*)'*** Vertex',ic2,'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      deallocate(chkflg_cord)
      
      do jc2=1,NBOUND
      do ic2=1,NFBOUN(jc2)
        if(.not.(chkflg_IFFACE(ic2,jc2))) then
          write(*,*)'*** Boundary face',ic2,'of boundary',jc2,
     &                                          'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      end do
      deallocate(chkflg_IFFACE)
      
      do ic2=1,ncell
        if(.not.(chkflg_lacell(ic2))) then
          write(*,*)'*** Material ID of cell',ic,'is not defined.'
          rtrnStat2=1;return
        end if
        if(.not.(chkflg_lvcell(ic2))) then
          write(*,*)'*** Vertices ID of cell',ic,'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      deallocate(chkflg_lacell)
      deallocate(chkflg_lvcell)

      rtrnStat2=0
      end subroutine checkFFRGridData
!
      end subroutine readFFRGrid_A
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine readFFRGrid(ifl,fnam,rtrnStat)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_FFRGFwords
      use FFRdata, only : nvrtx,ncell,SFBOUN,NBOUND,IBFACE,IFFACE,
     &                    NFBOUN,NBFS,lvcell,lacell,cord,kmesh
!----------------------------------------------------------------------------------------
!    SFBOUN : Boundary region name read from grid file.
!    NBOUND : Number of boundary regions read from grid file
!    boundIDMap : Indexmap from FFR-Internal boundary region ID to grid file boundary ID
!    numbname : Maximum number of bounary region for one boundary condition (usually 2)
!    IFFACE : Vertices ID constructing boundary face.
!    NFBOUN : Number of faces in a boundary region read form grid file.
!----------------------------------------------------------------------------------------
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifl
      character*80,intent(in) :: fnam
      integer,intent(out)     :: rtrnStat
!
! --- [local entities]
!
c      character(len=80),save :: fnam
      integer :: ios=0
c      integer :: ifl,ios=0
      logical :: flag
      integer :: ic,jc
c      integer :: ic,jc,rtrnStat
      character*80 :: fileFormatID
      character*80 :: datasetID,commentStr,datatypeID,keywordStr,fieldID
      integer :: datasetSize

      integer :: intTmp1,intTmp2,intTmp3,intTmp4,intTmp5
      real(8) :: dblTmp1,dblTmp2,dblTmp3
      character*80 :: tmpStr80_1,tmpStr80_2
!---------------------------------------------------------------------
!     material ID and material name arrrays are alloated as local
!     variable in this version. These variable should be allocated
!     as module variable and shared in program.
!---------------------------------------------------------------------
      integer :: numofMaterial
      integer,allocatable :: material_ID(:)
      character*80,allocatable :: materialName(:)
!---------------------------------------------------------------------
!     number of boundary region. In this version, boundary definition
!     is not considered boundary condition so periodic or touchinlet 
!     boundary is defined with two boundary regions. In this case
!     numofBoundary is same with NBOUND
!---------------------------------------------------------------------
      integer :: numofBoundary
!---------------------------------------------------------------------
!     Maximum number of faces in a boundary region. This value is used
!     to allocate NFBOUN array. But this is not used when NFBOUN is refered.
!---------------------------------------------------------------------
      integer :: NBCE
!---------------------------------------------------------------------
!     These logical arrays are used to check for substituting values to 
!     all array elements.
!---------------------------------------------------------------------
      logical,allocatable :: chkflg_cord(:)
      logical,allocatable :: chkflg_IFFACE(:,:)
      logical,allocatable :: chkflg_lacell(:),chkflg_lvcell(:)

!     confirm logical file id ifl is not used.
      inquire(unit=ifl,opened=flag)
      if(flag) then
        write(*,*)'*** Error at (readFFRGrid)'
        write(*,*)'*** Unit number ',ifl,'is alread used.'
        rtrnStat=1
        return
      end if
!     confirm the ffr-grid file exists.
      inquire(file=fnam,exist=flag)
      if(.not.(flag)) then
        write(*,*)'*** Error at (readFFRGrid)'
        write(*,*)'*** File ',trim(fnam),' does not exists.'
        rtrnStat=1
        return
      end if
      
      open(ifl,file=fnam,form='unformatted',
     &  status='unknown',action='read',iostat=ios)
      if(ios/=0) then
        write(*,*)'*** Cannot open FrontFlowRed Grid File:',
     &                                                trim(fnam)
        return
      end if

!     read file format ID
      read(ifl)fileFormatID
      if(fileFormatID==unformv2_FFGF) then
        call readFFRGridUGFv2(ios)
        if(ios/=0) then
          write(*,*)'*** Error in reading FFR-GRID'
          rtrnStat=1;return
        end if
      end if
      
      
      close(ifl)
      rtrnStat=0
      
      return
      
      contains
      
!///////////////////////////////////////////////////////////////////
!========================================================
      subroutine readFFRGridUGFv2(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      integer :: ios2
!
! --- read header field.
!
      read(ifl)datasetID
      if(datasetID/=newSet_FFGF) then
        write(*,*)'No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      read(ifl)commentStr
      read(ifl)datatypeID
      read(ifl)keywordStr
      read(ifl)commentStr
      read(ifl)datasetSize
      read(ifl)fieldID
      if(fieldID/=gridHead_FFGF) then
        write(*,*)'***No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      call readFFRGridHeaderField(ios2)
      if(ios2/=0) then
        write(*,*)'***Error in reading header field.'
        rtrnStat2=1;return
      end if

!     allocate memory region
      call allocFFRGridMemory(ios2)
      if(ios2/=0) then
        write(*,*)'***Error in allocating memory.'
        rtrnStat2=1;return
      end if

!     Initialize grid data
      call initFFRGridMemory(ios2)
      if(ios2/=0) then
        write(*,*)'***Error in initializing memory.'
        rtrnStat2=1;return
      end if
!
! --- read data field.
!
      do while(.true.)
        read(ifl)datatypeID
        if(datatypeID==fileend_FFGF) then
          rtrnStat2=0;  exit
        else if(datatypeID==newSet_FFGF) then
          write(*,*)'***Second dataset is detected.'
!
! --- In this version, FFR-GRID is construced by only one dataset
!
          rtrnStat2=1;  exit
        else if(datatypeID==customData_FFGF) then
          read(ifl)keywordStr
          if(keywordStr/=ffrGrid_FFGF) then
            write(*,*)'***Non grid datafield is deteced.'
!
! --- In this version, field skipping function is not implemented.
!
            rtrnStat2=1;  exit
          end if
          read(ifl)commentStr
          write(*,*)'    ###    Reading ',
     &            trim(adjustl(commentStr)),'...'
          read(ifl)datasetSize
          read(ifl)fieldID
          if(fieldID==elemType_FFGF) then
            call readFFRGridElemTypeField(ios2)
          else if(fieldID==bndryType_FFGF) then
            call readFFRGridBndryTypeField(ios2)
          else if(fieldID==vrtxCord_FFGF) then
            call readFFRGridVrtxCordField(ios2)
          else if(fieldID==bndryFace_FFGF) then
            call readFFRGridBndryFaceField(ios2)
          else if(fieldID==elemVrtx_FFGF) then
            call readFFRGridElemVrtxField(ios2)
          else
            write(*,*)'***Undefined datafield detected.',fieldID
            rtrnStat2=1;exit
          end if
          
        else
          write(*,*)'*** Undefined datatype is detected.',datatypeID
          rtrnStat2=1;  exit
        end if
      end do
      
!     check grid data
      if(rtrnStat2==0) then
        call checkFFRGridData(ios2)
        if(ios2/=0) then
          write(*,*)'*** Error in checking grid data.'
          rtrnStat2=1;return
        end if
      end if
      
      end subroutine readFFRGridUGFv2

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridHeaderField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      character*80 :: fieldEndID
      
      read(ifl)fieldVersion
      if(fieldVersion==1) then
!
! --- Read version 1 header field. FROM HERE -------------------------      
!
! --- number of vertices
        read(ifl)intTmp1
c *** only in FFR2VIZ
        nvrtx=intTmp1
        if(nvrtx<1) then
          write(*,*)
     &      '*** Number of vertices are wrong. Check configuration.'
          write(*,*)'*** Control:',nvrtx,'  Grid:',intTmp1
          rtrnStat2=1;return
        end if
! --- number of elements
        read(ifl)intTmp1
c *** only in FFR2VIZ
        ncell=intTmp1
        if(ncell<1) then
          write(*,*)
     &      '*** Number of elements are wrong. Check configuration.'
          write(*,*)'*** Control:',ncell,'  Grid:',intTmp1
          rtrnStat2=1;return
        end if
!       number of materials
        read(ifl)numofMaterial
!       number of boundary regions
        read(ifl)numofBoundary
        NBOUND=numofBoundary
!       number of boundary faces
        allocate(NFBOUN(1:NBOUND))
        NFBOUN=-1
        NBFS=0
        NBCE=0
        do ic2=1,NBOUND
          read(ifl,iostat=ios2)intTmp1,intTmp2
          if(ios2/=0) then
            write(*,*)'*** number of boundary region is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NBOUND)) then
            write(*,*)'*** Boundary region ID,',intTmp1,' is wrong.'
            rtrnStat2=1;return
          end if
          NFBOUN(intTmp1)=intTmp2
          NBFS=NBFS+NFBOUN(intTmp1)
        end do
        if(minval(NFBOUN)<0) then
          write(*,*)'*** Face number of boundary',minloc(NFBOUN),
     &                              ' is wrong:',minval(NFBOUN)
          rtrnStat2=1;return
        end if
!        NBFS=maxval(NFBOUN)
!  ---  Read version 1 header field. TO HERE -------------------------      
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=gridHeade_FFGF)) then
        write(*,*)'*** Abnormal termination of header field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0
      
      end subroutine readFFRGridHeaderField

!-----------------------------------------------------------------------
!========================================================
      subroutine allocFFRGridMemory(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: ios2,ic2
      character*80 :: fieldEndID
      
!     material ID
      allocate(material_ID(1:numofMaterial))
!     material name
      allocate(materialName(1:numofMaterial))
      
!     *In this version material ID number and material name are 
!      included from FFR-GRID file but not sheared (only as local variable).
!      Only lacell data is effective as global data.

!     boundary type
!     boundary name1
!     boundary name2
!     *In this version boundary type and name2 are not
!      included from FFR-GRID file. Only name1 is loaded to SFBOUN array.
      allocate(SFBOUN(0:NBOUND))
!      SFBOUN number 0 element is setted to ' '

!     vertices coordinate
!     *vertices array is allocated external subroutine and passed as 
!      dummy argument.

!     face number on a bounary region.
!     *face number list is allocated in subroutine readFFRGridHeaderField.
!     vertices ID constructs boudary face.
      allocate(IBFACE(1:NBFS))
      allocate(IFFACE(1:4,1:NBFS))

!     vertices list of every element
!     material ID of elements
!     *vertices list lvcell and material ID lacell are allocated in external 
!      subroutine and passd as dummy argument.
!     cell type ID of elements
!     *In this version, cell type data is not included from FFR-GRID.

!     These logical arrays are used to check for substituting values to 
!     all array elements.
      allocate(chkflg_cord(1:nvrtx))
      allocate(chkflg_IFFACE(1:NBFS,1:NBOUND))
      allocate(chkflg_lacell(1:ncell),chkflg_lvcell(1:ncell))

c *** only in FFR2VIZ
c     allocate main arrays
      allocate(cord(1:3,1:nvrtx))
      allocate(lacell(1:ncell))
      allocate(lvcell(1:8,1:ncell))
      allocate(kmesh(1:ncell))

      rtrnStat2=0
      end subroutine allocFFRGridMemory

!-----------------------------------------------------------------------
!========================================================
      subroutine initFFRGridMemory(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: ios2,ic2
      character*80 :: fieldEndID
      
!     material ID
      material_ID=0
!     material name
      materialName=' '

!     boundary type
!     boundary name1
!     boundary name2
      SFBOUN(:)=' '

!     vertices coordinate
      cord(:,:)=0.d0

!     face number on a bounary region.
      IBFACE(:)=0
      IFFACE(:,:)=0

!     vertices list of every element
      lvcell(:,:)=0
!     material ID of elements
      lacell(:)=0
      
!     cell type ID of elements


!     check flag
      chkflg_cord=.false.
      chkflg_IFFACE=.false.
      chkflg_lacell=.false.
      chkflg_lvcell=.false.

c *** only in FFR2VIZ
      kmesh=0

      rtrnStat2=0
      end subroutine initFFRGridMemory

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridElemTypeField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2,jc2
      integer :: nmat_l

      character*80 :: fieldEndID
      
      read(ifl)fieldVersion
      if(fieldVersion==1) then
!  ---  Read version 1 Element type field FROM HERE -------------------
        read(ifl)nmat_l
        if((nmat_l<0).or.(nmat_l>numofMaterial)) then
          write(*,*)
     &      '*** Material number in element type field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nmat_l
          read(ifl,iostat=ios2)intTmp1,tmpStr80_1
          if(ios2/=0) then
            write(*,*)'*** line number of element type is wrong.'
            rtrnStat2=1;return
          end if
          if(intTmp1==0) then
            write(*,*)'*** Material ID number is ZERO.'
            rtrnStat2=1;return
          end if
!         search same ID definition or not defined position
          do jc2=1,numofMaterial
            if((material_ID(jc2)==intTmp1).or.(material_ID(jc2)==0)) 
     &                                                            exit
          end do
          material_ID(jc2)=intTmp2; materialName(jc2)=tmpStr80_1
        end do
!  ---  Read version 1 Element type field TO HERE -------------------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=elemTypee_FFGF)) then
        write(*,*)'*** Abnormal termination of element type field.'
        rtrnStat2=1;return
      end if  
      
      rtrnStat2=0

      end subroutine readFFRGridElemTypeField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridBndryTypeField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nbound_l
      
      character*80 :: fieldEndID
      
      read(ifl)fieldVersion
      if(fieldVersion==1) then
!  ---  Read version 1 boundary type field. FROM HERE ------------------
        read(ifl)nbound_l
        if((nbound_l<0).or.(nbound_l>numofBoundary)) then
          write(*,*)
     &      '*** Boundary number in boundary type field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nbound_l
          read(ifl,iostat=ios2)intTmp1,intTmp2,tmpStr80_1,tmpStr80_2
          if(ios2/=0) then
            write(*,*)'*** line number of boundary type is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>numofBoundary)) then
            write(*,*)'*** Boundary ID number is wrong.'
            rtrnStat2=1;return
          end if
!         In this version, boundary condition ID (intTmp2) and 
!         boundary name2 (tmpStr80_2) is not used.
          SFBOUN(intTmp1)=tmpStr80_1
        end do
!  ---  Read version 1 boundary type field. TO HERE ------------------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=bndryTypee_FFGF)) then
        write(*,*)'*** Abnormal termination of boundary type field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridBndryTypeField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridVrtxCordField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nvertex_l
      character*80 :: fieldEndID
      
      
      read(ifl)fieldVersion
      if(fieldVersion==1) then
!  ---  Read version 1 vertices coordinate field. FROM HERE ---------
        read(ifl)nvertex_l
        if((nvertex_l<0).or.(nvertex_l>nvrtx)) then
          write(*,*)'*** ',
     &        'Vertices number in vertex coordinate field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nvertex_l
          read(ifl,iostat=ios2)intTmp1,dblTmp1,dblTmp2,dblTmp3
          if(ios2/=0) then
            write(*,*)'*** line number of vertex coordinate is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>nvrtx)) then
            write(*,*)'*** Vertex ID number is wrong.'
            rtrnStat2=1;return
          end if
          cord(1,intTmp1)=dblTmp1
          cord(2,intTmp1)=dblTmp2
          cord(3,intTmp1)=dblTmp3
          chkflg_cord(intTmp1)=.true.
        end do
!  ---  Read version 1 vertices coordinate field. TO HERE -----------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=vrtxCorde_FFGF)) then
        write(*,*)
     &      '*** Abnormal termination of vertices coordinates field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridVrtxCordField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridBndryFaceField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nface_l,nvrtxf_l,bndryTyp,bndryID
      character*80 :: fieldEndID
      
      read(ifl)fieldVersion
      if(fieldVersion==1) then
! ---  Read version 1 boundary face field. FROM HERE -----------
        read(ifl)bndryTyp
        read(ifl)bndryID
        read(ifl)nface_l,nvrtxf_l
        if((bndryID<0).or.(bndryID>numofBoundary)) then
          write(*,*)'*** Boundary ID in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if((nface_l<0).or.(nface_l>NFBOUN(bndryID))) then
          write(*,*)'*** Face number in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if(nvrtxf_l/=4) then
          write(*,*)'*** Vertices number in a face is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nface_l
          read(ifl,iostat=ios2)intTmp1,intTmp2,intTmp3,intTmp4,intTmp5
          if(ios2/=0) then
            write(*,*)'*** line number of boundary face is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NFBOUN(bndryID))) then
            write(*,*)'*** Face ID number is wrong.'
            rtrnStat2=1;return
          end if
          NBCE=NBCE+1
          IBFACE(NBCE)= bndryID
          IFFACE(1,NBCE)=intTmp2
          IFFACE(2,NBCE)=intTmp3
          IFFACE(3,NBCE)=intTmp4
          IFFACE(4,NBCE)=intTmp5
          chkflg_IFFACE(intTmp1,bndryID)=.true.
        end do
! ---  Read version 1 boundary face field. TO HERE -----------
      else
        write(*,*)'****Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=bndryFacee_FFGF)) then
        write(*,*)'*** Abnormal termination of boundary face field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridBndryFaceField

!-----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridElemVrtxField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nelem_l,vertxTmp(1:8)
      character*80 :: fieldEndID
c *** only in FFR2VIZ
      integer :: j,k

      read(ifl)fieldVersion
      if(fieldVersion==1) then
!  ---  Read version 1 element field. FROM HERE -----------------
        read(ifl)nelem_l
        if((nelem_l<0).or.(nelem_l>ncell)) then
          write(*,*)'*** Element ID in element field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nelem_l
          read(ifl,iostat=ios2)intTmp1,vertxTmp(1:8),intTmp2,intTmp3
          if(ios2/=0) then
            write(*,*)'*** line number of element is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>ncell)) then
            write(*,*)'*** Element ID number is wrong.'
            rtrnStat2=1;return
          end if
c *** only in FFR2VIZ
          k=0
          do j=1,8
             if(vertxTmp(j)>0) then
                k=k+1
                lvcell(k,intTmp1)=vertxTmp(j)
             end if
          end do
          select case(k)
          case(4)
             kmesh(intTmp1) = 1 ! tetra
          case(8)
             kmesh(intTmp1) = 2 ! hexa
          case(6)
             kmesh(intTmp1) = 3 ! prism
          case(5)
             kmesh(intTmp1) = 4 ! pyramid
          end select
          lacell(intTmp1)=intTmp2
!          In this version cell type  data is not included from FFR-GRID.
          chkflg_lacell(intTmp1)=.true.
          chkflg_lvcell(intTmp1)=.true.
        end do
!  ---  Read version 1 element field. FROM TO -----------------
      else
        write(*,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=elemVrtxe_FFGF)) then
        write(*,*)'*** Abnormal termination of element field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridElemVrtxField


!-----------------------------------------------------------------------
!========================================================
      subroutine checkFFRGridData(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: ios2,ic2,jc2
      
!     check all elements of arrays are defined.
      do ic2=1,nvrtx
        if(.not.(chkflg_cord(ic2))) then
          write(*,*)'*** Vertex',ic2,'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      deallocate(chkflg_cord)
      
      do jc2=1,NBOUND
      do ic2=1,NFBOUN(jc2)
        if(.not.(chkflg_IFFACE(ic2,jc2))) then
          write(*,*)'*** Boundary face',ic2,'of boundary',jc2,
     &                                          'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      end do
      deallocate(chkflg_IFFACE)
      
      do ic2=1,ncell
        if(.not.(chkflg_lacell(ic2))) then
          write(*,*)'*** Material ID of cell',ic,'is not defined.'
          rtrnStat2=1;return
        end if
        if(.not.(chkflg_lvcell(ic2))) then
          write(*,*)'*** Vertices ID of cell',ic,'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      deallocate(chkflg_lacell)
      deallocate(chkflg_lvcell)

      rtrnStat2=0
      end subroutine checkFFRGridData


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      end subroutine readFFRGrid
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!++++++Modified by T.Unemura 040531+++++++++++++++++++++++++++++++++++++
