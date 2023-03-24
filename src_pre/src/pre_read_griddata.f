!     subroutine read_grid_data
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_grid_data
     &    (mcell,mface,mvrtx,ncell,nvrtx,cord,lacell,
     &     lvcell,lfcell,lcface,lvface,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only       : lenfnm,getfil,gdformat,
     &                           ffrgridOutFlag,ffgFormatKey,
     &                           ifli,ifll,ifle,cntlnam
     &                           ,multiMaterialFlag
      use module_gf
      use module_boundary,only : nbcnd,boundName,numbname,ivbcnd,
     &                           lvbcnd,icbinr,lcbinr,kdintr,
     &                           kdbcnd,kdprdc,kdtchi,kdsld,
     &                           prdcAxis,prdcOffset,
     &                           boundIDMap,SFBOUN,NBOUND,
     &                           IFFACE,NFBOUN,NAME_I,
     &                           dum_bc
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mcell,mface,mvrtx,ncell
      integer,intent(inout) :: nvrtx
      real*8 ,intent(out)   :: cord(3,mvrtx)
      integer,intent(out)   :: lacell(  mcell)
      integer,intent(out)   :: lvcell(8,mcell)
      integer,intent(inout) :: lfcell(7,mcell)  
      integer,intent(inout) :: lcface(2,mface)  
      integer,intent(inout) :: lvface(4,mface)  
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      character(lenfnm)   :: fnam
      logical             :: fexist=.false.
      integer             :: k,n,ifl,ios,ierr1
      integer             :: no,kv,kc,kvs,kve,kcs,kce
      integer             :: nrec,nbdmn,nvbdmn,ncbdmn
      integer             :: i,j,m,mb,MXNODE,MXELEM,NB
      CHARACTER*60        :: FV_NAME,strID
      integer             :: major_version,minor_version,ifvreg
      character(len=10) :: featureName = ""
      character(len=0)  :: version = ""
      integer :: chos_mshr=1
      integer, parameter :: maxboun=100, maxdmn=100
! DEBUG
!      integer,parameter   :: lvfcel(4,6,4)=reshape( source=
!     &  (/1,3,2,0, 2,3,4,0, 3,1,4,0, 4,1,2,0, 0,0,0,0, 0,0,0,0,
!     &    1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0, 0,0,0,0,
!     &    1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2,0, 4,5,6,0, 0,0,0,0,
!     &    1,5,8,4, 2,3,7,6, 1,2,6,5, 3,4,8,7, 1,4,3,2, 5,6,7,8/),
!     &  shape=(/4,6,4/))
!      INTEGER,PARAMETER :: LVEDGL(2,4,2)=RESHAPE( SOURCE=
!     &            (/1,2, 2,3, 3,1, 0,0,
!     &              1,2, 2,3, 3,4, 4,1/),SHAPE=(/2,4,2/))
! --- only for hex:
!      INTEGER,PARAMETER :: lfedg(4,6)=RESHAPE( SOURCE=
!     &            (/1,5, 5,8, 8,4, 4,1,
!     &              2,3, 3,7, 7,6, 6,2,
!     &              1,2, 2,6, 6,5, 5,1,
!     &              3,4, 4,8, 8,7, 7,3,
!     &              1,4, 4,3, 3,2, 2,1,
!     &              5,6, 6,7, 7,8, 8,5/),SHAPE=(/4,6/))
!
!
      lacell=0
      lvcell=0
!
      ierror=0
      if(gdformat=='FV') then
!------------------------------------------------------
! --- read FieldView grid file greated by Gridgen  ----
!------------------------------------------------------
        call getfil(ifl,fnam,'vertex')
        fexist=.false.
        inquire(file=fnam,exist=fexist)
!
        if(.not.fexist) then
          write(ifle,'(1x,2a)') 'ERR: NOT found file: ',
     &                      TRIM(adjustl(fnam))
          stop 'in read_griddata'
        endif
!
        OPEN(ifl,FILE=fnam,FORM='FORMATTED',action='read',iostat=ios)
!
	if(ios/=0) then
	  write(ifle,'(1x,3a)')
     &    'MSG:  Cannot open file :',trim(fnam),'.'
	  ierror=1
	  return
	end if
!
!       check FV/UNS version
        READ(ifl,'(A60)')  FV_NAME
        IF(FV_NAME(1:9) .ne. 'FIELDVIEW') THEN
	  WRITE(ifle,'(1x,a)') 
     &      'MSG: This is NOT a FIELDVIEW ASCII file'
          ierror=1
	  return
        END IF
        rewind(ifl)
        read(FV_NAME,*) strID,major_version,minor_version
        if(major_version==3) then
           fnam=trim(fnam) // '.fvreg'
           ifvreg=ifl+1
           OPEN(ifvreg,FILE=fnam,FORM='FORMATTED',
     &                 action='read',iostat=ios)
           if(ios/=0) then
	     write(ifle,'(1x,3a)')
     &       'WRG:  Region file cannot open:',trim(fnam)
             ifvreg=-1             
           end if
           call readFVGrids_30(ifl,ifvreg,ios)
           if(ifvreg/=-1) then
              close(ifvreg)
           end if
        elseif(major_version==2 .and. minor_version==4) then
           call readFVGrids(ifl,ios)
        end if
!
	if(ios/=0) then
	  ierror=1
          return
	endif
        close(ifl)
!-------------------
!#ifdef STARCD
!-------------------
      else if(gdformat=='SC') then
      else if(gdformat=='GF') then
!--------------------------------------------------
! --- read GF grid file greated by Gridgen    -----
!--------------------------------------------------
!	call read_grid(mcell,mvrtx,ncell,nvrtx,cord,lacell,lvcell,ios)
!	if(ios/=0) then
!	  ierror=1
!	  return
!	end if
!	call read_boundary(mvrtx,ncell,ios)
!	if(ios/=0) then
!	  ierror=1
!	  return
!	end if
!
        if(ffgFormatKey=='A') then
          write(ifll,'(1x,a)') 'MSG: Begin read FFR ASCII grids' 
          call readFFRGrid_A
     &      (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell,ios)
        else if(ffgFormatKey=='U') then
          write(ifll,'(1x,a)') 'MSG: Begin read FFR UNFORMAT grids' 
          call readFFRGrid
     &      (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell,ios)
        else
          write(ifle,'(1x,a)') 'ERR: Not yet support FFR BINARY grids' 
          ierror=1;return
        end if
        write(ifll,'(1x,a)') 'MSG: End read FFR grids' 
!
        if(ios/=0) then
          ierror=1;return
        endif
!#ifdef GAMBIT
! DEBUG
      else if (gdformat=='GB') then
!-----------------------------------------------------------
! --- read GAMBIT grid file created by Gridgen -------------   
!-----------------------------------------------------------
         featureName = "FFR-GB"//char(0)
         version     = version  //char(0)

        call getfil(ifl,fnam,'fluent_grd')
        open(ifl,file=fnam,form='formatted',action='read',iostat=ios)
        if (ios/=0) then
          write(ifle,'(1x,3a)')
     &     'MSG:  Cannot open file :',trim(fnam)
          ierror=1
          return
        end if
!
        write(ifll,'(1x,a)') 'MSG: Begin read GAMBIT grids' 
        call readGAMBITgrids(
     &     mcell,mface,mvrtx,ncell,nvrtx,cord,lacell,
     &     lvcell,lfcell,lcface,lvface,
     &     ifl,ios)
        write(ifll,*) 'End read GAMBIT grids' 
        if(ios/=0) then
          ierror=1
          return
        end if
        close(ifl)
!
c      else if (gdformat=='GB-T') then
c!-----------------------------------------------------------
c! --- read GAMBIT grid file created by Gridgen -------------   
c!-----------------------------------------------------------
c        call getfil(ifl,fnam,'fluent_grd')
c        open(ifl,file=fnam,form='formatted',action='read',iostat=ios)
c        if (ios/=0) then
c          write(ifle,'(1x,3a)')
c     &     'MSG:  Cannot open file :',trim(fnam)
c          ierror=1
c          return
c        end if
c!
c        write(ifll,'(1x,a)') 'MSG: Begin read GAMBIT TETRA grids' 
c        call readGAMBITgrids_T(ifl,ios)
c        write(ifll,*) 'End read GAMBIT TETRA grids' 
c        if(ios/=0) then
c          ierror=1
c          return
c        end if
c        close(ifl)
c!
c      else if (gdformat=='GB-H') then
c!-----------------------------------------------------------
c! --- read GAMBIT grid file created by Gridgen -------------   
c!-----------------------------------------------------------
c        call getfil(ifl,fnam,'fluent_grd')
c        open(ifl,file=fnam,form='formatted',action='read',iostat=ios)
c        if (ios/=0) then
c          write(ifle,'(1x,2a)')
c     &     'MSG:  Cannot open file :',trim(fnam)
c          ierror=1
c          return
c        end if
c!
c        write(ifll,'(1x,a)') 'MSG: Begin read GAMBIT HEXA grids' 
c        call readGAMBITgrids_H(ifl,ios)
c        write(ifll,*) 'End read GAMBIT HEXA grids' 
c        if(ios/=0) then
c          ierror=1
c          return
c        end if
c        close(ifl)
c!
!#endif
!
!#ifdef SCRYU
      else if (gdformat=='SY') then
!----------------------------------------------
! --- ---- read SCRYU/Tetra grid file ---------      
!----------------------------------------------
         featureName = "FFR-SY"//char(0)
         version     = version  //char(0)

!
c
        call getfil(ifl,fnam,'scryu_grd')
        open(ifl,file=fnam,form='unformatted',action='read',iostat=ios)
        if (ios/=0) then
          write(ifle,'(1x,2a)')
     &     'MSG:  Cannot open file :',trim(fnam)
          ierror=1
          return
        end if
c
        write(ifll,'(1x,a)') 'Begin read SCRYU grids' 
        call readSCRYUgrids(ifl,ios)
        write(ifll,'(1x,a)') 'MSG: End read SCRYU grids'
        if(ios/=0) then
          ierror=1
          return
        end if
        close(ifl)
!
!#endif
      else if(gdformat=='NST') then
!-----------------------------------------------------------
! --- read NASTRAN grid file created by Gridgen -------------   
!-----------------------------------------------------------
         featureName = "FFR-NST"//char(0)
         version     = version   //char(0)

!
        call getfil(ifl,fnam,'nastran_grd')
        fexist=.false.
        inquire(file=fnam,exist=fexist)
        if(.not.fexist) then
          write(ifle,'(1x,2a)') 'ERR: NOT found file: ',
     &                      TRIM(adjustl(fnam))
          stop 'in read_griddata'
        endif
        open(ifl,file=fnam,form='formatted',action='read',iostat=ios)
        if (ios/=0) then
          write(ifle,'(1x,2a)')
     &     'MSG:  Cannot open file :',trim(fnam)
          ierror=1
          stop
        end if
        write(ifll,'(1x,a)') 'MSG: Begin read NASTRAN grids'
!
        call readNASTRANgrids(ifl,ios)
        write(ifll,'(1x,a)') 'MSG: End read NASTRAN grids'
        if(ios/=0) then
          ierror=1
          return
        end if
        close(ifl)
      else
!-------------------------------
! --- NOT support grid format
!-------------------------------
	write(ifle,'(1x,2a)')
     &   'ERR: Grid Format is NOT supported : ',trim(gdformat)
	ierror=1
	return
      end if
!
      DO mb=1,NBOUND
      SFBOUN(mb)=trim(adjustl(SFBOUN(mb)))
      enddo
!
      return
!
!//////////////////////////////////////////////////////////////////
      contains
!==================================================================
      subroutine readFVGrids_30(IUTFV,IUTRG,ierrorFVG)
!==================================================================
!     This subroutine read 'FieldView' grid format 
!     which is created by 'Gridgen'
!==================================================================
!
! --- [module arguments]
!
      integer,intent(in)     :: IUTFV,IUTRG
      integer,intent(inout)  :: ierrorFVG
!
! --- [local entities]
!
      integer      :: NGRIDS,FNVAR,NBFS,NE,NG,J,I
      CHARACTER*80 :: tempStrings
      CHARACTER*60 :: FV_NAME,strID1,strID2
      real(4)      :: tempReal,tempReal2,tempReal3
      integer      :: tempInt,boundID,boundNodes
      integer      :: nvrtx_tmp,ncell_tmp,readStat
      integer      :: iNV
      integer      :: iFS,NBFS_tmp=0
      integer      :: iNE
      INTEGER,allocatable :: IFFACE_tmp(:,:,:)
!
!      integer, parameter :: maxdmn=100
      integer            :: idum,idum1
      character*10, allocatable :: dmnnm(:)
      integer, allocatable      :: dmn1(:),dmn2(:),dmn0(:),dmn00(:)
      character(len=1)   :: chgNum='n'
      logical            :: lchgNum=.false.
!
! --- 
!
      iNV=0
      iFS=0
      iNE=0
      nvrtx_tmp=0
      ncell_tmp=0
      allocate(dmn1(maxdmn),dmn2(maxdmn),dmn0(maxdmn),dmn00(maxdmn))
      allocate(dmnnm(maxdmn))
      dmn0=0
      dmn00=0
      dmn1=0
      dmn2=0
      dmnnm=' '

      WRITE(ifll,'(1x,a)') 'MSG: Reading FV File...'
!
      READ(IUTFV,'(A60)') FV_NAME
      WRITE(ifll,'(1x,2a)') 'MSG: Reading field: ',trim(FV_NAME)
      IF(FV_NAME(1:9) .ne. 'FIELDVIEW') THEN
	WRITE(ifle,'(1x,a)') 
     &    'MSG: This is NOT a FIELDVIEW ASCII file'
	ierrorFVG=1
	return
      END IF
!
      READ(IUTFV,2000) tempStrings ! comment
      READ(IUTFV,2000) tempStrings ! comment
      READ(IUTFV,2000) tempStrings
      read(tempStrings,*) strID1, NGRIDS
      WRITE(ifll,'(1x,2a)') 'MSG: Reading field: ',trim(strID1)
      IF (trim(strID1).ne.'Grids') THEN
        WRITE(ifle,'(1x,2a)') 'MSG: ',
     & 'Section header exception happened at Grids'
	ierrorFVG=1
	return
      END IF
!
      IF (NGRIDS.LT.1) THEN
        WRITE(ifle,'(1x,a)') 'MSG: wrong num of grids'
	ierrorFVG=1
	return
      END IF
!
      READ(IUTFV,2000) tempStrings
      read(tempStrings,*) strID1, strID2, NBOUND
      WRITE(IFLL,'(1x,4a)') 
     &            'MSG: Reading field: ',trim(strID1),' ',trim(strID2) 
      IF(trim(strID1).ne.'Boundary'.or.trim(strID2).ne.'Table') THEN
        WRITE(ifle,'(1x,a)') 
     &  'MSG: Section header exception happened at Boundary Table'
	ierrorFVG=1
	return
      END IF
!
 2010 FORMAT(I2,I2,I2,A80)
      allocate(SFBOUN(0:NBOUND))
      allocate(NFBOUN(NBOUND))
      NFBOUN(:)=0
      DO mb=1,NBOUND
      READ(IUTFV,2010) tempInt,tempInt,tempInt, SFBOUN(mb)
      WRITE(IFLL,'(1x,2a,a)') 'MSG: ',
     &              'Reading boundary: ', trim(SFBOUN(mb))
      end do

!
! --- loop with NGRIDS
!
      do NG=1,NGRIDS
!
! --- 'Nodes'
!
      READ(IUTFV,2000) tempStrings
      read(tempStrings,*) strID1, tempInt
      WRITE(IFLL,'(1x,2a)') 'MSG: Reading field: ',trim(strID1)
      IF (trim(strID1).ne.'Nodes') THEN
	WRITE(ifle,'(1x,2a)') 'MSG: ',
     &                ' Section header exception happened at Nodes'
	ierrorFVG=1
	return
      END IF
!
      iNV=nvrtx_tmp
      nvrtx_tmp=nvrtx_tmp+tempInt
      if(nvrtx_tmp>mvrtx) then
         write(ifll,*) 
         write(ifle,*) 'ERR: mvrtx > ',nvrtx_tmp,' in ',cntlnam
         stop
      endif
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
      WRITE(IFLL,'(1x,4a)') 
     &            'MSG: Reading field: ',trim(strID1),' ',trim(strID2)
      IF (trim(strID1).ne.'Boundary' .or. trim(strID2).ne.'Faces') THEN
	WRITE(ifle,'(1x,2a)') 'MSG: ',
     &     'Section header exception happened at Boundary Faces'
	ierrorFVG=1
	return
      ENDIF
!
      if(NG==1) then
         iFS=0
         NBFS=NBFS_tmp
         allocate(IFFACE(4,NBFS,NBOUND))
         IFFACE=0
      else
         iFS=NBFS
         allocate(IFFACE_tmp(4,NBFS,NBOUND))
         IFFACE_tmp(1:4,1:NBFS,1:NBOUND)=IFFACE(1:4,1:NBFS,1:NBOUND)
         NBFS=NBFS+NBFS_tmp
         deallocate(IFFACE)
         allocate(IFFACE(4,NBFS,NBOUND))
         IFFACE=0
         IFFACE(1:4,1:iFS,1:NBOUND)=IFFACE_tmp(1:4,1:iFS,1:NBOUND)
         deallocate(IFFACE_tmp)
      end if
      DO j=iFS+1,NBFS
	READ(IUTFV,2000)  tempStrings
	READ(tempStrings, *) boundID,boundNodes
	NFBOUN(boundID)=NFBOUN(boundID)+1
	READ(tempStrings,*) boundID,boundNodes,
     & (IFFACE(i,NFBOUN(boundID),boundID),i=1,boundNodes)
!       shift vertex No.
        do i=1,4
        if(IFFACE(i,NFBOUN(boundID),boundID)/=0) then
           IFFACE(i,NFBOUN(boundID),boundID)=
     &     IFFACE(i,NFBOUN(boundID),boundID)+iNV
        end if
        end do
      end do
!
! --- 'Elements'
!
      
      READ(IUTFV,2000) tempStrings
      WRITE(IFLL,'(1x,2a)') 'MSG: Reading field: ',trim(tempStrings)
      IF (trim(tempStrings(1:8)).ne.'Elements') THEN
        WRITE(ifle,'(1x,a)') 'MSG: subroutine readFVGrids_30'
	WRITE(ifle,'(1x,a)')
     &  'Section header exception happened at Elements'
	ierrorFVG=1
	return
      END IF
!
      iNE=ncell_tmp
      dmn0(NG)=NG
      dmn1(NG)=iNE+1
      DO
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
	if(iNE.gt.mcell) then
	  WRITE(ifle,'(1x,2a)')
     &    'ERR: Element Number is larger than mcell in ',trim(cntlnam)
          WRITE(ifle,'(3a,2x,I12)') 
     &     'Increas mcell in ',cntlnam,' > ',iNE
	  ierrorFVG=1
	  return
	end if
!
	READ(tempStrings, *) I
      if(multiMaterialFlag)then
	IF(i == 1) THEN
! --- Tet mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(4,iNE), lvcell(1,iNE), lvcell(2,iNE), lvcell(3,iNE),
     &     lacell(iNE)
!       shift vertex No.
        lvcell(1:4,iNE)=lvcell(1:4,iNE)+iNV
	ELSEIF(i == 2) THEN
! --- Hex mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(2,iNE), lvcell(4,iNE), lvcell(3,iNE),
     &     lvcell(5,iNE), lvcell(6,iNE), lvcell(8,iNE), lvcell(7,iNE),
     &     lacell(iNE)
!       shift vertex No.
        lvcell(1:8,iNE)=lvcell(1:8,iNE)+iNV
	ELSEIF(i == 3) THEN
! --- Prism mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(4,iNE), lvcell(5,iNE), lvcell(2,iNE),
     &     lvcell(6,iNE), lvcell(3,iNE),
     &     lacell(iNE)
!       shift vertex No.
        lvcell(1:6,iNE)=lvcell(1:6,iNE)+iNV
	ELSEIF(i == 4) THEN
! --- Pyramid mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(2,iNE), lvcell(3,iNE), lvcell(4,iNE),
     &     lvcell(5,iNE),
     &     lacell(iNE)
!       shift vertex No.
        lvcell(1:5,iNE)=lvcell(1:5,iNE)+iNV
	END IF
      else
	IF(i == 1) THEN
! --- Tet mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(4,iNE), lvcell(1,iNE), lvcell(2,iNE), lvcell(3,iNE)
!	lacell(iNE)=NG
!       shift vertex No.
        lvcell(1:4,iNE)=lvcell(1:4,iNE)+iNV
	ELSEIF(i == 2) THEN
! --- Hex mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(2,iNE), lvcell(4,iNE), lvcell(3,iNE),
     &     lvcell(5,iNE), lvcell(6,iNE), lvcell(8,iNE), lvcell(7,iNE)
!	lacell(iNE)=NG
!       shift vertex No.
        lvcell(1:8,iNE)=lvcell(1:8,iNE)+iNV
	ELSEIF(i == 3) THEN
! --- Prism mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(4,iNE), lvcell(5,iNE), lvcell(2,iNE),
     &     lvcell(6,iNE), lvcell(3,iNE)
!	lacell(iNE)=NG
!       shift vertex No.
        lvcell(1:6,iNE)=lvcell(1:6,iNE)+iNV
	ELSEIF(i == 4) THEN
! --- Pyramid mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,iNE), lvcell(2,iNE), lvcell(3,iNE), lvcell(4,iNE),
     &     lvcell(5,iNE)
!	lacell(iNE)=NG
!       shift vertex No.
        lvcell(1:5,iNE)=lvcell(1:5,iNE)+iNV
	END IF
      endif
      end do
!
! --- end loop with NGRIDS
!
      end do !do ng=1,NGRIDS

!
! --- check 'ncell'
!
      if(nvrtx_tmp.ne.nvrtx) then
	WRITE(ifle,'(1x,a,2I12)') 
     &         'Number of vertex is not matched.',nvrtx_tmp,nvrtx
        WRITE(ifle,'(1x,a,I12,2a)') 
     &     'nvrtx= ',nvrtx_tmp,' => Reset [nvrtx] in ',cntlnam
	ierrorFVG=1
	return
      end if


      if(iNE.ne.ncell) then
	WRITE(ifle,*)
     &  'MSG: ',
     &  'Number of cell is not matched.',iNE,ncell
        WRITE(ifle,*) 'ncell= ',iNE,' => Reset [ncell] in ',cntlnam
	ierrorFVG=1
	return
      endif
!
!
! --- define material number
!
      if(IUTRG/=-1) then
      READ(IUTRG,2000) tempStrings ! comment
      READ(IUTRG,2000) tempStrings ! comment
      READ(IUTRG,'(A60)') FV_NAME
      WRITE(ifll,'(1x,2a)') 'MSG: Reading field: ',trim(FV_NAME)
      IF(FV_NAME(1:5) .ne. 'FVREG') THEN
	WRITE(ifle,'(1x,a)') 
     &    'MSG: This is NOT a FIELDVIEW REGION file'
        dmnnm='fluid'
!	ierrorFVG=1
!	return
      ELSE
        READ(IUTRG,2000) tempStrings ! DATASET_COORD_TYPE
        do ng=1,NGRIDS
           READ(IUTRG,2000) tempStrings
           IF (trim(tempStrings).ne.'REGION') THEN
              WRITE(ifle,'(1x,2a)') 'MSG: ',
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

      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
      do i=1,maxdmn
        idum=dmn0(i)   !idum=> region ID number
        if(idum/=0) then
          write(ifll,'(1X,a,I4,a,a10,a,I12,a,I12)') 
     &    'MSG : Material no.: IMAT_U(=lacell)=',i,' Region Name:['
     &    ,trim(dmnnm(idum)),'] cell# from',dmn1(i),' to',dmn2(i)
        endif
      enddo
      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
!
      write(ifll,'(1x,a)') '-'
      write(ifll,'(1x,a)') 
     & 'Attention: Input Material No. for multi-matrial Ver. of Gridgen'
      write(ifll,'(1x,a,4X,a)')
     &      'MSG: You should Change Material number !',
     &      ' Will you change Material number ? (n/y)'
      write(ifll,'(1x,a)')
     &'MSG: FrontFlowRed support + value for fluid Material no.[IMAT_U]'
      write(ifll,'(1x,a)')
     &'MSG: FrontFlowRed support - value for solid Material no.[IMAT_U]'
      
      read(*,*) chgNum
      lchgNum=.false.
      if(chgNum=='y'.or.chgNum=='Y') lchgNum=.true.
      do i=1,maxdmn
        idum1=dmn0(i)   !idum=> region ID number
        if(idum1/=0) then
          if(lchgNum) then
            write(ifll,'(1x,a,I4,a,a,a)') 
     &     'MSG: Mesher defined Material Number : ',i,
     &     ' Region Name:[',trim(dmnnm(idum1)),']'
            write(ifll,'(1x,a)') 
     &      'MSG: New number for this material: '
            read(*,*) idum1
            IF(idum1==0) STOP 'ERR: no-zero needed'
            dmn00(i)=idum1
          else
            idum1=i
            dmn00(i)=idum1
          endif
!          dmn00(i)=idum1
          do k=dmn1(i),dmn2(i)
            lacell(k)=idum1   !dmn0(i)
          end do
          write(ifll,'(1X,a,I4,a,a10,a,I12,a,I12)') 
     &    'MSG : IMAT_U(=lacell)=',idum1,' Region Name:['
     &    ,trim(dmnnm(idum1)),'] cell# from',dmn1(i),' to',dmn2(i)
          write(ifll,'(1X)')
          write(ifll,'(1X)')
        endif
      enddo
!
      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
      do i=1,maxdmn
        idum=dmn0(i)   !idum=> region ID number
        if(idum/=0) then
          write(ifll,'(1X,a,I4,a,a,a,I12,a,I12)') 
     &    "MSG : Material no.: IMAT_U(=lacell)=",
     &     dmn00(i)," Region Name:["
     &    ,trim(dmnnm(idum)),"] cell# from",dmn1(i)," to",dmn2(i)
        endif
      enddo
      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
!
      OPEN(61,FILE='material.fflow',FORM='FORMATTED',status='unknown',
     & iostat=ios)
      write(61,'(1X)')
      write(61,'(1x,a)') '**********************************'
      write(61,'(1X)')
      do i=1,maxdmn
        idum=dmn0(i)   !idum=> region ID number
        if(idum/=0) then
          write(61,'(1X,a,I4,a,a10,a,I12,a,I12)') 
     &    "MSG : Material no.: IMAT_U(=lacell)=",dmn00(i),
     &   " Region Name:["
     &    ,trim(dmnnm(idum)),"] cell# from",dmn1(i)," to",dmn2(i)
        endif
      enddo
      write(61,'(1X)')
      write(61,'(1x,a)') '**********************************'
      write(61,'(1X)')
      close(61)
!
      deallocate(dmn1,dmn2,dmn0,dmn00)
      deallocate(dmnnm)
!
      WRITE(ifll,'(1x,a)') '!!! FINISH READING  !!!.'
2000  FORMAT(A80)
      ierrorFVG=0
!
      end subroutine readFVGrids_30

!==================================================================
      subroutine readFVGrids(IUTFV,ierrorFVG)
!==================================================================
!     This subroutine read 'FieldView' grid format 
!     which is created by 'Gridgen'
!==================================================================
!
! --- [module arguments]
!
      integer,intent(in)     :: IUTFV
      integer,intent(inout)  :: ierrorFVG
!
! --- [local entities]
!
      integer      :: NGRIDS,FNVAR,NBFS,NE,J,I
      CHARACTER*80 :: tempStrings
      CHARACTER*60 :: FV_NAME
      real(4)      :: tempReal,tempReal2,tempReal3
      integer      :: tempInt,boundID,boundNodes
!
! --- 
!
      WRITE(ifll,'(1x,a)') 'MSG: Reading FV File...'
!
      READ(IUTFV,'(A60)') FV_NAME
      WRITE(ifll,'(1x,2a)') 'MSG: Reading field: ',trim(FV_NAME)
      IF(FV_NAME(1:9) .ne. 'FIELDVIEW') THEN
	WRITE(ifle,'(1x,a)') 
     &    'MSG: This is NOT a FIELDVIEW ASCII file'
	ierrorFVG=1
	return
      END IF
!
      READ(IUTFV,2000) tempStrings
      WRITE(ifll,'(1x,2a)') 
     &   'MSG: Reading field: ',trim(tempStrings)
      IF (tempStrings(:len_trim(tempStrings)).ne. 'CONSTANTS') THEN
        WRITE(ifle,'(1x,2a)') 'MSG: ', 
     &       'Section header exception happened at CONSTANTS'
        WRITE(ifle,'(1x,2a)') 
     &  '       1.miss type of "gdformat" in ',cntlnam
        WRITE(ifle,'(1x,2a)') 
     &   '       2.wrong charactor codes(CR<->CR+LF etc.',
     &   ' wrong FTP mode(ASCII/Binary) causes this error)'
        WRITE(ifle,'(1x,a)') '       3.wrong file, etc.'
	ierrorFVG=1
	return
      END IF
      READ(IUTFV,*) tempReal	!	TIME
      READ(IUTFV,*) tempReal	!	FSMACH
      READ(IUTFV,*) tempReal	!	ALPHA
      READ(IUTFV,*) tempReal	!	RE


      READ(IUTFV,2000) tempStrings
      WRITE(ifll,'(1x,2a)') 'MSG: Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'GRIDS') THEN
        WRITE(ifle,'(1x,2a)') 'MSG: ',
     & 'Section header exception happened at GRIDS'
	ierrorFVG=1
	return
      END IF
!
      READ(IUTFV,*) NGRIDS
      IF (NGRIDS.GT.1) THEN
        WRITE(ifle,'(1x,a)') 'MSG: Multi grid file is not supported'
	ierrorFVG=1
	return
      END IF
!
      READ(IUTFV,2000) tempStrings
      WRITE(IFLL,'(1x,2a)') 'MSG: Reading field: ',trim(tempStrings)
      IF(tempStrings.ne.'Boundary Table') THEN
        WRITE(ifle,'(1x,a)') 
     &  'MSG: Section header exception happened at Boundary Table'
	ierrorFVG=1
	return
      END IF
!
      READ(IUTFV,*) NBOUND
 2010 FORMAT(I1,A80)
      allocate(SFBOUN(0:NBOUND))
      DO mb=1,NBOUND
      READ(IUTFV,2010) tempInt, SFBOUN(mb)
      WRITE(IFLL,'(1x,2a,I8,a)') 'MSG: ',
     &              'Reading boundary: ', tempInt, SFBOUN(mb)
      end do
!
      READ(IUTFV,2000) tempStrings
      WRITE(IFLL,'(1x,2a)') 'MSG: Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'Variable Names') THEN
	WRITE(ifle,'(1x,a)')  
     &   'MSG: Section header exception happened at Variable Names'
	ierrorFVG=1
	return
      ENDIF
!
      READ(IUTFV,*) FNVAR
      IF (FNVAR.GE.1)
     &  WRITE(ifle,'(1x,2a)') 'MSG: ',
     &           'Sorry, this version does not support variables'
      DO j = 1 ,FNVAR	
	READ(IUTFV,2000) tempStrings	
	WRITE(ifle,'(1x,3a)') 'MSG: ',
     &                'WARRNING! IGNORE :',trim(tempStrings)
      end do
!
! --- 'Nodes'
!
      READ(IUTFV,2000) tempStrings
      WRITE(IFLL,'(1x,2a)') 'MSG: Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'Nodes') THEN
	WRITE(ifle,'(1x,2a)') 'MSG: ',
     &                ' Section header exception happened at Nodes'
	ierrorFVG=1
	return
      END IF
!
      READ(IUTFV,*) tempInt
      if(tempInt.ne.nvrtx) then
	WRITE(ifle,'(1x,a,2I12)') 
     &         'Number of vertex is not matched.',tempInt,nvrtx
        WRITE(ifle,'(1x,a,I12,2a)') 
     &     'nvrtx= ',tempInt,' => Reset [nvrtx] in ',cntlnam
	ierrorFVG=1
	return
      end if
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
      WRITE(IFLL,'(1x,2a)') 
     &              'MSG: Reading field: ',trim(tempStrings)
      IF (tempStrings.ne.'Boundary Faces') THEN
	WRITE(ifle,'(1x,2a)') 'MSG: ',
     &     'Section header exception happened at Boundary Faces'
	ierrorFVG=1
	return
      ENDIF
!
      READ(IUTFV,*) NBFS
      allocate(NFBOUN(NBOUND))
      allocate(IFFACE(4,NBFS,NBOUND))
      NFBOUN(:)=0
      IFFACE=0
      DO j=1,NBFS
	READ(IUTFV,2000)  tempStrings
	READ(tempStrings, *) boundID,boundNodes
	NFBOUN(boundID)=NFBOUN(boundID)+1
	READ(tempStrings,*) boundID,boundNodes,
     & (IFFACE(i,NFBOUN(boundID),boundID),i=1,boundNodes)
      end do
!
! --- 'Elements'
!
      READ(IUTFV,2000) tempStrings
      WRITE(IFLL,'(1x,2a)') 'MSG: Reading field: ',tempStrings
      IF (tempStrings.ne.'Elements') THEN
        WRITE(ifle,'(1x,a)') 'MSG: subroutine readFVGrids'
	WRITE(ifle,'(1x,a)')
     &  'Section header exception happened at Elements'
	ierrorFVG=1
	return
      END IF
!
      NE=0
      DO
        READ(IUTFV,2000) tempStrings
        IF (tempStrings .eq.'Variables') then
            write(ifll,'(1x,a)')  trim(tempStrings)
            exit
        endif
	NE=NE+1
	if(NE.gt.mcell) then
	  WRITE(ifle,'(1x,2a)')
     &    'ERR: Element Number is larger than mcell in ',trim(cntlnam)
          WRITE(ifle,'(3a,2x,I12)') 
     &     'Increas mcell in ',cntlnam,' > ',NE
	  ierrorFVG=1
	  return
	end if
!
	READ(tempStrings, *) I
      if(multiMaterialFlag)then
	IF(i == 1) THEN
! --- Tet mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(4,NE), lvcell(1,NE), lvcell(2,NE), lvcell(3,NE),
     &     lacell(NE)
	ELSEIF(i == 2) THEN
! --- Hex mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,NE), lvcell(2,NE), lvcell(4,NE), lvcell(3,NE),
     &     lvcell(5,NE), lvcell(6,NE), lvcell(8,NE), lvcell(7,NE),
     &     lacell(NE)
	ELSEIF(i == 3) THEN
! --- Prism mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,NE), lvcell(4,NE), lvcell(5,NE), lvcell(2,NE),
     &     lvcell(6,NE), lvcell(3,NE),
     &     lacell(NE)
	ELSEIF(i == 4) THEN
! --- Pyramid mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,NE), lvcell(2,NE), lvcell(3,NE), lvcell(4,NE),
     &     lvcell(5,NE),
     &     lacell(NE)
	END IF
      else
	IF(i == 1) THEN
! --- Tet mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(4,NE), lvcell(1,NE), lvcell(2,NE), lvcell(3,NE)
	lacell(NE)=1
	ELSEIF(i == 2) THEN
! --- Hex mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,NE), lvcell(2,NE), lvcell(4,NE), lvcell(3,NE),
     &     lvcell(5,NE), lvcell(6,NE), lvcell(8,NE), lvcell(7,NE)
	lacell(NE)=1
	ELSEIF(i == 3) THEN
! --- Prism mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,NE), lvcell(4,NE), lvcell(5,NE), lvcell(2,NE),
     &     lvcell(6,NE), lvcell(3,NE)
	lacell(NE)=1
	ELSEIF(i == 4) THEN
! --- Pyramid mesh
        READ(tempStrings, *) i, tempInt,
     &     lvcell(1,NE), lvcell(2,NE), lvcell(3,NE), lvcell(4,NE),
     &     lvcell(5,NE)
	lacell(NE)=1
	END IF
      endif
      end do
!
! --- check 'ncell'
!
      if(NE.ne.ncell) then
	WRITE(ifle,*)
     &  'MSG: ',
     &  'Number of cell is not matched.',NE,ncell
        WRITE(ifle,*) 'ncell= ',NE,' => Reset [ncell] in ',cntlnam
	ierrorFVG=1
	return
      endif
!
!
      WRITE(ifll,'(1x,a)') '!!! FINISH READING  !!!.'
2000  FORMAT(A80)
      ierrorFVG=0
!
      end subroutine readFVGrids
!
!
!#ifdef STARCD

!
! DEBUG
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine readGAMBITgrids(
     &     mcell,mface,mvrtx,ncell,nvrtx,cord,lacell,
     &     lvcell,lfcell,lcface,lvface,
     &     IUTFV,ierrorGB)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
! --- [dummy arguments]
!
      integer, intent(in)       :: IUTFV
      integer, intent(inout)    :: ierrorGB
!
      integer,intent(in)    :: mcell,mface,mvrtx,ncell
      integer,intent(inout) :: nvrtx
      real*8 ,intent(out)   :: cord(3,mvrtx)
      integer,intent(out)   :: lacell(  mcell)
      integer,intent(out)   :: lvcell(8,mcell)
      integer,intent(inout) :: lfcell(7,mcell)  
      integer,intent(inout) :: lcface(2,mface)  
      integer,intent(inout) :: lvface(4,mface)        
!
! --- [local entities]
!
      integer, allocatable      :: mfboun(:),nfpz(:)
      character*80, allocatable :: nmboun(:)
      character*10, allocatable :: dmnnm(:)
      integer, allocatable      :: dmn1(:),dmn2(:),dmn0(:),dmn00(:)
      integer, allocatable      :: kmesh(:)
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
      integer  :: NBFS,IMAT,IMAT1,NFCE,IS1,IS2,IS3,IC1,IC2,IS
      
      integer  :: IV1,IV2,IV3,IV4,IV5,IV6,IV,IIS,IERR
      integer  :: pnface,psface,nf,nc,nbnd,bndsfCnt
      integer  :: tmpI1,tmpI2,tmpI3,tmpI4,tmpI5,tmpI6,tmpI7
      integer  :: nnf4(2),nnf6(6),nnf(6),IVN0(4)
      integer  :: ICNUM(2),IVNUM(4),IVN1(4),IVN2(4),IVN3(4),IVN4(4)
      real(4)  :: tempReal,tempReal2,tempReal3,aa
      integer  :: idum,idum1
      integer   :: strptr
      character :: strID
      logical   :: gridgenflag
!      integer, parameter :: maxboun=100, maxdmn=100
!      integer, parameter :: maxboun=40, maxdmn=30
      integer   :: nvrtx_start,nvrtx_end,nvrtxx
      character(len=1) :: chgNum='n'
      logical :: lchgNum=.false.
!
      nvrtx_start=0
      nvrtx_end=0
!
      nbnd=0
      lacell=0
      lvcell(:,:)=0

!
! --- ----------------------------------------------------
      write(ifll,*) 'MSG: Reading GAMBIT grid file ...'
! --- ----------------------------------------------------
      boundID=0
      allocate(nmboun(1:maxboun)) ; nmboun=' '
      allocate(mfboun(maxboun))   ; mfboun=0
      allocate(nfpz(maxboun))     ; nfpz=0
!
      allocate(dmn1(maxdmn),dmn2(maxdmn),dmn0(maxdmn),dmn00(maxdmn))
      allocate(dmnnm(maxdmn))
!
      dmn0=0
      dmn00=0
      dmn1=0
      dmn2=0
      dmnnm=' '
      IMAT=0
      IMAT1=0
!
!     ###############################################
!     ##### read data regarding BC information ######
!     ###############################################
!
      read(IUTFV,'(A60)') FV_NAME
      write(ifll,*) 'MSG: First LINE in the FILE :',trim(FV_NAME)
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
      lvface=0
      lcface=0
      lfcell=0
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
              write(ifll,*) 'MSG: FACE number: NFCE=',NFCE
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
            end if
          elseif(tmpI1==12) then      ! CELLS block
            read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc10
            if(tmpc6(2:2)=='0') then
              read(tmpc10,'(z8)') tmpI1
              if(ncell/=tmpI1) then
                write(ifll,*) 
     &                'ERR: NCELL should be :',tmpI1,' in ',cntlnam
                stop
              endif
              write(ifll,*) 'MSG: CELL number: ncell=',tmpI1
! allocate
              allocate(kmesh(1:tmpI1),stat=ier)
              if(ier/=0) stop 'stop at allocating in read_griddata(GB)'
              kmesh=0

            endif
            if(tmpc6(2:2)/='0') then
!             2 more string in 12) section, in "regular cell section"
              IMAT=IMAT+1
              read(tempStrings,*) 
     &                       tmpc5,tmpc6,tmpc9,tmpc10,tmpc11,tmpc12
              read(tmpc6(2:3),'(z8)') tmpI1      ! Zone-ID
              read(tmpc9,'(z8)')  tmpI2
              read(tmpc10,'(z8)') tmpI3
              clen=len_trim(adjustl(tmpc12))
              do i=1,clen
                 strptr=scan(tmpc12,'()')
                 if(strptr>=1.and.strptr<=clen) then
                    tmpc12=tmpc12(1:strptr-1)//tmpc12(strptr+1:clen)
                 end if
              enddo
              read(tmpc12,'(z8)') tmpI4
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
                 do ! loop
                    if(i>tmpI3) exit
                    read(IUTFV,'(a)',ADVANCE='NO',EOR=101) tmpc7
 101                read(IUTFV,'(a)',ADVANCE='NO') tmpc7
                    do while(tmpc7==' ' .or. tmpc7=='(')
                       read(IUTFV,'(a)',ADVANCE='NO') tmpc7
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
              dmn0(IMAT)=tmpI1  !IMAT
              dmn1(IMAT)=tmpI2
              dmn2(IMAT)=tmpI3
              
            endif
          elseif(tmpI1==10) then
            read(tempStrings,*) tmpc5,tmpc6,tmpc9,tmpc10
            if(tmpc6(2:2)/='0') then
              read(tmpc9, '(z8)') tmpI1
              read(tmpc10,'(z8)') tmpI2
              do i=tmpI1,tmpI2
              nvrtxx=nvrtxx+1
              read(IUTFV,*) tempReal,tempReal2,tempReal3
              cord(1,i)=dble(tempReal )
              cord(2,i)=dble(tempReal2)
              cord(3,i)=dble(tempReal3)
              enddo
            endif
          elseif(tmpI1==39) then         !(maybe) boundary condition
            write(ifll,*) 'NOT YET SUPPORTED for Index 39'
          elseif(tmpI1==45) then         !boundary condition
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
          endif
        endif
      enddo
! --- --------------------------------------------------------------
      if(nvrtx/=nvrtxx) then
         write(ifll,*) 
     &        'WRN: NVRTX should be :',nvrtxx,' in ',cntlnam
         nvrtx=nvrtxx
      endif
!
      if(nvrtx>mvrtx) then
         write(ifll,*) 
         write(ifle,*) 'ERR: mvrtx > ',nvrtx,' in ',cntlnam
         stop
      endif

      if(nf.ne.NFCE) then
        write(ifle,*) 'ERR: reading error in face number',nf,NFCE
        stop
      endif
!
      if(nf>mface) then
        write(ifle,*) 'ERR: mface > ',nf,' in ',cntlnam
        stop ' read_grid_data '
      endif
!
      write(ifll,*) 'MSG: Total face number=',nf
      write(ifll,*) 'MSG: BC number without INTERIOR and FLUID=',nbnd
      if(dum_bc>0) then
        nbnd=nbnd+1
      endif
!
      isum=0
      do izn=1,maxboun
        if(nmboun(izn)/=' ') then
          isum=isum+nfpz(izn)
        endif
      enddo
      
      bndsfCnt=isum
      NBFS=isum
      write(ifll,*) 'MSG: BC FACES number= ',bndsfCnt
!
! --- --------------------------------------------------------------
!
      allocate(SFBOUN(0:nbnd))
      allocate(NFBOUN(nbnd))
      allocate(IFFACE(4,bndsfCnt,nbnd))
      NFBOUN=0
      IFFACE=0
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
      write(ifll,*) 'MSG: BC NAME :'
      do i=1,nbnd
      write(ifll,*) 'MSG: BC FACE',i,trim(SFBOUN(i)),NFBOUN(i)
      end do
!
      deallocate(nfpz)
!
!    ##################################################
!    ###### read data again for node coordinates ######
!    ##################################################
!
      boundID=0
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
                nf=0
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
                   IFFACE(1,nf,boundID)= tmpI2
                   IFFACE(2,nf,boundID)= tmpI3
                   IFFACE(3,nf,boundID)= tmpI4
                   IFFACE(4,nf,boundID)= 0
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
                   IFFACE(1,nf,boundID)= tmpI2
                   IFFACE(2,nf,boundID)= tmpI3
                   IFFACE(3,nf,boundID)= tmpI4
                   IFFACE(4,nf,boundID)= tmpI5
                end if
                end do
              elseif(elem_ty=='2') then
                write(ifll,*) '<Read element-type: linear>'
                write(ifll,*) 'NOT YET SUPPORTED'
              elseif(elem_ty=='3') then
                nf=0
                do i=psface,pnface
                nf=nf+1
                read(IUTFV,'(A60)',iostat=ier) tempStrings
                read(tempStrings,*) tmpc6,tmpc9,tmpc10,tmpc11,tmpc12
                read(tmpc6, '(z10)') tmpI2
                read(tmpc9, '(z10)') tmpI3
                read(tmpc10,'(z10)') tmpI4
                read(tmpc11,'(z10)') tmpI5
                read(tmpc12,'(z10)') tmpI6
                IFFACE(1,nf,boundID)= tmpI2
                IFFACE(2,nf,boundID)= tmpI3
                IFFACE(3,nf,boundID)= tmpI4
                IFFACE(4,nf,boundID)= 0
                enddo
                write(ifll,*) 
     &              '<Read element-type: triangular> for IFFACE'
              elseif(elem_ty=='4') then
                nf=0
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
                IFFACE(1,nf,boundID)=tmpI1
                IFFACE(2,nf,boundID)=tmpI2
                IFFACE(3,nf,boundID)=tmpI3
                IFFACE(4,nf,boundID)=tmpI4
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
      if(dum_bc>0) then
        NBOUND=boundID+1
      endif
      write(ifll,*) 'MSG : Face Number in BC:',nf,NFCE
      write(ifll,*) 'MSG : complete reading all the faces'
!
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
      if(dum_bc>0) then
        SFBOUN(NBOUND)='dummy'
        idum=0
        do izn=1,boundID
        if(trim(boundName(dum_bc,1))==trim(SFBOUN(izn))) then
          idum=izn
          exit
        endif
        enddo
        dum_bc=idum
        if(idum/=0) then
         NFBOUN(NBOUND)=NFBOUN(idum)
         IFFACE(:,1:NFBOUN(NBOUND),NBOUND)=IFFACE(:,1:NFBOUN(idum),idum)
        endif
      endif
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
      write(ifll,*) 'MSG : Completed connectivity '
! +++++++++++++++++++++++++
!
      deallocate(nmboun)
      deallocate(mfboun)
!
! --- define material number
!
      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
      do i=1,maxdmn
        idum=dmn0(i)   !idum=> region ID number
        if(idum/=0) then
          write(ifll,'(1X,a,I4,a,a10,a,I12,a,I12)') 
     &    'MSG : Material no.: IMAT_U(=lacell)=',i,' Region Name:['
     &    ,trim(dmnnm(idum)),'] cell# from',dmn1(i),' to',dmn2(i)
        endif
      enddo
      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
!
      write(ifll,'(1x,a)') '-'
      write(ifll,'(1x,a,4X,a)')
     &      'MSG: You Can Change Material number !',
     &      ' Will you change Material number ? (n/y)'
      write(ifll,'(1x,a)')
     &'MSG: FrontFlowRed support + value for fluid Material no.[IMAT_U]'
      write(ifll,'(1x,a)')
     &'MSG: FrontFlowRed support - value for solid Material no.[IMAT_U]'
      
      read(*,*) chgNum
      lchgNum=.false.
      if(chgNum=='y'.or.chgNum=='Y') lchgNum=.true.
      do i=1,maxdmn
        idum=dmn0(i)   !idum=> region ID number
        if(idum/=0) then
          if(lchgNum) then
            write(ifll,'(1x,a,I4,a,a,a)') 
     &     'MSG: Mesher defined Material Number : ',i,
     &     ' Region Name:[',trim(dmnnm(idum)),']'
            write(ifll,'(1x,a)') 
     &      'MSG: New number for this material: '
            read(*,*) idum1
            IF(idum1==0) STOP 'ERR: no-zero needed'
            dmn00(i)=idum1
          else
            idum1=i
            dmn00(i)=idum1
          endif
          do k=dmn1(i),dmn2(i)
            lacell(k)=idum1   !dmn0(i)
          end do
          write(ifll,'(1X,a,I4,a,a10,a,I12,a,I12)') 
     &    'MSG : IMAT_U(=lacell)=',idum1,' Region Name:['
     &    ,trim(dmnnm(idum)),'] cell# from',dmn1(i),' to',dmn2(i)
          write(ifll,'(1X)')
          write(ifll,'(1X)')
        endif
      enddo
!
      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
      do i=1,maxdmn
        idum=dmn0(i)   !idum=> region ID number
        if(idum/=0) then
          write(ifll,'(1X,a,I4,a,a,a,I12,a,I12)') 
     &    "MSG : Material no.: IMAT_U(=lacell)=",
     &     dmn00(i)," Region Name:["
     &    ,trim(dmnnm(idum)),"] cell# from",dmn1(i)," to",dmn2(i)
        endif
      enddo
      write(ifll,'(1X)')
      write(ifll,'(1x,a)') '**********************************'
      write(ifll,'(1X)')
!
      OPEN(61,FILE='material.fflow',FORM='FORMATTED',status='unknown',
     & iostat=ios)
      write(61,'(1X)')
      write(61,'(1x,a)') '**********************************'
      write(61,'(1X)')
      do i=1,maxdmn
        idum=dmn0(i)   !idum=> region ID number
        if(idum/=0) then
          write(61,'(1X,a,I4,a,a10,a,I12,a,I12)') 
     &    "MSG : Material no.: IMAT_U(=lacell)=",dmn00(i),
     &   " Region Name:["
     &    ,trim(dmnnm(idum)),"] cell# from",dmn1(i)," to",dmn2(i)
        endif
      enddo
      write(61,'(1X)')
      write(61,'(1x,a)') '**********************************'
      write(61,'(1X)')
      close(61)
!
      lfcell(:,:)=0
      lcface(:,:)=0
      lvface(:,:)=0
      deallocate(dmnnm)
!
      return
      end subroutine readGAMBITgrids

!#endif
!
!#ifdef SCRYU
!
! ***************************************************************
      subroutine readSCRYUgrids(IUTFV,ierrorGG)
! ***************************************************************
      implicit none
      integer, intent(in)  :: IUTFV
      integer, intent(inout) :: ierrorGG
!
      integer  :: i,j,nn,boundID,ierr
      character*12  :: LREGN
      integer  :: NBFS
      integer  :: L2D3D,NNODS,NELEM,NREGN,NREGI,IDVER,LFORT,NGR
     &           ,NDE,NE
      integer  :: nn0,nn1,nn2,nn3,nn4
      integer  :: Nsum,Fsum
!
      integer, allocatable :: LGR(:),MAT(:),NDNO(:),IE(:),IFA(:)
      integer, allocatable :: IEtot(:,:),IFAtot(:,:),nnv(:)
      integer, allocatable :: MFBOUN(:)
      real(4),allocatable  :: xcrd(:),ycrd(:),zcrd(:)
      character*1, allocatable :: IETYP(:)
      character*12, allocatable :: NMBOUN(:)
!
      return
      end subroutine readSCRYUgrids
!#endif
! --- -------------------------------------------------------
!==================================================================
      subroutine readNASTRANgrids(IUTFV,ierror)
!==================================================================
!     This subroutine read 'NASTRAN' grid format 
!==================================================================
!
! --- [module arguments]
!
      integer,intent(in)     :: IUTFV
      integer,intent(inout)  :: ierror
!
! --- [local entities]
!
      integer      :: NGRIDS,NELEMT,FNVAR,NBFS,NE,J,I
      CHARACTER*80 :: tempStrings
      CHARACTER*60 :: FV_NAME
      real(8)      :: dum1,dum2,dum3
      integer      :: tempInt,boundID,boundNodes,imat,ino
      integer      :: idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,
     &                idum9
      integer      :: idumBC=0,NdumBC=0,NBC=0
      real(4)      :: tempReal,tempReal2,tempReal3
      CHARACTER*8  :: cdum1
!
!
      end subroutine readNASTRANgrids
!
      end subroutine read_grid_data
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine storeGridData(nvrtx,mvrtx,mcell,ncell,
     &  cord,lacell,ierrorFVG)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only          : ifli,ifll,ifle,cntlnam
      use module_partitioner,only : givbcnd,glvbcnd,pMap,tempFlag
      use module_boundary,only : nbcnd,boundName,numbname,ivbcnd,
     &                           lvbcnd,icbinr,lcbinr,ivpair,
     &                        kdbuff,kdintr,kdbcnd,kdprdc,kdtchi,kdsld,
     &                           prdcAxis,prdcOffset,kdshutr,kdpors,
     &                           boundIDMap,SFBOUN,NBOUND,nobcnd,
     &                           IFFACE,NFBOUN,idis,NAME_I,set_rotmtrx1
      use module_param ,only   : scalfctr,smallNrm
      use module_material,only : nsold,nosld
      use module_material,only : nflud,nofld
!
      implicit none
!
      integer,intent(in)    :: nvrtx,mvrtx,mcell,ncell
      real*8 ,intent(in)    :: cord(3,mvrtx)
      integer,intent(in)    :: lacell(mcell)
      integer,intent(inout) :: ierrorFVG
!
! --- [local entities]
!
      real(8)             :: dum1,dum2,dum3,dumy(3),rot(3,3),dum4,dum(3)
      integer             :: IBIDa,IBIDb,kd
      integer             :: tempInt,K,J,I,M,IV,Int,nb,mb,mmb,IS
      integer             :: nbcndde
      integer,allocatable :: namflg(:,:)
      integer :: MMAT(-100:100)
!
! --- 
!

      MMAT=0
      do J=1,ncell
      I=lacell(J)
      MMAT(I)=1
      enddo
      do 1200 I=-100,100
      if(MMAT(I)==0) cycle
      do J=1,nflud
      if(nofld(J)==I) then
        goto 1200
      endif
      enddo
!
      do J=1,nsold
      if(nosld(J)==I) then
        goto 1200
      endif
      enddo
!
      write(ifle,'(1X,a,I6)') 
     &  'ERR: Cell Material is NOT defined in fflow.ctl',I
      stop 'stop in storeGridData'
 1200 enddo
!
      allocate(boundIDMap(nbcnd,numbname),namflg(nbcnd,numbname))
      boundIDMap=-1
!---------------------------
! --- serach undefined BC
!---------------------------
      nbcndde=nbcnd
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j)=='undefined') then
        nbcndde=nbcnd-1
        boundIDMap(nbcnd,1)=0
        boundIDMap(nbcnd,2)=-1
        SFBOUN(0)='undefined'
        NAME_I=0
      endif
      enddo
      enddo
!
      IF(nbcndde.LT.nbcnd-1) THEN
        WRITE(ifle,'(1x,2a)')
     &   ' ### ERR: NAME [undefined] BC too many in ',cntlnam
        ierrorFVG=1
        return
      ENDIF
!
      do mb=1,NBOUND
         write(ifll,'(1x,3a)') 'MSG: Boundary Condition Name: [',
     &      trim(adjustl(SFBOUN(mb))),'] defined in mesh file'
      enddo
!
      if(NBOUND.lt.nbcndde) then
	WRITE(ifle,'(2a,2I8)') 'MSG: ',
     &   'Number of boundary region is not matched',NBOUND,nbcnd
	ierrorFVG=1
        return
      endif
!-------------------------------------------------------------
! --- CHECK BC NAME between fflow.ctl and BC name in BC-file
!-------------------------------------------------------------
      allocate(tempFlag(NAME_I:NBOUND))
      namflg(:,:)=0
      tempFlag(:)=.false.
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j)/=' ') then
        do mb=NAME_I,NBOUND
        if(trim(adjustl(boundName(nb,j)))==trim(adjustl(SFBOUN(mb))))
     &  then
          tempFlag(mb)=.true.
          namflg(nb,j)=1
        endif
        enddo
        if(trim(adjustl(boundName(nb,j)))=='dummy') then
          namflg(nb,j)=1
        endif
      endif
      enddo
      enddo
!
      do mb=NAME_I,NBOUND
      if(.not.tempFlag(mb)) then
        write(ifle,'(2a,2x,2a)') 'ERR-2: NOT Defined BC name in ',
     &  cntlnam,' : ',trim(adjustl(SFBOUN(mb)))
        stop
      endif
      enddo
!
      do j=1,numbname
      do nb=1,nbcnd
      if(namflg(nb,j)==0.and.boundName(nb,j)/=' ') then
        write(ifle,*) 'ERR: NOT Found BC name in user-BC mesh file: ',
     &  trim(adjustl(boundName(nb,j)))
        stop
      endif
      enddo
      enddo
      deallocate(tempFlag,namflg)
!----------------------
! --- Map BC number
!----------------------
      do j=1,numbname
      do nb=1,nbcnd
      if(boundName(nb,j).ne.' ') then
        do mb=NAME_I,NBOUND
        if(trim(adjustl(boundName(nb,j))).eq.trim(adjustl(SFBOUN(mb))))
     &  then
          boundIDMap(nb,j)=mb
          exit
        endif
	enddo
      endif
      enddo
      enddo
!
      allocate(tempFlag(NAME_I:NBOUND))
      tempFlag=.false.
!
      do j=1,numbname
      do nb=1,nbcnd
      mb=boundIDMap(nb,j)
      if(boundIDMap(nb,j).ne.-1) then
        if((boundIDMap(nb,j)<0).or.(boundIDMap(nb,j)>NBOUND)) then
	  WRITE(ifle,*) 'MSG: ',
     &   'Wrong boundary ID.',boundIDMap(nb,j)
          ierrorFVG=1
	  return
        end if
        if(kdbcnd(0,nb)/=kdtchi) then
          if(tempFlag(mb)) then
	    WRITE(ifle,'(3a)')  
     &      'MSG:  Boundary [',
     &      trim(SFBOUN(boundIDMap(nb,j))),'] is defined twice.'
            ierrorFVG=1
            return
          endif
        endif
        tempFlag(mb)=.true.
      endif
      enddo
      enddo
!
      do mb=NAME_I,NBOUND
        do mmb=mb+1,NBOUND
        if(trim(adjustl(SFBOUN(mb)))==trim(adjustl(SFBOUN(mmb)))) then
          WRITE(ifle,'(3a)') ' ERR: Same BC name :[',
     &    trim(adjustl(SFBOUN(mb))),'] is used in grid file'
          ierrorFVG=1
          return
        endif
        ENDDO
      enddo
!
      do mb=NAME_I,NBOUND
      if(.not.(tempFlag(mb))) then
        WRITE(ifle,'(4a)') 
     & 'MSG: Boundary name: [',
     &   trim(SFBOUN(mb)),'] is NOT defined in ',cntlnam
        WRITE(ifle,'(2a)') 
     &     '    *** please check "boundary" in ',cntlnam
        WRITE(ifle,'(2a)') 
     &     '       1.miss type or missing in ',cntlnam
        WRITE(ifle,'(2a)') 
     &       '       2.wrong charactor codes(CR<->CR+LF etc.',
     &       ' wrong FTP mode(ASCII/Binary) causes this error)'
        WRITE(ifle,'(a)') '       3.wrong file, etc.'
	ierrorFVG=1
	return
      endif
      enddo
      deallocate(tempFlag)
!
!------------------------------------------------
! --- 
!------------------------------------------------
!
      allocate(tempFlag(0:nvrtx))
      allocate(givbcnd(0:NBOUND))
      givbcnd(:)=0
!
! --- 
!
      do 500 mb=1,NBOUND
      tempFlag(:)=.false.
      do 540 IS=1,NFBOUN(mb)
      do 510 Int=1,4
      IV=IFFACE(Int,IS,mb)
      tempFlag(IV)=.true.
 510  continue
 540  continue
      do iv=1,nvrtx
      if(tempFlag(iv)) then
        givbcnd(mb)=givbcnd(mb)+1
      endif
      enddo
 500  enddo
!
!
!
      do 550 mb=1,NBOUND
      write(ifll,'(a,I8,2x,a20,I8)')
     &   'MSG: BC NODE ',mb,trim(sfboun(mb)),givbcnd(mb)
      givbcnd(mb)=givbcnd(mb)+givbcnd(mb-1)
 550  continue
!
      allocate(glvbcnd(givbcnd(NBOUND)))
      glvbcnd(:)=0
!----------------------
! --- 
!----------------------
      tempInt=0
      do 600 mb=1,NBOUND
      tempFlag(:)=.false.
      do 640 IS=1,NFBOUN(mb)
      do 610 Int=1,4
      IV=IFFACE(Int,IS,mb)
      tempFlag(IV)=.true.
 610  continue
 640  continue
      do iv=1,nvrtx
      if(tempFlag(iv)) then
        tempInt=tempInt+1
        glvbcnd(tempInt)=iv
      endif
      enddo
      if(tempInt.ne.givbcnd(mb)) then
        WRITE(ifle,'(2a,3I8)')
     &  'MSG: ',
     &  'Number of boundary vertex is wrong.',mb,tempInt,givbcnd(mb)
	ierrorFVG=1
	stop
      endif
 600  continue
!
      allocate(icbinr(0:nbcnd))
      icbinr=0
      allocate(lcbinr(icbinr(nbcnd)))
!
!---------------------------------------------
! --- define periodic BC and touch-inlet BC --
!---------------------------------------------
!
      do 200 nb=1,nbcnd
      if(kdbcnd(0,nb)==kdprdc) then
!-------------------
! --- periodic BC: -
!-------------------
        if(idis(nb)==0) then
	  IBIDa=boundIDMap(nb,1)
 	  IBIDb=boundIDMap(nb,2)
	  if((givbcnd(IBIDa)-givbcnd(IBIDa-1))/=
     &       (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	     WRITE(ifle,'(1x,2a,2I4)')
     &      'MSG: Number of periodic boundary vertex is not matched.',
     &      (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &      (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	    ierrorFVG=1
            stop 'ERR: Stop at storeGridData subroutine'
	  endif
	  allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	  pMap(:)=-1
	  if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
	    do 210 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	    do 220 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
            if(((cord(2,glvbcnd(k))-cord(2,glvbcnd(j))-
     &           prdcOffset(nb,2))**2+
     &          (cord(3,glvbcnd(k))-cord(3,glvbcnd(j))-
     &           prdcOffset(nb,3))**2)
     &          .lt.smallNrm*scalfctr(nb)) then
	      pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &        =glvbcnd(k)
	      exit
	    endif
 220        continue
 210        continue
          elseif((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
	    do 230 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	    do 240 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
            if(((cord(1,glvbcnd(k))-cord(1,glvbcnd(j))-
     &          prdcOffset(nb,1))**2+
     &         (cord(3,glvbcnd(k))-cord(3,glvbcnd(j))-
     &         prdcOffset(nb,3))**2)<smallNrm*scalfctr(nb)) then
	      pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &        = glvbcnd(k)
  	      exit
            end if
 240        continue
 230        continue
          elseif((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
            do 250 j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
            do 260 k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
            if(((cord(2,glvbcnd(k))-cord(2,glvbcnd(j))-
     &            prdcOffset(nb,2))**2+
     &           (cord(1,glvbcnd(k))-cord(1,glvbcnd(j))-
     &           prdcOffset(nb,1))**2)<smallNrm*scalfctr(nb)) then
	      pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &        = glvbcnd(k)
	      exit
	    endif
 260        continue
 250        continue
          elseif((prdcAxis(nb)=='r').or.(prdcAxis(nb)=='R')) then
            rot(:,:)=0.d0
            kd=kdbcnd(0,nb)
            call set_rotmtrx1(nbcnd,kd,nb,rot)
            do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
!
            dum1=cord(1,glvbcnd(j))-prdcOffset(nb,1)
            dum2=cord(2,glvbcnd(j))-prdcOffset(nb,2)
            dum3=cord(3,glvbcnd(j))-prdcOffset(nb,3)
            dumy(1)=rot(1,1)*dum1
     &             +rot(2,1)*dum2
     &             +rot(3,1)*dum3
            dumy(2)=rot(1,2)*dum1
     &             +rot(2,2)*dum2
     &             +rot(3,2)*dum3
            dumy(3)=rot(1,3)*dum1
     &             +rot(2,3)*dum2
     &             +rot(3,3)*dum3
            dum4=1.D8
            do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
            dum1=(dumy(1)-(cord(1,glvbcnd(k))-prdcOffset(nb,1)))**2
     &          +(dumy(2)-(cord(2,glvbcnd(k))-prdcOffset(nb,2)))**2
     &          +(dumy(3)-(cord(3,glvbcnd(k))-prdcOffset(nb,3)))**2
            IF(dum1<smallNrm*scalfctr(nb)) then
              pMap(givbcnd(IBIDb-1)+j-givbcnd(IBIDa-1))
     &        = glvbcnd(k)
              exit
            endif
            enddo
            enddo
	  else
	    WRITE(ifle,'(1x,a)')  
     &      'MSG: Periodic axis setting is wrong.',
     &      prdcAxis(nb)
	    ierrorFVG=1
	    stop
	  endif
!
          do 270 j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
	  if((pMap(j)<1).or.(pMap(j)>nvrtx)) then
	    WRITE(ifle,'(1x,a,I8)')
     &      'ERR: Periodic boundary mapping error (no vertex).',pMap(j)
            WRITE(ifle,'(1x,a,I8)') 'MSG: BC no',nobcnd(nb)
	    ierrorFVG=1
	    stop
	  endif
 270      continue
!
          do 280 j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
          do 290 k=j+1,givbcnd(IBIDb)
	  if(pMap(j)==pMap(k)) then
	    WRITE(ifle,'(1x,2a)')  
     &      'ERR: ',
     &      'Periodic boundary mapping error (duplicated vertex).'
            WRITE(ifle,'(1x,2a)')  
     &      'MSG: decreasing value [&files/tolerance]'
	    ierrorFVG=1
	    stop
	  end if
 290      continue
 280      continue
!
          glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &    pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!
	  if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
	    do 300 j=1,givbcnd(IBIDa)-givbcnd(IBIDa-1)
	    if(((cord(2,glvbcnd(j+givbcnd(IBIDb-1)))-
     &         cord(2,glvbcnd(j+givbcnd(IBIDa-1)))-
     &         prdcOffset(nb,2))**2+
     &        (cord(3,glvbcnd(j+givbcnd(IBIDb-1)))-
     &         cord(3,glvbcnd(j+givbcnd(IBIDa-1)))-
     &         prdcOffset(nb,3))**2)>=smallNrm*scalfctr(nb)) then
               WRITE(ifle,'(1x,a)')
     &        'MSG: Periodic boundary mapping error (x-axis check).'
              ierrorFVG=1
              return
            endif
 300        continue
          elseif((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
	    do 310 j=1,givbcnd(IBIDa)-givbcnd(IBIDa-1)
	    if(((cord(1,glvbcnd(j+givbcnd(IBIDb-1)))-
     &         cord(1,glvbcnd(j+givbcnd(IBIDa-1)))-
     &         prdcOffset(nb,1))**2+
     &         (cord(3,glvbcnd(j+givbcnd(IBIDb-1)))-
     &         cord(3,glvbcnd(j+givbcnd(IBIDa-1)))- 
     &         prdcOffset(nb,3))**2)>=smallNrm*scalfctr(nb)) then
               WRITE(ifle,'(2a)') 'MSG: ',
     &        'Periodic boundary mapping error (y-axis check).'
              ierrorFVG=1
	      return
	    endif
 310        continue
          elseif((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
	    do 320 j=1,givbcnd(IBIDa)-givbcnd(IBIDa-1)
	    if(((cord(2,glvbcnd(j+givbcnd(IBIDb-1)))-
     &         cord(2,glvbcnd(j+givbcnd(IBIDa-1)))-
     &         prdcOffset(nb,2))**2+
     &        (cord(1,glvbcnd(j+givbcnd(IBIDb-1)))-
     &         cord(1,glvbcnd(j+givbcnd(IBIDa-1)))-
     &         prdcOffset(nb,1))**2)>=smallNrm*scalfctr(nb)) then
	      WRITE(ifle,'(2a)')  
     &        'MSG: ',
     &        'Periodic boundary mapping error (z-axis check).'
	      ierrorFVG=1
	      return
	    end if
 320        continue
	  end if
	  deallocate(pMap)
        elseif(idis(nb)==1) then
        endif
!
      elseif(kdbcnd(0,nb)==kdintr) then
!------------------------
! --- interface BC mesh:
!------------------------
        if(idis(nb)==0) then
          IBIDa=boundIDMap(nb,1)
 	  IBIDb=boundIDMap(nb,2)
          allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	  pMap(:)=-1
!
          if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &       (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	    WRITE(ifle,'(2a,2I10)') 'MSG: ',
     &      'Number of interface boundary vertex is not matched.',
     &      (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &      (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	    ierrorFVG=1
	    return
          else
            WRITE(ifll,'(a,I8)') 
     &      'MSG: Number of interface boundary NODE pair number =',
     &      givbcnd(IBIDa)-givbcnd(IBIDa-1)
	  endif
!
          do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	  do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
          dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &         (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &         (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
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
	    WRITE(ifle,'(2a)')
     &      'MSG: ',
     &      'Interface boundary mapping error (no vertex).'
	    ierrorFVG=1
	    stop
	  endif
          enddo
!
          do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
          do k=j+1,givbcnd(IBIDb)
	  if(pMap(j)==pMap(k)) then
	    WRITE(ifle,'(3a)')  
     &       'MSG: ',
     &       ' Interface boundary mapping error',
     &       ' (duplicated vertex on identical BC).'
	    ierrorFVG=1
	    stop
	  endif
          enddo
          enddo
!
          glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &    pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!
          deallocate(pMap)
        elseif(idis(nb)==1) then
          
        endif


      elseif(kdbcnd(0,nb)==kdbuff) then
!------------------------
! --- Buffle BC mesh:
!------------------------
        if(idis(nb)==0) then
          IBIDa=boundIDMap(nb,1)
 	  IBIDb=boundIDMap(nb,2)
          allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	  pMap(:)=-1
!
          if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &       (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	    WRITE(ifle,'(2a,2I10)') 'MSG: ',
     &      'Number of Buffle boundary vertex is not matched.',
     &      (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &      (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	    ierrorFVG=1
	    return
          else
            WRITE(ifll,'(a,2I8)') 
     &      'MSG: Number of Buffle boundary NODE pair number & nb=',
     &      givbcnd(IBIDa)-givbcnd(IBIDa-1),nb
	  endif
!
          do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	  do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
          dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &         (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &         (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
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
            
            WRITE(ifle,'(1X,a,I4)') 'MSG: nb=',nb
            WRITE(ifle,'(1X,a,I4)') 'pMap(j)= ',pMap(j)
            WRITE(ifle,'(1X,a,2G16.9)') 
     &      'MSG: scalfctr(nb)=  , smallNrm*scalfctr(nb)=',
     &      scalfctr(nb),smallNrm*scalfctr(nb)
	    WRITE(ifle,'(2a)')
     &      'MSG: ',
     &      'Buffle boundary mapping error (no vertex).'
            WRITE(ifle,'(2a)')
     &      'MSG: ',
     &      'Check &model/Buffle_shift'
	    ierrorFVG=1
	    stop
	  endif
          enddo
!
          do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
          do k=j+1,givbcnd(IBIDb)
	  if(pMap(j)==pMap(k)) then
            WRITE(ifle,'(1X,a,I4)') 'MSG: nb=' ,nb
            WRITE(ifle,'(1X,a,2G16.9)') 
     &   'MSG: scalfctr(nb)=  , smallNrm*scalfctr(nb)=',
     &   scalfctr(nb),smallNrm*scalfctr(nb)
	    WRITE(ifle,'(3a)')  
     &       'MSG: ',
     &       ' Buffle boundary mapping error',
     &       ' (duplicated vertex on identical BC).'
	    ierrorFVG=1
	    stop
	  endif
          enddo
          enddo
!
          glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &    pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!
          deallocate(pMap)
        elseif(idis(nb)==1) then
          
        endif

!--------------------------------------------------------------
      elseif(kdbcnd(0,nb)==kdshutr) then
!------------------------
! --- Shutter BC mesh:
!------------------------
        if(idis(nb)==0) then
          IBIDa=boundIDMap(nb,1)
 	  IBIDb=boundIDMap(nb,2)
          allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	  pMap(:)=-1
!
          if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &       (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	    WRITE(ifle,'(2a,2I10)') 'MSG: ',
     &      'Number of Shutter  boundary vertex is not matched.',
     &      (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &      (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	    ierrorFVG=1
	    return
          else
            WRITE(ifll,'(a,2I8)') 
     &      'MSG: Number of Shutter boundary NODE pair number & nb=',
     &      givbcnd(IBIDa)-givbcnd(IBIDa-1),nb
	  endif
!
          do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	  do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
          dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &         (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &         (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
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
            WRITE(ifle,'(1X,a,I4)') 'MSG: nb=',nb
            WRITE(ifle,'(1X,a,2G16.9)') 
     &   'MSG: scalfctr(nb)=  , smallNrm*scalfctr(nb)=',
     &   scalfctr(nb),smallNrm*scalfctr(nb)
	    WRITE(ifle,'(2a)')
     &      'MSG: ',
     &      'Shutter boundary mapping error (no vertex).'
	    ierrorFVG=1
	    stop
	  endif
          enddo
!
          do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
          do k=j+1,givbcnd(IBIDb)
	  if(pMap(j)==pMap(k)) then
            WRITE(ifle,'(1X,a,I4)') 'MSG: nb=' ,nb
            WRITE(ifle,'(1X,a,2G16.9)') 
     &   'MSG: scalfctr(nb)=  , smallNrm*scalfctr(nb)=',
     &   scalfctr(nb),smallNrm*scalfctr(nb)
	    WRITE(ifle,'(3a)')  
     &       'MSG: ',
     &       ' Shutter boundary mapping error',
     &       ' (duplicated vertex on identical BC).'
	    ierrorFVG=1
	    stop
	  endif
          enddo
          enddo
!
          glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &    pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!
          deallocate(pMap)
        elseif(idis(nb)==1) then
          
        endif


!--------------------------------------------------------------
      elseif(kdbcnd(0,nb)==kdpors) then
!------------------------
! --- Porous BC mesh:
!------------------------
        if(idis(nb)==0) then
          IBIDa=boundIDMap(nb,1)
 	  IBIDb=boundIDMap(nb,2)
          allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	  pMap(:)=-1
!
          if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &       (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	    WRITE(ifle,'(2a,2I10)') 'MSG: ',
     &      'Number of Porous boundary vertex is not matched.',
     &      (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &      (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	    ierrorFVG=1
	    return
          else
            WRITE(ifll,'(a,2I8)') 
     &      'MSG: Number of Porous boundary NODE pair number & nb=',
     &      givbcnd(IBIDa)-givbcnd(IBIDa-1),nb
	  endif
!
          do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	  do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
          dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &         (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &         (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
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
            WRITE(ifle,'(1X,a,I4)') 'MSG: nb=',nb
            WRITE(ifle,'(1X,a,2G16.9)') 
     &   'MSG: scalfctr(nb)=  , smallNrm*scalfctr(nb)=',
     &   scalfctr(nb),smallNrm*scalfctr(nb)
	    WRITE(ifle,'(2a)')
     &      'MSG: ',
     &      'Porous boundary mapping error (no vertex).'
	    ierrorFVG=1
	    stop
	  endif
          enddo
!
          do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
          do k=j+1,givbcnd(IBIDb)
	  if(pMap(j)==pMap(k)) then
            WRITE(ifle,'(1X,a,I4)') 'MSG: nb=' ,nb
            WRITE(ifle,'(1X,a,2G16.9)') 
     &   'MSG: scalfctr(nb)=  , smallNrm*scalfctr(nb)=',
     &   scalfctr(nb),smallNrm*scalfctr(nb)
	    WRITE(ifle,'(3a)')  
     &       'MSG: ',
     &       ' Porous boundary mapping error',
     &       ' (duplicated vertex on identical BC).'
	    ierrorFVG=1
	    stop
	  endif
          enddo
          enddo
!
          glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &    pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!
          deallocate(pMap)
        elseif(idis(nb)==1) then
          
        endif
!
      elseif(kdbcnd(0,nb)==kdsld) then
!----------------------
! --- sliding BC mesh:
!----------------------
        if(idis(nb)==0) then
          IBIDa=boundIDMap(nb,1)
 	  IBIDb=boundIDMap(nb,2)
          allocate(pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb)))
	  pMap(:)=-1
!
          if((givbcnd(IBIDa)-givbcnd(IBIDa-1)).ne.
     &       (givbcnd(IBIDb)-givbcnd(IBIDb-1))) then
	    WRITE(ifle,'(2a,2I8)') 'MSG: ',
     &      'Number of sliding boundary vertex is not matched.',
     &      (givbcnd(IBIDa)-givbcnd(IBIDa-1)),
     &      (givbcnd(IBIDb)-givbcnd(IBIDb-1))
	    ierrorFVG=1
	    stop
          else
            WRITE(ifll,'(a,I8)') 
     &      'MSG: Number of sliding boundary NODE pair number =',
     &      givbcnd(IBIDa)-givbcnd(IBIDa-1)
	  endif
!
          do j=givbcnd(IBIDa-1)+1,givbcnd(IBIDa)
	  do k=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)
!
          dum1=(cord(1,glvbcnd(k))-cord(1,glvbcnd(j)))**2+
     &         (cord(2,glvbcnd(k))-cord(2,glvbcnd(j)))**2+
     &         (cord(3,glvbcnd(k))-cord(3,glvbcnd(j)))**2
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
	    WRITE(ifle,'(2a,I8)')
     &      'MSG: ',
     &      'Sliding boundary mapping error (no vertex).',pMap(j)
	    ierrorFVG=1
	    stop
	  endif
          enddo
!
          do j=givbcnd(IBIDb-1)+1,givbcnd(IBIDb)-1
          do k=j+1,givbcnd(IBIDb)
	  if(pMap(j)==pMap(k)) then
	    WRITE(ifle,'(3a,3I8)')  
     &      'MSG: ',
     &      ' Sliding boundary mapping error',
     &      ' (duplicated vertex on identical BC).',j,k,pMap(j)
	   ierrorFVG=1
     	   stop
	  endif
          enddo
          enddo
!
          glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))=
     &    pMap(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
!
          deallocate(pMap)
        elseif(idis(nb)==1) then
        endif
      endif
 200  continue
!
!--------------------------------
! --- Creating ivbcnd; lvbcnd  --
!--------------------------------
!
      allocate(ivbcnd(0:nbcnd))
      allocate(ivpair(0:NBOUND))
      allocate(lvbcnd(givbcnd(NBOUND)))
      ivbcnd(:)=0
      lvbcnd(:)=0
      ivpair(:)=0
!
      do 400 nb=1,nbcndde
      IBIDa=boundIDMap(nb,1)
      IBIDb=boundIDMap(nb,2)
      ivbcnd(nb)=ivbcnd(nb-1)
      ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDa)-givbcnd(IBIDa-1)
      if((ivbcnd(nb)-ivbcnd(nb-1))/=(givbcnd(IBIDa)
     &  -givbcnd(IBIDa-1))) then
         WRITE(ifle,'(2a)') 'MSG: ',
     &  'Number of periodic boundary vertex is not matched.'
	ierrorFVG=1
	stop
      endif
!
      lvbcnd(ivbcnd(nb-1)+1:ivbcnd(nb))=
     &  glvbcnd(givbcnd(IBIDa-1)+1:givbcnd(IBIDa))
      
!
      if((kdbcnd(0,nb)==kdprdc.and.idis(nb)==0).or.
     &   (kdbcnd(0,nb)==kdsld .and.idis(nb)==0).or.
     &   (kdbcnd(0,nb)==kdbuff.and.idis(nb)==0).or.
     &   (kdbcnd(0,nb)==kdpors.and.idis(nb)==0).or.
     &   (kdbcnd(0,nb)==kdshutr.and.idis(nb)==0).or.
     &   (kdbcnd(0,nb)==kdintr.and.idis(nb)==0)) then
	 tempInt=ivbcnd(nb)
	 ivbcnd(nb)=ivbcnd(nb)+givbcnd(IBIDb)-givbcnd(IBIDb-1)
	 lvbcnd(tempInt+1:ivbcnd(nb))=
     &   glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
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
     &   glvbcnd(givbcnd(IBIDb-1)+1:givbcnd(IBIDb))
         ivpair(nb)=tempInt
      endif
 400  continue
!--------------------------
! --- end sub
!--------------------------
      deallocate(givbcnd)
      deallocate(glvbcnd)
      deallocate(tempFlag)
      ierrorFVG=0
!
      return
!
      end subroutine storeGridData
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine readFFRGrid_A
     & (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell,rtrnStat)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_FFRGFwords
      use module_io,only       : ifli,ifll,ifle
      use module_io,only       : getfil,ffgFormatKey
      use module_boundary,only : SFBOUN,NBOUND,boundIDMap,numbname,
     &                           IFFACE,NFBOUN
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
      integer,intent(in) :: mvrtx,mcell
      integer,intent(in) :: nvrtx,ncell
      real*8 ,intent(out) :: cord(3,mvrtx)
      integer,intent(out) :: lacell(  mcell)
      integer,intent(out) :: lvcell(8,mcell)
!
! --- [local entities]
!
      character(len=80),save :: fnam
      integer :: ifl,ios=0
      logical :: flag
      integer :: ic,jc,rtrnStat
      character*80 :: fileFormatID
      character*80 :: datasetID,commentStr,
     &                datatypeID,keywordStr,fieldID
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
!----------------------------------------------------------------------
      integer :: numofMaterial
      integer,allocatable :: material_ID(:)
      character*80,allocatable :: materialName(:)
!----------------------------------------------------------------------
!     number of boundary region. In this version, boundary definition
!     is not considered boundary condition so periodic or touchinlet 
!     boundary is defined with two boundary regions. In this case
!     numofBoundary is same with NBOUND
!----------------------------------------------------------------------
      integer :: numofBoundary
!----------------------------------------------------------------------
!     Maximum number of faces in a boundary region. This value is used
!     to allocate NFBOUN array. But this is not used when NFBOUN is refered.
!----------------------------------------------------------------------
      integer :: NBFS
!----------------------------------------------------------------------
!     These logical arrays are used to check for substituting values to 
!     all array elements.
!----------------------------------------------------------------------
      logical,allocatable :: chkflg_cord(:)
      logical,allocatable :: chkflg_IFFACE(:,:)
      logical,allocatable :: chkflg_lacell(:),chkflg_lvcell(:)
!
! --- 
!
      call getfil(ifl,fnam,'ffrgrid')
!
! --- confirm logical file id ifl is not used.
!
      inquire(unit=ifl,opened=flag)
      if(flag) then
        write(ifll,'(a)')'*** Error at (readFFRGrid)'
        write(ifll,'(a)')'*** Unit number ',ifl,'is alread used.'
        rtrnStat=1
        return
      end if
!
! --- confirm the ffr-grid file exists.
!
      inquire(file=fnam,exist=flag)
      if(.not.(flag)) then
        write(ifll,'(a)')'*** Error at (readFFRGrid)'
        write(ifll,'(3a)')'*** File ',trim(fnam),' does not exists.'
        rtrnStat=1
        return
      end if
!      
      open(ifl,file=fnam,form='formatted',
     &  status='unknown',action='read',iostat=ios)
      if(ios/=0) then
        write(ifll,'(2a)')'*** Cannot open FrontFlowRed Grid File:',
     &                                                trim(fnam)
        return
      end if
!
! --- read file format ID
!
      read(ifl,*) fileFormatID
      
      if(trim(adjustl(fileFormatID))==
     &            trim(adjustl(asciiv2_FFGF))) then
        call readFFRGridUGFv2(ios)
        if(ios/=0) then
          write(ifll,'(a)')'*** Error in reading FFR-GRID'
          rtrnStat=1;return
        end if
      end if
      
      
      close(ifl)
      rtrnStat=0
      
      return
      
      contains
      
!////////////////////////////////////////////////////////////////
!========================================================
      subroutine readFFRGridUGFv2(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      integer :: ios2
      

!     read header field.
      read(ifl,*)datasetID
      if(trim(adjustl(datasetID))/=trim(adjustl(newSet_FFGF))) then
        write(ifll,'(2a)')'No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      read(ifl,'(a)')commentStr
      read(ifl,*)datatypeID
      read(ifl,*)keywordStr
      read(ifl,'(a)')commentStr
      read(ifl,*,iostat=ios2)datasetSize
        if(ios2/=0) then
          write(ifll,'(a)')
     &    '*** Reading Error : FFR-GRID header datasize'
          rtrnStat2=1;return
        end if
      read(ifl,*)fieldID
      if(trim(adjustl(fieldID))/=trim(adjustl(gridHead_FFGF))) then
        write(ifll,*)'***No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      call readFFRGridHeaderField(ios2)
      if(ios2/=0) then
        write(ifll,*)'***Error in reading header field.'
        rtrnStat2=1;return
      end if

!     allocate memory region
      call allocFFRGridMemory(ios2)
      if(ios2/=0) then
        write(ifll,*)'***Error in allocating memory.'
        rtrnStat2=1;return
      end if

!     Initialize grid data
      call initFFRGridMemory(ios2)
      if(ios2/=0) then
        write(ifll,*)'***Error in initializing memory.'
        rtrnStat2=1;return
      end if

!     read data field.
      do while(.true.)
        read(ifl,*)datatypeID
        if(trim(adjustl(datatypeID))==trim(adjustl(fileend_FFGF))) then
          rtrnStat2=0;  exit
        else if(trim(adjustl(datatypeID))==
     &                        trim(adjustl(newSet_FFGF))) then
          write(ifll,*)'***Second dataset is detected.'
!         In this version, FFR-GRID is construced by only one dataset
          rtrnStat2=1;  exit
        else if(trim(adjustl(datatypeID))==
     &                      trim(adjustl(customData_FFGF))) then
          read(ifl,*)keywordStr
          if(trim(adjustl(keywordStr))/=
     &                        trim(adjustl(ffrGrid_FFGF))) then
            write(ifll,*)'***Non grid datafield is deteced.'
!          In this version, field skipping function is not implemented.
            rtrnStat2=1;  exit
          end if
          read(ifl,'(a)')commentStr
          write(ifll,*)'    ###    Reading ',
     &            trim(adjustl(commentStr)),'...'
          read(ifl,*,iostat=ios2)datasetSize
            if(ios2/=0) then
              write(ifll,*)'*** Reading Error : FFR-GRID datasize'
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
            write(ifll,*)'***Undefined datafield detected.',fieldID
            rtrnStat2=1;exit
          end if
          
        else
          write(ifll,*)'*** Undefined datatype is detected.',datatypeID
          rtrnStat2=1;  exit
        end if
      end do
      
!     check grid data
      if(rtrnStat2==0) then
        call checkFFRGridData(ios2)
        if(ios2/=0) then
          write(ifll,*)'*** Error in checking grid data.'
          rtrnStat2=1;return
        end if
      end if
      
      end subroutine readFFRGridUGFv2


!========================================================
      subroutine readFFRGridHeaderField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      character*80 :: fieldEndID
      
      read(ifl,*,iostat=ios2)fieldVersion
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID Header-1'
          rtrnStat2=1;return
        end if
      if(fieldVersion==1) then
!  ---  Read version 1 header field. FROM HERE -------------------------      
!       number of vertices
        read(ifl,*,iostat=ios2)intTmp1
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID header-2'
          rtrnStat2=1;return
        end if
        if(intTmp1/=nvrtx) then
          write(ifll,*)
     &      '*** Number of vertices are wrong. Check configuration.'
          write(ifll,*)'*** Control:',nvrtx,'  Grid:',intTmp1
          rtrnStat2=1;return
        end if
!       number of elements
        read(ifl,*,iostat=ios2)intTmp1
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID header-3'
          rtrnStat2=1;return
        end if
        if(intTmp1/=ncell) then
          write(ifll,*)
     &      '*** Number of elements are wrong. Check configuration.'
          write(ifll,*)'*** Control:',ncell,'  Grid:',intTmp1
          rtrnStat2=1;return
        end if
!       number of materials
        read(ifl,*,iostat=ios2)numofMaterial
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID header-4'
          rtrnStat2=1;return
        end if
!       number of boundary regions
        read(ifl,*,iostat=ios2)numofBoundary
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID header-5'
          rtrnStat2=1;return
        end if
        NBOUND=numofBoundary
!       number of boundary faces
        allocate(NFBOUN(1:NBOUND))
        NFBOUN=-1
        do ic2=1,NBOUND
          read(ifl,*,iostat=ios2)intTmp1,intTmp2
          if(ios2/=0) then
            write(ifll,*)'*** number of boundary region is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NBOUND)) then
            write(ifll,*)'*** Boundary region ID,',intTmp1,' is wrong.'
            rtrnStat2=1;return
          end if
          NFBOUN(intTmp1)=intTmp2
        end do
        if(minval(NFBOUN)<0) then
          write(ifll,*)'*** Face number of boundary',minloc(NFBOUN),
     &                              ' is wrong:',minval(NFBOUN)
          rtrnStat2=1;return
        end if
        NBFS=maxval(NFBOUN)
!  ---  Read version 1 header field. TO HERE -------------------------      
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(gridHeade_FFGF)))) then
        write(ifll,*)'*** Abnormal termination of header field.'
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
      allocate(IFFACE(1:4,1:NBFS,1:NBOUND))
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
      IFFACE(:,:,:)=0

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
        write(ifll,*)'*** Reading Error : FFR-GRID Elemtype-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 Element type field FROM HERE -------------------
        read(ifl,*,iostat=ios2)nmat_l
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID Elemtype-2'
          rtrnStat2=1;return
        end if
        if((nmat_l<0).or.(nmat_l>numofMaterial)) then
          write(ifll,*)
     &      '*** Material number in element type field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nmat_l
          read(ifl,*,iostat=ios2)intTmp1,tmpStr80_1
          if(ios2/=0) then
            write(ifll,*)'*** line number of element type is wrong.'
            rtrnStat2=1;return
          end if
          if(intTmp1==0) then
            write(ifll,*)'*** Material ID number is ZERO.'
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
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                  trim(adjustl(elemTypee_FFGF)))) then
        write(ifll,*)'*** Abnormal termination of element type field.'
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
        write(ifll,*)'*** Reading Error : FFR-GRID Boundarytype-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 boundary type field. FROM HERE ------------------
        read(ifl,*,iostat=ios2)nbound_l
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID Boundarytype-2'
          rtrnStat2=1;return
        end if
        if((nbound_l<0).or.(nbound_l>numofBoundary)) then
          write(ifll,*)
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
              write(ifll,*)'*** line number of boundary type is wrong.'
              rtrnStat2=1;return
            end if
            tmpStr80_2=' '
          end if
          if((intTmp1<1).or.(intTmp1>numofBoundary)) then
            write(ifll,*)'*** Boundary ID number is wrong.'
            rtrnStat2=1;return
          end if
!         In this version, boundary condition ID (intTmp2) and 
!         boundary name2 (tmpStr80_2) is not used.
          SFBOUN(intTmp1)=tmpStr80_1
        end do
!  ---  Read version 1 boundary type field. TO HERE ------------------
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(bndryTypee_FFGF)))) then
        write(ifll,*)'*** Abnormal termination of boundary type field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridBndryTypeField

!----------------------------------------------------------------------
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
        write(ifll,*)'*** Reading Error : FFR-GRID VertexCoordinate-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 vertices coordinate field. FROM HERE ---------
        read(ifl,*,iostat=ios2)nvertex_l
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID VertexCoordinate-2'
          rtrnStat2=1;return
        end if
        if((nvertex_l<0).or.(nvertex_l>nvrtx)) then
          write(ifll,*)'*** ',
     &        'Vertices number in vertex coordinate field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nvertex_l
          read(ifl,*,iostat=ios2)intTmp1,dblTmp1,dblTmp2,dblTmp3
          if(ios2/=0) then
            write(ifll,*)
     &        '*** line number of vertex coordinate is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>nvrtx)) then
            write(ifll,*)'*** Vertex ID number is wrong.'
            rtrnStat2=1;return
          end if
          cord(1,intTmp1)=dblTmp1
          cord(2,intTmp1)=dblTmp2
          cord(3,intTmp1)=dblTmp3
          chkflg_cord(intTmp1)=.true.
        end do
!  ---  Read version 1 vertices coordinate field. TO HERE -----------
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(vrtxCorde_FFGF)))) then
        write(ifll,*)
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
        write(ifll,*)'*** Reading Error : FFR-GRID BndryFace-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
! ---  Read version 1 boundary face field. FROM HERE -----------
        read(ifl,*,iostat=ios2)bndryTyp
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID BndryFace-2'
          rtrnStat2=1;return
        end if
        read(ifl,*,iostat=ios2)bndryID
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID BndryFace-3'
          rtrnStat2=1;return
        end if
        read(ifl,*,iostat=ios2)nface_l,nvrtxf_l
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID BndryFace-4'
          rtrnStat2=1;return
        end if
        if((bndryID<0).or.(bndryID>numofBoundary)) then
          write(ifll,*)
     &     '*** Boundary ID in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if((nface_l<0).or.(nface_l>NFBOUN(bndryID))) then
          write(ifll,*)
     &     '*** Face number in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if(nvrtxf_l/=4) then
          write(ifll,*)'*** Vertices number in a face is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nface_l
          read(ifl,*,iostat=ios2)
     &            intTmp1,intTmp2,intTmp3,intTmp4,intTmp5
          if(ios2/=0) then
            write(ifll,*)'*** line number of boundary face is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NFBOUN(bndryID))) then
            write(ifll,*)'*** Face ID number is wrong.'
            rtrnStat2=1;return
          end if
          IFFACE(1,intTmp1,bndryID)=intTmp2
          IFFACE(2,intTmp1,bndryID)=intTmp3
          IFFACE(3,intTmp1,bndryID)=intTmp4
          IFFACE(4,intTmp1,bndryID)=intTmp5
          chkflg_IFFACE(intTmp1,bndryID)=.true.
        end do
! ---  Read version 1 boundary face field. TO HERE -----------
      else
        write(ifll,*)'****Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(bndryFacee_FFGF)))) then
        write(ifll,*)'*** Abnormal termination of boundary face field.'
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

      read(ifl,*,iostat=ios2)fieldVersion
      if(ios2/=0) then
        write(ifll,*)'*** Reading Error : FFR-GRID Element-1'
        rtrnStat2=1;return
      end if
      if(fieldVersion==1) then
!  ---  Read version 1 element field. FROM HERE -----------------
        read(ifl,*,iostat=ios2)nelem_l
        if(ios2/=0) then
          write(ifll,*)'*** Reading Error : FFR-GRID Element-2'
          rtrnStat2=1;return
        end if
        if((nelem_l<0).or.(nelem_l>ncell)) then
          write(ifll,*)'*** Element ID in element field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nelem_l
          read(ifl,*,iostat=ios2)intTmp1,vertxTmp(1:8),intTmp2,intTmp3
          if(ios2/=0) then
            write(ifll,*)'*** line number of element is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>ncell)) then
            write(ifll,*)'*** Element ID number is wrong.'
            rtrnStat2=1;return
          end if
          lvcell(1:8,intTmp1)=vertxTmp(1:8)
          lacell(intTmp1)=intTmp2
! --- In this version cell type  data is not included from FFR-GRID.
          chkflg_lacell(intTmp1)=.true.
          chkflg_lvcell(intTmp1)=.true.
        end do
!  ---  Read version 1 element field. FROM TO -----------------
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,*,iostat=ios2)fieldEndID
      if((ios2/=0).or.(trim(adjustl(fieldEndID))/=
     &                        trim(adjustl(elemVrtxe_FFGF)))) then
        write(ifll,*)'*** Abnormal termination of element field.'
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
          write(ifll,*)'*** Vertex',ic2,'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      deallocate(chkflg_cord)
      
      do jc2=1,NBOUND
      do ic2=1,NFBOUN(jc2)
        if(.not.(chkflg_IFFACE(ic2,jc2))) then
          write(ifll,*)'*** Boundary face',ic2,'of boundary',jc2,
     &                                          'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      end do
      deallocate(chkflg_IFFACE)
      
      do ic2=1,ncell
        if(.not.(chkflg_lacell(ic2))) then
          write(ifll,*)'*** Material ID of cell',ic,'is not defined.'
          rtrnStat2=1;return
        end if
        if(.not.(chkflg_lvcell(ic2))) then
          write(ifll,*)'*** Vertices ID of cell',ic,'is not defined.'
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
      subroutine readFFRGrid
     & (mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell,rtrnStat)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_FFRGFwords
      use module_io,only       : ifli,ifll,ifle
      use module_io,only       : getfil,ffgFormatKey
      use module_boundary,only : SFBOUN,NBOUND,boundIDMap,numbname,
     &                           IFFACE,NFBOUN
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
      integer,intent(in) :: mvrtx,mcell
      integer,intent(in) :: nvrtx,ncell
      real*8 ,intent(out) :: cord(3,mvrtx)
      integer,intent(out) :: lacell(  mcell)
      integer,intent(out) :: lvcell(8,mcell)
!
! --- [local entities]
!
      character(len=80),save :: fnam
      integer :: ifl,ios=0
      logical :: flag
      integer :: ic,jc,rtrnStat
      character*80 :: fileFormatID
      character*80 :: datasetID,commentStr,
     &                datatypeID,keywordStr,fieldID
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
      integer :: NBFS
!---------------------------------------------------------------------
!     These logical arrays are used to check for substituting values to 
!     all array elements.
!---------------------------------------------------------------------
      logical,allocatable :: chkflg_cord(:)
      logical,allocatable :: chkflg_IFFACE(:,:)
      logical,allocatable :: chkflg_lacell(:),chkflg_lvcell(:)


      call getfil(ifl,fnam,'ffrgrid')
!     confirm logical file id ifl is not used.
      inquire(unit=ifl,opened=flag)
      if(flag) then
        write(ifll,*)'*** Error at (readFFRGrid)'
        write(ifll,*)'*** Unit number ',ifl,'is alread used.'
        rtrnStat=1
        return
      end if
!     confirm the ffr-grid file exists.
      inquire(file=fnam,exist=flag)
      if(.not.(flag)) then
        write(ifll,*)'*** Error at (readFFRGrid)'
        write(ifll,*)'*** File ',trim(fnam),' does not exists.'
        rtrnStat=1
        return
      end if
      
      open(ifl,file=fnam,form='unformatted',
     &  status='unknown',action='read',iostat=ios)
      if(ios/=0) then
        write(ifll,*)'*** Cannot open FrontFlowRed Grid File:',
     &                                                trim(fnam)
        return
      end if

!     read file format ID
      read(ifl)fileFormatID
      if(fileFormatID==unformv2_FFGF) then
        call readFFRGridUGFv2(ios)
        if(ios/=0) then
          write(ifll,*)'*** Error in reading FFR-GRID'
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
        write(ifll,*)'No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      read(ifl)commentStr
      read(ifl)datatypeID
      read(ifl)keywordStr
      read(ifl)commentStr
      read(ifl)datasetSize
      read(ifl)fieldID
      if(fieldID/=gridHead_FFGF) then
        write(ifll,*)'***No header field in file:',trim(fnam)
        rtrnStat2=1;return
      end if
      call readFFRGridHeaderField(ios2)
      if(ios2/=0) then
        write(ifll,*)'***Error in reading header field.'
        rtrnStat2=1;return
      end if

!     allocate memory region
      call allocFFRGridMemory(ios2)
      if(ios2/=0) then
        write(ifll,*)'***Error in allocating memory.'
        rtrnStat2=1;return
      end if

!     Initialize grid data
      call initFFRGridMemory(ios2)
      if(ios2/=0) then
        write(ifll,*)'***Error in initializing memory.'
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
          write(ifll,*)'***Second dataset is detected.'
!
! --- In this version, FFR-GRID is construced by only one dataset
!
          rtrnStat2=1;  exit
        else if(datatypeID==customData_FFGF) then
          read(ifl)keywordStr
          if(keywordStr/=ffrGrid_FFGF) then
            write(ifll,*)'***Non grid datafield is deteced.'
!
! --- In this version, field skipping function is not implemented.
!
            rtrnStat2=1;  exit
          end if
          read(ifl)commentStr
          write(ifll,*)'    ###    Reading ',
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
            write(ifll,*)'***Undefined datafield detected.',fieldID
            rtrnStat2=1;exit
          end if
          
        else
          write(ifll,*)'*** Undefined datatype is detected.',datatypeID
          rtrnStat2=1;  exit
        end if
      end do
      
!     check grid data
      if(rtrnStat2==0) then
        call checkFFRGridData(ios2)
        if(ios2/=0) then
          write(ifll,*)'*** Error in checking grid data.'
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
        if(intTmp1/=nvrtx) then
          write(ifll,*)
     &      '*** Number of vertices are wrong. Check configuration.'
          write(ifll,*)'*** Control:',nvrtx,'  Grid:',intTmp1
          rtrnStat2=1;return
        end if
! --- number of elements
        read(ifl)intTmp1
        if(intTmp1/=ncell) then
          write(ifll,*)
     &      '*** Number of elements are wrong. Check configuration.'
          write(ifll,*)'*** Control:',ncell,'  Grid:',intTmp1
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
        do ic2=1,NBOUND
          read(ifl,iostat=ios2)intTmp1,intTmp2
          if(ios2/=0) then
            write(ifll,*)'*** number of boundary region is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NBOUND)) then
            write(ifll,*)'*** Boundary region ID,',intTmp1,' is wrong.'
            rtrnStat2=1;return
          end if
          NFBOUN(intTmp1)=intTmp2
        end do
        if(minval(NFBOUN)<0) then
          write(ifll,*)'*** Face number of boundary',minloc(NFBOUN),
     &                              ' is wrong:',minval(NFBOUN)
          rtrnStat2=1;return
        end if
        NBFS=maxval(NFBOUN)
!  ---  Read version 1 header field. TO HERE -------------------------      
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=gridHeade_FFGF)) then
        write(ifll,*)'*** Abnormal termination of header field.'
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
      allocate(IFFACE(1:4,1:NBFS,1:NBOUND))

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

      rtrnStat2=0
      end subroutine allocFFRGridMemory

!----------------------------------------------------------------------
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
      IFFACE(:,:,:)=0

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

      rtrnStat2=0
      end subroutine initFFRGridMemory

!----------------------------------------------------------------------
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
          write(ifll,*)
     &      '*** Material number in element type field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nmat_l
          read(ifl,iostat=ios2)intTmp1,tmpStr80_1
          if(ios2/=0) then
            write(ifll,*)'*** line number of element type is wrong.'
            rtrnStat2=1;return
          end if
          if(intTmp1==0) then
            write(ifll,*)'*** Material ID number is ZERO.'
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
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=elemTypee_FFGF)) then
        write(ifll,*)'*** Abnormal termination of element type field.'
        rtrnStat2=1;return
      end if  
      
      rtrnStat2=0

      end subroutine readFFRGridElemTypeField

!----------------------------------------------------------------------
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
!  ---  Read version 1 boundary type field. FROM HERE -----------------
        read(ifl)nbound_l
        if((nbound_l<0).or.(nbound_l>numofBoundary)) then
          write(ifll,*)
     &      '*** Boundary number in boundary type field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nbound_l
          read(ifl,iostat=ios2)intTmp1,intTmp2,tmpStr80_1,tmpStr80_2
          if(ios2/=0) then
            write(ifll,*)'*** line number of boundary type is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>numofBoundary)) then
            write(ifll,*)'*** Boundary ID number is wrong.'
            rtrnStat2=1;return
          end if
!         In this version, boundary condition ID (intTmp2) and 
!         boundary name2 (tmpStr80_2) is not used.
          SFBOUN(intTmp1)=tmpStr80_1
        end do
!  ---  Read version 1 boundary type field. TO HERE ------------------
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=bndryTypee_FFGF)) then
        write(ifll,*)'*** Abnormal termination of boundary type field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridBndryTypeField

!----------------------------------------------------------------------
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
          write(ifll,*)'*** ',
     &        'Vertices number in vertex coordinate field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nvertex_l
          read(ifl,iostat=ios2)intTmp1,dblTmp1,dblTmp2,dblTmp3
          if(ios2/=0) then
            write(ifll,*)
     &     '*** line number of vertex coordinate is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>nvrtx)) then
            write(ifll,*)'*** Vertex ID number is wrong.'
            rtrnStat2=1;return
          end if
          cord(1,intTmp1)=dblTmp1
          cord(2,intTmp1)=dblTmp2
          cord(3,intTmp1)=dblTmp3
          chkflg_cord(intTmp1)=.true.
        end do
!  ---  Read version 1 vertices coordinate field. TO HERE -----------
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=vrtxCorde_FFGF)) then
        write(ifll,*)
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
          write(ifll,*)
     &      '*** Boundary ID in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if((nface_l<0).or.(nface_l>NFBOUN(bndryID))) then
          write(ifll,*)
     &   '*** Face number in boundary face field is wrong.'
          rtrnStat2=1;return
        end if
        if(nvrtxf_l/=4) then
          write(ifll,*)'*** Vertices number in a face is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nface_l
          read(ifl,iostat=ios2)intTmp1,intTmp2,intTmp3,intTmp4,intTmp5
          if(ios2/=0) then
            write(ifll,*)'*** line number of boundary face is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>NFBOUN(bndryID))) then
            write(ifll,*)'*** Face ID number is wrong.'
            rtrnStat2=1;return
          end if
          IFFACE(1,intTmp1,bndryID)=intTmp2
          IFFACE(2,intTmp1,bndryID)=intTmp3
          IFFACE(3,intTmp1,bndryID)=intTmp4
          IFFACE(4,intTmp1,bndryID)=intTmp5
          chkflg_IFFACE(intTmp1,bndryID)=.true.
        end do
! ---  Read version 1 boundary face field. TO HERE -----------
      else
        write(ifll,*)'****Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=bndryFacee_FFGF)) then
        write(ifll,*)'*** Abnormal termination of boundary face field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridBndryFaceField

!----------------------------------------------------------------------
!========================================================
      subroutine readFFRGridElemVrtxField(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: fieldVersion
      integer :: ios2,ic2
      integer :: nelem_l,vertxTmp(1:8)
      character*80 :: fieldEndID

      read(ifl)fieldVersion
      if(fieldVersion==1) then
!  ---  Read version 1 element field. FROM HERE -----------------
        read(ifl)nelem_l
        if((nelem_l<0).or.(nelem_l>ncell)) then
          write(ifll,*)'*** Element ID in element field is wrong.'
          rtrnStat2=1;return
        end if
        do ic2=1,nelem_l
          read(ifl,iostat=ios2)intTmp1,vertxTmp(1:8),intTmp2,intTmp3
          if(ios2/=0) then
            write(ifll,*)'*** line number of element is wrong.'
            rtrnStat2=1;return
          end if
          if((intTmp1<1).or.(intTmp1>ncell)) then
            write(ifll,*)'*** Element ID number is wrong.'
            rtrnStat2=1;return
          end if
          lvcell(1:8,intTmp1)=vertxTmp(1:8)
          lacell(intTmp1)=intTmp2
! --- In this version cell type  data is not included from FFR-GRID.
          chkflg_lacell(intTmp1)=.true.
          chkflg_lvcell(intTmp1)=.true.
        end do
!  ---  Read version 1 element field. FROM TO -----------------
      else
        write(ifll,*)'*** Ver.',fieldVersion,'data is not compatible.'
        rtrnStat2=1;return
      end if

!     check termination of data field.
      read(ifl,iostat=ios2)fieldEndID
      if((ios2/=0).or.(fieldEndID/=elemVrtxe_FFGF)) then
        write(ifll,*)'*** Abnormal termination of element field.'
        rtrnStat2=1;return
      end if
      
      rtrnStat2=0

      end subroutine readFFRGridElemVrtxField


!---------------------------------------------------------------------
!========================================================
      subroutine checkFFRGridData(rtrnStat2)
!========================================================
      integer,intent(out) :: rtrnStat2
      
      integer :: ios2,ic2,jc2
      
!     check all elements of arrays are defined.
      do ic2=1,nvrtx
        if(.not.(chkflg_cord(ic2))) then
          write(ifll,*)'*** Vertex',ic2,'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      deallocate(chkflg_cord)
      
      do jc2=1,NBOUND
      do ic2=1,NFBOUN(jc2)
        if(.not.(chkflg_IFFACE(ic2,jc2))) then
          write(ifll,*)'*** Boundary face',ic2,'of boundary',jc2,
     &                                          'is not defined.'
          rtrnStat2=1;return
        endif
      enddo
      enddo
      deallocate(chkflg_IFFACE)
!      
      do ic2=1,ncell
        if(.not.(chkflg_lacell(ic2))) then
          write(ifll,*)
     &       '*** Material ID of cell',ic,'is not defined.'
          rtrnStat2=1;return
        end if
        if(.not.(chkflg_lvcell(ic2))) then
          write(ifll,*)
     &       '*** Vertices ID of cell',ic,'is not defined.'
          rtrnStat2=1;return
        end if
      end do
      deallocate(chkflg_lacell)
      deallocate(chkflg_lvcell)
!
      rtrnStat2=0
      end subroutine checkFFRGridData
!
      end subroutine readFFRGrid
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$&$$$$$$$$$$$$$$$$$$$
