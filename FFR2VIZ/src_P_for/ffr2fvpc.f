!=============================================================================
      module ffr2fvpc_header
!=============================================================================
! ver.4.1.0
!=============================================================================
!--- Definitions
      integer,parameter :: mcomp               = 16
      integer,parameter :: MAX_FILE_NUMBER     = 1000
      integer,parameter :: MAX_INJECTOR        = 50
      integer,parameter :: FILE_NAME_LENGTH    = 15
      real*8,parameter  :: eps                 = 1.0e-30
      character(len=17) :: INFO_FILE_NAME      = "ffr2fvpc_info.txt"
      character(len=20) :: RESULT_FILE_NAME    = 'particle.fvp'
!
!--- Structures 
!
! - hold filename and related information list
      type t_file_info_list
        character(len=FILE_NAME_LENGTH) :: cpFileName
        integer :: iParcels(MAX_INJECTOR)
        integer :: iFIndex
        type(t_file_info_list),pointer :: next
      end type t_file_info_list

! - hold information writtern in "info_file_name" 
! -      and pointer to first element of file_info_list
      type t_ffrpt_info
        integer :: iFileNum
        integer :: iMaxParcels(MAX_INJECTOR)
        type(t_file_info_list),pointer :: first
      end type t_ffrpt_info

! - hold particle information list 
      type t_particle_list 
        integer :: iIndex
        real*8  :: dXYZ(3)
        real*8  :: dUVW(3)

        real*8  :: dTime
        real*8  :: dDura
        real*8  :: dDiam

        real*8  :: dInXYZ(3)
        real*8  :: dPcls
        real*8  :: dPMASS
        real*8  :: dTMP
        real*8  :: dPYS(mcomp)
        type(t_particle_list),pointer :: next
      end type t_particle_list

! - hold particle path begin/end position related to particle_list
      type t_path_list
        integer :: iNum
        type(t_particle_list),pointer :: l_begin
        type(t_particle_list),pointer :: l_end
      end type t_path_list

      end module ffr2fvpc_header
!=============================================================================

!=============================================================================
      program ffr2fvpc
!=============================================================================
      use ffr2fvpc_header
      implicit none

!--- Local variables
      integer,parameter :: iflr=10
      integer,parameter :: ifli=11
      integer,parameter :: iflp=12
      integer :: i, j, k, iFiles,kk,ii
      integer :: ierror
      integer :: idum1,idum2
      integer :: num_injection,num_path
      integer :: icom,ncomp
      character(len=256) :: cdum1
      character(len=20)  :: spenam(mcomp)
      character(len=5)   :: cdum3
      character(len=40)  :: cfmt
      real*8 :: dum_ms,dum_tmp,dum_pal
      real*8 :: XYZinj(3)
      real*8 :: ddum4(3),ddum5(3),ddum6,dum_D,dum_H
      real*8 :: dumm(mcomp),starttime
      type(t_ffrpt_info) :: theInfo
      type(t_file_info_list),pointer :: infoList,infoList_n,infoList_e
      type(t_particle_list),pointer  :: pL_n,currentPlistPtr
      type(t_path_list),allocatable  :: path(:)
!      integer, save,allocatable :: starttime
!
!--- Start
!

      write(*,'(a,/,a,/,a)')
     &      "+","++  particle path converter  v.4.1.0","+"
!
!--- Read & Store File Info List
!
      num_injection = 0
      num_path = 0
      do i = 1, MAX_INJECTOR
        theInfo%iMaxParcels(i) = -1
      enddo

      open(unit=ifli,file=INFO_FILE_NAME,iostat=ierror)
      if( ierror > 0 ) then
        write(*,*) "ERR: File Open Error at ",INFO_FILE_NAME
        stop
      endif
      read(ifli,*) cdum1, theInfo%iFileNum,RESULT_FILE_NAME
      write(*,'(a,a3,i4)') trim(cdum1)," = ",theInfo%iFileNum
      iFiles = theInfo%iFileNum

      i = 1
      allocate( infoList_e )
      nullify( infoList_e%next )
      read(ifli,*) cdum1
      infoList_e = T_FILE_INFO_LIST( cdum1,0,i,infoList_e%next )
      theInfo%first => infoList_e
      do i=2,iFiles
        read(ifli,*) cdum1
        allocate( infoList_n )
        nullify( infoList_n%next )
        infoList_n = T_FILE_INFO_LIST( cdum1,0,i,infoList_n%next )
        infoList_e%next => infoList_n
        infoList_e => infoList_n
      enddo
      close(ifli)
!
!--- Check the number of files
!
      if( MAX_FILE_NUMBER .lt. theInfo%iFileNum ) then
        write(*,*) "ERR: Too Many Files. ",
     &    MAX_FILE_NUMBER," < ",theInfo%iFileNum
        stop
      endif
!
!--- allocate path list
!
      write(*,'(a)') "search for max. particle index ..."
      !!!---- this loop can remove if max. particle ID get in some way
      infoList => theInfo%first !!!------- num_path ------------------
      do i = 1, theInfo%iFileNum
        open(unit=iflp,file=infoList%cpFileName,iostat=ierror)
        if( ierror > 0 ) then
          write(*,*) "ERR: File Open Error at ",infoList%cpFileName
          stop
        endif

        read(iflp,'(2i6,a6,g12.7,a)') idum1,idum2,cdum3,ddum6,cdum1
        num_injection = idum1
!        allocate(starttime(idum1))
        if( num_injection .gt. MAX_INJECTOR ) then
          write(*,*) "ERR: Too Many Injectors. MAX is ",MAX_INJECTOR
          stop
        endif

        read(iflp,*) ncomp
        if( ncomp > mcomp ) then
          print *,'This code is ONLY for ncomp=',mcomp
          stop 9999
        endif

        do icom = 1, ncomp
          read(iflp,*) spenam(icom)
        enddo

        do k = 1, num_injection
          read(iflp,'(i10,4g15.8)') idum1,XYZinj,starttime
!          starttime=ddum6
          if( idum1 <= 0 ) cycle
          
!
          do j = 1, idum1
            read(iflp,*) idum2,ddum4(1),ddum4(2),ddum4(3),
     &                   ddum5(1),ddum5(2),ddum5(3),
     &                   dum_D,dum_H,
     &                   dum_pal,dum_ms,dum_tmp,(dumm(kk),kk=1,ncomp)
          if( num_path < idum2 ) num_path=idum2
          enddo
        enddo
        close(iflp)
        infoList => infoList%next
      enddo !!!--------------------------- num_path ------------------

      allocate( path(num_path) )
      do i=1,num_path
        path(i)%iNum = 0
        nullify( path(i)%l_begin )
        nullify( path(i)%l_end   )
      enddo
!
!--- Read & Store Particle List & pointer list of Particle list
!
      infoList => theInfo%first
      do i = 1, theInfo%iFileNum
        write(*,'(a,a)') "Reading File:", infoList%cpFileName
        open(unit=iflp,file=infoList%cpFileName,iostat=ierror)
        if( ierror > 0 ) then
          write(*,*) "ERR: File Open Error at ",infoList%cpFileName
          stop
        endif

        read(iflp,'(2i6,a6,g12.7,a)') idum1,idum2,cdum3,ddum6,cdum1
        num_injection = idum1
        if( num_injection .gt. MAX_INJECTOR ) then
          write(*,*) "ERR: Too Many Injectors. MAX is ",MAX_INJECTOR
          stop
        endif

        read(iflp,*) ncomp
        if( ncomp > mcomp ) then
          print*,'This code is ONLY for ncomp=',mcomp
          stop 9999
        endif

        do icom = 1, ncomp
          read(iflp,*) spenam(icom)
        enddo

        do k = 1, num_injection
          read(iflp,'(i10,4g15.8)') idum1,XYZinj,starttime  !1111
!          starttime=ddum6
          if( idum1 <= 0 ) cycle
          infoList%iParcels(k) = idum1
          write(*,*) k," Parcels : ", infoList%iParcels(k)
          if( theInfo%iMaxParcels(k) < infoList%iParcels(k))
     &        theInfo%iMaxParcels(k) = infoList%iParcels(k)

          do j = 1, infoList%iParcels(k)
            read(iflp,*) idum1,ddum4(1),ddum4(2),ddum4(3),
     &                   ddum5(1),ddum5(2),ddum5(3),
     &                   dum_D,dum_H,
     &                  dum_pal,dum_ms,dum_tmp,(dumm(kk),kk=1,ncomp)

            do ii=1, 3
              if( abs(ddum4(ii)) < eps ) ddum4(ii) = 0.0
              if( abs(ddum5(ii)) < eps ) ddum5(ii) = 0.0
            enddo
            if( abs(dum_D) < eps ) dum_D = 0.0
            if( abs(dum_H) < eps ) dum_H = 0.0
            if( abs(dum_ms) < eps ) dum_ms = 0.0
            if( abs(dum_pal) < eps ) dum_pal = 0.0
            if( abs(dum_tmp) < eps ) dum_tmp = 0.0
            do ii=1, ncomp
              if( abs(dumm(ii)) < eps ) dumm(ii) = 0.0
            enddo
            allocate( pL_n )
            nullify( pL_n%next )
            pL_n = T_PARTICLE_LIST(idum1,ddum4,ddum5,ddum6,
     &                             dum_H,dum_D,XYZinj,
     &                             dum_pal,dum_ms,dum_tmp,dumm,
     &                             pL_n%next)

            if( path(idum1)%iNum .eq. 0 ) then 
              path(idum1)%l_begin => pL_n
              path(idum1)%l_end   => pL_n
            else
              path(idum1)%l_end%next => pL_n
              path(idum1)%l_end      => pL_n
            endif
            path(idum1)%iNum = path(idum1)%iNum+1 
          enddo
        enddo
        close(iflp)
        infoList => infoList%next
      enddo
!
! --- Write Output File : FieldView particle path format : header
!
      write(*,'(a)') "output file open ..."
      open(unit=iflr,file=trim(RESULT_FILE_NAME),iostat=ierror)
      if( ierror > 0 ) then
        write(*,*) "ERR: File Open Error at ",RESULT_FILE_NAME
        stop
      endif
      write(iflr,'(a,/,a)') 'FVPARTICLES 1 1', 'Variable Names'
      write(cdum3,'(i5)') ncomp+6
      write(iflr,'(a)') trim(adjustl(cdum3))
      write(iflr,'(a,/,a,/,a,/,a)') 'Duration', 
     &    'TimeHis','Diameter', 'Parcel'
      write(iflr,'(a,/,a)') 'PMASS', 'PTEM'
      do icom = 1, ncomp
        write(iflr,'(a)') trim(spenam(icom))
      enddo

      write(*,'(a,i8)') " particle path = ",num_path

      write(cdum3,'(i5)') ncomp+9

      write(cfmt,'(3a)') '(',trim(adjustl(cdum3)),'(1pe14.6))'


      
      do i = 1, num_path
        if( mod(i,500) .eq. 0 ) write(*,'(a,i8)') "Progress:",i
        if( path(i)%iNum <= 0 ) cycle



        currentPlistPtr => path(i)%l_begin
!        write(iflr,'(i6)') 1 + path(i)%iNum
        write(iflr,'(i6)') path(i)%iNum
        write(iflr,cfmt)
!     &    currentPlistPtr%dInXYZ(1),
!     &    currentPlistPtr%dInXYZ(2),
!     &    currentPlistPtr%dInXYZ(3),
!     &    0.0,
!     &    currentPlistPtr%dDura,
!     &    currentPlistPtr%dDiam,
!     &    currentPlistPtr%dPcls,
!     &    currentPlistPtr%dPMASS,
!     &    currentPlistPtr%dTMP,
!     &    (currentPlistPtr%dPYS(k),k=1,ncomp)

!        write(iflr,cfmt)
!     &    currentPlistPtr%dInXYZ(1),
!     &    currentPlistPtr%dInXYZ(2),
!     &    currentPlistPtr%dInXYZ(3),
!     &    currentPlistPtr%dTime,
!     &    currentPlistPtr%dDura,
!     &    currentPlistPtr%dDiam,
!     &    currentPlistPtr%dPcls,
!     &    currentPlistPtr%dPMASS,
!     &    currentPlistPtr%dTMP,
!     &    (currentPlistPtr%dPYS(k),k=1,ncomp)

        do while( associated(currentPlistPtr) )
          write(iflr,cfmt)
     &      currentPlistPtr%dXYZ(1),
     &      currentPlistPtr%dXYZ(2),
     &      currentPlistPtr%dXYZ(3),

     &      currentPlistPtr%dTime,
     &      currentPlistPtr%dDura,
     &      currentPlistPtr%dDiam,
     &      currentPlistPtr%dPcls,
     &      currentPlistPtr%dPMASS,
     &      currentPlistPtr%dTMP,
     &      (currentPlistPtr%dPYS(k),k=1,ncomp)
          currentPlistPtr => currentPlistPtr%next
        enddo
      enddo

      write(*,'(a)') "output file close ..."
      close( iflr )
      deallocate( path )
      write(*,'(a,/,a,/,a)') "+","++  terminated normally.","+"

      end program
!=============================================================================
