! **********************************************************************
      subroutine append_resultfv(fname,ifl,CPUnum,gridOnly,status)
! **********************************************************************
      use FVdata, only : hcell
      use FVheaderTag, only : FV_NODES,FV_FACES,FV_ELEMENTS,
     &                        FV_VARIABLES,FV_BNDRY_VARS
      use FFRdata, only : nvrtx,ncell,cord,NBFS,NBOUND,NFBOUN,
     &                    numVars,numBVars,lvcell,IBFACE,IFFACE,
     &                    kmesh,totdata
!
      implicit none
!
      character*80,intent(in) :: fname
      integer,intent(in) :: ifl,CPUnum
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
!
      integer :: grd,iMAT,iv
      integer :: i,j,k,l,ios,numSfBRG,nnum
      integer :: type,tempElems
      integer :: wallInfo4(1:4),wallInfo5(1:5),wallInfo6(1:6)
      integer :: headword
      integer :: ndCl(1:8)
      integer,allocatable :: ip(:)
      integer,allocatable :: ib(:)
      logical :: opn=.false.
!
! --- 
      allocate(ip(1:ncell))
      allocate(hcell(1:ncell),stat=ios)
      allocate(ib(1:NBFS),stat=ios)
      if(ios.ne.0) stop 'stop at allocating hcell in append_resultfv'
      ip=0
      hcell=0

!      inquire(file=fname,number=ios)
      inquire(file=fname,OPENED=opn)
!      if(ios>=0) then
      if(opn) then
        write(*,*)'(append_resultfv)File ',trim(fname),
     &            ' has been already opened as number ',ios
        stop 'ERR-1: append_resultfv'
        status=1; return
      end if
! --- 
      open(ifl,file=fname,form='unformatted',action='write',
     &                              position='append',iostat=ios)
      if(ios/=0) then
        write(*,*)'(append_resultfv) Cannot open ',trim(fname)
        status=1; return
      endif
!
!     --- OUTPUTING GRID COORDINATES
!
      write(ifl) FV_NODES, nvrtx
      write(ifl)(real(cord(1,i)),i=1,nvrtx),(real(cord(2,i)),i=1,nvrtx),
     &          (real(cord(3,i)),i=1,nvrtx)
!
!     --- OUTPUTING BOUNDARY-FACE DATA
!
      do i=1,NBOUND
         l=0
         ib=0
         do j=1,NBFS
         if(IBFACE(j)==i) then
            l=l+1
            ib(l)=j
         end if
         end do
         if(NFBOUN(i)/=l) then
            write(*,*)'(append_resultfv) Boundary faces is not match',
     &           l,NFBOUN(i)
            status=1; return
         end if
         write(ifl)FV_FACES,i,NFBOUN(i)
         write(ifl)((IFFACE(j,ib(k)),j=1,4),k=1,NFBOUN(i))
      end do
!
!     --- OUTPUTING ELEMENT DATA ???
!
      wallInfo4=0
      wallInfo5=0
      wallInfo6=0
!
      do type=1, 4
!        set up FVheader
!         1:tetra, 2:hex, 3:prsm, 4:prmd
!
!        *** note!
!        in real format, headword is identical for boundary
!         wall flags. int he following, this is not considered.
         ip=0
         tempElems=0
!        count cell num in types each other.
         do i = 1 , ncell
            if(kmesh(i)==type) then
               select case(type)
               case(1)
                  call ftn_encode_header(1,wallInfo4(1),headword)
               case(2)
                  call ftn_encode_header(2,wallInfo6(1),headword)
               case(3)
                  call ftn_encode_header(3,wallInfo5(1),headword)
               case(4)
                  call ftn_encode_header(4,wallInfo5(1),headword)
               case default
                  print *, '   ** unknown element type.'
               end select
               hcell(i)=headword
               tempElems=tempElems+1
!              pointer to cell num.
               ip(tempElems)=i
            end if
         end do

!        write out elemnts section
         if(tempElems>0) then
            select case(type)
            case(1)
!              make up lvcell order from FrontFlow to FIELDVIEW format.
               ndCl(1)=4
               ndCl(2)=1
               ndCl(3)=2
               ndCl(4)=3
               write(*,*) ' #### Tetra   elements=',tempElems
               write(ifl) FV_ELEMENTS,tempElems,0,0,0
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,4), 
     &              k=1,tempElems)
            case(2)
               ndCl(1)=1
               ndCl(2)=2
               ndCl(3)=4
               ndCl(4)=3
               ndCl(5)=5
               ndCl(6)=6
               ndCl(7)=8
               ndCl(8)=7
               write(*,*) ' #### Hexa    elements=',tempElems
               write(ifl) FV_ELEMENTS,0,tempElems,0,0
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,8), 
     &              k=1,tempElems)
            case(3)
               ndCl(1)=1
               ndCl(2)=4
               ndCl(3)=5
               ndCl(4)=2
               ndCl(5)=6
               ndCl(6)=3
               write(*,*) ' #### Prism   elements=',tempElems
               write(ifl) FV_ELEMENTS,0,0,tempElems,0
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,6), 
     &              k=1,tempElems)
            case(4)
               ndCl(1)=1
               ndCl(2)=2
               ndCl(3)=3
               ndCl(4)=4
               ndCl(5)=5
               write(*,*) ' #### Pyramid elements=',tempElems
               write(ifl) FV_ELEMENTS,0,0,0,tempElems
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,5), 
     &              k=1,tempElems)
            case default
               print *, '   ** unknown element type.'
            end select
         end if
!
      end do   !   do type=1, 4
!
!
!     --- OUTPUTING VARIABLES DATA
!
      write(ifl) FV_VARIABLES
      if(numVars>0) then
         write(ifl) ((real(totdata(i,j)),j=1,nvrtx),i=1,numVars)
      endif
!
!     --- OUTPUTING BC VARIABLES DATA
!
      write(ifl) FV_BNDRY_VARS
      if(numBVars>0) then
         write(*,*) 'ERR: FV Boundary Variables ',
     &              'does not support yet.'
!        write(ifl) bvar(:,:) ...
         status=1
         return
      elseif(numBVars==0) then
!        debug for FIELDVIEW 7.x
         numSfBRG=0
         do i=1,NBOUND
            numSfBRG=numSfBRG+NFBOUN(i)
         enddo
         do i=1,numSfBRG
            write(ifl) real(0.d0)
         end do
      endif
!      
      close(ifl)
      status=0

      deallocate(ip)
      deallocate(hcell)
      deallocate(ib)
!
      end subroutine append_resultfv

! **********************************************************************
      subroutine append_resultfv_ASCII(fname,ifl,CPUnum,gridOnly,status)
! **********************************************************************
      use FVdata, only : hcell
      use FVheaderTag, only : FV_NODES,FV_FACES,FV_ELEMENTS,
     &                        FV_VARIABLES,FV_BNDRY_VARS
      use FFRdata, only : nvrtx,ncell,cord,NBFS,NBOUND,NFBOUN,
     &                    numVars,numBVars,lvcell,IBFACE,IFFACE,
     &                    kmesh,totdata
!
      implicit none
!
      character*80,intent(in) :: fname
      integer,intent(in) :: ifl,CPUnum
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
!
      integer :: grd,iMAT,iv
      integer :: i,j,k,l,ios,numSfBRG,nnum
      integer :: type,tempElems
      integer :: wallInfo4(1:4),wallInfo5(1:5),wallInfo6(1:6)
      integer :: headword
      integer :: ndCl(1:8)
      integer,allocatable :: ip(:)
      integer,allocatable :: ib(:)
!
! --- 
      allocate(ip(1:ncell))
      allocate(hcell(1:ncell),stat=ios)
      allocate(ib(1:NBFS))
      if(ios.ne.0) stop 'stop at allocating hcell in append_resultfv'
      ip=0
      hcell=0

      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(append_resultfv_ASCII)File ',trim(fname),
     &            ' has been already opened as number ',ios
        stop 'ERR-2: append_resultfv_ASCII'
        status=1; return
      end if
! --- 
      open(ifl,file=fname,form='formatted',action='write',
     &                              position='append',iostat=ios)
!      open(ifl,file=fname,form='unformatted',action='write',
!     &                              position='append',iostat=ios)
      if(ios/=0) then
        write(*,*)'(append_resultfv_ASCII) Cannot open ',trim(fname)
        status=1; return
      endif
!
!     --- OUTPUTING GRID COORDINATES
!
      write(ifl,'(A)') 'Nodes'
      write(ifl,*) nvrtx
      do i=1,nvrtx
         write(ifl,*) real(cord(1,i)),real(cord(2,i)),real(cord(3,i))
      end do
!
!     --- OUTPUTING BOUNDARY-FACE DATA
!
      numSfBRG=0
      do i=1,NBOUND
         numSfBRG=numSfBRG+NFBOUN(i)
      end do
      write(ifl,'(A)') 'Boundary Faces'
      write(ifl,*) numSfBRG
      do i=1,NBOUND
         l=0
         ib=0
         do j=1,NBFS
         if(IBFACE(j)==i) then
           l=l+1
           ib(l)=j
         end if
         end do
         if(NFBOUN(i)/=l) then
            write(*,*)'(append_resultfv) Boundary faces is not match'
            status=1; return
         end if

         do k=1,NFBOUN(i)
            if((IFFACE(4,ib(k))==0) .or.
     &         (IFFACE(3,ib(k))==IFFACE(4,ib(k)))) then
               nnum=3
            else
               nnum=4
            endif
            write(ifl,*)i,nnum,(IFFACE(j,ib(k)),j=1,nnum)
         end do
      end do
!
!     --- OUTPUTING ELEMENT DATA ???
!
      wallInfo4=0
      wallInfo5=0
      wallInfo6=0
!
      write(ifl,'(A)') 'Elements'
      do i = 1 , ncell
        select case(kmesh(i))
!
        case(1) ! tet
!         make up lvcell order from FrontFlow to FIELDVIEW format.
          ndCl(1)=4
          ndCl(2)=1
          ndCl(3)=2
          ndCl(4)=3
!          write(ifl,'(6(I8,X))') 1,1,(lvcell(ndCl(j),i),j=1,4)
          write(ifl,*) 1,1,(lvcell(ndCl(j),i),j=1,4)
        case(2) ! hex
          ndCl(1)=1
          ndCl(2)=2
          ndCl(3)=4
          ndCl(4)=3
          ndCl(5)=5
          ndCl(6)=6
          ndCl(7)=8
          ndCl(8)=7
!          write(ifl,'(10(I8,X))')2,1,(lvcell(ndCl(j),i),j=1,8)
          write(ifl,*)2,1,(lvcell(ndCl(j),i),j=1,8)
        case(3) ! prsm
          ndCl(1)=1
          ndCl(2)=4
          ndCl(3)=5
          ndCl(4)=2
          ndCl(5)=6
          ndCl(6)=3
!          write(ifl,'(8(I8,X))') 3,1,(lvcell(ndCl(j),i),j=1,6)
          write(ifl,*) 3,1,(lvcell(ndCl(j),i),j=1,6)
        case(4) ! pyr
          ndCl(1)=1
          ndCl(2)=2
          ndCl(3)=3
          ndCl(4)=4
          ndCl(5)=5
!          write(ifl,'(7(I8,X))') 4,1,(lvcell(ndCl(j),i),j=1,5)
          write(ifl,*) 4,1,(lvcell(ndCl(j),i),j=1,5)
        end select
      end do

!
!
!     --- OUTPUTING VARIABLES DATA
!
      write(ifl,'(A)') 'Variables'
      if(numVars>0) then
         do i=1,numVars
            do j=1,nvrtx
               write(ifl,*)real(totdata(i,j))
            end do
         end do
      end if
!
!     --- OUTPUTING BC VARIABLES DATA
!
!      write(ifl,'(A)') 'Boundary Variables'
!      if(numBVars>0) then
!         write(*,*) 'ERR: FV Boundary Variables ',
!     &              'does not support yet.'
!!        write(ifl) bvar(:,:) ...
!         status=1
!         return
!      endif
!      
      close(ifl)
      status=0

      deallocate(ip)
      deallocate(hcell)
      deallocate(ib)
!
      end subroutine append_resultfv_ASCII

! **********************************************************************
      subroutine append_resultfv_SPLIT(fname_root,ifl,
     &                                 CPUnum,gridOnly,status)
! **********************************************************************
      use FVdata, only : hcell
      use FVheaderTag, only : FV_NODES,FV_FACES,FV_ELEMENTS,
     &                        FV_VARIABLES,FV_BNDRY_VARS
      use FFRdata, only : nvrtx,ncell,cord,NBFS,NBOUND,NFBOUN,
     &                    numVars,numBVars,lvcell,IBFACE,IFFACE,
     &                    kmesh,totdata
!
      implicit none
!
      character*80,intent(in) :: fname_root
      integer,intent(in) :: ifl,CPUnum
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
!
      integer :: grd,iMAT,iv
      integer :: i,j,k,l,ios,numSfBRG,nnum
      integer :: type,tempElems
      integer :: wallInfo4(1:4),wallInfo5(1:5),wallInfo6(1:6)
      integer :: headword
      integer :: ndCl(1:8)
      integer,allocatable :: ip(:)
      integer,allocatable :: ib(:)
      character*80 :: fname
!
! --- 
      allocate(ip(1:ncell))
      allocate(hcell(1:ncell),stat=ios)
      allocate(ib(1:NBFS))
      if(ios.ne.0) stop 'stop at allocating hcell in append_resultfv'
      ip=0
      hcell=0
!
!     ******** grid file *********
!
      fname=trim(adjustl(fname_root)) // '.gfv'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(append_resultfv)File ',trim(fname),
     &            ' has been already opened as number ',ios
        stop 'ERR-3: append_resultfv_SPLIT'
        status=1; return
      end if
! --- 
      open(ifl,file=fname,form='unformatted',action='write',
     &                              position='append',iostat=ios)
      if(ios/=0) then
        write(*,*)'(append_resultfv) Cannot open ',trim(fname)
        status=1; return
      endif
!
!     --- OUTPUTING GRID COORDINATES
!
      write(ifl) FV_NODES, nvrtx
      write(ifl)(real(cord(1,i)),i=1,nvrtx),(real(cord(2,i)),i=1,nvrtx),
     &          (real(cord(3,i)),i=1,nvrtx)
!
!     --- OUTPUTING BOUNDARY-FACE DATA
!
      do i=1,NBOUND
         l=0
         ib=0
         do j=1,NBFS
         if(IBFACE(j)==i) then
            l=l+1
            ib(l)=j
         end if
         end do
         if(NFBOUN(i)/=l) then
            write(*,*)'(append_resultfv) Boundary faces is not match'
            status=1; return
         end if
         write(ifl)FV_FACES,i,NFBOUN(i)
         write(ifl)((IFFACE(j,ib(k)),j=1,4),k=1,NFBOUN(i))
      end do
!
!     --- OUTPUTING ELEMENT DATA ???
!
      wallInfo4=0
      wallInfo5=0
      wallInfo6=0
!
      do type=1, 4
!        set up FVheader
!         1:tetra, 2:hex, 3:prsm, 4:prmd
!
!        *** note!
!        in real format, headword is identical for boundary
!         wall flags. int he following, this is not considered.
         ip=0
         tempElems=0
!        count cell num in types each other.
         do i = 1 , ncell
            if(kmesh(i)==type) then
               select case(type)
               case(1)
                  call ftn_encode_header(1,wallInfo4(1),headword)
               case(2)
                  call ftn_encode_header(2,wallInfo6(1),headword)
               case(3)
                  call ftn_encode_header(3,wallInfo5(1),headword)
               case(4)
                  call ftn_encode_header(4,wallInfo5(1),headword)
               case default
                  print *, '   ** unknown element type.'
               end select
               hcell(i)=headword
               tempElems=tempElems+1
!              pointer to cell num.
               ip(tempElems)=i
            end if
         end do

!        write out elemnts section
         if(tempElems>0) then
            select case(type)
            case(1)
!              make up lvcell order from FrontFlow to FIELDVIEW format.
               ndCl(1)=4
               ndCl(2)=1
               ndCl(3)=2
               ndCl(4)=3
               write(*,*) ' #### Tetra   elements=',tempElems
               write(ifl) FV_ELEMENTS,tempElems,0,0,0
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,4), 
     &              k=1,tempElems)
            case(2)
               ndCl(1)=1
               ndCl(2)=2
               ndCl(3)=4
               ndCl(4)=3
               ndCl(5)=5
               ndCl(6)=6
               ndCl(7)=8
               ndCl(8)=7
               write(*,*) ' #### Hexa    elements=',tempElems
               write(ifl) FV_ELEMENTS,0,tempElems,0,0
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,8), 
     &              k=1,tempElems)
            case(3)
               ndCl(1)=1
               ndCl(2)=4
               ndCl(3)=5
               ndCl(4)=2
               ndCl(5)=6
               ndCl(6)=3
               write(*,*) ' #### Prism   elements=',tempElems
               write(ifl) FV_ELEMENTS,0,0,tempElems,0
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,6), 
     &              k=1,tempElems)
            case(4)
               ndCl(1)=1
               ndCl(2)=2
               ndCl(3)=3
               ndCl(4)=4
               ndCl(5)=5
               write(*,*) ' #### Pyramid elements=',tempElems
               write(ifl) FV_ELEMENTS,0,0,0,tempElems
               write(ifl) (hcell(ip(k)),(lvcell(ndCl(j),ip(k)),j=1,5), 
     &              k=1,tempElems)
            case default
               print *, '   ** unknown element type.'
            end select
         end if
!
      end do   !   do type=1, 4
!
      close(ifl)
!
!
!     ******** result file *********
!
      fname=trim(adjustl(fname_root)) // '.rfv'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(append_resultfv)File ',trim(fname),
     &            ' has been already opened as number ',ios
        stop 'ERR-4: append_resultfv_SPLIT'
        status=1; return
      end if
! --- 
      open(ifl,file=fname,form='unformatted',action='write',
     &                              position='append',iostat=ios)
      if(ios/=0) then
        write(*,*)'(append_resultfv) Cannot open ',trim(fname)
        status=1; return
      endif
!
!     --- OUTPUTING GRID COORDINATES
!
      write(ifl) FV_NODES, nvrtx
!
!     --- OUTPUTING VARIABLES DATA
!
      write(ifl) FV_VARIABLES
      if(numVars>0) then
         write(ifl) ((real(totdata(i,j)),j=1,nvrtx),i=1,numVars)
      endif
!
!     --- OUTPUTING BC VARIABLES DATA
!
      write(ifl) FV_BNDRY_VARS
      if(numBVars>0) then
         write(*,*) 'ERR: FV Boundary Variables ',
     &              'does not support yet.'
!        write(ifl) bvar(:,:) ...
         status=1
         return
      elseif(numBVars==0) then
!        debug for FIELDVIEW 7.x
         numSfBRG=0
         do i=1,NBOUND
            numSfBRG=numSfBRG+NFBOUN(i)
         end do
         do i=1,numSfBRG
            write(ifl) real(0.d0)
         end do
      endif
!      
      close(ifl)

      status=0

      deallocate(ip)
      deallocate(hcell)
      deallocate(ib)
!
      end subroutine append_resultfv_SPLIT
