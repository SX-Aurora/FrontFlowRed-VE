c======================================================================c
      module tovtk
c
c @author     Takashi Matsushita
c @date       $Date: 2005/10/13 09:15:30 $
c @version    $Revision: 1.5 $
c
c Copyright:  (c) 2005 Takashi Matsushita, All rights reserved.
c License:    none
c Created:    27 September 2005
c Credits:    author(s) of append_resultfv.f
c
c @todo       nope
c @warnings   nope
c
c======================================================================c

      implicit none

c everything is private otherwise explicitly declared as public
      private

c----------------------------------------------------------------------c
c constants
c----------------------------------------------------------------------c
      character(len=*),
     &  parameter :: VTK_V2 = "# vtk DataFile Version 2.0"
      character(len=*),
     &  parameter :: VTK_ASCII = "ASCII"
      character(len=*),
     &  parameter :: VTK_UNSTRUCTURED_GRID = "DATASET UNSTRUCTURED_GRID"
      character(len=*),
     &  parameter :: VTK_POINT = "POINTS"
      character(len=*),
     &  parameter :: VTK_CELL = "CELLS"
      character(len=*),
     &  parameter :: VTK_CELL_TYPE = "CELL_TYPES"
      character(len=*),
     &  parameter :: VTK_CELL_DATA = "CELL_DATA"
      character(len=*),
     &  parameter :: VTK_VECTOR = "VECTORS"
      character(len=*),
     &  parameter :: VTK_SCALAR = "SCALARS"
      character(len=*),
     &  parameter :: VTK_LOOKUP_TABLE = "LOOKUP_TABLE"
      character(len=*),
     &  parameter :: VTK_POINT_DATA = "POINT_DATA"

      integer, parameter :: VTK_TETRA = 10
      integer, parameter :: VTK_HEXAHEDRON = 12
      integer, parameter :: VTK_WEDGE = 13
      integer, parameter :: VTK_PYRAMID = 14

c----------------------------------------------------------------------c
c public interfaces
c----------------------------------------------------------------------c
      public :: ffr_to_vtk

c----------------------------------------------------------------------c
c private interfaces
c----------------------------------------------------------------------c
      private :: open_vtk_file
      private :: close_vtk_file
      private :: write_vtk_header
      private :: write_vtk_points
      private :: write_vtk_cells
      private :: write_vtk_boundary
      private :: write_vtk_scalar
      private :: write_vtk_vector
      private :: ffr_to_vtk_ascii
      private :: ffr_to_vtk_binary

      contains
c======================================================================c
c implementations
c======================================================================c
c----------------------------------------------------------------------c
c private functions/subroutines
c----------------------------------------------------------------------c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function open_vtk_file(fileName) result (lun)
c----------------------------------------------------------------------
c @brief  opens a file for write
c @note   none
c @param  fileName [in] name of file to open
c @return logical unit number associated with opend file,
c         -1 if error occured
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c

        implicit none

        character(len=*), intent(in) :: fileName
        integer :: lun

        character(len=*), parameter :: FNAME = "open_vtk_file: "
        integer :: rc = 0

        write(*, *)'dbg> ', FNAME, trim(fileName), ' opening...'

        open(unit=lun,
     &       file=trim(fileName),
     &       status='unknown',
     &       form='formatted',
     &       action='write',
     &       iostat=rc)

        if (rc .ne. 0) then
          write(*, *)'err> ', FNAME, trim(fileName), ' open ', lun, rc
          lun = -1
        end if

        return
      end function open_vtk_file


c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function close_vtk_file(lun) result (rc)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c @brief  closes a file
c @note   none
c @param  lun [in] logical unit number of opened file
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c

        implicit none

        integer, intent(in) :: lun
        integer :: rc

        character(len=*), parameter :: FNAME = "close_vtk_file: "

        close(unit=lun, iostat=rc)

        if (rc .ne. 0) then
          write(*, *)'err> ', FNAME, ' close'
        end if

        return
      end function close_vtk_file

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function write_vtk_header(lun) result (rc)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c @brief  writes header part of vtk file
c @note   none
c @param  lun [in] logical unit number of opened file
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
!#undef LABEL
!#define LABEL 9010

        implicit none

        integer, intent(in) :: lun
        integer :: rc

        character(len=*), parameter :: FNAME = "write_vtk_header: "

        write(*, *)'dbg> ', FNAME, ' writing...'

        write(lun, fmt="(A)", iostat=rc, err=9010) trim(VTK_V2)
        write(lun, fmt="(A)", iostat=rc, err=9010) 'comment if any'
        write(lun, fmt="(A,/)", iostat=rc, err=9010) VTK_ASCII
        write(lun,
     &        fmt="(A)",
     &        iostat=rc,
     &        err=9010)
     &          VTK_UNSTRUCTURED_GRID

 9010   return
      end function write_vtk_header

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function write_vtk_points(lun) result (rc)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c @brief  writes coordinate points of nodes
c @note   none
c @param  lun [in] logical unit number of opened file
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
!#undef LABEL
!#define LABEL 9020

        use FFRdata, only : nvrtx, cord

        implicit none

        integer, intent(in) :: lun
        integer :: rc

        character(len=*), parameter :: FNAME = "write_vtk_points: "
c       character(len=*), parameter ::fmt_data = "(3E)"
! Modified by Y.Takahashi, 2021.11.16 ---------------------
!        character(len=*), parameter ::fmt_data = "(9E)"
        character(len=*), parameter ::fmt_data = "(9E20.7)"
! ---------------------------------------------------------
        integer :: ii

        write(*, *)'dbg> ', FNAME, ' writing...'

        write(lun,
     &        fmt="(A,X,I0,X,A)",
     &        iostat=rc,
     &        err=9020)
     &        VTK_POINT,
     &        nvrtx,
     &        "float"
        write(lun,
     &        fmt=fmt_data,
     &        iostat=rc,
     &        err=9020)
     &          (real(cord(1, ii)),
     &           real(cord(2, ii)),
     &           real(cord(3, ii)),
     &          ii=1, nvrtx)

 9020   return
      end function write_vtk_points

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function write_vtk_cells(lun) result (rc)
c----------------------------------------------------------------------
c @brief  writes cell data
c @note   none
c @param  lun [in] logical unit number of opened file
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
!#undef LABEL
!#define LABEL 9030

        use FFRdata, only : ncell, kmesh, lvcell

        implicit none

        integer, intent(in) :: lun
        integer :: rc

        character(len=*), parameter :: FNAME = "write_vtk_cells: "
c       character(len=*), parameter ::fmt_data = "(I0)"
c       integer, parameter :: newline = 1
        character(len=*), parameter ::fmt_data = "(1X,I0,$)"
        integer, parameter :: newline = 25

        integer, parameter :: nodes_per_tetra = 4
        integer, parameter :: nodes_per_hexa = 8
        integer, parameter :: nodes_per_wedge = 6
        integer, parameter :: nodes_per_pyramid = 5
        integer, parameter :: max_nodes_per_cell = 8

        integer :: ii, jj
        integer :: num_tetra=0, num_hexa=0, num_wedge=0, num_pyramid=0
        integer :: num_tokens, nodes
        integer :: nodes_map(1:max_nodes_per_cell)
        character(len=32) :: fmt_string

c count number of cells for each cell type
        do ii = 1, ncell

          select case(kmesh(ii))
            case(1)
              num_tetra = num_tetra + 1
            case(2)
              num_hexa = num_hexa + 1
            case(3)
              num_wedge = num_wedge + 1
            case(4)
              num_pyramid = num_pyramid + 1
            case default
              write(*, *)'err> ', FNAME, ' unknown type: ', kmesh(ii)
              stop
          end select

        end do

        write(*, *)'dbg> ', FNAME, '# of tetra cells = ', num_tetra
        write(*, *)'dbg> ', FNAME, '# of hexa cells = ', num_hexa
        write(*, *)'dbg> ', FNAME, '# of wedge cells = ', num_wedge
        write(*, *)'dbg> ', FNAME, '# of pyramid cells = ', num_pyramid

        num_tokens = num_tetra*(nodes_per_tetra + 1) +
     &               num_hexa*(nodes_per_hexa + 1) +
     &               num_wedge*(nodes_per_wedge + 1) +
     &               num_pyramid*(nodes_per_pyramid + 1)

        write(lun,
     &        fmt="(/,A,X,I0,X,I0)",
     &        iostat=rc,
     &        err=9030)
     &          VTK_CELL,
     &          ncell,
     &          num_tokens


c vtk has zero offset for node id
        lvcell = lvcell - 1

        write(*, *)'dbg> ', FNAME, ' writing...'

c dump node connectivity for each cell
        do ii = 1, ncell

          do jj = 1, max_nodes_per_cell
            nodes_map(jj) = jj
          end do

          select case(kmesh(ii))
            case(1)
              nodes = nodes_per_tetra
              fmt_string = "(5(I0,X))"
            case(2)
              nodes = nodes_per_hexa
              fmt_string = "(9(I0,X))"
            case(3)
              nodes_map(1) = 1
              nodes_map(2) = 3
              nodes_map(3) = 2
              nodes_map(4) = 4
              nodes_map(5) = 6
              nodes_map(6) = 5

              nodes = nodes_per_wedge
              fmt_string = "(7(I0,X))"
            case(4)
              nodes = nodes_per_pyramid
              fmt_string = "(6(I0,X))"
            case default
              write(*, *)'err> ', FNAME, ' unknown type: ', kmesh(ii)
              stop
          end select

          write(lun,
     &          fmt=fmt_string,
     &          err=9030)
     &            nodes,
     &            (lvcell(nodes_map(jj), ii), jj=1, nodes)
        end do


        write(lun,
     &        fmt="(/,A,X,I0)",
     &        iostat=rc,
     &        err=9030)
     &          VTK_CELL_TYPE,
     &          ncell

c dump cell type
        do ii = 1, ncell

          select case(kmesh(ii))
            case(1)
              nodes = VTK_TETRA
            case(2)
              nodes = VTK_HEXAHEDRON
            case(3)
              nodes = VTK_WEDGE
            case(4)
              nodes = VTK_PYRAMID
            case default
              write(*, *)'err> ', FNAME, ' unknown type: ', kmesh(ii)
              stop
          end select

          write(lun, fmt=fmt_data, iostat=rc, err=9030) nodes

          if (mod(ii, newline) .eq. 0) then
            write(lun, *, err=9030)
          end if

        end do

 9030   return
      end function write_vtk_cells

c**********************************************************************
      function write_vtk_boundary(lun) result (rc)
c**********************************************************************
c @brief  writes boundary data
c @note   none
c @param  lun [in] logical unit number of opened file
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
!#undef LABEL
!#define LABEL 9040

        use FFRdata, only : nbound,sfboun,nfboun,IBFACE,ifface,
     &                      nvrtx,NBFS

        implicit none

        integer, intent(in) :: lun
        integer :: rc

        character(len=*), parameter :: FNAME = "write_vtk_boundary: "
c       character(len=*), parameter :: fmt_data = "(I0,X)"
        character(len=*), parameter :: fmt_data = "(40(I0,X))"

        integer :: ii, jj, kk, node,ios=0,l,j,k
        integer, allocatable :: array(:)
        integer,allocatable :: ib(:)
!
! --- allocate
!
        allocate(array(1:nvrtx), stat=rc)
        allocate(ib(1:NBFS),stat=ios)
!
! --- 
!


        if (rc .ne. 0) then
          write(*, *)'err> ', FNAME, ' allocation'
          return
        end if

        write(lun,
     &        fmt="(/,/,A,X,I0)",
     &        iostat=rc,
     &        err=9040)
     &          VTK_POINT_DATA,
     &          nvrtx

        do ii=1,NBOUND
          l=0
          ib=0
          do j=1,NBFS
          if(IBFACE(j)==ii) then
            l=l+1
            ib(l)=j
          end if
          end do
          if(NFBOUN(ii)/=l) then
            write(*,*)'(append_resulten_VTK) ',
     &                 'Boundary faces is not match'
             return
          endif
          write(*, *)'dbg> ', FNAME, ' writing...', trim(sfboun(ii))
!
          write(lun,
     &          fmt="(/,A,X,A,X,A,/,A,X,A)",
     &          iostat=rc,
     &          err=9040)
     &            VTK_SCALAR,
     &            trim(sfboun(ii)),
     &            "char",
     &            VTK_LOOKUP_TABLE,
     &            "default"
!
          array = 0
          do kk=1,NFBOUN(ii)
          do jj=1,4
          node=IFFACE(jj,ib(kk))
          if(node/=0) then
            array(node)=1
          endif
          enddo
          enddo

          write(lun,
     &          fmt=fmt_data,
     &          iostat=rc,
     &          err=9040)
     &          array
        end do

        deallocate(array, stat=rc)
        if (rc .ne. 0) then
          write(*, *)'err> ', FNAME, ' de-allocation'
          return
        end if

 9040   return
      end function write_vtk_boundary

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function write_vtk_scalar(lun, byname, offset) result (rc)
c----------------------------------------------------------------------
c @brief  writes scalar data of each nodes
c @note   none
c @param  lun [in] logical unit number of opened file
c @param  byname [in] name of variable
c @param  offset [in] offset for totdata
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
!#undef LABEL
!#define LABEL 9050

        use FFRdata, only : nvrtx, totdata

        implicit none

        integer, intent(in) :: lun, offset
        character(len=*), intent(in) :: byname
        integer :: rc

        character(len=*), parameter :: FNAME = "write_vtk_scalar: "
c       character(len=*), parameter ::fmt_data = "(E)"
! Modified by Y.Takahashi, 2021.11.16 ---------------------
!        character(len=*), parameter ::fmt_data = "(9E)"
        character(len=*), parameter ::fmt_data = "(9E20.7)"
! ---------------------------------------------------------
        integer :: ii

        write(*, *)'dbg> ', FNAME, ' writing...', trim(byname)

        write(lun,
     &        fmt="(/,A,X,A,X,A)",
     &        iostat=rc,
     &        err=9050)
     &          VTK_SCALAR,
     &          trim(byname),
     &          "float"
        write(lun,
     &        fmt="(A,X,A)",
     &        iostat=rc,
     &        err=9050)
     &          VTK_LOOKUP_TABLE,
     &          "default"
        write(lun,
     &        fmt=fmt_data,
     &        iostat=rc,
     &        err=9050)
     &          (real(totdata(offset, ii)), ii=1, nvrtx)

 9050   return
      end function write_vtk_scalar

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function write_vtk_vector(lun, byname, offset) result (rc)
c----------------------------------------------------------------------
c @brief  writes scalar data of each nodes
c @note   none
c @param  lun [in] logical unit number of opened file
c @param  byname [in] name of variable
c @param  offset [in] offset for totdata
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
!#undef LABEL
!#define LABEL 9060

        use FFRdata, only : nvrtx, totdata

        implicit none

        integer, intent(in) :: lun, offset
        character(len=*), intent(in) :: byname
        integer :: rc

        character(len=*), parameter :: FNAME = "write_vtk_vector: "
c       character(len=*), parameter ::fmt_data = "(3E)"
! Modified by Y.Takahashi, 2021.11.16 ---------------------
!        character(len=*), parameter ::fmt_data = "(9E)"
        character(len=*), parameter ::fmt_data = "(9E20.7)"
! ---------------------------------------------------------
        integer :: ii

        write(*, *)'dbg> ', FNAME, ' writing...', trim(byname)

        write(lun,
     &        fmt="(/,A,X,A,X,A)",
     &        iostat=rc,
     &        err=9060)
     &          VTK_VECTOR,
     &          trim(byname),
     &          "float"
        write(lun,
     &        fmt=fmt_data,
     &        iostat=rc,
     &        err=9060)
     &          (real(totdata(offset, ii)),
     &           real(totdata(offset+1, ii)),
     &           real(totdata(offset+2, ii)),
     &           ii=1, nvrtx)

 9060   return
      end function write_vtk_vector


c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function ffr_to_vtk_ascii(fileName, gridOnly) result (rc)
c----------------------------------------------------------------------
c @brief  writes ascii vtk file
c @note   none
c @param  fileName [in] file name for output
c @param  gridOnly [in] flag to write grid data only or not
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
        use FFRdata, only : numVars, nameVars

        implicit none

        character(len=*), intent(in) :: fileName
        logical, intent(in) :: gridOnly
        integer :: rc

        integer :: lun, ii, pos

        rc = 1

        lun = open_vtk_file(fileName)

        if (lun .lt. 0) then
          return
        end if

        if (write_vtk_header(lun) .ne. 0) then
          return
        end if

        if (write_vtk_points(lun) .ne. 0) then
          return
        end if

        if (write_vtk_cells(lun) .ne. 0) then
          return
        end if

        if (write_vtk_boundary(lun) .ne. 0) then
          return
        end if

        if (.not. gridOnly) then
          ii = 1

          do while (ii <= numVars)

            pos = index(nameVars(ii), ';')

            if (pos .eq. 0) then
              rc = write_vtk_scalar(lun, nameVars(ii), ii)
              ii = ii + 1
            else
              rc = write_vtk_vector(lun,
     &                              nameVars(ii)(pos+1:),
     &                              ii)
              ii = ii + 3
            end if

            if (rc .ne. 0) then
              return
            end if

          end do

        end if

        if (close_vtk_file(lun) .lt. 0) then
          return
        end if

        rc = 0
        return
      end function ffr_to_vtk_ascii

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      function ffr_to_vtk_binary(fileName, gridOnly) result (rc)
c----------------------------------------------------------------------
c @brief  writes binary vtk file
c @note   none
c @param  fileName [in] file name for output
c @param  gridOnly [in] flag to write grid data only or not
c @param  rc [out] zero if successful, otherwise other than zero
c @return zero if successful, otherwise other than zero
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
        use FFRdata, only : nvrtx, cord,
     &                      ncell, kmesh, lvcell,
     &                      numVars, nameVars, totdata,
     &                      nbound, sfboun, nfboun, ifface,IBFACE,NBFS

c
c nvrtx:          number of nodes
c cord(3,nvrtx):  coordinate points of each node
c
c ncell:            number of cells
c kmesh(ncell):     array of cell type
c lvcell(8,ncell):  node connectivity of each cell (max. 8 nodes per cell)
c
c numVars:                number of variables
c nameVars(numVars):      name of variables
c totdata(numVars,nvrtx): data at each node
c NB: size of totdata changed from(:,1:nvrtx) to (:,0:nvrtx) at some point...
c
c nbound: number of boundary
c sfboun(nbound): name of boundary
c nfboun(nbound): number of faces in each boundary
c ifface(4, nfboun, nbound): list of nodes that belongs to faces which
c                            belongs to boundary, if ifface(4,i,j) = 0 ignore it
c

        implicit none

        character(len=*), intent(in) :: fileName
        logical, intent(in) :: gridOnly

        character(len=*), parameter :: FNAME = "ffr_to_vtk_binary: "
        integer :: ii, pos, header
        integer :: length, tuple, offset
        integer :: rc,l,j,k,i
        real(8), allocatable :: rarray(:,:)
        integer, allocatable :: iarray(:,:)
        integer,allocatable :: ib(:)
        
!
! --- 
!
        call write_vtk_header_bin(fileName)
        call write_vtk_points_bin(fileName, nvrtx, cord)
        call write_vtk_cells_bin(fileName, ncell, kmesh, lvcell)

c writes boundary data
        allocate(ib(1:NBFS))
        header=1
        do ii=1,NBOUND
          l=0
          ib=0
          do j=1,NBFS
          if(IBFACE(j)==ii) then
            l=l+1
            ib(l)=j
          end if
          end do
          allocate(iarray(1:4,1:nfboun(ii)), stat=rc)
          
          if (rc .ne. 0) then
            write(*, *)'err> ', FNAME, trim(fileName), ' allocation'
            return
          end if

          do k=1,NFBOUN(ii)
          do j=1,4
          iarray(j,k)=ifface(j,ib(k))
          enddo
          enddo
!          iarray =ifface(:,:,ii)

          call write_vtk_boundary_bin(fileName,
     &                                nfboun(ii),
     &                                iarray,
     &                                sfboun(ii),
     &                                nvrtx,
     &                                header)
          header = 0

          deallocate(iarray, stat=rc)
          if (rc .ne. 0) then
            write(*, *)'err> ', FNAME, ' de-allocation'
            return
          end if
        end do


c writes data at each node
        if (.not. gridOnly) then
          ii = 1

          tuple = size(totdata, 1)
          length = size(totdata) / tuple
c         write(*, *)'dbg> ', FNAME, buple, length, numVars, nvrtx

c need to ignore totdata(numVars,0)
          if (length .eq. nvrtx) then
            offset = 0
          else if (length .eq. nvrtx + 1) then
            offset = 1
          else
            write(*, *)'err> ', FNAME, ' todata: ', length, ' , ', nvrtx
            stop
          end if

          allocate(rarray(1:length, 1:tuple), stat=rc)
          if (rc .ne. 0) then
            write(*, *)'err> ', FNAME, trim(fileName), ' allocation'
            return
          end if

          rarray = transpose(totdata)

          do while (ii <= numVars)

            pos = index(nameVars(ii), ';')

            if (pos .eq. 0) then
              call write_vtk_scalar_bin(fileName,
     &                                  nvrtx,
     &                                  rarray(1,ii),
     &                                  offset,
     &                                  nameVars(ii))
              ii = ii + 1
            else
              call write_vtk_vector_bin(fileName,
     &                                  nvrtx,
     &                                  rarray(1,ii),
     &                                  offset,
     &                                  nameVars(ii))
              ii = ii + 3
            end if

          end do

          deallocate(rarray, stat=rc)
          if (rc .ne. 0) then
            write(*, *)'err> ', FNAME, ' de-allocation'
            return
          end if
        end if

        return

      end function ffr_to_vtk_binary

c----------------------------------------------------------------------c
c public functions/subroutines
c----------------------------------------------------------------------c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      subroutine ffr_to_vtk(fname, idum, CPUnum, gridOnly, rc, isAscii)
c----------------------------------------------------------------------
c @brief  writes vtk data
c @note   none
c @param  fname [in] file name for output
c @param  idum [in] dummy
c @param  CPUnum [in] dummy
c @param  gridOnly [in] flag to write grid data only or not
c @param  rc [out] zero if successful, otherwise other than zero
c @param  isAscii [in] .true. to write ascii data, .false. to write binary data
c @return none
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
        implicit none

        character(len=*), intent(in) :: fname
        integer, intent(in) :: idum, CPUnum
        logical, intent(in) :: gridOnly, isAscii
        integer ,intent(out) :: rc

        rc = 1

        if (isAscii) then
          if (ffr_to_vtk_ascii(fname, gridOnly) .ne. 0) then
            return
          end if
        else
          if (ffr_to_vtk_binary(fname, gridOnly) .ne. 0) then
            return
          end if
        end if

        rc = 0
        return
      end subroutine ffr_to_vtk

      end module tovtk
c eof
