      subroutine ftn_encode_header(elem_type, wall_info, header)
      use FVheaderTag, only : FV_TET_ELEM_ID,FV_HEX_ELEM_ID,
     &                        FV_PRISM_ELEM_ID,FV_PYRA_ELEM_ID,A_WALL
      implicit none
!  Arguments
      integer elem_type, wall_info(*), header

!  Local parameters (don t change these!)
      integer,parameter :: MAX_NUM_ELEM_FACES = 6
      integer,parameter :: BITS_PER_WALL = 3
      integer,parameter :: 
     &        ELEM_TYPE_BIT_SHIFT = MAX_NUM_ELEM_FACES*BITS_PER_WALL
!  Local variables
      integer i, nfaces

      if (elem_type .eq. FV_TET_ELEM_ID) then
         header = ishft (1, ELEM_TYPE_BIT_SHIFT)
         nfaces = 4
      else if (elem_type .eq. FV_HEX_ELEM_ID) then
         header = ishft (4, ELEM_TYPE_BIT_SHIFT)
         nfaces = 6
      else if (elem_type .eq. FV_PRISM_ELEM_ID) then
         header = ishft (3, ELEM_TYPE_BIT_SHIFT)
         nfaces = 5
      else if (elem_type .eq. FV_PYRA_ELEM_ID) then
         header = ishft (2, ELEM_TYPE_BIT_SHIFT)
         nfaces = 5
      else
         print *, 'ERROR:  Unknown element type'
         header = 0
         return
      endif

      do 3210 i = 1, nfaces
         if (wall_info(i) .lt. 0  .or.  wall_info(i) .gt. A_WALL) then
            print *, 'ERROR:  Bad wall value'
            header = 0
            return
         endif
         header = ior(header, ishft(wall_info(i), (i-1)*BITS_PER_WALL))
 3210 continue

      return
      end subroutine ftn_encode_header
!
