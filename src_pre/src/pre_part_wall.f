!
!
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine WALL(mvrtx,cord)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_boundary,only : kdbcnd,kxnone,nbcnd,ivbcnd,lvbcnd
      use module_io,only       : gdformat
      use module_partitioner
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mvrtx
      real*8 ,intent(inout)  :: cord  (3,mvrtx)
!
! --- [local entities]
!
      REAL*8  :: X0,Y0,Z0,DL,DLmin=1.d+24,XP,YP,ZP
      integer :: i,j,k,iw,iv,ipmin,ivW,iwmin,NWALL
      integer :: ips,ipe,ivl,nb,ierr1=0
!
! --- 
!
      allocate (WALLnode(INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating WALLnode(:) in WALL'
      allocate (WALLDIST(INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating WALLDIST(:) in WALL'
!
      if (WALLFLAG.eq.1) then
        open(21,file=WALFIL,status='unknown',form  ='unformatted')
        read (21) NWALL
        write (*,*) NWALL,INODTOT
        if (NWALL.ne.INODTOT) call ERROR_EXIT(5000,NWALL)
        read (21) (WALLnode(i),i=1,INODTOT)
        close(21)
        return
      endif
!
! --- Call shortest distance to wall from all vertex
!
      do 100 iv=1,INODTOT
      X0=cord(1,iv)
      Y0=cord(2,iv)
      Z0=cord(3,iv)
!
      if (mod(iv,10000).eq.0) then
        write (*,'(a, 2i12)') '    wall ', iv, INODTOT
      endif
!
! --- +------+
! --- | WALL | the shortest distance :
! --- +------+ 
!
! --- BC_WALL_tot: total wall vertex number
!
      do 140 iw=1,BC_WALL_tot
      ivW=BC_WALL(iw)
      XP=cord(1,ivW)
      YP=cord(2,ivW)
      ZP=cord(3,ivW)
      DL= dsqrt ((XP-X0)**2+(YP-Y0)**2+(ZP-Z0)**2)
      if(DL.lt.DLmin) then
        iwmin= ivW
        DLmin= DL
      endif
 140  continue
!
! --- iwmin: vertex number on wall 
! --- of the shortest distance
!
      WALLnode(iv)=iwmin
      WALLDIST(iv)=DLmin
!
 100  continue
!
      deallocate(WALLnode)
      deallocate(WALLDIST)
!
      return
!
      end subroutine WALL
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine ERROR_EXIT (IFLAG, nn)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use  module_partitioner
      implicit none
      integer     ,intent(in)  :: IFLAG, nn
      
      write (*,'(/,a)')                                                 
     &        "********** MESSAGE from Partitioner **********"
      if (IFLAG.ge. 1001 .and. IFLAG.lt.2000) then
        write (*,'(/,a)')                                               
     &        " ### ABORT : unexpected ZERO/minus in the orginal file"
        if (IFLAG.eq.1001) write (*,'(  a,/)')                          
     &        "     TOTAL NODE and/or ELEMENT NUMBER"
        if (IFLAG.eq.1002) write (*,'(  a,i8/)')                        
     &        "     BOUNDARY GROUP NUMBER (1:node, 2:elem, 3:suf)", nn
        if (IFLAG.eq.1003) write (*,'(  a,i8/)')                        
     &        "     BOUNDARY info ITEMs   (1:node, 2:elem, 3:suf)", nn
        if (IFLAG.eq.1004) write (*,'(  a,i8/)')                        
     &        "     ELEMENT type", nn
        if (IFLAG.eq.1005) write (*,'(  a,i8/)')                        
     &        "     ELEMENT connectivity in ", nn
        stop
      endif

! Erase last line feed '&', instead use -132 option -- by onishi
!      if (IFLAG.eq. 2001) then
!        write (*,'(/,a,i8/)')                                        &   
!     &        " ### ABORT : local node ID > N appears in ELEMENT", nn
!        stop
!      endif
!      ...

      if (IFLAG.eq. 2001) then
        write (*,'(/,a,i8/)')                                           
     &        " ### ABORT : local node ID > N appears in ELEMENT", nn
        stop
      endif

      if (IFLAG.eq.2002) then
        write (*,'(/,a  )')                                             
     &        " ### ABORT : local node ID > N appears in BOUNDARY"
        stop
      endif

      if (IFLAG.eq.5000) then
        write (*,'(/,a  )')                                             
     &        " ### ABORT : INVALID wall-law file"
        stop
      endif

      if (IFLAG.eq. 2) then
        write (*,'(/,a,/)')                                             
     &        " ### ABORT : ERROR in ORIGINAL GRID FILE : Parallel Info"
        stop
      endif

      if (IFLAG.eq.11) then
        write (*,'(/,a,/)')                                             
     &        " ### ABORT : ERROR in ORIGINAL GRID FILE"
        stop
      endif

      if (IFLAG.eq.12) then
        write (*,'(/,a,/)')                                             
     &        " ### ABORT : UNEXPECTED EOF in GRID/MeTiS FILE"
        stop
      endif

      if (IFLAG.eq.21) then
        write (*,'(/,a,/)')                                             
     &        " ### ABORT : ERROR in GRID/MeTiS FILE"
        stop
      endif

      if (IFLAG.eq.22) then
        write (*,'(/,a,/)')                                             
     &        " ### ABORT : UNEXPECTED EOF in GRID/MeTiS FILE"
        stop
      endif

      if (IFLAG.eq.31) then
        write (*,'(/,a,/)')                                             
     &        " ### ABORT : MeTiS file INCONSISTENCY"
        stop
      endif

      if (IFLAG.eq.32) then
        write (*,'(/,a,2i8/)')                                          
     &        " ### ABORT : INVALID PE  and node #", PETOT, nn
        stop
      endif

      if (IFLAG.eq.33) then
        write (*,'(/,a,i8/)')                                           
     &        " ### ABORT : INVALID element type", nn
        stop
      endif

      if (IFLAG.eq.5001) then
        write (*,'(/,a,i8/)')                                           
     &        " ### ABORT : UNSUPPORTED element type", nn
        stop
      endif

      if (IFLAG.eq.8001) then
        write (*,'(/,a,/)')                                             
     &        " ### ABORT : ERROR in CONTROL FILE"
        stop
      endif

      end
