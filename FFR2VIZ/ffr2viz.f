!======================================================================!
!                                                                      !
! Software Name : FFR2VIZ(FFRMOVIE)      (Ver. 1.0)                    !
!                                                                      !
! Main Program  : ffr2viz  (Multi Solver, HPC, Multi Domain)           !
!                                                                      !
!                          Written by Takeshi UNEMURA, Huilai ZHANG    !
!                                     Kazumi NISHIMURA                 !
!                                     2003/xx/xx                       !
!                                                                      !
! Software Name : FFR2VIZ(FFRMOVIE)      (Ver. 2.0)                    !
!                          Modified by Takeshi UNEMURA, Huilai ZHANG   !
!                                     Kazumi NISHIMURA, Keiji Ohnishi  !
!                                     2004/06/04                       !
!                                                                      !
!     Contact address : The University of Tokyo, FSIS Project          !
!                                                                      !
!======================================================================!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program ffr2viz
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      
      implicit none
!      
! --- 
!      
      integer :: iargc,argc
      integer :: i,ierror
      integer,parameter :: ifll=6,ifle=6
      character*56 :: argv
      character*56,allocatable :: arrayArg(:)
      
!
!
! --- Get Option index
!            
      argc=iargc()
      if(argc>=1) then
        allocate(arrayArg(1:argc))
        do 1010 i=1,argc
          call getarg(i,argv)
          arrayArg(i)=adjustl(argv)
 1010   continue
      end if

!
! --- Call main
!
      call main(argc,arrayArg,ierror)
      if( ierror.ne.0 ) goto 9999

      stop
!
! --- termination
!
 9999 continue
      write(ifle,*) '(ffr2viz)'
      write(ifle,*) '    ###    EXECUTION TERMINATED ABNORMALY'
      if( ifll.ne.ifle ) then
        write(ifll,*) '    ###    EXECUTION TERMINATED ABNORMALY'
      endif
      stop 9999
      end program ffr2viz

