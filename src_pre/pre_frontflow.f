!======================================================================!
!                                                                      !
! Software Name :                                                      !
!     FrontFlow/red   (Ver. 1.0) (Cell-Center,FACE-BASE)               !
!                                                                      !
!     Main Program : um0_FFRMAIN                                       !
!                                                                      !
!                          Written by Masayuki KAKEI,                  !
!                                     Huilai ZHANG, Osamu KITAMURA     !
!                                     Yoshinobu YAMADE, Masato IDA,    !
!                                     Takeshi UNEMURA, Eisuke YAMADA,  !
!                                     Nobuyuki TANIGUCHI               !
!                                     2003/03/25                       !
!                                                                      !
!     FontFlow/red    (Ver. 1.2) (EDGE-BASE,HPC)                       !
!                          Modified by                                 !
!                                     Huilai ZHANG                     !
!                                     Takeshi UNEMURA                  !
!                                     NAKAJIMA???                      !
!                                                                      !
!     Contact address : The University of Tokyo, FSIS Project          !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 2.0)                           !
!                                                                      !
! Main Program  : FRONTFLOW    (Vertex-Center,EDGE-BASE,HPC,           !
!                                  Multi-Domain)                       !
!                                                                      !
!                          Written by Huilai ZHANG                     !
!                                     Takeshi UNEMURA                  !
!                                     Nobuyuki TANIGUCHI               !
!                                     2003/03/25                       !
!                                                                      !
!     Contact address : The University of Tokyo, FSIS Project          !
!                                                                      !
!======================================================================!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program pre_FRONTFLOW
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!      USE MODULE_IO       ,ONLY : IFLI,IFLL,IFLE
!
!     Pre-treatment program for FFLOW
!
      implicit none
      integer :: mvrtx,mcell,mssfbc,mface,medge
      integer :: nvrtx,ncell
      integer :: ncomp,nrans,ncomp_suf,nphase,npotn
      integer :: nvrtxnw
      integer :: ierror
      integer :: IFLE=6,IFLL=6
!
      REAL  TIME_BEGAN, TIME_END   !jiang
!
!
!
! --- dimension size initialization
!
      mface =0
      medge =0
      mcell=0
      nvrtxnw=0
      mvrtx =0
      mssfbc=0
      ncell =0
      nvrtx =0
      ncomp =1
      nrans =1
      ncomp_suf=0
      nphase=ncomp
      npotn=0
      ierror=0
!-----------------------
! --- Title 
!-----------------------
      write (ifll,*)
      write (ifll,5000)
      write (ifll,5010)
! DEBUG
      write (ifll,5020)
      write (ifll,5030)
!
      write (ifll,5000)
!-----------------------
! --- read namelist file
!-----------------------
      call pre_namelist_admin
     & (mcell,mvrtx,mface,medge,mssfbc,ncell,nvrtx,ncomp,nrans,
     &  ncomp_suf,nphase,npotn,ierror)
      if( ierror.ne.0 ) goto 9999
!-----------------------
! --- main treatment
!-----------------------
      call pre_main
     & (mcell,mvrtx,mface,medge,mssfbc,
     &  ncell,nvrtx,nvrtxnw,ncomp,ncomp_suf,nphase,nrans,npotn,ierror)
      if( ierror.ne.0 ) goto 9999
!
! --- termination
!
      write(ifll,'(5x,a)') '%%%    Pre-Exec  terminated normally'
      write(ifll,5000)
 5000 format('|',108('='),'|')
 5010 format('|',46(' '),' Welcome to FrontFlow ',40(' '),'|')
 5020 format(' ',46(' '),' VERSION    = src_073.2')
 5030 format(' ',46(' '),' BUILD DATE = 2008-09-19-05:19:04')
      stop
 9999 continue
      write(ifle,'(a)') '(pre_FRONTFLOW)'
      write(ifle,'(a)') '    ###    EXECUTION TERMINATED ABNORMALLY'
      if( ifll.ne.ifle ) then
        write(ifll,'(a)') '    ###    EXECUTION TERMINATED ABNORMALLY'
      endif
      stop 9999
      end program pre_FRONTFLOW
