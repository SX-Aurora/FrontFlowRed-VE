      subroutine writeFLUENTHeader(fname_root,ifl,gridOnly,status)
!
      use FLUENTdata, only : F_VERSION
      use FFRdata, only : nvrtx,NFCE,ncell
!
      implicit none
!
      character(80),intent(in) :: fname_root
      integer,intent(in) :: ifl
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
      
      character(8) :: cdate
      integer :: i,ios
      character(80) :: fname
!
! --- -----------------------------------------------------
!     open FLUENT/cas format file and write header
! --- -----------------------------------------------------
      fname=trim(adjustl(fname_root)) // '.cas'
!      write(fname,*)trim(adjustl(fname_root)),'.cas'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeFLUENTHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeFLUENTHeader)Cannot create ',trim(fname)
        status=1; return
      end if
!
      write(ifl,112) '(0 "',F_VERSION,'")'
      write(ifl,110) ''
      write(ifl,110) '(0 "Machine Config:")'
      write(ifl,110) '(4 (60 0 0 1 2 4 4 4 8 8 4))'
      write(ifl,110) ''
      write(ifl,110) '(0 "Dimensions:")'
      write(ifl,110) '(2 3)'
      write(ifl,110) ''
      write(ifl,110) '(0 "Variables:")'
      write(ifl,110) '(37 ('
      write(ifl,111) '(case-config ('
      write(ifl,111) '(rp-acoustics? . #f)'
      write(ifl,111) '(rp-atm? . #f)'
      write(ifl,111) '(rp-axi? . #f)'
      write(ifl,111) '(rp-dpm-cache? . #f)'
      write(ifl,111) '(rp-unsteady? . #f)'
      write(ifl,111) '(rp-dual-time? . #f)'
      write(ifl,111) '(rp-amg? . #f)'
      write(ifl,111) '(rf-energy? . #f)'
      write(ifl,111) '(rp-hvac? . #f)'
      write(ifl,111) '(rp-inviscid? . #f)'
      write(ifl,111) '(rp-ke? . #t)'
      write(ifl,111) '(rp-kw? . #f)'
      write(ifl,111) '(rp-lam? . #f)'
      write(ifl,111) '(rp-les? . #f)'
      write(ifl,111) '(rp-lsf? . #f)'
      write(ifl,111) '(rp-net? . #f)'
      write(ifl,111) '(rp-react? . #f)'
      write(ifl,111) '(rp-sa? . #f)'
      write(ifl,111) '(rp-sa-des? . #f)'
      write(ifl,111) '(rp-seg? . #t)'
      write(ifl,111) '(rp-sge? . #f)'
      write(ifl,111) '(rp-spe? . #f)'
      write(ifl,111) '(rp-spe-part? . #f)'
      write(ifl,111) '(rp-spe-site? . #f)'
      write(ifl,111) '(rp-spe-surf? . #f)'
      write(ifl,111) '(rp-trb-scl? . #t)'
      write(ifl,111) '(rp-turb? . #t)'
      write(ifl,111) '(rp-absorbing-media? . #f)'
      write(ifl,111) '(rp-visc? . #t)'
      write(ifl,111) '(rp-v2f? . #f)'
      write(ifl,111) '(sg-cylindrical? . #f)'
      write(ifl,111) '(sg-disco? . #f)'
      write(ifl,111) '(sg-bee-gees? . #f)'
      write(ifl,111) '(sg-crev? . #f)'
      write(ifl,111) '(sg-dpm? . #f)'
      write(ifl,111) '(sg-dtrm? . #f)'
      write(ifl,111) '(sg-dynmesh? . #f)'
      write(ifl,111) '(sg-network? . #f)'
      write(ifl,111) '(sg-melt? . #f)'
      write(ifl,111) '(sg-mphase? . #f)'
      write(ifl,111) '(sg-p1? . #f)'
      write(ifl,111) '(sg-par-premix? . #f)'
      write(ifl,111) '(sg-pdf? . #f)'
      write(ifl,111) '(sg-pdf-transport? . #f)'
      write(ifl,111) '(sg-pollut? . #f)'
      write(ifl,111) '(sg-premixed? . #f)'
      write(ifl,111) '(sg-pull? . #f)'
      write(ifl,111) '(sg-rosseland? . #f)'
      write(ifl,111) '(sg-rsm? . #f)'
      write(ifl,111) '(sg-s2s? . #f)'
      write(ifl,111) '(sg-soot? . #f)'
      write(ifl,111) '(sg-swirl? . #f)'
      write(ifl,111) '(sg-udm? . #f)'
      write(ifl,111) '(sg-uds? . #f)'
      write(ifl,111) '(sg-vfr? . #f)'
      write(ifl,111) '(solar? . #f)'
      write(ifl,111) '(rp-3d? . #t)'
      write(ifl,111) '(rp-double? . #t)'
      write(ifl,111) '(rp-graphics? . #t)'
      write(ifl,111) '(rp-host? . #f)'
      write(ifl,111) '(rp-thread? . #f)'
      write(ifl,111) '(dpm-cache? . #t)))))'
      write(ifl,110) ''
      write(ifl,110) ''

      write(ifl,110) '(0 "Domain Variables:")'
      write(ifl,110) '(64 ())'
      write(ifl,110) ''

      close(ifl)

 110  format(A)
 111  format(A,1X,$)
 112  format(A,A,A)

! --- -----------------------------------------------------
!     open FLUENT/dat format file and write header
! --- -----------------------------------------------------
      if (.not.(gridOnly)) then

      fname=trim(adjustl(fname_root)) // '.dat'
!      write(fname,*)trim(adjustl(fname_root)),'.dat'
      inquire(file=fname,number=ios)
      if(ios>=0) then
        write(*,*)'(writeFLUENTHeader)File ',trim(fname),
     &            ' has been already opened as number ',ios
        status=1; return
      end if
!
      open(ifl,file=fname,form='formatted',action='write',iostat=ios)
      if(ios/=0) then
        write(*,*)'(writeFLUENTHeader)Cannot create ',trim(fname)
        status=1; return
      end if
!
      write(ifl,122) '(0 "',F_VERSION,'")'
      write(ifl,120) ''
      write(ifl,120) '(0 "Machine Config:")'
      write(ifl,120) '(4 (60 0 0 1 2 4 4 4 8 8 4))'
      write(ifl,120) ''
      write(ifl,120) '(0 "Grid size:")'
      write(ifl,121) '(33 (',ncell,NFCE,nvrtx,'))'
      write(ifl,120) ''
      write(ifl,120) '(0 "Variables:")'
      write(ifl,120) '(37 ('
      write(ifl,120) '(flow-time 0)'
      write(ifl,120) '(time-step 0)'
      write(ifl,120) '(periodic/pressure-derivative 0)'
      write(ifl,120) '(number-of-samples 0)'
      write(ifl,120) '(dpm/summary ())'
      write(ifl,120) '(operating-pressure 101325)'
      write(ifl,120) '(delta-time-sampled 0)))'
      write(ifl,120) ''

      close(ifl)

      end if

 120  format(A)
! Modified by Y.Takahashi, 2021.11.16
! 121  format(A,3I,A)
 121  format(A,3I8,A)
 122  format(A,A,A)

      status=0
!
      return
      end subroutine writeFLUENTHeader
