!
!      subroutine read_initial
!      subroutine read_source
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_initial
     & (NCV,NALLCV,NCVFAC,ncomp,nrans,
     & iters,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only     : lenfnm,ifle,getfil
      use module_model,only  : idrdp,mach0,comp,icaltb,rns_scl
      use module_rans,only   : akslw
      use module_dimnsn,only : xredc,yredc,zredc
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: NCV,NALLCV,NCVFAC,ncomp,nrans
      integer,intent(inout) :: iters
      integer,intent(out)   :: ierror
!
!
! -- [local entities]
!
      character(lenfnm) :: fnam
      integer :: l,n,ifl,ios,nrec,iter,itera
      integer :: ierr1
      integer :: NCVX,ncmpx,nrnsx,NCVFAX
      real*8  :: t0,y0(ncomp),p0,u0,v0,w0,aks0(nrans),vsgm0(2),rho0
      real*8  :: time,sum,vmx,sgm,r1,r2
!
!
!-< 1. Uniform flow field >-
!
      ierror=0
!
      if( iters.ge.0 ) goto 1000
!
      iters=0
!
      call nml_init
     &  (ncomp,nrans,t0,y0,p0,u0,v0,w0,rho0,vsgm0,aks0,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
 1000 continue
!
      return
!
 9999 continue
      write(ifle,'(a)') '(read_initial)'
!
      ierror=1
!
      end subroutine read_initial
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_source(mcell,ncell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only      : lenfnm,ifle,getfil,ifll
      use module_source,only  : nscnd
      USE module_usersub,ONLY : src_uvw,usrno,usryes
!
! 1. Input source file
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mcell,ncell
      integer,intent(out) :: ierror
!
!
! --- [local entities]
!
      character(lenfnm) :: fnam
      integer :: k,n,ifl,ios,ierr1
      integer :: no,kc,kcs,kce
      integer :: nrec,nsdmn,ncsdmn
      logical :: srcfil
!
!
!
      ierror=0
!
      if(src_uvw.eq.usryes) then
        write(ifll,'(2a)') 'MSG: User Subroutine [user/usrsrc.f] ',
     &  'is used for source terms by priority'
        return
      endif
!
      call getfil(ifl,fnam,'source')
      srcfil=.false.
      if( nscnd.gt.0 .and. fnam.eq.' ' ) then
        write(ifle,'(a)') ' ### error : data error'
        write(ifle,'(a)') 'lack of source file'
        goto 9999
      elseif( nscnd.lt.1 .and. fnam.ne.' ' ) then
        write(ifle,'(4a)') ' ### warning : ',
     &     'The source file ', trim(fnam),' is not used'
      elseif(nscnd.gt.0 .and. fnam.ne.' ' ) then
        inquire(file=fnam,exist=srcfil)
        if(.not.srcfil) then
          write(ifle,'(a)') ' ### error : file error'
          write(ifle,'(2a)') 
     &          'lack of source file',TRIM(adjustl(fnam))
          goto 9999
        endif
      elseif(fnam.eq.' ') then
        return
      endif
!
!
!-< 1. Input data >-
!
      return
!
 9999 continue
      write(ifle,'(a)') '(read_source)'
      ierror=1
!
      end subroutine read_source
!
