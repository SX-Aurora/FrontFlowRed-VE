!
!      subroutine nml_chkchr0
!      subroutine nml_chksiz
!      subroutine nml_errmsg0
!      subroutine nml_init
!      subroutine nml_listno
!      subroutine nml_sizes
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_chkchr0(ifle,cnam,cval,nn,ierr)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      integer     ,intent(in)  :: ifle
      character(*),intent(in)  :: cnam,cval
      integer     ,intent(in)  :: nn
      integer     ,intent(out) :: ierr
!
      ierr=0
      if( nn.gt.0 ) return
      write(ifle,*) '### error : data error'
      if( cval.eq.' ' ) then
        write(ifle,*) 'lack of data : ',cnam
      else
        write(ifle,*) cnam,' = ',cval
        if( nn.eq.0 ) then
          write(ifle,*) 'unknown word'
        else
          write(ifle,*) 'ambiguous word'
        endif
      endif
      write(ifle,*) '(nml_chkchr0)'
      ierr=1
      end subroutine nml_chkchr0
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_chksiz(ifle,cnam,jval,ival,mdnm,ierr)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      integer     ,intent(in)  :: ifle
      character(*),intent(in)  :: cnam
      integer     ,intent(in)  :: jval,ival
      character(*),intent(in)  :: mdnm
      integer     ,intent(out) :: ierr
!
      ierr=0
      if( ival.le.jval ) return
      write(ifle,*) '### error : parameter error'
      write(ifle,*) cnam,' =',jval
      write(ifle,*) 'it must be >=',ival
      write(ifle,*) 'reset ',cnam,' in ',mdnm
      write(ifle,*) '(nml_chksiz)'
      ierr=1
      end subroutine nml_chksiz
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_errmsg0(ifle,ios,nlst,ierr)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      integer     ,intent(in)  :: ifle,ios
      character(*),intent(in)  :: nlst
      integer     ,intent(out) :: ierr
      ierr=0
      if( ios.eq.0 ) return
      if( ios.gt.0 ) then
        write(ifle,*) '### error : read error'
      else
        write(ifle,*) '### error : end of file'
      endif
      write(ifle,*) 'namelist = ',nlst
      write(ifle,*) 'iostat =',ios
      write(ifle,*) '(nml_errmsg0)'
      ierr=1
      end subroutine nml_errmsg0

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_init
     &   (ncomp,nrans,t0,y0,p0,u0,v0,w0,rho0,vsgm0,aks0,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only : ifli,ifll,ifle
      use module_material,only : nofld,nflud
      use module_material,only : nosld,nsold
      use module_model,only    : idrdp,mach0,comp,incomp,icaltb,
     &                           rns_scl,sles,dles,lles
!
      implicit none
! 1. Input uniform initial condition
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ncomp,nrans
      real*8 ,intent(out) :: t0,y0(ncomp),p0,u0,v0,w0,vsgm0(2),rho0
      real*8 ,intent(out) :: aks0(nrans)
      integer,intent(out) :: ierror
!
! --- [namelist]
!
      integer,parameter :: mrans=10,mpotn=10
      integer,parameter :: mcomp=200

      integer :: IMAT_U,iflag,IMD,ICOM
      real*8  :: aks(mrans),potential(mpotn),potn_cond(mpotn)
      real*8  :: p,u,v,w,vsgm(2),t,ys(mcomp),dens
      namelist /initial/ IMAT_U,t,p,u,v,w,vsgm,ys,dens,aks,potential,
     &                   potn_cond
!
! --- [local entities]
!
      real*8, parameter :: undef=-huge(1.d0)
      integer,parameter :: iundef=-huge(1)
      real*8  :: sum
      integer :: i,j,ios,ierr1=0,nfld,nsld,nfsld,noinit(nflud+nsold)
!
!
!
      ierror=0
      nfld=0
      nsld=0
      nfsld=0
      noinit(:)=0
!
      call nml_chksiz(ifle,'mcomp',mcomp,ncomp,
     & '(nml_init)',ierr1)
      if( ierr1.ne.0 ) goto 9999
      call nml_chksiz(ifle,'mrans',mrans,nrans,
     & 'nml_init',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      rewind ifli
 350  continue
      IMAT_U=iundef
      t   =undef
      ys  =undef
      p   =undef
      u   =undef
      v   =undef
      w   =undef
      dens=undef
      vsgm=undef
      aks=undef
!
      read(ifli,initial,iostat=ios)
      if( ios.lt.0 ) goto 351
      call nml_errmsg0(ifle,ios,'initial',ierr1)
      if( ierr1.ne.0 ) then
        write(ifle,*) '### ERR: reading error in namelist/initial 1'
        goto 9999
      endif
!
      if(IMAT_U.eq.iundef) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : IMAT_U'
        goto 9999
      endif
!
      if(IMAT_U.gt.0) then
        do i=1,nflud
          if(IMAT_U.eq.nofld(i)) then
            nfld=nfld+1
          endif
        enddo
      elseif(IMAT_U.lt.0) then
        do i=1,nsold
          if(IMAT_U.eq.nosld(i)) then
            nsld=nsld+1
          endif
        enddo
      elseif(IMAT_U.eq.0) then
        write(ifle,*) '### ERR: IMAT_U=0 is NOT supported'
        write(ifle,*) '### IMAT_U > 0 or IMAT_U < 0'
        goto 9999
      endif
      nfsld=nfsld+1
!      if(nfsld.gt.nflud+nsold) then
!        write(ifle,*) '### ERR: too many set of namelist &initial'
!        goto 9999
!      endif
      noinit(nfsld)=IMAT_U
      do i=1,nfsld-1
      if(IMAT_U.eq.noinit(i))then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'this IMAT_U have been defined'
        goto 9999
      endif
      enddo
!
      if( t.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : t'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      elseif( t.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 't = ',t
        write(ifle,*) 'it must be > 0'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
!
      if(IMAT_U.lt.0) then
        write(ifll,*) ' ### ATT: Solid material'
        write(ifll,*) ' [&initial: ys,p,u,v,w,vsgm,dens,aks]',
     &                        'are NEGLECTED'
        goto 360
      endif
!
      if( dens.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : dens'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      elseif( dens.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'dens = ',dens
        write(ifle,*) 'it must be > 0'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
!
      sum=0.d0
      do 100 ICOM=1,ncomp
      sum=sum+ys(ICOM)
      if( ys(ICOM).eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : ys(ICOM)'
        write(ifle,*) 'ICOM =',ICOM
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      elseif( ys(ICOM).lt.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'ys(ICOM) = ',ys(ICOM)
        write(ifle,*) 'ICOM =',ICOM
        write(ifle,*) 'it must be >= 0'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
  100 continue
      if( sum.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'ys = ',(ys(i),i=1,ncomp)
        write(ifle,*) 'sum of them must be > 0'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
!
      if( p.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : p'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      elseif( p.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'p = ',p
        write(ifle,*) 'it must be > 0'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
!
!
      if( u.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : u'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
      if( v.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : v'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
      if( w.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : w'
        write(ifle,*) '[&initial] in fflow.ctl'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        goto 9999
      endif
!
      if( vsgm(1).eq.undef ) then
        vsgm=0.d0
      else
        if( vsgm(1).le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'vsgm(1) = ',vsgm(1)
          write(ifle,*) 'it must be > 0'
          write(ifle,*) '[&initial] in fflow.ctl'
          write(ifle,*) 'IMAT_U= ',IMAT_U
          goto 9999
        endif
        if( vsgm(2).eq.undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : vsgm(2)'
          write(ifle,*) '[&initial] in fflow.ctl'
          write(ifle,*) 'IMAT_U= ',IMAT_U
          goto 9999
        elseif( vsgm(2).le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'vsgm(2) = ',vsgm(2)
          write(ifle,*) 'it must be > 0'
          write(ifle,*) '[&initial] in fflow.ctl'
          write(ifle,*) 'IMAT_U= ',IMAT_U
          goto 9999
        endif
      endif
      if(rns_scl) then
        do IMD=1,nrans
        if( aks(IMD).eq.undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : aks(IMD) '
          write(ifle,*) '[&initial] in fflow.ctl'
          write(ifle,*) 'IMAT_U= ',IMAT_U
          write(ifle,*) 'IMD =',IMD
          goto 9999
        elseif( aks(IMD).lt.0. ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'aks(IMD) = ',aks(IMD)
          write(ifle,*) 'it must be >= 0'
          write(ifle,*) '[&initial] in fflow.ctl'
          write(ifle,*) 'IMAT_U= ',IMAT_U
          write(ifle,*) 'IMD =',IMD
          
          goto 9999
        endif
        enddo
      else
        aks=0.d0
      endif
!
 360  continue
!
      t0=t
      p0=p
      u0=u
      v0=v
      w0=w
      y0(1:ncomp)=ys(1:ncomp)
      vsgm0(1)=vsgm(1)
      vsgm0(2)=vsgm(2)
      if(rns_scl) then
        aks0(1:nrans)=aks(1:nrans)
      endif
      goto 350
 351  continue
!
      do j=1,nflud
      iflag=0
      do i=1,nfsld
      IMAT_U=noinit(i)
      if(nofld(j).eq.IMAT_U) then
        iflag=iflag+1
      endif
      enddo
      if(iflag.eq.0) then
          write(ifle,*) '### ERR-1: IMAT_U=',nofld(j),
     &                  ' Initial condition is NOT specified in',
     &                  ' namelist/initial'
        goto 9999
      endif
      enddo
!
      do j=1,nsold
      iflag=0
      do i=1,nfsld
      IMAT_U=noinit(i)
      if(IMAT_U.eq.nosld(j)) then
        iflag=iflag+1
      endif
      enddo
      if(iflag.eq.0) then
          write(ifle,*) '### ERR-2: IMAT_U=',nosld(j),
     &                  ' Initial condition is NOT specied in',
     &                  ' namelist/initial'
        goto 9999
      endif
      enddo
!
      if(nfld.lt.1.and.nsld.lt.1) then
        write(ifle,*) '### ERR: namelist/initial must be specified 2'
        goto 9999
      elseif(nfld.lt.1) then
        write(ifle,*) '### ERR: NO fluid field exists'
        write(ifle,*) '### ERR: This case is NOT supported yet'
!MHD_err        goto 9999
      elseif(nfld+nsld.ne.nflud+nsold) then
        write(ifle,*) 
     &            '### ERR-28: no. of namelist-initial NOT equal to',
     &                'namelist-fluid + namelist-solid'
        goto 9999
      endif
!
      return
!
 9999 continue
      write(ifle,*) '(nml_init)'
      ierror=1
      end subroutine nml_init
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_listno(nl,lst,wrd,nn)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
! 1. Return list no.
!      nn : =-1; ambiguous word
!         : = 0; unknown word
!         : > 0; identified word
!
! --- [dummy arguments]
!
      integer     ,intent(in)  :: nl
! for debug use -- by onishi
!      character, dimension(:), intent(in)  :: lst
      character(*),intent(in)  :: lst(nl)
      character(*),intent(in)  :: wrd
      integer     ,intent(out) :: nn
!
! --- [local entities]
!
      integer :: i,n,lng,ic
!
! add F95 trim functionality -- by onishi
!     on UNIX,Linux ... no need (right side of char is not filled with 'space')
!     on Windows    ... needed  (right side of char is filled with 'space')
!
      nn=0
      lng=min(len_trim(wrd),len_trim(lst(1)))
      if( wrd(:lng).eq.' ' ) return
      ic=0
      do 100 n=1,nl
      if( lst(n)(:lng).eq.' ' ) goto 100
      do 101 i=1,lng
      if( wrd(i:i).eq.' ' ) goto 102
      if( wrd(i:i).ne.lst(n)(i:i) ) goto 100
  101 continue
  102 continue
      ic=ic+1
      nn=n
  100 continue
      if( ic.gt.1 ) nn=-1
!
      end subroutine nml_listno
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_sizes(ifli,ifll,ifle,cntlnam,
     & mcello,mvrtxo,mfaceo,medgeo,mssfbco,
     & ncello,nvrtxo,ncompo,ncomp_sufo,lFC,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
! 1. Input sizes
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle,lFC
      character,intent(in):: cntlnam
      integer,intent(out) :: mcello,mfaceo,mvrtxo,medgeo,mssfbco
      integer,intent(out) :: ncello,nvrtxo,ncompo,ncomp_sufo
      integer,intent(inout) :: ierror
!
!
! --- [namelist]
!
      integer :: mcell,mface,ncell,nvrtx,ncomp,mvrtx,mssfbc,medge
      integer :: ncomp_surface,npotential,MXPRT
      namelist /sizes/ mcell,mface,medge,ncell,nvrtx,ncomp,
     &                 mvrtx,mssfbc,ncomp_surface,MXPRT!,npotential
!       mcell : upper limit of no. of cells
!       mface : upper limit of no. of faces
!       ncell : no. of cells
!       nvrtx : no. of vertices
!       ncomp : no. of components in mixed gas
!       nrans : no. of components in RANS model
!       ncomp_surface : no. of components in surface reaction
!       npotential    : no of potential number
! --- [local entities]
!
      integer,parameter :: iundef=-huge(1)
      integer :: i,ios
!
!
!
!-< 1. Input namelist >-
!
      mcell=iundef
      mface=iundef
      medge=iundef
      mvrtx=iundef
      mssfbc=iundef
      ncell=iundef
      nvrtx=iundef
      ncomp=iundef
      ncomp_surface=0
      npotential=0
!
      rewind ifli
      read (ifli,sizes,iostat=ios)
!      write(ifll,sizes)
      call nml_errmsg0(ifle,ios,'sizes',ierror)
      if( ierror.ne.0 ) goto 9999
!
      call errmsg1('mcell',mcell,ncell)
      call errmsg1('mface',mface,1)
      call errmsg1('medge',medge,1)
      call errmsg1('ncell',ncell,1)
      call errmsg1('mvrtx',mvrtx,nvrtx)
      call errmsg1('mssfbc',mssfbc,1)
      call errmsg1('ncomp',ncomp,1)
!      call errmsg1('ncomp_surface',ncomp,1)
!      call errmsg1('nrans',nrans,2)
      if( ierror.ne.0 ) goto 9999
!
      mcello=mcell
      mfaceo=mface
      medgeo=medge
      ncello=ncell
      nvrtxo=nvrtx
      ncompo=ncomp
      ncomp_sufo=ncomp_surface
!      npotno=npotential
!      if(lpotn) then
!        if(npotno==0) then
!          write(*,*) 'ERR: npotno=0'
!          write(*,*) 'MSG: Potential/=0'
!          stop 'MSG: stop at /sizes/'
!        endif
!      else
!        if(npotno/=0) then
!          write(*,*) 'MSG: Potential=0'
!          write(*,*) 'ERR: npotential=',npotno
!          stop 'MSG: set &model/Potential=1'
!        endif
!      endif
!      if(lFC>0) then
!        if(.not.lpotn) then
!          write(*,*) 
!     &  'ERR: &model/Potential/=1, NOT correct for  PEFC'
!          write(*,*) 'MSG: &model/Potential=1 for PEFC'
!          stop
!        endif
!        if(npotno<2) then
!          write(*,*) 
!     &  'ERR: npotential<2, NOT correct for  PEFC'
!          write(*,*) 'MSG: npotential>=2 for PEFC'
!          stop 
!        endif
!      endif
!
      mvrtxo=mvrtx
      mssfbco=mssfbc
!
      return
!
 9999 continue
      write(ifle,*) '(nml_sizes)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!==================================================
      subroutine errmsg1(vnam,idat,ib)
!==================================================
      implicit none
      character(*) :: vnam
      integer      :: idat,ib
      if(idat.eq.iundef) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : ',vnam
        ierror=1
        return
      endif
      if(idat.lt.ib) then
        write(ifle,*) '### error : data error'
        write(ifle,*) vnam,' =',idat
        write(ifle,*) 'it must be >=',ib
        ierror=1
      endif
      return
      end subroutine errmsg1
!
      end subroutine nml_sizes
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      logical function nml_comp_eq(str,v1,v2,s,f,e)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- compare values (real*8), equal (true) or not (false)
!---------------------------------------------------------------------
      implicit none
      character(*),intent(in) :: str    ! string of value
      real*8,intent(in) :: v1,v2        ! values to compare
      integer,intent(in) :: s           ! sequential number
      integer,intent(in) :: f           ! file number to write error script
      integer,intent(inout) :: e        ! total number of error
      nml_comp_eq = .false.             ! not equal
      if( v1==v2 ) then
        write(f,*) "### error : data error"
        write(f,*) "lack of data '",str,"' in no.",s
        e = 1+e
        nml_comp_eq = .true.            ! equal
      endif
      end function nml_comp_eq
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      logical function nml_comp_gt(str,v1,v2,s,f,e)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- compare values (real*8), v1 > v2 (true) or not (false)
!---------------------------------------------------------------------
      implicit none
      character(*),intent(in) :: str    ! string of value
      real*8,intent(in) :: v1,v2        ! values to compare
      integer,intent(in) :: s           ! sequential number
      integer,intent(in) :: f           ! file number to write error script
      integer,intent(inout) :: e        ! total number of error
      nml_comp_gt = .true.              ! v1 > v2
      if( v1<=v2 ) then
        write(f,*) "### error : data error"
        write(f,*) "'",str,"' =",v1," in no.",s
        write(f,*) "it must be > ",v2
        e = 1+e
        nml_comp_gt = .false.           ! v1 <= v2
      endif
      end function nml_comp_gt
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      logical function nml_comp_ge(str,v1,v2,s,f,e)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! ---- compare values (real*8), v1 >= v2 (true) or not (false)
!---------------------------------------------------------------------
      implicit none
      character(*),intent(in) :: str    ! string of value
      real*8,intent(in) :: v1,v2        ! values to compare
      integer,intent(in) :: s           ! sequential number
      integer,intent(in) :: f           ! file number to write error script
      integer,intent(inout) :: e        ! total number of error
      nml_comp_ge = .true.              ! v1 >= v2
      if( v1<v2 ) then
        write(f,*) "### error : data error"
        write(f,*) "'",str,"' =",v1," in no.",s
        write(f,*) "it must be >= ",v2
        e = 1+e
        nml_comp_ge = .false.           ! v1 < v2
      endif
      end function nml_comp_ge
!

