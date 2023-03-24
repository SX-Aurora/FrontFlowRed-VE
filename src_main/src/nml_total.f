!
!      subroutine nml_chkchr0
!      subroutine nml_chksiz
!      subroutine nml_errmsg0
!      subroutine nml_init
!      subroutine nml_listno
!      subroutine nml_sizesCHK
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_chkchr0(ifle,cnam,cval,nn,ierr)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
!
! --- [dummy arguments]
!
      integer     ,intent(in)  :: ifle
      character(*),intent(in)  :: cnam,cval
      integer     ,intent(in)  :: nn
      integer     ,intent(out) :: ierr
!
! --- [local entities]
!
      ierr=0
      if( nn.gt.0 ) return
      write(ifle,'(1X,A)') 'ERR: : data error'
      if( cval.eq.' ' ) then
        write(ifle,'(1X,2a)') 'MSG: lack of data : ',cnam
      else
        write(ifle,'(1X,4a)') 'ERR: ',cnam,' = ',cval
        if( nn.eq.0 ) then 
          write(ifle,'(1x,a)') 'ERR: unknown word'
        else
          write(ifle,*) 'ambiguous word'
        endif
      endif
      write(ifle,*) '(nml_chkchr0)'
      ierr=1
      end subroutine nml_chkchr0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_chksiz(ifle,cnam,jval,ival,mdnm,ierr)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
! --- [dummy arguments]
!
      integer     ,intent(in)  :: ifle
      character(*),intent(in)  :: cnam
      integer     ,intent(in)  :: jval,ival
      character(*),intent(in)  :: mdnm
      integer     ,intent(out) :: ierr
!
      ierr=0
      if( ival.le.jval ) return
      write(ifle,'(1X,A)') 'ERR: parameter error'
      write(ifle,'(1x,3a,I8)') 'MSG: ',cnam,' =',jval
      write(ifle,'(1x,a,I8)') 'MSG :it must be >=',ival
      write(ifle,'(1x,4a)') 'MSG: reset ',cnam,' in ',mdnm
      write(ifle,'(1x,a)') '(nml_chksiz)'
      ierr=1
      end subroutine nml_chksiz
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_errmsg0(ifle,ios,nlst,ierr)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
! --- [dummy arguments]
!
      integer     ,intent(in)  :: ifle,ios
      character(*),intent(in)  :: nlst
      integer     ,intent(out) :: ierr
      ierr=0
      if( ios.eq.0 ) return
      if( ios.gt.0 ) then
        write(ifle,'(1x,a)') 'ERR: read error'
      else
        write(ifle,'(1x,a)') 'MSG : end of file'
      endif
      write(ifle,*) 'MSG: namelist = ',nlst
      write(ifle,*) 'MSG : iostat =',ios
      write(ifle,*) 'MSG :(nml_errmsg0)'
      ierr=1
      end subroutine nml_errmsg0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_init
     &  (MAT_NO,t0,y0,p0,u0,v0,w0,rho0,vsgm0,aks0,potn0,
     &   t02,y02,u02,v02,w02,rho02,
     &   fai_iMHD,A_iMHD,eps_iMHD,eps0_iMHD,mu_iMHD,
     &   mu0_iMHD,sigma_iMHD,sigma0_iMHD,Current_iMHD,OMEGA_iMHD,
     &   ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- MXMATERIAL
      use module_dimension
      use module_constant
      use module_hpcutil,only  : my_rank
      use module_io,only       : ifli,ifll,ifle,cntlnam
      use module_model,only    : idrdp,mach0,comp,incomp,icaltb,
     &                           sles,dles,lles,ke,RSM,ical_MHD,SDES
      use module_material,only : nofld,nflud,nosld,nsold,rsld,BGCOMP
      use module_scalar,only   : rns_scl,ides,ical_cavi,ical_FC,pot_scl,
     &                           ip_pefc
      use module_Euler2ph,only : ieul2ph
      use module_VOF,     only : ical_vof
      use module_species,only  : spcnam,act
      use module_initial,only  : Ke_FC,Kc_FC,Pv_FC,Kap_FC,
     &                           sigma_FC,contang_FC,
     &                           surf_tenFC,
     &                           elecCondFC,ionCondFC,potn_c
!
      implicit none
!
! 1. Input uniform initial condition
!
! --- [dummy arguments]
!
      integer,intent(in)  :: MAT_NO(0:MXMAT)
      real*8 ,intent(out) :: t0(MXMAT),y0(MXMAT,mxcomp)
      real*8 ,intent(out) :: p0(MXMAT),u0(MXMAT),v0(MXMAT),w0(MXMAT)
      real*8 ,intent(out) :: vsgm0(MXMAT,2),rho0(MXMAT)
      real*8 ,intent(out) :: aks0(MXMAT,mxrans),potn0(MXMAT,mxpotn)
      real*8 ,intent(out) :: t02(MXMAT),y02(MXMAT,mxcomp)
      real*8 ,intent(out) :: u02(MXMAT),v02(MXMAT),w02(MXMAT)
      real*8 ,intent(out) :: rho02(MXMAT)
      real*8 ,intent(out) :: fai_iMHD(MXMAT,2)
      real*8 ,intent(out) :: A_iMHD(MXMAT,3,2)
      real*8 ,intent(out) :: eps_iMHD(MXMAT)
      real*8 ,intent(out) :: eps0_iMHD(MXMAT)
      real*8 ,intent(out) :: mu_iMHD(MXMAT)
      real*8 ,intent(out) :: mu0_iMHD(MXMAT)
      real*8 ,intent(out) :: sigma_iMHD(MXMAT)
      real*8 ,intent(out) :: sigma0_iMHD(MXMAT)
      real*8 ,intent(out) :: Current_iMHD(MXMAT,3,2)
      real*8 ,intent(out) :: OMEGA_iMHD(MXMAT)
      integer,intent(out) :: ierror
!
! --- [namelist]
!
      integer,parameter :: mrans=50,mcomp=200,mpotn=10
      integer:: IMAT_U,iflag,IMD,IMAT,j,I
      real*8 :: aks(mrans),potential(mpotn),potn_cond(mpotn)
      real*8 :: t,ys(mcomp)
      real*8 :: p,u,v,w,vsgm(2),dens
      real*8 :: t2,ys2(mcomp),u2,v2,w2,dens2,dum1,dum2
      real*8 :: MHD_init_FAI(2),MHD_init_A(3,2),MHD_eps,MHD_eps0,
     &          MHD_mu,MHD_mu0,MHD_sigma,MHD_sigma0,MHD_current(3,2),
     &          MHD_f
      real*8 :: KcFC,KeFC,PvFC,KappaFC,sigma,contact_angle
      integer:: surface_tension
      real*8 :: elec_cond,ion_cond
      
!      integer:: water_vapor
      namelist /initial/ IMAT_U,t,ys,p,u,v,w,vsgm,dens,aks,
     &                   potential,potn_cond
      namelist /E2P_init/ IMAT_U,t2,ys2,u2,v2,w2,dens2
      namelist /VOF_init/ IMAT_U,dens2
      namelist /CAVI_init/ IMAT_U,dens2
      namelist /FC_parameter/ IMAT_U,dens2,KcFC,KeFC,PvFC,
     &                   KappaFC,sigma,contact_angle,
     &                   Surface_tension,
     &                   elec_cond,ion_cond
      namelist /MHD_init/ IMAT_U,
     &                    MHD_init_FAI,
     &                    MHD_init_A,
     &                    MHD_eps,
     &                    MHD_eps0,
     &                    MHD_mu,
     &                    MHD_mu0,
     &                    MHD_sigma,
     &                    MHD_sigma0,
     &                    MHD_current,
     &                    MHD_f
      
!      namelist /MHD_init/ IMAT_U,MHD_init_FAI,MHD_init_A,
!     &                    MHD_eps,MHD_mu,MHD_sigma,
!     &                    MHD_current,MHD_f
!
! --- [local entities]
!
      real*8,parameter  :: undef=-huge(1.d0)
      integer,parameter :: iundef=-huge(1)
      integer :: IMAT_UU,IIMAT,ICOM,IBGG
      real*8  :: sum
      integer :: ios=0,ierr1=0,nfld,nsld,nfsld,noinit(MXMAT)
!
! --- 
!
      ierror=0
      nfld=0
      nsld=0
      nfsld=0
      noinit(:)=0
!
      ifli=10               !my_rank+10
      if(ifli<=0) ifli=1
      open(unit=ifli,file=cntlnam,status='unknown',
     &     action="read",iostat=ierror)
      if(ierror>0) then
        write(ifle,*) '*** can not open ',cntlnam,' file in nml_init',
     &            'ifli  =', ifli,
     &            'iostat=', ierror   
        goto 9999
      end if
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
      potential=undef
      potn_cond=undef
!
      read(ifli,initial,iostat=ios)
      if(ios.lt.0) goto 351
      call nml_errmsg0(ifle,ios,'initial',ierr1)
      if( ierr1.ne.0 ) then
        write(ifle,*) '### ERR-1: reading error'
        write(ifle,*) '[&initial] in fflow.ctl'
        call FFRABORT(1,'namelist/initial')
      endif
      if(IMAT_U.eq.iundef) then
        write(ifle,*) '### ERR-2 : data error'
        write(ifle,*) 'lack of data : IMAT_U'
        write(ifle,*) '[&initial] in fflow.ctl'
        call FFRABORT(1,'namelist/initial')
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
        write(ifle,*) '### ERR-5: IMAT_U=0 is NOT supported'
        write(ifle,*) '### IMAT_U > 0 or IMAT_U < 0'
        write(ifle,*) '[&initial] in fflow.ctl'
        call FFRABORT(1,'namelist/initial')
      endif
      nfsld=nfsld+1
!      if(nfsld.gt.nflud+nsold) then
!        write(ifle,*) '### ERR-7: too many set of namelist &initial'
!        call FFRABORT(1,'namelist/initial')
!      endif
      noinit(nfsld)=IMAT_U
      do i=1,nfsld-1
      if(IMAT_U.eq.noinit(i))then
        write(ifle,*) '### ERR-8 : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'this IMAT_U have been defined'
        write(ifle,*) '[&initial] in fflow.ctl'
        call FFRABORT(1,'namelist/initial')
      endif
      enddo
!
      if( t.eq.undef ) then
        write(ifle,*) '### ERR-9 : data error'
        write(ifle,*) 'lack of data : t'
        write(ifle,*) 'IMAT_U= ',IMAT_U
        call FFRABORT(1,'namelist/initial')
      elseif( t.le.0.d0 ) then
        write(ifle,*) '### ERR-10 : data error'
        write(ifle,*) 't = ',t
        write(ifle,*) 'it must be > 0' 
        write(ifle,*) 'IMAT_U= ',IMAT_U
        call FFRABORT(1,'namelist/initial')
      endif
!
      if(IMAT_U.lt.0) then
        write(ifll,*) ' ### ATT: Solid material'
        write(ifll,*) ' [&initial: ys,p,u,v,w,vsgm,dens,aks]',
     &                        'are NEGLECTED'
        goto 360
      endif
!
      if(dens.eq.undef) then
        write(ifle,*) '### ERR-16 : data error'
        write(ifle,*) 'lack of data : dens'
        call FFRABORT(1,'namelist/initial')
      elseif( dens.le.0.d0 ) then
        write(ifle,*) '### ERR-17 : data error'
        write(ifle,*) 'dens = ',dens
        write(ifle,*) 'it must be > 0'
        call FFRABORT(1,'namelist/initial')
      endif
!
      dum1=-1.d2
      IBGG=1
      sum=0.d0
      do 100 ICOM=1,ncomp
      sum=sum+ys(ICOM)
      if(ys(ICOM)>dum1) then
        dum1=ys(ICOM)
        IBGG=ICOM
      endif
      if(ys(ICOM).eq.undef) then
        write(ifle,*) '### ERR-11 : data error'
        write(ifle,*) 'lack of data : ys(ICOM)'
        write(ifle,*) 'ICOM =',ICOM
        call FFRABORT(1,'namelist/initial')
      elseif( ys(ICOM).lt.0.d0 ) then
        write(ifle,*) '### ERR-12 : data error'
        write(ifle,*) 'ys(ICOM) = ',ys(ICOM)
        write(ifle,*) 'ICOM =',ICOM
        write(ifle,*) 'it must be >= 0'
        call FFRABORT(1,'namelist/initial')
      endif
  100 continue
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      
      if(IMAT.gt.0) then
        IMAT_UU=nofld(IMAT)
      elseif(IMAT.lt.0) then
        IMAT_UU=nosld(-IMAT)
      endif
      if(IMAT<=0) cycle
      if(IMAT_UU/=IMAT_U) cycle
      
      if(act(BGCOMP(IMAT))==0) then
         write(ifle,'(1X,a,2I4)') 'ERR: material no= ',IMAT_UU,IMAT
         write(ifle,'(1X,a,I4)') 
     &  'ERR: bg_gas must be acvtive, species no=',BGCOMP(IMAT)
        call FFRABORT(1,'ERR: bg_gas ,call you FFR supportor')
      endif
      if(ys(BGCOMP(IMAT))<SML) then
        write(ifle,'(1X,a,I4)') 'ERR: material no= ',IMAT_UU
        write(ifle,'(3a)') 
     &  'ERR: Initial mass fraction of Back-Ground Gas :[',
     &   trim(spcnam(BGCOMP(IMAT))),
     &  '] MUST be > 0'
        write(ifle,'(2a,I4,3a)') 
     &    'MSG: Choise other component as Back Ground Gas',
     &    ', For example: ICOM=',IBGG,'  [',trim(spcnam(IBGG)),']'
        call FFRABORT(1,'ERR:')
      endif
      enddo

      if( sum.le.0.d0 ) then
        write(ifle,*) '### ERR-13 : data error'
        write(ifle,*) 'ys = ',(ys(i),i=1,ncomp)
        write(ifle,*) 'sum of them must be > 0'
        call FFRABORT(1,'namelist/initial')
      endif
!
      if( p.eq.undef ) then
        write(ifle,*) '### ERR-14 : data error'
        write(ifle,*) 'lack of data : p'
        call FFRABORT(1,'namelist/initial')
      elseif( p.le.0.d0 ) then
        write(ifle,*) '### ERR-15 : data error'
        write(ifle,*) 'p = ',p
        write(ifle,*) 'it must be > 0'
        call FFRABORT(1,'namelist/initial')
      endif
!
      if( u.eq.undef ) then
        write(ifle,*) '### ERR-18 : data error'
        write(ifle,*) 'lack of data : u'
        call FFRABORT(1,'namelist/initial')
      endif
      if( v.eq.undef ) then
        write(ifle,*) '### ERR-19 : data error'
        write(ifle,*) 'lack of data : v'
        call FFRABORT(1,'namelist/initial')
      endif
      if( w.eq.undef ) then
        write(ifle,*) '### ERR-20 : data error'
        write(ifle,*) 'lack of data : w'
        call FFRABORT(1,'namelist/initial')
      endif
!
      if( vsgm(1).eq.undef ) then
        vsgm=0.d0
      else
        if( vsgm(1).le.0.d0 ) then
          write(ifle,*) '### ERR-21 : data error'
          write(ifle,*) 'vsgm(1) = ',vsgm(1)
          write(ifle,*) 'it must be > 0'
          call FFRABORT(1,'namelist/initial')
        endif
        if( vsgm(2).eq.undef ) then
          write(ifle,*) '### ERR-22 : data error'
          write(ifle,*) 'lack of data : vsgm(2)'
          call FFRABORT(1,'namelist/initial')
        elseif( vsgm(2).le.0.d0 ) then
          write(ifle,*) '### ERR-23 : data error'
          write(ifle,*) 'vsgm(2) = ',vsgm(2)
          write(ifle,*) 'it must be > 0'
          call FFRABORT(1,'namelist/initial')
        endif
      endif
      if(rns_scl) then
        do IMD=1,nrans
        if( aks(IMD).eq.undef ) then
          write(ifle,*) '### ERR-24 : data error'
          write(ifle,*) 'lack of data : aks(IMD)'
          write(ifle,*) 'IMD =',IMD
          call FFRABORT(1,'namelist/initial')
        elseif( aks(IMD).lt.0. ) then
          write(ifle,*) '### ERR-25 : data error'
          write(ifle,*) 'aks(IMD) = ',aks(IMD)
          write(ifle,*) 'it must be >= 0'
          write(ifle,*) 'IMD =',IMD
          call FFRABORT(1,'namelist/initial')
        elseif(abs(aks(IMD)).lt.SML.and.
     &         icaltb==SDES.and.IMD==ides) then
          write(ifle,*) '### ERR-25 : Lack of Nu data error'
          write(ifle,*) 'aks(IMD) = ',aks(IMD)
          write(ifle,*) 'it must be > 0; or, set value Mu/rho'
          write(ifle,*) 'IMD =',IMD
          call FFRABORT(1,'namelist/initial')
        endif
        enddo
      else
        aks=0.d0
      endif
!pot_scl
      if(pot_scl) then
        do IMD=1,NPOTN
        if(potential(IMD).eq.undef) then
          write(ifle,*) '### ERR-26 : data error'
          write(ifle,*) 'lack of data : potential'
          write(ifle,*) 'IPOT =',IMD
          call FFRABORT(1,'namelist/initial')
        endif
        if(potn_cond(IMD).eq.undef) then
          write(ifle,*) '### ERR-26 : data error'
          write(ifle,*) 'lack of data : potn_cond'
          write(ifle,*) 'IPOT =',IMD
          call FFRABORT(1,'namelist/initial')
        endif
        enddo
      endif
!
 360  continue
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT.gt.0) then
        IMAT_UU=nofld(IMAT)
      elseif(IMAT.lt.0) then
        IMAT_UU=nosld(-IMAT)
      endif
      if(IMAT_UU.eq.IMAT_U.and.IMAT.gt.0) then
        t0(IIMAT)=t
        p0(IIMAT)=p
        u0(IIMAT)=u
        v0(IIMAT)=v
        w0(IIMAT)=w
        rho0(IIMAT)=dens
        do 200 ICOM=1,ncomp
        y0(IIMAT,ICOM)=ys(ICOM)
  200   continue
        vsgm0(IIMAT,1)=vsgm(1)
        vsgm0(IIMAT,2)=vsgm(2)
        if(rns_scl) then
          aks0(IIMAT,1:nrans)=aks(1:nrans)
        endif
        if(pot_scl) then
          potn0(IIMAT,1:npotn)=potential(1:npotn)
          potn_c(IIMAT,1:npotn)=potn_cond(1:npotn)
        endif
      elseif(IMAT_UU.eq.IMAT_U.and.IMAT.lt.0) then
        t0(IIMAT)=t
        p0(IIMAT)=0.d0
        u0(IIMAT)=0.d0
        v0(IIMAT)=0.d0
        w0(IIMAT)=0.d0
        rho0(IIMAT)=rsld(-IMAT)
        do 210 ICOM=1,ncomp
        y0(IIMAT,ICOM)=0.d0
  210   continue
        vsgm0(IIMAT,1)=0.d0
        vsgm0(IIMAT,2)=0.d0
        if(rns_scl) then
          aks0(IIMAT,1:nrans)=0.d0
        endif
        if(pot_scl) then
          potn0(IIMAT,1:npotn)=potential(1:npotn)
          potn_c(IIMAT,1:npotn)=potn_cond(1:npotn)
        endif
      endif
      enddo
      goto 350
 351  continue
!---------------------
! --- Euler two phase 
!---------------------
      if(ieul2ph>0) then
        rewind ifli
 450    continue
        IMAT_U=iundef
        t2   =undef
        ys2  =undef
        u2   =undef
        v2   =undef
        w2   =undef
        dens2=undef
        read(ifli,E2P_init,iostat=ios)
        if(ios.lt.0) goto 451
        call nml_errmsg0(ifle,ios,'E2P_init',ierr1)
        if(ierr1.ne.0) then
          write(ifle,*) '### ERR-1: reading error'
          call FFRABORT(1,'namelist/E2P_init')
          if(IMAT_U.eq.iundef) then
            write(ifle,*) '### ERR-2 : data error'
            write(ifle,*) 'lack of data : IMAT_U'
            call FFRABORT(1,'namelist/E2P_init')
          endif
        endif
        if(IMAT_U.eq.0) then
          write(ifle,*) '### ERR-5: IMAT_U=0 is NOT supported'
          write(ifle,*) '### IMAT_U > 0 or IMAT_U < 0'
          call FFRABORT(1,'namelist/E2P_init')
        endif
        if( t2.eq.undef ) then
          write(ifle,*) '### ERR-9 : data error'
          write(ifle,*) 'lack of data : t'
          call FFRABORT(1,'namelist/E2P_init')
        elseif( t2.le.0.d0 ) then
          write(ifle,*) '### ERR-10 : data error'
          write(ifle,*) 't = ',t2
          write(ifle,*) 'it must be > 0' 
          call FFRABORT(1,'namelist/E2P_init')
        endif
        if( dens2.eq.undef ) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'lack of data : dens'
          call FFRABORT(1,'namelist/E2P_init')
        elseif( dens2.le.0.d0 ) then
          write(ifle,*) '### ERR-17 : data error'
          write(ifle,*) 'dens2 = ',dens2
          write(ifle,*) 'it must be > 0'
          call FFRABORT(1,'namelist/E2P_init')
        endif
        sum=0.d0
        do 110 ICOM=1,ncomp
        sum=sum+ys2(ICOM)
        if( ys2(ICOM).eq.undef ) then
          write(ifle,*) '### ERR-11 : data error'
          write(ifle,*) 'lack of data : ys2(ICOM)'
          write(ifle,*) 'ICOM =',ICOM
          call FFRABORT(1,'namelist/E2P_init')
        elseif( ys2(ICOM).lt.0.d0 ) then
          write(ifle,*) '### ERR-12 : data error'
          write(ifle,*) 'ys2(ICOM) = ',ys2(ICOM)
          write(ifle,*) 'ICOM =',ICOM
          write(ifle,*) 'it must be >= 0'
          call FFRABORT(1,'namelist/E2P_init')
        endif
  110   continue
        if( sum.le.0.d0 ) then
          write(ifle,*) '### ERR-13 : data error'
          write(ifle,*) 'ys2 = ',(ys2(i),i=1,ncomp)
          write(ifle,*) 'sum of them must be > 0'
          call FFRABORT(1,'namelist/E2P_init')
        endif
!
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          IMAT_UU=nofld(IMAT)
        elseif(IMAT.lt.0) then
          IMAT_UU=nosld(-IMAT)
        endif
        if(IMAT_UU.eq.IMAT_U.and.IMAT.gt.0) then
          t02(IIMAT)=t2
          u02(IIMAT)=u2
          v02(IIMAT)=v2
          w02(IIMAT)=w2
          rho02(IIMAT)=dens2
          do ICOM=1,ncomp
          y02(IIMAT,ICOM)=ys2(ICOM)
          enddo
        elseif(IMAT_UU.eq.IMAT_U.and.IMAT.lt.0) then
          t02(IIMAT)=t2
          u02(IIMAT)=0.d0
          v02(IIMAT)=0.d0
          w02(IIMAT)=0.d0
          rho02(IIMAT)=dens2
          do ICOM=1,ncomp
          y02(IIMAT,ICOM)=0.d0
          enddo
        endif
        enddo
        goto 450
 451    continue
      endif  ! END Euler two phase
!
!---------------------
! --- VOF two phase 
!---------------------
!
      if(ical_vof==1) then
        rewind ifli
 550    continue
        IMAT_U=iundef
        dens2=undef
!
        ios=0

        read(ifli,VOF_init,iostat=ios)
        call nml_errmsg0(ifle,ios,'vof_init',ierr1)
        if(ios.lt.0) goto 551
        if( ierr1.ne.0 ) then
          write(ifle,*) '### ERR-1: reading error'
          call FFRABORT(1,'namelist/Vof_init')
          if(IMAT_U.eq.iundef) then
            write(ifle,*) '### ERR-2 : data error'
            write(ifle,*) 'lack of data : IMAT_U'
            call FFRABORT(1,'namelist/vof_init')
          endif
        endif
        if(IMAT_U.eq.0) then
          write(ifle,*) '### ERR-5: IMAT_U=0 is NOT supported'
          write(ifle,*) '### IMAT_U > 0 or IMAT_U < 0'
          call FFRABORT(1,'namelist/vof_init')
        endif
        if(dens2.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'lack of data : dens'
          call FFRABORT(1,'namelist/vof_init')
        elseif(dens2.le.0.d0) then
          write(ifle,*) '### ERR-17 : data error'
          write(ifle,*) 'dens2 = ',dens2
          write(ifle,*) 'it must be > 0'
          call FFRABORT(1,'namelist/vof_init')
        endif
!
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          IMAT_UU=nofld(IMAT)
        elseif(IMAT.lt.0) then
          IMAT_UU=nosld(-IMAT)
        endif
        if(IMAT_UU.eq.IMAT_U.and.IMAT.gt.0) then
          rho02(IIMAT)=dens2
        elseif(IMAT_UU.eq.IMAT_U.and.IMAT.lt.0) then
          rho02(IIMAT)=dens2
        endif
        enddo
!
        goto 550
 551    continue
      endif  ! END VOF two phase
!
      if(ical_cavi==1) then
        rewind ifli
 650    continue
        IMAT_U=iundef
        dens2=undef
        ios=0

        read(ifli,CAVI_init,iostat=ios)
        call nml_errmsg0(ifle,ios,'CAVI_init',ierr1)
        if(ios.lt.0) goto 651
        if( ierr1.ne.0 ) then
          write(ifle,*) '### ERR-1: reading error'
          call FFRABORT(1,'namelist/vof_init/CAVI_init/FC_parameter')
          if(IMAT_U.eq.iundef) then
            write(ifle,*) '### ERR-2 : data error'
            write(ifle,*) 'lack of data : IMAT_U'
            call FFRABORT(1,'namelist/vof_init/CAVI_init/FC_parameter')
          endif
        endif
        if(IMAT_U.eq.0) then
          write(ifle,*) '### ERR-5: IMAT_U=0 is NOT supported'
          write(ifle,*) '### IMAT_U > 0 or IMAT_U < 0'
          call FFRABORT(1,'namelist/vof_init/CAVI_init')
        endif
        if(dens2.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'lack of data : dens'
          call FFRABORT(1,'namelist/vof_init/CAVI_init')
        elseif(dens2.le.0.d0) then
          write(ifle,*) '### ERR-17 : data error'
          write(ifle,*) 'dens2 = ',dens2
          write(ifle,*) 'it must be > 0'
          call FFRABORT(1,'namelist/vof_init/CAVI_init')
        endif
!
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          IMAT_UU=nofld(IMAT)
        elseif(IMAT.lt.0) then
          IMAT_UU=nosld(-IMAT)
        endif
        if(IMAT_UU.eq.IMAT_U.and.IMAT.gt.0) then
          rho02(IIMAT)=dens2
        elseif(IMAT_UU.eq.IMAT_U.and.IMAT.lt.0) then
          rho02(IIMAT)=dens2
        endif
        enddo
!
        goto 650
 651    continue
      endif  ! END VOF two phase
!
      if(ical_FC>0) then
        rewind ifli
 750    continue
        IMAT_U=iundef
        dens2=undef
!        water_vapor=1
        KeFC=0.d0
        KcFC=0.d0
        PvFC=1.01325d5
        KappaFC=0.d0
        sigma=0.d0
        contact_angle=0.d0
!        electro_osmotic=0
        surface_tension=0
        elec_cond=0.d0
        ion_cond=0.d0
!
        ios=0
!
        read(ifli,FC_parameter,iostat=ios)
        if(ios.lt.0) goto 751
        write(ifll,*) 'MSG: &FC_parameter IMAT_U= ',IMAT_U
        call nml_errmsg0(ifle,ios,'FC_parameter',ierr1)
        if(ierr1.ne.0) then
          write(ifle,*) '### ERR-1: reading error'
          call FFRABORT(1,'namelist/vof_init/CAVI_init/FC_parameter')
          if(IMAT_U.eq.iundef) then
            write(ifle,*) '### ERR-2 : data error'
            write(ifle,*) 'lack of data : IMAT_U'
            call FFRABORT(1,'namelist/vof_init/CAVI_init/FC_parameter')
          endif
        endif
        if(IMAT_U.eq.0) then
          write(ifle,*) '### ERR-5: IMAT_U=0 is NOT supported'
          write(ifle,*) '### IMAT_U > 0 or IMAT_U < 0'
          call FFRABORT(1,'namelist/vof_init/CAVI_init')
        endif
        if(dens2.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'lack of data : dens'
          call FFRABORT(1,'namelist/vof_init/CAVI_init')
        elseif(dens2.le.0.d0) then
          write(ifle,*) '### ERR-17 : data error'
          write(ifle,*) 'dens2 = ',dens2
          write(ifle,*) 'it must be > 0'
          call FFRABORT(1,'namelist/vof_init/CAVI_init')
        endif
!
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          IMAT_UU=nofld(IMAT)
        elseif(IMAT.lt.0) then
          IMAT_UU=nosld(-IMAT)
        endif
        if(IMAT_UU.eq.IMAT_U.and.IMAT.gt.0) then
          rho02(IIMAT)=dens2
          if(ical_FC>0) then 
!            if(water_vapor>ncomp.or.water_vapor<1) then
!              call FFRABORT(1,'ERR: [water_vapor<ncomp]')
!            else
!              vap_no(IIMAT)=water_vapor 
              Ke_FC(IIMAT)=KeFC
              Kc_FC(IIMAT)=KcFC
              Pv_FC(IIMAT)=PvFC
              Kap_FC(IIMAT)=KappaFC
              sigma_FC(IIMAT)=sigma
              contang_FC(IIMAT)=contact_angle
!              osmoticFC(IIMAT)=electro_osmotic
              surf_tenFC(IIMAT)=surface_tension
              elecCondFC(IIMAT)=elec_cond
              ionCondFC(IIMAT)=ion_cond
              potn_c(IIMAT,ip_pefc(1))=elec_cond
              potn_c(IIMAT,ip_pefc(2))=ion_cond
!            endif
          endif 
        elseif(IMAT_UU.eq.IMAT_U.and.IMAT.lt.0) then
          rho02(IIMAT)=dens2
          if(ical_FC>0) then
            elecCondFC(IIMAT)=elec_cond
            ionCondFC(IIMAT)=ion_cond
            potn_c(IIMAT,ip_pefc(1))=elec_cond
            potn_c(IIMAT,ip_pefc(2))=ion_cond
          endif
        endif
        enddo
!
        goto 750
 751    continue
      endif  ! END VOF two phase

!
      if(ical_MHD==1.or.ical_MHD==2) then
        fai_iMHD(:,:)=0.d0
        A_iMHD(:,:,:)=0.d0
        eps_iMHD(:)=0.d0
        eps0_iMHD(:)=8.8541187818D-12
        mu_iMHD(:)=0.d0
        mu0_iMHD(:)=12.56637061D-7
        sigma_iMHD(:)=0.d0
        sigma0_iMHD(:)=0.d0
        Current_iMHD(:,:,:)=0.d0
!
        rewind ifli
 850    continue
        IMAT_U=iundef
        MHD_init_FAI=undef
        MHD_init_A=undef
        MHD_eps=undef
        MHD_eps0=8.8541187818D-12
        MHD_mu=undef
        MHD_mu0=12.56637061D-7
        MHD_sigma=undef
        MHD_sigma0=undef
        MHD_current=undef
        MHD_f=undef
!
        ios=0
        read(ifli,MHD_init,iostat=ios)
        if(ios.lt.0) goto 851
        call nml_errmsg0(ifle,ios,'MHD_init',ierr1)
        if( ierr1.ne.0 ) then
          write(ifle,*) '### ERR-1: reading error'
          write(ifle,*) '### ERR-2: Material No = ',IMAT_U
          call FFRABORT(1,'namelist/MHD_init')
          if(IMAT_U.eq.iundef) then
            write(ifle,*) '### ERR-2 : data error'
            write(ifle,*) 'lack of data : IMAT_U'
            call FFRABORT(1,'namelist/MHD_init')
          endif
        endif
!
        if(IMAT_U.eq.0) then
          write(ifle,*) '### ERR-5: IMAT_U=0 is NOT supported'
          write(ifle,*) '### IMAT_U > 0 or IMAT_U < 0'
          call FFRABORT(1,'namelist/MHD_init')
        endif
!
        if(MHD_init_FAI(1).eq.undef.or.MHD_init_FAI(2).eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_init_FAI'
          call FFRABORT(1,'namelist/MHD_init')
        endif
        DO I=1,3
        if(MHD_init_A(I,1).eq.undef.or.MHD_init_A(I,2).eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_init_FAI'
          call FFRABORT(1,'namelist/MHD_init')
        endif
        enddo
        DO I=1,3
        if(MHD_current(I,1).eq.undef.or.MHD_current(I,2).eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : initial: MHD_current'
          call FFRABORT(1,'namelist/MHD_init')
        endif
        enddo
        if(MHD_eps.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_eps'
          call FFRABORT(1,'namelist/MHD_init')
        endif
        if(MHD_eps0.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_eps0'
          call FFRABORT(1,'namelist/MHD_init')
        endif
        if(MHD_mu.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_mu'
          call FFRABORT(1,'namelist/MHD_init')
        endif
        if(MHD_mu0.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_mu0'
          call FFRABORT(1,'namelist/MHD_init')
        endif
        if(MHD_sigma.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_sigma'
          call FFRABORT(1,'namelist/MHD_init')
        endif
!        if(MHD_sigma0.eq.undef) then
!          write(ifle,*) '### ERR-16 : data error'
!          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
!          write(ifle,*) 'lack of data : MHD_sigma0'
!          call FFRABORT(1,'namelist/MHD_init')
!        endif
        if(MHD_f.eq.undef) then
          write(ifle,*) '### ERR-16 : data error'
          write(ifle,*) 'MSG: IMAT_U= ',IMAT_U
          write(ifle,*) 'lack of data : MHD_f'
          call FFRABORT(1,'namelist/MHD_init')
        endif
!
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          IMAT_UU=nofld(IMAT)
        elseif(IMAT.lt.0) then
          IMAT_UU=nosld(-IMAT)
        endif
        if(IMAT_UU.eq.IMAT_U.and.IMAT.gt.0) then
          fai_iMHD(IIMAT,1)=MHD_init_FAI(1)
          fai_iMHD(IIMAT,2)=MHD_init_FAI(2)
          DO I=1,3
          A_iMHD(IIMAT,I,1)=MHD_init_A(I,1)
          A_iMHD(IIMAT,I,2)=MHD_init_A(I,2)
          Current_iMHD(IIMAT,I,1)=MHD_current(I,1)
          Current_iMHD(IIMAT,I,2)=MHD_current(I,2)
          enddo
          eps_iMHD(IIMAT)=MHD_eps
          eps0_iMHD(IIMAT)=MHD_eps0
          mu_iMHD(IIMAT)=MHD_mu
          mu0_iMHD(IIMAT)=MHD_mu0
          sigma_iMHD(IIMAT)=MHD_sigma
          sigma0_iMHD(IIMAT)=MHD_sigma0
          OMEGA_iMHD(IIMAT)=MHD_f
        elseif(IMAT_UU.eq.IMAT_U.and.IMAT.lt.0) then
          fai_iMHD(IIMAT,1)=MHD_init_FAI(1)
          fai_iMHD(IIMAT,2)=MHD_init_FAI(2)
          DO I=1,3
          A_iMHD(IIMAT,I,1)=MHD_init_A(I,1)
          A_iMHD(IIMAT,I,2)=MHD_init_A(I,2)
          Current_iMHD(IIMAT,I,1)=MHD_current(I,1)
          Current_iMHD(IIMAT,I,2)=MHD_current(I,2)
          enddo
          eps_iMHD(IIMAT)=MHD_eps
          eps0_iMHD(IIMAT)=MHD_eps0
          mu_iMHD(IIMAT)=MHD_mu
          mu0_iMHD(IIMAT)=MHD_mu0
          sigma_iMHD(IIMAT)=MHD_sigma
          sigma0_iMHD(IIMAT)=MHD_sigma0
          OMEGA_iMHD(IIMAT)=MHD_f
        endif
        enddo
        goto 850
 851    continue
!
        dum1=OMEGA_iMHD(1)
        do IIMAT=2,NMAT
          dum2=OMEGA_iMHD(IIMAT)
          if(abs(dum1-dum2)>SML) then
            call FFRABORT(1,'All Material Hz MUST be same: MHD_f')
          else
            dum1=OMEGA_iMHD(IIMAT)
          endif
        enddo
      endif
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
          write(ifle,*) '### ERR-3: IMAT_U=',nofld(j),
     &                  ' Initial condition is NOT specified in',
     &                  ' namelist/initial'
          call FFRABORT(1,'namelist/initial')
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
          write(ifle,*) '### ERR-4: IMAT_U=',nosld(j),
     &                  ' Initial condition is NOT specied in',
     &                  ' namelist/solid'
          call FFRABORT(1,'namelist/initial')
      endif
      enddo
!
      if(nfld.lt.1.and.nsld.lt.1) then
        write(ifle,*) '### ERR-26: namelist/initial must be specified'
        call FFRABORT(1,'namelist/initial')
      elseif(nfld.lt.1) then
        write(ifle,*) '### WRN: NO fluid field exists'
        write(ifle,*) '### WRN-27: This case is NOT supported yet'
!MHD_err        call FFRABORT(2,'namelist/initial')
      elseif(nfld+nsld.ne.nflud+nsold) then
        write(ifle,*) 
     &            '### ERR-28: no. of namelist-initial NOT equal to',
     &                'namelist-fluid + namelist-solid'
        call FFRABORT(3,'namelist/initial')
      endif
!
 300  continue
!
      close(ifli)
!
      return
!
 9999 continue
      write(ifle,*) '(nml_init)'
      ierror=1
      end subroutine nml_init
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine nml_listno(nl,lst,wrd,nn)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. Return list no.
!      nn : =-1; ambiguous word
!         : = 0; unknown word
!         : > 0; identified word
!
      implicit none
! --- [dummy arguments]
!
      integer     ,intent(in)  :: nl
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
      if(wrd(:lng).eq.' ') return
      ic=0
      do 100 n=1,nl
      if(lst(n)(:lng).eq.' ') goto 100
      do 101 i=1,lng
      if(wrd(i:i).eq.' ' ) goto 102
      if(wrd(i:i).ne.lst(n)(i:i)) goto 100
  101 continue
  102 continue
      ic=ic+1
      nn=n
  100 continue
      if( ic.gt.1 ) nn=-1
!
      end subroutine nml_listno
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      logical function nml_comp_eq(str,v1,v2,s,f,e)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- compare values (real*8), equal (true) or not (false)
!---------------------------------------------------------------------
      implicit none
      character(*),intent(in) :: str    ! string of value
      real*8,intent(in) :: v1,v2        ! values to compare
      integer,intent(in) :: s           ! sequential number
      integer,intent(in) :: f           ! file number to write error script
      integer,intent(inout) :: e        ! total number of error
!
      nml_comp_eq = .false.             ! not equal
      if( v1==v2 ) then
        write(f,*) "### error : data error"
        write(f,*) "lack of data '",str,"' in no.",s
        e = 1+e
        nml_comp_eq = .true.            ! equal
      endif
      end function nml_comp_eq
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      logical function nml_comp_gt(str,v1,v2,s,f,e)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- compare values (real*8), v1 > v2 (true) or not (false)
!---------------------------------------------------------------------
      implicit none
      character(*),intent(in) :: str    ! string of value
      real*8,intent(in) :: v1,v2        ! values to compare
      integer,intent(in) :: s           ! sequential number
      integer,intent(in) :: f           ! file number to write error script
      integer,intent(inout) :: e        ! total number of error
!
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
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      logical function nml_comp_ge(str,v1,v2,s,f,e)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! ---- compare values (real*8), v1 >= v2 (true) or not (false)
!---------------------------------------------------------------------
      implicit none
      character(*),intent(in) :: str    ! string of value
      real*8,intent(in) :: v1,v2        ! values to compare
      integer,intent(in) :: s           ! sequential number
      integer,intent(in) :: f           ! file number to write error script
      integer,intent(inout) :: e        ! total number of error
!
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
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine checkExpire(yr,mt,dy,hr,mi,sc,rtrnStat)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none

      integer,intent(in) :: yr,mt,dy,hr,mi,sc
      integer,intent(out) :: rtrnStat


      integer :: today(8),remainD,remainH


      call date_and_time(VALUES=today)
      remainD=(yr-today(1))*365+(mt-today(2))*30+(dy-today(3))
      remainH=(hr-today(5))*3600+(mi-today(6))*60+(sc-today(7))

      if(remainD<0) then
        rtrnStat=1
      else if(remainD>0) then
        rtrnStat=0
      else if(remainD==0) then
        if(remainH<=0) then
          rtrnStat=1
        else if(remainH>0) then
          rtrnStat=0
        end if
      end if


 1000 format(2x,'MSG: Licnese will expire in :',
     &            I4.4,'.',I2.2,'.',
     &            I2.2,'  ',I2.2,':',I2.2,':',I2.2)
      if(rtrnStat==0) then
        write(*,1000)yr,mt,dy,hr,mi,sc
      else if(rtrnStat==1) then
        write(*,*)
     &'Trial period has expired. Please contact your supportor'
      end if
!
       return
       end subroutine checkExpire
