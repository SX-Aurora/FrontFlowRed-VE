!
!      module module_boundary 
!      module module_cgsolver
!      module module_chemcntl
!      module module_debug
!      module module_deltat
!      module module_dimnsn
!      module module_flags
!      module module_material
!      module module_gravity
!      module module_io
!      module module_les
!      module module_model
!      module module_movegrid
!      module module_nowtime
!      module module_output
!      module module_param
!      module module_rans
!      module module_simple
!      module module_source
!      module module_species
!      module module_time
!      module module_turbparm
!      module module_chemreac
!      module module_usersub
!      module module_hpc
!      module module_flow_sound
!      module module_partitioner
!      module module_hpc_input
!      module module_anim
!      module module_Euler2ph
!      module module_vof
!      module module_radslv
!      module module_radgroup
!      module module_radwork
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_boundary
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(21),parameter,private :: modnam='(module_boundary)'
!
!
!< module for boundary conditions & domains > 
!
!-< 1. boundary conditions >-
!
!--< 1.1 flags >--
!.and.idis(nb)==0
      integer,parameter :: kxnone=-999999
!
! kxnone : value of kdbcnd which indicate undefined 
!
      integer,parameter :: kdprdc=-1     !=> dis
      integer,parameter :: kdsymm=-2
      integer,parameter :: kdintr=-3     !=> dis
      integer,parameter :: kdilet=-4
      integer,parameter :: kdolet=-5
      integer,parameter :: kdtchi=-6     !=> dis
      integer,parameter :: kdfire=-7
      integer,parameter :: kdpres=-8
      integer,parameter :: kdstag=-9
      integer,parameter :: kdbuff=-10
      integer,parameter :: kdsld=-11     !=> dis
      integer,parameter :: kdcvd=-12
      integer,parameter :: kdshutr=-13
      integer,parameter :: kdpors=-14
      integer,parameter :: kdovst=-15
!
! kdprdc : value of kdbcnd(0,:) for periodic boundary
! kdsymm : value of kdbcnd(0,:) for symmetric boundary
! kdintr : value of kdbcnd(0,:) for interface boundary
! kdilet : value of kdbcnd(0,:) for inlet boundary
! kdolet : value of kdbcnd(0,:) for outlet boundary
! kdtchi : value of kdbcnd(0,:) for touch boundary
! kdfire : value of kdbcnd(0,:) for fire wall boundary
! kdpres : value of kdbcnd(0,:) for pressure boundary
!                                  (absolute pressure)
! kdstag : value of kdbcnd(0,:) for stagnation pressure boundary
!                                  (absolute pressure)
!
      logical,save :: lsldbc=.false.,lovstbc=.false.
!
!
      integer,parameter :: kvnslp=1, kvfslp=2, kvlglw=3,
     &                     kvrogf=4, kvmodl=5,kvEWF=6
      integer,parameter :: ktdirc=1, ktneum=2, kttrns=3,ktEWF=4
 4    integer,parameter :: kydirc=1, kyneum=2, kytrns=3,kyEWF=4
      integer,parameter :: kkdirc=1, kkneum=2, kklglw=3
!
      integer,parameter :: mpotn=10
!
! kv* : value of kdbcnd(1,:) for velocity 
!     : kvnslp / no-slip                        (=1)
!     : kvfslp / free-slip                      (=2)
!     : kvlglw / log-law                        (=3)
!     : kvrogf / rough-wall                     (=4)
! kt* : value of kdbcnd(2,:) for temperature
!     : ktdirc / T=const. (Dirichlet)           (=1)
!     : ktneum / dT/dn=const. (Neumann)         (=2)
!     : kttrns / heat transfer                  (=3)
! ky* : value of kdbcnd(3,:) for mass fraction
!     : kydirc / Ys=const. (Dirichlet)          (=1)
!     : kyneum / dYs/dn=const. (Neumann)        (=2)
!     : kytrns / mass transfer                  (=3)
! kk* : value of kdbcnd(4,:) for RANS model
!     : kkdirc / k=const. (Dirichlet)           (=1)
!     : kkneum / dk/dn=const. (Neumann)         (=2)
!     : kklglw / log-law                        (=3)
!
!--< 1.2 values >--
!
      integer,private :: ncomp=0, nrans=0
      integer :: iscomp=0, israns=0,distrb=0
      integer :: KSUF=0
!
! ncomp  : no. of Ys(i)  (chemical species)
! nrans  : no. of AKs(i) (dependent variables of RANS model)
! iscomp : start index of wdbcnd for Ys(i)
! israns : start index of wdbcnd for AKs(i)
!
      integer         :: nbcnd =0,NBOUND=0,NAME_I=1
      integer,private :: nbvals=0
      integer,save,allocatable         :: nobcnd (:)
      integer,save,allocatable         :: boundIDMap(:,:)
      CHARACTER*80,allocatable         :: SFBOUN(:)
      INTEGER,allocatable              :: IFFACE(:,:,:)
      INTEGER,allocatable              :: NFBOUN(:)	
      integer,save,allocatable         :: kdbcnd (:,:)
      integer,save,allocatable,private :: iwbcnd (:)
      real*8 ,save,allocatable         :: wdbcnd(:,:)
      real*8 ,save,allocatable,private :: wdtbcnd(:,:)
      real*8 ,save,allocatable         :: adbcnd (:,:)
      integer,save,allocatable         :: LBC_INDEX(:)
      integer,save,allocatable         :: idis(:)
      integer,save,allocatable         :: thinw(:)
      integer,save,allocatable         :: nblade(:,:)
      real*8,save,allocatable          :: radprop(:,:)
      integer,save,allocatable         :: radwalltype(:,:)
      logical,save,allocatable         :: conntBC(:)
!
! nbcnd       : no. of boundary conditions
! nbvals      : 1st dimension size of "wdbcnd"
! nobcnd(:)   : domain no.
! kdbcnd(1,:) : for velocity
! kdbcnd(2,:) : for temperature
! kdbcnd(3,:) : for mass fraction
! kdbcnd(4,:) : for RANS model
! iwbcnd      : pointer for "wdtbcnd"
! adbcnd(1,:) : heat transfer coefficient
! adbcnd(2,:) : mass transfer coefficient
! adbcnd(3,:) : factor for standard deviation of inlet velocity
!             : adbcnd(3,:)*(velocity) is standard deviation
! adbcnd(4,:) : factor for maximul fluctuation of inlet velocity
!             : adbcnd(4,:)*(velocity) is maximul fluctuation
! wdbcnd(0,:) : cycle
! wdbcnd(1,:) : p [Pa]
! wdbcnd(2,:) : u [m/s]
! wdbcnd(3,:) : v [m/s]
! wdbcnd(4,:) : w [m/s]
! wdbcnd(5,:) : T [K] or (lambda)*dT/dn [W/m^2]
! wdbcnd      : Ys [-] or (rho)D*d(Ys)/dn [kg/s/m^2]
! (6:5+ncomp,:)
! wdbcnd      : AKs or D*d(AKs)/dn
! (6+ncomp:5+ncomp+nrans,:)
! wdtbcnd     : time table for "wdbcnd"
!
!------- Jiang appended
! radprop(2,:):	wall radiative property, 
!       1: emissivity coef. 
!       2: absorbility coef.
! radwalltype(2,:): wall radiative type, and group flag
!	1: radiative type
!	1, diffuse-wall	,OK 
!	2, mirror-wall		,OK
!	3, directional-diffuse wall	,OK
!	4, user defined
!	2: groupflag,	1 or 0
!
!-< 2. boundary domains >-
!
      integer,save,allocatable :: ivbcnd(:)
      integer,save,allocatable :: lvbcnd(:)
      integer,save,allocatable :: ivpair(:)
!
! ivbcnd : pointer for "lvbcnd"
! lvbcnd : vertex no. in boundary domain
!
!
!-< 3. side no. of cell ( only for inner boundary ) >-
!
      integer,save,allocatable :: kdbinr(:,:)
      integer,save,allocatable :: icbinr(:)
      integer,save,allocatable :: lcbinr(:)
!
! kdbinr(1,:) : boundary condition no.
! kdbinr(2,:) : side no. of the boundary condition
! icbinr      : pointer for "lcbinr"
! lcbinr      : cell no. against side-1 of the inner boundary
!
!
!-< 99. for internal check >-
!
      integer,private :: iset1=0, iset2=0
      character*80,allocatable :: boundName(:,:)
      integer,save             :: numbname,undef_bcno=0
      character*1,allocatable  :: prdcAxis(:)
      real(8),allocatable      :: prdcOffset(:,:)
!
      integer,save             :: dum_bc=0
      integer,save,allocatable :: bcmat(:)
!///////////////////////////////////////////////////////////////////////
      contains
!
!
!< #1. check consistency of values in argument & module >---------------
!=================================================
      subroutine chkncomp(ival,sbnm)
!=================================================
      character(10),parameter :: subnam='(chkncomp)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.ncomp ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      stop 999
      end subroutine chkncomp
!
!< #2. check consistency of values in argument & module >---------------
!=================================================
      subroutine chknrans(ival,sbnm)
!=================================================
      character(10),parameter :: subnam='(chknrans)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.nrans ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      stop 999
      end subroutine chknrans
!
!< #3. interpolate values at present time making use of time table >----
!=================================================
      subroutine setnow(time)
!=================================================
      real*8,intent(in) :: time
      integer :: nwbcnd
      if( nbcnd.lt.1 ) return
      call modutl_chkset('>',1,iset1,modnam//'(setnow)')
      nwbcnd=iwbcnd(nbcnd)
      call modutl_setnow(time,nbcnd,nbvals,nwbcnd,
     & iwbcnd,wdtbcnd,wdbcnd)
      end subroutine setnow
!
!< #4. return rotational matrix of periodic boundary >------------------
!=================================================
      subroutine set_rotmtrx(nbcnx,rot)
!=================================================
      character(13) :: subnam='(set_rotmtrx)'
      integer,intent(in)  :: nbcnx
      real*8 ,intent(out) :: rot(3,3,nbcnx)
      integer :: i,j,k,nb,kd
!
      call modutl_chkset('>',1,iset1,modnam//subnam)
!
      if(nbcnx.ne.nbcnd.and.NAME_I==1) then
        write(*,*) '### program error -1- ',modnam,subnam
        stop 999
      endif
      do 100 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      if((kd==kdprdc.and.idis(nb)==0).or.
     &   (kd==kdsld.and.idis(nb)==0)) then
        k=0
        do 101 j=1,3
        do 101 i=1,3
        k=k+1
        rot(i,j,nb)=wdbcnd(k,nb)
  101   continue
      endif
  100 continue
      end subroutine set_rotmtrx
!=================================================
      subroutine set_rotmtrx1(nbcnx,kd,nb,rot)
!=================================================
      implicit none
      character(13) :: subnam='(set_rotmtrx)'
      integer,intent(in)  :: nbcnx,kd,nb
      real(8) ,intent(out) :: rot(3,3)
      integer :: i,j,k
!
      call modutl_chkset('>',1,iset1,modnam//subnam)
!
      if(nbcnx.ne.nbcnd) then
        write(*,*) '### program error -1- ',modnam,subnam
        stop 'module_boundary/set_rotmtrx'
      endif
!      do 100 nb=1,nbcnd
!      kd=kdbcnd(0,nb)
      if(kd==kdprdc.or.kd==kdtchi.or.kd==kdsld) then
        k=0
        do 101 j=1,3
        do 101 i=1,3
        k=k+1
        rot(i,j)=wdbcnd(k,nb)
  101   continue
      endif
!  100 continue
      end subroutine set_rotmtrx1
!
!< #5. store domain data >----------------------------------------------
!=================================================
!      subroutine strdomain(fnam,ifle,nbdmn,nvbdmn,ncbdmn,
!     & nobdmn,ivbdmn,lvbdmn,icbdmn,lcbdmn,ihpc,ierr)
!=================================================
!      character(11),parameter :: subnam='(strdomain)'
!      character(*),intent(in)  :: fnam
!      integer     ,intent(in)  :: ifle,nbdmn,nvbdmn,ncbdmn,ihpc
!      integer     ,intent(in)  :: nobdmn(nbdmn)
!      integer     ,intent(in)  :: ivbdmn(0:nbdmn),lvbdmn(nvbdmn)
!      integer     ,intent(in)  :: icbdmn(0:nbdmn),lcbdmn(ncbdmn)
!      integer     ,intent(out) :: ierr
!      integer :: nob(nbcnd)
!      integer :: i,j,k,m,n,kve,kce,nb,nd
!
! --- 
!
!      if(ihpc.eq.0) then
!        call modutl_chkset('>',1,iset1,modnam//subnam)
!        call modutl_chkset('=',1,iset2,modnam//subnam)
!      endif
!
!      ierr=0
!      kve=0
!      kce=0
!
!      do 100 nb=1,nbcnd
!      k=0
!      do 101 m=1,nbdmn
!      if(nobdmn(m).eq.nobcnd(nb)) k=m
!  101 continue
!      nob(nb)=k
!      if(k.lt.1) goto 100
!      if(kdbinr(1,nb).gt.0 .and.
!     &    icbdmn(k-1).ge.icbdmn(k) ) goto 9001
!      kve=kve+(ivbdmn(k)-ivbdmn(k-1))
!      kce=kce+(icbdmn(k)-icbdmn(k-1))
!  100 continue
!
!      allocate( ivbcnd(0:nbcnd), lvbcnd(kve),
!     &          icbinr(0:nbcnd), lcbinr(kce),
!     &          stat=ierr )
!
!      if( ierr.ne.0 ) then
!        write(ifle,*) '### error : allocation failed'
!        goto 9999
!      endif
!
!      ivbcnd=0
!      lvbcnd=0
!      icbinr=0
!      lcbinr=0
!      i=0
!      j=0
!
!      do 200 nb=1,nbcnd
!      if(nob(nb).lt.1) goto 209
!      nd=nob(nb)
!      do 201 m=ivbdmn(nd-1)+1,ivbdmn(nd)
!      i=i+1
!      lvbcnd(i)=lvbdmn(m)
!  201 continue
!      do 202 m=icbdmn(nd-1)+1,icbdmn(nd)
!      j=j+1
!      lcbinr(j)=lcbdmn(m)
!  202 continue
!  209 continue
!      ivbcnd(nb)=i
!      icbinr(nb)=j
!  200 continue
!
!      return
! 9001 continue
!      write(ifle,*) '### error : data error'
!      write(ifle,*) 'the boundary domain has no additional ',
!     &              'information about side of the face'
!      write(ifle,*) 'boundary domain no. =',nobcnd(n)
!      write(ifle,*) 'name of boundary file =',fnam(:len_trim(fnam))
! 9999 continue
!      write(ifle,*) modnam//subnam
!      ierr=1
!      end subroutine strdomain
!
!
!< #6. input namelist >---------------------------------------
!=============================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,lsuf,idens,lcomp,
     &      ncmpx,nrnsx,nbcndx,nphase,lrans,E2P,lsld,ierror)
!=============================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(in)      :: ncmpx,nrnsx,idens
      logical,intent(in)      :: lrans,E2P,lcomp
      logical,intent(in)      :: lsld,lsuf
      integer,intent(out)     :: nbcndx
      integer,intent(inout)   :: nphase
      integer,intent(out)     :: ierror
!
! --- [namelist]
!
      integer,parameter :: msite=10,mbulk=10
      integer,parameter :: lc=10,lsp=20
      integer,parameter :: mbctim=10
      integer,parameter :: mcomp=200, mrans=50,mneq=300
      character(lc) :: kind,MHD_kind,potential_kind(mpotn),particle_kind
      character(lc) :: vel,temp,conc,rans,MHD_A,MHD_F
      character(lsp):: site_species_name(msite,mcomp),
     &                 bulk_species_name(mbulk,mcomp)
      character(lc) :: site_name(msite),bulk_name(mbulk)
      real(8) :: site_dens(msite),bulk_thick(mbulk)
      real(8) :: site_species_molefraction(msite,mcomp),
     &           bulk_species_molefraction(mbulk,mcomp),
     &           bulk_species_dens(mbulk,mcomp)
!
      integer :: no,side,ntime,profile,open_air=0,turn_wall=0,
     &           heat_release
      real    :: time(mbctim)
      real    :: p(mbctim),u(mbctim),v(mbctim),w(mbctim)
      real    :: t(mbctim),ys(mcomp*mbctim)
      real    :: u2(mbctim),v2(mbctim),w2(mbctim),htc2,mtc2,vsgm2(2)
      real    :: t2(mbctim),ys2(mcomp*mbctim)
      real    :: aks(mrans*mbctim)
      real    :: cycle,htc,mtc,vsgm(2),rot(4)
      real*8  :: rpm,end_x,end_y,end_z,begin_x,begin_y,begin_z
      integer :: start_rot_step,rotup_step
      integer :: suf_reac_no(mneq)
      character*80 :: name,name2
      character*1 ::axis
      real    :: xoffset,yoffset,zoffset
      integer :: discont=0
      real(8) :: roughness,massflux,Mach_nu
      real(8) :: FAI_mhd(2),AAA_mhd(3,2)
      integer :: flux_balance,fix_mach
      integer :: n_blade1,n_blade2
      integer :: porous_BC=0
      real(8) :: C2_porous,perm_porous,thickness_porous=0,porosity
      integer :: surface_reaction=0,thinwall
      real(8) :: Potential_value(mpotn),V_cell
      real(8) :: rademi=0.0,radabsorb=0.0,radtran=0.0
      integer :: radtype=1,radgpflag=1 
      real*8  :: VOF_ang
!
      namelist /boundary/ no,side,kind,ntime,MHD_kind,
     &                    vel,temp,conc,rans,
     &                    MHD_A,MHD_F,AAA_mhd,FAI_mhd,
     &                    time,profile,p,u,v,w,t,ys,aks,
     &                    cycle,htc,mtc,vsgm,rot
     &                    ,name,name2,axis,xoffset,yoffset,zoffset,
     &                    u2,v2,w2,t2,ys2,htc2,mtc2,vsgm2,
     &                    open_air,
     &                    turn_wall,
     &                    rpm,end_x,end_y,end_z,
     &                    begin_x,begin_y,begin_z,
     &                    start_rot_step,rotup_step,
     &                    site_species_name,
     &                    bulk_species_name,
     &                    site_name,site_dens,
     &                    bulk_name,
     &                    bulk_thick,
     &                    suf_reac_no,
     &                    site_species_molefraction,
     &                    bulk_species_molefraction,
     &                    bulk_species_dens,heat_release,
     &                    discont,roughness,
     &                    massflux,flux_balance,
     &                    fix_mach,Mach_nu,
     &                    n_blade1,n_blade2,
     &                    porous_BC,C2_porous,perm_porous,porosity,
     &                    thickness_porous,
     &                    surface_reaction,
     &                    Potential_kind,Potential_value,
     &                    particle_kind,V_cell,
     &			  rademi,radabsorb,radtran,radtype,radgpflag,
     &                    VOF_ang,
     &                    thinwall
!
!
! --- [local entities]
!
      integer,parameter :: iundef=-huge(1)
      real   ,parameter :: undef=-huge(1.)
      logical :: llsld=.false.
      integer :: icom,iq,jq,nst,nm,ncm,nbk,nsite,nbulk,ity
!
      character(lc),parameter ::
     &   periodic  ='periodic  ', 
     &   symmetric ='symmetric ',
     &   inlet     ='inlet     ', 
     &   outlet    ='outlet    ', 
     &   wall      ='wall      ',
     &   interface ='interface ',
     &   touchinlet='touchinlet',
     &   firewall  ='firewall  ',
     &   static    ='static    ',
     &   stagnation='stagnation',
     &   buffle    ='buffle    ',
     &   sliding   ='sliding   ',
     &   cvdsurface='cvdsurface',
     &   shutter   ='shutter   ',
     &   porous    ='porous    ',
     &   overset   ='overset    '
! 
      character(lc),parameter ::
     &   noslip=    'no-slip   ', 
     &   freeslip=  'free-slip ',
     &   rough=     'rough-wall',
     &   satawall=   'sata-wall',   !u+=1/k*ln[(Yp+Ks)/Ks]
     &   model=     'model     ',
     &   EWF=     'EWF     ',
     &   dirichlet= 'Dirichlet ', 
     &   neumann=   'Neumann   ',
     &   loglaw=    'log-law   ', 
     &   transfer=  'transfer  '

!
      character(lc),parameter :: 
     &    kndlst(16)=(
     &   /symmetric,periodic,interface,inlet,outlet,wall,
     &    touchinlet,firewall,static,stagnation,buffle,sliding,
     &    cvdsurface,shutter,porous,overset/)
      character(lc),parameter :: 
     &    vknd(7)=(
     &   /noslip,freeslip,loglaw,rough,model,EWF,satawall/)
      character(lc),parameter :: 
     &    xknd(4)=(
     &   /dirichlet,neumann,transfer,EWF/)
      character(lc),parameter :: 
     &    rknd(3)=(
     &   /dirichlet,neumann,loglaw/)
!
      integer :: i,j,n,kt,nb,nwbcnd
      integer :: kd,kvel,ktmp,kcnc,krns
      integer :: yfix,idum
      integer :: ios=0,ierr1,ierr2
!
      call modutl_chkset('=',1,iset1,modnam//subnam)
!
!-< 1. Initial set >-
!
      ierror=0
!
      ncomp=ncmpx
      nrans=nrnsx
!
      call nml_chksiz(ifle,'mcomp',mcomp,ncomp,
     & modnam//subnam,ierr1)
      call nml_chksiz(ifle,'mrans',mrans,nrans,
     & modnam//subnam,ierr2)
      if( ierr1.ne.0 .or. ierr2.ne.0 ) goto 9999
!
!--< 1.1 count up sizes >--
!
      nbcnd =0
      nwbcnd=0
      rewind ifli
  100 continue
      ntime=1
      read(ifli,boundary,iostat=ios)
      if(ios.lt.0) goto 101
!
      nbcnd=nbcnd+1
!
      if(undef_bcno.ne.0) then
        write(ifle,*)
     &' ERR: name [undefined] must be defined as lasted BC in ',
     &  cntlnam
        write(ifle,*)
     &' ERR: Or only one [undefined] BC can be defined in ',cntlnam
        stop 'ERR: &namelist/boundary'
      endif
      if(name.eq.'undefined') then
        undef_bcno=nbcnd
      endif
      call nml_chksiz(ifle,'mbctim',mbctim,ntime,
     & modnam//subnam,ierr1)
      call nml_errmsg0(ifle,ios,'boundary',ierr2)
      if( ierr1.ne.0 .or. ierr2.ne.0 ) then
        write(ifle,boundary)
        write(ifle,*) ' ERR: Sequence no. of the namelist =',nbcnd
        goto 9999
      endif
      nwbcnd=nwbcnd+max(1,ntime)
      goto 100
  101 continue
!
!--< 1.2 allocate arrays >--
!
      nbcndx=nbcnd
      if(nbcndx.lt.1) then
        write(ile,'(2a)') 'ERR: No bundary BC read in ',cntlnam
        stop
      endif
      iscomp=5
      israns=5+ncomp
      nbvals=max(israns+nrans,9)
      allocate( nobcnd(nbcnd),
     &          LBC_INDEX(0:nbcnd),
     &          kdbcnd(0:4,nbcnd),
     &          iwbcnd(0:nbcnd),
     &          adbcnd(4,nbcnd),
     &          wdbcnd(0:nbvals,nbcnd),
     &          wdtbcnd(0:nbvals,nwbcnd),
     &          boundName(nbcnd,2),
     &          prdcAxis(nbcnd),prdcOffset(nbcnd,3),
     &          conntBC(nbcnd),
     &          idis(nbcnd),nblade(nbcnd,2),thinw(nbcnd),
     &          kdbinr(2,nbcnd), 
     &		radprop(2,nbcnd),radwalltype(2,nbcnd),!jiang

     &          stat=ierr1 )
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error : allocation failed'
        goto 9999
      endif
      nobcnd =0
      kdbcnd =kxnone
      iwbcnd =0
      adbcnd =0.d0
      wdbcnd =0.d0
      wdtbcnd=0.d0
      kdbinr =0
      numbname=2
      prdcAxis=' '
      prdcOffset=0.d0
      idis(:)=0
      thinw(:)=0
      conntBC(:)=.false.
      nblade(:,:)=0
!      bcmat(:)=0
      radprop=0.d0
      radwalltype=0
!
!
! --- < 2. Input namelist >-
!
! --- < 2.1 read namelist >--
!
      kind = ' '
!
      nb=0
      kt=0
      rewind ifli
! --- start read namelist file for every BC
 1000 continue
      vel   = ' '
      temp  = ' '
      conc  = ' ' 
      rans  = ' ' 
      no    = iundef
      side  = iundef
      ntime = 1
      time  = 0.d0
      cycle = undef
      profile=iundef
      p     = undef
      u     = undef
      v     = undef
      w     = undef
      t     = undef
      ys    = undef
      aks   = undef
      htc   = undef
      mtc   = undef
      vsgm  = undef
      rot   = 0.
      name  = ' '
      name2 = ' '
      axis  = ' '
      u2    = undef
      v2    = undef
      w2    = undef
      t2    = undef
      ys2   = undef
      htc2   = undef
      mtc2   = undef
      vsgm2  = undef
      xoffset = undef
      yoffset = undef
      zoffset = undef
      open_air=0
      thinwall=0
      turn_wall=0
      rpm=0.d0
      end_x=0.d0
      end_y=0.d0
      end_z=0.d0
      begin_x=0.d0
      begin_y=0.d0
      begin_z=1.d0
      start_rot_step=0
      rotup_step=0
!
! --- surface reaction
!
      site_species_name=''
      bulk_species_name=''
      site_name=''
      bulk_name=''
      site_dens=undef
      suf_reac_no=iundef
      num_site=1.d0
      discont=0
      n_blade1=0
      n_blade2=0
!
      surface_reaction=0
!
      read(ifli,boundary,iostat=ios)
      if(ios.lt.0 ) goto 1001
       call nml_errmsg0(ifle,ios,'boundary',ierr1)
      if( ierr1.ne.0 ) goto 9999
      nb=nb+1
!
! --- < 2.2 common procedure >--
!
! --- / global kind /
!
      call nml_listno(16,kndlst,kind,kd)
      call nml_chkchr0(ifle,'kind',kind,kd,ierr1)
      if( ierr1.ne.0 ) goto 9999
      kind=kndlst(kd)
!
      if( kind.eq.periodic.or. kind.eq.symmetric ) ntime=1
      kt=kt+ntime
      iwbcnd(nb)=kt
!
! --- / domain no. & side no. /
!
      if(kind.eq.periodic)  call chkintn('side',side)
      if(kind.eq.interface) call chkintn('side',side)
      if(ierror.ne.0 ) goto 9999
!
      if(no.eq.iundef) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : no'
        goto 9999
      endif
!
      do 200 n=1,nb-1
      if(no.ne.nobcnd(n)) goto 200
      ierr1=0
      if(side.eq.iundef) then
        ierr1=n
      elseif(kdbinr(1,n).gt.0) then
        ierr1=kdbinr(1,n)
      else
        kdbinr(1,n)=nb
        kdbinr(1,nb)=n
      endif
      if(ierr1.ne.0) then
        write(ifle,*) '  #### error : data error'
        write(ifle,*) 'no =',no
        write(ifle,*) 'this no. has been specified in previous line'
        write(ifle,*) 'sequence no. of the line =',ierr1
        goto 9999
      endif
  200 continue
      nobcnd(nb)=no
      boundName(nb,1)=name
      if(name.eq.' ') then
        write(ifle,*) ' EER: Boundary no. ',no,' must be specified'
      write(ifle,*) ' Please input the name corresponding to grid file'
        goto 9999
      endif
      boundName(nb,2)=name2
      prdcAxis(nb)=axis

      prdcOffset(nb,1)=dble(xoffset)
      prdcOffset(nb,2)=dble(yoffset)
      prdcOffset(nb,3)=dble(zoffset)
!
      if(side.ne.iundef) then
        if( side.lt.1 .or. side.gt.2 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'side =',side
          write(ifle,*) 'it must be 1 or 2'
          goto 9999
        endif
        if( kdbinr(1,nb).gt.0 ) then
          if( kdbinr(2,kdbinr(1,nb))+side.ne.3 ) then
            write(ifle,*) '### error : data error'
            write(ifle,*) 'side =',side
            if( side.eq.1 ) write(ifle,*) 'it must be 2'
            if( side.eq.2 ) write(ifle,*) 'it must be 1'
            goto 9999
          endif
        endif
        kdbinr(2,nb)=side
      endif
      if(kind/=outlet) then
        if(open_air/=0) then
        write(ifle,*) 
     &    ' ERR: open_air MUST be 0, While NOT [outlet] BC in ',
     &     cntlnam
          stop 'ERR: re-set open_air=0 '
        endif
      endif
      if(kind/=wall.and.
     &   kind/=cvdsurface.and.
     &   kind/=buffle.and.
     &   kind/=firewall.and.
     &   kind/=interface.and.
     &   kind/=outlet
     &   ) then
        if(turn_wall/=0) then
          write(ifle,*) 
     &  ' ERR: turn_wall MUST be 0, While NOT [wall] BC in ',cntlnam
          stop 'ERR: re-set turn_wall=0'
        endif
      endif
!
      radprop(1,nb)=rademi   !jiang
      radprop(2,nb)=radabsorb
      radwalltype(1,nb)=radtype
      radwalltype(2,nb)=radgpflag
!
      if( kind.eq.symmetric ) goto 1101
      if( kind.eq.periodic )  goto 1102
      if( kind.eq.inlet )     goto 1103
      if( kind.eq.outlet )    goto 1103
      if( kind.eq.touchinlet) goto 1104
      if( kind.eq.firewall)   goto 1105
      if( kind.eq.buffle)     goto 1105
      if( kind.eq.static)     goto 1106
      if( kind.eq.stagnation) goto 1106
      if( kind.eq.sliding)    goto 1107
      if( kind.eq.cvdsurface) goto 1105
      if( kind.eq.porous    ) goto 1105
      if( kind.eq.overset   ) goto 1108
      goto 1105
!
!--< 2.3 symmetric boundary >--
!
 1101 continue
      kdbcnd(0,nb)=kdsymm
      goto 1000
!
!--< 2.4 periodic boundary >--
!
 1102 continue
!
      if(axis==' ') then
        write(*,*) 'ERR: NOT define axis for periodic boundary nb=',nb
        write(*,*) 'MSG: axis="x", axis="y", axis="z" or axis="r"'
        stop 'STOP AT: &boundary/axis=""'
      elseif(axis=='x'.or.axis=='X'.or.
     &       axis=='y'.or.axis=='Y'.or.
     &       axis=='z'.or.axis=='Z'.or.
     &       axis=='r'.or.axis=='R') then
        if(discont==0) then
          if(xoffset==undef.or.yoffset==undef.or.zoffset==undef) then
            write(*,*) 'ERR: NOT define xoffset,yoffset,zoffset nb=',nb
            write(*,*) 'MSG: define them for periodic boundary'
            write(*,*) 'MSG: NOT necessary define them if discont=1 '
            stop 'STOP AT: &boundary'
          endif
        elseif(discont==1) then
          if(xoffset==undef.or.yoffset==undef.or.zoffset==undef) then
            write(*,*) 'WRN: NOT define xoffset,yoffset,zoffset nb=',nb
          endif
        endif
      else
        write(*,*) 
     &   'MSG: axis="x", axis="y", axis="z" or axis="r", nb=',nb
        stop 'STOP AT: &boundary'
      endif
!
      if(name2.eq.' ') then
        write(ifle,*) ' EER: Boundary no. ',no,' is periodic boundary',
     &  ' its counter-BCs name must be defined'
      write(ifle,*) ' Please input the name corresponding to grid file'
        goto 9999
      endif
      kdbcnd(0,nb)=kdprdc
      call romtrx(nb)
      if(discont==0.or.discont==1) then
        idis(nb)=discont
      else
         write(ifle,'(a)') 'EER: [discont] MUST be 0 or 1'
         write(ifle,'(a)') 
     &     'MSG: discont=0 : continuity periodic BC'
         write(ifle,'(a)') 
     &     'MSG: discont=1 : discontinuity periodic BC'
        stop'module_boundary/inputdata'
      endif
      goto 1000
!---------------------------------
!--< 2.5 inlet/outlet boundary >--
!---------------------------------
 1103 continue
      call chkntime
      call chkflt1(0,ntime,'u',u)
      call chkflt1(0,ntime,'v',v)
      call chkflt1(0,ntime,'w',w)
!
      if(open_air==5) then 
        call chkflt1(0,ntime,'t',t)
      else
        call chkflt1(1,ntime,'t',t)
      endif
      call chkflt2(1,ncomp,ntime,'ys',ys)
      if(E2P) then
        call chkflt1(0,ntime,'u2',u2)
        call chkflt1(0,ntime,'v2',v2)
        call chkflt1(0,ntime,'w2',w2)
        if(open_air==5) then
          call chkflt1(0,ntime,'t2',t2)
        else
          call chkflt1(1,ntime,'t2',t2)
        endif
        call chkflt2(1,ncomp,ntime,'ys2',ys2)
      endif
!
      if(kind.eq.inlet) then
!
        kdbcnd(0,nb)=kdilet
!        if(idens/=4) then
        if(flux_balance/=1) then
          call chkfltn('p(1)',p(1))
        endif
!
        if(profile.eq.iundef) then
          write(ifle,*) '### Inlet BC error : data error'
          write(ifle,*) 'Define profile=0,1 ro 2  in Inlet BC'
          write(ifle,*) 'profile=0: Uniform Vel. Profile'
          write(ifle,*) 'profile=1: User defined Vel. Profile'
          write(ifle,*) 'profile=2: Developed Vel. profile'
          goto 9999
        elseif(profile.eq.0) then
          distrb=0
        elseif(profile.eq.1) then
          distrb=1
        elseif(profile.eq.2) then
          distrb=2
        else
          write(ifle,*) 'Must be profile=0,1 or 2'
          goto 9999
        endif
      else
        kdbcnd(0,nb)=kdolet
        call chkflt1(1,ntime,'p',p)
      endif
      if(lrans) then
        call chkflt2(1,nrans,ntime,'aks',aks)
      endif
!
      if(kind.eq.inlet.and.vsgm(1).ne.undef) then
        call chkflt3('vsgm(1)',vsgm(1))
        call chkflt3('vsgm(2)',vsgm(2))
      endif
      if(E2P) then
        if(kind.eq.inlet.and.vsgm2(1).ne.undef) then
          call chkflt3('vsgm2(1)',vsgm2(1))
          call chkflt3('vsgm2(2)',vsgm2(2))
        endif
      endif
!
      call strdata(nb,1)
!
!
      if(ierror.ne.0) goto 9999
      goto 1000
!
!
!--< 2.? touch boundary >--
!
 1104 continue
!
      if(axis==' ') then
        write(*,*) 'ERR: NOT define axis for periodic boundary nb=',nb
        write(*,*) 'MSG: axis="x", axis="y", axis="z" or axis="r"'
        stop 'STOP AT: &boundary/axis=""'
      elseif(axis=='x'.or.axis=='X'.or.
     &       axis=='y'.or.axis=='Y'.or.
     &       axis=='z'.or.axis=='Z'.or.
     &       axis=='r'.or.axis=='R') then
      else
        write(*,*) 
     &   'MSG: axis="x", axis="y", axis="z" or axis="r", nb=',nb
        stop 'STOP AT: &boundary'
      endif
!
      if(name2.eq.' ') then
        write(ifle,*) ' EER: Boundary no. ',no,
     &  ' is touch-inlet boundary',
     &  ' its touching BC name must be specified by name2'
        write(ifle,*) ' Input the name corresponding to grid file'
        goto 9999
      endif     
!
      kdbcnd(0,nb)=kdtchi
!
      if(discont==0.or.discont==1) then
        idis(nb)=discont
      else
         write(ifle,'(a)') 'EER: [discont] MUST be 0 or 1'
         write(ifle,'(a)') 
     &            'MSG: discont=0 : continuity touch-inlet BC'
         write(ifle,'(a)') 
     &            'MSG: discont=1 : discontinuity touch-inlet BC'
        stop'module_boundary/inputdata'
      endif

!
      goto 1000
!
!--< 2.6 other boundary >-- wall BC (& interface BC)
!
 1105 continue
!
      thinw(nb)=0 
      if(thinwall==1) then 
        thinw(nb)=1 
      elseif(thinwall==0) then 
        thinw(nb)=0 
      else
        stop 'ERR: SHOULD: thinwall=0, or 1'
      endif 
!
      if(kind.eq.firewall) then
        if(surface_reaction==1) then
          call surfreac(nphase,nsite,msite,nbulk,mbulk,
     &                site_name,site_dens,
     &                bulk_name,lc,KSUF
     &                )
        endif
        call nml_listno(4,xknd,conc,kcnc)
        if( kcnc.gt.0 ) conc=xknd(kcnc)
        if(conc/=neumann) then
          write(ifle,'(1x,a,I4)') 
     &  " ### EER: firewall must be use [conc='Neumann'] on BC no=",
     &    no
          stop ': Stop at module_boundary'
        endif
        kdbcnd(0,nb)=kdfire
        call chkntime
        call chkflt1(0,ntime,'t',t)
        call chkfire(0,ncomp,ntime,'ys',ys)
        call chkfltn('p(1)',p(1))
      endif
!
      if(kind.eq.cvdsurface) then
        call nml_listno(4,xknd,temp,ktmp)
        if(ktmp.gt.0) temp=xknd(ktmp)
        if(temp/=neumann) then
!          write(ifle,*)
!     &  " ### EER: cvdsurface must be use [temp='Neumann']"
!          stop ': Stop at module_boundary'
        endif
        call nml_listno(4,xknd,conc,kcnc)
        if( kcnc.gt.0 ) conc=xknd(kcnc)
        if(conc/=neumann) then
          write(ifle,*)
     &  " ### EER: cvdsurface must be use [conc='Neumann']"
          stop ': Stop at module_boundary'
        endif
        kdbcnd(0,nb)=kdcvd
	KSUF=1
        call chkntime
!        call chk_t_neumann(1,ntime,'t',t)
        call chkfire(0,ncomp,ntime,'ys',ys)
        call chkfltn('p(1)',p(1))
!
        call surfreac(nphase,nsite,msite,nbulk,mbulk,
     &                site_name,site_dens,
     &                bulk_name,lc,KSUF
     &                )
!
!        nsite=0
!        do nst=1,msite
!        IF(site_name(nst)/=' ') then
!          nsite=nsite+1
!          nphase=nphase+1
!          if(site_dens(nst)==undef) then
!            write(ifle,'(a,I4)') 
!     &      'ERR* SITE density NOT defined at BC no=',no
!            stop 'ERR: SITE density [site_dens] must be defined'
!          endif
!        endif
!        enddo
!        nbulk=0
!        do nbk=1,mbulk
!        if(bulk_name(nbk)/=' ') then
!          nbulk=nbulk+1
!          nphase=nphase+1
!        endif
!        enddo        
      endif
!
      if(kind.eq.buffle) then
        if(surface_reaction==1) then
          call surfreac(nphase,nsite,msite,nbulk,mbulk,
     &                site_name,site_dens,
     &                bulk_name,lc,KSUF
     &                )
        endif
        if(discont/=0) stop 'ERR: Buffle BC only support discont=0'
        kdbcnd(0,nb)=kdbuff
        if(name2==' ') then 
          write(ifle,*) ' EER: Boundary no. ',no,
     &    ' is [buffle] boundary',
     &  ' pair BC name specification by [name2] is under-construction'
          write(ifle,*) 
     & ' Unset [name2] in [',cntlnam,']May be, you can not use buffle.'
          stop
        end if
        goto 1000 
      end if
!
      if(kind.eq.porous) then
        write(ifle,*) ' ERR: NOT Support [porous] BC.'
        stop 'MSG: Using [buffle] BC for porous BC'
        if(surface_reaction==1) then
          call surfreac(nphase,nsite,msite,nbulk,mbulk,
     &                site_name,site_dens,
     &                bulk_name,lc,KSUF
     &                )
        endif
        if(discont/=0) stop 'ERR: Buffle BC only support discont=0'
        kdbcnd(0,nb)=kdpors
        if(name2==' ') then 
          write(ifle,*) ' EER: Boundary no. ',no,
     &    ' is [porous] boundary',
     &  ' pair BC name specification by [name2] is under-construction'
          write(ifle,*) 
     & ' Unset [name2] in [',cntlnam,']May be, you can not use buffle.'
          stop
        end if
        if(C2_porous==undef.or.perm_porous==undef.or.
     &    thickness_porous==undef) then
          write(ifle,'(1X,A,1X,A)') 
     &  'ERR: NOT define C2_porous, perm_porous  and thickness_porous',
     &  "for Darcy's law"
          stop 'ERR: NOT define C2_porous, perm_porous'
        endif
        goto 1000 
      end if
!
!
      if(kind==shutter) then
        if(discont/=0) stop 'ERR: shutter BC only support discont=0'
        kdbcnd(0,nb)=kdshutr
        if(name2==' ') then 
          write(ifle,*) ' EER: Boundary no. ',no,
     &    ' is [buffle] boundary',
     &  ' pair BC name specification by [name2] is under-construction'
          write(ifle,*) 
     & ' Unset [name2] in [',cntlnam,']May be, you can not use buffle.'
          stop
        end if
      end if
!
      if(kind.eq.interface) then
        if(surface_reaction==1) then
          call surfreac(nphase,nsite,msite,nbulk,mbulk,
     &                site_name,site_dens,
     &                bulk_name,lc,KSUF
     &                )

        endif
        kdbcnd(0,nb)=kdintr
        if(name2.eq.' ') then
          write(ifle,*) ' EER: Boundary no. ',no,
     &    ' is [interface] boundary',
     &    ' its pair BC name must be specified by [name2]'
          write(ifle,*) 
     &    ' Set [name2] in [',cntlnam,
     &    '] , corresponding to grid file'
          goto 9999
        elseif(name2=='dummy') then
          dum_bc=nb
          if(discont==1) then
            write(ifle,'(1X,a,I4)') 
     &      'EER: [discont=1] NOT support name2="dummy" BC no= ',no
            stop 'ERR-MSG: define name2 for interface BC'
          endif
        endif  
        if(discont==0.or.discont==1) then
          idis(nb)=discont
!          if(discont==1) then
!            stop 'ERR: NOT support discontinuity interface BC'
!          endif
        else
          write(ifle,'(a)') 'EER: [discont] MUST be 0 or 1'
          write(ifle,'(a)') 
     &     'MSG: discont=0 : continuity periodic BC'
          write(ifle,'(a)') 
     &     'MSG: discont=1 : discontinuity periodic BC'
          stop 'module_boundary/inputdata'
        endif
      endif
!
!/ kind of each variable /
!
      call nml_listno(7,vknd,vel ,kvel)
!      call nml_chkchr0(ifle,'vel',vel,kvel,ierr1)
!
      call nml_listno(4,xknd,temp,ktmp)
      call nml_chkchr0(ifle,'temp',temp,ktmp,ierr1)
!
      call nml_listno(4,xknd,conc,kcnc)
!      call nml_chkchr0(ifle,'conc',conc,kcnc,ierr1)
!
      call nml_listno(3,rknd,rans,krns)
!      call nml_chkchr0(ifle,'rans',rans,krns,ierr1)
!

      if( ierr1.ne.0 ) goto 9999
      if( lrans ) then
      else
        rans=' '
      endif
!
      if( kvel.gt.0 ) vel =vknd(kvel)
      if( ktmp.gt.0 ) temp=xknd(ktmp)
      if( kcnc.gt.0 ) conc=xknd(kcnc)
      if( krns.gt.0 ) rans=rknd(krns)
!
	
      if( vel.eq.noslip )     kdbcnd(1,nb)=kvnslp
      if( vel.eq.freeslip )   kdbcnd(1,nb)=kvfslp
      if( vel.eq.loglaw )     kdbcnd(1,nb)=kvlglw
      if( vel.eq.rough  )     kdbcnd(1,nb)=kvrogf
      if( vel.eq.satawall  )  kdbcnd(1,nb)=kvsataW
      if( vel.eq.model  )     kdbcnd(1,nb)=kvmodl
      if( vel.eq.EWF  )     then
!        call FFRABORT(1,'ERR: vel=EWF is NOT support')
        kdbcnd(1,nb)=kvEWF
      endif
!
      if( temp.eq.dirichlet ) kdbcnd(2,nb)=ktdirc
      if( temp.eq.neumann )   kdbcnd(2,nb)=ktneum
      if( temp.eq.transfer )  kdbcnd(2,nb)=kttrns
      if( temp.eq.EWF )  then
!         call FFRABORT(1,'ERR: temp=EWF is NOT support')
         kdbcnd(2,nb)=ktEWF
      endif
!
      if( conc.eq.dirichlet ) kdbcnd(3,nb)=kydirc
      if( conc.eq.neumann )   kdbcnd(3,nb)=kyneum
      if( conc.eq.transfer )  kdbcnd(3,nb)=kytrns
      if( conc.eq.EWF )  then
!        call FFRABORT(1,'ERR: conc=EWF is NOT support')
        kdbcnd(3,nb)=kyEWF
      endif

!
      if( rans.eq.dirichlet ) kdbcnd(4,nb)=kkdirc
      if( rans.eq.neumann )   kdbcnd(4,nb)=kkneum
      if( rans.eq.loglaw )    kdbcnd(4,nb)=kklglw
!
!/ check & store data /
!
      if( kind.eq.interface ) then
        if( temp.eq.dirichlet ) then
          write(*,*) '### error : data error'
          write(*,*) 'temp =',temp
          write(*,*) 'it can not be = ',temp
          write(*,*) 'nb=' ,nb
          goto 9999
        endif
      endif
!
      call chkntime
!
      if(vel.eq.noslip.or.vel.eq.loglaw.or.vel==rough.or.vel==model
     & .or.vel==satawall) then
        if(turn_wall/=1) then
          if(u(1)==undef.or.v(1)==undef.or.w(1)==undef) then
            write(ifll,'(a)') 
     &'MSG: Velocity must be defined while use no-slip or log-law wall'
          endif
          call chkflt1(0,ntime,'u',u)
          call chkflt1(0,ntime,'v',v)
          call chkflt1(0,ntime,'w',w)
          if(E2P) then
            call chkflt1(0,ntime,'u2',u2)
            call chkflt1(0,ntime,'v2',v2)
            call chkflt1(0,ntime,'w2',w2)
!            call chkflt1(1,ntime,'t2',t2)
!            call chkflt2(1,ncomp,ntime,'ys2',ys2)
          endif
        endif
      else   ! slip wall BC
        call chkfltn('u(1)',u(1))
        call chkfltn('v(1)',v(1))
        call chkfltn('w(1)',w(1))
      endif
!
      if(temp==dirichlet.or.temp==transfer.or.temp==EWF) then
        idum=1
      elseif( temp.eq.neumann ) then
        idum=0
      else
        idum=-1
      endif
      if( kind.eq.interface .and.temp/=neumann) idum=-1
      if( idum.lt.0 ) then
        call chkfltn('t(1)',t(1))
        if(E2P) then
          call chkfltn('t2(1)',t2(1))
        endif
      else
        call chkflt1(idum,ntime,'t',t)
        if(E2P) then
          call chkflt1(idum,ntime,'t2',t2)
        endif
      endif
!
      if( conc.eq.dirichlet .or. conc.eq.transfer ) then
        yfix=1
      elseif( conc.eq.neumann ) then
        yfix=0
      else
        yfix=-1
      endif
      if( yfix.lt.0 ) then
        call chkfltn('ys(1)',ys(1))
        if(E2P) then
          call chkfltn('ys2(1)',ys2(1))
        endif
      else
        call chkflt2(yfix,ncomp,ntime,'ys',ys)
        if(E2P) then
          call chkflt2(yfix,ncomp,ntime,'ys2',ys2)
        endif
      endif
!
      if( rans.eq.dirichlet ) then
        call chkflt2(1,nrans,ntime,'aks',aks)
      elseif( rans.eq.neumann ) then
!        call chkfltn('aks(1)',aks(1))
        aks=0.
      elseif( rans.eq.loglaw ) then
        call chkfltn('aks(1)',aks(1))
      endif
!
      if( temp.eq.transfer ) call chkflt3('htc',htc)
      if( conc.eq.transfer ) call chkflt3('mtc',mtc)
      if(E2P) then
        if( temp.eq.transfer ) call chkflt3('htc2',htc2)
        if( conc.eq.transfer ) call chkflt3('mtc2',mtc2)
      endif
!
      call strdata(nb,yfix)
      if( ierror.ne.0 ) goto 9999
      goto 1000
!
! --- pressure and stagnation pressure
!
 1106 continue
      call chkntime
!      call chkflt1(0,ntime,'u',u)
!      call chkflt1(0,ntime,'v',v)
!      call chkflt1(0,ntime,'w',w)
      call chkflt1(1,ntime,'t',t)
      call chkflt2(1,ncomp,ntime,'ys',ys)
      call chkflt1(1,ntime,'p',p)
!
      if(kind.eq.static) then
        kdbcnd(0,nb)=kdpres
      elseif(kind.eq.stagnation) then
        kdbcnd(0,nb)=kdstag
      endif
      call strdata(nb,1)
      goto 1000
!
 1107 continue
!
      kdbcnd(0,nb)=kdsld
      lsldbc=.true.
      if(n_blade1<0) then
        write(ifle,*) ' ERR:n_blade1 must be lager then 0'
        stop
      elseif(n_blade1/=0) then
        nblade(nb,1)=n_blade1
      endif
!
      if(n_blade2<0) then
        write(ifle,*) ' ERR:n_blade2 must be lager then 0'
        stop
      elseif(n_blade2/=0) then
        nblade(nb,2)=n_blade2
      endif
!
      llsld=.true.
      if(name2.eq.' ') then
        write(ifle,*) ' EER: Boundary no. ',no,
     &  ' is [sliding] boundary,',
     &  ' In Sliding BC,2nd BC must be specified by name2'
        write(ifle,*) ' Input the name corresponding to grid file'
        stop
      endif
      call romtrx(nb)
!
      if(discont==1) then
        idis(nb)=discont
      elseif(discont==0) then
        write(ifle,'(a)') 'EER: [discont] MUST be 1 for sliding BC'
        stop'module_boundary/inputdata'
      elseif(discont==2) then
      else
         write(ifle,'(a)') 'EER: [discont] MUST be 0 or 1'
         write(ifle,'(a)') 
     &   'MSG: discont=0 : continuity sliding BC at initial condition'
         write(ifle,'(a)') 
     &'MSG: discont=1 : discontinuity sliding BC at initial condition'
        stop'module_boundary/inputdata'
      endif
      goto 1000
!
 1108 continue
!
      lovstbc=.true.
      if(kind==overset) then
        kdbcnd(0,nb)=kdovst
      endif
      goto 1000
!
 1001 continue
!
!
!-< 3. check for inner boundary >-
!
      do 300 n=1,nbcnd
      if( kdbinr(2,n).gt.0 .and. kdbinr(1,n).lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of namelist'
        write(ifle,*) 'no =',nobcnd(n)
        write(ifle,*) 'side =',3-kdbinr(2,n)
        goto 9999
      endif 
  300 continue
!
!-< 4. Initialize "wdbcnd" >-
!
      do 400 n=1,nbcnd
      kt=iwbcnd(n-1)+1
      do 401 i=1,nbvals
      wdbcnd(i,n)=wdtbcnd(i,kt)
  401 continue
  400 continue
!
      if(llsld.and.(.not.lsld)) then
        write(ifle,'(a)') 'WRN: Not define [rpm] '
        write(ifle,'(a)') 
     &  'WRN: [rpm] should be defined while using sliding BC'
      elseif((.not.llsld).and.lsld) then
        write(ifle,'(a)') 
     &  'WRN: Sliding BC must be defined while rpm is defined'
!        write(ifle,*)
!     &   'ERR: Sliding BC must be defined while rpm is defined'
!        write(ifle,*) 
!     &   'MSG: Please set [rpm=0.d0] or remove [rpm] in ',cntlnam
!        goto 9999
      endif
!
      if(lsuf.and.KSUF==0) then
        write(ifle,'(3a)') 'ERR: Surface reaction has been defined',
     &   " ,but no [kind='cvdsurface'] in &boundary in",cntlnam
!        stop 'at &boundary'
      endif

!
      return
!
 9999 continue
      write(ifle,*) modnam//subnam
      ierror=1
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #6.1 check real data >
!=================================================
      subroutine chkflt1(isw,nt,dnam,dval)
!=================================================
      integer     ,intent(in) :: isw,nt
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval(nt)
      integer :: i
      do 100 i=1,nt
      if( dval(i).eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) 'lack of data : ',dnam,'(i)'
        write(ifle,*) 'i =',i
        goto 9999
      endif
  100 continue
      if( isw.eq.0 ) return
      do 101 i=1,nt
      if( ( isw.eq.1 .and. dval(i).le.0.d0 ) .or. 
     &    ( isw.eq.2 .and. dval(i).lt.0.d0 ) ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) dnam,'(i) = ',dval(i)
      write(ifle,*) 'i =',i
      if( isw.eq.1 ) write(ifle,*) 'it must be > 0'
      if( isw.eq.2 ) write(ifle,*) 'it must be >= 0'
      goto 9999
      endif
  101 continue
      return
 9999 continue
      ierror=1
      end subroutine chkflt1
!
!=================================================
      subroutine chk_t_neumann(isw,nt,dnam,dval)
!=================================================
      implicit none
      integer     ,intent(in) :: isw,nt
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval(nt)
      integer                 :: i
      real*8,parameter        :: SML=1.d-15
!
      do 100 i=1,nt
      if( dval(i).eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) 'lack of data : ',dnam,'(i)'
        write(ifle,*) 'i =',i
        goto 9999
      endif
  100 continue
!
      if( isw.eq.0 ) return
!
      do 101 i=1,nt
      if(abs(dval(i)).gt.SML) then
        write(ifle,*) ' ### ERR: data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*)  dnam,'(i) = ',dval(i),'i =',i
        write(ifle,*) 'it must be = 0'
        goto 9999
      endif
  101 continue
      return
 9999 continue
      ierror=1
      end subroutine chk_t_neumann
!
!< #6.2 check real data >
!=================================================
      subroutine chkflt2(isw,nc,nt,dnam,dval)
!=================================================
      integer     ,intent(in) :: isw,nc,nt
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval(nc,nt)
      integer :: i,j
      do 100 i=1,nc
      do 100 j=1,nt
      if( dval(i,j).eq.undef ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) 'lack of data : ',dnam,'(i,j)'
      write(ifle,*) 'i,j =',i,j
      goto 9999
      endif
  100 continue
      if( isw.eq.0 ) return
      do 101 i=1,nc
      do 101 j=1,nt
      if( dval(i,j).lt.0.d0 ) then
      write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) dnam,'(i,j) = ',dval(i,j)
      write(ifle,*) 'i,j =',i,j
      write(ifle,*) 'it must be >= 0'
      goto 9999
      endif
  101 continue
      return
 9999 continue
      ierror=1
      end subroutine chkflt2
!
!< #6.3 check real data >
!=================================================
      subroutine chkflt3(dnam,dval) 
!=================================================
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval
      if( dval.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) 'lack of data : ',dnam
        ierror=1
      elseif( dval.le.0. ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) dnam,' = ',dval
        write(ifle,*) 'it must be > 0'
        ierror=1
      endif
      end subroutine chkflt3
!
!< #6.4 check real data >
!=================================================
      subroutine chkfltn(dnam,dval)
!=================================================
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval
      if( dval.eq.undef ) return
      write(ifle,*) '### error : data error'
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) dnam,' = ',dval
      write(ifle,*) 'it can not be specified'
      ierror=1
      end subroutine chkfltn
!
!< #6.5 check integer data >
!=================================================
      subroutine chkintn(dnam,dval)
!=================================================
      character(*),intent(in) :: dnam
      integer     ,intent(in) :: dval
!
      if( dval.eq.iundef ) return
      write(ifle,*) '### error : data error'
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) dnam,' = ',dval
      write(ifle,*) 'it can not be specified'
      ierror=1
      end subroutine chkintn

!=================================================
      subroutine chkfire(isw,nc,nt,dnam,dval)
!=================================================
      implicit none
      integer     ,intent(in) :: isw,nc,nt
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval(nc,nt)
      integer :: i,j
      do 100 i=1,nc
      do 100 j=1,nt
      if( dval(i,j).eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) 'lack of data : ',dnam,'(i,j)'
        write(ifle,*) 'i,j =',i,j
        goto 9999
      endif
  100 continue
      if(isw/=0 ) return
      do 101 i=1,nc
      do 101 j=1,nt
      if( dval(i,j).ne.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*)  dnam,'(i,j) = ',dval(i,j)
        write(ifle,*) 'i,j =',i,j
        write(ifle,*) 'it must be = 0'
        goto 9999
      endif
  101 continue
      return
 9999 continue
      ierror=1
      end subroutine chkfire
!
!< #6.6 check no. of time tables >
!=================================================
      subroutine chkntime
!=================================================
      integer :: i
      real    :: dum1
      if( ntime.lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) 'ntime =',ntime
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
      do 210 i=2,ntime
      if( time(i-1).ge.time(i) ) then
      write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) 'time(i-1),time(i) = ',time(i-1),time(i)
      write(ifle,*) 'time(i-1) must be < time(i)'
      write(ifle,*) 'i =',i
      goto 9999
      endif
  210 continue
      dum1=time(ntime)-time(1)
      if( cycle.ne.undef .and. cycle.le.dum1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) 'cycle = ',cycle
        write(ifle,*) 'it must be > time(ntime)-time(1)'
        write(ifle,*) 'time(ntime)-time(1) = ',dum1
        goto 9999
      endif
      return
 9999 continue
      ierror=1
      end subroutine chkntime
!
!< #6.7 store data >
!=================================================
      subroutine strdata(nb,yfix)
!=================================================
      integer,intent(in) :: nb,yfix
      integer :: i,j,i0,j0,k0
      real*8  :: ysw(mcomp),sum,fmax
                           wdbcnd(0,nb)=-1.d0
      if( cycle.ne.undef ) wdbcnd(0,nb)=dble(cycle)
      if( htc  .ne.undef ) adbcnd(1,nb)=dble(htc)
      if( mtc  .ne.undef ) adbcnd(2,nb)=dble(mtc)
      adbcnd(3,nb)=-1.d0
      if( vsgm(1).ne.undef ) then
        adbcnd(3,nb)=dble(vsgm(1))
        adbcnd(4,nb)=dble(vsgm(2))
      endif
      i0=iwbcnd(nb-1)
      do 100 i=1,ntime
                          wdtbcnd(0,i+i0)=dble(time(i))
      if( p(1).ne.undef ) wdtbcnd(1,i+i0)=dble(p(i))
      if( u(1).ne.undef ) wdtbcnd(2,i+i0)=dble(u(i))
      if( v(1).ne.undef ) wdtbcnd(3,i+i0)=dble(v(i))
      if( w(1).ne.undef ) wdtbcnd(4,i+i0)=dble(w(i))
      if( t(1).ne.undef ) wdtbcnd(5,i+i0)=dble(t(i))
      if( ys(1).eq.undef ) goto 109
      j0=ncomp*(i-1)
      sum=0.d0
      if( yfix.gt.0 ) then
        do 200 j=1,ncomp
        ysw(j)=dble(max(0.,ys(j+j0)))
        sum=sum+ysw(j)
  200   continue
        if( sum.le.0.d0 ) goto 9001
        sum=1.d0/sum
      else
        fmax=-1.d0
        do 210 j=1,ncomp
        ysw(j)=dble(ys(j+j0))
        sum=sum+ysw(j)
        if( abs(ysw(j)).gt.fmax ) then
          fmax=abs(ysw(j))
          k0=j
        endif
  210   continue
        if( abs(sum).gt.1.d-4 ) goto 9001
        ysw(k0)=ysw(k0)-sum
        sum=1.d0
      endif
      do 101 j=1,ncomp
      wdtbcnd(iscomp+j,i+i0)=ysw(j)*sum
  101 continue
  109 continue
      if( aks(1).ne.undef ) then
        do 102 j=1,nrans
        wdtbcnd(israns+j,i+i0)=dble(aks(j))
  102   continue
      endif
  100 continue
      return
 9001 continue
      write(ifle,*) '### error : data error'
      write(ifle,*) 'sum of ys(i-j) = ',sum
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      if( yfix.gt.0 ) write(ifle,*) 'it must be > 0'
      if( yfix.lt.1 ) write(ifle,*) 'it must be = 0'
      write(ifle,*) 'i,j =',j0+1,j0+ncomp
      ierror=1
      end subroutine strdata
!
!< #6.8 calculate rotational matrix for periodic boundary >
!=================================================
      subroutine romtrx(n)
!=================================================
      integer,intent(in)  :: n
      integer :: i
      real*8  :: th,aa,ax,ay,az,cost,sint,cosu
      ax=dble(rot(1))
      ay=dble(rot(2))
      az=dble(rot(3))
      aa=sqrt(ax*ax+ay*ay+az*az)
      th=(rot(4)/45.d0)*atan(1.d0)
      if( min(aa,abs(th)).gt.0.d0 ) goto 801
      do 100 i=2,8
      wdtbcnd(i,n)=0.d0
  100 continue
      wdtbcnd(1,n)=1.d0
      wdtbcnd(5,n)=1.d0
      wdtbcnd(9,n)=1.d0
      return
  801 continue
      ax=ax/aa
      ay=ay/aa
      az=az/aa
      cost=cos(th)
      sint=sin(th)
      cosu=1.d0-cost
      wdtbcnd(1,n)=ax*ax*cosu+   cost
      wdtbcnd(2,n)=ax*ay*cosu-az*sint
      wdtbcnd(3,n)=ax*az*cosu+ay*sint
      wdtbcnd(4,n)=ay*ax*cosu+az*sint
      wdtbcnd(5,n)=ay*ay*cosu+   cost
      wdtbcnd(6,n)=ay*az*cosu-ax*sint
      wdtbcnd(7,n)=az*ax*cosu-ay*sint
      wdtbcnd(8,n)=az*ay*cosu+ax*sint
      wdtbcnd(9,n)=az*az*cosu+   cost
      end subroutine romtrx
!=======================================================
      subroutine surfreac(nphase,nsite,msite,nbulk,mbulk,
     &                    site_name,site_dens,
     &                    bulk_name,lc,KSUF
     &                )
!=======================================================
      integer,intent(in)       :: msite,mbulk,lc
      integer,intent(inout)    :: nsite,nbulk,nphase,KSUF
      character(lc),intent(in) :: bulk_name(mbulk),site_name(msite)
      real(8),intent(inout)    :: site_dens(msite)
!=======================================================
!
        nsite=0
        do nst=1,msite
        IF(site_name(nst)/=' ') then
          nsite=nsite+1
          nphase=nphase+1
          if(site_dens(nst)==undef) then
            write(ifle,'(a,I4)') 
     &      'ERR* SITE density NOT defined at BC no=',no
            stop 'ERR: SITE density [site_dens] must be defined'
          endif
        endif
        enddo
        nbulk=0
        do nbk=1,mbulk
        if(bulk_name(nbk)/=' ') then
          nbulk=nbulk+1
          nphase=nphase+1
        endif
        enddo
        KSUF=1
      end subroutine surfreac
!
      end subroutine inputdata
!
      end module module_boundary


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_cgsolver
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(21),parameter,private :: modnam='(module_cgsolver)'
!
!< module for cg solver >
!
      real*8  :: aepscg=1.d-12
      real*8  :: repscg=1.d-12
      integer :: itercg=1000
      real*8  :: aepsbcg=1.d-12
      real*8  :: repsbcg=1.d-12,aepscg_p,repscg_p,
     &                    aepsbcg_MHD=1.d-12,repsbcg_MHD=1.d-12
      integer :: iterbcg=1000,itercg_p=100,iterbcg_MHD=100
      integer :: nostop=0
!
! aepscg  : absolute error of tolerance for iccg
! repscg  : relative error of tolerance for iccg
! itercg  : maximul iteration counts for iccg
! aepsbcg : absolute error of tolerance for bicgstab
! repsbcg : relative error of tolerance for bicgstab
! iterbcg : maximul iteration counts for bicgstab
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! [namelist]
      namelist /cgsolver/ aepscg,repscg,itercg,
     &                    aepsbcg,repsbcg,iterbcg,nostop,
     &                    aepscg_p,repscg_p,itercg_p,
     &                    iterbcg_MHD,aepsbcg_MHD,repsbcg_MHD
!
! [local entities]
      integer :: iset=0
      integer :: ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      rewind ifli
      read(ifli,cgsolver,iostat=ios)
      if( ios.lt.0 ) return
!      write(ifll,cgsolver)
      call nml_errmsg0(ifle,ios,'cgsolver',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if( max(aepscg,repscg).le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'aepscg,repscg = ',aepscg,repscg
        write(ifle,*) 'max. of them must be > 0'
        goto 9999
      endif
      if( itercg.lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'itercg =',itercg
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if( max(aepsbcg,repsbcg).le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'aepsbcg,repsbcg = ',aepsbcg,repsbcg
        write(ifle,*) 'max. of them must be > 0'
        goto 9999
      endif
      if( iterbcg.lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'iterbcg =',iterbcg
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_cgsolver

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_chemcntl
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(21),parameter,private :: modnam='(module_chemcntl)'
!
!
!< module for control data for chemical reaction >
!
      integer :: igT_iter=0,
     &                 Zeldovich_P=0,
     &                 Fuel_NO_P=0,
     &                 prompt_NO_P=0
      integer :: itintg=1000000
      integer :: msizfc=10
      real*8  :: divint=1.d-2
      real*8  :: abserr=1.d-5
      real*8  :: relerr=1.d-5
!
! itintg : maxmul no. of integration steps
! msizfc : factor of memory size
! divint : ratio of initial interval to total interval
! abserr : absolute error of tolerance
!        : for integration of ordinary differential equation
! relerr : relative error of tolerance
!        : for integration of ordinary differential equation
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! [namelist]
      namelist /chemcntl/ igT_iter,
     &                    Zeldovich_P,Fuel_NO_P,prompt_NO_P,
     &                    itintg,msizfc,divint,abserr,relerr
!
! [local entities]
      integer :: iset=0
      integer :: ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      rewind ifli
      read(ifli,chemcntl,iostat=ios)
      if( ios.lt.0 ) return
!      write(ifll,chemcntl)
      call nml_errmsg0(ifle,ios,'chemcntl',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if( itintg.lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'itintg =',itintg
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if( msizfc.lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'msizfc =',msizfc
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if( divint.le.0.d0 .or. divint.gt.1.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'divint = ',divint
        if( divint.le.0.d0 ) write(ifle,*) 'it must be > 0'
        if( divint.gt.1.d0 ) write(ifle,*) 'it must be <= 1'
        goto 9999
      endif
!
      if( abserr.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'abserr = ',abserr
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
      if( relerr.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'relerr = ',relerr
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_chemcntl

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_debug
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(18),parameter,private :: modnam='(module_debug)'
!
!< module for debugging >
!
      logical :: check_prep=.true.
!
      integer,private :: i
      integer,parameter,private :: mdebug=100
      integer :: idebug(mdebug)=(/(0,i=1,mdebug)/)
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
      integer :: flag(mdebug)
      namelist /debug/ flag
      integer :: iset=0
      integer :: i,ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      flag=0
      rewind ifli
      read (ifli,debug,iostat=ios)
      if( ios.lt.0 ) return
!      write(ifll,debug)
      call nml_errmsg0(ifle,ios,'debug',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      do 100 i=1,mdebug
      if( flag(i).gt.0 .and. flag(i).le.mdebug ) idebug(flag(i))=1
  100 continue
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_debug
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_deltat
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(19),parameter,private :: modnam='(module_deltat)'
!!!      include'interfacemmm'
!
!< module for data to estimate time increment >
!
      integer,parameter :: const=1, auto=2
      integer :: ideltt=const
      real*8  :: dtmax=1.d0, dtsafe=0.8d0
!
! const  : constant time increment
! auto   : auto time increment
! ideltt : flag for time increment
! dtmax  : maximum delta-t
! dtsafe : safety factor for delta-t
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in):: cntlnam
      integer,intent(out) :: ierror
!
! [namelist]
      character(10) :: option
      real*8        :: dt,safe
      namelist /deltat/ option,dt,safe
!
! [local entities]
      real*8,parameter :: undef=-huge(1.d0)
      character(8),parameter :: dtlst(2)=(/'constant','auto    ' /)
      integer :: iset=0
      integer :: ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      option=' '
      dt    =undef
      safe  =dtsafe
!
      rewind ifli
      read(ifli,deltat,iostat=ios)
!      write(ifll,deltat)
      call nml_errmsg0(ifle,ios,'deltat',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if( dt.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : dt'
        goto 9999
      elseif( dt.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'dt = ',dt
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if( safe.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'safe = ',safe
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      call nml_listno(2,dtlst,option,ideltt)
      call nml_chkchr0(ifle,'option',option,ideltt,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      dtmax =dt
      dtsafe=safe
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_deltat

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_dimnsn
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(19),parameter,private :: modnam='(module_dimnsn)'

!
!< module for dimension of computaional domain >
!
      logical :: xredc=.false., yredc=.false., zredc=.false.
!
! xredc : =.true.; x-direction is reduced
! yredc : =.true.; y-direction is reduced
! zredc : =.true.; z-direction is reduced
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! [namelist]
      logical :: x,y,z
      namelist /dimension/ x,y,z
!
! [local entities]
      integer :: iset=0
      integer :: ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      x=xredc
      y=yredc
      z=zredc
!
      rewind ifli
      read(ifli,dimension,iostat=ios)
      if( ios.lt.0 ) return
!      write(ifll,dimension)
      call nml_errmsg0(ifle,ios,'dimension',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      xredc=x
      yredc=y
      zredc=z
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_dimnsn

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_flags
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(18),parameter,private :: modnam='(module_flags)'
!
!
!< module for control flags >
!
!
!-< 1. Option to calculate convection term >-
!
!      integer,parameter :: up1st=1, up2nd=2,cnt2nd=3,up3rd=4
!      integer,parameter :: nolim=0, venkat=1
!      integer :: icnvvv=up1st,  icnvty=up1st,  icnvke=up1st
!      integer :: lmtrvv=venkat, lmtrty=venkat, lmtrke=venkat
!      real*8  :: rnck=0.3d0,bldfct=0.D0
!
! up1st  : 1st order upwind
! up2nd  : higher order
! nolim  : no limiter
! venkat : Venkatakrishnan limiter
! lmtrvv : flag for limiter of momentum equation
! lmtrty : flag for limiter of temp. & conc. equation
! lmtrke : flag for limiter of k-e equation
! icnvvv : flag for convection term of momentum equation
! icnvty : flag for convection term of temp. & conc. equation
! icnvke : flag for convection term of k-e equation
! rnck   : parameter for near-constant region in Vekatakrishnan limiter
!
!
!-< 2. Option to calculate viscous term >-
!
!      integer,parameter :: nocal=0, cal=1
!      integer :: ivscvv=cal, ivscty=cal, ivscke=cal
!
! nocal  : not calculate
!   cal  : calculate
! ivscvv : flag for viscous term of momentum equation
! ivscty : flag for diffusion term of temp. & conc. equation
! ivscke : flag for diffusion term of k-e equation
!
!
!-< 3. Option to integrate conv. & diff. equation >-
!
      integer,parameter :: eulere=0, euleri=1,Adams_Bashforth=2,
     &                     Adams_Moulton=3,Crank_Nicolson=4,
     &                     Runge_Kutta=5
      integer :: intgvv=euleri, intgty=euleri, intgke=euleri
!
! eulere : Euler explicit
! euleri : Euler implicit
! intgvv : flag for momentum equation
! intgty : flag for temp. & conc. equation
! intgke : flag for k-e equation
! intg*  : flag for time integration of conv. & diff. equation
!
!
!-< 4. Option to interpolate velocity >-
!
      integer,parameter :: face=0, center=1
      integer :: icalrv=face
!
! face   : interpolate velocity at cell face
! center : interpolate velocity at cell center
! icalrv : flag to interpolate velocity
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,lrans,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in):: cntlnam
      logical,intent(in)  :: lrans
      integer,intent(out) :: ierror
!
! [namelist]
!      logical       :: visc,diff_ty,diff_ke
!      character(10) :: cnv_vv,cnv_ty,cnv_ke
!      character(10) :: limitr_vv,limitr_ty,limitr_ke
      character(10) :: integ_vv,integ_ty,integ_ke,integ_MHD
      character(10) :: intrv
      namelist /flags/ 
     &                 integ_vv,integ_ty,integ_ke,
     &                 intrv,integ_MHD
!
! [local entities]
!      character( 3),parameter :: cnvlst(4)=(/'1st','2nd','c2d','3rd'/)
!      character( 8),parameter :: lmtlst(2)=(/'no      ','Venkatak'/)
      character( 8),parameter :: intlst(6)=(/'explicit','implicit',
     & 'adamsbfh','moultons','cranknic','rungekta'/)
      character( 6),parameter :: invlst(2)=(/'face  ','center'/)
      integer :: iset=0
      integer :: ios,ierr1,ierr2,ierr3,ierr4,ierr5,ierr6
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
!      visc     =lxi(ivscvv)
!      diff_ty  =lxi(ivscty)
!      diff_ke  =lxi(ivscke)
!      cnv_vv   =' '
!      cnv_ty   =' '
!      cnv_ke   =' '
!      limitr_vv=lmtlst(lmtrvv+1)
!      limitr_ty=lmtlst(lmtrty+1)
!      limitr_ke=lmtlst(lmtrke+1)
      integ_vv =' '
      integ_ty =' '
      integ_ke =' '
      intrv    =invlst(max(1,min(2,icalrv+1)))
!
      ierr1 = 0
      ierr2 = 0
      ierr3 = 0
      ierr4 = 0
      ierr5 = 0
      ierr6 = 0
!
      rewind ifli
      read(ifli,flags,iostat=ios)
!      write(ifll,flags)
      call nml_errmsg0(ifle,ios,'flags',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
!      call nml_listno(4,cnvlst,cnv_vv   ,icnvvv)
!      call nml_listno(4,cnvlst,cnv_ty   ,icnvty)
!      call nml_listno(2,lmtlst,limitr_vv,lmtrvv)
!      call nml_listno(2,lmtlst,limitr_ty,lmtrty)
      call nml_listno(6,intlst,integ_vv,intgvv)
      call nml_listno(6,intlst,integ_ty,intgty)
!      call nml_chkchr0(ifle,'cnv_vv'   ,cnv_vv   ,icnvvv,ierr1)
!      call nml_chkchr0(ifle,'cnv_ty'   ,cnv_ty   ,icnvty,ierr2)
!      call nml_chkchr0(ifle,'limitr_vv',limitr_vv,lmtrvv,ierr3)
!      call nml_chkchr0(ifle,'limitr_ty',limitr_ty,lmtrty,ierr4)
      call nml_chkchr0(ifle,'integ_vv' ,integ_vv ,intgvv,ierr5)
      call nml_chkchr0(ifle,'integ_ty' ,integ_ty ,intgty,ierr6)
      if( ierr1.ne.0 .or. ierr2.ne.0 .or.
     &    ierr3.ne.0 .or. ierr4.ne.0 .or.
     &    ierr5.ne.0 .or. ierr6.ne.0 ) goto 9999
      if( lrans ) then
        call nml_listno(6,intlst,integ_ke ,intgke)
        call nml_chkchr0(ifle,'integ_ke' ,integ_ke ,intgke,ierr3)
        if( ierr1.ne.0 .or. ierr2.ne.0 .or. ierr3.ne.0 ) goto 9999
      else
        intgke=1
      endif
      call nml_listno(2,invlst,intrv   ,icalrv)
      call nml_chkchr0(ifle,'intrv',intrv,icalrv,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
!      ivscvv=ixl(visc)
!      ivscty=ixl(diff_ty)
!      ivscke=ixl(diff_ke)
!
!      lmtrvv=lmtrvv-1
!      lmtrty=lmtrty-1
!      lmtrke=lmtrke-1
!
      intgvv=intgvv-1
      intgty=intgty-1
      intgke=intgke-1
!
      if(intgty.eq.Adams_Bashforth) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ty] is switch to be [implicit] '
        intgty=euleri
      endif
      if(intgke.eq.Adams_Bashforth) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ke] is switch to be [implicit] '
        if( lrans ) then
          intgke=euleri
        else
          intgke=eulere
        endif
      endif
      if(intgty.eq.Adams_Moulton) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ty] is switch to be [implicit] '
        intgty=euleri
      endif
      if(intgke.eq.Adams_Moulton) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ke] is switch to be [implicit] '
        if( lrans ) then
          intgke=euleri
        else
          intgke=eulere
        endif        
      endif
      if(intgty.eq.Crank_Nicolson) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ty] is switch to be [implicit] '
        intgty=euleri
      endif
      if(intgke.eq.Crank_Nicolson) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ke] is switch to be [implicit] '
        if( lrans ) then
          intgke=euleri
        else
          intgke=eulere
        endif        
      endif
      if(intgty.eq.Runge_Kutta) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ty] is switch to be [implicit] '
        intgty=euleri
      endif
      if(intgke.eq.Runge_Kutta) then
        write(ifle,2000) 
        write(ifle,*) '[integ_ke] is switch to be [implicit] '
        if( lrans ) then
          intgke=euleri
        else
          intgke=eulere
        endif        
      endif
 2000 FORMAT(/,2X,4X,3('?'),1X,'WARNING',1X,3('?'),/,6X)
!
      icalrv=icalrv-1
!
! --- Vertex Center uses 'icalrv=face'
      icalrv=face
!
!      if( .not.lrans ) ivscke=nocal
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
!
      end subroutine inputdata
!
      end module module_flags
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_material
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(18),parameter,private :: modnam='(module_material)'
!
! --- < 1. parameter for fluid domain >
!
! nflud : no. of fluid domains
! lclsd : flag for closed domain
!       : <1; open domain
!       : >0; closed domain
!
      integer,parameter :: suther=1, const=2,simpli=3,MK=4,chung=5
      real*8  :: sthrmu_i=18.d-6, sthrt_i=293.15d0, sthrc_i=117.d0
      real*8  :: prlmnr_i=0.72d0, sclmnr_i=1.d0
!
! suther : Sutherland's formula
! const  : constant
! simpli : simplified transport model
! isthr  : flag for viscosity
! sthrmu : viscosity at temp.=sthrt [Pa.s]
! sthrt  : reference temperature [K]
! sthrc  : Sutherland's constant [K]
! prlmnr : Prandtl no. for laminar flow [-]
! sclmnr : Schmidt no. for laminar flow [-]
!
!
      integer,parameter :: up1st=1, up2nd=2,cnt2nd=3,up3rd=4,
     &                     usi2nd=5,engcsv=6,c3bld=7
      integer,parameter :: nolim=0, venkat=1,slope=2
      real*8  :: rnck=0.3d0
!
! up1st  : 1st order upwind
! up2nd  : higher order
! nolim  : no limiter
! venkat : Venkatakrishnan limiter
! lmtrvv : flag for limiter of momentum equation
! lmtrty : flag for limiter of temp. & conc. equation
! lmtrke : flag for limiter of k-e equation
! icnvvv : flag for convection term of momentum equation
! icnvty : flag for convection term of temp. & conc. equation
! icnvke : flag for convection term of k-e equation
! rnck   : parameter for near-constant region in Vekatakrishnan limiter
!
! ---  Option to calculate viscous term >-
!
      integer,parameter :: nocal=0, cal=1
!
! nocal  : not calculate
!   cal  : calculate
! ivscvv : flag for viscous term of momentum equation
! ivscty : flag for diffusion term of temp. & conc. equation
! ivscke : flag for diffusion term of k-e equation
!
      integer :: nflud=1
      integer,save,allocatable :: lclsd(:)
      integer,save,allocatable :: nofld(:)
      integer,save,allocatable :: isthr(:)
      integer,save,allocatable :: PFstart(:)
      real*8 ,save,allocatable :: sthrmu(:),sthrt(:),sthrc(:)
      real*8 ,save,allocatable :: prlmnr(:),sclmnr(:)
      real*8 ,save,allocatable :: r_prlmnr(:),r_sclmnr(:)
      logical,save,allocatable :: cvdist(:)
      integer,save,allocatable :: icnvvv(:)
      integer,save,allocatable :: icnvty(:)
      integer,save,allocatable :: icnvke(:)
      integer,save,allocatable :: lmtrvv(:)
      integer,save,allocatable :: lmtrty(:)
      integer,save,allocatable :: lmtrke(:)
      real*8 ,save,allocatable :: bldf(:)
      integer,save,allocatable :: ivscvv(:)
      integer,save,allocatable :: ivscty(:)
      integer,save,allocatable :: ivscke(:)
! --- Sliding & rotation machine
      real*8 ,save,allocatable :: rotati(:)
      integer,save,allocatable :: ishaft(:)
      real*8 ,save,allocatable :: end(:,:),begin(:,:)
      integer,save,allocatable :: nsplit(:)
      integer,save             :: ical_sld=0
      integer,save,allocatable :: MAT_ALL(:)
      integer :: KMAT_S=0
!
      integer,save,allocatable :: DONOR_M(:)
! --- < 2. parameter for solid domain >
!
! nsold  : no. of solids [-]
! nosld  : solid no. [-]
! rsld   : density [kg/m^3]
! cpsld  : specific heat [J/kg/K]
! rmdsld : thermal conductivity [W/m/K]
! rotati : rpm [revolution per minute] with right-hand-law
! ishaft : axis of revolution,default=>[z]
!
      integer :: nsold=0
      integer,save,allocatable :: nosld(:)
      real*8 ,save,allocatable :: rsld(:)
      real*8 ,save,allocatable :: cpsld(:)
      real*8 ,save,allocatable :: rmdsld(:)
      integer,save,allocatable :: PSstart(:)
!
      real*8,save,allocatable  :: radmat(:,:)
!1: gabsorb, 2: pabsorb, 3: pscatter, 4:scabeta, 5: ptcdiamter
      INTEGER,save,allocatable :: radfludflag(:,:)
!1: gastype, 2: ptctype, 3: ptccal, 4: fgpflag
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. store data >-----------------------------------------------------
!=================================================
      subroutine strdata
     &  (ifle,nfludx,MATC_FLD,MATC_NO,jclos,ihpc,ierr)
!=================================================
      implicit none
      character(9),parameter :: subnam='(strdata)'
      integer,intent(in)  :: ifle,nfludx,ihpc
      integer,intent(in)  :: jclos(0:100),MATC_FLD(100),MATC_NO(200)
      integer,intent(out) :: ierr
      integer :: i,iset=0,IIMAT,IMAT
!
      if(ihpc.eq.0) then
        call modutl_chkset('=',1,iset,modnam//subnam)
      endif
      ierr=0
!
      allocate(lclsd(200),stat=ierr)
!
      if(ierr.ne.0) then
        write(ifle,*) '### error : allocation failed'
        goto 9999
      endif
!
      lclsd=0
      do 100 i=1,nfludx
      IIMAT=MATC_FLD(i)
      IMAT=MATC_NO(IIMAT)
      if(IMAT.lt.0) then
        write(ifle,*) '### error : FLUID has Negative IMAT No.',IMAT
        goto 9999
      endif
      if( jclos(IMAT).gt.0 ) then
! --- at last one outlet BC exists
        lclsd(IIMAT)=0
      else
! --- closed domain
        lclsd(IIMAT)=1
      endif
  100 continue
      return
 9999 continue
      write(ifle,*) modnam,subnam
      ierr=1
      end subroutine strdata
!
!============================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,
     &                     lsld,lrans,lovstbc,ierror)
!============================================================
!
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      integer,intent(in)    :: ifli,ifll,ifle
      character(*),intent(in)  :: cntlnam
      logical,intent(in)    :: lrans,lovstbc
      logical,intent(inout) :: lsld
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer                :: kwall=0
      character(3),parameter :: cond(2)=(/'yes','no '/)
!
! --- [namelist-fluid]
!
      integer       :: IMAT_U,Poisson_start
      character(10) :: muopt,wall_distance
      real*8        :: mu,t0,c,bldfct,mu_turb,mu_turb2
      real*8        :: Prandtl,Schmidt
      logical       :: visc,diff_ty,diff_ke,diff_MHD
      real*8        :: unrelax_P=1.0,unrelax_T=1.0,unrelax_Ys=1.0,
     &                  unrelax_V=1.d0,unrelax_RANS=1.d0,
     &                 unrelax_RC
      character(10) :: cnv_vv,cnv_ty,cnv_ke,cnv_MHD
      character(10) :: limitr_vv,limitr_ty,limitr_ke
      real*8        :: Prandtl2,Schmidt2,mu2,t02,c2,ts,hgs,hls
      real*8        :: rpm,end_x,end_y,end_z,begin_x,begin_y,begin_z
      integer       :: start_rot_step,rotup_step
      integer       :: N_SPLIT,unrelax_MHD
      real*8,parameter :: SML=1.d-15
      real(8)        :: angle_init=0.d0
      integer,parameter :: N_pors=80
      real(8)        :: C0_porous,C1_porous,perm_porous,C2_porous,
     &                  porosity
      character(LEN=5) :: porous_law
      integer        :: bg_gas=1
      real(8)        :: gabsorb=0.0,pabsorb=0.0,
     &			pscatter=0.0,scabeta=0.0,
     &			ptcdia=0.0,ptcpara(2),MUSCL_kapa
      INTEGER        :: scatype=1,gastype=1,ptccal=0,fgpflag=1
      INTEGER        :: Gauge_MHD_A,Gauge_MHD_Je,Gauge_MHD_J0
      INTEGER        :: IMAT_DONOR=0
!      real(8)        :: OMEGA_MHD
!
      namelist /fluid/ IMAT_U,muopt,
     &                 wall_distance,Poisson_start,
     &                 cnv_vv,cnv_ty,cnv_ke,cnv_MHD,
     &                 limitr_vv,limitr_ty,limitr_ke,
     &                 rnck,bldfct,MUSCL_kapa,
     &                 visc,diff_ty,diff_ke,diff_MHD,
     &                 unrelax_P,unrelax_T,unrelax_Ys,unrelax_V,
     &                 unrelax_RANS,unrelax_RC,
     &                 Prandtl ,Schmidt ,mu ,t0 ,c,mu_turb,mu_turb2,
     &                 Prandtl2,Schmidt2,mu2,t02,c2,ts,hgs,hls,
     &                 rpm,end_x,end_y,end_z,begin_x,begin_y,begin_z,
     &                 start_rot_step,rotup_step,angle_init,
     &                 porous_law,C0_porous,C1_porous,C2_porous,
     &                 perm_porous,porosity,bg_gas,
     &		       gabsorb,pabsorb,pscatter,scabeta,ptcdia,
     &		       gastype,scatype,ptccal,fgpflag,ptcpara,unrelax_MHD,
     &                 Gauge_MHD_A,Gauge_MHD_Je,Gauge_MHD_J0,
     &                 IMAT_DONOR
     &                 
!
!----------------------------------------------------------------------
! rpm : revolution of IMAT_U (revolution per minute)
! unrelax_P : relaxzation for pressure
! unrelax_T : relaxzation for enthalpy
! Poisson_start : 
!----------------------------------------------------------------------
! --- [namelist-solid]
!
      real*8  :: rho,cp,cndct
      namelist /solid/ IMAT_U,rho,cp,cndct,Poisson_start,unrelax_T,
     &                 diff_MHD,cnv_MHD,unrelax_MHD,
     &                 Gauge_MHD_A,Gauge_MHD_Je,Gauge_MHD_J0
!
! --- [local entities]
!
      character( 3),parameter :: cnvlst(8)=
     &         (/'1st','2nd','c2d','3rd','usi','eng','c3d','MSL'/)
      character( 8),parameter :: lmtlst(3)=(/'no      ','Venkatak',
     &                                      'slope   '/)
      integer,parameter   :: iundef=-huge(1)
      real*8 ,parameter   :: undef=-huge(1.d0)
      character(10) :: mulst(5)=
     & (/'Sutherland','constant  ','Simplified','MK        '
     &     ,'Chung     '
     &     /)

! erase 'parameter' attribute for handling mulst
!     because of avoid warning on Windows -- by onishi
!      character(10),parameter :: mulst(3)=
!     & (/'Sutherland','constant  ','Simplified'/)

      integer :: i,iflw,itrb,ios,ierr1,ierr2,ierr3,ierr4
      integer :: iset=0
      integer :: n,nd
      logical,external :: nml_comp_eq, nml_comp_gt
      real*8 :: rdum
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
! --- fluid material
!
!--< 1.1 count up sizes >--
!
      nd=0
      rewind ifli
!
      do
        read(ifli,fluid,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'fluid',ierr1)
        nd = nd+1
        if( ierr1/=0 ) then
          write(ifle,*) 'sequence no. of the [&fluid] =',nd
          goto 9999
        endif
      enddo
      if(nd.lt.1) then
        write(ifle,*) '### warning: no fluid specified'
!MHD_err        goto 9999
      endif
      nflud=nd
!
!--< 1.2 allocate arrays >--
!
      allocate( nofld (nflud),
     &          isthr (nflud),
     &          sthrmu(nflud),
     &          sthrt (nflud),
     &          sthrc (nflud),
     &          prlmnr(nflud),
     &          sclmnr(nflud),
     &          r_prlmnr(nflud),
     &          r_sclmnr(nflud),
     &          cvdist(nflud),
     &          PFstart(nflud),
     &          icnvvv(nflud),
     &          icnvty(nflud),
     &          icnvke(nflud),
     &          lmtrvv(nflud),
     &          lmtrty(nflud),
     &          lmtrke(nflud),
     &          bldf(nflud),
     &          ivscvv(nflud),
     &          ivscty(nflud),
     &          ivscke(nflud),
     &          rotati(nflud),
     &          ishaft(nflud),
     &          DONOR_M(nflud),
     &          nsplit(nflud),
     &          end(3,nflud),begin(3,nflud),
     &          radmat(5,nflud),radfludflag(4,nflud),!jiang
     &          stat=ierr1 )
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error-1 : allocation failed'
        goto 9999
      endif
      nofld=0
!
!-< 1.3 Input namelist >-
!
      sthrmu(:)=sthrmu_i
      sthrt(:) =sthrt_i
      sthrc(:) =sthrc_i
      prlmnr(:)=prlmnr_i
      sclmnr(:)=sclmnr_i
      r_prlmnr(:)=1.d0/prlmnr_i
      r_sclmnr(:)=1.d0/sclmnr_i
      cvdist(:)=.false.
      PFstart(:)=0
!
      icnvvv(:)=up1st
      icnvty(:)=up1st
      icnvke(:)=up1st
      lmtrvv(:)=nolim
      lmtrty(:)=nolim
      lmtrke(:)=nolim
      bldf(:)=0.D0
      ivscvv(:)=cal
      ivscty(:)=cal
      ivscke(:)=cal
      rotati(:)=0.D0
      ishaft(:)=0
      end(:,:)=0.D0
      begin(:,:)=0.D0
      DONOR_M(:)=0
!
      nd=0
      rewind ifli
      do nd=1,nflud
!
        IMAT_U =iundef
        muopt  =' '
        mu     =undef
        t0     =undef
        c      =undef
        Prandtl=prlmnr_i
        Schmidt=sclmnr_i
        wall_distance='no '
        rpm=undef
        end_x=0.d0
        end_y=0.d0
        end_z=0.d0
        begin_x=0.d0
        begin_y=0.d0
        begin_z=0.d0
        N_SPLIT=36
        IMAT_DONOR=0
!
        Poisson_start=iundef
        bldfct=0.d0
        visc     =.true.
        diff_ty  =.true.
        diff_ke  =.true.
        diff_MHD =.true.
        cnv_vv   =' '
        cnv_ty   =' '
        cnv_ke   =' '
        limitr_vv=lmtlst(slope+1)
        limitr_ty=lmtlst(slope+1)
        limitr_ke=lmtlst(slope+1)
        porous_law='no    '
!
        read(ifli,fluid,iostat=ios)
        if( ios<0 ) exit
        call nml_errmsg0(ifle,ios,'fluid',ierr1)
        if( ierr1/=0 ) goto 9999
!
        if( IMAT_U.eq.iundef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : IMAT_U in [&fluid]'
          goto 9999
        elseif( IMAT_U.le.0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'IMAT_U =',IMAT_U
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
!
        if(porous_law/='no    ') then
          if(IMAT_U<=N_pors) then
            write(ifle,'(1X,a)') 
     &       'ERR: Porous madia: must bc [IMAT_U>80]',IMAT_U
            stop 
          endif
        else
          if(IMAT_U>N_pors) then
             write(ifle,'(1X,a)') 
     &       'ERR: NON Porous madia(free fluid): must bc [IMAT_U<80]'
            stop 
          endif
        endif
!
        do n=1,nd-1
          if( IMAT_U.eq.nofld(n) ) then
            write(ifle,*) '### error : data error'
            write(ifle,*) 'IMAT_U =',IMAT_U
            write(ifle,*) 'this IMAT_U already has been specified'
            write(ifle,*) 'sequence no. of the line =',n
            goto 9999
          endif
        enddo
!        call nml_listno(4,mulst,muopt,isthr(nd))
        call nml_listno(5,mulst,muopt,isthr(nd))
        call nml_chkchr0(ifle,'muopt',muopt,isthr(nd),ierr1)
        if(ierr1.ne.0) goto 9999
!
        select case(isthr(nd))                ! transport model
          case(1)
            call read_Sutherland
          case(2)
            call read_constant
          case(3)
            call read_Simplified
        end select
!
        if(bldfct>1.d0.or.bldfct<0.d0) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'Blend Factor [bldfct] must between 1 and 0'
          goto 9999
        endif
!
        if(rpm/=undef) then
        rdum=(end_x-begin_x)*(end_x-begin_x)
     &      +(end_y-begin_y)*(end_y-begin_y)
     &      +(end_z-begin_z)*(end_z-begin_z)
        if(rdum.lt.SML) then
          write(ifle,'(a)') 'ERR: You must define revolution shaft '
          write(ifll,'(3a)') 
     &   'MSG: Vector of revolution shaft is defined as ',
     &                  'with : end_x,end_y,end_z',
     &                  'and  : begin_x,begin_y,begin_z'
          stop 'ERR: NOT define cood_end'
        endif
        if(N_SPLIT.lt.0)then
          write(ifle,'(a)') 'ERR: N_SPLIT must be > 0'
          stop 'ERR: NOT define [N_SPLIT]'
        endif
        if(abs(rpm).lt.SML) then
          ishaft(nd)=0
          rotati(nd)=0.d0
          lsld=.false.
          nsplit(nd)=0
        else
          ishaft(nd)=1
          rotati(nd)=rpm/60.d0*2.d0*3.1415926
          lsld=.true.
          end(1,nd)=end_x
          end(2,nd)=end_y
          end(3,nd)=end_z
          begin(1,nd)=begin_x
          begin(2,nd)=begin_y
          begin(3,nd)=begin_z
        endif
        endif
!
        call nml_listno(2,cond,wall_distance,kwall)
        call nml_chkchr0(ifle,'wall_distance',wall_distance,kwall,ierr1)
        if( ierr1.ne.0 ) goto 9999
!
        if(Poisson_start==iundef) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : Poisson_start in [&fluid]'
          write(ifle,*) 'Poisson_start is starting iter for Poisson'
          goto 9999
        elseif( Poisson_start<=-1 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'Poisson_start =',Poisson_start
          write(ifle,*) 'it must be >= 0'
          goto 9999
        endif
!
!        if(lovstbc) then
!          if(IMAT_DONOR==0) then
!            write(*,*) 'MSG: over set BC used'
!            write(*,*) 'ERR: NOT defined , IMAT_DONOR=',IMAT_DONOR
!            write(*,*) 
!     &        'MSG: Defining Donor Material no by  IMAT_DONOR in &fluid'
!          else
            DONOR_M(nd)=IMAT_DONOR
!          endif
!        endif
!
        call nml_listno(8,cnvlst,cnv_vv   ,icnvvv(nd))
        call nml_listno(8,cnvlst,cnv_ty   ,icnvty(nd))
        call nml_listno(3,lmtlst,limitr_vv,lmtrvv(nd))
        call nml_listno(3,lmtlst,limitr_ty,lmtrty(nd))
        call nml_chkchr0(ifle,'cnv_vv'   ,cnv_vv   ,icnvvv(nd),ierr1)
        call nml_chkchr0(ifle,'cnv_ty'   ,cnv_ty   ,icnvty(nd),ierr2)
        call nml_chkchr0(ifle,'limitr_vv',limitr_vv,lmtrvv(nd),ierr3)
        call nml_chkchr0(ifle,'limitr_ty',limitr_ty,lmtrty(nd),ierr4)
        if(ierr1.ne.0.or.ierr2.ne.0.or.
     &    ierr3.ne.0.or.ierr4.ne.0) goto 9999
!
        if(lrans) then
          call nml_listno(8,cnvlst,cnv_ke   ,icnvke(nd))
          call nml_listno(3,lmtlst,limitr_ke,lmtrke(nd))
          call nml_chkchr0(ifle,'cnv_ke'   ,cnv_ke   ,icnvke(nd),ierr1)
          call nml_chkchr0(ifle,'limitr_ke',limitr_ke,lmtrke(nd),ierr2)
          if( ierr1.ne.0 .or. ierr2.ne.0 .or. ierr3.ne.0 ) goto 9999
        else
          icnvke(nd)=1
          lmtrke(nd)=lmtrke(nd)+1
        endif
        ivscvv(nd)=ixl(visc)
        ivscty(nd)=ixl(diff_ty)
        ivscke(nd)=ixl(diff_ke)
!
        lmtrvv(nd)=lmtrvv(nd)-1
        lmtrty(nd)=lmtrty(nd)-1
        lmtrke(nd)=lmtrke(nd)-1
        if( .not.lrans ) ivscke(nd)=nocal
!
        nofld(nd)=IMAT_U
        sthrmu(nd)=mu
        sthrt(nd) =t0
        sthrc(nd) =c
        prlmnr(nd)=Prandtl
        sclmnr(nd)=Schmidt
        if(kwall.eq.1) then
          cvdist(nd)=.true.
        elseif(kwall.eq.2) then
          cvdist(nd)=.false.
        else
          write(ifle,*) "ERR: wall_distance='yes' or 'no' in ",cntlnam
          goto 9999
        endif
        PFstart(nd)=Poisson_start
!
      if(gabsorb.lt.0.0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'gabsorb =',gabsorb
        write(ifle,*) 'it should be >= 0'
        goto 9999
      end if
      if(pabsorb.lt.0.0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'pabsorb =',pabsorb
        write(ifle,*) 'it should be >= 0'
        goto 9999
      end if
      if(pscatter.lt.0.0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'pscatter =',pscatter
        write(ifle,*) 'it should be >= 0'
        goto 9999
      end if
      if(scatype.lt.1.OR.scatype.GT.5) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'scatype =',scatype
        write(ifle,*) 'it should be in (1 ~ 5)'
        goto 9999
      end if
      if(scatype.EQ.1.and.pscatter.GT.0.0) then
        if(abs(scabeta).gt.1d0) then
          write(ifle,*) '### error : data error'
	  write(ifle,*) 'scabeta =',scabeta
	  write(ifle,*) 'it should be in (-1 ~ 1)'
          goto 9999
        endif
      endif

      if(ptccal.EQ.1.and.ptcdia.LE.0.d0) then
	write(ifle,*) '### error : data error'
	write(ifle,*) 'ptcdia =',ptcdia
	write(ifle,*) 'particle diamater should be > 0'
	goto 9999
      endif

      if(gastype.ne.1.and.gastype.ne.2) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'gastype =',gastype
        write(ifle,*) 'it should be, 1: gray_gas, 2: real_gas'
        goto 9999
      endif

      radmat(1,nd)=gabsorb
      radmat(2,nd)=pabsorb
      radmat(3,nd)=pscatter
      radmat(4,nd)=scabeta
      radmat(5,nd)=ptcdia

      radfludflag(1,nd)=gastype
      radfludflag(2,nd)=scatype
      radfludflag(3,nd)=ptccal
      radfludflag(4,nd)=fgpflag

!
      enddo
!
      ical_sld=0
      do nd=1,nflud
      if(ishaft(nd)==1) then
        ical_sld=1
      endif
      enddo
!
	  
!
! --- solid material
!
      ierr1=0
      ios=0
      nd=0
      rewind ifli
  200 continue
      read(ifli,solid,iostat=ios)
      if( ios.lt.0 ) goto 201
      nd=nd+1
      call nml_errmsg0(ifle,ios,'solid',ierr1)
      if( ierr1.ne.0 ) then
        write(ifle,*) 'sequence no. of the [&solid] =',nd
        goto 9999
      endif
      goto 200
  201 continue
      if( nd.lt.1 ) goto 2001
      nsold=nd
!
!--< 2.1 allocate arrays >--
!
      allocate( nosld (nsold),
     &          rsld  (nsold),
     &          cpsld (nsold),
     &          rmdsld(nsold), stat=ierr1 )
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error-2 : allocation failed'
        goto 9999
      endif
      nosld =0
      rsld  =0.d0
      cpsld =0.d0
      rmdsld=0.d0
!
!-< 2.2 Input namelist >-
!
      nd=0
      rewind ifli
 2000 continue
!
      IMAT_U   =iundef
      rho  =undef
      cp   =undef
      cndct=undef
!
      read(ifli,solid,iostat=ios)
      if( ios.lt.0 ) goto 2001
      call nml_errmsg0(ifle,ios,'solid',ierr1)
      if( ierr1.ne.0 ) goto 9999
      nd=nd+1
!
      if( IMAT_U.eq.iundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : IMAT_U [&solid]'
        goto 9999
      elseif( IMAT_U.ge.0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'it must be < 0'
        goto 9999
      endif
      do 2100 n=1,nd-1
      if( IMAT_U.eq.nosld(n) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'this IMAT_U. already specified in previous line'
        write(ifle,*) 'sequence no. of the line =',n
        goto 9999
      endif
 2100 continue
!
      if( .not.nml_comp_eq("rho",rho,undef,nd,ifle,ierror) ) then
        if( nml_comp_gt("rho",rho,0.0d0,nd,ifle,ierror) )
     &    rsld(nd) = rho
      endif
      if( .not.nml_comp_eq("cp",cp,undef,nd,ifle,ierror) ) then
        if( nml_comp_gt("cp",cp,0.0d0,nd,ifle,ierror) )
     &    cpsld(nd) = cp
      endif
      if( .not.nml_comp_eq("cndct",cndct,undef,nd,ifle,ierror) ) then
        if( nml_comp_gt("cndct",cndct,0.0d0,nd,ifle,ierror) )
     &    rmdsld(nd) = cndct
      endif
!
      if( Poisson_start.eq.iundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : Poisson_start in [&solid]'
        write(ifle,*) 'Poisson_start is starting iter for Poisson'
        goto 9999
      elseif( Poisson_start.le.-1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'Poisson_start =',Poisson_start
        write(ifle,*) 'it must be >= 0'
        goto 9999
      endif
!
      nosld (nd)=IMAT_U
      rsld  (nd)=rho
      cpsld (nd)=cp
      rmdsld(nd)=cndct
!
      goto 2000
 2001 continue
!
!      call display
      allocate(MAT_ALL(nflud+nsold))
      KMAT_S=0
      do nd=1,nflud
        KMAT_S=KMAT_S+1
        MAT_ALL(KMAT_S)=nofld(nd)
      enddo
!
      do nd=1,nsold
        KMAT_S=KMAT_S+1
        MAT_ALL(KMAT_S)=nosld(nd)
      enddo
!
      if( ierror/=0 ) then      ! final error check
        write(ifle,*) "total error : ",ierror,modnam,subnam
      endif
      write(ifll,"(1x,80('-'))")
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      contains
!/////////////////////////////////////////////////////////////
!< #1.1 integer to logical > >>>
!============================================
      function lxi(i)
!============================================
      logical :: lxi
      integer,intent(in) :: i
                   lxi=.false.
      if( i.gt.0 ) lxi=.true.
      end function lxi
!
!< #1.2 logical to integer > >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!============================================
      function ixl(l)
!============================================
      integer :: ixl
      logical,intent(in) :: l
              ixl=0
      if( l ) ixl=1
      end function ixl
!
!============================================
      subroutine read_Sutherland                ! Sutherland
!============================================
      if( .not.nml_comp_eq("mu",mu,undef,nd,ifle,ierror) ) then
        if( nml_comp_gt("mu",mu,0.0d0,nd,ifle,ierror) )
     &    sthrmu(nd) = mu
      endif
      if( .not.nml_comp_eq("t0",t0,undef,nd,ifle,ierror) ) then
        if( nml_comp_gt("t0",t0,0.0d0,nd,ifle,ierror) )
     &    sthrt(nd) = t0
      endif
      if( .not.nml_comp_eq("c",c,undef,nd,ifle,ierror) )
     &  sthrc(nd) = c
      if( nml_comp_gt("Prandtl",Prandtl,0.0d0,nd,ifle,ierror) ) then
        prlmnr(nd) = Prandtl
        r_prlmnr(nd) = 1/prlmnr(nd)
      endif
      if( nml_comp_gt("Schmidt",Schmidt,0.0d0,nd,ifle,ierror) ) then
        sclmnr(nd) = Schmidt
        r_sclmnr(nd) = 1/sclmnr(nd)
      endif
      end subroutine read_Sutherland
!
!============================================
      subroutine read_constant                  ! constant
!============================================
      if( .not.nml_comp_eq("mu",mu,undef,nd,ifle,ierror) ) then
        if( nml_comp_gt("mu",mu,0.0d0,nd,ifle,ierror) )
     &    sthrmu(nd) = mu
      endif
      if( nml_comp_gt("Prandtl",Prandtl,0.0d0,nd,ifle,ierror) ) then
        prlmnr(nd) = Prandtl
        r_prlmnr(nd) = 1/prlmnr(nd)
      endif
      if( nml_comp_gt("Schmidt",Schmidt,0.0d0,nd,ifle,ierror) ) then
        sclmnr(nd) = Schmidt
        r_sclmnr(nd) = 1/sclmnr(nd)
      endif
      sthrt(nd) = sthrt_i
      sthrc(nd) = sthrc_i
      end subroutine read_constant
!
!============================================
      subroutine read_Simplified                ! simplified
!============================================
      prlmnr(nd) = 0.75                 ! always constant
      r_prlmnr(nd) = 1/prlmnr(nd)
      end subroutine read_Simplified
!
!============================================
      subroutine display
!============================================
      write(ifll,"(1x,80('-'))")
      do n=1,nflud
        write(ifll,"(1x,a,i3)")
     &      "Domain number              : ",n
        write(ifll,"(1x,2a)")
     &      "transport model            : ",mulst(isthr(n))
        if( isthr(n)==1 ) then          ! Sutherland
          call disp_Sutherland
        elseif( isthr(n)==2 ) then              ! constant
          call disp_constant
        elseif( isthr(n)==3 ) then              ! simplified
          call disp_Simplified
        endif
      enddo
      end subroutine display
!
!============================================
      subroutine disp_Sutherland
!============================================
      write(ifll,"(3x,a,e10.3,a)")
     &  "viscosity at 'T0' [Pa s] : ",sthrmu(n)," (mu0)"
      write(ifll,"(3x,a,f10.3,a)")
     &  "temperature [K]          : ",sthrt(n)," (T0)"
      write(ifll,"(3x,a,f10.3,a)")
     &  "constant value           : ",sthrc(n)," (C)"
      write(ifll,"(3x,a,f10.3,a)")
     &  "Schmidt number           : ",sclmnr(n)," (Sc)"
      write(ifll,"(3x,a,f10.3,a)")
     &  "Prandtl number           : ",prlmnr(n)," (Pr)"
      write(ifll,"(5x,2a)")
     &  "viscosity    : m = mu0*(T0+C)/(T+C)*(T/T0)^(3/2)",
     &  "  (T:temperature)"
      write(ifll,"(5x,2a)")
     &  "conductivity : l = Cp*m/Pr",
     &  "  (Cp:isopiestic specific heat)"
      write(ifll,"(5x,2a)")
     &  "diffusivity  : D = m/(r*Sc)",
     &  "  (r:mass density)"
      end subroutine disp_Sutherland
!
!============================================
      subroutine disp_constant
!============================================
      write(ifll,"(3x,a,e10.3,a)")
     &  "viscosity [Pa s]         : ",sthrmu(n)," (mu0)"
      write(ifll,"(3x,a,f10.3,a)")
     &  "Schmidt number           : ",sclmnr(n)," (Sc)"
      write(ifll,"(3x,a,f10.3,a)")
     &  "Prandtl number           : ",prlmnr(n)," (Pr)"
      write(ifll,"(5x,a)")
     &  "viscosity    : m = mu0"
      write(ifll,"(5x,2a)")
     &  "conductivity : l = Cp*m/Pr",
     &  "  (Cp:isopiestic specific heat)"
      write(ifll,"(5x,2a)")
     &  "diffusivity  : D = m/(r*Sc)",
     &  "  (r:mass density)"
      end subroutine disp_constant
!
!============================================
      subroutine disp_Simplified
!============================================
      write(ifll,"(3x,a,f10.3,a)")
     &  "Prandtl number           : ",prlmnr(n)," (Pr)"
      write(ifll,"(5x,a)")
     &  "viscosity    : m = a*Pr"
      write(ifll,"(5x,2a)")
     &  "conductivity : l = a*Cp",
     &  "  (Cp:isopiestic specific heat)"
      write(ifll,"(5x,2a)")
     &  "diffusivity  : D = a/(r*Le)",
     &  "  (r:mass density, Le:Lewis number)"
      write(ifll,"(10x,a)")
     &  "a = 2.58e-5*(T/298)^0.7"
      end subroutine disp_Simplified
!
      end subroutine inputdata
!
      end module module_material

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_gravity
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_gravity)'
!
!
!< module for gravity >
!
      real*8 :: ggg(3)=(/0.d0,0.d0,0.d0/)
      real*8 :: ramb=1.d0
      real*8 :: tamb=273.15d0, beta=0.21d-3,beta2=1.d0/273.d0
!
! ggg  : x,y,z components of gravity
! ramb : ambient density [kg/m^3]
! tamb : ambient temperature [K]
! beta : expansion coefficient [1/K]
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,lcomp,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      logical,intent(in)   :: lcomp
      integer,intent(out)  :: ierror
!
! [namelist]
      real*8 :: g(3)
      real*8 :: rho,t,rho2,t2
      integer :: ave_rho_T=0,p,ys(200)
      namelist /gravity/ g,rho,t,beta,beta2,rho2,t2,ave_rho_T,p,Ys
!
! [local entities]
      real*8,parameter :: undef=-huge(1.d0)
      real*8  :: betax
      integer :: iset=0
      integer :: ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      g   =ggg
      rho =ramb
      t   =tamb
!
      rewind ifli
      read(ifli,gravity,iostat=ios)
      if( ios.lt.0 ) return
!      write(ifll,gravity)
      call nml_errmsg0(ifle,ios,'gravity',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if( lcomp ) then
        if( rho.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'rho = ',rho
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
      else
        if( t.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 't = ',t
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
        if( beta.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'beta = ',beta
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
      endif
!
      ggg =g
      ramb=rho
      tamb=t
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_gravity

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_io
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(15),parameter,private :: modnam='(module_io)'
!
!
!< module for logical unit no. & file name >
!
!-< 1.logical unit no. for control data (input) >-
!
      integer,parameter :: ifli=1
      character(9),save :: cntlnam='fflow.ctl'
!
!-< 2.logical unit no. for log file & error message (output) >-
!
      integer,save      :: ifll=6, ifle=6
!
!-< 3.logical unit no. & name for external files >-
!
      integer,parameter              :: lenfnm=200
      integer,parameter,private      :: mfile=22
      integer          ,save,private :: iflno (-1:mfile)
      character(lenfnm),save,private :: filnam(-1:mfile)
      character(11)    ,save,private :: format(-1:mfile)
      character(7)     ,save,private :: status(-1:mfile)
      character(lenfnm),parameter,private :: blank=' '
!
! lenfnm : length of external file name
! mfile  : upper limit for no. of external files
! iflno  : logical unit no. for external file
! filnam : name of external file
! format : format identifier
! status : status identifier
!
!-< 99. for internal check >-
!
      integer,private :: iset=0
      character*4     :: gdformat   ! GF:GF format  FV:FieldView format
      real(8)         :: gdScale,angle_ssf,tolerance=1.d-7
      logical         :: GFrslt,FFrslt,FVrslt
!
!-------------------------------------------------------------
!     Added to output  FFR-GRID file
!     Flag for converting grid file from SC or FV to FFR
!-------------------------------------------------------------
      logical :: ffrgridOutFlag
!-------------------------------------------------------------
!     FFR grid format A: ASCII B: Bindary U: Unformatted
!-------------------------------------------------------------
      character*1 :: ffgFormatKey
!-------------------------------------------------------------
!     multi-material grid flag for FV-formatted grid
!-------------------------------------------------------------
      logical :: multiMaterialFlag
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(smallNrm,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      real*8,intent(inout) :: smallNrm
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      character(lenfnm) :: log,errmsg
      character(lenfnm) :: cell,vertex,boundary,initial
      character(lenfnm) :: p_initial,
     &                     p_restart,
     &                     p_result
      character(lenfnm) :: restart,result,history,force,anim
      character(lenfnm) :: rsltgf,rsltfv
      character(lenfnm) :: scInp,scVrt,scCel,scBnd
      character(lenfnm) :: fluent_grd
      character(lenfnm) :: scryu_grd
      character(lenfnm) :: nastran_grd
      character(lenfnm) :: cdcl,probe,fforce    !wallvar
      character(lenfnm) :: ffrgrid
      character*1       :: fgConv,ffrgridform
      character*3       :: multi_material
      namelist /files/ log,errmsg,
     &                 cell,vertex,boundary,
     &                 p_initial,
     &                 p_restart,
     &                 p_result,
     &                 initial,
     &                 restart,result,
     &                 history,force,anim,
     &                 gdformat,
     &                 gdScale,angle_ssf,tolerance,
     &                 scInp,scVrt,scCel,scBnd,  !
     &                 fluent_grd,               !
     &                 scryu_grd,                !
     &                 nastran_grd,              !
     &                 ffrgrid,ffrgridform,fgConv,
     &                 cdcl,probe,fforce
     &                 ,multi_material
!-----------------------
! --- [local entities]
!-----------------------
      character( 9),parameter :: formatted='formatted'
      character(11),parameter :: unformatted='unformatted'
      character( 3),parameter :: old='old', new='new'
      character( 7),parameter :: unknown='unknown'
      integer :: i,ios,no
      logical :: fexist=.false.
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
!-< 1. Input namelist >-
!
      log     =blank
      errmsg  =blank
      cell    =blank
      vertex  =blank
      boundary=blank
!      source  =blank
      initial =blank
      restart =blank
      result  =blank
      history =blank
      force   =blank
      anim    =blank
      gdformat='GF'
      rsltgf=blank
      rsltfv=blank
      GFrslt=.false.
      FFrslt=.false.
      FVrslt=.false.
      gdScale=1.d0
      angle_ssf=20.D0
      scInp  =blank
      scVrt  =blank
      scCel  =blank
      scBnd  =blank
      fluent_grd = blank
      scryu_grd = blank
      nastran_grd=blank
      ffrgrid=blank
      ffrgridOutFlag=.false.
      fgConv='n'
      ffrgridform='U'
      cdcl= blank
      probe= blank
      fforce=blank
      multi_material=blank
!
! --- inquire control file name
!
      inquire(file=cntlnam,exist=fexist)
      if(.not.fexist) then
        write(ifle,*) 'Control file [',cntlnam,'] NOT exist'
        stop ': Please check control file name'
      endif
      open(ifli,file=cntlnam,form='FORMATTED')
!----------------------------------------------
! --- control file (cntlnam='fflow.ctl' or fort.88)
!----------------------------------------------
      rewind ifli
      read (ifli,files,iostat=ios)
! ----------------------------------------------
      if( ios.ne.0 ) then
        if( ios.gt.0 ) write(*,*) '### error : read error'
        if( ios.lt.0 ) write(*,*) '### error : end of file'
        write(*,*) 'namelist = files'
        write(*,*) 'iostat =',ios
        goto 9999
      endif
      smallNrm=tolerance
      if(angle_ssf>90.d0.or.angle_ssf<0.d0) then
        write(ifle,*)'ERR: angle_ssf: ',angle_ssf
        write(ifle,*)'MSG: angle_ssf should be [0,90] (deg.)' 
      endif
!----------------------------------------------
! --- flag of 'if converting to FFR-GF'
!----------------------------------------------
      if((fgConv=='y').or.(fgConv=='Y')) ffrgridOutFlag=.true.
!
      if((ffrgridform=='a').or.(ffrgridform=='A')) then
        ffgFormatKey='A'
      else if((ffrgridform=='b').or.(ffrgridform=='B')) then
        ffgFormatKey='B'
      else if((ffrgridform=='u').or.(ffrgridform=='U')) then
        ffgFormatKey='U'
      else
        ffgFormatKey='U'
      end if
!-----------------------
! --- CHECK GRID FORMAT
!-----------------------
      if(gdformat=='FV') then
        call errmsg1(gdformat,vertex  ,'vertex')
      else if(gdformat=='SC') then
        call errmsg1(gdformat,scInp,'StarCD-Input')
        call errmsg1(gdformat,scVrt,'StarCD-Vertex')
        call errmsg1(gdformat,scBnd,'StarCD-Boundary')
        call errmsg1(gdformat,scCel,'StarCD-Cell')
! DEBUG
      else if(gdformat=='GB') then
        call errmsg1(gdformat,fluent_grd,'fluent_grd')
!      else if(gdformat=='GB-T'.or.gdformat=='GB-H') then
!        call errmsg1(gdformat,fluent_grd,'fluent_grd')
      else if(gdformat=='SY') then
        call errmsg1(gdformat,scryu_grd,'scryu_grd')
      else if(gdformat=='GF') then
        call errmsg1(gdformat,ffrgrid,'FFR_Grid')
      else if(gdformat=='NST') then
        call errmsg1(gdformat,nastran_grd,'nastran_grd')
      else
        write(*,*)'Assigned grid format keyword',
     &                        trim(gdformat),' is wrong.'
        write(*,*)'Check simulation control file.'
        ierror=1
        return
      endif
!
! --- 
!
! DEBUG
      if((ffrgridOutFlag).and.
     &  ((gdformat/='FV')
     &   .and.(gdformat/='SC')
     &   .and.(gdformat/='GB')
     &   .and.(gdformat/='NST')
     &  ))then
!      if((ffrgridOutFlag).and.
!     &  ((gdformat/='FV')
!     &   .and.(gdformat/='SC')
!     &   .and.(gdformat/='GB-H')
!     &   .and.(gdformat/='NST')
!     &  ))then

!--------------------------------------------------
! --- this version can convert to FFR-GRID from 
!     only Star-CD and FieldView grid.
!--------------------------------------------------
        ffrgridOutFlag=.false.
      end if
!
      if(ffrgridOutFlag) then
! --- when converting flag fgConv is YES, ffrgrid file name is required.
        call errmsg1(gdformat,ffrgrid,'FFR_Grid')
      end if
!
      if(result/=' ') FFrslt=.true.
      if(rsltgf/=' ') GFrslt=.true.
      if(rsltfv/=' ') FVrslt=.true.
      if( ierror.ne.0 ) goto 9999

      multiMaterialFlag=.false.
      if(multi_material=='yes') multiMaterialFlag=.true.
!
!
!-< 2. Set unit no., file name & some specifiers >-
!
      do 100 i=-1,mfile
      iflno (i)=i+12
      filnam(i)=blank
      format(i)=' '
      status(i)=' '
  100 continue
!
      
!
      no=-1
      call setspc(no,log,formatted,unknown)
      if( log.ne.blank ) ifll=iflno(no)
!
      no=0
      call setspc(no,errmsg,formatted,unknown)
      if( errmsg.ne.blank ) ifle=iflno(no)
!
      call setspc( 1,cell    ,unformatted,old)
      call setspc( 2,vertex  ,unformatted,old)
      call setspc( 3,boundary,unformatted,old)
      call setspc( 4,p_initial ,unformatted,old)
      call setspc( 5,initial ,unformatted,old)
      call setspc( 6,restart ,unformatted,unknown)
      call setspc( 7,result  ,unformatted,unknown)
      call setspc( 8,history ,unformatted,unknown)
      call setspc( 9, force  ,unformatted,unknown)
      call setspc( 10,rsltgf ,unformatted,unknown)
      call setspc( 11,rsltfv ,unformatted,unknown)
      call setspc( 12,scInp  ,  formatted,old)
      call setspc( 13,scVrt,    formatted,old)
      call setspc( 14,scBnd  ,  formatted,old)
      call setspc( 15,scCel  ,  formatted,old)
      call setspc( 16,fluent_grd,formatted,old)
      call setspc( 17,anim    ,unformatted,unknown)
      call setspc( 18,scryu_grd ,formatted,old)
      if(ffgFormatKey=='U') then
        call setspc( 19,ffrgrid,unformatted,old)
      else if(ffgFormatKey=='A') then
        call setspc( 19,ffrgrid,formatted,old)
      end if
      call setspc( 21,fforce,formatted,old)
      call setspc( 22,nastran_grd,formatted,old)
!
      return
!
 9999 continue
      write(*,*) modnam,subnam
      ierror=1
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1.1 error message >
!=================================================
      subroutine errmsg1(vnam0,vnam1,vnam2)
!=================================================
      character(*) :: vnam0,vnam1,vnam2
      if( vnam1.ne.' ' ) return
      write(*,'(1X,6a)') 'ERR: The grid format is defined as:',
     &          trim(vnam0),
     &   ', Set [',trim(vnam2),"='file_name'] in [&file] in ",
     &   cntlnam
      write(*,'(1X,4a)') 'lack of data : [',trim(vnam2),'] in ',cntlnam
      ierror=1
      end subroutine errmsg1
!
!< #1.2 set specifier >
!=================================================
      subroutine setspc(no,fin,fm,sta)
!=================================================
      integer,intent(in) :: no
      character(*),intent(in) :: fin,fm,sta
      filnam(no)=fin
      format(no)=fm
      status(no)=sta
      end subroutine setspc
!
      end subroutine inputdata
!
!
!< #2. open files >-----------------------------------------------------
!=================================================
      subroutine openfiles(ierror)
!=================================================
      integer,intent(out) :: ierror
      integer :: i,it,ios
      call modutl_chkset('=',2,iset,modnam//'(openfiles)')
      ierror=0
CCYYGF      do 100 i=1,mfile
CCYYGF      if( filnam(i).ne.blank ) then
CCYYGF      open(iflno(i),file=filnam(i),form=format(i),
CCYYGF     &     status=status(i),iostat=ios)
CCYYGF      if( ios.ne.0 ) then
CCYYGF        it=len_trim(filnam(i))
CCYYGF        write(*,*) '### error : file open error'
CCYYGF        write(*,*) 'file name = ',filnam(i)(:it)
CCYYGF        write(*,*) 'iostat =',ios
CCYYGF        ierror=1
CCYYGF        return
CCYYGF      endif
CCYYGF      endif
CCYYGF  100 continue
      end subroutine openfiles
!
!
!< #3. put logical unit no. & file name >-------------------------------
!=================================================
      subroutine getfil(ifl,fnam,cnam)
!=================================================
      character(8),parameter :: subnam='(getfil)'
      character(*),intent(in)  :: cnam
      integer     ,intent(out) :: ifl
      character(*),intent(out) :: fnam
      integer :: jfl
      call modutl_chkset('>',2,iset,modnam//subnam)
                               jfl=-999
      if( cnam.eq.'cell' )     jfl=1
      if( cnam.eq.'vertex' )   jfl=2
      if( cnam.eq.'boundary' ) jfl=3
      if( cnam.eq.'source' )   jfl=4
      if( cnam.eq.'initial' )  jfl=5
      if( cnam.eq.'restart' )  jfl=6
      if( cnam.eq.'result' )   jfl=7
      if( cnam.eq.'history' )  jfl=8
      if( cnam.eq.'force'   )  jfl=9
      if( cnam.eq.'rsltgf' )   jfl=10
      if( cnam.eq.'rsltfv' )   jfl=11
      if( cnam.eq.'scInp' )    jfl=12
      if( cnam.eq.'scVrt' )    jfl=13
      if( cnam.eq.'scBnd' )    jfl=14
      if( cnam.eq.'scCel' )    jfl=15
      if( cnam.eq.'fluent_grd')jfl=16
      if( cnam.eq.'anim')      jfl=17
      if( cnam.eq.'scryu_grd') jfl=18
      if( cnam.eq.'ffrgrid')   jfl=19
      if( cnam.eq.'fforce')   jfl=21
      if( cnam.eq.'nastran_grd')jfl=22
      if( jfl.eq.-999 ) then
        write(*,*) '### error : program error -1- ',modnam,subnam
        write(*,*) 'cnam = ',cnam
        stop 999
      endif
      ifl=iflno(jfl)
      fnam=filnam(jfl)
      end subroutine getfil
!
      end module module_io

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_les
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(16),parameter,private :: modnam='(module_les)'
!
!
!< module for LES model >
!
      real*8 :: csles=0.1d0
      real*8 :: apls =25.d0
      real(8) :: dlalpha=2.d0
      real(8) :: ign_iv=0
      integer:: ISTART=1000
      integer,save             :: iuvw_ave_rms_re
      integer,save             :: ip_ave,it_ave
      integer,save,allocatable :: icomp_ave(:),irans_ave(:)
      integer,save             :: it_rms,ip_rms
      integer,save,allocatable :: icomp_rms(:),irans_rms(:)
      integer,save             :: iuvwt_re
      integer,save,allocatable :: iuvw_rans_re(:)
!
      integer,save             :: iuvw_ave_rms_re_c
      integer,save             :: ip_ave_c,it_ave_c
      integer,save,allocatable :: icomp_ave_c(:),irans_ave_c(:)
      integer,save             :: it_rms_c,ip_rms_c
      integer,save,allocatable :: icomp_rms_c(:),irans_rms_c(:)
      integer,save             :: iuvwt_re_c
      integer,save,allocatable :: iuvw_rans_re_c(:)
!
! csles : Smagorinsky's constant
! ISTART: Start iteration for average value
! iu_ave: flag for if/not computing averaged value for LES (velocity)
! ip_ave: flag for if/not computing averaged value for LES (pressure)
! icomp_ave: flag for if/not computing averaged value for LES (speciese)
! irans_ave: flag for if/not computing averaged value for LES (scalar )
! iu_ave=0/1: 0=>not calculating ; 1=>calculating
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,nrans,ncomp,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle,nrans,ncomp
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      real*8 :: cs
      integer:: NSTART
      integer,parameter        :: msuf_out=11
      integer:: uvw_ave_rms_re=1,
     &          average_p=0,average_t=0,average_Mut=0,
     &          average_comp(200),average_rans(200),
     &          rms_t=0,rms_p=0,rms_Mut=0,
     &          rms_comp(200),rms_rans(200),
     &          re_uvw_t=0,
     &          re_uvw_rans(200),statis_file

      integer ::  re_uvw_comp(200)

      integer:: 
     &          surface_SDOT=0,
     &          GAS_WDOT=0,maxmin_uvw_p=0,
     &          depos_etch_speed=0,
     &          surface_mole_fraction=0,
     &          gas_mole_fraction=0,
     &          bulk_thickness=0
     &          
      


      character(LEN=80) :: RATE_EQ_SpecNAME(msuf_out)
      namelist /les/ cs,NSTART,dlalpha,ign_iv,
     &               uvw_ave_rms_re,maxmin_uvw_p,
     &               average_p,average_t,average_mut,
     &               average_comp,average_rans,
     &               rms_t,rms_p,rms_mut,
     &               rms_comp,rms_rans,
     &               re_uvw_t,
     &               re_uvw_rans,
     &               re_uvw_comp,
     &               surface_SDOT,gas_WDOT,depos_etch_speed,
     &               RATE_EQ_SpecNAME,surface_mole_fraction,
     &               gas_mole_fraction,
     &               bulk_thickness,statis_file

!
! --- [local entities]
!
      integer :: iset=0
      integer :: ios,ierr1
!
      allocate(icomp_ave(ncomp),
     &         irans_ave(nrans),
     6         icomp_rms(ncomp),
     &         irans_rms(nrans),
     &         iuvw_rans_re(nrans))
!
      allocate(icomp_ave_c(ncomp),
     &         irans_ave_c(nrans),
     6         icomp_rms_c(ncomp),
     &         irans_rms_c(nrans),
     &         iuvw_rans_re_c(nrans))
!
      iuvw_ave_rms_re=1
      ip_ave=0
      it_ave=0
      icomp_ave(:)=0
      irans_ave(:)=0
      ip_rms=0
      it_rms=0
      icomp_rms(:)=0
      irans_rms(:)=0
      iuvwt_re=0
      iuvw_rans_re(:)=0
!
      iuvw_ave_rms_re_c=1
      ip_ave_c=0
      it_ave_c=0
      icomp_ave_c(:)=0
      irans_ave_c(:)=0
      ip_rms_c=0
      it_rms_c=0
      icomp_rms_c(:)=0
      irans_rms_c(:)=0
      iuvwt_re_c=0
      iuvw_rans_re_c(:)=0
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      cs=csles
      NSTART=1000
      uvw_ave_rms_re=1
      average_p=0
      average_t=0
      average_comp(:)=0
      average_rans(:)=0
      rms_p=0;rms_t=0
      rms_comp(:)=0;rms_rans(:)=0
      re_uvw_t=0
      re_uvw_rans(:)=0
!
      rewind ifli
      read(ifli,les,iostat=ios)
      if( ios.lt.0 ) return
!      write(ifll,les)
      call nml_errmsg0(ifle,ios,'les',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if( cs.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'cs = ',cs
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      csles=cs
      ISTART=NSTART
      iuvw_ave_rms_re=uvw_ave_rms_re
      ip_ave=average_p
      it_ave=average_t
      icomp_ave(1:ncomp)=average_comp(1:ncomp)
      irans_ave(1:nrans)=average_rans(1:nrans)
      ip_rms=rms_p
      it_rms=rms_t
      icomp_rms(1:ncomp)=rms_comp(1:ncomp)
      irans_rms(1:nrans)=rms_rans(1:nrans)
      iuvwt_re=re_uvw_t
      iuvw_rans_re(1:nrans)=re_uvw_rans(1:nrans)
!
      if(iuvwt_re.eq.1.and.iuvw_ave_rms_re.eq.0) then 
        write(ifle,*) '### error : data error'
        write(ifle,*) '### reset [iuvw_ave_rms_re=1] in ',cntlnam,
     &  ' while [iuvwt_re=1]'
        stop ': at module_les'
      endif
      do i=1,nrans
      if(iuvw_rans_re(i)/=0.and.iuvw_ave_rms_re.eq.0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) '### reset [uvw_ave_rms_re=1] in ',cntlnam,
     &  'while [re_uvw_rans(i)=1]'
        stop ': at module_les'
      endif
      enddo
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_les

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_model
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(18),parameter,private :: modnam='(module_model)'
!
!
!< module for model flags >
!
!-< 1. Flow model >-
!
      integer,parameter :: incomp=0, mach0=1, comp=2
      integer :: idrdp=mach0,ICVREF=1
      integer :: ical_mvmsh
      integer,parameter :: PEFC=1,PAFC=2,MCFC=3,SOFC=4
      integer :: IFUCL=0
!
! incomp : incompressible
! mach0  : zero Mach no. approximation
! comp   : compressible
! idrdp  : flag for incompressible/compressible
!
!
!-< 2. Turbulence model >-
!
      integer,parameter :: noturb=0, ke=1,sles=2,dles=3, lles=4,
     &                     dns=5,RSM=6,ke_low=7,RNG=8,CHEN=9,
     &                     SDES=10,cnst=11,KE2S=12,KLES=13
      integer :: icaltb=ke
      logical :: rns_scl=.false.
!
      integer,parameter :: vertex_cen=1,cell_cen=2
      integer :: icon_cv=vertex_cen
!
! noturb : no model
! ke   : k-epsilon model
! sles   : LES model,Smagorinsky LES model
! dles   : Dynamic LES model
! lles   : Lagrangian Dynamic LES model
! dns    : Direct Numerical Simulation (same as 'noturb')
! RSM    : Reynolds Stress Equation Model
! 
! icaltb : flag for turbulence model
!
!
      integer,save :: ical_MHD,KMHD
      real*8 ,save :: Buffle_shft=0.d0
!
!-----------------------
! --- user function
!-----------------------
      integer,parameter :: uf_mm=200
      integer,save :: u_func(uf_mm)      
! 
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!================================================================
      subroutine inputdata
     &(ifli,ifll,ifle,lcomp,lrans,kemdl,lowke,lRNG,lCHEN,l2SKE,RANSMDL,
     & lSDES,lFC,lKLES,
     & RSMmdl,nrans,cntlnam,idens,ierror)
!================================================================
      implicit none
      character(11),parameter :: subnam='(inputdata)'
!
! 1. Input flags for model
!
! [dummy arguments]
!
      integer,intent(in)     :: ifli,ifll,ifle
      integer,intent(inout)  :: nrans,lFC,idens
      character(*),intent(in):: cntlnam
      logical,intent(out)    :: lcomp,lrans,kemdl,RSMmdl,lowke,
     &                          lRNG,lCHEN,lSDES,RANSMDL,l2SKE,lKLES
      integer,intent(out)    :: ierror
!
! --- [namelist]
!
      character(20) :: flow,trbmdl
      integer       :: cellref,cal_tys=1,cal_t=1,cal_reac=1,cal_surf=1,
     &                 cal_weak=0,cal_MHD=0,cal_vect=0,cal_mvmsh=0,
     &                 LBF_P=0,HPC_close=0
      real*8        :: pref,Buffle_shift=0.d0
      integer       :: NODE_CPU=1,moni_inter=100,iFIELD=0,EWT=0
      integer       :: density=0,cal_heat_reaction=1!,Potential=0
      character(20) :: Fuel_cell
      integer       :: cal_particle=0,cal_tet,cal_thermo_diff,rain_WDR
      character(LEN=6) :: CV
!
      integer       :: userflag(1:uf_mm)=0,iu
      namelist /user_function/ userflag
!
      namelist /model/ flow,trbmdl,
     &                cal_tys,cal_t,cal_reac,cal_surf,cal_weak,
     &                cal_MHD,cal_vect,NODE_CPU,moni_inter,cal_mvmsh,
     &                LBF_P,HPC_close,density,
     &                Fuel_cell,iFIELD,cal_heat_reaction,
     &                cal_particle,
     &                cal_tet,cal_thermo_diff,
     &                Buffle_shift,EWT,CV
!
!
! --- [local entities]
!
      character(14),parameter :: 
     &   flwlst(3)=(
     &    /'incompressible','zero-Mach     ','compressible  ' /)

      character( 4),parameter :: trblst(14)=(/
     & 'no  ','KE  ','SLES','DLES','LLES','DNS ','RSM ','LKE ',
     & 'RNG ','CHEN','SDES', 'cnst', '2SKE','KLES'/)

      character(4),parameter :: fulclist(6)=(/
     & 'PEFC','PAFC','MCFC','SOFC','AFC ','MDFC'/)
      character(6),parameter ::
     &   CVlst(2)=(/'vertex','cell  '/)
! --- PEFC: 
! --- PAFC: 
! --- MCFC: 
! --- SOFC: 
! --- AFC : 
      integer :: iset=0
      integer :: iflw,itrb,ios,ierr1,ierr2
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      flow  =' '
      CV=' '
      trbmdl=' '
      Fuel_cell=' '
      pref=101325.0
      cellref=1
      cal_tys=1
      cal_MHD=0
!      Potential=0
      density=0
!
      rewind ifli
      read(ifli,model,iostat=ios)
!      write(ifll,model)
      call nml_errmsg0(ifle,ios,'model',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      call nml_listno(3,flwlst,flow  ,idrdp)
      call nml_listno(14,trblst,trbmdl,icaltb)
      call nml_chkchr0(ifle,'flow'  ,flow  ,idrdp ,ierr1)
      call nml_chkchr0(ifle,'trbmdl',trbmdl,icaltb,ierr2)
      if( ierr1.ne.0 .or. ierr2.ne.0 ) goto 9999
      call nml_listno(2,CVlst,CV,icon_cv)
      call nml_chkchr0(ifle,'CV'  ,CV,icon_cv ,ierr1)
      if(icon_cv==cell_cen) then
        stop 'ERR: NOT support cell center in V3.1'
      endif
      if(ierr1.ne.0) goto 9999
!
      if(density==4) then
        idens=4
      endif
!
      if(Fuel_cell/=' ') then
        call nml_listno(6,fulclist,Fuel_cell,IFUCL)
        call nml_chkchr0(ifle,'Fuel_cell',Fuel_cell,IFUCL,ierr1)
        lFC=IFUCL
        write(ifle,'(1x,a)') "MSG: Fuel_cell= for NOT fuel-cell module"
        if(ierr1/=0) goto 9999
      endif
!
      idrdp=idrdp-1
      icaltb=icaltb-1
!
                          lcomp=.true.
      if(idrdp.eq.incomp) lcomp=.false.
!
!      if(.not.(idrdp==incomp.or.idrdp==mach0).and.cal_vect==1) then
!!      if(.not.(idrdp==comp).and.cal_vect==1) then
!        write(*,'(1X,2a)') 
!     &  'ERR: FFR Vector-Version ',
!     &   'ONLY support incompressible or zero Mach'
!          write(*,'(1X,a)') 
!     &  "MSG: set [flow='incomp']"
!          stop ' at [&model]' 
!      endif
!
      if(icaltb.eq.ke.or.icaltb.eq.RSM) then
!        lrans=.true.
        rns_scl=.true.
      else
!        lrans=.false.
        rns_scl=.false.
      endif
!
      kemdl=.false.
      RSMmdl=.false.
      lowke=.false.
      lRNG=.false.
      lCHEN=.false.
      lSDES=.false.
      lKLES=.false.
      RANSMDL=.false.
      if(icaltb.eq.ke) then
        nrans=2
        kemdl=.true.
        if(cal_vect==1) then
          write(*,'(1X,a)') 
     &  'ERR: FFR Vector Version NOT support RANS model'
          write(*,'(1X,a)') 
     &  "MSG: set [trbmdl='DNS'] or [trbmdl='SLES']"
!          stop ' at [&model]' 
        endif
        RANSMDL=.true.
      elseif(icaltb.eq.RSM) then
        nrans=7
        RSMmdl=.true.
        if(cal_vect==1) then
          write(*,'(1X,a)') 
     &  'ERR: FFR Vector Version NOT support RANS model'
          write(*,'(1X,a)') 
     &  "MSG: set [trbmdl='DNS'] or [trbmdl='SLES']"
          stop ' at [&model]' 
        endif
        RANSMDL=.true.
      elseif(icaltb.eq.ke_low) then
        nrans=2
        lowke=.true.
        if(cal_vect==1) then
          write(*,'(1X,a)') 
     &  'ERR: FFR Vector Version NOT support RANS model'
          write(*,'(1X,a)') 
     &  "MSG: set [trbmdl='DNS'] or [trbmdl='SLES']"
          stop ' at [&model]' 
        endif
        RANSMDL=.true.
      elseif(icaltb.eq.RNG) then
        nrans=2
        lRNG=.true.
        if(cal_vect==1) then
          write(*,'(1X,a)') 
     &  'ERR: FFR Vector Version NOT support RANS model'
          write(*,'(1X,a)') 
     &  "MSG: set [trbmdl='DNS'] or [trbmdl='SLES']"
          stop ' at [&model]' 
        endif
        RANSMDL=.true.
      elseif(icaltb.eq.CHEN) then
        nrans=2
        lCHEN=.true.
        if(cal_vect==1) then
          write(*,'(1X,a)') 
     &  'ERR: FFR Vector Version NOT support RANS model'
          write(*,'(1X,a)') 
     &  "MSG: set [trbmdl='DNS'] or [trbmdl='SLES']"
          stop ' at [&model]' 
        endif
        RANSMDL=.true.
      elseif(icaltb.eq.KE2S) then
        nrans=2
        l2SKE=.true.
        if(cal_vect==1) then
          write(*,'(1X,a)') 
     &  'ERR: FFR Vector Version NOT support RANS model'
          write(*,'(1X,a)') 
     &  "MSG: set [trbmdl='DNS'] or [trbmdl='SLES']"
          stop ' at [&model]' 
        endif
        RANSMDL=.true.
      elseif(icaltb.eq.SDES) then
        nrans=1
        lSDES=.true.
      elseif(icaltb==KLES) then
        nrans=1
        lKLES=.true.
      endif
!
      ICVREF=cellref
!
      if(.not.(cal_MHD==1.or.cal_MHD==0.or.cal_MHD==2))then
        write(ifle,'(1x,a)') 
     &   '    cal_MHD=1 : MHD (MagnetoHydroDynamics) Coupling '
        write(ifle,'(1x,a)') 
     &   '    cal_MHD=0 : No MHD Coupling '
        write(ifle,'(1x,a)') 
     &   '    cal_MHD=2 : Only MHD Solid '
        write(ifle,'(1x,a)') 
     &   '    ### EER: [cal_MHD] MUST be 1 or 0 '
        
        stop 'module_model'
      else
        ical_MHD=cal_MHD
        if(cal_MHD==1) then
          KMHD=1
        else
          KMHD=0
        endif
      endif
!
      if(.not.(cal_mvmsh==1.or.
     &         cal_mvmsh==0.or.
     &         cal_mvmsh==2.or.
     &         cal_mvmsh==3.or.
     &         cal_mvmsh==4.or.
     &         cal_mvmsh==5
     &  ))then
        write(ifle,'(6x,a)') 'EER: [cal_mvmsh] MUST be 0, 1 2 3 4 or 5'
!
        write(ifle,'(6x,2a)') 
     &   'MSG: cal_mvmsh=0 : ',
     &   'No [Moving] & [Removing/Adding] Mesh solver'
!
        write(ifle,'(6x,2a)') 
     &   'MSG: cal_mvmsh=1 : ',
     &   'Running [Moving] & [Removing/Adding] solver '
!
        write(ifle,'(6x,2a)') 
     &   'MSG: cal_mvmsh=2 : ',
     &   'Checking CV-GEOM of [Moving] mesh no-running solver'
!
        write(ifle,'(6x,4a)') 
     &   'MSG: cal_mvmsh=3 : ',
     &   'Viewing [Moving] & [Removing/Adding] mesh ',
     &   'no-checking CV-GEOM ',
     &   'and no-running solver'
!
        write(ifle,'(6x,3a)') 
     &   'MSG: cal_mvmsh=4 : ',
     &   'Checking CV-GEOM of [Moving] & [Removing/Adding] mesh ',
     &   'no-running solver'
!
        write(ifle,'(6x,2a)') 
     &   'MSG: cal_mvmsh=5 : ',
     &   'Running [Moving] & [Removing/Adding] mesh solver'
!        
        stop 'module_model'
      else
        ical_mvmsh=cal_mvmsh
      endif
!
!
      if(.not.(cal_tys.eq.2.or.cal_tys.eq.1.or.cal_tys.eq.0))then
        write(ifle,*) '    ### EER: [cal_tys] MUST be 1 or 0 '
        goto 9999
      endif
!
      if(IFUCL/=0) then
        write(ifll,'(1X,2a)') 'MSG: FULL-CELL Module: ',fulclist(IFUCL)
      endif
!
      Buffle_shft=0.d0
      if(abs(Buffle_shift)>1.d-15) then
        Buffle_shft=Buffle_shift
        write(ifll,*)
     &  'MSG: Buffle region z-coordinate been shifted :',Buffle_shift
        write(ifll,*) 'MSG: This flag is for Gridgen BUG'
      else
        write(ifll,*)
     &  'MSG: Buffle region z-coordinate NOT been shifted :'
        write(ifll,*) 'MSG: This flag is for Gridgen BUG'
      endif
!
!---------------------------
! --- 
!---------------------------
      userflag(:)=0
      u_func(:)=0
      rewind ifli
      read(ifli,user_function,iostat=ios)
      if( ios>0 ) goto 9999
      do iu=1,uf_mm
      if(userflag(iu)==1.and.iu==1) then
        write(ifll,'(1x,a,I4)') 'MSG: userflag(i) i= ',iu
        write(ifll,'(1x,a)') 'MSG: [atmospheric diffusion model] is ON'
        write(ifll,'(1x,2a)') 
     &     'MSG: [atmospheric diffusion model] is offered by ',
     &  'Center Research Institute of Electric Power Industry (CRIEPI)'
        u_func(iu)=1
      elseif(userflag(iu)==1.and.iu==2) then
        write(ifll,'(1x,a,I4)') 'MSG: userflag(i) i= ',iu
        write(ifll,'(1x,a)') 'MSG: [Multi region] flag is ON'
        write(ifll,'(1x,2a)') 
     &  'MSG: [Multi region] is offered by ',
     &  'Center Research Institute of Electric Power Industry (CRIEPI)'
        write(ifll,'(1x,a)') 
     & 'MSG: two region is defined by using [USER_UserDefine_region.f]'
        u_func(iu)=1
      elseif(userflag(iu)==1.and.iu==3) then
        write(ifll,'(1x,a,I4)') 'MSG: userflag(i) i= ',iu
        write(ifll,'(1x,a)') 'MSG: [Canopy region] flag is ON'
        write(ifll,'(1x,2a)') 
     &  'MSG: [Canopy region] is offered by ',
     &  'Center Research Institute of Electric Power Industry (CRIEPI)'
        write(ifll,'(1x,2a)') 
     &'MSG: Canopy region is ',
     &    'defined by using [USER_UserDefine_region.f]'
        u_func(iu)=1
      elseif(userflag(iu)==1) then
        write(ifle,'(1x,a,I4)') 
     &    'ERR: &user_function/userflag(i), i=',iu
        stop 'ERR: at &user_function'
      elseif(userflag(iu)/=1.and.userflag(iu)/=0) then
        stop
     &  'ERR-MSG: at &user_function/userflag(1:200) MUST be 0 or 1'
      else
        
      endif
      enddo      
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_model

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_movegrid
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(21),parameter,private :: modnam='(module_movegrid)'
!
!
!< module for moving grid >
!
      integer,parameter,private :: mgrid=10
!
      integer :: ngrid=0
      real*8  :: dtmin=0.d0
      real*8,save :: time(mgrid+1)
!
! ngrid : no. of grid tables
! dtmin : minimum time interval of grid tables
! time  : time of each grid
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! [namelist]
      real*8 :: cycle
      namelist /movegrid/ ngrid,time,cycle
!
! --- [local entities]
!
      integer :: iset=0
      integer :: i,inpt,ios,ierr1
      real*8  :: dum1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      ngrid= 1
      time = 0.d0
      cycle=-1.d0
      inpt = 1
!
      rewind ifli
      read(ifli,movegrid,iostat=ios)
      if( ios.lt.0 ) goto 100
!      write(ifll,movegrid)
      call nml_chksiz(ifle,'mgrid',mgrid,ngrid,modnam,ierr1)
      if( ierr1.ne.0 ) goto 9999
      call nml_errmsg0(ifle,ios,'movegrid',ierr1)
      if( ierr1.ne.0 ) goto 9999
      inpt=2

!      if(ngrid.le.1) then
!        return
!      endif

  100 continue
!
      if( ngrid.lt.inpt ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'ngrid =',ngrid
        write(ifle,*) 'it must be > 1'
        goto 9999
      endif
      dtmin=time(2)-time(1)
      do 110 i=2,ngrid
      dtmin=min(dtmin,time(i)-time(i-1))
      if( time(i-1).ge.time(i) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'time(i-1),time(i) = ',time(i-1),time(i)
        write(ifle,*) 'i =',i
        write(ifle,*) 'time(i-1) must be < time(i)'
        goto 9999
      endif
  110 continue
      dtmin=0.999d0*dtmin
!
      if( cycle.gt.0.d0 ) then
        dum1=time(ngrid)-time(1)
        if( cycle.le.dum1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'cycle = ',cycle
        write(ifle,*) 'it must be > time(ngrid)-time(1)'
        write(ifle,*) 'time(ngrid)-time(1) = ',dum1
        goto 9999
        endif
        time(ngrid+1)=time(1)+cycle
      else
        dum1=time(1)
        time(ngrid+1)=min(1.1d0*dum1,0.9d0*dum1,dum1-1.d0)
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_movegrid

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_nowtime
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!!      include'interfacemmm'
!
!< module for present time & iteration counts >
!
      integer :: iter=0
      real*8  :: time=0.d0
!
      end module module_nowtime

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_output
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(19),parameter,private :: modnam='(module_output)'
!
!
! 1. module for control data to output
!
      private
      public :: outchk,inputdata
!
!
      integer,parameter :: mspec=16, mfile=10
!
! mspec : upper limit for no. of specified output points
! mfile : no. of files
!
!
      character(LEN=9),parameter,public :: fnmlst(mfile)=(/
     &  'restart  ',
     &  'result   ',
     &  'history  ',
     &  'force    ',
     &  'anim     ',
     &  'cdcl     ',
     &  'probe    ',
     &  'fforce   ',
     &  'p_restart', 
     &  'p_result '/)

!      character(7),parameter,public :: fnmlst(mfile)=(/
!     & 'restart','result ','history','force  ','anim   ',
!     & 'cdcl   ','probe  ','fforce '/)
!
! fnmlst : list of file names
!
!
      integer,parameter,public :: none=-1, yes=1, no=0
!
! none : no output
! yes  : it is time to output
! no   : it is not time to output
!
!
      integer :: i
      integer,save,public :: kopt(mfile)
      real*8 ,save,public :: start_force
      real(8) ,save,public :: start_cdcl,start_prob,start_wallvar
      integer,save :: nspc(mfile)
      real*8 ,save :: tspc(mspec,mfile)
      real*8 ,save :: tps(mfile),tpe(mfile),tpd(mfile)
      logical :: set(mfile)=(/(.false.,i=1,mfile)/)
!
! kopt    : flag which indicates how to control
!         : =1; interval of iteration counts
!         : =2; specified iteration counts
!         : =3; interval of time
!         : =4; specified time
! nspc    : no. of specified output points
! tspc    : specified output points ( time or iteration counts )
! tps,tpe :
! & tpd   : start, end & interval of time or iteration counts
!
!
      integer :: iset=0
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!=================================================
      subroutine outchk(fnam,iter,time,delt,out)
!=================================================
      character(8),parameter :: subnam='(outchk)'
      character(*),intent(in)  :: fnam
      integer     ,intent(in)  :: iter
      real*8      ,intent(in)  :: time,delt
      integer     ,intent(out) :: out
      integer :: i,l,ns,nd,idum
      real*8  :: eps
      logical :: first(mfile)=(/(.true.,i=1,mfile)/)
!
!
!
      call modutl_chkset('>',1,iset,modnam//subnam)
      out=none
!
      call nml_listno(mfile,fnmlst,fnam,l)
      if( l.lt.1 ) then
        write(*,*) '### program error -1- ',modnam,subnam
        write(*,*) 'fnam = ',fnam
        stop 999
      endif
      if( .not.set(l) ) return
!
      if( kopt(l).gt.2 ) then
      eps=0.1d0*delt
      if( eps.le.0.d0 ) then
        write(*,*) '### program error -2- ',modnam,subnam
        write(*,*) 'eps = ',eps
        stop 999
      endif
      endif
!
      out=no
!
!/interval of iteration counts/
      if( kopt(l).eq.1 ) then
        if( iter.gt.nint(tpe(l)) ) return
        if( first(l) ) then
          ns=nint(tps(l))
          nd=nint(max(1.d0,tpd(l)))
          if( iter.gt.ns ) then
            ns=ns+((iter-ns)/nd)*nd
            if( iter.gt.ns ) ns=ns+nd
          endif
          tps(l)=ns
          tpd(l)=nd
          first(l)=.false.
        endif
        if( iter.lt.nint(tps(l)) ) return
        tps(l)=nint(tps(l)+tpd(l))
!/specified iteration counts/
      elseif( kopt(l).eq.2 ) then
        do 100 i=1,nspc(l)
        if( iter.eq.nint(tspc(i,l)) ) goto 101
  100   continue
        return
  101   continue
!
!/interval of time/
      elseif( kopt(l).eq.3 ) then
        if( time.gt.tpe(l)+eps ) return
        if( first(l) ) then
          tpd(l)=max(0.d0,tpd(l))
          if( time.gt.tps(l) .and. tpd(l).gt.0.d0 ) then
            tps(l)=tps(l)+int((time-tps(l))/tpd(l))*tpd(l)
            if( time.gt.tps(l) ) tps(l)=tps(l)+tpd(l)
          endif
          first(l)=.false.
        endif
        if( time.lt.tps(l)-eps ) return
        tps(l)=tps(l)+tpd(l)
!/specified time/
      else
        do 110 i=1,nspc(l)
        if( abs(time-tspc(i,l)).lt.eps ) goto 111
  110   continue
        return
  111   continue
      endif
!
      out=yes
      end subroutine outchk
!
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! [namelist]
      character(20) :: file
      character(10) :: type
      integer :: nspec,multi_result=0,BC_result=0
      real*8  :: start,end,inter
      real*8  :: spec(mspec)
      character(20)     :: filetype
!
      namelist /output/ file,type,start,end,inter,nspec,spec,
     &                  multi_result,filetype,BC_result
!
! [local entities]
!
      integer,parameter :: iundef=-huge(1)
      real*8 ,parameter :: undef=-huge(1.d0)
      character(7),parameter :: typlst(4)=(/
     & 'inter_i','spec_i ','inter_t','spec_t '/)
      integer :: lchk(mfile)
      integer :: i,nseq,nfl,ntp
      integer :: ios,ierr1,ierr2
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      type =' '
      start=undef
      inter=undef
      spec =undef
      nspec=iundef
!
      lchk=0
      nseq=0
      rewind ifli
 1000 continue
!
      file=' '
      end=huge(1)
!
      read(ifli,output,iostat=ios)
      if( ios.lt.0 ) goto 1001
      call nml_chksiz(ifle,'mspec',mspec,nspec,
     & '(module_output)',ierr1)
      if( ierr1.ne.0 ) goto 9999
!      write(ifll,output)
      call nml_errmsg0(ifle,ios,'output',ierr1)
      if( ierr1.ne.0 ) goto 9999
      nseq=nseq+1
!
      call nml_listno(mfile,fnmlst,file,nfl)
      call nml_listno(4,typlst,type,ntp)
      call nml_chkchr0(ifle,'file',file,nfl,ierr1)
      call nml_chkchr0(ifle,'type',type,ntp,ierr2)
      if( ierr1.ne.0 .or. ierr2.ne.0 ) goto 9999
!
      if( lchk(nfl).gt.0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'file = ',file
        write(ifle,*) 'this file already specified in previous line'
        write(ifle,*) 'sequence no. of the line =',lchk(nfl)
        goto 9999
      endif
      lchk(nfl)=nseq
!
      if( ntp.eq.1 .or. ntp.eq.3 ) then
        if( start.eq.undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : start'
          goto 9999
        endif
        if( end.lt.start ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'end = ',end
          write(ifle,*) 'it must be >= start'
          goto 9999
        endif
        if( inter.eq.undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : end'
          goto 9999
        elseif( inter.lt.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'inter = ',inter
          write(ifle,*) 'it must be >= 0'
          goto 9999
        endif
      else
        if( nspec.eq.iundef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : nspec'
          goto 9999
        elseif( nspec.lt.1 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'nspec =',nspec
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
        do 100 i=1,nspec
        if( spec(i).eq.undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : spec(i)'
          write(ifle,*) 'i =',i
          goto 9999
        endif
  100   continue
      endif
!
      set(nfl)=.true.
      kopt(nfl)=ntp
      if( ntp.eq.1 .or. ntp.eq.3 ) then
        tps(nfl)=start
        tpe(nfl)=end
        tpd(nfl)=inter
      else
        nspc(nfl)=nspec
        do 200 i=1,nspec
        tspc(i,nfl)=spec(i)
  200   continue
      endif
!
      goto 1000
 1001 continue
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_output

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_param 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!< module for named constants >
!
      real*8 :: yslw=1.d-20
      real*8,allocatable :: scalfctr(:)
! scalfctr(0) : max model size (x,y,z box size)
! scalfctr(:) : min face size for each boundary (min length on 1 face)
      real(8) :: smallNrm=1.0d-7
!
! yslw : lower limit of mass fraction
!
      end module module_param

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_rans
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_rans)'
!
!
!
!< 1. module for RANS model >
!
      integer,private :: nrans=0
!
      real*8,save,allocatable :: akslw(:)
      real*8,save,allocatable :: sgmaks(:)
!
! akslw  : lower limit of depedent variables
! sgmaks : parameter in diffusion term
!
!
!< 2. for k-epsilon model >
!
      real*8 :: cep1=1.45d0, cep2=1.90d0, cep3=1.45d0,cep4=-0.33d0
      real*8 :: cmu=0.09d0
      real*8 :: sgmk=1.4d0, sgme=1.3d0
      real(8) :: Schmidt_G=0.25d0, Schmidt_Z=0.5,
     &           Schmidt_G_M=1.d10,Schmidt_Z_M=0.7
!
! cep1 : parameter in generation term of eps. eq. [-]
! cep2 : parameter in generation term of eps. eq. [-]
! cep3 : parameter in generation term of eps. eq. [-]
! cmu  : parameter in turbulent viscosity [-]
! sgmk : parameter in diffusion term of k eq. [-]
! sgme : parameter in diffusion term of eps. eq. [-]
!
!
!-< 99. for internal check >-
!
      integer,private :: iset=0
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputke(ifli,ifll,ifle,cntlnam,nrnsx,kemdl,
     &   lowke,lRNG,lCHEN,l2SKE,lKLES,ierror)
!=================================================
      character(9),parameter :: subnam='(inputke)'
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle,nrnsx
      character(*),intent(in) :: cntlnam
      logical,intent(in)   :: kemdl,lowke,lRNG,lCHEN,l2SKE,lKLES
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      real*8 :: klow,epslow,sgm_C_dash,C_k,C_e
      namelist /kemodel/ klow,epslow,cep1,cep2,cep3,cep4,cep5,
     &                   cmu,sgmk,sgme,sgmm,sgmh,eta0,beta,
     &                   Schmidt_G,Schmidt_Z,
     &                   Schmidt_G_M,Schmidt_Z_M,sgm_C_dash,
     &                   C_k,C_e
!
! --- [local entities]
!
      real*8,parameter :: undef=-huge(1.d0)
      integer :: ios=0,ierr1
!
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
      ierror=0
      nrans=nrnsx
      if(lowke.or.kemdl.or.lRNG.or.lCHEN.or.l2SKE.or.nrans.gt.0) then
!        if( nrans.ne.2 ) then
!          write(*,*) '### program error -1- ',modnam,subnam
!          stop 999 ???
!        endif
      else
        return
      endif
!
      allocate( akslw(nrans), sgmaks(nrans), stat=ierr1 )
      if( ierr1.gt.0 ) goto 9999
!
      klow  =undef
      epslow=undef
      rewind ifli
      read(ifli,kemodel,iostat=ios)
!      write(ifll,kemodel)
      call nml_errmsg0(ifle,ios,'kemodel',ierr1)
      if( ios.gt.0 ) goto 9999
!
      if((kemdl.or.lowke.or.lRNG.or.lCHEN.or.l2SKE).and.nrans.gt.0) then
        if( klow.eq.undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : klow'
        elseif( klow.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'klow = ',klow
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
!
        if( epslow.eq.undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : epslow'
        elseif( epslow.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'epslow = ',epslow
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
!
        if( cep1.lt.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'cep1 = ',cep1
          write(ifle,*) 'it must be >= 0'
          goto 9999
        endif
!
        if( cep2.lt.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'cep2 = ',cep2
          write(ifle,*) 'it must be >= 0'
          goto 9999
        endif
!
        if( cep3.lt.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'cep3 = ',cep3
          write(ifle,*) 'it must be >= 0'
          goto 9999
        endif
!
        if( cmu.lt.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'cmu = ',cmu
          write(ifle,*) 'it must be >= 0'
          goto 9999
        endif
!
        if( sgmk.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'sgmk = ',sgmk
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
!
        if( sgme.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'sgme = ',sgme
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
!
        akslw (1)=klow
        akslw (2)=epslow
        sgmaks(1)=sgmk
        sgmaks(2)=sgme
      
      elseif(nrans.gt.0) then

      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputke
!
      end module module_rans

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_simple
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(19),parameter,private :: modnam='(module_simple)'
!
!
      integer :: nsmpl=1
      real*8  :: simeer=0.15D0
!
! nsmpl : no. of iteration counts in SIMPLE algorithm
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in):: cntlnam
      integer,intent(out) :: ierror
!
! [namelist]
      integer :: iter,PISO
      real*8  :: tolsimp
      namelist /simple/ iter,tolsimp,PISO
!
! [local entities]
      integer,parameter :: iundef=-huge(1)
      integer :: iset=0
      integer :: ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      iter=iundef
!
      rewind ifli
      read(ifli,simple,iostat=ios)
!      write(ifll,simple)
      call nml_errmsg0(ifle,ios,'simple',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if( iter.eq.iundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : iter'
        goto 9999
      elseif( iter.lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'iter =',iter
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      nsmpl=iter
      simeer=tolsimp
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_simple
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_source
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(19),parameter,private :: modnam='(module_source)'
!!!      include'interfacemmm'
!
!< module for source conditions & domains >
!
!-< 1. source conditions >-
!
!--< 1.1 flags >--
!
      integer,parameter :: kdpvol=1
      integer,parameter :: kdpmas=2
      integer,parameter :: kdmass=3
!
! kdpvol : value of kdscnd(:) for source term per unit volume
! kdpmas : value of kdscnd(:) for source term per unit mass
! kdmass : value of kdscnd(:) for sink/source of mass
!
!--< 1.2 values >--
!
      integer,private :: ncomp=0, nrans=0
      integer :: iscomp=0, israns=0
!
! ncomp  : no. of Ys(i) (chemical species)
! nrans  : no. of AKs(i) (dependent variables of RANS model)
! iscomp : start index of wdscnd for Ys(i)
! israns : start index of wdscnd for AKs(i)
!
      integer         :: nscnd =0
      integer,private :: nsvals=0
      integer,save,allocatable,private :: noscnd (:)
      integer,save,allocatable         :: kdscnd (:)
      integer,save,allocatable,private :: iwscnd (:)
      real*8 ,save,allocatable         :: wdscnd (:,:)
      real*8 ,save,allocatable,private :: wdtscnd(:,:)
!
! nscnd       : no. of source conditions
! nsvals      : 1st dimension size of "wdscnd"
! noscnd(:)   : domain no.
! iwscnd      : pointer for "wdtscnd"
! wdscnd(0,:) : cycle
! wdscnd(1,:) : r
! wdscnd(2,:) : u
! wdscnd(3,:) : v
! wdscnd(4,:) : w
! wdscnd(5,:) : T
! wdscnd      : Ys
! (6:5+ncomp,:)
! wdscnd      : AKs
! (6+ncomp:5+ncomp+nrans,:)
! wdtscnd     : time table for "wdscnd"
!
!
!-< 2. source domains >-
!
      integer,save,allocatable :: icscnd(:)
      integer,save,allocatable :: lcscnd(:)
!
! icscnd : pointer for "lcscnd"
! lcscnd : cell no. in source domain
!
!
!-< 99. for internal check >-
!
      integer,private :: iset1=0, iset2=0
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. check consistency of values in argument & module >---------------
!=================================================
      subroutine chkncomp(ival,sbnm)
!=================================================
      character(10),parameter :: subnam='(chkncomp)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.ncomp ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      stop 999
      end subroutine chkncomp
!
!< #2. check consistency of values in argument & module >---------------
!=================================================
      subroutine chknrans(ival,sbnm)
!=================================================
      character(10),parameter :: subnam='(chknrans)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.nrans ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      stop 999
      end subroutine chknrans
!
!< #3. interpolate values at present time making use of time table >----
!=================================================
      subroutine setnow(time)
!=================================================
      real*8,intent(in) :: time
      integer :: nwscnd
      if( nscnd.lt.1 ) return
      call modutl_chkset('>',1,iset1,modnam//'(setnow)')
      nwscnd=iwscnd(nscnd)
      call modutl_setnow(time,nscnd,nsvals,nwscnd,
     & iwscnd,wdtscnd,wdscnd)
      end subroutine setnow
!
!
!< #4. store domain data >----------------------------------------------
!=================================================
!      subroutine strdomain(ifle,mcell,nsdmn,nvsdmn,
!     & nosdmn,icsdmn,lcsdmn,ierr)
!=================================================
!      character(11),parameter :: subnam='(strdomain)'
!      integer,intent(in)  :: ifle,mcell,nsdmn,nvsdmn
!      integer,intent(in)  :: nosdmn(  mcell)
!      integer,intent(in)  :: icsdmn(0:mcell)
!      integer,intent(in)  :: lcsdmn(  mcell)
!      integer,intent(out) :: ierr
!      integer :: i,j,m,n
!      call modutl_chkset('>',1,iset1,modnam//subnam)
!      call modutl_chkset('=',1,iset2,modnam//subnam)
!      ierr=0
!      allocate( icscnd(0:nscnd),
!     &          lcscnd(nvsdmn), stat=ierr )
!      if( ierr.ne.0 ) then
!      write(ifle,*) '### error : allocation failed'
!      goto 9999
!      endif
!      icscnd(0)=0
!      i=0
!      do 100 n=1,nscnd
!      do 200 m=1,nsdmn
!      if( noscnd(n).eq.nosdmn(m) ) goto 201
!  200 continue
!      goto 100
!  201 continue
!      do 202 j=icsdmn(m-1)+1,icsdmn(m)
!      i=i+1
!      lcscnd(i)=lcsdmn(j)
!  202 continue
!      icscnd(n)=i
!  100 continue
!      return
! 9999 continue
!      write(ifle,*) modnam,subnam
!      ierr=1
!      end subroutine strdomain
!
!
!< #5. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,
     &                     ncmpx,nrnsx,lrans,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(in)   :: ncmpx,nrnsx
      logical,intent(in)   :: lrans
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      integer,parameter :: lc=10
      integer,parameter :: msrtim=10
      integer,parameter :: mcomp=200, mrans=50
      character(lc) :: kind
      integer :: no,ntime
      real    :: time(msrtim)
      real    :: r(msrtim),u(msrtim),v(msrtim),w(msrtim)
      real    :: t(msrtim),ys(mcomp*msrtim)
      real    :: aks(mrans*msrtim)
      real    :: cycle
      namelist /source/ no,kind,ntime,
     &                  time,r,u,v,w,t,ys,aks,
     &                  cycle
!
! --- [local entities]
!
      integer,parameter :: iundef=-huge(1)
      real   ,parameter :: undef=-huge(1.)
!
      character(lc),parameter :: mass='mass      ',
     &                      pervolume='pervolume ', 
     &                        permass='permass   '
!
      character(lc),parameter :: kndlst(3)=(/
     & mass, pervolume, permass /)
      real    :: dum1
      real*8  :: ysw(mcomp)
      integer :: i,j,n,i0,kt,ns,nwscnd,ncsdmn
      integer :: kd,kvel,ktmp,kcnc,krns
      integer :: ios,ierr1,ierr2,ierr3,ierr4
!
      call modutl_chkset('=',1,iset1,modnam//subnam)
      ierror=0
!
!
!-< 1. Initial set >-
!
      ncomp=ncmpx
      nrans=nrnsx
!
      call nml_chksiz(ifle,'mcomp',mcomp,ncomp,
     & modnam//subnam,ierr1)
      call nml_chksiz(ifle,'mrans',mrans,nrans,
     & modnam//subnam,ierr2)
      if( ierr1.ne.0 .or. ierr2.ne.0 ) goto 9999
!
!--< 1.1 count up sizes >--
!
      nscnd =0
      nwscnd=0
      rewind ifli
  100 continue
      ntime=1
      read(ifli,source,iostat=ios)
      if( ios.lt.0 ) goto 101
      nscnd=nscnd+1
      call nml_chksiz(ifle,'msrtim',msrtim,ntime,
     & modnam//subnam,ierr1)
      call nml_errmsg0(ifle,ios,'source',ierr2)
      if( ierr1.ne.0 .or. ierr2.ne.0 ) then
        write(ifle,*) 'sequence no. of the namelist =',nscnd
        goto 9999
      endif
      nwscnd=nwscnd+max(1,ntime)
      goto 100
  101 continue
!
!--< 1.2 allocate arrays >--
!
      iscomp=5
      israns=5+ncomp
      nsvals=israns+nrans
      allocate( noscnd(nscnd),
     &          kdscnd(nscnd),
     &          iwscnd(0:nscnd),
     &          wdscnd(0:nsvals,nscnd),
     &          wdtscnd(0:nsvals,nwscnd), stat=ierr1 )
      if( ierr1.ne.0 ) then
      write(ifle,*) '### error : allocation faild'
      goto 9999
      endif
      noscnd =0
      kdscnd =0
      iwscnd =0
      wdscnd =0.d0
      wdtscnd=0.d0
!
!
!-< 2. Input namelist >-
!
!--< 2.1 read namelist >--
!
      kind=pervolume
!
      ns=0
      kt=0
      rewind ifli
 1000 continue
      ntime = 1
      time  = 0.
      cycle = undef
      r     = undef
      u     = undef
      v     = undef
      w     = undef
      t     = undef
      ys    = undef
      aks   = undef
      read(ifli,source,iostat=ios)
      if( ios.lt.0 ) goto 1001
!      write(ifll,source)
      call nml_errmsg0(ifle,ios,'source',ierr1)
      if( ierr1.ne.0 ) goto 9999
      ns=ns+1
!
!--< 2.2 common procedure >--
!
!/ domain no. /
!
      if( no.eq.iundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : no'
        goto 9999
      endif
      do 200 n=1,ns-1
      if( no.eq.noscnd(n) ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) 'no =',no
      write(ifle,*) 'this no. already specified in previous line'
      write(ifle,*) 'sequence no. of the line =',n
      goto 9999
      endif
  200 continue
      noscnd(ns)=no
!
!/ global kind /
!
      call nml_listno(3,kndlst,kind,kd)
      call nml_chkchr0(ifle,'kind',kind,kd,ierr1)
      if( ierr1.ne.0 ) goto 9999
      kind=kndlst(kd)
!
!/ no. of times in time table /
!
      if( ntime.lt.1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'ntime =',ntime
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
      do 210 i=2,ntime
      if( time(i-1).ge.time(i) ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) 'time(i-1),time(i) = ',time(i-1),time(i)
      write(ifle,*) 'time(i-1) must be < time(i)'
      write(ifle,*) 'i =',i
      goto 9999
      endif
  210 continue
      dum1=time(ntime)-time(1)
      if( cycle.ne.undef .and. cycle.le.dum1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'cycle = ',cycle
        write(ifle,*) 'it must be > time(ntime)-time(1)'
        write(ifle,*) 'time(ntime)-time(1) = ',dum1
        goto 9999
      endif
      kt=kt+ntime
      iwscnd(ns)=kt
!
!/ check data /
!
      if( .not.lrans ) aks=0.
!
      if( kind.ne.mass ) then
        if( r(1).ne.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'r(1) can not be specified'
        goto 9999
        endif
        r=0.
        call chkflt0(1,1,ntime,'u',u)
        call chkflt0(1,1,ntime,'v',v)
        call chkflt0(1,1,ntime,'w',w)
        call chkflt0(1,1,ntime,'t',t)
        call chkflt0(2,ncomp,ntime,'ys',ys)
        call chkflt0(2,nrans,ntime,'aks',aks)
                              kdscnd(ns)=kdpvol
        if( kind.eq.permass ) kdscnd(ns)=kdpmas
      else
        call chkflt1(0,ntime,'r',r)
        call chkflt1(0,ntime,'u',u)
        call chkflt1(0,ntime,'v',v)
        call chkflt1(0,ntime,'w',w)
        call chkflt1(1,ntime,'t',t)
        call chkflt2(ncomp,ntime,'ys',ys)
        call chkflt2(nrans,ntime,'aks',aks)
        if( ierror.ne.0 ) goto 9999
        kdscnd(ns)=kdmass
      endif
!
!/ store data /
!
                           wdscnd(0,ns)=-1.d0
      if( cycle.ne.undef ) wdscnd(0,ns)=dble(cycle)
      i0=iwscnd(ns-1)
      do 220 i=1,ntime
      wdtscnd(0,i+i0)=dble(time(i))
      wdtscnd(1,i+i0)=dble(r(i))
      wdtscnd(2,i+i0)=dble(u(i))
      wdtscnd(3,i+i0)=dble(v(i))
      wdtscnd(4,i+i0)=dble(w(i))
      wdtscnd(5,i+i0)=dble(t(i))
      if( kind.eq.mass ) then
      call setys(i,ys,ysw,ierr1)
      if( ierr1.ne.0 ) goto 9999
      else
      ysw(1:ncomp)=dble(ys(ncomp*(i-1)+1:ncomp*i))
      endif
      do 221 j=1,ncomp
      wdtscnd(iscomp+j,i+i0)=ysw(j)
  221 continue
      do 222 j=1,nrans
      wdtscnd(israns+j,i+i0)=aks(j)
  222 continue
  220 continue
!
      goto 1000
!
 1001 continue
!
!
!-< 3. Initialize "wdscnd" >-
!
      do 300 n=1,nscnd
      kt=iwscnd(n-1)+1
      do 301 i=1,nsvals
      wdscnd(i,n)=wdtscnd(i,kt)
  301 continue
  300 continue
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      contains
!
!< #5.1 check real data >
!=================================================
      subroutine chkflt0(idm,nc,nt,dnam,dval)
!=================================================
      integer     ,intent(in)    :: idm,nc,nt
      character(*),intent(in)    :: dnam
      real        ,intent(inout) :: dval(nc,nt)
      integer :: i,j
      if( nc.lt.1 ) return
      if( dval(1,1).eq.undef ) then
        dval=0.
      else
        do 100 i=1,nc
        do 100 j=1,nt
        if( dval(i,j).eq.undef ) then
        write(ifle,*) '### error : data error'
        if( idm.gt.1 ) then
        write(ifle,*) 'lack of data : ',dnam,'(i,j)'
        write(ifle,*) 'i,j =',i,j
        else
        write(ifle,*) 'lack of data : ',dnam,'(i)'
        write(ifle,*) 'i =',j
        endif
        goto 9999
        endif
  100   continue
      endif
      return
 9999 continue
      ierror=1
      end subroutine chkflt0
!
!< #5.1 check real data >
!=================================================
      subroutine chkflt1(isw,nt,dnam,dval)
!=================================================
      integer     ,intent(in) :: isw,nt
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval(nt)
      integer :: i
      do 100 i=1,nt
      if( dval(1).eq.undef ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) 'lack of data : ',dnam,'(i)'
      write(ifle,*) 'i =',i
      goto 9999
      endif
  100 continue
      if( isw.eq.0 ) return
      do 101 i=1,nt
      if( ( isw.eq.1 .and. dval(i).le.0.d0 ) .or. 
     &    ( isw.eq.2 .and. dval(i).lt.0.d0 ) ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) dnam,'(i) = ',dval(i)
      write(ifle,*) 'i =',i
      if( isw.eq.1 ) write(ifle,*) 'it must be > 0'
      if( isw.eq.2 ) write(ifle,*) 'it must be >= 0'
      goto 9999
      endif
  101 continue
      return
 9999 continue
      ierror=1
      end subroutine chkflt1
!
!< #5.2 check real data >
!=================================================
      subroutine chkflt2(nc,nt,dnam,dval)
!=================================================
      integer     ,intent(in) :: nc,nt
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval(nc,nt)
      integer :: i,j
      do 100 i=1,nc
      do 100 j=1,nt
      if( dval(1,1).eq.undef ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) 'lack of data : ',dnam,'(i,j)'
      write(ifle,*) 'i,j =',i,j
      goto 9999
      endif
  100 continue
      do 101 i=1,nc
      do 101 j=1,nt
      if( dval(i,j).lt.0.d0 ) then
      write(ifle,*) '### error : data error'
      write(ifle,*) dnam,'(i,j) = ',dval(i,j)
      write(ifle,*) 'i,j =',i,j
      write(ifle,*) 'it must be >= 0'
      goto 9999
      endif
  101 continue
      return
 9999 continue
      ierror=1
      end subroutine chkflt2
!
!< #5.5 set mass fraction >
!=================================================
      subroutine setys(i,ys,ysw,ierr)
!=================================================
      integer,intent(in)  :: i
      real   ,intent(in)  :: ys (mcomp*msrtim)
      real*8 ,intent(out) :: ysw(mcomp)
      integer,intent(out) :: ierr
      integer :: j,i0
      real*8  :: sum
      i0=ncomp*(i-1)
      sum=0.d0
      do 100 j=1,ncomp
      ysw(j)=dble(max(0.,ys(j+i0)))
      sum=sum+ysw(j)
  100 continue
      if( sum.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'sum of ys(i-j) = ',sum
        write(ifle,*) 'it must be > 0'
        write(ifle,*) 'i,j =',i0+1,i0+ncomp
        ierr=1
        return
      endif
      sum=1.d0/sum
      do 101 j=1,ncomp
      ysw(j)=ysw(j)*sum
  101 continue
      end subroutine setys
!
      end subroutine inputdata
!
      end module module_source

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_time
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(17),parameter,private :: modnam='(module_time)'
!
!< module for start/end time >
!
      integer :: iters=0, itere=0
      logical :: steady
      
      real*8  :: timee=0.d0,toltim=0.001
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in):: cntlnam
      integer,intent(out) :: ierror
!
! [namelist]
      integer :: start,end,flowcon,MHD_cond
      real*8  :: end_time,toltime
      namelist /time/ start,end,end_time,flowcon,toltime,MHD_cond
!
! [local entities]
      integer,parameter :: iundef=-huge(1)
      integer :: iset=0
      integer :: ios,ierr1
!flowcon
!       =1: steady flow
!       =2: transient flow
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      start   =iundef
      end     =iundef
      end_time=huge(1.d0)
!
      rewind ifli
      read(ifli,time,iostat=ios)
!      write(ifll,time)
      call nml_errmsg0(ifle,ios,'time',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if( start.eq.iundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : start'
        goto 9999
      elseif( start.lt.-1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'start =',start
        write(ifle,*) 'it must be >= -1'
        goto 9999
      endif
!
      if( end.eq.iundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : end'
        goto 9999
      elseif( end.lt.start ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'end =',end
        write(ifle,*) 'it must be >= start'
        goto 9999
      endif
!
      if( end_time.lt.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'end_time = ',end_time
        write(ifle,*) 'it must be >= 0'
        goto 9999
      endif
!
      iters=start
      itere=end
      timee=end_time
!
      if(flowcon.eq.1.or.flowcon==3) then 
        steady=.true.
        if(toltime.gt.0.d0) then
          if(toltime.gt.1.d0) then
            write(ifle,*) '### error : data error'
            write(ifle,*) 'toltime must be smaller than one'
            stop
          else
            toltim=toltime
          endif
        else
          write(ifle,*) '### error : data error'
          write(ifle,*) 'toltime must be great than zero'
          stop
        endif
      elseif(flowcon.eq.2) then
        steady=.false.
        toltim=0.d0
      else
        write(ifle,*) '### error : data error'
        write(ifle,*) 'flowcon must be 1 or 2'
        stop
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_time
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_turbparm
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!< module for parameters of turbulent flow >
!
      real*8 :: prturb=0.9d0
      real*8 :: scturb=0.9d0
!
! prturb : turbulent Prandtl no.
! scturb : turbulent Schmidt no.
!
      real*8 :: akappa=0.42d0
      real*8 :: awfnc =5.5d0
!
! akappa : Karman constant
! awfnc  : constant of log-law profile near the wall
!
      end module module_turbparm
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_species
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_material,only : isthr,nflud
      implicit none
      character(20),parameter,private :: modnam='(module_species)'
!
!< module for species data >
!
!
!-< 0. Constants >-
!
      real*8,parameter :: gascns=8.314d0
 
!
! gascns : unversal gas constant [J/mol/K]
!
!
!-< 1. Name of species >-
!
      integer,private :: ns=0,ns_gas,ns_suf
!
      integer,parameter :: lenspc=20
      character(lenspc),save,allocatable :: spcnam(:)
!
!-< 2. Thermodynamic constant >-
!
      real*8,save,allocatable :: wm(:)
      real*8,save,allocatable :: r_wm(:)
      real*8,save,allocatable :: Ri(:)
      real*8,save,allocatable :: t3(:,:)
      real*8,save,allocatable :: a7(:,:,:)
      real*8,save,allocatable :: acpk(:,:)
      real*8,save,allocatable :: hform(:),href(:),tref_comp(:)
      real*8,save,allocatable :: Le(:)
      real*8,save,allocatable :: r_Le(:)
      real*8 :: hreff,tref=298.15d0
      logical,save,allocatable:: judg_t2(:)
      logical,save :: sw
!
!wm(:)     : molecular weight [kg/mol]
!acpk(:,:) : acpk(0) : heat of formation [J/kg]
!acpk(1)   : specific heat at constant pressure [J/kg/K]
!          : 4th-order polynomial of temperature
!r_wm(:)   : reciprocal of 'wm' [mol/kg]
!Ri(:)     : 'gascns'/'wm' [J/(kg K)]
!t3(:,:)   : point of temperature region [K]
!a7(:,:,:) : coefficient sets to specific heat at
!            constant pressure [J/(kg K)],
!            enthalpy [J/kg] & entropy [J/(kg K)]
!Le(:)     : Lewis Number of species [-]
!r_Le(:)   : reciprocal Lewis Number of species [-]
!judg_t2(:): judgment "t3(2,*)". every value is same or NOT
!sw        : switch of coefficient sets : T=> a7, F=> acpk
!
!
!-< 3. eular two phase if [ieul2ph=1] >
!
      real*8,save,allocatable :: wm2(:),acpk2(:,:)
!
!------------------------
!-< 4. surface reaction 
!------------------------
      character(lenspc),save,allocatable :: spcnamsf(:)
!-< 99. for internal check >-
!
      integer,private :: iset=0,icom
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. check consistency of values in argument & module >---------------
!=================================================
      subroutine chkncomp(ival,sbnm)
!=================================================
      implicit none
      character(10),parameter :: subnam='(chkncomp)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.ns_gas ) return
      call modutl_chkset('>',1,iset,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      stop ': at module_species/chkncomp'
      end subroutine chkncomp
!
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata_gas
     & (ifli,ifll,ifle,cntlnam,ncomp,E2P,lvof,ierror)
!=================================================
      implicit none
!
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle,ncomp
      character(*),intent(in) :: cntlnam
      logical,intent(in)   :: E2P,lvof
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      character(lenspc) :: name         ! species name
      real*8            :: weight       ! molecular weight [kg/mol]
      real*8            :: a(7,2)       ! coefficient set
      real*8            :: t(3)
                                        ! temperature [K] (1:upper limit,
                                        ! 2:branch point, 3:lower limit)
      real*8            :: cp(5)        ! coefficient set
      real*8            :: ho           ! enthalpy [J/kg]
      real*8            :: Lewis        ! Lewis number [-]
      integer           :: flag_a7=0,fa7,bg_gas
      real*8            :: weight2,cp2(5),ho2,tref2,diffusity,pref,
     &                     Eps_K,sigma
      integer            :: active=1
      character(200)     :: th_file,sp_file	! input file name
      character(200)     :: thfile,spfile
      integer,parameter  :: mcomp=200           ! max number of species
      integer            :: gasno		!
      real(8)            :: T_c,P_c,rho_c,K_ij(mcomp),OMEGA,dipole,kappa
      namelist /species/ name,weight,
     &                   cp, ho, tref,pref,
     &                   weight2,cp2,ho2,tref2,
     &                   flag_a7,a,t,Lewis,
     &                   th_file,sp_file,bg_gas,diffusity,
     &                   gasno,active,
     &                   Eps_K,sigma,
     &                   T_c,P_c,OMEGA,K_ij,Rho_c,dipole,kappa
!
!
! --- [local entities]
!
      real*8,parameter :: undef=-huge(1.d0)
      real*8  :: hk(5),tt,h00
      integer :: i,ios=0,iosf=0,ierr1,s,s2,nsno
      logical,external :: nml_comp_eq,nml_comp_gt
      logical :: Lel
!
      ns = ncomp
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
      ierr1=0
!
!-< 1. allocate array >-
!
      allocate( spcnam(ns),      ! name of species
     &          wm    (ns),      ! molecular weight [kg/mol]
     &          r_wm  (ns),      ! reciprocal of 'wm' [mol/kg]
     &          Ri    (ns),      ! 'gascns'/'wm' [J/(kg K)]
     &          wm2   (    ns),
     &          acpk2 (0:5,ns),
     &          stat=ierr1 )
      if(ierr1/=0) goto 9999
!------------------------------------------------
! --- judgment of coefficient set "a7" or "acpk"
!------------------------------------------------
 1002 continue
      rewind ifli
      nsno=0
      do
        read(ifli,species,iostat=ios)
        if(ios<0) then
          exit
        elseif(ios==0) then
          nsno=nsno+1
          if(nsno==1) then
            fa7=flag_a7
          endif
        else
          write(ifle,species)
          stop' Reading error in [&species]'
        endif
      enddo
      if(nsno.ne.ncomp) then
        write(ifle,*) '### ERR: ncomp /= nsno'
        write(ifle,*) '### MSG: ncomp= ',ncomp,' in ',cntlnam
        write(ifle,*) '### MSG: [&species] number : ',nsno
        write(ifle,*) '### reading error in "&species"'
        ierror = 1
        stop' Reading error in [&species]'
      endif
      ns_gas=nsno
      sw=.true.                        ! for "a7"
      if(fa7.eq.0)  sw=.false.         ! for "acpk"
!----------------------------------------------------
! --- reading species name, molecular weight [kg/mol]
!----------------------------------------------------
      spcnam=''
      rewind ifli
      do s=1,ns
        name=' '                      ! name of species
        read(ifli,species,iostat=ios)
        write(ifll,'(1x,a,I4,2X,a)') 
     &  'MSG: Species Number and Name: ',s,trim(name)
        if(name==' ') then            ! check "name", name of species
          write(ifle,*) '### reading error in "species"'
          write(ifle,*) 'lack of data "name" in no.',s
          ierror=1
        endif
        spcnam(s)=name
        if(s>=2) then
          do s2=1,s-1
            if(spcnam(s)==spcnam(s2)) then
              write(ifle,*) '### reading error in "species"'
              write(ifle,*) 'duplicated data : name(i)=name(j)=',
     &                       spcnam(s2)
              write(ifle,*) 'i,j =',s2,s
              ierror=1
            endif
          enddo
        endif
      enddo
!
      if(s-1/=ns) then                ! check number of species
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of namelists (species) =',s-1
        write(ifle,*) 'it must be =',ns
        ierror=1
      endif
!
!--------------------------------------------------------
! --- reading coefficient sets and display reading result
!--------------------------------------------------------
      if(sw) then
        allocate( a7(7,2,ns),    ! coefficient sets
     &            t3(  3,ns),    ! temperature [K]
     &            judg_t2(1),    ! judgment "t3(2,:)"
                                 ! . every value is same or NOT
     &            hform(ns),href(ns),tref_comp(ns),
     &            stat=ierr1)
!
! --- read again
!
        rewind ifli
        do s=1,ns
          read(ifli,species,iostat=ios)
        enddo
      else
!--------------
! --- a5 model
!--------------
        rewind ifli
        do s=1,ns
          read(ifli,species,iostat=ios)
        enddo
      endif
!----------------------
! --- final error check
!----------------------
      if(ierror/=0) then
        write(ifle,*) "total error : ",ierror,modnam,subnam
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata_gas
!
!=================================================
      subroutine inputdata_suf
     & (ifli,ifll,ifle,cntlnam,ncomp_suf,ierror)
!=================================================
      implicit none
!
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifli,ifll,ifle,ncomp_suf
      character(*),intent(in) :: cntlnam
      integer,intent(out)     :: ierror
!
! --- [namelist]
!
      character(lenspc) :: name         ! species name
      real*8            :: weight       ! molecular weight [kg/mol]
      real*8            :: a(7,2)       ! coefficient set
      real*8            :: t(3)
                                        ! temperature [K] (1:upper limit
                                        ! 2:branch point, 3:lower limit)
      real*8            :: ho           ! enthalpy [J/kg]
      real*8            :: Lewis        ! Lewis number [-]
      integer           :: flag_a7=0,fa7
      character(200)    :: th_file,sp_file ! input file name
      namelist /surface_species/ name,weight,tref,
     &                           flag_a7,a,t,Lewis
     &                           ,th_file,sp_file
!   
! --- [local entities]
!
      real*8,parameter :: undef=-huge(1.d0)
      real*8  :: hk(5),tt,h00
      integer :: i,ios=0,iosf=0,ierr1,s,s2,nsno
      logical,external :: nml_comp_eq,nml_comp_gt
      logical :: Lel
!
      iset=0
      ns = ncomp_suf
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
      ierr1=0
!
!-< 1. allocate array >-
!
      write(ifll,*)
      allocate(spcnamsf(ns),      ! name of species
     &         stat=ierr1 )
!------------------------------------------------
! --- judgment of coefficient set "a7" or "acpk"
!------------------------------------------------
      rewind ifli
      nsno=0
      do
        name=''
        read(ifli,surface_species,iostat=ios)
        if(ios<0) then
          exit
        elseif(ios==0) then
          nsno=nsno+1
          if(nsno==1) fa7=flag_a7
        else
          write(ifle,surface_species)
          stop' Reading error in [&surface_species]'
        endif
      enddo
      
      ns_suf=nsno
      if(ns_suf==0.and.ncomp_suf==0) return
      if(nsno/=ncomp_suf) then
        write(ifle,*) '### ERR: ncomp_suf /= nsno'
        write(ifle,*) '### MSG: ncomp_suf= ',ncomp_suf,' in ',cntlnam
        write(ifle,*) '### MSG: [&surface_species] number : ',nsno
        write(ifle,*) '### reading error in [&surface_species]'
        ierror=1
        stop ' Reading error in [&surface_species]'
      endif
      sw=.true.                        ! for "a7"
      if(fa7/=1) stop 'ERR:flag_a7 must be 1 for surface species'
!------------------------------------------------------
! --- reading species name, molecular weight [kg/mol]
!------------------------------------------------------
      spcnamsf=''
      rewind ifli
      do s=1,ns
        name=' '                       ! name of species
        read(ifli,surface_species,iostat=ios)
        write(ifll,'(1x,a,I4,2X,a)') 
     &  'MSG: Species Number and Name: ',s,trim(name)
        if(name==' ') then            ! check "name", name of species
          write(ifle,*) '### reading error in "species"'
          write(ifle,*) 'lack of data "name" in no.',s
          ierror=1
        endif
        spcnamsf(s)=name
        if(s>=2) then
          do s2=1,s-1
            if(spcnamsf(s)==spcnamsf(s2) ) then
              write(ifle,*) '### reading error in "species"'
              write(ifle,*) 'duplicated data : name(i)=name(j)=',
     &                       spcnamsf(s2)
              write(ifle,*) 'i,j =',s2,s
              ierror=1
            endif
          enddo
        endif
      enddo
!
      if( s-1/=ns ) then                  ! check number of species
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of namelists (species) =',s-1
        write(ifle,*) 'it must be =',ns
        ierror = 1
      endif
!
!----------------------
! --- final error check
!----------------------
      if( ierror/=0 ) then
        write(ifle,*) "total error : ",ierror,modnam,subnam
        stop
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      return
      end subroutine inputdata_suf
!
!=================================================
      subroutine str_retrieve(str,p,p_end,str_r,w)
!=================================================
      integer, intent(inout) :: p		! start position for search
      integer, intent(in) :: p_end		! end position of str (>p)
      character(p_end), intent(in) :: str	! string data
      character(p_end), intent(out) :: str_r	! retrieve data
      integer, intent(out) :: w			! width of str_r
      integer i,j,p_range(2)
      do i=p,p_end
        if(str(i:i)/=" ") then
          p_range(1) = i-1	! left edge position of retrieve data
          exit
        endif
      enddo
      do j=i,p_end
        if(str(j:j)==" ") then
          p_range(2) = j-1	! right edge position of retrieve data
          exit
        endif
      enddo
      str_r = str(p_range(1):p_range(2))
      p = p_range(2)+1
      w = p_range(2)-p_range(1)+1
      end subroutine str_retrieve
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!------------------------------------------------------------------------
! --- judgment of temperature region ( 1:higher region, 2:lower region )
!=================================================
      integer function judg_temp(ICOM,t)
!=================================================
      implicit none
        integer,intent(in) :: ICOM  ! classification number of species
        real*8, intent(in) :: t     ! temperature [K]
        if( t>t3(2,ICOM) ) then     ! split by temperature
          judg_temp = 1             ! higher
        else
          judg_temp = 2             ! lower
        endif
      end function judg_temp
!
!=================================================
      real*8 function c_pi(t,ICOM)  ! specific heat of species [J/(kg K)]
!=================================================
      implicit none
        integer,intent(in) :: ICOM  ! classification number of species
        real*8, intent(in) :: t     ! temperature [K]
        integer rj                  ! result of judgment
        rj = judg_temp(ICOM,t)
        c_pi =a7(1,rj,ICOM)
     &    +t*(a7(2,rj,ICOM)
     &    +t*(a7(3,rj,ICOM)
     &    +t*(a7(4,rj,ICOM)
     &    +t* a7(5,rj,ICOM) )))
        c_pi =c_pi*Ri(ICOM)
      end function c_pi
!
!=================================================
      real*8 function c_p(ncomp,t,yi) ! specific heat of mixture gas
                                          ! [J/(kg K)]
!=================================================
      implicit none
        integer,intent(in) :: ncomp       ! number of species
        real*8, intent(in) :: t,yi(ncomp) ! temperature [K],
                                          ! mass fraction of species
        integer icom                      ! classification number of species
        c_p=0
        do icom=1,ncomp
          c_p=c_p+yi(icom)*c_pi(t,icom)
        enddo
      end function c_p
!
!=================================================
      real*8 function h_t(icom,t)       ! (enthalpy of species)/(universal
                                        ! gas constant) [K]
!=================================================
      implicit none
        integer,intent(in) :: icom      ! classification number of species
        real*8, intent(in) :: t         ! temperature [K]
        integer rj                      ! result of judgment
        rj = judg_temp(icom,t)
        h_t =a7(6,rj,icom)
     &   +t*(a7(1,rj,icom)
     &   +t*(a7(2,rj,icom)*0.5d0
     &   +t*(a7(3,rj,icom)*1.d0/3.d0
     &   +t*(a7(4,rj,icom)*0.25d0
     &   +t*(a7(5,rj,icom)*0.2d0)))))
      end function h_t
!
!=================================================
      real*8 function enthalpy(ncomp,t,yi) ! enthalpy of mixture gas [J/kg]
!=================================================
        implicit none
        integer,intent(in) :: ncomp        ! number of species
        real*8, intent(in) :: t,yi(ncomp)  ! temperature [K],
                                           ! mass fraction of species
        integer icom                       ! classification number of species
        enthalpy=0.d0
        do icom=1,ncomp
          enthalpy=enthalpy
     &            +yi(icom)*(h_t(icom,t)*Ri(icom)-href(icom))
        enddo
      end function enthalpy
!
!=================================================
      real*8 function s_0(icom,t)       ! (entropy of species)/
                                        ! (universal gas constant) [-]
!=================================================
      implicit none
        integer,intent(in) :: icom      ! classification number of species
        real*8, intent(in) :: t         ! temperature [K]
        integer rj                      ! result of judgment
        rj = judg_temp(icom,t)
        s_0 = a7(7,rj,icom)
     &   +t*( a7(2,rj,icom)
     &   +t*( a7(3,rj,icom)*0.5d0
     &   +t*( a7(4,rj,icom)*0.33333d0
     &   +t*( a7(5,rj,icom)*0.25d0 ))))
     &   +    a7(1,rj,icom)*log(t)
      end function s_0
!
!=================================================
      real*8 function entropy(ncomp,t,yi) ! entropy of mixture gas [J/(kg K)]
!=================================================
      implicit none
        integer,intent(in) :: ncomp       ! number of species
        real*8, intent(in) :: t,yi(ncomp) ! temperature [K],
                                          ! mass fraction of species
        integer icom                      ! classification number of species
        entropy = 0.d0
        do icom=1,ncomp
          entropy = entropy+yi(icom)*s_0(icom,t)*Ri(icom)
        enddo
      end function entropy
!
!=================================================
      real*8 function sigma_yr(ncomp,yi)   ! [J/(kg K)]
!=================================================
      implicit none
        integer,intent(in) :: ncomp        ! number of species
        real*8, intent(in) :: yi(ncomp)    ! mass fraction of species
        integer icom                       ! classification number of species
        sigma_yr = 0
        do icom=1,ncomp
          sigma_yr = sigma_yr+yi(icom)*Ri(icom)
        enddo
      end function sigma_yr
!
      end module module_species
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_chemreac
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_species,only : a7
      implicit none
      character(17),parameter,private :: modnam='(module_chemreac)'
!
! --- < module for chemical reaction formula >
!
      integer,parameter :: igas=1,isurface=2
      integer,private   :: ns=0       ! number of species
      integer,private   :: ns_sf=0    ! number of species on SITE
      integer,parameter :: ovrall=1,  ! ovrall : overall reaction model
     &                     elemnt=2,  ! elemnt : elementary reaction model
     &                     ebrkup=3,  ! ebrkup : eddy break up model
     &                     ovalfr=4,  ! ovrall : overall reaction model + CO and SOOT
     &                     userreac=5,! userreac:user reaction model
     &                     stick_sf=6,! surface stick model
     &                     LH_ER=7,   ! surface
     &                     Bohm=8,    ! surface
     &                     elemarbi=9,! surface & Gas
     &                     BohmE=10,  ! surface by Ion-energy-dependent Bohm model
     &                     BohmY=11,   ! surface by Ion-enhanced Reac. Yield Bohm Model
     &                     Butler_Volmer=12,
     &                     Zeldovich=13,
     &                     oval_ebrkup=14,
     &                     SSFRRM=15

      integer :: nq=0,nneq            ! number of equation
      integer,save,allocatable:: ireq(:)      ! flag of reaction model
      real*8,save,allocatable :: vreq(:,:,:)  ! stoichiometric coefficients
      real*8,save,allocatable :: preq(:,:,:)  ! parameters for reaction rate
      real*8,save,allocatable :: mreq(:,:)    ! coefficients of 3rd body
      real*8,save,allocatable :: creq(:)      ! index number of mol
                                              ! concentration(overall)
      real*8,save,allocatable :: sub_vreq(:,:)! subtraction of "vreq" of each
                                              ! species in each equation
      real*8,save,allocatable :: sigma_vreq(:)! change of "sub_vreq" in
                                              ! each equation
      real*8,save,allocatable ::vreq_a7(:,:,:)! summation of "vreq" * "a7"
                                              ! in each equation
      integer,private :: iset=0
!
!///////////////////////////////////////////////////////////////////////
      contains
!---------------------------------------------------------
!< #1. check consistency of values in argument & module >-
!---------------------------------------------------------
!====================================
      subroutine chkncomp(ival,sbnm)
!====================================
      character(10),parameter :: subnam='(chkncomp)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
        if( ival==ns ) return
        call modutl_chkset('>',1,iset,modnam//subnam)
        write(*,*) '### program error -1- ',modnam//subnam
        write(*,*) 'sbnm=',sbnm
        stop
      end subroutine chkncomp
!----------------------------
!< #2. input name list    >--
!----------------------------
!=========================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,lsuf,
     &                     ncmpx,ncomp_sufx,lcomp,ierror)
!=========================================================
      use module_species,only : lenspc,spcnam,chk_spec=>chkncomp
      use module_species,only : spcnamsf
!
      implicit none
!
      character(11),parameter :: subnam='inputdata'
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle,ncmpx,ncomp_sufx
      character(*),intent(in) :: cntlnam
      logical,intent(in)   :: lcomp
      logical,intent(out)  :: lsuf
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      integer,parameter :: lc=8
      integer,parameter :: mcomp=200            ! number of species
      character(20)     :: model
      character(lenspc) :: nleft(mcomp+1)       ! name of species on left side
      character(lenspc) :: nright(mcomp+1)      ! name of species on right side
      real*8 :: cleft(mcomp+1),cright(mcomp+1) ! stoichiometric coefficient
      real*8  :: creac(mcomp+1)                 ! index number of mol
      real*8  :: pleft(mcomp+1),pright(mcomp+1) ! arbitrary power index
                                                ! concentration(overall)
      real*8  :: a,alpha,beta,e
      real*8  :: af,alphaf,ef,ab,alphab,eb      ! k = (af)*T^(alpahf)*
                                                ! EXP(-(ef)/(R*T))
      real*8  :: M(mcomp)                       ! coefficients of 3rd body
      real*8  :: cr1,cr2,r
      real*8,parameter :: WM_CO=28.d-3,WM_C=12.d-3
      real*8  :: Y_CO,nu_SOOT,nu_CO
      real*8  :: Y_SOOT,C_fuel,WM_fuel,nu_fuel
      real(8) :: a_LH_ER(mcomp),beta_LH_ER(mcomp),H_LH_ER(mcomp),
     &           L_LH_ER(mcomp),n_LH_ER(mcomp),m_LH_ER,k_LH_ER
!
      real(8) :: coeff_Bohm
      integer :: ignition_iter=0,stop_iter
      character(lc) :: kind,prs_kind
      integer :: unit_A_e
      real*8  :: af_low,alphaf_low,ef_low,
     &           ab_low,alphab_low,eb_low        ! for Lindemann etc.
      real*8  :: troe_a(2),troe_t1(2),troe_t2(2),troe_t3(2) ! Troe form
      real(8) :: I_ref_BV,alp_a_BV,alp_c_BV,
     &           r_BV_left(mcomp+1),r_BV_right(mcomp+1),
     &           ref_mol_frc_left(mcomp+1),
     &           ref_mol_frc_right(mcomp+1)
      namelist /chemreac/ model,nleft,nright,cleft,cright,
     &                                       pleft,pright,
     &                    creac,a,alpha,beta,e,
     &                    Y_SOOT,C_fuel,nu_fuel,WM_fuel,
     &                    af,alphaf,ef,
     &                    ab,alphab,eb,M,
     &                    cr1,cr2,r,ignition_iter,kind,
     &                    a_LH_ER,beta_LH_ER,H_LH_ER,L_LH_ER,
     &                    n_LH_ER,m_LH_ER,k_LH_ER,
     &                    coeff_Bohm,unit_A_e,
     &                    af_low,alphaf_low,ef_low,
     &                    ab_low,alphab_low,eb_low,
     &                    troe_a,troe_t1,troe_t2,troe_t3,
     &                    prs_kind,
     &                    I_ref_BV,
     &                    alp_a_BV,alp_c_BV,
     &                    r_BV_left,r_BV_right,
     &                    ref_mol_frc_left,ref_mol_frc_right,
     &   stop_iter
!
! --- [local entities]
!
      integer,parameter :: iundef=-huge(1)
      real*8 ,parameter :: undef=-huge(1.d0)
      character(LEN=11),parameter :: mdllst(15)=(/
     &             'overall    ','elementary ','eddybreakup',
     &             'fire       ','user       ','stick      ',
     &             'LH_ER      ','BOHM       ','arbi_elem  ',
     &             'E_BOHM     ','Y_BOHM     ','BV         ',
     &             'Zeldovich  ','eddy_over  ','SSFRRM     '/)
!
      character(lc),parameter :: 
     &    chemkd(2)=(/'gas     ','surface '/)
      character(lc),parameter :: gas='gas',surface='surface'
!
      integer :: ispc(mcomp)
      integer :: i,j,k,ios,imodl,s
      integer :: q,ierr1,kd
      real*8  :: dd2
      integer nol,nor              ! number of species in each side
      logical,external :: nml_comp_eq,nml_comp_gt,nml_comp_ge
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror = 0
!
!-< 1. Initial set >-
      nq=0                          ! count up equation
      nneq=0
      rewind ifli
      do
        read(ifli,chemreac,iostat=ios)
        if(ios<0) exit
        nq = nq+1
        call nml_errmsg0(ifle,ios,'chemreac',ierr1)
        if( ierr1/=0 ) then
          write(ifle,*) ' ### ERR: sequence no. of the namelist =',nq
          ierror = 1+ierror
        endif
      enddo
      nneq=nq
      if(nneq.gt.1) then
        if(.not.lcomp) then
          write(ifle,*) ' ### Warning: Chemical reaction must use ',
     &                     'compressible or low-Mach approximation'
        endif
      endif
      if(ncomp_sufx>0.and.nq<1) then
          write(ifle,'(2a)') 'ERR: NO surface reaction defined, ',
     &    'reset [&sizes/ncomp_surface=0]'
          stop
        endif
      if(nq<1) return                ! no reaction flow
      if(ierror.ne.0) goto 9999
!---------------------------
!--< 1.2 check interface >--
!---------------------------
      ns=ncmpx
      ns_sf=ncomp_sufx
      call nml_chksiz(ifle,'mcomp',mcomp,ns,    ! if ns>mcomp, ierr1=1
     &                      modnam//subnam,ierr1)
      if(ierr1/=0) then
        ierror=1
        goto 9999
      endif
      call chk_spec(ns,modnam//subnam)
!---------------------------
!--< 1.3 allocate arrays >--
!---------------------------
      allocate( ireq(           nq),          ! flag of reaction model
     &          vreq(ns+ns_sf,2,nq),          ! stoichiometric coefficients
     &          preq( 6,2,      nq),          ! parameters for reaction rate
     &          creq(ns+ns_sf     ),          ! index number of mol
                                              ! concentration(overall)
     &          mreq(ns+ns_sf,  nq),          ! coefficients of 3rd body
     &      sub_vreq(ns+ns_sf,  nq),          ! subtraction of "vreq" of
                                              ! each species in each equation
     &    sigma_vreq(           nq),          ! change of "sub_vreq" in each equation
     &       vreq_a7( 7,2,      nq),          ! summation of "vreq" * "a7"
                                              ! in each equation
     &          stat=ierr1 )
      if(ierr1/=0) then
        write(ifle,*) '### error : allocation failed'
        ierror = 1
      endif
!---------------------------
! --- initialization
!---------------------------
      ireq = 0
      vreq = 0
      preq = 0
      creq = 0
      mreq = -1
      sub_vreq = 0
      sigma_vreq = 0
      vreq_a7 = 0
!---------------------------
!-< 2. Input namelist >-
!---------------------------
!--< 2.1 read namelist >--
!---------------------------
!      model = ' '
!      rewind ifli
!      read(ifli,chemreac,iostat=ios)            ! read of namelist
!      call nml_listno(6,mdllst,model,imodl)
!      call nml_chkchr0(ifle,'model',model,imodl,ierr1)
!      if( ierr1/=0 ) then
!        write(ifle,*) "### error : unidentified model name"
!        ierror = 1
!        return
!      endif
!----------------------------------------------------------
! --- correspond to first flag for chemical reaction model
!----------------------------------------------------------
!      ireq(1:nq) = imodl
!
      rewind ifli
!
      do 200 q=1,nq
        ignition_iter=0
        model= ' '
        kind=' '
        nleft = ' '
        nright= ' '
        cleft = undef
        cright= undef
!/ overall /
        creac = undef    ! 
        a     = undef    ! Pre-exponential Factor (cm**3/mol-s)
        alpha = undef    !
        beta  = undef    !
        e     = undef    ! Activation Energy (kcal/mol)
!/ fire /
        Y_SOOT =undef    ! The fraction of fuel mass converted into SOOT
        C_fuel =undef    ! number of C in fuel Moleculor
        nu_fuel=undef    ! stoichiometric coefficient of fuel
        WM_fuel=undef    ! Molecular weight of fuel [kg/mol]
!/ elementary /
        af    = undef           ! k = (af)*T^(alpahf)*EXP(-(ef)/(R*T))
        alphaf= undef
        ef    = undef
        ab    = undef           ! backward reaction
        alphab= undef
        eb    = undef
        M     = undef           ! coefficients of 3rd body
!/ eddy break up /
        cr1   = undef
        cr2   = undef
        r     = undef
!-------------------------------------
        read(ifli,chemreac,iostat=ios)                  ! read of namelist
        call nml_errmsg0(ifle,ios,'chemreac',ierr1)     ! if ios>0, ierr1=1
        if( ierr1/=0 ) then
          write(ifle,*) "### error : reading namelist 'chemreac'"
          ierror = 1
          goto 9999
        endif
!-----------------------------
! --- identify reaction kind
!-----------------------------
        call nml_listno(2,chemkd,kind,kd)
        call nml_chkchr0(ifle,'kind',kind,kd,ierr1)
!
        if(ierr1/=0) stop '&chemreac/kind'
	if(kd==isurface) then
	  lsuf=.true.
	endif
!-------------------------
! --- identify model no.
!-------------------------
        call nml_listno(15,mdllst,model,imodl)
        call nml_chkchr0(ifle,'model',model,imodl,ierr1)
        if( ierr1.ne.0 ) goto 9999
        ireq(q)=imodl
!--------------------------------------------------
! --- check duplication between left and right side
!--------------------------------------------------
!        if(imodl/=elemnt.and.imodl/=userreac) then
!          call chk_lr(nleft,nright)
!        endif
! --- check duplication of left side
!        call chk_dup('nleft',nleft)
! --- check duplication of right side
!        call chk_dup('nright',nright)
!--------------------------------------------
! --- pick up number of species in each side
!--------------------------------------------
!        call pick_up_species
!        select case(imodl)
!          case(1)
!            call read_overall
!          case(2)
!            call read_elementary
!          case(3)
!            call read_eddy_break_up
!          case(4)
!            call read_overall_fire
!          case(5)
!            ! user reaction model
!          case(6)
!            ! surface stick reaction model
!        end select
 200  continue
!------------------
! --- chenk
!------------------
        if(ns_sf>0.and.(.NOT.lsuf)) then
          write(ifle,'(2a)') 'ERR: NO surface reaction defined, ',
     &    'reset [&sizes/ncomp_surface=0]'
          stop
        endif

!--------------------------------
! --- display of reading results
!--------------------------------
!
!      select case( imodl )
!        case(1)
!          call disp_overall
!        case(2)
!          call disp_elementary
!      end select
!
!----------------------
! --- final error check
!----------------------
!
      return
 9999 continue
      ierror=1
      if( ierror/=0 ) then
        write(ifle,*) "total error : ",ierror,modnam,subnam
        stop ': at module_chemreac'
      endif
!
      return
!
!///////////////////////////////////////////////////////////////////////
!
      contains
!--------------------------------------
! --- classification number of species
!--------------------------------------
!=================================================
       integer function c_no(name,str1,value,str2)
!=================================================
        character(lenspc),intent(in) :: name
        character,intent(in) :: str1,str2
        real*8,intent(in)  :: value
!
        c_no=chk_nam(ns,spcnam,name)
        if(c_no==0.and.ns_sf/=0) then
          c_no=chk_nam(ns_sf,spcnamsf,name)
          if(c_no/=ns+ns_sf+1) c_no=c_no+ns
        endif
! --- if "c_no" is identified(c_no>0), ierr1=0
        call nml_chkchr0(ifle,str1,name,c_no,ierr1)
        if( ierr1/=0 ) then
          write(ifle,*) "did not identify species name in equation",q
          ierror = 1
        endif
        if(value==iundef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'lack of data : ',str2,'(i)'
          write(ifle,*) 'i =',c_no
          ierror = 1
        elseif( value<1 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) str2,'(i) = ',value
          write(ifle,*) 'it must be > 0'
          write(ifle,*) 'i =',c_no
          ierror = 1
        endif
      end function c_no
!
!===========================================================
      subroutine pick_up_species  ! pick up species number
!===========================================================
      implicit none
      integer :: s
      nol = 0            ! number of species in left side
      nor = 0            ! number of species in right side
      do s=1,ns+ns_sf
        if(nleft(s)/=' ') then         ! left side
          j=c_no(nleft(s),'nleft',cleft(s),'cleft')
          if(j<=ns+ns_sf) then
!--------------------------------------------
! --- stoichiometric coefficient of left side
!--------------------------------------------
            vreq(j,1,q)=cleft(s)
            nol=nol+1
          elseif(j==ns+ns_sf+1) then
            mreq(1:ns+ns_sf,q)=M(1:ns+ns_sf)       ! 3rd body coefficient
          endif
          ispc(nol)=j
        endif
        if(nright(s)/=' ') then        ! right side
          j=c_no(nright(s),'nright',cright(s),'cright')
          if(j<=ns+ns_sf) then
!----------------------------------------------
! --- stoichiometric coefficient of right side
!----------------------------------------------
            vreq(j,2,q)=cright(s)
            nor=nor+1
          endif
        endif
      enddo
!
      do s=1,ns+ns_sf
        sub_vreq(s,q)=vreq(s,2,q)-vreq(s,1,q)
        sigma_vreq(q)=sigma_vreq(q)+sub_vreq(s,q)
      enddo
!
      if(nol<1 .or. nor<1) then               ! error check
        write(ifle,*) "### error-1 : data error"
        write(ifle,*) "lack of data in equation ",q
        ierror = 1
        write(ifle,*) "### error : nol & nor= ",nol,nor
        stop 'nol<1 or nor<1'
      elseif(imodl==2 .and. (nol>3 .or. nor>3)) then
        write(ifle,*) "### error-2 : data error"
        write(ifle,*) "too much data in equation ",q
        write(ifle,*) "elementary model affords 3 species in side."
        ierror = 1
      elseif(imodl==3 .and. nol/=2) then
        write(ifle,*) "### error-3 : data error"
        write(ifle,*) "no. of elements of nleft =",nol
         write(ifle,*) "it must be 2"
        write(ifle,*) "nleft(1) = name of fuel"
        write(ifle,*) "nleft(2) = name of oxidizing agent"
        ierror = 1
      endif
      end subroutine pick_up_species
!
!=================================================
      subroutine read_overall           ! imodl==1
!=================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      
      do i=1,nol
        if( nml_comp_eq("creac(i)",creac(i),undef,i,ifle,ierror) ) then
          if( .not.nml_comp_ge("creac(i)",creac(i),0,i,ifle,ierror) )
     &      exit
        endif
      enddo
      k=1
      do s=1,ns+ns_sf
        j=c_no(nleft(s),"nleft",cleft(s),"cleft")
! --- if j is identified(j>0), ierr1=0
        call nml_chkchr0(ifle,'nleft',nleft(s),j,ierr1)
        if( ierr1/=0 ) then
          write(ifle,*) "'nleft' is not identified in equation",q
          ierror = 1
        endif
! --- index number of mol concentration(overall)
        if( j<=ns+ns_sf ) creq(j) = creac(k)
        if( k<nol ) then
          k = k+1
        else
          exit
        endif
      enddo
      if(.not.nml_comp_eq("a",a,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("a",a,0.0d0,1,ifle,ierror) )
     &    preq(1,1,q)=a
      endif
      if(.not.nml_comp_eq("alpha",alpha,undef,1,ifle,ierror) )
     &  preq(2,1,q)=alpha
      if(.not.nml_comp_eq("beta",beta,undef,1,ifle,ierror) )
     &  preq(3,1,q)=beta
      if(.not.nml_comp_eq("e",e,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("e",e,0.0d0,1,ifle,ierror) )
     &    preq(4,1,q)=e
      endif
      end subroutine read_overall
!
!=======================================================
      subroutine read_overall_fire           ! imodl==1
!=======================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      
      character(lenspc) :: soot_nam,co_nam,o2_nam,co2_nam
      soot_nam='SOOT'
      co_nam='CO'
      o2_nam='O2'
      co2_nam='CO2'
!     
      do i=1,nol
        if(nml_comp_eq("creac(i)",creac(i),undef,i,ifle,ierror)) then
          if(.not.nml_comp_ge("creac(i)",creac(i),0.d0,i,ifle,ierror))
     &      exit
        endif
      enddo
      k=1
      do s=1,ns+ns_sf
        j=c_no(nleft(s),"nleft",cleft(s),"cleft")
        call nml_chkchr0(ifle,'nleft',nleft(s),j,ierr1)
        if(ierr1/=0) then
          write(ifle,*) "'nleft' is not identified in equation",q
          ierror = 1
        endif
        if(j<=ns+ns_sf) creq(j)=creac(k)
        if(k<nol) then
          k=k+1
        else
          exit
        endif
      enddo
!
      if(.not.nml_comp_eq("a",a,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("a",a,0.0d0,1,ifle,ierror) )
     &    preq(1,1,q)=a
      endif
!
      if(.not.nml_comp_eq("alpha",alpha,undef,1,ifle,ierror))
     &    preq(2,1,q)=alpha
!
      if(.not.nml_comp_eq("beta",beta,undef,1,ifle,ierror))
     &    preq(3,1,q)=beta
!
      if(.not.nml_comp_eq("e",e,undef,1,ifle,ierror)) then
        if(nml_comp_gt("e",e,0.0d0,1,ifle,ierror))
     &    preq(4,1,q)=e
      endif
!
      if(.not.nml_comp_eq("WM_fuel",WM_fuel,undef,1,ifle,ierror).and.
     &   .not.nml_comp_eq("C_fuel",C_fuel,undef,1,ifle,ierror).and.
     &   .not.nml_comp_eq("nu_fuel",nu_fuel,undef,1,ifle,ierror).and.
     &   .not.nml_comp_eq("Y_SOOT",Y_SOOT,undef,1,ifle,ierror)) then
!
        Y_CO=0.37d0*Y_SOOT
     &             +12.d0*C_fuel*0.0014d0/WM_fuel/nu_fuel
        nu_SOOT=WM_fuel/WM_C*Y_SOOT*nu_fuel
        nu_CO=  WM_fuel/WM_CO*Y_CO*nu_fuel
! --- SOOT Nu
        j=chk_nam(ns,spcnam,soot_nam)
        if(j==0.and.ns_sf/=0) then
          j=chk_nam(ns_sf,spcnamsf,soot_nam)
          if(j/=ns+ns_sf+1) j=j+ns
        endif
        vreq(j,2,q)=nu_SOOT
! --- CO Nu
        j= chk_nam(ns,spcnam,co_nam)
        if(j==0.and.ns_sf/=0) then
          j=chk_nam(ns_sf,spcnamsf,co_nam)
          if(j/=ns+ns_sf+1) j=j+ns
        endif
        vreq(j,2,q)=nu_CO
! --- O2 Nu
        j=chk_nam(ns,spcnam,o2_nam)
        if(j==0.and.ns_sf/=0) then
          j=chk_nam(ns_sf,spcnamsf,o2_nam)
          if(j/=ns+ns_sf+1) j=j+ns
        endif
        vreq(j,1,q)=vreq(j,1,q)-nu_CO-nu_SOOT
! --- CO2 Nu
        j=chk_nam(ns,spcnam,co2_nam)
        if(j==0.and.ns_sf/=0) then
          j=chk_nam(ns_sf,spcnamsf,co2_nam)
          if(j/=ns+ns_sf+1) j=j+ns
        endif
        vreq(j,2,q)=vreq(j,2,q)-nu_CO-nu_SOOT
        endif
!
       write(ifll,'(a)')' ==== stoichiometric coefficients ===='
        do s=1,ns+ns_sf
        write(ifll,'(a7,i3,a9,i3,10x,a8)')
     &           'Reac.= ',q,'   spec.no.=',s,spcnam(s)
        write(ifll,'(5x,e13.5,5x,e13.5)') vreq(s,1,q),vreq(s,2,q)
        enddo
      end subroutine read_overall_fire
!
!===================================================
      subroutine read_elementary        ! imodl==2
!===================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
!      
      if(.not.nml_comp_eq("af",af,undef,q,ifle,ierror)) then
        if(nml_comp_gt("af",af,0.0d0,q,ifle,ierror))
     &    preq(1,1,q) = af
      endif
      if(.not.nml_comp_eq("alphaf",alphaf,undef,q,ifle,ierror) )
     &    preq(2,1,q)=alphaf
      if(.not.nml_comp_eq("ef",ef,undef,q,ifle,ierror))
     &    preq(3,1,q)=ef
! --- =undef:calc. backward reaction with forward one (Both F and B)
! --- =0    :NOT consider backward reaction           (Only F)
! --- other :calc. backward reaction without forward one (Only B)
      preq(1,2,q) = ab
      if( preq(1,2,q)/=undef .and. preq(1,2,q)/=0.d0 ) then
          preq(2,2,q)=alphab
          preq(3,2,q)=eb
      endif
      do k=1,2
        do j=1,7
! --- summation of "vreq" * "a7" in each equation
          do s=1,ns
          vreq_a7(j,k,q)=vreq_a7(j,k,q)+sub_vreq(s,q)*a7(j,k,s)
          enddo
        enddo
      enddo
      end subroutine read_elementary
!
!=================================================
      subroutine read_eddy_break_up     ! imodl==3
!=================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      
      preq(1,1,q) = dble(ispc(1))
      preq(2,1,q) = dble(ispc(2))
      if(.not.nml_comp_eq("cr1",cr1,undef,q,ifle,ierror)) then
        if(nml_comp_gt("cr1",cr1,0.0d0,q,ifle,ierror))
     &    preq(3,1,q) = cr1
      endif
      dd2 = 0.d0
      if(cr2==undef) then
        cr2 = 0.d0
        dd2 = 10.d0
      elseif( cr2<=0.d0 ) then
        if( nml_comp_gt("cr2",cr2,0.0d0,q,ifle,ierror) ) return
      endif
      preq(4,1,q) = cr2
      if( .not.nml_comp_eq("r",r,undef,q,ifle,ierror) ) then
        if( nml_comp_gt("r",r,0.0d0,q,ifle,ierror) )
     &    preq(5,1,q) = r
      endif
!
      preq(6,1,q) = dd2
      end subroutine read_eddy_break_up
!
!< #2.1  check data >
!
!=================================================
      subroutine chk_dup(str,value)
!=================================================
      character(*),intent(in) :: str,value(mcomp+1)
      do s=1,(ns+ns_sf)-1
        if( value(s)==' ' ) cycle
        do j=s+1,(ns+ns_sf)
          if( value(j)==' ' ) cycle
          if( value(j)==value(s) ) then
            write(ifle,*) '### error : data error'
            write(ifle,*) 'duplicated data : ',str,
     &                    '(i)=',str,'(j)=',value(s)
            write(ifle,*) 'i,j =',s,j
            ierror = 1
            stop 'stop at chk_dup'
          endif
        enddo
      enddo
      end subroutine chk_dup
!
!---------------------------------------------------
! --- check duplication between left and right side
!===================================================
      subroutine chk_lr(left,right)
!===================================================
      character(lenspc),intent(in) :: left(mcomp+1)
                                   ! name of species on left side
      character(lenspc),intent(in) :: right(mcomp+1)
                                   ! name of species on right side
      do s=1,(ns+ns_sf)+1
        if(left(s)==' ') cycle
        do j=1,ns+ns_sf+1
          if( right(j)==' ' ) cycle
          if( left(s)==right(j) ) then
            write(ifle,*) '### error : data error'
            write(ifle,*) 'duplicated data : nleft(i)=nright(j)='
     &                    ,left(s)
            write(ifle,*) 'i,j =',s,j
            ierror = 1
            stop 'stop at chk_lr'
          endif
        enddo
      enddo
      end subroutine chk_lr
!
!-------------------------------------------------
! --- check species name
!=================================================
      integer function chk_nam(ncom,lst,wrd)
!=================================================
!  chk_nam : = 0; unknown name (did not declare in namelist 'species')
!          : > 0; identified name (nn: number of identified name)
!          : =ncom+1; 3rd body
        integer,intent(in)  :: ncom
        character(lenspc),intent(in)  :: lst(ncom)
        character(lenspc),intent(in)  :: wrd
        chk_nam = 0
        do s=1,ncom
          if( wrd==lst(s) ) then
            chk_nam = s
            exit
          elseif( wrd=='M' ) then
            chk_nam = ns+ns_sf+1
            exit
          endif
        enddo
      end function chk_nam
!
!-------------------------------------------------
! --- following is display part
!=================================================
      subroutine disp_overall
!=================================================
      call disp_head
      call disp_eq1(1)          ! left side
      write(ifll,"(a5,$)") arrow(preq(1,2,1))
      call disp_eq1(2)          ! right side
      call disp_coefficients(1,4)
      call disp_creq
      write(ifll,*)
      end subroutine disp_overall
!
!=================================================
      subroutine disp_elementary
!=================================================
      call disp_head
      do q=1,nq
        call disp_equation(q)
        write(ifll,*)
      enddo
!      call disp_vreq_a7
      end subroutine disp_elementary
!
!=================================================
      subroutine disp_vreq_a7
!=================================================
      do q=1,nq
        write(ifll,"(i3,$)") q
        do k=1,2
          write(ifll,"(e12.3,$)") (vreq_a7(j,k,q),j=1,7)
          if( k==1 ) write(ifll,"(/,a,$)") "   "
        enddo
        write(ifll,*)
      enddo
      end subroutine disp_vreq_a7
!
!=================================================
      subroutine disp_head
!=================================================
      write(ifll,"(80('-'))")
      write(ifll,"(2a)") "chemical reaction model     : ",
     &                   mdllst(ireq(1))
      write(ifll,"(a,i3)")"number of chemical equation : ",nq
      end subroutine disp_head
!
!=================================================
      subroutine disp_creq
!=================================================
      write(ifll,"(/,7x,a,$)") "index number     : "
      do s=1,ns
        if( creq(s)/=0 ) then
          write(ifll,"(2x,a6,$)") spcnam(s)
          write(ifll,"(a,f5.1,$)") "->",creq(s)
        endif
      enddo
      end subroutine disp_creq
!-------------------------------------------------
! --- display 1 equation for overall reaction
!=================================================
      subroutine disp_eq1(isw)
!=================================================
      integer,intent(in) :: isw  ! switch for side of equation
      k = 0
      do s=1,ns
        if( vreq(s,isw,1)/=0 ) then
          if( k>0 ) write(ifll,"(a3,$)") " + "
          write(ifll,"(f5.1,$)") vreq(s,isw,1)
          write(ifll,"(1x,a6,$)") spcnam(s)
          k = k+1
        endif
      enddo
      end subroutine disp_eq1
!
!-------------------------------------------------
! --- display 1 equation for elementary reaction
!=================================================
      subroutine disp_equation(q)
!=================================================
      integer,intent(in) :: q
                              ! classification number of equation
      write(ifll,"(i3,2x,a,$)") q,str_eqs(1,q)  ! left side
      write(ifll,"(a5,$)") arrow(preq(1,2,q))
      write(ifll,"(a,$)") str_eqs(2,q)          ! right side
      call disp_coefficients(q,3)
!      write(ifll,"(/,11x,9f4.0,$)") sub_vreq(1:ns,q)
      end subroutine disp_equation
!
!-------------------------------------------------
! --- display some coefficients
!=================================================
      subroutine disp_coefficients(q,n)
!=================================================
      integer,intent(in) :: q,n
      write(ifll,"(/,7x,a,$)") "forward reaction : "
      write(ifll,"(e12.3,$)") (preq(i,1,q),i=1,n)
      write(ifll,"(/,7x,a,$)") "backward reaction: "
      if( preq(1,2,q)==undef ) then
        write(ifll,"(a,$)") "<-- calculation"
      elseif( preq(1,2,q)==0 ) then
        write(ifll,"(a,$)") "<-- NOT consider"
      else
        write(ifll,"(e12.3,$)") (preq(i,2,q),i=1,n)
      endif
      if( mreq(1,q)/=-1 ) then
        write(ifll,"(/,7x,a,$)") "coefficients of M: "
        write(ifll,"(f5.1,$)") (mreq(s,q),s=1,ns)
      endif
      end subroutine disp_coefficients
!
!-------------------------------------------------
! --- display 1 equation's side
!=================================================
      character(24) function str_eqs(isw,q)
!=================================================
        integer,intent(in) :: isw ! switch
        integer,intent(in) :: q   ! classification number of equation
        str_eqs = ""
        k = 0
        do s=1,ns
          if( vreq(s,isw,q)/=0 ) then
            if( k==0 ) then
              str_eqs(1:6) = spcnam(s)
              k = k+1
              if( vreq(s,isw,q)>1 ) then
                str_eqs(7:9) = " + "
                str_eqs(10:15) = spcnam(s)
                k = k+1
              endif
            elseif( k==1 ) then
              str_eqs(7:9) = " + "
              str_eqs(10:15) = spcnam(s)
              k = k+1
              if( vreq(s,isw,q)>1 ) then
                str_eqs(16:18) = " + "
                str_eqs(19:24) = spcnam(s)
                k = k+1
              endif
            elseif( k==2 ) then
              str_eqs(16:18) = " + "
              str_eqs(19:24) = spcnam(s)
            endif
          endif
        enddo
        if( mreq(isw,q)/=-1 .and. k<3 ) then
          str_eqs(6*k+3*(k-1)+1:6*k+3*(k-1)+4) = " + M"
        endif
      end function str_eqs
!
!=================================================
      character(5) function arrow(arg)  ! display arrow
!=================================================
        real*8,intent(in) :: arg
        if( arg/=0 ) then    ! reversible reaction
          arrow = " <=> "
        else                 ! irreversible reaction
          arrow = "  => "
        endif
      end function arrow
      end subroutine inputdata
!
!-------------------------------------------------
! difference of
! (Gibbs energy on equation)/((universal gas constant)*(temperature)) [-]
!=================================================
      real*8 function dg_RT(t,q)
!=================================================
        use module_species,only : judg_temp
        real*8, intent(in) :: t   ! temperature [K]
        integer,intent(in) :: q   ! classification number of equation
        integer jt                ! judgment of temperature
        jt = judg_temp(1,t)
        dg_RT = vreq_a7(1,jt,q)*(1-log(t))
     &    +t*( -vreq_a7(2,jt,q)*0.5
     &    +t*( -vreq_a7(3,jt,q)/6
     &    +t*( -vreq_a7(4,jt,q)/12
     &    +t*( -vreq_a7(5,jt,q)*0.05 ))))
     &         +vreq_a7(6,jt,q)/t
     &         -vreq_a7(7,jt,q)
      end function dg_RT
!
!--------------------------------------------------------------------------
! --- difference of (Gibbs energy on equation)/(universal gas constant) [K]
!=================================================
      real*8 function dg_t(ns,t,q)
!=================================================
        use module_species,only : h_t,s_0
        integer,intent(in) :: ns  ! number of species
        real*8, intent(in) :: t   ! temperature [K]
        integer,intent(in) :: q   ! classification number of equation
        integer s                 ! classification number of species
        dg_t = 0
        do s=1,ns
          if( sub_vreq(s,q)==0 ) cycle
          dg_t = dg_t+sub_vreq(s,q)*(h_t(s,t)-t*s_0(s,t))
        enddo
      end function dg_t
!
!-------------------------------------------------
! --- revised mol concentration by 'mreq' [mol/m3]
!=================================================
      real*8 function revised_mc(ns,q,mc)
!=================================================
        integer,intent(in) :: ns      ! number of species
        integer,intent(in) :: q       ! classification number of equation
        real*8, intent(in) :: mc(ns)  ! mol concentration [mol/m3]
        integer s                     ! classification number of species
        revised_mc = 0
        do s=1,ns                     ! revised mol concentration
                                      ! by 'mreq' [mol/m3]
          revised_mc = revised_mc+mreq(s,q)*mc(s)
        enddo
      end function revised_mc
!
      end module module_chemreac
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine modutl_chkset(id,jset,iset,sbnm)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(*),intent(in)    :: id,sbnm
      integer     ,intent(in)    :: jset
      integer     ,intent(inout) :: iset
!
      if( id.ne.'=' ) then
        if( iset.ge.jset ) return
      else
        iset=iset+1
        if( iset.eq.jset ) return
      endif
      write(*,*) '### program error -1- ',sbnm,'(modutl_chkset)'
      write(*,*) 'id,jset,iset=',id,jset,iset
      stop 999
      end subroutine modutl_chkset


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine modutl_setnow(time,nset,nval,nwdt,iwd,wdt,wd)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. Interpolate values of variables with time
!    at present time from time tables
!
! [dummy arguments]
      real*8 ,intent(in)    :: time
      integer,intent(in)    :: nset,nval,nwdt
      integer,intent(in)    :: iwd(0:nset)
      real*8 ,intent(in)    :: wdt(0:nval,nwdt)
      real*8 ,intent(inout) :: wd (0:nval,nset)
!!!      include'interfacexxx'
!
! [local entities]
      integer :: i,k,n,ks,ke,kx,k1,k2,ncyc
      real*8  :: tim0,tim1,tim2,tcyc,fc1,fc2
!
!
!
! 1. tcyc=wd(0,n) is cycle of the time table (in),
!    wd(i,n) (i>0) are values at present time (out).
!    If tcyc is not positive, the time table is not cyclic.
! 2. ke=iwd(n) is the end of each time table,
!    ks=iwd(n-1)+1 is the start of that.
! 3. wdt(0,k) (ks<=k<=ke) are times of the time table,
!    wdt(i,k) (i>0) are values at time wdt(0,k) respectively.
! 4. If tcyc > 0 ( case of cyclic time table ),
!    the values at wdt(0,ks)+tcyc are wdt(i,ks) (i>0).
!    The condition wdt(0,ke) < wdt(0,ks)+tcyc is assumed.
! 5. If tcyc <= 0 ( case of not cyclic time table ),
!    the values at time < wdt(0,ks) are wdt(i,ks) and
!    the values at time > wdt(0,ke) are wdt(i,ke) (i>0).
! 6. If the no. of points in the time table (ke-ks+1) is less than 2,
!    the values do not vary with time and then
!    no procedure is performed.
!
      do 100 n=1,nset
      if( iwd(n).lt.iwd(n-1)+2 ) goto 100
      ks=iwd(n-1)+1
      ke=iwd(n)
      tim0=time
      kx=ke-1
      if( wd(0,n).gt.0.d0 ) then
        kx=ke
        tcyc=wd(0,n)
        tim1=wdt(0,ks)
        tim0=tim0-tim1
        ncyc=int(tim0/tcyc)
        tim0=tim0-dble(ncyc)*tcyc
        if( tim0.lt.0.d0 ) tim0=tim0+tcyc
        tim0=tim0+tim1
      endif
      k1=ks
      do 101 k=ks+1,kx
      if( wdt(0,k).lt.tim0 ) k1=k
  101 continue
      tim1=wdt(0,k1)
      if( k1.lt.ke ) then
        k2=k1+1
        tim2=wdt(0,k2)
      else
        k2=ks
        tim2=wdt(0,ks)+tcyc
      endif
      fc2=max(0.d0,min(1.d0,(tim0-tim1)/(tim2-tim1)))
      fc1=1.d0-fc2
      do 102 i=1,nval
      wd(i,n)=fc1*wdt(i,k1)+fc2*wdt(i,k2)
  102 continue
  100 continue
!
      return
!
      end subroutine modutl_setnow
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_gf
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!< module for GF >
      DATA IUT0    / 0 /
      DATA IUT5    / 5 /
      DATA IUT6    / 6 /
      DATA IWRITE  / 2 /
      DATA INAME   / 1 /
      DATA IDIM    / 3 /
      DATA IRESV   / 3 /
      DATA ICHECK  / 999999 /
C
      INTEGER IACT,NCOMFL,NCOMST
C
      PARAMETER ( MCOM = 10 )
      CHARACTER*60 COMFLE(MCOM)
      CHARACTER*60 COMSET(MCOM)
!
      end module module_gf
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_usersub
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_usersub)'
!
!
!
      integer,parameter,public :: usrno=0,usryes=1,usryesyes=2
      integer,public :: outusr,inlusr,iniusr,oulusr,walusr
      integer,public :: src_uvw,src_r,src_t,src_rans,src_fire
!
! usryes  : user subroutine is used by user
! usrno   : user subroutine is NOT used by user
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      implicit none
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy qrgument]
!
      
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! --- [local entities]
!
      character(6),parameter :: cond(3)=(/'yes   ','no    ','confrm'/)
      integer :: iset=0,co(3)=(/1,0,2/)
      integer :: ios,ierr1
!
! --- [namelist]
!
      integer       :: kout=1,kinl=1,ksrcu=1,kini=1,ksrcr=1,ksrct=1,
     &                 ksrcys=1,ksrcrans=1,kglb=1,kfire=1,koul=1,
     &                 kwal=1
      character(20) :: output,initial,initial_E2P,source_BC_E2P
      character(20) :: source_BC_CAVI
      character(20) :: inlet,outlet,wall,static
      character(20) :: source_uvw,
     &                 source_r,source_t,source_rans,
     &                 source_fire,MHD_initial,MHD_BC,MHD_CURRENT
      character(20) :: radpro,radpara,initial_particle,sata_rough_wall,
     &                 user_defined_region
      namelist /usrsub/output,inlet,initial,initial_E2P,source_uvw,
     &                 source_r,source_t,source_rans,
     &                 outlet,wall,static,source_BC_E2P,
     &                 source_fire,source_BC_CAVI,
     &                 MHD_initial,MHD_BC,MHD_CURRENT,
     &		       radpro,radpara,initial_particle,sata_rough_wall,
     &                 user_defined_region
!
      ierr1=0
! --- initialize flags
      outusr=0
      inlusr=0
      oulusr=0
      walusr=0
      src_uvw=0
      iniusr=0
      src_r=0
      src_t=0
      src_rans=0
      src_fire=0
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
      output=' '
      inlet=' '
      outlet=' '
      wall=' '
      source_uvw=' '
      initial=' '
      source_r=' '
      source_t=' '
      source_rans=' '
      source_fire=' '
!
      rewind ifli
      read(ifli,usrsub,iostat=ios)
      if( ios.lt.0 ) return
!
      call nml_errmsg0(ifle,ios,'usrsub',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      call nml_listno(3,cond,output,kout)
      call nml_chkchr0(ifle,'output',output,kout,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [output] into [&usrsub] in ',cntlnam
         goto 9999
      endif
!
      call nml_listno(3,cond,inlet,kinl)
      call nml_chkchr0(ifle,'inlet',inlet,kinl,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [inlet] into [&usrsub] in ',cntlnam
         goto 9999
      endif
!
      call nml_listno(3,cond,outlet,koul)
      call nml_chkchr0(ifle,'outlet',outlet,koul,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [outlet] into [&usrsub] in ',cntlnam
         stop ': at module_usersub/inputdata'
      endif
!
      call nml_listno(3,cond,wall,kwal)
      call nml_chkchr0(ifle,'wall',wall,kwal,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [wall] into [&usrsub] in ',cntlnam
         stop ': at module_usersub/inputdata'
      endif
!
      call nml_listno(3,cond,source_uvw,ksrcu)
      call nml_chkchr0(ifle,'source_uvw',source_uvw,ksrcu,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [source_uvw] into [&usrsub] in ',cntlnam
         goto 9999
      endif
!
      call nml_listno(3,cond,initial,kini)
      call nml_chkchr0(ifle,'initial',initial,kini,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [initial] into [&usrsub] in ',cntlnam
         goto 9999 
      endif
!
      call nml_listno(3,cond,source_r,ksrcr)
      call nml_chkchr0(ifle,'source_r',source_r,ksrcr,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [source_r] into [&usrsub] in ',cntlnam
         goto 9999
      endif
!
      call nml_listno(3,cond,source_t,ksrct)
      call nml_chkchr0(ifle,'source_t',source_t,ksrct,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [source_t] into [&usrsub] in ',cntlnam
         goto 9999
      endif
!
      call nml_listno(3,cond,source_rans,ksrcrans)
      call nml_chkchr0(ifle,'source_rans',source_rans,ksrcrans,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [source_rans] into [&usrsub] in ',cntlnam
         goto 9999
      endif
!
      call nml_listno(3,cond,source_fire,kfire)
      call nml_chkchr0(ifle,'source_fire',source_fire,kfire,ierr1)
      if( ierr1.ne.0 ) then
         write(ifle,*) 'add [source_fire] into [&usrsub] in ',cntlnam
         goto 9999
      endif
!----------------------
! --- define flags
!----------------------
      outusr=co(kout)
      inlusr=co(kinl)
      oulusr=co(koul)
      walusr=co(kwal)
      iniusr=co(kini)
      src_uvw=co(ksrcu)
      src_r=co(ksrcr)
      src_t=co(ksrct)
      src_rans=co(ksrcrans)
      src_fire=co(kfire)
!
!      if(outusr.eq.usryes) then
!        write(ifll,*) '    ###    User Subroutine [user/usrout.f]'
!      elseif(inlusr.eq.usryes) then
!        write(ifll,*) '    ###    User Subroutine [user/usrinl.f]'
!      elseif(srcusr.eq.usryes) then
!        write(ifll,*) '    ###    User Subroutine [user/usrsrc.f]'
!      elseif(iniusr.eq.usryes) then
!        write(ifll,*) '    ###    User Subroutine [user/usrini.f]'
!      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_usersub
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_hpc
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_hpc)'
!
      integer,save                  :: NPETOT,nbsldG
      integer,save,allocatable      :: nsldT(:,:),nsldbc(:),
     &                                 MAT_BCSLD(:,:)
      CHARACTER*80,save,allocatable :: nmsldL(:,:),nmsldG(:)
      
!
!/////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,nbcndx,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy qrgument]
!
      integer,intent(in)   :: ifli,ifll,ifle,nbcndx
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! --- [local entities]
!
      integer :: ierr1
      integer :: iset=0
!
! --- [namelist]
!
      integer :: ncpu,MONITOR
      namelist /hpc/ ncpu,MONITOR
!      
!
      ierr1=0
      ierror=0
      NPETOT=1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
      ncpu=1
!
      rewind ifli
      read(ifli,hpc,iostat=ios)
!      write(ifll,hpc)
      call nml_errmsg0(ifle,ios,'hpc',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      NPETOT=ncpu
!
      allocate(nsldT(0:nbcndx,NPETOT),nsldbc(NPETOT),stat=ierr1)
      allocate(nmsldL(nbcndx,NPETOT),nmsldG(nbcndx),stat=ierr1)
!      
      nsldT(:,:)=0
      nsldbc(:)=0
      nmsldL(:,:)=' '
      nmsldG(:)=' '
      if(ierr1.ne.0) then
         write(ifle,*) '### error : allocation failed: nsldT'
         goto 9999
      endif
!
      return
 9999 continue
      ierror=1
      write(ifle,*) modnam,subnam
      end subroutine inputdata
      end module module_hpc
!
!+++++Modified by T.Unemura 040430++++++++++++++++++++++++++++++++++++++
! !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       module module_flow_sound
! !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       character(20),parameter,private :: modnam='(module_flow_sound)'
! !
! ! --- [module arguments]
! !
!       logical,public           :: lsound,lshear
!       integer,save,allocatable :: isound(:),iname(:)
!       character*80,allocatable :: names(:)
!       integer :: nsond
! !
! !///////////////////////////////////////////////////////////////////////
!       contains
! !=================================================
!       subroutine inputdata
!      & (ifli,ifll,ifle,nbcnd,boundName,kdbcnd,kxnone,kdbuff,ierror)
! !=================================================
!       character(11),parameter :: subnam='(inputdata)'
! !
! !
! ! --- [dummy qrgument]
! !
!       integer,intent(in)    :: ifli,ifll,ifle,nbcnd,kxnone,kdbuff
!       integer,intent(in)    :: kdbcnd(0:4,nbcnd)
!       character(*),intent(in) :: boundName(nbcnd,2)
!       integer,intent(out)   :: ierror
! !
! ! --- [local entities]
! !
!       integer :: ierr1,ios,nb,nbs,kd
!       integer :: iset=0
!       integer :: ksound=0
!       character(3),parameter :: cond(2)=(/'yes','no '/)
!       character(7),parameter :: co(2)=(/'.true. ','.false.'/)
! !
! ! --- [namelist]
! !
!       character(20) :: flow_sound
!       character(80) :: name(50)
!       namelist /sound/ flow_sound,name
! !      
! !
!       ierr1=0
!       ierror=0
! !
! !
!       call modutl_chkset('=',1,iset,modnam//subnam)
! !
!       lsound=.false.
!       flow_sound=' '
!       name(:)=' '
! !
!       rewind ifli
!       read(ifli,sound,iostat=ios)
! 
!       call nml_errmsg0(ifle,ios,'sound',ierr1)
!       if( ierr1.ne.0 ) goto 9999
! !
!       call nml_listno(2,cond,flow_sound,ksound)
!       call nml_chkchr0(ifle,'flow_sound',flow_sound,ksound,ierr1)
!       if( ierr1.ne.0 ) goto 9999
! !
!       if(ksound.eq.1) then
!         lsound=.true.
!       elseif(ksound.eq.2) then
!         lsound=.false.
!       else
!         write(ifle,*) "ERR: flow_sound='yes' or 'no' in ',cntlnam
!         goto 9999
!       endif
!       if(.not.lsound) return
! !
!       do nb=1,nbcnd
!       do nbs=nb+1,nbcnd
!       if(name(nb).eq.name(nbs).and.name(nb).ne.' ') then
!         write(ifle,*) 'The name : ',name(nbs)(:len_trim(name(nbs))),
!      &  ' has been duplicated'
!         goto 9999
!       endif
!       enddo
!       enddo
! !
!       nsond=0
!       do nb=1,nbcnd
!       if(name(nb).ne.' ') then
!         nsond=nsond+1
!       endif
!       enddo
! !
!       allocate(names(nsond),iname(nsond))
!       nsond=0
!       names(:)=' '
!       do nb=1,nbcnd
!       if(name(nb).ne.' ') then
!         nsond=nsond+1
!         names(nsond)=name(nb)
!       endif
!       enddo
! !
!       allocate(isound(nbcnd))
!       isound(:)=0
!       iname(:)=0
!       do nbs=1,nsond
!       do nb=1,nbcnd
!       if(adjustl(names(nbs)).eq.adjustl(boundName(nb,1)))
!      &  then
!         isound(nb)=1
!         iname(nbs)=nb
!       endif
!       enddo
!       enddo
! !
!       do nbs=1,nsond
!        if(adjustl(names(nbs)).ne.' '.and.iname(nbs).eq.0) then
!          
!          write(ifle,*) '### error : namelist [sound]'
!          write(ifle,*)
!      &               names(nbs)(:len_trim(names(nbs))),
!      &              ' not found in boundary namelist'
!        goto 9999
!        endif
!       enddo
! !
!       do nbs=1,nsond
!         nb=iname(nbs)
!         kd=kdbcnd(0,nb)
!         if(kd.ne.kdbuff.and.kd.ne.kxnone) then
!           write(ifle,*) 
!      &    ' ERR: The boundary : ',names(nbs)(:len_trim(names(nbs))),
!      &    ' is NOT wall boundary'
!           goto 9999
!         endif
!       enddo
! !
!       return
!  9999 continue
!       ierror=1
!       write(ifle,*) modnam,subnam
!       end subroutine inputdata
!       end module module_flow_sound
! !
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_flow_sound
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_flow_sound)'
!
! --- [module arguments]
!
      logical :: soundObservFlag
      ! Flag of sound calculation enable/disable.
      integer :: numofSoundWall
      ! Number of sound souce walls.
      character*80,allocatable :: nameofSoundWall(:)
      ! List of sound wall name.
      integer,allocatable :: soundWall_bcNo(:)
      ! Boundary region index number of sound souce wall.
      real(8),allocatable :: soundWallCenterPosit(:,:)
      ! Position of center of sound source wall.
      integer :: numofSoundObserver
      ! Number of observer points.
      real(8),allocatable :: soundObservPosit(:,:)
      ! Positions of sound oberver points.
      integer,allocatable :: soundModelCode(:)
      ! Sound model code No. of the observer points.
      real(8),allocatable :: prsPrevSave(:,:),machPrevSave(:,:)
      ! Buffer for previous time step data.
      real(8),allocatable :: dflxPrevSave(:,:)
      ! Buffer for previous time step data.
      ! (Using for calculation of d/dtau)
      real(8) :: timePrevSave
      ! Time of the previous step.
!
!//////////////////////////////////////////////////////////////////////
      contains
!======================================================================
      subroutine inputdata
     & (ifli,ifll,ifle,cntlnam,
     6  nbcnd,boundName,kdbcnd,kxnone,kdbuff,ierror)
!======================================================================
      implicit none
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy qrgument] 
!
      integer,intent(in)    :: ifli,ifll,ifle,nbcnd,kxnone,kdbuff
      character(*),intent(in)  :: cntlnam
      integer,intent(in)    :: kdbcnd(0:4,nbcnd)
      character(*),intent(in) :: boundName(nbcnd,2)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer :: ierr1,ios
      integer :: iset=0
      integer :: ksound=0
      character(3),parameter :: cond(2)=(/'yes','no '/)
      character(7),parameter :: co(2)=(/'.true. ','.false.'/)
      character(6),parameter :: modelLabel(2)=(/'Curle ','FWH   '/)

      integer :: ic,jc,kc,nb,kd
      character*80,allocatable :: nameofSoundWallTmp(:)
!
! --- [namelist]
!
      character(20) :: enable_source,enable_observer
      character(80) :: wall_name
      real(8) :: position_x,position_y,position_z
      character*80 :: model
      namelist /sound_source/enable_source,wall_name
      namelist /sound_observer/enable_observer,position_x,position_y,
     &                         position_z,model
!
!
      ierr1=0
      ierror=0
      soundObservFlag=.true.
!
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
! --- Sound wall setting FROM HERE ------------------------------------------
!     Count up the number of enabled sound source walls.
      numofSoundWall=0
      rewind(ifli)
      do while(.true.)
        enable_source='yes'
        read(ifli,sound_source,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'sound_source',ierr1)
        if(ierr1/=0) then
          ierror=1
          soundObservFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_source,ksound)
        if(ksound==1) numofSoundWall=numofSoundWall+1
      end do
!     If there is no sound wall, set sound calculation flag to .false.
      if(numofSoundWall<=0) then
        soundObservFlag=.false.
        return
      end if
      
!     Prepare the list of wall name.
      allocate(nameofSoundWallTmp(1:numofSoundWall))
      
!     Rereading the namelist file.
      rewind(ifli)
      ic=0
      do while(.true.)
        enable_source='yes'
        wall_name=' '
        read(ifli,sound_source,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'sound_source',ierr1)
        if(ierr1/=0) then
          ierror=1
          soundObservFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_source,ksound)
        if(ksound==1) then
          ic=ic+1
          nameofSoundWallTmp(ic)=wall_name
        end if
      end do
!     Clear duplicated wall names and recount the number of walls.
      if(ic>=2) then
        do jc=1,ic-1
          if(nameofSoundWallTmp(jc)/=' ') then
            do kc=jc+1,ic
              if(nameofSoundWallTmp(kc)==nameofSoundWallTmp(jc))
     &            nameofSoundWallTmp(kc)=' '
            end do
          end if
        end do
        numofSoundWall=0
        do jc=1,ic
          if(nameofSoundWallTmp(jc)/=' ') 
     &            numofSoundWall=numofSoundWall+1
        end do
      end if
      allocate(nameofSoundWall(1:numofSoundWall))
      if(ic>=2) then
        kc=0
        do jc=1,ic
          if(nameofSoundWallTmp(jc)/=' ') then
            kc=kc+1
            nameofSoundWall(kc)=nameofSoundWallTmp(jc)
          end if
        end do
      else
        nameofSoundWall(:)=nameofSoundWallTmp(:)
      end if
      deallocate(nameofSoundWallTmp)

!     Search relation between sound source walls 
!      and boundary region index number.
      allocate(soundWall_bcNo(1:numofSoundWall))
      soundWall_bcNo(:)=0
      do jc=1,numofSoundWall
        do nb=1,nbcnd
          if(adjustl(nameofSoundWall(jc))==adjustl(boundName(nb,1))) 
     &      soundWall_bcNo(jc)=nb
        end do
      end do
!     Check the list.
      do jc=1,numofSoundWall
        if(soundWall_bcNo(jc)==0) then
          write(ifle,*)'### error : namelist [sound]'
          write(ifle,*)trim(nameofSoundWall(jc)),
     &                 ' not found in boundary namelist'
          ierror=1;return
        endif
      end do
!     Confirm the boundary conditions of sound source walls.
      do jc=1,numofSoundWall
        nb=soundWall_bcNo(jc)
        kd=kdbcnd(0,nb)
        if(kd.ne.kdbuff.and.kd.ne.kxnone) then
          write(ifle,*)
     &   ' ERR: The boundary : ',trim(nameofSoundWall(jc)),
     &   ' is NOT wall boundary'
          ierror=1;return
        endif
      end do
      
!      
      allocate(soundWallCenterPosit(1:numofSoundWall+1,1:3))
! --- Sound wall setting TO HERE ---------------------------------------
      
! --- Observer points setting FROM HERE --------------------------------
!     Count up the number of enabled obervers points
      numofSoundObserver=0
      rewind(ifli)
      do while(.true.)
        enable_observer='yes'
        read(ifli,sound_observer,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'sound_observer',ierr1)
        if(ierr1/=0) then
          ierror=1
          soundObservFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_observer,ksound)
        if(ksound==1) numofSoundObserver=numofSoundObserver+1
      end do
!     Even if there is no sound_oberver namelist,
!      Curle's equation is calculated.
      numofSoundObserver=numofSoundObserver+1
      
!     Preare arrays of observer points data.
      allocate(soundObservPosit(1:numofSoundObserver,1:3))
      allocate(soundModelCode(1:numofSoundObserver))
!     Reload name list file.
      rewind(ifli)
!       Data for default calculation of Curles equation
      ic=1
      soundObservPosit(1,:)=0.0d0
      soundModelCode(1)=1
      do while(.true.)
        enable_observer='yes'
        position_x=0.d0;position_y=0.d0;position_z=0.d0
        model='FWH'
        read(ifli,sound_observer,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'sound_observer',ierr1)
        if(ierr1/=0) then
          ierror=1
          soundObservFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_observer,ksound)
        if(ksound==1) then
          ic=ic+1
          soundObservPosit(ic,1)=position_x
          soundObservPosit(ic,2)=position_y
          soundObservPosit(ic,3)=position_z
          call nml_listno(2,modelLabel,model,kc)
          if((kc<=0).or.(kc>=3)) then
            ierror=1
            soundObservFlag=.false.
            write(ifle,*)'Unkown model equation:',trim(model)
            write(ifle,*) modnam,subnam
            return
          end if
          soundModelCode(ic)=kc
        end if
      end do
! --- Observer points setting TO HERE ----------------------------------


      end subroutine inputdata
      end module module_flow_sound
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_hpc_input
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_hpc_input)'
      integer             :: NPART,UCDFLAG,WALLFLAG
      character (len=80)  :: GRIDFIL,BCFIL,INPFIL,UCDFIL
      character (len=80)  :: WALFIL,DIRHED
      character (len=80)  :: HEADERg, HEADERb, HEADERc, HEADERw
!
!/////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata(ifli,ifll,ifle,NPETOT,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy qrgument]
!
      integer,intent(in)   :: ifli,ifll,ifle,NPETOT
      integer,intent(out)  :: ierror
!
! --- [local entities]
!
      integer :: ierr1
      integer,parameter :: hpcfil=11
!
! --- [namelist]
!
      character(80) :: mtsfil,wall,subdir
      character(80) :: hpc_boundary,hpc_comm,hpc_vertex,ucdfile
      character(80) :: hpc_source,hpc_result,hpc_anim
      character(80) :: hpc_initial,hpc_restart
      character(80)     :: hpc_p_result
      character(80)     :: hpc_p_restart
      character(80)     :: hpc_p_initial
      character(80)     :: hpc_probe        !!!modify
      character(80),parameter :: blank=' '
      integer :: NPE,walflg,ucdflg
      namelist /hpc_cntl/
     &        NPE,mtsfil,wall,subdir,
     &        hpc_vertex,hpc_boundary,hpc_comm,hpc_anim,
     &        hpc_source,hpc_result,hpc_initial,hpc_restart,
     &        walflg,ucdflg,ucdfile,
!     &        hpc_p_result,
     &        hpc_p_restart,
     &        hpc_p_initial
     &        ,hpc_probe     !!!modify
!      
!
      ierr1=0
      ierror=0
!
! --- namelist initialization
!
      NPE=1
      walflg=0
      ucdflg=0
      hpc_boundary=blank
      hpc_comm=blank
      hpc_vertex=blank
      mtsfil=blank
      wall=blank
      subdir=blank
      ucdfil=blank
      hpc_anim=blank
!
! --- variable initialization
!
      NPART=1
      UCDFLAG=0
      WALLFLAG=0

      GRIDFIL=blank
      BCFIL=blank

      INPFIL=blank
      UCDFIL=blank
      WALFIL=blank
      DIRHED=blank
      HEADERg=blank
      HEADERb=blank
      HEADERc=blank
      HEADERw=blank
!
! --- Read namelist
!
      rewind ifli
      read(ifli,hpc_cntl,iostat=ios)
      call nml_errmsg0(ifle,ios,'hpc_cntl',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      IF(NPE.NE.NPETOT) THEN
        write(ifle,*) 
     &   '*** ERR: NPE No. in CONT.HPC is NOT equal to ncpu in cort.1'
        WRITE(ifle,*) '### error : NPE must be equal to ',NPETOT
        write(ifle,*) '*** ATT: Are sure to use HPC'
        stop
      ELSEIF(NPE.LE.1) then
        WRITE(ifle,*) '### error : NPE must be greater than 1' 
        write(ifle,*) '*** ATT: Are sure to use HPC'
        stop
      ELSE
        NPART=NPE
      ENDIF
!
!
      UCDFLAG=ucdflg
      UCDFIL=ucdfile
      IF(UCDFLAG.eq.1) THEN
        if(UCDFIL.eq.blank) then
        WRITE(ifle,*) '### error : ucdfile not de defined if ucdfile=1'
          goto 9999
        endif
      ELSEIF(UCDFLAG.ne.1.and.UCDFLAG.ne.0) THEN
        WRITE(ifle,*) '### error : ucdflg must be 1 or 0'
        goto 9999
      ENDIF
!
      WALLFLAG=walflg
      WALFIL=wall
      IF(walflg.EQ.1) THEN
        IF(WALFIL.EQ.blank) THEN
        call errmsg(GRIDFIL,'vertex')
          WRITE(ifle,*) '### error :  wall not de defined if walflg=1'
          goto 9999
        ENDIF
      ELSEIF(WALLFLAG.ne.1.and.WALLFLAG.ne.0) THEN
        WRITE(ifle,*) '### error :  walflg must be 1 or 0'
        goto 9999
      ENDIF
!
      INPFIL=mtsfil
!
      call errmsg(INPFIL,'mtsfil')
      if(ierror.ne.0) goto 9999 
!
      DIRHED=subdir
      HEADERg=hpc_vertex
      HEADERb=hpc_boundary
      HEADERc=hpc_comm
      call errmsg(DIRHED,'subdir')
      if(ierror.ne.0) goto 9999 
      call errmsg(HEADERg,'hpc_vertex')
      if(ierror.ne.0) goto 9999 
      call errmsg(HEADERb,'hpc_boundary')
      if(ierror.ne.0) goto 9999 
      call errmsg(HEADERc,'hpc_comm')
      if(ierror.ne.0) goto 9999 
      HEADERw='work'
!
      return
 9999 continue
      ierror=1
      write(ifle,*) modnam,subnam
!
      contains
!
!=================================================
      subroutine errmsg(vnam1,vnam2)
!=================================================
      character(*) :: vnam1,vnam2
      if( vnam1.ne.' ' ) return
      write(*,*) '### error : data error'
      write(*,*) 'lack of data : ',vnam2
      ierror=1
      end subroutine errmsg
!
      end subroutine inputdata
!
      end module module_hpc_input

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_parameter
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
      integer :: NIFACE
      end module module_parameter
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_anim
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      integer,save             :: ianim=0
      integer,save             :: ianim_uvw,ianim_p
      integer,save             :: ianim_r,ianim_t
      integer,save,allocatable :: ianim_rans(:),ianim_comp(:)
! 
!/////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,nrans,ncomp,ierror)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle,nrans,ncomp
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      integer,parameter :: mcomp=200
      integer:: anim=0
      integer:: anim_uvw=0,anim_p=0,
     &          anim_r=0,anim_t=0,
     &          anim_rans(mcomp),anim_comp(mcomp)
      integer:: anim_GAS_WDOT(mcomp),anim_surface_SDOT(mcomp),
     &          anim_surface_Fraction(mcomp),
     &          anim_depos_etch_speed(mcomp)
      namelist/animation/
     &          anim,
     &          anim_uvw,anim_p,
     &          anim_r,anim_t,
     &          anim_rans,anim_comp,
     &          anim_GAS_WDOT,
     &          anim_surface_SDOT,
     &          anim_surface_Fraction,
     &          anim_depos_etch_speed
!
! --- [local entities]
!
      integer :: ios=0,ierr1=0
!
      allocate(ianim_rans(nrans),ianim_comp(ncomp))
!
      anim_rans(:)=0
      anim_comp(:)=0
!
      ianim=0
      ianim_uvw=0
      ianim_p=0
      ianim_r=0
      ianim_t=0
      ianim_rans(:)=0
      ianim_comp(:)=0
!
      anim_GAS_WDOT(:)=0
      anim_surface_SDOT(:)=0
      anim_surface_Fraction(:)=0
      anim_depos_etch_speed(:)=0
!
      rewind ifli
      read(ifli,animation,iostat=ios)
      if( ios.lt.0 ) return
      call nml_errmsg0(ifle,ios,'animation',ierr1)
      if(ierr1.ne.0) stop ': at module_anim'
!
      ianim=anim
      ianim_uvw=anim_uvw
      ianim_p=anim_p
      ianim_r=anim_r
      ianim_t=anim_t
      ianim_rans(1:nrans)=anim_rans(1:nrans)
      ianim_comp(1:ncomp)=anim_comp(1:ncomp)  
!
      end subroutine inputdata
      end module module_anim
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_Euler2ph
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     key-word for Euler 2 phase: ieul2ph, eul2ph
!
      implicit none
      character(21),parameter,private :: modnam='(module_Euler2ph)'
!
      integer,save :: ieul2ph=0,KE2P=0,NPHS=1
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,E2P,ierror)
!=================================================
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      logical,intent(out)  :: E2P
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      integer,parameter :: IPHMAX=10
      integer :: Euler2ph=0,N_phase=2
      character(20) :: phase_type(IPHMAX)
      
      namelist /Eul2ph/ Euler2ph,N_phase,phase_type
!
! --- [local entities]
!
      integer :: ierr1=0,ios=0
      
!
      Euler2ph=0
      N_phase=2
      rewind ifli
      read(ifli,Eul2ph,iostat=ios)
      if( ios.lt.0 ) return
      call nml_errmsg0(ifle,ios,'Eul2ph',ierr1)
      if( ios.gt.0 ) stop ': at module_Euler2ph'
!
      if(.not.(Euler2ph.eq.0.or.Euler2ph.eq.1)) then
        write(ifle,*) ' ### ERR : Euler2ph must be 0 or 1 at &Eul2ph'
        write(ifle,*) ' ### Reset [Euler2ph] in ',cntlnam
        stop ': at module_Euler2ph'
      endif
!
      if(Euler2ph==1) then
        if(N_phase<2) then
          write(ifle,'(1X,a)') 
     &   ' ### ERR : N_phase >= 2 for Euler Multi-phase flow'
          write(ifle,*) ' ### Reset [N_phase] in ',cntlnam
          stop ': at module_Euler2ph'
        else
          NPHS=N_phase-1
        endif
      else
        NPHS=1
      endif
!
      ieul2ph=Euler2ph
!
      if(ieul2ph.eq.1) then
        KE2P=1
        E2P=.true.
      else
        KE2P=0
        E2P=.false.
      endif
!
      ierror=0
      return
!
      end subroutine inputdata
!
      end module module_Euler2ph
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_vof
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     key-word for vof 2 phase: icalvof
!
      implicit none
      character(21),parameter,private :: modnam='(module_vof)'
!
      integer,save :: icalvof=0,change_phase=0
      integer,save :: LG=1,LS=2,intervof=1
      real(8),save :: coefcsf   
      integer,save :: iniflg
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!=======================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,lvof,ierror)
!=======================================================
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      logical,intent(out)  :: lvof
      integer,intent(out)  :: ierror
!
! --- [namelist]
!
      integer :: calvof=0
      integer :: interface=1,donor=1
      real*8  :: evap_h=0.d0
      character(10) :: cnv_vof,module
      real(8)  :: bldfct=0.8D0
      namelist /VOF/ calvof,interface,change_phase,evap_h,cnv_vof,
     &               module,donor,bldfct,
     &               coefcsf,iniflg
!
! --- [local entities]
!
      integer :: ierr1=0,ios=0
      
!
      calvof=0
      rewind ifli
      read(ifli,VOF,iostat=ios)
      if( ios.lt.0 ) return
      if( ios.gt.0 ) stop ': at &VOF error, stop at module_VOF'
!
      if(.not.(calvof.eq.0.or.calvof.eq.1)) then
        write(ifle,*) ' ### ERR : [calvof] must be 0 or 1 at &VOF'
        write(ifle,*) ' ### Reset [calvof] of &VOF in ',cntlnam
        stop ': at module_VOF'
      else
        write(ifll,'(a)') 
     &  'MSG : Volume Of Fraction [VOF] is defined with: '
        write(ifll,'(a)') '           VOF: F=[Liquid-Vol]/[Cell-Vol]'
      endif
!
      if(.not.(interface.eq.1.or.interface.eq.2)) then
        write(ifle,*) ' ### ERR : [interface] must be 1 or 2 !'
        write(ifle,*) ' ### MSG : interface=1 : liquid-gas interface'
        write(ifle,*) ' ### MSG : interface=2 : liquid-solid interface'
        stop ': at module_VOF'
      elseif(interface.eq.1) then
        write(ifle,*) ' ### MSG : VOF Used for liquid-gas interface'
        write(ifle,*) ' ### MSG : interface=1 : liquid-gas interface'
        write(ifle,*) ' ### MSG : interface=2 : liquid-solid interface'
      elseif(interface.eq.2) then
        write(ifle,*) ' ### MSG : VOF Used for liquid-solid interface'
        write(ifle,*) ' ### MSG : interface=1 : liquid-gas interface'
        write(ifle,*) ' ### MSG : interface=2 : liquid-solid interface'
      endif
!
      if(.not.(change_phase.eq.0.or.change_phase.eq.1)) then
      write(ifle,*) 
     &   ' ### ERR : [change_phase] must be 0 or 1 at &VOF'
      write(ifle,*) ' ### MSG : change_phase=0 : NO phase changing'
      write(ifle,*) ' ### MSG : change_phase=1 : Phase changing'
      stop ': at module_VOF'
      endif
!
      icalvof=calvof
      intervof=interface
!
      if(icalvof.eq.1) then
        lvof=.true.
      else
        lvof=.false.
      endif
!
      ierror=0
      return
!
      end subroutine inputdata
!
      end module module_VOF

!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_scalar
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
      character(21),parameter,private :: modnam='(module_scalar)'
!
      character(20),save,allocatable :: sclname(:)
      character(20),save,allocatable :: potname(:)
      integer,allocatable :: scl(:)
      integer :: icalke(2),icalrsm(7),icalvof(1)
      integer :: KRANS=0
      logical :: calgeq
      logical :: calxi
      logical :: caltemp
      logical :: calh
!
      integer,save :: idifflm=0
      integer,parameter,public :: ORG=0,consF=1,Blog=2,BlogM=3,
     &                            FltRans=4
!
      integer,save :: ical_POTN,Kpotn=0
!//////////////////////////////////////////////////////////////////////
      contains
!======================================================================
      subroutine inputdata
     & (ifli,ifll,ifle,cntlnam,
     &  lcomp,lrans,E2P,kemdl,lowke,lRNG,lCHEN,l2SKE,lSDES,lKLES,
     &  RSMmdl,RANSMDL,lvof,lFC,npotn,nrans,NPHS,ierror)
!======================================================================
      use module_model,only : u_func
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: ifli,ifll,ifle
      character(*),intent(in)  :: cntlnam
      integer,intent(inout) :: nrans,NPHS,lFC
      logical,intent(in)    :: E2P,lcomp,kemdl,RSMmdl,lvof,lowke,
     &                         lRNG,lCHEN,lSDES,RANSMDL,l2SKE,lKLES
      logical,intent(inout) :: lrans
      integer,intent(out)   :: ierror
      integer,intent(out)   :: npotn
!
! --- IF(ICAL(ISCA(ISC)))  CALL CALCSC(ISCA(ISC)
!
! --- [namelist]
!
      integer,parameter :: mrans=50
      integer,parameter :: mpotn=10
      character(20) :: scalar_name(mrans)
      character(20) :: Potential_name(mpotn)
      character(20) :: Geq,Xieq,Tempeq,heq
      character(80) :: flamelet_file,Lam_spd_file
      integer,save :: tur_vel,date_base=0,igT_iter,DENformat
!
      integer,parameter :: lc=10
      character(lc) :: flm_model
!
      character(lc),parameter ::
     &   ORG_les   =  'ORG_les   ',
     &   les_const =  'les_const ',
     &   C_les     =  'C_les     ',
     &   H_les     =  'H_les     ',
     &   rans_const=  'rans_const'
!
      character(lc),parameter ::
     &    flmknd(5)=(
     &   /ORG_les,les_const,C_les,H_les,rans_const/)
! --- Partial-pre-mixed:            -------------------------------
      real*8  :: Xi0,              ! Xi (constant value)
     &           fuel_rho,         ! Fuel density in unburn
     &           fuel_temp,        ! Fuel temperature in unburn
     &           Oxidant_rho,      ! Oxidant density in unburn
     &           Oxidant_temp,     ! Oxidant temperature in unburn
     &           Xi_max            ! Max Xi for Xi-Density profile

!
      namelist /scalar/ scalar_name  !,Geq,Xieq,Tempeq,heq
!
      namelist /flamelet/ Geq,Xieq,Tempeq,heq,tur_vel,
     &                    flamelet_file,Lam_spd_file,
     &           Xi0,              ! Xi (constant value)
     &           fuel_rho,         ! Fuel density in unburn
     &           fuel_temp,        ! Fuel temperature in unburn
     &           Oxidant_rho,      ! Oxidant density in unburn
     &           Oxidant_temp,     ! Oxidant temperature in unburn
     &           Xi_max,            ! Max Xi for Xi-Density profile
     &           date_base,
     &           flm_model,igT_iter,DENformat
! -----------------------------------------------------------------
!
      integer :: Cavi=0,Iter_cavi,Iter_poten
      real*8  :: T_ref=288.15,Pc=1944.61d6,T_sat,CeA_P,CeA_M
      integer :: Axis=1,Boundary_NO=-1,FSTR_OUT=0
      real*8  :: Power=1.d0
      namelist /Cavitation/ Cavi,T_ref,Pc,T_sat,
     &                      Iter_cavi,CeA_P,CeA_M,
     &                      Axis,Boundary_NO,FSTR_OUT,Power
!
      namelist /Potential/ Potential_name,Iter_poten
!
! --- [local entities]
!
      integer :: i_user,nransno=0,no=0,ios=0,ierr=0,nmodel=mrans
      integer :: icalgeq,icalxi,icaltemp=0,icalh,icalmdl
      integer :: i
      character(3),parameter :: flmltlst(2)=(/'no ','yes'/)
      character(3) :: text
      character(len=1) :: recal
!
      scalar_name(:)=' '
      calgeq  = .false.
      calxi   = .false.
      caltemp = .false.
      calh    = .false.
      Geq     = 'no '
      Xieq    = 'no '
      Tempeq  = 'no '
      heq     = 'no '
      date_base=0
      flm_model='ORG_les'
!
      rewind ifli
      read(ifli,scalar,iostat=ios)
      if(ios.gt.0) stop ': at reading error at &scalar'

      rewind ifli
      read(ifli,flamelet,iostat=ios)
      if(ios.gt.0) stop ': at reading error at &flamelet'

      rewind ifli
      read(ifli,Cavitation,iostat=ios)
      if(ios.gt.0) stop 'reading error at &Cavitation'
!
      call nml_listno(2,flmltlst,Geq,icalgeq)
      call nml_chkchr0(ifle,'Geq',Geq,icalgeq,ierror)
      if(ierror.ne.0) then
        stop 'reading error at &scalar about Geq'
      endif
      icalgeq=icalgeq-1
!
      call nml_listno(2,flmltlst,Xieq,icalxi)
      call nml_chkchr0(ifle,'Xieq',Xieq,icalxi,ierror)
      if(ierror.ne.0) then
        stop 'reading error at &scalar about Xieq'
      endif
      icalxi=icalxi-1
!
      call nml_listno(2,flmltlst,heq,icalh)
      call nml_chkchr0(ifle,'heq',heq,icalh,ierror)
      if(ierror.ne.0) then
        stop 'reading error at &scalar about heq'
      endif
      icalh=icalh-1
!
      icalmdl=0
      call nml_listno(5,flmknd,flm_model,icalmdl)
      call nml_chkchr0(ifle,'flm_model',flm_model,icalmdl,ierror)
      if(ierror.ne.0) then
        stop'reading error at &flamelet/flm_model'
      endif
      icalmdl=icalmdl-1
      
!
      if(icalgeq .eq.1) calgeq  = .true.
      if(icalxi  .eq.1) calxi   = .true.
      if(icaltemp.eq.1) caltemp = .true.
      if(icalh   .eq.1) calh    = .true.
!
      if(calxi) then
        idifflm=icalmdl
        if(idifflm==ORG) then  !=ORG
          write(ifll,'(1x,a)') 
     &   'MSG: ORG LES flamelet model is used for diffusion flame'
          if(date_base/=0) then
            stop 'ERR: date_base=0 for "ORG_les"'
          endif
        elseif(idifflm==consF) then  !=Blog
          write(ifll,'(1x,a)') 
     &   'Constant LES flamelet model is used for diffusion flame'
          if(date_base/=1) then
            stop'ERR: date_base=1 for "les_const"'
          endif
        elseif(idifflm==Blog) then  !=Blog
          write(ifll,'(1x,a)') 
     &   'Blog LES flamelet model is used for diffusion flame'
          if(date_base/=2) then
            stop'ERR: date_base=2 for "C_les"'
          endif
        elseif(idifflm==BlogM) then  !=BlogM
          write(ifll,'(1x,a)') 
     &   'Modified Blog LES flamelet model is used for diffusion flame'
          if(date_base/=2) then
            stop'ERR: date_base=2 for "H_les"'
          endif
        elseif(idifflm==FltRans) then  !=FltRans
          write(ifll,'(1x,a)') 
     &   'Constant RANS flamelet model is used for diffusion flame'
          if(date_base/=1) then
            stop'ERR: date_base=1 for "rans_const"'
          endif
        else
          write(ifle,'(1x,a)') 
     &     'MSG: flm_model="ORG_les": Defult (ORG) LES flamelet model'
          write(ifle,'(1x,a)') 
     &     'MSG: flm_model="les_const": Constant LES flamelet model'
          write(ifle,'(1x,a)') 
     &     'MSG: flm_model="C_les": Blog LES flamelet model'
          write(ifle,'(1x,a)') 
     &     'MSG: flm_model="H_les": Modified Blog LES flamelet model'
          write(ifle,'(1x,a)') 
     &     'MSG: flm_model="rans_const": Constant RANS flamelet model'
          write(ifle,'(1x,a)') 
     &     'MSG: Check &flamelet/flm_model '
          stop 'ERR:  &flamelet/flm_model'
        endif

        if(idifflm/=ORG) then
          if(date_base==0) then
             write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=1 for [chemtable_CH4_steady]'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=2 for [chemtable_CH4_unsteady]'
             stop'ERR: date_base=1 or date_base=2'
          elseif(date_base==1) then
            write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base'
            write(ifle,'(1x,a)') 
     &      'MSG: date_base=1 for [chemtable_CH4_steady]'
          elseif(date_base==2) then
            write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base'
            write(ifle,'(1x,a)') 
     &      'MSG: date_base=2 for [chemtable_CH4_unsteady]'
          else
            write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=1 for [chemtable_CH4_steady]'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=2 for [chemtable_CH4_unsteady]'
            stop'ERR: date_base=1 or date_base=2'
          endif
        else  !if(idifflm==ORG) then
          if(date_base/=0) then
            date_base=0
            stop'ERR: date_base=0 for flm_model="ORG_les"'
          endif
        endif
      endif
!
      if(calgeq.and.calxi.and.idifflm/=ORG) then
        write(ifle,'(1x,2a)') 
     &'MSG: ONLY flm_model="ORG_les" ',
     &    'support Partial pre-mixed flame model'
        stop'ERR: set [&flamelet/flm_model="ORG_les"]'
      elseif(calgeq.and.idifflm/=ORG) then
        write(ifle,'(1x,a)') 
     &   'MSG: ONLY flm_model="ORG_les" support pre-mixed flame model'
        stop'ERR: set [&flamelet/flm_model]'
      endif
!
      if(calxi) then
      if(.NOT.RANSMDL.and.idifflm==FltRans) then
        write(ifle,'(1x,a)') 
     &   'MSG: flm_model="rans_const": Constant RANS flamelet model'
        stop 'ERR: RANS model MUST be set for flm_model="rans_const"'
      elseif(RANSMDL.and.(idifflm/=FltRans.and.idifflm/=ORG)) then
        write(ifle,'(1x,a)') 
     &   'MSG: flm_model="rans_const": Constant RANS flamelet model'
        stop 'ERR: RANS model ONLY support flm_model="rans_const"'
      endif
      endif
!
      if(Cavi==1) then
        if(.not.lcomp) 
     &  stop 'ERR:Cavitation model must use [comp]'
      elseif(.not.(Cavi==0.or.Cavi==1)) then
        stop 'ERR: [Cavi] must be 0 or 1'
      endif
!-----------------------------------
! --- count scalar number 
!-----------------------------------
      no=0
!
      if(kemdl.or.lowke.or.lRNG.or.lCHEN.or.l2SKE) then
        no=no+2
      endif
!
      if(lSDES) then
        no=no+1
      endif
!
      if(RSMmdl) then
        no=no+7
      endif
!
      if(E2P) then
        no=no+(NPHS+1)
      endif
!
      if(lvof) then
        no=no+1
      endif
!
      if(lFC>0)then
        no=no+1
      endif
!
      if(Cavi==1) then
        no=no+2
      endif
!
      if(lKLES) then
        no=no+1
      endif
!
      if(calgeq ) no=no+1
      if(calxi  ) no=no+1
      if(caltemp) no=no+1
      if(calh   ) no=no+1
      if(calxi) then
        if(idifflm==ORG) then
          no=no+2   ! ixi_2,ixi_X
        elseif(idifflm==consF) then
          no=no+2   ! ixi_2,ixi_X
        elseif(idifflm==Blog) then
          no=no+4   ! ixi_2,ixi_X,ixi_C
        elseif(idifflm==BlogM) then
          no=no+4   ! ixi_2,ixi_X,ixi_C
        elseif(idifflm==FltRans) then
          no=no+2   ! ixi_2,ixi_X
        endif
      endif
!
      do i_user=1,mrans
      if(scalar_name(i_user).ne.' ') then
        no=no+1
      endif
      enddo
!
      if(u_func(1)==1) then
        no=no+1
      endif
!
      nransno=no
!
! --- allocating array
!
      allocate(sclname(nransno),scl(nransno),stat=ierr)
      if(ierr.ne.0) stop ': allocate error for [sclname]'  
      sclname(:)=' '
      scl(:)=0
!-----------------------------------
! --- repeat accont [no] and name
!-----------------------------------
      no=0
!
      nmodel=mrans
      if(kemdl.or.lowke.or.lRNG.or.lCHEN.or.l2SKE) then
        no=no+1
        sclname(no)='RANS_K'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='RANS_eps'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      IF(lSDES) then
        no=no+1
        sclname(no)='S-A_DES_Nu'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(RSMmdl) then
        no=no+1
        sclname(no)='RANS_eps'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='uu'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='vv'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='ww'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='uv'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='uw'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='vw'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(E2P) then
!        allocate(iaph(NPHS+1))
        do I=1,NPHS+1
        write(text,'(i3)') I
        text=adjustl(text)
        no=no+1
        sclname(no)='alph'//TRIM(text)
        nmodel=nmodel+1
        scl(no)=nmodel
        enddo
!
!        no=no+1
!        sclname(no)='alph2'
!        nmodel=nmodel+1
!        scl(no)=nmodel
      endif
!
      if(lvof) then
        no=no+1
        sclname(no)='VOF'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(lFC>0) then
        no=no+1
        sclname(no)='Saturation'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(Cavi==1) then
        no=no+1
        sclname(no)='Cavi_Y'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        no=no+1
        sclname(no)='VOID'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(lKLES) then
        no=no+1
        sclname(no)='LES_K'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(calgeq) then
        no=no+1
        sclname(no)='G'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(calxi) then
        no=no+1
        sclname(no)='Xi'
        nmodel=nmodel+1
        scl(no)=nmodel
!
        if(idifflm==ORG) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
!
          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
        elseif(idifflm==consF) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
!
          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
        elseif(idifflm==Blog) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel

          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel

          no=no+1
          sclname(no)='REAC_DEGRESS'
          nmodel=nmodel+1
          scl(no)=nmodel

          no=no+1
          sclname(no)='WDOT_C'
          nmodel=nmodel+1
          scl(no)=nmodel
        elseif(idifflm==BlogM) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel

          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel

          no=no+1
          sclname(no)='REAC_DEGRESS'
          nmodel=nmodel+1
          scl(no)=nmodel

          no=no+1
          sclname(no)='WDOT_C'
          nmodel=nmodel+1
          scl(no)=nmodel
        elseif(idifflm==FltRans) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
!
          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
        endif

      endif
!
      if(caltemp) then
        no=no+1
        sclname(no)='Temp'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      if(calh) then
        no=no+1
        sclname(no)='h_flmlt'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!--------------------------------------------------------
! --- User defined scalar must be put end of system model
!--------------------------------------------------------
      do i_user=1,mrans
      if(scalar_name(i_user).ne.' ') then
        no=no+1
        scl(no)=i_user
        sclname(no)=scalar_name(i_user)
      endif
      enddo
!
      if(u_func(1)==1) then
        no=no+1
        sclname(no)='c_dash^2'
        nmodel=nmodel+1
        scl(no)=nmodel
      endif
!
      write(ifll,2400) 
      write(ifll,2100)
      write(ifll,2150)
      write(ifll,2200)
      write(ifll,2000) 
      do no=1,nransno
        write(ifll,2300) no,TRIM(adjustl(sclname(no)))
        write(ifll,2400) 
 2000   format(1x,50('-'))
 2100   format(1x,4('#'),4x,'Scalar Information',20x,4('#'))
 2150   format(1x,4('#'),42x,4('#'))
 2200   format(1x,4('#'),4x,'Scalar No.',13x,'Scalar Name',4x,4('#'))
 2300   format(1x,4('#'),4x,I8.3,a20,10x,4('#'))
 2400   format(1x,50('#'))
      enddo
!
      nrans=nransno
      KRANS=0
      if(nrans.gt.0) then
        KRANS=1
      endif
      lrans=.false.
      if(nransno.gt.0) lrans=.true.
      if(nrans.gt.mrans) then
        write(ifll,'(a)') 'ERR : Salar number cannot great than ',
     &   mrans
        stop ': at module_scalar/inputdata'
      endif
!
!
!--------------------
! --- Potential scalar
!--------------------
      Potential_name(:)=' '
      rewind ifli
      read(ifli,Potential,iostat=ios)
      if(ios.gt.0) stop ': at reading error at &Potential'
!
      no=0
      if(lFC>0) then
        no=no+2
      endif
!-----------------
      do i_user=1,mpotn
      if(Potential_name(i_user).ne.' ') then
        no=no+1
      endif
      enddo
!
      npotn=no
!
      allocate(potname(npotn),stat=ierr)
      if(ierr.ne.0) stop ': allocate error for [potname]'  
      potname(:)=' '
!
!
!
!
      no=0
      if(lFC>0) then
        no=no+1
        potname(no)='FAI_S'
        no=no+1
        potname(no)='FAI_M'
      endif
!
      
!--------------------------------------------------------
! --- User defined scalar must be put end of system model
!--------------------------------------------------------
      do i_user=1,mpotn
      if(Potential_name(i_user).ne.' ') then
        no=no+1
        potname(no)=Potential_name(i_user)
      endif
      enddo
!
      write(ifll,*) 
      write(ifll,3400) 
      write(ifll,3100)
      write(ifll,3150)
      write(ifll,3200)
      write(ifll,3000) 
      do no=1,npotn
        write(ifll,3300) no,TRIM(adjustl(potname(no)))
        write(ifll,3400) 
 3000   format(1x,50('-'))
 3100   format(1x,4('#'),4x,'Potential Scalar Information',10x,4('#'))
 3150   format(1x,4('#'),42x,4('#'))
 3200   format(1x,4('#'),4x,
     &  'Potential No.',7x,'Potential Name',4x,4('#'))
 3300   format(1x,4('#'),4x,I8.3,a20,10x,4('#'))
 3400   format(1x,50('#'))
      enddo
!
      if(npotn/=0) Kpotn=1
!
      if(nransno/=0.or.npotn/=0) then
 1000 continue
      write(ifll,'(1X,a,2x,a)') 'MSG: ATTENTION Above Massage /',
     &   'type (y) & return'
      read(*,*) recal
      if(.not.((recal=='y'.or.recal=='Y'))) goto 1000
      endif
!
      ierror=0
      return
!
      end subroutine inputdata
      end module module_scalar
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_FFRGFwords
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character*80,parameter :: asciiv2_FFGF="#A_GF_V2"
      character*80,parameter :: unformv2_FFGF="#U_GF_V2"
      character*80,parameter :: binaryv2_FFGF="#B_GF_V2"
      character*80,parameter :: newSet_FFGF="#NEW_SET"
      character*80,parameter :: fileend_FFGF="#ENDFILE"
      character*80,parameter :: customData_FFGF="#CUSTOM_DATA"
      integer,parameter :: unkownDatasize_FFGF=0
!     For grids.
      character*80,parameter :: ffrGrid_FFGF="FFR_GRID_FIELD"
      character*80,parameter :: gridHead_FFGF="#FFR_GRID_HEADER"
      character*80,parameter :: gridHeade_FFGF="#FFR_GRID_HEADER_END"
      character*80,parameter :: elemType_FFGF="#ELEM_TYPE_LIST"
      character*80,parameter :: elemTypee_FFGF="#ELEM_TYPE_LIST_END"
      character*80,parameter :: bndryType_FFGF="#BNDRY_TYPE_LIST"
      character*80,parameter :: bndryTypee_FFGF="#BNDRY_TYPE_LIST_END"
      character*80,parameter :: vrtxCord_FFGF="#VRTX_CORD"
      character*80,parameter :: vrtxCorde_FFGF="#VRTX_CORD_END"
      character*80,parameter :: bndryFace_FFGF="#BNDRY_FACE"
      character*80,parameter :: bndryFacee_FFGF="#BNDRY_FACE_END"
      character*80,parameter :: elemVrtx_FFGF="#ELEM_VRTX"
      character*80,parameter :: elemVrtxe_FFGF="#ELEM_VRTX_END"
      character*80,parameter :: undefMatName="MAT"
      integer,parameter :: undefBndryType=0
      integer,parameter :: defaultCellTypeID=1
!     For fulid force.
      character*80,parameter :: ffrForce_FFGF="FFR_FORCE_FIELD"
      character*80,parameter :: forceHead_FFGF="#FFR_FORCE_HEADER"
      character*80,parameter :: forceHeade_FFGF="#FFR_FORCE_HEADER_END"
      character*80,parameter :: forceWall_FFGF="#FFR_FORCE_WALLS"
      character*80,parameter :: forceWalle_FFGF="#FFR_FORCE_WALLS_END"
      character*80,parameter :: forceModel_FFGF="#FFR_FORCE_MODEL"
      character*80,parameter :: forceModele_FFGF="#FFR_FORCE_MODEL_END"
      character*80,parameter :: forceObserv_FFGF="#FFR_FORCE_OBSERVER"
      character*80,parameter :: forceObserve_FFGF=
     &                                    "#FFR_FORCE_OBSERVER_END"
      character*80,parameter :: forceData_FFGF="#FFR_FORCE_DATA"
      character*80,parameter :: forceDatae_FFGF="#FFR_FORCE_DATA_END"


      end module module_FFRGFwords
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!Jiang Yuyan, 2005/10/17
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_rad
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_rad)'
!
!
!< module for radiation >
!
      character*10 :: radmodel='MC'  !'ZONE'/'MC'/'FVM'
      integer      :: radflag=0      !do or do not
      integer      :: radgpopt=1				
!use or use not 2-scale grids for fine mesh
      REAL*8       :: vfeps=1.0d-6
      REAL*8       :: feps=1.0d-6
      REAL*8       :: RLX=0.75d0
      INTEGER      :: DOMSN=4
      INTEGER      :: iModHF=1
      INTEGER      :: NRAY=1000
      INTEGER      :: RadMG=5000
      INTEGER      :: TeaRst=1

!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out)  :: ierror
!
! [namelist]
!
      character*10 :: rmodel,gasmodel
      integer      :: rflag,mgpnum,radjump,ngauss
      real*8       :: charlen
      namelist /radoption/ rflag,rmodel,radgpopt,mgpnum,domsn,
     &	feps,vfeps,RLX,nray,tearst,radjump,gasmodel,ngauss,charlen
!
! [local entities]
!
      real*8,parameter :: undef=-huge(1.d0)
      integer          :: iset=0
      integer          :: ios,ierr1
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      rewind ifli
      read(ifli,radoption,iostat=ios)
      if( ios.lt.0 ) return
      call nml_errmsg0(ifle,ios,'radoption',ierr1)
      if( ierr1.ne.0 ) goto 9999
      
      radflag=rflag
!
 1001 IF(radflag.eq.1) then
	if(rmodel(1:4).NE.'ZONE'.AND.
     &	   rmodel(1:2).NE.'MC'.AND.
     &     rmodel(1:3).NE.'FVM') then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'model = ',trim(rmodel)
          write(ifle,*) 'it must be <ZONE> <MC> <FVM>'
          goto 9999
        endif
      endif
! 
      radmodel = rmodel
      if(mgpnum.GT.0) RadMG=mgpnum
!
      if(radmodel(1:4).NE.'ZONE') iteaflag=0
      IF(iteaflag.EQ.1.OR.radmodel(1:4).NE.'ZONE') THEN
	iModHF=1
      ELSE
	iModHF=2
      ENDIF
	  
!
      return
!
 9999 continue

      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_rad
!
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_radsxf
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      INTEGER::	NRD,MRD,NRD0,MRD0,MID
      INTEGER:: NPALL,NDatExpt,NumTran,NumTranIn,NMatBoun
      CHARACTER*20 :: StackOpt
      REAL*8 :: AveVs(10)

	
      INTEGER,ALLOCATABLE :: FaceNum(:),NNPE(:),NNPS(:,:),
     &			   CNTIVITY(:,:,:),ElemKind(:)
	

      REAL*8, ALLOCATABLE  ::  Albedo(:),AreaVol(:),PROP(:),scabeta(:)
      INTEGER,ALLOCATABLE::  ETYPE(:)

      REAL*8, ALLOCATABLE::  CoefToRD(:),RaToExpt(:)
      INTEGER,ALLOCATABLE::  IndexToRD(:,:),IndexToExpt(:,:),
     &					   MatBounEnd(:),ExptID(:)
      REAL*8, ALLOCATABLE::  PRATIO(:),RaNeb(:)

!------ GAS
      REAL*8, ALLOCATABLE::  GAbsorb(:),Volume(:)
	
!------ WALL
      INTEGER,ALLOCATABLE::   WElem(:,:),WNEB(:),NodeBounID(:)
      REAL*8, ALLOCATABLE::	WGrayDeg(:),WArea(:),QW(:)
      REAL*8,ALLOCATABLE::    WABCD(:,:),WNV(:,:)
!for benchmark
      INTEGER,ALLOCATABLE::   WDAT(:,:),WPNUM(:,:)
      INTEGER :: NWDat
!------ PARTICLE
      REAL*8, ALLOCATABLE::	PAbsorb(:), PScatter(:)
!------ Ray
      INTEGER, ALLOCATABLE:: ReflectTime(:)
      REAL*8, ALLOCATABLE::  FinalPos(:,:)

!TEMPERATURE
      INTEGER,ALLOCATABLE::  VIndex(:),RevIndex(:)
c	REAL*8, ALLOCATABLE::  TT(:),TPtc(:),T(:),DIVQ(:),BB(:)

!------ Exchange-Area (Configure-factor)
      INTEGER, ALLOCATABLE:: RDId(:), RDN1(:),RDN2(:),RDIndex(:)
      REAL*8, ALLOCATABLE::  RDValue(:),RDV2(:)
!RDN1	 :   off-diagonal-line elements in Lower-part
!RDN2	 :   all off-diagonal-line elements


!------ BUCKET MESH
      INTEGER, ALLOCATABLE:: GBlockNum(:,:,:),GBlockIndex(:,:,:)
      INTEGER, ALLOCATABLE:: GBlockID(:)
					!HEXADEDRAL BUCKET
      INTEGER, ALLOCATABLE:: WBlockNum(:,:,:),WBlockIndex(:,:,:)
      INTEGER, ALLOCATABLE:: WBlockID(:)
      REAL*8, ALLOCATABLE::  RIN(:,:)
!SOLID ANGLE BUCKET IN A GLOBE SURFACE
      
!------ DUMMY ARRAYS
      INTEGER, ALLOCATABLE:: WK1(:),WK2(:),WK3(:),WK4(:),WK5(:),WK6(:)
      INTEGER, ALLOCATABLE:: WKMD1(:,:),WKMD2(:,:),WK7(:)
      REAL*8, ALLOCATABLE::  VWK1(:),VWK2(:),VWK3(:),VWKMD(:,:)
      REAL*8, ALLOCATABLE::  VWK3D(:,:,:),CTR(:,:)


      end module module_radsxf
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_radslv
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
	
      real*8,allocatable:: W1K1 (:)
      real*8,allocatable:: W1K2 (:)
      real*8,allocatable:: W1K3 (:)
      real*8,allocatable:: W1K4 (:)
      real*8,allocatable:: W1K5 (:)
      real*8,allocatable:: W1K6 (:)
      real*8,allocatable:: W1K7 (:)
      real*8,allocatable:: W1K8 (:)
      real*8,allocatable:: W1K9 (:)
      real*8,allocatable:: W1K10(:)
      real*8,allocatable:: W1K11(:)
      real*8,allocatable:: W1K12(:)

      INTEGER,allocatable:: WKN1(:)
      INTEGER,allocatable:: WKN2(:)
      INTEGER,allocatable:: WKN3(:)
      INTEGER,allocatable:: WKN4(:)
      INTEGER,allocatable:: WMK(:,:)
      end module module_radslv
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_radgroup
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	
      INTEGER:: IFlagGroup,NumGPTran
      real*8,allocatable:: GpProp(:)
      INTEGER,allocatable:: GpCNum(:),WSID(:,:)
!      
      REAL*8, ALLOCATABLE::  COEF(:)
      INTEGER,ALLOCATABLE::  NIndexToGp(:,:)
!      
      end module module_radgroup
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_radwork
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	
      real*8,allocatable  :: VWK1(:),VWK2(:),VWK3(:),VWK4(:)
      real*8,allocatable  :: VWKM1(:,:),VWKM2(:,:)
      real*8,allocatable  :: VWK3D(:,:,:)

      INTEGER,allocatable :: WK1(:),WK2(:),WK3(:),WK4(:),WK5(:),WK6(:)
      INTEGER,allocatable :: NWK1(:),NWK2(:),NWK3(:)
      INTEGER,allocatable :: WKM1(:,:),WKM2(:,:),WKM3(:,:)
      end module module_radwork
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

