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
!      module module_vector
!      module module_movegrid 
!      module module_nowtime  
!      module module_output
!      module module_param
!      module module_rans
!      module module_simple 
!      module module_source 
!      module module_time 
!      module module_turbparm 
!      module module_species 
!      module module_chemreac 
!      subroutine modutl_chkset 
!      subroutine modutl_setnow 
!      module module_gf 
!      module module_usersub 
!      module module_hpc 
!      module module_flow_sound 
!      module module_constant 
!      module module_dimension 
!      module module_metrix 
!      module module_anim 
!      module module_Euler2ph 
!      module module_vof 
!      module module_FFRGFwords 
!      module module_cdcl
!      module module_probe
!      module module_scalar 
!      module module_particle 
!      module_mphase_pc
!      module_rad
!      module_real_gas
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_boundary
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(LEN=21),parameter,private ::modnam='(module_boundary)'
!
!< module for boundary conditions & domains >
!
!-< 1. boundary conditions >-
!
!--< 1.1 flags >--
!
      integer,parameter :: kxnone=-999999
!
! kxnone : value of kdbcnd which indicate undefined (Wall BC) 
!
      integer,parameter :: kdprdc=-1
      integer,parameter :: kdsymm=-2
      integer,parameter :: kdintr=-3
      integer,parameter :: kdilet=-4
      integer,parameter :: kdolet=-5
      integer,parameter :: kdtchi=-6
      integer,parameter :: kdfire=-7
      integer,parameter :: kdpres=-8
      integer,parameter :: kdstag=-9
      integer,parameter :: kdbuff=-10
      integer,parameter :: kdsld=-11
      integer,parameter :: kdcvd=-12
      integer,parameter :: kdshutr=-13
      integer,parameter :: kdpors=-14
      integer,parameter :: kdovst=-15
!
      integer,parameter :: kdnature=-1
      integer,parameter :: kdinterf=-2
      integer,parameter :: kddefine=-3
      integer,parameter :: kdperiod=-4
!
! kdprdc : value of kdbcnd(0,:) for periodic boundary
! kdsymm : value of kdbcnd(0,:) for symmetric boundary
! kdintr : value of kdbcnd(0,:) for interface boundary   (wall)
! kdilet : value of kdbcnd(0,:) for inlet boundary       (inlet)
! kdolet : value of kdbcnd(0,:) for outlet boundary      (outlet)
! kdtchi : value of kdbcnd(0,:) for touch inlet boundary (inlet)
! kdfire : value of kdbcnd(0,:) for fire wall boundary   (wall)
! kdpres : value of kdbcnd(0,:) for pressure boundary
!                                  (absolute pressure)   (pressure)
! kdstag : value of kdbcnd(0,:) for stagnation pressure boundary
!                                  (absolute pressure)   (pressure)
! kdbuff : buffle wall BC                                (wall)
! kdsld  : sliding BC
! kdcvd  : value of kdbcnd(0,:) for CVD surface reaction 
!                                        wall boundary   (wall)
! kxnone : wall                                          (wall)
!
      logical,save :: lsldbc=.false.,lovstbc=.false.
!
!
      integer,parameter :: kvnslp=1, kvfslp=2, kvlglw=3,
     &                     kvrogf=4, kvmodl=5,kvEWF=6,kvsataW=7
      integer,parameter :: ktdirc=1, ktneum=2, kttrns=3,ktEWF=4
      integer,parameter :: kydirc=1, kyneum=2, kytrns=3,kyEWF=4
      integer,parameter :: kkdirc=1, kkneum=2, kklglw=3
      integer,parameter :: kmdirc=1, kmneum=2!, kmintr=3
      integer,parameter :: kpdirc=1, kpneum=2,kpintr=3,kpindr=4
      integer,parameter :: kpart_s=1,kpart_r=2,kpart_d=3
!
! kv* : value of kdbcnd(1,:) for velocity
!     : kvnslp / no-slip                        (=1)
!     : kvfslp / free-slip                      (=2)
!     : kvlglw / log-law                        (=3)
!     : kvrogf / rough-wall
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
! km* : value of kdbcnd(5,:) for MHD:A
!     : kmdirc / 
!     : kmneum / 
!     : value of kdbcnd(6,:) for MHD:FAI
!     : kmdirc / 
!     : kmneum / 
!     : value of kdbcnd(7,:) for MHD: BC type
!     : 

!--< 1.2 values >--
!
      integer         :: ical_buff=0,ical_poroBC=0
! ical_poroBC=0 : no porous jump
! ical_poroBC=1 : porous jump (only for buffle & outlet)
!                 Darcy's law
!                 -(mu/alpha*V+C2*rho*|V|V)
! poro(nb)      : =0 : no porous jump
!                 =1 : porous jump
! 
! surfreac(nb)  : =0 : no surface reaction  
!                 =1 : have surface reaction 
!                     (firewall,wall,interface,cvdsurface) 
!                     ical_suf=1
!                 =2 : PEFC (buffle)
!
      integer,private :: ncomp=0 , nrans=0
      integer         :: iscomp=0, israns=0,distrb=0
!
! ncomp  : no. of Ys(i) (chemical species)
! nrans  : no. of AKs(i) (dependent variables of RANS model)
! iscomp : start index of wdbcnd for Ys(i)
! israns : start index of wdbcnd for AKs(i)
      integer,parameter :: mpotn=10
!
      integer         :: nbcnd =0
      integer,private :: nbvals=0
      integer,save,allocatable          :: nobcnd (:)
!      integer,save,allocatable          :: bcndno (:)
      integer,save,allocatable          :: boundIDMap(:,:)
!      CHARACTER*80,allocatable          :: SFBOUN(:)
      integer,save,allocatable          :: kdbcnd (:,:)
      integer,save,allocatable          :: openout(:)
      integer,save,allocatable          :: poro(:)
      real(8),save,allocatable          :: C2pBC(:),alpha_coBC(:),
     &                                     thkposBC(:),prosty(:)

      integer,save,allocatable          :: rotwall(:)
      real(8) ,save,allocatable         :: endw(:,:),beginw(:,:)
      real(8) ,save,allocatable         :: omega(:)
      real(8) ,save,allocatable         :: strtimw(:)
      real(8) ,save,allocatable         :: rotupw(:)
      integer,save,allocatable,private  :: iwbcnd (:)
      real(8) ,save,allocatable         :: wdbcnd(:,:)
      real(8) ,save,allocatable         :: MHD_bnd(:,:,:)
      real(8) ,save,allocatable         :: wdbcnd2(:,:)
      real(8) ,save,allocatable,private :: wdtbcnd(:,:)
      real(8) ,save,allocatable,private :: wdtbcnd2(:,:)
      real(8) ,save,allocatable         :: adbcnd(:,:)
      real(8) ,save,allocatable         :: adbcnd2(:,:)
      integer, save,allocatable         :: LBC_INDEX(:)
      integer, save,allocatable         :: LBC_pair(:)
      integer, save,allocatable         :: MAT_BCIDX(:,:)
      real(8) ,save,allocatable         :: stagvel(:)
      real(8) ,save,allocatable         :: stagare(:)
      integer, save,allocatable         :: idis(:)
      integer ,save,allocatable         :: nblade(:,:)
      integer, save,allocatable         :: rotsld(:)
      real(8) ,save,allocatable         :: sizein(:)
      integer, save,allocatable         :: kdpotn(:,:)
      real(8) ,save,allocatable         :: vapotn(:,:)
      integer, save,allocatable         :: mem_cat(:)
!
      real*8 , save,allocatable         :: radprop(:,:)
      integer, save,allocatable         :: radwalltype(:)
!
      real(8) ,save,allocatable         :: vofang(:)
!
      real(8) ,save,allocatable         :: shut_strt(:)
      logical ,save,allocatable         :: pkdwall(:)
!
! nbcnd       : no. of boundary conditions
! nbvals      : 1st dimension size of "wdbcnd"
! nobcnd(:)   : domain no.
! kdbcnd(1,:) : for velocity
! kdbcnd(2,:) : for temperature
! kdbcnd(3,:) : for mass fraction
! kdbcnd(4,:) : for RANS model
! kdbcnd(5,:) : for MHD:A
! kdbcnd(6,:) : for MHD:FAI
! kdbcnd(7,:) : for MHD: TYPE
! openout(:)  : for fire open outlet BC (only for outlet)
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
!-< 2. boundary domains >-
!
      integer,save,allocatable :: ivbcnd(:)
      integer,save,allocatable :: lvbcnd(:)
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
!-----------------------
! --- surface reaction
!-----------------------
!nph
      integer,parameter             :: gasphase=1
      integer,parameter             :: surphase=2
      integer,parameter             :: blkphase=3
      logical,save,allocatable      :: IDX_SUFRAC(:,:)
      integer,save,allocatable      :: heat(:)
      integer,save                  :: nallcomp
      real(8),save,     allocatable :: phs_dns(:)
      real(8),save,     allocatable :: phs_thk(:)
      CHARACTER*80,save,allocatable :: phs_nam(:),phs_snm(:)
      integer,save,     allocatable :: phs_idx(:),phs_nbc(:),
     &                                 phs_typ(:),phs_com(:),
     &                                 phs_comr(:)
      real(8) ,save,allocatable     :: sigma_site(:),bi_site(:) 
      real(8) ,save,allocatable     :: m_indx(:),suf_dns(:)
      real(8) ,save,allocatable     :: n_bi(:,:),n_site(:,:)
      real(8) ,save,allocatable     :: num_site(:)
      real(8) ,save,allocatable     :: phs_inifrac(:),blk_dns(:)
!
      integer,save,allocatable      :: surfreac(:)
!
      integer,private :: iset1=0, iset2=0
      character*80,allocatable :: boundName(:,:)
      integer                  :: numbname,undef_bcno=0
      character*1,allocatable  :: prdcAxis(:)
      real(8),allocatable      :: prdcOffset(:,:)
!
      real(8),save,allocatable      :: rghnes(:)
      real(8),save,allocatable      :: masflx(:)
      real(8),save,allocatable      :: sumflx(:)
      real(8),save,allocatable      :: sumare(:)
      integer,save,allocatable      :: masflg(:)
      integer,save,allocatable      :: masbc(:)
      real(8),save,allocatable      :: dvel(:)
      logical,save                  :: imasflg=.false.
      integer,save,allocatable      :: machflg(:)
      real(8),save,allocatable      :: machno(:)
      logical,save                  :: ifixmach=.false.
!
      integer,save,allocatable      :: partkd(:)
!
! --- moving mesh
!
      integer                       :: NBOUND=0,NBFS=0
      INTEGER,allocatable           :: NFBOUN(:)
      CHARACTER*80,allocatable      :: SFBOUN(:)
      INTEGER,allocatable           :: IFFACE(:,:,:)
      INTEGER,allocatable           :: IBIDX(:)
!
      real(8),save :: FC_CELL_VOLT=1.d0
      integer,save :: FC_V_flag=0
!///////////////////////////////////////////////////////////////////////
      contains
!
!
!< #1. check consistency of values in argument & module >---------------
!=================================================
      subroutine chkncomp(ival,sbnm)
!=================================================
      implicit none
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      character(LEN=10),parameter :: subnam='(chkncomp)'
!
      if( ival.eq.ncomp ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      CALL FFRABORT(1,'module_boundary/chkncomp')
      end subroutine chkncomp
!
!
!< #2. check consistency of values in argument & module >---------------
!=================================================
      subroutine chknrans(ival,sbnm)
!=================================================
      implicit none
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      character(LEN=10),parameter :: subnam='(chknrans)'
      if( ival.eq.nrans ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      CALL FFRABORT(1,'module_boundary/chknrans')
      end subroutine chknrans
!
!
!< #3. interpolate values at present time making use of time table >----
!=================================================
      subroutine setnow(time)
!=================================================
      implicit none
      real(8),intent(in) :: time
      integer :: nwbcnd
      if( nbcnd.lt.1 ) return
      call modutl_chkset('>',1,iset1,modnam//'(setnow)')
      nwbcnd=iwbcnd(nbcnd)
      call modutl_setnow(time,nbcnd,nbvals,nwbcnd,
     & iwbcnd,wdtbcnd,wdbcnd)
      end subroutine setnow
!
!
!< #4. return rotational matrix of periodic boundary >------------------
!=================================================
      subroutine set_rotmtrx(nbcnx,kd,nb,rot)
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
        CALL FFRABORT(1,'module_boundary/set_rotmtrx')
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
      end subroutine set_rotmtrx
!
!< #5. store domain data >----------------------------------------------
!==================================================================
      subroutine strdomain(fnam,ifle,nbdmn,nvbdmn,ncbdmn,
     & nobdmn,ivbdmn,lvbdmn,icbdmn,lcbdmn,ierr)
!==================================================================
      implicit none
      character(*),intent(in)  :: fnam
      integer     ,intent(in)  :: ifle,nbdmn,nvbdmn,ncbdmn
      integer     ,intent(in)  :: nobdmn(nbdmn)
      integer     ,intent(in)  :: ivbdmn(0:nbdmn),lvbdmn(nvbdmn)
      integer     ,intent(in)  :: icbdmn(0:nbdmn),lcbdmn(ncbdmn)
      integer     ,intent(out) :: ierr
      integer :: nob(nbcnd)
      integer :: i,j,k,m,n,kve,kce,nb
      character(LEN=11),parameter :: subnam='(strdomain)'
!
! --- 
!
      call modutl_chkset('>',1,iset1,modnam//subnam)
      call modutl_chkset('=',1,iset2,modnam//subnam)
!
      ierr=0
      kve=0
      kce=0
!
      do 100 nb=1,nbcnd
      k=0
      do 101 m=1,nbdmn
      if(nobdmn(m).eq.nobcnd(nb)) k=m
  101 continue
      nob(nb)=k
      if( k.lt.1 ) goto 100
      if( kdbinr(1,nb).gt.0 .and.
     &    icbdmn(k-1).ge.icbdmn(k) ) goto 9001
      kve=kve+(ivbdmn(k)-ivbdmn(k-1))
      kce=kce+(icbdmn(k)-icbdmn(k-1))
  100 continue
!
      allocate( ivbcnd(0:nbcnd), lvbcnd(kve),
     &          icbinr(0:nbcnd), lcbinr(kce),
     &          stat=ierr)
!
      if( ierr.ne.0 ) then
        write(ifle,'(1x,a)') ' ### ERR : allocation failed'
        goto 9999
      endif
!
      ivbcnd=0
      lvbcnd=0
      icbinr=0
      lcbinr=0
      i=0
      j=0
!
      do 200 nb=1,nbcnd
      if(nob(nb).lt.1) goto 209
      k=nob(nb)
      do 201 m=ivbdmn(k-1)+1,ivbdmn(k)
      i=i+1
      lvbcnd(i)=lvbdmn(m)
  201 continue
      do 202 m=icbdmn(k-1)+1,icbdmn(k)
      j=j+1
      lcbinr(j)=lcbdmn(m)
  202 continue
  209 continue
      ivbcnd(nb)=i
      icbinr(nb)=j
  200 continue
!
      return
 9001 continue
      write(ifle,'(a)') 'ERR: data error'
      write(ifle,'(2a)') 
     &       'the boundary domain has no additional ',
     &       'information about side of the face'
      write(ifle,*) 'boundary domain no. =',nobcnd(n)
      write(ifle,*) 'name of boundary file =',fnam(:len_trim(fnam))
 9999 continue
      write(ifle,*) modnam//subnam
      ierr=1
      end subroutine strdomain
!
!< #6. input namelist >------------------------------------------------
!
!======================================================================
      subroutine inputdata
     &    (ifli,ifll,ifle,my_rank,cntlnam,mmcomp,ncmpx,ncomp_sufx,
     &     nphasex,nphasx,nrnsx,NPOTNx,nneq,spinam,
     &     lrans,kemdl,lowke,E2P,firewal,lrot,lsuf,isurface,lcomp,
     &     mmeq,ncompall,cofleft,net_vreq,blktyp,
     &     chem_bc,alpl,ical_sld,ical_suf,
     &     KSUF,NOMAT,imhd,lFC,lpotn,
     &     iprtcle,idens,ierror)
!
!======================================================================
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle,KSUF,mmeq,ncmpx,nphasex,
     &                        chem_bc(mmeq),ncomp_sufx,ncompall,mmcomp,
     &                        NOMAT,lFC,NPOTNx,iprtcle,idens
      character(*),intent(in) :: cntlnam,spinam(mmcomp)
      integer,intent(in)   :: my_rank,nrnsx,nneq,alpl(2),imhd
      integer,intent(inout):: blktyp(mmcomp)
      integer,intent(in)   :: isurface
      logical,intent(inout):: lrans,E2P,firewal,kemdl,lowke
      logical,intent(in)   :: lrot,lcomp
      logical,intent(in)   :: lsuf,lpotn
      integer,intent(inout):: ical_sld,ical_suf
      real*8, intent(in)   :: cofleft(mmeq*ncompall),
     &                        net_vreq(mmeq*ncompall)
      integer,intent(out)  :: ierror
      integer,intent(inout):: nphasx
!
      character(LEN=11),parameter :: subnam='(inputdata)'
!
! --- [namelist]
!
      integer,parameter  :: msite=11,mbulk=11 
      integer,parameter  :: lc=10,lenspc=20 
      integer,parameter  :: mbctim=10 
      integer,parameter  :: mcomp=200,mrans=50,mneq=300 
      real(8),parameter  :: SML=1.d-15,SML1=1.d-7 
      real(8)            :: rdum,dumy1,dumy2 
      character(lc) :: kind,MHD_kind,potential_kind(mpotn),dumc 
      character(lc) :: particle_kind
      character(lc) :: vel,temp,conc,rans,MHD_A,MHD_F
!
      character(lenspc):: site_species_name(msite,mcomp),
     &                    bulk_species_name(mbulk,mcomp)
      character(lc) :: site_name(msite),bulk_name(mbulk)
      real(8) :: site_dens(msite),bulk_thick(mbulk)
      real(8) :: 
     &              site_species_molefraction(msite,mcomp),
     &              bulk_species_molefraction(mbulk,mcomp),
     &              bulk_species_dens(mbulk,mcomp)
      integer :: suf_reac_no(mneq)=0,heat_release=0
!
      integer :: no,side,ntime,profile,open_air=0,turn_wall=0
      real    :: time(mbctim)
      real    :: p(mbctim),u(mbctim),v(mbctim),w(mbctim)
      real    :: t(mbctim)
      real*8  :: ys(mcomp*mbctim),ys2(mcomp*mbctim)
      real    :: u2(mbctim),v2(mbctim),w2(mbctim),htc2,mtc2,vsgm2(2)
      real    :: t2(mbctim)
      real*8  :: aks(mrans*mbctim)
      real    :: cycle,htc,mtc,vsgm(2),rot(4)
      real    :: rpm,end_x,end_y,end_z,begin_x,begin_y,begin_z
      real    :: xoffset,yoffset,zoffset
      integer :: start_rot_step,rotup_step
      real*8  :: VOF_ang
      character*80 :: name,name2
      character*1  :: axis,slash
      integer :: discont=0
      real(8) :: roughness,massflux,Mach_nu
      real(8) :: FAI_mhd(2),AAA_mhd(3,2)
      integer :: flux_balance,fix_mach=0
      integer :: n_blade1,n_blade2
!
      integer :: porous_BC=0
      real(8) :: C2_porous,perm_porous,thickness_porous,porosity
!
      integer :: surface_reaction=0
!
      real(8) :: rademi=0.d0,radabsorb=0.d0,radtran=0.d0
      integer :: radtype=1,radgpflag=1,thinwall=0
!
      real(8) :: Potential_value(mpotn),V_cell
!      real(8) :: shutter_speed,shutter_start
      logical :: inflow
!
      namelist /boundary/ no,side,kind,ntime,MHD_kind,
     &                    vel,temp,conc,rans,
     &                    MHD_A,MHD_F,AAA_mhd,FAI_mhd,
     &                    time,p,profile,u,v,w,t,ys,aks,
     &                    cycle,htc,mtc,vsgm,rot,
     &                    name,name2,axis,xoffset,yoffset,zoffset,
     &                    u2,v2,w2,t2,ys2,htc2,mtc2,vsgm2,
     &                    open_air,
     &                    turn_wall,
     &                    rpm,end_x,end_y,end_z,
     &                    begin_x,begin_y,begin_z,
     &                    start_rot_step,rotup_step,
     &                    surface_reaction,
     &                    site_name,
     &                    site_species_name,
     &                    site_species_molefraction,
     &                    site_dens,
     &                    bulk_name,
     &                    bulk_thick,
     &                    bulk_species_name,
     &                    bulk_species_molefraction,
     &                    bulk_species_dens,
     &                    heat_release,
     &                    suf_reac_no,
     &                    discont,roughness,
     &                    massflux,flux_balance,
     &                    Mach_nu,fix_mach,
     &                    n_blade1,n_blade2,
     &                    porous_BC,C2_porous,perm_porous,porosity,
     &                    thickness_porous,
     &                    Potential_kind,Potential_value,
     &                    particle_kind,
     &                    V_cell,
     &                    rademi,radabsorb,radtran,radtype,radgpflag,
     &                    VOF_ang,
!     &                    shutter_speed,shutter_start
     &                    thinwall
!
! --- [local entities] 
!
      integer,parameter :: iundef=-huge(1)
      real   ,parameter :: undef=-huge(1.)
      real*8 ,parameter :: dundef=-huge(1.d0)
      real*8    :: aph
      logical :: lcvdwal=.false. 
      integer :: icom,iq,jq,nst,nm,ncm,nbk,nsite,nbulk,ity,
     &           iflag
      character(lenspc):: chknam(mcomp),name_tmp
! 
      character(lc),parameter ::
     &   periodic=   'periodic  ',
     &   symmetric=  'symmetric ',
     &   inlet=      'inlet     ',
     &   outlet=     'outlet    ',
     &   wall=       'wall      ',
     &   interface=  'interface ',
     &   touchinlet= 'touchinlet',
     &   firewall=   'firewall  ',
     &   static=     'static    ',
     &   stagnation= 'stagnation',
     &   buffle=     'buffle    ',
     &   sliding=    'sliding   ',
     &   cvdsurface= 'cvdsurface',
     &   shutter   = 'shutter   ',
     &   porous    = 'porous    ',
     &   overset   ='overset    '
!
      character(lc),parameter ::
     &   noslip=     'no-slip   ',
     &   freeslip=   'free-slip ',
     &   rough=      'rough-wall',  !u+=1/k*ln[Y+]+B-(B-8.5+1/kln(Ks+))
     &   satawall=   'sata-wall',   !u+=1/k*ln[(Yp+Ks)/Ks]
     &   model=      'model     ',
     &   dirichlet=  'Dirichlet ',
     &   neumann=    'Neumann   ',
     &   EWF=        'EWF       ',
     &   loglaw=     'log-law   ',
     &   transfer=   'transfer  ',
     &   InterDiri=  'InterDiri '
!
      character(lc),parameter ::
     &   mhd_inter=  'mhd_inter ',
     &   nature   =  'nature    ',
     &   define   =  'define    ',
     &   mhd_period=  'mhd_period'
!     &   potn_inter= 'interface '
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
      character(10),parameter :: 
     &    phstyp(3)=(/'Gas-Phase ',
     &    'Site-Phase','Bulk-Phase'/)
      character(lc),parameter ::
     &    mhdknd(2)=(
     &   /dirichlet,neumann/)

      character(lc),parameter ::
     &    mhdtyp(4)=(
     &   /mhd_inter,nature,define,mhd_period/)
!
      character(lc),parameter ::
     &    potnlis(4)=(
     &   /Dirichlet,neumann,interface,InterDiri/)
!
      character(lc),parameter ::
     &    parknd(3)=(
     &   /'stay      ',
     &    'reflect   ',
     &    'dead      '/)
!
      integer :: i,j,k,n,kt,nb,nwbcnd
      integer :: kd,kvel,ktmp,kcnc,krns,kmhd_A,kmhd_F,kmhd
      integer :: kpotn,kpart
      integer :: yfix,idum
      integer :: ios,ierr1,ierr2
      integer :: pitch_one=0,mallcomp
!      real*8  :: dumy1,dumy2
!
      call modutl_chkset('=',1,iset1,modnam//subnam)
!
!---------------------------
!-< 1. Initial set >-
!---------------------------
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
      nphasx=1      !=1 for gas phase
      nallcomp=ncmpx
      icom=0
!
      rewind ifli
  100 continue
      site_name=''
      site_species_name=''
      bulk_name=''
      bulk_species_name=''
      ntime=1
      read(ifli,boundary,iostat=ios)
      if(ios.lt.0) goto 101
!
      nbcnd=nbcnd+1
!
      if(undef_bcno.ne.0) then
        write(ifle,'(2a)')
     &' ERR: name [undefined] must be defined as lasted BC in ',cntlnam
        write(ifle,'(2a)')
     &' ERR: Or only one [undefined] BC can be defined in ',cntlnam
        call FFRABORT(1,'&namelist/boundary')
      endif
      if(name.eq.'undefined') then
        undef_bcno=nbcnd
      endif
!
      call nml_chksiz(ifle,'mbctim',mbctim,ntime,
     & modnam//subnam,ierr1)
      call nml_errmsg0(ifle,ios,'boundary',ierr2)
      if( ierr1.ne.0 .or. ierr2.ne.0 ) then
        write(ifle,'(a,2I4)') 
     &         'ERR: ierr1, ierr2=',ierr1 , ierr2
        write(ifle,'(a,I4)') 
     &         'ERR: sequence no. of the namelist =',nbcnd
        goto 9999
      endif
      nwbcnd=nwbcnd+max(1,ntime)
!
! --- 
!
      nsite=0
      nm=0
      do nst=1,msite
      IF(site_name(nst)/=' ') then
        nsite=nsite+1
        do ncm=1,mcomp
        if(site_species_name(nst,ncm)/=' ') then
          nm=nm+1
          do i=1,lenspc
          slash=site_species_name(nst,ncm)(i:i)
          if(slash=='/') then
            name_tmp=site_species_name(nst,ncm)(1:i-1)
            exit
          endif
          enddo
          iflag=1
          do i=1,icom
          if(trim(name_tmp)==trim(chknam(i))) then
            iflag=0
          endif
          enddo

          if(iflag==1) then
            icom=icom+1
            chknam(icom)=trim(name_tmp)
          endif

        endif
        enddo
      endif
      enddo
      if(nsite>=msite) 
     & call FFRABORT
     &  (1,'ERR: SITE number in one substrate is limited in 10')
      nphasx=nphasx+nsite
      nallcomp=nallcomp+nm
      nbulk=0
      nm=0
      do nbk=1,mbulk
      if(bulk_name(nbk)/=' ') then
        nbulk=nbulk+1
        do ncm=1,mcomp
        if(bulk_species_name(nbk,ncm)/=' ') then
          nm=nm+1

          do i=1,lenspc
          slash=bulk_species_name(nbk,ncm)(i:i)
          if(slash=='/') then
            call FFRABORT(1,'ERR: Bulk species name have [/]')
          endif
          enddo
          name_tmp=trim(bulk_species_name(nbk,ncm))
          iflag=1
          do i=1,icom
          if(trim(name_tmp)==trim(chknam(i))) then
            iflag=0
          endif
          enddo
          
          if(iflag==1) then
            icom=icom+1
            chknam(icom)=trim(name_tmp)
          endif
        endif
        enddo
      endif
      enddo
      if(nbulk>=mbulk) call FFRABORT(1,'ERR: BULK number too many')
      nphasx=nphasx+nbulk
      nallcomp=nallcomp+nm
      goto 100
  101 continue
      mallcomp=nallcomp
!
! --- 
!
      if(nphasx/=nphasex) then
        write(ifle,*) 'nphasx/=nphasex =>',nphasx,nphasex
        call FFRABORT(1,'ERR: nphasx/=nphase')
      endif
!
      if(nsite==0.and.nbulk>0) then
        call FFRABORT(1,'ERR: Surface SITE NOT be defined')
      endif
!
      if(ncmpx+icom/=ncmpx+ncomp_sufx) then
        write(ifle,'(1x,a,/,1x,a,4I4)') 
     &   'ERR: Maybe SITE & BULK species name seting error', 
     &   ' Or: Having duplicated species name between phasees ',
     &   nallcomp,ncmpx+ncomp_sufx,ncmpx,ncomp_sufx
        call FFRABORT(1,'ERR: &boundary, ncomp_suf error')
      endif
!nph
      allocate(phs_idx(0:nphasex),
     &         phs_dns(nphasex),
     &         phs_nam(nphasex),
     &         phs_nbc(nphasex),
     &         phs_typ(nphasex),
     &         phs_thk(nphasex),
     &         sigma_site(nphasex),bi_site(nphasex),
     &         phs_snm(mallcomp),
     &         phs_com(mallcomp),
     &         phs_comr(mallcomp),
     &         num_site(mallcomp),
     &         phs_inifrac(mallcomp),
     &         blk_dns(mallcomp),
     &         stat=ierr1)
      if(ierr1.ne.0) 
     &   call FFRABORT(1,'ERR: &boundary allcating error')
      phs_idx=0
      phs_dns=0.d0
      phs_nam=''
      phs_nbc=0
      phs_typ=0   !gasphase  !=1
      phs_snm=''
      phs_inifrac=0.d0
      phs_thk=0.d0
      blk_dns=0.d0
      num_site=1.d0
!-----------------------
! --- gas phase
!-----------------------
      phs_idx(1)=ncmpx
      phs_dns(1)=0.d0
      phs_nam(1)='Gas'
      do i=1,ncmpx
      phs_snm(i)=spinam(i)
      enddo
      phs_nbc(1)=0
      phs_typ(1)=gasphase
!---------------------------
!--< 1.2 allocate arrays >--
!---------------------------
      iscomp=5
      israns=5+ncomp
      nbvals=max(israns+nrans,9)
      allocate( nobcnd(0:nbcnd),
!     &          bcndno(0:nbcnd),
     &          LBC_INDEX(0:nbcnd),
     &          LBC_pair(0:nbcnd), 
     &          MAT_BCIDX(0:nbcnd,2), 
     &          kdbcnd(0:7,nbcnd), 
     &          openout(nbcnd), 
     &          poro(nbcnd),C2pBC(nbcnd),alpha_coBC(nbcnd),
     &                           thkposBC(nbcnd),prosty(nbcnd),
     &          surfreac(nbcnd), 
     &          rotwall(nbcnd), 
     &          endw(3,nbcnd),beginw(3,nbcnd), 
     &          omega(nbcnd), 
     &          strtimw(nbcnd), 
     &          rotupw(nbcnd), 
     &          iwbcnd(0:nbcnd), 
     &          adbcnd(4,nbcnd),
     &          wdbcnd(0:nbvals,nbcnd),
     &          MHD_bnd(4,nbcnd,2),
     &          wdtbcnd(0:nbvals,nwbcnd),
     &          boundName(nbcnd,1:2),
     &          prdcAxis(nbcnd),prdcOffset(nbcnd,3),
     &          kdbinr(2,nbcnd),
     &          adbcnd2(4,nbcnd),
     &          wdbcnd2(0:nbvals,nbcnd),
     &          wdtbcnd2(0:nbvals,nwbcnd),
     &          stagvel(nbcnd),
     &          stagare(nbcnd),
     &          idis(nbcnd),
     &          nblade(nbcnd,2),
     &          rotsld(nbcnd),
     &          IDX_SUFRAC(nbcnd,max(nneq,1)),
     &          heat(nbcnd),
     &          n_bi(nbcnd,max(nneq,1)),
     &          n_site(nbcnd,max(nneq,1)),
     &          m_indx(max(nneq,1)),
     &          suf_dns(nbcnd),
     &          sizein(nbcnd),
     &          rghnes(nbcnd),
     &          masflx(nbcnd),
     &          masflg(nbcnd),
     &          sumflx(nbcnd),
     &          sumare(nbcnd),
     &          machflg(nbcnd),
     &          machno(nbcnd),
     &          masbc(nbcnd),
     &          partkd(nbcnd),
     &          dvel(nbcnd),
     &          kdpotn(mpotn,nbcnd),vapotn(mpotn,nbcnd),mem_cat(nbcnd),
     &          radprop(2,nbcnd),radwalltype(nbcnd),
     &          vofang(nbcnd),
     &          shut_strt(nbcnd),
     &          pkdwall(nbcnd),
     &          stat=ierr1 )
      if(ierr1.ne.0) then
        write(ifle,'(a)') 'ERR: allocation failed'
        goto 9999
      endif
!
      nobcnd =0
!      bcndno =0
      kdbcnd =kxnone
      openout=0
      poro(:)=0
      C2pBC(:)=0.d0
      thkposBC(:)=0.d0
      prosty(:)=0.d0
      alpha_coBC(:)=0.d0
      rotwall=0
      endw=0.d0
      beginw=0.d0
      omega=0.d0
      strtimw=0.d0
      rotupw=0.d0
      iwbcnd =0
      adbcnd =0.d0
      adbcnd2=0.d0
      wdbcnd =0.d0
      MHD_bnd =0.d0
      wdtbcnd=0.d0
      wdbcnd2 =0.d0
      wdtbcnd2=0.d0
      kdbinr =0
      numbname=2
      prdcAxis=' '
      prdcOffset=0.d0
      IDX_SUFRAC(:,:)=.false.
      heat(:)=0
      idis(:)=0
      nblade(:,:)=0
      rotsld(:)=0
      rghnes(:)=0.d0
      masflx(:)=0.d0
      masflg(:)=0
      machflg(:)=0
      machno(:)=0.d0
      kdpotn(:,:)=kpneum
      vapotn(:,:)=0.d0
      mem_cat(:)=0
      partkd(:)=kpart_s
!
      radprop(:,:)=0.d0
      radwalltype(:)=0
!
      vofang(:)=0.d0
!
      
      shut_strt(:)=0.d0
      pkdwall(:)=.false.
      surfreac(:)=0
!-------------------------
      
!-< 2. Input namelist >---
!-------------------------
!--< 2.1 read namelist >--
!-------------------------
      kind = ' '
      MHD_kind= ' '
!
      nb=0
      kt=0
      nphasx=1      !=1 for gas phase
      nallcomp=ncmpx
      rewind ifli
! --- start read namelist file for every BC
 1000 continue
      vel   = ' '
      temp  = ' '
      conc  = ' '
      rans  = ' '
      MHD_A = ' '
      MHD_F = ' '
      no    = iundef
      side  = iundef
      ntime = 1
      time  = 0.d0
      cycle = undef
      p     = undef
      profile=iundef
      u     = undef
      v     = undef
      w     = undef
      t     = undef
      ys    = dundef
      aks   = dundef
      htc   = undef
      mtc   = undef
!
      vsgm  = undef
      rot   = 0.
      name  = ' '
      name2 = ' '
      axis  = ' '
      u2    = undef
      v2    = undef
      w2    = undef
      t2    = undef
      ys2   = dundef
      htc2  = undef
      mtc2  = undef
      vsgm2 = undef
      xoffset = 0.0
      yoffset = 0.0
      zoffset = 0.0
      open_air=0
! --- sliding frame
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
! --- surface reaction
      site_species_name=''
      bulk_species_name=''
      site_name=''
      bulk_name=''
      site_dens=undef
      bulk_species_dens=undef
      suf_reac_no=iundef
      bulk_species_molefraction=undef
      site_species_molefraction=undef
      heat_release=0
      discont=0
      roughness=undef
      massflux=undef
      flux_balance=0
      Mach_nu=undef
      fix_mach=0
      n_blade1=1
      n_blade2=1
!
! --- MHD
!
      AAA_mhd=0.d0
      FAI_mhd=0.d0
!
! --- 
!
      porous_BC=0
      C2_porous=undef
      perm_porous=undef
      thickness_porous=undef
      porosity=undef
!-----------------------------
! --- 
!----------------------------- 
      surface_reaction=0
!
      Potential_kind(:)=' '
      Potential_value(:)=undef
!
      particle_kind=' '
!
      V_cell=undef
!
      VOF_ang=undef
!
!      shutter_speed=undef
!      shutter_start=undef
!
      read(ifli,boundary,iostat=ios)
      if(ios.lt.0) goto 1001
!      write(ifll,boundary)
      call nml_errmsg0(ifle,ios,'boundary',ierr1)
      if(ierr1.ne.0) goto 9999
      nb=nb+1
!----------------------------
!--< 2.2 common procedure >--
!----------------------------
!/ global kind /
!----------------------------
      call nml_listno(16,kndlst,kind,kd)
      call nml_chkchr0(ifle,'kind',kind,kd,ierr1)
      if( ierr1.ne.0 ) goto 9999
      kind=kndlst(kd)
!
      if(imhd>0) then
        call nml_listno(4,mhdtyp,MHD_kind,kmhd)
        call nml_chkchr0(ifle,'MHD_kind',MHD_kind,kmhd,ierr1)
        if( ierr1.ne.0 ) goto 9999
        MHD_kind=mhdtyp(kmhd)
        if(MHD_kind==mhd_inter) then
          kdbcnd(7,nb)=kdinterf
        elseif(MHD_kind==nature) then
          kdbcnd(7,nb)=kdnature
        elseif(MHD_kind==define) then
          kdbcnd(7,nb)=kddefine
        elseif(MHD_kind==mhd_period) then
          kdbcnd(7,nb)=kdperiod
        endif
        
        if(kind/=interface.and.MHD_kind==mhd_inter) then
          write(ifle,'(1x,a,I4)') 'ERR: boundary no= ',no
          call FFRABORT
     & (1,'ERR: IF kind/="interface" => MHD_kind="nature" or "define"')
        endif
        if(kind==interface.and.MHD_kind/=mhd_inter) then
          write(ifle,'(1x,a,I4)') 'ERR: boundary no= ',no
          call FFRABORT
     &    (1,'MSG: IF kind=="interface" => MHD_kind="mhd_inter"')
        endif
        IF(MHD_kind==define) then
          call nml_listno(2,mhdknd,MHD_A,kmhd_A)
          call nml_listno(2,mhdknd,MHD_F,kmhd_F)
          if( kmhd_A.gt.0 ) MHD_A=mhdknd(kmhd_A)
          if( kmhd_F.gt.0 ) MHD_F=mhdknd(kmhd_F)
          if(MHD_A.eq.dirichlet) kdbcnd(5,nb)=kmdirc
          if(MHD_A.eq.neumann)   kdbcnd(5,nb)=kmneum
          if(MHD_F.eq.dirichlet) kdbcnd(6,nb)=kmdirc
          if(MHD_F.eq.neumann)   kdbcnd(6,nb)=kmneum
          MHD_bnd(1,nb,1)=AAA_mhd(1,1)
          MHD_bnd(2,nb,1)=AAA_mhd(2,1)
          MHD_bnd(3,nb,1)=AAA_mhd(3,1)
          MHD_bnd(4,nb,1)=FAI_mhd(1)
          MHD_bnd(1,nb,2)=AAA_mhd(1,2)
          MHD_bnd(2,nb,2)=AAA_mhd(2,2)
          MHD_bnd(3,nb,2)=AAA_mhd(3,2)
          MHD_bnd(4,nb,2)=FAI_mhd(2)
        endif
      endif
!
      if(lFC>0.or.lpotn) then
        do I=1,npotnx
        kpotn=0
        dumc=Potential_kind(I)
        call nml_listno(4,potnlis,dumc,kpotn) 
        call nml_chkchr0
     &    (ifle,'Potential_kind',dumc,kpotn,ierr1)
        if( ierr1.ne.0 ) then
          write(ifle,'(1X,a,I4)') 'ERR: [Potential_kind] error at BC=',
     &   no
          call FFRABORT(1,'')
        endif
        if(kpotn==kpdirc.or.kpotn==kpneum) then
          if(Potential_value(I)==undef) then
            write(ifle,'(1X,a,2I4)') 
     &      'ERR: Define [Potential_value(I)] at BC=, I=',no,I
            call FFRABORT(1,'ERR: set [Potential_value(:)=?] ')
          endif
        endif
        kdpotn(I,nb)=kpotn
        vapotn(I,nb)=Potential_value(I)

        if(kpotn==kpintr) then
          if(kind/=buffle.and.kind/=interface) then
            write(ifle,'(1X,a,2I4)') 
     &      'ERR: Error Define [Potential_kind(I)] at BC=, I=',no,I
            write(ifle,'(1X,3a)') 
     &      'MSG: Potential_kind=interface only for ',
     &      'kind=buffle or ',
     &      'kind=interface'
            call FFRABORT(1,'MSG: reset [Potential_kind] ')
          endif
        endif
        if(kpotn==kpindr) then
          if(kind/=buffle.and.kind/=interface) then
            write(ifle,'(1X,a,2I4)') 
     &      'ERR: Error Define [Potential_kind(I)] at BC=, I=',no,I
            write(ifle,'(1X,3a)') 
     &      'MSG: Potential_kind="InterDiri" only for ',
     &      'kind=buffle or ',
     &      'kind=interface'
            call FFRABORT(1,'MSG: reset [Potential_kind] ')
          endif
          if(Potential_value(I)==undef) then
            write(ifle,'(1X,a,2I4)') 
     &      'ERR: To Define [Potential_value(I)] at BC=, I=',no,I
            write(ifle,'(1X,2a)') 'MSG: Potential I will refer to',
     &      ' number of potential defined by Potential_value(I)'
            call FFRABORT(1,'ERR-1: set [Potential_value(:)=?] ')
          elseif(int(Potential_value(I))==0) then
            write(ifle,'(1X,a,2I4)') 
     &      'ERR: To Define [Potential_value(I)] at BC=, I=',no,I
            write(ifle,'(1X,2a)') 'MSG: Potential I will refer to',
     &      ' number of potential defined by Potential_value(I)'
            call FFRABORT(1,'ERR-2: set [Potential_value(:)=?] ')
          endif
        endif
        enddo


        if(V_cell/=undef) then
          if(FC_V_flag>0) then
            call FFRABORT(1,'ERR: V_cell defined only once for PEFC')
          endif
          FC_V_flag=1
          FC_CELL_VOLT=V_cell
!          if(lFC>0) then
!            if(abs(V_cell-Potential_value(1))>1.d-10)then
!              call FFRABORT
!     &   (1,'ERR: V_cell defined NOT equal to potential value for PEFC')
!            endif
!          endif
        endif
      endif
!
      if(kind.eq.periodic.or.kind.eq.symmetric) ntime=1
      kt=kt+ntime
      iwbcnd(nb)=kt
!----------------------------
!/ domain no. & side no. /
!----------------------------
      if( kind.eq.periodic )  call chkintn('side',side)
      if( kind.eq.interface ) call chkintn('side',side)
      if( ierror.ne.0 ) goto 9999
!
      if(no.eq.iundef) then
        write(ifle,'(a)') 'ERR: data error'
        write(ifle,'(a)') 'lack of data : no'
        goto 9999
      endif
      do 200 n=1,nb-1
      if(no.ne.nobcnd(n) ) goto 200
      ierr1=0
      if( side.eq.iundef ) then
        ierr1=n
      elseif( kdbinr(1,n).gt.0 ) then
        ierr1=kdbinr(1,n)
      else
        kdbinr(1,n )=nb
        kdbinr(1,nb)=n
      endif
      if( ierr1.ne.0 ) then
      write(ifle,'(a)') 'ERR: data error'
      write(ifle,'(a,I4)') 'MSG: no =',no
      write(ifle,'(a)') 
     &   'MSG: this no. has been specified in previous line'
      write(ifle,'(a,I4)') 
     &   'MSG: sequence no. of the line =',ierr1
      goto 9999
      endif
  200 continue
      nobcnd(nb)=no
!      bcndno(no)=nb
      boundName(nb,1)=name
      if(name.eq.' ') then
        write(ifle,'(a,I4,a)') 
     &   'EER: Boundary no. ',no,' must be specified'
      write(ifle,'(a)') 
     &   'MSG: Please input the name corresponding to grid file'
        call FFRABORT(1,'module_boundary/inputdata')
      endif
      boundName(nb,2)=name2
      prdcAxis(nb)=axis
      prdcOffset(nb,1)=dble(xoffset)
      prdcOffset(nb,2)=dble(yoffset)
      prdcOffset(nb,3)=dble(zoffset)
!
      if(side.ne.iundef) then
        if(side.lt.1.or.side.gt.2) then
          write(ifle,'(a)') 'ERR: data error'
          write(ifle,'(a,I4)') 'side =',side
          write(ifle,'(a)') 'it must be 1 or 2'
          goto 9999
        endif
        if( kdbinr(1,nb).gt.0 ) then
          if( kdbinr(2,kdbinr(1,nb))+side.ne.3 ) then
            write(ifle,'(a)') '### error : data error'
            write(ifle,'(a,I4)') 'side =',side
            if( side.eq.1 ) write(ifle,'(a)') 'it must be 2'
            if( side.eq.2 ) write(ifle,'(a)') 'it must be 1'
            goto 9999
          endif
        endif
        kdbinr(2,nb)=side
      endif
!-------------------------
! --- 
!-------------------------
      if(kind/=outlet) then
        if(open_air/=0) then
        write(ifle,'(2a)') 
     &  'ERR: open_air MUST be 0, While NOT [outlet] BC in ',
     &   cntlnam
          call FFRABORT(1,'ERR: re-set open_air=0')
        endif
      endif
!-------------------------
! --- 
!-------------------------
      if(kind/=wall.and.
     &   kind/=cvdsurface.and.
     &   kind/=firewall.and.
     &   kind/=interface.and.
     &   kind/=outlet
     &    ) then
        if(turn_wall/=0) then
          write(ifle,'(2a)') 
     &    'ERR: turn_wall MUST be 0, While NOT [wall] BC in ',
     &    cntlnam
          call FFRABORT(1,'ERR: re-set turn_wall=0 ')
        endif
      endif
!--------------------------------------------------------
! --- rot wall for all wall-kind BC and outlet BC
!--------------------------------------------------------
      if(turn_wall==1) then
        rotwall(nb)=1
        omega(nb)=rpm/60.d0*2.d0*3.1415926d0
        endw(1,nb)=end_x
        endw(2,nb)=end_y
        endw(3,nb)=end_z
        beginw(1,nb)=begin_x
        beginw(2,nb)=begin_y
        beginw(3,nb)=begin_z
        strtimw(nb)=dble(max(start_rot_step,0))
        rotupw(nb)=dble(max(rotup_step,1))
        rdum=(end_x-begin_x)*(end_x-begin_x)
     &      +(end_y-begin_y)*(end_y-begin_y)
     &      +(end_z-begin_z)*(end_z-begin_z)
        if(rdum.lt.SML) then
          write(ifle,*) 'ERR: You must define revolution shaft '
          write(ifll,*) 
     &            'MSG: Vector of revolution shaft is defined as ',
     &            'with : end_x,end_y,end_z',
     &            'and  : begin_x,begin_y,begin_z'
        call FFRABORT(1,'ERR: NOT define revol. shaft at turning-wall')
        endif
      endif
!
!
!
      if(porous_BC==1) then
        if(C2_porous==undef.or.perm_porous==undef.or.
     &     thickness_porous==undef.or.porosity==undef) then
          write(ifle,'(1X,A,1X,A)') 
     &  "ERR: NOT define C2_porous, perm_porous porosity or",
     &  "  thickness_porous ,for Darcys law "
          call FFRABORT
     &   (1,'ERR: NOT define C2_porous, perm_porous')
        else
          C2pBC(nb)=C2_porous
          alpha_coBC(nb)=perm_porous
          thkposBC(nb)=thickness_porous
          prosty(nb)=porosity
        endif
      endif
!
! --- surface reaction
!
      if(   kind/=firewall
     & .and.kind/=buffle
     & .and.kind/=cvdsurface
     & .and.kind/=interface
     & .and.kind/=wall
     & ) then
        if(surface_reaction==1) then
          write(ifle,'(1x,a,I4)') 
     &    'ERR: Surface reaction setting error at BC=',no
          write(ifle,'(1x,2a)') 
     &     'MSG: Surface reaction only for ',
     &     'kind="firewall","interface","buffle","wall","cvdsurface"'
           call FFRABORT(1,'ERR: ')
        endif
      endif
!
!--------------------------------------------------------------------
!radwalltype:		1,diffuse-wall	,OK 
!radtype                2,mirror-wall	,OK
!			3,directional-diffuse wall,OK
!			4,user defined
!--------------------------------------------------------------------
!
      radprop(1,nb)=dble(rademi)
      radprop(2,nb)=dble(radabsorb)
      radwalltype(nb)=radtype
!
!
!
!
      inflow=kind.eq.inlet.or.kind.eq.outlet.or.kind.eq.touchinlet
      if(inflow) then
        masflg(nb)=flux_balance
        imasflg=imasflg.or.flux_balance>=1
        if(flux_balance>=1) then
          if(my_rank==0) then
            write(ifll,'(a)') 
     &      'MSG: Set [massflux] for [inlet] BC'
            write(ifll,'(3a)') 
     &      'MSG: (1)massflux>0 : inflow; ',
     &      'MSG: (2)massflux<0 : outflow ',
     &      'MSG: (3)massflux=0 : by Up-flow '
            write(ifll,'(1X,a)') 'MSG: flux_balance=1: massflux= KG/S'
            write(ifll,'(1X,a)') 'MSG: flux_balance=2: massflux= M^3/S'
            if(massflux==undef) then
              write(ifll,'(a)') 
     &        'ERR: Set [massflux] for [inlet] BC'
              write(ifll,'(a)') 
     &        'MSG: Set [massflux] for [inlet] BC'
              write(ifll,'(2a)') 
     &        'MSG: (1)massflux>0 : inflow; ',
     &        'MSG: (2)massflux<0 : outflow ',
     &        'MSG: (3)massflux=0 : by Up-flow '
              call FFRABORT(1,'ERR: NOT defined [massflux]')
            else
            endif
          endif
          masflx(nb)=-massflux
        endif
      endif
!
!
!
      if(fix_mach==1) then
        if(.NOT.
     &    ((kind==outlet.and.open_air==10).or.(kind==inlet)))
     &    then
          write(ifle,'(a)') 
     &     'ERR: fix_mach==1 support [open_air=10] BC'
          write(ifle,'(a)') 
     &     'ERR: and support [inlet] BC'
          call FFRABORT(1,'ERR: in &boundary')
        endif
      endif
!
      if(masflg(nb)>=1.and.fix_mach==1) then
        write(ifle,'(a)') 
     &    'ERR: fix_mach=1 NOT support flux_balance=1 or 2'
        call FFRABORT(1,'ERR: in &boundary')
      endif
!-------------------
! --- particle BC
!-------------------
      if(iprtcle==1.or.iprtcle==2) then
        call nml_listno(3,parknd,particle_kind,kpart)
        call nml_chkchr0
     &    (ifle,'particle_kind',particle_kind,kpart,ierr1)
        if(ierr1.ne.0) then
          write(ifle,'(1x,a,I4)') 'ERR: boundary no= ',no
          write(ifle,'(1x,4a)') 'ERR-MSG: set particle_kind=',parknd
          call FFRABORT(1,'ERR: set [particle_kind] in &boundary')
        endif
        partkd(nb)=kpart
        if(iprtcle==2) then
          if(kpart/=kpart_s) then
            write(ifle,'(1X,a,I4)') 
     &      'MSG: particle_kind set error at boundary no=',no
            write(ifle,'(1X,2a)') 
     &      'ERR-MSG: particle_kind="reflect"] ',
     &       'MUST be set, while [rain_WDR=1]'
     &      
            call FFRABORT(1,'ERR: &boundary/particle_kind')
          endif
        endif
      endif
      

! ------------------------------------------------
! --- 
! ------------------------------------------------
!
      if( kind.eq.symmetric ) goto 1101
      if( kind.eq.periodic )  goto 1102
      if( kind.eq.inlet )     goto 1103
      if( kind.eq.outlet )    goto 1103
      if( kind.eq.touchinlet) goto 1104  
      if( kind.eq.firewall)   goto 1105  !surface
      if( kind.eq.buffle)     goto 1105  !surface
      if( kind.eq.static)     goto 1106
      if( kind.eq.stagnation) goto 1106
      if( kind.eq.sliding)    goto 1107
      if( kind.eq.cvdsurface) goto 1105  !surface
      if( kind.eq.porous    ) goto 1105
      if( kind.eq.overset   ) goto 1108
!
      goto 1105  !wall  !surface
!-------------------------------
!--< 2.3 symmetric boundary >--
!-------------------------------
 1101 continue
      kdbcnd(0,nb)=kdsymm
      goto 1000
!-----------------------------
!--< 2.4 periodic boundary >--
!-----------------------------
 1102 continue
!      if(imhd>0) call FFRABORT(1,'ERR: MHD CANNOT use "periodic BC"')
      if(name2.eq.' ') then
        write(ifle,'(a,I4,2a)') 
     &   'EER: Boundary no. ',no,' is periodic boundary',
     &  ' its counter-BCs name must be defined'
      write(ifle,'(a)') 
     &    'MSG: Please input the name corresponding to grid file'
        call FFRABORT(1,'module_boundary/inputdata')
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
        call FFRABORT(1,'module_boundary/inputdata')
      endif
      goto 1000
!---------------------------------
!--< 2.5 inlet/outlet boundary >--
!---------------------------------
 1103 continue
!
      if((fix_mach==1).and.lcomp) then
          machflg(nb)=fix_mach
          machno(nb)=Mach_nu
          ifixmach=ifixmach.or.fix_mach==1
      elseif((fix_mach==1).and..NOT.lcomp) then
          write(ifle,'(a)') 
     &    'ERR: fix_mach==1 NOT support un-compressible'
          call FFRABORT(1,'ERR: in &boundary')
      else
          machflg(nb)=fix_mach
          machno(nb)=Mach_nu
      endif
!
      call chkntime
      call chkflt1(0,ntime,'u',u)
      call chkflt1(0,ntime,'v',v)
      call chkflt1(0,ntime,'w',w)
      if(open_air==5) then
        call chkflt1(0,ntime,'t',t)
      else
        call chkflt1(1,ntime,'t',t)
      endif
      call chkflt2(1,ncomp,ntime,'ys',ys)
!
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
!---------------
! --- INLET BC:
!---------------
      if(kind.eq.inlet) then
        kdbcnd(0,nb)=kdilet
!
        if(profile.eq.iundef) then
          write(ifle,'(a)') 'ERR: Inlet BC error : data error'
          write(ifle,'(a)') 'Define [profile] = 0,1 ro 2  in Inlet BC'
          write(ifle,'(a)') '[profile=0]: Uniform Vel. Profile'
          write(ifle,'(a)') '[profile=1]: User defined Vel. Profile'
          write(ifle,'(a)') '[profile=2]: Developed Vel. profile'
          CALL FFRABORT(1,'module_boundary/inputdata')
        elseif(profile.eq.0) then
          distrb=0
        elseif(profile.eq.1) then
          distrb=1
        elseif(profile.eq.2) then
          distrb=2
        else
          write(ifle,'(a)') 'MSG: Must be profile=0,1 or 2'
          CALL FFRABORT(1,'module_boundary/inputdata')
        endif
!
!        if((flux_balance==1.or.idens==4).and.lcomp) then
        if((flux_balance==1).and.lcomp) then
          call chkflt1(1,ntime,'p',p)
          if(ierror==1) then
!            write(ifle,'(a)') 
!     &      'MSG: Real Gas Model NEEDS defining Inflow pressure'
            write(ifle,'(2a)') 
     &      'MSG: While [flux_balance=1], ',
     &      'Compressible or zero model NEED defining Inflow pressure'
            CALL FFRABORT(1,'module_boundary/inputdata')
          endif
        else
          call chkfltn('p(1)',p(1))
        endif
!
!
      else
!---------------
! --- OUTLET BC
!---------------
        if(porous_BC==1) then
          poro(nb)=1
          ical_poroBC=1
        endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       openout(nb)=0      ! incompressible: dp/dx=0
                           ! compressible:   p=given (ONLY subsonic flow)
                           ! zero-mach :     dp/dx=0
        if(open_air==1) then 
          openout(nb)=1    ! (1) zero-mach open air 
                           ! (2) compressible open air 
                           ! (3) incompressible open air ???
        elseif(open_air==2) then 
          openout(nb)=2    ! suck-out BC for compresible 
        elseif(open_air==3) then
          openout(nb)=3    ! incompresible zero-pressure-out BC
        elseif(open_air==0) then
          openout(nb)=0    ! default
        elseif(open_air==4) then
          openout(nb)=4    ! sliding fan for incompressible 
        elseif(open_air==5) then
          openout(nb)=5    ! : Euler 2 Phase liquid surface 
                           ! : wall for liquid surface
                           ! : outlet for gas
        elseif(open_air==6) then 
          openout(nb)=6    ! : compresible inflow, 
                           !   p=F(p0,V,T) , DV/Dx=0.d0 
                           !   (p0=given)
        elseif(open_air==7) then 
          openout(nb)=7    ! :   
                           !    compresible outflow : ical_dens==4
                           !    dp/dx=0 or p=given
                           !    dv/dx=0
        elseif(open_air==8) then 
          openout(nb)=8    ! :  compresible inflow
                           !    (stagnation-pressure: p0=given)
                           !    Flux_balance=0 : NOT finshed
                           !    Flux_balance=1 : V=mass-rate; Dp/Dx=0 
                           !    V=F(p0,T); rho=F(p,T)
        elseif(open_air==9) then 
          openout(nb)=9    ! :  compresible outflow : 
                           !    all Mach number flow
                           !    (subsonic, transonic ,supersonic flow)
                           !                          DV/Dx=0.d0 
                           !               Mach>1.0 : Dp/Dx=0.d0
                           !               Mach<1.0 : p=given
        elseif(open_air==10) then 
          openout(nb)=10   ! :  compresible inflow : 
                           !                  (V=fix ,p=fix)
                           !                  fix_mach=0:  V= given 
                           !                  fix_mach=1:  Mach=given
                           !                         p=F(p0,V,T)
                           !                         (p0=given)
        else
          call FFRABORT(1,'ERR: open_air error')
        endif              ! No in-flow
!
! --- 
!       
!       'stagnation'       ! inflow BC:
                           ! compresible inflow : 
                           !               V=F(p0,T),(p0=given)
                           !               dp/dx=0
                           ! incompresible inflow : 
                           !               p=F(p0,T),(p0=given)
                           !               dv/dx=0
!                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        kdbcnd(0,nb)=kdolet
        call chkflt1(1,ntime,'p',p)
        kdbcnd(2,nb)=ktneum
        if( temp.eq.dirichlet ) kdbcnd(2,nb)=ktdirc
        if( temp.eq.neumann )   kdbcnd(2,nb)=ktneum
        if( temp.eq.transfer )  kdbcnd(2,nb)=kttrns
        
        kdbcnd(3,nb)=kyneum
!        if( conc.eq.dirichlet ) kdbcnd(3,nb)=kydirc
!        if( conc.eq.neumann )   kdbcnd(3,nb)=kyneum
!        if( conc.eq.transfer )  kdbcnd(3,nb)=kytrns
!        if( conc.eq.EWF )  kdbcnd(3,nb)=kyEWF
      endif
!
      
      if(lrans) then
        call chkflt2(1,nrans,ntime,'aks',aks)
      endif
!
      if(kind.eq.inlet.and.vsgm(1).ne.undef) then
        call chkflt3('vsgm(1)',vsgm(1))
        call chkflt3('vsgm(2)',vsgm(2))
      endif
!
      if(E2P) then
        if(kind.eq.inlet) then
          aph=aks(alpl(1))+aks(alpl(2))
          if(aks(alpl(1))>1.d0.or.aks(alpl(2))>1.d0) then
            write(ifle,'(1x,a,I4)') 'EER: Boundary no. ',no
            call 
     &   FFRABORT(1,'ERR-1:Volt Fraction is wrong setting in INLET')
          elseif(aks(alpl(1))<0.d0.or.aks(alpl(2))<0.d0) then
            write(ifle,'(1x,a,I4)') 'EER: Boundary no. ',no
            call  
     &   FFRABORT(1,'ERR-2:Volt Fraction is wrong setting in INLET')
          elseif(abs(aph-1.d0).gt.SML1) then
            write(ifle,'(1x,a,I4)') 'EER: Boundary no. ',no            
            call  
     &   FFRABORT(1,'ERR-3:Volt Fraction is wrong setting in INLET')
          endif
          dumy1=dble(u(1)*u(1)+v(1)*v(1)+w(1)*w(1))
          dumy2=dble(u2(1)*u2(1)+v2(1)*v2(1)+w2(1)*w2(1))
          if(aks(alpl(1))/=0.d0.and.dumy1==0.d0) then
           write(ifle,'(1x,a,I4)') 'EER: Boundary no. ',no
           call FFRABORT(1,'ERR:Volt Fraction and u,v,w in INLET')
          elseif(aks(alpl(2))/=0.d0.and.dumy2==0.d0) then
           write(ifle,'(1x,a,I4)') 'EER: Boundary no. ',no
           call FFRABORT(1,'ERR-1:Volt Fraction and u2,v2,w2 in INLET')
          elseif(aks(alpl(1))==1.d0.and.dumy2/=0.d0) then
           write(ifle,'(1x,a,I4)') 'EER: Boundary no. ',no
           call FFRABORT(1,'ERR-2:Volt Fraction and u2,v2,w2 in INLET')
          elseif(aks(alpl(2))==1.d0.and.dumy1/=0.d0) then
           write(ifle,'(1x,a,I4)') 'EER: Boundary no. ',no
           call FFRABORT(1,'ERR:Volt Fraction and u,v,w in INLET')
          endif
        endif
      endif
!
      call strdata(nb,1)
      if( ierror.ne.0 ) goto 9999
      if(E2P) then
        call strdata2(nb,1)
        if( ierror.ne.0 ) goto 9999
      endif
!      
      goto 1000
!--------------------------------
!--< 2.? touch inlet boundary >--
!--------------------------------
 1104 continue
!
!
      if(name2.eq.' ') then
        write(ifle,'(a,I4,2a)') 'EER: Boundary no. ',no,
     &  ' is touch-inlet boundary',
     &  ' its touching BC name must be specified by name2'
        write(ifle,'(a)') 
     &  ' Input the name corresponding to gridfile'
        call FFRABORT(1,'module_boundary/inputdata')
      elseif(adjustl(name2).eq.adjustl(name)) then
        write(ifle,'(2a)') 'EER: Boundary name2: ',adjustl(name2)
        write(ifle,'(a)') 'ERR: BC can not touch itself '
        call FFRABORT(1,'module_boundary/inputdata')
      endif
!
      kdbcnd(0,nb)=kdtchi
      call romtrx(nb)
      if(discont==0.or.discont==1) then
        idis(nb)=discont
      else
         write(ifle,'(a)') 
     &     'EER: [discont] MUST be 0 or 1'
         write(ifle,'(a)') 
     &     'MSG: discont=0 : continuity touch-inlet BC'
         write(ifle,'(a)') 
     &     'MSG: discont=1 : discontinuity touch-inlet BC'
        call FFRABORT(1,'module_boundary/inputdata')
      endif
!
!      masflg(nb)=flux_balance
!      imasflg=imasflg.or.flux_balance==1
!      if(flux_balance==1) then
!        if(my_rank==0) then
!          write(ifll,'(a)') 
!     &    'MSG: Set [massflux] for [touch-inlet] BC'
!          write(ifll,'(3a)') 
!     &    'MSG: (1)massflux>0 : inflow; ',
!     &    'MSG: (2)massflux<0 : outflow ',
!     &    'MSG: (3)massflux=0 : by Up-flow '
!          if(massflux==undef) then
 !           write(ifll,'(a)') 
!     &      'ERR: Set [massflux] for [touch-inlet] BC'
!            write(ifll,'(a)') 
!     &      'MSG: Set [massflux] for [touch-inlet] BC'
!            write(ifll,'(2a)') 
!     &      'MSG: (1)massflux>0 : inflow; ',
!     &      'MSG: (2)massflux<0 : outflow ',
!     &      'MSG: (3)massflux=0 : by Up-flow '
!            call FFRABORT(1,'ERR: NOT defined [massflux]')
!          else
!          endif
!        endif
!        masflx(nb)=-massflux
!      endif
!
      goto 1000
!
!--< 2.6 other boundary >--
!------------------
! --- firewall
!------------------
 1105 continue
      if(VOF_ang/=undef) then
        vofang(nb)=VOF_ang
      else
        vofang(nb)=0.d0
      endif
!
!
      if(kind.eq.firewall) then
        call nml_listno(4,xknd,conc,kcnc) 
        if( kcnc.gt.0 ) conc=xknd(kcnc) 
        if(conc/=neumann) then 
          write(ifle,'(1x,a,I4)') 
     &  " ### EER: firewall must be use [conc='Neumann'] on BC no=",
     &    no
          call FFRABORT(1,'Stop at module_boundary')
        endif 
        kdbcnd(0,nb)=kdfire
        firewal=.true. 
        call chkntime 
        call chkflt1(0,ntime,'t',t) 
        call chkfire(0,ncomp,ntime,'ys',ys) 
        call chkfltn('p(1)',p(1)) 
        if(surface_reaction==1) then
          call surface_reac(ifle,lenspc,
     &    nsite,nbulk,nphasx,nneq,nb,no,ntime,nbcnd,
     &    lsuf,lfc,
     &    cntlnam,
     &    mcomp,mbulk,msite,mneq,nphasex,
     &    suf_reac_no,
     &    bulk_name,site_name,
     &    site_species_name,bulk_species_name,
     &    site_dens,bulk_thick,
     &    site_species_molefraction,
     &    bulk_species_molefraction,
     &    bulk_species_dens,
     &    mallcomp,nallcomp,
     &    phs_dns,phs_thk,
     &    phs_nam,phs_snm,
     &    phs_idx,phs_nbc,phs_typ,phs_com,phs_comr,
     &    heat_release,heat,lcvdwal,surfreac
     &    )
        endif
      endif
!------------------
! --- surface
!------------------
      if(kind.eq.cvdsurface) then 
        kdbcnd(0,nb)=kdcvd 
        call surface_reac(ifle,lenspc,
     & nsite,nbulk,nphasx,nneq,nb,no,ntime,nbcnd,
     & lsuf,lfc,
     & cntlnam,
     & mcomp,mbulk,msite,mneq,nphasex,
     & suf_reac_no,
     & bulk_name,site_name,
     & site_species_name,bulk_species_name,
     & site_dens,bulk_thick,
     & site_species_molefraction,
     & bulk_species_molefraction,
     & bulk_species_dens,
     & mallcomp,nallcomp,
     & phs_dns,phs_thk,
     & phs_nam,phs_snm,
     & phs_idx,phs_nbc,phs_typ,phs_com,phs_comr,
     & heat_release,heat,lcvdwal,surfreac
     & )
        call nml_listno(4,xknd,temp,ktmp) 
        if(ktmp.gt.0) temp=xknd(ktmp)
        if(temp/=neumann) then
!          write(ifle,*)
!     &  " ### ERR: cvdsurface must be use [temp='Neumann']"
!          call FFRABORT(1,'Stop at module_boundary')
        endif
!
        call nml_listno(4,xknd,conc,kcnc)
        if( kcnc.gt.0 ) conc=xknd(kcnc)
        if(conc/=neumann) then
          write(ifle,'(2a)')
     &  " ### ERR: cvdsurface BC must be use [conc='Dirichlet'] in ",
     &    trim(cntlnam)
          call FFRABORT(1,'Stop at module_boundary')
        endif
        call chkntime
!        call chk_t_neumann(1,ntime,'t',t)
        call chkfire(0,ncomp,ntime,'ys',ys) !Ys=0 
        call chkfltn('p(1)',p(1)) 
!
        
      endif
!------------------
! --- buffle BC
!------------------
      if(kind==buffle) then
        if(porous_BC==1) then
          poro(nb)=1
          ical_poroBC=1
        endif
        ical_buff=1
        kdbcnd(0,nb)=kdbuff
        if(name2==' ') then
          write(ifle,*) ' EER: Boundary no. ',no,
     &    ' is [buffle] boundary',
     &   ' pair BC name specification by [name2] is under-construction'
          write(ifle,*) 'In ',cntlnam,
     &    ' Remove [name2]. May be, you can not use buffle.'
          call FFRABORT(1,'boundary/module')
        endif
        if(surface_reaction==1) then
          call surface_reac(ifle,lenspc,
     &    nsite,nbulk,nphasx,nneq,nb,no,ntime,nbcnd,
     &    lsuf,lfc,
     &    cntlnam,
     &    mcomp,mbulk,msite,mneq,nphasex,
     &    suf_reac_no,
     &    bulk_name,site_name,
     &    site_species_name,bulk_species_name,
     &    site_dens,bulk_thick,
     &    site_species_molefraction,
     &    bulk_species_molefraction,
     &    bulk_species_dens,
     &    mallcomp,nallcomp,
     &    phs_dns,phs_thk,
     &    phs_nam,phs_snm,
     &    phs_idx,phs_nbc,phs_typ,phs_com,phs_comr,
     &    heat_release,heat,lcvdwal,surfreac
     &    )
          surfreac(nb)=2
        endif
        goto 1000
      endif
!------------------------------
! --- Porous BC
!------------------------------
      if(kind==porous) then
        call FFRABORT(1,'ERR: Using [buffle] BC for porous BC')
        poro(nb)=1
        ical_poroBC=1
        kdbcnd(0,nb)=kdpors
        if(C2_porous==undef.or.perm_porous==undef.or.
     &     thickness_porous==undef.or.porosity==undef) then
          write(ifle,'(1X,A,1X,A)') 
     &  'ERR: NOT define C2_porous, perm_porous porosity ',
     &  "or  thickness_porous for Darcy s law"
          call FFRABORT
     &   (1,'ERR: NOT define C2_porous, perm_porous')
        else
          C2pBC(nb)=C2_porous
          alpha_coBC(nb)=perm_porous
          thkposBC(nb)=thickness_porous
          prosty(nb)=porosity
        endif
!
        if(surface_reaction==1) then
          call surface_reac(ifle,lenspc,
     &    nsite,nbulk,nphasx,nneq,nb,no,ntime,nbcnd,
     &    lsuf,lfc,
     &    cntlnam,
     &    mcomp,mbulk,msite,mneq,nphasex,
     &    suf_reac_no,
     &    bulk_name,site_name,
     &    site_species_name,bulk_species_name,
     &    site_dens,bulk_thick,
     &    site_species_molefraction,
     &    bulk_species_molefraction,
     &    bulk_species_dens,
     &    mallcomp,nallcomp,
     &    phs_dns,phs_thk,
     &    phs_nam,phs_snm,
     &    phs_idx,phs_nbc,phs_typ,phs_com,phs_comr,
     &    heat_release,heat,lcvdwal,surfreac
     &    )
          surfreac(nb)=1
        endif
        goto 1000
      endif


!------------------
! --- interface 
!------------------
      if(kind.eq.interface) then
        kdbcnd(0,nb)=kdintr
        if(name2.eq.' ') then
          write(ifle,*) ' EER: Boundary no. ',no,
     &    ' is [interface] boundary',
     &    ' its pair BC name must be specified by [name2]'
          write(ifle,*) 'In ',cntlnam,
     &    ' Set [name2] , corresponding to grid file'
          call FFRABORT(1,'boundary/module')
        endif
        if(discont==0.or.discont==1) then 
          idis(nb)=discont
!          if(discont==1) then 
!            call FFRABORT
!     &   (1,'ERR: NOT support discontinuity interface BC')
!          endif
        else
          write(ifle,'(a)') 'EER: [discont] MUST be 0 or 1'
          write(ifle,'(a)') 
     &     'MSG: discont=0 : continuity periodic BC'
          write(ifle,'(a)') 
     &     'MSG: discont=1 : discontinuity periodic BC'
          call FFRABORT(1,'module_boundary/inputdata')
        endif
        if(surface_reaction==1) then
          call surface_reac(ifle,lenspc,
     &    nsite,nbulk,nphasx,nneq,nb,no,ntime,nbcnd,
     &    lsuf,lfc,
     &    cntlnam,
     &    mcomp,mbulk,msite,mneq,nphasex,
     &    suf_reac_no,
     &    bulk_name,site_name,
     &    site_species_name,bulk_species_name,
     &    site_dens,bulk_thick,
     &    site_species_molefraction,
     &    bulk_species_molefraction,
     &    bulk_species_dens,
     &    mallcomp,nallcomp,
     &    phs_dns,phs_thk,
     &    phs_nam,phs_snm,
     &    phs_idx,phs_nbc,phs_typ,phs_com,phs_comr,
     &    heat_release,heat,lcvdwal,surfreac
     &    )
        endif
      endif
!
      if(kind==shutter) then 
!        if(shutter_start==undef) then
!          call FFRABORT(1,'ERR: NOT been defined ""shutter_start in fflow.ctl')
!        elseif(shutter_start<0.d0)  then
!          shut_strt(nb)=0.d0
!        else
!          shut_strt(nb)=shutter_start
!        endif
!        if(shutter_speed==undef) then
!        endif
        kdbcnd(0,nb)=kdintr 
        if(name2.eq.' ') then
          write(ifle,*) ' EER: Boundary no. ',no,
     &    ' is [shutter] boundary',
     &    ' its pair BC name must be specified by [name2]'
          write(ifle,*) 'In ',cntlnam,
     &    ' Set [name2] , corresponding to grid file'
          call FFRABORT(1,'boundary/module')
        endif 
      endif 
!----------------------------
! / kind of each variable /
!----------------------------
      call nml_listno(7,vknd,vel ,kvel)
      call nml_listno(4,xknd,temp,ktmp)
      call nml_listno(4,xknd,conc,kcnc)
      call nml_listno(3,rknd,rans,krns)
!
!      if(imhd>0) then
!        call nml_listno(2,mhdknd,MHD_A,kmhd_A)
!        call nml_listno(2,mhdknd,MHD_F,kmhd_F)
!      endif
!
      call nml_chkchr0(ifle,'temp',temp,ktmp,ierr1)
      if(ierr1.ne.0) goto 9999
      if(lrans) then
      else
        rans=' '
      endif
!
      if( kvel.gt.0 ) vel =vknd(kvel)
      if( ktmp.gt.0 ) temp=xknd(ktmp)
      if( kcnc.gt.0 ) conc=xknd(kcnc)
      if( krns.gt.0 ) rans=rknd(krns)
!
!      if(imhd>0) then
!        if( kmhd_A.gt.0 ) MHD_A=mhdknd(kmhd_A)
!        if( kmhd_F.gt.0 ) MHD_F=mhdknd(kmhd_F)
!      endif
!
      if( vel.eq.noslip )     kdbcnd(1,nb)=kvnslp
      if( vel.eq.freeslip )   kdbcnd(1,nb)=kvfslp
      if( vel.eq.loglaw )     kdbcnd(1,nb)=kvlglw
      if( vel.eq.rough  )     kdbcnd(1,nb)=kvrogf
      if( vel.eq.model  )     kdbcnd(1,nb)=kvmodl
      if( vel.eq.EWF  )     then
        call FFRABORT(1,'ERR: vel=EWF is NOT support')
        kdbcnd(1,nb)=kvEWF
      endif
      if( vel.eq.satawall  )  kdbcnd(1,nb)=kvsataW
      IF(vel.eq.rough.or.vel.eq.satawall) then
        if(roughness==undef) then
          if(vel.eq.rough) then
            write(ifle,'(1X,a)') 
     &      'WRN: NOT Define [roughness] while using [rough-wall]'
            write(ifle,'(1X,2a)') 
     &      'MSG: Default value [roughness=3.52D0 (~)]',
     &      'MSG: [roughness]=Ks*Ut/Nu'
            rghnes(nb)=3.52D0
          elseif(vel.eq.satawall) then
            write(ifle,'(1X,a)') 
     &      'WRN: NOT Define [roughness] while using [sata-wall]'
            write(ifle,'(1X,2a)') 
     &      'MSG: Default value [roughness=1.3D-5 (m)]',
     &      'MSG: [roughness]=Ks*Ut/Nu'
            rghnes(nb)=1.3D-5  ![m]
          endif
        endif
        if((roughness<3.52D0.or.roughness>200.d0)
     &      .and.kdbcnd(1,nb)==kvlglw) then
          call FFRABORT
     &  (1,'ERR:[roughness] : roughness>=3.52D0 or roughness<=200.d0')
        endif
        rghnes(nb)=roughness
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
      if( rans.eq.dirichlet ) then
!         if(my_rank==0) then
!           write(ifll,'(a)')
!     & 'MSG: k & e scalars are alway set in [loglaw] in High Re model'
!           write(ifll,'(a)')
!     & 'MSG: k & e scalars are alway set in default in Low Re model'
!         endif
         kdbcnd(4,nb)=kkdirc
      endif
      if( rans.eq.neumann   ) then
!         if(my_rank==0) then
!           write(ifll,'(a)')
!     & 'MSG: k & e scalars are alway set in [loglaw] in High Re model' 
!           write(ifll,'(a)')
!     & 'MSG: k & e scalars are alway set in default in Low Re model'
!         endif
         kdbcnd(4,nb)=kkneum
      endif
      if(rans.eq.loglaw) then
         if(kemdl.and.my_rank==0) then
           write(ifll,'(a)')
     &     'MSG: High Reynolds model , wall function (log-law)'
         elseif(lowke.and.my_rank==0) then
           write(ifll,'(a)')
     &     'MSG: Low Reynolds model , eps value alg for Near wall cell'
           write(ifll,'(a)') 'MSG: Select dirichlet or neumann'
         else
         endif
         kdbcnd(4,nb)=kklglw
      endif
!
!------------------------
! / check & store data /
!------------------------
!
      if(kind.eq.interface) then
        if( temp.eq.dirichlet ) then
          write(*,*) '### error : data error'
          write(*,*) 'temp =',temp
          write(*,*) 'it can not be = ',temp
          goto 9999
        endif
      endif
!
      call chkntime
!
      if(vel.eq.noslip.or.vel.eq.loglaw.or.vel==rough.or.
     &  vel==model.or.vel==satawall)
     &then
        if(turn_wall/=1) then
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
      else  ! slip wall BC
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
      if( kind.eq.interface.and.temp/=neumann ) idum=-1
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
      if(conc.eq.dirichlet .or. conc.eq.transfer ) then
        yfix=1
      elseif( conc.eq.neumann ) then
        yfix=0
      else
        yfix=-1
      endif
      if( yfix.lt.0 ) then
        call chkfltn2('ys(1)',ys(1))
        if(E2P) then
          call chkfltn2('ys2(1)',ys2(1))
        endif
      else
        call chkflt2(yfix,ncomp,ntime,'ys',ys)  !neumann yfix=0
        if(E2P) then
          call chkflt2(yfix,ncomp,ntime,'ys2',ys2)
        endif
      endif
!
      if(rans.eq.dirichlet ) then
        call chkflt2(1,nrans,ntime,'aks',aks)
        if(E2P) then
          aph=aks(alpl(1))+aks(alpl(2))
          if(abs(aph-1.d0).gt.SML1) then
            write(ifle,*)
     &      '### error : alph1+alph2 must be equal to 1.0 at BC'
            write(ifll,'(1X,A)') 'alph1/alph2 is reset Neumann BC '
            write(ifle,*) '  MSG : BC name :', TRIM(adjustl(name))
!            write(ifle,*) '  alph1+alph2= ',aph
!            call FFRABORT(1,'module_boundary/namelist &boundary')
          endif
        endif
      elseif(rans.eq.neumann) then
!        call chkfltn2('aks(1)',aks(1))
        aks=0.d0
      elseif(rans.eq.loglaw) then
        call chkfltn2('aks(1)',aks(1))
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
      if(ierror.ne.0) goto 9999
      if(E2P) then
        call strdata2(nb,yfix)
        if(ierror.ne.0) goto 9999
      endif
      goto 1000
!---------------------------------------
! --- pressure and stagnation pressure
!---------------------------------------
 1106 continue


      call chkntime
      call chkflt1(1,ntime,'t',t)
      call chkflt2(1,ncomp,ntime,'ys',ys)
      call chkflt1(1,ntime,'p',p)
!
      if(kind.eq.static) then
! --- t=> static temperature
! --- p=> static pressure
        kdbcnd(0,nb)=kdpres
      elseif(kind.eq.stagnation) then
! --- t=> static temperature
! --- p=> stagnation pressure
        kdbcnd(0,nb)=kdstag
!        call FFRABORT(1,'stagnation pressure BC not finished')
      endif
      call strdata(nb,1)
      goto 1000
!
!---------------------------------------
! --- Sliding BC
!---------------------------------------
!
 1107 continue
!
      kdbcnd(0,nb)=kdsld
      lsldbc=.true.
!
      nblade(nb,1)=1
      if(n_blade1<0) then
        write(ifle,*) ' ERR:n_blade1 must be lager then 0'
        write(ifle,*) ' MSG:Geom model is 1/n_blade1'
        call FFRABORT(1,'ERR:')
      elseif(n_blade1/=0) then
        nblade(nb,1)=n_blade1
        if(n_blade1>1) pitch_one=1
      endif
!
      nblade(nb,2)=1
      if(n_blade2<0) then
        write(ifle,*) ' ERR:n_blade2 must be lager then 0'
        write(ifle,*) ' MSG:Geom model is 1/n_blade2'
        call FFRABORT(1,'ERR:')
      elseif(n_blade2/=0) then
        nblade(nb,2)=n_blade2
        if(n_blade2>1) pitch_one=1
      endif

      if(name2.eq.' ') then
        write(ifle,'(a,I4,2a)') 'EER: Boundary no. ',no,
     &  ' is [sliding] boundary,',
     &  ' In Sliding BC,2nd BC must be specified by name2'
        write(ifle,'(a)') 
     &    'MSG: Input the name corresponding to grid file'
        call FFRABORT(1,'ERR: &boundary')
      endif    
      call romtrx(nb)
      if(discont==1) then
        idis(nb)=discont
      elseif(discont==0) then
        write(ifle,'(a)') 'EER: [discont] MUST be 1 for sliding BC'
!        idis(nb)=discont
        call FFRABORT(1,'ATT: re-run [prefflow]')
      else
         write(ifle,'(a)') 'EER: [discont] MUST be 0 or 1'
         write(ifle,'(a)') 
     &   'MSG: discont=0 : continuity sliding BC at initial condition'
         write(ifle,'(a)') 
     &'MSG: discont=1 : discontinuity sliding BC at initial condition'
        call FFRABORT(1,'module_boundary/inputdata')
      endif
      if(nblade(nb,1)>1.or.nblade(nb,2)>1) then
        idis(nb)=2
      endif
!
      goto 1000
!
 1108 continue

      lovstbc=.true.
      if(kind==overset) then
        kdbcnd(0,nb)=kdovst
      endif
      goto 1000
!
 1001 continue
!
      if(lFC>0) then
        if(FC_V_flag==0) then
          call FFRABORT(1,'ERR: V_cell NOT defined for PEFC')
        endif
      endif
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
!----------------------------
!-< 4. Initialize "wdbcnd" >-
!----------------------------
      do 400 n=1,nbcnd
      kt=iwbcnd(n-1)+1
      do 401 i=1,nbvals
      wdbcnd(i,n)=wdtbcnd(i,kt)
  401 continue
  400 continue
!
      if(E2P) then
        do n=1,nbcnd
        kt=iwbcnd(n-1)+1
        do i=1,nbvals
        wdbcnd2(i,n)=wdtbcnd2(i,kt)
        enddo
        enddo
      endif
!
!------------------------
! --- check sliding BC --
!------------------------
!
      if(lrot.and.(lsldbc.or.lovstbc)) then
        ical_sld=1
      endif
      if((lsldbc.or.lovstbc).and.(.not.lrot)) then
        ical_sld=2   !rpm=0 but have angle-rotation
        write(ifle,'(a)') 'WRN: Not define [rpm] '
        write(ifle,'(a)') 
     &   'WRN: [rpm] must be defined while using sliding BC'
      endif
      if((.not.lsldbc.and..not.lovstbc).and.lrot) then
        ical_sld=3
        write(ifle,'(a)')
     &   'WRN: Sliding BC should be defined while rpm is defined'
      endif
      if(pitch_one==1) then
        ical_sld=4
      endif
!
!-------------------------
! --- check touch inlet --
!-------------------------
!
      masbc(:)=0
      do n=1,nbcnd
      if(kdbcnd(0,n)==kdtchi) then
        do nb=1,nbcnd
        if(adjustl(boundName(n,2))==adjustl(boundName(nb,1))) then 
          masbc(n)=nb 
          exit
        endif
        enddo
        IF(masbc(n)==0) THEN
          call FFRABORT(1,'ERR: touch inlet BC: name2')
        endif
      endif
      enddo
!nph
!----------------------------------------------------------------------
! --- check cvd-wall and surface-reaction
!----------------------------------------------------------------------
      if((.not.lsuf).and.lcvdwal) then
        ical_suf=2
	write(ifle,'(2a)') 
     &  'WRN: No surface reaction defined in ',cntlnam
              write(ifle,'(3a)') 
     &  "MSG: set [&chemreac/kind='gas'] ",
     &  "for gas phase reaction in ",cntlnam
        write(ifle,'(3a)') 
     &  "MSG: set [&chemreac/kind='surface'] ",
     &  "for surface reaction in ",cntlnam
      endif
!
      if(lsuf.and.(.not.lcvdwal)) then
	write(ifle,'(2a)') 
     &  'ERR: No CVD wall defined in ',cntlnam
        write(ifle,'(2a)') 
     &  "MSG: [kind='cvdsurface'] in ",cntlnam
        write(ifle,'(a)') 
     &  "MSG: Re-run [prefflow]"
!        call FFRABORT(1,'Stop at module_boundary')
      endif
!
      if((.not.lcvdwal.and.KSUF==1).or.(lcvdwal.and.KSUF==0)) then
	write(ifle,'(a)') 
     &    'ERR : Surface reaction dimension error'
        write(ifle,'(a)') 'MSG : Re-run prefflow'
        call FFRABORT(1,'module_boundary/inputdata')
      endif
!
!----------------------------------------------------------------
! --- Surface reaction: List every Phase (gas, n-site, m-bulk)
!----------------------------------------------------------------
!
      if(lcvdwal.and.lsuf) then
        do i=1,nphasx
        phs_idx(i)=phs_idx(i)+phs_idx(i-1)
        enddo
!
        phs_com=0
        phs_comr=0
        do i=1,ncmpx
          phs_com(i)=i
          phs_comr(i)=i
        enddo
        do i=2,nphasx
        do j=phs_idx(i-1)+1,phs_idx(i)
        do k=ncmpx+1,ncmpx+ncomp_sufx   !(nallcomp/=ncmpx+ncomp_sufx)
          if(trim(phs_snm(j))==trim(spinam(k))) then
            if(phs_com(j)/=0.and.phs_com(j)/=k) then
              write(ifle,*)
     &        'ERR-1: Duplicated species name :',trim(spinam(k))
              call FFRABORT
     &        (1,'ERR: Duplicated species name in &surface_species')
            endif
            phs_com(j)=k
            phs_comr(k)=j
          endif
        enddo
        enddo
        enddo
!
        if(my_rank==0) then
          write(ifll,4998)
          do i=1,nphasx
          write(ifll,4998)
          write(ifll,'(2X,a,2X,a,I4)') 'Phase Information :',
     &     'Phs_no=',i
          write(ifll,4997)
     &    'BC_no=',nobcnd(phs_nbc(i)),
     &    'Phs_type=',phstyp(phs_typ(i)),'Phs_species_nu=',
     &    phs_idx(i)-phs_idx(i-1),'Phs_name=',trim(phs_nam(i))
! --- 
          do j=phs_idx(i-1)+1,phs_idx(i)
          write(ifll,4996) 'Phs_species_no= ',phs_com(j),
     &    'SIGAMA=',num_site(j),
     &    ' Phs_species_name= ',trim(phs_snm(j))
          enddo
          enddo
          write(ifll,4999)
          write(ifll,4995)
          write(ifll,4994)
!
        endif
!
        do j=2,nphasx
        do k=phs_idx(j-1)+1,phs_idx(j)
        do i=k+1,phs_idx(j)
        IF(trim(phs_snm(k))==trim(phs_snm(i))) then
          write(ifle,*)
     &    'ERR: Duplicated species name on BC:',trim(phs_snm(i))
          write(ifle,*)
     &    'MSG: Duplicated species name is ',
     &    'NOT allowed between phases in same substrate'
          call FFRABORT
     &    (1,'ERR: Duplicated species name in &boundary')
        endif
        enddo
        enddo
        enddo
!
       do k=ncmpx+1,ncmpx+ncomp_sufx  !nallcomp
          do i=k+1,ncmpx+ncomp_sufx   !nallcomp
              if(trim(spinam(k))==trim(spinam(i))) then
                write(ifle,*)
     &         'ERR-2: Duplicated species name:',trim(spinam(i))
                call FFRABORT
     &         (1,'ERR: Duplicated species name in &surface_species')
              endif
          enddo
       enddo
!
       do i=1,nphasx
         do j=phs_idx(i-1)+1,phs_idx(i)
          if(phs_com(j)==0) then
            write(ifle,'(3a)')
     &     'ERR: Speciee name defined in &boundary: ',trim(phs_snm(j)),
     &      ' is NOT defined in &surface_species'
            call FFRABORT(1,'ERR: Speciee name matching error')
          endif
         enddo
       enddo
!
       do i=ncmpx+1,ncmpx+ncomp_sufx  !nallcomp
         if(phs_comr(i)==0) then
           write(ifle,'(3a)')
     &     'ERR: Speciee name defined in &surface_species: ',
     6     trim(spinam(i)),
     &     ' is NOT defined in &boundary'
           call FFRABORT(1,'ERR: Speciee name matching error')
         endif
       enddo
!-----------------------------------------------------------------------
! --- calculation of index number m  !cofleft(ns*(q-1)+j)=cleft(s)
! --- sum of all the stoichiometric coefficients of reactants in suface
!-----------------------------------------------------------------------
!
        suf_dns(:)=0.d0
        n_bi(:,:)=1.d0         ! n_bi(nbcnd,iq)
        n_site(:,:)=1.d0       ! n_site(nbcnd,iq)

        do nb=1,nbcnd
        kd=kdbcnd(0,nb)
!        if(kd==kdcvd) then
        if(surfreac(nb)==1) then
          do i=1,nphasx
            ity=phs_typ(i)
            nbk=phs_nbc(i)
            if(nbk==nb.and.ity==surphase) then
              suf_dns(nb)=suf_dns(nb)+phs_dns(i)
            endif
          enddo
!
          do iq=1,nneq
          sigma_site(:)=0d0   ! sigma_site(nphasx)
          bi_site(:)=1.d0     ! bi_site(nphasx)
          if(IDX_SUFRAC(nb,iq)) then
            do i=1,nphasx
              ity=phs_typ(i)
              nbk=phs_nbc(i)
              if(ity==surphase.and.nbk==nb) then
                do j=phs_idx(i-1)+1,phs_idx(i)
                k=phs_com(j)
                n_bi(nb,iq)=n_bi(nb,iq)*
     &                  num_site(j)**cofleft(ncompall*(iq-1)+k)

                sigma_site(i)=sigma_site(i)
     &                          +net_vreq(ncompall*(iq-1)+k)
                bi_site(i)=bi_site(i)*
     &             num_site(j)**(-net_vreq(ncompall*(iq-1)+k))
                enddo
              endif
            enddo
            do i=1,nphasx
              ity=phs_typ(i)
              nbk=phs_nbc(i)
              if(ity==surphase.and.nbk==nb) then
                do j=phs_idx(i-1)+1,phs_idx(i)
                k=phs_com(j)
                n_site(nb,iq)=n_site(nb,iq)*bi_site(i)
     &                    *(phs_dns(i)**sigma_site(i))
                enddo
              endif
            enddo
!
            rdum=0.d0
            do i=1,nphasx
            ity=phs_typ(i)
            nbk=phs_nbc(i)
            if(ity==surphase.and.nbk==nb) then
               do j=phs_idx(i-1)+1,phs_idx(i)
                 k=phs_com(j)
                 rdum=rdum+net_vreq(ncompall*(iq-1)+k)*num_site(j)
               enddo
            endif
            enddo
            if(abs(rdum)>SML1) then
               write(ifll,'(2x,a,I4)')
     &  'WRN: SITE conservation is NOT keeped in reaction no :',iq
!               call FFRABORT(1,'ERR: ???')
            endif
          endif
          enddo
        endif
        enddo
!
!
          m_indx(:)=0.d0
          do iq=1,nneq
          do k=ncmpx+1,ncmpx+ncomp_sufx
          m_indx(iq)=m_indx(iq)+cofleft(ncompall*(iq-1)+k)
          enddo
          enddo

      endif
!
!
!
      blktyp(:)=0
      do i=1,nphasx
        ity=phs_typ(i)
        if(ity==gasphase) then
          blktyp(1:ncmpx)=1
        elseif(ity==surphase) then
          do j=phs_idx(i-1)+1,phs_idx(i)
          k=phs_com(j)
          blktyp(k)=2
          enddo
        elseif(ity==blkphase) then
          do j=phs_idx(i-1)+1,phs_idx(i)
          k=phs_com(j)
          blktyp(k)=3
          enddo
        endif
        enddo
!
 4995 format
     & (2X,'Phs_type=1: Gas-Phase; ',
     &  'Phs_type=2: Surface-Phase; Phs_type=3: Bulk-Phase')
 4996 format(2X,a,I4,' | ',a,F8.2,' | ',a,a)
 4997 format(2X,a,I4,/,2x,a,a,/,2x,a,I4,/,2x,a,a)
 4998 format(2X,108('-'))
 4999 format(2X,108('='))
 4994 format(2X,108('|'))
!
!----------------------------------------------------------------------
!
      return
!
 9999 continue
      write(ifle,*) modnam//subnam
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #6.1 check real data >
!=================================================
      subroutine chkflt1(isw,nt,dnam,dval)
!=================================================
      implicit none
      integer     ,intent(in) :: isw,nt
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval(nt)
      integer :: i
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
      if( ( isw.eq.1 .and. dval(i).le.0.d0 ) .or.
     &    ( isw.eq.2 .and. dval(i).lt.0.d0 ) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*)  dnam,'(i) = ',dval(i)
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
      real(8),parameter        :: SML=1.d-15
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
!< #6.2 check real data >
!=================================================
      subroutine chkflt2(isw,nc,nt,dnam,dval)
!=================================================
      implicit none
      integer     ,intent(in) :: isw,nc,nt
      character(*),intent(in) :: dnam
      real*8      ,intent(in) :: dval(nc,nt)
      integer :: i,j
      do 100 i=1,nc
      do 100 j=1,nt
      if( dval(i,j).eq.dundef ) then
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
        write(ifle,*)  dnam,'(i,j) = ',dval(i,j)
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
      implicit none
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
      implicit none
      character(*),intent(in) :: dnam
      real        ,intent(in) :: dval
!
      if( dval.eq.undef ) return
      write(ifle,*) '### error-1 : data error'
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) dnam,' = ',dval
      write(ifle,*) 'it can not be specified'
      ierror=1
      end subroutine chkfltn
!

!=================================================
      subroutine chkfltn2(dnam,dval)
!=================================================
      implicit none
      character(*),intent(in) :: dnam
      real*8      ,intent(in) :: dval
!
      if( dval.eq.dundef ) return
      write(ifle,*) '### error-1 : data error'
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      write(ifle,*) dnam,' = ',dval
      write(ifle,*) 'it can not be specified'
      ierror=1
      end subroutine chkfltn2
!< #6.5 check integer data >
!=================================================
      subroutine chkintn(dnam,dval)
!=================================================
      implicit none
      character(*),intent(in) :: dnam
      integer     ,intent(in) :: dval
      if( dval.eq.iundef ) return
      write(ifle,*) '### error-2 : data error'
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
      real*8      ,intent(in) :: dval(nc,nt)
      integer :: i,j
      do 100 i=1,nc
      do 100 j=1,nt
      if( dval(i,j).eq.dundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'BC no. =',nobcnd(nb)
        write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
        write(ifle,*) 'lack of data : ',dnam,'(i,j)'
        write(ifle,*) 'i,j =',i,j
        goto 9999
      endif
  100 continue
      if(isw.ne.0 ) return
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
      implicit none
      integer :: i
      real    :: dum1
      if( ntime.lt.1 ) then
        write(ifle,*) '### error : data error'
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
      implicit none
      integer,intent(in) :: nb,yfix
      integer :: i,j,i0,j0,k0
      real(8)  :: ysw(mcomp),sum,fmax
                           wdbcnd(0,nb)=-1.d0
      if( cycle.ne.undef ) wdbcnd(0,nb)=dble(cycle)
      if( htc  .ne.undef ) adbcnd(1,nb)=dble(htc)
      if( mtc  .ne.undef ) adbcnd(2,nb)=dble(mtc)
      adbcnd(3,nb)=-1.d0
      if(vsgm(1).ne.undef ) then
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
      if( ys(1).eq.dundef ) goto 109
      j0=ncomp*(i-1)
      sum=0.d0
      if(yfix.gt.0) then
        do 200 j=1,ncomp
        ysw(j)=max(0.d0,ys(j+j0))   !dble(max(0.d0,ys(j+j0)))
        sum=sum+ysw(j)
  200   continue
        if( sum.le.0.d0 ) goto 9001
        sum=1.d0/sum
      else
        fmax=-1.d0
        do 210 j=1,ncomp
        ysw(j)=ys(j+j0)
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
      if(aks(1).ne.dundef) then
        do 102 j=1,nrans
        wdtbcnd(israns+j,i+i0)=aks(j)
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
!=================================================
      subroutine strdata2(nb,yfix)
!=================================================
      implicit none
      integer,intent(in) :: nb,yfix
      integer :: i,j,i0,j0,k0
      real(8)  :: ysw(mcomp),sum,fmax
                           wdbcnd2(0,nb)=-1.d0
      if( cycle.ne.undef ) wdbcnd2(0,nb)=dble(cycle)
      if( htc2  .ne.undef ) adbcnd2(1,nb)=dble(htc2)
      if( mtc2  .ne.undef ) adbcnd2(2,nb)=dble(mtc2)
      adbcnd2(3,nb)=-1.d0
      if( vsgm2(1).ne.undef ) then
        adbcnd(3,nb)=dble(vsgm2(1))
        adbcnd(4,nb)=dble(vsgm2(2))
      endif
      i0=iwbcnd(nb-1)
      do 100 i=1,ntime
                          wdtbcnd2(0,i+i0)=dble(time(i))
      if( p(1).ne.undef ) wdtbcnd2(1,i+i0)=dble(p(i))
      if( u2(1).ne.undef ) wdtbcnd2(2,i+i0)=dble(u2(i))
      if( v2(1).ne.undef ) wdtbcnd2(3,i+i0)=dble(v2(i))
      if( w2(1).ne.undef ) wdtbcnd2(4,i+i0)=dble(w2(i))
      if( t2(1).ne.undef ) wdtbcnd2(5,i+i0)=dble(t2(i))
      if( ys2(1).eq.dundef ) goto 109
      j0=ncomp*(i-1)
      sum=0.d0
      if( yfix.gt.0 ) then
        do 200 j=1,ncomp
        ysw(j)=max(0.d0,ys2(j+j0))
        sum=sum+ysw(j)
  200   continue
        if( sum.le.0.d0 ) goto 9001
        sum=1.d0/sum
      else
        fmax=-1.d0
        do 210 j=1,ncomp
        ysw(j)=ys2(j+j0)
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
      wdtbcnd2(iscomp+j,i+i0)=ysw(j)*sum
  101 continue
  109 continue
      if( aks(1).ne.dundef ) then
        do 102 j=1,nrans
        wdtbcnd2(israns+j,i+i0)=aks(j)
  102   continue
      endif
  100 continue
      return
 9001 continue
      write(ifle,*) '### error : data error in E2P model'
      write(ifle,*) 'sum of ys2(i-j) = ',sum
      write(ifle,*) 'BC no. =',nobcnd(nb)
      write(ifle,*) 'BC kind :',TRIM(adjustl(kind))
      if( yfix.gt.0 ) write(ifle,*) 'it must be > 0'
      if( yfix.lt.1 ) write(ifle,*) 'it must be = 0'
      write(ifle,*) 'i,j =',j0+1,j0+ncomp
      ierror=1
      end subroutine strdata2

!
!< #6.8 calculate rotational matrix for periodic boundary >
!=================================================
      subroutine romtrx(n)
!=================================================
      implicit none
      integer,intent(in)  :: n
      integer :: i
      real(8)  :: th,aa,ax,ay,az,cost,sint,cosu
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
!nph
!=======================================================
      subroutine surface_reac(ifle,lenspc,
     & nsite,nbulk,nphasx,nneq,nb,no,ntime,nbcnd,
     & lsuf,lfc,
     & cntlnam,
     & mcomp,mbulk,msite,mneq,nphasex,
     & suf_reac_no,
     & bulk_name,site_name,
     & site_species_name,bulk_species_name,
     & site_dens,bulk_thick,
     & site_species_molefraction,
     & bulk_species_molefraction,
     & bulk_species_dens,
     & mallcomp,nallcomp,
     & phs_dns,phs_thk,
     & phs_nam,phs_snm,
     & phs_idx,phs_nbc,phs_typ,phs_com,phs_comr,
     & heat_release,heat,lcvdwal,surfreac
     & ) 
!=======================================================
!------------------------
! --- 
!------------------------
      integer,intent(in) :: nneq,ntime,nb,no,
     &                      nphasex,nbcnd,lfc
      integer,intent(inout) :: nsite,nbulk,nphasx,
     &                         nallcomp
      integer,intent(in) :: ifle,msite,lenspc,mcomp,mneq,
     &                      mbulk,mallcomp
      integer,intent(inout) :: suf_reac_no(mneq)
      integer,intent(inout) :: phs_idx(0:nphasex)
      integer,intent(inout) :: phs_nbc(nphasex)
      integer,intent(inout) :: phs_typ(nphasex)
      integer,intent(inout) :: phs_com(mallcomp), phs_comr(mallcomp)
      integer,intent(inout) :: heat(nbcnd),heat_release
      integer,intent(inout) :: surfreac(nbcnd)
!
      real(8),intent(in) :: site_dens(msite)
      real(8),intent(in) :: bulk_thick(mbulk)
      real(8),intent(in) :: 
     &              site_species_molefraction(msite,mcomp),
     &              bulk_species_molefraction(mbulk,mcomp),
     &              bulk_species_dens(mbulk,mcomp)
      real(8),intent(inout) :: phs_dns(nphasex)
      real(8),intent(inout) :: phs_thk(nphasex)
!
      character(*),intent(inout)  :: site_species_name(msite,mcomp)
      character(*),intent(in)  :: bulk_species_name(mbulk,mcomp)
      character(*),intent(in)  :: site_name(msite)
      character(*),intent(in)  :: bulk_name(mbulk)
      character(*),intent(in) :: cntlnam
      character(*),intent(inout) :: phs_nam(nphasex)
      character(*),intent(inout) :: phs_snm(mallcomp)
      logical,intent(in)  :: lsuf
      logical,intent(inout) :: lcvdwal
      
!
! --- 
!
      integer :: nst,iq,jq,idum,ios=0,nbk
      integer :: i,j,k
      integer :: nm,ncm
      integer,parameter :: iundef=-huge(1)
      real   ,parameter :: undef=-huge(1.)
      character*1  :: slash
!
      if(.not.lsuf) then 
	write(ifle,'(2a)')
     &  ' ### WRN: No surface reaction defined in ',cntlnam
      else
        nsite=0
        do nst=1,msite
        IF(site_name(nst)/='') then
          nsite=nsite+1
          nphasx=nphasx+1
          if(site_dens(nst)==undef) then
            write(ifle,'(a,I4)') 
     &      'ERR* SITE density NOT defined at BC no=',no
            call FFRABORT
     &      (1,'ERR: SITE density [site_dens] must be defined')
          else
            phs_dns(nphasx)=site_dens(nst)
            phs_thk(nphasx)=0.d0
          endif
          phs_nam(nphasx)=site_name(nst)
          phs_nbc(nphasx)=nb
          phs_typ(nphasx)=surphase
          nm=0
          do ncm=1,mcomp
          if(site_species_name(nst,ncm)/=' ') then
            ios=0
            idum=1
            nm=nm+1
            nallcomp=nallcomp+1
            do i=1,lenspc
            slash=site_species_name(nst,ncm)(i:i)
            if(slash=='/') then
              slash=site_species_name(nst,ncm)(i+1:i+1)
              read(slash,*,iostat=ios) idum
              if(ios>0) then
                write(ifle,'(2a)') 
     &             'WRN: site number NOT defined for ',
     &             trim(site_species_name(nst,ncm))
                write(ifle,'(2a)') '
     &             MSG: site number is defined to 1 for ',
     &             trim(site_species_name(nst,ncm))
                idum=1
              endif
              num_site(nallcomp)=dble(idum)
              site_species_name(nst,ncm)(i:lenspc)=''
              exit
            endif
            enddo
            if(nm==0) then
              write(ifle,'(a,I4,a)') 
     &        'ERR: NO site species name are defined at BC & SITE=',no,
     &        trim(site_species_name(nst,ncm))
              call FFRABORT(1,'ERR: NO site species defined')
            endif
            phs_snm(nallcomp)=site_species_name(nst,ncm)
            if(site_species_molefraction(nst,ncm)==undef) then
              call FFRABORT
     &           (1,'ERR: bulk_species_molefraction(:,:) not defined')
            else
              phs_inifrac(nallcomp)=site_species_molefraction(nst,ncm)
            endif
          endif
          enddo
          phs_idx(nphasx)=nm
        endif
        enddo
!
        do iq=1,nneq
        do jq=1,nneq
        if(suf_reac_no(jq)==iundef) cycle
        if(suf_reac_no(jq)>nneq.or.suf_reac_no(jq)<=0) then
          write(ifle,*) ' ### ERR: No such reaction defined'
          call FFRABORT(1,'&boundary/suf_reac_no')
        endif
        if(chem_bc(suf_reac_no(jq))/=isurface) then
          write(ifle,*) ' ### ERR: Reaction= ',suf_reac_no(jq),
     &                  ' is NOT surface reaction'
          call FFRABORT(1,'&boundary/suf_reac_no')
        endif
        if(iq==suf_reac_no(jq)) then
          IDX_SUFRAC(nb,iq)=.true.
        endif
        enddo
        enddo
        if(lFC>0) then
          if(suf_reac_no(1)==1) then
            mem_cat(nb)=1
          elseif(suf_reac_no(1)==2) then
            mem_cat(nb)=2
          endif
        else
          mem_cat(nb)=0
        endif
!
        nbulk=0
        do nbk=1,mbulk
        if(bulk_name(nbk)/=' ') then
          nbulk=nbulk+1
          nphasx=nphasx+1
          phs_nam(nphasx)=bulk_name(nbk)
          phs_nbc(nphasx)=nb
          phs_typ(nphasx)=blkphase
          phs_thk(nphasx)=bulk_thick(nbk)
          nm=0
          do ncm=1,mcomp
          if(bulk_species_name(nbk,ncm)/=' ') then
            nm=nm+1
            nallcomp=nallcomp+1
            phs_snm(nallcomp)=bulk_species_name(nbk,ncm)
            if(bulk_species_dens(nbk,ncm)==undef)then
              call FFRABORT
     &           (1,'ERR: bulk_species_dens(:,:) not defined')
            else
              blk_dns(nallcomp)=bulk_species_dens(nbk,ncm)
            endif
            if(bulk_species_molefraction(nbk,ncm)==undef)then
              call FFRABORT
     &           (1,'ERR: bulk_species_molefraction(:,:) not defined')
            else
              phs_inifrac(nallcomp)=bulk_species_molefraction(nbk,ncm)
            endif
          endif
          enddo
          if(nm==0) then 
            write(ifle,'(a,I4,a)') 
     &      'ERR: NO site species name are defined at BC & BULK=',no,
     &      trim(bulk_species_name(nbk,ncm)) 
            call FFRABORT(1,'ERR: NO bulk species defined') 
          endif
          phs_idx(nphasx)=nm 
        endif
        enddo
      endif
!
      if(heat_release/=0.and.heat_release/=1) then 
        write(ifle,'(a)') 
     &  'ERR: Surface reaction heat release flag is error'
        write(ifle,'(a)') 
     &  'MSG: heat_release=0: NOT considering formation heat'
        write(ifle,'(a)') 
     &  'MSG: heat_release=1: DO considering formation heat'
        call FFRABORT(1,'ERR:')
      else
        heat(nb)=heat_release
      endif
      lcvdwal=.true.
      surfreac(nb)=1
!
      end subroutine surface_reac
!
      end subroutine inputdata
!
      end module module_boundary
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_cgsolver
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(LEN=21),parameter,private :: modnam='(module_cgsolver)'
!
!
!< module for cg solver >
!
      real(8)  :: aepscg=1.d-12
      real(8)  :: repscg=1.d-12
      real(8)  :: aepsbcg=1.d-12
      real(8)  :: repsbcg=1.d-12
      real(8)  :: aepscg_p=1.d-18,repscg_p=1.d-18,
     &                    aepsbcg_MHD=1.d-12,repsbcg_MHD=1.d-12
      integer :: iterbcg=1000,itercg_p=100,iterbcg_MHD=100
      integer :: itercg=1000
      integer :: nostop=0
!
! aepscg  : absolute error of tolerance for iccg
! repscg  : relative error of tolerance for iccg
! itercg  : maximul iteration counts for iccg
! aepsbcg : absolute error of tolerance for bicgstab
! repsbcg : relative error of tolerance for bicgstab
! iterbcg : maximul iteration counts for bicgstab
! nostop  :  =0 => iters of ICCG should stop at itercg=1000
!            =1 => iters of ICCG should NOT stop
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >--------------------------------------------
!==================================================================
      subroutine inputdata(ifli,ifll,ifle,my_rank,cntlnam,ierror)
!==================================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle,my_rank
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
      character(LEN=11),parameter :: subnam='(inputdata)'
!
! --- [namelist]
!
      namelist /cgsolver/ nostop,
     &                    aepscg,repscg,
     &                    itercg,
     &                    
     &                    aepsbcg,repsbcg,
     &                    iterbcg,
     &                    
     &                    aepscg_p,repscg_p,
     &                    itercg_p,
     &                    
     &                    aepsbcg_MHD,repsbcg_MHD,
     &                    iterbcg_MHD
     &  
!
! --- [local entities]
!
      integer :: iset=0
      integer :: ios,ierr1
!

      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      rewind ifli
      read(ifli,cgsolver,iostat=ios)
      if( ios.lt.0 ) return
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
      if(nostop/=0.and.nostop/=1) then
        write(ifle,*) '### ERR : data error'
        write(ifle,*) 'nostop =',nostop
        write(ifle,*) 
     &      'it must be  [nostop=0] or [nostop=1] in ',cntlnam
        write(ifll,*) 
     &   '### MSG: [nostop=0] : ICCG should stop at ',
     &               itercg
        write(ifll,'(a,I8,a)') 
     &'### MSG: [nostop=1] : ICCG should NOT stop while reach at',
     &               itercg,' and go into next time step'
        goto 9999
      endif
      if(my_rank==0.and.nostop==0) then
        write(ifll,*) 
     &   '### MSG: [nostop=0] : ICCG should stop at ',
     &               itercg
      endif
      if(my_rank==0.and.nostop==1) then
         write(ifll,'(a,I8,a)') 
     &'### MSG: [nostop=1] : ICCG should NOT stop while reach at ',
     &               itercg,' and go into next time step'
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
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_chemcntl
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
!
!< module for control data for chemical reaction >private
!
      character(LEN=15),private,parameter :: modnam='module_chemcntl'
      integer,save  :: igT_iter=0,
     &                 Zeldovich_P=0,
     &                 Fuel_NO_P=0,
     &                 prompt_NO_P=0
      integer  :: itintg=1000000
      integer  :: msizfc=10
      real(8)  :: divint=1.d-2
      real(8)  :: abserr=1.d-5
      real(8)  :: relerr=1.d-5
      
!
! itintg : maxmul no. of integration steps
! msizfc : factor of memory size
! divint : ratio of initial interval to total interval
! abserr : absolute error of tolerance
!        : for integration of ordinary differential equation
! relerr : relative error of tolerance
!        : for integration of ordinary differential equation
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      implicit none
! --- MXMATERIAL
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      integer,intent(out) :: ierror
      character(*),intent(in) :: cntlnam
!
      character(LEN=11),parameter :: subnam='(inputdata)'
!
      integer :: iset=0
      integer :: ios,ierr1
!
! --- [namelist]
!
      namelist /chemcntl/ igT_iter,
     &                    Zeldovich_P,Fuel_NO_P,prompt_NO_P,
     &                    itintg,msizfc,divint,abserr,relerr
!
! --- [local entities]
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      rewind ifli
      read(ifli,chemcntl,iostat=ios)
      if( ios.lt.0 ) return
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
!
      if( relerr.le.0.d0 ) then
        write(ifle,*) ' ### ERR : data error'
        write(ifle,*) 'relerr = ',relerr
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      
      if(igT_iter<0) then
        igT_iter=0
        write(ifle,*) ' ### WRN : &chemcntl/igT_iter is reset to 0'
      endif
!
      if(Zeldovich_P<0) then
        Zeldovich_P=0
        write(ifle,*) ' ### WRN : &chemcntl/Zeldovich_P is reset to 0'
      endif
!
      if(Fuel_NO_P<0) then
        Fuel_NO_P=0
        write(ifle,*) ' ### WRN : &chemcntl/Fuel_NO_P is reset to 0'
      endif
!
      if(prompt_NO_P<0) then
        prompt_NO_P=0
        write(ifle,*) ' ### WRN : &chemcntl/prompt_NO_P is reset to 0'
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
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_debug
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
!
!< module for debugging >
!
      integer,parameter,private :: mdebug=100
      character(LEN=14),parameter,private :: modnam= '(module_debug)'
      logical :: check_prep=.true.
!
      integer,private :: i
      integer :: idebug(mdebug)=(/(0,i=1,mdebug)/)

!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      implicit none
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
      character(LEN=11),parameter :: subnam='(inputdata)'
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
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_deltat
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
!< module for data to estimate time increment >
!
      integer,parameter :: const=1,auto=2
      character(LEN=15),parameter,private :: modnam='(module_deltat)' 
      integer :: ideltt=const
      real(8)  :: dtmax=1.d0, dtsafe=0.8d0
      real(8)  :: maxcou,maxcou2
      integer  :: icvmxc
!
! const  : constant time increment
! auto   : auto time increment
! ideltt : flag for time increment
! dtmax  : maximum delta-t
! dtsafe : safety factor for delta-t
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=8),parameter :: dtlst(2)=(/'constant','auto    '/)
      character(LEN=10) :: option
!
      real(8),parameter :: undef=-huge(1.d0)
      integer :: iset=0
      integer :: ios,ierr1
!
! --- [namelist]
!
      real(8)        :: dt,safe
      namelist /deltat/ option,dt,safe
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

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_dimnsn
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(LEN=15),parameter,private :: modnam='(module_dimnsn)'
!
!
!< module for dimension of computaional domain >
!
      logical :: xredc=.false., yredc=.false., zredc=.false.
!
! xredc : =.true.; x-direction is reduced
! yredc : =.true.; y-direction is reduced
! zredc : =.true.; z-direction is reduced
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
      character(LEN=11),parameter :: subnam='(inputdata)'
!
! --- [namelist]
!
      logical :: x,y,z
      namelist /dimension/ x,y,z
!
! --- [local entities]
!
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
      
      if(ios.lt.0) return
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

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_flags
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(LEN=18),parameter,private :: modnam='(module_flags)'
      integer,parameter :: face=0, center=1
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
!      real(8) :: rnck=0.3d0,bldfct=0.8D0
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
!-< 3. Option to integrate conv. & diff. equation >-
!
      integer,parameter :: eulere=0, euleri=1,Adams_Bashforth=2,
     &                     Adams_Moulton=3,Crank_Nicolson=4,
     &                     Runge_Kutta=5
      integer :: intgvv =euleri, intgty=euleri, intgke=euleri
      integer :: intgMHD=euleri

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

      integer :: icalrv=face
!
! face   : interpolate velocity at cell face
! center : interpolate velocity at cell center
! icalrv : flag to interpolate velocity
!
!-< 5. Poisson Matrix >-
      integer :: SYMM=1
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=======================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,lrans,ierror)
!=======================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      logical,intent(in)  :: lrans
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character( 8),parameter :: intlst(6)=(/'explicit','implicit',
     & 'adamsbfh','moultons','cranknic','rungekta'/)
      character( 6),parameter :: invlst(2)=(/'face  ','center'/)
      integer :: iset=0
      integer :: ios,ierr1,ierr2,ierr3,ierr4,ierr5,ierr6,ierr7
!
! --- [namelist]
!
      character(LEN=10) :: integ_vv,integ_ty,integ_ke,integ_MHD
      character(LEN=10) :: intrv
      namelist /flags/
     &                 integ_vv,integ_ty,integ_ke,
     &                 intrv,integ_MHD
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
      integ_MHD=' '
      intrv    =invlst(max(1,min(2,icalrv+1)))
      ierr1 = 0
      ierr2 = 0
      ierr3 = 0
      ierr4 = 0
      ierr5 = 0
      ierr6 = 0
      ierr7 = 0
!!
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
      call nml_listno(6,intlst,integ_MHD,intgMHD)
!      call nml_chkchr0(ifle,'cnv_vv'   ,cnv_vv   ,icnvvv,ierr1)
!      call nml_chkchr0(ifle,'cnv_ty'   ,cnv_ty   ,icnvty,ierr2)
!      call nml_chkchr0(ifle,'limitr_vv',limitr_vv,lmtrvv,ierr3)
!      call nml_chkchr0(ifle,'limitr_ty',limitr_ty,lmtrty,ierr4)
      call nml_chkchr0(ifle,'integ_vv' ,integ_vv ,intgvv,ierr5)
      call nml_chkchr0(ifle,'integ_ty' ,integ_ty ,intgty,ierr6)
      call nml_chkchr0(ifle,'integ_MHD',integ_MHD,intgty,ierr7)
      if( ierr1.ne.0 .or. ierr2.ne.0 .or.
     &    ierr3.ne.0 .or. ierr4.ne.0 .or.
     &    ierr5.ne.0 .or. ierr6.ne.0 .or.ierr7/=0) goto 9999
!
      if(lrans) then
        call nml_listno(6,intlst,integ_ke ,intgke)
        call nml_chkchr0(ifle,'integ_ke' ,integ_ke ,intgke,ierr3)
        if( ierr1.ne.0 .or. ierr2.ne.0 .or. ierr3.ne.0 ) goto 9999
      else
        intgke=1
      endif
!
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
      if(intgMHD.eq.Adams_Bashforth) then
        write(ifle,2000)
        write(ifle,*) '[integ_MHD] is switch to be [implicit] '
        intgMHD=euleri
      endif
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
      if(intgMHD.eq.Adams_Moulton) then
        write(ifle,2000)
        write(ifle,*) '[integ_MHD] is switch to be [implicit] '
        intgMHD=euleri
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
      if(intgMHD.eq.Crank_Nicolson) then
        write(ifle,2000)
        write(ifle,*) '[integ_MHD] is switch to be [implicit] '
        intgMHD=euleri
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
      if(intgMHD.eq.Runge_Kutta) then
        write(ifle,2000)
        write(ifle,*) '[integ_MHD] is switch to be [implicit] '
        intgMHD=euleri
      endif
      if(intgty.eq.Runge_Kutta) then
        write(ifle,2000)
        write(ifle,*) '[integ_ty] is switch to be [implicit] '
        intgty=euleri
      endif
      if(intgke.eq.Runge_Kutta) then
        write(ifle,2000)
        write(ifle,*) '[integ_ke] is switch to be [implicit] '
        if(lrans) then
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
!      icalrv=face
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
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_material
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(LEN=18),parameter,private :: modnam='(module_material)'
!
! --- < 1. parameter for fluid domain >
!
! nflud : no. of fluid domains
! lclsd : flag for closed domain
!       : =1; given pressure: BC (pressure BC, comp-outlet BC)
!       : =0; closed domain: no-outlet BC, incomp-outlet BC, zero-outlet
!
!
      integer,parameter :: suther=1, const=2,simpli=3,MK=4,usrprp=5,
     &                     chung=6
      real(8)  :: sthrmu_i=18.d-6, sthrt_i=293.15d0, sthrc_i=117.d0
      real(8)  :: prlmnr_i=0.72d0, sclmnr_i=1.d0
      integer,save :: rg_mark=0
!
! suther : Sutherland's formula
! const  : constant
! simpli : simplified transport model
! chung  : Chung's dense fluid transport model
! isthr  : flag for viscosity
! sthrmu : viscosity at temp.=sthrt [Pa.s]
! sthrt  : reference temperature [K]
! sthrc  : Sutherland's constant [K]
! prlmnr : Prandtl no. for laminar flow [-]
! sclmnr : Schmidt no. for laminar flow [-]
!
      integer,parameter :: up1st=1, up2nd=2,cnt2nd=3,up3rd=4,
     &                     usi2nd=5,engcsv=6,c3bld=7,
     &                     mscl=8
      integer,parameter :: nolim=0, venkat=1,slope=2
      real(8)  :: rnck=0.3d0,bldfct=0.8D0
      real(8) ,save,allocatable :: MUSCL_k(:)
      real(8),parameter :: r3=1.d0/3.d0,SML=1.d-10
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
! MUSCL_k : =1/3 :3rd order accurate scheme 
!           -1: 2nd-order fully upwind
!           0 : 2nd-order upwind biased
!           1 : 2nd central scheme
! ---  Option to calculate viscous term >-
!
      integer,parameter :: nocal=0, cal=1
!      integer :: ivscvv=cal, ivscty=cal, ivscke=cal
! nocal  : not calculate
!   cal  : calculate
! ivscvv : flag for viscous term of momentum equation
! ivscty : flag for diffusion term of temp. & conc. equation
! ivscke : flag for diffusion term of k-e equation
!
      integer,save      :: ical_porous=0
      integer,parameter :: ino=1,ipower=2,idarcy=3,idarcy2=4
!
! ipower : POWER-LAW
! idarcy : Darcy-law
! idarcy2: Darcy-law + (dF/dt)/E+(dF/dx)/E^2
!

!
      integer :: nflud=1
      integer ,save,allocatable :: lclsd(:),iclosd(:)
      integer ,save,allocatable :: nofld(:)
      integer ,save,allocatable :: isthr(:)
      integer ,save,allocatable :: PFstart(:)
      real(8) ,save,allocatable :: sthrmu(:),sthrt(:),sthrc(:)
      real(8) ,save,allocatable :: turbmu(:)
      real(8) ,save,allocatable :: prlmnr(:),sclmnr(:)
      real(8) ,save,allocatable :: r_prlmnr(:),r_sclmnr(:)
      real(8) ,save,allocatable :: sthrmu2(:),sthrt2(:),sthrc2(:)
      real(8) ,save,allocatable :: prlmnr2(:),sclmnr2(:)
      real(8) ,save,allocatable :: r_prlmnr2(:),r_sclmnr2(:)
      real(8) ,save,allocatable :: turbmu2(:)
      real*8 ,save,allocatable :: tsat(:),hgsat(:),hlsat(:)
      logical ,save,allocatable :: cvdist(:)
      integer ,save,allocatable :: icnvvv(:)
      integer ,save,allocatable :: icnvty(:)
      integer ,save,allocatable :: icnvke(:)
      integer ,save,allocatable :: icnvMHD(:) !???zhang
      integer ,save,allocatable :: lmtrvv(:)
      integer ,save,allocatable :: lmtrty(:)
      integer ,save,allocatable :: lmtrke(:)
      integer ,save,allocatable :: ivscvv(:)
      integer ,save,allocatable :: ivscty(:)
      integer ,save,allocatable :: ivscke(:) 
      integer ,save,allocatable :: ivscMHD(:) !???zhang
      real(8) ,save,allocatable :: relax_fp(:)
      real(8) ,save,allocatable :: relax_C(:)
      real(8) ,save,allocatable :: relax_C_step(:)
      real(8) ,save,allocatable :: relax_frs(:)
      real(8) ,save,allocatable :: relax_fh(:)
      real(8) ,save,allocatable :: relax_fy(:)
      real(8) ,save,allocatable :: relax_fv(:)
      real(8) ,save,allocatable :: relax_frc(:)
      real(8) ,save,allocatable :: bldf(:)
      real(8) ,save,allocatable :: rotati(:)
      real(8) ,save,allocatable :: domegadt(:)
      real(8) ,save,allocatable :: rot_ang(:)
      real(8) ,save,allocatable :: rot_angold(:)
      real(8) ,save,allocatable :: sav_ang(:)
      real(8) ,save,allocatable :: rot_init(:)
      real(8) ,save,allocatable :: rotN0(:)
      real(8) ,save,allocatable :: strtim(:)
      real(8) ,save,allocatable :: rotup(:)
      integer ,save,allocatable :: ishaft(:)
      real(8) ,save,allocatable :: end(:,:),begin(:,:)
      integer ,save,allocatable :: nsplit(:)
      real*8  ,save,allocatable :: dpmxv(:),incmp_prs(:)
      integer ,save,allocatable :: ivscsldMHD(:)
      integer ,save,allocatable :: icnvsldMHD(:)
      integer ,save,allocatable :: ivscfldMHD(:)
      integer ,save,allocatable :: icnvfldMHD(:)
      real(8) ,save,allocatable :: relaxfldMHD(:)
      real(8) ,save,allocatable :: relaxsldMHD(:)
      integer ,save,allocatable :: GMHDA_fld(:)
      integer ,save,allocatable :: GMHDJ0_fld(:)
      integer ,save,allocatable :: GMHDJe_fld(:)
      integer ,save,allocatable :: GMHDA_sld(:)
      integer ,save,allocatable :: GMHDJ0_sld(:)
      integer ,save,allocatable :: GMHDJe_sld(:)
      integer ,save             :: ical_sld=0
      integer ,save             :: N_pors=80
      integer ,save,allocatable :: porous(:)
      real(8) ,save,allocatable :: C0(:),C1(:),C2p(:)
      real(8) ,save,allocatable :: alpha_co(:),porosty(:)
      integer ,save,allocatable :: BGCOMP(:)
!
! --- < 2. parameter for solid domain >
!
! nsold  : no. of solids [-]
! nosld  : solid no. [-]
! rsld   : density [kg/m^3]
! cpsld  : specific heat [J/kg/K]
! rmdsld : thermal conductivity [W/m/K]
! rotati : rpm [revolution per second [rps]]
! ishaft : revolution axis (only: x,y,z ; can not shift and can not rotating)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! ical_sld : =0: no-sliding-BC.and.no-over-set-BC  & no-rpm; 
!          : =1: have-sliding-BC.or.have-over-set-BC & have-rpm
!            =2: have-sliding-BC.or.have-over-set-BC & no-rpm
!          : =3: no-sliding-BC.and.no-over-set-BC   & have-rpm
!          : =4: have-sliding-BC & have-rpm & 1-pitch
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer :: nsold=0
      integer,save,allocatable  :: nosld(:)
      integer,save,allocatable  :: PSstart(:)
      real(8) ,save,allocatable :: rsld(:)
      real(8) ,save,allocatable :: cpsld(:)
      real(8) ,save,allocatable :: rmdsld(:)
      real(8) ,save,allocatable :: relax_sh(:)  
!
!      integer :: KMAT_S=0
!      integer,save,allocatable  :: MAT_ALL(:)
!      integer,save,allocatable  :: MAT_S(:)
!     
      real(8) ,save,allocatable :: relaxp(:)
      real(8) ,save,allocatable :: relaxC(:)
      real(8) ,save,allocatable :: relaxCi(:)
      real(8) ,save,allocatable :: relaxC_step(:)
      real(8) ,save,allocatable :: relaxh(:)
      real(8) ,save,allocatable :: relaxys(:)
      real(8) ,save,allocatable :: relaxv(:)
      real(8) ,save,allocatable :: relaxrc(:)
      real(8) ,save,allocatable :: relaxrs(:)

      real(8) ,save,allocatable :: OMGA (:)
      real(8) ,save,allocatable :: relax_MHD(:)
      integer ,save,allocatable :: GAUGE_A(:)
      integer ,save,allocatable :: GAUGE_J0(:)
      integer ,save,allocatable :: GAUGE_Je(:)
!      logical,save,allocatable  :: MHD_E(:)
!      logical,save,allocatable  :: MHD_M(:)
!
!     relaxp,relaxh: under-relaxation-factor for steady flow
!
!
! --- <3. parameter for radiation>
      real*8,save,allocatable  :: radmat(:,:)
!1: gabsorb, 2: pabsorb, 3: pscatter, 4:scabeta, 5: ptcdiamter
      INTEGER,save,allocatable :: radfludflag(:,:)
!1: gastype, 2: scatype, 3: ptccal, 4: fgpflag


!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. store data >-----------------------------------------------------
!====================================================
      subroutine deficld
     &  (NMAT,my_rank,ifle,nfludx,MXMAT,MAT_NO,
     &   idrdp,incomp,mach0,comp,jclos,Pstart,
     &   ioval,ical_vof,steady,iHPC_close,MHD_steady,
     &   ierr)
!====================================================
! --- lclsd(IIMAT)=0 : outlet BC 
! --- lclsd(IIMAT)=1 : ipfix=1
! --- lclsd(IIMAT)=2 : zero mach
! --- lclsd(IIMAT)=3 : compressible 
      implicit none
      integer,intent(in)  :: NMAT,ifle,nfludx,MXMAT,my_rank
      integer,intent(in)  :: idrdp,incomp,mach0,comp
      integer,intent(in)  :: jclos(0:100),MAT_NO(0:MXMAT)
!      integer,intent(in)  :: MATC_SLD(100)
      integer,intent(out) :: ierr
      integer,intent(in)  :: ioval,ical_vof,MHD_steady
      logical,intent(in)  :: steady,iHPC_close
      integer,intent(out) :: Pstart(MXMAT)
!
      character(LEN=9),parameter :: subnam='(strdata)'
      integer :: i,j,iset=0,IIMAT,IMAT,IMAT_U,IJMAT,imat_uu
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierr=0
!
      allocate(lclsd(0:MXMAT),iclosd(0:MXMAT),stat=ierr)
!
      if(ierr.ne.0) then
        write(ifle,*) '### error : allocation failed'
        goto 9999
      endif
!
      lclsd(1:MXMAT)=0
      iclosd(1:MXMAT)=0

      iclosd(0)=1
      if(idrdp.eq.comp) then
        lclsd(0)=3
      elseif(idrdp.eq.mach0) then
        lclsd(0)=2
      elseif(idrdp.eq.incomp) then
        lclsd(0)=1
      endif

      do 100 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT>0) then
        if(idrdp.eq.comp) then
          lclsd(IIMAT)=3
        elseif(idrdp.eq.mach0) then
          if(jclos(IMAT)==1.or.jclos(IMAT)==2) then !open
            lclsd(IIMAT)=0       !zhang8(=2,=0)
          else
            lclsd(IIMAT)=2
          endif
        elseif(idrdp.eq.incomp) then
          if(jclos(IMAT)==1.or.jclos(IMAT)==2) then  !open
            lclsd(IIMAT)=0   !0
          else 
            lclsd(IIMAT)=1   !1
          endif 
        ENDIF

        IF(jclos(IMAT)==4) then !sliding
          lclsd(:)=0
          exit
        endif

        if(ical_vof.eq.1) then
          lclsd(IIMAT)=0
        endif

        if(ioval==1) then
!          lclsd(IIMAT)=0
           lclsd(IIMAT)=1  !fire_shinonome
        endif

        if(jclos(IMAT).eq.0) then !closed
          iclosd(IIMAT)=1
        else
          iclosd(IIMAT)=0
        endif
      endif
 100  enddo

!      do 100 i=1,nfludx
!      IIMAT=MATC_FLD(i)
!      IMAT=MAT_NO(IIMAT)
!      if(IMAT.lt.0) then
!        write(ifle,*) '### error : FLUID has Negative IMAT No.',IMAT
!        goto 9999
!      endif
!      if(idrdp.eq.comp) then
!        lclsd(IIMAT)=3
!!        if(jclos(IMAT).eq.1.or.jclos(IMAT).eq.2) then
!!          lclsd(IIMAT)=0
!!        else
!!          lclsd(IIMAT)=2
!!        endif
!      elseif(idrdp.eq.mach0) then
!        if(jclos(IMAT)==1.or.jclos(IMAT)==2) then 
!          lclsd(IIMAT)=0       !zhang8(=2,=0)
!          lclsd(0)=0
!        else
!          lclsd(IIMAT)=2
!          lclsd(0)=2
!        endif
!      elseif(idrdp.eq.incomp) then
!!        if(iHPC_close) then
!!          lclsd(IIMAT)=1
!!          lclsd(0)=1
!!        else
!!          lclsd(IIMAT)=0
!!          lclsd(0)=0
!!        endif

!        if(jclos(IMAT)==1.or.jclos(IMAT)==2) then 
!          lclsd(IIMAT)=0 
!          lclsd(0)=0 
!        else 
!          lclsd(IIMAT)=1 
!          lclsd(0)=1 
!        endif 
!!        if(jclos(IMAT)==1.or.jclos(IMAT)==2) then
!!          lclsd(IIMAT)=0
!!        else
!!          lclsd(IIMAT)=1
!!        endif
!        IF(jclos(IMAT)==4) then !sliding
!          lclsd(:)=0
!          exit
!        endif
!      endif
!
!      if(ical_vof.eq.1) then
!        lclsd(IIMAT)=0
!      endif
!
!      if(ioval==1) then
!        lclsd(IIMAT)=0
!      endif
!
! --- jclos=0: closed-domain; jclos=1: Outlet-domain
!----------------------------------------------------
!      if(jclos(IMAT).eq.0) then
!        iclosd(IIMAT)=1
!      endif
!  100 continue
!
! --- 
!
!      if(steady) then 
        do i=1,nflud 
          IMAT_U=nofld(i)
          DO IIMAT=1,NMAT
          IMAT_UU=MAT_NO(IIMAT)
          IJMAT=0
          if(IMAT_UU==IMAT_U) then
            IJMAT=IIMAT
            exit
          endif
          ENDDO
          if(IJMAT==0) exit
          if(steady) then 
            relaxp(IJMAT)=relax_fp(i)
            relaxC(IJMAT)=relax_C(i)
            relaxCi(IJMAT)=relax_C(i)
            relaxC_step(IJMAT)=relax_C_step(i)
            relaxrs(IJMAT)=relax_frs(i)
            relaxh(IJMAT)=relax_fh(i)
            relaxys(IJMAT)=relax_fy(i)
            relaxv(IJMAT)=relax_fv(i)
            relaxrc(IJMAT)=relax_frc(i)
          else
            relaxp(IJMAT)=1.d0
            relaxC(IJMAT)=relax_C(i)
            relaxCi(IJMAT)=relax_C(i)
            relaxC_step(IJMAT)=relax_C_step(i)
            relaxh(IJMAT)=1.d0
            relaxys(IJMAT)=1.d0
            relaxv(IJMAT)=1.d0
            relaxrc(IJMAT)=1.d0
            relaxrs(IJMAT)=1.d0
          endif
          Pstart(IIMAT)=PFstart(i)
        enddo

        do i=1,nsold
          IMAT_U=nosld(i)
          DO IIMAT=1,NMAT
          IMAT_UU=MAT_NO(IIMAT)
          IJMAT=0
          if(IMAT_UU==IMAT_U) then
            IJMAT=IIMAT
            exit
          endif
          ENDDO
          if(IJMAT==0) exit
          if(steady) then 
            relaxh(IJMAT)=relax_sh(i)
          else
            relaxh(IJMAT)=1.d0!relax_sh(i)
          endif
            Pstart(IIMAT)=PSstart(i)
        enddo
!
      do i=1,nflud
        IMAT_U=nofld(i)
        IJMAT=0
        DO IIMAT=1,NMAT
        IMAT_UU=MAT_NO(IIMAT)
        IJMAT=0
        if(IMAT_UU==IMAT_U) then
          IJMAT=IIMAT
          exit
        endif
        ENDDO
        if(IJMAT==0) exit
        icnvMHD(IJMAT)=icnvfldMHD(i) 
        ivscMHD(IJMAT)=ivscfldMHD(i)
        GAUGE_A(IJMAT)=GMHDA_fld(i)
        GAUGE_J0(IJMAT)=GMHDJ0_fld(i)
        GAUGE_Je(IJMAT)=GMHDJe_fld(i)
        if(MHD_steady==1) then
          relax_MHD(IJMAT)=relaxfldMHD(i)
        else
          relax_MHD(IJMAT)=1.d0 
        endif
      enddo
      do i=1,nsold
        IMAT_U=nosld(i)
        
        DO IIMAT=1,NMAT
        IJMAT=0
        IMAT_UU=MAT_NO(IIMAT)
        if(IMAT_UU==IMAT_U) then 
          IJMAT=IIMAT
          exit
        endif
        ENDDO
        if(IJMAT==0) exit
        icnvMHD(IJMAT)=icnvsldMHD(i)
        ivscMHD(IJMAT)=ivscsldMHD(i)
        GAUGE_A(IJMAT)=GMHDA_sld(i)
        GAUGE_J0(IJMAT)=GMHDJ0_sld(i)
        GAUGE_Je(IJMAT)=GMHDJe_sld(i)
        if(MHD_steady==1) then
          relax_MHD(IJMAT)=relaxsldMHD(i)
        else
          relax_MHD(IJMAT)=1.d0
        endif
      enddo
!
      return
 9999 continue
      write(ifle,*) modnam,subnam
      ierr=1
      end subroutine deficld
!
!==========================================================
      subroutine inputdata(ifli,ifll,ifle,my_rank,NOMAT,
     &  cntlnam,lrans,E2P,lvof,lcavi,lrot,imhd,iMvmsh,lFC,
     &  MIX_KNI,CHUN,
     &  ncompall,ivector,rg_markx,ierror)
!==========================================================
!
      implicit none
!s
! --- [dummy arguments]
!
      integer,intent(in)    :: ifli,ifll,ifle,my_rank,iMvmsh
      integer,intent(in)    :: ncompall,lFC,ivector
      character(*),intent(in) :: cntlnam
      logical,intent(in)    :: lrans,E2P,lvof,lcavi
      logical,intent(inout) :: lrot
      integer,intent(inout) :: NOMAT,imhd,MIX_KNI,CHUN
      integer,intent(out)   :: ierror,rg_markx
!
! --- [local entities] 
!
      character(LEN=3),parameter  ::
     &   cnvlst(8)=(/'1st','2nd','c2d','3rd','usi','eng','c3d','MSL'/)
      character(LEN=8),parameter  :: lmtlst(3)=
     &              (/'no      ','Venkatak','slope   '/)
      character(LEN=10),parameter :: mulst(6)=
     &         (/'Sutherland','constant  ','Simplified','MK        ',
     &           'user      '
     &          ,'Chung     '
     &  /) 

!      character(LEN=10),parameter :: mulst(6)=
!     &         (/'Sutherland','constant  ','Simplified','MK        ',
!     &           'user      '/) 
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=3),parameter  :: cond(2)=(/'yes','no '/)
      integer,parameter      :: iundef=-huge(1)
      real(8) ,parameter     :: undef=-huge(1.d0)
      integer                :: kwall=0
      integer                :: nm
      real(8),parameter :: SML=1.d-15
!
      character(LEN=5),parameter  ::
     &       poros_lst(4)=(/'no   ','power','Darcy','PEFC '/)
!
! --- [local entities]
!
      real(8) :: rdum
      integer :: i,iflw,itrb,ios,ierr1,ierr2,ierr3,ierr4
      integer :: iset=0
      integer :: n,nd
      logical :: nml_comp_eq,nml_comp_gt
!
! --- [namelist-fluid]
!
!     IMAT_U=1~80    : fluid
!     IMAT_U=81~100  : porous
!     IMAT_U=-100~-1 : solid
!
      integer       :: IMAT_U,Poisson_start
      character(LEN=10) :: muopt,wall_distance
      real(8)        :: Prandtl,Schmidt,mu,t0,c,mu_turb,mu_turb2
      logical        :: visc,diff_ty,diff_ke,diff_MHD
      real(8)        :: unrelax_P=1.d0,unrelax_T=1.d0,unrelax_Ys=1.d0,
     &                  unrelax_V=1.d0,unrelax_RANS=1.d0,
     &                  unrelax_RC=1.d0,
     &                  unrelax_MHD=1.d0,
     &                  unrelax_sound=1.d0
      integer        :: unrelax_sound_step=-1
      character(LEN=10) :: cnv_vv,cnv_ty,cnv_ke,cnv_MHD
      character(LEN=10) :: limitr_vv,limitr_ty,limitr_ke
      real(8)        :: Prandtl2,Schmidt2,mu2,t02,c2,ts,hgs,hls
      real(8)        :: rpm,end_x,end_y,end_z,begin_x,begin_y,begin_z
      real(8)        :: angle_init=0.d0 !(deg ')
      integer        :: start_rot_step,rotup_step
      integer        :: N_SPLIT
      character(LEN=5) :: porous_law
      real(8)        :: C0_porous,C1_porous,perm_porous,C2_porous,
     &                  porosity,MUSCL_kapa
      integer        :: bg_gas=1
      real(8)        :: gabsorb=0.d0,pabsorb=0.d0,pscatter=0.d0,
     &                  scabeta=0.d0,ptcdia=0.d0,ptcpara(2)
      INTEGER        :: scatype=1,gastype=1,ptccal=0,fgpflag=0

      INTEGER        :: Gauge_MHD_A,Gauge_MHD_Je,Gauge_MHD_J0
      INTEGER        :: IMAT_DONOR=0
!     
!------------------------------------------------------------------------
!scatype
!         1,  Linear-anisotropic scatter / isotropic scatter, OK
!	  2,  large-diffuse particle scatter, OK
!	  3,  Rayleigh scatter, OK
!	  4,  Mie scatter
!	  5,  user defined scatter
!-----------------------------------------------------------------------
!porous_law
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
     &
     &                 gabsorb,pabsorb,pscatter,scabeta,ptcdia,
     &                 gastype,scatype,ptccal,fgpflag,ptcpara,
     &
     &                 unrelax_MHD,unrelax_sound,unrelax_sound_step,
     &                 Gauge_MHD_A,Gauge_MHD_Je,Gauge_MHD_J0,
     &                 IMAT_DONOR
     &                 
!------------------------
! --- [namelist-solid]
!------------------------
      real(8)  :: rho,cp,cndct
      namelist /solid/ IMAT_U,rho,cp,cndct,Poisson_start,unrelax_T,
     &                 diff_MHD,cnv_MHD,unrelax_MHD,
     &                 Gauge_MHD_A,Gauge_MHD_Je,Gauge_MHD_J0
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!----------------------------------------------
! --- fluid material
!----------------------------------------------
!--< 1.1 count up sizes >--
!
      nd=0
      rewind ifli
!
      do
        read(ifli,fluid,iostat=ios)
        if(ios.lt.0) exit
        call nml_errmsg0(ifle,ios,'fluid',ierr1)
        nd=nd+1
        if( ierr1.ne.0 ) then
          write(ifle,*) 'sequence no. of the namelist =',nd
          goto 9999
        endif
      enddo
  101 continue
      if( nd.lt.1 ) then
        write(ifle,*) '### warning: no fluid specified'
!MHD_err        call FFRABORT(1,'module_material/inputdata')
      endif
      nflud=nd
!
!--< 1.2 allocate arrays >--
!
      allocate( nofld (nflud),
     &          isthr (nflud),
     &          sthrmu(nflud),
     &          turbmu(nflud),
     &          sthrt (nflud),
     &          sthrc (nflud),
     &          prlmnr(nflud),
     &          sclmnr(nflud),
     &          r_prlmnr(nflud),
     &          r_sclmnr(nflud),
     &          sthrmu2(nflud),
     &          sthrt2 (nflud),
     &          sthrc2 (nflud),
     &          prlmnr2(nflud),
     &          sclmnr2(nflud),
     &          r_prlmnr2(nflud),
     &          r_sclmnr2(nflud),
     &          turbmu2(nflud),
     &          tsat(nflud),
     &          hgsat (nflud),
     &          hlsat (nflud),
     &          MUSCL_k(nflud),
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
     &          C0(nflud),C1(nflud),C2p(nflud),alpha_co(nflud),
     &          porosty(nflud),
     &          BGCOMP(nflud),
     &          porous(nflud),
     &          relax_fp(nflud),
     &          relax_C(nflud),
     &          relax_C_step(nflud),
     &          relax_fv(nflud),
     &          relax_frc(nflud),
     &          relax_fh(nflud),
     &          relax_fy(nflud),
     &          relax_frs(nflud),
     &          domegadt(0:nflud),
     &          rotati(0:nflud),
     &          rot_ang(0:nflud),
     &          rot_angold(0:nflud),
     &          sav_ang(0:nflud),
     &          rot_init(0:nflud),
     &          rotN0(0:nflud),
     &          strtim(0:nflud),
     &          rotup(0:nflud),
     &          ishaft(0:nflud),
     &          nsplit(0:nflud),
     &          end(3,0:nflud),begin(3,0:nflud),
     &          ivscfldMHD(nflud),
     &          icnvfldMHD(nflud),
     &          relaxfldMHD(nflud),

     &          GMHDA_fld(nflud),
     &          GMHDJe_fld(nflud),
     &          GMHDJ0_fld(nflud),

     &		radmat(5,nflud),radfludflag(4,nflud),
     &          stat=ierr1 )
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error-1 : allocation failed'
        call FFRABORT(1,'')
      endif
      nofld =0
!
      sthrmu(:)=sthrmu_i
      turbmu(:)=0.d0
      turbmu2(:)=0.d0
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
!      icnvMHD(:)=up1st
      icnvfldMHD(:)=up1st
      relaxfldMHD(:)=1.d0
      
!
      GMHDA_fld(:)=0
      GMHDJe_fld(:)=0
      GMHDJ0_fld(:)=0
!
      icnvke(:)=up1st
      lmtrvv(:)=nolim
      lmtrty(:)=nolim
      lmtrke(:)=nolim
      bldf(:)=0.8D0
      ivscvv(:)=cal
      ivscty(:)=cal
      ivscke(:)=cal
      ivscfldMHD(:)=cal
!      ivscMHD(:)=cal
      relax_fp(:)=1.d0
      MUSCL_k(:)=1.d0/3.d0
      relax_C(:)=1.d0
      relax_C_step(:)=1.d0
      relax_frs(:)=1.d0
      relax_fv(:)=1.d0
      relax_frc(:)=1.d0
      relax_fh(:)=1.d0
      relax_fy(:)=1.d0
!      OMGA_F(:)=0.D0
      rotati(:)=0.D0
      domegadt(:)=0.d0
      rot_ang(:)=0.d0
      rot_angold(:)=0.d0
      sav_ang(:)=0.d0
      rot_init(:)=0.d0
      rotN0(:)=0.D0
      strtim(:)=0.D0
      rotup(:)=0.d0
      ishaft(:)=0
      end(:,:)=0.D0
      end(3,:)=1.D0
      begin(:,:)=0.D0
      C0(:)=0.d0
      C1(:)=0.d0
      C2p(:)=0.d0
      alpha_co(:)=1.d0
      porous(:)=1
      porosty(:)=1.d0
      BGCOMP(:)=1
!
!-< 1.3 Input namelist >-
!
      nd=0
      rewind ifli
 1000 continue
!
      IMAT_U =iundef
      muopt  =' '
      mu     =undef
      mu_turb=0.d0
      mu_turb2=0.d0
      t0     =undef
      c      =undef
      Prandtl=prlmnr_i
      Schmidt=sclmnr_i
      mu2     =undef
      t02     =undef
      c2      =undef
      Prandtl2=prlmnr_i
      Schmidt2=sclmnr_i
      rpm=undef
      start_rot_step=0
      rotup_step=0
      angle_init=0.d0
      end_x=0.d0
      end_y=0.d0
      end_z=1.d0
      begin_x=0.d0
      begin_y=0.d0
      begin_z=0.d0
      N_SPLIT=36
      bg_gas=1
!
      wall_distance='no '
      Poisson_start=iundef
      bldfct=0.8d0
      MUSCL_kapa=1.d0/3.d0
      visc     =.true.
      diff_ty  =.true.
      diff_ke  =.true.
      diff_MHD =.true.
      cnv_vv   =' '
      cnv_ty   =' '
      cnv_ke   =' '
      cnv_MHD  =' '
      limitr_vv=lmtlst(slope+1)
      limitr_ty=lmtlst(slope+1)
      limitr_ke=lmtlst(slope+1)
!
      C0_porous=0.d0
      C1_porous=0.d0
      porous_law='no    '
      C2_porous=undef
      perm_porous=undef
      porosity=undef
!
!      OMEGA_MHD=0.d0
!
      Gauge_MHD_A=0
      Gauge_MHD_Je=0
      Gauge_MHD_J0=0
!
      read(ifli,fluid,iostat=ios)
      if( ios.lt.0 ) goto 1001
      call nml_errmsg0(ifle,ios,'fluid',ierr1)
      if( ierr1.ne.0 ) call FFRABORT(1,'')
      nd=nd+1
!
      if(IMAT_U.eq.iundef) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : IMAT_U'
        call FFRABORT(1,'')
      elseif( IMAT_U.le.0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'it must be > 0'
        call FFRABORT(1,'')
      endif
      do 1100 n=1,nd-1
      if( IMAT_U.eq.nofld(n) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'this IMAT_U already has been specified'
        write(ifle,*) 'sequence no. of the line =',n
        call FFRABORT(1,'')
      endif
 1100 continue
!
      if(porous_law/='no     ') then
        if(IMAT_U<=N_pors) then
          write(ifle,'(1X,a)') 
     &     'ERR: Porous madia: must bc [IMAT_U>N_pors]'
          call FFRABORT(1,'Porous madia model:')
        endif
        call nml_listno(4,poros_lst,porous_law,porous(nd))
        call nml_chkchr0(ifle,'porous_law',porous_law,porous(nd),ierr1)
        if( ierr1.ne.0 ) call FFRABORT(1,'ERR: Porous madia model:')
        ical_porous=1
        if(porous(nd)==ino) then
          
        elseif(porous(nd)==ipower) then
          C0(nd)=C0_porous
          C1(nd)=C1_porous
        elseif(porous(nd)==idarcy) then
!          if (C2_porous==undef.or.perm_porous==undef) then
!            write(ifle,'(1X,2a)') 
!     & 'MSG: Define Forchheimer coeff. [C2_porous] ',
!     & 'and permeability [perm_porous]'
!            call FFRABORT(1,'C2_porous,perm_porous')
!          else
!            C2p(nd)=C2_porous
!            alpha_co(nd)=perm_porous
!          endif
!
          if(porosity/=undef) then
            C2p(nd)=1.75d0/sqrt(150.d0*porosity**3)
            alpha_co(nd)=porosity**3*180.d0
            write(ifle,'(1X,2a)') 
     &      'MSG: C2=1.75d0/sqrt(150.d0*porosity**3)'
            write(ifle,'(1X,2a)') 'MSG: Calpha=porosity**3*180.d0'
          else
            if (C2_porous==undef.or.perm_porous==undef) then
              write(ifle,'(1X,2a)') 
     &        'MSG: Define Forchheimer coeff. [C2_porous] ',
     &        'and permeability [perm_porous]'
              call FFRABORT(1,'C2_porous,perm_porous')
            else
              C2p(nd)=C2_porous
              alpha_co(nd)=perm_porous
            endif
          endif
        elseif(porous(nd)==idarcy2) then
          if (perm_porous==undef) then
            write(ifle,'(1X,a)') 
     & 'ERR: permeability [perm_porous]'
            call FFRABORT(1,'MSG: set: perm_porous') 
          else 
            C2p(nd)=0.d0
            alpha_co(nd)=perm_porous
          endif

          if(porosity==undef) then 
            call FFRABORT(1,'MSG: porosity is NOT defined')
          else 
            porosty(nd)=porosity !1.75d0/sqrt(150.d0*porosity**3)
            alpha_co(nd)=perm_porous !porosity**3*180.d0
            if(porosity<1.d-25) call FFRABORT(1,'ERR:porosity<1.d-25')
          endif
        endif
      endif
!
      call nml_listno(6,mulst,muopt,isthr(nd))
      call nml_chkchr0(ifle,'muopt',muopt,isthr(nd),ierr1)
      if( ierr1.ne.0 ) call FFRABORT(1,'')
!
      select case(isthr(nd))                ! transport model
        case(1)
          call read_Sutherland
        case(2)
          call read_constant
        case(3)
          call read_Simplified
        case(4)
          call read_MK
        case(usrprp)
          call read_constant
        case(6)
          call read_Chung          
      end select
!
      if(E2P.or.lvof.or.lcavi) then
        select case(isthr(nd))
          case(1)
             call read_Sutherland2
          case(2)
            call read_constant2
          case(3)
            call read_Simplified2
        end select
      endif
!
      if( Prandtl.le.0.d0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'Prandtl =',Prandtl
        write(ifle,*) 'it must be > 0'
        call FFRABORT(1,'')
      endif
!
      if(E2P.or.lvof.or.lcavi) then
        if( Prandtl2.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'Prandtl2 =',Prandtl2
          write(ifle,*) 'it must be > 0'
          call FFRABORT(1,'')
        endif
      endif
!
      if(ncompall>1) then
        if(bg_gas>ncompall) then
          write(ifle,'(1x,a)') 
     & 'ERR: defined bg_gas as one of comps in &species'
          call FFRABORT(1,'ERR: Set [&fluid/bg_gas<ncomp]')
        endif
        if(bg_gas==0) then
          write(ifle,'(1x,a)') 
     & 'ERR: defined bg_gas as one of comps in &species'
          call FFRABORT(1,'ERR: Set [&fluid/bg_gas>0]')
        endif
        BGCOMP(nd)=bg_gas
      endif
!
      if(bldfct.gt.1.d0.or.bldfct.lt.0.d0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'Blend Factor (bldfct) must between 1 and 0'
        call FFRABORT(1,'')
      endif
!
!      if(rpm/=0.d0) then
      if(rpm/=undef) then
        if(Poisson_start/=0) then
          write(ifle,*) '### error : data error'
          write(ifle,*) '### error : Poisson_start must be 0 rpm/=0'
          call FFRABORT(1,'set Poisson_start=0 in fflow.ctl')
        endif
        rdum=(end_x-begin_x)*(end_x-begin_x)
     &      +(end_y-begin_y)*(end_y-begin_y)
     &      +(end_z-begin_z)*(end_z-begin_z)
        if(rdum.lt.SML) then
          write(ifle,*) 'ERR: You must define revolution shaft '
          write(ifll,*) 
     &         'MSG: Vector of revolution shaft is defined as ',
     &                  'with : end_x,end_y,end_z',
     &                  'and  : begin_x,begin_y,begin_z'
          call FFRABORT(1,'ERR: NOT define revolution shaft')
        endif
        if(N_SPLIT.lt.0)then
          write(ifle,*) 'ERR: N_SPLIT must be > 0'
          call FFRABORT(1,'ERR: NOT define [N_SPLIT]')
        endif
        if(abs(rpm).lt.SML) then
          ishaft(nd)=0
          rotati(nd)=0.d0
          rot_ang(nd)=angle_init*3.1415926d0/180.d0
          rot_angold(nd)=angle_init*3.1415926d0/180.d0 
          sav_ang(nd)=angle_init*3.1415926d0/180.d0    
          rot_init(nd)=angle_init*3.1415926d0/180.d0   
          rotN0(nd)=0.d0
          nsplit(nd)=0
          end(1,nd)=end_x
          end(2,nd)=end_y
          end(3,nd)=end_z
          begin(1,nd)=begin_x
          begin(2,nd)=begin_y
          begin(3,nd)=begin_z
          strtim(nd)=0.d0
          rotup(nd)=0.d0
        else
          ishaft(nd)=1
          rotati(nd)=0.d0
          rot_ang(nd)=angle_init*3.1415926d0/180.d0
          rot_angold(nd)=angle_init*3.1415926d0/180.d0
          sav_ang(nd)=angle_init*3.1415926d0/180.d0
          rot_init(nd)=angle_init*3.1415926d0/180.d0
          rotN0(nd)=rpm/60.d0*2.d0*3.1415926d0
          end(1,nd)=end_x
          end(2,nd)=end_y
          end(3,nd)=end_z
          begin(1,nd)=begin_x
          begin(2,nd)=begin_y
          begin(3,nd)=begin_z
          nsplit(nd)=N_SPLIT
          strtim(nd)=dble(max(start_rot_step,0))
          rotup(nd)=dble(max(rotup_step,1))
        endif
      endif
!
      if(unrelax_sound.gt.1.d0.or.unrelax_sound.lt.0.d0) then
        write(ifle,*) ' ERR: data error, unrelax_sound= ',unrelax_sound
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        goto 9999
      endif
      if(abs(unrelax_sound-1.d0)>1.d-10) then
        if(unrelax_sound_step==-1) then
          write(ifle,'(1x,2a)') 
     &   'ERR: set [unrelax_sound_step] ',
     &   'for under-relax number of sound speed'
          call FFRABORT(1,'ERR: in &fluid ')
        elseif(unrelax_sound_step<0) then
          write(ifle,'(1x,2a)') 
     &   'ERR: set [unrelax_sound_step>=0] ',
     &   'for under-relax number of sound speed in unrelax time step'
          call FFRABORT(1,'ERR: in &fluid ')
        endif
        unrelax_sound_step=max(1,unrelax_sound_step)
      endif
!
      if(unrelax_P.gt.1.d0.or.unrelax_P.lt.0.d0) then
        write(ifle,*) '### error : data error, unrelax_P= ',unrelax_P
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        goto 9999
      endif
! --- 
      if(unrelax_MHD.gt.1.d0.or.unrelax_MHD.lt.0.d0) then
        write(ifle,*) 
     &      '### error : data error, unrelax_MHD= ',unrelax_MHD
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        goto 9999
      endif
!
      if(unrelax_Ys.gt.1.d0.or.unrelax_Ys.lt.0.d0) then
        write(ifle,*) '### error : data error, unrelax_s= ',unrelax_Ys
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        goto 9999
      endif
!
      if(unrelax_T.gt.1.d0.or.unrelax_T.lt.0.d0) then
        write(ifle,*) '### error : data error, unrelax_T=',unrelax_T
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        goto 9999
      endif
!
      if(unrelax_V.gt.1.d0.or.unrelax_V.lt.0.d0) then
        write(ifle,*) '### error : data error, unrelax_V=',unrelax_V
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        goto 9999
      endif
!
      if(unrelax_RC.gt.1.d0.or.unrelax_RC.lt.0.d0) then
        write(ifle,*) '### error : data error, unrelax_RC=',unrelax_RC
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        write(ifle,*) 'Under-relaxation faceor Rhie-Chow ==> 0'
        goto 9999
      endif
!
      call nml_listno(2,cond,wall_distance,kwall)
      call nml_chkchr0(ifle,'wall_distance',wall_distance,kwall,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if(Poisson_start.eq.iundef) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : Poisson_start in namelist/fluid'
        write(ifle,*) 'Poisson_start is starting iter for Poisson'
        goto 9999
      elseif(Poisson_start.le.-1) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'Poisson_start =',Poisson_start
        write(ifle,*) 'it must be >= 0'
        goto 9999
      endif
!
      call nml_listno(8,cnvlst,cnv_vv   ,icnvvv(nd))
      if(ivector==1) then
        if(.not.(icnvvv(nd)==up1st
     &       .or.icnvvv(nd)==cnt2nd
     &       .or.icnvvv(nd)==up3rd)) then
          write(ifle,*) 'VECTOR VER. NOLY SUPPORT 1st, c2d,3rd for vel'
          call FFRABORT(1,'ERR:VECTOR VER. NOLY SUPPORT 1st, c2d,3rd')
        endif
      endif
      call nml_listno(8,cnvlst,cnv_ty   ,icnvty(nd))
      if(ivector==1) then
        if(.not.(icnvty(nd)==up1st
     &       .or.icnvty(nd)==cnt2nd
     &       .or.icnvty(nd)==up3rd)) then
          write(ifle,*) 
     &    'VECTOR VER. NOLY SUPPORT 1st, c2d,3rd for tmp & Ys'
          call FFRABORT(1,'ERR:VECTOR VER. NOLY SUPPORT 1st, c2d,3rd')
        endif
      endif
      if(icnvty(nd)==engcsv)  
     &    call FFRABORT(1,'ERR:Scheme [eng] only for velocity')
      call nml_listno(3,lmtlst,limitr_vv,lmtrvv(nd))
      call nml_listno(3,lmtlst,limitr_ty,lmtrty(nd))
      call nml_chkchr0(ifle,'cnv_vv'   ,cnv_vv   ,icnvvv(nd),ierr1)
      call nml_chkchr0(ifle,'cnv_ty'   ,cnv_ty   ,icnvty(nd),ierr2)

      call nml_chkchr0(ifle,'limitr_vv',limitr_vv,lmtrvv(nd),ierr3)
      call nml_chkchr0(ifle,'limitr_ty',limitr_ty,lmtrty(nd),ierr4)
      if(imhd>0) then
        call nml_listno(8,cnvlst,cnv_MHD  ,icnvfldMHD(nd))
        if(icnvfldMHD(nd)==engcsv) 
     &    call FFRABORT(1,'ERR:Scheme [eng] only for velocity')
        call nml_chkchr0(ifle,'cnv_MHD',cnv_MHD  ,icnvfldMHD(nd),ierr2)
        relaxfldMHD(nd)=unrelax_MHD
        GMHDA_fld(nd)=Gauge_MHD_A
        GMHDJ0_fld(nd)=Gauge_MHD_J0
        GMHDJe_fld(nd)=Gauge_MHD_Je
      endif
!      
      if( ierr1.ne.0 .or. ierr2.ne.0 .or.
     &    ierr3.ne.0 .or. ierr4.ne.0 ) goto 9999
!
      if( lrans ) then
        call nml_listno(8,cnvlst,cnv_ke   ,icnvke(nd))
        if(ivector==1) then
          if(.not.(icnvke(nd)==up1st
     &         .or.icnvke(nd)==cnt2nd
     &         .or.icnvke(nd)==up3rd)) then
            write(ifle,*) 
     &     'VECTOR VER. NOLY SUPPORT 1st, c2d,3rd for K,e & scalars'
          call FFRABORT(1,'ERR:VECTOR VER. NOLY SUPPORT 1st, c2d,3rd')
          endif
        endif
        if(icnvke(nd)==engcsv)  
     &    call FFRABORT(1,'ERR:Scheme [eng] only for velocity')
        call nml_listno(3,lmtlst,limitr_ke,lmtrke(nd))
        call nml_chkchr0(ifle,'cnv_ke'   ,cnv_ke   ,icnvke(nd),ierr1)
        call nml_chkchr0(ifle,'limitr_ke',limitr_ke,lmtrke(nd),ierr2)
        if(lvof.or.lcavi) icnvke(nd)=1
        if( ierr1.ne.0 .or. ierr2.ne.0 .or. ierr3.ne.0 ) goto 9999
      else
        icnvke(nd)=1
        lmtrke(nd)=lmtrke(nd)+1
      endif
      ivscvv(nd)=ixl(visc)
      ivscty(nd)=ixl(diff_ty)
      ivscke(nd)=ixl(diff_ke)
!
      MUSCL_k(nd)=MUSCL_kapa
      if(icnvvv(nd)==up3rd) then
        MUSCL_k(nd)=1.d0/3.d0
        write(ifll,'(1x,a)') 
     &   'MSG: MUSCL_kapa=1.0/3.0 MUST be set for 3rd upwind scheme '
        write(ifll,'(1x,a)') 
     &   'MSG: MUSCL_kapa=1.0/3.0 be reset in FFR'
      elseif(icnvvv(nd)==mscl) then
        if(ABS(MUSCL_kapa+1.d0)<SML) then
          write(ifll,'(1x,a)') 'MUSCL scheme is used'
          write(ifll,'(1x,a)') 'MUSCL_kapa=-1: 2nd-order fully upwind'
        elseif(ABS(MUSCL_kapa)<SML) then
          write(ifll,'(1x,a)') 'MUSCL scheme is used'
          write(ifll,'(1x,a)') 'MUSCL_kapa=0: 2nd-order upwind biased'
        elseif(ABS(MUSCL_kapa-1.d0)<SML) then
          write(ifll,'(1x,a)') 'MUSCL scheme is used'
          write(ifll,'(1x,a)') 'MUSCL_kapa=1: 2nd central scheme'
        else
          write(ifll,'(1x,a)') 'MUSCL scheme is used'
          write(ifll,*) 
     &         'User defined accurate scheme: MUSCL_kapa=',MUSCL_kapa
        endif
      endif

      ivscfldMHD(nd)=ixl(diff_MHD)
!
      lmtrvv(nd)=lmtrvv(nd)-1
      lmtrty(nd)=lmtrty(nd)-1
      lmtrke(nd)=lmtrke(nd)-1
      if(.not.lrans) ivscke(nd)=nocal
!
      nofld(nd)=IMAT_U
      if(isthr(nd)==1.or.isthr(nd)==2) then
        if(mu==undef) then 
          call FFRABORT(1,'ERR: [mu] NOT set')
        else
          sthrmu(nd)=mu
        endif 
      endif
      turbmu(nd)=mu_turb
      sthrt(nd) =t0 
      sthrc(nd) =c
      if( isthr(nd)/=3 ) prlmnr(nd)=Prandtl 
      prlmnr(nd)=Prandtl
      sclmnr(nd)=Schmidt
      if(E2P.or.lvof.or.lcavi.or.lFC>0) then 
        sthrmu2(nd)=mu2
        sthrt2(nd) =t02
        sthrc2(nd) =c2
        prlmnr2(nd)=Prandtl2
        sclmnr2(nd)=Schmidt2
        turbmu2(nd)=mu_turb2
        if(.not.lcavi) then
          tsat(nd) =ts
          hgsat(nd)=hgs
          hlsat(nd)=hls
        endif
      endif
      bldf(nd)=bldfct
      relax_fp(nd)=unrelax_P
      relax_C(nd)=unrelax_sound
      relax_C_step(nd)=unrelax_sound_step
      relax_fh(nd)=unrelax_T
      relax_fy(nd)=unrelax_Ys
      relax_fv(nd)=unrelax_V
      relax_frc(nd)=unrelax_RC
      relax_frs(nd)=unrelax_RANS
      if(kwall.eq.1) then
        cvdist(nd)=.true.
      elseif(kwall.eq.2) then
        cvdist(nd)=.false.
      else
        write(ifle,*) "ERR: wall_distance='yes' or 'no' in ",cntlnam
        call FFRABORT(1,'module_material/inputdata')
      endif
      PFstart(nd)=Poisson_start
!
      if(icnvty(nd)/=up1st) then
        write(ifle,'(1x,A)') 
     &  "WRN: cnv_ty MUST be [1st] for Unboundness problem"
        call FFRABORT(1,"MSG: Reset [cnv_ty= '1st']")
      endif
      if((icnvty(nd)==up3rd.or.icnvty(nd)==up2nd).and.
!     &   lmtrty(nd)==nolim) then
     &   lmtrty(nd)==slope) then
        write(ifle,'(1x,2A)') 
     &  "WRN: [limitr_ty= 'no'] ",
     &  "MUST be used with [3rd] & [2nd] for t/ys equation."
        call FFRABORT(1,"MSG: Reset [limitr_ty= 'slope']")
      endif
      if((icnvfldMHD(nd)==up3rd.or.icnvfldMHD(nd)==up2nd))then
        write(ifle,'(1x,2A)') 
     &  "WRN: MHD ",
     &  "MUST be  [1st] & [c2d] for MHD equation."
        call FFRABORT(1,"MSG: Reset MHD ")
      endif
!----------------------------
! --- radiation
!----------------------------
      if(gabsorb.lt.0.d0) then
        write(ifle,'(1x,a)') 'ERR: gas absorption coefficient (m^-1)'
        write(ifle,'(1x,a)') 'gabsorb =',gabsorb
        write(ifle,'(1x,a)') 'it should be >= 0.0'
        call FFRABORT(1,'ERR: &fluid/gabsorb')
      end if
      if(pabsorb.lt.0.d0) then
        write(ifle,'(1x,a)') 
     &     'ERR: particle absorption coefficient (m^-1)'
        write(ifle,'(1x,a)') 'pabsorb =',pabsorb
        write(ifle,'(1x,a)') 'it should be >= 0'
        call FFRABORT(1,'ERR: &fluid/pabsorb')
      end if
      if(pscatter.lt.0.d0) then
        write(ifle,'(1x,a)') 
     &  'ERR: particle scattering coefficient (m^-1)'
        write(ifle,'(1x,a)') 'pscatter =',pscatter
        write(ifle,'(1x,a)') 'it should be >= 0.d0'
        call FFRABORT(1,'ERR: &fluid/pscatter')
      end if
      if(scatype/=1.and.scatype/=2.and.scatype/=3.and.
     &   scatype/=4.and.scatype/=5) then
        write(ifle,'(1x,a)') 'ERR: defining particle scattering model '
        write(ifle,'(1x,a)') 'scatype =',scatype
        write(ifle,'(1x,a)') 'it should be in 1,2,3,4,5'
        write(ifle,'(1x,a)') 
     &   'scatype=1: Linear-anisotropic scatter / isotropic scatter'
        write(ifle,'(1x,a)') 
     &   'scatype=2: large-diffuse particle scatter'
        write(ifle,'(1x,a)') 
     &   'scatype=3: Rayleigh scatte'
        write(ifle,'(1x,a)') 
     &   'scatype=4: Mie scatter'
        write(ifle,'(1x,a)') 
     &   'scatype=5: user defined scatter'
        call FFRABORT(1,'ERR: &fluid/scatype')
      elseif(scatype==4) then
        write(ifle,'(1x,a)') 
     &   'scatype=4: Mie scatter'
        call FFRABORT(1,'ERR: &fluid/scatype=4: NOT finished')
      elseif(scatype==5) then
        write(ifle,'(1x,a)') 
     &   'scatype=5: user defined scatter'
        call FFRABORT(1,'ERR: &fluid/scatype=5: NOT finished')
      end if
!-----------------------------
! --- gastype=1 : gray gas
! --- gastype=2 : real gas
!-----------------------------
      if(gastype/=1.and.gastype/=2) then
        write(ifle,'(1x,a)') 'ERR: defining gas type'
        write(ifle,'(1x,a)') 'MSG: gastype=1: gray gas'
        write(ifle,'(1x,a)') 'MSG: gastype=1: real gas'
        write(ifle,'(1x,a)') 'gastype =',gastype
        write(ifle,'(1x,a)') 'it should be, 1: gray_gas, 2: real_gas'
        call FFRABORT(1,'ERR: &fluid/gastype')
      end if
!
      if(scatype==1.and.pscatter.GT.0.0) then
        if(abs(scabeta).gt.1d0) then
          write(ifle,'(1x,a)') 'ERR: coefficient for ',
     &      'Linear-anisotropic scatter / isotropic scatter'
	  write(ifle,'(1x,a)') 'scabeta =',scabeta
	  write(ifle,'(1x,a)') 'it should be in (-1 ~ 1)'
          call FFRABORT(1,'ERR: &fluid/scabeta')
	endif
      endif
!
      if(fgpflag/=0) then
        call FFRABORT(1,'ERR: NOT supported')
      endif
!
      if(ptccal.EQ.1.and.ptcdia.LE.0.d0) then
	write(ifle,'(1x,a)') 'ERR: particle averaging dia.: ptccal'
	write(ifle,'(1x,a)') 'ptcdia =',ptcdia
	write(ifle,'(1x,a)') 'particle diamater should be > 0'
	goto 9999
      end if
!
      if(gastype.ne.1.and.gastype.ne.2) then
        write(ifle,'(1x,a)') '### error : data error'
        write(ifle,'(1x,a)') 'gastype =',gastype
        write(ifle,'(1x,a)') 'it should be, 1: gray_gas, 2: real_gas'
        goto 9999
      end if
!
      radmat(1,nd) = gabsorb
      radmat(2,nd) = pabsorb
      radmat(3,nd) = pscatter
      radmat(4,nd) = scabeta
      radmat(5,nd) = ptcdia
!
      radfludflag(1,nd) = gastype
      radfludflag(2,nd) = scatype
      radfludflag(3,nd) = ptccal
      radfludflag(4,nd) = fgpflag
!
      rg_mark=0
      rg_markx=0
      if(gastype==2) then
        rg_mark=1
        rg_markx=1
      endif
!
      goto 1000
 1001 continue
!
      ical_sld=0
      lrot=.false.
      do nd=1,nflud
      if(ishaft(nd)==1) then
        ical_sld=1
        lrot=.true.
      endif
      enddo
!
!
!      if(nd>1.and.iMvmsh/=0) then
!        call FFRABORT(1,'ERR:Moving mesh NOT support mutil-material')
!      endif
!----------------------------------------------
! --- solid material
!----------------------------------------------
      nd=0
      rewind ifli
  200 continue
      read(ifli,solid,iostat=ios)
      if( ios.lt.0 ) goto 201
      nd=nd+1
      call nml_errmsg0(ifle,ios,'solid',ierr1)
      if( ierr1.ne.0 ) then
        write(ifle,*) 'sequence no. of the namelist =',nd
        call FFRABORT(1,'')
      endif
      goto 200
  201 continue
      nsold=nd
      NOMAT=nflud+nsold
      allocate(dpmxv(0:NOMAT),incmp_prs(0:NOMAT))
      dpmxv(:)=0.d0
      incmp_prs(:)=0.d0
      if( nd.lt.1 ) goto 2001
!
!--< 2.1 allocate arrays >--
!
      allocate( nosld (nsold),
     &          rsld  (nsold),
     &          cpsld (nsold),
     &          PSstart(nsold),
     &          relax_sh(nsold),
     &          ivscsldMHD(nsold),
     &          icnvsldMHD(nsold),
     &          relaxsldMHD(nsold),
     &          GMHDA_sld(nsold),
     &          GMHDJ0_sld(nsold),
     &          GMHDJe_sld(nsold),
     &          rmdsld(nsold), stat=ierr1)
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error-2 : allocation failed'
        call FFRABORT(1,'')
      endif
      nosld =0
      rsld  =0.d0
      cpsld =0.d0
      rmdsld=0.d0
      relax_sh=1.d0
      ivscsldMHD=cal
      icnvsldMHD=up1st
      relaxsldMHD=1.d0
      GMHDA_sld=0
      GMHDJ0_sld=0
      GMHDJe_sld=0
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
      Poisson_start=iundef
!
      read(ifli,solid,iostat=ios)
      if( ios.lt.0 ) goto 2001
      call nml_errmsg0(ifle,ios,'solid',ierr1)
      if( ierr1.ne.0 ) call FFRABORT(1,'')
      nd=nd+1
!
      if( IMAT_U.eq.iundef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'lack of data : IMAT_U'
        call FFRABORT(1,'')
      elseif( IMAT_U.ge.0 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'it must be < 0'
        call FFRABORT(1,'')
      endif
      do 2100 n=1,nd-1
      if( IMAT_U.eq.nosld(n) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'IMAT_U =',IMAT_U
        write(ifle,*) 'this IMAT_U. already specified in previous line'
        write(ifle,*) 'sequence no. of the line =',n
        call FFRABORT(1,'')
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
        write(ifle,*) 'lack of data : Poisson_start in namelist/fluid'
        write(ifle,*) 'Poisson_start is starting iter for Poisson'
        call FFRABORT(1,'')
      elseif( Poisson_start.le.-1 ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'Poisson_start =',Poisson_start
        write(ifle,*) 'it must be >= 0'
        call FFRABORT(1,'')
      endif
!
      if(unrelax_T.gt.1.d0.or.unrelax_T.lt.0.d0) then
        write(ifle,*) '### error : data error, unrelax_T=',unrelax_T
        write(ifle,*) 'Under-relaxation faceor must between 1 and 0'
        call FFRABORT(1,'')
      endif
!
      nosld (nd)=IMAT_U
      rsld  (nd)=rho
      cpsld (nd)=cp
      rmdsld(nd)=cndct
      PSstart(nd)=Poisson_start
      relax_sh(nd)=unrelax_T
      IF(imhd>0) then
        call nml_listno(8,cnvlst,cnv_MHD,icnvsldMHD(nd))
        GMHDA_sld(nd)=Gauge_MHD_A
        GMHDJ0_sld(nd)=Gauge_MHD_J0
        GMHDJe_sld(nd)=Gauge_MHD_Je

        if(icnvsldMHD(nd)==engcsv)  
     &    call FFRABORT(1,'ERR:Scheme [eng] only for velocity')
        call nml_chkchr0(ifle,'cnv_MHD',cnv_MHD,icnvsldMHD(nd),ierr2)
        relaxsldMHD(nd)=unrelax_MHD
      endif
      if(ierr2/=0) goto 9999
      if((icnvsldMHD(nd)==up3rd.or.icnvsldMHD(nd)==up2nd))then
        write(ifle,'(1x,2A)') 
     &  "WRN: MHD ",
     &  "MUST be  [1st] & [c2d] for MHD equation."
        call FFRABORT(1,"MSG: Reset MHD ")
      endif
      ivscsldMHD(nd)=ixl(diff_MHD)
!
      goto 2000
 2001 continue
!
!      allocate(MAT_ALL(nflud+nsold),MAT_S(nflud+nsold))
      allocate(relaxp(0:nflud+nsold),
     &         relaxC(0:nflud+nsold),
     &         relaxCi(0:nflud+nsold),
     &         relaxC_step(0:nflud+nsold),
     &         relaxh(0:nflud+nsold),
     &         relaxrs(0:nflud+nsold),
     &         relaxv(0:nflud+nsold),
     &         relaxrc(0:nflud+nsold),
     &         relax_MHD(0:nflud+nsold),
     &         GAUGE_A(0:nflud+nsold),
     &         GAUGE_J0(0:nflud+nsold),
     &         GAUGE_Je(0:nflud+nsold)
     &         )
      allocate(relaxys(0:nflud+nsold))
      
      allocate(OMGA(0:nflud+nsold)
     &        ,icnvMHD(0:nflud+nsold)
     &        ,ivscMHD(0:nflud+nsold)
     &         )
!
      relaxp(:)=1.d0
      relaxC(:)=1.d0
      relaxCi(:)=1.d0
      relaxC_step(:)=1.d0
      relaxrs(:)=1.d0
      relaxv(:)=1.d0
      relaxrc(:)=0.d0
      relaxh(:)=1.d0
      relaxys(:)=1.d0
      OMGA(:)=0.d0
      icnvMHD(:)=up1st
      ivscMHD(:)=cal
      relax_MHD(:)=1.d0
      GAUGE_A(:)=0
      GAUGE_J0(:)=0
      GAUGE_Je(:)=0
!      
!      nm=0
!      do nd=1,nflud
!        nm=nm+1
!        MAT_ALL(nm)=nofld(nd)
!      enddo
!
!      do nd=1,nsold
!        nm=nm+1
!        MAT_ALL(nm)=nosld(nd)
!      enddo
!      KMAT_S=nm
!
      call disp_head
      if( ierror/=0 ) then      ! final error check
        write(ifle,*) "total error : ",ierror,modnam,subnam
      endif
      if(my_rank==0) write(ifll,"(1x,80('-'))")
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      contains
!
!< #1.1 integer to logical > >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      function lxi(i)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      logical :: lxi
      integer,intent(in) :: i
                   lxi=.false.
      if( i.gt.0 ) lxi=.true.
      end function lxi
!
!< #1.2 logical to integer > >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      function ixl(l)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      integer :: ixl
      logical,intent(in) :: l
              ixl=0
      if( l ) ixl=1
      end function ixl
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_Sutherland       ! Sutherland
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(.not.nml_comp_eq("mu",mu,undef,nd,ifle,ierror)) then
        if(nml_comp_gt("mu",mu,0.0d0,nd,ifle,ierror)) then
          sthrmu(nd) = mu
        else
          if(my_rank==0) then
            write(ifle,'(a)') 
     &      'EER: mu MUST be great than 0.d0 for Sutherland'
            call FFRABORT(1,'ERR:')
          endif
        endif
      else
          if(my_rank==0) then
            write(ifle,'(a)') 'EER: mu MUST be defined for Sutherland'
            call FFRABORT(1,'ERR:')
          endif
      endif
!
      if(.not.nml_comp_eq("t0",t0,undef,nd,ifle,ierror)) then
        if(nml_comp_gt("t0",t0,0.0d0,nd,ifle,ierror)) then
          sthrt(nd) = t0
        else
          if(my_rank==0) then
            write(ifle,'(a)')  
     &      'EER: t0 MUST be great than 0.d0 for Sutherland'
            call FFRABORT(1,'ERR:')
          endif
        endif
      else
          if(my_rank==0) then
            write(ifle,'(a)') 'EER: t0 MUST be defined for Sutherland'
            call FFRABORT(1,'ERR:')
          endif
      endif
      if(.not.nml_comp_eq("c",c,undef,nd,ifle,ierror)) then
        sthrc(nd) = c
      else
          if(my_rank==0) then
            write(ifle,'(a)') 'EER: c MUST be defined for Sutherland'
            call FFRABORT(1,'ERR:')
          endif
      endif
!
      if(nml_comp_gt("Prandtl",Prandtl,0.0d0,nd,ifle,ierror)) then
        prlmnr(nd) = Prandtl
        r_prlmnr(nd) = 1.d0/Prandtl
      endif
      if(nml_comp_gt("Schmidt",Schmidt,0.0d0,nd,ifle,ierror)) then
        sclmnr(nd) = Schmidt
        r_sclmnr(nd) = 1.d0/Schmidt
      endif
      end subroutine read_Sutherland
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_constant           ! constant
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(.not.nml_comp_eq("mu",mu,undef,nd,ifle,ierror)) then
        if(nml_comp_gt("mu",mu,0.0d0,nd,ifle,ierror))
     &    sthrmu(nd) = mu
      endif
      if(nml_comp_gt("Prandtl",Prandtl,0.0d0,nd,ifle,ierror)) then
        prlmnr(nd) = Prandtl
        r_prlmnr(nd) = 1.d0/Prandtl
      endif
      if(nml_comp_gt("Schmidt",Schmidt,0.0d0,nd,ifle,ierror)) then
        sclmnr(nd) = Schmidt
        r_sclmnr(nd) = 1.d0/Schmidt
      endif
      sthrt(nd) = sthrt_i
      sthrc(nd) = sthrc_i
      end subroutine read_constant
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_Simplified       ! simplified
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      prlmnr(nd) = 0.75d0         ! always constant
      r_prlmnr(nd) = 1.d0/prlmnr(nd)
      end subroutine read_Simplified
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_MK
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      MIX_KNI=1
      end subroutine read_MK
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_Chung
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      CHUN=1
      end subroutine read_Chung
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_Sutherland2       ! Sutherland
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(.not.nml_comp_eq("mu2",mu2,undef,nd,ifle,ierror)) then
        if(nml_comp_gt("mu2",mu2,0.0d0,nd,ifle,ierror))
     &    sthrmu2(nd) = mu2
      endif
      if(.not.nml_comp_eq("t02",t02,undef,nd,ifle,ierror)) then
        if(nml_comp_gt("t02",t02,0.0d0,nd,ifle,ierror))
     &    sthrt2(nd) = t02
      endif
      if(.not.nml_comp_eq("c2",c2,undef,nd,ifle,ierror))
     &  sthrc2(nd)= c2
      if(nml_comp_gt("Prandtl2",Prandtl2,0.0d0,nd,ifle,ierror)) then
        prlmnr2(nd)= Prandtl2
        r_prlmnr2(nd)= 1.d0/Prandtl2
      endif
      if(nml_comp_gt("Schmidt2",Schmidt2,0.0d0,nd,ifle,ierror)) then
        sclmnr2(nd)= Schmidt2
        r_sclmnr2(nd)= 1.d0/Schmidt2
      endif
      end subroutine read_Sutherland2
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_constant2           ! constant
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(.not.nml_comp_eq("mu2",mu2,undef,nd,ifle,ierror)) then
        if(nml_comp_gt("mu2",mu2,0.0d0,nd,ifle,ierror))
     &    sthrmu2(nd) = mu2
      endif
      if(nml_comp_gt("Prandtl2",Prandtl2,0.0d0,nd,ifle,ierror)) then
        prlmnr2(nd)= Prandtl2
        r_prlmnr2(nd)= 1.d0/Prandtl2
      endif
      if(nml_comp_gt("Schmidt2",Schmidt2,0.0d0,nd,ifle,ierror)) then
        sclmnr2(nd) = Schmidt2
        r_sclmnr2(nd) = 1.d0/Schmidt2
      endif
      sthrt2(nd) = sthrt_i
      sthrc2(nd) = sthrc_i
      end subroutine read_constant2
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_Simplified2       ! simplified
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      prlmnr(nd) = 0.75d0         ! always constant
      r_prlmnr(nd) = 1.d0/prlmnr(nd)
      end subroutine read_Simplified2
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine disp_head
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(my_rank==0) write(ifll,"(1x,80('-'))")
      do n=1,nflud
        if(my_rank==0) write(ifll,"(1x,a,i3)")
     &      "Domain number              : ",n
        if(my_rank==0) write(ifll,"(1x,2a)")
     &      "transport model            : ",mulst(isthr(n))
        if( isthr(n)==1 ) then          ! Sutherland
          call disp_Sutherland
        elseif( isthr(n)==2 ) then              ! constant
          call disp_constant
        elseif( isthr(n)==3 ) then              ! simplified
          call disp_Simplified
        endif
      enddo
      end subroutine disp_head
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine disp_Sutherland
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(my_rank==0) write(ifll,"(3x,a,e10.3,a)")
     &  "viscosity at 'T0' [Pa s] : ",sthrmu(n)," (mu0)"
      if(my_rank==0) write(ifll,"(3x,a,f10.3,a)")
     &  "temperature [K]          : ",sthrt(n)," (T0)"
      if(my_rank==0) write(ifll,"(3x,a,f10.3,a)")
     &  "constant value           : ",sthrc(n)," (C)"
      if(my_rank==0) write(ifll,"(3x,a,f10.3,a)")
     &  "Schmidt number           : ",sclmnr(n)," (Sc)"
      if(my_rank==0) write(ifll,"(3x,a,f10.3,a)")
     &  "Prandtl number           : ",prlmnr(n)," (Pr)"
      if(my_rank==0) write(ifll,"(5x,2a)")
     &  "viscosity    : m = mu0*(T0+C)/(T+C)*(T/T0)^(3/2)",
     &  "  (T:temperature)"
      if(my_rank==0) write(ifll,"(5x,2a)")
     &  "conductivity : l = Cp*m/Pr",
     &  "  (Cp:isopiestic specific heat)"
      if(my_rank==0) write(ifll,"(5x,2a)")
     &  "diffusivity  : D = m/(r*Sc)",
     &  "  (r:mass density)"
      end subroutine disp_Sutherland
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine disp_constant
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(my_rank==0) write(ifll,"(3x,a,e10.3,a)")
     &  "viscosity [Pa s]         : ",sthrmu(n)," (mu0)"
      if(my_rank==0) write(ifll,"(3x,a,f10.3,a)")
     &  "Schmidt number           : ",sclmnr(n)," (Sc)"
      if(my_rank==0) write(ifll,"(3x,a,f10.3,a)")
     &  "Prandtl number           : ",prlmnr(n)," (Pr)"
      if(my_rank==0) write(ifll,"(5x,a)")
     &  "viscosity    : m = mu0"
      if(my_rank==0) write(ifll,"(5x,2a)")
     &  "conductivity : l = Cp*m/Pr",
     &  "  (Cp:isopiestic specific heat)"
      if(my_rank==0) write(ifll,"(5x,2a)")
     &  "diffusivity  : D = m/(r*Sc)",
     &  "  (r:mass density)"
      end subroutine disp_constant
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine disp_Simplified
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(my_rank==0) write(ifll,"(3x,a,f10.3,a)")
     &  "Prandtl number           : ",prlmnr(n)," (Pr)"
      if(my_rank==0) write(ifll,"(5x,a)")
     &  "viscosity    : m = a*Pr"
      if(my_rank==0) write(ifll,"(5x,2a)")
     &  "conductivity : l = a*Cp",
     &  "  (Cp:isopiestic specific heat)"
      if(my_rank==0) write(ifll,"(5x,2a)")
     &  "diffusivity  : D = a/(r*Le)",
     &  "  (r:mass density, Le:Lewis number)"
      if(my_rank==0) write(ifll,"(10x,a)")
     &  "a = 2.58e-5*(T/298)^0.7"
      end subroutine disp_Simplified
!
      end subroutine inputdata
!
      end module module_material

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_gravity
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(20),parameter,private :: modnam='(module_gravity)'
!
!
!< module for gravity >
!
      real(8) :: ggg(3)=(/0.d0,0.d0,0.d0/)
      real(8) :: ramb=1.d0
      real(8) :: tamb=273.15d0, beta=0.21d-3,      ! water
     &                          beta2=1.d0/273.d0 !air
      real(8) :: pamb=101325.d0
      real(8),allocatable :: ysamb(:)

      real(8) :: tamb2=273.15d0
      real(8) :: ramb2=1.d0
      integer:: buoyancy=0,i
      integer:: iave_rho_T
!
! ggg  : x,y,z components of gravity
! ramb : ambient density [kg/m^3]
! tamb : ambient temperature [K]
! beta : expansion coefficient [1/K]
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,ncomp,cntlnam,lcomp,
     &    E2P,ierror)
!=================================================
!
      implicit none
!
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle,ncomp
      character(*),intent(in) :: cntlnam
      logical,intent(in)  :: lcomp,E2P
      integer,intent(out) :: ierror
      
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      real(8),parameter :: undef=-huge(1.d0)
      real(8)  :: betax
      integer :: iset=0
      integer :: ios,ierr1
!
! --- [namelist]
!
      integer,parameter  :: mcomp=200
      real(8) :: g(3)
      real(8) :: rho,t,rho2,t2,P,ys(mcomp)
      integer :: ave_rho_T=0
      namelist /gravity/ g,rho,t,beta,beta2,rho2,t2,ave_rho_T,P,ys
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      g   =ggg
      rho =ramb
      t   =tamb
      p   =pamb
!
      rewind ifli
      read(ifli,gravity,iostat=ios)
      if( ios.lt.0 ) return
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
      if(.not.(ave_rho_T==0.or.ave_rho_T==1.or.
     &         ave_rho_T==2.or.ave_rho_T==3.or.
     &         ave_rho_T==4.or.
     &         ave_rho_T==20)
     &     ) then
          write(ifle,*) '### error : ave_rho_T error'
          write(ifle,*) 'ave_rho_T = ',ave_rho_T
          write(ifle,*) 'ave_rho_T=0,1,2,3, or 20 '
          write(ifle,'(1X,A)') 
     &  'MSG: ave_rho_T=0: using [rho] or [t] as reference parameter'
          write(ifle,'(6X,A,/,6x,A)') 
     &  'MSG: g*(dens-rho) for [zero] or [comp] flow ',
     &  'MSG: g*beta*(tmp-t)*dens for [incmop flow]'
          write(ifle,'(1X,2A)') 
     &  'MSG: ave_rho_T=1: Averaging all [rho] ',
     &  'or [tmp] as reference parameter'
          write(ifle,'(6X,A,/,6x,A)') 
     &  'MSG: g*(dens-rho_ave) for [zero] or [comp] flow ',
     &  'MSG: g*beta*(tmp-t_ave) for [incmop flow]'
          write(ifle,'(1X,A)') 
     & 'MSG: ave_rho_T=2: (g*dens) for [incomp] [comp] and [zero] flow'
          write(ifle,'(1X,A,/,6x,2A)') 
     &  'MSG: ave_rho_T=3: Averaging all [rho] as reference parameter',
     &    'MSG: g*(dens-rho_ave) for [incomp] flow ',
     &    'MSG: while, density=1,2,3'
          write(ifle,'(1X,A)') 
     &  'MSG: ave_rho_T=4:  g*(dens-rho_standard) for [comp] flow'
          write(ifle,'(1X,A,/,6x,A,/,6x,A)') 
     &  'MSG: For 2-pahse:',
     &  'MSG: ave_rho_T=20: g*(dens-rho) only for [incomp] flow  ',
     &  'MSG: ave_rho_T=0:  g*(t-tmp) only for [incomp] flow'
          
          goto 9999
      endif
!
      if(E2P) then
        if( rho.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'rho = ',rho
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
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
        if( rho2.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'rho2 = ',rho2
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
        if( t2.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 't2 = ',t2
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
        if( beta2.le.0.d0 ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'beta2 = ',beta2
          write(ifle,*) 'it must be > 0'
          goto 9999
        endif
      endif
!
      ggg =g
!
      ramb=rho
      tamb=t
      pamb=p
      ALLOCATE(ysamb(ncomp))
      do i=1,ncomp
      ysamb(i)=ys(i)
      enddo

      if(E2P) then
        ramb2=rho2
        tamb2=t2
      endif
      if(abs(ggg(1)).gt.0.d0.or.
     &   abs(ggg(2)).gt.0.d0.or.
     &   abs(ggg(3)).gt.0.d0) then
         buoyancy=1
      endif
      iave_rho_T=ave_rho_T
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_gravity

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_io
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
!
!< module for logical unit no. & file name >
!
!-< 1.logical unit no. for control data (input) >-
!
      integer,save      :: ifli=1
      character(LEN=9),save :: cntlnam='fflow.ctl'
!
!-< 2.logical unit no. for log file & error message (output) >-
!
      integer,save :: ifll=6, ifle=6
!
!-< 3.logical unit no. & name for external files >-
!
      character(LEN=11),parameter,private :: modnam='(module_io)'
      integer,parameter              :: lenfnm=200
      integer,parameter,private      :: mfile=29
      character(LEN=lenfnm),parameter,private :: blank=' '
      integer          ,save,private :: iflno (-1:mfile)
      character(LEN=lenfnm),save,private :: filnam(-1:mfile)
      character(LEN=11)    ,save,private :: format(-1:mfile)
      character(LEN=7)     ,save,private :: status(-1:mfile)
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
      integer,private  :: iset=0
      character(LEN=4) :: gdformat
      real(8) :: gdScale,angle_ssf,tolerance
!      logical :: GFrslt,FFrslt,FVrslt
!
!     GF:GF format  FV : Gridgen's FV format
!
      logical :: dflxTmpFFlag
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ierror)
!=================================================
      implicit none
!-------------------------
! --- [dummy arguments]
!-------------------------
      integer,intent(in)  :: ifli
      integer,intent(out) :: ierror
!-------------------------
! --- [local entities]
!-------------------------
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=9) ,parameter :: formatted='formatted'
      character(LEN=11),parameter :: unformatted='unformatted'
      character(LEN=3) ,parameter :: old='old', new='new'
      character(LEN=7) ,parameter :: unknown='unknown'
      integer :: i,ios=0,no
      logical :: fexist=.false.
!-------------------------
! --- [namelist]
!-------------------------
      character(lenfnm) :: log,errmsg
      character(lenfnm) :: cell,vertex,boundary,initial
      character(lenfnm) :: p_initial,
     &                     p_restart,
     &                     p_result
      character(lenfnm) :: restart,result,history,force,anim
      character(lenfnm) :: scInp,scVrt,scCel,scBnd
      character(lenfnm) :: fluent_grd
      character(lenfnm) :: scryu_grd
      character(lenfnm) :: nastran_grd
      character(lenfnm) :: cdcl,probe,fforce
      character(lenfnm) :: ffrgrid
      character(lenfnm) :: move_grid
      character*1 :: fgConv,ffrgridform
      character*3 :: multi_material
      namelist /files/ log,errmsg,
     &                 cell,vertex,boundary,
     &                 p_initial,
     &                 p_restart,
     &                 p_result,
     &                 initial,restart,result,
     &                 history,force,anim,
     &                 gdformat,
     &                 gdScale,angle_ssf,tolerance,
     &                 scInp,scVrt,scCel,scBnd,
     &                 fluent_grd,
     &                 scryu_grd,
     &                 nastran_grd,
     &                 ffrgrid,ffrgridform,fgConv,
     &                 cdcl,probe,fforce
     &                 ,multi_material
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!-------------------------
!-< 1. Input namelist >-
!-------------------------
      log     =blank
      errmsg  =blank
      cell    =blank
      vertex  =blank
      boundary=blank
!
      initial =blank
      restart =blank
      result  =blank
!
      p_initial=blank
      p_restart=blank
      p_result=blank
!
      history =blank
      force   =blank
      anim    =blank
      gdformat='GF'
      move_grid='move_grid'
!      rsltgf=blank
!      rsltfv=blank
!      GFrslt=.false.
!      FFrslt=.false.
!      FVrslt=.false.
      gdScale=1.d0
      scInp  ='star.inp'
      scVrt  ='star.vrt'
      scCel  ='star.cel'
      scBnd  ='star.bnd'
      fluent_grd = blank
      scryu_grd  = blank
      nastran_grd= blank
      ffrgrid=blank
      fgConv='n'
      ffrgridform='U'
      cdcl= blank
      probe=blank
      fforce=blank
!
!
!
      inquire(file=cntlnam,exist=fexist)
      if(.not.fexist) then
        write(ifle,*) 'Control file [',cntlnam,'] NOT exist'
        call FFRABORT(1,'Please check control file name')
      endif
      open(ifli,file=cntlnam,form='FORMATTED')
!----------------------------------------------
! --- control file (cntlnam=fflow.ctl or fort.1,fort.88)
!----------------------------------------------
      rewind ifli
      read (ifli,files,iostat=ios)
      if( ios.ne.0 ) then
        if( ios.gt.0 ) write(*,*)
     &         '### error : read error in module_total'
        if( ios.lt.0 ) write(*,*) '### error : end of file'
        write(*,*) 'namelist = files'
        write(*,*) 'iostat =',ios
        goto 9999
      endif
!
      if(gdformat=='FV') then
        call errmsg1(vertex  ,'vertex')
      else if(gdformat=='SC') then
        call errmsg1(scInp  ,'StarCD-Input')
        call errmsg1(scVrt  ,'StarCD-Vertex')
        call errmsg1(scBnd  ,'StarCD-Boundary')
        call errmsg1(scCel  ,'StarCD-Cell')
      else if(gdformat=='GB') then
        call errmsg1(fluent_grd  ,'fluent_grd')
      else if(gdformat=='SY') then
        call errmsg1(scryu_grd  ,'scryu_grd')
      else if(gdformat=='GF') then
        call errmsg1(ffrgrid ,'ffrgrid')
      else if(gdformat=='NST') then
        call errmsg1(nastran_grd,'nastran_grd')
      else
        call errmsg1(cell    ,'cell')
        call errmsg1(vertex  ,'vertex')
        call errmsg1(boundary,'boundary')
      end if
!      if(result.ne.' ') FFrslt=.true.
!      if(rsltgf.ne.' ') GFrslt=.true.
!      if(rsltfv.ne.' ') FVrslt=.true.
      if( ierror.ne.0 ) goto 9999
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
! --- log file: no=-1
!
      no=-1
      call setspc(-1,log,formatted,unknown)
      if( log.ne.blank ) ifll=iflno(no)
!
! --- error message file: no=0
!
      no=0
      call setspc(no,errmsg,formatted,unknown)
      if( errmsg.ne.blank ) ifle=iflno(no)
!
! ---
!
      call setspc( 1, cell    ,unformatted,old)
      call setspc( 2, vertex  ,unformatted,old)
      call setspc( 3, boundary,unformatted,old)
      call setspc( 4, p_initial ,unformatted,old)
      call setspc( 5, initial ,unformatted,old)
      call setspc( 6, restart ,unformatted,unknown)
      call setspc( 7, result  ,unformatted,unknown)
      call setspc( 8, history ,unformatted,unknown)
      call setspc( 9, force   ,unformatted,unknown)
      call setspc( 10,move_grid  ,unformatted,unknown)
!      call setspc( 11,restart_particle,unformatted,unknown)
      call setspc( 12,scInp,   formatted,old)
      call setspc( 13,scVrt,   formatted,old)
      call setspc( 14,scBnd,   formatted,old)
      call setspc( 15,scCel,   formatted,old)
      call setspc( 16,fluent_grd,formatted,old)
      call setspc( 17,anim    ,unformatted,unknown)
      call setspc( 18,scryu_grd,formatted,old)
      call setspc( 19,cdcl,   formatted,old)
      call setspc( 20,probe,  formatted,old)
      call setspc( 21,fforce,formatted,old)
      call setspc( 22,nastran_grd,formatted,old)
      call setspc( 28,p_restart,unformatted,unknown)
      call setspc( 29,p_result,unformatted,unknown)
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
      subroutine errmsg1(vnam1,vnam2)
!=================================================
      implicit none
      character(*) :: vnam1,vnam2
      if( vnam1.ne.' ' ) return
      write(*,*) '### error : data error'
      write(*,*) 'lack of data : ',vnam2
      ierror=1
      end subroutine errmsg1
!
!< #1.2 set specifier >
!=================================================
      subroutine setspc(no,fin,fm,sta)
!=================================================
      implicit none
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
      implicit none
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
      implicit none
      character(*),intent(in)  :: cnam
      integer     ,intent(out) :: ifl
      character(*),intent(out) :: fnam
      character(LEN=8),parameter :: subnam='(getfil)'
      integer :: jfl
!
      call modutl_chkset('>',2,iset,modnam//subnam)
                                jfl=-999
      if( cnam.eq.'cell'    )   jfl=1
      if( cnam.eq.'vertex'  )   jfl=2
      if( cnam.eq.'boundary')   jfl=3
      if( cnam.eq.'p_initial'  )   jfl=4
      if( cnam.eq.'initial' )   jfl=5
      if( cnam.eq.'restart' )   jfl=6
      if( cnam.eq.'result'  )   jfl=7
      if( cnam.eq.'history' )   jfl=8
      if( cnam.eq.'force'   )   jfl=9
      if( cnam.eq.'move_grid')  jfl=10
      if( cnam.eq.'scInp'   )   jfl=12
      if( cnam.eq.'scVrt'   )   jfl=13
      if( cnam.eq.'scBnd'   )   jfl=14
      if( cnam.eq.'scCel'   )   jfl=15
      if( cnam.eq.'fluent_grd') jfl=16
      if( cnam.eq.'anim')       jfl=17
      if( cnam.eq.'scryu_grd')  jfl=18
      if( cnam.eq.'cdcl')       jfl=19
      if( cnam.eq.'probe')      jfl=20
      if( cnam.eq.'fforce')     jfl=21
      if( cnam.eq.'nastran_grd')jfl=22
      if( cnam.eq.'eqrate_wk')  jfl=23
      if( cnam.eq.'statis')     jfl=24
      if( cnam.eq.'statis.tmp') jfl=25
      if( cnam.eq.'distance')   jfl=26
      if( cnam.eq.'geom')   jfl=27
      if( cnam.eq.'p_restart') jfl=28
      if( cnam.eq.'p_result')  jfl=29
      if( jfl.eq.-999 ) then
        write(*,*) '### error : program error -1- ',modnam,subnam
        write(*,*) 'cnam = ',cnam
        CALL FFRABORT(1,' module_io/getfil')
      endif
      ifl=iflno(jfl)
      fnam=filnam(jfl)
      end subroutine getfil
!
      end module module_io
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_les
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
!
      character(16),parameter,private :: modnam='(module_les)'
!
!< module for LES model >
!
      real(8) :: csles=0.1d0
      real(8) :: apls =25.d0
      real(8) :: dlalpha=2.d0
      real(8) :: ign_iv=0
      integer:: ISTART=1000
      integer,save             :: iuvw_ave_rms_re
      integer,save             :: ip_ave,it_ave,imu_ave
      integer,save,allocatable :: icomp_ave(:),irans_ave(:)
      integer,save             :: it_rms,ip_rms,imu_rms
      integer,save,allocatable :: icomp_rms(:),irans_rms(:)
      integer,save             :: iuvwt_re,imaxmin
      integer,save,allocatable :: iuvw_rans_re(:)
      integer,save,allocatable :: iuvwc_re(:)
!
      integer,save             :: iuvw_ave_rms_re_c
      integer,save             :: ip_ave_c,it_ave_c
      integer,save,allocatable :: icomp_ave_c(:),irans_ave_c(:)
      integer,save             :: it_rms_c,ip_rms_c
      integer,save,allocatable :: icomp_rms_c(:),irans_rms_c(:)
      integer,save             :: iuvwt_re_c
      integer,save,allocatable :: iuvw_rans_re_c(:)
      integer,save,allocatable :: iuvwc_re_c(:)
!
! csles : Smagorinsky's constant
! ISTART: Start iteration for average value
! iu_ave: flag for if/not computing averaged value for LES (velocity)
! ip_ave: flag for if/not computing averaged value for LES (pressure)
! icomp_ave: flag for if/not computing averaged value for LES (speciese)
! irans_ave: flag for if/not computing averaged value for LES (scalar )
! iu_ave=0/1: 0=>not calculating ; 1=>calculating
!
      integer,parameter        :: msuf_out=11
      integer,save             :: SDOT_suf=0,WDOT_gas=0,
     &                            depoeth_spd=0,
     &                            molefr_suf=0,
     &                            mol_rslt=0,
     &                            blk_thick=0
      integer,save,allocatable :: rat_E_S(:,:)       !(msuf_out,)
      integer,save             :: num_ratout=0
      integer,save             :: ista_momery=0,n_ave=0
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!===============================================================
      subroutine inputdata
     &          (ifli,ifll,ifle,mmcomp,mneq,ncompall,ncomp,nneq,
     &           my_rank,cntlnam,nrans,spinam,net_vreq,ierror)
!===============================================================
      implicit none
! --- MXMATERIAL
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifli,ifll,ifle,nrans,ncomp,mmcomp,
     &                           nneq,mneq,ncompall,my_rank
      character(*),intent(in) :: cntlnam,spinam(mmcomp)
      real*8, intent(in)      :: net_vreq(mneq*ncompall)
      integer,intent(out)     :: ierror
!
! --- [local entities]
!
      integer,parameter :: mcomp=200,mrans=50
      real*8,parameter  :: SML=1.d-10
      character(LEN=11),parameter :: subnam='(inputdata)'
      integer :: iset=0
      integer :: ios,ierr1,s,is,ieq
      real*8  :: dum
!
! --- [namelist]
!
      integer:: i,j
      real(8) :: cs=0.1d0
      integer:: NSTART=1000
      integer:: uvw_ave_rms_re=1,
     &          average_p=0,average_t=0,
     &          average_comp(mcomp),average_rans(mcomp),
     &          rms_t=0,rms_p=0,
     &          rms_comp(mcomp),rms_rans(mcomp),
     &          re_uvw_t=0,maxmin_uvw_p,
     &          re_uvw_rans(mrans),
     &          statis_file
      integer:: average_Mut=0,rms_Mut=0
      integer::  re_uvw_comp(mcomp)
!
      integer:: 
     &          surface_SDOT=0,
     &          GAS_WDOT=0,
     &          depos_etch_speed=0,
     &          surface_mole_fraction=0,gas_mole_fraction=0,
     &          bulk_thickness=0
      character(LEN=80) :: RATE_EQ_SpecNAME(msuf_out)
      character(LEN=1)  :: uder
      character(LEN=80) :: cdum
!
      namelist /les/ cs,NSTART,
     &               statis_file,
     &               dlalpha,ign_iv,
     &               uvw_ave_rms_re,
     &               average_p,maxmin_uvw_p,average_t,
     &               average_Mut,rms_Mut,
     &               average_comp,average_rans,
     &               rms_t,rms_p,
     &               rms_comp,rms_rans,
     &               re_uvw_t,
     &               re_uvw_comp,
     &               re_uvw_rans,
     &               surface_SDOT,gas_WDOT,depos_etch_speed,
     &               surface_mole_fraction,gas_mole_fraction,
     &               RATE_EQ_SpecNAME,
     &               bulk_thickness
!
      allocate(icomp_ave(ncomp),
     &         irans_ave(nrans),
     6         icomp_rms(ncomp),
     &         irans_rms(nrans),
     &         iuvw_rans_re(nrans)
     &         )
!
      allocate(iuvwc_re(ncomp))
!
      allocate(icomp_ave_c(ncomp),
     &         irans_ave_c(nrans),
     6         icomp_rms_c(ncomp),
     &         irans_rms_c(nrans),
     &         iuvwc_re_c(ncomp),
     &         iuvw_rans_re_c(nrans))
!
      allocate(rat_E_S(msuf_out,ncompall))
!
!
      iuvw_ave_rms_re=1
      ip_ave=0
      imu_ave=0
      it_ave=0
      icomp_ave(:)=0
      irans_ave(:)=0
      ip_rms=0
      imu_rms=0
      it_rms=0
      icomp_rms(:)=0
      irans_rms(:)=0
      iuvwt_re=0
      iuvwc_re(:)=0
      imaxmin=0
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
      iuvwc_re_c(:)=0
      iuvw_rans_re_c(:)=0
!
      ista_momery=0
!
      SDOT_suf=0
      blk_thick=0
      WDOT_gas=0
      depoeth_spd=0
      molefr_suf=0
      mol_rslt=0
      rat_E_S(:,:)=0
      num_ratout=0
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
 1000 continue
      rewind ifli
      cs=csles
      if(cs<0.d0) then
        cs=0.1d0
      endif
      if(dlalpha<0.d0) then
        dlalpha=2.d0
      endif
!
      if(ign_iv>3) then
        ign_iv=0
      endif
!
      NSTART=1000
      uvw_ave_rms_re=1
      average_p=0
      average_t=0
      average_Mut=0
      average_comp(:)=0
      average_rans(:)=0
      rms_p=0;rms_t=0;rms_Mut=0
      rms_comp(:)=0;rms_rans(:)=0
      re_uvw_t=0
      re_uvw_comp(:)=0
      maxmin_uvw_p=0
      re_uvw_rans(:)=0
      statis_file=1
!
      surface_SDOT=0
      bulk_thickness=0
      GAS_WDOT=0
      depos_etch_speed=0
      surface_mole_fraction=0
      gas_mole_fraction=0
      RATE_EQ_SpecNAME(:)=''
      dlalpha=2.d0
!
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
      if(statis_file==0) then
        ista_momery=1
      elseif(statis_file==1) then
        ista_momery=0
      else
        call FFRABORT(1,'ERR: statis_file/=0 or /=1 in &les')
      endif
      if(uvw_ave_rms_re==0.or.uvw_ave_rms_re==1) then
        iuvw_ave_rms_re=uvw_ave_rms_re
      else
        call FFRABORT(1,'ERR: uvw_ave_rms_re must bc 0 or 1')
      endif
      if(average_p==0.or.average_p==1) then
        ip_ave=average_p
      else
        call FFRABORT(1,'ERR: average_p must bc 0 or 1')
      endif
      if(average_Mut==0.or.average_Mut==1) then
        imu_ave=average_Mut
      else
        call FFRABORT(1,'ERR: average_Mut must bc 0 or 1')
      endif
      if(average_t==0.or.average_t==1) then
        it_ave=average_t
      else
        call FFRABORT(1,'ERR: average_t must bc 0 or 1')
      endif
      do i=1,ncomp
        if(average_comp(i)==0.or.average_comp(i)==1) then 
          icomp_ave(i)=average_comp(i)
        else
          call FFRABORT(1,'ERR: average_comp(i) must bc 0 or 1')
        endif
      enddo
      do i=1,nrans
        if(average_rans(i)==0.or.average_rans(i)==1) then
          irans_ave(i)=average_rans(i)
        else
          call FFRABORT(1,'ERR: average_rans(i) must bc 0 or 1')
        endif
      enddo
      if(rms_p==0.or.rms_p==1) then
        ip_rms=rms_p
      else
        call FFRABORT(1,'ERR: rms_p must bc 0 or 1')
      endif

      if(rms_Mut==0.or.rms_Mut==1) then
        imu_rms=rms_Mut
      else
        call FFRABORT(1,'ERR: rms_Mut must bc 0 or 1')
      endif

      if(rms_t==0.or.rms_t==1) then
        it_rms=rms_t
      else
        call FFRABORT(1,'ERR: rms_t must bc 0 or 1')
      endif
      do i=1,ncomp
      if(rms_comp(i)==0.or.rms_comp(i)==1) then
        icomp_rms(i)=rms_comp(i)
      else
        call FFRABORT(1,'ERR: rms_comp(i) must bc 0 or 1')
      endif
      enddo
      do i=1,nrans
      if(rms_rans(i)==0.or.rms_rans(i)==1) then
        irans_rms(i)=rms_rans(i)
      else
        call FFRABORT(1,'ERR: rms_rans(i) must bc 0 or 1')
      endif
      enddo
      if(re_uvw_t==0.or.re_uvw_t==1) then
        iuvwt_re=re_uvw_t
      else
        call FFRABORT(1,'ERR: re_uvw_t must bc 0 or 1')
      endif

      do i = 1, ncomp
         if(re_uvw_comp(i)==0.or.re_uvw_comp(i)==1) then
            iuvwc_re(i)=re_uvw_comp(i)
         else
            call FFRABORT(1,'ERR: re_uvw_t must bc 0 or 1')
         endif
      enddo


      if(maxmin_uvw_p==1.or.maxmin_uvw_p==0) then
        imaxmin=maxmin_uvw_p
      else
        call FFRABORT(1,'ERR: maxmin_uvw_p must bc 0 or 1')
      endif
      do i=1,nrans
        if(re_uvw_rans(i)==0.or.re_uvw_rans(i)==1) then
          iuvw_rans_re(i)=re_uvw_rans(i)
        else
          call FFRABORT(1,'ERR: re_uvw_rans(i) must bc 0 or 1')
        endif
      enddo
!
      if(iuvwt_re.eq.1.and.iuvw_ave_rms_re.eq.0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 
     &   ' ### MSG: reset [iuvw_ave_rms_re=1] while ',
     &  '[iuvwt_re=1] in ',cntlnam
        call FFRABORT(1,'module_les')
      endif
!

      do i=1,ncomp
         if(iuvwc_re(i)/=0.and.iuvw_ave_rms_re.eq.0) then
            write(ifle,*) '### error : data error'
            write(ifle,*) 
     &           '### reset [uvw_ave_rms_re=1], while ',
     &           '[re_uvw_rans(i)=1] in',cntlnam
            call FFRABORT(1,'module_les')
         endif
      enddo
!

      do i=1,nrans
      if(iuvw_rans_re(i)/=0.and.iuvw_ave_rms_re.eq.0) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 
     &    '### reset [uvw_ave_rms_re=1], while ',
     &  '[re_uvw_rans(i)=1] in',cntlnam
        call FFRABORT(1,'module_les')
      endif
      enddo
!
      if(surface_SDOT==0.or.surface_SDOT==1) then
        SDOT_suf=surface_SDOT
      else
        call FFRABORT(1,'ERR: surface_SDOT must bc 0 or 1')
      endif
      if(bulk_thickness==0.or.bulk_thickness==1) then
        blk_thick=bulk_thickness
      else
        call FFRABORT(1,'ERR: bulk_thickness must bc 0 or 1')
      endif

      if(depos_etch_speed==0.or.depos_etch_speed==1) then
        depoeth_spd=depos_etch_speed
      else
        call FFRABORT(1,'ERR: depos_etch_speed must bc 0 or 1')
      endif

      

      if(surface_mole_fraction==0.or.surface_mole_fraction==1) then
        molefr_suf=surface_mole_fraction
      else
        call FFRABORT(1,'ERR: surface_mole_fraction must bc 0 or 1')
      endif
      if(gas_mole_fraction==0) then
        mol_rslt=0
      elseif(gas_mole_fraction==1) then
        mol_rslt=1
        write(ifll,'(a)') 
     &   'MSG: result file will be output in mole fraction'
      elseif(gas_mole_fraction/=0.and.gas_mole_fraction/=1) then
        write(ifle,'(a)') 'ERR: gas_mole_fraction MUST be 0 or 1'
        call FFRABORT(1,'ERR: defining error for gas_mole_fraction')
      endif
      if(GAS_WDOT==0.or.GAS_WDOT==1) then
        WDOT_gas=GAS_WDOT
      else
        call FFRABORT(1,'ERR: GAS_WDOT must bc 0 or 1')
      endif
!
      if(nneq/=0) then
      do i=1,msuf_out
      if(RATE_EQ_SpecNAME(i)/=' ') then
        do j=1,80
        uder=RATE_EQ_SpecNAME(i)(j:j)
        ieq=0
        if(uder=='_') then
          cdum=RATE_EQ_SpecNAME(i)(1:j-1)
          read(cdum,*) ieq
          cdum=''
          cdum=RATE_EQ_SpecNAME(i)(j+1:80)
          exit
        endif
        enddo
        if(ieq<=0.or.ieq>nneq) then
          call FFRABORT
     &     (1,'ERR: Reading error in &les/RATE_EQ_SpecNAME(i)')
        endif
        is=0
        do s=1,ncompall
        if(trim(spinam(s))==trim(cdum)) then
          is=s
          exit
        endif
        enddo
        if(is==0) then
          write(ifle,'(2a)') 
     &    'ERR: Not found species name :' ,trim(cdum)
          call FFRABORT(1,'ERR: &les/RATE_EQ_SpecNAME(i)')
        endif
        dum=net_vreq(ncompall*(ieq-1)+is)
        if(abs(dum)<SML) then
          write(ifle,'(3a,I4)') 'ERR: not found species name: ',
     &    trim(cdum),' in eq.=',ieq
          call FFRABORT(1,'ERR: Defining error')
        else
          num_ratout=num_ratout+1
          if((num_ratout+1)>=msuf_out) then  !num_ratout<=10
            write(ifle,'(a)') 
     &      'ERR: in defining RATE_EQ_SpecNAME(i), i < or =10'
            call FFRABORT(1,'ERR: in &les/RATE_EQ_SpecNAME(i)')
          endif
          rat_E_S(i,is)=ieq  !*ncompall+is
        endif
        if(my_rank==0) then
          write(ifll,*)
          write(ifll,"(2x,108('+'))")
          write(ifll,'(2x,3a,I4,a)') 
     &   'MSG: Species production rate of ',
     &    trim(spinam(is)),' by reaction eq: ',
     &    ieq,' will be output'
        endif
      endif
      enddo
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_les
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_model 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(LEN=18),parameter,private :: modnam='(module_model)'
!
!
!< module for model flags >
!-------------------------
!-< 1. Flow model >-
!-------------------------
      integer,parameter :: incomp=0, mach0=1, comp=2, weak=3
      integer :: idrdp=mach0,ICVREF=1,ical_topt=1
      real(8) :: PREFM,expfact,Rh
      integer,parameter :: PEFC=1,PAFC=2,MCFC=3,SOFC=4,
     &                     AFC=5,MDFC=6
      integer :: IFUCL=0
      integer :: iheat=1
!-------------------------------
! --- canopy Drag Force coeff.
!-------------------------------
      real(8) :: drag_coef(21:50)
      
!
! incomp : incompressible
! mach0  : zero Mach no. approximation
! comp   : compressible
! idrdp  : flag for incompressible/compressible
!
!-------------------------
!-< 2. Turbulence model >-
!-------------------------
      integer,parameter :: noturb=0, ke=1,sles=2,dles=3, lles=4,
     &                     DNS=5,RSM=6,ke_low=7,RNG=8,CHEN=9,
     &                     SDES=10,cnst=11,KE2S=12,KLES=13
      integer :: icaltb=ke
!
      logical,save :: RANS_MDL=.false.
!
      integer,parameter :: vertex_cen=1,cell_cen=2
      integer :: icon_cv=vertex_cen
!
! noturb : no model
! ke     : k-epsilon model
! sles   : LES model,Smagorinsky LES model
! dles   : Dynamic LES model
! lles   : Lagrangian Dynamic LES model
! dns    : Direct Numerical Simulation (same as 'noturb')
! RSM    : Reynolds Stress Equation Model
!
! icaltb : flag for turbulence model
!
!
!-< 3. calculation flag>
      integer,save :: encomp=1
      logical,save :: ical_t=.true.,ical_reac=.true.,ical_surf=.true.
      logical,save :: ical_week=.false.
      logical,save :: ical_vect=.false.
!
      integer,save :: ical_mvmsh=0

!
      integer,save :: ical_MHD=0
      integer,save :: iLBF_P=0
      logical,save :: iHPC_close=.false.
      integer,save :: ical_dens=0
      integer,save :: ical_field=0
      integer,save :: ical_prt=0
      integer,save :: ical_WDR=0
      integer,save :: ical_thmoDIFF=0
      integer,save :: ical_tetra=0
!
      integer,save :: ical_EWT=0
!-----------------------
! --- user function
!-----------------------
      integer,parameter :: uf_mm=200
      integer,save :: u_func(uf_mm)
!---------------------------------------------------------------------------------------------------------
!     ical_prt=0: no particle calculation
!     ical_prt=1: Unsteady fluid field calculation + unsteady particle calculation
!     ical_prt=2: Firstly Steady fluid field calculation => unsteady particle calculation + heat transfor 
!     ical_prt=3: Firstly Steady fluid field calculation => Steady particle calculation (NOT finish)
!---------------------------------------------------------------------------------------------------------

!-< 4. openMP flag>
! --- openMP
! 
      integer :: nthrds=1   !NODE_CPU
!
!-< 5. openMP flag>
!
      integer,save :: monitor_stp
      
! encomp : =0 : hys_admin is not calculated
!          =1 : hys_admin is should be calculated
!          =2 : Firstly Steady fluid field calculation => 
!               then only hys is is calculated ,vel & prs are NOT calculated 
!///////////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=======================================================================
      subroutine inputdata
     &   (ifli,ifll,ifle,cntlnam,lcomp,l_temp,kemdl,RSMmdl,lowke,lKLES,
     &    l2SKE,lRNG,lCHEN,lSDES,RANSMDL,ivector,iMHD,
     &    iMvmsh,iBODY,lFC,ifld,
     &    iprtcle,ip_mdl,icomp,idens,ireac,ierror)
!=======================================================================
      implicit none
!
! 1. Input flags for model
!
! --- [dummy arguments]
!
      integer,intent(in)     :: ifli,ifll,ifle
      character(*),intent(in):: cntlnam
      logical,intent(out)    :: lcomp,kemdl,RSMmdl,lowke,lRNG,
     &                          lCHEN,l2SKE,lKLES,
     &                          lSDES,RANSMDL
      integer,intent(inout)  :: lFC,ifld,iprtcle,l_temp,ip_mdl,icomp
      integer,intent(out)    :: ierror,ivector,iMHD,iMvmsh,iBODY
      integer,intent(inout)  :: idens,ireac
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(14),parameter ::
     &   flwlst(3)=(
     &    /'incompressible','zero-Mach     ','compressible  '/)

      character( 4),parameter :: trblst(14)=(/
     & 'no  ','KE  ','SLES','DLES','LLES','DNS ','RSM ','LKE ',
     & 'RNG ','CHEN','SDES', 'cnst', '2SKE','KLES'/)

      character( 4),parameter :: fulclist(6)=(/
     & 'PEFC','PAFC','MCFC','SOFC','AFC ','MDFC'/)
      character(6),parameter ::
     &   CVlst(2)=(/'vertex','cell  '/)
! --- PEFC: 
! --- PAFC: 
! --- MCFC: 
! --- SOFC: 
! --- AFC : 
! --- MDFC:
      integer :: iset=0
      integer :: iflw,itrb,ios,ierr1,ierr2
!
! --- [namelist] 
!
      integer       :: cal_tys=1,cal_t=1,cal_reac=0,cal_surf=0,
     &                 cal_weak=0,cal_MHD=0,cal_vect=0,cal_mvmsh=0
      integer       :: NODE_CPU=1,moni_inter=1,LBF_P=0,HPC_close=0
      integer       :: density=0,iFIELD=0,cal_heat_reaction=1
      integer       :: cal_particle=0
      integer       :: cal_tet,cal_thermo_diff
      character(LEN=20) :: flow,trbmdl
      character(LEN=20) :: Fuel_cell
      real*8        :: Buffle_shift=0.d0
      integer       :: EWT=0  !Enhanced Wall Treatment
      character(LEN=6) :: CV
!---------------------------------------
! --- model namelist
!---------------------------------------
      namelist /model/ flow,trbmdl,
     &                 cal_tys,cal_t,cal_reac,cal_surf,cal_weak,
     &                 cal_MHD,cal_vect,NODE_CPU,moni_inter,cal_mvmsh,
     &                 LBF_P,HPC_close,density,
     &                 Fuel_cell,iFIELD,cal_heat_reaction,
     &                 cal_particle,
     &                 cal_tet,cal_thermo_diff,Buffle_shift,
     &                 EWT,CV
!
      integer       :: userflag(1:uf_mm)=0,iu
      namelist /user_function/ userflag
! --- 
!
      integer           :: rain_WDR=0
      integer,parameter :: iundef=-huge(1)
      integer,parameter :: M_injctr=20
      real*8            :: Rh_mm_Hour
!---------------------------------------
! --- rain namelist
!---------------------------------------
      namelist /rain/  rain_WDR,Rh_mm_Hour
!
!
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      flow  =' '
      CV=' '
      trbmdl=' '
      Fuel_cell=' '
      cal_tys=1
      cal_t=1
      cal_reac=0
      cal_surf=0
      cal_weak=0
      cal_MHD=0
      cal_vect=0
      cal_mvmsh=0
      moni_inter=1
      LBF_P=0
      HPC_close=0
      iFIELD=0
      cal_heat_reaction=1
      cal_particle=0
      cal_tet=0
      cal_thermo_diff=0
!
      EWT=0
!
      rain_WDR=0
      Rh_mm_Hour=iundef
!
      rewind ifli
      read(ifli,model,iostat=ios)
!      write(ifll,model)
      call nml_errmsg0(ifle,ios,'model',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
!---------------------------
! --- 
!---------------------------
      rewind ifli
      read(ifli,rain,iostat=ios)
!      call nml_errmsg0(ifle,ios,'rain',ierr1)
      if( ios>0 ) goto 9999
      if(rain_WDR==1) then
        if(Rh_mm_Hour==iundef) then
          write(ifle,'(1X,a)') 
     &  'MSG: rain_WDR=1, Rh_mm_Hour MUST be defined in mm/h'
          call FFRABORT(1,'ERR: reading error in &rain')
        else
          if(Rh_mm_Hour<0.d0) then
             call FFRABORT(1,'ERR: Rh_mm_Hour<0.d0 in &rain')
          else
             Rh=Rh_mm_Hour
          endif
        endif
      endif

!---------------------------
! --- 
!---------------------------
      userflag(:)=0
      u_func(:)=0
      rewind ifli
      read(ifli,user_function,iostat=ios)
      if( ios>0 ) goto 9999
      do iu=1,uf_mm
      if(userflag(iu)/=1.and.userflag(iu)/=0) 
     &   call FFRABORT(1,'ERR-MSG: userflag(:) MUST be 1 or 0')
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
        write(ifle,'(1x,a,I4)') 'ERR: &user_function/userflag(i), i=',iu
        call FFRABORT(1,'ERR: at &user_function')
      elseif(userflag(iu)/=1.and.userflag(iu)/=0) then
        call FFRABORT
     &  (1,'ERR-MSG: at &user_function/userflag(1:200) MUST be 0 or 1')
      else
        
      endif
      enddo
!---------------------------
!---------------------------
!---------------------------
      call nml_listno(3,flwlst,flow  ,idrdp)
      call nml_listno(14,trblst,trbmdl,icaltb)
      call nml_chkchr0(ifle,'flow'  ,flow  ,idrdp ,ierr1)
      call nml_chkchr0(ifle,'trbmdl',trbmdl,icaltb,ierr2)
      if(ierr1.ne.0.or.ierr2.ne.0) goto 9999
      call nml_listno(2,CVlst,CV,icon_cv)
      call nml_chkchr0(ifle,'CV'  ,CV,icon_cv ,ierr1)
      if(ierr1.ne.0) goto 9999
!
      if(Fuel_cell/=' ') then
        call nml_listno(6,fulclist,Fuel_cell,IFUCL)
        call nml_chkchr0(ifle,'Fuel_cell',Fuel_cell,IFUCL,ierr1)
        lFC=IFUCL
! --- 
         write(ifle,'(1x,a)') 
     &   "MSG: Fuel_cell= for NOT fuel-cell module"
        if(ierr1/=0) goto 9999
      endif
!
      idrdp=idrdp-1
      icaltb=icaltb-1
!
                          lcomp=.true.
      if(idrdp.eq.incomp) lcomp=.false.
!
!      if(icaltb.eq.ke.or.icaltb.eq.RSM) then
!        rns_scl=.true.
!      else
!        rns_scl=.false.
!      endif
!
      kemdl=.false.
      RSMmdl=.false.
      lowke=.false.
      lRNG=.false.
      lCHEN=.false.
      lSDES=.false.
      l2SKE=.false.
      lKLES=.false.
      RANSMDL=.false.
      if(icaltb.eq.ke) then
        kemdl=.true.
        RANSMDL=.true.
      elseif(icaltb.eq.RSM) then
        RSMmdl=.true.
        RANSMDL=.true.
      elseif(icaltb.eq.ke_low) then
        lowke=.true.
        RANSMDL=.true.
      elseif(icaltb.eq.RNG) then
        lRNG=.true.
        RANSMDL=.true.
      elseif(icaltb.eq.CHEN) then
        lCHEN=.true.
      elseif(icaltb.eq.KE2S)then 
        l2SKE=.true.
      elseif(icaltb.eq.SDES) then
        lSDES=.true.
      elseif(icaltb.eq.KLES) then
        lKLES=.true.
      endif
!
      if(u_func(2)==1) then
        if(icaltb/=KLES) then
          call FFRABORT(1,'ERR: u_func(2)=1 Only support trbmdl="KLES" ')
        endif
      endif
!
      if(icaltb==ke.or.icaltb==RSM.or.icaltb==ke_low.or.
     &   icaltb==RNG.or.icaltb==CHEN.or.icaltb==KE2S) then
        RANS_MDL=.true.
      else
        RANS_MDL=.false.
      endif
!
      if(userflag(1)==1.and..NOT.RANS_MDL) then
        call FFRABORT(1,'ERR: userflag(1)==1 ONLY support RANS model')
      endif
!
      if(EWT==1) then
        ical_EWT=1
        if(icaltb==RNG.or.icaltb==CHEN.or.
     &     icaltb==ke.or.icaltb==KE2S) then
        elseif(icaltb==ke_low) then
           call 
     &   FFRABORT(1,'ERR: EWT=1 NOT suible for Low-Re model')
        else
          call FFRABORT(1,'ERR: EWT=1 ONLY for RNG,CHEN or KE')
        endif
      else
        ical_EWT=0
      endif
!
      if(.not.(cal_tys==2.or.cal_tys==1.or.cal_tys==0))then
        write(ifle,*) '    ### EER: [cal_tys] MUST be 1, 0  or 2'
        call FFRABORT(1,'module_model')
      endif
      if(.not.(cal_t==1.or.cal_t==0.or.cal_t==2))then
        write(ifle,*) '    ### EER: [cal_t] MUST be 1 ,2 or 0 '
        call FFRABORT(1,'module_model')
      endif
      if(.not.(cal_reac==1.or.cal_reac.eq.0))then
        write(ifle,*) '    ### EER: [cal_reac] MUST be 1 or 0 '
        call FFRABORT(1,'module_model')
      endif
      if(.not.(LBF_P==1.or.LBF_P==0.or.
     &         LBF_P==2.or.LBF_P==3.or.LBF_P==4.or.LBF_P==5))then
        write(ifle,*) '    ### EER: [LBF_P] MUST be 0,1,2,3,4,5'
        call FFRABORT(1,'module_model')
      endif
      if(.not.(density==0.or.density==1.or.density==4.or.
     &         density==2.or.density==3.or.density==5
     &  ))then
        write(ifle,'(1x,a)') 'ERR: [Density] MUST be 0,1,2,3,4,5'
        write(ifle,'(1x,a)') 
     &    'MSG: Density=0: rho=const for [incompressible] flow'
        write(ifle,'(1x,2a)') 
     &    'MSG: Density=0: ',
     &     'rho=p/R/T/(SUM(Yi/Mi)) for [compressible] or [zero] flow'
        write(ifle,'(1x,2a)') 
     &    'MSG: Density=1: rho=[SUM(Yi/rhoi)]^-1 ',
     &    'for multi-component incompressible flow'
        write(ifle,'(1x,2a)') 
     &    'MSG: Density=2: rho=rho0/(1+beta*(T-T0)) ',
     &     'for tempture variation incompressible flow => NOT FINISHED'
        write(ifle,'(1x,a,/,17X,2a)') 
     &    'MSG: Density=3: rho=rho0/(1+beta*(T-T0)+SUM(beta_i*Yi)) ',
     &     'for tempture variation multi-component ',
     &      'incompressible flow => NOT FINISHED'
        write(ifle,'(1x,a,a)') 
     &    'MSG: Density=4: Redlich-Kwong state equation for ',
     &      'compressible flow'
        call FFRABORT(1,'module_model')
      endif
      if(.not.(HPC_close==1.or.HPC_close.eq.0))then
        write(ifle,*) '    ### EER: [HPC_close] MUST be 1 or 0 '
        call FFRABORT(1,'module_model')
      endif
      if(.not.(cal_surf==1.or.cal_surf.eq.0))then
        write(ifle,*) '    ### EER: [cal_surf] MUST be 1 or 0 '
        call FFRABORT(1,'module_model')
      endif
      if(.not.(cal_weak==1.or.cal_weak.eq.0)) then
        write(ifle,*) '    ### EER: [cal_weak] MUST be 1 or 0 '
        call FFRABORT(1,'module_model')
      endif
      if((idrdp==incomp.and.cal_t==2)
     &  .or.(idrdp==mach0.and.cal_t==2) 
     &  )then
        write(ifle,*) '    ### EER: [cal_t==2] error use ',
     &   'in incompressible flow'
        call FFRABORT(1,'module_model')
      endif
!
      if(idrdp/=incomp.and.cal_weak==1) then
        write(ifle,*) 
     &   '  ### EER: [cal_weak] MUST be set for [incomp] '
        call FFRABORT(1,'module_model')
      endif
!
      if(.not.(cal_vect==1.or.cal_vect==0))then
        write(ifle,*) '    cal_vect=1 : vector computer '
        write(ifle,*) '    cal_vect=0 : Scalar computer '
        write(ifle,*) '    ### EER: [cal_vect] MUST be 1 or 0 '
        call FFRABORT(1,'module_model')
      endif
!
      if(.not.(cal_mvmsh==1.or.
     &         cal_mvmsh==0.or.
     &         cal_mvmsh==2.or.
     &         cal_mvmsh==3.or.
     &         cal_mvmsh==4.or.
     &         cal_mvmsh==5
     &         ))then

        write(ifle,'(6x,a)') 'EER: [cal_mvmsh] MUST be 0, 1 2 3 4 or 5'
!
        write(ifle,'(6x,2a)') 
     &   'MSG: cal_mvmsh=0 : ',
     &   'No [Moving] & No [Removing/Adding] Mesh solver'
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
        call FFRABORT(1,'module_model')
      endif
!
      if(.not.(cal_MHD==1.or.cal_MHD==0.or.cal_MHD==2))then
        write(ifle,'(1x,a)') 
     &   '    cal_MHD=1 : MHD (MagnetoHydroDynamics) Coupling '
        write(ifle,'(1x,a)') 
     &   '    cal_MHD=2 : MHD Only (No fluid)'
        write(ifle,'(1x,a)') 
     &   '    cal_MHD=0 : No calculation of MHD, Fluid only'
        write(ifle,'(1x,a)') 
     &   '    ### EER: [cal_MHD] MUST be 0 ,1 or 2'
        call FFRABORT(1,'module_model')
      endif
!
      if(NODE_CPU<=0) then
        write(ifle,'(1x,a)') 
     &   '    NODE_CPU >0 : Node-cpu (PE) number for openMP'
        call FFRABORT(1,'module_model')
      endif
!
      if(IFUCL/=0) then
        write(ifll,'(1X,2a)') 'MSG: FULL-CELL Module: ',fulclist(IFUCL)        
      endif
!
      encomp=cal_tys 
      icomp=cal_tys 
!
      ical_t=cal_t==1   !zhang-cvd
      ical_topt=cal_t   !cal_t=1 : enthalpy-eq; cal_t=2 : dannetu 
!
      ical_reac=cal_reac==1
      ical_surf=cal_surf==1
      ical_week=cal_weak==1
      ical_MHD=cal_MHD
      iMHD=cal_MHD
      ical_vect=cal_vect==1
      ivector=cal_vect
      ical_mvmsh=cal_mvmsh
      iMvmsh=cal_mvmsh
      iLBF_P=LBF_P
      iBODY=LBF_P
      iHPC_close=HPC_close==1
      ical_dens=density
      if(cal_tys>=1.and.idrdp==incomp.and.density==0) then
        write(ifle,'(1x,a)') 
     &    'MSG: Density=0: rho=const for [incompressible] flow'
        write(ifle,'(1x,2a)') 
     &    'MSG: Density=0: ',
     &     'rho=p/R/T/(SUM(Yi/Mi)) for [compressible] or [zero] flow'
        write(ifle,'(1x,2a)') 
     &    'MSG: Density=1: rho=[SUM(Yi/rhoi)]^-1 using &initial/p ',
     &    'for multi-component incompressible flow'
        write(ifle,'(1x,2a)') 
     &    'MSG: Density=1 is reset in inner fflow ',
     &    'for multi-component incompressible flow'
        ical_dens=1
      endif
      if(ical_dens==1.or.ical_dens==2.or.ical_dens==3) then
        if(idrdp/=incomp) then
          write(ifle,'(1x,a)') 
     &    'MSG: Density=0: rho=const for [incompressible] flow'
          write(ifle,'(1x,2a)') 
     &    'MSG: Density=0: ',
     &     'rho=p/R/T/(SUM(Yi/Mi)) for [compressible] or [zero] flow'
          write(ifle,'(1x,2a)') 
     &    'MSG: Density=1: rho=[SUM(Yi/rhoi)]^-1 ',
     &    'for multi-component incompressible flow'
          write(ifle,'(1x,2a)') 
     &    'MSG: Density=2: rho=rho0/(1+beta*(T-T0)) ',
     &     'for tempture variation incompressible flow => NOT FINISHED'
          write(ifle,'(1x,3a)') 
     &    'MSG: Density=3: rho=rho0/(1+beta*(T-T0)+SUM(beta_i*Yi)) ',
     &     'for tempture variation multi-component ',
     &      'incompressible flow => NOT FINISHED'          
          call FFRABORT(1,'ERR: ')
        endif
      elseif(ical_dens==4) then
        idens=4
        if(idrdp/=comp) then
          write(ifle,'(1x,3a)') 
     &    'MSG: Density=4: Redlich-Kwong state equation for ',
     &      'compressible flow'
          call FFRABORT(1,'ERR: set [compressible] in &model ')
        else
          write(ifll,'(1x,2a)') 
     &    'MSG: Density=4: Redlich-Kwong state equation for ',
     &      'compressible flow'
        endif
      elseif(ical_dens==5) then
        
      endif
!
      nthrds=NODE_CPU
!
      if(ical_t) then
        if(encomp==0) encomp=1
      endif
!
      if(ical_reac) then
        ireac=1
      else
        ireac=0
      endif
!
      if(ical_reac.or.ical_surf) then
        if(encomp==0) encomp=1
        ical_t=.true.
      endif
!
      if(ical_t) then
        l_temp=1
      endif
!
      if(moni_inter<=0) then
        moni_inter=1
      endif
      monitor_stp=moni_inter
!
      if(iFIELD>0) then
        ifld=iFIELD
        ical_field=iFIELD
      else 
        ifld=0 
        ical_field=0 
      endif 
! 
      if(cal_heat_reaction==1) then 
        iheat=1
      else
        iheat=0
      endif
!
      if(.NOT.(rain_WDR==1.or.rain_WDR==0)) then
        call FFRABORT
     &    (1,'ERR-MSG: MUST set &model/rain_WDR=0 or rain_WDR=1')
      elseif(rain_WDR==1) then
        if(cal_particle==0) then
          call FFRABORT
     &    (1,'ERR-MSG: rain_WDR==1 MUST set &model/cal_particle=1')
        endif
      endif
!
!  if(rain_WDR==0) iprtcle=1
!  if(rain_WDR==1) iprtcle=2
!
      if(cal_particle==1.or.cal_particle==2.or.cal_particle==3) then 
        iprtcle=1
        ical_prt=cal_particle
        ip_mdl=cal_particle
        ical_WDR=rain_WDR
        if(ical_WDR==1) then
          iprtcle=2
        endif
      elseif(cal_particle==0)then 
        iprtcle=0
        ical_prt=0
        ip_mdl=0
        ical_WDR=0
        if(ical_WDR==1) then
          write(ifle,'(1x,a)')
     &    'setting up WDR model MUST set [&model/cal_particle>=1]'
          call FFRABORT
     &   (1,'ERR: stop at &model') 
        endif
      else
        call FFRABORT
     &   (1,'ERR: cal_particle=1,2,3 for partical TRACING')
      endif
!
      if(cal_thermo_diff==1) then 
        ical_thmoDIFF=1
      else
        ical_thmoDIFF=0
      endif
!
      if(cal_tet==1) then 
        expfact=0.d0
        ical_tetra=1
      else
        expfact=1.d0
        ical_tetra=0
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_model
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_vector
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(13),parameter,private :: modnam='(module_vector)'
!
!-----------------------------
!-< 1. Vector Computation >-
!-----------------------------
!
      integer,save :: ICVS_V,ICVE_V
      integer,save :: IDCS_V,IDCE_V
      integer,save :: ICVSIN_V,ICVEIN_V
      integer,save :: ICFS_V,ICFE_V
      integer,save,allocatable :: index_c(:)
      integer,save,allocatable :: index_f(:)
      integer,save :: MXUP,MXLW
      integer,save :: N2=128
      integer,save :: NUmax,NLmax
      integer,save :: NCOLOR=-1,NO
!
      integer,save :: IEMAX_V
!///////////////////////////////////////////////////////////////////////
      contains
!=======================================================================
      subroutine inputdata(icall,
     &     ifli,ifll,ifle,lvect,nthrds,iLBF_P,ical_prt,
     &     NCV,NCVIN,NALLCV,NCVFAC,NMAT,NSSFBC,MXMAT,MAXIE,
     &     MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,
     &     vctr,vctrp,NEIGHBOR,DNORM,SFAREA,SFCENT,
     &     MXCV_V,MXBND_V,MXCV_VP,MXBND_VP,
     &     IQ,MXALLCV,LVEDGE,MXCVFAC,imode,
     &     MEP,MXALLCV_P,MXCVFAC_P,
     &     ierror)
!=======================================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifli,ifll,ifle,MXMAT,nthrds,imode,
     &                           iLBF_P,ical_prt,MAXIE,icall
      integer,intent(in)      :: NCV,NCVIN,NALLCV,NCVFAC,NMAT,NSSFBC,
     &                           MEP,MXALLCV_P,MXCVFAC_P
      integer,intent(in)      :: MAT_INDEX(0:MXMAT)
      integer,intent(in)      :: MAT_CVEXT(0:MXMAT)
      integer,intent(in)      :: MAT_DCIDX(0:MXMAT)
      integer,intent(in)      :: MAT_CFIDX(0:MXMAT)
      integer,intent(out)     :: ierror
      logical,intent(in)      :: lvect
      integer,intent(in)      :: MXCV_V,MXBND_V,MXCV_VP,MXBND_VP,
     &                           MXALLCV,MXCVFAC
      integer,intent(in)      :: LVEDGE(2,MXCVFAC)

      integer,intent(inout)   :: vctr(MXCV_V,0:MXBND_V)
      integer,intent(inout)   :: vctrp(MXCV_VP,0:MXBND_VP)
      integer,intent(out)     :: NEIGHBOR(MEP,MXALLCV_P)
      real*8 ,intent(out)     :: DNORM   (MXCVFAC_P)

      real*8 ,intent(in)      :: SFAREA(4,  MXCVFAC)
      real*8 ,intent(in)      :: SFCENT(3,  MXCVFAC)
      integer,intent(inout)   :: IQ   (    MXALLCV)
!
      integer :: ICFL,ICVA,ICVB,ICVL,i,kl,ku,j,
     &           nnthrds=1,IIMAT,ICFS,ICFE,ICVS,ICVE,
     &           ICMAX,ICMIN,IB,IA
      REAL*8  :: dum1
      LOGICAL :: LDUMMYA, LDUMMYB
!-----------------------------------------------------
!      integer,external :: omp_get_num_threads,
!     &                    omp_get_max_threads,
!     &                    omp_get_thread_num,
!     &                    omp_get_num_procs
!      logical,external :: omp_get_dynamic,
!     &                    omp_in_parallel,
!     &                    omp_get_nested,
!     &                    omp_test_lock
!-----------------------------------------------------      
      ierror=0
      ICVS_V=1
      ICVE_V=1
      IDCS_V=1
      IDCE_V=1
      ICFS_V=1
      ICFE_V=1
      ICVSIN_V=1
      ICVEIN_V=1
!---------------------------
! --- 
!---------------------------
      if(iLBF_P==3.or.lvect) then
        IQ(:)=0
        do 100 IIMAT=1,NMAT 
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        IQ(ICVA)=IQ(ICVA)+1
        vctr(ICVA,IQ(ICVA))=-ICFL
        IQ(ICVB)=IQ(ICVB)+1
        vctr(ICVB,IQ(ICVB))= ICFL
        enddo
 100    enddo
        do 200 IIMAT=1,NMAT 
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        vctr(ICVL,0)=IQ(ICVL)
        enddo
 200    enddo
      endif
      if(ical_prt>=1) then
        do IIMAT=1,NMAT 
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        ICMAX=MAT_CVEXT(NMAT)
        ICMIN=MAT_CVEXT(NMAT-1)+1

        do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        IA=VCTRP(ICVA,0)+1
        IF(IA.LE.MEP) THEN
          vctrp(ICVA,IA)=-ICFL
          VCTRP(ICVA,0)=IA
          LDUMMYA=.TRUE.
          IF(ICVB.GE.ICMIN.AND.ICVB.LE.ICMAX) THEN 
            NEIGHBOR(IA,ICVA)=ICVB
          else
            NEIGHBOR(IA,ICVA)=-ICVB
          endif

        else
          call FFRABORT(1,'ERR: vctrp-1')
        endif
!
        IB=VCTRP(ICVB,0)+1
        IF(IB.LE.MEP) THEN
          IF(ICVB.GE.ICMIN.AND.ICVB.LE.ICMAX) THEN
            vctrp(ICVB,IB)=ICFL
            VCTRP(ICVB,0)=IB
            LDUMMYB=.TRUE.
            IF(ICVA.GE.ICMIN.AND.ICVA.LE.ICMAX) THEN 
              NEIGHBOR(IB,ICVB)=ICVA
            else
              NEIGHBOR(IB,ICVB)=-ICVA
            endif
          else
            LDUMMYB=.FALSE.
          endif
        else
          call FFRABORT(1,'ERR: vctrp-2')
        endif
        DNORM(ICFL)=SFCENT(1,ICFL)*SFAREA(1,ICFL)
     *             +SFCENT(2,ICFL)*SFAREA(2,ICFL)
     *             +SFCENT(3,ICFL)*SFAREA(3,ICFL)

        enddo
        enddo
      endif
!-------------------
! --- 
!-------------------
      if(lvect) then
!--------------------------------------------------------
        ICVS_V=1
        ICVE_V=NCV       !MAT_CVEXT(NMAT)
        IDCS_V=MAT_DCIDX(0)+1
        IDCE_V=MAT_DCIDX(NMAT)
        ICFS_V=MAT_CFIDX(0)+1
        ICFE_V=MAT_CFIDX(NMAT)
        ICVSIN_V=1
        ICVEIN_V=NCVIN   !MAT_INDEX(NMAT)
!
        DO ICVL=ICVS_V,ICVE_V
          ku=0
          kl=0
          do j=1,MXBND_V
          if(vctr(ICVL,j)>0) then
            ku=ku+1
            MXUP=max(MXUP,ku)
          elseif(vctr(ICVL,j)<0) then
            kl=kl+1
            MXLW=max(MXLW,kl)
          endif
          enddo
        enddo
!--------------------------------------------------------
!!$      call omp_set_dynamic(.false.)
!!$      call omp_set_num_threads(nthrds)
!!$omp   parallel default(shared)
!!$      nnthrds=omp_get_num_threads()
        if(nnthrds/=nthrds) then
          write(ifle,'(1x,a,2I4)') 
     &    'MSG: Reset NODE_CPU, NODE_CPU/num_threads=',
     &     nthrds,nnthrds
          write(ifle,'(1x,2a)') 
     &    'MSG: if using openMP parallel (NODE_CPU>1), ',
     &    'Recompile all source using openMP option'
          call FFRABORT(1,'ERR: &model:NODE_CPU /= num_threads')
        endif
!!$omp   end parallel 
!
        if(imode==1.and.icall==1) 
     &    ALLOCATE(index_c(0:nthrds),index_f(0:nthrds))
        index_c(1:nthrds)=int(ICVE_V/(nthrds))
        index_c(0)=0
        do i=1,nthrds
        index_c(i)=index_c(i-1)+index_c(i)
        enddo
        index_c(nthrds)=ICVE_V
!
        index_f(1:nthrds)=int(ICFE_V/(nthrds))
        index_f(0)=0
        do i=1,nthrds
        index_f(i)=index_f(i-1)+index_f(i)
        enddo
        index_f(nthrds)=ICFE_V
!--------------------------------------------------------
      elseif(iLBF_P==3) then
      elseif(ical_prt==1) then
      endif
!
!
      end subroutine inputdata
      end module module_vector
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_movegrid
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(21),parameter,private :: modnam='(module_movegrid)'
!
!
!< module for moving grid >
!
      integer,parameter,private :: mgrid=10
!
      integer       :: ngrid=0
      real(8)       :: dtmin=0.d0
      real(8),save  :: time(mgrid+1)
      integer,save  :: ical_mov=0
!piston moving mesh
      real(8),save  :: ratio,rps,LINER,the0,RRR,omeini,the_end,endpiston
      integer,save  :: pisBC,ipiston=0,mv_dir=3
!
! ngrid : no. of grid tables
! dtmin : minimum time interval of grid tables
! time  : time of each grid
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,lrot,ierror)
!=================================================
      implicit none
!
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      logical,intent(in)  :: lrot
      integer,intent(out) :: ierror
!
! --- [namelist]
!
      real(8) :: cycle
      integer :: cylinder,piston_top_BC,direction=3
      real(8) :: engine_rpm,connect_rod,crank_R,initial_omega,
     &                         liner_BDC,initial_CA,end_CA
      namelist /movegrid/ ngrid,time,cycle
      namelist /piston_moving/ engine_rpm,connect_rod,crank_R,
     &                         liner_BDC,initial_CA,piston_top_BC,
     &                         cylinder,direction,initial_omega,
     &                         end_CA
!
! --- [local entities]
!
      integer,parameter :: iundef=-huge(1)
      real   ,parameter :: undef=-huge(1.)
      character(LEN=11),parameter :: subnam='(inputdata)'
      integer :: iset=0
      integer :: i,inpt,ios,ierr1
      real(8) :: dum1,pi=3.1415926d0
!
      engine_rpm=undef
      connect_rod=undef
      crank_R=undef
      liner_BDC=undef
      initial_CA=undef
      piston_top_BC=iundef
      cylinder=0
      initial_omega=undef
      end_CA=undef

      rewind ifli
      read(ifli,piston_moving,iostat=ios)
      if(ios.lt.0) goto 101
      call nml_errmsg0(ifle,ios,'piston_moving',ierr1)
      if( ierr1/=0 ) goto 9999

      if(cylinder/=0.and.cylinder/=1) then
        call FFRABORT(1,'ERR: cylinder MUST be 1 or 0')
      endif
      if(engine_rpm==undef) then
        call FFRABORT(1,'ERR: engine_rpm MUST be defined')
      endif
      if(connect_rod==undef) then
        call FFRABORT(1,'ERR: connect_rod MUST be defined')
      endif
      if(crank_R==undef) then
        call FFRABORT(1,'ERR: crank_R MUST be defined')
      endif
      if(liner_BDC==undef) then
        call FFRABORT(1,'ERR: liner_BDC MUST be defined')
      endif
      if(initial_CA==undef) then
        call FFRABORT(1,'ERR: initial_CA MUST be defined')
      endif
      if(piston_top_BC==iundef) then
        call FFRABORT(1,'ERR: piston_top_BC MUST be defined')
      endif
      if(end_CA==undef) then
        call FFRABORT(1,'ERR: end_CA MUST be defined')
      elseif(end_CA < initial_CA) then
        call FFRABORT(1,'ERR: end_CA MUST be greater than initial_CA')
      endif
      if(initial_omega==undef) then 
        call FFRABORT(1,'ERR: initial_omega MUST be defined')
      endif
      mv_dir=direction
      if(mv_dir/=3.and.mv_dir/=2.and.mv_dir/=1) then
        call FFRABORT
     & (1,'ERR: mv_dir=1,2 or 3 for defining move direction:x,y or z')
      endif
      ratio=crank_R/(connect_rod+1.d-10)
      rps=engine_rpm/60.d0
      LINER=liner_BDC
      the0=initial_CA*pi/180.d0
      the_end=end_CA*pi/180.d0
      ipiston=cylinder
      RRR=crank_R
      pisBC=piston_top_BC
      omeini=initial_omega
      endpiston=(the_end-the0)/(2.d0*pi*rps+1.d-10)
!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
 101  continue
      ngrid= 1
      time = 0.d0
      cycle=-1.d0
      inpt = 1
!
      rewind ifli
      read(ifli,movegrid,iostat=ios)
      if( ios.lt.0 ) goto 100
      call nml_chksiz(ifle,'mgrid',mgrid,ngrid,modnam,ierr1)
      if( ierr1.ne.0 ) goto 9999
      call nml_errmsg0(ifle,ios,'movegrid',ierr1)
      if( ierr1.ne.0 ) goto 9999
      inpt=2

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
      ical_mov=0
      if(ngrid>1) then
        ical_mov=1
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

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_nowtime
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!< module for present time & iteration counts >
!
      integer :: iter=0
      real(8)  :: time=0.d0
!
      end module module_nowtime

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_output
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(LEN=15),parameter,private :: modnam='(module_output)'
!
!------------------------------------------
! 1. module for control data to output
!------------------------------------------
      private
      public :: outchk,inputdata
      public :: tps,mlt_rslt,fforceFiletype,BC_output,
     &          ProbeFiletype
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
!
! --- fnmlst : list of file names
!
!
      integer,parameter,public :: none=-1, yes=1, no=0
!
! none : no output
! yes  : it is time to output
! no   : it is not time to output
!
      integer :: i
      integer ,save,public :: kopt(mfile)
      real(8) ,save,public :: start_force
      real(8) ,save,public :: start_cdcl,start_prob,start_fforce
      integer ,save        :: nspc(mfile)
      real(8) ,save        :: tspc(mspec,mfile)
      real(8) ,save        :: tps(mfile),tpe(mfile),tpd(mfile)
      logical              :: set(mfile)=(/(.false.,i=1,mfile)/)
      integer ,save        :: mlt_rslt=0
      integer ,save :: fforceFiletype,BC_output,ProbeFiletype
!zhang0119
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
      implicit none
      character(*),intent(in)  :: fnam
      integer     ,intent(in)  :: iter
      real(8)     ,intent(in)  :: time,delt
      integer   ,intent(inout) :: out
      integer :: i,l,ns,nd,idum,knum
      real(8)  :: eps
      character(LEN=8),parameter :: subnam='(outchk)'
      logical :: first(mfile)=(/(.true.,i=1,mfile)/)
!
! ---
!
      call modutl_chkset('>',1,iset,modnam//subnam)
!
      knum=out
      out=none
!
      call nml_listno(mfile,fnmlst,fnam,l)
      if( l.lt.1 ) then
        write(*,*) '### program error -1- ',modnam,subnam
        write(*,*) 'fnam = ',fnam
        CALL FFRABORT(1,'module_output/outchk')
      endif
!
      if( .not.set(l) ) return
!
      if( kopt(l).gt.2 ) then
        eps=0.1d0*delt
        if( eps.le.0.d0 ) then
          write(*,*) '### program error -2- ',modnam,subnam
          write(*,*) 'eps = ',eps
          CALL FFRABORT(1,'module_output/outchk')
        endif
      endif
!
      out=no
!
! --- /interval of iteration counts/
!
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
        if(knum.eq.100) then
          if(mod(iter,nint(tpd(l))).eq.0) then
            out=yes
            return
          else
            return
          endif
        else
          if( iter.lt.nint(tps(l)) ) return
          tps(l)=nint(tps(l)+tpd(l))
        endif
!
! --- /specified iteration counts/
!
      elseif( kopt(l).eq.2 ) then
        do 100 i=1,nspc(l)
        if( iter.eq.nint(tspc(i,l)) ) goto 101
  100   continue
        return
  101   continue
!
! --- /interval of time/
!
      elseif( kopt(l).eq.3 ) then
        if( time.gt.tpe(l)+eps ) return
        if( first(l) ) then
          tpd(l)=max(0.d0,tpd(l))
          if( time.gt.tps(l) .and. tpd(l).gt.0.d0 ) then
            tps(l)=tps(l)+dble(int(time-tps(l))/tpd(l))*tpd(l)
            if( time.gt.tps(l) ) tps(l)=tps(l)+tpd(l)
          endif
          first(l)=.false.
        endif
        if( time.lt.tps(l)-eps ) return
        tps(l)=tps(l)+tpd(l)
!
! --- /specified time/
!
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
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
! --- [namelist]
!
      character(LEN=20) :: file
      character(LEN=10) :: type
      integer           :: nspec,multi_result=0,BC_result=0
      character(20)     :: filetype
      real(8)  :: start,end,inter
      real(8)  :: spec(mspec)
      
      namelist /output/ file,type,start,end,inter,nspec,spec,
     &                  multi_result,filetype,BC_result
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      integer,parameter :: iundef=-huge(1)
      real(8) ,parameter :: undef=-huge(1.d0)
      character(LEN=7),parameter :: typlst(4)=(/
     1                              'inter_i',
     2                              'spec_i ',
     3                              'inter_t',
     4                              'spec_t '/)
      character(LEN=5),parameter :: cfile(2)=(/'GF   ','Excel'/)
!      
      integer :: lchk(mfile)
      integer :: i,nseq,nfl,ntp
      integer :: ios,ierr1,ierr2
      integer :: kfforce=0
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
      multi_result=0
      filetype=' '
      BC_result=0
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
        if(nfl.eq.4) start_force=start
        if(nfl.eq.6) start_cdcl=start
        if(nfl.eq.7) start_prob=start
        if(nfl.eq.8) start_fforce=start
!        start_force=start
      else
        nspc(nfl)=nspec
        do 200 i=1,nspec
        tspc(i,nfl)=spec(i)
  200   continue
      endif
!
      if(multi_result==0) then
      elseif(multi_result==1) then
        mlt_rslt=yes
        write(ifll,'(2a)') 
     &   'MSG: result file will be output in form of '
     &          ,'.frontflow_****'
      elseif(multi_result/=0.and.multi_result/=1) then
        write(ifle,'(a)') 'ERR: multi_result MUST be 0 or 1'
        call FFRABORT(1,'ERR: defining error for multi_result')
      endif
      if(multi_result==1.and.nfl/=2) then
        call FFRABORT
     & (1,'ERR: [multi_result] only defined for [result] file')
      endif
!-----------------------
! --- fforce file type
!-----------------------
      if(nfl==8) then
        call nml_listno(2,cfile,filetype,kfforce)
        if(kfforce==2) then
          fforceFiletype=2
        else
          fforceFiletype=1
        endif
        if(trim(filetype)/=trim(cfile(1)).and.
     &     trim(filetype)/=trim(cfile(2)).and.
     &     trim(filetype)/=' ') then
          write(ifle,'(1x,a)') 
     &    "ERR: [filetype='GF'] or  [filetype='Excel']"
          call FFRABORT(1,'ERR: ')
        endif
!      if(nfl/=8.and.filetype/=' ') then
!        call FFRABORT
!     & (1,'ERR:  [filetype] only defined for [fforce] file')
!      endif
!-----------------------
! --- probe file type
!-----------------------
      else if(nfl==7) then
        call nml_listno(2,cfile,filetype,kfforce)
        if(kfforce==2) then
          ProbeFiletype=2
        else
          ProbeFiletype=1
        endif
        if(trim(filetype)/=trim(cfile(1)).and.
     &     trim(filetype)/=trim(cfile(2)).and.
     &     trim(filetype)/=' ') then
          write(ifle,'(1x,a)') 
     &    "ERR: [filetype='GF'] or  [filetype='Excel']"
          call FFRABORT(1,'ERR: Re-define filetype in [&probe]')
        endif
!
      endif
!
      if((nfl/=7.and.nfl/=8).and.filetype/=' ') then
        call FFRABORT
     & (1,'ERR:  [filetype] only defined for [fforce] or [probe] file')
      endif
      
!---------------
! --- BC output
!---------------
      if(BC_result==0) then
      elseif(BC_result==1) then
        BC_output=yes
        write(ifll,'(1x,a)')
     &   'MSG: BC result will be output in result file'
      elseif(BC_result/=0.and.BC_result/=1) then
        write(ifle,'(1x,a)') 'ERR: [BC_result] MUST be 0 or 1'
        call FFRABORT(1,'ERR: defining error for BC_result')
      else
        BC_output=no
      endif
      if(BC_result==1.and.nfl/=2) then
        call FFRABORT
     &  (1,'ERR: [BC_result] only defined for [result] file')
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

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_param
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- < module for named constants >
!
      real(8) :: yslw=1.d-20
!
! yslw : lower limit of mass fraction
!
!
! --- 
!    
      
      
      end module module_param
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_simple
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(LEN=15),parameter,private :: modnam='(module_simple)'
!
!
      integer :: nsmpl=1,I_PISO=0
      real(8)  :: simeer=0.15D0
      
!
! nsmpl : no. of iteration counts in SIMPLE algorithm
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,ierror)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
! --- [namelist]
!
      integer :: iter,PISO=0
      real(8)  :: tolsimp
      namelist /simple/ iter,tolsimp,PISO
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
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
      if(PISO==0) then
        I_PISO=0
      elseif(PISO==1) then
        write(ifll,*) 'MSG: PISO algorithm is USED'
        I_PISO=1
      else
        call FFRABORT(1,'ERR: &simple/PISO = 0 or 1')
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      end subroutine inputdata
!
      end module module_simple
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_source
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(LEN=15),parameter,private :: modnam='(module_source)'
!
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
      real(8) ,save,allocatable         :: wdscnd (:,:)
      real(8) ,save,allocatable,private :: wdtscnd(:,:)
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
      implicit none
      character(LEN=10),parameter :: subnam='(chkncomp)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.ncomp ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      CALL FFRABORT(1,'module_source/chkncomp')
      end subroutine chkncomp
!
!< #2. check consistency of values in argument & module >---------------
!=================================================
      subroutine chknrans(ival,sbnm)
!=================================================
      implicit none
      character(LEN=10),parameter :: subnam='(chknrans)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.nrans ) return
      call modutl_chkset('>',1,iset1,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      CALL FFRABORT(1,'module_source/chkncomp')
      end subroutine chknrans
!
!< #3. interpolate values at present time making use of time table >----
!=================================================
      subroutine setnow(time)
!=================================================
      implicit none
      real(8),intent(in) :: time
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
      subroutine strdomain(ifle,NCV,nsdmn,nvsdmn,
     & nosdmn,icsdmn,lcsdmn,ierr)
!=================================================
      implicit none
      integer,intent(in)  :: ifle,NCV,nsdmn,nvsdmn
      integer,intent(in)  :: nosdmn(  NCV)
      integer,intent(in)  :: icsdmn(0:NCV)
      integer,intent(in)  :: lcsdmn(  NCV)
      integer,intent(out) :: ierr
!
      character(LEN=11),parameter :: subnam='(strdomain)'
      integer :: i,j,m,n
!
      call modutl_chkset('>',1,iset1,modnam//subnam)
      call modutl_chkset('=',1,iset2,modnam//subnam)
      ierr=0
      allocate( icscnd(0:nscnd),
     &          lcscnd(nvsdmn), stat=ierr )
      if( ierr.ne.0 ) then
      write(ifle,*) '### error : allocation failed'
      goto 9999
      endif
      icscnd(0)=0
      i=0
      do 100 n=1,nscnd
      do 200 m=1,nsdmn
      if( noscnd(n).eq.nosdmn(m) ) goto 201
  200 continue
      goto 100
  201 continue
      do 202 j=icsdmn(m-1)+1,icsdmn(m)
      i=i+1
      lcscnd(i)=lcsdmn(j)
  202 continue
      icscnd(n)=i
  100 continue
      return
 9999 continue
      write(ifle,*) modnam,subnam
      ierr=1
      end subroutine strdomain
!
!
!< #5. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,
     &                     ncmpx,nrnsx,lrans,ierror)
!=================================================
      implicit none
! --- MXMATERIAL
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      integer,intent(in)  :: ncmpx,nrnsx
      logical,intent(in)  :: lrans
      integer,intent(out) :: ierror
!
! --- [namelist]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      integer,parameter :: lc=10
      integer,parameter :: msrtim=10
      integer,parameter :: mcomp=200, mrans=50
      character(LEN=lc) :: kind
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
      character(LEN=lc),parameter :: mass='mass      ',
     &                      pervolume='pervolume ', 
     &                        permass='permass   '
!
      character(LEN=lc),parameter :: kndlst(3)=(/
     & mass, pervolume, permass /)
      real    :: dum1
      real(8)  :: ysw(mcomp)
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
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
      integer,intent(in)  :: i
      real   ,intent(in)  :: ys (mcomp*msrtim)
      real(8) ,intent(out) :: ysw(mcomp)
      integer,intent(out) :: ierr
      integer :: j,i0
      real(8)  :: sum
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

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_time
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(17),parameter,private :: modnam='(module_time)'
!
!
!< module for start/end time >
!
      integer :: iters=0,itere=0,i_steady=1
      logical :: steady
      integer :: MHD_steady

      real(8)  :: timee=0.d0,toltim=0.001
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,iMHD,
     &  restart,ip_mdl,icomp,ierror)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ifli,ifll,ifle,iMHD,ip_mdl,icomp
      integer,intent(inout)  :: restart
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
! --- [namelist]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      integer :: start,end,flowcon,MHD_cond
      real(8) :: end_time,toltime
      namelist /time/ start,end,end_time,flowcon,toltime,MHD_cond
!
! --- [local entities]
!
      integer,parameter :: iundef=-huge(1)
      integer :: iset=0
      integer :: ios,ierr1
!flowcon
!       =1: steady flow
!       =2: transient flow
!
!MHD_cond
!       =1: steady condition
!       =2: transient condition
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
      elseif(start==-1) then
        restart=1
      endif
! --- 
      if(ip_mdl==2.and.(start==-1.or.start==0)) then
        write(ifle,'(1X,a)') 'ERR: NO fluid field '
        write(ifle,'(1X,a)') 
     & 'MSG: cal_particle=2 : fluid field MUST be calculated firstly'
        write(ifle,'(1X,a)') 
     & 'MSG: cal_particle=2 : and then particle tracing only'
        write(ifle,'(1X,2a)') 
     &  'MSG: RESET [start=-1] or [start=0] and [cal_particle=0] ',
     &      'for fluid field cal.'
        call FFRABORT
     & (1,'ERR: firstly cal. fluid field firstly ')
      endif
!
      if(icomp==2.and.(start==-1.or.start==0)) then
        write(ifle,'(1X,a)') 'ERR: NO fluid field '
        write(ifle,'(1X,a)') 
     & 'MSG: cal_tys=2 : fluid field MUST be calculated firstly'
        write(ifle,'(1X,2a)') 
     & 'MSG: cal_tys=2 : and then set [cal_tys=2] ',
     &   'for species transport only'
        write(ifle,'(1X,2a)') 
     &  'ERR: RESET [start=-1] or [start=0] and [cal_tys=0] ',
     &      'for fluid field cal.'
        call FFRABORT
     & (1,'ERR: After cal. fluid field  for species transport only')
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
      if(flowcon==1.or.flowcon==3) then
        steady=.true.
        i_steady=flowcon
        if(toltime.gt.0.d0) then
          if(toltime.gt.1.d0) then
            write(ifle,*) '### error : data error'
            write(ifle,*) 'toltime must be smaller than one'
            CALL FFRABORT(1,'module_time/inputdata')
          else
            toltim=toltime
          endif
        else
          write(ifle,*) '### error : data error'
          write(ifle,*) 'toltime must be great than zero'
          CALL FFRABORT(1,'module_time/inputdata')
        endif
      elseif(flowcon.eq.2) then
        steady=.false.
        toltim=0.d0
        i_steady=flowcon
      else
        write(ifle,*) '### error : data error'
        write(ifle,*) 'flowcon must be 1,2 or 3'
        CALL FFRABORT(1,'module_time/inputdata')
      endif
!
      if(imhd>0) then
        if(MHD_cond==1) then
! --- J-OMEGA steady methord
          MHD_steady=1
        elseif(MHD_cond==2) then
          MHD_steady=2
! --- Unsteady methord
        elseif(MHD_cond==3) then
! --- J-OMEGA unsteady methord
          MHD_steady=3
        else
          write(ifle,*) 'ERR: data error'
          write(ifle,*) 'MSG: flowcon must be 1 or 2'
          write(ifle,*) 
     &  'MSG: MHD_cond=1: Steady condition: J-OMEGA steady methord'
          write(ifle,*) 
     &  'MSG: MHD_cond=2: Unsteady condition: '
          write(ifle,*) 
     &  'MSG: MHD_cond=3: Steady condition: J-OMEGA unsteady methord'
          CALL FFRABORT(1,'module_time/inputdata')
        endif
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
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_turbparm
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!< module for parameters of turbulent flow >
!
      real(8) :: prturb=0.9d0
      real(8) :: scturb=0.9d0
!
! prturb : turbulent Prandtl no.
! scturb : turbulent Schmidt no.
!
      real(8) :: akappa=0.42d0,
     &           yplsm=11.63
      real(8) :: awfnc =5.5d0
      real(8) :: E_vel=9.d0
!
      real*8,parameter :: Rey_s=200.d0
      real*8  :: dRey=0.2d0*Rey_s
      real*8  :: Amu=34.48d0
!STAR:34.48  !FLUENT: 70.d0  !
      
!
! akappa : Karman constant
! awfnc  : constant of log-law profile near the wall
!
      end module module_turbparm
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_species
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_material,only : isthr,nflud
      implicit none
      character(LEN=16),parameter,private :: modnam='(module_species)'
!
!< module for species data >
      integer,parameter  :: mcomp=200  ! max number of species
!
!
!-< 0. Constants >-
!
      real(8),parameter :: gascns=8.314d0 ![J/mol/K]
      real(8),parameter :: gascns_1=8.314d0/(4.18605d0) ![Cal/mol/K]
!
! gascns : unversal gas constant [J/mol/K]
!
!
!-< 1. Name of species >-
!
      integer,parameter :: sp_gas=1,sp_stie=2
      integer,private :: ns=0,ns_gas,ns_suf
!
      integer,parameter :: lenspc=20
      character(lenspc),save,allocatable :: spcnam(:)
      INTEGER,save,allocatable :: spcno(:)
      integer,save,allocatable :: act(:)
!0:no emission gas
!1:H2O,   2:   CO2,    ... see user_mannual
!40:	user define gas/particle
!
!-< 2. Thermodynamic constant >-
!
      real(8),save,allocatable :: wm(:)
      real(8),save,allocatable :: r_wm(:)
      real(8),save,allocatable :: Ri(:)
      real(8),save,allocatable :: t3(:,:)
      real(8),save,allocatable :: a7(:,:,:)
      real(8),save,allocatable :: acpk(:,:)
      real(8),save,allocatable :: hform(:),href(:),
     &                            tref_comp(:)
      real(8),save,allocatable :: Le(:)
      real(8),save,allocatable :: diffu(:)
      real(8),save,allocatable :: eps_KK(:)
      real(8),save,allocatable :: sigmaa(:)
      real(8),save,allocatable :: r_Le(:)
      real(8)                  :: tref=298.15d0
      real(8)                  :: pref=101325.d0
      logical,save,allocatable :: judg_t2(:)
      logical,save             :: sw
      real(8),save,allocatable :: Tc(:),Pc(:),KIJ(:,:),
     &                            omga(:),Rhoc(:),
     &                            dpm(:),kpp(:)
!
      real(8),save :: t_bi1, t_bi2
!      
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
      real(8),save,allocatable :: hform2(:)      
      real(8),save,allocatable :: acpk2(:,:)
      real(8)                  :: tref2=298.15d0
!
!-< 99. for internal check >-
!
      integer,private :: iset=0,icom
      integer,save      :: I_H2=0,I_O2=0,I_H2O=0,
     &                     I_PROTON=0,I_H2OL=0,I_E=0
!
! --- real fluid
!
      real(8),parameter :: p1=0.42747d0,p2=0.08664d0,
     &                     p3=0.480d0,p4=1.574d0,p5=0.176d0 
!
      real(8),save,allocatable :: Aij0(:),alph(:),alph_d(:),alph_dd(:)
!
      real(8) :: delt_cp,delt_H,delt_C,VVV,TTT,PPP,dens
      
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. check consistency of values in argument & module >---------------
!=================================================
      subroutine chkncomp(ival,sbnm)
!=================================================
      implicit none
      character(LEN=10),parameter :: subnam='(chkncomp)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if( ival.eq.ns_gas ) return
      call modutl_chkset('>',1,iset,modnam//subnam)
      write(*,*) '### program error -1- ',modnam,subnam
      write(*,*) 'sbnm=',sbnm
      CALL FFRABORT(1,'module_species/chkncomp')
      end subroutine chkncomp
!
!< #1. input namelist >-------------------------------------------------
!
!========================================================
      subroutine inputdata_gas
     &           (ifli,ifll,ifle,mmcomp,ncomp,
     &            ncomp_suf,
     &            cntlnam,my_rank,
     &            spinam,E2P,lvof,lcavi,
     &            lFC,vapor,cool,MIX_KNI,CHUN,idens,
     &            ierror)
!========================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifli,ifll,ifle,ncomp,my_rank,ncomp_suf,
     &                           mmcomp,lFC,MIX_KNI,CHUN
      character(*),intent(in) :: cntlnam
      character(*),intent(out):: spinam(mmcomp)
      logical,intent(in)      :: E2P,lvof,lcavi
      integer,intent(in)      :: idens
      integer,intent(out)     :: ierror,vapor,cool
      

!
! --- [namelist]2
!
      character(lenspc)  :: name        ! species name
      real(8)            :: weight      ! molecular weight [kg/mol]
      real(8)            :: a(7,2)      ! coefficient set
      real(8)            :: t(3) 
                                        ! temperature [K] (1:upper limit,
                                        ! 2:branch point, 3:lower limit)
      real(8)            :: cp(5)       ! coefficient set
      real(8)            :: ho          ! enthalpy [J/kg]
      real(8)            :: Lewis       ! Lewis number [-]
      real(8)            :: weight2,cp2(5),ho2,tref2
      real(8)            :: diffusity,Eps_K,sigma
      integer            :: flag_a7=0,fa7
      character(200)     :: th_file,sp_file ! input file name
      character(200)     :: thfile,spfile
      integer            :: bg_gas=-1
      integer            :: gasno = 0
      integer            :: active=1
      real(8)            :: T_c,P_c,rho_c,K_ij(mcomp),OMEGA,dipole,kappa
      integer,save       :: real_mdl=0
!
      namelist /species/ name,weight,
     &                   cp ,ho ,tref,
     &                   cp2,ho2,tref2,
     &                   flag_a7,a,t,Lewis,
     &                   th_file,bg_gas,
     &                   diffusity,gasno,active,
     &                   Eps_K,sigma,
     &                   T_c,P_c,OMEGA,K_ij,Rho_c,
     &                   dipole,kappa
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      real(8),parameter :: undef=-huge(1.d0)
      real(8) :: hk(5),tt,h00
      integer :: i,ios=0,iosf=0,ierr1,s,s2,nsno,nofind,idum
      logical :: nml_comp_eq,nml_comp_gt,Lel,existth=.false.
!
!      integer,parameter :: lc=7
!      character(lc),parameter :: chemkd(2)=(/'gas    ','surface'/)
!      character(lc),parameter :: gas='gas',surface='surface'
      character(lenspc) :: spienm
      character(1) :: percnt
!
      ierror=0
      ns=ncomp+ncomp_suf
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
! flag_a7=0 : NASA five coefficient polynomial,heat-formation is necessary
! flag_a7=1 : Read file therm.dat, 
!             NASA 7 coefficient polynomial or input in fflow.ctl
!
!-< 1. allocate array >-
!
      allocate( spcnam(ns),      ! name of species
     &          wm    (ns),      ! molecular weight [kg/mol]
     &          r_wm  (ns),      ! reciprocal of 'wm' [mol/kg]
     &          Ri    (ns),      ! 'gascns'/'wm' [J/(kg K)]
     &          spcno (ns),	 ! gas/particle number, see user_manual,
     &          act(ns),
     &          eps_KK(ns),
     &          sigmaa(ns),
     &          Tc(ns),Pc(ns),Rhoc(ns),
     &          KIJ(ns,ns),omga(ns),dpm(ns),kpp(ns),
     &          Aij0(ns),alph(ns),alph_d(ns),alph_dd(ns),
     &          stat=ierr1 )
!      call memory_error
!------------------------------------------------
! --- judgment of coefficient set "a7" or "acpk"
!------------------------------------------------
 1002 continue
      flag_a7=0
      fa7=flag_a7
      nsno=0
      th_file = " "
      sp_file = " "
      thfile = " "
      spfile = " "
      act(:)=1
      eps_KK(:)=undef
      sigmaa(:)=undef
!
      Tc(:)=0.d0
      Pc(:)=0.d0
      omga(:)=0.d0
      KIJ(:,:)=0.d0
      Rhoc(:)=0.d0
      dpm(:)=0.d0
      kpp(:)=0.d0
!
      rewind ifli
      do
        Eps_K=undef
        sigma=undef
        T_c=undef
        P_c=undef
        OMEGA=undef
        K_ij(:)=undef
        Rho_c=undef
        dipole=undef
        kappa=undef
        active=1
        read(ifli,species,iostat=ios)
        if(ios<0) then
          exit
        elseif(ios==0) then
          nsno=nsno+1
          if(nsno==1) then
            fa7=flag_a7
            if(fa7==1) then
              if(th_file/=' ') then
                existth=.false.
                thfile=trim(th_file)
                inquire(file=thfile,exist=existth)
                if(.NOT.existth) then
                  write(ifle,'(2a)') 
     &            'ERR: Cannot found thermal file: ',trim(thfile)
                  write(ifle,'(2a,1x,a)') 
     &            'MSG: Copy thermal file ',trim(thfile),
     &             'into current directory'
            call FFRABORT(1,'NO thermal file') ! in first [&species]')
                endif
              endif
              if(sp_file/='') then
                spfile=trim(sp_file)
                existth=.false.
                inquire(file=spfile,exist=existth)
                if(.NOT.existth) then
                  write(ifle,'(2a)') 
     &            'ERR: Cannot found species file : ',trim(spfile)
                  write(ifle,'(2a)') 
     &            'MSG: Copy thermal file [species.dat] ',
     &             'into current directory'
                  call FFRABORT(1,'NO [species.dat] file')
                endif
              endif
            endif
          endif
          if(.not.(flag_a7==1.or.flag_a7==0)) then
            write(ifle,*) '### ERR: all flag_a7 should be 0 or 1'
            call FFRABORT(1,' Reading error-1 in [&species]')
          endif
          if(flag_a7/=fa7) then
         write(ifle,*) 
     &  '### ERR: NOT same flag_a7 value at species no=',nsno
            write(ifle,*) '### ERR: all flag_a7 should be 0 or 1'
            write(ifle,*) 
     &      '### MSG: flag_a7=1 : NASA seven coefficient polynomial '
            write(ifle,*) 
     &      '### MSG: flag_a7=0 : NASA five coefficient polynomial '
            call FFRABORT(1,' Reading error-2 in [&species]')
          endif
          if(active==0.or.active==1) then
             act(nsno)=active
          else
            call FFRABORT
     &  (1,'ERR: active=1 or 0 for defining active/passive scalar')
          endif
        else
          write(ifli,species)
          call FFRABORT(1,' Reading error-3 in [&species]')
        endif
!
        if(MIX_KNI==1) then
          write(ifle,'(1X)') 
          write(ifle,'(1X,a)') 'MSG: MIX KNI gas model is defined : '
          write(ifle,'(1X,2a)') 'MSG: species name :',trim(name)
          if(Eps_K==undef) then
            write(ifle,'(1X,a,I4)') 
     &      'ERR: Not Defining Eps_K in &species no',nsno
            write(ifle,'(1X,2a)') 
     &      'MSG: [Eps_K] is Epslon/K value for current species'
            call FFRABORT(1,'ERR: Epslon/K for MIX Kinetic Model')
          else
            write(ifll,*) 'MSG: Eps_K=',Eps_K
          endif
          if(sigma==undef) then
            write(ifle,'(1X,a,I4)') 
     &      'ERR: Not Defining sigma in &species no',nsno
            write(ifle,'(1X,2a)') 
     &      'MSG: [sigma] is Characteristic diameter ',
     &      'of the molecule in Angstroms'
            call FFRABORT(1,'ERR: sigma for MIX Kinetic Model')
          else
            write(ifll,*) 'MSG: sigma=',sigma
          endif
          
          Eps_KK(nsno)=Eps_K
          sigmaa(nsno)=sigma
          write(ifle,'(1X)') 
          
        endif
!
! --- 
!

        if(CHUN==1) then
          write(ifle,'(1X)') 
          write(ifle,'(1X,2a)') 'MSG: species name :',trim(name)
          real_mdl=0
          if(Eps_K/=undef.and.sigma/=undef) then
            real_mdl=1
            write(ifle,'(1X,a)') 'MSG: Eps_K & sigma input by user :'
          endif
          if(T_c/=undef.and.Rho_c/=undef.and.real_mdl/=1) then
            real_mdl=2
            write(ifle,'(1X,a)') 
     &          'MSG: Eps_K & sigma calculated with T_c & Rho_c :'
            write(ifle,'(1X,a)') 
     &          'MSG: Eps_K=T_c/1.2593D0 '
            write(ifle,'(1X,a)') 
     &          'MSG: sigma=0.809D0*(1.d6*weight/Rho_c)**(1.0D0/3.0D0)'
          endif
!
          if(real_mdl/=1.and.real_mdl/=2) then
            write(ifle,'(1X,a,I4)')  
     &      'ERR: Not Defining error in &species no',nsno
            call FFRABORT(1,'ERR: [Eps_K & sigma] or [T_c & Rho_c]')
          endif

          if(real_mdl==1) then
            Eps_KK(nsno)=Eps_K
            sigmaa(nsno)=sigma
          elseif(real_mdl==2) then
            Eps_KK(nsno)=T_c/1.2593D0
            sigmaa(nsno)=0.809D0*(1.d6*weight/Rho_c)**(1.0D0/3.0D0)
          else
            if(Eps_K==undef) then
              write(ifle,'(1X,a,I4)') 
     &        'ERR: Not Defining Eps_K in &species no',nsno
              write(ifle,'(1X,2a)') 
     &        'MSG: [Eps_K] is Epslon/K value for current species'
            endif
            if(sigma==undef) then
              write(ifle,'(1X,a,I4)') 
     &        'ERR: Not Defining sigma in &species no',nsno
              write(ifle,'(1X,2a)') 
     &        'MSG: [sigma] is Characteristic diameter ',
     &        'of the molecule in Angstroms'
            endif
!
            if(T_c==undef) then
              write(ifle,'(1X,a)') 'ERR: NOT defining T_c for real gas'
            else
              write(ifll,*) 'MSG: T_c=',T_c
            endif

            if(Rho_c==undef) then
              write(ifle,'(1X,a)') 
     &          'ERR: NOT defining Rho_c for real gas'
            else
              write(ifll,*) 'MSG: Rho_c=',Rho_c
            endif
!
            write(ifle,'(1X,a,I4)')  
     &      'ERR: Not Defining error in &species no',nsno
            call FFRABORT(1,'ERR: [Eps_K & sigma] or [T_c & Rho_c]')
          endif
!
          write(ifll,*) 'MSG: Eps_K=',Eps_KK(nsno)
          write(ifll,*) 'MSG: sigma=',sigmaa(nsno)
          if(kappa==undef) then
            write(ifle,'(1X,a)') 'ERR: NOT defining kappa for real gas'
            call FFRABORT(1,'ERR: kappa for real gas')
          else
            write(ifll,*) 'MSG: kappa=',kappa
          endif

          if(dipole==undef) then
          write(ifle,'(1X,a)') 
     &            'ERR: NOT defining dipole moment for real gas'
             call FFRABORT(1,'ERR: dipole for real gas')
          else
             write(ifll,*) 'MSG: dipole moment =',dipole
          endif
!
          if(T_c==undef) then
            write(ifle,'(1X,a)') 'ERR: NOT defining T_c for real gas'
            call FFRABORT(1,'ERR: T_c for real gas')
          else
            write(ifll,*) 'MSG: T_c=',T_c
          endif


          Tc(nsno)=T_c
          kpp(nsno)=kappa
          dpm(nsno)=dipole
        endif
!-----------------------
! --- real gas
!-----------------------
        if(idens==4) then
!
          if(CHUN==0) then
            call FFRABORT(1,'ERR: Real GAS MUST: &fluid/muopt="Chung"')
          endif
          write(ifle,'(1X,a)') 'MSG: real gas model is defined : '
!          write(ifle,'(1X,2a)') 'MSG: species name :',trim(name)

          if(T_c==undef) then
            write(ifle,'(1X,a)') 'ERR: NOT defining T_c for real gas'
            call FFRABORT(1,'ERR: T_c for real gas')
          else
            write(ifll,*) 'MSG: T_c=',T_c
          endif

          if(P_c==undef) then
            write(ifle,'(1X,a)') 'ERR: NOT defining P_c for real gas'
            call FFRABORT(1,'ERR: P_c for real gas')
          else
            write(ifll,*) 'MSG: P_c=',P_c
          endif


          if(OMEGA==undef) then
            write(ifle,'(1X,a)') 'ERR: NOT defining OMEGA for real gas'
            call FFRABORT(1,'ERR: OMEGA for real gas')
          else
            write(ifll,*) 'MSG: OMEGA=',OMEGA
          endif
!
          do i=1,ncomp
          if(K_ij(i)==undef) then
            write(ifle,'(1X,a,I4)') 'ERR: i=',i
            write(ifle,'(1X,a)') 
     &   'ERR: NOT defining K_ij(i) for real gas'
            call FFRABORT(1,'ERR: K_ij(i) for real gas')
          else
            KIJ(nsno,i)=K_ij(i)
            if(ncomp==1) KIJ(nsno,i)=0.d0
            write(ifll,*) 'MSG: i=,K_ij(i)=',i,KIJ(nsno,i)
          endif
          enddo

          Tc(nsno)=T_c    !
          Pc(nsno)=P_c
!          Rhoc(nsno)=Rho_c  
          omga(nsno)=OMEGA
!          dpm(nsno)=dipole
!          kpp(nsno)=kappa

        endif

      enddo
      write(ifle,'(1X)') 
!
!
!
      if(nsno.ne.ncomp) then
        write(ifle,*) ' ### ERR: ncomp /= nsno'
        write(ifle,*) ' ### MSG: ncomp= ',ncomp,' in file:',cntlnam
        write(ifle,*) ' ### MSG: [&species] number : ',nsno
        write(ifle,*) ' ### reading error in "&species"'
        write(ifle,*) ' ### MSG: Re-run prefflow '
        call FFRABORT(1,' Reading error in [&species]')
      endif
      ns_gas=nsno
      sw=.true.                        ! for "a7"
      if(fa7.eq.0)  sw=.false.         ! for "acpk"
!----------------------------------------------------
! --- reading species name, molecular weight [kg/mol]
!----------------------------------------------------
      spcnam=' '                      ! name of species
      spinam=' '     
      wm= 0                            ! molecular weight [kg/mol]
      r_wm= 0                          ! reciprocal of 'wm' [mol/kg]
      Ri= 0                            ! 'gascns'/'wm' [J/(kg K)]
      spcno=0
      
!      BGCOMP=0
      rewind ifli
      do s=1,ns_gas
        name=' '                       ! name of species
        weight = undef                 ! molecular weight [kg/mol]
        bg_gas=-1
        read(ifli,species,iostat=ios)
        call reading_judge             ! reading and judge of data file
        if(name==' ') then             ! check "name", name of species
          write(ifle,'(a)') 'ERR: reading error in [&species]'
          write(ifle,'(a,I4)') 'ERR: lack of species [name] in no.',s
          ierror=1
        endif
        spcnam(s)=name
	spinam(s)=name
        if(s>=2) then
          do s2=1,s-1
            if( spcnam(s)==spcnam(s2) ) then
              write(ifle,'(A)') 'ERR: reading error in [&species]'
              write(ifle,'(2A)') 
     &           'ERR: duplicated data : name(i)=name(j)=',
     &                       spcnam(s2)
              write(ifle,'(A,2I4)') 'i,j =',s2,s
              ierror = 1
            endif
          enddo
        endif
!------------------------------------------------
! --- check "weight", molecular weight [kg/mol]
!------------------------------------------------
        if(.not.nml_comp_eq("weight",weight,undef,s,ifle,ierror)) then
          if(nml_comp_gt("weight",weight,0.0d0,s,ifle,ierror)) then
            wm(s) = weight              ! molecular weight [kg/mol]
            r_wm(s) = 1.0d0/wm(s)       ! reciprocal of 'wm' [mol/kg]
            Ri(s) = gascns*r_wm(s)      ! 'gascns'/'wm' [J/(kg K)]
          endif
        endif
!-------------------------------------------------
! --- input and check backgrand gas: N2 or Ar etc
!-------------------------------------------------
        if(bg_gas/=-1) then
          call FFRABORT(1,'ERR: bg_gas will defined at &fluid')
        endif
!
        
!        if(ncomp==1) then
!          BGCOMP=1
!        elseIF(bg_gas==1.and.BGCOMP==0) THEN
!          BGCOMP=s
!          if(my_rank==0) then
!            write(ifll,'(2a)') 'MSG: BackGrand Gas is ',trim(spcnam(s))
!          endif
!        elseIF(bg_gas==1.and.BGCOMP/=0) THEN
!          WRITE(ifle,'(A)') 
!     &    'ERR: Only one gas can be defined as BackGround'
!          WRITE(ifle,'(2A)') 'MSG: Stop at :',trim(spcnam(s))
!          call FFRABORT(1,'ERR: ')
!        elseif(bg_gas/=1.and.bg_gas/=0) then
!          WRITE(ifle,'(2A)') 
!     &    'ERR: bg_gas MUST BE 0 OR 1 in ',trim(cntlnam)
!          call FFRABORT(1,'ERR: ')
!        endif
         if(gasno.lt.0.or.gasno.gt.40) then
	   WRITE(ifle,'(2A)') 
     &     'ERR: gasno should be in (0, 40) in',trim(cntlnam)
           WRITE(ifle,'(A)') '0: non-radiative media' 
	   WRITE(ifle,'(A)') '1~39: see manual, 40: user define' 
	   WRITE(ifle,'(2A)') 'if user define media >=2, you can specify',
     &	   'a number in 1~39 if without confliction.'
           call FFRABORT(1,'ERR: ')
	 endif
         spcno(s)=gasno
      enddo
!      if(BGCOMP==0) then
!        WRITE(ifle,'(A)') 'ERR: BackGround gas must be defined'
!        WRITE(ifle,'(3A)') 
!     &   'MSG: Set [bg_gas=1] in one [&species] for ',
!     &   'defining BackGround species in ',trim(cntlnam)
!        call FFRABORT(1,'ERR: ')
!      endif
!
      if( s-1/=ns_gas ) then                ! check number of species
        write(ifle,'(a)') '### error : data error'
        write(ifle,'(a,I4)') 'no. of namelists (species) =',s-1
        write(ifle,'(a,I4)') 'it must be =',ns
        ierror = 1
      endif
!
!--------------------------------------------------------
! --- reading coefficient sets and display reading result
!--------------------------------------------------------
      if(sw) then  !a7-model
        allocate( a7(7,2,ns),    ! coefficient sets
     &            t3(  3,ns),    ! temperature [K]
     &            judg_t2(1),    ! judgment "t3(2,:)"
                                 ! . every value is same or NOT
     &            hform(ns),href(ns),tref_comp(ns),
     &            diffu(ns),
     &            stat=ierr1)
!        call memory_error
!----------------
! --- read again
!----------------
        hform(:)=0.d0
        href(:)=0.d0        ![J/mol]
        tref_comp(:)=0.d0
        t3(:,:)=0.d0
        a7(:,:,:)=0.d0
        diffusity=undef
!
        rewind ifli
        do s=1,ns_gas
        a=undef                     ! coefficient
        t=0                         ! temperature [K]
        tref=undef
        read(ifli,species,iostat=ios)
        call reading_judge          ! reading and judge of data file
        nofind=0
        if(a(1,1)==undef) then
          spienm=spcnam(s)
          do i=1,lenspc
          percnt=spcnam(s)(i:i)
          if(percnt=='%') then
            spienm=trim(spcnam(s)(1:i-1))
            exit
          endif
          enddo
          if(th_file/=' ') then
            call from_files(ifle,s,a,t,LEN_TRIM(th_file),
     &       LEN_TRIM(sp_file),thfile,spfile,spienm,nofind)
          else
            write(ifle,'(a)') 
     &  'ERR: Not defined thermal file in first [&species]'
            write(ifle,'(a,/,a)') 
     &  'MSG: flag_a7=1 => NASA thermal file [therm.dat] will be used',
     &  '                  or input NASA 7 coefficients and t(1~3)'
            write(ifle,'(a,/,a)') 
     &  'MSG: flag_a7=0 => NASA five coefficient polynomial ',
     &  '                  and input parameter: cp, ho,tref,weight'
            call FFRABORT(1,'NO thermal file')
          endif
          if(nofind==1.and.a(1,1)==undef) then
            write(ifle,'(2X,4a)') 'ERR: NOT found Species name: ',
     &      trim(spcnam(s)),' in thermal file: ',trim(thfile)
            write(ifle,'(2X,2a)')
     &           "WRN: Please input a7 coeffient in ",trim(cntlnam)
            call 
     &      FFRABORT(1,'ERR: NOT found species name in thermal file')
          else
            call read_a7(s,ns,ifle,a,t,tref,undef,a7,t3,ierror)
            if(my_rank==0) then
              write(ifle,'(2x,4a)') 
     &        'MSG: t(1:3), a(1:14) read from ',
     &         trim(thfile), ' for species: ',trim(name)
            endif
          endif
        else
          call read_a7(s,ns,ifle,a,t,tref,undef,a7,t3,ierror)
          if(my_rank==0) then
            write(ifle,'(2x,4a)') 
     &      'MSG: t(1:3), a(1:14) read from ',
     &      trim(cntlnam), ' for species: ',trim(name)
          endif
        endif
!
        tref=dble(tref)
        href(s)=h_t(s,tref) ![---: H0/R]
!
        if(idens==4) then
          href(s)=h_tr(ns_gas,s,tref,Tc(s),Pc(s)) ![---: H0/R]
        endif
!
        tref_comp(s)=tref
        if(lFC>0) then
          if(diffusity==undef) then
            call FFRABORT(1,'ERR: define &species/diffusity for PEFC')
          else
            diffu(s)=diffusity
          endif
          if(trim(spinam(s))=='H2') then
            I_H2=s
          endif
          if(trim(spinam(s))=='O2') then
            I_O2=s
          endif
          if(trim(spinam(s))=='H2O-G') then
            I_H2O=s
            vapor=I_H2O
          endif
          if(trim(spinam(s))=='H+') then
            I_PROTON=s
          endif
          if(trim(spinam(s))=='E-') then
            I_E=s
          endif
          if(trim(spinam(s))=='H2O-L') then
            I_H2OL=s
            cool=I_H2OL
          endif
        endif
!
        enddo
!
        if(lFC>0) then
          if(I_H2==0) then
            call FFRABORT(1,'MSG: H2 for PEFC in &species')
          endif
          if(I_O2==0) then
            call FFRABORT(1,'MSG: O2 for PEFC in &species')
          endif
          if(I_H2O==0) then
            call FFRABORT(1,'MSG: H2O-G for PEFC in &species')
          endif
          if(I_PROTON==0) then
            call FFRABORT(1,'MSG: H+ for PEFC in &species')
          endif
          if(I_H2OL==0) then
            call FFRABORT(1,'MSG: H2O-L for PEFC in &species')
          endif
          if(I_E==0) then
            call FFRABORT(1,'MSG: E- for PEFC in &species')
          endif
        endif
!
!
        do s=1,ns_gas
          hform(s)=href(s)*Ri(s)    ![J/kg]  !Ri(:)     : 'gascns'/'wm' [J/(kg K)]
        enddo

        judg_t2(1)=.true.             ! every t3(2,*) is same
        do s=2,ns_gas
          if(t3(2,1)/=t3(2,s)) judg_t2(1)=.false.
        enddo
        call display_a7                    ! result of reading
!
        
      else  !a5-model
!--------------
! --- a5 model
!--------------
        allocate(acpk(0:5,ns) ,hform(ns), ! coefficient set
     &           acpk2(0:5,ns),hform2(ns),
     &           stat=ierr1)
!        call memory_error                  ! memory error check
        acpk=0.d0                          ! coefficient
        acpk2=0.d0
        tref=undef                         ! temperature [K]
        tref2=undef
        rewind ifli
        do s=1,ns_gas
          cp   = 0.d0                      ! coefficient
          cp(1)= undef
          ho   = undef                     ! enthalpy [J/kg]
          cp2   = 0.d0 
          cp2(1)= undef
          ho2   = undef
          read(ifli,species,iostat=ios)
          call reading_judge    ! reading and judge of data file
          call read_a5(s,ns,ifle,cp,tref,ho,undef,acpk,ierror)
          if(E2P.or.lvof.or.lcavi) then
            call read_a5(s,ns_gas,ifle,cp2,tref2,ho2,undef,acpk2,ierror)
          endif

          hform(s)=ho
          
          hform2(s)=ho2
        enddo
        call display_a5         ! result of reading
      endif
!
! --- Lewis number
!
      Lel=.false.
      do i=1,nflud
        if(isthr(i)==3) then
           Lel=.true.
        endif
      enddo
      if(Lel) then
        allocate(Le    (ns),
     &           r_Le  (ns),
     &           stat=ierr1)
!        call memory_error
        Le(1:ns) = 0.d0
        r_Le(1:ns) = 0.d0
        rewind ifli
        do s=1,ns_gas
          Lewis = undef
          read(ifli,species,iostat=ios)
          call reading_judge
          if(.not.nml_comp_eq("Lewis",Lewis,undef,s,ifle,ierror))
     &    then
            if( nml_comp_gt("Lewis",Lewis,0.0d0,s,ifle,ierror) )then
              Le(s) = Lewis
              r_Le(s) = 1.d0/Le(s)
            endif
          endif
        enddo
        call disp_Lewis   ! result of reading
      endif
!----------------------
! --- final error check
!----------------------
      if( ierror/=0 ) then
        write(ifle,'(a,I4,2a)') "total error : ",ierror,modnam,subnam
      endif
!
      if(sw) then  !a7-model
	t_bi1=t3(1,1)
	t_bi2=t3(3,1)
	do s=1,ns
	t_bi1=min(t_bi1,t3(1,s))
	t_bi2=max(t_bi2,t3(3,s))
	end do
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      return
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine reading_judge        ! reading and judge of data file
!=================================================
      if(s/=ns_gas.and.ios<0) then
        write(ifle,*) '### reading error in "species"'
        ierror = 1+ierror
      endif
      call nml_errmsg0(ifle,ios,'species',ierr1)
      if( ierr1/=0 ) ierror = 1+ierror
      end subroutine reading_judge
!
!
!=================================================
      subroutine display_a7             ! display of reading result
!=================================================
      if(my_rank==0) then
        write(ifll,*)
        write(ifll,5001) 'GAS-Phase Species List: '
        write(ifll,"(2x,a,i3)") "number of chemical species in GAS : ",
     &    ns_gas
 5001   format(2x,a,80('='))
      endif
      do s=1,ns_gas
        if(my_rank==0) write(ifll,'(i3,1x,a20)') s,spcnam(s)
                          ! name of species
        if(my_rank==0) 
     & write(ifll,'(13x,a,e13.5,a)') 'mi =',wm(s),' [kg/mol]'
                          ! molecular weight [kg/mol]
        if(my_rank==0) 
     & write(ifll,'(13x,a,e13.5,a)') 'Ri =',Ri(s),' [J/(K kg)]'
                          ! 'gascns'/'wm' [J/(kg K)]
        if(my_rank==0) 
     &  write(ifll,'(13x,a,3f8.1,a)') 'T =',(t3(i,s),i=1,3),' [K]'
                          ! point of temperature region [K]
        if(my_rank==0) write(ifll,'(13x,a,$)') 'a ='
        if(my_rank==0) write(ifll,'(7e13.5)') (a7(i,1,s),i=1,7)
                          ! coefficient set of upper region
!        if(my_rank==0) write(ifll,*)
        if(my_rank==0) write(ifll,'(16x,7e13.5)') (a7(i,2,s),i=1,7)
                          ! coefficient set of lower retion
        if(my_rank==0) write(ifll,*)
        if(my_rank==0) 
     &   write(ifll,'(13x,a,E13.5,1x,a,3X,E13.5,1x,a)') 
     &   'Formation enthalpy : ',hform(s),
     &   '[J/kg]',href(s)/4.186d0,'[kcal/mole]'
      enddo
      if(my_rank==0) then
      if(judg_t2(1))
     &  write(ifll,"(2x,5('='),a)") "every T[2] is same value!"
      endif
      end subroutine display_a7
!
!=================================================
      subroutine display_a5             ! display of reading result
!=================================================
      if(my_rank==0)
     &write(ifll,"(1x,a,i3)") "number of chemical species in GAS: ",
     &   ns_gas
      do s=1,ns_gas
        if(my_rank==0) write(ifll,'(i3,1x,a8,$)') s,spcnam(s)
                          ! name of species
        if(my_rank==0)
     & write(ifll,'(1x,a,e13.5,a)') 'mi =',wm(s),' [kg/mol]'
                          ! molecular weight [kg/mol]
        if(my_rank==0)
     & write(ifll,'(13x,a,e13.5,a)') 'Ri =',Ri(s),' [J/(K kg)]'
                          ! 'gascns'/'wm' [J/(kg K)]
        if(my_rank==0) write(ifll,'(13x,a,$)') 'cp ='
        if(my_rank==0) write(ifll,'(5e13.5)') (acpk(i,s),i=1,5)
                          ! coefficient set
        if(my_rank==0) 
     &   write(ifll,'(13x,a,E13.5,1x,a,3X,E13.5,1x,a)') 
     &   'Formation enthalpy : ',hform(s),
     &   '[J/kg]'!,href(s)/4.186d0,'[kcal/mole]'
        if(my_rank==0) write(ifll,*)
      enddo
      end subroutine display_a5
!
!-------------------------------------------------
! --- reading and check data
!=======================================================
      subroutine read_a7(s,ns,f,a,t,tref,undef,a7,t3,e)
!=======================================================
      integer,intent(in)  :: s
                  ! classification number of species
      integer,intent(in)  :: ns
                  ! number of species
      integer,intent(in)  :: f
                  ! file number to write error script
      real(8), intent(in)  :: a(7,2)
                  ! coefficient set
      real(8), intent(in)  :: t(3),tref
                  ! temperature [K]
                  ! (1:upper limit, 2:branch point, 3:lower limit)
      real(8), intent(in)  :: undef
                  ! judgment value
      real(8), intent(out) :: a7(7,2,ns)
                  ! coefficient sets
      real(8), intent(out) :: t3(3,ns)
                  ! temperature [K]
                  ! (1:upper limit, 2:branch point, 3:lower limit)
      integer,intent(out) :: e
                  ! total number of error
      integer :: i,j
      real(8)  :: tt
!--------------------------------
! --- check "a7", coefficient set
!--------------------------------
      A7_0:do i=1,2
        do j=1,7
          if(.not.nml_comp_eq("a",a(j,i),undef,s,f,e)) then
            a7(j,i,s)=a(j,i)
          else
            write(f,*) '"a" needs 14 values'
            exit A7_0
          endif
        enddo
      enddo A7_0
!-------------------------------------------------
! --- check "t3", point of temperature region [K]
! --- t(1):upper limit, t(2):branch point, t(3):lower limit
!-------------------------------------------------
      if( t(1)>t(2) .and. t(2)>t(3) .and. t(3)>0.d0 ) then
        t3(1:3,s)=t(1:3)
      else
        write(f,*) '### reading error in "species"'
        write(f,*) 'check "t" of no.',s
        write(f,'(a,3F14.5)') 'MSG: t(:)=',t(:)
        e=1+e
        call FFRABORT(1,'ERR: t(1:3) reading erorr')
      endif
!
      if( .not.nml_comp_eq("tref",tref,undef,s,f,e) ) then
        if( nml_comp_gt("tref",tref,0.0d0,s,f,e) )
     &    tt = tref
      endif
      if(e/=0) then 
        call FFRABORT(1,'ERR: t(1:3) reading erorr')
      endif
      end subroutine read_a7
!
!=================================================
      subroutine read_a5(s,ns,f,cp,tref,ho,undef,acpk,e)
                  ! reading and check data
!=================================================
      integer,intent(in)  :: s
                  ! classification number of species
      integer,intent(in)  :: ns
                  ! number of species
      integer,intent(in)  :: f
                  ! file number to write error script
      real(8), intent(in)  :: cp(5)
                  ! coefficient set
      real(8), intent(in)  :: tref
                  ! temperature [K]
      real(8), intent(in)  :: ho
                  ! enthalpy [J/kg]
      real(8), intent(in)  :: undef
                  ! judgment value
      real(8), intent(out) :: acpk(0:5,ns)
                  ! coefficient set
      integer,intent(out) :: e
                  ! total number of error
      integer i
      real(8) hk(5),tt
!
!----------------------------
! --- check "cp", coefficient
!----------------------------
      if(.not.nml_comp_eq("cp/cp2",cp(1),undef,s,f,e)) then
        if(nml_comp_gt("cp/cp2",cp(1),0.0d0,s,f,e)) then
          do i=1,5
            hk(i)=cp(i)/dble(i)
          enddo
        endif
      endif
!------------------------------------
! --- check "tref", temperature [K]
!------------------------------------
      if(.not.nml_comp_eq("tref/tref2",tref,undef,s,f,e)) then
!        if( nml_comp_gt("tref/tref2",tref,0.0d0,s,f,e)) then
        tt=tref
!        endif
      endif
!
! --- check "ho", formation enthalpy [J/kg]
!
      if(.not.nml_comp_eq("ho/ho2",ho,undef,s,f,e)) then
        acpk(0,s)=ho-((((hk(5) *tt
     &                  +hk(4))*tt
     &                  +hk(3))*tt
     &                  +hk(2))*tt
     &                  +hk(1))*tt
        do i=1,5
          acpk(i,s) = cp(i)             ! coefficient set
        enddo
      endif
      end subroutine read_a5
!=================================================
      subroutine disp_Lewis
!=================================================
      if(my_rank==0) write(ifll,'(2X,a)') "Lewis number of species:"
      if(my_rank==0) then
         do s=1,ns_gas
         write(ifll,'(2x,a8,2X,f8.3)') (spcnam(s)), Le(s)
         enddo
      endif
      if(my_rank==0) write(ifll,*)
      end subroutine disp_Lewis
!
      end subroutine inputdata_gas
!
!=========================================================================
      subroutine from_files
     &        (ifle,s,a77,t,len_th,len_sp,thfile,spfile,spienm,nofind)
!=========================================================================
      integer,intent(in) :: len_th,len_sp ! length of file name
      integer,intent(in) :: s,ifle
      real*8 ,intent(out):: a77(7,2),t(3)
      integer,intent(inout) :: nofind
      character(*),intent(in) :: thfile,spfile,spienm
!
      integer,parameter :: ird=44	  ! number to open file
      integer i,s2			  ! counter
      integer w,ios				  ! width of str_t
      character(120) :: str,str_t	  ! input data string
      character(len_th) :: th_		  ! file name of thermochemical data
      character(len_sp) :: sp_		  ! file name of species data
      character(3),parameter :: str_e3 = "e-3"
!
      th_=trim(thfile)
      sp_=trim(spfile)
      ios=0
!
      if(th_==' ') call FFRABORT(1,'ERR: th_file==???')
      open(ird,file=th_,iostat=ios,form='formatted')
      if(ios>0) then		! open error
        write(ifle,*) "cannot open file '",th_,"'"
        call FFRABORT(1,'ERR:')
      endif
      rewind ird
      do             ! species name search
        read(ird,'(a)',iostat=ios) str
        if(ios<0) then		! end of file
          write(ifle,'(2X,4a)') "WRN: cannot find '",trim(spienm),
     &                   "' in file :",trim(th_)
          nofind=1
          return
        elseif(str(80:80)=="1") then	! 1st line search
          do i=1,18
             if(str(i:i)==" ") exit		! 1st space search
          enddo
          if(spienm==str(1:i-1)) exit	! compare name
        endif
      enddo
      read(str(46:55),*) t(3)	! lower temperature [K]
      read(str(56:65),*) t(1)	! upper temperature [K]
      read(str(66:73),*) t(2)	! branch temperature [K]
      read(ird,'(a)',iostat=ios) str	! 2nd line
      if(str(80:80)=="2") then
         do i=13,80,15		! check spase after "E"
         if(str(i:i)==" ") str(i:i)="+"
         enddo
         do i=1,5
         read(str(15*(i-1)+1:15*i),*) a77(i,1)
         enddo
      endif
      read(ird,'(a)',iostat=ios) str	! 3rd line
      if(str(80:80)=="3") then
         do i=13,80,15		! check spase after "E"
         if(str(i:i)==" ") str(i:i)="+"
         enddo
         read(str( 1:15),*) a77(6,1)
         read(str(16:30),*) a77(7,1)
         do i=3,5
         read(str(15*(i-1)+1:15*i),*) a77(i-2,2)
         enddo
      endif
      read(ird,'(a)',iostat=ios) str	! 4th line
      if(str(80:80)=="4") then
         do i=13,80,15		! check spase after "E"
         if(str(i:i)==" ") str(i:i)="+"
         enddo
         do i=1,4
         read(str(15*(i-1)+1:15*i),*) a77(i+3,2)
         enddo
      endif
      close(ird)
!
      end subroutine from_files
!
!=================================================
      subroutine inputdata_suf
     & (ifli,ifll,ifle,cntlnam,my_rank,mmcomp,
     &  ncomp,ncomp_suf,spinam,ierror)
!=================================================
      implicit none
!
!
! --- [dummy arguments]
!
      integer,intent(in)      :: ifli,ifll,ifle,
     &                           ncomp_suf,ncomp,
     &                           my_rank,mmcomp
      character(*),intent(in) :: cntlnam
      character(*),intent(out):: spinam(mmcomp)
      integer,intent(out)     :: ierror
!
! --- [namelist]
!
      character(lenspc)  :: name        ! species name
      real(8)            :: weight      ! molecular weight [kg/mol]
      real(8)            :: a(7,2)      ! coefficient set
      real(8)            :: t(3)
                                        ! temperature [K] (1:upper limit,
                                        ! 2:branch point, 3:lower limit)
      real(8)            :: cp(5)       ! coefficient set
      real(8)            :: ho          ! enthalpy [J/kg]
      real(8)            :: Lewis       ! Lewis number [-]
      integer            :: flag_a7=0,fa7
      character(200)     :: th_file,sp_file ! input file name
      character(200)     :: thfile,spfile
      integer,parameter  :: mcomp=200       ! max number of species
      namelist /surface_species/ name,weight,tref,
     &                           flag_a7,a,t,Lewis,
     &                           th_file,sp_file
!
! flag_a7=1 : NASA seven coefficient polynomial
! flag_a7=0 : Not support for surface reaction
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      real(8),parameter :: undef=-huge(1.d0)
      real(8) :: hk(5),tt,h00
      integer :: i,ios=0,iosf=0,ierr1,s,s2,nsno,nofind
      logical :: nml_comp_eq,nml_comp_gt,Lel,existth=.false.
!
      character(lenspc) :: spienm
      character(1)      :: percnt
      integer,parameter :: lc=7
!
      if(ncomp_suf==0) return
      ierror=0
      ns=ncomp+ncomp_suf
      iset=0
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
!-< 1. allocate array >-
!------------------------------------------------
! --- judgment of coefficient set "a7" or "acpk"
!------------------------------------------------
      tref=undef                 ! temperature [k]
      flag_a7=0
      fa7=flag_a7
      th_file = ""
      sp_file = ""
      thfile = ""
      spfile = ""
      nsno=0
!
      rewind ifli
      do
        read(ifli,surface_species,iostat=ios)
        if(ios<0) then
          exit
        elseif(ios==0) then
          nsno=nsno+1
          if(nsno==1) then
            fa7=flag_a7
            if(fa7==1) then
              if(th_file/=' ') then
                existth=.false.
                thfile=trim(th_file)
                inquire(file=thfile,exist=existth)
                if(.NOT.existth) then
                  write(ifle,'(2a)') 
     &            'ERR: Cannot found thermal file: ',trim(thfile)
                  write(ifle,'(3a)') 
     &            'MSG: Copy thermal file  ',trim(thfile),
     &             'into current directory'
                  call FFRABORT(1,'NO thermal file')
                endif
              else
                call FFRABORT
     & (1,'ERR: [th_file] file name is NOT defined in &surface_species')
              endif
              if(sp_file/='') then
                spfile=trim(sp_file)
                existth=.false.
                inquire(file=spfile,exist=existth)
                if(.NOT.existth) then
                  write(ifle,'(2a)') 
     &            'ERR: Cannot found species file : ',trim(spfile)
                  write(ifle,'(2a)') 
     &            'MSG: Copy thermal file [species.dat] ',
     &             'into current directory'
                  call FFRABORT(1,'NO [species.dat] file')
                endif
              endif
            endif
          endif
          if(.not.(flag_a7==1.or.flag_a7==0)) then
            write(ifle,*) '### ERR: all flag_a7 should be 0 or 1'
            call FFRABORT(1,' Reading error in [&surface_species]')
          endif
          if(flag_a7/=fa7) then
            write(ifle,*)
     &     '### ERR: all flag_a7 should be 1 in &surface_species'
            write(ifle,*)
     &     '### MSG: flag_a7=1 : NASA seven coefficient polynomial '
            call FFRABORT(1,' Reading error in [&surface_species]')
          endif
        else
          call FFRABORT(1,' Reading error in [&surface_species]')
        endif
      enddo
      if(nsno==0) then
        
        return
      endif
      if(nsno.ne.ncomp_suf) then
        write(ifle,*) ' ### ERR: ncomp_suf /= nsno'
        write(ifle,*) ' ### MSG: ncomp_suf= ',ncomp_suf,' in file:',
     &                           cntlnam
        write(ifle,*) ' ### MSG: [&surface_species] number : ',nsno
        write(ifle,*) ' ### reading error in "&surface_species"'
        write(ifle,*) ' ### MSG: Re-run prefflow '
        call FFRABORT(1,' Reading error in [&surface_species]')
      endif
      ns_suf=nsno
!
      if(.NOT.sw) then
        write(ifle,'(2a)') 
     &  ' ### ERR: NASA-a7-form MUST be used in ',
     &   '&Surface_species and &species'
        call FFRABORT(1,'ERR: Surface_species MUST be flag_a7=1')
      else
        sw=.true.                        ! for "a7"
      endif
      if(fa7==0) call 
     %   FFRABORT(1,'ERR: Surface_species MUST be flag_a7=1')
!----------------------------------------------------
! --- reading species name, molecular weight [kg/mol]
!----------------------------------------------------
      rewind ifli
      do s=ns_gas+1,ns_gas+ns_suf
        name = ' '                     ! name of species
        weight = undef                 ! molecular weight [kg/mol]
        read(ifli,surface_species,iostat=ios)
        call reading_judge_suf         ! reading and judge of data file
        if(name==' ') then             ! check "name", name of species
          write(ifle,*) '### ERR: read in [&surface_species]'
          write(ifle,*) 'lack of data "name" in no.',s
          call FFRABORT(1,'ERR: read in [&surface_species]')
        endif
        spcnam(s)=name
	spinam(s)=name
        if(s>=2) then
          do s2=ns_gas+1,s-1
            if( spcnam(s)==spcnam(s2) ) then
              write(ifle,*) '### reading error in "surface_species"'
              write(ifle,*) 'duplicated data : name(i)=name(j)=',
     &                       spcnam(s2)
              write(ifle,*) 'i,j =',s2,s
              call 
     &        FFRABORT(1,'ERR: name duplicated in [&surface_species]')
            endif
          enddo
        endif
!------------------------------------------------
! --- check "weight", molecular weight [kg/mol]
!------------------------------------------------
        if(.not.nml_comp_eq("weight",weight,undef,s,ifle,ierror)) then
          if(nml_comp_gt("weight",weight,0.0d0,s,ifle,ierror)) then
            wm(s) = weight              ! molecular weight [kg/mol]
            r_wm(s) = 1.0d0/wm(s)       ! reciprocal of 'wm' [mol/kg]
            Ri(s) = gascns*r_wm(s)      ! 'gascns'/'wm' [J/(kg K)]
          endif
        endif
      enddo
!
      if( s-1/=ns ) then                ! check number of species
        write(ifle,*) '### error : data error'
        write(ifle,*) 'no. of namelists (species) =',s-1
        write(ifle,*) 'it must be =',ns
        ierror = 1
      endif
!
!----------------------------------------------------------
! --- reading coefficient sets and display reading result
!----------------------------------------------------------
!
!----------------
! --- read again
!----------------
!
      if(sw) then 
        rewind ifli
        do s=ns_gas+1,ns_gas+ns_suf
        a=undef                     ! coefficient
        t=0                         ! temperature [K]
        tref=undef
        read(ifli,surface_species,iostat=ios)
        call reading_judge_suf      ! reading and judge of data file
        nofind=0
        if(a(1,1)==undef) then
          spienm=spcnam(s)
          do i=1,lenspc
          percnt=spcnam(s)(i:i)
          if(percnt=='%') then
            spienm=trim(spcnam(s)(1:i-1))
            exit
          endif
          enddo
          if(th_file/=' ') then
            call from_files(ifle,s,a,t,LEN_TRIM(th_file),
     &            LEN_TRIM(sp_file),thfile,spfile,spienm,nofind)
          else
            write(ifle,'(a)') 
     &  'ERR: Not defined thermal file in first [&surface_species]'
            write(ifle,'(a,/,a)') 
     &  'MSG: flag_a7=1 => NASA thermal file [therm.dat] will be used', 
     &  '                  or input NASA 7 coefficients and t(1~3)'
            write(ifle,'(a)')
     &  'MSG: flag_a7=1 is only be used for [&surface_species]'
          endif
          if(nofind==1.and.a(1,1)==undef) then
            write(ifle,'(2X,4a)') 'ERR: NOT found Species name: ',
     &        trim(spcnam(s)),' in thermal file: ',trim(thfile)
            write(ifle,'(2X,2a)')
     &             "WRN: Please input a7 coeffient in ",trim(cntlnam)
            call 
     &        FFRABORT(1,'ERR: NOT found species name in thermal file')
          else
            call read_a7(s,ns,ifle,a,t,tref,undef,a7,t3,ierror)
            if(my_rank==0) then
              write(ifle,'(2x,4a)') 
     &        'MSG: Tref, t(1:3), a(1:14) read from ',
     &        trim(thfile), ' for species: ',trim(name)
            endif
          endif
        else
          call read_a7(s,ns,ifle,a,t,tref,undef,a7,t3,ierror)
          if(my_rank==0) then
            write(ifle,'(2x,4a)') 
     &      'MSG: Tref, t(1:3), a(1:14) read from ',
     &      trim(cntlnam), ' for species: ',trim(name)
          endif
        endif
!
        tref=dble(tref)
        href(s)=h_t(s,tref) ![---]
        tref_comp(s)=tref
        enddo
!
        do s=ns_gas+1,ns_gas+ns_suf
          hform(s)=href(s)*Ri(s)    ![J/kg]
        enddo
        judg_t2(1)=.true.             ! every t3(2,*) is same
        do s=ns_gas+2,ns_gas+ns_suf
        if(t3(2,1)/=t3(2,s)) judg_t2(1)=.false.
        enddo
        call display_a7                    ! result of reading
      else  !not support for a5-model
!--------------
! --- a5 model
!--------------
      endif
!
! --- Lewis number
!
      Lel=.false.
      do i=1,nflud
        if(isthr(i)==3) then
           Lel=.true.
        endif
      enddo
      if(Lel) then
!        call memory_error
        Le(ns_gas+1:ns_gas+ns_suf) = 0.d0
        r_Le(ns_gas+1:ns_gas+ns_suf) = 0.d0
        rewind ifli
        do s=ns_gas+1,ns_gas+ns_suf
          Lewis = undef
          read(ifli,surface_species,iostat=ios)
          call reading_judge_suf
          if(.not.nml_comp_eq("Lewis",Lewis,undef,s,ifle,ierror))
     &    then
            if( nml_comp_gt("Lewis",Lewis,0.0d0,s,ifle,ierror) )then
              Le(s) = Lewis
              r_Le(s) = 1.d0/Le(s)
            endif
          endif
        enddo
        call disp_Lewis   ! result of reading
      endif
!----------------------
! --- final error check
!----------------------
      if(ierror/=0) then
        write(ifle,*) "ERR: total error : ",ierror,modnam,subnam
        call FFRABORT(1,'ERR: reading [&surface_species]')
      endif
!
      return
!
 9999 continue
      write(ifle,*) modnam,subnam
      ierror=1
      return
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine reading_judge_suf        ! reading and judge of data file
!=================================================
      if(s/=(ns_gas+ns_suf).and.ios<0) then
        write(ifle,*) '### reading error in "surface_species"'
        ierror = 1+ierror
      endif
      call nml_errmsg0(ifle,ios,'surface_species',ierr1)
      if( ierr1/=0 ) ierror = 1+ierror
      end subroutine reading_judge_suf
!
!=================================================
      subroutine display_a7             ! display of reading result
!=================================================
      if(my_rank==0) then
        write(ifll,*)
        write(ifll,5001) 'Surface-Phase Species List: '
        write(ifll,"(2x,a,i3)") 
     &  "number of chemical species on SITE and BULK: ",
     &  ns_suf
 5001   format(2x,a,80('='))
       endif
      do s=ns_gas+1,ns_gas+ns_suf
        if(my_rank==0) write(ifll,'(i3,1x,a20)') s,spcnam(s)
                          ! name of species
        if(my_rank==0) 
     & write(ifll,'(13x,a,e13.5,a)') 'mi =',wm(s),' [kg/mol]'
                          ! molecular weight [kg/mol]
        if(my_rank==0) 
     & write(ifll,'(13x,a,e13.5,a)') 'Ri =',Ri(s),' [J/(K kg)]'
                          ! 'gascns'/'wm' [J/(kg K)]
        if(my_rank==0) 
     &  write(ifll,'(13x,a,3f8.1,a)') 'T =',(t3(i,s),i=1,3),' [K]'
                          ! point of temperature region [K]
        if(my_rank==0) write(ifll,'(13x,a,$)') 'a ='
        if(my_rank==0) write(ifll,'(7e13.5)') (a7(i,1,s),i=1,7)
                          ! coefficient set of upper region
!        if(my_rank==0) write(ifll,*)
        if(my_rank==0) write(ifll,'(16x,7e13.5)') (a7(i,2,s),i=1,7)
                          ! coefficient set of lower retion
        if(my_rank==0) write(ifll,*)
        if(my_rank==0) 
     &   write(ifll,'(13x,a,E13.5,1x,a,3X,E13.5,1x,a)') 
     &   'Formation enthalpy : ',hform(s),
     &   '[J/kg]',href(s)/4.186d0,'[kcal/mole]'
        if(my_rank==0) write(ifll,*)
      enddo
      if(my_rank==0) then
      if(judg_t2(1))
     &  write(ifll,"(2X,5('-'),a)") "every T[2] is same value!"
      endif
      end subroutine display_a7
!
!-------------------------------------------------
! --- reading and check data
!=================================================
      subroutine read_a7(s,ns,f,a,t,tref,undef,a7,t3,e)
!=================================================
      integer,intent(in)  :: s
                  ! classification number of species
      integer,intent(in)  :: ns
                  ! number of species
      integer,intent(in)  :: f
                  ! file number to write error script
      real(8), intent(in)  :: a(7,2)
                  ! coefficient set
      real(8), intent(in)  :: t(3),tref
                  ! temperature [K]
                  ! (1:upper limit, 2:branch point, 3:lower limit)
      real(8), intent(in)  :: undef
                  ! judgment value
      real(8), intent(out) :: a7(7,2,ns)
                  ! coefficient sets
      real(8), intent(out) :: t3(3,ns)
                  ! temperature [K]
                  ! (1:upper limit, 2:branch point, 3:lower limit)
      integer,intent(out) :: e
                  ! total number of error
      integer :: i,j
      real(8)  :: tt
!--------------------------------
! --- check "a7", coefficient set
!--------------------------------
      A7_0:do i=1,2
        do j=1,7
          if(.not.nml_comp_eq("a",a(j,i),undef,s,f,e)) then
            a7(j,i,s)=a(j,i)
          else
            write(f,*) '"a" needs 14 values'
            exit A7_0
          endif
        enddo
      enddo A7_0
!-------------------------------------------------
! --- check "t3", point of temperature region [K]
! --- t(1):upper limit, t(2):branch point, t(3):lower limit
!-------------------------------------------------
      if( t(1)>t(2) .and. t(2)>t(3) .and. t(3)>0.d0 ) then
        t3(1:3,s)=t(1:3)
      else
        write(f,*) '### reading error in "species"'
        write(f,*) 'check "t" of no.',s
        e=1+e
        call FFRABORT(1,'ERR: t(1:3) reading erorr')
      endif
      if( .not.nml_comp_eq("tref",tref,undef,s,f,e) ) then
        if( nml_comp_gt("tref",tref,0.0d0,s,f,e) )
     &    tt = tref
      endif
      if(e/=0) then 
        call FFRABORT(1,'ERR: t(1:3) reading erorr')
      endif
      end subroutine read_a7
!
!=================================================
      subroutine disp_Lewis
!=================================================
      if(my_rank==0) write(ifll,'(2X,a)') "Lewis number of species:"
      if(my_rank==0) then
         do s=ns_gas+1,ns_gas+ns_suf
         write(ifll,'(2x,a8,2X,f8.3)') (spcnam(s)), Le(s)
         enddo
      endif
      if(my_rank==0) write(ifll,*)
      end subroutine disp_Lewis
!
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
!------------------------------------------------------------------------
! --- judgment of temperature region ( 1:higher region, 2:lower region )
!=================================================
      integer function judg_temp(ICOM,t)
!=================================================
      implicit none
        integer,intent(in) :: ICOM  ! classification number of species
        real(8), intent(in) :: t    ! temperature [K]
        if(t>t3(2,ICOM)) then       ! split by temperature
          judg_temp=1               ! higher
        else
          judg_temp=2               ! lower
        endif
      end function judg_temp
!
!=================================================
      real(8) function c_pi(t,ICOM)  ! specific heat of species [J/(kg K)]
!=================================================
      implicit none
      integer,intent(in) :: ICOM   ! classification number of species
      real(8), intent(in) :: t     ! temperature [K]
      integer rj                   ! result of judgment
!
      rj = judg_temp(ICOM,t)
      c_pi =a7(1,rj,ICOM)
     &    +t*(a7(2,rj,ICOM)
     &    +t*(a7(3,rj,ICOM)
     &    +t*(a7(4,rj,ICOM)
     &    +t* a7(5,rj,ICOM))))
      c_pi =c_pi*Ri(ICOM)*dble(act(icom))
!
      end function c_pi
!
!=================================================
      real(8) function c_p(ncomp,t,yi)    ! specific heat of mixture gas
                                          ! [J/(kg K)]
!=================================================
      implicit none
        integer,intent(in) :: ncomp       ! number of species
        real(8), intent(in) :: t,yi(ncomp) ! temperature [K],
                                          ! mass fraction of species
        integer icom                      ! classification number of species
        
        c_p=0.d0
        do icom=1,ncomp
          c_p=c_p+yi(icom)*c_pi(t,icom)*dble(act(icom))
        enddo
      end function c_p
!
!=================================================
      real(8) function h_t(icom,t)       ! (enthalpy of species)/(universal
                                         !  gas constant) [K]
!=================================================
      implicit none
        integer,intent(in) :: icom      ! classification number of species
        real(8), intent(in) :: t        ! temperature [K]
        integer rj                      ! result of judgment
!
        rj = judg_temp(icom,t)
        h_t =a7(6,rj,icom)
     &   +t*(a7(1,rj,icom)
     &   +t*(a7(2,rj,icom)*0.5d0
     &   +t*(a7(3,rj,icom)*1.d0/3.d0
     &   +t*(a7(4,rj,icom)*0.25d0
     &   +t*(a7(5,rj,icom)*0.2d0)))))
        h_t=h_t*dble(act(icom))
      end function h_t
!
!=================================================
      real(8) function h_tr(ns_gas,icom,t,tc,pc)       
! (enthalpy of species)/(universal
!  gas constant) [K]
!=================================================
      implicit none
      integer,intent(in)   :: ns_gas,icom   
      real(8),intent(in)   :: t
      real(8),intent(out)  :: tc,pc
      integer :: rj,NR                   ! result of judgment
!
      real(8) :: yi(mcomp),dens,VVV,delt_cp,delt_H,delt_C,DVDP_R
!
      yi(1:mcomp)=0.d0
      yi(icom)=1.d0
      NR=2
      rj = judg_temp(icom,t)

      h_tr =a7(6,rj,icom)
     &   +t*(a7(1,rj,icom)
     &   +t*(a7(2,rj,icom)*0.5d0
     &   +t*(a7(3,rj,icom)*1.d0/3.d0
     &   +t*(a7(4,rj,icom)*0.25d0
     &   +t*(a7(5,rj,icom)*0.2d0)))))!+delt_H/Ri(icom)
      h_tr=h_tr*dble(act(icom))
      end function h_tr
!=================================================
      subroutine dhr(icom,TTT,PPP,delt_H,delt_cp)
!=================================================
      implicit none
      integer,intent(in)   :: icom
      real(8), intent(in)  :: TTT
      real(8), intent(out) :: delt_H,delt_cp
      real(8), intent(inout)  :: PPP
!
      real(8) :: MIX_W,AAA,lll,Aij0,BBB,alph,alph_d,alph_dd,AAA_D,
     &           AAA_DD,a_alph,Cij,Cji,A(4),RT1(2),RT2(2),RT3(2),
     &           DVDT,DVDP_R,dum1,dum2,dum3,a1,a2,a3,Tcc
      integer :: NR
!
      MIX_W=gascns*r_wm(icom)
      BBB=0.d0
      if(PPP<0.d0) PPP=Pc(icom)
      lll=p3+p4*omga(icom)-p5*omga(icom)*omga(icom)
      Aij0=p1*Ri(icom)**2*Tc(icom)**2/(Pc(icom))
      BBB=p2*Ri(icom)*Tc(icom)/(Pc(icom))
      alph=Aij0*(1.d0+lll*(1.d0-dsqrt(TTT/Tc(icom))))**2
      alph_d=-Aij0*lll/Tc(icom)*((1.d0+lll)
     &             *(TTT/Tc(icom))**(-0.5)-lll)
      alph_dd=0.5d0*Aij0*lll/Tc(icom)**2*(lll+1.d0)
     &             *(TTT/Tc(icom))**(-1.5d0)
      AAA=0.d0
      AAA_D=0.d0
      AAA_DD=0.d0
      a_alph=dsqrt(Aij0*Aij0*alph*alph)
      AAA=a_alph
      AAA_D=0.5d0
     &     *(alph*alph_d+alph_d*alph)
     &     *(Aij0*Aij0/(alph*alph))**0.5d0

      Cij=(alph_d*alph-alph*alph_d)
     &           /alph**2
      Cji=(alph_d*alph-alph*alph_d)
     &           /alph**2
      AAA_DD=0.25d0
     &        *(Aij0*Aij0/(alph*alph))**0.5d0
     & *(alph*(Cij*alph_d+2.d0*alph_dd)+alph*(Cji*alph_d+2.d0*alph_dd))
      a1=-MIX_W*TTT/PPP
      a2=AAA/PPP-BBB**2-MIX_W*TTT/PPP*BBB
      a3=-AAA*BBB/PPP
!
      dum1=0.d0
      dum2=0.d0
      dum3=0.d0
      CALL CUBIC(1.D0,A1,A2,A3,dum1,dum2,dum3,NR)
      VVV = DMAX1(dum1,dum2,dum3)   
      IF((NR.EQ.3).AND.(TTT.LE.Tc(icom))) THEN
        VVV = DMIN1(dum1,dum2,dum3)
      ENDIF

!
!      A(1)=1.d0
!      A(2)=-MIX_W*TTT/PPP
!      A(3)=AAA/PPP-BBB**2-MIX_W*TTT/PPP*BBB
!      A(4)=-AAA*BBB/PPP
!      dum1=1.d-10
!      CALL CUBIC_Nu(A,RT1,RT2,RT3,dum1,NR) 
!      VVV=RT1(1)
      dum1=-MIX_W/PPP
      dum2=AAA_D/PPP-MIX_W*BBB/PPP
      dum3=-BBB*AAA_D/PPP
      DVDT=-(dum1*VVV**2+dum2*VVV+dum3)/(3.d0*VVV**2+2.d0*a1*VVV+a2)
!
      dum1=MIX_W*TTT/PPP**2
      dum2=(BBB*MIX_W*TTT-AAA)/PPP**2
      dum3=BBB*AAA/PPP**2
      DVDP_R=-(3.d0*VVV**2+2.d0*a1*VVV+a2)/(dum1*VVV**2+dum2*VVV+dum3)
      delt_H=
     &       PPP*VVV
     &      -MIX_W*TTT
     &      +(AAA_D*TTT-AAA)*LOG((VVV+BBB)/VVV)/BBB 
      delt_cp=
     &        PPP*DVDT
     &       -MIX_W
     &       +TTT*AAA_DD*dlog((VVV+BBB)/VVV)/BBB
     &       -(AAA_D*TTT-AAA)*(1.d0/(VVV*(VVV+BBB)))*DVDT
      end subroutine dhr
!=================================================
!=================================================
      real(8) function enthalpy_r(ncomp,t,tr,yi)
!=================================================
        implicit none
        integer,intent(in) :: ncomp        ! number of species
        real(8), intent(in) :: t,tr,yi(ncomp)  ! temperature [K],
                                           ! mass fraction of species
        integer icom                       ! classification number of species
        real(8) :: dum1
        enthalpy_r=0.d0
        do icom=1,ncomp
        enthalpy_r=enthalpy_r
     &  +yi(icom)*((h_t(icom,t)-h_t(icom,tr))*Ri(icom))*dble(act(icom))
        enddo
      end function enthalpy_r
!=================================================
      real(8) function enthalpy(ncomp,t,yi) ! enthalpy of mixture gas [J/kg]
!=================================================
        implicit none
        integer,intent(in) :: ncomp        ! number of species
        real(8), intent(in) :: t,yi(ncomp) ! temperature [K],
                                           ! mass fraction of species
        integer icom                       ! classification number of species
        enthalpy=0.d0
        do icom=1,ncomp
          enthalpy=enthalpy
     &    +yi(icom)*(h_t(icom,t)*Ri(icom))*dble(act(icom))
        enddo
      end function enthalpy
!
!
!=================================================
      real(8) function h_ref(icom,t)   ! (enthalpy of species)/(universal
                                       !  gas constant) [K]
!=================================================
      implicit none
        integer,intent(in) :: icom    ! classification number of species
        real(8), intent(in) :: t      ! temperature [K]
        integer :: rj                 ! result of judgment
        real(8) :: ggg
!
        rj = judg_temp(icom,t)
        h_ref =a7(6,rj,icom)
     &   +t*(a7(1,rj,icom)
     &   +t*(a7(2,rj,icom)*0.5d0
     &   +t*(a7(3,rj,icom)*1.d0/3.d0
     &   +t*(a7(4,rj,icom)*0.25d0
     &   +t*(a7(5,rj,icom)*0.2d0)))))
        h_ref=(h_ref-href(icom))*dble(act(icom))
!        
      end function h_ref
!
!
!enthalpy_ref=[J/kg]
!=================================================
      real(8) function enthalpy_ref(ncomp,t,yi) ! enthalpy of mixture gas [J/kg]
!=================================================
        implicit none
        integer,intent(in) :: ncomp        ! number of species
        real(8), intent(in) :: t,yi(ncomp) ! temperature [K],
                                           ! mass fraction of species
        integer :: icom                    ! classification number of species
        real(8) :: hreff,tr
!
        enthalpy_ref=0.d0
        hreff=0.d0
        do icom=1,ncomp
        enthalpy_ref=enthalpy_ref
     &      +yi(icom)*(h_ref(icom,t)*Ri(icom))*dble(act(icom))
!
        enddo
      end function enthalpy_ref
!
!=================================================
      real(8) function s_0(icom,t)       ! (entropy of species)/
                                        ! (universal gas constant) [-]
!=================================================
      implicit none
        integer,intent(in) :: icom      ! classification number of species
        real(8), intent(in) :: t         ! temperature [K]
        integer rj   
! result of judgment
!        s_0=0.d0
!        if(act(icom)==0) return 
        rj = judg_temp(icom,t)
        s_0 = a7(7,rj,icom)
     &   +t*( a7(2,rj,icom)
     &   +t*( a7(3,rj,icom)*0.5d0
     &   +t*( a7(4,rj,icom)*0.33333d0
     &   +t*( a7(5,rj,icom)*0.25d0 ))))
     &   +    a7(1,rj,icom)*log(t)
        s_0 = s_0*dble(act(icom))
      end function s_0
!
!=================================================
      real(8) function entropy(ncomp,t,yi) ! entropy of mixture gas [J/(kg K)]
!=================================================
      implicit none
        integer,intent(in) :: ncomp       ! number of species
        real(8), intent(in) :: t,yi(ncomp) ! temperature [K],
                                          ! mass fraction of species
        integer icom                      ! classification number of species
        entropy = 0.d0
        do icom=1,ncomp
          entropy=entropy+yi(icom)*s_0(icom,t)*Ri(icom)*dble(act(icom))
        enddo
      end function entropy
!
!==========================================================
      real(8) function sigma_yr(ncomp,yi)   ! [J/(kg K)]
!==========================================================
      implicit none
        integer,intent(in) :: ncomp      ! number of species
        real(8), intent(in) :: yi(ncomp) ! mass fraction of species
        integer :: icom                  ! classification number of species
        sigma_yr=0.d0
        do icom=1,ncomp
          sigma_yr=sigma_yr+yi(icom)*Ri(icom)*dble(act(icom))
        enddo
      end function sigma_yr
!
      
C***********************************************************************
      SUBROUTINE CUBIC(A3,A2,A1,A0,Y1,Y2,Y3,NR)
C===============================================================
C     THIS SUBROUTINE FINDS THE REAL ROOTS OF THE CUBIC EQUATION
C          A3*Y**3+A2*Y**2+A1*Y*A0=0
C     REAL ROOTS ARE STORED IN Y1, Y2, Y3
C     NR IS THE NUMBER OF REAL ROOTS
C
C     PROGRAMMER - M. D. SCHUMAN - 5/18/70
C     MODIFICATION -  S. TSUDA   -11/08/07
C
C===============================================================
C+++++ DEFINE VARIABLES
      IMPLICIT NONE 
!
! --- [dummy arguments]
!
      real*8 ,intent(in)    :: A3,A2,A1,A0
      real*8 ,intent(out)   :: Y1,Y2,Y3
      integer,intent(out)   :: NR  
!

C----- I/O PARAMETERS
!      REAL*8 A1,A2,A3,A0
!      REAL*8 Y1,Y2,Y3

C----- LOCAL VARIABLES
      REAL*8 :: DUM,DUM2,A,B,X1,X2,X3,C,D,Q,R,Z 
C---------------------------------------------------------------
CHAPTER1 INITIAL SETTING
C---------------------------------------------------------------
      Y1 = 0.0D+0
      Y2 = 0.0D+0
      Y3 = 0.0D+0
C---------------------------------------------------------------
CHAPTER2 SOLVING CUBIC EQ.
C---------------------------------------------------------------
      IF(A3.NE.0.E0) THEN
         A  = A3
         B  = A2/3.0D+0
         X1 = B
         X2 = B
         X3 = B
         C  = A1/3.0D+0
         D  = A0
C+++++ PARAMETERS FOR DISCRIMINANT
         Q  = A*C-B**2                        
         R  = 0.5D+0*(3.0D+0*A*B*C-D*A**2)-B**3 
C+++++ FOR Q <= 0.0D+0; THREE ROOT
         IF(Q.LE.0.0D+0) THEN
            IF((Q**3+R**2).LE.0.0D+0) THEN
               Z    = DABS(R)/(DSQRT(-Q**3))
               DUM  = DATAN(DSQRT(1.0D+0-Z**2)/Z)/3.0D+0
               DUM2 = R/DABS(R)*2.0D+0*DSQRT(-Q)
C
               X1   = DUM2*DCOS(DUM)
               X2   = DUM2*DCOS(DUM+2.0944D+0)
               X3   = DUM2*DCOS(DUM+4.1888D+0)
               NR   = 3
            ELSE
               DUM  = 1.0D+0
               IF(DABS(R).GT.0.0D+0) DUM=(-Q)**3/(R*R)
               IF(DUM.GE.0.005D+0) THEN
                  DUM = DABS(R)/DSQRT(-Q**3)
                  Z   = DLOG(DUM+DSQRT(DUM**2-1.0D+0))/3.0D+0
                  X1  = R/DABS(R)*DSQRT(-Q)*(DEXP(Z)+DEXP(-Z))
                  NR  = 1
               ELSE
                  X1  = R/((DABS(R)**1.3333D+0*2.0D+0**0.3333D+0))
     &                 *((4.0D+0*R*R)**0.3333D+0-Q)
               ENDIF
               NR = 1
            ENDIF
C+++++ FOR Q > 0.0D+0; ONLY ONE ROOT 
         ELSE
            DUM = DABS(R)/DSQRT(Q**3)
            Z   = DLOG(DUM+DSQRT(DUM**2+1.0D+0))/3.0D+0
            X1  = R/DABS(R)*DSQRT(Q)*(DEXP(Z)-DEXP(-Z))
            NR  = 1
         ENDIF
C+++++ SETTING SOLUTIONS
         Y1 = (X1-B)/A
         Y2 = (X2-B)/A
         Y3 = (X3-B)/A
C-----------------------------------------------------------------------
CHAPTER3 SOLVING EQUATION EXCEPT FOR CUBIC EQUATION
C-----------------------------------------------------------------------
      ELSE
C+++++ QUADRATIC EQUATION; A3 = 0
         IF(A2.NE.0.0D+0) THEN
            DUM = A1**2-4.0D+0*A2*A0 ! DISCRIMINANT FORM
            NR  = 0
            IF(DUM.GE.0.0D+0) THEN
               Y1 = (-A1+DSQRT(DUM))/(2.0D+0*A2)
               Y2 = (-A1-DSQRT(DUM))/(2.0D+0*A2)
               NR = 2
            ENDIF
C+++++ LINEAR EQUATION; A3 = A2 = 0
            ELSE
               Y1 = -A0/A1       ! ONLY ONE SOLUTION
               NR = 1
            ENDIF
C
      ENDIF
C
C===== END OF OPERATION ===============================================
      RETURN
      END subroutine CUBIC
!
C***********************************************************************
      SUBROUTINE CUBIC_Nu(A,RT1,RT2,RT3,TOL,NOROOT)
C===============================================================
C     THIS SUBROUTINE FINDS THE REAL ROOTS OF THE CUBIC EQUATION
C          A3*Y**3+A2*Y**2+A1*Y+A0=0
C===============================================================
C+++++ DEFINE VARIABLES
      IMPLICIT NONE 
!
! --- [dummy arguments]
!
      real*8 ,intent(inout) :: A(4),TOL
      real*8 ,intent(out)   :: RT1(2),RT2(2),RT3(2)
      integer,intent(out)   :: NOROOT
!
C----- LOCAL VARIABLES
      REAL*8 :: ZERO,B(2),C(4),ZZ,Z,X,Y,ANG
      integer:: IT,I,N,NDA
!---------------------------------------------------------------
      ZERO=TOL/10.D0
!      IF(ABS(A(1))<ZERO) THEN
!        CALL QUADRT(a(2),RT1,RT2,RT3,TOL,NOROOT)
!        RETURN
!      ENDIF
!
      NOROOT=3
      IT=0
      DO 10 I=2,4
      A(I)=A(I)/A(1)
 10   ENDDO
      A(1)=1.D0
      IF(ABS(A(2))>ZERO) THEN
        NDA=3
        B(1)=1.D0
        B(2)=-A(2)*0.3333333D0
        CALL LINCNG(A,NDA,B,C)
        IT=1
      ELSE
        DO 20 I=1,4
        C(I)=A(I)
 20     ENDDO
      ENDIF
!
      X=(C(4)*C(4))*0.25+(C(3)**3)*0.037037D0
!
      IF(X>=0.D0) THEN
        X=SQRT(X)
        Y=-(C(4)*0.5D0)
        I=1
        RT1(I)=Y+X
 15     CONTINUE
        N=0
        IF(ABS(RT1(1))<0.D0)N=1
        IF(ABS(RT1(1))>ZERO) THEN
          RT1(I)=EXP((LOG(ABS(RT1(I))))*0.333333D0)
          IF(N==1) RT1(I)=-RT1(I) 
        ENDIF
        IF(I/=2) THEN
          I=2
          RT1(I)=Y-X
        ENDIF
        IF(I/=2) GOTO 15
        RT2(2)=((RT1(1)-RT1(2))*0.5D0)*1.732051D0
        RT1(1)=RT1(1)+RT1(2)
        RT1(2)=0.D0
        RT2(1)=-RT1(1)*0.5D0
        RT3(1)=RT2(1)
        RT3(2)=-RT2(2)
      ELSE
        ZZ=ABS(C(3))
        X=-(C(4)*0.5D0)/SQRT((ZZ**3)*0.037037D0)
        ANG=ACOS(X)
        Y=2.D0*SQRT(ZZ*0.333333333D0)
        ANG=ANG*0.333333333D0
        RT1(1)=Y*COS(ANG)
        RT2(1)=Y*COS(ANG+2.09440D0)
        RT3(1)=Y*COS(ANG+4.18879D0)
        RT1(2)=0.D0
        RT2(2)=0.D0
        RT3(2)=0.D0
      ENDIF
      IF(IT==0) RETURN
      RT1(1)=RT1(1)+B(2)
      RT2(1)=RT2(1)+B(2)
      RT3(1)=RT3(1)+B(2)
      RETURN
      END subroutine CUBIC_Nu
!
!***********************************************************************
      subroutine LINCNG(A,NDA,B,C)
!***********************************************************************
      IMPLICIT NONE 
!
      real*8 ,intent(inout) :: A(4),B(2),C(4)
      integer,intent(INout)   :: NDA
!----- LOCAL VARIABLES!
      integer:: IT,I,NZ
      real*8 :: B1PWR

      B1PWR=1.D0
      NZ=NDA+1
      DO 40 I=1,NZ
      C(I)=A(I)
 40   ENDDO
 6    CONTINUE
      DO 8 I=2,NZ
      C(I)=C(I)+C(I-1)*B(2)
 8    ENDDO
      C(NZ)=C(NZ)*B1PWR
      NZ=NZ-1
      B1PWR=B1PWR*B(1)
      IF(NZ<=1) THEN
      ELSE
        GOTO 6
      ENDIF
      C(1)=C(1)*B1PWR
      RETURN

      END subroutine LINCNG
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine rho_real_PR(IMODE,ICVL,ncomp,yi,TTT,PPP,dens,
     &                    VVV,delt_cp,delt_H,delt_C,DVDP_R)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_metrix,only  : XF=>rcomp
      implicit none 
!
      integer,intent(in)     :: ncomp,ICVL
      integer,intent(inout)  :: IMODE
      real*8 ,intent(inout)  :: TTT,PPP,dens
      real*8 ,intent(inout)  :: yi(ncomp)
      real*8 ,intent(out)    :: VVV,delt_cp,delt_H,delt_C,DVDP_R!delt_S
      
!---------------------  
! --- local 
!---------------------  
      integer :: icom,ii
      real*8           :: VVV0,
     &                    a_alph,Cij,Cji,
     &                    delt_S,MIX_W1,
     &                    MIX_W,AAA,BBB,AAA_D,AAA_DD,DVDT!,DVDP_R
!      
      real*8           :: a1,a2,a3,lll,dum1,dum2,dum3,func,dfunc
      integer          :: j,k,NR,i,itr=10
      real*8           :: tol=1.0d-8,A(4),RT1(2),RT2(2),RT3(2),Tcc
      real*8           :: dumm1,dumm2,dumm3,PPPC,TTTC
!-------------------------------------------------------------      
!
      end subroutine rho_real_PR 
!
      end module module_species 
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_chemreac
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_species,only : a7
      implicit none
      character(17),parameter,private :: modnam='(module_chemreac)'
      real*8,parameter :: SML=1.d-25
      real*8,parameter :: Faraday=96485.3415d0   !(C/mol)
!
! --- < module for chemical reaction formula >
!
      integer,save      :: ical_suf=0 ! ical_suf=1 : surface reaction
                                      ! ical_suf=2 : CVD wall defined
                                      !              but no surface reaction
!      integer,save      :: ical_P_R=0 ! =0 : NO particle reaction
!                                      ! =1 : particle reaction
!
      integer,save      :: iBV=0
      integer,save      :: IZ_O2=0,IZ_N2=0,IZ_NO=0,IZ_O=0
!
      integer,parameter :: igas=1,isurface=2
      integer,parameter :: falloff=1,activated=2
      integer,private   :: ns=0       ! number of species
      integer,private   :: ns_sf=0    ! number of species on SITE
      integer,private   :: ns_gs=0    ! number of species in GAS
!
      integer,parameter :: ovrall=1,  ! ovrall : overall reaction model(gas)
     &                     elemnt=2,  ! elemnt : elementary reaction model(gas/surface)
     &                     ebrkup=3,  ! ebrkup : eddy break up model(gas)
     &                     ovalfr=4,  ! ovrall : overall reaction model + CO and SOOT(gas)
     &                     userreac=5,! userreac:user reaction model(gas/surface)
     &                     stick_sf=6,! surface stick model
     &                     LH_ER=7,   ! surface
     &                     Bohm=8,    ! surface
     &                     elemarbi=9,! surface & Gas   (gas/surface)
     &                     BohmE=10,  ! surface by Ion-energy-dependent Bohm model
     &                     BohmY=11,  ! surface by Ion-enhanced Reac. Yield Bohm Model
     &                     Butler_Volmer=12,
     &                     Zeldovich=13,
     &                     oval_ebrkup=14,
     &                     SSFRRM=15
!
!     &                     P_CHAR_COM=15,
!     &                     P_CHAR_NO=16,
!     &                     P_HCN_O2_NO=17,
!     &                     P_HCN_NO_N2O2=18,
!     &                     P_prompt_NO=19,
!     &                     P_Zeld_NO=20
!
! BohmE=8  : Surface rection Bohm MODEL
! BohmE=10 : Ion-energy-dependent model
! BohmY=11 : Ion-enhanced Reac. Yield Model
! Butler_Volmer=12 : Butler_Volmer Model

      integer :: nq=0,nneq            ! number of equation
      integer,save,allocatable :: ireq(:)      ! flag of reaction model
      integer,save,allocatable :: M_3rd(:)     ! 
      integer,save,allocatable :: P_3rd(:,:)   ! 
      real(8),save,allocatable :: vreq(:,:,:)  ! stoichiometric coefficients
      real(8),save,allocatable :: preq(:,:,:)  ! parameters for reaction rate
      real(8),save,allocatable :: preq1(:,:)
      real(8),save,allocatable :: mreq(:,:)    ! coefficients of 3rd body (aki)
      real(8),save,allocatable :: creq(:,:)    ! index number of mol
!                                              ! concentration(overall)
      real(8),save,allocatable :: sub_vreq(:,:)! subtraction of "vreq" of each
!                                              ! species in each equation
      real(8),save,allocatable :: sigma_vreq(:)! change of "sub_vreq" in
                                               ! each equation
      real(8),save,allocatable :: vreq_a7(:,:,:)! summation of "vreq" * "a7"
                                                ! in each equation
      integer,save,allocatable :: kind_chem(:)  ! flag of if gas-phase or surface
      integer,save,allocatable :: kind_prs(:)   ! flag of pressure-dep. reaction
      integer,save,allocatable :: ig_iter(:)    ! ignition iter for every reaction
      integer,save,allocatable :: stp_iter(:)   ! 
      integer,save,allocatable :: UNIT_SI(:)    ! =1: cm-sec-mole-kcal-Kelvins
      real(8),save,allocatable :: const_R(:)
!
      real*8,save,allocatable :: lreq(:,:,:)    ! 
!                                                !parameters for low pressure(Lindemann)
      real*8,save,allocatable :: troe(:,:,:)    ! parameters for Troe form
      integer,save,allocatable :: l_sw(:,:)     ! flag of Lindemann form
!
      integer,private :: iset=0
      real(8),save,allocatable :: preq_LH_ER(:,:,:)
      real(8),save             :: m_LHER=2,k_LHER=1.d0
      integer,save,allocatable :: nosp_stick(:),nosp_Bohm(:)
      real(8),save,allocatable :: coefBohm(:)
      integer,save,allocatable :: Bohm_sub(:,:,:)
!
      real(8),save,allocatable :: BV_ref_I(:)
      real(8),save,allocatable :: BV_a(:,:)
      real(8),save,allocatable :: BV_para(:,:,:)
      real(8),save,allocatable :: BV_ref_molfrc(:,:,:)
!
      integer,save,allocatable :: LDCOM(:,:)
!///////////////////////////////////////////////////////////////////////
      contains
!---------------------------------------------------------
!< #1. check consistency of values in argument & module >-
!---------------------------------------------------------
!====================================
      subroutine chkncomp(ival,sbnm)
!====================================
      character(LEN=10),parameter :: subnam='(chkncomp)'
      integer     ,intent(in) :: ival
      character(*),intent(in) :: sbnm
      if(ival==ns) return
      call modutl_chkset('>',1,iset,modnam//subnam)
      write(*,*) '### program error -1- ',modnam//subnam
      write(*,*) 'sbnm=',sbnm
      stop
      end subroutine chkncomp
!---------------------------
!< #2. input name list >--
!---------------------------
!=======================================================
      subroutine inputdata(ifli,ifll,ifle,my_rank,mmcomp,lsuf,
     &           cntlnam,ncmpx,ncomp_sufx,
     &           lcomp,igT_iter,mneq,ncompall,
     &           chem_bc,blktyp,cofleft,net_vreq,coal_CHO,
     &           ireac,ierror)
!=======================================================
      use module_species,only : lenspc,spcnam,chk_spec=>chkncomp,
     &                          wm
!
      implicit none
!
!-------------------------
! --- [dummy arguments]
!-------------------------
      integer,intent(in)     :: ifli,ifll,ifle,ncmpx,my_rank,igT_iter,
     &                          mneq,ncompall,ncomp_sufx,ireac
      integer,intent(inout)  :: chem_bc(mneq)
      integer,intent(in)     :: mmcomp,blktyp(mmcomp)
      character(*),intent(in):: cntlnam
      logical,intent(in)     :: lcomp
      logical,intent(out)    :: lsuf
      integer,intent(out)    :: ierror,coal_CHO(mneq)
      real*8, intent(out)    :: cofleft(mneq*ncompall),
     &                          net_vreq(mneq*ncompall)
!
! --- [namelist]
!
      integer,parameter :: lc=8
      integer,parameter :: mcomp=200            ! number of species
      character(LEN=20) :: model
      character(lenspc) :: nleft(mcomp+1)        ! name of species on left side
      character(lenspc) :: nright(mcomp+1)       ! name of species on right side
      real(8)  :: cleft(mcomp+1),cright(mcomp+1) ! stoichiometric coefficient
      real(8)  :: pleft(mcomp+1),pright(mcomp+1) ! arbitrary power index
      real(8)  :: creac(mcomp+1)                 ! index number of mol-density
                                                 ! concentration(overall)
      real(8)  :: a,alpha,beta,e                 ! 1) overall reaction : 
                                                 ! k=A*T^alpha*P^beta*exp(-e/R/T)
                                                 ! 2) stick reaction :
                                                 ! si=A*B^alpha*exp(-e/R/T)
                                                 ! 3) user rection :
                                                 ! si=A*B^alpha*exp(-e/R/T)
      real(8)  :: af,alphaf,ef,ab,alphab,eb      ! elementary reaction : 
                                                 ! k = (af)*T^(alpahf)*EXP(-(ef)/(R*T))
      real*8  :: af_low,alphaf_low,ef_low,
     &           af_high,alphaf_high,ef_high     ! for Lindemann etc.
      real*8  :: troe_a(2),troe_t1(2),troe_t2(2),troe_t3(2) ! Troe form
      real(8)  :: M(mcomp)                       ! coefficients of 3rd body
!
      real(8)  :: a_LH_ER(mcomp),beta_LH_ER(mcomp),H_LH_ER(mcomp),
     &            L_LH_ER(mcomp),n_LH_ER(mcomp),m_LH_ER,k_LH_ER
!
      real(8)  :: coeff_Bohm
!
      real(8)  :: cr1,cr2,r
!
      real(8),parameter :: WM_CO=28.d-3,WM_C=12.d-3
      real(8)  :: Y_CO,nu_SOOT,nu_CO
      integer :: ignition_iter=0,stop_iter=10000000
      real(8) :: Y_SOOT,C_fuel,WM_fuel,nu_fuel
      character(lc) :: kind,prs_kind
      integer :: unit_A_e
      real(8) :: I_ref_BV,alp_a_BV,alp_c_BV,
     &           r_BV_left(mcomp+1),r_BV_right(mcomp+1),
     &           ref_mol_frc_left(mcomp+1),
     &           ref_mol_frc_right(mcomp+1)
!
      namelist /chemreac/ model,nleft,nright,cleft,cright,
     &                                       pleft,pright,
     &                    creac,a,alpha,beta,e,
     &                    Y_SOOT,C_fuel,nu_fuel,WM_fuel,
     &                    af,alphaf,ef,ab,alphab,eb,M,
     &                    cr1,cr2,r,ignition_iter,stop_iter,kind,
     &                    a_LH_ER,beta_LH_ER,H_LH_ER,L_LH_ER,
     &                    n_LH_ER,m_LH_ER,k_LH_ER,
     &                    coeff_Bohm,unit_A_e,
     &                    af_low,alphaf_low,ef_low,
     &                    af_high,alphaf_high,ef_high,
     &                    troe_a,troe_t1,troe_t2,troe_t3,
     &                    prs_kind,
     &                    I_ref_BV,
     &                    alp_a_BV,alp_c_BV,
     &                    r_BV_left,r_BV_right,
     &                    ref_mol_frc_left,ref_mol_frc_right
!
! --- [local entities]
!
!    unit_A_e=0 : SI unit: A=m^3*(n-1)/[mole^(n-1)*s*K^(-alpha)] ; E=J/mole 
!    unit_A_e=1 : cm unit: A=cm^3*(n-1)/[mole^(n-1)*s*K^(-alpha)]; E=cal/mole 
!    unit_A_e=2 : SI unit: m^3-kmol-s-J-K, input (E/R) for [ef] 
!
!
      character(LEN=11),parameter :: subnam='inputdata'
      integer,parameter  :: iundef=-huge(1)
      real(8) ,parameter :: undef=-huge(1.d0)
!
      character(LEN=11),parameter :: mdllst(15)=(/
     &             'overall    ','elementary ','eddybreakup',
     &             'fire       ','user       ','stick      ',
     &             'LH_ER      ','BOHM       ','arbi_elem  ',
     &             'E_BOHM     ','Y_BOHM     ','BV         ',
     &             'Zeldovich  ','eddy_over  ','SSFRRM     '
!
!     &             'P_CHAR_COMB','P_CHAR_NO  ','P_HCN_O2_NO',
!     &             'P_NO_N2O2  ','P_prompt_NO','P_Zeld_NO'
!
     &              /)
!
      character(LEN=9),parameter :: prslst(2)=(/'fall-off ',
     &                                          'activated'/)
!
!----------------------------------------------------------------
!
! overall        : gas
! elementary     : gas/surface
! eddybreakup    : gas
! fire           : gas
! user           : gas/surface
! stick          : surface
! LH_ER          : surface
! Bohm           : surface
! Zeldovich
!
!----------------------------------------------------------------
!
      character(lc),parameter :: 
     &    chemkd(2)=(/'gas     ','surface '/)
      character(lc),parameter :: gas='gas',surface='surface'
!
      integer :: ispc(mcomp)
      integer :: i,j,k,ios,imodl,s,kd,is,idum,idum1
      integer :: q,ierr1
      real(8) :: dd2,dum,dum1,dum2,dum3
      integer nol,nor              ! number of species in each side
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
!----------------------------------
! --- ireac
!----------------------------------
      if(ireac==0) return

!
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
      if(mcomp<ncmpx+ncomp_sufx) 
     & call FFRABORT(1,'ERR: mcomp<ncomp,ncomp_suf')
!--------------------
!-< 1. Initial set >-
!--------------------
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
!--------------------------
! --- check reaction model
!--------------------------
!        if(TRIM(adjustl(kind))==gas) then
!        elseif(TRIM(adjustl(kind))==surface) then
!          if(TRIM(adjustl(model))==)
!        else
!        endif
      enddo
!
      nneq=nq
!      if(nneq.gt.1) then
      if(nneq>0) then
        if(.not.lcomp) then
          write(ifle,*) ' ### Warning: Chemical reaction must use ',
     &                     'compressible or low-Mach approximation'
          call FFRABORT(1,' Change [flow] in [&model]')
        endif
      endif
      if(nq<1) return                ! no reaction flow
      if(ierror.ne.0) goto 9999
!---------------------------
!--< 1.2 check interface >--
!---------------------------
      ns=ncmpx+ncomp_sufx
      ns_gs=ncmpx
      ns_sf=ncomp_sufx
      call nml_chksiz(ifle,'mcomp',mcomp,ns,    ! if ns>mcomp, ierr1=1
     &                    modnam//subnam,ierr1)
      if(ierr1/=0) then
        ierror=1
        goto 9999
      endif
      call chk_spec(ns_gs,modnam//subnam)
!---------------------------
!--< 1.3 allocate arrays >--
!---------------------------
      allocate( ireq(     nq),    ! flag number of reaction model
     &         M_3rd(     nq),    ! 3rd flag
     &         P_3rd(   nq,2),    ! pressure 3rd flag
     &          vreq(ns,2,nq),    ! stoichiometric coefficients
     &          preq( 6,2,nq),    ! parameters for reaction rate
     &          preq1(6,nq),
     &          creq(ns,  nq),    ! index number of mol
                                  ! concentration(overall)
     &          mreq(ns,  nq),    ! coefficients of 3rd body
     &      sub_vreq(ns,  nq),    ! subtraction of "vreq" of
                                  ! each species in each equation
     &    sigma_vreq(     nq),    ! change of "sub_vreq" in each equation
     &       vreq_a7( 7,2,nq),    ! summation of "vreq" * "a7"
                                  ! in each equation
     &     kind_chem(     nq),
     &     kind_prs (     nq),
     &     ig_iter  (     nq),
     &     stp_iter  (     nq),
     &          lreq( 3,2,nq),    ! parameters for low pressure
     &           l_sw(  2,nq),    ! flag of Lindemann form
     &           troe(4,2,nq),    ! Troe form
     &     UNIT_SI  (     nq),
     &     const_R  (     nq),
     &    preq_LH_ER(5,ns,nq),
     &    nosp_stick(     nq),
     &     nosp_Bohm(     nq),
     &      coefBohm(     nq),
     &      Bohm_sub(ns,2,nq),
!
     &      BV_ref_I(nq),BV_a(2,nq),BV_para(ns,2,nq),
     &      BV_ref_molfrc(ns,2,nq),
     &      LDCOM(nq,2),
     &                stat=ierr1 )
   
!
      if(ierr1/=0) then
        write(ifle,*) '### error : allocation failed'
        ierror = 1
      endif
!---------------------------
! --- initialization
!---------------------------
      ireq = 0
      M_3rd=-1
      P_3rd=-1
      vreq = 0.d0
      preq = 0.d0
      preq1= 0.d0
      creq = 0.d0
      mreq = -1.d0
      sub_vreq = 0.d0
      sigma_vreq = 0.d0
      vreq_a7 = 0.d0
      kind_chem=igas
      kind_prs=0
      ig_iter=0
      stp_iter=10000000
      lreq(1:3,1:2,1:nq)=undef !0.d0
      l_sw(1:2,1:nq)=1      ! unimolecular
      troe(1:4,1:2,1:nq)=0.d0
      UNIT_SI=0
      preq_LH_ER=0.d0
      coefBohm(:)=1.d0
      Bohm_sub(:,:,:)=0
      nosp_Bohm(:)=1
      nosp_stick(:)=1
      BV_ref_I(:)=0.d0
      BV_a(:,:)=0.d0
      BV_para(:,:,:)=0.d0
!
      LDCOM(:,:)=0
!---------------------------
!-< 2. Input namelist >-
!---------------------------
!--< 2.1 read namelist >--
!---------------------------
!      model = ' '
!      rewind ifli
!      read(ifli,chemreac,iostat=ios)            ! read of namelist
!      call nml_listno(11,mdllst,model,imodl)
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
        stop_iter=10000000
        model= ' '
        kind=' '
        prs_kind=' '
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
        M     = 1.d0           ! coefficients of 3rd body
        pleft=undef
        pright=undef
        unit_A_e=0
! pressure dependance
        af_low     = undef       ! for Lindemann etc.
        alphaf_low = undef
        ef_low     = undef
        af_high    = undef       ! for Lindemann etc.
        alphaf_high= undef
        ef_high    = undef

        troe_a(1:2) = undef	 ! Troe form
        troe_t1(1:2) = undef
        troe_t2(1:2) = undef
        troe_t3(1:2) = undef        
!/ eddy break up /
        cr1   = undef
        cr2   = undef
        r     = undef
!/ LH_ER /
        a_LH_ER=undef
        beta_LH_ER=undef
        H_LH_ER=undef
        L_LH_ER=undef
        n_LH_ER=undef
        m_LH_ER=undef
        k_LH_ER=undef
!/Y_BOHM & E_BOHM/
        coeff_BOHM=undef
        kd=0
!/Butler_Volmer/
        I_ref_BV=undef
        alp_a_BV=undef
        alp_c_BV=undef
        r_BV_left=undef
        r_BV_right=undef
        ref_mol_frc_left=undef
        ref_mol_frc_right=undef
!/Zeldovich/
        coal_CHO(q)=0
!----------------------------------------------------
        read(ifli,chemreac,iostat=ios)                  ! read of namelist
        call nml_errmsg0(ifle,ios,'chemreac',ierr1)     ! if ios>0,ierr1=1
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
        if(ierr1/=0) call FFRABORT(1,'&chemreac/kind')
        kind_chem(q)=kd
        chem_bc(q)=kd
	if(kind_chem(q)==isurface) then
	  lsuf=.true.
	  ical_suf=1
!        elseif(kind_chem(q)==iparticle) then
!          ical_P_R=1
	endif
!-----------------------------------------------
! --- identify pressure-depend reaction kind 
!-----------------------------------------------
        kd=0
        call nml_listno(2,prslst,prs_kind,kd)
        if(ierr1/=0) call FFRABORT(1,'&chemreac/prs_kind')
        kind_prs(q)=kd

!-------------------------
! --- identify model no.
!-------------------------
        call nml_listno(15,mdllst,model,imodl)
        call nml_chkchr0(ifle,'model',model,imodl,ierr1)
        if( ierr1.ne.0 ) then
          write(ifle,'(1X,a)') 'ERR: &chemreac/model is unknown param.'
          write(ifle,*) 
     & " ### MSG: model='overall','elementary','eddybreakup','fire',",
     & " or 'user'"
          goto 9999
        endif
        ireq(q)=imodl
!----------------------------------
! --- check gas/surface model
!----------------------------------
        if(kind_chem(q)==igas) then 
          if(.not.(imodl==ovrall.or.
     &             imodl==oval_ebrkup.or.
     &             imodl==elemnt.or.
     &             imodl==ebrkup.or.
     &             imodl==ovalfr.or.
     &             imodl==elemarbi.or.
     &             imodl==Zeldovich.or.
     &             imodl==SSFRRM.or.
     &             imodl==userreac)
     &       ) then
            write(ifle,'(1X,a,1x,a,1x,a)') 'ERR: &chemreac/model :',
     &      trim(mdllst(imodl)),'is NOT for GAS reaction'
            write(ifle,'(1X,a,I4)') 'MSG: Reaction no=',q
            call FFRABORT(1,'ERR: GAS reaction')
          endif
        elseif(kind_chem(q)==isurface) then
          if(.not.(imodl==stick_sf.or.
     &             imodl==LH_ER.or.
     &             imodl==Bohm.or.
     &             imodl==BohmE.or.
     &             imodl==BohmY.or.
     &             imodl==elemarbi.or.
     &             imodl==elemnt.or.
     &             imodl==Butler_Volmer.or.
     &             imodl==userreac)
     &      ) then
            write(ifle,'(1X,a,1x,a,1x,a)') 'ERR: &chemreac/model :',
     &      trim(mdllst(imodl)),'is NOT for Surface reaction'
            write(ifle,'(1X,a,I4)') 'MSG: Reaction no=',q
            call FFRABORT(1,'ERR: Surface reaction')
          endif
!        elseif(kind_chem(q)==iparticle) then
!          if(.not.(imodl==P_CHAR_COM.or.
!     &             imodl==P_CHAR_NO.or.
!     &             imodl==P_HCN_O2_NO.or.
!     &             imodl==P_NO_N2O2.or.
!     &             imodl==P_prompt_NO.or.
!     &             imodl==P_Zeld_NO.or.
!     &             imodl==P_userreac)
!     &      ) then
!            write(ifle,'(1X,a,1x,a,1x,a)') 'ERR: &chemreac/model :',
!     &      trim(mdllst(imodl)),'is NOT for Particle reaction'
!            write(ifle,'(1X,a,I4)') 'MSG: Reaction no=',q
!            call FFRABORT(1,'ERR: Particle reaction')
!          endif
        else 
           write(ifle,'(1X,a,I4)') 'MSG: Reaction no=',q
          call FFRABORT(1,'ERR: NOT surport reaction model')
        endif 
!-------------------------
! --- ignition iter
!-------------------------
      if(ignition_iter.le.0) then
        ig_iter(q)=0
      else
        ig_iter(q)=ignition_iter
      endif
      if(stop_iter.le.0) then
        stp_iter(q)=0
      else
        stp_iter(q)=stop_iter
      endif
      if(ig_iter(q)<igT_iter) then
        write(ifle,*) ' ### WRN: Reaction no ',q
        write(ifle,*)
     & ' ### WRN: &chemreac/ignition_iter < &chemcntl/igT_iter'
        write(ifle,*)
     & ' ### WRN: Reaction should start on iter= ',igT_iter
      endif
!
! --- unit :UNIT_SI
!
      if(unit_A_e==0.or.unit_A_e==1.or.unit_A_e==2) then
        UNIT_SI(q)=unit_A_e
      else
        write(ifle,'(a)') 'ERR: unit_A_e MUST be 0,1 or 2'
        call FFRABORT(1,'err: unit_A_e')
      endif
!--------------------------------------------------
! --- check duplication between left and right side
!--------------------------------------------------
!        if(imodl==ovrall.or.imodl==ebrkup.or.imodl==ovalfr)
      if(imodl==ebrkup)
     &  then
          call chk_lr(nleft,nright)
      endif
!        call chk_lr(nleft,nright)
! --- check duplication of left side
        call chk_dup('nleft',nleft)
! --- check duplication of right side
        call chk_dup('nright',nright)
!---------------------------------------------------
! --- Check Ion-energy- Bohm & Ion-enhanced Bohm
!---------------------------------------------------
        if(imodl==BohmE.or.imodl==BohmY) then
          if(coeff_BOHM==undef) then
          if(imodl==BohmE) write(ifle,'(a)') 
     &    'ERR: corff. for Ion-energy-dependent BOHM model in eq.: ',q
          if(imodl==BohmE)  write(ifle,'(a)') 
     & 'ERR: corff. for Ion-enhanced Reac. BOHM Yield Model in eq.: ',q
          call FFRABORT(1,'ERR: Not input [coeff_BOHM] in [&chemreac]')
          endif
        else
          coefBohm(q)=coeff_BOHM
        endif
!--------------------------------------------
! --- pick up number of species in each side
!--------------------------------------------
        call pick_up_species
        select case(imodl)
          case(ovrall)              !1)
            call read_overall
          case(elemnt)              !2)
            call read_elementary
          case(ebrkup)              !3)
            call read_eddy_break_up
          case(ovalfr)              !4)
            call read_overall_fire
          case(oval_ebrkup)
            call read_eddybreakup_overall
          case(SSFRRM)
            call read_overall
          case(userreac)            !5)
            call read_user
          case(stick_sf)            !6)
            call read_stick
          case(LH_ER)               !7)
            call read_LH_ER
          case(Bohm)                !8)
            call read_stick
          case(elemarbi)
            call read_elementary
          case(BohmE)
            call read_user
          case(BohmY)
            call read_user
          case(Butler_Volmer)
             call read_Butler_Volmer
          case(Zeldovich)
             call read_Zeldovich
        end select
 200  continue
      if(ierror/=0) call FFRABORT(1,'ERR: reading error in &chemreac')
!------------------------------------------------
! --- 
!------------------------------------------------
      do q=1,nq
      do j=1,ncompall
      dum1=1.d10
      dum3=-1.d-10
      dum2=sub_vreq(j,q)   !vreq(j,1,q)
      if(dum2<-1.d-20) then
        if(dum2<dum1) then
          dum1=dum2
          LDCOM(q,1)=j
        endif
      endif
      if(dum2>1.d-20) then
        if(dum2>dum3) then
          dum3=dum2
          LDCOM(q,2)=j
        endif
      endif
      enddo
      enddo
!--------------------------------
! --- display of reading results
!--------------------------------
!
      do q=1,nq
      imodl=ireq(q)
      if(imodl==ovrall) then
        call disp_overall(q)
      elseif(imodl==ovalfr) then
        call disp_overall(q)
      elseif(imodl==elemnt) then
        call disp_elementary(q)
      else
        call disp_elementary(q)
      endif
      enddo

!
!----------------------
! --- final error check
!----------------------
!
      return
 9999 continue
      ierror=1
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
        real*8,intent(in) :: value
!
        c_no=chk_nam(ns,spcnam,name)
! --- if "c_no" is identified(c_no>0), ierr1=0
        call nml_chkchr0(ifle,str1,name,c_no,ierr1)
        if( ierr1/=0 ) then
          write(ifle,'(2X,2a)') 
     &     'ERR: Unknown species name: ',trim(name)
          write(ifle,'(2X,a,I4)') 
     &                "ERR: Unknown species name in equation",q
          call FFRABORT(1,'ERR:') 
        endif
        if(value==undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 
     & 'ERR: lack of data : ',str2,'(i)',' in reaction eq.=',q
          write(ifle,*) 'i =',c_no
          if(c_no==ns+1) write(ifle,*) 'ERR: Lack of 3rd M coeffient'
          call FFRABORT(1,'ERR:')
!        elseif( value<1 ) then
!          write(ifle,*) '### error : data error'
!          write(ifle,*) str2,'(i) = ',value
!          write(ifle,*) 'it must be > 0'
!          write(ifle,*) 'i =',c_no
!          ierror = 1
        endif
      end function c_no
!
!=================================================
       integer function c_nop(name,str1,value,str2)
!=================================================
        character(lenspc),intent(in) :: name
        character,intent(in) :: str1,str2
        real*8,intent(in)    :: value
!
        c_nop=chk_nam(ns,spcnam,name)
        call nml_chkchr0(ifle,str1,name,c_nop,ierr1)
        if( ierr1/=0 ) then
          write(ifle,*) "did not identify species name in equation",q
          ierror = 1
        endif
      end function c_nop
!
!===========================================================
      subroutine pick_up_species  ! pick up species number
!===========================================================
      implicit none
      integer      :: s,nostk,BYflag(2),no_elec,no_plsm,ii
      character(1) :: plasma,sharp,percnt
      real*8 :: wsum
      logical :: nml_comp_eq
      integer,save :: coal=0
!
      nol = 0            ! number of species in left side
      nor = 0            ! number of species in right side
!
      nostk=0
      BYflag=0
      no_elec=0
      no_plsm=0
!---------------------------------------
! --- check '#' in Sub-Reaction species
!---------------------------------------
      do 100 s=1,ns
      sharp=''
      if(nleft(s)/=' ') then
        sharp=nleft(s)(1:1)
        if(sharp=='#') then
          if(imodl/=BohmY) then
            write(ifle,'(a,I4,a)') 'ERR: Reaction No: ',q,
     &                ' is not Ion-enhanced Reac. Yield Bohm Model'
            write(ifle,'(a)') 
     &      'MSG: [#] can NOT be placed in front of species name'
            call FFRABORT
     &      (1,'ERR: Sub-Reaction can only be used by Y_BOHM model')
          endif
          BYflag(1)=BYflag(1)+1
        endif
      endif
      sharp=''
      if(nright(s)/=' ') then
        sharp=nright(s)(1:1)
        if(sharp=='#') then
          if(imodl/=BohmY) then
            write(ifle,'(a,I4,a)') 'ERR: Reaction No: ',q,
     &                ' is not Ion-enhanced Reac. Yield Bohm Model'
            write(ifle,'(a)') 
     &      'MSG: [#] can NOT be placed in front of species name'
            call FFRABORT
     &      (1,'ERR: Sub-Reaction can only be used by Y_BOHM model')
          endif
          BYflag(2)=BYflag(2)+1
        endif
      endif
 100  enddo
!------------------------------------------------
! --- check stick,Bohm,Y_Bohm,E_BOHM model 
! --- search species number
!------------------------------------------------
      do s=1,ns
        if(nleft(s)/=' ') then         ! left side
          sharp=nleft(s)(1:1)
          if(sharp=='#') then
            nleft(s)(1:1)=' '
            nleft(s)=adjustL(nleft(s))
          endif
          j=c_no(nleft(s),'nleft',cleft(s),'cleft')
          if(j<=ns) then
!
            if(imodl==stick_sf.and.j<=ns_gs) then
              nosp_stick(q)=j
              nostk=nostk+1
              if(nostk>=2) then
                write(ifle,'(2X,3a,I4)')
     &          'MSG: [stick] model only allow ',
     &          'exactly one gas-phase species reactant ',
     &          'in reaction equation: ',q
                call FFRABORT(1,'ERR: stick model using error')
              endif
            endif
!
            if(imodl==Bohm.or.imodl==BohmE.or.imodl==BohmY) then
              do i=1,lenspc
              plasma=nleft(s)(i:i)
              if(plasma=='+'.or.plasma=='-') then
                nosp_Bohm(q)=j
                no_plsm=no_plsm+1
                if(no_plsm>=2) then
                  write(ifle,'(2X,3a,I4)')
     &            'MSG: [Bohm] model only allow ',
     &            'exactly one gas-phase ion reactant',
     &             ' in reaction equation: '
     &             ,q
                  call FFRABORT(1,'ERR: Bohm model using error')
                endif
              endif
              if(plasma=='-') then
                write(ifle,'(a)') 
     &       'ERR: Surface reaction modelor nigative ion NOT finished'
                call FFRABORT(1,'MSG: Call your supportor')
              endif
              enddo
              if(trim(nleft(s))=='e'.or.trim(nleft(s))=='E') then
                no_elec=no_elec+1
              endif
              if(sharp=='#') then
                Bohm_sub(j,1,q)=1
              endif
            endif
!----------------------------------------------
! --- stoichiometric coefficient of left side
!----------------------------------------------
            vreq(j,1,q)=(cleft(s))
            cofleft(ns*(q-1)+j)=(cleft(s))
            nol=nol+1
          elseif(j==ns+1) then
            M_3rd(q)=1
            do is=1,ns_gs
            if(M(is)>=0.d0) then
              mreq(is,q)=M(is)       ! 3rd body coefficient
            endif
            enddo
            do ii=1,ns_gs
            if(M(ii)<0) then
              write(ifle,'(2X,a,I4)') 
     &        'MSG: at reaction reaction no=',q
               write(ifle,'(2X,2a)')
     &          'MSG: 3rd-body effect from ',trim(spcnam(ii)),
     &          ' Should be 1.0 in default value'
              call FFRABORT(1,'ERR: [M] value is NOT defined')
            endif
            enddo
          endif
          ispc(nol)=j
        endif
        if(nright(s)/=' ') then        ! right side
          sharp=nright(s)(1:1)
          if(sharp=='#') then
            nright(s)(1:1)=' '
            nright(s)=adjustL(nright(s))
          endif
          j=c_no(nright(s),'nright',cright(s),'cright')
          if(j<=ns) then
            if(sharp=='#') then
              Bohm_sub(j,2,q)=1
            endif
!----------------------------------------------
! --- stoichiometric coefficient of right side
!----------------------------------------------
            vreq(j,2,q)=(cright(s))
            nor=nor+1
          endif
        endif
      enddo
!----------------------------------------------
!
!----------------------------------------------
      do s=1,ns
        sub_vreq(s,q)=vreq(s,2,q)-vreq(s,1,q)
        net_vreq(ns*(q-1)+s)=sub_vreq(s,q)
      enddo
      do s=1,ns_gs
        sigma_vreq(q)=sigma_vreq(q)+sub_vreq(s,q)
      enddo
!
      if(nostk==0.and.imodl==stick_sf) then
        write(ifle,'(2X,3a,I4)')
     &       'MSG: [stick] model only allow ',
     &  'exactly one and at-least one gas-phase species reactant ',
     &       'in reaction equation: ',q
        call FFRABORT(1,'ERR: Can NOT used stick model')
      endif
!
      if(imodl==BohmY.and.(BYflag(1)==0.and.BYflag(2)==0)) then
        write(ifll,'(a,I4)') 
     &    'WRN: No Sub-Reaction species is specified in Reac.No= ',q
      elseif(imodl==BohmY.and.(BYflag(1)/=BYflag(2))) then
        write(ifll,'(2a,I4)') 
     &'WRN: Please Check Sub-Reaction if/not in element conservation',
     &  ' in  Reac.No= ',q
      endif
!
      if(no_plsm/=no_elec) then
        write(ifle,'(a,I4,4X,a,I4)') 'Plasma number= ',no_plsm,
     &  'Electron number',no_elec
        write(ifle,'(2a,I4)') 
     &  'ERR: Electron number NOT equal to Plasma number ',
     &    'in reaction eq.=',q
        call FFRABORT(1,'ERR: Plasma reaction error')
      endif
!
      wsum=0.d0
      do s=1,ns!_gs
        if(trim(spcnam(s))/='E'.and.trim(spcnam(s))/='e'.and.
     &     trim(spcnam(s))/='E-'.and.trim(spcnam(s))/='e-'
     &  ) then
          wsum=wsum+sub_vreq(s,q)*wm(s)
        endif
      enddo
!
      do s=1,ns
      if(nleft(s)=='Ca_Hb_Oc') then
        coal_CHO(q)=1
      endif
      enddo
!-------------------------
! --- coal conbustion
!-------------------------
      if(coal_CHO(q)==1) then
        if(.not.(imodl==oval_ebrkup.or.imodl==Zeldovich
     &     .or.imodl==SSFRRM)) then
          write(ifle,'(2X,2a)') 
     & " ERR: Coal combustion only',
     & ' supports ovrall-eddyBreakUp combined model"
          write(ifle,'(2X,a)') 'MSG: set [model=eddy_over] in '
          write(ifle,'(2X,a,I4)') 'MSG: reacton no= ',q
          call FFRABORT(1,'ERR: in &chemreac')
        endif
        coal=1
      endif 
!
      if(imodl==Zeldovich.and.coal/=1) then
        if(coal_CHO(q)/=1) then
          write(ifle,'(2X,a)') 
     &    "ERR: Zeldovich-NO only supports coal combustion"
          call FFRABORT(1,'ERR: in &chemreac')
        endif
      elseif(imodl==Zeldovich.and.q==1) then
        write(ifle,'(2X,a)') 
     &    "ERR: Zeldovich-NO model must be last &chemreac"
        call FFRABORT(1,'ERR: in &chemreac')
      endif
!
      if(abs(wsum)>1.d-15.and.coal_CHO(q)==0) then
        write(ifle,'(1x,a,I4)') 
     &  'WRN: Mass NOT equal to on both site ,Reaction eq= ',q
        write(ifle,'(1x,a,E14.6)') 'Difference mass is ',wsum
        call FFRABORT
     & (1,'MSG: check stoichiometric coefficient or molecular weight')
      endif
!
      if((imodl==ovrall.or.imodl==ovalfr).and.coal_CHO(q)==0) then
        if(nml_comp_eq("nu_fuel",nu_fuel,undef,1,ifle,ierror)) then
          write(ifle,'(1x,2a,I4)') 
     &  'ERR: Lack of Stoichiometric coefficient of ',
     &  'fuel for overall model ,reac. eq=',q
          call 
     &    FFRABORT(1,'MSG: Specifing  "nu_fuel" in [&chenreac]')
        else
          if(nu_fuel<1.d-10) then
            call 
     &      FFRABORT(1,'ERR: [nu_fuel] must be > 0.0 ')
          else
            do s=1,ns_gs
            vreq(s,2,q)=vreq(s,2,q)/nu_fuel
            vreq(s,1,q)=vreq(s,1,q)/nu_fuel
            enddo
          endif
        endif
      endif
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
      integer :: s
!
      do i=1,nol
        if( nml_comp_eq("creac(i)",creac(i),undef,i,ifle,ierror) ) then
          if( .not.nml_comp_ge("creac(i)",creac(i),0,i,ifle,ierror) )
     &      exit
        endif
      enddo
      k=1
      do s=1,ns
        j=c_no(nleft(s),"nleft",cleft(s),"cleft")
!
! --- if j is identified(j>0), ierr1=0
!
        call nml_chkchr0(ifle,'nleft',nleft(s),j,ierr1)
        if( ierr1/=0 ) then
          write(ifle,*) "'nleft' is not identified in equation",q
          ierror=1
        endif
!
! --- index number of mol concentration(overall)
!
        if( j<=ns ) then
          creq(j,q) = creac(k)
        endif
        if(k<nol) then
          k=k+1
        else
          exit
        endif
      enddo
      if(.not.nml_comp_eq("a",a,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("a",a,0.0d0,1,ifle,ierror) ) then
          preq(1,1,q)=a
        endif
      endif
      
      if(.not.nml_comp_eq("alpha",alpha,undef,1,ifle,ierror) )
     &  preq(2,1,q)=alpha
      if(.not.nml_comp_eq("beta",beta,undef,1,ifle,ierror) )
     &  preq(3,1,q)=beta
      if(.not.nml_comp_eq("e",e,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("e",e,0.0d0,1,ifle,ierror) ) then
          preq(4,1,q)=e
        endif
      endif
      end subroutine read_overall
!
!=======================================================
      subroutine read_overall_fire           ! imodl==1
!=======================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      
      character(lenspc) :: soot_nam,co_nam,o2_nam,co2_nam
      integer :: s
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
!
      do s=1,ns
        j=c_no(nleft(s),"nleft",cleft(s),"cleft")
        call nml_chkchr0(ifle,'nleft',nleft(s),j,ierr1)
        if(ierr1/=0) then
          write(ifle,*) "'nleft' is not identified in equation",q
          ierror=1
        endif
        if(j<=ns) then
          creq(j,q)=creac(k)
        endif
        if(k<nol) then
          k=k+1
        else
          exit
        endif
      enddo
!
      if(.not.nml_comp_eq("a",a,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("a",a,0.0d0,1,ifle,ierror) ) then
          preq(1,1,q)=a
!          dum=0.d0
!          do s=1,ncmpx          !ns
!          if(abs(vreq(s,1,q))>SML) then
!            dum=dum+creq(s,q)  ! vreq(s,1,q)
!          endif
!          enddo
!          dum1=0
!          do s=1,ncomp_sufx     !ns
!          if(abs(vreq(s,1,q))>SML) then
!            dum1=dum1+creq(s,q) !vreq(s,1,q)
!          endif
!          enddo
!          if(unit_A_e==0) then
!            preq(1,1,q)=a
!          elseif(unit_A_e==1) then
!            if(kind_chem(q)==igas) then
!              preq(1,1,q)=a*(10.d0**(-6*(dum-1.d0)))
!            else
!              preq(1,1,q)=a*(10.d0**(-6*dum-4*(dum1-1.d0)))
!            endif
!          else
!            if(my_rank==0) then
!              write(ifle,'(a)') 'ERR: unit_A_e MUST be 0 or 1'
!              write(ifle,'(a)') 
!     &        'MSG: unit_A_e=0: A=m^3*(n-1)/[mole^(n-1)*s*K^(-alpha)]'
!              write(ifle,'(12X,a)') '     E=J/mole'
!              write(ifle,'(a)')
!     &        'MSG: unit_A_e=1: A=cm^3*(n-1)/[mole^(n-1)*s*K^(-alpha)]'
!              write(ifle,'(12x,a)') '     E=cal/mole'
!              write(ifle,'(a,I4)') 'MSG: n= ',dum
!              write(ifle,'(2a)')
!     &       'MSG: Here, [n] is reactant species number ',
!     &       'considering stoichiometric coefficient '
!            endif
!            call FFRABORT(1,'ERR:')
!          endif
        endif
      endif

!
      if(.not.nml_comp_eq("alpha",alpha,undef,1,ifle,ierror))
     &    preq(2,1,q)=alpha
!
      if(.not.nml_comp_eq("beta",beta,undef,1,ifle,ierror))
     &    preq(3,1,q)=beta
!
      if(.not.nml_comp_eq("e",e,undef,1,ifle,ierror)) then
        if(nml_comp_gt("e",e,0.0d0,1,ifle,ierror)) then
!          if(unit_A_e==0) then
            preq(4,1,q)=e
!          elseif(unit_A_e==1) then
!            preq(4,1,q)=e*4.18605d0  ! 1[cal]=4.18605[J]
!          else
!            call FFRABORT(1,'ERR:unit_A_e MUST be 0 or 1')
!          endif
        endif
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
        vreq(j,2,q)=nu_SOOT
! --- CO Nu
        j= chk_nam(ns,spcnam,co_nam)
        vreq(j,2,q)=nu_CO
! --- O2 Nu
        j=chk_nam(ns,spcnam,o2_nam)
        vreq(j,1,q)=vreq(j,1,q)-nu_CO-nu_SOOT
! --- CO2 Nu
        j=chk_nam(ns,spcnam,co2_nam)
        vreq(j,2,q)=vreq(j,2,q)-nu_CO-nu_SOOT
        endif
!
        if(my_rank==0)
     &  write(ifll,'(a)')' ==== stoichiometric coefficients ===='
        do s=1,ns
        if(my_rank==0)write(ifll,'(a7,i3,a9,i3,10x,a8)')
     &           'Reac.= ',q,'   spec.no.=',s,spcnam(s)
        if(my_rank==0)
     &  write(ifll,'(5x,e13.5,5x,e13.5)') vreq(s,1,q),vreq(s,2,q)
        enddo
      end subroutine read_overall_fire
!
!===================================================
      subroutine read_elementary        ! imodl==2
!===================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      real*8 :: dum0
      integer :: s
!
      if(.not.nml_comp_eq("af",af,undef,q,ifle,ierror)) then
        if(nml_comp_gt("af",af,0.0d0,q,ifle,ierror)) then
          preq(1,1,q)=af
        endif
      endif
      if(.not.nml_comp_eq("alphaf",alphaf,undef,q,ifle,ierror)) then
        preq(2,1,q)=alphaf
      endif
      if(.not.nml_comp_eq("ef",ef,undef,q,ifle,ierror)) then
        preq(3,1,q)=ef
      endif
!
      if(kind_prs(q)==falloff) then
        if(.not.nml_comp_eq("af_low",af_low,undef,q,ifle,ierror)) then
          if( nml_comp_gt("af_low",af_low,0.0d0,q,ifle,ierror)) then
            lreq(1,1,q)=af_low
            P_3rd(q,1)=1
          endif
        endif
        if(.not.nml_comp_eq("alphaf_low",alphaf_low,
     &                                          undef,q,ifle,ierror))
     &      lreq(2,1,q)=alphaf_low
        if( .not.nml_comp_eq("ef_low",ef_low,undef,q,ifle,ierror))
     &      lreq(3,1,q)=ef_low
        if(ierror/=0) 
     &    call FFRABORT(1,'ERR: fall-off parameter error')
      endif
      if(kind_prs(q)==activated) then
        if(.not.nml_comp_eq("af_high",af_high,undef,q,ifle,ierror)) then
          if( nml_comp_gt("af_high",af_high,0.0d0,q,ifle,ierror)) then
            lreq(1,1,q)=preq(1,1,q)
            preq(1,1,q)=af_high
            P_3rd(q,1)=1
          endif
        endif
        if(.not.nml_comp_eq("alphaf_high",alphaf_high,
     &                            undef,q,ifle,ierror)) then
          lreq(2,1,q)=preq(2,1,q)
          preq(2,1,q)=alphaf_high
        endif
        if(.not.nml_comp_eq("ef_high",ef_high,undef,q,ifle,ierror)) then
          lreq(3,1,q)=preq(3,1,q)
          preq(3,1,q)=ef_high
        endif
        if(ierror/=0) 
     &    call FFRABORT(1,'ERR: activated parameter error')
      endif
      if(kind_prs(q)==activated.or.kind_prs(q)==falloff) then
        if(troe_a(1)/=undef) then		! Troe form
          troe(1,1,q) = troe_a(1)
          if(.not.nml_comp_eq("troe_t3(1)",
     &            troe_t3(1),undef,q,ifle,ierror))
     &      troe(2,1,q) = troe_t3(1)
          if(.not.nml_comp_eq("troe_t1(1)",
     &            troe_t1(1),undef,q,ifle,ierror))
     &      troe(3,1,q) = troe_t1(1)
          troe(4,1,q) = troe_t2(1)
        endif
      endif
! --------------------------------------------------------------------
! --- ab=: a_backward
! --- =undef:calc. backward reaction with forward one (Both F and B)
! --- =0    :NOT consider backward reaction           (Only F)
! --- other :calc. backward reaction without forward one (Only B)
      preq(1,2,q)=ab
      if(preq(1,2,q)/=undef.and.abs(preq(1,2,q))>SML) then
        preq(2,2,q)=alphab
        preq(3,2,q)=eb
      endif
!
! --------------------------------------------------------------------
      
! --------------------------------------------------------------------
      do k=1,2
        do j=1,7
! --- summation of "vreq" * "a7" in each equation
          do s=1,ns
          vreq_a7(j,k,q)=vreq_a7(j,k,q)+sub_vreq(s,q)*a7(j,k,s)
          enddo
        enddo
      enddo
!----------------------------
! --- arbitrary power index
!----------------------------
      if(imodl==elemarbi) then
        do s=1,ns
        if(nleft(s)/=' ') then         ! left side
          j=c_nop(nleft(s),'nleft',pleft(s),'pleft')
          if(j<=ns) then
! --- arbitrary power index of left side
            if(pleft(s)/=undef) vreq(j,1,q)=pleft(s)
          endif
        endif
!
        if(nright(s)/=' ') then        ! right side
          j=c_nop(nright(s),'nright',pright(s),'pright')
          if(j<=ns) then
! --- arbitrary power index of right side
            if(pright(s)/=undef) vreq(j,2,q)=pright(s)
          endif
        endif
        enddo
      endif
      end subroutine read_elementary
!
!=================================================
      subroutine read_Zeldovich     ! imodl==13
!=================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      integer :: s
!      
      do s=1,ns
        if(nleft(s)=='O2') then
          IZ_O2=chk_nam(ns,spcnam,nleft(s))
        endif
        if(nleft(s)=='N2') then
          IZ_N2=chk_nam(ns,spcnam,nleft(s))
        endif
        if(nright(s)=='NO') then
          IZ_NO=chk_nam(ns,spcnam,nright(s))
        endif
      enddo
        
      if(IZ_O2==0) then
        call FFRABORT(1,'MSG: IZ_O2=0 in Zeldovich')
      endif
      if(IZ_N2==0) then
        call FFRABORT(1,'MSG: IZ_N2=0 in Zeldovich')
      endif
      if(IZ_NO==0) then
        call FFRABORT(1,'MSG: IZ_NO=0 in Zeldovich')
      endif
      end subroutine read_Zeldovich

!===================================================
      subroutine read_stick        ! imodl==2
!===================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
!      
      if(.not.nml_comp_eq("a",a,undef,q,ifle,ierror) ) then
        preq(1,1,q)=a
      else
        write(ifle,'(2a,I4)') 
     &  'ERR: Prefactor [a] of Arrhenius NOT defined in &chemcntl ',
     &  'at eq.=',q
        call FFRABORT(1,'ERR: Maybe [af] be defined, ')
      endif
      if(.not.nml_comp_eq("alpha",alpha,undef,q,ifle,ierror)) then
        preq(2,1,q)=alpha
      else
        write(ifle,'(2a,I4)') 
     &  'ERR: Exponent [alpha: T^alpha] of Arrhenius NOT defined ',
     &  'in &chemcntl at eq.=',q
        call FFRABORT(1,'ERR: Maybe [alphaf] be defined, ')
      endif
      if(.not.nml_comp_eq("e",e,undef,q,ifle,ierror)) then
!        if(unit_A_e==0) then
          preq(3,1,q)=e
!        elseif(unit_A_e==1) then
!          preq(3,1,q)=e*4.18605d0  ! 1[cal]=4.18605[J]
!        else
!          call FFRABORT(1,'ERR:unit_A_e MUST be 0 or 1')
!        endif
!      else
!         write(ifle,'(2a,I4)') 
!     &   'ERR: Activation [E]  of Arrhenius NOT defined ',
!     &   'in &chemcntl at eq.=',q
!         call FFRABORT(1,'ERR: reading error in read_stick')
      endif
!      
      end subroutine read_stick
!
!=================================================
      subroutine read_eddy_break_up     ! imodl==3
!=================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      character(lenspc) :: soot_nam,co_nam,o2_nam,co2_nam
      soot_nam='SOOT'
      co_nam='CO'
      o2_nam='O2'
      co2_nam='CO2'
!      
      preq(1,1,q)=dble(ispc(1))
      preq(2,1,q)=dble(ispc(2))
      if(.not.nml_comp_eq("cr1",cr1,undef,q,ifle,ierror)) then
        if(nml_comp_gt("cr1",cr1,0.0d0,q,ifle,ierror))
     &    preq(3,1,q)=cr1
      endif
      dd2=0.d0
      if(cr2==undef) then
        cr2=0.d0
        dd2=10.d0
!      elseif(cr2<=0.d0) then
!	cr2=0.d0
!        dd2=10.d0
!        if(nml_comp_gt("cr2",cr2,0.0d0,q,ifle,ierror) ) return
      endif
      preq(4,1,q)=cr2
      if(.not.nml_comp_eq("r",r,undef,q,ifle,ierror)) then
        if(nml_comp_gt("r",r,0.0d0,q,ifle,ierror))
     &    preq(5,1,q)=r
      endif
      preq(6,1,q)=dd2
!-------------
! --- Fire
!-------------
      if(Y_SOOT<1.d-10) return
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
        vreq(j,2,q)=nu_SOOT
! --- CO Nu
        j= chk_nam(ns,spcnam,co_nam)
        vreq(j,2,q)=nu_CO
! --- O2 Nu
        j=chk_nam(ns,spcnam,o2_nam)
        vreq(j,1,q)=vreq(j,1,q)-nu_CO-nu_SOOT
! --- CO2 Nu
        j=chk_nam(ns,spcnam,co2_nam)
        vreq(j,2,q)=vreq(j,2,q)-nu_CO-nu_SOOT
        endif
!

      end subroutine read_eddy_break_up
!

!=================================================
      subroutine read_eddybreakup_overall
!read_eddy_break_up     ! imodl==3
!=================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      character(lenspc) :: soot_nam,co_nam,o2_nam,co2_nam
      soot_nam='SOOT'
      co_nam='CO'
      o2_nam='O2'
      co2_nam='CO2'
!      
      preq1(1,q)=dble(ispc(1))
      preq1(2,q)=dble(ispc(2))
      if(.not.nml_comp_eq("cr1",cr1,undef,q,ifle,ierror)) then
        if(nml_comp_gt("cr1",cr1,0.0d0,q,ifle,ierror))
     &    preq1(3,q)=cr1
      endif
      dd2=0.d0
!
      if(cr2==undef) then
        cr2=0.d0
        dd2=10.d0
!      elseif(cr2<=0.d0) then
!	cr2=0.d0
!        dd2=10.d0
!        if(nml_comp_gt("cr2",cr2,0.0d0,q,ifle,ierror) ) return
      endif
!
      preq1(4,q)=cr2
      if(.not.nml_comp_eq("r",r,undef,q,ifle,ierror)) then
        if(nml_comp_gt("r",r,0.0d0,q,ifle,ierror))
     &    preq1(5,q)=r
      endif
      preq1(6,q)=dd2
!
!-------------
! --- overall 
!-------------
!
      do i=1,nol
        if(nml_comp_eq("creac(i)",creac(i),undef,i,ifle,ierror)) then
          if(.not.nml_comp_ge("creac(i)",creac(i),0,i,ifle,ierror))
     &      exit
        endif
      enddo
      k=1
      do s=1,ns
        j=c_no(nleft(s),"nleft",cleft(s),"cleft")
!--------------------------------------
! --- if j is identified(j>0), ierr1=0
!--------------------------------------
        call nml_chkchr0(ifle,'nleft',nleft(s),j,ierr1)
        if( ierr1/=0 ) then
          write(ifle,*) "'nleft' is not identified in equation",q
          ierror=1
        endif
!-------------------------------------------------
! --- index number of mol concentration(overall)
!-------------------------------------------------
        if( j<=ns ) then
          creq(j,q) = creac(k)
        endif
        if(k<nol) then
          k=k+1
        else
          exit
        endif
      enddo
      if(.not.nml_comp_eq("a",a,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("a",a,0.0d0,1,ifle,ierror) ) then
          preq(1,1,q)=a
        endif
      endif
!
      if(.not.nml_comp_eq("alpha",alpha,undef,1,ifle,ierror))
     &  preq(2,1,q)=alpha
      if(.not.nml_comp_eq("beta",beta,undef,1,ifle,ierror))
     &  preq(3,1,q)=beta
      if(.not.nml_comp_eq("e",e,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("e",e,0.0d0,1,ifle,ierror) ) then
          preq(4,1,q)=e
        endif
      endif
      end subroutine read_eddybreakup_overall
!
!=================================================
      subroutine read_user           ! imodl==5
!=================================================
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
!
      if(.not.nml_comp_eq("a",a,undef,q,ifle,ierror) ) then
        preq(1,1,q)=a
      else
        write(ifle,'(2a,I4)') 
     &  'ERR: Prefactor [a] of Arrhenius NOT defined in &chemcntl ',
     &  'at eq.=',q
        call FFRABORT(1,'ERR: Maybe [af] be defined, ')
      endif
      if(.not.nml_comp_eq("alpha",alpha,undef,q,ifle,ierror)) then
         preq(2,1,q)=alpha
      else
         write(ifle,'(2a,I4)') 
     &  'ERR: Exponent [alpha: T^alpha] of Arrhenius NOT defined ',
     &  'in &chemcntl at eq.=',q
        call FFRABORT(1,'ERR: Maybe [alphaf] be defined, ')
      endif
      if(.not.nml_comp_eq("e",e,undef,q,ifle,ierror)) then
         preq(3,1,q)=e
      else
         write(ifle,'(2a,I4)') 
     &   'ERR: Activation [E]  of Arrhenius NOT defined ',
     &   'in &chemcntl at eq.=',q
         call FFRABORT(1,'ERR: reading error in user')
      endif
      end subroutine read_user
!
!===================================================
      subroutine read_LH_ER          ! imodl==7
!===================================================
!
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      integer :: s
!     
      do s=1,ns

      if(.not.nml_comp_eq("a_LH_ER",a_LH_ER,undef,q,ifle,ierror)) then
        preq_LH_ER(1,s,q)=a_LH_ER(s)
      endif

      if(.not.nml_comp_eq("beta_LH_ER",beta_LH_ER,undef,q,ifle,ierror))
     &  preq_LH_ER(2,s,q)=beta_LH_ER(s)

      if(.not.nml_comp_eq("H_LH_ER",H_LH_ER,undef,q,ifle,ierror))
     &  preq_LH_ER(3,s,q)=H_LH_ER(s)

      if(.not.nml_comp_eq("L_LH_ER",L_LH_ER,undef,q,ifle,ierror))
     &  preq_LH_ER(4,s,q)=L_LH_ER(s)

      preq_LH_ER(5,s,q)=1.d0
      if(.not.nml_comp_eq("n_LH_ER",n_LH_ER,undef,q,ifle,ierror))
     &  preq_LH_ER(5,s,q)=n_LH_ER(s)
      enddo


      preq(1,1,q)=2.d0
      if(.not.nml_comp_eq("m_LH_ER",m_LH_ER,undef,q,ifle,ierror))
     &  preq(1,1,q)=m_LH_ER

      
      if(.not.nml_comp_eq("k_LH_ER",k_LH_ER,undef,q,ifle,ierror))
     &  preq(2,1,q)=k_LH_ER
! --------------------------------------------------------------------
      end subroutine read_LH_ER
!=================================================
      subroutine read_Butler_Volmer
!===================================================
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
      integer :: s,j
      integer :: I_H2_l=0,I_O2_l=0,I_H2O_l,I_PROTON_l,I_H2OL_l
!
      if(.not.nml_comp_eq("I_ref_BV",I_ref_BV,undef,q,ifle,ierror) ) 
     &   then
        BV_ref_I(q)=I_ref_BV
      else
        write(ifle,'(2a,I4)') 
     &  'ERR: Ref current [I_ref_BV] of Butler_Volmer NOT defined',
     &  'at eq.=',q
        call FFRABORT(1,'ERR: [I_ref_BV] NOT be defined, ')
      endif
!
      if(.not.nml_comp_eq("alp_a_BV",alp_a_BV,undef,q,ifle,ierror) ) 
     &   then
        BV_a(1,q)=alp_a_BV
      else
        write(ifle,'(2a,I4)') 
     &  'ERR: Ref current [alp_a_BV] of Butler_Volmer NOT defined',
     &  'at eq.=',q
        call FFRABORT(1,'ERR: [alp_a_BV] NOT be defined, ')
      endif
!
      
      if(.not.nml_comp_eq("alp_c_BV",alp_c_BV,undef,q,ifle,ierror) ) 
     &   then
        BV_a(2,q)=alp_c_BV
      else
        write(ifle,'(2a,I4)') 
     &  'ERR: Ref current [alp_c_BV] of Butler_Volmer NOT defined',
     &  'at eq.=',q
        call FFRABORT(1,'ERR: [alp_c_BV] NOT be defined, ')
      endif
!
      do s=1,ns
      if(nleft(s)/=' ') then
        if(.not.nml_comp_eq("r_BV_left(s)",r_BV_left(s),
     &      undef,q,ifle,ierror) ) 
     &   then
          j=c_no(nleft(s),'r_BV_left(s)',r_BV_left(s),'r_BV_left')
          BV_para(j,1,q)=r_BV_left(s)
        else
          write(ifle,'(2a,I4)') 
     &    'ERR:  [r_BV_left] of Butler_Volmer NOT defined',
     &    'at eq.=',q
          call FFRABORT(1,'ERR: [r_BV_left] NOT be defined, ')
        endif
        if(.not.nml_comp_eq("ref_mol_frc_left(s)",ref_mol_frc_left(s),
     &      undef,q,ifle,ierror) ) 
     &   then
          j=c_no(nleft(s),
     &   'ref_mol_frc_left',ref_mol_frc_left(s),'ref_mol_frc_left')
          
          BV_ref_molfrc(j,1,q)=ref_mol_frc_left(s)
        else
          write(ifle,'(2a,I4)') 
     & 'ERR: [ref_mol_frc_left] of Butler_Volmer NOT defined',
     &    'at eq.=',q
          call FFRABORT(1,'ERR: [ref_mol_frc_left] NOT be defined, ')
        endif
        endif
!
! --- 
!
      if(nright(s)/=' ') then
        if(.not.nml_comp_eq("r_BV_right(s)",r_BV_right(s),
     &      undef,q,ifle,ierror) ) 
     &   then
          j=c_no(nright(s),'r_BV_right(s)',r_BV_right(s),'r_BV_right')
          BV_para(j,2,q)=r_BV_right(s)
        else
          write(ifle,'(2a,I4)') 
     &    'ERR: [r_BV_right] of Butler_Volmer NOT defined',
     &    'at eq.=',q
          call FFRABORT(1,'ERR: [r_BV_right] NOT be defined, ')
        endif
        if(.not.nml_comp_eq("ref_mol_frc_right(s)",ref_mol_frc_right(s),
     &      undef,q,ifle,ierror) ) 
     &   then
          j=c_no(nright(s),
     &   'ref_mol_frc_right',ref_mol_frc_right(s),'ref_mol_frc_right')
          BV_ref_molfrc(j,2,q)=ref_mol_frc_right(s)
        else
          write(ifle,'(2a,I4)') 
     & 'ERR:  [ref_mol_frc_right] of Butler_Volmer NOT defined',
     &    'at eq.=',q
          call FFRABORT(1,'ERR: [ref_mol_frc_right] NOT be defined, ')
        endif
      endif
      
      enddo
!
      iBV=1
!
      if(q==1) then
        do s=1,ns
        
        if(nleft(s)=='H2') then
          I_H2_l=chk_nam(ns,spcnam,nleft(s))
        endif
        enddo
        if(I_H2_l==0) then
          call FFRABORT(1,'MSG: Anode reaction set to reaction-1')
        endif
      endif
      if(q==2) then
        do s=1,ns
        if(trim(nleft(s))=='O2') then
          I_O2_l=chk_nam(ns,spcnam,nleft(s))
        endif
        if(trim(nright(s))=='H2O-G') then
          I_H2O_l=chk_nam(ns,spcnam,nright(s))
        endif
        if(trim(nleft(s))=='H+') then
          I_PROTON_l=chk_nam(ns,spcnam,nleft(s))
        endif
        enddo
        
        if(I_O2_l==0) then
          call FFRABORT(1,'ERR: Cathode reaction-2:O2')
        endif
        if(I_H2O_l==0) then
          call FFRABORT
     &   (1,'ERR: Cathode reaction-2:H2O-G')
        endif
        if(I_PROTON_l==0) then
          call FFRABORT
     &   (1,'ERR: Cathode reaction-2:H2O-G')
        endif
      endif
!
      end subroutine read_Butler_Volmer
!=================================================

!
!< #2.1  check data >
!
!=================================================
      subroutine chk_dup(str,value)
!=================================================
      character(*),intent(in) :: str,value(mcomp+1)
      integer :: s
!
      do s=1,ns-1
        if( value(s)==' ' ) cycle
        do j=s+1,ns
          if( value(j)==' ' ) cycle
          if( value(j)==value(s) ) then
            write(ifle,*) '### error : data error'
            write(ifle,'(2X,7a,I4)') 'duplicated data : ',str,
     &       '(i)=',str,'(j)=',trim(value(s)),'  eq= ',q
            write(ifle,*) 'i,j =',s,j
            ierror = 1
            call FFRABORT(1,'stop at chk_dup')
          endif
        enddo
      enddo
      end subroutine chk_dup
!
!---------------------------------------------------
! --- check duplication between left and right side
!===================================================
      subroutine chk_lr(left,right)
!===================================================]
      character(lenspc),intent(in) :: left(mcomp+1)
                                   ! name of species on left side
      character(lenspc),intent(in) :: right(mcomp+1)
                                   ! name of species on right side
      integer :: s
!
      do s=1,ns+1
        if(left(s)/=' '.and. left(s)/='M'
     &   .and.(left(s)/='e'.and.left(s)/='E')) then
          do j=1,ns+1
          if(right(j)/=' '.and.right(j)/='M'
     &    .and.(right(j)/='e'.and.right(j)/='E')) then
            
            if( left(s)==right(j) ) then
            write(ifle,'(a)') '### error : data error'
            write(ifle,'(3a,I4)') 
     &      'ERR: duplicated data : nleft(i)=nright(j)='
     &                    ,trim(left(s)),' in eq.=',q
            write(ifle,*) 'i,j =',s,j
            ierror = 1
            call FFRABORT(1,'stop at chk_lr')
            endif
          endif
          enddo
        endif
      enddo
      end subroutine chk_lr
!
!-------------------------------------------------
! --- check species name
!=================================================
      integer function chk_nam(ns,lst,wrd)
!=================================================
!  chk_nam : = 0; unknown name (did not declare in namelist 'species')
!          : > 0; identified name (nn: number of identified name)
!          : =ns+1; 3rd body
        integer,intent(in)  :: ns
        character(lenspc),intent(in)  :: lst(ns)
        character(lenspc),intent(in)  :: wrd
        integer :: s
        chk_nam = 0
        do s=1,ns
          if(wrd==lst(s)) then
            chk_nam = s
            exit
          elseif(wrd=='M') then
            chk_nam=ns+1
            exit
          endif
        enddo
      end function chk_nam
!
!-------------------------------------------------
! --- following is display part
!=================================================
      subroutine disp_overall(q)
!=================================================
      integer,intent(in) :: q
      call disp_head(q)
      call disp_eq1(1)          ! left side
      if(my_rank==0) write(ifll,'(a5,$)') arrow(preq(1,2,1))
      call disp_eq1(2)          ! right side
      call disp_coefficients(1,4)
      call disp_creq(q)
      if(my_rank==0) write(ifll,'(/,$)')
      end subroutine disp_overall
!
!=================================================
      subroutine disp_elementary(q)
!=================================================
      integer,intent(in) :: q
      call disp_head(q)
        call disp_equation(q)
        if(my_rank==0) write(ifll,'(a)') ' '
!      call disp_vreq_a7
      end subroutine disp_elementary
!
!=================================================
      subroutine disp_vreq_a7
!=================================================
      do q=1,nq
        if(my_rank==0) write(ifll,"(i3,$)") q
        do k=1,2
          if(my_rank==0)
     &    write(ifll,"(e12.3,$)") (vreq_a7(j,k,q),j=1,7)
          if( k==1 .and.my_rank==0) write(ifll,"(/,a,$)") "   "
        enddo
        if(my_rank==0) write(ifll,*)
      enddo
      end subroutine disp_vreq_a7
!
!=================================================
      subroutine disp_head(q)
!=================================================
      integer,intent(in) :: q
!
      if(my_rank==0) then
         write(ifll,"(2x,108('+'))")
         write(ifll,*)
      endif
      if(my_rank==0) 
     &    write(ifll,"(2X,2a)") "Chemical reaction model     : ",
     &                        trim(mdllst(ireq(q)))
      if(my_rank==0) 
     & write(ifll,"(2x,a,i3)")
     &  "Total number of chemical equation : ",nq
      end subroutine disp_head
!
!=================================================
      subroutine disp_creq(q)
!=================================================
      integer,intent(in) :: q
      integer :: s
!
      if(my_rank==0) write(ifll,"(/,7x,a,$)") "index number     : "
      do s=1,ns
        if( creq(s,q)/=0 ) then
          if(my_rank==0) write(ifll,"(2x,a6,$)") trim(spcnam(s))
          if(my_rank==0) write(ifll,"(a,f5.1,$)") " ->",creq(s,q)
        endif
      enddo
      end subroutine disp_creq
!-------------------------------------------------
! --- display 1 equation for overall reaction
!=================================================
      subroutine disp_eq1(isw)
!=================================================
      integer,intent(in) :: isw  ! switch for side of equation
      integer :: s

      k = 0
      do s=1,ns
        if( vreq(s,isw,1)/=0 ) then
          if( k>0 .and.my_rank==0) write(ifll,"(a3,$)") " + "
          if(my_rank==0) write(ifll,"(f5.1,$)") vreq(s,isw,1)
          if(my_rank==0) write(ifll,"(1x,a6,$)") trim(spcnam(s))
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
      character*80 :: sinp
! classification number of equation
      if(my_rank==0) then
        call str_eqs(1,q,sinp)
        write(ifll,"(1x,i3,2x,a,$)") q,trim(sinp) ! left side
      endif
      if(my_rank==0) then
        write(ifll,"(a5,$)") arrow(preq(1,2,q))
      endif
      if(my_rank==0) then
        call str_eqs(2,q,sinp)
       write(ifll,"(a,$)") trim(sinp)          ! right side
      endif
      call disp_coefficients(q,3)
!      write(ifll,'(/,11x,9f4.0,$)') sub_vreq(1:ns,q)
      end subroutine disp_equation
!
!-------------------------------------------------
! --- display some coefficients
!=================================================
      subroutine disp_coefficients(q,n)
!=================================================
      integer,intent(in) :: q,n
      if(my_rank==0) then 
        write(ifll,"(/,7x,a,$)") "forward  reaction : "
        write(ifll,"(e12.3,$)") (preq(i,1,q),i=1,n)
        if(UNIT_SI(q)==0) then
          write(ifll,"(2x,a,$)")
     &    'in SI (A=mole-m-K ; e=J/mole) unit'
        else
          write(ifll,"(2x,a,$)")
     &    'in (cm-sec-mole-kcal-Kelvins)'
        endif
        if(lreq(1,1,q)/=undef ) call disp_Lindemann(1,q,n)
        write(ifll,"(/,7x,a,$)") "backward reaction : "
        if(preq(1,2,q)==undef) then
          write(ifll,"(a,$)") "<-- using equilibrium constants"
        elseif(preq(1,2,q)==0) then
          write(ifll,"(a,$)") "<-- NOT consider"
        else
          write(ifll,"(e12.3,$)") (preq(i,2,q),i=1,n)
          if(UNIT_SI(q)==0) then
            write(ifll,"(2x,a,$)")
     &      'in SI (A=mole-m-K ; e=J/mole) unit'
          else
            write(ifll,"(2x,a,$)")
     &      'in (cm-sec-mole-kcal-Kelvins)'
          endif
          if(lreq(1,2,q)/=undef) call disp_Lindemann(2,q,n)
        endif
        if(M_3rd(q)>0) then
          call disp_M(q)
!          write(ifll,"(/,7x,a,$)") "coefficients of M: "
!          write(ifll,"(f5.1,$)") (mreq(s,q),s=1,ns_gs)
        endif
      endif
      end subroutine disp_coefficients
!
!=================================================
      subroutine disp_M(q)
!=================================================
      integer,intent(in) :: q  ! classification number of equation
      character(60) str_M
      integer :: s
!
      write(ifll,"(/,7x,a,$)") "coefficients of M: "
      do s=1,ns,10
        str_M = ""
        do i=1,10
          if(s+i-1>ns) exit
          write(str_M(i+4*(i-1):i+4*i),'(f5.1)') mreq(s+i-1,q)
        enddo
        if(s>9) then
          write(ifll,'(/,26x,a,$)') trim(str_M)
        else
          write(ifll,'(a,$)') trim(str_M)       ! 3rd body coefficient
        endif
      enddo
      end subroutine disp_M
!-------------------------------------------------
! --- display 1 equation's side
!=================================================
      subroutine str_eqs(isw,q,sinp)
!=================================================
        integer,intent(in) :: isw ! switch
        integer,intent(in) :: q   ! classification number of equation
        character(*),intent(inout):: sinp
!
        character(4) :: strrr
        integer :: lenc,lensum,idum,idum2
        integer :: s
!
        sinp = ""
        k=0
        lensum=0
        do s=1,ns
          if(lensum>80) call 
     &   FFRABORT(1,'ERR: in str_eqs, call your supportor')
          if(vreq(s,isw,q)>0 ) then
            k=k+1
            if(k/=1) then
               idum=lensum+1
               lensum=lensum+3
               sinp(idum:lensum)=' + '
            endif
            write(strrr,'(F4.2)') vreq(s,isw,q)
            idum=lensum+1
            lensum=lensum+4
            sinp(idum:lensum)=strrr
            lenc=len_trim(spcnam(s))
            idum=lensum+1
            lensum=lensum+lenc
            sinp(idum:lensum)=adjustl(spcnam(s))
          endif
        enddo

        if(M_3rd(q)==1) then
          idum=lensum+1
          lensum=lensum+4
          sinp(idum:lensum)=' + M'
        endif
      end subroutine str_eqs
!
!=================================================
      character(5) function arrow(arg)  ! display arrow
!=================================================
        real(8),intent(in) :: arg
        if(arg/=0) then    ! reversible reaction
          arrow = " <=> "
        else                 ! irreversible reaction
          arrow = "  => "
        endif
      end function arrow
!=================================================
      subroutine disp_Lindemann(arg,q,n)	! display Lindemann form
!=================================================
      integer, intent(in) :: arg	! reversible or irreversible
      integer, intent(in) :: q,n
      if(kind_prs(q)==falloff) then
        write(ifll,"(/,7x,a,$)") "low pressure(uni): "
      else
        write(ifll,"(/,7x,a,$)") "low pressure(bi) : "
      endif
      write(ifll,"(es12.3,$)") (lreq(i,arg,q),i=1,n)
      if(troe(1,arg,q)/=undef .and. troe(1,arg,q)/=0.d0) then	! Troe form
        write(ifll,"(/,7x,a,$)") "Troe form        : "
        if(troe(4,arg,q)/=undef) then
          write(ifll,"(es12.3,$)") (troe(i,arg,q),i=1,4)
        else
          write(ifll,"(es12.3,$)") (troe(i,arg,q),i=1,3)
          write(ifll,"(a,$)") "  no T**"
        endif
        if(0.d0==troe(2,arg,q) .or. 0.d0==troe(3,arg,q)) then
          write(ifll,"(/,2a,/,2a,$)")
     &      "order of Troe form      : ",
     &      "troe_a, troe_t3, troe_t1, troe_t2",
     &      "INPUT ERROR:",
     &      "troe_t3 and troe_t1 must not be 0."
          ierror = ierror+1
        endif
      endif
      end subroutine disp_Lindemann
!
      end subroutine inputdata
!
!-------------------------------------------------
! difference of
! (Gibbs energy on equation)/((universal gas constant)*(temperature)) [-]
!=================================================
      real(8) function dg_RT(t,q)
!=================================================
        use module_species,only : judg_temp
        real(8), intent(in) :: t   ! temperature [K]
        integer,intent(in)  :: q   ! classification number of equation
        integer jt                 ! judgment of temperature
        jt = judg_temp(1,t)
        dg_RT = vreq_a7(1,jt,q)*(1.d0-log(t))
     &    +t*( -vreq_a7(2,jt,q)*0.5d0
     &    +t*( -vreq_a7(3,jt,q)/6.d0
     &    +t*( -vreq_a7(4,jt,q)/12.d0
     &    +t*( -vreq_a7(5,jt,q)*0.05d0))))
     &         +vreq_a7(6,jt,q)/t
     &         -vreq_a7(7,jt,q)
      end function dg_RT
!
!--------------------------------------------------------------------------
! --- difference of (Gibbs energy on equation)/(universal gas constant) [K]
!=================================================
      real(8) function dg_t(ncom,t,q)
!=================================================
        use module_species,only : h_t,s_0
        integer,intent(in) :: ncom  ! number of species
        real(8),intent(in) :: t   ! temperature [K]
        integer,intent(in) :: q   ! classification number of equation
        integer :: s              ! classification number of species
        dg_t=0.d0
        do s=1,ncom
          if( sub_vreq(s,q)==0 ) cycle
          dg_t=dg_t+sub_vreq(s,q)*(h_t(s,t)-t*s_0(s,t))
        enddo
      end function dg_t
!
!-------------------------------------------------
! --- revised mol concentration by 'mreq' [mol/m3]
!=================================================
      real(8) function revised_mc(ncom,q,mc) 
!=================================================
        integer,intent(in) :: ncom       ! number of species
        integer,intent(in) :: q          ! classification number of equation
        real(8), intent(in) :: mc(ncom)  ! mol concentration [mol/m3]
        integer s                        ! classification number of species
        revised_mc = 0.d0
        do s=1,ncom                      ! revised mol concentration
                                         ! by 'mreq' [mol/m3]
          revised_mc = revised_mc+mreq(s,q)*mc(s)
        enddo
      end function revised_mc
!
      end module module_chemreac
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine modutl_chkset(id,jset,iset,sbnm)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
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
      write(*,'(1x,2a,2I10)') 'id,jset,iset= ',id,jset,iset
      CALL FFRABORT(1,'modutl_setnow')
      end subroutine modutl_chkset
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine modutl_setnow(time,nset,nval,nwdt,iwd,wdt,wd)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
! 1. Interpolate values of variables with time
!    at present time from time tables
!
! --- [dummy arguments]
!
      real(8) ,intent(in)    :: time
      integer,intent(in)    :: nset,nval,nwdt
      integer,intent(in)    :: iwd(0:nset)
      real(8) ,intent(in)    :: wdt(0:nval,nwdt)
      real(8) ,intent(inout) :: wd (0:nval,nset)
!
! --- [local entities]
!
      integer :: i,k,n,ks,ke,kx,k1,k2,ncyc
      real(8)  :: tim0,tim1,tim2,tcyc,fc1,fc2
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
        fc1=tim0/tcyc
        ncyc=int(fc1)
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
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_gf
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!< module for GF >
      implicit none
      integer  :: IUT0,IUT5,IUT6,IWRITE,IDIM,INAME,IRESV,ICHECK
      DATA IUT0    / 0 /
      DATA IUT5    / 5 /
      DATA IUT6    / 6 /
      DATA IWRITE  / 2 /
      DATA INAME   / 1 /
      DATA IDIM    / 3 /
      DATA IRESV   / 3 /
      DATA ICHECK  / 99 /
C
      INTEGER IACT,NCOMFL,NCOMST,MCOM
C
      PARAMETER ( MCOM = 20 )
      CHARACTER*60 COMFLE(MCOM)
      CHARACTER*60 COMSET(MCOM)
!
      end module module_gf
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_usersub
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_usersub)'
!
!
! --- [module arguments]
!
!
      integer,parameter,public :: usrno=0,usryes=1,usryesyes=2
      integer,public :: outusr,inlusr,iniusr,oulusr,walusr,stcusr
      integer,public :: src_uvw,src_r,src_t,src_rans,src_fire
      integer,public :: BC_MHD_user,ini_MHD_user,CURT_MHD_user
      integer,public :: iniusr_E2P,src_BC_E2P,src_BC_CAVI
      real(8) ,save  :: usr_dum1,usr_dum2,usr_dum3,
     &                  usr_dum4,usr_dum5,usr_dum6
      integer,public :: src_rad,para_rad
      integer,public :: src_p_ini,sata_ruf_wall,user_rgin
!
! usryes  : user subroutine is used by user
! usrno   : user subroutine is NOT used by user
!
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,firewal,imhd,ierror)
!=================================================
      use module_model,only : u_func 
      implicit none
!
! --- [dummy qrgument]
!
      integer,intent(in)  :: ifli,ifll,ifle,imhd
      character(*),intent(in) :: cntlnam
      logical,intent(in)  :: firewal
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=6),parameter :: 
     &                     cond(3)=(/'yes   ','no    ','confrm'/)
      integer :: iset=0,co(3)=(/1,0,2/)
      integer :: ios,ierr1
!
! --- [namelist]
!
      integer       :: kout=1,kinl=1,ksrcu=1,kini=1,ksrcr=1,ksrct=1,
     &                 ksrcys=1,ksrcrans=1,kglb=1,kfire=1,koul=1,
     &                 kwal=1,kMHD_i=1,kMHD_bc=1,kMHD_c=1,kstc=1,
     &                 kiniE2P=1,kBC_E2P=1,kBC_CAVI=1,
     &		       srcrad=1,pararad=1,kprt_ini=1,sata_wall=1,
     &                 User_region=1
      character(20) :: output,initial,initial_E2P,source_BC_E2P
      character(20) :: source_BC_CAVI
      character(20) :: inlet,outlet,wall,static
      character(20) :: source_uvw,
     &                 source_r,source_t,source_rans,
     &                 source_fire,MHD_initial,MHD_BC,MHD_CURRENT
      character(20) :: radpro,radpara
      character(20) :: initial_particle,sata_rough_wall
      character(20) :: user_defined_region
      namelist /usrsub/output,inlet,initial,source_uvw,initial_E2P,
     &                 source_r,source_t,source_rans,
     &                 outlet,wall,static,source_BC_E2P,
     &                 source_fire,source_BC_CAVI,
     &                 MHD_initial,MHD_BC,MHD_CURRENT,
     &		       radpro,radpara,
     &                 initial_particle,sata_rough_wall,
     &                 user_defined_region
!     &                 shutter
!
      ierr1=0 
!-----------------------
! --- initialize flags
!-----------------------
      outusr=0
      inlusr=0
      oulusr=0
      walusr=0
      stcusr=0
      src_uvw=0
      iniusr=0
      iniusr_E2P=0
      src_r=0
      src_t=0
      src_rans=0
      src_fire=0
      src_BC_E2P=0
      src_BC_CAVI=0
      src_rad=0
      src_p_ini=0
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
      output=' '
      inlet=' '
      outlet=' '
      wall=' '
      static=' '
      source_uvw=' '
      initial=' '
      initial_E2P=' '
      initial_particle=' '
      sata_rough_wall=' '
      source_BC_E2P=' '
      source_BC_CAVI=' '
      source_r=' '
      source_t=' '
      source_rans=' '
      source_fire=' '
      MHD_initial=' '
      MHD_BC=' '
      MHD_CURRENT=' '
!
      rewind ifli
      read(ifli,usrsub,iostat=ios)
      if( ios.lt.0 ) return
!
      call nml_errmsg0(ifle,ios,'usrsub',ierr1)
      if( ierr1.ne.0 ) call FFRABORT(1,'module_usersub/inputdata')
!
      call nml_listno(3,cond,output,kout)
      call nml_chkchr0(ifle,'output',output,kout,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [output] into [&usrsub]')
!
      call nml_listno(3,cond,inlet,kinl)
      call nml_chkchr0(ifle,'inlet',inlet,kinl,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [inlet] into [&usrsub] ')
!
      call nml_listno(3,cond,outlet,koul)
      call nml_chkchr0(ifle,'outlet',outlet,koul,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [inlet] into [&usrsub] ')
!
      call nml_listno(3,cond,wall,kwal)
      call nml_chkchr0(ifle,'wall',wall,kwal,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [inlet] into [&usrsub] ')
!
      call nml_listno(3,cond,source_uvw,ksrcu)
      call nml_chkchr0(ifle,'source_uvw',source_uvw,ksrcu,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [inlet] into [&usrsub] ')
!
      call nml_listno(3,cond,initial,kini)
      call nml_chkchr0(ifle,'initial',initial,kini,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [initial] into [&usrsub] ')
!
      call nml_listno(3,cond,initial_E2P,kiniE2P)
      call nml_chkchr0(ifle,'initial_E2P',initial_E2P,kiniE2P,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [initial_E2P] into [&usrsub] ')
!
      call nml_listno(3,cond,initial_particle,kprt_ini)
      call nml_chkchr0
     &  (ifle,'initial_particle',initial_particle,kprt_ini,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [initial_particle] into [&usrsub] ')
!
      call nml_listno(3,cond,sata_rough_wall,sata_wall)
      call nml_chkchr0
     &  (ifle,'sata_rough_wall',sata_rough_wall,sata_wall,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [sata_rough_wall] into [&usrsub] ')
!
      call nml_listno(3,cond,user_defined_region,User_region)
      call nml_chkchr0
     &(ifle,'user_defined_region',user_defined_region,User_region,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [user_defined_region] into [&usrsub] ')
!

      call nml_listno(3,cond,source_BC_E2P,kBC_E2P)
      call nml_chkchr0(ifle,'source_BC_E2P',source_BC_E2P,kBC_E2P,ierr1)
      if( ierr1.ne.0 )
     &  call FFRABORT(1,'add [source_BC_E2P] into [&usrsub] ')
!
      call nml_listno(3,cond,source_BC_CAVI,kBC_CAVI)
      call nml_chkchr0
     &    (ifle,'source_BC_CAVI',source_BC_CAVI,kBC_CAVI,ierr1)
      if( ierr1.ne.0 )
     &  call FFRABORT(1,'add [source_BC_CAVI] into [&usrsub] ')
!
      call nml_listno(3,cond,source_r,ksrcr)
      call nml_chkchr0(ifle,'source_r',source_r,ksrcr,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [source_r] into [&usrsub]')
!
      call nml_listno(3,cond,source_t,ksrct)
      call nml_chkchr0(ifle,'source_t',source_t,ksrct,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [source_t] into [&usrsub]')
!
      call nml_listno(3,cond,source_rans,ksrcrans)
      call nml_chkchr0(ifle,'source_rans',source_rans,ksrcrans,ierr1)
      if( ierr1.ne.0 )  
     &  call FFRABORT(1,'add [source_rans] into [&usrsub]')
!
      call nml_listno(3,cond,source_fire,kfire)
      call nml_chkchr0(ifle,'source_fire',source_fire,kfire,ierr1)
      if(ierr1.ne.0)  
     &  call FFRABORT(1,'add [source_fire] into [&usrsub]')
!
      call nml_listno(3,cond,static,kstc)
      call nml_chkchr0(ifle,'static',static,kstc,ierr1)
      if(ierr1.ne.0)  
     &  call FFRABORT(1,'add [static] into [&usrsub]')
!
      IF(imhd>0) then
        call nml_listno(3,cond,MHD_initial,kMHD_i)
        call nml_chkchr0(ifle,'MHD_initial',MHD_initial,kMHD_i,ierr1)
        if(ierr1.ne.0)  
     &  call FFRABORT(1,'add [MHD_initial] into [&usrsub]')
!
        call nml_listno(3,cond,MHD_BC,kMHD_bc)
        call nml_chkchr0(ifle,'MHD_BC',MHD_BC,kMHD_bc,ierr1)
        if(ierr1.ne.0)  
     &  call FFRABORT(1,'add [MHD_BC] into [&usrsub]')
!
        call nml_listno(3,cond,MHD_CURRENT,kMHD_c)
        call nml_chkchr0(ifle,'MHD_CURRENT',MHD_CURRENT,kMHD_c,ierr1)
        if(ierr1.ne.0)  
     &    call FFRABORT(1,'add [MHD_CURRENT] into [&usrsub]')
!
      endif
!
      call nml_listno(3,cond,radpro,srcrad)
      call nml_chkchr0(ifle,'radpro',radpro,srcrad,ierr1)
      if(ierr1.ne.0)  
     &  call FFRABORT(1,'add [radpro] into [&usrsub]')
      call nml_listno(3,cond,radpara,pararad)
      call nml_chkchr0(ifle,'radpara',radpara,pararad,ierr1)
      if(ierr1.ne.0)  
     &  call FFRABORT(1,'add [radpara] into [&usrsub]')
!----------------------
! --- define flags
!----------------------
      outusr=co(kout)
      inlusr=co(kinl)
      oulusr=co(koul)
      walusr=co(kwal)
      stcusr=co(kstc)
      iniusr=co(kini)
      iniusr_E2P=co(kiniE2P)
      src_p_ini=co(kprt_ini)
      sata_ruf_wall=co(sata_wall)
      user_rgin=co(User_region)

      src_uvw=co(ksrcu)
      src_r=co(ksrcr)
      src_t=co(ksrct)
      src_rans=co(ksrcrans)
      src_BC_E2P=co(kBC_E2P)
      src_BC_CAVI=co(kBC_CAVI)
      src_fire=co(kfire)
      
      IF(imhd>0) then
        BC_MHD_user=co(kMHD_bc)
        ini_MHD_user=co(kMHD_i)
        CURT_MHD_user=co(kMHD_c)
      endif
      if(firewal) then
        if(src_fire.eq.usrno) then
          write(ifle,*) ' ### ERR: fire wall must use user subroutine'
          write(ifle,*) " ### MSG: set [source_fire='yes'] "
          write(ifle,*) " ### MSG: compile user sub. [user_src_fire]"
          call FFRABORT(1,' module_usersub')
        endif
      endif
!
      if(u_func(2)==1.or.u_func(3)==1) then
        if(user_rgin/=usryes) then
          write(ifle,'(1X,a)') 
     &  'ERR: u_func(2)==1.or.u_func(3)==1 MUST set user define region'
          write(ifle,'(1X,a)') 'MSG: &usrsub/user_defined_region="yes"'
          write(ifle,'(1X,a)') 'MSG: Edit USER_UserDefine_region.f'
          call 
     & FFRABORT(1,'ERR-MSG: set [&usrsub/user_defined_region="yes"] ')
        endif
      endif
!

      src_rad=co(srcrad)
      para_rad=co(pararad)
!
      return
!
 9999 continue
      ierror=1
      write(ifle,*) modnam,subnam
      end subroutine inputdata
!
      end module module_usersub
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_hpc
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      integer :: MNVX_GL
      character(20),parameter,private :: modnam='(module_hpc)'
!
!///////////////////////////////////////////////////////////////////////
      contains
!============================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,NPE,ierror)
!============================================================
      implicit none
!
! --- [dummy qrgument]
!
      integer,intent(in)  :: ifli,ifll,ifle,NPE
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      integer :: ierr1,ios
      integer :: iset=0
!
! --- [namelist]
!
      logical :: lghpc
      integer :: ncpu,MONITOR
      namelist /hpc/ ncpu,MONITOR
!
!
      ierr1=0
      ierror=0
!
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
      ncpu=1
      MONITOR=1
      MNVX_GL=1
!
      rewind ifli
      read(ifli,hpc,iostat=ios)
      call nml_errmsg0(ifle,ios,'hpc',ierr1)
      if( ierr1.ne.0 ) 
     &   call FFRABORT(1,'ERR: Reading error in &hpc')
!
      lghpc=NPE==ncpu
      if(NPE.gt.1) then
         call hpcland(lghpc,1)
      endif
!      
      if(.not.lghpc) then   !NPE.ne.ncpu) then
        write(ifle,*) ' ### ERR : ncpu must be equal to',NPE
        call FFRABORT(1,'ERR: in module_hpc/inputdata')
      endif
!
      MNVX_GL=MONITOR
      if(MONITOR.lt.1) then
        write(ifle,*) ' ### ERR : MONITOR CELL NUMBER have be > 0',NPE
        call FFRABORT(1,'module_hpc/inputdata')
      endif
!
      return
 9999 continue
      ierror=1
      write(ifle,*) modnam,subnam
      end subroutine inputdata
      end module module_hpc
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_flow_sound1
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_flow_sound)'
!
! --- [module arguments]
!
      logical,public           :: lsound,lshear
      real(8) ,save,allocatable :: FX(:),FY(:),FZ(:)
      integer,save,allocatable :: isound(:),iname(:)
      character*80,allocatable :: names(:)
      integer :: nsond
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata
     & (ifli,ifll,ifle,cntlnam,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,ierror)
!=================================================
      implicit none
!
!
! --- [dummy qrgument]
!
      integer,intent(in) :: ifli,ifll,ifle,nbcnd,kxnone,kdfire,kdintr
     &                     ,kdbuff,kdcvd
      character(*),intent(in) :: cntlnam
      integer,intent(in)  :: kdbcnd(0:7,nbcnd)
      character(*),intent(in) :: boundName(nbcnd,2)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      character(LEN=3),parameter :: cond(2)=(/'yes','no '/)
      character(LEN=7),parameter :: co(2)=(/'.true. ','.false.'/)
      character(LEN=11),parameter :: subnam='(inputdata)'
      integer :: ierr1,ios,nb,nbs,kd
      integer :: iset=0
      integer :: ksound=0
!
! --- [namelist]
!
      character(20) :: flow_sound
      character(80) :: name(50)
      namelist /sound/ flow_sound,name
!
!
      ierr1=0
      ierror=0
!
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
      lsound=.false.
      flow_sound=' '
      name(:)=' '
!
      rewind ifli
      read(ifli,sound,iostat=ios)
!
      call nml_errmsg0(ifle,ios,'sound',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      call nml_listno(2,cond,flow_sound,ksound)
      call nml_chkchr0(ifle,'flow_sound',flow_sound,ksound,ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      if(ksound.eq.1) then
        lsound=.true.
      elseif(ksound.eq.2) then
        lsound=.false.
      else
        write(ifle,*) "ERR: flow_sound='yes' or 'no' in ",cntlnam
        call FFRABORT(1,'module_flow_sound/inputdata')
      endif
      if(.not.lsound) return
!
      do nb=1,nbcnd
      do nbs=nb+1,nbcnd
      if(name(nb).eq.name(nbs).and.name(nb).ne.' ') then
        write(ifle,*) 'The name : ',name(nbs)(:len_trim(name(nbs))),
     &  ' has been duplicated'
        call FFRABORT(1,'module_flow_sound/inputdata')
      endif
      enddo
      enddo
!
      nsond=0
      do nb=1,nbcnd
      if(name(nb).ne.' ') then
        nsond=nsond+1
      endif
      enddo
!
      allocate(names(nsond),iname(nsond))
      nsond=0
      names(:)=' '
      do nb=1,nbcnd
      if(name(nb).ne.' ') then
        nsond=nsond+1
        names(nsond)=name(nb)
      endif
      enddo
!
      allocate(isound(nbcnd))
      isound(:)=0
      iname(:)=0
      do nbs=1,nsond
      do nb=1,nbcnd
      if(adjustl(names(nbs)).eq.adjustl(boundName(nb,1))) then
        isound(nb)=1
        iname(nbs)=nb
      endif
      enddo
      enddo
!
      do nbs=1,nsond
       if(adjustl(names(nbs)).ne.' '.and.iname(nbs).eq.0) then
         write(ifle,*) '### error : namelist [sound]'
         write(ifle,*)
     &                 names(nbs)(:len_trim(names(nbs))),
     &                 ' not found in boundary namelist'
       call FFRABORT(1,'module_flow_sound/inputdata')
       endif
      enddo
!
      do nbs=1,nsond
       nb=iname(nbs)
       kd=kdbcnd(0,nb)
       if(.not.(kd.eq.kxnone.or.
     &          kd.eq.kdfire.or.
     &          kd.eq.kdintr.or.
     &          kd.eq.kdbuff.or.
     &          kd.eq.kdcvd)) then
         write(ifle,*)
     &   ' ERR: The boundary : ',names(nbs)(:len_trim(names(nbs))),
     &   ' is NOT wall boundary '
         call FFRABORT(1,'module_flow_sound/inputdata')
       endif
      enddo
!
!
      allocate(FX(nsond+1),FY(nsond+1),FZ(nsond+1))
!
      return
 9999 continue
      ierror=1
      write(ifle,*) modnam,subnam
      end subroutine inputdata
      end module module_flow_sound1
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_flow_sound
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_flow_sound)'
!
! --- [module arguments]
!
      logical                  :: soundObservFlag
      integer                  :: numofSoundWall
      character*80,allocatable :: nameofSoundWall(:)
      integer,allocatable      :: soundWall_bcNo(:)
      real(8),allocatable      :: soundWallCenterPosit(:,:)
      integer                  :: numofSoundObserver
      real(8),allocatable      :: soundObservPosit(:,:)
      integer,allocatable      :: soundModelCode(:)
      real(8),allocatable      :: prsPrevSave(:,:),machPrevSave(:,:)
      real(8),allocatable      :: dflxPrevSave(:,:)
      real(8),allocatable      :: totalSFaceArea(:)
      real(8)                  :: timePrevSave
!
!     soundObservFlag : Flag of sound calculation enable/disable.
!     numofSoundWall  : Number of sound souce walls.
!     nameofSoundWall : List of sound wall name.
!     soundWall_bcNo  : Boundary region index number of sound souce wall.
!     soundWallCenterPosit : Position of center of sound source wall.
!     numofSoundObserver : Number of observer points.
!     soundObservPosit: Positions of sound oberver points.
!     soundModelCode  : Sound model code No. of the observer points.
!     prsPrevSave     : Buffer for previous time step data.
!     dflxPrevSave    : Buffer for previous time step data.
!                       (Using for calculation of d/dtau)
!     timePrevSave    : Time of the previous step.
!
!///////////////////////////////////////////////////////////////////////
      contains
!=======================================================================
      subroutine inputdata
     & (ifli,ifll,ifle,cntlnam,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,ierror)
!=======================================================================
      implicit none
!
! --- [dummy qrgument]
!
      integer,intent(in)      :: ifli,ifll,ifle,
     &                           nbcnd,kxnone,kdfire,kdintr,kdbuff,
     &                           kdcvd
      character(*),intent(in)    :: cntlnam
      integer,intent(in)      :: kdbcnd(0:7,nbcnd)
      character(*),intent(in) :: boundName(nbcnd,2)
      integer,intent(out)     :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=3),parameter :: cond(2)=(/'yes','no '/)
      character(LEN=7),parameter :: co(2)=(/'.true. ','.false.'/)
      character(LEN=6),parameter :: modelLabel(2)=(/'Curle ','FWH   '/)
      integer :: ierr1,ios
      integer :: iset=0
      integer :: ksound=0

!
      integer :: ic,jc,kc,nb,kd
      character*80,allocatable :: nameofSoundWallTmp(:)
!
! --- [namelist]
!
      character(20) :: enable_source,enable_observer
      character(80) :: wall_name
      real(8)       :: position_x,position_y,position_z
      character*80  :: model
      namelist /sound_source/   enable_source,wall_name
      namelist /sound_observer/ enable_observer,position_x,position_y,
     &                          position_z,model
!
      ierr1=0
      ierror=0
      soundObservFlag=.true.
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!--------------------------------------------------------
! --- Sound wall setting FROM HERE ---------------------
! --- Count up the number of enabled sound source walls.
!--------------------------------------------------------
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
        endif
        call nml_listno(2,cond,enable_source,ksound)
        if(ksound==1) numofSoundWall=numofSoundWall+1
      end do
!
      if(numofSoundWall<=0) then
        soundObservFlag=.false.
        return
      end if
!-------------------------------------
! --- Prepare the list of wall name.
!-------------------------------------
      allocate(nameofSoundWallTmp(1:numofSoundWall))
!----------------------------------      
! --- Rereading the namelist file.
!----------------------------------
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
        endif
        call nml_listno(2,cond,enable_source,ksound)
        if(ksound==1) then
          ic=ic+1
          nameofSoundWallTmp(ic)=wall_name
        endif
      enddo
!--------------------------------------------------------------------
! --- Clear duplicated wall names and recount the number of walls.
!--------------------------------------------------------------------
      if(ic>=2) then
        do jc=1,ic-1
          if(nameofSoundWallTmp(jc)/=' ') then
            do kc=jc+1,ic
              if(nameofSoundWallTmp(kc)==nameofSoundWallTmp(jc))
     &            nameofSoundWallTmp(kc)=' '
            end do
          end if
        end do
!
        numofSoundWall=0
        do jc=1,ic
          if(nameofSoundWallTmp(jc)/=' ') 
     &            numofSoundWall=numofSoundWall+1
        enddo
      endif
!
      allocate(nameofSoundWall(1:numofSoundWall))
!
      if(ic>=2) then
        kc=0
        do jc=1,ic
          if(nameofSoundWallTmp(jc)/=' ') then
            kc=kc+1
            nameofSoundWall(kc)=nameofSoundWallTmp(jc)
          endif
        enddo
      else
        nameofSoundWall(:)=nameofSoundWallTmp(:)
      endif
      deallocate(nameofSoundWallTmp)
!--------------------------------------------------
! --- Search relation between sound source walls 
! --- and boundary region index number.
!--------------------------------------------------
      allocate(soundWall_bcNo(1:numofSoundWall))
      soundWall_bcNo(:)=0
      do jc=1,numofSoundWall
        do nb=1,nbcnd
          if(adjustl(nameofSoundWall(jc))==adjustl(boundName(nb,1))) 
     &      soundWall_bcNo(jc)=nb
        end do
      end do
!-------------------------
! --- Check the list.
!-------------------------
      do jc=1,numofSoundWall
        if(soundWall_bcNo(jc)==0) then
          write(ifle,*)'### error : namelist [sound]'
          write(ifle,*)' MSG: BC name: ', trim(nameofSoundWall(jc)),
     &                 ' not found in boundary namelist'
          call FFRABORT(1,'module_flow_sound/inputdata')
        endif
      end do
!------------------------------------------------------------
! --- Confirm the boundary conditions of sound source walls.
!------------------------------------------------------------
      do jc=1,numofSoundWall
        nb=soundWall_bcNo(jc)
        kd=kdbcnd(0,nb)
        if(.not.(kd.eq.kxnone.or.
     &           kd.eq.kdfire.or.
     &           kd.eq.kdintr.or.
     &           kd.eq.kdbuff.or.
     &           kd.eq.kdcvd)) then
          write(ifle,*)
     &   ' ERR: The boundary : ',trim(nameofSoundWall(jc)),
     &   ' is NOT wall boundary'
          call FFRABORT(1,'module_flow_sound/inputdata')
        endif
      enddo
!      
      allocate(soundWallCenterPosit(1:numofSoundWall+1,1:3))
      allocate(totalSFaceArea(1:numofSoundWall+1))
!-----------------------------------------------------------------
! --- Sound wall setting TO HERE ---------------------------------
!      
! --- Observer points setting FROM HERE --------------------------
! --- Count up the number of enabled obervers points
!-----------------------------------------------------------------
!
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
!--------------------------------------------------
! --- Even if there is no sound_oberver namelist,
!     Curle's equation is calculated.
!--------------------------------------------------
      numofSoundObserver=numofSoundObserver+1
!      
! --- Preare arrays of observer points data.
!
      allocate(soundObservPosit(1:numofSoundObserver,1:3))
      allocate(soundModelCode(1:numofSoundObserver))
!
! --- Reload name list file.
      rewind(ifli)
!-------------------------------------------------------
! --- Data for default calculation of Curles equation
!-------------------------------------------------------
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
            write(ifle,*)'ERR: Unkown Sound Model : ',trim(model)
            write(ifle,*) modnam,subnam
            return
          end if
          soundModelCode(ic)=kc
        end if
      end do
!--------------------------------------
! --- Observer points setting TO HERE 
!--------------------------------------
      end subroutine inputdata
      end module module_flow_sound
!
!zhang0119
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_fluidforce
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_fluidforce)'
!
! --- [module arguments]
!
      logical                  :: ForceFlag
      integer,save             :: Nforce
      character*80,allocatable :: NM_F(:)
      character*80,allocatable :: NM_F_WALL(:)
      integer,save,allocatable :: NO_F_WALL(:)
      integer,save,allocatable :: F_NUM(:),iwallvar(:)
      integer,save             :: total_F
      real(8),allocatable      :: end_mnt(:,:),begin_mnt(:,:),
     &                            unit_m(:,:)
!
      real(8),allocatable :: F_fric(:,:),F_pres(:,:),F_SUM(:,:)
      real(8),allocatable :: M_SUM(:,:)
!
!     ForceFlag   : Flag of fulid force 
!     Nforce      : Number of fulid force
!     NM_F : name of force (=FORCE_NO_1, when no=1) 
!     NM_F_WALL : wall-name in force-name
!     NO_F_WALL : wall-BC-number in force-name
!     F_NUM     : wall number in Force-No
!     total_F   : total force-wall number
!     end_mnt   : axis end point cood. of Moment-No
!     begin_mnt : axis beging point cood. of Moment-No
!
!///////////////////////////////////////////////////////////////////////
      contains
!=======================================================================
      subroutine inputdata
     & (ifli,ifll,ifle,cntlnam,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,kvnslp,kvlglw,kvfslp,kvmodl,kvEWF,
     &  ierror)
!=======================================================================
      implicit none
!
! --- [dummy qrgument]
!
      integer,intent(in)      :: ifli,ifll,ifle,
     &                           nbcnd,kxnone,kdfire,kdintr,kdbuff,
     &                           kdcvd,kvnslp,kvlglw,kvfslp,kvmodl,kvEWF
      character(*),intent(in) :: cntlnam
      integer,intent(in)      :: kdbcnd(0:7,nbcnd)
      character(*),intent(in) :: boundName(nbcnd,2)
      integer,intent(out)     :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=3),parameter :: cond(2)=(/'yes','no '/)
      character(LEN=7),parameter :: co(2)=(/'.true. ','.false.'/)
      character(LEN=5),parameter :: cfile(2)=(/'GF   ','Excel'/)
      integer :: ierr1,ios
      integer :: iset=0
      integer :: kfforce=0
!
      integer :: ic,jc,ib,kc,nb,kdv,nbx
!
! --- [namelist]
!
      integer,parameter :: MX_F_M=50
      character(80)     :: wall_name(1:MX_F_M),label
!      character(20)     :: filetype
      real*8            :: end_x,end_y,end_z,
     &                     begin_x,begin_y,begin_z
      integer,save      :: F_All(1:MX_F_M)=0
      real   ,parameter :: undef=-huge(1.)
      integer           :: IDUM1,IDUM2,I,II,NUMWALL
      integer           :: no
      character(8)      :: text
      integer           :: enable
      namelist /force_fluid/ no,enable,wall_name,label,
     &                       end_x,end_y,end_z,
     &                       begin_x,begin_y,begin_z
      real*8            :: dum1
!
      ierr1=0
      ierror=0
      ForceFlag=.false.
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!-----------------
! --- wall number
!-----------------
      NUMWALL=0
      do nb=1,nbcnd
      kdv=kdbcnd(1,nb)
      if(kdv==kvnslp.or.kdv==kvlglw.or.
     &    kdv==kvfslp.or.kdv==kvmodl.or.kdv==kvEWF) then
        NUMWALL=NUMWALL+1
      endif
      enddo
!-----------------
! read force & moment
!-----------------
      rewind(ifli)
      Nforce=0
      IDUM2=0
      do while(.true.)
        wall_name(:)=' '
        enable=-1
        ic=0
        read(ifli,force_fluid,iostat=ios)
        if(ios<0) exit
        if(enable==0) cycle
        if(enable==-1) then
          write(ifle,'(1x,2a)') 
     &   'ERR: Reset [enable] in [&force_fluid] in ',trim(cntlnam)
          write(ifle,'(1x,a)') 
     &    'MSG: [enable=0]: is inavailable'
          write(ifle,'(1x,a)') 
     &    'MSG: [enable=1]: is available'
          call FFRABORT(1,'ERR: in &force_fluid')
        endif
        call nml_errmsg0(ifle,ios,'force_fluid',ierr1)
        if(ierr1/=0) then
          ierror=1
          write(ifle,'(1x,a)') 
     &       'ERR: reading [&force_fluid] error'
          call FFRABORT(1,'ERR: in &force_fluid')
        endif
        ic=count(wall_name(:)/=' ')
        if(ic>0) then
          Nforce=Nforce+1
          ForceFlag=.true.
          if(TRIM(adjustl(wall_name(1)))=='ALL_WALL') then
            F_All(Nforce)=1
            IDUM2=IDUM2+NUMWALL
          else
            IDUM2=IDUM2+ic
          endif
        endif
      enddo
      total_F=IDUM2
!
      allocate(iwallvar(1:nbcnd))
      iwallvar(:)=1
!      
      if(.NOT.ForceFlag) return
!
! --- 
!
      allocate(end_mnt(3,max(1,Nforce)),begin_mnt(3,max(1,Nforce)))
      allocate(NM_F(max(1,Nforce)))
      allocate(F_NUM(0:max(1,Nforce)))
      allocate(NM_F_WALL(max(1,total_F)))
      allocate(NO_F_WALL(max(1,total_F)))
      allocate(F_fric(3,max(1,Nforce)))
      allocate(F_pres(3,max(1,Nforce)))
      allocate(F_SUM(3,max(1,Nforce)))
      allocate(M_SUM(4,max(1,Nforce)))
      allocate(unit_m(3,max(1,Nforce)))
      NM_F=' '
      end_mnt(:,:)=0.d0
      begin_mnt(:,:)=0.d0
      end_mnt(3,:)=1.d0
!--------------------
! --- re-read 
!--------------------
      Nforce=0
      F_NUM(:)=0
      rewind(ifli)
      do while(.true.)
        wall_name=' '
        label=' '
        no=0
        end_x=undef
        end_y=undef
        end_z=undef
        begin_x=undef
        begin_y=undef
        begin_z=undef
        read(ifli,force_fluid,iostat=ios)
        if(ios<0) exit
        ic=count(wall_name(:)/=' ')
        if(ic>0) then
          Nforce=Nforce+1
          F_NUM(Nforce)=F_NUM(Nforce-1)
          if(label==' ') then
            write(text,'(i8)') no
            NM_F(Nforce)='FORCE_NO_'//TRIM(text)
          else
            NM_F(Nforce)=trim(adjustl(label))
          endif
          if(F_All(Nforce)==1) then
            II=0
            do nb=1,nbcnd
            kdv=kdbcnd(1,nb)
            if(kdv==kvnslp.or.kdv==kvlglw.or.
     &         kdv==kvfslp.or.kdv==kvmodl.or.
     &         kdv==kvEWF)
     &      then
              II=II+1
              NM_F_WALL(II+F_NUM(Nforce))=
     &            trim(adjustl(boundName(nb,1)))
            endif
            enddo
            F_NUM(Nforce)=F_NUM(Nforce)+II
          else
            II=0
            do I=1,MX_F_M
            if(wall_name(I)/=' ')then
              II=II+1
              NM_F_WALL(II+F_NUM(Nforce))=trim(wall_name(I))
            endif
            enddo
            if(II/=ic) call FFRABORT(1,'ERR: Unknown error') 
            F_NUM(Nforce)=F_NUM(Nforce)+ic
            if(end_x==undef.or.end_y==undef.or.end_z==undef.or.
     &         begin_x==undef.or.begin_y==undef.or.begin_z==undef
     &       ) then
              write(ifle,'(1x,3a,I4)') 
     &        'ERR: NOT defined ',
     &        '[end_x,end_y,end_z,begin_x,begin_y,begin_z] in ',
     &        '[&force_fluid], no=',no
         call FFRABORT(1,'ERR: NOT defined [end_x]... or [begin_x]...')
            else
              dum1=(begin_x-end_x)**2
     &            +(begin_y-end_y)**2
     &            +(begin_z-end_z)**2
              if(dum1<1.D-15) then
               write(ifle,'(1x,a)') 'ERR: error in defining begin_x...'
                call FFRABORT(1,'ERR: chenk [&fluig_force]')
              else
                end_mnt(1,Nforce)=end_x
                end_mnt(2,Nforce)=end_y
                end_mnt(3,Nforce)=end_z
                begin_mnt(1,Nforce)=begin_x
                begin_mnt(2,Nforce)=begin_y
                begin_mnt(3,Nforce)=begin_z
              endif
            endif
          endif
        endif
      enddo
!-----------------------------
! --- 
!-----------------------------
      NO_F_WALL(:)=0
      do I=1,Nforce
        IDUM1=F_NUM(I-1)+1
        IDUM2=F_NUM(I)
        do II=IDUM1,IDUM2
          do nb=1,nbcnd
            kdv=kdbcnd(1,nb)
            if(kdv==kvnslp.or.kdv==kvlglw.or.kdv==kvfslp.or.
     &         kdv==kvmodl.or.kdv==kvEWF) 
     &      then
              if(trim(NM_F_WALL(II))==trim(adjustl(boundName(nb,1))))
     &        then
                NO_F_WALL(II)=nb
              endif
            else
              if(trim(NM_F_WALL(II))==trim(adjustl(boundName(nb,1))))
     &        then
                write(ifle,'(1x,4a)') 'ERR: NOT be Wall BC: ',
     &          trim(NM_F_WALL(II)),' in [&force_fluid] no=',
     &          trim(NM_F(I))
                call FFRABORT(1,'ERR: NOT wall BC name found')
              endif
            endif
          enddo
          if(NO_F_WALL(II)==0) 
     &    call FFRABORT(1,'ERR: NOT match BC name in [&force_fluid]')
        enddo
      enddo
!
      do I=1,Nforce
      unit_m(1,I)=end_mnt(1,I)-begin_mnt(1,I)
      unit_m(2,I)=end_mnt(2,I)-begin_mnt(2,I)
      unit_m(3,I)=end_mnt(3,I)-begin_mnt(3,I)
      dum1=dsqrt(unit_m(1,I)**2+unit_m(2,I)**2+unit_m(3,I)**2)
      unit_m(1,I)=unit_m(1,I)/dum1
      unit_m(2,I)=unit_m(2,I)/dum1
      unit_m(3,I)=unit_m(3,I)/dum1
      enddo
!
      end subroutine inputdata
      end module module_fluidforce
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine str2up(buf)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character, intent(inout):: buf*(*)
      integer :: i,n
      n=len(buf)
      do i=1,n
      if(ichar(buf(i:i))>=97 .and. ichar(buf(i:i))<=122)
     &        buf(i:i)=char(ichar(buf(i:i))-32)
      end do
      return
      end subroutine str2up
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine str2low(buf)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character, intent(inout):: buf*(*)
      integer :: i,n
      n=len(buf)
      do i=1,n
      if(ichar(buf(i:i))>=65 .and. ichar(buf(i:i))<=90)
     &        buf(i:i)=char(ichar(buf(i:i))+32)
      end do
      return
      end subroutine str2low
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_constant
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      integer,parameter :: ITWO=2,IZERO=0
!
      real(8) ,parameter :: SML=1.D-25
      real(8) ,parameter :: SML15=1.D-15
      real(8) ,parameter :: GREAT=1.D25
      real(8) ,parameter :: BIGG=1.D15
      real(8) ,parameter :: ZERO=0.D0
      real(8) ,parameter :: HALF=0.5D0
      real(8) ,parameter :: ONE=1.0D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0
      real(8) ,parameter :: FIVE=5.D0,SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0
      real(8) ,parameter :: NINE=9.D0,TEN=10.D0
      real(8) ,parameter :: pi=3.1415926d0
!
      end module module_constant
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_anim
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
      integer,save             :: ianim=0
      integer,save             :: ianim_uvw,ianim_p
      integer,save             :: ianim_r,ianim_t
      integer,save,allocatable :: ianim_rans(:),ianim_comp(:)
      integer,save,allocatable :: ianim_GAS_WDOT(:),
     &                            ianim_SDOT_suf(:),
     &                            ianim_molefr_suf(:),
     &                            ianim_depoeth_spd(:)
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,nrans,ncomp,mmcomp,
     &                     ncomp_suf,spinam,ierror)
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle,nrans,ncomp,ncomp_suf,mmcomp
      character(*),intent(in) :: cntlnam,spinam(mmcomp)
      integer,intent(out) :: ierror
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
      integer :: ios=0,ierr1=0,i
!
      ierror=0
      allocate(ianim_rans(nrans),ianim_comp(ncomp))
      allocate(ianim_GAS_WDOT(ncomp),
     &         ianim_SDOT_suf(ncomp+ncomp_suf),
     &         ianim_molefr_suf(ncomp+ncomp_suf),
     &         ianim_depoeth_spd(ncomp+ncomp_suf))
!
      anim_rans(:)=0
      anim_comp(:)=0
      anim_GAS_WDOT(:)=0
      anim_surface_SDOT(:)=0
      anim_surface_Fraction(:)=0
      anim_depos_etch_speed(:)=0
!
      ianim=0
      ianim_uvw=0
      ianim_p=0
      ianim_r=0
      ianim_t=0
      ianim_rans(:)=0
      ianim_comp(:)=0
      ianim_GAS_WDOT(:)=0
      ianim_SDOT_suf(:)=0
      ianim_molefr_suf(:)=0
      ianim_depoeth_spd(:)=0
!
      rewind ifli
      read(ifli,animation,iostat=ios)
      if( ios.lt.0 ) return
      call nml_errmsg0(ifle,ios,'animation',ierr1)
      if(ierr1.ne.0) call FFRABORT(1,'module_anim')
!
      ianim=anim
      ianim_uvw=anim_uvw
      ianim_p=anim_p
      ianim_r=anim_r
      ianim_t=anim_t
      ianim_rans(1:nrans)=anim_rans(1:nrans)
      ianim_comp(1:ncomp)=anim_comp(1:ncomp)
      ianim_GAS_WDOT(1:ncomp)=anim_GAS_WDOT(1:ncomp)
      do i=1,ncomp+ncomp_suf
      if(i<=ncomp.and.anim_surface_SDOT(i)==1) then
        ianim_SDOT_suf(i)=anim_surface_SDOT(i)
      else
        ianim_SDOT_suf(i)=anim_surface_SDOT(i)
      endif
      if(i<=ncomp.and.anim_surface_Fraction(i)==1) then
        write(ifle,'(2a)') 'ERR: Can NOT output animation for ',
     &  trim(spinam(i))
        write(ifle,'(a)') 
     6   'ERR: defining error in [anim_surface_Fraction]'
        CALL FFRABORT(1,'ERR:')
      else
        ianim_molefr_suf(i)=anim_surface_Fraction(i)
      endif
      if(i<=ncomp.and.anim_depos_etch_speed(i)==1) then
        write(ifle,'(2a)') 'ERR: Can NOT output animation for ',
     &  trim(spinam(i))
        write(ifle,'(a)') 
     6   'ERR: defining error in [anim_depos_etch_speed]'
        CALL FFRABORT(1,'ERR:')
      else
        ianim_depoeth_spd(i)=anim_depos_etch_speed(i)
      endif
      enddo
!
      end subroutine inputdata
      end module module_anim
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_Euler2ph
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     key-word for Euler 2 phase: ieul2ph
!
      implicit none
      character(21),parameter,private :: modnam='(module_Euler2ph)'
!
      integer,save :: ieul2ph=0,N_PHS=0
      integer,parameter :: IPHMAX=10
      integer,parameter :: kdphs_g=1,kdphs_l=2,kdphs_s=3
      integer :: kdphs(IPHMAX)=0
!
!/////////////////////////////////////////////////////////////////////
      contains
!
!===========================================================
      subroutine inputdata(ifli,ifll,ifle,KE2P,NPHS,cntlnam,E2P,ierror)
!===========================================================
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: ifli,ifll,ifle,KE2P,NPHS
      character(*),intent(in):: cntlnam
      logical,intent(out)    :: E2P
      integer,intent(out)    :: ierror
!
! --- [namelist]
!
      integer       :: Euler2ph=0,N_phase=2,I
      integer       :: kphs
      character(20) :: phase_type(IPHMAX)
      character(10),parameter ::
     &   GAS   ='GAS       ',
     &   LIQUID='LIQUID    ',
     &   SOLID ='SOLID     '
      character(10),parameter ::
     &    phstyp(3)=(/GAS,LIQUID,SOLID/)
!
      namelist /Eul2ph/ Euler2ph,N_phase,phase_type
!
! --- [local entities]
!
      integer :: ierr1=0,ios=0
!
      ierror=0
      Euler2ph=0
      phase_type(:)=' '
      rewind ifli
      read(ifli,Eul2ph,iostat=ios)
      if( ios.lt.0 ) return
      if( ios.gt.0 ) call FFRABORT(1,'module_Euler2ph')
!
      if(.not.(Euler2ph==0.or.Euler2ph==1.or.Euler2ph==2)) then
        write(ifle,'(1x,a)') 
     &   'ERR : Euler2ph must be 0/1/2 at &Eul2ph'

        write(ifle,'(1x,a)') 
     &  'MSG: Euler2ph=0:  No useing Euler two-phase model '

        write(ifle,'(1x,a)') 
     &   'MSG: Euler2ph=1:  Euler two-phase/one-pressure model '

        write(ifle,'(1x,a)') 
     &   'MSG: Euler2ph=2:  Euler two-phase/two-pressure model '

        write(ifle,'(1x,2a)') ' ### Reset [Euler2ph] in ',cntlnam
        
        call FFRABORT(1,'module_Euler2ph')
      endif
!
      if(N_phase<2.and.Euler2ph>0) then
        write(ifle,'(1X,a)') 
     & ' ### ERR : N_phase > 2 for Euler Multi-phase flow'
        write(ifle,*) ' ### Reset [N_phase] in ',cntlnam
        call FFRABORT(1,'Module_Euler2ph')
      endif
!
      if((KE2P==0.and.Euler2ph>0).or.(KE2P==1.and.Euler2ph==0)) then
        write(ifle,*) ' ### ERR : Euler2ph dimension error'
        write(ifle,*) ' ### MSG : Re-run prefflow'
        call FFRABORT(1,'module_Euler2ph/inputdata')
      endif
!
      if(Euler2ph>0) then
        if((NPHS+1)/=N_phase) then
          write(ifle,'(1x,a)') 'ERR : Euler2ph dimension error'
          write(ifle,'(1x,a)') 'MSG : [N_phase] have been changed'
          write(ifle,'(1x,a)') 'MSG : Re-run prefflow'
          call FFRABORT(1,'module_Euler2ph/inputdata')
        else
          N_PHS=NPHS
        endif
!
        do I=1,N_phase
        kphs=0
        call nml_listno(3,phstyp,phase_type(I),kphs)
        call nml_chkchr0(ifle,'phase_type',phase_type(I),kphs,ierr1)
        if(ierr1.ne.0) call FFRABORT(1,'ERR: ')
        phase_type(I)=phstyp(kphs)
        if(phase_type(I)==GAS) then
          kdphs(I)=kdphs_g
        elseif(phase_type(I)==LIQUID) then
          kdphs(I)=kdphs_l
        elseif(phase_type(I)==SOLID) then
          kdphs(I)=kdphs_s
        endif
        enddo
!
      endif
!
      ieul2ph=Euler2ph
!
      if(ieul2ph>0) then
        E2P=.true.
      else
        E2P=.false.
      endif
!
      return
!
      end subroutine inputdata
!
      end module module_Euler2ph
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_vof
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     key-word for VOF: ical_vof
!
      implicit none
      character(21),parameter,private :: modnam='(module_vof)'
!
      real(8),save :: evap_hh=0.d0,bldf_vof=0.8
      real(8),save :: coefcsf   
      integer,save :: iniflg=0
      integer,save :: ical_vof=0,change_phase=0
      integer,save :: LG=1,LS=2,intervof=1
      integer,parameter :: up1st=1,up2nd=2,cnt2nd=3,cicsam=4,hric=5
      integer,parameter :: grvty=1,byncy=2
      integer,parameter :: den1=1,den2=2
      integer,save :: cnvvof,modul,donr
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!=======================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,
     &  lvof,ivof,ierror)
!=======================================================
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle
      character(*),intent(in) :: cntlnam
      logical,intent(out) :: lvof
      integer,intent(out) :: ivof
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      character(LEN=3),parameter :: 
     &      cnvlst(5)=(/'1st','2nd','c2d','cic','hri'/)
      character(LEN=1),parameter :: modlst(3)=(/'g','b','n'/)
      integer :: ierr1=0,ios=0
!
! --- [namelist]
!
      integer       :: calvof=0
      integer       :: interface=1,donor=1
      character(LEN=10) :: cnv_vof,module
      real(8)  :: bldfct=0.8D0,evap_h
      
      namelist /VOF/ calvof,interface,evap_h,change_phase,cnv_vof,
     &               module,donor,bldfct,
     &               coefcsf,iniflg
!
      ierror=0
      calvof=0
      cnv_vof=' '
      cnvvof=up1st
      modul=grvty
      donr=den1
      
!
      rewind ifli
      read(ifli,VOF,iostat=ios)
      if( ios.lt.0 ) return
      if( ios.gt.0 ) call FFRABORT(1,'module_VOF')
!
      if(.not.(calvof.eq.0.or.calvof.eq.1)) then
        write(ifle,*) ' ### ERR : [calvof] must be 0 or 1 at &VOF'
        write(ifle,*) ' ### Reset [calvof] of &VOF in ',cntlnam
        write(ifll,*) 
     &     ' ### MSG : Liquid muse be defined in 1st phase'
        write(ifll,*) 
     &     ' ### MSG : VOF is defined as Liquid-VOF'
        call FFRABORT(1,'module_VOF')
      elseif(calvof.eq.1) then
        write(ifll,*) 
     &     ' ### MSG : Liquid muse be defined in 1st phase'
        write(ifll,*) 
     &     ' ### MSG : VOF is defined as Liquid-VOF'
        if(.not.(interface.eq.1.or.interface.eq.2)) then
          write(ifle,*) ' ### ERR : [interface] must be 1 or 2 !'
          write(ifle,*) ' ### MSG : interface=1 : liquid-gas interface'
        write(ifle,*) ' ### MSG : interface=2 : liquid-solid interface'
          call FFRABORT(1,'module_VOF')
        endif
        if(.not.(change_phase.eq.0.or.change_phase.eq.1)) then
          write(ifle,*) 
     &   ' ### ERR : [change_phase] must be 0 or 1 at &VOF'
          call FFRABORT(1,'module_VOF')
        endif
        if(.not.(donor.eq.den1.or.donor.eq.den2)) then
          write(ifle,*) 
     &  ' ### ERR : [donor]  must be 1 or 2 at &VOF'
        write(ifle,*) 
     & ' ### MSG : donor=1 => 1st phase is donor, 2nd phase is aceptor'
        write(ifle,*) 
     & ' ### MSG : donor=2 => 2nd phase is donor, 1st phase is aceptor'
        call FFRABORT(1,'module_VOF')
        endif
        call nml_listno(5,cnvlst,cnv_vof,cnvvof)
        call nml_chkchr0(ifle,'cnv_vof',cnv_vof,cnvvof,ierr1)
        if(ierr1.ne.0) call FFRABORT(1,'&VOF')
        call nml_listno(3,modlst,module,modul)
        call nml_chkchr0(ifle,'module',module,modul,ierr1)
        if(ierr1.ne.0) call FFRABORT(1,'&VOF')
!
        if(bldfct.gt.1.d0.or.bldfct.lt.0.d0) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'Blend Factor (bldfct) must between 1 and 0'
          call FFRABORT(1,'&VOF/bldfct')
        else
          bldf_vof=bldfct
        endif
      endif
!
      ical_vof=calvof
      intervof=interface
      donr=donor
!
      if(ical_vof.eq.1) then
        lvof=.true.
        ivof=1
      else
        lvof=.false.
        ivof=0
      endif
      evap_hh=evap_h
!
      return
!
      end subroutine inputdata
!
      end module module_VOF
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_FUEL
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
      character(21),parameter,private :: modnam='(module_FUEL)'
!
      real(8), save :: OPEN_V,osmo_Nd,Tau_diff
      integer, save :: No_Mem=0,No_AGDL=0,No_CGDL=0
      integer, save :: osmoticFC,col_no,vap_no
      real(8), save :: A_cat_th,C_cat_th
      real(8),parameter :: SML=1.d-15
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!=======================================================
      subroutine inputdata
     &  (ifli,ifll,ifle,ncompall,cntlnam,lFC,vapor,cool,ierror)
!=======================================================
!
      implicit none 
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle,lFC,ncompall,vapor,cool
      character(*),intent(in) :: cntlnam
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer :: ierr1=0,ios=0
!
! --- [namelist]
!
      integer       :: Mem_M_no,A_GDL_no,C_GDL_no,
     &                 electro_osmotic
      real*8        :: V_open,Nd_osmotic,Brug_tau
      real*8        :: A_catalyst_thick,C_catalyst_thick
      namelist /FUEL/ Mem_M_no,A_GDL_no,C_GDL_no,V_open, 
     &                electro_osmotic,Nd_osmotic,
     &                Brug_tau,
     &                A_catalyst_thick,C_catalyst_thick
!
      ierror=0
!
      if(lFC==0) return
      rewind ifli
      Mem_M_no=0
      A_GDL_no=0
      C_GDL_no=0
      V_open=0.d0
!      water_cool=7
      electro_osmotic=0
      A_catalyst_thick=0.d0
      C_catalyst_thick=0.d0
      read(ifli,FUEL,iostat=ios)
!      if( ios.lt.0 ) return
      if( ios.gt.0 ) call FFRABORT(1,'module_FUEL')
      if(Mem_M_no==0) then
        call FFRABORT(1,'ERR: NOT defining &FUEL/Mem_M_no')
      endif
      if(A_GDL_no==0) then
        call FFRABORT(1,'ERR: NOT defining &FUEL/A_GDL_no')
      endif
      if(C_GDL_no==0) then
        call FFRABORT(1,'ERR: NOT defining &FUEL/C_GDL_no')
      endif
      if(V_open==0.d0) then
        call FFRABORT(1,'ERR: NOT defining &FUEL/V_open')
      endif
      No_Mem=Mem_M_no
      No_AGDL=A_GDL_no
      No_CGDL=C_GDL_no
      OPEN_V=V_open
      if(cool>ncompall.or.cool<1) then
        call FFRABORT(1,'ERR: [PEFC-1]')
      else
        col_no=cool
      endif
      if(vapor>ncompall.or.vapor<1) then
        call FFRABORT(1,'ERR: [PEFC-1]')
      else
        vap_no=vapor
      endif
      osmoticFC=electro_osmotic
      osmo_Nd=max(0.d0,Nd_osmotic)
      Tau_diff=max(0.d0,Brug_tau)
      if(A_catalyst_thick<SML) then
        call FFRABORT(1,'ERR:NOT defined Anode catalyst thickness')
      else
        A_cat_th=A_catalyst_thick
      endif
      if(C_catalyst_thick<SML) then
        call FFRABORT(1,'ERR:NOT defined Cathode catalyst thickness')
      else
        C_cat_th=C_catalyst_thick
      endif
!
      return 
!
      end subroutine inputdata 
!
      end module module_FUEL

!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_scalar 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!
      character(21),parameter,private :: modnam='(module_scalar)'
!
      character(20),save,allocatable :: sclname(:)
      integer,save,allocatable :: scl(:)
!
      integer,save :: icalke(2),icalrsm(7),icalaph(2),ical_s
      integer,save :: ike(2),irsm(7),ivof,ides,icavi,ivold,ike_c,ikles
      logical,save :: rns_scl=.false.
      integer,parameter :: IPHMAX=10
      integer,save,allocatable :: iaph(:)
      integer,save :: icalvof,ical_FC
!
      character(20),save,allocatable :: potname(:)
      integer,save :: ip_pefc(2) 
      logical,save :: pot_scl=.false.,flag_flm_ys=.false.
! ---
      logical :: calgeq
      logical :: calxi
!
      logical :: caltemp
      logical :: calh

      integer :: igeq

      integer,save :: ixi
      integer,save :: ixi_C     !: Reaction degress
      integer,save :: ixi_2     !: Root meat square: rms(ixi)
      integer,save :: ixi_X     !: Scalar disspation rate
      integer,save :: ixi_WC
      integer,save :: idifflm=0
      integer,save :: comb_mdl=0
      integer,save :: frmt=1
      integer,parameter,public :: ORG=0,consF=1,Blog=2,BlogM=3,
     &                            FltRans=4
!
      logical,save :: ical_prp=.true.,ical_flmlt=.false.
      integer,save :: igT_iter_flm=0
!
      integer :: itemp
      integer :: ihflmlt
      integer,save :: istmdl,J
!
      real*8,allocatable :: fd_rho(:)
      real*8,allocatable :: fd_tmp(:)
      real*8,allocatable :: fd_lbv1(:)
      real*8,allocatable :: fd_lbv2(:)
      real*8 :: lbv_XXi(3),xi
      integer :: nrho
      integer :: ntmp
      integer :: nlbv1,nlbv2
!
      real*8, save :: XXi0,
     &                Xfuel_rho,
     &                Xfuel_temp,
     &                XOxidant_rho,
     &                XOxidant_temp,
     &                XXi_max
!
      integer,save :: ical_cavi=0,iterCVAI,iterPOTEN
      integer,save :: Axis_DIR=1,FSTR_F=0
      integer,save :: CAVI_B=-1
      real*8 ,save :: PWR=1.d0
!
      real*8,save  :: Pres_c=1944.61d6,P_str_v,Tref=288.15D0,
     &                T0_ref,KL,Rv=461.6D0, !(J/kg/K)
     &                a_cavi,b_cavi,c_cavi,T_satu,CeA_PP,CeA_MM
!
      character(LEN=80),save :: flmletfile='flamelet.data'
!
      integer,save :: n1,n2,n3,nvar
      integer,save :: nm,J_RHO
      real*8,save,allocatable        :: x1(:),x2(:),x3(:),xx(:)
      real*8,allocatable,save        :: table(:,:,:,:)
      real*8,allocatable,save        :: ORGtable(:,:),ORGspd(:,:)
      character(16),save,allocatable :: tags(:)
      character(16),save             :: tagss(2)
      integer,save :: NVAR_ORG=0,NVAR_COMP=3,iz,zmin,zmax
      character(80) ::  cdum
      real*8,save :: x1min,x1max
      real*8,save :: x2min,x2max
      real*8,save :: x3min,x3max
!
      integer,parameter :: mm=30
      real*8  :: POLYC(mm+1,2),dum1,dum2
      real*8,allocatable,save :: POLYC_YS(:,:) !(mm+1,ncomp)
      integer :: nn1,nn2,IFLAG
!
      integer,parameter :: mmCD=7,nnCD=7
      integer,save :: mCD=mmCD,nCD=nnCD
      real*8  :: POLYCDP(mmCD)
      real*8  :: re(nnCD),cdp(nnCD)
!//////////////////////////////////////////////////////////////////////
      contains
!======================================================================
      subroutine inputdata
     & (ifli,ifll,ifle,cntlnam,my_rank,
     &  lcomp,lrans,E2P,kemdl,lowke,l2SKE,lRNG,lCHEN,lSDES,lKLES,
     &  RSMmdl,RANSMDL,lvof,ncomp,nrans,NPHS,NPOTN,lcavi,lFC,lpotn,
     &  alpl,ikel,jvof,lflm,FLMLT,iCAVIT,iPEFC,ierror)
!======================================================================
      use module_model,only : u_func
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: ifli,ifll,ifle,my_rank
      character(*),intent(in) :: cntlnam
      integer,intent(in)    :: nrans,NPHS,NPOTN,ncomp
      logical,intent(in)    :: E2P,lcomp,kemdl,RSMmdl,lvof,lowke,
     &                         lRNG,lCHEN,lSDES,RANSMDL,l2SKE,lKLES
      logical,intent(inout) :: lcavi,lpotn
      integer,intent(inout) :: FLMLT,lFC,iCAVIT,iPEFC,lflm
      logical,intent(inout) :: lrans
      character(20)         :: Geq,Xieq,Tempeq,heq
      integer,intent(out)   :: ierror,alpl(2),ikel(2),jvof
!
! --- [namelist]
!
      integer,parameter :: lc=10
!
      character(lc),parameter ::
     &   ORG_les =  'ORG_les   ',
     &   les_const =  'les_const ',
     &   C_les=  'C_les     ',
     &   H_les=  'H_les     ',
     &   rans_const=  'rans_const'
!
      character(lc),parameter ::
     &    flmknd(5)=(
     &   /ORG_les,les_const,C_les,H_les,rans_const/)
!
      integer,parameter :: mrans=50
      integer,parameter :: mpotn=10
      integer,parameter :: mflame=50
      character(20) :: scalar_name(mrans)
      character(20) :: Potential_name(mpotn)
      character(80) :: flamelet_file,Lam_spd_file
      real*8,parameter :: undef=-huge(1.)
!      
! --- Partial-pre-mixed:            -------------------------------
!
      real*8  :: Xi0,              ! Xi (constant value)
     &           fuel_rho,         ! Fuel density in unburn
     &           fuel_temp,        ! Fuel temperature in unburn
     &           Oxidant_rho,      ! Oxidant density in unburn
     &           Oxidant_temp,     ! Oxidant temperature in unburn
     &           Xi_max            ! Max Xi for Xi-Density profile
! -----------------------------------------------------------------
      namelist /scalar/ scalar_name!,Geq,Xieq,Tempeq,heq
!
      integer :: tur_vel=1,date_base=0,lokflm=0,igT_iter=0,
     &           DENformat=-1
!
      
      character(lc) :: flm_model 
!     
      namelist /flamelet/ 
     &           Geq,
     &           Xieq,
     &           Tempeq,
     &           heq,
     &           tur_vel,
     &           flamelet_file,Lam_spd_file,
     &           Xi0,              ! Xi (constant value)
     &           fuel_rho,         ! Fuel density in unburn
     &           fuel_temp,        ! Fuel temperature in unburn
     &           Oxidant_rho,      ! Oxidant density in unburn
     &           Oxidant_temp,     ! Oxidant temperature in unburn
     &           Xi_max,           ! Max Xi for Xi-Density profile
     &           date_base,        ! =0 CHEMKIN datebase 
     &                             ! =1 FLAMEMASTER for steady model
     &                             ! =2 FLAMEMASTER for unsteady model
     &           flm_model,        ! 'ORG_les'
                                   ! 'les_const'
                                   ! 'C_les'
                                   ! 'H_les'
                                   ! 'rans_const'
     &           igT_iter,DENformat
!
      integer :: Cavi=0,Iter_cavi,Iter_poten,ICOM
      integer :: Axis=1,Boundary_NO=-1,FSTR_OUT=0
      real*8  :: Power=1.d0
      real*8  :: rho(mflame),tmp(mflame),lbv1(mflame),lbv2(mflame)
      real*8  :: lbv_Xi(3),dum1
      integer :: idum1
!
      real*8  :: T_ref=288.15,Pc=1944.61d6,T_sat,CeA_P,CeA_M
!
      namelist /Cavitation/ Cavi,T_ref,Pc,T_sat,
     &                      Iter_cavi,CeA_P,CeA_M,
     &                      Axis,Boundary_NO,FSTR_OUT,Power
!
      namelist /Potential/ Potential_name,Iter_poten
!
      namelist /flamelet_func/ nrho,rho,ntmp,tmp,nlbv1,nlbv2,
     &         lbv1,lbv2,lbv_Xi
! --- [local entities]
!
      integer :: i_user,nransno=0,no=0,ios=0,ierr=0,nmodel=mrans
      integer :: NPOTNno=0
      integer :: icalgeq,icalxi,icaltemp,icalh,icalmdl
      integer :: i
      character(3),parameter :: flmltlst(2)=(/'no ','yes'/) 
      character(3) :: text
      logical :: fexist=.false.
!
      ierror=0
      alpl(:)=0
      scalar_name(:)=' '
      calgeq  = .false.
      calxi   = .false.
      caltemp = .false.
      calh    = .false.
      Geq     = 'no '
      Xieq    = 'no '
      Tempeq  = 'no '
      heq     = 'no '
      igeq    = 1
      ixi     = 1
      itemp   = 1
      ihflmlt = 1
      CeA_PP  = 10.d0
      CeA_MM  = 1.d0
!
      Xi0=-1.d0
      fuel_rho=1.d0
      fuel_temp=298.15d0
      Oxidant_rho=1.205d0
      Oxidant_temp=298.15d0
      Xi_max=-1.d0
      flamelet_file=''
      Lam_spd_file=''
      flm_model='ORG_les'
      igT_iter=-1
      DENformat=-1
!BlogB=0
      date_base=0
! --- 
      nrho=0
      ntmp=0
      nlbv1=0
      nlbv2=0
      rho=undef
      tmp=undef
      lbv1=undef
      lbv2=undef
      lbv_Xi=-1.d0
!
      
!
      rewind ifli
      read(ifli,scalar,iostat=ios)
      if(ios.gt.0) call FFRABORT(1,'reading error at &scalar')
!
      rewind ifli
      read(ifli,flamelet,iostat=ios)
      if(ios.gt.0) call FFRABORT(1,'reading error at &flamelet')
!
      rewind ifli
      read(ifli,Cavitation,iostat=ios)
      if(ios.gt.0) call FFRABORT(1,'reading error at &Cavitation')
! 
!
      call nml_listno(2,flmltlst,Geq,icalgeq)
      call nml_chkchr0(ifle,'Geq',Geq,icalgeq,ierror)
      if(ierror.ne.0) then
        call FFRABORT(1,'reading error at &flamelet about Geq')
      endif
      icalgeq=icalgeq-1
!
! --- Set "icalxi" value
!
      call nml_listno(2,flmltlst,Xieq,icalxi)
      call nml_chkchr0(ifle,'Xieq',Xieq,icalxi,ierror)
      if(ierror.ne.0) then
        call FFRABORT(1,'reading error at &flamelet about Xieq')
      endif
      icalxi=icalxi-1
!
! --- Set "icaltemp" value
!
      call nml_listno(2,flmltlst,Tempeq,icaltemp)
      call nml_chkchr0(ifle,'Tempeq',Tempeq,icaltemp,ierror)
      if(ierror.ne.0) then
        call FFRABORT(1,'reading error at &flamelet about Tempeq')
      endif
      icaltemp=icaltemp-1
!-----------------------
! --- Set "icalh" value
!-----------------------
      call nml_listno(2,flmltlst,heq,icalh)
      call nml_chkchr0(ifle,'heq',heq,icalh,ierror)
      if(ierror.ne.0) then
        call FFRABORT(1,'reading error at &flamelet about heq')
      endif
      icalh=icalh-1
!
      icalmdl=0
      call nml_listno(5,flmknd,flm_model,icalmdl)
      call nml_chkchr0(ifle,'flm_model',flm_model,icalmdl,ierror)
      if(ierror.ne.0) then
        call FFRABORT(1,'reading error at &flamelet/flm_model')
      endif
      icalmdl=icalmdl-1
      idifflm=icalmdl
!----------------------------------------------
! --- set flamelet logical flags
!----------------------------------------------
      if(icalgeq .eq.1) calgeq  = .true.
      if(icalxi  .eq.1) calxi   = .true.
      if(icaltemp.eq.1) caltemp = .true.
      if(icalh   .eq.1) calh    = .true.
!
      if(calgeq.or.calxi.or.caltemp.or.calh) then
        if(igT_iter<0) then
!          write(ifle,*) 
!     &  'ERR: set &flamelet/igT_iter for flamelet model ignation step'
!          call FFRABORT(1,'MSG: set &flamelet/igT_iter ')
        else
           igT_iter_flm=igT_iter
        endif
      else
        igT_iter_flm=-1
      endif
!
      if(calxi.and.flm_model/='ORG_les') then
        if(DENformat==-1) then
          write(ifle,'(1X,a)') 
     &    'ERR: set &flamelet/DENformat for flamelet model'
          write(ifle,'(1X,a)') 
     &    'MSG: DENformat=1 for BABA format'
          write(ifle,'(1X,a)') 
     &    'MSG: DENformat=2 for WATANABE format'
          call FFRABORT(1,'ERR: NOT defined [DENformat]')
        elseif(DENformat/=1.and.DENformat/=2) then
          write(ifle,'(1X,a)') 
     &    'ERR: set &flamelet/DENformat=1 or =2 for flamelet model'
          call FFRABORT(1,'ERR: NOT defined [DENformat]')
        else
          frmt=DENformat
        endif
      endif
!      
      if(calxi) then
        idifflm=icalmdl
        if(idifflm==ORG) then  !=ORG
          write(ifll,'(1x,a)') 
     &   'MSG:ORG LES flamelet model is used for diffusion flame'
          ical_prp=.true.
          if(date_base/=0) then
            call FFRABORT(1,'ERR: date_base=0 for "ORG_les"')
          endif
        elseif(idifflm==consF) then  !=Blog
          write(ifll,'(1x,a)') 
     &   'Constant LES flamelet model is used for diffusion flame'
          ical_prp=.false.
          if(date_base/=1) then
            call FFRABORT(1,'ERR: date_base=1 for "les_const"')
          endif
        elseif(idifflm==Blog) then  !=Blog
          write(ifll,'(1x,a)') 
     &   'Blog LES flamelet model is used for diffusion flame'
          ical_prp=.false.
          if(date_base/=2) then
            call FFRABORT(1,'ERR: date_base=2 for "C_les"')
          endif
        elseif(idifflm==BlogM) then  !=BlogM
          write(ifll,'(1x,a)') 
     &   'Modified Blog LES flamelet model is used for diffusion flame'
          ical_prp=.false.
          if(date_base/=2) then
            call FFRABORT(1,'ERR: date_base=2 for "H_les"')
          endif
        elseif(idifflm==FltRans) then  !=FltRans
          write(ifll,'(1x,a)') 
     &   'Constant RANS flamelet model is used for diffusion flame'
          ical_prp=.false.
          if(date_base/=1) then
            call FFRABORT(1,'ERR: date_base=1 for "rans_const"')
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
          call FFRABORT(1,'ERR:  &flamelet/flm_model')
        endif
!
        if(idifflm/=ORG) then
          if(date_base==0) then
             write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base=0'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=1 for [chemtable_CH4_steady]'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=2 for [chemtable_CH4_unsteady]'
             call FFRABORT(1,'ERR: date_base=1 or date_base=2')
          elseif(date_base==1) then
            write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base=1'
            write(ifle,'(1x,a)') 
     &      'MSG: date_base=1 for [chemtable_CH4_steady]'
          elseif(date_base==2) then
            write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base=2'
            write(ifle,'(1x,a)') 
     &      'MSG: date_base=2 for [chemtable_CH4_unsteady]'
          else
            write(ifle,'(1x,a)') 
     &      'MSG: Flamelet date-base model: date_base'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=1 for [chemtable_CH4_steady]'
             write(ifle,'(1x,a)') 
     &      'MSG: date_base=2 for [chemtable_CH4_unsteady]'
            call FFRABORT(1,'ERR: date_base=1 or date_base=2')
          endif
        else  !if(idifflm==ORG) then
          if(date_base/=0) then
          date_base=0
          call FFRABORT(1,'ERR: date_base=0 for flm_model="ORG_les"')
          endif
        endif
      endif
!
      if(calgeq.and.calxi.and.idifflm/=ORG) then
        write(ifle,'(1x,2a)') 
     & 'MSG: ONLY flm_model="ORG_les" ',
     &   'support Partial pre-mixed flame model'
        call FFRABORT
     &    (1,'ERR: set [&flamelet/flm_model="ORG_les"]')
      elseif(calgeq.and.idifflm/=ORG) then
        write(ifle,'(1x,a)') 
     &   'MSG: ONLY flm_model="ORG_les" support pre-mixed flame model'
        call FFRABORT
     &    (1,'ERR: set [&flamelet/flm_model]')
      endif
!
      ical_cavi=Cavi
      if(Cavi==1) then
!        if(.not.lcomp) 
!     &  call FFRABORT(1,'ERR:Cavitation model must use [comp]')
        Pres_c=Pc
        Tref=T_ref
        a_cavi=3.353D-9*Pres_c
        b_cavi=-1.723D-6*Pres_c
        c_cavi=1.222D-3*Pres_c
        KL=2.d0*a_cavi*Tref+b_cavi
        T0_ref=(a_cavi*Tref*Tref+b_cavi*Tref+c_cavi)/KL-Tref
!        mu_ll=mu_l
!        mu_vv=mu_v
        T_satu=T_sat
        iterCVAI=max(0,Iter_cavi)
        iCAVIT=1
        CeA_PP=CeA_P
        CeA_MM=CeA_M
!
        CAVI_B=Boundary_NO 
        Axis_DIR=Axis 
        FSTR_F=FSTR_OUT
        PWR=Power
!
        if(Power>1.d0.or.Power<0.d0) then
          call FFRABORT(1,'ERR: Cavitation [Power] MUST BE [0~1]')
        endif
!
        if(.NOT.(FSTR_F==1.or.FSTR_F==0)) then
          call  FFRABORT
     &     (1,'ERR: FrontStar Force flag [FSTR_F] MUST BE [0 or 1]')
        endif
!
        if(.NOT.(Axis_DIR==1.or.Axis_DIR==2.or.Axis_DIR==3)) then
          call  FFRABORT
     &  (1,'ERR: Cavitation Dir MUST BE x,y or z By Axis_DIR=1,2 or 3')
        endif
!
      elseif(.not.(Cavi/=0.or.Cavi/=1)) then
        call FFRABORT(1,'ERR: [Cavi] must be 0 or 1')
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
      
!
      if(lFC>0)then
        no=no+1
      endif
!
      if(Cavi==1) then
        no=no+2
        lcavi=.true.
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
          no=no+2   ! ixi_2,ixi_
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
      if(ierr.ne.0) call FFRABORT(1,'allocate error for [sclname]')
!-------------------
! --- initialization
!-------------------
      sclname(:)=' '
      scl(:)=0

      icalke(:)=0
      ike(:)=1
      icalrsm(:)=0
      irsm(:)=1
      icalaph(:)=0   
      icalvof=0
      ivof=1
      ides=1
      icavi=1
      ical_s=1
      ical_FC=0
!
      igeq    = 1
      ixi     = 1
      itemp   = 1
      ihflmlt = 1
!
      ixi_C = -1
      ixi_2 = -1
      ixi_X = -1
      ixi_WC= -1
!
      ike_c=1
      ikles=1
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
        icalke(1)=1
        ike(1)=no
        ikel(1)=no

!
        no=no+1
        sclname(no)='RANS_eps'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalke(2)=1
        ike(2)=no
        ikel(2)=no
!
        if(l2SKE) then
          Re(1)=03.9d0;CDp(1)=1.250d0
          Re(2)=04.d0;CDp(2)=1.230d0
          Re(3)=10.d0;CDp(3)=0.714d0
          Re(4)=20.d0;CDp(4)=0.500d0
          Re(5)=30.d0;CDp(5)=0.450d0
          Re(6)=40.d0;CDp(6)=0.414d0
          Re(7)=89.d0;CDp(7)=0.286d0
          call STP(nCD,Re,cdp,mCD,POLYCDP(:),IFLAG)
          write(ifll,'(1X,a)') 'MSG: Poly for 2SKE model'
          do j=1,nCD
          xi=Re(j)
          call fitVALUE(xi,POLYCDP(:),mCD+1,dum1)
          write(*,'(1x,I3,3(1X,a,E16.5))') 
     &       j,'Re=',Re(j),'CDp=',cdp(j),'Cdp_p=',dum1
          enddo
        endif
      endif
!
      if(lSDES) then
        no=no+1
        sclname(no)='DES_Nu'
        nmodel=nmodel+1
        scl(no)=nmodel
        ides=no
      endif
!
      if(RSMmdl) then
        no=no+1
        sclname(no)='eps'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalrsm(1)=1
        irsm(1)=no
!
        no=no+1
        sclname(no)='uu'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalrsm(2)=1
        irsm(2)=no
!
        no=no+1
        sclname(no)='vv'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalrsm(3)=1
        irsm(3)=no
!
        no=no+1
        sclname(no)='ww'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalrsm(4)=1
        irsm(4)=no
!
        no=no+1
        sclname(no)='uv'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalrsm(5)=1
        irsm(5)=no
!
        no=no+1
        sclname(no)='uw'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalrsm(6)=1
        irsm(6)=no
!
        no=no+1
        sclname(no)='vw'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalrsm(7)=1
        irsm(7)=no
      endif
!
      allocate(iaph(IPHMAX))
      if(IPHMAX<NPHS+1) 
     &   call FFRABORT(1,'ERR:MultiPhase Max Phase < 10 ')
      if(E2P) then
        do I=1,NPHS+1
        write(text,'(i3)') I
        text=adjustl(text)
        no=no+1
!        
        sclname(no)='alph'//TRIM(text)
        nmodel=nmodel+1
        scl(no)=nmodel
        icalaph(I)=I
        iaph(I)=no
        alpl(I)=no
        enddo
!
!        no=no+1
!        sclname(no)='alph2'
!        nmodel=nmodel+1
!        scl(no)=nmodel
!        icalaph(2)=1
!        iaph(2)=no
!        alpl(2)=no
      else
        
      endif
!
      if(lvof) then
        no=no+1
        sclname(no)='VOF'
        nmodel=nmodel+1
        scl(no)=nmodel
        icalvof=1
        ivof=no
        jvof=no
!
      endif
!
      if(lFC>0) then
        no=no+1
        sclname(no)='Saturation'
        nmodel=nmodel+1
        scl(no)=nmodel
        ical_s=no
        ical_FC=lFC
      endif
!
      if(ical_cavi==1) then
        no=no+1
        sclname(no)='Cavi_Y'
        nmodel=nmodel+1
        scl(no)=nmodel
        icavi=no
!
        no=no+1
        sclname(no)='VOID'
        nmodel=nmodel+1
        scl(no)=nmodel
        ivold=no
      endif
!
      if(lKLES) then
        no=no+1
        sclname(no)='LES_K'
        nmodel=nmodel+1
        scl(no)=nmodel
        ikles=no
      endif
!
      if(calgeq) then
        no=no+1
        sclname(no)='G'
        nmodel=nmodel+1
        scl(no)=nmodel
        igeq=no
        FLMLT=1
        lflm=1
      endif
!
      if(calxi) then
        no=no+1
        sclname(no)='Xi'
        nmodel=nmodel+1
        scl(no)=nmodel
        ixi=no
        lflm=1
!
        if(idifflm==ORG) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_2=no
!
          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_X=no
        elseif(idifflm==consF) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_2=no
!
          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_X=no
        elseif(idifflm==Blog) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_2=no

          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_X=no

          no=no+1
          sclname(no)='REAC_DEGRESS'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_C=no

          no=no+1
          sclname(no)='WDOT_C'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_WC=no
        elseif(idifflm==BlogM) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_2=no

          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_X=no

          no=no+1
          sclname(no)='REAC_DEGRESS'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_C=no

          no=no+1
          sclname(no)='WDOT_C'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_WC=no
        elseif(idifflm==FltRans) then
          no=no+1
          sclname(no)='RMS_Xi'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_2=no
!
          no=no+1
          sclname(no)='Xi_DISS_RATE'
          nmodel=nmodel+1
          scl(no)=nmodel
          ixi_X=no
        endif
      endif
!
      if(calxi) then
        if(.NOT.RANSMDL.and.idifflm==FltRans) then
           write(ifle,'(1x,a)') 
     &   'MSG: flm_model="rans_const": Constant RANS flamelet model'
          call FFRABORT
     &   (1,'ERR: RANS model MUST be set for flm_model="rans_const"')
        elseif(RANSMDL.and.(idifflm/=FltRans.and.idifflm/=ORG)) then
          write(ifle,'(1x,a)') 
     &   'MSG: BlogB=4: Constant RANS flamelet model'
          call FFRABORT
     &   (1,'ERR: RANS model ONLY support flm_model="rans_const"')
        endif
      endif
!
      if(caltemp) then
        no=no+1
        sclname(no)='Temp'
        nmodel=nmodel+1
        scl(no)=nmodel
        itemp=no
      endif
!
      if(calh) then
        no=no+1
        sclname(no)='h_flmlt'
        nmodel=nmodel+1
        scl(no)=nmodel
        ihflmlt=no
      endif
!---------------------
! --- 
!---------------------
      if(calxi.or.calgeq) then
        ical_flmlt=.true.
!tur_vel=1:Damk\"oler  2:Yakhot
        if(tur_vel==1.or.tur_vel==2) then
          write(ifll,'(1x,a)') 'MSG: tur_vel=1: Damk\"oler model'
          write(ifll,'(1x,a)') 'MSG: tur_vel=2: Yakhot model'
          istmdl=tur_vel
        else
          write(ifle,'(1x,a)') 'ERR: turbulent flame speed model'
          write(ifle,'(1x,a)') 'MSG: tur_vel=1: Damk\"oler model'
          write(ifle,'(1x,a)') 'MSG: tur_vel=2: Yakhot model'
          call FFRABORT(1,'MSG:tur_vel/=1 or 2')
        endif
!
        if(calxi.and.calgeq) then 
          write(ifle,'(1x,a)')  
     &         'MSG: Partial pre-mixed flame model'
          write(ifle,'(1x,a)')  
     &         'MSG: Set [fuel_rho] for unburn Fuel density'
          write(ifle,'(1x,a)')  
     &         'MSG: Set [fuel_temp] for unburn Fuel temperature'
          write(ifle,'(1x,a)')  
     &         'MSG: Set [Oxidant_rho] for unburn Oxidant density'
          write(ifle,'(1x,a)') 
     &         'MSG: Set [Oxidant_temp] for unburn Oxidant temperature'
          write(ifle,'(1x,a)') 
     &         'MSG: Set [Xi_max] for Max Xi of Xi-Density profile'
!
          
          if(Xi_max<0.d0) then 
            write(ifle,'(1x,a)') 
     & 'MSG: [Xi_max]<0 : not define Maximun Xi of Xi-Density profile' 
!
            write(ifle,'(1x,3a)') 
     &      'MSG: Unburn Value: Fuel density,',
     &      'Oxidant density, Fuel temperature,  Oxidant density ',
     &      'should be decided by polynomial expression'
            XXi_max=Xi_max
          else 
            write(ifle,'(1x,a,E10.5)')  
     &         'MSG: Unburn [fuel_rho]=',fuel_rho
            write(ifle,'(1x,a,E10.5)')  
     &         'MSG: Unburn [fuel_temp]=',fuel_temp
            write(ifle,'(1x,a,E10.5)')  
     &         'MSG: Unburn [Oxidant_rho]=',Oxidant_rho
            write(ifle,'(1x,a,E10.5)') 
     &         'MSG: Unburn [Oxidant_temp]=',Oxidant_temp
            write(ifle,'(1x,a,E10.5)') 
     &         'MSG: [Xi_max]=',Xi_max
            Xfuel_rho     =fuel_rho
            Xfuel_temp    =fuel_temp
            XOxidant_rho  =Oxidant_rho
	    XOxidant_temp =Oxidant_temp
            XXi_max       =Xi_max
          endif
!
        elseif(.not.calxi.and.calgeq) then 
          write(ifle,'(1x,a,E10.5)') 
     &         'MSG: Unburn [Xi0]=',Xi0 
          XXi0=Xi0
          if(Xi0<0.d0) then 
            write(ifle,'(1x,2a)') 
     &      'ERR: Set [Xi0] for fuel concentration value of ',
     &         'pre-mixed flame model' 
            call FFRABORT(1,'ERR: in pre-mixed flame model') 
          else
          endif
        endif
!
        if(idifflm==ORG) then
!----------------------------------------
! --- ORG flamelet LES combustion model
!----------------------------------------
          flmletfile=trim(flamelet_file)
!
          if(flmletfile=='') then
            call FFRABORT
     &     (1,'Set flamelet date base file &flamelet/flamelet_file')
          endif
!       
          inquire(file=flmletfile,exist=fexist)
!
          if (.not.fexist) then
            write(ifle,*) 
     &      'MSG: FLAMELET DATA FILE file [',
     &      trim(flmletfile),'] NOT exist'
            call FFRABORT(1,'Please check FLAMELET DATA FILE name')
          endif
!
          open(50,file=flmletfile,iostat=ios)
!
          if(ios.gt.0) then
! --- sa 
            write(ifle,'(2X,2a)') 
     &      'ERR: Reading flamelet file: ',trim(flmletfile)
            call FFRABORT(1,'Error on reading flamelet data')
          endif
! --- Read a number of Density(Inverse) Polynomial Coefficients
          rewind(50)
          read(50,flamelet_func,iostat=ios,ERR=130) 
          if(ios.gt.0) then 
            write(ifle,*) 'ERR: Reading File error:',trim(flmletfile)
            call FFRABORT(1,'reading error at &flamelet_func')
          endif
          if(nrho==0.and.ntmp==0) then
!
! -- CHEMKIN
!
            close(50)
            open(50,file=flmletfile,iostat=ios,
     &        form='FORMATTED') 
            rewind(50)
            read(50,*,ERR=120) idum1
            nvar=idum1
            NVAR_ORG=nvar
            allocate(tags(nvar))
            read(50,*) (tags(i),i=1,nvar) !tags
            if(my_rank==0) then
              do i=1,nvar
              write(ifll,'(1x,a,I4,4X,a)') 'MSG: ',I,trim(tags(I))
              enddo
            endif
            idum1=0
!            
            do 
               read(50,*,iostat=ios) dum1
            if(ios<0) exit
            idum1=idum1+1
            enddo
!
            iz=idum1
            write(ifll,'(1x,a,4X,I8)') 'MSG: Colume number=',idum1
            allocate(ORGtable(idum1,nvar))
            allocate(x1(idum1),xx(idum1))
            rewind(50)
            read(50,*,ERR=120) nvar
            read(50,*) (tags(i),i=1,nvar)
            do i=1,idum1
            read(50,*,ERR=200) (ORGtable(i,j),j=1,nvar)
            enddo
            close(50)
            if(my_rank==0) then
              write(cdum,*) "(15x,",nvar,"a15",")"
              write(ifll,cdum) (tags(i),i=1,nvar)
              write(cdum,*) "(4x,","I8,",nvar,"E15.6",")"
              do i=1,idum1
              write(ifll,cdum) i,(ORGtable(i,j),j=1,nvar)
              enddo
            endif
!
            x1(:)=ORGtable(:,1)
            xx(:)=ORGtable(:,3)
            call STP(idum1,x1,xx,mm,POLYC(:,1),IFLAG)
            do j=1,idum1
            xi=x1(j)
            call fitVALUE(xi,POLYC(:,1),mm+1,dum1)
            if(my_rank==0) then
              write(ifll,'(1x,I6,3(1X,a,E16.5))') 
     &         j,'Z=',x1(j),'rho=',xx(j),'rho_poly=',dum1
            endif
            enddo
            nrho=mm+1
            allocate(fd_rho(0:nrho-1),stat=ierr)
            if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_rho") 
            fd_rho=0.0d0
            do i=0,nrho-1
            fd_rho(i)=POLYC(nrho-i,1)
            enddo
!
            xx(:)=ORGtable(:,2)
            call STP(idum1,x1,xx,mm,POLYC(:,1),IFLAG)
            do j=1,idum1
            xi=x1(j)
            call fitVALUE(xi,POLYC(:,1),mm+1,dum1)
            if(my_rank==0) then
              write(*,'(1x,I6,3(1X,a,E16.5))') 
     &         j,'Z=',x1(j),'tmp=',xx(j),'tmp_poly=',dum1
            endif
            enddo
            ntmp=mm+1
            allocate(fd_tmp(0:ntmp-1),stat=ierr)
            if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_tmp") 
            fd_tmp=0.0d0
            do i=0,ntmp-1
            fd_tmp(i)=POLYC(ntmp-i,1)
            enddo
!---------------------------
! --- species vailazation 
!---------------------------POLYC_YS(:,:)
            if(nvar-NVAR_COMP/=ncomp) then
              write(ifle,'(1X,a)') 'set &species and ncomp in fflow.ctl'
              call FFRABORT(1,'ERR: nvar-3/=ncomp')
            endif

            if(nvar-NVAR_COMP==ncomp) then
              flag_flm_ys=.true.
              allocate(POLYC_YS(mm+1,ncomp),stat=ierr)
              if(ierr/=0) call FFRABORT(1,"ERR: allocating POLYC_YS") 
              do ICOM=1,NCOMP
              xx(:)=ORGtable(:,NVAR_COMP+icom)
              call STP(idum1,x1,xx,mm,POLYC_YS(:,ICOM),IFLAG)

!              do j=1,idum1
!              xi=x1(j)
!              call fitVALUE(xi,POLYC_YS(:,ICOM),mm+1,dum1)
!              write(*,'(1x,I6,3(1X,a,E16.5))') 
!     &         j,'Z=',x1(j),'tmp=',xx(j),'tmp_poly=',dum1
!              enddo
!              stop 4433
              enddo
            endif
!
! ---
!
            if(calgeq) then
              nlbv1=mm+1
              nlbv2=mm+1
              ierr=0
              allocate(fd_lbv1(0:nlbv1-1),stat=ierr)
              if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_lbv1-1") 
              allocate(fd_lbv2(0:nlbv2-1),stat=ierr)
              if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_lbv2-1") 
              fd_lbv1=0.0d0
              fd_lbv2=0.0d0


              flmletfile=trim(Lam_spd_file)
              if(flmletfile=='') then
                call FFRABORT
     &          (1,'Set flamelet date base file &flamelet/Lam_spd_file')
              endif
              inquire(file=flmletfile,exist=fexist)
              if (.not.fexist) then
                write(ifle,*) 
     &          'MSG: FLAMELET DATA FILE file [',
     &          trim(flmletfile),'] NOT exist'
                call FFRABORT(1,'Please check FLAMELET DATA FILE name')
              endif
              open(50,file=flmletfile,iostat=ios,
     &        form='FORMATTED') 
              rewind(50)
              read(50,*) (tagss(i),i=1,2)
!            
              idum1=0
              do 
              read(50,*,iostat=ios) dum1
              if(ios<0) exit
              idum1=idum1+1
              enddo
              allocate(ORGspd(idum1,2))
              rewind(50)
              read(50,*) (tagss(i),i=1,2)
              do i=1,idum1
              read(50,*,ERR=200) (ORGspd(i,j),j=1,2)
              enddo
              close(50)
              nvar=2
              write(cdum,*) "(15x,",nvar,"a15",")"
              write(ifll,cdum) (tagss(i),i=1,2)
              write(cdum,*) "(4x,","I8,",nvar,"E15.6",")"
              do i=1,idum1
              write(ifll,cdum) i,(ORGspd(i,j),j=1,2)
              enddo
!
              dum2=-1.d0
              J=0
              do i=1,idum1
                dum1=ORGspd(i,2)
                if(dum1>dum2) then
                  dum2=dum1
                  J=I
                endif
              enddo
              lbv_XXi(1)=ORGspd(1,1)
              lbv_XXi(2)=ORGspd(J,1)
              lbv_XXi(3)=ORGspd(idum1,1)
              write(ifll,'(1X,a,3F15.4)')'2 POLYNOMIAL RANG:',lbv_XXi(:)
!
              deallocate(x1,xx)
              allocate(x1(idum1),xx(idum1))
              x1=0.d0
              xx=0.d0
              do I=1,J
              x1(I)=ORGspd(I,1)
              xx(I)=ORGspd(I,2)
              enddo

              POLYC(:,1)=0.d0
              call STP(J,x1,xx,mm,POLYC(:,1),IFLAG)
              do j=1,J
              xi=x1(j)
              call fitVALUE(xi,POLYC(:,1),mm+1,dum1)
              write(*,'(1x,I6,3(1X,a,E16.5))') 
     &           j,'Z=',x1(j),'spd1=',xx(j),'spd_poly1=',dum1
              enddo
              do i=0,nlbv1-1
              fd_lbv1(i)=POLYC(nlbv1-i,1)
              enddo
              
              
              x1=0.d0
              xx=0.d0
              do I=J+1,idum1
              x1(I)=ORGspd(I,1)
              xx(I)=ORGspd(I,2)
              enddo
              POLYC(:,1)=0.d0
              nn1=idum1-J
              call STP(nn1,x1(J+1:idum1),xx(J+1:idum1),
     &                mm,POLYC(:,1),IFLAG)
              do j=J+1,idum1
              xi=x1(j)
              call fitVALUE(xi,POLYC(:,1),mm+1,dum1)
              write(*,'(1x,I6,3(1X,a,E16.5))') 
     &           j,'Z=',x1(j),'spd2=',xx(j),'spd_poly2=',dum1
              enddo
              do i=0,nlbv2-1
              fd_lbv2(i)=POLYC(nlbv2-i,1)
              enddo
            endif
!----------------------------------
! --- poly
!----------------------------------
          else  !ios=0
            close(50) 
            if(nrho>mflame.or.nrho>mflame.or.
     &         nlbv1>mflame.or.nlbv2>mflame) then
              call FFRABORT
     &  (1,'nrho>mflame.or.nrho>mflame.or.nlbv1>mflame.or.nlbv2>mflame')
            endif
            allocate(fd_rho(0:nrho-1),stat=ierr)
            if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_rho") 
            allocate(fd_tmp(0:ntmp-1),stat=ierr)
            if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_tmp") 
!
!
            fd_rho=0.0d0
            fd_tmp=0.0d0
!
! --- Read Density(Inverse) Polynomial Coefficients
!
            do i=0,nrho-1
            fd_rho(i)=rho(i+1)
            enddo
! --- Read a number of Temperature Polynomial Coefficients
            do i=0,ntmp-1
            fd_tmp(i)=tmp(i+1)
            enddo
!

            if(calgeq) then
!
              allocate(fd_lbv1(0:nlbv1-1),stat=ierr)
              if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_lbv1") 
              allocate(fd_lbv2(0:nlbv2-1),stat=ierr)
              if(ierr/=0) call FFRABORT(1,"ERR: allocating fd_lbv2")

              fd_lbv1=0.0d0
              fd_lbv2=0.0d0

              do I=1,3
              if(lbv_Xi(I)<0.d0.and.I/=3) then
                write(ifle,*) 'ERR: NOT defined lbv_Xi(i), i=',I
                write(ifle,*) 'ERR: Reading file',trim(flamelet_file)
                call FFRABORT(1,'ERR: Reading error')
              endif
              lbv_XXi(I)=lbv_Xi(I)
              enddo
!
              do i=0,nlbv1-1
              if(lbv1(i+1)==undef) then
                write(ifle,*) 'ERR: NOT defined lbv1(i), i=',I
                write(ifle,*) 'ERR: Reading file',trim(flamelet_file)
                call FFRABORT(1,'ERR: Reading error')
              endif
              fd_lbv1(i)=lbv1(i+1) 
              enddo 
!
              if(lbv_Xi(3)>0.d0) then
                do i=0,nlbv2-1
                if(lbv2(i+1)==undef) then
                  write(ifle,*) 'ERR: NOT defined lbv2(i), i=',I
                  write(ifle,*) 'ERR: Reading file',trim(flamelet_file)
                  call FFRABORT(1,'ERR: Reading error')
                endif
                fd_lbv2(i)=lbv2(i+1) 
                enddo 
              endif
!
            endif
          endif
        elseif(calxi) then
!------------------------------------------
! --- Other flamelet LES combustion model  
!------------------------------------------
          comb_mdl=date_base
          flmletfile=trim(adjustl(flamelet_file))
          if(flmletfile=='') then
            call FFRABORT
     &       (1,'Set flamelet date base file &flamelet/flamelet_file')
          else
            write(ifll,'(1x,2a)') 'File name :',
     &         trim(adjustl(flamelet_file))
          endif

          inquire(file=flmletfile,exist=fexist)
          if (.not.fexist) then
              write(ifle,*) 
     &       'MSG: FLAMELET DATA FILE file [',
     &        trim(flmletfile),'] NOT exist'
              call FFRABORT(1,'Please check FLAMELET DATA FILE name')
          endif
          CALL checkExpire(2022,9,30,23,59,59,lokflm)
            if(lokflm/=0) then
              write(ifle,'(1x,a)') 
     &     'MSG: DENTYUU-KEN SUPPORTTED DATE-BASE FILE HAS BEEN Expired'
              call FFRABORT(1,'Flamelet date-base file')
            endif
!            rewind(50) 
!            open(50,file=flmletfile,iostat=ios) 
            if(comb_mdl==1) then
              if(idifflm==FltRans.or.idifflm==consF) then
              else
                call FFRABORT(1,'MSG: comb_mdl=1: BlogB=1 or BlogB=4')
              endif
            elseif(comb_mdl==2) then
              if(idifflm==Blog.or.idifflm==BlogM) then
              else
                call FFRABORT(1,'MSG: comb_mdl=2: BlogB=2 or BlogB=3')
              endif
            endif
!
!          if(comb_mdl==1) then
!             open(50,file=flmletfile,form='unformatted',iostat=ios)
!          elseif(comb_mdl==2) then
!             open(50,file=flmletfile,form='unformatted',iostat=ios)
!          endif
!watnabe-sann
           if(frmt==2) then
             open(50,file=trim(flmletfile),form='unformatted',
     &         status='old',iostat=ios) 
             read(50,ERR=100) n1,n2,n3,nvar
             allocate(x1(n1),xx(n1),x2(n2),x3(n3))
             allocate(table(n1,n2,n3,nvar))
             allocate(tags(nvar))
             read(50) x1,x2,x3
             read(50) tags
             read(50) table
             close(50)
!BABA-sann
           else
             open(50,file=trim(flmletfile),form='unformatted',
     &         status='old',iostat=ios) 
             read(50,ERR=100) n1,n2,n3,nvar
             allocate(x1(n1),xx(n1),x2(n2),x3(n3))
             allocate(table(n1,n2,n3,nvar))
             allocate(tags(nvar))
             read(50) x1,x2,x3
             read(50) table
             read(50) tags
             close(50)
           endif
!
            write(ifll,'(1x,a)') 'MSG: Flamelet date-base table :'
            do i=1,nvar
            write(ifll,'(1x,a,I4,4X,a)') 'MSG: ',I,trim(tags(I))
            enddo
!
            x1min=minval(x1)
            x1max=maxval(x1)
            x2min=minval(x2)
            x2max=maxval(x2)
            x3min=minval(x3)
            x3max=maxval(x3)
!
            write(*,'(A,I5)') 'nVars ',nvar
            write(*,'(A, 2ES15.5)') 
     &       'chemtable: x1', minval(x1), maxval(x1)
            write(*,'(A, 2ES15.5)')
     &        'chemtable: x2', minval(x2), maxval(x2)
            write(*,'(A, 2ES15.5)') 
     &       'chemtable: x3', minval(x3), maxval(x3)
            do i=1, nvar 
            write(*,'(1x,A,A, 2ES15.5)')  'chemtable:  ', tags(i),             
     &                               minval(table(:,:,:,i)),             
     &                               maxval(table(:,:,:,i))
            enddo 
!
            DO J=1,nvar 
            idum1=1
            if('RHO'== tags(J)) then 
              J_RHO=J 
              do I=1,n1 
              xx(I)=table(I,1,idum1,J_RHO) 
              enddo 
              nm=1 
              dum2=-1.d0
              do I=1,n1
                dum1=xx(I)
                if(dum1>dum2) then
                  nm=I
                  dum2=dum1
                endif
              enddo
              EXIT
            endif
            enddo
!
            nm=30
            nn1=nm
            nn2=n1-nn1
!
            call STP(nn1,x1,xx,mm,POLYC(:,1),IFLAG)
            call STP(nn2,x1(nn1+1:n1),xx(nn1+1:n1),mm,POLYC(:,2),IFLAG)
!--------------
! --- check
!--------------
            write(ifll,'(1X,a)') ' Poly for RHO'
            dum2=0.d0
            do j=1,nn1
            xi=x1(j)
            call fitVALUE(xi,POLYC(:,1),mm+1,dum1)
            dum2=dum2+(dum1-xx(j))**2
            write(*,'(1x,I6,3(1X,a,E16.5))') 
     &         j,'Z=',x1(j),'rho=',xx(j),'rho_p=',dum1
            enddo
            write(*,'(1x,A,3x,E16.5)') 'RMS=',dum2

            dum2=0.d0
            do j=nn1+1,n1
            xi=x1(j)
            call fitVALUE(xi,POLYC(:,2),mm+1,dum1)
            dum2=dum2+(dum1-xx(j))**2
            write(*,'(1x,I6,3(1X,a,E16.5))') 
     &       j,'Z=',x1(j),'rho=',xx(j),'rho_p=',dum1
            enddo
            write(*,'(1x,A,3x,E16.5)') 'RMS=',dum2
            write(ifll,'(1X,a)') ' end Poly for '
!
          endif
        endif
!----------------------------------------------------------
! --- User defined scalar must be put end of system model 
!----------------------------------------------------------
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
          ike_c=no
        endif
!
!
!      do no=1,nransno
!        write(ifll,*) '    ### Scalar Information               ###' 
!        write(ifll,*) '    ### Scalar Sequential Order    : ',no
!        write(ifll,*) '    ### Scalar Name                : ',
!     &                     TRIM(adjustl(sclname(no)))
!        write(ifll,*) '    ### Scalar Usr&Sys Defined No. : ',scl(no) 
!        write(ifll,2000)
! 2000   format(3x,45('-'),2x)
!      enddo
!
      if(nrans.ne.nransno) then
        write(ifle,*) ' ### ERR: nrans= ',nrans,' nransno= ',nransno
        write(ifle,*) ' ### MSG: Please contact your supportor'
        write(ifle,*) ' ### MSG: Or Re-run prefflow '
        call FFRABORT (1,'[nrans] NOT equal to [nransno]')
      endif
!
      if(nransno.gt.0) then
        lrans=.true.
        rns_scl=.true.
      else
        lrans=.false.
        rns_scl=.false.
      endif
!
      if(nransno.gt.mrans) then
        write(ifll,*) '    ### ERR : Salar number cannot great than ',
     &   mrans
        call FFRABORT (1,'module_scalar/inputdata')
      endif
!
!
!
!
!--------------------
! --- Potential scalar
!--------------------
      Potential_name(:)=' '
      rewind ifli
      read(ifli,Potential,iostat=ios)
      if(ios.gt.0) stop ': at reading error at &Potential'
      iterPOTEN=max(0,Iter_poten)
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
      NPOTNno=no
!
      if(NPOTNno/=NPOTN) then
        write(ifle,*) ' ### ERR: NPOTN= ',NPOTN,' NPOTNno= ',NPOTNno
        write(ifle,*) ' ### MSG: Please contact your supportor'
        write(ifle,*) ' ### MSG: Or Re-run prefflow '
        call FFRABORT (1,'[NPOTN] NOT equal to [NPOTNno]')        
      endif
      allocate(potname(npotn),stat=ierr)
      if(ierr.ne.0) stop ': allocate error for [potname]'  
!---------------------
! --- initialization
!---------------------
      potname(:)=' '
      ip_pefc(:)=1
!
      no=0
      if(lFC>0) then
        no=no+1
        potname(no)='FAI_S'
        ip_pefc(1)=no
!
        no=no+1
        potname(no)='FAI_M'
        ip_pefc(2)=no
      endif
!-----------------
      do i_user=1,mpotn
      if(Potential_name(i_user).ne.' ') then
        no=no+1
        potname(no)=Potential_name(i_user)
      endif
      enddo
!
      if(NPOTNno.gt.0) then
         pot_scl=.true.
         lpotn=.true.
      else
         pot_scl=.false.
         lpotn=.false.
      endif
!
      if(lFC>0) then
        iPEFC=lFC
      endif
!
      return
!
 100  CALL FFRABORT(1,'ERR: reading [flamelet_file] file error')
 120  CALL FFRABORT(1,'ERR: reading [read(50,ERR=120) nvar] error')
 130  CALL FFRABORT(1,'ERR: reading [flamelet_func] error')
 200  CALL FFRABORT(1,'ERR: reading [reading ORGtable(i,j)] error')

      end subroutine inputdata
!
!=====================================================================
      subroutine fitVALUE(xi,fldat,nfl,phi)
!=====================================================================
      integer,intent(in)   :: nfl
      real*8,intent(in)    :: xi
      real*8,intent(inout) :: phi
      real*8,intent(in)    :: fldat(nfl)
      integer :: ifl
!
      phi=0.0d0
      do ifl=1,nfl
      phi=phi+fldat(ifl)*xi**dble(nfl-ifl)
      enddo
!
      end subroutine fitVALUE
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine STP(M,X,Y,N,C,IFLAG)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      implicit none

!
! --- [dummy arguments]
!
      integer,intent(inout) :: M,IFLAG
      integer,intent(in)    :: N
      real*8, intent(inout) :: X(M),Y(M),C(N+1)
!
! --- [local entities]
!
!      integer,parameter :: MA=20
      integer           :: I,J,K,N1,N2,NNPN,N1PN1,N1N1,IFLAG1
      integer           :: I1J,L,JNP2,IN1,IJ,IM1JP1,NC,N1N2
      real*8            :: XPWR
      real*8,allocatable :: A(:)
!
!
      N1N2=N1*N2
!
      N1=N+1
      N2=N1+1

      ALLOCATE(A(N1*N2))

      NNPN=N*N1
      N1PN1=N1+N1
      N1N1=NNPN+N1
      
      I1J=1
      do 30 J=1,N
      A(I1J)=0.d0
      I1J=I1J+N1
 30   enddo

      L=NNPN+1
      do 40 J=1,N1PN1
      A(L)=0.d0
      L=L+1
 40   enddo

      do 10  K=1,M
      XPWR=X(K)
      I1J=N1+1
      do 20 J=2,N1
      JNP2=N1N1+J
      A(I1J)=A(I1J)+XPWR
      A(JNP2)=A(JNP2)+XPWR*Y(K)
      XPWR=XPWR*X(K)
      I1J=I1J+N1
 20   enddo

      do 50 I=2,N1
      IN1=NNPN+I
      A(IN1)=A(IN1)+XPWR
      XPWR=XPWR*X(K)
 50   enddo

      A(N1N1+1)=A(N1N1+1)+Y(K)

 10   enddo

      A(1)=dble(M)
      do 60 I=2,N1
      IJ=I
      do 80 J=1,N
      IM1JP1=IJ+N
      A(IJ)=A(IM1JP1)
      IJ=IJ+N1
 80   enddo
 60   enddo

      NC=-N2
      
      CALL MATI(IFLAG,IFLAG,N1,NC,A,N1,N1N2,XPWR)
      K=N1N1+N1
      DO 70 I=1,N1
      C(I)=A(K)
      K=K-1
 70   enddo

      return
!
      end subroutine STP
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine MATI(ISOL,IDSOL,NR,NC,A,MRA,MA,DET)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
       implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: ISOL,IDSOL,NR,NC,MRA
      integer,intent(in) :: MA
      real*8, intent(inout) :: A(MA),DET

      
!
! --- [local entities]
!
      integer :: KWA(NR)
      integer :: IR,IC,IBMP,JBMP,KBMP,NES,MDIV,IRIC,MAD,MSER,KSER
      integer :: MZ,II,JJ,KK
      integer :: NET,iflag,I,IP,J,JR,K,L
      real*8  :: dum1,PIV,PSTO
!
      IR=NR
      ISOL=1
      IDSOL=1

      IF(NR<=0.or.(IR-MRA)>0) then
        ISOL=3
        return
      endif

      IC=IABS(NC)
      IF((IC-IR)<0) IC=IR
      IBMP=1
      JBMP=MRA
      KBMP=JBMP+IBMP
      NES=IR*JBMP
      NET=IC*JBMP
      
      IF(NC==0) then
        ISOL=3
        return
      endif
      
      IF(NC<0) then
        MDIV=JBMP+1
        IRIC=IR-IC
      else
        MDIV=1
      endif

      MAD=MDIV
      MSER=1
      KSER=IR
      MZ=1
      DET=1.d0
 15   continue

      IFLAG=0
      PIV=0.d0
      I=MSER
 17   IF(I<=KSER) then
        if((ABS(A(I))-PIV)>0.d0) then
          PIV=ABS(A(I))
          IP=I
        ENDIF
        I=I+IBMP
        goto 17
      endif

      IF(PIV==0.d0) then
        DET=0.d0
        ISOL=2
        IDSOL=1
        return
      endif

      IF(NC<0) then
        I=IP
        J=MSER
      else
        I=IP-((IP-1)/JBMP)*JBMP
        J=MSER-((MSER-1)/JBMP)*JBMP
        JJ=MSER/KBMP+1
        II=JJ+(IP-MSER)
        KWA(JJ)=II
      endif


      IF(IP<MSER) then
        ISOL=3
        return
      endif


 27   IF(IP>MSER) then
        IF((J-NET)<=0) then
          PSTO=A(I)
          A(I)=A(J)
          A(J)=PSTO
          I=I+JBMP
          J=J+JBMP
          goto 27
        endif
        DET=-DET
      endif

      PSTO=A(MSER)
      DET=DET*PSTO
      IF(DET==0.d0) then
        IDSOL=1
        ISOL=2
        return
      endif
      PSTO=1.d0/PSTO

      A(MSER)=1.d0

      I=MDIV

 37   IF((I-NET)<=0) then
        A(I)=A(I)*PSTO
        I=I+JBMP
        goto 37
      endif

 47   IF((MZ-KSER)<=0) then
        IF((MZ-MSER)/=0) then
          I=MAD
          J=MDIV
          PSTO=A(MZ)
           IF(PSTO/=0.d0) then
              A(MZ)=0.d0
 57           IF((J-NET)<=0) then
                A(I)=A(I)-A(J)*PSTO
                J=J+JBMP
                I=I+JBMP
                GOTO 57
              endif
           endif
        endif
        MAD=MAD+IBMP
        MZ=MZ+IBMP
        goto 47
      endif

      KSER=KSER+JBMP

      IF(KSER<=NES) then
        MSER=MSER+KBMP
        IF(NC>=0) then
          MDIV=MDIV+IBMP
          MZ=((MSER-1)/JBMP)*JBMP+1
          MAD=1
          IFLAG=15
        endif
        IF(IFLAG/=15) then
          MDIV=MDIV+KBMP
          IF(IRIC/=0) then
            MZ=((MSER-1)/JBMP)*JBMP+1
          else
            MZ=MSER+IBMP
          endif
          MAD=MZ+JBMP
          IFLAG=15
        endif
      endif

      IF(IFLAG==15) goto 15
      IF(NC<0) return

      JR=IR

 67   if(JR/=0) then
         if(JR==0.or.KWA(JR)<JR) then
            ISOL=3
            return
         endif

         IF(KWA(JR)>JR) then
           K=(JR-1)*JBMP
           J=K+IR
           L=(KWA(JR)-1)*JBMP+IR
 77        IF(J>=K) then
             IF(J>K) then
               PSTO=A(L)
               A(L)=A(J)
               A(J)=PSTO
               J=J-IBMP
               L=L-IBMP
             endif
             GOTO 77
           endif
           JR=JR-1
         endif
         goto 67
      endif

      return

      end subroutine  MATI

      end module module_scalar
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_rans
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(LEN=13),parameter,private :: modnam='(module_rans)'
!
! --- < 1. module for RANS model >
!
      integer,private :: nrans=0
!
      real(8),save,allocatable :: akslw(:)
      real(8),save,allocatable :: sgmaks(:)
      real*8,save,allocatable  :: akshg(:)
      real*8,save,allocatable  :: sgmaksm(:)
!
!
! --- akslw  : lower limit of depedent variables 
! --- sgmaks : parameter in diffusion term (sigama for eddy diff)
! --- akshg  : high limit of depedent variables 
! --- sgmaksm: parameter in diffusion term (sigama for molecule diff)
!
!< 2. for k-epsilon model >
!
      real(8) :: cep1=1.44d0, cep2=1.92d0, cep3=1.3d0,cep4=-0.33d0,
     &           cep5=0.25d0
      real(8) :: cmu=0.09d0
      real(8) :: sgmk=1.0d0, sgme=1.3d0,sgmh=0.09d0,sgmm=0.09d0,
     &           eta0=4.38d0,beta=0.012d0
      real(8) :: akdes,cb1,cb2,sgmDES,cw1,cw2,cw3,cv1,cDES
      real(8) :: Schmidt_G=0.25d0, Schmidt_Z=0.5,
     &           Schmidt_G_M=1.d10,Schmidt_Z_M=0.7
!
      real(8) :: L2s=2.d-3
!
      real(8) :: Ck,Ce
!
! --- cep1 : parameter in generation term of eps. eq. [-]
! --- cep2 : parameter in generation term of eps. eq. [-]
! --- cep3 : parameter in generation term of eps. eq. [-]
! --- cmu  : parameter in turbulent viscosity [-]
! --- sgmk : parameter in diffusion term of k eq. [-]
! --- sgme : parameter in diffusion term of eps. eq. [-]
!
!-< 99. for internal check >-
!
      integer,private :: iset=0
!
!///////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namelist >-------------------------------------------------
!=======================================================
      subroutine inputke(ifli,ifll,ifle,cntlnam,
     &          nrnsx,kemdl,lowke,l2SKE,lRNG,lCHEN,lSDES,lKLES,E2P,
     &                   calgeq,calxi,caltemp,calh,ierror)
!=======================================================
!
      use module_scalar,only : icalke,icalrsm,icalaph,ical_s,ike,irsm,
     &                         ivof,ides,icavi,ivold,ike_c,ixi,igeq,
     &                         itemp,ihflmlt,iaph
      use module_model,only  : u_func
!
      implicit none
!
! --- MXMATERIAL
!
!
! ---- [dummy arguments]
!
      integer,intent(in)  :: ifli,ifll,ifle,nrnsx
      character(*),intent(in) :: cntlnam
      logical,intent(in)  :: kemdl,lowke,lRNG,lCHEN,lSDES,l2SKE,
     &                       calgeq,calxi,caltemp,calh,E2P,lKLES
      integer,intent(out) :: ierror
!
! --- [namelist] 
!
      real(8) :: klow,epslow,sgm_C_dash=1.d0
      real(8) :: C_k=0.05d0,C_e=1.d0
!
      namelist /kemodel/ klow,epslow,cep1,cep2,cep3,cep4,cep5,
     &                   cmu,sgmk,sgme,sgmm,sgmh,eta0,beta,
     &                   Schmidt_G,Schmidt_Z,
     &                   Schmidt_G_M,Schmidt_Z_M,
     &                   sgm_C_dash,C_k,C_e
     &                   
!
      real(8) :: akppa
      namelist /DESmodel/ akppa,cb1,cb2,sgmDES,cw1,cw2,cw3,cv1,cDES
!
! ---- [local entities]
!
      character(LEN=9),parameter :: subnam='(inputke)'
      real(8),parameter :: undef=-huge(1.d0)
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
!          CALL FFRABORT(1,'module_rans/inputke')
!        endif
      else
        return
      endif
!
      allocate( akslw(nrans), sgmaks(nrans), stat=ierr1 )
      akslw(:)=1.d-20
      sgmaks(:)=1.d0
      
!
      if(E2P) then
        akslw(iaph(1))=0.d0
        akslw(iaph(2))=0.d0
        sgmaks(iaph(1))=1.d0
        sgmaks(iaph(2))=1.d0
      endif
      if( ierr1.ne.0 ) goto 9999
!
      if(kemdl.or.l2SKE) then
        cep1=1.44d0
        cep2=1.92d0
        cep3=1.3d0
        cep4=-0.33d0
        cmu=0.09d0
        sgmk=1.0d0
        sgme=1.3d0
        if(l2SKE) then
           L2s=2.d0
        endif
      elseif(lowke) then
        cep1=1.00d0
        cep2=1.83d0
        cep3=1.3d0
        cep4=-0.33d0
        cmu=0.084d0
        sgmk=1.69d0
        sgme=1.3d0
        
      elseif(lRNG) then
        cep1=1.42d0
        cep2=1.68d0
        cep3=1.42d0
        cep4=-0.387d0
        cmu=0.085d0
        sgmk=0.719d0
        sgme=0.719d0
        sgmh=0.09d0
        sgmm=0.09d0
        eta0=4.38d0
        beta=0.012d0
      elseif(lCHEN) then 
        cmu=0.09d0
        sgmk=0.75d0
        sgme=1.15d0
        sgmh=0.09d0
        sgmm=0.09d0
        cep1=1.15d0
        cep2=1.9d0
        cep3=1.4d0
        cep4=-0.33d0
        cep5=0.25d0
      elseif(lSDES) then
        akdes=0.41d0
        cb1=0.1355d0
        cb2=0.622d0
        sgmDES=0.6667d0
        cw1=cb1/akdes**2+(1.d0+cb2)/sgmDES
        cw2=0.3d0
        cw3=2.d0
        cv1=7.1d0
        cDES=0.65d0
      elseif(lKLES) then
        Ck=0.05d0
        Ce=1.d0
      endif
!
      
      allocate(akshg(nrans),sgmaksm(nrans),stat=ierr1)
      akshg(:)=1.d30
      sgmaksm(:)=1.d0
!
      if(ierr1.ne.0) goto 9999
!
! --- T.Tominaga Add 2004/12/22 START
      if(calgeq) then
         akslw  (igeq) = 0.d0
         akshg  (igeq) = 1.d0

         sgmaks (igeq) = Schmidt_G
         sgmaksm(igeq) = Schmidt_G_M
      endif
      if(calxi) then
         akslw  (ixi) = 0.0d0
         akshg  (ixi) = 1.0d0

         sgmaks (ixi) = Schmidt_Z
         sgmaksm(ixi) = Schmidt_Z_M
      endif
      if(caltemp) then
         akslw  (itemp) = 0.0d0
         akshg  (itemp) = 1.0d0
         sgmaks (itemp) = 0.5d0
         sgmaksm(itemp) = 0.7d0
      endif
      if(calh) then
         akslw  (ihflmlt) = 0.0d0
         akshg  (ihflmlt) = 1.0d0
         sgmaks (ihflmlt) = 0.7d0
         sgmaksm(ihflmlt) = 0.75d0
      endif
!      
      if(u_func(1)==1) then
        sgmaks(ike_c) = sgm_C_dash 
      endif
!
      klow  =undef
      epslow=undef
      sgm_C_dash=1.d0
      rewind ifli
      read(ifli,kemodel,iostat=ios)
      call nml_errmsg0(ifle,ios,'kemodel',ierr1)
      if(ios.gt.0) call FFRABORT(1,'&kemodel')
!
      if( klow.eq.undef ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'ERR: lack of data : klow'
        call FFRABORT(1,'&kemodel')
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
        call FFRABORT(1,'&kemodel')
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
      if( cep5.lt.0.d0.and.(lCHEN) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'cep5 = ',cep5
        write(ifle,*) 'it must be >= 0'
        goto 9999
      endif
!
      if( cep5.lt.0.d0 ) then
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
      if( sgmh.le.0.d0.and.(lRNG.or.lCHEN) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'sgmh = ',sgmh
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if(sgmm.le.0.d0.and.(lRNG.or.lCHEN)) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'sgmm = ',sgmm
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if( eta0.le.0.d0.and.(lRNG) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'eta0 = ',eta0
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if( beta.le.0.d0.and.(lRNG) ) then
        write(ifle,*) '### error : data error'
        write(ifle,*) 'beta = ',beta
        write(ifle,*) 'it must be > 0'
        goto 9999
      endif
!
      if(kemdl.or.lowke.or.lRNG.or.lCHEN) then
        akslw (ike(1))=klow
        akslw (ike(2))=epslow
        sgmaks(ike(1))=sgmk
        sgmaks(ike(2))=sgme
      endif
!
      if(lKLES) then
        if(C_k<0.d0) call FFRABORT(1,'ERR: &kemodel/C_k for KLES model')
        if(C_e<0.d0) call FFRABORT(1,'ERR: &kemodel/C_e for KLES model')
      else
        Ck=C_k
        Ce=C_e
      endif
!
      if(lSDES) then
        akppa=undef
        cb1=undef
        cb2=undef
        sgmDES=undef
        cw2=undef
        cw3=undef
        cv1=undef
        cDES=undef
        rewind ifli
        read(ifli,DESmodel,iostat=ios)
        call nml_errmsg0(ifle,ios,'DESmodel',ierr1)
        if(ios.gt.0) call FFRABORT(1,'&DESmodel')
        if(ios<0) then
          akdes=0.41d0
          cb1=0.1355d0
          cb2=0.622d0
          sgmDES=0.6667d0
          cw1=cb1/akdes**2+(1.d0+cb2)/sgmDES
          cw2=0.3d0
          cw3=2.d0
          cv1=7.1d0
          cDES=0.65d0
        else
          if(akppa==undef) then
            akdes=0.41d0
          else
            akdes=akppa
          endif
          if(cb1==undef) then
            cb1=0.1355d0
          endif
          if(cb2==undef) then
            cb2=0.622d0
          endif
          if(sgmDES==undef) then
            sgmDES=0.6667d0
          endif
          cw1=cb1/akdes**2+(1.d0+cb2)/sgmDES
          if(cw2==undef) then
            cw2=0.3d0
          endif
          if(cw3==undef) then
            cw3=2.d0
          endif
          if(cv1==undef) then
            cv1=7.1d0
          endif
          if(cDES==undef) then
            cDES=0.65d0
          endif
        endif
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

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_FFRGFwords
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
!------------------------------------------------------------
! --- Data format identifier --------------------------------
!------------------------------------------------------------
      character*80,parameter :: asciiv2_FFGF="#A_GF_V2"
      character*80,parameter :: unformv2_FFGF="#U_GF_V2"
      character*80,parameter :: binaryv2_FFGF="#B_GF_V2"
!------------------------------------------------------------
! --- Dataset identifier ------------------------------------
!------------------------------------------------------------
      character*80,parameter :: newSet_FFGF="#NEW_SET"
      character*80,parameter :: fileend_FFGF="#ENDFILE"
      character*80,parameter :: customData_FFGF="#CUSTOM_DATA"
      integer,parameter :: unkownDatasize_FFGF=0
!------------------------------------------------------------
! --- FFR-GF Datafield Keywords
!     For grids.
! -----------------------------------------------------------
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
!     For time series data.
      character*80,parameter :: ffrTseries_FFGF="FFR_TSERIES_FIELD"
      character*80,parameter :: tseriesHead_FFGF="#FFR_TSERIES_HEADER"
      character*80,parameter :: tseriesHeade_FFGF=
     &                                    "#FFR_TSERIES_HEADER_END"
      character*80,parameter :: tseriesObserv_FFGF=
     &                                    "#FFR_TSERIES_OBSERVER"
      character*80,parameter :: tseriesObserve_FFGF=
     &                                    "#FFR_TSERIES_OBSERVER_END"
      character*80,parameter :: tseriesData_FFGF="#FFR_TSERIES_DATA"
      character*80,parameter :: tseriesDatae_FFGF=
     &                                    "#FFR_TSERIES_DATA_END"
!--------------------
! --- For text data.
!--------------------
      character(80),parameter :: ffrText_FFGF="FFR_TEXT_FIELD"
      character(80),parameter :: textData_FFGF="#FFR_TEXT_DATA"
      character(80),parameter :: textDatae_FFGF="#FFR_TEXT_DATA_END"
!     For Axis data
      character(80),parameter :: ffrAxis_FFGF="FFR_AXIS_FIELD"
      character(80),parameter :: axisData_FFGF="#FFR_AXIS_DATA"
      character(80),parameter :: axisDatae_FFGF="#FFR_AXIS_DATA_END"
!     For Position data
!      character(80),parameter :: ffrAxis_FFGF="FFR_POSITION_FIELD"
!      character(80),parameter :: PosiData_FFGF="#FFR_POSITION_DATA"
!      character(80),parameter :: PosiDatae_FFGF="#FFR_POSITION_DATA_END"

      end module module_FFRGFwords
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_cdcl
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_cdcl)'
!
! --- [module arguments]
!
      logical :: cdcloutputFlag
! Flag of Cd/Cl value calculation.
      integer :: numofCdClWall
! Number of walls.
      character*80,allocatable :: nameofCdClWall(:)
! List of wall name.
      integer,allocatable :: CdClWall_bcNo(:)
! Boundary region index number of wall.
      real(8),allocatable :: projAreaS(:)
! wall area projected to inlet velocity direction
      real(8) :: inflowVelCdCl(1:3)
!      
!///////////////////////////////////////////////////////////////////////
      contains
!=======================================================================
      subroutine inputdata
     & (ifli,ifll,ifle,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,ierror)
!=======================================================================
      implicit none
!
!
! --- [dummy qrgument]
!
      integer,intent(in)    :: ifli,ifll,ifle,nbcnd,kxnone,kdfire,kdintr
     &                              ,kdbuff,kdcvd
      integer,intent(in)    :: kdbcnd(0:7,nbcnd)
      character(*),intent(in) :: boundName(nbcnd,2)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=3),parameter :: cond(2)=(/'yes','no '/)
      character(LEN=7),parameter :: co(2)=(/'.true. ','.false.'/)
      integer :: ierr1,ios
      integer :: iset=0
      integer :: kcdcl=0

      integer :: ic,jc,kc,nb,kd
      character*80,allocatable :: nameofCdClWallTmp(:)
!
! --- [namelist]
!
      character(20) :: enable_wall
      character(80) :: wall_name
      real(8) :: inflow_u,inflow_v,inflow_w
      namelist /cdcl_output/enable_wall,wall_name,
     &                        inflow_u,inflow_v,inflow_w
!
!
      ierr1=0
      ierror=0
      cdcloutputFlag=.true.
!
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
! --- Cd/Cl wall setting FROM HERE ------------------------------------------
!     Count up the number of enabled walls.
      numofCdClWall=0
      rewind(ifli)
      do while(.true.)
        enable_wall='yes'
        read(ifli,cdcl_output,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'cdcl_output',ierr1)
        if(ierr1/=0) then
          ierror=1
          cdcloutputFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_wall,kcdcl)
        if(kcdcl==1) numofCdClWall=numofCdClWall+1
      end do
!     If there is no wall, set calculation flag to .false.
      if(numofCdClWall<=0) then
        cdcloutputFlag=.false.
        return
      end if
      
!     Prepare the list of wall name.
      allocate(nameofCdClWallTmp(1:numofCdClWall))
      
!     Rereading the namelist file.
      rewind(ifli)
      ic=0
      inflow_u=0.d0
      inflow_v=0.d0
      inflow_w=0.d0
      do while(.true.)
        enable_wall='yes'
        wall_name=' '
        read(ifli,cdcl_output,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'cdcl_output',ierr1)
        if(ierr1/=0) then
          ierror=1
          cdcloutputFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_wall,kcdcl)
        if(kcdcl==1) then
          ic=ic+1
          nameofCdClWallTmp(ic)=wall_name
        end if
      enddo
!
!     inlet flow velocity for calculating CD/CL
!
      inflowVelCdCl(1)=inflow_u
      inflowVelCdCl(2)=inflow_v
      inflowVelCdCl(3)=inflow_w
!     Clear duplicated wall names and recount the number of walls.
      if(ic>=2) then
        do jc=1,ic-1
          if(nameofCdClWallTmp(jc)/=' ') then
            do kc=jc+1,ic
              if(nameofCdClWallTmp(kc)==nameofCdClWallTmp(jc))
     &            nameofCdClWallTmp(kc)=' '
            end do
          end if
        end do
        numofCdClWall=0
        do jc=1,ic
          if(nameofCdClWallTmp(jc)/=' ') 
     &            numofCdClWall=numofCdClWall+1
        end do
      end if
      allocate(nameofCdClWall(1:numofCdClWall))
      if(ic>=2) then
        kc=0
        do jc=1,ic
          if(nameofCdClWallTmp(jc)/=' ') then
            kc=kc+1
            nameofCdClWall(kc)=nameofCdClWallTmp(jc)
          end if
        end do
      else
        nameofCdClWall(:)=nameofCdClWallTmp(:)
      end if
      deallocate(nameofCdClWallTmp)

!     Search relation between walls 
!      and boundary region index number.
      allocate(CdClWall_bcNo(1:numofCdClWall))
      CdClWall_bcNo(:)=0
      do jc=1,numofCdClWall
        do nb=1,nbcnd
          if(adjustl(nameofCdClWall(jc))==adjustl(boundName(nb,1))) 
     &      CdClWall_bcNo(jc)=nb
        end do
      end do
!     Check the list.
      do jc=1,numofCdClWall
        if(CdClWall_bcNo(jc)==0) then
          write(ifle,*)'### error : namelist [cdcl]'
          write(ifle,*)trim(nameofCdClWall(jc)),
     &                 ' not found in boundary namelist'
          call FFRABORT(1,'module_cdcl/inputdata')
        endif
      end do
!     Confirm the boundary conditions of walls.
      do jc=1,numofCdClWall
        nb=CdClWall_bcNo(jc)
        kd=kdbcnd(0,nb)
        if(.not.(kd.eq.kxnone.or.
     &           kd.eq.kdfire.or.
     &           kd.eq.kdintr.or.
     &           kd.eq.kdbuff.or.
     &           kd.eq.kdcvd)) then
          write(ifle,*)
     &   ' ERR: The boundary : ',trim(nameofCdClWall(jc)),
     &   ' is NOT wall boundary'
          call FFRABORT(1,'module_flow_cdcl/inputdata')
        endif
      end do
      
!      
! ---  wall setting TO HERE ---------------------------------------

      end subroutine inputdata
      end module module_cdcl
!

!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_probe
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: 
     &                        modnam='(module_probe)'
!
! --- [module arguments]
!
      logical :: probeFlag 
      ! Flag of probe data output enable/disable.
      integer :: numofProbes
      ! Number of output points.
      real(8),allocatable :: probePosit(:,:)
      real(8),allocatable :: probePositIn(:,:)
      integer,allocatable :: probeCVIndex(:,:),NO_YS(:)
      integer,allocatable :: ave_step(:)
      
      ! (:,1) : CV number   (:.2) : cpu number
      ! Positions of output points.
      integer,allocatable :: probeDataType(:,:)
      ! output data type and number of components. (1:numofProbes,1:2)
      !    (:,1)        (:,2)
      !  1:Pressure     1
      !  2:Velocity     3
      !  3:Temparature  1
      !  4:Density      1
      !  5:pp0          1
      character*80,allocatable :: probeDataLabel(:)
      ! output data label (space character is not allowed in label string)
      real(8),allocatable :: probeData(:,:)
      ! output data
!      integer :: ProbeFiletype
!      ! Output File Format
!      !  1: GF(default)
!      !  2: Excel(space paused)

!
!///////////////////////////////////////////////////////////////////////
      contains
!=======================================================================
      subroutine inputdata
     & (ifli,ifll,ifle,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,ierror)
!=======================================================================
      implicit none
!
!
! --- [dummy qrgument]
!
      integer,intent(in)    :: ifli,ifll,ifle,nbcnd,kxnone,kdfire,kdintr
     &                              ,kdbuff,kdcvd
      integer,intent(in)    :: kdbcnd(0:7,nbcnd)
      character(*),intent(in) :: boundName(nbcnd,2)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='(inputdata)'
      character(LEN=3),parameter :: cond(2)=(/'yes','no '/)
      character(LEN=7),parameter :: co(2)=(/'.true. ','.false.'/)
      character(LEN=5),parameter :: cfile(2)=(/'GF   ','Excel'/)
      integer :: ierr1,ios
      integer :: iset=0
      integer :: kprobe=0
      integer :: ic,jc,kc,nb,kd
!
! --- [namelist]
!
      character(LEN=20) :: enable_probe
      real(8) :: position_x,position_y,position_z
      character(LEN=3)   :: identifier
      character(LEN=80)  :: label
      integer            :: ys_no=0
      integer            :: average_step=1
      namelist /probe/enable_probe,position_x,position_y,
     &                position_z,identifier,label,ys_no
     &                ,average_step
      
!
!
      ierr1=0
      ierror=0
      probeFlag=.true.
!      filetype='GF'
!
!
      call modutl_chkset('=',1,iset,modnam//subnam)
!
! --- Probe setting FROM HERE --------------------------------
!     Count up the number of enabled probes
      numofProbes=0
      rewind(ifli)
      do while(.true.)
        enable_probe='yes'
        read(ifli,probe,iostat=ios)
        if(ios<0) then
          exit
        endif
        call nml_errmsg0(ifle,ios,'probe',ierr1)
        if(ierr1/=0) then
          ierror=1
          probeFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_probe,kprobe)
        if(kprobe==1) numofProbes=numofProbes+1
      end do
!     If there is no probe, set calculation flag to .false.
      if(numofProbes<=0) then
        probeFlag=.false.
        return
      end if
!----------------------------------      
!     Preare arrays of probe data.
!----------------------------------
      allocate(probePosit(1:numofProbes,1:3))
      allocate(probePositIn(1:numofProbes,1:3))
      allocate(probeCVIndex(1:numofProbes,1:2))
      allocate(probeDataType(1:numofProbes,1:2))
      allocate(probeDataLabel(1:numofProbes))
      allocate(NO_YS(1:numofProbes))
      allocate(ave_step(1:numofProbes))
!----------------------------------
!     Reload name list file.
!----------------------------------
      ic=0
      rewind(ifli)
      do while(.true.)
        enable_probe='yes'
        position_x=0.d0;position_y=0.d0;position_z=0.d0
        label=' '
        read(ifli,probe,iostat=ios)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'probe',ierr1)
        if(ierr1/=0) then
          ierror=1
          probeFlag=.false.
          write(ifle,*) modnam,subnam
          return
        end if
        call nml_listno(2,cond,enable_probe,kprobe)
        if(kprobe==1) then
          ic=ic+1
          probePosit(ic,1)=position_x
          probePosit(ic,2)=position_y
          probePosit(ic,3)=position_z
          probeDataLabel(ic)=trim(adjustl(label))
          if((identifier=='prs').or.(identifier=='PRS')) then
            probeDataType(ic,1)=1
            probeDataType(ic,2)=1
            if(probeDataLabel(ic)==' ')
     &            probeDataLabel(ic)='Pressure'
          else if((identifier=='vel').or.(identifier=='VEL')) then
            probeDataType(ic,1)=2
            probeDataType(ic,2)=3
            if(probeDataLabel(ic)==' ')
     &            probeDataLabel(ic)='Velocity_uvw'
          else if((identifier=='tmp').or.(identifier=='TMP')) then
            probeDataType(ic,1)=3
            probeDataType(ic,2)=1
            if(probeDataLabel(ic)==' ')
     &            probeDataLabel(ic)='Temperature'
          else if((identifier=='den').or.(identifier=='DEN')) then
            probeDataType(ic,1)=4
            probeDataType(ic,2)=1
            if(probeDataLabel(ic)==' ')
     &            probeDataLabel(ic)='Density'
          else if((identifier=='pp0').or.(identifier=='PP0')) then
            probeDataType(ic,1)=5
            probeDataType(ic,2)=1
            if(probeDataLabel(ic)==' ')
     &            probeDataLabel(ic)='Thermal_pressure'
          else if((identifier=='ys').or.(identifier=='YS')) then
            probeDataType(ic,1)=6
            probeDataType(ic,2)=1
            if(probeDataLabel(ic)==' ')
     &            probeDataLabel(ic)='Mass_Fraction'
            if(ys_no==0) then
              call FFRABORT(1,'ERR: NOT defined Ys number')
            else
              NO_YS(ic)=ys_no
            endif
          else if((identifier=='avp').or.(identifier=='AVP')) then
            if(average_step<=0) then
              call FFRABORT(1,'ERR: wrong average step')
            else
              ave_step(ic)=average_step
            endif
            probeDataType(ic,1)=7
            probeDataType(ic,2)=average_step
            if(probeDataLabel(ic)==' ')
     &            probeDataLabel(ic)='Pressure(averaged)'
          end if
!
!          call nml_listno(2,cfile,filetype,kprobe)
!          if(kprobe==2) then
!             ProbeFiletype=2
!          else
!             ProbeFiletype=1
!          end if
!
        end if
      end do
! --- Observer points setting TO HERE ----------------------------------

! --- allocate output data buffer --------------------------------------
      allocate(probeData(1:numofProbes,
     &                  1:maxval(probeDataType(1:numofProbes,2))))

      end subroutine inputdata
      end module module_probe
!


!
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_mphase_pc
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! idrift     :    =0  drift flux model
!                 =1  two-fluid model
! i_contm    :    =1  for the contaminated liquid
!                 =2  for the slightly contaminated liquid
!                 =3  for the pure liquid
! i_pc     :      =0  no phase change
!                 =1  with phase change
! i_gsrc     :    =0  no gas source(aeration)
!                 =1  with gas source(aeration)
! dia        :    the diameter of bubble
! sigma      :    surface tension between gas & liquid[N/m]
! fn         :    a factor for calculating cd of a swarm of bubbles
!
      integer    :: ip_cnv=0           !! considering the convective term or not
      integer    :: ip_vis=0           !! considering the viscous term or not
      integer    :: ipresbf=1          !! for the calculation of body force 
      integer    :: ipresrhs=1         !! for the calculation of rhs
      integer    :: iaout=0
      integer    :: i_vtlim=0          !! flag for whether vt is limited or not
      integer    :: irchow=2           !! irchow=0 : no Rhie-Chow interpolation
                                       !! irchow=1 : Rhie-Chow interpolation for liquid only                                       !! irchow=2 : Rhie-Chow interpolation for liquid & gas
      integer    :: i_timnew=0         !! limit to time-step
      integer    :: iht_old=0          !! =1 -> the old one is used for calculating ht coef.
!
!! ak
!
      integer    :: idrift=1           !! two-fluid model
      integer    :: i_contm=1          !! for the contaminated liquid
      integer    :: isalp=2            !! one-eq or two-eq for the vol. fractions
      integer    :: imueconst=0        !! constant eddy viscosity
      integer    :: ike2p=1
      integer    :: itke=1
      integer    :: iadf=1
!
!! pv
!
      integer    :: i_poisson=30       !! flag for Poisson equation
      integer    :: ipbc=2             !! flag for the pressure BC
      integer    :: iv_pea=0           !! =1 : partial elimination algorithm
      integer    :: itvmin=1           !! min iteration no. for the interphase coupling
      integer    :: itvmax=10          !! max iteration no. for the interphase coupling
      real(8)    :: veps=2.0d-4        !! velocity tolerance for the interphase coupling
      integer    :: ivdf=1
!
!! temp
!
      integer    :: i_pc=1             !! =0: no phase change; =1: with phase change
      integer    :: i_cnds=0           !! =0: no condensation; =1: with condensation
      integer    :: icaltmp=1
      integer    :: ibchtb=1           !! calculating the boiling heat transfer coef
      integer    :: ihtbwrn=1          !! output of wrn msg when calculating htb
      integer    :: itdf=1
      integer    :: iydf=1
      integer    :: istintf=1          !! implicitly treat the interface or not
      real(8)    :: dtgup=5.0          !! 
      integer    :: ihpres=1           !! considering the variation of pressure or not
!
!! gradient
!
      integer    :: igradp=1           !!      
      integer    :: igradv=1           !!      
      integer    :: igradt=0           !!      
      integer    :: igrada=1           !!      
      integer    :: igradk=0           !!      
      integer    :: igrado=0           !!      
!
!! gridcells
!
!!      integer    :: icase=-1           !! =-1 : F-type model; -2 : Cube: others: input
!!      integer    :: ncvx0=7008         !! the number of cells for icase=[others]
!
!! parameters
!
      real(8)    :: d_b0=0.005         !! the diameter of gas bubbles
      real(8)    :: d_d0=0.001         !! the diameter of liquid droplets
      real(8)    :: sigma=0.072        !! air-water surface tension
      real(8)    :: fn=1.75            !! parameter for the bubble swarm effect
      real(8)    :: alp_1=0.30         !! lower limit to alpha for 2-phase flow pattern
      real(8)    :: alp_2=0.75         !! upper limit to alpha for 2-phase flow pattern
      real(8)    :: rat_pc=0.5         !! the ratio for evaporation
      real(8)    :: omega=1.0          !! the under relaxation fator for pressure correction(0.0 - 1.0)
!
!! gassrc
!
      integer    :: i_gsrc=0           !! no aeration
      real(8)    :: xl=-1.0d+20        !! the range of gas generation region
      real(8)    :: yl=-1.0d+20        !! the range of gas generation region
      real(8)    :: zl=-1.0d+20        !! the range of gas generation region
      real(8)    :: xu= 1.0d+20        !! the range of gas generation region
      real(8)    :: yu= 1.0d+20        !! the range of gas generation region
      real(8)    :: zu= 1.0d+20        !! the range of gas generation region
      real(8)    :: gsrc=1.0d-3        !! gas source(+)/sink(-) [m^3/(m^3.sec)]
!
!! limit
!
      integer    :: ilim_a=1           !! limit to alpha
      integer    :: ilim_m=1           !! limit to mass_i
      integer    :: ilim_cd=0          !! limit to cd
      integer    :: ilim_hl=0          !! limit to h(enthalpy)
      integer    :: ilim_hg=1          !! limit to h(enthalpy)
      integer    :: ilim_tmpl=0        !! limit to temp
      integer    :: ilim_tmpg=0        !! limit to temp
      integer    :: niter0=-10000000   !!
!
!! dbg
!
      integer    :: idbg_fric=0
      integer    :: idbg_ke=0          !! wr k-e
      integer    :: idbg_m=0           !! wr mass_i
      integer    :: idbg_nrans=0       !! wr nrans
      integer    :: idbg_sam=0         !! wr sol_alp
      integer    :: idbg_kesrc=0       !! ke src
      integer    :: idbg_eul=0         !! ke eulerian
      integer    :: idbg_alp=0         !! output of alpha for debug
      integer    :: idbg_rho=0         !! output of rho for debug
      integer    :: idbg_tmp=0         !! output of temp for debug
      integer    :: idbg_vt=0          !! output of v_term for debug
      integer    :: idbg_vel=0         !! output of vel for debug(cal convective term)
      integer    :: idbg_prs=0         !! output of parameters in prs_admin
      integer    :: i_coord=0          !! output of z-coordinates
!
      real(8)    :: qtemps,qt_wfrate,qtemp,qtemp1,qtemp2
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine input_param
!=================================================
      implicit none
!
! --- [dummy arguments]
!
      integer  :: iu_2
      integer  :: ios,ierror
      logical  :: l_exist
!
! --- [namelist]
!
      namelist /flag/ idrift,i_contm,i_pc,i_cnds,i_gsrc,i_poisson,
     &                ip_cnv,ip_vis,isalp,imueconst,i_vtlim,itvmin,
     &                itvmax,iaout,icaltmp,ike2p,itke,ibchtb,ihtbwrn
      namelist /parameters/ d_b0,d_d0,sigma,fn,veps,alp_1,alp_2
      namelist /gassrc/ xl,yl,zl,xu,yu,zu,gsrc
      namelist /limit/ ilim_a,ilim_m,ilim_cd,
     &                 iht_old
      namelist /debug/ idbg_fric,idbg_ke,idbg_m,idbg_nrans,
     &                 idbg_sam,idbg_kesrc,idbg_eul,idbg_alp,
     &                 idbg_rho,idbg_tmp,idbg_vt,idbg_vel,idbg_prs,
     &                 i_coord
!
! --------------------------------------------------------------
      inquire(file='fort.2', exist=l_exist, iostat=ios)
! --------------------------------------------------------------
!
      
      if(l_exist) then
        iu_2=2
!
! --------------------------------------------------------------
        open(unit=iu_2,file='fort.2',iostat=ios)
! --------------------------------------------------------------
!
        read(iu_2,flag,iostat=ios)
!
        read(iu_2,parameters,iostat=ios)
!
        read(iu_2,gassrc,iostat=ios)
!
        read(iu_2,limit,iostat=ios)
!
        read(iu_2,debug,iostat=ios)
!
        close(iu_2)
!
      end if
!
      if(i_gsrc == 0) gsrc=0.0
!
	return
!
      end subroutine input_param
!
      end module module_mphase_pc                     
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_rad
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character(20),parameter,private :: modnam='(module_rad)'
      INTEGER,private :: I,J,K,L 
!
!
!< module for radiation >
!
      character*10::	radmodel='FVM'	!'ZONE'/'MC'/'FVM'
      character*10::	gasmodel='SLW'	!'SLW'/'WSGG'/'EWBCK'/'SNBCK'
      integer::		radflag=0  	!do or do not
      REAL*8 ::		vfeps=1.0d-6
      REAL*8 ::		feps=1.0d-5
      REAL*8 ::		RLX=0.75d0
      REAL*8 ::		teaeps=1.0d-5
      REAL*8 ::		charlen=0.d0	!characteristic lenth = 3.6V/F

      INTEGER::		iteaflag=0
      INTEGER::		iModHF=1
      INTEGER::		nray=0
      INTEGER::		radgpopt=1
      INTEGER::		mgpnum=6000
      INTEGER::		tearst=0
      INTEGER::		NDIV=8
      INTEGER::		DOMSN=1
      INTEGER::         RadJump=1
      INTEGER::         ngauss=4

c	radflag=0  	!do or do not
c	vfeps=1.0E-8	!relative residual for smoothing exchange-area
c	feps=1.0E-8	!absolute residual for smoothing exchange-area
c	RLX=0.75	!relaxition coef. for GS interation
c	NDIV=1		!partition of solid angle for FVM/DOM/DTM models
c	teaeps=1.0E-5	!residual for solving Jiang's Model for 
c			!			anisotropic scattering media
c	beta=0		!coef. for anisotropic scattering, =0 means 
c			!			isotropic media
c	iteaflag=0	!flag for calculate total-exchange-area
c	iModHF		!mode for heat transfer calculate,
	  		!1: direct cal., 2: by solving equation


!parameter available flag
      INTEGER :: ParaFlag(0:40)
!for slw
      REAL*8 :: almn(0:3,0:3,0:3,40),blmn(0:3,0:3,0:3,40),
     &	 cminmax(2,40)
!for wsgg
      REAL*8 :: awsgg(4,6), bwsgg(4,4,5), cwsgg(4,4,3,5)
      REAL*8 :: dwsgg(6,3), TotCoef(100)
!I----n, J----m, K----l
!H2O
      DATA   (((almn(I,J,K,1),I=0,3),J=0,3),K=0,3)
     &		/1.6103d0,-4.0931d0, 5.1435d0, -2.0857d0,
     &		-0.81812d0, 15.5525d0,-21.819d0, 9.8775d0,
     &		2.6001d0, -21.204d0, 31.0828d0, -14.279d0,
     &		-1.3171d0, 9.6524d0, -14.474d0, 6.6747d0,
     &		0.440187d0, -0.63348d0, 0.871627d0, -0.38798d0,
     &		-0.82164d0,5.0239d0,-5.9818d0,2.6355d0,
     &		1.5149d0,-7.8032d0,9.8642d0,-4.1931d0,
     &		-0.81023d0, 3.727d0, -4.874d0, 1.9868d0,
     &		0.106647d0, -0.43116d0, 0.689598d0, -0.29831d0,
     &		-0.38573d0, 1.8865d0, -2.9712d0, 1.2834d0,
     &		0.578351d0, -2.6218d0, 4.2698d0, -1.7929d0,
     &		-0.28014d0, 1.1785d0, -1.9568d0, 0.787249d0,
     &		8.25027d-3, -3.28556d-2, 6.81563d-2, -3.04815d-2,
     &		-3.10578d-2, 0.123369d0, -0.26154d0, 0.117452d0,
     &		4.39319d-2, -0.15792d0, 0.350948d0, -0.15308d0,
     &		-2.03699d-2, 6.61142d-2, -0.15283d0, 6.34035d-2/
       DATA   (((blmn(I,J,K,1),I=0,2),J=0,3),K=0,3)
     &		/4.72d0, -8.5482d0, 5.2394d0,
     &		-0.84969d0, 0.312478d0, -0.13804d0,
     &		-3.47243d-2, 4.02461d-2,-5.80104d-2,
     &		5.79830d-4, 3.94125d-3, -5.29017d-3,
     &		-8.9615d0, 16.9547d0, -10.76d0,
     &		1.5861d0, -2.0166d0, 1.46d0,
     &		4.3473d-2, -0.67133d0, 0.633231d0,
     &		2.87067d-3, -7.0683d-2, 6.2371d-2,
     &		9.1461d0, -17.327d0, 11.1864d0,
     &		-1.3975d0, 1.9965d0, -1.6935d0,
     &		8.46419d-2, 0.599994d0, -0.70054d0,
     &		7.14719d-3, 6.62086d-2, -6.87294d-2,
     &		-3.5504d0, 6.624d0, -4.3058d0,
     &		0.485392d0, -0.7071d0, 0.689109d0,
     &		-6.77456d-2, -0.18179d0, 0.269308d0,
     &		-5.92726d-3, -2.04694d-2, 2.56411d-2/	
      DATA   cminmax(1,1) / 3.0d-5 /
      DATA   cminmax(2,1) / 60.d0 /
!CO2
      DATA   (((almn(I,J,K,2),I=0,3),J=0,3),K=0,3)
     &		/2.45702d0, -5.45334d0, 6.53751d0, -2.52344d0,
     &		-4.0232d0, 15.67297d0, -24.3247d0, 11.33757d0,
     &		7.54549d0, -23.8023d0, 39.51896d0, -19.1137d0,
     &		-3.63104d0, 11.9078d0, -20.3606d0, 9.97877d0,
     &		7.65678d-2, 2.36184d0, -3.95061d0, 2.17482d0,
     &		0.2901819d0, -12.0041d0, 22.44342d0, -13.0467d0,
     &		-0.64282d0, 21.5003d0, -40.8667d0, 23.66762d0,
     &		0.3942158d0, -11.5818d0, 22.05176d0, -12.6536d0,
     &		-3.30582d-2, 0.4367742d0, -0.725331d0, 0.4138566d0,
     &		0.3672993d0, -3.52466d0, 6.74885d0, -3.96295d0,
     &		-0.69811d0, 6.60703d0, -12.9667d0, 7.58713d0,
     &		0.3831158d0, -3.65683d0, 7.19415d0, -4.16496d0,
     &		-1.87927d-3, 1.92123d-2, -3.25863d-2, 1.98493d-2,
     &		2.85033d-2, -0.223537d0, 0.4402715d0, -0.26267d0,
     &		-5.49594d-2, 0.4370937d0, -0.881494d0, 0.521958d0,
     &		3.04198d-2, -0.247793d0, 0.4990777d0, -0.291566d0/
      DATA   cminmax(1,2) / 3.0d-5 /
      DATA   cminmax(2,2) / 600.d0 /
!absorb coefficient / emissive coefficient
      DATA   ((awsgg(I,J),I=1,3),J=1,6)
     &	  /	0.3966d0,  15.64d0,   394.3d0,
     &		0.4098d0,  6.325d0,   120.5d0,
     &		0.4496d0,  7.113d0,   119.7d0,
     &		0.4303d0,  7.055d0,   178.1d0,
     &		0.4201d0,  6.516d0,   131.9d0,
     &		1.2531d0,  8.4258d0,  87.064d0 /
!weight of emissivity
      DATA   (((bwsgg(I,J,k),I=1,4),J=1,3),K=1,5)
     &	  /	0.4334d-1,  2.62d-4,   -1.56d-7,  2.565d-11,
     &		-0.4814d-1, 2.822d-4,  -1.794d-7, 3.274d-11,
     &		0.5492d-1,  0.1087d-4, -0.35d-7,  0.9123d-11,
     &		5.977d-1,   -5.119d-4, 3.042d-7,  -5.564d-11,
     &		0.5677d-1,  3.333d-4,  -1.967d-7, 2.718d-11,
     &		1.8d-1,     -2.334d-4, 1.008d-7,  -1.454d-11,
     &		6.324d-1,   -8.358d-4, 6.135d-7,  -13.03d-11,
     &		-0.2016d-1, 7.145d-4,  -5.212d-7, 9.868d-11,
     &		3.5d-1,     -5.04d-4,  2.425d-7,  -3.888d-11,
     &		5.15d-1,    -2.303d-4, 0.9779d-7, -1.494d-11,
     &		0.7749d-1,  3.399d-4,  -2.297d-7, 3.77d-11,
     &		1.907d-1,   -1.824d-4, 0.5608d-7, -0.5122d-11,
     &		6.508d-1,   -5.551d-4, 3.029d-7,  -5.353d-11,
     &		-0.2504d-1, 6.112d-4,  -3.882d-7, 6.528d-11,
     &		2.718d-1,   -3.118d-4, 1.221d-7,  -1.612d-11 /
      DATA   ((dwsgg(I,J),I=1,6),J=1,3) 
     &	  /	1.6879d-1,  2.5682d-4, 9.5161d-8,
     &		-3.166d-10,1.4834d-13,-2.156d-17, 
     &		4.9577d-2, 9.3954d-4, -1.6416d-6,
     &		1.1478d-9, -3.76d-13, 4.7503d-17, 
     &		2.789d-1, -5.1265d-4, 6.732d-7, 
     &		-5.1488d-10,1.8887d-13, -2.5856d-17 /
!wdight of absorbility for I=4,5
      DATA   ((((cwsgg(I,J,K,L),I=1,4),J=1,4),K=1,3),L=4,5)
     &	/	0.55657,    -0.62824d-3, 0.31876d-6,   -0.52922d-10,
     &		0.32964d-3, 0.27744d-6,  -0.26105d-9,  0.37807d-13,
     &		-0.53441d-6,0.33753d-9,  -0.10348d-12, 0.26027d-16,
     &		0.12381d-9, -0.90223d-13,0.38675d-16,  -0.99306d-20,
     &	    0.16676d-1, 0.15769d-3,  -0.10937d-6,  0.19588d-10,
     &		0.5091d-3,  -0.76773d-6, 0.40784d-9,   -0.69622d-13,
     &	    0.3762d-7,  0.18729d-9,  -0.15889d-12, 0.30781d-16,
     &		-0.3251d-10,-0.26171d-13,0.29848d-16,  -0.58387d-20,
     &		0.28689d-1, 0.20697d-3,  -0.17474d-6,  0.37238d-10,
     &		0.24221d-3, -0.55686d-6, 0.34884d-9,   -0.67887d-13,
     &		-0.19492d-6,0.36102d-9,  -0.2148d-12,  0.41305d-16,
     &		0.41721d-10,-0.73d-13,  0.431d-16,    -0.83182d-20,
     &		0.59324d0,    -0.61741d-3, 0.29248d-6,   -0.45823d-10,
     &		0.35739d-3, 0.22122d-6,  -0.2638d-9,   0.45951d-13,
     &		-0.71313d-6,0.46181d-9,  -0.70858d-13, 0.38038d-17,
     &		0.17806d-9, -0.11654d-12,0.19939d-16,  -0.13486d-20,
     &		-0.35664d-1,0.21502d-3,  -0.13648d-6,  0.24284d-10,
     &		0.51605d-3, -0.70037d-6, 0.3868d-9,    -0.70429d-13,
     &		0.12245d-6, 0.99434d-10, -0.15598d-12, 0.37664d-16,
     &		-0.57563d-10,-0.10109d-13,0.35273d-16, -0.89872d-20,
     &		0.12951,    0.5452d-4,   -0.80049d-7,  0.17813d-10,
     &		0.1521d-3,  -0.3775d-6,  0.21019d-9,   -0.36011d-13,
     &		-0.13165d-6,0.20719d-9,  -0.9672d-13,  0.14807d-16,
     &		0.26872d-10, -0.34803d-13,0.14336d-16, -0.19754d-20  /
!//////////////////////////////////////////////////////////////////////
      contains
!
!< #1. input namdlist >------------------------------------------------
!=================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,rg_markx,ierror)
!=================================================
      implicit none
      character(11),parameter :: subnam='(inputdata)'
      
!----------------------
! [dummy arguments]
!----------------------
      integer,intent(in)      :: ifli,ifll,ifle,rg_markx
      character(*),intent(in) :: cntlnam
      
      integer,intent(out)     :: ierror
!----------------------
! [namelist]
!----------------------
      character*10 :: rmodel
      integer :: rflag
      namelist /radoption/ 
     &                   rflag,
     &                   rmodel,
     &                   gasmodel,

     &                   feps,  ! NOT
     &                   vfeps,  ! NOT
     &                   RLX,  ! NOT
     &                   nray,  ! NOT
     &                   radgpopt,  ! NOT

     &                   domsn,
     &                   mgpnum,  ! NOT
     &                   tearst,  ! NOT
     &                   radjump,
     &                   ngauss,
     &                   charlen

!----------------------
! [local entities]
!----------------------
      real*8,parameter :: undef=-huge(1.d0)
      integer :: iset=0
      integer :: ios,ierr1
!
!-----------------------------------------------------------------------
! --- 
!-----------------------------------------------------------------------
      call modutl_chkset('=',1,iset,modnam//subnam)
      ierror=0
!
      rewind ifli
      read(ifli,radoption,iostat=ios)
      if( ios.lt.0 ) return
      call nml_errmsg0(ifle,ios,'radoption',ierr1)
      if( ierr1.ne.0 ) goto 9999

      radflag=rflag
      if(radflag.ne.1) return
 1001 IF(radflag.eq.1) then
        if((rmodel(1:4).ne.'ZONE').and.(rmodel(1:2).ne.'MC').and.
     &    (rmodel(1:3).ne.'FVM')) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 'model = ',trim(rmodel)
          write(ifle,*) 'it must be <ZONE> <MC> <FVM>'
          goto 9999
        endif
      ENDIF
!
      radmodel=rmodel
      if(radmodel(1:4).NE.'ZONE') iteaflag=0
      IF(iteaflag.EQ.1.OR.radmodel(1:4).NE.'ZONE') THEN
  	iModHF=1
      ELSE
  	iModHF=2
      ENDIF
!	 
      IF(radmodel(1:3).eq.'FVM') THEN
	IF( (DOMSN.NE.2).AND.(DOMSN.NE.4).AND.(DOMSN.NE.6).AND.
     &	(DOMSN.NE.8) ) THEN
	  write(ifle,*) 'ERR: Direction zone defining error'
	  write(ifle,*) 'domsn = ',domsn
	  write(ifle,*) 'it should be <2,4,6, or 8>'
	call FFRABORT(1,'ERR: &radoption/domsn')
      ENDIF
!
      NDIV=DOMSN*(DOMSN+2)
!
      ENDIF
!
      if(
     &   gasmodel(1:3).ne.'SLW'.AND.
     &	 gasmodel(1:4).ne.'WSGG'.AND.
!     &	 gasmodel(1:5).ne.'EWBCK'.AND.
!     &	 gasmodel(1:5).ne.'SNBCK'
     &	 gasmodel(1:4).ne.'GRAY'
     &   ) then
	write(ifle,*) '### error : Real_Gas_Model ERROR'
        write(ifle,*) 'Available gas models are : '
        write(ifle,*) '<SLW> <WSGG> <GRAY> '
        call FFRABORT(1,'ERR: ')
      endif
!
      if(rg_markx==1) then
        if(
     &   gasmodel(1:3).ne.'SLW'.AND.
     &	 gasmodel(1:4).ne.'WSGG'
     &   ) then
          write(ifle,*) 
     &   'MSG: [gastype=2] defined in &fluid for real gas'
          write(ifle,*) 'Available gas models are : '
          write(ifle,*) "gasmodel='WSGG' ,or 'SLW'"
          call FFRABORT(1,'ERR: [GRAY]')
        endif
      elseif(rg_markx==0) then
        if(gasmodel(1:4)/='GRAY') then
          write(ifle,*) 
     &   'MSG: [gastype=1] defined in &fluid for gray gas'
          write(ifle,*) 'Available gas models are : '
          write(ifle,*) "gasmodel='GRAY' "
          call FFRABORT(1,'ERR: [WSGG] or [SLW]')
        endif
      endif
!
      if(gasmodel(1:3).eq.'SLW'.and.charlen.le.0.d0) then
	write(ifle,'(1x,a)') 'ERR: Please give characteristic lenth'
        write(ifle,'(1x,a)') '      charlen = 3.6V/F'
        write(ifle,'(1x,a)') 'where V and F are total',
     &		 ' volume and surface area'
        call FFRABORT(1,'ERR: [WSGG] or [SLW]')
      endif

      if(gasmodel(1:4).eq.'SLW') then
	ParaFlag(0:2)=1
        if(NGauss<3) then
          write(ifle,'(1x,a)') 'MSG: intergration for SLW'
!          call FFRABORT(1,'ERR: NGauss<3')
        endif
        
      endif
!
      if(gasmodel(1:4).eq.'WSGG') then
        if(NGauss==0) then
          call FFRABORT(1,'ERR: NGauss=0')
        elseif(NGauss>4) THEN
          call FFRABORT(1,'ERR: NGauss>4')
        endif
	ParaFlag(0:2)=1
      endif
!
      IF(gasmodel(1:4).eq.'GRAY') THEN
        if(NGauss/=1) then
           call FFRABORT(1,'ERR: NGauss/=1')
        endif
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
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_size
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      integer,save         :: 
      implicit none
      public :: nml_sizes

!//////////////////////////////////////////////////////////////////////
!
!
      contains
!
C***********************************************************************
      SUBROUTINE nml_sizes(ifli,ifll,ifle,cntlnam,MXPRTX,iprtcle,ierror)
C===============================================================
      IMPLICIT NONE 
!
! --- [dummy arguments]
!
      integer,intent(out)     :: MXPRTX,ierror
      integer,intent(in)      :: ifli,ifll,ifle,iprtcle
      character(*),intent(in) :: cntlnam
!
      integer :: mcell,mface,ncell,nvrtx,mvrtx,ncomp,mssfbc,medge
      integer :: ncomp_surface,MXPRT
      integer :: ios,ierr1
      namelist /sizes/ mcell,mface,medge,ncell,nvrtx,ncomp,
     &                 mvrtx,mssfbc,ncomp_surface,MXPRT
!
      if(iprtcle==0) then
        MXPRTX=1
        return
      endif
!
      MXPRT=-1
      rewind ifli
      read (ifli,sizes,iostat=ios)
      call nml_errmsg0(ifle,ios,'sizes',ierr1)
      if(ierr1.ne.0) call FFRABORT(1,'ERR: read error &sizes')
      if(MXPRT==-1) then
        write(ifle,'(1X,a)') 'Maximum Particle number by MXPRT'
        call FFRABORT(1,'ERR: NOT defined [MXPRT] in &sizes')
      endif
      MXPRTX=MXPRT
      
      return
      end subroutine nml_sizes
!
      end module module_size




!
!ical_dens==4
!flamelet
!PISO
!LBF_P=3,5
!chung
!overset
