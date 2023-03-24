!      
!      SUBROUTINE FFR_TRACING 
!      SUBROUTINE FFR_TRACING01 
!      SUBROUTINE FFR_TRACING03 
!      SUBROUTINE FFR_TRACING02_F 
!      SUBROUTINE FFR_TRACING_NEWUVW 
!      SUBROUTINE NUKIYAMA_TANAZAWA  !initial  ! num. fix 
!      SUBROUTINE NUKIYAMA_TANAZAWA_VOLUMEBASE ! volume fix ,num=>free  
!      SUBROUTINE CALC_NUKIYAMATANAZAWA_AB 
!      subroutine calc_gamma 
!      SUBROUTINE add_random_walk 
!      SUBROUTINE MT_RANDOM 
!      SUBROUTINE READ_PARTICLE_INI 
!      SUBROUTINE BREAKUP_NORMAL_TAB 
!      SUBROUTINE read_particle_restart 
!      SUBROUTINE LOCATE_PARTICLE_TO_CELL 
!      SUBROUTINE BREAKUP_PILCH 
!      SUBROUTINE BREAKUP_IMPROVED_TAB 
!      SUBROUTINE write_particle_restart 
!      SUBROUTINE RE_INJ_P 
!      SUBROUTINE P_result_S 
!      NOx subroutine
!      subroutine Combustion_Char
!      subroutine Zeldovich_NO_P
!      real*8 function Zeldovich_Func
!      real*8 function Particle_PDF
!      real*8 function Zeldovich_func_PDF
!      subroutine prompt_NO 
!      subroutine FUEL_NO 
!      subroutine FUEL_NO2 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_particle 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_species,only  : a7,wm,r_wm,Ri,gascns,href,hform,
     &                           tref_comp,judg_temp,act
      use module_chemreac,only : vreq
      character(20),parameter,private :: 
     &                        modnam='(module_particle)' 
!
! --- [module arguments]
!
      real*8,save,allocatable   :: mass_p_total(:)
      integer,save,allocatable   :: lack_p(:),re_inj(:)
!
      integer, save     :: ITER_P 
      integer,parameter :: M_injctr=20,mcomp=200 
      integer           :: N_injctr=0,Num_injector=10000
      integer, save,allocatable :: NINALL(:)
      integer, save,allocatable :: NPALL(:)
      integer, save,allocatable :: NINEND(:)
      real*8,  save,allocatable :: SPHEIGHT(:),VOLUMEPC(:)
!-------------------------
! Logical Unit Number 
!--------------------------------------
      integer, parameter :: iPartLU=19 
!--------------------------------------
! Setting File Name
!-------------------------
      character*13, parameter :: cPSETTINGFNAM="particle.ctl"
      character*12, parameter :: cPTINIFNAM="sprayini.txt"
!-------------------------
! Initial Conditions
!-------------------------
      real*8, save,allocatable :: dUVWNEW(:)
      real*8, save,allocatable :: dSMD(:)
      real*8, save,allocatable :: dInjVolume(:)
      real*8, save,allocatable :: RHOP0(:)
      real*8, save,allocatable :: particle_comp(:,:)
!
!-----------------------------------------------------------------
! Spray Injection Time and Number of Injection Parcels at a STEP
!-----------------------------------------------------------------
      integer, save,allocatable :: iSprayInjStart(:)
      integer, save,allocatable :: iSprayInjEnd(:)
      integer, save,allocatable :: iSprayInjTime(:)
      real*8,  save,allocatable :: starttime(:)
      integer, save,allocatable :: iSprayInjNum(:)
      integer, save,allocatable :: iSprayInjInter(:)
      integer, save,allocatable :: iSprayInjTotal(:)
      logical, save,allocatable :: INJ_FLAG(:)
!-------------------------
! Output Control
!-------------------------
!      integer, save :: iOUTITER
!
!-------------------------
! Injection Point
!-------------------------
      real*8, save,allocatable :: dXP0(:)
      real*8, save,allocatable :: dYP0(:)
      real*8, save,allocatable :: dZP0(:)
!
      real*8, save,allocatable :: endx(:)
      real*8, save,allocatable :: endy(:)
      real*8, save,allocatable :: endz(:)
      real*8, save,allocatable :: axis_inj(:,:)!axis_inj(3,:)
!
!-------------------------
! Breakup Model 
!-------------------------
!
      character*4,save,allocatable :: cBreakupModel(:)
      logical, save,allocatable :: lBreakupNONE(:)
      logical, save,allocatable :: lBreakupPILCH(:)
      logical, save,allocatable :: lBreakupITAB(:)
      logical, save,allocatable :: lBreakupNTAB(:)
!
!-------------------------
! Sub cycling 
!-------------------------
      integer,save :: ITER_INNER
!-----------
! Restart 
!-----------
!      logical,save :: lPTRESTART
!---------------------------------------------
! Parameter for Nukiyama-Tanazawa Function 
!---------------------------------------------
      REAL*8, save,allocatable :: dAlpha(:)
      REAL*8, save,allocatable :: dBeta(:)
!---------------------------------------------
! --- Spray Geometry Settings
!---------------------------------------------
      character*4,save,allocatable:: cNozzleType(:)
      REAL*8, save,allocatable :: dSPRAYANGLE(:)
      REAL*8, save,allocatable :: dSPANGLEMAG(:)
      REAL*8, save,allocatable :: dFANANGLE(:)
      REAL*8, save,allocatable :: dSLITWIDTH(:)
      REAL*8, save,allocatable :: dSLITdir(:)
!-------------------------
! Injector TYPE 
!-------------------------
      integer, save :: imanual=1,islit=2,ihole=3,i_coal=4,
     &                 inlet=5,user=6
      integer,save,allocatable :: injkd(:)
!-------------------------
! EVAPRATION MODEL 
!-------------------------
      integer, save :: iflag_eva=0
      integer, save :: ievaN=0,ievaF=1,iAbSi=2,iLaKn=3,ievaU=4
      integer, save :: ihetD=1,ihetMass=2,ihetLaKn=3,ihetU=4
      integer, save :: imovD=1,imovMilr=2,imovSomm=3,imovU=4,imovMrsi=5,
     &                 imovSub
!
!AbSi=2 :Abramzon  Sirignano model
!LaKn=3 : Langmuir-Knudsen model
!
      integer,save,allocatable :: evakd(:,:)
      integer,save,allocatable :: heatexkd(:)
      integer,save,allocatable :: movkd(:)
      REAL*8,save,allocatable  :: Latent(:,:)
      REAL*8,save,allocatable  :: PP_0(:,:)
      REAL*8,save,allocatable  :: TT_0(:,:)
      REAL*8,save,allocatable  :: EVAP_COMP(:,:)
      REAL*8,save,allocatable  :: PP_sat(:,:),TT_vap(:,:),TT_boil(:,:)
!----------------------------
! --- Breakup Model   
!----------------------------
      real*8,save,allocatable :: dPDstion(:)
      real*8,save,allocatable :: dPDstvel(:)
      real*8,save,allocatable :: dLiquidMu(:)
      real*8,save,allocatable :: dLiquidStens(:)
      real*8,save,allocatable :: dCf(:)
      real*8,save,allocatable :: dCk(:)
      real*8,save,allocatable :: dCd(:)
      real*8,save,allocatable :: dCb(:)
      real*8,save,allocatable :: dCv(:)
      real*8,save,allocatable :: dK(:)
      real*8,save,allocatable :: ddcd(:)
!
      integer, save :: iconst=1,iRosR=2,inoml=3,inuki=4,ikuro=5
      integer,save,allocatable :: dropkd(:)
      real*8,save,allocatable  :: dRR_a(:)
!----------------------------
!--- particle_information   
!----------------------------
      real*8, save :: dMinRelaxTime
      integer,save :: iMaxPassedCell
      real*8, save :: dMaxRandomWalk
      real*8, save :: dAveRandomWalk
      integer,save :: iTotalCountOfRW
!----------------------------
! --- temp array 
!----------------------------
      integer, allocatable,save :: IPSND(:),IPQUENE(:),
     &                             NPNEW(:),NSTA(:),NRES(:),NRER(:),
     &                             NP0(:),NDEL(:),NPLEFT(:),NDEL1(:),
     &                             NP_REINJ(:),NP_DEL(:)
      integer, allocatable,save :: IPCPU(:),W1K(:)
      REAL*8 , allocatable,save :: TRESIDUAL(:),HIST(:),WK(:,:)
      REAL*8 , allocatable,save :: HHIST(:,:)
!
      integer, allocatable,save :: NPEXCHAGOUT(:),NPEXCHAGIN(:) 
      INTEGER                   :: NPOUTALL,NPINALL

      integer, allocatable,save :: JPASS(:)
      integer, allocatable,save :: JUDGE(:),BCKND(:)
!---------------------
! --- coal
!---------------------
      integer,save      :: IFLAG_COAL=0
      integer,save :: total_start
      
      INTEGER,PARAMETER :: VS_CHO=1,V_CHO=2,
     &                     V_H2O=3,VS_H2O=4,
     &                     V_HCN=5,VS_HCN=6,
     &                     DIMV=6
      integer,parameter :: icoal=1,iliquid=2,iglass=3,isolid=4
      integer,parameter :: idensity=1,ivolume=2
      real*8,parameter  :: M_CO2=44.010d-3,
     &                     M_H2O=18.016d-3,
     &                     M_O2=32.d-3
!     &                     Tref=298.15d0
      integer,save      :: II_CHO,II_N2,II_CO2,II_CO,II_H2O,
     &                     II_NO,II_N,II_O,II_H,II_OH,II_HCN,
     &                     II_CN,II_C,II_S,II_ASH,II_O2
!
      integer,allocatable ,save :: ifuel(:)
      integer,allocatable ,save :: itype(:)
      real*8, allocatable ,save :: Q_fac(:),
     &                             MF_H2O_p(:),MF_vola_p(:),
     &                             MF_fixC_p(:),MF_ash_p(:)
      real*8, allocatable ,save :: MF_C_e(:),MF_H_e(:),
     &                             MF_O_e(:),MF_N_e(:),MF_S_e(:)
!
      real*8, allocatable ,save :: a_C(:),b_H(:),c_O(:)
      real*8, allocatable ,save :: Q_MF_vola(:),
     &                             Q_MF_fixC(:),
     &                             Q_MF_H2O(:),
     &                             Q_MF_ASH(:),
     &                             net_MF_HCN(:),
     &                             net_MF_Nchar(:)
      real*8, allocatable ,save :: DH_CHO(:),HR_cha(:),
     &                             LH_CHO(:),LH_H2O(:)
!
      real*8, allocatable ,save :: M_CHO(:),H_CHO(:)
      real*8, allocatable ,save :: xx_O2(:)
      real*8, allocatable ,save :: H_CO2a(:),H_H2Oa(:),H_O2a(:)
!
      real*8,parameter          :: sigm_setbolt=5.67D-8  !W/(m^2-K^4)
      real*8, allocatable ,save :: eps_rad(:)
!
      real*8, allocatable ,save :: Tini_p(:)
      real*8, allocatable ,save :: Cp_pa(:)  !CL
      real*8, allocatable ,save :: htc_p(:)  !J/(m^2 s K)
      real*8, allocatable ,save :: heat_f(:),Tref_p(:)
!
      real*8, allocatable ,save :: Avap_CHO(:),Evap_CHO(:)
      real*8, allocatable ,save :: Avap_H2O(:),Evap_H2O(:)
      real*8, allocatable ,save :: Avap_HCN(:),Evap_HCN(:)
!
      integer,allocatable ,save :: 
     &                             ical_la_coal(:),
     &                             ical_evap(:)
!
      real*8, allocatable ,save :: Acomb_C(:),Ecomb_C(:)
      real*8, allocatable ,save :: FF_char(:),eeta_N_NO(:) 
      real*8, allocatable ,save :: p00_CHAR(:)
      real*8, allocatable ,save :: AA_re_NO(:),EE_re_NO(:),AA_E(:) 
      real*8, allocatable ,save :: AA_promp(:),EE_promp(:) 
!
      real*8,allocatable, save :: 
     &           AA_no1(:),
     &           ff_index(:),
     &           EE_no1(:),
     &           AA_no2(:),
     &           EE_no2(:)
!
!
!
      integer, allocatable ,save :: NO_BC(:),P_MF(:)
      integer, allocatable ,save :: bc_nb(:)
! --- 
      REAL*8,  save,allocatable  :: PMASSC(:)
! ---
      integer, allocatable ,save :: partreac(:),INJ_RNUM(:)
      integer,save,allocatable   :: IDX_PRTRAC(:,:)
      
!---------------------------
! --- particle chemreac 
!---------------------------
      integer,save      :: ical_P_R=0 ! =0 : NO particle reaction
                                      ! =1 : particle reaction
      integer,parameter :: 
     &                     P_CHAR_COM   =1,
     &                     P_CHAR_NO    =2,
     &                     P_HCN_NO_N2O2=3,
     &                     P_elemnt     =4,
     &                     P_ovrall     =5,
     &                     P_user       =6
!
      integer :: nqP=0,nneqP                     ! number of equation
!
      integer,save,allocatable :: ig_iterP(:)    ! ignition iter for every reaction
      integer,save,allocatable :: stp_iterP(:)   ! stop reaction for every reaction
!
      integer,save,allocatable :: ireqP(:)       ! reaction model
      real(8),save,allocatable :: vreqP(:,:,:)   ! stoichiometric coefficients
      integer,save,allocatable :: slidP(:,:,:)   ! stoichiometric coefficients
      real(8),save,allocatable :: preqP(:,:,:)   ! parameters for reaction rate
      real(8),save,allocatable :: creqP(:,:)     ! index number of mol
      real(8),save,allocatable :: sub_vreqP(:,:) ! subtraction of "vreq" of each
!                                                ! species in each equation
      integer,save,allocatable :: sub_slidP(:,:) 
      real(8),save,allocatable :: sigma_vreqP(:) ! change of "sub_vreq" in
                                                 ! each equation
!
      real(8),save,allocatable :: nu_intnst(:),rain_vel(:)
      real(8),save,allocatable :: Rh_d(:),func_d(:),pcl_WDR(:)
!///////////////////////////////////////////////////////////////////////
      contains
!=======================================================================
      subroutine inputdata(ifli,ifll,ifle,cntlnam,restart,
     &                     iprtcle,iinj,ncomp,mneq,l_temp,lcomp,
     &                     my_rank,coal_CHO,nneq,
     &                     ireac,igT_iter,
     &                     spinam,ierror)
!=======================================================================  
      USE module_dimension,only : HPC_dead,MXCVFAC,
     &                            A_dead,BC_dead,P_Reced,MASS_dead
      use module_boundary,only  : nobcnd,LBC_INDEX,nbcnd
      use module_metrix,only    : lreac,lreacP
      use module_model,only     : Rh
!
      implicit none
!
! --- [dummy qrgument] 
!
      integer,intent(in)      :: ifli,ifll,ifle,iprtcle,ncomp,
     &                           restart,
     &                           l_temp,my_rank,coal_CHO(mneq),nneq,mneq
      
! --- I_CHO,I_H2O,I_C,I_ash,I_O2
      character(*),intent(in) :: cntlnam
      character(*),intent(in) :: spinam(ncomp)
      integer,intent(out)     :: ierror,iinj
      integer,intent(in)      :: ireac,igT_iter
      logical,intent(in)      :: lcomp
!
! --- [namelist] !injector
!
      integer,parameter  :: mmneq=300
      real*8 :: smd,vel,inj_mass_rate,x,y,z,rho_p,end_x,end_y,end_z
      real*8 :: ys_p(mcomp),cd_coef
      real*8 :: width,fanangle,sprayspread,angle
      integer:: inj_start,inj_num,subcycling,slit_dir
      integer:: inj_steps,inj_inter,inj_MAX
      integer:: nin=0,ierr1=0
      real*8 :: alpha,beta
      real*8 :: liquidmu,stens,cf,ck,cd,cb,cv,k
      real   ,parameter :: undef=-huge(1.)
      integer,parameter :: iundef=-huge(1)
      character(len=4)  :: breakup
      character(len=64) :: rsfile
      integer:: ios=0,icom,nb,NO
      character(4),parameter :: manual='manu',
     &                          slit='slit',
     &                          hole='hole',
     &                          coal_i='coal',
     &                          ilet='ilet',
     &                          usri='user'
      character(4),parameter :: 
     &                          injlst(6)=(/
     &                          manual,
     &                          slit,
     &                          hole,
     &                          coal_i,
     &                          ilet,usri
     &                          /)
      
!
! --- [namelist] !coal 
!
      character(len=6) :: fuel_type
      character(len=7) :: const_type
      character(len=4) :: nozzletype
      character(len=4) :: HeatExch_model
      character(len=4) :: movement_model
      character(len=4) :: drop_distribution
!
      real*8 :: a_RR
!
      real*8 :: Q_factor,HR_char,
     &          MF_H2O_proximate,
     &          MF_vola_proximate,
     &          MF_fixC_proximate,
     &          MF_ash_proximate,
     &          COALCV,Latent_H_CHO,Latent_H_H2O
      real*8 :: MF_C_ultimate,MF_H_ultimate,MF_O_ultimate,
     &          MF_N_ultimate,MF_S_ultimate
      real*8 :: H_CO2,H_H2O,H_O2,kappa
!
      character(6),parameter :: coal  ='coal  ',
     &                          liquid='liquid',
     &                          glass ='glass ',
     &                          solid ='solid '
      character(6),parameter :: fullst(4)=(/coal,liquid,glass,solid/) 
!
      character(7),parameter :: dens='density',
     &                          volu='volume '
      character(7),parameter :: cotlst(2)=(/dens,volu/)
!
      character(4),parameter :: evaF='evaF'
      character(4),parameter :: AbSi='AbSi'
      character(4),parameter :: LaKn='LaKn'
      character(4),parameter :: evaU='evaU'
      character(4),parameter :: evalst(4)=(/evaF,AbSi,LaKn,evaU/)

      character(4),parameter :: hetD='Dflt'
      character(4),parameter :: hetMass='Mass'
      character(4),parameter :: hetLaKn='LaKn'
      character(4),parameter :: hetU='User'
      character(4),parameter :: hetlst(4)=(/hetD,hetMass,hetLaKn,hetU/)

      character(4),parameter :: movD='Dflt'
      character(4),parameter :: movMilr='Milr'
      character(4),parameter :: movSomm='Somm'
      character(4),parameter :: movU='User'
      character(4),parameter :: movMrsi='Mrsi'
      character(4),parameter :: movSub='SubM'
      character(4),parameter :: movlst(6)=(/movD,movMilr,
     &                          movSomm,movU,movMrsi,movSub/)
!
      character(4),parameter :: cnst='cnst'
      character(4),parameter :: RosR='RosR'  !Rosin-Rammler distribution
      character(4),parameter :: noml='noml'  !normal distribution
      character(4),parameter :: nuki='Nuki'
      character(4),parameter :: kuro='kuro'
      character(4),parameter :: drplst(5)=(/cnst,RosR,noml,nuki,kuro/)
!
      integer :: fkd,s
      real*8 :: cp_p,T_ini_p,rad_emissivity,heat_fraction,T_ref_p
      real*8 :: Av_CHO,Ev_CHO,Av_H2O,Ev_H2O
!
      real*8 :: Av_HCN,Ev_HCN
      integer:: cal_evap_coal=0,
     &          cal_evap,
     &          q,rj,output_MF

      
!--------------
! --- local 
!--------------
      real*8 :: M_C,M_H,M_O,M_N
      integer :: bc_no=0,I,iq,jq
!----------------------------
! --- local evapration
!----------------------------
      real*8 :: Latent_heat(mcomp),
     &          P_sat(mcomp),
     &          T_vap(mcomp),
     &          T_boil(mcomp)
      integer ::          evap_ys(mcomp)
      character(len=4) :: evap_model(mcomp)
      real*8  :: P_0(mcomp),T_0(mcomp)
      integer :: cal_reaction
!
      integer :: particle_reac_no(mmneq)=iundef
!----------------------------------------- 
! --- Namelists for Particle Settings 
!----------------------------------------- 
      namelist /injector_number/
     &          Num_injector,total_start
!
      namelist /particle_model/ 
     &         smd,vel,inj_mass_rate,!inj_vol_mass,
     &         x,y,z,
     &         end_x,end_y,end_z,
     &         ys_p,
     &         inj_start,inj_steps,inj_num,inj_inter,inj_MAX,
     &         alpha,beta,rho_p,
     &         subcycling,
     &         breakup,
     &         nozzletype,width,fanangle,
     +         sprayspread,angle,slit_dir,
     %         liquidmu,stens,cf,ck,cd,cb,cv,k,
     &         cd_coef,
!
     &         fuel_type,const_type,
     &         evap_model,HeatExch_model,movement_model,
     &         drop_distribution,a_RR,
!
     &         Q_factor,COALCV,
     &         HR_char,Latent_H_CHO,Latent_H_H2O,
     &         MF_H2O_proximate,MF_vola_proximate,
     &         MF_fixC_proximate,MF_ash_proximate,
     &         MF_C_ultimate,MF_H_ultimate,MF_O_ultimate,
     &         MF_N_ultimate,MF_S_ultimate,
     &         H_CO2,H_H2O,H_O2,
!
     &         cp_p,T_ini_p,rad_emissivity,heat_fraction,T_ref_p,
!
     &         Av_CHO,Ev_CHO,Av_H2O,Ev_H2O,Av_HCN,Ev_HCN,
!
     &         cal_evap_coal,
     &         cal_evap,
!
     &         bc_no,
     &         output_MF,
!
     &         Latent_heat,P_sat,T_vap,T_boil,evap_ys,
     &         P_0,T_0,
!
     &         cal_reaction,particle_reac_no

!
!
!-----------------------
!     local  entities   
!-----------------------
      REAL*8 :: dum1,dum2,dum3,dl,t
!
      if(iprtcle==0) return
!
!-----------------------------
! ---- read particle_model
!-----------------------------
!
      Num_injector=iundef
      total_start=iundef
!
      rewind ifli
      read(ifli,injector_number,iostat=ios,err=200)
      if(Num_injector==iundef) then
        call FFRABORT(1,'ERR: Set &injector_number/[Num_injector]')
      endif
      if(total_start==iundef) then
        write(ifle,'(1X,a)') 'ERR: Set total injection time step'
        call FFRABORT(1,'ERR: Set &injector_number/[total_start]')
      endif
      if(ios<0) call FFRABORT(1,'ERR: &injector_number reading error')
!      
!
      rewind ifli
      N_injctr=0
      do I=1,Num_injector
        read(ifli,particle_model,iostat=ios,err=100)
        if(ios<0) exit
        call nml_errmsg0(ifle,ios,'&particle_model',ierr1)
        if(ierr1/=0) 
     &  call FFRABORT(1,'ERR: &particle_model reading error')
        N_injctr=N_injctr+1
      enddo
      iinj=N_injctr
!---------------------------
! --- initial reaction arry
!---------------------------
      call particle_chem
     &          (1,0,N_injctr,ifli,ifll,ifle,iprtcle,my_rank,ncomp,
     &           cntlnam,lcomp,igT_iter,ireac,ierror)
!
      if(Num_injector>N_injctr) then
        write(ifle,'(2X,2a)')
     &  'ERR: Num_injector in &injector_number',
     &  '> number of &particle_model'
        call FFRABORT
     & (1,'ERR: Redefining Num_injector in &injector_number')
      endif
!
      allocate(lreacP(nneqP),stat=ierror)
!-----------------------------
! --- allocate arry
!-----------------------------
      allocate(
     &         latent(ncomp,N_injctr),
     &         EVAP_COMP(ncomp,N_injctr), 
     &         PP_sat(ncomp,N_injctr), 
     &         TT_vap(ncomp,N_injctr), 
     &         TT_boil(ncomp,N_injctr), 
     &         TT_0(ncomp,N_injctr), 
     &         PP_0(ncomp,N_injctr), 
     &         stat=ierr1)
      if(ierr1/=0) 
     &   call FFRABORT(1,'ERR: &particle_model allcating latent(:,:)')
!-----------------------------
! --- allocate arry
!-----------------------------
      allocate(
     6         NINALL(0:N_injctr),NPALL(N_injctr),NINEND(0:N_injctr),
     &         SPHEIGHT(0:N_injctr),VOLUMEPC(0:N_injctr),
     &         dUVWNEW(N_injctr),
     &         particle_comp(N_injctr,ncomp),
     &         dSMD(N_injctr),
     &         dInjVolume(N_injctr),
     &         RHOP0(N_injctr),
     &         mass_p_total(N_injctr),
     &         lack_p(N_injctr),
     &         re_inj(N_injctr),
     &         iSprayInjStart(N_injctr),
     &         iSprayInjInter(N_injctr),
     &         iSprayInjEnd(N_injctr),
     &         iSprayInjTime(N_injctr),
     &         starttime(N_injctr),
     &         iSprayInjNum(N_injctr),
     &         iSprayInjTotal(N_injctr),
     &         INJ_FLAG(N_injctr),
     &         dXP0(N_injctr),
     &         dYP0(N_injctr),
     &         dZP0(N_injctr),
     &         endx(N_injctr),
     &         endy(N_injctr),
     &         endz(N_injctr),
     &         axis_inj(3,N_injctr),
     &         cBreakupModel(N_injctr),
     &         lBreakupNONE(N_injctr),
     &         lBreakupPILCH(N_injctr),
     &         lBreakupITAB(N_injctr),
     &         lBreakupNTAB(N_injctr),
     &         dAlpha(N_injctr),
     &         dBeta(N_injctr),
     &         cNozzleType(N_injctr),
     &         dSPRAYANGLE(N_injctr),
     &         dSPANGLEMAG(N_injctr),
     &         dFANANGLE(N_injctr),
     &         dSLITWIDTH(N_injctr),
     &         dSLITdir(N_injctr),
     &         dLiquidMu(N_injctr),
     &         dLiquidStens(N_injctr),
     &         dCf(N_injctr),
     &         dCk(N_injctr),
     &         dCd(N_injctr),
     &         ddCd(N_injctr),
     &         dCb(N_injctr),
     &         dCv(N_injctr),
     &         dK(N_injctr),
     &         injkd(N_injctr),
     &         evakd(N_injctr,ncomp),
     &         heatexkd(N_injctr),
     &         movkd(N_injctr),
     &         dropkd(N_injctr),
     &         dRR_a(N_injctr),
     &         nu_intnst(N_injctr),rain_vel(N_injctr),
     &         Rh_d(N_injctr),func_d(N_injctr),pcl_WDR(N_injctr),
     &         stat=ierr1)
      if(ierr1/=0)
     &   call FFRABORT(1,'ERR: &particle_model allcating error-0')
!-------------------------------
! ---
!-------------------------------
      NINALL=0
      NPALL=0
      NINEND=0
      SPHEIGHT=0.d0
      VOLUMEPC=0.d0
      dUVWNEW=0.d0
      particle_comp=0.d0
      dSMD=0.d0
      dInjVolume=0.d0
      RHOP0=0.d0
      iSprayInjStart=0.d0
      iSprayInjInter=0.d0
      iSprayInjEnd=0.d0
      iSprayInjTime=0.d0
      starttime=0.d0
      iSprayInjStart=0
      iSprayInjInter=0
      iSprayInjEnd=0
      iSprayInjTime=0
      iSprayInjNum=0
      iSprayInjTotal=0
!
      injkd(:)=imanual
      evakd(:,:)=ievaN
      heatexkd(:)=ihetD
      movkd(:)=imovD
      dropkd(:)=iRosR
      dRR_a(:)=3.d0
!      Latent_L(:)=0.d0
      Latent(:,:)=0.d0
      PP_0(:,:)=0.d0
      TT_0(:,:)=0.d0
      EVAP_COMP(:,:)=0.d0
      PP_sat(:,:)=0.d0
      TT_vap(:,:)=0.d0
      TT_boil(:,:)=0.d0
!      
      pcl_WDR(:)=1.d0
!
      allocate(ifuel(N_injctr),
     &         itype(N_injctr),
     &         Q_fac(N_injctr),
     &         a_C(N_injctr),
     &         b_H(N_injctr),
     &         c_O(N_injctr),
     &         MF_H2O_p(N_injctr),
     &         MF_vola_p(N_injctr),
     &         MF_fixC_p(N_injctr),
     &         MF_ash_p(N_injctr),
     &         MF_C_e(N_injctr),
     &         MF_H_e(N_injctr),
     &         MF_O_e(N_injctr),
     &         MF_N_e(N_injctr),
     &         MF_S_e(N_injctr),
     &         Q_MF_vola(N_injctr),
     &         Q_MF_fixC(N_injctr),
     &         Q_MF_H2O(N_injctr),
     &         net_MF_HCN(N_injctr),
     &         net_MF_Nchar(N_injctr),
     &         Q_MF_ASH(N_injctr),

     &         DH_CHO(N_injctr),
     &         HR_cha(N_injctr),
     &         LH_CHO(N_injctr),
     &         LH_H2O(N_injctr),

     &         M_CHO(N_injctr),
     &         H_CHO(N_injctr),

     &         xx_O2(N_injctr),
     &         H_CO2a(N_injctr),
     &         H_H2Oa(N_injctr),
     &         H_O2a(N_injctr),

     &         eps_rad(N_injctr),

     &         Tini_p(N_injctr),
     &         Tref_p(N_injctr),
     &         Cp_pa(N_injctr),

     &         heat_f(N_injctr),

     &         Avap_CHO(N_injctr),
     &         Evap_CHO(N_injctr),
     &         Avap_H2O(N_injctr),
     &         Evap_H2O(N_injctr),
     &         Avap_HCN(N_injctr),
     &         Evap_HCN(N_injctr),

!     &         ical_CHAR_COMB(N_injctr),
     &         ical_la_coal(N_injctr),
     &         ical_evap(N_injctr),
!
     &         NO_BC(N_injctr),
     &         bc_nb(N_injctr),
     &         P_MF(N_injctr),
!
     &         partreac(N_injctr),INJ_RNUM(N_injctr),
     &         IDX_PRTRAC(N_injctr,max(1,nneqP)),
!
     &         stat=ierr1)

      if(ierr1/=0) 
     &   call FFRABORT(1,'ERR: &particle_model allcating error-1')
!
!
      ddcd(:)=1.d0
      Tref_p(:)=285.15
!------------------------------
! --- read particle model 
!------------------------------
      ierror=0
      nin=0
      NO_BC(:)=0
      bc_nb(:)=0
!
      partreac(:)=0
      INJ_RNUM(:)=0
      IDX_PRTRAC(:,:)=0
!-------------------------------
! --- read particle model
!-------------------------------
      rewind ifli
      do I=1,N_injctr       !Num_injector
        nin=nin+1
        ys_p(:)=undef
        smd=undef
        vel=undef
!
        rho_p=undef
        cd_coef=1.d0        
!
        inj_mass_rate=undef
        x=undef
        y=undef
        z=undef
        inj_start=iundef
        inj_steps=iundef
        inj_num=iundef
        inj_inter=1
        inj_MAX=0
        breakup=''
        nozzletype=''
        fuel_type=''
        const_type=''
        HeatExch_model=''        
        movement_model=''
        drop_distribution=''
        a_RR=undef
!        
        width=undef       
        fanangle=undef    
        sprayspread=undef 
        angle=undef
        slit_dir=iundef
!
        subcycling=iundef
!
        liquidmu=undef
        stens=undef
!
        Q_factor =undef
        HR_char  =undef
        MF_H2O_proximate=undef
        MF_vola_proximate=undef
        MF_fixC_proximate=undef
        MF_ash_proximate=undef
        MF_C_ultimate=undef
        MF_H_ultimate=undef
        MF_O_ultimate=undef
        MF_N_ultimate=undef
        MF_S_ultimate=undef
        cp_p=undef
        T_ini_p=300.d0
        rad_emissivity=undef
        heat_fraction=undef
        Av_CHO=undef
        Ev_CHO=undef
        Av_H2O=undef
        Ev_H2O=undef
        Av_HCN=undef
        Ev_HCN=undef
!        
!        cal_CHAR_COMB=0
        cal_evap_coal=0
        cal_evap=0
!
        bc_no=0
        output_MF=0
!
        Latent_heat(:)=undef
        P_sat(:)=undef
        T_vap(:)=undef
        T_boil(:)=undef
        evap_ys(:)=iundef
        evap_model(:)=''    
        P_0(:)=undef
        T_0(:)=undef
!
        cal_reaction=0
!
!        parcel_WDR=undef
!
        read(ifli,particle_model,iostat=ios,err=100)
        if(ios<0) exit
!
        if(cal_reaction>0) then
          if(nneqP==0) then
            call FFRABORT(1,'ERR: NO &particle_chemreac is define')
          endif
        endif
!----------------------------------
        if(ical_P_R==1) then
          if(cal_reaction<0) then
            write(ifle,*) 
     &     'ERR: [cal_reaction>0 or =0] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          else
            partreac(NIN)=cal_reaction
!
            if(cal_reaction>0) then
              do iq=1,nneqP
              do jq=1,mmneq
              if(particle_reac_no(jq)==iundef) cycle
              if(particle_reac_no(jq)>nneqP.or.particle_reac_no(jq)<=0) 
     &        then
                write(ifle,*) ' ### ERR: No such reaction defined'
                call FFRABORT(1,'&particle_chemreac/particle_reac_no')
              endif
              if(iq==particle_reac_no(jq)) then
                INJ_RNUM(NIN)=INJ_RNUM(NIN)+1  
                IDX_PRTRAC(NIN,INJ_RNUM(NIN))=iq    

              endif
              enddo
              enddo
            endif
          endif
        else
          if(cal_reaction>0) then
            write(ifle,*) 
     &     'ERR: cal_reaction>0 in &particle_model'
            write(ifle,*) 
     &     'MSG: NO particle reaction is set in &chemreac'
            call FFRABORT(1,'namelist/particle_model')
          endif
        endif
!
!        if(nneqP/=0.and.cal_reaction/=0)then
!          call particle_chem
!     &    (1,0,N_injctr,ifli,ifll,ifle,iprtcle,my_rank,ncomp,
!     &     cntlnam,lcomp,igT_iter,ireac,ierror)
!        endif
!----------------------------------
! --- 
!
        if(nozzletype=='') then
          write(ifle,*) 
     &     'ERR: lack of [nozzletype] in &particle_model'
          write(ifle,*) 
     &     'MSG: nozzletype="manu","slit","hole", "coal","ilet"'
          call FFRABORT(1,'namelist/particle_model')
        else
           fkd=0
           call nml_listno(6,injlst,nozzletype,fkd)
           call nml_chkchr0(ifle,'nozzletype',nozzletype,fkd,ierr1)
           if(ierr1/=0)
     &       call FFRABORT(1,'ERR:nozzletype is NOT SURPPORTED') 

           cNozzleType(nin)=nozzletype
           injkd(nin)=fkd
!
           if((cNozzleType(nin)/='slit')
     &     .and.(cNozzleType(nin)/='hole')
     &     .and.(cNozzleType(nin)/='manu')
     &     .and.(cNozzleType(nin)/='coal')
     &     .and.(cNozzleType(nin)/='ilet')
     &     .and.(cNozzleType(nin)/='user')
     &      ) then
             write(*,*) "ERR: NozzleType must be ",
     &                  '[manu],[slit],[hole],[coal],[ilet]' 
             write(*,*) "Resseted NozzleType to slit."
             call FFRABORT(1,'ERR:nozzletype is NOT SURPPORTED')
!
           endif
           if(injkd(nin)==imanual) then
           elseif(injkd(nin)==islit) then
           elseif(injkd(nin)==ihole) then
           elseif(injkd(nin)==i_coal) then
           elseif(injkd(nin)==inlet) then
           endif
        endif
!
        if(drop_distribution=='') then
          write(ifle,'(1x,a)') 
     &   'lack of data : [drop_distribution] in &particle_model'
          write(ifle,'(1x,a)') 
     &     'MSG: drop_distribution="cnst","RosR" or "noml"'
          write(ifle,'(1x,3a)') 
     &     'MSG: drop_distribution="cnst": constant DIA',
     &     '     drop_distribution="RosR": Rosin-Rammler',
     &     '     drop_distribution="noml": normal'
          call FFRABORT(1,'namelist/particle_model')
        else
          fkd=0
          call nml_listno(5,drplst,drop_distribution,fkd)
          call nml_chkchr0
     &         (ifle,'drop_distribution',drop_distribution,fkd,ierr1)
          if(ierr1/=0) then
              write(ifle,'(1x,a)') 
     &        'ERR: [drop_distribution] is error'
              write(ifle,'(1x,4a)') 
     &      'MSG: drop_distribution="cnst": constant DIA',
     &      '     drop_distribution="RosR": Rosin-Rammler',
     &      '     drop_distribution="noml": normal',
     &      '     drop_distribution="Nuki": Nukiyama-Tanasawa'
              call FFRABORT(1,'ERR: drop_distribution is error')
          endif
          dropkd(nin)=fkd
          if(fkd==iconst) then
            write(ifll,'(1x,a)') 'MSG: drop-size model: constant'
          elseif(fkd==iRosR) then
            if(a_RR==undef) then
              call FFRABORT(1,'ERR: NOT define a_RR for Rosin-Rammler')
            else
              dRR_a(NIN)=a_RR
            endif
            write(ifll,'(1x,a)') 'MSG: drop-size model: Rosin-Rammler'
          elseif(fkd==inoml) then
            write(ifll,'(1x,a)') 'MSG: drop-size model: normal'
          elseif(fkd==inuki) then
            write(ifll,'(1x,a)') 
     &             'MSG: drop-size model: Nukiyama-Tanasawa'
          endif
        endif
!
        if(fuel_type=='') then
          write(ifle,*) 'lack of data : [fuel_type] in &particle_model'
          write(ifle,*) 
     &     'MSG: fuel_type="coal","liquid","glass","solid"'
          call FFRABORT(1,'namelist/particle_model')
        else
          fkd=0
          call nml_listno(4,fullst,fuel_type,fkd)
          call nml_chkchr0(ifle,'fuel_type',fuel_type,fkd,ierr1)
          if(ierr1/=0) then
            write(ifle,*) 
     &     'MSG: fuel_type="coal","liquid"'
            call FFRABORT(1,'ERR:fuel_type is NOT SURPPORTED')
          endif
          ifuel(nin)=fkd
        endif
!
        if(iprtcle==2) then
          if(ifuel(nin)/=iliquid) then
            write(ifle,'(1x,a)') 
     &      'ERR-MSG: fuel_type="liquid" for rain_WDR=1'
            write(ifle,'(1x,a,I6)') 
     &      'MSG: error in particle_model no=',nin
            call FFRABORT(1,'ERR:fuel_type error')
          endif
        endif

!
        if(const_type=='') then
          write(ifle,*) 'lack of data: [const_type] in &particle_model'
          write(ifle,*) 
     &     'MSG: const_type="density","volume"'
          call FFRABORT(1,'namelist/particle_model')
        else
          fkd=0
          call nml_listno(2,cotlst,const_type,fkd)
          call nml_chkchr0(ifle,'const_type',const_type,fkd,ierr1)
          if(ierr1/=0) then
            write(ifle,*) 
     &       'MSG: const_type="density","volume"'
            call FFRABORT(1,'ERR:const_type is error')
          endif
          itype(nin)=fkd
        endif
!
!
        if(fuel_type=='liquid') then
          do icom=1,NCOMP
          if(evap_model(icom)=='') then
            write(ifll,*) 
     &       'MSG: evaporation model for species: ',
     &        ICOM ,' NOT be defined'
            evakd(nin,ICOM)=ievaN
          else
            iflag_eva=1
            fkd=0
            call nml_listno(4,evalst,evap_model(ICOM),fkd)
            call nml_chkchr0
     &       (ifle,'evap_model',evap_model(icom),fkd,ierr1)
            if(ierr1/=0) then
              write(ifle,*) 
     &        'ERR: evap_model(?) ?= ',ICOM, 
     &            'is NOT supported','INI= ',NIN
              write(ifll,*) 
     &        'MSG: evap_model=',evalst(:),'are supported'
              call FFRABORT(1,'ERR: evap_model is error')
            endif
            evakd(nin,ICOM)=fkd
          endif
          enddo
!
          if(cal_evap<0) then
            call FFRABORT(1,'ERR:cal_evap<0')
          else
            ical_evap(nin)=cal_evap
          endif
!
          DO ICOM=1,NCOMP
          if(evap_ys(ICOM)==iundef) then
            write(ifle,*) 
     &      'ERR: evap_ys(ICOM),NOT defined, ICOM=',
     &                       ICOM
            write(ifle,*) 'INIJECTOR No= ',NIN
            write(ifle,*) 'Evaporation comp define '
            call FFRABORT(1,'STOP at namelist/particle_model')
          else
            if(evap_ys(ICOM)==1) then
              EVAP_COMP(ICOM,nin)=dble(evap_ys(ICOM))
              if(evakd(nin,ICOM)==ievaN)then
                write(ifle,*) 
     &         'ERR: evap_model(?) ICOM=? ',
     &         ICOM, 'NOT been defined','INI= ',NIN
                write(ifle,*) 'MSG: evap_ys(ICOM),ICOM=',
     &                       evap_ys(ICOM),ICOM
                write(ifle,*) 'INIJECTOR No= ',NIN
                call FFRABORT
     &          (1,'STOP at namelist/particle_model')
              endif
            elseif(evap_ys(ICOM)/=0) then
              write(ifle,*) 'ERR: evap_ys(ICOM),ICOM=',
     &                       evap_ys(ICOM),ICOM
              write(ifle,*) 'MSG: evap_ys(ICOM)=0 or 1'
              write(ifle,*) 'INIJECTOR No= ',NIN
              call FFRABORT(1,'namelist/particle_model')
            endif
          endif
!
          if(evap_ys(ICOM)/=1) cycle
!
          if(Latent_heat(ICOM)==undef) then
            write(ifle,*) 
     6   'ERR: [Latent_heat(?)] in &particle_model,?=',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
            Latent(ICOM,nin)=Latent_heat(ICOM)
          endif

          if(P_0(ICOM)==undef) then
            write(ifle,*) 
     6   'ERR: lack of [P_0(?)] in &particle_model,?=',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
             PP_0(ICOM,nin)=P_0(ICOM)
          endif

          if(T_0(ICOM)==undef) then
            write(ifle,*) 
     6   'ERR: lack of [T_0(?)] in &particle_model,?=',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
             TT_0(ICOM,nin)=T_0(ICOM)
          endif
! --- 
          if(T_vap(ICOM)==undef) then
            write(ifle,*) 
     6   'ERR: lack of [T_vap(?)] in &particle_model,?=',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
            TT_vap(ICOM,nin)=T_vap(ICOM)
          endif

          if(T_boil(ICOM)==undef) then
            write(ifle,*) 
     6   'ERR: lack of [T_boil(?)] in &particle_model,?=',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
            TT_boil(ICOM,nin)=T_boil(ICOM)
          endif
! --- 
          enddo
        endif
! --- 
        if(HeatExch_model=='') then
          write(ifle,*) 
     &       'ERR: lack of data : [HeatExch_model] in &particle_model'
          write(ifle,*) 
     &       'MSG: HeatExch_model="????":'
          call FFRABORT(1,'namelist/particle_model')
        else
          fkd=0
          call nml_listno(4,hetlst,HeatExch_model,fkd)
          call nml_chkchr0(ifle,'HeatExch_model',HeatExch_model,fkd,ierr1)
          if(ierr1/=0) then
            write(ifle,*) 
     &      'MSG: HeatExch_model="????","..."'
            call FFRABORT(1,'ERR: HeatExch_model is error')
          endif
          heatexkd(nin)=fkd
        endif
! ---- 
        if(Movement_model=='') then
          write(ifle,*) 
     &       'ERR: lack of data : [Movement_model] in &particle_model'
          write(ifle,*) 
     &       'MSG: Movement_model="????":'
          call FFRABORT(1,'namelist/particle_model')
        else
          fkd=0
          call nml_listno(6,movlst,movement_model,fkd)
          call nml_chkchr0(ifle,'movement_model',movement_model,fkd,ierr1)
          if(ierr1/=0) then
            write(ifle,*) 
     &      'MSG: movement_model="????","..."'
            call FFRABORT(1,'ERR: movement_model is error')
          endif
          movkd(nin)=fkd
        endif
! ---- 
        if(breakup=='') then
          write(ifle,*) 'lack of data : [breakup] in &particle_model'
          write(ifle,*) 'MSG: breakup="NONE","PILC","NTAB","ITAB"'
          write(ifle,*) 'MSG: breakup="NONE" for solid particle'
          call FFRABORT(1,'namelist/particle_model')
        else
          
        endif
!
        if(smd==undef) then
          write(ifle,*) 'lack of data : [smd] in &particle_model'
          call FFRABORT(1,'namelist/particle_model')
        endif

        if(vel==undef) then
          write(ifle,*) 'lack of data : [vel] in &particle_model'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(rho_p==undef) then
          write(ifle,*) 'lack of data : [rho_p] in &particle_model'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(cd_coef==undef) then
          write(ifle,*) 'lack of data : [cd_coef] in &particle_model'
          call FFRABORT(1,'namelist/particle_model')
        elseif(cd_coef>1.d0) then
          call FFRABORT(1,'ERR: New FFR version : cd_coef<=1.d0')
        endif
!
        if(inj_mass_rate==undef) then
          write(ifle,'(1X,2a)') 
     & 'ERR: lack of data : [inj_mass_rate] in &particle_model'
          call FFRABORT(1,'namelist/particle_model')
        endif

        if(x==undef.or.y==undef.or.z==undef) then
          write(ifle,*) 'lack of data : [x,y,z] in &particle_model'
          call FFRABORT(1,'namelist/particle_model')
        endif

        if(inj_start==iundef) then
          write(ifle,*) 'lack of data : [inj_start] in &particle_model'
          write(ifle,*) 'MSG: injection start time step'
          call FFRABORT(1,'namelist/particle_model')
        endif

        if(inj_steps==iundef) then
          write(ifle,*) 'lack of data : [inj_steps] in &particle_model'
          write(ifle,*) 'MSG: injection injection time-step number'
          call FFRABORT(1,'namelist/particle_model')
        endif

        if(inj_num==iundef) then
          write(ifle,*) 'lack of data : [inj_num] in &particle_model'
          write(ifle,*) 'MSG: particle number in every time step'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(subcycling==iundef) then
          write(ifle,*) 
     &   'lack of data : [subcycling] in &particle_model'
          write(ifle,*) 
     &   'MSG: particle subcycling number in every time step'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(width==undef) then
          write(ifle,*) 'lack of data : [width] in &particle_model' 
          write(ifle,*) 'MSG: [width] is Width when slit-injector'
          write(ifle,*) 'MSG: [width] is ',
     &        'Diameter when hole-injector or manual'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(slit_dir==iundef) then
          write(ifle,*) 'lack of data : [slit_dir] in &particle_model' 
          write(ifle,*) 'MSG: [slit_dir] is slit-Width direction'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(fanangle==undef) then
          write(ifle,*) 'lack of data : [fanangle] in &particle_model'
          write(ifle,*) 'MSG: fanangle '
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(sprayspread==undef.and.nozzletype=='slit') then
          write(ifle,*) 
     &  'lack of data : [sprayspread] in &particle_model'
          write(ifle,*) 'MSG: sprayspread only for [slit] injector'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
        if(angle==undef) then
          write(ifle,*) 'lack of data : [angle] in &particle_model'
          write(ifle,*) 'MSG: angle injection direction angle fron Y'
          call FFRABORT(1,'namelist/particle_model')
        endif
!
!
        dUVWNEW(nin)=vel
        particle_comp(nin,:)=0.d0
        dSMD(nin)=smd
        ddcd(nin)=cd_coef
        RHOP0(nin)=rho_p
        dAlpha(nin)=alpha
        dBeta(nin)=beta
!                           :inj_mass_rate=(kg/s)
!                            
!

        if(fuel_type=='coal') then
          dInjVolume(nin)=inj_mass_rate     !(kg/s)
        elseif(fuel_type=='liquid') then
          dInjVolume(nin)=inj_mass_rate     !(kg/s)
        else
          dInjVolume(nin)=inj_mass_rate     !(kg/s)
        endif
!
        dXP0(nin) =x
        dYP0(nin) =y
        dZP0(nin) =z
        endx(nin)=end_x
        endy(nin)=end_y
        endz(nin)=end_z
        dum1=end_x-x
        dum2=end_y-y
        dum3=end_z-z
        dl=dsqrt(dum1**2+dum2**2+dum3**2)
        axis_inj(1,nin)=dum1/dl
        axis_inj(2,nin)=dum2/dl
        axis_inj(3,nin)=dum3/dl
!

        INJ_FLAG(N_injctr)=.false.
        iSprayInjStart(nin)=inj_start
        iSprayInjTime(nin)=inj_steps
        iSprayInjNum(nin)=inj_num
        iSprayInjTotal(nin)=inj_MAX
        iSprayInjInter(nin)=inj_inter
        if(iSprayInjStart(nin).le.0) then
          write(*,*) "WRN: inj_start is less than 1."
          write(*,*) "MSG: Reseted inj_start to 1."
          iSprayInjStart(nin)=1
        endif
        starttime(nin)=0.d0
        iSprayInjEnd(nin)=iSprayInjStart(nin)+iSprayInjTime(nin)-1
!
        if(iSprayInjNum(nin).lt.0) then
          write(*,*) "WRN: inj_num must be positive value."
          write(*,*) "MSG: Reseted inj_num to zero."
          iSprayInjNum(nin)=0
        endif
!
        
        dSLITWIDTH(nin)=width
        dSLITdir(nin)=dble(slit_dir)
        dSPRAYANGLE(nin)=angle
        dSPANGLEMAG(nin)=sprayspread
        dFANANGLE(nin)=fanangle
!
        ITER_INNER=subcycling
!
        cBreakupModel(nin)=breakup
        if(cBreakupModel(nin).eq."NONE") THEN
          lBreakupNONE(nin)= .TRUE.
          lBreakupPILCH(nin)= .FALSE.
          lBreakupNTAB(nin)= .FALSE.
          lBreakupITAB(nin)= .FALSE.
        else if(cBreakupModel(nin).eq."PILC") THEN
          lBreakupNONE(nin)= .FALSE.
          lBreakupPILCH(nin)= .TRUE.
          lBreakupNTAB(nin)= .FALSE.
          lBreakupITAB(nin)= .FALSE.
        else if(cBreakupModel(nin).eq."NTAB") THEN
          lBreakupNONE(nin)= .FALSE.
          lBreakupPILCH(nin)= .FALSE.
          lBreakupNTAB(nin)= .TRUE.
          lBreakupITAB(nin)= .FALSE.
        else if(cBreakupModel(nin).eq."ITAB") THEN
          lBreakupNONE(nin)= .FALSE.
          lBreakupPILCH(nin)= .FALSE.
          lBreakupNTAB(nin)= .FALSE.
          lBreakupITAB(nin)= .TRUE.
        else
          write(*,*) 
     &   "Warning!!Breakup Model must be NONE,NTAB,PILCH."
          write(*,*) "Breakup Model will be set to NONE."
          cBreakupModel(nin)="NONE"
          lBreakupNONE(nin)=.TRUE.
          lBreakupPILCH(nin)=.FALSE.
          lBreakupNTAB(nin)=.FALSE.
          lBreakupITAB(nin)=.FALSE.
        endif
!
!        lPTRESTART=.FALSE.
!        if(restart.eq.'yes') lPTRESTART=.TRUE.
!
        if(.not.lBreakupNONE(nin)) then
          dLiquidMu(nin)=0.77D-03
          dLiquidStens(nin)=0.024d0
          dCf(nin)=0.3333333333333333333333
          dCk(nin)=8.0
          dCd(nin)=5.0
          dCb(nin)=0.5
          dCv(nin)=1.0
          dK(nin)=3.33333333333333333333333
!
!
!
          dLiquidMu(nin)= liquidmu
          dLiquidStens(nin)=stens
          dCf(nin)=cf
          dCk(nin)=ck
          dCb(nin)=cb
          dCd(nin)=cd
          dCv(nin)=cv
          dK(nin)=k
        endif
!
        if(iprtcle==2) then
          if(injkd(nin)/=inlet) then
            write(ifle,'(1x,a)')
     &       'ERR: &/particle_model/nozzletype.ne."ilet"'
            call FFRABORT(1,'MSG: rain model need nozzletype="ilet"')
          endif
!
!          if(parcel_WDR==undef) then
!            write(ifle,'(1x,a)')
!     &       'ERR: &/particle_model/parcel_WDR for rain_WDR'
!            call FFRABORT(1,'MSG: rain model need parcel_WDR')
!          endif
!          pcl_WDR(NIN)=parcel_WDR
        endif
!
        if(ifuel(nin)==icoal) then
          IF(Q_factor==undef) THEN
            write(ifle,*) 
     &      'lack of data : [Q_factor] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(HR_char==undef) THEN
            write(ifle,*) 
     &      'lack of data : [HR_char] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_H2O_proximate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_H2O_proximate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_vola_proximate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_vola_proximate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_fixC_proximate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_fixC_proximate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_ash_proximate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_ash_proximate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_C_ultimate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_C_ultimate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_H_ultimate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_H_ultimate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_O_ultimate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_O_ultimate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_N_ultimate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_N_ultimate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
          IF(MF_S_ultimate==undef) THEN
            write(ifle,*) 
     &      'lack of data : [MF_S_ultimate] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          ENDIF
!
          IF(nin==1) THEN
          do S=1,NCOMP
            if(trim(spinam(s))=='Ca_Hb_Oc') then
              II_CHO=s
            endif
            if(trim(spinam(s))=='O2') then
              II_O2=s
            endif
            if(trim(spinam(s))=='N2') then
              II_N2=s
            endif
            if(trim(spinam(s))=='CO2') then
              II_CO2=s
            endif
            if(trim(spinam(s))=='CO') then
              II_CO=s
            endif
            if(trim(spinam(s))=='H2O') then
              II_H2O=s
            endif
            if(trim(spinam(s))=='NO') then
              II_NO=s
            endif
            if(trim(spinam(s))=='N') then
              II_N=s
            endif
            if(trim(spinam(s))=='O') then
              II_O=s
            endif
            if(trim(spinam(s))=='H') then
              II_H=s
            endif
            if(trim(spinam(s))=='OH'.or.trim(spinam(s))=='NH3') then
!            if(trim(spinam(s))=='NH3') then
              II_OH=s
            endif
            if(trim(spinam(s))=='HCN') then
              II_HCN=s
            endif
            if(trim(spinam(s))=='CN') then
              II_CN=s
            endif
            if(trim(spinam(s))=='C') then
              II_C=s
            endif
            if(trim(spinam(s))=='S') then
              II_S=s
            endif
            if(trim(spinam(s))=='ASH') then
              II_ASH=s
            endif
          enddo
          if(II_CHO==0) then
            call FFRABORT(1,'ERR:II_CHO=0)')
          endif
          if(II_O2==0) then
            call FFRABORT(1,'ERR:II_O2=0)')
          endif
          if(II_N2==0) then
            call FFRABORT(1,'ERR:II_N2=0)')
          endif
          if(II_CO2==0) then
            call FFRABORT(1,'ERR:II_CO2=0)')
          endif
          if(II_CO==0) then
            call FFRABORT(1,'ERR:II_CO=0)')
          endif
          if(II_H2O==0) then
            call FFRABORT(1,'ERR:II_H2O=0)')
          endif
          if(II_N==0) then
            call FFRABORT(1,'ERR:II_N=0)')
          endif
          if(II_O==0) then
            call FFRABORT(1,'ERR:II_O=0)')
          endif
          if(II_NO==0) then
            call FFRABORT(1,'ERR:II_NO=0)')
          endif
          if(II_H==0) then
            call FFRABORT(1,'ERR:II_H=0)')
          endif
          if(II_OH==0) then
            call FFRABORT(1,'ERR:II_OH=0)')
          endif
          if(II_HCN==0) then
            call FFRABORT(1,'ERR:II_HCN=0)')
          endif
          if(II_CN==0) then
            call FFRABORT(1,'ERR:II_CN=0)')
          endif
          if(II_C==0) then
            call FFRABORT(1,'ERR:II_C=0)')
          endif
          if(II_S==0) then
            call FFRABORT(1,'ERR:II_S=0)')
          endif
          if(II_ASH==0) then
            call FFRABORT(1,'ERR:II_ASH=0)')
          endif
          ENDIF
!
          Q_fac(nin)=Q_factor
          MF_H2O_p(nin)=MF_H2O_proximate
          MF_vola_p(nin)=MF_vola_proximate
          MF_fixC_p(nin)=MF_fixC_proximate
          MF_ash_p(nin)=MF_ash_proximate
          dum1=!MF_H2O_p(nin)
     &        +MF_vola_p(nin)
     &        +MF_fixC_p(nin)
     &        +MF_ash_p(nin)
          IF(abs(1.d0-dum1)>1.d-10) THEN 
            WRITE(IFLE,'(1X,a)') 'ERR: Proximate percent /= 100%'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:CHO= ',MF_vola_p(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:C= ',MF_fixC_p(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:Ash= ',MF_ash_p(nin)*100.d0,' %'         
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &        'MSG: Proximate percent = CHO+C+ASH =',
     &         dum1*100.d0,'% /=100%'
            CALL FFRABORT(1,'ERR: stop at module_particle')
          ENDIF
! --- WATER-CONTANTED PROXIMATE VALUE
          dum1=MF_H2O_proximate
     &        +MF_vola_proximate
     &        +MF_fixC_proximate
     &        +MF_ash_proximate
          if(abs(dum1)<1.d-10) then
            call FFRABORT(1,'ERR: sum of proximate value = 0%')
          endif
          MF_H2O_p(nin)=MF_H2O_proximate/dum1
          MF_vola_p(nin)=MF_vola_proximate/dum1
          MF_fixC_p(nin)=MF_fixC_proximate/dum1
          MF_ash_p(nin)=MF_ash_proximate/dum1
! 
          MF_C_e(nin)=MF_C_ultimate
          MF_H_e(nin)=MF_H_ultimate
          MF_O_e(nin)=MF_O_ultimate
          MF_N_e(nin)=MF_N_ultimate
          MF_S_e(nin)=MF_S_ultimate
          dum1=MF_C_e(nin)
     &        +MF_H_e(nin)
     &        +MF_O_e(nin)
     &        +MF_N_e(nin)
     &        +MF_S_e(nin)
     &        +MF_ash_p(nin)
!
          IF(abs(1.d0-dum1)>1.d-10) THEN 
            WRITE(IFLE,'(1X,a)') 'ERR: Ultimate percent /= 100%'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:C= ',MF_C_e(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:H= ',MF_H_e(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:O= ',MF_O_e(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:N= ',MF_N_e(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:S= ',MF_S_e(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &       'MSG:ash= ',MF_ash_p(nin)*100.d0,' %'
            WRITE(IFLE,'(1X,a,E15.7,a)') 
     &        'MSG:Ultimate percent =C+H+O+N+S+ASH= ',
     &        dum1*100.d0,'% /=100%'
!            CALL FFRABORT(1,'ERR: stop at module_particle')
          ENDIF
!
          Q_MF_vola(nin)=Q_fac(nin)*MF_vola_p(nin)
          write(ifll,'(1x,a,E15.7)') 
     &     'MSG: Q-H2O Mass fraction of CHO=',Q_MF_vola(nin)
          Q_MF_fixC(nin)=MF_fixC_p(nin)-(Q_MF_vola(nin)-MF_vola_p(nin))
          write(ifll,'(1x,a,E15.7)') 
     &     'MSG: Q-H2O Mass fraction of FixC=',Q_MF_fixC(nin)
          Q_MF_H2O(nin)=MF_H2O_p(nin)
          write(ifll,'(1x,a,E15.7)') 
     &     'MSG: Q-H2O Mass fraction of H2O=',Q_MF_H2O(nin)
          Q_MF_ASH(nin)=MF_ash_p(nin)
          write(ifll,'(1x,a,E15.7)') 
     &     'MSG: Q-H2O Mass fraction of ASH=',Q_MF_ASH(nin)
!
          kappa=Q_MF_vola(nin)/(Q_MF_vola(nin)+Q_MF_fixC(nin))
          dum1=(Q_MF_vola(nin)+Q_MF_fixC(nin))
          net_MF_HCN(nin)=dum1*MF_N_e(nin)*kappa
          write(ifll,'(1x,a,E15.7)') 
     &     'MSG: NET Mass fraction of HCN=',net_MF_HCN(nin)
          net_MF_Nchar(nin)=dum1*MF_N_e(nin)*(1.d0-kappa)!kappa
          write(ifll,'(1x,a,E15.7)') 
     &     'MSG: NET Mass fraction of charN=',net_MF_Nchar(nin)
!
          dum1=(Q_MF_vola(nin)/(Q_MF_vola(nin)+Q_MF_fixC(nin)))
          M_H=MF_H_e(nin)/dum1
          M_O=MF_O_e(nin)/dum1
!          M_N=MF_N_e(nin)/dum1
          M_C=1.d0-M_H-M_O!-M_N
          dum1=(M_C/12.d0+M_H+M_O/16.d0)
          a_C(nin)=(M_C/12.d0)/dum1
          b_H(nin)=M_H/dum1
          c_O(nin)=(M_O/16.d0)/dum1
!
          dum1=a_C(nin)+b_H(nin)+c_O(nin)
          IF(abs(1.d0-dum1)>1.d-10) THEN 
            WRITE(IFLE,'(1X,a,3E15.7)')
     &      'MSG: CaHbOc: a,b,c=',a_C(nin),b_H(nin),c_O(nin)
            WRITE(IFLE,'(1X,a,E15.7,a)')
     &      'MSG: a+b+c=',dum1,' /=1.d0'
            CALL FFRABORT(1,'ERR: stop at module_particle')
          endif
          dum1=a_C(nin)   
          a_C(nin)=a_C(nin)/dum1!a_C(nin)
          b_H(nin)=b_H(nin)/dum1!a_C(nin)
          c_O(nin)=c_O(nin)/dum1!a_C(nin)
          if(my_rank==0) then
            write(ifll,'(2X,a,3f10.3)') 'MSG: CaHbOc => a,b,c=',
     &                         a_C(nin),b_H(nin),c_O(nin)
          endif
!----------------------------
! --- 
!----------------------------
          dum1=1.d0-MF_H2O_p(nin)-MF_ash_p(nin)-Q_MF_vola(nin)
          HR_cha(nin)=HR_char        !*dum1
!
          dum2=Q_MF_vola(nin)
          LH_CHO(nin)=Latent_H_CHO*dum2
          latent(II_CHO,nin)=Latent_H_CHO
!
          dum3=MF_H2O_p(nin)
          LH_H2O(nin)=Latent_H_H2O*dum3
          latent(II_H2O,nin)=Latent_H_H2O
!
          DH_CHO(nin)=COALCV          !J/KG
     &               -HR_cha(nin)*dum1
     &               +LH_CHO(nin)!*dum2
     &               +LH_H2O(nin)!*dum3
!BABA
          DH_CHO(nin)=DH_CHO(nin)/dum2
!
!
          M_CHO(nin)=(12.d0+b_H(nin)/a_C(nin)+16.d0*c_O(nin)/a_C(nin))
     &              /1000.d0
          DH_CHO(nin)=DH_CHO(nin)*M_CHO(nin)  !J/mol
          wm(II_CHO)=M_CHO(nin)
          r_wm(II_CHO)=1.d0/M_CHO(nin)
          Ri(II_CHO)=gascns*r_wm(II_CHO)
!
! --- 2C+O2=>2CO+Q1
! --- 2CO+O2=>2CO2+Q2
! --- 2C+2O2=>2CO2+2Q3  => Q3=HR_cha(nin)
!
          dum1=2.d0*HR_cha(nin)*wm(II_CHO)-
     &         (2.d0*href(II_CO2)-2.d0*href(II_CO)-href(II_O2))
          HR_cha(nin)=dum1/WM(II_C)       !J/kg
          dum1=-(href(II_CO)-0.5d0*href(II_O2)-href(II_C))/WM(II_C)
          HR_cha(nin)=dum1
!     
          do q=1,nneq
          if(coal_CHO(q)==1) then
            vreq(II_CHO,1,q)=1.d0 !a_C(nin)
            vreq(II_CHO,2,q)=0.d0
            vreq(II_O2,1,q)=0.5d0*(a_C(nin)+0.5d0*b_H(nin)-c_O(nin))
            vreq(II_O2,2,q)=0.d0
            vreq(II_CO,1,q)=0.d0
            vreq(II_CO,2,q)=a_C(nin)
            vreq(II_H2O,1,q)=0.d0
            vreq(II_H2O,2,q)=0.5d0*b_H(nin)
          endif
          enddo
          
          if(my_rank==0) then
            write(ifll,'(2X,a,f10.3)') 'MSG: M_CaHbOc= ',M_CHO(nin)
          endif
!
          H_CO2a(nin)=H_CO2    !/M_CO2  !(J/mol) 
          H_H2Oa(nin)=H_H2O    !/M_H2O
          H_O2a(nin)=H_O2      !/M_O2
          xx_O2(nin)=(2.d0*a_C(nin)+b_H(nin)/2.-c_O(nin))/2.d0
          H_CHO(nin)=a_C(nin)*H_CO2a(nin)
     &              +b_H(nin)/2.*H_H2Oa(nin)
     &              -xx_O2(nin)*H_O2a(nin)
     &              +DH_CHO(nin)           ![J/mol] 
!--------------------------
! --- thermal properity
!--------------------------

          if(nin==1.and.II_CHO/=0) then
            t=tref_comp(II_CHO)
            a7(6,1,II_CHO)=H_CHO(nin)/8.314d0
     &                   -T*(a7(1,1,II_CHO)
     &                   +T*(a7(2,1,II_CHO)*0.5d0
     &                   +T*(a7(3,1,II_CHO)*1.d0/3.d0
     &                   +T*(a7(4,1,II_CHO)*0.25d0
     &                   +T*(a7(5,1,II_CHO)*0.2d0)))))
            a7(6,2,II_CHO)=H_CHO(nin)/8.314d0
     &                   -T*(a7(1,2,II_CHO)
     &                   +T*(a7(2,2,II_CHO)*0.5d0
     &                   +T*(a7(3,2,II_CHO)*1.d0/3.d0
     &                   +T*(a7(4,2,II_CHO)*0.25d0
     &                   +T*(a7(5,2,II_CHO)*0.2d0)))))
            rj = judg_temp(II_CHO,tref_comp(II_CHO))
            dum1=a7(6,rj,II_CHO)
     &      +t*(a7(1,rj,II_CHO)
     &      +t*(a7(2,rj,II_CHO)*0.5d0
     &      +t*(a7(3,rj,II_CHO)*1.d0/3.d0
     &      +t*(a7(4,rj,II_CHO)*0.25d0
     &      +t*(a7(5,rj,II_CHO)*0.2d0)))))
            href(II_CHO)=dum1
!dum1*dble(act(II_CHO))    ![---]
            hform(II_CHO)=href(II_CHO)*Ri(II_CHO)  ![J/kg]
          endif
!
          dum3=0.d0
          do icom=1,ncomp
          if(ys_p(icom)==undef) then
            write(ifle,*) '### ERR-11 : data error'
          write(ifle,*) 'lack of data : ys_p(ICOM) in &particle_model'
            write(ifle,*) 'ICOM =',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
            if(icom==II_CHO) then
              particle_comp(nin,icom)=Q_MF_vola(nin)-net_MF_HCN(nin)
            elseif(icom==II_H2O) then
               particle_comp(nin,icom)=Q_MF_H2O(nin)
            elseif(icom==II_C) then
              particle_comp(nin,icom)=Q_MF_fixC(nin)-net_MF_Nchar(nin)
            elseif(icom==II_ASH) then
              particle_comp(nin,icom)=Q_MF_ASH(nin)
            elseif(icom==II_HCN) then
              particle_comp(nin,icom)=net_MF_HCN(nin)
            elseif(icom==II_N) then
              particle_comp(nin,icom)=net_MF_Nchar(nin)
            else
              particle_comp(nin,icom)=ys_p(icom)
            endif
            dum3=dum3+particle_comp(nin,icom)
          endif
          enddo
          if(abs(dum3-1.d0)<1.d-10) then
            do icom=1,ncomp
            particle_comp(nin,icom)=particle_comp(nin,icom)/dum3
            enddo
          else
! --- 
            write(ifle,'(1X,a,E15.7,a)') 
     &    'MSG: Particle YS=',dum3, '/=1.d0'
            call FFRABORT(1,'ERR: SUM(YS) of net mass-fraction /= 1')
          endif
!
!          if(cal_CHAR_COMB<0) then
!            call FFRABORT(1,'ERR:cal_CHAR_COMB or 1')
!          else
!            ical_CHAR_COMB(nin)=cal_CHAR_COMB
!          endif
!
          if(cal_evap_coal<0) then
            call FFRABORT(1,'ERR:cal_evap_coal=0 or 1')
          else
            ical_la_coal(nin)=cal_evap_coal
          endif
!
!
        elseif(ifuel(nin)==iliquid) then
!----------------------
! --- liquid fuel 
!----------------------
          do icom=1,ncomp
          if(ys_p(icom)==undef) then 
            write(ifle,*) '### ERR-11 : data error'
            write(ifle,*) 'lack of data : ys(ICOM) in &particle_model'
            write(ifle,*) 'ICOM =',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
            particle_comp(nin,icom)=ys_p(icom)
          endif
          enddo
          if(liquidmu==undef) then
            write(ifle,'(1X,a)') 
     &       'lack of data : [liquidmu] in &particle_model'
            write(ifle,*) 'MSG: Mu for liquid particle'
            call FFRABORT(1,'namelist/particle_model')
          endif
          if(stens==undef) then
            write(ifle,*) 'lack of data : [stens] in &particle_model'
            write(ifle,*) 
     &        'MSG: stens for liquid particle surface tension'
            call FFRABORT(1,'namelist/particle_model')
          endif
!----------------------
! --- 
!----------------------
        elseif(ifuel(nin)==iglass.or.ifuel(nin)==isolid) then
          do icom=1,ncomp
          if(ys_p(icom)==undef) then 
            write(ifle,*) '### ERR-11 : data error'
            write(ifle,*) 'lack of data : ys(ICOM) in &particle_model'
            write(ifle,*) 'ICOM =',ICOM
            call FFRABORT(1,'namelist/particle_model')
          else
            particle_comp(nin,icom)=ys_p(icom)
          endif
          enddo
        endif
!-------------------------------------
! --- 
!-------------------------------------
        if(l_temp==1) then
          if(cp_p==undef) then
            write(ifle,*) 'lack of data : [cp_p] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(T_ini_p==undef) then
            write(ifle,*) 'lack of data : [T_ini_p] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(rad_emissivity==undef) then
            write(ifle,*) 
     &      'lack of data : [rad_emissivity] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!

          if(heat_fraction==undef) then
            write(ifle,*) 
     &      'lack of data : [heat_fraction] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(Av_CHO==undef.and.ifuel(nin)==icoal) then
            write(ifle,*)
     &      'lack of data : [Av_CHO] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(Ev_CHO==undef.and.ifuel(nin)==icoal) then
            write(ifle,*) 
     &      'lack of data : [Ev_CHO] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(Av_HCN==undef.and.ifuel(nin)==icoal) then
            write(ifle,*) 
     &      'lack of data : [Av_HCN] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(Ev_HCN==undef.and.ifuel(nin)==icoal) then
            write(ifle,*) 
     &      'lack of data : [Ev_HCN] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(Av_H2O==undef.and.ifuel(nin)==icoal) then
            write(ifle,*) 
     &      'lack of data : [Av_H2O] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          if(Ev_H2O==undef.and.ifuel(nin)==icoal) then
            write(ifle,*) 
     &      'lack of data : [Ev_H2O] in &particle_model'
            call FFRABORT(1,'namelist/particle_model')
          endif
!
          Tini_p(nin)=T_ini_p
          Cp_pa(nin)=cp_p
          eps_rad(nin)=rad_emissivity
!          cond_p(nin)=conduct_p
          heat_f(nin)=min(1.d0,max(heat_fraction,0.d0))
          Tref_p(nin)=T_ref_p
!
          Avap_CHO(nin)=Av_CHO
          Evap_CHO(nin)=Ev_CHO
          Avap_H2O(nin)=Av_H2O
          Evap_H2O(nin)=Ev_H2O
          Avap_HCN(nin)=Av_HCN
          Evap_HCN(nin)=Ev_HCN
        endif
!
        if(bc_no==0) then
        else
          do nb=1,nbcnd
          NO=nobcnd(nb)
          if(bc_no==NO) then
            bc_nb(nin)=nb
            NO=0
            exit
          endif
          enddo
          if(NO/=0) then
            call FFRABORT(1,'ERR: bc_no NOT existed')
          else
            NO_BC(nin)=bc_no
          endif
        endif
!
        P_MF(nin)=output_MF 
!
      enddo
!
      
      if(N_injctr==0) then
        call FFRABORT
     & (1,'ERR:&particle_model namelist NOT denfined in fflow.ctl')
      endif
!
! --- coal
!
      HPC_dead=-MXCVFAC-100
      P_Reced=-MXCVFAC-1
      A_dead=-MXCVFAC-2
      BC_dead=-MXCVFAC-3
      MASS_dead=-MXCVFAC-4
      return
!
 100  call FFRABORT(1,'Err100: Reading namelist [particle_model]')
 200  call FFRABORT(1,'Err200: Reading namelist [injector_number]')
!
      end subroutine inputdata

!
!=======================================================
      subroutine particle_chem
     &          (imode,NIN,N_injctr,ifli,ifll,ifle,
     &           iprtcle,my_rank,ncomp,
     &           cntlnam,lcomp,igT_iter,ireac,ierror)
!=======================================================
      use module_species,only : lenspc,spcnam,
     &                          chk_spec=>chkncomp,
     &                          wm
!
      implicit none
!
!
! --- [dummy arguments]
!
      integer,intent(in)     :: ifli,ifll,ifle,my_rank,igT_iter,
     &                          ireac,N_injctr
      integer,intent(in)     :: imode,NIN,iprtcle,ncomp
      character(*),intent(in):: cntlnam
      logical,intent(in)     :: lcomp
      integer,intent(out)    :: ierror
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
      real(8)  :: af,alphaf,ef                   ! elementary reaction : 
                                                 ! k = (af)*T^(alpahf)*EXP(-(ef)/(R*T))
      integer :: ignition_iter=0,stop_iter=10000000 
      integer :: nol,nor   
      real(8),parameter :: undef=-huge(1.)
      integer,parameter :: iundef=-huge(1)
!----------------------------------------------
      real*8 :: Acomb_char,Ecomb_char,p0
      real*8 :: F_char,eta_N_NO
      real*8 :: A_no1,f_index,E_no1,A_no2,E_no2
      real*8 :: A_re_NO,E_re_NO,A_E
      real*8 :: A_promp,E_promp
!----------------------------------------------
      real(8)  :: nu_fuel
!
      namelist /particle_chemreac/ 
     &                    model,
     &                    nleft,nright,
     &                    cleft,cright,
     &                    pleft,pright,
     &                    a,alpha,beta,e,
     &                    creac,nu_fuel,
     &                    af,alphaf,ef,
     &                    ignition_iter,stop_iter,
!
     &                    Acomb_char,Ecomb_char,P0,
!
     &                    F_char,eta_N_NO,
!
     &                    A_no1,f_index,E_no1,A_no2,E_no2,
!
     &                    A_re_NO,E_re_NO,A_E,
!
     &                    A_promp,E_promp
!
! --- [local entities]
!
      character(LEN=11),parameter :: subnam='inputdata'
!
      character(LEN=12),parameter :: mdllst(6)=(/
     &             'P_CHAR_COMB ','P_CHAR_NO   ',
     &             'P_NO_N2O2   ',
     &             'P_elementary','P_overall   ',
     &             'P_user      '/)
!
!----------------------------------------------------------------
!
      integer :: i,j,k,ios,imodl,s,kd,is,idum,idum1
!      
      integer :: q,ierr1,ns
      real(8) :: dd2,dum,dum1,dum2,dum3
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
!----------------------------------
! --- ireac 
!----------------------------------
!      if(ireac==0) return
!
      ierror=0
!
!---------------------
!-< 1. Initial set >-
!---------------------
!      if(imode==1) then
        nqP=0  
        nneqP=0
        rewind ifli
!
        do
        ierr1=0
        read(ifli,particle_chemreac,iostat=ios)
        if(ios<0) exit
        nqP= nqP+1
        call nml_errmsg0(ifle,ios,'particle_chemreac',ierr1)
        if( ierr1/=0 ) then
          write(ifle,'(1X,a,I8)') 
     &     'ERR: sequence no. of the namelist =',nqP
          call FFRABORT(1,'&particle_chemreac')
        endif
        enddo
!
        nneqP=nqP
        if(nneqP>0) then
          if(.not.lcomp) then
            write(ifle,*) ' ### Warning: Chemical reaction must use ',
     &                     'compressible or low-Mach approximation'
            call FFRABORT(1,' Change [flow] in [&model]')
          endif
        endif
!

        if(nqP>0) then
          ical_P_R=1
        else
          ical_P_R=0
        endif
!
        if(nqP==0) return  ! no reaction flow
        if(ierror.ne.0) goto 9999
!
!---------------------------
!--< 1.2 check interface >--
!---------------------------
!
        ns=ncomp
        call nml_chksiz(ifle,'mcomp',mcomp,ns,    ! if ns>mcomp, ierr1=1
     &                    modnam//subnam,ierr1)
        if(ierr1/=0) then
          ierror=1
          call FFRABORT(1,"ERR: ns>mcomp in ")
        endif
!---------------------------
!--< 1.3 allocate arrays >--
!---------------------------
        allocate( ireqP(     nqP),  ! flag number of reaction model
     &            vreqP(ns,2,nqP),  ! stoichiometric coefficients
     &            slidP(ns,2,nqP), 

     &            preqP( 6,2,nqP),  ! parameters for reaction rate
     &            creqP(ns,  nqP),  ! index number of mol
                                    ! concentration(overall)
     &        sub_vreqP(ns,  nqP),  ! subtraction of "vreq" of
     &        sub_slidP(ns,  nqP),  ! 
                                    ! each species in each equation
     &      sigma_vreqP(     nqP),  ! change of "sub_vreq" in each equation
     &         ig_iterP(     nqP),
     &        stp_iterP(     nqP),
     &               stat=ierr1 )
        if(ierr1/=0) then
          write(ifle,*) '### error : allocation failed'
          ierror = 1
        endif
!
      allocate(
     &        Acomb_C(nneqP),
     &        Ecomb_C(nneqP),
     &        p00_CHAR(nneqP),

     &        FF_char(nneqP),
     &        eeta_N_NO(nneqP),

     &        AA_no1(nneqP),
     &        EE_no1(nneqP),
     &        AA_no2(nneqP),
     &        EE_no2(nneqP),
     &        ff_index(nneqP),

     &        AA_re_NO(nneqP),
     &        EE_re_NO(nneqP),
     &        AA_E(nneqP),

     &        AA_promp(nneqP),
     &        EE_promp(nneqP),

     &             stat=ierr1 )      
!
      Acomb_C(:)=0.d0
      Ecomb_C(:)=0.d0
      p00_CHAR(:)=0.d0
      FF_char(:)=0.d0
      eeta_N_NO(:)=0.d0
      AA_no1(:)=0.d0
      EE_no1(:)=0.d0
      AA_no2(:)=0.d0
      EE_no2(:)=0.d0
      ff_index(:)=0.d0
      AA_re_NO(:)=0.d0
      EE_re_NO(:)=0.d0
      AA_E(:)=0.d0
      AA_promp(:)=0.d0
      EE_promp(:)=0.d0      
!
!      endif 
!---------------------------
! --- initialization 
!---------------------------
      ireqP=0
      vreqP=0.d0
      slidP=0
      preqP=0.d0
      creqP=0.d0
      sub_vreqP=0.d0
      sigma_vreqP=0.d0
      ig_iterP=0
      stp_iterP=10000000
!--------------------------
! --- P_CHAR_COM
!--------------------------
      Acomb_char=undef
      Ecomb_char=undef
      p0=undef
!-------------------------
! --- P_CHAR_NO
!-------------------------
      F_char=undef
      eta_N_NO=undef
!---------------------
! --- FUEL NO
!---------------------
      A_no1=undef
      f_index=undef
      E_no1=undef
      A_no2=undef
      E_no2=undef
!---------------------
! --- Fuel NO reduction 
!---------------------
      A_re_NO=undef
      E_re_NO=undef
      A_E=undef
!---------------------
! --- 
!---------------------
      A_promp=undef
      E_promp=undef
!---------------------------
!-< 2. Input namelist >-
!---------------------------
!--< 2.1 read namelist >--
!---------------------------
      rewind ifli

      do 200 q=1,nqP

        ignition_iter=0
        stop_iter=10000000

        model= ' '

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
!
        nu_fuel=undef    ! stoichiometric coefficient of fuel
!/ elementary /
        af    = undef    ! k = (af)*T^(alpahf)*EXP(-(ef)/(R*T))
        alphaf= undef
        ef    = undef

        pleft=undef
        pright=undef
!----------------------------------------------------
        read(ifli,particle_chemreac,iostat=ios)    ! read of namelist
        call nml_errmsg0(ifle,ios,'particle_chemreac',ierr1)     ! if ios>0,ierr1=1
        if( ierr1/=0 ) then
          write(ifle,*) "ERR: reading namelist 'particle_chemreac'"
          ierror = 1
          goto 9999
        endif
!-----------------------------
! --- identify reaction kind
!-----------------------------
!-------------------------
! --- identify model no.
!-------------------------
        call nml_listno(6,mdllst,model,imodl)
        call nml_chkchr0(ifle,'model',model,imodl,ierr1)
        if( ierr1.ne.0 ) then
          write(ifle,'(1X,a)') 
     &    'ERR: &particle_chemreac/model is unknown param.'
          write(ifle,*) 
     &    " ### MSG: model='overall',
     &    'elementary','eddybreakup','fire',",
     &    " or 'user'"
          call FFRABORT(1,'ERR: &particle_chemreac/model')
        endif
        ireqP(q)=imodl
!----------------------------------
! --- check gas/surface model
!----------------------------------
!-------------------------
! --- ignition iter
!-------------------------
      if(ignition_iter.le.0) then
        ig_iterP(q)=0
      else
        ig_iterP(q)=ignition_iter
      endif
      if(stop_iter.le.0) then
        stp_iterP(q)=0
      else
        stp_iterP(q)=stop_iter
      endif
!
      if(ig_iterP(q)<igT_iter) then
        write(ifle,*) ' ### WRN: Reaction no ',q
        write(ifle,*)
     & ' ### WRN: &particle_chemreac/ignition_iter < &chemcntl/igT_iter'
        write(ifle,*)
     & ' ### WRN: Reaction should start on iter= ',igT_iter
      endif
!--------------------------------------------------
! --- check duplication between left and right side
!--------------------------------------------------
      call chk_lr(nleft,nright)
! --- check duplication of left side
      call chk_dup('nleft',nleft)
! --- check duplication of right side
      call chk_dup('nright',nright)
!--------------------------------------------
! --- pick up number of species in each side
!--------------------------------------------
      call pick_up_species
!
      select case(imodl)
        case(P_CHAR_COM) 
          call read_P_CHAR_COM
        case(P_CHAR_NO)
          call read_P_CHAR_NO
        case(P_HCN_NO_N2O2)    
          call read_P_HCN_NO_N2O2
        case(P_ovrall)   
          call read_overall
        case(P_elemnt)
          call read_elementary
        case(P_user)      
          call read_user
      end select
 200  continue
      if(ierror/=0) 
     & call FFRABORT(1,'ERR: reading error in &particle_chemreac')
!------------------------------------------------
! --- 
!------------------------------------------------
!--------------------------------
! --- display of reading results 
!--------------------------------
!
      do q=1,nqP
      imodl=ireqP(q)
      if(imodl==P_ovrall) then
        call disp_overall(q)
      elseif(imodl==P_elemnt) then
        call disp_elementary(q)
      else
        call disp_elementary(q)
      endif
      enddo
!
! 300  call FFRABORT(1,'Err300: Reading namelist [injector_number]')
!----------------------
! --- final error check
!----------------------
!
      return
 9999 continue
      ierror=1
!
      return
!///////////////////////////////////////////////////////////////////////
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
!
        if(value==undef ) then
          write(ifle,*) '### error : data error'
          write(ifle,*) 
     &    'ERR: lack of data : ',str2,'(i)',' in reaction eq.=',q
          write(ifle,*) 'i =',c_no
          if(c_no==ns+1) write(ifle,*) 'ERR: Lack of 3rd M coeffient'
          call FFRABORT(1,'ERR:')
        endif
      end function c_no
!
!===================================================
       integer function c_nop(name,str1,value,str2)
!===================================================
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
      integer      :: s,no_elec,no_plsm,ii,i
      character(1) :: plasma,sharp,percnt
      real*8 :: wsum
      logical :: nml_comp_eq
      integer,save :: coal=0
      character(lenspc) :: spienm
!
      nol = 0            ! number of species in left side
      nor = 0            ! number of species in right side
!------------------------------------------------
! --- check stick,Bohm,Y_Bohm,E_BOHM model 
! --- search species number
!------------------------------------------------
      do s=1,ns
        if(nleft(s)/=' ') then
          do i=1,lenspc
          percnt=nleft(s)(i:i)
          ii=0
          if(percnt=='@') then
            nleft(s)=trim(nleft(s)(1:i-1))
            ii=1
            exit
          endif
          enddo
          j=c_no(nleft(s),'nleft',cleft(s),'cleft')
          if(j<=ns) then
            vreqP(j,1,q)=(cleft(s))
            if(ii==1) slidP(j,1,q)=1 
            nol=nol+1
          endif
        endif
! --- 
        if(nright(s)/=' ') then 
          do i=1,lenspc
          percnt=nright(s)(i:i) 
          ii=0
          if(percnt=='@') then
            nright(s)=trim(nright(s)(1:i-1)) 
            ii=1
            exit
          endif
          enddo
          j=c_no(nright(s),'nright',cright(s),'cright')
          if(j<=ns) then 
            vreqP(j,2,q)=(cright(s))
            if(ii==1) slidP(j,2,q)=1 
            nor=nor+1
          endif
        endif
      enddo
!----------------------------------------------
! --- 
!----------------------------------------------
      do s=1,ns
        sub_vreqP(s,q)=vreqP(s,2,q)-vreqP(s,1,q)
      enddo
!
      do s=1,ns
        sigma_vreqP(q)=sigma_vreqP(q)+sub_vreqP(s,q)
      enddo
!
!
!      if(nol<1 .or. nor<1) then               ! error check
!        write(ifle,*) "### error-1 : data error"
!        write(ifle,*) "lack of data in equation ",q
!        ierror = 1
!        write(ifle,*) "### error : nol & nor= ",nol,nor
!        stop 'nol<1 or nor<1'
!      elseif(imodl==2 .and. (nol>3 .or. nor>3)) then
!        write(ifle,*) "### error-2 : data error"
!        write(ifle,*) "too much data in equation ",q
!        write(ifle,*) "elementary model affords 3 species in side."
!        ierror = 1
!      elseif(imodl==3 .and. nol/=2) then
!        write(ifle,*) "### error-3 : data error"
!        write(ifle,*) "no. of elements of nleft =",nol
!         write(ifle,*) "it must be 2"
!        write(ifle,*) "nleft(1) = name of fuel"
!        write(ifle,*) "nleft(2) = name of oxidizing agent"
!        ierror = 1
!      endif
      
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
!          creqP(j,q) = creac(k)
        endif
        if(k<nol) then
          k=k+1
        else
          exit
        endif
      enddo
      if(.not.nml_comp_eq("a",a,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("a",a,0.0d0,1,ifle,ierror) ) then
          preqP(1,1,q)=a
        endif
      endif
  
      if(.not.nml_comp_eq("alpha",alpha,undef,1,ifle,ierror) )
     &  preqP(2,1,q)=alpha
      if(.not.nml_comp_eq("beta",beta,undef,1,ifle,ierror) )
     &  preqP(3,1,q)=beta
      if(.not.nml_comp_eq("e",e,undef,1,ifle,ierror) ) then
        if(nml_comp_gt("e",e,0.0d0,1,ifle,ierror) ) then
          preqP(4,1,q)=e
        endif
      endif
      end subroutine read_overall
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
          preqP(1,1,q)=af
        endif
      endif
      if(.not.nml_comp_eq("alphaf",alphaf,undef,q,ifle,ierror)) then
        preqP(2,1,q)=alphaf
      endif
      if(.not.nml_comp_eq("ef",ef,undef,q,ifle,ierror)) then
        preqP(3,1,q)=ef
      endif
! --------------------------------------------------------------------
! --------------------------------------------------------------------
      
!----------------------------
! --- arbitrary power index
!----------------------------
      if(imodl==P_elemnt) then
        do s=1,ns
        if(nleft(s)/=' ') then         ! left side
          j=c_nop(nleft(s),'nleft',pleft(s),'pleft')
          if(j<=ns) then 
            if(pleft(s)/=undef) vreqP(j,1,q)=pleft(s)
          endif
        endif
!
        if(nright(s)/=' ') then        ! right side
          j=c_nop(nright(s),'nright',pright(s),'pright')
          if(j<=ns) then
            if(pright(s)/=undef) vreqP(j,2,q)=pright(s)
          endif
        endif
        enddo
      endif
      end subroutine read_elementary
!
!=================================================
      subroutine read_user
!=================================================
      logical :: nml_comp_eq,nml_comp_gt,nml_comp_ge
!
      if(.not.nml_comp_eq("a",a,undef,q,ifle,ierror) ) then
        preqP(1,1,q)=a
      else
        write(ifle,'(2a,I4)') 
     &  'ERR: Prefactor [a] of Arrhenius NOT defined in &particle_chemreac',
     &  'at eq.=',q
        call FFRABORT(1,'ERR: Maybe [a] be defined, ')
      endif

      if(.not.nml_comp_eq("alpha",alpha,undef,q,ifle,ierror)) then
         preqP(2,1,q)=alpha
      else
         write(ifle,'(2a,I4)') 
     &  'ERR: Exponent [alpha: T^alpha] of Arrhenius NOT defined ',
     &  'in &chemcntl at eq.=',q
        call FFRABORT(1,'ERR: Maybe [alpha] be defined, ')
      endif

      if(.not.nml_comp_eq("e",e,undef,q,ifle,ierror)) then
         preqP(3,1,q)=e
      else
         write(ifle,'(2a,I4)') 
     &   'ERR: Activation [E]  of Arrhenius NOT defined ',
     &   'in &chemcntl at eq.=',q
         call FFRABORT(1,'ERR: reading error in user')
      endif
      end subroutine read_user
!
!=================================================
      subroutine read_P_CHAR_COM
!=================================================
      if(Acomb_char==undef) then
        write(ifle,'(1x,a,I4)') 
     &    'ERR: [Acomb_char] is NOT defined, reaction no= ',q
        write(ifle,'(1x,a)') 
     &  'MSG: [P_CHAR_COMB] model is defined '
        call FFRABORT(1,'STOP : &particle_chemreac')
      endif

      if(Ecomb_char==undef) then
        write(ifle,'(1x,a,I4)') 
     &    'ERR: [Ecomb_char] is NOT defined, reaction no= ',q
        write(ifle,'(1x,a)') 
     &  'MSG: [P_CHAR_COMB] model is defined '
        call FFRABORT(1,'STOP : &particle_chemreac')
      endif

      if(p0==undef) then
        write(ifle,'(1x,a,I4)') 
     &    'ERR: [p0] is NOT defined, reaction no= ',q
        write(ifle,'(1x,a)') 
     &  'MSG: [P_CHAR_COMB] model is defined '
        call FFRABORT(1,'STOP : &particle_chemreac')
      endif
      Acomb_C(q)=Acomb_char
      Ecomb_C(q)=Ecomb_char
      p00_CHAR(q)=P0
      return
      end subroutine read_P_CHAR_COM
!######################################################
!======================================================
      subroutine read_P_CHAR_NO
!======================================================
      if(F_char==undef) then
        write(ifle,'(1x,a,I4)') 
     &    'ERR: [F_char] is NOT defined, reaction no= ',q
        write(ifle,'(1x,a)') 
     &  'MSG: [P_CHAR_NO] model is defined '
        call FFRABORT(1,'STOP : &particle_chemreac')
      else
        FF_char(q)=F_char
      endif
!
      if(eta_N_NO==undef) then
        write(ifle,'(1x,a,I4)') 
     &    'ERR: [eta_N_NO] is NOT defined, reaction no= ',q
        write(ifle,'(1x,a)') 
     &  'MSG: [P_CHAR_NO] model is defined '
        call FFRABORT(1,'STOP : &particle_chemreac')
      else
        eeta_N_NO(q)=eta_N_NO
      endif
      return
      end subroutine read_P_CHAR_NO
!######################################################
!======================================================
      subroutine read_P_HCN_O2_NO
!======================================================
!
            if(A_no1==undef) then
              write(ifle,'(1x,a,I4)') 
     &        'ERR: [A_no1] is NOT defined, reaction no= ',q
              write(ifle,'(1x,a)') 
     &        'MSG: [P_HCN_O2_NO] model is defined '
               call FFRABORT(1,'STOP : &particle_chemreac')
            else
              AA_no1(q)=A_no1
            endif
            if(E_no1==undef) then
              write(ifle,'(1x,a,I4)') 
     &        'ERR: [A_no1] is NOT defined, reaction no= ',q
              call FFRABORT(1,'STOP : &particle_chemreac')
            else
              EE_no1(q)=E_no1
            endif
!
            if(A_no2==undef) then
              write(ifle,'(1x,a,I4)') 
     &        'ERR: [A_no1] is NOT defined, reaction no= ',q
              call FFRABORT(1,'STOP : &particle_chemreac')
            else
              AA_no2(q)=A_no2
            endif
            if(E_no2==undef) then
              write(ifle,'(1x,a,I4)') 
     &        'ERR: [A_no2] is NOT defined, reaction no= ',q
              call FFRABORT(1,'STOP : &particle_chemreac')
            else
              EE_no2(q)=E_no2
            endif
            if(f_index==undef) then
              write(ifle,'(1x,a,I4)') 
     &        'ERR: [f_index] is NOT defined, reaction no= ',q
              call FFRABORT(1,'STOP : &particle_chemreac')
            else
              ff_index(q)=f_index
            endif
!      
      end subroutine read_P_HCN_O2_NO
!######################################################

!======================================================
      subroutine read_P_HCN_NO_N2O2
!======================================================
      if(A_re_NO==undef) then
        write(ifle,'(1x,a,I4)') 
     &        'ERR: [A_re_NO] is NOT defined, reaction no= ',q
        call FFRABORT(1,'STOP : &particle_chemreac')
      else
        AA_re_NO(q)=A_re_NO
      endif
!
          if(E_re_NO==undef) then
            write(ifle,'(1x,a,I4)') 
     &        'ERR: [E_re_NO] is NOT defined, reaction no= ',q
            call FFRABORT(1,'STOP : &particle_chemreac')
          else
            EE_re_NO(q)=E_re_NO
          endif

          if(A_E==undef) then
            write(ifle,'(1x,a,I4)') 
     &        'ERR: [A_E] is NOT defined, reaction no= ',q
!            call FFRABORT(1,'ERR:A_E not defined')
          else
            AA_E(q)=A_E
          endif
      end subroutine read_P_HCN_NO_N2O2
!
!======================================================
      subroutine read_P_prompt_NO
!======================================================
      return
      if(A_promp==undef) then
        write(ifle,'(1x,a,I4)') 
     &        'ERR: [A_promp] is NOT defined, reaction no= ',q
        call FFRABORT(1,'STOP : &particle_chemreac')
      else
        AA_promp(q)=A_promp
      endif
      if(E_promp==undef) then
        write(ifle,'(1x,a,I4)') 
     &        'ERR: [E_promp] is NOT defined, reaction no= ',q
        call FFRABORT(1,'STOP : &particle_chemreac')
      else
        EE_promp(q)=E_promp
      endif
      return
      end subroutine read_P_prompt_NO
!
!======================================================
      subroutine read_P_Zeld_NO
!======================================================
      end subroutine read_P_Zeld_NO
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
!===================================================
      character(lenspc),intent(in) :: left(mcomp+1)
                                   ! name of species on left side
      character(lenspc),intent(in) :: right(mcomp+1)
                                   ! name of species on right side
      integer :: s
!
      do s=1,ns
        if(left(s)/=' '.and. left(s)/='M') then
          do j=1,ns
          if(right(j)/=' '.and.right(j)/='M') then
            if( left(s)==right(j) ) then
            write(ifle,'(a)') '### error : data error'
            write(ifle,'(3a,I4)') 
     &      'ERR: duplicated data : nleft(i)=nright(j)='
     &                    ,trim(left(s)),' in eq.=',q
            write(ifle,*) 'i,j =',s,j
            ierror=1
            call FFRABORT(1,'STOP AT: chk_lr')
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
      if(my_rank==0) write(ifll,'(a5,$)') arrow(preqP(1,2,1))
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
      end subroutine disp_elementary
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
     &                        trim(mdllst(ireqP(q)))
      if(my_rank==0) 
     & write(ifll,"(2x,a,i3)")
     &  "Total number of chemical equation : ",nqP
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
        if( creqP(s,q)/=0 ) then
          if(my_rank==0) write(ifll,"(2x,a6,$)") trim(spcnam(s))
          if(my_rank==0) write(ifll,"(a,f5.1,$)") " ->",creqP(s,q)
        endif
      enddo
      end subroutine disp_creq
!
!-------------------------------------------------
! --- display 1 equation for overall reaction 
!=================================================
      subroutine disp_eq1(isw)
!=================================================
      integer,intent(in) :: isw  ! switch for side of equation
      integer :: s
      k = 0
      do s=1,ns
        if(vreqP(s,isw,1)/=0) then
          if( k>0 .and.my_rank==0) write(ifll,"(a3,$)") " + "
          if(my_rank==0) write(ifll,"(f5.1,$)") vreqP(s,isw,1)
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
      if(my_rank==0) then
        call str_eqs(1,q,sinp)
        write(ifll,"(1x,i3,2x,a,$)") q,trim(sinp) ! left side
      endif
      if(my_rank==0) then
        write(ifll,"(a5,$)") arrow(preqP(1,2,q))
      endif
      if(my_rank==0) then
        call str_eqs(2,q,sinp)
       write(ifll,"(a,$)") trim(sinp)          ! right side
      endif
      call disp_coefficients(q,3)
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
        write(ifll,"(e12.3,$)") (preqP(i,1,q),i=1,n)
        write(ifll,"(/,7x,a,$)") "backward reaction : "
        if(preqP(1,2,q)==undef) then
          write(ifll,"(a,$)") "<-- using equilibrium constants"
        elseif(preqP(1,2,q)==0) then
          write(ifll,"(a,$)") "<-- NOT consider"
        else
          write(ifll,"(e12.3,$)") (preqP(i,2,q),i=1,n)
        endif
      endif
      end subroutine disp_coefficients
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
          if(vreqP(s,isw,q)>0 ) then
            k=k+1
            if(k/=1) then
               idum=lensum+1
               lensum=lensum+3
               sinp(idum:lensum)=' + '
            endif
            write(strrr,'(F4.2)') vreqP(s,isw,q)
            idum=lensum+1
            lensum=lensum+4
            sinp(idum:lensum)=strrr
            lenc=len_trim(spcnam(s))
            idum=lensum+1
            lensum=lensum+lenc
            sinp(idum:lensum)=adjustl(spcnam(s))
          endif
        enddo
!
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

      end subroutine particle_chem
!
!      end subroutine inputdata
      end module module_particle
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FFR_TRACING(IT,ITS,DT0,ffexit,time,reaction,
     1           MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_CFIDX,MAT_CAL,
     2           LVEDGE,SFCENT,CVVOLM,WIFACE,CVCENT,SFAREA,
     3           GRDC,RHO,RMU,VEL,TMP,YYS,WDOT,
     4           DNORM,NEIGHBOR,
     5           INSIDE,JOUT,PMASS,DIA,XYZP,UVWP,PARCEL,PYS,TMPP,RHOP,
     &           DYDTP,DHDTP,HIS_P,
     6           HEIGHT,FU,FV,FW,IERROR,
     7           MOMCELL,EVAP_Y,EVAP_H,LBC_SSF,vctrp,
     &           P_rad,radinten,
     8           IPCELL,NEIBCELL,aks,MASSYS,CP,
     7           pp0,rds,rmd,hhs,cr,FIELD_U,WDRBND
     9           ) 
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments] 
!
      USE module_dimension 
      USE module_hpcutil 
      USE module_particle     
      USE module_constant,only : SML
      use module_io,only       : ifle,ifll,gdScale
      use module_trace ,only   : CC,HEIGHT0,HEIGHT1,VELN
      use module_time,    only : itere
      use module_particle,only : ifuel,icoal,iliquid,particle_comp,
     &                           iglass,isolid,total_start
      use module_particle,only : II_CHO,II_H2O,II_C,II_ash
      use module_particle,only : IFLAG_COAL,VS_CHO,V_CHO,V_H2O,VS_H2O,
     &                           VS_HCN,V_HCN
      use module_metrix,only   : ys,rcomp,eva_comp,sinj_comp
      use module_species,only  : gascns,wm,r_wm,Ri,act,sw
      use module_model,  only  : encomp,ical_t,ical_reac,ical_prt
      use module_model,  only  : ical_WDR
      use module_chemcntl,only : igT_iter,Zeldovich_P,Fuel_NO_P,
     &                           prompt_NO_P
!      use module_chemreac,only : ical_P_R
      use module_trace ,only   : WK1_P,PRT_VF
      use module_output,only   : outchk,fnmlst,none,yes,mlt_rslt
      use module_time,    only : iters
      USE module_usersub, ONLY : src_p_ini,usryes
      use module_metrix,only   : WDOTP   =>d2work2
      use module_metrix,only   : tempp   =>d2work3
!      use module_metrix,only   : WDOTP_P =>d2work1
      
      use module_rad   , only  : NDIV
      use module_rad   , only  : radflag
      use module_model,   only : monitor_stp
      use module_boundary,only : nbcnd
      use module_scalar,only   : calxi,idifflm,ORG,ical_flmlt
!
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)    :: IT,ITS
      logical,intent(in)    :: ffexit
      REAL*8,INTENT(INOUT)  :: DT0
      REAL*8,INTENT(IN)     :: time
      INTEGER,INTENT(IN)    :: MAT_NO   (      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_INDEX(      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CFIDX(      0:MXMAT)
      LOGICAL,INTENT(IN)    :: MAT_CAL  (      0:MXMAT)
      LOGICAL,INTENT(IN)    :: reaction
      INTEGER,INTENT(IN)    :: MAT_CV   (      MXALLCV)
!
      INTEGER,INTENT(IN)    :: LVEDGE   (    2,MXCVFAC)
      REAL*8 ,INTENT(IN)    :: SFCENT   (    3,MXCVFAC)
      REAL*8 ,INTENT(IN)    :: CVVOLM   (      MXALLCV)
      REAL*8 ,INTENT(IN)    :: WIFACE   (    2,MXCVFAC)
      REAL*8 ,INTENT(IN)    :: CVCENT   (    3,MXALLCV)
      REAL*8 ,INTENT(IN)    :: SFAREA   (    4,MXCVFAC)
      INTEGER,INTENT(IN)    :: LBC_SSF  (      MXSSFBC)
!
      REAL*8,INTENT(INOUT)  :: GRDC     (      MXALLCV,3,3)
      REAL*8,INTENT(IN)     :: RHO      (      MXALLCV,2)
      REAL*8,INTENT(IN)     :: RMU      (      MXALLCV)
      REAL*8,INTENT(IN)     :: VEL      (      MXALLCV,3,2)
      REAL*8,INTENT(INOUT)  :: WDOT     (      MXALLCV,MXCOMP)
      REAL*8,INTENT(IN)     :: YYS      (      MXALLCV,MXCOMP)
      REAL*8,INTENT(IN)     :: TMP      (      MXALLCV)
      real*8,intent(in)     :: aks      (      MXALLCVR,MXRANS)
!
      REAL*8,INTENT(INOUT)  :: RHOP     (      MXPRT)
      REAL*8,INTENT(INOUT)  :: TMPP     (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: DNORM    (      MXCVFAC_P)
      INTEGER,INTENT(IN)    :: NEIGHBOR (MEP  ,MXALLCV_P)
!
      INTEGER,INTENT(INOUT) :: INSIDE   (      MXPRT)
      INTEGER,INTENT(INOUT) :: JOUT     (      MXPRT)
      
      REAL*8,INTENT(INOUT)  :: PMASS    (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: PYS      (      MXPRT,MXCOMP,2)
      REAL*8,INTENT(INOUT)  :: PARCEL   (      MXPRT)
      REAL*8,INTENT(INOUT)  :: DIA      (      MXPRT)
      REAL*8,INTENT(INOUT)  :: XYZP     (      MXPRT,3)
      REAL*8,INTENT(INOUT)  :: UVWP     (      MXPRT,3)
!
      REAL*8,INTENT(INOUT)  :: HEIGHT   (  MEP,MXPRT)
      REAL*8,INTENT(INOUT)  :: FU       (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: FV       (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: FW       (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: DHDTP    (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: DYDTP    (      MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: HIS_P    (      MXPRT)
      INTEGER,INTENT(INOUT) :: IERROR

      REAL*8,INTENT(INOUT)  :: MOMCELL  (      MXALLCV_P,4)
      REAL*8,INTENT(INOUT)  :: EVAP_Y   (      MXALLCV_P,MXCOMP)
      REAL*8,INTENT(INOUT)  :: EVAP_H   (      MXALLCV_P,2)

      INTEGER,INTENT(INOUT) :: IPCELL   (      MXALLCV_P)
      INTEGER,INTENT(INOUT) :: NEIBCELL (      MXALLCV_P)
      integer,intent(inOUT) :: vctrP     (MXCV_VP,0:MXBND_VP)
      REAL*8,INTENT(INOUT)  :: P_rad    (MXALLCV_RAD,2)
      REAL*8,INTENT(IN)     :: radinten (NDNG,MXALLCV_RAD) 
      REAL*8,INTENT(INOUT)  :: MASSYS(   MXALLCV,MXcomp)
      REAL*8,INTENT(INOUT)  :: cp    (        MXALLCV)
      REAL*8 ,INTENT(IN)    :: pp0   (        MXALLCV)
      REAL*8 ,INTENT(IN)    :: RDS   (        MXALLCV,MXcomp)
      REAL*8,INTENT(IN)     :: RMD   (        MXALLCV)
      REAL*8 ,INTENT(INOUT) :: hhs   (        MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cr    (        MXALLCV)
      real*8 ,intent(inout) :: FIELD_U(       MXCV_D,NFLID)
      real*8 ,intent(inout) :: WDRBND(        MXCV_WDR,N_inj)
!
! --- [local entities]
!
!      REAL*8   :: ReynoldsPP(MXPRT),PMASS0(MXPRT)
      REAL*8   :: Nusselt,Prandtl,ReynoldsPP
      INTEGER      :: IERR,MCELL
      INTEGER      :: NP,MFACE,IPCV
      integer, save,allocatable :: NPSUM(:),NPINJ(:,:)

!
!------------------------------------------------
! --- VARIABLES FOR FORCE FEEDBACK CALCULATION 
!------------------------------------------------
!-------------------------------------------------------------
!     NOTES:
!     MEVENT: TOTAL NUMBER OF EVENTS
!     ONE EVENT: ONE PATICALE PASS THROUGH ONE CELL
!     JPASS(:) :: CELL NUMBER WHERE PARTICLE PASSING THROUGH
!-------------------------------------------------------------
! --- HPC - MPI 
!--------------------
!----------------------------------------
! for random walk calc. (k-e) 
!----------------------------------------
!
      INTEGER :: I,J,K,L,M,N,N0,N00,K1,K2,IR,IL,J2,IE00,INB
      INTEGER :: IFA,MOUT,KE,KE0,IRK,NTIME,NSTA1,NP01,
     &           NEIB,NCPU
!
!-------------------------------------------------------
!-< 1. SET UP SOME ARRAYS >-	PROGRAM PARTICLETRACING 
!-------------------------------------------------------
!
      REAL*8       :: PI,D_p,Kc,Kd
      INTEGER,save :: ITER0
      CHARACTER*19 :: OUTNAMES
!-----------------
! --- WORK AREA
!-----------------
!
!-------------------------------------------------------
! For Breakup Calculation
!-------------------------------------------------------
      REAL*8 :: dTbu,dTemp,dum1,dum2,dum3,Pg_O2,dumi
      real*8 :: dum4,zzz
      REAL*8 :: dDiamReduc,dDIANEW,dVolume
      REAL*8 :: DIM_U,RHO_U,T_U
      LOGICAL:: lBreakup,pat_reac
      integer:: iBreakupCount,out,IIII,NNN,idiv
!-------------------------------------------------------
! For MPI Output 
!-------------------------------------------------------
      integer :: ITER_IN_COUNT,istat,IP,ICFL,ICFLL,INI,IPS,IPE
      integer :: IPP,IPPP,niter=4
      real*8  :: DT_IN, dTSave,DT00
!-----------------------
! --- [local entities] 
!-----------------------
      REAL*8 :: UVWP0(3),DT(4),XYZP0(3),pdum
      REAL*8 :: VNORMAL,TIMEMIN,REYNOLDSP,SLIPV,FCD,US,VS,WS
      REAL*8 :: DIA0,RHO0,XMAX,YMAX,ZMAX
!
      REAL*8 :: XP0,YP0,ZP0,BETAP0,ALFA0,U0,Kv
!
      INTEGER           :: MAXBUF,IMODE
      INTEGER,PARAMETER :: NDIS  =20
      INTEGER,save :: iflag=0,i_initial=0,iflag1=0,reinj=0
      INTEGER           :: I_HPC=0,NRR
      INTEGER :: NUMINJCT, NPALLTOTAL,icom,NREINJ=0,NPMX,
     &           NpAllTotal2
      REAL*8  :: dTempRandomWalk
      REAL*8  :: AREP,diag,SORC,QRAD,H_CHAR_COM,THR,RADI,wegt
      REAL*8  :: comb_k,fuel_char_NO,X_HCN,X_NO,X_O2,dn1,dn2
      integer :: IIMAT,IMAT,ICVS,ICVE,ICVL
      integer :: ifl=0,ierr1=0,N_O2,NB
      integer :: nbx,icpu,idum,idum2
!
      REAL*8  ::  DP_M,TP_M,RHOP_M,MASSP_M,VELP_M(1:3)
      REAL*8  ::  G_T_M,G_VEL_M(1:3),rmu_G_M,G_RHO_M,G_RAMD_M,G_CP_M
      real*8  ::  Pr_bulk,Re_d,htc
!
!----------------------------------------------------------------------
!
      MAXBUF    =MXPRT
      CC(:)     =0.d0 
      HEIGHT0(:)=0.d0 
      HEIGHT1(:)=0.d0 
      VELN(:)   =0.d0 
!
!-----------------------------------
! --- first step ; restart step 
!-----------------------------------
!
      if(iflag==0) then
        IFLAG_COAL=0
        do INI=1,N_injctr
        if(ifuel(INI)==icoal.or.ifuel(INI)==iglass) then
           IFLAG_COAL=icoal
        endif
        enddo
!
        if(ical_P_R==1) then
!          ALLOCATE(WDOTP(MXCV,MXcomp),stat=ierr1) 
!          if(ierr1/=0) call 
!     &    FFRABORT(1,'allocating WDOTP error in FFR_TRACING')
!          WDOTP(:,:)=0.d0
        endif

        allocate(TRESIDUAL(MXPRT),stat=ierr)
        if(ierr/=0) call 
     &  FFRABORT(1,'ERR:allocate error-1 in FFR_TRACING') 
!
        allocate(PMASSC(MXPRT),stat=ierr)
        if(ierr/=0) call 
     &    FFRABORT(1,'ERR:allocate PMASSC in FFR_TRACING') 
        PMASSC(:)=0.d0
!
        if(NPE>1) then 
          allocate(IMPORT_indexIP(0:NEIBPETOT))
          allocate(EXPORT_indexIP(0:NEIBPETOT))
          allocate(IPSND(MXPRT)) 
          allocate(IPQUENE(MXPRT)) 
          allocate(NPEXCHAGOUT(NEIBPETOT)) 
          ALLOCATE(NPEXCHAGIN(NEIBPETOT)) 
!
          IMPORT_indexIP=0
          EXPORT_indexIP=0
          IPSND(:)=0
          IPQUENE(:)=0
          NPEXCHAGOUT(:)=0
          NPEXCHAGIN(:)=0
          NPOUTALL=0
          NPINALL=0
        endif
!
        ALLOCATE(NPSUM(NPE),NPINJ(NPE,N_injctr))
        NPSUM(:)=0
        NPINJ(:,:)=0
        ALLOCATE(IBC_P_NO(1:NPE,N_injctr,nbcnd,2))
        ALLOCATE(IBC_P_TOL(N_injctr,nbcnd,2))
        IBC_P_NO(:,:,:,:)=0
        IBC_P_TOL(:,:,:)=0

!
        NREINJ=0
        NPMX=0
        do INI=1,N_injctr
        NREINJ=NREINJ+iSprayInjNum(INI)
        NPMX=MAX(NPMX,iSprayInjNum(INI))
        enddo
!
        ALLOCATE(NP_REINJ(NREINJ),NP_DEL(NREINJ),
     &    WK(2,NPMX),W1K(NPMX),
     &    HHIST(2,NPMX),HIST(NREINJ))
        ALLOCATE(NRES(0:N_injctr),NRER(N_injctr)) 
! ---   
        NP_REINJ(:)=0
        NP_DEL(:)=0
        NRES(:)=0
        NRER(:)=0
!
        do INI=1,N_injctr 
        NRES(INI)=iSprayInjNum(INI) 
        enddo
!
        do INI=1,N_injctr
        NRES(INI)=NRES(INI-1)+iSprayInjNum(INI) 
        enddo
!
        ALLOCATE(JPASS(MEVENT))
        ALLOCATE(BCKND(MXCVFAC_P))
        ALLOCATE(JUDGE(MXALLCV_P))
        allocate(IPCPU(MXPRT)) 
        ALLOCATE(NPNEW(N_injctr),NP0(N_injctr),
     &        NDEL(N_injctr),NDEL1(N_injctr))
        ALLOCATE(NPLEFT(N_injctr))
        ALLOCATE(NSTA(N_injctr))
!
        ALLOCATE(PRT_VF(MXALLCV_P),stat=ierr1) 
        if(ierr1/=0) call 
     &   FFRABORT(1,'allocating PRT_VF error in FFR_TRACING')
        PRT_VF(:)=0.d0
!
        IPCPU(:)=-1
        TRESIDUAL=0.d0
        JPASS=0
        JUDGE=0
        JOUT=HPC_dead
        HIS_P=0.d0
        NPNEW=0
        NP0=0
        NDEL=0
        NDEL1=0
        NPLEFT=0
        NSTA=0
        NP_REINJ=0
        NRER=0
        BCKND=0
        PYS=0.d0
!
        iflag=1
!
      endif    !if(iflag==0) then
!
! --- save 
!
      dTSave=DT0
      PI=1.0D0
      PI=4.0*DATAN(PI)
      IE00=LVEDGE(1,MAT_CFIDX(0)+1)
!
!--------------------------- 
! --- initial condition 
!--------------------------- 
!
      IF(IT==ITS+1.or.i_initial==0) THEN
        i_initial=1 
!---------------------------------------------------------------------
        if(my_rank.eq.root)
     +  write(ifll,*) "MSG:[FFR_TRACING]:Reading particle.ctl..."
!
        ALLOCATE(dPDstion(MXPRT))
        ALLOCATE(dPDstvel(MXPRT))
        dPDstion=0.d0
        dPDstvel=0.d0
!
!---------------------------------------------------------------------
! --- Make Particle Initial Condition or Read from Initial File.
!---------------------------------------------------------------------
!
        write(ifll,*)
     +    "MSG:[FFR_TRACING]:Setting up spray initial condition..."
!
        CALL READ_PARTICLE_INI (NPSUM,NCV,MXALLCV,
     &    MXPRT,PI,MXCOMP,MXCVFAC,MXSSFBC,
     &    SFCENT,LBC_SSF,SFAREA,LVEDGE,CVCENT,
     *    DIA,XYZP,UVWP,RHOP,PARCEL,DT0,TMPP,PYS,
     &    PMASS,INSIDE,JOUT,HIS_P)
!
!------------------------------------------------- 
! --- Broadcast Initial Conditon to Other Proc. 
!------------------------------------------------- 
! 
        if(my_rank==root) 
     +    write(ifll,*) 
     +     "MSG:[FFR_TRACING]:Broadcasting initial condition..." 
!-------------------------------------------------------------
!
        NPALL(:)=0
!
        ITER_P=ITS
!
!
!------------------------------------------------------------
! --- GENERATE NEIGHBOR(MEP,MXALLCV), NVCRSIGN(MEP,MXALLCV) 
!------------------------------------------------------------
!
        if(iters.gt.total_start) then 
          CALL read_particle_restart(IT,time,
     &    NPALL,INSIDE,JOUT,NINALL,NSTA,NRES,
     &    PMASS,DIA,XYZP,UVWP,PARCEL,RHOP,TMPP,
     &    HEIGHT,FU,FV,FW,PYS(:,:,1),HIS_P,
     &    ierror)
        endif 
!
        CALL FFR_TRACING01(NMAT,MXALLCV,MXCVFAC,MEP,
     1           MAT_NO,MAT_CFIDX,MAT_CAL,MAT_CVEXT, 
     2           LVEDGE,SFAREA,SFCENT,CVCENT, 
     3           VCTRP,NEIGHBOR,DNORM,IERROR,BCKND,
     4           MXSSFBC,MXCVFAC_P,LBC_SSF) 
!
        if(my_rank.eq.root) write(ifll,*) 
     +    "MSG:[FFR_TRACING]:End of Initialization."
!
        NEIBCELL(:)=0
        IPCELL(:)=0
!
        IF(NPE>1) THEN 
          CALL PARTICLE_SEND_RECV01(MXALLCV,NCV,NEIBCELL,IPCELL)
        ENDIF
!
      ENDIF  !IF(IT==ITS+1) 
!
      ITER_P=ITER_P+1
      if(IT<total_start) return
!
      do IPCV=1,MXALLCV
      if(IPCELL(IPCV)>NEIBPETOT) then
        write(*,*) 'ERR: IPCELL(IPCV)= ',IPCELL(IPCV)
        call FFRABORT(1,'ERR-2: JOUT(IP)>NEIBPETOT')
      endif 
      enddo 
!
      if(reaction.and.ical_P_R==1) then
        ALLOCATE(tempp(MXPRT,2),stat=ierr1) 
        if(ierr1/=0) call 
     &    FFRABORT(1,'allocating tempp error in FFR_TRACING')
        ALLOCATE(WDOTP(MXCV,MXcomp),stat=ierr1) 
        if(ierr1/=0) call 
     &    FFRABORT(1,'allocating WDOTP error in FFR_TRACING')
        WDOTP(:,:)=0.d0
      endif
!
!
!---------------------------------------------------
! --- 
!---------------------------------------------------
!      if(ical_t) then 
!        if(sw) then 
!          call cal_t2cp(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
!     &          tmp(:,1),yys(:,:,1),pp0(:,1),rho(:,1),
!     &           MASSYS,cp,MASSYS)
!        else
!          call cal_t2hcp(1,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
!     &         tmp(:,1),yys(:,:,1),hhs,MASSYS,cofd,cp,MASSYS)
!        endif
!      endif
!
!----------------------------------------
! --- particle reaction flag: prt_reac
!----------------------------------------
!
      pat_reac=
     1         IFLAG_COAL==icoal.and.
     2         (encomp>=1.or.ical_flmlt).and.
     3         IT.gt.igT_iter.and.
     5         ical_reac
      if(II_O2/=0.and.II_NO/=0.and.pat_reac) then
        MASSYS(:,:)=0.d0 
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.lt.0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        MASSYS(ICVL,II_O2)=YYS(ICVL,II_O2)  !(kg)
     &                 *RHO(ICVL,1)*CVVOLM(ICVL)
        MASSYS(ICVL,II_NO)=YYS(ICVL,II_NO)
     &                 *RHO(ICVL,1)*CVVOLM(ICVL)
        enddo
        enddo
      endif
!
!---------------------------------------------------
! --- initial end 
!---------------------------------------------------
!---------------------------------------------------
! ---
!---------------------------------------------------
!
      DO INB=1,VCTRP(IE00,0)
      ICFLL=VCTRP(IE00,INB)
      ICFL=abs(ICFLL)
      CC(INB)=-DNORM(ICFL)
      DO I=1,3
      CC(INB)=CC(INB)+SFAREA(I,ICFL)*CVCENT(I,IE00)
      ENDDO
      CC(INB)=CC(INB)*SIGN(1,ICFLL)
      ENDDO
! 
      DO INB=1,VCTRP(IE00,0)
      HEIGHT0(INB)=CC(INB)
      ENDDO
!
!
!----------------------------------------------
!
      IF(NPE>1) then  !zhang-1
        do INI=1,N_injctr
        IPS=NINALL(INI-1)+1
        IPE=NINALL(INI-1)+NPALL(INI) 
        do IP=IPS,IPE
        if(JOUT(IP)==0) then 
          IPCV=INSIDE(IP) 
          IF(IPCV==0.or.IPCV>MXALLCV) then 
            write(ifll,*) 'MSG--: ',IP,IPCV 
            call FFRABORT(1,'ERR:IPCV==0')
          endif
          IF(IPCELL(IPCV)>=1.or.IPCV>NCV) then 
!          IF(IPCELL(IPCV)>=1) then 
            JOUT(IP)=IPCELL(IPCV) 
            if(JOUT(IP)>NEIBPETOT) then
              write(ifll,*) 'MSG: ',IPCV,IP,JOUT(IP)
              call FFRABORT(1,'ERR-0: JOUT(IP)>NEIBPETOT')
            endif
          endif
        elseif(JOUT(IP)>NEIBPETOT) then
          write(ifll,*) 'ERR: at:IP=, JOUT(IP)',IP,JOUT(IP)
          write(ifll,*) 'MSG: call your FFR supporter'
          call FFRABORT(1,'ERR-1: JOUT(IP)>NEIBPETOT')
        endif
        enddo
        enddo
      endif
!
!---------------------------------------------------------------------------
! --- CALCULATE THE VELOCITY GRADIENT, PRESSURE GRADIENT ... AT CELL CENTER 
!---------------------------------------------------------------------------
!
      CALL GRAD_CELL(3,10,
     &  MAT_CV,MAT_CVEXT,MAT_NO,MAT_CAL,MAT_CFIDX,vctrp,
     &  LVEDGE,SFAREA,SFCENT,WIFACE,CVVOLM,VEL(:,:,2),GRDC)
      NP0(:)=NPALL(:)
!------------------------
! --- 
!------------------------
      do INI=1,N_injctr 
      IF(ITER_P==iSprayInjStart(INI)) then
        starttime(INI)=time
      endif
      enddo
!
!--------------------- ---------------------------
! --- WRITE(*,*) "CALL INITIALIZATION  PARTICLES "
!--------------------- ---------------------------
!
      NPNEW(:)=0
      INJ_FLAG(:)=.false.
      do INI=1,N_injctr
      INJ_FLAG(INI)=
     & mod((ITER_P-iSprayInjStart(INI)),iSprayInjInter(INI))==0
      IF((ITER_P.GE.iSprayInjStart(INI)).and.
     &   (ITER_P.LE.iSprayInjEnd(INI)).and.
     &   INJ_FLAG(INI)
     &    .and.(NINALL(INI-1)+NPALL(INI)<NINALL(INI))
     &   ) THEN
        NPNEW(INI)=iSprayInjNum(INI)
      ENDIF
      NSTA(INI)=NSTA(INI)+NPNEW(INI)
      enddo
      
!
      if(.false.) then
      do INI=1,N_injctr
      IPS=NINALL(INI-1)
      DO IP=1,NPNEW(INI)
        NSTA1=IPS+NSTA(INI)+IP
        NP01=IPS+NP0(INI)+IP
        PARCEL(NP01)=PARCEL(NSTA1)
        DIA(NP01)=DIA(NSTA1)
        XYZP(NP01,1)=XYZP(NSTA1,1)
        XYZP(NP01,2)=XYZP(NSTA1,2)
        XYZP(NP01,3)=XYZP(NSTA1,3)
        UVWP(NP01,1)=UVWP(NSTA1,1)
        UVWP(NP01,2)=UVWP(NSTA1,2)
        UVWP(NP01,3)=UVWP(NSTA1,3)
        PMASS(NP01,1)=PMASS(NSTA1,1)
        PMASS(NP01,2)=PMASS(NSTA1,2)
        RHOP(NP01)=RHOP(NSTA1)
        TMPP(NP01,1)=TMPP(NSTA1,1)
        TMPP(NP01,2)=TMPP(NSTA1,2)
        PYS(NP01,:,1)=PYS(NSTA1,:,1)
        PYS(NP01,:,2)=PYS(NSTA1,:,2)
        DHDTP(NP01,1)=DHDTP(NSTA1,1)
        DHDTP(NP01,2)=DHDTP(NSTA1,2)
        DYDTP(NP01,:)=DYDTP(NSTA1,:)
        INSIDE(NP01)=INSIDE(NSTA1)
        HEIGHT(:,NP01)=HEIGHT(:,NSTA1)
!
!---------------------------------------------
! ++ Set Initial Value to Breakup Params ++ 
!---------------------------------------------
	dPDstion(NP01)=0.d0
	dPDstvel(NP01)=0.d0
      ENDDO
      enddo
      endif
!----------------------------
! --- 
!----------------------------

      do INI=1,N_injctr 
      NPALL(INI)=NP0(INI)+NPNEW(INI) 
      INJ_FLAG(INI)= 
     &  INJ_FLAG(INI).and.NPALL(INI)
     &  >=NINALL(INI)-NINALL(INI-1)
      NPALL(INI)=MIN(NPALL(INI),
     &           NINALL(INI)-NINALL(INI-1)) 
      enddo
!---------------------------
! --- 
!---------------------------
      do INI=1,N_injctr 
      IPS=NINALL(INI-1)+1+NP0(INI)
      IPE=NINALL(INI-1)+NPALL(INI)
      do IP=IPS,IPE
      IF(JOUT(IP)==HPC_dead) cycle
      IE00=min(max(INSIDE(IP),IE00),NCV)
      CALL LOCATE_PARTICLE_TO_CELL(1,IP,INI,
     +            INSIDE,JOUT,VCTRP,XYZP,DIA,RHOP,HEIGHT,UVWP,
     +            HEIGHT0,CVCENT,SFAREA,SFCENT,NEIGHBOR,IE00,
     +            DNORM,IPCELL,BCKND,LVEDGE,WDRBND,PARCEL,DT0)
      if(src_p_ini==usryes) then
        T_U=TMPP(IP,1)
        rcomp(:)=PYS(IP,:,1)
        DIM_U=DIA(IP)
        RHO_U=RHOP(IP)
        XYZP0(1)=XYZP(IP,1)
        XYZP0(2)=XYZP(IP,2)
        XYZP0(3)=XYZP(IP,3)
        UVWP0(1)=UVWP(IP,1)
        UVWP0(2)=UVWP(IP,2)
        UVWP0(3)=UVWP(IP,3)
        call USER_ini_particle
     &      (IP,INI,IT,time,dTSave,ncomp,
     &      DIM_U,RHO_U,XYZP0,UVWP0,rcomp,T_U)
        XYZP(IP,1)=XYZP0(1)
        XYZP(IP,2)=XYZP0(2)
        XYZP(IP,3)=XYZP0(3)
        UVWP(IP,1)=UVWP0(1)
        UVWP(IP,2)=UVWP0(2)
        UVWP(IP,3)=UVWP0(3)
        TMPP(IP,1)=T_U
        PYS(IP,:,1)=rcomp(:)
        RHOP(IP)=RHO_U
        DIA(IP)=DIM_U
      endif
      enddo
      enddo
!
!----------------------------------------------
      do 1231 INI=1,N_injctr 
        nb=bc_nb(INI)                 !       NO_BC(INI)
        if(NRER(INI)<=0) cycle 
        IPPP=1
        NNN=0
        do 1220
        NNN=NNN+1
        if(NNN>=10000) call FFRABORT(1,'ERR: NO LIMIT DO loop')
        N0=0
        N00=0
        do 1230 ICPU=1,NPE
        IDUM=IBC_P_NO(ICPU,INI,nb,2)
        N0=N0+IDUM
        if(IDUM==0) cycle
        if(NPINJ(ICPU,INI)>IDUM) cycle
        N00=N00+1
        if(IPPP>NRER(INI)) goto 1210
        IP=NP_REINJ(NRES(INI-1)+IPPP)
        
        NPINJ(ICPU,INI)=NPINJ(ICPU,INI)+1
        if(IP==0) then
          cycle
        endif
        IPCPU(IP)=ICPU-1
        
        IPPP=IPPP+1
 1230   enddo
        if(N0==0.or.N00==0) goto 1210
 1220   enddo
 1210   continue
 1231 enddo
!
      do INI=1,N_injctr 
        if(NRER(INI)<=0) cycle 
!        if(injkd(INI)/=i_coal.and.injkd(INI)/=inlet) cycle !zhang-1
        if(NPE>1) then  !zhang-1
          DO ICPU=1,NPE
          IDUM=NPINJ(ICPU,INI)
          call hpcisum(IDUM)
          NPINJ(ICPU,INI)=IDUM
          enddo
        endif
        nb=bc_nb(INI)
        DO ICPU=1,NPE
        if(NPINJ(ICPU,INI)>=IBC_P_TOL(INI,nb,1)) then
           NPINJ(ICPU,INI)=0
        endif
        enddo
      enddo
! --- 

      do INI=1,N_injctr 
        if(NRER(INI)>0) then
          do IPP=1,NRER(INI)
          IP=NP_REINJ(NRES(INI-1)+IPP) 
          IF(IP==0) cycle
          if(JOUT(IP)>=0) cycle 

          if(IPCPU(IP)/=my_rank)  cycle 
          JOUT(IP)=0
!
          call RE_INJ_P(IP,INI,MXPRT,PI,MXCOMP,MXCVFAC,MXSSFBC,
     &      SFCENT,LBC_SSF,SFAREA,JOUT,IPCPU,
     *      DIA,XYZP,UVWP,RHOP,PARCEL,NDIS,dTSave,TMPP,PYS,
     &      PMASS,HIS_P,INSIDE,LVEDGE)
!
          IE00=min(INSIDE(IP),NCV)
          CALL LOCATE_PARTICLE_TO_CELL(2,IP,INI,
     +            INSIDE,JOUT,VCTRP,XYZP,DIA,RHOP,HEIGHT,UVWP,
     +            HEIGHT0,CVCENT,SFAREA,SFCENT,NEIGHBOR,IE00,
     +            DNORM,IPCELL,BCKND,LVEDGE,WDRBND,PARCEL,DT0)
          if(src_p_ini==usryes) then 
            rcomp(:)=PYS(IP,:,1)
            T_U=TMPP(IP,1)
            DIM_U=DIA(IP)
            RHO_U=RHOP(IP)
            XYZP0(1)=XYZP(IP,1)
            XYZP0(2)=XYZP(IP,2)
            XYZP0(3)=XYZP(IP,3)
            UVWP0(1)=UVWP(IP,1)
            UVWP0(2)=UVWP(IP,2)
            UVWP0(3)=UVWP(IP,3)
            call USER_ini_particle
     &      (IP,INI,IT,time,dTSave,ncomp,
     &      DIM_U,RHO_U,XYZP0,UVWP0,rcomp,T_U)
            XYZP(IP,1)=XYZP0(1)
            XYZP(IP,2)=XYZP0(2)
            XYZP(IP,3)=XYZP0(3)
            UVWP(IP,1)=UVWP0(1)
            UVWP(IP,2)=UVWP0(2)
            UVWP(IP,3)=UVWP0(3)
            RHOP(IP)=RHO_U
            DIA(IP)=DIM_U
            TMPP(IP,1)=T_U
          endif
          enddo
        endif
      enddo
!     
! --- 
!     
      do INI=1,N_injctr
      IPS=NINALL(INI-1)!+1
      IPE=NINALL(INI-1)+NPALL(INI)
      DO IP=IPS+NP0(INI)+1,IPE
!      if(JOUT(IP)>=0.or.JOUT(IP)==HPC_dead) cycle
      FU(IP,2)=0.d0
      FV(IP,2)=0.d0
      FW(IP,2)=0.d0
      ENDDO
      enddo
C----------------------------------------------------
C --- WRITE(*,*) "END OF INITAILIZATION PARTICLES" 
C----------------------------------------------------
      do INI=1,N_injctr 
      NUMINJCT=0
      if(NPNEW(INI)>0) THEN 
        NUMINJCT=NUMINJCT+NPNEW(INI)   !NRER(INI)
      endif
      if(NRER(INI)>0) then
        NUMINJCT=NUMINJCT+NRER(INI)
      endif 
!!!!      NUMINJCT=NUMINJCT-NDEL(INI)
!
      if(NPE>1) call hpcisum(NUMINJCT)   !zhang-1   
      
      if(mod(IT,monitor_stp)==0) then
        if(my_rank.eq.root)
     &    write(ifll,fmt="(3X,'MSG:Injector=',I4,1X,I8,1X,E15.6,
     &   ' (SMD) particles were injected.')")
     &  INI,NUMINJCT/NPE,dSMD(INI)
      endif
      enddo

!
!----------------------------------------------------
! --- TRACE PARTICLES IN THE FLOW FIELD 
!----------------------------------------------------
!
      dMinRelaxTime=1.d0
!      dMaxPtReynolds=0.d0
      iMaxPassedCell=0
      DT_IN=DT0/DBLE(ITER_INNER)
      iBreakupCount=0
!
      if(mod(IT,monitor_stp)==0) then
        if(my_rank.eq.root) then
        write(ifll,fmt="(4X,25('+'),' Particle caluclation : 
     &   ',I4,25('+'))")
     +        ITER_P
!        
!        write(ifll,fmt="(10X,'Fluid ITER_P=',I4)") ITER_P
        write(ifll,fmt="(10X,'Fluid ITER_P=',I6)") ITER_P
        write(ifll,fmt="(10X,'Fluid Time Step=',E12.5)") DT0
        write(ifll,fmt="(10X,'Particle Time Step=',E12.5)") DT_IN
      endif
      endif
!
      MOMCELL(:,:)=0.d0
!
!
      dMaxRandomWalk=0.d0
      dAveRandomWalk=0.d0
      iTotalCountOfRW=0
!
!
!---------------------------------------------------------
! --- Sub cycle time-step 
!---------------------------------------------------------
      if(radflag==1) then
        P_rad(:,1)=0.d0
        P_rad(:,2)=0.d0
      endif
!
      PRT_VF(:)=0.d0
!
!
!
      DO 1000 ITER_IN_COUNT=1,ITER_INNER 
!
!
!      
      PMASS(:,2)=PMASS(:,1)
      PYS(:,:,2)=PYS(:,:,1)
      TMPP(:,2)=TMPP(:,1)
!
      JPASS(:)=0
      JUDGE(:)=0
!
!
      do 1500 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      if(ifuel(INI)==iliquid) then
        DO 200 IP=IPS,IPE
        if(JOUT(IP)/=0) cycle
        DT0=DT_IN
        lBreakup=.FALSE.
!-----------------------------
! --- Breakup Calculation    +
!-----------------------------
        if(.not.lBreakupNONE(INI)) THEN
          if(lBreakupNTAB(INI).or.lBreakupITAB(INI)) THEN
            if(lBreakupNTAB(INI)) THEN
              UVWP0(1)=UVWP(IP,1)
              UVWP0(2)=UVWP(IP,2)
              UVWP0(3)=UVWP(IP,3)
              CALL BREAKUP_NORMAL_TAB(INI,IT,ITS,
     1        DT0,GRDC,IP,INSIDE,
     2        RHOP,DIA,XYZP,UVWP0,
     3        PARCEL,VEL,RHO,RMU,CVCENT,lBreakup,dTbu,dDiamReduc, 
     4        CVVOLM,aks)
              if(lBreakup.and..false.) then
                JPASS(:)=0
                JUDGE(:)=0
                KE=INSIDE(IP)
                KE0=KE
                UVWP0(1)=UVWP(IP,1)
                UVWP0(2)=UVWP(IP,2)
                UVWP0(3)=UVWP(IP,3)
                DO INB=1,VCTRP(KE,0)
                HEIGHT1(INB)=HEIGHT(INB,IP)
                ENDDO
                FU(IP,1)=FU(IP,2)
                FV(IP,1)=FV(IP,2)
                FW(IP,1)=FW(IP,2)
!
                CALL FFR_TRACING02_F(INI,
     *          IP,XYZP,UVWP,HEIGHT,FU,FV,FW,RHOP,
     *          INSIDE,KE,KE0,UVWP0,VCTRP,HEIGHT1,SFAREA,SFCENT,
     *          dTbu,GRDC,CVCENT,DIA,
     *          VEL,RHO,RMU,NEIGHBOR,DNORM,BCKND,
     *          JPASS,JUDGE,JOUT,MOMCELL,PMASS(:,1),PARCEL,
     *          IPCELL,TRESIDUAL,CVVOLM,aks,LVEDGE,WDRBND)
                DT0=DT0-dTbu
              endif
            else
              UVWP0(1)=UVWP(IP,1)
              UVWP0(2)=UVWP(IP,2)
              UVWP0(3)=UVWP(IP,3)
              CALL BREAKUP_IMPROVED_TAB
     &        (IT,ITS,INI,
     &        DT0,GRDC,IP,INSIDE,
     &        RHOP,DIA,XYZP,UVWP0,
     &        PARCEL,VEL,RHO,RMU,CVCENT,lBreakup,dTbu,dDiamReduc,
     &        CVVOLM,aks)
            endif
!
            if(lBreakup) then
!----------------------
! Change Diameter
!----------------------
              call utl_random(dTemp)
              CALL NUKIYAMA_TANAZAWA_VOLUMEBASE(INI,dDIANEW,
     +               DIA(IP)/dDiamReduc,dTemp)
!
              if(dDIANEW .gt. DIA(IP)) dDIANEW = DIA(IP)
              dVolume=DIA(IP)*DIA(IP)*DIA(IP)*PARCEL(IP)
              PARCEL(IP)=dVolume/(dDIANEW*dDIANEW*dDIANEW)
              iBreakupCount = iBreakupCount+1
!
              DIA(IP)=dDIANEW
              PMASS(IP,1)=PI*DIA(IP)*DIA(IP)*DIA(IP)*RHOP(IP)/6.0
              dPDstion(IP)=0.d0
              dPDstvel(IP)=0.d0
!--------------------------------------------------
! If N-th particle goes out, skip the calculation.
!--------------------------------------------------
              if(JOUT(IP)>0) THEN
                TRESIDUAL(IP)=TRESIDUAL(IP)+DT0
                goto 200
              ENDIF
            endif
!
          else if(lBreakupPILCH(INI)) THEN
            UVWP0(1)=UVWP(IP,1)
            UVWP0(2)=UVWP(IP,2)
            UVWP0(3)=UVWP(IP,3)
            CALL BREAKUP_PILCH(INI,IT,ITS,
     &         DT0,GRDC,N,INSIDE,
     2         RHOP,DIA,XYZP,UVWP0,
     3         PARCEL,VEL,RHO,RMU,CVCENT,lBreakup,dTbu,dDiamReduc,
     4         CVVOLM,aks)
!
            if(lBreakup) then
              dDIANEW=DIA(IP)/dDiamReduc
              dVolume=DIA(IP)*DIA(IP)*DIA(IP)*PARCEL(IP)
              PARCEL(IP)=dVolume/(dDIANEW*dDIANEW*dDIANEW)
              DIA(IP)=dDIANEW
              PMASS(IP,1)=PI*DIA(IP)*DIA(IP)*DIA(IP)*RHOP(IP)/6.d0
              iBreakupCount=iBreakupCount+1
            endif
          endif
        endif
!
!++++++++++++++++++++++++++ Breakup End +++++++++++++++++++++++++
!
        KE=INSIDE(IP)
        KE0=KE
        UVWP0(1)=UVWP(IP,1)
        UVWP0(2)=UVWP(IP,2)
        UVWP0(3)=UVWP(IP,3)
        DO K1=1,VCTRP(KE,0)
        HEIGHT1(K1)=HEIGHT(K1,IP)
        ENDDO
        FU(IP,1)=FU(IP,2)
        FV(IP,1)=FV(IP,2)
        FW(IP,1)=FW(IP,2)
!
        CALL FFR_TRACING02_F(INI,
     *  IP,XYZP,UVWP,HEIGHT,FU,FV,FW,RHOP,
     *  INSIDE,KE,KE0,UVWP0,VCTRP,HEIGHT1,SFAREA,SFCENT,
     *  DT0,GRDC,CVCENT,DIA,
     *  VEL,RHO,RMU,NEIGHBOR,DNORM,BCKND,
     *  JPASS,JUDGE,JOUT,MOMCELL,PMASS(:,1),PARCEL,
     *  IPCELL,TRESIDUAL,CVVOLM,aks,LVEDGE,WDRBND)
!
 200    enddo
      elseif(ifuel(INI)==icoal.or.
     &       ifuel(INI)==iglass.or.
     &       ifuel(INI)==isolid) then 

        DO 210 IP=IPS,IPE
        if(JOUT(IP)/=0) cycle
        DT0=DT_IN
        KE=INSIDE(IP)
        KE0=KE
        UVWP0(1)=UVWP(IP,1)
        UVWP0(2)=UVWP(IP,2)
        UVWP0(3)=UVWP(IP,3)
        DO K1=1,VCTRP(KE,0)
        HEIGHT1(K1)=HEIGHT(K1,IP)
        ENDDO
        FU(IP,1)=FU(IP,2)
        FV(IP,1)=FV(IP,2)
        FW(IP,1)=FW(IP,2)
!
        CALL FFR_TRACING02_F(INI,
     *  IP,XYZP,UVWP,HEIGHT,FU,FV,FW,RHOP,
     *  INSIDE,KE,KE0,UVWP0,VCTRP,HEIGHT1,SFAREA,SFCENT,
     *  DT0,GRDC,CVCENT,DIA,
     *  VEL,RHO,RMU,NEIGHBOR,DNORM,BCKND,
     *  JPASS,JUDGE,JOUT,MOMCELL,PMASS(:,1),PARCEL,
     *  IPCELL,TRESIDUAL,CVVOLM,aks,LVEDGE,WDRBND)
 210    enddo
      endif
 1500 enddo
!-----------------------------------------
!     HPC - MPI 
!-----------------------------------------
!    WHEN A PARTICLE GOES INTO AN OUTER CELL, 
!        IT WILL BE SENT TO THE NEIGHBOR DOMAIN 
!    WHERE THE CELL IS A INNER CELL.
!-----------------------------------------
!
!
! --- 
!
      IF(NPE>1) THEN !zhang-
        IPSND(:)=0
        I_HPC=0
        NPEXCHAGOUT(:)=0
        NPEXCHAGIN(:)=0
        NPOUTALL=0
        do INI=1,N_injctr
        IPS=NINALL(INI-1)+1
        IPE=NINALL(INI-1)+NPALL(INI)
        DO IP=IPS,IPE

        if(JOUT(IP)<=0) cycle
        NCPU=JOUT(IP)           !IPCELL(JOUT(IP))
        
        IF(NCPU>0) THEN 
          NPEXCHAGOUT(NCPU)=NPEXCHAGOUT(NCPU)+1
          NPOUTALL=NPOUTALL+1
          IPSND(NPOUTALL)=IP
          IPQUENE(NPOUTALL)=NPEXCHAGOUT(NCPU)
          I_HPC=I_HPC+1
        ENDIF
        ENDDO
        enddo
! --- 
        CALL PARTICLE_SEND_RECV02(NPEXCHAGIN,
     &       NPEXCHAGOUT) 
!---------------------------------
! --- 
!---------------------------------
        NPINALL=0 
        DO NCPU=1,NEIBPETOT 
          NPINALL=NPINALL+NPEXCHAGIN(NCPU)
        ENDDO
!
!----------------------------------------------
! --- 
!----------------------------------------------
!
        CALL hpcimax(I_HPC)
!
        if(I_HPC>0) then  !zhang-1
          IMODE=1
          CALL PARTICLE_MPI(IMODE,
     &       MXPRT,MXALLCV_P,
     &       MXCVFAC_P,N_injctr,
     &       XYZP,UVWP,DIA,TMPP,PYS,
     *       FU,FV,FW,RHOP,PARCEL,PMASS(:,1),PMASSC,INSIDE,
     *       NEIBCELL,TRESIDUAL,
     *       NPEXCHAGOUT,NPEXCHAGIN,IPSND,IPQUENE,
     &       dPDstion,dPDstvel,
     &       NINALL,JOUT,IPCELL,!IPCPU,
     6       NPOUTALL,NPINALL
     &  )
        endif

C-------------------------------------------------------
C     TREAT THE BOUNDARY CONDITION (JOUT = -1) AGAIN
C     IGNORE COLLISION EFFECT
C-------------------------------------------------------
      endif
!
!------------------------------------
! --- 
!------------------------------------
      DYDTP(:,:)=0.d0
      DHDTP(:,:)=0.d0
      if(reaction.and.ical_P_R==1) then
        WDOTP(:,:)=0.d0
      endif
!----------------------------------
! --- Spray liquid evapration
!----------------------------------
      if(
     &   iflag_eva==1.and.
     &   ((encomp>=1.or.ical_flmlt).or.(calxi.and.idifflm/=ORG))
     &   ) then
      endif
!------------------------------------
! --- CHO, HCN & H2O evaporation 
!------------------------------------
      if(
     &  (encomp>=1.or.ical_flmlt).and.
     &  IFLAG_COAL==icoal.and.
     &  ical_reac.and.
     &  IT.gt.igT_iter
     &  ) then
      endif
!----------------------------------------
! --- 
!----------------------------------------
      if(reaction.and.ical_P_R==1) then
      endif
!--------------------------------------------
! --- Char (Solid C) combustion (2C+O2=2CO)
! --- Char-NO :2N+O2=2NO
!--------------------------------------------
!      if(
!     &  encomp>=1.and.
!     &  IFLAG_COAL==icoal.and.
!     &  ical_reac.and.
!     &  IT.gt.igT_iter
!     &  ) then
!         call Combustion_Char( 
!     &  IT,ITER_P,time,DT_IN,N_injctr,PI,
!     &  NINALL,INSIDE,JOUT,PARCEL,DIA,TMPP,NPALL,      !WK1_P,
!     &  DYDTP,DHDTP,PMASS,PYS,RHO,
!     &  WDOTP,tmp,yys,EVAP_Y,EVAP_H,CVVOLM,MASSYS,
!     &  ierror)
!      endif
!
!---------------------------------
! --- FUEL-NO: 2HCN+O2=>2NO
!              2HCN+2NO=>2N2+O2
!---------------------------------
!
      if(
     &  (encomp>=1.or.ical_flmlt).and.
     &  IFLAG_COAL==icoal.and.
     &  ical_reac.and.
     &  IT.gt.igT_iter
     &  ) then
        if(IT>=Fuel_NO_P) then
        endif
      endif 
!---------------------------------
! --- Prompt NO 
!---------------------------------
      if(
     &  (encomp>=1.or.ical_flmlt).and.
     &  ical_reac.and.
     &  it.gt.igT_iter.and.
     &  IFLAG_COAL==icoal
     &  ) then
        if(it>=prompt_NO_P) then
        endif
        if(it>=Zeldovich_P) then
        endif
      endif
!--------------------------------------------
! --- enthalpy transport & particle Temp
!--------------------------------------------
      if(ical_t) then
      endif
!
      if(ical_t) then
        QRAD=0.d0
        dum3=0.d0
        DT0=DT_IN
        do 3000 INI=1,N_injctr
        IPS=NINALL(INI-1)+1
        IPE=NINALL(INI-1)+NPALL(INI)
! --- Glass only
        if(ifuel(INI)==iglass) then
          DO IP=IPS,IPE
            if(JOUT(IP)/=0) cycle
            IPCV=INSIDE(IP)
            if(IPCV>NCV) cycle
            DP_M=DIA(IP)
            TP_M=TMPP(IP,1)
            RHOP_M=RHOP(IP)
            MASSP_M=PMASS(IP,1)
            VELP_M(1:3)=UVWP(IP,1:3)
!
            G_T_M=TMP(IPCV)
            G_VEL_M(1:3)=VEL(IPCV,1:3,1)
            rmu_G_M=RMU(IPCV)
            G_RHO_M=RHO(IPCV,1)
            G_RAMD_M=0.03   !rmdc(IPCV)
            G_CP_M=1000.0  !cp(IPCV)
            call USER_PARTICLE_TEMP1(MXPRT,IP,N_injctr,DT0,
     &           DP_M,TP_M,RHOP_M,MASSP_M,VELP_M,
     &           TMPP(:,1),DIA,
     &           G_T_M,G_VEL_M,rmu_G_M,G_RHO_M,G_RAMD_M,G_CP_M)
            RHOP(IP)=RHOP_M
            PMASS(IP,1)=MASSP_M
            TMPP(IP,1)=TP_M
            DIA(IP)=DP_M
          enddo
          cycle
        endif
! --- Glass only end
 3000   enddo
      endif
!--------------------------
! --- 
!--------------------------
      if(encomp>=1.or.ical_flmlt) then
        dum2=0.d0
        do 4000 INI=1,N_injctr
        IPS=NINALL(INI-1)+1
        IPE=NINALL(INI-1)+NPALL(INI) 
        if(ifuel(INI)/=icoal.and.
     &     ifuel(INI)/=iglass.and.
     &     ifuel(INI)/=isolid.and.
     &     ifuel(INI)/=iliquid
     &     ) cycle  
        DT0=DT_IN 
!
        DO 4100 IP=IPS,IPE 
        if(JOUT(IP)/=0) cycle
        do ICOM=1,NCOMP
        PYS(IP,ICOM,1)=max(0.d0,
     &  (PYS(IP,ICOM,1)*PMASS(IP,1)+DT0*DYDTP(IP,icom)))!/PMASS(IP,1))
        enddo
 4100   enddo
!
        DO 4200 IP=IPS,IPE
        if(JOUT(IP)/=0) cycle
        dum1=0.d0
        do ICOM=1,NCOMP
        dum1=dum1+PYS(IP,ICOM,1)
        enddo
        PMASS(IP,1)=dum1
        PYS(IP,:,1)=PYS(IP,:,1)/(dum1+SML)
 4200   enddo
!
 4000   enddo

!        if(ITER_IN_COUNT==1) then
!        print*,'XXXXXXXXXX-1',PMASS(1:2,1)
!        print*,'XXXXXXXXXX-2',DIA(1:2)
!        endif
!---------------------------------
! --- UPDATE MASS & DIA 
!---------------------------------
!        IF(IFLAG_COAL==icoal) then
         call UPDATE_P(
     &   IT,ITER_P,time,DT_IN,DIMV,N_injctr,PI,
     &   NINALL,JOUT,NPALL,INSIDE,DYDTP,PARCEL,
     &   PMASS,DIA,RHOP,HIS_P,XYZP,PRT_VF)
        ENDIF
!      endif
!
 1000 enddo
!
      
!      print*,'particle release mass=',mass_p_total(1)-PMASS(1,1)*PARCEL(1) 
      
! ----------------------------------------------------
!
!      if(
!!     &  IFLAG_COAL==icoal.and..false.  !1111
!     &  IFLAG_COAL==icoal.or.
!     &  ) then 
!
      reinj=0
      lack_p(:)=0
      re_inj(:)=0
      do INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      DO IP=IPS,IPE
      if(JOUT(IP)==0) then
        lack_p(INI)=lack_p(INI)+1
      endif
      enddo
      enddo
!
      do INI=1,N_injctr
      if(lack_p(INI)>iSprayInjTotal(INI)) then
        re_inj(INI)=1
        reinj=1
      endif
      enddo


      if(reinj==1) then
        NP_REINJ(:)=0
!-------------------------
! --- particle on BC 
!-------------------------
!
        do 5400 INI=1,N_injctr 
!        if(re_inj(INI)==0) cycle 
        NRER(INI)=0 
        
!        if(injkd(INI)==i_coal.or.injkd(INI)==inlet) then
        if(ifuel(INI)==icoal.and.
     &     ifuel(INI)==iglass.or.
     &     ifuel(INI)==isolid.or.
     &     ifuel(INI)==iliquid) then 
!
        if(.NOT.INJ_FLAG(INI)) cycle 
        if(ITER_P>iSprayInjEnd(INI)) cycle 
        IPS=NINALL(INI-1)+1
        IPE=NINALL(INI-1)+NPALL(INI)
        if(NINALL(INI-1)+NPALL(INI)<NINALL(INI)) cycle
        DO 5410 IP=IPS,IPE 
        if(NRER(INI)<iSprayInjNum(INI)) then
          if(JOUT(IP)==HPC_dead) cycle
          if(JOUT(IP)<0.or.JOUT(IP)==BC_dead 
     &    .or.JOUT(IP)==A_dead.or.JOUT(IP)==mass_dead) then 
            NRER(INI)=NRER(INI)+1
            NP_REINJ(NRES(INI-1)+NRER(INI))=IP
!            JOUT(IP)=P_Reced
          endif
        endif
 5410   enddo
        endif
 5400   enddo
!
!-------------------------
! --- mass zero 
!-------------------------
        do 5300 INI=1,N_injctr
!        if(re_inj(INI)==0) cycle 
        if(ifuel(INI)==icoal.or.
     &     ifuel(INI)==iglass.or.
     &     ifuel(INI)==isolid.or.
     &     ifuel(INI)==iliquid) then
          if(.NOT.INJ_FLAG(INI)) cycle
          if(ITER_P>iSprayInjEnd(INI)) cycle
          if(NRER(INI)>=iSprayInjNum(INI)) cycle 

          IPS=NINALL(INI-1)+1
          IPE=NINALL(INI-1)+NPALL(INI)
          if(NINALL(INI-1)+NPALL(INI)<NINALL(INI)) cycle 
          if(ifuel(INI)==icoal
     &    ) then
             DO IP=IPS,IPE
             if(JOUT(IP)==HPC_dead) cycle
             if(NRER(INI)<iSprayInjNum(INI)) then
               if(abs(1.d0-PYS(IP,II_ASH,1))<SML) then 
                 NRER(INI)=NRER(INI)+1
                 NP_REINJ(NRES(INI-1)+NRER(INI))=IP
               endif
             endif
             enddo
          elseif(ifuel(INI)==iliquid.or.
     &           ifuel(INI)==iglass.or.
     &           ifuel(INI)==isolid
     &       ) then
             if(IT<ical_evap(INI).or.ical_evap(INI)==0) cycle
             DO IP=IPS,IPE
             if(JOUT(IP)==HPC_dead) cycle
             if(NRER(INI)<iSprayInjNum(INI)) then
               dum1=0.D0
               dum2=0.D0
               DO ICOM=1,NCOMP
               DUM1=DUM1+PYS(IP,ICOM,1)
               DUM2=DUM2+PYS(IP,ICOM,1)*(1.D0-EVAP_COMP(ICOM,INI))
               ENDDO
               if(abs(1.d0-DUM2)<SML) then 
                 ! evap-comp = 0%
                 NRER(INI)=NRER(INI)+1
                 NP_REINJ(NRES(INI-1)+NRER(INI))=IP
               endif
             endif
             enddo
          endif
        else
          
        endif
 5300   enddo
!
!-------------------------
! --- History time 
!-------------------------
!
        NDEL1(:)=NRER(:)
!        
        NP_DEL(:)=0
        HIST(:)=0.d0
        HHIST(:,:)=0.d0
        do 5200 INI=1,N_injctr
        if(re_inj(INI)==0) cycle 
        NP_DEL(:)=0
!        if(ifuel(INI)/=icoal.and.ifuel(INI)/=iglass) cycle
!        if(injkd(INI)==i_coal.or.injkd(INI)==inlet) then
        if(.NOT.INJ_FLAG(INI)) cycle
        if(ITER_P>iSprayInjEnd(INI)) cycle
        do 5210 IPP=1,iSprayInjNum(INI)
        if(NRER(INI)<iSprayInjNum(INI)) then
          IPS=NINALL(INI-1)+1
          IPE=NINALL(INI-1)+NPALL(INI)
          dum1=-1.d10
          IPPP=0
          if(IPS<=IPE) then 
            DO 5220 IP=IPS,IPE 
            if(JOUT(IP)==HPC_dead) cycle
            dum2=HIS_P(IP)
            
!!!!!!!!!!!!!!!!!!!            if(JOUT(IP)==0) then
            if(JOUT(IP)==P_Reced.or.JOUT(IP)>0) cycle
              IF(dum2>dum1) then
                dum1=dum2
                IPPP=IP
              endif
!!!!!!            endif
 5220       enddo
!
            if(IPPP/=0) then
              NDEL1(INI)=NDEL1(INI)+1
              NP_DEL(NDEL1(INI))=JOUT(IPPP)
!
              NRER(INI)=NRER(INI)+1
              NP_REINJ(NRES(INI-1)+NRER(INI))=IPPP
              HIST(NRES(INI-1)+NRER(INI))=dum1
              JOUT(IPPP)=P_Reced 
            endif
          endif
        endif
 5210   enddo
! --- 
        if(NDEL1(INI)==0) cycle 
        I=NRES(INI-1)+NDEL1(INI)
        J=NRER(INI)
        N0=0
        DO K=I,J
         N0=N0+1
         JOUT(NP_REINJ(K))=NP_DEL(N0)
        enddo
!        endif
 5200   enddo 
!
        if(NPE>1) then
          do INI=1,N_injctr 
          idum=NRER(INI) 
          call hpcimax(idum) 
          NRER(INI)=idum
          enddo
        endif
!----------------------------------------------
        do INI=1,N_injctr 
        HHIST(:,:)=0.d0
        WK=0.d0
        W1K=-1
        if(.NOT.INJ_FLAG(INI)) cycle
        IMAT=NRES(INI-1)+NRER(INI) 
        IIMAT=NRES(INI-1)+1
        if(NRER(INI)==0) cycle 
        idum=NRER(INI)
        HHIST(1,1:idum)=HIST(IIMAT:IMAT)
        HHIST(2,1:idum)=NP_REINJ(IIMAT:IMAT)
        N=iSprayInjNum(INI) 
        call HPC_MAXLOC(N,HHIST(1:2,1:N),WK(1:2,1:N),W1K(1:N)) 
        NP_REINJ(IIMAT:IMAT)=0
        DO I=1,N
        if(W1K(I)<=0) cycle 
        NP_REINJ(NRES(INI-1)+I)=W1K(I)
!        JOUT(W1K(I))=P_Reced   !HPC_dead  !BABA  !p_n=11
        JOUT(W1K(I))=HPC_dead  !BABA    !!p_n=40
        ENDDO
        enddo
!
      endif
!---------------------------------------
! --- time-averaged rate for gas phase
!---------------------------------------
      do 7000 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT.lt.0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      DT0=dTSave
!
      if(reaction.and.ical_P_R==1) then
        do ICOM=1,ncomp
        do ICVL=ICVS,ICVE
        EVAP_Y(ICVL,ICOM)=EVAP_Y(ICVL,ICOM)/dTSave  !(kg/s/m^3)
        dum1=WDOTP(ICVL,ICOM)/dTSave    ! (kg/s/m^3)
        WDOT(ICVL,ICOM)=WDOT(ICVL,ICOM)+dum1
        enddo
        enddo
      else
        do ICOM=1,ncomp
        do ICVL=ICVS,ICVE
        EVAP_Y(ICVL,ICOM)=EVAP_Y(ICVL,ICOM)/dTSave  !(kg/s/m^3)
        enddo
        enddo
      endif
!
      do ICVL=ICVS,ICVE
      EVAP_H(ICVL,1)=EVAP_H(ICVL,1)/dTSave            ! (J/s/m^3)
      EVAP_H(ICVL,2)=EVAP_H(ICVL,2)/dTSave            ! (J/s/m^3)
      PRT_VF(ICVL)=PRT_VF(ICVL)/dTSave/CVVOLM(ICVL)
      ENDDO
      if(radflag==1) then
        do ICVL=ICVS,ICVE
        P_rad(ICVL,1)=P_rad(ICVL,1)/dTSave
        P_rad(ICVL,2)=P_rad(ICVL,2)/dTSave
        enddo
      endif
!
      do I=1,4
      do ICVL=ICVS,ICVE
      MOMCELL(ICVL,I)=MOMCELL(ICVL,I)/dTSave      !N
      enddo
      enddo
 7000 enddo
!
!------------------------------------
! --- output particle information 
!------------------------------------
      NPALLTOTAL=0
      DO INI=1,N_injctr
      NPALLTOTAL=NPALLTOTAL+NPALL(INI)
      enddo
!
      NpAllTotal2=0
      do ini=1,n_injctr
      IPS=NINALL(ini-1)+1
      IPE=NINALL(ini-1)+NPALL(ini)
      do IP=IPS,IPE
      if(JOUT(IP)==HPC_dead.or.JOUT(IP)==BC_dead) cycle
      NpAllTotal2=NpAllTotal2+1
      enddo
      enddo
!
      if(NPE>1) then   !zhang-1
        call hpcisum(NPALLTOTAL)
        call hpcisum(NPALLTOTAL2)
!
        if(NPALLTOTAL.gt.0) then
          call hpcrmin(dMinRelaxTime)
          call hpcimax(iMaxPassedCell)
          call hpcisum(iBreakupCount)
          call hpcrmax(dMaxRandomWalk)
          call hpcrsum(dAveRandomWalk)
          call hpcisum(iTotalCountOfRW)
!
          NPALLTOTAL=NPALLTOTAL/NPE
          if(my_rank .eq. root) then
            write(ifll,fmt="(' Total Parcels       :',I6)") NPALLTOTAL
           write(ifll,fmt="(' Total Parcels2       :',I6)") NPALLTOTAL2
            write(ifll,fmt="(' Min Relaxation Time :',E15.6)")
     +          dMinRelaxTime
            if(dMinRelaxTime .lt. dTSave/DBLE(ITER_INNER))
     +        write(ifll,*) "[Warning]: Relaxation Time exceeds dt!"
            write(ifll,fmt="(' Max Passed Cell     :',I3)") 
     &        iMaxPassedCell
            DO INI=1,N_injctr
            if(.not.lBreakupNONE(INI)) then 
              write(ifll,fmt="(' Num of Parcels broken up:',I4)")
     +            iBreakupCount
            endif
            enddo
!
            if(iTotalCountOfRW .gt. 0) then
              write(ifll,fmt="(' MaxContributionOfRandomWalk :',
     &        E15.6,' %')")
     +        dMaxRandomWalk
              write(ifll,fmt="(' AveContributionOfRandomWalk :',
     &        E15.6,' %')")
     +        dAveRandomWalk/DBLE(iTotalCountOfRW+SML)
            endif
            write(ifll,fmt="(77('+'))")
          endif
        endif
      endif
!
      DT0=dTSave

      if(reaction.and.ical_P_R==1) then
        DEALLOCATE(WDOTP,tempp)
      endif
!---------------------
! --- temp memory 
!---------------------
      if(IFLAG_COAL==icoal) then
        PMASS(:,2)=PMASSC(:)
      endif
!
      RETURN !zreturn

      END SUBROUTINE FFR_TRACING
!
!*********************************************************************
!
! Making Cell Nieghbor lists and connections.
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FFR_TRACING01(NMAT,MXALLCV,MXCVFAC,MEP,
     1           MAT_NO,MAT_CFIDX,MAT_CAL,MAT_CVEXT,
     2           LVEDGE,SFAREA,SFCENT,CVCENT,
     3           VCTRP,NEIGHBOR,DNORM,IERROR,BCKND,
     4           MXSSFBC,MXCVFAC_P,LBC_SSF)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only       : ifle,ifll
      USE module_dimension,only : MXCV_VP,MXBND_VP
      USE MODULE_BOUNDARY,ONLY  : NBCND,MAT_BCIDX,LBC_INDEX
      use module_trace ,ONLY    : CC
      use module_particle ,only : N_injctr
      
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)  :: NMAT,MXALLCV,MXCVFAC,MEP,MXCVFAC_P
      INTEGER,INTENT(IN)  :: MAT_NO(      0:NMAT)
      INTEGER,INTENT(IN)  :: MAT_CFIDX(   0:NMAT)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:NMAT)
      LOGICAL,INTENT(IN)  :: MAT_CAL(     0:NMAT)
      INTEGER,INTENT(IN)  :: LVEDGE(   2,MXCVFAC)
      REAL*8,INTENT(IN)   :: SFAREA(   4,MXCVFAC)
      REAL*8,INTENT(IN)   :: SFCENT(   3,MXCVFAC)
      REAL*8,INTENT(IN)   :: CVCENT(   3,MXALLCV)
!
      INTEGER,INTENT(INOUT)  :: VCTRP(MXCV_VP,0:MXBND_VP)
      INTEGER, INTENT(IN) :: NEIGHBOR(MEP,MXALLCV)
      REAL*8,INTENT(INOUT)   :: DNORM(MXCVFAC) 
      INTEGER, INTENT(INOUT) :: IERROR 
!
      INTEGER,INTENT(IN)  :: MXSSFBC
      INTEGER,INTENT(IN)  :: LBC_SSF(MXSSFBC)
      INTEGER,INTENT(OUT) :: BCKND(MXCVFAC_P)
! BCKND
! --- [LOCAL ENTITIES]
! 
      INTEGER :: IIMAT,IMAT,ICFL,ICFLL,ICFS,ICFE,ICVA,ICVB,IA,IB,IFA
      INTEGER :: ICMAX,ICMIN,I,J,K,ICV,IDC,IBFL,IBFS,IBFE,ICVL,IC
      REAL*8  :: TEMP,TEMP1
      INTEGER :: IAMAX
      INTEGER :: K2,nb
      LOGICAL :: LDUMMYA, LDUMMYB
!       
!
      ICMAX=MAT_CVEXT(NMAT)
      ICMIN=MAT_CVEXT(0)+1
!
      DO 20 nb=1,NBCND
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        DO 40 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          BCKND(ICFL)=nb
  40    CONTINUE
  20  CONTINUE
!
      DO ICVL=ICMIN,ICMAX
      DO IC=1,VCTRP(ICVL,0)
      ICFLL=VCTRP(ICVL,IC)
      ICFL=abs(ICFLL)
      CC(IC)=-DNORM(ICFL)
      DO I=1,3
      CC(IC)=CC(IC)+SFAREA(I,ICFL)
     *             *CVCENT(I,ICVL)
      ENDDO
      CC(IC)=CC(IC)*SIGN(1,ICFLL)
!
      if(CC(IC).LT.0) then 
        VCTRP(ICVL,IC)=-VCTRP(ICVL,IC) 
        ICFLL=VCTRP(ICVL,IC)
        ICFL=abs(ICFLL)
        CC(IC)=-DNORM(ICFL)
        DO I=1,3
        CC(IC)=CC(IC)+SFAREA(I,ICFL)
     *               *CVCENT(I,ICVL)
        ENDDO
        CC(IC)=CC(IC)*SIGN(1,ICFLL)
      endif
!
      IF(CC(IC).LT.0) THEN
        WRITE(ifll,*) "ERROR: ICVL=,IC=,CC(IC)=",ICVL,IC,CC(IC),
     + (CVCENT(I,ICVL),I=1,3),DNORM(ICFL),(SFCENT(I,ICFL),I=1,3)
        WRITE(ifll,*) "MSG: Redefining initial particle coordinates"
!        call FFRABORT(1,'ERR: FFR_TRACING01')
      ENDIF
      ENDDO
      ENDDO
!
      RETURN
      END  SUBROUTINE FFR_TRACING01   
!
!********************************************************************
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FFR_TRACING03(INI,
     *  IP,HEIGHT,IPCV,XYZP0,INSIDE,VCTRP,SFAREA,
     *  DT0,CVCENT,DNORM)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only       : ifle,ifll
      USE module_dimension 
      use module_trace ,only   : CC,HEIGHT2=>HEIGHT0,
     &                           HEIGHT1,VELN
      use module_particle,only : N_injctr
!
      IMPLICIT NONE
!
!------------------------
! --- [DUMMY ARGUMENTS]
!------------------------
      INTEGER,intent(in)     :: IP,IPCV,INI
      REAL*8 ,intent(in)     :: DT0 
      INTEGER,intent(in)     :: VCTRP(MXCV_VP ,0:MXBND_VP) 
      REAL*8 ,intent(in)     :: SFAREA(4,    MXCVFAC)   
      REAL*8 ,intent(in)     :: CVCENT(3,    MXALLCV_P) 
      REAL*8 ,intent(in)     :: DNORM(       MXCVFAC_P) 
      REAL*8 ,intent(inout)  :: HEIGHT(MEP,  MXPRT)    
      INTEGER,intent(inout)  :: INSIDE(      MXPRT)
      REAL*8 ,intent(inOUT)  :: XYZP0(3)
!
!     <<<<<<< VARIABLES FOR COLLISION CALCULATION >>>>>>> 
!
!-----------------------
! --- [LOCAL ENTITIES]
!-----------------------
      INTEGER :: K1,K2,ICFL,ICFLL,INB,I
!
      REAL*8 :: XYZP1(3),UVWP0(3)
      REAL*8 :: VNORMAL,TIMEMIN
!
      XYZP1(1)=CVCENT(1,IPCV)
      XYZP1(2)=CVCENT(2,IPCV)
      XYZP1(3)=CVCENT(3,IPCV)
!
      UVWP0(1)=(XYZP0(1)-XYZP1(1))/DT0
      UVWP0(2)=(XYZP0(2)-XYZP1(2))/DT0
      UVWP0(3)=(XYZP0(3)-XYZP1(3))/DT0
      INSIDE(IP)=IPCV
!
      DO INB=1,VCTRP(IPCV,0)
      ICFLL=VCTRP(IPCV,INB)
      ICFL=abs(ICFLL)
      CC(INB)=-DNORM(ICFL) 
      DO I=1,3 
      CC(INB)=CC(INB)+SFAREA(I,ICFL)*XYZP1(I)
      ENDDO
      CC(INB)=CC(INB)*SIGN(1,ICFLL)
      HEIGHT2(INB)=CC(INB)
      ENDDO  	  
!
      TIMEMIN=5.d0*DT0
      DO INB=1,VCTRP(IPCV,0)
      VNORMAL=0.d0
      ICFLL=VCTRP(IPCV,INB)
      ICFL=abs(ICFLL)
      DO I=1,3
      VNORMAL=VNORMAL-UVWP0(I)*SFAREA(I,ICFL)*SIGN(1,ICFLL)
      ENDDO
      IF(HEIGHT2(INB).GT.0.AND.HEIGHT2(INB).LT.VNORMAL*TIMEMIN) THEN
        TIMEMIN=HEIGHT2(INB)/VNORMAL
      ENDIF
      VELN(INB)=VNORMAL 
      ENDDO 
!
      IF(TIMEMIN.GT.DT0) THEN
        DO INB=1,VCTRP(IPCV,0) 
        HEIGHT2(INB)=HEIGHT2(INB)-VELN(INB)*DT0 
        ENDDO 
      ELSE
        DO INB=1,VCTRP(IPCV,0)
        HEIGHT2(INB)= HEIGHT2(INB)-0.95*VELN(INB)*TIMEMIN
        XYZP0(1)=XYZP1(1)+0.95*UVWP0(1)*TIMEMIN
        XYZP0(2)=XYZP1(2)+0.95*UVWP0(2)*TIMEMIN
        XYZP0(3)=XYZP1(3)+0.95*UVWP0(3)*TIMEMIN
        ENDDO 
      ENDIF

      DO INB=1,VCTRP(IPCV,0)
      HEIGHT(INB,IP)=HEIGHT2(INB)
      ENDDO     
!        
      RETURN
!	
      END SUBROUTINE FFR_TRACING03
!
!********************************************************************!
!
! Tracing particle by first order explicit.
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE  FFR_TRACING02_F(INI,
     *  IP,XYZP,UVWP,HEIGHT,FU,FV,FW,RHOP,
     *  INSIDE,IPCV,KE0,UVWP0,VCTRP,HEIGHT1,SFAREA,SFCENT,
     *  DT0,GRDC,CVCENT,DIA,
     *  VEL,RHO,RMU,NEIGHBOR,DNORM,BCKND,
     *  JPASS,JUDGE,JOUT,MOMCELL,PMASS,PARCEL,
     *  IPCELL,TRESIDUAL,CVVOLM,aks,LVEDGE,WDRBND)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!      use module_metrix,only : sta1,sta2,req1,req2

      use module_io,only        : ifle,ifll
      use module_dimension
      use module_hpcutil
      USE module_particle ,ONLY : dMinRelaxTime,
     &                            iMaxPassedCell,N_injctr,
     &                            iglass,bc_nb
      use module_trace ,only    : CC,VELN
      use module_gravity ,only  : ggg
      use module_constant
      use module_boundary,only : partkd,kpart_s,kpart_r,kpart_d,pkdwall
      use module_model,only  : ical_WDR
!
!----------------------------
! for parallel code debuging 
!----------------------------
! --- [DUMMY ARGUMENTS] 
!----------------------------
!
      IMPLICIT NONE
!
!
      INTEGER,intent(in)   :: IP,KE0,INI
      INTEGER,intent(inout):: IPCV
      REAL*8, intent(inout):: XYZP    (MXPRT,3)
      REAL*8, intent(inout):: UVWP    (MXPRT,3)
      REAL*8, intent(inout):: HEIGHT  (MEP,MXPRT)
      REAL*8, intent(inout):: FU      (MXPRT,2)
      REAL*8, intent(inout):: FV      (MXPRT,2)
      REAL*8, intent(inout):: FW      (MXPRT,2)
      real*8, intent(in)   :: RHOP    (MXPRT)
      INTEGER,intent(inout):: INSIDE  (MXPRT)
      REAL*8, intent(inout):: UVWP0   (3)
      INTEGER,intent(in)   :: VCTRP(MXCV_VP,0:MXBND_VP)
      REAL*8, intent(in)   :: HEIGHT1 (MEP)
      REAL*8, intent(in)   :: SFAREA  (4,MXCVFAC)
      REAL*8, intent(in)   :: DT0
      REAL*8, intent(in)   :: GRDC    (MXALLCV_P,3,3)
      REAL*8, intent(in)   :: CVCENT  (3,MXALLCV_P)
      REAL*8, intent(in)   :: DIA     (MXPRT) 
      REAL*8, intent(in)   :: VEL     (MXALLCV,3,2)
      real*8, intent(in)   :: RHO     (MXALLCV,2)
      REAL*8, intent(in)   :: RMU     (MXALLCV)
      REAL*8, INTENT(IN)   :: CVVOLM  (MXALLCV)
      real*8 ,intent(in)   :: aks     (MXALLCVR,MXRANS)
!
      INTEGER,intent(in)   :: NEIGHBOR(MEP,MXALLCV_P)
      REAL*8, intent(in)   :: DNORM   (    MXCVFAC_P)
!
      INTEGER,intent(inout):: JOUT    (    MXPRT)
      INTEGER,intent(inout):: JUDGE   (    MXCVFAC_P)
      REAL*8, intent(inout):: MOMCELL (    MXALLCV_P,4)
      REAL*8, intent(in)   :: PMASS   (    MXPRT)
      REAL*8, intent(in)   :: PARCEL  (    MXPRT)
      INTEGER,intent(in)   :: IPCELL  (    MXALLCV_P)
      REAL*8, intent(inout):: TRESIDUAL(   MXPRT)
      REAL*8 ,INTENT(IN)   :: SFCENT   (    3,MXCVFAC)
      INTEGER,INTENT(in)   :: BCKND(MXCVFAC_P)
!
      INTEGER,INTENT(INOUT) :: JPASS(MEVENT)
      INTEGER,INTENT(IN)    :: LVEDGE   (    2,MXCVFAC)
      real*8 ,intent(inout) :: WDRBND(        MXCV_WDR,N_inj)
!-------------------------
! --- [LOCAL ENTITIES]
!-------------------------
      INTEGER :: I,J,K,L,M,K1,K2,K3,NALL,INB,ICVB,kdp,nb
      INTEGER :: ICFL,ICFLL,MOUT,IRK,NTIME,KEN,KENEIB
      REAL*8 :: dum1,dum2,dl,U1,U2,U3,R1,R2,R3,ierr=0
!
!
!            <<<<<<< WORK AREA >>>>>>>
!
      REAL*8,PARAMETER :: third=1.d0/3.d0
!
      REAL*8 :: 
     1          UVW0(3),XYZP0(3),
     2          UVW1(3),XYZP1(3),DTIME
      REAL*8 :: VNORMAL,TIME,TIMEMIN,TIMEMIN1
      REAL*8 :: VN1,TEMP
      REAL*8 :: FXCELL,FYCELL,FZCELL
!
      real*8 :: PARTTIME
      real*8 :: dRelaxationTime,UVWFluc(3)
      integer:: KENB,ICVL
!
      TRESIDUAL(IP)=DT0
      TIMEMIN1=0.d0
      NALL=1
!
      IF(IPCELL(IPCV)>=1) THEN !zhang-
        JOUT(IP)=IPCELL(IPCV) 
        INSIDE(IP)=IPCV
        GOTO 50
      ENDIF
!      if(IPCV>NCV) then
!        INSIDE(IP)=NCV   !IPCV
!        JOUT(IP)=A_dead  !0   1111
!        GOTO 52
!      endif
      JPASS(NALL)=IPCV

!
 100  CONTINUE
!
      TIMEMIN=5.d0*TRESIDUAL(IP)
!
      DO 200 INB=1,VCTRP(IPCV,0)
      VNORMAL=0.d0
      ICFLL=VCTRP(IPCV,INB)
      ICFL=abs(ICFLL)
      DO I=1,3
      VNORMAL=VNORMAL-UVWP0(I)*SFAREA(I,ICFL)*SIGN(1,ICFLL)
      ENDDO
!
      IF(HEIGHT(INB,IP).GT.0
     &  .AND.HEIGHT(INB,IP).LT.VNORMAL*TIMEMIN) 
     & THEN
        TIMEMIN=HEIGHT(INB,IP)/(VNORMAL+SML)
        MOUT=INB
      ENDIF
      VELN(INB)=VNORMAL
 200  ENDDO 
!
!
      IF(TIMEMIN.GT.TRESIDUAL(IP)) THEN
!
        CALL FFR_TRACING_NEWUVW(INI,
     *  IP,XYZP,UVWP,HEIGHT,FU,FV,FW,RHOP,
     *  INSIDE,IPCV,UVWP0,VCTRP,SFAREA,
     *  TRESIDUAL(IP),GRDC,CVCENT,DIA,VEL,RHO,RMU,DNORM,
     *  MOMCELL,PMASS,PARCEL,JOUT,IPCELL,
     *  CVVOLM,aks)
!
        TRESIDUAL(IP)=0.d0
!--------------------------------------------------------
!    IF THIS CELL IS A OUTER CELL (other CPU region)
!--------------------------------------------------------
!
        DO INB=1,VCTRP(IPCV,0)
        ICFLL=VCTRP(IPCV,INB)
        ICFL=abs(ICFLL)
        CC(INB)=-DNORM(ICFL)
        DO I=1,3
        CC(INB)=CC(INB)+SFAREA(I,ICFL)*XYZP(IP,I)
        ENDDO
        CC(INB)=CC(INB)*SIGN(1,ICFLL)
        ENDDO
!
        DO INB=1,VCTRP(IPCV,0)
        HEIGHT(INB,IP)=CC(INB)
        ENDDO
!
        GO TO 50
!
      ELSE
!-----------------------------------------------
! --- TRACE THE PARTICLE IN THE NEIGHBOR CELL   
!-----------------------------------------------
        DTIME=TIMEMIN
        TRESIDUAL(IP)=TRESIDUAL(IP)-DTIME
!
        IF(NEIGHBOR(MOUT,IPCV)<0) THEN
          ICFL=ABS(VCTRP(IPCV,MOUT))
          nb=BCKND(ICFL)
!          if(nb==bc_nb(INI)) goto 50
          if(nb==0) CALL FFRABORT(1,'ERR:BCKND-1')
          kdp=partkd(nb)
          if(kdp==kpart_d) then
            if(ical_WDR==1.and.pkdwall(nb).and.JOUT(IP)==0
     &        .and.nb/=bc_nb(INI)) then
              ICVL=LVEDGE(1,ICFL)
              dum1=PI*DIA(IP)**3/6.0D0*PARCEL(IP)/SFAREA(4,ICFL)
              WDRBND(ICVL,INI)=WDRBND(ICVL,INI)+dum1
              JOUT(IP)=BC_dead
            else
              JOUT(IP)=BC_dead
            endif
            ICVB=abs(NEIGHBOR(MOUT,IPCV))
            INSIDE(IP)=ICVB
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)+U1*dum1
            XYZP(IP,2)=XYZP(IP,2)+U2*dum1
            XYZP(IP,3)=XYZP(IP,3)+U3*dum1
            UVWP(IP,:)=0.d0
          elseif(kdp==kpart_s) then
            if(ical_WDR==1.and.pkdwall(nb).and.JOUT(IP)==0
     &        .and.nb/=bc_nb(INI)) then
              ICVL=LVEDGE(1,ICFL)
              dum1=PI*DIA(IP)**3/6.0D0*PARCEL(IP)/SFAREA(4,ICFL)
              WDRBND(ICVL,INI)=WDRBND(ICVL,INI)+dum1
              JOUT(IP)=BC_dead
            else
              JOUT(IP)=-nb   !-ICFL    !-1 !HIBI
            endif

!            if(.not.pkdwall(nb))  then
!              JOUT(IP)=BC_dead
!            endif
            ICVB=abs(NEIGHBOR(MOUT,IPCV))
            INSIDE(IP)=ICVB
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)+U1*dum1
            XYZP(IP,2)=XYZP(IP,2)+U2*dum1
            XYZP(IP,3)=XYZP(IP,3)+U3*dum1
            UVWP(IP,:)=0.d0
          elseif(kdp==kpart_r) then
             
!            JOUT(IP)=-ICFL*INT(INT(IPCV/NCV))      !0
            if(.not.pkdwall(nb))  then
              JOUT(IP)=BC_dead
            endif
            JOUT(IP)=0

            INSIDE(IP)=IPCV
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)

!            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)-U1*dum1
            XYZP(IP,2)=XYZP(IP,2)-U2*dum1
            XYZP(IP,3)=XYZP(IP,3)-U3*dum1
            if(dum2<0.d0) then
              goto 50
            endif

            R1=dl*dum2*SFAREA(1,ICFL)
            R2=dl*dum2*SFAREA(2,ICFL)
            R3=dl*dum2*SFAREA(3,ICFL)
            U1=UVWP(IP,1)-R1
            U2=UVWP(IP,2)-R2
            U3=UVWP(IP,3)-R3
            UVWP(IP,1)=U1-R1
            UVWP(IP,2)=U2-R2
            UVWP(IP,3)=U3-R3
          endif
          GOTO 50
        ENDIF

!
        CALL FFR_TRACING_NEWUVW(INI,
     *      IP,XYZP,UVWP,HEIGHT,FU,FV,FW,RHOP,
     *      INSIDE,IPCV,UVWP0,VCTRP,SFAREA,
     *      DTIME,GRDC,CVCENT,DIA,VEL,RHO,RMU,DNORM,
     *      MOMCELL,PMASS,PARCEL,JOUT,IPCELL,
     *      CVVOLM,aks)
!
        FU(IP,1)=FU(IP,2)
        FV(IP,1)=FV(IP,2)
        FW(IP,1)=FW(IP,2)
!
!
        NALL=NALL+1
        NALL=MIN(NALL,MEVENT)  !????
        IPCV=NEIGHBOR(MOUT,IPCV)
        JPASS(NALL)=IPCV
        INSIDE(IP)=IPCV
        if(IPCV==0) call FFRABORT(1,'IPCV=000000000')
!
        DO INB=1,VCTRP(IPCV,0)
        ICFLL=VCTRP(IPCV,INB)
        ICFL=abs(ICFLL)
        CC(INB)=-DNORM(ICFL)
        DO I=1,3
        CC(INB)=CC(INB)+SFAREA(I,ICFL)*XYZP(IP,I)
        ENDDO
        CC(INB)=CC(INB)*SIGN(1,ICFLL)
        ENDDO
!
        DO INB=1,VCTRP(IPCV,0)
        HEIGHT(INB,IP)=CC(INB)
        IF(HEIGHT(INB,IP).lt.0.d0) then
          HEIGHT(INB,IP)=-1.d0*CC(INB)
        endif
        ENDDO
!----------------------------------------------
!    IF THIS CELL IS A OUTER CELL (HPC-MPI) 
!----------------------------------------------
        IF(IPCELL(IPCV).GE.1) THEN 
          JOUT(IP)=IPCELL(IPCV)
          INSIDE(IP)=IPCV
          GOTO 50
        ENDIF 
!
        TIMEMIN1=TIMEMIN
!
        GOTO 100
        ENDIF
 50   CONTINUE
!
      IF(IPCELL(IPCV).GE.1) THEN 
        JOUT(IP)=IPCELL(IPCV)
        INSIDE(IP)=IPCV
      ENDIF 
 52   CONTINUE
!--------------------------------
! Count Number of Passed Cell.
!--------------------------------
      IF(NALL.gt.iMaxPassedCell) iMaxPassedCell=NALL
!
      DO K1=1,NALL
      KEN=JPASS(K1)
      if(KEN==0) cycle 
      JUDGE(KEN)=0
      ENDDO
!********************************************************************
! --- 
!********************************************************************
      RETURN
      END SUBROUTINE FFR_TRACING02_F
C
C********************************************************************
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE  FFR_TRACING_NEWUVW(INI,
     *  IP,XYZP,UVWP,HEIGHT,FU,FV,FW,RHOP,
     *  INSIDE,IPCV,UVWP0,VCTRP,SFAREA,
     *  DT0,GRDC,CVCENT,DIA,VEL,RHO,RMU,DNORM,
     *  MOMCELL,PMASS,PARCEL,JOUT,IPCELL,
     *  CVVOLM,aks)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only       : ifle,ifll
      use module_dimension
      USE module_particle
      use module_trace ,only   : CC,VELN
      use module_gravity,only  : ggg
      use module_constant
      use module_particle,only : ifuel,icoal,iliquid,iglass
      use module_particle,only : movkd,imovD,imovMilr,imovSomm,imovU
!------------------------
! --- [DUMMY ARGUMENTS] 
!------------------------
      IMPLICIT NONE
!
      INTEGER,intent(in)   :: IP,INI
      INTEGER,intent(inout):: IPCV
      REAL*8, intent(inout):: XYZP    (      MXPRT,3)
      REAL*8, intent(inout):: UVWP    (      MXPRT,3)
      REAL*8, intent(inout):: HEIGHT  (  MEP,MXPRT)
      REAL*8, intent(inout):: FU      (      MXPRT,2)
      REAL*8, intent(inout):: FV      (      MXPRT,2)
      REAL*8, intent(inout):: FW      (      MXPRT,2)
      real*8, intent(in)   :: RHOP    (      MXPRT)
      INTEGER,intent(inout):: INSIDE  (      MXPRT)
      INTEGER,intent(inout):: JOUT  (      MXPRT)
      REAL*8, intent(inout):: UVWP0   (      3)
      INTEGER,intent(in)   :: VCTRP    (MXCV_VP,0:MXBND_VP)
      REAL*8, intent(in)   :: SFAREA  (4,    MXCVFAC)
      REAL*8, intent(in)   :: DT0
      REAL*8, intent(in)   :: GRDC    (      MXALLCV,3,3)
      REAL*8, intent(in)   :: CVCENT  (3,    MXALLCV)
      REAL*8, intent(in)   :: DIA     (MXPRT)
      REAL*8, intent(in)   :: VEL     (      MXALLCV,3,2)
      real*8, intent(in)   :: RHO     (      MXALLCV,2)
      REAL*8, intent(in)   :: RMU     (      MXALLCV)
      REAL*8, intent(in)   :: DNORM   (      MXCVFAC_P)
      INTEGER,INTENT(INOUT):: IPCELL  (      MXALLCV_P)
!
      REAL*8, intent(inout):: MOMCELL (      MXALLCV_P,4)
      REAL*8, intent(in)   :: PMASS   (      MXPRT)
      REAL*8, intent(in)   :: PARCEL  (      MXPRT)
      REAL*8, INTENT(IN)   :: CVVOLM  (      MXALLCV)
      real*8 ,intent(in)   :: aks     (      MXALLCVR,MXRANS)
!------------------------
! --- [LOCAL ENTITIES]
!------------------------
      INTEGER :: K1,K2,K3
!
!
C            <<<<<<< WORK AREA >>>>>>>
C--------------------------
C Gravity Force Setting
C--------------------------
!
      REAL*8 :: UVW0F(3),XYZP0(3),dum1,dum2,PA_MS,FDD
      REAL*8 :: VNORMAL
      REAL*8 :: RE_P,SLIPV,FCD,FTAU,US,VS,WS,F1,FTAU_r
      REAL*8 :: FXCELL,FYCELL,FZCELL
!
      real*8 :: dRelaxationTime,UVWFluc(3),y,body(3)
      integer:: KENB,I,INB
      REAL*8 :: dTempRandomWalk
      REAL*8 :: aaa,bbb,U_blow,Re_b,D3,VOLUME,DP,b1,b2,b3,b4,rmda,fai
!------------------------------------------------------------
! --- INTERPOLATE THE PHYSICAL VARIABLES IN THE PRESENT CELL
!------------------------------------------------------------
      DO I=1,3
      XYZP0(I)=XYZP(IP,I)+UVWP0(I)*DT0
      ENDDO
!
      DO INB=1,VCTRP(IPCV,0)
      HEIGHT(INB,IP)=HEIGHT(INB,IP)-VELN(INB)*DT0
      IF(HEIGHT(INB,IP).LT.0) THEN
        CALL FFR_TRACING03(INI,
     &         IP,HEIGHT,IPCV,XYZP0,INSIDE,VCTRP,SFAREA,
     *         DT0,CVCENT,DNORM)
        GO TO 120
      ENDIF
      ENDDO
 120  CONTINUE
!
      IF(IPCELL(IPCV).GE.1) THEN 
        JOUT(IP)=IPCELL(IPCV)
        INSIDE(IP)=IPCV
        return
      ENDIF 
!--------------------------------
! --- 
!--------------------------------
      UVW0F(1)=VEL(IPCV,1,2)
     *         +GRDC(IPCV,1,1)*(XYZP0(1)-CVCENT(1,IPCV))
     *         +GRDC(IPCV,2,1)*(XYZP0(2)-CVCENT(2,IPCV))
     *         +GRDC(IPCV,3,1)*(XYZP0(3)-CVCENT(3,IPCV))
      UVW0F(2)=VEL(IPCV,2,2)
     *         +GRDC(IPCV,1,2)*(XYZP0(1)-CVCENT(1,IPCV))
     *         +GRDC(IPCV,2,2)*(XYZP0(2)-CVCENT(2,IPCV))
     *         +GRDC(IPCV,3,2)*(XYZP0(3)-CVCENT(3,IPCV))
      UVW0F(3)=VEL(IPCV,3,2)
     *         +GRDC(IPCV,1,3)*(XYZP0(1)-CVCENT(1,IPCV))
     *         +GRDC(IPCV,2,3)*(XYZP0(2)-CVCENT(2,IPCV))
     *         +GRDC(IPCV,3,3)*(XYZP0(3)-CVCENT(3,IPCV))
!---------------------------------------------------------------------
! --- Add Gaussian Noise to Velocity components (Random Walk Model)
!---------------------------------------------------------------------
      CALL add_random_walk(IP,IPCV,GRDC,CVVOLM,UVWFluc,aks)
      UVW0F(1)=UVW0F(1)+UVWFluc(1)
      UVW0F(2)=UVW0F(2)+UVWFluc(2)
      UVW0F(3)=UVW0F(3)+UVWFluc(3)
C
      US=UVWP0(1)-UVW0F(1)
      VS=UVWP0(2)-UVW0F(2)
      WS=UVWP0(3)-UVW0F(3)
C
      DP=DIA(IP)
      SLIPV=SQRT(US*US+VS*VS+WS*WS)
      RE_P=DP*SLIPV*RHO(IPCV,1)/RMU(IPCV)
!
      U_blow=SQRT(UVW0F(1)**2+UVW0F(2)**2+UVW0F(3)**2)
      Re_b=DP*U_blow*RHO(IPCV,1)/RMU(IPCV)
!
      FTAU_r=18.D0*RMU(IPCV)/(RHOP(IP)*DP*DP)
      PA_MS=PMASS(IP)
!-------------------------------------
! --- FCD=24.d0/RE_P*F1
! --- FDD=FCD*RE_P/24.d0*FTAU_r
!-------------------------------------
      if(movkd(INI)==imovD) then
        F1=1.0D0+0.15d0*RE_P**0.687d0
        if(RE_P>=1.d3) F1=0.44d0*RE_P/24.d0
        FDD=F1*FTAU_r*ddcd(INI)
      elseif(movkd(INI)==imovMilr) then
        aaa=0.06d0+0.077d0*exp(-0.4d0*RE_P)
        bbb=0.4d0+0.77d0*exp(-0.04d0*RE_P)
        F1=(1.0D0+0.0545d0*RE_P
     &           +0.1d0*RE_P**0.5*max(0.d0,(1.d0-0.03d0*RE_P)))/
     &           (1.d0+aaa*(abs(Re_b))**bbb)
        if(RE_P>=1.d3) F1=0.44d0*RE_P/24.d0
        FDD=F1*FTAU_r*ddcd(INI)
      elseif(movkd(INI)==imovSomm) then
        F1=1.0D0+1.d0/6.d0**RE_P**(2.d0/3.d0)
        if(RE_P>=1.d3) F1=0.44d0*RE_P/24.d0
        FDD=F1*FTAU_r*ddcd(INI)
      elseif(movkd(INI)==imovU) then
        call USER_particle_mov_f1(INI,IP,DT0,
     &       SLIPV,U_blow,Re_P,RE_B,FTAU_r,
     &       PA_MS,UVW0F,UVWP0,RHO(IPCV,1),RHOP(IP),F1)
        FDD=F1*FTAU_r*ddcd(INI)
      elseif(movkd(INI)==imovMrsi) then !Morsi-Alexander
        fai=1.d0
        b1=exp(2.3288-6.4581*fai+2.4486*fai**2)
        b2=0.0964d0+0.5565*fai
        b3=exp(4.9050d0-13.8944d0*fai+18.4222d0*fai**2-10.2599d0*fai**3)
        b4=exp(1.4681d0+12.2584d0*fai-20.7322d0*fai**2+15.8855d0*fai**3)
        F1=1.d0+b1*RE_P**b2
        FCD=24.d0/RE_P*F1+b3*RE_P/(b4+RE_P)
        FDD=FCD*RE_P/24.d0*FTAU_r
      elseif(movkd(INI)==imovSub) then !Sub-micron particle
        rmda=1.d-6 ! molecular mean free path
        FDD=FTAU_r/(1.d0+2.d0*rmda/DP*
     &      (1.257d0+0.4d0*exp(1.1d0*DP/(2.d0*rmda))))
      endif

      dum1=1.d0+FDD*DT0
      body(:)=ggg(:)*(RHOP(IP)-RHO(IPCV,2))/RHOP(IP)
      UVWP0(1)=(UVWP(IP,1)+(body(1)+FDD*UVW0F(1))*DT0)/dum1
      UVWP0(2)=(UVWP(IP,2)+(body(2)+FDD*UVW0F(2))*DT0)/dum1
      UVWP0(3)=(UVWP(IP,3)+(body(3)+FDD*UVW0F(3))*DT0)/dum1

      FU(IP,2)=FDD*(UVWP0(1)-UVW0F(1))
      FV(IP,2)=FDD*(UVWP0(2)-UVW0F(2))
      FW(IP,2)=FDD*(UVWP0(3)-UVW0F(3))
!BABA-KUROSE
!      FU(IP,2)=FDD*UVWP0(1)
!      FV(IP,2)=FDD*UVWP0(2)
!      FW(IP,2)=FDD*UVWP0(3)
!
      PA_MS=PMASS(IP)*PARCEL(IP)
      FXCELL=(FU(IP,2))*PA_MS
      FYCELL=(FV(IP,2))*PA_MS
      FZCELL=(FW(IP,2))*PA_MS
      MOMCELL(IPCV,1)=MOMCELL(IPCV,1)+FXCELL*DT0   !N
      MOMCELL(IPCV,2)=MOMCELL(IPCV,2)+FYCELL*DT0
      MOMCELL(IPCV,3)=MOMCELL(IPCV,3)+FZCELL*DT0
      MOMCELL(IPCV,4)=MOMCELL(IPCV,4)+FDD*DT0*PA_MS
!
      DO I=1,3
      UVWP(IP,I)=UVWP0(I)
      XYZP(IP,I)=XYZP0(I)
      ENDDO
!
      END SUBROUTINE FFR_TRACING_NEWUVW
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C subroutine for making
C NUKIYAMA-TANAZAWA distribution from random number dTemp.
C A random number dTemp must be 0.0 to 1.0
C Changing dDY and iIter is not recommended.
C
C 2005.05.13 Jun ARAI
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine NUKIYAMA_TANAZAWA(INI,dDIAM, dSMD, dTemp)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only       : ifle,ifll
      USE module_particle,only : dAlpha,dBeta
!
      implicit none
!
! --- [module arguments]
!
      INTEGER,INTENT(IN)  :: INI
      real*8, INTENT(OUT) :: dDIAM
      real*8, INTENT(IN)  :: dSMD
      real*8, INTENT(IN)  :: dTemp
!
      real*8, parameter :: dDY = 0.001
      integer, parameter:: iIter = 2600
      real*8, save :: ALPHA, BETA
      real*8, save :: A,B
!
! --- LOCAL ENTITIES
!
!       dY: Diameter non-dimensionalized by Sauter Mean Diameter
!       dX: Integral of NUKIYAMA-TANAZAWA function
!
      real*8 :: dX,dY,dTemp2
      integer:: i
      logical,save :: lInit=.TRUE.
!
      if(lInit) then
        ALPHA = dAlpha(INI)
        BETA  = dBeta(INI)
        CALL CALC_NUKIYAMATANAZAWA_AB(INI,ALPHA,BETA,A,B)
        lInit = .FALSE.
      endif
!
      dX=0.d0
      dY=0.d0
!
      do i=1,iIter
      dY=dY+dDY
      dX=dX+A*dY**ALPHA*DEXP(-1.0*B*dY**BETA)*dDY
      if(dX.ge.dTemp) goto 100
      enddo
 100  continue
!
      call utl_random(dTemp2)
      dDIAM=(dY-dTemp2*dDY)*dSMD
        
      IF(dDIAM .le. 0.d0) GOTO 100
        
      return
      end subroutine NUKIYAMA_TANAZAWA
!
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C subroutine for making
C NUKIYAMA-TANAZAWA distribution from random number dTemp.
C Use this routine, if only the volume of each parcel is constant.
C A random number dTemp must be 0.0 to 1.0
C Changing dDY and iIter is not recommended.
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine NUKIYAMA_TANAZAWA_VOLUMEBASE(INI,dDIAM,dSMD,dTemp)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only       : ifle,ifll
      USE module_particle,only : dAlpha,dBeta
!      
      implicit none
!
! --- [module arguments]
!
      INTEGER,INTENT(IN) :: INI
      real*8, INTENT(OUT):: dDIAM
      real*8, INTENT(IN) :: dSMD
      real*8, INTENT(IN) :: dTemp
C
      real*8, parameter :: dDY = 0.002
      integer, parameter:: iIter = 2000
C----------------------------------------------------------------------
C       LOCAL ENTITIES
C
C       dY: Diameter non-dimensionalized by Sauter Mean Diameter
C       dX: Integral of NUKIYAMA-TANAZAWA function (Volume Based)
C       dNormalize : A constant for normalizing the function.
C----------------------------------------------------------------------
      real*8, save :: ALPHA, BETA
      real*8, save :: A,B

      real*8	:: dX,dY,dTemp2
      integer	:: i
      logical,save :: lInit=.TRUE.
      real*8, save :: dNormalize
!-------------------------------------------------
! --- 
!-------------------------------------------------
      if(lInit) then
        ALPHA=dAlpha(INI)
        BETA=dBeta(INI)
        CALL CALC_NUKIYAMATANAZAWA_AB(INI,ALPHA,BETA,A,B)
        dX=0.d0
        dY=0.d0
        do i=1,iIter
        dY=dY+dDY
        dX=dX+A*dY**ALPHA*DEXP(-1.0*B*dY**BETA)*dY**3*dDY
        enddo
        dNormalize=dX
        lInit = .FALSE.
      endif
!
      dX=0.d0
      dY=0.d0
!
      do i=1,iIter
      dY=dY+dDY
      dX=dX+A*dY**ALPHA*DEXP(-1.0*B*dY**BETA)*dY**3*dDY/dNormalize
      if(dX.ge.dTemp) goto 100
      enddo
 100  continue
!
      call utl_random(dTemp2)
      dDIAM=(dY-dTemp2*dDY)*dSMD
!        
      IF(dDIAM.le.0.d0) GOTO 100
!        
      return
!
      end subroutine NUKIYAMA_TANAZAWA_VOLUMEBASE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
C Calculate A and B in the NUKIYAMA-TANAZAWA Distribution
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine CALC_NUKIYAMATANAZAWA_AB(INI,ALPHA,BETA,A,B)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_io,only       : ifle,ifll
      USE module_hpcutil,	only :my_rank,root
!
      implicit none
!
      INTEGER,INTENT(IN) :: INI
      real*8, intent(in) :: ALPHA,BETA
      real*8, intent(out):: A, B
!
      integer :: i
      real*8  :: temp1, temp2, temp3
!
      call calc_gamma((ALPHA+4.d0)/BETA,temp1)
      call calc_gamma((ALPHA+3.d0)/BETA,temp2)
      call calc_gamma((ALPHA+1.d0)/BETA,temp3)
!
      A= (temp1/temp2)**(ALPHA+1.d0)*BETA/temp3
      B= (temp1/temp2)**BETA
!
      if(my_rank.eq.root) then
        write(ifll,*) "+++ NUKIYAMA-TANAZAWA DISTRIBUTION +++"
        write(ifll,*) "ALPHA=",alpha,",BETA =",beta,", B=",B,", A=",A
        write(ifll,fmt="(50('+'))")
      endif
!
      end subroutine CALC_NUKIYAMATANAZAWA_AB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C Subroutine for calculating rough approximation of GAMMA Fucntion
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine calc_gamma(z,gm)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only       : ifle,ifll
      implicit none
!
      real*8, intent(in) :: z
      real*8, intent(out):: gm
!
      integer ::iz,i
      real*8  ::zm,gzm,temp
      real*8, parameter :: eulerconst=0.57721566490153286060651209d0
!
      gm=1.d0
      zm=z
!
      if(z.gt.0) then
        iz=INT(z)
        do i=1,iz-1
        zm=zm-1.d0
        gm=gm*zm
        enddo
        if(zm.eq.1.d0) return
      else
        iz=-1*INT(z)
        do i=1,iz
        gm=gm/zm
        zm=zm+1.d0
        enddo
      endif
!
      gzm=DEXP(-1.d0*eulerconst*zm)/zm
      do i=1,50000
      temp=DEXP(zm/DBLE(i))/(1.d0+zm/DBLE(i))
      gzm=gzm*temp
      enddo
!
      gm=gm*gzm
!
      return
!
      end subroutine calc_gamma
!
!
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine add_random_walk(N,ICVL,GRDC,CVVOLM,UVWFluc,aks)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      use module_dimension
      use module_io,only       : ifle,ifll
      use module_model, only : idrdp,incomp,mach0,icaltb,KE2S,
     &                            ke,sles,dles,lles,noturb,RSM,ke_low
      use module_scalar,only : ike
      use module_time,  only : steady
!
      implicit none
!------------------------
! --- Dummy Auguments 
!------------------------
      integer, intent(IN) :: ICVL,N
      REAL*8,  INTENT(IN) :: CVVOLM(   MXALLCV)
      REAL*8,  INTENT(IN) :: GRDC  (   MXALLCV,3,3)
      real*8  ,intent(in) :: aks   (   MXALLCVR,MXRANS)
      real*8 , intent(OUT):: UVWFluc(3)
!------------------------
! --- Local Entities
!------------------------
      real*8 ,parameter :: r1pn=1.d0/3.d0, Cnut=5.d-2, Ceps=1.d0
      real*8 ,parameter :: twothirds=2.d0/3.d0
      real*8 ,parameter :: twopi=2.d0*3.1415926535897932d0
      real*8  :: SIJ11,SIJ22,SIJ33,SIJ12,SIJ13,SIJ23,del
      real*8  :: dspd,dKs,dTemp(2)
      integer :: i
      real*8,save :: MAXKS=0.d0
!-----------------------------------------------------------
! --- Call random seed on calculating the first particle.
!-----------------------------------------------------------
      if(N.eq.1) CALL RANDOM_SEED()
      dKs=0.d0
!
! --- < 1. Calculate dissipation function >- 2*SijSij
!
!      if(icaltb.eq.sles.or.(.not.steady)) then
      if(icaltb==sles) then
! --- LES Smagorinsky model
        SIJ11=grdc(ICVL,1,1)
        SIJ22=GRDC(ICVL,2,2)
        SIJ33=GRDC(ICVL,3,3)
        SIJ12=GRDC(ICVL,1,2)+GRDC(ICVL,2,1)
        SIJ13=GRDC(ICVL,1,3)+GRDC(ICVL,3,1)
        SIJ23=GRDC(ICVL,2,3)+GRDC(ICVL,3,2)
        dspd=(2.d0*SIJ11*SIJ11
     &       +2.d0*SIJ22*SIJ22
     &       +2.d0*SIJ33*SIJ33
     &       +(SIJ12*SIJ12+SIJ13*SIJ13+SIJ23*SIJ23))
        del=CVVOLM(ICVL)**r1pn
        dKs=(Cnut/Ceps)*del**2*DSQRT(dspd)
      elseif(icaltb==ke.or.icaltb==KE2S) then
! --- K-E model
        dKs=aks(ICVL,ike(1))
        continue
      elseif(icaltb.eq.ke_low) then
        dKs=aks(ICVL,ike(1))
        continue
      elseif(icaltb.eq.RSM) then
        continue
      endif
C----------------------------------------------------------------------------
C Calculate Gaussian Noise from 2 random numbers by inverse function method.
C (Normally distributed random number)
C----------------------------------------------------------------------------
      call utl_random(dTemp(1))
      call utl_random(dTemp(2))
      if(dTemp(1) .le. 0.d0) dTemp(1) = 1.d-10
      if(dTemp(2) .le. 0.d0) dTemp(2) = 1.d-10
      UVWFluc(1)=dsqrt(-2.d0*DLOG(dTemp(1)))*DCOS(twopi*dTemp(2))
      UVWFluc(2)=dsqrt(-2.d0*DLOG(dTemp(1)))*DSIN(twopi*dTemp(2))
      call utl_random(dTemp(1))
      call utl_random(dTemp(2))
      if(dTemp(1).le.0.d0) dTemp(1)=1.d-10
      if(dTemp(2).le.0.d0) dTemp(2)=1.d-10
      UVWFluc(3)=dsqrt(-2.d0*DLOG(dTemp(1)))*DCOS(twopi*dTemp(2))
C
      do i=1,3
      UVWFluc(i)=UVWFluc(i)*dsqrt(dKs*twothirds)
      enddo
C
      if(dKs .gt. MAXKS) then
        MAXKS = dKs
      endif
      return
      end subroutine add_random_walk
C********************************************************************
C
C       Pseudo Random Number Generation unit (MT-Method)
C
C********************************************************************
* A C-program for MT19937: Real number version
*   genrand() generates one pseudorandom real number (double)
* which is uniformly distributed on [0,1]-interval, for each
* call. sgenrand(seed) set initial values to the working area
* of 624 words. Before genrand(), sgenrand(seed) must be
* called once. (seed is any 32-bit integer except for 0).
* Integer generator is obtained by modifying two lines.
*   Coded by Takuji Nishimura, considering the suggestions by
* Topher Cooper and Marc Rieffel in July-Aug. 1997.
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Library General Public
* License as published by the Free Software Foundation; either
* version 2 of the License, or (at your option) any later
* version.
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
* See the GNU Library General Public License for more details.
* You should have received a copy of the GNU Library General
* Public License along with this library; if not, write to the
* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
* 02111-1307  USA
*
* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
* When you use this, send an email to: matumoto@math.keio.ac.jp
* with an appropriate reference to your work.
*
************************************************************************
* Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
*
*   genrand()      -> double precision function grnd()
*   sgenrand(seed) -> subroutine sgrnd(seed)
*                     integer seed
*
* This program uses the following non-standard intrinsics.
*   ishft(i,n): If n>0, shifts bits in i by n positions to left.
*               If n<0, shifts bits in i by n positions to right.
*   iand (i,j): Performs logical AND on corresponding bits of i and j.
*   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
*   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
*
************************************************************************
* subroutine mt_random(randomnumber)
* if this routine is called at first time, initialization will be done
* by using system time. MOD(system time, sec*min*hour)
*
************************************************************************
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine mt_random(randomnumber)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
      use module_io,only       : ifle,ifll
      USE module_hpcutil,	only :my_rank,NPE
      implicit none
*
* Period parameters
      integer, parameter :: N     =  624
      integer, parameter :: N1    =  N+1
      integer, parameter :: M     =  397
      integer, parameter :: MATA  = -1727483681
*                                    constant vector a
      integer, parameter :: UMASK = -214748648 !-2147483648
*                                    most significant w-r bits
      integer, parameter :: LMASK =  214743647
*                                    least significant r bits
* Tempering parameters
      integer, parameter :: TMASKB = -165838656  !-1658038656
      integer, parameter :: TMASKC = -272236544
*
      integer, save :: mt(0:N-1), mti,seed
      real*8  :: grnd,randomnumber
      integer :: iTimeSys(2),iTimes(9),y,kk
      integer :: TSHFTU,TSHFTS,TSHFTT,TSHFTL
      logical, save :: lInit=.TRUE.
*
      integer,save :: mag01(0:1)
      data mti/N1/
*                     mti==N+1 means mt[N] is not initialized
      data mag01/0, MATA/
*                        mag01(x) = x * MATA for x=0,1
C
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
C-------------------
C Initialization
C-------------------
      if(lInit) then
!        CALL DTIME(iTimeSys)
!        CALL GMTIME(iTimeSys(2),iTimes)
!        seed = mod(iTimeSys(2),iTimes(1)*iTimes(2)*iTimes(3)+1)
!        if(seed.eq.0) seed = 4357
        seed = 4357
        if(NPE>1) then
          write(ifll,*) "MSG:[mt_random]:Random Seed =",seed,
     +               ",my_rank =",my_rank
        else
          write(ifll,*) "MSG:[mt_random]:Random Seed =",seed
        endif
        
        mt(0)= iand(seed,-1)
        do 900 mti=1,N-1
          mt(mti) = iand(69069*mt(mti-1),-1)
  900   continue
        lInit = .FALSE.
      endif
C
C
*
      if(mti.ge.N) then  !generate N words at one time
        if(mti.eq.N+1) then  !if not initialized
          seed = 4357
          mt(0)= iand(seed,-1)
          do 950 mti=1,N-1
            mt(mti)=iand(69069 * mt(mti-1),-1)
  950     continue
*                              a default initial seed is used
        endif
*
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
           y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
           mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
*
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
*
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif
*
      randomnumber=grnd
      return
      end subroutine mt_random
C
C********************************************************************
C
C Subroutine For Fan Spray Nozzle (Volume based DDM)
C
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE READ_PARTICLE_INI(NPSUM,NCV,MXALLCV,
     &  MXPRT,PI,MXCOMP,MXCVFAC,MXSSFBC,
     &  SFCENT,LBC_SSF,SFAREA,LVEDGE,CVCENT,
     *  DIA,XYZP,UVWP,RHOP,PARCEL,DT0,TMPP,PYS,
     &  PMASS,INSIDE,JOUT,HIS_P) 
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_io,only      : ifle,ifll 
      USE module_particle
      USE module_hpcutil,only : my_rank,root,IBC_P_NO,NPE,IBC_P_TOL,
     &                          NEIBPETOT 
      use module_model,  only : ical_t,encomp,ical_reac,ical_WDR,Rh
      use module_particle,only: Tini_p,Cp_pa,heat_f,eps_rad,
     &                          Tref_p
      use module_particle,only: ifuel,icoal,iliquid,isolid
      USE module_particle,only: IFLAG_COAL,iglass          ! debag OSHIMA 2021.10.25  by NuFD
!     USE module_particle,only: IFLAG_COAL,icoal,iglass    << debag delete OSHIMA 2021.10.25 by NuFD 
      use module_constant,only: SML
      use module_trace ,only  : WK1_P
      use module_boundary,only : nobcnd,LBC_INDEX,nbcnd
      use module_particle,only : starttime,PMASSC
      use module_scalar,only   : calxi,idifflm,ORG 
      use module_particle,only : dropkd,iconst
      use module_scalar,only   : ical_flmlt
!
      IMPLICIT NONE
!
!--------------------
! Dummy Auguments 
!--------------------
!
      INTEGER, INTENT(IN)    :: MXPRT,MXCOMP,MXCVFAC,MXSSFBC,NCV,MXALLCV
      INTEGER, INTENT(INOUT) :: NPSUM(NPE)
      REAL*8,  INTENT(IN)    :: PI,DT0
      real*8 , intent(in)    :: SFCENT(3,MXCVFAC)
      INTEGER ,intent(in)    :: LBC_SSF( MXSSFBC)
      REAL*8,  INTENT(INOUT) :: DIA     (MXPRT)
      REAL*8,  INTENT(INOUT) :: XYZP    (MXPRT,3)
      REAL*8,  INTENT(INOUT) :: UVWP    (MXPRT,3)
      REAL*8,  INTENT(INOUT) :: RHOP    (MXPRT)
      REAL*8,  INTENT(INOUT) :: PARCEL  (MXPRT)
      REAL*8,  INTENT(INOUT) :: PYS     (MXPRT,MXCOMP,2)
      REAL*8,  INTENT(INOUT) :: TMPP    (MXPRT,2)
      REAL*8,  INTENT(INOUT) :: PMASS   (MXPRT,2)
      REAL*8 , INTENT(IN)    :: SFAREA(4,MXCVFAC)
      INTEGER, INTENT(INOUT) :: INSIDE(  MXPRT)
      INTEGER, INTENT(INOUT) :: JOUT  (  MXPRT)
      REAL*8 , INTENT(INOUT) :: HIS_P(   MXPRT)
      INTEGER, INTENT(IN)    :: LVEDGE(2,MXCVFAC)
      REAL*8 ,INTENT(IN)     :: CVCENT(3,MXALLCV)
!
!--------------------
! --- Local Entities 
!--------------------
!
      INTEGER :: I,ifl,INI,ierr1=0,NPALL1,IP,IPS,IPE,NO 
      REAL*8  :: DIAMIN,DIAMAX,D3,dFLOWRATE,dum1,dum2,dum3,dum4,dum5=0.d0
      REAL*8  :: VOLUME,TEMP1,TEMP2,TEMP3,ANGLE1,ANGLE2
      REAL*8  :: unitt(3),dum0,unit(3),dumm(3),SML1=1.d-15
      INTEGER :: icom,IBFS,IBFE,nb,IBFL,ICFL,IDUM,ICPU,ISUM,
     &           IPP,n0,ICV
      INTEGER,save :: INJCPU=0,ICONT=0,rain_BCNO
      REAL*8,save  :: rain_area=0.d0
      real*8,  save,allocatable :: nuP_m2(:),nuP_IBFL(:)
!
! --- 
!
      mass_p_total(:)=0.d0
      IBC_P_NO=0
      IBC_P_TOL=0
      NINALL(:)=0
      NPALL1=0
!
!------------------------------------------
! --- 
!------------------------------------------
!
      if(ical_WDR==1) then
        allocate(nuP_m2(N_injctr),stat=ierr1)
        if(ierr1/=0) call FFRABORT 
     &    (1,'ERR: nuP_m2 allocate error in READ_PARTICLE_INI')
!
        do INI=1,N_injctr
        NB=bc_nb(INI)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        dum1=0.d0
        do IBFL=IBFS,IBFE 
        ICFL=LBC_SSF(IBFL)
        dum1=dum1+SFAREA(4,ICFL) 
        enddo
        nuP_m2(INI)=dum1
        enddo
!
! --- 
!
        IDUM=0
        do nb=1,nbcnd 
        NO=nobcnd(nb) 
        do INI=1,N_injctr
        if(injkd(INI)/=inlet) cycle
        if(IDUM/=0) then
          if(NO_BC(INI)/=IDUM) then
            call FFRABORT(1,'ERR: rain BC MUST BE only one')
          endif
        endif
        if(NO_BC(INI)==NO) then 
          IDUM=NO_BC(INI)
        elseif(NO_BC(INI)==0)then
          call FFRABORT(1,'ERR: NOT defined BC no')
        endif
        enddo
        enddo
!
        NB=IDUM
        rain_BCNO=IDUM
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        dum1=0.d0
        do IBFL=IBFS,IBFE 
        ICFL=LBC_SSF(IBFL)
        dum1=dum1+SFAREA(4,ICFL) 
        enddo
        rain_area=dum1
!
        nu_intnst(:)=0.d0
        dum2=0.d0
        do INI=1,N_injctr
        dumm(1)=8000.d0*exp(-4.1d0*Rh**(-0.21)*1.D3*dSMD(INI)) !Nd
        dumm(2)=dumm(1)*PI*dSMD(INI)**3/6.d0  !m^3/m^3
!        
        dum1=140.d0*dSMD(INI)**0.5 !Vt
        rain_vel(INI)=dum1
        dumm(3)=dum1*dumm(2)   !m^3/m^2-s
        func_d(INI)=dumm(3)
        Rh_d(INI)=dumm(3)      !m^3/m^2-s
        dum2=dum2+dumm(3)      !m^3/m^2-s
        nu_intnst(INI)=dum1*dumm(1)  !1/m^2-s
        enddo
!
        func_d(:)=func_d(:)/(dum2+SML)
        dum3=Rh*1.d-3/3600.d0/dum2
!
        do INI=1,N_injctr
        nu_intnst(INI)=nu_intnst(INI)*dum3  !1/m^2-s

        Rh_d(INI)=Rh_d(INI)*dum3   !!m^3/m^2-s
        write(ifll,*)
        write(ifll,'(1x,a,I6,a,E15.4)') 
     &       'MSG: INJ_NO= ',INI,' Rain drop DMD(SMD)= ',dSMD(INI)

        write(ifll,'(1x,a,I16,a)') 
     &   'MSG: number density= ',INT(nu_intnst(INI)),' 1/m^2-s'
!
        dum1=(iSprayInjInter(INI)*DT0*rain_area*nu_intnst(INI))
        pcl_WDR(INI)=dum1/iSprayInjNum(INI)

        write(ifll,'(1x,a,E15.4,a,I4,a)') 
     &     'MSG: Drop number =',dum1,' by every ',
     &     iSprayInjInter(INI),' time steps'

        write(ifll,'(1x,a,I12,a,I4,a)') 
     &     'MSG: parcel number =',iSprayInjNum(INI),' by every ',
     &     iSprayInjInter(INI),' time steps'

        write(ifll,'(1x,a,E15.4,a)') 
     &     'MSG: Drop number =',pcl_WDR(INI),' in every parcel'

        enddo
!
      endif  !if(ical_WDR==1) then
!----------------------------------------
! ---
!----------------------------------------
      do 2000 INI=1,N_injctr
      if(iSprayInjTotal(INI)==0) then
        NINALL(INI)=iSprayInjNum(INI) 
     &  *max(1,1+int((iSprayInjEnd(INI)-iSprayInjStart(INI))
     &  /iSprayInjInter(INI)))
        if(iSprayInjTime(INI)==0) NINALL(INI)=0
        IF(iSprayInjTotal(INI)==0) NINALL(INI)=0
      else
        NINALL(INI)=min(
     &  iSprayInjNum(INI)
     &  *max(1,1+int((iSprayInjEnd(INI)-iSprayInjStart(INI))
     &  /iSprayInjInter(INI))),
     &  iSprayInjTotal(INI))
      endif
!
      SPHEIGHT(INI)=dSLITWIDTH(INI)*0.5d0
     &        /DTAN(dFANANGLE(INI)*0.5d0*PI/180.d0)
      NPALL1=NPALL1+NINALL(INI)
!
      if(NINALL(INI)==0) then
        VOLUMEPC(INI)=0.d0
      else
        VOLUMEPC(INI)=dInjVolume(INI)/DBLE(NINALL(INI))
      endif
 2000 enddo
!
      if(NPALL1.gt.MXPRT) then
        write(ifle,'(1X,a,2I10)') 'ERR: NPALL1>MXPRT',NPALL1,MXPRT
        goto 100
      endif

!---------------------------------------
! --- 
!---------------------------------------
      do INI=1,N_injctr
      NINALL(INI)=NINALL(INI-1)+NINALL(INI)
      enddo
!--------------------
! --- counter BC no
!--------------------
      do INI=1,N_injctr 
      if(.NOT.(injkd(INI)==i_coal.or.injkd(INI)==inlet)) cycle 
      do nb=1,nbcnd
      NO=nobcnd(nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(NO_BC(INI)==NO) then
        do IBFL=IBFS,IBFE
        IBC_P_NO(my_rank+1,INI,nb,1)=
     &           IBC_P_NO(my_rank+1,INI,nb,1)+1
        enddo
      endif
      enddo
      enddo
!-------------------------
! --- total BC counter
!-------------------------
      do INI=1,N_injctr 
      if(.NOT.(injkd(INI)==i_coal.or.injkd(INI)==inlet)) cycle 
      do nb=1,nbcnd 
      no=nobcnd(nb)
      if(NO_BC(INI)==no) then 
        IDUM=IBC_P_NO(my_rank+1,INI,nb,1) 
        call hpcisum(IDUM) 
        IBC_P_TOL(INI,nb,1)=IDUM 
! ---  
        do ICPU=1,NPE 
        IDUM=IBC_P_NO(ICPU,INI,nb,1) 
        call hpcibcast(ICPU-1,IDUM) 
        IBC_P_NO(ICPU,INI,nb,1)=IDUM 
        enddo 
      endif
      enddo
      enddo
!------------------------------------- 
! --- Particle by BC%
!-------------------------------------
      do INI=1,N_injctr
      if(.NOT.(injkd(INI)==i_coal.or.injkd(INI)==inlet)) cycle 
      do nb=1,nbcnd 
      NO=nobcnd(nb) 
      if(NO_BC(INI)==NO) then 
        dum1=dble(IBC_P_NO(my_rank+1,INI,nb,1))
     &      /dble(SML+IBC_P_TOL(INI,nb,1))
        idum=INT(dum1*(NINALL(INI)-NINALL(INI-1))) 
        IBC_P_NO(my_rank+1,INI,nb,1)=idum 
        IBC_P_TOL(INI,nb,1)=NINALL(INI)-NINALL(INI-1)
      endif
      enddo
      enddo
!-------------------------------------
! --- 
!-------------------------------------
      do 110 INI=1,N_injctr 
      if(.NOT.(injkd(INI)==i_coal.or.injkd(INI)==inlet)) cycle 
      do 113 nb=1,nbcnd 
      NO=nobcnd(nb) 
      if(NO_BC(INI)==NO) then 
        ISUM=0
        ICONT=0
        do 112 ICPU=1,NPE 
        IDUM=IBC_P_NO(ICPU,INI,nb,1) 
        call hpcibcast(ICPU-1,IDUM) 
        IBC_P_NO(ICPU,INI,nb,1)=IDUM 
        ISUM=ISUM+IDUM
        if(IDUM/=0) ICONT=ICPU
 112    enddo
        if(ICONT/=0) then 
          IDUM=IBC_P_TOL(INI,nb,1)-ISUM
          IBC_P_NO(ICONT,INI,nb,1)=
     &     IBC_P_NO(ICONT,INI,nb,1)+IDUM
        endif
      endif
 113  enddo
 110  enddo
!----------------------------------
      do INI=1,N_injctr 
        if(.NOT.(injkd(INI)==i_coal.or.injkd(INI)==inlet)) cycle 
        do nb=1,nbcnd 
        NO=nobcnd(nb) 
        if(NO_BC(INI)==NO) then 
          do ICPU=1,NPE 
          IDUM=IBC_P_NO(ICPU,INI,nb,1)
          if(IDUM/=0) then
            IBC_P_TOL(INI,nb,2)=IBC_P_TOL(INI,nb,2)+10
          endif
          enddo
          IBC_P_TOL(INI,nb,2)=
     &    MIN(IBC_P_TOL(INI,nb,2),iSprayInjTotal(INI))
!
          ISUM=0
          ICONT=0
          do ICPU=1,NPE
          IDUM=IBC_P_NO(ICPU,INI,nb,1)
          if(IDUM/=0) then 
            dum1=dble(IBC_P_NO(ICPU,INI,nb,1))
     &          /dble(IBC_P_TOL(INI,nb,1))
            IBC_P_NO(ICPU,INI,nb,2)=INT(dum1*IBC_P_TOL(INI,nb,2))
            ISUM=ISUM+IBC_P_NO(ICPU,INI,nb,2)
            ICONT=ICPU
          endif
          enddo
          if(ICONT/=0) then
           IDUM=(NINALL(INI)-NINALL(INI-1))-ISUM
           IBC_P_NO(ICONT,INI,nb,2)=IBC_P_NO(ICONT,INI,nb,2)+IDUM
          endif
        endif
        enddo
      enddo
!----------------------------------
! --- 
!----------------------------------
      do INI=1,N_injctr 
       if(injkd(INI)==i_coal.or.injkd(INI)==inlet) then 
        nb=bc_nb(INI)   
        IPS=NINALL(INI-1)+1 
        IPE=NINALL(INI)
!
        IPP=IPS
!
        NPSUM(:)=0 
        do 550
        n0=0
        do 551 ICPU=1,NPE 
          IDUM=IBC_P_NO(ICPU,INI,nb,2)
!          
          N0=n0+IDUM
          if(IDUM==0) cycle 
!
          if(NPSUM(ICPU)<=IDUM) then 
            NPSUM(ICPU)=NPSUM(ICPU)+1 
            IPCPU(IPP)=ICPU-1
            IPP=IPP+1
            if(IPP>IPE) goto 500
          endif
 551    enddo
        if(N0==0)   goto 500
 550    enddo 
 500    continue
       else
        IPCPU(:)=root   !??????
       endif
      enddo
!
!------------------------------------------
! --- 
!------------------------------------------
      ifl=iPartLU
!
      open(ifl,FILE=cPTINIFNAM,STATUS="REPLACE",err=20)
      write(ifl,fmt="(
     +        ' X Y Z U V W RHO DIA PARCEL')")
!
      do 1000 INI=1,N_injctr
      write(ifl,fmt="(I10)") NPALL1
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI)
      DIAMIN=1.d0
      DIAMAX=0.d0
      if(injkd(INI)==ihole) then
!---------------
! --- ihole
!---------------
         unit(1)=endx(INI)-dXP0(INI) 
         unit(2)=endy(INI)-dYP0(INI) 
         unit(3)=endz(INI)-dZP0(INI) 
         dum4=sqrt(unit(1)**2+unit(2)**2+unit(3)**2)
         unit(1)=unit(1)/dum4
         unit(2)=unit(2)/dum4
         unit(3)=unit(3)/dum4
         DO IP=IPS,IPE
         if(IPCPU(IP)/=my_rank) cycle
         JOUT(IP)=0
         call utl_random(TEMP1)
         if(dropkd(INI)/=iconst) then
           call drop_distrib(INI,DIA(IP),dSMD(INI),TEMP1)
         elseif(dropkd(INI)==iconst) then
           DIA(IP)=dSMD(INI)
         else

         endif
!         if(.not.lBreakupNONE(INI)) then
!           call utl_random(TEMP1)
!           CALL NUKIYAMA_TANAZAWA_VOLUMEBASE
!     &          (INI,DIA(IP),dSMD(INI),TEMP1)
!         else
!           DIA(IP)=dSMD(INI)
!         endif
         DIAMIN=DMIN1(DIA(IP),DIAMIN)
         DIAMAX=DMAX1(DIA(IP),DIAMAX)
         D3=DIA(IP)*DIA(IP)*DIA(IP)
         VOLUME=PI*D3/6.0D0
         RHOP(IP)=RHOP0(INI)
         PMASS(IP,1)=VOLUME*RHOP(IP)
         PARCEL(IP)=DT0*dInjVolume(INI)/PMASS(IP,1) 
     &   /iSprayInjNum(INI)*iSprayInjInter(INI)
         mass_p_total(INI)=mass_p_total(INI)+PMASS(IP,1)*PARCEL(IP)
!
         call utl_random(TEMP2)
         call utl_random(TEMP3)
         XYZP(IP,1)=dXP0(INI)   !+dSLITWIDTH(INI)*TEMP2*DCOS(TEMP3*2.d0*PI)
         XYZP(IP,2)=dYP0(INI)
         XYZP(IP,3)=dZP0(INI)   !+dSLITWIDTH(INI)*TEMP2*DSIN(TEMP3*2.d0*PI)
!dSPRAYANGLE(INI)*PI/180.d0
         ANGLE1=abs(TEMP3)*2.d0*PI!/180.d0  !alpha
         ANGLE2=abs(TEMP2)*dFANANGLE(INI)*0.5d0*PI/180.d0  !beta
!!+dSPRAYANGLE(INI)*PI/180.d0

!-----------------------------------------------------
!         dum1=dUVWNEW(INI)*DCOS(ANGLE2)
!         dum2=dUVWNEW(INI)*DSIN(ANGLE1)*DSIN(ANGLE2)
!         dum3=dUVWNEW(INI)*DCOS(ANGLE1)*DSIN(ANGLE2)
!-----------------------------------------------------
!         dumm(1)=dUVWNEW(INI)*unit(1)*DCOS(ANGLE2)
!         dumm(2)=dUVWNEW(INI)*unit(2)*DCOS(ANGLE2)
!         dumm(3)=dUVWNEW(INI)*unit(3)*DCOS(ANGLE2)
!-----------------------------------------------------OK
!         UVWP(IP,1)=      dUVWNEW(INI)*DSIN(ANGLE1)*DSIN(ANGLE2)
!         UVWP(IP,2)=-1.d0*dUVWNEW(INI)*DCOS(ANGLE2)
!         UVWP(IP,3)=      dUVWNEW(INI)*DCOS(ANGLE1)*DSIN(ANGLE2)
!-----------------------------------------------------
!         UVWP(IP,1)=dUVWNEW(INI)*DSIN(ANGLE1)*DCOS(ANGLE2)
!         UVWP(IP,2)=-1.d0*dUVWNEW(INI)*DCOS(ANGLE1)*DCOS(ANGLE2)
!         UVWP(IP,3)=-1.d0*dUVWNEW(INI)*DSIN(ANGLE2)
!-----------------------------------------------------
         UVWP(IP,1)=dUVWNEW(INI)*DSIN(ANGLE2)*DCOS(ANGLE1)
         UVWP(IP,2)=dUVWNEW(INI)*DSIN(ANGLE2)*DSIN(ANGLE1)
         UVWP(IP,3)=-1.d0*dUVWNEW(INI)*DCOS(ANGLE2)
!-----------------------------------------------------
         write(ifl,fmt="(I6,9(1X,E12.5))") IP,XYZP(IP,1),XYZP(IP,2),
     +         XYZP(IP,3),UVWP(IP,1),UVWP(IP,2),UVWP(IP,3),RHOP(IP),
     +         DIA(IP),PARCEL(IP) 

         INSIDE(IP)=NCV   !MIN(NCV,LVEDGE(1,ICFL))    
         enddo
      elseif(injkd(INI)==islit) then 
!----------------------
! --- islit 
!----------------------
         dumm(1)=dble(1.d0-dSLITdir(INI))/dble(SML1+1.d0-dSLITdir(INI))
         dumm(2)=dble(2.d0-dSLITdir(INI))/dble(SML1+2.d0-dSLITdir(INI))
         dumm(3)=dble(3.d0-dSLITdir(INI))/dble(SML1+3.d0-dSLITdir(INI))
         dumm(1)=1.d0-dumm(1)
         dumm(2)=1.d0-dumm(2)
         dumm(3)=1.d0-dumm(3)
         unit(1)=endx(INI)-dXP0(INI) 
         unit(2)=endy(INI)-dYP0(INI) 
         unit(3)=endz(INI)-dZP0(INI) 
         dum4=sqrt(unit(1)**2+unit(2)**2+unit(3)**2)
         unit(1)=unit(1)/dum4
         unit(2)=unit(2)/dum4
         unit(3)=unit(3)/dum4
         DO IP=IPS,IPE 
         if(IPCPU(IP)/=my_rank) cycle
         JOUT(IP)=0
         if(dropkd(INI)/=iconst) then
           call drop_distrib(INI,DIA(IP),dSMD(INI),TEMP1)
         elseif(dropkd(INI)==iconst) then
           DIA(IP)=dSMD(INI)
         else

         endif         
!         if(.not.lBreakupNONE(INI)) then
!           call utl_random(TEMP1)
!           CALL NUKIYAMA_TANAZAWA_VOLUMEBASE
!     &          (INI,DIA(IP),dSMD(INI),TEMP1)
!         else
!           DIA(IP)=dSMD(INI)
!         endif
         DIAMIN=DMIN1(DIA(IP),DIAMIN)
         DIAMAX=DMAX1(DIA(IP),DIAMAX)
         D3=DIA(IP)*DIA(IP)*DIA(IP)
         VOLUME=PI*D3/6.0D0
         RHOP(IP)=RHOP0(INI)
         PMASS(IP,1)=VOLUME*RHOP(IP)
         mass_p_total(INI)=mass_p_total(INI)+PMASS(IP,1)*PARCEL(IP)
!
         PARCEL(IP)=DT0*dInjVolume(INI)/PMASS(IP,1) 
     &   /iSprayInjNum(INI)*iSprayInjInter(INI)
!
         call utl_random(TEMP2)
         call utl_random(TEMP3)
         ANGLE1=dFANANGLE(INI)*(TEMP2-0.5d0)*PI/180.d0
         ANGLE2=(dSPANGLEMAG(INI)*(TEMP3-0.5d0)
!     &         +dSPRAYANGLE(INI)
     &         )*PI/180.d0

         XYZP(IP,1)=dXP0(INI)+SPHEIGHT(INI)*DTAN(ANGLE1)*dumm(1)
         XYZP(IP,2)=dYP0(INI)+SPHEIGHT(INI)*DTAN(ANGLE1)*dumm(2)
         XYZP(IP,3)=dZP0(INI)+SPHEIGHT(INI)*DTAN(ANGLE1)*dumm(3)
!
!         XYZP(IP,1)=dXP0(INI)+SPHEIGHT(INI)*DTAN(ANGLE1) 
!         XYZP(IP,2)=dYP0(INI)
!         XYZP(IP,3)=dZP0(INI)
!
         UVWP(IP,1)=dUVWNEW(INI)*DSIN(ANGLE1)*DCOS(ANGLE2)*dumm(1)
     &              +(dUVWNEW(INI)*DCOS(ANGLE1)*DCOS(ANGLE2)
     &              +dUVWNEW(INI)*DSIN(ANGLE2))*unit(1)
!
         UVWP(IP,2)=dUVWNEW(INI)*DSIN(ANGLE1)*DCOS(ANGLE2)*dumm(2)
     &              +(dUVWNEW(INI)*DCOS(ANGLE1)*DCOS(ANGLE2)
     &              +dUVWNEW(INI)*DSIN(ANGLE2))*unit(2)
!
         UVWP(IP,3)=dUVWNEW(INI)*DSIN(ANGLE1)*DCOS(ANGLE2)*dumm(3)
     &              +(dUVWNEW(INI)*DCOS(ANGLE1)*DCOS(ANGLE2)
     &              +dUVWNEW(INI)*DSIN(ANGLE2))*unit(3)

!         UVWP(IP,1)=dUVWNEW(INI)*DSIN(ANGLE1)*DCOS(ANGLE2)
!         UVWP(IP,2)=-1.d0*dUVWNEW(INI)*DCOS(ANGLE1)*DCOS(ANGLE2)
!         UVWP(IP,3)=-1.d0*dUVWNEW(INI)*DSIN(ANGLE2)
!
         INSIDE(IP)=NCV   !MIN(NCV,LVEDGE(1,ICFL))         
         write(ifl,fmt="(I6,9(1X,E12.5))") IP,XYZP(IP,1),XYZP(IP,2),
     +         XYZP(IP,3),UVWP(IP,1),UVWP(IP,2),UVWP(IP,3),RHOP(IP),
     +         DIA(IP),PARCEL(IP)
         enddo
      elseif(injkd(INI)==i_coal.or.injkd(INI)==inlet) then
!----------------------
! --- i_coal,inlet
!----------------------
        do nb=1,nbcnd 
        NO=nobcnd(nb) 
        if(NO_BC(INI)==NO) then 
          NO=0
          exit
        elseif(NO_BC(INI)==0)then
          exit
        endif
        enddo
!
        IF(NO/=0)then 
          write(ifle,'(1X,a,I4.4,a)') 
     &      'ERR: bc_no=',NO,' is NOT existed' 
          call FFRABORT(1,'ERR:bc_no in &particle_model') 
        endif 
!-------------------------
! --- 
!-------------------------
        unit(1)=endx(INI)-dXP0(INI) 
        unit(2)=endy(INI)-dYP0(INI) 
        unit(3)=endz(INI)-dZP0(INI) 
        dum4=sqrt(unit(1)**2+unit(2)**2+unit(3)**2)
        unit(1)=unit(1)/dum4
        unit(2)=unit(2)/dum4
        unit(3)=unit(3)/dum4
!
        ICONT=0    
        INJCPU=0   
!-------------------------
! --- 
!-------------------------
        if(ical_WDR==1) then
          if(dropkd(INI)/=iconst) then 
            write(ifle,*) 'INJECTOR number=' ,INI 
            write(ifle,*) 
     &      "drop_distribution='cnst' for rain model"
            call FFRABORT(1,'ERR: rain model must ') 
          endif
          nuP_m2(INI)=dble(IPE-IPS+1)/nuP_m2(INI)
          NB=bc_nb(INI)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          allocate(nuP_IBFL(IBFE-IBFS+1),stat=ierr1)
          if(ierr1/=0) call FFRABORT 
     & (1,'ERR: nuP_IBFL allocate error in READ_PARTICLE_INI') 
          nuP_IBFL(:)=0.d0
        endif
!
        DO IP=IPS,IPE
        if(IPCPU(IP)/=my_rank) cycle
        JOUT(IP)=0
        call utl_random(TEMP1)
        if(dropkd(INI)/=iconst) then
          call drop_distrib(INI,DIA(IP),dSMD(INI),TEMP1)
        elseif(dropkd(INI)==iconst) then
          DIA(IP)=dSMD(INI)
        else
        endif
!
        D3=DIA(IP)*DIA(IP)*DIA(IP)
        VOLUME=PI*D3/6.0D0
        RHOP(IP)=RHOP0(INI)
        PMASS(IP,1)=VOLUME*RHOP(IP)

        RHOP(IP)=RHOP0(INI)
!
!
        XYZP(IP,1)=dXP0(INI)
        XYZP(IP,2)=dYP0(INI)
        XYZP(IP,3)=dZP0(INI)
!
!
        if(NO_BC(INI)/=0) then
          nb=bc_nb(INI)
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
 3200     call utl_random(TEMP2)
          IBFL=int(IBFS+TEMP2*(IBFE-IBFS))
          if(IBFL<IBFS.or.IBFL>IBFE) then 
            IBFL=IBFS
          endif
          ICFL=LBC_SSF(IBFL) 
          if(ical_WDR==1) then
            IDUM=IBFL-IBFS+1
            dum1=nuP_m2(INI)*SFAREA(4,ICFL)
            if(nuP_IBFL(IDUM)>dum1) then
              goto 3200
            endif
            nuP_IBFL(IDUM)=nuP_IBFL(IDUM)+1.d0

          endif
!
          unit(1)=-SFAREA(1,ICFL) !HIBI
          unit(2)=-SFAREA(2,ICFL)
          unit(3)=-SFAREA(3,ICFL)
!
          if(ical_WDR==1) then
            PARCEL(IP)=pcl_WDR(INI)
            dUVWNEW(INI)=rain_vel(INI)
            UVWP(IP,1)=dUVWNEW(INI)*unit(1)
            UVWP(IP,2)=dUVWNEW(INI)*unit(2)
            UVWP(IP,3)=dUVWNEW(INI)*unit(3)
          else
            PARCEL(IP)=DT0*dInjVolume(INI)/PMASS(IP,1) 
     &      /iSprayInjNum(INI)*iSprayInjInter(INI)
!
            UVWP(IP,1)=dUVWNEW(INI)*unit(1)
            UVWP(IP,2)=dUVWNEW(INI)*unit(2)
            UVWP(IP,3)=dUVWNEW(INI)*unit(3)
          endif
!
          mass_p_total(INI)=mass_p_total(INI)+PMASS(IP,1)*PARCEL(IP) 

          write(ifl,fmt="(I6,9(1X,E12.5))") IP,XYZP(IP,1),XYZP(IP,2),
     +     XYZP(IP,3),UVWP(IP,1),UVWP(IP,2),UVWP(IP,3),RHOP(IP),
     +     DIA(IP),PARCEL(IP)


!
          unitt(1)=SFAREA(1,ICFL)
          unitt(2)=SFAREA(2,ICFL)
          unitt(3)=SFAREA(3,ICFL)
!
          call utl_random(TEMP1)
          call utl_random(TEMP2)
          call utl_random(TEMP3)
!
          dum0=TEMP1**2+TEMP2**2+TEMP3**2+SML 
          TEMP1=TEMP1/dum0
          TEMP2=TEMP2/dum0
          TEMP3=TEMP3/dum0
!
          dum0=SQRT(SFAREA(4,ICFL)/PI)
          dum1=dum0*TEMP1
          dum2=dum0*TEMP2
          dum3=dum0*TEMP3
!
          dum0=unitt(1)*dum1+unitt(2)*dum2+unitt(3)*dum3 
          dum1=dum1-dum0*unitt(1)
          dum2=dum2-dum0*unitt(2)
          dum3=dum3-dum0*unitt(3)
!
          XYZP(IP,1)=SFCENT(1,ICFL)+dum1
          XYZP(IP,2)=SFCENT(2,ICFL)+dum2
          XYZP(IP,3)=SFCENT(3,ICFL)+dum3
          INSIDE(IP)=MIN(NCV,LVEDGE(1,ICFL))
          if(IP==IPE) then
            dumm(1)=SFCENT(1,ICFL)
            dumm(2)=SFCENT(2,ICFL)
            dumm(3)=SFCENT(3,ICFL)
            INJCPU=my_rank
            ICONT=ICONT+1
          endif
        else
        endif
        enddo
!
        if(NINALL(INI)-NINALL(INI-1)==0) cycle
        call hpcisum(INJCPU)
        if(INJCPU>NPE) then
          write(ifll,'(2X,I4)') INJCPU
          call FFRABORT(1,'ERR-2: READ_PARTICLE_INI')
        endif
!
        call hpcrbcast_ary(INJCPU,dumm,3)
        dXP0(INI)=dumm(1)
        dYP0(INI)=dumm(2)
        dZP0(INI)=dumm(3)
!
        call hpcisum(ICONT)
        if(ICONT/=1) then
          write(ifll,'(2X,a,I4)') 'ERR-1: Call your supporter, ICONT=', 
     &                          ICONT
          call FFRABORT(1,'ERR: READ_PARTICLE_INI') 
        endif
        if(ical_WDR==1) deallocate(nuP_IBFL)  
      elseif(injkd(INI)==user) then  
!-------------------- 
! --- imanual 
!-------------------- 
        call USER_PARTICLE_INITIAL(INI,IPS,IPE, 
     &  NCV,MXALLCV,MXPRT,MXCOMP,MXCVFAC,MXSSFBC, 
     &  SFCENT,LBC_SSF,SFAREA,LVEDGE,CVCENT, 
     *  DIA,XYZP,UVWP,RHOP,PARCEL,DT0,TMPP,PYS, 
     &  PMASS,INSIDE,JOUT,HIS_P)  
      elseif(injkd(INI)==imanual) then
!--------------------
! --- imanual
!--------------------
        IP=1
        do ICV=1,NCV
        if(mod(ICV,3)==0) then
          IP=IP+1
          if(IPCPU(IP)/=my_rank) cycle
          if(IP>MXPRT) exit
          XYZP(IP,1)=CVCENT(1,ICV)
          XYZP(IP,2)=CVCENT(2,ICV)
          XYZP(IP,3)=CVCENT(3,ICV)
          JOUT(IP)=0
          DIA(IP)=dSMD(INI)
          D3=DIA(IP)*DIA(IP)*DIA(IP)
          VOLUME=PI*D3/6.0D0
          RHOP(IP)=RHOP0(INI)
          PMASS(IP,1)=VOLUME*RHOP(IP)
          PARCEL(IP)=1.d0
          UVWP(IP,1)=0.d0
          UVWP(IP,2)=0.d0
          UVWP(IP,3)=0.d0
          INSIDE(IP)=NCV
        endif
        enddo
      endif
!
!----------------------------------------------
! --- initial temperature 
!----------------------------------------------
!
      if(ical_t) then
        DO IP=IPS,IPE
        TMPP(IP,1)=Tini_p(INI)
        enddo
      ELSE
        DO IP=IPS,IPE
        TMPP(IP,1)=300.D0
        enddo
      endif
!
      if((encomp>=1.or.ical_flmlt)) then
        if(ifuel(INI)==icoal.or.
     &     ifuel(INI)==iglass.or.
     &     ifuel(INI)==isolid) then
          DO IP=IPS,IPE
          do icom=1,MXCOMP
          PYS(IP,icom,1)=particle_comp(INI,icom)
          enddo
          enddo
        elseif(ifuel(INI)==iliquid) then
          do icom=1,MXCOMP
          DO IP=IPS,IPE
          PYS(IP,icom,1)=particle_comp(INI,icom)
          enddo
          enddo
        endif
      elseif(calxi.and.idifflm/=ORG) then
        DO IP=IPS,IPE
        do icom=1,MXCOMP
        PYS(IP,icom,1)=particle_comp(INI,icom)
        enddo
        enddo
      endif
!
      if(IFLAG_COAL==icoal) then
        PMASSC(:)=PYS(:,II_C,1)*PMASS(:,1)
      endif

      close(ifl,err=200)
C------------------------------------------
C  Particle Initial Value Setting End.
C  Initial Data Output
C------------------------------------------
      dFLOWRATE=dInjVolume(INI)/1000.d0
      write(ifll,4020)
      write(ifll,4000)
      write(ifll,4001) N_injctr,INI
      write(ifll,4020)
      write(ifll,fmt="(4X,'Initial Velocity   :',E12.4)") dUVWNEW(INI)
      write(ifll,fmt="(4X,'SauterMeanDiameter :',E12.4)") dSMD(INI)
      write(ifll,fmt="(4X,'Dens.Liquid[kg/m^3]:',E12.4)") RHOP0(INI)
      write(ifll,fmt="(4X,'Maximum Diameter   :',E12.4)") DIAMAX
      write(ifll,fmt="(4X,'Minimum Diameter   :',E12.4)") DIAMIN
      write(ifll,fmt="(4X,'Num of Total Parcel:',I6)") NINALL(INI)
      write(ifll,fmt="(4X,'Total Inject mass rate [kg/s]:',E12.4)")
     &             dInjVolume(INI)
      write(ifll,fmt="(4X,'Flow Rate[kg/msec]:',E12.4)") dFLOWRATE
      write(ifll,fmt="(4X,'Start of Injection :',I5)") 
     &    iSprayInjStart(INI)
      write(ifll,fmt="(4X,'End of Injection   :',I5)") iSprayInjEnd(INI)
      write(ifll,fmt="(4X,'Injecting Steps    :',I5)") 
     &    iSprayInjTime(INI)
      write(ifll,fmt="(4X,'Inject Pc at a STEP:',I5)") iSprayInjNum(INI)
      write(ifll,fmt="(4X,'Injection Point X  :',E12.4)")dXP0(INI)
      write(ifll,fmt="(4X,'Injection Point Y  :',E12.4)")dYP0(INI)
      write(ifll,fmt="(4X,'Injection Point Z  :',E12.4)")dZP0(INI)
      write(ifll,fmt="(4X,'Breakup Model      :',A4)") 
     &    cBreakupModel(INI)
      write(ifll,fmt="(4X,'Num of Sub-cycling :',I5)") ITER_INNER
      write(ifll,'(4X,a)') '------ NUKIYAMA-TANAZAWA CONDITIONS ------'
      write(ifll,fmt="(4X,'ALPHA(Dist. Funct.):',E12.4)") dAlpha(INI)
      write(ifll,fmt="(4X,'BETA (Dist. Funct.):',E12.4)") dBeta(INI)
      write(ifll,'(4X,a)') '------- SPRAY GEOMETORY CONDITIONS -------'
      write(ifll,fmt="(4X,'NOZZLE TYPE        :',A4)") cNozzleType(INI)
      write(ifll,fmt="(4X,'Width (Diameter)[m]:',E12.4)") 
     &    dSLITWIDTH(INI)
      write(ifll,fmt="(4X,'FAN(CONE)ANGLE[deg]:',E12.4)") dFANANGLE(INI)
      write(ifll,fmt="(4X,'SPREADANGLE   [deg]:',E12.4)") 
     &    dSPANGLEMAG(INI)
      write(ifll,fmt="(4X,'INJECT AXIS ANGLE  [deg]:',E12.4)") 
     &                 dSPRAYANGLE(INI)
      if(.not.lBreakupNONE(INI)) then
        write(ifll,'(4X,a)') '------ TAB MODEL CONDITIONS ------'
        write(ifll,fmt="(4X,A4,E12.4)") "lmu :",dLiquidMu(INI)
        write(ifll,fmt="(4X,A4,E12.4)") "stns:",dLiquidStens(INI)
        write(ifll,fmt="(4X,A4,E12.4)") "Cf  :",dCf(INI)
        write(ifll,fmt="(4X,A4,E12.4)") "Ck  :",dCk(INI)
        write(ifll,fmt="(4X,A4,E12.4)") "Cb  :",dCb(INI)
        write(ifll,fmt="(4X,A4,E12.4)") "Cd  :",dCd(INI)
        write(ifll,fmt="(4X,A4,E12.4)") "Cv  :",dCv(INI)
        write(ifll,fmt="(4X,A4,E12.4)") "K   :",dK(INI)
      endif

 1000 ENDDO

!
   50 FORMAT(6E12.4)
 4020 FORMAT(2X,'|',108('='),'|')      
 4000 FORMAT(2X,'|',20X,'Particle Initial Condition',61X,'|')
 4001 FORMAT(2X,'|',20X,
     &  'Total Injector number=',I8,20X,'Inj. No=',I8,21X,'|')
 4002 FORMAT(2X,'|',20X,'Injector number: ',I8,20X,'|')
!
C-------------------------------------------------- 
C --- Particle Initial Data output end
C--------------------------------------------------
C --- Check NUKIYAMA-TANAZAWA Dist. if you need.
C--------------------------------------------------
!
      RETURN
!
   10 CALL FFRABORT(1,'ERR: Opening Particle_info.txt')
   20 CALL FFRABORT(1,'ERR: Opening '//trim(cPTINIFNAM))
  100 CALL FFRABORT(1,'ERR: Exceed Max Particle Number')
  200 CALL FFRABORT(1,'ERR: Closing '//trim(cPTINIFNAM))
!
      END SUBROUTINE READ_PARTICLE_INI
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE RE_INJ_P(IP,INI,MXPRT,PI,MXCOMP,MXCVFAC,MXSSFBC,
     &  SFCENT,LBC_SSF,SFAREA,JOUT,IPCPU,
     *  DIA,XYZP,UVWP,RHOP,PARCEL,NDIS,DT0,TMPP,PYS,
     &  PMASS,HIS_P,INSIDE,LVEDGE)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      USE module_dimension,only : NCV
      use module_io,only      : ifle,ifll
      USE module_hpcutil,only : my_rank,root
      use module_model,  only : ical_t,encomp
      use module_particle,only: Tini_p,Cp_pa,heat_f,eps_rad,
     &                          Tref_p,bc_nb,mass_p_total
      use module_particle,only: ifuel,icoal,iliquid,iglass,isolid
      USE module_particle,only: IFLAG_COAL,icoal,dUVWNEW,NO_BC,imanual,
     &                          particle_comp,VS_CHO,VS_H2O,VS_H2O,
     &                          VS_HCN,net_MF_HCN,Q_MF_vola,
     &                          iSprayInjInter,iSprayInjInter,
     &                          iSprayInjNum,dSMD,dXP0,dYP0,dZP0,
     &                          endx,endy,endz,RHOP0,injkd,NINALL,
     &                          N_injctr,NPALL,ihole,islit,i_coal,
     &                          Q_MF_H2O,dInjVolume,inlet
     &                          
      use module_constant,only: SML
      use module_trace ,only  : WK1_P
      use module_boundary,only : nobcnd,LBC_INDEX
      use module_scalar,only   : ical_flmlt
!
      IMPLICIT NONE
!
!--------------------
! Dummy Auguments 
!--------------------
!
      INTEGER, INTENT(IN)   :: MXPRT,MXCOMP,MXCVFAC,MXSSFBC
      INTEGER, INTENT(IN)   :: NDIS,IP,INI
      REAL*8, INTENT(IN)    :: PI,DT0
      real*8 ,intent(in)    :: SFCENT(3,MXCVFAC)
      INTEGER ,intent(in)   :: LBC_SSF( MXSSFBC)
      REAL*8, INTENT(INOUT) :: DIA     (MXPRT)
      REAL*8, INTENT(INOUT) :: XYZP    (MXPRT,3)
      REAL*8, INTENT(INOUT) :: UVWP    (MXPRT,3)
      REAL*8, INTENT(INOUT) :: RHOP    (MXPRT)
      REAL*8, INTENT(INOUT) :: PARCEL  (MXPRT)
      REAL*8,INTENT(INOUT)  :: PYS(     MXPRT,MXCOMP,2)
      REAL*8,INTENT(INOUT)  :: TMPP  (  MXPRT,2)
      REAL*8,INTENT(INOUT)  :: PMASS (  MXPRT)
      REAL*8,INTENT(INOUT)  :: HIS_P (  MXPRT)
      REAL*8 ,INTENT(IN)    :: SFAREA(4,MXCVFAC)
      integer,intent(inout) :: JOUT(    MXPRT)
      integer,intent(inout) :: IPCPU(    MXPRT)
      INTEGER,INTENT(INOUT) :: INSIDE(  MXPRT)
      INTEGER,INTENT(IN)    :: LVEDGE(  2,MXCVFAC)
!
!--------------------
! --- Local Entities 
!--------------------
!
      REAL*8  :: D3,dum1,dum2,dum3,dum4
      REAL*8  :: VOLUME
      REAL*8  :: unitt(3),dum0,unit(3)
      REAL*8,save  :: TEMP1,TEMP2,TEMP3
      INTEGER :: icom,IBFS,IBFE,nb,IBFL,ICFL
!
! --- 
!
!
      RETURN
!
      END SUBROUTINE RE_INJ_P


C
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BREAKUP_NORMAL_TAB(INI,IT,ITS,
     1      DT0,GRDC,N,INSIDE,
     &      RHOP,DIA,XYZP,UVWP0,
     2      PARCEL,VEL,RHO,RMU,CVCENT,lBreakup,dTbu,dDiamReduc,
     3      CVVOLM,aks)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! [DUMMY ARGUMENTS]
!

      USE module_dimension
      use module_io,only       : ifle,ifll
      USE module_particle ,ONLY : dPDstion,
     &                            dPDstvel,
     &                            dLiquidMu,
     &                            dLiquidStens,dCf,dCk,dCd,dCb,dK,dCv,
     &                            N_injctr
!
      IMPLICIT NONE

! Dummy Auguments
!
      INTEGER,INTENT(IN)  :: N,INI
      INTEGER,INTENT(IN)  :: IT,ITS
      REAL*8,INTENT(IN)   :: DT0
      REAL*8,INTENT(IN)   :: RHO (      MXALLCV,2)
      REAL*8,INTENT(IN)   :: RHOP(      MXPRT)
      REAL*8,INTENT(INOUT):: DIA (      MXPRT) 
      REAL*8,INTENT(INOUT):: UVWP0(3)
      REAL*8,INTENT(INOUT):: PARCEL(    MXPRT)
      REAL*8,INTENT(IN)   :: XYZP(      MXPRT,3)
      REAL*8,INTENT(IN)   :: RMU (      MXALLCV)
      REAL*8,INTENT(IN)   :: VEL (      MXALLCV,3,2)
      REAL*8,intent(in)   :: CVCENT(3  ,MXALLCV)
      REAL*8,intent(in)   :: GRDC(      MXALLCV,3,3)
      INTEGER,intent(in)  :: INSIDE(    MXPRT)
      LOGICAL,intent(out) :: lBreakup
      REAL*8, intent(out) :: dTbu
      REAL*8, intent(out) :: dDiamReduc
      REAL*8, INTENT(IN)  :: CVVOLM(    MXALLCV)
      real*8 ,intent(in)  :: aks(       MXALLCVR,MXRANS)
!
! Local Entities
!
      REAL*8 :: dWeber, dWec, dOmega, dTdrcp,dVolume,dY0,dYd0
      REAL*8 :: dAmp,dCOSwdt,dSINwdt,dOMGrcp,dPhi
      REAL*8 :: dTemp1,dTemp2,dTemp3,dDTbu
      REAL*8 :: US,VS,WS,SLIPV,dDIANEW
      REAL*8 :: UVW0(3),UVWFluc(3),Vnorm
      INTEGER:: KE, i
      integer:: iIncrement(2)
      REAL*8, parameter :: PI=3.14159265358979323846264
      LOGICAL,save      :: lCallFlag=.FALSE.
!
      if(.not.lCallFlag) then
        write(ifll,*) "Subroutine BREAKUP_NORMAL_TAB called"
        lCallFlag = .TRUE.
      endif

      KE=INSIDE(N)
      lBreakup=.FALSE.
      dDiamReduc = 1.d0
C----------------------------
C Calculate Slip Velocity
C----------------------------
      UVW0(1)= VEL(KE,1,2)
     *            +GRDC(KE,1,1)*(XYZP(N,1)-CVCENT(1,KE))
     *            +GRDC(KE,2,1)*(XYZP(N,2)-CVCENT(2,KE))
     *            +GRDC(KE,3,1)*(XYZP(N,3)-CVCENT(3,KE))
      UVW0(2)= VEL(KE,2,2)
     *            +GRDC(KE,1,2)*(XYZP(N,1)-CVCENT(1,KE))
     *            +GRDC(KE,2,2)*(XYZP(N,2)-CVCENT(2,KE))
     *            +GRDC(KE,3,2)*(XYZP(N,3)-CVCENT(3,KE))
      UVW0(3)= VEL(KE,3,2)
     *            +GRDC(KE,1,3)*(XYZP(N,1)-CVCENT(1,KE))
     *            +GRDC(KE,2,3)*(XYZP(N,2)-CVCENT(2,KE))
     *            +GRDC(KE,3,3)*(XYZP(N,3)-CVCENT(3,KE))
C----------------------------------------------------------------
C Add Gaussian Noise to Velocity components (Random Walk Model)
C----------------------------------------------------------------
      CALL add_random_walk(N,KE,GRDC,CVVOLM,UVWFluc,aks)
!
      UVW0(1)=UVW0(1)+UVWFluc(1)
      UVW0(2)=UVW0(2)+UVWFluc(2)
      UVW0(3)=UVW0(3)+UVWFluc(3)
C
C
      US=UVWP0(1)-UVW0(1)
      VS=UVWP0(2)-UVW0(2)
      WS=UVWP0(3)-UVW0(3)
C      
      SLIPV=SQRT(US*US + VS*VS + WS*WS)
C--------------------------------
C Calculate TAB Model Constants
C--------------------------------
      dY0=dPDstion(N)
      dYd0=dPDstvel(N)
      dWeber=RHO(KE,2)*(SLIPV*SLIPV)*DIA(N)*0.5d0/dLiquidStens(INI)
      dWec=dCf(INI)*dWeber/(dCk(INI)*dCb(INI))
      dTdrcp=2.d0*dCd(INI)*dLiquidMu(INI)/(RHOP(N)*DIA(N)*DIA(N))
      dOmega=8.d0*dCk(INI)*
     &   dLiquidStens(INI)/(RHOP(N)*DIA(N)*DIA(N)*DIA(N))
     1         -dTdrcp*dTdrcp
C----------------------------------------------------------------
C If the square of Omega is negative, distortion is negligible.
C----------------------------------------------------------------
      if(dOmega .gt. 0.d0) then
        dOmega=DSQRT(dOmega)
      else 
	return
      endif
C----------------------
C Preparation
C----------------------
      dCOSwdt = DCOS(dOmega*DT0)
      dSINwdt = DSIN(dOmega*DT0)
      dOMGrcp = 1.d0/dOmega
C--------------------------------------------------------
C Calculate time advanced Nondimensional Distortion.
C--------------------------------------------------------
      dPDstion(N) = dWec+DEXP(-1.d0*DT0*dTdrcp)
     1                   *((dY0-dWec)*dCOSwdt
     2                     +dOMGrcp*(dYd0+(dY0-dWec)*dTdrcp)
     3                      *dSINwdt)
C-----------------------------------------------------------------
C Calculate time advanced Nondimensional Distortion Velocity.
C-----------------------------------------------------------------
      dPDstvel(N) = (dWec-dPDstion(N))*dTdrcp
     1             + dOmega*DEXP(-1.d0*DT0*dTdrcp)
     2               *(dOMGrcp*(dYd0+(dY0-dWec)*dTdrcp)*dCOSwdt
     3                 -(dY0-dWec)*dSINwdt)

C------------------------------------------------
C Calculate Amplitude for Undumped Oscillation
C And if dWec + dAmp < 1, breakup is impossible.
C------------------------------------------------
      dAmp=DSQRT((dY0-dWec)**2+(dYd0*dOMGrcp)**2)
      if((dWec+dAmp).le.1.d0) return
      if((dWec-1.d0).gt.dAmp) THEN
        dTbu=0.5d0*DT0
        goto 150
      endif
!
      if(dPDstion(N) .ge. 1.0) THEN
        dTbu=0.5d0*DT0
        goto 150
      endif
C--------------------------------
C Calculate Breakup Time dTbu
C--------------------------------
      dPhi=DATAN(dOMGrcp*dYd0/(dY0-dWec)) 
      dTbu=(DACOS((dWec-1.d0)/dAmp)+dPhi)*dOMGrcp 
      if(dTbu.lt.0.d0) THEN 
        write(ifll,*) "Error!! Negative Breakup time."
        write(ifll,1000) "Wec,Omega,Amp,Phi,:",dWec, dOmega, dAmp, dPhi
 1000   format(A19,4(E12.5))
      endif
 150  continue
C-------------------------------------------------
C --- write(*,*) "breakuptime:",dTbu
C-------------------------------------------------
C Judge breakup or not, and treatment for breakup
C-------------------------------------------------
      if(dTbu<DT0) THEN
C---------------------
C Preparation
C---------------------
        dCOSwdt = DCOS(dOmega*dTbu)
        dSINwdt = DSIN(dOmega*dTbu)
C---------------------------------------------------------------
C Calculate Nondimensional Distortion at Breakup.
C---------------------------------------------------------------
        dPDstion(N) = dWec+DEXP(-1.d0*dTbu*dTdrcp)
     1                     *((dY0-dWec)*dCOSwdt
     2                       +dOMGrcp*(dYd0+(dY0-dWec)*dTdrcp)
     3                        *dSINwdt)
C---------------------------------------------------------------
C Calculate Nondimensional Distortion Velocity at Breakup.
C---------------------------------------------------------------
        dPDstvel(N) = (dWec-dPDstion(N))*dTdrcp
     1               + dOmega*DEXP(-1.d0*dTbu*dTdrcp)
     2                 *(dOMGrcp*(dYd0+(dY0-dWec)*dTdrcp)*dCOSwdt
     3                   -(dY0-dWec)*dSINwdt)
C---------------------------------------------------------------
C Calculate New Diameter by Energy Conservation
C---------------------------------------------------------------
        dDiamReduc = 1.d0+8.d0*dK(INI)/20.d0
     1    +RHOP(N)*DIA(N)*DIA(N)*DIA(N)*dPDstvel(N)*dPDstvel(N)
     2     /(8.d0*dLiquidStens(INI))*(6.d0*dK(INI)-5.d0)/120.d0
!---------------------------------------------------------------
C +++++++++++++++ Add Normal Velocity +++++++++++++++
!---------------------------------------------------------------
        call utl_random(dTemp1)
        call utl_random(dTemp2)
        dTemp1 = dTemp1*2.d0*PI
        dTemp2 = dTemp2*2.d0*PI
C
        Vnorm =  dCv(INI)*dCb(INI)*DIA(N)*0.5d0*dPDstvel(N)
        UVWP0(1)=UVWP0(1)+Vnorm*DSIN(dTemp1)*DCOS(dTemp2)
        UVWP0(2)=UVWP0(2)+Vnorm*DCOS(dTemp1)*DCOS(dTemp2)
        UVWP0(3)=UVWP0(3)+Vnorm*DSIN(dTemp1)
!
        if(dDiamReduc.lt.1.0) then
          write(ifll,*) "Err. Child Particle is larger than parent"
        endif
        if(dDiamReduc.lt.0.0) then
          write(ifll,*) "Err. Negative Diameter"
          call FFRABORT(1,'ERR: Normal TAB')
        endif
        lBreakup=.TRUE.
      END IF
C
C ++++++++++++++++ Breakup Treatment End ++++++++++++++++
C
      END SUBROUTINE BREAKUP_NORMAL_TAB
C
C
C Read Restart file for Particle Calculation
! 
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE read_particle_restart(
     &  IT,time,NPALL,INSIDE,JOUT,NINALL,NSTA,NRES,
     &  PMASS,DIA,XYZP,UVWP,PARCEL,RHOP,TMPP,
     &  HEIGHT,FU,FV,FW,PYS,HIS_P,
     &  ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension 
      use module_io,only       : ifll,ifle,getfil
      use module_hpcutil, only : NPE,my_rank,INITin_P,root
      use module_particle,only : iPartLU,lBreakupNONE
      use module_particle,only : dPDstion,dPDstvel,N_injctr,
     &                           IFLAG_COAL,icoal,iglass,
     &                           VS_CHO,V_CHO,V_H2O,VS_H2O,
     &                           VS_HCN,V_HCN
      use module_time,    only : iters
      use module_trace ,only   : WK1_P,IWK_P
!
!  1. Write restart file
!
      implicit none
!
!--- [dummy arguments]
!
      integer,intent(in)   :: IT
      real*8 ,intent(in)   :: time
      integer,intent(INOUT)  :: NPALL (N_injctr),
     %                          NINALL(0:N_injctr),
     &                          NSTA(N_injctr),
     &                          NRES(0:N_injctr)
      INTEGER,INTENT(INOUT):: INSIDE(        MXPRT)
      INTEGER,INTENT(INOUT):: JOUT  (        MXPRT)
      REAL*8 ,INTENT(INOUT):: PMASS(         MXPRT,2)
      REAL*8 ,INTENT(INOUT):: PARCEL(        MXPRT)
      REAL*8 ,INTENT(INOUT):: DIA(           MXPRT)
      REAL*8 ,INTENT(INOUT):: RHOP(          MXPRT)
      REAL*8 ,INTENT(INOUT):: XYZP(        MXPRT,3)
      REAL*8 ,INTENT(INOUT):: UVWP(        MXPRT,3)
      REAL*8 ,INTENT(INOUT):: HEIGHT(    MEP,MXPRT)
      REAL*8 ,INTENT(INOUT):: FU(          MXPRT,2)
      REAL*8 ,INTENT(INOUT):: FV(          MXPRT,2)
      REAL*8 ,INTENT(INOUT):: FW(          MXPRT,2)
      REAL*8,INTENT(INOUT) :: PYS(  MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT) :: TMPP (       MXPRT)
      REAL*8,INTENT(INOUT) :: HIS_P(       MXPRT)
      
      integer,intent(OUT):: ierror
!
!--- [local entities]
!
      character(LEN=128) :: fnam
      character(len=8)   :: text
      character(LEN=8)   :: vname
      integer :: I,J,ifl=0,IP,iterr,INI,INII,IPS,
     & IPE,IDUM,IDUMX,IPSX,IPEX
      real*8  :: timer
      logical :: inifil
      integer,allocatable :: NPALLX(:),NINALLX(:),NSTAX(:),NRESX(:)
      integer,save :: iflag=0,ierr1=0
!
!
      ALLOCATE(NPALLX(N_injctr),
     &         NINALLX(0:N_injctr),
     &         NSTAX(N_injctr),
     &         NRESX(N_injctr))
      NPALLX=0
      NINALLX=0
      NSTAX=0
      NRESX=0
!
      ierror=0
!
      if(ifl==0) then
        call getfil(ifl,fnam,'p_initial')
      endif
      inifil=.false.
      if(NPE>1) then
         fnam=trim(INITin_P)!//trim(fnam)
      else
      if(fnam.eq.' ') then
        write(ifle,*) '### error : p_initial file name not defined'
        write(ifle,*) 
     &   'lack of initial_particle file name in control file'
        call FFRABORT(1,'ERR:')
      elseif(fnam.ne.' ') then
        inquire(file=fnam,exist=inifil)
        if(.not.inifil) then
          write(ifle,'(1x,a)') 
     &      'ERR: particle restart file reading error'
          write(ifle,'(1x,2a)') 
     &      'lack of initial_particle file in work Directory',
     &       TRIM(adjustl(fnam))
          write(ifle,'(1x,2a)') 
     &     'MSG: or: changing ',
     &     '&injector_number/total_start > &time/start'
          call FFRABORT(1,'ERR:')
        endif
      endif
      endif
!
!      if(NPE>1) then
!        fnam=INITin_P
!      endif
!
      open(ifl,file=trim(fnam),form='unformatted',
     &     status='OLD',err=100)
!
      do INI=1,N_injctr
      read(ifl) iterr,timer,INII,NPALLX(INI),NINALLX(INI),
     &          NSTAX(INI),NRESX(INI)
!
      IDUMX=NINALLX(INI)-NINALLX(INI-1)
      IDUM=NINALL(INI)-NINALL(INI-1)
!
      if(IDUMX/=IDUM) then
        iflag=1
      endif
      NPALL(INI)=NPALLX(INI)
      NSTA(INI)=NSTAX(INI)
      NRES(INI)=NRESX(INI)
!
!
      enddo      
!
      if(iterr/=iters) then
        write(ifll,'(1X,a,I4,I4)')
     & 'START TIME STEP /= RESTART FILE START STEP',iterr,iters
        call FFRABORT(1,'ERR: TIME-STEP NOT MACHED')
      endif
!
         

      do INI=1,N_injctr
      if(NPE*my_rank.eq.0) then
        write(ifll,*) "MSG:[read_particle_restart]:file:",trim(fnam)
        write(ifll,*) "MSG:[read_particle_restart]:IT=",iterr
        write(ifll,*) "MSG:[read_particle_restart]:time=",timer
        write(ifll,*) "MSG:[read_particle_restart]:NPALL=",NPALL(INI)
      endif
      enddo
!
      IDUM=MAX(MEP,NCOMP,3)
      allocate(WK1_P(MXPRT,IDUM),IWK_P(MXPRT),stat=ierr1)
      if(ierr1/=0) 
     &   call FFRABORT(1,'ERR: allocate error in particle restart')
!
      do INI=1,N_injctr
      IPSX=NINALLX(INI-1)+1
      IPEX=NINALLX(INI-1)+NPALLX(INI)
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI)
      IDUMX=NPALLX(INI)-1       !NINALLX(INI)-NINALLX(INI-1)
      IDUM=NINALL(INI)-NINALL(INI-1)
!
      WK1_P(:,1:3)=0.d0
      read(ifl) vname
      read(ifl) ((WK1_P(IP,J),IP=IPSX,IPEX),J=1,3)
      if(iflag==1) then
        do J=1,3
        XYZP(IPS:IPS+IDUMX,J)=WK1_P(IPSX:IPSX+IDUMX,J)
        enddo
      else
        do J=1,3
        XYZP(IPS:IPS+IDUMX,J)=WK1_P(IPS:IPS+IDUMX,J)
        enddo
      endif
!
      WK1_P(:,1:3)=0.d0
      read(ifl) vname
      read(ifl) ((WK1_P(IP,J),IP=IPSX,IPEX),J=1,3)
      if(iflag==1) then
        do J=1,3
        UVWP(IPS:IPS+IDUMX,J)=WK1_P(IPSX:IPSX+IDUMX,J)
        enddo
      else
        do J=1,3
        UVWP(IPS:IPS+IDUMX,J)=WK1_P(IPS:IPS+IDUMX,J)
        enddo
      endif
!
      WK1_P(:,1)=0.d0
      read(ifl) vname
      read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
      if(iflag==1) then
        DIA(IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,1)
      else
        DIA(IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,1)
      endif
!
      WK1_P(:,1)=0.d0
      read(ifl) vname
      read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
      if(iflag==1) then
        PARCEL(IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,1)
      else
        PARCEL(IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,1)
      endif
!
      WK1_P(:,1)=0.d0
      read(ifl) vname
      read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
      if(iflag==1) then
        RHOP(IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,1)
      else
        RHOP(IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,1)
      endif
!
      WK1_P(:,1)=0.d0
      read(ifl) vname
      read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
      if(iflag==1) then
        PMASS(IPS:IPS+IDUMX,1)=WK1_P(IPSX:IPSX+IDUMX,1)
      else
        PMASS(IPS:IPS+IDUMX,1)=WK1_P(IPS:IPS+IDUMX,1)
      endif
!
      WK1_P(:,1)=0.d0
!      read(ifl) vname
      read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
      if(iflag==1) then
        PMASS(IPS:IPS+IDUMX,2)=WK1_P(IPSX:IPSX+IDUMX,1)
      else
        PMASS(IPS:IPS+IDUMX,2)=WK1_P(IPS:IPS+IDUMX,1)
      endif
!
      WK1_P(:,1)=0.d0
      read(ifl) vname
      read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
      if(iflag==1) then
        HIS_P(IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,1)
      else
        HIS_P(IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,1)
      endif
!
      IWK_P(:)=0
      read(ifl) vname
      read(ifl) (IWK_P(IP),IP=IPSX,IPEX)
      if(iflag==1) then
        INSIDE(IPS:IPS+IDUMX)=IWK_P(IPSX:IPSX+IDUMX)
      else
        INSIDE(IPS:IPS+IDUMX)=IWK_P(IPS:IPS+IDUMX)
      endif
!
      WK1_P(:,1:MEP)=0.d0
      read(ifl) vname
      read(ifl) ((WK1_P(IP,J),J=1,MEP),IP=IPSX,IPEX)
      if(iflag==1) then
        do J=1,MEP
        HEIGHT(J,IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,J)
        enddo
      else
        do J=1,MEP
        HEIGHT(J,IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,J)
        enddo
      endif
!
      WK1_P(:,1:2)=0.d0
      read(ifl) vname
      read(ifl) ((WK1_P(IP,J),IP=IPSX,IPEX),J=1,2)
      if(iflag==1) then
        do J=1,2
        FU(IPS:IPS+IDUMX,J)=WK1_P(IPSX:IPSX+IDUMX,J)
        enddo
      else
        do J=1,2
        FU(IPS:IPS+IDUMX,J)=WK1_P(IPS:IPS+IDUMX,J)
        enddo
      endif
!
      WK1_P(:,1:2)=0.d0
      read(ifl) vname
      read(ifl) ((WK1_P(IP,J),IP=IPSX,IPEX),J=1,2)
      if(iflag==1) then
        do J=1,2
        FV(IPS:IPS+IDUMX,J)=WK1_P(IPSX:IPSX+IDUMX,J)
        enddo
      else
        do J=1,2
        FV(IPS:IPS+IDUMX,J)=WK1_P(IPS:IPS+IDUMX,J)
        enddo
      endif
!
      WK1_P(:,1:2)=0.d0
      read(ifl) vname
      read(ifl) ((WK1_P(IP,J),IP=IPSX,IPEX),J=1,2)
      if(iflag==1) then
        do J=1,2
        FW(IPS:IPS+IDUMX,J)=WK1_P(IPSX:IPSX+IDUMX,J)
        enddo
      else
        do J=1,2
        FW(IPS:IPS+IDUMX,J)=WK1_P(IPS:IPS+IDUMX,J)
        enddo
      endif
!
      if(.not.lBreakupNONE(INI)) then
        WK1_P(:,1)=0.d0
        read(ifl,err=200) vname
        read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
        if(iflag==1) then
          dPDstion(IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,1)
        else
          dPDstion(IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,1)
        endif
!
        WK1_P(:,1)=0.d0
        read(ifl) vname
        read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
        if(iflag==1) then
          dPDstvel(IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,1)
        else
          dPDstvel(IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,1)
        endif
      endif
!
      WK1_P(:,1:NCOMP)=0.d0
      read(ifl) ((WK1_P(IP,J),IP=IPSX,IPEX),J=1,NCOMP)
      if(iflag==1) then
        do J=1,NCOMP
        PYS(IPS:IPS+IDUMX,J)=WK1_P(IPSX:IPSX+IDUMX,J)
        enddo
      else
        do J=1,NCOMP
        PYS(IPS:IPS+IDUMX,J)=WK1_P(IPS:IPS+IDUMX,J)
        enddo
      endif
!
      WK1_P(:,1)=0.d0
      read(ifl) (WK1_P(IP,1),IP=IPSX,IPEX)
      if(iflag==1) then
        TMPP(IPS:IPS+IDUMX)=WK1_P(IPSX:IPSX+IDUMX,1)
      else
        TMPP(IPS:IPS+IDUMX)=WK1_P(IPS:IPS+IDUMX,1)
      endif
!
      IWK_P(:)=0
!      read(ifl) vname
      read(ifl) (IWK_P(IP),IP=IPSX,IPEX)
      if(iflag==1) then
        JOUT(IPS:IPS+IDUMX)=IWK_P(IPSX:IPSX+IDUMX)
      else
        JOUT(IPS:IPS+IDUMX)=IWK_P(IPS:IPS+IDUMX)
      endif
!
      if(IFLAG_COAL==icoal) then
!        read(ifl) (WK1_P(IP,VS_CHO),IP=IPS,IPE)
!        read(ifl) (WK1_P(IP,V_CHO),IP=IPS,IPE)
!        read(ifl) (WK1_P(IP,VS_H2O),IP=IPS,IPE)
!        read(ifl) (WK1_P(IP,V_H2O),IP=IPS,IPE)
!        read(ifl) (WK1_P(IP,VS_HCN),IP=IPS,IPE)
!        read(ifl) (WK1_P(IP,V_HCN),IP=IPS,IPE)
      endif
      enddo
!
      close(ifl)
      DEALLOCATE(WK1_P,IWK_P)
      DEALLOCATE(NPALLX,NINALLX,NSTAX,NRESX)
      return
  100 CALL FFRABORT(1,'Err: Opening Restart File:'//trim(fnam))
  200 CALL FFRABORT(1,'Err: Reading Restart File:BreakupModel')
      end subroutine read_particle_restart
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE LOCATE_PARTICLE_TO_CELL(imode,IP,INI,
     +            INSIDE,JOUT,VCTRP,XYZP,DIA,RHOP,HEIGHT,UVWP,
     +            HEIGHT0,CVCENT,SFAREA,SFCENT,NEIGHBOR,IE00,
     +            DNORM,IPCELL,BCKND,LVEDGE,WDRBND,PARCEL,DT0)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
!-------------------------
! --- [module arguments]
!-------------------------
      USE module_hpcutil 
      USE module_dimension 
      use module_io,only        : ifle,ifll
      use module_trace ,only    : CC,VELN
      use module_particle ,only : N_injctr,bc_nb
      use module_constant
      use module_boundary,only  : partkd,kpart_s,kpart_r,kpart_d,pkdwall
      use module_model,only  : ical_WDR
!
      implicit none
!----------------------
! --- DUMMY AUGUMENT
!----------------------
!      integer, intent(in)   :: NP0(N_injctr),NPALL(N_injctr),IE00
!      integer, intent(in)   :: NINALL(0:N_injctr)
      integer, intent(in)   :: IE00,IP,imode,INI
      real*8,  intent(in)   :: DT0
      integer, intent(in)   :: VCTRP      (MXCV_VP,0:MXBND_VP)
      integer, intent(in)   :: NEIGHBOR  (MEP  ,MXALLCV_P)
      integer, intent(in)   :: IPCELL    (      MXALLCV_P)
      real*8,  intent(inout):: XYZP      (MXPRT,3)
      real*8,  intent(in)   :: CVCENT    (3    ,MXALLCV)
      real*8,  intent(in)   :: HEIGHT0   (      MEP)
      real*8,  intent(in)   :: DNORM     (      MXCVFAC_P)
      real*8,  intent(in)   :: SFAREA    (4    ,MXCVFAC)
      real*8,  intent(in)   :: DIA       (      MXPRT)
      real*8,  intent(in)   :: RHOP      (      MXPRT)
      real*8,  intent(inout):: HEIGHT    (MEP  ,MXPRT)
      integer, intent(inout):: INSIDE    (      MXPRT)
      integer, intent(inout):: JOUT      (      MXPRT)
      REAL*8 ,INTENT(IN)    :: SFCENT   (    3,MXCVFAC)
      REAL*8,INTENT(INOUT)  :: UVWP     (      MXPRT,3)
      INTEGER,INTENT(in)    :: BCKND(MXCVFAC_P)
      INTEGER,INTENT(IN)    :: LVEDGE   (    2,MXCVFAC)
      real*8 ,intent(inout) :: WDRBND(        MXCV_WDR,N_inj)
      REAL*8,INTENT(IN)  :: PARCEL   (      MXPRT)
C------------------------
C      LOCAL ENTITIES 
C------------------------
      integer :: N, KE,K1,K2,MOUT,ICFL,ICFLL,IPCV,I,INB
      integer :: ICVB,nb,kdp,IDUM,ICVL
      real*8  :: XYZP0(3),UVWP0(3),TIMEMIN,VNORMAL
      REAL*8  :: dum1,dum2,dl,U1,U2,U3,R1,R2,R3
!--------------------------------------
! --- 
!--------------------------------------
!      do INI=1,N_injctr
!      IPS=NINALL(INI-1)+1+NP0(INI)
!      IPE=NINALL(INI-1)+NPALL(INI)
!      DO 50 IP=IPS,IPE
      IPCV=INSIDE(IP)
      INSIDE(IP)=IE00
      IF(IPCELL(IPCV)>=1) then
        JOUT(IP)=IPCELL(IPCV) 
        INSIDE(IP)=INSIDE(IP)
        return
      endif

      XYZP0(1)=CVCENT(1,IE00)
      XYZP0(2)=CVCENT(2,IE00)
      XYZP0(3)=CVCENT(3,IE00)
      UVWP0(1)=(XYZP(IP,1)-CVCENT(1,IE00))/DT0
      UVWP0(2)=(XYZP(IP,2)-CVCENT(2,IE00))/DT0
      UVWP0(3)=(XYZP(IP,3)-CVCENT(3,IE00))/DT0
!
      DO INB=1,VCTRP(IE00,0)
      HEIGHT(INB,IP)=HEIGHT0(INB)
      ENDDO
!
      IDUM=0
 100  continue
!
      IDUM=IDUM+1
      TIMEMIN=5.d0*DT0
      DO INB=1,VCTRP(IPCV,0)
      VNORMAL=0.d0
      ICFLL=VCTRP(IPCV,INB)
      ICFL=abs(ICFLL)
      DO I=1,3
      VNORMAL=VNORMAL-UVWP0(I)*SFAREA(I,ICFL)*SIGN(1,ICFLL)
      ENDDO
      IF(HEIGHT(INB,IP).GT.0.d0
     &  .AND.HEIGHT(INB,IP)/TIMEMIN.LT.VNORMAL) THEN
        TIMEMIN=HEIGHT(INB,IP)/VNORMAL
        MOUT=INB
      ENDIF
      VELN(INB)=VNORMAL
      ENDDO 
!---------------------------------------------------------------
! --- INTERPOLATE THE PHYSICAL VARIABLES IN THE PRESENT CELL 
!---------------------------------------------------------------
!      if(.false.) then
      IF(TIMEMIN.GT.DT0+SML) THEN
        DO INB=1,VCTRP(IPCV,0)
        HEIGHT(INB,IP)=HEIGHT(INB,IP)-VELN(INB)*DT0
        IF(HEIGHT(INB,IP).LT.0) THEN
          ICFL=ABS(VCTRP(IPCV,INB))
          nb=BCKND(ICFL)
          if(nb==0) return  !zhang2
          if(nb==0) then
            JOUT(IP)=A_dead         !MXALLCV+1
!            call FFRABORT(1,'ERR-1 A_dead') !zhang2
            WRITE(ifll,*) " OUT OF DOMAIN, IP= ",
     &                     my_rank,IP,TIMEMIN,DT0
            write(ifll,'(2X,a)') 
     &           'MSG: This particle is losed and reset in Inlet' 
            goto 150
          endif
          kdp=partkd(nb)            !BABA1
!          if(nb==bc_nb(INI)) goto 150
          if(kdp==kpart_d) then
            if(ical_WDR==1.and.pkdwall(nb).and.JOUT(IP)==0
     &        .and.nb/=bc_nb(INI)) then
              ICVL=LVEDGE(1,ICFL)
              dum1=PI*DIA(IP)**3/6.0D0*PARCEL(IP)/SFAREA(4,ICFL)
              WDRBND(ICVL,INI)=WDRBND(ICVL,INI)+dum1
              JOUT(IP)=BC_dead
            endif
            JOUT(IP)=BC_dead
            ICVB=abs(NEIGHBOR(INB,IPCV))
            INSIDE(IP)=ICVB
!
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)+U1*dum1
            XYZP(IP,2)=XYZP(IP,2)+U2*dum1
            XYZP(IP,3)=XYZP(IP,3)+U3*dum1
            UVWP(IP,:)=0.d0 
          elseif(kdp==kpart_s) then
            if(ical_WDR==1.and.pkdwall(nb).and.JOUT(IP)==0
     &        .and.nb/=bc_nb(INI)) then
              ICVL=LVEDGE(1,ICFL)
              dum1=PI*DIA(IP)**3/6.0D0*PARCEL(IP)/SFAREA(4,ICFL)
              WDRBND(ICVL,INI)=WDRBND(ICVL,INI)+dum1
              JOUT(IP)=BC_dead
            endif
!
            JOUT(IP)=-nb   !-ICFL    !-1 !HIBI
!            if(.not.pkdwall(nb).and.nb/=bc_nb(INI))  then
!              JOUT(IP)=BC_dead
!            endif
            ICVB=abs(NEIGHBOR(INB,IPCV))
            INSIDE(IP)=ICVB
!
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)+U1*dum1
            XYZP(IP,2)=XYZP(IP,2)+U2*dum1
            XYZP(IP,3)=XYZP(IP,3)+U3*dum1
            UVWP(IP,:)=0.d0 
          elseif(kdp==kpart_r) then
!            JOUT(IP)=-ICFL*INT(INT(IPCV/NCV))      !0
             
            if(.not.pkdwall(nb).and.nb/=bc_nb(INI))  then
              JOUT(IP)=BC_dead
            endif
            JOUT(IP)=0
            INSIDE(IP)=IPCV
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
!            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)-U1*dum1
            XYZP(IP,2)=XYZP(IP,2)-U2*dum1
            XYZP(IP,3)=XYZP(IP,3)-U3*dum1
            if(dum2<0.d0) then
              goto 150
            endif

            R1=dl*dum2*SFAREA(1,ICFL)
            R2=dl*dum2*SFAREA(2,ICFL)
            R3=dl*dum2*SFAREA(3,ICFL)
            U1=UVWP(IP,1)-R1
            U2=UVWP(IP,2)-R2
            U3=UVWP(IP,3)-R3
            UVWP(IP,1)=U1-R1
            UVWP(IP,2)=U2-R2
            UVWP(IP,3)=U3-R3
          endif
          GO TO 150 
        ENDIF 
        ENDDO 
        JOUT(IP)=-ICFL*INT(INT(IPCV/NCV))  !zhang
        INSIDE(IP)=IPCV !IPCV 
        GO TO 150
      ELSE
!--------------------------------------------------
C	  TRACE THE PARTICLE IN THE NEIGHBOR CELL
!--------------------------------------------------
        IF(NEIGHBOR(MOUT,IPCV)<=0) THEN 
          ICFL=ABS(VCTRP(IPCV,MOUT))
          nb=BCKND(ICFL)
          if(nb==0) CALL FFRABORT(1,'ERR:BCKND-3')
          kdp=partkd(nb)  !BABA1
!          if(nb==bc_nb(INI)) goto 150
          if(kdp==kpart_d) then

            if(ical_WDR==1.and.pkdwall(nb).and.JOUT(IP)==0
     &        .and.nb/=bc_nb(INI)) then
              ICVL=LVEDGE(1,ICFL)
              dum1=PI*DIA(IP)**3/6.0D0*PARCEL(IP)/SFAREA(4,ICFL)
              WDRBND(ICVL,INI)=WDRBND(ICVL,INI)+dum1
              JOUT(IP)=BC_dead
            else
              JOUT(IP)=BC_dead
            endif
            ICVB=abs(NEIGHBOR(INB,IPCV))
            INSIDE(IP)=ICVB
!
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)+U1*dum1
            XYZP(IP,2)=XYZP(IP,2)+U2*dum1
            XYZP(IP,3)=XYZP(IP,3)+U3*dum1
            UVWP(IP,:)=0.d0 
          elseif(kdp==kpart_s) then
            if(ical_WDR==1.and.pkdwall(nb).and.JOUT(IP)==0
     &        .and.nb/=bc_nb(INI)) then
              ICVL=LVEDGE(1,ICFL)
              dum1=PI*DIA(IP)**3/6.0D0*PARCEL(IP)/SFAREA(4,ICFL)
              WDRBND(ICVL,INI)=WDRBND(ICVL,INI)+dum1
              JOUT(IP)=BC_dead
            else
              JOUT(IP)=-nb   !-ICFL  !1 !HIBI
            endif
!            if(.not.pkdwall(nb).and.nb/=bc_nb(INI))  then
!               JOUT(IP)=BC_dead
 !           endif
            ICVB=abs(NEIGHBOR(MOUT,IPCV))
            INSIDE(IP)=ICVB
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)
            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
            dum1=dum1/(dum2+SML)
            XYZP(IP,1)=XYZP(IP,1)+U1*dum1
            XYZP(IP,2)=XYZP(IP,2)+U2*dum1
            XYZP(IP,3)=XYZP(IP,3)+U3*dum1
            UVWP(IP,:)=0.d0
          elseif(kdp==kpart_r) then
!            JOUT(IP)=-ICFL*INT(INT(IPCV/NCV))      !0
            if(.not.pkdwall(nb).and.nb/=bc_nb(INI))  then
              JOUT(IP)=BC_dead
            endif
            JOUT(IP)=0
            INSIDE(IP)=IPCV
            R1=SFCENT(1,ICFL)-XYZP(IP,1)
            R2=SFCENT(2,ICFL)-XYZP(IP,2)
            R3=SFCENT(3,ICFL)-XYZP(IP,3)
            dl=dsqrt(UVWP(IP,1)**2+UVWP(IP,2)**2+UVWP(IP,3)**2)+SML
            U1=UVWP(IP,1)/dl
            U2=UVWP(IP,2)/dl
            U3=UVWP(IP,3)/dl
            dum1=R1*SFAREA(1,ICFL)+R2*SFAREA(2,ICFL)+R3*SFAREA(3,ICFL)

            dum2=U1*SFAREA(1,ICFL)+U2*SFAREA(2,ICFL)+U3*SFAREA(3,ICFL)
!            dum1=dum1/(dum2+SML)
            
            XYZP(IP,1)=XYZP(IP,1)-U1*dum1
            XYZP(IP,2)=XYZP(IP,2)-U2*dum1
            XYZP(IP,3)=XYZP(IP,3)-U3*dum1
            if(dum2<0.d0) then
              goto 150
            endif
            R1=dl*dum2*SFAREA(1,ICFL)
            R2=dl*dum2*SFAREA(2,ICFL)
            R3=dl*dum2*SFAREA(3,ICFL)
            U1=UVWP(IP,1)-R1
            U2=UVWP(IP,2)-R2
            U3=UVWP(IP,3)-R3
            UVWP(IP,1)=U1-R1
            UVWP(IP,2)=U2-R2
            UVWP(IP,3)=U3-R3
          endif
          GO TO 150
        ENDIF
        IPCV=NEIGHBOR(MOUT,IPCV)
        IF(IPCELL(IPCV).GE.1) THEN
          INSIDE(IP)=IPCV
          JOUT(IP)=IPCELL(IPCV)
          GO TO 150 
        ENDIF 

        JOUT(IP)=-ICFL*INT(INT(IPCV/NCV))  !zhang

        INSIDE(IP)=IPCV  !zhang
        DO INB=1,VCTRP(IPCV,0)
        ICFLL=VCTRP(IPCV,INB)
        ICFL=abs(ICFLL)
        CC(INB)=-DNORM(ICFL)
        DO I=1,3
        CC(INB)=CC(INB)+SFAREA(I,ICFL)*XYZP0(I)
        ENDDO
        CC(INB)=CC(INB)*SIGN(1,ICFLL) !
        ENDDO
!
        DO INB=1,VCTRP(IPCV,0)
        HEIGHT(INB,IP)=CC(INB)
        ENDDO
        if(IDUM>=100) then !goto 150
           write(ifll,*) 'ERR: Can NOT FOUND NEXT CV FOR IP= ',IP 
           JOUT(IP)=A_dead !MXALLCV+1
!            call FFRABORT(1,'ERR-2 A_dead') !zhang2
!          call FFRABORT(1,'ERR: Contact your FFR supporter')
           goto 150
        endif
        GO TO 100  !zhang456
      ENDIF
!      endif
!
  150 CONTINUE
      if(IPCV==0) call 
     &   FFRABORT(1,'IPCV=000000000 LOCATE_PARTICLE_TO_CELL')
      IF(IPCELL(IPCV)>=1) then
        INSIDE(IP)=IPCV
        JOUT(IP)=IPCELL(IPCV) 
        return
      endif
  50  CONTINUE
!      enddo

      END SUBROUTINE LOCATE_PARTICLE_TO_CELL
!
!
C Pilch and Erdman breakup model. by Hisashi Ito.
CCSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BREAKUP_PILCH(INI,IT,ITS,
     1      DT0,GRDC,N,INSIDE,
     &      RHOP,DIA,XYZP,UVWP0,
     2      PARCEL,VEL,RHO,RMU,CVCENT,lBreakup,dTbu,dDiamReduc,
     3      CVVOLM,aks)
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! [DUMMY ARGUMENTS]
!
      use module_dimension
      use module_particle
      use module_io,only       : ifle,ifll
!
      IMPLICIT NONE

! Dummy Auguments
!
      INTEGER,INTENT(IN)  :: N,INI
      INTEGER,INTENT(IN)  :: IT,ITS
      REAL*8,INTENT(IN)   :: DT0
      REAL*8,INTENT(IN)   :: RHO(  MXALLCV,2)
      REAL*8,INTENT(IN)   :: RHOP(      MXPRT)
      REAL*8,INTENT(INOUT):: DIA(       MXPRT) 
      REAL*8,INTENT(IN)   :: UVWP0(3)
      REAL*8,INTENT(INOUT):: PARCEL(    MXPRT)
      REAL*8,INTENT(IN)   :: XYZP(      MXPRT,3)
      REAL*8,INTENT(IN)   :: RMU(MXALLCV)
      REAL*8,INTENT(IN)   :: VEL(MXALLCV,3,2)
      REAL*8,intent(in)   :: CVCENT(3,MXALLCV)
      REAL*8,intent(in)   :: GRDC(MXALLCV,3,3)
      INTEGER,intent(in)  :: INSIDE(MXPRT)
      LOGICAL,intent(out) :: lBreakup
      REAL*8, intent(out) :: dTbu
      REAL*8, intent(out) :: dDiamReduc
      REAL*8, INTENT(IN)  :: CVVOLM(MXALLCV)
      real*8 ,intent(in)  :: aks(MXALLCVR,MXRANS)
!
! Local Entities
!
      REAL*8 :: dWeber, dWec, Oh
      REAL*8 :: dTnondim, Vd, DIAs
      REAL*8 :: DIANEW, dRATE
      REAL*8 :: US,VS,WS,SLIPV
      REAL*8 :: UVW0(3),UVWFluc(3)
      INTEGER:: KE, i
      REAL*8, parameter :: PI=3.14159265358979323846264
      REAL*8, parameter :: B1=0.375d0
      REAL*8, parameter :: B2=0.2274d0
      LOGICAL,save      :: lCallFlag
      DATA lCallFlag/.FALSE./
!
!
      IF(.not.lCallFlag) THEN
        WRITE(ifll,*) 
     &   "MSG: [BREAKUP_PILCH] :PILCH BREAKUP MODEL is ON."
        lCallFlag = .TRUE.
      END IF

      KE=INSIDE(N)
      lBreakup=.TRUE.
      dDiamReduc = 1.d0
C
C Calculate Slip Velocity
C
        UVW0(1) = VEL(KE,1,2)
     *            +GRDC(KE,1,1)*(XYZP(1,N)-CVCENT(1,KE))
     *            +GRDC(KE,2,1)*(XYZP(2,N)-CVCENT(2,KE))
     *            +GRDC(KE,3,1)*(XYZP(3,N)-CVCENT(3,KE))
        UVW0(2) = VEL(KE,2,2)
     *            +GRDC(KE,1,2)*(XYZP(1,N)-CVCENT(1,KE))
     *            +GRDC(KE,2,2)*(XYZP(2,N)-CVCENT(2,KE))
     *            +GRDC(KE,3,2)*(XYZP(3,N)-CVCENT(3,KE))
        UVW0(3) = VEL(KE,3,2)
     *            +GRDC(KE,1,3)*(XYZP(1,N)-CVCENT(1,KE))
     *            +GRDC(KE,2,3)*(XYZP(2,N)-CVCENT(2,KE))
     *            +GRDC(KE,3,3)*(XYZP(3,N)-CVCENT(3,KE))
C
C Add Gaussian Noise to Velocity components (Random Walk Model)
C
        CALL add_random_walk(N,KE,GRDC,CVVOLM,UVWFluc,aks)
        UVW0(1)=UVW0(1)+UVWFluc(1)
        UVW0(2)=UVW0(2)+UVWFluc(2)
        UVW0(3)=UVW0(3)+UVWFluc(3)
C
C
         US = UVWP0(1) - UVW0(1)
         VS = UVWP0(2) - UVW0(2)
         WS = UVWP0(3) - UVW0(3)
C      
         SLIPV = SQRT(US*US + VS*VS + WS*WS)
C
C Calculate Reitz and Erdman Model Constants
C
      dWeber = RHO(KE,2)*(SLIPV*SLIPV)*DIA(N)/dLiquidStens(INI)
      Oh = dLiquidMu(INI)/((RHOP(N)*DIA(N)*dLiquidStens(INI))**0.5d0)
      dWec = 12.d0*(1.d0+1.077d0*(Oh**1.6d0))
C
C If  dWeber < dWec ,breakup is impossible
C
      IF(dWeber .LE. dWec) THEN
c        WRITE(ifll,*) 
!     &   "The Weber number is smaller than the critical Weber number."
c        WRITE(ifll,*) "So breakup is impossible."
        lBreakup = .FALSE.
        return
      END IF
C
C Calculate the dimensionless total breakup times dTnondim
C
      IF((dWeber .GE. 12.d0) .AND. (dWeber .LE. 18.d0)) THEN
        dTnondim = 6.d0*(dWeber-12.d0)**(-0.25d0)
      ELSE IF((dWeber .GT. 18.d0) .AND. (dWeber .LE. 45.d0)) THEN
        dTnondim = 2.45d0*(dWeber-12.d0)**(0.25d0)
      ELSE IF((dWeber .GT. 45.d0) .AND. (dWeber .LE. 350.d0)) THEN
        dTnondim = 14.1d0*(dWeber-12.d0)**(-0.25d0)
      ELSE IF((dWeber .GT. 350.d0) .AND. (dWeber .LE. 2670.d0)) THEN
        dTnondim = 0.766d0*(dWeber-12.d0)**(0.25d0)
      ELSE IF(dWeber .GT. 2670.d0) THEN
        dTnondim = 5.5d0
      ELSE
c        WRITE(ifll,*) "The Weber number is smaller than 12."
c        WRITE(ifll,*) "Breakup is impossible"
        lBreakup = .FALSE.
        RETURN
      END IF
C
C Calculate the total breakup times dT from dTnondim
C
      dTbu = dTnondim*DIA(N)*((RHOP(N)/RHO(KE,2))**0.5d0)/SLIPV
C
C Evaluate the stable droplet diameter DIAs
C
      Vd = SLIPV*((RHO(KE,2)/RHOP(N))**0.5d0)
     +          *(B1*dTnondim+B2*((dTnondim)**2))
      DIAs = dWec*dLiquidStens(INI)*((1-(Vd/SLIPV))**(-2.d0))
     +           /(RHO(KE,2)*(SLIPV**2.d0))
C
C Calculate the reduction rate of diameter dRATE
C
      dRATE = (DIAs-DIA(N))/dTbu
      IF (dRATE .GT. 0.d0) THEN
c        WRITE(ifll,*) "Error. Child Particle is larger than parent."
        lBreakup=.FALSE.
        return
      END IF
C
C Calcurate the New diamater
C

C Check DT_IN and dTbu
        IF (DT0 .GT. dTbu) THEN
c          WRITE(*,*) "Warning! DT_IN > dTbu"
          dDiamReduc = DIA(N)/DIAs
          return
        END IF
C
      DIANEW=dRATE*DT0+DIA(N)

C Check the value of DIANEW
        IF (DIANEW .LT. 0.d0) THEN
          WRITE(ifll,*) "Error. The new diameter(DIANEW) is negetive"
          stop
        END IF
C
C Calculate dDiamReduc
C
        dDiamReduc = DIA(N)/DIANEW
        return
      END SUBROUTINE BREAKUP_PILCH
C
!C********************************************************************
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C IMPROVED TAB breakup model. by Jun ARAI
C J.H.PARK ET AL. , Atomization and Sprays, Vol.12, pp.387-401, 2002
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BREAKUP_IMPROVED_TAB(IT,ITS,INI,
     1      DT0,GRDC,N,INSIDE,
     &      RHOP,DIA,XYZP,UVWP0,
     2      PARCEL,VEL,RHO,RMU,CVCENT,lBreakup,dTbu,dDiamReduc,
     3      CVVOLM,aks)
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! [DUMMY ARGUMENTS]
!
      USE module_particle  
      use module_dimension
      use module_io,only       : ifle,ifll
!
      IMPLICIT NONE

! Dummy Auguments
!
      INTEGER,INTENT(IN)  :: N,INI
      INTEGER,INTENT(IN)  :: IT,ITS
      REAL*8,INTENT(IN)   :: DT0
      REAL*8,INTENT(IN)   :: RHO (     MXALLCV,2)
      REAL*8,INTENT(IN)   :: RHOP(     MXPRT)
      REAL*8,INTENT(INOUT):: DIA (     MXPRT) 
      REAL*8,INTENT(INOUT):: UVWP0(3)
      REAL*8,INTENT(INOUT):: PARCEL(   MXPRT)
      REAL*8,INTENT(IN)   :: XYZP(   3,MXPRT)
      REAL*8,INTENT(IN)   :: RMU(      MXALLCV)
      REAL*8,INTENT(IN)   :: VEL(      MXALLCV,3,2)
      REAL*8,intent(in)   :: CVCENT(3, MXALLCV)
      REAL*8,intent(in)   :: GRDC(     MXALLCV,3,3)
      INTEGER,intent(in)  :: INSIDE(   MXPRT)
      LOGICAL,intent(out) :: lBreakup
      REAL*8, intent(out) :: dTbu
      REAL*8, intent(out) :: dDiamReduc
      REAL*8, INTENT(IN)   :: CVVOLM(MXALLCV)
      real*8 ,intent(in)  :: aks(MXALLCVR,MXRANS)
!
! Local Entities
!
      REAL*8 :: dWeber,dOmega,dVolume,dYC(5),dVC(5),dRhor,dMur
      REAL*8 :: dKY(4),dKV(4),radi, dRep, dTRK(5)
      REAL*8 :: dA1,dA2,dA3,dA4,dTemp1,dTemp2
      REAL*8 :: US,VS,WS,SLIPV,dDIANEW
      REAL*8 :: UVW0(3),UVWFluc(3),Vnorm
      INTEGER:: KE, IRK, ITERIN, I
      REAL*8, parameter :: PI=3.141592653589793238d0
      LOGICAL,save      :: lCallFlag=.FALSE.
!
!
      if(.not.lCallFlag) then
        write(ifll,*)
     +    "MSG:[BREAKUP_IMPROVED_TAB]: ITAB BREAKUP MODEL is ON."
        lCallFlag = .TRUE.
      endif
C
      KE=INSIDE(N)
      lBreakup=.FALSE.
      dDiamReduc = 1.d0
      radi = DIA(N)*0.5d0
      dTbu = 0.d0
C
C Calculate Slip Velocity
C
        UVW0(1) = VEL(KE,1,2)
     *            +GRDC(KE,1,1)*(XYZP(1,N)-CVCENT(1,KE))
     *            +GRDC(KE,2,1)*(XYZP(2,N)-CVCENT(2,KE))
     *            +GRDC(KE,3,1)*(XYZP(3,N)-CVCENT(3,KE))
        UVW0(2) = VEL(KE,2,2)
     *            +GRDC(KE,1,2)*(XYZP(1,N)-CVCENT(1,KE))
     *            +GRDC(KE,2,2)*(XYZP(2,N)-CVCENT(2,KE))
     *            +GRDC(KE,3,2)*(XYZP(3,N)-CVCENT(3,KE))
        UVW0(3) = VEL(KE,3,2)
     *            +GRDC(KE,1,3)*(XYZP(1,N)-CVCENT(1,KE))
     *            +GRDC(KE,2,3)*(XYZP(2,N)-CVCENT(2,KE))
     *            +GRDC(KE,3,3)*(XYZP(3,N)-CVCENT(3,KE))
C
C Add Gaussian Noise to Velocity components (Random Walk Model)
C
        CALL add_random_walk(N,KE,GRDC,CVVOLM,UVWFluc,aks)
        UVW0(1)=UVW0(1)+UVWFluc(1)
        UVW0(2)=UVW0(2)+UVWFluc(2)
        UVW0(3)=UVW0(3)+UVWFluc(3)
C
C
         US = UVWP0(1) - UVW0(1)
         VS = UVWP0(2) - UVW0(2)
         WS = UVWP0(3) - UVW0(3)
C
         SLIPV = SQRT(US*US + VS*VS + WS*WS)
         if(SLIPV .le. 0.d0) return
C
C Calculate TAB Model Constants
C
      DTRK(1)= DT0*SLIPV/radi
      ITERIN = INT(DTRK(1))
C
C If ITERIN is too large, the deformation is negligible.
C
      if(ITERIN .gt. 50) return
C
      DTRK(1) = DTRK(1)/DBLE(ITERIN)
      DTRK(2)= DTRK(1)*0.5d0
      DTRK(3)= DTRK(1)*0.5d0
      DTRK(4)= DTRK(1)
      DTRK(5)= DTRK(1)/6.d0
      dYC(1) = dPDstion(N)
      dVC(1) = dPDstvel(N)
      dWeber = RHO(KE,2)*(SLIPV*SLIPV)*radi/dLiquidStens(INI)
C
C if the weber number is too small, deformation is negligible.
C
      if(dWeber .lt. 0.5d0) return
C
      dOmega = dCk(INI)*dLiquidStens(INI)/(RHOP(N)*radi*radi*radi)
      dRep = radi*SLIPV*RHO(KE,2)/RMU(KE)
      dRhor = RHOP(N)/RHO(KE,2)
      dMur  = dLiquidMu(INI)/RMU(KE)
C
      dA1=-1.d0*dCd(INI)*dMur/(dRep*dRhor)
      dA2=(2.d0*dCf(INI)-dCk(INI)/dWeber)/dRhor
      dA3=dCf(INI)*dCb(INI)/dRhor
      dA4=dCf(INI)/(dCb(INI)*dRhor)
C
      do I=1, ITERIN
        do IRK=2, 5
          dKV(IRK-1)
     +      = dA1*dVC(IRK-1) + (dA2+dA3*dYC(IRK-1))*dYC(IRK-1) + dA4
          dVC(IRK) = dVC(1) + DTRK(IRK)*dKV(IRK-1)
          dKY(IRK-1) = dVC(IRK)
          dYC(IRK) = dYC(1) + DTRK(IRK)*dKY(IRK-1)
        enddo
C
        dPDstvel(N) = dVC(1)
     +     + DTRK(5)*(dKV(1) + 2.d0*dKV(2) + 2.d0*dKV(3) + dKV(4))
        dPDstion(N) = dYC(1)
     +     + DTRK(5)*(dKY(1) + 2.d0*dKY(2) + 2.d0*dKY(3) + dKY(4))
C------------------------
C Breakup Criterion
C------------------------
        dYC(5)=1.d0+dCb(INI)*dPDstion(N)
        if((2.d0*dYC(5)**5 + 1.d0/dYC(5) - 4.d0/dYC(5)**4)
     +     .gt. 2.5107d0*dWeber) THEN
C------------------------------------------------
C Calculate New Diameter by Energy Conservation
C------------------------------------------------
          dDiamReduc = 1.0+8.0*dK(INI)/20.0
     1      +RHOP(N)*DIA(N)*DIA(N)*DIA(N)*dPDstvel(N)*dPDstvel(N)
     2       /(8.0*dLiquidStens(INI))*(6.0*dK(INI)-5.0)/120.0
C------------------------------------------------
C +++++++++++++++ Add Normal Velocity +++++++++++
C------------------------------------------------
          call utl_random(dTemp1)
          call utl_random(dTemp2)
          dTemp1 = dTemp1*2.d0*PI
          dTemp2 = dTemp2*2.d0*PI
C
          Vnorm =  dCv(INI)*dCb(INI)*DIA(N)*0.5d0*dPDstvel(N)
          UVWP0(1)=UVWP0(1)+Vnorm*DSIN(dTemp1)*DCOS(dTemp2)
          UVWP0(2)=UVWP0(2)+Vnorm*DCOS(dTemp1)*DCOS(dTemp2)
          UVWP0(3)=UVWP0(3)+Vnorm*DSIN(dTemp1)
C
C++++++++++++++++++++++++++++++++++++++++++++++++++
C
          if(dDiamReduc.lt.1.0) then
            write(ifll,*) "Err. Child Particle is larger than parent"
          endif
          if(dDiamReduc.lt.0.0) then
            write(ifll,*) "Err. Negative Diameter"
            CALL FFRABORT(1,'Err:Normal TAB')
          endif
          lBreakup=.TRUE.
          goto 100
        END IF
      ENDDO
 100  CONTINUE
C      if(IT .gt. 300) WRITE(ifll,*) "DEBUG02"

      return
C
      END SUBROUTINE BREAKUP_IMPROVED_TAB
C
C
C Writing Restart file for Particle Calculation, Jun ARAI,2005-10-11
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE write_particle_restart(
     &  IT,time,delt,ffexit,outputfile,MP,MEP,INSIDE,
     &  PMASS,DIA,XYZP,UVWP,PARCEL,RHOP,TMPP,
     &  HEIGHT,FU,FV,FW,PYS,HIS_P,JOUT,
     &  ierror)
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension,only : MXCOMP,NCOMP,MXPRT
      use module_hpcutil,  only : NPE,my_rank,RESTRTout_P
      USE module_particle  
      use module_output,only    : outchk,fnmlst,none,yes
      use module_io,only        : ifle,ifll,getfil
      USE module_particle,only  : IFLAG_COAL,icoal,iglass,
     &                            VS_CHO,V_CHO,V_H2O,VS_H2O,
     &                            VS_HCN,V_HCN
!      use module_trace ,only    : WK1_P
!
!  1. Write restart file
!
      implicit none
!
!--- [dummy arguments]
!
      integer,intent(in) :: IT
      logical,intent(in) :: ffexit,outputfile
      real*8 ,intent(in) :: time,delt
      INTEGER,INTENT(IN) :: MP,MEP
!
      INTEGER,INTENT(IN) :: INSIDE(        MP)
      INTEGER,INTENT(IN) :: JOUT  (        MP)

      REAL*8 ,INTENT(IN) :: PMASS(         MP,2)
      REAL*8 ,INTENT(IN) :: PARCEL(        MP)
      REAL*8 ,INTENT(IN) :: DIA(           MP)
      REAL*8 ,INTENT(IN) :: RHOP(          MP)
      REAL*8 ,INTENT(IN) :: XYZP(          MP,3)
      REAL*8 ,INTENT(IN) :: UVWP(          MP,3)
      REAL*8 ,INTENT(IN) :: HEIGHT(        MEP,MP)
      REAL*8 ,INTENT(IN) :: FU(            MP,2)
      REAL*8 ,INTENT(IN) :: FV(            MP,2)
      REAL*8 ,INTENT(IN) :: FW(            MP,2)
      REAL*8 ,INTENT(IN) :: PYS(  MXPRT,MXCOMP,2)
      REAL*8,INTENT(IN)  :: TMPP(       MXPRT)
      REAL*8,INTENT(IN)  :: HIS_P(       MXPRT)
!      REAL*8,INTENT(IN) :: PHHH(  MXPRT,2)
!      
      integer,intent(OUT):: ierror
!
!--- [local entities]
!
      integer,parameter  :: vnum=14
      character(LEN=128) :: fnam
      character(len=8)   :: text
      character(LEN=8)   :: vname(vnum)
      integer :: I,J,ifl,IP,ios,out,INI,IPS,IPE
!
!
      ierror=0
      vname = ' '
      vname(1 ) = 'XYZ_P'
      vname(2 ) = 'UVW_P'
      vname(3 ) = 'DIA_P'
      vname(4 ) = 'PARCEL_P'
      vname(5 ) = 'RHO_P'
      vname(6 ) = 'P_MASS'
      vname(7 ) = 'INSIDE'
      vname(8 ) = 'HEIGHT'
      vname(9 ) = 'FU'
      vname(10) = 'FV'
      vname(11) = 'FW'
      vname(12) = 'DPDSTION'
      vname(13) = 'DPDSTVEL'
      vname(14) = 'HISTORY'
!
      call outchk(fnmlst(9),IT,time,delt,out)
      if(out.eq.none) return
      if(.not.ffexit .and. out.ne.yes.and..not.outputfile) return
!
      call getfil(ifl,fnam,'p_restart')
      if(NPE>1) then
        fnam=RESTRTout_P
      else
      endif
      fnam = trim(trim(fnam))
      write(text,'(i8)') IT
      text=adjustl(text)
      fnam=TRIM(TRIM(fnam)//'_'//TRIM(text))
      
      open(ifl,file=trim(fnam),form='unformatted',
     &     status='unknown',iostat=ios)
      do INI=1,N_injctr
      write(ifl) IT,time,INI,NPALL(INI),NINALL(INI),
     &          NSTA(INI),NRES(INI)
      enddo
      do INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
!
      write(ifl) vname(1)
      write(ifl) ((XYZP(IP,J),IP=IPS,IPE),J=1,3)
      write(ifl) vname(2)
      write(ifl) ((UVWP(IP,J),IP=IPS,IPE),J=1,3)
      write(ifl) vname(3)
      write(ifl) (DIA(IP),IP=IPS,IPE)
      write(ifl) vname(4)
      write(ifl) (PARCEL(IP),IP=IPS,IPE)
      write(ifl) vname(5)
      write(ifl) (RHOP(IP),IP=IPS,IPE)
      write(ifl) vname(6)
      write(ifl) (PMASS(IP,1),IP=IPS,IPE)
      write(ifl) (PMASS(IP,2),IP=IPS,IPE)
      write(ifl) vname(14)
      write(ifl) (HIS_P(IP),IP=IPS,IPE)
      write(ifl) vname(7)
      write(ifl) (INSIDE(IP),IP=IPS,IPE)
      write(ifl) vname(8)
      write(ifl) ((HEIGHT(J,IP),J=1,MEP),IP=IPS,IPE)
      write(ifl) vname(9)
      write(ifl) ((FU(IP,J),IP=IPS,IPE),J=1,2)
      write(ifl) vname(10)
      write(ifl) ((FV(IP,J),IP=IPS,IPE),J=1,2)
      write(ifl) vname(11)
      write(ifl) ((FW(IP,J),IP=IPS,IPE),J=1,2)
      if(.not.lBreakupNONE(INI)) then
        write(ifl) vname(12)
        write(ifl) (dPDstion(IP),IP=IPS,IPE)
        write(ifl) vname(13)
        write(ifl) (dPDstvel(IP),IP=IPS,IPE)
      endif
      write(ifl) ((PYS(IP,J,1),IP=IPS,IPE),J=1,NCOMP)
      write(ifl) (TMPP(IP),IP=IPS,IPE)
      write(ifl) (JOUT(IP),IP=IPS,IPE)
      if(IFLAG_COAL==icoal) then
!        write(ifl) (WK1_P(IP,VS_CHO),IP=IPS,IPE)
!        write(ifl) (WK1_P(IP,V_CHO),IP=IPS,IPE)
!        write(ifl) (WK1_P(IP,VS_H2O),IP=IPS,IPE)
!        write(ifl) (WK1_P(IP,V_H2O),IP=IPS,IPE)
!        write(ifl) (WK1_P(IP,VS_HCN),IP=IPS,IPE)
!        write(ifl) (WK1_P(IP,V_HCN),IP=IPS,IPE)
      endif
      enddo
      close(ifl)
      return
      end subroutine write_particle_restart
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Latent_P(
     &  IT,ITER_P,time,DT_IN,N_injctr,
     &  NINALL,INSIDE,JOUT,PARCEL,TMPP,NPALL,    !WK1_P,
     &  PMASS,PYS,
     &  DYDTP,DHDTP,EVAP_Y,EVAP_H,CVVOLM,MASSYS,
     &  ierror)
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension 
      USE module_hpcutil     
      USE module_constant,only : SML
      use module_io,only       : ifle,ifll
      USE module_particle,only : Evap_CHO,Avap_CHO,ical_la_coal
      use module_particle,only : IFLAG_COAL,VS_CHO,V_CHO,V_H2O,VS_H2O,
     &                           VS_HCN,V_HCN
      use module_particle,only : II_CHO,II_H2O,II_C,II_ash
      use module_particle,only : ifuel,icoal,iliquid,isolid,
     &                           particle_comp
      use module_particle,only : Avap_H2O,Evap_H2O,LH_H2O,Avap_HCN,
     &                           Evap_HCN,iglass
      use module_particle,only : II_HCN,LH_CHO,latent
      use module_species,only  : gascns
!
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)    :: IT,ITER_P,N_injctr
      REAL*8,INTENT(IN)     :: time,DT_IN
      INTEGER,INTENT(IN)    :: INSIDE   (      MXPRT)
      INTEGER,INTENT(IN)    :: NINALL   (    0:N_injctr)
      INTEGER,INTENT(IN)    :: NPALL    (      N_injctr)
      INTEGER,INTENT(INout) :: JOUT     (      MXPRT)
      REAL*8,INTENT(IN)     :: PARCEL   (      MXPRT)
      REAL*8,INTENT(IN)     :: TMPP     (      MXPRT)
      REAL*8,INTENT(IN)     :: CVVOLM   (      MXALLCV)
      REAL*8,INTENT(INOUT)  :: PYS      (      MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: PMASS    (      MXPRT)
!
      REAL*8,INTENT(INOUT)  :: DYDTP    (MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: DHDTP    (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: EVAP_Y   (      MXALLCV_P,MXCOMP)
      REAL*8,INTENT(INOUT)  :: EVAP_H   (      MXALLCV_P,2)
      REAL*8,INTENT(INOUT)  :: MASSYS   (      MXALLCV,MXcomp)

      INTEGER,INTENT(INOUT) :: ierror
!
! --- [local entities]
!
      integer :: IP,ICFL,INI,IPS,IPE
      INTEGER :: IPCV
      REAL*8  :: dum1,dum2,dum3,pdum,DT0,Kv,vol
      REAL*8,parameter :: f_latent_pt=1.d0 
      REAL*8,save  :: f_latent_gas=1.d0-f_latent_pt
!      
      return
      end subroutine Latent_P
!
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE UPDATE_P(
     &   IT,ITER_P,time,DT_IN,DIMV,N_injctr,PI,
     &   NINALL,JOUT,NPALL,INSIDE,DYDTP,PARCEL,
     &   PMASS,DIA,RHOP,HIS_P,XYZP,PRT_VF)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension   
      USE module_hpcutil    
      
      USE module_particle,only : ifuel,icoal,iglass,
     &                           itype,idensity,ivolume,mass_p_total
      use module_model,  only  : encomp,ical_t,ical_reac
      use module_chemcntl,only : igT_iter
      USE module_particle,only : injkd,i_coal,inlet
      USE module_constant,only : SML
      
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)    :: IT,ITER_P,DIMV,N_injctr
      REAL*8,INTENT(IN)     :: time,DT_IN,PI
      INTEGER,INTENT(IN)    :: INSIDE   (      MXPRT)
      INTEGER,INTENT(IN)    :: NINALL   ( 0:N_injctr)
      INTEGER,INTENT(IN)    :: NPALL    (   N_injctr)
      INTEGER,INTENT(INOUT) :: JOUT     (      MXPRT)
      REAL*8,INTENT(INOUT)  :: PARCEL   (      MXPRT)
      REAL*8,INTENT(INOUT)  :: DIA      (      MXPRT)
      REAL*8,INTENT(INOUT)  :: PMASS    (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: RHOP     (      MXPRT)
      REAL*8,INTENT(IN)     :: DYDTP    (MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: HIS_P    (      MXPRT)
      REAL*8,INTENT(IN)     :: XYZP     (      MXPRT,3)
      REAL*8,INTENT(INOUT)  :: PRT_VF   (      MXALLCV_P)
!
! --- [local entities]
!
      integer :: IP,ICFL,INI,IPS,IPE,nbx4=4
      INTEGER :: IPCV,icom,idum,nbx=8,nbx1=3,nbx2=4,idum2=0
      REAL*8  :: DT0,dum1,dum2,dum3
      REAL*8  :: dmin,dmin1,dmax,dmax1
      REAL*8,parameter ::  third=1.d0/3.d0
      REAL*8,save :: dum_eva=0.d0
!
      do 6000 INI=1,N_injctr

!        if(injkd(INI)/=i_coal.and.injkd(INI)/=inlet.and.) cycle
!        if(ifuel(INI)==iglass) cycle
!        IF(ifuel(INI)/=icoal) cycle
         dum_eva=0.d0
         dmin1=1.d10
         dmax1=-1.d10
         idum=0
         idum2=0

          DT0=DT_IN
          IPS=NINALL(INI-1)+1
          IPE=NINALL(INI-1)+NPALL(INI)
          if(itype(INI)==idensity) then
            DO 6100 IP=IPS,IPE
            if(JOUT(IP)==HPC_dead) cycle
            IPCV=INSIDE(IP)
            dum1=0.d0
            dum1=PMASS(IP,1)/RHOP(IP)
            DIA(IP)=(6.d0*dum1/PI)**third 
            if(dum1<SML) JOUT(IP)=MASS_dead
!
!            if(JOUT(IP)==0) then   !  hibi-san    
!               idum2=idum2+1
!               HIS_P(IP)=HIS_P(IP)+DT_IN !HIBI
!            endif
            HIS_P(IP)=HIS_P(IP)+DT_IN   !hibi-san
            PRT_VF(IPCV)=PRT_VF(IPCV)+dum1*PARCEL(IP)*DT_IN 
!------------------------------------
! hibi-san
!------------------------------------
!!            if(-JOUT(IP)==nbx) then !HIBI  case-2(ann-1)
!!            if(-JOUT(IP)==nbx1.or.-JOUT(IP)==nbx2) then !HIBI  case-1
!!            if(XYZP(IP,1)>3.d0) then !HIBI  case-2(ann-2)
!            if(-JOUT(IP)==nbx4) then !HIBI  case-2(jiko-ji)
!!              JOUT(IP)=A_dead  !HIBI  case-2(ann-2)
!              dum_eva=dum_eva+HIS_P(IP)
!              idum=idum+1
!              if(HIS_P(IP)<dmin1) then
!                dmin1=HIS_P(IP)!dum_eva
!              endif
!              if(HIS_P(IP)>dmax1) then
!                dmax1=HIS_P(IP)!dum_eva
!              endif
!            endif

 6100       ENDDO


            

!            print*,'INJ live particle= ',idum2
!            print*,'nb= ',nbx
!            print*,'INJ= ',INI,'arre  nu= ',idum
!            print*,'INJ= ',INI,'ave time= ',dum_eva/dble(idum+1.d-10)
!            print*,'INJ= ',INI,'min time= ',dmin1
!            print*,'INJ= ',INI,'max time= ',dmax1

          elseif(itype(INI)==ivolume) then
            DO 6200 IP=IPS,IPE
            if(JOUT(IP)==HPC_dead) cycle
            IPCV=INSIDE(IP)
            dum1=DIA(IP)**3*PI/6.d0
            RHOP(IP)=PMASS(IP,1)/(dum1+SML)
            if(PMASS(IP,1)<SML) JOUT(IP)=MASS_dead
!            if(JOUT(IP)==0) HIS_P(IP)=HIS_P(IP)+DT_IN !HIBI
            HIS_P(IP)=HIS_P(IP)+DT_IN
            PRT_VF(IPCV)=PRT_VF(IPCV)+dum1*PARCEL(IP)*DT_IN 
 6200       ENDDO
          endif

 6000   ENDDO

      return
      end subroutine UPDATE_P 
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE P_result_S(
     &   IT,ITER_P,time,dTSave,DIMV,N_injctr,ffexit,outputfile,
     &   MAT_CAL,MAT_NO,MAT_INDEX,MAT_CVEXT,
     &   NINALL,NPALL,HIS_P,JOUT,
     &   XYZP,UVWP,DIA,PARCEL,PMASS,TMPP,PYS)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension   
      USE module_hpcutil    
      use module_output,only   : outchk,fnmlst,none,yes,mlt_rslt 
      use module_species,only  : spcnam 
      use module_io,only       : getfil,ifle,ifll,gdScale 
      USE module_particle,only : dYP0,dZP0,dXP0,P_MF,starttime 
      use module_particle,only : iPartLU,IPSND
      use module_particle,only : PMASSC
!
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!     
      INTEGER,INTENT(IN)    :: IT,ITER_P,DIMV,N_injctr
      REAL*8,INTENT(IN)     :: time,dTSave
      logical,intent(in)    :: ffexit,outputfile
      INTEGER,INTENT(IN)    :: MAT_NO   (      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_INDEX(      0:MXMAT)
      LOGICAL,INTENT(IN)    :: MAT_CAL  (      0:MXMAT)
      INTEGER,INTENT(IN)    :: NINALL   ( 0:N_injctr)
      INTEGER,INTENT(IN)    :: NPALL    (   N_injctr)
      INTEGER,INTENT(IN)    :: JOUT     (      MXPRT)
!
      REAL*8,INTENT(INOUT)     :: XYZP     (      MXPRT,3)
      REAL*8,INTENT(IN)     :: UVWP     (      MXPRT,3)
      REAL*8,INTENT(INOUT)  :: DIA      (      MXPRT)
      REAL*8,INTENT(IN)     :: PARCEL   (      MXPRT)
      REAL*8,INTENT(IN)     :: PMASS    (      MXPRT,2)
      REAL*8,INTENT(IN)     :: HIS_P    (      MXPRT)
      REAL*8,INTENT(IN)     :: TMPP     (      MXPRT,2)
      REAL*8,INTENT(IN)     :: PYS      (      MXPRT,MXCOMP,2)
!
      REAL*8 :: dum1,dum2,dum3,dumm
      integer:: K,out,N0,IPS,IPE,IP,INI,icom,idum
      CHARACTER*6  :: CHA2,CHA3,for1
      CHARACTER*15 :: OUTNAME 
      CHARACTER(LEN=128) :: pt_result=''
      character*128:: formt
      integer :: ifl=0
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
      call outchk(fnmlst(10),IT,time,dTSave,out)
      if(ffexit.or.out==yes.or.outputfile) then
        if(ifl==0) then
          call getfil(ifl,pt_result,'p_result')
        endif
!
        if(mlt_rslt==yes) then
          write(CHA2,fmt="(I6.6)") IT
          OUTNAME=trim(pt_result)//CHA2
        else
          OUTNAME=trim(pt_result)
        endif
        N0=3+3+5+ncomp+1
!
        if(NPE>1) THEN
C----------------------
C MPI Merged Output 
C----------------------
          CALL MERGEOUT_PTRESULT_MPI(
     %  ncomp,N0,OUTNAME,MXPRT,NPALL,time,N_injctr,NINALL,
     &  XYZP,UVWP,DIA,PARCEL,HIS_P,PMASS(:,1),PMASSC,
     &  TMPP(:,1),PYS(:,:,1),JOUT,IPSND,
     &  iPartLU,dYP0,dZP0,dXP0,P_MF,starttime
     &   )
        else
C----------------------
C Single Output 
C----------------------
          write(for1,'(I3)') ncomp
          write(formt,'(a,a,a)')
     &   '(I10,6(3X,E15.8),5(3X,E15.8),',trim(for1),'(3X,E15.8))'
!
          OPEN(ifl,FILE=OUTNAME,STATUS="REPLACE")
          write(ifl,fmt="(2I6,' TIME=',E12.7,
     +        ' X Y Z U V W DIA HIS PARCEL,PMASS,TMP_P,Ys_P')")
     &    N_injctr,N0,time
          WRITE(ifl,'(I3)') ncomp
          do icom=1,ncomp
          WRITE(ifl,'(a)') trim(spcnam(icom)) 
          enddo
          do INI=1,N_injctr
          dum1=dXP0(INI)/gdScale
          dum2=dYP0(INI)/gdScale
          dum3=dZP0(INI)/gdScale
          IPS=NINALL(INI-1)+1
          IPE=NINALL(INI-1)+NPALL(INI)
          idum=NPALL(INI)
          if(IPS>IPE) idum=0
          write(ifl,2501) idum,dum1,dum2,dum3,starttime(INI)
          if(P_MF(INI)==1) then
            DO IP=IPS,IPE 
            if(JOUT(IP)==mass_dead) then
              XYZP(IP,1)=dXP0(INI)
              XYZP(IP,2)=dYP0(INI)
              XYZP(IP,3)=dZP0(INI)
            endif
            WRITE(ifl,formt) IP,
     &        (XYZP(IP,K)/gdScale,K=1,3),
     &        (UVWP(IP,K),K=1,3),
     &         DIA(IP),
     &         HIS_P(IP),  !PARCEL(IP),HIS_P(IP)
     &         PARCEL(IP),
     &         PMASS(IP,1),
     &         TMPP(IP,1),
     &        (PYS(IP,K,1),K=1,ncomp)
            ENDDO
          else
            call FFRABORT(1,'ERR: Particle result file option error')
            DO IP=IPS,IPE
            if(JOUT(IP)==mass_dead) then
              XYZP(IP,1)=dXP0(INI)
              XYZP(IP,2)=dYP0(INI)
              XYZP(IP,3)=dZP0(INI)
            endif
            WRITE(ifl,formt) IP,
     &        (XYZP(IP,K)/gdScale,K=1,3),
     &        (UVWP(IP,K),K=1,3),
     &         DIA(IP),    !7
     &         HIS_P(IP),  !8
     &         PARCEL(IP), !9
     &         PMASS(IP,1),!10
     &         TMPP(IP,1), !11
     &        (PMASS(IP,1)*PYS(IP,K,1),K=1,ncomp)
            ENDDO
          endif
          enddo
! 2500     format(I10,6(3X,E15.8),4(3X,E15.8),16(3X,E15.8)) 
! 2500     format(I10,6(3X,E15.8),4(3X,E15.8),9(3X,E15.8)) 
 2501     FORMAT(I10,4E15.8) 
          CLOSE(ifl)
        ENDIF
      ENDIF
!
      return
      end SUBROUTINE P_result_S
!
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Evaporation_P(
     &  IT,ITER_P,time,DT_IN,N_injctr,
     &  NINALL,INSIDE,JOUT,PARCEL,TMPP,NPALL,
     &  PMASS,PYS,
     &  DYDTP,DHDTP,EVAP_Y,EVAP_H,CVVOLM,MASSYS,
     &  yys,pp0,tmp,rho,rhop,rmu,vel,uvwp,rds,DIA,
     &  rmd,
     &  ierror)
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE module_dimension 
      USE module_hpcutil     
      use module_io,only       : ifle,ifll
      USE module_particle,only : evakd,ievaF,iAbSi,iLaKn,ievaU,
     &                           ical_evap,dLiquidMu
      use module_particle,only : ifuel,icoal,iliquid,isolid,
     &                           particle_comp,
     &                           iglass
!
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)    :: IT,ITER_P,N_injctr
      REAL*8,INTENT(IN)     :: time,DT_IN
      INTEGER,INTENT(IN)    :: INSIDE   (      MXPRT)
      INTEGER,INTENT(IN)    :: NINALL   (    0:N_injctr)
      INTEGER,INTENT(IN)    :: NPALL    (      N_injctr)
      INTEGER,INTENT(INout) :: JOUT     (      MXPRT)
      REAL*8,INTENT(IN)     :: PARCEL   (      MXPRT)
      REAL*8,INTENT(IN)     :: TMPP     (      MXPRT)
      REAL*8,INTENT(IN)     :: CVVOLM   (      MXALLCV)
      REAL*8,INTENT(INOUT)  :: PYS      (      MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: PMASS    (      MXPRT,2)
!
      REAL*8,INTENT(INOUT)  :: DYDTP    (MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: DHDTP    (      MXPRT,2)
      REAL*8,INTENT(INOUT)  :: EVAP_Y   (      MXALLCV_P,MXCOMP)
      REAL*8,INTENT(INOUT)  :: EVAP_H   (      MXALLCV_P,2)
      REAL*8,INTENT(INOUT)  :: MASSYS   (      MXALLCV,MXcomp)

      INTEGER,INTENT(INOUT) :: ierror
!
      REAL*8,INTENT(IN)     :: RHO      (      MXALLCV)
      REAL*8,INTENT(IN)     :: RMU      (      MXALLCV)
      REAL*8,INTENT(IN)     :: VEL      (      MXALLCV,3)
      REAL*8,INTENT(IN)     :: YYS      (      MXALLCV,MXCOMP)
      REAL*8,INTENT(IN)     :: TMP      (      MXALLCV)
      REAL*8,INTENT(IN)     :: RHOP     (      MXPRT)
      REAL*8,INTENT(IN)     :: pp0      (      MXALLCV)
      REAL*8,INTENT(IN)     :: UVWP     (      MXPRT,3)
      REAL*8,INTENT(IN)     :: RDS      (      MXALLCV,MXcomp)
      REAL*8,INTENT(IN)     :: RMD      (      MXALLCV)
      REAL*8, intent(in)    :: DIA      (MXPRT)
!
! --- [local entities]
!
      integer :: IP,ICFL,INI,IPS,IPE,ICOM
!      
      ierror=0
!
      return
      end subroutine Evaporation_P
!
!
!
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine drop_distrib(INI,DIA,dSMD,random)
!CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only       : ifle,ifll
      USE module_constant,only : SML
!      USE module_particle
      USE module_particle,only : dropkd,iconst,iRosR,inoml,inuki,ikuro
!
      IMPLICIT NONE
!
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)     :: INI
      REAL*8 ,INTENT(IN)     :: dSMD,random
      REAL*8 ,INTENT(OUT)    :: DIA
!
! --- [local entities]
!
      if(dropkd(INI)==iRosR) then
        call FFRABORT(1,'ERR: normal distribution NOT support')
      elseif(dropkd(INI)==inoml) then
        call FFRABORT(1,'ERR: normal distribution NOT support')
      elseif(dropkd(INI)==ikuro) then
        call FFRABORT(1,'ERR: normal distribution NOT support')
      elseif(dropkd(INI)==inuki) then
        CALL NUKIYAMA_TANAZAWA_VOLUMEBASE(INI,DIA,dSMD,random)
      endif
!
!
      
      return
!
!//////////////////////////////////////////////////////////////////////
      contains
!
!============================================================= 
      subroutine RosR_dis(INI,DIA,dSMD,random)
!=============================================================
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C subroutine for making
C Rosin-Rammler inverse distribution from random number random.
C A random number X must be 0.0 to 1.0
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      USE module_particle,only : dRR_a
!
      implicit none

      INTEGER,INTENT(IN) :: INI
      real*8,intent(in)  :: dSMD
      real*8,intent(in)  :: random
      real*8,intent(out) :: DIA

      if(dRR_a(INI)/=0.d0.and.random>0.d0.and.random<1.d0 ) then
        DIA=dSMD*(-log(1.d0-random))**(1.d0/dRR_a(INI))
      else
        DIA=0.d0
      endif

      return      
      return
      end subroutine RosR_dis
!
!============================================================= 
      subroutine normal_dis(INI,DIA,dSMD,random)
!=============================================================
      INTEGER,INTENT(IN) :: INI
      real*8,intent(in)  :: dSMD
      real*8,intent(in)  :: random
      real*8,intent(out) :: DIA
      return
      end subroutine normal_dis
!
      end subroutine drop_distrib
!
!

!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
      subroutine chem_admin_particle(
     &    IT,ITER_P,time,DT0,N_injctr,NINALL,NPALL,
     &    INSIDE,JOUT,PARCEL,DIA,TMPP,WDOTP,
     &    DYDTP,DHDTP,PMASS,PYS,RHOP,MASSYS,EVAP_H,EVAP_Y,
     &    MAT_NO,MAT_INDEX,MAT_CVEXT,MAT_CAL,LBC_SSF,LVEDGE,
     &    SFAREA,SFCENT,CVCENT,CVVOLM,
     &    vel,pp0,rho,aks,tmp,yys,wdot,tempp,
     &    zzs,cp,cr,FIELD_U,ierror)
!S$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle,lenfnm,getfil
      USE module_particle,only : partreac,IDX_PRTRAC
      USE module_particle,only : nneqP,INJ_RNUM,IFLAG_COAL,icoal
      USE module_particle,only : ireqP
      USE module_particle,only : stp_iterP,ig_iterP
      use module_particle,only : P_CHAR_COM,P_CHAR_NO,
     &                           P_HCN_NO_N2O2,
     &                           P_elemnt,P_ovrall,P_user
      use module_metrix,only   : lreacP
      use module_chemcntl,only : igT_iter,Fuel_NO_P,prompt_NO_P,
     6                           Zeldovich_P
      use module_species,only  : r_wm
      
!
      implicit none
!
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IT,ITER_P
      integer,intent(in)    :: N_injctr
      real*8 ,intent(in)    :: DT0,time
      integer,intent(in)    :: NINALL (0:N_injctr)
      INTEGER,INTENT(IN)    :: NPALL  (  N_injctr)
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_INDEX( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   ( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  ( 0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF(   MXSSFBC)
      INTEGER,INTENT(IN)    :: LVEDGE (2, MXCVFAC)
      REAL*8 ,INTENT(IN)    :: CVCENT (3, MXALLCV)
      REAL*8 ,INTENT(IN)    :: CVVOLM(    MXALLCV)
      REAL*8 ,INTENT(IN)    :: SFAREA (4, MXCVFAC)
      REAL*8 ,INTENT(IN)    :: SFCENT (3, MXCVFAC)
!
      INTEGER,INTENT(IN)    :: INSIDE (      MXPRT)
      INTEGER,INTENT(IN)    :: JOUT   (      MXPRT)
      REAL*8,INTENT(IN)     :: PARCEL (      MXPRT)
      REAL*8,INTENT(IN)     :: TMPP   (      MXPRT,2)
      REAL*8,INTENT(IN)     :: DIA    (      MXPRT)
      REAL*8,INTENT(INOUT)  :: DYDTP  (      MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: DHDTP  (      MXPRT,2)
      REAL*8,INTENT(IN)     :: PMASS  (      MXPRT,2)
      REAL*8,INTENT(IN)     :: PYS    (      MXPRT,MXCOMP,2)
      REAL*8,INTENT(INOUT)  :: RHOP   (      MXPRT)
      REAL*8,INTENT(INOUT)  :: WDOTP  (      MXCV,MXcomp)
!
      REAL*8, INTENT(IN)    :: VEL    (   MXALLCV,3,2)
      REAL*8 ,INTENT(IN)    :: pp0    (   MXALLCV)
      REAL*8, INTENT(IN)    :: RHO    (   MXALLCV,2)
      real*8, intent(in)    :: aks    (   MXALLCVR,MXRANS)
      REAL*8, INTENT(IN)    :: TMP    (   MXALLCV)
      REAL*8, INTENT(IN)    :: YYS    (   MXALLCV,MXCOMP)
      REAL*8, INTENT(INOUT) :: WDOT   (   MXALLCV,MXCOMP)
      REAL*8 ,INTENT(INOUT) :: zzs    (   MXALLCV,MXcomp)
      REAL*8,INTENT(INOUT)  :: cp     (   MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cr     (   MXALLCV)
      real*8 ,intent(inout) :: FIELD_U(   MXCV_D,NFLID)
      REAL*8,INTENT(INOUT)  :: MASSYS (   MXALLCV,MXcomp)
      REAL*8,INTENT(INOUT)  :: EVAP_H (   MXALLCV_P,2)
      REAL*8,INTENT(INOUT)  :: EVAP_Y (   MXALLCV_P,MXCOMP)
      REAL*8,INTENT(INOUT)  :: tempp  (   MXPRT,2)
!      REAL*8,INTENT(INOUT)  :: WDOTP_P(   MXPRT,MXcomp)
      integer,intent(inout) :: ierror

!
! --- [local entities]
!
      integer :: IIMAT,IMAT,IRC,ICVS,ICVE,ICVL,ICOM
      integer :: IBFS,IBFE,IBFL,ICFL,IDC,idum1,idum2
      integer :: iphg,IRCC 
      integer :: IPS,IPE,INI
      integer :: nb,kd,nbx,ntp,icoms,icome,iph,ipcom,icomm
      real*8,save  :: dum1,dum2
!      real*8,parameter :: PI=4.0*DATAN(1.d0)

!---------------------------------------------------------------------
      return
      end subroutine chem_admin_particle
!
!
!
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE USER_PARTICLE_INITIAL(INI,IPS,IPE,
     &  NCV,MXALLCV,MXPRT,MXCOMP,MXCVFAC,MXSSFBC,
     &  SFCENT,LBC_SSF,SFAREA,LVEDGE,CVCENT,
     *  DIA,XYZP,UVWP,RHOP,PARCEL,DT0,TMPP,PYS,
     &  PMASS,INSIDE,JOUT,HIS_P) 
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      USE module_particle,only : mass_p_total,iSprayInjNum,
     &                           iSprayInjInter,dInjVolume,dSLITWIDTH,
     &                           dXP0,dYP0,dZP0,iconst,dropkd,dSMD,RHOP0
!
      IMPLICIT NONE 
!
!--------------------
! Dummy Auguments 
!--------------------
!
      INTEGER, INTENT(IN)    :: INI,IPS,IPE
      INTEGER, INTENT(IN)    :: MXPRT,MXCOMP,MXCVFAC,MXSSFBC,NCV,MXALLCV
      REAL*8,  INTENT(IN)    :: DT0
      INTEGER ,intent(in)    :: LBC_SSF( MXSSFBC)
      INTEGER, INTENT(IN)    :: LVEDGE(2,MXCVFAC)
      REAL*8 , INTENT(IN)    :: SFAREA(4,MXCVFAC)
      REAL*8 , INTENT(IN)    :: CVCENT(3,MXALLCV)
      real*8 , intent(in)    :: SFCENT(3,MXCVFAC)
!
      INTEGER, INTENT(INOUT) :: INSIDE(  MXPRT)
      INTEGER, INTENT(INOUT) :: JOUT  (  MXPRT)
      REAL*8,  INTENT(INOUT) :: DIA     (MXPRT)
      REAL*8,  INTENT(INOUT) :: XYZP    (MXPRT,3)
      REAL*8,  INTENT(INOUT) :: UVWP    (MXPRT,3)
      REAL*8,  INTENT(INOUT) :: RHOP    (MXPRT)
      REAL*8,  INTENT(INOUT) :: PARCEL  (MXPRT)
      REAL*8,  INTENT(INOUT) :: PYS     (MXPRT,MXCOMP,2)
      REAL*8,  INTENT(INOUT) :: TMPP    (MXPRT,2)
      REAL*8,  INTENT(INOUT) :: PMASS   (MXPRT,2)
      REAL*8 , INTENT(INOUT) :: HIS_P   (MXPRT)
!
!--------------------
! --- Local Entities 
!--------------------
!
      INTEGER      :: I,ifl,ierr1=0,NPALL1,IP,NO,ICV,ICVL
      REAL*8,save  :: dum1,dum2,dum3,dum4,TEMP1,TEMP2,TEMP3
      REAL*8       :: PI,VOLUME,D3
      INTEGER,save :: NCVX,NLOC=0
!
      real*8,save,allocatable :: locpos(:)
!
! -------------------------------------------------------------------------
!
      PI=1.0D0
      PI=4.0*DATAN(PI)
!
! --- user section
!
      if(INI==1) then
        NCVX=NCV
        allocate(locpos(NCVX),stat=ierr1)
        if(ierr1/=0) call FFRABORT 
     &    (1,'ERR: locpos allocate error in USER_PARTICLE_INI')
        locpos(:)=0
        NLOC=0
        do ICV=1,NCV
        dum1=CVCENT(1,ICV)
        dum2=CVCENT(2,ICV)
        dum3=CVCENT(3,ICV)
!
!        if(dum1<-6.8d0.and.               dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-1
!        if(dum1>-6.2d0.and.dum1<-5.8.and.dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-2
!        if(dum1>-5.2d0.and.dum1<-4.8.and.dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-3
!        if(dum1>-4.2d0.and.dum1<-3.8.and.dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-4
!        if(dum1>-3.2d0.and.dum1<-2.8.and.dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-5
!        if(dum1>-2.2d0.and.dum1<-1.8.and.dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-6
!        if(dum1>-1.2d0.and.dum1<-0.8.and.dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-7
!        if(dum1>-0.2d0.and.dum1< 0.2.and.dum3>-0.15d0.and.dum3<0.0d0) then  !case-2-8
!
!        if(dum1<-6.8d0.and.               dum3>0.05d0.and.dum3<0.20d0) then !case-1-1
        if(dum1>-6.2d0.and.dum1<-5.8.and.dum3>0.05d0.and.dum3<0.20d0) then  !case-1-2
!        if(dum1>-5.2d0.and.dum1<-4.8.and.dum3>0.05d0.and.dum3<0.20d0) then  !case-1-3
!        if(dum1>-4.2d0.and.dum1<-3.8.and.dum3>0.05d0.and.dum3<0.20d0) then  !case-1-4
!        if(dum1>-3.2d0.and.dum1<-2.8.and.dum3>0.05d0.and.dum3<0.20d0) then  !case-1-5
!        if(dum1>-2.2d0.and.dum1<-1.8.and.dum3>0.05d0.and.dum3<0.20d0) then  !case-1-6
!        if(dum1>-1.2d0.and.dum1<-0.8.and.dum3>0.05d0.and.dum3<0.20d0) then  !case-1-7
!        if(dum1>-0.2d0.and.dum1< 0.2.and.dum3>0.05d0.and.dum3<0.20d0) then  !case-1-8
          NLOC=NLOC+1
          locpos(NLOC)=ICV
        endif 
        enddo 
      endif 
!
! --- 
!
      DO 1000 IP=IPS,IPE
!
!      if(IPCPU(IP)/=my_rank) cycle
! --- 
!
      if(dropkd(INI)/=iconst) then
        call drop_distrib(INI,DIA(IP),dSMD(INI),TEMP1)
      elseif(dropkd(INI)==iconst) then
        DIA(IP)=dSMD(INI)
      else
      endif 
!
      JOUT(IP)=0
      D3=DIA(IP)*DIA(IP)*DIA(IP)
      VOLUME=PI*D3/6.D0
      RHOP(IP)=RHOP0(INI)
      PMASS(IP,1)=VOLUME*RHOP(IP)
      mass_p_total(INI)=mass_p_total(INI)+PMASS(IP,1)*PARCEL(IP)
      PARCEL(IP)=DT0*dInjVolume(INI)/PMASS(IP,1) 
     &   /iSprayInjNum(INI)*iSprayInjInter(INI)
!
      call utl_random(TEMP2)
      call utl_random(TEMP3)
      ICV=int(1+TEMP2*NLOC)
      if(ICV<1.or.ICV>NLOC) then 
        ICV=locpos(1)
      endif
      ICVL=locpos(ICV)
      INSIDE(IP)=ICVL  
      
!------------------------
! --- initial cood
!------------------------
      XYZP(IP,1)=CVCENT(1,ICVL)
      XYZP(IP,2)=CVCENT(2,ICVL)
      XYZP(IP,3)=CVCENT(3,ICVL)
!------------------------
! --- initial velosity
!------------------------
      UVWP(IP,1)=0.d0
      UVWP(IP,2)=0.d0
      UVWP(IP,3)=0.d0
!
 1000 enddo

      return

      END SUBROUTINE USER_PARTICLE_INITIAL
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE USER_PARTICLE_TEMP1(MXPRT,IP,N_injctr,DELT_M,
     &           DP_M,TP_M,RHOP_M,MASSP_M,VELP_M,
     &           TMPP,DIAP,
     &           G_T_M,G_VEL_M,rmu_G_M,G_RHO_M,G_RAMD_M,G_CP_M
     &      )
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      IMPLICIT REAL*8 (A-H,O-Z)
!!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)    :: MXPRT,IP,N_injctr
      REAL*8,INTENT(IN)     :: DELT_M
      REAL*8,INTENT(INOUT)  :: DP_M,TP_M,RHOP_M,MASSP_M,VELP_M(3)
      REAL*8,INTENT(IN)     :: TMPP(MXPRT),DIAP(MXPRT)
!
      REAL*8,INTENT(IN)     :: G_T_M,G_VEL_M(3),rmu_G_M,G_RHO_M
      REAL*8,INTENT(IN)     :: G_RAMD_M,G_CP_M
!      
!
! --- [local entities]
!
      PARAMETER (ND=10)
      PARAMETER (N=10)
      REAL*8,save,allocatable  ::  T(:,:),XP(:,:)
      DIMENSION XPNEW(ND),TN(ND),XPP(N),TT(N)
      DIMENSION RHO(ND),CP(ND),RAMDA(ND)
      DIMENSION VELP(3),VELGAS(3)
!      REAL*8 :: 
      
      INTEGER,save  ::  iflagini=0
      INTEGER :: NSUB=1,NNSUB
!------------------------------------------------------------
!
      RHO_BATCH=1400.0
      CP_BATCH=1100.0
      RAMDA_BATCH=1.11
      RHO_GLASS=2600.0
      CP_GLASS=1300.0
      RAMDA_GLASS=30.0
!
! --- -------------------------
!
!      RHO_GAS=1.0
!      CP_GAS=1000.0
!      RAMDA_GAS=0.03
!      RMU_GAS=1.0D-4
!      TAIR=2000.0
!
!      VELP(1)=0.0
!      VELP(2)=0.0
!      VELP(3)=-1.0
!      VELGAS(1)=0.0
!      VELGAS(2)=0.0
!      VELGAS(3)=-2.0

!
! --- -------------------------
      RCO2=0.161
      RAIR=0.4D-3
      RENERGY=2230.0D+3
C
      A_PSAT=0.5D-3
      B_PSAT=0.0
      A_DIFF=2.0D-4
      B_DIFF=0.0
C
      TMELT=1200.0
!--------------------------------


!      DELT=0.00001
      QRAD=0.0
!      DP=2.0D-4
!      N=10
!      WRITE(6,*) 'particle diameter(mm)?'
!      READ(5,*) DP
!      DP=0.1d0
!      DP=DP*0.001
!
!
!
      if(iflagini==0) then
        ALLOCATE(T(N,MXPRT),XP(N,MXPRT))
        iflagini=1
        DO IIP=1,MXPRT
        DO I=1,N
         T(I,IIP)=TMPP(IIP)
        END DO
!
        DX=0.5*DIAP(IIP)/N
!
        DO I=1,N
         XP(I,IIP)=DX*I
        END DO
        enddo
      endif
C
!----------------------------------------
      DELT=DELT_M
      DELT_SUB=DELT/NSUB

      
      
      TAIR=G_T_M-273.15
      VELGAS(1:3)=G_VEL_M(1:3)
      VELP(1:3)=VELP_M(1:3)
      RMU_GAS=rmu_G_M
      RHO_GAS=G_RHO_M
      RAMDA_GAS=G_RAMD_M
      CP_GAS=G_CP_M
!----------------------------------------
C

C
      SUMQHEAT=0.0
      SUMCO2=0.0
      SUMENERGY=0.0
!      NTIME=0
!      TIME=0.0
!      NSTEP=0
!      WRITE(8,'(2A)') ' TIME AVE_TEMP TOT_MASS QHEAT VOLATIL_MAS '
!     1               ,' SUMQHEAT,SUMENERGY,SUMCO2 '
C
!  AVE_TEMP  ,TOT_MASS,RHO,DIA,
!  GAS_TEMP,VEL,UVWP,RHO,

  100 CONTINUE

      do 2000 NNSUB=1,NSUB
         

!      TIME=TIME+DELT_SUB
!      NSTEP=NSTEP+1
      XPP(:)=XP(:,IP)
      TT(:)=T(:,IP)
      CALL PARTICLE_HEAT
     1           (DELT_SUB,XPP,TT,QRAD,RHO,CP,RAMDA,N,VELP,VELGAS
     2           ,RHO_BATCH,CP_BATCH,RAMDA_BATCH
     3           ,RHO_GLASS,CP_GLASS,RAMDA_GLASS
     4           ,RHO_GAS,CP_GAS,RAMDA_GAS,RMU_GAS
     5           ,RCO2,RAIR,RENERGY,TMELT,TAIR
     6           ,A_PSAT,B_PSAT,A_DIFF,B_DIFF
     7           ,CO2,AIR,ENERGY,QHEAT
     8           ,VOLATIL_MAS,TN,XPNEW,AVE_TEMP,TOT_MASS)
      XP(:,IP)=XPP(:)
      T(:,IP)=TT(:)
C
      IF (VOLATIL_MAS.LT.1.0D-90) VOLATIL_MAS=0.0
      DO I=1,N
         T(I,IP)=TN(I)    !OK
         XP(I,IP)=XPNEW(I)   !OK
      END DO
      DP=XPNEW(N)*2.d0
C
      SUMQHEAT=SUMQHEAT+QHEAT*DELT_SUB   
      SUMCO2=SUMCO2+CO2*DELT_SUB
      SUMENERGY=SUMENERGY+ENERGY*DELT_SUB

 2000 enddo
!
!      WRITE(8,'(10G13.5)') TIME,AVE_TEMP,TOT_MASS,QHEAT,VOLATIL_MAS
!     1                    ,SUMQHEAT,SUMENERGY,SUMCO2
C
!
!
      DP_M=DP
      TP_M=AVE_TEMP+273.15
      RHOP_M=TOT_MASS/(DP**3*4.1315926d0/6.d0)
      MASSP_M=TOT_MASS
      
!
      return
!
      end SUBROUTINE USER_PARTICLE_TEMP1
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      SUBROUTINE PARTICLE_HEAT
     1           (DELT,XP,T,QRAD,RHO,CP,RAMDA,N,VELP,VELGAS
     2           ,RHO_BATCH,CP_BATCH,RAMDA_BATCH
     3           ,RHO_GLASS,CP_GLASS,RAMDA_GLASS
     4           ,RHO_GAS,CP_GAS,RAMDA_GAS,RMU_GAS
     5           ,RCO2,RAIR,RENERGY,TMELT,TAIR
     6           ,A_PSAT,B_PSAT,A_DIFF,B_DIFF
     7           ,CO2,AIR,ENERGY,QHEAT
     8           ,VOLATIL_MASS,TN,XPNEW,AVE_TEMP,TOT_MASS)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ND=20)
C
      DIMENSION XP(1),T(1),TN(1),RHO(1),CP(1),RAMDA(1),XPNEW(1)
      DIMENSION VELP(3),VELGAS(3)
C
      DIMENSION A(ND),B(ND),C(ND),R(ND),X(ND)
      DIMENSION XC(ND),DX(ND)
      DIMENSION OLDMASS(ND),OLDVOL(ND)
      DIMENSION ANEWMASS(ND),ANEWVOL(ND)
      DIMENSION IGLASS(ND)
C
      DATA TR,PATM,RCONST   / 273.15D-0, 101325.0D-0, 8314.3D-3 /
      DATA WNa2O            / 61.98D-3 /
      DATA PI               / 3.14159265D-0 /
C
C ----------------------- statment function
c      PSATX(TT)=0.5D-0
c      DIFF(TT)=2.0D-4
      PSATX(TT)=A_PSAT+B_PSAT*TT
      DIFF(TT)=A_DIFF+B_DIFF*TT
C
C
C------------------ 物性値の設定
C
      DP=2.0D-0*XP(N)
      VOL1=0.0D-0
      VOL2=4.0D-0/3.0D-0*PI*XP(1)**3
      DO I=1,N
         OLDVOL(I)=VOL2-VOL1
         VOL1=VOL2
         IF (I.LT.N) VOL2=4.0D-0/3.0D-0*PI*XP(I+1)**3
         IF (T(I).GT.TMELT) THEN
            IGLASS(I)=1
            RHO(I)=  RHO_GLASS
            CP(I)=   CP_GLASS
            RAMDA(I)=RAMDA_GLASS
         ELSE
            IGLASS(I)=0
            RHO(I)=  RHO_BATCH
            CP(I)=   CP_BATCH
            RAMDA(I)=RAMDA_BATCH
         END IF
         OLDMASS(I)=OLDVOL(I)*RHO(I)
         XPNEW(I)=XP(I)
      END DO
C
C----------------- 熱伝達量の計算
C
      AREA=PI*DP**2
      UREL=DSQRT( (VELP(1)-VELGAS(1))**2
     1           +(VELP(2)-VELGAS(2))**2
     2           +(VELP(3)-VELGAS(3))**2 )
      Re=UREL*DP*RHO_GAS/RMU_GAS
      Pr=RMU_GAS*CP_GAS/RAMDA_GAS
      ANu=2.0D-0 + 0.6D-0*Re**0.5D-0*Pr**(1.0D-0/3.0D-0)
      H=ANu*RAMDA_GAS/DP
      QHEAT=H*AREA*(TAIR-T(N))
C
C------------------------------------------------ 温度計算
C
      XC(1)=XP(1)*0.5D-0
      DX(1)=XP(1)
      DO I=2,N
         XC(I)=0.5D-0*(XP(I)+XP(I-1))
         DX(I)=XP(I)-XP(I-1)
      END DO
C
C----------------------------- I=1
      RHOCT=RHO(1)*CP(1)/DELT
      RAMI=RAMDA(1)
      RAMP1=2.0D-0/(1.0D-0/RAMDA(1)+1.0D-0/RAMDA(2))
      A(1)=0.0D-0
      B(1)=RHOCT+RAMP1/(DX(1)+0.5D-0*DX(2))/DX(1)
     1          +RAMI/XC(1)/DX(1)
      C(1)=     -RAMP1/(DX(1)+0.5D-0*DX(2))/DX(1)
     1          -RAMI/XC(1)/DX(1)
      R(1)=RHOCT*T(1)
C----------------------------- I=2,N-1
      DO I=2,N-1
         RHOCT=RHO(I)*CP(I)/DELT
         RAMI=RAMDA(I)
         RAMM1=2.0D-0/(1.0D-0/RAMDA(I)+1.0D-0/RAMDA(I-1))
         RAMP1=2.0D-0/(1.0D-0/RAMDA(I)+1.0D-0/RAMDA(I+1))
         A(I)=     -2.0D-0*RAMM1/(DX(I)+DX(I-1))/DX(I)
     1             +RAMI/XC(I)/DX(I)
         B(I)=RHOCT+2.0D-0*RAMP1/(DX(I)+DX(I+1))/DX(I)
     1             +2.0D-0*RAMM1/(DX(I)+DX(I-1))/DX(I)
         C(I)=     -2.0D-0*RAMP1/(DX(I)+DX(I+1))/DX(I)
     1             -RAMI/XC(I)/DX(I)
         R(I)=RHOCT*T(I)
      END DO
C----------------------------- I=N
         RHOCT=RHO(N)*CP(N)/DELT
         RAMI=RAMDA(N)
         RAMM1=2.0D-0/(1.0D-0/RAMDA(N)+1.0D-0/RAMDA(N-1))
         A(N)=     -2.0D-0*RAMM1/(DX(N)+DX(N-1))/DX(N)
     1             +RAMI/XC(N)/DX(N)
         B(N)=RHOCT+2.0D-0*RAMM1/(DX(N)+DX(N-1))/DX(N)
     1             -RAMI/XC(N)/DX(N)
         C(N)=0.0D-0
         VOL=4.0D-0/3.0D-0*PI*(XP(N)**3-XP(N-1)**3)
         R(N)=RHOCT*T(N)+(QHEAT+QRAD)/VOL
C
C----------------- 3重対角行列を解く
C
      CALL TRISOL(A,B,C,R,X,N,IFLG)
      DO I=1,N
         TN(I)=R(I)
      END DO
C
C------------------------------------------------ 温度計算 終了
C
C
C----------- CO2,AIR,ENERGY,AVE_TEMP,TOT_MASS,XPNEW の計算
C
      CO2=0.0D-0
      AIR=0.0D-0
      ENERGY=0.0D-0
      AVE_TEMP=0.0D-0
      TOT_MCP=0.0D-0
      TOT_MASS=0.0D-0
      DO I=1,N
         IF (IGLASS(I).EQ.0 .AND. TN(I).GT.TMELT) THEN
            CO2_MESH=RCO2*OLDMASS(I)
            AIR_MESH=RAIR*OLDMASS(I)
            ENERGY_MESH=RENERGY*OLDMASS(I)
            CO2=CO2+CO2_MESH
            AIR=AIR+AIR_MESH
            ENERGY=ENERGY+ENERGY_MESH
            ANEWMASS(I)=OLDMASS(I)-CO2_MESH-AIR_MESH
            ANEWVOL(I)=ANEWMASS(I)/RHO_GLASS
            TOT_MCP=TOT_MCP+ANEWMASS(I)*CP_GLASS
            AVE_TEMP=AVE_TEMP+TN(I)*ANEWMASS(I)*CP_GLASS
         ELSE
            ANEWMASS(I)=OLDMASS(I)
            ANEWVOL(I)=OLDVOL(I)
            TOT_MCP=TOT_MCP+ANEWMASS(I)*CP(I)
            AVE_TEMP=AVE_TEMP+TN(I)*ANEWMASS(I)*CP(I)
         END IF
         TOT_MASS=TOT_MASS+ANEWMASS(I)
      END DO
      CO2=CO2/DELT
      AIR=AIR/DELT
      ENERGY=ENERGY/DELT
      AVE_TEMP=AVE_TEMP/TOT_MCP
c      go to 55
C
      VOL=ANEWVOL(1)
      XPNEW(1)=(3.0D-0/4.0D-0*VOL/PI)**(1.0D-0/3.0D-0)
      DO I=2,N
         VOL=VOL+ANEWVOL(I)
         XPNEW(I)=(3.0D-0/4.0D-0*VOL/PI)**(1.0D-0/3.0D-0)
      END DO
C
C----------------- Na2O 揮発量(VOLATIL_MASS)の計算
C
      IF (IGLASS(N).EQ.1) THEN
         TT=TN(N)+TR
         PSAT=PSATX(TT)
         COEF_DIF=DIFF(TT)
         DP=2.0D-0*XPNEW(N)
         Re=UREL*DP*RHO_GAS/RMU_GAS
         Sc=RMU_GAS/RHO_GAS/COEF_DIF
         Sh=2.0D-0 + 0.6D-0*Re**0.5D-0*Sc**(1.0D-0/3.0D-0)
         Hd=Sh*COEF_DIF/DP
         AREA=PI*DP**2
         Csurf=PSAT*PATM/(RCONST*TT)*WNA2O
         VOLATIL_MASS=VOLATIL_MASS+Hd*AREA*Csurf*DELT
      END IF
   55 continue
C
C
      RETURN
      END
C
C
C--------------------------------------------sub ３重対角行列
C
C
      SUBROUTINE TRISOL(A,B,C,R,X,N,IFLG)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(1),B(1),C(1),R(1),X(1)
C
      IFLG=0
      IF (DABS(B(N)).LT.1.0D-30) THEN
         IFLG=1
         RETURN
      END IF
C
      A(N) = A(N)/B(N)
      R(N) = R(N)/B(N)
      DO 100 I=2,N
         II = -I + N + 2
         BN = 1./(B(II-1) - A(II)*C(II-1))
         IF(II.GT.2) A(II-1) = A(II-1) * BN
         R(II-1) = (R(II-1)-C(II-1)*R(II)) * BN
100   CONTINUE
      X(1) = R(1)
      DO 200 I=2,N
         X(I) = R(I) - A(I)*X(I-1)
200   CONTINUE
      DO 300 I=1,N
         R(I) = X(I)
300   CONTINUE
C
C
      RETURN
      END
