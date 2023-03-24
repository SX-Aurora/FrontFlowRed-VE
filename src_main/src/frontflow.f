!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 3.1)                           !
!                                                                      !
!                 This code is re-developed by NuFD CO.LTD             !
!     Main Program : FRONTFLOW                                         !
!                          Both Cell-Center and Vertex-Center          !
!                                                                      !
!                          Developer: Huilai   ZHANG                   !
!                          Adviser  :                                  !
!                                     Nobuyuki OSHIMA                  !
!                                     Makoto   TSUBOKURA               !
!                                     Takuji   NAKASHIMA               !
!                                     Toshio   KOBAYASHI               !
!                                     Yuichi   ITOH                    !
!                                     2008/06/15                       !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 3.0)                           !
!                                                                      !
!     Main Program : FRONTFLOW   (Vector Routine (for vector computer) !
!                                 Lagrangian 2P [Particle Tracking]    !
!                                 LES Flamelet conbustion model        ! 
!                                 Radiation [FVM,MCM & ZONE]           ! 
!                                 Dynamical SGS modle                  !  
!                                                                      !
!                          Written by Huilai   ZHANG                   !
!                                     Yuyan    Jiang                   !
!                                     Takuji   TOMINAGA                !
!                                     Nobuyuki OSHIMA                  !
!                                     Makoto   TSUBOKURA               !
!                                     2007/06/01                       !
!                                                                      !
!     Contact address : The University of Tokyo, RSS21 Project         !
!                                                                      !
!======================================================================!
!
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 3.0beta)                       !
!                                                                      !
!     Main Program : FRONTFLOW   (Vector Routine (for vector computer) !
!                                 Lagrangian 2P [Particle Tracking]    !
!                                 LES Flamelet conbustion model        ! 
!                                 Radiation [FVM,MCM & ZONE]           ! 
!                                 Dynamical SGS modle                  !  
!                                                                      !
!                          Written by Huilai   ZHANG                   !
!                                     Yuyan    Jiang                   !
!                                     Takuji   TOMINAGA                !
!                                     Takeshi  UNEMURA                 !
!                                     Nobuyuki OSHIMA                  !
!                                     Makoto   TSUBOKURA               !
!                                     2007/06/01                       !
!                                                                      !
!     Contact address : The University of Tokyo, RSS21 Project         !
!                                                                      !
!======================================================================!
! 
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 2.8a)                          !
!                                                                      !
!     Main Program : FRONTFLOW    (Vector Routine, DES,                !
!                                      Serval RANS Models)             !
!                                                                      !
!                          Written by Huilai   ZHANG                   !
!                                     Takuji   TOMINAGA                !
!                                     Keiji    ONISHI                  !
!                                     Eisuke   YAMADA,                 !
!                                     Takafumi SUGINAKA                !
!                                     Yuyan    Jiang                   !
!                                     Takeshi  UNEMURA                 !
!                                     Nobuyuki OSHIMA                  !
!                                     Makoto   TSUBOKURA               !
!                                     2006/06/15                       !
!                                                                      !
!     Contact address : The University of Tokyo, RSS21 Project         !
!                                                                      !
!======================================================================!
! 
!======================================================================!
!                                                                      !
! Software Name : FrontFlow/red   (Ver. 2.0)                           !
!                                                                      !
!     Main Program : FRONTFLOW    (Vertex-Center,EDGE-BASE,HPC,        !
!                                  Multi-Domain[thermal solid coupled])!
!                                                                      !
!                          Written by Huilai ZHANG                     !
!                                     Takeshi UNEMURA                  !
!                                     Nobuyuki TANIGUCHI               !
!                                     2005/05/26                       !
!                                                                      !
!     Contact address : The University of Tokyo, FSIS Project          !
!                                                                      !
!======================================================================!
!
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
!                                     2004/03/25                       !
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
!      MODULE     MODULE_ARRAY_FFLOW
!      MODULE     MODULE_ELEMENT_FFLOW
!      PROGRAM    FRONTFLOW
!      SUBROUTINE INPUT_FFLOW
!      SUBROUTINE MAIN_FFLOW
!      SUBROUTINE SUF_ADMIN
!      SUBROUTINE HYS_ADMIN
!      SUBROUTINE RANS_ADMIN
!      SUBROUTINE VEL_ADMIN
!      SUBROUTINE PRS_ADMIN
!      SUBROUTINE RVA_ADMIN
!      SUBROUTINE ACS_MONITOR
!      SUBROUTINE ERROR_STP
!      SUBROUTINE ERROR_SMAC
!      SUBROUTINE INITVAR
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine INPUT_FFLOW(Dmns,ivector,iMHD,iMvmsh,ivof,iBODY,FLMLT,
     &           iCAVIT,iPEFC,iprtcle,iinj,ifld,MXPRTX,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      USE MODULE_HPCUTIL
      USE module_dimension,only : nrans,ncomp,ncomp_suf,!ncomp_prt
     &                            KE2P,KSUF,
     &                            ncompall,nphase,NPHS,NPOTN
      USE module_dimension,ONLY : NCV,NCVIN,NCVFAC,MXALLCV,
     &                            NMAT,NSSFBC,MXMAT
      use module_metrix   ,only : cofleft =>W1K8
      use module_metrix   ,only : net_vreq=>W1K9
      use module_chemreac ,only : UNIT_SI,nneq,const_R,M_3rd,P_3rd
      use module_species  ,only : gascns,gascns_1
!
      USE MODULE_IO       ,ONLY : IFLI,IFLL,IFLE,OPENFILES,cntlnam
      USE MODULE_IO       ,ONLY : INPUT_IO       =>INPUTDATA
      USE MODULE_MODEL    ,ONLY : INPUT_MODEL    =>INPUTDATA
      USE MODULE_FLAGS    ,ONLY : INPUT_FLAGS    =>INPUTDATA
      USE MODULE_TIME     ,ONLY : INPUT_TIME     =>INPUTDATA
      USE MODULE_DELTAT   ,ONLY : INPUT_DELTAT   =>INPUTDATA
      USE MODULE_SIMPLE   ,ONLY : INPUT_SIMPLE   =>INPUTDATA
      USE MODULE_GRAVITY  ,ONLY : INPUT_GRAVITY  =>INPUTDATA
      USE MODULE_CGSOLVER ,ONLY : INPUT_CGSOLVER =>INPUTDATA
      USE MODULE_RANS     ,ONLY : INPUT_ransMODEL=>INPUTKE
      USE MODULE_LES      ,ONLY : INPUT_LES      =>INPUTDATA
      USE MODULE_SPECIES  ,ONLY : INPUT_SPECIES_gas  =>INPUTDATA_gas
      USE MODULE_SPECIES  ,ONLY : INPUT_SPECIES_suf  =>INPUTDATA_suf
      USE MODULE_CHEMCNTL ,ONLY : INPUT_CHEMCNTL =>INPUTDATA
      USE MODULE_CHEMREAC ,ONLY : INPUT_CHEMREAC =>INPUTDATA
      USE MODULE_BOUNDARY ,ONLY : INPUT_BOUNDARY =>INPUTDATA
      USE MODULE_SOURCE   ,ONLY : INPUT_SOURCE   =>INPUTDATA
      USE module_material ,ONLY : INPUT_MAT      =>INPUTDATA
      USE MODULE_MOVEGRID ,ONLY : INPUT_MOVEGRID =>INPUTDATA
      USE MODULE_OUTPUT   ,ONLY : INPUT_OUTPUT   =>INPUTDATA
      USE MODULE_DIMNSN   ,ONLY : INPUT_DIMNSN   =>INPUTDATA
      USE MODULE_DEBUG    ,ONLY : INPUT_DEBUG    =>INPUTDATA
      USE MODULE_USERSUB  ,ONLY : INPUT_USER     =>INPUTDATA
      USE MODULE_HPC      ,ONLY : INPUT_HPC      =>INPUTDATA
      use module_hpc_input,only : hpc_input      =>inputdata
      use module_flow_sound,only: INPUT_sound    =>inputdata
      use module_flow_sound1,only: INPUT_sound1  =>inputdata
      use module_anim     ,only : input_anim     =>inputdata
      use module_size     ,only : nml_sizes
      USE MODULE_BOUNDARY ,ONLY : nbcnd,boundName,kdbcnd,kxnone,
     &                            kdfire,kdintr,kdbuff,kdcvd,
     &                            kvnslp,kvlglw,kvfslp,kvmodl,kvEWF,
     &                            kdshutr,lsldbc,lovstbc
      use module_Euler2ph ,ONLY : input_eul2ph   =>inputdata
      use module_vof      ,ONLY : input_vof      =>inputdata
      use module_scalar   ,ONLY : input_scalar   =>inputdata
      use module_cdcl     ,only : input_cdcl     =>inputdata
      use module_probe    ,only : input_prob     =>inputdata
      use module_fluidforce,only: input_fforce   =>inputdata
      use module_FUEL      ,only: input_FUEL     =>inputdata
      use module_particle  ,only: input_particle =>inputdata
      use module_particle  ,only: particle_chem

      use module_rad,only       : input_radiation =>inputdata
!
      USE module_material ,ONLY : ical_sld
      use module_chemreac ,ONLY : ical_suf,nneq,igas,isurface,
     &                            kind_chem,preq,lreq,
     &                            ovalfr,ovrall,stick_sf,
     &                            BohmE,BohmY,Bohm,userreac,
     &                            creq,ireq,vreq,elemarbi,
     &                            kind_chem,
     &                            kind_prs,falloff
      USE MODULE_CHEMCNTL ,ONLY : igT_iter
      use module_model,only     : nthrds,encomp,ical_t,mach0,idrdp
      use module_scalar   ,ONLY : calgeq,calxi,caltemp,calh,
     &                            igeq,ixi,itemp,ihflmlt,ical_prp,
     &                            ical_flmlt
!
      implicit none
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(OUT)    :: IERROR
      INTEGER,INTENT(INOUT)  :: ivector,iMHD,iMvmsh,iBODY,FLMLT,iCAVIT,
     &                          ifld,iPEFC,iprtcle,iinj,ivof,MXPRTX
      character(len=80),intent(out) :: Dmns
      integer,parameter      :: mcomp=200,mrans=50,mneq=300 
!see module_boundary
      INTEGER                :: chem_bc(mneq)=0
!
! --- [local entities]
!
      logical :: lcomp,lrans,E2P=.false.,kemdl,RSMmdl,firewal=.false.,
     &           lvof=.false.,
     &           lrot=.false.,
     6           RANSMDL=.false.,
!     &           lsldbc=.false.,
!     &           lovstbc=.false.,
     &           lsuf=.false.,lowke=.false.,
     &           lRNG=.false.,lCHEN=.false.,lmhd=.false.,
     &           lSDES=.false.,lcavi=.false.,lpotn=.false.,
     &           l2SKE=.false.,lKLES=.false.
      logical :: u_rgn(1:200)=.false.
      integer :: lFC=0,vapor=1,cool=7,coal_CHO(1:mneq)=0
      integer :: rg_markx=0,restart=0,ip_mdl=0,icomp=0
      integer :: idens=0,ireac=0,lflm=0
      integer :: MIX_KNI=0
      integer :: CHUN=0
      integer :: alpl(2),ikel(2),jvof,
     &           NOMAT=1
      integer :: l_temp=0
      character(len=20) :: spinam(mcomp)
      integer :: blktyp(mcomp)
      real*8  :: dum,dum1,dum2
      real*8,parameter :: SML=1.d-10
      integer :: q,s
      real*8 ,parameter :: undef = -huge(1.d0)
      integer :: nphasex=0
      
!      integer :: I_CHO=0,I_H2O=0,I_C=0,I_ash=0,I_O2
!
! --- ------------------------------------------------------------------------
!
! license check
!
      character(len=10) :: featureName = "FFR"
      character(len=0)  :: version = ""
!

      CALL checkExpire(2022,9,30,23,59,59,ierror)
      if(ierror/=0) call FFRABORT(1,'FRONTFLOW')
!----------------------------------------
! --- ncompall=ncomp+ncomp_suf;nallcomp
!----------------------------------------
      ALLOCATE(cofleft(1:ncompall*mneq))
      ALLOCATE(net_vreq(1:ncompall*mneq))
      cofleft=0.d0
      net_vreq=0.d0
!
      if(mcomp<ncompall) then
        call FFRABORT(1,'ERR: mcomp<ncompall, call your supportor')
      endif
!
! --- < file names & open files >-
!
      ifli=10   !huilai
      if(ifli<=0) ifli=1
      open(unit=ifli,file=cntlnam,status='unknown',
     &     action="read",iostat=ierror)
      if(ierror>0) then
        write(ifle,*) '*** can not open [',cntlnam,'] file',
     &            'ifli  =', ifli,
     &            'iostat=', ierror   
        goto 9999
      end if
!
      call input_io(ifli,ierror)
      if(ierror.ne.0) goto 9999
      call openfiles(ierror)
      if(ierror.ne.0) goto 9999
!----------------------------
! --- < model flags >-
!----------------------------
      call input_model
     &  (ifli,ifll,ifle,cntlnam,lcomp,l_temp,kemdl,RSMmdl,lowke,lKLES,
     &   l2SKE,lRNG,lCHEN,lSDES,RANSMDL,ivector,
     &   iMHD,iMvmsh,iBODY,lFC,ifld,
     &   iprtcle,ip_mdl,icomp,idens,ireac,ierror)
      if( ierror.ne.0 ) goto 9999
!
      call nml_sizes(ifli,ifll,ifle,cntlnam,MXPRTX,iprtcle,ierror)
      if( ierror.ne.0 ) goto 9999
!----------------------------
!-< Euler two Phase flow : >
!----------------------------
      call input_eul2ph(ifli,ifll,ifle,KE2P,NPHS,cntlnam,E2P,ierror)
      if( ierror.ne.0 ) goto 9999
!------------------------------------------
!-< VOF two Phase flow : defining [nrans]>
!------------------------------------------
      call input_vof(ifli,ifll,ifle,cntlnam,lvof,ivof,ierror)
!---------------------------------------
!-< scalar equation : defining [nrans]>
!---------------------------------------
      call input_scalar
     &     (ifli,ifll,ifle,cntlnam,my_rank,
     &      lcomp,lrans,E2P,kemdl,lowke,l2SKE,lRNG,lCHEN,lSDES,lKLES,
     &      RSMmdl,RANSMDL,lvof,ncomp,nrans,
     &      NPHS,NPOTN,lcavi,lFC,lpotn,
     &      alpl,ikel,jvof,lflm,FLMLT,iCAVIT,iPEFC,ierror)
      if(lflm==1) then
        encomp=0
        l_temp=1
        ical_t=.true.
        idrdp=mach0
      endif
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< data for every material >-
!---------------------------------------
      call input_MAT(ifli,ifll,ifle,my_rank,NOMAT,cntlnam,
     &     lrans,E2P,lvof,lcavi,lrot,imhd,iMvmsh,lFC,
     &     MIX_KNI,CHUN,
     &     ncompall,ivector,rg_markx,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
! --- < some flags >-
!---------------------------------------
      call input_flags(ifli,ifll,ifle,cntlnam,lrans,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< start/end of time integrarion >-
!---------------------------------------
      call input_time(ifli,ifll,ifle,cntlnam,iMHD,
     &  restart,ip_mdl,icomp,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< data to estimate time increment >-
!---------------------------------------
      call input_deltat(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< control data for SIMPLE >-
!---------------------------------------
      call input_simple(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< parameter for cg solver >-
!---------------------------------------
      call input_cgsolver(ifli,ifll,ifle,my_rank,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< parameter of k-e model >-
!---------------------------------------
      call input_ransmodel(ifli,ifll,ifle,cntlnam,
     &     nrans,kemdl,lowke,l2SKE,lRNG,lCHEN,lSDES,lKLES,E2P,
     &     calgeq,calxi,caltemp,calh,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< gas species data >- 
!---------------------------------------
      call input_species_gas(ifli,ifll,ifle,mcomp,
     &         ncomp,ncomp_suf,
     &         cntlnam,my_rank,
     &         spinam,E2P,lvof,lcavi,
     &         lFC,vapor,cool,MIX_KNI,CHUN,idens,
     &         ierror)
      if( ierror.ne.0 ) goto 9999
!
!---------------------------------------
!-< surface species data >-
!---------------------------------------
      call input_species_suf
     &   (ifli,ifll,ifle,cntlnam,my_rank,mcomp,
     &    ncomp,ncomp_suf,spinam,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< data for gravity >-
!---------------------------------------
      call input_gravity(ifli,ifll,ifle,ncomp,cntlnam,lcomp,E2P,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< chemical reaction >-
!---------------------------------------
      call input_chemcntl(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!

      call input_chemreac(ifli,ifll,ifle,my_rank,mcomp,lsuf,
     &    cntlnam,ncomp,ncomp_suf,lcomp,igT_iter,
     &    mneq,ncompall,chem_bc,blktyp,
     &    cofleft,net_vreq,coal_CHO,ireac,ierror)
      if( ierror.ne.0 ) goto 9999
!----------------------------
!-< parameter of LES model >-
!----------------------------
      call input_les
     &  (ifli,ifll,ifle,mcomp,mneq,ncompall,ncomp,nneq,
     &   my_rank,cntlnam,nrans,spinam,net_vreq,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< boundary condition >-
!---------------------------------------
      call input_boundary
     &  (ifli,ifll,ifle,my_rank,cntlnam,mcomp,ncomp,ncomp_suf,
     &   nphase,nphasex,nrans,NPOTN,nneq,spinam,
     &   lrans,kemdl,lowke,E2P,firewal,lrot,lsuf,isurface,lcomp,
     &   mneq,ncompall,cofleft,net_vreq,blktyp,
     &   chem_bc,alpl,ical_sld,ical_suf,KSUF,NOMAT,imhd,lFC,
     &   lpotn,iprtcle,idens,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
! --- cm unit: A=cm^3*(n-1)/[mole^(n-1)*s*K^(-alpha)]; E=kcal/mole
!---------------------------------------
      if(nneq>0) then
        const_R(:)=gascns
        do q=1,nneq
        if(UNIT_SI(q)==1) then
          const_R(q)=gascns_1
          if(ireq(q)==stick_sf.or.
     &     ireq(q)==BohmE.or.
     &     ireq(q)==BohmY.or.
     &     ireq(q)==Bohm.or.
     &     ireq(q)==userreac
     &     ) cycle
          dum=0.d0
          do s=1,ncomp+ncomp_suf
          if(abs(vreq(s,1,q))>SML) then
            dum1=vreq(s,1,q)
            if(ireq(q)==ovrall.or.
     &        ireq(q)==ovalfr) dum1=creq(s,q)
            if(blktyp(s)==1) then
              dum=dum+6.d0*dum1
            elseif(blktyp(s)==2) then
              dum=dum+4.d0*dum1
            elseif(blktyp(s)==3) then
            endif
          endif
          enddo
!
! --- 
!
          dum2=0.d0
          do s=1,ncomp+ncomp_suf
          if(abs(vreq(s,2,q))>SML) then
            dum1=vreq(s,2,q)
            if(ireq(q)==ovrall.or.
     &        ireq(q)==ovalfr) dum1=creq(s,q)
            if(blktyp(s)==1) then
              dum2=dum2+6.d0*dum1
            elseif(blktyp(s)==2) then
              dum2=dum2+4.d0*dum1
            elseif(blktyp(s)==3) then
            endif
          endif
          enddo
!
          if(M_3rd(q)==1.and.
     &       kind_chem(q)==igas.and.
     &       P_3rd(q,1)/=1) then
            dum=dum+6.d0
            dum2=dum2+6.d0
          endif
!
          if(M_3rd(q)==1.and.kind_chem(q)==isurface) then
            dum=dum+4.d0
            dum2=dum2+4.d0 
          endif
!
          if(kind_chem(q)==igas) then
            if(kind_prs(q)==0) then
              preq(1,1,q)=preq(1,1,q)*(10.d0**(-dum+6.d0))
              if(preq(1,2,q)/=undef) then
                preq(1,2,q)=preq(1,2,q)*(10.d0**(-dum2+6.d0))
              endif
!ogasa 
!              if(q.eq.7 .or. q.eq.8) then
!                preq(1,2,q)=preq(1,2,q)*(10.d0**(-dum2+6.d0))
!              endif

            else
              if(kind_prs(q)==falloff) then
                preq(1,1,q)=preq(1,1,q)*(10.d0**(-dum+6.d0))
                lreq(1,1,q)=lreq(1,1,q)*(10.d0**(-dum))
                if(preq(1,2,q)/=undef) then 
                  preq(1,2,q)=preq(1,2,q)*(10.d0**(-dum2+6.d0))
                endif
                if(lreq(1,2,q)/=undef) then
                  lreq(1,2,q)=lreq(1,2,q)*(10.d0**(-dum2))
                endif
              else
                lreq(1,1,q)=lreq(1,1,q)*(10.d0**(-dum+6.d0))
                preq(1,1,q)=preq(1,1,q)*(10.d0**(-dum+12.d0))
                
                if(lreq(1,2,q)/=undef) then
                  lreq(1,2,q)=lreq(1,1,q)*(10.d0**(-dum2+6.d0))
                endif
                if(preq(1,2,q)/=undef) then 
                  preq(1,2,q)=preq(1,1,q)*(10.d0**(-dum2+12.d0))
                endif
              endif
            endif
          elseif(kind_chem(q)==isurface) then
            preq(1,1,q)=preq(1,1,q)*(10.d0**(-dum+4.d0))
          endif
          if(my_rank==ROOT) then
            write(*,'(2x,a,I4,a,E12.5,a)') 
     &   'MSG: Reaction number= ',q,' A=',preq(1,1,q),
     &   ' in unit of M-mol-J-K'
          endif
        endif
        enddo
      endif
!-----------------
!-< source term >-
!-----------------
      call input_source
     &   (ifli,ifll,ifle,cntlnam,ncomp,nrans,lrans,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< control data for moving grid >-
!---------------------------------------
      call input_movegrid(ifli,ifll,ifle,cntlnam,lrot,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< control data to output >-
!---------------------------------------
      call input_output(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< dimension of computational domain >-
!---------------------------------------
      call input_dimnsn(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< debug flags >-
!---------------------------------------
      call input_debug(ifli,ifll,ifle,cntlnam,ierror)
      if( ierror.ne.0 ) goto 9999
!----------------------
!-< User subroutine >-
!----------------------
      call input_user(ifli,ifll,ifle,cntlnam,firewal,imhd,ierror)
      if( ierror.ne.0 ) goto 9999
!------------------
!-< HPC in fort.1>-
!------------------
      call input_hpc(ifli,ifll,ifle,cntlnam,NPE,ierror)
      if( ierror.ne.0 ) goto 9999
!--------------------
!-< HPC in CONT.HPC>-
!--------------------
      if(NPE.gt.1) then
         call hpc_input(my_rank,NPE,ifli,ifll,ifle,cntlnam,Dmns,ierror)
      endif
      if(ierror.ne.0) goto 9999
!---------------------------------------
!-< flow sound>-
!---------------------------------------
!
      call INPUT_sound
     & (ifli,ifll,ifle,cntlnam,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!-< animation >-
!---------------------------------------
      call input_anim(ifli,ifll,ifle,cntlnam,nrans,ncomp,mcomp,
     &                ncomp_suf,spinam,ierror)
      if( ierror.ne.0 ) goto 9999
!
!---------------------------------------------
!-<CD/CL value>-
!---------------------------------------------
      call INPUT_cdcl
     & (ifli,ifll,ifle,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,ierror)
      if(ierror.ne.0) goto 9999
!---------------------------------------------
!-<Probe data>-
!---------------------------------------------
      call INPUT_prob
     & (ifli,ifll,ifle,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,ierror)
      if(ierror.ne.0) goto 9999
!
!--------------------
!-< Wall variables >-
!--------------------
      call input_fforce     !zhang0119
     & (ifli,ifll,ifle,cntlnam,nbcnd,boundName,kdbcnd,kxnone,
     &  kdfire,kdintr,kdbuff,kdcvd,kvnslp,kvlglw,kvfslp,kvmodl,kvEWF,
     &  ierror)
      if(ierror.ne.0) goto 9999
!------------------
! --- FUEL
!------------------
      call input_FUEL(ifli,ifll,ifle,ncompall,cntlnam,
     &     lFC,vapor,cool,ierror)
      if(ierror.ne.0) goto 9999
!
!------------------
! --- particle 
!------------------
!

      call input_particle(ifli,ifll,ifle,cntlnam,restart,
     &     iprtcle,iinj,ncomp,mneq,l_temp,lcomp,
     &     my_rank,coal_CHO,nneq,
     &     ireac,igT_iter,
     &     spinam,ierror) 
!
!---------------------------------------
!-< data for radiation heat transfer >-		
!---------------------------------------
      call input_radiation(ifli,ifll,ifle,cntlnam,rg_markx,ierror)
!
      close(ifli)
!
      deALLOCATE(cofleft,net_vreq)
!
      return
!
 9999 continue
!
      if(my_rank.eq.ROOT) write(ifle,*) '(INPUT_FFLOW)'
      ierror=1
      end subroutine INPUT_FFLOW
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine main_FFLOW(
     &  LVEDGE,LVRTCV,LFUTAU,
     &  LBC_SSF,LCYCSF,FRSTCV,UTAU,DISINL,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,LCV_CV,
     &  Pstart,mat_cal,
     &  SFAREA,SFCENT,CVVOLM,CVVOL0,CVCENT,wiface,DISALL,xta,
     &  rva,rho,prs,pp0,tmp,yys,vel,velf,aks,rmut,rmu,rmd,rds,
     &  wdot,ccc,
     &  rmx,hhh,hhs,cps,cp,cr,diag,kdbt,kdby,
     &  rvx,rvd,grdc,
!     &  vlasum,dvdd,vctr,vctrp,FIELD_U,
     &  vlasum,dvdd,vctrp,FIELD_U,
     &  prsbnd,tmpbnd,yysbnd,velbnd,aksbnd,
     &  htcbnd,mtcbnd,radbnd,HTFLUX,
     &  SUFBND,SDOT,MOLFRC,MOLCEN,SITDEN,BLKTHK,!WDRBND,
     &  LCYCOLD,wifsld,OPPANG,locmsh,
     &  ipfix,itery,iterk,iterv,iterp,repsy,repsk,repsv,repsp,
     &  aepsy,aepsk,aepsv,aepsp,erry,errk,errv,errp,
     &  errsdy,errsdk,errsdv,
!MHD
     &  reps_FAI,aeps_FAI,err_FAI,
     &  reps_A,aeps_A,err_A,iter_FAI,iter_A,errsdA,errsdFAI,
!Particle
     1  DNORM,NEIGHBOR,
     2  INSIDE,JOUT,PMASS,DIA,XYZP,UVWP,PARCEL,RHOP,PYS,TMPP,
     3  HEIGHT,FU,FV,FW,
     4  MOMCELL,EVAP_Y,EVAP_H,IPCELL,NEIBCELL,DYDTP,DHDTP,HIS_P,
! Euler two phase
     &  VEL2,PRS2,tmp2,RVA2,RHO2,RMUT2,RMU2,
     &  YYS2,RDS2,WDOT2,RMD2,ccc2,dvdd2,
     &  HEATL_C,HEATG_C,HLS,HGS,TMP_SAT,FRIC_C,MASS_I,UTAU2,hhh2,
     &  VELBND2,TMPBND2,YYSBND2,
     &  HTCBND2,MTCBND2,RADBND2,
! potential scalar
     &  iptfix,POTNAL,POFLUX,POTNBC,
     &  iter_POTN,reps_POTN,aeps_POTN,err_POTN,
!
     &  itery2,iterv2,repsy2,repsv2,
     &  aepsy2,aepsv2,erry2,errv2,
     &  errsdy2,errsdv2,
!
     &  RadHeatFlux,P_rad,PAbsorb,PScatter,
     &  RadHeat,
     &  ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension 
      use module_hpcutil
      use module_constant
      use module_io,      only : ifll,ifle,ifli
      use module_nowtime, only : iter,time
      use module_movegrid,only : ipiston,endpiston
      use module_model,   only : idrdp,incomp,mach0,weak,icaltb,
     &                           ke,sles,dles,lles,dns,RSM,encomp,
     &                           ical_t,ical_topt,ical_reac,comp,
     &                           ical_surf,ical_week,ical_vect,
     &                           nthrds,monitor_stp,ical_MHD,
     &                           ical_mvmsh,iLBF_P,ical_dens,
     &                           PEFC,PAFC,MCFC,SOFC,AFC,MDFC,
     &                           iheat,ical_prt
      use module_time,    only : iters,itere,timee,steady,toltim,
     &                           i_steady
      use module_simple,  only : nsmpl,simeer,I_PISO
      use module_boundary,only : bc_setnow => setnow,
     &                           MAT_BCIDX,kdpres,kdolet,wdbcnd,kdbcnd,
     &                           nbcnd,kdcvd,heat,nobcnd,surfreac,
     &                           lovstbc ,imasflg,ifixmach
      use module_source,  only : src_setnow => setnow 
      use module_flags,   only : intgvv,intgty,intgke,eulere,euleri,
     &                           Adams_Bashforth,Adams_Moulton,
     &                           Runge_Kutta,Crank_Nicolson,icalrv,
     &                           face
      use module_species,only  : sw,spcnam,wm
      use module_usersub, only : usrno,usryes,outusr
      use module_param
      use module_deltat,  only : ideltt,const
      use module_chemreac,only : nneq,ig_iter,ical_suf
      use module_chemcntl,only : igT_iter
      use module_Euler2ph,only : ieul2ph
      use module_VOF     ,only : ical_vof
      use module_scalar,  only : icalke,icalrsm,icalaph,
     &                           ike,irsm,iaph,ivof,
     &                           rns_scl,ical_cavi,ivold,icavi,
     &                           calgeq,calxi,caltemp,calh,
     &                           ical_FC,ical_s,
     &                           pot_scl,iterPOTEN
      use module_material,only : lclsd,iclosd,nofld,
     &                           nosld,nflud,ical_sld,incmp_prs
      USE module_movegrid,only : ical_mov
      use module_metrix,only   : MHD_A,MHD_FAI,MHD_RVA,SIGMA,MHDbnd,
     &                           MHD_CRNT,MHD_CRT0,tempMHD,MHD_FAI0
      use module_metrix,only   : movflg,cord,lacell,lvcell,
     &                           lvface,lbface,lcface,lfcell,LEFACE,
     &                           listbc,
     &                           area,volume,gface,gcell
      use module_metrix,only   : bdyfrc,vflc,yplusf,rva_s
      use module_metrix,only   : t_dir,ANGBND
      use module_boundary,only : 
     &                           kdintr,ktneum,kttrns, 
     &                           nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX
      use module_rad   , only  : radflag,radmodel,radjump
      use module_radsxf, only  : sumcoef,radinten,GWAbsorb
!     &                          RadHeatFlux,P_rad,
!     &                           PAbsorb,PScatter,
!     &                           radinten,GWAbsorb,RadHeat
      use module_particle,only : ITER_P,DIMV,N_injctr,NINALL,NPALL,
     &                           ical_P_R,nneqP
      use module_gravity  ,only : ggg
      use module_metrix,only    : msk,vctr
      use MODULE_ARRAY_FFLOW ,only : iterR,repsR,aepsR,errR
      use module_metrix,only    : SHUTFL
      use module_scalar,only    : idifflm,ORG,Blog,BlogM,consF,FltRans,
     &                            ical_prp,ical_flmlt
!
      use module_metrix,only    : WDRBND
!
      
!
      implicit none
!
!-----------------------
! --- [dummy arguments] 
!-----------------------
      integer,intent(inout)  :: ierror
      integer,INTENT(INOUT)  :: ipfix(MXMAT),iptfix(MXMAT)
      integer,INTENT(INOUT)  :: itery(0:ncomp),
     &                          iterk(mxrans),
     &                          iterv(3),
     &                          iterp
      real*8,INTENT(INOUT)   :: repsy(0:ncomp),
     &                          repsk(mxrans),
     &                          repsv(3),
     &                          repsp
      real*8,INTENT(INOUT)   :: aepsy(0:ncomp),
     &                          aepsk(mxrans),
     &                          aepsv(3),
     &                          aepsp
      real*8,INTENT(INOUT)   :: erry(0:ncomp),
     &                          errk(mxrans), 
     &                          errv(3),
     &                          errp
      real*8,INTENT(INOUT)   :: errsdy(0:ncomp),
     &                          errsdk(mxrans),
     &                          errsdv(3)
      real*8, intent(inout)  :: reps_FAI(2),
     &                          aeps_FAI(2),
     &                          err_FAI(2),
     &                          reps_A(3,2),
     &                          aeps_A(3,2),
     &                          err_A(3,2),
     &                          errsdA(3,2),errsdFAI(2)
      integer,INTENT(INOUT)  :: iter_FAI(2),iter_A(3,2)
!      
!----------------------
! --- Euler two phase
!----------------------
!
!      integer,INTENT(INOUT)  :: itery2(0:ncomp,NPHS),
!     &                          iterv2(3,NPHS)
!      real*8,INTENT(INOUT)   :: repsy2(0:ncomp,NPHS), 
!     &                          repsv2(3,NPHS)
!      real*8,INTENT(INOUT)   :: aepsy2(0:ncomp,NPHS) , 
!     &                          aepsv2(3,NPHS)
!      real*8,INTENT(INOUT)   :: erry2(0:ncomp,NPHS)  ,  
!     &                          errv2(3,NPHS)
!      real*8,INTENT(INOUT)   :: errsdy2(0:ncomp,NPHS),
!     &                          errsdv2(3,NPHS)
!
      integer,INTENT(INOUT)  :: itery2(0:ncomp),
     &                          iterv2(3)
      real*8,INTENT(INOUT)   :: repsy2(0:ncomp), 
     &                          repsv2(3)
      real*8,INTENT(INOUT)   :: aepsy2(0:ncomp) , 
     &                          aepsv2(3)
      real*8,INTENT(INOUT)   :: erry2(0:ncomp)  ,  
     &                          errv2(3)
      real*8,INTENT(INOUT)   :: errsdy2(0:ncomp),
     &                          errsdv2(3)
!
      real*8 ,intent(inout)  :: err_POTN(MXPOTN)
      real*8 ,intent(inout)  :: reps_POTN(MXPOTN),aeps_POTN(MXPOTN)
      integer,intent(inout)  :: iter_POTN(MXPOTN)
!---------------
! --- 
!---------------
      INTEGER,INTENT(INOUT)   :: LVEDGE(      2,MXCVFAC)
      INTEGER,INTENT(INOUT)   :: LVRTCV(        MXALLCV)
      INTEGER,INTENT(INOUT)   :: LFUTAU(        MXCV)
      INTEGER,INTENT(INOUT)   :: MAT_NO(        0:MXMAT)
      INTEGER,INTENT(INOUT)   :: MAT_CV(        MXALLCV)
      INTEGER,INTENT(INOUT)   :: LCV_CV(        MXALLCV)
      INTEGER,INTENT(INOUT)   :: MAT_INDEX(     0:MXMAT)
      INTEGER,INTENT(INOUT)   :: MAT_CVEXT(     0:MXMAT)
      INTEGER,INTENT(INOUT)   :: MAT_DCIDX(     0:MXMAT)
      INTEGER,INTENT(INOUT)   :: MAT_CFIDX(     0:MXMAT)
      INTEGER,INTENT(INOUT)   :: Pstart(          MXMAT)
      logical,INTENT(INOUT)   :: mat_cal(       0:MXMAT)
!
      INTEGER,INTENT(INOUT)   :: LCYCSF(        MXSSFBC)
      INTEGER,INTENT(INOUT)   :: LBC_SSF(       MXSSFBC)
!
      REAL*8 ,INTENT(INOUT)   :: SFAREA(      4,MXCVFAC)
      REAL*8 ,INTENT(INOUT)   :: SFCENT(      3,MXCVFAC)
      REAL*8 ,INTENT(INOUT)   :: CVVOLM(        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: CVVOL0(        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: CVCENT(      3,MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: WIFACE(        MXCVFAC)
!
      REAL*8 ,INTENT(INOUT)   :: XTA   (        MXCVFAC)
      REAL*8 ,INTENT(INOUT)   :: FRSTCV(        MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: DISALL(        MXCV)
      REAL*8 ,INTENT(INOUT)   :: DISINL(        MXSSFBC)
!
      REAL*8 ,INTENT(INOUT)   :: RVA   (        MXCVFAC,2)
      REAL*8 ,INTENT(INOUT)   :: RHO   (        MXALLCV,2)
      REAL*8 ,INTENT(INOUT)   :: PRS   (        MXALLCV,2)
      REAL*8 ,INTENT(INOUT)   :: PP0   (        MXALLCV,2)
      REAL*8 ,INTENT(INOUT)   :: TMP   (        MXALLCV,2)
      REAL*8 ,INTENT(INOUT)   :: YYS   (        MXALLCV,MXCOMP,2)
      REAL*8 ,INTENT(INOUT)   :: VEL   (        MXALLCV,3,2)
      REAL*8 ,INTENT(INOUT)   :: VELF  (        MXCVFAC_B,3,2)
      REAL*8 ,INTENT(INOUT)   :: AKS   (        MXALLCVR,MXRANS,2)
      REAL*8 ,INTENT(INOUT)   :: RMUT  (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: RMX   (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: RMU   (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: RMD   (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: RDS   (        MXALLCV,MXCOMP)
      REAL*8 ,INTENT(INOUT)   :: WDOT  (        MXALLCV,MXCOMP)
      REAL*8 ,INTENT(INOUT)   :: CCC   (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: hhh   (        MXALLCV,2)
      REAL*8 ,INTENT(INOUT)   :: hhs   (        MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT)   :: cps   (        MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT)   :: cp    (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: cr    (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: PRSBND(        MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: TMPBND(        MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: YYSBND(        MXSSFBC,MXCOMP)
      REAL*8 ,INTENT(INOUT)   :: VELBND(        MXSSFBC,3)
      REAL*8 ,INTENT(INOUT)   :: AKSBND(        MXSSFBCR,MXRANS)
      REAL*8 ,INTENT(INOUT)   :: HTCBND(        MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: HTFLUX(        MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: MTCBND(        MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: RADBND(        MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: SUFBND(        MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(INOUT)   :: SDOT  (        MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(INOUT)   :: MOLFRC(      MXSSFBC_SUF,MXCOMPALL,2)
      REAL*8 ,INTENT(INOUT)   :: MOLCEN(        MXSSFBC_SUF,MXCOMPALL)
      REAL*8 ,INTENT(INOUT)   :: SITDEN(        MXSSFBC_SUF,MXPHASE,2)
      REAL*8 ,INTENT(INOUT)   :: BLKTHK(        MXSSFBC_SUF,MXPHASE,2)

      REAL*8 ,INTENT(INOUT)   :: UTAU  (        0:MXSSFBC)
      REAL*8 ,INTENT(INOUT)   :: GRDC  (        MXALLCV,3,3)
      REAL*8 ,INTENT(INOUT)   :: DIAG  (        MXALLCV)
      REAL*8 ,INTENT(INOUT)   :: VLASUM(        MXCV)
      REAL*8 ,INTENT(INOUT)   :: dvdd  (        MXCV, 3,2)
      REAL*8 ,INTENT(INOUT)   :: rvx   (        MXALLCV,3)
      REAL*8 ,INTENT(INOUT)   :: rvd   (        MXCVFAC,3)
      INTEGER,INTENT(INOUT)   :: kdbt  (        MXCVFAC)
      INTEGER,INTENT(INOUT)   :: kdby  (        MXCVFAC)
!
! --- Euler 2 Phase (if ieul2ph=1)
!
!      REAL*8,INTENT(INOUT)  :: VEL2     (       MXALLCV2,3,2,NPHS)
!      REAL*8,INTENT(INOUT)  :: tmp2     (       MXALLCV2,2,NPHS)
!      REAL*8,INTENT(INOUT)  :: RVA2     (       MXCVFAC2,2,NPHS)
!      REAL*8,INTENT(INOUT)  :: RHO2     (       MXALLCVC,2,NPHS)
!      REAL*8,INTENT(INOUT)  :: RMUT2    (       MXALLCV2  ,NPHS)
!      REAL*8,INTENT(INOUT)  :: RMU2     (       MXALLCV2  ,NPHS)
!      REAL*8,INTENT(INOUT)  :: YYS2     (       MXALLCV2,MXCOMP,2,NPHS)
!      REAL*8,INTENT(INOUT)  :: RDS2     (       MXALLCV2,MXCOMP,NPHS)
!      REAL*8,INTENT(INOUT)  :: WDOT2    (       MXALLCV2,MXCOMP,NPHS)
!      REAL*8,INTENT(INOUT)  :: RMD2     (       MXALLCV2,NPHS)
!      REAL*8,INTENT(INOUT)  :: ccc2     (       MXALLCV2,NPHS)
!      REAL*8,INTENT(INOUT)  :: dvdd2    (       MXCV2,3,2,NPHS)
!
!      REAL*8,INTENT(INOUT)  :: UTAU2    (     0:MXSSFBC2,NPHS)
!      REAL*8,INTENT(INOUT)  :: VELBND2  (       MXSSFBC2,3,NPHS)
!      REAL*8,INTENT(INOUT)  :: TMPBND2  (       MXSSFBC2,NPHS)
!      REAL*8,INTENT(INOUT)  :: YYSBND2  (       MXSSFBC2,MXCOMP,NPHS)
!      REAL*8,INTENT(INOUT)  :: HTCBND2  (       MXSSFBC2,NPHS)
!      REAL*8,INTENT(INOUT)  :: MTCBND2  (       MXSSFBC2,NPHS)
!      REAL*8,INTENT(INOUT)  :: RADBND2  (       MXSSFBC2,NPHS)
!
      REAL*8,INTENT(INOUT)  :: VEL2     (       MXALLCV2,3,2)
      REAL*8,INTENT(INOUT)  :: PRS2     (       MXALLCV2,2)
      REAL*8,INTENT(INOUT)  :: tmp2     (       MXALLCV2,2)
      REAL*8,INTENT(INOUT)  :: RVA2     (       MXCVFAC2,2)
      REAL*8,INTENT(INOUT)  :: RHO2     (       MXALLCVC,2)
      REAL*8,INTENT(INOUT)  :: RMUT2    (       MXALLCV2  )
      REAL*8,INTENT(INOUT)  :: RMU2     (       MXALLCV2  )
      REAL*8,INTENT(INOUT)  :: YYS2     (       MXALLCV2,MXCOMP,2)
      REAL*8,INTENT(INOUT)  :: RDS2     (       MXALLCV2,MXCOMP)
      REAL*8,INTENT(INOUT)  :: WDOT2    (       MXALLCV2,MXCOMP)
      REAL*8,INTENT(INOUT)  :: RMD2     (       MXALLCV2)
      REAL*8,INTENT(INOUT)  :: ccc2     (       MXALLCV2)
      REAL*8,INTENT(INOUT)  :: hhh2     (       MXALLCV2,2)
      REAL*8,INTENT(INOUT)  :: dvdd2    (       MXCV2,3,2)
!
      REAL*8,INTENT(INOUT)  :: UTAU2    (     0:MXSSFBC2)
      REAL*8,INTENT(INOUT)  :: VELBND2  (       MXSSFBC2,3)
      REAL*8,INTENT(INOUT)  :: TMPBND2  (       MXSSFBC2)
      REAL*8,INTENT(INOUT)  :: YYSBND2  (       MXSSFBC2,MXCOMP)
      REAL*8,INTENT(INOUT)  :: HTCBND2  (       MXSSFBC2)
      REAL*8,INTENT(INOUT)  :: MTCBND2  (       MXSSFBC2)
      REAL*8,INTENT(INOUT)  :: RADBND2  (       MXSSFBC2)
!
      REAL*8,INTENT(INOUT)  :: HEATL_C  (       MXALLCV2)
      REAL*8,INTENT(INOUT)  :: HEATG_C  (       MXALLCV2)     
      REAL*8,INTENT(INOUT)  :: HLS      (       MXALLCV2)
      REAL*8,INTENT(INOUT)  :: HGS      (       MXALLCV2)
      REAL*8,INTENT(INOUT)  :: TMP_SAT  (       MXALLCV2)
      REAL*8,INTENT(INOUT)  :: FRIC_C   (       MXALLCV2)
      REAL*8,INTENT(INOUT)  :: MASS_I   (       MXALLCV2)
!
      integer,intent(inout) :: LCYCOLD  (       MXSSFBC_SLD)
      integer,intent(in)    :: locmsh   (       MXSSFBC)
      real*8 ,intent(inout) :: wifsld   (       MXSSFBC_SLD)
      real*8 ,intent(inout) :: OPPANG   (       MXSSFBC_SLD)
!      integer,intent(in)    :: vctr     (       MXCV_V,0:MXBND_V)
      integer,intent(in)    :: vctrp    (       MXCV_VP,0:MXBND_VP)
      real*8 ,intent(inout) :: FIELD_U  (       MXCV_D,NFLID)
!
      REAL*8,INTENT(INOUT)  :: POTNAL   (       MXALLCVP,MXPOTN,2)
      REAL*8,INTENT(INOUT)  :: POFLUX   (       MXCVFACP,MXPOTN)
      REAL*8,INTENT(INOUT)  :: POTNBC   (       MXSSFBCP,MXPOTN)
!
      integer,INTENT(INOUT) :: NEIGHBOR(MEP,MXALLCV_P)	 	
!      INTEGER,INTENT(INOUT) :: NVCRSIGN(MEP,MXALLCV_P)  
!      INTEGER,INTENT(INOUT) :: LFCV(0:MEP,MXALLCV_P)    
      real*8,INTENT(INOUT)  :: XYZP(MXPRT,3)
      real*8,INTENT(INOUT)  :: UVWP(MXPRT,3)
      real*8,INTENT(INOUT)  :: PARCEL(MXPRT)
      real*8,INTENT(INOUT)  :: DIA(MXPRT)
      real*8,INTENT(INOUT)  :: PMASS(MXPRT,2)
      real*8,INTENT(INOUT)  :: PYS(MXPRT,MXCOMP,2)
      real*8,INTENT(INOUT)  :: HIS_P(MXPRT)
      real*8,INTENT(INOUT)  :: TMPP(MXPRT,2)
      real*8,INTENT(INOUT)  :: HEIGHT(MEP,MXPRT)
      real*8,INTENT(INOUT)  :: FU(MXPRT,2)
      real*8,INTENT(INOUT)  :: FV(MXPRT,2)
      real*8,INTENT(INOUT)  :: FW(MXPRT,2)
      REAL*8,INTENT(INOUT)  :: RHOP(MXPRT)
      integer,INTENT(INOUT) :: INSIDE(MXPRT)
      integer,INTENT(INOUT) :: JOUT(MXPRT)
      REAL*8,INTENT(INOUT)  :: DYDTP (MXPRT,MXCOMP)
      REAL*8,INTENT(INOUT)  :: DHDTP (MXPRT,2)
      REAL*8,INTENT(IN)     :: DNORM(MXCVFAC_P)
      REAL*8 ,INTENT(INOUT) :: MOMCELL(MXALLCV_P,4)
      REAL*8 ,INTENT(INOUT) :: EVAP_Y(MXALLCV_P,MXCOMP)
      REAL*8 ,INTENT(INOUT) :: EVAP_H(MXALLCV_P,2)
      INTEGER,INTENT(INOUT) :: IPCELL(MXALLCV_P)
      INTEGER,INTENT(INOUT) :: NEIBCELL(MXALLCV_P)
!
      REAL*8 ,INTENT(INout)   :: P_rad      (MXALLCV_RAD,2)
      REAL*8 ,INTENT(INOUT)   :: RadHeatFlux(MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: RadHeat    (MXALLCV_RAD)
!      REAL*8 ,INTENT(INOUT)   :: sumcoef    (MXALLNG)
      REAL*8 ,INTENT(INOUT)   :: PAbsorb    (MXALLCV_RAD)
      REAL*8 ,INTENT(INOUT)   :: PScatter   (MXALLCV_RAD)
!      REAL*8 ,INTENT(INOUT)   :: GWAbsorb   (MXALLNG)
!      REAL*8 ,INTENT(INOUT)   :: radinten  (NDNG,MXALLCV_RAD)
!
! --- [local entities]
!
      integer,save :: ictl=0
      integer :: ismpl,msmpl=0,iphs=1
      integer :: MNVX,iact0,IIMAT,IMAT,ICVS,ICVE,IMAT_SIZE,IMAT_U
      real*8  :: times,deltt,rsdc,p_grdc=1.d0,volctr
      logical :: ffexit,cvgnc,stpcvg,explicit,reaction=.false.,
     &           stopfile_s,stopfile_l,outputfile
      integer :: nb,kd,myid,ipiso,npiso=3
      logical :: inoheat=.false.,lvect=.false.
      real*8  :: errabA(3,2)
!
      real*8  :: DT0
!      integer :: NPALL
      integer :: mph_1=1,mph_2=2
! --- 

!
      integer :: no(12)
      integer :: i,j,k,low,ii,icv,icvmn,ICVL,KMAT,l,IDCS,IDCE
!
      integer :: iterp2
      real*8  :: aepsp2,repsp2,errp2
      logical :: mark_raddo
      real*8  :: diff,avet
!
      real*8  :: gggg=0.d0
!--------------------------
! --- ZERO Initialization
!--------------------------
!
      rva   =0.d0
      rho   =1.d0
      prs   =1.d0
      pp0   =1.d0
      tmp   =1.d0
      yys   =1.d0
      vel   =0.d0
      velf  =0.d0
      aks   =0.d0
      rmut  =0.d0
      rmx   =0.d0
      rmu   =0.d0
      rmd   =0.d0
      rds   =0.d0
      wdot  =0.d0   ![kg/(s*m^3)]
      
      ccc   =1.d0
      hhh   =0.d0
      hhs   =0.d0
      cps   =0.d0
      cp    =0.d0
      cr    =0.d0
!
      prsbnd=1.d0
      tmpbnd=1.d0
      yysbnd=1.d0
      velbnd=0.d0
      aksbnd=0.d0
      htcbnd=0.d0
      HTFLUX=0.d0
      mtcbnd=0.d0
      radbnd=0.d0
!
      SUFBND=0.d0   
      SDOT=0.d0     ! [kg/(s*m^2)] for both site and bulk
      MOLFRC=0.d0   ! [%]
      MOLCEN=0.d0   ! [site/bulk: mol/m^2,   or gas:mol/m^3]
      SITDEN=0.d0   ! [mol/m^2]
      BLKTHK=0.d0   ! [m]
!
      WDRBND=0.d0
!      LCYCOLD=0
!      wifsld=1.d0
!
      UTAU=1.d0
      DIAG=0.d0
      grdc=0.d0
!
      vlasum=0.d0
      dvdd=0.d0
      rvx=0.d0
!
      itery=0
      iterk=0
      iterv=0
      iterp=0
      ipfix=0
      iptfix=0
      repsy=0.d0
      repsk=0.d0
      repsv=0.d0
      repsp=0.d0
      aepsy=0.d0
      aepsk=0.d0
      aepsv=0.d0
      aepsp=0.d0
      erry=0.d0
      errk=0.d0
      errv=0.d0
      errsdy=0.d0
      errsdk=0.d0
      errsdv=0.d0
      rsdc=0.d0
!
      reps_FAI=0.d0
      aeps_FAI=0.d0
      err_FAI=0.d0
      reps_A=0.d0
      aeps_A=0.d0
      err_A=0.d0
      iter_FAI=0;iter_A=0
      errsdA=0.d0;errsdFAI=0.d0
!
      iter_POTN=0
      reps_POTN=0
      aeps_POTN=0
      err_POTN=0.d0
!-------------------
! --- Euler 2 Phase
!-------------------
      VEL2=0.d0
      PRS2=1.d0
      tmp2=0.d0
      RVA2=0.d0
      RHO2=1.d0
      RMUT2=0.d0            !1.254d-2
      RMU2=0.d0
      YYS2=1.d0
      RDS2=0.d0             !2.41d-2
      WDOT2=0.d0
      ccc2=1.d0
      dvdd2=0.d0
      RMD2=0.d0             !2.41d-2
      UTAU2=0.d0
      VELBND2=0.d0
      TMPBND2=0.d0
      YYSBND2=0.d0
      HTCBND2=0.d0
      MTCBND2=0.d0
      RADBND2=0.d0
      HEATL_C=1.d0
      HEATG_C=1.d0
      HLS=4.19d5 
      hgs=2.67d6 
      tmp_sat=1.00d2 
      TMP_SAT=0.d0
      FRIC_C=0.d0
      MASS_I=0.d0
!
      itery2=0
      iterv2=0
!
      repsy2=0.d0
      repsv2=0.d0
!
      aepsy2=0.d0
      aepsv2=0.d0
!
      erry2  =0.d0
      errv2  =0.d0
      errsdy2=0.d0
      errsdv2=0.d0
!
      MHD_A=0.d0
      MHD_FAI=0.d0
      MHD_FAI0=0.d0
      MHD_RVA=0.d0
      SIGMA=0.d0
      MHDbnd=0.d0
      MHD_CRNT=0.d0
      MHD_CRT0=0.d0
      tempMHD=0.d0
!
      POTNAL(:,:,:)=0.d0
      POFLUX(:,:)=0.d0
      POTNBC(:,:)=0.d0
!
      EVAP_Y(:,:)=0.d0
      EVAP_H(:,:)=0.d0
!
      t_dir(:,:)=0.d0
      ANGBND(:)=0.d0
      gggg=0.d0
!
      iterp2=0
      repsp2=0.d0;aepsp2=0.d0;errp2=0.d0
!
      P_RAD(:,:)=0.d0
      RadHeatFlux=0.d0
      RadHeat=0.d0
!
      PMASS(:,:)=0.d0
      DIA=0.d0
      XYZP=0.d0
      UVWP=0.d0
      PARCEL=0.d0
      PYS=0.d0
      TMPP=0.d0
      RHOP=0.d0
      HIS_P=0.d0
      HEIGHT=0.d0
      MOMCELL=0.d0
      EVAP_Y=0.d0
      EVAP_H=0.d0
!
      if(ipiston==1) then
        timee=endpiston
      endif
!
      if(ical_cavi==1) then
      endif
!
      if(MXCV_D==MXALLCV) then
        FIELD_U=0.d0
      endif
!
!-------------------------------
! --- Reset Logical Block
!-------------------------------
! --- 1) Explicit Euler Method:
!-------------------------------
!
      explicit=
     &         intgvv.eq.eulere.OR.intgvv.eq.Adams_Bashforth.or.
     &         intgvv.eq.Runge_Kutta
      if(explicit.and.idrdp.eq.incomp) then
        if(nsmpl.gt.1) then
          if(my_rank.eq.ROOT) write(ifle,5000) 
 5000     FORMAT(/,2X,4X,3('?'),1X,'WARNING',1X,3('?'),/,6X,
     &    'SIMPLE-iter must be 1 when using explicit time-integral',/)
          if(my_rank.eq.ROOT) write(ifle,*)
     &    '    SIMPLE-iter=1 has been reseted'
          nsmpl=1
        endif
      endif
      if(lovstbc) then
        nsmpl=1
      endif
!
      if(I_PISO==1) then
        npiso=1  !3
      elseif(I_PISO==0) then
        npiso=1
      endif
!
!
      if(explicit) then
        if(intgke==eulere.or.intgty==eulere) then
        if(i_steady==3) then
          call FFRABORT(1,'ERR: flowcon=3 must use implicit')
        endif
        endif
      endif
!
      if(i_steady==3) then
        nsmpl=1
      endif
!
! --------------------------------------------------------------------
! --------------------------------------------------------------------
! --- monitor ICV
! ----------------
      call initvar(MNVX,MAT_CV)
! --------------------------------------------------------------------
! --- < 4. Input data >-
! --------------------------------------------------------------------
      if(ical_t.or.idrdp==mach0.or.
     &  (nneq.ge.1.and.(ical_reac.or.ical_surf))) then
        ical_t=.true.
      elseif(idrdp.eq.comp.and.ical_topt==1) then
        ical_t=.true.
        ical_topt=1
      elseif(idrdp.eq.comp.and.ical_topt==2) then
        ical_t=.false.
        ical_topt=2   ! dannetu jyouken
      endif
!-----------------------------
! --- Volume center factor
!-----------------------------
      if(icalrv==face) then
        volctr=0.d0
      else
        volctr=1.d0
      endif
!
!-----------------------0000
! --- 
!-----------------------
      IF(ical_vect.and.ical_MHD/=0) then
        call FFRABORT(1,'ERR: VECTOR Ver. NOT SUPPORT for MHD')
      endif
!-----------------------
! --- 
!-----------------------
!      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcimin(lclsd(IIMAT))
!        CALL hpcimin(iclosd(IIMAT))
!        enddo
!        DO IIMAT=1,NMAT 
!        CALL hpcimin(lclsd(IIMAT))
!        CALL hpcimin(iclosd(IIMAT))
!        enddo
!      endif
!
      if(ical_vof==1) then
      endif
! ---------------------------------
! --- < 4.1 initial condition >----
! ---------------------------------
      call read_initial(
     & iters,times,
     & CVVOLM,CVCENT,DISALL,
     & MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,mat_cfidx,
     & LBC_SSF,LCYCSF,mat_cal,LVEDGE,wiface,
     & dvdd,pp0(:,1),prs(:,1),tmp(:,1),yys(:,:,1),vel(:,:,1),
     & velf(:,:,1),rva(:,1),
     & aks(:,:,1),rho(:,1),hhh(:,1),
     & MHD_A,MHD_FAI,SIGMA,MHD_CRT0,MHD_CRNT,
     & dvdd2,rho2(:,1),tmp2(:,1),yys2(:,:,1),
     & vel2(:,:,1),PRS2(:,1),rva2(:,1),hhh2(:,1),
     & MOLFRC,SITDEN(:,:,1),BLKTHK(:,:,1),wifsld,LCYCOLD,OPPANG,
     & cord,lacell,lvcell,
     & POTNAL(:,:,1),WDRBND,
     & ierror)
       if(ierror.ne.0) goto 9999
!
! --- 
!
      call dc_symprv
     &   (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &    LCYCOLD,wifsld,OPPANG,
     &    SFAREA,SFCENT,vel(:,:,1),0,1)
!
      if(ieul2ph>0) then
      endif
! -----------------------------------
! --- enthalpy (hhh) initialzation
! -----------------------------------
      if(ical_t.and.ical_prp) then
        if(iters.le.0) then
          if(sw) then
            call cal_t2h(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &          tmp(:,1),yys(:,:,1),hhs,hhh(:,1),rho(:,1),pp0(:,1)) !pp0
            if(ieul2ph>0) then
            endif
          else
            if(ical_vof==1) then
            elseif(ical_cavi==1) then
            else
              call cal_t2hcp(mph_1,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &          tmp(:,1),yys(:,:,1),hhs,cps,hhh(:,1),cp,cr,prs)
              if(ieul2ph>0) then
              endif
            endif
          endif
        else
          if(ical_vof.eq.1) then 
          elseif(ical_cavi==1) then
          else
            if(iheat==1) then
              call cal_h2t(mph_1,NCV,
     &        MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,hhh(:,1),
     &        yys(:,:,1),tmp(:,1),prs(:,1),rho(:,1))
              if(ieul2ph>0) then
              endif
            endif
          endif
        endif
      endif
! -----------------------------------------0000
! --- calculate density and sound speed^2 
! -----------------------------------------
      if(idrdp/=incomp.or.
     &  (idrdp==incomp.and.ical_dens==1).or.
     &  (idrdp==incomp.and.ical_dens==2).or.
     &  (ical_dens==4).or.
     &  (idrdp==incomp.and.ical_dens==5).or.
     &  (idrdp==incomp.and.ical_dens==3)
     &  .or.(ical_cavi==1)
     &  ) then 
        if(calgeq.or.calxi) then 
          if(idifflm==ORG) then  !ORG  idifflm==ORG
            call cal_flamelet
     &          (iter,MAT_CVEXT,MAT_DCIDX,MAT_CAL,MAT_NO,aks,rho,tmp,YYS)
          else
          endif
        else
          call bc_setbnd(mph_1,iter,time,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,DISINL,SFAREA,SFCENT,locmsh,
     & tmp(:,1),rho(:,1),vel(:,:,1),yys(:,:,1),
     &  aks(:,:,1),prs(:,1),XTA,ccc,
     &  prsbnd,aksbnd,
     &  velbnd,tmpbnd,yysbnd,htcbnd,mtcbnd,ANGBND,SHUTFL)
!
          call cal_rhoccc !1)
     & (mph_1,NALLCV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &  DIAG,hhs,cps,cp,cr,
     &  aks(:,:,1),tmp(:,1),yys(:,:,1),pp0(:,1),prs,rho,
     &  rmu,rmd,
     &  RHO2(:,1),VEL(:,:,1),ccc,
     &  yysbnd,prsbnd,tmpbnd,LBC_SSF,lvedge,1)
        endif

        if(ieul2ph>0) then
          call FFRABORT(1,'ERR:Euler 2-Phase NOT support E2P')
        endif
!
      endif
!
      if(iters.lt.1) then
        if(ieul2ph>0) then
        endif
        call cal_rva(
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &  MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,wiface,SFCENT,
     &  LCYCOLD,wifsld,OPPANG,
     &  xta,rho(:,1),vel(:,:,1),rvx,rva(:,1))
      endif
!
      call reset_vars(0,
     &  LVEDGE,LCYCSF,
     &  MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,mat_cal,MAT_DCIDX,
     &  pp0,aks,
     &  rva,yys,vel,hhh,tmp,prs,
     &  rva2,yys2,vel2,hhh2,tmp2,prs2)
!
!---------------------------------
!--< 4.2 source domain >--
!---------------------------------
      call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &               LVEDGE,LBC_SSF,LCYCSF,msk)
!
      if(NPE.eq.1.and..not.ical_vect) then
        call read_source(NCV,ierror)
      endif
      if( ierror.ne.0 ) goto 9999 
!
!-< II. Time marching >-
!
      iter=iters
      time=times
!
 1000 continue
      iter=iter+1
!------------------------------0000
! --- < 2. Reset time level >-
!------------------------------
      if(ical_mvmsh==2.or.
     &   ical_mvmsh==3.or.
     &   ical_mvmsh==4
     &   ) goto 1300
!
!------------------------------------------
! --- < 1. Preliminary set >-
!------------------------------------------
! --- < 1.0 Set material calculation flag>
!------------------------------------------
      reaction=.false.
      if((nneq>=1.or.nneqP>=1).and.(ical_reac.or.ical_surf.or.ical_P_R==1)
     &  .and.iter.gt.igT_iter) then
        reaction=.true.
      endif
      call mat_flag(iter,Pstart,mat_cal)
!------------------------------------------
! --- < 1.1 transport properties >--      
!     set_trnsprp(rmx) ==> set_deltat     
!     set_trnsprp(fric_c) ==> set_deltat2 
!------------------------------------------
      if(ical_vof==1) then
      elseif(ical_cavi==1) then
      else
        if(ical_prp) then  !ORG idifflm/=ORG
          call set_trnsprp(mph_1,
     &    LVEDGE,LBC_SSF,LCYCSF,
     &    MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,wifsld,LCYCOLD,
     &    DIAG,hhs,cps,cp,cr,
     &    tmp(:,1),yys(:,:,1),aks(:,:,1),prs(:,1),rho(:,1),
     &    rmu,rmd,rds,rmx)
          if(ieul2ph>0) then
          endif
        endif
      endif
!----------------------------------
! --- < 1.2 turbulent viscosity >--
!----------------------------------
      if(.not.explicit.and..not.(ideltt.eq.const)) then
        call cal_eddy_viscosity(iter,deltt,time,
     &    LVEDGE,LBC_SSF,LCYCSF,LFUTAU,
     &    MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,MAT_INDEX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,DISALL,FRSTCV,GRDC,OPPANG,
     &    rho(:,1),vel(:,:,1),rmu,rmut,rmut2,aks(:,:,1),tmp(:,1),utau,
     &    wifsld,LCYCOLD,vctr,FIELD_U,
     &    yplusf,vflc)
      endif
!-----------------------------
! --- < 1.3 time increment >--
!------------------------------------------------------------
! --- two cases are set: constant delt-t, and auto-delt-t 
!------------------------------------------------------------
!
 1300 continue
      call set_deltat(mph_1,iter,
     &  LVEDGE,LBC_SSF,LCYCSF,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  CVVOLM,CVCENT,SFAREA,
     &  cp,cr,wifsld,LCYCOLD,vctr,ccc,wiface,
     &  rva,rho(:,1),rmx,rmut,deltt)
      if(ieul2ph>0) then
      endif
!-------------------------------0000
! --- < 2. Reset time level >-
!-------------------------------
      time=time+deltt
!---------------------------------------------
! --- < 2.1 metrics in case of moving grid >--
!---------------------------------------------
      CVVOL0=CVVOLM
!
      call metric_admin(2,iter,time,
     &  LVEDGE,LVRTCV,LFUTAU,
     &  LBC_SSF,LCYCSF,FRSTCV,DISINL,
     &  SFAREA,SFCENT,CVVOLM,CVVOL0,CVCENT,wiface,DISALL,xta,
     &  MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LCYCOLD,wifsld,OPPANG,
     &  movflg,cord,lacell,lvcell,
     &  lvface,lbface,lcface,lfcell,LEFACE,listbc,locmsh,
     6  area,volume,gface,gcell,
     &  deltt,ierror)
      if(ierror.ne.0 ) goto 9999
      if(ical_mvmsh==2.or.
     &    ical_mvmsh==3.or.
     &    ical_mvmsh==4
     &   ) goto 1400

      if(ical_sld/=0) 
     &  call CONTROL_RPM(time,deltt,iter,MAT_NO)
!---------------------------------------------
! --- Define injector CV number 
!---------------------------------------------
      if(ictl==0) then
        call cal_injcoord(deltt,time,iter,ictl,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  CVVOLM,CVCENT)
      endif
!-----------------------------------
! --- < 2.2 dependent variables >--
!-----------------------------------
      tmp(  :,2)=tmp(  :,1)
      hhh(  :,2)=hhh(  :,1)
      yys(:,:,2)=yys(:,:,1)
      vel(:,:,2)=vel(:,:,1)
      velf(:,:,2)=velf(:,:,1)
      rho(  :,2)=rho(  :,1)
      prs(  :,2)=prs(  :,1)
      pp0(  :,2)=pp0(  :,1)
      if(ieul2ph==0) then
        rva(  :,2)=rva(  :,1)
        dvdd(:,:,2)=dvdd(:,:,1)
      endif
      if(ical_FC==PEFC) then
        dvdd(:,:,1)=0.d0
      endif
      aks(:,:,2)=aks(:,:,1)
!
      wdot(:,:)=0.d0
!      EVAP_Y(:,:)=0.d0
!      EVAP_H(:,:)=0.d0
!
      MOLFRC(:,:,2)=MOLFRC(:,:,1)
      SITDEN(:,:,2)=SITDEN(:,:,1)
      BLKTHK(:,:,2)=BLKTHK(:,:,1)
!--------------
! --- ieul2ph
!--------------
      if(ieul2ph>0) then
      endif
      if(ical_MHD/=0) then
        MHD_A(:,:,1,2)=MHD_A(:,:,1,1)
        MHD_A(:,:,2,2)=MHD_A(:,:,2,1)
        MHD_FAI(:,1,2)=MHD_FAI(:,1,1)
        MHD_FAI(:,2,2)=MHD_FAI(:,2,1)
      endif
      if(ical_FC==PEFC) then
        POTNAL(:,:,2)=POTNAL(:,:,1)
      endif
!
      if(MXCV_D==MXALLCV) then
        FIELD_U=0.d0
      endif
!
      call src_setnow(time)
      call bc_setnow(time)
!---------------------0000
! --- single phase BC
!---------------------(IBFL,1)
      call bc_setbnd(mph_1,iter,time,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,DISINL,SFAREA,SFCENT,locmsh,
     &  tmp(:,1),rho(:,1),vel(:,:,1),yys(:,:,1),aks(:,:,1),
     &  prs(:,1),XTA,ccc,
     &  prsbnd,aksbnd,
     &  velbnd,tmpbnd,yysbnd,htcbnd,mtcbnd,ANGBND,SHUTFL)
!
      if(ieul2ph>0) then
      endif
!
! --- pressure BC
!
      incmp_prs(:)=1.d10
      do nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      kd=kdbcnd(0,nb)
      if(kd==kdpres.or.kd==kdolet) then
        if(incmp_prs(IIMAT)>=wdbcnd(1,nb)) then
          incmp_prs(IIMAT)=wdbcnd(1,nb)
        endif
      endif
      enddo
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        call hpcrmin(incmp_prs(IIMAT))
!        incmp_prs(0)=1.d10
!        enddo
        DO IIMAT=1,NMAT
        call hpcrmin(incmp_prs(IIMAT))
!        incmp_prs(0)=1.d10
        enddo
        
      endif
!
      if(idrdp.eq.mach0) then
        call bc_pp0
     &  (LVEDGE,mat_cal,
     &   MAT_CV,MAT_CVEXT,MAT_NO,LBC_SSF,
     &   prsbnd,pp0(:,1))
      endif
!
      call bc_prs
     & (nflud,ipfix,deltt,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,FRSTCV,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,LBC_SSF,mat_cal,
     &  LVEDGE,LCYCSF,vel,prsbnd,tmpbnd,rho(:,1),ccc,yys(:,:,1),
     &  wifsld,LCYCOLD,
     &  prs(:,1),pp0(:,1),tmp(:,1),rva(:,1))
!
      if(encomp>=1) then
        call bc_tys
     &    (mph_1,LVEDGE,LCYCSF,LBC_SSF,mat_cal,MAT_DCIDX,MAT_NO,
     &     LCYCOLD,wifsld,
     &     rva(:,1),tmpbnd,yysbnd,vel,rho(:,1),ccc,
     &     tmp(:,1),yys(:,:,1),hhh(:,1),pp0(:,1))
      endif
!
      if(idrdp/=incomp.or.
     &  (idrdp==incomp.and.ical_week).or.
     &  (idrdp==incomp.and.ical_dens==1).or.
     &  (idrdp==incomp.and.ical_dens==2).or.
     &  (ical_dens==4).or.
     &  (idrdp==incomp.and.ical_dens==5).or.
     &  (idrdp==incomp.and.ical_dens==3)
     &  .or.(ical_cavi==1)
     &   ) then
        if(calgeq.or.calxi) then  !OK

           if(idifflm==ORG) then  !ORG  idifflm==ORG
             call cal_flamelet
     &          (iter,MAT_CVEXT,MAT_DCIDX,MAT_CAL,MAT_NO,aks,rho,tmp,YYS)
           else
           endif
        else
          call cal_rhoccc !2)
     &    (mph_1,NALLCV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &     DIAG,hhs,cps,cp,cr,
     &   aks(:,:,1),tmp(:,1),yys(:,:,1),pp0(:,1),prs,rho,
     &   rmu,rmd,
     &   RHO2(:,1),VEL(:,:,1),ccc,
     &            yysbnd,prsbnd,tmpbnd,LBC_SSF,lvedge,2)
        endif
      endif
!
      if(imasflg.or.ifixmach) then
        call bc_setbnd(mph_1,iter,time,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,mat_cal,DISINL,SFAREA,SFCENT,locmsh,
     & tmp(:,1),rho(:,1),vel(:,:,1),yys(:,:,1),aks(:,:,1),prs(:,1),XTA,ccc,
     &  prsbnd,aksbnd,
     &  velbnd,tmpbnd,yysbnd,htcbnd,mtcbnd,ANGBND,SHUTFL)
      endif
!
      call bc_vel
     & (mph_1,deltt,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,SFCENT,
     &   LCYCOLD,wifsld,OPPANG,locmsh,
     &   xta,rva(:,1),rho(:,1),velbnd,prsbnd,prs(:,1),ccc,aksbnd,
     &   vel,tmp(:,1),yys(:,:,1),10)
!
      call bc_rva   
     & (mph_1,time,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &  LCYCOLD,wifsld,OPPANG,
     &  CVCENT,SFCENT,wiface,locmsh,aks(:,:,1),
     &  SFAREA,rho(:,1),vel,xta,velbnd,prsbnd,aksbnd,
     &  prs(:,1),rva(:,1),1)
      if(ieul2ph>0) then
      endif
!------------------
! --- two phase BC
!------------------
      if(ieul2ph>0) then
      endif
!-------------------------
!-< 3. Time integration >-
!-------------------------
! --- In every time-step, 'corrector number':ismpl=1,nsmpl
! --- If we are computing an unsteady flow or transient flow, 
!     and time accuracy is required, correcting iterations must be 
!     continued until the entire system of all non-linear eqs. are
!     satisfied to within a narrow tolerance.
!
      msmpl=0
      do 1100 ismpl=1,nsmpl

!----------------------------------------------------------------
!--< 3.2 turbulent viscosity >-- only for implicit time-integral
!----------------------------------------------------------------
      if(ical_prt==2) goto 1600 
!
      IF(ical_MHD/=2) then
        call cal_eddy_viscosity(iter,deltt,time,
     &   LVEDGE,LBC_SSF,LCYCSF,LFUTAU,
     &   MAT_NO,MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,MAT_INDEX,
     &   SFAREA,SFCENT,wiface,CVCENT,CVVOLM,DISALL,FRSTCV,GRDC,OPPANG,
     &   rho(:,1),vel(:,:,1),rmu,rmut,rmut2,aks(:,:,1),tmp(:,1),utau,
     &   wifsld,LCYCOLD,vctr,FIELD_U,
     &   yplusf,vflc)
      endif
!      
!----------------------------------------------------------------
! --- 
!----------------------------------------------------------------
      if(pot_scl.and.iter>iterPOTEN) then 
        call potential_admin(deltt,iter,iterPOTEN,time,ismpl,nsmpl,
     &   ipfix,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,FRSTCV,disall,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  wifsld,LCYCOLD,FIELD_U,
     &  POTNAL,POFLUX,POTNBC,
     &  kdbt,RMX,diag,rvx(:,1),rho(:,1),rho2(:,1),iptfix,aks(:,:,1),
     &  grdc(:,:,1),vctr,rva(:,1),
     &  yys(:,:,1),tmp(:,1),prs(:,1),OPPANG,
     &  iter_POTN,reps_POTN,aeps_POTN,err_POTN
     &       )
      endif
!
!------------------------------------------
!--< 3.6 RANS model & scalar transport >--
!------------------------------------------
!
      if(rns_scl) then
        if(ieul2ph>0) then
        endif
!
        if(encomp/=2) then
        call rans_admin(
     &  deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,FRSTCV,disall,xta,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  rva,rmu,rmut,rho(:,1),tmp,vel(:,:,1),prs(:,1),!????
     &  yys(:,:,1),dvdd(:,:,1),EVAP_Y,
     &  rva2,rmu2,rmut2,rho2(:,1),tmp2(:,1),vel2(:,:,1),
     &  yys2(:,:,1),
     &  aksbnd,aks,MASS_I,
     &  grdc,diag,kdbt,rvx,
     &  rvd,RMX,cp,cr,RMD,
!     &  LCYCOLD,wifsld,OPPANG,vctr,FIELD_U,POTNAL(:,:,1),
     &  LCYCOLD,wifsld,OPPANG,FIELD_U,POTNAL(:,:,1),
     &  pp0(:,1),vflc,yplusf,rva_s,utau,LFUTAU,VELBND,
     &  iterk,repsk,aepsk,errk,ierror)
        ENDIF

        endif
        if( ierror.ne.0 ) goto 9999
!----------------------------------
 1600 continue
!----------------------------------
      if(ical_MHD==1.or.ical_MHD==2) then
      endif
!----------------------------------
!--< 3.3-0 Chemical Reaction rate
!----------------------------------
!
      if(ical_prt==2) goto 1601
!
      if(reaction) then   !if(reaction.and..not.steady) then
        if(ical_FC==PEFC) then
        else
          wdot(:,:)=0.d0
        endif
        SDOT(:,:)=0.d0
        call chem_admin_unstdy(iter,deltt,
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,MAT_CFIDX,MAT_CV,
     &  LBC_SSF,LCYCSF,LVEDGE,wiface,wifsld,LCYCOLD,OPPANG,vctr,
     &  SFAREA,SFCENT,CVCENT,CVVOLM,vel(:,:,1),
     &  prs(:,1),rho(:,1),pp0(:,1),aks(:,:,1),tmp(:,1),yys(:,:,1),wdot,
     &  SDOT,MOLFRC,MOLCEN,SITDEN(:,:,1),
     &  hhs,cp,cr,FIELD_U,rvx(:,1),grdc,
     &  ierror)
        if( ierror.ne.0 ) goto 9999
      endif
!
!      if(NFLID>=1) then
!        FIELD_U(1:NCV,1)=my_rank+1
!      endif
!--------------------------------
! --- Suface reaction mechanism
!--------------------------------
!      
      if(reaction.and.ical_suf==1) then
      endif
!
!------------------------------------------------------------
!--< 3.3-1 Radiation heat transfer by FVM/DOM/DTM method
!		--- Developed by Jiang Yuyan, 2006/4.�`5.
!------------------------------------------------------------
!	
      if(radflag==1.and.radmodel(1:3).eq.'FVM') then
        if(iter-iters.ge.2.and.mod(iter-iters,radjump).eq.0) then
	   mark_raddo=.true.
	   if(iter-iters.GT.radjump) then
             diff=0.d0
             avet=0.d0
             do ICV=1,NCV
!here, sumcoef is old value of tmp
             diff=diff+(tmp(ICV,1)-sumcoef(ICV))		
     &		 *(tmp(ICV,1)-sumcoef(ICV))
             avet=avet+tmp(ICV,1)
             
             enddo
             diff=sqrt(diff)
             if(avet.LE.1.0e-2) then
               mark_raddo=.false.
             else if(diff/avet.lt.1.0e-4) then
               mark_raddo=.false.
             endif
	   endif
           mark_raddo=.true.
	   if(NPE.GT.1) call hpclor(mark_raddo)
	   if(mark_raddo) then
	   call radfvm_admin 
     &  (ismpl,iphs,deltt,iter,p_grdc,time,ictl,
     &	LVEDGE,LCYCSF,LBC_SSF,
     &	SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &	rva(:,1),rho,vel(:,:,1),prs,pp0,
     &	kdbt,vctr,cps,cp,diag,
     &	MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &	LCYCOLD,wifsld,
     &  P_rad,RadHeatFlux,RadHeat,sumcoef,PAbsorb,PScatter,GWAbsorb,
     &  radinten,
     &	tmpbnd,yysbnd,htcbnd,mtcbnd,radbnd,tmp,yys,
     &	iterR,repsR,aepsR,errR,ierror)
           endif
        endif
      endif
!-----------------------------------------------------
! --- PARTICLE TRACING 
!-----------------------------------------------------
!
 1601 continue
!
      if(ical_prt>=1.and.ismpl==1) then 
        if(ical_prt==2.and.ical_t) then
          call cal_h2t(1,NCV,
     &    MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
     &    hhh(:,1),yys(:,:,1),tmp(:,1),prs(:,1),rho(:,1))
        endif

        DT0=deltt
        EVAP_Y(:,:)=0.d0
        EVAP_H(:,:)=0.d0

        CALL FFR_TRACING(ITER,ITERS,DT0,ffexit,time,reaction,
     1           MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_CFIDX,MAT_CAL,
     2           LVEDGE,SFCENT,CVVOLM,WIFACE,CVCENT,SFAREA,
     3           GRDC,RHO,RMU,VEL,TMP(:,1),YYS(:,:,1),WDOT,
     4           DNORM,NEIGHBOR,
     5           INSIDE,JOUT,PMASS,DIA,XYZP,UVWP,PARCEL,PYS,TMPP,RHOP,
     &           DYDTP,DHDTP,HIS_P,
     6           HEIGHT,FU,FV,FW,IERROR,
     7           MOMCELL,EVAP_Y,EVAP_H,LBC_SSF,vctrp,
     &           P_rad,radinten,
     8           IPCELL,NEIBCELL,aks,cps,cp,
     7           pp0,rds,rmd,hhs,cr,FIELD_U,WDRBND
     &           ) 
        if(iter==iters+1) wdot(:,:)=0.d0
        if(ical_prt==2) goto 1602
      endif 
!---------------------------------------
!--< 3.3 temperature & mass fraction >--
!---------------------------------------
!--------------------
! --- Euler 2 phase  
!--------------------
      if(ieul2ph>0.and.encomp>=1) then
      endif
!----------------------------------------------------
! --- Single Phase or first phase of Euler 2 phase
!----------------------------------------------------
      if(encomp>=1) then 
        iphs=mph_1
        call hys_admin(ismpl,iphs,deltt,iter,p_grdc,time,ictl,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,
     &  rva(:,1),rmd,rds,rmut,rho,vel(:,:,1),prs,pp0,wdot,rmu,
     &  utau,frstcv,
     &  EVAP_Y,EVAP_H,
     &  hhh,hhs,cps,cp,cr,rvx(:,1),rvx(:,2),rvx(:,3),
     &  grdc,diag,rmx,kdbt,kdby,rvd,
!!!!!!!!!     &  aks(:,:,1),vctr,dvdd(:,:,1),ccc,
     &  aks(:,:,1),dvdd(:,:,1),ccc,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  SUFBND,SDOT,LCYCOLD,wifsld,OPPANG,FIELD_U,
     &  tmpbnd,yysbnd,htcbnd,mtcbnd,radbnd,HTFLUX,
     &  tmp,yys,MASS_I,
     &  RadHeat,RadHeatFlux,pabsorb,pscatter,
     &  itery,repsy,aepsy,erry,ierror)
!
        call bc_tys
     &    (mph_1,LVEDGE,LCYCSF,LBC_SSF,mat_cal,MAT_DCIDX,MAT_NO,
     &     LCYCOLD,wifsld,
     &     rva(:,1),tmpbnd,yysbnd,vel,rho(:,1),ccc,
     &     tmp(:,1),yys(:,:,1),hhh(:,1),pp0(:,1))
      endif
!
      if(idrdp/=incomp
     &  .or.(ical_cavi==1)
     &  ) then
        if(ieul2ph>0) then
          call FFRABORT(1,'ERR:E2P NOT support Zero and comp')
        endif
        if(calgeq.or.calxi) then
          if(idifflm==ORG) then  !ORG  idifflm==ORG
             call cal_flamelet
     &          (iter,MAT_CVEXT,MAT_DCIDX,MAT_CAL,MAT_NO,aks,rho,tmp,YYS)
          else
          endif
        else
          call cal_rhoccc !3)
     &      (mph_1,NALLCV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &      DIAG,hhs,cps,cp,cr,
     &      aks(:,:,1),tmp(:,1),yys(:,:,1),pp0(:,1),prs,rho,
     &      rmu,rmd,
     &      RHO2(:,1),VEL(:,:,1),ccc,
     &      yysbnd,prsbnd,tmpbnd,LBC_SSF,lvedge,3)
        endif
      endif
!------------------------------------------------------------
!--< 3.4 Euler 2 phase flow model. 
!------------------------------------------------------------
! --- set flag for boundary condition: kdbt=>kdbp; kdby=>kdbv
!------------------------------------------------------------VVVV
      call bc_kdbpv(LBC_SSF,LCYCSF,mat_cal,MAT_NO,kdbt,kdby)
!------------------------------------------------------------
! --- Euler 2 phase velocity
!------------------------------------------------------------
      if(ieul2ph>0.and.encomp/=2) then
      endif
!--------------------------
!--< 3.4 velocity field >--
!--------------------------
      IF(ical_MHD/=2.and.encomp/=2) then
        IF(ical_MHD==1) then
          DO I=1,3
          rvx(:,I)=MHD_A(:,I,1,2) 
          enddo
        endif
        iphs=mph_1
        if(ieul2ph>0)
     &   call bc_kdbpv_2p(iphs,LBC_SSF,LCYCSF,mat_cal,kdby)
!
        call vel_admin(
     &  iphs,iter,ismpl,deltt,nsmpl,rsdc,time,ictl,volctr,
     &  LVEDGE,LCYCSF,LFUTAU,LBC_SSF,LCYCOLD,wifsld,OPPANG,
     &  SFAREA,SFCENT,CVCENT,CVVOLM,CVVOL0,wiface,FRSTCV,DISALL,
     &  vel2,rho2,prs2,aks,fric_c,dvdd2,MASS_I,SDOT,xta,pp0(:,1),
     &  rmu,rmut,tmp,rva,vel,rho,prs,yys(:,:,1),utau,dvdd,
     &  grdc,diag,kdby,rvd,velf,
     &  rmx,cp,cr,rvx,bdyfrc,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  FIELD_U,MOMCELL,
     &  angbnd,t_dir,gggg,
     &  iterv,repsv,aepsv,errv,ierror)
!--------------------------------
!--< 3.4 Mass flux in CV face >--
!--------------------------------
        if(ieul2ph>0.and.encomp/=2) then
        endif
!
        if(encomp/=2) then
          iphs=mph_1
          call rva_admin(iphs,iter,ismpl,time,deltt,nsmpl,volctr,
     &    LVEDGE,LCYCSF,LFUTAU,LBC_SSF,
     &    SFAREA,SFCENT,CVCENT,CVVOLM,CVVOL0,wiface,FRSTCV,DISALL,
     &    rva,vel,rho,prs,xta,tmp(:,1),
     &    grdc,rvx,aks(:,:,1),bdyfrc,velf,
     &    MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &    prsbnd,
     &    LCYCOLD,wifsld,OPPANG,locmsh,FIELD_U,
     &    ierror)
        endif
!----------------------------------------------------------
! --- <pressure, correcting pressure, correcting velocity> 
!----------------------------------------------------------
        if(ieul2ph==2.and.encomp/=2) then
        endif
!
        iphs=mph_1
        iterp=0
        do ipiso=1,npiso
        if(ipiso>1) then
!          rvd(:,1)=rva(:,1)   !piso
!          rva(:,1)=rva(:,2)
          if(ieul2ph>0) then
            call FFRABORT(1,'ERR: piso NOT support E2P')
          endif
          call bc_rva  !OK !zh000 
     &    (mph_1,time,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &    LCYCOLD,wifsld,OPPANG,
     &    CVCENT,SFCENT,wiface,locmsh,aks(:,:,1),
     &    SFAREA,rho(:,1),vel,xta,velbnd,prsbnd,aksbnd,
     &    prs(:,1),rva(:,1),2)
        endif
        if(encomp/=2) then

          call prs_admin(iphs,ipiso,
     &    iter,ismpl,deltt,time,nsmpl,ipfix,p_grdc,ictl,volctr,
     &    LVEDGE,LCYCSF,LBC_SSF,
     &    SFAREA,SFCENT,CVCENT,CVVOLM,CVVOL0,wiface,pp0(:,1),
     &    ccc,rva,vel,rho,prs,yys(:,:,1),tmp,velf,
     &    ccc2,rva2,vel2,rho2,prs2,yys2(:,:,1),tmp2,
     &    aks,SDOT,mass_i,
     &    EVAP_Y,
     &    grdc,diag,kdbt,RMX,cr,bdyfrc,rvx,
     &    MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &    LCYCOLD,wifsld,OPPANG,FIELD_U,
     &    iterp,repsp,aepsp,errp,ierror)
        endif
        enddo
!        if(npiso>1) then   !piso
!          rva(:,1)=rvd(:,1)+rva(:,2)
!        endif
!
        if( ierror.ne.0 ) goto 9999 
!-----------------------------------------------
! --- Calculating temperature using enthalpy 
!-----------------------------------------------
        if(idrdp.ne.incomp.and..NOT.ical_flmlt) then 
          call cal_temp(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,yys(:,:,1),
     &                pp0(:,1),prs(:,1),rho(:,1),tmp(:,1),p_grdc)
        endif
!----------------------------
!--< 3.5 reset variables >-- 
!----------------------------
!----------------------------
! --- single phase BC 
!----------------------------
        call bc_prs !OK
     & (nflud,ipfix,deltt,
     &  SFAREA,SFCENT,CVVOLM,CVCENT,FRSTCV,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,LBC_SSF,mat_cal,
     &  LVEDGE,LCYCSF,vel,prsbnd,tmpbnd,rho(:,1),ccc,
     &  yys(:,:,1),wifsld,LCYCOLD,
     &  prs(:,1),pp0(:,1),tmp(:,1),rva(:,1))
!
!        call bc_tys
!     &    (mph_1,LVEDGE,LCYCSF,LBC_SSF,mat_cal,MAT_DCIDX,MAT_NO,
!     &     LCYCOLD,wifsld,
!     &     rva(:,1),tmpbnd,yysbnd,vel,rho(:,1),ccc,tmp(:,1),
!     &     yys(:,:,1),hhh(:,1))  !OK
!
        if(idrdp/=incomp
     &    .or.(ical_cavi==1)
     &    ) then
          if(calgeq.or.calxi) then
!            if(idifflm==ORG) then  !ORG  idifflm==ORG
!              call cal_flamelet
!     &          (iter,MAT_CVEXT,MAT_DCIDX,MAT_CAL,MAT_NO,aks,rho,tmp,YYS)
!            else
!            endif
          else
!            call cal_rhoccc  !4)
!     &       (mph_1,NALLCV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
!     &        DIAG,hhs,cps,cp,cr,
!     &        aks(:,:,1),tmp(:,1),yys(:,:,1),pp0(:,1),prs,rho,
!     &        rmu,rmd,
!     &        RHO2(:,1),VEL(:,:,1),ccc,
!     &            yysbnd,prsbnd,tmpbnd,LBC_SSF,lvedge,4)
          endif
        endif
!
        call bc_vel  !OK
     & (mph_1,deltt,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,SFCENT,
     &  LCYCOLD,wifsld,OPPANG,locmsh,
     &  xta,rva(:,1),rho(:,1),velbnd,prsbnd,prs(:,1),ccc,aksbnd,
     &  vel,tmp(:,1),yys(:,:,1),30)
!
        call bc_rva  !OK !zh000 
     &  (mph_1,time,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &  LCYCOLD,wifsld,OPPANG,
     &  CVCENT,SFCENT,wiface,locmsh,aks(:,:,1),
     &  SFAREA,rho(:,1),vel,xta,velbnd,prsbnd,aksbnd,
     &  prs(:,1),rva(:,1),2)
        if(ieul2ph>0) then
        endif
      endif
!------------------
! --- two phase BC 
!------------------
      if(ieul2ph>0) then 
      endif  !endif(ical_MHD/=2)
!-----------------------------------------
!--< 3.6 RANS model & scalar transport >--
!-----------------------------------------
!
!----------------------------------------------
!--< 3.1 reset variables => solid material  >--
!----------------------------------------------
      call reset_vars(0,
     &  LVEDGE,LCYCSF,
     &  MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,mat_cal,MAT_DCIDX,
     &  pp0,aks,
     &  rva,yys,vel,hhh,tmp,prs,
     &  rva2,yys2,vel2,hhh2,tmp2,prs2)
!----------------------------
!-----------------------------
! --- Simple method iteration
!-----------------------------
!
 1602 continue
      msmpl=msmpl+1

!
      if(.not.steady.and.
     &  (mod(iter,monitor_stp)==0.or.
     &   ffexit.or.outputfile)) then
        call acs_monitor
     &   (msmpl,iter,itery,iterk,iterv,iterp,MNVX,
     &    repsy,aepsy,repsk,aepsk,repsv,aepsv,repsp,aepsp,
     &    erry,errk,errv,errp,errsdy,errsdk,errsdv,
     &    reps_FAI,aeps_FAI,err_FAI,
     &    reps_A,aeps_A,err_A,iter_FAI,iter_A,errsdA,errsdFAI,
     &    iter_POTN,reps_POTN,aeps_POTN,err_POTN,POTNAL(:,:,1),
     &    time,deltt,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &    MHD_A,MHD_FAI,
     &    CVVOLM,tmp(:,1),vel(:,:,1),yys(:,:,1),prs(:,1),aks(:,:,1),
     &    ccc,rho(:,1))
        if(ieul2ph>0) then
        endif
      endif
!
!----------------------------------------
! --- SIMPLE Method Converging criteria:
!----------------------------------------
!
      call error_smac(iter,iters,
     &                ncomp,mxrans,nrans,erry,errk,errv,errp,
     &                itery,iterk,
     &                err_FAI,err_A,
     &                erry2,errv2,
     &                cvgnc)
!      if(cvgnc.and.(ical_sld==0)) goto 1200
      if(cvgnc) goto 1200
!
 1100 continue
 1200 continue
!-----------------------------------------------------
! --- Time step Converging criteria for steady flow:
!-----------------------------------------------------
      if(steady) then
        call error_stp(deltt,
     &           errsdy,errsdk,errsdv,errsdA,errsdFAI,
     &           MAT_CV,MAT_CVEXT,MAT_INDEX,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,vel,aks,MHD_A,MHD_FAI,errabA,
     &           stpcvg)
        if(mod(iter,monitor_stp)==0.or.ffexit.or.outputfile) then 
            call acs_monitor
     &   (msmpl,iter,itery,iterk,iterv,iterp,MNVX,
     &    repsy,aepsy,repsk,aepsk,repsv,aepsv,repsp,aepsp,
     &    erry,errk,errv,errp,errsdy,errsdk,errsdv,
     &    reps_FAI,aeps_FAI,err_FAI,
     &    reps_A,aeps_A,err_A,iter_FAI,iter_A,errsdA,errsdFAI,
     &    iter_POTN,reps_POTN,aeps_POTN,err_POTN,POTNAL(:,:,1),
     &    time,deltt,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &    MHD_A,MHD_FAI,
     &    CVVOLM,tmp(:,1),vel(:,:,1),yys(:,:,1),prs(:,1),aks(:,:,1)
     &   ,ccc,rho(:,1))
          if(ieul2ph>0) then
          endif
        endif
      else
        stpcvg=.false.
      endif
 1400 continue
!
!---------------
!-< 4. Output >-
!---------------
!
      stopfile_l=.false.
      stopfile_s=.false.
      ffexit=.false.
      if((iter.ge.itere).or.(time.ge.timee)) then 
          ffexit=.true.
      elseif(stpcvg) then
          ffexit=.true.
      else
          inquire(file='STOP',exist=stopfile_l)
          inquire(file='stop',exist=stopfile_s)
          if(stopfile_l.or.stopfile_s) then    
            ffexit=.true.
          else
            ffexit=.false.
          endif
      endif
      outputfile=.false.
      inquire(file='output',exist=outputfile)
!
      if(NPE.gt.1) then
        call hpcland(ffexit,1)
        call hpcland(outputfile,1)
      endif
!
!--------------------------
! --- < 4.3 result file >--
!--------------------------
!
      call write_result(
     &      iter,time,deltt,ffexit,outputfile,vlasum,
     &      pp0(:,1),prs(:,1),rho(:,1),tmp(:,1),yys(:,:,1),vel(:,:,1),
     &      aks(:,:,1),rmut,rva(:,1),WDOT,ccc,
     &      MHD_A(:,:,:,1),MHD_FAI(:,:,1),MHD_CRNT,MHD_CRT0,tempMHD,
     &      MHD_A(:,:,:,2),
     &      rho2(:,1),tmp2(:,1),rmut2,vel2(:,:,1),yys2(:,:,1),prs2(:,1),
     &      SUFBND,SDOT,MOLFRC(:,:,1),BLKTHK(:,:,1),
     &      mat_cal,lcycold,wifsld,oppang,wiface,
     &      LCV_CV,LVEDGE,LBC_SSF,LCYCSF,
     &      MAT_CV,MAT_NO,MAT_CVEXT,MAT_INDEX,CVCENT,
     &      rmu,utau,SFAREA,SFCENT,FRSTCV,CVVOLM,
     &      cord,lacell,lvcell,
     &      POTNAL(:,:,1),POFLUX,POTNBC,FIELD_U,DISALL,WDRBND,
     &      ierror)
!----------------------
      if(ical_prt>=1) then  !1111
        call P_result_S(
     &   iter,ITER_P,time,deltt,DIMV,N_injctr,ffexit,outputfile,
     &   MAT_CAL,MAT_NO,MAT_INDEX,MAT_CVEXT,
     &   NINALL,NPALL,HIS_P,JOUT,
     &   XYZP,UVWP,DIA,PARCEL,PMASS,TMPP,PYS)
      endif
!----------------------
      if( ierror.ne.0 ) goto 9999
      if((ical_mvmsh==2.or.ical_mvmsh==3.or.ical_mvmsh==4)
     &   .and.mod(iter,monitor_stp)==0) then
        call acs_monitor
     &   (msmpl,iter,itery,iterk,iterv,iterp,MNVX,
     &    repsy,aepsy,repsk,aepsk,repsv,aepsv,repsp,aepsp,
     &    erry,errk,errv,errp,errsdy,errsdk,errsdv,
     &    reps_FAI,aeps_FAI,err_FAI,
     &    reps_A,aeps_A,err_A,iter_FAI,iter_A,errsdA,errsdFAI,
     &    iter_POTN,reps_POTN,aeps_POTN,err_POTN,POTNAL(:,:,1),
     &    time,deltt,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &    MHD_A,MHD_FAI,
     &    CVVOLM,tmp(:,1),vel(:,:,1),yys(:,:,1),prs(:,1),aks(:,:,1),
     &    ccc ,rho(:,1))
        goto 1500
      endif
!----------------------
!--< 4.1 statistics >--
!----------------------
        call cal_statis(iter,vel(:,:,1),vlasum,
     &           CVCENT,MAT_NO,MAT_CVEXT,MAT_INDEX,
     &           pp0(:,1),prs(:,1),rho(:,1),tmp(:,1),
     &           yys(:,:,1),aks(:,:,1),rmut)
!------------------------
!--< 4.2 restart file >--
!------------------------
      call write_restart(
     &  iter,time,deltt,ffexit,outputfile,dvdd,
     &  pp0(:,1),prs(:,1),rho(:,1),hhh(:,1),yys(:,:,1),vel(:,:,1),
     &  aks(:,:,1),rva(:,1),
     &  MHD_A(:,:,:,1),MHD_FAI(:,:,1),MHD_CRNT,
     &  MOLFRC(:,:,1),SITDEN(:,:,1),BLKTHK(:,:,1),
     &  rho2(:,1),tmp2(:,1),yys2(:,:,1),vel2(:,:,1),prs2(:,1),hhh2(:,1),
     &  rva2(:,1),dvdd2,
     &  cord,lacell,lvcell,
     &  POTNAL(:,:,1),WDRBND,
     &  ierror)
      if( ierror.ne.0 ) goto 9999
!---------------------------------------
!--< 4.2-2 restart file for particle >--
!---------------------------------------
      if(ical_prt>=1) then  !.and.iter>=ical_prt) then
        call write_particle_restart(
     &  iter,time,deltt,ffexit,outputfile,MXPRT,MEP,INSIDE,
     &  PMASS,DIA,XYZP,UVWP,PARCEL,RHOP,TMPP,
     &  HEIGHT,FU,FV,FW,PYS,HIS_P,JOUT,
     &  ierror)
      endif
      if( ierror.ne.0 ) goto 9999
!--------------------------
! --- < 4.3 result file >--
!--------------------------
!------------------------------------
! --- < 4.4 result animation file >--
!------------------------------------
      call write_anim(iter,time,deltt,ffexit,LCV_CV,LVEDGE,LBC_SSF,
     &                pp0(:,1),prs(:,1),rho(:,1),tmp(:,1),yys(:,:,1),
     &                vel(:,:,1),aks(:,:,1),
     &                MAT_NO,MAT_CVEXT,MAT_INDEX,CVCENT,
     &                SUFBND,SDOT,MOLFRC(:,:,1),WDOT,
     &                ierror)
!----------------------------
!--< 4.4 sound force file >--
!----------------------------
      call write_force_history
     & (iter,time,deltt,ffexit,
     &  MAT_NO,LVEDGE,LBC_SSF,SFAREA,SFCENT,
     &  vel,prs,340.d0,ierror)
!
      call write_cdcl
     & (iter,time,deltt,ffexit,
     &  LVEDGE,LBC_SSF,SFAREA,SFCENT,
     &  prs(:,:),rho(:,:),ierror)
! 
      if(NPE.eq.1) then
        call write_probe(iter,time,deltt,ffexit,
     &    MAT_CVEXT,MAT_INDEX,CVCENT,
     &    prs,pp0,rho,tmp,yys,vel,aks,
     &    ierror)
      else
        call write_probe_para(iter,time,deltt,ffexit,
     &    MAT_CVEXT,MAT_INDEX,CVCENT,
     &    prs,pp0,rho,tmp,yys,vel,aks,
     &    ierror)
      endif
!
      call write_fluid_force (iter,time,deltt,ffexit,
     &  MAT_NO,LVEDGE,LBC_SSF,SFAREA,SFCENT,
     &  FRSTCV,vel,prs,rmut,rho,utau,ierror)
!
!
!--------------------------
! --- Write user output >--
!--------------------------
 1500 continue
       if(outusr.eq.usryes) then
         if(ical_prp) then
           if(sw) then
             call cal_t2cp(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,prs,rho,rmu,rmd,cps,cp,cr)
           else
             call cal_t2hcp(1,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,hhs,cps,rvx(:,1),cp,cr,prs)
           endif
         endif
         call user_output(
     &        MXMAT,MXCOMP,MXRANS,MXCV,MXALLCV,
     &        MXALLCVR,NMAT,NCOMP,NRANS,
     &        MXALLCV2,MXALLCV_RAD,
     &        NCV,nflud,
     &        nofld,MAT_NO,MAT_CVEXT,MAT_DCIDX,LCV_CV,
     &        NPE,my_rank,root,
     &        CVCENT,CVVOLM,disall,
     &        iter,time,deltt,spcnam,wm,
     &        aks(:,:,1),
     &        prs(:,1),rho(:,1),tmp(:,1),yys(:,:,1),
     &        vel(:,:,1),rmut,RVA(:,1),
     &        prs2(:,1),rho2(:,1),tmp2(:,1),yys2(:,:,1),
     &        vel2(:,:,1),rmut2,RVA2(:,1),
     &        nbcnd,nobcnd,MXSSFBC,MXCVFAC,MXCVFAC2,
     &        LBC_INDEX,LBC_SSF,LVEDGE,LCYCSF,SFAREA,SFCENT,locmsh,
     &        HTFLUX,htcbnd,mtcbnd,
     &        rmu,utau,cp,cps,
     &        radflag,radheatflux
     &        )     
      endif
!
      call unlink('output')
      if( .not.ffexit ) goto 1000
!
      if( ierror.ne.0 ) goto 9999
!
      if(my_rank.eq.ROOT) then 
        write(ifll,3020)
        write(ifll,3010)
        write(ifll,3020)
      endif
!
      return
!
 9999 continue
      if(my_rank.eq.ROOT) then 
        write(ifle,*) 'main_FFLOW'
        write(ifll,3025)
        write(ifle,3015)
        write(ifll,3025)
        if( ifll.ne.ifle ) then
          write(ifll,'(1x,a)') '### execution terminated abnormaly'
        endif
 3010   FORMAT(2X,'|',40X,'execution terminated normaly',40X,'|')
 3015   FORMAT(2X,'|',40X,'execution terminated abnormaly',38X,'|')
 3020   FORMAT(2X,'|',108('*'),'|')
 3025   FORMAT(2X,'|',108('#'),'|')
 3030   FORMAT(2X,'|',40X,a,36X,'|')
      endif
!---------------------
! --- abnormaly stop 
!---------------------
      CALL FFRABORT(1,'main_FFLOW')
!
      end subroutine main_FFLOW
!
!velbnd
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine rans_admin(
     &  deltt,iter,time,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,CVVOL0,FRSTCV,disall,xta,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  rva,rmu,rmut,rho,tmp,vel,prs,yys,mdot,EVAP_Y,  
     &  rva2,rmu2,rmut2,rho2,tmp2,vel2,yys2,
     &  aksbnd,aks,MASS_I,
     &  grdc,diag,kdbk,dtau,
     &  rvd,dsclk,cp,cofd,RMD,
!     &  LCYCOLD,wifsld,OPPANG,vctr,FIELD_U,POTNAL,
     &  LCYCOLD,wifsld,OPPANG,FIELD_U,POTNAL,
     &  pp0,vflc,yplusf,rva_s,utau,LFUTAU,VELBND,
     &  iterk,repsk,aepsk,errk,ierror) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!       dsclk  <=  RMX
!       cofd   <=  cr
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only      : ifle
      use module_flags,only   : intgke,eulere
      use module_material,only: icnvke,lmtrke,ivscke,cal,nofld,
     &                          porosty,sthrmu2,relaxrs
      use module_rans,only    : sgmaks,akslw,sgmaksm,akshg
      use module_rans,only    : sgmDES,cb2,Ck,Ce
      USE module_usersub,ONLY : src_rans,usrno,usryes,src_BC_E2P,
     &                          src_BC_CAVI
! --- 
      use module_scalar,only  : icalke,icalrsm,icalaph,
     &                          ike,irsm,iaph,ivof,ides,sclname,ike_c,
     &                          rns_scl,calgeq,calxi,caltemp,calh,
     &                          igeq,ixi,itemp,ihflmlt,
     &                          ical_cavi,icavi,ivold,
     &                          Axis_DIR,FSTR_F,CAVI_B,PWR,
     &                          ical_FC,ical_s,
     &                          ixi_C,ixi_WC,ixi_2,ixi_X,idifflm,
     &                          ORG,Blog,BlogM,consF,FltRans,Axis_DIR,
     &                          x3min,x3max,comb_mdl,
     &                          x1min,x1max,x2min,x2max,
     &                          ikles
      use module_model,only   : icaltb,ke,RSM,ke_low,RNG,CHEN,SDES,KE2S,
     &                          sles,ical_vect,expfact,ical_dens,u_func,
     &                          KLES
      use module_Euler2ph,only: ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      use module_vof     ,only: ical_vof,intervof,LG,LS,change_phase,
     &                          cnvvof,cicsam
      use module_metrix,only  : dkdt  =>d2work1,
     &                          diagk =>d2work2,
     &                          dfk   =>d2work4,
     &                          dflxk =>d2work3,
     &                          akstmp=>d1work1,
     &                          flmspd=>W1K11
      use module_metrix,only  : errrans
      use module_metrix,only  : rcomp,sinj_comp,eva_comp,vctr
      use module_boundary,only : kdbcnd,boundName,nobcnd,
     &                           LBC_INDEX,nbcnd,MAT_BCIDX
      use module_model,only    : idrdp,mach0,comp,ical_prt,
     &                           PEFC,PAFC,MCFC,SOFC,AFC,MDFC
      use module_initial ,only : rho0,rho02,Ke_FC,Kc_FC,Pv_FC,
     &                           Kap_FC,sigma_FC,contang_FC,
     &                           surf_tenFC,p0
      use module_gravity ,only : ggg
      use module_FUEL    ,only : No_Mem,No_AGDL,No_CGDL,vap_no
      use module_species ,only : r_wm
      use module_gravity ,only : ggg,buoyancy
      use module_time,    only : i_steady
      use module_particle,only : iflag_eva
      
!
! 1. Reynolds Averaged Navier-Stokes model
!
      implicit none
!
! --- [dummy arguments]
! 
      real*8 ,intent(in)    :: deltt,time
      integer,intent(in)    :: iter
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO(       0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV(       MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF   (   MXSSFBC)
      real*8 ,intent(in)    :: SFAREA    (4, MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3, MXCVFAC)
      real*8 ,intent(in)    :: wiface    (   MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (   MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (   MXALLCV)
      real*8 ,intent(in)    :: FRSTCV    (   MXssfbc)
      real*8 ,intent(in)    :: DISALL    (   MXCV)
      real*8 ,intent(inout)    :: rva       (   MXCVFAC,2)
      real*8 ,intent(in)    :: rmu       (   MXALLCV)
      real*8 ,intent(in)    :: rmut      (   MXALLCV)
      real*8 ,intent(inout) :: rho       (   MXALLCV  ,2)
      real*8 ,intent(in)    :: tmp       (   MXALLCV,2)
      real*8 ,intent(inout) :: vel       (   MXALLCV,3)

      real*8 ,intent(in)    :: rva2       (   MXCVFAC2,2)
      real*8 ,intent(in)    :: rmu2       (   MXALLCV2)
      real*8 ,intent(in)    :: rmut2      (   MXALLCV2)
      real*8 ,intent(in)    :: rho2       (   MXALLCVC,2)
      real*8 ,intent(in)    :: tmp2       (   MXALLCV2)
      real*8 ,intent(in)    :: vel2       (   MXALLCV2,3)

      real*8 ,intent(in)    :: aksbnd(       MXssfbcR,MXrans)
      real*8 ,intent(inout) :: aks   (       MXALLCVR,MXRANS,2)
      real*8 ,intent(inout) :: rvd       (   MXCVFAC,3)
!      real*8 ,intent(inout) :: tau       (   MXALLCV)
      real*8 ,intent(inout) :: dsclk     (   MXALLCV)
      real*8 ,intent(inout) :: cofd      (   MXALLCV)
      real*8 ,intent(inout) :: grdc      (   MXALLCV,3,3)
      real*8 ,intent(inout) :: dtau      (   MXALLCV,3)
      real*8 ,intent(inout) :: diag      (   MXALLCV)
      INTEGER,INTENT(INOUT) :: kdbk      (   MXCVFAC)
      real*8 ,intent(inout) :: mass_i    (   MXALLCV2)
      integer,intent(out)   :: ierror,iterk(MXrans)
      real*8 ,intent(out)   :: repsk(mxrans),aepsk(MXrans),errk(MXrans)
      integer,intent(in)    :: LCYCOLD   (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG    (MXSSFBC_SLD)
!      integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
      real*8 ,intent(in)    :: pp0(MXALLCV)
      real*8 ,intent(in)    :: vflc(MXCV_F)
      real*8 ,intent(in)    :: yplusf(MXCV_F)
      real*8 ,intent(in)    :: rva_s(MXCVFAC_F)
      real*8 ,intent(in)    :: prs(MXALLCV)
      REAL*8 ,INTENT(INOUT) :: mdot(MXCV,3)
      real*8 ,intent(in)    :: YYS ( MXALLCV,MXCOMP)
      real*8 ,intent(in)    :: YYS2( MXALLCV2,MXCOMP)
      REAL*8 ,INTENT(IN)    :: POTNAL(MXALLCVP,MXPOTN)
      REAL*8 ,INTENT(IN)    :: utau   (0:MXSSFBC)
      REAL*8 ,INTENT(INOUT) :: RMD   (        MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cp    (        MXALLCV)
      REAL*8 ,INTENT(INOUT) :: XTA    (  MXCVFAC)
      integer,intent(in)    :: LFUTAU (     MXCV)
      REAL*8 ,INTENT(INOUT) :: VELBND(       MXSSFBC,3)
      REAL*8 ,INTENT(INOUT) :: EVAP_Y(MXALLCV_P,MXCOMP)
!
! --- [local entities]
!
      real*8  :: duma(1),dll
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp,ndf,nq,ierr1=0
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IIMAT,ICVS,ICVE,ICVL
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: IDCS,IDCE,ICFL,ICFS,ICFE,no
      integer :: IMODE,ICVLA,ICVLB
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      real*8  :: vdlt,rdeltt,big,dflx,dum1,dum2,dum3,dum4,dum5,dum6
      real*8  :: wi1,wi2,grx,gry,grz,dx,dy,dz,gf0,gf1,gf2,gf3
      real*8  :: grdf(3),dxx,dyx,dzx,evapmdl,dl,dlvect,dumd(3)
      real*8  :: repski,aepski,void1,void2
      real*8  :: u,v,w,u2,v2,w2,p,dens,dens2,T_WALL,T,T2,mass
      integer :: iterki,IMAT_U,IBFS,IBFE,IBFL
      logical :: viscous,noviscs,othrs
      integer :: cnv_rans
      character*80 :: BC_name
      character :: sn
      real*8  :: Ramda_l,Ramda_g,dpds,Dc,dumsave=0.d0,dumm=0.d0
      real*8  :: Cz=1.d0,Cx=2.d0
      integer :: ndiag
!
!      real*8,allocatable  :: errabs(:)
!
      ALLOCATE (dkdt (MXALLCVR,MXrans),stat=ierr1)
      ALLOCATE (diagk(MXALLCVR,MXrans),stat=ierr1)
      ALLOCATE (dfk  (MXALLCVR,MXrans),stat=ierr1)
      ALLOCATE (dflxk(MXCVFACR,MXrans),stat=ierr1)
      ALLOCATE (akstmp(MXALLCVR),stat=ierr1)
      if(ierr1/=0) call FFRABORT(1,'allocating error in rans_admin')
!--------------------------------------------------------------------
! --- 
!--------------
      vdlt=1.d0/deltt
      rdeltt=1.d0/deltt
      if(i_steady==3) then
        rdeltt=1.d0    !1212
      endif
      kdbk=0.d0
      diagk=0.d0
      dkdt=0.d0
      dfk =0.d0
!
!-< 1. Set boundary condition >-
! --- kdbk(ICF) k-eps BC 
!              =0:periodic/CV-face
!              =1:other
!              =2:sysmetric/kkneum
!-------------------
! --- BC condition
!-------------------
      call bc_rans(LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     &             rva(:,1),rva2(:,1),aksbnd,aks,kdbk)
!
      if(ieul2ph>0) then
      endif
!
!---------------------------
!-< 2. Calculate tendency >-
!---------------------------
!------------------------
!--< 2.1 source term >-- 
!------------------------
!      tau=0.d0
!
!-----------------------------------------
! --- Scalar User source term (subroutine)
!-----------------------------------------
!
      if(src_rans.eq.usryes) then
      endif
!
!-----------------------------------------
! --- VOID User source term (subroutine) 
!-----------------------------------------
!
      if(src_BC_E2P.eq.usryes.and.ieul2ph>0) then
      endif
!
!----------------------------
! --- K-E model source term 
!----------------------------
!
      if(icaltb==ke.or.
     &   icaltb==ke_low.or.
     &   icaltb==RNG.or.
     &   icaltb==CHEN.or.
     &   icaltb==KE2S.or.
     &   icaltb==KLES.or.
     &   icaltb==SDES) 
     &  then
        if(ieul2ph>0) call FFRABORT(1,'ERR:E2P NOT Support KE model')
        call keps_src(
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  DISALL,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  grdc,dtau,rmut,rmu,rho,vel,tmp(:,1),aks,dkdt,
     &  LCYCOLD,wifsld,utau,LFUTAU,
     &  diagk,vctr)
      endif
!--------------
! --- C_dush^2
!--------------
      if(u_func(1)==1) then
        call C_dash_src(
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  DISALL,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  grdc,dtau,rmut,rmu,rho,vel,tmp(:,1),aks,yys,dkdt,
     &  LCYCOLD,wifsld,utau,LFUTAU,
     &  diagk,vctr)
      endif
!
      if(icaltb==KLES) then
      endif
!

!---------------------------
! --- RSM model  source term
!---------------------------
!
      if(icaltb.eq.RSM) then

      endif
!-----------------------------------
! --- Cavitation model  source term 
!-----------------------------------
      if(ical_cavi==1) then 
      endif
!
! --- 
!
      if(ical_FC==PEFC) then
        call FC_src(0,MXrans,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  rho,rho2(:,1),tmp(:,1),aks,dkdt,diagk,prs,yys,mDOT,
     &  LCYCOLD,wifsld)
!        if(buoyancy==1) then
!        do IIMAT=1,NMAT
!        IMAT=MAT_NO(IIMAT)
!        if(.not.mat_cal(IIMAT).or.IMAT.le.0) cycle
!        IMAT_U=nofld(IMAT)
!      if((IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem)) 
!     &  cycle
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_CVEXT(IIMAT)
!        DO ICVL=ICVS,ICVE
!        dum1=aks(ICVL,ical_s,1)
!        do IV=1,3
!        dvdt(ICVL,IV)=dvdt(ICVL,IV)
!     &               +dum1*ggg(IV)*(rho(ICVL,1)-rho2(ICVL,1))
!        enddo
!        enddo
!        enddo

!        endif
      endif
!--------------------------------------
! --- Euler tow phase flow source term
! --- Vaporation: mass_i(:)>0
! --- Condensation : mass_i(:)<0
!--------------------------------------
      if(ieul2ph>0) then
      endif
!--------------------------------------
! --- VOF equation source term
!--------------------------------------
      if(ical_vof==1)then
      endif
      if(ical_vof==1.and.change_phase.eq.1) then
      endif
!----------------------------------
! --- K-E Model wall function 
!----------------------------------
      if(icaltb==KE2S.or.icaltb==ke.or.
     &   icaltb==ke_low.or.icaltb==RNG.or.icaltb==CHEN) 
     &  then
        big=1.d20/deltt
        if(ical_vect) then  !NNOTVECT
          call bc_wallke_VECT(big,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &    SFAREA,SFCENT,CVCENT,CVVOLM,rmu,utau,rho,vel,aks,dkdt,diagk)
        else
          call bc_wallke(big,LVEDGE,LBC_SSF,LCYCSF,mat_cal,MAT_NO,
     &    SFAREA,SFCENT,CVCENT,CVVOLM,FRSTCV,
     &    rmu,rho,vel,aks,dkdt,diagk)
        endif
      endif
!--------------------------------------- !NOT finished 
! --- diffusion combustion source term 
!---------------------------------------
      if(calxi.and.idifflm/=ORG) then !NOT finished =>ixi  !ORG  idifflm/=ORG
        if(ical_prt>=1.and.iflag_eva==1) then !Evaperation source
          do IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle 
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICOM=1,ncomp
          do ICVL=ICVS,ICVE
          dkdt(ICVL,ixi)=dkdt(ICVL,ixi)+EVAP_Y(ICVL,ICOM)
          enddo
          enddo
          enddo
        endif
      endif
!--------------------------------------------
! --- 
!--------------------------------------------
      if(calxi) then            !NOT finished =>ixi_2
        if(idifflm==FltRans) then
          akstmp(:)=aks(:,ixi,1)
          call grad_cell(1,14,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,akstmp,grdc(:,:,1)) 
!
          dum1=1.d0/sgmaks(ixi_2)
          dum2=1.d0/sgmaksm(ixi_2)
          do IIMAT=1,NMAT               !eddy diffusion 
          if(.not.mat_cal(IIMAT)) cycle 
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE 
          dum3=rho(ICVL,1)*tmp(ICVL,2)*
!     &         (rmut(ICVL)*dum1)*
     &         
     &         (grdc(ICVL,1,1)**2
     &         +grdc(ICVL,2,1)**2
     &         +grdc(ICVL,3,1)**2)
          dum4=aks(ICVL,ixi_X,1)*rho(ICVL,1)
          dum4=Cx*aks(ICVL,ike(2),1)/
     &   (aks(ICVL,ike(1),1)+SML)*rho(ICVL,1)*aks(ICVL,ixi_2,1)
          dkdt(ICVL,ixi_2)=dkdt(ICVL,ixi_2)+(dum3-dum4)
!     &                    +max(0.d0,dum3-dum4) 
          dum3=Cx*aks(ICVL,ike(2),1)/
     &   (aks(ICVL,ike(1),1)+SML)*rho(ICVL,1)
          diagk(ICVL,ixi_2)=diagk(ICVL,ixi_2)+dum3
          enddo
          enddo
        endif
      endif
!-------------------------------------------
! --- 
!-------------------------------------------
      if(calxi) then !NOT finished => ixi_C
        if(idifflm==BlogM) then
          call FFRABORT(1,'ERR: NOT finished for Blog=3')
          do IIMAT=1,NMAT       ! !C sourve term
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1 
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          dkdt(ICVL,ixi_C)=
     &         dkdt(ICVL,ixi_C)+rho(ICVL,1)*aks(ICVL,ixi_WC,1)
          diagk(ICVL,ixi_C)=0.d0
          enddo
          enddo
        endif
!
        if(idifflm==Blog) then
          do IIMAT=1,NMAT       !C sourve term
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          dkdt(ICVL,ixi_C)=
!     &         dkdt(ICVL,ixi_C)+rho(ICVL,1)*aks(ICVL,ixi_WC,1)
     &         dkdt(ICVL,ixi_C)+aks(ICVL,ixi_WC,1)
          diagk(ICVL,ixi_C)=0.d0
          enddo
          enddo
        endif
      endif
!---------------------------------------
! --- Other scalar's wall function ???
!---------------------------------------
!
!      othrs=.false.
!      do IMD=1,nrans
!      if(othrs) then
!      endif
!      enddo
!
!----------------------------------------------
!--< 2.2 time difference term & source term >--
!----------------------------------------------
!
      do 250 IMD=1,nrans
      do 200 IIMAT=1,NMAT    !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) goto 200
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
!
      if(ical_vof==1.and.IMD.eq.ivof) then
      elseif
     &  ((icaltb==KE2S.or.icaltb==ke.or.icaltb==
     &    ke_low.or.icaltb==RNG.or.icaltb==CHEN)
     &  .and.(IMD.eq.ike(1).or.IMD.eq.ike(2))) then
        if(ical_dens==4) then
        else
          do ICVL=ICVS,ICVE
          dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          diagk(ICVL,IMD)=dum1
     &      +CVVOLM(ICVL)*(rho(ICVL,1)*diagk(ICVL,IMD))
          dkdt(ICVL,IMD)=dum1*(aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
     &                +CVVOLM(ICVL)*(dkdt(ICVL,IMD))
          enddo
        endif
      elseif(u_func(1)==1.and.ike_c==IMD)then
        do ICVL=ICVS,ICVE
          dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          diagk(ICVL,IMD)=dum1
     &      +CVVOLM(ICVL)*(rho(ICVL,1)*diagk(ICVL,IMD))
          dkdt(ICVL,IMD)=dum1*(aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
     &                +CVVOLM(ICVL)*(dkdt(ICVL,IMD))
        enddo
      elseif(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) then
      elseif(ical_cavi==1.and.(IMD==icavi.or.IMD==ivold)) then
      elseif(ical_FC==PEFC.and.IMD==ical_s) then
!        IMAT_U=nofld(IMAT)
!        if(.NOT.(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem)) 
!     &  cycle
        do ICVL=ICVS,ICVE
        dum1=rho02(IIMAT)*CVVOL0(ICVL)*rdeltt
        dum2=porosty(IMAT)
        diagk(ICVL,IMD)=dum1+CVVOLM(ICVL)*diagk(ICVL,IMD)
        dkdt(ICVL,IMD)=dum2*dum1*(aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
     &                +dum2*CVVOLM(ICVL)*(dkdt(ICVL,IMD))
        enddo
      elseif(calxi.and.
     &  (IMD==ixi.or.IMD==ixi_C.or.IMD==ixi_2.or.IMD==ixi_X)) then
        if(IMD==ixi_X) cycle
        if(IMD==ixi.or.
     &    (idifflm==FltRans.and.IMD==ixi_2).or.
     &    (idifflm==BlogM.and.IMD==ixi_C).or.
     &    (idifflm==Blog .and.IMD==ixi_C)
     &     ) then
          do ICVL=ICVS,ICVE
          dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
          diagk(ICVL,IMD)=dum1+CVVOLM(ICVL)*diagk(ICVL,IMD)
          dkdt(ICVL,IMD)=dum1*(aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
     &                +CVVOLM(ICVL)*(dkdt(ICVL,IMD))
          enddo
        endif
      else
        do ICVL=ICVS,ICVE
        dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
        diagk(ICVL,IMD)=dum1+CVVOLM(ICVL)*diagk(ICVL,IMD)
        dkdt(ICVL,IMD)=dum1*(aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
     &                +CVVOLM(ICVL)*(dkdt(ICVL,IMD))
        enddo
      endif
!
      if(i_steady==3.and..false.) then
        do ICVL=ICVS,ICVE
        dum1=rho(ICVL,2)*CVVOL0(ICVL)*rdeltt
        dkdt(ICVL,IMD)=dkdt(ICVL,IMD)-dum1*aks(ICVL,IMD,2)
        enddo
      endif
 200  enddo
!
!
!--------------------------
!--< 2.3 diffusion term >--
!--------------------------
!
!      if(ivscke.ne.cal) goto 219
!
      if(ical_vof==1.and.IMD.eq.ivof) then
      endif
      if(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) then
        dflxk(:,IMD)=0.d0
        cycle
      endif
      if(ical_cavi==1.and.(IMD==icavi.or.IMD==ivold)) then
        dflxk(:,IMD)=0.d0
        cycle
      endif
      if(calxi.and.
     &  (IMD==ixi.or.IMD==ixi_C.or.IMD==ixi_2.or.IMD==ixi_X)) then
        if(IMD==ixi_X) then
          dflxk(:,IMD)=0.d0
          cycle
        endif
        if(IMD==ixi_2.and.idifflm/=FltRans) then
          dflxk(:,IMD)=0.d0
          cycle
        endif
      endif
!
      akstmp(:)=aks(:,IMD,1)
      call grad_cell(1,14,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,akstmp,grdc(:,:,1))
!
      call dc_symprv
     &(1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & LCYCOLD,wifsld,OPPANG,
     & SFAREA,SFCENT,grdc(:,:,1),1,1)
!
      if(icaltb==SDES.and.IMD==ides) then
        dum1=1.d0/sgmDES
        do 221 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dfk(ICVL,IMD)=(rmu(ICVL)+rmut(ICVL))*dum1 
        enddo
  221   enddo
      elseif(icaltb==KLES.and.IMD==ikles) then
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dfk(ICVL,IMD)=(rmut(ICVL))
        enddo
        enddo
      elseif(ical_FC==PEFC.and.IMD==ical_s) then
! --- surface tension (diffusion term)
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT<=0) cycle
        IMAT_U=nofld(IMAT)
        if(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem) then
          if(surf_tenFC(IIMAT)==0) cycle 
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          dum3=porosty(IMAT)
          do ICVL=ICVS,ICVE
          dum1=aks(ICVL,IMD,1)
          Ramda_l=dum1*(2.d0-dum1)
          Ramda_g=1.d0-Ramda_l
          dum2=sthrmu2(IMAT)
          dum4=sigma_FC(IIMAT)*
     &       cos(contang_FC(IIMAT)*3.14159/180.d0)*
     &       dsqrt(dum3/Kap_FC(IIMAT))
          dum5=2.d0*2.12d0*(1.d0-dum1)
     &      -1.417d0-3.d0*1.262*(1.d0-dum1)**2
          dpds=dum4*dum5
          Dc=-Kap_FC(IIMAT)*dum1**3*Ramda_g/dum2*dpds
          if(Dc<1.d-10) then
            Dc=Kap_FC(IIMAT)*(dum1+0.01d0)*rho02(IIMAT)*9.8d0/dum2*
     &         (3.7d0*1.73d-4*(exp(-3.7d0*(dum1-0.494))
     &     +exp(3.7*(dum1-0.494))))
            dfk(ICVL,IMD)=Dc
          else
            dfk(ICVL,IMD)=Dc*rho02(IIMAT)
          endif
          enddo
          if(IMAT_U==No_Mem)then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            dum3=porosty(IMAT)
            do ICVL=ICVS,ICVE
            p=prs(ICVL)
            T=tmp(ICVL,1)
            rcomp(:)=YYS(ICVL,:)
            dum2=Pv_FC(IIMAT)
            dum3=aks(ICVL,ical_s,1)
            call ramd(rcomp,p,T,dum2,dum3,dum4,r_wm,ncomp,vap_no,dum1)
            if(dum1<=2.d0) then
              dum3=1.d0
            elseif(dum1>2.d0.and.dum1<=3.d0) then
              dum3=1.d0+2.d0*(dum1-2.d0)
            elseif(dum1>3.d0.and.dum1<=4.d0) then
              dum3=3.d0-1.38d0*(dum1-3.d0)
            else
              dum3=2.563d0-0.33d0*dum1+0.0264*dum1**2-0.000671*dum1**3
            endif
            dum2=rho02(IIMAT)*1.d-10*exp(2416.d0*(1.d0/303.d0-1.d0/T))
            dfk(ICVL,IMD)=dfk(ICVL,IMD)+dum2*dum3   !Very small=10^-10
            enddo
          endif
        endif
        enddo
      elseif(calgeq.and.IMD==igeq) then
        dum1=1.d0/sgmaks(IMD)
        dum2=1.d0/sgmaksm(IMD)
        do IIMAT=1,NMAT   !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1 
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dfk(ICVL,IMD)=rmu(ICVL)*dum2+rmut(ICVL)*dum1
        enddo
        enddo
      elseif(calxi.and.
     &   (IMD==ixi.or.IMD==ixi_C.or.IMD==ixi_2)
     &   ) then
        if(IMD==ixi.or.
     &    (idifflm==FltRans.and.IMD==ixi_2).or.
     &    (idifflm==BlogM.and.IMD==ixi_C).or.
     &    (idifflm==Blog .and.IMD==ixi_C)
     &     ) then
          if(IMD==ixi.and.idifflm==ORG) then  !ORG  idifflm==ORG 
            dum1=1.d0/sgmaks(IMD)
            dum2=1.d0/sgmaksm(IMD)
            do IIMAT=1,NMAT   !ICV=1,NCV
            if(.not.mat_cal(IIMAT)) cycle
            IMAT=MAT_NO(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1 
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE 
            dfk(ICVL,IMD)=rmu(ICVL)*dum2+rmut(ICVL)*dum1  !NOT Finished
            enddo
            enddo
          elseif((IMD==ixi.or.IMD==ixi_2).and.idifflm==FltRans) then
            dum1=1.d0/sgmaks(IMD)
            dum2=1.d0/sgmaksm(IMD)
            do IIMAT=1,NMAT   !ICV=1,NCV 
            if(.not.mat_cal(IIMAT)) cycle 
            IMAT=MAT_NO(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1 
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
            dfk(ICVL,IMD)=rmut(ICVL)*dum1
!rmut(ICVL)*dum1+rho(ICVL,1)*tmp(ICVL,2)
            enddo
            enddo
          elseif((IMD==ixi).and.
     &    (idifflm==consF.or.idifflm==Blog.or.idifflm==BlogM)) then
            dum1=1.d0/sgmaks(IMD)
            dum2=1.d0/sgmaksm(IMD)
            do IIMAT=1,NMAT   !ICV=1,NCV 
            if(.not.mat_cal(IIMAT)) cycle 
            IMAT=MAT_NO(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1 
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
            dfk(ICVL,IMD)=rmd(ICVL)/cp(ICVL)+rmut(ICVL)*dum1
            enddo 
            enddo 
          elseif(IMD==ixi_C.and.(idifflm==Blog.or.idifflm==BlogM)) then
            dum1=1.d0/sgmaks(IMD)
            dum2=1.d0/sgmaksm(IMD)
            do IIMAT=1,NMAT   !ICV=1,NCV 
            if(.not.mat_cal(IIMAT)) cycle 
            IMAT=MAT_NO(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1 
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
!            dfk(ICVL,IMD)=rho(ICVL,1)*tmp(ICVL,2)+rmut(ICVL)*dum1
            dfk(ICVL,IMD)=rmd(ICVL)/cp(ICVL)+rmut(ICVL)*dum1 !same as Z
            enddo
            enddo
          endif
        endif
      else
        dum1=1.d0/sgmaks(IMD)
        do 211 IIMAT=1,NMAT   !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 211
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dfk(ICVL,IMD)=rmu(ICVL)+rmut(ICVL)*dum1
        enddo
  211   enddo
      endif
!
      if(NMAT>1) then
        akstmp(:)=dfk(:,IMD)
        call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,akstmp)
        dfk(:,IMD)=akstmp(:)
      endif
!
      do 212 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) goto 212
      IMAT=MAT_NO(IIMAT)
      if(IMAT<=0) cycle
      IMAT_U=nofld(abs(IMAT))
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
      ICVA=LVEDGE(1,ICFL)
      ICVB=LVEDGE(2,ICFL)
      wi1=wiface(ICFL)
      wi2=1.d0-wiface(ICFL)
      dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
      dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
      dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
      dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
      dx=dx/dl
      dy=dy/dl
      dz=dz/dl
      dlvect=abs(dx*SFAREA(1,ICFL)
     &          +dy*SFAREA(2,ICFL)
     &          +dz*SFAREA(3,ICFL))
      grx=wi1*grdc(ICVA,1,1)+wi2*grdc(ICVB,1,1)
      gry=wi1*grdc(ICVA,2,1)+wi2*grdc(ICVB,2,1)
      grz=wi1*grdc(ICVA,3,1)+wi2*grdc(ICVB,3,1)
      gf1=(aks(ICVB,IMD,1)-aks(ICVA,IMD,1))/(dlvect*dl)
      gf2=grx*SFAREA(1,ICFL)+gry*SFAREA(2,ICFL)+grz*SFAREA(3,ICFL)
      gf3=(grx*dx+gry*dy+grz*dz)
      dum1=wi1*dfk(ICVA,IMD)+wi2*dfk(ICVB,IMD)
!      dum2=dfk(ICVA,IMD)
!      dum3=dfk(ICVB,IMD)
!      dum1=dum2*dum3/(dum2*wi1+dum3*wi2+SML)
!-------------------(1)
!      dflxk(ICFL,IMD)=dum1*(gf1+expfact*(gf2-gf3))
      dflxk(ICFL,IMD)=dum1*(gf1+0.d0*(gf2-gf3))
!------------------- antei(2)
!      grdf(1)=SFAREA(1,ICFL)*gf1+grx*(SFAREA(1,ICFL)-dx/dlvect)
!      grdf(2)=SFAREA(2,ICFL)*gf1+gry*(SFAREA(2,ICFL)-dy/dlvect)
!      grdf(3)=SFAREA(3,ICFL)*gf1+grz*(SFAREA(3,ICFL)-dz/dlvect)
!      dflxk(ICFL,IMD)=
!     &     (wi1*dfk(ICVA,IMD)
!     &     +wi2*dfk(ICVB,IMD))
!     &     *(SFAREA(1,ICFL)*grdf(1)
!     &      +SFAREA(2,ICFL)*grdf(2)
!     &      +SFAREA(3,ICFL)*grdf(3))
!-------------------
      enddo
  212 enddo  !
!----------------------------------
! --- Spalart and Allmaras model 
!----------------------------------
      IF(icaltb==SDES.and.IMD==ides) then
        dum2=1.d0/sgmDES
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        gf1=GRDC(ICVL,1,1)**2
     &     +GRDC(ICVL,2,1)**2
     &     +GRDC(ICVL,3,1)**2
        dum1=cb2*rho(ICVL,1)*dum2*gf1
        dkdt(ICVL,IMD)=dkdt(ICVL,IMD)+dum1*CVVOLM(ICVL)
        enddo
        enddo
      endif
!
  250 enddo
!------------------------------------
! --- BC for dflxk of all scalar 
!------------------------------------
      call bc_ransdif
     &  (LBC_SSF,LCYCSF,mat_cal,aksbnd,rva(:,1),rva2(:,1),dflxk)
!
      do 260 IMD=1,nrans
      do 220 IIMAT=1,NMAT     !ICF=1,NCVFAC
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      IMAT_U=nofld(abs(IMAT))
      if(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) cycle
      if(ical_cavi==1.and.(IMD==icavi.or.IMD==ivold)) cycle
      if(ical_vof==1.and.(IMD==ivof)) cycle
      if(calxi.and.(IMD==ixi_X.or.
     &  (IMD==ixi_2.and.idifflm/=FltRans))) cycle
      ICFS=MAT_CFIDX(IIMAT-1)+1
      ICFE=MAT_CFIDX(IIMAT)
      do ICFL=ICFS,ICFE
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        dum1=dflxk(ICFL,IMD)*SFAREA(4,ICFL)
        dkdt(ICVA,IMD)=dkdt(ICVA,IMD)+dum1
        dkdt(ICVB,IMD)=dkdt(ICVB,IMD)-dum1
      enddo
  220 continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,dkdt(:,IMD))
      ENDIF
!
  219 continue
!
!---------------------------
!--< 2.4 convection term >--
!---------------------------
      akstmp(:)=aks(:,IMD,1)
!---------------------------
! --- flamelet model
!---------------------------
      if (calgeq.and.IMD.eq.igeq) then
!-------------------------------------
! --- Set Turbulent Burning Velocity
!-------------------------------------
! --- Define a Laminar Burning Velocity from chemical conditions
!-----------------------------------------------------------------
        call cal_bv(MAT_CVEXT,MAT_CAL,MAT_INDEX,aks,pp0,flmspd)
!-----------------------------------------------------------------
! --- Define a Turbulent Burning Velocity from flow conditions
!-----------------------------------------------------------------
        if(icaltb.eq.sles) then
          call cal_trbbv(MAT_CVEXT,MAT_CAL,vflc,flmspd)
!--------------------------------------------
! --- Dumping function for Burning Velocity
!--------------------------------------------
          call cal_dmpbv(MAT_CVEXT,MAT_CAL,yplusf,flmspd)
        endif
!--------------------------------------------
! --- Cal. (Advection + Propagation) Flux 
!--------------------------------------------
	call cal_prpgf
     &          (MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,mat_cal
     &          ,LVEDGE,LBC_SSF,LCYCSF,LCYCOLD,OPPANG
     &          ,SFAREA,SFCENT,wiface,wifsld,CVVOLM
     &          ,rva,akstmp,flmspd,grdc,vctr,rva_s)
!
        call conv_term 
     &  (icnvke,lmtrke,deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,vctr,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &  SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &  rvd(:,1:2),vel,
     &  grdc(:,:,1),rva_s,akstmp,dkdt(:,IMD),
     &  3)
      elseif (calxi.and.
     &  (IMD==ixi.or.IMD==ixi_C.or.IMD==ixi_2.or.IMD==ixi_X))
     &  then
        if(IMD==ixi_X.or.(idifflm/=FltRans.and.IMD==ixi_2)) cycle
!
        if(IMD==ixi.or.
     &    (idifflm==FltRans.and.IMD==ixi_2).or.
     &    (idifflm==BlogM.and.IMD==ixi_C).or.
     &    (idifflm==Blog .and.IMD==ixi_C)
     &     ) then
          call conv_term
     &    (icnvke,lmtrke,deltt,
     &    LVEDGE,LBC_SSF,LCYCSF,vctr,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &    rvd(:,1:2),vel,
     &    grdc(:,:,1),rva,akstmp,dkdt(:,IMD),
     &    5)
        endif
      else
        if(ical_vof==1.and.IMD==ivof) then
        elseif(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) then
        elseif(ical_FC==PEFC.and.IMD==ical_s) then
          call conv_term_FC
     &    (icnvke,lmtrke,deltt,
     &    LVEDGE,LBC_SSF,LCYCSF,vctr,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &    rvd(:,1:2),vel,dtau,POTNAL,
     &    grdc(:,:,1),rva,akstmp,dkdt(:,IMD),rho(:,1),rho2(:,1),
     &    tmp(:,1),prs,yys,aks,FIELD_U,
     &    4)
          do IIMAT=1,NMAT
          if(.not.mat_cal(IIMAT)) cycle
          IMAT=MAT_NO(IIMAT)
          if(IMAT<0) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            dkdt(ICVS:ICVE,IMD)=0.d0
          else
!            IMAT_U=nofld(IMAT)
!            if(IMAT_U==No_AGDL.or.IMAT_U==No_CGDL.or.IMAT_U==No_Mem) cycle
!            ICVS=MAT_CVEXT(IIMAT-1)+1
!            ICVE=MAT_CVEXT(IIMAT)
!            dkdt(ICVS:ICVE,IMD)=0.d0
          endif
          enddo
        elseif(ical_cavi==1.and.(IMD==icavi.or.IMD==ivold)) then
        else
          call conv_term
     &    (icnvke,lmtrke,deltt,
     &    LVEDGE,LBC_SSF,LCYCSF,vctr,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,MAT_DCIDX,
     &    SFAREA,SFCENT,wiface,CVCENT,CVVOLM,
     &    rvd(:,1:2),vel,
     &    grdc(:,:,1),rva,akstmp,dkdt(:,IMD),
     &    5)
        endif
      endif
!
!---------------------------
!-< 3. Update k & epsilon >-
!---------------------------
!--< 3.1 Euler explicit   >-
!---------------------------
      if(intgke.eq.eulere
     &  ) then
!
        do 300 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) goto 300
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dum1=1.d0/(diagk(ICVL,IMD)+SML)
        dkdt(ICVL,IMD)=dum1*dkdt(ICVL,IMD)
        enddo
 300    enddo
!--------------------------
!--< 3.2 Euler implicit >--
!--------------------------
      else
!
        ndf=1 !nrans
	nq=1  !nrans
        rvd=0.d0
        dsclk=0.d0
        duma=0.d0
        cofd=1.d0
!?????????????????????????????????????????
! --- scalar d(s) NOT use diffusion trem 
!        dfk(:,IMD)=0.d0 
!?????????????????????????????????????????
        sn=sclname(IMD)
        if(ical_vof==1.and.IMD.eq.ivof) then
          IMODE=3
          dfk(:,IMD)=0.d0
        elseif(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) then
        elseif(ical_cavi==1.and.IMD==icavi) then
          dfk(:,IMD)=0.d0
          IMODE=4   !3
        else 
          IMODE=4   !3
        endif
!
        if(ieul2ph==0) then
          if(ical_vect) then !NNOTVECT
            rvd=0.d0 
            ndf=1 
            nq=1
            dsclk=0.d0
            duma=0.d0
            cofd=1.d0
            ndiag=1
            call solve_cnvdif_vect
     &           (.true.,1,ndf,nq,sn,time,ndiag,
     &      LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,
     &      MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,
     &      MAT_NO,mat_cal,MAT_CFIDX,
     &      rho,kdbk,rva,rvd(:,1:2),dfk(:,IMD),diagk(:,IMD),
     &      LCYCOLD,wifsld,OPPANG,
     &      cofd,dsclk,duma,dkdt(:,IMD),deltt,vctr,
     &      iterki,repski,aepski,iter,ierror,IMODE)
            repsk(IMD)=repski
            aepsk(IMD)=aepski
            iterk(IMD)=iterki
          else
            call solve_cnvdif
     &      (.true.,1,ndf,nq,sn,time,relaxrs,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,SFCENT,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    rho,kdbk,rva,rvd(:,1:2),dfk(:,IMD),diagk(:,IMD),
     &    LCYCOLD,wifsld,OPPANG,CVVOLM,
     &    cofd,dsclk,duma,dkdt(:,IMD),deltt,
     &    iterki,repski,aepski,iter,ierror,IMODE)
          endif
        endif
!
        repsk(IMD)=repski
        aepsk(IMD)=aepski
        iterk(IMD)=iterki
        if( ierror.ne.0 ) goto 9999
      endif
 260  enddo
!
      if(.false..and.ical_vect) then   !only for KE,RNG,CHEN model
        if(ieul2ph>0.or.ical_cavi==1.or.ical_FC==PEFC.or.ical_vof==1) 
     &  then
          call 
     &   FFRABORT(1,'ERR: VECT-Ver ONLY support KE, RNG CHEN and SDES')
        endif
        rvd=0.d0
        ndf=nrans
	nq=nrans
        dsclk=0.d0
        duma=0.d0
        cofd=1.d0
        ndiag=nrans
        call solve_cnvdif_vect
     &      (.true.,1,ndf,nq,sn,time,ndiag,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,wiface,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    rho,kdbk,rva,rvd(:,1:2),dfk(:,IMD),diagk(:,IMD),
!     &    rho,kdbk,rva,rvd(:,1:2),dfk,diagk,
     &    LCYCOLD,wifsld,OPPANG,
     &    cofd,dsclk,duma,dkdt,deltt,vctr,
     &    iterk,repsk,aepsk,iter,ierror,IMODE)
      endif
!-------------------------------------------------
!--< 3.3 reset k & epsilon at (n+1) time level >--
!-------------------------------------------------
      dum1=0.d0
      do 330 IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
!
      if(IMAT<0) then
        aks(ICVS:ICVE,:,1)=0.d0
        cycle
      endif
!
      IMAT_U=nofld(IMAT)
      do 310 IMD=1,nrans
!      
      if(ical_cavi==1.and.(IMD==ivold.or.IMD==icavi)
     &  .or.(ical_FC==PEFC.and.IMD==ical_s)
     &  .or.(ical_vof==1.and.IMD.eq.ivof)
     &  )then
!
        if(IMD==ivold) cycle
        if(ical_cavi==1.or.IMD==icavi) then
        endif
!
        if(ical_FC==PEFC.and.IMD==ical_s) then
          do ICVL=ICVS,ICVE
          dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
          aks(ICVL,IMD,1)=min(1.d0,max(dum1,0.d0))
          enddo
        endif
!
        if(ical_vof==1.and.IMD.eq.ivof) then
        endif
!        
      elseif(ieul2ph>0.and.(IMD==iaph(1).or.IMD==iaph(2))) then
      elseif(calgeq.and.IMD==igeq) then 
        do ICVL=ICVS,ICVE
        dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
        aks(ICVL,IMD,1)=min(akshg(IMD),max(dum1,akslw(IMD)))
        enddo
      elseif(calxi.and.
     &      (IMD==ixi.or.IMD==ixi_C.or.IMD==ixi_2.or.IMD==ixi_X)
     &   ) then 
        if(IMD==ixi.or.
     &    (idifflm==FltRans.and.IMD==ixi_2).or.
     &    (idifflm==BlogM.and.IMD==ixi_C).or.
     &    (idifflm==Blog .and.IMD==ixi_C)
     &     ) then
          if(comb_mdl==2) then
            if(IMD==ixi_C) then
              do ICVL=ICVS,ICVE
              dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
              aks(ICVL,IMD,1)=min(x3max,max(dum1,x3min))
              enddo
            elseif(IMD==ixi) then
              do ICVL=ICVS,ICVE
              dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
              aks(ICVL,IMD,1)=min(x1max,max(dum1,x1min))
              enddo
            elseif(IMD==ixi_2) then
              do ICVL=ICVS,ICVE
              dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
              aks(ICVL,IMD,1)=min(x2max,max(dum1,x2min))
              enddo
            endif
          elseif(comb_mdl==1) then
            if(IMD==ixi_X) then
              do ICVL=ICVS,ICVE
              dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
              aks(ICVL,IMD,1)=min(x3max,max(dum1,x3min))
              enddo
            elseif(IMD==ixi) then
              do ICVL=ICVS,ICVE
              dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
              aks(ICVL,IMD,1)=min(x1max,max(dum1,x1min))
              enddo
            elseif(IMD==ixi_2) then
              do ICVL=ICVS,ICVE
              dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
              aks(ICVL,IMD,1)=min(x2max,max(dum1,x2min))
              enddo
            endif
          else
            do ICVL=ICVS,ICVE
            dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
            aks(ICVL,IMD,1)=min(akshg(IMD),max(dum1,akslw(IMD)))
            enddo
          endif
!        elseif(idifflm==BlogM.and.IMD==ixi_2) then
!          do ICVL=ICVS,ICVE
!          dum1=aks(ICVL,IMD,1)+dkdt(ICVL,IMD)*relaxrs(IIMAT)
!          aks(ICVL,IMD,1)=min(akshg(IMD),max(dum1,akslw(IMD)))
!          enddo
        elseif(IMD==ixi_2) then
          cycle
        elseif(IMD==ixi_X) then
          cycle
        endif
      else
        do ICVL=ICVS,ICVE
        aks(ICVL,IMD,1)=max(akslw(IMD),aks(ICVL,IMD,1)
     &                 +relaxrs(IIMAT)*dkdt(ICVL,IMD))
        enddo
      endif
  310 enddo
 330  enddo
!--------------------------
! --- Convergence error:
!--------------------------
!
      errk(:)=ZERO
      errrans(:)=ZERO
      Do 420 IIMAT=1,NMAT    !ICV=1,NCVIN
      if(.not.mat_cal(IIMAT)) goto 420
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do 421 IMD=1,nrans
      do ICVL=ICVS,ICVE
      errk(IMD)=errk(IMD)+dkdt(ICVL,IMD)*dkdt(ICVL,IMD)
      errrans(IMD)=errrans(IMD)+aks(ICVL,IMD,1)*aks(ICVL,IMD,1)
      enddo
 421  enddo
 420  enddo
!
      if(ical_vof==1) then
        errrans(ivof)=errrans(ivof)+errk(ivof)
      endif
!
      if(NPE>1) then
         do  IMD=1,nrans
         CALL hpcrsum(errrans(IMD))
         CALL hpcrsum(errk(IMD))
        enddo
      endif
!
      do 425 IMD=1,nrans
      errk(IMD)=dsqrt(errk(IMD))/(dsqrt(errrans(IMD))+SML)
  425 enddo
!
!
      if(ieul2ph>0) then
      endif
!
      if(ical_vof==1) then
      endif
!
      if(ical_FC==PEFC) then
        akstmp(:)=aks(:,ical_s,1)
        call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,akstmp)
        aks(:,ical_s,1)=akstmp(:)
      endif
!
!
      if(ical_cavi==1) then
      endif
!
      if(icaltb==ke.or.
     &   icaltb==ke_low.or.icaltb==RNG.or.icaltb==CHEN.or.icaltb==KE2S) 
     &  then
        Do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) then
          do IMD=1,nrans
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          if((icaltb==ke.or.
     &        icaltb==ke_low.or.
     &        icaltb==RNG.or.
     &        icaltb==KE2S.or.
     &        icaltb==CHEN).and.
     &        (IMD.eq.ike(1).or.IMD.eq.ike(2))) then
            aks(ICVS:ICVE,ike(1),1)=0.d0
            aks(ICVS:ICVE,ike(2),1)=0.d0
          endif
          enddo
!
          if(u_func(1)==1) then
            aks(ICVS:ICVE,ike_c,1)=0.d0
          endif
        endif
        enddo
      endif
!
      if(ieul2ph>0) then
      endif
!---------------------------------------
! --- 
!---------------------------------------
      if(calxi) then   
        if(idifflm/=FltRans) then
          akstmp(:)=aks(:,ixi,1)
          call grad_cell(1,14,
     &    MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &    LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,akstmp,grdc(:,:,1))
          if(comb_mdl==1) then
            Cz=2.d0
          else
            Cz=1.d0
          endif
          do IIMAT=1,NMAT 
          IMAT=MAT_NO(IIMAT)
          if(IMAT<0) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
          dum3=grdc(ICVL,1,1)**2+grdc(ICVL,2,1)**2+grdc(ICVL,3,1)**2
          dum4=Cz*CVVOLM(ICVL)**(2.d0/3.d0)
          aks(ICVL,ixi_2,1)=dum3*dum4!*rho(ICVL,1)
          enddo
          enddo
        endif
!-------------------------------
! --- Scalar dissipation rate ixi_X
!-------------------------------ike(1)
        if(idifflm==FltRans) then
          do IIMAT=1,NMAT 
          IMAT=MAT_NO(IIMAT)
          if(IMAT<0) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
          dum3=aks(ICVL,ike(2),1)/(aks(ICVL,ike(1),1)+SML)
          dum4=max(Cx*dum3*aks(ICVL,ixi_2,1),x3min)
!          aks(ICVL,ixi_X,1)=min(dum4,akshg(ixi_X))
          aks(ICVL,ixi_X,1)=min(dum4,x3max)

!          dum3=aks(ICVL,ixi,1)
!          dum4=min(aks(ICVL,ixi_2,1),dum3-dum3**2) 
!          aks(ICVL,ixi_2,1)=dum4
          enddo
          enddo
        else
          dum1=1.d0/sgmaks(ixi)
          dum2=1.d0/sgmaksm(ixi)
          do IIMAT=1,NMAT 
          IMAT=MAT_NO(IIMAT)
          if(IMAT<0) cycle
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
          dum4=rmu(ICVL)*dum2+rmut(ICVL)*dum1
          dum3=grdc(ICVL,1,1)**2+grdc(ICVL,2,1)**2+grdc(ICVL,3,1)**2
          aks(ICVL,ixi_X,1)=Cx*dum3*dum4
          enddo
          enddo
        endif
      endif

      DEALLOCATE (dkdt,diagk,dfk,dflxk,akstmp)
!
      
      return
!
 9999 continue
      if(my_rank.eq.ROOT) write(ifle,*) '(rans_admin)'
      ierror=1
!
      end subroutine rans_admin
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine acs_monitor
     & (msmpl,iter,itery,iterk,iterv,iterp,MNVX,
     &  repsy,aepsy,repsk,aepsk,repsv,aepsv,repsp,aepsp,
     &  erry,errk,errv,errp,errsdy,errsdk,errsdv,
     &  reps_FAI,aeps_FAI,err_FAI,
     &  reps_A,aeps_A,err_A,iter_FAI,iter_A,errsdA,errsdFAI,
     &  iter_POTN,reps_POTN,aeps_POTN,err_POTN,POTNAL,
     &  time,deltt,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &  MHD_A,MHD_FAI,
     &  CVVOLM,tmp,vel,yys,prs,aks,ccc,rho)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only     : ifli,ifll,ifle
      use module_simple,only : nsmpl,simeer
      use module_time,only   : steady,toltim
      use module_hpc,only    : MNVX_GL
      use module_scalar,only : icalke,icalrsm,icalaph,
     &                         ike,irsm,iaph,ides,
     &                         rns_scl,calgeq,calxi,ixi,igeq,
     &                         ical_cavi,icavi,ivold,ivof,ical_FC,
     &                         pot_scl,ical_s,ical_flmlt,ikles
      use module_material,only : ical_sld,rotati,ishaft,end,begin,
     &                           rot_ang,OMGA
      use module_boundary,only : kdsld,nbcnd,kdbcnd,boundName,MAT_BCIDX,
     &                           LBC_INDEX
      use module_species,only  : spcnam
      use module_metrix,only   : ymin=>rcomp,ymax=>ys
      use module_model,only    : ical_vect,nthrds,ical_MHD,ke,ke_low,
     &                           RNG,CHEN,SDES,icaltb,KE2S,ical_mvmsh,
     &                           ical_dens,KLES
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_VOF     ,only : ical_vof
      use module_deltat  ,only : maxcou,maxcou2
      use module_time,    only : i_steady
      use module_model,   only : vertex_cen,cell_cen,icon_cv,ical_dens
!
      implicit none
!
! --- [dummy arguments] 
!
      integer,intent(in) :: iter,msmpl,MNVX
      integer,intent(in) :: itery(0:ncomp),iterk(mxrans),iterv(3),iterp
      real*8 ,intent(in) :: repsy(0:ncomp),aepsy(0:ncomp)
      real*8 ,intent(in) :: repsk(mxrans),aepsk(mxrans)
      real*8 ,intent(in) :: repsv(3),aepsv(3),repsp,aepsp
      real*8 ,intent(in) :: erry(0:ncomp),errk(mxrans),errv(3),errp
      real*8 ,intent(in) :: errsdy(0:ncomp),errsdk(mxrans),errsdv(3)
      real*8 ,intent(in) :: time,deltt
      INTEGER,INTENT(IN) :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN) :: MAT_INDEX(   0:MXMAT)
      INTEGER,INTENT(IN) :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN) :: MAT_DCIDX(   0:MXMAT)
      INTEGER,INTENT(IN) :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN) :: mat_cal  (   0:MXMAT)
      real*8 ,intent(in) :: tmp  (       MXALLCV)
      real*8 ,intent(in) :: ccc  (       MXALLCV)
      real*8 ,intent(in) :: rho  (       MXALLCV)
      real*8 ,intent(in) :: CVVOLM(      MXALLCV)
      real*8 ,intent(in) :: vel  (       MXALLCV,3)
      real*8 ,intent(in) :: yys  (       MXALLCV,MXcomp)
      real*8 ,intent(in) :: prs  (       MXALLCV)
      real*8 ,intent(in) :: aks  (       MXALLCVR,MXrans)
      real*8 ,intent(in) :: MHD_A  (     MXCV_MHD,3,2,2)
      real*8 ,intent(in) :: MHD_FAI(     MXCV_MHD,2,2)
      real*8, intent(in) :: reps_FAI(2),
     &                      aeps_FAI(2),
     &                      err_FAI(2),
     &                      reps_A(3,2),
     &                      aeps_A(3,2),
     &                      err_A(3,2)
      real*8 ,intent(in) :: errsdA(3,2),errsdFAI(2)
      integer,INTENT(IN) :: iter_FAI(2),iter_A(3,2)
!
      real*8 ,intent(in) :: err_POTN(MXPOTN)
      real*8 ,intent(in) :: reps_POTN(MXPOTN),aeps_POTN(MXPOTN)
      integer,intent(in) :: iter_POTN(MXPOTN)
      REAL*8,INTENT(IN)  :: POTNAL   (       MXALLCVP,MXPOTN)
      
!
!     err_FAI,err_A
!
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: icvmxc,IVGL
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: IIMAT,ICVS,ICVE,ICVL,idumy=0,idum,IPOT
      real*8  :: tmin,tmax,pmin,pmax,umin(3),umax(3),dumy=0.D0,dum1
      real*8  :: AMHDmin(3,2),FMHDmin(2),AMHDmax(3,2),FMHDmax(2)
      real*8  :: velmgn,courant=0.d0,maxcoux=0.d0
      character(len=14) :: max_c,CV_type
      
      real*8  :: akmin(2),akmax(2)
      real*8  :: velloc(3),aksloc(2),prsloc,tmploc,yysloc
      real*8  :: alpha,rot
      character(len=6) :: opt,optmve
      real*8  :: ptnmin(8),ptnmax(8),ptnloc(8),err_P(8),reps_P(8),
     &           aeps_P(8),dden4
      integer :: iter_P(8)
!
      IF(ical_vect) then
         opt='VECTOR'
      else
         opt='SCALAR'
      endif
      if(icon_cv==vertex_cen) then
        CV_type='AT VERTEX   : '
      else
        CV_type='AT CELL     : '
      endif
      IF(ical_mvmsh/=0) then
         optmve='MOVING'
      else
         optmve='NOMOVE'
      endif
      velloc=0.d0
      aksloc=0.d0
      prsloc=0.d0
      tmploc=0.d0
      yysloc=0.d0
!
      tmin=GREAT
      pmin=GREAT
      
      tmax=-GREAT
      pmax=-GREAT
!
!
      do 30 i=1,3
      umin(i)=GREAT
      umax(i)=-GREAT
  30  continue
!
      IF(ical_MHD/=0) then
        do i=1,2
        do J=1,3
        AMHDmin(J,i)=GREAT
        AMHDmax(J,i)=-GREAT
        enddo
        FMHDmin(i)=GREAT
        FMHDmax(i)=-GREAT
        enddo
      endif
!
      do 10 ICOM=1,ncomp
      ymin(ICOM)=GREAT
      ymax(ICOM)=-GREAT
  10  continue
!
      do 20 IMD=1,2   !nrans
      akmin(IMD)=GREAT
      akmax(IMD)=-GREAT
  20  continue
!
      if(pot_scl) then
        ptnmin(:)=GREAT
        ptnmax(:)=-GREAT
        ptnloc(:)=0.d0
        do IPOT=1,min(8,NPOTN)
        reps_P(IPOT)=reps_POTN(IPOT)
        aeps_P(IPOT)=aeps_POTN(IPOT)
        iter_P(IPOT)=iter_POTN(IPOT)
        enddo
      endif
!
      if(ical_mvmsh==2.or.ical_mvmsh==3.or.ical_mvmsh==4) goto 1000
      IF(ical_vect) then
        DO ICVL=ICVS_V,ICVE_V
          tmin=min(tmin,tmp(ICVL))
          pmin=min(pmin,prs(ICVL))
          tmax=max(tmax,tmp(ICVL))
          pmax=max(pmax,prs(ICVL))
        enddo
!
        do i=1,3
        DO ICVL=ICVS_V,ICVE_V
        umin(i)=min(umin(i),vel(ICVL,i))
        umax(i)=max(umax(i),vel(ICVL,i))
        enddo
        enddo
!
        do ICOM=1,ncomp
        DO ICVL=ICVS_V,ICVE_V
        ymin(ICOM)=min(ymin(ICOM),yys(ICVL,ICOM))
        ymax(ICOM)=max(ymax(ICOM),yys(ICVL,ICOM))
        enddo
        enddo
!
        if(rns_scl) then
          if(icaltb==ke.or.
     &       icaltb==ke_low.or.
     &       icaltb==RNG.or.
     &       icaltb==KE2S.or.
     &       icaltb==CHEN)
     &    then
            do IMD=1,2
            DO ICVL=ICVS_V,ICVE_V
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ike(IMD)))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ike(IMD)))
            enddo
            enddo
          elseif(icaltb==SDES) then
            IMD=ides
            DO ICVL=ICVS_V,ICVE_V
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ides))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ides))
            enddo
          elseif(icaltb==KLES) then
            IMD=ikles
            DO ICVL=ICVS_V,ICVE_V
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ikles))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ikles))
            enddo
          endif
        endif
      else
        do 100 IIMAT=1,NMAT   !ICV=1,NCVIN
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        DO ICVL=ICVS,ICVE
        tmin=min(tmin,tmp(ICVL))
        if(IMAT>0) pmin=min(pmin,prs(ICVL))
        tmax=max(tmax,tmp(ICVL))
        if(IMAT>0) pmax=max(pmax,prs(ICVL))
        do 120 i=1,3
        umin(i)=min(umin(i),vel(ICVL,i))
        umax(i)=max(umax(i),vel(ICVL,i))
 120    enddo
        do 130 ICOM=1,ncomp
        ymin(ICOM)=min(ymin(ICOM),yys(ICVL,ICOM))
        ymax(ICOM)=max(ymax(ICOM),yys(ICVL,ICOM))
 130    enddo

        if(ical_dens==4) then
        endif
        if(rns_scl.and.ical_dens/=4) then
          if(icaltb==ke.or.
     &       icaltb==ke_low.or.
     &       icaltb==RNG.or.
     &       icaltb==KE2S.or.
     &       icaltb==CHEN)
     &    then
            do 140 IMD=1,2
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ike(IMD)))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ike(IMD)))
 140        continue
          elseif(icaltb==SDES) then
            IMD=1
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ides))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ides))
          elseif(icaltb==KLES) then
            IMD=1
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ikles))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ikles))
          elseif(ical_flmlt) then
            if(calgeq) then
              IMD=igeq
              akmin(IMD)=min(akmin(IMD),aks(ICVL,IMD))
              akmax(IMD)=max(akmax(IMD),aks(ICVL,IMD))
            elseif(calxi) then
              IMD=ixi
              akmin(IMD)=min(akmin(IMD),aks(ICVL,IMD))
              akmax(IMD)=max(akmax(IMD),aks(ICVL,IMD))
            elseif(calgeq.and.calxi) then
              IMD=1
              akmin(IMD)=min(akmin(IMD),aks(ICVL,ixi))
              akmax(IMD)=max(akmax(IMD),aks(ICVL,ixi))
              IMD=2
              akmin(IMD)=min(akmin(IMD),aks(ICVL,igeq))
              akmax(IMD)=max(akmax(IMD),aks(ICVL,igeq))
            endif
          elseif(ical_cavi==1) then 
            IMD=1
            akmin(IMD)=min(akmin(IMD),aks(ICVL,icavi))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,icavi))
            IMD=2
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ivold))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ivold))
          elseif(ical_FC>0) then
            IMD=1
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ical_s))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ical_s))
          elseif(ical_vof==1) then
            IMD=1
            akmin(IMD)=min(akmin(IMD),aks(ICVL,ivof))
            akmax(IMD)=max(akmax(IMD),aks(ICVL,ivof))
          endif
        endif
!
        if(pot_scl) then
          do IPOT=1,min(NPOTN,8)
          ptnmin(IPOT)=min(ptnmin(IPOT),POTNAL(ICVL,IPOT))
          ptnmax(IPOT)=max(ptnmax(IPOT),POTNAL(ICVL,IPOT))
          enddo
        endif
        enddo
 100    enddo
        IF(ICAL_MHD/=0) then
          do IIMAT=1,NMAT   !ICV=1,NCVIN
          IMAT=MAT_NO(IIMAT)
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          DO ICVL=ICVS,ICVE
          do j=1,2
          do i=1,3
          AMHDmin(i,j)=min(AMHDmin(i,j),MHD_A(ICVL,i,j,1))
          AMHDmax(i,j)=max(AMHDmax(i,j),MHD_A(ICVL,i,j,1))
          enddo
          FMHDmin(j)=min(FMHDmin(j),MHD_FAI(ICVL,j,1))
          FMHDmax(j)=max(FMHDmax(j),MHD_FAI(ICVL,j,1))
          enddo
          enddo
          enddo
        endif
      endif
!
      if(NPE.gt.1) then
        call hpcrmin_0(tmin)
        call hpcrmin_0(pmin)
        call hpcrmax_0(tmax)
        call hpcrmax_0(pmax)
        do i=1,3
        call hpcrmin_0(umin(i))
        call hpcrmax_0(umax(i))
        enddo
        do ICOM=1,ncomp
        call hpcrmin_0(ymin(ICOM))
        call hpcrmax_0(ymax(ICOM))
        enddo
        if(rns_scl) then
          do IMD=1,2    !nrans
          call hpcrmin_0(akmin(IMD))
          call hpcrmax_0(akmax(IMD))
          enddo
        endif
        if(pot_scl) then
          do IPOT=1,min(NPOTN,8)
          call hpcrmin_0(akmin(IPOT))
          call hpcrmax_0(akmax(IPOT))
          enddo
        endif
      endif
!
! --- 
!
!     if(.not.steady) then
        maxcoux=ZERO
        icvmxc=1
        do 400 IIMAT=1,NMAT   !iv=1,NCVIN
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT<0) cycle

        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
!
        dden4=0.d0
        if(ical_dens==4) then
          dden4=1.d0
        endif
        dden4=0.d0

        DO ICVL=ICVS,ICVE 
          velmgn=vel(ICVL,1)*vel(ICVL,1)+
     &           vel(ICVL,2)*vel(ICVL,2)+
     &           vel(ICVL,3)*vel(ICVL,3)
          velmgn=dsqrt(velmgn+dden4*CCC(ICVL))
          courant=velmgn*deltt/(CVVOLM(ICVL)**0.3333D0+SML)
!
          if(courant.gt.maxcoux) then
            maxcoux=courant
            icvmxc=ICVL
          endif
        enddo
 400    continue 
!
        if(maxcoux>maxcou) then  !max_c='MAX.COURANT'
          max_c='MAX.CV_COURANT'
          maxcou=maxcoux
        else
          max_c='MAX.SF_COURANT'
          !maxcou=maxcou
        endif
!     endif
!
! --- HPC: MNVX_GL will be read in namelist
!
      if(NPE.gt.1) then
        IVGL=NOD_IDg(MAT_CV(icvmxc))
        call hpcmaxloc(IVGL,maxcou)
        if(MNVX.gt.0) then
          do i=1,3
          velloc(i)=vel(MNVX,i)
          enddo
          prsloc=prs(MNVX)
          tmploc=tmp(MNVX)
!
          if(ical_dens==4) then
             aksloc(1)=rho(MNVX)
             aksloc(2)=sqrt(ccc(MNVX))
          endif
          if(rns_scl.and.ical_dens/=4) then
            if(icaltb==ke.or.
     &       icaltb==ke_low.or.
     &       icaltb==RNG.or.
     &       icaltb==KE2S.or.
     &       icaltb==CHEN)
     &      then
              aksloc(1)=aks(MNVX,ike(1))
              aksloc(2)=aks(MNVX,ike(2))
            elseif(icaltb==SDES) then
              aksloc(1)=aks(MNVX,ides)
              aksloc(2)=0.d0
            elseif(icaltb==KLES) then
              aksloc(1)=aks(MNVX,ikles)
              aksloc(2)=0.d0
            elseif(ical_flmlt) then
              if(calgeq) then
                aksloc(1)=aks(MNVX,igeq)
              elseif(calxi) then
                aksloc(1)=aks(MNVX,ixi)
              elseif(calgeq.and.calxi) then
                aksloc(1)=aks(MNVX,ixi)
                aksloc(2)=aks(MNVX,igeq)
              endif
            elseif(ical_cavi==1) then
              aksloc(1)=aks(MNVX,icavi)
              aksloc(2)=aks(MNVX,ivold)
            elseif(ical_FC>0) then
              aksloc(1)=aks(MNVX,ical_s)
            elseif(ical_vof==1) then
              aksloc(1)=aks(MNVX,ivof)
            endif
          endif
          yysloc=yys(MNVX,1)
          if(pot_scl) then
            do IPOT=1,min(NPOTN,8)
            ptnloc(IPOT)=POTNAL(MNVX,IPOT)
            enddo
          endif
        endif
        call hpcrsum_0(prsloc)
        call hpcrsum_0(tmploc)
        call hpcrsum_0(aksloc(1))
        call hpcrsum_0(aksloc(2))
        call hpcrsum_0(yysloc)
        call hpcrsum_0(velloc(1))
        call hpcrsum_0(velloc(2))
        call hpcrsum_0(velloc(3))
        if(pot_scl) then
          do IPOT=1,min(NPOTN,8)
          call hpcrsum(ptnloc(IPOT))
          enddo
        endif
      else
        IVGL=MAT_CV(icvmxc)
        do i=1,3
        velloc(i)=vel(MNVX,i)
        enddo
        prsloc=prs(MNVX)
        tmploc=tmp(MNVX)
        if(ical_dens==4) then
          aksloc(1)=rho(MNVX)
          aksloc(2)=sqrt(ccc(MNVX))
        endif

        if(rns_scl.and.ical_dens/=4) then
          if(icaltb==ke.or.
     &       icaltb==ke_low.or.
     &       icaltb==RNG.or.
     &       icaltb==KE2S.or.
     &       icaltb==CHEN)
     &    then
            aksloc(1)=aks(MNVX,ike(1))
            aksloc(2)=aks(MNVX,ike(2))
          elseif(icaltb==SDES) then
            aksloc(1)=aks(MNVX,ides)
            aksloc(2)=0.d0
          elseif(icaltb==KLES) then
            aksloc(1)=aks(MNVX,ikles)
            aksloc(2)=0.d0
          elseif(ical_flmlt) then
            if(calgeq) then
              aksloc(1)=aks(MNVX,igeq)
            elseif(calxi) then
              aksloc(1)=aks(MNVX,ixi)
            elseif(calgeq.and.calxi) then
              aksloc(1)=aks(MNVX,ixi)
              aksloc(2)=aks(MNVX,igeq)
            endif
          elseif(ical_cavi==1) then
            aksloc(1)=aks(MNVX,icavi)
            aksloc(2)=aks(MNVX,ivold)
          elseif(ical_FC>0) then
            aksloc(1)=aks(MNVX,ical_s)
          elseif(ical_vof==1) then
            aksloc(1)=aks(MNVX,ivof)
          endif
        endif
        yysloc=yys(MNVX,1)
        if(pot_scl) then
          do IPOT=1,min(NPOTN,8)
          ptnloc(IPOT)=POTNAL(MNVX,IPOT)
          enddo
        endif
      endif  
!
! --- OUTPUT:
!=====================================================================
 1000 continue
      if(my_rank.eq.ROOT) then
        write(ifll,3001) iter,opt,optmve,
     &     NPE,nthrds,i_steady,
     &    time,deltt,max_c,maxcou,CV_type,IVGL
        if(ical_dens==4) then
          write(ifll,3009)
        endif
        if(rns_scl.and.ical_dens/=4) then
          if(icaltb==ke.or.
     &       icaltb==ke_low.or.
     &       icaltb==RNG.or.
     &       icaltb==KE2S.or.
     &       icaltb==CHEN)
     &    then
            write(ifll,3005)
          elseif(icaltb==SDES.or.icaltb==KLES) then
            write(ifll,3005)
          elseif(ical_flmlt) then
            write(ifll,3004)
          elseif(ical_cavi==1) then
            write(ifll,3006)
          elseif(ical_FC>0) then
            write(ifll,3007)
          elseif(ical_vof==1) then
            write(ifll,3008)
          else
            write(ifll,3005)
          endif
        elseif(ical_dens/=4) then
           write(ifll,3005)
        endif
!
        write(ifll,3020)
!
        write(ifll,3035) (velloc(i),i=1,3),prsloc,tmploc,
     &            (aksloc(i),i=1,2),yysloc
!
        write(ifll,3040) (umin(i),i=1,3),pmin,tmin,
     &            (akmin(i),i=1,2),(ymin(i),i=1,1)
!
        write(ifll,3042) (umax(i),i=1,3),pmax,tmax,
     &            (akmax(i),i=1,2),(ymax(i),i=1,1)
!
        write(ifll,3020)
        write(ifll,3030) msmpl
        write(ifll,3036) 
        write(ifll,3032) 
!
        if(mxrans.ge.2) then
          idum=0
          do i=1,ncomp
            if(itery(i)>idum) idum=itery(i)
          enddo
          
          write(ifll,3044) (iterv(i),i=1,3),iterp,itery(0),
     &    (iterk(i),i=1,2),idum!(itery(i),i=1,1)
!
          dum1=-1.d0
          do i=1,ncomp
            if(repsy(i)>dum1) dum1=repsy(i)
          enddo
          write(ifll,3046) (repsv(i),i=1,3),repsp,repsy(0),
     &    (repsk(i),i=1,2),dum1!(repsy(i),i=1,1)
!
          dum1=-1.d0
          do i=1,ncomp
            if(aepsy(i)>dum1) dum1=aepsy(i)
          enddo
          write(ifll,3048) (aepsv(i),i=1,3),aepsp,aepsy(0),
     &    (aepsk(i),i=1,2),dum1!(aepsy(i),i=1,1)
!
          write(ifll,3020)
          write(ifll,3049) simeer
!
          IF(msmpl.GE.nsmpl.and.nsmpl.ge.2) write(ifll,3054)
!
          dum1=-1.d0
          do i=1,ncomp
            if(erry(i)>dum1) dum1=erry(i)
          enddo
          write(ifll,3050) (errv(i),i=1,3),errp,erry(0),
     &    (errk(i),i=1,2),dum1!(erry(i),i=1,1)
!
          if(steady) then
            dum1=-1.d0
            do i=1,ncomp
              if(errsdy(i)>dum1) dum1=errsdy(i)
            enddo
            write(ifll,3020)
            write(ifll,3060) toltim
            write(ifll,3058) (errsdv(i),i=1,3),errsdy(0),
     &     (errsdk(i),i=1,2),dum1!(errsdy(i),i=1,1)
          endif
        elseif(mxrans.le.1) then
          idum=0
          do i=1,ncomp
            if(itery(i)>idum) idum=itery(i)

          enddo
          write(ifll,3044) (iterv(i),i=1,3),iterp,itery(0),
     &    (iterk(i),i=1,1),idumy,idum!(itery(i),i=1,1)

!
          dum1=-1.d0
          do i=1,ncomp
            if(repsy(i)>dum1) dum1=repsy(i)
          enddo
          write(ifll,3046) (repsv(i),i=1,3),repsp,repsy(0),
     &    (repsk(i),i=1,1),dumy,dum1!(repsy(i),i=1,1)

!
          dum1=-1.d0
          do i=1,ncomp
            if(aepsy(i)>dum1) dum1=aepsy(i)
          enddo
          write(ifll,3048) (aepsv(i),i=1,3),aepsp,aepsy(0),
     &    (aepsk(i),i=1,1),dumy,dum1!(aepsy(i),i=1,1)
!
          write(ifll,3020)
          write(ifll,3049) simeer
!
          IF(msmpl.GE.nsmpl.and.nsmpl.ge.2) write(ifll,3054)

          dum1=-1.d0
          do i=1,ncomp
            if(erry(i)>dum1) dum1=erry(i)
          enddo
          write(ifll,3050) (errv(i),i=1,3),errp,erry(0),
     &    (errk(i),i=1,1),dumy,dum1!(erry(i),i=1,1)
!
          if(steady) then
            dum1=-1.d0
            do i=1,ncomp
              if(errsdy(i)>dum1) dum1=errsdy(i)
            enddo
            write(ifll,3020)
            write(ifll,3060) toltim
            dum1=-1.d0
            do i=1,ncomp
              if(errsdy(i)>dum1) dum1=errsdy(i)
            enddo
            write(ifll,3058) (errsdv(i),i=1,3),errsdy(0),
     &      (errsdk(i),i=1,1),dumy,dum1!(errsdy(i),i=1,1)
          endif
        endif
!
        write(ifll,4020)
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          IF(IMAT.lt.0) cycle
          if(ishaft(IMAT)==1) then
            alpha=rot_ang(IMAT)*180.d0/3.14159d0
            rot=rotati(IMAT)*60.d0/(2.d0*3.1415926d0)
            write(ifll,3100) 
            write(ifll,3110) IMAT,rot,alpha
          endif
        enddo
!
        IF(ical_MHD/=0) then
          write(ifll,4021)
          write(ifll,4005)   
       
          write(ifll,3040) ((AMHDmin(i,j),i=1,3),j=1,2),
     &            FMHDmin(1),FMHDmin(2)
!
          write(ifll,3042) ((AMHDmax(i,j),i=1,3),j=1,2),
     &            FMHDmax(1),FMHDmax(2)
          dum1=time*OMGA(1)*180.d0/3.1415926d0
          write(ifll,4032) dum1
          write(ifll,4044) (iter_A(i,1),i=1,3),
     &    (iter_A(i,2),i=1,3),iter_FAI(1),iter_FAI(2)

          write(ifll,4046) (reps_A(i,1),i=1,3),
     &    (reps_A(i,2),i=1,3),reps_FAI(1),reps_FAI(2)

          write(ifll,4048) (aeps_A(i,1),i=1,3),
     &    (aeps_A(i,2),i=1,3),aeps_FAI(1),aeps_FAI(1)

          write(ifll,4049) 

          write(ifll,4050) (err_A(i,1),i=1,3),
     &    (err_A(i,2),i=1,3),err_FAI(1),err_FAI(2)

          IF(steady) then
            write(ifll,4060) 
            write(ifll,4058) (errsdA(i,1),i=1,3),
     &      (errsdA(i,2),i=1,3)
          endif
        endif


        IF(pot_scl) then
          
          write(ifll,4021)
          write(ifll,4006)
          write(ifll,3020)
       
          write(ifll,3035) (ptnloc(I),I=1,8)
          write(ifll,3040) (ptnmin(I),I=1,8)
          write(ifll,3042) (ptnmax(I),I=1,8)
          write(ifll,3036) 

          write(ifll,4044) (iter_P(i),i=1,8)
          write(ifll,4046) (reps_P(i),i=1,8)
          write(ifll,4048) (aeps_P(i),i=1,8)
          write(ifll,4020)
        endif
        write(ifll,4020)
      endif
!
      return
!
!=====================================================================
 3001 FORMAT(
     1 /,2X,'|',108('='),'|',
     2 /,2X,'|',2X,'TIME STEP   : ',I10,6x,a6,6X,a6,6X,
     2 'MPI: ',I6,2X,'openMP: ',I6,3x,'FLOW TYPE:',I3,9X,'|',
     3 /,2X,'|',2X,'TIME : ',E10.4,4X,
     3            'TIME INTERVAL: ',E10.4,4X,
     3             a14,': ',E10.4,2X,
     3             a14,I10,4X,'|',
     4 /,2X,'|',108('-'),'|')
 3002 FORMAT(
     1 /,2X,'|',108('='),'|',
     2 /,2X,'|',2X,'TIME STEP   : ',I10,82X,'|',
     3 /,2X,'|',2X,'TIME : ',E10.4,4X,
     3            'TIME INTERVAL: ',E10.4,64X,'|',
     4 /,2X,'|',108('-'),'|')

 3004 FORMAT(2X,'|',2X,'Velocity Field :',
     &  1X,'U',11X,'V',11X,'W',11X,'P',11X,'T',
     & 11X,'Z',9X,'EPS',6X,'COMP_1',4X,'|')
 3005 FORMAT(2X,'|',2X,'Velocity Field :',
     &  1X,'U',11X,'V',11X,'W',11X,'P',11X,'T',
     & 11X,'K',9X,'EPS',6X,'COMP_1',4X,'|')
 3006 FORMAT(2X,'|',2X,'Velocity Field :',
     &  1X,'U',11X,'V',11X,'W',11X,'P',11X,'T',
     & 11X,'Y',8X,'VOID',6X,'COMP_1',4X,'|')
 3007 FORMAT(2X,'|',2X,'Velocity Field :',
     &  1X,'U',11X,'V',11X,'W',11X,'P',11X,'T',
     & 11X,'S',9X,'EPS',6X,'COMP_1',4X,'|')
 3008 FORMAT(2X,'|',2X,'Velocity Field :',
     &  1X,'U',11X,'V',11X,'W',11X,'P',11X,'T',
     & 11X,'VOF',7X,'EPS',6X,'COMP_1',4X,'|')
 3009 FORMAT(2X,'|',2X,'Velocity Field :',
     &  1X,'U',11X,'V',11X,'W',11X,'P',11X,'T',
     & 11X,'RHO',7X,' Cs',6X,'COMP_1',4X,'|')

 3010 FORMAT(2X,'|',8X,8(2X,F10.3),4X,'|')
 3020 FORMAT(2X,'|',108('-'),'|')
 3025 FORMAT(2X,'|',108('.'),'|')
! 3030 FORMAT(2X,'|',2X,'OUTER-ITER NO. : ',I10,79X,'|')
 3030 FORMAT(2X,'|',2X,'OUTER-ITER NO. : ',I10,67X,
     &       'ALL COMP',4X,'|')
 3032 FORMAT(2X,'|',2X,'CG/ICCG :',97X,'|')
 3035 FORMAT(2X,'|',2X,'MONI :',8(2X,E10.4),4X,'|')
 3036 FORMAT(2X,'|',108X,'|')
 3040 FORMAT(2X,'|',2X,'MIN  :',8(2X,E10.4),4X,'|')
 3042 FORMAT(2X,'|',2X,'MAX  :',8(2X,E10.4),4X,'|')
 3044 FORMAT(2X,'|',2X,'ITER :',8(2X,I10)  ,4X,'|')
 3046 FORMAT(2X,'|',2X,'RES  :',8(2X,E10.4),4X,'|')
 3048 FORMAT(2X,'|',2X,'ABS  :',8(2X,E10.4),4X,'|')
 3049 FORMAT(2X,'|',2X,'OUTER-ITER ERROR :',20X,
     &       'CNVERGENCE TOLERANCE :',
     &       4X,E10.4,32X,'|')
 3050 FORMAT(2X,'|',2X,'ERR  :',8(2X,E10.4),4X,'|')
 3054 FORMAT(2X,'|',2X,
     &     ' ???WARNING: OUTER-ITER GREAT THEN USER DEFINED ',58X,'|')
 3058 FORMAT(2X,'|',2X,'ERR  :',3(2X,E10.4),3X,6('-'),
     &3X,4(2X,E10.4),4X,'|')
 3060 FORMAT(2X,'|',2X,'STEADY FLOW ERROR :',19X,
     &       'CNVERGENCE TOLERANCE :',
     &       4X,E10.4,32X,'|')
 3100 FORMAT(2X,'|',2X,'ROTATION MATERIAL INFO : ',81X,'|')
! 3110 FORMAT(2X,'|',2X,'MATRIAL NUMBER :',I4,
!     &      ' CURRENT ANGLE:',F12.5,2X,'<=>',2x,'MATRIAL NUMBER :',I4,
!     &       ' CURRENT ANGLE:',F12.5,5X,'|')
 3110 FORMAT(2X,'|',2X,'MATRIAL NUMBER :',I4,4x,'|',
     &       ' CURRENT RPM:',F12.5,4X,' CURRENT ANGLE:',F12.5,25X,'|')
 3120 FORMAT(2X,'|',2X,'SLIDING BC NAME ',A,1X,'&',1X,A,63x,'|')
 3130 FORMAT(2X,'|',2X,'MATRIAL NUMBER :',I4,
     &       ' CURRENT RPM  :',F12.5,2X,'<=>',2x,'MATRIAL NUMBER :',I4,
     &       ' CURRENT RPM  :',F12.5,5X,'|')
 4005 FORMAT(2X,'|',2X,'MHD FIELD :',4X,'Ar1',9X,'Ar2',9X,'Ar3',9X,
     &   'Ai1',9X,'Ai2',
     & 9X,'Ai3',8X,'FAIr',8X,'FAIi',4X,'|')
 4006 FORMAT(2X,'|',2X,'POTENTIAL :',4X,'PT1',9X,'PT2',9X,'PT3',9X,
     &   'PT4',9X,'PT5',
     & 9X,'PT6',9X,'PT7',9X,'PT8',4X,'|')
 4021 FORMAT(2X,'|',108('~'),'|')
 4032 FORMAT(2X,'|',2X,'CG/ICCG :',23X,'OMEGA*T=',(1X,E10.4),55X,'|')
 4044 FORMAT(2X,'|',2X,'ITER :',8(2X,I10)  ,4X,'|')
 4046 FORMAT(2X,'|',2X,'RES  :',8(2X,E10.4),4X,'|')
 4048 FORMAT(2X,'|',2X,'ABS  :',8(2X,E10.4),4X,'|')
 4049 FORMAT(2X,'|',2X,'OUTER-ITER EEROR :',88X,'|')
 4050 FORMAT(2X,'|',2X,'ERR  :',8(2X,E10.4),4X,'|')
 4058 FORMAT(2X,'|',2X,'ERR  :',6(2X,E10.4),6X,6('-'),
     &       6X,6('-'),4X,'|')
 4060 FORMAT(2X,'|',2X,'STEADY FLOW ERROR :',87X,'|')
 4020 FORMAT(2X,'|',108('='),'|')
!
!=====================================================================
      end subroutine acs_monitor
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine acs_monitor2
     &    (time,deltt,itery2,iterk,iterv2,
     &    iterp2,repsp2,aepsp2,errp2,
     &    repsy2,aepsy2,repsk,aepsk,repsv2,aepsv2,
     &    erry2,errk,errv2,errsdy2,errsdk,errsdv2,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &    CVVOLM,tmp2,vel2,yys2,prs2,aks)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_io,only     : ifli,ifll,ifle
      use module_simple,only : nsmpl
      use module_time,only   : steady,toltim
      use module_hpc,only    : MNVX_GL
      use module_scalar,only : icalke,icalrsm,icalaph,
     &                         ike,irsm,iaph,ides,
     &                         rns_scl
      use module_material,only : ical_sld,rotati,ishaft,end,begin,
     &                           rot_ang,OMGA
      use module_boundary,only : kdsld,nbcnd,kdbcnd,boundName,MAT_BCIDX,
     &                           LBC_INDEX
      use module_species,only  : spcnam
      use module_metrix,only   : ymin=>rcomp,ymax=>ys
      use module_model,only    : ical_vect,nthrds,ical_MHD,ke,ke_low,
     &                           RNG,CHEN,SDES,icaltb,ical_mvmsh
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_deltat  ,only : maxcou,maxcou2
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: iterp2
      real*8 ,intent(in) :: repsp2,aepsp2,errp2
      real*8 ,intent(in) :: time,deltt
      integer,intent(in) :: itery2(0:ncomp),iterk(mxrans),iterv2(3)
      real*8 ,intent(in) :: repsy2(0:ncomp),aepsy2(0:ncomp)
      real*8 ,intent(in) :: repsk(mxrans),aepsk(mxrans)
      real*8 ,intent(in) :: repsv2(3),aepsv2(3)
      real*8 ,intent(in) :: erry2(0:ncomp),errk(mxrans),errv2(3)
      real*8 ,intent(in) :: errsdy2(0:ncomp),errsdk(mxrans),errsdv2(3)
      INTEGER,INTENT(IN) :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN) :: MAT_INDEX(   0:MXMAT)
      INTEGER,INTENT(IN) :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN) :: MAT_DCIDX(   0:MXMAT)
      INTEGER,INTENT(IN) :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN) :: mat_cal  (   0:MXMAT)
      real*8 ,intent(in) :: CVVOLM(      MXALLCV)
      real*8 ,intent(in) :: tmp2 (       MXALLCV2)
      real*8 ,intent(in) :: prs2 (       MXALLCV2)
      real*8 ,intent(in) :: vel2 (       MXALLCV2,3)
      real*8 ,intent(in) :: yys2 (       MXALLCV2,MXcomp)
      real*8 ,intent(in) :: aks  (       MXALLCVR,MXrans)
!
!     err_FAI,err_A
!
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: icvmxc,IVGL
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: IIMAT,ICVS,ICVE,ICVL,idumy=0,idum=0
      real*8  :: tmin,tmax,pmin,pmax,umin(3),umax(3),dumy=0.D0,dum1
      real*8  :: velmgn,courant=0.d0
      real*8  :: akmin(2),akmax(2)
      real*8  :: velloc(3),aksloc(2),prsloc,tmploc,yysloc
      real*8  :: alpha,rot
      character(len=6) :: opt='E2P'
!
      velloc=0.d0
      aksloc=0.d0
      prsloc=0.d0
      tmploc=0.d0
      yysloc=0.d0
!
      tmin=GREAT
      pmin=GREAT
!      
      tmax=-GREAT
      pmax=-GREAT
!
      do 30 i=1,3
      umin(i)=GREAT
      umax(i)=-GREAT
  30  continue
!
      do 10 ICOM=1,ncomp
      ymin(ICOM)=GREAT
      ymax(ICOM)=-GREAT
  10  continue
!
      do 20 IMD=1,2   !nrans
      akmin(IMD)=GREAT
      akmax(IMD)=-GREAT
  20  continue
!
      IF(ical_vect) then
        call FFRABORT(1,'ERR: VECTOR Ver. NOT SUPPORT E2P')
      else
        do 100 IIMAT=1,NMAT   !ICV=1,NCVIN
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        DO ICVL=ICVS,ICVE
        tmin=min(tmin,tmp2(ICVL))
        pmin=min(pmin,prs2(ICVL))
        tmax=max(tmax,tmp2(ICVL))
        pmax=max(pmax,prs2(ICVL))
        do 120 i=1,3
        umin(i)=min(umin(i),vel2(ICVL,i))
        umax(i)=max(umax(i),vel2(ICVL,i))
 120    enddo
        do 130 ICOM=1,ncomp
        ymin(ICOM)=min(ymin(ICOM),yys2(ICVL,ICOM))
        ymax(ICOM)=max(ymax(ICOM),yys2(ICVL,ICOM))
 130    enddo
        if(rns_scl) then
          akmin(1)=min(akmin(1),aks(ICVL,iaph(1)))
          akmax(1)=max(akmax(1),aks(ICVL,iaph(1)))
          akmin(2)=min(akmin(2),aks(ICVL,iaph(2)))
          akmax(2)=max(akmax(2),aks(ICVL,iaph(2)))
        endif
        enddo
 100    enddo
      endif
!
      if(NPE.gt.1) then
        call hpcrmin_0(tmin)
        call hpcrmin_0(pmin)
        call hpcrmax_0(tmax)
        call hpcrmax_0(pmax)
        do i=1,3
        call hpcrmin_0(umin(i))
        call hpcrmax_0(umax(i))
        enddo
        do ICOM=1,ncomp
        call hpcrmin_0(ymin(ICOM))
        call hpcrmax_0(ymax(ICOM))
        enddo
        do IMD=1,2    !nrans
        call hpcrmin_0(akmin(IMD))
        call hpcrmax_0(akmax(IMD))
        enddo
      endif
!
! --- 
!
!     if(.not.steady) then
!        maxcou=ZERO
!        icvmxc=1
!        do 400 IIMAT=1,NMAT   !iv=1,NCVIN
!        if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        DO ICVL=ICVS,ICVE
!        velmgn=vel2(ICVL,1)*vel2(ICVL,1)+
!     &         vel2(ICVL,2)*vel2(ICVL,2)+
!     &         vel2(ICVL,3)*vel2(ICVL,3)
!        velmgn=dsqrt(velmgn)
!        courant=velmgn*deltt/(CVVOLM(ICVL)**0.3333D0+SML)
!        if(courant.gt.maxcou) then
!          maxcou=courant
!          icvmxc=ICVL
!        endif
!        enddo
! 400    continue
!     endif
!
! --- HPC: MNVX_GL will be read in namelist
!
!=====================================================================
 1000 continue
      if(my_rank.eq.ROOT) then
        write(ifll,3001) 
        write(ifll,3005)
        write(ifll,3020)
!
!        write(ifll,3035) (velloc(i),i=1,3),prsloc,tmploc,
!     &            (aksloc(i),i=1,2),yysloc
!
        write(ifll,3040) (umin(i),i=1,3),pmin,tmin,
     &            (akmin(i),i=1,2),(ymin(i),i=1,1)
!
        write(ifll,3042) (umax(i),i=1,3),pmax,tmax,
     &            (akmax(i),i=1,2),(ymax(i),i=1,1)
!
        write(ifll,3020)
!        write(ifll,3036) 
!        write(ifll,3032) 
!
        idum=0
        do i=1,ncomp
          if(itery2(i)>idum) idum=itery2(i)
        enddo
        write(ifll,3044) (iterv2(i),i=1,3),iterp2,itery2(0),
     &    (iterk(iaph(i)),i=1,2),idum
!
        dum1=-1.d0
        do i=1,ncomp
          if(repsy2(i)>dum1) dum1=repsy2(i)
        enddo
        write(ifll,3046) (repsv2(i),i=1,3),repsp2,repsy2(0),
     &    (repsk(i),i=1,2),dum1
!
        dum1=-1.d0
        do i=1,ncomp
          if(aepsy2(i)>dum1) dum1=aepsy2(i)
        enddo
        write(ifll,3048) (aepsv2(i),i=1,3),aepsp2,aepsy2(0),
     &    (aepsk(iaph(i)),i=1,2),dum1
!
        write(ifll,3020)
!
        dum1=-1.d0
        do i=1,ncomp
          if(erry2(i)>dum1) dum1=erry2(i)
        enddo
        write(ifll,3050) (errv2(i),i=1,3),errp2,erry2(0),
     &    (errk(iaph(i)),i=1,2),dum1
!
        if(steady) then
          dum1=-1.d0
          do i=1,ncomp
            if(errsdy2(i)>dum1) dum1=errsdy2(i)
          enddo
          write(ifll,3020)
          write(ifll,3060) toltim
          write(ifll,3058) (errsdv2(i),i=1,3),errsdy2(0),
     &     (errsdk(iaph(i)),i=1,2),dum1
        endif
        write(ifll,4020)
      endif
!
      return
!
!=====================================================================
 3001 FORMAT(
!     1 /,2X,'|',108('='),'|',
     2 2X,'|',2X,'E2P',103x,'|',
     4 /,2X,'|',108('-'),'|')
 3005 FORMAT(2X,'|',2X,'Velocity Field :',
     &  1X,'U',11X,'V',11X,'W',11X,'P',11X,'T',
     & 9X,'AP1',9X,'AP2',6X,'COMP_1',4X,'|')
 3010 FORMAT(2X,'|',8X,8(2X,F10.3),4X,'|')
 3020 FORMAT(2X,'|',108('-'),'|')
 3025 FORMAT(2X,'|',108('.'),'|')
 3030 FORMAT(2X,'|',2X,'OUTER-ITER NO. : ',I10,79X,'|')
 3032 FORMAT(2X,'|',2X,'CG/ICCG :',97X,'|')
! 3035 FORMAT(2X,'|',2X,'MONI :',8(2X,E10.4),4X,'|')
 3036 FORMAT(2X,'|',108X,'|')
 3040 FORMAT(2X,'|',2X,'MIN  :',8(2X,E10.4),4X,'|')
 3042 FORMAT(2X,'|',2X,'MAX  :',8(2X,E10.4),4X,'|')
 3044 FORMAT(2X,'|',2X,'ITER :',8(2X,I10)  ,4X,'|')
 3046 FORMAT(2X,'|',2X,'RES  :',8(2X,E10.4),4X,'|')
 3048 FORMAT(2X,'|',2X,'ABS  :',8(2X,E10.4),4X,'|')
 3049 FORMAT(2X,'|',2X,'OUTER-ITER ERROR :',20X,
     &       'CNVERGENCE TOLERANCE :',
     &       4X,E10.4,32X,'|')
 3050 FORMAT(2X,'|',2X,'ERR  :',8(2X,E10.4),4X,'|')
 3054 FORMAT(2X,'|',2X,
     &     ' ???WARNING: OUTER-ITER GREAT THEN USER DEFINED ',58X,'|')
 3058 FORMAT(2X,'|',2X,'ERR  :',3(2X,E10.4),3X,6('-'),
     &3X,4(2X,E10.4),4X,'|')
 3060 FORMAT(2X,'|',2X,'STEADY FLOW ERROR :',19X,
     &       'CNVERGENCE TOLERANCE :',
     &       4X,E10.4,32X,'|')
 3100 FORMAT(2X,'|',2X,'ROTATION MATERIAL INFO : ',81X,'|')
! 3110 FORMAT(2X,'|',2X,'MATRIAL NUMBER :',I4,
!     &      ' CURRENT ANGLE:',F12.5,2X,'<=>',2x,'MATRIAL NUMBER :',I4,
!     &       ' CURRENT ANGLE:',F12.5,5X,'|')
 3110 FORMAT(2X,'|',2X,'MATRIAL NUMBER :',I4,4x,'|',
     &       ' CURRENT RPM:',F12.5,4X,' CURRENT ANGLE:',F12.5,25X,'|')
 3120 FORMAT(2X,'|',2X,'SLIDING BC NAME ',A,1X,'&',1X,A,63x,'|')
 3130 FORMAT(2X,'|',2X,'MATRIAL NUMBER :',I4,
     &       ' CURRENT RPM  :',F12.5,2X,'<=>',2x,'MATRIAL NUMBER :',I4,
     &       ' CURRENT RPM  :',F12.5,5X,'|')
 4021 FORMAT(2X,'|',108('~'),'|')
 4032 FORMAT(2X,'|',2X,'CG/ICCG :',23X,'OMEGA*T=',(1X,E10.4),55X,'|')
 4044 FORMAT(2X,'|',2X,'ITER :',8(2X,I10)  ,4X,'|')
 4046 FORMAT(2X,'|',2X,'RES  :',8(2X,E10.4),4X,'|')
 4048 FORMAT(2X,'|',2X,'ABS  :',8(2X,E10.4),4X,'|')
 4049 FORMAT(2X,'|',2X,'OUTER-ITER EEROR :',88X,'|')
 4050 FORMAT(2X,'|',2X,'ERR  :',8(2X,E10.4),4X,'|')
 4058 FORMAT(2X,'|',2X,'ERR  :',6(2X,E10.4),6X,6('-'),
     &       6X,6('-'),4X,'|')
 4060 FORMAT(2X,'|',2X,'STEADY FLOW ERROR :',87X,'|')
 4020 FORMAT(2X,'|',108('='),'|')
!
!=====================================================================
      end subroutine acs_monitor2

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine error_smac(iter,iters,
     &                      ncomp,mxrans,nrans,erry,errk,errv,errp,
     &                      itery,iterk,
     &                      err_FAI,err_A,
     &                      erry2,errv2,
     &                      cvgnc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
      use module_simple, only : simeer
      use module_model,  only : ical_MHD
      use module_Euler2ph,only : ieul2ph
      use module_scalar,only : icavi,ivold,ical_cavi
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)  :: iter,iters,ncomp,nrans,mxrans
      integer,intent(in)  :: itery(0:ncomp),iterk(mxrans)
      real*8, intent(in)  :: erry(0:ncomp),errk(mxrans),errv(3),errp
      real*8, intent(in)  :: erry2(0:ncomp),errv2(3)
      real*8, intent(in)  :: err_FAI(2),err_A(3,2)
      logical,intent(out) :: cvgnc
!
! --- [local entities]
!
      real*8 :: errmax
      integer:: icom,IMD,i,j
      integer,save :: i_e2p=0
      real*8,parameter  :: err_stand=1.d3
      integer,parameter :: iter_stand=5,iter_E2P=3
!     
      cvgnc=.false.
      errmax=0.d0
!
!      do 100 icom=0,ncomp
!      if(itery(icom)>0) then
!        errmax=max(errmax,erry(icom))
!      endif
! 100  continue

      if(nrans>0) then
        do 120 IMD=1,nrans
        if(ical_cavi==1.and.(IMD==icavi.or.IMD==ivold)) cycle
        errmax=max(errmax,errk(IMD))
 120    continue
      endif
!
      do 140 i=1,3
      errmax=max(errmax,errv(i))
 140  continue
      errmax=max(errmax,errp)
!
      if(ieul2ph>0) then
      endif
!
      IF(ical_MHD/=0) then
        do i=1,2
        errmax=max(errmax,err_FAI(i))
        do j=1,3
        errmax=max(errmax,err_A(j,i))
        enddo
        enddo
      endif
!
      if(NPE.gt.1) then
        call hpcrmax(errmax)
      endif
!
      if(ieul2ph>0) then
      else
        if(errmax>err_stand.and.iter>(iters+iter_stand)) then
          call FFRABORT(1,'SMAC/SIMPLE NOT CONVERGED')
        endif
      endif
!
      if(errmax.lt.simeer) then
        cvgnc=.true.
      else
        cvgnc=.false.
      endif
!
      if(NPE.gt.1) then
        call hpcland(cvgnc,2)
      endif
!
      return
!
      end subroutine error_smac


!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine error_stp(
     &       deltt,
     &       errsdy,errsdk,errsdv,errsdA,errsdFAI,
     &       MAT_CV,MAT_CVEXT,MAT_INDEX,MAT_DCIDX,MAT_NO,mat_cal,
     &       tmp,yys,vel,aks,MHD_A,MHD_FAI,errabA,
     &       stpcvg)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
      use module_scalar,  only : icalke,icalrsm,icalaph,
     &                           ike,irsm,iaph,
     &                           rns_scl
      use module_model,   only : ical_vect,ical_mhd
!
      use module_time,  only   : toltim
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
!
      implicit none
! --- [dummy arguments]
!

      real*8, intent(in)  :: deltt
      INTEGER,INTENT(IN)  :: MAT_CV   ( MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_DCIDX(0:MXMAT)
      real*8, intent(in)  :: tmp(       MXALLCV,2)
      real*8, intent(in)  :: yys(       MXALLCV,MXcomp,2)
      real*8, intent(in)  :: vel(       MXALLCV,3,2)
      real*8, intent(in)  :: aks(       MXALLCVR,MXrans,2)
      real*8 ,intent(in)  :: MHD_A     (MXCV_MHD,3,2,2)
      real*8 ,intent(in)  :: MHD_FAI   (MXCV_MHD,2,2)
      real*8, intent(out) :: errsdy(0:ncomp),errsdk(mxrans),errsdv(3)
      real*8, intent(inout) :: errsdA(3,2),errsdFAI(2)
      real*8 ,intent(in)  :: errabA(3,2)
      logical,intent(out) :: stpcvg
!
! --- [local entities]
!
      real*8 :: errmax
      integer:: IMD,i,ICV,ICOM,IMHD
      integer:: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      real*8 :: erraby(0:ncomp),errabk(mxrans),errabv(3)
      real*8 :: courant
!

      erraby=0.d0
      errabk=0.d0
      errabv=0.d0
      errsdy=0.d0
      errsdk=0.d0
      errsdv=0.d0
!      errsdA=0.d0
!     
! --- 
!
      IF(ical_vect) then  !NNOTVECT
        do ICVL=ICVS_V,ICVE_V
        errsdy(0)=errsdy(0)+
     &          (tmp(ICVL,2)-tmp(ICVL,1))*(tmp(ICVL,2)-tmp(ICVL,1))
        erraby(0)=erraby(0)+
     &          (tmp(ICVL,1)*tmp(ICVL,1))
        enddo
!
        do ICOM=1,ncomp
        do ICVL=ICVS_V,ICVE_V
        errsdy(ICOM)=errsdy(ICOM)+
     &          (yys(ICVL,ICOM,2)-yys(ICVL,ICOM,1))
     &         *(yys(ICVL,ICOM,2)-yys(ICVL,ICOM,1))
        erraby(ICOM)=erraby(ICOM)+
     &          (yys(ICVL,ICOM,1)*yys(ICVL,ICOM,1)
     &          +yys(ICVL,ICOM,2)*yys(ICVL,ICOM,2))
        enddo
        enddo
!
        if(rns_scl) then
          do IMD=1,nrans
          do ICVL=ICVS_V,ICVE_V
          errsdk(IMD)=errsdk(IMD)+
     &            (aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
     &           *(aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
          errabk(IMD)=errabk(IMD)+
     &            (aks(ICVL,IMD,1))
     &           *(aks(ICVL,IMD,1))
          enddo
          enddo
        endif
! 
        do i=1,3
        do ICVL=ICVS_V,ICVE_V
        errsdv(i)=errsdv(i)+
     &          (vel(ICVL,i,2)-vel(ICVL,i,1))
     &         *(vel(ICVL,i,2)-vel(ICVL,i,1))
        errabv(i)=errabv(i)+
     &          (vel(ICVL,i,1))
     &         *(vel(ICVL,i,1))
        enddo
        enddo
! 
      else
        do 200 IIMAT=1,NMAT     !ICV=1,NCVIN
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        do ICVL=ICVS,ICVE
!
        errsdy(0)=errsdy(0)+
     &          (tmp(ICVL,2)-tmp(ICVL,1))*(tmp(ICVL,2)-tmp(ICVL,1))
        erraby(0)=erraby(0)+
     &          (tmp(ICVL,1)*tmp(ICVL,1))
!
        do 220 ICOM=1,ncomp
        errsdy(ICOM)=errsdy(ICOM)+
     &          (yys(ICVL,ICOM,2)-yys(ICVL,ICOM,1))
     &         *(yys(ICVL,ICOM,2)-yys(ICVL,ICOM,1))
        erraby(ICOM)=erraby(ICOM)+
     &          (yys(ICVL,ICOM,1)*yys(ICVL,ICOM,1)
     &          +yys(ICVL,ICOM,2)*yys(ICVL,ICOM,2))
 220    enddo
!
        if(rns_scl) then
          do 240 IMD=1,nrans
          errsdk(IMD)=errsdk(IMD)+
     &            (aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
     &           *(aks(ICVL,IMD,2)-aks(ICVL,IMD,1))
          errabk(IMD)=errabk(IMD)+
     &            (aks(ICVL,IMD,1))
     &           *(aks(ICVL,IMD,1))
 240      enddo
        endif
! 
        do 260 i=1,3
        errsdv(i)=errsdv(i)+
     &          (vel(ICVL,i,2)-vel(ICVL,i,1))
     &         *(vel(ICVL,i,2)-vel(ICVL,i,1))
        errabv(i)=errabv(i)+
     &          (vel(ICVL,i,1))
     &         *(vel(ICVL,i,1))
 260    enddo
!
!        if(ical_mhd>0) then
!          DO IMHD=1,2
!          do 270 i=1,3
!          errsdA(i,IMHD)=errsdA(i,IMHD)+
!     &            (MHD_A(ICVL,i,IMHD,2)-MHD_A(ICVL,i,IMHD,1))
!     &           *(MHD_A(ICVL,i,IMHD,2)-MHD_A(ICVL,i,IMHD,1))
!          errabA(i,IMHD)=errabA(i,IMHD)+
!     &            (MHD_A(ICVL,i,IMHD,1))
!     &           *(MHD_A(ICVL,i,IMHD,1))
! 270      enddo
!          enddo
!        endif
!
        enddo
! 

200     enddo

      endif
!
      if(NPE.gt.1) then
        do 225 ICOM=0,ncomp
        call hpcrsum(errsdy(ICOM))
        call hpcrsum(erraby(ICOM))
 225    continue
!
        do 245 IMD=1,nrans
        call hpcrsum(errsdk(IMD))
        call hpcrsum(errabk(IMD))
 245    continue
! 
        do 265 i=1,3
        call hpcrsum(errsdv(i))
        call hpcrsum(errabv(i))
 265    continue
!
        IF(ical_MHD>0) then
        DO IMHD=1,2
        do 275 i=1,3
        call hpcrsum(errsdA(i,IMHD))
        call hpcrsum(errabA(i,IMHD))
 275    continue
        enddo
        endif
      endif
!
! --- 
!
      errsdy(0)=dSQRT(errsdy(0))/(dSQRT(erraby(0))+SML)
      do 320 ICOM=1,ncomp
      errsdy(ICOM)=dSQRT(errsdy(ICOM))/(dSQRT(erraby(ICOM))+SML)
 320  continue
!
      do 340 IMD=1,nrans
      errsdk(IMD)=dSQRT(errsdk(IMD))/(dSQRT(errabk(IMD))+SML)
 340  continue
! 
      do 360 i=1,3
      errsdv(i)=dSQRT(errsdv(i))/(dSQRT(errabv(i))+SML)
 360  continue
!
      IF(ical_mhd>0) then
        DO IMHD=1,2
        do 370 i=1,3
        errsdA(i,IMHD)=dSQRT(errsdA(i,IMHD))/(dSQRT(errabA(i,IMHD))+SML)
 370    enddo
        enddo
      endif
!
! --- 
!
      errmax=0.d0
      do 100 icom=0,ncomp
      errmax=max(errmax,errsdy(icom))
 100  continue
      do 120 IMD=1,nrans
      errmax=max(errmax,errsdk(IMD))
 120  continue
      do 140 i=1,3
      errmax=max(errmax,errsdv(i))
 140  continue
      if(ical_mhd>0) then
        DO IMHD=1,2
        do 150 i=1,3
        errmax=max(errmax,errsdA(i,IMHD))
 150    enddo
        enddo
      endif
!
      if(NPE.gt.1) then
        call hpcrmax(errmax)
      endif
!
      if(errmax.lt.toltim) then
        stpcvg=.true.
      else
        stpcvg=.false.
      endif
!
      if(NPE.gt.1) then
        call hpcland(stpcvg,3)
      endif
!
      return
      end subroutine error_stp
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine initvar(MNVX,MAT_CV)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_hpcutil
      use module_dimension
      use module_hpc,only    : MNVX_GL
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(inout)  :: MNVX
      integer,intent(in)     :: MAT_CV(MXALLCV)
!
! --- [local entities]
!
      integer :: IV,ICVL
!
!
!
      if(MNVX_GL>NCVIN) THEN
        CALL FFRABORT
     &   (1,'ERR: Reset [MONITOR] in &hpc less then MAX. vertex number')
      ENDIF
      MNVX=0
      if(NPE.gt.1) then
        do 100 ICVL=1,NCVIN
        if(MNVX_GL==NOD_IDg(MAT_CV(ICVL))) then
          MNVX=ICVL
        endif
 100    continue
      else
        MNVX=MAT_CV(MNVX_GL)
      endif
!
      return
      end subroutine initvar
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine title
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only      : ifll
!
      implicit none
! --- [dummy arguments]
!
!
! --- [local entities]
!
      INTEGER :: val(8)
!
      write(ifll,*)
      write(ifll,1000)
      write(ifll,1050)
! --- 
      write(ifll,1100)
      write(ifll,1200)
! DEBUG
      write(ifll,1210)
      write(ifll,1220)
!
      CALL DATE_AND_TIME(VALUES=val)
      write(ifll,1300) val(1),val(2),val(3)
      write(ifll,1400) val(5),val(6),val(7)
! --- 
      write(ifll,1050)
      write(ifll,1000)
!
 1000 FORMAT(2X,'####',101('#'),'####')
 1050 FORMAT(2X,'    ',101(' '),'    ')
 1100 FORMAT(2X,'    ',37X,' WELCOME TO FrontFlow/red  ',37X,'    ')
 1200 FORMAT(2X,'    ',37X,' FrontFlow/red Version 2.0 ',37X,'    ')
 1210 FORMAT(2X,'    ',37X,' VERSION    = src_073.2')
 1220 FORMAT(2X,'    ',37X,' BUILD DATE = 2008-09-19-05:19:04')
 1300 FORMAT
     &(2X,'    ',38X,'Year:',I4.4,2X,'Mon :',
     &    I2.2,2X,'Date:',I2.2,36X,'    ')
 1400 FORMAT
     &(2X,'    ',38X,'Hour:',I4.4,2X,'Min :',
     &    I2.2,2X,'Sec :',I2.2,36X,'    ')
!
      return
      end subroutine title
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine mat_flag(iter,Pstart,mat_cal)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
!
      implicit none
! --- [DUMMY ARGUMENTS]
!
      INTEGER,intent(in) :: iter 
      INTEGER,intent(in) :: Pstart (MXMAT)
      logical,intent(out):: mat_cal(0:MXMAT)
!
! --- [local entities] 
!
      INTEGER :: IIMAT
!
      do 100 IIMAT=1,NMAT
      mat_cal(IIMAT)=.false.
      if(iter.ge.Pstart(IIMAT)) then
        mat_cal(IIMAT)=.true.
      elseif(iter.eq.0) then
        mat_cal(IIMAT)=.true.
      endif
 100  continue
!
      mat_cal(0)=.false.
!
      return
!
      end subroutine mat_flag
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
      subroutine prs_admin
     & (iphs,ipiso,
     &  iter,ismpl,deltt,time,nsmpl,ipfix,p_grdc,ictl,volctr,
     &  LVEDGE,LCYCSF,LBC_SSF,
     &  SFAREA,SFCENT,CVCENT,CVVOLM,CVVOL0,wiface,pp0,
     &  ccc,rva,vel,rho,prs,yys,tmp,velf,
     &  ccc2,rva2,vel2,rho2,prs2,yys2,tmp2,
     &  aks,SDOT,mass_i,
     &  EVAP_Y,
     &  grdc,diag,kdbp,dp,dr,bdyfrc,rvx,
     &  MAT_NO,MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_CFIDX,mat_cal,
     &  LCYCOLD,wifsld,OPPANG,FIELD_U,
     &  iterp,repsp,aepsp,errp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
!-------------------------------------------------
!     dp    <=   RMX(:)
!     dr    <=   cr(:)
!-------------------------------------------------
!
! --- [module arguments]
!
      use module_dimension
      use module_constant
      use module_hpcutil
!
      use module_io,only       : ifle,ifll
      use module_model,only    : idrdp,incomp,mach0,comp,ical_week,
     &                           ical_vect,nthrds,iLBF_P
      use module_material,only : nflud,lclsd,!,KMAT_S,MAT_S,
     &                           nofld,nsold,nosld,relaxp,ical_sld,
     &                           rotati,ishaft,end,begin,rot_ang,
     &                           dpmxv,
     &   ical_porous,porous,porosty,idarcy2,relaxv,relaxC
      use module_model,only    : ICVREF,PREFM,PEFC
      use module_boundary,only : kdbcnd,kdilet,kvlglw,kvnslp,kxnone,
     &                           kdintr,boundName,kdfire,kdbuff,kdcvd,
     &                           kdolet,kdprdc,kdsymm,kdtchi,kdshutr,
     &                           distrb,kdpres,kdstag,
     &                           LBC_INDEX,nbcnd,MAT_BCIDX,openout,
     &                           kdpres,kdstag,kdsld,idis,LBC_pair,
     &                           kdovst,lovstbc,nobcnd,
     &                           phs_idx,phs_com,surfreac,ical_poroBC
      use module_Euler2ph,only : ieul2ph,kdphs_g,kdphs_l,kdphs_s,kdphs
      USE module_usersub,ONLY  : src_r,usrno,usryes,src_fire
      use module_movegrid,only : ngrid
      use module_initial ,only : icvprs,rho0,rho02
      use module_VOF     ,only : ical_vof
      use module_chemreac,ONLY : ical_suf
      use module_scalar,  only : ivof,iaph,ical_cavi,ivold,ical_s
      use module_species, only : wm
      use module_flags  , only : SYMM
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_metrix,only   : rcomp,eva_comp,sinj_comp,matsld
      use module_metrix,only   : tmpfac=>d2vect
      use module_gravity ,only : beta,ggg
      use module_scalar,  only : ical_FC
      use module_time,    only : i_steady
      use module_metrix,only   : bdyf,vctr
      use module_model,   only : ical_prt,ical_dens
      use module_scalar  ,ONLY : calgeq,calxi
!      use module_metrix,only   : cnvq=>W1K9
!      use module_metrix,only   : cnvv=>W1K8
!
! 1.  Update velocity field
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,ismpl,nsmpl,iphs,ipiso
      real*8 ,intent(inout) :: deltt,time,p_grdc,volctr
      integer,intent(inout) :: ictl
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_NO(      0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CV(      MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(   0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)    
      integer,intent(inout) :: ipfix     (    MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)    :: CVCENT    (3,MXALLCV)
      real*8 ,intent(in)    :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)    :: CVVOL0    (  MXALLCV)
      real*8 ,intent(in)    :: wiface    (  MXCVFAC)
      real*8 ,intent(inout) :: ccc       (  MXALLCV)
      real*8 ,intent(inout) :: rva       (  MXCVFAC ,2)
      real*8 ,intent(inout) :: vel       (  MXALLCV ,3,2)
      real*8 ,intent(inout) :: velf      (  MXCVFAC_B ,3,2)
      real*8 ,intent(in)    :: yys       (  MXALLCV,MXcomp)
      real*8 ,intent(inout) :: rho       (  MXALLCV ,2)
      real*8 ,intent(inout) :: prs       (  MXALLCV ,2) 
      real*8 ,intent(inout) :: pp0       (  MXALLCV)
      real*8 ,intent(in)    :: tmp       (  MXALLCV,2)
!
      real*8 ,intent(in)    :: tmp2       (  MXALLCV2,2)
      real*8 ,intent(inout) :: ccc2       (  MXALLCV2)
      real*8 ,intent(inout) :: rva2       (  MXCVFAC2,2)
      real*8 ,intent(inout) :: vel2       (  MXALLCV2,3,2)
      real*8 ,intent(in)    :: yys2       (  MXALLCV2,MXcomp)
      real*8 ,intent(inout) :: rho2       (  MXALLCVC,2)
      real*8 ,intent(inout) :: prs2       (  MXALLCV2,2) 
!
      real*8 ,intent(in)    :: aks       (  MXALLCVR,mxrans,2)
      REAL*8 ,INTENT(IN)    :: SDOT      (  MXSSFBC_SUF,MXCOMPALL)
      integer,intent(in)    :: LCYCOLD   (  MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld    (  MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG    (  MXSSFBC_SLD)
!
      real*8 ,intent(inout) :: grdc      (  MXALLCV ,3,3)
      real*8 ,intent(inout) :: diag      (  MXALLCV)
      REAL*8 ,INTENT(INOUT) :: dp        (  MXALLCV)
      REAL*8 ,INTENT(INOUT) :: dr        (  MXALLCV)
      INTEGER,INTENT(IN)    :: kdbp      (  MXCVFAC)
      integer,intent(out)   :: ierror,iterp
      real*8 ,intent(out)   :: repsp,aepsp,errp
!     integer,intent(in)    :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(inout) :: bdyfrc(MXCV_B,3,2)
      real*8 ,intent(in)    :: mass_i    (   MXALLCV2)
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
      real*8 ,intent(in)    :: EVAP_Y(MXALLCV_P,MXCOMP)
      REAL*8 ,INTENT(INOUT) :: rvx       (  MXALLCV,3)
!
! --- [local entities]
!
      real*8  :: dpAB,dl,dll,dlvect,errabp,vol,ru,rv,rw,dxx,dyx,dzx
      real*8  :: dum1,dum2,dum3,dum4,dum5,dum6,dumcom,dumsim,
     &           rrdeltt,rdeltt,grdsf,dump5
      real*8  :: sinj_r,sinj_u,sinj_v,sinj_w,sinj_t
      real*8  :: grx,gry,grz
      real*8  :: u,v,w,dens,suu,suv,suw,utauu,t,p
      real*8  :: vof1A,vof2A,vof1B,vof2B,dumA,dumB,dumC,gf1
      integer :: i,j,k,l,m,n,ierr1=0,iclosed=1,idum
      integer :: ICOM,IMD,ICH,IFLD,ICTP,IDCA,ICVB1,IDCB,IDCBO
      integer :: IMAT,IIMAT,KMAT,IMAT_U,ICVS,ICVE,ICVL,IWL
      integer :: ICFS,ICFE,ICFL,ICFO,ICVBO
      integer :: IO,IN,IP,IS,IW,IV
      integer :: IC,IE,ICF,ICV,IBF,NB,KD,IDC,ICVFIX
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
      integer :: ICVLA,ICVLB,IBFS,IBFE,IMODE,IDCS,IDCE,IBFL
      integer :: ICODE
      integer :: imvflg=0
      logical :: outlet_3=.false.,lg_sld=.false.
      real*8  :: gf0,wi1,wi2,zerovof,dum,rvaerr,aph1,aph2
      real*8  :: eva_mass,unitmass,T_WALL,
     %           vel_av(3),dx2,dy2,dz2,ds1,ds2,ds3,
     &           eva_heat,eva_T,rvaicfp,rvaicf,ad_comp,rvasld
      character*80 :: BC_name
      integer :: IIMAT1,IIMAT2,iscal
      integer :: ICMAT,myid,ICMAT1,no
      integer :: IIMATS(2),IMATS(2),ISLD,ISLD2
      integer :: ipcom,icoms,icome,iph
      real*8  :: unit(3,2),rbb(3,3,2),fbb(3,3,2)
      real*8  :: sf11,sf12,sf13,ds21,ds22,ds23,dx,dy,dz,th(3)
      real*8  :: dxo,dyo,dzo,costh1,costh2,sinth1,sinth2,avsf
      real*8  :: org_x1,org_y1,org_x2,org_y2,dum_typ,dum_2
      real*8  :: flgmov
      real*8,save :: dum_eva=0.d0
      integer :: ICV_D,ICV_A
!
! --- rva_out array
!
      flgmov=1.d0
      dr=0.d0
      dp=0.d0
      rdeltt=1.d0/deltt
      rrdeltt=deltt
      diag=0.d0
      if(i_steady==3) then      ! prs_admin 
        rdeltt=1.d0 
        rrdeltt=1.d0  !1.d0
      endif
!
      dump5=0.d0
      IF(iLBF_P==5)then 
        dump5=1.d0
      endif
!----------------------------------------------------------
! --- source term in continuity equation: dp (Ax=b, dp=>b)
!----------------------------------------------------------
      if(iphs==2) goto 1000
      rvaerr=0.d0
      if(src_r.eq.usryes) then
         do 155 IIMAT=1,NMAT   !ICN=1,NCV
         if(.not.mat_cal(IIMAT)) goto 155
         IMAT=MAT_NO(IIMAT)
         ICVS=MAT_CVEXT(IIMAT-1)+1
         ICVE=MAT_CVEXT(IIMAT)
         if(IMAT.gt.0) then
           IMAT_U=nofld(IMAT)
           DO 156 ICVL=ICVS,ICVE
           ICV=MAT_CV(ICVL)
           if(NPE.gt.1) ICV=NOD_IDg(ICV)
           u=vel(ICVL,1,1)
           v=vel(ICVL,2,1)
           w=vel(ICVL,3,1)
           t=300.d0
           p=prs(ICVL,1) 
           do icom=1,ncomp
             rcomp(icom)=yys(ICVL,icom)
           enddo
           dens=rho(ICVL,1)
           sinj_u=0.d0
           sinj_v=0.d0
           sinj_w=0.d0
           sinj_r=0.d0
           sinj_t=0.d0
           sinj_comp(:)=0.d0
           vol=CVVOLM(ICVL)
           dxx=CVCENT(1,ICVL)
           dyx=CVCENT(2,ICVL)
           dzx=CVCENT(3,ICVL)
!
           call user_src_r(4,ictl,
     &     deltt,iter,time,ICV,IMAT_U,iphs,ncomp,dxx,dyx,dzx,
     &          u,v,w,T,p,dens,rcomp,vol,
     &          sinj_r,sinj_u,sinj_v,sinj_w,sinj_t,sinj_comp)
           if(abs(sinj_r).gt.SML) then   !for sinj_r>0 and <0
             dp(ICVL)=dp(ICVL)+sinj_r    !sinj_r=[kg/(s*m^3)]
             rvaerr=rvaerr+sinj_r!*vol
           endif
 156       continue
         elseif(IMAT.lt.0) then
           IMAT_U=nosld(-IMAT)
         endif
 155     continue
      endif
!
      if(ical_prt>=1) then 
         do IIMAT=1,NMAT
         if(.not.mat_cal(IIMAT)) cycle
         IMAT=MAT_NO(IIMAT)
         ICVS=MAT_CVEXT(IIMAT-1)+1
         ICVE=MAT_INDEX(IIMAT)
         if(IMAT.gt.0) then
           do ICVL=ICVS,ICVE
           ad_comp=0.d0
           do icom=1,ncomp
           ad_comp=ad_comp+EVAP_Y(ICVL,ICOM)  !
           enddo
           dum_eva=ad_comp !+ad_comp*CVVOLM(ICVL)
           dp(ICVL)=dp(ICVL)+dum_eva
           enddo
         endif
         enddo
      endif
!---------------------------------------
! --- Fire wall evaporation source term
!---------------------------------------
      if(src_fire.eq.usryes) then
        do nb=1,nbcnd
        BC_name=boundName(nb,1)
        IIMAT=MAT_BCIDX(nb,1)
        IMAT=MAT_NO(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        if(IMAT.lt.0) cycle
        IMAT_U=nofld(IMAT)
        kd=kdbcnd(0,nb)
        no=nobcnd(nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd.eq.kdfire) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL) 
          ICVL=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          u=vel(ICVL,1,1)
          v=vel(ICVL,2,1)
          w=vel(ICVL,3,1)
          T=tmp(ICVL,1)
          T_WALL=tmp(IDC,1)
          if(idrdp.eq.mach0.or.idrdp.eq.comp) then
            p=prs(ICVL,1)
          else
            p=prs(ICVL,1)+pp0(ICVL)
          endif
          dens=rho(ICVL,1)
          do icom=1,ncomp
            rcomp(icom)=yys(ICVL,icom)
          enddo
          eva_comp(:)=0.d0
          eva_mass=0.d0
          eva_heat=0.d0
          eva_T=tmp(ICVL,1)
!
          call user_src_fire
     &   (no,BC_name,deltt,iter,time,IMAT_U,iphs,ncomp,
     &    u,v,w,p,dens,rcomp,T_WALL,
     &    T,eva_T,eva_comp,eva_mass,eva_heat)
!
          unitmass=0.d0
          if(eva_mass.gt.0.d0) then     !unitmass=[kg/(s*m^3)]
            unitmass=eva_mass*SFAREA(4,ICFL)/CVVOLM(ICVL)
            dp(ICVL)=dp(ICVL)+unitmass  !dp(ICVL)=[kg/(s*m^3)]
          else
          endif
          enddo
        endif
        enddo
      endif
!
!------------------------------------------
! --- surface reaction (sticking function) 
!------------------------------------------
!
      if(ical_suf==1) then
      endif
!
!--------------------------------------------
! --- User source term for Euler two Phase
!--------------------------------------------
!
 1000 continue
      if(ieul2ph>0) then
      endif
!
!----------------------------------------
      if(iLBF_P==5) then
!        call dc_symprv
!     &    (1,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
!     &     LCYCOLD,wifsld,OPPANG,
!     &     SFAREA,SFCENT,bdyfrc(:,:,1),0,1)
        call grad_cell(3,7,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,bdyfrc(:,:,1),grdc)
      endif
!----------------------------------------
! --- dlt-density and Moving Mesh source 
!----------------------------------------
!
      if(ical_vect) then  !NOTHPC
        if(ieul2ph>0) then
          call FFRABORT(1,'ERR:Vector NOT support E2P')
        endif
!CIDR NODEP
!        if(idrdp==mach0) then
        do ICVL=ICVS_V,ICVE_V
        dp(ICVL)=(rho(ICVL,2)*CVVOL0(ICVL)
     &           -rho(ICVL,1)*CVVOLM(ICVL))*rdeltt
     &           -CVVOLM(ICVL)*dp(ICVL)
        enddo
!        endif
!
        do IE=1,MAXIE
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V 
        IF(ABS(vctr(ICVL,IE))>0) then
        dp(ICVL)=dp(ICVL)
     &    +sign(1,vctr(ICVL,IE))
     &    *rva(ABS(vctr(ICVL,IE)),1)
        ENDIF
        enddo
        ENDDO
      else
        if(ical_cavi==1.or.ical_FC==PEFC.or.ical_dens==4) then 
          if(ical_cavi==1) then  !cavit
            iscal=ivold
          elseif(ical_FC==PEFC) then
            iscal=ical_s
          endif
!
          if(ical_cavi==1) then 
          elseif(ical_FC==PEFC) then
            dp=0.d0
            do IIMAT=1,NMAT
            IMAT=MAT_NO(IIMAT)
            if(.not.mat_cal(IIMAT)) cycle   
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE 
            dp(ICVL)=
     &         -(CVVOL0(ICVL)*aks(ICVL,iscal,2)-
     &           CVVOLM(ICVL)*aks(ICVL,iscal,1))*rdeltt
     &          +CVVOLM(ICVL)*dp(ICVL)/rho(ICVL,1)
            enddo
            enddo

          endif
!
          if(ical_cavi==1.or.ical_FC==PEFC) then 
            do IIMAT=1,NMAT                 !ICF=1,NCVFAC 
            if(.not.mat_cal(IIMAT)) cycle 
            ICFS=MAT_CFIDX(IIMAT-1)+1 
            ICFE=MAT_CFIDX(IIMAT)
            dum1=rho0(IIMAT)
            dum2=rho02(IIMAT)
            do ICFL=ICFS,ICFE
            ICVLA=LVEDGE(1,ICFL)
            ICVLB=LVEDGE(2,ICFL)
            wi1=wiface(ICFL)
            wi2=1.d0-wiface(ICFL)
            dum=wi1*rho(ICVLA,1)+wi2*rho(ICVLB,1)
            dp(ICVLA)=dp(ICVLA)-rva(ICFL,1) !/dum  !cavit
            dp(ICVLB)=dp(ICVLB)+rva(ICFL,1) !/dum
            enddo
            enddo
          endif
!
          if(ical_dens==4) then 
          endif
!
        elseif(ical_vof==1) then
        elseif(ieul2ph>0) then
        else
          IF(iLBF_P==5) THEN 
            do IIMAT=1,NMAT 
            ICVS=MAT_CVEXT(IIMAT-1)+1 
            ICVE=MAT_CVEXT(IIMAT) 
            do ICVL=ICVS,ICVE 
            dp(ICVL)=dp(ICVL)-
     &      deltt*(grdc(ICVL,1,1)+grdc(ICVL,2,2)+grdc(ICVL,3,3))
!     &     deltt*(bdyfrc(ICVL,1,1)+bdyfrc(ICVL,2,1)+bdyfrc(ICVL,3,1))!/CVVOLM(ICVL)
            enddo
            enddo
          ENDIF
!
          do IIMAT=1,NMAT
            IMAT=MAT_NO(IIMAT)
            if(.not.mat_cal(IIMAT)) cycle 
            ICVS=MAT_CVEXT(IIMAT-1)+1 
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE 
            dp(ICVL)=
!zhang5 (1)
!     &          (CVVOL0(ICVL)*prs(ICVL,2)-
!     &           CVVOLM(ICVL)*prs(ICVL,1))*rdeltt/ccc(ICVL)
!zhang5
     &               (rho(ICVL,2)*CVVOL0(ICVL)
     &               -rho(ICVL,1)*CVVOLM(ICVL))*rdeltt*flgmov 
     &             +CVVOLM(ICVL)*dp(ICVL)
            enddo 
          enddo 
!zhang5  (2)
          if(idrdp==comp.and..false.) then
            dumA=0.d0
            dumB=1.d0-dumA
            rvx(:,1)=0.d0
            rvx(:,2)=0.d0
            do IIMAT=1,NMAT 
            IMAT=MAT_NO(IIMAT)
            if(.not.mat_cal(IIMAT).or.IMAT<0) cycle 
            ICFS=MAT_CFIDX(IIMAT-1)+1 
            ICFE=MAT_CFIDX(IIMAT)
            do ICFL=ICFS,ICFE 
            ICVA=LVEDGE(1,ICFL) 
            ICVB=LVEDGE(2,ICFL) 
            wi1=wiface(ICFL) 
            wi2=1.d0-wiface(ICFL) 
            ru=(wi1*vel(ICVA,1,1)+wi2*vel(ICVB,1,1)) 
            rv=(wi1*vel(ICVA,2,1)+wi2*vel(ICVB,2,1)) 
            rw=(wi1*vel(ICVA,3,1)+wi2*vel(ICVB,3,1)) 
            dum3=rho(ICVA,1)*wi1+rho(ICVB,1)*wi2
            dum1=SFAREA(4,ICFL)* 
     &      (SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw)
            dum2=max(0.d0,dum1)*rho(ICVA,1)
     &          +min(0.d0,dum1)*rho(ICVB,1)
            dum4=dum1*dum3
            rvx(ICVA,1)=rvx(ICVA,1)+(dum2*dumA+dum4*dumB)!+rva(ICFL,1) 
            rvx(ICVB,1)=rvx(ICVB,1)-(dum2*dumA+dum4*dumB)!-rva(ICFL,1) 
!            rvx(ICVA,2)=rvx(ICVA,2)+dum1
!            rvx(ICVB,2)=rvx(ICVB,2)-dum1
            enddo
            enddo
!
            do IIMAT=1,NMAT
            if(.not.mat_cal(IIMAT)) cycle
            IMAT=MAT_NO(IIMAT)
            if(IMAT.lt.0) cycle
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
!            dp(ICVL)=dp(ICVL)+(rho(ICVL,1)*rvx(ICVL,2)-rvx(ICVL,1))  !bad
            dp(ICVL)=dp(ICVL)-rvx(ICVL,1)        !(B 1)
!            dp(ICVL)=dp(ICVL)+(-rho(ICVL,1)*rvx(ICVL,2))  !bad
            enddo 
            enddo 
          else
            do IIMAT=1,NMAT
            if(.not.mat_cal(IIMAT)) cycle 
            ICFS=MAT_CFIDX(IIMAT-1)+1
            ICFE=MAT_CFIDX(IIMAT)
            do ICFL=ICFS,ICFE
            ICVLA=LVEDGE(1,ICFL)
            ICVLB=LVEDGE(2,ICFL)
            dp(ICVLA)=dp(ICVLA)-rva(ICFL,1) 
            dp(ICVLB)=dp(ICVLB)+rva(ICFL,1) 
            enddo
            enddo
          endif
!
!
!          if(lovstbc) then
!            do nb=1,nbcnd
!            kd=kdbcnd(0,nb)
!            if(kd==kdovst) then
!              IBFS=LBC_INDEX(nb-1)+1
!              IBFE=LBC_INDEX(nb)
!              do IBFL=IBFS,IBFE 
!              ICFL=LBC_SSF(IBFL) 
!              ICVP=LCYCSF(IBFL)
!              ICVLA=LVEDGE(1,ICFL)
!              ICVLB=LVEDGE(2,ICFL)
!              dum1=SFAREA(4,ICFL)*
!     &        (SFAREA(1,ICFL)*(vel(ICVLA,1,1))
!     &        +SFAREA(2,ICFL)*(vel(ICVLA,2,1))
!     &        +SFAREA(3,ICFL)*(vel(ICVLA,3,1)))
!              dum2=SFAREA(4,ICFL)*
!     &        (SFAREA(1,ICFL)*(vel(ICVLB,1,1))
!     &        +SFAREA(2,ICFL)*(vel(ICVLB,2,1))
!     &        +SFAREA(3,ICFL)*(vel(ICVLB,3,1)))
!              dp(ICVLA)=dp(ICVLA)-dum1-dum2
!              enddo
!            endif
!            enddo
!          endif
!
        endif
      endif
!
!--------------------------
! --- CV Variance of mass 
!--------------------------
!
      if(i_steady==3) then
        if(idrdp.eq.comp) then
          do IIMAT=1,NMAT
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          diag(ICVS:ICVE)=
     &    CVVOLM(ICVS:ICVE)/(ccc(ICVS:ICVE))/relaxp(IIMAT)
          dp(ICVS:ICVE)=dp(ICVS:ICVE)/relaxp(IIMAT)
          enddo
        else
          diag=0.d0
          call closed(ipfix)
          do IIMAT=1,NMAT
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          dp(ICVS:ICVE)=dp(ICVS:ICVE)!/relaxp(IIMAT)
          enddo
        endif
        if(ical_week) then
          do IIMAT=1,NMAT
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          diag(ICVS:ICVE)=
     &    CVVOLM(ICVS:ICVE)/(ccc(ICVS:ICVE))/relaxp(IIMAT)
          enddo
        endif
      else
        if(idrdp.eq.comp) then
          if(ical_cavi==1) then !cavit
            diag(:NCV)=CVVOLM(:NCV)/(ccc(:NCV))*rdeltt  !/rho(:NCV,1)
          elseif(ical_dens==4) then 
            diag(:NCV)=CVVOLM(:NCV)/(ccc(:NCV))*rdeltt
          else 
            diag(:NCV)=CVVOLM(:NCV)/(ccc(:NCV))*rdeltt
          endif
        else
          diag=0.d0
          call closed(ipfix)
        endif

        if(ical_week) then
          diag(:NCV)=rdeltt*CVVOLM(:NCV)/(ccc(:NCV))
        endif

        if(ieul2ph==2) then
          call bc_kdbp(iphs,LBC_SSF,LCYCSF,mat_cal,kdbp)
        endif 
      endif
!
!-----------------------------------
! --- < 6. Solve Poisson equation >-
!     ical_sld==1=>SYMM=2 OK
!     ical_sld==2=>SYMM=1 X
!-----------------------------------
!
      IMODE=5
      if(ical_sld==1.or.ical_sld==2.or.ical_sld==4) then
        SYMM=1
      else
        SYMM=1
      endif
!      if(idrdp.eq.comp) then
      if(idrdp.eq.comp.or.lovstbc) then
        SYMM=2          !=2
      endif
      if(ieul2ph==1) then
        SYMM=5
      elseif(ieul2ph==2) then
        SYMM=1
      endif
      if(ical_vect) then 
        SYMM=2          !=2
      endif
      if(ical_FC==PEFC) then
        SYMM=5
      endif
      if(ical_vof==1) then
        SYMM=5
      endif
!
!      if(ical_porous==1.or.ical_poroBC==1) then 
!        SYMM=2
!        IMODE=11
!      endif
!
      if(SYMM==1) then
        if(ical_vect) then
          call FFRABORT(1,'ERR: call solve_poisson_vect')
        else
          call solve_poisson
     &    (iter,nflud,rrdeltt,time,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,SFCENT,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,
     &    MAT_NO,mat_cal,MAT_CFIDX,
     &    ipfix,kdbp,diag,dp,rho(:,1),
     &    LCYCOLD,wifsld,OPPANG,
     &    iterp,repsp,aepsp,IMODE,ierr1)
        endif
      elseif(SYMM==2) then 
        if(ical_vect) then                   !NOVECT 
          call solve_poisson_unsymm_vect
     &    (iter,nflud,rrdeltt,time,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    ipfix,kdbp,diag,dp,rho(:,1),ccc,rva(:,1),
     &    LCYCOLD,wifsld,SFCENT,!vctr,
     &    iterp,repsp,aepsp,ierr1)
        else
          if(idrdp.eq.comp) then 
            IMODE=10 
          endif
          if(ical_dens==4) then 
            IMODE=10
          endif 
          call solve_poisson_unsymm
     &    (iter,nflud,rrdeltt,time,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,SFCENT,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    ipfix,kdbp,diag,dp,rho(:,1),rva(:,1),ccc,tmp,vel(:,:,1),
     &    LCYCOLD,wifsld,
     &    iterp,repsp,aepsp,IMODE,ierr1)
        endif
      elseif(SYMM==3) then
        call solve_poisson_1d
     &  (deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  ipfix,kdbp,diag,dp,rho(:,1),
     &  iterp,repsp,aepsp,ierr1)
      elseif(SYMM==5) then
         if(ieul2ph==1) then
           IMODE=5
         elseif(ical_FC==PEFC) then
           IMODE=8
         elseif(ical_vof==1) then
           IMODE=9 
         endif
         
!        call solve_poisson_e2p
!     &    (iter,nflud,deltt,time,
!     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
!     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
!     &    ipfix,kdbp,diag,dp,rho(:,1),rho2(:,1),aks(:,:,1),
!     &    LCYCOLD,wifsld,
!     &    iterp,repsp,aepsp,IMODE,ierr1)
         call solve_poisson_unsymm_e2p 
     &    (iter,nflud,deltt,time,
     &    LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &    MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &    ipfix,kdbp,diag,dp,rho(:,1),rho2(:,1),aks(:,:,1),rva(:,1),ccc,
     &    LCYCOLD,wifsld,
     &    iterp,repsp,aepsp,IMODE,ierr1)
      elseif(SYMM==4) then
        dp(:)=0.d0
        if(my_rank==root) then
          write(ifll,'(1X,A)') 'MSG: Pressure Poisson is Solved'
        endif
      endif
      
!
      if(ierr1.ne.0) goto 9999
!
!
!----------------------------------
! --- Boundary Condition on face: 
!----------------------------------
! --- kdbp(ICF) pressure BC
!              =0:periodic/CV-face
!              =1:outlet
!              =2:sysmetric/inlet/wall
!              =3:(inlet)
!--------------------------------
! --- Boundary Condition for dp 
!--------------------------------
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,dp)
!
!------------------------ 
! --- outlet BC for dp 
!------------------------
      do 500 nb=1,nbcnd 
      IIMAT=MAT_BCIDX(nb,1) 
      if(.not.mat_cal(IIMAT)) cycle 
      kd=kdbcnd(0,nb)
      if((kd==kdolet.and.openout(nb)/=5
     &  .and.openout(nb)/=8
!     &  .and.openout(nb)/=9
     &  .and.openout(nb)/=7)
     &  .or.   kd==kdpres.or. 
     &   (kd==kdstag.and.idrdp/=comp)) 
     &  then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        dp(IDC)=0.d0
        enddo
      elseif(kd==kdsld.and.idis(nb)>=1) then 
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICVB=LVEDGE(1,ICFP)
          ICVBO=LVEDGE(1,ICFO)
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
          dp(IDC)=wi1*dp(ICVB)+wi2*dp(ICVBO)
          enddo
      endif
  500 enddo
!
      if(iLBF_P==2.and..false.) then
        do nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        kd=kdbcnd(0,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd==kxnone)then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          dum2=(bdyfrc(ICVA,1,iphs))*(SFCENT(1,ICFL)-CVCENT(1,ICVA))
     &        +(bdyfrc(ICVA,2,iphs))*(SFCENT(2,ICFL)-CVCENT(2,ICVA))
     &        +(bdyfrc(ICVA,3,iphs))*(SFCENT(3,ICFL)-CVCENT(3,ICVA))
          dum1=vel(ICVA,1,1)*SFAREA(1,ICFL)
     &        +vel(ICVA,2,1)*SFAREA(2,ICFL)
     &        +vel(ICVA,3,1)*SFAREA(3,ICFL)
!          dum2=0.5d0*dum1*dum1*rho(ICVA,1)
          dp(ICVB)=dp(ICVA)+sign(dum2,-dum1)
          enddo
        endif
      enddo
      endif
!
! --- 
!
!-------------------------------
!--< 7. Modify velocity field >-
!------------------------------- 3333 
!
      if(iLBF_P==1) then              !important 
        call grad_cell_body(1,18,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,dp,bdyfrc(:,:,1),
     &  grdc(:,:,1))
      elseif(iLBF_P==0.or.iLBF_P==2.or.iLBF_P==4.or.iLBF_P==5) then
        call grad_cell(1,18,
     &       MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &       LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,dp,grdc(:,:,1))
      elseif(iLBF_P==3) then          !6666
        call grad_cell_least(1,18,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,rva(:,1),
     &  dp,grdc(:,:,1))
      endif
!
!----------------------------------------
! --- OUTLET, INLET and Touch-InletBC: >-
!----------------------------------------
!
      do 325 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 325
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if((kd.eq.kdolet
     &.or.kd.eq.kdilet.or.kd.eq.kdtchi).and.idrdp==comp) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &       grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &       grdc(ICV,3,1)*SFAREA(3,ICFL)
        grx=grdc(ICV,1,1)-dum1*SFAREA(1,ICFL)
        gry=grdc(ICV,2,1)-dum1*SFAREA(2,ICFL)
        grz=grdc(ICV,3,1)-dum1*SFAREA(3,ICFL)

        grdc(ICV,1,2)=dum1*SFAREA(1,ICFL)
        grdc(ICV,2,2)=dum1*SFAREA(2,ICFL)
        grdc(ICV,3,2)=dum1*SFAREA(3,ICFL)
        grdc(IDC,1,1)=dum1*SFAREA(1,ICFL)
        grdc(IDC,2,1)=dum1*SFAREA(2,ICFL)
        grdc(IDC,3,1)=dum1*SFAREA(3,ICFL)
        enddo

        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        grdc(ICV,1,1)=grdc(ICV,1,2)
        grdc(ICV,2,1)=grdc(ICV,2,2)
        grdc(ICV,3,1)=grdc(ICV,3,2)
        enddo
        
      elseif(kd.eq.kxnone.or.kd.eq.kdfire.or.
     &       kd.eq.kdintr.or.kd.eq.kdcvd) then
!        if(iLBF_P==2.or.iLBF_P==1.or.iLBF_P==3) cycle
        if(iLBF_P==2.or.iLBF_P==1) cycle
!        cycle
! --- wall 
        IF(ical_vect) then  !important-1
          do IBFL=IBFS,IBFE
          dum1=grdc(LVEDGE(1,LBC_SSF(IBFL)),1,1)
     &             *SFAREA(1,LBC_SSF(IBFL))+
     &         grdc(LVEDGE(1,LBC_SSF(IBFL)),2,1)
     &             *SFAREA(2,LBC_SSF(IBFL))+
     &         grdc(LVEDGE(1,LBC_SSF(IBFL)),3,1)
     &             *SFAREA(3,LBC_SSF(IBFL))
          tmpfac(IBFL,1)=dum1*SFAREA(1,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
          grdc(LVEDGE(1,LBC_SSF(IBFL)),1,2)=tmpfac(IBFL,1)
          enddo
!
          do IBFL=IBFS,IBFE
          dum1=grdc(LVEDGE(1,LBC_SSF(IBFL)),1,1)
     &             *SFAREA(1,LBC_SSF(IBFL))+
     &         grdc(LVEDGE(1,LBC_SSF(IBFL)),2,1)
     &             *SFAREA(2,LBC_SSF(IBFL))+
     &         grdc(LVEDGE(1,LBC_SSF(IBFL)),3,1)
     &             *SFAREA(3,LBC_SSF(IBFL))
          tmpfac(IBFL,1)=dum1*SFAREA(2,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
          grdc(LVEDGE(1,LBC_SSF(IBFL)),2,2)=tmpfac(IBFL,1)
          enddo
!
          do IBFL=IBFS,IBFE
          dum1=grdc(LVEDGE(1,LBC_SSF(IBFL)),1,1)
     &             *SFAREA(1,LBC_SSF(IBFL))+
     &         grdc(LVEDGE(1,LBC_SSF(IBFL)),2,1)
     &             *SFAREA(2,LBC_SSF(IBFL))+
     &         grdc(LVEDGE(1,LBC_SSF(IBFL)),3,1)
     &             *SFAREA(3,LBC_SSF(IBFL))
          tmpfac(IBFL,1)=dum1*SFAREA(3,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
          grdc(LVEDGE(1,LBC_SSF(IBFL)),3,2)=tmpfac(IBFL,1)
          enddo
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          grdc(ICV,1,1)=grdc(ICV,1,2)
          grdc(ICV,2,1)=grdc(ICV,2,2)
          grdc(ICV,3,1)=grdc(ICV,3,2)
          enddo

        else
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=(grdc(ICV,1,1))*SFAREA(1,ICFL)+
     &         (grdc(ICV,2,1))*SFAREA(2,ICFL)+
     &         (grdc(ICV,3,1))*SFAREA(3,ICFL)
          grdc(ICV,1,2)=dum1*SFAREA(1,ICFL)
          grdc(ICV,2,2)=dum1*SFAREA(2,ICFL)
          grdc(ICV,3,2)=dum1*SFAREA(3,ICFL)
          enddo
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          grdc(ICV,1,1)=grdc(ICV,1,2)
          grdc(ICV,2,1)=grdc(ICV,2,2)
          grdc(ICV,3,1)=grdc(ICV,3,2)
          enddo
        endif
      elseif(kd==kdprdc.and.idis(nb)==1) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        dum1=grdc(ICV,1,1)*SFAREA(1,ICFL)+
     &       grdc(ICV,2,1)*SFAREA(2,ICFL)+
     &       grdc(ICV,3,1)*SFAREA(3,ICFL)
        grdc(IDC,1,1)=dum1*SFAREA(1,ICFL)
        grdc(IDC,2,1)=dum1*SFAREA(2,ICFL)
        grdc(IDC,3,1)=dum1*SFAREA(3,ICFL)
        enddo
      elseif(kd==kdsld.and.idis(nb)>=1) then
      endif
 325  continue
!-------------------------------------------------------------
! --- < 8. Correcting pressure and New velocity: (New=n+1) >--
!-------------------------------------------------------------
      dumsim=1.d0
      if(iter.eq.1.and.ismpl.eq.1) then
        dumsim=0.d0
      endif
!
                        dumcom=0.d0
      if(idrdp.eq.comp) dumcom=1.d0
!      if(ical_dens==4) dumcom=0.d0
!
      if(ical_vect) then
         IIMAT=1
!CIDR NODEP
         DO ICVL=ICVS_V,ICVE_V
         dum2=rho(ICVL,1)+dr(ICVL)+dumcom*dp(ICVL)/(ccc(ICVL)+SML)
         prs(ICVL,1)=prs(ICVL,1)+dp(ICVL)*relaxp(IIMAT)
         rho(ICVL,1)=dum2
         if(rho(ICVL,1).lt.0) rho(ICVL,1) = rho(ICVL,2)    !!!okabe
         vel(ICVL,1,1)=(rho(ICVL,1)*vel(ICVL,1,1)
     &                      -rrdeltt*grdc(ICVL,1,1))/(dum2+SML)
         vel(ICVL,2,1)=(rho(ICVL,1)*vel(ICVL,2,1)
     &                      -rrdeltt*grdc(ICVL,2,1))/(dum2+SML)
         vel(ICVL,3,1)=(rho(ICVL,1)*vel(ICVL,3,1)
     &                      -rrdeltt*grdc(ICVL,3,1))/(dum2+SML)
         enddo
!
      else
        if(ieul2ph>0) then
        elseif(
     &         ical_FC==PEFC.or.
     &         ical_vof==1
     &         ) then
          dum_2=1.d0
          if(iLBF_P==2.or.iLBF_P==3) dum_2=0.d0 
          do IIMAT=1,NMAT
            IMAT=MAT_NO(IIMAT)
            if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle 
            dum1=porosty(IMAT) 
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
            prs(ICVL,1)=dum_2*prs(ICVL,1)+dp(ICVL)*relaxp(IIMAT)
            do IV=1,3
            vel(ICVL,IV,1)=vel(ICVL,IV,1)
     &          -rrdeltt*dum1*grdc(ICVL,IV,1)/rho(ICVL,1)
            enddo
            enddo
          enddo
        elseif(ical_cavi==1) then
          dum_2=1.d0
!          if(iLBF_P==2.or.iLBF_P==3) dum_2=0.d0 
          if(iLBF_P==2) dum_2=0.d0 
          do IIMAT=1,NMAT
            IMAT=MAT_NO(IIMAT)
            if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle 
            dum1=porosty(IMAT) 
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            do ICVL=ICVS,ICVE
            prs(ICVL,1)=dum_2*prs(ICVL,1)+dp(ICVL)*relaxp(IIMAT)
            dum2=rho(ICVL,1)+dr(ICVL)
     &          +dumcom*dp(ICVL)/(ccc(ICVL)+SML)*relaxp(IIMAT)
            do IV=1,3
            vel(ICVL,IV,1)=(rho(ICVL,1)*vel(ICVL,IV,1)
     &          -rrdeltt*grdc(ICVL,IV,1))/(dum2+SML)
            enddo
!!!!!!!!            rho(ICVL,1)=dum2 
            enddo
          enddo
        elseif(ical_dens==4) then
        else !ical_cavi==1
          dum_2=1.d0
!          if(iLBF_P==2.or.iLBF_P==3) dum_2=0.d0    !7777 
          if(iLBF_P==2) dum_2=0.d0    !7777
          do 600 IIMAT=1,NMAT 
          IMAT=MAT_NO(IIMAT) 
          if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle 
          ICVS=MAT_CVEXT(IIMAT-1)+1 
          ICVE=MAT_CVEXT(IIMAT) 
          dum1=porosty(IMAT)  
          if(iLBF_P==5) then
            do ICVL=ICVS,ICVE  
            prs(ICVL,1)=dum_2*prs(ICVL,1)+dp(ICVL)*relaxp(IIMAT)
            dum2=rho(ICVL,1)+dr(ICVL)+dumcom*dp(ICVL)/(ccc(ICVL)+SML)
            do IV=1,3
            vel(ICVL,IV,1)=(rho(ICVL,1)*vel(ICVL,IV,1)
     &             +dump5*rrdeltt*bdyfrc(ICVL,IV,1)
     &             -rrdeltt*grdc(ICVL,IV,1))/(dum2+SML)
            enddo
            rho(ICVL,1)=dum2
            enddo
          else
            
            do ICVL=ICVS,ICVE  
            prs(ICVL,1)=dum_2*prs(ICVL,1)+dp(ICVL)*relaxp(IIMAT)
            dum2=rho(ICVL,1)+dr(ICVL)
     &          +dumcom*dp(ICVL)/(ccc(ICVL)+SML)*relaxp(IIMAT)
!
            do 602 IV=1,3 !zhang5 (3) 
            vel(ICVL,IV,1)=vel(ICVL,IV,1)!*rho(ICVL,1)

!     &       -rrdeltt*grdc(ICVL,IV,1)/rho(ICVL,1)*relaxp(IIMAT) !(dum2+SML)
     &       -rrdeltt*grdc(ICVL,IV,1)/dum2!*relaxp(IIMAT) !(dum2+SML)
  602       enddo
            rho(ICVL,1)=dum2 
            enddo
          endif
  600     enddo
        endif
      endif 
!      if(ical_porous==1) then
!        do IIMAT=1,NMAT
!        IMAT=MAT_NO(IIMAT)
!        if(IMAT.le.0) cycle
!        IMAT_U=nofld(IMAT)
!        if(porous(IMAT)==idarcy2) then
!            ICVS=MAT_CVEXT(IIMAT-1)+1
!            ICVE=MAT_CVEXT(IIMAT)
!            DO ICVL=ICVS,ICVE
!            dp(ICVL)=dp(ICVL)*porosty(IMAT)
!            enddo
!        endif
!        enddo
!      endif

!----------------
! --- zhang-cvd
!----------------
      if(ical_vect) then
      else
      if(idrdp.eq.comp) then
        if(ical_dens==4) then
        else
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        DO ICV=ICVS,ICVE
        if(rho(ICV,1)<0.d0) then
          write(ifle,'(1x,a,E10.4,I8)') 'ERR: Negative density',
     &    rho(ICV,1),ICV
           write(ifle,'(1X,4(1x,/,a),I8)') 
     &     'MSG: Check 1) Initial pressure, too low pressure',
     &     '           2) Reaction Mechanism maybe in unbalance',
     &     '           3) Decreasing [unrelax_P] in fflow.ctl',
     &     '     ICV=',ICV
          call FFRABORT
     & (1,'ERR: Negative density; check initial-condition or reactions')
        endif
        enddo
        enddo
        endif
      endif
      endif
!
!----------------------------------------------
! --- New Mass flux on Face: m[n+1]=m[n]-dp*S 
!----------------------------------------------
!--------------
! --- CV-face:
!--------------
      if(ical_vect) then   !NNOTVECT
!CIDR NODEP
        IIMAT=1
        DO ICFL=ICFS_V,ICFE_V  !
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dum1=wi1*ccc(ICVA)+wi2*ccc(ICVB)
        dum4=wi1*dp(ICVA)+wi2*dp(ICVB)
        dum2=wi1*rho(ICVA,1)+wi2*rho(ICVB,1)
        dum3=rva(ICFL,1)/(dum1+SML)/dum2*dum4*relaxp(IIMAT)
        rva(ICFL,1)=rva(ICFL,1)-rrdeltt*SFAREA(4,ICFL)*!relaxp(IIMAT)*
     &   (dp(ICVB)-dp(ICVA))/
     &         abs((CVCENT(1,ICVB)-CVCENT(1,ICVA))*SFAREA(1,ICFL)
     &            +(CVCENT(2,ICVB)-CVCENT(2,ICVA))*SFAREA(2,ICFL)
     &            +(CVCENT(3,ICVB)-CVCENT(3,ICVA))*SFAREA(3,ICFL))
     &            +dum3*dumcom
!**2
!     &      *dsqrt((CVCENT(1,ICVB)-CVCENT(1,ICVA))**2
!     &            +(CVCENT(2,ICVB)-CVCENT(2,ICVA))**2
!     &            +(CVCENT(3,ICVB)-CVCENT(3,ICVA))**2)
        enddo
      elseif(ieul2ph>0) then
      elseif(ical_vof==1) then
      elseif(ical_cavi==1) then
      elseif(ical_dens==4) then 
      else
        if(iLBF_P==4) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          dum1=porosty(IMAT) 
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
          dpAB=dp(ICVB)-dp(ICVA)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &            +dy*SFAREA(2,ICFL)
     &            +dz*SFAREA(3,ICFL))
          dx=dx/dl
          dy=dy/dl
          dz=dz/dl
          dum1=wi1*rho(ICVA,1)+wi2*rho(ICVB,1)
          grdsf=porosty(IMAT)*rrdeltt*SFAREA(4,ICFL)*dpAB/dlvect  
          dum2=rrdeltt*dpAB/dlvect  

          rva(ICFL,1)=rva(ICFL,1)-grdsf     !+dum2  
          velf(ICFL,1,1)=velf(ICFL,1,1)-dum2/dum1*dx!SFAREA(1,ICFL)
          velf(ICFL,2,1)=velf(ICFL,2,1)-dum2/dum1*dy!SFAREA(2,ICFL)
          velf(ICFL,3,1)=velf(ICFL,3,1)-dum2/dum1*dz!SFAREA(3,ICFL)
          enddo
          enddo
        elseif(iLBF_P==5) then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
          dum1=porosty(IMAT) 
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)

          do ICFL=ICFS,ICFE
          ICVA=LVEDGE(1,ICFL)
          ICVB=LVEDGE(2,ICFL)
          wi1=wiface(ICFL)
          wi2=1.d0-wiface(ICFL)
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
          dpAB=dp(ICVB)-dp(ICVA)
          dlvect=abs(dx*SFAREA(1,ICFL)
     &              +dy*SFAREA(2,ICFL)
     &              +dz*SFAREA(3,ICFL))
          dum1=wi1*ccc(ICVA)+wi2*ccc(ICVB)
          dum4=wi1*dp(ICVA)+wi2*dp(ICVB)
          dum2=wi1*rho(ICVA,1)+wi2*rho(ICVB,1)
          dum3=rva(ICFL,1)/(dum1+SML)/dum2*dum4*relaxp(IIMAT)
          ru=(wi1*bdyfrc(ICVA,1,1)+wi2*bdyfrc(ICVB,1,1))
          rv=(wi1*bdyfrc(ICVA,2,1)+wi2*bdyfrc(ICVB,2,1))
          rw=(wi1*bdyfrc(ICVA,3,1)+wi2*bdyfrc(ICVB,3,1))
          dum2=SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw
          grdsf=porosty(IMAT)
     &         *rrdeltt*SFAREA(4,ICFL)*dpAB/dlvect         !*relaxv(IIMAT)
          rva(ICFL,1)=rva(ICFL,1)-grdsf+dum3*dumcom        !+dum2  
     &               +rrdeltt*dum2*dump5
          enddo
          enddo
        else  ! zhang5 (4) 
          do IIMAT=1,NMAT  
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle 
          dum1=porosty(IMAT)
          ICFS=MAT_CFIDX(IIMAT-1)+1 
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE 
          ICVA=LVEDGE(1,ICFL) 
          ICVB=LVEDGE(2,ICFL) 
          wi1=wiface(ICFL) 
          wi2=1.d0-wiface(ICFL) 
          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA) 
          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA) 
          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA) 
          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML) 
          ru=(wi1*vel(ICVA,1,1)+wi2*vel(ICVB,1,1))
          rv=(wi1*vel(ICVA,2,1)+wi2*vel(ICVB,2,1))
          rw=(wi1*vel(ICVA,3,1)+wi2*vel(ICVB,3,1))
          dum1=wi1*ccc(ICVA)+wi2*ccc(ICVB) 
          dum4=SFAREA(4,ICFL)*
     &    (SFAREA(1,ICFL)*ru+SFAREA(2,ICFL)*rv+SFAREA(3,ICFL)*rw)
          dum5= 
     &                dp(ICVA)*max(0.d0,dum4)
     &               +dp(ICVB)*min(0.d0,dum4)
          dum5=dum5/dum1
          dpAB=dp(ICVB)-dp(ICVA) 
          dum4=wi1*dp(ICVA)+wi2*dp(ICVB) 
          dum2=wi1*rho(ICVA,1)+wi2*rho(ICVB,1) 
          vof1A=rho(ICVA,1)+dr(ICVA) 
     &          +dumcom*dp(ICVA)/(ccc(ICVA)+SML)*relaxp(IIMAT)
          vof1B=rho(ICVB,1)+dr(ICVB) 
     &          +dumcom*dp(ICVB)/(ccc(ICVB)+SML)*relaxp(IIMAT)
          dum3=wi1*vof1A+wi2*vof1B 
          
          dlvect=abs(dx*SFAREA(1,ICFL) 
     &              +dy*SFAREA(2,ICFL) 
     &              +dz*SFAREA(3,ICFL))
          grdsf=rrdeltt*SFAREA(4,ICFL)*dpAB/dlvect*relaxp(IIMAT) 
!     &               *porosty(IMAT)*dum2/dum3
          dum5=rva(ICFL,1)/(dum1+SML)/dum2*dum4*relaxp(IIMAT)
          dum6=dum3/dum2
          dum5=0.d0
          rva(ICFL,1)=rva(ICFL,1)+dum5*dumcom-grdsf
          enddo
          enddo
!zhang5 (5)
          if(.false.)then
          do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(.not.mat_cal(IIMAT).or.IMAT.lt.0) cycle 
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          dum2=rho(ICVL,1)+dr(ICVL)
     &      +dumcom*dp(ICVL)/(ccc(ICVL)+SML)*relaxp(IIMAT)
          rho(ICVL,1)=dum2
          enddo
          enddo
          endif
!zhang5  (6)
          if(idrdp==comp.and..false.) then
            dp=0.d0
            do IIMAT=1,NMAT
            if(.not.mat_cal(IIMAT)) cycle 
            ICFS=MAT_CFIDX(IIMAT-1)+1
            ICFE=MAT_CFIDX(IIMAT)
            do ICFL=ICFS,ICFE
            ICVLA=LVEDGE(1,ICFL)
            ICVLB=LVEDGE(2,ICFL)
            dp(ICVLA)=dp(ICVLA)-rva(ICFL,1) 
            dp(ICVLB)=dp(ICVLB)+rva(ICFL,1) 
            enddo
            enddo
!
            do IIMAT=1,NMAT
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            DO ICVL=ICVS,ICVE
            dum1=rho(ICVL,1)
     &                +deltt*dp(ICVL)/CVVOLM(ICVL)
            do IV=1,3
            vel(ICVL,IV,1)=vel(ICVL,IV,1)
     &       -rrdeltt*grdc(ICVL,IV,1)/(dum1+SML)
!     &       -rrdeltt*grdc(ICVL,IV,1)/(rho(ICVL,1)+SML)
            enddo
            rho(ICVL,1)=dum1
            enddo
            enddo
          endif
!
!          do IIMAT=1,NMAT
!          IMAT=MAT_NO(IIMAT)
!          if(.not.mat_cal(IIMAT).or.IMAT<0) cycle
!          dum1=porosty(IMAT) 
!          ICFS=MAT_CFIDX(IIMAT-1)+1
!          ICFE=MAT_CFIDX(IIMAT)
!          do ICFL=ICFS,ICFE
!          ICVA=LVEDGE(1,ICFL)
!          ICVB=LVEDGE(2,ICFL)
!          wi1=wiface(ICFL)
!          wi2=1.d0-wiface(ICFL)
!          dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
!          dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
!          dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
!          dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)
!          dpAB=dp(ICVB)-dp(ICVA)
!          dlvect=abs(dx*SFAREA(1,ICFL)
!     &              +dy*SFAREA(2,ICFL)
!     &              +dz*SFAREA(3,ICFL))
!          dum1=wi1*ccc(ICVA)+wi2*ccc(ICVB)
!          dum4=wi1*dp(ICVA)+wi2*dp(ICVB)
!          dum2=wi1*rho(ICVA,1)+wi2*rho(ICVB,1)
!          dum3=rva(ICFL,1)/(dum1+SML)/dum2*dum4*relaxp(IIMAT)
!          grdsf=porosty(IMAT)
!     &         *rrdeltt*SFAREA(4,ICFL)*dpAB/dlvect         !*relaxv(IIMAT)
!          rva(ICFL,1)=rva(ICFL,1)-grdsf+dum3*dumcom        !+dum2  
!          enddo
!          enddo
        endif
      endif
!-------------------------
! --- rva in slinding BC  
!-------------------------
      if(ical_sld/=0.and..false.) then
!      if(ical_sld/=0) then
        do nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        kd=kdbcnd(0,nb)
        if(kd==kdovst) then
          IBFS=LBC_INDEX(nb-1)+1
          IBFE=LBC_INDEX(nb)
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVP=LCYCSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          
          enddo
        elseif(kd==kdsld.and..false.) then
          do  ISLD=1,2
          IF(ISLD==1) THEN
            ISLD2=2
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=IBFS+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2-1
            elseif(idis(nb)>=1) then
              IBFS=LBC_INDEX(nb-1)+1
              IBFE=LBC_pair(nb)
            endif
          ELSE
            ISLD2=1
            if(idis(nb)==0) then
              IBFS=LBC_INDEX(nb-1)+1+(LBC_INDEX(nb)-LBC_INDEX(nb-1))/2
              IBFE=LBC_INDEX(nb)
            elseif(idis(nb)>=1) then
              IBFS=LBC_pair(nb)+1
              IBFE=LBC_INDEX(nb)
            endif
          ENDIF
          IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
          if(.not.mat_cal(IIMATS(ISLD))) cycle
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          dum1=dsqrt(unit(1,ISLD)**2+unit(2,ISLD)**2+unit(3,ISLD)**2)
          unit(1,ISLD)=unit(1,ISLD)/dum1
          unit(2,ISLD)=unit(2,ISLD)/dum1
          unit(3,ISLD)=unit(3,ISLD)/dum1
          th(ISLD)=rot_ang(IMATS(ISLD))
          CALL rotth(unit(:,ISLD),th(ISLD),fbb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=fbb(j,i,ISLD)
          enddo
          enddo
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          dum1=dsqrt(unit(1,ISLD2)**2+unit(2,ISLD2)**2+
     &                unit(3,ISLD2)**2)
          unit(1,ISLD2)=unit(1,ISLD2)/dum1
          unit(2,ISLD2)=unit(2,ISLD2)/dum1
          unit(3,ISLD2)=unit(3,ISLD2)/dum1
          th(ISLD2)=rot_ang(IMATS(ISLD2))
          CALL rotth(unit(:,ISLD2),th(ISLD2),fbb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=fbb(j,i,ISLD2)
          enddo
          enddo
!
          costh1=cos(th(ISLD))
          sinth1=sin(th(ISLD))
          costh2=cos(th(ISLD2))
          sinth2=sin(th(ISLD2))
          org_x1=begin(1,IMATS(ISLD))*costh1
     &          -begin(2,IMATS(ISLD))*sinth1
          org_y1=begin(1,IMATS(ISLD))*sinth1
     &          +begin(2,IMATS(ISLD))*costh1
          org_x2=begin(1,IMATS(ISLD2))*costh2
     &          -begin(2,IMATS(ISLD2))*sinth2
          org_y2=begin(1,IMATS(ISLD2))*sinth2
     &          +begin(2,IMATS(ISLD2))*costh2
!
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
!
          ICVA=LVEDGE(1,ICFL)
          IDCA=LVEDGE(2,ICFL)
!
          ICVB=LVEDGE(1,ICFP)
          IDCB=LVEDGE(2,ICFP)

          ICVBO=LVEDGE(1,ICFO)
          IDCBO=LVEDGE(2,ICFO)
!
          sf11=SFAREA(1,ICFL)
          sf12=SFAREA(2,ICFL)
          sf13=SFAREA(3,ICFL)
          avsf=SFAREA(4,ICFL)
          dx=CVCENT(1,IDCA)-CVCENT(1,ICVA)
          dy=CVCENT(2,IDCA)-CVCENT(2,ICVA)
          dz=CVCENT(3,IDCA)-CVCENT(3,ICVA)
          dum1=dsqrt(dx*dx+dy*dy+dz*dz)
          dll=(dx*sf11+dy*sf12+dz*sf13)**2/dum1
          gf1=avsf*rrdeltt*(dp(IDCA)-dp(ICVA))/dll!*relaxp(IIMAT)
          rva(ICFL,1)=-gf1+rva(ICFL,1)
          enddo
          enddo
        endif
        enddo
      endif
!-----------------------------
! --- Mass Convergence error:
!-----------------------------
 7001 continue
      if(ical_vect) then  !NNOTVECT
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V
        diag(ICVL)=0.d0
        dp(ICVL)=0.d0
        enddo
!
        do IE=1,MAXIE  !vctr(ICVL,0)
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V  !index_c(myid)+1,index_c(myid+1)
        IF(ABS(vctr(ICVL,IE))>0) then
        diag(ICVL)=diag(ICVL)
     &       +abs(rva(ABS(vctr(ICVL,IE)),1))
        ENDIF
        enddo
        enddo
!
        errp=0.d0
        errabp=0.d0
!
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V !index_c(myid)+1,index_c(myid+1)
        errp=errp+(dp(ICVL)*dp(ICVL))
        errabp=errabp+diag(ICVL)*diag(ICVL)
        enddo
!
      else
        if(ieul2ph==1) then
        else
          diag=0.d0
          do 700 IIMAT=1,NMAT   !ICF=1,NCVFAC
          if(.not.mat_cal(IIMAT)) goto 700
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          diag(ICVLA)=diag(ICVLA)+ABS(rva(ICFL,1))
          diag(ICVLB)=diag(ICVLB)+ABS(rva(ICFL,1))
          enddo
  700     continue
!
          dp=0.d0
          do IIMAT=1,NMAT     !ICF=1,NCVFAC
          if(.not.mat_cal(IIMAT)) cycle
          ICFS=MAT_CFIDX(IIMAT-1)+1
          ICFE=MAT_CFIDX(IIMAT)
          do ICFL=ICFS,ICFE
          ICVLA=LVEDGE(1,ICFL)
          ICVLB=LVEDGE(2,ICFL)
          dp(ICVLA)=dp(ICVLA)-rva(ICFL,1)
          dp(ICVLB)=dp(ICVLB)+rva(ICFL,1)
          enddo
          enddo
!
          errp=0.d0
          errabp=0.d0
          Do 725 IIMAT=1,NMAT   !ICV=1,NCVIN
          if(.not.mat_cal(IIMAT)) goto 725
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_INDEX(IIMAT)
          do ICVL=ICVS,ICVE
          errp=errp+(dp(ICVL)*dp(ICVL))
          errabp=errabp+diag(ICVL)*diag(ICVL)
          enddo
  725     continue
        endif
      endif
!
      if(NPE.gt.1) then
        call hpcrsum(errp)
        call hpcrsum(errabp)
      endif
      errp=dsqrt(errp)/(dsqrt(errabp)+SML)
!
!
!
!----------------------------------------
! --- compressible flow
!----------------------------------------
      if(idrdp.eq.comp) then
        do 300 IIMAT=1,NMAT !ICV=1,NCV
        if(.not.mat_cal(IIMAT)) goto 300
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        pp0(ICVS:ICVE)=prs(ICVS:ICVE,1)
        pp0(IDCS:IDCE)=prs(IDCS:IDCE,1)
 300    enddo
      
      elseif(idrdp.eq.incomp.or.idrdp.eq.mach0) then
!      if(idrdp.eq.incomp.or.idrdp.eq.mach0.or.idrdp.eq.comp) then
!-------------------------------------------
! --- relative pressure (max volume)
!-------------------------------------------
        dpmxv(:)=0.d0
        if(imvflg.eq.0.or.(ngrid.gt.1)) then
          icvprs(:)=1
          ICMAT1=0
          if(my_rank.eq.root) then
            do IIMAT=1,NMAT
            if(lclsd(IIMAT)==2.or.
     &         lclsd(IIMAT)==0.or.lclsd(IIMAT)==1) then
              IMAT=MAT_NO(IIMAT)
              if(IMAT.gt.0) then
                dum1=-1.d0
                icvprs(IIMAT)=0
                dpmxv(IIMAT)=0.d0
                ICVS=MAT_CVEXT(IIMAT-1)+1
                ICVE=MAT_CVEXT(IIMAT)
                do ICVL=ICVS,ICVE
                if(CVVOLM(ICVL).gt.dum1) then
                  dum1=CVVOLM(ICVL)
                  icvprs(IIMAT)=ICVL
                  ICMAT1=ICVL
                endif
                enddo
              endif
            endif
            enddo
          endif
! --- sliding BC
!          matsld(:)=0
!          ICMAT=MXMAT+1
          ICMAT=icvprs(1)
          lg_sld=.false.
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
          if(kd==kdovst.or.kd==kdsld.or.
     &       kd==kdbuff.OR.kd==kdshutr) then
!              IIMAT1=MAT_BCIDX(nb,1)
!              IIMAT2=MAT_BCIDX(nb,2)
!              matsld(IIMAT1)=1
!              matsld(IIMAT2)=1
!              IF(ICMAT.GE.IIMAT1) THEN
!                 ICMAT=IIMAT1
!              ENDIF
!              IF(ICMAT.GE.IIMAT2) THEN
!                 ICMAT=IIMAT2
!              ENDIF
              lg_sld=.true.
          endif
          enddo

          IF(NPE>1) THEN
             call hpclor(lg_sld)
          ENDIF

!          IF(lg_sld) THEN
!            call hpcibcast(root,ICMAT)
!            icvprs(:)=ICMAT
!          ENDIF

!
          imvflg=1
!-------------------------------------
! --- 
!-------------------------------------
          do nb=1,nbcnd
          kd=kdbcnd(0,nb)
          IIMAT1=MAT_BCIDX(nb,1)
          IIMAT2=MAT_BCIDX(nb,2)
          if((kd==kdolet.and.openout(nb)==3).or.
     &       (kd==kdolet.and.openout(nb)==1).or.
!     &       (kd==kdolet.and.openout(nb)==4).or.
     &       (kd==kdpres).or.(kd==kdstag)) then 
            outlet_3=.true.
          endif
          enddo
          if(NPE>1) then
            call hpclor(outlet_3)
          endif
        endif

!-------------------------------------
! --- 
!-------------------------------------

        if(.not.outlet_3) then
          if(my_rank.eq.root) then
            do IIMAT=1,NMAT
            if(lclsd(IIMAT)==2.or.
     &         lclsd(IIMAT)==0.or.lclsd(IIMAT)==1) then
              dpmxv(IIMAT)=prs(icvprs(IIMAT),1)
            endif
            enddo
            if(lg_sld) then
              dum1=prs(icvprs(1),1)
              dpmxv(:)=dum1
            endif
          endif
!-----------------------------
! --- relative pressure 
!-----------------------------
          if(NPE.gt.1) then
!            DO KMAT=1,KMAT_S
!            IIMAT=MAT_S(KMAT)
!            if(lclsd(IIMAT)==2.or.
!     &         lclsd(IIMAT)==0.or.lclsd(IIMAT)==1) then
!              CALL hpcrbcast(root,dpmxv(IIMAT))
!            endif
!            enddo
            DO IIMAT=1,NMAT
            if(lclsd(IIMAT)==2.or.
     &         lclsd(IIMAT)==0.or.lclsd(IIMAT)==1) then
              CALL hpcrbcast(root,dpmxv(IIMAT))
            endif
            enddo
          endif
!
          if(ical_vect) then   ! NNOTVECT 
            DO ICVL=ICVS_V,ICVE_V 
            prs(ICVL,1)=prs(ICVL,1)-dpmxv(1)
            enddo
            DO ICVL=IDCS_V,IDCE_V
            prs(ICVL,1)=prs(ICVL,1)-dpmxv(1)
            ENDDO
          else
            
            do IIMAT=1,NMAT
            if(lclsd(IIMAT)==2.or.
     &         lclsd(IIMAT)==0.or.lclsd(IIMAT)==1) then
            if(.not.mat_cal(IIMAT)) cycle
            IMAT=MAT_NO(IIMAT)
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            IDCS=MAT_DCIDX(IIMAT-1)+1
            IDCE=MAT_DCIDX(IIMAT)
            if(IMAT.gt.0) then 
              do ICVL=ICVS,ICVE
              prs(ICVL,1)=prs(ICVL,1)-dpmxv(IIMAT)
              enddo
              do ICVL=IDCS,IDCE
              prs(ICVL,1)=prs(ICVL,1)-dpmxv(IIMAT)
              enddo
            endif
            endif
            enddo
          endif
        endif
!
        do IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          if(IMAT.lt.0) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_CVEXT(IIMAT)
            prs(ICVS:ICVE,1)=0.d0
          endif
        enddo
      endif
!
      if(ieul2ph==1) then
      endif
!
      if(iLBF_P==2.and..false.) then
        do nb=1,nbcnd
        IIMAT=MAT_BCIDX(nb,1)
        if(.not.mat_cal(IIMAT)) cycle
        kd=kdbcnd(0,nb)
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(kd.eq.kxnone .or. kd.eq.kdfire .or.
     &     kd.eq.kdintr .or.kd.eq.kdcvd) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum1=vel(ICV,1,1)*SFAREA(1,ICFL)
     &        +vel(ICV,2,1)*SFAREA(2,ICFL)
     &        +vel(ICV,3,1)*SFAREA(3,ICFL)
          vel(ICV,1,1)=vel(ICV,1,1)-dum1*SFAREA(1,ICFL)
          vel(ICV,2,1)=vel(ICV,2,1)-dum1*SFAREA(2,ICFL)
          vel(ICV,3,1)=vel(ICV,3,1)-dum1*SFAREA(3,ICFL)
          enddo
        endif
        enddo
      endif
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXALLCV,NCV,prs(:,1))
      ENDIF
!
      return
!
 9999 continue
!
      if(my_rank.eq.ROOT) write(ifle,*) '(prs_admin)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!
! --- < procedure for closed domain >
!
!=================================================
      subroutine closed(ipfix)
!=================================================
      use module_material ,only : iclosd
      use module_metrix,only   : dvdc,dpdt,rsdc
      implicit none
      integer,intent(out) :: ipfix(MXMAT)
      real*8  :: drdt,dt,dum1
      real*8  :: dvdcm,rsdcm
!
!-------------------------------------------------
!
      dvdc(:)=0.d0
      dpdt(:)=0.d0
      rsdc(:)=0.d0
      if(ical_vect) then   !NNOTVECT
        if(lclsd(1)==1) then
          dum1=-1.d0
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V
          dvdc(1)=dvdc(1)+CVVOLM(ICVL)/ccc(ICVL)
          rsdc(1)=rsdc(1)+dp(ICVL)
          if(CVVOLM(ICVL).gt.dum1) then
            dum1=CVVOLM(ICVL)
            ipfix(1)=ICVL
          endif
          enddo
        elseif(lclsd(1)==2) then
          IIMAT=1
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V
          dvdc(IIMAT)=dvdc(IIMAT)+CVVOLM(ICVL)/ccc(ICVL)
          rsdc(IIMAT)=rsdc(IIMAT)+dp(ICVL)
          enddo
!
          do nb=1,nbcnd
            IIMAT1=MAT_BCIDX(nb,1)
            IIMAT2=MAT_BCIDX(nb,2)
            kd=kdbcnd(0,nb)
            if(kd==kdsld.or.kd==kdbuff.OR.
     &         kd==kdshutr.or.kd==kdovst) then
              dvdcm=dvdc(IIMAT1)+dvdc(IIMAT2)
              rsdcm=rsdc(IIMAT1)+rsdc(IIMAT2)
              dvdc(IIMAT1)=dvdcm
              dvdc(IIMAT2)=dvdcm
              rsdc(IIMAT1)=rsdcm
              rsdc(IIMAT2)=rsdcm
              dvdc(0)=0.d0
              rsdc(0)=0.d0
            endif
          enddo
        endif
      else
        do 100 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          if(lclsd(IIMAT)==1) then
            dum1=-1.d0
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_INDEX(IIMAT)
            do 101 ICVL=ICVS,ICVE
            dvdc(IIMAT)=dvdc(IIMAT)+CVVOLM(ICVL)/ccc(ICVL)
            rsdc(IIMAT)=rsdc(IIMAT)+dp(ICVL)
            if(CVVOLM(ICVL).gt.dum1) then
              dum1=CVVOLM(ICVL)
              ipfix(IIMAT)=ICVL
            endif
 101        enddo
          elseif(lclsd(IIMAT)==2) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_INDEX(IIMAT)
            do ICVL=ICVS,ICVE
            dvdc(IIMAT)=dvdc(IIMAT)+CVVOLM(ICVL)/ccc(ICVL)
            rsdc(IIMAT)=rsdc(IIMAT)+dp(ICVL)
            enddo
            do nb=1,nbcnd
            IIMAT1=MAT_BCIDX(nb,1)
            IIMAT2=MAT_BCIDX(nb,2)
            kd=kdbcnd(0,nb)
! --- ---
            if(kd==kdsld.or.kd==kdbuff.OR.
     &        kd==kdshutr.or.kd==kdovst) then
              dvdcm=dvdc(IIMAT1)+dvdc(IIMAT2)
              rsdcm=rsdc(IIMAT1)+rsdc(IIMAT2)
              dvdc(IIMAT1)=dvdcm
              dvdc(IIMAT2)=dvdcm
              rsdc(IIMAT1)=rsdcm
              rsdc(IIMAT2)=rsdcm
              dvdc(0)=0.d0
              rsdc(0)=0.d0
            endif
            enddo
          endif
        endif
 100    enddo
      endif
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          if(lclsd(IIMAT)==1.or.lclsd(IIMAT)==2) then
!            call hpcrsum(dvdc(IIMAT))
!            call hpcrsum(rsdc(IIMAT))
!            dvdc(0)=0.d0
!            rsdc(0)=0.d0
!          endif
!        enddo
!
        DO IIMAT=1,NMAT
        if(lclsd(IIMAT)==1.or.lclsd(IIMAT)==2) then
          call hpcrsum(dvdc(IIMAT))
          call hpcrsum(rsdc(IIMAT))
!          dvdc(0)=0.d0
!          rsdc(0)=0.d0
        endif
        enddo
      endif
!
      do IIMAT=1,NMAT
      if(.not.mat_cal(IIMAT)) cycle
      if(lclsd(IIMAT)==1.or.lclsd(IIMAT)==2) then
        dpdt(IIMAT)=rsdc(IIMAT)/(dvdc(IIMAT)+SML)
      endif
      enddo
!
      if(ical_vect) then  !NNOTVECT
        IIMAT=1
        if(lclsd(1)==1) then
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V
          drdt=dpdt(IIMAT)/ccc(ICVL)
          dp(ICVL)=dp(ICVL)-CVVOLM(ICVL)*drdt
          enddo
        else
!CIDR NODEP
          DO ICVL=ICVS_V,ICVE_V
          dp(ICVL)=dp(ICVL)-CVVOLM(ICVL)*dpdt(IIMAT)/ccc(ICVL)
          enddo
        endif
        if(iclosd(IIMAT)==1.or.iclosd(IIMAT)==2) then
          if(idrdp.eq.mach0) then
!CIDR NODEP
            DO ICVL=ICVS_V,ICVE_V
            pp0(ICVL)=pp0(ICVL)+dpdt(IIMAT)*rrdeltt  !rrdeltt=deltt
            dr(ICVL)=dpdt(IIMAT)/ccc(ICVL)*rrdeltt
            enddo
          else
            dt=0.d0
          endif
        endif
      else
        do 200 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          if(lclsd(IIMAT)==1) then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_INDEX(IIMAT)
            do ICVL=ICVS,ICVE
            drdt=dpdt(IIMAT)/ccc(ICVL)
            dp(ICVL)=dp(ICVL)-CVVOLM(ICVL)*drdt
            enddo
          elseif(lclsd(IIMAT)==2)then
            ICVS=MAT_CVEXT(IIMAT-1)+1
            ICVE=MAT_INDEX(IIMAT)
            do ICVL=ICVS,ICVE
            drdt=dpdt(IIMAT)/ccc(ICVL)
            dp(ICVL)=dp(ICVL)-CVVOLM(ICVL)*drdt
            enddo
          endif
        endif
 200    enddo
!
        do 300 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        if(IMAT.gt.0) then
          if(iclosd(IIMAT)==1.or.iclosd(IIMAT)==2) then
            if(idrdp.eq.mach0) then
              dt=deltt
              ICVS=MAT_CVEXT(IIMAT-1)+1
              ICVE=MAT_INDEX(IIMAT)
              do 103 ICVL=ICVS,ICVE
              drdt=dpdt(IIMAT)/ccc(ICVL)
              pp0(ICVL)=pp0(ICVL)+dpdt(IIMAT)*dt
              dr(ICVL)=drdt*dt
 103          continue
            else
              dt=0.d0
            endif
          endif
        endif
 300    enddo
      endif
!
      end subroutine closed
!
      end subroutine prs_admin
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine suf_admin(
     &  iter,deltt,
     &  MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,LBC_SSF,LVEDGE,
     &  SFAREA,SFCENT,CVCENT,CVVOLM,
     &  vel,prs,rho,tmp,yys,wdot,FIELD_U,
     &  SDOT,MOLFRC,MOLCEN,SUFBND,SITDEN,BLKTHK,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     dZ/dt=sigma*sdot/rho-Z/rho*d(rho)/dt-Z/sigma*d(sigma)/dt
!     1) steady: 
!
!
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_constant
      use module_io,only       : ifll,ifle
      use module_boundary,only :nbcnd,kdbcnd,MAT_BCIDX,LBC_INDEX,kdcvd,
     &                          phs_idx,phs_dns,
     &                          phs_nbc,phs_typ,phs_snm,phs_nam,
     &                          gasphase,surphase,blkphase,phs_com,
     &                          phs_inifrac,num_site,blk_dns,phs_thk,
     &                          surfreac
      use module_species,only  : wm,r_wm
!     
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter
      real*8 ,intent(in)    :: deltt
      INTEGER,INTENT(IN)    :: MAT_CVEXT( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX( 0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO   ( 0:MXMAT)
      logical,INTENT(IN)    :: mat_cal  ( 0:MXMAT)
      INTEGER,INTENT(IN)    :: LBC_SSF(   MXSSFBC)
      INTEGER,INTENT(IN)    :: LVEDGE (2, MXCVFAC)
      REAL*8 ,INTENT(IN)    :: CVCENT (3, MXALLCV)
      REAL*8 ,INTENT(IN)    :: CVVOLM(    MXALLCV)
      REAL*8 ,INTENT(IN)    :: SFAREA (4, MXCVFAC)
      REAL*8 ,INTENT(IN)    :: SFCENT (3, MXCVFAC) 
      real*8 ,intent(in)    :: rho     (  MXALLCV) 
      real*8 ,intent(in)    :: tmp     (  MXALLCV) 
      real*8 ,intent(in)    :: vel     (  MXALLCV,3) 
      real*8 ,intent(in)    :: prs     (  MXALLCV) 
      real*8 ,intent(in)    :: yys     (  MXALLCV,MXcomp) 
      real*8 ,intent(in)    :: wdot    (  MXALLCV,MXcomp) 
      REAL*8 ,INTENT(IN)    :: SDOT    (  MXSSFBC_SUF,MXCOMPALL) 
      REAL*8 ,INTENT(INOUT) :: MOLFRC  (  MXSSFBC_SUF,MXCOMPALL,2) 
      REAL*8 ,INTENT(IN)    :: MOLCEN  (  MXSSFBC_SUF,MXCOMPALL) 
      REAL*8 ,INTENT(INOUT) :: SUFBND  (  MXSSFBC_SUF,MXCOMPALL) 
      REAL*8 ,INTENT(INOUT) :: SITDEN  (  MXSSFBC_SUF,MXPHASE,2) 
      REAL*8 ,INTENT(INOUT) :: BLKTHK  (  MXSSFBC_SUF,MXPHASE,2) 
      real*8 ,intent(inout) :: FIELD_U(MXCV_D,NFLID)
      integer,intent(out)   :: ierror
!
! --- [local entities]
!
      integer :: nb,kd,nbx,ntp,icom,icoms,icome,iph,ipcom
      integer :: IBFS,IBFE,IBFL,ICFL,ICVL
      integer :: IIMAT,IMAT
      real*8  :: dum1,dum2,dum3,dum4,thichness
!     
      ierror=0
!
      return
      end subroutine suf_admin
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
