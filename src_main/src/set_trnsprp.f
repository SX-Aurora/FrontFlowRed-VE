!
!     subroutine set_trnsprp
!     subroutine set_vof_t2cph
!     subroutine set_trnsprp_t2h_cavi
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine set_trnsprp(mph,
     &  LVEDGE,LBC_SSF,LCYCSF,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,wifsld,LCYCOLD,
     &  hhh,hhs,cps,cp,cr,
     &  tmp,yys,aks,prs,rho,
     &  rmu,rmd,rds,rmx)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments] 
!
      use module_dimension 
      use module_material ,only : isthr,suther,const,simpli,MK,usrprp,
     &                            sthrmu,sthrt,sthrc,
     &                            prlmnr,sclmnr,nflud,
     &                            sthrmu2,sthrt2,sthrc2,
     &                            prlmnr2,sclmnr2,porosty,
     &                            chung
      use module_species,  only : wm,chkncomp,sw,diffu,gascns
      use module_material, only : rmdsld,nofld
      use module_Euler2ph,only  : ieul2ph
      use module_model,only     : ical_vect,nthrds,ical_thmoDIFF
      use module_vector,only    : ICVS_V,ICVE_V,
     &                            ICFS_V,ICFE_V,
     &                            ICVSIN_V,ICVEIN_V,
     &                            IDCS_V,IDCE_V,index_c,index_f
      use module_scalar,only    : ical_cavi,Tref,ivold,
     &                            icavi,T0_ref,ical_FC
      use module_FUEL    ,only : No_Mem,No_AGDL,No_CGDL,OPEN_V,Tau_diff
!      use module_metrix,only   : zz_temp=>d2work1!,zs_temp=>d2work2
      use module_species  ,only : r_wm
!      use module_metrix,only   : rcomp,TH_DIFF,r_wmm=>ys,wmm=>eva_comp,
!     &                           Xi=>sinj_comp
      use module_metrix,only   : W1K9,W1K8
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mph
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)

      REAL*8 ,INTENT(INOUT) :: hhh   (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: hhs   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cps   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp    (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cr    (       MXALLCV)
      real*8 ,intent(in)    :: tmp   (       MXALLCV)
      real*8 ,intent(in)    :: prs   (       MXALLCV)
      real*8 ,intent(in)    :: yys   (       MXALLCV,MXcomp)
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(out)   :: rmu   (       MXALLCV)
      real*8 ,intent(out)   :: rmd   (       MXALLCV)
      real*8 ,intent(out)   :: rds   (       MXALLCV,MXcomp)
      real*8 ,intent(out)   :: rmx   (       MXALLCV)
      real*8 ,intent(in)    :: aks   (       MXALLCVR,mxrans)
      real*8 ,intent(in)    :: rho   (       MXALLCV)
!
! --- [local entities]
!
      real*8  :: vwm    (   mxcomp)
      real*8  :: zs     (   mxcomp),SML=1.D-15
      real*8  :: tt,zz,eps,dum1,dum2
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,ICOM,myid,IMAT_U
      real*8,allocatable :: 
     &    sthrmu_t(:),sthrt_t(:),sthrc_t(:),prlmnr_t(:),sclmnr_t(:)
      real*8,parameter :: A=2.58d-5
             ! constant for simplified transport model [kg/(m*s)]
      real*8,parameter :: T0=298.d0
             ! reference temperature for simplified transport model [K]

!
      call chkncomp(ncomp,'(set_trnsprp)')
!
!
!
      ALLOCATE(sthrmu_t(nflud),
     &         sthrt_t(nflud),
     &         sthrc_t(nflud),
     &         prlmnr_t(nflud),
     &         sclmnr_t(nflud))
!
! 
!
      if(ieul2ph>0.and.mph==2) then
        sthrmu_t(:)=sthrmu2(:)
        sthrt_t(:)=sthrt2(:)
        sthrc_t(:)=sthrc2(:)
        prlmnr_t(:)=prlmnr2(:)
        sclmnr_t(:)=sclmnr2(:)
      else
        sthrmu_t(:)=sthrmu(:)
        sthrt_t(:)=sthrt(:)
        sthrc_t(:)=sthrc(:)
        prlmnr_t(:)=prlmnr(:)
        sclmnr_t(:)=sclmnr(:)
      endif
!
!-< 1. Fluid part >-
!
!
!--< 1.1 specific heat at constant pressure >--
!
      if(sw) then
        call cal_t2cp(MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,prs,rho,rmu,rmd,cps,cp,cr)
      else
        call cal_t2hcp(mph,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,
     &           tmp,yys,hhs,cps,hhh,cp,cr,prs)
      endif 
! 
      if(ical_vect.and..false.) then 
!      if(ical_vect) then 
        IMAT=1 
!CIDR NODEP 
        
        DO ICVL=ICVS_V,ICVE_V !index_c(myid)+1,index_c(myid+1)
        rmu(ICVL)=sthrmu_t(IMAT)
        rmd(ICVL)=cp(ICVL)*rmu(ICVL)/prlmnr_t(IMAT)
        enddo
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V !index_c(myid)+1,index_c(myid+1)
        rmx(ICVL)=rmu(ICVL)
        enddo
!
        cr(:)=0.d0
        do icom=1,ncomp
!CIDR NODEP
        do ICVL=ICVS_V,ICVE_V
        dum1=yys(ICVL,icom)*r_wm(icom)
        rds(ICVL,icom)=dum1
        cr(ICVL)=cr(ICVL)+dum1
        enddo
        enddo
!
        do icom=1,ncomp
        do ICVL=ICVS_V,ICVE_V
        rds(ICVL,icom)=rmu(ICVL)*cr(ICVL)/sclmnr_t(IMAT)
     &     *(1.d0-yys(ICVL,icom))/max(1.d-20,cr(ICVL)-rds(ICVL,icom))
        enddo
        enddo
!
        call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmx)
        call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmu)
        DEALLOCATE(sthrmu_t,sthrt_t,sthrc_t,prlmnr_t,sclmnr_t)
        return
      endif
!------------------------------------------
!--< 1.2 viscosity >-- Scalar computation
!------------------------------------------
      do IIMAT=1,NMAT 
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        if(IMAT>0) then       ! fluid
          IMAT_U=nofld(IMAT)
          select case(isthr(IMAT))
            case(suther)                ! Sutherland's formula
              call model_Sutherland(ICVS,ICVE)
            case(const)                 ! constant viscosity
              if(ical_vect) then
                call model_Constant_V(ICVS,ICVE)
              else
                call model_Constant(ICVS,ICVE)
              endif
            case(simpli)  ! ical_vect ! simplified transport model
              call model_Simpli(ICVS,ICVE,prlmnr_t,
     &            tmp,cp,rmu,rmd,rds,
     &            IMAT,MXALLCV,nflud,MXcomp)
            case(MK)                    ! Mix_Kinetic model 
              call model_Mix_Kinetic(ICVS,ICVE)  
            case(usrprp)
              call USER_PRO(IMAT_U,ICVS,ICVE)
            case(chung)

          end select
        else                    ! solid
          call cal_solid(ICVS,ICVE)
        endif
      enddo
!
!
!
      if(ical_FC>0) then
        do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        if(IMAT<0)cycle
        IMAT_U=nofld(IMAT)
        if(IMAT_U==No_Mem.or.IMAT_U==No_AGDL.or.IMAT_U==No_CGDL) then 
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          do ICVL=ICVS,ICVE
          dum1=tmp(ICVL)
          dum2=prs(ICVL)
          dum1=(dum1/300.d0)**1.5d0*(1.01325d5/dum2)
          do icom=1,ncomp
          rds(ICVL,icom)=dum1*diffu(icom)
          enddo
          enddo
        endif
        enddo
      endif
!
!-< 3. Scale of diffusivity >-
!
      do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT>0) then       ! fluid
          do ICVL=ICVS,ICVE
            rmx(ICVL)=max(rmu(ICVL),rmd(ICVL)/(cp(ICVL)+SML))
            do icom=1,ncomp
            rmx(ICVL)=max(rmx(ICVL),rds(ICVL,icom))
            enddo
          enddo
        else                    ! solid
          rmx(ICVS:ICVE)=rmd(ICVS:ICVE)/cp(ICVS:ICVE)
        endif
      enddo
!
!-< 4. Set value at dummy cells >-
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmu)
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmd)
!
      do ICOM=1,ncomp
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rds(:,ICOM))
      enddo
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmx)
!
      DEALLOCATE(sthrmu_t,sthrt_t,sthrc_t,prlmnr_t,sclmnr_t)
      return
!
!///////////////////////////////////////////////////////////////////////
      contains
!==============================================================================
      subroutine model_Sutherland(do_start,do_end)      ! Sutherland's formula
!==============================================================================
      use module_material ,only : sthrmu,sthrt,sthrc,r_prlmnr
      integer,intent(in) :: do_start,do_end
!
      dum1 = sthrmu_t(IMAT)*(sthrt_t(IMAT)+sthrc_t(IMAT))
      do ICVL=do_start,do_end
        rmu(ICVL) = dum1/(tmp(ICVL)                 ! viscosity [Pa s]
     &     +sthrc_t(IMAT))*sqrt((tmp(ICVL)/sthrt_t(IMAT))**3)
        rmd(ICVL)=cp(ICVL)*rmu(ICVL)/prlmnr_t(IMAT) ! thermal conductivity [J/(m s K)]
        call diffusivity_density(ICVL)              ! (mass diffusivity of 
                                                    ! species)*(density) [kg/(m s)]
      enddo
      end subroutine model_Sutherland
!
!===========================================================================
      subroutine model_Constant(do_start,do_end)    ! constant viscosity
!===========================================================================
      use module_material ,only : sthrmu,r_prlmnr 
      integer,intent(in) :: do_start,do_end
      do ICVL=do_start,do_end
! --- constant
        rmu(ICVL) = sthrmu_t(IMAT)
      ! viscosity [Pa s]
        rmd(ICVL) = cp(ICVL)*rmu(ICVL)/prlmnr_t(IMAT)
      ! thermal conductivity [J/(m s K)]
!
!--------------------------------------------------------------------
! --- FLUENT
!
!        rmu(ICVL)=3.17d-6+2.04D-8*tmp(ICVL)-3.52D-12*tmp(ICVL)**2
!        rmd(ICVL)=3.8D-2+5.41D-4*tmp(ICVL)
!     &           -2.51D-7*tmp(ICVL)**2+8.57D-11*tmp(ICVL)**3
!--------------------------------------------------------------------
        call diffusivity_density(ICVL)          ! (mass diffusivity of 
                                                ! species)*(density) [kg/(m s)]
      enddo
!
      end subroutine model_Constant
!
!===========================================================================
      subroutine model_Constant_V(do_start,do_end)    ! constant viscosity
!===========================================================================
      use module_material ,only : sthrmu,r_prlmnr 
      integer,intent(in) :: do_start,do_end
!
      real*8,parameter   :: eps=1.d-20
      real*8             :: zs(mxcomp),zz
!
! -------------------------------------------------------------
!
      do ICVL=do_start,do_end
! --- constant
        rmu(ICVL) = sthrmu_t(IMAT)
      ! viscosity [Pa s]
        rmd(ICVL) = cp(ICVL)*rmu(ICVL)/prlmnr_t(IMAT)
      ! thermal conductivity [J/(m s K)]
!
      enddo
!------------------------------
! --- diffusivity_density
! (mass diffusivity of species)*(density) [kg/(m s)]
!------------------------------
      W1K8(:)=0.d0
      do icom=1,ncomp
      do ICVL=do_start,do_end
      W1K8(ICVL)=W1K8(ICVL)+yys(ICVL,icom)*r_wm(icom)
      enddo
      enddo
!
      do icom=1,ncomp
      do ICVL=do_start,do_end
      dum1=yys(ICVL,icom)*r_wm(icom)
      rds(ICVL,icom)=rmu(ICVL)*W1K8(ICVL)/sclmnr_t(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,W1K8(ICVL)-dum1)
      enddo
      enddo
!
      end subroutine model_Constant_V
!


!
!===========================================================================
      subroutine USER_PRO(IMAT_U,do_start,do_end)        ! constant viscosity
!===========================================================================
      use module_io,only        : ifle,ifll
      use module_material ,only : sthrmu,r_prlmnr 
      use module_metrix,only    : mu_i,k_i,thmdif,difu
      use module_metrix,only   : cp_i=>rcomp,ys,TH_DIFF,wmm=>eva_comp,
     &                           Xi=>sinj_comp
!
      integer,intent(in) :: IMAT_U,do_start,do_end
!
      real*8  :: mu_usr,rmd_usr,difu_usr
      real*8  :: tu,pu,rhou,cpui,cpu,T_MW
      integer :: i_mu,i_rmd,i_diffu,i_thm_dif
!
!
      do ICVL=do_start,do_end
      wmm(:)=wm(:)
      ys(:)=yys(ICVL,:)
      dum1=0.d0
      T_MW=0.d0
      do icom=1,ncomp
      dum1=dum1+ys(ICOM)*r_wm(ICOM)
      enddo
      dum1=1.d0/dum1
      do icom=1,ncomp
      Xi(ICOM)=(ys(ICOM)*r_wm(ICOM))*dum1
      enddo
      do icom=1,ncomp
      T_MW=T_MW+Xi(icom)/wmm(icom)
      enddo
      T_MW=1.d0/T_MW
!
! --- constant 
!
      tu=tmp(ICVL)
      pu=prs(ICVL)
      rhou=rho(ICVL)
      cp_i(:)=cps(ICVL,:)
      cpu=cp(ICVL)
!
      mu_usr=-1.d0
      rmd_usr=-1.d0
      difu(:)=-1.d0
      thmdif(:)=-1.d0
!
      i_mu=0
      i_rmd=0
      i_diffu=0
      i_thm_dif=0
!
      call USER_PROPERTY(ICVL,IMAT_U,MXcomp,wmm,
!     &           i_mu,i_rmd,i_diffu,i_thm_dif,
     &           T_MW,tu,pu,rhou,Ys,Xi,cp_i,cpu,
     &           mu_usr,rmd_usr,difu,thmdif)
      if(mu_usr<0.d0) then 
        write(ifle,*) 'ERR: [mu_usr] is not defined in USER_PROPERTY'
        call FFRABORT(1,'ERR: in USER_PROPERTY')
      else
        rmu(ICVL)=mu_usr
      endif
!
      if(rmd_usr<0.d0) then 
        write(ifle,*) 'ERR: [rmd_usr] is not defined in USER_PROPERTY'
        call FFRABORT(1,'ERR: in USER_PROPERTY')
      else
        rmd(ICVL)=rmd_usr 
      endif
!
      do ICOM=1,ncomp 
      if(difu(ICOM)<0.d0) then 
        write(ifle,*) 
     &    'ERR: [diffusity(icom)] is not defined in USER_PROPERTY'
        call FFRABORT(1,'ERR: in USER_PROPERTY')
      else 
        rds(ICVL,:)=difu(:)*rhou
      endif
!
      if(ical_thmoDIFF==1) then 
        if(thmdif(ICOM)<0.d0) then 
          write(ifle,*) 
     &   'ERR: [thermo_diffusity(icom)] is not defined in USER_PROPERTY'
          call FFRABORT(1,'ERR: in USER_PROPERTY')
        else
          TH_DIFF(ICVL,:)=thmdif(:)*rhou
        endif
      endif
      enddo
!
!        call diffusivity_density(ICVL)    ! (mass diffusivity of 
                                           ! species)*(density) [kg/(m s)] 
      enddo
      end subroutine USER_PRO
!
!===========================================================================
      subroutine model_Simpli(do_start,do_end,prlmnr_t,tmp,cp,
     &             rmu,rmd,rds,
     &             IMAT,MXALLCV,nflud,MXcomp)
  !simplified transport model
!===========================================================================
      use module_material ,only : prlmnr
      use module_species,only : r_Le

      integer,intent(in)    :: do_start,do_end,IMAT,nflud,MXALLCV,
     &                         MXcomp
      real*8 ,intent(in)    :: tmp   (       MXALLCV)
      real*8 ,intent(in)    :: prlmnr_t(nflud)
      REAL*8 ,INTENT(IN)    :: cp    (       MXALLCV)
      real*8 ,intent(out)   :: rmu   (       MXALLCV)
      real*8 ,intent(out)   :: rmd   (       MXALLCV)
      real*8 ,intent(out)   :: rds   (       MXALLCV,MXcomp)
      
!
      real*8,parameter :: A=2.58d-5,dex=0.7D0 
      real*8,parameter :: T0=298.d0
      integer :: ICVL
!
      if(ical_vect) then
!CIDR NODEP 
        do ICVL=do_start,do_end
        dum1=A*(tmp(ICVL)/T0)**dex 
        W1K8(ICVL)=dum1
        rmu(ICVL)=dum1*prlmnr_t(IMAT) ! viscosity [Pa s] 
        rmd(ICVL)=dum1*cp(ICVL)       ! thermal conductivity [J/(m s K)] 
        enddo

        do icom=1,ncomp
        do ICVL=do_start,do_end
        rds(ICVL,icom)=W1K8(ICVL)*r_Le(icom)! (mass diffusivity of 
                                        ! species)*(density) [kg/(m s)] 
        enddo
        enddo
      else
        do ICVL=do_start,do_end
        dum1 = A*(tmp(ICVL)/T0)**0.7D0 
        rmu(ICVL) = dum1*prlmnr_t(IMAT) ! viscosity [Pa s] 
        rmd(ICVL) = dum1*cp(ICVL)       ! thermal conductivity [J/(m s K)] 
        do icom=1,ncomp 
        rds(ICVL,icom) = dum1*r_Le(icom)! (mass diffusivity of 
                                        ! species)*(density) [kg/(m s)] 
        enddo
        enddo
      endif
      end subroutine model_Simpli
!
!===========================================================================
      subroutine model_Mix_Kinetic(do_start,do_end)        ! constant viscosity
!===========================================================================
      use module_material ,only : sthrmu,r_prlmnr
      use module_species  ,only : eps_KK,sigmaa
      use module_metrix,only    : mu_i,k_i
      use module_metrix,only   : rcomp,TH_DIFF,r_wmm=>ys,wmm=>eva_comp,
     &                           Xi=>sinj_comp
!
      implicit none
!
!
      integer,intent(in) :: do_start,do_end
      integer imixkinetic,j
      real*8, parameter :: Amu=1.16145d0,Bmu=0.14874d0,Cmu=0.52487d0,
     &                     Dmu=0.77320d0,Emu=2.16178d0,Fmu=2.43787d0
      real*8, parameter :: AD=1.060636d0,BD=0.15610d0, CD=0.193d0,
     &                     DD=0.47635d0, ED=1.03587d0, FD=1.52996d0,
     &                     GD=1.76474d0, HD=3.89411d0
      real*8, parameter :: ASTR=1.11d0,BSTR=1.15D0,CSTR=0.86d0
      real*8, parameter :: k_B=1.3806503D-23
      real*8, parameter :: mu_c=2.6693d-6,k_c=1.9891D-4,Dij_c=1.8583D-2
      real*8, parameter :: mu_u=0.1d0,k_u=418.67979999999994d0
      real*8, parameter :: C_F1=0.511d0,C_F2=0.489D0,C_F3=0.659d0,
     &                     C_F=-2.59D-7  !-2.59D-8
      real*8 omega_k,t_star,
     &       dum_k,dum_mu,K_ij,dum1_flnt,dum2_flnt
      real*8 :: t_mw,dum1,dum2,dum3,dum4,dum5,dum6,dij,T,eps_KKij,
     &          Omega_mu,Omega_D
      Omega_mu(T)=Amu*T**(-Bmu)+Cmu*exp(-Dmu*T)+Emu*exp(-Fmu*T)
      Omega_D(T)=AD*T**(-BD)+CD*exp(-DD*T)+ED*exp(-FD*T)
     &                    +GD*exp(-HD*T)
!---------------------------------------------------------------------------
! --- 
!---------------------------------------------------------------------------
      wmm(:)=wm(:)*1.d3
!--------------------
! --- 
!--------------------
      do ICVL=do_start,do_end 
        rcomp(:)=yys(ICVL,:)
        dum1=0.d0
        do icom=1,ncomp
          dum1=dum1+rcomp(ICOM)*r_wm(ICOM)
        enddo
        dum1=1.d0/dum1
        do icom=1,ncomp
          Xi(ICOM)=(rcomp(ICOM)*r_wm(ICOM))*dum1
        enddo

        T_MW=0.d0
        do icom=1,ncomp
        t_star=tmp(ICVL)/eps_KK(ICOM)
! --- MIX kinetic Method :
        mu_i(icom)=mu_c*sqrt(wmm(icom)*tmp(ICVL))/Omega_mu(t_star)/
     &             sigmaa(icom)**2  !conv from g/cm-s to Pa-s
! --- MIX kinetic Method :
        k_i(icom)=k_c*sqrt(tmp(ICVL)/wmm(icom))/Omega_mu(t_star)/
     &             sigmaa(icom)**2*k_u  ! conv from [cal/s/cm/K] to [W/m/K]
! --- Eucken's Law
!        K_i(icom)=mu_i(icom)*(cps(ICVL,icom)+5.d0/4.d0*gascns/wmm(icom))
!
        K_i(icom)=mu_i(icom)*(cps(ICVL,icom)+5.d0/4.d0*gascns/wm(icom))
        T_MW=T_MW+Xi(icom)/wmm(icom)
        end do
!
        T_MW=1.d0/T_MW
!
        dum_mu=0.d0
        dum_k=0.d0
        DO ICOM=1,ncomp
          dum1=0.d0
          dum2=0.d0
          dum3=0.d0
          dum5=0.d0
          dum6=0.d0
          dum1_flnt=0.d0
          dum2_flnt=0.d0
          do J=1,ncomp

          dum1=dum1+Xi(J)/sqrt(8.d0)/sqrt(1.d0+wm(ICOM)/wm(j))
     &           *(1.d0+sqrt(k_i(ICOM)/k_i(j))*wm(j)/wm(ICOM))**2
          dum2=dum2+Xi(J)/sqrt(8.d0)/sqrt(1.d0+wm(ICOM)/wm(j))
     &  *(1.d0+sqrt(mu_i(ICOM)/mu_i(j))*(wm(j)/wm(ICOM))**0.25d0)**2
!
          if(j==ICOM) cycle
!--- Diffusion : Di
          eps_KKij=k_B*sqrt(eps_KK(ICOM)*eps_KK(J))
          t_star=k_B*tmp(ICVL)/eps_KKij
          dum4=(0.5d0*(sigmaa(ICOM)+sigmaa(J)))**2
     &        *Omega_D(t_star)*prs(ICVL)
          DIJ=Dij_c*tmp(ICVL)**(1.5)*
     &        sqrt(1.d0/wmm(ICOM)+1.d0/wmm(J))/(dum4+1.d-20)
!
          K_ij=1.d0
          dum3=dum3+K_ij*wmm(j)*wmm(ICOM)/T_MW**2*DIJ       !thermo-diffusion
!
          dum5=dum5+Xi(J)/(DIJ+SML)
          dum6=dum6+Xi(J)*wmm(J)
!
          dum1_flnt=dum1_flnt+wmm(J)**C_F1*Xi(J)
          dum2_flnt=dum2_flnt+wmm(J)**C_F2*Xi(J)
          enddo
!
          dum_mu=dum_mu+Xi(ICOM)*mu_i(ICOM)/(dum2+1.d-20) 
!          dum_k=dum_k+Xi(ICOM)*k_i(ICOM)/(dum1+1.d-20)
          dum_k=dum_k+Xi(ICOM)*k_i(ICOM)/(dum2+1.d-20)
!-----------------------------------------------------------
! --- [6] :
          rds(ICVL,icom)=dum6/T_MW/(dum5+1.d-20)*rho(ICVL)
!-----------------------------------------------------------
! --- [CFD-ACE] & [FLUENT] : 
!          rds(ICVL,icom)=(1.d0-Xi(icom))/(dum5+1.d-20)*rho(ICVL)
!-----------------------------------------------------------
          if(ical_thmoDIFF==1) then
! --- [6] :
!            TH_DIFF(ICVL,icom)=dum3*rho(ICVL)
! --- [2] FLUENT 

            TH_DIFF(ICVL,icom)=C_F*tmp(ICVL)**C_F3
     &     *(wmm(ICOM)**C_F1*Xi(ICOM)/dum1_flnt-rcomp(ICOM))
     &     *(dum1_flnt/dum2_flnt)
          endif 
        END DO 
        rmd(ICVL)=dum_k  
        rmu(ICVL)=dum_mu 
!
!        thermal conductivity [J/(m s K)]
!        call diffusivity_density(ICVL)          ! (mass diffusivity of 
                                                 ! species)*(density) [kg/(m s)]
      enddo
      end subroutine model_Mix_Kinetic
!
!
!===========================================================================
      subroutine diffusivity_density(ICVL) 
!===========================================================================
! --- (mass diffusivity of species)*(density) [kg/(m s)] 
      use module_material ,only : r_sclmnr
      integer,intent(in) :: ICVL     ! classification number of cell
      real*8,parameter   :: eps=1.d-20
      real*8             :: zs(mxcomp),zz
      zz=0.d0
      do icom=1,ncomp
        zs(icom)=yys(ICVL,icom)*r_wm(icom)
        zz=zz+zs(icom)
      enddo
      do icom=1,ncomp
        rds(ICVL,icom)=rmu(ICVL)*zz/sclmnr_t(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,zz-zs(icom))
                       ! species)*(density) [kg/(m s)]
      enddo
!
      end subroutine diffusivity_density
!
!===========================================================================
      subroutine cal_solid(do_start,do_end)   ! calculation of solid part 
!===========================================================================
      use module_material, only : rmdsld
      integer,intent(in) :: do_start,do_end
      do ICVL=do_start,do_end
        rmu(ICVL) = -1.d0           ! viscosity [Pa s]
        rmd(ICVL) = rmdsld(-IMAT)  ! thermal conductivity [J/(m s K)]
        !rmd(ICVL) = 0.001421639D0*tmp(ICVL)+0.9099076D0  ! Sumitomo version!!!yohnishi

        rds(ICVL,1:ncomp) = 0.d0   ! (mass diffusivity of species)*(density) 
                                   !                              [kg/(m s)]
      enddo
      end subroutine cal_solid
!
      end subroutine set_trnsprp
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine set_vof_t2cph(icode,
     &  LVEDGE,LBC_SSF,LCYCSF,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,wifsld,LCYCOLD,
     &  hhh,hhs,cps,cp,cr,
     &  aks,tmp,yys,rho,rmu,rmd,rds,rmx)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_io,only        : ifle,ifll
      use module_material ,only : isthr,suther,const,simpli,
     &                            sthrmu,sthrt,sthrc,
     &                            prlmnr,sclmnr,nflud,
     &                            sthrmu2,sthrt2,sthrc2,
     &                            prlmnr2,sclmnr2
      use module_species,  only : sw,gascns,acpk,acpk2,wm,chkncomp,
     &                            r_wm
      use module_material, only : rmdsld,rsld
      use module_initial , only : rho0,rho02
      use module_scalar,   only : ivof
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC),icode
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      REAL*8 ,INTENT(INOUT) :: hhh   (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: hhs   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cps   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(OUT)   :: cp    (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cr    (       MXALLCV)
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(in)    :: aks   (       MXALLCVR,mxrans)
      real*8 ,intent(in)    :: tmp   (       MXALLCV)
      real*8 ,intent(in)    :: yys   (       MXALLCV,MXcomp)
      real*8 ,intent(out)   :: rho   (       MXALLCV)
      real*8 ,intent(out)   :: rmu   (       MXALLCV)
      real*8 ,intent(out)   :: rmd   (       MXALLCV)
      real*8 ,intent(out)   :: rds   (       MXALLCV,MXcomp)
      real*8 ,intent(out)   :: rmx   (       MXALLCV)
!
! --- [local entities]
!
      real*8  :: zs(mxcomp),SML=1.D-15,eps=1.d-20
      real*8  :: tt,zz,dum1,dum2,dum3,dum4,dum5,dum6,dumA,dumB
      real*8  :: cp1,cp2,cps1,cps2,vof1,vof2,hhs1,hhs2,hhhx1,hhhx2
      real*8  :: ahk(0:5,mxcomp),ahk2(0:5,mxcomp)
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,ICOM,IDCS,IDCE,n
!
! --- 
!
      call chkncomp(ncomp,'(set_trnsprp)')
!
      do 100 ICOM=1,ncomp
        ahk(0,ICOM)=acpk(0,ICOM)
        do 105 n=1,5
        ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
  105   continue
  100 continue
      do 200 ICOM=1,ncomp
        ahk2(0,ICOM)=acpk2(0,ICOM)
        do 205 n=1,5
        ahk2(n,ICOM)=acpk2(n,ICOM)/dble(n)
  205   continue
  200 continue
!
!----------------
!-< 1. density >-
!----------------
!      do IIMAT=1,NMAT
!        if( .not.mat_cal(IIMAT) ) cycle
!        IMAT = MAT_NO(IIMAT)
!        ICVS = MAT_CVEXT(IIMAT-1)+1
!        ICVE = MAT_CVEXT(IIMAT)
!        IDCS = MAT_DCIDX(IIMAT-1)+1
!        IDCE = MAT_DCIDX(IIMAT)
!        if( IMAT>0 ) then       ! fluid
!          call cal_fluid(IIMAT,ICVS,ICVE)
!          call cal_fluid(IIMAT,IDCS,IDCE)
!        else                    ! solid
!          do ICVL=ICVS,ICVE
!            rho(ICVL) = rsld(-IMAT)
!          enddo
!        endif
!      enddo
!----------------------------------------------
!--< 1.1 specific heat at constant pressure >--
!----------------------------------------------
      if(sw) then
        WRITE(ifle,*) ' ### ERR : Not support a7 model for VOF '
        call FFRABORT(1,'set_vofprp')
      endif
!---------------------
!--< 1.2 viscosity >--
!---------------------
      do IIMAT=1,NMAT
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        if(IMAT>0) then                 ! fluid
          select case(isthr(IMAT))
            case(suther)                ! Sutherland's formula
              call model_Sutherland(IIMAT,IMAT,ICVS,ICVE)
            case(const)                 ! constant viscosity
              call model_Constant(IIMAT,IMAT,ICVS,ICVE)
            case(simpli)                ! simplified transport model
              WRITE(ifle,*) 
     &        ' ### ERR :  simplified transport model ',
     &        'NOT supported for VOF '
              call FFRABORT(1,'set_vofprp')
          end select
        else                    ! solid
          call cal_solid(IMAT,ICVS,ICVE)
        endif
      enddo
!-----------------------------
!-< 3. Scale of diffusivity >-
!-----------------------------
      if(icode.eq.1) then
      do IIMAT=1,NMAT
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        if( IMAT>0 ) then       ! fluid
          do ICVL=ICVS,ICVE
            rmx(ICVL) = max(rmu(ICVL),rmd(ICVL)/(cp(ICVL)+SML))
            do icom=1,ncomp
            rmx(ICVL) = max(rmx(ICVL),rds(ICVL,icom))
            enddo
          enddo
        else                    ! solid
          rmx(ICVS:ICVE) = rmd(ICVS:ICVE)/cp(ICVS:ICVE)
        endif
      enddo
      endif
!---------------------------------
!-< 4. Set value at dummy cells >-
!---------------------------------
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmu)
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmd)
!
      do ICOM=1,ncomp
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rds(:,ICOM))
      enddo
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmx)
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,cp)
!
      return
!
!///////////////////////////////////////////////////////////////////////
      contains
!==============================================================================
      subroutine model_Sutherland(IIMAT,IMAT,do_start,do_end)  ! Sutherland's formula
!==============================================================================
      use module_material ,only : sthrmu,sthrt,sthrc,r_prlmnr
      integer,intent(in) :: do_start,do_end,IMAT,IIMAT
!
      dumA=sthrmu(IMAT)*(sthrt(IMAT)+sthrc(IMAT))
      dumB=sthrmu2(IMAT)*(sthrt2(IMAT)+sthrc2(IMAT))
      do ICVL=do_start,do_end
        vof1=aks(ICVL,ivof)
        vof2=1.d0-aks(ICVL,ivof)
        dum1=
     &  dumA/(tmp(ICVL)+sthrc(IMAT))*sqrt((tmp(ICVL)/sthrt(IMAT))**3)
        dum2=
     &  dumB/(tmp(ICVL)+sthrc2(IMAT))*sqrt((tmp(ICVL)/sthrt2(IMAT))**3)
        rmu(ICVL)=dum1*dum2/(dum2*vof1+dum1*vof2) ! viscosity [Pa s]
        tt=tmp(ICVL)
        cp1=0.d0
        cp2=0.d0
        zz =0.d0
        hhhx1=0.d0
        hhhx2=0.d0
        do 100 icom=1,ncomp
        hhs1=((((
     &      ahk(5,ICOM) *tt
     &     +ahk(4,ICOM))*tt
     &     +ahk(3,ICOM))*tt
     &     +ahk(2,ICOM))*tt
     &     +ahk(1,ICOM))*tt
     &     +ahk(0,ICOM)
        hhhx1=hhhx1+hhs1*yys(ICVL,ICOM)
        cps1=(((
     &      acpk(5,ICOM) *tt
     &     +acpk(4,ICOM))*tt
     &     +acpk(3,ICOM))*tt
     &     +acpk(2,ICOM))*tt
     &     +acpk(1,ICOM)
        cp1=cp1+yys(ICVL,ICOM)*cps1
!
        hhs2=((((
     &      ahk2(5,ICOM) *tt
     &     +ahk2(4,ICOM))*tt
     &     +ahk2(3,ICOM))*tt
     &     +ahk2(2,ICOM))*tt
     &     +ahk2(1,ICOM))*tt
     &     +ahk2(0,ICOM)
        hhhx2=hhhx2+hhs2*yys(ICVL,ICOM)
        cps2=(((
     &      acpk2(5,ICOM) *tt
     &     +acpk2(4,ICOM))*tt
     &     +acpk2(3,ICOM))*tt
     &     +acpk2(2,ICOM))*tt
     &     +acpk2(1,ICOM)
        cp2=cp2+yys(ICVL,ICOM)*cps2
!
        hhs(ICVL,ICOM)=hhs1*vof1+hhs2*vof2
        cps(ICVL,ICOM)=cps1*vof1+cps2*vof2
        zs(icom) = yys(ICVL,icom)*r_wm(icom)
        zz = zz+zs(icom)
 100    continue
        hhh(ICVL)=hhhx1*vof1+hhhx2*vof2
        cp(ICVL)=(cp1*vof1*rho0(IIMAT)+cp2*vof2*rho02(IIMAT))/rho(ICVL)
        dum3=cp1*dum1/prlmnr(IMAT)     
        dum4=cp2*dum2/prlmnr2(IMAT)
        rmd(ICVL)=dum3*dum4/(dum4*vof1+dum3*vof2) ! thermal conductivity [J/(m s K)]
        do icom=1,ncomp
        dum5=dum1*zz/prlmnr(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,zz-zs(icom))
        dum6=dum2*zz/prlmnr2(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,zz-zs(icom))
        rds(ICVL,icom)=dum5*dum6/(dum6*vof1+dum5*vof2)! (mass diffusivity of 
                                                      ! species)*(density) [kg/(m s)]
        
        enddo
      enddo
      end subroutine model_Sutherland
!===========================================================================
      subroutine model_Constant(IIMAT,IMAT,do_start,do_end)    ! constant viscosity
!===========================================================================
      use module_material ,only : sthrmu,r_prlmnr
      integer,intent(in) :: do_start,do_end,IMAT,IIMAT
      dumA=sthrmu(IMAT)
      dumB=sthrmu2(IMAT)
      dum1=sthrmu(IMAT)
      dum2=sthrmu2(IMAT)
      do ICVL=do_start,do_end
        vof1=aks(ICVL,ivof)
        vof2=1.d0-aks(ICVL,ivof)
!        rmu(ICVL)=dumA*dumB/(dumB*vof1+dumA*vof2) ! viscosity [Pa s]
        rmu(ICVL)=vof1*dumA+vof2*dumB
        tt=tmp(ICVL)
        cp1=0.d0
        cp2=0.d0
        zz =0.d0
        do 100 icom=1,ncomp
        hhhx1=0.d0
        hhhx2=0.d0
        hhs1=((((
     &      ahk(5,ICOM) *tt
     &     +ahk(4,ICOM))*tt
     &     +ahk(3,ICOM))*tt
     &     +ahk(2,ICOM))*tt
     &     +ahk(1,ICOM))*tt
     &     +ahk(0,ICOM)
        hhhx1=hhhx1+hhs1*yys(ICVL,ICOM)
        cps1=(((
     &      acpk(5,ICOM) *tt
     &     +acpk(4,ICOM))*tt
     &     +acpk(3,ICOM))*tt
     &     +acpk(2,ICOM))*tt
     &     +acpk(1,ICOM)
        cp1=cp1+yys(ICVL,ICOM)*cps1
!
        hhs2=((((
     &      ahk2(5,ICOM) *tt
     &     +ahk2(4,ICOM))*tt
     &     +ahk2(3,ICOM))*tt
     &     +ahk2(2,ICOM))*tt
     &     +ahk2(1,ICOM))*tt
     &     +ahk2(0,ICOM)
        hhhx2=hhhx2+hhs2*yys(ICVL,ICOM)
        cps2=(((
     &      acpk2(5,ICOM) *tt
     &     +acpk2(4,ICOM))*tt
     &     +acpk2(3,ICOM))*tt
     &     +acpk2(2,ICOM))*tt
     &     +acpk2(1,ICOM)
        cp2=cp2+yys(ICVL,ICOM)*cps2
!
        hhs(ICVL,ICOM)=hhs1*vof1+hhs2*vof2
        cps(ICVL,ICOM)=cps1*vof1+cps2*vof2
        zs(icom) = yys(ICVL,icom)*r_wm(icom)
        zz = zz+zs(icom)
 100    continue
        hhh(ICVL)=hhhx1*vof1+hhhx2*vof2
        cp(ICVL)=(cp1*vof1*rho0(IIMAT)+cp2*vof2*rho02(IIMAT))/rho(ICVL) 
        dum3=cp1*dum1/prlmnr(IMAT)  
        dum4=cp2*dum2/prlmnr2(IMAT) 
        rmd(ICVL)=dum3*dum4/(dum4*vof1+dum3*vof2) ! thermal conductivity [J/(m s K)] 
        do icom=1,ncomp
        dum5=dum1*zz/prlmnr(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,zz-zs(icom))
        dum6=dum2*zz/prlmnr2(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,zz-zs(icom))
        rds(ICVL,icom)=dum5*dum6/(dum6*vof1+dum5*vof2+SML)! (mass diffusivity of 
        enddo                                         ! species)*(density) [kg/(m s)] 
      enddo
      end subroutine model_Constant
!
!========================================================================================
      subroutine cal_solid(IMAT,do_start,do_end)             ! calculation of solid part
!========================================================================================
      use module_material, only : rmdsld
      integer,intent(in) :: do_start,do_end,IMAT
      do ICVL=do_start,do_end
        rmu(ICVL) = 1.d0                   ! viscosity [Pa s]
        rmd(ICVL) = rmdsld(-IMAT)          ! thermal conductivity [J/(m s K)]
        rds(ICVL,1:ncomp) = 1.d0           ! (mass diffusivity of species)*(density) [kg/(m s)]
      enddo
      end subroutine cal_solid
!
!=================================================
      subroutine cal_fluid(IIMAT,do_start,do_end)
!=================================================
      integer,intent(in) :: do_start,do_end,IIMAT
      do ICVL=do_start,do_end  ! [kg/m^3]
      rho(ICVL)=aks(ICVL,ivof)*rho0(IIMAT)
     &         +(1.d0-aks(ICVL,ivof)*rho02(IIMAT))
      enddo
      end subroutine cal_fluid
!
      end subroutine set_vof_t2cph
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine set_trnsprp_t2h_cavi(icode,
     &  LVEDGE,LBC_SSF,LCYCSF,
     &  MAT_CV,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,wifsld,LCYCOLD,
     &  hhh,hhs,cps,cp,cr,
     &  tmp,yys,aks,rho,rmu,rmd,rds,rmx)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments] 
!
      use module_dimension
      use module_io,only        : ifle,ifll
      use module_material ,only : isthr,suther,const,simpli,
     &                            sthrmu,sthrt,sthrc,
     &                            prlmnr,sclmnr,nflud,
     &                            sthrmu2,sthrt2,sthrc2,
     &                            prlmnr2,sclmnr2
      use module_species,  only : wm,chkncomp,sw,acpk,acpk2,wm,r_wm
      use module_material, only : rmdsld
      use module_Euler2ph,only  : ieul2ph
      use module_model,only     : ical_vect,nthrds
      use module_vector,only    : ICVS_V,ICVE_V,
     &                            ICFS_V,ICFE_V,
     &                            ICVSIN_V,ICVEIN_V,
     &                            IDCS_V,IDCE_V,index_c,index_f
      use module_scalar,only    : ical_cavi,Tref,ivold,
     &                            icavi,T0_ref
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: LVEDGE    (2, MXCVFAC),icode
      integer,intent(in)    :: LBC_SSF   (   MXSSFBC)
      integer,intent(in)    :: LCYCSF    (   MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV    (   MXALLCV)
      INTEGER,INTENT(IN)    :: MAT_CVEXT (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX (   0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO    (   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal   (   0:MXMAT)
      REAL*8 ,INTENT(INOUT) :: hhh   (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: hhs   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cps   (       MXALLCV,MXcomp)
      REAL*8 ,INTENT(INOUT) :: cp    (       MXALLCV)
      REAL*8 ,INTENT(INOUT) :: cr    (       MXALLCV)
      real*8 ,intent(in)    :: tmp   (       MXALLCV)
      real*8 ,intent(in)    :: yys   (       MXALLCV,MXcomp)
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(out)   :: rmu   (       MXALLCV)
      real*8 ,intent(out)   :: rmd   (       MXALLCV)
      real*8 ,intent(out)   :: rds   (       MXALLCV,MXcomp)
      real*8 ,intent(out)   :: rmx   (       MXALLCV)
      real*8 ,intent(in)    :: aks   (       MXALLCVR,mxrans)
      real*8 ,intent(in)    :: rho   (       MXALLCV)
!
! --- [local entities]
!
      real*8  :: zs(mxcomp),SML=1.D-15,eps=1.d-20
      real*8  :: tt,zz,dum1,dum2,dum3,dum4,dum5,dum6,dumA,dumB
      real*8  :: cp1,cp2,cps1,cps2,vof1,vof2,hhs1,hhs2,hhhx1,hhhx2
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,ICOM,myid,n,IDCS,IDCE
      real*8,parameter :: A=2.58d-5
             ! constant for simplified transport model [kg/(m*s)]
      real*8,parameter :: T0=298.d0
             ! reference temperature for simplified transport model [K]
      real*8  :: ahk(0:5,mxcomp),ahk2(0:5,mxcomp)
      
!
      call chkncomp(ncomp,'(set_trnsprp)')
!
      do 100 ICOM=1,ncomp
      ahk(0,ICOM)=acpk(0,ICOM)
      do 105 n=1,5
      ahk(n,ICOM)=acpk(n,ICOM)/dble(n)
  105 enddo
  100 enddo
      do 200 ICOM=1,ncomp
      ahk2(0,ICOM)=acpk2(0,ICOM)
      do 205 n=1,5
      ahk2(n,ICOM)=acpk2(n,ICOM)/dble(n)
  205 enddo
  200 enddo
!
      if(sw) then
        WRITE(ifle,*) ' ### ERR : Not support a7 model for CAVITATION'
        call FFRABORT(1,'set_vofprp')
      endif

!------------------------------------------
!--< 1.2 viscosity >-- Scalar computation
!------------------------------------------
      do IIMAT=1,NMAT
        IMAT = MAT_NO(IIMAT)
        ICVS = MAT_CVEXT(IIMAT-1)+1
        ICVE = MAT_CVEXT(IIMAT)
        IDCS=MAT_DCIDX(IIMAT-1)+1
        IDCE=MAT_DCIDX(IIMAT)
        if(IMAT>0) then       ! fluid
          call model_cavit(ICVS,ICVE)
          call model_cavit(IDCS,IDCE)
        else                  ! solid
          call cal_solid(ICVS,ICVE)
          call cal_solid(IDCS,IDCE)
        endif
      enddo
!
!-< 3. Scale of diffusivity >-
!
      if(icode.eq.1) then
      do IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT>0) then       ! fluid
          do ICVL=ICVS,ICVE
            rmx(ICVL)=max(rmu(ICVL),rmd(ICVL)/(cp(ICVL)+SML))
            do icom=1,ncomp
            rmx(ICVL)=max(rmx(ICVL),rds(ICVL,icom))
            enddo
          enddo
        else                    ! solid
          rmx(ICVS:ICVE)=rmd(ICVS:ICVE)/cp(ICVS:ICVE)
        endif
      enddo
      endif
!
!-< 4. Set value at dummy cells >-
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmu)
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmd)
!
      do ICOM=1,ncomp
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rds(:,ICOM))
      enddo
!
      call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,rmx)
!
      return
!
!///////////////////////////////////////////////////////////////////////
      contains
!===========================================================================
      subroutine model_cavit(do_start,do_end)    
!===========================================================================
      use module_material ,only : sthrmu,r_prlmnr
      
      integer,intent(in) :: do_start,do_end
      real*8  :: alpha1,alpha2,SMR1,SMR2,XMYL,alpha,
     &           XMYA,SSD,T0SD,T,mu_ll,mu_vv,T2
!
      mu_ll=sthrmu(IMAT)
      mu_vv=sthrmu2(IMAT)
      dum1=sthrmu(IMAT)
      dum2=sthrmu2(IMAT)
      do ICVL=do_start,do_end
!
      T=tmp(ICVL)
      T0SD  = 416.7D0
      SSD   = 861.D0
      XMYA  = mu_vv*(T/T0SD)**1.5D0*(T0SD+SSD)/(T+SSD)
!
      T2    = Tref-273.15D0  ! Kelvin ==> Cerucus degree
      SMR1  = 0.03368D0
      SMR2  = 0.00022099D0
      XMYL  = mu_ll*(1.D0+SMR1*T2+SMR2*T2**2)**(-1.D0)

      
      alpha = aks(ICVL,ivold)
      rmu(ICVL)=(1.d0-alpha)*(1.d0+2.5d0*alpha)*XMYL+alpha*XMYA
      
      rmu(ICVL)=(1.d0-alpha)*(1.d0+2.5d0*alpha)*mu_ll+alpha*mu_vv
      
      dum1=cp(ICVL)*XMYL/prlmnr(IMAT)
      dum2=cp(ICVL)*XMYA/prlmnr2(IMAT)
      rmd(ICVL)=(1.d0-alpha)*(1.d0+2.5d0*alpha)*dum1+alpha*dum2
!
      alpha1=1.d0-aks(ICVL,ivold)
      alpha2=aks(ICVL,ivold)
      tt=tmp(ICVL)
      cp1=0.d0
      cp2=0.d0
      zz =0.d0
      hhhx1=0.d0
      hhhx2=0.d0
      do 100 icom=1,ncomp
        hhs1=((((
     &      ahk(5,ICOM) *tt
     &     +ahk(4,ICOM))*tt
     &     +ahk(3,ICOM))*tt
     &     +ahk(2,ICOM))*tt
     &     +ahk(1,ICOM))*tt
     &     +ahk(0,ICOM)
        hhhx1=hhhx1+hhs1*yys(ICVL,ICOM)
        cps1=(((
     &      acpk(5,ICOM) *tt
     &     +acpk(4,ICOM))*tt
     &     +acpk(3,ICOM))*tt
     &     +acpk(2,ICOM))*tt
     &     +acpk(1,ICOM)
        cp1=cp1+yys(ICVL,ICOM)*cps1
!
        hhs2=((((
     &      ahk2(5,ICOM) *tt
     &     +ahk2(4,ICOM))*tt
     &     +ahk2(3,ICOM))*tt
     &     +ahk2(2,ICOM))*tt
     &     +ahk2(1,ICOM))*tt
     &     +ahk2(0,ICOM)
        hhhx2=hhhx2+hhs2*yys(ICVL,ICOM)
        cps2=(((
     &      acpk2(5,ICOM) *tt
     &     +acpk2(4,ICOM))*tt
     &     +acpk2(3,ICOM))*tt
     &     +acpk2(2,ICOM))*tt
     &     +acpk2(1,ICOM)
        cp2=cp2+yys(ICVL,ICOM)*cps2
!
        hhs(ICVL,ICOM)=hhs1*alpha1+hhs2*alpha2
        cps(ICVL,ICOM)=cps1*alpha1+cps2*alpha2
        zs(icom) = yys(ICVL,icom)*r_wm(icom)
        zz = zz+zs(icom)
 100  enddo
      hhh(ICVL)=hhhx1*alpha1+hhhx2*alpha2
      cp(ICVL)=(cp1*alpha1+cp2*alpha2)
      do icom=1,ncomp
      dum5=dum1*zz/prlmnr(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,zz-zs(icom))
      dum6=dum2*zz/prlmnr2(IMAT)
     &          *(1.d0-yys(ICVL,icom))/max(eps,zz-zs(icom))
      rds(ICVL,icom)=dum5*dum6/(dum6*alpha1+dum5*alpha2+SML)
                                        ! (mass diffusivity of 
      enddo
      enddo
      end subroutine model_cavit
!
!===========================================================================
      subroutine cal_solid(do_start,do_end)   ! calculation of solid part 
!===========================================================================
      use module_material, only : rmdsld 
      integer,intent(in) :: do_start,do_end 
      do ICVL=do_start,do_end 
        rmu(ICVL) = 1.d0           ! viscosity [Pa s] 
        rmd(ICVL) = rmdsld(-IMAT)  ! thermal conductivity [J/(m s K)] 
        rds(ICVL,1:ncomp) = 0.d0   ! (mass diffusivity of species)*(density) 
                                   !                              [kg/(m s)] 
      enddo
      end subroutine cal_solid
!
      end subroutine set_trnsprp_t2h_cavi
