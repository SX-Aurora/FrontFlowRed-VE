!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine USER_PROPERTY(ICVL,IMAT_U,ncomp,wm,
     &           T_MW,temp,prs,rho,Ys,Xi,cp_i,cp,
     &           mu_usr,rmd_usr,diffusivity,thmo_diffusivity) 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      use module_species  ,only : eps_KK,sigmaa
!     Viscosity            :  (Ps-s)
!     Thermal Conductivity :  (w/m-k)
!     Diffusivity          :  (m^2/s)
!     thmo_diffusivity     :  ()
!
!
!
!---------------------------------------------------------------------
      IMPLICIT NONE 
!
! --- [DUMMY ARGUMENTS] 
!
      integer,intent(in)      :: ICVL,IMAT_U,ncomp 
      real*8 ,intent(in)      :: wm(ncomp)
      real*8 ,intent(in)      :: T_MW
      real*8 ,intent(in)      :: temp 
      real*8 ,intent(in)      :: prs 
      real*8 ,intent(in)      :: rho 
      real*8 ,intent(in)      :: Ys(ncomp)
      real*8 ,intent(in)      :: Xi(ncomp)
      real*8 ,intent(in)      :: cp_i(ncomp)
      real*8 ,intent(in)      :: cp
      real*8 ,intent(inout)   :: mu_usr
      real*8 ,intent(inout)   :: rmd_usr
      real*8 ,intent(inout)   :: diffusivity(ncomp)
      real*8 ,intent(inout)   :: thmo_diffusivity(ncomp)
!
! --- [LOCAL ENTITIES]
!
      integer :: ICOM,J
      real*8, parameter :: C_F1=0.511d0,C_F2=0.489D0,C_F3=0.659d0,
     &                     C_F=-2.59D-7
      real*8, parameter :: k_B=1.3806503D-23,Dij_c=1.8583D-2
      real*8, parameter :: AD=1.060636d0,BD=0.15610d0, CD=0.193d0,
     &                     DD=0.47635d0, ED=1.03587d0, FD=1.52996d0,
     &                     GD=1.76474d0, HD=3.89411d0,
     & SML=1.d-20
      
      real*8 :: eps_KK(2)!=207.6,eps_KK(2)=38
      real*8 :: sigmaa(2)!=4.084,sigmaa(2)=2.92

!
! --- user array 
!
      real*8 :: dum1_flnt,dum2_flnt,zz,zs(2),dum4,dum5,dum6,t_star,
     &          eps_KKij,Omega_D,T,DIJ,wmm(2)
!     
      Omega_D(T)=AD*T**(-BD)+CD*exp(-DD*T)+ED*exp(-FD*T)
     &                    +GD*exp(-HD*T)
!
      eps_KK(1)=207.6
      eps_KK(2)=38
      sigmaa(1)=4.084
      sigmaa(2)=2.92
! --- 
!     
      wmm(:)=wm(:)*1.d3
      dum1_flnt=0.d0 
      dum2_flnt=0.d0 
      do J=1,ncomp 
      dum1_flnt=dum1_flnt+wm(J)**C_F1*Xi(J) 
      dum2_flnt=dum2_flnt+wm(J)**C_F2*Xi(J) 
      enddo 
! 
      mu_usr=3.17d-6+2.04D-8*temp-3.52D-12*temp**2 
      rmd_usr=3.8D-2+5.41D-4*temp
     &           -2.51D-7*temp**2+8.57D-11*temp**3 
!
!
!------------------------------------------------------------------
      do icom=1,ncomp
      dum6=0.d0    
      dum5=0.d0   
      do J=1,ncomp
      if(J==icom) cycle
      eps_KKij=k_B*sqrt(eps_KK(ICOM)*eps_KK(J))
      t_star=k_B*temp/eps_KKij
      dum4=(0.5d0*(sigmaa(ICOM)+sigmaa(J)))**2
     &        *Omega_D(t_star)*prs
      DIJ=Dij_c*temp**(1.5)*
     &        sqrt(1.d0/wmm(ICOM)+1.d0/wmm(J))/(dum4+1.d-20)
      dum6=dum6+Xi(J)*wmm(J)
      dum5=dum5+Xi(J)/(DIJ+SML)
      enddo
      diffusivity(icom)=dum6/(T_MW*1.d3)/(dum5+1.d-20)!*rho
!      diffusivity(icom)=(1.d0-Xi(icom))/(dum5+1.d-20)!*rho(ICVL)
      enddo
!------------------------------------------------------------------
!diffusivity(ICOM)
!      zz=0.d0
!      do icom=1,ncomp
!      zs(icom)=ys(icom)/wm(icom)
!      zz=zz+zs(icom)
!      enddo
!
!      do icom=1,ncomp
!      diffusivity(ICOM)=mu_usr*zz/0.72
!     &     *(1.d0-ys(icom))/max(1.d-20,zz-zs(icom))/rho
!      ! species)*(density) [kg/(m s)]
!      enddo
!------------------------------------------------------------------
!      do ICOM=1,ncomp 
!      diffusivity(ICOM)=7.234D-8*temp+4.659D-10*temp**2 
!     &           -8.016d-14*temp**3
!
!      thmo_diffusivity(ICOM)=C_F*temp**C_F3
!     &     *(wm(ICOM)**C_F1*Xi(ICOM)/dum1_flnt-Ys(ICOM)) 
!     &     *(dum1_flnt/dum2_flnt) 
!      enddo
!
      return
      end subroutine USER_PROPERTY
