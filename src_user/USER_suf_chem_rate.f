!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_suf_chem_rate(
     &           no_bc,IRC,BC_name,ncomp,ncomp_suf,ugc,wm,Ps,Temp,
     &           dens,
     &           stk_a,stk_b,stk_c,
     &           Mass_Fraction,
     &           area,vol,
     &           mol_fraction,
     &           mol_concentration,
     &           stoich_coef_for,
     &           stoich_coef_rev,
     &           stkrat)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- CALCULATE USER Surface REACTION RATE ---------------------------
!
! Vables coming into this subroutine:
!
!     no_bc       : No. of CVD wall BC in fflow.ctl file
!     BC_name     : CVD wall BC name
!
!     ncomp       : Total Species number of GAS 
!                   defined at [&species] in fflow.ctl
!
!     ncomp_suf   : Total Species number of Surface (SITE and BULK)
!                   defined at [&&surface_species] in fflow.ctl  
!
!     ugc         : R=Univeral Gas Constant (R=8.314d0) [J/mol/K]
!     wm(ncomp)   : Mi=molecular weight [kg/mol]
!
!                 ! Probability function : 
!                 ! Si=A*T^alpha*P^beta*exp(-E/R/T)
!                 ! a     : Prefactor defined     at &chemreac in fflow.ctl
!                 ! alpha : Exponent of T defined at &chemreac in fflow.ctl
!                 ! e     : Activation defined    at &chemreac in fflow.ctl
!     stk_a       : stk_a => a
!     stk_b       : stk_b => alpha
!     stk_c       : stk_c => e
!
!     Ps(ncomp)   : Pi=Partial pressure [Pa] of each species
!     Temp        : T=Temperature [K] of 1st CV
!     mol_fraction     : mole fraction (%)
!     mol_concentration: mole concentration (mol/m^3, or mol/m^2)
!     stoich_coef_for  : forward reaction stoichometric coefficients (v')
!     stoich_coef_rev  : reverse reaction stoichometric coefficients (v'')
!     
! Vables returning main routine
!     stkrat      : stiking rate coefficient ;[mol/m^2/s]
!     
!                 ! Heltz-Knudssen surface sticking model:
!                 ! stkrat==Ri=Si*Pi/sqrt(2*pi*Mi*R*T);
!
!                 ! stkrat > 0 => production
!                 !        < 0 => disappearance
!
! Attention : Stoichometric coefficients should be multipled 
!             onto [stkrat] for eq.=IRC in main-code as below
!             to calculate "stick rate [mol/m^2/s]"
!             --------------------------
!             do icom=1,ncomp+ncomp_suf
!             stkrat(icom)=stkrat(icom)*(stoich_coef_rev(icom)-stoich_coef_for(icom))
!             enddo                           
!             -------------------------
!
! Example   : standard sticking model's in main code
!             --------------------------
!             dum1=stk_a(IRC)*Temp**stk_b(IRC)*exp(-stk_c(IRC)/ugc/Temp)
!             Si=max(min(dum1,1.d0),-1.d0)
!             dum2=
!    &        dum1*mol_concentration(icom_gas(IRC))**stoich_coef_for(icom_gas(IRC))
!             stkrat(icom)=dum2
!             --------------------------
!             here :
!                [icom_gas(IRC)] is a only-one gas species number in fflow.ctl
!             which takes part in the reaction of IRC, 
!             because sticking model just allow exactly-only one gas phase 
!             species as reactant, however other reactants can be SITE species
!---------------------------------------------------------------------
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      integer,intent(in)      :: no_bc
      integer,intent(in)      :: IRC
      integer,intent(in)      :: ncomp,ncomp_suf
      character(*),intent(in) :: BC_name
      real*8 ,intent(in)      :: ugc,Temp,wm(ncomp+ncomp_suf)
      real*8 ,intent(in)      :: Ps(ncomp+ncomp_suf)
      real*8 ,intent(in)      :: dens,Mass_Fraction(ncomp),area,vol
      real*8 ,intent(in)      :: mol_fraction(ncomp+ncomp_suf)
      real*8 ,intent(in)      :: mol_concentration(ncomp+ncomp_suf)
      real*8 ,intent(in)      :: stoich_coef_for(ncomp+ncomp_suf)
      real*8 ,intent(in)      :: stoich_coef_rev(ncomp+ncomp_suf)
      real*8 ,intent(in)      :: stk_a,stk_b,stk_c
      real*8 ,intent(out)     :: stkrat(ncomp+ncomp_suf)
!
! --- [LOCAL ENTITIES]
!
      real*8,parameter  :: pi=3.1415926d0,SML=1.d-20
      integer:: icom
      real*8 :: dum1,dum2,dum3
!------------------------------------------
! --- One step reaction for SiH4=>Si+2H2
!------------------------------------------
!      real*8,parameter  :: k1=1.75d3,k2=4.d4,k3=10.d0,Patoat=9.8692D-6
!      integer,parameter :: icom_SiH4=1,icom_H2=2,
!     &                     icom_Si=3
!      real*8 :: k,T,P_atm(3)
!
!------------------------------------------
! --- 3 step reaction :
!------------------------------------------
      integer,parameter :: icom_SiH4=1,icom_H2=2,
     &                     icom_SiH2=3,icom_Si=4
      real*8,parameter  :: KB=7.826D0
      real*8 :: gamma,TT
!
!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- One step reaction for SiH4(g)=>Si(s)+2H2(g)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!--------------------------
! --- Reaction eq. No=1 :
!--------------------------
!      K(T)=1.25D9*exp(-18500.d0/T)
!      stkrat(:)=0.d0
!      P_atm(1:ncomp)=Ps(1:ncomp)*Patoat
!      if(IRC==1) then
!        dum1=1.d0+k1*P_atm(icom_H2)+k2*P_atm(icom_SiH4)
!        dum2=k3*K(Temp)*P_atm(icom_SiH4)/dum1
!        stkrat(icom_Si)  =dum2
!        stkrat(icom_SiH4)=dum2
!        stkrat(icom_H2)  =dum2
!      endif
!      return
!
!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- 3 step reaction :
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      gamma(TT)=0.0537d0*exp(-9400.d0/TT) 
      stkrat(:)=0.d0
!-----------------------------------------------
! --- Reaction eq. No=2 :SiH4(g)=>Si(s)+2H2(g)
!-----------------------------------------------
      if(IRC==2) then
        dum1=gamma(Temp)*mol_concentration(icom_SiH4)
     &      *sqrt(KB*Temp/(2.d0*pi*wm(icom_SiH4)))
        stkrat(icom_SiH4)=dum1
        stkrat(icom_Si)=dum1
        stkrat(icom_H2)=dum1
      endif
!-----------------------------------------------
! --- Reaction eq. No=3 :SiH2(g)=>Si(s)+H2(g)
!-----------------------------------------------
      if(IRC==3) then
        dum1=0.1d0*mol_concentration(icom_SiH2)
     &      *sqrt(KB*Temp/(2.d0*pi*wm(icom_SiH2)))
        stkrat(icom_SiH2)=dum1
        stkrat(icom_Si)=dum1
        stkrat(icom_H2)=dum1
      endif
!---------------------------------------------------------------------
!
      return
      end subroutine user_suf_chem_rate
