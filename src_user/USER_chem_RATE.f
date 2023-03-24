!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE USER_CHEM_RATE(
     &                     MXMAT,MXCOMP,MXALLCV,NMAT,NCOMP,nflud,
     &                     mxcomp_suf,ncomp_suf,
     &                     nofld,MAT_NO,MAT_CVEXT,MAT_DCIDX,CVCENT,
     &                     IRC,nnrc,vreq,kind_chem,spcnam,wm,ugc,
     &                     moleX,yi,tmp,dumya,dumyb,wdot)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- CALCULATE USER REACTION RATE of GAS Phase ----------------------
!
!---------------------------------------------------------------------
! --- DIMENSION 
!---------------------------------------------------------------------
!
!     MXMAT    : material number (array's dimension)
!     MXCOMP   : species number  (array's dimension)
!     MXALLCV  : CV number       (array's dimension)
!     NMAT     : material number (do-loop up-limit)
!     NCOMP    : species number  (do-loop up-limit)
!     nflud    : fluid number    
!
!---------------------------------------------------------------------
! --- CONNECTIVITY
!---------------------------------------------------------------------
!
!     MAT_CVEXT (0:MXMAT)     : Pointer for CV's do-loop
!     MAT_DCIDX (0:MXMAT)     : Pointer for dummy-CV's do-loop
!     MAT_NO    (0:MXMAT)     : Pointer for material no.
!     nofld     (  nflud)     : Pointer for fluid-material no
!     CVCENT    (3,MXALLCV)   : CV center coordinates 
!
!---------------------------------------------------------------------
! --- PHYSICAL VABLES
!---------------------------------------------------------------------
!
!     moleX (mxallcv,mxcomp)  : mol concentration . [mol/m^3]
!     yi    (mxallcv,mxcomp)  : mass fraction 
!     tmp   (mxallcv)         : temprature    [k]
!     wdot  (mxallcv,mxcomp)  : net value of reaction rate [kg/(m^3*s)]
!     dumya (mxallcv)         : warking array
!     dumyb (mxallcv)         : warking array
!
!---------------------------------------------------------------------
! --- Reaction 
!---------------------------------------------------------------------
!
!     IRC              : reaction sequence no. in fflow.ctl
!     nnrc             : total reaction number.
!     vreq(icom,i,IRC) : stoichiometric coefficients.
!                        icom : species sequence no. in fflow.ctl.
!                        i    : =1  =>  left side of reaction equation IRC.
!                               =2  =>  right side of reaction equation IRC.
!                        IRC  : reaction sequence no. in fflow.ctl.
!     kind_chem(IRC)   : kind of reaction. (gas-phase or surface)
!                        =1 => gas-phase.
!                        =2 => surface.
!     spcnam(icom)     : species name.
!     wm(icom)         : molecular weight of icom. [kg/mol]
!     ugc              : univeral gas constant.{=8.314d0[j/mol/k]}
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!
! --- [DUMMY ARGUMENTS]
!
      INTEGER,INTENT(IN)     :: MXMAT,MXCOMP,MXALLCV,NMAT,NCOMP,nflud,
     &                          mxcomp_suf,ncomp_suf
      INTEGER,INTENT(IN)     :: MAT_NO    (0:MXMAT)
      INTEGER,INTENT(IN)     :: MAT_CVEXT (0:MXMAT)
      INTEGER,INTENT(IN)     :: MAT_DCIDX (0:MXMAT)
      INTEGER,INTENT(IN)     :: nofld     (  nflud)
      REAL*8 ,INTENT(IN)     :: CVCENT    (3,MXALLCV)
!
      integer,intent(in)     :: IRC,nnrc
      real*8 ,intent(in)     :: vreq(mxcomp+mxcomp_suf,2,nnrc)
      integer,intent(in)     :: kind_chem(nnrc)
      character(*),intent(in):: spcnam(mxcomp+mxcomp_suf)
      real*8 ,intent(in)     :: wm(mxcomp+mxcomp_suf)
      real*8 ,intent(in)     :: ugc
!
      real*8 ,intent(in)     :: moleX     (mxallcv,mxcomp)
      real*8 ,intent(in)     :: yi        (mxallcv,mxcomp)
      real*8 ,intent(in)     :: tmp       (mxallcv)
      real*8 ,intent(inout)  :: dumya     (mxallcv)
      real*8 ,intent(inout)  :: dumyb     (mxallcv)
      real*8 ,intent(inout)  :: wdot      (mxallcv,mxcomp)
!
! --- [LOCAL ENTITIES]
!
      integer,save :: iflag=0
      real*8  :: dum1,dum2,dum3,dum4,dum5,dum6,dum7
      INTEGER :: ICVL,ICV,ICOM
      INTEGER :: IIMAT,IMAT,IMAT_U,ICVS,ICVE
      REAL*8 ,PARAMETER :: SML=1.D-20
!
! --- user array 
!
      integer,save :: ndim=2
!                                                   changed ('05.05.17)
      integer :: ierr=0,idum,i,j,ios=0,ifl=30,ifl2=31,IER   
!
      real*8  :: dmin
      integer :: inear,jnear
!                                                   changed ('05.05.17)
      real*8,save,allocatable  :: u_rate(:,:)
!
      integer,parameter ::  kmax=26,rmax=81
!                                                   changed ('05.05.17)
      integer,parameter ::  icom_CL2=2,
     &                      icom_E=7,
     &                      icom_CL=3,
     &                      icom_N2=5,
     &                      icom_N=4,
     &                      left=1,right=2
      real*8,save       ::  rate(rmax,kmax,2)
!
      real*8,dimension(kmax) :: zzz
      real*8,dimension(rmax) :: rrr
!
!---------------------------------------------------------------------
!
      dumya(:)=0.d0
      dumyb(:)=0.d0
!
      if(iflag==0) then
!                                                   changed ('05.05.17)
	ALLOCATE (u_rate(MXALLCV,2),stat=ierr)
!
        u_rate=0.d0
        if(ierr/=0) then
          call FFRABORT(1,'ALLOCATE error at u_rate(:,:)')
        endif
!-------------------------------------------------------------------
! --- (IRC=9) E + CL2 --> E + 2CL   and  (IRC=16) E + N2 --> E + 2N
!-------------------------------------------------------------------
        open(ifl,file='cldist.dat',form='formatted',status='old',
     &       action='read',iostat=ios)
        read(ifl,*) idum,idum   !=> kmax,rmax
        read(ifl,*) idum,idum
        read(ifl,*) (zzz(j),j=kmax,1,-1)
        read(ifl,*) (rrr(i),i=1,rmax)
!                                                   changed ('05.05.17)
        open(ifl2,file='ndist.dat',form='formatted',status='old',
     &       action='read',iostat=ios)
        read(ifl2,*) idum,idum   !=> kmax,rmax
        read(ifl2,*) idum,idum
        read(ifl2,*) (zzz(j),j=kmax,1,-1)
        read(ifl2,*) (rrr(i),i=1,rmax)
!
          do j=kmax,1,-1
        do i=1,rmax
!                                                   changed ('05.05.17)
            read(ifl,*) rate(i,j,1)
            read(ifl2,*) rate(i,j,2)
!
          enddo
        enddo
        close(ifl)
!                                                   changed ('05.05.17)
        close(ifl2)
!                                                   changed ('05.05.17)
        rrr=0.01d0*dble(rrr)
        zzz=0.01d0*dble(zzz)
        do j=1,kmax
          zzz(j) = zzz(j) + 0.17
        enddo
!
! --- serach near cood.
!
        do IIMAT=1,NMAT

        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        DUM1=dsqrt(CVCENT(1,ICVL)**2+CVCENT(2,ICVL)**2)
        DUM3=CVCENT(3,ICVL)

        dmin=1.d0
        if(DUM3>0.17d0.and.DUM1<0.16d0) then
          do i=1,rmax
          do j=1,kmax
          DUM4=rrr(i)
          DUM5=zzz(j)
          DUM6=(DUM1-DUM4)**2+(DUM3-DUM5)**2
          if(DUM6<dmin) then
            dmin=DUM6
            inear=i
            jnear=j
          endif
          enddo
          enddo
!                                                   changed ('05.05.17)
          u_rate(ICVL,1)=rate(inear,jnear,1)  ! =>attention
          u_rate(ICVL,2)=rate(inear,jnear,2)
!
        endif
        enddo
        enddo
        iflag=1
      endif
!
      if(IRC==9) then     !E + CL2 --> E + 2CL
        do 200 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          IMAT_U=nofld(IMAT)     ! IMAT_U: material number in fflow.ctl
          if(IMAT<0) cycle       ! IMAT < 0 : solid material
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          dum2=vreq(ICOM_CL2,right,IRC)-vreq(ICOM_CL2,left,IRC)
          dum3=vreq(ICOM_CL,right,IRC)-vreq(ICOM_CL,left,IRC)
          do 300 ICVL=ICVS,ICVE
!                                                   changed ('05.05.17)
            dum1=u_rate(ICVL,1)*moleX(ICVL,ICOM_CL2)
     &           **vreq(ICOM_CL2,left,IRC)
!
!-----------------
! -----  CL2  ----
!-----------------
            WDOT(ICVL,ICOM_CL2)=WDOT(ICVL,ICOM_CL2)
     &                        +wm(ICOM_CL2)*dum2*dum1
!-----------------
! -----  CL  -----
!-----------------
            WDOT(ICVL,ICOM_CL)=WDOT(ICVL,ICOM_CL)
     &                       +wm(ICOM_CL)*dum3*dum1
!-------------------
! -----  E  --------
!-------------------
            WDOT(ICVL,ICOM_E)=0.d0
 300      enddo
 200    enddo
      endif
!
!                                                   added ('05.05.17)
      if(IRC==16) then     !E + N2 --> E + 2N
        do 600 IIMAT=1,NMAT
          IMAT=MAT_NO(IIMAT)
          IMAT_U=nofld(IMAT)     ! IMAT_U: material number in fflow.ctl
          if(IMAT<0) cycle       ! IMAT < 0 : solid material
          ICVS=MAT_CVEXT(IIMAT-1)+1
          ICVE=MAT_CVEXT(IIMAT)
          dum2=vreq(ICOM_N2,right,IRC)-vreq(ICOM_N2,left,IRC)
          dum3=vreq(ICOM_N,right,IRC)-vreq(ICOM_N,left,IRC)
          do 500 ICVL=ICVS,ICVE
!                                                   changed ('05.05.17)
            dum1=u_rate(ICVL,2)*moleX(ICVL,ICOM_N2)
     &           **vreq(ICOM_N2,left,IRC)
!
!-----------------
! -----  N2   ----
!-----------------
            WDOT(ICVL,ICOM_N2)=WDOT(ICVL,ICOM_N2)
     &                          +wm(ICOM_N2)*dum2*dum1
!-----------------
! -----  N   -----
!-----------------
            WDOT(ICVL,ICOM_N)=WDOT(ICVL,ICOM_N)
     &                         +wm(ICOM_CL)*dum3*dum1
!-------------------
! -----  E  --------
!-------------------
            WDOT(ICVL,ICOM_E)=0.d0
 500      enddo
 600    enddo
      endif
      return
      end subroutine USER_CHEM_RATE
!
