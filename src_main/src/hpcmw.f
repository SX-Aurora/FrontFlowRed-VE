!
!      subroutine solve_poisson_hpc
!      subroutine ICCG12_HPC
!      subroutine ICCG12_HPC_VECT
!      subroutine mSORT
!      SUBROUTINE CMICCG
!      subroutine SMC
!      SUBROUTINE VICCG
!      SUBROUTINE AXMLT
!      SUBROUTINE DECOMP
!      SUBROUTINE LUSUBB
!      SUBROUTINE SEND_RECIVE
!      
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        subroutine  solve_poisson_1D(deltt,
     &  LVEDGE,LBC_SSF,LCYCSF,SFAREA,CVCENT,CVVOLM,wiface,
     &  MAT_CV,MAT_INDEX,MAT_CVEXT,MAT_DCIDX,MAT_NO,mat_cal,MAT_CFIDX,
     &  ipfix,kdbf,diag,rhs,rho,
     &  iterp,repsp,aepsp,ierror)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     NCV  : Total number in local CPU domain 
!            (=>NP in HPCMW Gr.)
!            (NCV=NCVIN when one cpu)
!     NCVIN: Inner CV number in local CPU domain
!            (NCVIN<NCV when multi-CPU)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension 
      use module_constant
      use module_hpcutil
      use module_io       ,only : ifll,ifle
      use module_cgsolver ,only : aepscg,repscg,itercg
      use module_material ,only : lclsd
      use module_metrix,only    : msk
      use module_metrix,only    : aalw
      use module_metrix,only    : aaup  =>W2K1
      use module_metrix,only    : BB
      use module_metrix,only    : INL
      use module_metrix,only    : IAL
      use module_metrix,only    : INU
      use module_metrix,only    : IAU   =>IW2K2
      use module_metrix,only    : indexL=>iwork1
      use module_metrix,only    : indexU=>iwork2
      use module_metrix,only    : aacol=>aax,NCOL1,NCOL2

      use module_metrix,only    : itemL =>iwork3
      use module_metrix,only    : itemU =>iwork4
      use module_metrix,only    : AL    =>W1K7
      use module_metrix,only    : AU    =>W1K8
      use module_metrix,only    : D     =>W1K1
      use module_vector,only    : ICVS_V,ICVE_V,
     &                            ICFS_V,ICFE_V,
     &                            ICVSIN_V,ICVEIN_V,
     &                            IDCS_V,IDCE_V,N2
!
      implicit none
! --- [dummy qrgument]
!
      real*8 ,intent(in)    :: deltt
      integer,intent(in)    :: ipfix (    MXMAT) 
      integer,intent(in)    :: kdbf  (  MXCVFAC)
      integer,intent(in)    :: LVEDGE(2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF( MXSSFBC)
      integer,intent(in)    :: LCYCSF(  MXSSFBC)
      INTEGER,INTENT(IN)    :: MAT_CV(     MXCV)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_DCIDX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)
      integer,intent(in)    :: MAT_CFIDX(0:MXMAT)
      real*8 ,intent(in)    :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)    :: wiface(  MXCVFAC)
      real*8 ,intent(in)    :: CVCENT(3,MXALLCV)
      real*8 ,intent(inout) :: diag  (  MXALLCV)
      real*8 ,intent(inout) :: rhs   (  MXALLCV)
      real*8 ,intent(in)    :: CVVOLM(  MXALLCV)
      real*8 ,intent(inout) :: rho    (  MXALLCV)
      integer,intent(out)   :: ierror,iterp
      real*8 ,intent(out)   :: repsp,aepsp

!
! --- [local entities]
!
      real*8  :: dx,dy,dz,alf,aepsq,repsq,alfn,dl,ER,pref,rmax
      real*8  :: wifa,wifb,dva,dvb,coef,avegra,avegrb
      integer :: IOPT,nl,nu,iterq,ierr=0
      integer :: i,j,k,l,m,n,kk
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IMODE,ICOM,IMD,ICM,IFLAG,IL,NSOLVR
!
      integer               :: NPL,NPU
      integer               :: ERROR
!
! --- +-------------------+
! --- | PARAMETER CONTROL |
! --- +-------------------+
!
      allocate(aaup(MXCV,MAXIE),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating aaup(:,:) error-1 in solve_poisson_1D'
        call FFRABORT(1,'solve_poisson_1D')
      endif
      allocate (IAU(MXCV,MAXIE),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating IAU(:,:) error-2 in solve_poisson_1D'
        call FFRABORT(1,'solve_poisson_1D')
      endif

!
      ierror=0
!
!
! --- < 1. Make up coefficient matrix >-
!
! --- < 1.1 all over the domain >--
!
!      call bc_prdmsk(MAT_NO,MAT_CVEXT,MAT_DCIDX,mat_cal,
!     &               LVEDGE,LBC_SSF,LCYCSF,msk)
!
! --- 
!
      D(:)=0.d0
      INL=0
      INU=0
      aalw=0.d0
      aaup=0.d0
      IAL=0
      IAU=0
      BB(1:MXCV)=rhs(1:MXCV)
!
      do ICFL=ICFS_V,ICFE_V
      if(kdbf(ICFL).lt.2) then
        ICVA=LVEDGE(1,ICFL)
        ICVB=LVEDGE(2,ICFL)
! --- regular-face/periodic and outlet BC
          alf=SFAREA(4,ICFL)
     &    *deltt*dsqrt((CVCENT(1,ICVB)-CVCENT(1,ICVA))**2
     &                +(CVCENT(2,ICVB)-CVCENT(2,ICVA))**2
     &                +(CVCENT(3,ICVB)-CVCENT(3,ICVA))**2)
     &           /(((CVCENT(1,ICVB)-CVCENT(1,ICVA))*SFAREA(1,ICFL)
     &             +(CVCENT(2,ICVB)-CVCENT(2,ICVA))*SFAREA(2,ICFL)
     &             +(CVCENT(3,ICVB)-CVCENT(3,ICVA))*SFAREA(3,ICFL))**2)
!----------------------------------------
        if(kdbf(ICFL).ne.1) then
! --- regular & periodic face
          ICVA=msk(ICVA)
          ICVB=msk(ICVB)
          D(ICVA)=D(ICVA)+alf
          D(ICVB)=D(ICVB)+alf
!
          nl=min(ICVA,ICVB)
          nu=max(ICVA,ICVB)
!
          INU(nl)=INU(nl)+1
          IAU(nl,INU(nl))=nu
!
          INL(nu)=INL(nu)+1
          IAL(nu,INL(nu))=nl
!
          aaup(nl,INU(nl))=-alf
          aalw(nu,INL(nu))=-alf
        elseif( kdbf(ICFL).eq.1 ) then
! --- outlet BC:kdbf(ICFL)=1; Inlet:kdbf(ICFL)=3 !
          D(ICVA)=D(ICVA)+alf
        endif
      endif
      enddo
!
! --- < 1.2 Clear solid part >--
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,D)
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
! --- < 1.3 Pressure fix point >-
!
!      do 130 IIMAT=1,NMAT
!      if(lclsd(IIMAT).eq.1) then
!        n=ipfix(IIMAT)
!        BB(n)=0.d0
!        D(n)=1.d0
!
!        do 131 i=1,INU(n)
!          m=IAU(n,i)
!          do 132 j=1,INL(m)
!          if( IAL(m,j).eq.n ) aalw(m,j)=0.d0
!  132     continue
!  131   continue
!
!        do 133 i=1,INL(n)
!          m=IAL(n,i)
!          do 134 j=1,INU(m)
!          if( IAU(m,j).eq.n ) aaup(m,j)=0.d0
! 134      continue
! 133    continue
!
!        do j=1,INL(n)
!        aalw(n,j)=0.d0
!        enddo
!        do j=1,INU(n)
!        aaup(n,j)=0.d0
!        enddo
!        INL(n)=0
!        INU(n)=0
!      endif
!  130 continue
!
! --- < (1) CG of compressed array version, BY Zhang >
!
! --- Reordering
!
!
! --- CG solver: (NSOLVR=1:By Zhang ;NSOLVR=2:By IBM ;NSOLVR=3: By FFLOW
!
      do ICVL=ICVS_V,ICVE_V
        NL=INL(ICVL)
        NCOL1=0
        NCOL2=0
        aacol=0.d0
        do 145 k=1,NL
        NCOL1(k)=IAL(ICVL,k)
        aacol(k)=aalw(ICVL,k)
 145    enddo
        call mSORT(NCOL1,NCOL2,NL)
        do 150 k=NL,1,-1
        IAL(ICVL,NL-k+1)=NCOL1(NCOL2(k))
        aalw(ICVL,NL-k+1)=aacol(NCOL2(k))
 150    enddo
!
        NU=INU(ICVL)
        NCOL1=0
        NCOL2=0
        aacol=0.d0
        do 155 k=1,NU
        NCOL1(k)=IAU(ICVL,k)
        aacol(k)=aaup(ICVL,k)
 155    enddo
        call mSORT(NCOL1,NCOL2,NU)
        do 157 k=NU,1,-1
        IAU(ICVL,NU-k+1)= NCOL1(NCOL2(k))
        aaup(ICVL,NU-k+1)=aacol(NCOL2(k))
 157    enddo
      enddo
!
      allocate (indexL(0:MXCV), indexU(0:MXCV),stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error-3 in solve_poisson_1D'
        call FFRABORT(1,'solve_poisson_1D')
      endif
!
      indexL=0
      indexU=0
!
      do ICVL=ICVS_V,ICVE_V
        indexL(ICVL)=indexL(ICVL-1)+INL(ICVL)
        indexU(ICVL)=indexU(ICVL-1)+INU(ICVL)
      enddo
!
        
      NPL=indexL(NCV)
      NPU=indexU(NCV)
! --- 
      allocate (itemL(NPL),itemU(NPU),
     &          AL(NPL),AU(NPU),
     &          stat=ierr)
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error-4 in solve_poisson_1D'
        call FFRABORT(1,'solve_poisson_1D')
      endif
      AU(:)=0.d0
! --- 
      do ICVL=ICVS_V,ICVE_V
        do 175 k=1,INL(ICVL)
          kk=k+indexL(ICVL-1)
          itemL(kk)=IAL(ICVL,k)
          AL(kk)=  aalw(ICVL,k)
 175    continue
        do 178 k=1,INU(ICVL)
          kk=k+indexU(ICVL-1)
          itemU(kk)=IAU(ICVL,k)
          AU(kk)=  aaup(ICVL,k)
 178    continue
      enddo
!
! --- +------------------+
! --- | ITERATIVE solver |
! --- +------------------+
!

      aepsq=aepscg
      repsq=repscg
      iterq=itercg
      ERROR=0
!
      deallocate(aaup,IAU)
!
      NSOLVR=1
      if(NSOLVR==1) then
        CALL ICCG12_HPC_VECT
     &    (MXCV,NCVIN,NCV,NPL,NPU,MXMAT,NMAT,
     &     MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &     D,AL,indexL,itemL,AU,indexU,itemU,BB,
     &     repsq,aepsq,iterq,ERROR)
      else
      endif


      aepsp=aepsq
      repsp=repsq
      iterp=iterq
!
      if(ERROR.ne.0) then
          write(ifle,*) 
     &  '*** ERR:  Pressure Poisson iter. (iterq) too large: ',
     &   iterq,repsq,aepsq
          write(ifle,*) 'At my_rank= ',my_rank
          call FFRABORT(1,'solve_poisson_1D')
      endif
!
      do ICVL=ICVS_V,ICVE_V
        rhs(ICVL)=BB(ICVL)
      enddo
!
      deallocate(itemL,itemU,indexL,indexU,AL,AU)
!
      return
!
!
 9999 continue
      write(ifle,*) '(solve_poisson_1D)'
      ierror=1
      return
      end subroutine solve_poisson_1D
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine mSORT (STEM,INUM,NN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      integer,intent(in)    :: NN
      integer,intent(in)    :: STEM(NN)
      integer,intent(inout) :: INUM(NN)
      integer :: ii,jj,ITEM
!
      do ii = 1,NN
        INUM(ii)= ii
      enddo

      do ii= 1,NN-1
!CDIR NOVECTOR
      do jj= 1,NN-ii
        if (STEM(INUM(jj)) .lt. STEM(INUM(jj+1))) then
          ITEM      = INUM(jj+1)
          INUM(jj+1)= INUM(jj)
          INUM(jj)  = ITEM
        endif
      enddo
      enddo
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      return
      end

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine ICCG12_HPC_VECT
     &(MXNP,NCVIN,NP,NPL,NPU,MXMAT,NMAT,
     & MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     & D,AL,indexL,itemL,AU,indexU,itemU,B,
     & reps,aeps,itr,ERROR)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     This subroutine is ICCG solver for serial CPU and HPC
!     Unsymm ICCG solver                                
!     N : inner vertex number
!     NP: Total vertex number
!     
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_hpcutil
      use module_cgsolver,only : nostop
      use module_io,    only : ifll,ifle
      use module_metrix,only    : P=>W1K2
      use module_metrix,only    : Q=>W1K3
      use module_metrix,only    : R=>W1K4
      use module_metrix,only    : Z=>W1K5
      use module_metrix,only    : X=>W1K6
!      use module_material,only  : KMAT_S,MAT_ALL,MAT_S
      use module_vector,only    : ICVS_V,ICVE_V,
     &                            ICFS_V,ICFE_V,
     &                            ICVSIN_V,ICVEIN_V,
     &                            IDCS_V,IDCE_V
!
      implicit none
!
! --- [dummy arguments]
!
      real*8 ,intent(inout)  ::  reps,aeps
      integer,intent(inout)  ::  itr,NCVIN,NP,NPL,NPU,MXNP
      integer,intent(inout)  ::  ERROR
      integer,intent(in)     ::  MXMAT,NMAT
      integer,intent(in)     ::  MAT_INDEX(0:MXMAT)
      integer,intent(in)     ::  MAT_CVEXT(0:MXMAT)
      integer,intent(in)     ::  MAT_NO(   0:MXMAT)
      logical,INTENT(IN)     ::  mat_cal(  0:MXMAT)
!
      real*8 ,intent(inout)  ::  D(0:MXNP)
      real*8 ,intent(inout)  ::  B(MXNP)
      real*8 ,intent(inout)  ::  AU(NPU)
      real*8 ,intent(inout)  ::  AL(NPL)
!
      integer,intent(in)     ::  indexU(0:MXNP)
      integer,intent(in)     ::  itemU(    NPU)
      integer,intent(in)     ::  indexL(0:MXNP)
      integer,intent(in)     ::  itemL(    NPL)
!
! --- [local entities]
!
      integer  :: LsU,LeU,LsL,LeL,L,Lnod,maxitr=300
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0,ONE=1.d0
      integer  :: i,j,k,it1,it2,jer,itrS,NCV
      real*8   :: BETA,X1,X2,Y,ALPHA,repsS,
     &            S,C1,C2,C3,rmax,rmax0
      integer  :: IMAT,IIMAT,ICVS,ICVE,ICVL

!
! --- allocating array:
! --- INIT. 
!
      rmax=1.d0
!CIDR NODEP
      do 50 ICVL=1,NCVIN
        rmax=min(rmax,D(ICVL))
 50   continue
!
! --- MPI: RMAX
!     
      IF(NPE.gt.1) THEN
        CALL hpcrmin(RMAX)
      ENDIF
!
      if(RMAX.le.ZERO.and.my_rank.eq.root) then
        write(ifle,*) '### interface error -1- (ICCG12_HPC_VECT)'
        call FFRABORT(1,'ICCG12_HPC_VECT')
      endif
!
      ERROR=0
      NCV=NP
      jer=0
      it1=1
      it2=itr
      itr=0
      rmax0=ONE
      maxitr=it2
!
!CIDR NODEP
      do ICVL=1,NCV
        P(ICVL)=ZERO
        Q(ICVL)=ZERO
        Z(ICVL)=ONE
        R(ICVL)=ZERO
        X(ICVL)=ZERO
      enddo 
      RMAX=ZERO
!CIDR NODEP
      DO ICVL=1,NCV
        RMAX=max(RMAX,ABS(B(ICVL)))
      ENDDO
      RMAX0=RMAX
!
      if(NPE.gt.1) then
        CALL hpcrsum(RMAX)
      endif
      if(RMAX<=max(aeps,reps*rmax0)) goto 999
      rmax0=rmax
!
! --- < 1. diagonal SCALING >
!CIDR NODEP
      DO ICVL=1,NCV
        D(ICVL)=1.d0/(Dsqrt(D(ICVL)))
        Z(ICVL)=D(ICVL)
      ENDDO
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,D)
      ENDIF
!
      do ICVL=1,NCV
        LsU=indexU(ICVL-1)+1
        LeU=indexU(ICVL)
!CIDR NODEP
        DO L=LsU,LeU
          Lnod=itemU(L)
          AU(L)=AU(L)*D(ICVL)*D(Lnod)
        ENDDO
!
        LsL=indexL(ICVL-1)+1
        LeL=indexL(ICVL)
!CIDR NODEP
        DO L=LsL,LeL
          Lnod=itemL(L)
          AL(L)=AL(L)*D(ICVL)*D(Lnod)
        ENDDO
      enddo
!
      rmax=0.d0
!CIDR NODEP
      DO ICVL=1,NCV   !ICVS_V,ICVE_V   !ICV=1,NP
        B(ICVL)=B(ICVL)*D(ICVL)
        D(ICVL)=1.D0
        rmax=rmax+abs(B(ICVL))
      ENDDO
!
      if(NPE.gt.1) then
        CALL hpcrsum(RMAX)
      endif
!
      if(RMAX<=max(aeps,reps*rmax0)) goto 999
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,D)
      ENDIF
!---------------------------------------------
! --- < 3. incomplete Cholesky DECOMPOSION >
!---------------------------------------------
!
      DO 10 ICVL=1,NCV    !ICVS_V,ICVE_V    !ICV=1,NP
      if(d(ICVL).le.0.d0) then
        write(ifle,*) '### interface error -2- (ICCG12_HPC_VECT)'
        write(ifle,*) 'ICVL=, D(ICVL)= ',ICVL,d(ICVL)
        call FFRABORT(1,'ICCG12_HPC_VECT')
      endif
      D(ICVL)=1.d0/(SML+D(ICVL))
      LsL=indexU(ICVL-1)+1
      LeL=indexU(ICVL)
      DO L=LsL,LeL
      Lnod=itemU(L)
      d(Lnod)=d(Lnod)-d(ICVL)*AU(L)*AU(L)
      ENDDO
 10   CONTINUE
!---------------------------------------------
! --- < 4. Initial set for iteration >
!---------------------------------------------
  888 continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,X)
      ENDIF
!---------------
! --- Q(:)=A*X
!---------------
      DO 30 ICVL=1,NCVIN    !ICVSIN_V,ICVEIN_V   !ICV=1,N
      S=D(ICVL)*X(ICVL)
      LsL=indexL(ICVL-1)+1
      LeL=indexL(ICVL)
!CIDR NODEP
      DO L=LsL,LeL
        Lnod=itemL(L)
        S=S+AL(L)*X(Lnod)
      ENDDO
!
      LsU=indexU(ICVL-1)+1
      LeU=indexU(ICVL)
!CIDR NODEP
      DO L=LsU,LeU
        Lnod=itemU(L)
        S=S+AU(L)*X(Lnod)
      ENDDO
      Q(ICVL)=S
 30   CONTINUE
!---------------
! --- R(:)
!---------------
      rmax=ZERO
!CIDR NODEP
      DO 100 ICVL=1,NCV !ICVSIN_V,ICVEIN_V   !ICV=1,N
        R(ICVL)=B(ICVL)-Q(ICVL)
        rmax=rmax+abs(R(ICVL))
 100  CONTINUE
!------------------------------
! --> MPI_SUM => rmax
!------------------------------
      if(NPE.gt.1) then
        CALL hpcrsum(RMAX)
      endif
      if(RMAX<=max(aeps,reps*rmax0)) goto 999
!---------------
! --- P(:)=LDLT
!---------------
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,P)
      ENDIF
!
      DO 110 ICVL=1,NCV    !ICVSIN_V,ICVEIN_V   !ICV=1,N
      S=R(ICVL)
      LsL=indexL(ICVL-1)+1
      LeL=indexL(ICVL)
      DO L=LsL,LeL
        Lnod=itemL(L)
        S=S-AL(L)*P(Lnod)
      ENDDO
      P(ICVL)=D(ICVL)*S
 110  CONTINUE
!
      DO 130 ICVL=NCV,1,-1   !ICVEIN_V,ICVSIN_V,-1   !ICV=N,1,-1
      S=0.D0
      LsU=indexU(ICVL-1)+1
      LeU=indexU(ICVL)
      DO L=LsU,LeU
        Lnod=itemU(L)
        S=S+AU(L)*P(Lnod)
      ENDDO
      P(ICVL)=P(ICVL)-D(ICVL)*S
 130  CONTINUE
!
      C1=0.D0
!CIDR NODEP
      DO 200 ICVL=1,NCV   !ICVS_V,ICVE_V   !ICV=1,N
        C1=C1+R(ICVL)*P(ICVL)
 200  CONTINUE
!------------------------------
! --> MPI_SUM => C1
!------------------------------
      if(NPE.gt.1) then
        CALL hpcrsum(C1)
      endif
!
      if(ABS(C1).lt.SML) jer=1
!
      if(jer.gt.0) then
!CIDR NODEP
        do 330 ICVL=1,NCV   !ICVSIN_V,ICVEIN_V   !ICV=1,N
          x(ICVL)=x(ICVL)+p(ICVL)
  330   continue
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
! --- ||||||||||||||||||||||| ITERATION START |||||||||||||||||||||||||
!
      DO 300 itr=it1,it2
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXNP,NP,P)
      ENDIF
!---------------
! --- Q=A*P
!---------------
      DO 310 ICVL=1,NCV   !ICVSIN_V,ICVEIN_V   !ICV=1,N
      S=D(ICVL)*P(ICVL)
!
      LsL=indexL(ICVL-1)+1
      LeL=indexL(ICVL)
!CIDR NODEP
      DO L=LsL,LeL
        Lnod=itemL(L)
        S=S+AL(L)*P(Lnod)
      ENDDO
!
      LsU=indexU(ICVL-1)+1
      LeU=indexU(ICVL)
!CIDR NODEP
      DO L=LsU,LeU
        Lnod=itemU(L)
        S=S+AU(L)*P(Lnod)
      ENDDO
      Q(ICVL)=S
 310  CONTINUE
!------------------
! --> MPI => Q(:)
!------------------
      C2=0.D0
!CIDR NODEP
      DO 400 ICVL=1,NCVIN   !ICVSIN_V,ICVEIN_V   !ICV=1,N
        C2=C2+P(ICVL)*Q(ICVL)
 400  CONTINUE
!------------------
! --> MPI_SUM => C2
!------------------
      if(NPE.gt.1) then
        CALL hpcrsum(C2)
      endif
!
      IF(C2.eq.zero) THEN
        jer=1
        it1=itr
        goto 888
      ENDIF
!
!
      ALPHA=C1/(SML+C2)
!
      X1=0.D0
      X2=0.D0
!CIDR NODEP
      DO 500 ICVL=1,NCVIN    !ICVSIN_V,ICVEIN_V    !ICV=1,N
        Y=X(ICVL)
        X(ICVL)=X(ICVL)+ALPHA*P(ICVL)
        R(ICVL)=R(ICVL)-ALPHA*Q(ICVL)
        X1=X1+Y*Y
        X2=X2+(X(ICVL)-Y)*(X(ICVL)-Y)
 500  CONTINUE
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,R)
      endif
      rmax=0.d0
!CIDR NODEP
      DO 550 ICVL=ICVSIN_V,ICVEIN_V   !ICV=1,N
        rmax=max(rmax,abs(R(ICVL)))
 550  CONTINUE
!
! --> MPI => X1,X2,rmax
!
      if(NPE.gt.1) then
        CALL hpcrmax(RMAX)
      endif
!
      if( RMAX.le.min(aeps,reps*rmax0)) then
         ERROR=0
         goto 999
      endif
!------------------
! --- LDLT
!------------------
!
! --> MPI => Q(:)
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,Q)
      ENDIF
!
      DO 510 ICVL=ICVSIN_V,ICVEIN_V   !ICV=1,N
      S=R(ICVL)
      LsL=indexL(ICVL-1)+1
      LeL=indexL(ICVL)
!CIDR NODEP
      DO L=LsL,LeL
        Lnod=itemL(L)
        S=S-AL(L)*Q(Lnod)
      ENDDO
      Q(ICVL)=D(ICVL)*S
 510  CONTINUE
!------------------
! --> MPI => Q(:)
!------------------
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,Q)
      ENDIF
!
      DO 530 ICVL=NCVIN,1,-1   !ICVEIN_V,ICVSIN_V,-1   !ICV=N,1,-1
      S=0.D0
      LsU=indexU(ICVL-1)+1
      LeU=indexU(ICVL)
!CIDR NODEP
      DO L=LsU,LeU
        Lnod=itemU(L)
        S=S+AU(L)*Q(Lnod)
      ENDDO
      Q(ICVL)=Q(ICVL)-D(ICVL)*S
 530  CONTINUE
!
      C3=ZERO
!CIDR NODEP
      DO 540 ICVL=ICVSIN_V,ICVEIN_V   !ICV=1,N
        C3=C3+R(ICVL)*Q(ICVL)
 540  CONTINUE
!---------------------
! --> MPI_SUM => C3
!---------------------
      if(NPE.gt.1) then
        CALL hpcrsum(C3)
      endif
!
      if(ABS(C3).lt.SML) then
        jer=1
        it1=itr
        goto 888
      endif
!
      BETA=C3/(SML+C1)
      C1=C3
!
!CIDR NODEP
      DO 770 ICVL=ICVSIN_V,ICVEIN_V   !ICV=1,N
        P(ICVL)=Q(ICVL)+BETA*P(ICVL)
 770  CONTINUE
!
 300  CONTINUE
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &   ' ### WRN: ICCG Solver is NOT converged: ',
     & 'iccg_Iter= ', maxitr,
     &  ' NOT stop and go into next time step'
         endif
         goto 999
      endif
!
      ERROR=1
!
 999  CONTINUE
!
      aeps=max(aeps,rmax)
      reps=max(reps,rmax/(rmax0))
!
      do ICVL=ICVSIN_V,ICVEIN_V
        B(ICVL)=X(ICVL)*Z(ICVL)
      enddo
!------------------
! --> MPI => B(:)
!------------------
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXNP,NP,B)
      ENDIF
!
      if(NPE.gt.1) then
        CALL hpcimax(itr)
        call hpcrmax(reps)
      endif
!
      RETURN
!
      END subroutine ICCG12_HPC_VECT
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CMICCG(calsld,LL,IL,NL,LU,IU,NU,NCVIN,
     &           MAT_INDEX,MAT_CVEXT,MAT_NO,ical_sld,NMAT,MXMAT,
     &           MXALLCV,MXCV,N,N2,INO,NO,IT,IJT,IW,IW2K_VECT)
!***********************************************************************
!*  CUTHILL-MCKEE METHOD FOR VECTOR ICCG.                              *
!*  COPYRIGHT:  YASUNORI USHIRO AND TSUTOMU OGUNI  SEP. 1 1991  VER. 1 *
!***********************************************************************
      use module_hpcutil
      use module_vector  ,ONLY : NUmax,NLmax
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)       :: MXALLCV,MXCV,NL,NU,NCVIN,N
      integer,intent(in)       :: ical_sld,NMAT,MXMAT
      integer,intent(inout)    :: IL(MXCV)
      integer,intent(inout)    :: LL(MXALLCV,NL)
      integer,intent(inout)    :: IU(MXCV)
      integer,intent(inout)    :: LU(MXCV,NU)
      integer,intent(inout)    :: IT(MXCV)
      integer,intent(inout)    :: IJT(0:MXCV)
      integer,intent(inout)    :: INO(0:MXCV)
      integer,intent(inout)    :: NO
      integer,intent(in)       :: N2
      INTEGER,INTENT(IN)       :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)       :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)       :: MAT_NO   (0:MXMAT)
      integer,intent(inout)    :: IW2K_VECT(MXCV,2)
      logical,INTENT(IN)       :: calsld
      integer,intent(inout)    :: IW(MXCV)
!      real*8 ,intent(inout)    :: BB(MXCV)
!      real*8 ,intent(inout)    :: D(0:MXCV)
!
! --- [local entities]
!
      integer :: KCT,KC,KCT0,KMIN,KMAX,INEW,IOLD,NN,NMX
      integer :: I,J,K,JC,JN,II,N1,N3,JJ,NCV,INEW1
      integer :: IIMAT,ICVS,ICVE,ICVL,NCVINN
      real*8  :: RMAXS,rmaxs0 ! Added by Y.Takahashi, 2021/11/13
      !
!      if(calsld.and.(ical_sld==1.or.
!     &   ical_sld==2.or.ical_sld==4)) then
        IW2K_VECT(:,:)=0
!        do IIMAT=1,NMAT
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        do ICVL=ICVS,ICVE
!        IW2K_VECT(ICVL,1)=2
!        enddo
!        enddo
!      endif
!-----------------------------------
! --- 
!-----------------------------------
      NCV    = N
      INO(0) = 0
      INO(1) = 0
      NO     = 1
      KCT    = 0
      IW(:)=0
      NCVINN=NCVIN
      if(calsld.and.(ical_sld==1.or.ical_sld==2.or.ical_sld==4)) then
        NCVINN=NCV
      endif
!---------------------
!-- RCM ordering
!---------------------
!
      IW(:)=0
      IT(:)=0
      DO 10 I=1,NCVINN
!      IF(IL(I)==0.and.IW2K_VECT(I,1)==1) THEN
      IF(IL(I).EQ.0) THEN
        INO(NO)=INO(NO)+1
        KCT=KCT+1
        IT(KCT)=I
        IW(I)=-1
        IW2K_VECT(I,2)=1
      ELSE
        IW(I)=0
      ENDIF
   10 ENDDO
!
      KC   = KCT
      KCT0 = 0
!
 100  CONTINUE
!
      KMIN=NCVINN
      KMAX=1
      DO J=1,KC
        JC=IT(KCT0+J)
        JN=IU(JC)
        DO I=1,JN
          II=LU(JC,I)
!          IF(II==0.or.IW2K_VECT(I,1)==0) cycle  !6666
          IW(II)=IW(II)+1
          KMIN=MIN0(II,KMIN)
          KMAX=MAX0(II,KMAX)
        ENDDO
      ENDDO
!
      NO=NO+1
      KCT0=KCT
!
      DO 60 I=KMIN,MIN(KMAX,NCVINN)
!        IF(IW(I)==IL(I).and.IL(I)/=0.and.IW2K_VECT(I,1)==1) THEN
        IF(IW(I)==IL(I).and.IL(I)/=0) THEN
          KCT=KCT+1
          IW(I)=-1
          IT(KCT)=I
          IW2K_VECT(I,2)=1
        ENDIF
 60   ENDDO
      KC=KCT-KCT0
      INO(NO)=KCT
      IF (KC.NE.0.AND.NO.LE.NCVINN) GO TO 100
!
!-----------
!--- add
!-----------
      IF(calsld.and.(ical_sld==1.or.
     &       ical_sld==2.or.ical_sld==4)) then
!CIDR NODEP
        KCT0=KCT
        do I=1,NCV
        if(IW2K_VECT(I,2)/=1) then
!        if(IW2K_VECT(I,1)/=1) then
          KCT=KCT+1
          IT(KCT)=I
        endif
        enddo
!
        KC=KCT-KCT0
        IF(KC.NE.0) then
          NO=NO+1
          INO(NO)=KCT
        endif
!
        DO INEW=1,NCV
        IOLD=IT(INEW)
        IJT(IOLD)=INEW
        ENDDO
!
      else
!CIDR NODEP
        DO I=NCVIN+1,N
        IT(I)=I
        IJT(I)=I
        enddo
!CIDR NODEP
        DO INEW=1,NCVIN
        IOLD=IT(INEW)
        IJT(IOLD)=INEW
        ENDDO
      endif
!
! --- ---------------------------
!
      NMX=N
      N1 = NMX / N2
      N3 = NMX - N1 * N2
      JJ=0
      DO 200 J=1,N1
        DO 190 I=1,N2
        JJ = JJ + 1
        IW(JJ)=NMX+I
  190   ENDDO
  200 ENDDO
      DO 210 I=1,N3
        JJ=JJ+1
        IW(JJ)=NMX+I
  210 ENDDO
!
! --- ---------------------------
!
      DO J=1,NL
!CDIR NODEP
        DO I=1,N
        IF (LL(I,J).EQ.0) THEN
          LL(I,J)=IW(I)
        ELSE
          LL(I,J)=IJT(LL(I,J))
        ENDIF
        enddo
      enddo
!
      DO J=1,NU
!CDIR NODEP
        DO 260 I=1,N
        IF (LU(I,J).EQ.0) THEN
         LU(I,J)=IW(I)
        ELSE
         LU(I,J)=IJT(LU(I,J))
        ENDIF
  260   CONTINUE
      enddo
!
! --- ---------------------------
!
      DO 110 J=1,NL
!CDIR NODEP
        DO 90 I=1,N
          IW(I)=LL(IT(I),J)
 90     enddo
!CDIR NODEP
        DO I=1,N
          LL(I,J) = IW(I)
        enddo
 110  enddo
!
      DO 140 J=1,NU
!CDIR NODEP
        DO 120 I=1,N
          IW(I) = LU(IT(I),J)
 120    enddo
!CDIR NODEP
        DO 130 I=1,N
          LU(I,J) = IW(I)
 130    enddo
 140  enddo
!
! --- ---------------------------
!
!CDIR NODEP
      DO 150 I=1,N
        IW(I) = IL(IT(I))
 150  enddo
!CDIR NODEP
      DO 160 I=1,N
        IL(I) = IW(I)
 160  enddo
!CDIR NODEP
      DO I=1,N
        IW(I) = IU(IT(I))
      enddo
	  RMAXS  = 0.     ! add debag  OSHIMA 2021.10.25 by NuFD
	  rmaxs0 = 0.     ! add debag  OSHIMA 2021.10.25 by NuFD
!CDIR NODEP
      DO I=1,N
        IU(I) = IW(I)
      enddo
!6666 ------------------------------------------
!      DO I=1,N
!        IW(I) = IW2K_VECT(IT(I),1)
!      enddo
!      DO I=1,N
!        IW2K_VECT(I,1)=IW(I)
!      enddo
!      DO I=1,N
!        IW(I) = IW2K_VECT(IT(I),2)
!      enddo
!      DO I=1,N
!        IW2K_VECT(I,2)=IW(I)
!      enddo
!
      IW(:)=0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      IW(ICVL)=1
      enddo
      enddo
!--------------------------------------------
      NLmax= -N
      NUmax= -N
      do i= 1, N
        NLmax= max(NLmax, IL(i))
        NUmax= max(NUmax, IU(i))
      enddo
      RETURN
!
      END SUBROUTINE CMICCG
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE VICCG(IMODE1,NCOLOR,calsld,IVDIM,
     &           MXALLCV,MXCV,N,NCVIN,NU,NL,N2,NO,
     &                 MAT_INDEX,MAT_CVEXT,MAT_NO,ical_sld,NMAT,MXMAT,
     &                 AL,LL,IL,
     &                 AU,LU,IU,
     &                 B,X,LN,D,Z,
     +                 IT,IJT,W,IW2K_VECT,TEMP,IW,
     &                 ITR,EPS,aeps,IER
     &   )
***********************************************************************
*  ICCG method on a vector computer.                                  *
*  Copyright:  Yasunori Ushiro and Tsutomu Oguni  SEP. 1 1992  VER. 2 *
***********************************************************************
!
      use module_hpcutil
      use module_cgsolver,only : nostop
      use module_io,only       : ifll,ifle
      use module_vector  ,ONLY : NUmax,NLmax
!
      implicit none
!      include 'mpif.h'
!
! --- [ dummy arguments ]
!
      integer,intent(in)       :: MXALLCV,MXCV,N,NCVIN,NU,NL,N2,NO,
     &                            NCOLOR,IMODE1,IVDIM
      integer,intent(in)       :: ical_sld,NMAT,MXMAT
      integer,intent(inout)    :: IER,ITR
      real*8 ,intent(inout)    :: EPS,aeps
      real*8 ,intent(inout)    :: AL(MXALLCV,0:NL),AU(MXCV,0:NU)
      integer,intent(in)       :: LN(0:MXCV)
      integer,intent(in)       :: LL(MXALLCV,NL),LU(MXCV,NU)
      integer,intent(in)       :: IL(MXCV),IU(MXCV)
!
      integer,intent(in)       :: IJT(0:MXCV)
      integer,intent(in)       :: IT(MXCV)
!
      real*8 ,intent(inout)    :: W(MXCV+N2,IVDIM)
      real*8 ,intent(inout)    :: D(0:MXCV)
      real*8 ,intent(inout)    :: B(MXCV)
      real*8 ,intent(inout)    :: Z(MXCV)
      real*8 ,intent(inout)    :: X(MXCV+N2)
      real*8 ,intent(inout)    :: TEMP(MXCV)
      INTEGER,INTENT(IN)       :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)       :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)       :: MAT_NO (  0:MXMAT)
      integer,intent(in)       :: IW2K_VECT(MXCV,2)
      logical,INTENT(IN)       :: calsld
      integer,intent(inout)       :: IW(MXCV)
!
! --- [local entities]
!
      integer :: maxitr=300
      integer :: I,IS,IE,J,L,K,IIMAT,ICVS,ICVE,ICVL
      real*8  :: R1,R2,R3,RNORM,BETA,BNORM,ALP,ERR
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0
      integer,parameter :: ID=1,IR=2,IP=3,IQ=4
      real*8  :: RMAXS,rmaxs0
      integer :: ierr,NCVINN
      real*8  :: dum1
!
      maxitr=ITR
      NCVINN=NCVIN
      IF(calsld.and.(ical_sld==1.or.ical_sld==2.or.ical_sld==4)) then
        NCVINN=N   !CVIN
      endif
!
!W(I,1)=>
!W(I,2)=>R(ICVL)
!W(I,3)=>P(ICVL)
!W(I,4)=>Q(ICVL)
!
      W=0.d0
!--------------------
C  Change B,X,D 
!--------------------
!CDIR NODEP
      if(NCOLOR<0) then 
        DO I=1,N
        W(I,1)=B(IT(I))
        W(I,2)=X(IT(I))
        W(I,3)=D(IT(I))
        W(I,4)=dble(IW(IT(I)))
        enddo
!CDIR NODEP
        DO I=1,N
        B(I)=W(I,1)
        X(I)=W(I,2)
        D(I)=W(I,3)
        IW(I)=INT(W(I,4))  !IW(IT(I))
        enddo
      else
        call FFRABORT(1,'ERR: gfgdg')
        DO I=1,N
        W(I,1)=B(IT(I))
        W(I,4)=dble(IW(IT(I)))
        enddo
!CDIR NODEP
        DO I=1,N
        B(I)=W(I,1)
        IW(I)=INT(W(I,4))  !IW(IT(I))
        enddo
      endif
!-------------------
C  Change AL,AU
!-------------------
!
      if(NCOLOR<0) then 
        DO 10 J=1,NLmax
!CDIR NODEP
        DO I=1,N
          W(I,1)=AL(IT(I),J)
        enddo
!CDIR NODEP
        DO I=1,N
          AL(I,J)=W(I,1)
        enddo
   10   enddo
!
        DO 16 J=1,NUmax
!CDIR NODEP
        DO I=1,N
          W(I,1)=AU(IT(I),J)
        enddo
!CDIR NODEP
        DO I=1,N
          AU(I,J)=W(I,1)
        enddo
   16   enddo
      endif
!
!-------------------------------------
!      X(N+1:)=0.d0
!-------------------------------------
!
!CDIR NODEP
      W(1:N,ID)=D(1:N)
!CDIR NODEP
!!!        D(:)=0.d0    !7777
      W(N+1:MXCV+N2,ID)=0.0D0
!CDIR NODEP
      DO I=1,N
      W(I,ID)=1.d0/dsqrt(SML+W(I,ID))
      Z(I)=W(I,ID)
      ENDDO
!
      IF(IMODE1==1) THEN
        DO L=1,NUmax
!CDIR NODEP
        do I=1,N
        AU(I,L)=AU(I,L)*W(I,ID)*W(LU(I,L),ID)
        enddo
        ENDDO
        DO L=1,NLmax
!CDIR NODEP
        do I=1,N
        AL(I,L)=AL(I,L)*W(I,ID)*W(LL(I,L),ID)
        enddo
        ENDDO
!
!CDIR NODEP
        DO I=1,N
        B(I)=B(I)*W(I,ID)
        W(I,ID)=1.D0
        ENDDO
      ELSE
        DO L=1,NUmax
!CDIR NODEP
        do I=1,N
        AU(I,L)=AU(I,L)*W(I,ID)*W(LU(I,L),ID)
     &         *sign(1.d0,D(I))
        enddo
        ENDDO
        DO L=1,NLmax
!CDIR NODEP
        do I=1,N
        AL(I,L)=AL(I,L)*W(I,ID)*W(LL(I,L),ID)
     &         *sign(1.d0,D(I))
        enddo
        ENDDO
!
!CDIR NODEP
        DO I=1,N
        B(I)=B(I)*W(I,ID)*sign(1.d0,D(I))
        W(I,ID)=1.D0
        ENDDO
      ENDIF
!
!-------------------------------------
!CDIR NODEP  !1
      DO I=1,MXCV+N2
        X(I)   =0.0D0
        W(I,IR)=0.0D0
        W(I,IP)=0.0D0
        W(I,IQ)=0.0D0
      enddo
!
! --- 
!
!--------------------------------------------
!  Incomlete Cholesky decomposition 
!--------------------------------------------
       CALL DECOMP_IC(NCVIN,N2,B,NUmax,NLmax,ical_sld,
     &     calsld,
     &     AL,LL,NL,IL,
     &     AU,LU,NU,IU,
     &     IT,IJT,IW2K_VECT,TEMP,IW,
     &     MXALLCV,MXCV,N,LN,NO,W(:,ID),EPS,IER)

!       CALLDECOMP_LU(NCVIN,N2,B,NUmax,NLmax,ical_sld,
!     &           calsld,
!     &           AL,LL,NL,IL,
!     &           AU,LU,NU,IU,
!     &           IT,IJT,IW2K_VECT,
!     &           MXALLCV,MXCV,N,LN,NO,W(:,ID),EPS,IER)
!----------
C  R=A*X
!----------
      CALL AXMLT(1,IMODE1,
     &           MXALLCV,NCVIN,N2,NUmax,NLmax,calsld,ical_sld,
     &           AL,LL,NL,AU,LU,NU,MXCV,N,TEMP,IW,
     &           IT,IJT,IW2K_VECT,X,W(1,IR))
!--------------------------------7777
      BNORM = 0.0D0
      DO I=1,NCVINN
      BNORM=BNORM+B(I)**2*dble(IW(I))
      enddo
      IF(NPE.gt.1) CALL hpcrsum(BNORM)
!--------------------------------
!      rmaxs=0.d0
!      DO I=1,NCVINN
!      rmaxs=(max(rmaxs,abs(B(I))*dble(IW(I))))
!      ENDDO
!      IF(NPE.gt.1) call hpcrmax(RMAXS)
!      rmaxs0=rmaxs
!--------------------------------
!CDIR NODEP
      DO I=1,N
        W(I,IR)=(B(I)-W(I,IR))!*dble(IW2K_VECT(I,1))
      enddo
!
!-------------------------------------
C  P=(LU)**(-1)*R 
!-------------------------------------
      CALL LUSUBB(N2,W(1,ID),AL,LL,NL,AU,LU,NU,NCVIN,calsld,
     &            ical_sld,NUmax,NLmax,IW2K_VECT,TEMP,IW,
     &            MXALLCV,MXCV,N,LN,NO,IT,IJT,W(1,IR),W(1,IP))
!
      R1=0.0D0

      DO I=1,NCVINN
      R1=R1+W(I,IR)*W(I,IP)*dble(IW(I))
      enddo
      
      IF(NPE.gt.1) CALL hpcrsum(R1)

!------------------------------------------------------------------
C  Iterations      
!------------------------------------------------------------------
      DO 200 L=1,ITR
!-------------
C  Q= A*P
!-------------
      CALL AXMLT(2,IMODE1,
     &           MXALLCV,NCVIN,N2,NUmax,NLmax,calsld,ical_sld,
     &           AL,LL,NL,AU,LU,NU,MXCV,N,TEMP,IW,
     &           IT,IJT,IW2K_VECT,W(1,IP),W(1,IQ))
      R2=0.0D0

      DO I=1,NCVINN
      R2=R2+W(I,IP)*W(I,IQ)*dble(IW(I))
      enddo
!
      IF(NPE.gt.1) CALL hpcrsum(R2)
!
      IF(DABS(R2).LT.SML) THEN
        WRITE(*,'(1x,a,I8,F14.4)') 'ZERO IN R2 ', L,DABS(R2)
        IER=2
        GOTO 900
      elseif(DABS(R2)>1.d10) then
        call FFRABORT(1,'ERR: VICCG: R2 > 1.E10')
      ENDIF
!
      ALP=R1/R2
!
!CDIR NODEP
      DO I=1,N
        X(I)   =(X(I)   +ALP*W(I,IP))!*dble(IW2K_VECT(I,1))
        W(I,IR)=(W(I,IR)-ALP*W(I,IQ))!*dble(IW2K_VECT(I,1))
      enddo
!-------------------------------------------7777
      RNORM=0.0D0
      DO I=1,NCVINN
      RNORM=RNORM+W(I,IR)**2*dble(IW(I))
      enddo
      IF(NPE.gt.1) CALL hpcrsum(RNORM)
      ERR=DSQRT(RNORM/BNORM)
!-------------------------------------------
!      rmaxs=0.d0
!      DO I=1,NCVINN
!      rmaxs=max(rmaxs,abs(W(I,IR))*dble(IW(I)))
!      ENDDO
!      IF(NPE.gt.1) CALL hpcrmax(RMAXS)
!      ERR=RMAXS
!-------------------------------------------
!
!-----------------
C  Q=(LU)**(-1)*R 
!-----------------
!
      CALL LUSUBB(N2,W(1,ID),AL,LL,NL,AU,LU,NU,NCVIN,calsld,
     &            ical_sld,NUmax,NLmax,IW2K_VECT,TEMP,IW,
     &            MXALLCV,MXCV,N,LN,NO,IT,IJT,W(1,IR),W(1,IQ))
!
!aeps EPS
      IF(ERR.LT.max(aeps,EPS*rmaxs0)) THEN !max(aeps,reps*rmaxs0)
        IER=0
        GOTO 900
      ELSE
        R3=0.0D0
        DO I=1,NCVINN
        R3=R3+W(I,IR)*W(I,IQ)*dble(IW(I))
        enddo
!
        if(NPE.gt.1) CALL hpcrsum(R3)
!
        IF(DABS(R1).LT.SML) THEN
          WRITE(ifle,*) 'Zero in R1 ', L,DABS(R1),my_rank
          IER = 3
          GO TO 900
        ENDIF
!
        BETA=R3/R1
        R1=R3
!CDIR NODEP
        DO I=1,N  !import
        W(I,IP)=(W(I,IQ)+BETA*W(I,IP))!*dble(IW2K_VECT(I,1))
        enddo
      ENDIF
!
  200 CONTINUE
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &    ' ### WRN: VICCG Solver is NOT converged: ',
     & 'VICCG_Iter= ', maxitr,
     & ' NOT stop and go into next time step'
         endif
         IER=0
         goto 900
      endif
!
      IER=1
!
 900  CONTINUE
!
!CDIR NODEP
      DO I=1,N
      W(I,3)=X(IJT(I))*Z(IJT(I))
      enddo
!CDIR NODEP
      DO I=1,N
      X(I)=W(I,3)
      enddo
!
      ITR = L
      EPS = ERR
!
!      if(NPE.gt.1) then
!        CALL SOLVER_SEND_RECV(1,MXCV,N,X(1:))
!      endif
!
      RETURN
!
      END SUBROUTINE VICCG
!


CCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE AXMLT(ICALL,IMODE1,
     &           MXALLCV,NCVIN,N2,NUmax,NLmax,calsld,ical_sld,
     &                 AL,LL,NL,AU,LU,NU,MXCV,N,TEMP,IW,
     &                 IT,IJT,IW2K_VECT,P,Q)
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!----------------------------------
! --- Y=A*X  
!----------------------------------
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MXALLCV,MXCV,N,NU,NL,N2,NCVIN,
     &                         NUmax,NLmax,ical_sld,IMODE1,
     &                         ICALL
      integer,intent(in)    :: LL(MXALLCV,NL),LU(MXCV,NU)
!      real*8 ,intent(inout) :: D(0:MXCV)
      real*8 ,intent(in)    :: AL(MXALLCV,0:NL)
      real*8 ,intent(in)    :: AU(MXCV,0:NU)
      real*8 ,intent(inout) :: P(MXCV+N2),Q(MXCV+N2)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      real*8 ,intent(inout) :: TEMP(MXCV)
      logical,INTENT(IN)    :: calsld
      integer,intent(inout)       :: IW(MXCV)
!------------------------
! --- [local entities]
!------------------------
      integer :: I,J,IS,IE,K,II  
!--------------------
C  Q= D*P + AL*P
!--------------------
!
!CDIR NODEP
      DO I=1,N
        Q(I)=(P(I)+AL(I,1)*P(LL(I,1)))*dble(IW(I))
      enddo
!
      DO 80 J=2,NLmax 
!CDIR NODEP
        DO I=1,N
          Q(I)=(Q(I)+AL(I,J)*P(LL(I,J)))*dble(IW(I))
        enddo
   80 enddo
!
!--------------------
C  Q=Q + AU*P 
!--------------------
!
      DO 100 J=1,NUmax
!CDIR NODEP
        DO I=1,N
        Q(I)=(Q(I)+AU(I,J)*P(LU(I,J)))*dble(IW(I))
        enddo
  100 enddo
!
!-----------------------------------------------------------
      RETURN
!-----------------------------------------------------------
      END SUBROUTINE AXMLT
C
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE DECOMP_IC(NCVIN,N2,B,NUmax,NLmax,ical_sld,
     &           calsld,
     &           AL,LL,NL,IL,
     &           AU,LU,NU,IU,
     &           IT,IJT,IW2K_VECT,TEMP,IW,
     &           MXALLCV,MXCV,N,LN,NO,DD,EPS,IER)
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
! Incomlete decomposition
!
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: IER,NCVIN,NUmax,NLmax
      integer,intent(in)    :: ical_sld
      real*8 ,intent(in)    :: EPS
      integer,intent(in)    :: MXALLCV,MXCV,N,NL,N2,NO
      real*8 ,intent(in)    :: AL(MXALLCV,0:NL)
      integer,intent(in)    :: LN(0:MXCV)
      integer,intent(in)    :: LL(MXALLCV,NL)

      integer,intent(in)    :: NU
      integer,intent(in)    :: LU(MXCV,NU)
      real*8 ,intent(in)    :: AU(MXCV,0:NU)
      integer,intent(in)    :: IL(MXCV)
      integer,intent(in)    :: IU(MXCV)

      real*8 ,intent(inout) :: DD(MXCV+N2)
!      real*8 ,intent(inout) :: D(0:MXCV)
      real*8 ,intent(inout) :: B(MXCV)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      real*8 ,intent(inout) :: TEMP(MXCV)
      logical,INTENT(IN)    :: calsld
      integer,intent(inout)       :: IW(MXCV)
!
! --- [local entities]
!
      integer :: I,IS,IE,J,K,Lnod,L
      real*8  :: EPS1,DW,dum1
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      EPS1=EPS*1.0D-3
      IER=0
!------------------------------------------------------
      DD(N+1:MXCV+N2)=0.d0

      DO 100 J=1,NO
        IS=LN(J-1)+1
        IE=LN(J)
        IF((IE-IS+1)>0) then
!
!----------------------------------------
!
!CDIR NODEP
          DO I=IS,IE
          DD(I)=1.D0-AL(I,1)**2*DD(LL(I,1))!*dble(IW(I))
          enddo
!
          DO 60 K=2,NLmax-1
!CDIR NODEP
          DO I=IS,IE
          DD(I)=DD(I)-AL(I,K)**2*DD(LL(I,K))!*dble(IW(I))
          enddo
   60     enddo
!
!CDIR NODEP
          DO I=IS,IE
          DW=DD(I)-AL(I,NLmax)**2*DD(LL(I,NLmax))!*dble(IW(I))
          IF (DABS(DW).LT.EPS1) THEN
            IER = 1
            WRITE(*,*) '(SUBR. VICCG) Singular at step = ', I
            RETURN
          ENDIF
          DD(I)=1.0D0/DW!*dble(IW2K_VECT(I,1))
          enddo
        elseif(.false.) then !7777
!
!----------------------------------------
!
!CDIR NOVECTOR
          DO I=IS,IE
          DD(I)=1.D0-AL(I,1)**2*DD(LL(I,1))
          enddo
          DO K=2,NLmax-1
!CDIR NOVECTOR
          DO I=IS,IE
            DD(I)=DD(I)-AL(I,K)**2*DD(LL(I,K))
          enddo
          enddo
!CDIR NOVECTOR
          DO I=IS,IE
          DW=DD(I)-AL(I,NLmax)**2*DD(LL(I,NLmax))
          IF (DABS(DW).LT.EPS1) THEN
            IER = 1
            WRITE(*,*) '(SUBR. VICCG) Singular at step = ', I
            RETURN
          ENDIF
          DD(I)=1.0D0/DW!*IW(I)
          enddo
        endif
  100 enddo

      END SUBROUTINE DECOMP_IC
!
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE DECOMP_LU(NCVIN,N2,B,NUmax,NLmax,ical_sld,
     &           calsld,
     &           AL,LL,NL,IL,
     &           AU,LU,NU,IU,
     &           IT,IJT,IW2K_VECT,
     &           MXALLCV,MXCV,N,LN,NO,DD,EPS,IER)
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
! Incomlete decomposition
!
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: IER,NCVIN,NUmax,NLmax
      integer,intent(in)    :: ical_sld
      real*8 ,intent(in)    :: EPS
      integer,intent(in)    :: MXALLCV,MXCV,N,NL,N2,NO
      real*8 ,intent(in)    :: AL(MXALLCV,0:NL)
      integer,intent(in)    :: LN(0:MXCV)
      integer,intent(in)    :: LL(MXALLCV,NL)

      integer,intent(in)    :: NU
      integer,intent(in)    :: LU(MXCV,NU)
      real*8 ,intent(in)    :: AU(MXCV,0:NU)
      integer,intent(in)    :: IL(MXCV)
      integer,intent(in)    :: IU(MXCV)

      real*8 ,intent(inout) :: DD(MXCV+N2)
!      real*8 ,intent(inout) :: D(0:MXCV)
      real*8 ,intent(inout) :: B(MXCV)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
!      real*8 ,intent(inout) :: TEMP(MXCV)
      logical,INTENT(IN)    :: calsld
!
! --- [local entities]
!
      integer :: I,IS,IE,J,K,Lnod,L
      real*8  :: EPS1,DW,dum1
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      EPS1=EPS*1.0D-3
      IER=0
!------------------------------------------------------
      DD(N+1:MXCV+N2)=0.d0
!      return

!      
!      DO I=1,N
!      DD(I)=1.d0/(DD(I))
!      do 231 j=1,IU(I)   !IQ(ICVL,1)+1,IQ(ICVL,2)
!      k=LU(I,j)
!      dum1=0.d0
!      do 232 l=1,IL(k)  !l=1,IQ(k,1)
!      IF(LL(k,l)==I) dum1=DD(I)*AL(k,l)
! 232  enddo
 !     DD(k)=DD(k)-dum1*AU(I,j)
! 231  enddo
!      ENDDO
!      return
!
! --- ---------------------------
!

      DO 100 J=1,NO
        IS=LN(J-1)+1
        IE=LN(J)
        IF((IE-IS+1)>0) then
!
!----------------------------------------
!
!CDIR NODEP
          DO I=IS,IE
!          DD(I)=1.D0-AL(I,1)**2*DD(LL(I,1))
          DD(I)=1.D0-AL(I,1)*AU(I,1)
     &         *SQRT(DD(LL(I,1))*DD(LU(I,1)))
          enddo
!
          DO 60 K=2,NLmax-1
!CDIR NODEP
          DO I=IS,IE
!          DD(I)=DD(I)-AL(I,K)**2*DD(LL(I,K))
          DD(I)=DD(I)-AL(I,K)*AU(I,K)
     &         *SQRT(DD(LL(I,K))*DD(LU(I,K)))
          enddo

   60     enddo
!
!CDIR NODEP
!
          DO I=IS,IE
!          DW=DD(I)-AL(I,NLmax)**2*DD(LL(I,NLmax))
          DW=DD(I)-AL(I,NLmax)*AU(I,NLmax)*
     &       SQRT(DD(LL(I,NLmax))*DD(LU(I,NLmax)))
          IF (DABS(DW).LT.EPS1) THEN
            IER = 1
            WRITE(*,*) '(SUBR. VICCG) Singular at step = ', I
            RETURN
          ENDIF
          DD(I)=1.0D0/DW
          enddo
        endif
  100 enddo
!
      END SUBROUTINE DECOMP_LU
!
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE LUSUBB(N2,DD,AL,LL,NL,AU,LU,NU,NCVIN,calsld,
     &              ical_sld,NUmax,NLmax,IW2K_VECT,TEMP,IW,
     &              MXALLCV,MXCV,N,LN,NO,IT,IJT,R,Q)
CCCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MXALLCV,MXCV,N2,N,NL,NU,
     &                         NO,NCVIN,
     &                         NUmax,NLmax
      real*8 ,intent(inOUT) :: AL(MXALLCV,0:NL)
      real*8 ,intent(in)    :: AU(MXCV,0:NU)
      real*8 ,intent(in)    :: DD(MXCV+N2)
      real*8 ,intent(inout) :: Q(MXCV+N2)
      real*8 ,intent(inout) :: R(MXCV+N2)
      integer,intent(in)    :: LL(MXALLCV,NL),LU(MXCV,NU)
      integer,intent(in)    :: LN(0:MXCV)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
!      real*8 ,intent(inout) :: D(0:MXCV)
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      real*8 ,intent(inout) :: TEMP(MXCV)
      integer,intent(in)    :: ical_sld
      logical,intent(in)    :: calsld
      integer,intent(in)    :: IW(MXCV)
!
! --- [local entities]
!
      integer :: I,IS,IE,J,K,ISUM
      real*8  :: SUM
!
!-----------------
C  Q=(L)**(-1)*R
!-----------------
      DO 100 J=1,NO
        IS=LN(J-1)+1
        IE=LN(J)
        IF((IE-IS+1)>0) then
!CDIR NODEP
          DO 20 I=IS,IE
          Q(I)=(R(I)-AL(I,1)*Q(LL(I,1)))*dble(IW(I))
 20       enddo
!
          DO 60 K=2,NLmax-1
!CDIR NODEP
          DO I=IS,IE
          Q(I)=(Q(I)-AL(I,K)*Q(LL(I,K)))*dble(IW(I))
          enddo
   60     enddo
!
!CDIR NODEP
          DO I=IS,IE
          Q(I)=(DD(I)*(Q(I)-AL(I,NLmax)*Q(LL(I,NLmax))))
     &         *dble(IW(I))
          enddo
        elseif(.false.) then
!CDIR NOVECTOR
          DO I=IS,IE
          Q(I)=(R(I)-AL(I,1)*Q(LL(I,1)))!*dble(IW2K_VECT(I,1))!*dble(IW(I))
          enddo
!
          DO K=2,NLmax-1
!CDIR NOVECTOR
          DO I=IS,IE
          Q(I)=(Q(I)-AL(I,K)*Q(LL(I,K)))!*dble(IW2K_VECT(I,1))!*dble(IW(I))
          enddo
          enddo
!
!CDIR NOVECTOR
          DO I=IS,IE
          Q(I)=(DD(I)*(Q(I)-AL(I,NLmax)*Q(LL(I,NLmax))))
     &        !*dble(IW2K_VECT(I,1))!*dble(IW(I))
          enddo
        endif
  100 enddo
!-----------------
C  Q= (U)**(-1)*Q 
!-----------------
      CALL SEND_RECIVE(calsld,ical_sld,
     & IW2K_VECT,TEMP,MXCV,N,1,MXCV+N2,IT,IJT,Q)
!------------------------------------------------------
      DO 200 J=NO,1,-1         !NO-1,1,-1
        IS=LN(J-1)+1
        IE=LN(J)
        DO 140 K=1,NUmax
        IF((IE-IS+1)>0) then
!CDIR NODEP
          DO 120 I=IS,IE
          Q(I)=(Q(I)-DD(I)*AU(I,K)*Q(LU(I,K)))*dble(IW(I))
 120      enddo
        elseif(.false.) then
!CDIR NOVECTOR
          DO I=IS,IE
          Q(I)=(Q(I)-DD(I)*AU(I,K)*Q(LU(I,K)))*dble(IW(I))
          enddo
        endif
 140    enddo
  200 enddo
!------------------------------------------------------
      CALL SEND_RECIVE(calsld,ical_sld,
     & IW2K_VECT,TEMP,MXCV,N,1,MXCV+N2,IT,IJT,Q)
!------------------------------------------------------
      RETURN

      END SUBROUTINE LUSUBB
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SEND_RECIVE(calsld,ical_sld,
     & IW2K_VECT,TEMP,MXCV,N,INI,MXCVN2,IT,IJT,WORK)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: ical_sld,MXCV,N,INI,MXCVN2
      logical,intent(in)    :: calsld
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      real*8 ,intent(inOUT) :: TEMP(MXCV)
      real*8 ,intent(inOUT) :: WORK(INI:MXCVN2)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
!      integer,intent(in)    :: IW(MXCV)
!
! --- [local entities]
!
      integer :: I
!
!
!
      if(NPE.gt.1) then
        DO I=1,N
          TEMP(I)=WORK(IJT(I))
      ENDDO
        call SOLVER_SEND_RECV0(1,1,MXCV,N,TEMP)
        DO I=1,N
          WORK(I)=TEMP(IT(I))
        enddo
      endif
!
      return
!
      end SUBROUTINE SEND_RECIVE
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE VICCG_B(m,cn,nq,ndiag,IMODE1,NCOLOR,calsld,IVDIM,
     &           MXALLCV,MXCV,N,NCVIN,NU,NL,N2,NO,
     &                 MAT_INDEX,MAT_CVEXT,MAT_NO,ical_sld,NMAT,MXMAT,
     &                 AL,LL,IL,
     &                 AU,LU,IU,
     &                 B,X,LN,D,Z,
     +                 IT,IJT,W,IW2K_VECT,TEMP,IW,
     &                 ITR,EPS,IER
     &   )
***********************************************************************
*  ICCG method on a vector computer.                                  *
*  Copyright:  Yasunori Ushiro and Tsutomu Oguni  SEP. 1 1992  VER. 2 *
***********************************************************************
!
      use module_hpcutil
      use module_cgsolver,only : nostop
      use module_io,only       : ifll,ifle
      use module_vector  ,ONLY : NUmax,NLmax
      use module_model,   only : monitor_stp
      use module_nowtime, only : iter,time
!
      implicit none
!
! --- [ dummy arguments ]
!
      character,intent(in)     :: cn
      integer,intent(in)       :: MXALLCV,MXCV,N,NCVIN,NU,NL,N2,NO,
     &                            NCOLOR,IMODE1,IVDIM,m
      integer,intent(in)       :: nq,ndiag
      integer,intent(in)       :: ical_sld,NMAT,MXMAT
      integer,intent(inout)    :: IER,ITR
      real*8 ,intent(inout)    :: EPS
      real*8 ,intent(inout)    :: AL(MXALLCV,0:NL),AU(MXCV,0:NU)
      integer,intent(in)       :: LN(0:MXCV)
      integer,intent(inout)    :: LL(MXALLCV,NL),LU(MXCV,NU)
      integer,intent(in)       :: IL(MXCV),IU(MXCV)
!
      integer,intent(in)       :: IJT(0:MXCV)
      integer,intent(in)       :: IT(MXCV)
!
      real*8 ,intent(inout)    :: W(MXCV+N2,IVDIM)
      real*8 ,intent(inout)    :: D(0:MXCV)
      real*8 ,intent(inout)    :: B(MXCV)
      real*8 ,intent(inout)    :: Z(MXCV)
      real*8 ,intent(inout)    :: X(MXCV+N2)
      real*8 ,intent(inout)    :: TEMP(MXCV)
      INTEGER,INTENT(IN)       :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)       :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)       :: MAT_NO (  0:MXMAT)
      integer,intent(in)       :: IW2K_VECT(MXCV,2)
      logical,INTENT(IN)       :: calsld
      integer,intent(inout)    :: IW(MXCV)            
!
! --- [local entities]
!
      integer :: maxitr=300
      integer :: I,IS,IE,J,L,K,IIMAT,ICVS,ICVE,ICVL
      real*8  :: RR1,RR2,RR3,RNORM,BNORM,ALP,ERR
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0
      integer,parameter :: P0=1,P1=2,Q0=3,Q1=4,R0=5,R1=6,R2=7,ID=8
      real*8  :: RMAXS,rmaxs0,s,c1,ttt,sss,uuu,vvv,alpha,beta,omega
      integer :: ierr,NCVINN
      real*8  :: dum1,dum2
      logical,save :: multi_M
!----------------------------
! --- 
!----------------------------
      maxitr=ITR
!
      multi_M=calsld.and.
     &   (ical_sld==1.or.ical_sld==2.or.ical_sld==4)
      NCVINN=NCVIN
      IF(multi_M) then
        NCVINN=N        !CVIN
      endif
!----------------------------
      ttt=0.d0
      sss=0.d0
      uuu=0.d0
      vvv=0.d0
      alpha=0.d0
      beta=0.d0
      omega=0.d0
!
      W=0.d0
      X=0.d0
!--------------------
C  Change B,X,D 
!--------------------
!CDIR NODEP
      if(m==1.or.ndiag>1) then
        DO I=1,N
        W(I,1)=B(IT(I))
        W(I,2)=X(IT(I))
        W(I,3)=D(IT(I))
        W(I,4)=dble(IW(IT(I)))
        enddo
!CDIR NODEP
        DO I=1,N
        B(I)=W(I,1)
        X(I)=W(I,2)
        D(I)=W(I,3)
        IW(I)=INT(W(I,4))  !IW(IT(I))
        enddo
      else
        DO I=1,N
        W(I,1)=B(IT(I))
!!!!!!!!!!!!        W(I,4)=dble(IW(IT(I)))
        enddo
!CDIR NODEP
        DO I=1,N
        B(I)=W(I,1)
!!!!!!!!!!!!        IW(I)=INT(W(I,4))  !IW(IT(I))
        enddo
      endif
!-------------------
C  Change AL,AU
!-------------------
      if(m==1) then
        DO 10 J=1,NLmax
!CDIR NODEP
        DO I=1,N
        W(I,1)=AL(IT(I),J)
        enddo
!CDIR NODEP
        DO I=1,N
        AL(I,J)=W(I,1)
        enddo
   10   enddo
!
        DO 16 J=1,NUmax
!CDIR NODEP
        DO I=1,N
        W(I,1)=AU(IT(I),J)
        enddo
!CDIR NODEP
        DO I=1,N
        AU(I,J)=W(I,1)
        enddo
   16   enddo
      endif
!-------------------------------------
      X(N+1:)=0.d0 
!-------------------------------------
!CDIR NODEP
      W(1:N,ID)=D(1:N)
!CDIR NODEP
      W(N+1:MXCV+N2,ID)=0.D0
!CDIR NODEP
      DO I=1,N
      W(I,ID)=1.d0/dsqrt(W(I,ID))
      Z(I)=W(I,ID)
      ENDDO
!--------------------------
! --- diagonal scaling 
!--------------------------
      if(m==1) then
        DO L=1,NUmax
!CDIR NODEP
        do I=1,N
        AU(I,L)=AU(I,L)*W(I,ID)*W(LU(I,L),ID)*sign(1.d0,D(I))
        enddo
        ENDDO
!
        DO L=1,NLmax
!CDIR NODEP
        do I=1,N
        AL(I,L)=AL(I,L)*W(I,ID)*W(LL(I,L),ID)*sign(1.d0,D(I))
        enddo
        ENDDO
      endif
!
!CDIR NODEP
      DO I=1,N
      B(I)=B(I)*W(I,ID)*sign(1.d0,D(I))         !D(I)
      W(I,ID)=1.D0
      ENDDO
!
!-------------------------------------
!CDIR NODEP  !1
!      DO I=1,MXCV+N2
      X(:)    =0.0D0
      W(:,1:7)=0.0D0
!      enddo
!
! --- 
!
!--------------------------------------------
!  incomplete LU decomposition
!--------------------------------------------
      CALL DECOMP_LU_B(NCVIN,N2,NUmax,NLmax,ical_sld,
     &     multi_M,
     &     AL,LL,NL,IL,
     &     AU,LU,NU,IU,
     &     IT,IJT,IW2K_VECT,TEMP,
     &     MXALLCV,MXCV,N,LN,NO,W(1,ID),EPS,IER)
!      
!----------
C  R=A*X 
!----------
!
      CALL AXMLT_B(1,IMODE1,
     &           MXALLCV,NCVIN,N2,NUmax,NLmax,multi_M,ical_sld,
     &           AL,LL,NL,AU,LU,NU,MXCV,N,IW,
     &           IT,IJT,IW2K_VECT,X,W(1,R0))  !R0
!
!
!CDIR NODEP
      DO I=1,NCVINN
      W(I,R0)=(B(I)-W(I,R0))*dble(IW(I))!5555
      ENDDO
!
      sss=0.d0
      if(multi_M) then
        DO I=1,NCVINN
        dum2=sign(1.d0,W(I,R0))
     &       *max(1.d0,abs(W(I,R0)))
        sss=sss+dum2*W(I,R0)*dble(IW(I))
        enddo
!
!        do IIMAT=1,NMAT
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        do I=ICVS,ICVE
!        dum2=sign(1.d0,W(IJT(I),R0))
!     &       *max(1.d0,abs(W(IJT(I),R0)))
!        sss=sss+dum2*W(IJT(I),R0)
!        enddo
!        enddo
!
        DO I=1,NCVINN
        W(I,P0)=W(I,R0)*dble(IW(I))!5555 
        W(I,R1)=W(I,R0)*dble(IW(I))!5555 
        dum1=sign(1.d0,W(I,R0))*max(1.d0,abs(W(I,R0)))
        W(I,R0)=dum1*dble(IW(I))!5555
        enddo
!        DO I=1,NCVINN
!        W(I,P0)=W(I,R0)
!        W(I,R1)=W(I,R0)
!        dum1=sign(1.d0,W(I,R0))*max(1.d0,abs(W(I,R0)))
!        W(I,R0)=dum1
!        ENDDO

      else
        DO I=1,NCVINN
        W(I,P0)=W(I,R0)
        W(I,R1)=W(I,R0)
        dum1=sign(1.d0,W(I,R0))*max(1.d0,abs(W(I,R0)))
        sss=sss+dum1*W(I,R0)*dble(IW(I))
        W(I,R0)=dum1*dble(IW(I))!5555
        ENDDO
      endif
!
      IF(NPE.gt.1) CALL hpcrsum(sss)
!CDIR NODEP
!-------------------------------------
C  P=(LU)**(-1)*R 
!-------------------------------------
      CALL LUSUBB_B(
     &    m,N2,W(1,ID),AL,LL,NL,AU,LU,NU,NCVIN,multi_M,
     &    ical_sld,NUmax,NLmax,IW2K_VECT,TEMP,IW,
     &    MXALLCV,MXCV,N,LN,NO,IT,IJT,W(1,R1),W(1,P1))
!                                             R1,P1
      DO I=1,NCVINN
      X(I)=(X(I)+W(I,P1))*dble(IW(I))!5555
      ENDDO
!
!------------------------------------------------------------------
C  Iterations      
!------------------------------------------------------------------
      DO 200 L=1,ITR
      CALL LUSUBB_B(
     &    m,N2,W(1,ID),AL,LL,NL,AU,LU,NU,NCVIN,multi_M,
     &    ical_sld,NUmax,NLmax,IW2K_VECT,TEMP,IW,
     &    MXALLCV,MXCV,N,LN,NO,IT,IJT,W(1,P0),W(1,P1))
!-------------
C  Q= A*P
!-------------
      CALL AXMLT_B(2,IMODE1,
     &           MXALLCV,NCVIN,N2,NUmax,NLmax,multi_M,ical_sld,
     &           AL,LL,NL,AU,LU,NU,MXCV,N,IW,
     &           IT,IJT,IW2K_VECT,W(1,P1),W(1,Q1))
!                                  P1,Q1
      vvv=0.0D0
!      IF(multi_M) then
!        do IIMAT=1,NMAT
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        do I=ICVS,ICVE
!        vvv=vvv+W(IJT(I),R0)*W(IJT(I),Q1)
!        enddo
!        enddo
!      else
        DO I=1,NCVINN
        vvv=vvv+W(I,R0)*W(I,Q1)*dble(IW(I))
        enddo
!      endif
!
      IF(NPE.gt.1) CALL hpcrsum(vvv)
!
      alpha=sss/(vvv+SML)
!
!CDIR NODEP
      DO I=1,N
        W(I,Q0)=W(I,R1)-alpha*W(I,Q1)
      enddo
!
!-----------------
C  Q=(LU)**(-1)*R
!-----------------
      CALL LUSUBB_B(
     &    m,N2,W(1,ID),AL,LL,NL,AU,LU,NU,NCVIN,multi_M,
     &    ical_sld,NUmax,NLmax,IW2K_VECT,TEMP,IW,
     &    MXALLCV,MXCV,N,LN,NO,IT,IJT,W(1,Q0),W(1,R1))
!
      CALL AXMLT_B(3,IMODE1,
     &           MXALLCV,NCVIN,N2,NUmax,NLmax,multi_M,ical_sld,
     &           AL,LL,NL,AU,LU,NU,MXCV,N,IW,
     &           IT,IJT,IW2K_VECT,W(1,R1),W(1,R2))
!                                             Q0,R1
      uuu=0.d0
      vvv=0.d0
!      if(multi_M) then
!        do IIMAT=1,NMAT
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        do I=ICVS,ICVE
!        uuu=uuu+W(IJT(I),R2)*W(IJT(I),Q0)
!        vvv=vvv+W(IJT(I),R2)**2
!        enddo
!        enddo
!      else
        DO I=1,NCVINN
        uuu=uuu+W(I,R2)*W(I,Q0)*dble(IW(I))
        vvv=vvv+W(I,R2)**2*dble(IW(I))
        enddo
!      endif
      if(NPE.gt.1) then
        CALL hpcrsum(uuu)
        CALL hpcrsum(vvv)
      endif
      omega=uuu/(vvv+SML)
!
      DO I=1,NCVINN
      X(I)=X(I)+(alpha*W(I,P1)+omega*W(I,R1))*dble(IW(I))!5555
      enddo
!
      DO I=1,NCVINN
      W(I,R1)=(W(I,Q0)-omega*W(I,R2))*dble(IW(I))!5555
      enddo
!
      rmaxs=0.d0
!      if(multi_M) then
!        do IIMAT=1,NMAT
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        do I=ICVS,ICVE
!        rmaxs=max(rmaxs,abs(W(IJT(I),R1)))
!        enddo
!        enddo
!      else
        DO I=1,NCVINN
        rmaxs=max(rmaxs,dble(IW(I))*abs(W(I,R1)))
        enddo
!      endif
!
      IF(NPE.gt.1) CALL hpcrmax(rmaxs)
!
      IF(rmaxs.LT.EPS) THEN
        IER=0
        GOTO 900
      ELSEIF(rmaxs>1.d10) then
        call FFRABORT(1,'ERR: VICCG_B: rmaxs > 1.E10')
      ENDIF
!
      ttt=0.d0
!      if(multi_M) then
!        do IIMAT=1,NMAT
!        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
!        do I=ICVS,ICVE
!        ttt=ttt+W(IJT(I),R1)*W(IJT(I),R0)
!        enddo
!        enddo
!      else
        DO I=1,NCVINN
        ttt=ttt+W(I,R1)*W(I,R0)*dble(IW(I))
        enddo
!      endif
!
      IF(NPE.gt.1) call hpcrsum(ttt)
!
      beta=ttt*alpha/(sss*omega+SML)
!
      do 470 I=1,NCVINN    !i=1,NCVIN
      W(I,P0)=(W(I,R1)+beta*(W(I,P0)-omega*W(I,Q1)))*dble(IW(I))!5555
  470 enddo
!
      sss=ttt
!
  200 CONTINUE
!
      if(nostop==1) then
!        if(mod(iter,monitor_stp)==0) then
!         if(my_rank==root) then
!           write(ifle,'(2a,I4,2a)')
!     &    ' ### WRN: VICCG_B Solver is NOT converged: ',
!     & 'Bicg_Iter= ', maxitr,
!     & ' NOT stop and go into next time step char= ',
!     & cn
!         endif
!         endif
         IER=0
         goto 900
      endif
!
      IER=1
!
 900  CONTINUE
!
!CDIR NODEP
      DO I=1,N 
      W(I,3)=X(IJT(I))*Z(IJT(I))
      enddo
!CDIR NODEP
      DO I=1,N
      X(I)=W(I,3)
      enddo
      
!
      ITR = L
      EPS = rmaxs
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,N,X(1:))
      endif
      RETURN
!
      END SUBROUTINE VICCG_B
!
CCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE AXMLT_B(ICALL,IMODE1,
     &           MXALLCV,NCVIN,N2,NUmax,NLmax,multi_M,ical_sld,
     &                 AL,LL,NL,AU,LU,NU,MXCV,N,IW,
     &                 IT,IJT,IW2K_VECT,X,P)
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!----------------------------------
! --- P=A*X  
!----------------------------------
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MXALLCV,N,NU,NL,N2,NCVIN,
     &                         NUmax,NLmax,ical_sld,IMODE1,
     &                         ICALL,MXCV
      integer,intent(in)    :: LL(MXALLCV,NL),LU(MXCV,NU)
      real*8 ,intent(in)    :: AL(MXALLCV,0:NL)
      real*8 ,intent(in)    :: AU(MXCV,0:NU)
      real*8 ,intent(inout) :: P(MXCV+N2),X(MXCV+N2)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      logical,INTENT(IN)    :: multi_M
      integer,intent(in)    :: IW(MXCV)      
!------------------------
! --- [local entities]
!------------------------
      integer :: I,J,IS,IE,K,II  
!------------------------
      DO I=1,N
      P(I)=(X(I)+AL(I,1)*X(LL(I,1)))*dble(IW(I))
      enddo
!
      DO 80 J=2,NLmax
!CDIR NODEP
      DO I=1,N
      P(I)=(P(I)+AL(I,J)*X(LL(I,J)))*dble(IW(I))
      enddo
   80 enddo
!
!--------------------
C  Q=Q + AU*P 
!--------------------
!
      DO 100 J=1,NUmax
!CDIR NODEP
      DO I=1,N
      P(I)=(P(I)+AU(I,J)*X(LU(I,J)))*dble(IW(I))
      enddo
  100 enddo
!
!-----------------------------------------------------------
      RETURN
!-----------------------------------------------------------
      END SUBROUTINE AXMLT_B
!
!
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE DECOMP_LU_B(NCVIN,N2,NUmax,NLmax,ical_sld,
     &           multi_M,
     &           AL,LL,NL,IL,
     &           AU,LU,NU,IU,
     &           IT,IJT,IW2K_VECT,TEMP,
     &           MXALLCV,MXCV,N,LN,NO,DD,EPS,IER)
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
! Incomlete decomposition
!
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: IER,NCVIN,NUmax,NLmax
      integer,intent(in)    :: ical_sld
      real*8 ,intent(in)    :: EPS
      integer,intent(in)    :: MXALLCV,MXCV,N,NL,N2,NO
      real*8 ,intent(in)    :: AL(MXALLCV,0:NL)
      integer,intent(in)    :: LN(0:MXCV)
      integer,intent(in)    :: LL(MXALLCV,NL)

      integer,intent(in)    :: NU
      integer,intent(in)    :: LU(MXCV,NU)
      real*8 ,intent(in)    :: AU(MXCV,0:NU)
      integer,intent(in)    :: IL(MXCV)
      integer,intent(in)    :: IU(MXCV)
      real*8 ,intent(inout) :: DD(MXCV+N2)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      logical,INTENT(IN)    :: multi_M
      real*8 ,intent(inout)    :: TEMP(MXCV)      
!
! --- [local entities]
!
      integer :: I,IS,IE,J,K,Lnod,L
      real*8  :: EPS1,DW,dum1
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      EPS1=EPS*1.0D-3
      IER=0
      return
!-----------------------------------------
      DD(N+1:MXCV+N2)=0.d0
!
!      DO I=1,N
!!      DD(I)=1.d0/(DD(I))
!      do 231 j=1,IU(I)     !IQ(ICVL,1)+1,IQ(ICVL,2)
!      k=LU(I,j)
!      dum1=0.d0
!      do 232 l=1,IL(k)     !l=1,IQ(k,1)
!      IF(LL(k,l)==I) dum1=DD(I)*AL(LU(I,j),l)
! 232  enddo
!      DD(LU(I,j))=DD(LU(I,j))-dum1*AU(I,j)
! 231  enddo
!      ENDDO
!      return
!
!
      DO 100 J=1,NO
        IS=LN(J-1)+1
        IE=LN(J)
        IF((IE-IS+1)>0) then
!
!----------------------------------------
!
!CDIR NODEP
!          DO I=IS,IE 
!          DD(I)=1.d0/(DD(I))
!          enddo 
!
          DO I=IS,IE 
          DD(I)=DD(I)-AL(I,1)*DD(LL(I,1))*AU(I,1) !'incomp' 
!          DD(I)=DD(I)-AL(I,1)**2*DD(LL(I,1))
!          DD(I)=DD(I)-AL(I,1)*DD(LL(I,1))*DD(I) 
!          DD(I)=DD(I)-AL(I,1)*DD(LL(I,1))     
          enddo 
!
          DO 60 K=2,NLmax-1 
!CDIR NODEP
          DO I=IS,IE   !SSSSS
          DD(I)=DD(I)-AL(I,K)*DD(LL(I,K))*AU(I,K) 
!          DD(I)=DD(I)-AL(I,K)**2*DD(LL(I,K))
!          DD(I)=DD(I)-AL(I,K)*DD(LL(I,K))*DD(I) 
!          DD(I)=DD(I)-AL(I,K)*DD(LL(I,K))
          enddo
   60     enddo
!
!CDIR NODEP
!
          DO I=IS,IE   !SSSSS
          DW=DD(I)-AL(I,NLmax)*DD(LL(I,NLmax))*AU(I,NLmax) 
!          DW=DD(I)-AL(I,NLmax)**2*DD(LL(I,NLmax))
!          DW=DD(I)-AL(I,NLmax)*DD(LL(I,NLmax))*DD(I) 
!          DW=DD(I)-AL(I,NLmax)*DD(LL(I,NLmax))
          IF (DABS(DW).LT.EPS1) THEN
            IER = 1
            WRITE(*,*) '(SUBR. VICCG_B) Singular at step = ', I
            RETURN
          ENDIF
          DD(I)=1.0D0/DW
          enddo
        endif
  100 enddo
!
!      DO I=1,N  !5555
!      DD(I)=1.d0/DD(I)
!      enddo
!      CALL SEND_RECIVE_B(multi_M,ical_sld,
!     & IW2K_VECT,TEMP,MXCV,N,1,MXCV+N2,IT,IJT,DD)
!
      END SUBROUTINE DECOMP_LU_B
!
CCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE LUSUBB_B(
     &      m,N2,DD,AL,LL,NL,AU,LU,NU,NCVIN,multi_M,
     &      ical_sld,NUmax,NLmax,IW2K_VECT,TEMP,IW,
     &      MXALLCV,MXCV,N,LN,NO,IT,IJT,R,Q)
CCCCSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MXALLCV,MXCV,N2,N,NL,NU,
     &                         NO,NCVIN,m,
     &                         NUmax,NLmax
      real*8 ,intent(inOUT) :: AL(MXALLCV,0:NL)
      real*8 ,intent(in)    :: AU(MXCV,0:NU)
      real*8 ,intent(in)    :: DD(MXCV+N2)
      real*8 ,intent(inout) :: Q(MXCV+N2)
      real*8 ,intent(inout) :: R(MXCV+N2)
      integer,intent(in)    :: LL(MXALLCV,NL),LU(MXCV,NU)
      integer,intent(in)    :: LN(0:MXCV)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      real*8 ,intent(inout) :: TEMP(MXCV)
      integer,intent(in)    :: ical_sld
      logical,intent(in)    :: multi_M
      integer,intent(in)    :: IW(MXCV)            
!
! --- [local entities]
!
      integer :: I,IS,IE,J,K,ISUM
      real*8  :: SUM
!
!-----------------
C  Q=(L)**(-1)*R 
!-----------------
      DO 100 J=1,NO
        IS=LN(J-1)+1
        IE=LN(J)
        IF((IE-IS+1)>0) then
!CDIR NODEP
          DO 20 I=IS,IE
          Q(I)=(R(I)-AL(I,1)*Q(LL(I,1)))*dble(IW(I))
 20       enddo
!
          DO 60 K=2,NLmax-1
!CDIR NODEP
          DO I=IS,IE
          Q(I)=(Q(I)-AL(I,K)*Q(LL(I,K)))*dble(IW(I))
          enddo
   60     enddo
!
!CDIR NODEP
          DO I=IS,IE
          Q(I)=(DD(I)*(Q(I)-AL(I,NLmax)*Q(LL(I,NLmax))))
     &         *dble(IW(I))
          enddo
        endif
  100 enddo
!-----------------
C  Q= (U)**(-1)*Q 
!-----------------
!
      CALL SEND_RECIVE_B(multi_M,ical_sld,
     & IW2K_VECT,TEMP,MXCV,N,1,MXCV+N2,IT,IJT,Q)
!-----------------------------------------------------------
      DO 200 J=NO,1,-1         !NO-1,1,-1
        IS=LN(J-1)+1
        IE=LN(J)
        DO 140 K=1,NUmax
        IF((IE-IS+1)>0) then
!CDIR NODEP
          DO 120 I=IS,IE 
          Q(I)=(Q(I)-DD(I)*AU(I,K)*Q(LU(I,K)))*dble(IW(I))
 120      enddo
        endif
 140    enddo
 200  enddo
!-----------------------------------------------------------
      CALL SEND_RECIVE_B(multi_M,ical_sld,
     & IW2K_VECT,TEMP,MXCV,N,1,MXCV+N2,IT,IJT,Q)
!
!-----------------------------------------------------------
      RETURN
!
      END SUBROUTINE LUSUBB_B
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SEND_RECIVE_B(multi_M,ical_sld,
     & IW2K_VECT,TEMP,MXCV,N,INI,MXCVN2,IT,IJT,WORK)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: ical_sld,MXCV,N,INI,MXCVN2
      logical,intent(in)    :: multi_M
      integer,intent(in)    :: IW2K_VECT(MXCV,2)
      real*8 ,intent(inOUT) :: TEMP(MXCV)
      real*8 ,intent(inOUT) :: WORK(INI:MXCVN2)
      integer,intent(in)    :: IJT(0:MXCV)
      integer,intent(in)    :: IT(MXCV)
!
! --- [local entities]
!
      integer :: I
!
!
!
      if(NPE.gt.1) then
        DO I=1,N
          TEMP(I)=WORK(IJT(I))!*dble(IW(IJT(I)))
        ENDDO
        call SOLVER_SEND_RECV0(1,1,MXCV,N,TEMP)
        DO I=1,N
          WORK(I)=TEMP(IT(I))
        enddo
      endif
!
      return
!
      end SUBROUTINE SEND_RECIVE_B
