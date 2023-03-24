!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine user_output(
     &         MXMAT,MXCOMP,MXRANS,MXCV,MXALLCV,
     &         MXALLCVR,NMAT,NCOMP,NRANS,
     &         MXALLCV2,MXALLCV_RAD,
     &         NCV,nflud,
     &         nofld,MAT_NO,MAT_CVEXT,MAT_DCIDX,LCV_CV,
     &         NPE,my_rank,root,
     &         CVCENT,CVVOLM,disall,
     &         iter,time,deltt,spcnam,wm,
     &         aks,
     &         prs,rho,tmp,yys,vel,rmut,RVA,
     &         prs2,rho2,tmp2,yys2,vel2,rmut2,RVA2,
     &         nbcnd,nobcnd,MXSSFBC,MXCVFAC,MXCVFAC2,
     &         LBC_INDEX,LBC_SSF,LVEDGE,LCYCSF,SFAREA,SFCENT,locmsh,
     &         HTFLUX,htcbnd,mtcbnd,
     &         rmu,utau,cp,cps,
     &         radflag,radqw)     
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     User subroutine for output result
!
!=====================================================================
!---------------------------------------------------------------------
! --- DIMENSION 
!---------------------------------------------------------------------
!
!     MXMAT    : Material number      (array's dimension)
!     MXCOMP   : Species number       (array's dimension)
!     MXALLCV  : CV number            (array's dimension)
!     MXALLCVR : CV number for scalar (array's dimension)
!     NMAT     : Material number      (do-loop up-limit)
!     NCOMP    : Total Species number (do-loop up-limit)
!     NRANS    : Total Scalar transfer equation's number
!                                     (do-loop up-limit)
!     nflud    : fluid number         (array's dimension)
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
!     LCV_CV    (  MXALLCV)   : Globle CV No. to Local material No.
!
!---------------------------------------------------------------------
! --- Computation info
!---------------------------------------------------------------------
!
!     iter      : Time step iter No.
!     time      : time
!     deltt     : time interval
!     NPE       : Total cpu number 
!     my_rank   : My cpu number 
!     spcnam(MXCOMP) : name of species
!
!---------------------------------------------------------------------
! --- PHYSICAL VABLES
!---------------------------------------------------------------------
!
!     prs       : Pressure at  Control Volume center
!     rho       : Density  at  Control Volume center
!     tmp       : Temperature at  Control Volume center
!     yys       : Species mass fraction  at  Control Volume center
!     vel       : Velocity  at  Control Volume center
!     aks       : Scalar value (K,e etc.) at  Control Volume center 
!
!=====================================================================
      use module_usersub,only     : usr_dum1,usr_dum2,usr_dum3,
     &                              usr_dum4,usr_dum5,usr_dum6
      use module_boundary,only    : MAT_BCIDX,boundName
      use module_fluidforce,only  : ForceFlag,iwallvar
      use module_turbparm,only    : prturb,scturb
      use module_material ,only   : prlmnr
!
      implicit none
!
! --- [dummy arguments] 
!
      integer,intent(in)     :: MXMAT,MXCOMP,MXRANS,MXCV,MXALLCV,
     &                           MXALLCVR
      integer,intent(in)     :: nbcnd,MXSSFBC,MXCVFAC,MXCVFAC2
      integer,intent(in)     :: MXALLCV_RAD
      integer,intent(in)     :: MXALLCV2
      integer,intent(in)     :: NMAT,NCOMP,NRANS,NCV,nflud
      INTEGER,INTENT(IN)     :: MAT_NO    (0:MXMAT)
      INTEGER,INTENT(IN)     :: MAT_CVEXT (0:MXMAT)
      INTEGER,INTENT(IN)     :: MAT_DCIDX (0:MXMAT)
      INTEGER,INTENT(IN)     :: nofld     (  nflud)
      integer,intent(in)     :: LCV_CV    (MXALLCV)
      integer,intent(in)     :: LCYCSF   (   MXSSFBC)
      integer,intent(in)     :: NPE,my_rank,root
!
      REAL*8 ,INTENT(IN)     :: CVCENT  (3,MXALLCV)
      REAL*8 ,INTENT(IN)     :: CVVOLM  (  MXALLCV)
      REAL*8 ,INTENT(IN)     :: DISALL  (     MXCV)
!      
      integer,intent(in)     :: iter
      real*8 ,intent(in)     :: time,deltt
      character(*),intent(in):: spcnam(mxcomp)
      REAL*8 ,INTENT(IN)     :: WM(mxcomp)
      real*8 ,intent(in)     :: aks   (   MXALLCVR,MXRANS)
!   
      real*8 ,intent(in)     :: prs   (   MXALLCV)
      real*8 ,intent(in)     :: rho   (   MXALLCV)
      real*8 ,intent(in)     :: tmp   (   MXALLCV)
      real*8 ,intent(in)     :: yys   (   MXALLCV,MXCOMP)
      real*8 ,intent(inout)     :: vel   (   MXALLCV,3)  !????
      real*8 ,intent(in)     :: rmut  (   MXALLCV)
      real*8 ,intent(in)     :: RVA   (   MXCVFAC)
!
      real*8 ,intent(in)     :: prs2  (   MXALLCV2)
      real*8 ,intent(in)     :: rho2  (   MXALLCV2)
      real*8 ,intent(in)     :: tmp2  (   MXALLCV2)
      real*8 ,intent(in)     :: yys2  (   MXALLCV2,MXCOMP)
      real*8 ,intent(in)     :: vel2  (   MXALLCV2,3)
      real*8 ,intent(in)     :: rmut2 (   MXALLCV2)
      real*8 ,intent(in)     :: RVA2  (   MXCVFAC2)
!
      integer,intent(in)     :: nobcnd(0:nbcnd)
      integer,intent(in)     :: LBC_INDEX(0:nbcnd)
      integer,intent(in)     :: LBC_SSF( MXSSFBC)
      integer,intent(in)     :: LVEDGE(2,MXCVFAC)
      real*8 ,intent(in)     :: SFAREA(4,MXCVFAC)
      real*8 ,intent(in)     :: SFCENT(3,MXCVFAC)
      integer,intent(in)     :: locmsh(  MXSSFBC)
      integer,intent(in)     :: radflag
      real*8 ,intent(in)     :: radqw (MXALLCV_RAD)
!
      real*8 ,intent(in)     :: rmu   (   MXALLCV)
      REAL*8 ,INTENT(IN)     :: utau  ( 0:MXSSFBC)
!
      real*8 ,intent(in)     :: cps   (   MXALLCV,MXcomp)
      real*8 ,intent(in)     :: cp    (   MXALLCV)
!
      real*8 ,intent(in)     :: htcbnd    (    MXSSFBC)
      real*8 ,intent(in)     :: HTFLUX    (    MXSSFBC)
      real*8 ,intent(in)     :: mtcbnd    (    MXSSFBC)
!---------------------------------------------------------------
      integer,save :: iflag=0,iflagg=0
      integer      :: ICOM,nb,IIMAT,kd,IBFS,IBFE,ICFL,IBFL
      real*8       :: dum1,dum2,dum3,dum4,dum(35)
!

      integer,parameter :: outbcID=2,iter_s=1000,iter_ss=10010
      real(8),parameter :: outbcfact=1.d0
      real(8),save :: totalArea=-1.0d0    
      integer      :: nb1,ICV,IDC,ICVS,ICVE,ICVL,ios
      real(8)      :: af1,af2,v11,v12,v13,v21,v22,v23,dumm=0.d0
!
!
!
      real(8)::yp,rnu,utaul,yplus,xxxx,yyyy,zzzz
      integer   ::no,i
      real(8):: heat,minx,x,xx,mdiff,dd1,dd2,wi1,wi2,dd0
!---------------------------------------------------------------
!---------------------------------------------------------------


      if(mod(iter,5)==0) then
      open(2345,file="thrmflx.dat",form='formatted',status='replace')
      do nb=1,nbcnd
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      no=nobcnd(nb)
      if(no==1) then   !wall
        do IBFL=IBFS,IBFE
           ICFL=LBC_SSF(IBFL)
           ICV=LVEDGE(1,ICFL)
           IDC=LVEDGE(2,ICFL)
           xxxx=CVCENT(1,IDC)
           if(CVCENT(3,ICV).le.0.5)then
             write(2345,1001) xxxx,-HTFLUX(IBFL)*1.0d-6
           endif
!           if(IBFL==IBFS+100) print*,'USER-1',-HTFLUX(IBFL)
        enddo
      endif
      enddo
 1001 format(2(e13.5))
      close(2345)

      endif
      

      return





      if(iflag==0) then
         open(60,FILE='smoke',
     &       form='formatted',status='unknown',iostat=ios)
         iflag=1
      endif
      if(mod(iter,50)==0) then
        do nb=1,nbcnd
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        if(nb==9) then 
          dum1=0.d0
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum2=(vel(ICV,1)*SFAREA(1,ICFL)
     &         +vel(ICV,2)*SFAREA(2,ICFL)
     &         +vel(ICV,3)*SFAREA(3,ICFL))*SFAREA(4,ICFL)
          dum2=max(dum2,0.d0)
          dum1=dum1+yys(ICV,2)*rho(ICV)*dum2
          enddo
        elseif(nb==8) then
          dum4=0.d0
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          dum2=(vel(ICV,1)*SFAREA(1,ICFL)
     &         +vel(ICV,2)*SFAREA(2,ICFL)
     &         +vel(ICV,3)*SFAREA(3,ICFL))*SFAREA(4,ICFL)
          dum4=dum4+yys(ICV,2)*rho(ICV)*dum2
          enddo
        endif
        enddo
        write(60,'(1x,I8,a,3E14.6)') iter,'RATE=',dum1,dum4,dum1/dum4
      endif
!      

      return
!
      end subroutine user_output
