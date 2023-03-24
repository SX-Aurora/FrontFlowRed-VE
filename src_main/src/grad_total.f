!
!     subroutine grad_cell
!     subroutine grad_cell_body
!     subroutine grad_cell_least
!     subroutine AGMOR
!     subroutine BMAQR
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine grad_cell(ndim,imode,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,phi,grdc) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1.  Calculate gradient at cell center
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
      use module_metrix,only   : bdyf
      use MODULE_ARRAY_FFLOW,only   : MAT_INDEX
      
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim,imode
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real*8 ,intent(in)  :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)  :: wiface    (  MXCVFAC)
      real*8 ,intent(in)  :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)  :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)  :: phi       (  MXALLCV,ndim)
      real*8 ,intent(out) :: grdc      (  MXALLCV,3,ndim)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IDIM,iv
      integer :: IMAT_U,myid,IE
      real*8  :: dum1,dum2,dum3,wi1,wi2,grdf(3),dx,dy,dz,dl
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-30,ZERO=0.D0
!
      do IDIM=1,ndim
      grdc(1:MXALLCV,1:3,IDIM)=0.d0
      enddo
!------------------------------------------------
! --- Cell center gradient using Gauss' theorem
! --- Gauss' theorem: du/dxi=SUMj(ujSji)/Vol;i=1,2,3(x,y,z);j=1,N 
!------------------------------------------------
      if(NPE.gt.1) then
        call SOLVER_SEND_RECV(ndim,MXALLCV,NCV,phi)
      endif
!----------------------------------
! --- 
!----------------------------------
! DEBUG
      if(ical_vect) then
        do IDIM=1,ndim
!
        do IE=1,MAXIE   !vctr(ICVL,0)
!CIDR NODEP
#ifdef SX_TUNE
!NEC$ retain(phi)
!NEC$ retain(vctr)
!NEC$ gather_reorder
!NEC$ vob
!NEC$ vovertake
        DO ICVL=ICVS_V,ICVE_V
        ICFL=vctr(ICVL,IE)
        if(ICFL==0) cycle
        if(ICFL<0) then
          grdc(ICVL,1,IDIM)=grdc(ICVL,1,IDIM)
     &       +SFAREA(1,-ICFL)
     &       *(wiface(-ICFL)*phi(ICVL,IDIM)
     &       +(1.d0-wiface(-ICFL))*phi(LVEDGE(2,-ICFL),IDIM))
     &       *SFAREA(4,-ICFL)
          grdc(ICVL,2,IDIM)=grdc(ICVL,2,IDIM)
     &       +SFAREA(2,-ICFL)
     &       *(wiface(-ICFL)*phi(ICVL,IDIM)
     &       +(1.d0-wiface(-ICFL))*phi(LVEDGE(2,-ICFL),IDIM))
     &       *SFAREA(4,-ICFL)
          grdc(ICVL,3,IDIM)=grdc(ICVL,3,IDIM)
     &       +SFAREA(3,-ICFL)
     &       *(wiface(-ICFL)*phi(ICVL,IDIM)
     &       +(1.d0-wiface(-ICFL))*phi(LVEDGE(2,-ICFL),IDIM))
     &       *SFAREA(4,-ICFL)
!
        elseif(ICFL>0) then
          grdc(ICVL,1,IDIM)=grdc(ICVL,1,IDIM)
     &       -SFAREA(1,ICFL)
     &       *((1.d0-wiface(ICFL))*phi(ICVL,IDIM)
     &       +wiface(ICFL)*phi(LVEDGE(1,ICFL),IDIM))
     &       *SFAREA(4,ICFL)
          grdc(ICVL,2,IDIM)=grdc(ICVL,2,IDIM)
     &       -SFAREA(2,ICFL)
     &       *((1.d0-wiface(ICFL))*phi(ICVL,IDIM)
     &       +wiface(ICFL)*phi(LVEDGE(1,ICFL),IDIM))
     &       *SFAREA(4,ICFL)
          grdc(ICVL,3,IDIM)=grdc(ICVL,3,IDIM)
     &       -SFAREA(3,ICFL)
     &       *((1.d0-wiface(ICFL))*phi(ICVL,IDIM)
     &       +wiface(ICFL)*phi(LVEDGE(1,ICFL),IDIM))
     &       *SFAREA(4,ICFL)
        endif
        enddo
#else
        DO ICVL=ICVS_V,ICVE_V
        ICFL=vctr(ICVL,IE)
        if(ICFL==0) cycle
        if(ICFL<0) then
          grdc(ICVL,1,IDIM)=grdc(ICVL,1,IDIM)
     &       +SFAREA(1,-ICFL)
     &       *(wiface(-ICFL)*phi(ICVL,IDIM)
     &       +(1.d0-wiface(-ICFL))*phi(LVEDGE(2,-ICFL),IDIM))
     &       *SFAREA(4,-ICFL)
          grdc(ICVL,2,IDIM)=grdc(ICVL,2,IDIM)
     &       +SFAREA(2,-ICFL)
     &       *(wiface(-ICFL)*phi(ICVL,IDIM)
     &       +(1.d0-wiface(-ICFL))*phi(LVEDGE(2,-ICFL),IDIM))
     &       *SFAREA(4,-ICFL)
          grdc(ICVL,3,IDIM)=grdc(ICVL,3,IDIM)
     &       +SFAREA(3,-ICFL)
     &       *(wiface(-ICFL)*phi(ICVL,IDIM)
     &       +(1.d0-wiface(-ICFL))*phi(LVEDGE(2,-ICFL),IDIM))
     &       *SFAREA(4,-ICFL)
!
        elseif(ICFL>0) then
          grdc(ICVL,1,IDIM)=grdc(ICVL,1,IDIM)
     &       -SFAREA(1,ICFL)
     &       *((1.d0-wiface(ICFL))*phi(ICVL,IDIM)
     &       +wiface(ICFL)*phi(LVEDGE(1,ICFL),IDIM))
     &       *SFAREA(4,ICFL)
          grdc(ICVL,2,IDIM)=grdc(ICVL,2,IDIM)
     &       -SFAREA(2,ICFL)
     &       *((1.d0-wiface(ICFL))*phi(ICVL,IDIM)
     &       +wiface(ICFL)*phi(LVEDGE(1,ICFL),IDIM))
     &       *SFAREA(4,ICFL)
          grdc(ICVL,3,IDIM)=grdc(ICVL,3,IDIM)
     &       -SFAREA(3,ICFL)
     &       *((1.d0-wiface(ICFL))*phi(ICVL,IDIM)
     &       +wiface(ICFL)*phi(LVEDGE(1,ICFL),IDIM))
     &       *SFAREA(4,ICFL)
        endif
        enddo
#endif /** SX_TUNE **/
        enddo
!
!CIDR NODEP
        DO ICVL=ICVS_V,ICVE_V
        grdc(ICVL,1,IDIM)=grdc(ICVL,1,IDIM)/CVVOLM(ICVL)
        grdc(ICVL,2,IDIM)=grdc(ICVL,2,IDIM)/CVVOLM(ICVL)
        grdc(ICVL,3,IDIM)=grdc(ICVL,3,IDIM)/CVVOLM(ICVL)
        ENDDO
        ENDDO
      else
        do 200 IDIM=1,ndim
        do 100 IIMAT=1,NMAT            !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dum1=(wi1*phi(ICVLA,IDIM)
     &       +wi2*phi(ICVLB,IDIM))
     &       *SFAREA(4,ICFL)

!
        do 101 iv=1,3
        dum2=dum1*SFAREA(IV,ICFL)
        grdc(ICVLA,IV,IDIM)=grdc(ICVLA,IV,IDIM)+dum2
        grdc(ICVLB,IV,IDIM)=grdc(ICVLB,IV,IDIM)-dum2
 101    enddo
        enddo
 100    enddo
!
        do 410 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
!        ICVE=MAT_INDEX(IIMAT)
        ICVE=MAT_CVEXT(IIMAT)
        do 111 iv=1,3
        do 110 ICVL=ICVS,ICVE
        dum1=1.d0/(CVVOLM(ICVL)+SML)
        grdc(ICVL,iv,IDIM)=grdc(ICVL,iv,IDIM)*dum1
!     &             -bdyf(ICVL,iv)*dum1  !phi(ICVL,IDIM)
 110    enddo
 111    enddo
 410    enddo
 200    enddo
      endif
!----------------------------------
! --- MPI 
!----------------------------------
      if(NPE.gt.1) then
        DO idim=1,ndim
        call SOLVER_SEND_RECV(3,MXALLCV,NCV,grdc(:,:,idim))
        ENDDO
      endif
!
      return
      end subroutine grad_cell
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine grad_cell_dir(ndim,imode,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,phi,grdc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1.  Calculate gradient at cell center
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim,imode
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real*8 ,intent(in)  :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)  :: wiface    (  MXCVFAC)
      real*8 ,intent(in)  :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)  :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)  :: phi       (  MXALLCV,3,ndim)
      real*8 ,intent(out) :: grdc      (  MXALLCV,ndim)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IDIM,iv
      integer :: IMAT_U,myid,IE
      real*8  :: dum1,dum2,wi1,wi2,grdf(3),dx,dy,dz,dl,dum3,dum4
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
      do IDIM=1,ndim
      grdc(1:MXALLCV,IDIM)=0.d0
      enddo
!------------------------------------------------------------------
! --- Cell center gradient using Gauss' theorem
! --- Gauss' theorem: du/dxi=SUMj(ujSji)/Vol;i=1,2,3(x,y,z);j=1,N 
!------------------------------------------------------------------
!      if(NPE.gt.1) then
!        do IDIM=1,ndim
!        call SOLVER_SEND_RECV(3,MXALLCV,NCV,phi(:,1:3,IDIM))
!        enddo
!      endif

! --- 
!----------------------------------
!
      do 200 IDIM=1,ndim
      do 100 IIMAT=1,NMAT            !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=(1.d0-wiface(ICFL))

        grdf(1)=SFAREA(1,ICFL)
        grdf(2)=SFAREA(2,ICFL)
        grdf(3)=SFAREA(3,ICFL)

        do 101 IV=1,3
        dum1=wi1*phi(ICVLA,IV,IDIM)+wi2*phi(ICVLB,IV,IDIM)
        dum4=dum1*grdf(IV)*SFAREA(4,ICFL)!*SFAREA(iv,ICFL)
        grdc(ICVLA,IDIM)=grdc(ICVLA,IDIM)+dum4
        grdc(ICVLB,IDIM)=grdc(ICVLB,IDIM)-dum4
 101    enddo

        enddo
 100  enddo
!

      do 410 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 110 ICVL=ICVS,ICVE
        dum1=1.d0/(CVVOLM(ICVL)+SML)
        grdc(ICVL,IDIM)=grdc(ICVL,IDIM)*dum1
 110    enddo
 410  enddo
 200  enddo
!----------------------------------
! --- MPI
!----------------------------------okabe
      if(NPE.gt.1) then
        do IDIM=1,ndim
        call SOLVER_SEND_RECV(1,MXALLCV,NCV,grdc(:,IDIM))
        enddo
      endif
!
      return
      end subroutine grad_cell_dir
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine grad_cell_dir1(ndim,imode,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,phi,grdc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1.  Calculate gradient at cell center
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim,imode
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real*8 ,intent(in)  :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)  :: wiface    (  MXCVFAC)
      real*8 ,intent(in)  :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)  :: SFCENT    (3,MXCVFAC)
      real(8),intent(in)  :: CVCENT    (3,MXALLCV)
      real*8 ,intent(in)  :: phi       (  MXALLCV,ndim)
      real*8 ,intent(out) :: grdc      (  MXALLCV,ndim)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICV,IDC,ICVP,IDCP
      integer :: IDIM,iv
      integer :: IMAT_U,myid,IE
      real*8  :: dum1,dum2,wi1,wi2,grdf(3),dx,dy,dz,dl,dum3,dum4
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
      do IDIM=1,ndim
      grdc(1:MXALLCV,IDIM)=0.d0
      enddo
!
      do 200 IDIM=1,ndim
      do 100 IIMAT=1,NMAT            !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=(1.d0-wiface(ICFL))
        dx=CVCENT(1,ICVLB)-CVCENT(1,ICVLA)
        dy=CVCENT(2,ICVLB)-CVCENT(2,ICVLA)
        dz=CVCENT(3,ICVLB)-CVCENT(3,ICVLA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz+SML)  
        grdf(1)=SFAREA(1,ICFL)
        grdf(2)=SFAREA(2,ICFL)
        grdf(3)=SFAREA(3,ICFL)

!        do 101 IV=1,3
        dum1=(phi(ICVLB,IDIM)-phi(ICVLA,IDIM))/dl!*grdf(IV)
     &    *SFAREA(4,ICFL)
        grdc(ICVLA,IDIM)=grdc(ICVLA,IDIM)+dum1
        grdc(ICVLB,IDIM)=grdc(ICVLB,IDIM)-dum1
! 101    enddo

        enddo
 100  enddo
!

      do 410 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 110 ICVL=ICVS,ICVE
        dum1=1.d0/(CVVOLM(ICVL)+SML)
        grdc(ICVL,IDIM)=grdc(ICVL,IDIM)*dum1
 110    enddo
 410  enddo



 200  enddo


!
      return
      end subroutine grad_cell_dir1

!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine grad_cell_body(ndim,imode,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,phi,bdyfrc,grdc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1.  Calculate gradient at cell center
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim,imode
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real*8 ,intent(in)  :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)  :: wiface    (  MXCVFAC)
      real*8 ,intent(in)  :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)  :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)  :: phi       (  MXALLCV,ndim)
      real*8 ,intent(out) :: grdc      (  MXALLCV,3,ndim)
      real*8 ,intent(in)  :: bdyfrc(MXCV_B,3)
      real*8 ,intent(in)  :: CVCENT    (3, MXALLCV)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IDIM,iv
      integer :: IMAT_U,myid,IE
      real*8  :: dum1,dum2,wi1,wi2,grdf(3),dx,dy,dz,dl,dx2,dy2,dz2
      real*8  :: dum3,dum4,dum5
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
      do IDIM=1,ndim
      grdc(:,:,IDIM)=0.d0
      enddo
!-----------------------------------------------------------------
! --- Cell center gradient using Gauss' theorem
! --- Gauss' theorem: du/dxi=SUMj(ujSji)/Vol;i=1,2,3(x,y,z);j=1,N 
!-----------------------------------------------------------------
      if(NPE.gt.1) then
        call SOLVER_SEND_RECV(ndim,MXALLCV,NCV,phi)
      endif
!----------------------------------
! --- 
!----------------------------------
      if(ical_vect) then
        call FFRABORT(1,'ERR:  grad_cell_body')
      else
        do 200 IDIM=1,ndim
        do 100 IIMAT=1,NMAT    !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        wi1=wiface(ICFL)
        wi2=1.d0-wiface(ICFL)
        dx=-(SFCENT(1,ICFL)-CVCENT(1,ICVLA))
        dy=-(SFCENT(2,ICFL)-CVCENT(2,ICVLA))
        dz=-(SFCENT(3,ICFL)-CVCENT(3,ICVLA))
        dx2=-(SFCENT(1,ICFL)-CVCENT(1,ICVLB))
        dy2=-(SFCENT(2,ICFL)-CVCENT(2,ICVLB))
        dz2=-(SFCENT(3,ICFL)-CVCENT(3,ICVLB))
!----------------------------------------------------------------
!        grdf(1)=dx*wi1+dx2*wi2
!        grdf(2)=dy*wi1+dy2*wi2
!        grdf(3)=dz*wi1+dz2*wi2
        grdf(1)=wi1*dx*bdyfrc(ICVLA,1)+wi2*dx2*bdyfrc(ICVLB,1)
        grdf(2)=wi1*dy*bdyfrc(ICVLA,2)+wi2*dy2*bdyfrc(ICVLB,2)
        grdf(3)=wi1*dz*bdyfrc(ICVLA,3)+wi2*dz2*bdyfrc(ICVLB,3)
!        dum4=grdf(1)*dum1+grdf(2)*dum2+grdf(3)*dum3
        dum5=(wi1*phi(ICVLA,IDIM)+wi2*phi(ICVLB,IDIM)!-dum4
     &        )*SFAREA(4,ICFL)
!----------------------------------------------------------------
        do 101 iv=1,3
        dum4=dum5*SFAREA(IV,ICFL)!-grdf(iv)*SFAREA(4,ICFL)
        grdc(ICVLA,IV,IDIM)=grdc(ICVLA,IV,IDIM)+dum4
        grdc(ICVLB,IV,IDIM)=grdc(ICVLB,IV,IDIM)-dum4
 101    enddo
        enddo
 100    enddo
!
        do 410 IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do 111 iv=1,3
        do 110 ICVL=ICVS,ICVE
        dum1=1.d0/(CVVOLM(ICVL)+SML)
!        grdc(ICVL,iv,IDIM)=grdc(ICVL,iv,IDIM)*dum1
        grdc(ICVL,iv,IDIM)=grdc(ICVL,iv,IDIM)*dum1+bdyfrc(ICVL,iv)
 110    enddo
 111    enddo
 410    enddo
 200    enddo
      endif
!----------------------------------
! --- MPI 
!----------------------------------
      if(NPE.gt.1) then
        DO idim=1,ndim
        call SOLVER_SEND_RECV(3,MXALLCV,NCV,grdc(:,:,idim))
        ENDDO
      endif
!
      return
      end subroutine grad_cell_body
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine grad_cell_least(ndim,imode,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,rva,phi,grdc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_model,only    : ical_vect,nthrds
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V,index_c,index_f
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim,imode
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real*8 ,intent(in)  :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)  :: wiface    (  MXCVFAC)
      real*8 ,intent(in)  :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)  :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(in)  :: phi       (  MXALLCV,ndim)
      real*8 ,intent(out) :: grdc      (  MXALLCV,3,ndim)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
      real*8 ,intent(in)  :: CVCENT    (3, MXALLCV)
      real*8 ,intent(in)  :: rva       (   MXCVFAC)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,LL=1
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP,ICFLL
      integer :: IDIM,iv
      integer :: IMAT_U,myid,IE,IC
      real*8  :: dum1,dum2,wi1,wi2,grdf(3),dx,dy,dz,dl,C(3),dlvect
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8,allocatable :: A(:,:),Q(:,:),B(:)
!
      do IDIM=1,ndim
      grdc(:,:,IDIM)=0.d0
      enddo
      allocate(A(MAXIE,3),Q(MAXIE,MAXIE),B(MAXIE))
      A(:,:)=0.d0
      Q(:,:)=0.d0
      B(:)=0.d0
      C(:)=0.d0
!----------------------------------
! --- Cell center gradient using Gauss' theorem
! --- Gauss' theorem: du/dxi=SUMj(ujSji)/Vol;i=1,2,3(x,y,z);j=1,N 
!----------------------------------
      if(NPE.gt.1) then
        call SOLVER_SEND_RECV(ndim,MXALLCV,NCV,phi)
      endif
!----------------------------------
! --- 
!----------------------------------
      if(ical_vect) then
        call FFRabort
     &  (1,'ERR: Least square methord NOT support Vector CPU ')
      else
        do 200 IDIM=1,ndim
        do 100 IIMAT=1,NMAT    !ICF=1,NCVFAC 
        if(.not.mat_cal(IIMAT)) cycle
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        M=vctr(ICVL,0)
        do IE=1,M
        ICFL=vctr(ICVL,IE)
        ICFLL=ABS(ICFL)
        IC=(3-ICFL/ICFLL)/2
!        IC=(3+ICFL/ICFLL)/2
        ICVB=LVEDGE(IC,ICFLL)
        ICVA=ICVL
        dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
        dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
        dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
        dl=dsqrt(dx*dx+dy*dy+dz*dz)
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl
        dlvect=abs(dx*SFAREA(1,ICFLL) 
     &            +dy*SFAREA(2,ICFLL) 
     &            +dz*SFAREA(3,ICFLL)) 
!        dum1=dble(ICFL/ICFLL)  !*SFAREA(4,ICFLL) 
        dum1=-dble(ICFL/ICFLL)  !*SFAREA(4,ICFLL) 
!        B(IE)=dlvect*(phi(ICVB,IDIM)-phi(ICVA,IDIM))/dl!*SFAREA(4,ICFLL) 
        B(IE)=(phi(ICVB,IDIM)-phi(ICVA,IDIM))/(dl*dlvect)!*SFAREA(4,ICFLL) 
        A(IE,1:3)=SFAREA(1:3,ICFLL)*dum1 
        enddo
        CALL AGMOR(A,MAXIE,M,3,B,Q,LL,C) 
        if(LL==0) call FFRABORT(1,'ERR:LL=0 in grad_cell_least') 
        do 101 iv=1,3
        grdc(ICVL,IV,IDIM)=B(iv)
 101    enddo
        enddo
 100    enddo
 200    enddo
!
      endif
!----------------------------------
! --- MPI
!----------------------------------
      if(NPE.gt.1) then
        DO idim=1,ndim
        call SOLVER_SEND_RECV(3,MXALLCV,NCV,grdc(:,:,idim))
        ENDDO
      endif
!
      deallocate(A,Q,B)
      return
      end subroutine grad_cell_least 
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine AGMOR(A,MAX,M,N,B,Q,L,C) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MAX,M,N  !M=N; N=3
      real*8 ,intent(in)    :: A(MAX,N)
      real*8 ,intent(inout) :: Q(MAX,MAX)
      real*8 ,intent(out)   :: B(MAX),C(N)
      integer,intent(inout)    :: L
!
! --- [local entities]
!
      real*8  :: D
      integer :: I,J


      CALL BMAQR(A,MAX,M,N,Q,L) 
      IF(L==0) RETURN
      DO 20 I=1,N
      D=0.D0 
      DO 10 J=1,M
      D=D+Q(J,I)*B(J)
 10   ENDDO
      C(I)=D
 20   ENDDO
!
      B(N)=C(N)/A(N,N)
      DO 40 I=N-1,1,-1
      D=0.D0
      DO 30 J=I+1,N
      D=D+A(I,J)*B(J)
 30   ENDDO
      B(I)=(C(I)-D)/A(I,I)
 40   ENDDO
!
      RETURN
      end subroutine AGMOR
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine BMAQR(A,MAX,M,N,Q,L)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 !
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil      
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: MAX,M,N  !M=N; N=3
      real*8 ,intent(inout) :: A(MAX,N)
      real*8 ,intent(inout) :: Q(MAX,MAX)
      integer,intent(inout) :: L
!
! --- [local entities]
!
      real*8  :: APPHA,T,U
      integer :: I,J,K,NN
!
      IF(M<N) THEN
        L=0
        RETURN
      ENDIF
!
      DO 10 I=1,M
      DO J=1,M
      Q(I,J)=0.D0
      IF(I==J) Q(I,J)=1.D0
      ENDDO
 10   ENDDO
!
      NN=N
      IF(M==N) NN=M-1
      DO 200 K=1,NN
      U=0.D0
      DO 20 I=K,M
      IF(ABS(A(I,K))>U)  U=ABS(A(I,K))
 20   ENDDO
      APPHA=0.D0
      DO 30 I=K,M
      T=A(I,K)/U
      APPHA=APPHA+T*T
 30   ENDDO
      IF(A(K,K)>0.D0) U=-U
      APPHA=U*SQRT(APPHA)
      IF(ABS(APPHA)+1.D0==1.D0) THEN
        L=0
        RETURN
      ENDIF
      U=SQRT(2.D0*APPHA*(APPHA-A(K,K)))
      IF(U+1.D0/=1.D0) THEN
        A(K,K)=(A(K,K)-APPHA)/U
        DO 50 I=K+1,M
        A(I,K)=A(I,K)/U
 50     ENDDO
        DO 80 J=1,M
        T=0.D0
        DO 60 L=K,M
        T=T+A(L,K)*Q(L,J)
 60     ENDDO
        DO 70 I=K,M
        Q(I,J)=Q(I,J)-2.D0*T*A(I,K)
 70     ENDDO

 80     ENDDO
        DO 110 J=K+1,N
        T=0.D0
        DO 90 L=K,M
        T=T+A(L,K)*A(L,J)
 90     ENDDO
        DO 100 I=K,M
        A(I,J)=A(I,J)-2.D0*T*A(I,K)
 100    ENDDO
 110    ENDDO
        A(K,K)=APPHA
        DO 120 I=K+1,M
        A(I,K)=0.D0
 120    ENDDO

      ENDIF
 200  ENDDO
      L=1
      DO 210 I=1,M-1
      DO 220 J=I+1,M
      T=Q(I,J)
      Q(I,J)=Q(J,I)
      Q(J,I)=T
 220  ENDDO
 210  ENDDO
!
      RETURN
!
      END subroutine BMAQR
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hat_bar(ndim,imode,dlalpha,R_deltV,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &  LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &  phi,grdc,mphi)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_les,only   : ign_iv
      use module_dimension
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim,imode
      real(8),intent(in)  :: dlalpha
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real(8),intent(in)  :: SFAREA    (4,MXCVFAC)
      real(8),intent(in)  :: wiface    (  MXCVFAC)
      real(8),intent(in)  :: CVVOLM    (  MXALLCV)
      real(8),intent(in)  :: SFCENT    (3,MXCVFAC)
      real(8),intent(in)  :: CVCENT    (3,MXALLCV)
      real(8),intent(in)  :: phi       (  MXALLCV,ndim)
      real(8),intent(out) :: grdc      (  MXALLCV,3,ndim)
      real(8),intent(out) :: mphi      (  MXALLCV,ndim)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      integer,intent(in)  :: LCYCOLD(  MXSSFBC_SLD)
      real*8 ,intent(in)  :: wifsld(   MXSSFBC_SLD)
      real*8 ,intent(in)  :: OPPANG(   MXSSFBC_SLD)
      real*8 ,intent(in)  :: R_deltV(  MXALLCV)
!
! --- [local entities]
!
      integer :: icvl,icvla,icvlb,nd,iv,n,ierr=0,IIMAT,IMAT,ICVS,ICVE
      real(8) :: delh,dum1
      real(8),parameter :: r1pn=1.d0/3.d0
      integer :: ICFS,ICFE,ICFL
      integer :: flag_g=1
!
!------------------------------------
! --- D^2(FAI)/D^2x ()
!------------------------------------
!
      if(flag_g==1) then
        call grad_cell(ndim,1,
     &     MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &     LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,phi,grdc)
!
        if(ndim==3) then
          call dc_vgrad(MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,
     &              LCYCOLD,wifsld,OPPANG,
     &              grdc)
        else
          do iv=1,3
          call dc_symprs(1,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &               mat_cal,grdc(:,iv,1))
          enddo
        endif
!
        call grad_cell_dir(ndim,1,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,grdc,mphi)
      else
!----------------------------------------------
! --- ! --- D^2(FAI)/D^2x (integer methord)
!----------------------------------------------
        call grad_cell_dir1(ndim,1,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,CVCENT,phi,mphi)
      endif
!
!----------------------------------------------
! 
!----------------------------------------------
      do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        if(IMAT<0) cycle
        do nd=1,ndim
        do ICVL=ICVS,ICVE
        delh=R_deltV(ICVL)*(CVVOLM(ICVL))**r1pn     !*dlalpha
        mphi(ICVL,nd)=phi(ICVL,nd)+delh**2/24.d0*mphi(ICVL,nd) 
        enddo
        enddo
      enddo
!
      return
!     
      end subroutine hat_bar
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine top_hat(ndim,imode,dlalpha,
     &  MAT_CV,MAT_CVEXT,MAT_NO,mat_cal,MAT_CFIDX,vctr,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,cvcent,
     &  LBC_SSF,LCYCSF,LCYCOLD,wifsld,OPPANG,
     &  phi,grdc,mphi)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_les,only   : ign_iv
      use module_dimension
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim,imode
      real(8),intent(in)  :: dlalpha
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real(8),intent(in)  :: SFAREA    (4,MXCVFAC)
      real(8),intent(in)  :: wiface    (  MXCVFAC)
      real(8),intent(in)  :: CVVOLM    (  MXALLCV)
      real(8),intent(in)  :: SFCENT    (3,MXCVFAC)
      real(8),intent(in)  :: CVCENT    (3,MXALLCV)
      real(8),intent(in)  :: phi       (  MXALLCV,ndim)
      real(8),intent(out) :: grdc      (  MXALLCV,3,ndim)
      real(8),intent(out) :: mphi      (  MXALLCV,ndim)
      integer,intent(in)  :: vctr(MXCV_V,0:MXBND_V)
      integer,intent(in)  :: LBC_SSF(  MXSSFBC)
      integer,intent(in)  :: LCYCSF (  MXSSFBC)
      integer,intent(in)  :: LCYCOLD(  MXSSFBC_SLD)
      real*8 ,intent(in)  :: wifsld(   MXSSFBC_SLD)
      real*8 ,intent(in)  :: OPPANG(   MXSSFBC_SLD)
!
! --- [local entities]
!
      integer :: icvl,icvla,icvlb,nd,iv,IDIM,ierr=0,IIMAT,IMAT,ICVS,ICVE
      integer :: ICFS,ICFE,ICFL
      real(8) :: delh,dum1,dum2
      real(8),parameter :: r1pn=1.d0/3.d0
      integer,parameter :: m_dor=3
!
!------------------------------------
! --- 
!------------------------------------
!
      do IDIM=1,ndim
      grdc(1:MXALLCV,1:3,IDIM)=0.d0
      enddo

      do 200 IDIM=1,ndim
        do 100 IIMAT=1,NMAT            !ICF=1,NCVFAC
        if(.not.mat_cal(IIMAT)) cycle
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=phi(ICVLA,IDIM)*CVVOLM(ICVLA)
     &      +phi(ICVLB,IDIM)*CVVOLM(ICVLB)
        dum2=CVVOLM(ICVLA)+CVVOLM(ICVLB)
        grdc(ICVLA,1,IDIM)=grdc(ICVLA,1,IDIM)+dum1
        grdc(ICVLB,1,IDIM)=grdc(ICVLB,1,IDIM)+dum1
        grdc(ICVLA,2,IDIM)=grdc(ICVLA,2,IDIM)+dum2
        grdc(ICVLB,2,IDIM)=grdc(ICVLB,2,IDIM)+dum2
        grdc(ICVLA,3,IDIM)=grdc(ICVLA,3,IDIM)+1.d0
        grdc(ICVLB,3,IDIM)=grdc(ICVLB,3,IDIM)+1.d0
        enddo
 100    enddo
!
        do iv=1,m_dor
        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        dum2=grdc(ICVL,3,IDIM)
        dum1=phi(ICVL,IDIM)*CVVOLM(ICVL)*dum2
        grdc(ICVL,1,IDIM)=grdc(ICVL,1,IDIM)+dum1
        grdc(ICVL,2,IDIM)=grdc(ICVL,2,IDIM)*dum2
        enddo
        enddo
        enddo

        do IIMAT=1,NMAT
        if(.not.mat_cal(IIMAT)) cycle
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do ICVL=ICVS,ICVE
        mphi(ICVL,IDIM)=grdc(ICVL,1,IDIM)/grdc(ICVL,2,IDIM)
        enddo
        enddo

 200  enddo
!
      end subroutine top_hat

!==========================================================================!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine corr_area(ndim,
     &  MAT_CV,MAT_CVEXT,MAT_NO,MAT_CFIDX,
     &  LVEDGE,SFAREA,SFCENT,wiface,CVVOLM,grdc)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1.  Calculate gradient at cell center
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_metrix,only    : bdyf
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: ndim
      INTEGER,INTENT(IN)  :: MAT_CV   (   MXALLCV)
      INTEGER,INTENT(IN)  :: MAT_CVEXT(   0:MXMAT)
      INTEGER,INTENT(IN)  :: MAT_NO   (   0:MXMAT)
!      logical,INTENT(IN)  :: mat_cal  (   0:MXMAT)
      integer,intent(in)  :: MAT_CFIDX(   0:MXMAT)
      integer,intent(in)  :: LVEDGE    (2,MXCVFAC)
      real*8 ,intent(in)  :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)  :: wiface    (  MXCVFAC)
      real*8 ,intent(in)  :: CVVOLM    (  MXALLCV)
      real*8 ,intent(in)  :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(out) :: grdc      (  MXALLCV,3)
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,IDCS,IDCE,ICFL,ICFS,ICFE
      integer :: ICVLA,ICVLB,ICVA,ICVB,ICV,IDC,ICVP,IDCP
      integer :: IDIM,iv
      integer :: IMAT_U,myid,IE
      real*8  :: dum1,dum2,wi1,wi2,grdf(3),dx,dy,dz,dl
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-30,ZERO=0.D0
!
      grdc(1:MXALLCV,1:3)=0.d0
!----------------------------------
! --- Cell center gradient using Gauss' theorem
! --- Gauss' theorem: du/dxi=SUMj(ujSji)/Vol;i=1,2,3(x,y,z);j=1,N 
!----------------------------------
      do 100 IIMAT=1,NMAT            !ICF=1,NCVFAC
        ICFS=MAT_CFIDX(IIMAT-1)+1
        ICFE=MAT_CFIDX(IIMAT)
        do ICFL=ICFS,ICFE
        ICVLA=LVEDGE(1,ICFL)
        ICVLB=LVEDGE(2,ICFL)
        dum1=SFAREA(4,ICFL)
        do 101 iv=1,3
        dum2=dum1*SFAREA(IV,ICFL)
        grdc(ICVLA,IV)=grdc(ICVLA,IV)+dum2
        grdc(ICVLB,IV)=grdc(ICVLB,IV)-dum2
 101    enddo
        enddo
 100  enddo
!
      
      do 410 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_CVEXT(IIMAT)
        do IV=1,3
        do ICVL=ICVS,ICVE
        bdyf(ICVL,IV)=grdc(ICVL,IV)!/(CVVOLM(ICVL)+SML)
        enddo
        enddo
 410  enddo

      return
      end subroutine corr_area
!
