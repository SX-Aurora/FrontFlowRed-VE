!
!      subroutine bc_check()
!      subroutine bc_kdbpv1
!      subroutine dc_metric
!      subroutine utl_boxmlr
!      subroutine utl_random
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_check
     & (mface,mcell,nface,mvrtx,nvrtx,
     & lbface,lcface,lvface,lacell,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only       : ifle
      use module_model,only    : icaltb,rns_scl
      use module_boundary,only : kdbcnd,kdprdc,kdsymm,kdintr,kdsld,
     &                           kdilet,kdolet,kdtchi,kxnone,kdpres,
     &                           kdstag,idis,kdbuff,kdshutr,kdpors,
     &                           kdovst
!
! 1. Check validity of boundary condition
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mface,mcell,nface,nvrtx,mvrtx
      integer,intent(in)  :: lbface(2,mface)
      integer,intent(in)  :: lcface(2,mface)
      integer,intent(in)  :: lvface(4,mface)
      integer,intent(in)  :: lacell(  mcell)
      integer,intent(out) :: ierror
!
! --- [local entities]
!
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
!
!
!
      ierror=0
!
      do 1000 IS=1,nface
      nb=lbface(1,IS)
      if(nb.gt.0) then
        kd=kdbcnd(0,nb)
        if(kd.ne.kdsymm.and.kd.ne.kdprdc.and.kd.ne.kdtchi.and.
     &     kd.ne.kdsld.and.kd/=kdovst) then
!
          IC1=lcface(1,IS)
          kdv=kdbcnd(1,nb)
          kdt=kdbcnd(2,nb)
          kdy=kdbcnd(3,nb)
          kdk=kdbcnd(4,nb)
          kdp=kxnone
!
!-< 1. fluid face >-
!
          IMAT=lacell(IC1)
          if(IMAT.gt.0) then
            if(kd.ne.kdilet
     &    .and.kd.ne.kdolet
     &    .and.kd.ne.kdpres
     &    .and.kd.ne.kdstag
     &    .and.kd/=kdbuff
     &    .and.kd/=kdshutr
     &    .and.kd/=kdpors
     &    ) then
              call errmsg1(1,kdv,IC1,'vel',nb)
              call errmsg1(1,kdt,IC1,'temp',nb)
              call errmsg1(1,kdy,IC1,'conc',nb)
              if(rns_scl) call errmsg1(1,kdk,IC1,'RANS',nb)
              if(ierror.ne.0 ) goto 9999
            endif
          endif
!
!-< 2. solid face >-
!
          if(IMAT.lt.0) then
            if(kd.eq.kdilet) kdv=kxnone+1
            if(kd.eq.kdolet) kdp=kxnone+1
            call errmsg1(2,kdt,IC1,'temp',nb)
            if(kd.ne.kdintr
     &    .and.kd/=kdbuff
     &    .and.kd/=kdshutr
     &    .and.kd/=kdpors
     &    ) then
              call errmsg2(2,kdv,IC1,'vel',nb)
              call errmsg2(2,kdp,IC1,'pressure',nb)
              call errmsg2(2,kdy,IC1,'conc',nb)
              if(rns_scl) call errmsg2(2,kdk,IC1,'RANS',nb)
            endif
            if(ierror.ne.0) goto 9999
          endif
!
        endif
      endif
 1000 continue
!
      return
!
 9999 continue
      ie=3+min(1,lvface(4,IS))
      write(ifle,*) 'vertices of the face 1=',(lvface(i,IS),i=1,ie)
      write(ifle,*) 'boundary condition no. =',nb
      write(ifle,*) '(bc_check)'
      ierror=1
!
!///////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine errmsg1(ifs,kd,n1,vnam,nb)
!=================================================
      integer     ,intent(in) :: ifs,kd,n1,nb
      character(*),intent(in) :: vnam
      character(5),parameter  :: cfs(2)=(/'fluid','solid'/)
!
      if( kd.ne.kxnone ) return
!
      write(ifle,*) '### error : data error'
      write(ifle,*) 'boundary condition for [',vnam,'] is not ',
     &              'specified on [',cfs(ifs),'] face'
      write(ifle,*) 'cell no. =',n1
      write(ifle,*) 'BC no=',nb
      ierror=1
      end subroutine errmsg1
!=================================================
      subroutine errmsg2(ifs,kd,n1,vnam,nb)
!=================================================
      integer     ,intent(in) :: ifs,kd,n1,nb
      character(*),intent(in) :: vnam
      character(5),parameter  :: cfs(2)=(/'fluid','solid'/)
!
      if( kd.eq.kxnone ) return
!
      write(ifle,*) '### error : data error'
      write(ifle,*) 'boundary condition for [',vnam,'] can not ',
     &              'be specified on [',cfs(ifs),'] face'
      write(ifle,*) 'cell no. =',n1
      write(ifle,*) 'BC no=',nb
      ierror=1
      end subroutine errmsg2
      end subroutine bc_check
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine bc_kdbpv1(mface,nface,lbface,kdbp)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_boundary,only : kdbcnd,kdprdc,kdsymm,kdilet,kdolet
!
! 1. Set boundary condition flag for pressure & velocity
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mface,nface
      integer,intent(in)  :: lbface(2,mface)
      integer,intent(out) :: kdbp  (  mface)
!
! --- [local entities]
!
      integer :: nb,kd,m1
      integer :: IC,IS,IV,IE,IVA,IVB,IC1,IC2
      integer :: ILS,ILV,ILE,IVV,IEE
      integer :: ICTP,ICOM,IMAT,IMD,IFL
      integer :: IBS,IPS,IDC
!
      kdbp=0
      do 1000 IS=1,nface
      nb=lbface(1,IS)
      if(nb.gt.0) then
        nb=lbface(1,IS)
        kd=kdbcnd(0,nb)
!
        if(kd.eq.kdsymm) then
          kdbp(IS)=2
        elseif(kd.ne.kdprdc) then
          if(kd.eq.kdolet) then
            kdbp(IS)=1
          else
            kdbp(IS)=2
          endif
        endif
!
      elseif(lbface(1,IS).lt.0) then
        kdbp(IS)=999
! --- have defined: kdbp(IS)=0
!      else
!       kdbp(IS)=0
      endif
 1000 continue
!
      return
      end subroutine bc_kdbpv1
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_metric
     & (NALLCV,NCVFAC,NSSFBC,nedge,medge,mssfbc,mvrtx,
     &  LVEDGE,LBCSSF,LCYCSF,SFCENT,SFAREA,CVCENT,CVVOLM,ICALL,icode)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_boundary,only : set_rotmtrx,nbcnd,
     &                           kdbcnd,kdprdc,kdsymm,kdolet,kdilet,
     &                           kdsld,idis,kdintr,kdbuff,kdshutr,
     &                           kdpors,kdovst
      use module_io,only       : ifle
!
! 1. Set metric component at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: NALLCV,NCVFAC,NSSFBC,nedge,ICALL,icode
      integer,intent(in)    :: medge,mssfbc,mvrtx
      integer,intent(in)    :: LVEDGE    (2,medge)
      integer,intent(in)    :: LBCSSF    (  mssfbc)
      integer,intent(in)    :: LCYCSF    (  mssfbc)
      real*8 ,intent(inout) :: SFCENT    (3,medge)
      real*8 ,intent(inout) :: SFAREA    (4,medge)
      real*8 ,intent(inout) :: CVCENT    (3,mvrtx)
      real*8 ,intent(inout) :: CVVOLM    (  mvrtx)
!
! --- [local entities] 
!
      real*8  :: aa(3,3,nbcnd)
      real*8  :: r1(3),r2(3),prd,dx,dy,dz,dlvect
      real*8  :: alf,dl,dum1
      real*8,parameter :: SML=1.d-15
      integer :: i,j,k,l,m,n,kdv,kdt,kdy,kdk,kdp
      integer :: ICOM,IMD,ICH,IFLD,IMAT,ICTP
      integer :: IC,IS,IV,IE,ICF,ICV,IBF,NB,KD,IDC
      integer :: ICVA,ICVB,IVA,IVB,IC1,IC2,IBFP,ICFP,ICVP,IDCP
!
!-< 1. Set values at dummy cell >-
!
      aa=0.d0
      call set_rotmtrx(nbcnd,aa)
!
      DO 1000 IBF=1,nssfbc
      nb=LBCSSF(IBF)
      if(nb.eq.0) goto 1000
      kd=kdbcnd(0,abs(nb))
      if(nb.lt.0.and.(kd==kdprdc.and.idis(abs(nb))==0)) goto 1000
!      if(nb.lt.0.and.(kd==kdsld.and.idis(abs(nb))==0)) goto 1000
      nb=abs(nb)
      kd=kdbcnd(0,nb)
      ICF=IBF+nedge
      ICV=LVEDGE(1,ICF)
      IDC=LVEDGE(2,ICF)
!
!--< 1.1 periodic boundary >--
!
      if(kd==kdprdc) then
        if(idis(nb)==0) then
          IBFP=LCYCSF(IBF)
          ICFP=nedge+IBFP
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          CVVOLM(IDC)=CVVOLM(ICVP)
          CVVOLM(IDCP)=CVVOLM(ICV)
          do 101 i=1,3
          r1(i)=(CVCENT(i,ICV)-SFCENT(i,ICF))
          r2(i)=(CVCENT(i,ICVP)-SFCENT(i,ICFP))
  101     continue
          do 102 i=1,3
          CVCENT(i,IDC)=SFCENT(i,ICF)
     &                +(aa(i,1,nb)*r2(1)
     &                 +aa(i,2,nb)*r2(2)
     &                 +aa(i,3,nb)*r2(3))
          CVCENT(i,IDCP)=SFCENT(i,ICFP)
     &                +(aa(1,i,nb)*r1(1)
     &                 +aa(2,i,nb)*r1(2)
     &                 +aa(3,i,nb)*r1(3))
  102     continue
        elseif(idis(nb)==1) then
          IBFP=LCYCSF(IBF)
          ICFP=nedge+IBFP
          ICVP=LVEDGE(1,ICFP)
          CVVOLM(IDC)=CVVOLM(ICV)
          do i=1,3
          r1(i)=(CVCENT(i,ICV)-SFCENT(i,ICF))
          enddo
          prd=r1(1)*SFAREA(1,ICF)
     &       +r1(2)*SFAREA(2,ICF)
     &       +r1(3)*SFAREA(3,ICF)
          prd=prd+prd

          do i=1,3
          CVCENT(i,IDC)=SFCENT(i,ICF)+(r1(i)-prd*SFAREA(i,ICF))
          enddo
        endif
      elseif(kd==kdovst) then
!         IBFP=LCYCSF(IBF)
!         ICFP=nedge+IBFP
!         ICVP=LVEDGE(1,ICFP)
!         r1(:)=SFCENT(:,ICFP)-CVCENT(:,ICVP) 
!
        r1(:)=SFCENT(:,ICF)-CVCENT(:,ICV)
        
        CVVOLM(IDC)=CVVOLM(ICV)        
        prd=dsqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
        do i=1,3
          CVCENT(i,IDC)=SFCENT(i,ICF)+prd*SFAREA(i,ICF)
        enddo
      elseif(kd==kdsld) then
!      elseif(kd==kdsld.or.kd==kdovst) then
!
! --- sliding BC zhang8
!
!        if(idis(nb)==0) then
!          IBFP=LCYCSF(IBF)
!          ICFP=nedge+IBFP
!          ICVP=LVEDGE(1,ICFP)
!          IDCP=LVEDGE(2,ICFP)
!          CVVOLM(IDC) =CVVOLM(ICVP)
!          CVVOLM(IDCP)=CVVOLM(ICV)
!          do i=1,3
!          CVCENT(i,IDC) =SFCENT(i,ICF) !  CVCENT(i,ICVP)
!          CVCENT(i,IDCP)=SFCENT(i,ICFP)!  CVCENT(i,ICV)
!          enddo
!        elseif(idis(nb)==1) then 
          IBFP=LCYCSF(IBF)
          ICFP=nedge+IBFP
          ICVP=LVEDGE(1,ICFP)
          CVVOLM(IDC) =CVVOLM(ICV)
!          do i=1,3
!          r1(i)=(CVCENT(i,ICV)-SFCENT(i,ICF))
!          enddo
          r1(:)=SFCENT(:,ICFP)-CVCENT(:,ICVP) 
!          prd=r1(1)*SFAREA(1,ICF)
!     &       +r1(2)*SFAREA(2,ICF)
!     &       +r1(3)*SFAREA(3,ICF)
!          prd=prd+prd
          prd=dsqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
          do i=1,3
!          CVCENT(i,IDC)=SFCENT(i,ICF)+(r1(i)-prd*SFAREA(i,ICF))
          CVCENT(i,IDC)=SFCENT(i,ICF)+prd*SFAREA(i,ICF)
          enddo
!        endif
!
!--< 1.2 symmetric boundary >--
!
      elseif(kd.eq.kdsymm) then
        CVVOLM(IDC)=CVVOLM(ICV)
        do 201 i=1,3
        r1(i)=CVCENT(i,ICV)-SFCENT(i,ICF)
  201   continue
        prd=r1(1)*SFAREA(1,ICF)
     &     +r1(2)*SFAREA(2,ICF)
     &     +r1(3)*SFAREA(3,ICF)
        prd=prd+prd
        do 202 i=1,3
        CVCENT(i,IDC)=SFCENT(i,ICF)+(r1(i)-prd*SFAREA(i,ICF))
  202   continue
      elseif(kd.eq.kdolet) then
!
!--< 1.3 outlet boundary >--
!
        CVVOLM(IDC)=0.d0
        do 302 i=1,3
        CVCENT(i,IDC)=SFCENT(i,ICF)
  302   continue
      elseif(kd==kdbuff.or.kd==kdshutr.or.kd==kdpors) then
        IBFP=LCYCSF(IBF)
        if(IBFP==0) cycle
        ICFP=nedge+IBFP
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        CVVOLM(IDC)=CVVOLM(ICVP)
        CVVOLM(IDCP)=CVVOLM(ICV)
        do i=1,3
        CVCENT(i,IDC)=CVCENT(i,ICVP)
        CVCENT(i,IDCP)=CVCENT(i,ICV)
        enddo
      elseif(kd==kdintr) then
        if(idis(nb)==0) then
          IBFP=LCYCSF(IBF)
          if(IBFP==0) cycle
          ICFP=nedge+IBFP
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          CVVOLM(IDC)=CVVOLM(ICVP)
          CVVOLM(IDCP)=CVVOLM(ICV)
          do i=1,3
          CVCENT(i,IDC)=CVCENT(i,ICVP)
          CVCENT(i,IDCP)=CVCENT(i,ICV)
          enddo
        elseif(idis(nb)==1) then
          IBFP=LCYCSF(IBF)
          ICFP=nedge+IBFP
          ICVP=LVEDGE(1,ICFP)
          CVVOLM(IDC) =CVVOLM(ICV)
          r1(:)=SFCENT(:,ICFP)-CVCENT(:,ICVP) 
          prd=dsqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
          do i=1,3
          CVCENT(i,IDC)=SFCENT(i,ICF)+prd*SFAREA(i,ICF)
          enddo
        endif
!       
!--< 1.3 other boundary (incluing "interface")>-- 
!
      else
        CVVOLM(IDC)=0.d0
        do 301 i=1,3
        CVCENT(i,IDC)=SFCENT(i,ICF)
  301   continue
      endif
!
 1000 continue
!
! --- 
!
      do 2000 IBF=1,nssfbc   !ICF=1,NCVFAC
      ICF=IBF+nedge
      nb=LBCSSF(IBF)
      if(nb.eq.0) goto 2000
      ICVA=LVEDGE(1,ICF)
      ICVB=LVEDGE(2,ICF)
      if(CVCENT(1,ICVB).eq.0.d0.and.
     &   CVCENT(2,ICVB).eq.0.d0.and.
     &   CVCENT(3,ICVB).eq.0.d0) then
        write(ifle,*) 
     &   ' ### ERR: Dummy Cell coordinates are NOT defined'
        write(ifle,*) ' Please Contact FrontFlow developer',nb
        stop 'stop at dc_metric'
      endif
 2000 continue
!---------------------
! --- check all CV
!---------------------
      ICFP=0
      do 3000 ICF=1,NCVFAC
      if(ICF.gt.NEDGE) goto 3000
      if(SFCENT(1,ICF).eq.0.d0.and.
     &   SFCENT(2,ICF).eq.0.d0.and.
     &   SFCENT(3,ICF).eq.0.d0) then
        ICFP=ICFP+1
      endif
 3000 continue
      if(ICFP.ge.2) then
        write(ifle,*) 
     &  ' ### ERR: CV Sub Face coordinates are NOT defined'
        write(ifle,*) ' Please Contact FrontFlow developer'
        stop 'stop at dc_metric'
      endif
!
      ICFP=0
      do 3010 ICF=1,NCVFAC
      if(abs(SFAREA(4,ICF)).le.SML) then
        ICFP=ICFP+1
      endif
 3010 continue
      if(ICFP.gt.0) then
        write(ifle,*) 
     &   ' ### ERR: CV Sub Face Area is ZERO N=',ICFP
        write(ifle,*) ' Please Contact FrontFlow developer'
        stop 'stop at dc_metric'
      endif
!
! --- 
!
      do 4000 ICF=1,NCVFAC
      ICVA=LVEDGE(1,ICF)
      ICVB=LVEDGE(2,ICF)
      if(ICVB.lt.ICVA)then
        write(ifle,*) ' ### ERR: ICVB < ICVA'
        write(ifle,*) ' Please contact FrontFlow Developer'
        stop 'stop at dc_metric'
      endif
      dx=CVCENT(1,ICVB)-CVCENT(1,ICVA)
      dy=CVCENT(2,ICVB)-CVCENT(2,ICVA)
      dz=CVCENT(3,ICVB)-CVCENT(3,ICVA)
      alf=(dx*SFAREA(1,ICF)
     &    +dy*SFAREA(2,ICF)
     &    +dz*SFAREA(3,ICF))
!      if(alf.le.0.d0) then
!        write(ifle,'(a)')
!     &  'WRN: Opposite Dir. between Edge vector and SF vector'
!        write(ifle,'(5x,a,I4,a,I10,E10.4)')
!     &   'ICALL=',ICALL,' At SF=',ICF,alf
!        write(ifle,'(5x,a,I10,a,I10)') 
!     &   'SF between ICVA=',ICVA,'and ICVB=',ICVB
!      endif
      if(abs(alf)<1.d-10) then
        write(ifle,'(1X,a)') 'ERROR: Center of cell is too near'
        write(ifle,'(1X,a,2x,I8,4E16.8)') 
     &              'MSG: ICF= ',ICF,SFAREA(1:4,ICF)
        write(ifle,'(1X,a,I8,3E16.8)') 
     &              'MSG: ICVA=, Coord.=',ICVA,CVCENT(:,ICVA)
        write(ifle,'(1X,a,I8,3E16.8)')  
     &       'MSG: ICVB=, Coord.=',ICVB,CVCENT(:,ICVB)
        
        write(ifle,'(1X,a)') 'MSG: 2D MESH is NOT used by FFR'
        write(ifle,'(1X,a,F16.4)') 'Distance= ',abs(alf)
        stop 'ERR: Mesh too bad'
      endif
 4000 continue
!
      return
!
      end subroutine dc_metric
!
