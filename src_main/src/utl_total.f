!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      
!      subroutine utl_bcgstb
!      subroutine utl_bcgstb_MAT
!
!      subroutine utl_iccg_hpc      =Xa
!      subroutine utl_iccg_hpc_MAT  =Xa
!      subroutine utl_iccg          =Xa
!      subroutine utl_iccg_vect
!      subroutine utl_bcgstb_vect
!      subroutine utl_iccg_hpc_vect
!      
!      subroutine utl_gausel(mko,nko,n,ip,a,b,x)
!      subroutine utl_gausel_nop(mko,nko,n,a,b,x)
!      subroutine utl_getpvt(mko,nko,n,a,ip,ierr)
!      subroutine utl_getpvt_nop(mko,nko,n,a,ierr)
!      subroutine utl_random(r)
!      
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
      subroutine utl_bcgstb(
     &   cn,m,MXALLCV,MXCV,MAXIE1,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,ip,aa,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aeps,reps,itr,ier,IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
!
!  1. Bi-CGstab solver with diagonal scaling
!     & incomplete LU decomposition
!
!     nn   : no. of elements
!     MAXIE   : maximul no. of off-diagonal elements
!     IQ   : no. of off-diagonal elements
!     ip   : address of off-diagonal element
!     aa   : coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence  (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence  (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!     MAT_INDEX=>NCVIN
!     MAT_CVEXT=>NCV
!     MAT_DCIDX=NCV~NALLCV
!
! --- [module arguments]
!
      use module_hpcutil
      use module_cgsolver,only : nostop
      use module_io,only       : ifli,ifll,ifle
      use module_metrix,only   : D  =>W1K1
      use module_metrix,only   : p0 =>W1K2
      use module_metrix,only   : p1 =>W1K3
      use module_metrix,only   : q0 =>W1K4
      use module_metrix,only   : q1 =>W1K5
      use module_metrix,only   : r0 =>W1K6
!
      use module_metrix,only   : r1 =>W1K7
      use module_metrix,only   : r2 =>W1K8
      use module_metrix,only   : x  =>W1K9
      use module_metrix,only   : z  =>W1K10
      use module_metrix,only   : a  =>W2K1
!      use module_material,only : KMAT_S,MAT_ALL,MAT_S
      use module_metrix,only   : t,s,u,v,alpha,beta,omega
      use module_model,only    : ical_vect,nthrds
      use module_dimension,only : MXCOMP
!
      implicit none
!
! --- [dummy arguments]
!
      character,intent(in)  :: cn
      integer,intent(in)    :: m,MXALLCV,MXCV,MXMAT,NMAT,IQMXCV
      integer,intent(in)    :: NALLCV,NCV,NCVIN,MAXIE,IMODE,MAXIE1
      integer,intent(in)    :: IQ(IQMXCV,2)
      integer,intent(in)    :: ip(MXALLCV,MAXIE)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: aa(MXALLCV,0:MAXIE)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
      real*8 ,parameter  :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0
      real*8 :: c1,dum1,dum2,rmaxs,rmaxs0
!
      integer:: maxitr=300,IITER=0
      integer:: i,j,k,l,it1,it2,jer,m1,ICV,ierr,KMAT,idum
      integer:: IMAT,IIMAT,ICVS,ICVE,ICVL
!
      ierr=0
      IITER=0
!
      if(.not.ical_vect) then
        allocate(a(MXCV,0:MAXIE),stat=ierr)
      endif
!
      t(:)=0.d0
      s(:)=0.d0
      u(:)=0.d0
      v(:)=0.d0
      alpha(:)=0.d0
      beta(:)=0.d0
      omega(:)=0.d0
      a=0.d0
!
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in utl_bcgstb'
        call FFRABORT(1,'utl_bcgstb')
      endif
!
!-< 1. Preliminary set >-
!
!
      rmaxs=1.d0
      do 100 IIMAT=1,NMAT    !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) goto 100
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      DO ICVL=ICVS,ICVE
      rmaxs=min(rmaxs,abs(aa(ICVL,0)))
      enddo
  100 continue
!
      IF(NPE.gt.1) THEN
        CALL hpcrmin(RMAXS)
      ENDIF
!
      if(RMAXS.lt.0.d0) then
        write(ifle,'(a,E14.4,3X,a,3X,I4)') 
     &  '### interface error -1- (utl_bcgstb)',rmaxs,cn,m
        CALL FFRABORT(1,'utl_bcgstb')
      elseif(RMAXS.eq.0.d0) then
        deallocate(a)
        return
      endif
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
!
      d (:)=0.d0
      p0(:)=0.d0
      p1(:)=0.d0
      q0(:)=0.d0
      q1(:)=0.d0
      r0(:)=0.d0
      r1(:)=0.d0
      r2(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0
!
!-< 2. Preconditioning >-
!
!
!--< 2.1 diagonal scaling >--
!
      do 2000 IIMAT=1,NMAT    !
      if(.not.mat_cal(IIMAT)) goto 2000
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do 200 ICVL=ICVS,ICVE   !ICV=1,NCV
      d(ICVL)=1.d0/(dsqrt(abs(aa(ICVL,0))))
      z(ICVL)=d(ICVL)
      a(ICVL,0)=sign(1.d0,aa(ICVL,0))
  200 continue
 2000 continue
!
! --- Low and Upper region
!
      rmaxs=0.d0
      do 2500 IIMAT=1,NMAT    !
      if(.not.mat_cal(IIMAT)) goto 2500 
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do 210 ICVL=ICVS,ICVE   !ICV=1,NCV
      do j=1,IQ(ICVL,2)
      a(ICVL,j)=aa(ICVL,j)*d(ICVL)*d(ip(ICVL,j))*a(ICVL,0)
      enddo
  210 continue
!
      dum1=-1.d0
      do 220 ICVL=ICVS,ICVE    !ICV=1,NCV
      bb(ICVL)=bb(ICVL)*d(ICVL)*a(ICVL,0)
      x(ICVL)=0.d0
      rmaxs=max(rmaxs,abs(bb(ICVL)))
  220 continue
!
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      d(ICVS:ICVE)=1.d0
      a(ICVS:ICVE,0)=0.d0
 2500 continue
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
      rmaxs0=RMAXS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- one iteration at least

      if(cn=='y'.or.IMODE==7.or.cn=='V') then
      else
        if(RMAXS.le.aeps) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
! --- all below by NCVIN
!
      do 3000 IIMAT=1,NMAT    !
      if(.not.mat_cal(IIMAT)) goto 3000
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
!
!--< 2.2 incomplete LU decomposition >--
!
      dum1=0.d0
      do 230 ICVL=ICVS,ICVE    !ICV=1,NCVIN
      if(d(ICVL).le.0.d0) then
        dum1=d(ICVL)
      endif
      d(ICVL)=1.d0/(d(ICVL))
! --- upper region
      do 231 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      k=ip(ICVL,j)
      c1=0.d0
! --- low
      do 232 l=1,IQ(k,1)
      if(ip(k,l).eq.ICVL) c1=d(ICVL)*a(k,l)
  232 continue
      d(k)=d(k)-c1*a(ICVL,j)
  231 continue
  230 continue
!
 3000 continue
!
      if(NPE.gt.1) then
        call hpcrmin(dum1)
      endif
!
      if(dum1.lt.0.d0) then
        write(ifle,'(3X,A,I4)') 
     &   '### CPU NUMBER MY_RANK= ',MY_RANK
        write(ifle,'(a,E14.4,3X,a,3X,I4)')
     &  '### interface error -2- (utl_bcgstb)',dum1,cn,m
        CALL FFRABORT(1,'utl_bcgstb')
      endif
!
! --- 
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,D)
      endif
!
!-< 3. Initial set for iteration >-
!
  888 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X)
      endif
!
! --- 
!
      rmaxs=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 300 ICVL=ICVS,ICVE
      dum1=x(ICVL)
! --- upper and low region:
      do 301 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*x(ip(ICVL,j))
  301 continue
!
      r0(ICVL)=bb(ICVL)-dum1
      rmaxs=max(rmaxs,abs(r0(ICVL)))
  300 continue
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if(cn=='y') then
        if(MXCOMP>1.and.IITER==0) then
        else
          if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
        endif
      elseif(cn=='V') then
        if(itr>1) then
          if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
        endif
      else
        if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      s(:)=0.d0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 310 ICVL=ICVS,ICVE     !i=1,NCVIN
        p0(ICVL)=r0(ICVL)
        r1(ICVL)=r0(ICVL)
        dum1=sign(1.d0,r0(ICVL))*max(1.d0,abs(r0(ICVL)))
        s(IIMAT)=s(IIMAT)+dum1*r0(ICVL)
        r0(ICVL)=dum1
  310 continue
      enddo
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(s(IIMAT))
!          s(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(s(IIMAT))
        enddo
      endif
      c1=0.d0
      do IIMAT=1,NMAT
        c1=c1+abs(s(IIMAT))
      enddo
      if(NPE.gt.1) then
        CALL hpcrsum(c1)
      endif
!
      if(c1.lt.1.d-50) jer=1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      if( jer.gt.0 ) then
        do IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        do 320 ICVL=ICVS,ICVE     !i=1,NCVIN
          dum1=0.d0
! --- low region:
          do 321 j=1,IQ(ICVL,1)
!            k=ip(ICVL,j)
!            if(k>ICVE.or.k<ICVS) cycle   !zhang???
            dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  321     continue
          p1(ICVL)=d(ICVL)*(r1(ICVL)-dum1)
  320   continue
        enddo
!
        IF(NPE.gt.1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
        endif
!
        do IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        do 322 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper region: 
        do 323 j=IQ(ICVL,1)+1,IQ(ICVL,2)
!            k=ip(ICVL,j)
!            if(k>ICVE.or.k<ICVS) cycle   !zhang???
        dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  323   continue
        p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
        x(ICVL)=x(ICVL)+p1(ICVL)
  322   continue
        enddo
!
!
        jer=0
        it1=it1+1
        itr=it1
        IITER=1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. Iteration >-
!
! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      do 1000 itr=it1,it2
!
        IF(NPE.gt.1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
        endif
!------------------------------------------------------1
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 400 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=0.d0
! --- low 
      do 401 j=1,IQ(ICVL,1)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  401 enddo
      p1(ICVL)=d(ICVL)*(p0(ICVL)-dum1)
  400 enddo
      enddo
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif
!
! --- 
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 402 ICVL=ICVE,ICVS,-1    !i=NCVIN,1,-1
      dum1=0.d0
! --- upper
      do 403 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  403 continue
      p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
  402 enddo
      enddo
!------------------------------------------------------1
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif
!
      v(:)=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 410 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=p1(ICVL)
! --- low + upper
      do 411 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  411 enddo
      q1(ICVL)=dum1
      v(IIMAT)=v(IIMAT)+r0(ICVL)*dum1
  410 enddo
      enddo
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(v(IIMAT))
!          v(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(v(IIMAT))
        enddo
      endif      
!
      c1=0.0d0
      do IIMAT=1,NMAT
      c1=c1+abs(v(IIMAT))
      enddo
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c1)
      ENDIF
      if(c1==0.d0) then
        jer=1
        it1=itr
        IITER=1
        goto 888
      endif
!
      do IIMAT=1,NMAT
      if(v(IIMAT)==0.d0)  v(IIMAT)=SML
      alpha(IIMAT)=s(IIMAT)/(v(IIMAT)) 
      enddo
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,R1)
      endif
!---------------------------------------------------2
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 420 ICVL=ICVS,ICVE    !i=1,NCVIN
        q0(ICVL)=r1(ICVL)-alpha(IIMAT)*q1(ICVL)
  420 continue
!      enddo
!
!      do IIMAT=1,NMAT
!      IMAT=MAT_NO(IIMAT)
!      ICVS=MAT_CVEXT(IIMAT-1)+1
!      ICVE=MAT_INDEX(IIMAT)
!      if(.not.mat_cal(IIMAT)) cycle
      do 430 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum2=0.d0
! --- low
        do 431 j=1,IQ(ICVL,1)
          dum2=dum2+a(ICVL,j)*r1(ip(ICVL,j))
  431   continue
        r1(ICVL)=d(ICVL)*(q0(ICVL)-dum2)
  430 continue
      enddo
!
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,R1)
      endif
!
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 432 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper
        do 433 j=IQ(ICVL,1)+1,IQ(ICVL,2)
        dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  433   continue
        r1(ICVL)=r1(ICVL)-d(ICVL)*dum1
  432 continue
      enddo
!---------------------------------------------------2
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,R1)
      endif
!
      u(:)=0.d0
      v(:)=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 440 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=r1(ICVL)
        do 441 j=1,IQ(ICVL,2)
! --- low+upper
        dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  441   continue
        r2(ICVL)=dum1
        u(IIMAT)=u(IIMAT)+dum1*q0(ICVL)
        v(IIMAT)=v(IIMAT)+dum1*dum1
  440 continue
      enddo
!--------------------------------------------3
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(u(IIMAT))
!          CALL hpcrsum(v(IIMAT))
!          u(0)=0.d0
!          v(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(u(IIMAT))
        CALL hpcrsum(v(IIMAT))
        enddo
      endif
      c1=0.d0
      do IIMAT=1,NMAT
        c1=c1+abs(v(IIMAT))
      enddo
      if(NPE.gt.1) then
        CALL hpcrsum(c1)
      endif
!
      if(c1==0.d0) then
        jer=1
        it1=itr
        IITER=1
        goto 888
      endif
!
      do IIMAT=1,NMAT
      if(v(IIMAT)==0.d0)  v(IIMAT)=SML
      omega(IIMAT)=u(IIMAT)/(v(IIMAT))
      enddo
!
      rmaxs=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 450 ICVL=ICVS,ICVE   !i=1,NCVIN
        x(ICVL)=x(ICVL)+(alpha(IIMAT)*p1(ICVL)+omega(IIMAT)*r1(ICVL))
        r1(ICVL)=q0(ICVL)-omega(IIMAT)*r2(ICVL)
        rmaxs=max(rmaxs,abs(r1(ICVL)))
  450 enddo
      enddo
!--------------------------------------------
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!
      if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
!
      t(:)=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 460 ICVL=ICVS,ICVE   !i=1,NCVIN
        t(IIMAT)=t(IIMAT)+r0(ICVL)*r1(ICVL)
  460 continue
      enddo
!--------------------------------------------
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(t(IIMAT))
!          t(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(t(IIMAT))
        enddo
      ENDIF
!
      c1=0.d0
      do IIMAT=1,NMAT
      c1=c1+abs(s(IIMAT)*omega(IIMAT))
      enddo
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c1)
      ENDIF
      if(c1==0.d0) then
        jer=1
        it1=itr
        IITER=1
        goto 888
      endif
!
      do IIMAT=1,NMAT
      beta(IIMAT)=t(IIMAT)*alpha(IIMAT)/(SML+s(IIMAT)*omega(IIMAT))
      enddo
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 470 ICVL=ICVS,ICVE    !i=1,NCVIN
      p0(ICVL)=r1(ICVL)+beta(IIMAT)*(p0(ICVL)-omega(IIMAT)*q1(ICVL))
  470 continue
      enddo
      s(:)=t(:)
!
 1000 continue
!
      if(nostop==1) then
        if(my_rank==root) then
           write(ifle,'(1X,3a,I4,a)')
     &    ' ### WRN1: BiCGSTAB Solver is NOT converged: ',trim(cn),
     & ' Bicg_Iter= ', maxitr,' NOT stop and go into next time step'
        endif
        goto 999
      endif
!
      ier=1
!
  999 continue
!
      aeps=max(aeps,rmaxs)
      if(rmaxs0.lt.SML) then
        reps=max(reps,1.d0)
      else
        reps=max(reps,rmaxs/(rmaxs0+SML))
      endif
!
      do 650 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(mat_cal(IIMAT)) then
        do 500 ICVL=ICVS,ICVE
        bb(ICVL)=x(ICVL)*z(ICVL)
  500   continue
      endif
 650  continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      if(.not.ical_vect) then
        deallocate(a)
      endif
!
      return
      end subroutine utl_bcgstb
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_bcgstb1(
     &   cn,m,MXALLCV,MXCV,MAXIE1,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,ip,aa,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aeps,reps,itr,ier,IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. Bi-CGstab solver with diagonal scaling
!     & incomplete LU decomposition
!
!     nn   : no. of elements
!     MAXIE   : maximul no. of off-diagonal elements
!     IQ   : no. of off-diagonal elements
!     ip   : address of off-diagonal element
!     aa   : coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence  (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence  (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!     MAT_INDEX=>NCVIN
!     MAT_CVEXT=>NCV
!     MAT_DCIDX=NCV~NALLCV
!
! --- [module arguments]
!
      use module_hpcutil
      use module_cgsolver,only : nostop
      use module_io,only       : ifli,ifll,ifle
      use module_metrix,only   : D  =>W1K1
      use module_metrix,only   : p0 =>W1K2
      use module_metrix,only   : p1 =>W1K3
      use module_metrix,only   : q0 =>W1K4
      use module_metrix,only   : q1 =>W1K5
      use module_metrix,only   : r0 =>W1K6
!
      use module_metrix,only   : r1 =>W1K7
      use module_metrix,only   : r2 =>W1K8
      use module_metrix,only   : x  =>W1K9
      use module_metrix,only   : z  =>W1K10
      use module_metrix,only   : a  =>W2K1
!      use module_material,only : KMAT_S,MAT_ALL,MAT_S
      use module_metrix,only   : t,s,u,v,alpha,beta,omega
      use module_model,only    : ical_vect,nthrds
      use module_dimension,only : MXCOMP
!
      implicit none
!
! --- [dummy arguments]
!
      character,intent(in)  :: cn
      integer,intent(in)    :: m,MXALLCV,MXCV,MXMAT,NMAT,IQMXCV
      integer,intent(in)    :: NALLCV,NCV,NCVIN,MAXIE,IMODE,MAXIE1
      integer,intent(in)    :: IQ(IQMXCV,2)
      integer,intent(in)    :: ip(MXALLCV,MAXIE)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: aa(MXALLCV,0:MAXIE)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
      real*8 ,parameter  :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0
      real*8 :: c1,dum1,dum2,rmaxs,rmaxs0
!
      integer:: maxitr=300
      integer:: i,j,k,l,it1,it2,jer,m1,ICV,ierr,KMAT,idum
      integer:: IMAT,IIMAT,ICVS,ICVE,ICVL
!
!
      ierr=0
!
      if(.not.ical_vect) then
        allocate(a(MXCV,0:MAXIE),stat=ierr)
      endif
!
      t(:)=0.d0
      s(:)=0.d0
      u(:)=0.d0
      v(:)=0.d0
      alpha(:)=0.d0
      beta(:)=0.d0
      omega(:)=0.d0
      a=0.d0
!
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in utl_bcgstb'
        call FFRABORT(1,'utl_bcgstb')
      endif
!
!-< 1. Preliminary set >-
!
      rmaxs=1.d0
      do 100 IIMAT=1,NMAT    !ICV=1,NCV
      if(.not.mat_cal(IIMAT)) goto 100
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      DO ICVL=ICVS,ICVE
      rmaxs=min(rmaxs,abs(aa(ICVL,0)))
      enddo
  100 continue
!
      IF(NPE.gt.1) THEN
        CALL hpcrmin(RMAXS)
      ENDIF
!
      if(RMAXS.lt.0.d0) then
        write(ifle,'(a,E14.4,3X,a,3X,I4)') 
     &  '### interface error -1- (utl_bcgstb)',rmaxs,cn,m
        CALL FFRABORT(1,'utl_bcgstb')
      elseif(RMAXS.eq.0.d0) then
        deallocate(a)
        return
      endif
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
!
      d (:)=0.d0
      p0(:)=0.d0
      p1(:)=0.d0
      q0(:)=0.d0
      q1(:)=0.d0
      r0(:)=0.d0
      r1(:)=0.d0
      r2(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0
!
!-< 2. Preconditioning >-
!
!
!--< 2.1 diagonal scaling >--
!
      do 2000 IIMAT=1,NMAT    !
      if(.not.mat_cal(IIMAT)) goto 2000
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      do 200 ICVL=ICVS,ICVE   !ICV=1,NCV
      d(ICVL)=1.d0/(dsqrt(abs(aa(ICVL,0))))
      z(ICVL)=d(ICVL)
      a(ICVL,0)=sign(1.d0,aa(ICVL,0))
  200 continue
 2000 continue
!
! --- Low and Upper region
!
      rmaxs=0.d0
      do 2500 IIMAT=1,NMAT    !
      if(.not.mat_cal(IIMAT)) goto 2500
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do 210 ICVL=ICVS,ICVE   !ICV=1,NCV
      do j=1,IQ(ICVL,2)
      a(ICVL,j)=aa(ICVL,j)*d(ICVL)*d(ip(ICVL,j))*a(ICVL,0)
      enddo
  210 continue
!
      do 220 ICVL=ICVS,ICVE    !ICV=1,NCV
      bb(ICVL)=bb(ICVL)*d(ICVL)*a(ICVL,0)
      x(ICVL)=0.d0
      rmaxs=max(rmaxs,abs(bb(ICVL)))
  220 continue
!
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      d(ICVS:ICVE)=1.d0
      a(ICVS:ICVE,0)=0.d0
 2500 continue
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- one iteration at least
      if(cn=='y') then
      else
        if(RMAXS.le.aeps) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rmaxs0=RMAXS
!
! --- all below by NCVIN
!
      do 3000 IIMAT=1,NMAT    !
      if(.not.mat_cal(IIMAT)) goto 3000
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
!
!--< 2.2 incomplete LU decomposition >--
!
      dum1=0.d0
      do 230 ICVL=ICVS,ICVE    !ICV=1,NCVIN
      if(d(ICVL).le.0.d0) then
        dum1=d(ICVL)
      endif
      d(ICVL)=1.d0/(d(ICVL))
! --- upper region
      do 231 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      k=ip(ICVL,j)
      c1=0.d0
! --- low
      do 232 l=1,IQ(k,1)
      if(ip(k,l).eq.ICVL) c1=d(ICVL)*a(k,l)
  232 continue
      d(k)=d(k)-c1*a(ICVL,j)
  231 continue
  230 continue
!
 3000 continue
!
      if(NPE.gt.1) then
        call hpcrmin(dum1)
      endif
!
      if(dum1.lt.0.d0) then
        write(ifle,'(3X,A,I4)') 
     &  '### CPU NUMBER MY_RANK= ',MY_RANK
        write(ifle,'(a,E14.4,3X,a,3X,I4)')
     &  '### interface error -2- (utl_bcgstb)',dum1,cn,m
        CALL FFRABORT(1,'utl_bcgstb')
      endif
!
! --- 
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,D)
      endif
!
!-< 3. Initial set for iteration >-
!
  888 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X)
      endif
!
! --- 
!
      rmaxs=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 300 ICVL=ICVS,ICVE
      dum1=x(ICVL)
! --- upper and low region:
      do 301 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*x(ip(ICVL,j))
  301 continue
!
      r0(ICVL)=bb(ICVL)-dum1
      rmaxs=max(rmaxs,abs(r0(ICVL)))
  300 continue
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(cn=='y') then
        if(MXCOMP>1) then
        else
          if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
        endif
      else
        
        if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      s(:)=0.d0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 310 ICVL=ICVS,ICVE     !i=1,NCVIN
        p0(ICVL)=r0(ICVL)
        r1(ICVL)=r0(ICVL)
        dum1=sign(1.d0,r0(ICVL))*max(1.d0,abs(r0(ICVL)))
        s(IIMAT)=s(IIMAT)+dum1*r0(ICVL)
        r0(ICVL)=dum1
  310 continue
      enddo
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(s(IIMAT))
!          s(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(s(IIMAT))
        enddo
      endif
      c1=0.d0
      do IIMAT=1,NMAT
        c1=c1+abs(s(IIMAT))
      enddo
      if(NPE.gt.1) then
        CALL hpcrsum(c1)
      endif
!
      if(c1.lt.1.d-50) jer=1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      if( jer.gt.0 ) then
        do IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        do 320 ICVL=ICVS,ICVE     !i=1,NCVIN
          dum1=0.d0
! --- low region:
          do 321 j=1,IQ(ICVL,1)
!            k=ip(ICVL,j)
!            if(k>ICVE.or.k<ICVS) cycle   !zhang???
            dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  321     continue
          p1(ICVL)=d(ICVL)*(r1(ICVL)-dum1)
  320   continue
        enddo
!
        IF(NPE.gt.1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
        endif
!
        do IIMAT=1,NMAT
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        if(.not.mat_cal(IIMAT)) cycle
        do 322 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper region:
        do 323 j=IQ(ICVL,1)+1,IQ(ICVL,2)
!            k=ip(ICVL,j)
!            if(k>ICVE.or.k<ICVS) cycle   !zhang???
        dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  323   continue
        p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
        x(ICVL)=x(ICVL)+p1(ICVL)
  322   continue
        enddo
!
!
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. Iteration >-
!
! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      do 1000 itr=it1,it2
!
        IF(NPE.gt.1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
        endif
!------------------------------------------------------1
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 400 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=0.d0
! --- low
      do 401 j=1,IQ(ICVL,1)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  401 enddo
      p1(ICVL)=d(ICVL)*(p0(ICVL)-dum1)
  400 enddo
      enddo
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif
!
! --- 
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 402 ICVL=ICVE,ICVS,-1    !i=NCVIN,1,-1
      dum1=0.d0
! --- upper
      do 403 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  403 continue
      p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
  402 enddo
      enddo
!------------------------------------------------------1
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif
!
      v(:)=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 410 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=p1(ICVL)
! --- low + upper
      do 411 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  411 enddo
      q1(ICVL)=dum1
      v(IIMAT)=v(IIMAT)+r0(ICVL)*dum1
  410 enddo
      enddo
!
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(v(IIMAT))
!          v(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(v(IIMAT))
        enddo
      endif      
!
      c1=0.0d0
      do IIMAT=1,NMAT
      c1=c1+abs(v(IIMAT))
      enddo
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c1)
      ENDIF
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      do IIMAT=1,NMAT
      if(v(IIMAT)==0.d0)  v(IIMAT)=SML
      alpha(IIMAT)=s(IIMAT)/(v(IIMAT)) 
      enddo
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,R1)
      endif
!---------------------------------------------------2
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 420 ICVL=ICVS,ICVE    !i=1,NCVIN
        q0(ICVL)=r1(ICVL)-alpha(IIMAT)*q1(ICVL)
  420 continue
!      enddo
!
!      do IIMAT=1,NMAT
!      IMAT=MAT_NO(IIMAT)
!      ICVS=MAT_CVEXT(IIMAT-1)+1
!      ICVE=MAT_INDEX(IIMAT)
!      if(.not.mat_cal(IIMAT)) cycle
      do 430 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum2=0.d0
! --- low
        do 431 j=1,IQ(ICVL,1)
          dum2=dum2+a(ICVL,j)*r1(ip(ICVL,j))
  431   continue
        r1(ICVL)=d(ICVL)*(q0(ICVL)-dum2)
  430 continue
      enddo
!
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,R1)
      endif
!
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 432 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper
        do 433 j=IQ(ICVL,1)+1,IQ(ICVL,2)
          dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  433   continue
        r1(ICVL)=r1(ICVL)-d(ICVL)*dum1
  432 continue
      enddo
!---------------------------------------------------2
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXALLCV,NCV,R1)
      endif
!
      u(:)=0.d0
      v(:)=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 440 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=r1(ICVL)
        do 441 j=1,IQ(ICVL,2)
! --- low+upper
        dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  441   continue
        r2(ICVL)=dum1
        u(IIMAT)=u(IIMAT)+dum1*q0(ICVL)
        v(IIMAT)=v(IIMAT)+dum1*dum1
  440 continue
      enddo
!--------------------------------------------3
      if(NPE.gt.1) then
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(u(IIMAT))
!          CALL hpcrsum(v(IIMAT))
!          u(0)=0.d0
!          v(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(u(IIMAT))
        CALL hpcrsum(v(IIMAT))
        enddo
      endif
      c1=0.d0
      do IIMAT=1,NMAT
        c1=c1+abs(v(IIMAT))
      enddo
      if(NPE.gt.1) then
        CALL hpcrsum(c1)
      endif
!
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      do IIMAT=1,NMAT
      if(v(IIMAT)==0.d0)  v(IIMAT)=SML
      omega(IIMAT)=u(IIMAT)/(v(IIMAT))
      enddo
!
      rmaxs=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 450 ICVL=ICVS,ICVE   !i=1,NCVIN
        x(ICVL)=x(ICVL)+(alpha(IIMAT)*p1(ICVL)+omega(IIMAT)*r1(ICVL))
        r1(ICVL)=q0(ICVL)-omega(IIMAT)*r2(ICVL)
        rmaxs=max(rmaxs,abs(r1(ICVL)))
  450 enddo
      enddo
!--------------------------------------------
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!
      if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
!
      t(:)=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 460 ICVL=ICVS,ICVE   !i=1,NCVIN
        t(IIMAT)=t(IIMAT)+r0(ICVL)*r1(ICVL)
  460 continue
      enddo
!--------------------------------------------
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          CALL hpcrsum(t(IIMAT))
!          t(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrsum(t(IIMAT))
        enddo
      ENDIF
!
      c1=0.d0
      do IIMAT=1,NMAT
      c1=c1+abs(s(IIMAT)*omega(IIMAT))
      enddo
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c1)
      ENDIF
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      do IIMAT=1,NMAT
      beta(IIMAT)=t(IIMAT)*alpha(IIMAT)/(SML+s(IIMAT)*omega(IIMAT))
      enddo
!
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      do 470 ICVL=ICVS,ICVE    !i=1,NCVIN
      p0(ICVL)=r1(ICVL)+beta(IIMAT)*(p0(ICVL)-omega(IIMAT)*q1(ICVL))
  470 continue
      enddo
      s(:)=t(:)
!
 1000 continue
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &    ' ### WRN2: BiCGSTAB Solver is NOT converged: ',
     & 'Bicg_Iter= ', maxitr,' NOT stop and go into next time step'
         endif
         goto 999
      endif
!
      ier=1
!
  999 continue
!
      aeps=max(aeps,rmaxs)
      if(rmaxs0.lt.SML) then
        reps=max(reps,1.d0)
      else
        reps=max(reps,rmaxs/rmaxs0)
      endif
!
      do 650 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(mat_cal(IIMAT)) then
        do 500 ICVL=ICVS,ICVE
        bb(ICVL)=x(ICVL)*z(ICVL)
  500   continue
      endif
 650  continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      if(.not.ical_vect) then
        deallocate(a)
      endif
!
      return
      end subroutine utl_bcgstb1
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_bcgstb_MAT
     &  (cn,MXALLCV,MXCV,MAXIE1,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,ip,aa,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aeps,reps,itr,ier,IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. Bi-CGstab solver with diagonal scaling
!     & incomplete LU decomposition
!
!     nn   : no. of elements
!     MAXIE   : maximul no. of off-diagonal elements
!     IQ   : no. of off-diagonal elements
!     ip   : address of off-diagonal element
!     aa   : coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!     MAT_INDEX=>NCVIN
!     MAT_CVEXT=>NCV
!     MAT_DCIDX=NCV~NALLCV
!
! --- [module arguments]
!
      use module_hpcutil
      use module_metrix,only : D  =>W1K1
      use module_metrix,only : p0 =>W1K2
      use module_metrix,only : p1 =>W1K3
      use module_metrix,only : q0 =>W1K4
      use module_metrix,only : q1 =>W1K5
      use module_metrix,only : r0 =>W1K6
      use module_metrix,only : r1 =>W1K7
      use module_metrix,only : r2 =>W1K8
      use module_metrix,only : x  =>W1K9
      use module_metrix,only : z  =>W1K10
      use module_metrix,only : a  =>W2K1
      use module_io,only       : ifli,ifll,ifle
!      use module_material,only : KMAT_S,MAT_ALL,MAT_S
      use module_cgsolver,only : nostop
      use module_model,only    : ical_vect,nthrds
      use module_dimension,only : MXCOMP
!
      implicit none
!
! --- [dummy arguments]
!
      character,intent(in)  :: cn
      integer,intent(in)    :: MXALLCV,MXCV,MAXIE1,MXMAT,NMAT,IQMXCV
      integer,intent(in)    :: NALLCV,NCV,NCVIN,MAXIE,IMODE
      integer,intent(in)    :: IQ(IQMXCV,2)
      integer,intent(in)    :: ip(MXALLCV,MAXIE)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: aa(MXALLCV,0:MAXIE)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0

      real*8 :: c1,t,s,u,v,rmaxs,rmaxs0
      real*8 :: alpha,beta,omega,dum1
      integer:: maxitr=300
      integer:: i,j,k,l,it1,it2,jer,m1,ICV,ierr,KMAT
      integer:: IMAT,IIMAT,ICVS,ICVE,ICVL
!
!
      t=0.d0
      s=0.d0
      u=0.d0
      v=0.d0
      alpha=0.d0
      beta=0.d0
      omega=0.d0
      ierr=0
!
      if(.not.ical_vect) then
        allocate(a(MXCV,0:MAXIE),
     &          stat=ierr)
      endif
!
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in utl_bcgstb'
        call FFRABORT(1,'utl_bcgstb')
      endif

!
!-< 1. Preliminary set >-
!
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      rmaxs=1.d0
      DO ICVL=ICVS,ICVE
      rmaxs=min(rmaxs,abs(aa(ICVL,0)))
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmin(RMAXS)
      ENDIF
!
      if(RMAXS.lt.0.d0) then
         write(*,'(a,E14.4,a)') 
     &  '### interface error -1- (utl_bcgstb)',rmaxs,cn
        CALL FFRABORT(1,'utl_bcgstb')
      elseif(RMAXS.eq.0.d0) then
        deallocate(a)
        return
      endif
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
!
      d (:)=0.d0
      p0(:)=0.d0
      p1(:)=0.d0
      q0(:)=0.d0
      q1(:)=0.d0
      r0(:)=0.d0
      r1(:)=0.d0
      r2(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0
!-----------------------------
!-< 2. Preconditioning >-
!-----------------------------
!-----------------------------
!--< 2.1 diagonal scaling >--
!-----------------------------
!---------------------------------------------------------------
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!---------------------------------------------------------------
      do 200 ICVL=ICVS,ICVE
      d(ICVL)=1.d0/(dsqrt(abs(aa(ICVL,0))))
      z(ICVL)=d(ICVL)
      a(ICVL,0)=sign(1.d0,aa(ICVL,0))
  200 continue
!-----------------------------
! --- Low and Upper region
!-----------------------------
      do 210 ICVL=ICVS,ICVE
      do j=1,IQ(ICVL,2)
      a(ICVL,j)=aa(ICVL,j)*d(ICVL)*d(ip(ICVL,j))*a(ICVL,0)
      enddo
  210 continue
!
      rmaxs=0.d0
      do 220 ICVL=ICVS,ICVE
      bb(ICVL)=bb(ICVL)*d(ICVL)*a(ICVL,0)
      x(ICVL)=0.d0
      rmaxs=max(rmaxs,abs(bb(ICVL)))
  220 continue
!
      d(ICVS:ICVE)=1.d0
      a(ICVS:ICVE,0)=0.d0
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- one iteration at least
      if(cn=='y') then
      else
        if(RMAXS.le.aeps) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rmaxs0=RMAXS
!---------------------------------------
! --- all below by NCVIN
!---------------------------------------
!--< 2.2 incomplete LU decomposition >--
!---------------------------------------
      do 230 ICVL=ICVS,ICVE
      if(d(ICVL).le.0.d0) then
        write(*,'(a,E14.4,a)') 
     &  '### interface error -2- (utl_bcgstb)',d(ICVL),cn
        call FFRABORT(1,'interface error -2- utl_bcgstb')
      endif
      d(ICVL)=1.d0/(d(ICVL))
! --- upper region
      do 231 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      k=ip(ICVL,j)
      c1=0.d0
! --- low
      do 232 l=1,IQ(k,1)
      if(ip(k,l).eq.ICVL) c1=d(ICVL)*a(k,l)
  232 continue
      d(k)=d(k)-c1*a(ICVL,j)
  231 continue
  230 continue
!
!-< 3. Initial set for iteration >-
!
  888 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X)
      endif
!
! --- 
!
      rmaxs=0.d0
      do 300 ICVL=ICVS,ICVE     !ICV=1,NCVIN
      dum1=x(ICVL)
! --- upper and low region:
      do 301 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*x(ip(ICVL,j))
  301 continue
      r0(ICVL)=bb(ICVL)-dum1
      rmaxs=max(rmaxs,abs(r0(ICVL)))
  300 continue
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(cn=='y') then
        if(MXCOMP>1) then
        else
          if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
        endif
      else
        if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      s=0.d0
      do 310 ICVL=ICVS,ICVE     !i=1,NCVIN
        p0(ICVL)=r0(ICVL)
        r1(ICVL)=r0(ICVL)
        dum1=sign(1.d0,r0(ICVL))*max(1.d0,abs(r0(ICVL)))
        s=s+dum1*r0(ICVL)
        r0(ICVL)=dum1
  310 continue
!
      if(NPE.gt.1) then
          CALL hpcrsum(s)
      endif
      c1=abs(s)
!
      if(c1==0.d0) jer=1
!
      if(jer.gt.0) then
        do 320 ICVL=ICVS,ICVE     !i=1,NCVIN
          dum1=0.d0
! --- low region:
          do 321 j=1,IQ(ICVL,1)
            dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  321     continue
          p1(ICVL)=d(ICVL)*(r1(ICVL)-dum1)
  320   continue
!
        IF(NPE.gt.1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
        endif
!
        do 322 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper region:
        do 323 j=IQ(ICVL,1)+1,IQ(ICVL,2)
        dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  323   continue
        p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
        x(ICVL)=x(ICVL)+p1(ICVL)
  322   continue
!
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. Iteration >-
!
! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      do 1000 itr=it1,it2
!
        IF(NPE.gt.1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
        endif
!
      do 400 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=0.d0
! --- low
      do 401 j=1,IQ(ICVL,1)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  401 continue
      p1(ICVL)=d(ICVL)*(p0(ICVL)-dum1)
  400 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif
!
! --- 
!
      do 402 ICVL=ICVE,ICVS,-1    !i=NCVIN,1,-1
      dum1=0.d0
! --- upper
      do 403 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  403 continue
      p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
  402 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif
!
      v=0.d0
      do 410 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=p1(ICVL)
! --- low + upper
      do 411 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  411 continue
      q1(ICVL)=dum1
      v=v+r0(ICVL)*dum1
  410 continue
!
      if(NPE.gt.1) then
        CALL hpcrsum(v)
      endif      
!
      c1=abs(v)
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      alpha=s/(v) 
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
      endif
!420
      do ICVL=ICVS,ICVE    !i=1,NCVIN
        q0(ICVL)=r1(ICVL)-alpha*q1(ICVL)
      enddo
!
      do ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=0.d0
! --- low
        do 431 j=1,IQ(ICVL,1)
          dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  431   continue
        r1(ICVL)=d(ICVL)*(q0(ICVL)-dum1)
      enddo
!
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
      endif
!
      do 432 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper
        do 433 j=IQ(ICVL,1)+1,IQ(ICVL,2)
          dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  433   continue
        r1(ICVL)=r1(ICVL)-d(ICVL)*dum1
  432 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
      endif
!
      u=0.d0
      v=0.d0
      do 440 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=r1(ICVL)
        do 441 j=1,IQ(ICVL,2)
! --- low+upper
        dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  441   continue
        r2(ICVL)=dum1
        u=u+dum1*q0(ICVL)
        v=v+dum1*dum1
  440 continue
!
      if(NPE.gt.1) then
          CALL hpcrsum(u)
          CALL hpcrsum(v)
      endif
      c1=abs(v)
!
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      omega=u/(v)
!
      rmaxs=0.d0
      do 450 ICVL=ICVS,ICVE
        x(ICVL)=x(ICVL)+(alpha*p1(ICVL)+omega*r1(ICVL))
        r1(ICVL)=q0(ICVL)-omega*r2(ICVL)
        rmaxs=max(rmaxs,abs(r1(ICVL)))
  450 continue
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!
      if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
!
      t=0.d0
      do 460 ICVL=ICVS,ICVE   !i=1,NCVIN
        t=t+r0(ICVL)*r1(ICVL)
  460 continue
!
      IF(NPE.gt.1) THEN
        CALL hpcrsum(t)
      ENDIF
!
      c1=abs(s*omega)
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c1)
      ENDIF
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      beta=t*alpha/(SML+s*omega)
!
      do 470 ICVL=ICVS,ICVE
      p0(ICVL)=r1(ICVL)+beta*(p0(ICVL)-omega*q1(ICVL))
  470 continue
      s=t
!
 1000 continue
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &  ' ### WRN3: BiCGSTAB Solver is NOT converged: ',
     &  'Bicg_Iter= ', maxitr,' NOT stop and go into next time step'
         endif
         goto 999
      endif
!
      ier=1
!
  999 continue
!
      aeps=max(aeps,rmaxs)
      if(rmaxs0.lt.SML) then
        reps=max(reps,1.d0)
      else
        reps=max(reps,rmaxs/rmaxs0)
      endif
!
        do 500 ICVL=ICVS,ICVE   !i=1,NCV
        bb(ICVL)=x(ICVL)*z(ICVL)
  500   continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      if(.not.ical_vect) then
        deallocate(a)
      endif
!
      return
      end subroutine utl_bcgstb_MAT
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_bcgstb_MAT1
     &  (cn,MXALLCV,MXCV,MAXIE1,MXMAT,IQMXCV,NMAT,
     &   NALLCV,NCV,NCVIN,MAXIE,IQ,ip,aa,
     &   MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &   bb,aeps,reps,itr,ier,IMODE)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. Bi-CGstab solver with diagonal scaling
!     & incomplete LU decomposition
!
!     nn   : no. of elements
!     MAXIE   : maximul no. of off-diagonal elements
!     IQ   : no. of off-diagonal elements
!     ip   : address of off-diagonal element
!     aa   : coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!     MAT_INDEX=>NCVIN
!     MAT_CVEXT=>NCV
!     MAT_DCIDX=NCV~NALLCV
!
! --- [module arguments]
!
      use module_hpcutil
      use module_metrix,only : D  =>W1K1
      use module_metrix,only : p0 =>W1K2
      use module_metrix,only : p1 =>W1K3
      use module_metrix,only : q0 =>W1K4
      use module_metrix,only : q1 =>W1K5
      use module_metrix,only : r0 =>W1K6
      use module_metrix,only : r1 =>W1K7
      use module_metrix,only : r2 =>W1K8
      use module_metrix,only : x  =>W1K9
      use module_metrix,only : z  =>W1K10
      use module_metrix,only : a  =>W2K1
      use module_io,only       : ifli,ifll,ifle
!      use module_material,only : KMAT_S,MAT_ALL,MAT_S
      use module_cgsolver,only : nostop
      use module_model,only    : ical_vect,nthrds
      use module_dimension,only: MXCOMP
!
      implicit none
!
! --- [dummy arguments]
!
      character,intent(in)  :: cn
      integer,intent(in)    :: MXALLCV,MXCV,MAXIE1,MXMAT,NMAT,IQMXCV
      integer,intent(in)    :: NALLCV,NCV,NCVIN,MAXIE,IMODE
      integer,intent(in)    :: IQ(IQMXCV,2)
      integer,intent(in)    :: ip(MXALLCV,MAXIE)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO(   0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: aa(MXALLCV,0:MAXIE)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0

      real*8 :: c1,t,s,u,v,rmaxs,rmaxs0
      real*8 :: alpha,beta,omega,dum1
      integer:: maxitr=300,IITER=0,IELOC=0
      integer:: i,j,k,l,it1,it2,jer,m1,ICV,ierr,KMAT
      integer:: IMAT,IIMAT,ICVS,ICVE,ICVL,ICVE1
      
!
!
      t=0.d0
      s=0.d0
      u=0.d0
      v=0.d0
      alpha=0.d0
      beta=0.d0
      omega=0.d0
      ierr=0
      IITER=0
!
      if(.not.ical_vect) then
        allocate(a(MXCV,0:MAXIE),
     &          stat=ierr)
      endif
!
      if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in utl_bcgstb'
        call FFRABORT(1,'utl_bcgstb')
      endif

!
!-< 1. Preliminary set >-
!
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      rmaxs=1.d0
      DO ICVL=ICVS,ICVE
      rmaxs=min(rmaxs,abs(aa(ICVL,0)))
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmin(RMAXS)
      ENDIF
!
      if(RMAXS.lt.0.d0) then
         write(*,'(a,E14.4,a)') 
     &  '### interface error -1- (utl_bcgstb)',rmaxs,cn
        CALL FFRABORT(1,'utl_bcgstb')
      elseif(RMAXS.eq.0.d0) then
        deallocate(a)
        return
      endif
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
!
      d (:)=0.d0
      p0(:)=0.d0
      p1(:)=0.d0
      q0(:)=0.d0
      q1(:)=0.d0
      r0(:)=0.d0
      r1(:)=0.d0
      r2(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0
!-----------------------------
!-< 2. Preconditioning >-
!-----------------------------
!-----------------------------
!--< 2.1 diagonal scaling >---
!-----------------------------
!---------------------------------------------------------------
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!      ICVE=MAT_INDEX(NMAT)
!---------------------------------------------------------------
      do 200 ICVL=ICVS,ICVE
      d(ICVL)=1.d0/(dsqrt(abs(aa(ICVL,0))))
      z(ICVL)=d(ICVL)
      a(ICVL,0)=sign(1.d0,aa(ICVL,0))
  200 continue
!-----------------------------
! --- Low and Upper region
!-----------------------------
      do 210 ICVL=ICVS,ICVE
      do j=1,IQ(ICVL,2)
      a(ICVL,j)=aa(ICVL,j)*d(ICVL)*d(ip(ICVL,j))*a(ICVL,0)
      enddo
  210 continue
!
!
      do 220 ICVL=ICVS,ICVE
      bb(ICVL)=bb(ICVL)*d(ICVL)*a(ICVL,0)
      x(ICVL)=0.d0
  220 continue
!
      rmaxs=0.d0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      rmaxs=max(rmaxs,abs(bb(ICVL)))
      enddo
      enddo


      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      d(ICVS:ICVE)=1.d0
      a(ICVS:ICVE,0)=0.d0
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
      rmaxs0=RMAXS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- one iteration at least
      if(cn=='y'.or.cn=='V'.or.IMODE==6) then
!      if(cn=='y'.or.cn=='V') then
      else
        if(RMAXS.le.aeps) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------
! --- all below by NCVIN 
!---------------------------------------
!--< 2.2 incomplete LU decomposition >--
!---------------------------------------
!6666
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      
!      ICVS=1
!      ICVE=MAT_CVEXT(NMAT)
      do 230 ICVL=ICVS,ICVE
      if(d(ICVL).le.0.d0) then
        write(*,'(a,E14.4,a)') 
     &  '### interface error -2- (utl_bcgstb)',d(ICVL),cn
        call FFRABORT(1,'interface error -2- utl_bcgstb')
      endif
      d(ICVL)=1.d0/(d(ICVL))
! --- upper region
      do 231 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      k=ip(ICVL,j)
      c1=0.d0
! --- low
      do 232 l=1,IQ(k,1)
      if(ip(k,l).eq.ICVL) c1=d(ICVL)*a(k,l)
  232 continue
      d(k)=d(k)-c1*a(ICVL,j)
  231 continue
  230 continue
      enddo
!
!-< 3. Initial set for iteration >-
!
  888 continue
!9999-1
!      IF(NPE.gt.1) then
!        CALL SOLVER_SEND_RECV(1,MXCV,NCV,D)
!      endif
!
! --- 
!6666
!      rmaxs=0.d0
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!      do IIMAT=1,NMAT
!      ICVS=MAT_CVEXT(IIMAT-1)+1
!      ICVE=MAT_INDEX(IIMAT)
      do 300 ICVL=ICVS,ICVE     !ICV=1,NCVIN
      dum1=x(ICVL)
! --- upper and low region:
      do 301 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*x(ip(ICVL,j))
  301 continue
      r0(ICVL)=bb(ICVL)-dum1
!      rmaxs=max(rmaxs,abs(r0(ICVL)))
  300 continue
!      enddo
!9999-2
!      IF(NPE.gt.1) then
!        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R0)
!!        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
!      endif

      rmaxs=0.d0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE     !ICV=1,NCVIN
      rmaxs=max(rmaxs,abs(r0(ICVL)))
      enddo
      enddo      
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      if(cn=='y') then
        if(MXCOMP>1.and.IITER==0) then
        else
          if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
        endif
      elseif(cn=='V'.or.IMODE==6) then
!      elseif(cn=='V') then
        if(itr>1) then
          if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
        endif
      else
        if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!6666
      s=0.d0
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)

      do IIMAT=1,NMAT
      ICVE1=MAT_INDEX(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
!6666
      do 310 ICVL=ICVS,ICVE     !i=1,NCVIN
        p0(ICVL)=r0(ICVL)
        r1(ICVL)=r0(ICVL)
        dum1=sign(1.d0,r0(ICVL))*max(1.d0,abs(r0(ICVL)))
        if(ICVL<=ICVE1) s=s+dum1*r0(ICVL)
!        s=s+dum1*r0(ICVL)
        r0(ICVL)=dum1
  310 continue
      enddo
!
      if(NPE.gt.1) then
          CALL hpcrsum(s)
      endif
      c1=abs(s)
!
      if(c1==0.d0) jer=1
!
      
      if(jer.gt.0) then
        ICVS=1
        ICVE=MAT_CVEXT(NMAT)
        do 320 ICVL=ICVS,ICVE     !i=1,NCVIN
          dum1=0.d0
! --- low region:
          do 321 j=1,IQ(ICVL,1)
            dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  321     continue
          p1(ICVL)=d(ICVL)*(r1(ICVL)-dum1)
  320   continue
!9999-3
!        IF(NPE.gt.1) then
!          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
!        endif
!
        do 322 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper region:
        do 323 j=IQ(ICVL,1)+1,IQ(ICVL,2)
        dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  323   continue
        p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
        x(ICVL)=x(ICVL)+p1(ICVL)
  322   continue
!
        jer=0
        it1=it1+1
        itr=it1
        IITER=1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. Iteration >-
!
! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      do 1000 itr=it1,it2
!9999-4
        IF(NPE.gt.1) then
          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P0)
        endif
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!      ICVE1=MAT_INDEX(NMAT)
      do 400 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=0.d0
! --- low
      do 401 j=1,IQ(ICVL,1)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  401 continue
      p1(ICVL)=d(ICVL)*(p0(ICVL)-dum1)
  400 continue
!9999-5
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif
!
! --- 
!
      do 402 ICVL=ICVE,ICVS,-1    !i=NCVIN,1,-1
      dum1=0.d0
! --- upper
      do 403 j=IQ(ICVL,1)+1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  403 continue
      p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
  402 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
      endif

!6666 important
      v=0.d0
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)

      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      ICVE1=MAT_INDEX(IIMAT)

      do 410 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=p1(ICVL)
! --- low + upper
      do 411 j=1,IQ(ICVL,2)
      dum1=dum1+a(ICVL,j)*p1(ip(ICVL,j))
  411 enddo
      q1(ICVL)=dum1
      if(ICVL<=ICVE1) v=v+r0(ICVL)*dum1
!!!!      v=v+r0(ICVL)*dum1
  410 enddo
      enddo

!
      if(NPE.gt.1) then
        CALL hpcrsum(v)
      endif      
!
      c1=abs(v)
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      alpha=s/(v) 
!9999-6
!      IF(NPE.gt.1) then
!        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
!      endif
!420
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do ICVL=ICVS,ICVE    !i=1,NCVIN
        q0(ICVL)=r1(ICVL)-alpha*q1(ICVL)
      enddo
!9999-7 important
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,Q0)
      endif
!
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!      ICVE=MAT_INDEX(IIMAT)
!      do IIMAT=1,NMAT
!      ICVS=MAT_CVEXT(IIMAT-1)+1
!      ICVE=MAT_INDEX(IIMAT)
!!      ICVE=MAT_CVEXT(IIMAT)
      do ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=0.d0
! --- low
        do 431 j=1,IQ(ICVL,1)
          dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  431   continue
        r1(ICVL)=d(ICVL)*(q0(ICVL)-dum1)
      enddo
!      enddo
!
!9999-8 important
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
      endif
!
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!      ICVE=MAT_INDEX(IIMAT)
!      do IIMAT=1,NMAT
!      ICVS=MAT_CVEXT(IIMAT-1)+1
!      ICVE=MAT_INDEX(IIMAT)
      do 432 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper
        do 433 j=IQ(ICVL,1)+1,IQ(ICVL,2)
          dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  433   continue
        r1(ICVL)=r1(ICVL)-d(ICVL)*dum1
  432 continue
!      enddo
!9999-9 important
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
      endif
!
      u=0.d0
      v=0.d0
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)         !important
      ICVE1=MAT_INDEX(IIMAT)
      do 440 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=r1(ICVL)
        do 441 j=1,IQ(ICVL,2)
! --- low+upper
        dum1=dum1+a(ICVL,j)*r1(ip(ICVL,j))
  441   continue
        r2(ICVL)=dum1
        if(ICVL<=ICVE1) then
          u=u+dum1*q0(ICVL)
          v=v+dum1*dum1
        endif
  440 continue
      enddo
!
!
      if(NPE.gt.1) then
          CALL hpcrsum(u)
          CALL hpcrsum(v)
      endif
      c1=abs(v)
!
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      omega=abs(u)/(v)
!
      dum1=0.d0
      rmaxs=0.d0
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!      
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      ICVE1=MAT_INDEX(IIMAT)
      do 450 ICVL=ICVS,ICVE
        x(ICVL)=x(ICVL)+(alpha*p1(ICVL)+omega*r1(ICVL))
        dum1=dum1+x(ICVL)
        r1(ICVL)=q0(ICVL)-omega*r2(ICVL)
        if(ICVL<=ICVE1) rmaxs=max(rmaxs,abs(r1(ICVL)))
!        rmaxs=max(rmaxs,abs(r1(ICVL)))
  450 enddo
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!
      if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
!
!      IIMAT=1
!      ICVS=1
!      ICVE=MAT_CVEXT(NMAT)
      t=0.d0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do 460 ICVL=ICVS,ICVE   !i=1,NCVIN
      t=t+r0(ICVL)*r1(ICVL)
  460 enddo
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrsum(t)
      ENDIF
!
      c1=abs(s*omega)
!      IF(NPE.gt.1) THEN
!        CALL hpcrsum(c1)
!      ENDIF
      if(c1==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      beta=t*alpha/(s*omega)
!
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 470 ICVL=ICVS,ICVE
      p0(ICVL)=r1(ICVL)+beta*(p0(ICVL)-omega*q1(ICVL))
  470 continue
      s=t
!
 1000 continue
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &  ' ### WRN4: BiCGSTAB Solver is NOT converged: ',
     &  'Bicg_Iter= ', maxitr,' NOT stop and go into next time step'
         endif
         goto 999
      endif
!
      ier=1
!
  999 continue
!
      aeps=max(aeps,rmaxs)
      if(rmaxs0.lt.SML) then
        reps=max(reps,1.d0)
      else
        reps=max(reps,rmaxs/rmaxs0)
      endif
!
      IIMAT=1
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 500 ICVL=ICVS,ICVE   !i=1,NCV
      bb(ICVL)=x(ICVL)*z(ICVL)
  500 continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,BB)
      ENDIF
!
      if(.not.ical_vect) then
        deallocate(a)
      endif
!
      return
      end subroutine utl_bcgstb_MAT1

!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_iccg_hpc(iter,
     &  MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &  MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &  IL,ip,a,bb,aeps,reps,itr,IMODE,ier)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. CG solver with diagonal scaling
!     & incomplete Cholesky decomposition
!
!     NCV   : no. of elements
!     mm   : maximul no. of off-diagonal elements
!     IL   : no. of off-diagonal elements
!     ip   : address of off-diagonal element
!     au   : upper triangular part of coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!
      use module_hpcutil
      use module_constant
      use module_cgsolver,only : nostop
      use module_io,only     : ifll,ifle
      use module_metrix,only : d=> W1K1
      use module_metrix,only : p=> W1K2
      use module_metrix,only : q=> W1K3
      use module_metrix,only : r=> W1K4
      use module_metrix,only : x=> W1K5
      use module_metrix,only : z=> W1K6
!      use module_material,only : KMAT_S,MAT_ALL,MAT_S
      use module_metrix,only : c1,c2,c3,alpha,beta,
     &                         rmaxs0,rmaxs
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,MXALLCV,MXCV,MAXIE,MXMAT,IMODE
      integer,intent(in)    :: NMAT,NALLCV,NCV,NCVIN
      integer,intent(in)    :: IL(MXCV)
      integer,intent(in)    :: ip(MXALLCV,MAXIE)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: a(MXALLCV,0:MAXIE)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
!
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,ICV,KMAT
      integer :: i,j,k,l,it1,it2,jer,ierr,maxitr=300
      real*8  :: alf,sm,c11,c22,c33
      real*8 ,parameter :: XHALF=0.50D0
      logical :: ICHECK
!
!      allocate(c1(0:MXMAT),c2(0:MXMAT),c3(0:MXMAT),
!     &         alpha(MXMAT),beta(MXMAT))
!-< 1. preliminary set >-
!
      rmaxs(:)=1.d0
      do 100 IIMAT=1,NMAT    !i=1,NCV
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT).or.IMAT<0) goto 100
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      DO ICVL=ICVS,ICVE
      rmaxs(IIMAT)=min(rmaxs(IIMAT),a(ICVL,0))
      enddo
  100 continue
!
      IF(NPE>1) THEN
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcrmin(RMAXS(IIMAT))
!        RMAXS(0)=1.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrmin(RMAXS(IIMAT))
        enddo
      ENDIF

      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(RMAXS(IIMAT).le.0.d0) then
        write(*,*) '### interface error -1- (utl_iccg)',
     &          IIMAT,RMAXS(IIMAT)
        CALL FFRABORT(1,'utl_iccg')
      endif
      enddo
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
!
      d(:)=0.d0
      p(:)=0.d0
      q(:)=0.d0
      r(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0
!
!-< 2. Preconditioning >-
!
!--< 2.1 diagonal scaling >--
!
      do 2000 IIMAT=1,NMAT    !i=1,NCV
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 2000
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      DO ICVL=ICVS,ICVE
      d(ICVL)=1.d0/(dsqrt(a(ICVL,0)))
      z(ICVL)=d(ICVL)
      enddo
 2000 continue
!
      rmaxs(:)=0.d0
      do 2400 IIMAT=1,NMAT    !i=1,NCV
      
      if(.not.mat_cal(IIMAT)) goto 2400
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
!
      do 210 ICVL=ICVS,ICVE    !i=1,NCV
      do j=1,IL(ICVL)
      a(ICVL,j)=a(ICVL,j)*d(ICVL)*d(ip(ICVL,j))
      enddo
  210 continue
!
      do 220 ICVL=ICVS,ICVE    !i=1,NCV
      bb(ICVL)=bb(ICVL)*d(ICVL)
      rmaxs(IIMAT)=(max(rmaxs(IIMAT),abs(bb(ICVL))))
  220 continue
!
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      d(ICVS:ICVE)=1.d0
      x(ICVS:ICVE)=0.d0
 2400 continue
!
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcrmax(RMAXS(IIMAT))
!        RMAXS(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrmax(RMAXS(IIMAT))
        enddo
      ENDIF
!
      ICHECK=.true.
      do IIMAT=1,NMAT 
      ICHECK=ICHECK.and.(RMAXS(IIMAT).le.aeps)
      enddo
      IF(NPE.gt.1) call hpcland(ICHECK,3)
      rmaxs0=1.d0
      if(ICHECK) goto 999
      rmaxs0=rmaxs
!
!--< 2.2 incomplete Cholesky decomposition >--
!
      c11=1.d0
      do 3000 IIMAT=1,NMAT    
      if(.not.mat_cal(IIMAT)) goto 3000
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do 230 ICVL=ICVS,ICVE    !i=1,NCVIN
      if(d(ICVL).le.0.d0) then
        c11=d(ICVL)
      endif
      d(ICVL)=1.d0/(d(ICVL))
      do 231 j=1,IL(ICVL)
      k=ip(ICVL,j)
      d(k)=d(k)-d(ICVL)*a(ICVL,j)*a(ICVL,j)
  231 continue
  230 continue
 3000 continue
      if(NPE>1) then
        call hpcrmin(c11)
      endif
      if(c11.lt.0.d0) then
        write(*,*) '### interface error -2- (utl_iccg)'
        CALL FFRABORT(1,'utl_iccg')
      endif
!
!-< 3. Initial set for iteration >-
!
  888 continue
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X)
      endif
!
      do 4000 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 4000
      r(ICVS:ICVE)=x(ICVS:ICVE)
 4000 continue
!
      rmaxs(:)=0.d0
      do 4100 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 4100
      do 301 ICVL=ICVS,ICVE  !i=1,NCV
      do 301 j=1,IL(ICVL)
      k=ip(ICVL,j)
      r(ICVL)=r(ICVL)+a(ICVL,j)*x(k)
      r(k)=r(k)+a(ICVL,j)*x(ICVL)
  301 continue
 4100 continue
!MAT_CVEXT
      do 4200 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 4200
      do 302 ICVL=ICVS,ICVE    !i=1,NCV
      r(ICVL)=bb(ICVL)-r(ICVL)
      rmaxs(IIMAT)=max(rmaxs(IIMAT),abs(r(ICVL)))
      a(ICVL,0)=0.d0
  302 continue
 4200 continue
!
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcrmax(RMAXS(IIMAT))
!        RMAXS(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrmax(RMAXS(IIMAT))
        enddo
      ENDIF

      ICHECK=.true.
      DO IIMAT=1,NMAT
      ICHECK=ICHECK.and.
     &   (RMAXS(IIMAT).le.max(aeps,reps*rmaxs0(IIMAT)))
      enddo
      IF(NPE.gt.1) call hpcland(ICHECK,3)
      if(ICHECK) goto 999


!MAT_INDEX
      do 5000 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 5000
      do 310 ICVL=ICVS,ICVE    !i=1,NCV
      p(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
      do 311 j=1,IL(ICVL)
      k=ip(ICVL,j)
      a(k,0)=a(k,0)+a(ICVL,j)*p(ICVL)
  311 continue
  310 continue
 5000 continue
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P)
      endif
!
      c1(:)=0.d0
      do 5100 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 5100
      do 320 ICVL=ICVE,ICVS,-1    !i=NCV,1,-1
        sm=0.d0
        do 321 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*p(ip(ICVL,j))
  321   continue
        p(ICVL)=p(ICVL)-d(ICVL)*sm
        c1(IIMAT)=c1(IIMAT)+r(ICVL)*p(ICVL)
  320 continue
 5100 continue
!
! --- c1
!
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          call hpcrsum(c1(IIMAT))
!          c1(0)=0.d0
!        ENDDO
        DO IIMAT=1,NMAT
        call hpcrsum(c1(IIMAT))
        enddo
      endif
      c11=0.d0
      do IIMAT=1,NMAT
      c11=c11+c1(IIMAT)
      enddo
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c11)
      ENDIF
      if(c11.eq.0.d0) jer=1
!
      if(jer.gt.0) then
        do 7500 IIMAT=1,NMAT
        IMAT=MAT_NO(IIMAT)
        ICVS=MAT_CVEXT(IIMAT-1)+1
        ICVE=MAT_INDEX(IIMAT)
        if(.not.mat_cal(IIMAT)) goto 7500
        do 330 ICVL=ICVS,ICVE   !i=1,NCV
        x(ICVL)=x(ICVL)+p(ICVL)
  330   continue
 7500   continue
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. iteration >-
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      do 1000 itr=it1,it2
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P)
      endif
!
      do 7000 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 7000
      do 400 ICVL=ICVS,ICVE
        q(ICVL)=p(ICVL)
  400 enddo
 7000 enddo
!
      c2(:)=0.d0
      do 7700 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 7700
      do 401 ICVL=ICVS,ICVE
      do 401 j=1,IL(ICVL)
        k=ip(ICVL,j)
        q(ICVL)=q(ICVL)+a(ICVL,j)*p(k)
        q(k)=q(k)+a(ICVL,j)*p(ICVL)
  401 continue
      do 402 ICVL=ICVS,ICVE   !i=1,NCV
        c2(IIMAT)=c2(IIMAT)+p(ICVL)*q(ICVL)
  402 continue
 7700 continue
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,Q) !increase
      endif
!---------
! --- c2  
!---------
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          call hpcrsum(c2(IIMAT))
!          c2(0)=0.d0
!        ENDDO
        DO IIMAT=1,NMAT
        call hpcrsum(c2(IIMAT))
        enddo
      endif
      c22=0.d0
      do IIMAT=1,NMAT
      c22=c22+c2(IIMAT)
      enddo
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c22)
      ENDIF
      if(c22==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!------------------
! --- alpha
!------------------
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) c2(IIMAT)=1.d0
      alpha(IIMAT)=c1(IIMAT)/(c2(IIMAT))
      enddo
!
      rmaxs(:)=0.d0
      do 7100 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 7100
      do 410 ICVL=ICVS,ICVE
        x(ICVL)=x(ICVL)+alpha(IIMAT)*p(ICVL)
        r(ICVL)=r(ICVL)-alpha(IIMAT)*q(ICVL)
        rmaxs(IIMAT)=max(rmaxs(IIMAT),abs(r(ICVL)))
!        a(ICVL,0)=0.d0
  410 continue
 7100 continue
      a(:,0)=0.d0
!
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!        IIMAT=MAT_S(KMAT)
!        CALL hpcrmax(RMAXS(IIMAT))
!        RMAXS(0)=0.d0
!        enddo
        DO IIMAT=1,NMAT
        CALL hpcrmax(RMAXS(IIMAT))
        enddo
      ENDIF
!
      ICHECK=.true.
      DO IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) cycle
      ICHECK=ICHECK.and.
     &   (RMAXS(IIMAT).le.min(aeps,reps*rmaxs0(IIMAT)))
      enddo
      IF(NPE.gt.1) call hpcland(ICHECK,3)
      if(ICHECK) goto 999
!
!      IF(NPE.gt.1) THEN
!        CALL SOLVER_SEND_RECV(1,MXCV,NCV,r)
!      ENDIF
!
      do 7150 IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 7150
      do 420 ICVL=ICVS,ICVE   !i=1,NCV
        q(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
        do 421 j=1,IL(ICVL)
          k=ip(ICVL,j)
          a(k,0)=a(k,0)+a(ICVL,j)*q(ICVL)
  421   continue
  420 continue
 7150 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,Q)    !important
      endif
!
      c3(:)=0.d0
      do 7250 IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 7250
      do 430 ICVL=ICVE,ICVS,-1   !i=NCV,1,-1
        sm=0.d0
        do 431 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*q(ip(ICVL,j))
  431   continue
        q(ICVL)=q(ICVL)-d(ICVL)*sm
        c3(IIMAT)=c3(IIMAT)+r(ICVL)*q(ICVL)
  430 continue
 7250 continue
!---------
! --- c3
!---------
      IF(NPE.gt.1) THEN
!        DO KMAT=1,KMAT_S
!          IIMAT=MAT_S(KMAT)
!          call hpcrsum(c3(IIMAT))
!          c3(0)=0.d0
!        ENDDO
        DO IIMAT=1,NMAT
        call hpcrsum(c3(IIMAT))
        enddo
      endif
!
      c33=0.d0
      do IIMAT=1,NMAT
      c33=c33+c3(IIMAT)
      enddo
!
      if(NPE.gt.1) then
        CALL hpcrsum(c33)
      endif
!
      if(c33==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!------------
! --- beta
!------------
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0) c1(IIMAT)=1.d0
      beta(IIMAT)=c3(IIMAT)/(c1(IIMAT))
      c1(IIMAT)=c3(IIMAT)
      enddo
      c11=c33
!
      do 7200 IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_CVEXT(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 7200
      do 440 ICVL=ICVS,ICVE    !i=1,NCV
      p(ICVL)=q(ICVL)+beta(IIMAT)*p(ICVL)
  440 continue
 7200 continue
!
 1000 continue
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &  ' ### WRN: ICCG Solver is NOT converged: ',
     & 'iccg_Iter= ', maxitr,' NOT stop and go into next time step'
         endif
         goto 999
      endif
!
      ier=1
!
  999 continue
!
      DO IIMAT=1,NMAT 
      IMAT=MAT_NO(IIMAT)
      if(IMAT<0)rmaxs0=1.d0
      aeps=max(aeps,rmaxs(IIMAT))
      reps=max(reps,rmaxs(IIMAT)/(SML+rmaxs0(IIMAT)))
      enddo
!
      do 650 IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(mat_cal(IIMAT)) then
        do 500 ICVL=ICVS,ICVE   !i=1,NCV
        bb(ICVL)=x(ICVL)*z(ICVL)
  500   continue
      endif
 650  continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXCV,NCV,BB)
      ENDIF
      return
      end subroutine utl_iccg_hpc
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_iccg_hpc_MAT(iter,ical_sld,
     &  MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &  MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &  IL,ip,a,bb,aeps,reps,itr,IMODE,ier)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. CG solver with diagonal scaling
!     & incomplete Cholesky decomposition
!
!     NCV   : no. of elements
!     mm   : maximul no. of off-diagonal elements
!     IL   : no. of off-diagonal elements
!     ip   : address of off-diagonal element
!     au   : upper triangular part of coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!
      use module_hpcutil
      use module_cgsolver,only : nostop
      use module_io,only       : ifll,ifle
      use module_metrix,only   : d=> W1K1
      use module_metrix,only   : p=> W1K2
      use module_metrix,only   : q=> W1K3
      use module_metrix,only   : r=> W1K4
      use module_metrix,only   : x=> W1K5
      use module_metrix,only   : z=> W1K6
!      use module_material,only : KMAT_S,MAT_ALL,MAT_S
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,MXALLCV,MXCV,MAXIE,MXMAT,
     &                         ical_sld,IMODE
      integer,intent(in)    :: NMAT,NALLCV,NCV,NCVIN
      integer,intent(in)    :: IL(MXCV)
      integer,intent(in)    :: ip(MXALLCV,MAXIE)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: a(MXALLCV,0:MAXIE)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
!
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,ICV,KMAT,ICVS1,ICVE1
      integer :: i,j,k,l,it1,it2,jer,ierr,maxitr=300
      real*8  :: sm,alpha,beta,c11,c22,c33
      real*8  :: rmaxs0,rmaxs,dum1
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0
!
!-< 1. preliminary set >-
!
      rmaxs=1.d0
!
      do 100 IIMAT=1,NMAT    !i=1,NCV
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      DO ICVL=ICVS,ICVE
      rmaxs=min(rmaxs,a(ICVL,0))
      enddo
 100  enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmin(RMAXS)
      ENDIF
!
      if(RMAXS.le.0.d0) then
        write(*,*) '### interface error -1- (utl_iccg_hpc_MAT)',RMAXS
        CALL FFRABORT(1,'utl_iccg_hpc_MAT')
      endif
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
!
      d(:)=0.d0
      p(:)=0.d0
      q(:)=0.d0
      r(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0
!
!-< 2. Preconditioning >-
!
!--< 2.1 diagonal scaling >--
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      DO ICVL=ICVS,ICVE
      d(ICVL)=1.d0/(dsqrt(a(ICVL,0)))
      z(ICVL)=d(ICVL)
      enddo
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
!
      do 215 ICVL=ICVS,ICVE    !i=1,NCV 
      do 210 j=1,IL(ICVL)
      k=ip(ICVL,j)
      a(ICVL,j)=a(ICVL,j)*d(ICVL)*d(k)
 210  enddo
 215  enddo
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 220 ICVL=ICVS,ICVE    !i=1,NCV
      bb(ICVL)=bb(ICVL)*d(ICVL)
 220  enddo
!
      rmaxs=0.d0
      do IIMAT=1,NMAT    !i=1,NCV
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do 225 ICVL=ICVS,ICVE    !i=1,NCV
      rmaxs=(max(rmaxs,abs(bb(ICVL))))
 225  enddo
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
      if(IMODE==7) then
        
      else
        if(RMAXS.le.aeps) goto 999
      endif
      rmaxs0=rmaxs
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      d(ICVS:ICVE)=1.d0
      x(ICVS:ICVE)=0.d0
!
!--< 2.2 incomplete Cholesky decomposition >--
!
      do IIMAT=1,NMAT 
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
!      ICVE1=MAT_INDEX(IIMAT)
      do 230 ICVL=ICVS,ICVE 
      if(d(ICVL).le.0.d0) then
        write(*,*) '### interface error -2- (utl_iccg_hpc_MAT)'
        write(*,*) '### ICVL=, d(ICVL)=',ICVL,d(ICVL)
        CALL FFRABORT(1,'utl_iccg_hpc_MAT')
      endif
      d(ICVL)=1.d0/(d(ICVL))
      do 231 j=1,IL(ICVL)
      k=ip(ICVL,j)
      d(k)=d(k)-d(ICVL)*a(ICVL,j)*a(ICVL,j)
  231 enddo
  230 enddo
      enddo
!
!-< 3. Initial set for iteration >-
!
  888 continue
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X)
      endif
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      r(ICVS:ICVE)=x(ICVS:ICVE)
!
! --- fixed 
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 301 ICVL=ICVS,ICVE  !i=1,NCV 
      do j=1,IL(ICVL)
      k=ip(ICVL,j)
      r(ICVL)=r(ICVL)+a(ICVL,j)*x(k)
      r(k)=r(k)+a(ICVL,j)*x(ICVL)
      enddo
  301 continue
! --- fixed
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 302 ICVL=ICVS,ICVE                 !i=1,NCV
      r(ICVL)=bb(ICVL)-r(ICVL)
  302 enddo
      a(:,0)=0.d0
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      rmaxs=0.d0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do 305 ICVL=ICVS,ICVE    !i=1,NCV
      rmaxs=max(rmaxs,abs(r(ICVL)))
 305  enddo
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF

      if(IMODE==7) then
        if(itr>1) then
          if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
        endif
      else
        if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
      endif
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 310 ICVL=ICVS,ICVE    !i=1,NCV
        p(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
        do 311 j=1,IL(ICVL)
          k=ip(ICVL,j)
          a(k,0)=a(k,0)+a(ICVL,j)*p(ICVL)
  311   continue
  310 continue
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P)
      endif
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 320 ICVL=ICVE,ICVS,-1    !i=NCV,1,-1
        sm=0.d0
        do 321 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*p(ip(ICVL,j))
  321   continue
        p(ICVL)=p(ICVL)-d(ICVL)*sm
  320 continue
!
! --- c11
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      c11=0.d0
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      c11=c11+r(ICVL)*p(ICVL)
      enddo
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c11)
      ENDIF
      if(c11.eq.0.d0) jer=1
!
      if(jer.gt.0) then
        ICVS=1
        ICVE=MAT_CVEXT(NMAT)
        do 330 ICVL=ICVS,ICVE   !i=1,NCV
        x(ICVL)=x(ICVL)+p(ICVL)
  330   continue
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. iteration >-
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      do 1000 itr=it1,it2
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P)
      endif
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 400 ICVL=ICVS,ICVE
        q(ICVL)=p(ICVL)
  400 continue
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 401 ICVL=ICVS,ICVE
      do 401 j=1,IL(ICVL)
        k=ip(ICVL,j)
        q(ICVL)=q(ICVL)+a(ICVL,j)*p(k)
        q(k)=q(k)      +a(ICVL,j)*p(ICVL)
  401 continue
!
      if(NPE.gt.1) then  !OK
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,Q)
      endif
!
      c22=0.d0
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do 402 ICVL=ICVS,ICVE   !i=1,NCV
        c22=c22+p(ICVL)*q(ICVL)
  402 continue
      enddo
!
! --- c2
!
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c22)
      ENDIF
      if(c22==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
! --- alpha
!
      alpha=c11/c22
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 410 ICVL=ICVS,ICVE
        x(ICVL)=x(ICVL)+alpha*p(ICVL)
        r(ICVL)=r(ICVL)-alpha*q(ICVL)
  410 continue
      a(:,0)=0.d0
!
      rmaxs=0.d0
      do IIMAT=1,NMAT
      IMAT=MAT_NO(IIMAT)
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
        rmaxs=max(rmaxs,abs(r(ICVL)))
      ENDDO
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAXS)
      ENDIF
!
      if(RMAXS.le.min(aeps,reps*rmaxs0)) goto 999
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 420 ICVL=ICVS,ICVE   !i=1,NCV
        q(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
        do 421 j=1,IL(ICVL)
          k=ip(ICVL,j)
          a(k,0)=a(k,0)+a(ICVL,j)*q(ICVL)
  421   continue
  420 continue
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,Q)    !important
      endif
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 430 ICVL=ICVE,ICVS,-1   !i=NCV,1,-1
        sm=0.d0
        do 431 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*q(ip(ICVL,j))
  431   continue
        q(ICVL)=q(ICVL)-d(ICVL)*sm
  430 continue
!
      c33=0.d0
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      do ICVL=ICVS,ICVE
      c33=c33+r(ICVL)*q(ICVL)
      enddo
      enddo
!
! --- c3
!
      if(NPE.gt.1) then
        CALL hpcrsum(c33)
      endif
      if(c33==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
! --- beta 
!
      beta=c33/c11
      c11=c33
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 440 ICVL=ICVS,ICVE
      p(ICVL)=q(ICVL)+beta*p(ICVL)
  440 continue
!
 1000 continue
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &   ' ### WRN: ICCG Solver is NOT converged: ',
     & 'iccg_Iter= ', maxitr,' NOT stop and go into next time step'
         endif
         goto 999
      endif
!
      ier=1
!
  999 continue
!
      aeps=max(aeps,rmaxs)
      reps=max(reps,rmaxs/(rmaxs0+SML))
!
      ICVS=1
      ICVE=MAT_CVEXT(NMAT)
      do 500 ICVL=ICVS,ICVE   !i=1,NCV
      bb(ICVL)=x(ICVL)*z(ICVL)

  500 continue
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXCV,NCV,BB)
      ENDIF
      return
      end subroutine utl_iccg_hpc_MAT
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_iccg(
     &  iter,MXCV,NCV,NMAT,mm,MXMAT,MXALLCV,
     &  MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &  IL,ip,a,bb,aeps,reps,itr,IMODE,ier)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. CG solver with diagonal scaling
!     & incomplete Cholesky decomposition
!
!     NCV   : no. of elements
!     mm   : maximul no. of off-diagonal elements
!     IL   : no. of off-diagonal elements
!     ip   : address of off-diagonal element
!     au   : upper triangular part of coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!
!
! --- [module arguments]
!
      use module_metrix,only : d=> W1K1
      use module_metrix,only : p=> W1K2
      use module_metrix,only : q=> W1K3
      use module_metrix,only : r=> W1K4
      use module_metrix,only : x=> W1K5
      use module_metrix,only : z=> W1K6
      use module_cgsolver,only : nostop
      use module_io,only     : ifll,ifle
      use module_hpcutil,only: my_rank
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,NCV,mm,MXCV,MXMAT,NMAT,MXALLCV
      integer,intent(in)    :: IL(MXCV),IMODE
      integer,intent(in)    :: ip(MXALLCV,mm)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: a(MXALLCV,0:mm)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,ICVS,ICVE,ICVL,ICV,maxitr=300
      integer  i,j,k,l,it1,it2,jer,ierr,itrmx
      real*8   c1,c2,c3,sm,alpha,beta
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0
      real*8,allocatable :: rmax0(:),rmax(:)
!
      allocate(rmax0(MXMAT),rmax(MXMAT))
!
      
!-< 1. preliminary set >-
!

      itrmx=0
      rmax(:)=1.d0
      do 100 IIMAT=1,NMAT    !i=1,NCV
      if(.not.mat_cal(IIMAT)) goto 100
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      DO ICVL=ICVS,ICVE
      rmax(IIMAT)=min(rmax(IIMAT),a(ICVL,0))
      enddo
!
      if(rmax(IIMAT).le.0.d0) then
       write(*,*) '### interface error -1- (utl_iccg)',IIMAT,rmax(IIMAT)
        CALL FFRABORT(1,'utl_iccg 1')
      endif
  100 continue
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
      rmax0(:)=1.d0
!
      d(:)=0.d0
      p(:)=0.d0
      q(:)=0.d0
      r(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0
!
!-< 2. Preconditioning >-
!
!--< 2.1 diagonal scaling >--
!
      IIMAT=1
      do 2000 IIMAT=1,NMAT     !i=1,NCV
      if(.not.mat_cal(IIMAT)) goto 2000
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
 
      DO ICVL=ICVS,ICVE
      d(ICVL)=1.d0/(dsqrt(a(ICVL,0)))
      z(ICVL)=d(ICVL)
      enddo
      do 210 ICVL=ICVS,ICVE     !i=1,NCV
      do j=1,IL(ICVL)
      a(ICVL,j)=a(ICVL,j)*d(ICVL)*d(ip(ICVL,j))
      enddo
  210 continue
 
      rmax(IIMAT)=0.d0
      do 220 ICVL=ICVS,ICVE     !i=1,NCV
      bb(ICVL)=bb(ICVL)*d(ICVL)
      x(ICVL)=0.d0
      rmax(IIMAT)=max(rmax(IIMAT),abs(bb(ICVL)))
      d(ICVL)=1.d0
  220 continue
 
       if(rmax(IIMAT).le.aeps) then
         goto 2000
       endif
       rmax0(IIMAT)=rmax(IIMAT)
!
!--< 2.2 incomplete Cholesky decomposition >--
!
      do 230 ICVL=ICVS,ICVE
      if(d(ICVL).le.0.d0) then
        write(ifle,'(3X,A,I4)') 
     &  '### CPU NUMBER MY_RANK= ',MY_RANK
        write(*,*) '### ICVL=, d(ICVL)=',ICVL,d(ICVL)
        write(*,*) '### interface error -2- (utl_iccg)'
        CALL FFRABORT(1,'utl_iccg')
      endif
      d(ICVL)=1.d0/(d(ICVL))
      do 231 j=1,IL(ICVL)
      k=ip(ICVL,j)
      d(k)=d(k)-d(ICVL)*a(ICVL,j)*a(ICVL,j)
  231 continue
  230 continue
!
 2000 continue
!
!-< 3. Initial set for iteration >-
!
      do 600 IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(.not.mat_cal(IIMAT)) goto 600
      if(rmax(IIMAT).le.aeps) goto 600
!
  888 continue
!
      rmax(IIMAT)=0.d0
      do ICVL=ICVS,ICVE
      r(ICVL)=x(ICVL)
      enddo
!
      do 301 ICVL=ICVS,ICVE  !i=1,NCV
      do 301 j=1,IL(ICVL)
      k=ip(ICVL,j)
      r(ICVL)=r(ICVL)+a(ICVL,j)*x(k)
      r(k)=r(k)+a(ICVL,j)*x(ICVL)
  301 continue
!
      do 302 ICVL=ICVS,ICVE    !i=1,NCV
      r(ICVL)=bb(ICVL)-r(ICVL)
      rmax(IIMAT)=max(rmax(IIMAT),abs(r(ICVL)))
      a(ICVL,0)=0.d0
  302 continue
!
      if(rmax(IIMAT).le.max(aeps,reps*rmax0(IIMAT))) 
     &   goto 600
!
      do 310 ICVL=ICVS,ICVE    !i=1,NCV
      p(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
      do 311 j=1,IL(ICVL)
      k=ip(ICVL,j)
      a(k,0)=a(k,0)+a(ICVL,j)*p(ICVL)
  311 continue
  310 continue
!
      c1=0.d0
      do 320 ICVL=ICVE,ICVS,-1    !i=NCV,1,-1
        sm=0.d0
        do 321 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*p(ip(ICVL,j))
  321   continue
        p(ICVL)=p(ICVL)-d(ICVL)*sm
        c1=c1+r(ICVL)*p(ICVL)
  320 continue
!
      if(c1.eq.0.d0) jer=1
!
      if(jer.gt.0) then
        do 330 ICVL=ICVS,ICVE   !i=1,NCV
        x(ICVL)=x(ICVL)+p(ICVL)
  330   continue
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. iteration >-
!
      do 1000 itr=it1,it2
      if(itr>itrmx) itrmx=itr
!
      do 400 ICVL=ICVS,ICVE
        q(ICVL)=p(ICVL)
  400 continue
      do 401 ICVL=ICVS,ICVE
      do 401 j=1,IL(ICVL)
        k=ip(ICVL,j)
        q(ICVL)=q(ICVL)+a(ICVL,j)*p(k)
        q(k)=q(k)+a(ICVL,j)*p(ICVL)
  401 continue
      c2=0.d0
!      
      do 402 ICVL=ICVS,ICVE   !i=1,NCV
        c2=c2+p(ICVL)*q(ICVL)
  402 continue
!
      if(c2==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      alpha=c1/(c2)
      rmax(IIMAT)=0.d0
      do 410 ICVL=ICVS,ICVE
        x(ICVL)=x(ICVL)+alpha*p(ICVL)
        r(ICVL)=r(ICVL)-alpha*q(ICVL)
        rmax(IIMAT)=max(rmax(IIMAT),abs(r(ICVL)))
        a(ICVL,0)=0.d0
  410 continue
!
      if(rmax(IIMAT).le.min(aeps,reps*rmax0(IIMAT)) ) goto 600
!
      do 420 ICVL=ICVS,ICVE   !i=1,NCV
        q(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
        do 421 j=1,IL(ICVL)
          k=ip(ICVL,j)
          a(k,0)=a(k,0)+a(ICVL,j)*q(ICVL)
  421   continue
  420 continue

      c3=0.d0
      do 430 ICVL=ICVE,ICVS,-1   !i=NCV,1,-1
        sm=0.d0
        do 431 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*q(ip(ICVL,j))
  431   continue
        q(ICVL)=q(ICVL)-d(ICVL)*sm
        c3=c3+r(ICVL)*q(ICVL)
  430 continue
!
      if(c3==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      beta=c3/(c1)
      c1=c3
      do 440 ICVL=ICVS,ICVE    !i=1,NCV
        p(ICVL)=q(ICVL)+beta*p(ICVL)
  440 continue
!
 1000 continue
!
      if(nostop==1) then
        write(ifle,'(a,I4,a,I4a)')
     &  ' ### WRN: ICCG Solver is NOT converged: IIMAT= ',IIMAT,
     & ' iccg_Iter= ', maxitr,' NOT stop and go into next time step'
       cycle
      else
        write(ifle,'(2a,I4,2(a,I4))')
     &  ' ### ERR: ICCG Solver is NOT converged: ',
     & 'iccg_Iter= ', maxitr,' IIMAT=',IIMAT,' my_rank=',my_rank
        ier=1
      endif
!
 600  continue
!
  999 continue
!
      do IIMAT=1,NMAT
      if(mat_cal(IIMAT)) then
        aeps=max(aeps,rmax(IIMAT))
        reps=max(reps,rmax(IIMAT)/(rmax0(IIMAT)))
      endif
      enddo
!
      
      do 650 IIMAT=1,NMAT
      ICVS=MAT_CVEXT(IIMAT-1)+1
      ICVE=MAT_INDEX(IIMAT)
      if(mat_cal(IIMAT)) then
        do 500 ICVL=ICVS,ICVE   !i=1,NCV
        bb(ICVL)=x(ICVL)*z(ICVL)
  500   continue
      endif
 650  continue
!
      itr=itrmx
!
      deallocate(rmax0,rmax)
      return
      end subroutine utl_iccg
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_gausel(mko,nko,n,ip,a,b,x)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- direct solver
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mko,nko,n
      integer,intent(in)  :: ip(mko,n)
      real*8 ,intent(in)  :: a(mko,n,n),b(mko,n)
      real*8 ,intent(out) :: x(mko,n)
!
! --- [local entities]
!
      integer :: i,j,l
!
!
!/forward elimination/
      do 100 i = 1, n
      do 101 l = 1, nko
      x(l,i) = b(l,ip(l,i))
  101 continue
      do 120 j = 1, i-1
      do 121 l = 1, nko
      x(l,i) = x(l,i) - a(l,ip(l,i),j)*x(l,j)
  121 continue
  120 continue
  100 continue
!
!/backward substitution/
      do 200 i = n, 1, -1
      do 220 j = i+1, n
      do 221 l = 1, nko
      x(l,i) = x(l,i) - a(l,ip(l,i),j)*x(l,j)
  221 continue
  220 continue
      do 202 l = 1, nko
      x(l,i) = x(l,i)*a(l,ip(l,i),i)
  202 continue
  200 continue
!
      return
      end subroutine utl_gausel

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_gausel_nop(mko,nko,n,a,b,x)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mko,nko,n
      real*8 ,intent(in)  :: a(mko,n,n),b(mko,n)
      real*8 ,intent(out) :: x(mko,n)
!
! --- [local entities]
!
      integer :: i,j,l
!
!
!/forward elimination/
      do 100 i = 1, n
      do 101 l = 1, nko
      x(l,i) = b(l,i)
  101 continue
      do 120 j = 1, i-1
      do 121 l = 1, nko
      x(l,i) = x(l,i) - a(l,i,j)*x(l,j)
  121 continue
  120 continue
  100 continue
!
!/backward substitution/
      do 200 i = n, 1, -1
      do 220 j = i+1, n
      do 221 l = 1, nko
      x(l,i) = x(l,i) - a(l,i,j)*x(l,j)
  221 continue
  220 continue
      do 202 l = 1, nko
      x(l,i) = x(l,i)*a(l,i,i)
  202 continue
  200 continue
!
      return
      end subroutine utl_gausel_nop

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_getpvt(mko,nko,n,a,ip,ierr)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)    :: mko,nko,n
      real*8 ,intent(inout) :: a(mko,n,n)
      integer,intent(out)   :: ip(mko,n)
      integer,intent(out)   :: ierr(mko)
!
! --- [local entities]
!
      real*8  :: wk(mko,n),amx(mko),cf0(mko),cf1(mko)
      integer :: m(mko)
      integer :: i,j,k,l,im
!
!
      do 10 l = 1, nko
      ierr(l) = 0
   10 continue
!
! /default pivot setting/
      do 100 j = 1, n
      do 101 l = 1, nko
      ip(l,j) = j
  101 continue
  100 continue
!
! /get pivot and matrix/
      do 200 k = 1, n
!
      do 201 l = 1, nko
      m(l)   = k
      amx(l) = abs(a(l,ip(l,k),k))
 201  continue
      do 220 i = k+1, n
      do 221 l = 1, nko
      if( abs(a(l,ip(l,i),k)).gt.amx(l) ) then
        m(l)   = i
        amx(l) = abs(a(l,ip(l,i),k))
      endif
  221 continue
  220 continue
      do 210 l = 1, nko
      if( m(l).ne.k ) then
        im         = ip(l,k)
        ip(l,k)    = ip(l,m(l))
        ip(l,m(l)) = im
      endif
      if( abs(a(l,ip(l,k),k)).gt.0.d0 ) then
        cf0(l) = 1.d0/a(l,ip(l,k),k)
      else
        cf0(l) = 1.d0
        ierr(l) = 1
      endif
      a(l,ip(l,k),k) = cf0(l)
  210 continue
      do 240 i = k+1, n
      do 241 l = 1, nko
      cf1(l) = a(l,ip(l,i),k)*cf0(l)
      a(l,ip(l,i),k) = cf1(l)
  241 continue
      do 260 j = k+1, n
      do 261 l = 1, nko
      wk(l,j) = a(l,ip(l,i),j) - a(l,ip(l,k),j)*cf1(l)
  261 continue
  260 continue
      do 280 j = k+1, n
      do 281 l = 1, nko
      a(l,ip(l,i),j) = wk(l,j)
  281 continue
  280 continue
  240 continue
!
  200 continue
!
      return
      end subroutine utl_getpvt

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_getpvt_nop(mko,nko,n,a,ierr)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [dummy arguments]
!
      integer,intent(in)    :: mko,nko,n
      real*8 ,intent(inout) :: a(mko,n,n)
      integer,intent(out)   :: ierr(mko)
!
! --- [local entities]
!
      real*8,parameter :: dd=1.d-70
      integer :: i,j,k,l
!
!
      do 10 l = 1, nko
      ierr(l) = 0
   10 continue
!
!/LU decomposition/
      do 200 k = 1, n
!
      do 210 l = 1, nko
      if( abs(a(l,k,k)).gt.0.d0 ) then
        a(l,k,k) = 1.0d0/a(l,k,k)
      else
        a(l,k,k) = 1.0d0
        ierr(l) = 1
      endif
  210 continue
!
      do 240 i = k+1, n
      do 241 l = 1, nko
      a(l,i,k) = a(l,i,k)*a(l,k,k)
  241 continue
      do 260 j = k+1, n
      do 261 l = 1, nko
      a(l,i,j) = a(l,i,j) - a(l,k,j)*a(l,i,k)
  261 continue
  260 continue
  240 continue
!
  200 continue
!
      return
      end subroutine utl_getpvt_nop
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_random(r)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! 1. random number generator ( 0<r<1 )
!
      implicit none
      real*8,intent(out) :: r
      integer,parameter :: k=1899, m=2**23
      integer,parameter :: mr=2**12, md=m/mr
      logical,save :: first=.true.
      integer :: ii=1, l=0, u=0, x
      INTEGER :: val(8)
      if( first ) then
        CALL DATE_AND_TIME(VALUES=val)
        ii=int(val(8))
        l=mod(ii,mr)
        u=(ii-l)/mr
        first=.false.
      endif
      x=k*l
      l=mod(x,mr)
      u=mod(k*u+(x-l)/mr,md)
      r=dble(u*mr+l)/dble(m)
      end subroutine utl_random
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine utl_iccg_hpc_vect(iter,
     &  MXALLCV,MXCV,MAXIE,MXMAT,NMAT,NALLCV,NCV,NCVIN,
     &  MAT_INDEX,MAT_CVEXT,MAT_NO,mat_cal,
     &  IL,ip,a,bb,aeps,reps,itr,ier)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil
      use module_cgsolver,only : nostop
      use module_model   ,only : nthrds
      use module_io,only       : ifll,ifle
      use module_metrix,only   : d=> W1K1
      use module_metrix,only   : p=> W1K2
      use module_metrix,only   : q=> W1K3
      use module_metrix,only   : r=> W1K4
      use module_metrix,only   : x=> W1K5
      use module_metrix,only   : z=> W1K6
!      use module_material,only : KMAT_S,MAT_ALL,MAT_S
      use module_vector,only   : ICVS_V,ICVE_V,
     &                           ICFS_V,ICFE_V,
     &                           ICVSIN_V,ICVEIN_V,
     &                           IDCS_V,IDCE_V
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: iter,MXALLCV,MXCV,MAXIE,MXMAT
      integer,intent(in)    :: NMAT,NALLCV,NCV,NCVIN
      integer,intent(in)    :: IL(MXCV)
      integer,intent(in)    :: ip(MXALLCV,MAXIE)
      INTEGER,INTENT(IN)    :: MAT_INDEX(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_CVEXT(0:MXMAT)
      INTEGER,INTENT(IN)    :: MAT_NO (  0:MXMAT)
      logical,INTENT(IN)    :: mat_cal(  0:MXMAT)      
      real*8 ,intent(inout) :: a(MXALLCV,0:MAXIE)
      real*8 ,intent(inout) :: bb(MXCV)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
!
! --- [local entities]
!
!
      integer :: ICVL
      integer :: my_s,my_e,i
      integer :: j,k,l,it1,it2,jer,ierr,maxitr=300
      real*8  :: sm,alpha,beta,c11,c22,c33
      real*8  :: rmax0,rmax,dum1
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-50,ZERO=0.D0
!-------------------
! --- openMP
!-------------------
!
!-< 1. preliminary set >-
!
      rmax=1.d0
!
      DO ICVL=ICVSIN_V,ICVEIN_V
        rmax=min(rmax,a(ICVL,0))
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmin(RMAX)
      ENDIF
!
      if(RMAX.le.0.d0) then
        write(*,*) '### interface error -1- (utl_iccg_hpc_vect)',RMAX
        CALL FFRABORT(1,'utl_iccg_hpc_vect')
      endif
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
      maxitr=it2
!
      DO ICVL=ICVSIN_V,ICVEIN_V
      d(ICVL)=0.d0
      p(ICVL)=0.d0
      q(ICVL)=0.d0
      r(ICVL)=0.d0
      x(ICVL)=0.d0
      z(ICVL)=0.d0
      enddo
!
!-< 2. Preconditioning >-
!
!--< 2.1 diagonal scaling >--
!
      DO ICVL=ICVSIN_V,ICVEIN_V
      d(ICVL)=1.d0/(dsqrt(a(ICVL,0)))
      z(ICVL)=d(ICVL)
      enddo
!
      do 215 ICVL=ICVSIN_V,ICVEIN_V    !i=1,NCV
      do 210 j=1,IL(ICVL)
      k=ip(ICVL,j)
      a(ICVL,j)=a(ICVL,j)*d(ICVL)*d(k)
 210  enddo
 215  enddo
!
      do 220 ICVL=ICVSIN_V,ICVEIN_V    !i=1,NCV
      bb(ICVL)=bb(ICVL)*d(ICVL)
 220  enddo
!
      rmax=0.d0
      do 225 ICVL=ICVSIN_V,ICVEIN_V    !i=1,NCV
      rmax=(max(rmax,abs(bb(ICVL))))
 225  enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAX)
      ENDIF
      if(RMAX.le.aeps) then
        goto 999
      endif
      rmax0=rmax
!
      do ICVL=ICVSIN_V,ICVEIN_V
      d(ICVL)=1.d0
      x(ICVL)=0.d0
      enddo
!
!--< 2.2 incomplete Cholesky decomposition >--
!
      !NOomp-1
      do 230 ICVL=ICVSIN_V,ICVEIN_V    !i=1,NCVIN
      if(d(ICVL).le.0.d0) then
        write(*,*) '### interface error -2- (utl_iccg_hpc_vectMAT)'
        write(*,*) '### ICVL=, d(ICVL)=',ICVL,d(ICVL)
        CALL FFRABORT(1,'utl_iccg_hpc_vect')
      endif
      d(ICVL)=1.d0/(d(ICVL))
      do 231 j=1,IL(ICVL)
      k=ip(ICVL,j)
      d(k)=d(k)-d(ICVL)*a(ICVL,j)*a(ICVL,j)
  231 enddo
  230 enddo
!
!-< 3. Initial set for iteration >-
!
  888 continue
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X)
      endif
!
      DO ICVL=ICVSIN_V,ICVEIN_V
      r(ICVL)=x(ICVL)
      enddo
!
! --- fixed
!
!----- openMP -2- ?------------------------------
      do 301 ICVL=ICVSIN_V,ICVEIN_V  !i=1,NCV
      do j=1,IL(ICVL)
      k=ip(ICVL,j)
      r(ICVL)=r(ICVL)+a(ICVL,j)*x(k)
      r(k)=r(k)+a(ICVL,j)*x(ICVL)
      enddo
  301 enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! --- fixed
!
      do 302 ICVL=ICVSIN_V,ICVEIN_V                 !i=1,NCV
      r(ICVL)=bb(ICVL)-r(ICVL)
  302 enddo
!
      rmax=0.d0
      do 305 ICVL=ICVSIN_V,ICVEIN_V  !i=1,NCV
      rmax=max(rmax,abs(r(ICVL)))
 305  enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAX)
      ENDIF
      if( RMAX.le.max(aeps,reps*rmax0)) goto 999
!----- openMP -3- ?   ------------------------------
!      !NOomp-2
      a(:,0)=0.d0
      do 310 ICVL=ICVSIN_V,ICVEIN_V
        p(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
        do 311 j=1,IL(ICVL)
          k=ip(ICVL,j)
          a(k,0)=a(k,0)+a(ICVL,j)*p(ICVL)
  311   enddo
  310 enddo
!----- openMP -3- ------------------------------
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P)
      endif
      do 320 ICVL=ICVEIN_V,ICVSIN_V,-1    !i=NCV,1,-1
        sm=0.d0
        do 321 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*p(ip(ICVL,j))
  321   enddo
        p(ICVL)=p(ICVL)-d(ICVL)*sm
 320  enddo
!
! --- c11
!
      c11=0.d0
      do ICVL=ICVSIN_V,ICVEIN_V 
      c11=c11+r(ICVL)*p(ICVL)
      enddo
!
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c11)
      ENDIF
      if(c11.eq.0.d0) jer=1
!
      if(jer.gt.0) then
        do 330 ICVL=ICVSIN_V,ICVEIN_V    !i=1,NCV
        x(ICVL)=x(ICVL)+p(ICVL)
  330   enddo
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. iteration >-
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      do 1000 itr=it1,it2
!
      if(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P)
      endif
!
      do 400 ICVL=ICVSIN_V,ICVEIN_V
        q(ICVL)=p(ICVL)
  400 enddo
!
      do 401 ICVL=ICVSIN_V,ICVEIN_V
      do j=1,IL(ICVL)
        k=ip(ICVL,j)
        q(ICVL)=q(ICVL)+a(ICVL,j)*p(k)
        q(k)=q(k)      +a(ICVL,j)*p(ICVL)
      enddo
  401 enddo
!
      if(NPE.gt.1) then  !OK
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,Q)
      endif
!
      c22=0.d0
      do 402 ICVL=ICVSIN_V,ICVEIN_V   !i=1,NCV
        c22=c22+p(ICVL)*q(ICVL)
  402 enddo
!
! --- c2
!
      IF(NPE.gt.1) THEN
        CALL hpcrsum(c22)
      ENDIF
      if(c22==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
! --- alpha
!
      alpha=c11/c22
!
      do 410 ICVL=ICVSIN_V,ICVEIN_V
        x(ICVL)=x(ICVL)+alpha*p(ICVL)
        r(ICVL)=r(ICVL)-alpha*q(ICVL)
        a(ICVL,0)=0.d0
  410 enddo
!
      rmax=0.d0
!
      do ICVL=ICVSIN_V,ICVEIN_V
        rmax=max(rmax,abs(r(ICVL)))
      ENDDO
!
      IF(NPE.gt.1) THEN
        CALL hpcrmax(RMAX)
      ENDIF
!
      if(RMAX.le.min(aeps,reps*rmax0)) goto 999
!--------------------
! --- forward
!--------------------
      do 420 ICVL=ICVSIN_V,ICVEIN_V   !i=1,NCV
        q(ICVL)=d(ICVL)*(r(ICVL)-a(ICVL,0))
        do 421 j=1,IL(ICVL)
          k=ip(ICVL,j)
          a(k,0)=a(k,0)+a(ICVL,j)*q(ICVL)
  421   ENDDO
  420 ENDDO
!
      IF(NPE.gt.1) then
        CALL SOLVER_SEND_RECV(1,MXCV,NCV,Q)    !important
      endif
!--------------------
! --- backward 
!--------------------
      do 430 ICVL=ICVEIN_V,ICVSIN_V,-1   !i=NCV,1,-1
        sm=0.d0
        do 431 j=1,IL(ICVL)
          sm=sm+a(ICVL,j)*q(ip(ICVL,j))
  431   ENDDO
        q(ICVL)=q(ICVL)-d(ICVL)*sm
  430 ENDDO
!
      c33=0.d0
      do ICVL=ICVSIN_V,ICVEIN_V
      c33=c33+r(ICVL)*q(ICVL)
      enddo
!
! --- c3
!
      if(NPE.gt.1) then
        CALL hpcrsum(c33)
      endif
      if(c33==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
! --- beta 
!
      beta=c33/c11
      c11=c33
!
      do 440 ICVL=ICVSIN_V,ICVEIN_V
      p(ICVL)=q(ICVL)+beta*p(ICVL)
  440 enddo
!
 1000 continue
!
      if(nostop==1) then
         if(my_rank==root) then
           write(ifle,'(2a,I4,a)')
     &   ' ### WRN: ICCG Solver is NOT converged: ',
     & 'iccg_Iter= ', maxitr,' NOT stop and go into next time step'
         endif
         goto 999
      endif
!
      ier=1
!
  999 continue
!
      aeps=max(aeps,rmax)
      reps=max(reps,rmax/(rmax0+SML))
!
      do 500 ICVL=ICVSIN_V,ICVEIN_V   !i=1,NCV
      bb(ICVL)=x(ICVL)*z(ICVL)
  500 enddo
!
      IF(NPE.GT.1) THEN
        CALL SOLVER_SEND_RECV (1,MXCV,NCV,BB)
      ENDIF

      return
      end subroutine utl_iccg_hpc_vect
!
