!
!      subroutine dc_symprs()
!      subroutine dc_symprv()
!      subroutine dc_symvel
!      subroutine dc_symvect
!      subroutine dc_vgrad()
!      subroutine dc_symprs_vof
!      subroutine dc_symMHD
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_symprs(nko,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &                     mat_cal,aa) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments] 
!
      use module_dimension
      use module_boundary,only : kdbcnd,kdprdc,LBC_INDEX,nbcnd,kdtchi,
     &                           MAT_BCIDX,kdintr,kdsld,kdcvd,idis,
     &                           kdbuff,kdshutr,lbc_pair
     &                           ,kdovst
      use module_model,only    : ical_vect
      use module_metrix,only   : tmpfac=>d2vect,tmpsld
      use module_metrix,only   : SHUTFL
!
! 1. Set scalar at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: wifsld ( MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(inout) :: aa     ( MXALLCV,nko)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,I
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,ICVBO,ICFO
      real*8  :: wi1,wi2,dum1
      integer :: ISLD,ISLD2,IBFS1,IBFE1,IBFS2,IBFE2
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc) then
! --- Periodic & interface BC:
        if(idis(nb)==0) then
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do 201 k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          aa(IDCP,k)=aa(ICV,k)
  201     continue
 1250     continue
        elseif(idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
!
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          do k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          enddo
          enddo
        endif
      elseif(kd==kdbuff) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do k=1,nko
!          dum1=0.5d0*(aa(ICVP,k)+aa(ICV,k))
!          aa(IDC,k)=dum1
!          aa(IDCP,k)=dum1
!          aa(ICV,k)=dum1
!          aa(ICVP,k)=dum1
          aa(IDC,k)=aa(ICVP,k)
          aa(IDCP,k)=aa(ICV,k)
          enddo
        enddo
      elseif(kd==kdshutr) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          if(SHUTFL(IBFL)==0) then
            do k=1,nko
            aa(IDC,k)=aa(ICVP,k)
            aa(IDCP,k)=aa(ICV,k)
            enddo
          else
            do k=1,nko
            aa(IDC,k)=aa(ICV,k)   !aa(ICVP,k)
            aa(IDCP,k)=aa(ICVP,k)   !aa(ICV,k)
            enddo
          endif
        enddo
      elseif(kd==kdovst) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICVP=LCYCSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICV=LVEDGE(1,ICFL)
          do k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          enddo
        enddo
      elseif(kd==kdsld) then
       if(ical_vect) then 
        do k=1,nko
        do ISLD=1,2
        ISLD2=3-ISLD
        IBFS1=(LBC_INDEX(nb-1)+1)*(2-ISLD)
        IBFE1=(LBC_pair(nb))*(2-ISLD)
        IBFS2=(LBC_pair(nb)+1)*(2-ISLD2)
        IBFE2=(LBC_INDEX(nb))*(2-ISLD2)
        DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
!        ICFL=LBC_SSF(IBFL)
!        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICFO=LCYCOLD(IBFL)
        ICVB=LVEDGE(1,ICFP)
        ICVBO=LVEDGE(1,ICFO)
        
        tmpsld(IBFL,1)=  !aa(IDC,k)
     &           wifsld(IBFL)*aa(ICVB,k)
     &        +(1.d0-wifsld(IBFL))*aa(ICVBO,k)
        enddo

        DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
        aa(LVEDGE(2,LBC_SSF(IBFL)),k)=tmpsld(IBFL,1)
        enddo
!
        enddo
        enddo
       else
        do 1200 ISLD=1,2
        ISLD2=3-ISLD
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
!
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICVA=LVEDGE(1,ICFL)

        ICFP=LCYCSF(IBFL)
        ICFO=LCYCOLD(IBFL)

        ICVB=LVEDGE(1,ICFP)
        ICVBO=LVEDGE(1,ICFO)
        wi1=wifsld(IBFL)
        wi2=1.d0-wi1
        do k=1,nko
        aa(IDC,k)=wi1*aa(ICVB,k)+wi2*aa(ICVBO,k)  !?zhang???? 
!        aa(IDC,k)=aa(ICVA,k)
!        aa(IDC,k)=aa(ICVB,k)
        enddo
        enddo
 1200   enddo
       endif
      else
! --- Other BC:
        if(ical_vect) then 
          do k=1,nko  !huilai
!NEC$ vector
          do IBFL=IBFS,IBFE
          tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),k)
          ENDDO
!NEC$ vector
          do IBFL=IBFS,IBFE
          aa(LVEDGE(2,LBC_SSF(IBFL)),k)=tmpfac(IBFL,1)
          ENDDO
          enddo
        else
          do 1150 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 101 k=1,nko
          aa(IDC,k)=aa(ICV,k)
  101     continue
 1150     continue
        endif
!
        if(kd==kdintr) then
          do IBFL=IBFS,IBFE
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFP)
          IDC=LVEDGE(2,ICFP)
          do k=1,nko
          aa(IDC,k)=aa(ICV,k)
          enddo
          enddo
        endif
      endif
 1000 continue
!
      return
!
      end subroutine dc_symprs
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_symprv
     &(nko,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & LCYCOLD,wifsld,OPPANG,
     & SFAREA,SFCENT,aa,imode,ivel)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     imode==0 : velocity
!     imode==1 : grdc
!     imode==2 : rvx
!
! --- [module arguments]
!
      USE module_dimension
      use module_hpcutil,only  : my_rank
      use module_boundary,only : set_rotmtrx,nbcnd,LBC_INDEX,
     &                           kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdintr,kdsld,MAT_BCIDX,rotsld,
     &                           idis,LBC_pair,kdbuff,kdshutr
     &                           ,kdovst
      use module_material,only : rotati,ishaft,end,begin,nsplit,
     &                           rot_ang
      use module_nowtime ,only : iter,time
      use module_material,only : ical_sld
      use module_model,only    : ical_vect
      use module_metrix,only   : tmpfac=>d2vect,tmpsld
      use module_metrix,only   : SHUTFL
!
! 1. Set vector at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko,imode,ivel
      integer,intent(in)    :: MAT_NO    (  0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(inout) :: aa        (  MXALLCV,3,nko)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
!
! --- [local entities]
!
      real*8  :: prd
      real*8  :: rot(3,3)
      real*8  :: shft(3),ufix(3,2)
      real*8  :: unit(3,2),th(2),bb(3,3,2)
      real*8  :: rbb(3,3,2)
      real*8  :: urot(3,2)
      real*8  :: vr(3,2),r(3,2)
      real*8  :: v0(3,2),dr(2),vell(3),wi1,wi2,dum1,dum2,dum3
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2)
      integer :: i,j,ICVBO,ICFO,ibfs1,ibfe1,ibfs2,ibfe2
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k,ISLD,ISLD2
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc)then  !zhang???
         rot=0.d0
         call set_rotmtrx(nbcnd,kd,nb,rot)
!--< 1.1 periodic boundary >--
        do 1150 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 101 k=1,nko
!
        aa(IDC,1,k)=rot(1,1)*aa(ICVP,1,k)
     &             +rot(1,2)*aa(ICVP,2,k)
     &             +rot(1,3)*aa(ICVP,3,k)
        aa(IDC,2,k)=rot(2,1)*aa(ICVP,1,k)
     &             +rot(2,2)*aa(ICVP,2,k)
     &             +rot(2,3)*aa(ICVP,3,k)
        aa(IDC,3,k)=rot(3,1)*aa(ICVP,1,k)
     &             +rot(3,2)*aa(ICVP,2,k)
     &             +rot(3,3)*aa(ICVP,3,k)
!
        if(idis(nb)==0) then   !zhang???
          aa(IDCP,1,k)=rot(1,1)*aa(ICV,1,k)
     &                +rot(2,1)*aa(ICV,2,k)
     &                +rot(3,1)*aa(ICV,3,k)
          aa(IDCP,2,k)=rot(1,2)*aa(ICV,1,k)
     &                +rot(2,2)*aa(ICV,2,k)
     &                +rot(3,2)*aa(ICV,3,k)
          aa(IDCP,3,k)=rot(1,3)*aa(ICV,1,k)
     &                +rot(2,3)*aa(ICV,2,k)
     &                +rot(3,3)*aa(ICV,3,k)
        endif
 101    continue
 1150   continue
      elseif(kd==kdbuff) then 
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
!        
        do k=1,nko 
        aa(IDC,1,k)=aa(ICVP,1,k)
        aa(IDC,2,k)=aa(ICVP,2,k)
        aa(IDC,3,k)=aa(ICVP,3,k)
        aa(IDCP,1,k)=aa(ICV,1,k)
        aa(IDCP,2,k)=aa(ICV,2,k)
        aa(IDCP,3,k)=aa(ICV,3,k)
        enddo
        enddo
      elseif(kd==kdshutr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        if(SHUTFL(IBFL)==0) then
          do k=1,nko
          aa(IDC,1,k)=aa(ICVP,1,k)
          aa(IDC,2,k)=aa(ICVP,2,k)
          aa(IDC,3,k)=aa(ICVP,3,k)
          aa(IDCP,1,k)=aa(ICV,1,k)
          aa(IDCP,2,k)=aa(ICV,2,k)
          aa(IDCP,3,k)=aa(ICV,3,k)
          enddo
        else
          do k=1,nko
          aa(IDC,1,k)=aa(ICV,1,k)   !aa(ICVP,1,k)
          aa(IDC,2,k)=aa(ICV,2,k)   !aa(ICVP,2,k)
          aa(IDC,3,k)=aa(ICV,3,k)   !aa(ICVP,3,k)
          aa(IDCP,1,k)=aa(ICVP,1,k)   !aa(ICV,1,k)
          aa(IDCP,2,k)=aa(ICVP,2,k)   !aa(ICV,2,k)
          aa(IDCP,3,k)=aa(ICVP,3,k)   !aa(ICV,3,k)
          enddo
        endif
        enddo
      elseif(kd==kdovst) then 
        IIMATS(1)=MAT_BCIDX(nb,1)
        IIMATS(2)=MAT_BCIDX(nb,2)
        IMATS(1)=MAT_NO(IIMATS(1))
        IMATS(2)=MAT_NO(IIMATS(2))
! 
        th(1)=rot_ang(IMATS(1))
        th(2)=rot_ang(IMATS(2))
! 
        unit(1,2)=end(1,IMATS(2))-begin(1,IMATS(2))
        unit(2,2)=end(2,IMATS(2))-begin(2,IMATS(2))
        unit(3,2)=end(3,IMATS(2))-begin(3,IMATS(2))
        CALL rotth(unit(:,2),th(2),bb(:,:,2))
! 
        unit(1,1)=end(1,IMATS(1))-begin(1,IMATS(1))
        unit(2,1)=end(2,IMATS(1))-begin(2,IMATS(1))
        unit(3,1)=end(3,IMATS(1))-begin(3,IMATS(1))
        CALL rotth(unit(:,1),th(1),bb(:,:,1))
! 
        do i=1,3
        do j=1,3
        rbb(i,j,1)=bb(j,i,1)
        rbb(i,j,2)=bb(j,i,2)
        enddo
        enddo
! 
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LCYCSF(IBFL)
        r(1,1)=SFCENT(1,ICFL)-begin(1,IMATS(1))
        r(2,1)=SFCENT(2,ICFL)-begin(2,IMATS(1))
        r(3,1)=SFCENT(3,ICFL)-begin(3,IMATS(1))
        call AXB_UNIT_C(unit(:,1),r(:,1),vr(:,1))
        dr(1)=r(1,1)*unit(1,1)
     &          +r(2,1)*unit(2,1)
     &          +r(3,1)*unit(3,1)
        r(1,1)=r(1,1)-dr(1)*unit(1,1)
        r(2,1)=r(2,1)-dr(1)*unit(2,1)
        r(3,1)=r(3,1)-dr(1)*unit(3,1)
        dr(1)=dsqrt(r(1,1)*r(1,1)
     &             +r(2,1)*r(2,1)
     &             +r(3,1)*r(3,1))
        v0(1,1)=dr(1)*rotati(IMATS(1))*vr(1,1)
        v0(2,1)=dr(1)*rotati(IMATS(1))*vr(2,1)
        v0(3,1)=dr(1)*rotati(IMATS(1))*vr(3,1)
!
        dum1=SFCENT(1,ICFL)-begin(1,IMATS(2)) !7777
        dum2=SFCENT(2,ICFL)-begin(2,IMATS(2)) 
        dum3=SFCENT(3,ICFL)-begin(3,IMATS(2)) 
!
        r(1,2)=rbb(1,1,1)*dum1
     &        +rbb(1,2,1)*dum2
     &        +rbb(1,3,1)*dum3
        r(2,2)=rbb(2,1,1)*dum1
     &        +rbb(2,2,1)*dum2
     &        +rbb(2,3,1)*dum3
        r(3,2)=rbb(3,1,1)*dum1
     &        +rbb(3,2,1)*dum2
     &        +rbb(3,3,1)*dum3
!
        call AXB_UNIT_C(unit(:,2),r(:,2),vr(:,2)) 
!
        dr(2)=r(1,2)*unit(1,2)
     &       +r(2,2)*unit(2,2)
     &       +r(3,2)*unit(3,2)
!
        r(1,2)=r(1,2)-dr(2)*unit(1,2)
        r(2,2)=r(2,2)-dr(2)*unit(2,2)
        r(3,2)=r(3,2)-dr(2)*unit(3,2)
        dr(2)=dsqrt(r(1,2)*r(1,2)
     &             +r(2,2)*r(2,2)
     &             +r(3,2)*r(3,2))
        v0(1,2)=dr(2)*rotati(IMATS(2))*vr(1,2)
        v0(2,2)=dr(2)*rotati(IMATS(2))*vr(2,2)
        v0(3,2)=dr(2)*rotati(IMATS(2))*vr(3,2)        
!
        IF(imode==1) THEN
          v0(:,:)=0.D0
        ENDIF
!
        do k=1,nko
        vell(1)=aa(ICVP,1,k)
        vell(2)=aa(ICVP,2,k)
        vell(3)=aa(ICVP,3,k)
! --- rot2=> fix2
        ufix(1,2)=
     &            rbb(1,1,2)*(vell(1))!+v0(1,2))
     &           +rbb(1,2,2)*(vell(2))!+v0(2,2))
     &           +rbb(1,3,2)*(vell(3))!+v0(3,2))
     &           +v0(1,2)
        ufix(2,2)=
     &            rbb(2,1,2)*(vell(1))!+v0(1,2))
     &           +rbb(2,2,2)*(vell(2))!+v0(2,2))
     &           +rbb(2,3,2)*(vell(3))!+v0(3,2))
     &           +v0(2,2)
        ufix(3,2)=
     &            rbb(3,1,2)*(vell(1))!+v0(1,2))
     &           +rbb(3,2,2)*(vell(2))!+v0(2,2))
     &           +rbb(3,3,2)*(vell(3))!+v0(3,2))
     &           +v0(3,2)
! --- fix2 => rot1
        urot(1,1)= rbb(1,1,1)*(ufix(1,2))
     &            +rbb(2,1,1)*(ufix(2,2))
     &            +rbb(3,1,1)*(ufix(3,2))
     &            -v0(1,1)
        urot(2,1)= rbb(1,2,1)*(ufix(1,2))
     &            +rbb(2,2,1)*(ufix(2,2))
     &            +rbb(3,2,1)*(ufix(3,2))
     &            -v0(2,1)
        urot(3,1)= rbb(1,3,1)*(ufix(1,2))
     &            +rbb(2,3,1)*(ufix(2,2))
     &            +rbb(3,3,1)*(ufix(3,2))
     &            -v0(3,1)
        
        aa(IDC,1,k)=urot(1,1)
        aa(IDC,2,k)=urot(2,1)
        aa(IDC,3,k)=urot(3,1)
        enddo !enddo k=1,...
        enddo
!
!        if(ivel==1) then 
!          do IBFL=IBFS,IBFE
!          ICFL=LBC_SSF(IBFL)
!          ICV=LVEDGE(1,ICFL)
!          IDC=LVEDGE(2,ICFL)
!          ICVP=LCYCSF(IBFL)
!          do k=1,nko
!          dum1=(SFAREA(1,ICFL)*aa(ICVP,1,k)
!     &         +SFAREA(2,ICFL)*aa(ICVP,2,k)
!     &         +SFAREA(3,ICFL)*aa(ICVP,3,k))
!          if(dum1>0.d0) then
!            aa(IDC,1,k)=aa(ICV,1,k)
!            aa(IDC,2,k)=aa(ICV,2,k)
!            aa(IDC,3,k)=aa(ICV,3,k)
!          else
!          endif
!          enddo
!          enddo
!        endif
      elseif(kd==kdsld.and.idis(nb)==1) then
!        call set_rotmtrx(nbcnd,kd,nb,rot)
        if(ical_vect) then 
          do k=1,nko
          do 1460 ISLD=1,2
          ISLD2=3-ISLD
          IBFS1=(LBC_INDEX(nb-1)+1)*(2-ISLD)
          IBFE1=(LBC_pair(nb))*(2-ISLD)
          IBFS2=(LBC_pair(nb)+1)*(2-ISLD2)
          IBFE2=(LBC_INDEX(nb))*(2-ISLD2)
          IMATS(ISLD)=MAT_NO(MAT_BCIDX(nb,ISLD))
          th(ISLD)=rot_ang(MAT_NO(MAT_BCIDX(nb,ISLD)))
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
          rbb(1,1,ISLD)=bb(1,1,ISLD)
          rbb(1,2,ISLD)=bb(2,1,ISLD)
          rbb(1,3,ISLD)=bb(3,1,ISLD)
          rbb(2,1,ISLD)=bb(1,2,ISLD)
          rbb(2,2,ISLD)=bb(2,2,ISLD)
          rbb(2,3,ISLD)=bb(3,2,ISLD)
          rbb(3,1,ISLD)=bb(1,3,ISLD)
          rbb(3,2,ISLD)=bb(2,3,ISLD)
          rbb(3,3,ISLD)=bb(3,3,ISLD)

          IMATS(ISLD2)=MAT_NO(MAT_BCIDX(nb,ISLD2))
          th(ISLD2)=rot_ang(IMATS(ISLD2))
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
          rbb(1,1,ISLD2)=bb(1,1,ISLD2)
          rbb(1,2,ISLD2)=bb(2,1,ISLD2)
          rbb(1,3,ISLD2)=bb(3,1,ISLD2)
          rbb(2,1,ISLD2)=bb(1,2,ISLD2)
          rbb(2,2,ISLD2)=bb(2,2,ISLD2)
          rbb(2,3,ISLD2)=bb(3,2,ISLD2)
          rbb(3,1,ISLD2)=bb(1,3,ISLD2)
          rbb(3,2,ISLD2)=bb(2,3,ISLD2)
          rbb(3,3,ISLD2)=bb(3,3,ISLD2)
! --- imode


         if(imode==0) then

          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)= !r(1,ISLD2)=
     &      rbb(1,1,ISLD)*(SFCENT(1,ICFL)
     &                     -begin(1,IMATS(ISLD2)))
     &     +rbb(1,2,ISLD)*(SFCENT(2,ICFL)
     &                     -begin(2,IMATS(ISLD2)))
     &     +rbb(1,3,ISLD)*(SFCENT(3,ICFL)
     &                     -begin(3,IMATS(ISLD2))) 
          tmpsld(IBFL,2)= !r(2,ISLD2)=
     &      rbb(2,1,ISLD)*(SFCENT(1,ICFL)
     &                     -begin(1,IMATS(ISLD2)))
     &     +rbb(2,2,ISLD)*(SFCENT(2,ICFL)
     &                     -begin(2,IMATS(ISLD2)))
     &     +rbb(2,3,ISLD)*(SFCENT(3,ICFL)
     &                     -begin(3,IMATS(ISLD2))) 
          tmpsld(IBFL,3)= !r(3,ISLD2)=
     &      rbb(3,1,ISLD)*(SFCENT(1,ICFL)
     &                     -begin(1,IMATS(ISLD2)))
     &     +rbb(3,2,ISLD)*(SFCENT(2,ICFL)
     &                     -begin(2,IMATS(ISLD2)))
     &     +rbb(3,3,ISLD)*(SFCENT(3,ICFL)
     &                     -begin(3,IMATS(ISLD2)))
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,8)=   !dr(ISLD2)
     &             tmpsld(IBFL,1)*unit(1,ISLD2)
     &            +tmpsld(IBFL,2)*unit(2,ISLD2)
     &            +tmpsld(IBFL,3)*unit(3,ISLD2)          
          enddo

!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,7)=   !vr(3,ISLD2)=
     &                (unit(1,ISLD2)*tmpsld(IBFL,2)
     &                -unit(2,ISLD2)*tmpsld(IBFL,1))
          tmpsld(IBFL,6)=   !vr(2,ISLD2)=
     &                (unit(3,ISLD2)*tmpsld(IBFL,1)
     &                -unit(1,ISLD2)*tmpsld(IBFL,3))
          tmpsld(IBFL,5)=   !vr(1,ISLD2)=
     &                (unit(2,ISLD2)*tmpsld(IBFL,3)
     &                -unit(3,ISLD2)*tmpsld(IBFL,2))
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          dum1=
     &    dsqrt(tmpsld(IBFL,5)**2+tmpsld(IBFL,6)**2+tmpsld(IBFL,7)**2)
          tmpsld(IBFL,5)=tmpsld(IBFL,5)/dum1
          tmpsld(IBFL,6)=tmpsld(IBFL,6)/dum1
          tmpsld(IBFL,7)=tmpsld(IBFL,7)/dum1
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)    !r(1,ISLD2)
     &       =tmpsld(IBFL,1)-tmpsld(IBFL,8)*unit(1,ISLD2)
          tmpsld(IBFL,2)    !r(2,ISLD2)
     &       =tmpsld(IBFL,2)-tmpsld(IBFL,8)*unit(2,ISLD2)
          tmpsld(IBFL,3)    !r(3,ISLD2)
     &       =tmpsld(IBFL,3)-tmpsld(IBFL,8)*unit(3,ISLD2)
          enddo

          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)  !dr(ISLD2)
          tmpsld(IBFL,8)=dsqrt(
     &                    tmpsld(IBFL,1)**2
     &                   +tmpsld(IBFL,2)**2
     &                   +tmpsld(IBFL,3)**2)
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)!       v0(1,ISLD2)
     &      =tmpsld(IBFL,8)*rotati(IMATS(ISLD2))*tmpsld(IBFL,5)
          tmpsld(IBFL,2)!       v0(2,ISLD2)
     &      =tmpsld(IBFL,8)*rotati(IMATS(ISLD2))*tmpsld(IBFL,6)
          tmpsld(IBFL,3)!       v0(3,ISLD2)
     &      =tmpsld(IBFL,8)*rotati(IMATS(ISLD2))*tmpsld(IBFL,7)
          enddo
!
         else
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          tmpsld(IBFL,1)=0.d0
          tmpsld(IBFL,2)=0.d0
          tmpsld(IBFL,3)=0.d0
          enddo
         endif



!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICVP=LVEDGE(1,ICFP)
          ICVBO=LVEDGE(1,ICFO)
          tmpsld(IBFL,5)=  !ufix(1,ISLD2)=
     &                  rbb(1,1,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,1,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,1,k))
     &                 +rbb(1,2,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,2,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,2,k))
     &                 +rbb(1,3,ISLD2)
     &                        *(wifsld(IBFL)*aa(ICVP,3,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,3,k))
     &                 +tmpsld(IBFL,1)!     v0(1,ISLD2)
          tmpsld(IBFL,6)=! ufix(2,ISLD2)
     &                  rbb(2,1,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,1,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,1,k))
     &                 +rbb(2,2,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,2,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,2,k))
     &                 +rbb(2,3,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,3,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,3,k))
     &                 +tmpsld(IBFL,2)!     v0(2,ISLD2)
          tmpsld(IBFL,7)=! ufix(3,ISLD2)
     &                  rbb(3,1,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,1,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,1,k))
     &                 +rbb(3,2,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,2,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,2,k))
     &                 +rbb(3,3,ISLD2)
     &                       *(wifsld(IBFL)*aa(ICVP,3,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,3,k))
     &                 +tmpsld(IBFL,3)!     v0(3,ISLD2)
          enddo
!
         if(imode==0) then
! --- fix2 => rot1
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)=  !r(1,ISLD)=
     &              (SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))

     &            -((SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &              *unit(1,ISLD)
     &             +(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &              *unit(2,ISLD)
     &             +(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &              *unit(3,ISLD))
     &              *unit(1,ISLD)


          tmpsld(IBFL,2)=  !r(2,ISLD)=
     &              (SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))

     &            -((SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &              *unit(1,ISLD)
     &             +(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &              *unit(2,ISLD)
     &             +(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &              *unit(3,ISLD))
     &              *unit(2,ISLD)


          tmpsld(IBFL,3)=  !r(3,ISLD)=
     &              (SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))

     &            -((SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &              *unit(1,ISLD)
     &             +(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &              *unit(2,ISLD)
     &             +(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &              *unit(3,ISLD))
     &              *unit(3,ISLD)
          enddo

          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,4)=dsqrt   !dr(ISLD)=
     &                  (tmpsld(IBFL,1)**2
     &                  +tmpsld(IBFL,2)**2
     &                  +tmpsld(IBFL,3)**2)
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,3)=  ! vr(3,ISLD)=
     &           (unit(1,ISLD)
     &          *(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &           -unit(2,ISLD)
     &          *(SFCENT(1,ICFL)-begin(1,IMATS(ISLD))))

          tmpsld(IBFL,2)=! vr(2,ISLD)=
     &           (unit(3,ISLD)
     &          *(SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &          -unit(1,ISLD)
     &          *(SFCENT(3,ICFL)-begin(3,IMATS(ISLD))))

          tmpsld(IBFL,1)=! vr(1,ISLD)=
     &           (unit(2,ISLD)
     &          *(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &          -unit(3,ISLD)
     &          *(SFCENT(2,ICFL)-begin(2,IMATS(ISLD))))
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          dum1=
     &    dsqrt(tmpsld(IBFL,1)**2+tmpsld(IBFL,2)**2+tmpsld(IBFL,3)**2)
          tmpsld(IBFL,1)=tmpsld(IBFL,1)/dum1
          tmpsld(IBFL,2)=tmpsld(IBFL,2)/dum1
          tmpsld(IBFL,3)=tmpsld(IBFL,3)/dum1
          enddo          
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          
          tmpsld(IBFL,1) !v0(1,ISLD)
     &    =tmpsld(IBFL,4)*rotati(IMATS(ISLD))*tmpsld(IBFL,1)
          tmpsld(IBFL,2) !v0(2,ISLD)
     &    =tmpsld(IBFL,4)*rotati(IMATS(ISLD))*tmpsld(IBFL,2)
          tmpsld(IBFL,3) !v0(3,ISLD)
     &    =tmpsld(IBFL,4)*rotati(IMATS(ISLD))*tmpsld(IBFL,3)
          enddo

         else
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          tmpsld(IBFL,1)=0.d0
          tmpsld(IBFL,2)=0.d0
          tmpsld(IBFL,3)=0.d0
          enddo
         endif


          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
           aa(IDC,1,k)=  !urot(1,ISLD)= 
     &                  rbb(1,1,ISLD)*tmpsld(IBFL,5)
     &                 +rbb(2,1,ISLD)*tmpsld(IBFL,6)
     &                 +rbb(3,1,ISLD)*tmpsld(IBFL,7)
     &                 -tmpsld(IBFL,1) ! v0(1,ISLD)
           aa(IDC,2,k)= !urot(2,ISLD)= 
     &                  rbb(1,2,ISLD)*tmpsld(IBFL,5)
     &                 +rbb(2,2,ISLD)*tmpsld(IBFL,6)
     &                 +rbb(3,2,ISLD)*tmpsld(IBFL,7)
     &                 -tmpsld(IBFL,2) !  v0(2,ISLD)
           aa(IDC,3,k)=  !urot(3,ISLD)= 
     &                  rbb(1,3,ISLD)*tmpsld(IBFL,5)
     &                 +rbb(2,3,ISLD)*tmpsld(IBFL,6)
     &                 +rbb(3,3,ISLD)*tmpsld(IBFL,7)
     &                 -tmpsld(IBFL,3) !  v0(3,ISLD)
!

          ENDDO
!

 1460     enddo
          enddo
        else
!------------------
! --- scalar ------
!------------------
!          call set_rotmtrx(nbcnd,kd,nb,rot)
          do 1200 ISLD=1,2
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
          IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
          th(ISLD)=rot_ang(IMATS(ISLD))
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD)=bb(j,i,ISLD)
          enddo
          enddo
!
          IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
          IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
          th(ISLD2)=rot_ang(IMATS(ISLD2))
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=bb(j,i,ISLD2)
          enddo
          enddo
!
          if(imode==0) then 
            DO IBFL=IBFS,IBFE
            wi1=wifsld(IBFL)
            wi2=1.d0-wi1
!
            ICFL=LBC_SSF(IBFL)
            ICFP=LCYCSF(IBFL)
            ICFO=LCYCOLD(IBFL)
            ICV=LVEDGE(1,ICFL)
            ICVP=LVEDGE(1,ICFP)
            ICVBO=LVEDGE(1,ICFO)
            IDC=LVEDGE(2,ICFL)
!A
            r(1,ISLD)=SFCENT(1,ICFL)-begin(1,IMATS(ISLD))
            r(2,ISLD)=SFCENT(2,ICFL)-begin(2,IMATS(ISLD))
            r(3,ISLD)=SFCENT(3,ICFL)-begin(3,IMATS(ISLD))
!
            call AXB_UNIT_C(unit(:,ISLD),r(:,ISLD),vr(:,ISLD))
!
            dr(ISLD)=r(1,ISLD)*unit(1,ISLD)
     &              +r(2,ISLD)*unit(2,ISLD)
     &              +r(3,ISLD)*unit(3,ISLD)
            r(1,ISLD)=r(1,ISLD)-dr(ISLD)*unit(1,ISLD)
            r(2,ISLD)=r(2,ISLD)-dr(ISLD)*unit(2,ISLD)
            r(3,ISLD)=r(3,ISLD)-dr(ISLD)*unit(3,ISLD)
            dr(ISLD)=dsqrt(r(1,ISLD)*r(1,ISLD)
     &                    +r(2,ISLD)*r(2,ISLD)
     &                    +r(3,ISLD)*r(3,ISLD))
            v0(1,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(1,ISLD)
            v0(2,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(2,ISLD)
            v0(3,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(3,ISLD)
!
            
            dum1=SFCENT(1,ICFL)-begin(1,IMATS(ISLD2)) !7777
            dum2=SFCENT(2,ICFL)-begin(2,IMATS(ISLD2)) 
            dum3=SFCENT(3,ICFL)-begin(3,IMATS(ISLD2)) 
!
            r(1,ISLD2)=rbb(1,1,ISLD)*dum1
     &                +rbb(1,2,ISLD)*dum2
     &                +rbb(1,3,ISLD)*dum3
            r(2,ISLD2)=rbb(2,1,ISLD)*dum1
     &                +rbb(2,2,ISLD)*dum2
     &                +rbb(2,3,ISLD)*dum3
            r(3,ISLD2)=rbb(3,1,ISLD)*dum1
     &                +rbb(3,2,ISLD)*dum2
     &                +rbb(3,3,ISLD)*dum3
!
            call AXB_UNIT_C(unit(:,ISLD2),r(:,ISLD2),vr(:,ISLD2)) 
!
            dr(ISLD2)=r(1,ISLD2)*unit(1,ISLD2)
     &               +r(2,ISLD2)*unit(2,ISLD2)
     &               +r(3,ISLD2)*unit(3,ISLD2)
!
            r(1,ISLD2)=r(1,ISLD2)-dr(ISLD2)*unit(1,ISLD2)
            r(2,ISLD2)=r(2,ISLD2)-dr(ISLD2)*unit(2,ISLD2)
            r(3,ISLD2)=r(3,ISLD2)-dr(ISLD2)*unit(3,ISLD2)
            dr(ISLD2)=dsqrt(r(1,ISLD2)*r(1,ISLD2)
     &                     +r(2,ISLD2)*r(2,ISLD2)
     &                     +r(3,ISLD2)*r(3,ISLD2))
            v0(1,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(1,ISLD2)
            v0(2,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(2,ISLD2)
            v0(3,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(3,ISLD2)
!
            do k=1,nko
!
! --- rot2=> fix2  zhang8:rbb
!
            vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
            vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
            vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
            ufix(1,ISLD2)=
     &            rbb(1,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &           +rbb(1,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &           +rbb(1,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &           +v0(1,ISLD2)
            ufix(2,ISLD2)=
     &                  rbb(2,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(2,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(2,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(2,ISLD2)
            ufix(3,ISLD2)=
     &                  rbb(3,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(3,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(3,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(3,ISLD2)
! --- fix2 => rot1
            urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                   +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                   +rbb(3,1,ISLD)*(ufix(3,ISLD2))
     &                   - v0(1,ISLD)
            urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                   +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                   +rbb(3,2,ISLD)*(ufix(3,ISLD2))
     &                   - v0(2,ISLD)
            urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                   +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                   +rbb(3,3,ISLD)*(ufix(3,ISLD2))
     &                   - v0(3,ISLD)
!
            aa(IDC,1,k)=urot(1,ISLD)
            aa(IDC,2,k)=urot(2,ISLD)
            aa(IDC,3,k)=urot(3,ISLD)
!zhang????
          
            enddo
            enddo
          elseif(imode==1) then
! --- V0=0.d0
            DO IBFL=IBFS,IBFE
            wi1=wifsld(IBFL)
            wi2=1.d0-wi1
            ICFL=LBC_SSF(IBFL)
            IDC=LVEDGE(2,ICFL)
            ICFP=LCYCSF(IBFL)
            ICVP=LVEDGE(1,ICFP)
            ICFO=LCYCOLD(IBFL)
            ICVBO=LVEDGE(1,ICFO)
            do k=1,nko
            vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
            vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
            vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
! --- rot2=> fix2  zhang8:rbb
            ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))
     &                   +rbb(1,2,ISLD2)*(vell(2))
     &                   +rbb(1,3,ISLD2)*(vell(3))
            ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))
     &                   +rbb(2,2,ISLD2)*(vell(2))
     &                   +rbb(2,3,ISLD2)*(vell(3))
            ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))
     &                   +rbb(3,2,ISLD2)*(vell(2))
     &                   +rbb(3,3,ISLD2)*(vell(3))
! --- fix2 => rot1
            urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                   +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                   +rbb(3,1,ISLD)*(ufix(3,ISLD2))
            urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                   +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                   +rbb(3,2,ISLD)*(ufix(3,ISLD2))
            urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                   +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                   +rbb(3,3,ISLD)*(ufix(3,ISLD2))
            aa(IDC,1,k)=urot(1,ISLD)
            aa(IDC,2,k)=urot(2,ISLD)
            aa(IDC,3,k)=urot(3,ISLD)
!zhang????
            enddo
            enddo
          endif
 1200     enddo
        endif
      elseif(kd==kdsld.and.idis(nb)==2) then
        do 1300 ISLD=1,2
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
        IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
        th(ISLD)=rot_ang(IMATS(ISLD))
        unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
        unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
        unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
        CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD)=bb(j,i,ISLD)
        enddo
        enddo
!new
        IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
        IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
!        th(ISLD2)=rot_ang(IMATS(ISLD2))
        unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
        unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
        unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
!        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
!        do i=1,3
!        do j=1,3
!        rbb(i,j,ISLD2)=bb(j,i,ISLD2)
!        enddo
!        enddo
!
        if(imode==0) then
          DO IBFL=IBFS,IBFE
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
!
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)
          ICVBO=LVEDGE(1,ICFO)
          IDC=LVEDGE(2,ICFL)
!
          r(1,ISLD)=SFCENT(1,ICFL)-begin(1,IMATS(ISLD))
          r(2,ISLD)=SFCENT(2,ICFL)-begin(2,IMATS(ISLD))
          r(3,ISLD)=SFCENT(3,ICFL)-begin(3,IMATS(ISLD))
!
          call AXB_UNIT_C(unit(:,ISLD),r(:,ISLD),vr(:,ISLD))
!
          dr(ISLD)=r(1,ISLD)*unit(1,ISLD)
     &            +r(2,ISLD)*unit(2,ISLD)
     &            +r(3,ISLD)*unit(3,ISLD)
          r(1,ISLD)=r(1,ISLD)-dr(ISLD)*unit(1,ISLD)
          r(2,ISLD)=r(2,ISLD)-dr(ISLD)*unit(2,ISLD)
          r(3,ISLD)=r(3,ISLD)-dr(ISLD)*unit(3,ISLD)
          dr(ISLD)=dsqrt(r(1,ISLD)*r(1,ISLD)
     &                  +r(2,ISLD)*r(2,ISLD)
     &                  +r(3,ISLD)*r(3,ISLD))
          v0(1,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(1,ISLD)
          v0(2,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(2,ISLD)
          v0(3,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(3,ISLD)
!new
          th(ISLD2)=OPPANG(IBFL)
          CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=bb(j,i,ISLD2)
          enddo
          enddo
!
          dum1=SFCENT(1,ICFL)-begin(1,IMATS(ISLD2)) !7777
          dum2=SFCENT(2,ICFL)-begin(2,IMATS(ISLD2)) 
          dum3=SFCENT(3,ICFL)-begin(3,IMATS(ISLD2)) 
          r(1,ISLD2)=rbb(1,1,ISLD)*dum1
     &              +rbb(1,2,ISLD)*dum2
     &              +rbb(1,3,ISLD)*dum3
          r(2,ISLD2)=rbb(2,1,ISLD)*dum1
     &              +rbb(2,2,ISLD)*dum2
     &              +rbb(2,3,ISLD)*dum3
          r(3,ISLD2)=rbb(3,1,ISLD)*dum1
     &              +rbb(3,2,ISLD)*dum2
     &              +rbb(3,3,ISLD)*dum3
!
          call AXB_UNIT_C(unit(:,ISLD2),r(:,ISLD2),vr(:,ISLD2))
!
          dr(ISLD2)=r(1,ISLD2)*unit(1,ISLD2)
     &             +r(2,ISLD2)*unit(2,ISLD2)
     &             +r(3,ISLD2)*unit(3,ISLD2)
!
          r(1,ISLD2)=r(1,ISLD2)-dr(ISLD2)*unit(1,ISLD2)
          r(2,ISLD2)=r(2,ISLD2)-dr(ISLD2)*unit(2,ISLD2)
          r(3,ISLD2)=r(3,ISLD2)-dr(ISLD2)*unit(3,ISLD2)
          dr(ISLD2)=dsqrt(r(1,ISLD2)*r(1,ISLD2)
     &                   +r(2,ISLD2)*r(2,ISLD2)
     &                   +r(3,ISLD2)*r(3,ISLD2))
          v0(1,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(1,ISLD2)
          v0(2,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(2,ISLD2)
          v0(3,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(3,ISLD2)
!
          do k=1,nko
!
! --- rot2=> fix2  zhang8:rbb
!
          vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
          vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
          vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
          ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(1,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(1,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(1,ISLD2)
          ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(2,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(2,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(2,ISLD2)
          ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(3,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(3,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(3,ISLD2)
! --- fix2 => rot1
          urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,1,ISLD)*(ufix(3,ISLD2))
     &                 - v0(1,ISLD)
          urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,2,ISLD)*(ufix(3,ISLD2))
     &                 - v0(2,ISLD)
          urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,3,ISLD)*(ufix(3,ISLD2))
     &                 - v0(3,ISLD)
!
          aa(IDC,1,k)=urot(1,ISLD)
          aa(IDC,2,k)=urot(2,ISLD)
          aa(IDC,3,k)=urot(3,ISLD)
          enddo
          enddo
        elseif(imode==1) then
! --- V0=0.d0
          DO IBFL=IBFS,IBFE
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          ICFO=LCYCOLD(IBFL)
          ICVBO=LVEDGE(1,ICFO)
!new
          th(ISLD2)=OPPANG(IBFL)
          CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=bb(j,i,ISLD2)
          enddo
          enddo
!
          do k=1,nko
          vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
          vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
          vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
! --- rot2=> fix2  zhang8:rbb
          ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))
     &                 +rbb(1,2,ISLD2)*(vell(2))
     &                 +rbb(1,3,ISLD2)*(vell(3))
          ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))
     &                 +rbb(2,2,ISLD2)*(vell(2))
     &                 +rbb(2,3,ISLD2)*(vell(3))
          ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))
     &                 +rbb(3,2,ISLD2)*(vell(2))
     &                 +rbb(3,3,ISLD2)*(vell(3))
! --- fix2 => rot1
          urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,1,ISLD)*(ufix(3,ISLD2))
          urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,2,ISLD)*(ufix(3,ISLD2))
          urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,3,ISLD)*(ufix(3,ISLD2))
          aa(IDC,1,k)=urot(1,ISLD)
          aa(IDC,2,k)=urot(2,ISLD)
          aa(IDC,3,k)=urot(3,ISLD)
!zhang????
          enddo
          enddo
        endif
 1300   enddo
        
      elseif(kd.eq.kdsymm) then
!--< 1.2 symmetric boundary >--
	if(ical_vect) then 
	  do k=1,nko
!	  do IBFL=IBFS,IBFE
!          prd=2.d0*(SFAREA(1,LBC_SSF(IBFL))*aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
!     &             +SFAREA(2,LBC_SSF(IBFL))*aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
!     &             +SFAREA(3,LBC_SSF(IBFL))*aa(LVEDGE(1,LBC_SSF(IBFL)),3,k))
!          aa(IDC,1,k)=aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)-prd*SFAREA(1,LBC_SSF(IBFL))
!          aa(IDC,2,k)=aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)-prd*SFAREA(2,LBC_SSF(IBFL))
!          aa(IDC,3,k)=aa(LVEDGE(1,LBC_SSF(IBFL)),3,k)-prd*SFAREA(3,LBC_SSF(IBFL))
!          enddo
          
          do IBFL=IBFS,IBFE
          prd=2.d0*(SFAREA(1,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &             +SFAREA(2,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &             +SFAREA(3,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),3,k))
          tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &                  -prd*SFAREA(1,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
            aa(LVEDGE(2,LBC_SSF(IBFL)),1,k)=tmpfac(IBFL,1)
          enddo
          do IBFL=IBFS,IBFE
            prd=2.d0*(SFAREA(1,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &               +SFAREA(2,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &               +SFAREA(3,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),3,k))
            tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &          -prd*SFAREA(2,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
            aa(LVEDGE(2,LBC_SSF(IBFL)),2,k)=tmpfac(IBFL,1)
          enddo
          do IBFL=IBFS,IBFE
            prd=2.d0*(SFAREA(1,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &               +SFAREA(2,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &               +SFAREA(3,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),3,k))
            tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),3,k)
     &          -prd*SFAREA(3,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
            aa(LVEDGE(2,LBC_SSF(IBFL)),3,k)=tmpfac(IBFL,1)
          enddo
          enddo
	else
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 201 k=1,nko
          prd=2.d0*(SFAREA(1,ICFL)*aa(ICV,1,k)
     &             +SFAREA(2,ICFL)*aa(ICV,2,k)
     &             +SFAREA(3,ICFL)*aa(ICV,3,k))
          aa(IDC,1,k)=aa(ICV,1,k)-prd*SFAREA(1,ICFL)
          aa(IDC,2,k)=aa(ICV,2,k)-prd*SFAREA(2,ICFL)
          aa(IDC,3,k)=aa(ICV,3,k)-prd*SFAREA(3,ICFL)
  201     continue
 1250     continue
 	endif

! --- interface BC
      elseif(kd.eq.kdintr) then
        do 1550 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 501 k=1,nko
        aa(IDC,1,k)=aa(ICV,1,k)
        aa(IDC,2,k)=aa(ICV,2,k)
        aa(IDC,3,k)=aa(ICV,3,k)
!
        aa(IDCP,1,k)=aa(ICVP,1,k)
        aa(IDCP,2,k)=aa(ICVP,2,k)
        aa(IDCP,3,k)=aa(ICVP,3,k)
 501    continue
 1550   continue

!--< 1.3 other boundary >--
      elseif(kd/=kdtchi) then
        if(ical_vect) then 
          do k=1, nko
!NEC$ ivdep
          do IBFL=IBFS,IBFE
          tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
          enddo
!NEC$ ivdep
          do IBFL=IBFS,IBFE
          aa(LVEDGE(2,LBC_SSF(IBFL)),1,k)=tmpfac(IBFL,1)
          enddo
!NEC$ ivdep
          do IBFL=IBFS,IBFE
          tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
          enddo
!NEC$ ivdep
          do IBFL=IBFS,IBFE
          aa(LVEDGE(2,LBC_SSF(IBFL)),2,k)=tmpfac(IBFL,1)
          enddo
!NEC$ ivdep
          do IBFL=IBFS,IBFE
          tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),3,k)
          enddo
!NEC$ ivdep
          do IBFL=IBFS,IBFE
          aa(LVEDGE(2,LBC_SSF(IBFL)),3,k)=tmpfac(IBFL,1)
          enddo
!          
          enddo

        else
          do 1350 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 301 k=1,nko
          aa(IDC,1,k)=aa(ICV,1,k)
          aa(IDC,2,k)=aa(ICV,2,k)
          aa(IDC,3,k)=aa(ICV,3,k)
  301     continue
 1350     continue
        endif
      endif
 1000 continue
!
! --- --< 1.4 touch-inlet BC >--
!
!      do 2000 nb=1,nbcnd
!      IIMAT=MAT_BCIDX(nb,1)
!      if(.not.mat_cal(IIMAT)) goto 2000
!      kd=kdbcnd(0,nb)
!      if(kd.eq.kdtchi) then
!        call set_rotmtrx(nbcnd,kd,nb,rot)
!        IBFS=LBC_INDEX(nb-1)+1
!        IBFE=LBC_INDEX(nb)
!        do 1450 IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        ICV=LVEDGE(1,ICFL)
!        IDC=LVEDGE(2,ICFL)
!        ICFP=LCYCSF(IBFL)
!        ICVP=LVEDGE(1,ICFP)
!        IDCP=LVEDGE(2,ICFP)
!        do 401 k=1,nko
!        aa(IDC,1,k)= rot(1,1)*aa(ICVP,1,k)
!     &              +rot(1,2)*aa(ICVP,2,k)
!     &              +rot(1,3)*aa(ICVP,3,k)
!        aa(IDC,2,k)= rot(2,1)*aa(ICVP,1,k)
!     &              +rot(2,2)*aa(ICVP,2,k)
!     &              +rot(2,3)*aa(ICVP,3,k)
!        aa(IDC,3,k)= rot(3,1)*aa(ICVP,1,k)
!     &              +rot(3,2)*aa(ICVP,2,k)
!     &              +rot(3,3)*aa(ICVP,3,k)
! 401    continue
! 1450   continue
!      endif
! 2000 continue
!
      return
      end subroutine dc_symprv
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_symvel
     &(nko,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & LCYCOLD,wifsld,OPPANG,
     & SFAREA,SFCENT,aa,imode,ivel)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     imode==0 : velocity
!     imode==1 : grdc
!
! --- [module arguments]
!
      USE module_dimension
      use module_hpcutil,only  : my_rank
      use module_boundary,only : set_rotmtrx,nbcnd,LBC_INDEX,
     &                           kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdintr,kdsld,MAT_BCIDX,idis,
     &                           LBC_pair,rotsld,kdbuff
     &                           ,kdovst
      use module_material,only : rotati,ishaft,end,begin,nsplit,
     &                           rot_ang
      use module_nowtime,only  : iter,time
      use module_material,only : ical_sld
      use module_metrix,only   : tmpsld
      use module_model,only    : ical_vect
!
! 1. Set vector at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko,imode,ivel
      integer,intent(in)    :: MAT_NO    (  0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(inout) :: aa        (  MXALLCV,3,nko)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
!
! --- [local entities]
!
      real*8  :: prd
      real*8  :: shft(3),ufix(3,2)
      real*8  :: unit(3,2),th(2),bb(3,3,2),vell(3)
      real*8  :: rbb(3,3,2)
      real*8  :: urot(3,2)
      real*8  :: vr(3,2),r(3,2)
      real*8  :: v0(3,2),dr(2),wi1,wi2,dum1,dum2,dum3
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2)
      integer :: i,j,ICVBO,ICFO,ibfs1,ibfe1,ibfs2,ibfe2
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k,ISLD,ISLD2
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
! --- 
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      if(kd==kdovst) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        IIMATS(1)=MAT_BCIDX(nb,1)
        IIMATS(2)=MAT_BCIDX(nb,2)
        IMATS(1)=MAT_NO(IIMATS(1))
        IMATS(2)=MAT_NO(IIMATS(2))

        th(1)=rot_ang(IMATS(1))
        th(2)=rot_ang(IMATS(2))

        unit(1,2)=end(1,IMATS(2))-begin(1,IMATS(2))
        unit(2,2)=end(2,IMATS(2))-begin(2,IMATS(2))
        unit(3,2)=end(3,IMATS(2))-begin(3,IMATS(2))
        CALL rotth(unit(:,2),th(2),bb(:,:,2))

        unit(1,1)=end(1,IMATS(1))-begin(1,IMATS(1))
        unit(2,1)=end(2,IMATS(1))-begin(2,IMATS(1))
        unit(3,1)=end(3,IMATS(1))-begin(3,IMATS(1))
        CALL rotth(unit(:,1),th(1),bb(:,:,1))

        do i=1,3
        do j=1,3
        rbb(i,j,1)=bb(j,i,1)
        rbb(i,j,2)=bb(j,i,2)
        enddo
        enddo
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LCYCSF(IBFL)
        r(1,1)=SFCENT(1,ICFL)-begin(1,IMATS(1))
        r(2,1)=SFCENT(2,ICFL)-begin(2,IMATS(1))
        r(3,1)=SFCENT(3,ICFL)-begin(3,IMATS(1))
        call AXB_UNIT_C(unit(:,1),r(:,1),vr(:,1))
        dr(1)=r(1,1)*unit(1,1)
     &          +r(2,1)*unit(2,1)
     &          +r(3,1)*unit(3,1)
        r(1,1)=r(1,1)-dr(1)*unit(1,1)
        r(2,1)=r(2,1)-dr(1)*unit(2,1)
        r(3,1)=r(3,1)-dr(1)*unit(3,1)
        dr(1)=dsqrt(r(1,1)*r(1,1)
     &             +r(2,1)*r(2,1)
     &             +r(3,1)*r(3,1))
        v0(1,1)=dr(1)*rotati(IMATS(1))*vr(1,1)
        v0(2,1)=dr(1)*rotati(IMATS(1))*vr(2,1)
        v0(3,1)=dr(1)*rotati(IMATS(1))*vr(3,1)
!
        dum1=SFCENT(1,ICFL)-begin(1,IMATS(2)) !7777
        dum2=SFCENT(2,ICFL)-begin(2,IMATS(2)) 
        dum3=SFCENT(3,ICFL)-begin(3,IMATS(2)) 
!
        r(1,2)=rbb(1,1,1)*dum1
     &        +rbb(1,2,1)*dum2
     &        +rbb(1,3,1)*dum3
        r(2,2)=rbb(2,1,1)*dum1
     &        +rbb(2,2,1)*dum2
     &        +rbb(2,3,1)*dum3
        r(3,2)=rbb(3,1,1)*dum1
     &        +rbb(3,2,1)*dum2
     &        +rbb(3,3,1)*dum3
!
        call AXB_UNIT_C(unit(:,2),r(:,2),vr(:,2)) 
!
        dr(2)=r(1,2)*unit(1,2)
     &       +r(2,2)*unit(2,2)
     &       +r(3,2)*unit(3,2)
!
        r(1,2)=r(1,2)-dr(2)*unit(1,2)
        r(2,2)=r(2,2)-dr(2)*unit(2,2)
        r(3,2)=r(3,2)-dr(2)*unit(3,2)
        dr(2)=dsqrt(r(1,2)*r(1,2)
     &             +r(2,2)*r(2,2)
     &             +r(3,2)*r(3,2))
        v0(1,2)=dr(2)*rotati(IMATS(2))*vr(1,2)
        v0(2,2)=dr(2)*rotati(IMATS(2))*vr(2,2)
        v0(3,2)=dr(2)*rotati(IMATS(2))*vr(3,2)        
!
        IF(imode==1) THEN
          v0(:,:)=0.D0
        ENDIF

        do k=1,nko
        vell(1)=aa(ICVP,1,k)
        vell(2)=aa(ICVP,2,k)
        vell(3)=aa(ICVP,3,k)
! --- rot2=> fix2
        ufix(1,2)=
     &            rbb(1,1,2)*(vell(1))!+v0(1,2))
     &           +rbb(1,2,2)*(vell(2))!+v0(2,2))
     &           +rbb(1,3,2)*(vell(3))!+v0(3,2))
     &           +v0(1,2)
        ufix(2,2)=
     &            rbb(2,1,2)*(vell(1))!+v0(1,2))
     &           +rbb(2,2,2)*(vell(2))!+v0(2,2))
     &           +rbb(2,3,2)*(vell(3))!+v0(3,2))
     &           +v0(2,2)
        ufix(3,2)=
     &            rbb(3,1,2)*(vell(1))!+v0(1,2))
     &           +rbb(3,2,2)*(vell(2))!+v0(2,2))
     &           +rbb(3,3,2)*(vell(3))!+v0(3,2))
     &           +v0(3,2)
! --- fix2 => rot1
        urot(1,1)= rbb(1,1,1)*(ufix(1,2))
     &            +rbb(2,1,1)*(ufix(2,2))
     &            +rbb(3,1,1)*(ufix(3,2))
     &            -v0(1,1)
        urot(2,1)= rbb(1,2,1)*(ufix(1,2))
     &            +rbb(2,2,1)*(ufix(2,2))
     &            +rbb(3,2,1)*(ufix(3,2))
     &            -v0(2,1)
        urot(3,1)= rbb(1,3,1)*(ufix(1,2))
     &            +rbb(2,3,1)*(ufix(2,2))
     &            +rbb(3,3,1)*(ufix(3,2))
     &            -v0(3,1)

        aa(IDC,1,k)=urot(1,1)
        aa(IDC,2,k)=urot(2,1)
        aa(IDC,3,k)=urot(3,1)
        enddo
        
        enddo
        
      elseif(kd.eq.kdsld.and.idis(nb)==1) then
       if(ical_vect) then 
          do k=1,nko
          do 1450 ISLD=1,2
          ISLD2=3-ISLD
          IBFS1=(LBC_INDEX(nb-1)+1)*(2-ISLD)
          IBFE1=(LBC_pair(nb))*(2-ISLD)
          IBFS2=(LBC_pair(nb)+1)*(2-ISLD2)
          IBFE2=(LBC_INDEX(nb))*(2-ISLD2)
          IMATS(ISLD)=MAT_NO(MAT_BCIDX(nb,ISLD))
          th(ISLD)=rot_ang(MAT_NO(MAT_BCIDX(nb,ISLD)))
          unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
          unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
          unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
          CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
          rbb(1,1,ISLD)=bb(1,1,ISLD)
          rbb(1,2,ISLD)=bb(2,1,ISLD)
          rbb(1,3,ISLD)=bb(3,1,ISLD)
          rbb(2,1,ISLD)=bb(1,2,ISLD)
          rbb(2,2,ISLD)=bb(2,2,ISLD)
          rbb(2,3,ISLD)=bb(3,2,ISLD)
          rbb(3,1,ISLD)=bb(1,3,ISLD)
          rbb(3,2,ISLD)=bb(2,3,ISLD)
          rbb(3,3,ISLD)=bb(3,3,ISLD)

          IMATS(ISLD2)=MAT_NO(MAT_BCIDX(nb,ISLD2))
          th(ISLD2)=rot_ang(IMATS(ISLD2))
          unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
          unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
          unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
          CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
          rbb(1,1,ISLD2)=bb(1,1,ISLD2)
          rbb(1,2,ISLD2)=bb(2,1,ISLD2)
          rbb(1,3,ISLD2)=bb(3,1,ISLD2)
          rbb(2,1,ISLD2)=bb(1,2,ISLD2)
          rbb(2,2,ISLD2)=bb(2,2,ISLD2)
          rbb(2,3,ISLD2)=bb(3,2,ISLD2)
          rbb(3,1,ISLD2)=bb(1,3,ISLD2)
          rbb(3,2,ISLD2)=bb(2,3,ISLD2)
          rbb(3,3,ISLD2)=bb(3,3,ISLD2)
! --- imode
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)= !r(1,ISLD2)=
     &               rbb(1,1,ISLD)*(SFCENT(1,ICFL)
     &              -begin(1,IMATS(ISLD2)))
     &              +rbb(1,2,ISLD)*(SFCENT(2,ICFL)
     &              -begin(2,IMATS(ISLD2)))
     &              +rbb(1,3,ISLD)*(SFCENT(3,ICFL)
     &              -begin(3,IMATS(ISLD2))) 
          tmpsld(IBFL,2)= !r(2,ISLD2)=
     &               rbb(2,1,ISLD)*(SFCENT(1,ICFL)
     &              -begin(1,IMATS(ISLD2)))
     &              +rbb(2,2,ISLD)*(SFCENT(2,ICFL)
     &              -begin(2,IMATS(ISLD2)))
     &              +rbb(2,3,ISLD)*(SFCENT(3,ICFL)
     &              -begin(3,IMATS(ISLD2))) 
          tmpsld(IBFL,3)= !r(3,ISLD2)=
     &               rbb(3,1,ISLD)*(SFCENT(1,ICFL)
     &              -begin(1,IMATS(ISLD2)))
     &              +rbb(3,2,ISLD)*(SFCENT(2,ICFL)
     &              -begin(2,IMATS(ISLD2)))
     &              +rbb(3,3,ISLD)*(SFCENT(3,ICFL)
     &              -begin(3,IMATS(ISLD2)))
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,8)=   !dr(ISLD2)
     &             tmpsld(IBFL,1)*unit(1,ISLD2)
     &            +tmpsld(IBFL,2)*unit(2,ISLD2)
     &            +tmpsld(IBFL,3)*unit(3,ISLD2)          
          enddo

!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,7)=   !vr(3,ISLD2)=
     &                (unit(1,ISLD2)*tmpsld(IBFL,2)
     &                -unit(2,ISLD2)*tmpsld(IBFL,1))
          tmpsld(IBFL,6)=   !vr(2,ISLD2)=
     &                (unit(3,ISLD2)*tmpsld(IBFL,1)
     &                -unit(1,ISLD2)*tmpsld(IBFL,3))
          tmpsld(IBFL,5)=   !vr(1,ISLD2)=
     &                (unit(2,ISLD2)*tmpsld(IBFL,3)
     &                -unit(3,ISLD2)*tmpsld(IBFL,2))
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          dum1=
     &    dsqrt(tmpsld(IBFL,5)**2+tmpsld(IBFL,6)**2+tmpsld(IBFL,7)**2)
          tmpsld(IBFL,5)=tmpsld(IBFL,5)/dum1
          tmpsld(IBFL,6)=tmpsld(IBFL,6)/dum1
          tmpsld(IBFL,7)=tmpsld(IBFL,7)/dum1
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)    !r(1,ISLD2)
     &       =tmpsld(IBFL,1)-tmpsld(IBFL,8)*unit(1,ISLD2)
          tmpsld(IBFL,1)    !r(2,ISLD2)
     &       =tmpsld(IBFL,2)-tmpsld(IBFL,8)*unit(2,ISLD2)
          tmpsld(IBFL,1)    !r(3,ISLD2)
     &       =tmpsld(IBFL,3)-tmpsld(IBFL,8)*unit(3,ISLD2)
          enddo

          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,8)=dsqrt(
     &                    tmpsld(IBFL,1)**2
     &                   +tmpsld(IBFL,2)**2
     &                   +tmpsld(IBFL,3)**2)
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)!       v0(1,ISLD2)
     &               =tmpsld(IBFL,8)*rotati(IMATS(ISLD2))*tmpsld(IBFL,5)
          tmpsld(IBFL,2)!       v0(2,ISLD2)
     &               =tmpsld(IBFL,8)*rotati(IMATS(ISLD2))*tmpsld(IBFL,6)
          tmpsld(IBFL,3)!       v0(3,ISLD2)
     &               =tmpsld(IBFL,8)*rotati(IMATS(ISLD2))*tmpsld(IBFL,7)
          enddo
!

!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICVP=LVEDGE(1,ICFP)
          ICVBO=LVEDGE(1,ICFO)
          tmpsld(IBFL,5)=
     &                  rbb(1,1,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,1,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,1,k))
     &                 +rbb(1,2,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,2,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,2,k))
     &                 +rbb(1,3,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,3,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,3,k))
     &                 +tmpsld(IBFL,1)!     v0(1,ISLD2)
          tmpsld(IBFL,6)=
     &                  rbb(2,1,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,1,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,1,k))
     &                 +rbb(2,2,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,2,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,2,k))
     &                 +rbb(2,3,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,3,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,3,k))
     &                 +tmpsld(IBFL,2)!     v0(2,ISLD2)
          tmpsld(IBFL,7)=
     &                  rbb(3,1,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,1,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,1,k))
     &                 +rbb(3,2,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,2,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,2,k))
     &                 +rbb(3,3,ISLD2)
     &                 *(wifsld(IBFL)*aa(ICVP,3,k)
     &                  +(1.d0-wifsld(IBFL))*aa(ICVBO,3,k))
     &                 +tmpsld(IBFL,3)!     v0(3,ISLD2)
          enddo

! --- fix2 => rot1
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,1)=  !r(1,ISLD)=
     &              (SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &            -((SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &              *unit(1,ISLD)
     &             +(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &              *unit(2,ISLD)
     &             +(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &              *unit(3,ISLD))
     &              *unit(1,ISLD)
          tmpsld(IBFL,2)=  !r(2,ISLD)=
     &              (SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &            -((SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &              *unit(1,ISLD)
     &             +(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &              *unit(2,ISLD)
     &             +(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &              *unit(3,ISLD))
     &              *unit(2,ISLD)
          tmpsld(IBFL,3)=  !r(3,ISLD)=
     &              (SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &            -((SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &              *unit(1,ISLD)
     &             +(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &              *unit(2,ISLD)
     &             +(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &              *unit(3,ISLD))
     &              *unit(3,ISLD)
          enddo

          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,4)=dsqrt   !dr(ISLD)=
     &                  (tmpsld(IBFL,1)**2
     &                  +tmpsld(IBFL,2)**2
     &                  +tmpsld(IBFL,3)**2)
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          tmpsld(IBFL,3)=  ! vr(3,ISLD)=
     &           (unit(1,ISLD)
     &          *(SFCENT(2,ICFL)-begin(2,IMATS(ISLD)))
     &           -unit(2,ISLD)
     &          *(SFCENT(1,ICFL)-begin(1,IMATS(ISLD))))

          tmpsld(IBFL,2)=! vr(2,ISLD)=
     &           (unit(3,ISLD)
     &          *(SFCENT(1,ICFL)-begin(1,IMATS(ISLD)))
     &          -unit(1,ISLD)
     &          *(SFCENT(3,ICFL)-begin(3,IMATS(ISLD))))

          tmpsld(IBFL,1)=! vr(1,ISLD)=
     &           (unit(2,ISLD)
     &          *(SFCENT(3,ICFL)-begin(3,IMATS(ISLD)))
     &          -unit(3,ISLD)
     &          *(SFCENT(2,ICFL)-begin(2,IMATS(ISLD))))
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          dum1=
     &    dsqrt(tmpsld(IBFL,1)**2+tmpsld(IBFL,2)**2+tmpsld(IBFL,3)**2)
          tmpsld(IBFL,1)=tmpsld(IBFL,1)/dum1
          tmpsld(IBFL,2)=tmpsld(IBFL,2)/dum1
          tmpsld(IBFL,3)=tmpsld(IBFL,3)/dum1
          enddo     
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          
          tmpsld(IBFL,1) !v0(1,ISLD)
     &       =tmpsld(IBFL,4)*rotati(IMATS(ISLD))*tmpsld(IBFL,1)
          tmpsld(IBFL,2) !v0(2,ISLD)
     &       =tmpsld(IBFL,4)*rotati(IMATS(ISLD))*tmpsld(IBFL,2)
          tmpsld(IBFL,3) !v0(3,ISLD)
     &       =tmpsld(IBFL,4)*rotati(IMATS(ISLD))*tmpsld(IBFL,3)
          enddo
!
          DO IBFL=IBFS1+IBFS2,IBFE1+IBFE2
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
           aa(IDC,1,k)=  !urot(1,ISLD)= 
     &                  rbb(1,1,ISLD)*tmpsld(IBFL,5)
     &                 +rbb(2,1,ISLD)*tmpsld(IBFL,6)
     &                 +rbb(3,1,ISLD)*tmpsld(IBFL,7)
     &                 -tmpsld(IBFL,1) ! v0(1,ISLD)
           aa(IDC,2,k)= !urot(2,ISLD)= 
     &                  rbb(1,2,ISLD)*tmpsld(IBFL,5)
     &                 +rbb(2,2,ISLD)*tmpsld(IBFL,6)
     &                 +rbb(3,2,ISLD)*tmpsld(IBFL,7)
     &                 -tmpsld(IBFL,2) !  v0(2,ISLD)
           aa(IDC,3,k)=  !urot(3,ISLD)= 
     &                  rbb(1,3,ISLD)*tmpsld(IBFL,5)
     &                 +rbb(2,3,ISLD)*tmpsld(IBFL,6)
     &                 +rbb(3,3,ISLD)*tmpsld(IBFL,7)
     &                 -tmpsld(IBFL,3) !  v0(3,ISLD)
!
          ENDDO
!

 1450     enddo
          enddo
        
       else
        do 1200 ISLD=1,2
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
!
        IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
        IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
        th(ISLD)=rot_ang(IMATS(ISLD))
        unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
        unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
        unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
        CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD)=bb(j,i,ISLD)
        enddo
        enddo
!
        IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
        IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
        th(ISLD2)=rot_ang(IMATS(ISLD2))
        unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
        unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
        unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD2)=bb(j,i,ISLD2)
        enddo
        enddo
!
        DO IBFL=IBFS,IBFE
!
        wi1=wifsld(IBFL)
        wi2=1.d0-wi1
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        ICFO=LCYCOLD(IBFL)

        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LVEDGE(1,ICFP)
        ICVBO=LVEDGE(1,ICFO)

        r(1,ISLD)=SFCENT(1,ICFL)-begin(1,IMATS(ISLD))
        r(2,ISLD)=SFCENT(2,ICFL)-begin(2,IMATS(ISLD))
        r(3,ISLD)=SFCENT(3,ICFL)-begin(3,IMATS(ISLD))
!
        call AXB_UNIT_C(unit(:,ISLD),r(:,ISLD),vr(:,ISLD))
!
        dr(ISLD)=r(1,ISLD)*unit(1,ISLD)
     &          +r(2,ISLD)*unit(2,ISLD)
     &          +r(3,ISLD)*unit(3,ISLD)
        r(1,ISLD)=r(1,ISLD)-dr(ISLD)*unit(1,ISLD)
        r(2,ISLD)=r(2,ISLD)-dr(ISLD)*unit(2,ISLD)
        r(3,ISLD)=r(3,ISLD)-dr(ISLD)*unit(3,ISLD)
        dr(ISLD)=dsqrt(r(1,ISLD)*r(1,ISLD)
     &                +r(2,ISLD)*r(2,ISLD)
     &                +r(3,ISLD)*r(3,ISLD))
        v0(1,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(1,ISLD)
        v0(2,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(2,ISLD)
        v0(3,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(3,ISLD)

        dum1=SFCENT(1,ICFL)-begin(1,IMATS(ISLD2)) !7777
        dum2=SFCENT(2,ICFL)-begin(2,IMATS(ISLD2)) 
        dum3=SFCENT(3,ICFL)-begin(3,IMATS(ISLD2))
        r(1,ISLD2)=rbb(1,1,ISLD)*dum1
     &            +rbb(1,2,ISLD)*dum2
     &            +rbb(1,3,ISLD)*dum3
        r(2,ISLD2)=rbb(2,1,ISLD)*dum1
     &            +rbb(2,2,ISLD)*dum2
     &            +rbb(2,3,ISLD)*dum3
        r(3,ISLD2)=rbb(3,1,ISLD)*dum1
     &            +rbb(3,2,ISLD)*dum2
     &            +rbb(3,3,ISLD)*dum3
!
        call AXB_UNIT_C(unit(:,ISLD2),r(:,ISLD2),vr(:,ISLD2))
!
        dr(ISLD2)=r(1,ISLD2)*unit(1,ISLD2)
     &           +r(2,ISLD2)*unit(2,ISLD2)
     &           +r(3,ISLD2)*unit(3,ISLD2)
        r(1,ISLD2)=r(1,ISLD2)-dr(ISLD2)*unit(1,ISLD2)
        r(2,ISLD2)=r(2,ISLD2)-dr(ISLD2)*unit(2,ISLD2)
        r(3,ISLD2)=r(3,ISLD2)-dr(ISLD2)*unit(3,ISLD2)
        dr(ISLD2)=dsqrt(r(1,ISLD2)*r(1,ISLD2)
     &                 +r(2,ISLD2)*r(2,ISLD2)
     &                 +r(3,ISLD2)*r(3,ISLD2))
        v0(1,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(1,ISLD2)
        v0(2,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(2,ISLD2)
        v0(3,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(3,ISLD2)
!
        do k=1,nko
        vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
        vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
        vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
! --- rot2=> fix2  zhang8:rbb
        ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &               +rbb(1,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &               +rbb(1,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &               +v0(1,ISLD2)
        ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &               +rbb(2,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &               +rbb(2,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &               +v0(2,ISLD2)
        ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &               +rbb(3,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &               +rbb(3,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &               +v0(3,ISLD2)
! --- fix2 => rot1
        urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &               +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &               +rbb(3,1,ISLD)*(ufix(3,ISLD2))
     &               - v0(1,ISLD)
        urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &               +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &               +rbb(3,2,ISLD)*(ufix(3,ISLD2))
     &               - v0(2,ISLD)
        urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &               +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &               +rbb(3,3,ISLD)*(ufix(3,ISLD2))
     &               - v0(3,ISLD)
        aa(IDC,1,k)=urot(1,ISLD)
        aa(IDC,2,k)=urot(2,ISLD)
        aa(IDC,3,k)=urot(3,ISLD)
        enddo
        enddo
 1200   enddo
       endif
      elseif(kd.eq.kdsld.and.idis(nb)==2) then
        do 1300 ISLD=1,2
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
!
        IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
        IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
        th(ISLD)=rot_ang(IMATS(ISLD))
        unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
        unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
        unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
        CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD)=bb(j,i,ISLD)
        enddo
        enddo
!new
        IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
        IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
!        th(ISLD2)=rot_ang(IMATS(ISLD2))
        unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
        unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
        unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
!        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
!        do i=1,3
!        do j=1,3
!        rbb(i,j,ISLD2)=bb(j,i,ISLD2)
!        enddo
!        enddo
!
        DO IBFL=IBFS,IBFE
!
        wi1=wifsld(IBFL)
        wi2=1.d0-wi1
        ICFL=LBC_SSF(IBFL)
        ICFP=LCYCSF(IBFL)
        ICFO=LCYCOLD(IBFL)

        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICVP=LVEDGE(1,ICFP)
        ICVBO=LVEDGE(1,ICFO)

        r(1,ISLD)=SFCENT(1,ICFL)-begin(1,IMATS(ISLD))
        r(2,ISLD)=SFCENT(2,ICFL)-begin(2,IMATS(ISLD))
        r(3,ISLD)=SFCENT(3,ICFL)-begin(3,IMATS(ISLD))
!
        call AXB_UNIT_C(unit(:,ISLD),r(:,ISLD),vr(:,ISLD))
!
        dr(ISLD)=r(1,ISLD)*unit(1,ISLD)
     &          +r(2,ISLD)*unit(2,ISLD)
     &          +r(3,ISLD)*unit(3,ISLD)
        r(1,ISLD)=r(1,ISLD)-dr(ISLD)*unit(1,ISLD)
        r(2,ISLD)=r(2,ISLD)-dr(ISLD)*unit(2,ISLD)
        r(3,ISLD)=r(3,ISLD)-dr(ISLD)*unit(3,ISLD)
        dr(ISLD)=dsqrt(r(1,ISLD)*r(1,ISLD)
     &                +r(2,ISLD)*r(2,ISLD)
     &                +r(3,ISLD)*r(3,ISLD))
        v0(1,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(1,ISLD)
        v0(2,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(2,ISLD)
        v0(3,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(3,ISLD)
!new
        th(ISLD2)=OPPANG(IBFL)
        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD2)=bb(j,i,ISLD2)
        enddo
        enddo
!
        dum1=SFCENT(1,ICFL)-begin(1,IMATS(ISLD2)) !7777
        dum2=SFCENT(2,ICFL)-begin(2,IMATS(ISLD2)) 
        dum3=SFCENT(3,ICFL)-begin(3,IMATS(ISLD2)) 
        r(1,ISLD2)=rbb(1,1,ISLD)*dum1
     &            +rbb(1,2,ISLD)*dum2
     &            +rbb(1,3,ISLD)*dum3
        r(2,ISLD2)=rbb(2,1,ISLD)*dum1
     &            +rbb(2,2,ISLD)*dum2
     &            +rbb(2,3,ISLD)*dum3
        r(3,ISLD2)=rbb(3,1,ISLD)*dum1
     &            +rbb(3,2,ISLD)*dum2
     &            +rbb(3,3,ISLD)*dum3
!
        call AXB_UNIT_C(unit(:,ISLD2),r(:,ISLD2),vr(:,ISLD2))
!
        dr(ISLD2)=r(1,ISLD2)*unit(1,ISLD2)
     &           +r(2,ISLD2)*unit(2,ISLD2)
     &           +r(3,ISLD2)*unit(3,ISLD2)
        r(1,ISLD2)=r(1,ISLD2)-dr(ISLD2)*unit(1,ISLD2)
        r(2,ISLD2)=r(2,ISLD2)-dr(ISLD2)*unit(2,ISLD2)
        r(3,ISLD2)=r(3,ISLD2)-dr(ISLD2)*unit(3,ISLD2)
        dr(ISLD2)=dsqrt(r(1,ISLD2)*r(1,ISLD2)
     &                 +r(2,ISLD2)*r(2,ISLD2)
     &                 +r(3,ISLD2)*r(3,ISLD2))
        v0(1,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(1,ISLD2)
        v0(2,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(2,ISLD2)
        v0(3,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(3,ISLD2)
!
        do k=1,nko
        vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
        vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
        vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
! --- rot2=> fix2  zhang8:rbb
        ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1)+v0(1,ISLD2))
     &               +rbb(1,2,ISLD2)*(vell(2)+v0(2,ISLD2))
     &               +rbb(1,3,ISLD2)*(vell(3)+v0(3,ISLD2))
!     &               +v0(1,ISLD2)
        ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1)+v0(1,ISLD2))
     &               +rbb(2,2,ISLD2)*(vell(2)+v0(2,ISLD2))
     &               +rbb(2,3,ISLD2)*(vell(3)+v0(3,ISLD2))
!     &               +v0(2,ISLD2)
        ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1)+v0(1,ISLD2))
     &               +rbb(3,2,ISLD2)*(vell(2)+v0(2,ISLD2))
     &               +rbb(3,3,ISLD2)*(vell(3)+v0(3,ISLD2))
!     &               +v0(3,ISLD2)
! --- fix2 => rot1
        urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &               +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &               +rbb(3,1,ISLD)*(ufix(3,ISLD2))
     &               - v0(1,ISLD)
        urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &               +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &               +rbb(3,2,ISLD)*(ufix(3,ISLD2))
     &               - v0(2,ISLD)
        urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &               +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &               +rbb(3,3,ISLD)*(ufix(3,ISLD2))
     &               - v0(3,ISLD)
        aa(IDC,1,k)=urot(1,ISLD)
        aa(IDC,2,k)=urot(2,ISLD)
        aa(IDC,3,k)=urot(3,ISLD)
!
        enddo
        enddo
 1300   enddo
      elseif(kd==kdbuff) then
        cycle
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)

        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do k=1,nko
!
        aa(IDC,1,k)=aa(ICV,1,k)
        aa(IDC,2,k)=aa(ICV,2,k)
        aa(IDC,3,k)=aa(ICV,3,k)

        aa(IDCP,1,k)=aa(ICVP,1,k)
        aa(IDCP,2,k)=aa(ICVP,2,k)
        aa(IDCP,3,k)=aa(ICVP,3,k)
!
        enddo
        enddo
      endif
 1000 continue
!
      return
      end subroutine dc_symvel
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_vgrad
     &  (MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,
     &   LCYCOLD,wifsld,OPPANG,vv)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_boundary,only : set_rotmtrx,nbcnd,LBC_INDEX,
     &                           kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdintr,kdsld,MAT_BCIDX,rotsld,
     &                           idis,LBC_pair,kdbuff,kdshutr
     &                           ,kdovst
      use module_material,only : rotati,ishaft,end,begin,ical_sld,
     &                           rot_ang
      use module_nowtime,only  : iter,time
      use module_model,only    : ical_vect
      use module_metrix,only   : tmpfac=>d2vect
      use module_metrix,only   : SHUTFL
      use module_metrix,only   : tmpsld
!
! 1. Set gradient of vector at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      INTEGER,INTENT(IN)    :: MAT_NO    (  0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(inout) :: vv        ( MXALLCV,3,3)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
!
! --- [local entities]
!
      real*8  :: x1(3,3),x2(3,3),y(3),xx(3,3),aa(3,3)
      real*8  :: g1(3,3),g2(3,3),wi1,wi2
      real*8  :: prd,unit(3,2),th(2),bb(3,3,2),rbb(3,3,2),vell(3,3)
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2),ICVBO,ICFO,ierr1=0
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,i,j,ISLD,ISLD2
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ibfs1,ibfe1,ibfs2,ibfe2
!
!-< 1. Set values at dummy cell >-
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdprdc) then
        call set_rotmtrx(nbcnd,kd,nb,aa)
!--< 1.1 periodic boundary >--
        if(idis(nb)==0) then
          do 1150 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do 101 j=1,3
          do 101 i=1,3
          x1(i,j)=aa(i,1)*vv(ICV,1,j)   !=>A^-1*TA
     &           +aa(i,2)*vv(ICV,2,j)
     &           +aa(i,3)*vv(ICV,3,j)
          x2(i,j)=aa(1,i)*vv(ICVP,1,j)  !=>A*TB
     &           +aa(2,i)*vv(ICVP,2,j)
     &           +aa(3,i)*vv(ICVP,3,j)
 101      continue
          do 102 j=1,3
          do 102 i=1,3
          vv(IDCP,i,j)=x1(i,1)*aa(1,j)  !=>A^-1*TA*A
     &                +x1(i,2)*aa(2,j)
     &                +x1(i,3)*aa(3,j)
          vv(IDC,i,j)= x2(i,1)*aa(j,1)  !=>A*TB*A^-1
     &                +x2(i,2)*aa(j,2)
     &                +x2(i,3)*aa(j,3)
 102      continue
 1150     continue
        elseif(idis(nb)==1) then !zhang???
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          do j=1,3
          do i=1,3
!          x1(i,j)=aa(i,1)*vv(ICV,1,j)   !=>A^-1*TA
!     &         +aa(i,2)*vv(ICV,2,j)
!     &         +aa(i,3)*vv(ICV,3,j)
          x2(i,j)=aa(1,i)*vv(ICVP,1,j)  !=>A*TB
     &         +aa(2,i)*vv(ICVP,2,j)
     &         +aa(3,i)*vv(ICVP,3,j)
          enddo
          enddo
!
          do j=1,3
          do i=1,3
!          vv(IDCP,i,j)=x1(i,1)*aa(1,j)  !=>A^-1*TA*A
!     &              +x1(i,2)*aa(2,j)
!     &              +x1(i,3)*aa(3,j)
          vv(IDC,i,j)= x2(i,1)*aa(j,1)  !=>A*TB*A^-1
     &              +x2(i,2)*aa(j,2)
     &              +x2(i,3)*aa(j,3)
          enddo
          enddo
          enddo
        endif
      elseif(kd.eq.kdbuff) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do j=1,3
          do i=1,3
          vv(IDCP,i,j)=vv(ICV,i,j)
          vv(IDC,i,j)=vv(ICVP,i,j)
          enddo
          enddo
        enddo
      elseif(kd.eq.kdshutr) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          if(SHUTFL(IBFL)==0) then
            do j=1,3
            do i=1,3
            vv(IDCP,i,j)=vv(ICV,i,j)
            vv(IDC,i,j)=vv(ICVP,i,j)
            enddo
            enddo
          else
            do j=1,3
            do i=1,3
            vv(IDCP,i,j)=vv(ICVP,i,j)   !vv(ICV,i,j)
            vv(IDC,i,j)=vv(ICV,i,j)     !vv(ICVP,i,j)
            enddo
            enddo
          endif
        enddo
      elseif(kd==kdovst) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          do j=1,3
          do i=1,3
          vv(IDC,i,j)=0.d0
          enddo
          enddo
        enddo
      elseif(kd.eq.kdsld.and.idis(nb)==1) then
!       if(ical_vect) then
!       else
        do 1220 ISLD=1,2
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
!
        IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
        IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
        th(ISLD)=rot_ang(IMATS(ISLD))
        unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
        unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
        unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
        CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD)=bb(j,i,ISLD)
        enddo
        enddo
!
        IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
        IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
        th(ISLD2)=rot_ang(IMATS(ISLD2))
        unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
        unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
        unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD2)=bb(j,i,ISLD2)
        enddo
        enddo
!
        DO 1160 IBFL=IBFS,IBFE
        wi1=wifsld(IBFL)
        wi2=1.d0-wi1
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
!
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        ICFO=LCYCOLD(IBFL)
        ICVBO=LVEDGE(1,ICFO)
        vell(1,1)=wi1*vv(ICVP,1,1)+wi2*vv(ICVBO,1,1)
        vell(1,2)=wi1*vv(ICVP,1,2)+wi2*vv(ICVBO,1,2)
        vell(1,3)=wi1*vv(ICVP,1,3)+wi2*vv(ICVBO,1,3)
        vell(2,1)=wi1*vv(ICVP,2,1)+wi2*vv(ICVBO,2,1)
        vell(2,2)=wi1*vv(ICVP,2,2)+wi2*vv(ICVBO,2,2)
        vell(2,3)=wi1*vv(ICVP,2,3)+wi2*vv(ICVBO,2,3)
        vell(3,1)=wi1*vv(ICVP,3,1)+wi2*vv(ICVBO,3,1)
        vell(3,2)=wi1*vv(ICVP,3,2)+wi2*vv(ICVBO,3,2)
        vell(3,3)=wi1*vv(ICVP,3,3)+wi2*vv(ICVBO,3,3)
!
! --- 2=>G  (x2=> Tg=B2*T2*B2^-1)
        do j=1,3
        do i=1,3
        x2(i,j)=bb(1,i,ISLD2)*vell(1,j)
     &         +bb(2,i,ISLD2)*vell(2,j)
     &         +bb(3,i,ISLD2)*vell(3,j)
        enddo
        enddo
        do j=1,3
        do i=1,3
        g2(i,j)=x2(i,1)*bb(j,1,ISLD2)
     &         +x2(i,2)*bb(j,2,ISLD2)
     &         +x2(i,3)*bb(j,3,ISLD2)
        enddo
        enddo
! --- G=>2  (x1=>  T2=B2^-1*Tg*B2)
        do j=1,3
        do i=1,3
        x2(i,j)=bb(i,1,ISLD)*g2(1,j)
     &         +bb(i,2,ISLD)*g2(2,j)
     &         +bb(i,3,ISLD)*g2(3,j)
        enddo
        enddo
        do j=1,3
        do i=1,3
        vv(IDC,i,j)=x2(i,1)*bb(1,j,ISLD)
     &             +x2(i,2)*bb(2,j,ISLD)
     &             +x2(i,3)*bb(3,j,ISLD)
!
        enddo
        enddo
 1160   enddo
 1220   enddo
!       endif
      elseif(kd.eq.kdsld.and.idis(nb)==2) then
        do 1300 ISLD=1,2
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
!
        IIMATS(ISLD)=MAT_BCIDX(nb,ISLD)
        IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
        th(ISLD)=rot_ang(IMATS(ISLD))
        unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
        unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
        unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
        CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
!        do i=1,3
!        do j=1,3
!        rbb(i,j,ISLD)=bb(j,i,ISLD)
!        enddo
!        enddo
!new
        IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
        IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
!        th(ISLD2)=rot_ang(IMATS(ISLD2))
        unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
        unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
        unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
!        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
!        do i=1,3
!        do j=1,3
!        rbb(i,j,ISLD2)=bb(j,i,ISLD2)
!        enddo
!        enddo
!
        DO IBFL=IBFS,IBFE
        wi1=wifsld(IBFL)
        wi2=1.d0-wi1
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICV=LVEDGE(1,ICFL)
!
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        ICFO=LCYCOLD(IBFL)
        ICVBO=LVEDGE(1,ICFO)
!
        th(ISLD2)=OPPANG(IBFL)
        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
!
        vell(1,1)=wi1*vv(ICVP,1,1)+wi2*vv(ICVBO,1,1)
        vell(1,2)=wi1*vv(ICVP,1,2)+wi2*vv(ICVBO,1,2)
        vell(1,3)=wi1*vv(ICVP,1,3)+wi2*vv(ICVBO,1,3)
        vell(2,1)=wi1*vv(ICVP,2,1)+wi2*vv(ICVBO,2,1)
        vell(2,2)=wi1*vv(ICVP,2,2)+wi2*vv(ICVBO,2,2)
        vell(2,3)=wi1*vv(ICVP,2,3)+wi2*vv(ICVBO,2,3)
        vell(3,1)=wi1*vv(ICVP,3,1)+wi2*vv(ICVBO,3,1)
        vell(3,2)=wi1*vv(ICVP,3,2)+wi2*vv(ICVBO,3,2)
        vell(3,3)=wi1*vv(ICVP,3,3)+wi2*vv(ICVBO,3,3)
!
! --- 2=>G  (x2=> Tg=B2*T2*B2^-1)
        do j=1,3
        do i=1,3
        x2(i,j)=bb(1,i,ISLD2)*vell(1,j)
     &         +bb(2,i,ISLD2)*vell(2,j)
     &         +bb(3,i,ISLD2)*vell(3,j)
        enddo
        enddo
        do j=1,3
        do i=1,3
        g2(i,j)=x2(i,1)*bb(j,1,ISLD2)
     &         +x2(i,2)*bb(j,2,ISLD2)
     &         +x2(i,3)*bb(j,3,ISLD2)
        enddo
        enddo
! --- G=>2  (x1=>  T2=B2^-1*Tg*B2)
        do j=1,3
        do i=1,3
        x2(i,j)=bb(i,1,ISLD)*g2(1,j)
     &         +bb(i,2,ISLD)*g2(2,j)
     &         +bb(i,3,ISLD)*g2(3,j)
        enddo
        enddo
        do j=1,3
        do i=1,3
        vv(IDC,i,j)=x2(i,1)*bb(1,j,ISLD)
     &             +x2(i,2)*bb(2,j,ISLD)
     &             +x2(i,3)*bb(3,j,ISLD)
        enddo
        enddo
        enddo
 1300   enddo

      elseif(kd.eq.kdintr) then
! --- interface
        do 1450 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 402 j=1,3
        do 402 i=1,3
        vv(IDC,i,j)=vv(ICV,i,j)
        vv(IDCP,i,j)=vv(ICVP,i,j)
 402    continue
 1450   continue
!      elseif(kd.eq.kdtchi) then
!--< 1.4 touch-inlet boundary >--
!        do 1450 IBFL=IBFS,IBFE
!        ICFL=LBC_SSF(IBFL)
!        IDC=LVEDGE(2,ICFL)
!        ICFP=LCYCSF(IBFL)
!        ICVP=LVEDGE(1,ICFP)
!        IDCP=LVEDGE(1,ICFP)
!        do 401 j=1,3
!        do 401 i=1,3
!        x1(i,j)=aa(i,1)*vv(IDCP,1,j)
!     &         +aa(i,2)*vv(IDCP,2,j)
!     &         +aa(i,3)*vv(IDCP,3,j)
! 401    continue
!        do 402 j=1,3
!        do 402 i=1,3
!        vv(IDC,i,j)= x1(i,1)*aa(j,1)
!     &              +x1(i,2)*aa(j,2)
!     &              +x1(i,3)*aa(j,3)
! 402    continue
! 1450   continue
      elseif(kd.eq.kdsymm) then
!--< 1.2 symmetric boundary >--x
        do 1250 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
          do 201 j=1,3
          do 201 i=1,3
          x1(i,j)=SFAREA(i,ICFL)*SFAREA(j,ICFL)
 201      continue
          do 202 j=1,3
          y(j)=vv(ICV,j,1)*SFAREA(1,ICFL)
     &        +vv(ICV,j,2)*SFAREA(2,ICFL)
     &        +vv(ICV,j,3)*SFAREA(3,ICFL)
          do 202 i=1,3
          x2(i,j)=(x1(i,1)*vv(ICV,1,j)
     &            +x1(i,2)*vv(ICV,2,j)
     &            +x1(i,3)*vv(ICV,3,j))
     &           +(x1(1,j)*vv(ICV,i,1)
     &            +x1(2,j)*vv(ICV,i,2)
     &            +x1(3,j)*vv(ICV,i,3))
  202     continue
          prd=4.d0*(SFAREA(1,ICFL)*y(1)
     &             +SFAREA(2,ICFL)*y(2)
     &             +SFAREA(3,ICFL)*y(3))
          do 203 j=1,3
          do 203 i=1,3
          vv(IDC,i,j)=vv(ICV,i,j)-(x2(i,j)+x2(i,j))+prd*x1(i,j)
  203     continue
 1250     continue
!
!--< 1.3 other boundary >--
      else
        if(ical_vect) then 
          do j=1,3
          do i=1,3
!          do IBFL=IBFS,IBFE
!          ICFL=LBC_SSF(IBFL)
!          ICV=LVEDGE(1,ICFL)
!          IDC=LVEDGE(2,ICFL)
!          vv(IDC,i,j)=vv(ICV,i,j)
!          enddo
!NEC$ vector
          do IBFL=IBFS,IBFE
          tmpfac(IBFL,1)=vv(LVEDGE(1,LBC_SSF(IBFL)),i,j)
          enddo
!NEC$ vector
          do IBFL=IBFS,IBFE
          vv(LVEDGE(2,LBC_SSF(IBFL)),i,j)=tmpfac(IBFL,1)
          enddo
          enddo
          enddo
        else
          do IBFL=IBFS,IBFE
          do j=1,3
          do i=1,3
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          vv(IDC,i,j)=vv(ICV,i,j)
          enddo
          enddo
          enddo
        endif

!
      endif
 1000 continue
!
      return
      end subroutine dc_vgrad
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine rotth(unit,th,bb)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      real*8 ,intent(inout)  :: unit(3),th
      real*8 ,intent(out)    :: bb(3,3)
      real*8  :: aa,ax,ay,az,cost,sint,cosu
      real*8,parameter :: SML=1.d-8
!
      ax=unit(1)
      ay=unit(2)
      az=unit(3)
      aa=sqrt(ax*ax+ay*ay+az*az)
      unit(1)=ax/aa
      unit(2)=ay/aa
      unit(3)=az/aa
      if(min(aa,abs(th)).gt.SML) goto 200
      bb(:,:)=0.d0
      bb(1,1)=1.d0
      bb(2,2)=1.d0
      bb(3,3)=1.d0
      return
!
 200  continue
      ax=ax/aa
      ay=ay/aa
      az=az/aa
      cost=cos(th)
      sint=sin(th)
      cosu=1.d0-cost
      bb(1,1)=ax*ax*cosu+   cost  !1
      bb(2,1)=ax*ay*cosu-az*sint  !2
      bb(3,1)=ax*az*cosu+ay*sint  !3
      bb(1,2)=ay*ax*cosu+az*sint  !4
      bb(2,2)=ay*ay*cosu+   cost  !5
      bb(3,2)=ay*az*cosu-ax*sint  !6
      bb(1,3)=az*ax*cosu-ay*sint  !7
      bb(2,3)=az*ay*cosu+ax*sint  !8
      bb(3,3)=az*az*cosu+   cost  !9
      return
      end subroutine rotth
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE AXB_UNIT_C(A,B,C)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!---------------------------------------
!     A X B = C 
!---------------------------------------
      real*8 ,intent(in)   :: A(3),B(3)
      real*8 ,intent(OUT)  :: C(3)
      real*8 :: D
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
!
!
      C(3)=(A(1)*B(2)-A(2)*B(1))
      C(2)=(A(3)*B(1)-A(1)*B(3))
      C(1)=(A(2)*B(3)-A(3)*B(2))
      D=DSQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))+SML
      C(3)=C(3)/D
      C(2)=C(2)/D
      C(1)=C(1)/D
!
      RETURN
      END subroutine AXB_UNIT_C
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_symprs_vof(nko,LVEDGE,LBC_SSF,LCYCSF,mat_cal,aks)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_boundary,only : kdbcnd,kdprdc,LBC_INDEX,nbcnd,kdtchi,
     &                           MAT_BCIDX,kdintr,kxnone,kdfire,kdsymm,
     &                           kdintr,kdbuff,kdsld,kdcvd,idis
     &                           ,kdovst
      use module_scalar,  only : ivof
!
! 1. Set scalar at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(inout) :: aks       (  MXALLCVR,mxrans)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd.eq.kdprdc.or.kd==kdsld) then
! --- Periodic & interface BC:
        if(idis(nb)==0) then
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          aks(IDC,ivof)=aks(ICVP,ivof)
          aks(IDCP,ivof)=aks(ICV,ivof)
 1250     continue
        elseif(idis(nb)>=1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          aks(IDC,ivof)=aks(ICVP,ivof)
          enddo
        endif
      elseif(kd.eq.kxnone.or.
     &       kd.eq.kdfire.or.
     &       kd.eq.kdsymm.or.
     &       kd.eq.kdintr.or.
!     &       kd.eq.kdbuff.or.
     &       kd.eq.kdcvd) then
! --- Other BC:
        do 1150 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        aks(IDC,ivof)=aks(ICV,ivof)
 1150   continue
      endif
 1000 continue
!
      return
!
      end subroutine dc_symprs_vof 
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_symMHD(nko,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &                     mat_cal,aa)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_boundary,only : kdbcnd,kdprdc,LBC_INDEX,nbcnd,kdtchi,
     &                           MAT_BCIDX,kdintr,kdsld,kdcvd,idis
     &                           ,kdovst
      use module_model,only    : ical_vect
      use module_metrix,only   : tmpfac=>d2vect
!
! 1. Set scalar at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(inout) :: aa        (  MXALLCV,nko)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,I
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,ICVBO,ICFO
      real*8  :: wi1,wi2
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc) then
! --- Periodic & interface BC:
        if(idis(nb)==0) then
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do 201 k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          aa(IDCP,k)=aa(ICV,k)
  201     continue
 1250     continue
        elseif(idis(nb)==1) then
          do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)

          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          do k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          enddo
          enddo
        endif
      elseif(kd==kdsld) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        IDC=LVEDGE(2,ICFL)
        ICVA=LVEDGE(1,ICFL)

        ICFP=LCYCSF(IBFL)
        ICFO=LCYCOLD(IBFL)
!!!        ICVP=LVEDGE(1,ICFP)


        ICVB=LVEDGE(1,ICFP)
        ICVBO=LVEDGE(1,ICFO)
        wi1=wifsld(IBFL)
        wi2=1.d0-wi1
        do k=1,nko
        aa(IDC,k)=wi1*aa(ICVB,k)+wi2*aa(ICVBO,k)  !?zhang????
!!!!        aa(IDC,k)=aa(ICVP,k)
        enddo
        enddo
      elseif(kd==kdintr) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          aa(IDCP,k)=aa(ICV,k)
          enddo
        enddo
      else
! --- Other BC:
        if(ical_vect) then 
          do k=1,nko  !huilai
!NEC$ vector
          do IBFL=IBFS,IBFE
          tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),k)
          ENDDO
!NEC$ vector
          do IBFL=IBFS,IBFE
          aa(LVEDGE(2,LBC_SSF(IBFL)),k)=tmpfac(IBFL,1)
          ENDDO
          enddo
        else
          do 1150 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 101 k=1,nko
          aa(IDC,k)=aa(ICV,k)
  101     continue
 1150     continue
        endif
!
      endif
 1000 continue
!
      return
!
      end subroutine dc_symMHD 
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_sym_curt_MHD
     &           (nko,LVEDGE,LBC_SSF,LCYCSF,wifsld,LCYCOLD,
     &            mat_cal,aa)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_boundary,only : kdbcnd,kdprdc,LBC_INDEX,nbcnd,kdtchi,
     &                           MAT_BCIDX,kdintr,kdsld,kdcvd,idis
     &                           ,kdovst
!
! 1. Set scalar at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: wifsld(  MXSSFBC_SLD)
      integer,intent(in)    :: LCYCOLD( MXSSFBC_SLD)
      real*8 ,intent(inout) :: aa        (  MXALLCV,nko)
!
! --- [local entities]
!
      integer :: IMAT,IIMAT,I
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ICVA,ICVB,ICVBO,ICFO
      real*8  :: wi1,wi2
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) goto 1000
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc) then
! --- Periodic & interface BC:
        if(idis(nb)==0) then
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do 201 k=1,nko
          aa(IDC,k)=aa(ICVP,k)
          aa(IDCP,k)=aa(ICV,k)
  201     continue
 1250     continue
        elseif(idis(nb)==1) then
          call FFRABORT(1,'dc_sym_curt_MHD-1')
        endif
      elseif(kd==kdsld) then
         call FFRABORT(1,'dc_sym_curt_MHD-2')
      elseif(kd==kdintr) then
        do IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do k=1,nko
          aa(IDCP,k)=aa(ICVP,k)
          aa(IDC,k)=aa(ICV,k)
          enddo
        enddo
      else
! --- Other BC:
          do 1150 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 101 k=1,nko
          aa(IDC,k)=aa(ICV,k)
  101     continue
 1150     continue

      endif
 1000 continue
!
      return
!
      end subroutine dc_sym_curt_MHD 

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_vgrad_MHD
     &  (MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,SFAREA,
     &   LCYCOLD,wifsld,OPPANG,vv)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments] 
!
      use module_dimension 
      use module_boundary,only : set_rotmtrx,nbcnd,LBC_INDEX,
     &                           kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdintr,kdsld,MAT_BCIDX,rotsld,
     &                           idis,LBC_pair,
     &                           kmdirc,kmneum,kdprdc,
     &                           kdnature,kdinterf,kddefine
     &                           ,kdovst
      use module_material,only : rotati,ishaft,end,begin,ical_sld,
     &                           rot_ang
      use module_nowtime,only  : iter,time
      use module_model,only    : ical_vect
      
!
! 1. Set gradient of vector at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      INTEGER,INTENT(IN)    :: MAT_NO    (  0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(inout) :: vv        (  MXALLCV,3,3)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)

!
! --- [local entities] 
!
      real*8  :: x1(3,3),x2(3,3),y(3),xx(3,3),aa(3,3)
      real*8  :: g1(3,3),g2(3,3),wi1,wi2
      real*8  :: prd,unit(3,2),th(2),bb(3,3,2),rbb(3,3,2),vell(3,3)
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2),ICVBO,ICFO,ierr1=0
      integer :: nb,kd,kdm,kdF,kdA,ISLD,ISLD2,i,j
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
      integer :: ibfs1,ibfe1,ibfs2,ibfe2
!
!-< 1. Set values at dummy cell >-
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      IMAT=MAT_NO(IIMAT)
      if(.not.mat_cal(IIMAT)) cycle
      kdm=kdbcnd(7,nb)
      kdF=kdbcnd(6,nb)
      kdA=kdbcnd(5,nb)
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc) then 
          do IBFL=IBFS,IBFE 
          ICFL=LBC_SSF(IBFL) 
          ICV=LVEDGE(1,ICFL) 
          IDC=LVEDGE(2,ICFL) 
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          IDCP=LVEDGE(2,ICFP)
          do j=1,3
          do i=1,3
          x1(i,j)=aa(i,1)*vv(ICV,1,j)
     &           +aa(i,2)*vv(ICV,2,j)
     &           +aa(i,3)*vv(ICV,3,j)
          x2(i,j)=aa(1,i)*vv(ICVP,1,j)
     &           +aa(2,i)*vv(ICVP,2,j)
     &           +aa(3,i)*vv(ICVP,3,j)
          enddo
          enddo
!
          do j=1,3
          do i=1,3
          vv(IDCP,i,j)=x1(i,1)*aa(1,j)
     &                +x1(i,2)*aa(2,j)
     &                +x1(i,3)*aa(3,j)
          vv(IDC,i,j)= x2(i,1)*aa(j,1)
     &                +x2(i,2)*aa(j,2)
     &                +x2(i,3)*aa(j,3)
          enddo
          enddo
!
          enddo
      elseif(kd==kdintr) then
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do j=1,3
        do i=1,3
        vv(IDC,i,j)=vv(ICVP,i,j)
        vv(IDCP,i,j)=vv(ICV,i,j)
        enddo
        enddo
        enddo
!
      elseif(kd.eq.kdsymm) then
        do 1250 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        do j=1,3
        do i=1,3
        x1(i,j)=SFAREA(i,ICFL)*SFAREA(j,ICFL)
        enddo
        enddo
        do j=1,3
        y(j)=vv(ICV,j,1)*SFAREA(1,ICFL)
     &        +vv(ICV,j,2)*SFAREA(2,ICFL)
     &        +vv(ICV,j,3)*SFAREA(3,ICFL)
        do i=1,3
        x2(i,j)=(x1(i,1)*vv(ICV,1,j)
     &            +x1(i,2)*vv(ICV,2,j)
     &            +x1(i,3)*vv(ICV,3,j))
     &           +(x1(1,j)*vv(ICV,i,1)
     &            +x1(2,j)*vv(ICV,i,2)
     &            +x1(3,j)*vv(ICV,i,3))
        enddo
        enddo
        prd=4.d0*(SFAREA(1,ICFL)*y(1)
     &             +SFAREA(2,ICFL)*y(2)
     &             +SFAREA(3,ICFL)*y(3))
        do j=1,3
        do i=1,3
        vv(IDC,i,j)=vv(ICV,i,j)-(x2(i,j)+x2(i,j))+prd*x1(i,j)
        enddo
        enddo
1250    enddo
      else
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        do J=1,3
        do I=1,3
        vv(IDC,i,j)=vv(ICV,i,j)
        enddo
        enddo
        enddo
      endif

!
 1000 enddo
      return
      end subroutine dc_vgrad_MHD


!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine dc_symprvMHD
     &(nko,MAT_NO,LVEDGE,LBC_SSF,LCYCSF,mat_cal,
     & LCYCOLD,wifsld,OPPANG,
     & SFAREA,SFCENT,aa,imode)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!     imode==0 : velocity
!     imode==1 : grdc
!     imode==2 : rvx
!
! --- [module arguments]
!
      USE module_dimension
      use module_hpcutil,only  : my_rank
      use module_boundary,only : set_rotmtrx,nbcnd,LBC_INDEX,
     &                           kdbcnd,kdprdc,kdsymm,kdilet,kdolet,
     &                           kdtchi,kdintr,kdsld,MAT_BCIDX,rotsld,
     &                           idis,LBC_pair,kdbuff,kdshutr
     &                           ,kdovst
      use module_material,only : rotati,ishaft,end,begin,nsplit,
     &                           rot_ang
      use module_nowtime ,only : iter,time
      use module_material,only : ical_sld
      use module_model,only    : ical_vect
      use module_metrix,only   : tmpfac=>d2vect
      use module_metrix,only   : SHUTFL
!
! 1. Set vector at dummy cell
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: nko,imode
      integer,intent(in)    :: MAT_NO    (  0:MXMAT)
      integer,intent(in)    :: LVEDGE    (2,MXCVFAC)
      integer,intent(in)    :: LBC_SSF   (  MXSSFBC)
      integer,intent(in)    :: LCYCSF    (  MXSSFBC)
      logical,INTENT(IN)    :: mat_cal   (  0:MXMAT)
      real*8 ,intent(in)    :: SFAREA    (4,MXCVFAC)
      real*8 ,intent(in)    :: SFCENT    (3,MXCVFAC)
      real*8 ,intent(inout) :: aa        (  MXALLCV,3,nko)
      integer,intent(in)    :: LCYCOLD    (MXSSFBC_SLD)
      real*8 ,intent(in)    :: wifsld     (MXSSFBC_SLD)
      real*8 ,intent(in)    :: OPPANG     (MXSSFBC_SLD)
!
! --- [local entities]
!
      real*8  :: prd
      real*8  :: rot(3,3)
      real*8  :: shft(3),ufix(3,2)
      real*8  :: unit(3,2),th(2),bb(3,3,2)
      real*8  :: rbb(3,3,2)
      real*8  :: urot(3,2)
      real*8  :: vr(3,2),r(3,2)
      real*8  :: v0(3,2),dr(2),vell(3),wi1,wi2,dum1,dum2,dum3
      integer :: IMAT,IIMAT,IIMATS(2),IMATS(2)
      integer :: i,j,ICVBO,ICFO
      integer :: nb,kd,kdv,kdt,kdy,kdk,kdp,k,ISLD,ISLD2
      integer :: IBFS,IBFE,IBFL,IBFP,ICFL,ICFP,ICV,IDC,ICVP,IDCP
!
! --- 
!
      do 1000 nb=1,nbcnd
      IIMAT=MAT_BCIDX(nb,1)
      if(.not.mat_cal(IIMAT)) cycle
      kd=kdbcnd(0,nb)
      IBFS=LBC_INDEX(nb-1)+1
      IBFE=LBC_INDEX(nb)
      if(kd==kdprdc)then  !zhang???
        rot=0.d0
        call set_rotmtrx(nbcnd,kd,nb,rot)
!--< 1.1 periodic boundary >--
        do 1150 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 101 k=1,nko
!
        aa(IDC,1,k)=rot(1,1)*aa(ICVP,1,k)
     &             +rot(1,2)*aa(ICVP,2,k)
     &             +rot(1,3)*aa(ICVP,3,k)
        aa(IDC,2,k)=rot(2,1)*aa(ICVP,1,k)
     &             +rot(2,2)*aa(ICVP,2,k)
     &             +rot(2,3)*aa(ICVP,3,k)
        aa(IDC,3,k)=rot(3,1)*aa(ICVP,1,k)
     &             +rot(3,2)*aa(ICVP,2,k)
     &             +rot(3,3)*aa(ICVP,3,k)
!
          aa(IDCP,1,k)=rot(1,1)*aa(ICV,1,k)
     &                +rot(2,1)*aa(ICV,2,k)
     &                +rot(3,1)*aa(ICV,3,k)
          aa(IDCP,2,k)=rot(1,2)*aa(ICV,1,k)
     &                +rot(2,2)*aa(ICV,2,k)
     &                +rot(3,2)*aa(ICV,3,k)
          aa(IDCP,3,k)=rot(1,3)*aa(ICV,1,k)
     &                +rot(2,3)*aa(ICV,2,k)
     &                +rot(3,3)*aa(ICV,3,k)
 101    continue
 1150   continue
      elseif(kd==kdbuff) then 
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
!        
        do k=1,nko 
        aa(IDC,1,k)=aa(ICVP,1,k)
        aa(IDC,2,k)=aa(ICVP,2,k)
        aa(IDC,3,k)=aa(ICVP,3,k)
        aa(IDCP,1,k)=aa(ICV,1,k)
        aa(IDCP,2,k)=aa(ICV,2,k)
        aa(IDCP,3,k)=aa(ICV,3,k)
        enddo
        enddo
      elseif(kd==kdshutr) then
        call FFRABORT(1,'MSG: dc_symprvMHD=>kd==kdshutr')
        do IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        if(SHUTFL(IBFL)==0) then
          do k=1,nko
          aa(IDC,1,k)=aa(ICVP,1,k)
          aa(IDC,2,k)=aa(ICVP,2,k)
          aa(IDC,3,k)=aa(ICVP,3,k)
          aa(IDCP,1,k)=aa(ICV,1,k)
          aa(IDCP,2,k)=aa(ICV,2,k)
          aa(IDCP,3,k)=aa(ICV,3,k)
          enddo
        else
          do k=1,nko
          aa(IDC,1,k)=aa(ICV,1,k)   !aa(ICVP,1,k)
          aa(IDC,2,k)=aa(ICV,2,k)   !aa(ICVP,2,k)
          aa(IDC,3,k)=aa(ICV,3,k)   !aa(ICVP,3,k)
          aa(IDCP,1,k)=aa(ICVP,1,k)   !aa(ICV,1,k)
          aa(IDCP,2,k)=aa(ICVP,2,k)   !aa(ICV,2,k)
          aa(IDCP,3,k)=aa(ICVP,3,k)   !aa(ICV,3,k)
          enddo
        endif
        enddo
      elseif(kd==kdsld.and.idis(nb)==1) then  !ical_sld==1=>rotsld(nb)==1
        call FFRABORT(1,'MSG: dc_symprvMHD=kd==kdsld.and.idis(nb)==1')
        call set_rotmtrx(nbcnd,kd,nb,rot)
        do 1200 ISLD=1,2
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
        IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
        th(ISLD)=rot_ang(IMATS(ISLD))
        unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
        unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
        unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
        CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD)=bb(j,i,ISLD)
        enddo
        enddo
!
        IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
        IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
        th(ISLD2)=rot_ang(IMATS(ISLD2))
        unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
        unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
        unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
        CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD2)=bb(j,i,ISLD2)
        enddo
        enddo
!
        if(imode==0) then
          DO IBFL=IBFS,IBFE
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
!
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)
          ICVBO=LVEDGE(1,ICFO)
          IDC=LVEDGE(2,ICFL)
!
          r(1,ISLD)=SFCENT(1,ICFL)-begin(1,IMATS(ISLD))
          r(2,ISLD)=SFCENT(2,ICFL)-begin(2,IMATS(ISLD))
          r(3,ISLD)=SFCENT(3,ICFL)-begin(3,IMATS(ISLD))
!
          call AXB_UNIT_C(unit(:,ISLD),r(:,ISLD),vr(:,ISLD))
!
          dr(ISLD)=r(1,ISLD)*unit(1,ISLD)
     &            +r(2,ISLD)*unit(2,ISLD)
     &            +r(3,ISLD)*unit(3,ISLD)
          r(1,ISLD)=r(1,ISLD)-dr(ISLD)*unit(1,ISLD)
          r(2,ISLD)=r(2,ISLD)-dr(ISLD)*unit(2,ISLD)
          r(3,ISLD)=r(3,ISLD)-dr(ISLD)*unit(3,ISLD)
          dr(ISLD)=dsqrt(r(1,ISLD)*r(1,ISLD)
     &                  +r(2,ISLD)*r(2,ISLD)
     &                  +r(3,ISLD)*r(3,ISLD))
          v0(1,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(1,ISLD)
          v0(2,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(2,ISLD)
          v0(3,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(3,ISLD)
!
          
          dum1=SFCENT(1,ICFL)-begin(1,IMATS(ISLD2)) !7777
          dum2=SFCENT(2,ICFL)-begin(2,IMATS(ISLD2)) 
          dum3=SFCENT(3,ICFL)-begin(3,IMATS(ISLD2)) 
!
          r(1,ISLD2)=rbb(1,1,ISLD)*dum1
     &              +rbb(1,2,ISLD)*dum2
     &              +rbb(1,3,ISLD)*dum3
          r(2,ISLD2)=rbb(2,1,ISLD)*dum1
     &              +rbb(2,2,ISLD)*dum2
     &              +rbb(2,3,ISLD)*dum3
          r(3,ISLD2)=rbb(3,1,ISLD)*dum1
     &              +rbb(3,2,ISLD)*dum2
     &              +rbb(3,3,ISLD)*dum3
!
          call AXB_UNIT_C(unit(:,ISLD2),r(:,ISLD2),vr(:,ISLD2)) 
!
          dr(ISLD2)=r(1,ISLD2)*unit(1,ISLD2)
     &             +r(2,ISLD2)*unit(2,ISLD2)
     &             +r(3,ISLD2)*unit(3,ISLD2)
!
          r(1,ISLD2)=r(1,ISLD2)-dr(ISLD2)*unit(1,ISLD2)
          r(2,ISLD2)=r(2,ISLD2)-dr(ISLD2)*unit(2,ISLD2)
          r(3,ISLD2)=r(3,ISLD2)-dr(ISLD2)*unit(3,ISLD2)
          dr(ISLD2)=dsqrt(r(1,ISLD2)*r(1,ISLD2)
     &                   +r(2,ISLD2)*r(2,ISLD2)
     &                   +r(3,ISLD2)*r(3,ISLD2))
          v0(1,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(1,ISLD2)
          v0(2,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(2,ISLD2)
          v0(3,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(3,ISLD2)
!
          do k=1,nko
!
! --- rot2=> fix2  zhang8:rbb
!
          vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
          vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
          vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
          ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(1,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(1,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(1,ISLD2)
          ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(2,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(2,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(2,ISLD2)
          ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(3,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(3,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(3,ISLD2)
! --- fix2 => rot1
          urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,1,ISLD)*(ufix(3,ISLD2))
     &                 - v0(1,ISLD)
          urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,2,ISLD)*(ufix(3,ISLD2))
     &                 - v0(2,ISLD)
          urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,3,ISLD)*(ufix(3,ISLD2))
     &                 - v0(3,ISLD)
!
          aa(IDC,1,k)=urot(1,ISLD)
          aa(IDC,2,k)=urot(2,ISLD)
          aa(IDC,3,k)=urot(3,ISLD)
!zhang????
          
          enddo
          enddo
        elseif(imode==1) then
! --- V0=0.d0
          DO IBFL=IBFS,IBFE
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          ICFO=LCYCOLD(IBFL)
          ICVBO=LVEDGE(1,ICFO)
          do k=1,nko
          vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
          vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
          vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
! --- rot2=> fix2  zhang8:rbb
          ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))
     &                 +rbb(1,2,ISLD2)*(vell(2))
     &                 +rbb(1,3,ISLD2)*(vell(3))
          ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))
     &                 +rbb(2,2,ISLD2)*(vell(2))
     &                 +rbb(2,3,ISLD2)*(vell(3))
          ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))
     &                 +rbb(3,2,ISLD2)*(vell(2))
     &                 +rbb(3,3,ISLD2)*(vell(3))
! --- fix2 => rot1
          urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,1,ISLD)*(ufix(3,ISLD2))
          urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,2,ISLD)*(ufix(3,ISLD2))
          urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,3,ISLD)*(ufix(3,ISLD2))
          aa(IDC,1,k)=urot(1,ISLD)
          aa(IDC,2,k)=urot(2,ISLD)
          aa(IDC,3,k)=urot(3,ISLD)
!zhang????
          enddo
          enddo
        endif
 1200   enddo
      elseif(kd==kdsld.and.idis(nb)==2) then
        call FFRABORT(1,'MSG: dc_symprvMHD=kd==kdsld.and.idis(nb)==2')
        do 1300 ISLD=1,2
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
        IMATS(ISLD)=MAT_NO(IIMATS(ISLD))
        th(ISLD)=rot_ang(IMATS(ISLD))
        unit(1,ISLD)=end(1,IMATS(ISLD))-begin(1,IMATS(ISLD))
        unit(2,ISLD)=end(2,IMATS(ISLD))-begin(2,IMATS(ISLD))
        unit(3,ISLD)=end(3,IMATS(ISLD))-begin(3,IMATS(ISLD))
        CALL rotth(unit(:,ISLD),th(ISLD),bb(:,:,ISLD))
        do i=1,3
        do j=1,3
        rbb(i,j,ISLD)=bb(j,i,ISLD)
        enddo
        enddo
!
        IIMATS(ISLD2)=MAT_BCIDX(nb,ISLD2)
        IMATS(ISLD2)=MAT_NO(IIMATS(ISLD2))
        unit(1,ISLD2)=end(1,IMATS(ISLD2))-begin(1,IMATS(ISLD2))
        unit(2,ISLD2)=end(2,IMATS(ISLD2))-begin(2,IMATS(ISLD2))
        unit(3,ISLD2)=end(3,IMATS(ISLD2))-begin(3,IMATS(ISLD2))
!
        if(imode==0) then
          DO IBFL=IBFS,IBFE
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
!
          ICFL=LBC_SSF(IBFL)
          ICFP=LCYCSF(IBFL)
          ICFO=LCYCOLD(IBFL)
          ICV=LVEDGE(1,ICFL)
          ICVP=LVEDGE(1,ICFP)
          ICVBO=LVEDGE(1,ICFO)
          IDC=LVEDGE(2,ICFL)
!
          r(1,ISLD)=SFCENT(1,ICFL)-begin(1,IMATS(ISLD))
          r(2,ISLD)=SFCENT(2,ICFL)-begin(2,IMATS(ISLD))
          r(3,ISLD)=SFCENT(3,ICFL)-begin(3,IMATS(ISLD))
!
          call AXB_UNIT_C(unit(:,ISLD),r(:,ISLD),vr(:,ISLD))
!
          dr(ISLD)=r(1,ISLD)*unit(1,ISLD)
     &            +r(2,ISLD)*unit(2,ISLD)
     &            +r(3,ISLD)*unit(3,ISLD)
          r(1,ISLD)=r(1,ISLD)-dr(ISLD)*unit(1,ISLD)
          r(2,ISLD)=r(2,ISLD)-dr(ISLD)*unit(2,ISLD)
          r(3,ISLD)=r(3,ISLD)-dr(ISLD)*unit(3,ISLD)
          dr(ISLD)=dsqrt(r(1,ISLD)*r(1,ISLD)
     &                  +r(2,ISLD)*r(2,ISLD)
     &                  +r(3,ISLD)*r(3,ISLD))
          v0(1,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(1,ISLD)
          v0(2,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(2,ISLD)
          v0(3,ISLD)=dr(ISLD)*rotati(IMATS(ISLD))*vr(3,ISLD)
!new
          th(ISLD2)=OPPANG(IBFL)
          CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=bb(j,i,ISLD2)
          enddo
          enddo
!
          dum1=SFCENT(1,ICFL)-begin(1,IMATS(ISLD2)) !7777
          dum2=SFCENT(2,ICFL)-begin(2,IMATS(ISLD2)) 
          dum3=SFCENT(3,ICFL)-begin(3,IMATS(ISLD2)) 
          r(1,ISLD2)=rbb(1,1,ISLD)*dum1
     &              +rbb(1,2,ISLD)*dum2
     &              +rbb(1,3,ISLD)*dum3
          r(2,ISLD2)=rbb(2,1,ISLD)*dum1
     &              +rbb(2,2,ISLD)*dum2
     &              +rbb(2,3,ISLD)*dum3
          r(3,ISLD2)=rbb(3,1,ISLD)*dum1
     &              +rbb(3,2,ISLD)*dum2
     &              +rbb(3,3,ISLD)*dum3
!
          call AXB_UNIT_C(unit(:,ISLD2),r(:,ISLD2),vr(:,ISLD2))
!
          dr(ISLD2)=r(1,ISLD2)*unit(1,ISLD2)
     &             +r(2,ISLD2)*unit(2,ISLD2)
     &             +r(3,ISLD2)*unit(3,ISLD2)
!
          r(1,ISLD2)=r(1,ISLD2)-dr(ISLD2)*unit(1,ISLD2)
          r(2,ISLD2)=r(2,ISLD2)-dr(ISLD2)*unit(2,ISLD2)
          r(3,ISLD2)=r(3,ISLD2)-dr(ISLD2)*unit(3,ISLD2)
          dr(ISLD2)=dsqrt(r(1,ISLD2)*r(1,ISLD2)
     &                   +r(2,ISLD2)*r(2,ISLD2)
     &                   +r(3,ISLD2)*r(3,ISLD2))
          v0(1,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(1,ISLD2)
          v0(2,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(2,ISLD2)
          v0(3,ISLD2)=dr(ISLD2)*rotati(IMATS(ISLD2))*vr(3,ISLD2)
!
          do k=1,nko
!
! --- rot2=> fix2  zhang8:rbb
!
          vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
          vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
          vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
          ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(1,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(1,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(1,ISLD2)
          ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(2,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(2,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(2,ISLD2)
          ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))!+v0(1,ISLD2))
     &                 +rbb(3,2,ISLD2)*(vell(2))!+v0(2,ISLD2))
     &                 +rbb(3,3,ISLD2)*(vell(3))!+v0(3,ISLD2))
     &                 +v0(3,ISLD2)
! --- fix2 => rot1
          urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,1,ISLD)*(ufix(3,ISLD2))
     &                 - v0(1,ISLD)
          urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,2,ISLD)*(ufix(3,ISLD2))
     &                 - v0(2,ISLD)
          urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,3,ISLD)*(ufix(3,ISLD2))
     &                 - v0(3,ISLD)
!
          aa(IDC,1,k)=urot(1,ISLD)
          aa(IDC,2,k)=urot(2,ISLD)
          aa(IDC,3,k)=urot(3,ISLD)
          enddo
          enddo
        elseif(imode==1) then
! --- V0=0.d0
          DO IBFL=IBFS,IBFE
          wi1=wifsld(IBFL)
          wi2=1.d0-wi1
          ICFL=LBC_SSF(IBFL)
          IDC=LVEDGE(2,ICFL)
          ICFP=LCYCSF(IBFL)
          ICVP=LVEDGE(1,ICFP)
          ICFO=LCYCOLD(IBFL)
          ICVBO=LVEDGE(1,ICFO)
!new
          th(ISLD2)=OPPANG(IBFL)
          CALL rotth(unit(:,ISLD2),th(ISLD2),bb(:,:,ISLD2))
          do i=1,3
          do j=1,3
          rbb(i,j,ISLD2)=bb(j,i,ISLD2)
          enddo
          enddo
!
          do k=1,nko
          vell(1)=wi1*aa(ICVP,1,k)+wi2*aa(ICVBO,1,k)
          vell(2)=wi1*aa(ICVP,2,k)+wi2*aa(ICVBO,2,k)
          vell(3)=wi1*aa(ICVP,3,k)+wi2*aa(ICVBO,3,k)
! --- rot2=> fix2  zhang8:rbb
          ufix(1,ISLD2)=rbb(1,1,ISLD2)*(vell(1))
     &                 +rbb(1,2,ISLD2)*(vell(2))
     &                 +rbb(1,3,ISLD2)*(vell(3))
          ufix(2,ISLD2)=rbb(2,1,ISLD2)*(vell(1))
     &                 +rbb(2,2,ISLD2)*(vell(2))
     &                 +rbb(2,3,ISLD2)*(vell(3))
          ufix(3,ISLD2)=rbb(3,1,ISLD2)*(vell(1))
     &                 +rbb(3,2,ISLD2)*(vell(2))
     &                 +rbb(3,3,ISLD2)*(vell(3))
! --- fix2 => rot1
          urot(1,ISLD)= rbb(1,1,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,1,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,1,ISLD)*(ufix(3,ISLD2))
          urot(2,ISLD)= rbb(1,2,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,2,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,2,ISLD)*(ufix(3,ISLD2))
          urot(3,ISLD)= rbb(1,3,ISLD)*(ufix(1,ISLD2))
     &                 +rbb(2,3,ISLD)*(ufix(2,ISLD2))
     &                 +rbb(3,3,ISLD)*(ufix(3,ISLD2))
          aa(IDC,1,k)=urot(1,ISLD)
          aa(IDC,2,k)=urot(2,ISLD)
          aa(IDC,3,k)=urot(3,ISLD)
!zhang????
          enddo
          enddo
        endif
 1300   enddo
        
      elseif(kd.eq.kdsymm) then
!--< 1.2 symmetric boundary >--
        call FFRABORT(1,'MSG: dc_symprvMHD=kd.eq.kdsymm)')
	if(ical_vect) then 
	  do k=1,nko
          do IBFL=IBFS,IBFE
          prd=2.d0*(SFAREA(1,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &             +SFAREA(2,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &             +SFAREA(3,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),3,k))
            tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &          -prd*SFAREA(1,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
            aa(LVEDGE(2,LBC_SSF(IBFL)),1,k)=tmpfac(IBFL,1)
          enddo
          do IBFL=IBFS,IBFE
            prd=2.d0*(SFAREA(1,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &               +SFAREA(2,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &               +SFAREA(3,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),3,k))
            tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &          -prd*SFAREA(2,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
            aa(LVEDGE(2,LBC_SSF(IBFL)),2,k)=tmpfac(IBFL,1)
          enddo
          do IBFL=IBFS,IBFE
            prd=2.d0*(SFAREA(1,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),1,k)
     &               +SFAREA(2,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),2,k)
     &               +SFAREA(3,LBC_SSF(IBFL))
     &          *aa(LVEDGE(1,LBC_SSF(IBFL)),3,k))
            tmpfac(IBFL,1)=aa(LVEDGE(1,LBC_SSF(IBFL)),3,k)
     &          -prd*SFAREA(3,LBC_SSF(IBFL))
          enddo
          do IBFL=IBFS,IBFE
            aa(LVEDGE(2,LBC_SSF(IBFL)),3,k)=tmpfac(IBFL,1)
          enddo
          enddo
	else
          do 1250 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 201 k=1,nko
          prd=2.d0*(SFAREA(1,ICFL)*aa(ICV,1,k)
     &             +SFAREA(2,ICFL)*aa(ICV,2,k)
     &             +SFAREA(3,ICFL)*aa(ICV,3,k))
          aa(IDC,1,k)=aa(ICV,1,k)-prd*SFAREA(1,ICFL)
          aa(IDC,2,k)=aa(ICV,2,k)-prd*SFAREA(2,ICFL)
          aa(IDC,3,k)=aa(ICV,3,k)-prd*SFAREA(3,ICFL)
  201     continue
 1250     continue
 	endif

! --- interface BC
      elseif(kd.eq.kdintr) then
        do 1550 IBFL=IBFS,IBFE
        ICFL=LBC_SSF(IBFL)
        ICV=LVEDGE(1,ICFL)
        IDC=LVEDGE(2,ICFL)
        ICFP=LCYCSF(IBFL)
        ICVP=LVEDGE(1,ICFP)
        IDCP=LVEDGE(2,ICFP)
        do 501 k=1,nko
        aa(IDC,1,k)=aa(ICVP,1,k)
        aa(IDC,2,k)=aa(ICVP,2,k)
        aa(IDC,3,k)=aa(ICVP,3,k)
!
        aa(IDCP,1,k)=aa(ICV,1,k)
        aa(IDCP,2,k)=aa(ICV,2,k)
        aa(IDCP,3,k)=aa(ICV,3,k)
 501    continue
 1550   continue

!--< 1.3 other boundary >--
      elseif(kd/=kdtchi) then
        if(ical_vect) then
          call FFRABORT(1,'MSG: dc_symprvMHD=kd/=kdtchi & ical_vect')
        else
          do 1350 IBFL=IBFS,IBFE
          ICFL=LBC_SSF(IBFL)
          ICV=LVEDGE(1,ICFL)
          IDC=LVEDGE(2,ICFL)
          do 301 k=1,nko
          aa(IDC,1,k)=aa(ICV,1,k)
          aa(IDC,2,k)=aa(ICV,2,k)
          aa(IDC,3,k)=aa(ICV,3,k)
  301     continue
 1350     continue
        endif
      endif
 1000 continue
!
! --- --< 1.4 touch-inlet BC >--
!
      return
      end subroutine dc_symprvMHD
!
!
!
