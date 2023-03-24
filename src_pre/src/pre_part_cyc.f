!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine CYCLE
     & (ncell,nvrtx,nface,nedge,nssfbc,NCV,NCVFAC,NALLCV,IBCCYC,IBCTCI,
     &  IBCSLD,IBCINT,
     &  medge,mssfbc,mface,mcell,
     &  lbface,lcface,lvface)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_boundary,only : kdbcnd,kdprdc,kdtchi,kdsld,kdintr,idis,
     &                           kdbuff,kdshutr,kdpors
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: ncell,nvrtx,nface,nedge,nssfbc
      integer,intent(in) :: medge,mssfbc,mface,mcell
      integer,intent(in) :: NCV,NCVFAC,NALLCV,IBCCYC,IBCTCI,IBCSLD,
     &                      IBCINT
!
      integer,intent(in) :: lbface(2,mface)
      integer,intent(in) :: lcface(2,mface)
      integer,intent(in) :: lvface(4,mface)
!
! --- [local entities]
!
      integer :: NEcur,icyc,in1,in2,in3,in4,NSUF,isuf,icel,icelC,nn
      integer :: icyc1,icyc2,icyc3,icyc4,ip1,ip2,ip3,iq1,iq2,iq3
      integer :: in1C,in2C,in3C,isufp,isufpair,isufo
      integer :: isufp1,isufp2
      integer :: ivL,IS,IS2,IS3,nb,kd
      integer :: inA(4),inB(4),IF1,IF2,IC1,IC2,IC3
!---------------
! --- 
!---------------
      isufCYCLtot=2*(IBCCYC+IBCINT)
      
      allocate (ISUF_NODE(max(isufCYCLtot,1),4))             !cyc
      allocate (SUFtoCELL(max(isufCYCLtot,1),2))             !cyc
      allocate (ISUF_NODE_P(2,max(IBCCYC+IBCINT,1),4))       !cyc
!
      isufPAIRtot=IBCTCI
      allocate (ISUF_NODE_pair(max(isufPAIRtot,1),4))        !touch-inlet
      allocate (SUFtoCELL_pair(max(isufPAIRtot,1),3))        !touch-inlet
      allocate (ISUF_NODE_O(3,max(IBCTCI,1),4))              !touch-inlet
!
      ISUF_NODE=0
      SUFtoCELL=0
      ISUF_NODE_P=0
      ISUF_NODE_pair=0
      SUFtoCELL_pair=0
      ISUF_NODE_O=0
!---------------------------------
! --- periodic AND touch-inlet ---
!---------------------------------
      isufpair=0
      isuf=0
      do 100 IS=1,nface
      nb=lbface(1,IS)
      if(nb.le.0) goto 100
      kd=kdbcnd(0,nb)
      IS2=lbface(2,IS)
      if((kd==kdprdc.and.idis(nb)==0).or.
     &   (kd==kdbuff.and.idis(nb)==0).or.
     &   (kd==kdpors.and.idis(nb)==0).or.
     &   (kd==kdintr.and.idis(nb)==0).or.
     &   (kd==kdshutr.and.idis(nb)==0)
     &) then
!
! --- main face
!
        IC1=lcface(1,IS)
        IC2=lcface(1,IS2)
!
        isuf=isuf+1
        SUFtoCELL(isuf,1)=IC1
        SUFtoCELL(isuf,2)=IC2
        do 150 ivL=1,4
        ISUF_NODE(isuf,ivL)=lvface(ivL,IS)
 150    continue
! --- vice face
        isuf=isuf+1
        SUFtoCELL(isuf,1)=IC2
        SUFtoCELL(isuf,2)=IC1
        do 155 ivL=1,4
        ISUF_NODE(isuf,ivL)=lvface(ivL,IS2)
 155    continue
!
        if(isuf.gt.isufCYCLtot) then
          write(*,*) '*** ERR: isuf.gt.isufCYCLtot',isuf,isufCYCLtot
          stop
        endif
!
      elseif(kd==kdtchi.and.idis(nb)==0) then
! --- 
        IS3=lbface(2,IS2)
        IC1=lcface(1,IS)
        IC2=lcface(1,IS2)
        if(IS3.ne.0) IC3=lcface(1,IS3)
!
        isufpair=isufpair+1
        SUFtoCELL_pair(isufpair,1)=IC1
        SUFtoCELL_pair(isufpair,2)=IC2
        if(IS3.ne.0) SUFtoCELL_pair(isufpair,3)=IC3
        do 260 ivL=1,4
        ISUF_NODE_pair(isufpair,ivL)=lvface(ivL,IS)
 260    continue

!
        if(isufpair.gt.isufPAIRtot) then
          write(*,*) '*** ERR: isuf.gt.isufCYCLtot',isuf,isufCYCLtot
          stop
        endif
      endif
!
 100  continue
!
      if(isufpair.ne.isufPAIRtot) then 
        write(*,*) '*** ERR: isufpair.ne.isufPAIRtot',
     &                  isufpair,'',isufPAIRtot
        stop
      else
        write(* ,'(  "    ###    isufPAIRtot =",2i12)') 
     &           isufPAIRtot,isufpair
      endif
      if(isuf.ne.isufCYCLtot) then
        write(*,*) '*** ERR: isuf.ne.isufCYCLtot',isuf,'',isufCYCLtot
        stop
      else
        write(* ,'(  "    ###    isufCYCLtot =",2i12)') 
     &                           isufCYCLtot,isuf
      endif
!
! --- 
!
      isufp=0
      isufp2=0
      isufo=0
      do 200 IS=1,nface
      nb=lbface(1,IS)
      if(nb.eq.0) goto 200
      kd=kdbcnd(0,ABS(nb))
      if((kd==kdprdc.and.idis(abs(nb))==0).or.
     &   (kd==kdbuff.and.idis(abs(nb))==0).or.
     &   (kd==kdpors.and.idis(abs(nb))==0).or.
     &   (kd==kdintr.and.idis(abs(nb))==0).or.
     &   (kd==kdshutr.and.idis(abs(nb))==0)
     &    ) then
        if(nb.gt.0) then
          isufp=isufp+1
          do 250 ivL=1,4
          ISUF_NODE_P(1,isufp,ivL)=lvface(ivL,IS)
 250      continue
          IS2=lbface(2,IS)
          do 255 ivL=1,4
          ISUF_NODE_P(2,isufp,ivL)=lvface(ivL,IS2)
 255      continue
        elseif(nb.lt.0) then
          isufp2=isufp2+1
        endif
      elseif(kd==kdtchi.and.idis(abs(nb))==0) then
! --- touch-inlet + interface BC
        IS2=lbface(2,IS)
        IS3=lbface(2,IS2)
        if(IS3.eq.0) goto 200
        isufo=isufo+1
        do 357 ivL=1,4
        ISUF_NODE_O(1,isufo,ivL)=lvface(ivL,IS)
 357    continue
        do 350 ivL=1,4
        ISUF_NODE_O(2,isufo,ivL)=lvface(ivL,IS2)
 350    continue
        do 355 ivL=1,4
        ISUF_NODE_O(3,isufo,ivL)=lvface(ivL,IS3)
 355    continue
      endif
 200  continue
!
      if(isufp.ne.isufp2) then
        write(*,*) '*** ERR: isufp NOT equal to isufp2',isufp,isufp2
        stop
      else if(isufp.ne.IBCCYC+IBCINT) then
        write(*,*) '*** ERR: isufp NOT equal to IBCCYC+IBCINT',isufp,
     &   IBCCYC+IBCINT
        stop
      else if(isufo.ne.IBCTCI) then
        write(*,*) '*** Warnig: isufo NOT equal to IBCTCI',isufo,
     &   IBCTCI
      endif
      NSUF_NODE_P=isufp
      NSUF_NODE_O=isufo
!
      return
!
      end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
