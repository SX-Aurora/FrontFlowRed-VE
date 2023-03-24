!
!***
!*** subroutine COMM
!*** subroutine re_comm
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine COMM(mcell,mface,nface,lvcell,lvface,lcface,lbface)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_io,only         : ifle,ifll
      use module_parameter,only  : NIFACE
      use module_boundary,only   : kdbcnd,nbcnd,idis,boundName,
     &                             kdprdc,kdintr,kdsld,kdtchi,kdbuff,
     &                             kdshutr,kdpors
      use module_partitioner,only: LBC_SSF=> WK6
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)     :: mcell,mface,nface
      integer,intent(inout)  :: lvcell(8,mcell)
      integer,intent(inout)  :: lvface(4,mface)
      integer,intent(inout)  :: lcface(2,mface)
      integer,intent(inout)  :: lbface(2,mface)
!
! --- [local entities]
!
      integer :: iv,ic,ierr1=0
      integer :: ivinL,ivexL,ivL,ivLG,ivexLG,ivexLGS,ivinLGS
      integer :: icinL,icL,icLG,ivLGw
      integer :: ivW,inWin,icpu0,icpuW
      integer :: icpuc,ncpuc
      integer :: icpu,icou,icel,i,j,k,ncpu,ncpu2,ncon,local,long,nn
      integer :: iext
      integer :: icpu1,icpu2,icpu3,icpu4,icpu5,icpu6,icpu7,icpu8
      integer :: icup1,icup2,icup3,icup4,icup5,icup6,icup7,icup8
      integer :: in1,in2,in3,in4,in5,in6,in7,in8
      integer :: im1,im2,im3,im4,im5,im6,im7,im8
      integer :: ip1,ip2,ip3,ip4
      integer :: isuf,ityp,isufpair
      integer :: icel1,icel2,icel3
      integer :: jj,js,je
      integer :: iv1,iv2,iw,iw3,ii,kmax1,kmax2
      integer :: OpenStatus
      integer :: BCSLD(nbcnd),LBC_IND(0:nbcnd),nb,nbx,ISB,NBSSF
      integer :: IS,IS1,IS2,IC1,IC2,IBFL,IBFS,IBFE,kd
      character*80 :: nameBC1,nameBC2
!
      write (ifll,'(a)') '  - LOCAL NODE'
!
!----------------------
! --- List sliding BC
!----------------------
!
      LBC_IND=0
      do nb=1,nbcnd
      kd=kdbcnd(0,nb)
!      if((kd==kdprdc.and.idis(nb)==1).or.
!     &   (kd==kdsld.and.idis(nb)==1).or.
!     &   (kd==kdsld.and.idis(nb)==0).or.
!     &   (kd==kdintr.and.idis(nb)==1).or.
!     &   (kd==kdtchi.and.idis(nb)==1)
!     &   ) then
        LBC_IND(nb)=LBC_IND(nb-1)
        ISB=0
        do IS=NIFACE,nface
        NBSSF=lbface(1,IS)
        if(nb==NBSSF) then
          ISB=ISB+1
        endif
        enddo
!
        do IS=NIFACE,nface
        NBSSF=lbface(1,IS)
        if(nb==-NBSSF) then
          ISB=ISB+1
        endif
        enddo
        LBC_IND(nb)=LBC_IND(nb)+ISB
!      endif
      enddo
!
      allocate(LBC_SSF(max(LBC_IND(nbcnd),1)),stat=ierr1) 
      if(ierr1.ne.0) 
     &     stop 'stop at allocating LBC_SSF(:) in COMM'
      LBC_SSF=0
!
      DO nb=1,nbcnd
      kd=kdbcnd(0,nb)
!      if((kd==kdprdc.and.idis(nb)==1).or.
!     &   (kd==kdsld.and.idis(nb)==1).or.
!     &   (kd==kdsld.and.idis(nb)==0).or.
!     &   (kd==kdintr.and.idis(nb)==1).or.
!     &   (kd==kdtchi.and.idis(nb)==1)
!     &   ) then
        ISB=LBC_IND(nb-1)
        do IS=NIFACE,nface
        NBSSF=lbface(1,IS)
        if(nb==NBSSF) then
          ISB=ISB+1
          LBC_SSF(ISB)=IS
        ENDIF
        enddo
!
        do IS=NIFACE,nface
        NBSSF=lbface(1,IS)
        if(nb==-NBSSF) then
          ISB=ISB+1
          LBC_SSF(ISB)=IS
        endif
        enddo
!      endif
      enddo
!
! --- +-----------------+
! --- | LOCAL vertex ID |
! --- +-----------------+
!
! --- NODE_ID(i,2): NODE_ID(i,2)='CPU NO.'=my_rank+1
! --- NODE_ID(i,1): NODE_ID(iv_gloabl,1)=iv_local
! --- NODtotL(icpu): Total inner Vertex no. in icpu 
!
      allocate (NODtotL(PETOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocating NODtotL(:) in COMM'
!
      NODtotL(:)=0
      do 50 icpu=1,PETOT
      ivinL=0
      do 60 iv=1,INODTOT
      if(NODE_ID(iv,2).eq.icpu) then
        ivinL=ivinL+1
        NODE_ID(iv,1)=ivinL
      endif
 60   continue
      NODtotL(icpu)=ivinL
 50   continue
!
      write (ifll,'(a)') '  - LOCAL CELL'
!
! --- +------------+
! --- | LOCAL CELL |
! --- +------------+
!
      allocate (ELMstack(0:PETOT), ELMtotL(PETOT), WK1(ICELTOT))
!
      ELMstack=0
      ELMtotL=0
!------------------------------------------------
! --- Decompositon Cell into CPU
!------------------------------------------------
      do 100 icpu=1,PETOT
      icL=0
      WK1=0
      do 110 ic=1,ICELTOT
      ityp=ICELTYP(ic)
! --- 
      if(ityp.eq.1) then
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
! --- 
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu)then
!
          icL=icL+1
          WK1(ic)=1
        endif
      elseif(ityp.eq.2) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
! --- 
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu.or.
     &     icpu5.eq.icpu)then
!
          icL=icL+1
          WK1(ic)=1
        endif

      elseif(ityp.eq.3) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
        icpu6=NODE_ID(in6,2)
! --- 
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu.or.
     &     icpu5.eq.icpu.or.
     &     icpu6.eq.icpu)then
!
          icL=icL+1
          WK1(ic)=1
        endif
      elseif(ityp.eq.4) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        in7=lvcell(7,ic)
        in8=lvcell(8,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
        icpu6=NODE_ID(in6,2)
        icpu7=NODE_ID(in7,2)
        icpu8=NODE_ID(in8,2)
! --- 
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu.or.
     &     icpu5.eq.icpu.or.
     &     icpu6.eq.icpu.or.
     &     icpu7.eq.icpu.or.
     &     icpu8.eq.icpu)then
!
          icL=icL+1
          WK1(ic)=1
        endif
      endif
 110  continue
!------------------------------------------------
! --- cyclc BC and interface BC
!------------------------------------------------
      if(isufCYCLtot.gt.0) then
        do 120 isuf=1,isufCYCLtot-1,2
!        do 120 isuf=1,isufCYCLtot
        in1=ISUF_NODE(isuf,1)
        in2=ISUF_NODE(isuf,2)
        in3=ISUF_NODE(isuf,3)
        in4=ISUF_NODE(isuf,4)
        in5=ISUF_NODE(isuf+1,1)
        in6=ISUF_NODE(isuf+1,2)
        in7=ISUF_NODE(isuf+1,3)
        in8=ISUF_NODE(isuf+1,4)
        icel1=SUFtoCELL(isuf,1)
        icel2=SUFtoCELL(isuf,2)
!
        if(in4.ne.0) then
! --- 
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu4=NODE_ID(in4,2)
          icpu5=NODE_ID(in5,2)
          icpu6=NODE_ID(in6,2)
          icpu7=NODE_ID(in7,2)
          icpu8=NODE_ID(in8,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu4.eq.icpu.or.
     &       icpu5.eq.icpu.or.
     &       icpu6.eq.icpu.or.
     &       icpu7.eq.icpu.or.
     &       icpu8.eq.icpu
     &    ) then
!
             if(WK1(icel1).eq.0) then
               icL=icL+1
               WK1(icel1)=1
             endif
             if(WK1(icel2).eq.0) then
               icL=icL+1
               WK1(icel2)=1
             endif
          endif
        else
! --- 
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu5=NODE_ID(in5,2)
          icpu6=NODE_ID(in6,2)
          icpu7=NODE_ID(in7,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu5.eq.icpu.or.
     &       icpu6.eq.icpu.or.
     &       icpu7.eq.icpu
     &) then
             if(WK1(icel1).eq.0) then
               icL=icL+1
               WK1(icel1)=1
             endif
             if(WK1(icel2).eq.0) then
               icL=icL+1
               WK1(icel2)=1
             endif
          endif
        endif
!
 120    continue
      endif
!------------------------
! --- touch-inlet BC
!------------------------
      if(isufPAIRtot.gt.0) then
        do 600 isufpair=1,isufPAIRtot
        in1=ISUF_NODE_pair(isufpair,1)
        in2=ISUF_NODE_pair(isufpair,2)
        in3=ISUF_NODE_pair(isufpair,3)
        in4=ISUF_NODE_pair(isufpair,4)
        icel1=SUFtoCELL_pair(isufpair,1)
        icel2=SUFtoCELL_pair(isufpair,2)
        icel3=SUFtoCELL_pair(isufpair,3)
        if(in4.ne.0) then
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu4=NODE_ID(in4,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu4.eq.icpu) then
             if(WK1(icel1).eq.0) then
               icL=icL+1
               WK1(icel1)=1
             endif
             if(WK1(icel2).eq.0) then
               icL=icL+1
               WK1(icel2)=1
             endif
             if(icel3.ne.0) then
               if(WK1(icel3).eq.0) then
                 icL=icL+1
                 WK1(icel3)=1
               endif
             endif
          endif
        else
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu) then
             if(WK1(icel1).eq.0) then
               icL=icL+1
               WK1(icel1)=1
             endif
             if(WK1(icel2).eq.0) then
               icL=icL+1
               WK1(icel2)=1
             endif
             if(icel3.ne.0) then
               if(WK1(icel3).eq.0) then
                 icL=icL+1
                 WK1(icel3)=1
               endif
             endif
          endif
        endif
 600    continue
      endif
!---------------------------
! --- SLIDING BC
!---------------------------
      BCSLD(:)=0
      do 660 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      IBFS=LBC_IND(nb-1)+1
      IBFE=LBC_IND(nb)
      if(kd==kdsld.and.idis(nb)==0) then
!      if(kd==kdsld) then
        do IBFL=IBFS,IBFE
        IS1=LBC_SSF(IBFL)
        IS2=lbface(2,IS1)
        IC1=lcface(1,IS1)
        IC2=lcface(1,IS2)
!
        in1=lvface(1,IS1)
        in2=lvface(2,IS1)
        in3=lvface(3,IS1)
        in4=lvface(4,IS1)
!
        im1=lvface(1,IS2)
        im2=lvface(2,IS2)
        im3=lvface(3,IS2)
        im4=lvface(4,IS2)
!
        if(in4.ne.0) then
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu4=NODE_ID(in4,2)
          icup1=NODE_ID(im1,2)
          icup2=NODE_ID(im2,2)
          icup3=NODE_ID(im3,2)
          icup4=NODE_ID(im4,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu4.eq.icpu.or.
     &       icup1.eq.icpu.or.
     &       icup2.eq.icpu.or.
     &       icup3.eq.icpu.or.
     &       icup4.eq.icpu)then
             BCSLD(nb)=1
             GOTO 660
          endif
        ELSE
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icup1=NODE_ID(im1,2)
          icup2=NODE_ID(im2,2)
          icup3=NODE_ID(im3,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icup1.eq.icpu.or.
     &       icup2.eq.icpu.or.
     &       icup3.eq.icpu)then
             BCSLD(nb)=1
             GOTO 660
          endif
        endif
        enddo
      endif
 660  continue
!
      DO nb=1,nbcnd
      if(BCSLD(nb)==1) then
        IBFS=LBC_IND(nb-1)+1
        IBFE=LBC_IND(nb)
        do IBFL=IBFS,IBFE
        IS1=LBC_SSF(IBFL)
        IS2=lbface(2,IS1)
        IC1=lcface(1,IS1)
        IC2=lcface(1,IS2)
        if(WK1(ic1)==0) then
           icL=icL+1
           WK1(ic1)=1
        endif
        if(WK1(ic2)==0) then
           icL=icL+1
           WK1(ic2)=1
        endif
        enddo
      endif
      ENDDO
!----------------------------
! --- discontinuous BC
!----------------------------
!     kdprdc BC
!        lbface(1,IS1)= nb
!        lbface(1,IS2)=-nb
!     kdsld BC
!        lbface(1,IS1)= nb
!        lbface(1,IS2)=-nb
!     kdintr BC
!        lbface(1,IS1)= nb
!        lbface(1,IS2)=-nb
!     kdtchi BC
!        lbface(1,IS1)= nb
!        lbface(1,IS2)= nb2
!----------------------------
      BCSLD(:)=0
      do 681 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      IBFS=LBC_IND(nb-1)+1
      IBFE=LBC_IND(nb)

      if(idis(nb)==1) then
        if(kd==kdsld.or.
     &     kd==kdprdc.or.
     &     kd==kdbuff.or.
     &     kd==kdshutr.or.
     &     kd==kdpors.or.
     &     kd==kdintr.or.
     &     kd==kdtchi) then
          nameBC2=boundName(nb,2)
          do IBFL=IBFS,IBFE
          IS1=LBC_SSF(IBFL)
          IC1=lcface(1,IS1)
          in1=lvface(1,IS1)
          in2=lvface(2,IS1)
          in3=lvface(3,IS1)
          in4=lvface(4,IS1)
          if(in4.ne.0) then
            icpu1=NODE_ID(in1,2)
            icpu2=NODE_ID(in2,2)
            icpu3=NODE_ID(in3,2)
            icpu4=NODE_ID(in4,2)
            if(icpu1.eq.icpu.or.
     &         icpu2.eq.icpu.or.
     &         icpu3.eq.icpu.or.
     &         icpu4.eq.icpu) then
               BCSLD(nb)=1
               if(kd==kdtchi) then
                 do nbx=1,nbcnd
                 if(nbx==nb) cycle
                 do j=1,2
                 nameBC1=boundName(nbx,j)
                 if(trim(nameBC2)==trim(nameBC1)) then
                   BCSLD(nbx)=1
                   GOTO 680
                 endif
                 enddo
                 enddo
               endif
               GOTO 680
            endif
          else
            icpu1=NODE_ID(in1,2)
            icpu2=NODE_ID(in2,2)
            icpu3=NODE_ID(in3,2)
            if(icpu1.eq.icpu.or.
     &         icpu2.eq.icpu.or.
     &         icpu3.eq.icpu) then
               BCSLD(nb)=1
               if(kd==kdtchi) then
                 do nbx=1,nbcnd
                 if(nbx==nb) cycle
                 do j=1,2
                 nameBC1=boundName(nbx,j)
                 if(trim(nameBC2)==trim(nameBC1)) then
                   BCSLD(nbx)=1
                   GOTO 680
                 endif
                 enddo
                 enddo
               endif
               GOTO 680
            endif
          endif
          enddo
        endif
      endif
 680  continue
 681  enddo
!
      DO nb=1,nbcnd
      if(BCSLD(nb)==1) then
        IBFS=LBC_IND(nb-1)+1
        IBFE=LBC_IND(nb)
        do IBFL=IBFS,IBFE
        IS1=LBC_SSF(IBFL)
        IC1=lcface(1,IS1)
        if(WK1(ic1)==0) then
           icL=icL+1
           WK1(ic1)=1
        endif
        enddo
      endif
      ENDDO
!
      ELMstack(icpu)=icL
!
 100  continue
!
!------------------------------------------------------
! --- 
!------------------------------------------------------
!
      do 130 icpu=1,PETOT
      ELMstack(icpu)=ELMstack(icpu-1)+ELMstack(icpu)
 130  continue
!
      nn=ELMstack(PETOT)
      allocate(ELMitem(nn))
!------------------------
! --- ELMitem(icLG)=icel
!------------------------
      ELMitem=0
!
      do 140 icpu=1,PETOT
      icLG=ELMstack(icpu-1)
      WK1=0
!
      do 160 ic=1,ICELTOT
      ityp=ICELTYP(ic)
! --- 
      if(ityp.eq.1) then
!
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
!
! --- 
!
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu)then
!
          icLG=icLG+1
          WK1(ic)=1
          ELMitem(icLG)=ic
        endif
      elseif(ityp.eq.2) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
!
! --- 
!
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu.or.
     &     icpu5.eq.icpu)then
!
          icLG=icLG+1
          WK1(ic)=1
          ELMitem(icLG)=ic
        endif
!
      elseif(ityp.eq.3) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
        icpu6=NODE_ID(in6,2)
!
! --- 
!
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu.or.
     &     icpu5.eq.icpu.or.
     &     icpu6.eq.icpu)then
!
          icLG=icLG+1
          WK1(ic)=1
          ELMitem(icLG)=ic
        endif

      elseif(ityp.eq.4) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        in7=lvcell(7,ic)
        in8=lvcell(8,ic)
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
        icpu6=NODE_ID(in6,2)
        icpu7=NODE_ID(in7,2)
        icpu8=NODE_ID(in8,2)
!
! --- 
!
        if(icpu1.eq.icpu.or.
     &     icpu2.eq.icpu.or.
     &     icpu3.eq.icpu.or.
     &     icpu4.eq.icpu.or.
     &     icpu5.eq.icpu.or.
     &     icpu6.eq.icpu.or.
     &     icpu7.eq.icpu.or.
     &     icpu8.eq.icpu)then
!
          icLG=icLG+1
          WK1(ic)=1
          ELMitem(icLG)=ic
        endif
      endif
! --- 
 160  continue
!--------------------------------
! --- cyclc BC and interface BC
!--------------------------------
      if(isufCYCLtot.gt.0) then
!        do 150 isuf=1,isufCYCLtot
        do 150 isuf=1,isufCYCLtot-1,2
        in1=iSUF_NODE(isuf,1)
        in2=iSUF_NODE(isuf,2)
        in3=iSUF_NODE(isuf,3)
        in4=iSUF_NODE(isuf,4)
        in5=ISUF_NODE(isuf+1,1)
        in6=ISUF_NODE(isuf+1,2)
        in7=ISUF_NODE(isuf+1,3)
        in8=ISUF_NODE(isuf+1,4)
        icel1=SUFtoCELL(isuf,1)
        icel2=SUFtoCELL(isuf,2)
!
        if(in4.ne.0) then
! --- 
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu4=NODE_ID(in4,2)
          icpu5=NODE_ID(in5,2)
          icpu6=NODE_ID(in6,2)
          icpu7=NODE_ID(in7,2)
          icpu8=NODE_ID(in8,2)
!
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu4.eq.icpu.or.
     &       icpu5.eq.icpu.or.
     &       icpu6.eq.icpu.or.
     &       icpu7.eq.icpu.or.
     &       icpu8.eq.icpu
     &   ) then
!
!
            if(WK1(icel1).eq.0) then
              icLG=icLG+1
              WK1(icel1)=1
              ELMitem(icLG)= icel1
            endif
            if(WK1(icel2).eq.0) then
              icLG=icLG+1
              WK1(icel2)=1
              ELMitem(icLG)= icel2
            endif
          endif
        else
! --- 
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu5=NODE_ID(in5,2)
          icpu6=NODE_ID(in6,2)
          icpu7=NODE_ID(in7,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu5.eq.icpu.or.
     &       icpu6.eq.icpu.or.
     &       icpu7.eq.icpu
     &    ) then
             if(WK1(icel1).eq.0) then
               icLG=icLG+1
               WK1(icel1)=1
               ELMitem(icLG)= icel1
             endif
             if(WK1(icel2).eq.0) then
               icLG=icLG+1
               WK1(icel2)=1
               ELMitem(icLG)= icel2
             endif
          endif
        endif
!
 150    enddo
      endif
!--------------------------------
! --- 
!--------------------------------
      if(isufPAIRtot.gt.0) then
        do 610 isufpair=1,isufPAIRtot
        in1=ISUF_NODE_pair(isufpair,1)
        in2=ISUF_NODE_pair(isufpair,2)
        in3=ISUF_NODE_pair(isufpair,3)
        in4=ISUF_NODE_pair(isufpair,4)
        icel1=SUFtoCELL_pair(isufpair,1)
        icel2=SUFtoCELL_pair(isufpair,2)
        icel3=SUFtoCELL_pair(isufpair,3)
        if(in4.ne.0) then
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu4=NODE_ID(in4,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu4.eq.icpu)then
            if(WK1(icel1).eq.0) then
              icLG=icLG+1
              WK1(icel1)=1
              ELMitem(icLG)= icel1
            endif
            if(WK1(icel2).eq.0) then
              icLG=icLG+1
              WK1(icel2)=1
              ELMitem(icLG)= icel2
            endif
            if(icel3.ne.0) then
              if(WK1(icel3).eq.0) then
                icLG=icLG+1
                WK1(icel3)=1
                ELMitem(icLG)= icel3
              endif
            endif
          endif
        else
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu) then
            if(WK1(icel1).eq.0) then
              icLG=icLG+1
              WK1(icel1)=1
              ELMitem(icLG)= icel1
            endif
            if(WK1(icel2).eq.0) then
              icLG=icLG+1
              WK1(icel2)=1
              ELMitem(icLG)= icel2
            endif
            if(icel3.ne.0) then
              if(WK1(icel3).eq.0) then
                icLG=icLG+1
                WK1(icel3)=1
                ELMitem(icLG)= icel3
              endif
            endif
          endif
        endif
 610    continue
      endif
!--------------------------------
! --- SLIDING BC
!--------------------------------
      BCSLD(:)=0
      do 670 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      IBFS=LBC_IND(nb-1)+1
      IBFE=LBC_IND(nb)
      if(kd.eq.kdsld.and.idis(nb)==0) then
!      if(kd.eq.kdsld) then
        do IBFL=IBFS,IBFE
        IS1=LBC_SSF(IBFL)
        IS2=lbface(2,IS1)
        IC1=lcface(1,IS1)
        IC2=lcface(1,IS2)
!
        in1=lvface(1,IS1)
        in2=lvface(2,IS1)
        in3=lvface(3,IS1)
        in4=lvface(4,IS1)
!
        im1=lvface(1,IS2)
        im2=lvface(2,IS2)
        im3=lvface(3,IS2)
        im4=lvface(4,IS2)
!
        if(in4.ne.0) then
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icpu4=NODE_ID(in4,2)
          icup1=NODE_ID(im1,2)
          icup2=NODE_ID(im2,2)
          icup3=NODE_ID(im3,2)
          icup4=NODE_ID(im4,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icpu4.eq.icpu.or.
     &       icup1.eq.icpu.or.
     &       icup2.eq.icpu.or.
     &       icup3.eq.icpu.or.
     &       icup4.eq.icpu)then
             BCSLD(nb)=1
             GOTO 670
          endif
        ELSE
          icpu1=NODE_ID(in1,2)
          icpu2=NODE_ID(in2,2)
          icpu3=NODE_ID(in3,2)
          icup1=NODE_ID(im1,2)
          icup2=NODE_ID(im2,2)
          icup3=NODE_ID(im3,2)
          if(icpu1.eq.icpu.or.
     &       icpu2.eq.icpu.or.
     &       icpu3.eq.icpu.or.
     &       icup1.eq.icpu.or.
     &       icup2.eq.icpu.or.
     &       icup3.eq.icpu)then
             BCSLD(nb)=1
             GOTO 670
          endif
        endif
        enddo
      endif
 670  continue
!
      DO nb=1,nbcnd
      if(BCSLD(nb)==1) then
        IBFS=LBC_IND(nb-1)+1
        IBFE=LBC_IND(nb)
        do IBFL=IBFS,IBFE
        IS1=LBC_SSF(IBFL)
        IS2=lbface(2,IS1)
        IC1=lcface(1,IS1)
        IC2=lcface(1,IS2)
        if(WK1(ic1)==0) then
          icLG=icLG+1
          WK1(ic1)=1
          ELMitem(icLG)=ic1
        endif
        if(WK1(ic2)==0) then
          icLG=icLG+1
          WK1(ic2)=1
          ELMitem(icLG)=ic2
        endif
        enddo
      endif
      ENDDO
!
!----------------------------
! --- discontinuous BC
!----------------------------
      BCSLD(:)=0
      do 690 nb=1,nbcnd
      kd=kdbcnd(0,nb)
      IBFS=LBC_IND(nb-1)+1
      IBFE=LBC_IND(nb)
      if(idis(nb)==1) then
        if(kd==kdsld.or.
     &     kd==kdprdc.or.
     &     kd==kdbuff.or.
     &     kd==kdshutr.or.
     &     kd==kdpors.or.
     &     kd==kdtchi.or.
     &     kd==kdintr) then
          nameBC1=boundName(nb,1)
          nameBC2=boundName(nb,2)
          do IBFL=IBFS,IBFE
          IS1=LBC_SSF(IBFL)
          IC1=lcface(1,IS1)
          in1=lvface(1,IS1)
          in2=lvface(2,IS1)
          in3=lvface(3,IS1)
          in4=lvface(4,IS1)
          if(in4.ne.0) then
            icpu1=NODE_ID(in1,2)
            icpu2=NODE_ID(in2,2)
            icpu3=NODE_ID(in3,2)
            icpu4=NODE_ID(in4,2)
            if(icpu1.eq.icpu.or.
     &         icpu2.eq.icpu.or.
     &         icpu3.eq.icpu.or.
     &         icpu4.eq.icpu) then
               BCSLD(nb)=1
               if(kd==kdtchi) then
                 do nbx=1,nbcnd
                 if(nbx==nb) cycle
                 do j=1,2
                 nameBC1=boundName(nbx,j)
                 if(trim(nameBC2)==trim(nameBC1)) then
                   BCSLD(nbx)=1
                   GOTO 690
                 endif
                 enddo
                 enddo
               endif
               GOTO 690
            endif
          else
            icpu1=NODE_ID(in1,2)
            icpu2=NODE_ID(in2,2)
            icpu3=NODE_ID(in3,2)
            if(icpu1.eq.icpu.or.
     &         icpu2.eq.icpu.or.
     &         icpu3.eq.icpu) then
               BCSLD(nb)=1
               if(kd==kdtchi) then
                 do nbx=1,nbcnd
                 if(nbx==nb) cycle
                 do j=1,2
                 nameBC1=boundName(nbx,j)
                 if(trim(nameBC2)==trim(nameBC1)) then
                   BCSLD(nbx)=1
                   GOTO 690
                 endif
                 enddo
                 enddo
               endif
               GOTO 690
            endif
          endif
          enddo
        endif
      endif
 690  enddo
!
      DO nb=1,nbcnd
      if(BCSLD(nb)==1) then
        IBFS=LBC_IND(nb-1)+1
        IBFE=LBC_IND(nb)
        do IBFL=IBFS,IBFE
        IS1=LBC_SSF(IBFL)
        IC1=lcface(1,IS1)
        if(WK1(ic1)==0) then
           icLG=icLG+1
           WK1(ic1)=1
           ELMitem(icLG)=ic1
        endif
        enddo
      endif
      enddo
!
 140  continue
!
      deallocate (WK1)
!
! === 
!
      write (ifll,'(a)') '  - NEIGHBORING PE'
!
! --- +----------------+
! --- | NEIGHBORING PE |
! --- +----------------+
!
      allocate 
     & (NEIBPETOT(PETOT),NEIBPE(PETOT,PETOT),NEIBPEinv(PETOT,PETOT))
!
      NEIBPETOT=0
      NEIBPE   =0
      NEIBPEinv=0
!
      do 170 icpu=1,PETOT
!-------------------------------------------------------------------
! --- | Local ID of icpu from ELMstack(icpu-1)+1 to ELMstack(icpu) |
!-------------------------------------------------------------------
      do 180 icLG=ELMstack(icpu-1)+1, ELMstack(icpu)
!
      ic=ELMitem(icLG)
      ityp=ICELTYP(ic)
      
! --- 
      if(ityp.eq.1) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
!
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
!
        if(icpu1.ne.icpu) then
          NEIBPEinv(icpu ,icpu1)=1
          NEIBPEinv(icpu1,icpu )=1
        endif
!
        if(icpu2.ne.icpu) then
          NEIBPEinv(icpu ,icpu2)=1
          NEIBPEinv(icpu2,icpu )=1
        endif
!
        if(icpu3.ne.icpu) then
          NEIBPEinv(icpu ,icpu3)=1
          NEIBPEinv(icpu3,icpu )=1
        endif
!
        if(icpu4.ne.icpu) then
          NEIBPEinv(icpu ,icpu4)=1
          NEIBPEinv(icpu4,icpu )=1
        endif
      elseif(ityp.eq.2) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
!
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
!
        if(icpu1.ne.icpu) then
          NEIBPEinv(icpu ,icpu1)=1
          NEIBPEinv(icpu1,icpu )=1
        endif
!
        if(icpu2.ne.icpu) then
          NEIBPEinv(icpu ,icpu2)=1
          NEIBPEinv(icpu2,icpu )=1
        endif
!
        if(icpu3.ne.icpu) then
          NEIBPEinv(icpu ,icpu3)=1
          NEIBPEinv(icpu3,icpu )=1
        endif
!
        if(icpu4.ne.icpu) then
          NEIBPEinv(icpu ,icpu4)=1
          NEIBPEinv(icpu4,icpu )=1
        endif
!
        if(icpu5.ne.icpu) then
          NEIBPEinv(icpu ,icpu5)=1
          NEIBPEinv(icpu5,icpu )=1
        endif
      elseif(ityp.eq.3) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
!
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
        icpu6=NODE_ID(in6,2)
!
        if(icpu1.ne.icpu) then
          NEIBPEinv(icpu ,icpu1)=1
          NEIBPEinv(icpu1,icpu )=1
        endif
!
        if(icpu2.ne.icpu) then
          NEIBPEinv(icpu ,icpu2)=1
          NEIBPEinv(icpu2,icpu )=1
        endif
!
        if(icpu3.ne.icpu) then
          NEIBPEinv(icpu ,icpu3)=1
          NEIBPEinv(icpu3,icpu )=1
        endif
!
        if(icpu4.ne.icpu) then
          NEIBPEinv(icpu ,icpu4)=1
          NEIBPEinv(icpu4,icpu )=1
        endif
        if(icpu5.ne.icpu) then
          NEIBPEinv(icpu ,icpu5)=1
          NEIBPEinv(icpu5,icpu )=1
        endif
        if(icpu6.ne.icpu) then
          NEIBPEinv(icpu ,icpu6)=1
          NEIBPEinv(icpu6,icpu )=1
        endif
      elseif(ityp.eq.4) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        in7=lvcell(7,ic)
        in8=lvcell(8,ic)
!
        icpu1=NODE_ID(in1,2)
        icpu2=NODE_ID(in2,2)
        icpu3=NODE_ID(in3,2)
        icpu4=NODE_ID(in4,2)
        icpu5=NODE_ID(in5,2)
        icpu6=NODE_ID(in6,2)
        icpu7=NODE_ID(in7,2)
        icpu8=NODE_ID(in8,2)
!
        if(icpu1.ne.icpu) then
          NEIBPEinv(icpu ,icpu1)=1
          NEIBPEinv(icpu1,icpu )=1
        endif
!
        if(icpu2.ne.icpu) then
          NEIBPEinv(icpu ,icpu2)=1
          NEIBPEinv(icpu2,icpu )=1
        endif
!
        if(icpu3.ne.icpu) then
          NEIBPEinv(icpu ,icpu3)=1
          NEIBPEinv(icpu3,icpu )=1
        endif
!
        if(icpu4.ne.icpu) then
          NEIBPEinv(icpu ,icpu4)=1
          NEIBPEinv(icpu4,icpu )=1
        endif
        if(icpu5.ne.icpu) then
          NEIBPEinv(icpu ,icpu5)=1
          NEIBPEinv(icpu5,icpu )=1
        endif
        if(icpu6.ne.icpu) then
          NEIBPEinv(icpu ,icpu6)=1
          NEIBPEinv(icpu6,icpu )=1
        endif
        if(icpu7.ne.icpu) then
          NEIBPEinv(icpu ,icpu7)=1
          NEIBPEinv(icpu7,icpu )=1
        endif
        if(icpu8.ne.icpu) then
          NEIBPEinv(icpu ,icpu8)=1
          NEIBPEinv(icpu8,icpu )=1
        endif
      endif
!
 180  continue
!
 170  continue
!
! --- 
!
      do 190 icpu1=1,PETOT
      ncpu=0
      do 195 icpu2=1,PETOT
      if(NEIBPEinv(icpu1,icpu2).ne.0) then
        ncpu=ncpu+1
        NEIBPE(icpu1,ncpu)=icpu2
      endif
 195  continue
      NEIBPETOT(icpu1)=ncpu
 190  continue
!
      NEIBPEinv= 0
!
! --- 
!
      do 200 icpu1=1,PETOT
      do 210 ncpu=1,NEIBPETOT(icpu1)
      icpu2=NEIBPE(icpu1,ncpu)
      do 220 ncpuc=1,NEIBPETOT(icpu2)
      if(NEIBPE(icpu2,ncpuc).eq.icpu1) then
        NEIBPEinv(icpu1,ncpu)=ncpuc
      endif
 220  continue
 210  continue
 200  continue
!
! --- +---------------------------+
! --- | LOCAL NODE & IMPORT TABLE |
! --- +---------------------------+
!
      allocate (NODstack(0:PETOT))
      NODstack=0
! --- 
      allocate (WK1(INODTOT))
      allocate (WK2(INODTOT))
! --- 
      do 300 icpu=1,PETOT
      WK1=0
      do 240 icLG=ELMstack(icpu-1)+1, ELMstack(icpu)
      ic=ELMitem(icLG)
      ityp=ICELTYP(ic)
      if(ityp.eq.1) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        WK1(in1)=1
        WK1(in2)=1
        WK1(in3)=1
        WK1(in4)=1
      elseif(ityp.eq.2) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        WK1(in1)=1
        WK1(in2)=1
        WK1(in3)=1
        WK1(in4)=1
        WK1(in5)=1
      elseif(ityp.eq.3) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        WK1(in1)=1
        WK1(in2)=1
        WK1(in3)=1
        WK1(in4)=1
        WK1(in5)=1
        WK1(in6)=1
      elseif(ityp.eq.4) then
! --- 
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        in7=lvcell(7,ic)
        in8=lvcell(8,ic)
        WK1(in1)=1
        WK1(in2)=1
        WK1(in3)=1
        WK1(in4)=1
        WK1(in5)=1
        WK1(in6)=1
        WK1(in7)=1
        WK1(in8)=1
      endif
! --- 
 240  continue
!-----+----------------------------------------+
! --- | NODtotL(icpu) : total inter vertex no. |
!-----+----------------------------------------+
      WK2=0
      ivexL=0
      do 250 iv=1,INODTOT
      if(WK1(iv).eq.1.and.NODE_ID(iv,2).ne.icpu) then
        ivexL=ivexL+1
        WK2(ivexL)=iv
      endif
 250  continue
!-----+--------------------------------+
! --- | NODstack(icpu): inter + exteri |
!-----+--------------------------------+
      NODstack(icpu)=ivexL+NODtotL(icpu)
!
      open(11,file=WORKFIL(icpu),status='unknown',form= 'unformatted',
     &     iostat = OpenStatus)
      if(OpenStatus > 0) then
        write (ifll,'(a)') '*** Cannot open file: Please confirm that 
     & the HPC control sub-directries(ie. hpc_0000) exist.'
        STOP
      endif
      rewind (11)
      write  (11)  ivexL
      write  (11)  (WK2(j),j=1,ivexL)
      close  (11)
 300  continue
!
! --- make array: NODitem(ivLG)=iv
!
      do 310 icpu=1,PETOT
      NODstack(icpu)=NODstack(icpu-1)+NODstack(icpu)
 310  continue
!------------------------------------------
! --- | NODitem(ivLG)=iv; interior-vertex
!------------------------------------------
      nn=NODstack(PETOT)
      allocate(NODitem(nn))
!
! --- 
!
      do 330 icpu=1,PETOT
      ivLG=NODstack(icpu-1)
      do 320 iv=1,INODTOT
      if(NODE_ID(iv,2).eq.icpu) then
        ivLG=ivLG+1
        NODitem(ivLG)=iv
      endif
 320  continue
!---------------------------------------------
! --- (NODitem(ivLG)=iv: interior+1=>exterior
!---------------------------------------------
      open(11,file=WORKFIL(icpu),status='unknown', 
     &                           form= 'unformatted')
      read (11) ii
      read (11) (NODitem(j),j=ivLG+1,NODstack(icpu))
      close (11)
 330  continue
!
      deallocate(WK1,WK2)
!
      allocate(WK1(0:PETOT),WK2(PETOT),WK3(PETOT))
!
! --- 
!
      do 400 icpu=1,PETOT
      WK1=0
      WK2=0
      do 410 ncpu=1,NEIBPETOT(icpu)
      icpu2=NEIBPE(icpu,ncpu)
      WK2(icpu2)=ncpu
 410  continue
!----------------------------------------------------------------
! --- plus exterior vertex: NODstack(icpu-1) + NODtotL(icpu)
! --- WK1(ncpu): exterior part's pointer
!----------------------------------------------------------------
      ivLG=NODstack(icpu-1)+NODtotL(icpu)
!
      do 420 ivexLG=ivLG+1,NODstack(icpu)
      iv=NODitem(ivexLG)
      ncpu=WK2(NODE_ID(iv,2))
      WK1(ncpu)=WK1(ncpu)+1
 420  continue
!
      do 430 ncpu=1,NEIBPETOT(icpu)
      WK1(ncpu)=WK1(ncpu-1)+WK1(ncpu)
 430  continue
!---------------------------------------------
! --- total exterior vertex on icpu: nn
!---------------------------------------------
      nn=WK1(NEIBPETOT(icpu))
      allocate (WK4(nn))
! 
      WK3=0
      WK4=0
!---------------------------------------------
! --- interior vertex; ivexLG
!---------------------------------------------
      do 440 ivexLG=ivLG+1,NODstack(icpu)
      iv=NODitem(ivexLG)
      ncpu=WK2(NODE_ID(iv,2))
      WK3(ncpu)=WK3(ncpu)+1
      ivL=ivexLG-NODstack(icpu-1)
      ivexLGS=WK1(ncpu-1)+WK3(ncpu)
      WK4(ivexLGS)=ivL
 440  continue
!
      open (11,file=WORKFIL(icpu),status='unknown',
     &                          form= 'unformatted')
      rewind(11)
      write(11)  NEIBPETOT(icpu)
      write(11) (WK1(ncpu),ncpu=1,NEIBPETOT(icpu))
      write(11) (WK4(ivexLGS),ivexLGS=1,WK1(NEIBPETOT(icpu)))
      WK1(0)=0
      close (11)        
!
      deallocate (WK4)
!
 400  continue
!
      deallocate(WK1,WK2,WK3)
!===
!
      write (ifll,'(a)') '  - EXPORT'
!
! --- +--------------+
! --- | EXPORT TABLE |
! --- +--------------+
!
      allocate (WK1(0:PETOT), WK2(0:PETOT))
!
! --- +---------------------------------------------------------------+
!
      do 500 icpu=1,PETOT
!
      open(11,file=WORKFIL(icpu),status='unknown',  
     &                        form= 'unformatted')
      read(11) ncpu

      allocate(IMPORT_index(0:ncpu))
               IMPORT_index= 0

      read(11) (IMPORT_index(j),j=1,ncpu)
      IMPORT_index(0)=0
      ivexLGS=IMPORT_index(ncpu)

      allocate(IMPORT_item(ivexLGS))

      read(11)(IMPORT_item(j),j=1,ivexLGS)
      close (11)
!
! --- 
!
      WK2=0
      do 510 ncpu=1,NEIBPETOT(icpu)
      icpuc=NEIBPE(icpu,ncpu)
      ncpuc=NEIBPEinv(icpu,ncpu)
! --- 
      open (11, file= WORKFIL(icpuc), status='unknown',
     &                              form= 'unformatted')
      read(11) jj
      read(11) (WK1(j),j=1,jj)
      WK1(0)=0
      WK2(ncpu)=WK1(ncpuc)-WK1(ncpuc-1)
      close (11)
 510  continue
!
      do 520 ncpu=1,NEIBPETOT(icpu)
      WK2(ncpu)=WK2(ncpu-1)+WK2(ncpu)
 520  continue
!
      allocate (EXPORT_index(0:NEIBPETOT(icpu)))
      allocate (EXPORT_item(WK2(NEIBPETOT(icpu))))
      EXPORT_index=0
      EXPORT_item =0
!
      do 530 ncpu=1,NEIBPETOT(icpu)
      icpuc=NEIBPE(icpu,ncpu)
      ncpuc=NEIBPEinv(icpu,ncpu)
      open (11,file=WORKFIL(icpuc),status='unknown',   
     &                          form= 'unformatted')
      read(11) jj
      read(11) (WK1(j),j=1,jj)
      WK1(0)=0
      EXPORT_index(ncpu)=WK1(ncpuc)-WK1(ncpuc-1)
!
      allocate (WK3(WK1(jj)))
      read(11) (WK3(j),j=1,WK1(jj))
!
      jS=WK1(ncpuc-1)+1
      jE=WK1(ncpuc)
!
      do 540 ivexLGS=WK1(ncpuc-1)+1,WK1(ncpuc)
      ivinLGS=ivexLGS-WK1(ncpuc-1)+WK2(ncpu-1)
      ivL=WK3(ivexLGS)
      ivexLG=NODstack(icpuc-1)+ivL
      iv=NODitem(ivexLG)
      ivinL=NODE_ID(iv,1)
      EXPORT_item(ivinLGS)=ivinL
 540  continue
!
      deallocate(WK3)
      close(11)
 530  continue
!
      do 550 ncpu=1,NEIBPETOT(icpu)
      EXPORT_index(ncpu)=EXPORT_index(ncpu-1)+EXPORT_index(ncpu)
 550  continue
!----------------------------
! --- COMMout(icpu)
!----------------------------
      open (11,file=COMMout(icpu),status= 'unknown',form='formatted')
      kmax1= IMPORT_index(NEIBPETOT(icpu))
      kmax2= EXPORT_index(NEIBPETOT(icpu))
      write (11,'(a)')   '#NEIBPEtot'
      write (11,'(6i12)') NEIBPETOT(icpu)
      write (11,'(a)')   '#NEIBPE'         
      write (11,'(6i12)') (NEIBPE(icpu,ncpu)-1,ncpu=1,NEIBPETOT(icpu))  !-1
      write (11,'(a)')   '#IMPORT index'
      write (11,'(6i12)') (IMPORT_index(ncpu),ncpu=1,NEIBPETOT(icpu))
      write (11,'(a)')   '#IMPORT items'
      write (11,'(6i12)') (IMPORT_item (ncpu),ncpu=1,kmax1)
      write (11,'(a)')   '#EXPORT index'
      write (11,'(6i12)') (EXPORT_index(ncpu),ncpu=1,NEIBPETOT(icpu))
      write (11,'(a)')   '#EXPORT items'
      write (11,'(6i12)') (EXPORT_item (ncpu),ncpu=1,kmax2)
      write (11,'(a)')   '#INTERNAL NODE'
      write (11,'(6i12)')  NODtotL(icpu)
      write (11,'(a)')   '#TOTAL NODE (NO WALL)'
      write (11,'(6i12)')  NODstack(icpu)-NODstack(icpu-1)
      write (11,'(a)')   '#GLOBAL NODE ID (NO WALL)'
      write (11,'(6i12)') (NODitem(i),
     &                     i=NODstack(icpu-1)+1,NODstack(icpu))
      close (11)
!
      deallocate (IMPORT_index, IMPORT_item)
      deallocate (EXPORT_index, EXPORT_item)
!
 500  continue
!
      deallocate (WK1,WK2,LBC_SSF)
!
      return
!
      end subroutine COMM
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine re_comm(mssfbc)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_boundary,only    : nbcnd
      use module_hpc,only         : NPETOT
      use module_io,only          : ifle,ifll
      use module_partitioner,only : LCV_CV => WK1
      use module_partitioner,only : MAT_CV => WK5
      use module_partitioner,only : NOD_IDg=> WK3
      use module_partitioner,only : MAT_CVG=> WK4
      use module_partitioner,only : DISALL => WALLDIST
      use module_partitioner,only : LFUTAU => wk6
      use module_material,only    : ical_sld,rotati,ishaft
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: mssfbc
!
! --- [local entities]
!
      integer :: my_rank,ifl,icpu,ifchl,kmax1,kmax2,k,ios
      integer :: INTNODTOT,NODTOT,ICVL,ICV,NEIBPETOT_S
      integer :: NCVX,ierr=0
      integer :: ierr1=0
      character(len=80),save :: fname
!-----------------------
! --- Serial Resorting
!-----------------------
      open(11,file='MAT_CVG',FORM='unformatted',status='unknown',
     &        iostat=ios)
      read (11) NCVX
      allocate (MAT_CVG(NCVX))
      read (11)(MAT_CVG(ICV),ICV=1,NCVX)
      CLOSE(11)
!      if(ical_sldx==1) then
!        allocate(LSLDG(mssfbc,2),stat=ierr1)
!        if(ierr1.ne.0) 
!     &  stop 'stop at allocating LSLDG in list_output_geom'
!        nbsldx=0
!        nsldx(:)=0
!        LSLDG(:,:)=0
!        read(11) nbsldx,(nsldx(nbx),nbx=1,nbsldx)
!        do nbx=1,nbsldx
!        read(11) (LSLDG(isd,1),LSLDG(isd,2),isd=1,nsldx(nbx))
!        enddo
!        deallocate(LSLDG)
!      endif
!      close(11)
!--------------------
! --- 
!--------------------
!      open (11,FILE='distance',FORM='UNFORMATTED',
!     &            action='read',iostat=ios)
!      read (11) NCVX
!      allocate (DISALL(NCVX),LFUTAU(NCVX),stat=ierr)
!      if(ierr/=0) stop 'stop at allocating DISALL in '
!      read (11,ERR=1200) (DISALL(MAT_CVG(ICVL)),ICVL=1,NCVX)
!      read (11,ERR=1200) (LFUTAU(MAT_CVG(ICVL)),ICVL=1,NCVX)
!      close(11)
!---------------------------------
! --- Reading file: RODRICV(icpu)
!---------------------------------
!      if(ical_sld==1) then
!        nslid(:)=0
!        do icpu=1,NPETOT
!        do nbx=1,nsldbc(icpu)
!          nslid(icpu)=nslid(icpu)+nsldT(nbx,icpu)
!        enddo
!        enddo
!        do icpu=1,NPETOT
!        nslid(icpu)=nslid(icpu-1)+nslid(icpu)
!        enddo
!        allocate(LSLDL(nslid(NPETOT),2),stat=ierr1)
!        if(ierr1.ne.0)
!     &  stop 'stop at allocating LSLDL in re_comm'
!        LSLDL(:,:)=0
!      endif
!
      do 100 icpu=1,NPETOT
      allocate(LCV_CV(NALLCV_hpc(icpu)))
      allocate(MAT_CV(NALLCV_hpc(icpu)))
!
      ifl=50
      open(ifl,file=RODRICV(icpu),FORM='unformatted',status='unknown',
     &            iostat=ios)
      read(ifl) (LCV_CV(ICV),ICV=1,NALLCV_hpc(icpu))
      read(ifl) (MAT_CV(ICV),ICV=1,NALLCV_hpc(icpu))
!      if(ical_sldx==1) then
!        nbsldx=0
!        nsldx(:)=0
!        read(ifl) nbsldx,(nsldx(nbx),nbx=1,nbsldx)
!        do nbx=1,nbsldx
!        read(ifl) (LSLDL(isd,1),LSLDL(isd,2),
!     &             isd=nslid(icpu-1)+1,nsldx(nbx)+nslid(icpu-1))
!        enddo
!      endif
      close(ifl)
!---------------------------------------
! --- Reading comm file: COMMout(icpu)
!---------------------------------------
      ifchl=11
      open(ifchl,file=COMMout(icpu),status= 'unknown',form='formatted')
      read (ifchl,'(a80)' ,ERR=98) LINE
      read (ifchl,'(6i12)',ERR=98) NEIBPETOT_S
      read (ifchl,'(a80)' ,ERR=98) LINE
      allocate(WK2(NEIBPETOT_S))
      read (ifchl,'(6i12)',ERR=98) (WK2(k),k=1,NEIBPETOT_S)
      read (ifchl,'(a80)' ,ERR=98) LINE
      allocate(IMPORT_index(NEIBPETOT_S))
      read (ifchl,'(6i12)',ERR=98) (IMPORT_index(k),k=1,NEIBPETOT_S)
      read (ifchl,'(a80)' ,ERR=98) LINE
      kmax1= IMPORT_index(NEIBPETOT_S)
      allocate(IMPORT_item(kmax1))
      read (ifchl,'(6i12)',ERR=98) (IMPORT_item(k),k=1,kmax1)
      read (ifchl,'(a80)' ,ERR=98) LINE
      allocate(EXPORT_index(NEIBPETOT_S))
      read (ifchl,'(6i12)',ERR=98) (EXPORT_index(k),k=1,NEIBPETOT_S)
      kmax2= EXPORT_index(NEIBPETOT_S)
      allocate(EXPORT_item(kmax2))
      read (ifchl,'(a80)' ,ERR=98) LINE
      read (ifchl,'(6i12)',ERR=98) (EXPORT_item(k),k=1,kmax2)
      read (ifchl,'(a80)' ,ERR=98) LINE
      read (ifchl,'(6i12)',ERR=98) INTNODTOT
      read (ifchl,'(a80)' ,ERR=98) LINE
      read (ifchl,'(6i12)',ERR=98) NODTOT
      read (ifchl,'(a80)' ,ERR=98) LINE
      allocate(NOD_IDg(NODTOT))
      read (ifchl,'(6i12)',ERR=98) (NOD_IDg(k),k=1,NODTOT)
      close(ifchl)
!----------------------------------------------
! --- Rewriteing comm file: COMMout(icpu)
!----------------------------------------------
      open (ifchl,file=COMMout(icpu),status= 'unknown',form='formatted')
      kmax1= IMPORT_index(NEIBPETOT_S)
      kmax2= EXPORT_index(NEIBPETOT_S)
      write (ifchl,'(a)')   '#NEIBPEtot'
      write (ifchl,'(6i12)') NEIBPETOT_S
      write (ifchl,'(a)')   '#NEIBPE'         
      write (ifchl,'(6i12)') (WK2(k),k=1,NEIBPETOT_S)
      write (ifchl,'(a)')   '#IMPORT index'
      write (ifchl,'(6i12)') (IMPORT_index(k),k=1,NEIBPETOT_S)
      write (ifchl,'(a)')   '#IMPORT items'
      write (ifchl,'(6i12)') (LCV_CV(IMPORT_item(k)),k=1,kmax1)
      write (ifchl,'(a)')   '#EXPORT index'
      write (ifchl,'(6i12)') (EXPORT_index(k),k=1,NEIBPETOT_S)
      write (ifchl,'(a)')   '#EXPORT items'
      write (ifchl,'(6i12)') (LCV_CV(EXPORT_item(k)),k=1,kmax2)
      write (ifchl,'(a)')   '#INTERNAL NODE'
      write (ifchl,'(6i12)')  INTNODTOT
      write (ifchl,'(a)')   '#TOTAL NODE (NO WALL)'
      write (ifchl,'(6i12)')  NODTOT
      write (ifchl,'(a)')   '#GLOBAL NODE ID (NO WALL)'
      write (ifchl,'(6i12)') (NOD_IDg(ICV),ICV=1,NODTOT)
      close(ifchl)
!-----------------------------------------
! --- HPC wall distance file :'wall_hpc'
!-----------------------------------------
!      my_rank=icpu-1
!      HEADER1='wall_hpc'
!      call DEFINE_FILE_NAME(DIRHED,HEADER1,my_rank,walldis(icpu))
!      fname=walldis(icpu)
!      open (11,FILE=fname,FORM='UNFORMATTED',
!     &            action='write',iostat=ios)
!      write(11) (DISALL(NOD_IDg(MAT_CV(ICVL))),ICVL=1,NCV_hpc(icpu))
!      close(11)
!-----------------
! --- deallocate
!-----------------
      deallocate(WK2,LCV_CV,MAT_CV,
     &   IMPORT_index,IMPORT_item,EXPORT_index,EXPORT_item,NOD_IDg)
!
 100  continue
!
!      deallocate(DISALL,LFUTAU,MAT_CVG)
      deallocate(MAT_CVG)
!
      return
 98   stop 'read comm file'
 1200 stop ' distance file read error in re_comm'
      end subroutine re_comm
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
