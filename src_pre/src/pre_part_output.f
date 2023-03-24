!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine OUTPUT_hpc_grid
     &           (nface,ncelb,NCV,nedge,nssfbc,NCVFAC,NALLCV,
     &            ncell,mcell,mface,medge,mvrtx,
     &            lfcell,lcface,LEFACE,lvface,cord,lvcell,lacell,
     &            lbface)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_io,only       : ifli,ifll,ifle,gdformat
      use module_boundary,only : nbcnd,nobcnd,ivbcnd,lvbcnd,kdbcnd,
     &                           icbinr,lcbinr,
     &                           kdprdc,kxnone,kdsymm,kdtchi,kdsld,
     &                           kdintr,kdilet,kdolet,undef_bcno,
     &                           boundIDMap,SFBOUN,NBOUND,idis,kdbuff,
     &                           ivpair,kdshutr,kdpors
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: nface,ncelb,NCV,nedge,nssfbc,NCVFAC,NALLCV
      integer,intent(in)  :: ncell
      integer,intent(in)  :: mcell,mface,medge,mvrtx
      integer,intent(in)  :: lfcell(7,mcell)
      integer,intent(in)  :: lcface(2,mface)
      integer,intent(in)  :: LEFACE(5,MFACE)
      integer,intent(in)  :: lvface(4,mface)
      real*8 ,intent(in)  :: cord  (3,mvrtx)
      integer,intent(in)  :: lvcell(8,mcell)
      integer,intent(in)  :: lacell(  mcell)
      integer,intent(in)  :: lbface(2,mface)
!
! --- [local entities]
!
      integer :: ic,icL,ivL,ivLG,icpu,ib,ibin,ibwl,ibsm,ibcy
      integer :: ie,ibcy1,ibcy2,ip0,ip1,ips,ipe,ivv,iv
      integer :: nn0,nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,NN
      integer :: icpu1,icpu2,icpu3,icpu4
      integer :: icpu11,icpu12,icpu13,icpu14
      integer :: icpu21,icpu22,icpu23,icpu24
      integer :: iv11,iv12,iv13,iv14,iv21,iv22,iv23,iv24
      integer :: iv31,iv32,iv33,iv34
      integer :: in1,in2,in3,in4,in5,in6,in7,in8
      integer :: iv1,iv2,iv3,iv4
      integer :: NNX,NNA,NNB,i,j,k,N1,N0,IS,ILV
      integer :: nbufa,nbufb,nb,mb,icLG,ierr1=0
      integer :: ityp,NSUF_LOCAL,isufp,isufo,nob,mb1,mb2,kd,isuf,
     &           icel1,icel2
      integer :: in11,in12,in13,in14,in15,in16,in17,in18
      integer :: in21,in22,in23,in24,in25,in26,in27,in28
      integer :: iwk1,iwk2,IS2
!
      real*8, allocatable :: XYZ(:,:)
!
! --- 
!
      in1=0
      in2=0
      in3=0
      in4=0
      in5=0
      in6=0
      in7=0
      in8=0

      nn0=0
      nn1=1
      nn2=2
      nn3=3
      nn4=4
      nn5=5
      nn6=6
      nn7=7
      nn8=8
!
      HEADER1='NODE'
!
      allocate (WK2(0:INODTOT),stat=ierr1)
      WK2=0
      if(ierr1.ne.0) stop 'stop at allocation -1- in OUTPUT'
!
! --- DO-LOOP 100 START
!
      do 100 icpu=1,PETOT
!
      HEADER2='ELEMENT'
!
! --- +-----------------+
! --- | WRITE MESH file |
! --- +-----------------+
!
      open (11,file=GRIDout(icpu), status= 'unknown')
!
      write(11,'(a)') TRIM(adjustL(gdformat))
      write(11,'(a)') TRIM(adjustL(HEADER1))
!
! --- Total Local Nodes (Inter, exter, wall)
!
      NNX=NODstack(icpu)-NODstack(icpu-1)
      NNA=ELMstack(icpu)-ELMstack(icpu-1)
      NNB=NODtotL(icpu)
      nvrtx_hpc(icpu)=NNX
      allocate (XYZ(NNX,3),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -2- in OUTPUT'
!
      XYZ=0.d0
      ivL=0
! --- vertex coordinates
      do 110 j=NODstack(icpu-1)+1,NODstack(icpu)
      ivL=ivL+1
      iv=NODitem(j)
      XYZ(ivL,1)=cord(1,iv)
      XYZ(ivL,2)=cord(2,iv)
      XYZ(ivL,3)=cord(3,iv)
 110  continue
!
      write(11,'(6i12)') nn3, NNX
      write(11,'(5(1pe14.6))') ((XYZ(i,k),k=1,3),i=1,NNX)
      deallocate (XYZ)
!----------------------------------------
! --- Avoiding Same gloabl vertex number
!----------------------------------------
      WK2=0
      ivL=0
      do 111 j=NODstack(icpu-1)+1,NODstack(icpu)
      iv=NODitem(j)
      if(WK2(iv).eq.0) then
        ivL=ivL+1
        WK2(iv)=ivL
      else
        write(ifle,*) ' ### ERROR-1 : One Vertex appears twice in icpu'
        write(ifle,*) ' ICPU = ',icpu
        write(ifle,*) NNB,NNA,NNX,WK2(iv),ivL,iv
        stop
      endif
 111  continue
      if(ivL.ne.NNX) then
        write(ifle,*) ' ### ERROR-2 :',ivL,NNX
        stop
      endif
!----------------------------------------
! --- Local CONNECTIVITY
!----------------------------------------
      write(11,'(a)') TRIM(adjustL(HEADER2))
!
      NNX=ELMstack(icpu)-ELMstack(icpu-1)
      allocate(IW1(NNX,8),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -3- in OUTPUT'
      ncell_hpc(icpu)=NNX
      write(ifle,*) '    ###    NCELL_H= ',NCELL_HPC(ICPU)
      IW1=0
      icL=0
      do 120 icLG=ELMstack(icpu-1)+1,ELMstack(icpu)
      icL=icL+1
      ic=ELMitem(icLG)
      ityp=ICELTYP(ic)
      if(ityp.eq.1) then
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        IW1(icL,1)=WK2(in1)
        IW1(icL,2)=WK2(in2)
        IW1(icL,3)=WK2(in3)
        IW1(icL,4)=WK2(in4)
        IW1(icL,5)=0
        IW1(icL,6)=0
        IW1(icL,7)=0
        IW1(icL,8)=0
      elseif(ityp.eq.2) then
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        IW1(icL,1)=WK2(in1)
        IW1(icL,2)=WK2(in2)
        IW1(icL,3)=WK2(in3)
        IW1(icL,4)=WK2(in4)
        IW1(icL,5)=WK2(in5)
        IW1(icL,6)=0
        IW1(icL,7)=0
        IW1(icL,8)=0
      elseif(ityp.eq.3) then
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        IW1(icL,1)=WK2(in1)
        IW1(icL,2)=WK2(in2)
        IW1(icL,3)=WK2(in3)
        IW1(icL,4)=WK2(in4)
        IW1(icL,5)=WK2(in5)
        IW1(icL,6)=WK2(in6)
        IW1(icL,7)=0
        IW1(icL,8)=0
      elseif(ityp.eq.4) then
        in1=lvcell(1,ic)
        in2=lvcell(2,ic)
        in3=lvcell(3,ic)
        in4=lvcell(4,ic)
        in5=lvcell(5,ic)
        in6=lvcell(6,ic)
        in7=lvcell(7,ic)
        in8=lvcell(8,ic)
        IW1(icL,1)=WK2(in1)
        IW1(icL,2)=WK2(in2)
        IW1(icL,3)=WK2(in3)
        IW1(icL,4)=WK2(in4)
        IW1(icL,5)=WK2(in5)
        IW1(icL,6)=WK2(in6)
        IW1(icL,7)=WK2(in7)
        IW1(icL,8)=WK2(in8)
      endif
      NN=NODstack(icpu)-NODstack(icpu-1)
      if(WK2(in1).gt.NN.or.
     &   WK2(in2).gt.NN.or.
     &   WK2(in3).gt.NN.or.
     &   WK2(in4).gt.NN.or.
     &   WK2(in5).gt.NN.or.
     &   WK2(in6).gt.NN.or.
     &   WK2(in7).gt.NN.or.
     &   WK2(in8).gt.NN) then
      print*,'  ####  ERR-3: Local VERTEX NO, Greater than NCV at icpu'
        print*,'     ',
     &  icpu,NN,': ',WK2(in1),WK2(in2),WK2(in3),WK2(in4),
     &               WK2(in5),WK2(in6),WK2(in7),WK2(in8)
        stop
      endif
 120  continue
      write(11,'(6i12)') nn8, NNX
      write(11,'(6i12)') ((IW1(i,k),k=1,8),i=1,NNX)

!
      deallocate(IW1)
!-----------------------
! --- Material Type
!-----------------------
      HEADER2='MATERIAL_TYPE'
      write(11,'(a)') TRIM(adjustL(HEADER2))
!
      NNX=ELMstack(icpu)-ELMstack(icpu-1)
      allocate(WK1(NNX),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -4- in OUTPUT'
      WK1=0
      icL=0
      do 130 icLG=ELMstack(icpu-1)+1,ELMstack(icpu)
      icL=icL+1
      ic=ELMitem(icLG)
      WK1(icL)=lacell(ic)
 130  continue
!
      write(11,'(6i12)') nn1,NNX
      write(11,'(6i12)') (WK1(icL),icL=1,NNX)
!
      deallocate (WK1)
!
      write(11,'(a)') '#ENDFILE'
!
      close (11)
!-------------------------------------------------------------
! --- Calculating local :
! --- 1)nface_hpc(icpu),2)ncelb_hpc(icpu),3)nedge_hpc(icpu)
! --- 4)NCV_hpc(icpu),5)NALLCV_hpc(icpu),6)NCVFAC_hpc(icpu)
! --- 7)nssfbc_hpc(icpu);8)IEMAX_hpc(icpu)
!-------------------------------------------------------------
      allocate(WK3(0:ncelb),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -5- in OUTPUT'
      allocate(WK4(0:nface),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -6- in OUTPUT'
      allocate(WK5(0:nedge),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -7- in OUTPUT'
!
      WK3=0
      WK4=0
      WK5=0
      nface_hpc(icpu)=0
      ncelb_hpc(icpu)=0
      NCV_hpc(icpu)=0
      nedge_hpc(icpu)=0
!
      do 140 icLG=ELMstack(icpu-1)+1,ELMstack(icpu)
      ic=ELMitem(icLG)
      do 142 i=1,lfcell(7,ic)
      IS=lfcell(i,ic)
      if(IS>nface) then
        write(ifle,'(1x,a,2I12)') 'ERR: IS>nface',IS,nface
        stop 'Stop at pre_part_output'
      endif
      WK4(IS)=1
 142  continue
 140  continue
!
      do 144 IS=1,nface
      if(WK4(IS).eq.1) then
        nface_hpc(icpu)=nface_hpc(icpu)+1
        do 146 i=1,2
        ic=lcface(i,IS)
        WK3(ic)=1
 146    continue
        do 147 i=1,LEFACE(5,IS)
        ie=LEFACE(i,IS)
        WK5(ie)=1
 147    continue
      endif
 144  continue
!
      do 148 ic=1,ncelb
      if(WK3(ic).eq.1) then
        ncelb_hpc(icpu)=ncelb_hpc(icpu)+1
      endif
 148  continue
!
      do 150 ie=1,nedge
      if(WK5(ie).eq.1) then
        nedge_hpc(icpu)=nedge_hpc(icpu)+1
      endif
 150  continue
!---------------------------------------------------------------
! --- Calculating local NCV_hpc(icpu): not include wall vertex
!---------------------------------------------------------------
      NCV_hpc(icpu)=NODstack(icpu)-NODstack(icpu-1)
!
! --- 
!
      deallocate (WK3)
      deallocate (WK4)
      deallocate (WK5)
!
 100  continue
!///////////////////////////////////////////////////////////////////
!
      deallocate (WK2)
!
! --- +---------------+
! --- | WRITE BC file |
! --- +---------------+
!
      allocate (WK1(INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -8- in OUTPUT'
      allocate (WK2(0:INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -9- in OUTPUT'
      allocate (WK3(0:INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -10- in OUTPUT'
      allocate (WK4(INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -11- in OUTPUT'
      allocate (WK5(0:INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -12- in OUTPUT'
      allocate (WK6(0:INODTOT),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -13- in OUTPUT'
      allocate (IW1(0:INODTOT,2),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -14- in OUTPUT'
      allocate (IW2(2,mface),stat=ierr1)
      if(ierr1.ne.0) stop 'stop at allocation -14- in OUTPUT'
!
      write(ifll,*)'    ###    WRITE AND CHECK BC LIST -------------'
!---------------------
! --- FOR all CPU 
!---------------------
      do 400 icpu=1,PETOT
!
      write(ifll,*)
      write(ifll,*)
     &'    ###    ||||||||    ICPU = ',ICPU,'    ||||||||'
!
      HEADER2='Boundary_Table'
!
      open (11,file=BCout(icpu),status='unknown')
!
      write(11,'(a)')    TRIM(adjustL(gdformat))
      write(11,'(a)')    HEADER2
      write(11,'(6i12)') NBOUND
!
      WK1=0
      ivL=0
      do 220 j=NODstack(icpu-1)+1,NODstack(icpu)
      ivL=ivL+1
      iv=NODitem(j)
      WK1(iv)=ivL
 220  continue
!------------------------
! --- All Pairs BC face
!------------------------
      if (NSUF_NODE_P.gt.0.or.NSUF_NODE_O.gt.0) then
        WK2=0 
        WK2(0)=1
        WK3=0
        WK3(0)=1
        NSUF_LOCAL=0
!....................................
! --- < 1. Cycle and Sliding BC >
!....................................
        do 300 isufp=1,NSUF_NODE_P
        iv11=ISUF_NODE_P(1,isufp,1)
        iv12=ISUF_NODE_P(1,isufp,2)
        iv13=ISUF_NODE_P(1,isufp,3)
        iv14=ISUF_NODE_P(1,isufp,4)
!
        iv21=ISUF_NODE_P(2,isufp,1)
        iv22=ISUF_NODE_P(2,isufp,2)
        iv23=ISUF_NODE_P(2,isufp,3)
        iv24=ISUF_NODE_P(2,isufp,4)
!
        if(iv14.ne.0.and.iv24.ne.0) then
!
          icpu11=NODE_ID(iv11,2)
          icpu12=NODE_ID(iv12,2)
          icpu13=NODE_ID(iv13,2)
          icpu14=NODE_ID(iv14,2)
          icpu21=NODE_ID(iv21,2)
          icpu22=NODE_ID(iv22,2)
          icpu23=NODE_ID(iv23,2)
          icpu24=NODE_ID(iv24,2)
!
!          if( WK1(iv11)/=0.and.    !????zhang
!     &        WK1(iv12)/=0.and.   
!     &        WK1(iv13)/=0.and.   
!     &        WK1(iv14)/=0.and.   
!     &        WK1(iv21)/=0.and.   
!     &        WK1(iv22)/=0.and.   
!     &        WK1(iv23)/=0.and.   
!     &        WK1(iv24)/=0) then
           if(icpu11.eq.icpu.or.
     &        icpu12.eq.icpu.or.
     &        icpu13.eq.icpu.or.
     &        icpu14.eq.icpu.or.
     &        icpu21.eq.icpu.or.
     &        icpu22.eq.icpu.or.
     &        icpu23.eq.icpu.or.
     &        icpu24.eq.icpu) then
            NSUF_LOCAL=NSUF_LOCAL+1
            WK2(iv11)=1
            WK2(iv12)=1
            WK2(iv13)=1
            WK2(iv14)=1
            WK3(iv21)=1
            WK3(iv22)=1
            WK3(iv23)=1
            WK3(iv24)=1
          endif
        elseif(iv14.eq.0.and.iv24.eq.0) then
          icpu11=NODE_ID(iv11,2)
          icpu12=NODE_ID(iv12,2)
          icpu13=NODE_ID(iv13,2)
          icpu21=NODE_ID(iv21,2)
          icpu22=NODE_ID(iv22,2)
          icpu23=NODE_ID(iv23,2)
!          if( WK1(iv11)/=0.and.
!     &        WK1(iv12)/=0.and.   
!     &        WK1(iv13)/=0.and.   
!     &        WK1(iv21)/=0.and.   
!     &        WK1(iv22)/=0.and.   
!     &        WK1(iv23)/=0
!     &   ) then
           if(icpu11.eq.icpu.or.
     &        icpu12.eq.icpu.or.
     &        icpu13.eq.icpu.or.
     &        icpu21.eq.icpu.or.
     &        icpu22.eq.icpu.or.
     &        icpu23.eq.icpu) then
            NSUF_LOCAL=NSUF_LOCAL+1
            WK2(iv11)=1
            WK2(iv12)=1
            WK2(iv13)=1
            WK3(iv21)=1
            WK3(iv22)=1
            WK3(iv23)=1
          endif
        else
          write(ifle,*)
     & '    ###    ERR-4 : PERIDICAL BC FACE NOT SAME TYPE '
          stop
        endif
 300    continue
!....................................
! --- < 2. Touch-inlet BC >
!....................................
        do 350 isufo=1,NSUF_NODE_O
        iv11=ISUF_NODE_O(1,isufo,1)
        iv12=ISUF_NODE_O(1,isufo,2)
        iv13=ISUF_NODE_O(1,isufo,3)
        iv14=ISUF_NODE_O(1,isufo,4)
        iv21=ISUF_NODE_O(2,isufo,1)
        iv22=ISUF_NODE_O(2,isufo,2)
        iv23=ISUF_NODE_O(2,isufo,3)
        iv24=ISUF_NODE_O(2,isufo,4)
        iv31=ISUF_NODE_O(3,isufo,1)
        iv32=ISUF_NODE_O(3,isufo,2)
        iv33=ISUF_NODE_O(3,isufo,3)
        iv34=ISUF_NODE_O(3,isufo,4)
        if(iv14.ne.0) then
          icpu11=NODE_ID(iv11,2)
          icpu12=NODE_ID(iv12,2)
          icpu13=NODE_ID(iv13,2)
          icpu14=NODE_ID(iv14,2)
          if(icpu11.eq.icpu.or.
     &       icpu12.eq.icpu.or.
     &       icpu13.eq.icpu.or.
     &       icpu14.eq.icpu) then
            WK2(iv21)=1
            WK2(iv22)=1
            WK2(iv23)=1
            WK2(iv24)=1
            WK3(iv31)=1
            WK3(iv32)=1
            WK3(iv33)=1
            WK3(iv34)=1
          endif
        elseif(iv14.eq.0.and.iv24.eq.0) then
          icpu11=NODE_ID(iv11,2)
          icpu12=NODE_ID(iv12,2)
          icpu13=NODE_ID(iv13,2)
          if(icpu11.eq.icpu.or.
     &       icpu12.eq.icpu.or.
     &       icpu13.eq.icpu) then
            WK2(iv21)=1
            WK2(iv22)=1
            WK2(iv23)=1
            WK3(iv31)=1
            WK3(iv32)=1
            WK3(iv33)=1
          endif
        endif
 350    continue
!
        write(ifll,*)
     &'    ###    | PERIODIC/INTERFACE BC |'
        write(ifll,*)
     &'    ###    TOTAL PERIODIC/INTERFACE FACE PAIRS : ', 
     &     NSUF_LOCAL
!
        nbufa=0
        nbufb=0
        do 1200 IS=1,nface
        if(lcface(2,IS).le.ncell) goto 1200
        do 1201 i=1,4
        ILV=lvface(i,IS)
        if(WK2(ILV).ne.1) goto 1210
 1201   continue
        nbufa=nbufa+1
        goto 1200
 1210   continue
!
        do 1211 i=1,4
        ILV=lvface(i,IS)
        if(WK3(ILV).ne.1) goto 1220
 1211   continue
        nbufb=nbufb+1
        
 1220   continue
!
 1200   continue
!
!        if(nbufb.ne.nbufa) then
        if(nbufb.ne.nbufa.and.min(nbufa,nbufb).lt.NSUF_LOCAL) then
          write(ifle,'( 
     &   "   ###    ERR-5 : FACE-A NOT EQUAL FACE-B ",2i12)')
     &    NBUFA,NBUFB
!          stop  !?????zhang
        elseif(nbufa==0.and.nbufb/=0) then
          write(ifle,'("ERR: PERIODIC or Interface no.=0")')
          stop
        else
          write(ifll,
     &'("     ###    TOTAL PERIODIC/INTERFACE FACE-A VS. FACE-B : ",
     &  2I12)')
     &     nbufa,nbufb
        endif
!------------------------------------------
! --- check cyclic and interface BC:
!------------------------------------------
        ibcy=0
        do 260 ib=1,BC_CYCL_tot
        iv1=BC_CYCL(ib,1)
        iv2=BC_CYCL(ib,2)
        if(WK2(iv1).eq.1.and.WK3(iv2).eq.1) then
          if(WK1(iv1).eq.0.or.WK1(iv2).eq.0) then
            write(ifle,*) 
     &    '    ***  ERR-6-1: Periodical pair Vertex not in same CPU'
            write(ifle,*) WK1(iv1),'',WK1(iv2),'',ibcy,iv1,iv2
            stop
          else
            ibcy=ibcy+1
          endif
        endif
 260    continue
!
        write (ifll,*) 
     &   '    ###    TOTAL CYCLIC BC  VERTEX NO.= ',IBCY
!
        if(BC_CYCL_tot/=0.and.ibcy==0) then
          write(*,*) 'WRN: Not found CYCLIC BC VERTEX'
        endif
!
        ibcy=0
        do ib=1,BC_INTR_tot
        iv1=BC_INTR(ib,1)
        iv2=BC_INTR(ib,2)
!
        if((WK2(iv1)==1.and.WK3(iv2)==1).or.
     &     (WK2(iv2)==1.and.WK3(iv1)==1)) then
          if(WK1(iv1).eq.0.or.WK1(iv2).eq.0) then
            write(ifle,*) 
     &    '    ***  ERR-6-2: Interface pair Vertex not in same CPU'
            write(ifle,*) WK1(iv1),'',WK1(iv2),'',ibcy,iv1,iv2
            stop
          else
            ibcy=ibcy+1
          endif
        endif
        enddo
        write(ifll,*) 
     &  '    ###    TOTAL INTERFACE/Buffle BC VERTEX NO.= ',IBCY
        if(BC_INTR_tot/=0.and.ibcy==0) then
!          stop 'ERR: Not found INTERFACE BC VERTEX' !zhang
        endif
      endif
!
!----------------------------------------------------------------------
! --- All BC ----------------------------------------------------------
!----------------------------------------------------------------------
! 
      N0=0
      N1=0
!
      do 500 nb=1,nbcnd
      nob=nobcnd(nb)
!---------------------------------
! --- Check CYCLIC & INTERFACE BC 
!---------------------------------
      IW2(:,1:4)=0
      kd=kdbcnd(0,nb)
      if((kd==kdprdc.and.idis(nb)==0)
     &  .or.(kd==kdbuff.and.idis(nb)==0)
     &  .or.(kd==kdpors.and.idis(nb)==0)
     &  .or.(kd==kdintr.and.idis(nb)==0)
     &  .or.(kd==kdshutr.and.idis(nb)==0)
     &   ) then
        N0=N0+2
        ip1=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)
!
        if(kd==kdprdc.and.idis(nb)==0) then
          WK5=0;WK5(0)=1
          WK6=0;WK6(0)=1
          do 450 i=ip0+1,ip0+ip1
          IV1=lvbcnd(i)
          IV2=lvbcnd(i+ip1)
          if (WK1(IV1).ne.0.and.WK2(IV1).ne.0) then
            Wk5(IV1)=1
          endif
          if (WK1(IV2).ne.0.and.WK3(IV2).ne.0) then
            Wk6(IV2)=1
          endif
 450      continue
        elseif(kd==kdintr.and.idis(nb)==0) then
          WK5=0;WK5(0)=1
          WK6=0;WK6(0)=1
          do i=ip0+1,ip0+ip1
          IV1=lvbcnd(i)
          IV2=lvbcnd(i+ip1) 
          if (WK1(IV1).ne.0.and.WK2(IV1).ne.0.or.
     &        WK1(IV1).ne.0.and.WK3(IV1).ne.0
     &      ) then
            Wk5(IV1)=1
          endif
          if (WK1(IV2).ne.0.and.WK3(IV2).ne.0.OR.
     &        WK1(IV2).ne.0.and.WK2(IV2).ne.0
     &      ) then
            Wk6(IV2)=1
          endif
          enddo
        elseif(kd==kdbuff.and.idis(nb)==0) then
          WK5=0;WK5(0)=1
          WK6=0;WK6(0)=1
          do i=ip0+1,ip0+ip1
          IV1=lvbcnd(i)
          IV2=lvbcnd(i+ip1) 
          if (WK1(IV1).ne.0.and.WK2(IV1).ne.0.or.
     &        WK1(IV1).ne.0.and.WK3(IV1).ne.0
     &      ) then
            Wk5(IV1)=1
          endif
          if (WK1(IV2).ne.0.and.WK3(IV2).ne.0.OR.
     &        WK1(IV2).ne.0.and.WK2(IV2).ne.0
     &      ) then
            Wk6(IV2)=1
          endif
          enddo
        elseif(kd==kdshutr.and.idis(nb)==0) then
          WK5=0;WK5(0)=1
          WK6=0;WK6(0)=1
          do i=ip0+1,ip0+ip1
          IV1=lvbcnd(i)
          IV2=lvbcnd(i+ip1) 
          if (WK1(IV1).ne.0.and.WK2(IV1).ne.0.or.
     &        WK1(IV1).ne.0.and.WK3(IV1).ne.0
     &      ) then
            Wk5(IV1)=1
          endif
          if (WK1(IV2).ne.0.and.WK3(IV2).ne.0.OR.
     &        WK1(IV2).ne.0.and.WK2(IV2).ne.0
     &      ) then
            Wk6(IV2)=1
          endif
          enddo
        elseif(kd==kdpors.and.idis(nb)==0) then
          WK5=0;WK5(0)=1
          WK6=0;WK6(0)=1
          do i=ip0+1,ip0+ip1
          IV1=lvbcnd(i)
          IV2=lvbcnd(i+ip1) 
          if (WK1(IV1).ne.0.and.WK2(IV1).ne.0.or.
     &        WK1(IV1).ne.0.and.WK3(IV1).ne.0
     &      ) then
            Wk5(IV1)=1
          endif
          if (WK1(IV2).ne.0.and.WK3(IV2).ne.0.OR.
     &        WK1(IV2).ne.0.and.WK2(IV2).ne.0
     &      ) then
            Wk6(IV2)=1
          endif
          enddo

        endif
!--------------------------
! --- Check BC face number
!--------------------------
        nbufa=0
        nbufb=0
        IW1(:,1)=0
        IW1(0,1)=1
        IW1(:,2)=0
        IW1(0,2)=1
        do 2200 IS=1,nface
        if(lcface(2,IS).le.ncell) cycle
        NN=1
        do 2201 i=1,4
        ILV=lvface(i,IS)
        NN=NN*WK5(ILV)
 2201   enddo
        IF(NN.ne.0) then
          IW1(lvface(1,IS),1)=1
          IW1(lvface(2,IS),1)=1
          IW1(lvface(3,IS),1)=1
          IW1(lvface(4,IS),1)=1
          nbufa=nbufa+1
          IW2(1,nbufa)=IS
        else
          goto 2210
        endif
        goto 2200
 2210   continue
        NN=1
        do 2211 i=1,4
        ILV=lvface(i,IS)
        NN=NN*WK6(ILV)
 2211   continue
        IF(NN.ne.0) THEN
          IW1(lvface(1,IS),2)=1
          IW1(lvface(2,IS),2)=1
          IW1(lvface(3,IS),2)=1
          IW1(lvface(4,IS),2)=1
          nbufb=nbufb+1
          IW2(2,nbufb)=IS
        ELSE
          GOTO 2220
        ENDIF
 2220   continue
 2200   continue
!
        IF(nbufb.ne.nbufa) then
          write(ifle,*) 
     &    'ERR-7 : FACE-A /= FACE-B : ',
     &           NBUFA,NBUFB,'AT BC NO.: ',NOB
          stop
!        elseIF(nbufa==0) then
!          write(ifle,*) 'ERR: BC NO.=',NOB
!          STOP 'ERR: Not found PERIODIC/INTERFACE BC face'
        else
          write(ifll,'(a,2I12,a,I8)')
     &     '    ###    PERIODIC/INTERFACE BC FACE-A VS FACE-B : '
     &     ,NBUFA,NBUFB,' at BC NO.: ',NOB
        endif
!----------------
! --- CYCLIC BC
!----------------
        ivv=0
        WK4=0
        do 550 i=ip0+1,ip0+ip1
        IV=lvbcnd(i)
        if(IW1(IV,1)/=0) then
!        if(WK2(IV).ne.0.and.WK1(IV)/=0) then
!        if(WK2(IV)==1) then
!        if(Wk5(IV)==1) then
          N1=N1+1
          ivL=WK1(IV)
          ivv=ivv+1
          WK4(ivv)=ivL
!          locmsh
        endif
 550    continue
        ibcy1=ivv
        mb1=boundIDMap(nb,1)
        write(11,'(a)') TRIM(adjustL(SFBOUN(mb1)))
        write(11,'(3i12)') nn1,ivv,nbufa
        write(11,'(6i12)') (WK4(i),i=1,ivv)
        if(nbufa>0) then
          do i=1,nbufa
          write(11,'(4i12)') (WK1(lvface(1:4,(IW2(1,i)))))
          enddo
        endif
!
        ivv=0
        WK4=0
        do 570 i=ip0+1,ip0+ip1
        IV=lvbcnd(i+ip1)
        if(IW1(IV,2)/=0) then
!        if(WK3(IV).ne.0.and.WK1(IV)/=0) then
!        if(WK3(IV)==1) then
!        if(Wk6(IV)==1) then
          N1=N1+1
          ivL=WK1(IV)
          ivv=ivv+1
          WK4(ivv)=ivL
        endif
 570    continue
        ibcy2=ivv
        mb2=boundIDMap(nb,2)
        write(11,'(a)')    TRIM(adjustL(SFBOUN(mb2)))
        write(11,'(3i12)') nn1,ivv,nbufb
        write(11,'(6i12)') (WK4(i),i=1,ivv)
        if(nbufb>0) then
         do i=1,nbufb
         write(11,'(4i12)') (WK1(lvface(1:4,(IW2(2,i)))))
         enddo
        endif
        write(ifll,
     & '("    ###    PERIODIC/INTERFACE BC VERTEX-A VS. VERTEX-B :",
     &      2I8)') IBCY1,IBCY2
!
        if(ibcy1/=ibcy2) then
          write(ifle,*)'  ####  ERR-8: ', 
     &    ' PERIODIC BC NUMBER NOT EQUAL AT ICPU= ',ICPU
          write(ifle,*) 'IBCY1= ',IBCY1,' IBCY2= ',IBCY2
          stop
!        elseif(ibcy1==0) then
!          stop 'PERIODIC/INTERFACE VERTEX No.=0'
        endif
      elseif(kd==kdsld.and.idis(nb)==0) then
        N0=N0+2
        ivv=0
        WK4=0
        ip1=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)
        do i=ip0+1,ip0+ip1
        IV=lvbcnd(i)
        if(WK1(IV).ne.0) then
          N1=N1+1
          ivL=WK1(IV)
          ivv=ivv+1
          WK4(ivv)=ivL
        endif
        enddo
        ibcy1=ivv


        mb1=boundIDMap(nb,1)
        write(11,'(a)') TRIM(adjustL(SFBOUN(mb1)))
        write(11,'(3i12)') nn1,ivv,0
        write(11,'(6i12)') (WK4(i),i=1,ivv)
!
        ivv=0
        WK4=0
        do i=ip0+1,ip0+ip1
        IV=lvbcnd(i+ip1)
        if(WK1(IV).ne.0) then
          N1=N1+1
          ivL=WK1(IV)
          ivv=ivv+1
          WK4(ivv)=ivL
        endif
        enddo
        ibcy2=ivv
        mb2=boundIDMap(nb,2)
        write(11,'(a)')    TRIM(adjustL(SFBOUN(mb2)))
        write(11,'(3i12)')  nn1,ivv,0
        write(11,'(6i12)') (WK4(i),i=1,ivv)
        if(ibcy2/=ibcy1) stop ' BC sliding in hpc_out'
      elseif(kd==kdtchi.and.idis(nb)==0) then
!---------------------
! --- touch inlet BC:
!---------------------
        N0=N0+1
        ips=ivbcnd(nb-1)+1
        ipe=ivbcnd(nb)
        ivv=0
        WK4=0
        do 575 i=ips,ipe
        IV=lvbcnd(i)
        if(WK1(IV).ne.0) then
          N1=N1+1
          ivv=ivv+1
          ivL=WK1(iv)
          WK4(ivv)=ivL
        endif
 575    continue
        mb=boundIDMap(nb,1)
        write(11,'(a)') TRIM(adjustL(SFBOUN(mb)))
        write(11,'(3i12)') nn1,ivv,0
        write(11,'(6i12)') (WK4(i),i=1,ivv)
      else
!---------------------
! --- Other BC 
!---------------------
        if((kd==kdprdc.and.idis(nb)==1).or.
     &     (kd==kdintr.and.idis(nb)==1).OR.
     &     (kd==kdBUFF.and.idis(nb)==1).OR.
     &     (kd==kdpors.and.idis(nb)==1).OR.
     &     (kd==kdshutr.and.idis(nb)==1)
     &    ) then
          N0=N0+2
          ivv=0
          WK4=0
          ip0=ivbcnd(nb-1)+1
          ip1=ivpair(nb)
          do i=ip0,ip1
          IV=lvbcnd(i)
          if(WK1(IV).ne.0) then
            N1=N1+1
            ivL=WK1(IV)
            ivv=ivv+1
            WK4(ivv)=ivL
          endif
          enddo
          mb1=boundIDMap(nb,1)
          write(11,'(a)') TRIM(adjustL(SFBOUN(mb1)))
          write(11,'(3i12)') nn1,ivv,0
          write(11,'(6i12)') (WK4(i),i=1,ivv)
!
          ivv=0
          WK4=0
          ip0=ivpair(nb)+1
          ip1=ivbcnd(nb)
          do i=ip0,ip1
          IV=lvbcnd(i)
          if(WK1(IV).ne.0) then
            N1=N1+1
            ivL=WK1(IV)
            ivv=ivv+1
            WK4(ivv)=ivL
          endif
          enddo
          mb2=boundIDMap(nb,2)
          write(11,'(a)')    TRIM(adjustL(SFBOUN(mb2)))
          write(11,'(3i12)')  nn1,ivv,0
          write(11,'(6i12)') (WK4(i),i=1,ivv) 
!
        elseif (kd==kdsld.and.idis(nb)==1) then
          N0=N0+2
          ivv=0
          WK4=0
          ip0=ivbcnd(nb-1)+1
          ip1=ivpair(nb)
          do i=ip0,ip1
          IV=lvbcnd(i)
          if(WK1(IV)/=0) then 
            N1=N1+1
            ivL=WK1(IV)
            ivv=ivv+1
            WK4(ivv)=ivL 
          endif
          enddo
          if((ivpair(nb)-ivbcnd(nb-1))==ivv) then
          else
            ivv=0
          endif
          mb1=boundIDMap(nb,1) 
          write(11,'(a)') TRIM(adjustL(SFBOUN(mb1)))
          write(11,'(3i12)') nn1,ivv,0
          write(11,'(6i12)') (WK4(i),i=1,ivv)
          
!
          ivv=0
          WK4=0
          ip0=ivpair(nb)+1
          ip1=ivbcnd(nb)
          do i=ip0,ip1
          IV=lvbcnd(i)
          if(WK1(IV)/=0) then
            N1=N1+1
            ivL=WK1(IV)
            ivv=ivv+1
            WK4(ivv)=ivL
          endif
          enddo
          if((ivbcnd(nb)-ivpair(nb))==ivv) then
          else
            ivv=0
          endif

          mb2=boundIDMap(nb,2)
          write(11,'(a)')    TRIM(adjustL(SFBOUN(mb2)))
          write(11,'(3i12)')  nn1,ivv,0
          write(11,'(6i12)') (WK4(i),i=1,ivv) 
!
          
        else
          N0=N0+1
          ips=ivbcnd(nb-1)+1
          ipe=ivbcnd(nb)
          ivv=0
          WK4=0
          do 560 i=ips,ipe
          IV=lvbcnd(i)
          if(WK1(IV).ne.0) then
            N1=N1+1
            ivv=ivv+1
            ivL=WK1(iv)
            WK4(ivv)=ivL
          endif
 560      continue
          mb=boundIDMap(nb,1)
          write(11,'(a)') TRIM(adjustL(SFBOUN(mb)))
          write(11,'(3i12)') nn1,ivv,0
          write(11,'(6i12)') (WK4(i),i=1,ivv)

        endif

      endif
!
 500  continue
!
      IVBCTOT_HPC(icpu)=N1
      write (11,'(a)') '#ENDFILE'
      close (11)
!
      IF(NBOUND.ne.N0.and.undef_bcno.eq.0) then
        write(ifle,*) '  #### ERR-9: Boundary number NOT matched'
        write(ifle,*) '  Read Serial BC No = ',NBOUND
        write(ifle,*) '  Write HPC BC No = ',N0
        stop
!      elseif(NBOUND.ne.N0-1.and.undef_bcno.ne.0) then
!        write(ifle,*) '  #### ERR-10: Boundary number NOT matched'
!        write(ifle,*) '  Read Serial BC No = ',NBOUND
!        write(ifle,*) '  Write HPC BC No = ',N0
!        stop
      endif
!
 400  continue

      deallocate (IW1,IW2,WK6,WK5,WK4,WK3,WK2,WK1)
      deallocate (SFBOUN,boundIDMap)
!
      return
!
      end subroutine OUTPUT_hpc_grid
