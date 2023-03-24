!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine DATA
     & (ncell,nvrtx,nface,nedge,nssfbc,NCV,NCVFAC,NALLCV,
     &  mvrtx,mcell,medge,mssfbc,
     &  cord,lvcell,lacell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_boundary,only : nbcnd,ivbcnd,lvbcnd,kdbcnd,
     &                           icbinr,lcbinr,
     &                           kdprdc,kxnone,kdsymm,kdsld,
     &                           kdintr,kdilet,kdolet,kdtchi,kvnslp,
     &                           kvlglw,idis,kdbuff,kdshutr,kdpors
      use module_io,only : ifle
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: ncell,nvrtx,nface,nedge,nssfbc
      integer,intent(in) :: mcell,mvrtx,medge,mssfbc
      integer,intent(in) :: NCV,NCVFAC,NALLCV
      real*8 ,intent(in) :: cord  (3,mvrtx)
      integer,intent(in) :: lvcell(8,mcell)
      integer,intent(in) :: lacell(  mcell)
!
! --- [local entities]
!
      integer,parameter :: nv2typ(8)=(/0,0,0,1,2,3,0,4/)
      integer           :: I,IC,kclv,ILV,nb,ips,ipe,ip1,ip0,IV1,IV2,ivl
      integer           :: ierr1,IV,kd,kdv
!
      ierr1=0
!
! --- 
!
      INODTOT=nvrtx
      allocate(NODE_ID(INODTOT,2),stat=ierr1 )
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error : allocating NODE_ID(:,:)'
        stop 'subroutine DATA'
      endif
!
      NODE_ID(:,:)=0
!
      ICELTOT=ncell
!
      allocate(ICELTYP(ncell),stat=ierr1 )
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error : allocating ICELTYP(:)'
        stop 'subroutine DATA'
      endif
      ICELTYP=0
!
      do 120 ic=1,ICELTOT
      kclv=0
      do 125 ILV=1,8
      if (lvcell(ILV,ic).gt.0) kclv=ILV
 125  continue
      ICELTYP(ic)=nv2typ(kclv)
 120  continue
!
!
      BC_INLT_tot=0
      BC_WALL_tot=0
      BC_FREE_tot=0
      BC_SYMT_tot=0
      BC_CYCL_tot=0
      BC_BODY_tot=0
      BC_MWAL_tot=0
      BC_TCIN_tot=0
      BC_SLID_tot=0
      BC_INTR_tot=0
!
      do 1000 nb=1,nbcnd
      ips=ivbcnd(nb-1)
      ipe=ivbcnd(nb)
      kd=kdbcnd(0,nb)
      kdv=kdbcnd(1,nb)
      if((kdv==kvnslp.or.kdv==kvlglw)
     &  .and.(kd/=kdintr).and.(kd/=kdbuff)) then
        BC_WALL_tot=BC_WALL_tot+ipe-ips
      elseif( kd==kdprdc.and.idis(nb)==0) then
        BC_CYCL_tot=BC_CYCL_tot+(ipe-ips)/2
      elseif( kd==kdsld.and.idis(nb)==0 ) then
        BC_SLID_tot=BC_SLID_tot+(ipe-ips)/2
      elseif( kd==kdsymm ) then
        BC_SYMT_tot=BC_SYMT_tot+ipe-ips
      elseif( kd==kdilet ) then
        BC_INLT_tot=BC_INLT_tot+ipe-ips
      elseif( kd==kdolet ) then
        BC_FREE_tot=BC_FREE_tot+ipe-ips
      elseif( kd==kdtchi.and.idis(nb)==0) then
        BC_TCIN_tot=BC_TCIN_tot+ipe-ips
      elseif((kd==kdintr.or.kd==kdbuff.or.kd==kdshutr.or.kd==kdpors)
     &    .and.idis(nb)==0 ) then
        BC_INTR_tot=BC_INTR_tot+ipe-ips
        
      endif
!
 1000 continue
!
      if(BC_INLT_tot.gt.0) then 
        allocate(BC_INLT(BC_INLT_tot))
        allocate(BC_IV3D(BC_INLT_tot,3))
      endif
      if(BC_WALL_tot.gt.0) allocate(BC_WALL(BC_WALL_tot))
      if(BC_SYMT_tot.gt.0) allocate(BC_SYMT(BC_SYMT_tot))
      if(BC_CYCL_tot.gt.0) allocate(BC_CYCL(BC_CYCL_tot,2))
      if(BC_INTR_tot.gt.0) allocate(BC_INTR(BC_INTR_tot,2))
      if(BC_BODY_tot.gt.0) allocate(BC_BODY(BC_BODY_tot))
      if(BC_FREE_tot.gt.0) allocate(BC_FREE(BC_FREE_tot))
      if(BC_MWAL_tot.gt.0) then
        allocate(BC_MWAL(BC_MWAL_tot))
        allocate(BC_MV3D(BC_MWAL_tot,3))
      endif
      if(BC_TCIN_tot.gt.0)  allocate(BC_TCIN(BC_TCIN_tot))
!
      BC_INLT_tot=0
      BC_WALL_tot=0
      BC_FREE_tot=0
      BC_SYMT_tot=0
      BC_CYCL_tot=0
      BC_BODY_tot=0
      BC_MWAL_tot=0
      BC_TCIN_tot=0
      BC_SLID_tot=0
      BC_INTR_tot=0
!
      do 2000 nb=1,nbcnd
      ips=ivbcnd(nb-1)+1
      ipe=ivbcnd(nb)
      kd=kdbcnd(0,nb)
      kdv=kdbcnd(1,nb)
      if((kdv==kvnslp.or.kdv==kvlglw).and.
     &      (kd/=kdintr).and.(kd/=kdbuff)
     & .and.(kd/=kdshutr).and.(kd/=kdpors)) then
        do 2100 ivl=ips,ipe
        BC_WALL_tot=BC_WALL_tot+1
        IV=lvbcnd(ivl)
        BC_WALL(BC_WALL_tot)=IV
 2100   continue
      elseif(kd==kdprdc.and.idis(nb)==0) then
        ip1=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)        
        do 2300 ivl=ip0+1,ip0+ip1
        IV1=lvbcnd(ivl)
        IV2=lvbcnd(ivl+ip1)
        BC_CYCL_tot=BC_CYCL_tot+1
        BC_CYCL(BC_CYCL_tot,1)=IV1
        BC_CYCL(BC_CYCL_tot,2)=IV2
 2300   continue
      elseif(kd==kdsld.and.idis(nb)==0) then
        ip1=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)        
        do ivl=ip0+1,ip0+ip1
        enddo
      elseif( kd==kdsymm ) then
        do 2400 ivl=ips,ipe
        BC_SYMT_tot=BC_SYMT_tot+1
        IV=lvbcnd(ivl)
        BC_SYMT(BC_SYMT_tot)=IV
 2400   continue
      elseif( kd==kdilet ) then
        do 2500 ivl=ips,ipe
        BC_INLT_tot=BC_INLT_tot+1
        IV=lvbcnd(ivl)
        BC_INLT(BC_INLT_tot)=IV
        BC_IV3D(BC_INLT_tot,:)=0
 2500   continue
      elseif( kd==kdolet ) then
        do 2600 ivl=ips,ipe
        BC_FREE_tot=BC_FREE_tot+1
        IV=lvbcnd(ivl)
        BC_FREE(BC_FREE_tot)=IV
 2600   continue
      elseif( kd==kdtchi.and.idis(nb)==0 ) then
        do 2800 ivl=ips,ipe
        BC_TCIN_tot=BC_TCIN_tot+1
        IV=lvbcnd(ivl)
        BC_TCIN(BC_TCIN_tot)=IV
 2800   continue
      elseif( (kd==kdintr.or.kd==kdbuff.or.kd==kdshutr.or.kd==kdpors)
     &    .and.idis(nb)==0) then
        ip1=(ivbcnd(nb)-ivbcnd(nb-1))/2
        ip0=ivbcnd(nb-1)
        do 2700 ivl=ip0+1,ip0+ip1
        IV1=lvbcnd(ivl)
        IV2=lvbcnd(ivl+ip1)
        BC_INTR_tot=BC_INTR_tot+1
        BC_INTR(BC_INTR_tot,1)=IV1
        BC_INTR(BC_INTR_tot,2)=IV2
 2700   continue
      endif
!
 2000 continue
!
      end subroutine DATA






