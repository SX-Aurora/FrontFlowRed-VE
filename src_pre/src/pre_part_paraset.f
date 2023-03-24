!
!     subroutine PARASET
!     subroutine DEFINE_FILE_NAME
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine PARASET(medge,nedge,mssfbc,IEMAX,
     &           LBCSSF,LCYCSF,LVEDGE)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_partitioner
      use module_io,only : ifle,ifll
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: medge,nedge,mssfbc,IEMAX
      integer,intent(in) :: LVEDGE(2,medge)
      integer,intent(in) :: LBCSSF    ( mssfbc)
      integer,intent(in) :: LCYCSF    ( mssfbc)
!
! --- [local entities]
!
      integer :: ierr1=0,ii,i
      integer :: LENGTH,my_rank,NPARTMAX
      character(len=80) :: PARM_FILE
!
!
! --- METIS input (Creat METIS Input file => STOP)
!
      call CREATE_METIS_INPUT(medge,nedge,mssfbc,IEMAX,
     &           LBCSSF,LCYCSF,LVEDGE)
!
! --- Partitioning Graph-Mesh by METIS
!
      call metis(INPFIL,NPART)
!
! --- 
!
      IF(NPART.GT.1.AND.NPART.LT.10) THEN
        write(METISFIL,'(i1)') NPART
      ELSEIF(NPART.GE.10.AND.NPART.LT.100) THEN

        write(METISFIL,'(i2)') NPART
      ELSEIF(NPART.GE.100.AND.NPART.LT.1000) THEN
        write(METISFIL,'(i3)') NPART
      ELSE
        WRITE(ifll,*) '*** The number of partitions greater than 1000'
        WRITE(ifll,*) 'Please check subroutine PARASET'
        STOP
      ENDIF
!
      INPFIL=adjustL(INPFIL)
      LENGTH=len_trim(INPFIL)
      METISFIL=adjustl(METISFIL)
      METISFIL=INPFIL(1:LENGTH)//'.part'//'.'//TRIM(METISFIL)
      open (85,file=METISFIL,status='unknown')
      NPARTMAX=-100
      do 100 i=1,INODTOT
      read(85,*,err=998,end=999) ii
      NODE_ID(i,2)=ii+1
      NPARTMAX=max(NPARTMAX,ii+1)
 100  continue
      close(85)
! --- NPARTMAX='CPU number'
      PETOT=NPARTMAX
!
      IF(PETOT.NE.NPART) then
        write (ifll,*) 'ERROR :'
        write (ifll,*) 
     & 'number of partitions is NOT equal to CPU No.',NPART, PETOT
        stop
      ENDIF
!
      if (PETOT.lt.1) then
        call ERROR_EXIT(32,0)
      elseif(PETOT.gt.INODTOT) then
        call ERROR_EXIT(32,INODTOT)
      endif
!
! --- Defining file name:
!
      allocate (GRIDout(PETOT), BCout(PETOT), COMMout(PETOT), 
     &          WORKFIL(PETOT),hpcdmn(PETOT),GEOMout(PETOT),
     &          RODRICV(PETOT),walldis(PETOT),SLIDNG(PETOT),
     &	        RADout(PETOT),	!jiang
     &          stat=ierr1)
!
      do 200 my_rank=0,PETOT-1
      call DEFINE_FILE_NAME(DIRHED,HEADERg,my_rank, GRIDout(my_rank+1))
      call DEFINE_FILE_NAME(DIRHED,HEADERb,my_rank, BCout  (my_rank+1))
      call DEFINE_FILE_NAME(DIRHED,HEADERc,my_rank, COMMout(my_rank+1))
      call DEFINE_FILE_NAME(DIRHED,HEADERw,my_rank, WORKFIL(my_rank+1))
!jiang
      call DEFINE_FILE_NAME(DIRHED,HEADERw,my_rank, RADout(my_rank+1))
      PARM_FILE='dimen.parm'
      call 
     &  DEFINE_FILE_NAME(DIRHED,PARM_FILE,my_rank, hpcdmn(my_rank+1))
 200  continue
!
! --- 
!
      allocate(nvrtx_hpc(PETOT))
      allocate(ncell_hpc(PETOT))
      allocate(ncelb_hpc(PETOT))
      allocate(nface_hpc(PETOT))
      allocate(nedge_hpc(PETOT))
      allocate(NCV_hpc(PETOT))
      allocate(NALLCV_hpc(PETOT))
      allocate(NCVFAC_hpc(PETOT))
      allocate(nssfbc_hpc(PETOT))
      allocate(ncomp_hpc(PETOT))
      allocate(ncomp_suf_hpc(PETOT))
      allocate(nphase_hpc(PETOT))
      allocate(nrans_hpc(PETOT))
      allocate(npotn_hpc(PETOT))
      allocate(NMAT_hpc(PETOT))
      allocate(IEMAX_hpc(PETOT))
      allocate(NBCINL_hpc(PETOT))
      allocate(NBCWAL_hpc(PETOT))
      allocate(NBCCYC_hpc(PETOT))
      allocate(NBCOUT_hpc(PETOT))
      allocate(NBCSYM_hpc(PETOT))
      allocate(NBCTCI_hpc(PETOT))
      allocate(NBCSLD_hpc(PETOT))
      allocate(NBCINT_hpc(PETOT))
      allocate(IVBCTOT_HPC(PETOT))
!
      nvrtx_hpc=0
      ncell_hpc=0
      ncelb_hpc=0
      nface_hpc=0
      nedge_hpc=0
      NCV_hpc=0
      NALLCV_hpc=0
      NCVFAC_hpc=0
      nssfbc_hpc=0
      ncomp_hpc=0
      ncomp_suf_hpc=0
      nphase_hpc=0
      nrans_hpc=0
      NMAT_hpc=0
      IEMAX_hpc=0
      NBCINL_hpc=0
      NBCWAL_hpc=0
      NBCCYC_hpc=0
      NBCOUT_hpc=0
      NBCSYM_hpc=0      

!
      return
!
 997  continue
        call ERROR_EXIT (8001,0)
 998  continue
        call ERROR_EXIT (21,0)
 999  continue
        call ERROR_EXIT (22,0)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      end subroutine PARASET
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine DEFINE_FILE_NAME(DIRHED,HEADER,my_rank,filname)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      implicit none
!
! --- [dummy arguments]
!
      character(80),intent(inout) :: filname,HEADER,DIRHED
      integer,     intent(in)    :: my_rank
!
! --- [local entities]
!
      character (len= 4)         :: SUBindex
      integer                    :: ID,LENGTH,LENGT
      
!
! --- 
!
      HEADER=adjustL(HEADER)
      LENGTH=len_trim(HEADER)
      DIRHED=adjustL(DIRHED)
      LENGT=len_trim(DIRHED)
      write(SUBindex,'(i4.4)') my_rank
!
      filname=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADER(1:LENGTH)
!
      return
      end subroutine DEFINE_FILE_NAME
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine CREATE_METIS_INPUT(medge,nedge,mssfbc,IEMAX,
     &           LBCSSF,LCYCSF,LVEDGE)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_partitioner
      use module_io      ,only : ifle,ifll
      use module_boundary,only : nbcnd,LBC_INDEX,kdbcnd,conntBC,
     &                           idis,kdintr,kdbuff,kdshutr,kdpors
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in) :: medge,nedge,IEMAX,mssfbc
      integer,intent(in) :: LVEDGE(2,medge)
      integer,intent(in) :: LBCSSF    ( mssfbc)
      integer,intent(in) :: LCYCSF    ( mssfbc)
!
! --- 
!
      integer  :: i,j,k,in1,in2,ik1,ik2,ie,iv,NEIBMAXmetis
      integer  :: ierr1=0
      character(len=20)  :: iedge
      integer  :: ICFL,ICFP,IBFL,IBFS,IBFE,NB,kd
      integer,save  :: IEDGTOT_M
!
! --- WK1(inod)=WK1(innod)+1
!
      allocate (WK1(INODTOT),stat=ierr1)
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error : allocating WK1(INODTOT)'
        stop 'subroutine PARASET'
      endif
!
      WK1(:)=0
      IEDGTOT=nedge
      IEDGTOT_M=nedge
      do 100 ie=1,IEDGTOT
      in1=LVEDGE(1,ie)
      in2=LVEDGE(2,ie)
      WK1(in1)=WK1(in1)+1
      WK1(in2)=WK1(in2)+1
 100  continue
!
      DO nb=1,nbcnd
      kd=kdbcnd(0,NB)
      if(conntBC(nb)) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        DO IBFL=IBFS,IBFE
        IEDGTOT_M=IEDGTOT_M+1
        ICFL=LBCSSF(IBFL)
        ICFP=LCYCSF(IBFL)
        in2=LVEDGE(1,ICFP)
        in1=LVEDGE(1,ICFL)
        WK1(in2)=WK1(in2)+1
        WK1(in1)=WK1(in1)+1
        enddo
      endif
      enddo
!
! --- 'NEIBMAXmetis' is max edge number
!
      NEIBMAXmetis=-100
      do 150 iv=1,INODTOT
      NEIBMAXmetis=max(NEIBMAXmetis,WK1(iv))
 150  continue
!
! --- IW1(inod1,NEIBMAXmetis)=inod2
!
      allocate(IW1(INODTOT,NEIBMAXmetis),stat=ierr1)
      if( ierr1.ne.0 ) then
        write(ifle,*) '### error : allocating IW1(INODTOT,NEIBMAXmetis)'
        stop 'subroutine PARASET'
      endif
      WK1=0
      IW1=0
!
! --- neighboring NODEs : IW1(inod1,WK1(inod1))=inod2
!
      do 250 ie=1,IEDGTOT
      in1=LVEDGE(1,ie)
      in2=LVEDGE(2,ie)
      WK1(in1)=WK1(in1)+1
      WK1(in2)=WK1(in2)+1
      IW1(in1,WK1(in1))=in2
      IW1(in2,WK1(in2))=in1
 250  continue
!
!
!
      DO nb=1,nbcnd
      kd=kdbcnd(0,NB)
      if(conntBC(nb)) then
        IBFS=LBC_INDEX(nb-1)+1
        IBFE=LBC_INDEX(nb)
        DO IBFL=IBFS,IBFE
        ICFL=LBCSSF(IBFL)
        ICFP=LCYCSF(IBFL)
        in2=LVEDGE(1,ICFP)
        in1=LVEDGE(1,ICFL)
        WK1(in2)=WK1(in2)+1
        WK1(in1)=WK1(in1)+1
        IW1(in1,WK1(in1))=in2
        IW1(in2,WK1(in2))=in1
        enddo
      endif
      enddo      
!
! --- INPFIL=test.m
!

      write(iedge,*) "(",IEMAX,"i12",")"
      open  (22,file=INPFIL,status='unknown')
      write (22,'(2i12)') INODTOT, IEDGTOT_M
!
      do 500 iv=1,INODTOT
      write (22,iedge) (IW1(iv,k),k=1,WK1(iv))
 500  continue
      close (22)
!
      deallocate (WK1)
      deallocate (IW1)
!
      return
!
      end subroutine CREATE_METIS_INPUT

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
