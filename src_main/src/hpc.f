!     The subroutines included in this file is used for real HPC
!
!      subroutine HPC_INIT
!      subroutine FFRABORT
!      subroutine FFTIME
!      subroutine SOLVER_SEND_RECV
!      subroutine hpcimax
!      subroutine hpcimax_0
!      subroutine hpcrmax
!      subroutine hpcrmax_0
!      subroutine hpcimin
!      subroutine hpcimin_0
!      subroutine hpcrmin
!      subroutine hpcrmin_0
!      subroutine hpcisum
!      subroutine hpcisum_0
!      subroutine hpcrsum
!      subroutine hpcrsum_0
!      subroutine hpcibcast
!      subroutine hpcrbcast
!      subroutine hpcland
!      subroutine hpclor
!      subroutine hpcmaxloc
!      subroutine read_griddata_hpc
!      subroutine read_hpc_GF
!      subroutine read_hpc
!      subroutine comm
!      SUBROUTINE MERGEOUT_PTRESULT_MPI
!      subroutine hpcrbcast_stream
!      subroutine hpc_igather
!      subroutine hpc_dgather
!      subroutine hpcrbcast_ary
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine HPC_INIT(NCPU,my_cpu)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     initializing MPI
!
      use module_hpcutil
!
      implicit none
      integer,intent(out)  :: NCPU,my_cpu
      integer  :: ierr,ITERHPC
      real*8   :: RESID
!
      NCPU=1
      my_cpu=0
!
      NPE=1
      PETOT=1
      ROOT=0
      my_rank=0
     
      RESID=1.D-6
      ITERHPC=1000
!
      call MPI_INIT      (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr )
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr )

      call MPI_COMM_DUP  (MPI_COMM_WORLD, SOLVER_COMM, ierr)
!
      NPE=PETOT
      NCPU=NPE
      my_cpu=my_rank
!
      return
!
      end subroutine HPC_INIT
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine FFRABORT(ICODE,subname)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
      use module_io,    only : ifll,ifle
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout)   ::  ICODE
      character(*),intent(in) ::  subname
!
! --- [local entities]
!
      integer  :: ierr
!      integer  :: COMM
!
!
      IF(ICODE.EQ.1) THEN 
        WRITE(ifle,'(1X,2a)') 'STOP Subroutine name: ',
     &                 subname(:len_trim(subname))
        WRITE(ifle,'(1X,a,I8)') 'STOP at MY_RANK= ',MY_RANK
!
        call MPI_ABORT(SOLVER_COMM,ierr)
      ELSEIF(ICODE.EQ.0) THEN 
        if(my_rank==root) write(ifll,3000)
        call MPI_FINALIZE (ierr)
        stop
      ELSE
        WRITE(ifle,*) 'STOP ERROR CODE NO.= ',ICODE,
     &                  ' MY_RANK= ',MY_RANK
        call MPI_ABORT(SOLVER_COMM,ierr)
      ENDIF
 3000 FORMAT(2X,'|',108('#'),'|')
!
      end subroutine FFRABORT
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine FFTIME(ITIME)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
      use module_io,    only : ifll,ifle
!
      implicit none
! --- [dummy arguments]
!
      integer,intent(in)  ::  ITIME
!
! --- [local entities]
!
      real*8,save   :: mpi_start=0.d0,mpi_end
      integer  :: NDT,NHOUR,NMIN,NSECEND
      integer  :: ierr
      INTEGER  :: val(8)
!
      CALL DATE_AND_TIME(VALUES=val)
      IF(ITIME.eq.0) then
        mpi_start=MPI_WTIME()
      else if(ITIME.eq.999) then
         mpi_end=MPI_WTIME()
         if(my_rank.eq.ROOT) then 
           NDT=INT(mpi_end-mpi_start)
           NHOUR=NDT/3600
           NMIN=MOD(NDT,3600)/60
           NSECEND=MOD((MOD(NDT,3600)),60)
           write(ifll,3020) NHOUR,NMIN,NSECEND,
     &              val(1),val(2),val(3),val(5),val(6),val(7)
         endif
      else
         mpi_end=MPI_WTIME()
         if(my_rank.eq.ROOT) then 
           NDT=INT(mpi_end-mpi_start)
           NHOUR=NDT/3600
           NMIN=MOD(NDT,3600)/60
           NSECEND=MOD((MOD(NDT,3600)),60)
           write(ifll,3010) NHOUR,NMIN,NSECEND,
     &              val(1),val(2),val(3),val(5),val(6),val(7)
         endif
      endif
 3010 FORMAT(2X,'|',10X,' Used time : ',I4.4,':',I2.2,':',I2.2,
     & 2X,' | Date in Year:',I4.4,2X,'Mon :',I2.2,2X,'Day :',I2.2, 
     & 2X,'Hour:',I2.2,2X,'Min :',I2.2,2X,'Sec :',I2.2,8X,'|')
 3020 FORMAT(2X,'|',10X,' Wall time : ',I4.4,':',I2.2,':',I2.2,
     & 2X,' | Stop in Year:',I4.4,2X,'Mon :',I2.2,2X,'Day :',I2.2, 
     & 2X,'Hour:',I2.2,2X,'Min :',I2.2,2X,'Sec :',I2.2,8X,'|')
      return 
      end subroutine FFTIME 
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine SOLVER_SEND_RECV(ndim,MXNP,NP,X)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil
      use module_metrix,only : WS,WR,sta1,sta2,req1,req2
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: MXNP,NP,ndim
      real*8, intent(inout):: X(MXNP,ndim)
!
! --- [local entities]
!
!      integer,save, allocatable :: sta1(:,:)
!      integer,save, allocatable :: sta2(:,:)
!      integer,save, allocatable :: req1(:  )
!      integer,save, allocatable :: req2(:  )
!      real*8 ,allocatable       :: WS(:),WR(:)
      integer                   :: neib,istart,inum,k,idim
      integer,save              :: NFLAG=0
      integer  :: ierr
!      data NFLAG/0/             
!
! --- INIT.
!
! DEBUG
!      if (NFLAG.eq.0) then
!        allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
!        allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
!        allocate (req1(NEIBPETOT))
!        allocate (req2(NEIBPETOT))
!        allocate(WS(EXPORT_index(NEIBPETOT)),
!     &           WR(IMPORT_index(NEIBPETOT)))
!        NFLAG=1
!      endif
! --- 
!      allocate(WS(EXPORT_index(NEIBPETOT)),
!     &         WR(IMPORT_index(NEIBPETOT)))
      
!
      do 100 idim=1,ndim
      WS=0.d0
      WR=0.d0
!
! --- SEND
!
      do neib= 1, NEIBPETOT
        istart= EXPORT_index(neib-1)
        inum  = EXPORT_index(neib  ) - istart
        do k=istart+1,istart+inum
           WS(k)=X(EXPORT_item(k),idim)
        enddo


        call MPI_ISEND (WS(istart+1), inum, MPI_DOUBLE_PRECISION,       
     &                  NEIBPE(neib), 0, SOLVER_COMM,                   
     &                  req1(neib), ierr)


      enddo
      
!
! --- RECEIVE 
!
      do neib= 1, NEIBPETOT
        istart=IMPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart
        call MPI_IRECV (WR(istart+1), inum, MPI_DOUBLE_PRECISION,
     &                  NEIBPE(neib), 0, SOLVER_COMM,
     &                  req2(neib), ierr)
      enddo
!
! --- WAIT
!
      call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
!   
      do neib=1,NEIBPETOT
        istart=IMPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart
      do k=istart+1,istart+inum
        X(IMPORT_item(K),idim)= WR(k)
      enddo
      enddo
!
      call MPI_WAITALL (NEIBPETOT,req1,sta1,ierr)
 100  continue
!
!      deallocate(WS,WR)
!
      end subroutine solver_send_recv
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine SOLVER_SEND_RECV0(ndim,ini,MXNP,NP,X)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil
      use module_metrix,only : WS,WR,sta1,sta2,req1,req2
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: ini,MXNP,NP,ndim
      real*8, intent(inout):: X(ini:MXNP,ndim)
!
! --- [local entities]
!
      integer                   :: neib,istart,inum,k,idim
      integer,save              :: NFLAG=0
      integer  :: ierr
!
! --- INIT.
!
#ifdef SX_TUNE
      do 100 idim=1,ndim
      WS=0.d0
      WR=0.d0
!
! --- SEND
!
      do neib= 1, NEIBPETOT
        istart= EXPORT_index(neib-1)
        inum  = EXPORT_index(neib  ) - istart
        do k=istart+1,istart+inum
           WS(k)=X(EXPORT_item(k),idim)
        enddo
        call MPI_SEND (WS(istart+1), inum, MPI_DOUBLE_PRECISION,
     &                  NEIBPE(neib), 0, SOLVER_COMM,
     &                  ierr)
      enddo
!
! --- RECEIVE
!
      do neib= 1, NEIBPETOT
        istart=IMPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart
        call MPI_IRECV (WR(istart+1), inum, MPI_DOUBLE_PRECISION,
     &                  NEIBPE(neib), 0, SOLVER_COMM,
     &                  req2(neib), ierr)
      enddo
!
! --- WAIT
!
      call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
!
      do neib=1,NEIBPETOT
        istart=IMPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart
      do k=istart+1,istart+inum
        X(IMPORT_item(K),idim)= WR(k)
      enddo
      enddo
!
 100  continue
#else
      do 100 idim=1,ndim
      WS=0.d0
      WR=0.d0
!
! --- SEND
!
      do neib= 1, NEIBPETOT
        istart= EXPORT_index(neib-1)
        inum  = EXPORT_index(neib  ) - istart
        do k=istart+1,istart+inum
           WS(k)=X(EXPORT_item(k),idim)
        enddo
        call MPI_ISEND (WS(istart+1), inum, MPI_DOUBLE_PRECISION,       
     &                  NEIBPE(neib), 0, SOLVER_COMM,                   
     &                  req1(neib), ierr)
      enddo
!
! --- RECEIVE
!
      do neib= 1, NEIBPETOT
        istart=IMPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart
        call MPI_IRECV (WR(istart+1), inum, MPI_DOUBLE_PRECISION,
     &                  NEIBPE(neib), 0, SOLVER_COMM,
     &                  req2(neib), ierr)
      enddo
!
! --- WAIT
!
      call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
!   
      do neib=1,NEIBPETOT
        istart=IMPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart
      do k=istart+1,istart+inum
        X(IMPORT_item(K),idim)= WR(k)
      enddo
      enddo
!
      call MPI_WAITALL (NEIBPETOT,req1,sta1,ierr)
 100  continue
#endif /** SX_TUNE **/
!
!      deallocate(WS,WR)
!
      end subroutine solver_send_recv0

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcimax(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(inout)  :: rin
      integer                :: rout
      integer  :: ierr
!
      CALL MPI_allREDUCE (rin,rout,1,MPI_INTEGER,
     &                      MPI_MAX,SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcimax
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcimax_0(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(inout)  :: rin
      integer                :: rout
      integer  :: ierr
!
      CALL MPI_REDUCE (rin,rout,1,MPI_INTEGER,
     &                      MPI_MAX,0,SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcimax_0
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrmax(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      real*8 ,intent(inout)  :: rin
      real*8                 :: rout
      integer  :: ierr
!
      CALL MPI_allREDUCE (rin,rout, 1, MPI_DOUBLE_PRECISION,
     &                      MPI_MAX, SOLVER_COMM, ierr)
      rin=rout
      return
      end subroutine hpcrmax
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrmax_0(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      real*8 ,intent(inout)  :: rin
      real*8                 :: rout
      integer  :: ierr
!
      CALL MPI_REDUCE(rin,rout, 1, MPI_DOUBLE_PRECISION,
     &                      MPI_MAX, 0,SOLVER_COMM, ierr)
      rin=rout
      return
      end subroutine hpcrmax_0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcimin(rin) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(inout)  :: rin
      integer                :: rout
      integer  :: ierr
!
      CALL MPI_allREDUCE (rin,rout, 1, MPI_INTEGER,
     &                      MPI_MIN, SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcimin
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcimin_0(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(inout)  :: rin
      integer                :: rout
      integer  :: ierr
!
      CALL MPI_REDUCE (rin,rout, 1, MPI_INTEGER,
     &                      MPI_MIN, 0,SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcimin_0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrmin(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      real*8 ,intent(inout)  :: rin
      real*8                 :: rout
      integer  :: ierr
!
      CALL MPI_allREDUCE (rin,rout, 1, MPI_DOUBLE_PRECISION,
     &                      MPI_MIN, SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcrmin
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrmin_0(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      real*8 ,intent(inout)  :: rin
      real*8                 :: rout
      integer                :: ierr
!
      CALL MPI_REDUCE (rin,rout, 1, MPI_DOUBLE_PRECISION,
     &                      MPI_MIN, 0,SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcrmin_0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcisum(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(inout)  :: rin
      integer                :: rout
      integer                :: ierr
!
      CALL MPI_allREDUCE (rin,rout, 1, MPI_INTEGER,
     &                      MPI_SUM, SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcisum
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcisum_0(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(inout)  :: rin
      integer                :: rout
      integer                :: ierr
!
      CALL MPI_REDUCE (rin,rout, 1, MPI_INTEGER,
     &                      MPI_SUM, 0,SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcisum_0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrsum(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      real*8 ,intent(inout)  :: rin
      real*8                 :: rout
      integer                :: ierr
!
      CALL MPI_allREDUCE (rin,rout, 1, MPI_DOUBLE_PRECISION,
     &                      MPI_SUM, SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcrsum
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrsum_0(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      real*8 ,intent(inout)  :: rin
      real*8                 :: rout
      integer                :: ierr
!
      CALL MPI_REDUCE (rin,rout, 1, MPI_DOUBLE_PRECISION,
     &                      MPI_SUM, 0,SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcrsum_0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrasum(rin,idim)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(in)    :: idim
      real*8,intent(inout)  :: rin(idim)
      real*8                :: rout(idim)
      integer               :: ierr
!
      CALL MPI_allREDUCE (rin,rout,idim, MPI_DOUBLE_PRECISION,
     &                      MPI_SUM, SOLVER_COMM, ierr)
      rin=rout
!
      return
      end subroutine hpcrasum
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrasum_0(rin,idim,idst)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(in)    :: idim,idst
      real*8,intent(inout)  :: rin(idim)
!
      real*8                :: rout(idim)
      integer               :: ierr
!
      rout=0.d0
      CALL MPI_REDUCE(rin,rout,idim, MPI_DOUBLE_PRECISION,
     &                      MPI_SUM,idst,SOLVER_COMM, ierr)
      rin(:)=rout(:)
!
      return
      end subroutine hpcrasum_0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcibcast(icpu,RIN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- one-to-all
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(inout)  :: RIN
      integer,intent(in)     :: icpu
      integer                :: ierr
!
      call MPI_BCAST (RIN,1,MPI_INTEGER,icpu,SOLVER_COMM,ierr)
!
      return
      end subroutine hpcibcast
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrbcast(icpu,RIN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- one-to-all
      use module_hpcutil
!
! --- [dummy arguments]
!???
      implicit none
      integer,intent(in)     :: icpu
      real*8 ,intent(inout)  :: RIN
      integer                :: ierr
!
      call MPI_BCAST(RIN,1,MPI_DOUBLE_PRECISION,icpu,SOLVER_COMM,ierr)
!
      return
      end subroutine hpcrbcast
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrbcast_ary(icpu,RIN,NDIM)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- one-to-all
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(in)     :: icpu,NDIM
      real*8 ,intent(inout)  :: RIN(NDIM)
      integer                :: ierr
!
      call MPI_BCAST
     &  (RIN,NDIM,MPI_DOUBLE_PRECISION,icpu,SOLVER_COMM,ierr)
!
      return
      end subroutine hpcrbcast_ary
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrbcast_array(icpu,NDIM,N,WORK)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- one-to-all
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(in)     :: icpu
      integer,intent(in)     :: NDIM,N
      real*8 ,intent(inout)  :: WORK(N,NDIM)
      integer                :: ierr,i
!
      do i=1,NDIM
      call MPI_BCAST(WORK(:,i),N,MPI_DOUBLE_PRECISION,
     &  icpu,SOLVER_COMM,ierr)
      enddo
!
      return
      end subroutine hpcrbcast_array
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcland(RIN,imode)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- 'logical and' for all cpu
!
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      logical,intent(inout)   :: RIN
      integer,intent(in)      :: imode
      logical                 :: ROUT
      integer                 :: ierr
!
!      ierr=0
!
      CALL MPI_allREDUCE (rin,rout,1,MPI_LOGICAL,
     &                     MPI_LAND,SOLVER_COMM,ierr)
!     &                     MPI_BAND,SOLVER_COMM,ierr)
      rin=rout

!
      return
      end subroutine hpcland
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpclor(RIN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- 'logical or' for all cpu
!
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      logical,intent(inout)   :: RIN
      logical                 :: ROUT!=.false.
      integer                 :: ierr
!
      CALL MPI_allREDUCE (rin,rout, 1, MPI_LOGICAL,
     &                      MPI_LOR, SOLVER_COMM, ierr)
!     &                      MPI_BOR, SOLVER_COMM, ierr)
      rin=rout
!MPI_LOR
      return
      end subroutine hpclor
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcmaxloc(INDX,RMAX)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- 'MAX + Location' for all cpu
!
      use module_hpcutil
!
! --- [dummy arguments]
!
      integer,intent(inout) :: INDX
      real*8 ,intent(inout) :: RMAX
      integer               :: ierr
!
! --- [local entities]
!
      real*8  :: RIN(2),out(2)
!
      RIN(1)=RMAX
      RIN(2)=INDX
      ierr=0
      CALL MPI_allREDUCE (rin,out, 1, MPI_2DOUBLE_PRECISION,
     &                      MPI_MAXLOC,SOLVER_COMM,ierr)
      RMAX=out(1)
      INDX=out(2)
!
      return
      end subroutine hpcmaxloc


!!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      subroutine read_griddata_hpc(cord,lacell,lvcell,ierror)
!!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!! --- [module arguments]
!!
!      use module_io      ,only : ifle,ifll
!      use module_io,only       : gdformat
!      use module_dimension
!      use module_gf
!!
!      implicit none
!!
!! --- [dummy arguments]
!!
!      real*8 ,intent(out) :: cord(3,mxvrtx)
!      integer,intent(out) :: lacell(mxcelb)
!      integer,intent(out) :: lvcell(8,mxcell)
!      integer,intent(out) :: ierror
!!
!      integer :: ios=0
!      ierror=0
!! --- 
!      if(gdformat=='GF') then
!!	call read_hpc_GF(cord,lacell,lvcell,ios)
!	if(ios/=0) then
!          write(ifle,*) 'ERR: reading read_hpc_GF'
!	  ierror=1
!	  call FFRABORT(1,'read_griddata_hpc')
!	endif
!      else
!! --- Other mesher grid data has same format in HPC
!        call read_hpc(cord,lacell,lvcell,ios)
!        call comm       
!	if(ios/=0) then
!          write(ifle,*) 'ERR: reading read_hpc'
!	  ierror=1
!          call FFRABORT(1,'read_griddata_hpc')
!	endif        
!      end if
!!
!      return
!!
!      end subroutine  read_griddata_hpc
!
!!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      subroutine read_hpc(cord,lacell,lvcell,ierror)
!!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!     Get following parameter : ncell,nvrtx
!!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!! --- [module arguments]
!!
!      use module_dimension
!      use module_hpcutil
!      use module_io      ,only : lenfnm,ifle,getfil,ifll
!      use module_boundary,only : nbcnd,boundName,numbname,ivbcnd,
!     &                           lvbcnd,icbinr,lcbinr,
!     &                           kdbcnd,kdprdc,prdcAxis,prdcOffset,
!     &                           boundIDMap,SFBOUN,NBOUND
!      implicit none
!!
!! --- [dummy arguments]
!!
!      integer,intent(out)     :: lacell(  mxcell)
!      integer,intent(out)     :: lvcell(8,mxcell)
!      real*8 ,intent(out)     :: cord  (3,mxvrtx)
!      integer,intent(inout)   :: ierror
!!
!! --- [local entities]
!!
!      character(len=80)   :: HEADER1,HEADER2
!      integer,allocatable :: prdcMap(:)
!      integer :: nn0,nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8
!      integer :: NNX,i,j,k,ivv,icL
!      integer :: ic,iv
!      integer :: ierr1,ihpc,ifchl
!      integer,allocatable :: nobdmn(:)
!      integer,allocatable :: ivbdmn(:)
!      integer,allocatable :: lvbdmn(:)
!      integer,allocatable :: icbdmn(:)
!      integer,allocatable :: lcbdmn(:)
!      real(8),parameter   :: smallNrm=1.0d-13
!      integer             :: boundIDa,boundIDb
!      integer             :: tempInt,M,Int,nb,mb,IS
      
!
!
!      ifchl=11!my_rank+11
!      open (ifchl,file=GRIDin, status= 'unknown')
!      read (ifchl,*,ERR=9998) HEADER1
!      read (ifchl,*,ERR=9998) HEADER2
!      read (ifchl,*,ERR=9998) nn3, NNX
!      if(NNX.ne.NVRTX_L) then
!        write(ifle,*) ' #### ERR: NNX NOT equal to NVRTX_L',NNX,NVRTX_L
!        WRITE(ifle,*) ' ####    Please run prefflow program again'
!        call FFRABORT(1,'read_hpc')
!      else
!        nvrtx=NNX
!      endif
!      read (ifchl,'(5e14.6)') ((cord(k,i),k=1,nn3),i=1,NNX)
!      cord=DBLE(cord)
!! ---  Local CONNECTIVITY
!      read (ifchl,*,ERR=9998) HEADER2
!      read (ifchl,*,ERR=9998) nn8, NNX
!      if(NNX.ne.NCELL_L) then
!        write(ifle,*) ' #### ERR: NNX NOT equal to NCELL_L',NNX,NCELL_L
!        WRITE(ifle,*) ' ####    Please run prefflow program again'
!        call FFRABORT(1,'read_hpc')
!      else
!        ncell=NNX
!      endif
!      read (ifchl,'(6i12)',ERR=9998) ((lvcell(k,i),k=1,nn8),i=1,NNX)
!! --- Material Type
!      read (ifchl,*,ERR=9998) HEADER2
!      read (ifchl,*,ERR=9998) nn1,NNX
!      read (ifchl,'(6i12)',ERR=9998) (lacell(icL),icL=1,NNX)
!      read (ifchl,*,ERR=9998) HEADER2
!!
!      close (ifchl)
!!
!! --- Read BC file
!!
!      open (ifchl,file=BCin,status='unknown')
!      read (ifchl,*,ERR=9998) HEADER1
!      read (ifchl,*,ERR=9998) HEADER2
!      read (ifchl,*,ERR=9998) NBOUND
!
!      allocate(SFBOUN(NBOUND))
!      allocate(boundIDMap(nbcnd,numbname))
!      allocate( nobdmn(NBOUND),
!     &          ivbdmn(0:NBOUND), lvbdmn(MXBCTOT),
!     &          icbdmn(0:NBOUND), lcbdmn(NBOUND),
!     &          stat=ierr1 )
!      allocate(ivbcnd(0:nbcnd), lvbcnd(MXBCTOT))
!      nobdmn=0
!      ivbdmn=0
!      lvbdmn=0
!      ivbcnd=0
!      lvbcnd=0
!
!      NBCTOT=0
!      do 100 mb=1,NBOUND
!      read (ifchl,*,ERR=9998) SFBOUN(mb)
!      read (ifchl,*,ERR=9998) nn1,ivv
!      ivbdmn(mb)=ivbdmn(mb-1)+ivv
!      read (ifchl,'(6i12)',ERR=9998) 
!     &            (lvbdmn(i),i=ivbdmn(mb-1)+1,ivbdmn(mb))
!      NBCTOT=NBCTOT+ivv
! 100  continue
!
!      close(ifchl)
!
!      IF(NBCTOT.NE.NBCTOT_L) then
!        write(ifle,*) ' #### ERR: BC Vertex Number error'
!        write(ifle,*) ' NBCTOT, NBCTOT_L= ',NBCTOT, NBCTOT_L
!        WRITE(ifle,*) ' ####    Please run prefflow program again'
!        call FFRABORT(1,'read_hpc')
!      endif
!
!
!      boundIDMap=-1
!      do j=1,numbname
!      do nb=1,nbcnd
!      if(boundName(nb,j).ne.' ') then
!        do mb=1,NBOUND
!        if(adjustl(boundName(nb,j)).eq.adjustl(SFBOUN(mb))) then
!          boundIDMap(nb,j)=mb
!          exit
!        endif
!	end do
!      end if
!      end do
!      end do
!
! --- Reording lvbdmn
!
!      do 200 nb=1,nbcnd
!      if(kdbcnd(0,nb).eq.kdprdc) then
!	boundIDa=boundIDMap(nb,1)
! 	boundIDb=boundIDMap(nb,2)
!	allocate(prdcMap(ivbdmn(boundIDb-1)+1:ivbdmn(boundIDb)))
!	prdcMap=-1
!	if((prdcAxis(nb)=='x').or.(prdcAxis(nb)=='X')) then
!	  do 210 j=ivbdmn(boundIDa-1)+1,ivbdmn(boundIDa)
!	  do 220 k=ivbdmn(boundIDb-1)+1,ivbdmn(boundIDb)
!          if(((cord(2,lvbdmn(k))-cord(2,lvbdmn(j))-prdcOffset(nb,2))**2
!     &       +(cord(3,lvbdmn(k))-cord(3,lvbdmn(j))-prdcOffset(nb,3))**2
!     &        ).lt.smallNrm) then
!	    prdcMap(ivbdmn(boundIDb-1)+j-ivbdmn(boundIDa-1))
!     &      = lvbdmn(k)
!	    exit
!	  end if
! 220      continue
! 210      continue
!	else if((prdcAxis(nb)=='y').or.(prdcAxis(nb)=='Y')) then
!	  do 230 j=ivbdmn(boundIDa-1)+1,ivbdmn(boundIDa)
!	  do 240 k=ivbdmn(boundIDb-1)+1,ivbdmn(boundIDb)
!	  if(((cord(1,lvbdmn(k))-cord(1,lvbdmn(j))- 
!     &         prdcOffset(nb,1))**2+
!     &        (cord(3,lvbdmn(k))-cord(3,lvbdmn(j))- 
!     &         prdcOffset(nb,3))**2)<smallNrm) then
!	    prdcMap(ivbdmn(boundIDb-1)+j-ivbdmn(boundIDa-1))
!     &      = lvbdmn(k)
!  	    exit
!          end if
! 240      continue
! 230      continue
!        else if((prdcAxis(nb)=='z').or.(prdcAxis(nb)=='Z')) then
!          do 250 j=ivbdmn(boundIDa-1)+1,ivbdmn(boundIDa)
!          do 260 k=ivbdmn(boundIDb-1)+1,ivbdmn(boundIDb)
!	  if(((cord(2,lvbdmn(k))-cord(2,lvbdmn(j))- 
!     &         prdcOffset(nb,2))**2+ 
!     &        (cord(1,lvbdmn(k))-cord(1,lvbdmn(j))- 
!     &         prdcOffset(nb,1))**2)<smallNrm) then
!	    prdcMap(ivbdmn(boundIDb-1)+j-ivbdmn(boundIDa-1))
!     &      = lvbdmn(k)
!	    exit
!	  end if
! 260      continue
! 250      continue
!	else
!	  WRITE(ifle,*)
!     &    '(readhpc)Periodic axis setting is wrong.',
!     &    prdcAxis(nb)
!	  call FFRABORT(1,'read_hpc')
!	end if
!!
!        lvbdmn(ivbdmn(boundIDb-1)+1:ivbdmn(boundIDb))=
!     &  prdcMap(ivbdmn(boundIDb-1)+1:ivbdmn(boundIDb))
!!
!	deallocate(prdcMap)
!      end if
! 200  continue
!
! --- 
!
!      do 400 nb=1,nbcnd
!      boundIDa=boundIDMap(nb,1)
!      boundIDb=boundIDMap(nb,2)
!      ivbcnd(nb)=ivbcnd(nb-1)
!      ivbcnd(nb)=ivbcnd(nb)+ivbdmn(boundIDa)-ivbdmn(boundIDa-1)
!      if((ivbcnd(nb)-ivbcnd(nb-1))/=(ivbdmn(boundIDa) 
!     &  -ivbdmn(boundIDa-1))) then
!	WRITE(ifle,*)  
!     &  '(readhpc)',
!     &  'Number of periodic boundary vertex is not matched.'
!	call FFRABORT(1,'read_hpc')
!      end if
!      lvbcnd(ivbcnd(nb-1)+1:ivbcnd(nb))=
!     &lvbdmn(ivbdmn(boundIDa-1)+1:ivbdmn(boundIDa))
!      if(kdbcnd(0,nb)==kdprdc) then
!	tempInt=ivbcnd(nb)
!	ivbcnd(nb)=ivbcnd(nb)+ivbdmn(boundIDb)-ivbdmn(boundIDb-1)
!	lvbcnd(tempInt+1:ivbcnd(nb))= 
!     &  lvbdmn(ivbdmn(boundIDb-1)+1:ivbdmn(boundIDb))
!      end if
! 400  continue
!
!      deallocate(nobdmn,ivbdmn,lvbdmn,icbdmn,lcbdmn)
!
!      return

! 9998 continue
!      write(ifle,*)'ERR : read grid data in read_hpc'
!      call FFRABORT(1,'read_hpc')
! 9999 continue
!      write(ifle,*) '(read_hpc)'
!      ierror=1
!
!      end subroutine read_hpc
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine comm
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Get following parameter: NCVIN,NCV
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_dimension
      use module_hpcutil
      use module_io,only : ifle,ifll
      use module_metrix,only : WS,WR,sta1,sta2,req1,req2
!
! --- [local entities]
!
      implicit none
      integer :: nn,nn0,nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,ifchl
      integer :: i,j,k,ICV
!
!
! --- Communication file input
!
      ifchl=11!my_rank+11
      open (ifchl, file= COMMin, form='formatted', status='unknown')
!
! --- NEIBPE
!
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) NEIBPETOT
      allocate (NEIBPE(NEIBPETOT))
      allocate (IMPORT_index(0:NEIBPETOT))
      allocate (EXPORT_index(0:NEIBPETOT))
      !
!
      NEIBPE= 0
      IMPORT_index= 0
      EXPORT_index= 0
!      IMPORT_indexIP= 0
!      EXPORT_indexIP= 0
!2
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) (NEIBPE(k), k= 1, NEIBPETOT) !1
!
! --- IMPORT
!3 
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) (IMPORT_index(k),k=1,NEIBPETOT) !2
      nn= IMPORT_index(NEIBPETOT)
      allocate (IMPORT_item(nn))
      read (ifchl,'(a80)' ,ERR=9998) LINE
!4
      read (ifchl,'(6i12)',ERR=9998) (IMPORT_item(k), k= 1, nn) !3
!
! --- EXPORT
!4
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) (EXPORT_index(k),k=1,NEIBPETOT) !4
      nn= EXPORT_index(NEIBPETOT)
!     
      allocate (EXPORT_item(nn))
!      
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) (EXPORT_item(k), k= 1, nn) !5
!
! --- NODE info.
! --- inter vertex
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) INTNODTOT
!
      if(INTNODTOT.ne.NCVIN_L) then
        write(ifle,*) ' #### ERR: NCVIN NOT equal to NCVIN_L',NCVIN,
     &                  NCVIN_L
        WRITE(ifle,*) ' ####    Please run prefflow program again'
        call FFRABORT(1,'comm')
      else
        NCVIN=INTNODTOT
      endif
!
! --- inter+exter vertex
!
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) NODTOT
      if(NODTOT.ne.NCV_L) then
        write(ifle,*) ' #### ERR: NCV NOT equal to NCV_L',NCV,NCV_L
        WRITE(ifle,*) ' ####    Please run prefflow program again'
        call FFRABORT(1,'comm')
      else
        NCV=NODTOT
      endif
!
      allocate (NOD_IDg(NODTOT))
      read (ifchl,'(a80)' ,ERR=9998) LINE
      read (ifchl,'(6i12)',ERR=9998) (NOD_IDg(ICV),ICV=1,NODTOT) !6
!
! --- ELEMENT info.
!      read (ifchl,'(a80)' ,ERR=9998) LINE
!      read (ifchl,'(6i12)',ERR=9998) ELMtot
!
!      allocate (ELM_IDg(ELMtot))
!      read (ifchl,'(a80)' ,ERR=9998) LINE
!      read (ifchl,'(6i12)',ERR=9998) (ELM_IDg(k), k= 1, ELMtot) !7
! --- ADDITIONAL info.
      LENT = max (NODTOT, IMPORT_index(NEIBPETOT),
     &                       EXPORT_index(NEIBPETOT))
      NODtotDIFF= INTNODTOT - NODTOT
!
      close (ifchl)
!
! DEBUG
      allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
      allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
      allocate (req1(NEIBPETOT))
      allocate (req2(NEIBPETOT))
      allocate(WS(max(1,EXPORT_index(NEIBPETOT))),
     &         WR(MAX(1,IMPORT_index(NEIBPETOT))))
!
      return
!
 9998 call FFRABORT(1,'comm')
!
      end subroutine comm
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_igather(ncpu,send_data,recv_buff)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(in)    :: ncpu
      integer,intent(in)    :: send_data
      integer,intent(inout) :: recv_buff(1:ncpu)
      integer :: ierr
!
      call mpi_allgather(send_data,1,MPI_INTEGER,
     &                recv_buff,1,MPI_INTEGER,SOLVER_COMM,ierr)
      end subroutine hpc_igather

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_dgather(ncpu,send_data,recv_buff)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
!      
      integer,intent(in)    :: ncpu
      real(8),intent(in)    :: send_data
      real(8),intent(inout) :: recv_buff(1:ncpu)
      integer :: ierr
!
      call mpi_allgather(send_data,1,MPI_DOUBLE_PRECISION,
     &                   recv_buff,1,MPI_DOUBLE_PRECISION,
     &                   SOLVER_COMM,ierr)
!
      end subroutine hpc_dgather
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_dsend(ndim,send_buff,dest,tagID)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
!
      integer,intent(in) :: ndim,dest,tagID
      real(8),intent(in) :: send_buff(1:ndim)
      integer :: ierr
!      
      call mpi_send(send_buff,ndim,MPI_DOUBLE_PRECISION,dest,tagID,
     &                  SOLVER_COMM,ierr)
!      
      end subroutine hpc_dsend

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_drecv(ndim,recv_buff,source,tagID)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
!
      integer,intent(in)    :: ndim,source,tagID
      real(8),intent(inout) :: recv_buff(1:ndim)
      integer :: ierr
      INTEGER :: STATUS(MPI_STATUS_SIZE)
!      
      call mpi_recv(recv_buff,ndim,MPI_DOUBLE_PRECISION,source,tagID,
     &                  SOLVER_COMM,STATUS,ierr)
!     
      end subroutine hpc_drecv
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcbarrier(IERR)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(inout) :: ierr
!
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      end subroutine hpcbarrier
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PARTICLE_MPI(IMODE,
     &                 MP,MCELL,
     &                 MFACE,N_injctr,
     &                 XYZP,UVWP,DIA,TMPP,PYS,
     *                 FU,FV,FW,RHOP,PARCEL,PMASS,PMASSC,INSIDE,
     *                 NEIBCELL,TRESIDUAL,
     *                 NPEXCHAGOUT,NPEXCHAGIN,IPSND,IPQUENE,
     &                 dPDstion,dPDstvel,
     &                 NINALL,JOUT,IPCELL,!IPCPU,
     6                 NPOUTALL,NPINALL
     &   )
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C----------------------
      use module_hpcutil
      USE module_dimension 
!      use module_particle ,only : NPALL
!
!----------------------
      IMPLICIT NONE
!----------------------
! --- [Dummy Auguments]
C----------------------
      INTEGER,intent(in)   :: MP,MCELL,MFACE,N_injctr,IMODE
      REAL*8, intent(inout):: XYZP(MP,3)
      REAL*8, intent(inout):: UVWP(MP,3)
      REAL*8, intent(inOUT):: DIA(MP) 
      REAL*8,INTENT(INOUT) :: PYS(MP,MXCOMP,2)
      REAL*8, intent(inOUT):: TMPP(MP) 
      REAL*8, intent(inout):: FU(MP,2)
      REAL*8, intent(inout):: FV(MP,2)
      REAL*8, intent(inout):: FW(MP,2)

      real*8, intent(inOUT):: RHOP( MP )
      REAL*8, intent(inOUT):: PARCEL(MP)
      INTEGER,intent(inout):: INSIDE(MP)
      REAL*8, intent(inOUT):: PMASS(MP)
      REAL*8, intent(inOUT):: PMASSC(MP)
      
      REAL*8, intent(inout):: TRESIDUAL(MP)
      INTEGER,INTENT(INOUT):: NEIBCELL(MCELL)
      INTEGER,INTENT(INOUT):: IPCELL(MCELL)
      INTEGER,INTENT(INOUT):: IPSND(MP),IPQUENE(MP)
      INTEGER,INTENT(INOUT):: NPEXCHAGOUT(NEIBPETOT),
     *                        NPEXCHAGIN(NEIBPETOT)
      REAL*8, intent(inout):: dPDstion(MP)
      REAL*8, intent(inout):: dPDstvel(MP)
      INTEGER,INTENT(INOUT):: NINALL(0:N_injctr)
!      INTEGER,INTENT(INOUT):: IPCPU(MP)
!
      INTEGER,INTENT(IN)::
     &                        NPOUTALL,
     &                        NPINALL
      INTEGER,intent(inout)::JOUT(MP)      
!----------------
! --- local 
!----------------
      real*8 ,save,allocatable  :: BUFSND(:),BUFRCV(:)
      INTEGER :: INI,IPS,IPE,NNDEL,MXBUFO,MXBUFI,MAXBUF
!
C************************************************************************
      INTEGER            :: iSNDRCVSIZE=20,NPIN
      INTEGER            :: isca=20
      INTEGER, PARAMETER :: MAXDOM=100
      real*8             :: PI
C
C************************************************************************
C
!---------------------------------------
! --- [local entities]
!---------------------------------------
      integer,save, allocatable :: sta1(:,:)
      integer,save, allocatable :: sta2(:,:)
      integer,save, allocatable :: req1(:  )
      integer,save, allocatable :: req2(:  )
      integer                   :: NCPU,neib,istart,inum,MSGLEN,k,ICOM,
     &                             IPP
      integer,save              :: NFLAG=0,NCON,I_RANK
      integer:: IERR,I,IP,NSTART,IPCV,MSGTYP,ISEND,IRECV,N
      integer:: req
!------------------------
! --- INIT.
!------------------------
      if(IMODE==1) then
        iSNDRCVSIZE=isca+NCOMP
      else
        iSNDRCVSIZE=1
      endif
!
      PI=4.d0*DATAN(1.d0)
      if(NFLAG==0) then
        allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (req1(NEIBPETOT))
        allocate (req2(NEIBPETOT))
        sta1=0
        sta2=0
        req1=0
        req2=0
        NFLAG=1
        call hpcisum(NFLAG)
      endif
!-------------------
! --- 
!-------------------
      EXPORT_indexIP(0)=0 
      IMPORT_indexIP(0)=0 

! --- 
!
      DO NCPU=1,NEIBPETOT 
      EXPORT_indexIP(NCPU)=EXPORT_indexIP(NCPU-1) 
     *               +NPEXCHAGOUT(NCPU)*iSNDRCVSIZE
      IMPORT_indexIP(NCPU)=IMPORT_indexIP(NCPU-1) 
     *               +NPEXCHAGIN(NCPU)*iSNDRCVSIZE
      ENDDO 
!
      IERR=0
!
      MXBUFO=iSNDRCVSIZE*(NPOUTALL)
      MXBUFI=(NPINALL)*iSNDRCVSIZE
      call  hpcimax(MXBUFO)
      call  hpcimax(MXBUFI)
!
      MAXBUF=MAX(MAXBUF,MP,MXBUFO,MXBUFI)
!      
      allocate(BUFSND(MAXBUF),BUFRCV(MAXBUF))
      BUFSND=0.d0
      BUFRCV=0.d0
!

!---------------------------------------------
! --- CHECK THE INTERNAL ARRAY SIZE
!---------------------------------------------
!
!--------------------------------------
! POST ALL THE EXPECTED RECEIVES 
!--------------------------------------
      CALL MPI_BARRIER(SOLVER_COMM,IERR) !1111
!--------------------------------------
! --- SET UP THE SEND BUFFER
!--------------------------------------
      if(IMODE==2) then
        if(NPOUTALL>0) then 
          DO I=1,NPOUTALL 
          IP=IPSND(I)
          NCPU=JOUT(IP)                 !I_RANK=IPCPU(IP) 
          NSTART=EXPORT_indexIP(NCPU-1)+(IPQUENE(I)-1)*iSNDRCVSIZE
          BUFSND(NSTART+1)=dble(IP)
          JOUT(IP)=HPC_dead
          IPSND(I)=0
          enddo
        endif
      else
!
      if(NPOUTALL>0) then 
        DO I=1,NPOUTALL
        IP=IPSND(I)
        NCPU=JOUT(IP)
        
        JOUT(IP)=HPC_dead
        NSTART=EXPORT_indexIP(NCPU-1)+(IPQUENE(I)-1)*iSNDRCVSIZE

        IPCV=INSIDE(IP)
        IF(NSTART+iSNDRCVSIZE>MAXBUF) then
          write(*,*) 'MSG:-0',my_rank,'',NCPU-1
          write(*,*) 'MSG:-1',IPQUENE(I),EXPORT_indexIP(NCPU-1)
          write(*,*) 'MSG:-2',NSTART,iSNDRCVSIZE,MAXBUF
          call FFRABORT(1,'ERR_MSG-1:MAXBUF too small')
        endif
!

        NCON=1
        BUFSND(NSTART+1)  = dble(NEIBCELL(IPCV))  !ICV_NB
        NCON=NCON+1      
        BUFSND(NSTART+2)  = dble(IP)
        NCON=NCON+1      
        BUFSND(NSTART+3)  = DIA(IP)
        NCON=NCON+1      
        BUFSND(NSTART+4)  = PARCEL(IP)
        NCON=NCON+1      
        BUFSND(NSTART+5)  = RHOP(IP)
        NCON=NCON+1                    
        BUFSND(NSTART+6)  = XYZP(IP,1)
        NCON=NCON+1      
        BUFSND(NSTART+7)  = XYZP(IP,2)
        NCON=NCON+1      
        BUFSND(NSTART+8)  = XYZP(IP,3)
        NCON=NCON+1      
        BUFSND(NSTART+9)  = UVWP(IP,1)
        NCON=NCON+1      
        BUFSND(NSTART+10) = UVWP(IP,2)
        NCON=NCON+1      
        BUFSND(NSTART+11) = UVWP(IP,3)
        NCON=NCON+1      
        BUFSND(NSTART+12) = FU(IP,2)
        NCON=NCON+1      
        BUFSND(NSTART+13) = FV(IP,2)
        NCON=NCON+1      
        BUFSND(NSTART+14) = FW(IP,2)
        NCON=NCON+1      
        BUFSND(NSTART+15) = TRESIDUAL(IP)
        NCON=NCON+1      
        BUFSND(NSTART+16) = dPDstion(IP)
        NCON=NCON+1      
        BUFSND(NSTART+17) = dPDstvel(IP)
        NCON=NCON+1      
        BUFSND(NSTART+18) = TMPP(IP)
        NCON=NCON+1      
        BUFSND(NSTART+19) = PMASS(IP)
        NCON=NCON+1      
        BUFSND(NSTART+20) = PMASSC(IP)
!

        DO ICOM=1,NCOMP
        NCON=NCON+1      
        BUFSND(NSTART+NCON) =PYS(IP,icom,1)
        enddo
        if(NCON/=iSNDRCVSIZE) then
          call FFRABORT
     &   (1,'ERR: [NCON/=iSNDRCVSIZE] Call your supportor ')
        endif
        ENDDO
      endif
      endif
!
!
C-----------------------------
C --- SEND THE PARTICLES 
C-----------------------------
!
!

      DO 220 NEIB=1,NEIBPETOT
      MSGTYP=0
      ISEND=NEIBPE(NEIB)
      MSGLEN=iSNDRCVSIZE*NPEXCHAGOUT(NEIB)
      NSTART=EXPORT_indexIP(NEIB-1)+1
      IF(MSGLEN/=0) then
        IF(NSTART+MSGLEN-1.GT.MAXBUF) THEN
          CALL FFRABORT(1,'ERR_MSG-2:MAXBUF too small') 
        ENDIF
        CALL MPI_ISEND(BUFSND(NSTART),MSGLEN,MPI_DOUBLE_PRECISION,
     &      ISEND,0,SOLVER_COMM,req1(neib), ierr)
!      else
!        ISEND=MPI_PROC_NULL
!        MSGLEN=1
!        NSTART=1
!        BUFSND(NSTART)=1.d0
      endif
!      CALL MPI_ISEND(BUFSND(NSTART),MSGLEN,MPI_DOUBLE_PRECISION,
!     &      ISEND,0,SOLVER_COMM,req1(neib), ierr)
  220 CONTINUE
!
      DO 110 NEIB=1,NEIBPETOT
      IRECV =NEIBPE(NEIB)
      MSGLEN=iSNDRCVSIZE*NPEXCHAGIN(NEIB)
      NSTART=IMPORT_indexIP(NEIB-1)+1
!
      if(MSGLEN/=0) then
        IF(NSTART+MSGLEN-1.GT.MAXBUF) THEN
          CALL FFRABORT(1,'ERR_MSG-3:MAXBUF too small')
        ENDIF
        CALL MPI_IRECV(BUFRCV(NSTART),MSGLEN,MPI_DOUBLE_PRECISION,
     &     IRECV,0,SOLVER_COMM,req2(NEIB),IERR) !SOLVER_COMM
!      else
!        IRECV=MPI_PROC_NULL
!        MSGLEN=1
!        NSTART=1
!        BUFSND(NSTART)=1.d0
      endif
!      CALL MPI_IRECV(BUFRCV(NSTART),MSGLEN,MPI_DOUBLE_PRECISION,
!     &     IRECV,0,SOLVER_COMM,req2(NEIB),IERR) !SOLVER_COMM
  110 ENDDO
!
C------------------
! --- WAIT 
!------------------
      call MPI_WAITall (NEIBPETOT, req2, sta2, ierr)
      call MPI_WAITall (NEIBPETOT, req1, sta1, ierr)
C------------------------------------
C     ADD THE RECEIVED PARTICLES 
C------------------------------------
      if(IMODE==2) then
        if(NPINALL>=1) THEN 
          DO IPP=1,NPINALL
          NSTART=(IPP-1)*iSNDRCVSIZE
          IP=INT(BUFRCV(NSTART+1))
          JOUT(IP)=-MXCVFAC-1
          IPSND(IP)=-1
          enddo
        endif
      else
      if(NPINALL>=1) THEN 
        DO IPP=1,NPINALL
        NSTART=(IPP-1)*iSNDRCVSIZE
! --- ---
        if(NSTART+iSNDRCVSIZE>MXBUFI) 
     &     call FFRABORT(1,'ERR: >MXBUFI') 
!
        IP=INT(BUFRCV(NSTART+2))
        JOUT(IP)=0
        INSIDE(IP)=INT(BUFRCV(NSTART+1)) 
        DIA(IP)=BUFRCV(NSTART+3)

        PARCEL(IP)=BUFRCV(NSTART+4)
        RHOP(IP)=BUFRCV(NSTART+5)
        XYZP(IP,1)=BUFRCV(NSTART+6)
        XYZP(IP,2)=BUFRCV(NSTART+7)
        XYZP(IP,3)=BUFRCV(NSTART+8)
        UVWP(IP,1)=BUFRCV(NSTART+9)
        UVWP(IP,2)=BUFRCV(NSTART+10)
        UVWP(IP,3)=BUFRCV(NSTART+11)
!
        FU(IP,2)=BUFRCV(NSTART+12)
        FV(IP,2)=BUFRCV(NSTART+13)
        FW(IP,2)=BUFRCV(NSTART+14)
!
        TRESIDUAL(IP)=BUFRCV(NSTART+15)
        dPDstion(IP)=BUFRCV(NSTART+16)
        dPDstvel(IP)=BUFRCV(NSTART+17)

        TMPP(IP)=BUFRCV(NSTART+18)
        PMASS(IP)=BUFRCV(NSTART+19)

        PMASSC(IP)=BUFRCV(NSTART+20)

        NCON=isca
        DO ICOM=1,NCOMP
        NCON=NCON+1
        PYS(IP,icom,1)=BUFRCV(NSTART+NCON)
        enddo

        IF(NCON>iSNDRCVSIZE) call FFRABORT(1,'ERR: NCON>iSNDRCVSIZE')
        ENDDO
      ENDIF
      endif
!
      deallocate(BUFSND,BUFRCV) !1111
!
      return
!
      END SUBROUTINE PARTICLE_MPI
C
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine PARTICLE_SEND_RECV02(NPEXCHAGIN,NPEXCHAGOUT)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)  :: NPEXCHAGOUT(NEIBPETOT)
      integer,intent(out) :: NPEXCHAGIN(NEIBPETOT)
!
! --- [local entities]
!
      integer,save, allocatable :: sta1(:,:)
      integer,save, allocatable :: sta2(:,:)
      integer,save, allocatable :: req1(:  )
      integer,save, allocatable :: req2(:  )
      real*8,save ,allocatable       :: WS(:),WR(:)
      integer                   :: neib,istart,inum,k,idim
      integer,save              :: NFLAG=0
      integer  :: ierr
!
! --- INIT.
!
      if (NFLAG.eq.0) then
        allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (req1(NEIBPETOT))
        allocate (req2(NEIBPETOT))
        allocate(WS(NEIBPETOT),WR(NEIBPETOT))
        NFLAG=1 
      endif
! --- 

!
      WS=0.d0
      WR=0.d0
!
! --- SEND
!
      
      do neib=1,NEIBPETOT
        WS(neib)=dble(NPEXCHAGOUT(NEIB))+0.1d0
        inum=1
        call MPI_ISEND (WS(neib), inum, MPI_DOUBLE_PRECISION,
     &                  NEIBPE(neib), 0, SOLVER_COMM,
     &                  req1(neib), ierr)
      enddo
!
! --- RECEIVE 
!
      do neib= 1, NEIBPETOT
        inum  =1
        call MPI_IRECV (WR(neib), inum, MPI_DOUBLE_PRECISION,
     &                  NEIBPE(neib), 0, SOLVER_COMM,
     &                  req2(neib), ierr)
      enddo
!
! --- WAIT
!
      call MPI_WAITALL (NEIBPETOT,req1,sta1,ierr)
      call MPI_WAITALL (NEIBPETOT,req2,sta2, ierr)
!
      do neib=1,NEIBPETOT
        NPEXCHAGIN(NEIB)= INT(WR(neib))
      enddo
!
!      deallocate(WS,WR)
!
      return
      end subroutine PARTICLE_SEND_RECV02
      
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine PARTICLE_SEND_RECV01(MXNP,NP,NEIBCELL,IPCELL)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil
!
      implicit none
!
! --- [dummy arguments]
!
      integer,intent(in)   :: MXNP,NP
      integer,intent(inout):: NEIBCELL(MXNP),IPCELL(MXNP)
!
! --- [local entities]
!
      integer,save, allocatable :: sta1(:,:)
      integer,save, allocatable :: sta2(:,:)
      integer,save, allocatable :: req1(:  )
      integer,save, allocatable :: req2(:  )
      real*8 ,allocatable       :: WS(:),WR(:)
      integer                   :: ncpu,istart,inum,k,idim,neib
      integer,save              :: NFLAG
      integer  :: ierr
      data NFLAG/0/
!
! --- INIT.
!
      if (NFLAG.eq.0) then
        allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (req1(NEIBPETOT))
        allocate (req2(NEIBPETOT))
        NFLAG=1 
      endif
! --- 
      allocate(WS(EXPORT_index(NEIBPETOT)),
     &         WR(IMPORT_index(NEIBPETOT)))
!
C
      WS=0.d0
      WR=0.d0
!
! --- SEND
!
      do neib= 1, NEIBPETOT
        istart= EXPORT_index(neib-1)
        inum  = EXPORT_index(neib  ) - istart
        do k=istart+1,istart+inum
           WS(k)= EXPORT_item(k)+0.1
        enddo
        call MPI_ISEND (WS(istart+1), inum, MPI_DOUBLE_PRECISION,       
     &                  NEIBPE(neib), 0, SOLVER_COMM,                   
     &                  req1(neib), ierr)
      enddo
!
      
! --- RECEIVE
!
      do neib= 1, NEIBPETOT
        istart=IMPORT_index(neib-1)
        inum  =IMPORT_index(neib)-istart
        call MPI_IRECV (WR(istart+1), inum, MPI_DOUBLE_PRECISION,
     &                  NEIBPE(neib), 0, SOLVER_COMM,
     &                  req2(neib), ierr)
      enddo
!
! --- WAIT
!
      call MPI_WAITALL(NEIBPETOT,req1,sta1,ierr)
      call MPI_WAITALL(NEIBPETOT,req2,sta2,ierr)
!
! --- 
!
      do ncpu=1,NEIBPETOT 
      istart=IMPORT_index(ncpu-1)
      inum=IMPORT_index(ncpu)-istart
      do k=istart+1,istart+inum
      NEIBCELL(IMPORT_item(K))=INT(WR(k))  !neib_number
      IPCELL(IMPORT_item(K))=ncpu
!NEIBPE(ncpu)
!ncpu
      enddo
      enddo
!
!
      deallocate(WS,WR)
!
      return
      end subroutine PARTICLE_send_recv01
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrbcast_stream(icpu,RIN,NUM)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- one-to-all
      use module_hpcutil
!
! --- [dummy arguments]
!
      implicit none
      integer,intent(in)     :: icpu
      integer,intent(in)     :: NUM
      real*8 ,intent(inout)  :: RIN(NUM)
      integer                :: ierr 
!
      call MPI_BCAST(RIN,NUM,MPI_DOUBLE_PRECISION,
     +               icpu,SOLVER_COMM,ierr)
!
      return
!
      end subroutine hpcrbcast_stream
!
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MERGEOUT_PTRESULT_MPI(
     &   ncomp,N0,OUTNAME,MP,NPALL,TIME,N_injctr,NINALL,
     +   XYZP,UVWP,DIA,PARCEL,HIS_P,PMASS,PMASSC,
     &   TMPP,PYS,JOUT,IPSND,
     &   iPartLU,dYP0,dZP0,dXP0,P_MF,starttime
     &)
CSS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     
      USE module_dimension,only : MXCOMP,HPC_dead
      USE module_hpcutil
      use module_species,only  : spcnam
      use module_io,only       : getfil,ifle,ifll,gdScale
!      
      IMPLICIT NONE
!
! --- [dummy arguments]
!
      CHARACTER(*),INTENT(IN) :: OUTNAME
      INTEGER,INTENT(IN) :: iPartLU
      REAL*8, INTENT(IN) :: dYP0(N_injctr),
     &                      dZP0(N_injctr),
     &                      dXP0(N_injctr),
     &                      starttime(N_injctr)
      INTEGER,INTENT(IN) :: P_MF(N_injctr)
      INTEGER,INTENT(IN) :: MP,N_injctr,N0,ncomp
      INTEGER,INTENT(IN) :: NPALL(N_injctr)
      INTEGER,INTENT(IN) :: NINALL(0:N_injctr)
      REAL*8, INTENT(IN) :: TIME
      REAL*8, INTENT(IN) :: XYZP(MP,3)
      REAL*8, INTENT(IN) :: UVWP(MP,3)
      REAL*8, INTENT(IN) :: DIA(MP)
      REAL*8, INTENT(IN) :: PARCEL(MP)
      REAL*8, INTENT(IN) :: HIS_P(MP)
      REAL*8, INTENT(IN) :: PMASS(MP)
      REAL*8, INTENT(IN) :: PMASSC(MP)
      REAL*8, INTENT(IN) :: PYS(MP,MXCOMP)
      REAL*8, INTENT(IN) :: TMPP(MP)
      INTEGER,INTENT(IN) :: JOUT(MP)
      INTEGER,INTENT(INout) :: IPSND(MP)
!
! --- [local entities]
!
      integer :: ierr,I,J,K,L,M,N,INI,N00
      integer,save :: IFLAG=0,NPMAXDIM=0
      integer :: icom,IP,IPS,IPE,idum,NNN,NPCPU_TOL,IV
      character*128:: formt
      real*8  :: dum1,dum2,dum3
!
      integer,allocatable,save  :: MPI_INT_SNDBUF(:,:)
      integer,allocatable,save  :: MPI_INT_RCVBUF(:,:)
!
      real*8 ,allocatable,save  :: MPI_REAL_SNDBUF(:)
      real*8 ,allocatable,save  :: MPI_REAL_RCVBUF(:)
!
      real*8,allocatable,save  :: TEMP_P(:,:)
!
      integer,allocatable,save :: NPALL_TOTAL(:)
      integer,allocatable,save :: NPALL_MAX(:)
!
      CHARACTER*6  :: for1
!
      integer :: NPCPU_MAX=0,icpu,MIP,NINJ
      integer,save :: writeflg=0,ierr1=0
!
!
!
!
C--------------------------------------------------
C Get the number of particles in other process.
C--------------------------------------------------
      if(IFLAG==0) then
         IFLAG=1
         ALLOCATE(MPI_INT_SNDBUF(2,N_injctr),stat=ierr1)
         if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE MPI_INT_SNDBUF')

         ALLOCATE(MPI_INT_RCVBUF(1:NPE*2,N_injctr),stat=ierr1)
         if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE MPI_INT_RCVBUF')

         ALLOCATE(NPALL_TOTAL(N_injctr),stat=ierr1)
         if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE NPALL_TOTAL')

         ALLOCATE(NPALL_MAX(0:N_injctr),stat=ierr1)
         if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE NPALL_MAX')

         ALLOCATE(MPI_REAL_SNDBUF(MP),stat=ierr1)
         if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE MPI_REAL_SNDBUF')

!         ALLOCATE(MPI_REAL_RCVBUF(MP*NPE),stat=ierr1)
!         if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE MPI_REAL_RCVBUF')

         ALLOCATE(TEMP_P(MP,N0),stat=ierr1)
         if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE TEMP_P')

         MPI_INT_SNDBUF=0
         MPI_INT_RCVBUF=0
         MPI_REAL_SNDBUF=0.d0
         MPI_REAL_RCVBUF=0.d0
      endif
!
      NPALL_MAX(:)=0
      NPALL_TOTAL(:)=0
!
      do INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      DO IP=IPS,IPE
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      NPALL_TOTAL(INI)=NPALL_TOTAL(INI)+1
      enddo
      enddo
!
      do INI=1,N_injctr
!
        MPI_INT_SNDBUF(1,INI)=my_rank
        MPI_INT_SNDBUF(2,INI)=NPALL_TOTAL(INI) 
        CALL MPI_GATHER(MPI_INT_SNDBUF(:,INI),
     +                  2,
     +                  MPI_INTEGER,
     &                  MPI_INT_RCVBUF(:,INI),
     +                  2,
     +                  MPI_INTEGER,
     +                  root,
     +                  SOLVER_COMM,
     +                  ierr)
      enddo 

!--------------------------------------
! --- 
!--------------------------------------
      if(my_rank.eq.root) then
        NPALL_TOTAL(:)=0
        NPALL_MAX(:)=-1
        NPCPU_MAX=-1
        NPCPU_TOL=0
!
        do INI=1,N_injctr 
        do icpu=1,NPE
        NNN=MPI_INT_RCVBUF(2*icpu,INI)
        NPALL_MAX(INI)=MAX(NPALL_MAX(INI),NNN)
        NPALL_TOTAL(INI)=NPALL_TOTAL(INI)+NNN
        enddo
        NPCPU_TOL=NPCPU_TOL+NPALL_TOTAL(INI)
        enddo
      endif
!
!
!
      NPMAXDIM=0
      do INI=1,N_injctr
      IDUM=NPALL_MAX(INI)
      CALL MPI_BCAST(IDUM,1,
     +               MPI_INTEGER,root,SOLVER_COMM,ierr)
      NPALL_MAX(INI)=IDUM
      NPMAXDIM=MAX(NPALL_MAX(INI),NPMAXDIM)
      enddo
      NPMAXDIM=NPE*NPMAXDIM
!
      ALLOCATE(MPI_REAL_RCVBUF(NPMAXDIM),stat=ierr1)
      if(ierr1/=0) call FFRABORT(1,'ERR: ALLOCATE MPI_REAL_RCVBUF')
!--------------------------------
! --- 
!--------------------------------

C------------------------
C Gather XYZP data 1,2,3
C------------------------
!      if(my_rank.eq.root) write(ifll,'(2X,a)') "Gathering MPI datas ..."
!
      N00=0
      do IV=1,3
      N00=N00+1
      do 1000 INI=1,N_injctr 
      IPS=NINALL(INI-1)+1 
      IPE=NINALL(INI-1)+NPALL(INI) 
      MPI_REAL_SNDBUF(:)=0.d0 
      N=0
      do IP=IPS,IPE                   !   N=1,NPALL(INI)
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle 
      N=N+1
      MPI_REAL_SNDBUF(N)=XYZP(IP,IV)
      enddo
!
!
      if(NPALL_MAX(INI)/=0) then
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,
     +                  ierr)
      endif
!
      if(my_rank==root) then 
        if(INI==1) M=0 
        do icpu=1,NPE
        do N=1,MPI_INT_RCVBUF(2*icpu,INI)
        M=M+1
        L=(icpu-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
      endif
 1000 enddo
      enddo
C-------------------------
C Gather UVWP data 4,5,6
C-------------------------
      
      do 2100 i=1,3
      N00=N00+1
      do 2000 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      MPI_REAL_SNDBUF(:)=0.d0
      N=0
      do IP=IPS,IPE    ! N=1,NPALL(INI)
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=UVWP(IP,i)
      enddo
      if(NPALL_MAX(INI)/=0) then
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     &                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,
     +                  ierr)
      endif
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
        endif
 2000 enddo
 2100 enddo

C-------------------------
C Gather DIA data 7
C-------------------------
      MPI_REAL_SNDBUF(:)=0.d0
      N00=N00+1
      do 4000 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      N=0
      do IP=IPS,IPE    ! do N=1,NPALL(INI)
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=DIA(IP)
      enddo
      if(NPALL_MAX(INI)/=0) then      
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,
     +                  ierr)
      endif 
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
      endif
 4000 enddo
      MIP=M
!
!-------------------------
! Gather PARCEL data 8
!-------------------------
!

      MPI_REAL_SNDBUF(:)=0.d0
      N00=N00+1
      do 4500 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      N=0
      do IP=IPS,IPE    !
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=HIS_P(IP)
      enddo
      if(NPALL_MAX(INI)/=0) then
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,ierr)
      endif
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
      endif
 4500 enddo

!-------------------------
! Gather PARCEL data 9
!-------------------------

      MPI_REAL_SNDBUF(:)=0.d0
      N00=N00+1
      do 5000 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      N=0
      do IP=IPS,IPE    !
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=PARCEL(IP)
      enddo
      if(NPALL_MAX(INI)/=0) then      
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,ierr)
      endif
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
      endif
 5000 enddo

!
!-------------------------
! Gather PMASS data 10
!-------------------------
      MPI_REAL_SNDBUF(:)=0.d0
      N00=N00+1
      do 5500 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      N=0
      do IP=IPS,IPE    !
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=PMASS(IP)
      enddo
      if(NPALL_MAX(INI)/=0) then
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,ierr)
      endif
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
      endif
 5500 enddo

!-------------------------
! Gather TMPP data 11
!-------------------------
      MPI_REAL_SNDBUF(:)=0.d0
      N00=N00+1
      do 6000 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      N=0
      do IP=IPS,IPE    !
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=TMPP(IP)
      enddo
      if(NPALL_MAX(INI)/=0) then
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,ierr)
      endif
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
      endif
 6000 enddo

!-------------------------
! Gather PYS data 12-17
!-------------------------
      do ICOM=1,NCOMP
      MPI_REAL_SNDBUF(:)=0.d0
      N00=N00+1
      do 6500 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      N=0
      do IP=IPS,IPE    !
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=PYS(IP,ICOM)
      enddo
      if(NPALL_MAX(INI)/=0) then
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,ierr)
      endif
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=MPI_REAL_RCVBUF(L)
        enddo
        enddo
      endif
 6500 enddo
      enddo
      
!-------------------------
! Gather inj_no data
!-------------------------
      MPI_REAL_SNDBUF(:)=0.d0
      N00=N00+1
      NINJ=N00
      do 7000 INI=1,N_injctr
      IPS=NINALL(INI-1)+1
      IPE=NINALL(INI-1)+NPALL(INI)
      N=0
      do IP=IPS,IPE    !
      if(JOUT(IP)>0.or.JOUT(IP)==HPC_dead) cycle
      N=N+1
      MPI_REAL_SNDBUF(N)=dble(INI)
      enddo
      if(NPALL_MAX(INI)/=0) then
        CALL MPI_GATHER(MPI_REAL_SNDBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  MPI_REAL_RCVBUF,
     +                  NPALL_MAX(INI),
     +                  MPI_DOUBLE_PRECISION,
     +                  root,
     +                  SOLVER_COMM,ierr)
      endif
      if(my_rank.eq.root) then
        if(INI==1) M=0
        do j=1,NPE
        do N=1,MPI_INT_RCVBUF(2*j,INI)
        M=M+1
        L=(j-1)*NPALL_MAX(INI)+N
        TEMP_P(M,N00)=(MPI_REAL_RCVBUF(L))
        enddo
        enddo
      endif
 7000 enddo
C-------------------------
C Data gathering end
C-------------------------
      if(my_rank.eq.root) then
        write(for1,'(I3)') ncomp
          write(formt,'(a,a,a)')
     &   '(I10,6(3X,E15.8),5(3X,E15.8),',trim(for1),'(3X,E15.8))'
        OPEN(iPartLU,FILE = OUTNAME,STATUS="REPLACE")
        write(iPartLU,fmt="(2I6,' TIME=',E12.7,
     +        ' X Y Z U V W DIA HIS PARCEL,PMASS,TMP_P,Ys_P')")
     &    N_injctr,N0,time
        WRITE(iPartLU,'(I3)') ncomp
!
        do icom=1,ncomp
        WRITE(iPartLU,'(a)') trim(spcnam(icom))
        enddo
!
        NPALL_MAX=0
        IPSND=0
        do INI=1,N_injctr
        DO M=1,MIP
        IF(INT(TEMP_P(M,NINJ))==INI) then
          NPALL_MAX(INI)=NPALL_MAX(INI)+1
          IDUM=NPALL_MAX(INI-1)+NPALL_MAX(INI)
          IPSND(IDUM)=M 
        endif
        enddo
        NPALL_MAX(INI)=NPALL_MAX(INI-1)+NPALL_MAX(INI)
        enddo

        do INI=1,N_injctr
          dum1=dXP0(INI)/gdScale
          dum2=dYP0(INI)/gdScale
          dum3=dZP0(INI)/gdScale
          IPS=NPALL_MAX(INI-1)+1
          IPE=NPALL_MAX(INI)
          idum=NPALL_MAX(INI)-NPALL_MAX(INI-1)
          write(iPartLU,2501) idum,dum1,dum2,dum3,starttime(INI)
          DO IP=IPS,IPE
            M=IPSND(IP)
            WRITE(iPartLU,formt) M,
     &        (TEMP_P(M,K)/gdScale,K=1,3),
     &        (TEMP_P(M,K),K=4,6),
     &         TEMP_P(M,7),
     &         TEMP_P(M,8),
     &         TEMP_P(M,9),
     &         TEMP_P(M,10),
     &         TEMP_P(M,11),
     &        (TEMP_P(M,11+K),K=1,ncomp)
          ENDDO
        enddo
 2501   FORMAT(I10,4E15.8) 
        CLOSE(iPartLU)
      endif
!
      deALLOCATE(MPI_REAL_RCVBUF)
      RETURN
!
      END SUBROUTINE MERGEOUT_PTRESULT_MPI
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE hpc_MAXLOC(N,ARR,WK,W1K) 
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil 
!
! --- [dummy arguments] 
!
      implicit none 
!
      INTEGER ,intent(inout)  :: N
      real*8  ,intent(inout)  :: ARR(2,N),WK(2,N)
      INTEGER ,intent(inout)  :: W1K(N)
!
      integer  :: ierr=0
!
! --- 
!
      call MPI_REDUCE(ARR,WK,N,MPI_2DOUBLE_PRECISION,
     &                MPI_MAXLOC,ROOT,SOLVER_COMM,ierr) 
! 
      W1K(1:N)=INT(WK(2,1:N))
      call MPI_BCAST 
     &  (W1K,N,MPI_INTEGER,ROOT,SOLVER_COMM,ierr)  
!
      END SUBROUTINE hpc_MAXLOC 
!
