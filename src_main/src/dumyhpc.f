!     The subroutines included in this file is used for dummy HPC
!
!      subroutine HPC_INIT
!      subroutine FFRABORT
!      subroutine FFTIME
!      subroutine SOLVER_SEND_RECV
!      subroutine hpcimax
!      subroutine hpcimax_0
!      subroutine hpcrmax
!      subroutine hpcrmax_0icpu,RIN,NUM)
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
!      subroutine hpclor
!      subroutine hpcland
!      subroutine hpcmaxloc
!      subroutine read_griddata_hpc
!      subroutine comm
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine HPC_INIT(NCPU,my_cpu)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
      integer,intent(out)  :: NCPU,my_cpu
!
      NPE=1
      PETOT=1
      ROOT=0
      my_rank=0
!      NEIBPETOT=0
      RESID=1.D-6
      ITERHPC=1000
      NCPU=1
      my_cpu=0
!
      end subroutine HPC_INIT
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine FFRABORT(ICODE,subname)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
      use module_io,    only : ifll,ifle
!
! --- [dummy arguments]
!
      integer,intent(in)      ::  ICODE
      character(*),intent(in) ::  subname
!
! --- [local entities]
!
        IF(ICODE.EQ.1) THEN 
!
          WRITE(ifle,*) 'STOP Subroutine name: ',
     &               subname(:len_trim(subname))
          STOP 9999
        ELSEIF(ICODE.EQ.0) THEN
          write(ifll,3000)
          STOP
        else
          WRITE(ifle,*) 'STOP ERROR CODE NO.= ',ICODE,
     &                  ' MY_RANK= ',MY_RANK
          STOP
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
! --- [dummy arguments]
!
      integer,intent(in)  ::  ITIME
!
! --- [local entities]
!
      integer  :: start=0,end
      integer  :: NDT,NHOUR,NMIN,NSECEND
      INTEGER  :: val(8)
      integer :: TIME
!
      CALL DATE_AND_TIME(VALUES=val)
      IF(ITIME.eq.0) then
         start=TIME()
      else if(ITIME.eq.999) then
         end=TIME()
         if(my_rank.eq.ROOT) then 
           NDT=(end-start)
           NHOUR=NDT/3600
           NMIN=MOD(NDT,3600)/60
           NSECEND=MOD((MOD(NDT,3600)),60)
           write(ifll,3020) NHOUR,NMIN,NSECEND,
     &              val(1),val(2),val(3),val(5),val(6),val(7)
         endif
      else
         end=TIME()
         if(my_rank.eq.ROOT) then 
           NDT=(end-start)
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
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine SOLVER_SEND_RECV(ndim,MXNP,NP,X)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil
!
! --- [dummy arguments]
!
      integer,intent(in)   :: MXNP,NP,ndim
      real*8, intent(inout):: X(MXNP)
!
      return
!
      end subroutine solver_send_recv
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine SOLVER_SEND_RECV0(ndim,ini,MXNP,NP,X)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      use module_hpcutil
!
! --- [dummy arguments]
!
      integer,intent(in)   :: MXNP,NP,ndim,ini
      real*8, intent(inout):: X(MXNP)
!
      return
!
      end subroutine solver_send_recv0
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcimax(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      integer,intent(inout)  :: rin
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
      integer,intent(inout)  :: rin
      return
      end subroutine hpcimax_0
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrmax(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      real*8 ,intent(inout)  :: rin
!
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
      real*8 ,intent(inout)  :: rin
      return
      end subroutine hpcrmax_0
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcimin(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      integer,intent(inout)  :: rin
      return
      end subroutine hpcimin
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcimin_0(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!icpu,RIN,NUM)
! --- [dummy arguments]
!
      integer,intent(inout)  :: rin
      return
      end subroutine hpcimin_0
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrmin(rin)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      real*8 ,intent(inout)  :: rin
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
      real*8 ,intent(inout)  :: rin
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
      integer,intent(inout)  :: rin
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
      integer,intent(inout)  :: rin
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
      real*8 ,intent(inout)  :: rin
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
      real*8 ,intent(inout)  :: rin
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
      integer,intent(in) :: idim
      real*8,intent(in)  :: rin(idim)

!
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
      integer,intent(in) :: idim,idst
      real*8,intent(in)  :: rin(idim)
!
      return
      end subroutine hpcrasum_0
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcibcast(icpu,RIN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      integer,intent(inout)  :: icpu,RIN
      return
      end subroutine hpcibcast
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrbcast(icpu,RIN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      use module_hpcutil
!
! --- [dummy arguments]
!
      integer,intent(inout)  :: icpu
      real*8 ,intent(inout)  :: RIN
!
      return
      end subroutine hpcrbcast
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
!
      return
      end subroutine hpcrbcast_array



!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrbcast_ary(icpu,NDIM,N,WORK)
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
      return
      end subroutine hpcrbcast_ary
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcland(RIN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- 'logical and' for all cpu
!
      use module_hpcutil
!
! --- [dummy arguments]
!
      logical,intent(inout)   :: RIN
      return
      end subroutine hpcland
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
      return
      end subroutine hpcmaxloc
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpclor(RIN)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- 'logical and' for all cpu
!
      use module_hpcutil
!
! --- [dummy arguments]
!
      logical,intent(inout)   :: RIN
      return
      end subroutine hpclor
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine read_griddata_hpc()
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      return
      end subroutine read_griddata_hpc
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine comm
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      return
      end subroutine comm
!
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_igather(ncpu,send_data,recv_buff)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      
      integer,intent(in) :: ncpu
      integer,intent(in) :: send_data
      integer,intent(inout) :: recv_buff(1:ncpu)
      integer :: ierr
      
      return
      end subroutine hpc_igather

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_dgather(ncpu,send_data,recv_buff)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      
      integer,intent(in) :: ncpu
      real(8),intent(in) :: send_data
      real(8),intent(inout) :: recv_buff(1:ncpu)
      integer :: ierr
      
      return
      end subroutine hpc_dgather
      
      
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_dsend(ndim,send_buff,dest,tagID)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [dummy arguments]
!
      implicit none
!
      integer,intent(in) :: ndim,dest,tagID
      real(8),intent(in) :: send_buff(1:ndim)
!      
      end subroutine hpc_dsend

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpc_drecv(ndim,recv_buff,source,tagID)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [dummy arguments]
!
      implicit none
!
      integer,intent(in) :: ndim,source,tagID
      real(8),intent(inout) :: recv_buff(1:ndim)

      end subroutine hpc_drecv

!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- Dummy Particle
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcbarrier()
      end subroutine hpcbarrier
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PARTICLE_MPI()
      end SUBROUTINE PARTICLE_MPI
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine PARTICLE_SEND_RECV02()
      end subroutine PARTICLE_SEND_RECV02
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine PARTICLE_SEND_RECV01()
      end subroutine PARTICLE_SEND_RECV01
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MERGEOUT_PTRESULT_MPI()
      end SUBROUTINE MERGEOUT_PTRESULT_MPI
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hpcrbcast_stream()
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      end subroutine hpcrbcast_stream
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
! --- 
!
      
! 
      W1K(1:N)=INT(ARR(2,1:N))
!
      END SUBROUTINE hpc_MAXLOC 
