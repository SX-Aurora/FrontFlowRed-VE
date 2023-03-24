! *** *****************************************************************
      subroutine append_resultstar(fname_root,ifl,
     &     gridOnly,gridFormat,status)
! *** *****************************************************************
      use FFRdata, only : nvrtx,ncell,NBOUND,SFBOUN,NFBOUN,
     &                    NBFS,IBFACE,IFFACE,
     &                    lvcell,lacell,cord,kmesh,totdata,
     &                    numVars,nameVars,ncompFFR,ieul2ph,nrnsxFFR
      use FFRreaddata, only : iuvw_ave_rms_rex,iuvwt_rex,
     &                        ip_avex,it_avex,ip_rmsx,it_rmsx,
     &                        icomp_avex,icomp_rmsx
!      use ProgData
      implicit none
!
      character(80),intent(in) :: fname_root
      integer,intent(in) :: ifl
      logical,intent(in) :: gridOnly
      character(2),intent(in)  :: gridFormat
      integer,intent(out) :: status
      integer :: i,j,k,IV,IC,NF,ios
      integer :: len,cnt1,cnt2,cntv
      character(80) :: fname
      character(32) :: scalarExt,scalarTemp
      character(32) :: vectorExt
      logical :: vectorflag
      integer :: kclv,prevNode,MatchNode(1:8),NMAT
      integer,allocatable :: MAT_NO(:)
!
      end subroutine append_resultstar
