! *** *****************************************************************
      subroutine append_resultscry(fname,ifl,gridOnly,status)
! *** *****************************************************************
      use SCRYUheaderTag, only : L2D3D,LFORT
      use SCRYUdata, only : LVCT,NDATA,IRECL,LNX,LCORD,NBNN,NTTS,
     &                      LKAI,NTTE,NRTY,TITLE,NGFAX,LRGN,LNAM,IRTY,
     &                      OverlapStart_n,OverlapEnd,
     &                      IPTYP,IPMAT,sumnod
      use SCRYUpre, only : NNODS,NELEM,NREGN,NDE,NDNO,
     &                     IEtot,fnode,FTYP,IETYP,NKND
      use FFRdata, only : nvrtx,ncell,NBOUND,SFBOUN,NFBOUN,
     &                    NBFS,IBFACE,IFFACE,
     &                    lvcell,lacell,cord,kmesh,totdata,
     &                    numVars,nameVars,ncompFFR,ieul2ph,nrnsxFFR,
     &                    npotnxFFR,NFLIDxFFR
      use FFRreaddata, only : iuvw_ave_rms_rex,iuvwt_rex,
     &                        ip_avex,it_avex,ip_rmsx,it_rmsx,
     &                        icomp_avex,icomp_rmsx
!      use ProgData
      implicit none
!
      character*80,intent(in) :: fname
      integer,intent(in) :: ifl
      logical,intent(in) :: gridOnly
      integer,intent(out) :: status
      integer :: i,j,k,ii,jj,ios
      integer :: skp0
      integer :: skp1,skp2,skp3,skp4,skp5,skp6,skp7,skp8,skp9
!
      integer,allocatable :: WNDNO(:)
!
! ---
      return
      end subroutine append_resultscry
