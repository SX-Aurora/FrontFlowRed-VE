!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE USER_UserDefine_region(IMAT,ICVL,xyz,regionNo,Cd)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      IMPLICIT none
!
! --- [dummy arguments]
!
      integer,intent(in)    :: IMAT,ICVL
      REAL*8 ,intent(in)    :: xyz(3)
      REAL*8 ,intent(out)   :: Cd
      integer,intent(out)   :: regionNo
!
! --- local arguments
!
!      REAL*8,parameter :: xmin=1.2d0,xmax=2.2d0,ymin=-0.2d0,ymax=0.2d0
      REAL*8,parameter :: xmin=1.5d0,xmax=1.6d0,ymin=-0.05d0,ymax=0.05d0
!      REAL*8,parameter :: xmin=-1.d0,xmax=1.0d0,ymin=-1.d0,ymax=1.0d0,
!     &                    zmin=-0.001d0,zmax=1.0d0
!,zmin,zmax
!
      REAL*8 :: dum1,dum2,dum3
!
!##################################################################
! Explain:
!     userflag(2)=1
!     regionNo=1: KLES region : 1.5 order SGS model
!             =2: SLES region : Smogorinsky model
!             =3: DLES region : Dynamic-SGS 
! 
!     userflag(3)=1
!     regionNo=21~50 : canopy region 
!                      Cd [Drag Force coeff.]
!##################################################################
!------------------------------------------------
! --- multi-region for different LES model
!------------------------------------------------
!      if(xyz(1)>xmin.and.xyz(1)<xmax.and.xyz(2)>ymin.and.xyz(2)<ymax) 
!     &  then
!!        regionNo=2     ! Smogorinsky model
!        regionNo=3     ! Dynamic Smogorinsky model
!      else
!        regionNo=1     ! 1.5 order SGS model
!      endif
!------------------------
! --- canopy region
!------------------------
      if(xyz(1)>xmin.and.xyz(1)<xmax.and.xyz(2)>ymin.and.xyz(2)<ymax) 
     &  then
        regionNo=21 
        Cd=1000.d0
      endif
!------------------------
! --- canopy region
!------------------------
!      if(xyz(1)>xmin.and.xyz(1)<xmax.and.xyz(2)>ymin.and.xyz(2)<ymax.and.
!     &   xyz(3)>zmin.and.xyz(3)<zmax)then
!        regionNo=21 
!        Cd=100.d0
!      endif
!
      return
!
      end SUBROUTINE USER_UserDefine_region

        
