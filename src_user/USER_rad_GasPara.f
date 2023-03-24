!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE User_rad_GasPara(gasmodel,ParaFlag,almn,blmn,
     &	                       awsgg,bwsgg,cwsgg)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!       gasmodel : For Real Gas Model ('SLW' or 'WSGG')
!       ParaFlag : =1 : user specifies 'almn, blmn' 
!                  =0 : Current user subroutine is NOT used 
!       almn    : coefficient of model's expression (for 'SLW')
!       blmn    : coefficient of model's expression (for 'SLW')
!       awsgg   : coefficient of model's expression (for 'WSGG')
!       bwsgg   : coefficient of model's expression (for 'WSGG')
!       cwsgg   : coefficient of model's expression (for 'WSGG')
!
!
!
!
        IMPLICIT NONE 
!
! --- [dummy arguments]
!
        REAL*8,intent(in)       :: almn(0:3,0:3,0:3,40)
        REAL*8,intent(in)       :: blmn(0:3,0:3,0:3,40)
	CHARACTER*10,intent(in) :: gasmodel
	INTEGER,intent(out)     :: ParaFlag(0:40)
!
! --- for slw
!	REAL*8,intent(out) :: almn(0:3,0:3,0:3,40),blmn(0:3,0:3,0:3,40)
     
! --- for wsgg
	REAL*8,intent(out) :: awsgg(4,6), bwsgg(4,4,5), cwsgg(4,4,3,5)
!
!--------------------------------------------------------------------------
!		- 2 -   WSGG model
!--------------------------------------------------------------------------
!
	ParaFlag     = 0
	ParaFlag(1)  = 1	!H2O,CO2
	ParaFlag(2)  = 1
	ParaFlag(5)  = 1	!CO
	ParaFlag(6)  = 2	!CH4
	ParaFlag(27) = 2	!C2H6
!awsgg	----	absorption coefficient
	awsgg      = 0.d0
	awsgg(1,1) = 0.d0
	awsgg(2,1) = 0.69d0
	awsgg(3,1) = 7.4d0
	awsgg(4,1) = 80.d0
	awsgg(1,2) = 3.85d0
	awsgg(2:4,2) = 0.d0
!bwsgg	----	weight
	bwsgg       =0.d0
	bwsgg(1,1,1)=0.364d0
	bwsgg(2,1,1)=4.74d-5
!	
	bwsgg(1,2,1)=0.266d0
	bwsgg(2,2,1)=7.19d-5
!	
	bwsgg(1,3,1)=0.252d0
	bwsgg(2,3,1)=-7.41d-5
!
	bwsgg(1,4,1)=0.118d0
	bwsgg(2,4,1)=-4.52d-5
!
	RETURN
!--------------------------------------------
!----->
!
	END SUBROUTINE User_rad_GasPara
