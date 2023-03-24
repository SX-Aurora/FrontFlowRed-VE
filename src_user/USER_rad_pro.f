!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE user_rad_pro
     &  (ICV,ncomp,spcno,xyz,t,p,rho,yys,gab,pab,sca)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!       ICV :  for radiation properity setting
!       ncomp : species component number 
!       xyz   : CV center coordinate
!       t     : Gas tempturature
!       p     : pressure
!       rho   : gas density
!       yys   : Mas fraction for species
!       gab   : gas absorbtion coefficient
!       pab   : particle absorbtion coefficient
!       sca   : particle scattering coefficient 
!
!
!
	IMPLICIT none
!
! --- [dummy arguments]
!
        integer,intent(in)   :: ICV,ncomp
        integer,intent(in)   :: spcno(ncomp)
        REAL*8 ,intent(in)   :: yys(ncomp),xyz(3),p,rho,t
        REAL*8 ,intent(out)  :: gab,pab,sca
!
	REAL*8	rr,ff,tc
!	INTEGER	spcno(ncomp),id
!	
	
	RETURN

!----->

	END SUBROUTINE user_rad_pro

