	
C======================================================================C
C                                                                      C
C Subroutine COSVALUE()						       C
C                                                                      C
C                                       2003/12/01                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

!-----------------------------------------------------------------------------
!	Subroutine List
!							/////// -- usage
!     subroutine	DistPorj()
!     subroutine	FaceEqu()
!     subroutine	FaceAreaVect()
!     subroutine	PloymFaceArea()
!     subroutine	CosValue()
!     subroutine	SinValue()
!     subroutine	IntersectNode()
!     subroutine	GetGridIndex()
!     subroutine	ZuobiaoRot()
!     subroutine	RotVector()
!     subroutine	GetTanLine()
!     subroutine	CheckRightHandTurn()
!     subroutine	WallFaceEquVect()
!     subroutine	ElemVolume()
!     subroutine	GetSolidAngle()
!     subroutine	PointInGlobe()
!     subroutine	GetDistination()
!     subroutine	GetDistance()
!     subroutine	BuildZuoBiao()
!     subroutine	AreaVector()
!     subroutine	AreaVector2()
!     subroutine	FaceAreaEqu()
!     subroutine	SurfaceExtract()
!     subroutine	GetToken()
!     subroutine	GetArraySize()
!     subroutine	SymmTreatRD()
!     subroutine	SmoothLTRI()
!     subroutine	SmoothLTRI2()
!     subroutine	SmoothBand()
!     subroutine	WallGroupEASum()
!     subroutine	GetSegment()
!     subroutine	GetNodeNeb()
!     subroutine	GetPrdPair()
!										/////// -- shape function
!	subroutine  quadri4shapefun()
!	subroutine  tet4shapefun()
!	subroutine	tri3shapefun()
!	subroutine  pyramid5shapefun()
!	surrouinte	prism6shapefun()
!	subroutine	hex8shapefun()
!	subroutine  qua4pos()
!	subroutine	tri3pos()
!	subroutine  Tet4Pos()
!	subroutine  Pyramid5Pos()
!	subroutine  Prism6Pos()
!	subroutine  Hex8Pos()
!	subroutine  Pyramid6Arrange()
!	subroutine  Prism6Arrange()
!	subroutine  Hex8Arrange()
!	subroutine  TetVolume()
!	subroutine  PyramidVolume()
!	subroutine  GaussTri3()
!	subroutine  GaussSibian()
!	subroutine  GaussFaceQita
!	subroutine  GaussTet4()
!	subroutine  GaussHex8()
!	subroutine  GaussVolQita()
!	subroutine	isinface()
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE DistProj(A,B,C,D,X,Y,Z,XNB,YNB,ZNB,tyx,tyy,tyz,dist)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT NONE  
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
C This subroutine is used to get the distance from a point to a face
C and the projection in the face
C face			 : Ax+By+Cz+D = 0
C Point			 : P(x,y,z)
C Projection point : (tyx,tyy,tyz)
C Distance		 : Dist
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      real*8,intent(INOUT) :: A,B,C,D,X,Y,Z,
     &                        XNB,YNB,ZNB,tyx,tyy,tyz,dist
!
      real*8 :: x1,y1,z1,x2,y2,z2,r1,r2
!
	dist = DABS(A*X+B*Y+C*Z+D)/DSQRT(A*A+B*B+C*C)
!
	x1 = X-dist*XNB
	y1 = Y-dist*YNB
	z1 = Z-dist*ZNB

	x2 = X+dist*XNB
	y2 = Y+dist*YNB
	z2 = Z+dist*ZNB

	r1 = DABS(A*x1+B*y1+C*z1+D)
	r2 = DABS(A*x2+B*y2+C*z2+D)
	
	IF(r1.GT.r2) THEN
		tyx = x2
		tyy = y2
		tyz = z2
	ELSE
		tyx = x1
		tyy = y1
		tyz = z1
		!WRITE(*,*) 'r1=', r1
	END IF
	
	!WRITE(IUT0,*) 'r=', r1,r2,xnb,ynb,znb
	

	return

	END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FaceEqu(ABCD,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT NONE
      real*8,intent(INOUT) :: ABCD(4)
      real*8,intent(INOUT) :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
! --- 
      real*8 :: Vector(3),vm

!	
		
C This subroutine is used to get the face equation in form:
C Ax+By+Cz+D = 0
C And face normal vector	
!	

      ABCD(1) = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
      ABCD(2) = (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1)
      ABCD(3) = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      ABCD(4) = -(ABCD(1)*x1+ABCD(2)*y1+ABCD(3)*z1)
	
      vm = SQRT(ABCD(1)*ABCD(1) + ABCD(2)*ABCD(2) + ABCD(3)*ABCD(3))
      ABCD(1) = ABCD(1)/vm
      ABCD(2) = ABCD(2)/vm
      ABCD(3) = ABCD(3)/vm
      ABCD(4) = ABCD(4)/vm
      RETURN
      END 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE FaceAreaVect(Vector,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,A)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      
      IMPLICIT NONE 
!
      real*8,intent(INOUT) :: Vector(3)
      real*8,intent(INOUT) :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,A
!
      real*8 :: ABCD(4),ax,ay,az,bx,by,bz,xnb,ynb,znb
!
!
      AX  = X2-X1
      AY  = Y2-Y1
      AZ  = Z2-Z1
      BX  = X3-X1
      BY  = Y3-Y1
      BZ  = Z3-Z1

      XNB = AY*BZ-AZ*BY
      YNB = AZ*BX-AX*BZ
      ZNB = AX*BY-AY*BX
      A = (XNB*XNB+YNB*YNB+ZNB*ZNB)
      A = 0.5*DSQRT(A)

      IF(A.LE.1.0E-16) THEN
	vector(1) = 0.0
	vector(2) = 0.0
	vector(3) = 0.0
      ELSE
	vector(1)  = 0.5*XNB/A
	vector(2)  = 0.5*YNB/A
	vector(3)  = 0.5*ZNB/A
      END IF
      RETURN
      END  SUBROUTINE FaceAreaVect
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PolymFaceArea(XYZ,M,N,Area)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT NONE    
      INTEGER,intent(INOUT) :: M,N
      real*8, intent(INOUT) :: XYZ(3,M),Area
!
      INTEGER :: i,i2,i3
      real*8  :: ax,ay,az,bx,by,bz,xnb,ynb,znb,a

      area = 0
      DO I = 1, N-2
	I2 = I+1	!- ((I+1)/N)*N
	I3 = I+2	!- ((I+2)/N)*N
	AX  = XYZ(1,I2)-XYZ(1,1)
	AY  = XYZ(2,I2)-XYZ(2,1)
	AZ  = XYZ(3,I2)-XYZ(3,1)
	BX  = XYZ(1,I3)-XYZ(1,1)
	BY  = XYZ(2,I3)-XYZ(2,1)
	BZ  = XYZ(3,I3)-XYZ(3,1)

	XNB = AY*BZ-AZ*BY
	YNB = AZ*BX-AX*BZ
	ZNB = AX*BY-AY*BX
	A = (XNB*XNB+YNB*YNB+ZNB*ZNB)
	A = 0.5*DSQRT(A)
	area = area + A
	
	END DO	
		
	RETURN
	END 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE COSVALUE(X1,Y1,Z1,X2,Y2,Z2,value)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	IMPLICIT NONE  !REAL*8(A-H,O-Z)
!
	real*8,intent(INOUT) :: X1,Y1,Z1,X2,Y2,Z2,value
	real*8 :: d1,d2
!
C	This subroutine is used to calculate the cosine of the angle
C	of the two lines
	
	d1 = SQRT(x1*x1 + y1*y1 + z1*z1)
	d2 = SQRT(x2*x2 + y2*y2 + z2*z2)
	
	value = (x1*x2 + y1*y2 +z1*z2)
	IF(d1*d2.GT.1.0E-16) THEN
	  value = value/(d1*d2)
	ELSE
	  value = 0.0
	END IF
!
	RETURN
!
	END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE SINVALUE(ax,ay,az,bx,by,bz,value)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	IMPLICIT NONE !REAL*8(A-H,O-Z)
	real*8,intent(INOUT) :: ax,ay,az,bx,by,bz,value
!
	real*8 :: d1,d2,v1,v2,v3
C	This subroutine is used to calculate the sine of the angle
C	of the two lines
	
	d1 = SQRT(ax*ax + ay*ay + az*az)
	d2 = SQRT(bx*bx + by*by + bz*bz)
	
	v1 = ay*bz - az*by
	v2 = az*bx - ax*bz
	v3 = ax*by - ay*bx
	value = SQRT(v1*v1 + v2*v2 + v3*v3)

	IF(d1*d2.GT.1.0E-12) THEN
          value = value/(d1*d2)
	ELSE
	  value = 0.0
	END IF

	RETURN


	END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE  IntersectNode
     &    (A,B,C,D,X0,Y0,Z0,Dist,Eta,Seta,XW,YW,ZW,
     &	   idirect,ierr)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT NONE
      real*8,intent(INOUT) :: A,B,C,D,X0,Y0,Z0,Dist,Eta,Seta,XW,YW,ZW
      INTEGER,intent(INOUT) :: idirect,ierr
!
      real*8 :: r2



	!DIRECT = 1	:	Forward direction
	!DIRECT =-1 :	Backward direction
		
c	X1 = X0 - sin(Eta)*cos(Seta)*Dist
c	Y1 = Y0 - sin(Eta)*sin(Seta)*Dist
c	Z1 = Z0 - cos(Eta)*Dist

	XW = X0 + sin(Eta)*cos(Seta)*Dist*idirect
	YW = Y0 + sin(Eta)*sin(Seta)*Dist*idirect
	ZW = Z0 + cos(Eta)*Dist*idirect

c	r1 = DABS(A*x1+B*y1+C*z1+D)
	r2 = DABS(A*XW+B*YW+C*ZW+D)/(A*A+B*B+C*C)

	
	IF(r2.GT.1.0E-4) THEN
          ierr = 1
	ELSE
          ierr = 0
	END IF

	RETURN

	END 

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE  GetGridIndex(N,Grid,Point,I)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
	IMPLICIT NONE
	INTEGER,intent(INOUT) :: N
	real*8,intent(INOUT) :: Grid(N),Point
	INTEGER :: I
	
        I = 1
	DO 10 WHILE (Point.GT.Grid(I+1).AND.I.LT.N)
        I = I + 1
 10     CONTINUE

	RETURN

	END



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE  
     &  ZuobiaoRot(Seta,Eta,NX,NY,NZ,TX1,TY1,TZ1,TX2,TY2,TZ2)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.14159265)
	REAL*8 Seta,Eta,NX,NY,NZ,TX1,TY1,TZ1,TX2,TY2,TZ2
!direction cosine
!new OZ, normal direction
	unit = DSQRT(NX*NX+NY*NY+NZ*NZ)
	zcosx = NX / unit
	zcosy = NY / unit
	zcosz = NZ / unit
			
	
!new OX, tangent line in plane, 1 
	unit = DSQRT(TX1*TX1+TY1*TY1+TZ1*TZ1)
	xcosx = TX1 / unit
	xcosy = TY1 / unit
	xcosz = TZ1 / unit


!new OY, tangent line in plane, 2 
	unit = DSQRT(TX2*TX2+TY2*TY2+TZ2*TZ2)
	ycosx = TX2 / unit
	ycosy = TY2 / unit
	ycosz = TZ2 / unit

	xnew = sin(Eta)*cos(Seta)
	ynew = sin(Eta)*sin(Seta)
	znew = cos(Eta)
	
	xori = xcosx*xnew + ycosx*ynew + zcosx*znew
	yori = xcosy*xnew + ycosy*ynew + zcosy*znew
	zori = xcosz*xnew + ycosz*ynew + zcosz*znew

	a = DSQRT(xori*xori+yori*yori+zori*zori)
	xori = xori/a
	yori = yori/a
	zori = zori/a

	Eta = acos(zori)

	a = DSQRT(xori*xori+yori*yori)
	IF(a.GT.0) THEN
		Seta = acos(xori/a)
	ELSE
		Seta = 0
	END IF
	IF(yori.lt.0.0) Seta = -Seta
	RETURN

	END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE  WallRelativeVect
     &    (V,NX,NY,NZ,TX1,TY1,TZ1,TX2,TY2,TZ2)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- 
     	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.14159265)
	REAL*8 Seta,Eta,NX,NY,NZ,TX1,TY1,TZ1,TX2,TY2,TZ2
	REAL*8 V(3)

!
!This is to get the vector expressed with respect to the wall normal
!direction
!new OZ, normal direction
!
	unit = DSQRT(NX*NX+NY*NY+NZ*NZ)
	zcosx = NX / unit
	zcosy = NY / unit
	zcosz = NZ / unit
!new OX, tangent line in plane, 1 
	unit = DSQRT(TX1*TX1+TY1*TY1+TZ1*TZ1)
	xcosx = TX1 / unit
	xcosy = TY1 / unit
	xcosz = TZ1 / unit
!new OY, tangent line in plane, 2 
	unit = DSQRT(TX2*TX2+TY2*TY2+TZ2*TZ2)
	ycosx = TX2 / unit
	ycosy = TY2 / unit
	ycosz = TZ2 / unit
	xori = xcosx*V(1) + ycosx*V(2) + zcosx*V(3)
	yori = xcosy*V(1) + ycosy*V(2) + zcosy*V(3)
	zori = xcosz*V(1) + ycosz*V(2) + zcosz*V(3)
	a = DSQRT(xori*xori+yori*yori+zori*zori)
	V(1) = xori/a
	V(2) = yori/a
	V(3) = zori/a
	RETURN

	END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE RotVector(AX,AY,AZ,BX,BY,BZ,vx,vy,vz,A)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT REAL*8(A-H,O-Z)
!This subroutine is used to vector product (rotation) of the two lines
!0->1 and 0->2
! --- 
      XNB = AY*BZ-AZ*BY
      YNB = AZ*BX-AX*BZ
      ZNB = AX*BY-AY*BX
      A = (XNB*XNB+YNB*YNB+ZNB*ZNB)
      A = 0.5*DSQRT(A)
      IF(A.LE.1.0E-12) THEN
	vx = 0.0
	vy = 0.0
	vz = 0.0
      ELSE
	vx  = 0.5*XNB/A
	vy  = 0.5*YNB/A
	vz  = 0.5*ZNB/A
      END IF
      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GetTanLine(WABCD,WNV,WTV1,WTV2,M,NDOF)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 :: WABCD(4,M),WNV(NDOF,M),WTV1(NDOF,M),WTV2(NDOF,M)
      REAL*8 :: P1(3),P2(3)


      DO 100 I = 1,M
!WTV1
	J = 1
      DO 50 WHILE(WABCD(J,I).EQ.0)
	J = J+1
 50   CONTINUE

      a = WABCD(4,I)
      DO 60 K = 1,3
	IF(J.NE.K) THEN
	  P1(K) = 1.0
	  a = a + WABCD(K,I)*P1(K)
      END IF
	
 60   CONTINUE
      P1(J) = a/WABCD(J,I)
		
      a = WABCD(4,I)
      DO 65 K = 1,3
	IF(J.NE.K) THEN
	  P2(K) = 2.0
          a=a+WABCD(K,I)*P2(K)
        END IF
	
 65   CONTINUE
      P2(J) = a/WABCD(J,I)
		
	
      DO 70 J = 1,3
	WTV1(J,I) = P1(J)-P2(J)
 70   CONTINUE	
	
      a = SQRT(WTV1(1,I)*WTV1(1,I)+WTV1(2,I)*WTV1(2,I)+
     &	WTV1(3,I)*WTV1(3,I))
      DO 80 J = 1,3
	WTV1(J,I) = WTV1(J,I)/a
 80   CONTINUE
!WTV2
      CALL RotVector(WNV(1,I),WNV(2,I),WNV(3,I),
     &	   WTV1(1,I),WTV1(2,I),WTV1(3,I),	
     &	   WTV2(1,I),WTV2(2,I),WTV2(3,I),b)
 100  CONTINUE

      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CheckRightHandTurn(IUT0,CTR,FaceNode,RHTIsOutward,N2D)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT REAL*8(A-H,O-Z)
		
      REAL*8 :: CTR(3),FaceNode(N2D,3),NV(3),FaceCTR(3),TV(3)
      CHARACTER*20 :: RHTIsOutward
	

!Check Whether the Rotate-Vector of Right_Hand_Turn of Vetices is ourward

      FaceCTR = 0
      DO 20 IDOF = 1,3
      DO 10 J = 1,3
	FaceCTR(IDOF) = FaceCTR(IDOF)+FaceNode(J,IDOF)
 10   CONTINUE
	FaceCTR(IDOF) = FaceCTR(IDOF)/3.0
	TV(IDOF) = FaceCTR(IDOF)-CTR(IDOF)
 20   CONTINUE

      AX = FaceNode(2,1)-FaceNode(1,1)
	AY = FaceNode(2,2)-FaceNode(1,2)
	AZ = FaceNode(2,3)-FaceNode(1,3)
	BX = FaceNode(3,1)-FaceNode(1,1)
	BY = FaceNode(3,2)-FaceNode(1,2)
	BZ = FaceNode(3,3)-FaceNode(1,3)
	CALL RotVector(AX,AY,AZ,BX,BY,BZ,NV(1),NV(2),NV(3),A)

	CALL COSVALUE(NV(1),NV(2),NV(3),TV(1),TV(2),TV(3),value)
      IF(value.GT.0) THEN
	RHTIsOutward = 'Yes'
	WRITE(IUT0,*) '      Connectivity Is anti-clockwise'
      ELSE
	RHTIsOutward = 'No'
	WRITE(IUT0,*) '      Connectivity Is clockwise'
      END IF
      RETURN
      END 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE WallFaceEquVect(WABCD,WNV,WElem,WArea,XYZ,ELEM,WNEB,
     &			ElemID,NNPE,NDOF,NP,NE,NW,N2D,N3D,NKIND)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER WElem(N2D+2,NW),WNEB(NW)
      INTEGER ELEM(N3D,NE),ElemID(NE),NNPE(NKIND)
      REAL*8	WABCD(4,NW),WNV(NDOF,NW),WArea(NW)
      REAL*8	XYZ(NDOF,NP)
      REAL*8	ABCD(4),NV(3),GCTR(3),FCTR(3),TV(3)
      CHARACTER*20 RHTIsOutward
!get normal vector and face area
      DO 100 I = 1,NW
!normal vector
        IP1 = WElem(2,I)
	IP2 = WElem(3,I)
	IP3 = WElem(4,I)

	CALL FaceEqu(ABCD,XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &					  XYZ(1,IP2),XYZ(2,IP2),XYZ(3,IP2),
     &					  XYZ(1,IP3),XYZ(2,IP3),XYZ(3,IP3))
	DO 50 J = 1,4
	  WABCD(J,I) = ABCD(J)
 50	CONTINUE
	CALL FaceAreaVect(NV,XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &						 XYZ(1,IP2),XYZ(2,IP2),XYZ(3,IP2),
     &						 XYZ(1,IP3),XYZ(2,IP3),XYZ(3,IP3),Area)	
	DO 60 J = 1,3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WNV should be the OUTWARD normal direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	WNV(J,I) = NV(J)
			
 60	CONTINUE
	WArea(I) = Area
	DO 70 J = 2,WElem(1,I)-2
          IP1 = WElem(J+2,I)
	  IP2 = WElem(J+3,I)
	  IP3 = WElem(2,I)
          CALL FaceAreaVect(NV,XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &	    XYZ(1,IP2),XYZ(2,IP2),XYZ(3,IP2),
     &	    XYZ(1,IP3),XYZ(2,IP3),XYZ(3,IP3),Area)
          WArea(I) = WArea(I) + Area
 70     CONTINUE	
 100	CONTINUE
!make the normal vectors direct outward
	DO 200 IW = 1,NW
          FCTR = 0.0
	  DO 150 J=1,WElem(1,IW)
	  IP = WElem(J+1,IW)
	  FCTR(:) = FCTR(:)+XYZ(:,IP)
 150      CONTINUE
	  FCTR(:) = FCTR(:)/WElem(1,IW)
          IE = WNEB(IW)
	  INNPE = NNPE(ElemID(IE))
	  GCTR(:) = 0.0
	  DO 10 I = 1,INNPE
	    IP = ELEM(I,IE)
	    GCTR(:) = GCTR(:) + XYZ(:,IP)
10	  CONTINUE
	  GCTR(:) = GCTR(:)/INNPE
	  TV(:) = FCTR(:)-GCTR(:)
	  CALL COSVALUE
     &    (WNV(1,IW),WNV(2,IW),WNV(3,IW),TV(1),TV(2),TV(3),
     &				  value)

          IF(value.LT.0) THEN
!inward
            WNV(:,IW) = -WNV(:,IW)
	  END IF
200	CONTINUE

	RETURN

!--------> end of WallFaceEquVect()

	END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ElemVolume(XYZ,ELEM,V,N2D,N3D,NDOF,NFACE,NP,NE,MP,ME,
     &		  NKIND,ElemID,FaceNum,NNPE,NNPS,Local,AveV)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT REAL*8(A-H,O-Z)

      REAL*8  XYZ(NDOF,MP),V(NE)
	INTEGER ELEM(N3D,ME),NNPE(NKIND)
	INTEGER ElemID(NE), FaceNum(NKIND),NNPS(NKIND,NFACE),
     &		Local(NKIND,NFACE,N2D)
	
	REAL*8 CTR(3),P1(3),P2(3),P3(3)

C
C
C      CALCULATE ELEMENT'S VOLUME
C         ( 3-D ; SINGLE PRECISION, REFERING ELEMENT BY TYPE LIST )
C
C
C     ARGUMENT LISTINGS
C       (1) INPUT

C          XYZ    (NDOF,IP) ; X,Y,Z-DIR. COORDINATE  OF NODE
C          ELEM (I,IE) ; NODE TABLE,
C          NE          ; NUMBER OF TOTAL ELEMENTS
C          NP          ; NUMBER OF TOTAL    NODES
C          N3D           ; NUMBER OF NODES ASSIGNED TO ONE ELEMENT
C
C       (2) OUTPUT
C          NV    (NDOF,IE) ; X,Y,Z COMPONENT OF THE NORMAL VECTOR OF SURFACE IE
C          A     (IE) ; AREA OF SURFACE IE
C
C
	
	
	AveV = 0
	DO 200 IE = 1 , NE
	  V(IE) = 0

	  DO 20 IDOF = 1,NDOF
		CTR(IDOF) = 0.0
		DO 10 J = 1, NNPE(ElemID(IE))
		  IP = ELEM(J,IE)
		  CTR(IDOF) = CTR(IDOF) + XYZ(IDOF,IP)
10		CONTINUE
		CTR(IDOF) = CTR(IDOF)/NNPE(ElemID(IE))
20	  CONTINUE

	  ID = ElemID(IE)
	  DO 50 JFACE = 1,FaceNum(ID)
		DO 50 JNODE = 1, NNPS(ID,JFACE)-2
			
			IP1 = ELEM(Local(ID,JFACE,JNODE+1),IE)
			IP2 = ELEM(Local(ID,JFACE,JNODE+2),IE)
			IP3 = ELEM(Local(ID,JFACE,1),IE)
              
			DO 30 IDOF = 1,NDOF
				P1(IDOF) = XYZ(IDOF,IP1)
				P2(IDOF) = XYZ(IDOF,IP2)
				P3(IDOF) = XYZ(IDOF,IP3)
30			CONTINUE
		
		AX  = (CTR(2)-P2(2))*(CTR(3)-P3(3))
     &		  -(CTR(2)-P3(2))*(CTR(3)-P2(3))
		AX  = (CTR(1)-P1(1))*AX

		AY  = (CTR(1)-P2(1))*(CTR(3)-P3(3))
     &		  -(CTR(1)-P3(1))*(CTR(3)-P2(3))
		AY  = (CTR(2)-P1(2))*AY

		AZ  = (CTR(1)-P2(1))*(CTR(2)-P3(2))
     &		  -(CTR(1)-P3(1))*(CTR(2)-P2(2))
		AZ  = (CTR(3)-P1(3))*AZ

          VVV = (AX - AY + AZ)/6.0

		V(IE) = V(IE) + ABS(VVV)
50		CONTINUE
		AveV = AveV+V(IE)
  200 CONTINUE
 
	AveV = AveV/NE
      RETURN
      END

C======================================================================C
	SUBROUTINE  GetSolidAngle(Seta,Eta,VX,VY,VZ,CX,CY,CZ)
C======================================================================C
     	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.141592653589793)
	REAL*8 Seta,Eta,VX,VY,VZ,CX,CY,CZ


	!Seta:	circumferential angle
	!eta:	cone angle

	!cone angle
	
	unit = DSQRT((VX-CX)*(VX-CX)+(VY-CY)*(VY-CY)+(VZ-CZ)*(VZ-CZ))
	IF(unit.GT.0.0)  THEN
		cosz = (VZ-CZ) / unit
		eta = acos(cosz)
	ELSE
		eta = 0
	END IF

	
	!circumferential angle
	unit = DSQRT((VX-CX)*(VX-CX)+(VY-CY)*(VY-CY))
	IF(unit.GT.0.0)  THEN
		cosx = (VX-CX) / unit
		cosy = (VY-CY) / unit
		seta = acos(cosx)
		IF(cosy.lt.0) seta = 2.0*PAI-seta
	ELSE
		seta = 0.0
	END IF

	RETURN

	END


	SUBROUTINE  PointInGlobe(V0,CTR,VW,ROUT,Eta,Seta,IFLAG)

     	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.141592653589793)
	REAL*8 V0(3),CTR(3),VW(3),ROUT,Eta,Seta


	fdao = 1.0/ROUT

	
	CALL GetDistance(CTR(1),CTR(2),CTR(3),V0(1),V0(2),V0(3),r1)


	IFLAG = 1

	d1 = 2.0*ROUT
	d2 = 0
	CALL GetDistination(V0(1),V0(2),V0(3),TX1,TY1,TZ1,d1,Eta,Seta)
	CALL GetDistance(CTR(1),CTR(2),CTR(3),TX1,TY1,TZ1,r1)
	CALL GetDistination(V0(1),V0(2),V0(3),TX2,TY2,TZ2,d2,Eta,Seta)
	CALL GetDistance(CTR(1),CTR(2),CTR(3),TX2,TY2,TZ2,r2)

	IF(r2.GT.ROUT) THEN
		CALL COSVALUE(V0(1)-CTR(1),V0(2)-CTR(2),V0(3)-CTR(3),
     &				  TX1-V0(1),TY1-V0(2),TZ1-V0(3),value)
		IF(value.GT.-1.0E-2) THEN
				!V0 is out of the globe, and the ray directs outward, 
				!will not intersect any wall face
			IFLAG = 0
			RETURN
		END IF
	END IF	

	
	d3 = 0.5*(d1+d2)

	itime = 1
	
50	CALL GetDistination(V0(1),V0(2),V0(3),VW(1),VW(2),VW(3),d3,Eta,Seta)
	CALL GetDistance(CTR(1),CTR(2),CTR(3),VW(1),VW(2),VW(3),r3)
	
	itime = itime + 1
	error = abs(r3-ROUT)*fdao

	!write(*,*) 'err=',error,itime
	


	IF(error.LE.5.0E-03) GOTO 100

		IF(r2.LT.ROUT.AND.r3.LT.ROUT) THEN		!----v2----v3----P----v1
	
			r2 = max(r2,r3)
			d2 = max(d2,d3)
			d3 = (d1+d2)*0.5
	
						
		ELSE IF(r2.GT.ROUT.AND.r3.GT.ROUT) THEN	!----P----v2----v3----v1
			
			r1 = min(r2,r3)
			d1 = min(d2,d3)
			r2 = 0.0
			d2 = 0.0
			d3 = (d2+d1)*0.5
	
		ELSE
			
			d1 = max(d2,d3)
			d2 = min(d2,d3)
			r1 = max(r2,r3)
			r2 = min(r2,r3)
			d3 = 0.5*(d2+d3)

		END IF

		GOTO 50

100	CONTINUE




	RETURN


	END
C======================================================================C
	SUBROUTINE  GetDistination(SX,SY,SZ,TX,TY,TZ,D,Eta,Seta)
C======================================================================C
     	IMPLICIT REAL*8(A-H,O-Z)
	
	TX = SX + sin(Eta)*cos(Seta)*D
	TY = SY + sin(Eta)*sin(Seta)*D
	TZ = SZ + cos(Eta)*D


	RETURN
	END	
C======================================================================C
	SUBROUTINE  GetDistance(SX,SY,SZ,TX,TY,TZ,D)
C======================================================================C
     	IMPLICIT REAL*8(A-H,O-Z)
	
	D = SQRT((SX-TX)*(SX-TX)+(SY-TY)*(SY-TY)+(SZ-TZ)*(SZ-TZ))


	RETURN


	END	
C======================================================================C
	SUBROUTINE BuildZuoBiao(NV,TV1,TV2,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8	ABCD(4),NV(NDOF),TV1(NDOF),TV2(NDOF)
	REAL*8	P1(3),P2(3)


		
	!TV1
	
	unit = SQRT(NV(1)*NV(1)+NV(2)*NV(2)+NV(3)*NV(3))
	IF(unit.EQ.0.0) THEN
		WRITE(*,*) 'Error in BuildZuoBiao():'
		WRITE(*,*) 'Zero normal vector !'
		STOP
	ELSE
		NV = NV/unit
	END IF

	J = 1
	
	DO 50 WHILE(ABS(NV(J)).LE.1.0E-12)
		J = J+1
50	CONTINUE
	
	VMODE = 0.0
	DO 60 K = 1,3
		IF(J.NE.K) THEN
			TV1(K) = 1.0
			VMODE = VMODE + NV(K)*TV1(K)
		END IF
	
60	CONTINUE
	

	TV1(J) = -VMODE/NV(J)
	
	VMODE = SQRT(TV1(1)*TV1(1)+TV1(2)*TV1(2)+TV1(3)*TV1(3))
	TV1 = TV1 / VMODE

	!TV2
	CALL RotVector(NV(1),NV(2),NV(3),
     &			   TV1(1),TV1(2),TV1(3),	
     &			   TV2(1),TV2(2),TV2(3),b)
	


100	CONTINUE

	RETURN



	END







C======================================================================C
C                                                                      C
C Subroutine AreaVector()												 C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2003/12/01                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C



	SUBROUTINE AreaVector(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,vx,vy,vz,A)
	IMPLICIT REAL*8(A-H,O-Z)
	!This subroutine is used to vector product (rotation) of the two lines
	!0->1 and 0->2


	AX  = X1-X0
      AY  = Y1-Y0
	AZ  = Z1-Z0
	
	BX  = X2-X0
      BY  = Y2-Y0
	BZ  = Z2-Z0

C
      XNB = AY*BZ-AZ*BY
      YNB = AZ*BX-AX*BZ
      ZNB = AX*BY-AY*BX
      A = (XNB*XNB+YNB*YNB+ZNB*ZNB)
      A = 0.5*DSQRT(A)


	IF(A.LE.1.0E-5) THEN
		vx = 0.0
		vy = 0.0
		vz = 0.0
	ELSE
		vx  = 0.5*XNB/A
		vy  = 0.5*YNB/A
		vz  = 0.5*ZNB/A
	END IF
	
	

	RETURN
	END





C======================================================================C
C                                                                      C
C Subroutine AreaVector2()											 C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2003/12/01                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C



	SUBROUTINE AreaVector2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,vx,vy,vz,A,ratio,
     &					  iflag)
	IMPLICIT REAL*8(A-H,O-Z)
	!This subroutine is used to vector product (rotation) of the two lines
	!0->1 and 0->2


	AX  = X1-X0
      AY  = Y1-Y0
	AZ  = Z1-Z0
	
	BX  = X2-X0
      BY  = Y2-Y0
	BZ  = Z2-Z0

C
      XNB = AY*BZ-AZ*BY
      YNB = AZ*BX-AX*BZ
      ZNB = AX*BY-AY*BX
      A = (XNB*XNB+YNB*YNB+ZNB*ZNB)
      A = 0.5*DSQRT(A)
	b=A/ratio

	if(iflag.EQ.0) return

	IF(b.LE.1.0E-10) THEN
		vx = 0.0
		vy = 0.0
		vz = 0.0
	ELSE
		vx  = 0.5*XNB/A
		vy  = 0.5*YNB/A
		vz  = 0.5*ZNB/A
	END IF
	
	

	RETURN
	END






C======================================================================C
C                                                                      C
C Subroutine FaceAreaEqu()											 C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/09                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C


	SUBROUTINE FaceAreaEqu(Vector,ABCD,XYZ,NP,MP,NDOF)
	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XYZ(NDOF,MP),ABCD(4),VECTOR(3)
	
		
C	This subroutine is used to get the face equation in form:
C		Ax+By+Cz+D = 0
C	And face normal vector	
	
	!	
	x1 = XYZ(1,1)
	x2 = XYZ(1,2)
	x3 = XYZ(1,3)
	y1 = XYZ(2,1)
	y2 = XYZ(2,2)
	y3 = XYZ(2,3)
	z1 = XYZ(3,1)
	z2 = XYZ(3,2)
	z3 = XYZ(3,3)

	ABCD(1) = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
	ABCD(2) = (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1)
	ABCD(3) = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
	ABCD(4) = -(ABCD(1)*x1+ABCD(2)*y1+ABCD(3)*z1)


	!
	AX  = X2-X1
      AY  = Y2-Y1
	AZ  = Z2-Z1
	BX  = X3-X1
      BY  = Y3-Y1
	BZ  = Z3-Z1

      XNB = AY*BZ-AZ*BY
      YNB = AZ*BX-AX*BZ
      ZNB = AX*BY-AY*BX
      A = (XNB*XNB+YNB*YNB+ZNB*ZNB)
      A = 0.5*DSQRT(A)
	IF(A.LE.1.0E-16) THEN
		vector(1) = 0.0
		vector(2) = 0.0
		vector(3) = 0.0
	ELSE
		vector(1)  = 0.5*XNB/A
		vector(2)  = 0.5*YNB/A
		vector(3)  = 0.5*ZNB/A
	END IF
	
	
	RETURN


	END 









C======================================================================C
C                                                                      C
C Subroutine SurfaceExtract()									         C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/06/09                     C
C                                       2006/03/25                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

	SUBROUTINE	SurfaceExtract(IUT0,XYZ,ELEM,
     &				ElemID,FaceNum,NNPS,Local,NNPE,WPCH,
     &				NP,NE,NW,NWT,MP,ME,MW,NDOF,N2D,N3D,NFACE,
     &				MK,NEB,NBOUN)


	use module_radsxf, only : NodeBounId

	USE MODULE_RADWORK, ONLY : WK1,WK2,WK3,WK4,WK5
	USE MODULE_RADWORK, ONLY : VWK3D => VWKM1
	USE MODULE_RADWORK, ONLY : PATCH => WKM1

	
 	IMPLICIT REAL*8(A-H,O-Z)
CC
	REAL*8  XYZ(NDOF,MP)
	INTEGER	ELEM(N3D,ME)
	INTEGER ElemID(ME),WPCH(N2D+2,MW)
	INTEGER FaceNum(MK),NNPS(MK,NFACE),Local(MK,NFACE,N2D)
	INTEGER NNPE(MK),NEB(MW)
	REAL*8	CTR(3)
	


	ALLOCATE (WK1(MP),WK2(MP),stat=ierr)
	IF(ierr.NE.0) STOP 'allocate work arrays in SurfaceExtract()'

C	 THIS CODE IS USED TO EXTRACT THE SURFACE NODE AND PATCHES FROM A 
C	GLOBAL MESH.
	

	!reset
	WK1 = 0
	WK2 = 0


	!to judge wether this is a surface mesh
	
c	MaxFaceNum = 0
c	DO 50 I = 1, NE
c		IF(ElemID(I).GT.0) THEN
c			MaxFaceNum = MAX(FaceNum(ElemID(I)),MaxFaceNum)
c		END IF
c50	CONTINUE
		!It is already a surface mesh, return
c	IF(MaxFaceNum.LE.1) RETURN
		


	!WK1 counts the time of appearance of a node as the minimum one in one element, 
	ic=0
	DO 100 IE = 1,NE
		
		IK = ElemID(IE)
		DO 90 IFACE = 1, FaceNum(IK)
			

			MINID = MP+1
			DO 80 J= 1, NNPS(IK,IFACE)
				ID = Local(IK,IFACE,J)
				ID = ELEM(ID,IE)
				MINID = MIN(ID,MINID)

80			CONTINUE
			WK1(MINID) = WK1(MINID) + 1
			ic = ic +1
90		CONTINUE
		
100	CONTINUE
	
	
	!WK2 marks the start point in WK3 and WK3D
	WK2(1) = 0
	DO 200 IP = 2,NP
		WK2(IP) = WK2(IP-1) + WK1(IP-1)
200	CONTINUE
	WK1 = 0
	
	MM=WK2(NP)
	ALLOCATE(VWK3D(NDOF,MM),PATCH(N2D,MM),WK3(MM),WK4(MM),WK5(MM),
     &		  stat=ierr)
	IF(ierr.NE.0) STOP 'allocate work arrays in SurfaceExtract()'
	WK3 = 0
	WK4 = 0
	WK5 = 0
	VWK3D = 0

	
	!Extract the faces of each element, and mark the internal ones
	!WK3 stacks the face id, WK3D stacks its center point
	dvalue = (XYZ(1,1)-XYZ(1,2))*(XYZ(1,1)-XYZ(1,2))+
     &		 (XYZ(2,1)-XYZ(2,2))*(XYZ(2,1)-XYZ(2,2))+
     &		 (XYZ(3,1)-XYZ(3,2))*(XYZ(3,1)-XYZ(3,2))
	dvalue = SQRT(dvalue)

	

	IPch = 0
	Identical = 0
	DO 300 IE = 1,NE
		
		IK = ElemID(IE)
		DO 290 IFACE = 1, FaceNum(IK) 
		

			MINID = MP+1
			CTR(1) = 0
			CTR(2) = 0
			CTR(3) = 0
			DO 220 J= 1, NNPS(IK,IFACE)
				ID = Local(IK,IFACE,J)
				ID = ELEM(ID,IE)
				MINID = MIN(ID,MINID)
				CTR(1) = CTR(1) + XYZ(1,ID)/NNPS(IK,IFACE)
				CTR(2) = CTR(2) + XYZ(2,ID)/NNPS(IK,IFACE)
				CTR(3) = CTR(3) + XYZ(3,ID)/NNPS(IK,IFACE)
220			CONTINUE

		

			Mark = -1
			I = 0
			DO 250 WHILE(I.LT.WK1(MINID).AND.Mark.EQ.-1)
				
				I = I + 1
				ID = I + WK2(MINID)
			
				dis = (VWK3D(1,ID)-CTR(1))*(VWK3D(1,ID)-CTR(1))+
     &				  (VWK3D(2,ID)-CTR(2))*(VWK3D(2,ID)-CTR(2))+
     &				  (VWK3D(3,ID)-CTR(3))*(VWK3D(3,ID)-CTR(3))
				dis = SQRT(dis)

	
				IF(dis.LE.1.0E-12*dvalue)	Mark = ID	!identical face
					
250			CONTINUE

c		 if(mark.GT.0) then
c			WRITE(IUT0,*) 'IE=',IE, 'Iface=',iface,'mark=',mark
c c    &			,minid,wk1(minid)
c		end if

			IF(Mark.EQ.-1) THEN		
				!it is new face, stack it
				IPch = IPch + 1
				
		
				WK4(IPch) = NNPS(IK,IFACE)
						!WK4 marks the nnpes of the patch
				WK5(IPch) = IE
						!WK5 records the neighboring element id

							
				DO 260 J = 1, NNPS(IK,IFACE)
					K = Local(IK,IFace,J)
					PATCH(J,IPch) = ELEM(K,IE)
260				CONTINUE


				WK1(MINID) = WK1(MINID) + 1

				ID = WK2(MINID) + WK1(MINID)
				WK3(ID) = IPch
				VWK3D(1,ID) = CTR(1)
				VWK3D(2,ID) = CTR(2)
				VWK3D(3,ID) = CTR(3)

											
			ELSE
			
				!it is inner face, delete from the stack array
				Identical = Identical + 1

				JP = WK3(Mark)
				WK4(JP) = -WK4(JP) 
											!------iPch = 1
											!			  2		
											!			 ...
											!			 JP		---->	-1	XXXX
											!			 ...			
			
				!//
c				ID = WK1(MINID) + WK2(MINID)
c				WK3(Mark) = WK3(ID)
c				VWK3D(1,Mark) = VWK3D(1,ID)
c				VWK3D(2,Mark) = VWK3D(2,ID)
c				VWK3D(3,Mark) = VWK3D(3,ID)

c				WK1(MINID) = WK1(MINID) - 1

											!-----------wk2(I)
											!			 +1		
											!			 ...
											!			 +mark		<----      
											!			 ...			|
											!			 +wk1(I)-1		|
											!			 +wk1(I)	-----
		
			
			END IF	
	
290		CONTINUE

300	CONTINUE



	!project the wall patches into the global surface patches
	
	DO 400 IW = 1,NW
	
		MINID = MP+1
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0
		DO 320 J= 1, WPCH(1,IW)
			ID = WPCH(J+1,IW)
			MINID = MIN(ID,MINID)
			
			CTR(1) = CTR(1) + XYZ(1,ID)/WPCH(1,IW)
			CTR(2) = CTR(2) + XYZ(2,ID)/WPCH(1,IW)
			CTR(3) = CTR(3) + XYZ(3,ID)/WPCH(1,IW)
320		CONTINUE

		Mark = -1
		I = 0
		DO 350 WHILE(I.LT.WK1(MINID).AND.Mark.EQ.-1)
				
			I = I + 1
			ID = I + WK2(MINID)
		
			dis = (VWK3D(1,ID)-CTR(1))*(VWK3D(1,ID)-CTR(1))+
     &			  (VWK3D(2,ID)-CTR(2))*(VWK3D(2,ID)-CTR(2))+
     &			  (VWK3D(3,ID)-CTR(3))*(VWK3D(3,ID)-CTR(3))
			dis = SQRT(dis)
			
			IF(dis.LE.1.0E-10*dvalue) THEN
				Mark = ID			!identical face
				JP = WK3(Mark)
				IF(WK4(JP).EQ.0.OR.ABS(WK4(JP)).NE.WPCH(1,IW)) THEN
					WRITE(IUT0,*) 'Error abort in SurfaceExtract():'
					WRITE(IUT0,*) 'the patches are not identical'
					STOP
				ELSE
					WK4(JP) = 100+IW
				END IF
				NEB(IW) = WK5(JP)

			END IF
					
350		CONTINUE
		
		IF(Mark.EQ.-1) THEN
		WRITE(IUT0,*) 'Error abort in SurfaceExtract():'
		WRITE(IUT0,*) 'find no corresponding patch in the global_data'
		STOP
		END IF

400	CONTINUE


	
	
	!reclaim the interface or slide boundary
	DO 450 I=1,IPch
	  IF(WK4(I).LT.0) THEN
		Mark=0
		ID1=NodeBounID(PATCH(1,I))
		IF(ID1.EQ.0) CYCLE

		DO J = 2, -WK4(I)
		  ID2 = NodeBounID(PATCH(J,I))
		  IF(ID2.EQ.0) THEN
			Mark=1
			EXIT
		  ELSE IF(ID2.GT.0.AND.ID2.NE.ID1) THEN
			Mark=1
			EXIT
		  END IF
		END DO
		IF(Mark.EQ.0) WK4(I)=-WK4(I)

	  END IF
450	CONTINUE


	!count the surface patches and mark the surface node
	WK1 = 0
	IE = NW
	DO 500 I = 1, IPch
		IF(WK4(I).GT.0.AND.WK4(I).LT.100) THEN
			!append new patches
			IE = IE + 1
			WPCH(1,IE) = WK4(I)
			DO 430 J = 1, WK4(I)
				WPCH(J+1,IE) = PATCH(J,I)
430			CONTINUE
			NEB(IE) = WK5(I)
			WPCH(WK4(I)+2,IE) = NBOUN+1
		ELSE IF(WK4(I).GT.100)	THEN	
			!rearrange the connectivity of old patches
			IW = WK4(I)-100
			DO 440 J = 1, WPCH(1,IW)
				WPCH(J+1,IW) = PATCH(J,I)
440			CONTINUE
			
		END IF

500	CONTINUE
	
	NWT = IE
	
	
			
c	WRITE(IUT0,*) '  --------- < Extracted Surface Mesh > ------------'
	WRITE(IUT0,'(a,I8)') '  Total_Surface_Cell    = ', NWT
	WRITE(IUT0,'(a,I8)') '  Initial_Surface_Cell  = ', NW
	WRITE(IUT0,'(a,I8)') '  Appended_Surface_Cell = ', NWT-NW
c	WRITE(IUT0,*) '  -------------------------------------------------'


	DEALLOCATE (VWK3D,PATCH,WK1,WK2,WK3,WK4,WK5)


	RETURN

!------> end of SurfaceExtract()
	END







C======================================================================C
C                                                                      C
C Subroutine SurfaceExtract()									         C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/06/09                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	

	SUBROUTINE	SurfaceExtract_old(IUT0,XYZ,ELEM,VWK3D,WK1,WK2,WK3,WK4,
     &				WK5,ElemID,FaceNum,NNPS,Local,NNPE,WPCH,PATCH,
     &				NP,NE,NW,NWT,MP,ME,MM,MW,NDOF,N2D,N3D,NFACE,
     &				MK,NEB,NBOUN)

	
 	IMPLICIT REAL*8(A-H,O-Z)
CC
	REAL*8  XYZ(NDOF,MP),VWK3D(NDOF,MM)
	INTEGER	ELEM(N3D,ME),WK1(MM),WK2(MM),WK3(MM),WK4(MM),WK5(MM)
	INTEGER ElemID(ME),PATCH(N2D,MM),WPCH(N2D+2,MW)
	INTEGER FaceNum(MK),NNPS(MK,NFACE),Local(MK,NFACE,N2D)
	INTEGER NNPE(MK),NEB(MW)
	REAL*8	CTR(3)
	

C	 THIS CODE IS USED TO EXTRACT THE SURFACE NODE AND PATCHES FROM A 
C	GLOBAL MESH.
	

	!reset
	WK1 = 0
	WK2 = 0
	WK3 = 0
	WK4 = 0
	WK5 = 0
	VWK3D = 0
	


	!to judge wether this is a surface mesh
	
c	MaxFaceNum = 0
c	DO 50 I = 1, NE
c		IF(ElemID(I).GT.0) THEN
c			MaxFaceNum = MAX(FaceNum(ElemID(I)),MaxFaceNum)
c		END IF
c50	CONTINUE
		!It is already a surface mesh, return
c	IF(MaxFaceNum.LE.1) RETURN
		


	!WK1 counts the time of appearance of a node as the minimum one in one element, 
	ic=0
	DO 100 IE = 1,NE
		
		IK = ElemID(IE)
		DO 90 IFACE = 1, FaceNum(IK)
		
			MINID = MP+1
			DO 80 J= 1, NNPS(IK,IFACE)
				ID = Local(IK,IFACE,J)
				ID = ELEM(ID,IE)
				MINID = MIN(ID,MINID)

80			CONTINUE
			WK1(MINID) = WK1(MINID) + 1
			ic = ic +1
90		CONTINUE
		
100	CONTINUE
	
	
	!WK2 marks the start point in WK3 and WK3D
	WK2(1) = 0
	DO 200 IP = 2,NP
		WK2(IP) = WK2(IP-1) + WK1(IP-1)
200	CONTINUE
	WK1 = 0


	
	!Extract the faces of each element, and mark the internal ones
	!WK3 stacks the face id, WK3D stacks its center point
	dvalue = (XYZ(1,1)-XYZ(1,2))*(XYZ(1,1)-XYZ(1,2))+
     &		 (XYZ(2,1)-XYZ(2,2))*(XYZ(2,1)-XYZ(2,2))+
     &		 (XYZ(3,1)-XYZ(3,2))*(XYZ(3,1)-XYZ(3,2))
	dvalue = SQRT(dvalue)

	

	IPch = 0
	Identical = 0
	DO 300 IE = 1,NE
		
		IK = ElemID(IE)
		DO 290 IFACE = 1, FaceNum(IK) 
		

			MINID = MP+1
			CTR(1) = 0
			CTR(2) = 0
			CTR(3) = 0
			DO 220 J= 1, NNPS(IK,IFACE)
				ID = Local(IK,IFACE,J)
				ID = ELEM(ID,IE)
				MINID = MIN(ID,MINID)
				CTR(1) = CTR(1) + XYZ(1,ID)/NNPS(IK,IFACE)
				CTR(2) = CTR(2) + XYZ(2,ID)/NNPS(IK,IFACE)
				CTR(3) = CTR(3) + XYZ(3,ID)/NNPS(IK,IFACE)
220			CONTINUE

		

			Mark = -1
			I = 0
			DO 250 WHILE(I.LT.WK1(MINID).AND.Mark.EQ.-1)
				
				I = I + 1
				ID = I + WK2(MINID)
			
				dis = (VWK3D(1,ID)-CTR(1))*(VWK3D(1,ID)-CTR(1))+
     &				  (VWK3D(2,ID)-CTR(2))*(VWK3D(2,ID)-CTR(2))+
     &				  (VWK3D(3,ID)-CTR(3))*(VWK3D(3,ID)-CTR(3))
				dis = SQRT(dis)

	
				IF(dis.LE.1.0E-12*dvalue)	Mark = ID	!identical face
					
250			CONTINUE

c		 if(mark.GT.0) then
c			WRITE(IUT0,*) 'IE=',IE, 'Iface=',iface,'mark=',mark
c c    &			,minid,wk1(minid)
c		end if

			IF(Mark.EQ.-1) THEN		
				!it is new face, stack it
				IPch = IPch + 1
						
				WK4(IPch) = NNPS(IK,IFACE)
						!WK4 marks the nnpes of the patch
				WK5(IPch) = IE
						!WK5 records the neighboring element id

							
				DO 260 J = 1, NNPS(IK,IFACE)
					K = Local(IK,IFace,J)
					PATCH(J,IPch) = ELEM(K,IE)
260				CONTINUE


				WK1(MINID) = WK1(MINID) + 1

				ID = WK2(MINID) + WK1(MINID)
				WK3(ID) = IPch
				VWK3D(1,ID) = CTR(1)
				VWK3D(2,ID) = CTR(2)
				VWK3D(3,ID) = CTR(3)

											
			ELSE
			
				!it is inner face, delete from the stack array
				Identical = Identical + 1

				JP = WK3(Mark)
				WK4(JP) = 0
											!------iPch = 1
											!			  2		
											!			 ...
											!			 JP		---->	-1	XXXX
											!			 ...			
			
				!//
c				ID = WK1(MINID) + WK2(MINID)
c				WK3(Mark) = WK3(ID)
c				VWK3D(1,Mark) = VWK3D(1,ID)
c				VWK3D(2,Mark) = VWK3D(2,ID)
c				VWK3D(3,Mark) = VWK3D(3,ID)

c				WK1(MINID) = WK1(MINID) - 1

											!-----------wk2(I)
											!			 +1		
											!			 ...
											!			 +mark		<----      
											!			 ...			|
											!			 +wk1(I)-1		|
											!			 +wk1(I)	-----
		
			
			END IF	
	
290		CONTINUE

300	CONTINUE



	!project the wall patches into the global surface patches
	DO 400 IW = 1,NW
	
		MINID = MP+1
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0
		DO 320 J= 1, WPCH(1,IW)
			ID = WPCH(J+1,IW)
			MINID = MIN(ID,MINID)
			CTR(1) = CTR(1) + XYZ(1,ID)/WPCH(1,IW)
			CTR(2) = CTR(2) + XYZ(2,ID)/WPCH(1,IW)
			CTR(3) = CTR(3) + XYZ(3,ID)/WPCH(1,IW)
320		CONTINUE
	

		Mark = -1
		I = 0
		DO 350 WHILE(I.LT.WK1(MINID).AND.Mark.EQ.-1)
				
			I = I + 1
			ID = I + WK2(MINID)
		
			dis = (VWK3D(1,ID)-CTR(1))*(VWK3D(1,ID)-CTR(1))+
     &			  (VWK3D(2,ID)-CTR(2))*(VWK3D(2,ID)-CTR(2))+
     &			  (VWK3D(3,ID)-CTR(3))*(VWK3D(3,ID)-CTR(3))
			dis = SQRT(dis)

	
			IF(dis.LE.1.0E-12*dvalue) THEN
				Mark = ID	!identical face
				JP = WK3(Mark)
				IF(WK4(JP).EQ.0.OR.WK4(JP).NE.WPCH(1,IW)) THEN
					WRITE(IUT0,*) 'Error abort in SurfaceExtract():'
					WRITE(IUT0,*) 'the patches are not identical'
					STOP
				ELSE
					WK4(JP) = -IW
				END IF
				NEB(IW) = WK5(JP)

			END IF
					
350		CONTINUE
		
		IF(Mark.EQ.-1) THEN
		WRITE(IUT0,*) 'Error abort in SurfaceExtract():'
		WRITE(IUT0,*) 'find no corresponding patch in the global_data'
		STOP
		END IF

400	CONTINUE


	
	!count the surface patches and mark the surface node

	WK1 = 0
	IE = NW
	DO 500 I = 1, IPch
		IF(WK4(I).GT.0)	THEN
			!append new patches
			IE = IE + 1
			WPCH(1,IE) = WK4(I)
			DO 430 J = 1, WK4(I)
				WPCH(J+1,IE) = PATCH(J,I)
430			CONTINUE
			NEB(IE) = WK5(I)
			WPCH(WK4(I)+2,IE) = NBOUN+1
		ELSE IF(WK4(I).LE.-1)	THEN	
			!rearrange the connectivity of old patches
			IW = -WK4(I)
			DO 440 J = 1, WPCH(1,IW)
				WPCH(J+1,IW) = PATCH(J,I)
440			CONTINUE
			
		END IF

500	CONTINUE
	
	NWT = IE
	!IF(NWT.GT.NW)	NBOUN = NBOUN + 1
	
	
			
c	WRITE(IUT0,*) '  --------- < Extracted Surface Mesh > ------------'
	WRITE(IUT0,*) '  Total_PATCH     = ', NWT
	WRITE(IUT0,*) '  Origen_PATCH    = ', NW
	WRITE(IUT0,*) '  Appended_PATCH  = ', NWT-NW
c	WRITE(IUT0,*) '  -------------------------------------------------'


	RETURN

!------> end of SurfaceExtract()
	END







C======================================================================C
C                                                                      C
C Subroutine GetToken()									             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/05                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

!C
!C  Utility Subroutine for processing data
!C
	SUBROUTINE GetToken(NF,TOKEN,LENGTH,NC)
	
	CHARACTER(LEN=LENGTH) TOKEN
	INTEGER NF,N

	N = 0
	DO 10 I = 1, LENGTH
		
		IF(token(i:i).NE.'') THEN
			N = N + 1
		ELSE
			GOTO 20 
		END IF

10	CONTINUE
	


20	NC = N
	

c	WRITE(*,*) 'The last character = ', TOKEN(n:n)
c	WRITE(*,*) 'NF = ', NF



	RETURN

!---------->
	END 





	
C======================================================================C
C                                                                      C
C Subroutine GetHead()									             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/05                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

!C
!C  Utility Subroutine for processing data
!C
	SUBROUTINE GetHead(TOKEN,LENGTH,NC)
	
	CHARACTER(LEN=LENGTH) TOKEN
	INTEGER NF,N

	N = 0
	DO 10 I = 1, LENGTH
		
		IF(token(i:i).EQ.' ') THEN
			N = N + 1
		ELSE
			GOTO 20 
		END IF

10	CONTINUE
	


20	NC = N+1
	

c	WRITE(*,*) 'The last character = ', TOKEN(n:n)
c	WRITE(*,*) 'NF = ', NF



	RETURN

!---------->
	END 




	
C======================================================================C
C                                                                      C
C Subroutine GetArraySize()								             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/05                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

	SUBROUTINE	GetArraySize(IUT0,NNN,MRD0,NRD0,NRD,MID,NE,NW,VMIN,
     &				AveGCoef,AveWGray,AveVol,confidence,Model,
     &				NRAY,StackOpt)
	IMPLICIT REAL*8(A-H,O-Z)
	INCLUDE  'pre_common.f'
		
	PARAMETER (PAI	= 3.14159265)

	CHARACTER*20  StackOpt,Model
	REAL*8	eta(50)
	

	

C0----------------------------------------------------------------------
	DATA   eta / 1.82116, 1.64486, 1.53440, 1.45216, 1.38586,
     &			 1.32988, 1.28118, 1.23790, 1.19880, 1.16306,
     &			 1.13008, 1.09938, 1.07062, 1.04352, 1.01788,
     &			 0.99352, 0.97028, 0.94804, 0.92670, 0.90618,
     &			 0.88640, 0.86728, 0.84878, 0.83084, 0.81342,
     &			 0.79648, 0.77998, 0.76390, 0.74820, 0.73286,
     &			 0.71786, 0.70318, 0.68880, 0.67470, 0.66086,
     &			 0.64726, 0.63390, 0.62076, 0.60784, 0.59512,
     &			 0.58258, 0.57022, 0.55804, 0.54602, 0.53416,
     &			 0.52244, 0.51086, 0.49942, 0.48812, 0.47694 /

	
	CoefErr = 0.1


C1---------------------	ZONE METHOD	------------------------------
	



	IF(MODEL(1:4).EQ.'ZONE') THEN
		!stack as a 1-D array anologous to a
		!Lower-Triangle-Array
		NNN = NRD
		StackOpt = 'LTRI'
		MRD0 = NRD*(NRD+1)/2
				
		VMIN = 1.0/NRD*VMIN
		
		MID = 1

		
		WRITE(IUT0,*) '  NRD =', NRD
		WRITE(IUT0,*) '  Exchange-Area Save Mode=', StackOpt(1:4)
			
		RETURN

	END IF





C1---------------------	MONTE-CARLO METHOD	------------------------------

	I1 = (100 - CONFIDENCE*100)
	I1 = MIN(I1,50)
	I1 = MAX(I1,1)
	eta1 = eta(I1)

	IF(AveGCoef.GT.0) THEN
		Dist0 = -1.0*LOG(1.0-0.95)/AveGCoef
		NNN0 = 4.0/3.0*PAI*Dist0**3/AveVol
		
	ELSE
		NNN0 = NW
	END IF
	
	NNN = MIN(NNN0,NRD)

	
c	IF(MC%RaySpecOpt(1:4).EQ.'AUTO')	THEN
	IF(NRAY.LE.0) THEN
	NRAY = 1.0*NNN/CoefErr/10
	NRAY = MIN(10000,NRAY)
	NRAY = MAX(100,NRAY)
	END IF
	

c	IF(NNN.GT.NRD/4) THEN
		!stack as a 1-D array anologous to a
		!Lower-Triangle-Array
		NNN = NRD
		StackOpt = 'LTRI'
		MRD0 = NRD*(NRD+1)/2
				
		VMIN = 1.0/NRD*VMIN
		
		MID = 1

c	ELSE
		!stack as a 1-D array anologous with band NNN
		!and element indices
c		StackOpt = 'BAND'
c
c		VMIN = 1.0/NNN*VMIN
c		 
c		!enlarge
c		MRD0 = NRD*NNN
c		MID  = MRD0
c
c	END IF
	
c	WRITE(IUT0,*) '  NRAY=',NRAY
c	WRITE(IUT0,*) '  MRD0=',MRD0
c	WRITE(IUT0,*) '  NRD =', NRD
c	WRITE(IUT0,*) '  Exchange-Area Save Mode=', TRIM(StackOpt)
	!WRITE(IUT0,*) '    VMIN=',VMIN

		!V1MIN--	the lower limit of main READs values, 
		!			the READs above this value will appear in the coef. matrix of the equation
		!			It means the ratio to the average.
		!V2MIN--	the lower limit of small READs values, 
		!			the READs above this value will appear in the equation. as source terms,
		!			It means the ratio to the average.
		!			the source terms will be calculated individually during the solving iteration
		!			the READs under this value will appear in the equation as one group source term




	RETURN
!------>

	END 





C======================================================================C
C                                                                      C
C Subroutine SymmTreatRD_1()								             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/05                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

	SUBROUTINE	SymmTreatRD1(Dij,Dji,Ai,Aj)
	IMPLICIT REAL*8(A-H,O-Z)	
	
	!this symmetry_treatment scheme is proposed by
	!R.I.Loehrcke, J.S.Dolaghan, and P.J.Burns, Smoothing Monte-Carlo Exchange Factors,
	!ASME-T J. Heat Transfer, V.117, P524-526,1995.
	
	!conf:		confidence
	!Dij:		exchange area
	!
	!
	Fij = Dij/Ai
	Fji = Dji/Aj
	IF(Fij.EQ.0) THEN
		Cij = 1
		Cji = 0
	ELSE IF(Fji.EQ.0) THEN
		Cji = 1
		Cij = 0
	ELSE
		Cij = (1-Fij)/Fij
		Cji = (1-Fji)/Fji
	END IF
	
	
	c1 = Cji/(Cij+Cji)
	!c2 = (cij*cij)/(cij*cij+cji*cji)
	
	Dij = c1*Dij + (1-c1)*Dji
	Dji = Dij


	END SUBROUTINE
	




	

C======================================================================C
C                                                                      C
C Subroutine SymmTreatRD_2()								             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/05                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

	SUBROUTINE	SymmTreatRD(Dij,Dji,Ai,Aj)
	IMPLICIT REAL*8(A-H,O-Z)
	
	!this symmetry_treatment scheme is proposed by
	!T. Omori, J.H. Yang, S. Kato, S. Murakami, Radiative Heat Transfer Analysis Method for
	!Coupled Simulation of Convection and Radiation in Large-Scale and Complicated Enclusures
	!SHASE of Japan, No.88, 103-113, 2003.1
	
	!am:		exponent of Ai and Aj used in area-weighted averaging
	!Dij:		exchange area
	am = 1

	!value = (Dij*(Aj**am)+Dji*(Ai**am))/(Aj**am+Ai**am)
	Dij = (Dij*Aj+Dji*Ai)/(Aj+Ai)

	Dji = Dij


C-----------------------------------------------------------------------------------
C	NOTE:	this method is proved to be better,			Jiang Yuyan. 2005/08/22/
C

	END SUBROUTINE
C======================================================================C
	SUBROUTINE	SmoothLTRI_1(IUT0,NN,MM,StaPoint,D,A,X,WK,B,
     &					   feps,vfeps,RLX,iterate)
C======================================================================C

	
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8 D(MM),A(NN),WK(NN),X(NN),B(NN)
	INTEGER StaPoint(NN+1)
     
		
C----------------------------------------------------------------------------------C
c	this is a smoothing program for exchange-area
c	the final value obeys the conservation-law and reciprocity-relation:
c			sum(Dij)j = Ai
c			Dij = Dji
c	this program is correspondent to StackOpt='LTRI'
c	
c	this symmetry_treatment scheme is proposed by
c	R.I.Loehrcke, J.S.Dolaghan, and P.J.Burns, Smoothing Monte-Carlo Exchange Factors,
c	ASME-T J. Heat Transfer, V.117, P524-526,1995.
c	
c	&
c
c	M.E. Larsen and J.R. Howell, Least-Squares Smoothing of Direct-Exchange Areas in 
c	Zonal Analysis, ASME-T J. Heat Transfer, V.108, P.239-242, 1986.
C----------------------------------------------------------------------------------C
C
C		EQUATION:		[R] X = B
C
C				Rij = Wij,		Rii = Wii + sum(Wij)j
C				Bi = Ai - sum(Dij)j
C				Wij = Dij^2
C
C							SYMMETRY:	Rij=Rji
C		Correction:
C				Dij = Dij + Wij(Xi+Xj)
C----------------------------------------------------------------------------------C
	iterate = 0
	residence1 = 0

	!Stack Wij
10	WK = 0.0
	B(1:NN) = A(1:NN)
	
	iflag = 0
	tv = 0
	
	DO 100 I = 1,NN
		ID = StaPoint(I)
		DO J=1,I
			WK(I) = WK(I) + D(ID+J)*D(ID+J)
			B(I)  = B(I) - D(ID+J)
		END DO
		DO J=I+1,NN
			ID1 = StaPoint(J)
			WK(I) = WK(I) + D(ID1+I)*D(ID1+I)
			B(I)  = B(I) - D(ID1+I)
		END DO

		if(iterate.eq.0) residence1= residence1 + abs(B(I))	
		tv = tv + sqrt(wk(i))
100	CONTINUE
	if(iterate.eq.0.and.tv.gt.0) residence1 = residence1/tv
	!write(iut0,*) 'residence=',residence1
	IF(residence1.LT.vfeps) GOTO 1999
	
	!get correction 
	err = 1000.0
	itr = 0
	iflag = 1
	X = 0

	DO 1000 WHILE(err.GT.feps)
		itr = itr + 1
		err = 0

		DO 200 I = 1,NN
			
			
			bb = B(I)
			ID = StaPoint(I)
			DO J=1,I-1
				bb = bb - (D(ID+J)*D(ID+J))*X(J)
			END DO
			
			DO J=I+1,NN
				ID1 = StaPoint(J)
				bb = bb - (D(ID1+I)*D(ID1+I))*X(J)
			END DO
			vlast = X(I)
			vpre = bb/(D(ID+I)*D(ID+I)+WK(I))
			X(I) = X(I)*(1-rlx) + vpre*rlx
			err = err + ABS(X(I)-vlast)
	!if(abs(x(i)).gt.10000) then
c		write(iut0,*) I,ID,x(I),xlast
c		write(iut0,*) bb,b(i),d(id+i),D(ID+I)*D(ID+I)+WK(I)
c		pause
	!end if

200		CONTINUE
		err = err/NN
		WRITE(IUT0,*) '  STEP=',itr, 'ERR=',err, ' / ',feps
1000	CONTINUE


	!correct the Exchange-Area
	B = A
	residence2 = 0
	minus = 0
	DO 2000 I=1,NN
		!correct
		ID = StaPoint(I)
		DO J=1,I-1
			D(ID+J) = D(ID+J)+(D(ID+J)*D(ID+J))*(X(I)+X(J))
			IF(D(ID+J).LT.0) THEN
				
				!WRITE(IUT0,*) 'ERROR Aborb in SmoothLTRI():'
				!WRITE(IUT0,*) D(ID+J)
				!WRITE(IUT0,*) 'Exchange-area < 0 at(',I,',',J,')'
				!STOP
				minus = minus + 1
			END IF
		END DO

		!check conservation
		ID = StaPoint(I)
		DO J=1,I
			B(I)  = B(I) - D(ID+J)
		END DO
		DO J=I+1,NN
			ID1 = StaPoint(J)
			B(I)  = B(I) - D(ID1+I)
		END DO

		residence2 = residence2 + abs(B(I))

2000	CONTINUE
	
	WRITE(IUT0,*) '  Found Negative Values:', minus
	if(tv.GT.0) residence2 = residence2/tv
	
	
	iterate = iterate+1
	if(residence2.GT.vfeps) GOTO 10
	
1999	IF(iflag.EQ.1) THEN
	 correct = residence2/residence1*100
	 WRITE(IUT0,*) '  Total_Residence:'
	 WRITE(IUT0,*) '  Before Smoothing = ',residence1
	 WRITE(IUT0,*) '  After  Smoothing = ',residence2,correct,'%'
	 WRITE(IUT0,*) '  Correction times = ',iterate
	ELSE
	 WRITE(IUT0,*) '  Total_Residence = ',residence1, 
     &	 '  ...need no correction'
	END IF

	RETURN

!-------> end of SmoothLTRI1()


	END


C======================================================================C
	SUBROUTINE	SmoothLTRI_2(IUT0,NRD,MRD,RDIndex,RDValue,PROP,VWK2,
     3						feps,vfeps,RLX,iterate)
C======================================================================C
	
	
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8 VWK2(NRD)
	REAL*8 RDValue(MRD),PROP(NRD)
	INTEGER RDIndex(NRD+1)
     

C----------------------------------------------------------------------------------C
c	this is a smoothing program for exchange-area
c	the final value obeys the conservation-law and reciprocity-relation:
c			sum(Dij)j = 1
c			Dij = Dji
c	this program is correspondent to StackOpt='LTRI'
c	
c	this symmetry_treatment scheme is proposed by
c	T. Omori, J.H. Yang, S. Kato, S. Murakami, Radiative Heat Transfer Analysis Method for
c	Coupled Simulation of Convection and Radiation in Large-Scale and Complicated Enclusures
c	SHASE of Japan, No.88, 103-113, 2003.1
c		
C----------------------------------------------------------------------------------C
	
	err = vfeps+1
	istep = 0


	DO 1000 WHILE(err.GT.vfeps)

		istep = istep + 1

		!sum
		err = 0
		VWK2(:) = 0.0
		tv = 0
		DO 100 IRD = 1,NRD
			
			!Fii is save in RDSelf()
			
			!for left part J<I, Dij  (Dij = Ai*Fij)
			DO J=1,IRD
				ID = RDIndex(IRD) + J		
				VWK2(IRD) = VWK2(IRD)+ RDValue(ID)
				
			END DO

			!for left part J<I, Dij -> Dji
			DO J=IRD+1,NRD			!Fii is save in RDSelf()
				ID = RDIndex(J) + IRD
				VWK2(IRD) = VWK2(IRD)+ RDValue(ID)
			END DO
			
			err = err + ABS((VWK2(IRD)-PROP(IRD)))/PROP(IRD)
100		CONTINUE
		err = err/NRD
		if(tv.GT.0) err = err/tv

		!get new values
		DO 300 IRD=1,NRD
		  DO 300 J = 1, IRD
			a1 = PROP(IRD)*PROP(J)/(PROP(J)+PROP(IRD))
			ID = RDIndex(IRD)+J
			a2 = 1.0/VWK2(IRD)+1.0/VWK2(J)
			RDValue(ID) = a1*a2*RDValue(ID)
			

300		CONTINUE
		
		write(IUT0,'(I3,a,E12.3,a,E12.3)') istep,': Resi=  ',
     &			err,' / ',vfeps

1000	CONTINUE
	
	iterate = istep

	RETURN

!-------> end of SmoothLTRI2()


	END

C======================================================================C
        SUBROUTINE	SmoothBAND(IUT0,NE,NW,NRD,MRD0,MRD,MID,
     &					RDId,RDN1,RDN2,RDIndex,RDValue,
     &					PROP,VMIN,feps,vfeps,RLX,iterate)
C======================================================================C     

	USE MODULE_RADWORK, ONLY : VWK =>VWK4
	USE MODULE_RADWORK, ONLY : B => VWK2
	USE MODULE_RADWORK, ONLY : X => VWK3

	USE MODULE_RADWORK, ONLY : WK =>WK1
	USE MODULE_RADWORK, ONLY : WK2,WK3

	
	IMPLICIT REAL*8(A-H,O-Z)

	REAL*8  RDValue(MRD0),PROP(NRD)
	INTEGER RDId(MID),RDN1(NRD),RDN2(NRD),RDIndex(NRD+1)
	     	

	ALLOCATE(VWK(MRD0),B(NRD),X(NRD),WK(MRD0),WK2(NRD+1),WK3(NRD),
     &	stat=ierr)
	IF(ierr.NE.0) STOP 'allocate arrays in SmoothBAND()'

C----------------------------------------------------------------------------
C	Make starting points
C

	!num of Dij, going
	WK2 = 0
	DO 100 I = 1, MRD
		IRD = RDId(I)
		WK2(IRD) =  WK2(IRD) + 1
100	CONTINUE
	
	!num of Dji,coming
	WK3 = 0
	DO 110 IRD = 2, NRD+1
		WK3(IRD) = RDIndex(IRD+1)-RDIndex(IRD)
		M = MAX(WK3(IRD),WK2(IRD))
		WK2(IRD) =  WK2(IRD-1) + M		!WK2 marks new start-point
110	CONTINUE
	MRD = WK2(NRD+1)


	RDN1 = 0		!left large off-diagnal element
	RDN2 = 0		!total large off-diagnal element
	
	VWK = 0.0
	WK = 0

C----------------------------------------------------------------------------
C	Get Reciprocity Average
C

	DO 1000 IRD = 1,NRD
		tvalue = 0.0
		IMRD = RDIndex(IRD)+1
		VWK(IMRD) = RDValue(IMRD)
		
		DO 200 J = 2,RDIndex(IRD+1)-RDIndex(IRD)
			IMRD = IMRD + 1
			JRD = RDId(IMRD)

			!if J2 > IRD, take symmetric average
			IF(JRD.GT.IRD) THEN
				mark = 0
				K = 1
				NUM = RDIndex(JRD+1)-RDIndex(JRD)
				DO 150 WHILE(mark.EQ.0.AND.K.LT.NUM)
					K = K + 1
					JMRD = RDIndex(JRD) + K
					IRDP = RDId(JMRD)
					IF(IRDP.EQ.IRD) THEN
						a1 = RDValue(IMRD)
						a2 = RDValue(JMRD)
						CALL SymmTreatRD(a1,a2,PROP(IRD),PROP(JRD))
						ID = WK2(IRD)+J
						VWK(ID) = a1
						ID = WK2(JRD)+K
						VWK(ID) = a1
								!reverse the value sign, for informing that
								!they were treated
						RDValue(IMRD) = -a1		
						RDValue(JMRD) = -a1

						!found the counter-part
						mark = 1
								
					END IF

150				CONTINUE
				!found no counter-part
				IF(mark.EQ.0) THEN
					a1 = RDValue(IMRD)
					CALL SymmTreatRD(a1,0.d0,PROP(IRD),PROP(JRD))
					ID = WK2(IRD)+J
					VWK(ID) = a1
								!reverse the value sign, for informing that
								!it was treated
					RDValue(IMRD) = -a1

					!append it to jrd
					WK3(JRD) = WK3(JRD) + 1
					ID = WK2(JRD)+WK3(JRD)
					IF(ID.GT.WK2(JRD+1)) GOTO 4000
					VWK(ID) = a1

				END IF
			
			!if J2<IRD and Value > 0, it has not been treat and has no counter-part
			ELSE IF(RDValue(IMRD).GT.0)	THEN
				a1 = RDValue(IMRD)
				CALL SymmTreatRD(a1,0.d0,PROP(IRD),PROP(JRD))
				ID = WK2(IRD)+J
				VWK(ID) = a1
							!reverse the value sign, for informing that
							!it was treated
				RDValue(IMRD) = -a1

				!append it to JRD
				WK3(JRD) = WK3(JRD) + 1
				ID = WK2(JRD)+WK3(JRD)
				IF(ID.GT.WK2(JRD+1)) GOTO 4000
				VWK(ID) = a1
			END IF

200		CONTINUE
				

1000	CONTINUE

	RDId = WK
	RDValue = VWK
	RDIndex = WK2
	
C----------------------------------------------------------------------------
C	Smoothing
C



C----------------------------------------------------------------------------------C
c	this is a smoothing program for exchange-area
c	the final value obeys the conservation-law and reciprocity-relation:
c			sum(Dij)j = Ai
c			Dij = Dji
c	this program is correspondent to StackOpt='LTRI'
c	
c	this symmetry_treatment scheme is proposed by
c	R.I.Loehrcke, J.S.Dolaghan, and P.J.Burns, Smoothing Monte-Carlo Exchange Factors,
c	ASME-T J. Heat Transfer, V.117, P524-526,1995.
c	
c	&
c
c	M.E. Larsen and J.R. Howell, Least-Squares Smoothing of Direct-Exchange Areas in 
c	Zonal Analysis, ASME-T J. Heat Transfer, V.108, P.239-242, 1986.
C----------------------------------------------------------------------------------C
C
C		EQUATION:		[R] X = B
C
C				Rij = Wij,		Rii = Wii + sum(Wij)j
C				Bi = Ai - sum(Dij)j
C				Wij = Fij^2
C
C							SYMMETRY:	Rij=Rji
C		Correction:
C				Dij = Dij + Wij(Xi+Xj)
C----------------------------------------------------------------------------------C
	
	!Stack Wij
1111	VWK = 0.0
	residence1 = 0
	B = 0
	iflag = 0
	tv = 0
	DO 2100 I = 1,NRD
		err = PROP(I)
		DO J=RDIndex(I)+1,RDIndex(I+1)
			VWK(J) = RDValue(J)*RDValue(J)
			B(I) = B(I) + VWK(J)
			err  = err - RDValue(J)
		END DO
		VWK(RDIndex(I)+1) = VWK(RDIndex(I)+1) + B(I)
		
		tv = tv + sqrt(B(I))
		residence1 = residence1 + abs(err)	

2100	CONTINUE
	if(tv.gt.0) residence1=residence1/tv
	if(residence1.LT.vfeps) GOTO 3999
	
	!get correction 
	iterate = iterate+1
	iflag = 1
	err = 1000.0
	itr = 0
	DO 3000 WHILE(err.GT.feps)
		itr = itr + 1
		err = 0

		DO 2200 I = 1,NRD
			bb = B(I)
			
			DO J=RDIndex(I)+2,RDIndex(I+1)
				ID = RDId(J)
				bb = bb - VWK(J)*X(ID)
			END DO
			
			vlast = X(I)
			ID = RDIndex(I)+1
			vpre = bb/VWK(ID)
			X(I) = X(I)*(1-rlx) + vpre*rlx
			err = err + ABS(X(I)-vlast)
2200		CONTINUE
		err = err/NRD

		WRITE(IUT0,*) '  STEP=',itr, 'ERR=',err, ' / ',feps
3000	CONTINUE


	!correct the Exchange-Area
	VWK = 0.0
	DO 3100 I=1,NRD
		!correct
		DO J=RDIndex(I)+1,RDIndex(I+1)
			ID = RDId(J)
			VWK(J) = VWK(J)+(RDValue(J)*RDValue(J))*X(ID)
		END DO
3100	CONTINUE
	
	DO 3200 I=1,MRD
		!correct
		RDValue(I) = RDValue(I) + VWK(J)
3200	CONTINUE


	residence2 = 0
	DO 3300 I=1,NRD
		!check conservation
		err = PROP(I)
		DO J=RDIndex(I)+1,RDIndex(I+1)
			err  = err - RDValue(J)
		END DO
			
		residence2= residence2 + abs(err)	

3300	CONTINUE
	if(tv.gt.0) residence2=residence2/tv
	if(residence2.GT.vfeps) GOTO 1111

3999	IF(iflag.EQ.1) THEN
	 correct = residence2/residence1*100
	 WRITE(IUT0,*) '  Total_Residence:'
	 WRITE(IUT0,*) '  Before Smoothing = ',residence1
	 WRITE(IUT0,*) '  After  Smoothing = ',residence2,correct,'%'
	 WRITE(IUT0,*) '  Correction Times = ',iterate
	ELSE
	 WRITE(IUT0,*) '  Total_Residence = ',residence1, 
     &	 'need no correction'
	END IF

	DEALLOCATE(VWK,B,X,WK,WK2,WK3)



	RETURN

	

4000	WRITE(IUT0,*) 'ERROR ABORT IN  SmoothBAND():'
	WRITE(IUT0,*) 'The record number exceeds it original setting'
	STOP
!-------> end of SmoothBAND()

	END


C======================================================================C
	SUBROUTINE WallGroupEASum(IUT0,IUT1,NE,NW,NRD0,NRD,MRD,MID,N2D,
     &					NBOUN,StackOpt,RDIndex,RDValue,RDId,PROP,
     &					EASum,WElem,VIndex)
C======================================================================C
	

	USE module_radwork, ONLY : WK1

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  RDValue(MRD),PROP(NRD0),EASum(0:NBOUN,0:NBOUN),
     &		SumProp(0:NBOUN)
	INTEGER RDId(MID),RDIndex(NRD+1),WElem(N2D+2,NW),VIndex(NRD0)
	CHARACTER*20	StackOpt
	

	ALLOCATE(WK1(NRD0))

	MaxNum = 0

	WK1(:) = 0
	SumProp(:) = 0
	
	DO 10 I=1,NE
		ID = VIndex(I)
		IF(ID.GT.0) SumProp(0) = SumProp(0) + PROP(I)

10	CONTINUE

	DO 20 IW=1,NW
		IIT = IW+NE
		ID = VIndex(IIT)
		IN = WElem(1,IW)
		IF(ID.GT.0) THEN
			WK1(ID) = WElem(IN+2,IW)
			IF(WK1(ID).GT.0) 
     &			SumProp(WK1(ID)) = SumProp(WK1(ID)) + PROP(IIT)
		END IF
		
		MaxNum = MAX(MaxNum,WElem(IN+2,IW))
		
20	CONTINUE
	
	

	EASum = 0.0

	IF(StackOpt(1:4).EQ.'LTRI')	THEN

		DO 100 IRD = 1,NRD
			
			!for left part J<I, Dij  (Dij = Ai*Fij)
			DO J=1,IRD
				ID = RDIndex(IRD) + J		
				I1 = WK1(IRD)
				I2 = WK1(J)
				EASum(I1,I2) = EASum(I1,I2)+ RDValue(ID)		!/PROP(IRD)
			END DO

			!for left part J<I, Dij -> Dji
			DO J=IRD+1,NRD			!Fii is save in RDSelf()
				ID = RDIndex(J) + IRD
				I1 = WK1(IRD)
				I2 = WK1(J)
				EASum(I1,I2) = EASum(I1,I2)+ RDValue(ID)		!/PROP(IRD)
			END DO

100		CONTINUE

	ELSE

		DO 200 I = 1,NRD
			I1 = WK1(I)
			DO J=RDIndex(I)+1,RDIndex(I+1)
				I2 = WK1(RDId(J))
				EASum(I1,I2) = EASum(I1,I2) + RDValue(J)		!/PROP(I)
			
			END DO
		
200		CONTINUE
		
	
	END IF


	!fff = sum(rdvalue)

	DO 250 I=0,MaxNum
		DO 250 J=0,MaxNum
			IF(SumProp(I).GT.0) EASum(I,J) = EASum(I,J)/SumProp(I)
250	CONTINUE

	OPEN(IUT1,FILE='BExchangeArea-Summary.dat')
	
	WRITE(IUT1,*) 'Exchange-Area between Gas and Walls, <0: is Gas>'
	
	!NOTE:	The wall numbers are assumed to be continuous	
	WRITE(IUT1,'(I16, 20I12)') (I,I=0,MaxNum)
	DO I=0,MaxNum
		WRITE(IUT1,'(I4, 20F12.5)') I,(EASum(I,J),J=0,MaxNum)
	END DO

	CLOSE(IUT1)


	DEALLOCATE(WK1)	


	RETURN


!---------> end of WallGroupEASum()
	END 

C======================================================================C
	SUBROUTINE GetSegment(DV,DD,ratio)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8	DV(3),DD(3)
	IDOF = 1

	dmag = sqrt(DV(1)*DV(1) + DV(2)*DV(2) + DV(3)*DV(3))
	DV = DV/dmag

	vmax = abs(DV(1)/dd(1))
	DO 30 I = 2, 3
		va = abs(DV(I)/dd(I))
		IF(va.GT.vmax) THEN
		IDOF = I
		vmax = va
		END IF
30	CONTINUE
	
	dseg = 1.0/vmax/ratio

	
	RETURN


	END 

C======================================================================C	
	SUBROUTINE GetNodeNeb(IUT0,WElem,ELEM,
     &			ElemKind,NNPE,XYZ,NPIN,NBOUND,
     &			ME,MP,NP,NE,NW,NWT,N2D,N3D,NDOF,NKind)
C======================================================================C

	USE MODULE_RADSXF, ONLY : WNEB,RaNeb,VIndex,NodeBounID
	USE MODULE_RADWORK, ONLY : WK1

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8	XYZ(NDOF,MP)
	INTEGER WElem(N2D+2,NWT)
	INTEGER ELEM(N3D,ME),ElemKind(NE),NNPE(NKind)
	!INTEGER NNPS(NKIND,NFACE),CNTIVITY(NKIND,NFACE,N2D),FaceNum(NKIND)
	
	!mark participating nodes to allocate arrays
	

	ALLOCATE(WK1(NP))

	WK1=0
	DO IW=1,NW
		
		IF(WElem(WElem(1,IW)+2,IW).EQ.NBOUND+1) CYCLE		
				!interface import node, no need to cal. heat divergence

		DO J=1,WElem(1,IW)
			IP=WElem(J+1,IW)
			ID=NodebounID(IP)
			IF(ID.GT.0.AND.ID.NE.NBOUND+1) WK1(IP)=WNEB(IW)
		END DO
	
	END DO


	DO IE=1,NE
	  Mark=0
	  DO J=1,NNPE(ElemKind(IE))
		IP=ELEM(J,IE)
		IF(NodeBounID(IP).EQ.NBOUND+1) THEN
			Mark=1			!contains interface import and export vertex
			EXIT
		END IF
	  END DO
	  IF(Mark.EQ.0) CYCLE

	  DO J=1,NNPE(ElemKind(IE))
		IP=ELEM(J,IE)
		IF(NodeBounID(IP).EQ.0) THEN
			WK1(IP)=-IE
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!export vertex, (except other boundary vertex)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		END IF
	  END DO
	END DO


	DO IP=1,NP
	  IE=WK1(IP)
	  
	  IF(IE.GT.0) THEN	
	    !boundary node, its neighbor is internal node
		distmin = 1.0E30
		markdo = 0
		DO J=1,NNPE(ElemKind(IE))
			JP=ELEM(J,IE)
			IF(NodeBounID(JP).GT.0) CYCLE	
			
			!internal node
			dist=(XYZ(1,IP)-XYZ(1,JP))**2+(XYZ(2,IP)-XYZ(2,JP))**2
     &			 +(XYZ(3,IP)-XYZ(3,JP))**2
			IF(dist.LE.distmin) THEN
				distmin = dist
				WK1(IP)=JP		!the nearest intenal node
				markdo = 1
			END IF
		END DO
		
		IF(markdo.EQ.0) STOP 'can not find internal neighbor'
			  
	  ELSE IF(IE.LT.0) THEN		
	    !export node, its neighbor is boundary node (import node)
		IE=-IE
		distmin = 1.0E30
		markdo = 0
		DO J=1,NNPE(ElemKind(IE))
			JP=ELEM(J,IE)
			IF(NodeBounID(JP).NE.NBOUND+1) CYCLE	
			
			!import node
			dist=(XYZ(1,IP)-XYZ(1,JP))**2+(XYZ(2,IP)-XYZ(2,JP))**2
     &			 +(XYZ(3,IP)-XYZ(3,JP))**2
			IF(dist.LE.distmin) THEN
				distmin = dist
				WK1(IP)=-JP		!the nearest import node
				markdo=1
			END IF
		END DO
	
		IF(markdo.EQ.0) STOP 'can not find import neighbor'
	
	  END IF

	END DO

	
	DEALLOCATE(WNEB)
	ALLOCATE(WNEB(NP))
	WNEB(1:NP)=WK1(1:NP)

	
	DO I=1,NP			!see RDVElemToNode() or NodeToGPProp()

	  IF(WNEB(I).GT.0) THEN			
		!!!!!!!only need for value transfer from internal nodes to surface
		RaNeb(I)= RaNeb(I)/RaNeb(WNEB(I))
		
	  ELSE
		RaNeb(I)= 1.0
	  END IF

	END DO

c	open(1,file='neb1.dat')
c	write(1,*) np,nbound
c	do I=1,np
c		write(1,*) i,wneb(i),raneb(i),nodebounid(i)
c	end do
c	close(1)


	DO I=1,NP			!see RDVElemToNode() or NodeToGPProp()
	  WNEB(I)=ABS(WNEB(I))
	END DO

	
	DEALLOCATE(WK1)



	RETURN

!-------> 

	END
C======================================================================C
	SUBROUTINE GetPrdPair(IUT0,NBOUND,prdData,RadWallProp,
     &		MP,NP,ME,NE,NWT,NW,NDOF,N2D,N3D,NKIND,NFACE,
     &		XYZ,ELEM,FanWei,AveVol)
C======================================================================C
	
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID
	USE MODULE_RADSXF, ONLY : FaceNum,NNPE,NNPS,CNTIVITY,ElemKind
	USE MODULE_RADSXF, ONLY : WNV,WELEM,WGrayDeg,CTR

	USE module_radwork, ONLY : WK2

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8	radwallprop(NBOUND),prdData(NBOUND,5)
	REAL*8	XYZ(NDOF,MP),FanWei(2,NDOF)
	INTEGER MFind(NBOUND)
	INTEGER ELEM(N3D,ME)
	REAL*8	V0(3),VT(3),VW(3),DV(3),Seta,Eta



	!optimized Block Array size
	TotVol=(FanWei(2,1)-FanWei(1,1))*(FanWei(2,1)-FanWei(1,1))
     &		*(FanWei(2,1)-FanWei(1,1))
	NB= (TotVol/AveVol)**0.33*1.5
	wholedist = 0
	DO  I=1,3
		wholedist = max(wholedist,FanWei(2,I)-FanWei(1,I))
	END DO

	MFind = 0
	MarkDir =-1
	iben=0

	ALLOCATE(WK2(ME))

	DO 1000 IBOUN=1,NBOUND
		IPair = prdData(IBOUN,1)
		IF(IPair.EQ.0.OR.MFind(IBOUN).EQ.1) CYCLE
		
		WRITE(IUT0,*) '   Boundary :', IBOUN

		!stack wall Element
		CALL StackWallPrd(IUT0,XYZ,WElem,MP,NP,NWT,N2D,N3D,NDOF,
     &			   FanWei,WNV,NPBW,NB,0.8d0,NBOUND,IPair,prdData)


		!find wall-face pair
		DO IW=1,NWT
			JB=WElem(WELEM(1,IW)+2,IW)
			IF(JB.NE.IBOUN) CYCLE
			
			IIT = NE+IW

			V0(1)=CTR(1,IIT)+prdData(JB,3)			!x-offset
			V0(2)=CTR(2,IIT)+prdData(JB,4)			!y-offset
			V0(3)=CTR(3,IIT)+prdData(JB,5)			!z-offset
			
			Dist0 = 10.0*wholedist
		
			Istime=0
			
200			Istime= Istime+1
			IDOF=prdData(JB,2)
		
			DV = 0.d0
			DV(IDOF) = MarkDir*Dist0

c			DV(1) = MarkDir*Dist0*WNV(1,IW)
c			DV(2) = MarkDir*Dist0*WNV(2,IW)
c			DV(3) = MarkDir*Dist0*WNV(3,IW)
			
			VT(1) = V0(1) + DV(1)
			VT(2) = V0(2) + DV(2)
			VT(3) = V0(3) + DV(3)
			
			CALL GetSolidAngle(Seta,Eta,VT(1),VT(2),VT(3),
     &						   V0(1),V0(2),V0(3))

			!looking for the conterpart in periodic wall surface
			
			IDWall = 0
			CALL FindArrivalWall(DistW,Dist0,V0,VT,VW,DV,Seta,Eta,
     &		  IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NWT,iray,ierr,iben,
     &		  XYZ,FanWei)
			
			
			IF(IDWall.LE.0.AND.Istime.EQ.1) THEN
				MarkDir = -MarkDir
				GOTO 200
			END IF

			IF(IDWall.LE.0) THEN
				WRITE(IUT0,*) '   Error Abort in GetPrdPair()::'
				WRITE(IUT0,*) '   Find no periodic conterpart !!!'
				WRITE(IUT0,'(3x,a,I3,a,I6)') 'Boun_No.:',IBOUN, 
     &				'IW=',IW
				WRITE(IUT0,'(3x,a,3F12.3)') 
     &				'CTR=',CTR(1,IIT),CTR(2,IIT),CTR(3,IIT)
				STOP
			END IF

			WGrayDeg(IW)= -IDWall
			WGrayDeg(IDWall)= -IW

			!WRITE(IUT0,*) iw,idwall, ctr(:,ne+iw),ctr(:,ne+idwall)

		END DO

		MFind(IBOUN)=1
		MFind(-IPair)=1
		DEALLOCATE( WBlockNum,WBlockIndex,WBlockID)
1000	CONTINUE

	
	DEALLOCATE(WK2)


	RETURN

!-------> 

	END





C======================================================================C
C                                                                      C
C Subroutine Quadri4ShapeFun()								         C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2004/07/05                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE Quadri4ShapeFun(NODE,N3D,POINT,fun,NDOF)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D)
	REAL*8  V1(3),V2(3),V3(3),V4(3),POINT(3)
	REAL*8	S1(3),S2(3),fun(4),kesai,eta

C		THIS Subroutine IS USED TO GET THE SHAPE FUNCTION OF A quadrilateral.

	
							!THE VECTICES ARRANGEMENT:
							!	
							!     4 --------- 3 
							!     |		|	  |
							!	  |	----P---- |
							!     |		|	  |
							!     1 --------- 2
							!			  -->coef

	DO 10 I = 1,3
		V1(I) = NODE(I,1)
		V2(I) = NODE(I,2)
		V3(I) = NODE(I,3)
		V4(I) = NODE(I,4)
10	CONTINUE

	kesai = SiBian4(V1,V2,V4,V3,POINT)
	eta = SiBian4(V1,V4,V2,V3,POINT)


	fun(1) = (1-kesai)*(1-eta)
	fun(2) = kesai*(1-eta)
	fun(3) = kesai*eta
	fun(4) = (1-kesai)*eta


	RETURN

!----->
		
	END






C======================================================================C
C                                                                      C
C Function SiBian4()										             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/18                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	REAL*8 FUNCTION SiBian4(V1,V2,V3,V4,P)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  V1(3),V2(3),V3(3),V4(3),P(3)
	REAL*8	S1(3),S2(3)

C		THIS FUNCTION IS USED TO GET THE SHAPE FUNCTION OF A quadrangle.
C		THE RETURN VALUE IS THE COEFFICIENT OF THE INTERPOLATION RELATION
C		BETWEEN THE LINE 1-->2 AND 3-->4
C	    THE COEF IS DETERMINED BY CALCULATING THE ANGLE BETWEEN P-S1 AND S2-S1,
C		IF THE SINE IS ZERO, THEN P LIES IN S1-->S2.
C		REMINDING THAT 1-->S1 = COEF*(1-->2) [THE SAME FOR S2], WE THEN GET THE
C		VALUE COEF.

	
							!THE VECTICES ARRANGEMENT:
							!	
							!     3 ----S2--- 4 
							!     		|	 
							!	  		P	
							!     		|	 
							!     1 ----S1--- 2
							!			  -->coef


	VCOS = 0
	VSIN = 1
	COEF1 = 1
	COEF0 = 0
	cleftmax = 0
	crightmin = 1
	

	
	
	CALL AreaVector(V1(1),V1(2),V1(3),V3(1),V3(2),V3(3),
     $			P(1),P(2),P(3),vx0,vy0,vz0,AAA)
	CALL AreaVector(V2(1),V2(2),V2(3),V4(1),V4(2),V4(3),
     $			P(1),P(2),P(3),vx1,vy1,vz1,AAA)
c	CALL SINVALUE(V1(1)-P(1),V1(2)-P(2),V1(3)-P(3),
c     $			  V3(1)-P(1),V3(2)-P(2),V3(3)-P(3),VSIN0)
c	CALL SINVALUE(V2(1)-P(1),V2(2)-P(2),V2(3)-P(3),
c     $			  V4(1)-P(1),V4(2)-P(2),V4(3)-P(3),VSIN1)
	CALL COSVALUE(V1(1)-P(1),V1(2)-P(2),V1(3)-P(3),
     $			  V3(1)-P(1),V3(2)-P(2),V3(3)-P(3),VCOS0)
	CALL COSVALUE(V2(1)-P(1),V2(2)-P(2),V2(3)-P(3),
     $			  V4(1)-P(1),V4(2)-P(2),V4(3)-P(3),VCOS1)
	CALL COSVALUE(vx1,vy1,vz1,vx0,vy0,vz0,VCOS)
	
	!write(*,*) 'vcos', vcos0,vcos1

	IF(1-abs(VCOS0).LE.1.0E-6) THEN
		SiBian4 = 0
		RETURN
	END IF

	IF(1-abs(VCOS1).LE.1.0E-6) THEN
		SiBian4 = 1
		RETURN
	END IF

	IF(VCOS.EQ.1.0) THEN
		WRITE(*,*) 'ERROR in SiBian4(): the shapefun is < 0 or > 1' 
		STOP
	END IF
	
	IF(VCOS0.EQ.0.0) THEN
		DO 10 I= 1,3
			S1(I) = V1(I)-0.5*(V2(I)-V1(I))
			S2(I) = V3(I)-0.5*(V4(I)-V3(I))
10		CONTINUE
		CALL COSVALUE(S1(1)-P(1),S1(2)-P(2),S1(3)-P(3),
     $		  S2(1)-P(1),S2(2)-P(2),S2(3)-P(3),VCOS0)
      END IF	
	vx = vx0
	vy = vy0
	vz = vz0

	iter = 0
100	iter = iter + 1
	
	vxp = vx
	vyp = vy
	vzp = vz
	COEF = COEF1

	DO 20 I= 1,3
		S1(I) = V1(I)+COEF1*(V2(I)-V1(I))
		S2(I) = V3(I)+COEF1*(V4(I)-V3(I))
20	CONTINUE
	
	CALL AreaVector(S1(1),S1(2),S1(3),S2(1),S2(2),S2(3),
     $			P(1),P(2),P(3),vx,vy,vz,AAA)

	sx1 = S1(1)-P(1)
	sy1 = S1(2)-P(2)
	sz1 = S1(3)-P(3)
	sx2 = S2(1)-P(1)
	sy2 = S2(2)-P(2)
	sz2 = S2(3)-P(3)
	sx0 = S1(1)-S2(1)
	sy0 = S1(2)-S2(2)
	sz0 = S1(3)-S2(3)
	r0 = sx0*sx0 + sy0*sy0 + sz0*sz0
	r1 = (sx1*sx1 + sy1*sy1 + sz1*sz1)/r0
	r2 = (sx2*sx2 + sy2*sy2 + sz2*sz2)/r0

	IF(r1.LE.1E-3.OR.r2.LE.1.0E-3) THEN
		SiBian4 = COEF 
		RETURN
	END IF


	CALL SinValue(S1(1)-P(1),S1(2)-P(2),S1(3)-P(3),
     $			S2(1)-P(1),S2(2)-P(2),S2(3)-P(3),VSIN)	
	
	IF(VSIN.GT.1.0E-4) THEN
			
		CALL COSVALUE(vxp,vyp,vzp,vx,vy,vz,VCOS)
		CALL COSVALUE(vx0,vy0,vz0,vx,vy,vz,VCOS0)
		
		IF(VCOS.GT.0) THEN	!both are in one side

			IF(VCOS0.GT.0) THEN		!both are in LEFT side
									!0---c1---c0---p----1
									!0---c0---c1---p----1
				cleftmax = max(cleftmax,coef0,coef1)
				COEF = (COEF0+crightmin)*0.5
				
			ELSE					!both are in RIGHT side
									!0---p---c1---c0----1
									!0---p---c0---c1----1
				crightmin = min(crightmin,coef0,coef1)
				COEF = (COEF0+cleftmax)*0.5
			END IF
			
		ELSE 
									!across the P
									!0---c0---p---c1----1
									!0---c1---p---c0----1
				COEF = (COEF1+COEF0)*0.5
	
		END IF
		
		
		COEF0 = COEF1
		COEF1 = COEF
	
		!write(*,*) 'icase = ',icase,'vsin=',vsin
		!write(*,*) 'coef=',coef0,coef1,coef
	!pause
	
		GOTO 100

	END IF
	
	SiBian4 = COEF

	RETURN

!----->
		
	END
C======================================================================C
	SUBROUTINE Qua4Pos(NODE,N3D,POINT,fun,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D),POINT(3),fun(4)

C		THIS FUNCTION IS USED TO GET THE POSITION IN A quadrangle BY SHAPE FUNCTION.

	
	COEF1 = fun(1)
	COEF2 = fun(2)

	fun(1) = (1-COEF1)*(1-COEF2)
	fun(2) = COEF1*    (1-COEF2)
	fun(3) = COEF1*    COEF2
	fun(4) = (1-COEF1)*COEF2


	DO 50 I = 1, 3
	  POINT(I) = 0.0
	  DO 50 J = 1,4
		 POINT(I)=POINT(I)+fun(J)*NODE(I,J)
50	CONTINUE

	RETURN

!----->
		
	END






C======================================================================C
C                                                                      C
C Subroutine Tri3ShapeFun()								             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/18                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE Tri3ShapeFun(NODE,N3D,POINT,fun,NDOF)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D)
	REAL*8  V1(3),V2(3),V3(3),POINT(3),fun(3)
	
	

C		THIS Subroutine IS USED TO GET THE SHAPE FUNCTION OF A Triangle.
C		THE RETURN VALUE IS THE SHAPE FUNCTION FOR AREA RATIO OF P-V2-V3 TO VTHIS-V2-V3

	
							!THE VECTICES ARRANGEMENT:
							!	                 V3 
							!		           // |
							!		         / /  |
							!		       /  /   |
							!		     /	 P	  |
							!		   /	 `    |
							!	     /	       `  |
							!     V1 -----------`V2 
							!			  

	DO 10 I = 1,3
		V1(I) = NODE(I,1)
		V2(I) = NODE(I,2)
		V3(I) = NODE(I,3)
		
10	CONTINUE
	
	CALL AreaVector(V1(1),V1(2),V1(3),
     &				V2(1),V2(2),V2(3),
     &				V3(1),V3(2),V3(3),
     &				tvx,tvy,tvz,AreaTotal)

	CALL AreaVector(POINT(1),POINT(2),POINT(3),
     &				V2(1),V2(2),V2(3),
     &				V3(1),V3(2),V3(3),
     &				tvx,tvy,tvz,Area)
	
	fun(1) = Area/AreaTotal

	CALL AreaVector(POINT(1),POINT(2),POINT(3),
     &				V1(1),V1(2),V1(3),
     &				V3(1),V3(2),V3(3),
     &				tvx,tvy,tvz,Area)
	
	fun(2) = Area/AreaTotal

	CALL AreaVector(POINT(1),POINT(2),POINT(3),
     &				V1(1),V1(2),V1(3),
     &				V2(1),V2(2),V2(3),
     &				tvx,tvy,tvz,Area)
	
	fun(3) = Area/AreaTotal


	
c	fun(1) = SanJiao3(V1,V2,V3,POINT)
c	fun(2) = SanJiao3(V2,V1,V3,POINT)
c	fun(3) = SanJiao3(V3,V1,V2,POINT)


	RETURN

!----->
		
	END







C======================================================================C
C                                                                      C
C Function SanJiao3()										             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2004/07/05                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	REAL*8 FUNCTION SanJiao3(VTHIS,V2,V3,POINT)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  VTHIS(3),V2(3),V3(3),POINT(3)
	

C		THIS FUNCTION IS USED TO GET THE SHAPE FUNCTION OF A Triangle.
C		THE RETURN VALUE IS THE AREA RATIO OF P-V2-V3 TO VTHIS-V2-V3

	
							!THE VECTICES ARRANGEMENT:
							!	                 V3 
							!		           // |
							!		         / /  |
							!		       /  /   |
							!		     /	 P	  |
							!		   /	 `    |
							!	     /	       `  |
							!     VTHIS ---------`V2 
							!			  

	CALL AreaVector(VTHIS(1),VTHIS(2),VTHIS(3),
     &				V2(1),V2(2),V2(3),
     &				V3(1),V3(2),V3(3),
     &				tvx,tvy,tvz,AreaTotal)

	CALL AreaVector(Point(1),Point(2),Point(3),
     &				V2(1),V2(2),V2(3),
     &				V3(1),V3(2),V3(3),
     &				tvx,tvy,tvz,Area)
C
	
	SanJiao3 = Area/AreaTotal

	RETURN

!----->
		
	END
C======================================================================C
	SUBROUTINE Tri3Pos(NODE,N3D,POINT,fun,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D), POINT(3),fun(3)
	
	

C		THIS Subroutine IS USED TO GET THE POSITION in A Triangle by SHAPE FUNCTION.

	
							!THE VECTICES ARRANGEMENT:
							!	                 V3 
							!		           // |
							!		         / /  |
							!		       /  /   |
							!		     /	 P	  |
							!		   /	 `    |
							!	     /	       `  |
							!     V1 -----------`V2 
							!			  
	fun = fun/2.0
	fun(3) = 1.0-fun(1)-fun(2)

	DO 20 I = 1,3
		POINT(I) = fun(1)*NODE(I,1)+fun(2)*NODE(I,2)
     &			+fun(3)*NODE(I,3)
20	CONTINUE


	RETURN

!----->
		
	END







C======================================================================C
C                                                                      C
C Subroutine Tet4ShapeFun()								             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/18                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE Tet4ShapeFun(NODE,N3D,POINT,fun,NDOF)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D)
	REAL*8  V1(3),V2(3),V3(3),V4(3),POINT(3),fun(4)
	
	

C		THIS Subroutine IS USED TO GET THE SHAPE FUNCTION OF A Tetrahedral.
C
C														  

	DO 10 I = 1,3
		V1(I) = NODE(I,1)
		V2(I) = NODE(I,2)
		V3(I) = NODE(I,3)
		V4(I) = NODE(I,4)
10	CONTINUE
	
	CALL TetVolume(V1,V2,V3,V4,VTotal)
	
	CALL TetVolume(POINT,V2,V3,V4,VPart)
	fun(1) = VPart/VTotal

	CALL TetVolume(POINT,V1,V3,V4,VPart)
	fun(2) = VPart/VTotal

	CALL TetVolume(POINT,V1,V2,V4,VPart)
	fun(3) = VPart/VTotal

	CALL TetVolume(POINT,V1,V2,V3,VPart)
	fun(4) = VPart/VTotal


	sumfun = fun(1)+fun(2)+fun(3)+fun(4)

	IF(sumfun.GT.1.05) THEN
		WRITE(*,*) 'WARNING,Sumary of shape function > 1'
		WRITE(*,*) 'The point may be out of element'
	END IF


	RETURN

!----->
		
	END






C======================================================================C
C                                                                      C
C Subroutine Pyramid5ShapeFun()							             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/19                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE Pyramid5ShapeFun(NODE,N3D,POINT,fun,ElemID,NNPS,Local,
     &							 NFACE,N2D,NDOF,NKIND)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D),kesai,eta,S1(3),S2(3),S3(3),S4(3)
	REAL*8  VPex(3),V1(3),V2(3),V3(3),V4(3),POINT(3),fun(5)
	INTEGER	Local(NKIND,NFACE,N2D),IndexID(5),NNPS(NKIND,NFACE),
     &	ElemID,IdPex
	

C		THIS Subroutine IS USED TO GET THE SHAPE FUNCTION OF A Pyramid.
C
C														  


	DO 10 JFACE = 1,5					!bottom quadrileteral

		IF(NNPS(ElemID,JFACE).EQ.4) THEN
			IndexID(1) = Local(ElemID,JFACE,1)
			IndexID(2) = Local(ElemID,JFACE,2)
			IndexID(3) = Local(ElemID,JFACE,3)
			IndexID(4) = Local(ElemID,JFACE,4)
			GOTO 15
		END IF
10	CONTINUE
	
15	DO 20 I = 1,5
		IF(IndexID(I).EQ.0) IdPex = I
20	CONTINUE
	!write(*,*) IdPex
	DO 30 I = 1,4
		!write(*,*) IndexID(I)
		IF(IndexID(I).GT.5) THEN
			WRITE(*,*) 'ERROR in Pyramid3DShapeFun, IP > 5'
			STOP
		END IF
		
30	CONTINUE
	

	VPex(1) = NODE(1,IdPex)
	VPex(2) = NODE(2,IdPex)
	VPex(3) = NODE(3,IdPex) 

	DO 40 I = 1,3

		V1(I) = NODE(I,IndexID(1))
		V2(I) = NODE(I,IndexID(2))
		V3(I) = NODE(I,IndexID(3))
		V4(I) = NODE(I,IndexID(4))
40	CONTINUE
	
	CALL TetVolume(VPex,V1,V2,V3,VTotal)
	CALL TetVolume(POINT,V1,V2,V3,VPart)
	fun(IdPex) = VPart/VTotal
	
	
	!
	DO 50 I = 1,3
		S1(I) = V1(I) + (VPex(I)-V1(I))*fun(IdPex)
		S2(I) = V2(I) + (VPex(I)-V2(I))*fun(IdPex)
		S3(I) = V3(I) + (VPex(I)-V3(I))*fun(IdPex)
		S4(I) = V4(I) + (VPex(I)-V4(I))*fun(IdPex)
	!write(*,*) S1(I),S2(I),S3(I),S4(I),point(I)
50	CONTINUE

	kesai = SiBian4(S1,S2,S4,S3,POINT)
	eta = SiBian4(S1,S4,S2,S3,POINT)

	!write(*,*) 'kesai', kesai,eta

	fun(IndexID(1)) = (1-kesai)*(1-eta)*(1-fun(IdPex))
	fun(IndexID(2)) = kesai*(1-eta)*(1-fun(IdPex))
	fun(IndexID(3)) = kesai*eta*(1-fun(IdPex))
	fun(IndexID(4)) = (1-kesai)*eta*(1-fun(IdPex))

	sumfun = fun(1)+fun(2)+fun(3)+fun(4)+fun(5)
	IF(sumfun.GT.1.05) THEN
		WRITE(*,*) 'WARNING,Sumary of shape function > 1'
		WRITE(*,*) 'The point may be out of element'
	END IF

	RETURN

!----->
		
	END




C======================================================================C
C                                                                      C
C Subroutine Prism6ShapeFun()								             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/19                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE Prism6ShapeFun(NODE,N3D,POINT,fun,ElemID,NNPS,Local,
     &						   NFACE,N2D,NDOF,NKIND)

	IMPLICIT REAL*8(A-H,O-Z)
	INTEGER ElemID
	REAL*8  NODE(NDOF,N3D),kesai,eta, S1(3),S2(3),S3(3)
	REAL*8  V1(3),V2(3),V3(3),V4(3),V5(3),V6(3),POINT(3),fun(6)
	INTEGER	Local(NKIND,NFACE,N2D),IndexID(6),NNPS(NKIND,NFACE)
	

C		THIS Subroutine IS USED TO GET THE SHAPE FUNCTION OF A Prism.
C
C	
C											4--------5		--- COEF =1
C											|`	   / |
C											|  `6/   |
C											|   |    |
C											1---|----2		--- COEF =0									  
C											 ` 	|  /
C											   `3/
C
C

	
	


	DO 10 JFACE = 1,5					!first triangle

		IF(NNPS(ElemID,JFACE).EQ.3) THEN
			IndexID(1) = Local(ElemID,JFACE,1)
			IndexID(2) = Local(ElemID,JFACE,2)
			IndexID(3) = Local(ElemID,JFACE,3)
			GOTO 15
		END IF
10	CONTINUE
	


15	DO 30 JFACE = 1,5
		IF(NNPS(ElemID,JFACE).EQ.4) THEN		!quadrileterals
			DO 25 KNODE = 1,4
			  ID = Local(ElemID,JFACE,KNODE)
			  DO 23 MID = 1,3 
				IF(ID.EQ.IndexID(MID)) THEN		!node in the first triangle
				  NextID = Local(ElemID,JFACE,MOD(KNODE+1,4))
				  Mark = 0						
				  DO 20 I = 1,3
					IF(NextID.EQ.IndexID(I)) Mark = 1	
												!next node is also in the first triangle
20				  CONTINUE
				  IF(Mark.EQ.0) IndexID(MID+3) = NextID
									!if the next node is not in the first triangle,
									!it must be in the second triangle, and
									! MID-->NextID is in an edge
				END IF
23			  CONTINUE
25			CONTINUE
		END IF
30	CONTINUE

	

	DO 40 I=1,6
		
	
		IF(IndexID(I).LE.0) THEN
			WRITE(*,*) 'ERROR in Prism6ShapeFun: Node ID is not found'
			WRITE(*,*) 'Please check the Local() setting'
			STOP
		END IF
40	CONTINUE
	
	

	DO 50 I = 1,3

		V1(I) = NODE(I,IndexID(1))
		V2(I) = NODE(I,IndexID(2))
		V3(I) = NODE(I,IndexID(3))
		V4(I) = NODE(I,IndexID(4))
		V5(I) = NODE(I,IndexID(5))
		V6(I) = NODE(I,IndexID(6))
	
50	CONTINUE



	COEF = GetPointSect(V1,V2,V3,V4,V5,V6,POINT)

	!write(*,*) 'COEF=',COEF

	DO 100 I= 1,3
		S1(I) = V1(I)+COEF*(V4(I)-V1(I))
		S2(I) = V2(I)+COEF*(V5(I)-V2(I))
		S3(I) = V3(I)+COEF*(V6(I)-V3(I))
100	CONTINUE

	
200	CALL AreaVector(S1(1),S1(2),S1(3),
     &				S2(1),S2(2),S2(3),
     &				S3(1),S3(2),S3(3),
     &				tvx,tvy,tvz,AreaTotal)

	CALL AreaVector(POINT(1),POINT(2),POINT(3),
     &				S2(1),S2(2),S2(3),
     &				S3(1),S3(2),S3(3),
     &				tvx,tvy,tvz,Area)

	ratio = Area/AreaTotal
	fun(IndexID(1))= (1-COEF)*ratio
	fun(IndexID(4))= COEF*ratio
	
	!write(*,*) 'ratio1=',ratio


	CALL AreaVector(POINT(1),POINT(2),POINT(3),
     &				S1(1),S1(2),S1(3),
     &				S3(1),S3(2),S3(3),
     &				tvx,tvy,tvz,Area)

	ratio = Area/AreaTotal
	fun(IndexID(2))= (1-COEF)*ratio
	fun(IndexID(5))= COEF*ratio

	!write(*,*) 'ratio2=',ratio

	CALL AreaVector(POINT(1),POINT(2),POINT(3),
     &				S1(1),S1(2),S1(3),
     &				S2(1),S2(2),S2(3),
     &				tvx,tvy,tvz,Area)

	ratio = Area/AreaTotal
	fun(IndexID(3))= (1-COEF)*ratio
	fun(IndexID(6))= COEF*ratio

	
	!write(*,*) 'ratio3=',ratio
	!sumfun = fun(1)+fun(2)+fun(3)+fun(4)+fun(5)+fun(6)

	IF(sum(fun).GT.1.01) THEN
		WRITE(*,*) 'WARNING in Prism6Shapefun: '
		WRITE(*,*) '  Sumary of shape function > 1'
		WRITE(*,*) '  The point may be out of element'
	END IF

	RETURN

!----->
		
	END









C======================================================================C
C                                                                      C
C SUBROUTINE IsInFace()									             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/19                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE IsInFace(V1,V2,V3,P,value)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  V1(3),V2(3),V3(3),P(3),ABCD(4)
	

C		THIS FUNCTION IS USED TO Build a surface Ax+By+C+D=0 occupied by V1,V2,V3
C		AND get the excess value of the point P
C		If the P is in the face, the value should be zero.

	
	ABCD(1) = (V2(2)-V1(2))*(V3(3)-V1(3))-(V3(2)-V1(2))*(V2(3)-V1(3))
	ABCD(2) = (V3(1)-V1(1))*(V2(3)-V1(3))-(V2(1)-V1(1))*(V3(3)-V1(3))
	ABCD(3) = (V2(1)-V1(1))*(V3(2)-V1(2))-(V3(1)-V1(1))*(V2(2)-V1(2))
	ABCD(4) = -(ABCD(1)*V1(1)+ABCD(2)*V1(2)+ABCD(3)*V1(3))

C
	vmax = 0.0
	DO 10 I = 1,4
		vmax = max(abs(ABCD(I)),vmax)
10	CONTINUE

	value = (ABCD(1)*P(1)+ABCD(2)*P(2)+ABCD(3)*P(3)+ABCD(4))/vmax
	
	!write(*,*) 'abcd=',ABCD(1),ABCD(2),ABCD(3),ABCD(4),value

	RETURN

!----->
	
	END






C======================================================================C
C                                                                      C
C FUNCTION GetPointSect()									             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/19                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	REAL*8 FUNCTION GetPointSect(V1,V2,V3,V4,V5,V6,POINT)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  S1(3),S2(3),S3(3)
	REAL*8  V1(3),V2(3),V3(3),V4(3),V5(3),V6(3),POINT(3)


		
C		THIS function is to get a section containing P,
C		and returns a coefficient COEF as shape function
C		
C
C	
C											4--------5		--- COEF =1
C											|`	   / |
C											|  `6/   |
C											|   |    |
C											S1--|---S2		--- COEF								  
C											|`   P / |
C											| ` S3/  |
C											|   |    |
C											1---|----2		--- COEF =0									  
C											 ` 	|  /
C											   `3/
C
C



	CALL IsInFace(V1,V2,V3,Point,value1)
	CALL IsInFace(V4,V5,V6,Point,value2)
	vmax = max(abs(value1),abs(value2))

	IF(abs(value1)/vmax.LE.1.E-2) THEN
		GetPointSect = 0			!is in bottom face
		RETURN
	END IF
	
	IF(abs(value2)/vmax.LE.1.E-2) THEN
		GetPointSect = 1			!is in top face
		RETURN
	END IF

	IF(value1*value2.GT.1.0E-1) THEN
		WRITE(*,*) 'WARNING in GetPointSect(): Point is out'
		write(*,*) value1,value2
		!STOP
	END IF

	
	
	clowmax = 0
	chighmin = 1
	
		
	COEF0 = -0.1
	DO 10 I= 1,3
		S1(I) = V1(I)+COEF0*(V4(I)-V1(I))
		S2(I) = V2(I)+COEF0*(V5(I)-V2(I))
		S3(I) = V3(I)+COEF0*(V6(I)-V3(I))
10	CONTINUE
	CALL IsInFace(S1,S2,S3,Point,value1)
	value = value1

	COEF1 = 1.1
	iter = 0
100	iter = iter + 1
	
	!IF((iter.GT.20).AND.(iter/20.EQ.0)) COEF1 = 1 + 0.1*(iter/20)
	valueP=value
	COEF = COEF1
	value = value1
	IF(iter.GT.100) GOTO 1000
		
	DO 50 I= 1,3
		S1(I) = V1(I)+COEF1*(V4(I)-V1(I))
		S2(I) = V2(I)+COEF1*(V5(I)-V2(I))
		S3(I) = V3(I)+COEF1*(V6(I)-V3(I))
50	CONTINUE
	
	CALL IsInFace(S1,S2,S3,Point,value)
	

	IF(abs(value).GT.1.0E-4) THEN
		
		IF(value*valueP.GT.0) THEN		!both sections are in one side

			IF(value1*valueP.GT.0) THEN	 !section-face is too low
									
				clowmax = max(clowmax,coef0,coef1)
				COEF = (COEF0+chighmin)*0.5
				
			ELSE					!section-face is too high

				chighmin = min(chighmin,coef0,coef1)
				COEF = (COEF0+clowmax)*0.5
			END IF
			
		ELSE 
									!across the Point
									
				COEF = (COEF1+COEF0)*0.5
	
		END IF
		
		
		COEF0 = COEF1
		COEF1 = COEF
	
		!write(*,*) 'iter = ',iter,'value=',value,'coef=',coef
		!write(*,*) 'coef=',coef0,coef1,coef
	!pause
		

		GOTO 100

	END IF

1000	GetPointSect = COEF

	RETURN



	END 

	

	






C======================================================================C
C                                                                      C
C Subroutine Hex8ShapeFun()								             C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/01/19                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE Hex8ShapeFun(NODE,N3D,POINT,fun,ElemID,NNPS,Local,
     &						NFACE,N2D,NDOF,NKIND)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D),kesai,eta, S1(3),S2(3),S3(3)
	REAL*8  V1(3),V2(3),V3(3),V4(3),V5(3),V6(3),V7(3),V8(3),
     &		POINT(3),fun(8)
	INTEGER	Local(NKIND,NFACE,N2D),IndexID(8),NNPS(NKIND,NFACE),ElemID
	

C		THIS Subroutine IS USED TO GET THE SHAPE FUNCTION OF A hexahedron
C
C	
C											8--------7		
C										   /|	   / |
C										  5-------6  |
C										  |	|     |  |
C					y					  |	4-----|--3		
C										  |/	  | / 
C					|					  1-------2  
C					|
C				   /----- x
C				  /
C				 z
							!first face (bottom face)
	IndexID(1) = Local(ElemID,1,1)
	IndexID(2) = Local(ElemID,1,2)
	IndexID(3) = Local(ElemID,1,3)
	IndexID(4) = Local(ElemID,1,4)

	
	
15	DO 30 JFACE = 2,6
		DO 25 KNODE = 1,4
		  ID = Local(ElemID,JFACE,KNODE)
		  DO 23 MID = 1,4 
			IF(ID.EQ.IndexID(MID)) THEN		!node in the first face
			  NextID = Local(ElemID,JFACE,MOD(KNODE+1,4))
			  Mark = 0		
			  DO 20 I = 1,4
				IF(NextID.EQ.IndexID(I)) Mark = 1	
											!next node is also in the first face
20			  CONTINUE
			  IF(Mark.EQ.0) IndexID(MID+4) = NextID
								!if the next node is not in the first face,
								!it must be in the top face, and
								! MID-->NextID is in an edge
			END IF
23		  CONTINUE
25		CONTINUE
30	CONTINUE

	

	DO 40 I=1,8
		
		IF(IndexID(I).LE.0) THEN
			WRITE(*,*) 'ERROR in Hex8ShapeFun: Node ID is not found'
			WRITE(*,*) 'Please check the Local() setting'
			STOP
		END IF
40	CONTINUE

	

	DO 50 I = 1,3

		V1(I) = NODE(I,IndexID(1))
		V2(I) = NODE(I,IndexID(2))
		V3(I) = NODE(I,IndexID(3))
		V4(I) = NODE(I,IndexID(4))
		V5(I) = NODE(I,IndexID(5))
		V6(I) = NODE(I,IndexID(6))
		V7(I) = NODE(I,IndexID(7))
		V8(I) = NODE(I,IndexID(8))
	
50	CONTINUE


	COEF1 = GetPointSect(V1,V5,V8,V2,V6,V7,POINT)		!X--dir
	!write(*,*) 'coef1=',coef1
	COEF2 = GetPointSect(V1,V2,V3,V5,V6,V7,POINT)		!Y--dir
	!write(*,*) 'coef2=',coef2
	COEF3 = GetPointSect(V1,V2,V6,V4,V3,V7,POINT)		!Z--dir
	!write(*,*) 'coef2=',coef2

	fun(IndexID(1)) = (1-COEF1)*(1-COEF2)*(1-COEF3)
	fun(IndexID(2)) = COEF1*    (1-COEF2)*(1-COEF3)
	fun(IndexID(3)) = COEF1*    (1-COEF2)*    COEF3
	fun(IndexID(4)) = (1-COEF1)*(1-COEF2)*    COEF3

	fun(IndexID(5)) = (1-COEF1)*    COEF2*(1-COEF3)
	fun(IndexID(6)) = COEF1*        COEF2*(1-COEF3)
	fun(IndexID(7)) = COEF1*        COEF2*    COEF3
	fun(IndexID(8)) = (1-COEF1)*    COEF2*    COEF3


	!sumfun = fun(1)+fun(2)+fun(3)+fun(4)+fun(5)+fun(6)

	IF(sum(fun).GT.1.01) THEN
		WRITE(*,*) 'WARNING in Hex8Shapefun: '
		WRITE(*,*) '  Sumary of shape function > 1'
		WRITE(*,*) '  The point may be out of element'
	END IF

	RETURN


!----->
		
	END

C======================================================================C
	SUBROUTINE Tet4Pos(NODE,N3D,V0,fun,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D)
	REAL*8  V1(3),V2(3),V3(3),V4(3),V0(3),fun(4)
	
	

C		THIS Subroutine IS USED TO GET THE POSITION in the elem by SHAPE FUNCTION.
C
C														  

	
	fun = fun/3.0 
	fun(4) = 1.0-fun(1)-fun(2)-fun(3)

	DO 20 I = 1,3
		V0(I) = fun(1)*NODE(I,1)+fun(2)*NODE(I,2)
     &			+fun(3)*NODE(I,3)+fun(4)*NODE(I,4)
20	CONTINUE


	RETURN

!----->
		
	END

C======================================================================C
	SUBROUTINE Pyramid5Pos(V,N3D,POINT,fun,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  V(NDOF,N3D),kesai,eta, POINT(3),fun(5)

	

C		THIS Subroutine IS USED TO GET THE POSITION in the element by SHAPE FUNCTION.
C
C														  
	fun(5) =fun(1)
	kesai = fun(2)
	eta   = fun(3)

	fun(1) = (1-kesai)*(1-eta)*(1-fun(5))
	fun(2) = kesai*(1-eta)*(1-fun(5))
	fun(3) = kesai*eta*(1-fun(5))
	fun(4) = (1-kesai)*eta*(1-fun(5))




	DO 50 I =  1,3
		
	  POINT(I)=V(I,5)*fun(5)
     &		   +V(I,1)*fun(1)+V(I,2)*fun(2)+V(I,3)*fun(3)+V(I,4)*fun(4)

50	CONTINUE	


	RETURN

!----->
		
	END
C======================================================================C
	SUBROUTINE Pyramid5Arrange(NODE,V,ElemID,NNPS,Local,
     &						 NFACE,N2D,N3D,NKIND,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D),kesai,eta,S1(3),S2(3),S3(3),S4(3)
	REAL*8  V(3,N3D)
	INTEGER	Local(NKIND,NFACE,N2D),IndexID(5),NNPS(NKIND,NFACE),
     &		ElemID,IdPex
	

C		THIS Subroutine IS USED TO ReArrange the CONNECTIVITY OF A Pyramid Element.
C
C														  


	DO 10 JFACE = 1,5					!get bottom quadrileteral

		IF(NNPS(ElemID,JFACE).EQ.4) THEN
			IndexID(1) = Local(ElemID,JFACE,1)
			IndexID(2) = Local(ElemID,JFACE,2)
			IndexID(3) = Local(ElemID,JFACE,3)
			IndexID(4) = Local(ElemID,JFACE,4)
			GOTO 15
		END IF
10	CONTINUE
	
15	DO 20 I = 1,5
		IF(IndexID(I).EQ.0) IdPex = I
20	CONTINUE
	!write(*,*) IdPex
	DO 30 I = 1,4
		!write(*,*) IndexID(I)
		IF(IndexID(I).GT.5) THEN
			WRITE(*,*) 'ERROR in Pyramid5Arrange, IP > 5'
			STOP
		END IF
		
30	CONTINUE
	
	V(1,5) = NODE(1,IdPex)
	V(2,5) = NODE(2,IdPex)
	V(3,5) = NODE(3,IdPex) 

	DO 40 I = 1,3

		V(I,1) = NODE(I,IndexID(1))
		V(I,2) = NODE(I,IndexID(2))
		V(I,3) = NODE(I,IndexID(3))
		V(I,4) = NODE(I,IndexID(4))
40	CONTINUE
	

	RETURN

!----->
		
	END
C======================================================================C
	SUBROUTINE Prism6Pos(NODE,N3D,POINT,fun,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8  NODE(NDOF,N3D),kesai,eta, POINT(3),fun(6)
		

C		THIS Subroutine IS USED TO GET THE POSITION by SHAPE FUNCTION OF A Prism.
C
C	
C											4--------5		--- COEF =1
C											|`	   / |
C											|  `6/   |
C											|   |    |
C											1---|----2		--- COEF =0									  
C											 ` 	|  /
C											   `3/
C
C

	
	kesai = 0.5*fun(1)
	eta = 0.5*fun(2)
	COEF = fun(3)
	
	fun(1)=(1-kesai-eta)*(1-COEF)
	fun(2)=kesai*(1-COEF)
	fun(3)=eta*(1-COEF)
	fun(4)=(1-kesai-eta)*COEF
	fun(5)=kesai*COEF
	fun(6)=eta*COEF


	DO 60 I = 1,3
		POINT(I) = 0.0
		DO 60 J= 1,6

		POINT(I)=POINT(I)+fun(J)*NODE(I,J)

60	CONTINUE

	RETURN

!----->
		
	END

C======================================================================C
	SUBROUTINE Prism6Arrange(NODE,V,ElemID,NNPS,Local,
     &					 NFACE,N2D,N3D,NKIND,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	INTEGER ElemID
	REAL*8  NODE(NDOF,N3D),kesai,eta, S1(3),S2(3),S3(3)
	REAL*8  V(3,N3D),POINT(3),fun(6)
	INTEGER	Local(NKIND,NFACE,N2D),IndexID(6),NNPS(NKIND,NFACE)
	

C		THIS Subroutine IS USED TO ReArrange the CONNECTIVITY A Prism Element.
C
	

	DO 10 JFACE = 1,5					!first triangle

		IF(NNPS(ElemID,JFACE).EQ.3) THEN
			IndexID(1) = Local(ElemID,JFACE,1)
			IndexID(2) = Local(ElemID,JFACE,2)
			IndexID(3) = Local(ElemID,JFACE,3)
			GOTO 15
		END IF
10	CONTINUE
	

15	DO 30 JFACE = 1,5
		IF(NNPS(ElemID,JFACE).EQ.4) THEN		!quadrileterals
			DO 25 KNODE = 1,4
			  ID = Local(ElemID,JFACE,KNODE)
			  DO 23 MID = 1,3 
				IF(ID.EQ.IndexID(MID)) THEN		!node in the first triangle
				  KN1 =  KNODE+1
				  IF(KN1.GT.4) KN1=KN1-4
				  NextID = Local(ElemID,JFACE,KN1)
				  Mark = 0						
				  DO 20 I = 1,3
					IF(NextID.EQ.IndexID(I)) Mark = 1	
												!next node is also in the first triangle
20				  CONTINUE
				  IF(Mark.EQ.0) IndexID(MID+3) = NextID
									!if the next node is not in the first triangle,
									!it must be in the second triangle, and
									! MID-->NextID is in an edge
				END IF
23			  CONTINUE
25			CONTINUE
		END IF
30	CONTINUE


	DO 40 I=1,6
		
		IF(IndexID(I).LE.0) THEN
			WRITE(*,*) 'ERROR in Prism6Arrange: Node ID is not found'
			WRITE(*,*) 'Please check the Local() setting'
			STOP
		END IF
40	CONTINUE
		

	DO 50 I = 1,3

		V(I,1) = NODE(I,IndexID(1))
		V(I,2) = NODE(I,IndexID(2))
		V(I,3) = NODE(I,IndexID(3))
		V(I,4) = NODE(I,IndexID(4))
		V(I,5) = NODE(I,IndexID(5))
		V(I,6) = NODE(I,IndexID(6))
	
50	CONTINUE



	RETURN

!----->
		
	END

C======================================================================C
	SUBROUTINE Hex8Pos(NODE,N3D,POINT,fun,NDOF)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D), POINT(3),fun(8)
	
	

C		THIS Subroutine IS USED TO GET THE SHAPE FUNCTION OF A hexahedron
C
C	
C											8--------7		
C										   /|	   / |
C										  5-------6  |
C										  |	|     |  |
C					y					  |	4-----|--3		
C										  |/	  | / 
C					|					  1-------2  
C					|
C				   /----- x
C				  /
C				 z
							!first face (bottom face)


	COEF1 = fun(1)
	COEF2 = fun(2)
	COEF3 = fun(3)


	fun(1) = (1-COEF1)*(1-COEF2)*(1-COEF3)
	fun(2) = COEF1*    (1-COEF2)*(1-COEF3)
	fun(3) = COEF1*    (1-COEF2)*    COEF3
	fun(4) = (1-COEF1)*(1-COEF2)*    COEF3

	fun(5) = (1-COEF1)*    COEF2*(1-COEF3)
	fun(6) = COEF1*        COEF2*(1-COEF3)
	fun(7) = COEF1*        COEF2*    COEF3
	fun(8) = (1-COEF1)*    COEF2*    COEF3
	
	DO 50 I = 1, 3
	  POINT(I) = 0.0
	  DO 50 J = 1,8
		 POINT(I)=POINT(I)+fun(J)*NODE(I,J)
50	CONTINUE



	RETURN


!----->
		
	END
C======================================================================C
	SUBROUTINE Hex8Arrange(NODE,V,ElemID,NNPS,Local,
     &					 NFACE,N2D,N3D,NKIND,NDOF)
C======================================================================C

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(NDOF,N3D),V(NDOF,N3D),kesai,eta, S1(3),S2(3),S3(3)
	INTEGER	Local(NKIND,NFACE,N2D),IndexID(8),NNPS(NKIND,NFACE),ElemID
	

C		THIS Subroutine IS USED TO ReArrange the CONECTIVITY of A hexahedron Element
C
C	
C										8--------7		
C									   /|	   / |
C									  5-------6  |
C									  |	|     |  |
C				y					  |	4-----|--3		
C									  |/	  | / 
C				|					  1-------2  
C				|
C			   /----- x
C			  /
C			 z
							!first face (bottom face)



	IndexID(1) = Local(ElemID,1,1)
	IndexID(2) = Local(ElemID,1,2)
	IndexID(3) = Local(ElemID,1,3)
	IndexID(4) = Local(ElemID,1,4)
	
	

15	DO 30 JFACE = 2,6
		DO 25 KNODE = 1,4
		  ID = Local(ElemID,JFACE,KNODE)
		  !write(*,*) jface,knode,id

		  DO 23 MID = 1,4 
			IF(ID.EQ.IndexID(MID)) THEN		!node in the first face
			  KN1 =  KNODE+1
			  IF(KN1.GT.4) KN1=KN1-4
			  NextID = Local(ElemID,JFACE,KN1)
			  Mark = 0	
			  DO 20 I = 1,4
				IF(NextID.EQ.IndexID(I)) Mark = 1	
											!next node is also in the first face
20			  CONTINUE
	
			  IF(Mark.EQ.0) IndexID(MID+4) = NextID
								!if the next node is not in the first face,
								!it must be in the top face, and
								! MID-->NextID is in an edge
			END IF
23		  CONTINUE
25		CONTINUE
30	CONTINUE

	
	DO 40 I=1,8
		IF(IndexID(I).LE.0) THEN
			WRITE(*,*) 'ERROR in Hex8Arrange: Node ID is not found'
			WRITE(*,*) 'Please check the Local() setting'
			STOP
		END IF
40	CONTINUE


	DO 50 I = 1,3

		V(I,1) = NODE(I,IndexID(1))
		V(I,2) = NODE(I,IndexID(2))
		V(I,3) = NODE(I,IndexID(3))
		V(I,4) = NODE(I,IndexID(4))
		V(I,5) = NODE(I,IndexID(5))
		V(I,6) = NODE(I,IndexID(6))
		V(I,7) = NODE(I,IndexID(7))
		V(I,8) = NODE(I,IndexID(8))
	
50	CONTINUE


	RETURN


!----->
		
	END





C======================================================================C
C                                                                      C
C Subroutine TetVolume()										  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/1/17                      C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C

	

      SUBROUTINE TetVolume(P0,P1,P2,P3,VVV)
      IMPLICIT REAL*8(A-H,O-Z)

	
	REAL*8 P0(3),P1(3),P2(3),P3(3),VVV

C
C
C      CALCULATE a Tetrahedral'S VOLUME
C         ( 3-D ; SINGLE PRECISION, REFERING ELEMENT BY TYPE LIST )

	
		a1 = P0(1)-P1(1)
		a2 = P0(2)-P1(2)
		a3 = P0(3)-P1(3)
		
		b1 = P0(1)-P2(1)
		b2 = P0(2)-P2(2)
		b3 = P0(3)-P2(3)

		c1 = P0(1)-P3(1)
		c2 = P0(2)-P3(2)
		c3 = P0(3)-P3(3)


		AX  = (b2*c3-c2*b3)*a1
		AY  = (b1*c3-c1*b3)*a2
		AZ  = (b1*c2-c1*b2)*a3

          VVV = abs((AX - AY + AZ))/6.0

      RETURN


      END



	


C======================================================================C
C                                                                      C
C Subroutine PyramidVolume()									  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/09/02                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
      SUBROUTINE PyramidVolume(NODE,VVV)
      IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8 NODE(3,5),VVV
	REAL*8 V0(3),V1(3),V2(3),V3(3)
C
C
C      CALCULATE a Tetrahedral'S VOLUME
C         ( 3-D ; SINGLE PRECISION, REFERING ELEMENT BY TYPE LIST )

	
	
	V0(:) = NODE(:,1)

	V1(:) = NODE(:,2)
	V2(:) = NODE(:,3)
	V3(:) = NODE(:,4)
	CALL TetVolume(V0,V1,V2,V3,va)

	!V1(:) = NODE(:,2)
	V2(:) = NODE(:,5)
	!V3(:) = NODE(:,4)
	CALL TetVolume(V0,V1,V2,V3,vb)

	VVV =  vb+ va

      RETURN


      END







C======================================================================C
C                                                                      C
C Subroutine GaussTri3()										  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/08/26                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE	GaussTri3(value,NODE,N3D,V0,fun,alfa,area,iautoflag,
     &					  NVwo,NVta,icswo,icsta)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(3,N3D),NVwo(3),NVta(3),ds(4)
	REAL*8  VC(3), V0(3),fun(4),jacob


C	This subroutine gives the GAUSS approximation of the area-integration
C	of a Triangle. 4 quadrature points, the residence is O(h^5)
C
C
C
C
	!R0 = SQRT(0.75)			!circumradius, 0.5/sin(60deg)
	At = 0.75*SQRT(3.0)*0.75	!3/4* 3^0.5 * R0^2
	jacob = 2.0*Area			!Jacob coef. is constant
	
	VC(:) = (NODE(1,:)+NODE(2,:)+NODE(3,:))/3.0

	IF(iautoflag.EQ.1) THEN

		DO I = 1,3
		  CALL GetDistance(NODE(1,I),NODE(2,I),NODE(3,I),
     &		V0(1),V0(2),V0(3),ds(I))
		  
		END DO
		
		CALL GetDistance(VC(1),VC(2),VC(3),V0(1),V0(2),V0(3),ds(4))
	

		DO I =1,4
			IF(ds(I).GT.0.0)	THEN
				fun(I) = exp(-ds(I)*alfa)/(ds(I)*ds(I))
			ELSE
				fun(I) = 0
			END IF
			
		END DO

	END IF
	
	IF(icswo.EQ.1) THEN

		DO I = 1,3
		  CALL CosValue(NVwo(1),NVwo(2),NVwo(3),
     &				  V0(1)-NODE(1,I),V0(2)-NODE(2,I),V0(3)-NODE(3,I),
     &				  cs)
		  cs = max(cs,0.d0)
		  fun(I) = fun(I)*cs
		END DO

		CALL CosValue(NVwo(1),NVwo(2),NVwo(3),
     &				  V0(1)-VC(1),V0(2)-VC(2),V0(3)-VC(3),
     &				  cs)
		cs = max(cs,0.d0)
		fun(4) = fun(4)*cs

	END IF

	IF(icsta.EQ.1) THEN

		DO I = 1,3
		  CALL CosValue(NVta(1),NVta(2),NVta(3),
     &				  NODE(1,I)-V0(1),NODE(2,I)-V0(2),NODE(3,I)-V0(3),
     &				  cs)
		  cs = max(cs,0.d0)
		  fun(I) = fun(I)*cs
		END DO

		CALL CosValue(NVta(1),NVta(2),NVta(3),
     &				  VC(1)-V0(1),VC(2)-V0(2),VC(3)-V0(3),
     &				  cs)
		cs = max(cs,0.d0)
		fun(4) = fun(4)*cs

	END IF


	
	value = 1.0/12*(fun(1)+fun(2)+fun(3))
	value = value+ 0.75*fun(4)
	value = value*Jacob*At


	RETURN

	END








C======================================================================C
C                                                                      C
C Subroutine GaussSibian()									  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/08/26                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE	GaussSibian(value,NODE,N3D,V0,fun,alfa,iautoflag,
     &						NVwo,NVta,icswo,icsta)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(3,N3D),VC(3),A(4),ds(4),NVwo(3),NVta(3)
	REAL*8  V0(3), V1(3),V2(3),V3(3),V4(3),fun(4),jacob,S(4)
	REAL*8	PNODE(3,4)


C	This subroutine gives the GAUSS approximation of the area-integration
C	of a Quatriangle. 4 quadrature points, the residence is O(h^6)
C
C
C
C
	
	As = 1.0					!4*h^2, h=0.5
	ps = SQRT(1.0/3.0)*0.5


	S(1) = 0.5-Ps
	S(2) = 0.5-Ps

	CALL Qua4Pos(NODE,N3D,V1,S,3)
	PNode(:,1) = V1(:)

	S(1) = 0.5-Ps
	S(2) = 0.5+Ps
	CALL GetDistance(V2(1),V2(2),V2(3),V0(1),V0(2),V0(3),ds(2))
	PNode(:,2) = V2(:)

	S(1) = 0.5+Ps
	S(2) = 0.5-Ps
	CALL Qua4Pos(NODE,N3D,V3,S,3)
	PNode(:,3) = V3(:)

	S(1) = Ps+0.5
	S(2) = Ps+0.5
	CALL GetDistance(V4(1),V4(2),V4(3),V0(1),V0(2),V0(3),ds(4))
	PNode(:,4) = V4(:)


	IF(iautoflag.EQ.1) THEN
		
		CALL GetDistance(V1(1),V1(2),V1(3),V0(1),V0(2),V0(3),ds(1))
		CALL GetDistance(V2(1),V2(2),V2(3),V0(1),V0(2),V0(3),ds(2))
		CALL GetDistance(V3(1),V3(2),V3(3),V0(1),V0(2),V0(3),ds(3))
		CALL GetDistance(V4(1),V4(2),V4(3),V0(1),V0(2),V0(3),ds(4))
		
		vc(:) = 0.25*(v1(:)+v2(:)+v3(:)+v4(:))

		vc(:) = 0.25*(Node(:,1)+Node(:,2)+Node(:,3)+Node(:,4))

		DO I =1,4
			IF(ds(I).GT.0.0)	THEN
				fun(I) = exp(-ds(I)*alfa)/(ds(I)*ds(I))
			ELSE
				fun(I) = 0
			END IF
			
		END DO

	END IF


	
	if(icswo.EQ.1) THEN
		DO I = 1,4
		CALL CosValue(NVwo(1),NVwo(2),NVwo(3),V0(1)-PNode(1,I),
     &				  V0(2)-PNode(2,I),V0(3)-PNode(3,I),cs)
		cs = max(cs,0.d0)
		fun(I) = fun(I)*cs
		END DO

	END IF


	if(icsta.EQ.1) THEN
		DO I = 1,4
		CALL CosValue(NVta(1),NVta(2),NVta(3),PNode(1,I)-V0(1),
     &				  PNode(2,I)-V0(2),PNode(3,I)-V0(3),cs)
		cs = max(cs,0.d0)
		fun(I) = fun(I)*cs
		END DO

	END IF

	!center point
	VC(1) = (NODE(1,1)+NODE(1,2)+NODE(1,3)+NODE(1,4))*0.25
	VC(2) = (NODE(2,1)+NODE(2,2)+NODE(2,3)+NODE(2,4))*0.25
	VC(3) = (NODE(3,1)+NODE(3,2)+NODE(3,3)+NODE(3,4))*0.25
		
	DO I = 1,4
		I1 = I
		IA = (I1+1)		
		IF(IA.GT.4) IA=IA-4
		IB = (I1+3)
		IF(IB.GT.4) IB=IB-4

		PNODE(:,1) = NODE(:,I1)
		PNODE(:,2) = (NODE(:,I1)+NODE(:,IA))/2.0
		PNODE(:,3) = (NODE(:,I1)+NODE(:,IB))/2.0
		PNODE(:,4) = VC(:)

		CALL PolymFaceArea(PNODE,4,4,A(I))
			
	END DO

	
	value = fun(1)*A(1)+fun(2)*A(2)+fun(3)*A(3)+fun(4)*A(4)
	

	RETURN

	END 









C======================================================================C
C                                                                      C
C Subroutine GaussFaceQita()									  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/09/06                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE	GaussFaceQita(value,NODE,N3D,NP,V0,fun,alfa,iautoflag,
     &						NVwo,NVta,icswo,icsta)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(3,N3D),A(NP),ds(NP),NVwo(3),NVta(3)
	REAL*8  V0(3), fun(NP),jacob,S(4)
	REAL*8	PNODE(3,4),CTR(3),VC(3,NP)


C	This subroutine gives the GAUSS approximation of the area-integration
C	of a Quatriangle. 4 quadrature points, the residence is O(h^6)
C
C
C
C
	CTR = 0
	DO I= 1,NP
		CTR(:) = CTR(:) + NODE(:,I)
	END DO
	CTR(:) = CTR(:)/NP


	DO I = 1,NP
		I1 = I
		IA = (I1+1)		
		IF(IA.GT.NP) IA=IA-NP
		IB = (I1+NP-1)
		IF(IB.GT.NP) IB=IB-NP

		PNODE(:,1) = NODE(:,I1)
		PNODE(:,2) = (NODE(:,I1)+NODE(:,IA))/2.0
		PNODE(:,3) = (NODE(:,I1)+NODE(:,IB))/2.0
		PNODE(:,4) = CTR(:)

		VC(:,I) = (PNODE(:,1)+PNODE(:,2)+PNODE(:,3)+PNODE(:,4))/4.0

		CALL PolymFaceArea(PNODE,4,4,A(I))
			
	END DO


	IF(iautoflag.EQ.1) THEN
		
		DO I = 1, NP
			CALL GetDistance(VC(1,I),VC(2,I),VC(3,I),
     &						 V0(1),V0(2),V0(3),ds(I))
		
			IF(ds(I).GT.0.0)	THEN
				fun(I) = exp(-ds(I)*alfa)/(ds(I)*ds(I))
			ELSE
				fun(I) = 0
			END IF
			
		END DO

	END IF
	
	if(icswo.EQ.1) THEN
		DO I = 1,NP
		CALL CosValue(NVwo(1),NVwo(2),NVwo(3),V0(1)-VC(1,I),
     &				  V0(2)-VC(2,I),V0(3)-VC(3,I),cs)
		cs = max(cs,0.d0)
		fun(I) = fun(I)*cs
		END DO

	END IF


	if(icsta.EQ.1) THEN
		DO I = 1,4
		CALL CosValue(NVta(1),NVta(2),NVta(3),VC(1,I)-V0(1),
     &				  VC(2,I)-V0(2),VC(3,I)-V0(3),cs)
		cs = max(cs,0.d0)
		fun(I) = fun(I)*cs
		END DO

	END IF

	value = 0.0
	DO I=1,NP
		value = value + fun(I)*A(I)
	END DO

	RETURN

	END 






C======================================================================C
C                                                                      C
C Subroutine GaussTet4()										  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/08/26                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE	GaussTet4(value,NODE,N3D,V0,fun,
     &					  alfa,volume,NV,iautoflag,icosflag)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(3,N3D),ds(8),NV(3),cs(8)
	REAL*8  VC(3,4), V0(3),fun(8),jacob


C	This subroutine gives the GAUSS approximation of the volume-integration
C	of a Tetrahedron. 8 quadrature points, the residence is O(h^6)
C
C
C
C

	DO I = 1,4
		I1 = I
		I2 = (I1+1)		
		IF(I2.GT.4) I2=I2-4
		I3 = (I1+2)
		IF(I3.GT.4) I3=I3-4
		
		VC(1,I) = (NODE(1,I1)+NODE(2,I1)+NODE(3,I1))/3.0
		VC(2,I) = (NODE(1,I2)+NODE(2,I2)+NODE(3,I2))/3.0
		VC(3,I) = (NODE(1,I3)+NODE(2,I3)+NODE(3,I3))/3.0
		
	END DO
		

	IF(iautoflag.EQ.1) THEN

		DO I =1,4
		  CALL GetDistance(NODE(1,I),NODE(2,I),NODE(3,I),
     &		V0(1),V0(2),V0(3),ds(I))
		END DO
     		
		DO I = 1,4
			CALL GetDistance(VC(1,I),VC(2,I),VC(3,I),
     &		V0(1),V0(2),V0(3),ds(I+4))
	
		END DO
		

		DO I =1,8
			IF(ds(I).GT.0.0)	THEN
				fun(I) = exp(-ds(I)*alfa)/(ds(I)*ds(I))
			ELSE
				fun(I) = 0
			END IF
		END DO

	END IF
	

	IF(icosflag.EQ.1) THEN
		DO I =1,4
			CALL CosValue(V0(1)-NODE(1,I),V0(2)-NODE(2,I),
     &				V0(3)-NODE(3,I),NV(1),NV(2),NV(3),cs(I))
		END DO
     		
		DO I = 1,4
			CALL CosValue(V0(1)-VC(1,I),V0(2)-VC(2,I),V0(3)-VC(3,I),
     &					  NV(1),NV(2),NV(3),cs(I+4))
	
		END DO
		

		DO I =1,8
			cs(I) = max(cs(I),0.d0)
			fun(I) = fun(I)*cs(I)
		END DO

	END IF
			
	value = (fun(1)+fun(2)+fun(3)+fun(4))/40.0
	value = value+ (fun(5)+fun(6)+fun(7)+fun(8))*9.0/40.0
	value = value*Volume

	RETURN

	END








C======================================================================C
C                                                                      C
C Subroutine GaussHex8()										  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/08/26                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE	GaussHex8(value,NODE,N3D,V0,CTR,fun,CNTIVITY,N2D,
     &					  NFACE,NKIND,IKIND,alfa,NV,iautoflag,icosflag)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(3,N3D),FC(3,6),CTR(3),ds(6),V(6),NV(3)
	REAL*8  V0(3), fun(6),Pyramid(3,5),cs(6)
	INTEGER	CNTIVITY(NKIND,NFACE,N2D)


C	This subroutine gives the GAUSS approximation of the volume-integration
C	of a Hexagon. 4 quadrature points, the residence is O(h^6)
C
C
C
C
	
c	CTR(:) = 0
c	DO I = 1,8
c		CTR(:) = CTR(:) + NODE(:,I)
c	END DO
c	CTR(:) = CTR(:)/8.0
	Pyramid (:,1) = CTR(:)

	DO I = 1,6
	  FC(:,I) = 0.0
	  DO J = 1,4	
		IP = CNTIVITY(IKIND,I,J)
		FC(:,I) = FC(:,I) + NODE(:,IP)
	  END DO
	  FC(:,I) = FC(:,I)/4.0
	END DO



	IF(iautoflag.EQ.1) THEN
		
				
		DO I = 1,6
		  
		  CALL GetDistance(FC(1,I),FC(2,I),FC(3,I),
     &					   V0(1),V0(2),V0(3),ds(I))
		  IF(ds(I).GT.0.0)	THEN
			fun(I) = exp(-ds(I)*alfa)/(ds(I)*ds(I))
		  ELSE
			fun(I) = 0
		  END IF
		END DO
		
	END IF
	

	IF(icosflag.EQ.1) THEN
			
		DO I = 1,6
		  CALL CosValue(V0(1)-FC(1,I),V0(2)-FC(2,I),V0(3)-FC(3,I),
     &					NV(1),NV(2),NV(3),cs(I))
		  cs(I) = max(0.d0,cs(I))
		  fun(I) = fun(I)*cs(I)
		END DO
		
	END IF

	
	DO I = 1,6
		DO J = 1,4	
			IP = CNTIVITY(iKIND,I,J)
		    Pyramid(:,J+1) = NODE(:,IP)
		END DO
		CALL PyramidVolume(Pyramid,V(I))
	END DO

	value = 0
	DO I = 1,6
		value = value + fun(I)*V(I)
	END DO


	RETURN

	END 





C======================================================================C
C                                                                      C
C Subroutine GaussVolQita()									  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/09/02                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE	GaussVolQita_bkp(value,NODE,N3D,NP,V0,CTR,fun,
     &					  alfa,volume,NV,iautoflag,icosflag)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(3,N3D),ds(NP+1),NV(3)
	REAL*8  CTR(3),VC(3), V0(3),fun(NP),cs(NP)


C	This subroutine gives the GAUSS approximation of the volume-integration
C	of a Prism or Pyramid element, or any other. 
C	
C	The method should be improved
C
C
C
C


	IF(iautoflag.EQ.1) THEN

		DO I =1,NP
		  CALL GetDistance(NODE(1,I),NODE(2,I),NODE(3,I),
     &		V0(1),V0(2),V0(3),ds(I))
		END DO
		
		CALL GetDistance(VC(1),VC(2),VC(3),V0(1),V0(2),V0(3),ds(NP+1))
		     		
		DO I =1,NP+1
			IF(ds(I).GT.0.0)	THEN
				fun(I) = exp(-ds(I)*alfa)/(ds(I)*ds(I))
			ELSE
				fun(I) = 0
			END IF
		END DO

	END IF
	
	value = 0
	DO I = 1,NP+1
		value = value+fun(I)
	END DO
	value = value/(NP+1)*Volume


	RETURN

	END





C======================================================================C
C                                                                      C
C Subroutine GaussVolQita()									  	     C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                                       2005/09/02                     C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C                                                                      C
C  PSE GROUP                                                           C
C======================================================================C
	SUBROUTINE	GaussVolQita(value,NODE,N3D,V0,CTR,fun,NP,CNTIVITY,
     &				  NNPS,MFACE,MKIND,N2D,NFACE,IKIND,alfa,volume,NV,
     &				  iautoflag,icosflag)

	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8  NODE(3,N3D),ds(NP+1),NV(3),V(NFACE)
	REAL*8  CTR(3),FC(3,NFACE), V0(3),fun(NFACE),cs(NFACE)
	REAL*8	Pyramid(3,5),t1(3),t2(3),t3(3),t4(3)
	INTEGER	CNTIVITY(MKIND,MFACE,N2D),NNPS(MKIND,MFACE)


C	This subroutine gives the GAUSS approximation of the volume-integration
C	of a Prism or Pyramid element, or any other. 
C	
C	The method should be improved
C
C
C
C

	DO I = 1,NFACE
		  
	  FC(:,I) = 0.0
	  DO J = 1,	NNPS(IKIND,I)
		IP = CNTIVITY(IKIND,I,J)
		FC(:,I) = FC(:,I) + NODE(:,IP)
	  END DO
	  FC(:,I) = FC(:,I)/NNPS(IKIND,I)
	END DO

	
	IF(iautoflag.EQ.1) THEN
				
		DO I = 1,NFACE
		  
		  CALL GetDistance(FC(1,I),FC(2,I),FC(3,I),
     &					   V0(1),V0(2),V0(3),ds(I))
		  IF(ds(I).GT.0.0)	THEN
			fun(I) = exp(-ds(I)*alfa)/(ds(I)*ds(I))
		  ELSE
			fun(I) = 0
		  END IF
		END DO
		
	END IF
	

	IF(icosflag.EQ.1) THEN
			
		DO I = 1,NFACE
		  CALL CosValue(V0(1)-FC(1,I),V0(2)-FC(2,I),V0(3)-FC(3,I),
     &					NV(1),NV(2),NV(3),cs(I))
		  cs(I) = max(0.d0,cs(I))
		  fun(I) = fun(I)*cs(I)
		END DO
		
	END IF

	
	DO I = 1,NFACE
	  IF(NNPS(IKIND,I).EQ.3) THEN
		T1(:) = CTR(:)
		
		IP = CNTIVITY(IKIND,I,1)
		T2(:) = NODE(:,IP)

		IP = CNTIVITY(IKIND,I,2)
		T3(:) = NODE(:,IP)
		
		IP = CNTIVITY(IKIND,I,3)
		T4(:) = NODE(:,IP)
		
		
		CALL TetVolume(T1,T2,T3,T4,V(I))
	  ELSE
		
		Pyramid(:,1) = CTR(:)
		DO J = 1,NNPS(IKIND,I)
			IP = CNTIVITY(iKIND,I,J)
		    Pyramid(:,J+1) = NODE(:,IP)
		END DO
		CALL PyramidVolume(Pyramid,V(I))
		
	  END IF
	END DO

	value = 0
	DO I = 1,NFACE
		value = value + fun(I)*V(I)
	END DO


	RETURN


	END

