!
!	Subroutine List 

!	(subroutines for MonteCarlo Method)
!		subroutine	GasRayTracer()
!		subroutine	WallRayTracer()
!		subroutine	GasRayTracerInhom()
!		subroutine	WallRayTracerInhom()
!		subroutine	FindArrivalElem()
!		subroutine	FindArrivalWall()
!		subroutine	MarchOn()
!		subroutine	Random1()
!		subroutine	Random2() 
!		subroutine	Random3()
!		subroutine	Random4()
!		subroutine	ScatterAngleIso()
!		subroutine	ScatterAngleLAS()
!		subroutine	ScatterAngleTANI()
!		subroutine	EmitAngleWall()
!		subroutine	GetBlockExtin()
!
!	(subroutines for ZONE Method)
!		subroutine	WallZone()	
!		subroutine	GasZone()
!		subroutine	SelfZone()
!	
!	(subroutines for Double-Grid Scheme)
!	subroutine  makegroup()
!	subroutine	symmgrouprdv()
!	
!	(subroutines for Making Bucket Grid)
!	subroutine stackvolume()
!	subroutine stackwall()
!	subroutine stackpatch()
!	subroutine stackvolumegroup()
!	subroutine stackwallgroup()
!	subroutine stackwallprd()
C======================================================================C
C                                                                      C
C SUBROUTINE GasRayTracer()					       C
C                                                                      C
C                              Ver1:    2005/04/11~		       C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C			                                               C
C======================================================================C
	SUBROUTINE  GasRayTracer(IUT0,VMIN,NRAY,NP,NE,NW,MP,ME,NWT,
     2		NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD,MID,
     3		IRD,IMRD,AveWArea,AveVol,NTRay,
     4		AvePAbsorb,AvePScatter,AveGAbsorb,NGP,
     5		StackOpt,IFlagGp,ModCal,GlobeCTR,FanWei,XYZ,ELEM)

	
	USE MODULE_RADSXF, ONLY : FaceNum,NNPE,NNPS,CNTIVITY,ElemKind
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE MODULE_RADSXF, ONLY : GAbsorb,PAbsorb,PScatter,WGrayDeg
	USE MODULE_RADSXF, ONLY : WElem,WABCD,WNV,WArea,Volume,WNeb
	USE MODULE_RADSXF, ONLY : RDId,RDN1,RDN2,RDIndex,RDValue
	USE MODULE_RADSXF, ONLY : EType,VIndex,PROP,Albedo,CTR,scabeta
	USE MODULE_RADGROUP, ONLY : GpCNum,GpProp
	USE module_radwork, ONLY : NumAbsorb => WK1
	USE module_radwork, ONLY : WK2
	
c	USE MODULE_RADSXF, ONLY : GAbsorb,PAbsorb,PScatter,WGrayDeg

	IMPLICIT REAL*8(A-H,O-Z)
!	IMPLICIT NONE
!
	INTEGER,INTENT(IN)    :: IUT0,NRAY

	REAL*8,PARAMETER :: PAI=3.141592653589793,SHIGMA = 5.67e-8
	
	INTEGER ELEM(N3D,ME)
	REAL*8	XYZ(NDOF,MP),GlobeCTR(3),FanWei(2,NDOF)

	REAL*8 AveExtinCoef,AvePAbsorb,AvePScatter,AveGAbsorb,AveWArea
	CHARACTER*20 StackOpt

!StackOpt=	  1 :	get exchange-area directly
!		  2 :	get exchange-area based on reciprocity
!Private
       REAL*8	NODE(NDOF,N3D),NODEP(NDOF,N3D),fun(N3D)
       REAL*8	V0(3),VT(3),VW(3),VLast(3),DV(3)
       REAL*8	NV(3),TV1(3),TV2(3)
	INTEGER NTRay

c	real*8	tmp(3000,50)


	ALLOCATE(NumAbsorb(NRD0),WK2(NRD0))
	WK2 = 0
	
	AveExtinCoef = AveGAbsorb + AvePAbsorb + AvePScatter
	AveAlbedo = 0
	IF(AveExtinCoef.GT.0.and.ModCal.EQ.0) THEN
		AveAlbedo = AvePScatter/AveExtinCoef
	END IF
	IF(AveExtinCoef.EQ.0.d0) AveAlbedo=1.10
	iben = 0
	ireset=0
	wholedist = 0
	DO  I=1,3
		wholedist = max(wholedist,FanWei(2,I)-FanWei(1,I))
	END DO
	

	SEEDX = 0				!position
	SEEDY = 0
	SEEDZ = 0
	SEED1 = 0				!angle
	SEED2 = 0
	SEED3 = 0				!distance
	SEED4 = 0
	SEED5 = 0				!albedo
	SEED6 = 0
	NTRAY0=NTRAY

	DO 1000 IE = 1, NE
		
	IF(VIndex(IE).LT.1) CYCLE
		
!IF(MOD(VIndex(IE),100).EQ.0) WRITE(IUT0,*) 'IE = ',IE
!boundary condition element

C--------------------
C	Pre-treat
C--------------------
    

	DO 100 J =1,NNPE(ElemKind(IE))
		IP = ELEM(J,IE)
		DO 100 IDOF = 1,NDOF
		  NODE(IDOF,J) = XYZ(IDOF,IP)
100	CONTINUE


!Rearrange the conectivity
	IF(NNPE(ElemKind(IE)).EQ.4) THEN
		NODEP = NODE

	ELSE IF(NNPE(ElemKind(IE)).EQ.5) THEN
		CALL Pyramid5Arrange(NODE,NODEP,ElemKind(IE),
     &		NNPS,CNTIVITY,NFACE,N2D,N3D,NKIND,NDOF)

	ELSE IF(NNPE(ElemKind(IE)).EQ.6) THEN
		CALL Prism6Arrange(NODE,NODEP,ElemKind(IE),
     &		NNPS,CNTIVITY,NFACE,N2D,N3D,NKIND,NDOF)

	ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
		CALL Hex8Arrange(NODE,NODEP,ElemKind(IE),
     &		NNPS,CNTIVITY,NFACE,N2D,N3D,NKIND,NDOF)

	ELSE IF(NNPE(ElemKind(IE)).GT.8) THEN
	WRITE(IUT0,*) '!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!'
	WRITE(IUT0,*) 'IE=',IE,'NNPE=',NNPE(ElemKind(IE))
	WRITE(IUT0,*) 'Shape function for this kind of ELEM is not'
	WRITE(IUT0,*) 'available!'
	STOP
	END IF

		
C--------------------
C	Ray Tracing
C--------------------
	
10		NumAbsorb = 0
		
		
		RayRatio=MIN(MAX(1.0d0, AveVol/Volume(IE)),10.0d0)
		IF(IFlagGp.EQ.1) THEN
		  IGP=VIndex(IE)
	  !IF(GpCNum(IGP).GT.1) RayRatio=MAX(1.0/GpCNum(IGP),0.01)
	  IF(GpCNum(IGP).GT.1) RayRatio=MAX(PROP(IE)/GpPROP(IGP),0.1)
		END IF
		
		NRAY1=NRAY*RayRatio
		NTRay=NTRay+NRAY1

		!WRITE(IUT0,'("IE=",I5,", NRAY=",I6) ') IE,NRAY1
		IF(MOD(IE,100).EQ.1) 
     &		WRITE(IUT0,'("     IE=",I8," /[",I8,"]") ') IE-1,NE

		DO 500 IRAY = 1, NRAY1

c		if(mod(iray,5000).EQ.0.OR.IRAY.EQ.NRAY1) 
c     &				WRITE(IUT0,*) '   IRAY = ',IRAY

!emission start position V0
 123		CALL RANDOM1(SEEDX,fun(1))			
		CALL RANDOM2(SEEDY,fun(2))
		CALL RANDOM3(SEEDZ,fun(3))

		IF(NNPE(ElemKind(IE)).EQ.4) THEN
			CALL Tet4Pos(NODEP,N3D,V0,fun,NDOF)
		ELSE IF(NNPE(ElemKind(IE)).EQ.5) THEN
			CALL Pyramid5Pos(NODEP,N3D,V0,fun,NDOF)
		ELSE IF(NNPE(ElemKind(IE)).EQ.6) THEN
			CALL Prism6Pos(NODEP,N3D,V0,fun,NDOF)
		ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
			CALL Hex8Pos(NODEP,N3D,V0,fun,NDOF)
		END IF
		VLast = V0

!emission angle
		CALL RANDOM4(SEED1,SEED2,RanSeta)
		CALL RANDOM4(SEED1,SEED2,RanEta)
		Seta = RanSeta*2.0*PAI
		Eta = ACOS(1-2*RanEta)
	
			
!c------------this is a loop for ray emission and tracing

!emission possible distance
 150		CALL RANDOM4(SEED3,SEED4,RanDist)
		IF(AveExtinCoef.GT.0) THEN
			Dist0 = -1.0*LOG(1.0-RanDist)/AveExtinCoef
		ELSE
			Dist0 = 10.0*wholedist
		END IF

200		IF(Seta.LT.0) Seta = Seta+2*PAI
		
		DV(1) = sin(Eta)*cos(Seta)*Dist0
		DV(2) = sin(Eta)*sin(Seta)*Dist0
		DV(3) = cos(Eta)*Dist0
		VT(1) = V0(1) + DV(1)
		VT(2) = V0(2) + DV(2)
		VT(3) = V0(3) + DV(3)
!looking for an intersecting wall surface if any
		CALL FindArrivalWall(DistW,Dist0,V0,VT,VW,DV,Seta,Eta,
     &		  IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NWT,iray,ierr,iben,
     &		  XYZ,FanWei)
	
		IF(ierr.GT.0) THEN
		 ireset=ireset+1
		 GOTO 123
	        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! found no surface (idwall=0), or it is too far, (idwall<0)
!		
         IF(IDWall.LE.0) THEN
                CALL FindArrivalElem(iray,MatchID,VT,N2D,N3D,NDOF,
     &		NB,NPBG,NP,NE,MP,ME,NFACE,NKIND,XYZ,ELEM,FanWei)
	
			
		IF(MatchID.LE.0) THEN
                  ireset=ireset+1
		  GOTO 123
		END IF
!particle scattering for participating media
!for non-scattering gas, please set Albedo = 0
                CALL RANDOM4(SEED5,SEED6,RanAlbedo)
		IF(RanAlbedo.GE.AveAlbedo) THEN
		  NumAbsorb(MatchID) = NumAbsorb(MatchID) + 1
!scattered
		  ELSE
			V0 = VT			!New emission position
			CALL RANDOM4(SEED1,SEED2,RanEta1)
			CALL RANDOM4(SEED1,SEED2,RanSeta1)
			
			!scattering	angle
			
			!Eta1 = ACOS(1-2*RanEta1)
			!Seta1 = RanSeta1*2.0*PAI
			!reset travel	
			CALL ScatterAngleLAS(Seta,Eta,RanSeta1,RanEta1,
     &			EType(MatchID),scabeta(MatchID))
				
			GOTO 150		!distance reset

	          END IF
			
		  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intersect with a surface patch
!
		  ELSE
			
			CALL RANDOM4(SEED5,SEED6,RanWGray)
			
			JRD = VIndex(IDWALL+NE)
			!absorbed
		IF(JRD.GT.0.AND.RanWGray.LE.1.0-Albedo(IDWALL+NE)
     &			.and.DistW.GT.0.0) THEN
			NumAbsorb(IDWall+NE) = NumAbsorb(IDWall+NE) + 1
			
!periodic boundary wall
			ELSE IF(WGrayDeg(IDWALL).LT.-0.9) THEN	
!the value is the ID of wall-surf	
				
				IPair = ABS(WGrayDeg(IDWALL))

				Dist0 = Dist0 - DistW
				
		V0(1) = VW(1) + (CTR(1,NE+IPair)-CTR(1,NE+IDWALL))
		V0(2) = VW(2) + (CTR(2,NE+IPair)-CTR(2,NE+IDWALL))
		V0(3) = VW(3) + (CTR(3,NE+IPair)-CTR(3,NE+IDWALL))
		Vlast = V0
		GOTO 200

!reflected
		ELSE
!reflected by wall
		IF(DistW.GT.0.0) THEN
!performed one true step of reflection
		  Dist0 = Dist0 - DistW
		  Vlast = V0
		  V0 = VW
                ELSE
!the present step is idle, resume it
!this case probably arise for the ray in corner.
                  V0 = Vlast
		END IF
!Isotropic reflecting face
		CALL RANDOM4(SEED1,SEED2,RanSeta)
		CALL RANDOM4(SEED1,SEED2,RanEta)
		NV(1) = WNV(1,IDWall)
		NV(2) = WNV(2,IDWall)
		NV(3) = WNV(3,IDWall)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WNV is the outward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CALL BuildZuoBiao(NV,TV1,TV2,3)
		CALL ReflectAngleMC(Seta,Eta,RanSeta,RanEta,NV,
     &					EType(IDWALL+NE))
		GOTO 200
		END IF
                END IF

500		CONTINUE




			
C-----------------------------
C	Calcuate R-E-A-D values
C-----------------------------

	
	IF(StackOpt(1:4).EQ.'BAND') GOTO  699


!----------------------------------Stack Method 1------------------------------------------
!---
!---	Stack it into a 1-D array, anologous to a Lower-Trianglar-Array (NRD0,NRD0)/2
!---
!------------------------------------------------------------------------------------

		
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************!
         IF(IFlagGp.NE.1) THEN	!not grouped

           IRD = VIndex(IE)
	   ratio =  1.0*NumAbsorb(IE)/NRAY1
	   RDValue(RDIndex(IRD)+IRD) = PROP(IE)*ratio
	   selfratio = ratio
	   sumofratio = ratio
           DO 600 J = 1,IE-1
	    JRD = VIndex(J)
		 IF(JRD.GT.0) THEN
		ratio1 = 1.0*NumAbsorb(J)/NRAY1	!present value
		sumofratio = sumofratio + ratio1
		ratio1 = ratio1*PROP(IE)
		ratio2 = RDValue(RDIndex(IRD)+JRD)	
!call last value
		CALL  SymmTreatRD(ratio1,ratio2,PROP(IE),PROP(J))
		RDValue(RDIndex(IRD)+JRD) = ratio1
		 ELSE IF(NumAbsorb(J).GT.0) THEN
			WRITE(IUT0,*) 'ID = ',J
			WRITE(IUT0,*) 'An element received ray, '
			WRITE(IUT0,*) 'but does not appear in equation'
			STOP
		 END IF	
		 
600		CONTINUE

				
		DO 620 J = IE+1,NRD0
	 
		 JRD = VIndex(J)
		 IF(JRD.GT.0) THEN
			 ratio = 1.0*NumAbsorb(J)/NRAY1
!present value deposit to (JRD,IRD)
			 RDValue(RDIndex(JRD)+IRD) = ratio*PROP(I)
			 sumofratio = sumofratio + ratio
		 ELSE IF(NumAbsorb(J).GT.0) THEN
			WRITE(IUT0,*) 'ID = ',J,'JRD=',JRD
			WRITE(IUT0,*) 'An element received ray, '
		WRITE(IUT0,*) 'but does not appear in equation'
			STOP
		 END IF		 
620		CONTINUE
		

	
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !**********************************************************!
	  ELSE		!Grouped
		

		IGP = VIndex(IE)
		IRD = VIndex(IE)
		sumofratio = 0
		selfratio =  1.0*NumAbsorb(IE)/NRAY1
c		ID = IGP*NGP+IGP
c		RDValue(ID) =RDValue(ID)+PROP(IRD)*ratio

		DO 640 J = 1,NE+NW
		 
		 IF(NumAbsorb(J).EQ.0) CYCLE

		 JGP = VIndex(J)
		 ratio1 = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio1

		 IF(JGP.GT.0) THEN
			 ID = (IGP-1)*NGP+JGP
			 RDValue(ID)=RDValue(ID)+PROP(IE)*ratio1
		 ELSE
			WRITE(IUT0,*) 'ID = ',J, 'JGP=',JGP
			WRITE(IUT0,*) 'An element received ray, '
			WRITE(IUT0,*) 'but does not appear in equation'
			STOP
		 END IF
		 
640		CONTINUE


	  END IF

	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c633	  WRITE(IUT0,'(a,E16.8)') '     Sum of values=', sumofratio
c	  WRITE(IUT0,'(a,E16.8)') '     Self absorption=',selfratio

	  GOTO 1000
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!----------------------------------Method-2------------------------------------------
!---
!---	Stack view-factors in a 1-D array, with band NNN and Element Index in RDId(:)
!---
!------------------------------------------------------------------------------------
	
699		IRD = VIndex(IE)
		RDIndex(IRD) = IMRD
		IMRD = IMRD + 1
		ratio = 1.0*NumAbsorb(IE)/NRAY1
		RDId(IMRD) = IRD
		
		RDValue(IMRD) = ratio*PROP(IE)
		sumofratio = ratio

		IF(IMRD.GT.MRD-MRD/NRD0) THEN
			WRITE(IUT0,*) 'Error in GasRayTracer():'
	WRITE(IUT0,*) 'Stack maybe overflow, please enlarge MRD= '
     &			,MRD
			STOP
		END IF
				
		
		DO 700 J = 1,IE-1
		 ratio = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio

		 JRD = VIndex(J)
		 IF(ratio.GT.0) THEN
			IF(JRD.GT.0) THEN
				IMRD = IMRD + 1
				RDId(IMRD) = JRD
				RDValue(IMRD) = ratio*PROP(IE)
			 ELSE
				WRITE(IUT0,*) 'ID = ',J
		WRITE(IUT0,*) 'An element received ray, '
		WRITE(IUT0,*) 'but does not appear in equation'
				STOP
			 END IF	
		 END IF

700		CONTINUE
		
		DO 710 J = IE+1,NRD0
		 ratio = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio

		 JRD = VIndex(J)
		 IF(ratio.GT.0) THEN
			IF(JRD.GT.0) THEN
				IMRD = IMRD + 1
				RDId(IMRD) = JRD
				RDValue(IMRD) = ratio*PROP(IE)
			ELSE
		WRITE(IUT0,*) 'ID = ',J
		WRITE(IUT0,*) 'An element received ray, '
		WRITE(IUT0,*) 'but does not appear in equation'
		STOP
			END IF
	     END IF
710		CONTINUE


999		CONTINUE
c		WRITE(IUT0,*) '   Total Record Num=', IMRD-RDIndex(IRD)
c		WRITE(IUT0,*) '   Sum of values=',sumofratio	
		
1000  CONTINUE

	WRITE(IUT0,'(a,I10)') '    Total tracing rays = ',NTRay-NTRAY0
	WRITE(IUT0,'(a,I5)') '     Expanded-search times =', iben
	WRITE(IUT0,'(a,I5)') '     Ray reset times       =', ireset


	DEALLOCATE(NumAbsorb,WK2)

	RETURN

!-------> end of GasRayTracer()

	END 
C======================================================================C
	SUBROUTINE  GasRayTracerInhom(IUT0,VMIN,NRAY,NP,NE,NW,MP,ME,NWT,
     2		NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD,MID,IRD,
     3		IMRD,dseg,AveWArea,AveVol,NTRay,
     4		AvePAbsorb,AvePScatter,AveGAbsorb,NGP,
     5		StackOpt,IFlagGp,ModCal,GlobeCTR,FanWei,
     6		XYZ,ELEM)
C======================================================================C
	USE MODULE_RADSXF, ONLY : FaceNum,NNPE,NNPS,CNTIVITY,ElemKind
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE MODULE_RADSXF, ONLY : GAbsorb,PAbsorb,PScatter,WGrayDeg
	USE MODULE_RADSXF, ONLY : WElem,WABCD,WNV,WArea,Volume,WNEB
	USE MODULE_RADSXF, ONLY : RDId,RDN1,RDN2,RDIndex,RDValue
	USE MODULE_RADSXF, ONLY : EType,VIndex,PROP,Albedo,CTR,scabeta
	USE MODULE_RADGROUP, ONLY : GpCNum,GpProp
	USE module_radwork, ONLY : NumAbsorb => WK1
	USE module_radwork, ONLY : VWK3D,WK2

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)

	INTEGER ELEM(N3D,ME)
	REAL*8	XYZ(NDOF,MP),GlobeCTR(3),FanWei(2,NDOF)

	REAL*8 AveExtinCoef,AvePAbsorb,AvePScatter,AveGAbsorb,AveWArea
	CHARACTER*20 StackOpt

!StackOpt=	  1 :	get exchange-area directly
!		  2 :	get exchange-area based on reciprocity


	!Private
	REAL*8	NODE(NDOF,N3D),NODEP(NDOF,N3D),fun(N3D)
      REAL*8	V0(3),VT(3),VW(3),VLast(3),DV(3)
	REAL*8	NV(3),TV1(3),TV2(3),ddd(3)
	INTEGER NTRay

c	real*8	tmp(3000,50)


	!ALLOCATE(NumAbsorb(NRD0),WK2(NRD0),VWK3D(NB,NB,NB))
	ALLOCATE(NumAbsorb(NRD0),WK2(NRD0))
	WK2 = 0	

	!CALL GetBlockExtin(IUT0,NB,NPBG)

	

	AveExtinCoef = AveGAbsorb + AvePAbsorb + AvePScatter
	
	iben = 0
	ireset=0
	wholedist = 0
	DO  I=1,3
		wholedist = max(wholedist,FanWei(2,I)-FanWei(1,I))
	END DO

	DO I = 1,3
		ddd(I) = (FanWei(2,I)-FanWei(1,I))/NB
	END DO
	
	

	SEEDX = 0				!position
	SEEDY = 0
	SEEDZ = 0
	SEED1 = 0				!angle
	SEED2 = 0
	SEED3 = 0				!distance
	SEED4 = 0
	SEED5 = 0				!albedo
	SEED6 = 0
	NTRAY0=NTRAY
		
	DO 1000 IE = 1, NE

	
		
		IF(VIndex(IE).LT.1) CYCLE
		
		
C--------------------
C	Pre-treat
C--------------------
    

		DO 100 J =1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			DO 100 IDOF = 1,NDOF
			  NODE(IDOF,J) = XYZ(IDOF,IP)
100		CONTINUE


		!Rearrange the conectivity
		IF(NNPE(ElemKind(IE)).EQ.4) THEN
			NODEP = NODE

		ELSE IF(NNPE(ElemKind(IE)).EQ.5) THEN
			CALL Pyramid5Arrange(NODE,NODEP,ElemKind(IE),
     &			NNPS,CNTIVITY,NFACE,N2D,N3D,NKIND,NDOF)

		ELSE IF(NNPE(ElemKind(IE)).EQ.6) THEN
			CALL Prism6Arrange(NODE,NODEP,ElemKind(IE),
     &			NNPS,CNTIVITY,NFACE,N2D,N3D,NKIND,NDOF)

		ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
			CALL Hex8Arrange(NODE,NODEP,ElemKind(IE),
     &			NNPS,CNTIVITY,NFACE,N2D,N3D,NKIND,NDOF)

		ELSE IF(NNPE(ElemKind(IE)).GT.8) THEN
	WRITE(IUT0,*) '!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!'
	WRITE(IUT0,*) 'IE=',IE,'NNPE=',NNPE(ElemKind(IE))
	WRITE(IUT0,*) 'Shape function for this kind of ELEM is not'
	WRITE(IUT0,*) 'available!'
			STOP
		END IF


			
C--------------------
C	Ray Tracing
C--------------------
	
10		NumAbsorb = 0
		
		RayRatio=MIN(MAX(1.0d0, AveVol/Volume(IE)),10.0d0)
		IF(IFlagGp.EQ.1) THEN
		  IGP=VIndex(IE)
	  IF(GpCNum(IGP).GT.1) RayRatio=MAX(PROP(IE)/GpPROP(IGP),0.1)
		END IF

		NRAY1=NRAY*RayRatio
		NTRay=NTRay+NRAY1

		!WRITE(IUT0,'("IE=",I5,", NRAY=",I6) ') IE,NRAY1
		IF(MOD(IE,100).EQ.1) 
     &		WRITE(IUT0,'("     IE=",I8," /[",I8,"]") ') IE-1,NE

		DO 500 IRAY = 1, NRAY1


c		if(mod(iray,5000).EQ.0.OR.IRAY.EQ.NRAY1) 
c     &				WRITE(IUT0,*) '   IRAY = ',IRAY

			!emission start position V0
123			CALL RANDOM1(SEEDX,fun(1))	
			CALL RANDOM2(SEEDY,fun(2))
			CALL RANDOM3(SEEDZ,fun(3))

			IF(NNPE(ElemKind(IE)).EQ.4) THEN
				CALL Tet4Pos(NODEP,N3D,V0,fun,NDOF)

			ELSE IF(NNPE(ElemKind(IE)).EQ.5) THEN
				CALL Pyramid5Pos(NODEP,N3D,V0,fun,NDOF)

			ELSE IF(NNPE(ElemKind(IE)).EQ.6) THEN
				CALL Prism6Pos(NODEP,N3D,V0,fun,NDOF)

			ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
				CALL Hex8Pos(NODEP,N3D,V0,fun,NDOF)

			END IF

			VLast = V0

			!emission angle
			CALL RANDOM4(SEED1,SEED2,RanSeta)
			CALL RANDOM4(SEED1,SEED2,RanEta)
			Seta = RanSeta*2.0*PAI
			Eta = ACOS(1-2*RanEta)
!c------------this is a loop for ray emission and tracing
			
			!emission possible distance
150			CALL RANDOM4(SEED3,SEED4,RanDist)
		ExtinCoef = GAbsorb(IE)+PAbsorb(IE)+PScatter(IE)
			Dist0 = -1.0*LOG(1.0-RanDist)/ExtinCoef
			
			IE1 = IE
			Vlast = V0

!looking for a final volume or intersecting wall surface if any
			MarkDaoda = 0
			IDWALL = 0
			
			DO 333 WHILE(MarkDaoda.EQ.0)
				
				IF(Seta.LT.0) Seta = Seta+2*PAI
				
!CALL  GetSegment(DV,ddd,dseg,2.0)
								
		ExtinCoef1 = GAbsorb(IE1)+PAbsorb(IE1)+PScatter(IE1)
	        IF(ExtinCoef1.GT.0.d0) THEN
		ExtinCoef0 = ExtinCoef
		ExtinCoef = ExtinCoef1
		Dist0 = Dist0*ExtinCoef0/ExtinCoef
		DistGo = min(Dist0,dseg)
!possible travel distance
		ELSE
		DistGo = dseg
		END IF
!looking for an intersecting wall surface if any
		IF(IDWALL.EQ.0) THEN
		  DV(1) = sin(Eta)*cos(Seta)*Dist0
!*wholedist*1.e5	
		  DV(2) = sin(Eta)*sin(Seta)*Dist0
		  DV(3) = cos(Eta)*Dist0

	          VT(1) = V0(1) + DV(1)
		  VT(2) = V0(2) + DV(2)
                  VT(3) = V0(3) + DV(3)
	  CALL FindArrivalWall(DistW,Dist0,V0,VT,VW,DV,Seta,Eta,
     &		IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NWT,iray,ierr,iben,
     &			    XYZ,FanWei)
				

	  IF(IDWALL.EQ.0) THEN
	!Error abort
c		WRITE(IUT0,*) 'ERROR in GasRayTracerInhom():'
c		WRITE(IUT0,*) 'iray =',iray
c		WRITE(IUT0,*) v0,dv
c		WRITE(IUT0,*) 'The ray intersect no wall'
c    		WRITE(IUT0,*) 'please check the mesh or code !!!'
c		STOP
c		ireset=ireset+1
c		goto 123
		DistW=Dist0+DistGo+wholedist
	  END IF

				
	ELSE
	DistW = SQRT((VW(1)-V0(1))*(VW(1)-V0(1))+ 
     &	(VW(2)-V0(2))*(VW(2)-V0(2)) +
     &	(VW(3)-V0(3))*(VW(3)-V0(3))	)
	END IF

	DV(1) = sin(Eta)*cos(Seta)*DistGo			!Dist0
	DV(2) = sin(Eta)*sin(Seta)*DistGo			!Dist0
	DV(3) = cos(Eta)*DistGo					!Dist0
	VT(1) = V0(1) + DV(1)
	VT(2) = V0(2) + DV(2)
	VT(3) = V0(3) + DV(3)

							
!march on until encountering an element with different property
!	(i.e., until the optical-extinction coef. is changed)
! ExtinCoef keeps unchanged, but Dist0,DistW and V0,VT will be
! changed
        imarchon=0
!arrive at a volume element
	IF(DistW.GT.DistGo.AND.imarchon.LE.2) THEN
	  CALL FindArrivalElem(iray,MatchID,VT,N2D,N3D,NDOF,
     &    NB,NPBG,NP,NE,MP,ME,NFACE,NKIND,XYZ,ELEM,FanWei)

        IF(MatchID.LE.0) THEN
!Error abort
c		WRITE(IUT0,*) 'ERROR in GasRayTracerInhom():'
c		WRITE(IUT0,*) 'iray =',iray
c		WRITE(IUT0,*) v0,dv
c	WRITE(IUT0,*) 'The ray goes to unknown place'
c		WRITE(IUT0,*) 'please check the mesh or code !!!'
c	STOP
	ireset=ireset+1
	goto 123
	END IF
	ELSE
	  MatchID=WNEB(ABS(IDWALL))
	END IF
				

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	intersect a surface patch
!
	IF(imarchon.EQ.3.OR.
     &	(IDWALL.NE.0.AND.DistW.LE.DistGo))THEN
	IDWALL = ABS(IDWALL)
!intersect with the wall surface
	CALL RANDOM4(SEED5,SEED6,RanWGray)
	JRD = VIndex(IDWALL+NE)	
!absorbed
	IF(JRD.GT.0.AND.
     &	RanWGray.LE.1.0-Albedo(IDWALL+NE)) THEN
	NumAbsorb(IDWall+NE) = NumAbsorb(IDWall+NE)+1
	MarkDaoda = 1
!periodic boundary wall
	ELSE IF(WGrayDeg(IDWALL).LT.-0.9) THEN	
!the value is the ID of wall-surf	
	IPair = ABS(WGrayDeg(IDWALL))
	Dist0 = Dist0 - DistW
        V0(1)=VW(1)+(CTR(1,NE+IPair)-CTR(1,NE+IDWALL))
	V0(2)=VW(2)+(CTR(2,NE+IPair)-CTR(2,NE+IDWALL))
	V0(3)=VW(3)+(CTR(3,NE+IPair)-CTR(3,NE+IDWALL))
	Vlast=V0
!reflected
	ELSE
!reflected by wall
cIF(DistW.GT.0.0) THEN
					
!performed one true step of reflection
	IE2 = WNEB(IDWALL)
	ExtinCoef2 = GAbsorb(IE2)+PAbsorb(IE2)
     &	 + PScatter(IE2)
	IF(ExtinCoef2.GT.0) THEN
	  ratio=2.*ExtinCoef/(ExtinCoef+ExtinCoef2)
	  Dist0 = Dist0 - DistW*ratio
	END IF
	Vlast = V0
	V0 = VW
	IE1 = IE2
c	ELSE
c	!the present step is idle, resume it
c	!this case probably arise for the ray in corner.
c	V0 = Vlast
cEND IF
!Isotropic reflecting face
	CALL RANDOM4(SEED1,SEED2,RanSeta)
	CALL RANDOM4(SEED1,SEED2,RanEta)
c				Seta = RanSeta*2.0*PAI
c				Eta = ACOS(SQRT(1.0-RanEta))

	NV(1) = WNV(1,IDWall)
	NV(2) = WNV(2,IDWall)
	NV(3) = WNV(3,IDWall)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WNV is the outward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c				CALL BuildZuoBiao(NV,TV1,TV2,3)
c				
c	CALL ZuobiaoRot(Seta,Eta,NV(1),NV(2),NV(3),
c     &		TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
						
	CALL ReflectAngleMC(Seta,Eta,RanSeta,RanEta,
     &	NV,EType(IDWALL+NE))
	IDWALL = 0	!reset
	END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	arrive at an element
!
!NO-Participating element, pass it WITHOUT distance-reduction
	ELSE IF(VIndex(MatchID).EQ.0) THEN	
!not a participating element
	V0 = VT
	IE2 = MatchID
	Vlast = V0
	IE1 = MatchID
				
!participating element, pass it WITH distance-reduction
        ELSE IF(Dist0.GT.1.01*DistGo) THEN
!continue the process
	V0 = VT
	IE2 = MatchID
	ExtinCoef2 =GAbsorb(IE2)+PAbsorb(IE2)+PScatter(IE2)
	ratio = 2.0*ExtinCoef/(ExtinCoef+ExtinCoef2)
	Dist0 = Dist0 - DistGo*ratio
	Vlast = V0
	IE1 = MatchID
!final position, travel terminate at this element
	ELSE
!arrive at a volume
!particle scattering for participating media
!for non-scattering gas, please set Albedo = 0
	CALL RANDOM4(SEED3,SEED4,RanAlbedo)
        JRD =VIndex(MatchID)
					
!absorbed
	IF(RanAlbedo.GE.Albedo(MatchID)) THEN
	NumAbsorb(MatchID) = NumAbsorb(MatchID) + 1
	MarkDaoda = 1

!scattered
	ELSE
	V0 = VT			!New emission position
	CALL RANDOM4(SEED1,SEED2,RanEta1)
	CALL RANDOM4(SEED1,SEED2,RanSeta1)
!scattering angle 
!Eta1 = ACOS(1-2*RanEta1)
!Seta1 = RanSeta1*2.0*PAI
!reset travel
	CALL ScatterAngleLAS(Seta,Eta,RanSeta1,
     &	RanEta1,EType(MatchID),scabeta(MatchID))

	CALL RANDOM4(SEED5,SEED6,RanDist)
	IE1 = MatchID
	ExtinCoef = GAbsorb(IE1)+PAbsorb(IE1)
     &		+ PScatter(IE1)
	Dist0 = -1.0*LOG(1.0-RanDist)/ExtinCoef		
!reset travel
	Vlast = V0
	IDWALL = 0	!reset
	END IF
        END IF
 333	CONTINUE


 500	CONTINUE
C-----------------------------
C	Calcuate R-E-A-D values
C-----------------------------

	
        IF(StackOpt(1:4).EQ.'BAND') GOTO  699


!----------------------------------Stack Method 1------------------------------------------
!---
!---	Stack it into a 1-D array, anologous to a Lower-Trianglar-Array (NRD0,NRD0)/2
!---
!------------------------------------------------------------------------------------
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************!
	IF(IFlagGp.NE.1) THEN		!not grouped

	IRD = VIndex(IE)
	ratio =  1.0*NumAbsorb(IE)/NRAY1
	selfratio=ratio
	RDValue(RDIndex(IRD)+IRD) = PROP(IE)*ratio

	sumofratio = ratio

	DO 600 J = 1,IE-1
		 
	 JRD = VIndex(J)
	 IF(JRD.GT.0) THEN
	   ratio1 = 1.0*NumAbsorb(J)/NRAY1		!present value
           sumofratio = sumofratio + ratio1
	   ratio1 = ratio1*PROP(IE)
	   ratio2 = RDValue(RDIndex(IRD)+JRD)	!call last value
           CALL  SymmTreatRD(ratio1,ratio2,PROP(IE),PROP(J))
		RDValue(RDIndex(IRD)+JRD) = ratio1
           ELSE IF(NumAbsorb(J).GT.0) THEN
		WRITE(IUT0,*) 'ID = ',J
		WRITE(IUT0,*) 'An element received ray, '
		WRITE(IUT0,*) 'but does not appear in equation'
		STOP
           END IF	
		 
600	   CONTINUE
	   DO 620 J = IE+1,NRD0
           JRD = VIndex(J)
	   IF(JRD.GT.0) THEN
	     ratio = 1.0*NumAbsorb(J)/NRAY1
!present value deposit to (JRD,IRD)
             RDValue(RDIndex(JRD)+IRD) = ratio*PROP(IE)
             sumofratio = sumofratio + ratio
	     ELSE IF(NumAbsorb(J).GT.0) THEN
	     WRITE(IUT0,*) 'ID = ',J
             WRITE(IUT0,*) 'An element received ray, '
	     WRITE(IUT0,*) 'but does not appear in equation'
             STOP
	   END IF		 
620	   CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************!
           ELSE	!Grouped
	    IGP = VIndex(IE)
	    IRD = VIndex(IE)
		sumofratio = 0
		selfratio =  1.0*NumAbsorb(IE)/NRAY1
c		ID = IGP*NGP+IGP
c		RDValue(ID) =RDValue(ID)+PROP(IRD)*ratio

		DO 640 J = 1,NE+NW
		 
		 IF(NumAbsorb(J).EQ.0) CYCLE

		 JGP = VIndex(J)
		 ratio1 = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio1

		 IF(JGP.GT.0) THEN
			 ID = (IGP-1)*NGP+JGP
			 RDValue(ID)=RDValue(ID)+PROP(IE)*ratio1
		 ELSE
			WRITE(IUT0,*) 'ID = ',J, 'JGP=',JGP
			WRITE(IUT0,*) 'An element received ray, '
			WRITE(IUT0,*) 'but does not appear in equation'
			STOP
		 END IF
		 
640		CONTINUE


	  END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c633	  WRITE(IUT0,'(a,E16.8)') '     Sum of values=', sumofratio
c	  WRITE(IUT0,'(a,E16.8)') '     Self absorption=',selfratio

          GOTO 1000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------Method-2------------------------------------------
!---
!---	Stack view-factors in a 1-D array, with band NNN and Element Index in RDId(:)
!---
!------------------------------------------------------------------------------------
	
 699	IRD = VIndex(IE)
	RDIndex(IRD) = IMRD
	IMRD = IMRD + 1
	ratio = 1.0*NumAbsorb(IE)/NRAY1
	RDId(IMRD) = IRD
		
	RDValue(IMRD) = ratio*PROP(IE)
	sumofratio = ratio

	IF(IMRD.GT.MRD-MRD/NRD0) THEN
	WRITE(IUT0,*) 'Error in GasRayTracer():'
	WRITE(IUT0,*) 'Stack maybe overflow, please enlarge MRD= '
     &	,MRD
	STOP
	END IF
				
		
		DO 700 J = 1,IE-1
		 ratio = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio

		 JRD = VIndex(J)
		 IF(ratio.GT.0) THEN
			IF(JRD.GT.0) THEN
				IMRD = IMRD + 1
				RDId(IMRD) = JRD
				RDValue(IMRD) = ratio*PROP(IE)
			 ELSE
			WRITE(IUT0,*) 'ID = ',J
			WRITE(IUT0,*) 'An element received ray, '
			WRITE(IUT0,*) 'but does not appear in equation'
			STOP
			 END IF	
		 END IF

700		CONTINUE
		
		DO 710 J = IE+1,NRD0
		 ratio = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio

		 JRD = VIndex(J)
		 IF(ratio.GT.0) THEN
			IF(JRD.GT.0) THEN
				IMRD = IMRD + 1
				RDId(IMRD) = JRD
				RDValue(IMRD) = ratio*PROP(IE)
			ELSE
			WRITE(IUT0,*) 'ID = ',J
			WRITE(IUT0,*) 'An element received ray, '
			WRITE(IUT0,*) 'but does not appear in equation'
			STOP
			END IF
	     END IF
710		CONTINUE


999		CONTINUE
		
c		WRITE(IUT0,*) '   Total Record Num=', IMRD-RDIndex(IRD)
c		WRITE(IUT0,*) '   Sum of values=',sumofratio	
		
1000  CONTINUE

	WRITE(IUT0,'(a,I10)') '    Total tracing rays = ',NTRay-NTRAY0
	WRITE(IUT0,'(a,I5)') '     Expanded-search times =', iben
	WRITE(IUT0,'(a,I5)') '     Ray reset times       =', ireset


	!DEALLOCATE(NumAbsorb,WK2,VWK3D)
	DEALLOCATE(NumAbsorb,WK2)

	RETURN

!-------> end of GasRayTracerInhom()

	END 
C======================================================================C
	SUBROUTINE  WallRayTracer(IUT0,VMIN,NRAY,NP,NE,NW,MP,ME,NWT,
     2		NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD,
     3		MID,IRD,IMRD,AveWArea,AveVol,NTRay,
     4		AvePAbsorb,AvePScatter,AveGAbsorb,NGP,
     5		StackOpt,IFlagGp,ModCal,GlobeCTR,FanWei,XYZ,ELEM)
C======================================================================C
	USE MODULE_RADSXF, ONLY : FaceNum,NNPE,NNPS,CNTIVITY,ElemKind
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE MODULE_RADSXF, ONLY : GAbsorb,PAbsorb,PScatter,WGrayDeg
	USE MODULE_RADSXF, ONLY : WElem,WABCD,WNV,WArea,Volume,WNEB
	USE MODULE_RADSXF, ONLY : RDId,RDN1,RDN2,RDIndex,RDValue
	USE MODULE_RADSXF, ONLY : EType,VIndex,PROP,Albedo,CTR,scabeta
	USE MODULE_RADGROUP, ONLY : GpCNum,GpProp
	USE module_radwork, ONLY : NumAbsorb => WK1
	USE module_radwork, ONLY : WK2
c	USE MODULE_RADSXF, ONLY : GAbsorb,PAbsorb,PScatter,WGrayDeg

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)
	
	INTEGER ELEM(N3D,ME)
	REAL*8	XYZ(NDOF,MP),GlobeCTR(3),FanWei(2,NDOF)
	
	REAL*8 AveExtinCoef,AvePAbsorb,AvePScatter,AveGAbsorb,AveWArea
	CHARACTER*20 StackOpt
					
!NW:	number of wall surfaces
!NWT:	number of total surface patches, (NW+1:NWT) are the 
!	surfaces that are not specified as wall

!Private
	REAL*8	NODE(NDOF,N3D),NODEP(NDOF,N3D),fun(N3D)
        REAL*8	V0(3),VT(3),VW(3),DV(3),VLast(3)
	REAL*8	NV(3),TV1(3),TV2(3)
	INTEGER NTRay

	ALLOCATE(NumAbsorb(NRD0),WK2(NRD0))
	WK2 = 0

	AveExtinCoef = AveGAbsorb + AvePAbsorb + AvePScatter
	AveAlbedo = 0.d0
	IF(AveExtinCoef.GT.0.and.ModCal.EQ.0) THEN
		AveAlbedo = AvePScatter/AveExtinCoef
	END IF
	IF(AveExtinCoef.EQ.0.d0) AveAlbedo=1.10

	iben=0
	ireset=0
	wholedist = 0
	DO  I=1,3
	wholedist = max(wholedist,FanWei(2,I)-FanWei(1,I))
	END DO
!
	SEEDX = 0		!position
	SEEDY = 0
	SEEDZ = 0
	SEED1 = 0		!seta,eta
	SEED2 = 0
	SEED3 = 0		!distance
	SEED4 = 0		
	SEED5 = 0		!scattering albedo
	SEED6 = 0
	NTRAY0=NTRAY	
    
	DO 1000 IW = 1, NW
	IF(VIndex(IW+NE).LT.1) CYCLE
!WRITE(IUT0,*) 'IW = ',IW
!IF(MOD(VIndex(IW+NE),100).EQ.0) WRITE(IUT0,*) 'IW = ',IW
!boundary condition element
C-----------------------
C	Pre-treat
C-----------------------
	DO 100 J =1,WElem(1,IW)
	IP = WElem(J+1,IW)
	DO 100 IDOF = 1,NDOF
        NODE(IDOF,J) = XYZ(IDOF,IP)
100	CONTINUE
C--------------------
C	Ray Tracing of wall
C--------------------
		NumAbsorb = 0
		
		RayRatio=MIN(MAX(1.0d0, AveWArea/WArea(IW)),10.0d0)
		IF(IFlagGp.EQ.1) THEN
		  IGP=VIndex(NE+IW)
		  !IF(GpCNum(IGP).GT.1) RayRatio=1.0/GpCNum(IGP)
		  IF(GpCNum(IGP).GT.1) 
     &		  RayRatio=MAX(PROP(NE+IW)/GpPROP(IGP),0.1)
		END IF
				
		NRAY1=NRAY*RayRatio
		NTRay=NTRay+NRAY1

		!WRITE(IUT0,'("IW=",I5,", NRAY=",I6) ') IW,NRAY1
		
		IF(MOD(IW,100).EQ.1) 
     &		  WRITE(IUT0,'("     IW=",I8," /[",I8,"]") ') IW-1,NW
     
	
		DO 500 IRAY = 1, NRAY1

c			if(mod(iray,5000).EQ.0.OR.IRAY.EQ.NRAY) 
c     &				WRITE(IUT0,*) '   IRAY = ',IRAY
		
			!emission start position V0
123			CALL RANDOM1(SEEDX,fun(1))		
			CALL RANDOM2(SEEDY,fun(2))
			CALL RANDOM3(SEEDZ,fun(3))

			IF(WElem(1,IW).EQ.3) THEN
				CALL Tri3Pos(NODE,N3D,V0,fun,NDOF)

			ELSE IF(WElem(1,IW).EQ.4) THEN
				CALL Qua4Pos(NODE,N3D,V0,fun,NDOF)

			ELSE IF(WElem(1,IW).GT.4) THEN
		WRITE(IUT0,*) '!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!'
		WRITE(IUT0,*) 'IW=',IW,'NNPS=',WElem(1,IW)
		WRITE(IUT0,*) 'Shape function for this kind of face is'
		WRITE(IUT0,*) 'not available!'
				STOP
			END IF

			VLast = V0

			!emission angle
			CALL RANDOM4(SEED1,SEED2,RanSeta)
			CALL RANDOM4(SEED1,SEED2,RanEta)
c			Seta = RanSeta*2.0*PAI
c			Eta = ACOS(SQRT(1-RanEta))
	
			NV(1) = WNV(1,IW)
			NV(2) = WNV(2,IW)
			NV(3) = WNV(3,IW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WNV is the outward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
c		CALL BuildZuoBiao(NV,TV1,TV2,3)
c		CALL ZuobiaoRot(Seta,Eta,NV(1),NV(2),NV(3),
c     &		TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
			

		CALL EmitAngleWall(Seta,Eta,RanSeta,RanEta,NV,
     &		EType(IW+NE))
			
	
!c------------this is a loop for ray emission and tracing
!emission possible distance
 150		CALL RANDOM4(SEED3,SEED4,RanDist)
		IF(AveExtinCoef.GT.0) THEN
		Dist0 = -1.0*LOG(1.0-RanDist)/AveExtinCoef
		ELSE
		Dist0 = 10.0*wholedist
		END IF
			

201		DV(1) = sin(Eta)*cos(Seta)*Dist0
		DV(2) = sin(Eta)*sin(Seta)*Dist0
		DV(3) = cos(Eta)*Dist0
		VT(1) = V0(1) + DV(1)
		VT(2) = V0(2) + DV(2)
		VT(3) = V0(3) + DV(3)
			
!looking for an intersecting wall surface if any
		CALL FindArrivalWall(DistW,Dist0,V0,VT,VW,DV,Seta,Eta,
     &		  IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NWT,iray,ierr,iben,
     &		  XYZ,FanWei)
		IF(ierr.GT.0) THEN
		ireset=ireset+1
		goto 123
		END IF
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! found no surface (idwall=0), or it is too far, (idwall<0)
!	
	IF(IDWall.LE.0) THEN
	  CALL FindArrivalElem(iray,MatchID,VT,N2D,N3D,NDOF,
     &	  NB,NPBG,NP,NE,MP,ME,NFACE,NKIND,XYZ,ELEM,FanWei)
	IF(MatchID.LE.0) THEN
	ireset=ireset+1
	goto 123
        END IF
			
!particle scattering for participating media
!for non-scattering gas, please set Albedo = 0
	  CALL RANDOM4(SEED5,SEED6,RanAlbedo)
!absorbed
	  IF(RanAlbedo.GE.AveAlbedo) THEN
            NumAbsorb(MatchID) = NumAbsorb(MatchID) + 1
!scattered
	  ELSE
	V0 = VT
        CALL RANDOM4(SEED1,SEED2,RanEta1)
	CALL RANDOM4(SEED1,SEED2,RanSeta1)
				
!Isotropic scattering				
!Eta1 = ACOS(1-2*RanEta1)
!Seta1 = RanSeta1*2.0*PAI
!reset travel
	CALL ScatterAngleLAS(Seta,Eta,RanSeta1,RanEta1,
     &		EType(MatchID),scabeta(MatchID))

!CALL RANDOM4(SEED3,SEED4,RanDist)

	GOTO 150	!distance reset

        END IF
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intersect with a surface patch
!
        ELSE

	CALL RANDOM4(SEED5,SEED6,RanWGray)
	JRD = VIndex(IDWALL+NE)
!absorbed
        IF(JRD.GT.0.AND.RanWGray.LE.1.0-Albedo(IDWALL+NE)
     &	.and.DistW.GT.0.0) THEN
			
!absorbed by wall
c				FinalPos(IRay,1) = VW(1)
c				FinalPos(IRay,2) = VW(2)
c				FinalPos(IRay,3) = VW(3)
	NumAbsorb(IDWall+NE) = NumAbsorb(IDWall+NE) + 1
			
!periodic boundary wall
	ELSE IF(WGrayDeg(IDWALL).LT.-0.9) THEN	
				
	IPair = ABS(WGrayDeg(IDWALL))
	Dist0 = Dist0 - DistW
	V0(1) = VW(1) + (CTR(1,NE+IPair)-CTR(1,NE+IDWALL))
	V0(2) = VW(2) + (CTR(2,NE+IPair)-CTR(2,NE+IDWALL))
	V0(3) = VW(3) + (CTR(3,NE+IPair)-CTR(3,NE+IDWALL))
	Vlast = V0
	GOTO 201
!reflected
	ELSE
!reflected by wall
	IF(DistW.GT.0.0) THEN
!performed one true step of reflection
	Dist0 = Dist0 - DistW
	Vlast = V0
        V0 = VW
	ELSE
!the present step is idle, resume it
	V0 = Vlast
        END IF

!Isotropic reflecting face
	CALL RANDOM4(SEED1,SEED2,RanSeta)
	CALL RANDOM4(SEED1,SEED2,RanEta)
	NV(1) = WNV(1,IDWall)
	NV(2) = WNV(2,IDWall)
        NV(3) = WNV(3,IDWall)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WNV is the outward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
c			CALL BuildZuoBiao(NV,TV1,TV2,3)
c				
c			CALL ZuobiaoRot(Seta,Eta,NV(1),NV(2),NV(3),
c     &			TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
				
	CALL ReflectAngleMC(Seta,Eta,RanSeta,RanEta,NV,
     &					EType(IDWALL+NE))
	GOTO 201				
	END IF
        END IF

 500	CONTINUE




			
C-----------------------------
C	Calcuate R-E-A-D values
C-----------------------------
        IIT = NE + IW
	IF(StackOpt(1:4).EQ.'BAND') GOTO  699
		
!----------------------------------Stack Method 1------------------------------------------
!---
!---	Stack it into a 1-D array, anologous to a Lower-Trianglar-Array (NRD0,NRD0)/2
!---
!------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************!
        IF(IFlagGp.NE.1) THEN		!not grouped

	IRD = VIndex(IIT)
	ratio = 1.0*NumAbsorb(IIT)/NRAY1
!PROP(IRD) = WArea(IW)*WGrayDeg(IW)
	RDValue(RDIndex(IRD)+IRD) =ratio*PROP(IIT)
	selfratio=ratio
	sumofratio = ratio
	DO 600 J = 1,IIT-1
	 JRD = VIndex(J)
	 IF(JRD.GT.0) THEN
	 ratio1 = 1.0*NumAbsorb(J)/NRAY1
	 sumofratio = sumofratio + ratio1
	 ratio1 = ratio1*PROP(IIT)
	 ratio2 = RDValue(RDIndex(IRD)+JRD)
	 CALL  SymmTreatRD(ratio1,ratio2,PROP(IIT),PROP(J))
	 RDValue(RDIndex(IRD)+JRD) = ratio1
         ELSE IF(NumAbsorb(J).GT.0) THEN
	WRITE(IUT0,*) 'ID = ',J
	WRITE(IUT0,*) 'An element received ray, '
	WRITE(IUT0,*) 'but does not appear in equation'
	STOP
         END IF	
		 
 600	CONTINUE
	
		
	DO 620 J = IIT+1,NRD0
	 
	 JRD = VIndex(J)
	 IF(JRD.GT.0) THEN			
!deposit to (JRD,IRD)
	 ratio = 1.0*NumAbsorb(J)/NRAY1
	 RDValue(RDIndex(JRD)+IRD) = ratio*PROP(IIT)
	 sumofratio = sumofratio + ratio
         ELSE IF(NumAbsorb(J).GT.0) THEN
	WRITE(IUT0,*) 'ID = ',J
	WRITE(IUT0,*) 'An element received ray, '
	WRITE(IUT0,*) 'but does not appear in equation'
	STOP
        END IF
 620	CONTINUE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************!
        ELSE	!Grouped
	IGP = VIndex(IIT)
	IRD = VIndex(IIT)
	sumofratio = 0
	selfratio =  1.0*NumAbsorb(IIT)/NRAY1
c		ID = IGP*NGP+IGP
c		RDValue(ID) =RDValue(ID)+PROP(IRD)*ratio
	DO 640 J = 1,NE+NW
	 IF(NumAbsorb(J).EQ.0) CYCLE
	 JGP = VIndex(J)
	 ratio1 = 1.0*NumAbsorb(J)/NRAY1
	 sumofratio = sumofratio + ratio1
	 IF(JGP.GT.0) THEN
	 ID = (IGP-1)*NGP+JGP
	 RDValue(ID)=RDValue(ID)+PROP(IIT)*ratio1
         ELSE
	WRITE(IUT0,*) 'ID = ',J, 'JGP=',JGP
	WRITE(IUT0,*) 'An element received ray, '
	WRITE(IUT0,*) 'but does not appear in equation'
	STOP
         END IF
		 
 640	CONTINUE
        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c633	  WRITE(IUT0,'(a,E16.8)') '     Sum of values=', sumofratio
c	  WRITE(IUT0,'(a,E16.8)') '     Self absorption=',selfratio

	  GOTO 1000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------Method-2------------------------------------------
!---
!---	Stack view-factors in a 1-D array, with band NNN and Element Index in RDId(:)
!---
!------------------------------------------------------------------------------------
	
 699	IRD = VIndex(IIT)
	RDIndex(IRD) = IMRD
	IMRD = IMRD + 1
		ratio=1.0*NumAbsorb(IIT)/NRAY1
		RDId(IMRD) = IRD
		!PROP(IRD) = WArea(IW)*WGrayDeg(IW)
		RDValue(IMRD) = PROP(IIT)*ratio
		sumofratio = ratio

		IF(IMRD.GT.MRD-MRD/NRD0) THEN
			WRITE(IUT0,*) 'Error in WallRayTracer():'
	WRITE(IUT0,*) 'Stack maybe overflow, please enlarge MRD= '
     &			,MRD
			STOP
		END IF
	

		
		DO 700 J = 1,IIT-1
		 ratio = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio

		 JRD = VIndex(J)
		 IF(ratio.GT.0.0) THEN
			IF(JRD.GT.0) THEN
				IMRD = IMRD + 1
				RDId(IMRD) = JRD
				RDValue(IMRD) = ratio*PROP(IIT)
			ELSE 
			WRITE(IUT0,*) 'ID = ',J
			WRITE(IUT0,*) 'An element received ray, '
			WRITE(IUT0,*) 'but does not appear in equation'
				STOP
			END IF
		 END IF	

700		CONTINUE
		
		DO 710 J = IIT+1,NRD0
		 ratio = 1.0*NumAbsorb(J)/NRAY1
		 sumofratio = sumofratio + ratio

		 JRD = VIndex(J)
		 IF(ratio.GT.0.0) THEN
			IF(JRD.GT.0) THEN
				IMRD = IMRD + 1
				RDId(IMRD) = JRD
				RDValue(IMRD) = ratio*PROP(IIT)
			ELSE IF(NumAbsorb(J).GT.0) THEN
			WRITE(IUT0,*) 'ID = ',J
			WRITE(IUT0,*) 'An element received ray, '
			WRITE(IUT0,*) 'but does not appear in equation'
			STOP
			END IF	
		  END IF
710		CONTINUE


999		continue
c		WRITE(IUT0,*) '   Total Record Num=', IMRD-RDIndex(IRD)
c		WRITE(IUT0,*) '   Sum of values=',sumofratio	
		
1000  CONTINUE

	WRITE(IUT0,'(a,I10)') '    Total tracing rays = ',NTRay-NTRAY0
	WRITE(IUT0,'(a,I5)') '     Expanded-search times =', iben
	WRITE(IUT0,'(a,I5)') '     Ray reset times       =', ireset


	DEALLOCATE(NumAbsorb,WK2)


	RETURN

!-------> end of WallRayTracer()

	END 
C======================================================================C
	SUBROUTINE  WallRayTracerInhom(IUT0,VMIN,NRAY,NP,NE,NW,MP,ME,NWT,
     2		NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD,MID,IRD,
     3		IMRD,dseg,AveWArea,AveVol,NTRay,
     4		AvePAbsorb,AvePScatter,AveGAbsorb,NGP,
     5		StackOpt,IFlagGp,ModCal,GlobeCTR,FanWei,XYZ,ELEM)
C======================================================================C

	USE MODULE_RADSXF, ONLY : FaceNum,NNPE,NNPS,CNTIVITY,ElemKind
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE MODULE_RADSXF, ONLY : GAbsorb,PAbsorb,PScatter,WGrayDeg
	USE MODULE_RADSXF, ONLY : WElem,WABCD,WNV,WArea,Volume,WNEB
	USE MODULE_RADSXF, ONLY : RDId,RDN1,RDN2,RDIndex,RDValue
	USE MODULE_RADSXF, ONLY : EType,VIndex,PROP,Albedo,CTR,scabeta
	USE MODULE_RADGROUP, ONLY : GpCNum,GpProp
	USE module_radwork, ONLY : NumAbsorb => WK1
	USE module_radwork, ONLY : VWK3D,WK2


	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)
	
	INTEGER ELEM(N3D,ME)
	REAL*8	XYZ(NDOF,MP),GlobeCTR(3),FanWei(2,NDOF)

	REAL*8 AveExtinCoef,AvePAbsorb,AvePScatter,AveGAbsorb,AveWArea
	CHARACTER*20 StackOpt

!NW:	number of wall surfaces
!NWT:	number of total surface patches, (NW+1:NWT) are the 
!		surfaces that are not specified as wall

	!Private
	REAL*8	NODE(NDOF,N3D),NODEP(NDOF,N3D),fun(N3D)
      REAL*8	V0(3),VT(3),VW(3),DV(3),VLast(3)
	REAL*8	NV(3),TV1(3),TV2(3)
	INTEGER NTRay
	

	!ALLOCATE(NumAbsorb(NRD0),WK2(NRD0),VWK3D(NB,NB,NB))
	ALLOCATE(NumAbsorb(NRD0),WK2(NRD0))

	WK2 = 0	

	!CALL GetBlockExtin(IUT0,NB,NPBG)


	AveExtinCoef = AveGAbsorb + AvePAbsorb + AvePScatter
	iben=0
	ireset=0
	wholedist = 0
	DO  I=1,3
	wholedist = max(wholedist,FanWei(2,I)-FanWei(1,I))
	END DO


	SEEDX = 0		!position
	SEEDY = 0
	SEEDZ = 0
	SEED1 = 0		!angle
	SEED2 = 0
	SEED3 = 0		!distance
	SEED4 = 0		
	SEED5 = 0		!scattering albedo
	SEED6 = 0
	NTRAY0=NTRAY
		    
	DO 1000 IW = 1, NW
		
		
		IF(VIndex(IW+NE).LT.1) CYCLE
		
!IF(MOD(VIndex(IW+NE),100).EQ.0) WRITE(IUT0,*) 'IW = ',IW
			!boundary condition element
		
		


C-----------------------
C	Pre-treat
C-----------------------
    

		DO 100 J =1,WElem(1,IW)
			IP = WElem(J+1,IW)
			DO 100 IDOF = 1,NDOF
			  NODE(IDOF,J) = XYZ(IDOF,IP)
100		CONTINUE



				
C--------------------
C	Ray Tracing of wall
C--------------------
		NumAbsorb = 0

		RayRatio=MIN(MAX(1.0d0, AveWArea/WArea(IW)),10.0d0)
		IF(IFlagGp.EQ.1) THEN
		  IGP=VIndex(NE+IW)
		  !IF(GpCNum(IGP).GT.1) RayRatio=1.0/GpCNum(IGP)
		  IF(GpCNum(IGP).GT.1) 
     &		  RayRatio=MAX(PROP(NE+IW)/GpPROP(IGP),0.1)
		END IF
		NRAY1=NRAY*RayRatio
		NTRay=NTRay+NRAY1

		!WRITE(IUT0,'("IW=",I5,", NRAY=",I6) ') IW,NRAY1
		IF(MOD(IW,100).EQ.1) 
     &		  WRITE(IUT0,'("     IW=",I8," /[",I8,"]") ') IW-1,NW


	    
		DO 500 IRAY = 1, NRAY1

c			if(mod(iray,5000).EQ.0.OR.IRAY.EQ.NRAY) 
c     & 			WRITE(IUT0,*) '   IRAY = ',IRAY
		
			!emission start position V0
123			CALL RANDOM1(SEEDX,fun(1))	
			CALL RANDOM2(SEEDY,fun(2))
			CALL RANDOM3(SEEDZ,fun(3))

			IF(WElem(1,IW).EQ.3) THEN
				CALL Tri3Pos(NODE,N3D,V0,fun,NDOF)

			ELSE IF(WElem(1,IW).EQ.4) THEN
				CALL Qua4Pos(NODE,N3D,V0,fun,NDOF)

			ELSE IF(WElem(1,IW).GT.4) THEN
		WRITE(IUT0,*) '!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!'
		WRITE(IUT0,*) 'IW=',IW,'NNPS=',WElem(1,IW)
		WRITE(IUT0,*) 'Shape function for this kind of face is'
		WRITE(IUT0,*) 'not available!'
		STOP
		END IF

		VLast = V0

		!emission angle
		CALL RANDOM4(SEED1,SEED2,RanSeta)
		CALL RANDOM4(SEED1,SEED2,RanEta)
			
			
			
c			Seta = RanSeta*2.0*PAI
c			Eta = ACOS(SQRT(1-RanEta))
	
		NV(1) = WNV(1,IW)
		NV(2) = WNV(2,IW)
		NV(3) = WNV(3,IW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WNV is the outward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
c		CALL BuildZuoBiao(NV,TV1,TV2,3)
c		CALL ZuobiaoRot(Seta,Eta,NV(1),NV(2),NV(3),
c     &		TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
			

	CALL EmitAngleWall(Seta,Eta,RanSeta,RanEta,NV,
     &			EType(IW+NE))
!c------------this is a loop for ray emission and tracing
!emission possible distance
 150	CALL RANDOM4(SEED3,SEED4,RanDist)
	IE = WNEB(IW)
	ExtinCoef = GAbsorb(IE)+PAbsorb(IE)+PScatter(IE)
	IF(ExtinCoef.GT.0) THEN
	Dist0 = -1.0*LOG(1.0-RanDist)/ExtinCoef
	ELSE
	Dist0 = -1.0*LOG(1.0-RanDist)/AveExtinCoef
	END IF
	IE1 = IE
	Vlast = V0
!looking for a final volume or intersecting wall surface if any
	MarkDaoda = 0
	IDWALL = 0
	DO 333 WHILE(MarkDaoda.EQ.0)
	IF(Seta.LT.0) Seta = Seta+2*PAI
	ExtinCoef1 = GAbsorb(IE1)+PAbsorb(IE1)+PScatter(IE1)
	IF(ExtinCoef1.GT.0) THEN
	IF(ExtinCoef.EQ.0) THEN
	Dist0 = -1.0*LOG(1.0-RanDist)/ExtinCoef1
	ExtinCoef = ExtinCoef1
	ExtinCoef0 = ExtinCoef
	ELSE
	ExtinCoef0 = ExtinCoef
	ExtinCoef = ExtinCoef1
        Dist0 = Dist0*ExtinCoef0/ExtinCoef
	ENDIF
	DistGo = min(Dist0,dseg)	!possible travel distance
	ELSE
	DistGo = dseg
	END IF
!looking for an intersecting wall surface if any
	IF(IDWALL.EQ.0) THEN
	  DV(1) = sin(Eta)*cos(Seta)*Dist0		!*wholedist*1.e5	
	  DV(2) = sin(Eta)*sin(Seta)*Dist0
	  DV(3) = cos(Eta)*Dist0

	  VT(1) = V0(1) + DV(1)
	  VT(2) = V0(2) + DV(2)
	  VT(3) = V0(3) + DV(3)

	  CALL FindArrivalWall(DistW,Dist0,V0,VT,VW,DV,Seta,Eta,
     &		IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NWT,iray,ierr,iben,
     &	        XYZ,FanWei)
				
				
	 IF(IDWALL.EQ.0) THEN
c					!Error abort
c					WRITE(IUT0,*) 'ERROR in GasRayTracerInhom():'
c					WRITE(IUT0,*) 'iray =',iray
c					WRITE(IUT0,*) v0,dv
c					WRITE(IUT0,*) 'The ray intersect no wall'
c    					WRITE(IUT0,*) 'please check the mesh or code !!!'
c					STOP
c					ireset=ireset+1
c					goto 123
           DistW=Dist0+DistGo+wholedist
	 END IF
				
	ELSE
	DistW = SQRT((VW(1)-V0(1))*(VW(1)-V0(1))+ 
     &	(VW(2)-V0(2))*(VW(2)-V0(2)) +
     &	(VW(3)-V0(3))*(VW(3)-V0(3))	)
	ENDIF
	DV(1) = sin(Eta)*cos(Seta)*DistGo			!Dist0
	DV(2) = sin(Eta)*sin(Seta)*DistGo			!Dist0
	DV(3) = cos(Eta)*DistGo					!Dist0

	VT(1) = V0(1) + DV(1)
	VT(2) = V0(2) + DV(2)
	VT(3) = V0(3) + DV(3)

!march on until encountering an element with different property
!	(i.e., until the optical-extinction coef. is changed)
! ExtinCoef keeps unchanged, but Dist0,DistW and V0,VT may be
! changed
        imarchon=0
c				IF(IDWALL.NE.0.AND.DistW.GT.DistGo) THEN
c				  CALL MarchOn(iray,IE1,IDWALL,ExtinCoef1,V0,VT,DistW,
c     &				Dist0,DistGo,seta,eta,N2D,N3D,NDOF,NB,NPBG,
c     &				NP,NE,MP,ME,NFACE,NKIND,imarchon)
!write(*,*) 'imarchon=',imarchon
c				END IF
										
!arrive at a volume element
!IF(DistW.GT.DistGo.AND.imarchon.LE.2) THEN
	IF(DistW.GT.DistGo.AND.imarchon.LE.2) THEN
	  CALL FindArrivalElem(iray,MatchID,VT,N2D,N3D,NDOF,
     &	  NB,NPBG,NP,NE,MP,ME,NFACE,NKIND,XYZ,ELEM,FanWei)

        IF(MatchID.LE.0) THEN
!Error abort
c					WRITE(IUT0,*) 'ERROR in GasRayTracerInhom():'
c					WRITE(IUT0,*) 'iray =',iray
c					WRITE(IUT0,*) v0,dv
c					WRITE(IUT0,*) 'The ray goes to unknown place'
c    					WRITE(IUT0,*) 'please check the mesh or code !!!'
c					STOP
	ireset=ireset+1
	goto 123
        END IF
	ELSE
        MatchID=WNEB(ABS(IDWALL))
	END IF
C-----------------------------2006/04/01---------------------------------------C
c
c				DV(1) = sin(Eta)*cos(Seta)*DistGo
c				DV(2) = sin(Eta)*sin(Seta)*DistGo
c				DV(3) = cos(Eta)*DistGo
c
c				VT(1) = V0(1) + DV(1)
c				VT(2) = V0(2) + DV(2)
c				VT(3) = V0(3) + DV(3)
c
c				
c				!looking for an intersecting wall surface if any
c				IF(IDWALL.EQ.0) THEN
c
c				 CALL FindArrivalWall(DistW,Dist0,V0,VT,VW,DV,Seta,Eta,
c     &				IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NWT,iray,ierr,iben)
c				ELSE
c				  DistW = SQRT((VW(1)-V0(1))*(VW(1)-V0(1)) + 
c     &						   (VW(2)-V0(2))*(VW(2)-V0(2)) +
c     &						   (VW(3)-V0(3))*(VW(3)-V0(3))	)
c				END IF
c
c
c				!arrive at a volume element
c				CALL FindArrivalElem(iray,MatchID,VT,N2D,N3D,NDOF,
c     &					NB,NPBG,NP,NE,MP,ME,NFACE,NKIND)
c	
c
c				IF(MatchID.LT.0.AND.IDWALL.LT.0) THEN
c					!reset
cc					WRITE(IUT0,*) 'ERROR in WallRayTracerInhom():'
cc					WRITE(IUT0,*) 'iray =',iray
cc					WRITE(IUT0,*) v0,dv
cc					WRITE(IUT0,*) 'The ray goes to unknown place'
cc   					WRITE(IUT0,*) 'please check the mesh or code !!!'
cc					STOP
c					ireset=ireset+1
c					goto 123
c				END IF
C--------------------------------2006/04/01---------------------------------------C
c	if(iw.eq.809) write(*,*) iray,idwall,matchid,imarchon,
C    &				distw,distgo,dist0,extincoef,extincoef1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	intersect a surface patch
!
	IF(imarchon.EQ.3.OR. 
     &	(IDWALL.NE.0.AND.DistW.LE.DistGo))THEN
!	IF(IDWALL.NE.0.AND.DistW.LE.DistGo)THEN
!intersect with the wall surface
	IDWALL = ABS(IDWALL)
	CALL RANDOM4(SEED5,SEED6,RanWGray)
	JRD = VIndex(IDWALL+NE)
!absorbed
        IF(JRD.GT.0.AND.
     &	RanWGray.LE.1.0-Albedo(IDWALL+NE)) THEN
	NumAbsorb(IDWall+NE) = NumAbsorb(IDWall+NE)+1
	MarkDaoda = 1
!periodic boundary wall
	ELSE IF(WGrayDeg(IDWALL).LT.-0.9) THEN	
!the value is the ID of wall-surf	
	IPair = ABS(WGrayDeg(IDWALL))
	Dist0 = Dist0 - DistW
						
	V0(1)=VW(1)+(CTR(1,NE+IPair)-CTR(1,NE+IDWALL))
	V0(2)=VW(2)+(CTR(2,NE+IPair)-CTR(2,NE+IDWALL))
	V0(3)=VW(3)+(CTR(3,NE+IPair)-CTR(3,NE+IDWALL))
	Vlast=V0
					
!reflected
	ELSE
!reflected by wall
	IF(DistW.GT.0.0) THEN
!performed one true step of reflection
	IE2 = WNEB(IDWALL)
	ExtinCoef2 = GAbsorb(IE2)+PAbsorb(IE2)
     &	 + PScatter(IE2)
	IF(ExtinCoef2.GT.0) THEN
	ratio=2.*ExtinCoef/(ExtinCoef+ExtinCoef2)
	Dist0 = Dist0 - DistW*ratio
	END IF
	Vlast = V0
	V0 = VW
	IE1 = IE2
        ELSE
!the present step is idle, resume it
!this case probably arise for the ray in corner.
	V0 = Vlast
        END IF
!Isotropic reflecting face
	CALL RANDOM4(SEED1,SEED2,RanSeta)
	CALL RANDOM4(SEED1,SEED2,RanEta)
c				Seta = RanSeta*2.0*PAI
c				Eta = ACOS(SQRT(1.0-RanEta))

	NV(1) = WNV(1,IDWall)
	NV(2) = WNV(2,IDWall)
	NV(3) = WNV(3,IDWall)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WNV is the outward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	CALL BuildZuoBiao(NV,TV1,TV2,3)
c	CALL ZuobiaoRot(Seta,Eta,NV(1),NV(2),NV(3),
c     &	TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
!
        CALL ReflectAngleMC(Seta,Eta,RanSeta,RanEta,
     &	NV,EType(IDWALL+NE))
	IDWALL = 0		!reset
	END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	arrive at an element
!
!NO-Participating element, pass it WITHOUT distance-reduction
	ELSE IF(VIndex(MatchID).EQ.0) THEN	
!not a participating element
	V0 = VT
	IE2 = MatchID
	Vlast = V0
	IE1 = MatchID
!participating element, pass it WITH distance-reduction
	ELSE IF(Dist0.GT.1.01*DistGo) THEN	
!continue the process
	V0 = VT
	IE2 = MatchID
	ExtinCoef2 =GAbsorb(IE2)+PAbsorb(IE2)+PScatter(IE2)
	ratio = 2.0*ExtinCoef/(ExtinCoef+ExtinCoef2)
	Dist0 = Dist0 - DistGo*ratio
	Vlast = V0
	IE1 = MatchID
	
!final position, travel terminate at this element
	ELSE
!arrive at a volume
!particle scattering for participating media
!for non-scattering gas, please set Albedo = 0
	CALL RANDOM4(SEED5,SEED6,RanAlbedo)
	JRD = VIndex(MatchID)
!absorbed
	IF(RanAlbedo.GE.Albedo(MatchID)) THEN
	NumAbsorb(MatchID) = NumAbsorb(MatchID) + 1
	MarkDaoda = 1
!scattered
	ELSE
	V0 = VT			!New emission position
	CALL RANDOM4(SEED1,SEED2,RanEta1)
	CALL RANDOM4(SEED1,SEED2,RanSeta1)
!Scattering angle
!Eta1 = ACOS(1-2*RanEta1)
!Seta1 = RanSeta1*2.0*PAI
!reset travel
	CALL ScatterAngleLAS(Seta,Eta,RanSeta1,
     &	RanEta1,EType(MatchID),scabeta(MatchID))

	CALL RANDOM4(SEED3,SEED4,RanDist)
	IE1 = MatchID
	ExtinCoef = GAbsorb(IE1)+PAbsorb(IE1)
     &		+PScatter(IE1)
	Dist0 = -1.0*LOG(1.0-RanDist)/ExtinCoef		
	Vlast = V0
	IDWALL = 0		!reset
	END IF
	END IF
 333	CONTINUE
 500	CONTINUE
C-----------------------------
C	Calcuate R-E-A-D values
C-----------------------------
        IIT = NE + IW
	IF(StackOpt(1:4).EQ.'BAND') GOTO  699
!----------------------------------Stack Method 1------------------------------------------
!---
!---	Stack it into a 1-D array, anologous to a Lower-Trianglar-Array (NRD0,NRD0)/2
!---
!------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************!
        IF(IFlagGp.NE.1) THEN		!not grouped
	IRD = VIndex(IIT)
	ratio = 1.0*NumAbsorb(IIT)/NRAY1
!PROP(IRD) = WArea(IW)*WGrayDeg(IW)
	RDValue(RDIndex(IRD)+IRD) =ratio*PROP(IIT)
	selfratio =ratio
	sumofratio = ratio
	DO 600 J = 1,IIT-1
        JRD = VIndex(J)
	IF(JRD.GT.0) THEN
	 ratio1 = 1.0*NumAbsorb(J)/NRAY1
	 sumofratio = sumofratio + ratio1
	 ratio1 = ratio1*PROP(IIT)
	 ratio2 = RDValue(RDIndex(IRD)+JRD)
	 CALL  SymmTreatRD(ratio1,ratio2,PROP(IIT),PROP(J))
	 RDValue(RDIndex(IRD)+JRD) = ratio1
         ELSE IF(NumAbsorb(J).GT.0) THEN
	 WRITE(IUT0,*) 'ID = ',J
	 WRITE(IUT0,*) 'An element received ray, '
	 WRITE(IUT0,*) 'but does not appear in equation'
	STOP
        END IF	
 600    CONTINUE
	DO 620 J = IIT+1,NRD0
	 JRD = VIndex(J)
	 IF(JRD.GT.0) THEN			
	 ratio = 1.0*NumAbsorb(J)/NRAY1
	 RDValue(RDIndex(JRD)+IRD) = ratio*PROP(IIT)
	 sumofratio = sumofratio + ratio
         ELSE IF(NumAbsorb(J).GT.0) THEN
	 WRITE(IUT0,*) 'ID = ',J
	 WRITE(IUT0,*) 'An element received ray, '
	 WRITE(IUT0,*) 'but does not appear in equation'
	 STOP
	 END IF
 620     CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************!
         ELSE		!Grouped
	IGP = VIndex(IIT)
	IRD = VIndex(IIT)
	sumofratio = 0
	selfratio =  1.0*NumAbsorb(IIT)/NRAY1
c		ID = IGP*NGP+IGP
c		RDValue(ID) =RDValue(ID)+PROP(IRD)*ratio
        DO 640 J = 1,NE+NW
	 IF(NumAbsorb(J).EQ.0) CYCLE
	 JGP = VIndex(J)
	 ratio1 = 1.0*NumAbsorb(J)/NRAY1
	 sumofratio = sumofratio + ratio1
	 IF(JGP.GT.0) THEN
	 ID = (IGP-1)*NGP+JGP
	 RDValue(ID)=RDValue(ID)+PROP(IIT)*ratio1
         ELSE
         WRITE(IUT0,*) 'ID = ',J, 'JGP=',JGP
		WRITE(IUT0,*) 'An element received ray, '
		WRITE(IUT0,*) 'but does not appear in equation'
		STOP
	 END IF
		 
 640	CONTINUE


	  END IF
!
	  GOTO 1000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------Method-2------------------------------------------
!---
!---	Stack view-factors in a 1-D array, with band NNN and Element Index in RDId(:)
!---
!------------------------------------------------------------------------------------
	
 699	IRD = VIndex(IIT)
	RDIndex(IRD) = IMRD
	IMRD = IMRD + 1
	ratio=1.0*NumAbsorb(IIT)/NRAY1
	RDId(IMRD) = IRD
!PROP(IRD) = WArea(IW)*WGrayDeg(IW)
	RDValue(IMRD) = PROP(IIT)*ratio
	sumofratio = ratio
	IF(IMRD.GT.MRD-MRD/NRD0) THEN
	WRITE(IUT0,*) 'Error in WallRayTracer():'
	WRITE(IUT0,*) 'Stack maybe overflow, please enlarge MRD= '
     &			,MRD
	STOP
	END IF
	DO 700 J = 1,IIT-1
	 ratio = 1.0*NumAbsorb(J)/NRAY1
	 sumofratio = sumofratio + ratio
	 JRD = VIndex(J)
	 IF(ratio.GT.0.0) THEN
	IF(JRD.GT.0) THEN
	IMRD = IMRD + 1
	RDId(IMRD) = JRD
	RDValue(IMRD) = ratio*PROP(IIT)
        ELSE 
	WRITE(IUT0,*) 'ID = ',J
	WRITE(IUT0,*) 'An element received ray, '
	WRITE(IUT0,*) 'but does not appear in equation'
	STOP
	END IF
        END IF	
 700    CONTINUE
	DO 710 J = IIT+1,NRD0
	 ratio = 1.0*NumAbsorb(J)/NRAY1
	 sumofratio = sumofratio + ratio
	 JRD = VIndex(J)
	 IF(ratio.GT.0.0) THEN
	IF(JRD.GT.0) THEN
	 IMRD = IMRD + 1
	 RDId(IMRD) = JRD
	 RDValue(IMRD) = ratio*PROP(IIT)
        ELSE IF(NumAbsorb(J).GT.0) THEN
	WRITE(IUT0,*) 'ID = ',J
	WRITE(IUT0,*) 'An element received ray, '
	WRITE(IUT0,*) 'but does not appear in equation'
	STOP
        END IF	
	END IF
 710	CONTINUE
 999	continue
 1000   CONTINUE

	WRITE(IUT0,'(a,I10)') '    Total tracing rays = ',NTRay-NTRAY0
	WRITE(IUT0,'(a,I5)') '     Expanded-search times =', iben
	WRITE(IUT0,'(a,I5)') '     Ray reset times       =', ireset

!DEALLOCATE(NumAbsorb,WK2,VWK3D)
	DEALLOCATE(NumAbsorb,WK2)

	RETURN

!-------> end of WallRayTracerInhom()
	END 
C======================================================================C
	SUBROUTINE	FindArrivalElem(iray,MatchID,VT,
     &		N2D,N3D,NDOF,NB,NPBG,NP,NE,MP,ME,NFACE,NKIND,
     &		XYZ,ELEM,FanWei)
C======================================================================C
	
	USE MODULE_RADSXF, ONLY : ElemKind,FaceNum,NNPE,NNPS,CNTIVITY
	USE MODULE_RADSXF, ONLY : Volume,CTR
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE module_radwork, ONLY : WK2
	

	IMPLICIT REAL*8(A-H,O-Z)

	INTEGER ELEM(N3D,ME)
	REAL*8	XYZ(NDOF,MP),FanWei(2,NDOF)
	
	!Private
	REAL*8  VT(3),P1(3),P2(3),P3(3)
	REAL*8	FaceCTR(NDOF),FaceNV(NDOF),ABCD(4)
	REAL*8	Vertex(NDOF,N3D),FaceVertex(NDOF,N2D),vmemo(8)
	INTEGER LocalID(3,3),OBJ(5000)
	

	FDA = 0.1
	!WK2 = 0
	!OBJ = 0


	iwk = 0
	dx = (FanWei(2,1)-FanWei(1,1))/NB
	dy = (FanWei(2,2)-FanWei(1,2))/NB
	dz = (FanWei(2,3)-FanWei(1,3))/NB
	kx = 1
	ky = 1
	kz = 1

	jx = (VT(1)-FanWei(1,1))/dx+1
	jy = (VT(2)-FanWei(1,2))/dy+1
	jz = (VT(3)-FanWei(1,3))/dz+1

	MatchID=0

	IF((jx.GT.NB).OR.(jy.GT.NB).OR.(jz.GT.NB)) THEN
		MatchID = 0
		RETURN
	ENDIF

	IF((jx.LT.1).OR.(jy.LT.1).OR.(jz.LT.1)) THEN
		MatchID = 0
		RETURN
	END IF

	imax = (VT(1)-FanWei(1,1)+FDA*dx)/dx+1
	imin = (VT(1)-FanWei(1,1)-FDA*dx)/dx+1
	jmax = (VT(2)-FanWei(1,2)+FDA*dy)/dy+1
	jmin = (VT(2)-FanWei(1,2)-FDA*dy)/dy+1
	kmax = (VT(3)-FanWei(1,3)+FDA*dz)/dz+1
	kmin = (VT(3)-FanWei(1,3)-FDA*dz)/dz+1

	imax = min(NB,imax)
	imin = max(1,imin)
	jmax = min(NB,jmax)
	jmin = max(1,jmin)
	kmax = min(NB,kmax)
	kmin = max(1,kmin)

	

	!stack objects
	num = GBlockNum(jx,jy,jz)
	ib = 0
	DO 100 I =1,NUM
		
		ID = GBlockIndex(jx,jy,jz)+I 
		IE = GBlockId(ID)
		IF(WK2(IE).EQ.0) THEN
			ib = ib + 1
			OBJ(IB) = IE
			WK2(IE) = IB
		END IF

100	CONTINUE

	

	DO 220 ix = imin, imax
	  DO 220 iy = jmin, jmax
	    DO 220 iz = kmin, kmax

		   num = GBlockNum(ix,iy,iz)
		   DO 210 I =1,NUM
				ID = GBlockIndex(ix,iy,iz)+I 
				IE = GBlockId(ID)
				IF(WK2(IE).EQ.0) THEN
				ib = ib + 1
				OBJ(IB) = IE
				WK2(IE) = IB
				END IF
210		   CONTINUE
220	CONTINUE

	

	!seach the matching object
	NUMALL = IB
	!IF(numall.GT.200) write(*,*) 'IB=',ib
	IB = 0
	MatchID = 0
	DO 690 WHILE(IB.LT.NUMALL.AND.MatchID.LE.0)
			
			IB = IB + 1
			IE = OBJ(IB)
	
			!
			VVV = 0.0
			ID = ElemKind(IE)
			DO 650 JFACE = 1,FaceNum(ID)
			  DO 650 JNODE = 1, NNPS(ID,JFACE)-2
			IP1 = ELEM(CNTIVITY(ID,JFACE,JNODE+1),IE)
			IP2 = ELEM(CNTIVITY(ID,JFACE,JNODE+2),IE)
			IP3 = ELEM(CNTIVITY(ID,JFACE,1),IE)
              
			DO IDOF = 1,NDOF
				P1(IDOF) = XYZ(IDOF,IP1)
				P2(IDOF) = XYZ(IDOF,IP2)
				P3(IDOF) = XYZ(IDOF,IP3)
			END DO
		
			AX  = (VT(2)-P2(2))*(VT(3)-P3(3))
     &				-(VT(2)-P3(2))*(VT(3)-P2(3))
			AX  = (VT(1)-P1(1))*AX

			AY  = (VT(1)-P2(1))*(VT(3)-P3(3))
     &				-(VT(1)-P3(1))*(VT(3)-P2(3))
			AY  = (VT(2)-P1(2))*AY

			AZ  = (VT(1)-P2(1))*(VT(2)-P3(2))
     &				-(VT(1)-P3(1))*(VT(2)-P2(2))
			AZ  = (VT(3)-P1(3))*AZ

			VVV = VVV + ABS(AX - AY + AZ)
			
650		CONTINUE

		vvv = vvv/6.0
		dvolume=abs(VVV-Volume(IE))/Volume(IE)
		
		IF(dvolume.LT.1.0d-6) THEN
			MatchID = IE
c		ELSE IF(dvolume.LT.1.0d-4) THEN
c		
c			MatchID = -IE

		END IF

690	CONTINUE

	MatchID = ABS(MatchID)
	

1000	DO I = 1,NUMALL
		IE = OBJ(I)
		WK2(IE) = 0
	END DO


	RETURN


!-------> end of FindArrivalElem()

	END
C======================================================================C
	SUBROUTINE	FindArrivalElem_old(iray,MatchID,VT,
     &		N2D,N3D,NDOF,NB,NPBG,NP,NE,MP,ME,NFACE,NKIND,
     &		XYZ,ELEM,FanWei)
C======================================================================C
	
	USE MODULE_RADSXF, ONLY : ElemKind,FaceNum,NNPE,NNPS,CNTIVITY
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE module_radwork, ONLY : WK2
	

	IMPLICIT REAL*8(A-H,O-Z)

	INTEGER ELEM(N3D,ME)
	REAL*8	XYZ(NDOF,MP),FanWei(2,NDOF)


	!Private
	REAL*8  VT(3)
	REAL*8	FaceCTR(NDOF),FaceNV(NDOF),ABCD(4),CTR(NDOF)
	REAL*8	Vertex(NDOF,N3D),FaceVertex(NDOF,N2D),vmemo(8)
	INTEGER LocalID(3,3),OBJ(5000)
	

	FDA = 0.1
	!WK2 = 0
	!OBJ = 0


	iwk = 0
	dx = (FanWei(2,1)-FanWei(1,1))/NB
	dy = (FanWei(2,2)-FanWei(1,2))/NB
	dz = (FanWei(2,3)-FanWei(1,3))/NB
	kx = 1
	ky = 1
	kz = 1

	jx = (VT(1)-FanWei(1,1))/dx+1
	jy = (VT(2)-FanWei(1,2))/dy+1
	jz = (VT(3)-FanWei(1,3))/dz+1
	

	IF((jx.GT.NB).OR.(jy.GT.NB).OR.(jz.GT.NB)) THEN
		MatchID = -1
		RETURN
	ENDIF

	IF((jx.LT.1).OR.(jy.LT.1).OR.(jz.LT.1)) THEN
		MatchID = -1
		RETURN
	END IF

	imax = (VT(1)-FanWei(1,1)+FDA*dx)/dx+1
	imin = (VT(1)-FanWei(1,1)-FDA*dx)/dx+1
	jmax = (VT(2)-FanWei(1,2)+FDA*dy)/dy+1
	jmin = (VT(2)-FanWei(1,2)-FDA*dy)/dy+1
	kmax = (VT(3)-FanWei(1,3)+FDA*dz)/dz+1
	kmin = (VT(3)-FanWei(1,3)-FDA*dz)/dz+1

	imax = min(NB,imax)
	imin = max(1,imin)
	jmax = min(NB,jmax)
	jmin = max(1,jmin)
	kmax = min(NB,kmax)
	kmin = max(1,kmin)

	

	!stack objects
	num = GBlockNum(jx,jy,jz)
	ib = 0
	DO 100 I =1,NUM
		
		ID = GBlockIndex(jx,jy,jz)+I 
		IE = GBlockId(ID)
		IF(WK2(IE).EQ.0) THEN
			ib = ib + 1
			OBJ(IB) = IE
			WK2(IE) = IB
		END IF

100	CONTINUE

	

	DO 220 ix = imin, imax
	  DO 220 iy = jmin, jmax
	    DO 220 iz = kmin, kmax

		   num = GBlockNum(ix,iy,iz)
		   DO 210 I =1,NUM
				ID = GBlockIndex(ix,iy,iz)+I 
				IE = GBlockId(ID)
				IF(WK2(IE).EQ.0) THEN
				ib = ib + 1
				OBJ(IB) = IE
				WK2(IE) = IB
				END IF
210		   CONTINUE
220	CONTINUE

	

	!seach the matching object
	NUMALL = IB
	!IF(numall.GT.200) write(*,*) 'IB=',ib
	IB = 0
	MatchID = -1
	DO 690 WHILE(IB.LT.NUMALL.AND.MatchID.EQ.-1)
			
			IB = IB + 1
			IE = OBJ(IB)
			
			MarkNotIt = 0
			MarkInEdgeLine = 0
			JFACE = 0
			CTR = 0.0
			value1 = 0.0
		
			DO 600 INNPE = 1, NNPE(ElemKind(IE))
				DO 600 IDOF = 1, NDOF
				IPK = ELEM(INNPE,IE)
				Vertex(IDOF,INNPE) = XYZ(IDOF,IPK)
				CTR(IDOF)=CTR(IDOF)+XYZ(IDOF,IPK)
600			CONTINUE
			CTR = CTR/NNPE(ElemKind(IE))

	

			DO 670 WHILE((JFACE.LT.FaceNum(ElemKind(IE)))
     &			.AND.(MarkNotIt.EQ.0))
				
				JFACE =  JFACE + 1
				NumPoint = NNPS(ElemKind(IE),JFACE)
				
!get face center and normal vector (and equation)
				DO 610 IDOF = 1, NDOF
				  FaceCTR(IDOF) = 0
				  DO 605 K = 1, NumPoint
			IPK = CNTIVITY(ElemKind(IE),JFACE,K)
			FaceCTR(IDOF) = FaceCTR(IDOF)+Vertex(IDOF,IPK)
			FaceVertex(IDOF,K)= Vertex(IDOF,IPK)
605	              CONTINUE
		        FaceCTR(IDOF)=FaceCTR(IDOF)/NumPoint
			
610		CONTINUE


		Call FaceAreaEqu(FaceNV,ABCD,FaceVertex,NumPoint,
     &					 N2D,NDOF)
!cheak whether FaceNV is outward
		vx1 = FaceCTR(1)-CTR(1)
		vy1 = FaceCTR(2)-CTR(2)
		vz1 = FaceCTR(3)-CTR(3)
		CALL COSVALUE(vx1,vy1,vz1,
     &		  FaceNV(1),FaceNV(2),FaceNV(3),
     &		  value)		
!if inward, reverse it
		if(value.LT.-1.0E-10) FaceNV(:) = -FaceNV(:)
!check cosine of the angle
		vx1 = VT(1)-FaceCTR(1)
		vy1 = VT(2)-FaceCTR(2)
		vz1 = VT(3)-FaceCTR(3)
		vmode = SQRT(vx1*vx1+vy1*vy1+vz1*vz1)
		IF(vmode.GT.0) THEN
		vx1 = vx1/vmode
		vy1 = vy1/vmode
		vz1 = vz1/vmode
	        END IF
		CALL COSVALUE(vx1,vy1,vz1,
     &		  FaceNV(1),FaceNV(2),FaceNV(3),
     &		  value)
!write(*,*) 'facenum',Jface,numpoint,vx1,vy1,vz1
		IF(abs(value) < 1.0E-4)	THEN
			!			  
		        !	    /|     V1
			!		/|    /
			!		/|  /
			!		/|/  
			!		/|-------->FaceNV
			!		/|
		        !
			!the angle between v1 and FaceNV is about 90 deg,
			!the point is located at the surface of the elem
		MarkInEdgeLine = 1
		ELSE
	!IF((JFACE.GT.1).AND.(value*value1.lt.0.0)) MarkNotIt=1
		IF(value*value1.LE.-1.0E-10) MarkNotIt = 1
	!the node is not containned by the element
		value1 = value
		END IF
!value1 = value
		vmemo(JFACE) = value
	
 670		CONTINUE
	
		IF(MarkNotIt.EQ.0) THEN
		MatchID = IE
!GOTO 1000
c				vx1 = CTR(1)-VT(1)
c				vx2 = CTR(2)-VT(2)
c				vx3 = CTR(3)-VT(3)
c				PEDIST	= SQRT(vx1*vx1+vx2*vx2+vx3*vx3)
c				
c				IF(PEDIST.GT.MouseNebDist(IE)) THEN
c					MatchID(IP)= -1
c				END IF
		ENDIF
	
690	CONTINUE
	
1000	DO I = 1,NUMALL
	IE = OBJ(I)
	WK2(IE) = 0
	END DO

	RETURN

!-------> end of FindArrivalElem()

	END
C======================================================================C
	SUBROUTINE	MarchOn(iray,IDE,IDWALL,xshcoef,V0,VT,DistW,Dist0,
     &		DistGo,seta,eta,N2D,N3D,NDOF,NB,NPBG,NP,NE,MP,ME,NFACE,
     &		NKIND,ires,XYZ,ELEM,FanWei)
C======================================================================C
	
	USE MODULE_RADSXF, ONLY : ElemKind,FaceNum,NNPE,
     &						  NNPS,CNTIVITY
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE module_radwork, ONLY : WK2
	USE module_radwork, ONLY : ExtinCoef => VWK3D
	

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(NSEG=2)

	INTEGER ELEM(N3D,ME)
	REAL*8	XYZ(NDOF,MP),FanWei(2,NDOF)

	!Private
	REAL*8  VT(3),V0(3),ddd(3),DV(3),DistW,Dist0
	REAL*8  V1(3)
	REAL*8	Vertex(NDOF,N3D),FaceVertex(NDOF,N2D),vmemo(8)
	INTEGER LocalID(3,3),OBJ(100)
	


!------ires
!   = 0,		do nothing
!   = 1,        go directly to a new place, with dist0 reduction
!   = 2,        directly terminated at the desination volume cell
!   = 3,        directly terminated at the wall	- IDWall



	iwk = 0
	dx = (FanWei(2,1)-FanWei(1,1))/NB
	dy = (FanWei(2,2)-FanWei(1,2))/NB
	dz = (FanWei(2,3)-FanWei(1,3))/NB
	ddd(1) = dx
	ddd(2) = dy
	ddd(3) = dz

	jx = (V0(1)-FanWei(1,1))/dx+1
	jy = (V0(2)-FanWei(1,2))/dy+1
	jz = (V0(3)-FanWei(1,3))/dz+1
	
	VV1 = ExtinCoef(jx,jy,jz)
	
	dvalue=abs(vv1-xshcoef)
	if((vv1+xshcoef).gt.0.d0) dvalue=dvalue/(vv1+xshcoef)
	IF(dvalue.GT.0.05)	THEN
		!it is an inhomogeneous bucket
		ires=0
		DV(1) = sin(Eta)*cos(Seta)*DistGo
		DV(2) = sin(Eta)*sin(Seta)*DistGo
		DV(3) = cos(Eta)*DistGo
		VT(1) = V0(1) + DV(1)
		VT(2) = V0(2) + DV(2)
		VT(3) = V0(3) + DV(3)
		return
	END IF

c2-------------------------------------------------------------	
	
	ivnear=0
	iwnear=0

	IF(DistW.GT.Dist0) THEN
		ivnear=1
		DistM = Dist0
	ELSE
		iwnear=1	
		DistM = DistW
	END IF

	DV(1) = sin(Eta)*cos(Seta)*DistM
	DV(2) = sin(Eta)*sin(Seta)*DistM
	DV(3) = cos(Eta)*DistM

	IDOF = 1
	vmax = abs(DV(1)/ddd(1))
	vmax = max(vmax,abs(DV(2)/ddd(2)))
	vmax = max(vmax,abs(DV(3)/ddd(3)))

	MSEG = vmax*NSEG+1
	dseg = 1.0/MSEG
	
	i0 = 0
	j0 = 0
	k0 = 0
	ib = 0
	ib0 = ib

	ISEG = 0


c3-------------------------------------------------------------
	
	ISEG = 0
	Mark = 0
	VT = V0

	IF(VV1.EQ.0.d0) THEN
	  
	  DO 400 WHILE(ISEG.LT.MSEG+1.AND.Mark.EQ.0)
	
		ISEG = ISEG+1
		tseg=min((ISEG-1)*dseg,1.0d0)
		V1(1) = V0(1) + tseg*DV(1)
		V1(2) = V0(2) + tseg*DV(2)
		V1(3) = V0(3) + tseg*DV(3)

		i = (V1(1)-FanWei(1,1))/dx+1
		j = (V1(2)-FanWei(1,2))/dy+1
		k = (V1(3)-FanWei(1,3))/dz+1
		
		IF(i.GT.NB.OR.i.LT.1) CYCLE
		IF(j.GT.NB.OR.j.LT.1) CYCLE
		IF(k.GT.NB.OR.k.LT.1) CYCLE
		
		VV2=ExtinCoef(i,j,k)
		
		IF(VV2.GT.1.0d-3) THEN
			Mark = 1
		ELSE
			VT(:)=V1(:)				!continue
			
		END IF
		
400	  CONTINUE
		
	ELSE

	  DO 410 WHILE(ISEG.LT.MSEG+1.AND.Mark.EQ.0)
	
		ISEG = ISEG+1 
		tseg=min((ISEG-1)*dseg,1.0d0)
		V1(1) = V0(1) + tseg*DV(1)
		V1(2) = V0(2) + tseg*DV(2)
		V1(3) = V0(3) + tseg*DV(3)

		i = (V1(1)-FanWei(1,1))/dx+1
		j = (V1(2)-FanWei(1,2))/dy+1
		k = (V1(3)-FanWei(1,3))/dz+1
		
		IF(i.GT.NB.OR.i.LT.1) CYCLE
		IF(j.GT.NB.OR.j.LT.1) CYCLE
		IF(k.GT.NB.OR.k.LT.1) CYCLE
		
		VV2=ExtinCoef(i,j,k)
		dvalue=abs(vv2-vv1)/vv1

		IF(dvalue.GT.0.05) THEN
			Mark = 1
		ELSE
			VT(:)=V1(:)			!continue
		END IF

410	  CONTINUE

	END IF

	
	
	IF(Mark.EQ.1) THEN	!encountered different elements
		
		dvalue= (vt(1)-v1(1))*(vt(1)-v1(1))+
     &			(vt(2)-v1(2))*(vt(2)-v1(2))+
     &			(vt(3)-v1(3))*(vt(3)-v1(3))		
		dvalue=sqrt(dvalue)
		IF(dvalue/DistM.LE.5.0d-3) THEN
			Mark=0		!very close
		ELSE
			DistM=DistM-dvalue
		END IF
	END IF

	
	IF(Mark.EQ.0) THEN
		!may go to the destination directly
		IF(iwnear.EQ.1) THEN	!arrived at wall
			ires=3
			DistW =0.0
			IF(xshcoef.GT.0.d0) Dist0=Dist0-DistW
						
		ELSE			!terminate at an element
			ires=2
			if(xshcoef.GT.0.d0) Dist0 = 0.0
			DistW = DistW-DistM
			DistGo=0.0
		END IF

	ELSE

		ires = 1
		V0 = VT

		DistW=DistW-DistM
		IF(xshcoef.GT.0.d0) Dist0=Dist0-DistM
		DistGo=min(DistGo,Dist0)
		
		DV(1) = sin(Eta)*cos(Seta)*DistGo
		DV(2) = sin(Eta)*sin(Seta)*DistGo
		DV(3) = cos(Eta)*DistGo
					
		VT(1) = V0(1) + DV(1)
		VT(2) = V0(2) + DV(2)
		VT(3) = V0(3) + DV(3)
	END IF
	
	!starting point proceed to new place
	V0 = VT



	RETURN


!-------> end of MarchOn()

	END



C======================================================================C
	SUBROUTINE	FindArrivalWall(DistW,Dist0,V0,V1,VW,DV,Seta,Eta,
     &			  IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NW,iray,ierr,iben,
     &			  XYZ,FanWei)
C======================================================================C

	USE MODULE_RADSXF, ONLY : WABCD,WNV,WElem,WArea
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID,GBlockNum
	USE module_radwork, ONLY : WK2


	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,NSEG = 5)
	
	REAL*8	XYZ(NDOF,MP),FanWei(2,NDOF)

	!Private
	REAL*8  V0(3),VW(3),V1(3),DV(3),ddd(3),vst(3)
     	INTEGER FaceID(5000),ista(3),iend(3),ithis(3)

	!ierr ---	0:		normal return
	!			1:		out of range
	!			2:		should arrived at wall, but failed to find it
	

	

	FDA = 0.1
	ierr = 0
	IDWall = 0
	iwk = 0
	iacrossflag = 0
	!FaceID = 0
	!WK2(1:NW) = 0

c2-------------------------------------------------------------
c	IF((jx.GT.NB).OR.(jy.GT.NB).OR.(jz.GT.NB)) THEN
c			!WRITE(IUT0,*) '  num is out of range,(jx,jy,jz) =',jx,jy,jz
c		ierr = 1
c		RETURN
c	ENDIF
c	IF((jx.LT.1).OR.(jy.LT.1).OR.(jz.LT.1)) THEN
c			!WRITE(IUT0,*) '  num is out of range,(jx,jy,jz) =',jx,jy,jz
c		ierr = 1
c		RETURN
c	END IF
c	IF(GBlockNum(jx,jy,jz).EQ.0) THEN
c			!WRITE(IUT0,*) '  origen is out of range'
c		ierr = 1
c		RETURN
c	END IF


c1-------------------------------------------------------------	


	dx = (FanWei(2,1)-FanWei(1,1))/NB
	dy = (FanWei(2,2)-FanWei(1,2))/NB
	dz = (FanWei(2,3)-FanWei(1,3))/NB
	ddd(1) = dx
	ddd(2) = dy
	ddd(3) = dz
	
	ista(1) = (V0(1)-FanWei(1,1))/dx+1
	ista(2) = (V0(2)-FanWei(1,2))/dy+1
	ista(3) = (V0(3)-FanWei(1,3))/dz+1
	
	iend(1) = (V1(1)-FanWei(1,1))/dx+1
	iend(2) = (V1(2)-FanWei(1,2))/dy+1
	iend(3) = (V1(3)-FanWei(1,3))/dz+1

	

	mmk=0
	IF((iend(1).GT.NB).OR.(iend(2).GT.NB).OR.
     &  (iend(3).GT.NB)) THEN
		iacrossflag = 1
		mmk=1
	ELSE IF((iend(1).LT.1).OR.(iend(2).LT.1).OR.
     &  (iend(3).LT.1)) THEN
		iacrossflag = 1
		mmk=2
	ELSE IF(GBlockNum(iend(1),iend(2),iend(3)).EQ.0) THEN
		iacrossflag = 1
	ENDIF

	
c2-------------------------------------------------------------	
	IDOF = 1
	vmax = abs(DV(1)/ddd(1))
	DO 30 I = 2, 3
		va = abs(DV(I)/ddd(I))
		IF(va.GT.vmax) THEN
		IDOF = I
		vmax = va
		END IF
30	CONTINUE
	
	
	MSEG = vmax*NSEG+1
	dseg = 1.0/MSEG
	
	
	i0 = 0
	j0 = 0
	k0 = 0
	ib = 0
	ib0 = ib
	
	ISEG = 0
	DO 1400 WHILE(ISEG.LT.MSEG+1.AND.IDWALL.LE.0)
	
		ISEG = ISEG+1  
		x1 = V0(1) + (ISEG-1)*DV(1)*dseg
		y1 = V0(2) + (ISEG-1)*DV(2)*dseg
		z1 = V0(3) + (ISEG-1)*DV(3)*dseg

		i = (x1-FanWei(1,1))/dx+1
		j = (y1-FanWei(1,2))/dy+1
		k = (z1-FanWei(1,3))/dz+1
		
		IF(i.GT.NB.OR.i.LT.1) GO TO 1400
		IF(j.GT.NB.OR.j.LT.1) GO TO 1400
		IF(k.GT.NB.OR.k.LT.1) GO TO 1400
		IF(i0.EQ.i.AND.j0.EQ.j.and.k0.EQ.k) GOTO 1400
		
	
		i0=i
		j0=j
		k0=k

		!left
		ib0 = ib
		DO 50 JJ =1,WBlockNum(i,j,k)
			II = WBlockIndex(i,j,k)+JJ
			IW = WBlockId(II)
			
			IF(WK2(IW).EQ.0) THEN
				ib = ib+1
				WK2(IW) = ib
				FaceID(ib) = IW 
			END IF
		
50		CONTINUE

	
			
		IA = ib0
		DO 690  WHILE(IA.LT.IB.AND.IDWALL.EQ.0)
			
			IA = IA + 1
			IW = FaceID(IA)

			MarkNotIt = 0
			MarkInEdgeLine = 0

					
			! wnv is outward
			CALL  COSVALUE(WNV(1,IW),WNV(2,IW),WNV(3,IW),
     &				   DV(1),DV(2),DV(3),vcos)

			!could not intersect this wall face
			IF(vcos.LE.0)	CYCLE		!GOTO 690
			
	
			!project node in face
		  CALL DistProj(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &		  WABCD(4,IW),V0(1),V0(2),V0(3),
     &		  WNV(1,IW),WNV(2,IW),WNV(3,IW),X2,Y2,Z2,DistN)
		CALL COSVALUE( X2-V0(1),Y2-V0(2),Z2-V0(3),
     &			   DV(1),DV(2),DV(3),value2)
		IF(value2.LE.0)	CYCLE			!GOTO 690
		DistW = DistN/vcos
		CALL IntersectNode(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			WABCD(4,IW),V0(1),V0(2),V0(3),DistW,Eta,Seta,
     &			VW(1),VW(2),VW(3),1,istatus)
		IF(istatus.GT.0) CYCLE			!GOTO 690
		aaa = 0.0		
		DO 610 K=1,WElem(1,IW)
		K1 = K+1
		IF(K1.GT.WElem(1,IW)) K1=K1-WElem(1,IW)
		IP0 = WElem(K+1,IW)
		IP1 = WElem(K1+1,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &		XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &		VW(1),VW(2),VW(3),vx,vy,vz,bbb,1.d0,0)
		aaa = aaa + bbb
610		CONTINUE
		darea = abs(aaa-WArea(IW))/WArea(IW)
		IF(darea.LE.1.0d-6) THEN
		IDWall = IW
		ierr = 0
                ELSE IF(darea.LE.1.0d-2) THEN
		IDWall = -IW
		END IF
690	CONTINUE	
1400	CONTINUE
	numall = ib
	IDWALL = ABS(IDWALL)
	IF(IDWALL.GT.0) THEN	!found a possible intersecting surface patch
	IF(DistW.GT.Dist0) THEN
	IDWall = -IDWALL			!too far away, do not intersect
	ierr = 0
        END IF
	DO I = 1,NUMALL
	IW = FaceID(I)
	WK2(IW) = 0
	END DO
	RETURN
	END IF
c2-------------------------------------------------------------
	IF(IDWall.EQ.0.AND.iacrossflag.EQ.1) THEN
	!enlarge scanning area
	iben = iben+1

		imin = min(ista(1),iend(1))-1
		imax = max(ista(1),iend(1))+1
		jmin = min(ista(2),iend(2))-1
		jmax = max(ista(2),iend(2))+1
		kmin = min(ista(3),iend(3))-1
		kmax = max(ista(3),iend(3))+1

	    imin = max(1,imin)
		imax = min(NB,imax)
		jmin = max(1,jmin)
		jmax = min(NB,jmax)
		kmin = max(1,kmin)
		kmax = min(NB,kmax)

		ib = numall
		DO 2200 ix = imin, imax
		  DO 2200 iy = jmin, jmax
			DO 2200 iz = kmin, kmax
				DO 2190 I=1,WBlockNum(ix,iy,iz)
					ID = WBlockIndex(ix,iy,iz)+I
					IW = WBlockId(ID)
					IF(WK2(IW).EQ.0) THEN
						ib = ib+1
						WK2(IW) = ib
						FaceID(ib) = IW 
					END IF
2190				CONTINUE
2200		CONTINUE

		
		!search again
		IA = NUMALL
		DO 2690  WHILE(IA.LT.IB.AND.IDWALL.EQ.0)
			
			IA = IA+1
			IW = FaceID(IA)

			MarkNotIt = 0
			MarkInEdgeLine = 0
			
			CALL  COSVALUE(WNV(1,IW),WNV(2,IW),WNV(3,IW),
     &				   DV(1),DV(2),DV(3),vcos)

			!could not intersect this wall face
			IF(vcos.LE.0) CYCLE
					
		CALL DistProj(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &		  WABCD(4,IW),V0(1),V0(2),V0(3),
     &		  WNV(1,IW),WNV(2,IW),WNV(3,IW),X2,Y2,Z2,DistN)

		
			CALL COSVALUE( X2-V0(1),Y2-V0(2),Z2-V0(3),
     &					   DV(1),DV(2),DV(3),value2)
			IF(value2.LE.0) CYCLE
			
			DistW = DistN/vcos
						
		CALL IntersectNode(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			WABCD(4,IW),V0(1),V0(2),V0(3),DistW,Eta,Seta,
     &			VW(1),VW(2),VW(3),1,istatus)
		IF(istatus.GT.0) CYCLE		!not it
		aaa = 0.0
		DO 2610 K=1,WElem(1,IW)
		K1 = K+1
		IF(K1.GT.WElem(1,IW)) K1=K1-WElem(1,IW)
		IP0 = WElem(K+1,IW)
		IP1 = WElem(K1+1,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &		XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &		VW(1),VW(2),VW(3),vx,vy,vz,bbb,1.d0,0)
		aaa = aaa + bbb
2610		CONTINUE
		darea = abs(aaa-WArea(IW))/WArea(IW)
		IF(darea.LE.1.0d-6) THEN
		IDWall = IW
		ierr = 0
	        ELSE IF(darea.LE.1.0d-2) THEN
		IDWall = -IW
		END IF
2690		CONTINUE
		IF(IDWall.EQ.0) THEN
		ierr = 2
		END IF
	END IF

	IDWALL= ABS(IDWALL)
	IF(IDWALL.GT.0.AND.DistW.GT.Dist0) THEN
		IDWall = -IDWALL !too far away, do not intersect
		ierr = 0
	END IF

	!reset the work array
	DO I = 1,MAX(IB,NUMALL)
		IW = FaceID(I)
		WK2(IW) = 0
	END DO

	RETURN
!-------> end of FindArrivalWall()
	END

C======================================================================C
	SUBROUTINE FindArrivalWall_old
     &            (DistW,Dist0,V0,V1,VW,DV,Seta,Eta,
     &		  IDWall,N2D,N3D,NDOF,NB,NPBW,MP,NW,iray,ierr,iben,
     &			  XYZ,FanWei)
C======================================================================C

	USE MODULE_RADSXF, ONLY : WABCD,WNV,WElem,WArea
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID,GBlockNum
	
	USE module_radwork, ONLY : WK2


	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,NSEG = 4)

	REAL*8	XYZ(NDOF,MP),FanWei(2,NDOF)

	REAL*8  V0(3),VW(3),V1(3),DV(3),ddd(3),vst(3)
     	INTEGER FaceID(5000),ista(3),iend(3),ithis(3)

	!ierr --0:		normal return
	!	1:		out of range
	!	2:		should arrived at wall, but failed to find it
	



	FDA = 0.1
	ierr = 0
	IDWall = 0
	iwk = 0
	iacrossflag = 0
	!FaceID = 0
	!WK2(1:NW) = 0

c2-------------------------------------------------------------
c	IF((jx.GT.NB).OR.(jy.GT.NB).OR.(jz.GT.NB)) THEN
c			!WRITE(IUT0,*) '  num is out of range,(jx,jy,jz) =',jx,jy,jz
c		ierr = 1
c		RETURN
c	ENDIF
c	IF((jx.LT.1).OR.(jy.LT.1).OR.(jz.LT.1)) THEN
c			!WRITE(IUT0,*) '  num is out of range,(jx,jy,jz) =',jx,jy,jz
c		ierr = 1
c		RETURN
c	END IF
c	IF(GBlockNum(jx,jy,jz).EQ.0) THEN
c			!WRITE(IUT0,*) '  origen is out of range'
c		ierr = 1
c		RETURN
c	END IF


c1-------------------------------------------------------------	

	dx = (FanWei(2,1)-FanWei(1,1))/NB
	dy = (FanWei(2,2)-FanWei(1,2))/NB
	dz = (FanWei(2,3)-FanWei(1,3))/NB
	ddd(1) = dx
	ddd(2) = dy
	ddd(3) = dz
	
	ista(1) = (V0(1)-FanWei(1,1))/dx+1
	ista(2) = (V0(2)-FanWei(1,2))/dy+1
	ista(3) = (V0(3)-FanWei(1,3))/dz+1
	
	iend(1) = (V1(1)-FanWei(1,1))/dx+1
	iend(2) = (V1(2)-FanWei(1,2))/dy+1
	iend(3) = (V1(3)-FanWei(1,3))/dz+1

	
! --- 7
	IF((iend(1).GT.NB).OR.(iend(2).GT.NB).OR.(iend(3).GT.NB)) THEN
	iacrossflag = 1
	ELSE IF((iend(1).LT.1).OR.(iend(2).LT.1).
     &            OR.(iend(3).LT.1)) THEN
	iacrossflag = 1
	ELSE IF(GBlockNum(iend(1),iend(2),iend(3)).EQ.0) THEN
	iacrossflag = 1
	ENDIF
c2-------------------------------------------------------------	
	IDOF = 1
	vmax = abs(DV(1)/ddd(1))
	DO 30 I = 2, 3
	va = abs(DV(I)/ddd(I))
	IF(va.GT.vmax) THEN
	IDOF = I
	vmax = va
	END IF
30	CONTINUE
!	
	MSEG = vmax*NSEG+1
	dseg = 1.0/MSEG
!	
	i0 = 0
	j0 = 0
	k0 = 0
	ib = 0
	ib0 = ib
	
	ISEG = 0
	DO 1400 WHILE(ISEG.LT.MSEG+1.AND.IDWALL.EQ.0)
!
	ISEG = ISEG+1  
	x1 = V0(1) + (ISEG-1)*DV(1)*dseg
	y1 = V0(2) + (ISEG-1)*DV(2)*dseg
	z1 = V0(3) + (ISEG-1)*DV(3)*dseg

	i = (x1-FanWei(1,1))/dx+1
	j = (y1-FanWei(1,2))/dy+1
	k = (z1-FanWei(1,3))/dz+1
		
	IF(i.GT.NB.OR.i.LT.1) GO TO 1400
	IF(j.GT.NB.OR.j.LT.1) GO TO 1400
	IF(k.GT.NB.OR.k.LT.1) GO TO 1400
	IF(i0.EQ.i.AND.j0.EQ.j.and.k0.EQ.k) GOTO 1400
!left
	ib0 = ib
	DO 50 JJ =1,WBlockNum(i,j,k)
	II = WBlockIndex(i,j,k)+JJ
	IW = WBlockId(II)
	IF(WK2(IW).EQ.0) THEN
	ib = ib+1
	WK2(IW) = ib
	FaceID(ib) = IW 
	END IF
 50	CONTINUE
	IA = ib0
	DO 690  WHILE(IA.LT.IB.AND.IDWALL.EQ.0)
	IA = IA + 1
	IW = FaceID(IA)
	MarkNotIt = 0
	MarkInEdgeLine = 0
! wnv is outward
	CALL  COSVALUE(WNV(1,IW),WNV(2,IW),WNV(3,IW),
     &			   DV(1),DV(2),DV(3),vcos)
			!could not intersect this wall face
			IF(vcos.LE.0)	CYCLE		!GOTO 690
					

			
			!project node in face
		CALL DistProj(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			  WABCD(4,IW),V0(1),V0(2),V0(3),
     &			  WNV(1,IW),WNV(2,IW),WNV(3,IW),X2,Y2,Z2,DistN)

			
			!
		CALL COSVALUE( X2-V0(1),Y2-V0(2),Z2-V0(3),
     &				   DV(1),DV(2),DV(3),value2)
		IF(value2.LE.0)	CYCLE			!GOTO 690
		DistW = DistN/vcos
		CALL IntersectNode(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &		WABCD(4,IW),V0(1),V0(2),V0(3),DistW,Eta,Seta,
     &		VW(1),VW(2),VW(3),1,istatus)
		IF(istatus.GT.0) CYCLE				!GOTO 690
		MarkNotIt = 0
		MarkInEdgeLine = 0
		K = 0
		vx1 = 0.0
		vy1 = 0.0
		vz1 = 0.0
		IP0 = WElem(2,IW)
		IP1 = WElem(3,IW)
		IP2 = WElem(4,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &			XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &			XYZ(1,IP2),XYZ(2,IP2),XYZ(3,IP2),
     &			vx,vy,vz,aaa,1.d0,0)
		DO 610 WHILE((K.LT.WElem(1,IW)).AND.(MarkNotIt.EQ.0))
		K =  K + 1
		K1 = MOD(K,WElem(1,IW))+1
		IP0 = WElem(K+1,IW)
		IP1 = WElem(K1+1,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &			XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &			VW(1),VW(2),VW(3),vx,vy,vz,bbb,aaa,1)
		vmodule = DSQRT(vx*vx + vy*vy +vz*vz)
		value = 0
		IF(vmodule < 1.0E-3)	THEN
!the angle is less than 1 deg,
!the point is at the edge-line of the patch
		MarkInEdgeLine = 1
		ELSE
		CALL COSVALUE(vx1,vy1,vz1,vx,vy,vz,value)
		IF(value.LT.-1.0E-4) MarkNotIt = 1
!the vectors are opposite,
!this means the point is not bounded by 
!the element/patch
		vx1 = vx
		vy1 = vy
		vz1 = vz
		END IF
610	CONTINUE
	IF((MarkNotIt.EQ.0)) THEN	!found the intersecting wall patch
	IDWall = IW
	ierr = 0
        ENDIF
	
690	CONTINUE	

1400	CONTINUE

	numall = ib

	IF(IDWALL.GT.0) THEN	!found a possible intersecting surface patch
	
		IF(DistW.GT.Dist0) THEN
			IDWall = -IDWALL			!too far away, do not intersect
			ierr = 0
		END IF
	
		DO I = 1,NUMALL
			IW = FaceID(I)
			WK2(IW) = 0
		END DO
			
		RETURN

	END IF




c2-------------------------------------------------------------
	
	IF(IDWall.EQ.0.AND.iacrossflag.EQ.1) THEN
		
		!enlarge scanning area
		iben = iben+1

		imin = min(ista(1),iend(1))-1
		imax = max(ista(1),iend(1))+1
		jmin = min(ista(2),iend(2))-1
		jmax = max(ista(2),iend(2))+1
		kmin = min(ista(3),iend(3))-1
		kmax = max(ista(3),iend(3))+1

	    imin = max(1,imin)
		imax = min(NB,imax)
		jmin = max(1,jmin)
		jmax = min(NB,jmax)
		kmin = max(1,kmin)
		kmax = min(NB,kmax)

		ib = numall
		DO 2200 ix = imin, imax
		  DO 2200 iy = jmin, jmax
			DO 2200 iz = kmin, kmax
				DO 2190 I=1,WBlockNum(ix,iy,iz)
					ID = WBlockIndex(ix,iy,iz)+I
					IW = WBlockId(ID)
					IF(WK2(IW).EQ.0) THEN
						ib = ib+1
						WK2(IW) = ib
						FaceID(ib) = IW 
					END IF
2190				CONTINUE
2200		CONTINUE

		
		!search again
		IA = NUMALL
		DO 2690  WHILE(IA.LT.IB.AND.IDWALL.EQ.0)
			
			IA = IA+1
			IW = FaceID(IA)

			MarkNotIt = 0
			MarkInEdgeLine = 0
			
			CALL  COSVALUE(WNV(1,IW),WNV(2,IW),WNV(3,IW),
     &				   DV(1),DV(2),DV(3),vcos)


			!could not intersect this wall face
			IF(vcos.LE.0) GOTO 2690
					
		CALL DistProj(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			  WABCD(4,IW),V0(1),V0(2),V0(3),
     &		  WNV(1,IW),WNV(2,IW),WNV(3,IW),X2,Y2,Z2,DistN)

			

		CALL COSVALUE( X2-V0(1),Y2-V0(2),Z2-V0(3),
     &				   DV(1),DV(2),DV(3),value2)
		IF(value2.LE.0) GOTO 2690
			
		DistW = DistN/vcos
						
		CALL IntersectNode(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			WABCD(4,IW),V0(1),V0(2),V0(3),DistW,Eta,Seta,
     &			VW(1),VW(2),VW(3),1,istatus)
		IF(istatus.GT.0) GOTO 2690			!not it
		MarkNotIt = 0
		MarkInEdgeLine = 0
		K = 0
		vx1 = 0.0
		vy1 = 0.0
		vz1 = 0.0
		IP0 = WElem(2,IW)
		IP1 = WElem(3,IW)
		IP2 = WElem(4,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &			XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &			XYZ(1,IP2),XYZ(2,IP2),XYZ(3,IP2),
     &			vx,vy,vz,aaa,1.d0,0)
		DO 2610 WHILE((K.LT.WElem(1,IW)).AND.(MarkNotIt.EQ.0))
		K =  K + 1
		K1 = MOD(K,WElem(1,IW))+1
		IP0 = WElem(K+1,IW)
		IP1 = WElem(K1+1,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &				XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &				VW(1),VW(2),VW(3),vx,vy,vz,bbb,aaa,1)
		vmodule = DSQRT(vx*vx + vy*vy +vz*vz)
		value = 0
		IF(vmodule < 1.0E-3)	THEN
!the angle is less than 1 deg,
!the point is at the edge-line of the patch
		MarkInEdgeLine = 1
		ELSE
		CALL COSVALUE(vx1,vy1,vz1,vx,vy,vz,value)
		IF(value.LT.-1.0E-4) MarkNotIt = 1
!the vectors are opposite,
!this means the point is not bounded by 
!the element/patch
		vx1 = vx
		vy1 = vy
		vz1 = vz
		END IF
2610		CONTINUE
		IF((MarkNotIt.EQ.0)) THEN	
!found the intersecting wall patch
		IDWall = IW
		ierr = 0
	        ENDIF
2690		CONTINUE
	        IF(IDWall.EQ.0) THEN
		ierr = 2
		END IF
	END IF
!reset the work array
	DO I = 1,MAX(IB,NUMALL)
	IW = FaceID(I)
	WK2(IW) = 0
	END DO

	RETURN

!-------> end of FindArrivalWall()

	END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE FindArrivalWall2(DistW,Dist0,V0,V1,VW,DV,Seta,Eta,
     &		  MeshRange,IDWall,WABCD,WNV,XYZ,WElem,WBlockNum,
     &		  WBlockIndex,WBlockID,GBlockNum,N2D,N3D,NDOF,NB,NPBW,
     &		  MP,NW,WK,iray,ierr)
!
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793)

	REAL*8  XYZ(NDOF,MP),V0(3),VW(3),V1(3),DV(3),ddd(3),vst(3)
	REAL*8  WABCD(4,NW),WNV(NDOF,NW),MeshRange(2,NDOF)
	INTEGER WElem(N2D+1,NW),WK(NW)
	INTEGER WBlockNum(NB,NB,NB),WBlockIndex(NB,NB,NB),
     &		WBlockID(NPBW),GBlockNum(NB,NB,NB)
     	INTEGER FaceID(500)
     			

!ierr ---
!	0:		normal return
!	1:		out of range
!	2:		should arrived at wall, but failed to find it
	

!---------------------------------------------------------------------------------------------!
!NOTE:
!THIS METHOD CAN ONLY APPLYED TO MESH WITHOUT INNER SURFACE OR MULTI-LAYER WALLS
!---------------------------------------------------------------------------------------------!



	FDA = 0.1
	FaceID = 0
	WK(1:NW) = 0
	iwk = 0
	ierr = 0
	IDWall = -1

	
c1-------------------------------------------------------------
	
	dx = (MeshRange(2,1)-MeshRange(1,1))/NB
	dy = (MeshRange(2,2)-MeshRange(1,2))/NB
	dz = (MeshRange(2,3)-MeshRange(1,3))/NB
	kx = 1
	ky = 1
	kz = 1

	jx = (V0(1)-MeshRange(1,1))/dx+1
	jy = (V0(2)-MeshRange(1,2))/dy+1
	jz = (V0(3)-MeshRange(1,3))/dz+1



c3-------------------------------------------------------------	

	iacrossflag = 0

	kx = (V1(1)-MeshRange(1,1))/dx+1
	ky = (V1(2)-MeshRange(1,2))/dy+1
	kz = (V1(3)-MeshRange(1,3))/dz+1


	IF((kx.GT.NB).OR.(ky.GT.NB).OR.(kz.GT.NB)) THEN
		iacrossflag = 1
	ELSE IF((kx.LT.1).OR.(ky.LT.1).OR.(kz.LT.1)) THEN
		iacrossflag = 1
	ELSE IF(GBlockNum(kx,ky,kz).EQ.0) THEN
		iacrossflag = 1
	ENDIF



!c4---------------------------------------------------------------
100	imin = min(jx,kx)
	imax = max(jx,kx)
	jmin = min(jy,ky)
	jmax = max(jy,ky)
	kmin = min(jz,kz)
	kmax = max(jz,kz)
	
	imax = min(NB,imax)
	imin = max(1,imin)
	jmax = min(NB,jmax)
	jmin = max(1,jmin)
	kmax = min(NB,kmax)
	kmin = max(1,kmin)

	ib = 0
	DO 200 ix = imin, imax
	  DO 200 iy = jmin, jmax
	    DO 200 iz = kmin, kmax
			
			DO 190 I=1,WBlockNum(ix,iy,iz)
			ID = WBlockIndex(ix,iy,iz)+I
			IW = WBlockId(ID)
			IF(WK(IW).EQ.0) THEN
				ib = ib+1
				WK(IW) = ib
				FaceID(ib) = IW 
			END IF
190			CONTINUE
200	CONTINUE

	numall = ib
	write(*,*) numall

	IF(numall.EQ.0) THEN		!find no wall
		IF(iacrossflag.EQ.1) THEN
			!should have, but failed to get
			IDWall = -1
			ierr = 2
			RETURN
		ELSE
			!arrive at inner place
			IDWall = -1
			ierr = 0
			RETURN
		END IF
	END IF


!c5---------------------------------------------------------------



	IDWall = -1
	
	DO 690 IB = NUMALL,1,-1
  
			IW = FaceID(IB)

			MarkNotIt = 0
			MarkInEdgeLine = 0
			
			CALL  COSVALUE(WNV(1,IW),WNV(2,IW),WNV(3,IW),
     &					   DV(1),DV(2),DV(3),vcos)


			!could not intersect this wall face
			IF(vcos.LE.0) GOTO 690
					
		CALL DistProj(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			  WABCD(4,IW),V0(1),V0(2),V0(3),
     &			  WNV(1,IW),WNV(2,IW),WNV(3,IW),X1,Y1,Z1,DistN)

			

		CALL COSVALUE( X1-V0(1),Y1-V0(2),Z1-V0(3),
     &				   DV(1),DV(2),DV(3),value2)
		IF(value2.LE.0) GOTO 690
			
		DistW = DistN/vcos
						
		CALL IntersectNode(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			WABCD(4,IW),V0(1),V0(2),V0(3),DistW,Eta,Seta,
     &			VW(1),VW(2),VW(3),1,istatus)
		IF(istatus.GT.0) GOTO 690
!not it
		
			MarkNotIt = 0
			MarkInEdgeLine = 0
			K = 0
			vx1 = 0.0
			vy1 = 0.0
			vz1 = 0.0
			

			IP0 = WElem(2,IW)
			IP1 = WElem(3,IW)
			IP2 = WElem(4,IW)
			
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &			XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &			XYZ(1,IP2),XYZ(2,IP2),XYZ(3,IP2),
     &			vx,vy,vz,aaa,1.d0,0)

		DO 610 WHILE((K.LT.WElem(1,IW)).AND.(MarkNotIt.EQ.0))
			K =  K + 1
			K1 = MOD(K,WElem(1,IW))+1
			IP0 = WElem(K+1,IW)
			IP1 = WElem(K1+1,IW)
	
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &				XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &			VW(1),VW(2),VW(3),vx,vy,vz,bbb,aaa,1)
		vmodule = DSQRT(vx*vx + vy*vy +vz*vz)
		value = 0
		IF(vmodule < 1.0E-3)	THEN
			!the angle is less than 1 deg,
			!the point is at the edge-line of the patch
		MarkInEdgeLine = 1
		ELSE
		CALL COSVALUE(vx1,vy1,vz1,vx,vy,vz,value)
!value = +-1
		IF(value.LT.-1.0E-4) MarkNotIt = 1
!the vectors are opposite,
!this means the point is not bounded by 
!the element/patch
		vx1 = vx
		vy1 = vy
		vz1 = vz
		END IF
 610		CONTINUE
		IF((MarkNotIt.EQ.0)) THEN	
!found the intersecting wall patch
		IDWall = IW
		ierr = 0
		GOTO 1000
		ENDIF
 690	CONTINUE
	
 1000	IF(DistW.GT.Dist0) THEN
	IDWall = -1
	ierr = 0
	RETURN
	END IF
!
	IF(IDWall.EQ.-1.AND.iacrossflag.EQ.1) THEN
	ierr = 2
        ELSE
	ierr = 0
	END IF

2000	RETURN

!-------> end of FindArrivalWall()

	END
!
	SUBROUTINE	
     & GetIntersectWall(DistW,Dist0,V0,V1,VW,DV,Seta,Eta,
     & MeshRange,IDWall,WABCD,WNV,XYZ,WElem,WBlockNum,
     & WBlockIndex,WBlockID,GBlockNum,N2D,N3D,NDOF,NB,NPBW,
     & MP,NW,WK,istat)

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793)

	REAL*8  XYZ(NDOF,MP),V0(3),VW(3),V1(3),DV(3),ddd(3),vst(3)
	REAL*8  WABCD(4,NW),WNV(NDOF,NW),MeshRange(2,NDOF)
	INTEGER WElem(N2D+1,NW),WK(NW)
	INTEGER WBlockNum(NB,NB,NB),WBlockIndex(NB,NB,NB),
     &		WBlockID(NPBW),GBlockNum(NB,NB,NB)
     	INTEGER FaceID(500)
     			

!ierr ---
!	0:		normal return
!	1:		out of range
!	2:		should arrived at wall, but failed to find it
	

!---------------------------------------------------------------------------------------------!
!NOTE:
!		THIS METHOD CAN ONLY APPLYED TO MESH WITHOUT INNER SURFACE OR MULTI-LAYER WALLS
!---------------------------------------------------------------------------------------------!

	FDA = 0.1
	ierr = 0
	IDWall = -1
	iwk = 0

	istat = 0
	
c1-------------------------------------------------------------
	
	dx = (MeshRange(2,1)-MeshRange(1,1))/NB
	dy = (MeshRange(2,2)-MeshRange(1,2))/NB
	dz = (MeshRange(2,3)-MeshRange(1,3))/NB

	kx = (V1(1)-MeshRange(1,1))/dx+1
	ky = (V1(2)-MeshRange(1,2))/dy+1
	kz = (V1(3)-MeshRange(1,3))/dz+1

	

	IF((kx.GT.NB).OR.(ky.GT.NB).OR.(kz.GT.NB)) THEN
		return
	ELSE IF((kx.LT.1).OR.(ky.LT.1).OR.(kz.LT.1)) THEN
		return
	ELSE IF(GBlockNum(kx,ky,kz).EQ.0) THEN
		return
	ENDIF


	DO 190 I=1,WBlockNum(kx,ky,kz)
		ID = WBlockIndex(kx,ky,kz)+I
		IW = WBlockId(ID)
		IF(WK(IW).EQ.0) THEN
			ib = ib+1
			WK(IW) = ib
			FaceID(ib) = IW 
		END IF
190	CONTINUE

	numall = ib


c2---------------------------------------------------------------------

	IDWall = -1
	
	DO 690 IB = 1,NUMALL
  
			IW = FaceID(IB)

			MarkNotIt = 0
			MarkInEdgeLine = 0
			
			CALL  COSVALUE(WNV(1,IW),WNV(2,IW),WNV(3,IW),
     &					   DV(1),DV(2),DV(3),vcos)


			!could not intersect this wall face
			IF(vcos.LE.0) GOTO 690
					
		CALL DistProj(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			  WABCD(4,IW),V0(1),V0(2),V0(3),
     &			  WNV(1,IW),WNV(2,IW),WNV(3,IW),X1,Y1,Z1,DistN)

			
		CALL COSVALUE( X1-V0(1),Y1-V0(2),Z1-V0(3),
     &				   DV(1),DV(2),DV(3),value2)
		IF(value2.LE.0) GOTO 690
			
		DistW = DistN/vcos
						
		CALL IntersectNode(WABCD(1,IW),WABCD(2,IW),WABCD(3,IW),
     &			WABCD(4,IW),V0(1),V0(2),V0(3),DistW,Eta,Seta,
     &			VW(1),VW(2),VW(3),1,istatus)
		IF(istatus.GT.0) GOTO 690	!not it
		MarkNotIt = 0
		MarkInEdgeLine = 0
		K = 0
		vx1 = 0.0
		vy1 = 0.0
		vz1 = 0.0
		IP0 = WElem(2,IW)
			IP1 = WElem(3,IW)
			IP2 = WElem(4,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &				XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &				XYZ(1,IP2),XYZ(2,IP2),XYZ(3,IP2),
     &				vx,vy,vz,aaa,1.d0,0)
		DO 610 WHILE((K.LT.WElem(1,IW)).AND.(MarkNotIt.EQ.0))
			K =  K + 1
			K1 = MOD(K,WElem(1,IW))+1
			IP0 = WElem(K+1,IW)
			IP1 = WElem(K1+1,IW)
		CALL AreaVector2(XYZ(1,IP0),XYZ(2,IP0),XYZ(3,IP0),
     &			XYZ(1,IP1),XYZ(2,IP1),XYZ(3,IP1),
     &			VW(1),VW(2),VW(3),vx,vy,vz,bbb,aaa,1)
		vmodule = DSQRT(vx*vx + vy*vy +vz*vz)
		value = 0
		IF(vmodule < 1.0E-3)	THEN
!the angle is less than 1 deg,
!the point is at the edge-line of the patch
		MarkInEdgeLine = 1
		ELSE
		CALL COSVALUE(vx1,vy1,vz1,vx,vy,vz,value)
		IF(value.LT.-1.0E-4) MarkNotIt = 1
!the vectors are opposite,
!this means the point is not bounded by 
!the element/patch
		vx1 = vx
		vy1 = vy
		vz1 = vz
		END IF
 610		CONTINUE
		IF((MarkNotIt.EQ.0)) THEN	
!found the intersecting wall patch
		IDWall = IW
		istat = 0
		RETURN
		ENDIF
690	CONTINUE

!-------> end of FindArrivalWall()

	END
C======================================================================C

	!Period <= 6.0E5
	SUBROUTINE  RANDOM1(SEED,RAN)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 RAN,SEED

	!seed0 = 5249347
	IF(SEED.EQ.0) SEED = 5249347.0
100	SEED = DMOD(23.0*SEED,10000001.0d0)
	RAN = SEED*1.0E-7
	
	END



	!Period <= 2^31
	SUBROUTINE  RANDOM2(SEED,RAN)
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 RAN,SEED

	!seed0 = 5249347
	IF(SEED.EQ.0) SEED = 5249347.0
100	SEED = DMOD(131075.0*SEED,2147483649.0d0)
	RAN = SEED/2147483649.0
	
	END



	!Period <= 2^31
	SUBROUTINE  RANDOM3(SEED,RAN)
	IMPLICIT REAL*8(A-H,O-Z)	
	REAL*8 RAN,SEED

	!seed0 = 9
	IF(SEED.EQ.0) SEED = 9.0
100	SEED = DMOD(69069.0*SEED,2147483648.0d0)
	RAN = SEED/2147483648.0

	END



	!Period > 2^31
	SUBROUTINE  RANDOM4(SEED1,SEED2,RAN)
	IMPLICIT REAL*8(A-H,O-Z)	
	REAL*8 RAND,RAN,SEED1,SEED2

	!seed1 = 9, seed2 = 11
	IF(SEED1.EQ.0) SEED1 = 9.0
	IF(SEED2.EQ.0) SEED2 = 11.0
100	RAND = SEED1*65539.0 + SEED2*65539.0
	RAND = DMOD(RAND,2147483648.0d0)
	RAN = RAND/2147483648.0

	SEED1 = SEED2
	SEED2 = RAND

	END


C======================================================================C
	SUBROUTINE  ScatterAngleIso(Seta,Eta,Seta1,Eta1)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8	NV(3),TV1(3),TV2(3)

!this is a subroutine to give the emission angle of the scattered
!light bundle, 
!seta, eta	:	incidence angle (circum,cone)
!seta1, eta1	:	scattering angle (circum,cone)


	!to make a coordinate system for the incidence direction
	NV(1) = sin(Eta)*cos(Seta)
	NV(2) = sin(Eta)*sin(Seta)
	NV(3) = cos(Eta)
	CALL BuildZuoBiao(NV,TV1,TV2,3)
	
	!get the angles (eta1,seta1) in origenal coodinates
	CALL ZuobiaoRot(Seta1,Eta1,NV(1),NV(2),NV(3),
     &		   TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))	
	

	Seta = Seta1
	Eta = Eta1

	RETURN

	END

C======================================================================C
       SUBROUTINE  ScatterAngleLAS
     & (Seta0,Eta0,RanSeta,RanEta,ptype,Beta)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.1415926)
	REAL*8	NV(3),TV1(3),TV2(3)
	INTEGER ptype

!this is a subroutine to give the emission angle of the scattered
!light bundle, 
!PTYPE:		1,Linear-anisotropic scatter / isotropic scatter, 
!		2,large-diffuse particle scatter
!		3,Rayleigh scatter





!=========================================================================!
!	PTYPE = 1	:	Isotropic / Linear Anisoptropic scattering, (MIE theory)
!=========================================================================!

!the phase function is defined by a linear Legendre function, that,
!
!		fai(eta)  = 1 + beta*cos(eta)
!
!(if beta = 0, it is isotropic scattering)
!
!and, RanEta = integrate_[0.5*fai(eta)*sin(eta)]	| (0--eta)
!here, the solid angle = sin(eta)*deta*dseta
!
!RanEta = {(1+0.25beta) - [cos(eta) + 0.25*beta*cos(2eta)]}/2.0
!--> get Eta
	
!----------------------- ATTENTION!! --------------------------
!For this kind of scatter phase function, there should be
!		-1.0 <= beta <= 1.0
!--------------------------------------------------------------

!seta0, eta0,		incident angle (circum,cone)
!seta, eta,			scattering bundle angle (circum,cone), ==> seta0,eta0


	IF(PTYPE.EQ.1) THEN

	  seta = 2.0*PAI*RanSeta
	
	  IF(beta.EQ.0) THEN
		!Isotropic scattering
		Eta = ACOS(1.0-2.0*RanEta)
		
		seta0 = seta
		eta0 = eta
		RETURN

	  ELSE
!unisotropic scattering, get eta's value by dichotomy scheme
!the phase function increases monotonously with angle eta
		
		eta1 = 0.0
		eta2 = PAI

		eta = (eta1+eta2)/2
	v3 = ((1+0.25*beta)-(cos(eta)+0.25*beta*cos(2.0*eta)))/2.0

		err = abs(v3-RanEta)
	
		itr = 0
		DO 100 WHILE(err.GT.1.0E-4)
		  itr = itr + 1
		  IF(v3.LT.RanEta) THEN		!too small
			eta1 = eta
			eta = (eta + eta2)/2.0
		  ELSE
			eta2 = eta
			eta = (eta + eta1)/2.0
		  END IF
	  v3 = ((1+0.25*beta)-(cos(eta)+0.25*beta*cos(2.0*eta)))/2.0
		  err = abs(v3-RanEta)
		  
		  !write(*,*) itr,err

100		CONTINUE

	  END IF
	  
	  GOTO 1000
	END IF
!=========================================================================!
!	PTYPE = 2 :	Large-Diffuse Particle Scattering
!=========================================================================!
!this is a subroutine to give the emission angle of the scattered
!light bundle, 
!the phase function is defined at the textbook of Taniguchi et al. (1994, Page.18, 110)
!
!		fai(eta)  = 8/(3*PAI)*(sin(eta) - eta*cos(eta))
!
!
!and, RanEta = integrate_[0.5*fai(eta)*sin(eta)]	| (0--eta)
!
!RanEta = 2/(3Pai)*[eta - 3/4*sin(2eta) + eta*cos(2eta)/2]		---> get eta
!seta0, eta0,		incident angle (circum,cone)
!seta, eta,			scattering bundle angle (circum,cone), ==> seta0,eta0

	IF(PTYPE.EQ.2) THEN
	seta = 2.0*PAI*RanSeta
!unisotropic scattering, get eta's value by dichotomy scheme
!the phase function increases monotonously with angle eta
	eta1 = 0.0
	eta2 = PAI
	eta = (eta1+eta2)/2
	v3 = (eta-0.75*sin(2.*eta)+0.5*eta*cos(2.0*eta))*2.0/3.0/Pai
	err = abs(v3-RanEta)
	DO 200 WHILE(err.GT.1.0E-4)
	  IF(v3.LT.RanEta) THEN		!too small
	eta1 = eta
	eta = (eta + eta2)/2.0
	  ELSE
	eta2 = eta
	eta = (eta + eta1)/2.0
	  END IF
        v3 = (eta-0.75*sin(2.*eta)+0.5*eta*cos(2.*eta))
     &		   *2./3./Pai
	err = abs(v3-RanEta)
		
200	CONTINUE

	GOTO 1000
	END IF
	
!=========================================================================!
!	PTYPE = 3	:  Rayleigh Scattering for Small Particles
!=========================================================================!
!this is a subroutine to give the emission angle of the scattered
!light bundle, 
!the phase function is defined at the textbook of Taniguchi et al. (1994, Page.18, 110)
!
!		fai(eta)  = 3/4*(1 + cos(eta)**2)
!
!
!and, RanEta = integrate_[0.5*fai(eta)*sin(eta)]	| (0--eta)
!
!RanEta = {1.0-3/16[(5*cos(eta) + 1/3*cos(3eta)]}/2.0	---> get eta
!seta0, eta0,		incident angle (circum,cone)
!seta, eta,	scattering bundle angle (circum,cone), ==> seta0,eta0
!	
	IF(PTYPE.EQ.3) THEN
	seta = 2.0*PAI*RanSeta
!unisotropic scattering, get eta's value by dichotomy scheme
!the phase function increases monotonously with angle eta
	eta1 = 0.0
	eta2 = PAI

	eta = (eta1+eta2)/2
	v3 = 0.5-0.09375*(5.0*cos(eta)+cos(3*eta)/3.0)
			
	err = abs(v3-RanEta)
	DO 300 WHILE(err.GT.1.0E-4)
	  IF(v3.LT.RanEta) THEN		!too small
		eta1 = eta
		eta = (eta + eta2)/2.0
	  ELSE
		eta2 = eta
		eta = (eta + eta1)/2.0
	  END IF
	  v3 =  0.5-0.09375*(5.0*cos(eta)+cos(3*eta)/3.0)
	  err = abs(v3-RanEta)
300		CONTINUE
	END IF
!to make a coordinate system for the incidence direction
 1000	NV(1) = sin(Eta0)*cos(Seta0)
	NV(2) = sin(Eta0)*sin(Seta0)
	NV(3) = cos(Eta0)
	CALL BuildZuoBiao(NV,TV1,TV2,3)
!get the angles (eta1,seta1) in origenal coodinates
	CALL ZuobiaoRot(Seta,Eta,NV(1),NV(2),NV(3),
     &		   TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))	
	Seta0 = Seta
	Eta0 = Eta
	RETURN
 
	END
C======================================================================C
	SUBROUTINE  ReflectAngleMC(Seta,Eta,RanSeta,RanEta,NVout,WTYPE)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.1415926)
	REAL*8	NVout(3),TV1(3),TV2(3),NVin(3),Vin(3)
	INTEGER Wtype

!this is a subroutine to give the reflecting direction of an incident
!light bundle to a wall surface.


!WTYPE:		1,	diffuse-wall 
!		2,  mirror-wall
!		3,	directional-diffuse wall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NVout is the outward normal direction
!NVin is the inward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	NVin(:)=-NVout(:)
	IF(WTYPE.EQ.1) THEN	
	Seta = RanSeta*2.0*PAI
	Eta = ACOS(SQRT(1.0-RanEta))
!-------
!Taniguchi'book, P.49~50
!For diffuse surface :	f(Eta) = 2*cos(eta)*sin(eta)* deta
!RanEta = integrate(f(eta))|0~eta = 1.0-cos(eta)**2
	CALL BuildZuoBiao(NVin,TV1,TV2,3)
	CALL ZuobiaoRot(Seta,Eta,NVin(1),NVin(2),NVin(3),
     &			TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
				

	ELSE IF(WTYPE.EQ.2) THEN

		Vin(1) = sin(Eta)*cos(Seta)
		Vin(2) = sin(Eta)*sin(Seta)
		Vin(3) = cos(Eta)
				
	f1 = Vin(1)*NVout(1)+Vin(2)*NVout(2)+Vin(3)*NVout(3)
!f2 = NVout(1)*NVout(1) + NVout(2)*NVout(2)+NVout(3)*NVout(3)
	f2 = 1.0

	Vin(:) = Vin(:)-2.0*f1/f2*NVout(:)
	f3=SQRT(Vin(1)*Vin(1)+Vin(2)*Vin(2)+Vin(3)*Vin(3))
	IF(f3.EQ.0) STOP 'Error in ReflectAngleMC: mirror reflect'
	Vin(:) = Vin(:)/f3
	
!							    
!				Vin		|	  Vin`
!				      `		|  	  /
!				       `	|    /
!				        `	|   /
!				         `	!  /
!					      `	| /
!                                   ____________|/_________
!						|
!						|
!						|
!						NVout
		Eta = acos(Vin(3))
		f4 = DSQRT(Vin(1)*Vin(1)+Vin(2)*Vin(2))
		IF(f4.GT.0) THEN
			Seta = acos(Vin(1)/f4)
		ELSE
			Seta = 0
		END IF
		IF(Vin(2).LT.0.0) Seta = -Seta
	ELSE IF(WTYPE.EQ.3) THEN	
		Seta = RanSeta*2.0*PAI
!Eta = ACOS(SQRT(1.0-RanEta))
!-------
!Taniguchi'book, P.49~50
!For directional-diffuse surface :	f(Eta) = 2*cos(eta)*sin(eta)*cos(eta)* deta
!RanEta = integrate(f(eta))|0~eta = (2/3-cos(eta)/2-cos(3*eta)/6)*1.5
		eta1 = 0.0
		eta2 = 0.5*PAI
		eta = 0.5*(eta1+eta2)
		v3 = 1.0-0.75*cos(eta)-0.25*cos(3.0*eta)
		err = abs(v3-RanEta)
		DO 300 WHILE(err.GT.1.0E-4)
		  IF(v3.LT.RanEta) THEN		!too small
			eta1 = eta
			eta = (eta + eta2)/2.0
		  ELSE
			eta2 = eta
			eta = (eta + eta1)/2.0
		  END IF
		  v3 =  1.0-0.75*cos(eta)-0.25*cos(3.0*eta)
		  err = abs(v3-RanEta)
		
300		CONTINUE
		CALL BuildZuoBiao(NVin,TV1,TV2,3)
		CALL ZuobiaoRot(Seta,Eta,NVin(1),NVin(2),NVin(3),
     &			TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
	END IF
	RETURN
!---------> end of reflectAngleMC()
	END
C======================================================================C
	SUBROUTINE  EmitAngleWall(Seta,Eta,RanSeta,RanEta,NVout,WTYPE)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.1415926)
	REAL*8	NVout(3),TV1(3),TV2(3),NVin(3),Vin(3)
	INTEGER Wtype

!this is a subroutine to give the reflecting direction of an incident
!light bundle to a wall surface.


!WTYPE:		1,	diffuse-wall
!			2,  mirror-wall, collimated emission
!			3,	directional-diffuse wall


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NVout is the outward normal direction
!NVin is the inward normal direction
!Please pay high attention!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	NVin(:)=-NVout(:)
	IF(WTYPE.EQ.1) THEN	
		Seta = RanSeta*2.0*PAI
		Eta = ACOS(SQRT(1.0-RanEta))
!-------
!Taniguchi'book, P.49~50
!For diffuse surface :	f(Eta) = 2*cos(eta)*sin(eta)* deta
!RanEta = integrate(f(eta))|0~eta = 1.0-cos(eta)**2
		CALL BuildZuoBiao(NVin,TV1,TV2,3)
		CALL ZuobiaoRot(Seta,Eta,NVin(1),NVin(2),NVin(3),
     &			TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
	ELSE IF(WTYPE.EQ.2) THEN
		Eta = acos(NVin(3))
		f4 = DSQRT(NVin(1)*NVin(1)+NVin(2)*NVin(2))
		IF(f4.GT.0) THEN
			Seta = acos(NVin(1)/f4)
		ELSE
			Seta = 0
		END IF
		IF(NVin(2).LT.0.0) Seta = -Seta
	ELSE IF(WTYPE.EQ.3) THEN	
		Seta = RanSeta*2.0*PAI
!Eta = ACOS(SQRT(1.0-RanEta))
!-------
!Taniguchi'book, P.49~50
!For directional-diffuse surface :	f(Eta) = 2*cos(eta)*sin(eta)*cos(eta)* deta
!RanEta = integrate(f(eta))|0~eta = (2/3-cos(eta)/2-cos(3*eta)/6)*1.5
		eta1 = 0.0
		eta2 = 0.5*PAI
		eta = 0.5*(eta1+eta2)
		v3 = 1.0-0.75*cos(eta)-0.25*cos(3.0*eta)
		err = abs(v3-RanEta)
		DO 300 WHILE(err.GT.1.0E-4)
		  IF(v3.LT.RanEta) THEN		!too small
	            eta1 = eta
		    eta = (eta + eta2)/2.0
	          ELSE
			eta2 = eta
			eta = (eta + eta1)/2.0
		  END IF
		  v3 =  1.0-0.75*cos(eta)-0.25*cos(3.0*eta)
		  err = abs(v3-RanEta)
300		CONTINUE
		CALL BuildZuoBiao(NVin,TV1,TV2,3)
		CALL ZuobiaoRot(Seta,Eta,NVin(1),NVin(2),NVin(3),
     &			TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))
	END IF
	RETURN
 
!---------> end of reflectAngleMC()

	END
!
	SUBROUTINE  ScatterAngleTANI(Seta0,Eta0,RanSeta,RanEta,Beta)

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI = 3.1415926)
	REAL*8	NV(3),TV1(3),TV2(3)


!this is a subroutine to give the emission angle of the scattered
!light bundle, 
!the phase function is defined at the textbook of Taniguchi et al. (1994, Page.18, 110)
!
!		fai(eta)  = 8/(3*PAI)*(sin(eta) - eta*cos(eta))
!
!
!and, RanEta = integrate_[fai(eta)]	| (0--eta)
!
!RanEta = 4/(3Pai)*[eta/2 - 3/8*sin(2eta) + eta*cos(2eta)/4]		---> get eta
	
!----------------------- ATTENTION!! --------------------------
!the integration method is different form above. It seems not valid but
!we still apply it.
!--------------------------------------------------------------

	
!seta0, eta0,	incident angle (circum,cone)
!seta, eta,	scattering bundle angle (circum,cone), ==> seta0,eta0


	seta = 2.0*PAI*RanSeta
	
	
	!unisotropic scattering, get eta's value by dichotomy scheme
	!the phase function increases monotonously with angle eta
		
	eta1 = 0.0
	eta2 = PAI
	
	eta = (eta1+eta2)/2
	v3=
     &   (0.5*eta-0.375*sin(2.*eta)+0.25*eta*cos(2.0*eta))*4.0/3.0/Pai
! --- O		
	err = abs(v3-RanEta)

	DO 100 WHILE(err.GT.1.0E-4)
		
		  IF(v3.LT.RanEta) THEN		!too small
			eta1 = eta
			eta = (eta + eta2)/2.0
		  ELSE
			eta2 = eta
			eta = (eta + eta1)/2.0
		  END IF

		  v3 = (0.5*eta-0.375*sin(2.*eta)+0.25*eta*cos(2.*eta))
     &		   *4./3./Pai
		  err = abs(v3-RanEta)
	
100	CONTINUE

	!to make a coordinate system for the incidence direction
	NV(1) = sin(Eta0)*cos(Seta0)
	NV(2) = sin(Eta0)*sin(Seta0)
	NV(3) = cos(Eta0)
	CALL BuildZuoBiao(NV,TV1,TV2,3)
	
	!get the angles (eta1,seta1) in original coodinates
	CALL ZuobiaoRot(Seta,Eta,NV(1),NV(2),NV(3),
     &		   TV1(1),TV1(2),TV1(3),TV2(1),TV2(2),TV2(3))	
	

	Seta0 = Seta
	Eta0 = Eta

	RETURN
 
 !-------->

	END
C======================================================================C	

	SUBROUTINE   GetBlockExtin(IUT0,NB,NBPG)

	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE MODULE_RADSXF, ONLY : GAbsorb,PAbsorb,PScatter
	USE module_radwork, ONLY : VWK3D

	IMPLICIT REAL*8(A-H,O-Z)


	VWK3D(:,:,:) = 0.d0
	DO 100 K=1,NB
	  DO 100 J=1,NB
	 	DO 100 I=1,NB

		   ID = GBlockIndex(i,j,k)
		   NUM = GBlockNum(I,J,K)
		   DO JJ = 1, NUM
			  IE=GBlockID(ID+JJ)
			  VWK3D(I,J,K)=VWK3D(I,J,K)
     &			  +GAbsorb(IE)+PAbsorb(IE)+PScatter(IE)
		   END DO
		   IF(NUM.GT.0) VWK3D(I,J,K)=VWK3D(I,J,K)/NUM

100	CONTINUE
	RETURN
	END
C======================================================================C
	SUBROUTINE  WallZone(IUT0,NP,NE,MP,ME,NW,NWT,
     2		NKIND,NFACE,N2D,N3D,NDOF,NRD,NRD0,MRD,
     3		XYZ,ELEM,RADIUS,XiaoJuli,AbsorbCoef,IFlagGP)

	USE MODULE_RADSXF, ONLY : ElemKind,FaceNum,NNPE,NNPS,CNTIVITY
	USE MODULE_RADSXF, ONLY : Volume,WElem,WArea,WABCD,WNV,CTR
	USE MODULE_RADSXF, ONLY : RDN1,RDN2,RDIndex,RDValue,VIndex,PROP

	IMPLICIT REAL*8(A-H,O-Z)

	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)

	INTEGER ELEM(N3D,ME)
	REAL*8 XYZ(NDOF,MP),RADIUS(NRD0)
	REAL*8 AbsorbCoef

	
	!Private
	REAL*8	NODE(NDOF,N3D),NODEP(NDOF,N3D),fun(N3D),ExArea(20)
      REAL*8	V0(3),VT(3),VW(3),VLast(3),DV(3)
	REAL*8	NV(3),NV1(3),NV2(3),VE(3,20)
	REAL*8	CTRme(3)



	DO 1000 IW = 1, NW
		IIT = NE + IW
		IF(VIndex(IIT).LT.1) CYCLE
		!WRITE(IUT0,*) 'IW = ',IW
		IF(MOD(IW,100).EQ.1) 
     &		WRITE(IUT0,'("     IW=",I8," /[",I8,"]") ') IW-1,NW
						
C--------------------
C	Pre-treat
C--------------------
    
		!CTR = 0.0
		DO 10 J =1,WELEM(1,IW)
			IP = WELEM(J+1,IW)
			DO 10 IDOF = 1,NDOF
			  NODE(IDOF,J) = XYZ(IDOF,IP)
			  
10		CONTINUE

		!Set Emitting Points
		IF(WELEM(1,IW).EQ.3) THEN		!trianglar
			NG = 4
			DO I = 1,3
				VE(:,I) = NODE(:,I)
			END DO
			VE(:,4) = (NODE(:,1)+NODE(:,2)+NODE(:,3))/3.0
						
		ELSE IF(WELEM(1,IW).EQ.4) THEN	!qudrilateral
			NG = 4
			ps = SQRT(1.0/3.0)*0.5

			fun(1) = 0.5-Ps
			fun(2) = 0.5-Ps
			CALL Qua4Pos(NODE,N3D,V0,fun,NDOF)
			VE(:,1) = V0(:)

			fun(1) = 0.5+Ps
			fun(2) = 0.5-Ps
			CALL Qua4Pos(NODE,N3D,V0,fun,NDOF)
			VE(:,2) = V0(:)

			fun(1) = 0.5-Ps
			fun(2) = 0.5+Ps
			CALL Qua4Pos(NODE,N3D,V0,fun,NDOF)
			VE(:,3) = V0(:)

			fun(1) = Ps+0.5
			fun(2) = Ps+0.5
			CALL Qua4Pos(NODE,N3D,V0,fun,NDOF)
			VE(:,4) = V0(:)
		
		ELSE
			NG = WELEM(1,IW)
			DO I = 1,WELEM(1,IW)
			  VE(:,I) = CTR(:,IIT)
			  VE(:,I) = VE(:,I) + NODE(:,I)
			  I1 = (I+1)		!-((I+1)/NG)*NG
			  IF(I1.GT.NG) I1=I1-NG

			  I2 = (I+NG-1)		!-((I+NG-1)/NG)*NG
			  IF(I2.GT.NG) I2=I2-NG

			  VE(:,I) = VE(:,I) + NODE(:,I1)
			  VE(:,I) = VE(:,I) + NODE(:,I2)

			  VE(:,I) = VE(:,I)/4.0
			END DO
		
		END IF




		
C---------------------------------------
C	Calculate Direct-Exchange Area
C---------------------------------------
	
!-----1:	Gas Objects

		NV1(:) = - WNV(:,IW)

		DO 500 JE = 1, NE
		
		  IF(VIndex(JE).LT.1) CYCLE
		  
		  !NV is dummy
		  
		  CALL GetDistance(CTR(1,IIT),CTR(2,IIT),CTR(3,IIT),
     &			   CTR(1,JE),CTR(2,JE),CTR(3,JE),dist)
		  		  
		  
		  MarkAbs = 0
		  IF(dist/RADIUS(JE).LE.XiaoJuli) THEN
			MarkAbs = 1
!-----------------------------------------------------------------
!The two elements are close or the target element is very large,
!It is need to get Gauss-Polymal-Integration to improve the accuracy
!of the calculation.
!-----------------------------------------------------------------
			
			DO 20 J =1,NNPE(ElemKind(JE))
				IP = ELEM(J,JE)
				DO 20 IDOF = 1,NDOF
					NODEP(IDOF,J) = XYZ(IDOF,IP)
20			CONTINUE

		  END IF

		  MarkEmi = 0
c		  IF(dist/RADIUS(IIT).LE.XiaoJuli) MarkEmi = 1
			!-----------------------------------------------------------------
			!The two elements are close or the emitting element is very large,
			!It is need to get Gauss-Polymal-Integration to improve the accuracy
			!of the calculation.
			!-----------------------------------------------------------------

		  IF(MarkEmi.EQ.1.AND.MarkAbs.EQ.1) THEN
			!Double Gauss Integration
			
			!absorbing element			-JE
			DO 100 IG = 1,NG
				V0(:) = VE(:,IG)
				CTRme(:) = CTR(:,JE)
				
			IF(NNPE(ElemKind(JE)).EQ.4) THEN
			CALL GaussTet4(ExArea(IG),NODEP,N3D,V0,fun,
     &			  AbsorbCoef,VOLUME(JE),NV1,1,1)
		
			ELSE IF(NNPE(ElemKind(JE)).EQ.8) THEN
		CALL GaussHex8(ExArea(IG),NODEP,N3D,V0,CTRme,fun,
     &			  CNTIVITY,N2D,NFACE,NKIND,ElemKind(JE),
     &			  AbsorbCoef,NV1,1,1)

			ELSE
			ik = ElemKind(JE)
		CALL GaussVolQita(ExArea(IG),NODEP,N3D,V0,CTRme,fun,
     &		  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &		  FaceNum(ik),IK,AbsorbCoef,VOLUME(JE),NV1,1,1)
				END IF 
  100			CONTINUE
!Emitting element		-IIT,	V0,NV1,NV2 are dummys
			CTRme(:) = CTR(:,IIT)
						
			IF(WElem(1,IW).EQ.3) THEN
			CALL GaussTri3(VALUE,NODE,N3D,V0,ExArea,
     &				  AbsorbCoef,WArea(IW),0,NV1,NV2,0,0)
		
			ELSE IF(WElem(1,IW).EQ.4) THEN
			CALL GaussSibian(VALUE,NODE,N3D,V0,ExArea,
     &		  			 AbsorbCoef,0,NV1,NV2,0,0)
	
			ELSE
			CALL GaussFaceQita(VALUE,NODE,N3D,
     &			  WElem(1,IW),V0,ExArea,AbsorbCoef,
     &			  0,NV1,NV2,0,0)
	
			END IF
			
		  ELSE IF(MarkEmi.EQ.1) THEN
!Gauss Integration for Emissive Element,	-IIT, NV2 is dummy
			V0(:) =  CTR(:,JE)
			CTRme(:) =  CTR(:,IIT)
					
			IF(WElem(1,IW).EQ.3) THEN
				CALL GaussTri3(VALUE,NODE,N3D,V0,fun,
     &				  AbsorbCoef,WArea(IW),1,NV1,NV2,1,0)
		
			ELSE IF(WElem(1,IW).EQ.4) THEN
				CALL GaussSibian(VALUE,NODE,N3D,V0,fun,
     &		  			 AbsorbCoef,1,NV1,NV2,1,0)
	
			ELSE
				CALL GaussFaceQita(VALUE,NODE,N3D,
     &				  WElem(1,IW),V0,fun,AbsorbCoef,
     &				  1,NV1,NV2,1,0)
	
			END IF
			
			VALUE = VALUE*Volume(JE)

		  ELSE IF(MarkAbs.EQ.1) THEN
!Gauss Integration for Absorbing Element,	-JE,
			V0(:) =  CTR(:,IIT)
			CTRme(:) =  CTR(:,JE)
			
			IF(NNPE(ElemKind(JE)).EQ.4) THEN
				CALL GaussTet4(VALUE,NODEP,N3D,V0,fun,
     &				  AbsorbCoef,VOLUME(JE),NV1,1,1)
		
			ELSE IF(NNPE(ElemKind(JE)).EQ.8) THEN
			CALL GaussHex8(VALUE,NODEP,N3D,V0,CTRme,fun,
     &			  CNTIVITY,N2D,NFACE,NKIND,ElemKind(JE),
     &			  AbsorbCoef,NV1,1,1)

			ELSE
			ik = ElemKind(JE)
			CALL GaussVolQita(VALUE,NODEP,N3D,V0,CTRme,fun,
     &			  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &			  FaceNum(ik),IK,AbsorbCoef,VOLUME(JE),NV1,1,1)
			END IF
			
			VALUE = VALUE*WArea(IW)

		  ELSE
			!Direct Integration
			VALUE = exp(-dist*AbsorbCoef)/(dist*dist)
			VALUE = VALUE*WArea(IW)*Volume(JE)

	CALL CosValue(NV1(1),NV1(2),NV1(3), CTR(1,JE)-CTR(1,IIT),
     &		  CTR(2,JE)-CTR(2,IIT),	CTR(3,JE)-CTR(3,IIT),
     &		  cosv)
	cosv = max(cosv,0.0)
	VALUE = VALUE*cosv
	  END IF
!Exchange area
	  VALUE = VALUE/PAI*AbsorbCoef	
	  IF(IFlagGP.NE.1) THEN
		  IRD = VIndex(IIT)
		  JRD = VIndex(JE)
		  value2 = RDValue(RDIndex(IRD)+JRD)
		  CALL SymmTreatRD(VALUE,VALUE2,PROP(IIT),PROP(JE))
			  RDValue(RDIndex(IRD)+JRD) = VALUE
		  ELSE		!grouped mode
			  IGP = VIndex(IIT)
			  JGP = VIndex(JE)
			  ID = (IGP-1)*NRD+JGP
			  RDValue(ID) = RDValue(ID)+VALUE
		  END IF
		  		
500		CONTINUE


	
!-----1:	Wall Objects

		NV1(:) = -WNV(:,IW)

		DO 800 JW = 1, NW
		
		  JJT =  NE+JW
		  NV2(:) = -WNV(:,JW)
			!here NV is INWARD normal vector

		  IF(VIndex(JJT).LT.1) CYCLE
		  
		  CALL GetDistance(CTR(1,JJT),CTR(2,JJT),CTR(3,JJT),
     &			   CTR(1,IIT),CTR(2,IIT),CTR(3,IIT),dist)
		  
		  MarkAbs = 0
		  IF(dist/RADIUS(JJT).LE.XiaoJuli) THEN
			MarkAbs = 1
!-----------------------------------------------------------------
!The two elements are close or the target element is very large,
!It is need to get Gauss-Polymal-Integration to improve the accuracy
!of the calculation.
!-----------------------------------------------------------------
		DO 520 J =1,WElem(1,JW)
		IP = WElem(J+1,JW)
		DO 520 IDOF = 1,NDOF
		NODEP(IDOF,J) = XYZ(IDOF,IP)
 520		CONTINUE
	        END IF
	        MarkEmi = 0	
!c 		IF(dist/RADIUS(IIT).LE.XiaoJuli) MarkEmi = 1
!-----------------------------------------------------------------
!The two elements are close or the emitting element is very large,
!It is need to get Gauss-Polymal-Integration to improve the accuracy
!of the calculation.
!-----------------------------------------------------------------
	  IF(MarkEmi.EQ.1.AND.MarkAbs.EQ.1) THEN
!Double Gauss Integration
		DO 600 IG = 1,NG	!absorbing element, JW
		V0(:) = VE(:,IG)
		IF(WElem(1,JW).EQ.3) THEN
		CALL GaussTri3(ExArea(IG),NODEP,N3D,V0,fun,
     &		AbsorbCoef,WArea(JW),1,NV2,NV1,1,1)
		ELSE IF(WElem(1,JW).EQ.4) THEN
		CALL GaussSibian(ExArea(IG),NODEP,N3D,V0,fun,
     &	        AbsorbCoef,1,NV2,NV1,1,1)
		ELSE
		CALL GaussFaceQita(ExArea(IG),NODEP,N3D,
     &		  WElem(1,JW),V0,fun,AbsorbCoef,
     &		  1,NV2,NV1,1,1)
		END IF
600			CONTINUE
!Emitting element		-IW,	V0,NV1,NV2 are dummys
c			V0(:) = CTR(:,IIT)
c			CTRme(:) = CTR(:,IE)
		IF(WElem(1,IW).EQ.3) THEN
		CALL GaussTri3(VALUE,NODE,N3D,V0,ExArea,
     &			  AbsorbCoef,WArea(IW),0,NV2,NV1,0,0)
		ELSE IF(WElem(1,IW).EQ.4) THEN
		CALL GaussSibian(VALUE,NODE,N3D,V0,ExArea,
     &		AbsorbCoef,0,NV2,NV1,0,0)
	
			ELSE
				CALL GaussFaceQita(VALUE,NODE,N3D,
     &				  WElem(1,IW),V0,ExArea,AbsorbCoef,
     &				  0,NV2,NV1,0,0)
	
			END IF
			
		  ELSE IF(MarkEmi.EQ.1) THEN
!Gauss Integration for Emissive Element,	-IW
			V0(:) =  CTR(:,JJT)
					
			IF(WElem(1,IW).EQ.3) THEN
			CALL GaussTri3(VALUE,NODE,N3D,V0,fun,
     &			  AbsorbCoef,WArea(IW),1,NV1,NV2,1,1)
		
			ELSE IF(WElem(1,IW).EQ.4) THEN
				CALL GaussSibian(VALUE,NODE,N3D,V0,fun,
     &		  			 AbsorbCoef,1,NV1,NV2,1,1)
	
			ELSE
				CALL GaussFaceQita(VALUE,NODE,N3D,
     &				  WElem(1,IW),V0,fun,AbsorbCoef,
     &				  1,NV1,NV2,1,1)
	
			END IF
			
			VALUE = VALUE*WArea(JW)

		  ELSE IF(MarkAbs.EQ.1) THEN
!Gauss Integration for Absorbing Element,	-IIT,	NV2 is dummy
			V0(:) =  CTR(:,IIT)
			
			IF(WElem(1,JW).EQ.3) THEN
			CALL GaussTri3(VALUE,NODEP,N3D,V0,fun,
     &			  AbsorbCoef,WArea(JW),1,NV2,NV1,1,1)
		
			ELSE IF(WElem(1,JW).EQ.4) THEN
			CALL GaussSibian(VALUE,NODEP,N3D,V0,fun,
     &			  		 AbsorbCoef,1,NV2,NV1,1,1)
	
			ELSE
			CALL GaussFaceQita(VALUE,NODEP,N3D,
     &			  WElem(1,JW),V0,fun,AbsorbCoef,
     &			  1,NV2,NV1,1,1)
	
			END IF

			VALUE = VALUE*WArea(IW)

		  ELSE
!Direct Integration
			VALUE = exp(-dist*AbsorbCoef)/(dist*dist)
			VALUE = VALUE*WArea(JW)*WArea(IW)
	CALL CosValue(NV1(1),NV1(2),NV1(3), CTR(1,JJT)-CTR(1,IIT),
     &  CTR(2,JJT)-CTR(2,IIT),CTR(3,JJT)-CTR(3,IIT),
     &					  cosv1)
	cosv1 = max(cosv1,0.0)
	CALL CosValue(NV2(1),NV2(2),NV2(3), CTR(1,IIT)-CTR(1,JJT),
     &		  CTR(2,IIT)-CTR(2,JJT),CTR(3,IIT)-CTR(3,JJT),
     &		  cosv2)
	cosv2 = max(cosv2,0.0)
	VALUE = VALUE*cosv1*cosv2
        END IF
!save Exchange-area
		  IF(IFlagGP.NE.1) THEN
			  IRD = VIndex(IIT)
			  JRD = VIndex(JJT)
			  IF(IRD.EQ.JRD) THEN
				RDValue(RDIndex(IRD)+IRD) = VALUE
			  ELSE IF(JRD.LT.IRD) THEN
				value2 = RDValue(RDIndex(IRD)+JRD)
	CALL  SymmTreatRD(VALUE,VALUE2,PROP(IIT),PROP(JJT))
				RDValue(RDIndex(IRD)+JRD) = VALUE

			  ELSE
				RDValue(RDIndex(JRD)+IRD) = VALUE
			  END IF
			ELSE
			  
			  IGP = VIndex(IIT)
			  JGP = VIndex(JJT)
			  ID = (IGP-1)*NRD+JGP
			  RDValue(ID) = RDValue(ID)+VALUE
			END IF
		
 800		CONTINUE
 1000  CONTINUE
       RETURN
       END 
C======================================================================C
	SUBROUTINE  GasZone(IUT0,NP,NE,MP,ME,NW,NWT,
     2		NKIND,NFACE,N2D,N3D,NDOF,NRD,NRD0,MRD,
     3		XYZ,ELEM,RADIUS,XiaoJuli,AbsorbCoef,IFlagGP)


	USE MODULE_RADSXF, ONLY : ElemKind,FaceNum,NNPE,NNPS,CNTIVITY
	USE MODULE_RADSXF, ONLY : Volume,WElem,WArea,WABCD,WNV,CTR
	USE MODULE_RADSXF, ONLY : RDN1,RDN2,RDIndex,RDValue,VIndex,PROP

	IMPLICIT REAL*8(A-H,O-Z)

	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)

	INTEGER ELEM(N3D,ME)
	REAL*8 XYZ(NDOF,MP),RADIUS(NRD0)
	REAL*8 AbsorbCoef

	
	!Private
	REAL*8	NODE(NDOF,N3D),NODEP(NDOF,N3D),fun(N3D),ExArea(20)
      REAL*8	V0(3),VT(3),VW(3),VLast(3),DV(3)
	REAL*8	NV(3),NV1(3),NV2(3),VE(3,20)
	REAL*8	CTRme(3)




	DO 1000 IE = 1, NE
		
		IF(VIndex(IE).LT.1) CYCLE
!WRITE(IUT0,*) 'IE = ',IE
		IF(MOD(IE,100).EQ.1) 
     &	WRITE(IUT0,'("     IE=",I8," /[",I8,"]") ') IE-1,NE	
! --- 
C--------------------
C	Pre-treat
C--------------------
    
		!CTR = 0.0
		DO 10 J =1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			DO 10 IDOF = 1,NDOF
			  NODE(IDOF,J) = XYZ(IDOF,IP)
			  !CTR(IDOF)=CTR(IDOF) + XYZ(IDOF,IP)
10		CONTINUE
		!CTR = CTR / NNPE(ElemKind(IE))

		!Set Emitting Points
		IF(NNPE(ElemKind(IE)).EQ.4) THEN
			NG = 8
			DO I = 1,4
				VE(:,I) = NODE(:,I)
			END DO
			DO I = 1,4
				I1 = (I+1)
				IF(I1.GT.4)	I1 = I1-4
				I2 = (I+2)
				IF(I2.GT.4)	I2 = I2-4
		VE(:,I+4) = (NODE(:,I)+NODE(:,I1)+NODE(:,I2))/3.0
			END DO
			
		ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
			NG = 6
			DO I = 1,6
			  VE(:,I) = 0.0
			  DO J = 1,4
				IK = ElemKind(IE)
				IP = CNTIVITY(IK,I,J)
				VE(:,I) = VE(:,I) + NODE(:,IP)
			  END DO
			  VE(:,I) = VE(:,I)/4.0
			END DO

		ELSE
			NG = NNPE(ElemKind(IE))+1
			VE(:,1) = CTR(:,IE)
			DO I = 1,NNPE(ElemKind(IE))
				VE(:,I+1) = NODE(:,I)
			END DO
		
		END IF




		
C---------------------------------------
C	Calculate Direct-Exchange Area
C---------------------------------------
	
!-----1:	Gas Objects

		DO 500 JE = 1, NE
		
		  IF(VIndex(JE).LT.1) CYCLE
!		  IF(JE.EQ.IE) CYCLE
		  
		  !NV is dummy
		  
		  CALL GetDistance(CTR(1,IE),CTR(2,IE),CTR(3,IE),
     &		   CTR(1,JE),CTR(2,JE),CTR(3,JE),dist)
		  		  
		  
		  MarkAbs = 0
		  IF(dist/RADIUS(JE).LE.XiaoJuli) THEN
			MarkAbs = 1
!-----------------------------------------------------------------
!The two elements are close or the target element is very large,
!It is need to get Gauss-Polymal-Integration to improve the accuracy
!of the calculation.
!-----------------------------------------------------------------
			
			DO 20 J =1,NNPE(ElemKind(JE))
				IP = ELEM(J,JE)
				DO 20 IDOF = 1,NDOF
					NODEP(IDOF,J) = XYZ(IDOF,IP)
20			CONTINUE

		  END IF

		  MarkEmi = 0	
c		  IF(dist/RADIUS(IE).LE.XiaoJuli) MarkEmi = 1
!-----------------------------------------------------------------
!The two elements are close or the emitting element is very large,
!It is need to get Gauss-Polymal-Integration to improve the accuracy
!of the calculation.
!-----------------------------------------------------------------
		  		  
		  IF(MarkEmi.EQ.1.AND.MarkAbs.EQ.1) THEN
			!Double Gauss Integration
			!NV is always dummy

			!absorbing element			-JE
			DO 100 IG = 1,NG
				V0(:) = VE(:,IG)
				CTRme(:) = CTR(:,JE)

				IF(NNPE(ElemKind(JE)).EQ.4) THEN
			CALL GaussTet4(ExArea(IG),NODEP,N3D,V0,fun,
     &			  AbsorbCoef,VOLUME(JE),NV,1,0)
		
			ELSE IF(NNPE(ElemKind(JE)).EQ.8) THEN
		CALL GaussHex8(ExArea(IG),NODEP,N3D,V0,CTRme,fun,
     &		  CNTIVITY,N2D,NFACE,NKIND,ElemKind(JE),
     &			  AbsorbCoef,NV,1,0)
				
				ELSE
		ik = ElemKind(JE)
		CALL GaussVolQita(ExArea(IG),NODEP,N3D,V0,CTRme,fun,
     &		  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &		  FaceNum(ik),IK,AbsorbCoef,VOLUME(JE),NV,1,0)
		
				END IF

100			CONTINUE
!Emitting element		-IE,	V0,CTRme are dummys
			CTRme(:) = CTR(:,IE)
						
			IF(NNPE(ElemKind(IE)).EQ.4) THEN
		CALL GaussTet4(VALUE,NODE,N3D,V0,ExArea,
     &		  AbsorbCoef,VOLUME(IE),NV,0,0)
		
		ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
		CALL GaussHex8(VALUE,NODE,N3D,V0,CTRme,ExArea,
     &		  CNTIVITY,N2D,NFACE,NKIND,ElemKind(IE),
     &		  AbsorbCoef,NV,0,0)
		ELSE
!v0,v1 are dummy
		ik = ElemKind(IE)
		CALL GaussVolQita(VALUE,NODE,N3D,V0,CTRme,ExArea,
     &		  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &		  FaceNum(ik),IK,AbsorbCoef,VOLUME(IE),NV,0,0)
		END IF
	  ELSE IF(MarkEmi.EQ.1) THEN
!Gauss Integration for Emissive Element,	-IE
	V0(:) =  CTR(:,JE)
	CTRme(:) =  CTR(:,IE)
	IF(NNPE(ElemKind(IE)).EQ.4) THEN
	CALL GaussTet4(VALUE,NODEP,N3D,V0,fun,
     &				  AbsorbCoef,VOLUME(IE),NV,1,0)
		
	ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
		CALL GaussHex8(VALUE,NODE,N3D,V0,CTRme,fun,
     &		  CNTIVITY,N2D,NFACE,NKIND,ElemKind(IE),
     &			  AbsorbCoef,NV,1,0)
	ELSE
		ik = ElemKind(IE)
		CALL GaussVolQita(VALUE,NODE,N3D,V0,CTRme,fun,
     &		  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &		  FaceNum(ik),IK,AbsorbCoef,VOLUME(IE),NV,1,0)
	END IF
	VALUE = VALUE*Volume(JE)
        ELSE IF(MarkAbs.EQ.1) THEN
!Gauss Integration for Absorbing Element,	-JE
	V0(:) =  CTR(:,IE)
	CTRme(:) =  CTR(:,JE)
	IF(NNPE(ElemKind(JE)).EQ.4) THEN
		CALL GaussTet4(VALUE,NODEP,N3D,V0,fun,
     &		  AbsorbCoef,VOLUME(JE),NV,1,0)
	ELSE IF(NNPE(ElemKind(JE)).EQ.8) THEN
	CALL GaussHex8(VALUE,NODEP,N3D,V0,CTRme,fun,
     &	  CNTIVITY,N2D,NFACE,NKIND,ElemKind(JE),
     &		  AbsorbCoef,NV,1,0)

			ELSE
			ik = ElemKind(JE)
			CALL GaussVolQita(VALUE,NODEP,N3D,V0,CTRme,fun,
     &			  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &			  FaceNum(ik),IK,AbsorbCoef,VOLUME(JE),NV,1,0)
			END IF
			
			VALUE = VALUE*Volume(IE)

		  ELSE
			!Direct Integration
			VALUE = exp(-dist*AbsorbCoef)/(dist*dist)
			VALUE = VALUE*Volume(IE)*Volume(JE)

		  END IF
		  
		  !Exchange area
		  VALUE = VALUE/PAI*AbsorbCoef*AbsorbCoef
	  
		  IF(IFlagGP.NE.1) THEN
			  IRD = VIndex(IE)
			  JRD = VIndex(JE)
			  IF(JE.EQ.IE) THEN
				RDValue(RDIndex(IRD)+IRD) = VALUE
			  ELSE IF(JE.LT.IE) THEN
		value2 = RDValue(RDIndex(IRD)+JRD)
		CALL  SymmTreatRD(VALUE,VALUE2,PROP(IE),PROP(JE))
		RDValue(RDIndex(IRD)+JRD) = VALUE

			  ELSE
				RDValue(RDIndex(JRD)+IRD) = VALUE
			  END IF
		  ELSE
			  IGP = VIndex(IE)
			  JGP = VIndex(JE)
			  ID = (IGP-1)*NRD+JGP
			  RDValue(ID) = RDValue(ID)+VALUE
		  END IF	
500		CONTINUE




	
!-----1:	Wall Objects

		DO 800 IW = 1, NW
		
		  IIT =  NE+IW
		  NV1(:) = -WNV(:,IW)
			!here NV is INWARD normal vector

		  IF(VIndex(IIT).LT.1) CYCLE
		  
		  CALL GetDistance(CTR(1,IE),CTR(2,IE),CTR(3,IE),
     &			   CTR(1,IIT),CTR(2,IIT),CTR(3,IIT),dist)
		  
		  MarkAbs = 0
		  IF(dist/RADIUS(IIT).LE.XiaoJuli) THEN
			MarkAbs = 1
!-----------------------------------------------------------------
!The two elements are close or the target element is very large,
!It is need to get Gauss-Polymal-Integration to improve the accuracy
!of the calculation.
!-----------------------------------------------------------------
			
			DO 520 J =1,WElem(1,IW)
				IP = WElem(J+1,IW)
				DO 520 IDOF = 1,NDOF
					NODEP(IDOF,J) = XYZ(IDOF,IP)
520			CONTINUE

		  END IF

		  MarkEmi = 0	
c		  IF(dist/RADIUS(IE).LE.XiaoJuli) MarkEmi = 1
!-----------------------------------------------------------------
!The two elements are close or the emitting element is very large,
!It is need to get Gauss-Polymal-Integration to improve the accuracy
!of the calculation.
!-----------------------------------------------------------------
				  
		  IF(MarkEmi.EQ.1.AND.MarkAbs.EQ.1) THEN
!Double Gauss Integration
!aborbting element			-JE,NV2 is dummy
!NV1(:) = -WNV(:,IW)
			DO 600 IG = 1,NG
				V0(:) = VE(:,IG)
				
				IF(WElem(1,IW).EQ.3) THEN
			CALL GaussTri3(ExArea(IG),NODEP,N3D,V0,fun,
     &			  AbsorbCoef,WArea(IW),1,NV1,NV2,1,0)
		
			ELSE IF(WElem(1,IW).EQ.4) THEN
			CALL GaussSibian(ExArea(IG),NODEP,N3D,V0,fun,
     &		  			 AbsorbCoef,1,NV1,NV2,1,0)
	
		ELSE
			CALL GaussFaceQita(ExArea(IG),NODEP,N3D,
     &			  WElem(1,IW),V0,fun,AbsorbCoef,
     &			  1,NV1,NV2,1,0)
	
		END IF
600			CONTINUE
!Emitting element		-IE,	V0,V1,NV1 are dummys
!V0(:) = CTR(:,IIT)
!CTRme(:) = CTR(:,IE)
		IF(NNPE(ElemKind(IE)).EQ.4) THEN
		CALL GaussTet4(VALUE,NODE,N3D,V0,ExArea,
     &				  AbsorbCoef,VOLUME(IE),NV1,0,0)
		
			ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
		CALL GaussHex8(VALUE,NODE,N3D,V0,CTRme,ExArea,
     &		  CNTIVITY,N2D,NFACE,NKIND,ElemKind(IE),
     &				  AbsorbCoef,NV1,0,0)

			ELSE
				!v0,v1 are dummy
				ik = ElemKind(IE)
		CALL GaussVolQita(VALUE,NODE,N3D,V0,CTRme,ExArea,
     &		  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &		  FaceNum(ik),IK,AbsorbCoef,VOLUME(IE),NV1,0,0)
				
	END IF
			
		  	 
		  
	 ELSE IF(MarkEmi.EQ.1) THEN
	!Gauss Integration for Emissive Element,	-IE
	V0(:) =  CTR(:,IIT)
	CTRme(:) = CTR(:,IE)
	IF(NNPE(ElemKind(IE)).EQ.4) THEN
		CALL GaussTet4(VALUE,NODEP,N3D,V0,fun,
     &				  AbsorbCoef,VOLUME(IE),NV1,1,1)
		
	ELSE IF(NNPE(ElemKind(IE)).EQ.8) THEN
	CALL GaussHex8(VALUE,NODE,N3D,V0,CTRme,fun,
     &	  CNTIVITY,N2D,NFACE,NKIND,ElemKind(IE),
     &		  AbsorbCoef,NV1,1,1)
	ELSE
	ik = ElemKind(IE)
	CALL GaussVolQita(VALUE,NODE,N3D,V0,CTRme,fun,
     &	  NNPE(ik),CNTIVITY,NNPS,NFACE,NKIND,N2D,
     &	  FaceNum(ik),IK,AbsorbCoef,VOLUME(IE),NV1,1,1)
	END IF
	VALUE = VALUE*WArea(IW)
        ELSE IF(MarkAbs.EQ.1) THEN
!Gauss Integration for Absorbing Element,	-IW, NV2 is dummy
	V0(:) =  CTR(:,IE)
	IF(WElem(1,IW).EQ.3) THEN
	CALL GaussTri3(VALUE,NODEP,N3D,V0,fun,
     &		  AbsorbCoef,WArea(IW),1,NV1,NV2,1,0)
		
	ELSE IF(WElem(1,IW).EQ.4) THEN
	CALL GaussSibian(VALUE,NODEP,N3D,V0,fun,
     &			 AbsorbCoef,1,NV1,NV2,1,0)
	
	ELSE
	CALL GaussFaceQita(VALUE,NODEP,N3D,
     &		  WElem(1,IW),V0,fun,AbsorbCoef,
     &		  1,NV1,NV2,1,0)
	
			END IF

			VALUE = VALUE*Volume(IE)

		  ELSE
			!Direct Integration
			VALUE = exp(-dist*AbsorbCoef)/(dist*dist)
			VALUE = VALUE*Volume(IE)*WArea(IW)

	CALL CosValue(NV1(1),NV1(2),NV1(3), CTR(1,IE)-CTR(1,IIT),
     &		  CTR(2,IE)-CTR(2,IIT),	CTR(3,IE)-CTR(3,IIT),
     &			  cosv)
			cosv = max(cosv,0.0)

			VALUE = VALUE*cosv

		  END IF
		  
		  !Exchange area
		  VALUE = VALUE/PAI*AbsorbCoef
		  
		  !save Exchange-area
		  IF(IFlagGP.NE.1) THEN
			  JRD = VIndex(IIT)
			  RDValue(RDIndex(JRD)+IRD) = VALUE
		  ELSE
			  IGP = VIndex(IE)
			  JGP = VIndex(IIT)
			  ID = (IGP-1)*NRD+JGP
			  RDValue(ID) = RDValue(ID)+VALUE
		  END IF
		
800		CONTINUE



		
1000  CONTINUE

	


	RETURN

!-------> end of GasZone()

	END 
C======================================================================C
	SUBROUTINE  SelfZone
     &   (NRD,NRD0,MRD,RDValue,RDIndex,RevIndex,PROP)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)

	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)

	
	INTEGER RDIndex(NRD+1),RevIndex(NRD)

	REAL*8 RDValue(MRD),PROP(NRD0)


	DO 100 IRD = 1,NRD
			
		VVV = 0
			
		!for left part J<I, Dij  (Dij = Ai*Fij)
		DO J=1,IRD-1
			ID = RDIndex(IRD) + J		
			VVV = VVV+ RDValue(ID)
			
		END DO

!for left part J<I, Dij -> Dji
		DO J=IRD+1,NRD	
!Fii is save in RDSelf()
			ID = RDIndex(J) + IRD
			VVV = VVV+ RDValue(ID)
		END DO
		ID = RDIndex(IRD) + IRD
		RDValue(ID) = MAX(PROP(RevIndex(IRD))-VVV,0.0)
100	CONTINUE
	
	
	RETURN

!-------> end of SelfZone()

	END 

C======================================================================C
	SUBROUTINE  MakeGroup(IUT0,MG,NRD,NRD0,NGP,NBOUN,
     &		MP,NP,ME,NE,NWT,NW,NDOF,N2D,N3D,NKIND,NFACE,NMAT,
     &		XYZ,ELEM,ElemKind,FaceNum,NNPE,NNPS,CNTIVITY,
     &		MatID,FanWei,GlobeCTR,WGFlag,IFlagGroup)
C======================================================================C

	use module_material, only:  radmat,radfludflag
	
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID
	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID
	USE MODULE_RADSXF, ONLY : WElem,WArea,Volume,CTR
	USE MODULE_RADSXF, ONLY : VIndex,PROP,Albedo
	USE module_radgroup, ONLY : GpCNum,GpProp,WSID
	USE module_radwork, ONLY : IndexToGp => WK1
c	USE MODULE_RADSXF, ONLY : RDId,RDN1,RDN2,RDIndex,RDValue

	
	IMPLICIT REAL*8(A-H,O-Z)
	
	INTEGER ELEM(N3D,ME),ElemKind(NE),WGFlag(NBOUN)
	INTEGER FaceNum(NKIND),NNPE(NKIND),NNPS(NKIND,NFACE),
     &		CNTIVITY(NKIND,NFACE,N2D),MatID(NE)
	REAL*8  XYZ(NDOF,MP),FanWei(2,NDOF),GlobeCTR(3)
	INTEGER WSNO(NBOUN),IGasGpFlag,IFlagGroup
	

	INTEGER NUM1(2,NMAT),NUM2(2,NBOUN)

C------------------------------------------------------------------------------------------C
C
C	This subroutine can group up the adjacent cells or surfaces, in order to reduce 
C	the computation time and memory for too fine grids.
C	NOTES:
C			1. a group does not contain a face and a volume cell simutaneously
C			2. the surface which is marked by WGFlag(IK)=0 is not grouped
C			3. faces belong to different boundaries do not enter one group
C			---- this makes the groups contain only insignificant elements
C			
C------------------------------------------------------------------------------------------C
C
C----------------------------- 1.  Make Bucket Mesh --------------------------- 
C

	NB1 = (0.75*MG)**0.3333
	NB2 = (0.50*MG)**0.3333
	MPBG = NE*1.1
	MPBW = NW*1.1
	WRITE(IUT0,*) '  BUCKET ARRAY FOR GAS AND WALL:', NB1,NB2


	!Gas Element
	ALLOCATE(GBlockNum(NB1,NB1,NB1),GBlockIndex(NB1,NB1,NB1))
	ALLOCATE(GBlockID(MPBG))

	CALL StackVolumeGroup(IUT0,CTR,ELEM,ME,N3D,NDOF,NE,NRD0,NNPE,
     &		 NKIND,ElemKind,GBlockID,GBlockNum,GBlockIndex,
     &		 FanWei,MPBG,NPBG,NB1)
	
	
	!Wall Element
	ALLOCATE(WBlockNum(NB2,NB2,NB2),WBlockIndex(NB2,NB2,NB2))
	ALLOCATE(WBlockID(MPBW))

	CALL StackWallGroup(IUT0,CTR,WElem,NW,NE,NRD0,N2D,N3D,NDOF,
     &	   WBlockID,WBlockNum,WBlockIndex,FanWei,MPBW,NPBW,NB2)



C
C--------------- 2. Cell (Element/Surface) --> Group index ----------------- 
C		
	ALLOCATE(IndexToGP(NRD0),GpCNum(NRD), GpProp(NRD),
     &		 WSID(NBOUN,10000))
	
	WSNO=0
	WSID=0
	NGP=0
	IndexToGP = 0
	GpCNum = 0
	GpProp = 0.d0
			
	!gas
	NUM1 = 0
	
	DO 100 K=1,NB1
  	  DO 100 J=1,NB1
		DO 100 I=1,NB1
			maxf=0
			WSNO(:)=0
			DO IB=1,GBlockNum(I,J,K)
				ID = GBlockIndex(I,J,K)+IB
				IE = GBlockID(ID)
				IRD = Vindex(IE)
				IK = MatID(IE)
				IF(IK.LE.0) CYCLE
				
				NUM1(1,IK)=NUM1(1,IK)+1
				
			IF(IRD.GT.0.AND.RadFludFlag(4,IK).EQ.1) THEN
					WSNO(IK) = WSNO(IK) + 1
					WSID(IK,WSNO(IK))=IE
					maxf=MAX(maxf,IK)
				END IF
			END DO
	
			IA=0
			DO IK=1,maxf
				IF(WSNO(IK).EQ.0) CYCLE
				IA=IA+1
				NUM1(2,IK)=NUM1(2,IK)+1
				DO JD=1,WSNO(IK)
					IE = WSID(IK,JD)
					KGP = NGP+IA
					IndexToGP(IE)=KGP
				GpCNum(KGP)=GpCNum(KGP)+1
				GpProp(KGP)=GpProp(KGP)+PROP(IE)
				END DO
			END DO
			NGP=NGP+IA
 100	CONTINUE
	

	mgas=0
	DO I=1,NGP
		mgas=max(mgas,GpCNum(I))
	END DO
	NGP1=NGP
	IF(NGP1.EQ.0) WRITE(IUT0,*) '    No Fluid Cell is grouped !'
!wall
	NUM2 = 0
	DO 200 K=1,NB2
  	  DO 200 J=1,NB2
		DO 200 I=1,NB2
		maxb=0
		WSNO(:)=0
		DO IB=1,WBlockNum(I,J,K)
			ID = WBlockIndex(I,J,K)+IB
			IW = WBlockID(ID)
			IIT = NE+IW
			IRD = Vindex(IIT)
			IK = WElem(WElem(1,IW)+2,IW)
			IF(IK.LE.0) CYCLE

			NUM2(1,IK)=NUM2(1,IK)+1

			IF(IRD.GT.0.AND.WGFlag(IK).EQ.1) THEN
				WSNO(IK) = WSNO(IK) + 1
				WSID(IK,WSNO(IK))=IIT
				maxb=MAX(maxb,IK)
			END IF
		END DO
		
		IA=0
		DO IK=1,maxb
			IF(WSNO(IK).EQ.0) CYCLE
			IA=IA+1
			NUM2(2,IK)=NUM2(2,IK)+1
			DO JD=1,WSNO(IK)
				IIT = WSID(IK,JD)
				KGP = NGP+IA
				IndexToGP(IIT)=KGP
				GpCNum(KGP)=GpCNum(KGP)+1
				GpProp(KGP)=GpProp(KGP)+PROP(IIT)
			END DO
		END DO
		
		NGP=NGP+IA

200	CONTINUE

	mwall=0
	DO I=NGP1+1,NGP
		mwall=max(mwall,GpCNum(I))
	END DO
	IF(NGP-NGP1.EQ.0) 
     &    WRITE(IUT0,*) '    No Wall Surface is grouped !'

	DEALLOCATE(GBlockNum,GBlockIndex,WBlockNum,WBlockIndex)
	DEALLOCATE(GBlockID,WBlockID)
	DEALLOCATE(WSID)
	
	N1=NGP
	IFlagGroup=1

	!-------------- except treatment : no group cells
	IF(N1.EQ.0) THEN
	WRITE(IUT0,*)  	'No cell or surface has been grouped !!!!'
	WRITE(IUT0,*)  	'Grouping is resumed, RETURN !!!'
	IFlagGroup=0
	DEALLOCATE(IndexToGP,GpCNum, GpProp)
	RETURN
	END IF



	!append the remained not grouped elements
	DO I=1,NRD0
		IRD= VIndex(I)
		IF(IRD.GT.0.AND.IndexToGP(I).EQ.0) THEN
			NGP=NGP+1
			IndexToGP(I)=NGP
			GpCNum(NGP)=1
			GpProp(NGP)=PROP(I)
		END IF
	END DO

	!-------------- except treatment : no group cells
	IF(NGP.GT.0.8*NRD) THEN
	WRITE(IUT0,'(a,I5,a,I5)') '  Origen Cell=', NRD, 'By Group=', 
     &							NGP
	WRITE(IUT0,*)' Grouping can not improve efficicency,resumed!'
		
	IFlagGroup=0
	DEALLOCATE(IndexToGP,GpCNum, GpProp)
	RETURN

	END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	VIndex(:)=IndexToGP(:)
	DEALLOCATE(IndexToGP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Monitor 
	WRITE(IUT0,'(a,I5)')  	'  Total Cells = ', NGP
	WRITE(IUT0,'(a,I5)')  	'    >  Group Cells     = ', N1
	WRITE(IUT0,'(a,I5)')  	'    >  Fluid Groups    = ', NGP1
	
	
c	IF(NUM1(1,-1).GT.0) THEN
c	 WRITE(IUT0,'(a,I2,a,I7,a,I5)') 	
c     &	'      <-1>  ', NUM1(1,-1),':',NUM1(2,-1)
c	END IF

	DO I =1,NMAT
	WRITE(IUT0,'(a,I2,a,I7,a,I5)') 	
     & '      <',I,'>  ', NUM1(1,I),':',
     &	NUM1(2,I)
	END DO
	WRITE(IUT0,'(a,I5)')  	'    Largest fluid group  = ', mgas

	WRITE(IUT0,'(a,I5)')  	'    >  Wall Groups     = ', N1-NGP1
	DO I =1,NBOUN-1
	WRITE(IUT0,'(a,I3,a,I5,a,I5)') 
     &	'      <',I,'>  ', NUM2(1,I),':',
     &	NUM2(2,I)
	END DO
	WRITE(IUT0,'(a,I5)')  	'    Largest wall group   = ', mwall

	WRITE(IUT0,'(a,I5)')  	'    >  Ungrouped Cells = ', NGP-N1

	RETURN

!-------> end of MakeGroup()


	END 
C======================================================================C
	SUBROUTINE  SymmGroupRDV(IUT0,IUT1,NRD0,MRD0,NRD,MRD)
C======================================================================C
	USE module_radgroup, ONLY : GpCNum,GpProp
	USE MODULE_RADSXF, ONLY : RDIndex,RDValue
	!USE module_radwork, ONLY : WK1


	IMPLICIT REAL*8(A-H,O-Z)
	
	DO 100 I=1,NRD
	  DO J=1,I-1
		ID1=(I-1)*NRD + J
		ID2=(J-1)*NRD + I
		V1 = RDValue(ID1)
		V2 = RDValue(ID2)
		CALL  SymmTreatRD(V1,V2,GpProp(I),GpProp(J))
		RDValue(ID1)=V1
		RDValue(ID2)=V2

	  END DO
100	CONTINUE

	
	OPEN(IUT1,FILE='rdtmp.bin',
     &	 form='unformatted',status='unknown',iostat=ios)
	
	DO 200 I=1,NRD
		ID=(I-1)*NRD
		WRITE(IUT1) (RDValue(J),J=ID+1,ID+I)
	
200	CONTINUE
	
	CLOSE(IUT1)
	

	
	RDIndex = 0
	DO 600 I = 2, NRD+1
		RDIndex(I) = RDIndex(I-1) + I-1
600	CONTINUE
	MRD = RDIndex(NRD+1)+100
	MRD0 = MRD
	
	DEALLOCATE(RDValue)
	ALLOCATE(RDValue(MRD0))



	OPEN(IUT1,FILE='rdtmp.bin',
     &	 form='unformatted',status='unknown',iostat=ios)

	DO I=1,NRD
		ID=RDIndex(I)
		READ(IUT1) (RDValue(J),J=ID+1,ID+I)
	
	END DO
	
	CLOSE(IUT1)
	
	
	RETURN

!-------> end of SymmGroupRDV()


	END 

C======================================================================C
	SUBROUTINE StackVolume(IUT0,XYZ,ELEM,MP,ME,N3D,NDOF,NP,NE,NNPE,
     &			 NKIND,NMAT,ElemID,MatID,MeshRange,NPB,MB,FDA)
C======================================================================C


	USE MODULE_RADSXF, ONLY : GBlockNum,GBlockIndex,GBlockID

     	
	IMPLICIT REAL*8(A-H,O-Z)
	!PARAMETER( MaxLimit=2000, idstep=1 )
C	
	REAL*8  XYZ(NDOF,MP),MeshRange(2,NDOF)
	INTEGER ELEM(N3D,ME),NNPE(20),ElemID(NE),MatID(NE)
	CHARACTER*10  StackMode
	REAL*8	CTR(3)
	
	

10	dx = (MeshRange(2,1)-MeshRange(1,1))/MB
	dy = (MeshRange(2,2)-MeshRange(1,2))/MB
	dz = (MeshRange(2,3)-MeshRange(1,3))/MB
	maxnum = 0
	minnum = 1.0E5

	ALLOCATE(GBlockNum(MB,MB,MB),GBlockIndex(MB,MB,MB))
	GBlockNum(:,:,:) = 0
	GBlockIndex(:,:,:) = 0

		
	
	!write(*,*) 'NE=',NE,MB
	
	!count  
	DO 400 IE = 1, NE
		ID = MatID(IE)
		IF(ID.GT.NMAT.OR.ID.LE.0) CYCLE

		xmax=MeshRange(1,1)
		xmin=MeshRange(2,1)
		ymax=MeshRange(1,2)
		ymin=MeshRange(2,2)
		zmax=MeshRange(1,3)
		zmin=MeshRange(2,3)
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0

		DO 390 J = 1, NNPE(ElemID(IE))
			IP = ELEM(J,IE)
			xmax = 	max(xmax, XYZ(1,IP))
			xmin = 	min(xmin, XYZ(1,IP))
			ymax = 	max(ymax, XYZ(2,IP))
			ymin = 	min(ymin, XYZ(2,IP))
			zmax = 	max(zmax, XYZ(3,IP))
			zmin = 	min(zmin, XYZ(3,IP))

			CTR(1) = CTR(1) + XYZ(1,IP)
			CTR(2) = CTR(2) + XYZ(2,IP)
			CTR(3) = CTR(3) + XYZ(3,IP)
390		CONTINUE
		CTR(1) = CTR(1)/NNPE(ElemID(IE))
		CTR(2) = CTR(2)/NNPE(ElemID(IE))
		CTR(3) = CTR(3)/NNPE(ElemID(IE))

		xmax = 	xmax + dx*FDA
		xmin = 	xmin - dx*FDA
		ymax = 	ymax + dy*FDA
		ymin = 	ymin - dy*FDA
		zmax = 	zmax + dz*FDA
		zmin = 	zmin - dz*FDA
	
		ixmax = (xmax-MeshRange(1,1))/dx+1
		ixmin = (xmin-MeshRange(1,1))/dx+1
		iymax = (ymax-MeshRange(1,2))/dy+1
		iymin = (ymin-MeshRange(1,2))/dy+1
		izmax = (zmax-MeshRange(1,3))/dz+1
		izmin = (zmin-MeshRange(1,3))/dz+1
		
		ictr = (CTR(1)-MeshRange(1,1))/dx+1
		jctr = (CTR(2)-MeshRange(1,2))/dy+1
		kctr = (CTR(3)-MeshRange(1,3))/dz+1

		ixmax = min(MB,ixmax)
!		ixmax = min(ictr+idstep,ixmax)
		ixmin = max(1,ixmin)
!		ixmin = max(ictr-idstep,ixmin)
		iymax = min(MB,iymax)
!		iymax = min(jctr+idstep,iymax)
		iymin = max(1,iymin)
!		iymin = max(jctr-idstep,iymin)
		izmax = min(MB,izmax)
!		izmax = min(kctr+idstep,izmax)
		izmin = max(1,izmin)
!		izmin = max(kctr-idstep,izmin)
		
		num = (ixmax-ixmin)*(iymax-iymin)*(izmax-izmin)
	
		DO 395 iz = izmin, izmax
		 DO 395 iy = iymin, iymax
		  DO 395 ix = ixmin, ixmax
			GBlockNum(ix,iy,iz) = GBlockNum(ix,iy,iz) + 1
395		CONTINUE

400	CONTINUE

	

	
	!put GBlockIndex
	NPB = 0
	maxnum = GBlockNum(1,1,1)
	DO 450 K = 1, MB
		DO 450 J = 1, MB
			DO 450 I = 1, MB
			GBlockIndex(i,j,k) = NPB
			NPB = NPB + GBlockNum(i,j,k)
			maxnum = max(GBlockNum(i,j,k), maxnum)
			GBlockNum(i,j,k) = 0
		
450	CONTINUE
	

	WRITE(IUT0,*) '  maximum number in a GAS bucket box = ', maxnum
	!WRITE(IUT0,*) '  minimum number in a block = ', minnum
	WRITE(IUT0,*) '  Stacked', NPB, ' times in total'

	MPB=NPB
	ALLOCATE(GBlockID(NPB),stat=ierr)
	IF(ierr.NE.0) STOP 'allocate array in StackVolume()'


	IF(NPB.GT.20000000) THEN
		WRITE(IUT0,*) '     !!! WARNING in StackVolume()::'
		WRITE(IUT0,*) '     Too many stack element!', NPB
	END IF


	DO 490 IE = 1, NE
		ID = MatID(IE)
		IF(ID.GT.NMAT.OR.ID.LE.0) CYCLE

		xmax=MeshRange(1,1)
		xmin=MeshRange(2,1)
		ymax=MeshRange(1,2)
		ymin=MeshRange(2,2)
		zmax=MeshRange(1,3)
		zmin=MeshRange(2,3)
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0

		DO 460 J = 1, NNPE(ElemID(IE))
			IP = ELEM(J,IE)
			xmax = 	max(xmax, XYZ(1,IP))
			xmin = 	min(xmin, XYZ(1,IP))
			ymax = 	max(ymax, XYZ(2,IP))
			ymin = 	min(ymin, XYZ(2,IP))
			zmax = 	max(zmax, XYZ(3,IP))
			zmin = 	min(zmin, XYZ(3,IP))

			CTR(1) = CTR(1) + XYZ(1,IP)
			CTR(2) = CTR(2) + XYZ(2,IP)
			CTR(3) = CTR(3) + XYZ(3,IP)
460		CONTINUE
		CTR(1) = CTR(1)/NNPE(ElemID(IE))
		CTR(2) = CTR(2)/NNPE(ElemID(IE))
		CTR(3) = CTR(3)/NNPE(ElemID(IE))

		xmax = 	xmax + dx*FDA
		xmin = 	xmin - dx*FDA
		ymax = 	ymax + dy*FDA
		ymin = 	ymin - dy*FDA
		zmax = 	zmax + dz*FDA
		zmin = 	zmin - dz*FDA
	
		ixmax = (xmax-MeshRange(1,1))/dx+1
		ixmin = (xmin-MeshRange(1,1))/dx+1
		iymax = (ymax-MeshRange(1,2))/dy+1
		iymin = (ymin-MeshRange(1,2))/dy+1
		izmax = (zmax-MeshRange(1,3))/dz+1
		izmin = (zmin-MeshRange(1,3))/dz+1
		
		ictr = (CTR(1)-MeshRange(1,1))/dx+1
		jctr = (CTR(2)-MeshRange(1,2))/dy+1
		kctr = (CTR(3)-MeshRange(1,3))/dz+1

		ixmax = min(MB,ixmax)
!		ixmax = min(ictr+idstep,ixmax)
		ixmin = max(1,ixmin)
!		ixmin = max(ictr-idstep,ixmin)
		iymax = min(MB,iymax)
!		iymax = min(jctr+idstep,iymax)
		iymin = max(1,iymin)
!		iymin = max(jctr-idstep,iymin)
		izmax = min(MB,izmax)
!		izmax = min(kctr+idstep,izmax)
		izmin = max(1,izmin)
!		izmin = max(kctr-idstep,izmin)

		DO 465 iz = izmin, izmax
		 DO 465 iy = iymin, iymax
		  DO 465 ix = ixmin, ixmax

			GBlockNum(ix,iy,iz) = GBlockNum(ix,iy,iz) + 1
		num = GBlockNum(ix,iy,iz) + GBlockIndex(ix,iy,iz)
			GBlockId(num) = IE

465		CONTINUE


490	CONTINUE


!------> end of StackVolume()

	END 

C======================================================================C
	SUBROUTINE StackWall(IUT0,XYZ,WElem,MP,NP,NW,N2D,N3D,NDOF,
     &				 MeshRange,WNV,NPB,MB,FDA,wsmk,NBOUND)
C======================================================================C
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID


	IMPLICIT REAL*8(A-H,O-Z)
	!PARAMETER( MaxLimit=2000, idstep=1 )
C	
	REAL*8  XYZ(NDOF,MP),MeshRange(2,NDOF),WNV(NDOF,NW)
	INTEGER WElem(N2D+2,NW),WSMK(-NBOUND:NBOUND)
	CHARACTER*10  StackMode
	REAL*8	CTR(3),Gap
	
	
10	dx = (MeshRange(2,1)-MeshRange(1,1))/MB
	dy = (MeshRange(2,2)-MeshRange(1,2))/MB
	dz = (MeshRange(2,3)-MeshRange(1,3))/MB
	Gap = MAX(MAX(dx,dy),dz)*1.2
	maxnum = 0
	minnum = 1.0E5
	

	ALLOCATE(WBlockNum(MB,MB,MB),WBlockIndex(MB,MB,MB))

	WBlockNum(:,:,:) = 0
	WBlockIndex(:,:,:) = 0

	

	
	!write(*,*) 'NE=',NE,MB
	
	!count 
	DO 400 IW = 1, NW
		IB = WELEM(WElem(1,IW)+2,IW)
		IF(WSMK(IB).EQ.0) CYCLE

		xmax=MeshRange(1,1)
		xmin=MeshRange(2,1)
		ymax=MeshRange(1,2)
		ymin=MeshRange(2,2)
		zmax=MeshRange(1,3)
		zmin=MeshRange(2,3)
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0

		DO 390 J = 1, WElem(1,IW)
			IP = WElem(J+1,IW)
			xmax = 	max(xmax, XYZ(1,IP))
			xmin = 	min(xmin, XYZ(1,IP))
			ymax = 	max(ymax, XYZ(2,IP))
			ymin = 	min(ymin, XYZ(2,IP))
			zmax = 	max(zmax, XYZ(3,IP))
			zmin = 	min(zmin, XYZ(3,IP))

			CTR(1) = CTR(1) + XYZ(1,IP)
			CTR(2) = CTR(2) + XYZ(2,IP)
			CTR(3) = CTR(3) + XYZ(3,IP)
390		CONTINUE
		CTR(1) = CTR(1)/WElem(1,IW)
		CTR(2) = CTR(2)/WElem(1,IW)
		CTR(3) = CTR(3)/WElem(1,IW)

		Gx = ABS(Gap*WNV(1,IW))
		Gy = ABS(Gap*WNV(2,IW))
		Gz = ABS(Gap*WNV(3,IW))

		xmax = 	xmax + (dx + Gx)*FDA
		xmin = 	xmin - (dx + Gx)*FDA
		ymax = 	ymax + (dy + Gy)*FDA
		ymin = 	ymin - (dy + Gy)*FDA
		zmax = 	zmax + (dz + Gz)*FDA
		zmin = 	zmin - (dz + Gz)*FDA
	
		ixmax = (xmax-MeshRange(1,1))/dx+1
		ixmin = (xmin-MeshRange(1,1))/dx+1
		iymax = (ymax-MeshRange(1,2))/dy+1
		iymin = (ymin-MeshRange(1,2))/dy+1
		izmax = (zmax-MeshRange(1,3))/dz+1
		izmin = (zmin-MeshRange(1,3))/dz+1
		
		ictr = (CTR(1)-MeshRange(1,1))/dx+1
		jctr = (CTR(2)-MeshRange(1,2))/dy+1
		kctr = (CTR(3)-MeshRange(1,3))/dz+1

		ixmax = min(MB,ixmax)
!		ixmax = min(ictr+idstep,ixmax)
		ixmin = max(1,ixmin)
!		ixmin = max(ictr-idstep,ixmin)
		iymax = min(MB,iymax)
!		iymax = min(jctr+idstep,iymax)
		iymin = max(1,iymin)
!		iymin = max(jctr-idstep,iymin)
		izmax = min(MB,izmax)
!		izmax = min(kctr+idstep,izmax)
		izmin = max(1,izmin)
!		izmin = max(kctr-idstep,izmin)
		
		num = (ixmax-ixmin)*(iymax-iymin)*(izmax-izmin)
	
		DO 395 iz = izmin, izmax
		 DO 395 iy = iymin, iymax
		  DO 395 ix = ixmin, ixmax
			WBlockNum(ix,iy,iz) = WBlockNum(ix,iy,iz) + 1
395		CONTINUE

400	CONTINUE

	

	
	!put WBlockIndex
	NPB = 0
	maxnum = WBlockNum(1,1,1)
	DO 450 K = 1, MB
		DO 450 J = 1, MB
			DO 450 I = 1, MB
			WBlockIndex(i,j,k) = NPB
			NPB = NPB + WBlockNum(i,j,k)
			maxnum = max(WBlockNum(i,j,k), maxnum)
			WBlockNum(i,j,k) = 0
		
450	CONTINUE
	

	WRITE(IUT0,*) '  maximum number in a Surface bucket box = ', 
     &			  maxnum
	!WRITE(IUT0,*) '  minimum number in a block = ', minnum
	WRITE(IUT0,*) '  Stacked', NPB, ' times in total'


	MPB=NPB
	ALLOCATE(WBlockID(NPB),stat=ierr)
	IF(ierr.NE.0) STOP 'allocate array in StackWall()'


	IF(NPB.GT.10000000) THEN
		WRITE(IUT0,*) '     !!! WARNING: in StackWall():'
		WRITE(IUT0,*) '     Too many stacks', NPB
	END IF


	DO 490 IW = 1, NW
		IB = WELEM(WElem(1,IW)+2,IW)
		IF(WSMK(IB).EQ.0) CYCLE

		xmax=MeshRange(1,1)
		xmin=MeshRange(2,1)
		ymax=MeshRange(1,2)
		ymin=MeshRange(2,2)
		zmax=MeshRange(1,3)
		zmin=MeshRange(2,3)
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0

		DO 460 J = 1, WElem(1,IW)
			IP = WElem(J+1,IW)
			xmax = 	max(xmax, XYZ(1,IP))
			xmin = 	min(xmin, XYZ(1,IP))
			ymax = 	max(ymax, XYZ(2,IP))
			ymin = 	min(ymin, XYZ(2,IP))
			zmax = 	max(zmax, XYZ(3,IP))
			zmin = 	min(zmin, XYZ(3,IP))

			CTR(1) = CTR(1) + XYZ(1,IP)
			CTR(2) = CTR(2) + XYZ(2,IP)
			CTR(3) = CTR(3) + XYZ(3,IP)
460		CONTINUE
		CTR(1) = CTR(1)/WElem(1,IW)
		CTR(2) = CTR(2)/WElem(1,IW)
		CTR(3) = CTR(3)/WElem(1,IW)

		Gx = ABS(Gap*WNV(1,IW))
		Gy = ABS(Gap*WNV(2,IW))
		Gz = ABS(Gap*WNV(3,IW))

		xmax = 	xmax + (dx + Gx)*FDA
		xmin = 	xmin - (dx + Gx)*FDA
		ymax = 	ymax + (dy + Gy)*FDA
		ymin = 	ymin - (dy + Gy)*FDA
		zmax = 	zmax + (dz + Gz)*FDA
		zmin = 	zmin - (dz + Gz)*FDA
	
		ixmax = (xmax-MeshRange(1,1))/dx+1
		ixmin = (xmin-MeshRange(1,1))/dx+1
		iymax = (ymax-MeshRange(1,2))/dy+1
		iymin = (ymin-MeshRange(1,2))/dy+1
		izmax = (zmax-MeshRange(1,3))/dz+1
		izmin = (zmin-MeshRange(1,3))/dz+1
		
		ictr = (CTR(1)-MeshRange(1,1))/dx+1
		jctr = (CTR(2)-MeshRange(1,2))/dy+1
		kctr = (CTR(3)-MeshRange(1,3))/dz+1

		ixmax = min(MB,ixmax)
!		ixmax = min(ictr+idstep,ixmax)
		ixmin = max(1,ixmin)
!		ixmin = max(ictr-idstep,ixmin)
		iymax = min(MB,iymax)
!		iymax = min(jctr+idstep,iymax)
		iymin = max(1,iymin)
!		iymin = max(jctr-idstep,iymin)
		izmax = min(MB,izmax)
!		izmax = min(kctr+idstep,izmax)
		izmin = max(1,izmin)
!		izmin = max(kctr-idstep,izmin)

		DO 465 iz = izmin, izmax
		 DO 465 iy = iymin, iymax
		  DO 465 ix = ixmin, ixmax

			WBlockNum(ix,iy,iz) = WBlockNum(ix,iy,iz) + 1
		num = WBlockNum(ix,iy,iz) + WBlockIndex(ix,iy,iz)
			WBlockID(num) = IW

465		CONTINUE


490	CONTINUE


!------> end of StackWall()

	END 
C======================================================================C
	SUBROUTINE StackPatch(IUT0,XYZ,WallElem,MP,N2D,NDOF,NP,NE,NW,
     &				 BlockNodeID,WBlockNum,WBlockIndex,
     &			MeshRange,MPB,NPB,MB,FDA,ROUT,RIN,GlobeCTR)
C======================================================================C     	

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER(PAI=3.141592653589793)
C	
	REAL*8  XYZ(NDOF,MP),MeshRange(2,NDOF),RIN(MB,MB),ROUT
	INTEGER WallElem(N2D+1,NW)
	INTEGER BlockNodeID(MPB), WBlockNum(MB,MB),
     &		WBlockIndex(MB,MB)
	CHARACTER*10  StackMode
	REAL*8	CTR(3),GlobeCTR(3),vs(N2D+1)
	
	

10	deta  = PAI/MB
	dseta = 2.0*PAI/MB
	
	
	DO 310 I =  1, MB
		DO 310 J =  1, MB
			WBlockNum(i,j) = 0
			WBlockIndex(i,j) = 0
310	CONTINUE	
	
	GlobeCTR = 0.0
	DO 100 I = 1, NP
	  DO 100 J = 1,NDOF
		GlobeCTR(J)=GlobeCTR(J)+XYZ(J,I)
100	CONTINUE
	GlobeCTR = GlobeCTR/NP


	rr=(MeshRange(2,1)
     & -MeshRange(1,1))*(MeshRange(2,1)-MeshRange(1,1))
     & +(MeshRange(2,2)-MeshRange(1,2))*(MeshRange(2,2)-MeshRange(1,2))
     & +(MeshRange(2,3)-MeshRange(1,3))*(MeshRange(2,3)-MeshRange(1,3))
	radius = sqrt(rr)

	ROUT = 0.0
	RIN  = radius*1.0E10
	RIN_MIN = radius*1.0E10

	deta = PAI/MB
	dseta = 2*PAI/MB
		
	!count 
	DO 400 IW = 1, NW
		

		etamax = -PAI
		etamin = 2.0*PAI
		setamax = -5.0*PAI
		setamin = 5.0*PAI
		
		DO 320 J = 1, WallElem(1,IW)
			IP = WallElem(J+1,IW)
	radius =  (XYZ(1,IP)-GlobeCTR(1))*(XYZ(1,IP)-GlobeCTR(1))
     &		 +(XYZ(2,IP)-GlobeCTR(2))*(XYZ(2,IP)-GlobeCTR(2))
     &		 +(XYZ(3,IP)-GlobeCTR(3))*(XYZ(3,IP)-GlobeCTR(3))
			radius = SQRT(radius)
			ROUT = MAX(ROUT,radius)
			RIN_MIN = MIN(RIN_MIN,radius)
			
	CALL GetSolidAngle(seta,eta,XYZ(1,IP),XYZ(2,IP),XYZ(3,IP),
     &			   GlobeCTR(1),GlobeCTR(2),GlobeCTR(3))
		!seta :	circumferential angle
		!eta  : cone angle (zenith angle) 
			
			vs(j) = seta
			setamin =  min(setamin,seta)
			setamax =  max(setamax,seta)
			etamax  =  max(etamax, eta)
			etamin  =  min(etamin, eta)
			

320		CONTINUE
		
		IF(setamax-setamin.GT.PAI) THEN
					
			setamax = -5.0*PAI
			setamin = 5.0*PAI
			DO 325 J =1,WallElem(1,IW)
				if(vs(j).GT.PAI) vs(j) = vs(j)-2*PAI
				setamin =  min(setamin,vs(J))
				setamax =  max(setamax,vs(J))
325			CONTINUE
		END IF
		
		IF(setamax.LT.0) setamax = setamax+2*PAI
		IF(setamin.LT.0) setamin = setamin+2*PAI
		etamax = etamax + 0.25*deta
		setamax = setamax + 0.25*dseta
		etamin = etamin - 0.25*deta
		setamin = setamin - 0.25*dseta


c		IF((setamax-setamin).GT.PAI) THEN
c			a = setamax
c			setamax = setamin+2.0*PAI
c			setamin = a
c		ELSE IF(setamin.LT.0) THEN
c			setamax = setamax+2.0*PAI
c			setamin = setamin+2.0*PAI
c		END IF

		ietamax = etamax/deta+1
		ietamin = etamin/deta+1
		ietamax = min(MB,ietamax)
		ietamin = max(1,ietamin)

		
		isetamax = setamax/dseta+1
		isetamin = setamin/dseta+1
		
		IF(isetamax.LT.isetamin) THEN
			
			DO 330 iseta = isetamin, MB
				DO 330 ieta = ietamin, ietamax

  		WBlockNum(ieta,iseta) = WBlockNum(ieta,iseta) + 1
		RIN(ieta,iseta) = min(RIN(ieta,iseta),radius)

330			CONTINUE			
			DO 340 iseta = 1, isetamax
				DO 340 ieta = ietamin, ietamax

  		WBlockNum(ieta,iseta) = WBlockNum(ieta,iseta) + 1
		RIN(ieta,iseta) = min(RIN(ieta,iseta),radius)

340			CONTINUE

		ELSE
			isetamax = min(MB,isetamax)
			isetamin = max(1,isetamin)

		DO 380 iseta = isetamin, isetamax
		DO 380 ieta = ietamin, ietamax
  		WBlockNum(ieta,iseta) = WBlockNum(ieta,iseta) + 1
		RIN(ieta,iseta) = min(RIN(ieta,iseta),radius)
 380		CONTINUE

		END IF 
		

400	CONTINUE


	!
	ROUT = ROUT*1.05
	RIN = RIN*0.95
	RIN_MIN = RIN_MIN*0.95

	
	!put blockindex

	NPB = 0
	maxnum = WBlockNum(1,1)
	DO 450 I = 1, MB
		DO 450 J = 1, MB
			WBlockIndex(i,j) = NPB
			NPB = NPB + WBlockNum(i,j)
			maxnum = max(WBlockNum(i,j), maxnum)
			
			WBlockNum(i,j) = 0
			
450	CONTINUE
	
	WRITE(IUT0,*) ' '
	WRITE(IUT0,*) 
! --- 
     & '  maximum number in a WALL bucket box = ', maxnum
	WRITE(IUT0,*) '  Stacked', NPB, ' times in total'
	WRITE(IUT0,*) '  Outer radias of Globe = ', ROUT
	WRITE(IUT0,*) '  Minimum Inner radias of Globe = ', RIN_MIN

	IF(NPB.GT.MPB) THEN
		WRITE(IUT0,*) 'ERROR in StackPatch():'
		WRITE(IUT0,*) 'Please enlarge MPB'
		STOP
	END IF



	!count 
	DO 500 IW = 1, NW
		

		etamax = -PAI
		etamin = 2.0*PAI
		setamax = -5.0*PAI
		setamin = 5.0*PAI

		DO 420 J = 1, WallElem(1,IW)
			IP = WallElem(J+1,IW)
						
	CALL GetSolidAngle(seta,eta,XYZ(1,IP),XYZ(2,IP),XYZ(3,IP),
     &			   GlobeCTR(1),GlobeCTR(2),GlobeCTR(3))
	!seta :	circumferential angle
	!eta  : cone angle
			
			vs(j) = seta
			setamin =  min(setamin,seta)
			setamax =  max(setamax,seta)
			etamin  =  min(etamin, eta)
			etamax  =  max(etamax, eta)

420		CONTINUE

		IF(setamax-setamin.GT.PAI) THEN

			setamax = -5.0*PAI
			setamin = 5.0*PAI
		
			DO 425 J =1,WallElem(1,IW)
				if(vs(j).GT.PAI) vs(j) = vs(j)-2.0*PAI
				setamin =  min(setamin,vs(J))
				setamax =  max(setamax,vs(J))
425			CONTINUE
		

		END IF
		
		IF(setamax.LT.0) setamax = setamax+2*PAI
		IF(setamin.LT.0) setamin = setamin+2*PAI
		etamax = etamax + 0.25*deta
		setamax = setamax + 0.25*dseta
		etamin = etamin - 0.25*deta
		setamin = setamin - 0.25*dseta


c		IF((setamax-setamin).GT.PAI) THEN
c			a = setamax
c			setamax = setamin+2.0*PAI
c			setamin = a
c		ELSE IF(setamin.LT.0) THEN
c			setamax = setamax+2.0*PAI
c			setamin = setamin+2.0*PAI
c		END IF

		ietamax = etamax/deta+1
		ietamin = etamin/deta+1
		ietamax = min(MB,ietamax)
		ietamin = max(1,ietamin)

		
		isetamax = setamax/dseta+1
		isetamin = setamin/dseta+1
		
		IF(isetamax.LT.isetamin) THEN
			
			DO 430 iy = isetamin, MB
				DO 430 ix = ietamin, ietamax
			WBlockNum(ix,iy) = WBlockNum(ix,iy) + 1
			num = WBlockNum(ix,iy) + WBlockIndex(ix,iy)
			BlockNodeId(num) = IW

430			CONTINUE			
			DO 440 iy = 1, isetamax
				DO 440 ix = ietamin, ietamax

	  		WBlockNum(ix,iy) = WBlockNum(ix,iy) + 1
			num = WBlockNum(ix,iy) + WBlockIndex(ix,iy)
			BlockNodeId(num) = IW

440			CONTINUE



		ELSE
			isetamax = min(MB,isetamax)
			isetamin = max(1,isetamin)
			DO 480 iy = isetamin, isetamax
			DO 480 ix = ietamin, ietamax
	  		WBlockNum(ix,iy) = WBlockNum(ix,iy) + 1
			num = WBlockNum(ix,iy) + WBlockIndex(ix,iy)
			BlockNodeId(num) = IW
480			CONTINUE

		END IF 
		
		
500	CONTINUE


	DO 600 I=1,MB
	  DO 600 J = 1,MB
		IF(WBlockNum(I,J).EQ.0) RIN(I,J) = ROUT
		IF(WBlockNum(I,J).GT.100) THEN

		WRITE(*,*) I,J,WBlockNum(I,J),i*deta,j*dseta

	END IF
600	CONTINUE

	RETURN

!------>

	END
C======================================================================C
	SUBROUTINE StackVolumeGroup(IUT0,CTR,ELEM,ME,N3D,NDOF,NE,NRD,
     &			 NNPE,NKIND,ElemID,BlockElemID,BlockElemNum,
     &			 BlockIndex,MeshRange,MPB,NPB,MB)
C======================================================================C     	
	IMPLICIT REAL*8(A-H,O-Z)
	!PARAMETER( MaxLimit=2000, idstep=1 )
C	
	REAL*8  CTR(NDOF,NRD),MeshRange(2,NDOF)
	INTEGER ELEM(N3D,ME),NNPE(20),ElemID(NE)
	INTEGER BlockElemID(MPB), BlockElemNum(MB,MB,MB),
     &		BlockIndex(MB,MB,MB)
	CHARACTER*10  StackMode
	
	
	

10	dx = (MeshRange(2,1)-MeshRange(1,1))/MB
	dy = (MeshRange(2,2)-MeshRange(1,2))/MB
	dz = (MeshRange(2,3)-MeshRange(1,3))/MB
	maxnum = 0
	minnum = 1.0E5

	BlockElemNum(:,:,:) = 0
	BlockIndex(:,:,:) = 0

	
	
	!write(*,*) 'NE=',NE,MB
	
	!count  
	DO 400 IE = 1, NE
		
		
		ix = (CTR(1,IE)-MeshRange(1,1))/dx+1
		iy = (CTR(2,IE)-MeshRange(1,2))/dy+1
		iz = (CTR(3,IE)-MeshRange(1,3))/dz+1

		ix = min(MB,ix)
		ix = max(1,ix)
		iy = min(MB,iy)
		iy = max(1,iy)
		iz = min(MB,iz)
		iz = max(1,iz)
		
		BlockElemNum(ix,iy,iz) = BlockElemNum(ix,iy,iz) + 1

400	CONTINUE


	
	!put blockindex
	NPB = 0
	maxnum = BlockElemNum(1,1,1)
	DO 450 K = 1, MB
		DO 450 J = 1, MB
			DO 450 I = 1, MB
			BlockIndex(i,j,k) = NPB
			NPB = NPB + BlockElemNum(i,j,k)
			maxnum = max(BlockElemNum(i,j,k), maxnum)
			BlockElemNum(i,j,k) = 0
		
450	CONTINUE
	
	WRITE(IUT0,*) '  Gas Bucket For Grouping ...'
c	WRITE(IUT0,*) '  maximum number in a GAS bucket box = ', maxnum
c	WRITE(IUT0,*) '  Stacked', NPB, ' times in total'


	DO 490 IE = 1, NE
		ix = (CTR(1,IE)-MeshRange(1,1))/dx+1
		iy = (CTR(2,IE)-MeshRange(1,2))/dy+1
		iz = (CTR(3,IE)-MeshRange(1,3))/dz+1

		ix = min(MB,ix)
		ix = max(1,ix)
		iy = min(MB,iy)
		iy = max(1,iy)
		iz = min(MB,iz)
		iz = max(1,iz)
		
		BlockElemNum(ix,iy,iz) = BlockElemNum(ix,iy,iz) + 1
		num = BlockElemNum(ix,iy,iz) + BlockIndex(ix,iy,iz)
		BlockElemId(num) = IE

490	CONTINUE


!------> end of stackvolumegroup()

	END 
C======================================================================C
	SUBROUTINE StackWallGroup(IUT0,CTR,WElem,NW,NE,NRD,N2D,N3D,
     &			 NDOF,BlockWallID,BlockWallNum,BlockIndex,
     &			 MeshRange,MPB,NPB,MB)
C======================================================================C     	
	IMPLICIT REAL*8(A-H,O-Z)
!PARAMETER( MaxLimit=2000, idstep=1 )
C	
	REAL*8  CTR(NDOF,NRD),MeshRange(2,NDOF),WNV(NDOF,NW)
	INTEGER WElem(N2D+2,NW)
	INTEGER BlockWallID(MPB), BlockWallNum(MB,MB,MB),
     &		BlockIndex(MB,MB,MB)
	CHARACTER*10  StackMode
	
	
	
10	dx = (MeshRange(2,1)-MeshRange(1,1))/MB
	dy = (MeshRange(2,2)-MeshRange(1,2))/MB
	dz = (MeshRange(2,3)-MeshRange(1,3))/MB
	Gap = MAX(MAX(dx,dy),dz)*0.25
	
	BlockWallNum(:,:,:) = 0
	BlockIndex(:,:,:) = 0

	!count 
	DO 400 IW = 1, NW
		
		IIT=NE+IW
		ix = (CTR(1,IIT)-MeshRange(1,1))/dx+1
		iy = (CTR(2,IIT)-MeshRange(1,2))/dy+1
		iz = (CTR(3,IIT)-MeshRange(1,3))/dz+1

		ix = min(MB,ix)
		ix = max(1,ix)
		iy = min(MB,iy)
		iy = max(1,iy)
		iz = min(MB,iz)
		iz = max(1,iz)
		BlockWallNum(ix,iy,iz) = BlockWallNum(ix,iy,iz) + 1

400	CONTINUE


	
	!put blockindex
	NPB = 0
	maxnum = BlockWallNum(1,1,1)
	DO 450 K = 1, MB
		DO 450 J = 1, MB
			DO 450 I = 1, MB
			BlockIndex(i,j,k) = NPB
			NPB = NPB + BlockWallNum(i,j,k)
			maxnum = max(BlockWallNum(i,j,k), maxnum)
			BlockWallNum(i,j,k) = 0
		
450	CONTINUE
	
	WRITE(IUT0,*) '  Surface Bucket For Grouping ...'
c	WRITE(IUT0,*) '  maximum number in a Surface bucket box = ', 
c     &			     maxnum
c	WRITE(IUT0,*) '  Stacked', NPB, ' times in total'


	DO 490 IW = 1, NW

		IIT=NE+IW
		ix = (CTR(1,IIT)-MeshRange(1,1))/dx+1
		iy = (CTR(2,IIT)-MeshRange(1,2))/dy+1
		iz = (CTR(3,IIT)-MeshRange(1,3))/dz+1

		ix = min(MB,ix)
		ix = max(1,ix)
		iy = min(MB,iy)
		iy = max(1,iy)
		iz = min(MB,iz)
		iz = max(1,iz)

		BlockWallNum(ix,iy,iz) = BlockWallNum(ix,iy,iz) + 1
		num = BlockWallNum(ix,iy,iz) + BlockIndex(ix,iy,iz)
		BlockWallId(num) = IW

490	CONTINUE


!------> end of stackwallgroup()

	END 
C======================================================================C
	SUBROUTINE StackWallPrd(IUT0,XYZ,WElem,MP,NP,NW,N2D,N3D,NDOF,
     &			MeshRange,WNV,NPB,MB,FDA,NBOUN,IBOUN,
     &					prdData)
C======================================================================C
	USE MODULE_RADSXF, ONLY : WBlockNum,WBlockIndex,WBlockID


	IMPLICIT REAL*8(A-H,O-Z)
	!PARAMETER( MaxLimit=2000, idstep=1 )
C	
	REAL*8  XYZ(NDOF,MP),MeshRange(2,NDOF),WNV(NDOF,NW)
	REAL*8  prdData(NBOUN,5)
	INTEGER WElem(N2D+2,NW)
	CHARACTER*10  StackMode
	REAL*8	CTR(3),Gap
	
	
10	dx = (MeshRange(2,1)-MeshRange(1,1))/MB
	dy = (MeshRange(2,2)-MeshRange(1,2))/MB
	dz = (MeshRange(2,3)-MeshRange(1,3))/MB
	Gap = MAX(MAX(dx,dy),dz)*1.2
	maxnum = 0
	minnum = 1.0E5
	

	ALLOCATE(WBlockNum(MB,MB,MB),WBlockIndex(MB,MB,MB))

	WBlockNum(:,:,:) = 0
	WBlockIndex(:,:,:) = 0

	
	
	!count 
	DO 400 IW = 1, NW
		
		NK =  WElem(WElem(1,IW)+2,IW)
		
		IF(NK.NE.IBOUN) CYCLE

		xmax=MeshRange(1,1)
		xmin=MeshRange(2,1)
		ymax=MeshRange(1,2)
		ymin=MeshRange(2,2)
		zmax=MeshRange(1,3)
		zmin=MeshRange(2,3)
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0

		DO 390 J = 1, WElem(1,IW)
			IP = WElem(J+1,IW)
			xmax = 	max(xmax, XYZ(1,IP))
			xmin = 	min(xmin, XYZ(1,IP))
			ymax = 	max(ymax, XYZ(2,IP))
			ymin = 	min(ymin, XYZ(2,IP))
			zmax = 	max(zmax, XYZ(3,IP))
			zmin = 	min(zmin, XYZ(3,IP))

			CTR(1) = CTR(1) + XYZ(1,IP)
			CTR(2) = CTR(2) + XYZ(2,IP)
			CTR(3) = CTR(3) + XYZ(3,IP)
390		CONTINUE
		CTR(1) = CTR(1)/WElem(1,IW)
		CTR(2) = CTR(2)/WElem(1,IW)
		CTR(3) = CTR(3)/WElem(1,IW)

		Gx = ABS(Gap*WNV(1,IW))
		Gy = ABS(Gap*WNV(2,IW))
		Gz = ABS(Gap*WNV(3,IW))

		xmax = 	xmax + (dx + Gx)*FDA
		xmin = 	xmin - (dx + Gx)*FDA
		ymax = 	ymax + (dy + Gy)*FDA
		ymin = 	ymin - (dy + Gy)*FDA
		zmax = 	zmax + (dz + Gz)*FDA
		zmin = 	zmin - (dz + Gz)*FDA
	
		ixmax = (xmax-MeshRange(1,1))/dx+1
		ixmin = (xmin-MeshRange(1,1))/dx+1
		iymax = (ymax-MeshRange(1,2))/dy+1
		iymin = (ymin-MeshRange(1,2))/dy+1
		izmax = (zmax-MeshRange(1,3))/dz+1
		izmin = (zmin-MeshRange(1,3))/dz+1
		
		ictr = (CTR(1)-MeshRange(1,1))/dx+1
		jctr = (CTR(2)-MeshRange(1,2))/dy+1
		kctr = (CTR(3)-MeshRange(1,3))/dz+1

		ixmax = min(MB,ixmax)
!		ixmax = min(ictr+idstep,ixmax)
		ixmin = max(1,ixmin)
!		ixmin = max(ictr-idstep,ixmin)
		iymax = min(MB,iymax)
!		iymax = min(jctr+idstep,iymax)
		iymin = max(1,iymin)
!		iymin = max(jctr-idstep,iymin)
		izmax = min(MB,izmax)
!		izmax = min(kctr+idstep,izmax)
		izmin = max(1,izmin)
!		izmin = max(kctr-idstep,izmin)
		
		num = (ixmax-ixmin)*(iymax-iymin)*(izmax-izmin)
	
		DO 395 iz = izmin, izmax
		 DO 395 iy = iymin, iymax
		  DO 395 ix = ixmin, ixmax
			WBlockNum(ix,iy,iz) = WBlockNum(ix,iy,iz) + 1
395		CONTINUE

400	CONTINUE

	
!put WBlockIndex
	NPB = 0
	maxnum = WBlockNum(1,1,1)
	DO 450 K = 1, MB
		DO 450 J = 1, MB
			DO 450 I = 1, MB
			WBlockIndex(i,j,k) = NPB
			NPB = NPB + WBlockNum(i,j,k)
			maxnum = max(WBlockNum(i,j,k), maxnum)
			WBlockNum(i,j,k) = 0
450	CONTINUE
	MPB=NPB
	ALLOCATE(WBlockID(NPB),stat=ierr)
	IF(ierr.NE.0) STOP 'allocate array in StackWallPrd()'


	IF(NPB.GT.10000000) THEN
		WRITE(IUT0,*) '     !!! WARNING: in StackWallPrd():'
		WRITE(IUT0,*) '     Too many stacks', NPB
	END IF


	DO 490 IW = 1, NW
		NK =  WElem(WElem(1,IW)+2,IW)
		IF(NK.NE.IBOUN) CYCLE
		xmax=MeshRange(1,1)
		xmin=MeshRange(2,1)
		ymax=MeshRange(1,2)
		ymin=MeshRange(2,2)
		zmax=MeshRange(1,3)
		zmin=MeshRange(2,3)
		CTR(1) = 0
		CTR(2) = 0
		CTR(3) = 0

		DO 460 J = 1, WElem(1,IW)
			IP = WElem(J+1,IW)
			xmax = 	max(xmax, XYZ(1,IP))
			xmin = 	min(xmin, XYZ(1,IP))
			ymax = 	max(ymax, XYZ(2,IP))
			ymin = 	min(ymin, XYZ(2,IP))
			zmax = 	max(zmax, XYZ(3,IP))
			zmin = 	min(zmin, XYZ(3,IP))

			CTR(1) = CTR(1) + XYZ(1,IP)
			CTR(2) = CTR(2) + XYZ(2,IP)
			CTR(3) = CTR(3) + XYZ(3,IP)
460		CONTINUE
		CTR(1) = CTR(1)/WElem(1,IW)
		CTR(2) = CTR(2)/WElem(1,IW)
		CTR(3) = CTR(3)/WElem(1,IW)

		Gx = ABS(Gap*WNV(1,IW))
		Gy = ABS(Gap*WNV(2,IW))
		Gz = ABS(Gap*WNV(3,IW))

		xmax = 	xmax + (dx + Gx)*FDA
		xmin = 	xmin - (dx + Gx)*FDA
		ymax = 	ymax + (dy + Gy)*FDA
		ymin = 	ymin - (dy + Gy)*FDA
		zmax = 	zmax + (dz + Gz)*FDA
		zmin = 	zmin - (dz + Gz)*FDA
	
		ixmax = (xmax-MeshRange(1,1))/dx+1
		ixmin = (xmin-MeshRange(1,1))/dx+1
		iymax = (ymax-MeshRange(1,2))/dy+1
		iymin = (ymin-MeshRange(1,2))/dy+1
		izmax = (zmax-MeshRange(1,3))/dz+1
		izmin = (zmin-MeshRange(1,3))/dz+1
		
		ictr = (CTR(1)-MeshRange(1,1))/dx+1
		jctr = (CTR(2)-MeshRange(1,2))/dy+1
		kctr = (CTR(3)-MeshRange(1,3))/dz+1

		ixmax = min(MB,ixmax)
!		ixmax = min(ictr+idstep,ixmax)
		ixmin = max(1,ixmin)
!		ixmin = max(ictr-idstep,ixmin)
		iymax = min(MB,iymax)
!		iymax = min(jctr+idstep,iymax)
		iymin = max(1,iymin)
!		iymin = max(jctr-idstep,iymin)
		izmax = min(MB,izmax)
!		izmax = min(kctr+idstep,izmax)
		izmin = max(1,izmin)
!		izmin = max(kctr-idstep,izmin)

		DO 465 iz = izmin, izmax
		 DO 465 iy = iymin, iymax
		  DO 465 ix = ixmin, ixmax
		WBlockNum(ix,iy,iz) = WBlockNum(ix,iy,iz) + 1
		num = WBlockNum(ix,iy,iz) + WBlockIndex(ix,iy,iz)
		WBlockID(num) = IW

465		CONTINUE


490	CONTINUE


!------> end of StackWallPrd()

	END 

