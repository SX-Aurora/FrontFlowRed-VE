
!
!	Subroutine List 
!	surbroutine  getproperty()
!	surbroutine  readrestart()
!	surbroutine  writerestart()
!	surbroutine  outputview()
!	(part of subroutines for Equation Solvers)
!	subroutine MakeEnergyEqu()
!	subroutine MakeEnergyEquGp()
!	subroutine T4Solver()
!	subroutine SORLTRI()
!	subroutine ICCGLTRI()
!	subroutine GSLTRI2()
!	subroutine SorBand()
!	subroutine ICCGBand()
!	subroutine GSBand()
!	subroutine GetHeatFlux()
!	subroutine RDVElemToNode()
!	subroutine GetNodeIndex()
!	subroutine NodeToGPProp
!	subroutine UserDefRadPro
!-----------------------------------------------------------------------------
C======================================================================C
C                                                                      C
C SUBROUTINE GetProperty()			 		       C
C                                                                      C
!		2005/07/08~					       C
C                                                                      C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C			                                               C
C======================================================================C
	SUBROUTINE GetProperty(IUT0,MODEL,NE,NW,WRadProp,WType,
     &	 MatID,NMAT,NWT,NBOUN,N2D,NRD0,NRD,WGrayDeg,scabeta,
     &	 GAbsorb,PAbsorb,PScatter,WElem,VIndex,Prop,WArea,Volume,
     &	 Albedo,CTR,AveVs,HomFlag,EType,userflag)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
      use module_material, only:  radmat,radfludflag

!radmat(:,flud)	 1:	 gabsorb	
!				 2:	 pabsorb	
!				 3:	 pscatter	
!				 4:	 scabeta
!				 5:	 ptcdia		(no use heat)
!radfludflag(:,flud)	 1:	 gastype
!				 2:	 scatype
!				 3:	 ptc cal. flag
!				 4:	 fgpflag
!						
	IMPLICIT REAL*8(A-H,O-Z)
	INCLUDE 'pre_common.f'
!	
!	
	INTEGER WElem(N2D+2,NWT),VIndex(NRD0),MatID(NE)
	INTEGER WType(NBOUN),HomFlag,EType(NRD0),userflag
	REAL*8 GAbsorb(NE),PAbsorb(NE),PScatter(NE),WGrayDeg(NWT)
	REAL*8 WArea(NWT),Volume(NE),PROP(NRD0),Albedo(NRD0)
	REAL*8 WRadProp(NBOUN),CTR(3,NRD0),AveVs(10),scabeta(NE)

!avevs:		1,		ave_extincoef
!			2,		ave_wgraydeg
!			3,		ave_gabsorb
!			4,		ave_pabsorb
!			5,		ave_pscatter
	CHARACTER*10    MODEL
!		
	GAbsorb = 0.
	PAbsorb = 0.
	PScatter = 0.
	WGrayDeg = 0.
	AveVs = 0
	ETYPE = 1
	gmax=0.0
	gmin=0.0
	Mark=0
	
	IF(userflag.eq.1) THEN
          DO IE = 1,NE
	  ID=MatID(IE)
!IF(ID.LT.1.OR.ID.GT.2) WRITE(*,*) IE,ID
            IF(ID.LE.0.OR.ID.GT.NMAT) CYCLE

		CALL UserDefRadPro(CTR(1,IE),CTR(2,IE),CTR(3,IE),
     &			GAbsorb(IE),PAbsorb(IE),PScatter(IE))
	
		EType(IE) = RadFludFlag(2,ID)
		IF(ETYPE(IE).GT.1.AND.PScatter(IE).GT.0) Mark=1
!RadFludFlag(2,:)	scattering type
!EType =1,	Linear-anisotropic scatter / isotropic scatter, OK
!	2,  large-diffuse particle scatter, OK
!	3,	Rayleigh scatter, OK
!	4,	Mie scatter
!	5,  user defined scatter
!-----------------------------------------		
          scabeta(IE)=radmat(4,ID)
	  
	  END DO
	ELSE
	  DO IE = 1,NE
		ID=MatID(IE)
		IF(ID.LE.0.OR.ID.GT.NMAT) CYCLE

		IF(RadFludFlag(3,ID).EQ.1.AND.PScatter(IE).GT.0) THEN
		 WRITE(IUT0,*) 'Error Abort in RadSxf() '
		 WRITE(IUT0,*) 'Because you set <ptccal = 1>'
		 WRITE(IUT0,*) 'Please Set <radwxflag=1>'
	 WRITE(IUT0,*) 'and edit the subroutine UserDefRadPro()'
	 WRITE(IUT0,*) 'to calculate the scattering parameters'
		 STOP
		END IF
		
		IF(radfludflag(1,ID).eq.2) THEN
		 WRITE(IUT0,*) 'Error Abort in RadSxf() '
		 WRITE(IUT0,*) 'Because you use REAL_GAS'
	 WRITE(IUT0,*) 'Please edit the subroutine UserDefRadPro()'
		 WRITE(IUT0,*) 'to calculate the gas parameters '
		 WRITE(IUT0,*) 'or set radmodel=<FVM>'
		 STOP
		END IF


		GAbsorb(IE) = RadMat(1,ID)
		PAbsorb(IE) = RadMat(2,ID)
		PScatter(IE) = RadMat(3,ID)
		EType(IE) = RadFludFlag(2,ID)
		IF(ETYPE(IE).GT.1.AND.PScatter(IE).GT.0) Mark=1
!
!radmat(4,:)	scattering type
!EType =	1,	Linear-anisotropic scatter / isotropic scatter, OK
!		2,  large-diffuse particle scatter, OK
!		3,	Rayleigh scatter, OK
!		4,	Mie scatter
!		5,  user defined scatter


		scabeta(IE)=radmat(4,ID)
	
	  END DO
	END IF
	
	AveVs(3)=SUM(GAbsorb(1:NE))/NE
	AveVs(4)=SUM(PAbsorb(1:NE))/NE
	AveVs(5)=SUM(PScatter(1:NE))/NE
	AveVs(1)=SUM(AveVs(3:5))

	
	
	DO IW = 1,NWT
		IB = WElem(1,IW)+2
		IB = WElem(IB,IW)
		IB = ABS(IB)		!for period pair

		WGrayDeg(IW) = WRadProp(IB)			
!NOTE:
!period boundary, WRadProp(IB)=-1.0, see pre_radiation.f

		ETYPE(IW+NE) = WType(IB)
		IF(ETYPE(IW+NE).GT.1)	Mark=1
		IF(WGrayDeg(IW).LT.1.d0) Mark=1
     
!WTYPE:		1,	diffuse-wall	,OK 
!		2,  mirror-wall		,OK
!		3,	directional-diffuse wall	,OK
!		4,  user defined

	END DO
	AveVs(2)=SUM(WGrayDeg(1:NWT))/NWT


!error abort
	IF(Mark.EQ.1.AND.Model(1:4).EQ.'ZONE') THEN
	 WRITE(IUT0,*) 'Error Abort in RadSxf() '
	 WRITE(IUT0,*) 'Either the wall_reflect or the gas_scatter is '
	 WRITE(IUT0,*) 'anisotropic, so that <ZONE> can not be used'
	 WRITE(IUT0,*) 'Please change model to <MC> or <FVM>'
	 STOP
	END IF
	
	
	
C4-----------------------------------------------------------------------
	VIndex = 0
	IRD = 0
!	
	IF(MODEL(1:4).EQ.'ZONE') THEN
          DO 100 I= 1,NE
          aa = GAbsorb(I)+PAbsorb(I)+PScatter(I)
	  IF(aa.GT.0.0) THEN
	    IRD = IRD + 1
	    VIndex(I) = IRD
	    PROP(I) = 4.0*aa*Volume(I)
	    Albedo(I)=0
	  ENDIF
100	  CONTINUE
!
          DO 200 I= 1,NW
          IB = WElem(1,I)+2
	  IB = WElem(IB,I)
	  IF(IB.LE.0) CYCLE		
!for period/interface/slide boundry
!only host boundary (ib>0) is participating
	  IRD = IRD + 1
	  VIndex(NE+I) = IRD
	  PROP(NE+I) = WArea(I)
	  Albedo(NE+I)=0.0
200	  CONTINUE
	ELSE
	  Albedo(:)=1.0d0
	  gmax=GAbsorb(1)+PAbsorb(1)
	  gmin=gmax
	  DO 110 I= 1,NE
	  aa = (GAbsorb(I)+PAbsorb(I)+PScatter(I))
	  gmax=max(aa,gmax)
	  gmin=min(aa,gmin)
		  
	  IF(aa.GT.0.0) THEN
	    IRD = IRD + 1
	    VIndex(I) = IRD
	    PROP(I) = 4.0*(GAbsorb(I)+PAbsorb(I))*Volume(I)
	    Albedo(I)=PScatter(I)/aa
	  END IF
110	  CONTINUE
	
	  DO 210 I= 1,NW
	  IB = WElem(1,I)+2
	  IB = WElem(IB,I)
	  IF(IB.LE.0) CYCLE
!for period/interface/slide boundry
! only host boundary (ib>0) is participating
	  IF(WGrayDeg(I).GT.0.0) THEN
	    IRD = IRD + 1
	    VIndex(NE+I) = IRD
	    PROP(NE+I) = WArea(I)*WGrayDeg(I)
	    Albedo(NE+I)=1.0-WGrayDeg(I)
	  ENDIF
210	CONTINUE
	ENDIF
!
	NRD=IRD
!	
	IF(gmax.EQ.gmin) THEN
	  HomFlag = 1
        ELSE
	  HomFlag = 0
	END IF
		
	IF(NRD.EQ.0) THEN
	WRITE(IUT0,*) 'Error in GetProperty():'
	WRITE(IUT0,*) 'No element could absorb ray because'
	WRITE(IUT0,*) 'both gas and wall have zero absorption coef.'
	STOP
	END IF
	

	RETURN

2000	WRITE(*,*) 'ERROR ABORT in GetPorperty()'
	WRITE(*,*) FileProp, 'is not existant'
	STOP

!-------> end of GetProperty()

	END

C======================================================================C
	SUBROUTINE ReadRestart(IUT2,NRD0,MRD0,NRD,MRD,StackOpt)
C======================================================================C
	USE MODULE_RADSXF, ONLY : RDVALUE,RDIndex,RDN1,RDN2,RDID
	USE MODULE_RADSXF, ONLY : PROP

	IMPLICIT REAL*8(A-H,O-Z)
	INTEGER NRD0,MRD0,NRD,MRD
	CHARACTER*20 StackOpt

	OPEN(IUT2,FILE='rdtmp.bin',		!CONVERT='BIG_ENDIAN',
     &	 form='unformatted',status='unknown',iostat=ios)
		read(IUT2) StackOpt
		read(IUT2) NRD,MRD
		read(IUT2) (RDIndex(I),I=1,NRD+1)
		read(IUT2) (RDValue(I),I=1,MRD)
		IF(StackOpt(1:4).EQ.'BAND') THEN
			read(IUT2) (RDN1(I),I=1,NRD)
			read(IUT2) (RDN2(I),I=1,NRD)
			read(IUT2) (RDId(I),I=1,MRD)
		END IF
		CLOSE(IUT2)


	RETURN

!--------->

	END 

C======================================================================C
	SUBROUTINE WriteRestart(IUT2,NRD0,MRD0,NRD,MRD,StackOpt)
C======================================================================C		
	USE MODULE_RADSXF, ONLY : RDVALUE,RDIndex,RDN1,RDN2,RDID
	USE MODULE_RADSXF, ONLY : PROP

	IMPLICIT REAL*8(A-H,O-Z)
	INTEGER NRD0,MRD0,NRD,MRD
	CHARACTER*20 StackOpt

		OPEN(IUT2,FILE='rdtmp.bin',		!CONVERT='BIG_ENDIAN',
     &	 form='unformatted',status='unknown',iostat=ios)
		write(IUT2) StackOpt
		write(IUT2) NRD,MRD
		write(IUT2) (RDIndex(I),I=1,NRD+1)
		write(IUT2) (RDValue(I),I=1,MRD)
		IF(StackOpt(1:4).EQ.'BAND') THEN
			write(IUT2) (RDN1(I),I=1,NRD)
			write(IUT2) (RDN2(I),I=1,NRD)
			write(IUT2) (RDId(I),I=1,MRD)
		END IF
		CLOSE(IUT2)


	RETURN

!--------->

	END 

C======================================================================C
	SUBROUTINE OutputView(XYZ,ELEM,WElem,WNEB,
     &		    FaceNum,NNPE,NNPS,CNTIVITY,ElemKind,
     &			NRD0,NRD,NE,ME,NW,NWT,NP,MP,NDOF,N2D,N3D,
     &			NKIND,NFACE)
C======================================================================C
	
	IMPLICIT REAL*8(A-H,O-Z)
	
	REAL*8 	XYZ(NDOF,MP)
	INTEGER ELEM(N3D,ME),ElemKind(NE)
	INTEGER WElem(N2D+2,NWT),WNEB(NW)
	CHARACTER*80 FileRes
	CHARACTER*20 Scheme
	INTEGER FaceNum(NKIND),NNPE(NKIND),NNPS(NKIND,NFACE),
     &		CNTIVITY(NKIND,NFACE,N2D)

	


C---------------------------- mesh grid data ---------------------------------------
	OPEN(1,FILE='View-Mesh.dat')

	!
	WRITE(1,'(I8)') NP
	DO 20 I = 1, NP
		WRITE(1,'(I18, 3G16.5)') I,XYZ(1,I),
     &			XYZ(2,I),XYZ(3,I)
20	CONTINUE
	
	WRITE(1,*) NKIND
	DO 30 I = 1,NKIND
		WRITE(1,'(2I8)') NNPE(I),FaceNum(I)
30	CONTINUE
	DO 40 I = 1,NKIND
		IF(FaceNum(I).GT.0) THEN
		DO 35 J = 1,FaceNum(I)
		WRITE(1,'(5I8)') NNPS(I,J),
     &		(CNTIVITY(I,J,K),K=1,NNPS(I,J))
35		CONTINUE
		END IF
40	CONTINUE
		
	WRITE(1,'(I8)') NE
	DO 60 I = 1, NE
		ID = ElemKind(I)		  
		WRITE(1,'(10I8)') ID,I,
     &		  (Elem(J,I),J=1,NNPE(ID))
60	CONTINUE


	!surface element		
	WRITE(1,'(I8)') NW
	DO 300 I = 1, NW
	WRITE(1,'(10I8)') (WElem(J,I),J=1,WElem(1,I)+2), WNEB(I)
		!welem:
		!(1,:) = n		:	num of vertex, n
		!(2~n+1,:)		:	vertex no.
		!(n+2,:)		:	Wall Group ID

300	CONTINUE
	

C---------------------------- result data ---------------------------------------

	CLOSE(1)


	RETURN

!---------> end of OutputView()
	END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE MakeEnergyEqu(VMIN,StackOpt,NRD0,NRD,MRD,MID)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	USE MODULE_RADSXF, ONLY :  RDId,RDN1,RDN2,RDIndex,RDValue
	USE MODULE_RADSXF, ONLY :  VIndex,PROP

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)
	CHARACTER*20 StackOpt
	

C----------------------------------------------------------------------------C
C	Make A * X = B
C	Energy conservation equation:
C	PiiXi = sum(Dij*Xj) + Bi	-----
C			
C		
C	Aii*Xi + sum(Aij*Xj) = Bi	<----
C	Aij = -Dij		(for i!=j)
C	Aii = Pii-Dii
C	Bi = HSource
!
	IF(StackOpt(1:4).EQ.'LTRI') THEN
!		
!Calculate Bi

	DO 400 II=1,NRD0

!BB(I) = HSource(I)
	  I=VIndex(II)
	  IF(I.LE.0) CYCLE
	  ID = RDIndex(I) 		
	  RDValue(ID+I) = PROP(II)-RDValue(ID+I) !Aii
!Calculate Aij
          DO 300 J = 1,I-1
!unknown element
          RDValue(ID+J) = - RDValue(ID+J)
!---------------ATTENTION!!!-------------------!
!it is negative for Aij(i!=j)
!----------------------------------------------!
 300	  CONTINUE
		     			      
400	CONTINUE
		
	ELSE
!Calculate Bi
          DO 1400 II=1,NRD0
            I=VIndex(II)
	    IF(I.LE.0) CYCLE
              ID = RDIndex(I)
!BB(I) = HSource(I)
              RDValue(ID+1)=PROP(II)-RDValue(ID+1) !Aii
!Aij
            DO 1400 J = RDIndex(I)+2,RDIndex(I+1)
		RDValue(J) = -RDValue(J)
						
1400	    CONTINUE
		
		!Normalize the equation with Aii
		DO 1600 I = 1,NRD
			ID = RDIndex(I)
			Aii = RDValue(ID+1)
			DO J = ID+1, RDIndex(I+1)
				RDValue(J) = RDValue(J)/Aii
			END DO
			!BB(I) = BB(I)/Aii
1600		CONTINUE
		
		!move big coefs to front
		!this is a special treatment for ICCG SOLVER
		RDN1 = 0		!left big coefs
		RDN2 = 0		!total big coefs, RDN2-RDN1 = right big coefs
		DO 1700 I = 1,NRD
			ISTA = RDIndex(I) + 1
			DO 1700 J = RDIndex(I)+2, RDIndex(I+1)
				
			  IF(RDValue(J).GT.VMIN) THEN
				!big coef
				RDN2(I) = RDN2(I) + 1
				ID1 = ISTA + RDN2(I)
				VA = RDValue(ID1)
				IA = RDId(ID1)
				RDValue(ID1) = RDValue(J)
				RDId(ID1) = RDId(J)
				RDValue(J) = VA
				RDId(J) = IA

				IF(RDId(ID1).LT.I) THEN		
					!left off-diagnal elements
					RDN1(I) = RDN1(I) + 1
					IDL = ISTA + RDN1(I)
					VB = RDValue(IDL)
					IB = RDId(IDL)
					RDValue(IDL) = RDValue(ID1)
					RDId(IDL) = RDId(ID1)
					RDValue(ID1)=VB
					RDId(ID1)=IB
				END IF
			  END IF
		
1700		CONTINUE


	END IF

	
	RETURN

!-------> end of MakeEnergyEqu()

	END



C======================================================================C
	SUBROUTINE MakeEnergyEquGP(IUT0,NRD0,NRD,MRD)
C======================================================================C

	USE MODULE_RADGROUP, ONLY : GpProp
	USE MODULE_RADSXF, ONLY : RDValue,RDIndex

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)

	
		
C----------------------------------------------------------------------------C
C	Make A * X = B
C	Energy conservation equation:
C			PiiXi = sum(Dij*Xj) + Bi	-----
C											 |
C											 |
C			Aii*Xi + sum(Aij*Xj) = Bi	<----
C			Aij = -Dij		(for i!=j)
C			Aii = Pii-Dii
C			Bi = HSource



	!Calculate Bi
	DO 400 I=1,NRD
	
		!BB(I) = HSource(I)
		ID = RDIndex(I) 		
  		RDValue(ID+I) = GpPROP(I)-RDValue(ID+I)		!Aii
	
		!Calculate Aij
		DO 300 J = 1,I-1
			!unknown element
			
			RDValue(ID+J) = - RDValue(ID+J)
!---------------ATTENTION!!!-------------------!
!it is negative for Aij(i!=j)
!----------------------------------------------!
300		CONTINUE

400	CONTINUE
	
	
	RETURN

!-------> end of MakeEnergyEquGP()

	END
C======================================================================C
	SUBROUTINE MakeEnergyEqu2(VIndex,RDValue,
     &            RDIndex,RDN1,RDN2,RDId,
     &		 TT,T,BB,PROP,WK,VMIN,
     &		 StackOpt,NRD0,NRD,MRD,MID)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)
	INTEGER RDId(MID),RDN1(NRD),RDN2(NRD),RDIndex(NRD+1)
	INTEGER VIndex(NRD0),WK(NRD)
	REAL*8  RDValue(MRD),PROP(NRD)
	REAL*8  TT(NRD),T(NRD0),BB(NRD)
	CHARACTER*20 StackOpt
	

	!Make A * T = B

	
	WK = 0
	IV = 0
	DO 100 I =1,NRD0
		ID = VIndex(I)
		IF(ID.GT.0) THEN
			IF(T(I).LE.0) THEN
				IV = IV+1
				WK(ID) = IV
			END	IF
			TT(ID) = SHIGMA*T(I)**4
		END IF
100	CONTINUE
	IVRD = IV		!!!!!------->


	BB = 0
	IF(StackOpt(1:4).EQ.'LTRI') THEN
		
		!Calculate Bi
		DO 400 I=1,NRD

		  I1 = WK(I)
		  IF(I1.GT.0) THEN
		
				!unknown element
			ID = RDIndex(I)
			DO 200 J = 1,I-1
				J1 = WK(J)
				IF(J1.EQ.0) THEN
					!known element, stack to BB
					energy = RDValue(ID+J)*TT(J)
					BB(I1) = BB(I1) + energy
				END IF
200			CONTINUE
			
			
			DO 210 J = I+1,NRD
				J1 = WK(J)
				IF(J1.EQ.0) THEN
					!known element, stack to BB
					ID1 = RDIndex(J)
					energy = RDValue(ID1+I)*TT(J)
					BB(I1) = BB(I1) + energy
				END IF
210			CONTINUE
			
		  END IF
400		CONTINUE
		
		!Calculate Aij
		IMRD = 0
		DO 500 I =1,NRD
		  I1 = WK(I)

		  IF(I1.GT.0) THEN
			!unknown elements stack it to matrix
			
			ID = RDIndex(I)
			RDIndex(I1) = IMRD
			DO 300 J = 1,I-1
				J1 = WK(J)
				IF(J1.GT.0) THEN
					!unknown element
				RDValue(IMRD+J1) = - RDValue(ID+J)
!---------------ATTENTION!!!-------------------!
!it is negative for Aij(i!=j)
!----------------------------------------------!
				END IF
300			CONTINUE
			
			RDValue(IMRD+I1) = PROP(I)-RDValue(ID+I)
			IMRD = IMRD + I1
			
		  END IF
				      
500		CONTINUE

		RDIndex(IVRD+1) = IMRD
	
!---------------ATTENTION!!!-------------------!
!DO NOT TRY TO NORMALIZE THE EQUATION, BECAUSE 
!WE ONLY HAVE LOWER-HALF OF THE MATRIX
!----------------------------------------------!

	ELSE

		!Calculate Bi
		DO 1400 I=1,NRD
		  I1 = WK(I)
		  IF(I1.GT.0) THEN
				!unknown element
			DO 1200 J = RDIndex(I)+2,RDIndex(I+1)
				ID = RDId(J)
				J1 = WK(ID)
				IF(J1.EQ.0) THEN
					!known element, stack to BB
					energy = RDValue(J)*TT(ID)
					BB(I1) = BB(I1) + energy
				END IF
1200			CONTINUE
						
		  END IF
1400		CONTINUE
		
		!Calculate Aij
		IMRD = 0
		DO 1500 I =1,NRD
		  I1 = WK(I)
		  IF(I1.GT.0) THEN
			!unknown elements stack it to matrix
			ID = RDIndex(I)
			RDIndex(I1) = IMRD
			
			IMRD = IMRD + 1
			RDValue(IMRD) = PROP(I)-RDValue(ID+1)
			RDId(IMRD) = I1
						
			DO 1300 J = ID+2,RDIndex(I+1)
				J1 = WK(J)
				IF(J1.GT.0) THEN
					!unknown element
					IMRD = IMRD + 1
					RDValue(IMRD) = - RDValue(J)
					RDId(IMRD) = J1
!---------------ATTENTION!!!-------------------!
!it is negative for Aij(i!=j)
!----------------------------------------------!
				END IF
1300			CONTINUE
							
		  END IF
		  
1500		CONTINUE

		RDIndex(IVRD+1) = IMRD

		!Normalize the equation with Aii
		DO 1600 I = 1,IVRD
			ID = RDIndex(I)
			Aii = RDValue(ID+1)
			DO J = ID+1, RDIndex(I+1)
				RDValue(J) = RDValue(J)/Aii
			END DO
			BB(I) = BB(I)/Aii
1600		CONTINUE

		!move big coefs to front
		!this is a special treatment for ICCG SOLVER
		RDN1 = 0		!left big coefs
		RDN2 = 0		!total big coefs, RDN2-RDN1 = right big coefs
		DO 1700 I = 1,IVRD
			ISTA = RDIndex(I) + 1
			DO 1700 J = RDIndex(I)+2, RDIndex(I+1)
				
			  IF(RDValue(J).GT.VMIN) THEN
				!big coef
				RDN2(I) = RDN2(I) + 1
				ID1 = ISTA + RDN2(I)
				VA = RDValue(ID1)
				IA = RDId(ID1)
				RDValue(ID1) = RDValue(J)
				RDId(ID1) = RDId(J)
				RDValue(J) = VA
				RDId(J) = IA
  
				IF(RDId(ID1).LT.I) THEN		
					!left off-diagnal elements
					RDN1(I) = RDN1(I) + 1
					IDL = ISTA + RDN1(I)
					VB = RDValue(IDL)
					IB = RDId(IDL)
					RDValue(IDL) = RDValue(ID1)
					RDId(IDL) = RDId(ID1)
					RDValue(ID1)=VB
					RDId(ID1)=IB
				END IF
			  END IF
		
1700		CONTINUE


	END IF


	!New index of T -> TT
	DO 2000 I = 1,NRD0
		ID = VIndex(I)
		IF(ID.GT.0) THEN
			IV = WK(ID)
			VIndex(I) = IV
		END IF

2000	CONTINUE
									!		|
	NRD = IVRD						!<------|
	MRD = IMRD
	RDIndex(NRD+1) = MRD

	RETURN

!-------> end of makeenergyequ2()

	END

C======================================================================C
	SUBROUTINE T4Solver(IUT0,NRD0,MRD,NRD,MID,NVRD,T,TT,BB,VWK,
     &			RDValue,RDN1,RDN2,RDIndex,RDId,VIndex,RevIndex,
     &			aeps,reps,RLX,AlgorithmOpt,StackOpt)
C======================================================================C	

		
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)

	INTEGER RDId(MID),RDN1(NRD),RDN2(NRD),RDIndex(NRD+1)
	INTEGER VIndex(NRD0),RevIndex(NVRD)
	REAL*8  RDValue(MRD)
	REAL*8  TT(NRD),T(NRD0),BB(NRD),VWK(NRD)
	
	CHARACTER*20 AlgorithmOpt,StackOpt

	itrmax = 300

!aeps	- the absolute error tolerance for the convergence of equation solution
!reps	- the relative error tolerance for the convergence of equation solution
!RLX	- relaxtion factor
	
	
	
	IF(StackOpt.EQ.'LTRI')	THEN

		IF(AlgorithmOpt(1:3).EQ.'SOR') THEN

		CALL SORLTRI(IUT0,RDIndex,RDValue,TT,BB,RevIndex,
     &				 aeps,reps,RLX,NRD,MRD,NVRD,itr)
	
		ELSE IF(AlgorithmOpt(1:2).EQ.'CG') THEN
			
		CALL ICCGLTRI(IUT0,NRD,MRD,NVRD,RevIndex,
     &		RDIndex,RDValue,BB,TT,aeps,reps,itrmax,itr,ier)
		ELSE
			
	CALL GauseLTRI(NVRD,NRD,MRD,RDValue,BB,TT,RDIndex,RevIndex)
		END IF


	ELSE

		IF(AlgorithmOpt(1:3).EQ.'SOR') THEN

		CALL SORBAND(IUT0,RDIndex,RDValue,TT,BB,RDId,RevIndex,
     &				 aeps,reps,RLX,NRD,MRD,MID,NVRD,itr)

		ELSE IF(AlgorithmOpt(1:2).EQ.'CG') THEN
			
	CALL BICGSTAB_BAND(IUT0,NRD,MRD,NVRD,RevIndex,RDId,RDN1,
     &		RDN2,RDIndex,RDValue,BB,TT,aeps,reps,itrmax,itr,ier)

		ELSE
			
		CALL GauseBAND(NVRD,NRD,MRD,RDValue,RDId,BB,TT,RDIndex,
     &				   RevIndex)

		END IF


	END IF
	

	!Transfer values to T
	DO 100 I =1,NRD0
		ID = VIndex(I)
		IF(ID.GT.0) THEN
			TT(ID) = (TT(ID)/SHIGMA)**0.25
			T(I) = TT(ID)
		END IF
			
100	CONTINUE



	RETURN

!-------> T4Solver()

	END


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine BICGSTAB_BAND(IUT0,NRD,MRD,NVRD,RevIndex,
     2		   RDId,IL,IQ,RDIndex,RDValue,bb,TT,aeps,reps,itrmax,
     3		   itr,ier)
	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. Bi-CGstab solver with diagonal scaling
!     & incomplete LU decomposition
!
!     nn   : no. of elements
!     IEMAX   : maximul no. of off-diagonal elements
!     IQ   : no. of off-diagonal elements
!     IL   : no. of off-diagonal elements at left part (lower part)
!     ip   : address of off-diagonal element
!     aa   : coefficient matrix
!     bb   : right hand side vector (in)
!          : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence (out)
!     itr  : maximul iteration counts (in)
!          : iteration counts at convergence (out)
!     ier  : return code
!
! --- [module arguments]
!
      
      use module_radslv,only : d  =>W1K1
      use module_radslv,only : p0 =>W1K2
      use module_radslv,only : p1 =>W1K3
      use module_radslv,only : q0 =>W1K4
      use module_radslv,only : q1 =>W1K5
      use module_radslv,only : r0 =>W1K6
      use module_radslv,only : r1 =>W1K7
      use module_radslv,only : r2 =>W1K8
      use module_radslv,only : x  =>W1K9
      use module_radslv,only : z  =>W1K10
      use module_radslv,only : diag_a  =>W1K11
	use module_radslv,only : a0  =>W1K12
	

	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (MXMAT=1,NMAT=1)
!
c      implicit none
!
! --- [dummy arguments]
!	
      real*8 ,intent(inout) :: bb(NRD)
      real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
	INTEGER	RDId(MRD), IL(NRD),IQ(NRD),RDIndex(NRD),RevIndex(NVRD)
	REAL*8	RDValue(MRD), RDLeft(NRD),RDSelf(NRD)

!
! --- [local entities]
!
      real*8 ,parameter :: XHALF=0.50D0,SML=1.D-25,ZERO=0.D0
      real*8 :: c1,t(0:MXMAT),s(0:MXMAT),u(0:MXMAT),v(0:MXMAT),
     &          rmaxs,rmaxs0
      real*8 :: alpha(MXMAT),beta(MXMAT),omega(MXMAT),dum1
      integer:: i,j,k,l,it1,it2,jer,m1,ICV,ierr,KMAT
      integer:: IMAT,IIMAT,ICVS,ICVE,ICVL

	DATA ifle / 0 / 
!
!

	
	ALLOCATE (d(NRD),p0(NRD),p1(NRD),q0(NRD),q1(NRD),r0(NRD),
     &		  r1(NRD),r2(NRD),x(NRD),z(NRD),diag_a(NRD),a0(NRD),
     &		  stat=ierr)

	if(ierr.ne.0) then
        write(ifle,*) 'allocating array error in utl_bcgstb'
	  stop
      endif

      t(:)=0.d0
      s(:)=0.d0
      u(:)=0.d0
      v(:)=0.d0
      alpha(:)=0.d0
      beta(:)=0.d0
      omega(:)=0.d0
      ierr=0

	!-----------------------------------!
	!	Jiang 2005/07/04
	!-----------------------------------!
	ICVS = 1
	ICVE = NRD
	IIMAT = 1
	!a0(:) <---- aa(:,0)
	!diag_a(:) <---- a(:,0)
	!RDValue(:) <---- a(:,j!=0)

	a0 = 1.0




!
!-< 1. Preliminary set >-
!
      rmaxs=1.d0
      DO 100 IIMAT=1,NMAT    !ICV=1,NCV
c      if(.not.mat_cal(IIMAT)) goto 100
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
      DO ICVL=ICVS,ICVE
      rmaxs=min(rmaxs,abs(a0(ICVL)))
      END DO
  100 CONTINUE
!
c      IF(NPE.gt.1) THEN
c        CALL hpcrmin(RMAXS)
c      ENDIF
!
c      if(RMAXS.lt.0.d0) then
c        WRITE(IUT0,*) '### interface error -1- (utl_bcgstb)',rmaxs
c        CALL ABORT(1,'utl_bcgstb')
c      elseif(RMAXS.eq.0.d0) then
c        deallocate(r1,r2,x,z,a)
c        return
c      endif
!
      ier=0
      jer=0
      it1=1
      it2=itr
      itr=0
!
      d (:)=0.d0
      p0(:)=0.d0
      p1(:)=0.d0
      q0(:)=0.d0
      q1(:)=0.d0
      r0(:)=0.d0
      r1(:)=0.d0
      r2(:)=0.d0
      x (:)=0.d0
!
!-< 2. Preconditioning >-
!
!
!--< 2.1 diagonal scaling >--
!
      DO 2000 IIMAT=1,NMAT    !
c      if(.not.mat_cal(IIMAT)) goto 2000
c      ICVS=MAT_INTER(IIMAT-1)+1
c      ICVE=MAT_INTER(IIMAT)
      DO 200 ICVL=ICVS,ICVE   !ICV=1,NCV
      d(ICVL)=1.d0/(sqrt(abs(a0(ICVL)))+SML)
      z(ICVL)=d(ICVL)
      diag_a(ICVL)=sign(1.d0,a0(ICVL))
  200 CONTINUE
 2000 CONTINUE


!
! --- Low and Upper region
!

	!------------------------------------------------!
	!we set d(:) = diag_a(:) = 1
	!it is not necessary to change RDValue(:) 
	!					---- Jiang 2005/07/04
	!------------------------------------------------!

c      rmaxs=0.d0
c      DO 2500 IIMAT=1,NMAT    !
c      if(.not.mat_cal(IIMAT)) goto 2500
c      ICVS=MAT_INTER(IIMAT-1)+1
c      ICVE=MAT_INTER(IIMAT)
c      DO 210 ICVL=ICVS,ICVE   !ICV=1,NCV
c      DO 210 j=1,IQ(ICVL)
c      RDValue(JJ)=aRDValue(JJ)*d(ICVL)*d(RDId(JJ))*diag_a(ICVL)
c  210 CONTINUE
!
      DO 220 ICVL=ICVS,ICVE    !ICV=1,NCV
      bb(ICVL)=bb(ICVL)*d(ICVL)*diag_a(ICVL)
      x(ICVL)=0.d0
      rmaxs=max(rmaxs,abs(bb(ICVL)))
  220 CONTINUE

!
c      ICVS=MAT_INTER(IIMAT-1)+1
c      ICVE=MAT_INTER(IIMAT)
      d(ICVS:ICVE)=1.d0
      diag_a(ICVS:ICVE)=0.d0
 2500 CONTINUE
!
c      IF(NPE.gt.1) THEN
c        CALL hpcrmax(RMAXS)
c      ENDIF
!
      if(RMAXS.le.aeps) goto 999
!
      rmaxs0=RMAXS
!      
! --- all below by NCVIN 
!
      DO 3000 IIMAT=1,NMAT    !
c      if(.not.mat_cal(IIMAT)) goto 3000
c      ICVS=MAT_INTER(IIMAT-1)+1
c      ICVE=MAT_INTER(IIMAT)
!
!--< 2.2 incomplete LU decomposition >--
!
      dum1=0.d0
      DO 230 ICVL=ICVS,ICVE    !ICV=1,NCVIN
      if(d(ICVL).le.0.d0) then
        dum1=d(ICVL)
      endif
      d(ICVL)=1.d0/(d(ICVL)+SML)
! --- upper region
      DO 231 j=IL(ICVL)+1,IQ(ICVL)
	JJ = j + RDIndex(ICVL)
      k=RDId(JJ)
      c1=0.d0
! --- low
      DO 232 l=1,IL(k)
	KK = l + RDIndex(k)
      if(RDId(KK).eq.ICVL) c1=d(ICVL)*RDValue(KK)
  232 CONTINUE
      d(k)=d(k)-c1*RDValue(JJ)
  231 CONTINUE
  230 CONTINUE
!
 3000 CONTINUE
!
c      if(NPE.gt.1) then
c        call hpcrmin(dum1)
c      endif
!
      if(dum1.lt.0.d0) then
        write(*,*) '### CPU NUMBER MY_RANK= ',MY_RANK
        write(*,*) '### interface error -2- (mc_bcgstb)',dum1
c        CALL ABORT(1,'utl_bcgstb')
      endif
!
! --- 
!
c      IF(NPE.gt.1) then
c        CALL SOLVER_SEND_RECV(1,MXCV,NCV,D)
c      endif
!
!-< 3. Initial set for iteration >-
!
  888 CONTINUE
!
c      IF(NPE.gt.1) then
c        CALL SOLVER_SEND_RECV(1,MXCV,NCV,X)
c      endif
!
! --- 
!
      rmaxs=0.d0
      DO 4000 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 4000
      DO 300 ICVL=ICVS,ICVE     !ICV=1,NCVIN
      dum1=x(ICVL)
! --- upper and low region:
      DO 301 j=1,IQ(ICVL)
	JJ = j + RDIndex(ICVL)
      dum1=dum1+RDValue(JJ)*x(RDId(JJ))
  301 CONTINUE
!
      r0(ICVL)=bb(ICVL)-dum1
      rmaxs=max(rmaxs,abs(r0(ICVL)))
  300 CONTINUE
 4000 CONTINUE
!
c      IF(NPE.gt.1) THEN
c        CALL hpcrmax(RMAXS)
c      ENDIF
!
      if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
!
      s(:)=0.d0
      DO 5000 IIMAT=1,NMAT
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 5000
      DO 310 ICVL=ICVS,ICVE     !i=1,NCVIN
        p0(ICVL)=r0(ICVL)
        r1(ICVL)=r0(ICVL)
        dum1=sign(1.d0,r0(ICVL))*max(1.d0,abs(r0(ICVL)))
        s(IIMAT)=s(IIMAT)+dum1*r0(ICVL)
        r0(ICVL)=dum1
  310 CONTINUE
 5000 CONTINUE
!
c      if(NPE.gt.1) then
c        DO KMAT=1,KMAT_S
c          IIMAT=MAT_S(KMAT)
c          CALL hpcrsum(s(IIMAT))
c        END DO
c      endif
      c1=0.d0
      DO IIMAT=1,NMAT
        c1=c1+abs(s(IIMAT))
      END DO
c      if(NPE.gt.1) then
c        CALL hpcrsum(c1)
c      endif
!
      if(c1.lt.1.d-20) jer=1
!
      if( jer.gt.0 ) then
        DO 6000 IIMAT=1,NMAT
c        ICVS=MAT_INDEX(IIMAT-1)+1
c        ICVE=MAT_INDEX(IIMAT)
c        if(.not.mat_cal(IIMAT)) goto 6000
        DO 320 ICVL=ICVS,ICVE     !i=1,NCVIN
          dum1=0.d0
! --- low region:
          DO 321 j=1,IL(ICVL)
		  JJ = j + RDIndex(ICVL)
            dum1=dum1+RDValue(JJ)*p1(RDId(JJ))
  321     CONTINUE
          p1(ICVL)=d(ICVL)*(r1(ICVL)-dum1)
  320   CONTINUE
 6000   CONTINUE
!
c        IF(NPE.gt.1) then
c          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
c        endif
!
        DO 6100 IIMAT=1,NMAT
c        ICVS=MAT_INDEX(IIMAT-1)+1
c        ICVE=MAT_INDEX(IIMAT)
c        if(.not.mat_cal(IIMAT)) goto 6100
        DO 322 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper region:
        DO 323 j=IL(ICVL)+1,IQ(ICVL)
	  JJ = j + RDIndex(ICVL)
        dum1=dum1+RDValue(JJ)*p1(RDId(JJ))
  323   CONTINUE
        p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
        x(ICVL)=x(ICVL)+p1(ICVL)
  322   CONTINUE
 6100   CONTINUE
!
!
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif
!
!-< 4. Iteration >-
!
! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
      DO 1000 itr=it1,it2
!
c        IF(NPE.gt.1) then
c          CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
c        endif
!
      DO 7000 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7000
      DO 400 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=0.d0
! --- low
      DO 401 j=1,IL(ICVL)
	JJ = j + RDIndex(ICVL)
      dum1=dum1+RDValue(JJ)*p1(RDId(JJ))
  401 CONTINUE
      p1(ICVL)=d(ICVL)*(p0(ICVL)-dum1)
  400 CONTINUE
 7000 CONTINUE
!
c      IF(NPE.gt.1) then
c        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
c      endif
!
! --- 
!
      DO 7100 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7100
      DO 402 ICVL=ICVE,ICVS,-1    !i=NCVIN,1,-1
      dum1=0.d0
! --- upper
      DO 403 j=IL(ICVL)+1,IQ(ICVL)
	JJ = j + RDIndex(ICVL)
      dum1=dum1+RDValue(JJ)*p1(RDId(JJ))
  403 CONTINUE
      p1(ICVL)=p1(ICVL)-d(ICVL)*dum1
  402 CONTINUE
 7100 CONTINUE
!
c      IF(NPE.gt.1) then
c        CALL SOLVER_SEND_RECV(1,MXCV,NCV,P1)
c      endif
!
      v(:)=0.d0
      DO 7200 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7200
      DO 410 ICVL=ICVS,ICVE    !i=1,NCVIN
      dum1=p1(ICVL)
! --- low + upper
      DO 411 j=1,IQ(ICVL)
	JJ = j + RDIndex(ICVL)
      dum1=dum1+RDValue(JJ)*p1(RDId(JJ))
  411 CONTINUE
      q1(ICVL)=dum1
      v(IIMAT)=v(IIMAT)+r0(ICVL)*dum1
  410 CONTINUE
 7200 CONTINUE
!
c      if(NPE.gt.1) then
c        DO KMAT=1,KMAT_S
c          IIMAT=MAT_S(KMAT)
c          CALL hpcrsum(v(IIMAT))
c        END DO
c      endif      
!
      c1=0.0d0
      DO IIMAT=1,NMAT
      c1=c1+abs(v(IIMAT))
      END DO
c      IF(NPE.gt.1) THEN
c        CALL hpcrsum(c1)
c      ENDIF
      if(c1.lt.1.d-20) then
        jer=1
        it1=itr
        goto 888
      endif
!
      DO IIMAT=1,NMAT
      alpha(IIMAT)=s(IIMAT)/(v(IIMAT)+SML) 
      END DO
!
c      IF(NPE.gt.1) then
c        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
c      endif
!
      DO 7300 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7300
      DO 420 ICVL=ICVS,ICVE    !i=1,NCVIN
        q0(ICVL)=r1(ICVL)-alpha(IIMAT)*q1(ICVL)
  420 CONTINUE
!
      DO 430 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=0.d0
! --- low
        DO 431 j=1,IL(ICVL)
		JJ = j + RDIndex(ICVL)
          dum1=dum1+RDValue(JJ)*r1(RDId(JJ))
  431   CONTINUE
        r1(ICVL)=d(ICVL)*(q0(ICVL)-dum1)
  430 CONTINUE
 7300 CONTINUE
!
!
c      IF(NPE.gt.1) then
c        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
c      endif
!
!
      DO 7400 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7400
      DO 432 ICVL=ICVE,ICVS,-1   !i=NCVIN,1,-1
        dum1=0.d0
! --- upper
        DO 433 j=IL(ICVL)+1,IQ(ICVL)
		JJ = j + RDIndex(ICVL)
          dum1=dum1+RDValue(JJ)*r1(RDId(JJ))
  433   CONTINUE
        r1(ICVL)=r1(ICVL)-d(ICVL)*dum1
  432 CONTINUE
 7400 CONTINUE
!
c      IF(NPE.gt.1) then
c        CALL SOLVER_SEND_RECV(1,MXCV,NCV,R1)
c      endif
!
      u(:)=0.d0
      v(:)=0.d0
      DO 7500 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7500
      DO 440 ICVL=ICVS,ICVE    !i=1,NCVIN
        dum1=r1(ICVL)
        DO 441 j=1,IQ(ICVL)
! --- low+upper
	  JJ = j + RDIndex(ICVL)
        dum1=dum1+RDValue(JJ)*r1(RDId(JJ))
  441   CONTINUE
        r2(ICVL)=dum1
        u(IIMAT)=u(IIMAT)+dum1*q0(ICVL)
        v(IIMAT)=v(IIMAT)+dum1*dum1
  440 CONTINUE
 7500 CONTINUE
!
c      if(NPE.gt.1) then
c        DO KMAT=1,KMAT_S
c          IIMAT=MAT_S(KMAT)
c          CALL hpcrsum(u(IIMAT))
c          CALL hpcrsum(v(IIMAT))
c        END DO
c      endif
      c1=0.d0
      DO IIMAT=1,NMAT
        c1=c1+abs(v(IIMAT))
      END DO
c      if(NPE.gt.1) then
c        CALL hpcrsum(c1)
c      endif
!
      if(c1.lt.1.d-20) then
        jer=1
        it1=itr
        goto 888
      endif
!
      DO IIMAT=1,NMAT
      omega(IIMAT)=u(IIMAT)/(v(IIMAT)+SML)
      END DO
!
      rmaxs=0.d0
      DO 7600 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7600
      DO 450 ICVL=ICVS,ICVE   !i=1,NCVIN
        x(ICVL)=x(ICVL)+(alpha(IIMAT)*p1(ICVL)+omega(IIMAT)*r1(ICVL))
        r1(ICVL)=q0(ICVL)-omega(IIMAT)*r2(ICVL)
        rmaxs=max(rmaxs,abs(r1(ICVL)))
  450 CONTINUE
 7600 CONTINUE
!
c      IF(NPE.gt.1) THEN
c        CALL hpcrmax(RMAXS)
c      ENDIF
!
      if(RMAXS.le.max(aeps,reps*rmaxs0)) goto 999
!
      t(:)=0.d0
      DO 7700 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7700
      DO 460 ICVL=ICVS,ICVE   !i=1,NCVIN
        t(IIMAT)=t(IIMAT)+r0(ICVL)*r1(ICVL)
  460 CONTINUE
 7700 CONTINUE
!
c      IF(NPE.gt.1) THEN
c        DO KMAT=1,KMAT_S
c          IIMAT=MAT_S(KMAT)
c          CALL hpcrsum(t(IIMAT))
c        END DO
c      ENDIF
!
      c1=0.d0
      DO IIMAT=1,NMAT
      c1=c1+abs(s(IIMAT)*omega(IIMAT))
      END DO
c      IF(NPE.gt.1) THEN
c        CALL hpcrsum(c1)
c      ENDIF
      if(c1.lt.1.d-20 ) then
        jer=1
        it1=itr
        goto 888
      endif
!
      DO IIMAT=1,NMAT
      beta(IIMAT)=t(IIMAT)*alpha(IIMAT)/(s(IIMAT)*omega(IIMAT)+SML)
      END DO
!
      DO 7800 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(.not.mat_cal(IIMAT)) goto 7800
      DO 470 ICVL=ICVS,ICVE    !i=1,NCVIN
      p0(ICVL)=r1(ICVL)+beta(IIMAT)*(p0(ICVL)-omega(IIMAT)*q1(ICVL))
  470 CONTINUE
 7800 CONTINUE
      s(:)=t(:)
!
 1000 CONTINUE
!
      ier=1
!
  999 CONTINUE
!
      aeps=max(aeps,rmaxs)
      reps=max(reps,rmaxs/(rmaxs0+SML))
!
      DO 650 IIMAT=1,NMAT
c      IMAT=MAT_NO(IIMAT)
c      ICVS=MAT_INDEX(IIMAT-1)+1
c      ICVE=MAT_INDEX(IIMAT)
c      if(mat_cal(IIMAT)) then
        DO 500 ICVL=ICVS,ICVE   !i=1,NCV
        bb(ICVL)=x(ICVL)*z(ICVL)
  500   CONTINUE
c      endif
 650  CONTINUE

!
c      IF(NPE.GT.1) THEN
c        CALL SOLVER_SEND_RECV (1,MXCV,NCV,BB)
c      ENDIF
!
      DEALLOCATE (d,p0,p1,q0,q1,r0,r1,r2,x,z,diag_a,a0)
!
      RETURN

!-------> end of Subroutine Bicgstab_band()

      END

!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine ICCGLTRI(IUT0,NRD,MRD,NVRD,RevIndex,
     &			RDIndex,RDValue,BB,TT,aeps,reps,itrmax,itr,ier)
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  1. CG solver with diagonal scaling
!     & incomplete Cholesky decomposition
!
!     NCV   : no. of elements
!     mm   : maximul no. of off-diagonal elements
!     au   : upper triangular part of coefficient matrix
!     bb   : right hand side vector (in)
!     TT   : solution vector (out)
!     aeps : absolute error of tolerance for convergence (in)
!          : absolute error of residual at convergence (out)
!     reps : relative error of tolerance for convergence (in)
!          : relative error of residual at convergence (out)
!     itrmax  : maximul iteration counts (in)
!     itr     : iteration counts at convergence (out)
!     ier  : return code
!
!
! --- [module arguments]
!
      use module_radslv,only : d=> W1K1
      use module_radslv,only : p=> W1K2
      use module_radslv,only : q=> W1K3
      use module_radslv,only : r=> W1K4
      use module_radslv,only : x=> W1K5
      use module_radslv,only : z=> W1K6
	use module_radslv,only : AA=> W1K7
	use module_radslv,only : B=> W1K8
	use module_radslv,only : AAStaID=> WKN1
	use module_radslv,only : WK=> WKN2

!
      IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (MXMAT=1,NMAT=1)
!
c      implicit none
!
! --- [dummy arguments]
!	
      real*8 ,intent(inout) :: bb(NRD),TT(NRD)
      INTEGER,intent(inout) :: RDIndex(NRD+1),RevIndex(NVRD)
	REAL*8,intent(inout)  :: RDValue(MRD)
	real*8 ,intent(inout) :: aeps,reps
      integer,intent(inout) :: itr
      integer,intent(out)   :: ier
	
	
	MVRD = NVRD*(NVRD+1)/2 + 1

	ALLOCATE (d(NRD),p(NRD),q(NRD),r(NRD),x(NRD),z(NRD),
     &		  AA(MVRD),AAStaID(NVRD+1),B(NVRD),WK(NRD),
     &		  stat=ierr)
	if(ierr.ne.0) stop 'stop at allocation -1- in ICCGLTRI()'
	AA(:) = RDValue(:)

!
!-< 0. transfer arrays >-
!
!		RDValue	--> AA
!		BB		--> B
!		RDIndex	--> AAStaID
	
	
	WK(:) = 0
	DO I=1,NVRD
		IRD=RevIndex(I)
		WK(IRD) = I
	END DO

	!B -->B
	DO I=1,NVRD
		IRD=RevIndex(I)
		B(I) = BB(IRD)
		DO JRD=1,IRD
			IF(WK(JRD).EQ.0) THEN
				ID = RDIndex(IRD)+JRD
				B(I)=B(I)-RDValue(ID)*TT(JRD)
			END IF
		END DO
		DO JRD=IRD+1,NRD
			IF(WK(JRD).EQ.0) THEN
				ID = RDIndex(JRD)+IRD
				B(I)=B(I)-RDValue(ID)*TT(JRD)
			END IF
		END DO
	END DO

	!RDValue,RDInde --> AA,AAStaID
	IVRD=0
	DO I=1,NVRD
		IRD=RevIndex(I)
		AAStaID(I)= IVRD
		DO JRD=1,IRD
			IF(WK(JRD).GT.0) THEN
				IVRD = IVRD + 1
				ID = RDIndex(IRD)+JRD
				AA(IVRD)=RDValue(ID)
			END IF
		END DO
	END DO
	AAStaID(NVRD+1)= IVRD



!-< 1. preliminary set >-
	
!
!
      ier=0
      jer=0
      it1=1
      it2=itrmax
      itr=0
      maxitr=it2
      
!
      
	d(:)=0.d0
      p(:)=0.d0
      q(:)=0.d0
      r(:)=0.d0
      x(:)=0.d0
      z(:)=0.d0


!
!-< 2. Preconditioning >-
!
!--< 2.1 diagonal scaling >--
!
      
      DO IRD=1,NVRD

		ID = AAStaID(IRD)+IRD
		d(IRD)=1.d0/(dsqrt(AA(ID)))
		z(IRD)=d(IRD)
      END DO
!
      DO IRD=1,NVRD
		DO JRD=1,IRD-1
			ID = AAStaID(IRD)+JRD
			AA(ID)=AA(ID)*d(IRD)*d(JRD)
		END DO
	END DO
 
!
      rmax=0.d0
      DO 220 IRD=1,NVRD    !i=1,NCV
		B(IRD)=B(IRD)*d(IRD)
		x(IRD)=0.d0
		rmax=max(rmax,abs(B(IRD)))
		d(IRD)=1.d0			!!!
  220 CONTINUE
!
      IF(rmax.le.aeps) THEN
        goto 2000
      END IF
      rmax0=rmax


!
!--< 2.2 incomplete Cholesky decomposition >--
!
      DO 230 IRD=1,NVRD
c		IF(d(IRD).le.0.d0) then			!no need, if pre-scaling is done
c			WRITE(IUT0,*) 'Error abort in ICCGLTRI()'
c			WRITE(IUT0,*) 'Aii = 0'
c			STOP
c		END IF
		d(IRD)=1.d0/(d(IRD))
		DO 231 JRD=1,IRD-1
			ID = AAStaID(IRD)+JRD
			d(JRD)=d(JRD)-d(IRD)*AA(ID)*AA(ID)
  231		CONTINUE
  230 CONTINUE
!
 2000 CONTINUE


!
!-< 3. Initial set for iteration >-
!
!
  888 CONTINUE
!
      rmax=0.d0
      DO IRD=1,NVRD
		r(IRD)=x(IRD)
      END DO
!
      DO 301 IRD=1,NVRD  !i=1,NCV
        DO 301 JRD=1,IRD-1
		ID = AAStaID(IRD)+JRD
		r(IRD)=r(IRD)+AA(ID)*x(JRD)
		r(JRD)=r(JRD)+AA(ID)*x(IRD)
  301 CONTINUE
!
      DO 302 IRD=1,NVRD    !i=1,NCV
		r(IRD)=B(IRD)-r(IRD)
		rmax=max(rmax,abs(r(IRD)))
		ID = AAStaID(IRD)+IRD
		AA(ID)=0.d0
  302 CONTINUE
!
      if(rmax.le.max(aeps,reps*rmax0)) goto 600
!
      DO 310 IRD=1,NVRD    !i=1,NCV
		ID = AAStaID(IRD)+IRD
		p(IRD)=d(IRD)*(r(IRD)-AA(ID))
		DO 311 JRD=1,IRD-1
			ID1 = AAStaID(IRD)+JRD
			AA(ID)=AA(ID)+AA(ID1)*p(IRD)
  311		CONTINUE
  310 CONTINUE
!
      c1=0.d0
      DO 320 IRD=NVRD,1,-1    !i=NCV,1,-1
		sm=0.d0
		DO 321 JRD=1,IRD-1
			ID1 = AAStaID(IRD)+JRD
			sm=sm+AA(ID1)*p(JRD)
  321		CONTINUE
		p(IRD)=p(IRD)-d(IRD)*sm
		c1=c1+r(IRD)*p(IRD)
  320 CONTINUE
!
      if(c1.eq.0.d0) jer=1
!
      if(jer.gt.0) then
        DO 330 IRD=1,NVRD   !i=1,NCV
		x(IRD)=x(IRD)+p(IRD)
  330   CONTINUE
        jer=0
        it1=it1+1
        itr=it1
        if(it1.le.it2) goto 888
      endif



!
!-< 4. iteration >-
!
      DO 1000 itr=it1,it2
      if(itr>itrmx) itrmx=itr
!
      DO 400 IRD=1,NVRD
		q(IRD)=p(IRD)
  400 CONTINUE

      DO 401 IRD=1,NVRD
	  DO 401 JRD=1,IRD-1
		ID = AAStaID(IRD)+JRD
		q(IRD)=q(IRD)+AA(ID)*p(JRD)
		q(IRD)=q(IRD)+AA(ID)*p(IRD)
  401 CONTINUE
      c2=0.d0
!      
      DO 402 IRD=1,NVRD   !i=1,NCV
		c2=c2+p(IRD)*q(IRD)
  402 CONTINUE
!
      if(c2==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      alpha=c1/c2
      rmax=0.d0
      DO 410 IRD=1,NVRD
        x(IRD)=x(IRD)+alpha*p(IRD)
        r(IRD)=r(IRD)-alpha*q(IRD)
        rmax=max(rmax,abs(r(IRD)))
        ID = AAStaID(IRD)+IRD
	  AA(ID)=0.d0
  410 CONTINUE
!
      if(rmax.le.min(aeps,reps*rmax0) ) goto 600
!
      DO 420 IRD=1,NVRD   !i=1,NCV
	  ID=AAStaID(IRD)+IRD
        q(IRD)=d(IRD)*(r(IRD)-AA(ID))
        DO 421 JRD=1,IRD-1
          ID1=AAStaID(IRD)+JRD
          AA(ID)=AA(ID)+AA(ID1)*q(IRD)
  421   CONTINUE
  420 CONTINUE

      c3=0.d0
      DO 430 IRD=NVRD,1,-1   !i=NCV,1,-1
        sm=0.d0
        DO 431 JRD=1,IRD-1
		ID=AAStaID(IRD)+JRD
          sm=sm+AA(ID)*q(JRD)
  431   CONTINUE
        q(IRD)=q(IRD)-d(IRD)*sm
        c3=c3+r(IRD)*q(IRD)
  430 CONTINUE
!
      if(c3==0.d0) then
        jer=1
        it1=itr
        goto 888
      endif
!
      beta=c3/(c1)
      c1=c3
      DO 440 IRD=1,NVRD    !i=1,NCV
        p(IRD)=q(IRD)+beta*p(IRD)
  440 CONTINUE
!
 1000 CONTINUE
!
 600  CONTINUE
!
!
!
      
      DO 500 IRD=1,NVRD
		B(IRD)=x(IRD)*z(IRD)
  500 CONTINUE
      
	DO 550 IRD=1,NVRD
		I = RevIndex(IRD)
		TT(I)=B(IRD)
  550 CONTINUE

 !
      itr=itrmx
!
      
	DEALLOCATE (d,p,q,r,x,z,AA,AAStaID,B)      

	
	RETURN
      
!-------> end of Subroutine ICCGLTRI()
	
	END 

C======================================================================C
      SUBROUTINE GauseLTRI(NVRD,NRD,MRD,RDValue,BB,TT,RDIndex,RevIndex)
C======================================================================C
! --- direct solver
!
	use module_radslv,only : X=> W1K1
	use module_radslv,only : AA=> W1K7
	use module_radslv,only : B=> W1K8
	use module_radslv,only : AAStaID=> WKN1
	use module_radslv,only : WK=> WKN2

      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(in)  :: RDIndex(NRD+1),RevIndex(NVRD)
      real*8 ,intent(in)  :: RDValue(MRD),BB(NRD)
      real*8 ,intent(out) :: TT(NRD)
!
! --- [local entities]
!
      
	MVRD = NVRD*(NVRD+1)/2 + 1

	ALLOCATE( AA(MVRD),AAStaID(NVRD+1),B(NVRD),WK(NRD),X(NRD),
     &		  stat=ierr)

	WK(:) = 0
	DO I=1,NVRD
		IRD=RevIndex(I)
		WK(IRD) = I
	END DO

	!B -->B
	DO I=1,NVRD
		IRD=RevIndex(I)
		B(I) = BB(IRD)
		DO JRD=1,IRD
			IF(WK(JRD).EQ.0) THEN
				ID = RDIndex(IRD)+JRD
				B(I)=B(I)-RDValue(ID)*TT(JRD)
			END IF
		END DO
		DO JRD=IRD+1,NRD
			IF(WK(JRD).EQ.0) THEN
				ID = RDIndex(JRD)+IRD
				B(I)=B(I)-RDValue(ID)*TT(JRD)
			END IF
		END DO
	END DO

	!RDValue,RDInde --> AA,AAStaID
	IVRD=0
	DO I=1,NVRD
		IRD=RevIndex(I)
		AAStaID(I)= IVRD
		DO JRD=1,IRD
			IF(WK(JRD).GT.0) THEN
				IVRD = IVRD + 1
				ID = RDIndex(IRD)+JRD
				AA(IVRD)=RDValue(ID)
			END IF
		END DO
	END DO
	AAStaID(NVRD+1)= IVRD


!
!
!/forward elimination/
      DO 100 I = 1, NVRD
		ID = AAStaID(I)+I			!aii
		diagV=1.0/AA(ID)				
		DO 120 J = I+1, NVRD
			
			ID1= AAStaID(J)+I			!aji/aii
			coef = AA(ID1)*diagV	

			DO 121 K = I+1,J		
										!k=j,NVRD is not need to cal. because of array symmetry
				ID2 = AAStaID(K)+I		!Aik, symmetric, Aki

				ID3 = AAStaID(J)+K		!Ajk

				AA(ID3) = AA(ID3) - coef*AA(ID2)

  121			CONTINUE
			B(J) = B(J) - coef*B(I)
	
  120		CONTINUE
  100 CONTINUE
!
!/backward substitution/
      DO 200 I = NVRD, 1, -1
		v1 = B(I)
		DO 220 J = I+1, NVRD		
			ID = AAStaID(J)+I		!v1 = b - Aij*Xj, J=i+1,N
			v1 = v1 - AA(ID)*X(J)
  220		CONTINUE
		ID = AAStaID(I)+I
		X(I) = v1/AA(ID)			!X = V1/Aii
	
  200 CONTINUE


!
!
      DO 300 I = 1, NVRD
		IRD = RevIndex(I)
		TT(IRD)= X(I)
  
  300 CONTINUE
!
      
	
	
	DEALLOCATE( AA,AAStaID,B,WK,X)

	RETURN

!--------> end of GauseLTRI

      END
C======================================================================C
      SUBROUTINE GauseBAND(NVRD,NRD,MRD,RDValue,RDId,BB,TT,
     &					 RDIndex,RevIndex)
C======================================================================C
! --- direct solver
!
	use module_radslv,only : X=> W1K1
	use module_radslv,only : AA=> W1K7
	use module_radslv,only : B=> W1K8
	use module_radslv,only : AAStaID=> WKN1
	use module_radslv,only : WK=> WKN2

      IMPLICIT REAL*8(A-H,O-Z)
      integer RDIndex(NRD+1),RevIndex(NVRD),RDId(MRD)
      real*8  RDValue(MRD)
      real*8  TT(NRD),BB(NRD)
!
! --- [local entities]
!
      
	MVRD = NVRD*(NVRD+1)/2 + 1

	ALLOCATE( AA(MVRD),AAStaID(NVRD+1),B(NVRD),WK(NRD),X(NRD),
     &		  stat=ierr)

	WK(:) = 0
	DO I=1,NVRD
		IRD=RevIndex(I)
		WK(IRD) = I
	END DO

	!B -->B
	DO I=1,NVRD
		IRD=RevIndex(I)
		B(I) = BB(IRD)
		DO JRD=1,IRD-1
			IF(WK(IRD).EQ.0) THEN
				ID = RDIndex(IRD)+JRD
				B(I)=B(I)-RDValue(ID)*TT(JRD)
			END IF
		END DO
		DO JRD=IRD+1,NRD
			IF(WK(IRD).EQ.0) THEN
				ID = RDIndex(JRD)+IRD
				B(I)=B(I)-RDValue(ID)*TT(JRD)
			END IF
		END DO
	END DO

	!RDValue,RDInde --> AA,AAStaID
	IVRD=0
	DO I=1,NVRD
		IRD=RevIndex(I)
		AAStaID(I)= IVRD
		DO JRD=1,IRD
			IF(WK(IRD).EQ.1) THEN
				IVRD = IVRD + 1
				ID = RDIndex(IRD)+JRD
				AA(IVRD)=RDValue(ID)
			END IF
		END DO
	END DO
	AAStaID(NVRD+1)= IVRD


!
!
!/forward elimination/
      DO 100 I = 1, NVRD
		ID = AAStaID(I)+I			!aii
		diagV=1.0/AA(ID)				
		DO 120 J = I+1, NVRD
			
			ID1= AAStaID(J)+I			!aji/aii
			coef = AA(ID1)*diagV	

			DO 121 K = I+1,J		
										!k=j,NVRD is not need to cal. because of array symmetry
				ID2 = AAStaID(K)+I		!Aik, symmetric, Aki

				ID3 = AAStaID(J)+K		!Ajk

				AA(ID3) = AA(ID3) - coef*AA(ID2)		!Ajk = Ajk-Aji/Aii*Aik
														!j,k = i+1,N
  121			CONTINUE
			B(J) = B(J) - coef*B(I)
  120		CONTINUE
  100 CONTINUE
!
!/backward substitution/
      DO 200 I = NVRD, 1, -1
		v1 = B(I)
		DO 220 J = i+1, NVRD		
			ID = AAStaID(J)+I		!v1 = b - Aij*Xj, J=i+1,N
			v1 = v1 - AA(ID)*X(J)
  220		CONTINUE
		ID = AAStaID(I)+I
		X(I) = v1/AA(ID)			!X = V1/Aii
  200 CONTINUE


!
!
      DO 300 IRD = 1, NRD
		TT(IRD)= X(WK(IRD))
  
  300 CONTINUE
!
      	
	
	DEALLOCATE( AA,AAStaID,B,WK,X)

	RETURN

!--------> end of GauseBAND()

      END

C======================================================================C
	SUBROUTINE SORBAND(IUT0,RowStaPoint,Aij,XX,BB,ElemId,RevIndex,
     &			aeps,reps,RLX,NN,MM,MID,NVRD,itr)
C======================================================================C
	
	IMPLICIT REAL*8(A-H,O-Z)
	INTEGER ElemID(MID),RowStaPoint(NN+1),RevIndex(NVRD)
	REAL*8  Aij(MM)
	REAL*8  BB(NN),XX(NN)
	
	!rlx = 0.75

C---------------------------------------------------------------
C
C	This is a solver by Successive Overrelaxation Method 
C		A * X = B
C
C---------------------------------------------------------------

!-------
	err1 = aeps + 1
	err2 = reps + 1

	itr = 0
	DO 1000 WHILE(err1.GT.aeps.OR.err2.GT.reps)
		
		err1 = 0
		err2 = 0
		itr = itr + 1

		DO 200 IVRD = 1,NVRD

		  IN = RevIndex(IVRD)
		  
		  aa = 0
		  DO 100 IM = RowStaPoint(IN)+2,RowStaPoint(IN+1)
			JN = ElemID(IM)
			aa = aa + Aij(IM)*XX(JN)
100		  CONTINUE
		
	      !vsum = sum(Aij(RowStaPoint(IN)+2:RowStaPoint(IN+1)))
	
		  vnew = (BB(IN)-aa)/Aij(RowStaPoint(IN)+1)
		  vold = XX(IN)
		  XX(IN) = (1.0-RLX)*vold + RLX*vnew
		  
		  vv =  ABS(XX(IN)-vold)
		  err1 = err1 + vv*vv
		  IF(XX(IN).NE.0) err2 = err2 + (vv/XX(IN))*(vv/XX(IN))

200		CONTINUE
		err1 = sqrt(err1)/NVRD
		err2 = sqrt(err2)/NVRD

		WRITE(IUT0,*) 'Iteration = ',itr
		WRITE(IUT0,*) '   aeps = ',err1, ' / ', aeps
		WRITE(IUT0,*) '   reps = ',err2, ' / ', reps
		

1000	CONTINUE

		
	RETURN

!------> end of SorBAND()
	END 

C======================================================================C
	SUBROUTINE SORLTRI(IUT0,RowStaPoint,Aij,XX,BB,RevIndex,
     &				   aeps,reps,RLX,NN,MM,NVRD,itr)
C======================================================================C

	IMPLICIT REAL*8(A-H,O-Z)
	INTEGER RowStaPoint(NN+1),RevIndex(NVRD)
	REAL*8  Aij(MM)
	REAL*8  BB(NN),XX(NN)
	

C---------------------------------------------------------------
C
C	This is a solver by Successive Overrelaxation Method
C		A * X = B
C
C---------------------------------------------------------------


!-------
	IF(RLX.LE.0) RLX=0.75
	
	err1 = aeps + 1
	err2 = reps + 1

	itr = 0
	DO 1000 WHILE(err1.GT.aeps.OR.err2.GT.reps)
		
		err1 = 0
		err2 = 0
		itr = itr + 1

		DO 200 IVRD = 1,NVRD
		  
		  I = RevIndex(IVRD)

		  aa = 0.0
		  ID = RowStaPoint(I)
		  DO 100 J = 1,I-1
			aa = aa + Aij(ID+J)*XX(J)
100		  CONTINUE

		  
		  DO 150 J = I+1,NN
			ID1 = RowStaPoint(J)
			aa = aa + Aij(ID1+I)*XX(J)
150		  CONTINUE
		
	      !vsum = sum(Aij(RowStaPoint(IN)+2:RowStaPoint(IN+1)))
	
		  vnew = (BB(I)-aa)/Aij(ID+I)
		  vold = XX(I)
		  XX(I) = (1.0-RLX)*vold + RLX*vnew
		  
		  vv = ABS(XX(I)-vold)
		  err1 = err1 + vv*vv
		  
		  IF(XX(I).NE.0) err2 = err2 + (vv/XX(I))*(vv/XX(I))
		
200		CONTINUE
		err1 = sqrt(err1)/NVRD
		err2 = sqrt(err2)/NVRD

		!WRITE(IUT0,*) 'Iteration = ',itr
		!WRITE(IUT0,*) '   aeps = ',err1, ' / ', aeps
		!WRITE(IUT0,*) '   reps = ',err2, ' / ', reps

1000	CONTINUE

	WRITE(IUT0,*) 'Iteration = ',itr
	WRITE(IUT0,*) '   aeps = ',err1, ' / ', aeps
	WRITE(IUT0,*) '   reps = ',err2, ' / ', reps


	RETURN
!---------> end of SorLTRI()

	END

C======================================================================C
	SUBROUTINE GetHeatFlux(NW,NE,NRD0,NRD,MRD,MID,N2D,TT,DIVQ,QW,WNEB,
     &			   VIndex,RDValue,RDId,RDIndex,BB,WArea,Volume,
     &				   WElem,StackOpt,QSUM)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (PAI = 3.141592653589793,SHIGMA = 5.67e-8)
	!
	INTEGER RDId(MID),RDIndex(NRD+1),VIndex(NRD0),WNEB(NW)
	INTEGER WElem(N2D+2,NW),numw(50)
	REAL*8  RDValue(MRD),WArea(NW),Volume(NE)
	REAL*8  TT(NRD),DIVQ(NRD0),QW(NW),BB(NRD)
	REAL*8	QSUM(2,50)	!1:		total heat flux
				!2:		average heat flux
	CHARACTER*20 StackOpt

C-----------------------------------------------------------------------------------------C
C	Energy conservation equation for unknown element:
C			PiiXi = sum(Dij*Xj) + Bi	-----
C											 |
C											 |
C			Aii*Xi + sum(Aij*Xj) = Bi	<----
C			Aij = -Dij		(for i!=j)
C			Aii = Pii-Dii
C			Bi = HSource			
C
C	DivQ for known element:
C			DivQ = Aii*shigma*Ti^4 + sum(Aij*shigma*Tj^4) - Bi
C			
C			Positive DivQ means Qgo > Qcome


	
	DIVQ = 0.0

	st = 0
	DO i=1,ne
		st = st + TT(i)**2
	end DO
	st = sqrt(st/ne)
	write(*,*) '	----------------------------'
	write(*,*) '		ave_t=',st
	write(*,*) '	----------------------------'
	
	IF(StackOpt(1:4).EQ.'LTRI') THEN
		DO 1000 II = 1,NRD0
			I = VIndex(II)
			IF(I.GT.0) THEN
				ID = RDIndex(I)
			DivQ(II) = SHIGMA*RDValue(ID+I)*TT(I)**4
	
				DO J=1,I-1
		 DivQ(II) = DivQ(II)+RDValue(ID+J)*SHIGMA*TT(J)**4
    	
									!Aij < 0
				END DO
				
				DO J=I+1,NRD
					ID1 = RDIndex(J)
		 DivQ(II) = DivQ(II)+RDValue(ID1+I)*SHIGMA*TT(J)**4

									!Aji = Aij < 0
				END DO

								
			END IF

		
1000		CONTINUE

	ELSE

		DO 2000 II = 1,NRD0
		I = VIndex(II)
		IF(I.GT.0) THEN
		ID = RDIndex(I)
		DivQ(II) = SHIGMA*RDValue(ID+1)*TT(I)**4 - BB(I)
		DO J=ID+2,RDIndex(I+1)
		J1 = RDId(J)
		DivQ(II) = DivQ(II)+RDValue(J)*SHIGMA*TT(J1)**4
		END DO
				

	END IF
2000		CONTINUE

	END IF

	DO I=1,NE
		DivQ(I)=DivQ(I)/Volume(I)
	END DO
	DO I=1,NW
		DivQ(NE+I)=DivQ(NE+I)/WArea(I)
	END DO




	!get heat flux summary of each wall group
	QSUM = 0.0
	NUMW = 0
	DO IW = 1,NW
		IN = WElem(1,IW)
		ID = WElem(IN+2,IW)
		IIT = IW+NE
		IF(ID.GT.0) THEN
			QSUM(1,ID) = QSUM(1,ID) + DivQ(IIT)*WArea(IW)
			NUMW(ID) = NUMW(ID) + 1
		END IF
	END DO

	OPEN(1,FILE='Summary.dat',POSITION='APPEND')
	WRITE(1,*) 'HEAT FLUX ...'
	DO ID = 1,50
		IF(NUMW(ID).GT.0) THEN
			QSUM(2,ID) = QSUM(1,ID)/NUMW(ID)
		WRITE(1,*) ID,':   ', NUMW(ID),QSUM(1,ID),QSUM(2,ID)
		END IF
	END DO
	CLOSE(1)
	
	

	RETURN
!---------> end of GetHeatFlux()

	END




C======================================================================C
	
	SUBROUTINE RDVElemToNode_old
     &       (IUT0,MM,ME,MP,NP,NE,NW,NWT,NRD0,NRD,
     &		MRD,MRD0,MID,N2D,N3D,NKind,StackOpt,NPIN,
     &		ELEM)
C======================================================================C
	
	USE MODULE_RADWORK, only : NRDValue	=> VWK1
	USE MODULE_RADWORK, only : NPROP	=> VWK2
	USE MODULE_RADWORK, only : RevIndex => WK1
	USE MODULE_RADWORK, only : NRDIndex	=> WK5
	USE MODULE_RADWORK, only : NVIndex	=> WK6

	USE MODULE_RADSXF, only : PRATIO,RaNeb,AreaVol
	USE MODULE_RADSXF, only : RDValue,RDId,RDIndex,RDN1,RDN2
	USE MODULE_RADSXF, only : WElem,ElemKind,NNPE,Volume,WArea
	USE MODULE_RADSXF, only : PROP,VIndex,NodeBounID

	IMPLICIT REAL*8(A-H,O-Z)
	
	INTEGER ELEM(N3D,ME)
	CHARACTER*20 StackOpt,StackOpt2
	

C
C	1.	Index for participating nodes, allocate working arrays
C

	!mark participating nodes to allocate arrays
	ALLOCATE(RevIndex(NRD),NVIndex(NP),NPROP(NP),NRDIndex(NP+1))
	ALLOCATE(PRATIO(NP),AreaVol(NP),RaNeb(NP),stat=ierr)


	NVIndex(1:NP)=0
	CALL GetNodeIndex(ME,NE,NW,NWT,NP,N2D,N3D,NKind,WElem,
     &		ELEM,ElemKind,NNPE,VIndex,NVIndex,NRD2,NPIN,NPBI)
	
	RevIndex=0
	DO I=1,NRD0
		IRD=VIndex(I)
		IF(IRD.GT.0) RevIndex(IRD)=I
	END DO
		
	StackOpt2='LTRI'
	NRDIndex = 0
	DO I = 2, NRD2+1
		NRDIndex(I) = NRDIndex(I-1) + I-1
	END DO
	MRD2 = NRDIndex(NRD2+1)
		

	ALLOCATE(NRDValue(MRD2),stat=ierr)
	if(ierr.ne.0) stop 'stop at allocation -1- in RDVElemToNode()'


C
C   2	. Node Equivalent Property
C

	!NPROP is used for radiation heat transfer calculation
	!	and hence the some nodes has no value

	AreaVol(:)=0.d0
	NPROP(:)=0.d0
	DO IW = 1, NW
		  NN=0
		  DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			IF(NVIndex(IP).LT.0) THEN
				NN=NN+1
			END IF
		  END DO
				   
		  DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			IF(NVIndex(IP).LT.0) THEN
				NPROP(IP)=NPROP(IP)+PROP(NE+IW)/NN
!AreaVol(IP)=AreaVol(IP)+WArea(IW)/NN	
			END IF
		  END DO
	END DO
	

	DO IE=1,NE
		  NN=0
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
		IF(NVIndex(IP).LE.NPIN.AND.NVIndex(IP).GT.0) THEN
	          NN=NN+1
		END IF
		  END DO
		  		  
		  DO J=1,NNPE(ElemKind(IE))
		IP=ELEM(J,IE)
		IF(NVIndex(IP).LE.NPIN.AND.NVIndex(IP).GT.0) THEN
		  NPROP(IP)=NPROP(IP)+PROP(IE)/NN
		  !AreaVol(IP)=AreaVol(IP)+Volume(IE)/NN
		END IF
	  END DO

	END DO



!AreaVol, PRATIO is used to give real heat flux divergence and heat flux
1234	PRATIO(:)=0.d0
	
	DO IW = 1, NW
		  DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			PRATIO(IP)=PRATIO(IP)+PROP(NE+IW)/WElem(1,IW)
			AreaVol(IP)=AreaVol(IP)-WArea(IW)/WElem(1,IW)
		  END DO
	END DO
	DO IE=1,NE
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
			IF(AreaVol(IP).GE.0.0d0) THEN
		  PRATIO(IP)=PRATIO(IP)+PROP(IE)/NNPE(ElemKind(IE))
		  AreaVol(IP)=AreaVol(IP)+Volume(IE)/NNPE(ElemKind(IE))
		END IF
		  END DO
	END DO
	
	

!RaNeb is to transfer heat flux divergence from internal nodes to the nodes in surface 
!	(transfering imaginary temperatures from interface(import)_nodes to export_nodes
!	for data communication needs no RaNeb)
	RaNeb(:)=0.d0
	DO IE=1,NE
	DO J=1,NNPE(ElemKind(IE))
	IP=ELEM(J,IE)
	RaNeb(IP)=RaNeb(IP)+PROP(IE)/NNPE(ElemKind(IE))
	END DO
	END DO

	

C
C   3	. Transfer values from cell to vertices
C

!---------------------------------------------------------------------------------
!	
!	Transfer values from cell to vertices
!
!---------------------------------------------------------------------------------
	NRDValue(:)=0.0
	IEEND=VIndex(NE)
	!NPROP(:) = 0
	
	IF(StackOpt(1:4).EQ.'LTRI') THEN
		!Transfor for LRTI mode,	 LTRI	=>	LTRI

		!emitting cell is internal element 
		DO 200 IE=1,NE
		  IRD = VIndex(IE)
		  IF(IRD.LE.0) CYCLE	!not a participating element
		  
		  JNNPE=0
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
			INRD=NVIndex(IP)
			IF(INRD.GT.0) JNNPE=JNNPE+1
		  END DO

		  IF(JNNPE.EQ.0) CYCLE
	  
		  DO 150 J=1,NNPE(ElemKind(IE))
			
			IP=ELEM(J,IE)
			INRD=NVIndex(IP)
			IF(INRD.LE.0) CYCLE 	
!transfer Value of INTERNAL element to INTERNAL node only
			
!prop => Nprop, Volume -> AreaVol
!NPROP(IP) = NPROP(IP)+PROP(IE)/JNNPE
			
!AreaVol(IP) = AreaVol(IP)+Volume(IE)/JNNPE
c			NBB(IP) = NBB(IP) + BB(IE)/JNNPE
c			NT(IP) = NT(IP) + T(IE)
c			NCX(IP)=NCX(IP)+1
					    
!transfer values to nodes, Dj->inrd,	(j=1,NRD2)
			DO 145 JRD=1,IEEND
			  ID = RDIndex(MAX(IRD,JRD))+MIN(JRD,IRD)
			  JE = RevIndex(JRD)

			  KNNPE=0
			  DO K=1,NNPE(ElemKind(JE))
				IP=ELEM(K,JE)
				JNRD=NVIndex(IP)
				IF(JNRD.GT.0) KNNPE=KNNPE+1	
!because JE is internal elem, its value transfers only to internal nodes
			  END DO
			  
		  IF(KNNPE.LE.0) CYCLE
			  
		  VA=RDValue(ID)/KNNPE/JNNPE		!!!!!!!!!
		  DO K=1,NNPE(ElemKind(JE))
			IP=ELEM(K,JE)
			JNRD=NVIndex(IP)
			IF(JNRD.GT.0) THEN
			IF(JNRD.LE.INRD) THEN
			ID3=NRDIndex(INRD)+JNRD
			ELSE
			ID3=NRDIndex(JNRD)+INRD
                        END IF
			NRDValue(ID3)=NRDValue(ID3)+VA
		    END IF
			  END DO

145			CONTINUE

			DO 146 JRD=IEEND+1,NRD
				ID = RDIndex(JRD)+IRD
				JW = RevIndex(JRD)-NE
				KNNPS=0
				DO K=1,WElem(1,JW)
					IP=WELEM(K+1,JW)
					JNRD=NVIndex(IP)
					IF(JNRD.LT.0) KNNPS=KNNPS+1
!because JW is boundary face,its value transfers only to boundary nodes
				END DO
				  
				IF(KNNPS.LE.0) CYCLE
				  
				VA=RDValue(ID)/KNNPS/JNNPE		!!!!!!!!!
				DO K=1,WElem(1,JW)
					IP=WELEM(K+1,JW)
					JNRD=NVIndex(IP)
					IF(JNRD.LT.0) THEN		!boundary nodes
					JNRD=-JNRD
					IF(JNRD.LE.INRD) THEN
						ID3=NRDIndex(INRD)+JNRD
					ELSE
						ID3=NRDIndex(JNRD)+INRD
					END IF
					NRDValue(ID3)=NRDValue(ID3)+VA
					END IF
				END DO
146			CONTINUE			

150		  CONTINUE
200		CONTINUE
	
		
		!emitting cell is boundary faces
		DO 300 IW=1,NWT
		  IRD = VIndex(NE+IW)
		  IF(IRD.LE.0) CYCLE	!not a participating face
		  
		  JNNPS=0
		  DO J=1,WElem(1,IW)
			IP=WELEM(J+1,IW)
			INRD=NVIndex(IP)
			IF(INRD.LT.0) JNNPS=JNNPS+1
!boundary face's valve transfers to boundary nodes noly
		  END DO
		  IF(JNNPS.EQ.0) CYCLE
		  
		  DO 250 J=1,WElem(1,IW)
			
			IP=WELEM(J+1,IW)
			INRD=NVIndex(IP)
			

			IF(INRD.GE.0) CYCLE
!transfer Value of boundary face to boundary nodes only
			
			INRD=-INRD
			!prop => Nprop
			!NPROP(IP) = NPROP(IP)+PROP(IW+NE)/JNNPS
		
			!AreaVol(IP) = AreaVol(IP)+WArea(IW)/JNNPS
c			NBB(IP) = NBB(IP) + BB(IRD)/JNNPS
c			NT(IP) = NT(IP) + T(NE+IW)
c			NCX(IP)=NCX(IP)+1

!transfer values to nodes, Dj->inrd,	(j=1,NRD2)
			DO 245 JRD=1,IEEND
			    ID = RDIndex(IRD)+JRD
			    JE = RevIndex(JRD)
			
			    !absorbing cell is internal element
			    KNNPE=0
				DO K=1,NNPE(ElemKind(JE))
					IP=ELEM(K,JE)
					JNRD=NVIndex(IP)
			IF(JNRD.GT.0) KNNPE=KNNPE+1	
!because JE is internal elem, its value transfers only to internal nodes
			END DO
				  
				IF(KNNPE.LE.0) CYCLE
				  
				VA=RDValue(ID)/KNNPE/JNNPS		!!!!!!!!!
				DO K=1,NNPE(ElemKind(JE))
					IP=ELEM(K,JE)
					JNRD=NVIndex(IP)
					IF(JNRD.GT.0) THEN
						
							!jnrd is internal node, JNRD always < INRD
				IF(JNRD.LE.INRD) THEN	
								
					ID3=NRDIndex(INRD)+JNRD
				ELSE
					ID3=NRDIndex(JNRD)+INRD
				END IF
				NRDValue(ID3)=NRDValue(ID3)+VA
			    END IF
			END DO
245			CONTINUE

			DO 246 JRD=IEEND+1,NRD
			    ID = RDIndex(MAX(IRD,JRD))+MIN(IRD,JRD)
			    JW = RevIndex(JRD)-NE
				KNNPS=0
				DO K=1,WElem(1,JW)
					IP=WELEM(K+1,JW)
					JNRD=NVIndex(IP)
					IF(JNRD.LT.0) KNNPS=KNNPS+1
						!because JW is boundary face,its value transfers only to boundary nodes
				END DO
				  
				IF(KNNPS.LE.0) CYCLE
				  
				VA=RDValue(ID)/KNNPS/JNNPS		!!!!!!!!!
				DO K=1,WElem(1,JW)
					IP=WELEM(K+1,JW)
					JNRD=NVIndex(IP)
					IF(JNRD.LT.0) THEN		!boundary nodes
						JNRD=-JNRD
					IF(JNRD.LE.INRD) THEN
						ID3=NRDIndex(INRD)+JNRD
					ELSE
						ID3=NRDIndex(JNRD)+INRD
					END IF
					NRDValue(ID3)=NRDValue(ID3)+VA
				    END IF
				END DO
246			CONTINUE			

250		  CONTINUE
300		CONTINUE


		DO 400 I=1,NRD2
			DO J=1,I-1
				ID=NRDIndex(I)+J
				NRDValue(ID)=NRDValue(ID)*0.5
			END DO
400		CONTINUE





	
	!transfer for StackOpt='BAND' mode,		BAND => LTRI
	ELSE
	  
!emitting cell is internal element 
		DO 1200 IE=1,NE
		  IRD = VIndex(IE)
		  IF(IRD.LE.0) CYCLE	!not a participating element
		  
		  JNNPE=0
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
			INRD=NVIndex(IP)
			IF(INRD.GT.0) JNNPE=JNNPE+1
		  END DO

		  IF(JNNPE.EQ.0) CYCLE
		
			  
		  DO 1150 J=1,NNPE(ElemKind(IE))
			
			IP=ELEM(J,IE)
			INRD=NVIndex(IP)
		
			IF(INRD.LE.0) CYCLE	
				!transfer Value of INTERNAL element to INTERNAL node only
			
			
			!prop => Nprop
			!NPROP(IP) = NPROP(IP)+PROP(IE)/JNNPE
			
			!AreaVol(IP) = AreaVol(IP)+Volume(IE)/JNNPE
c			NBB(IP) = NBB(IP) + BB(IRD)/JNNPE
c			NT(IP) = NT(IP) + T(IE)
c			NCX(IP)=NCX(IP)+1
		 
!transfer values to nodes, Dj->inrd,	(j=1,NRD2)
			DO 1149 JMRD=RDIndex(IRD)+1,RDIndex(IRD+1)
			  JRD = RDId(JMRD)
			  JE = RevIndex(JRD)
						  
			  KNNPE=0
			  DO K=1,NNPE(ElemKind(JE))
				IP=ELEM(K,JE)
				JNRD=NVIndex(IP)
				IF(JNRD.GT.0) KNNPE=KNNPE+1	
!because JE is internal elem, its value transfers only to internal nodes
			  END DO
			  
			  IF(KNNPE.LE.0) CYCLE
			  
			  VA=RDValue(JMRD)/KNNPE/JNNPE		!!!!!!!!!
			  DO K=1,NNPE(ElemKind(JE))
				IP=ELEM(K,JE)
				JNRD=NVIndex(IP)
				IF(JNRD.GT.0) THEN
					IF(JNRD.LE.INRD) THEN
						ID3=NRDIndex(INRD)+JNRD
					ELSE
						ID3=NRDIndex(JNRD)+INRD
					END IF
					NRDValue(ID3)=NRDValue(ID3)+VA
			   END IF
			  END DO

1149			CONTINUE			
1150		  CONTINUE
1200		CONTINUE
	
		
		!emitting cell is boundary faces
		DO 1300 IW=1,NWT
		  IRD = VIndex(NE+IW)
		  IF(IRD.LE.0) CYCLE	!not a participating face
		  
		  JNNPS=WElem(1,IW)
		  DO J=1,WElem(1,IW)
			IP=WELEM(J+1,IW)
			INRD=NVIndex(IP)
			IF(INRD.LT.0) JNNPS=JNNPS+1		
!boundary face's valve transfers to boundary nodes noly
		  END DO
		  IF(JNNPS.EQ.0) CYCLE
		  
		  DO 1250 J=1,WElem(1,IW)
			
			IP=WELEM(J+1,IW)
			INRD=NVIndex(IP)
			IF(INRD.GE.0) CYCLE
			INRD=-INRD
!transfer Value of boundary face to boundary nodes only
			
			!prop => Nprop
			!NPROP(IP) = NPROP(IP)+PROP(IW+NE)/JNNPS
						
			!AreaVol(IP) = AreaVol(IP)+WArea(IW)/JNNPS
c			NBB(IP) = NBB(IP) + BB(IRD)/JNNPS
c			NT(IP) = NT(IP) + T(NE+IW)
c			NCX(IP)=NCX(IP)+1

!transfer values to nodes, Dj->inrd,	(j=1,NRD2)
			DO 1249 JMRD=RDIndex(IRD)+1,RDIndex(IRD+1)
			  JRD = RDId(JMRD)
			  JEW = RevIndex(JRD)
			
			  IF(JEW.LE.NE)  THEN		!absorbing cell is internal element
				  KNNPE=0
				  DO K=1,NNPE(ElemKind(JEW))
					IP=ELEM(K,JEW)
					JNRD=NVIndex(IP)
					IF(JNRD.GT.0) KNNPE=KNNPE+1
						!because JE is internal elem, its value transfers only to internal nodes
				  END DO
				  
				  IF(KNNPE.LE.0) CYCLE
				  
				  VA=RDValue(ID)/KNNPE/JNNPS		!!!!!!!!!
				  DO K=1,NNPE(ElemKind(JEW))
					IP=ELEM(K,JEW)
					JNRD=NVIndex(IP)
					IF(JNRD.GT.0) THEN
							!jnrd is internal node, JNRD always < INRD
					IF(JNRD.LE.INRD) THEN	
					ID3=NRDIndex(INRD)+JNRD
					ELSE
					ID3=NRDIndex(JNRD)+INRD
				END IF
					NRDValue(ID3)=NRDValue(ID3)+VA
				   END IF
				  END DO
			  ELSE
				  JW = JEW-NE
				  KNNPS=0
				  DO K=1,WElem(1,JW)
					IP=WELEM(K+1,JW)
					JNRD=NVIndex(IP)
					IF(JNRD.LT.0) KNNPS=KNNPS+1
!because JW is boundary face,its value transfers only to boundary nodes
				  END DO
				  
				  IF(KNNPS.LE.0) CYCLE
				  
				  VA=RDValue(ID)/KNNPS/JNNPS
				  DO K=1,WElem(1,JW)
					IP=WELEM(K+1,JW)
					JNRD=NVIndex(IP)
					IF(JNRD.LT.0) THEN
						JNRD=-JNRD
						IF(JNRD.LE.INRD) THEN
						ID3=NRDIndex(INRD)+JNRD
						ELSE
						ID3=NRDIndex(JNRD)+INRD
						END IF
					NRDValue(ID3)=NRDValue(ID3)+VA
				   END IF
				  END DO

			  END IF

1249			CONTINUE			
1250		  CONTINUE
1300		CONTINUE

	END IF


C
C   4	. Conservation check
C
!---------------------------------------------------------------------
!	conservation check ...
!---------------------------------------------------------------------
	v1 = sum(volume(1:ne))
	w1 = sum(warea(1:nw))
	r1 = sum(RDValue(1:mrd))
	r2 = sum(NRDValue(1:mrd2))
	v2 = 0
	w2 = 0
	DO I=1,NP
		IF(AreaVol(I).GT.0) THEN
			v2 = v2+AreaVol(I)
		ELSE
			w2 = w2-AreaVol(I)
		END IF
	END DO

	exc1=0
	DO II=1,NRD0
		I=VIndex(II)
		IF(I.LE.0) CYCLE
		vv=0
		DO J=1,I
			ID=RDIndex(I)+J
			vv=vv+RDValue(ID)
		END DO
		DO J=I+1,NRD
			ID=RDIndex(J)+I
			vv=vv+RDValue(ID)
		END DO
		exc1=exc1+ abs(PROP(II)-vv)
	END DO
	
	exc2=0
	DO II=1,NP
		I=ABS(NVIndex(II))
		IF(I.LE.0) CYCLE
		vv=0
		DO J=1,I
			ID=NRDIndex(I)+J
			vv=vv+NRDValue(ID)
		END DO
		DO J=I+1,NRD2
			ID=NRDIndex(J)+I
			vv=vv+NRDValue(ID)
		END DO
		exc2=exc2+ (NPROP(II)-vv)
	
	END DO
	vv1=sum(prop(1:NRD0))
	vv2=sum(nprop(1:NP))

	write(IUT0,*) '  Conservation of Total-Exchange-Area ::'
	write(IUT0,*) '  Cell summary values'
cwrite(IUT0,'("  Volume=",F8.3," Area=",F8.3," RD=",F8.3)') v1,w1,r1 
	write(IUT0,'("   TEA= ",F10.6, ", Residence= ",F10.6)') r1,exc1
	write(IUT0,*) '  Vertex summary values'
cwrite(IUT0,'("  Volume=",F8.3," Area=",F8.3," RD=",F8.3)') v2,w2,r2
	write(IUT0,'("   TEA= ",F10.6, ", Residence= ",F10.6)') r2,exc2
	write(IUT0,'("   NRD= ",I6, ", NP= ",I6)') NRD2,NP

C
C   5	. Replace RDValue
C
!---------------------------------------------------------------------
!	replace RDValue ...
!---------------------------------------------------------------------
	MRD  = MRD2
	MRD0 = MRD
	NRD  = NRD2
	NRD0 = NP
	MID  = 1
	StackOpt='LTRI'

	DEALLOCATE(RDValue,RDIndex,RDId,PROP,VIndex)

	ALLOCATE(RDValue(MRD),RDIndex(NRD+1),PROP(NP),
     &		 RDId(MID),VIndex(NRD0),stat=ierr)
	if(ierr.ne.0) stop 'stop at allocation -2- in RDVElemToNode()'

	RDValue(:)	=	NRDValue(:)
	RDIndex(:)	=	NRDIndex(:)
	AreaVol(:)  =   ABS(AreaVol(:))

	DO I=1,NP
		VIndex(I) =	ABS(NVIndex(I))
		PROP(I)	 =	NPROP(I)
		IF(PROP(I).GT.0) THEN
			PRATIO(I)=PRATIO(I)/PROP(I)
		ELSE
			PRATIO(I)=0.d0		!no participation nodes
		END IF
		
	END DO




	
	DEALLOCATE(RevIndex,NRDValue,NRDIndex,NPROP,NVIndex)


	RETURN
!---------> end of RDVElemToNode_old()

	END

 






 
C======================================================================C
        SUBROUTINE GetNodeIndex(ME,NE,NW,NWT,NP,N2D,N3D,NKind,WElem,
     &		ELEM,ElemKind,NNPE,VIndex,NVIndex,NRD2,NPIN,NPBI)
C======================================================================C

	USE MODULE_RADSXF,ONLY : NodeBounID

	IMPLICIT REAL*8(A-H,O-Z)
	
	INTEGER WElem(N2D+2,NWT)
	INTEGER ELEM(N3D,ME),ElemKind(NE),NNPE(NKind)
	INTEGER VIndex(NE+NW),NVIndex(NP)


	NVIndex=0
	NodeBounID = 0



	DO 100 IE=1,NE
		IRD=VIndex(IE)
		IF(IRD.GT.0) THEN
			IK = ElemKind(IE)
			DO J=1,NNPE(IK)
				IP=Elem(J,IE)
				NVIndex(IP)=1				! -internal participating nodes
			END DO
		END IF
100	CONTINUE
	
	
	DO 110 IW=1,NW
		IIT=IW+NE
		IRD=VIndex(IIT)
		IF(IRD.GT.0) THEN
			NPS = WElem(1,IW)
			DO J=2,NPS+1
				IP=WElem(J,IW)
				IK2=WElem(WElem(1,IW)+2,IW)
		IF(NVIndex(IP).LT.0.AND.NVIndex(IP).NE.-IK2) THEN
					NVIndex(IP)=-100000
					NodeBounID(IP)=0
						!corner nodes, except
				ELSE
					NVIndex(IP)=-IK2
					NodeBounID(IP)=IK2
						!surface nodes
				END IF
								
			END DO
		END IF
110	CONTINUE
	
	!count internal and boundary node number
	NPIN=0
	NPBI=0
	DO IP=1,NP
		IF(NVIndex(IP).LT.0) THEN	!boundary nodes
			IF(NVIndex(IP).GT.-100000) THEN
				NPBI=NPBI+1
				NVIndex(IP)=-NPBI
			ELSE
				NVIndex(IP)=0
			END IF
		ELSE IF(NVIndex(IP).GT.0) THEN	!internal nodes
			NPIN=NPIN+1
			NVIndex(IP)=NPIN
		END IF
	END DO
	NRD2=NPBI+NPIN

	DO IP=1,NP
		IF(NVIndex(IP).LT.0) THEN	!boundary nodes
			NVIndex(IP)=NVIndex(IP)-NPIN
		END IF
	END DO





	RETURN
!---------> end of GetNodeIndex()

	END








C======================================================================C
	SUBROUTINE NodeToGPProp
     &  (IUT0,MM,ME,MP,NP,NE,NW,NWT,N2D,N3D,NKind,
     &			NRD,ELEM)
C======================================================================C
	

	USE MODULE_RADGROUP, only : Coef,NIndexToGP,GpProp,NumGPTran
	
	USE MODULE_RADSXF, only :   PRATIO,RaNeb,NodeBounID,AreaVol
	USE MODULE_RADSXF, only :   PROP,VIndex,WArea,Volume,WElem
	USE MODULE_RADSXF, only :   ElemKind,NNPE

	USE MODULE_RADWORK, only :  NVIndex	=> WK2
	USE MODULE_RADWORK, only :  VWK1
	USE MODULE_RADWORK, only :  NPROP	=> VWK2
	


	IMPLICIT REAL*8(A-H,O-Z)
	
	INTEGER IDUM(NKIND),ELEM(N3D,ME)
	CHARACTER*20 StackOpt,StackOpt2
	
	!mark participating nodes to allocate arrays
	ALLOCATE(PRATIO(NP),AreaVol(NP),RaNeb(NP))
	ALLOCATE(NVIndex(NP),VWK1(MAX(NP,NE+NW)),NPROP(NP))
	

c
c
	CALL GetNodeIndex(ME,NE,NW,NWT,NP,N2D,N3D,NKind,WElem,
     &		ELEM,ElemKind,NNPE,VIndex,NVIndex,NRD2,NPIN,NPBI)
	
	NVIndex(:) = ABS(NVIndex(:))

C
C	1. Node Equivalent Property
C
	

	PRATIO(:)=0.d0
	AreaVol(:)=0.d0
	RaNeb(:)=0.d0

	!NPROP is used for radiation heat transfer calculation
	!	and hence the some nodes has no value
	NumGPTran=0
	DO IW = 1, NW
		  NN=0
		  DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			IF(NVIndex(IP).GT.NPIN) THEN
				NN=NN+1
			END IF
		  END DO
		  NumGPTran=NumGPTran+NN
		   
		  DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			IF(NVIndex(IP).GT.NPIN) THEN
				NPROP(IP)=NPROP(IP)+PROP(NE+IW)/NN
				!AreaVol(IP)=AreaVol(IP)+WArea(IW)/NN
			END IF
		  END DO
	END DO
		

	DO IE=1,NE
		  NN=0
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
		IF(NVIndex(IP).LE.NPIN.AND.NVIndex(IP).GT.0) THEN
			  NN=NN+1
			END IF
		  END DO
		  NumGPTran=NumGPTran+NN
		  
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
		IF(NVIndex(IP).LE.NPIN.AND.NVIndex(IP).GT.0) THEN
			  NPROP(IP)=NPROP(IP)+PROP(IE)/NN
			  !AreaVol(IP)=AreaVol(IP)+Volume(IE)/NN
			END IF
		  END DO

	END DO

	

!AreaVol,PRATIO is used to give real heat flux divergence and heat flux
	VWK1(:)=0.d0
	DO IW = 1, NW
		DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			VWK1(IP)=1.0	!mark it as boundary nodes
			PRATIO(IP)=PRATIO(IP)+PROP(NE+IW)/WElem(1,IW)
			AreaVol(IP)=AreaVol(IP)+WArea(IW)/WElem(1,IW)
		END DO
	END DO
	DO IE=1,NE
		DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
			IF(VWK1(IP).EQ.0.0d0) THEN
		  PRATIO(IP)=PRATIO(IP)+PROP(IE)/NNPE(ElemKind(IE))
		  AreaVol(IP)=AreaVol(IP)+Volume(IE)/NNPE(ElemKind(IE))
			END IF
		END DO
	END DO



!RaNeb is to transfer heat flux divergence from internal nodes to the nodes in surface 
!	(transfering imaginary temperatures from interface(import)_nodes to export_nodes
!	for data communication needs no RaNeb)
	
	DO IE=1,NE
	DO J=1,NNPE(ElemKind(IE))
	IP=ELEM(J,IE)
        RaNeb(IP)=RaNeb(IP)+PROP(IE)/NNPE(ElemKind(IE))
	END DO
	END DO




C
C   2	. Node To Group data transfer index and coefficient
C

	
	ALLOCATE(COEF(NumGPTran),NIndexToGP(2,NumGPTran))


	!surface node
	NumGPTran=0
	DO IW = 1, NW

		  IGP=VIndex(NE+IW)
		  IF(IGP.LE.0) CYCLE

		  NN=0
		  DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			IF(NVIndex(IP).GT.NPIN) THEN
				NN=NN+1
			END IF
		  END DO
		  
		  JP=0
		  DO J=2,WElem(1,IW)+1
			IP=WElem(J,IW)
			IF(NVIndex(IP).GT.NPIN) THEN
				JP=JP+1
				NIndexToGP(1,NumGPTran+JP)=IP
				NIndexToGP(2,NumGPTran+JP)=IGP
			COEF(NumGPTran+JP)=PROP(NE+IW)/NN/GpPROP(IGP)
			END IF
		  END DO
		  NumGPTran=NumGPTran+NN
	END DO
		
	!internal node
	DO IE=1,NE
		  
		  IGP=VIndex(IE)
		  IF(IGP.LE.0) CYCLE
		  
		  NN=0
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
		IF(NVIndex(IP).LE.NPIN.AND.NVIndex(IP).GT.0) THEN
			  NN=NN+1
			END IF
		  END DO
		  
		  JP=0
		  DO J=1,NNPE(ElemKind(IE))
			IP=ELEM(J,IE)
		IF(NVIndex(IP).LE.NPIN.AND.NVIndex(IP).GT.0) THEN
				JP=JP+1
				IGP=VIndex(IE)
				NIndexToGP(1,NumGPTran+JP)=IP
				NIndexToGP(2,NumGPTran+JP)=IGP
			COEF(NumGPTran+JP)=PROP(IE)/NN/GpPROP(IGP)
			END IF
		  END DO
		  
		  NumGPTran=NumGPTran+NN
	END DO


	!check conservation
	VWK1(:)=0.d0
	DO I=1,NumGPTran
		IGP=NIndexToGP(2,I)
		VWK1(IGP)=VWK1(IGP)+COEF(I)
	END DO

	DO IGP=1,NRD
		IF(abs(VWK1(IGP)-1.0).GT.1.0e-3) 
     &	WRITE(IUT0,*) 'Coef violated conservation', IGP,VWK1(IGP)
	END DO


C
C   3	. post treatment
C

	!
	DO I=1,NP
		IF(NVIndex(I).GT.0) THEN
			PRATIO(I)=PRATIO(I)/NPROP(I)
		ELSE
			PRATIO(I)=0.d0	!!!	No participating nodes
		END IF
	END DO
	
	DEALLOCATE(PROP)
	ALLOCATE(PROP(NP))
	PROP(:)=NPROP(:)

	DEALLOCATE(NVIndex,VWK1,NPROP)

	RETURN
!---------> end of NodeToGPIndex()

	END






C======================================================================C	
	SUBROUTINE RDVElemToNode(IUT0,MM,ME,MP,NP,NE,NW,NWT,N2D,
     $  N3D,NKind,
     &	NDOF,NRD,NRD0,NBOUN,NMAT,MAT_NO,ELEM,MatID,XYZ,
     &	IFlagGroup,ICPU)
C======================================================================C
	

        USE MODULE_RADSXF, only :   CoefToRD,IndexToRD,
     &                              NumTran,NumTranIn
	USE MODULE_RADSXF, only :   RaToExpt,IndexToExpt,
     &                              ExptID,NDatExpt
	USE MODULE_RADSXF, only :   PROP,VIndex,WArea,
     &                              Volume,WElem,WNEB
	USE MODULE_RADSXF, only :   ElemKind,NNPE
	USE MODULE_RADSXF, only :   AreaVol,MatBounEnd,NPALL
	use module_radsxf, only :   rdindex,rdvalue

	USE MODULE_RADGROUP, ONLY : GpProp

	USE MODULE_RADWORK, only :  WK1
	USE MODULE_RADWORK, only :  VWK1,VWK4
	USE MODULE_RADWORK, only :  NVIndex => WK2
	USE MODULE_RADWORK, only :  NPROP   => VWK3

	
!----- only for benchmark
	use module_radsxf, only :   NWDat,WDAT,WPNUM
	
	
	IMPLICIT REAL*8(A-H,O-Z)
	
	INTEGER ELEM(N3D,ME),MATID(NE),MAT_NO(-100:100)
	REAL*8	XYZ(NDOF,MP)
	
	CHARACTER*20 StackOpt,StackOpt2
	CHARACTER*80 tmpchar
	


	ALLOCATE(MatBounEnd(0:NMAT+NBOUN),NVIndex(3*NP),
     &	 WK1(MAX(NP,NRD0)),VWK1(MAX(NP,NRD0)),VWK4(MAX(NP,NRD0)),
     &		 stat=ierr)

	
C
C	1. COUNT NPALL (internal nodes and surface nodes, 
C		  list by different materials and boundaries) 

	NPALL=0
	NVIndex    =0
	MatBounEnd=0
	minv=np
	maxv=0

	DO IM=1,NMAT
	  
	  WK1=0
	  DO IE=1,NE
	    
		ID = MAT_NO(MatID(IE))
		IF(ID.NE.IM) CYCLE

		DO J=1,NNPE(ElemKind(IE))
		
			IP=ELEM(J,IE)
			WK1(IP)=1
		END DO
	  END DO

		minv=np
	maxv=0
	  DO I=1,NP
		IF(WK1(I).GT.0) THEN
			NPALL=NPALL+1
			NVIndex(NPALL)=I
	minv=min(minv,I)
	maxv=max(maxv,I)
		END IF
	  END DO
	  MatBounEnd(IM)=NPALL
	
c	write(*,*) 'imat=',im,npall-MatBounEnd(IM-1),minv,maxv

	END DO

	DO IB=1,NBOUN

	  WK1=0
	  DO IW = 1, NW
		IB1=WELEM(WElem(1,IW)+2,IW)

		IF(IB1.EQ.IB) THEN
			DO J=2,WElem(1,IW)+1
				IP=WElem(J,IW)
				WK1(IP)=1
			END DO
		END IF
	  END DO
  
	  DO I=1,NP
		IF(WK1(I).GT.0) THEN
			NPALL=NPALL+1
			NVIndex(NPALL)=I
		END IF
	  END DO
	  MatBounEnd(NMAT+IB)=NPALL
	END DO

	WRITE(IUT0,*) '  NPALL= ', NPALL



C
C  2.  NPROP, AreaVol
C
	ALLOCATE(NPROP(NPALL),AreaVol(NPALL))
	NPROP   =0.d0
	AreaVol =0.d0
	NPALL = 0


	DO IM=1,NMAT
	  
	  WK1=0
	  VWK1=0
	  VWK4=0
	  DO IE=1,NE
	    
		ID = MAT_NO(MatID(IE))
		IF(ID.NE.IM) CYCLE
		
		DO J=1,NNPE(ElemKind(IE))
		
		IP=ELEM(J,IE)
		VWK1(IP)=VWK1(IP)+PROP(IE)/NNPE(ElemKind(IE))
		VWK4(IP)=VWK4(IP)+Volume(IE)/NNPE(ElemKind(IE))
		WK1(IP)=1
		END DO
	  END DO

	  DO I=1,NP
		IF(WK1(I).GT.0) THEN
			NPALL=NPALL+1
			NPROP(NPALL)=VWK1(I)
			AreaVol(NPALL)=VWK4(I)
		END IF
	  END DO
	END DO


	DO IB=1,NBOUN

	  WK1=0
	  VWK1=0
	  VWK4=0
	  DO IW = 1, NW
		
		IB1=WELEM(WElem(1,IW)+2,IW)
		IF(IB1.EQ.IB) THEN
			DO J=2,WElem(1,IW)+1
				IP=WElem(J,IW)
				WK1(IP)=1
			VWK1(IP)=VWK1(IP)+PROP(NE+IW)/WElem(1,IW)
			VWK4(IP)=VWK4(IP)+WArea(IW)/WElem(1,IW)
			END DO
		END IF
	  END DO
  
	  
	  DO I=1,NP
		IF(WK1(I).GT.0) THEN
			NPALL=NPALL+1
			NPROP(NPALL)=VWK1(I)
			AreaVol(NPALL)=VWK4(I)
		END IF
	  END DO
	END DO
	



C
C   3.  Value Index from node to Radiation Equation
C
	NumTran=0

	!internal node
	DO IE=1,NE
		IRD=VIndex(IE)
		IF(IRD.GT.0) NumTran=NumTran+NNPE(ElemKind(IE))
	END DO

	DO IW = 1, NW
		IRD=VIndex(NE+IW)
		IF(IRD.GT.0) NumTran=NumTran+WElem(1,IW)

	END DO
	
	NumTranIn = 0
	DO IW = 1, NW
		IRD=VIndex(NE+IW)
		IB = WElem(WElem(1,IW)+2,IW)
	IF(IRD.GT.0.AND.IB.EQ.NBOUN) NumTranIN=NumTranIN+WElem(1,IW)

	END DO
	NumTranIn=NumTran-NumTranIN

	WRITE(IUT0,'(a,I8,a,I8)') '  NumTran= ', NumTran,
     &			  ', NumTranIn= ', NumTranIN


	
	ALLOCATE(CoefToRD(NumTran),IndexToRD(2,NumTran))
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF(IFlagGroup.NE.1) THEN
	
	  !Internal Node
	  NumTran=0
	  DO IM=1,NMAT
	  
		WK1=0
		DO I=MatBounEnd(IM-1)+1,MatBounEnd(IM)
			IP=NVIndex(I)
			WK1(IP)=I
		END DO
		  
	  
		DO IE=1,NE

			IF(VIndex(IE).LE.0) CYCLE
			ID = MAT_NO(MatID(IE))
			IF(ID.NE.IM) CYCLE
						
			
			JP=0
			DO J=1,NNPE(ElemKind(IE))
				IP=ELEM(J,IE)
				JP=JP+1
				IRD=VIndex(IE)
				IndexToRD(1,NumTran+JP)=WK1(IP)
				IndexToRD(2,NumTran+JP)=IRD
			CoefToRD(NumTran+JP)=1.0/NNPE(ElemKind(IE))
			END DO
					
			NumTran=NumTran+NNPE(ElemKind(IE))
		END DO

	  END DO

  	  !surface node
	  DO IB=1,NBOUN
	  
		WK1=0
		DO I=MatBounEnd(NMAT+IB-1)+1,MatBounEnd(NMAT+IB)
			IP=NVIndex(I)
			WK1(IP)=I
		END DO
	  
		DO IW=1,NW

			IF(VIndex(NE+IW).LE.0) CYCLE
			IF(WElem(WElem(1,IW)+2,IW).NE.IB) CYCLE
		
			JP=0
			DO J=1,WElem(1,IW)
				IP=WElem(J+1,IW)
				JP=JP+1
				IRD=VIndex(NE+IW)
				IndexToRD(1,NumTran+JP)=WK1(IP)
				IndexToRD(2,NumTran+JP)=IRD
				CoefToRD(NumTran+JP)=1.0/WElem(1,IW)
			END DO
					
			NumTran=NumTran+WElem(1,IW)
		END DO

	  END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	ELSE

	  !Internal Node
	  NumTran=0
	  DO IM=1,NMAT
	  
		WK1=0
		DO I=MatBounEnd(IM-1)+1,MatBounEnd(IM)
			IP=NVIndex(I)
			WK1(IP)=I
		END DO
		  
	  
		DO IE=1,NE

			IRD = VIndex(IE)
			IF(IRD.LE.0) CYCLE
			
			ID = MAT_NO(MatID(IE))
			IF(ID.NE.IM) CYCLE
						
			aaa = PROP(IE)/NNPE(ElemKind(IE))/GpProp(IRD)
			JP=0
			DO J=1,NNPE(ElemKind(IE))
				IP=ELEM(J,IE)
				JP=JP+1
				IndexToRD(1,NumTran+JP)=WK1(IP)
				IndexToRD(2,NumTran+JP)=IRD
				CoefToRD(NumTran+JP)=aaa
			END DO
					
			NumTran=NumTran+NNPE(ElemKind(IE))
		END DO

	  END DO

  	  !surface node
	  DO IB=1,NBOUN
	  
		WK1=0
		DO I=MatBounEnd(NMAT+IB-1)+1,MatBounEnd(NMAT+IB)
			IP=NVIndex(I)
			WK1(IP)=I
		END DO
	  
		DO IW=1,NW

			IRD=VIndex(NE+IW)
			IF(IRD.LE.0) CYCLE
			IF(WElem(WElem(1,IW)+2,IW).NE.IB) CYCLE
		
			aaa = PROP(NE+IW)/WElem(1,IW)/GpProp(IRD)
			JP=0
			DO J=1,WElem(1,IW)
				IP=WElem(J+1,IW)
				JP=JP+1
				IndexToRD(1,NumTran+JP)=WK1(IP)
				IndexToRD(2,NumTran+JP)=IRD
				CoefToRD(NumTran+JP)=aaa
			END DO
					
			NumTran=NumTran+WElem(1,IW)
		END DO

	  END DO
		
		
	END IF



C
C   4	. conservation check
C

	VWK1(:)=0.d0
	DO I=1,NumTran
		IRD=IndexToRD(2,I)
		VWK1(IRD)=VWK1(IRD)+CoefToRD(I)
	END DO

	DO IRD=1,NRD
		IF(abs(VWK1(IRD)-1.0).GT.1.0e-3) 
     &	WRITE(IUT0,*) 'Coef violated conservation', IRD,VWK1(IRD)
	END DO

	

C
C   5. Data for communication in parallel computation
C
	NDatExpt=0
	IF(ICPU.LT.0) GOTO 8888

	WK1 = 0
	DO IW=1,NW
		IRD = VIndex(NE+IW)
		IF(IRD.LE.0) CYCLE
		IF(WElem(WElem(1,IW)+2,IW).EQ.NBOUN) THEN
			WK1(IRD) = 1
			IE = WNEB(IW)
			IF(IE.LE.0) STOP 'strange neighbor id'
			IRD = VIndex(IE)
			IF(IRD.GT.0) WK1(IRD)=1
		END IF
	END DO

	VWK4 = 0.d0
	DO I=1,NRD
		IF(WK1(I).EQ.1)  THEN
		DO J=1,I-1
		ID = RDIndex(I)+J
		IF(WK1(J).EQ.0) VWK4(I) = VWK4(I) - RDValue(ID)
		END DO
		DO J=I+1,NRD
		ID = RDIndex(J)+I
		IF(WK1(J).EQ.0) VWK4(I) = VWK4(I) - RDValue(ID)
		END DO
		END IF
	END DO

	
	WK1 = 0
	DO IW=1,NW

		IF(VIndex(NE+IW).LE.0) CYCLE
		
		IF(WElem(WElem(1,IW)+2,IW).EQ.NBOUN) THEN	!domain interface
		
			IE = WNEB(IW)
			IF(IE.LE.0) STOP 'strange neighbor id'

			DO J=1,WElem(1,IW)
				IP=WElem(J+1,IW)
				WK1(IP)=-1
			END DO
		ELSE	!real boundary
			DO J=1,WElem(1,IW)
				IP=WElem(J+1,IW)
				IF(WK1(IP).EQ.0) WK1(IP)=1
			END DO
		END IF
	END DO

c	DO IE=1,NE
c		IF(VIndex(NE).LE.0) CYCLE
c
c		NN = 0
c		Mark = 0
c		DO J=1,NNPE(ElemKind(IE))
c			IP = ELEM(J,IE)
c			IF(WK1(IP).LT.0) THEN
c				Mark=1
c			ELSE IF(WK1(IP).EQ.0) THEN
c				NN=NN+1
c			END IF
c		END DO
c		IF(Mark.EQ.1) NDatExpt=NDatExpt+NN
c	END DO
	
	
	DO IW=1,NW

		IF(VIndex(NE+IW).LE.0) CYCLE
		IF(WElem(WElem(1,IW)+2,IW).LT.NBOUN) CYCLE
		
		IE = WNEB(IW)
		
		NN = 0
		DO J=1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			IF(WK1(IP).EQ.0) NN=NN+1
		END DO
		NDatExpt=NDatExpt+NN
		IF(VIndex(IE).GT.0) NDatExpt=NDatExpt+NN
	END DO
	

	WRITE(IUT0,*) '  NDatExpt = ', NDatExpt


	!!!!!!!!!!!!!!!!!!!!!!
	ALLOCATE(ExptID(NRD))
	ALLOCATE(RaToExpt(NDatExpt),IndexToExpt(2,NDatExpt))
	
	
	NDatExpt=0
	ExptID = 0
	VWK1 = 0.d0

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	IF(IFlagGroup.NE.1) THEN

c	  DO IE=1,NE
c		IF(VIndex(IE).LE.0) CYCLE
c		
c		NN = 0
c		Mark = 0
c		DO J=1,NNPE(ElemKind(IE))
c			IP = ELEM(J,IE)
c			IF(WK1(IP).LT.0) THEN
c				Mark=1
c			ELSE IF(WK1(IP).EQ.0) THEN
c				NN=NN+1
c			END IF
c		END DO
c		IF(Mark.EQ.0) CYCLE
c
c		IRD = VIndex(IE)
c		ExptID(IRD) = 1
c		JP=0
c		DO J=1,NNPE(ElemKind(IE))
c			IP = ELEM(J,IE)
c			IF(WK1(IP).EQ.0) THEN
c				JP=JP+1
c				IndexToExpt(1,NDatExpt+JP)=IP
c				IndexToExpt(2,NDatExpt+JP)=IRD
c				RaToExpt(NDatExpt+JP)=1.0/NN
c				VWK1(IP) = VWK1(IP)+PROP(IE)/NN
c			END IF
c		END DO
c		NDatExpt = NDatExpt +NN
c		
c	  END DO

	  NDatExpt1=NDatExpt
	  
	  DO IW=1,NW

		IF(VIndex(NE+IW).LE.0) CYCLE
		IF(WElem(WElem(1,IW)+2,IW).LT.NBOUN) CYCLE
		
		IE = WNEB(IW)

		NN = 0
		DO J=1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			IF(WK1(IP).EQ.0) NN=NN+1
		END DO
		
		JP=0
		IRD = VIndex(NE+IW)
		ExptID(IRD) = 2
		DO J=1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			IF(WK1(IP).EQ.0) THEN
				JP=JP+1
				IndexToExpt(1,NDatExpt+JP)=IP
				IndexToExpt(2,NDatExpt+JP)=IRD
				RaToExpt(NDatExpt+JP)=1.0/NN
				VWK1(IP) = VWK1(IP)+VWK4(IRD)/NN
			END IF
		END DO
		NDatExpt = NDatExpt +NN

		IRD = VIndex(IE)
		IF(IRD.LE.0) CYCLE
		ExptID(IRD) = 1
		JP=0
		DO J=1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			IF(WK1(IP).EQ.0) THEN
				JP=JP+1
				IndexToExpt(1,NDatExpt+JP)=IP
				IndexToExpt(2,NDatExpt+JP)=IRD
				RaToExpt(NDatExpt+JP)=1.0/NN
				VWK1(IP) = VWK1(IP)+VWK4(IRD)/NN
			END IF
		END DO
		NDatExpt = NDatExpt +NN

	  END DO
	  
	  DO I=1,NDatExpt
		 IP=IndexToExpt(1,I)
		 RaToExpt(I)=RaToExpt(I)/VWK1(IP)
	
	  END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ELSE

c	  DO IE=1,NE
c		IF(VIndex(IE).LE.0) CYCLE
c		
c		NN = 0
c		Mark = 0
c		DO J=1,NNPE(ElemKind(IE))
c			IP = ELEM(J,IE)
c			IF(WK1(IP).LT.0) THEN
c				Mark=1
c			ELSE IF(WK1(IP).EQ.0) THEN
c				NN=NN+1
c			END IF
c		END DO
c		IF(Mark.EQ.0) CYCLE
c
c		IRD = VIndex(IE)
c		aaa = PROP(IE)/NN/GpProp(IRD)
c		ExptID(IRD) = 1
c		JP=0
c		DO J=1,NNPE(ElemKind(IE))
c			IP = ELEM(J,IE)
c			IF(WK1(IP).EQ.0) THEN
c				JP=JP+1
c				IndexToExpt(1,NDatExpt+JP)=IP
c				IndexToExpt(2,NDatExpt+JP)=IRD
c				RaToExpt(NDatExpt+JP)=aaa
c				VWK1(IP) = VWK1(IP)+PROP(IE)/NN
c			END IF
c		END DO
c		NDatExpt = NDatExpt +NN
c
c	  END DO
	  
	  
	  DO IW=1,NW

		IF(VIndex(NE+IW).LE.0) CYCLE
		IF(WElem(WElem(1,IW)+2,IW).LT.NBOUN) CYCLE
		
		IE = WNEB(IW)

		NN = 0
		DO J=1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			IF(WK1(IP).EQ.0) NN=NN+1
		END DO
		
		JP=0
		IRD = VIndex(NE+IW)
		aaa = PROP(NE+IW)/NN/GpProp(IRD)
		ExptID(IRD) = 2
		DO J=1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			IF(WK1(IP).EQ.0) THEN
				JP=JP+1
				IndexToExpt(1,NDatExpt+JP)=IP
				IndexToExpt(2,NDatExpt+JP)=IRD
				RaToExpt(NDatExpt+JP)=aaa
				VWK1(IP) = VWK1(IP)+VWK4(IRD)/NN
			END IF
		END DO
		NDatExpt = NDatExpt +NN

		IRD = VIndex(IE)
		IF(IRD.LE.0) CYCLE
		ExptID(IRD) = 1
		JP=0
		aaa = PROP(IE)/NN/GpProp(IRD)
		DO J=1,NNPE(ElemKind(IE))
			IP = ELEM(J,IE)
			IF(WK1(IP).EQ.0) THEN
				JP=JP+1
				IndexToExpt(1,NDatExpt+JP)=IP
				IndexToExpt(2,NDatExpt+JP)=IRD
				RaToExpt(NDatExpt+JP)=aaa
				VWK1(IP) = VWK1(IP)+VWK4(IRD)/NN
			END IF
		END DO
		NDatExpt = NDatExpt +NN


	  END DO
	  
	  DO I=1,NDatExpt
		 IP=IndexToExpt(1,I)
		 RaToExpt(I)=RaToExpt(I)/VWK1(IP)
	  END DO
	  
	END IF


8888	CONTINUE

	GOTO 9000
	
C----- only for benchmark

	open(1,file='hotspec.txt')
	read(1,*) tmpchar
	read(1,*) NWDAT
	read(1,*) tmpchar
	read(1,*) lenth
	read(1,*) tmpchar
	read(1,*) IB
	read(1,*) tmpchar
	read(1,*) yz0
	close(1)
	x0=-0.5

	ALLOCATE(WDAT(NPALL,4),WPNUM(NWDAT,4))
	
	
	WDAT = 0
      dx =1.0/(NWDat-1.0) 
	dyz=lenth/(NWDat-1.0)
      WPNUM=0

      DO I=MatBounEND(NMAT+IB-1)+1,MatBounEND(NMAT+IB)
       IP = NVIndex(I)
       IF(ABS(XYZ(1,IP)-x0).LT.dx*0.25) THEN
          j=0
          do k=1,NWDat
            a1 = x0-dx*0.5 + dx*(k-1)
            a2 = x0+dx*0.5 + dx*(k-1)
            if(XYZ(2,IP).gt.a1.and.XYZ(2,IP).lt.a2) then
              j=k
              exit
            end if
          end do
          
		if(j.EQ.0.OR.j.GT.nwdat) stop 'error j index'
          IF(ABS(XYZ(3,IP)).LE.1.01*dyz) THEN
			WDAT(I,1)=J
			WPNUM(J,1)=WPNUM(J,1)+1
		END IF

		j=0
		do k=1,NWDat
            a1 = yz0-dyz*0.5 + dyz*(k-1)
            a2 = yz0+dyz*0.5 + dyz*(k-1)
            if(XYZ(3,IP).gt.a1.and.XYZ(3,IP).lt.a2) then
           
              j=k
              exit
            end if
          end do
          
		if(j.EQ.0.OR.j.GT.nwdat) stop 'error j index'
		WDAT(I,3)=J
		WPNUM(J,3)=WPNUM(J,3)+1
		
       END IF
      END DO

	DO I= 1,MatBounEND(1)
	  IP = NVIndex(I)
	  IF(ABS(XYZ(1,IP)-x0-dx).LT.dx*0.25) THEN
          j=0
          do k=1,NWDat
            a1 = x0-dx*0.5 + dx*(k-1)
            a2 = x0+dx*0.5 + dx*(k-1)
            if(XYZ(2,IP).gt.a1.and.XYZ(2,IP).lt.a2) then
              j=k
              exit
            end if
          end do
		if(j.EQ.0.OR.j.GT.nwdat)  stop 'error j index'
		IF(ABS(XYZ(3,IP)).LE.1.01*dyz) THEN
			WDAT(I,2)=J
			WPNUM(J,2)=WPNUM(J,2)+1
		END IF
		
		j=0
          do k=1,NWDat
            a1 = yz0-dyz*0.5 + dyz*(k-1)
            a2 = yz0+dyz*0.5 + dyz*(k-1)
            if(XYZ(3,IP).gt.a1.and.XYZ(3,IP).lt.a2) then
              j=k
              exit
            end if
          end do
		if(j.EQ.0.OR.j.GT.nwdat) stop 'error j index'
		WDAT(I,4)=J
		WPNUM(J,4)=WPNUM(J,4)+1
	  END IF
	END DO


	WRITE(*,*) '  HotWall groups ::'
	DO I=1,NWDAT
	WRITE(*,'(5I6)')I,WPNum(I,1),WPNum(I,2),WPNum(I,3),WPNum(I,4)
	END DO



9000	DEALLOCATE(VIndex,PROP)
	ALLOCATE(VIndex(NPALL),PROP(NPALL))

	PROP(1:NPALL) = NPROP(1:NPALL)
	VIndex(1:NPALL) = NVIndex(1:NPALL)





	DEALLOCATE(WK1,VWK1,VWK4,NVIndex,NPROP)




	RETURN
!---------> end of RDVElemToNode_new()

	END






	
C======================================================================C
	SUBROUTINE UserDefRadPro(x,y,z,gab,pab,sca)
C======================================================================C
	IMPLICIT REAL*8(A-H,O-Z)

	rr = SQRT((x+4.5)**2+(y-0.0)**2+(z-0.0)**2)
	  	
	IF(rr.LE.2.50) THEN

		gab = 0.2
		pab = 0.3*(2.5-rr)/2.50
		!sca=0.1
	ELSE
		gab = 0.1
	END IF
	
	RETURN

!----->

	END
