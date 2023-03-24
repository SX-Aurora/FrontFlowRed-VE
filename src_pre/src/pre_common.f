  

     	TYPE WALLBOUN
		INTEGER::	WHomFlag,WVarFlag
		INTEGER::	WTYPE,BTYPE		
						
						!WType = 1	:		gray wall (black wall included)
						!		 2  :		mirror-reflect
						!		 3  :		user-define	
						
						!BType = 1	:		specified T
						!		 2  :		unknown T
						!		 3  :		convection, unknown T
		REAL*8::	T,GrayDeg,q
	END TYPE


	TYPE GASPROP
		INTEGER::	HaveGasFlag,GHomFlag,GIsoFlag,GVarFlag
		INTEGER::	GTYPE
						!GType = 1	:		gray gas
						!		 2  :		non-gray gas, user-define absorption model
		REAL*8::	AbsorbCoef
	END TYPE


	TYPE PTCPROP
		INTEGER::	HavePtcFlag,PHomFlag,PIsoFlag,PVarFlag
		INTEGER::	PTYPE
						!PType = 1	:		isotropic scatter
						!		 2  :		unisotropic scatter
		REAL*8::	AbsorbCoef,ScatterCoef
		REAL*8::	beta
						!for Legendre Type unisotropic scatter, 
						!beta = 0, is isotropic
	END TYPE

	
					!HaveGasFlag:		Have Gas?
					!HavePtcFlag:		Have Paticles?
					!HaveGasFlag:		Have Wall?
					!GHomFlag,GIsoFlag,GVarFlag,GGrayFlag
					!Is Gas:
					!	homogeneous,isotropic,variable with temperature and gray-gas?
					!P,W ...


	TYPE FILENAMES
		CHARACTER*80::	FileMesh, FileProp, FileInit, 
     &					FileRes,FileConfig
						!FileMesh:		Mesh File Name
						!FileProp:		Property File Name (only need for inhomogeneous media or wall)
						!FileInit:		Initial Field File Name (only need for)
						!FileRes:		For result output
						!FileConfig:	For output Configuration-factor (/ Exchange-Area)
		CHARACTER*20::	MeshFormat		! FV  / GF	/ADV  /SELF
	END TYPE


	TYPE GSOLVER
		
		INTEGER::		EquiFlag
						!EquiFlag = 1, For equilibrium condition, solve the gas temperature
						

		CHARACTER*20::	Model
					!ModelOpt =			'MC'	:	Monte-Carlo Method
					!					'DOM'	:	Discrete Ordinates Method
					!					'FVM'	:	Finite Volume Method
					!					'ZONE'	:	Zone Method
					!					'NET'	:	Net-Radiation Method
					
		CHARACTER*20::	Algorithum
					!!!
					!AlgorithmOpt = 'SOR'	:	Successive Overrelaxation Method 
					!			    'CG'    :	BI-CGSTAB method with preconditioning
					!				'GS'	:	Gaussian Elimination Method
		REAL*8::	aeps,reps,feps,vfeps
					!aeps:			tolerance for average residence of equation T
					!reps:			tolerance for relative residence of equation T
					!feps:			tolerance for residence of exchange-area equation
					!vfeps:			tolerance for error of exchange-area
					
		REAL*8::    RLX
					!RLX:			relaxition in SOR solver

		REAL*8::	VMIN
					!a value for narrowing the matrix band


	END TYPE




	TYPE MONTECARLO
		CHARACTER*20::	RaySpecOpt,EmitType
						!RaySpecOPt:	ray number specify method
						!				'AUTO' / 'SELF'
						!EmitType:		'GAS' / 'WALL' / 'BOTH'
						
		INTEGER::		NRAY
		
		REAL*8::		CONFIDENCE
						!For estimation of the ray-number variance, give the confidence for error approximation
						!the distribution of propability is a canonical one
		REAL*8::		CoefErr
						!average error of ray-number absorption coefficient

	END TYPE


	TYPE ZONEMETHOD
	
		REAL*8::		XiaoJuli
						!For refined calculation of the close elements

	END TYPE



	TYPE CPUTIMES

		REAL	TIME_BEGAN, TIME_END
		REAL	TIME_TRACE1,TIME_TRACE2
		REAL	TIME_SLV1,TIME_SLV2

	END TYPE



	TYPE MISSION

		INTEGER IRestaRead
							!IRstRead =0:	No read
							!		  =1:	Read Restart Data
		INTEGER IRestaWrite
							!IRstWrite=0:	No Write
							!		  =1:	Write Restart Data
		INTEGER JobType
							!		  =1:	Make Exchange-Area Only
							!		  =2:	Calculate temperature and heat flux
		INTEGER HaveInitFlag
							!		  =0:	No Init Field
							!		  =1:	HaveInitField

	END TYPE


	TYPE TEA
		REAL*8::		teaeps
						!residual for Total-Exchange-Area calculation				
			
		INTEGER::		ModTea
						!modtea=1	:	cal. Total-Exchange-Area directly
						!modtea=2	:	cal. Direct-Exchange-Area by MC and solve Total-Exchange-Area by equation

		INTEGER::		MethodCalTea
						!MethodCalTea=1	:		Hottel's traditional method,  
						!			 =2	:		Jiang's method,  

		INTEGER::		IFlagTEA,IModHF,NDIV
						!IFlagTEA = 1,	flag for calculation of the Total-Exchange-Area
						!IModHF	=1,	calculate the heat flux directly
						!       =2, calculate the heat flux with scattering/reflection process.
						!NDIV  :	Number of division of the solid angle

	END TYPE TEA
