
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!---- programs are developed for radiation, Jiang Yuyan, 2005/5-12
!
!  contral_panel ::
!      subroutine advance_rad()
!      subroutine advance_rad_hpc()
!	 subroutine pre_part_radout()
!	 subroutine radsxf()







!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine advance_rad
     & (mvrtx,mcell,mface,nvrtx,ncell,nface,cord,lacell,lvcell,
     &  lbface,lvface)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_partitioner
      use module_rad
      use module_boundary,only:  nbcnd,ivbcnd,lvbcnd,boundIDMap,
     &                           kxnone,kdbcnd,kdprdc,kdsymm,kdintr,
     &                           kdilet,kdolet,kdtchi,kdfire,kdpres,
     &				   kdstag,kdbuff,kdsld,NBOUND,numbname,
     &                           IFFACE,NFBOUN,SFBOUN,boundName,idis,
     &                          radprop,radwalltype,prdcAxis,prdcOffset
      use module_material, only:  radmat,radfludflag,nflud

     
      !/////////
	USE MODULE_RADSXF, only:	NodeBounID
      
      
c      integer,parameter :: kdprdc=-1     !=> dis
c      integer,parameter :: kdsymm=-2
c      integer,parameter :: kdintr=-3     !=> dis
c      integer,parameter :: kdilet=-4
c      integer,parameter :: kdolet=-5
c      integer,parameter :: kdtchi=-6     !=> dis
c      integer,parameter :: kdfire=-7
c      integer,parameter :: kdpres=-8
c      integer,parameter :: kdstag=-9
c      integer,parameter :: kdbuff=-10
c      integer,parameter :: kdsld=-11     !=> dis
c      integer,parameter :: kdcvd=-12
     
C	cord(3,mvrtx)	cordinates
C	lvcell(8,mcell)	element 
C	lacell(mcell)	mat_id
C	lvface(4,mface) vertices of face
C	lbface(2,mface) boundary no., 1, 2:period
C	nbcnd	 : no. of boundary conditions
C	CNTIVITY:	connectivity of element cells
C    SFBOUN : Boundary region name read from grid file.
C    NBOUND : Number of boundary regions read from grid file
C    boundIDMap : Indexmap from FFR-Internal boundary region ID <in fflow.ctl> to 
C                 grid file boundary ID
C    numbname : Maximum number of bounary region for one boundary 
C               condition (usually 2)
C    IFFACE : Vertices ID constructing boundary face.
C    NFBOUN : Number of faces in a boundary region read form grid file.
C    kdbcnd :	saves boundary name and Temperature ...
C    
C---------------------------------------------------------------------
     
!
      implicit none

!-------connectivity      
      integer,parameter :: NKIND=4
      integer,parameter :: lvfcel(4,6,4)=reshape( source=
     &  (/1,3,2,0, 2,3,4,0, 3,1,4,0, 4,1,2,0, 0,0,0,0, 0,0,0,0,
     &    1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0, 0,0,0,0,
     &    1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2,0, 4,5,6,0, 0,0,0,0,
     &    1,5,8,4, 2,3,7,6, 1,2,6,5, 3,4,8,7, 1,4,3,2, 5,6,7,8/),
     &  shape=(/4,6,4/) )
!
! --- [dummy arguments]
!
      integer,intent(in) :: mvrtx,mcell,mface
      integer,intent(in) :: nvrtx,ncell,nface
      real*8 ,intent(inout) :: cord(3,mvrtx)
!
      integer,intent(in) :: lacell(  mcell)
      integer,intent(in) :: lvcell(8,mcell)
      integer,intent(in) :: lvface(4,mface)
      integer,intent(in) :: lbface(2,mface)
      

      
!-------- radiation, Jiang
!radmat(:,nflud)
!1: gabsorb, 2: pabsorb, 3: pscatter, 4:scabeta, 5: ptcdiamter
!radfludtype(:,nflud)
!1: gastype, 2: ptctype, 3: ptccal, 4: fgpflag
								
!local arrays
     
      integer :: lvf(NKIND,6,4),RWallType2(nbcnd+1),RWgpflag2(nbcnd+1)
      real*8 :: RWallProp2(nbcnd+1),prdData(nbcnd+1,5)
	integer::	ikind,i,j,k,ib,ibound,idum(10)
				!1:		conterpart name
				!2:		axis
				!3-5:	offset at x,y,z direction
     
c      prdcAxis(nb)=axis
c      prdcOffset(nb,1)=dble(xoffset)
c      prdcOffset(nb,2)=dble(yoffset)
c      prdcOffset(nb,3)=dble(zoffset)    
!
! --- [local entities]
!
     
     
	RWgpflag2(:) = 1
     
      
!
! ------- pre-treat for consistency with radsxf()
!
      do 100 ikind=1,4
      	do 100 j=1,6
      	  do 100 i=1,4
      		lvf(ikind,j,i)=lvfcel(i,j,ikind)
  100  continue 
  
  	  do 110 ib = 1,nbcnd
  	  	ikind = kdbcnd(0,ib)
  	  	RWallType2(ib)=1				!diffuse-reflection
  	  	RWallProp2(ib)=1.0
  	  	prdData(ib,:)= 0
  	  	IF(ikind.EQ.kdsymm) THEN
  	  		RWallProp2(ib) = 0.0
  	  	ELSE IF(ikind.EQ.kdilet.OR.
     &     		ikind.EQ.kdolet.OR.
     &     		ikind.EQ.kdtchi) THEN
  	  		RWallProp2(ib) = 1.0
  	  	ELSE IF(ikind.EQ.kdprdc) THEN	
 !for periodic boundary, give the offset distance and axis
 !prdData(ib,1)= boundIDMap(ib,2)
  	 IF(prdcAxis(ib)(1:1).EQ.'X'.OR.prdcAxis(ib)(1:1).EQ.'x') THEN
		prdData(ib,2)= 1
         ELSE IF(prdcAxis(ib)(1:1).EQ.'Y'.OR.prdcAxis(ib)(1:1).EQ.'y') 
     &  then

				prdData(ib,2)= 2
		 ELSE
				prdData(ib,2)= 3
		 END IF
  	  		prdData(ib,3)= prdcOffset(ib,1)
  	  		prdData(ib,4)= prdcOffset(ib,2)
  	  		prdData(ib,5)= prdcOffset(ib,3)
  	  		RWallProp2(ib) = -1.0		!mark as periodic face
  	  	ELSE
  	  		RWallProp2(ib)= radprop(1,ib)
  	  		RWallType2(ib)= radwalltype(1,ib)
			RWgpflag2(ib)=radwalltype(2,ib)
		
  	  	END IF
  110 continue
	
	!only for partition-domain-interface, here they are dummys, no use
	RWallProp2(nbcnd+1) = 1.0
	RWallType2(nbcnd+1) = 1
	RWgpflag2(nbcnd+1)  = 0
	prdData(nbcnd+1,:) = 0


	!Get NodeBounID
	ALLOCATE(NodeBounID(nvrtx))
     
      CALL radsxf(-1,mvrtx,mcell,nvrtx,ncell,nflud,
     &		cord,lacell,lvcell,lvf,NKIND,
     &  		RWallProp2,RWallType2,RWgpflag2,prdData)
  
      
!
!      deallocate(IFFACE,NFBOUN)
! --- 
!
      return
      end subroutine advance_rad
















!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine advance_rad_hpc
     & (icpu,mvrtx,mcell,nvrtx,ncell,cord,lacell,lvcell)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
	use module_io, only : ifll,ifle
      use module_partitioner
      use module_rad,	  only:    radmodel,radflag
      use module_boundary,only:  nbcnd,ivbcnd,lvbcnd,boundIDMap,
     &                           kxnone,kdbcnd,kdprdc,kdsymm,kdintr,
     &                           kdilet,kdolet,kdtchi,kdfire,kdpres,
     &				   kdstag,kdbuff,kdsld,NBOUND,numbname,
     &                           IFFACE,NFBOUN,SFBOUN,boundName,idis,
     &                         radprop,radwalltype,prdcAxis,prdcOffset
      use module_material, only:  radmat,radfludflag,nflud

      !/////////
	USE module_radsxf,only: NodeBounID
	USE module_radsxf,only: WK7
      
      
c      integer,parameter :: kdprdc=-1     !=> dis
c      integer,parameter :: kdsymm=-2
c      integer,parameter :: kdintr=-3     !=> dis
c      integer,parameter :: kdilet=-4
c      integer,parameter :: kdolet=-5
c      integer,parameter :: kdtchi=-6     !=> dis
c      integer,parameter :: kdfire=-7
c      integer,parameter :: kdpres=-8
c      integer,parameter :: kdstag=-9
c      integer,parameter :: kdbuff=-10
c      integer,parameter :: kdsld=-11     !=> dis
c      integer,parameter :: kdcvd=-12
     
C	cord(3,mvrtx)	cordinates
C	lvcell(8,mcell)	element 
C	lacell(mcell)	mat_id
C	lvface(4,mface) vertices of face
C	lbface(2,mface)
C	nbcnd	 : no. of boundary conditions
C	CNTIVITY:	connectivity of element cells
C    SFBOUN : Boundary region name read from grid file.
C    NBOUND : Number of boundary regions read from grid file
C    boundIDMap : Indexmap from FFR-Internal boundary region ID <in fflow.ctl> to 
C                 grid file boundary ID
C    numbname : Maximum number of bounary region for one boundary 
C               condition (usually 2)
C    IFFACE : Vertices ID constructing boundary face.
C    NFBOUN : Number of faces in a boundary region read form grid file.
C    kdbcnd :	saves boundary name and Temperature ...
C    
C---------------------------------------------------------------------
     
!
      implicit none

!-------connectivity      
      integer,parameter :: NKIND=4
      integer,parameter :: lvfcel(4,6,4)=reshape( source=
     &  (/1,3,2,0, 2,3,4,0, 3,1,4,0, 4,1,2,0, 0,0,0,0, 0,0,0,0,
     &    1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0, 0,0,0,0,
     &    1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2,0, 4,5,6,0, 0,0,0,0,
     &    1,5,8,4, 2,3,7,6, 1,2,6,5, 3,4,8,7, 1,4,3,2, 5,6,7,8/),
     &  shape=(/4,6,4/) )
!
! --- [dummy arguments]
!
      integer,intent(in) :: mvrtx,mcell,icpu
      integer,intent(in) :: nvrtx,ncell
      real*8 ,intent(inout) :: cord(3,mvrtx)
      integer,intent(in) :: lacell(  mcell)
      integer,intent(in) :: lvcell(8,mcell)
      
				
!local arrays
      
	character*80	tmpStr80
      integer :: lvf(NKIND,6,4),RWallType2(nbcnd+1),RWgpflag2(nbcnd+1)
      real*8 ::  RWallProp2(nbcnd+1),prdData(nbcnd+1,5)
	integer::	ikind,i,j,k,ib,ibound,nn1,ivv,nbs
	integer::	idum(4),MM,mface,nface,iface,mb,nb
			!prdData(NBOUND+1,5) -- 1:		conterpart name
			!						2:		axis
			!						3-5:	offset at x,y,z direction
     


	ALLOCATE (NodeBounID(nvrtx),WK7(nvrtx))
	
	NodeBounID(1:nvrtx) = 0
	
!-------------------
! --- Read BC file
!-------------------
      open (31,file=BCout(icpu),status='unknown')
      read (31,*) tmpStr80
      read (31,*) tmpStr80
      read (31,*) k
!
	IF(radmodel(1:3).NE.'FVM') THEN
	  do ib=1,nbcnd
		if(kdbcnd(0,ib).eq.kdprdc) then
		  write(ifle,*) 'Error abort :'
		  write(ifle,*) '  Only <FVM> model can be used for'
	  write(ifle,*) '  Parallel_computation && Periodic_boundary'
		  STOP
		end if
	  end do
	END IF
     
      !write(ifll,*) k
      do 100 mb=1,nbcnd
		read (31,'(a80)') tmpStr80
		!SFBOUN(mb) = trim(adjustl(tmpStr80))
		
		read (31,*) nn1,ivv,nbs
		!write(ifll,*) nn1,ivv,nbs
		DO I=1,(ivv-1)/6+1
		read (31,*) (WK7(k),k=(I-1)*6+1,MIN(i*6,ivv))
		END DO

		DO I=1,ivv
			IF(NodeBounID(WK7(I)).EQ.0) THEN
				NodeBounID(WK7(I))=mb
			ELSE
				NodeBounID(WK7(I))=-100
	!2nd time appearance, this is a corner node, do not use
	!it to specify boundary id
			END IF
		END DO

		if(nbs>0) then
		  do i=1,nbs
		  read(31,*) (idum(k),k=1,4)
		  enddo
		endif

	    if(kdbcnd(0,mb).eq.kdintr.or.kdbcnd(0,mb).eq.kdsld) then
			read (31,'(a80)') tmpStr80
				
			read (31,*) nn1,ivv,nbs
			DO I=1,(ivv-1)/6+1
			read (31,*) (WK7(k),k=(I-1)*6+1,MIN(i*6,ivv))
			END DO

			if(nbs>0) then
			  do i=1,nbs
			  read(31,*) (idum(k),k=1,4)
			  enddo
			endif	
	    end if
 100  continue
      
	close(31)

	DEALLOCATE(WK7)


! ------- pre-treat for consistency with radsxf()
!
      do 110 ikind=1,4
      	do 110 j=1,6
      	  do 110 i=1,4
      		lvf(ikind,j,i)=lvfcel(i,j,ikind)
  110  continue 
  
	

	RWgpflag2(:) = 1

  	do 120 ib = 1,nbcnd
  	  	ikind = kdbcnd(0,ib)
  	  	RWallType2(ib)=1		!diffuse-reflection
  	  	RWallProp2(ib)=1.0
  	  	prdData(ib,:)= 0
	
  	  	IF(ikind.EQ.kdsymm) THEN
  	  		RWallProp2(ib) = 0.0
			RWallType2(ib) = 2		!mirror surface
  	  	ELSE IF(ikind.EQ.kdilet.OR.
     &     		ikind.EQ.kdolet.OR.
     &     		ikind.EQ.kdtchi) THEN
  	  		
			RWallProp2(ib) = 1.0
			
  	  	ELSE IF(ikind.EQ.kdprdc) THEN	
!for periodic boundary, give the offset distance and axis
!prdData(ib,1)=boundIDMap(ib,2)
  	 IF(prdcAxis(ib)(1:1).EQ.'X'.OR.prdcAxis(ib)(1:1).EQ.'x') THEN
				prdData(ib,2)= 1
         ELSE IF(prdcAxis(ib)(1:1).EQ.'Y'.OR.prdcAxis(ib)(1:1).EQ.'y') 
     &	THEN
				prdData(ib,2)= 2
		 ELSE
				prdData(ib,2)= 3
		 END IF
  	  	 prdData(ib,3)= prdcOffset(ib,1)
  	  	 prdData(ib,4)= prdcOffset(ib,2)
  	  	 prdData(ib,5)= prdcOffset(ib,3)

  	  	 RWallProp2(ib) = -1.0		!mark as periodic face
  	  	ELSE
  	  	
			RWallProp2(ib)= radprop(1,ib)
  	  		RWallType2(ib)= radwalltype(1,ib)
			RWgpflag2(ib)= radwalltype(2,ib)
  	  	!write(ifll,*) radprop(1,ib), radwalltype(ib),ib
		END IF
  120 continue
	
	!for partition-domain-interface,treat as black-wall
	RWallProp2(nbcnd+1) = 1.0
	RWallType2(nbcnd+1)= 1
	RWgpflag2(nbcnd+1)= 0			!interface not grouped
	prdData(nbcnd+1,:) = 0

!transfor from cell to vertex

           
      CALL radsxf(icpu,mvrtx,mcell,nvrtx,ncell,nflud,
     &		cord,lacell,lvcell,lvf,NKIND,
     &  		RWallProp2,RWallType2,RWgpflag2,prdData)
      
!
!      deallocate(IFFACE,NFBOUN)
! --- 
!
	
	deallocate(SFBOUN,boundIDMap)

      return
      end subroutine advance_rad_hpc






!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine pre_part_radout (icpu)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- [module arguments]
!
      use module_io,only : ifll,ifle
	use module_partitioner
	
!radiation, appended by Jiang yuyan, 2005/10/28
	use module_radsxf, only  : RDValue,RDIndex,PROP,RDN1,RDN2,RDId,
     &				   VIndex,AreaVol,NodeBounID,
     &				   NRD,MRD,NRD0,MRD0,MID,StackOpt
!
      implicit none
!
! --- [dummy arguments]
	integer,intent(in) :: icpu
!
! --- [local entities]
!
      integer :: i,j,IS,boundID,my_rank,iv,ic,mb
      integer :: ifl,ios=0
      integer :: NBFS,NBOUND_TEMP
      character(len=80),save :: fnam




	ifl=50
      open(ifl,file=RADout(ICPU),		!CONVERT='BIG_ENDIAN',
     &         FORM='unformatted',status='unknown',iostat=ios)
      if(ios/=0) then
        write(ifll,*)'*** Cannot create Rad File:',RADout(ICPU)
        return
      end if
!
      
!appended by Jiang Yuyan for Radiative Heat Transfer
!2005/10/28
	WRITE(ifl) StackOpt
	WRITE(ifl) NRD,NRD0,MRD,MID

	WRITE(ifl) (NodeBounID(I),I=1,NRD0)
	WRITE(ifl) (VIndex(I),I=1,NRD0)
	WRITE(ifl) (AreaVol(I),I=1,NRD)
	WRITE(ifl) (PROP(I),I=1,NRD)
	
	WRITE(ifl) (RDIndex(I),I=1,NRD+1)
	WRITE(ifl) (RDValue(I),I=1,MRD)
	IF(StackOpt(1:4).EQ.'BAND') THEN
		WRITE(ifl) (RDN1(I),I=1,NRD)
	    WRITE(ifl) (RDN2(I),I=1,NRD)
		WRITE(ifl) (RDId(I),I=1,MRD)
	END IF
      close(ifl)
!
!      deallocate(IFFACE,NFBOUN)
! --- 
!
      return
      end subroutine pre_part_radout


C-------------------------------------------------------------------------------
C/
C/	this subroutine inherits from the code RADSXF with part of the functions
C/	 2005/10/25
C/	(RADSXF is a Research_Oriented_Radiation_Code developed by YY.JIANG)
C-------------------------------------------------------------------------------

C======================================================================C
C                                                                      C
C PROGRAM RadXF()														 C
C                                                                      C
C                                       WRITTEN BY YY.JIANG            C
C                              Ver1:    2005/04/11~2005/07/15/         C
C                              Ver2:    2005/07/21~2005/09/01/         C
C                              Ver3:    2005/10/01~2005/11/01/	       C
C                                                                      C
C			                                               C
C PARTS OF THE CODE WERE DEVELOPED BY YY.JIANG BEFORE 2005/03/01       C
C ALL RIGHTS RESERVED, COPYRIGHT(C). UNIVERSITY OF TOKYO, FSIS PROJECT C
C======================================================================C
!
	SUBROUTINE RADSXF(ICPU,MP,ME,NP,NE,nfluid,
     &			XYZ,MatID,ELEM,LOCAL,NKIND,
     &			radwallprop,radwalltype,rwgpflag,prdData)

    !///////// imodhf			
	use module_io, only : iut0 => ifll
	use module_io, only : ifle
	
	use module_rad	!in pre_module_total.f

	USE MODULE_RADSXF
	USE MODULE_RADGROUP
	use module_boundary,only:  NBOUND,IFFACE,
     &         NFBOUN,SFBOUN,boundIDMap
	use module_boundary,only:  nbcnd,kdbcnd,kdprdc,kdintr,kdsld,
     &				 prdcAxis,prdcOffset,idis
	use module_material, only: radmat,radfludflag	
!in radmodule_total.f

	
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER ( N3D  = 8,		    N2D =  4,
     2			NFACE= 6,			NDOF = 3,
     3			FDA  = 0.15,		NSEG = 2,
     5			PAI	= 3.14159265,	SHIGMA = 5.67e-8)

C		
C	----------------------------------------------------------------------------------------- C
C			A GENERAL-PURPOSE CODE															  C
C			DEVELOPED FOR RADIATIVE HEAT TRANSFER COMBINED WITH FLUID FLOW OR HEAT CONDUCTION C
C										--- BY YY JIANG										  C
C									FSIS Proj. RSS21 Proj. IIS,The University of Tokyo		  C
C	----------------------------------------------------------------------------------------- C
C	 THIS IS A PROGRAM FOR THE CALCULATION OF THE RADIATION HEAT TRANSFER COUPLED
C    WITH FLUID FLOW, BOTH THE SHAPE AND THE MESH GRID COULD BE ARBITRARY.
C
C      MAIN FEATURES:
C		(1) MONTE-CARLO METHOD, ZONE / NET-RADIATION METHOD, FINITE VOLUME / DISCRETE ORDINATES METHOD
C		(2) THE GRID MAY BE TETRAHEDRAL OR HEXAHEDRAL, PRISM, PYRAMID ..., OR A HYBRID ONE
C		(3)	BOTH 2-D AND 3-D PROBLEMS CAN BE SOLVED,
C		(4) GAS MAY BE TRANSPARENT OR PARTICIPENT,
C		(5) WALL REFLECTION MAY BE DIFFUSIVE, MIRRORY, OR WHATEVER DEFINED BY USER
C		(6)	PARTICLE RADIATION AND SCATTERING IS INCLUDED
C		(7) THE RADIATIVE PROPERTIES MAY BE OPTIONALLY,
C					ISOTROPIC / UN-ISOTROPIC
C					HOMERGENOURS / NON-HOMERGENOURS
C					CONSTANT / VARIABLE WITH TEMPERATURE
C					GRAY / NON-GRAY (REAL GAS WITH ABSORBING BAND)
C
C	  IT IS EXTREMELY IMPORTANT TO MAKE A REASONABLE CHOICE WITH RESPECT TO
C	THE MODEL FOR THE PROPERTY OF REAL GAS AND/OR PARTICLES. THESE PARTS ARE 
C	ALWAYS OPEN FOR CUSTOM DEVELOPMENTS. 
C				/[ particlemodel.f ]
C	
C	  HOWEVER YOU ARE NOT APPECIATED TO MODIFY THE CODE IN THE FOLLOWING FILES
C	UNLESS YOU ARE QUITE FAMILIAR WITH WHAT IS GOING IN THIS PROGRAM. AN HASTY
C	MODIFICATION MAY HARM BOTH THE LOGIC AND EFFICIENCY OF THE CODE.
C				/[ function.f bucket.f shapefun.f solver.f 
C				   montecarlo.f zone.f fvm.f ]
C
C	  THIS FIGURE SHOWS THE IMPLICIT COORDINATES WE USED, THE VERTICE' CONNECTVITY
C	IN AN ELENMENT MAY BE CLOCKWISE OR ANTI-CLOCKWISE. THE CODE CAN JUDGE THERM CELL
C	BY CELL AUTOMATICCALLY.
C
C											 Y
C											 |
C											 |
C										     |_________________
C										 	/|				  /	
C										   / |		U		 /|
C										  /	 |			    / |
C									     /___|_____________/  |
C									     |	 |		  	   |  |
C									     |	 |		S	   |  |
C									     | W |			   | E|
C									     |	 |			   |  |
C									     |	 |	  N		   |  |
C										 |	 /-------------|--/----------> X
C										 |  /			   | /
C										 | /		D	   |/
C									     |/________________/
C										 /
C										/
C									   /
C									  /
C									 Z
C			
	
	
	!------ FLUID		(-->transfer)
	INTEGER ELEM(N3D,ME),MatID(NE),LOCAL(NKIND,NFACE,N2D)
	REAL*8	XYZ(NDOF,MP),GlobeCTR(NDOF),FanWei(2,NDOF)
	integer :: MAT_NO(-100:100)

	!------ WALL		(-->transfer)
	INTEGER	radwalltype(nbcnd+1),rwgpflag(nbcnd+1),
     &		wsmk(-nbcnd-1:nbcnd+1)
	REAL*8	radwallprop(nbcnd+1),prdData(nbcnd+1,5)
!prdData(nbcnd+1,5) -- 1:		conterpart name
!		2:		axis
!		3-5:	offset at x,y,z direction

	!Private Options
	CHARACTER*20 RHTIsOutward,tmpchar
					!InfoOutOpt,ViewOutOpt,SampleOpt, RestartOutOpt,OutwardViewIsClock :
					!						'Yes' | 'No'
					!"RightHandTurnIsOutward", is necessary for wall face searching
					!MissionOpt-----
					! PASSIVE	: Heat transfer through convection or conduction is dominant
					! POSITIVE  :    ""					radiation is dominant
						
					!!!
					!StackOpt =			'BAND'	:	stack as a 1-D Array with element indeices(M*NNN)
					!								suitable for large-scale mesh with optical-thick gas
					!					'LTRI'	:	stack as a 1-D Array, which is anologous to
					!								a Lower-Trianglar-Array (M*M/2), suitable for 
					!								small-scale mesh or mesh with optical-thin gas
					

	!OTHERS
	REAL*8	V1(3),V2(3),V3(3),V4(3),FaceNode(N2D,3)
	REAL*8	EASUM(0:nbcnd+1,0:nbcnd+1),QSUM(2,nbcnd+1),ddd(3)
	INTEGER NE,NP,NW,NWT,GHomFlag,NUMP(-nbcnd-1:nbcnd+1)
	INTEGER*2 ISTAGE
	INTEGER RevBounMap(NBOUND+1)
	REAL	TIME_TRACE1,TIME_TRACE2,TIME_BEGAN,TIME_END
	
C
     
	DATA  IUT1 / 11 /
	DATA  IUT2 / 12 /
	DATA  IUT3 / 13 /
      

C***********************************************************************************C
C----
C----								1. INITIALIZE
C----
C***********************************************************************************C

			
	CALL CPU_TIME(TIME_BEGAN)
	ISTAGE = 1


		

	WRITE(IUT0,*) '=================================================='
	WRITE(IUT0,*) '-                                                -'
	WRITE(IUT0,*) '-           PROGRAM < RADSXF > START             -'
	WRITE(IUT0,*) '-                                                -'
	WRITE(IUT0,*) '=================================================='	


C==============================================================================C	
C
C	Pretreat Mesh
C
C

	WRITE(IUT0,9900) '<<--',ISTAGE,'-->> PRE-TREAT MESH'
	ISTAGE=ISTAGE+1


	ALLOCATE(FaceNum(NKIND),NNPE(NKIND),NNPS(NKIND,NFACE),
     &		 CNTIVITY(NKIND,NFACE,N2D),ElemKind(NE))
	

	MM=MAX(ME,MP)
	!MPB = MAX(MAX(20*ME,2000000),MM*N3D)

	FaceNum=0
	NNPE=0
	NNPS=0

	FaceNum(1)=4
	FaceNum(2)=5
	FaceNum(3)=5
	FaceNum(4)=6

	NNPE(1)=4
	NNPE(2)=5
	NNPE(3)=6
	NNPE(4)=8

	NNPS(1,1:4)=3
	NNPS(2,1)=4
	NNPS(2,2:5)=3
	NNPS(3,1:3)=4
	NNPS(3,4:5)=3
	NNPS(4,1:6)=4

	DO 10 I=1,N2D
	 DO 10 J=1,NFACE
	  DO 10 K=1,NKIND
		CNTIVITY(K,J,I) = LOCAL(K,J,I)
10	CONTINUE
	 
C
C    boundIDMap : Indexmap from FFR-Internal boundary region ID <in fflow.ctl> to 
C                 grid file boundary ID
C
C		boundIDMap :    fflow.ctl ->	mesh
C		RevBounMap :    mesh	  ->	fflow.ctl 
C
C
C		
C		nbcnd,IFFACE,NFBOUN,SFBOUN		
C					------   The num and turn is the same as input mesh 
C
C		radwalltype,rwgpflag,radwallprop,prdData 
C					------   the turn is as that in <fflow.ctl>
C							 and hence the RevBounMap transfers the Mesh_Boun_No -> FFR_BOUN_No		
C
C
	RevBounMap=0
	DO I=1,nbcnd
		IB = boundIDMap(I,1)		!fflow.ctl -> mesh
		RevBounMap(IB)=I			!mesh -> fflow.ctl
		IB2 = boundIDMap(I,2)
		IF(IB2.GT.0) RevBounMap(IB2) = -I
	END DO
	RevBounMap(NBOUND+1)=nbcnd+1
	!write(*,*) 'NBOUND', nbcnd,nbound

	MarkPeriod=0
	wsmk=1
	DO I=1,nbcnd
		!ib1 = boundIDMap(I,1)
		ib2 = boundIDMap(I,2)
		IF(kdbcnd(0,i).EQ.kdprdc) THEN
		  !RevBounMap(ib2) = -I
		  prdData(I,1)=-I
		  MarkPeriod=1
		END IF
		if(kdbcnd(0,i).EQ.kdsld.or.kdbcnd(0,i).EQ.kdintr) then
			wsmk(-i) = 0
			if(radwallprop(i).LE.0.d0) wsmk(i) = 0
			
		end if
	END DO

	!error abort
	IF(MarkPeriod.EQ.1.AND.RadModel(1:4).EQ.'ZONE') THEN
	 WRITE(IUT0,*) 'Error Abort in RadSxf() '
       WRITE(IUT0,*) 'Zone method not applicable for periodic boundary'
	 WRITE(IUT0,*) 'Please change model to <MC> or <FVM>'
	 STOP
	END IF
	

	!get ELEMKIND
	DO 100 I=1,NE
		INNPE=0
		LOOP_A : DO J=1,8
			IF(ELEM(J,I).GT.0) THEN
				INNPE=INNPE+1
			ELSE
				EXIT LOOP_A
			END IF
		END DO LOOP_A
		
		CHECK_NNPE: SELECT CASE (INNPE)
		CASE (4)
			ElemKind(I) = 1
		CASE (5)
			ElemKind(I) = 2
		CASE (6)
			ElemKind(I) = 3
		CASE (8)
			ElemKind(I) = 4
		END SELECT CHECK_NNPE
	
100	CONTINUE


	!Make WElem
	NUMP(:)=0
		
	IF(ICPU.LT.0)	THEN	!Global Mode, Serial Computation

	  NW = 0
	  MW = 0
	  DO IB = 1,NBOUND
		ibfflow=RevBounMap(IB)
		!IF(ibfflow.EQ.0) CYCLE
		NW=NW+NFBOUN(IB)
	  END DO
	  MW=NW
	  NWT=NW
	  
	  ALLOCATE(WELEM(N2D+2,MW),WNEB(MW),stat=ierr)	  


	  NodeBounID(:)=0
	  IW = 0
	  DO 200 IB = 1,NBOUND			!bound in mesh
		ibfflow=RevBounMap(IB)
		!IF(ibfflow.EQ.0) CYCLE		
					!NOTE::
					!slave pair of interface/slide boundary is not included 
					!(at this time, RevBounMap(IB))
					!but slave pair of periodic boundary is included
		
		DO J=1,NFBOUN(IB)
		  I2D=0
		  IW=IW+1
		  IN1 = 0
		  DO K=1,N2D
			IN2 = IFFACE(K,J,IB)
		    IF(IN2.GT.0.AND.IN2.NE.IN1) THEN
		  		I2D=I2D+1
		  		WElem(K+1,IW)=IN2
				IN1 = IN2

				!IP = WElem(K+1,IW)
c				IF(NodeBounID(IP).EQ.0) THEN
c					NodeBounID(IP)=IB			!1st appearence
c				ELSE IF(NodeBounID(IP).NE.IB) THEN
c					!2nd appearance, and belongs to different boundaries, it is a corner node
c					NodeBounID(IP)=-100
c				END IF
			ELSE
		  		EXIT
		  	END IF
		  END DO
		  WElem(1,IW)=I2D
	
		  IBfflow=RevBounMap(IB)	!mesh -> fflow.ctl
		  WElem(I2D+2,IW)=	IBfflow
		  NUMP(IBfflow)=NUMP(IBfflow)+1
		END DO


200	  CONTINUE

	  !extract all surface patches, append the surfaces that
	  !are not specified as wall
	  !get WNEB
	 
	  CALL  SurfaceExtract(IUT0,XYZ,ELEM,ElemKind,FaceNum,NNPS,
     &			CNTIVITY,NNPE,WELEM,NP,NE,NW,NWT,MP,ME,MW,
     &			NDOF,N2D,N3D,NFACE,NKIND,WNEB,nbcnd)

									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
									!NOTE: nbcnd should not be changed
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  	


	WRITE(IUT0,*) '  SURFACE CELL NUMBER OF EACH BOUNDARY :: '
		
		DO IB=1,nbcnd
		    WRITE(IUT0,'("   <",I3," >  :  ", I5)') IB,
     &				NUMP(IB)
			IF(NUMP(-IB).GT.0) 
     &		WRITE(IUT0,'("           : <pair>",I5)') NUMP(-IB)
		END DO
		
	
	ELSE	!Partition Mode, Parallel Computation

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!--------------------------------------------!
		!periodic boundary is excepted
		!--------------------------------------------!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!extract all surface patches, append the surfaces that
		!are not specified as wall
		
		MW=NE*1.5
		NW=0
		NWT=0
				
		ALLOCATE (WKMD1(N2D+2,MW),WNEB(MW))
	 
				
		WKMD1 = 0

		!get WNEB and appending surface patches
		CALL  SurfaceExtract(IUT0,XYZ,ELEM,
     &			ElemKind,FaceNum,NNPS,CNTIVITY,NNPE,
     &			WKMD1,NP,NE,NW,NWT,MP,ME,MW,
     &			NDOF,N2D,N3D,NFACE,NKIND,WNEB,nbcnd)
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
									!NOTE: nbcnd should not be changed
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		ALLOCATE (WELEM(N2D+2,NWT))
		WElem(:,:)=0
		NW = NWT
          
		iadd = 0
		DO I=1,NWT
			DO J=1,N2D+2
				WElem(J,I)=WKMD1(J,I)
			END DO
				
			Mark=0	
			DO J=1,WElem(1,I)
				IN=WElem(J+1,I)
				
				IF(NodeBounID(IN).EQ.0) THEN
	!contains undefined node, it is an interface patch
					WElem(WElem(1,I)+2,I)=nbcnd+1
					iadd=iadd+1
					Mark=1
					Exit
				END IF
				
			END DO
	
			IF(Mark.EQ.0) THEN	!physical boundary patch
				DO J=1,WElem(1,I)
					IN=WElem(J+1,I)
					IF(NodeBounID(IN).GT.0) THEN 
				WElem(WElem(1,I)+2,I)=NodeBounID(IN)
						EXIT
					END IF
				END DO
				DO J=1,WElem(1,I)
					IN=WElem(J+1,I)
					ID1 = WElem(WElem(1,I)+2,I)
					ID2 = NodeBounID(IN)
				IF(ID2.GE.0.AND.ID1.NE.ID2) THEN 
		  WRITE(IUT0,*) 'Error: boundary id is different!'
		  WRITE(IUT0,*) 'Patch ID=',ID1,'Node ID=',ID2
				  STOP
				END IF
				END DO
			END IF
			IB=WElem(WElem(1,I)+2,I)
			NUMP(IB)=NUMP(IB)+1
	    END DO

		DEALLOCATE (WKMD1)


	WRITE(IUT0,*) '  SURFACE CELL NUMBER OF EACH BOUNDARY :: '
		DO IB=1,nbcnd
		WRITE(IUT0,'("   <",I3," >  :  ", I5)') IB,NUMP(IB) 
		IF(NUMP(-IB).GT.0) 
     &		WRITE(IUT0,'("           : <pair>",I5)') NUMP(-IB)
		END DO
	
		IB=nbcnd+1
		IF(NUMP(IB).GT.0) 
     &	WRITE(IUT0,'("   <hpc_domain_interface>  :  ", I5)') NUMP(IB) 
		
	END IF
      

	
	

C==============================================================================C	
C
C	Allocate working arrays and calculate some parameters
C
C	
	
	!Build Important Parameters
		
			!NW:	number of wall patches
			!NWT:	number of total surface patches
			!NRD0:	All element
	
	NRD0 = NE+NW	!R-E-A-D arrays, for all elements
		
	!allocate arrays
	ALLOCATE (WNV(NDOF,NWT),WArea(NWT),WABCD(4,NWT),WGrayDeg(NWT))
	ALLOCATE (VOLUME(NE),GAbsorb(NE),PAbsorb(NE),PScatter(NE))
	ALLOCATE (VIndex(NRD0),Albedo(NRD0),PROP(NRD0),scabeta(NE))
	ALLOCATE (ETYPE(NRD0))

	
	!get the face equation, normal and tangent vectors of wall faces
	!the normal vector MUST BE OUTWARD,
	!this surbroutine can get the OUTWARD normal vector automatically
	
	CALL WallFaceEquVect(WABCD,WNV,WElem,WArea,XYZ,ELEM,WNEB,
     &		 ElemKind,NNPE,NDOF,MP,ME,NWT,N2D,N3D,NKIND)
	
	Mark=0
	DO IW=1,NW
		
		ib = WELEM(WElem(1,IW)+2,IW)
		if(wsmk(ib).EQ.0) CYCLE
	if(kdbcnd(0,ib).NE.kdsld.AND.kdbcnd(0,ib).NE.kdintr) CYCLE
				
		IE=WNEB(IW)
		IMAT = MatID(IE)
		IF(IMAT.LE.0) THEN
			WNV(1,IW)=-WNV(1,IW)
			WNV(2,IW)=-WNV(2,IW)
			WNV(3,IW)=-WNV(3,IW)
			mark=ib
		END IF
		!if(ib.eq.6) write(*,*) iw,ie,imat

	END DO
		
	!get volumes of element
	CALL ElemVolume(XYZ,ELEM,Volume,N2D,N3D,NDOF,NFACE,NP,NE,MP,ME,
     &			NKIND,ElemKind,FaceNum,NNPE,NNPS,CNTIVITY,
     &				AveVol)


	!get center of the element
	ALLOCATE(CTR(3,NRD0))
	DO 250 IE =1,NE
		IK = ElemKind(IE)
		CTR(:,IE) = 0.0
		DO J=1,NNPE(IK)
			IP = ELEM(J,IE)
			CTR(:,IE) = CTR(:,IE) + XYZ(:,IP)
		END DO
		CTR(:,IE) = CTR(:,IE)/NNPE(IK)
250	CONTINUE
		
	AveWarea=0
	DO 260 IW =1,NW
		IIT = NE+IW
		CTR(:,IIT) = 0.0
		DO J=1,WElem(1,IW)
			IP = WELEM(J+1,IW)
			CTR(:,IIT) = CTR(:,IIT) + XYZ(:,IP)
		END DO
		CTR(:,IIT) = CTR(:,IIT)/WElem(1,IW)
		AveWarea=AveWArea+WArea(IW)
260	CONTINUE
	AveWarea=AveWarea/NW

	
	!get meshrange
	FanWei(1,1)=XYZ(1,1)
	FanWei(1,2)=XYZ(2,1)
	FanWei(1,3)=XYZ(3,1)
	FanWei(2,1)=XYZ(1,1)
	FanWei(2,2)=XYZ(2,1)
	FanWei(2,3)=XYZ(3,1)
      DO 270 IP = 1, NP
		FanWei(2,1) = max(XYZ(1,IP), FanWei(2,1))
		FanWei(2,2) = max(XYZ(2,IP), FanWei(2,2))
		FanWei(2,3) = max(XYZ(3,IP), FanWei(2,3))
		FanWei(1,1) = min(XYZ(1,IP), FanWei(1,1))
		FanWei(1,2) = min(XYZ(2,IP), FanWei(1,2))
		FanWei(1,3) = min(XYZ(3,IP), FanWei(1,3))
270   CONTINUE
	FanWei(2,1) = FanWei(2,1) + (FanWei(2,1)-FanWei(1,1))*FDA/2.
	FanWei(2,2) = FanWei(2,2) + (FanWei(2,2)-FanWei(1,2))*FDA/2.
	FanWei(2,3) = FanWei(2,3) + (FanWei(2,3)-FanWei(1,3))*FDA/2.
	FanWei(1,1) = FanWei(1,1) - (FanWei(2,1)-FanWei(1,1))*FDA/2.
	FanWei(1,2) = FanWei(1,2) - (FanWei(2,2)-FanWei(1,2))*FDA/2.
	FanWei(1,3) = FanWei(1,3) - (FanWei(2,3)-FanWei(1,3))*FDA/2.


c	CALL OutputView(XYZ,ELEM,WElem,WNEB,
c     &			    FaceNum,NNPE,NNPS,CNTIVITY,ElemKind,
c     &			NRD0,NRD,NE,ME,NW,NWT,NP,MP,NDOF,N2D,N3D,
c     &			NKIND,NFACE)








C==============================================================================C	
C
C	Radiation Property
C
C
	WRITE(IUT0,9900) '<<--',ISTAGE,'-->> GET PROPERTY'
	ISTAGE=ISTAGE+1


	!this module is consistent as list_output_geom()
	MAT_NO=0
      DO IE=1,NE
		IMAT=MatID(IE)
		MAT_NO(IMAT)=1
	END DO
!
      NMAT=0
      DO IMAT=-100,100
		if(MAT_NO(IMAT).eq.1) then
		  NMAT=NMAT+1
		  MAT_NO(IMAT)=NMAT
		endif

	END DO
	NSOLID = NMAT-NFLUID
	write(IUT0,'(a,I3)')  '   MATERIAL = ', NMAT
	write(IUT0,'(a,I3)')  '      FLUID = ', NFLUID
	write(IUT0,'(a,I3)')  '      SOLID = ', NSOLID


        print*,'UUUUUUUU',NRD0
	!property
	CALL GetProperty(IUT0,radmodel,NE,NW,RadwallProp,RadwallType,
     &	 MatID,NMAT,NWT,nbcnd+1,N2D,NRD0,NRD,WGrayDeg,scabeta,
     &	 GAbsorb,PAbsorb,PScatter,WElem,VIndex,Prop,WArea,Volume,
     &	 Albedo,CTR,AveVs,GHomFlag,ETYPE,radwxflag)
	
	



C==============================================================================C	
C
C	Matching countpart for periodic boundaries
C
C
	IF(MarkPeriod.EQ.1) THEN
		WRITE(IUT0,9900) '<<--',ISTAGE,
     &		'-->> SET SURFACE-PAIR FOR PERIODIC BOUNDARY'
		ISTAGE=ISTAGE+1
	
		CALL  GetPrdPair(IUT0,nbcnd+1,prdData,RadWallProp,
     &		MP,NP,ME,NE,NWT,NW,NDOF,N2D,N3D,NKIND,NFACE,
     &		XYZ,ELEM,FanWei,AveVol)

	END IF



     


C==============================================================================C	
C
C	Estimate the memory size for depositing Exchange-area
C
C
	WRITE(IUT0,9900) '<<--',ISTAGE,'-->> GET ARRAY SIZE AND STACK-MODE'
	ISTAGE=ISTAGE+1
	
	!determine the coef. matrix memory size MRD0
	!and the Lower limits of the matrix coef.
	confid = 0.95
	VMIN=0.1
	
	CALL	GetArraySize(IUT0,NNN,MRD0,NRD0,NRD,MID,NE,NW,VMIN,
     &	 AveVs(1),AveVs(2),AveVol,confid,radmodel,NRAY,StackOpt)


	IF(RadGPOpt.EQ.0.AND.NRD.GT.1.5*RadMG) THEN
		
		WRITE(IUT0,*) '  ############################################'
		WRITE(IUT0,*) '   The mesh is too fine. it is preferred'
		WRITE(IUT0,*) '   to use double-scale scheme by setting'
		WRITE(IUT0,*) '   radgpopt = 1'
		WRITE(IUT0,*) '  ############################################'
		IF(NRD.GT.2**16-1) RadGPOpt=1
	END IF




C==============================================================================C	
C
C	Make coarse grids by Grouping the cells and surfaces
C	for large-scale mesh.  (>=5000)
C				-------- 2006/03/30
C		
C	
	IFlagGroup=0
	IF(RadGPOpt.EQ.1.AND.NRD.GT.1.5*RadMG) THEN 

	
		!-----Stack Gas Elem and Wall Elem for Quick Ray-Positioning
	WRITE(IUT0,9900) '<<--',ISTAGE,'-->> GROUP THE ELEMENTS'
	ISTAGE=ISTAGE+1
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!		
		!	<<	VIndex = IndexToGP    >>
		!		
		CALL MakeGroup(IUT0,RadMG,NRD,NRD0,NGP,nbcnd+1,
     &		MP,NP,ME,NE,NWT,NW,NDOF,N2D,N3D,NKIND,NFACE,NFLUID,
     &		XYZ,ELEM,ElemKind,FaceNum,NNPE,NNPS,CNTIVITY,
     &		MatID,FanWei,GlobeCTR,rwgpflag,IFlagGroup)

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF(IFlagGroup.EQ.1)  THEN
			NRD=NGP
			StackOpt = 'LTRI'
		END IF
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
			
	END IF



C==============================================================================C	
C
C	Allocate arrays for Exchange-Area,  << RDValue,RDIndex,RDN1,RDN2,RDId >>
C		
C		
	IF(StackOpt(1:4).EQ.'LTRI')	THEN
		
		ALLOCATE (RDN1(1),RDN2(1),RDId(1))
		ALLOCATE (RDIndex(NRD+1))
	
		RDIndex = 0
		DO 600 I = 2, NRD+1
			RDIndex(I) = RDIndex(I-1) + I-1
600		CONTINUE
		MRD = RDIndex(NRD+1)+100
		MRD0 = MRD
	
		IF(IFlagGroup.EQ.1) MRD0=NRD*NRD
	ELSE

		ALLOCATE (RDN1(NRD),RDN2(NRD),RDId(MID))
		ALLOCATE (RDIndex(NRD+1))

	END IF
	
	MRD=MRD0
	ALLOCATE (RDValue(MRD0),stat=ierr)
	IF(ierr.NE.0) STOP 'at Allocate RDValue in RadSXF()'

	IRD = 0
	IMRD = 0
	RDN1 = 0
	RDN2 = 0
	RDValue = 0.d0
	RDId = 0
	
	dseg = (0.75*AveVol/PAI)**0.33*NSEG
	
	WRITE(IUT0,*) '  NRAY=',NRAY
	WRITE(IUT0,*) '  NRD =', NRD
	WRITE(IUT0,*) '  MRD0=',MRD0
	WRITE(IUT0,*) '  Exchange-Area Save Mode=', TRIM(StackOpt)




C==============================================================================C	
C
C	Read Exchange-Area data from rdtmp.bin
C
C

	IF(TeaRst.EQ.0) THEN
	!-----Stack Gas Elem and Wall Elem for Quick Ray-Positioning
	WRITE(IUT0,9900) '<<--',ISTAGE,'-->> READ EXCHANGE-AREA DATA'
	ISTAGE=ISTAGE+1
	CALL ReadRestart(IUT2,NRD0,MRD0,NRD,MRD,StackOpt)
	GOTO 8888
	END IF

	
C***********************************************************************************C
C----
C----						2. CALCULATE EXCHANGE AREA (by MC / ZONE)
C----								
C***********************************************************************************C


C==============================================================================C	
C
C	Calculate the exchange-area for zone and monte-carlo method
C
C



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!			MONTE-CARLO METHOD
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF(RadModel(1:2).EQ.'MC')	THEN


!-----Stack Gas Elem and Wall Elem for Quick Ray-Positioning
		WRITE(IUT0,9900) '<<--',ISTAGE,'-->> STACK BUCKETS'
		ISTAGE=ISTAGE+1

		!optimized Block Array size
	TotVol=(FanWei(2,1)-FanWei(1,1))*(FanWei(2,1)-FanWei(1,1))
     &		   *(FanWei(2,1)-FanWei(1,1))
		NB= (TotVol/AveVol)**0.33*1.5


		!Gas Element
	CALL StackVolume(IUT0,XYZ,ELEM,MP,ME,N3D,NDOF,NP,NE,NNPE,
     &		NKIND,NMAT,ElemKind,MatID,FanWei,NPBG,NB,0.15d0)

	
		!Wall Element
		CALL StackWall(IUT0,XYZ,WElem,MP,NP,NWT,N2D,N3D,NDOF,
     &			   FanWei,WNV,NPBW,NB,0.8d0,wsmk,nbcnd+1)


			
!-----CALCULATE EXCHANGE-AREA VALUES BY M-C METHOD
	
		WRITE(IUT0,9900) '<<--',ISTAGE,'-->> ',
     &		'CALCULATE EXCHANGE-AREA BY M-C METHOD'
		ISTAGE=ISTAGE+1

		CALL CPU_TIME(TIME_TRACE1)
			
		
		mark3 = 0		!always=0
		NTRay=0			!total traced rays

		IF(GHomFlag.EQ.1) THEN
	  !Homogeneous media
	  WRITE(IUT0,*) '  Homogeneous Media '
		  
	  WRITE(IUT0,*) '  ---- Ray Emitted from Gas Elements... '
	  CALL GasRayTracer(IUT0,VMIN,NRAY,NP,NE,NW,MP,ME,NWT,
     2	NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD0,MID,IRD,
     3		IMRD,AveWArea,AveVol,NTRay,AveVs(3),AveVs(4),
     4		AveVs(5),NGP,StackOpt,IFlagGroup,mark3,GlobeCTR,FanWei,
     5		XYZ,ELEM)

		  WRITE(IUT0,*) '  ---- Ray Emitted from Wall Elements... '
	  CALL WallRayTracer(IUT0,VMIN,NRAY,NP,NE,NW,MP,ME,NWT,
     2	NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD0,MID,IRD,
     3	IMRD,AveWArea,AveVol,NTRay,AveVs(3),AveVs(4),
     4	AveVs(5),NGP,StackOpt,IFlagGroup,mark3,GlobeCTR,FanWei,
     5	XYZ,ELEM)
			  
		ELSE
		  !Inhomogeneous media
		  WRITE(IUT0,*) '  Inhomogeneous Media '
		  
		  
	  WRITE(IUT0,*) '  ---- Ray Emitted from Gas Elements... '
		  CALL GasRayTracerInhom(IUT0,VMIN,NRAY,NP,NE,NW,MP,
     2		ME,NWT,NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD0,
     3		MID,IRD,IMRD,dseg,AveWArea,AveVol,NTRay,
     4	AveVs(3),AveVs(4),AveVs(5),NGP,StackOpt,IFlagGroup,mark3,
     5        GlobeCTR,FanWei,XYZ,ELEM)

	  WRITE(IUT0,*) '  ---- Ray Emitted from Wall Elements... '
		  CALL WallRayTracerInhom(IUT0,VMIN,NRAY,NP,NE,NW,MP,
     2		ME,NWT,NB,NPBG,NPBW,NKIND,NFACE,N2D,N3D,NDOF,NRD0,MRD0,
     3		MID,IRD,IMRD,dseg,AveWArea,AveVol,NTRay,
     4	AveVs(3),AveVs(4),AveVs(5),NGP,StackOpt,IFlagGroup,mark3,
     5        GlobeCTR,FanWei,XYZ,ELEM)


		  IF(StackOpt(1:4).EQ.'BAND') THEN
			  MRD = IMRD
	  		  RDIndex(NRD+1) = MRD
			  WRITE(IUT0,*) '  Total Record Values = ', MRD
		  END IF
		  		  
		END IF

		DEALLOCATE(GBlockNum,GBlockIndex,WBlockNum,WBlockIndex)
		DEALLOCATE(GBlockID,WBlockID)
		
		CALL CPU_TIME(TIME_TRACE2)

	END IF
	
		

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!			ZONE METHOD
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF(radmodel(1:4).EQ.'ZONE')	THEN
		WRITE(IUT0,9900) '<<--',ISTAGE,'-->> ',
     &			'CALCULATE EXCHANGE-AREA BY ZONE METHOD'
		ISTAGE=ISTAGE+1



		ALLOCATE(VWK1(NRD0))		!VWK1 is Radius of the element

		!get equivalent radius
		DO IE =1,NE
			VWK1(IE) = (Volume(IE)/PAI*0.75)**0.3333
		END DO
		DO IW =1,NW
			IIT = NE+IW
			VWK1(IIT) = SQRT(WArea(IW)/PAI)
		END DO
		

		!calculate the exchange area
		CALL CPU_TIME(TIME_TRACE1)

		XiaoJuLi=0.0
		WRITE(IUT0,*) '  ---- Gas Elements ... '
		CALL GasZone(IUT0,NP,NE,MP,ME,NW,NWT,
     2			NKIND,NFACE,N2D,N3D,NDOF,NRD,NRD0,MRD,
     3			XYZ,ELEM,VWK1,XiaoJuli,
     5			AveVs(3)+AveVs(4),IFlagGroup)
		
		WRITE(IUT0,*) '  ---- Wall Elements	... '
		CALL WallZone(IUT0,NP,NE,MP,ME,NW,NWT,
     2			NKIND,NFACE,N2D,N3D,NDOF,NRD,NRD0,MRD,
     3			XYZ,ELEM,VWK1,XiaoJuli,
     5			AveVs(3)+AveVs(4),IFlagGroup)

		!self direct-exchange-area
		!CALL SelfZone(NRD,MRD,RDValue,RDIndex,PROP)

		CALL CPU_TIME(TIME_TRACE2)
	
		DEALLOCATE(VWK1)
	END IF


	t2 = TIME_TRACE2-TIME_TRACE1
	WRITE(IUT0,9910) 
     & ' TIME TOOK IN CALCULATING EXCHANGE-AREA =',t2,
     &	'[s]'


	

C==============================================================================C	
C
C	Smooth exchange-area by use of reprocity and consevation principles
C
C
	WRITE(IUT0,9900) '<<--',ISTAGE,'-->> SMOOTH THE EXCHANGE-AREA'
	ISTAGE=ISTAGE+1

	IF(IFlagGroup.EQ.1) THEN
		!RDValue is reallocated in this subroutine
		CALL SymmGroupRDV(IUT0,IUT3,NRD0,MRD0,NRD,MRD)
	END IF


	ALLOCATE (VWK1(NRD),VWK2(NRD))


	IF(IFlagGroup.EQ.1)  THEN
		VWK1(:)=GpProp(:)
	ELSE
		DO I=1,NRD0
			IRD=VIndex(I)
			IF(IRD.GT.0) VWK1(IRD)=PROP(I)
		END DO
	END IF

	
	IF(StackOpt(1:4).EQ.'LTRI') THEN 

	CALL SmoothLTRI_2(IUT0,NRD,MRD,RDIndex,RDValue,VWK1,VWK2,
     &					feps,vfeps,RLX,itr)
	
	ELSE

		CALL SmoothBAND(IUT0,NE,NW,NRD,MRD0,MRD,MID,
     &				RDId,RDN1,RDN2,RDIndex,RDValue,
     &				VWK1,VMIN,feps,vfeps,RLX,iterate)

	END IF

	DEALLOCATE (VWK1,VWK2)


	!Get the Exchange-Area between each wall groups
	WRITE(IUT0,9900) '<<--',ISTAGE,
     &			'-->> CAL. EXCHANGE-AREA BETWEEN WALLS'
	ISTAGE=ISTAGE+1
	CALL WallGroupEASum(IUT0,IUT0+5,NE,NWT,NRD0,NRD,MRD,MID,N2D,
     &		nbcnd+1,StackOpt,RDIndex,RDValue,RDId,PROP,EASum,
     &		WElem,VIndex)






C==============================================================================C	
C
C	Energy Equation
C
C

2222	WRITE(IUT0,9900) '<<--',ISTAGE,
     &			'-->> MAKE ENERGY EQUATION FOR <FrontFlow/RED>'
	ISTAGE=ISTAGE+1

	!Make C * T = B
	IF(IFlagGroup.EQ.1) THEN
          CALL MakeEnergyEquGP(IUT0,NRD0,NRD,MRD)
	ELSE
          CALL MakeEnergyEqu(VMIN,StackOpt,NRD0,NRD,MRD,MID)
	END IF





	
C==============================================================================C	
C
C	output temporary data of exchange-area
C
C
	IF(TeaRst.EQ.1) CALL WriteRestart(IUT2,NRD0,MRD0,NRD,MRD,StackOpt)
	






C==============================================================================C	
C
C	Make the index from vertex to cell / cell-group
C
C
8888	WRITE(IUT0,9900) '<<--',ISTAGE,
     &		'-->> DATA TRANSFER INDEX BETWEEN NODE AND CELL'
	ISTAGE=ISTAGE+1

	CALL RDVElemToNode(IUT0,MM,ME,MP,NP,NE,NW,NWT,N2D,N3D,NKind,
     &	NDOF,NRD,NRD0,nbcnd+1,NMAT,MAT_NO,ELEM,MatID,XYZ,IFlagGroup,
     &	ICPU)
					!NOTE: nbcnd+1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	NRD0 = NP
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



c	IF(IFlagGroup.EQ.1) THEN
c	  
c	  !Make the index from vertex to cell-group
c	  CALL NodeToGPProp(IUT0,MM,ME,MP,NP,NE,NW,NWT,N2D,N3D,
c     &			NKind,NRD,ELEM)
c
c	ELSE
c	
c	  !Transfer the Exchange-Area data from cells to vertices
c	  CALL RDVElemToNode(IUT0,MM,ME,MP,NP,NE,NW,NWT,NRD0,
c     &			NRD,MRD,MRD0,MID,N2D,N3D,NKind,StackOpt,NPIN,
c     &			ELEM)
c	END IF
c	!Get neighboring nodes of boundary vertex
c	CALL GetNodeNeb(IUT0,WElem,ELEM,
c     &				ElemKind,NNPE,XYZ,NPIN,nbcnd,
c     &			    ME,MP,NP,NE,NW,NWT,N2D,N3D,NDOF,NKind)



	
C***********************************************************************************C
C----
C----								3. FINISH
C----
C***********************************************************************************C

	DEALLOCATE (WElem,WNV,WABCD,WArea,NodeBounID,WNEB)
	DEALLOCATE (Volume,GAbsorb,WGrayDeg,PAbsorb,PScatter,Albedo)
	DEALLOCATE (scabeta,CTR,ETYPE)
	DEALLOCATE(FaceNum,NNPE,NNPS,CNTIVITY,ElemKind)
	IF(IFlagGroup.EQ.1) DEALLOCATE (GpCNum,GpProp)

		

	CALL CPU_TIME(TIME_END)

	t1 = TIME_END-TIME_BEGAN

	WRITE(IUT0,9900) '<<--',ISTAGE,	'-->> FINISHED'
	ISTAGE=ISTAGE+1
	
	WRITE(IUT0,9910) ' TOTAL CPU TIME = ', t1, ' (s)'
	IF(TeaRst.EQ.1) THEN
	IF(RadModel(1:3).NE.'FVM')
     & WRITE(IUT0,9910) 
     &  ' TIME IN CALCULATING EXCHANGE-AREA = ',t2,'(s)'
	IF(RadModel(1:2).EQ.'MC') 
     & WRITE(IUT0,'(a,I12)') '  Total Traced Rays =', NTRay
	END IF

	WRITE(IUT0,*) '=================================================='
	WRITE(IUT0,*) '-                                                -'
	WRITE(IUT0,*) '-         RADSXF FINISHED SUCCESSFULLY           -'
	WRITE(IUT0,*) '-                                                -'
	WRITE(IUT0,*) '=================================================='	
	



	RETURN
	
9900	FORMAT(1x,a,I2,a)
9910	FORMAT(1x,a,F9.2,a)

!--------- End of subroutine radsxf()

	
	END SUBROUTINE RADSXF



