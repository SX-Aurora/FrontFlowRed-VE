!
!      module module_include
!      module module_hpc_input
!      module module_hpcutil
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      module module_include
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! --- array dimension
!
!      include 'dimen.h_hpc'
!
!      end module module_include
!
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_hpc_input
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      character(20),parameter,private :: modnam='(module_hpc_input)'
      integer             :: NPART,UCDFLAG,WALLFLAG
      character (len=80)  :: GRIDFIL,BCFIL,INPFIL,UCDFIL
      character (len=80)  :: WALFIL,DIRHED
      character (len=80)  :: HEADERg, HEADERb, HEADERc, HEADERw
      character (len=80)  :: HEADERs,HEADERi,HEADERr,HEADERt,HEADERf
      character (len=80)  :: HEADER_Pr,HEADER_Pi
      character (len=80)  :: INITin_P,RESTRTout_P
      
      character (len=80)  :: HEADERa,HEADERm
      character (len=80)  :: GRIDin,BCin,COMMin
      character (len=80)  :: SRCin,INITin,RESLTout,RESTRTout,RESLTgrid
      character (len=80)  :: FORCEout,GEOMin,WALLin
      character (len=80)  :: statis,animfile,sld_wk,suf_wk
      character (len=80) :: WLDISout
      character (len=80) :: dflxTmpFname1,dflxTmpFname2
      character (len=80)  :: HEADER_pro,PROBEout     !!!modify
!
!/////////////////////////////////////////////////////////////////////
      contains
!=================================================
      subroutine inputdata(my_rank,NPETOT,ifli,ifll,ifle,cntlnam,
     &                  Dmns,ierror)
!=================================================
      character(11),parameter :: subnam='(inputdata)'
!
! --- [dummy qrgument]
!
      integer,intent(in)     :: ifli,ifll,ifle,my_rank,NPETOT
      character,intent(in)   :: cntlnam
      integer,intent(out)    :: ierror
      character(len=80),intent(out) :: Dmns
!
! --- [local entities]
!
      integer :: ierr1=0
      integer           :: hpcfil=11
      integer           :: ID,LENGTH,LENGT
      character(len= 4) :: SUBindex
      character(80)     :: mtsfil,boundary,vertex,wall,subdir
      character(80)     :: hpc_boundary,hpc_comm,hpc_vertex,ucdfile
      character(80)     :: hpc_anim
      character(80)     :: hpc_source,hpc_result
      character(80)     :: hpc_initial,hpc_restart,hpc_force
      character(80),parameter :: blank=' '
      integer           :: NPE,walflg,ucdflg,ios
      character(80)     :: hpc_p_result
      character(80)     :: hpc_p_restart
      character(80)     :: hpc_p_initial
      character(80)     :: hpc_probe       !!!modify
!
! --- [namelist]
!
      namelist /hpc_cntl/
     &        NPE,mtsfil,boundary,vertex,wall,subdir,
     &        hpc_vertex,hpc_boundary,hpc_comm,hpc_force,
     &        hpc_source,hpc_result,hpc_initial,hpc_restart,
     &        hpc_anim,
     &        walflg,ucdflg,ucdfile,
!     &        hpc_p_result,
     &        hpc_p_restart,
     &        hpc_p_initial
     &        ,hpc_probe         !!!modify
!
! --- Start:
!
      ierror=0
!
! --- namelist initialization
!
      NPE=1
      walflg=0
      ucdflg=0
      hpc_boundary=blank
      hpc_comm=blank
      hpc_vertex=blank
      hpc_source=blank
      hpc_result=blank
      hpc_initial=blank
      hpc_restart=blank
      hpc_force=blank
      hpc_anim=blank
      mtsfil=blank
      boundary=blank
      vertex=blank
      wall=blank
      subdir=blank
      ucdfil=blank
      hpc_probe=blank
!
! --- variable initialization
!
      NPART=1
      UCDFLAG=0
      WALLFLAG=0
      GRIDFIL=blank
      BCFIL=blank
      INPFIL=blank
      UCDFIL=blank
      WALFIL=blank
      DIRHED=blank
      HEADERg=blank
      HEADERb=blank
      HEADERc=blank
      HEADERw=blank
      HEADERs=blank
      HEADERi=blank
      HEADERr=blank

      HEADER_Pr=blank
      HEADER_Pi=blank

      HEADERt=blank
      HEADERm=blank
      HEADERf=blank
      HEADERa=blank
      SRCin=blank
      INITin=blank
      RESLTout=blank
      RESLTgrid=blank
      RESTRTout=blank

      INITin_P=blank
      RESTRTout_P=blank

      animfile=blank
      FORCEout=blank
      GEOMin=blank
      WALLin=blank
      statis=blank
      sld_wk=blank
      suf_wk=blank

      HEADER_pro=blank  !!!modify
      PROBEout=blank    !!!modify
      
!++++Modified by T. Unemura 030910++++++++++++++++++++++++++++++++++++++
      WLDISout=blank
      dflxTmpFname1=blank
      dflxTmpFname2=blank
!
! --- Read namelist
!
      hpcfil=11!NPETOT  !+my_rank+1
!      if(hpcfil>99.or.hpcfil<1) then
!        write(ifle,*) '*** Warning: over-flow FORTRAN max unit number',
!     &            'in opening hpc.cntl file.'
!        goto 9999
!      end if
!
!      open (hpcfil,file='CONT.HPC',status='unknown')
!
      rewind ifli
      read(ifli,hpc_cntl,iostat=ios)
!      close(hpcfil)
!      write(ifll,hpc_cntl)
      call nml_errmsg0(ifle,ios,'hpc_cntl',ierr1)
      if( ierr1.ne.0 ) goto 9999
!
      IF(NPE.NE.NPETOT) THEN
        write(ifle,*) 
     &   '*** ERR: NPE No. in CONT.HPC is NOT equal to ncpu in cort.1'
        WRITE(ifle,*) '### error : NPE must be equal to ',NPETOT
        call FFRABORT(1,'module_hpc_input/inputdata')
      ELSEIF(NPE.LE.1) then
        WRITE(ifle,*) '### error : NPE must be greater than 1' 
        write(ifle,*) '*** ATT: Are sure to use HPC'
        call FFRABORT(1,'module_hpc_input/inputdata')
      ELSE
        NPART=NPE
      ENDIF
!
      UCDFLAG=ucdflg
      UCDFIL=ucdfile
      IF(UCDFLAG.eq.1) THEN
        if(UCDFIL.eq.blank) then
        WRITE(ifle,*) '### error : ucdfile not de defined if ucdfile=1'
          goto 9999
        endif
      ELSEIF(UCDFLAG.ne.1.and.UCDFLAG.ne.0) THEN
        WRITE(ifle,*) '### error : ucdflg must be 1 or 0'
        goto 9999
      ENDIF
!
      WALLFLAG=walflg
      WALFIL=wall
      IF(walflg.EQ.1) THEN
        IF(WALFIL.EQ.blank) THEN
        call errmsg(GRIDFIL,'vertex')
          WRITE(ifle,*) '### error :  wall not de defined if walflg=1'
          goto 9999
        ENDIF
      ELSEIF(WALLFLAG.ne.1.and.WALLFLAG.ne.0) THEN
        WRITE(ifle,*) '### error :  walflg must be 1 or 0'
        goto 9999
      ENDIF
!
      GRIDFIL=vertex
      BCFIL=boundary
      INPFIL=mtsfil
!
      call errmsg(INPFIL,'mtsfil')
!      call errmsg(GRIDFIL,'vertex')
!      call errmsg(BCFIL,'boundary')
!
      DIRHED=subdir
      HEADERg=hpc_vertex
      HEADERb=hpc_boundary
      HEADERc=hpc_comm
      HEADERs=hpc_source
      HEADERi=hpc_initial
      HEADERr=hpc_restart

      HEADER_Pr=hpc_p_restart
      HEADER_Pi =hpc_p_initial

      HEADERt=hpc_result
      HEADERm='move_grid'
      HEADERf=hpc_force
      HEADERw='geom'
      HEADERa=hpc_anim
      
      HEADER_pro=hpc_probe     !!!modify

      call errmsg(DIRHED,'subdir')
      call errmsg(HEADERg,'hpc_vertex')
      call errmsg(HEADERb,'hpc_boundary')
      call errmsg(HEADERc,'hpc_comm')
      call errmsg(HEADERa,'hpc_anim')
!      call errmsg(HEADERs,'hpc_source')
!      call errmsg(HEADERi,'hpc_initial')
      call errmsg(HEADERr,'hpc_restart')
      call errmsg(HEADERt,'hpc_result')
      call errmsg(HEADER_pro,'hpc_probe')     !!!modify
!
! --- 
!
      write(SUBindex,'(i4.4)') my_rank
      DIRHED=adjustL(DIRHED)
      LENGT=len_trim(DIRHED)
!
      HEADERg=adjustL(HEADERg)
      LENGTH=len_trim(HEADERg)
      GRIDin=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERg(1:LENGTH)
!
      HEADERb=adjustL(HEADERb)
      LENGTH=len_trim(HEADERb)
      BCin=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERb(1:LENGTH)
!
      HEADERc=adjustL(HEADERc)
      LENGTH=len_trim(HEADERc)
      COMMin=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERc(1:LENGTH)
!
      Dmns=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'dimen.parm'
!
      HEADERs=adjustL(HEADERs)
      LENGTH=len_trim(HEADERs)      
      SRCin=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERs(1:LENGTH)
!
      HEADERi=adjustL(HEADERi)
      LENGTH=len_trim(HEADERi)      
      INITin=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERi(1:LENGTH)
!
      HEADERr=adjustL(HEADERr)
      LENGTH=len_trim(HEADERr)      
      RESTRTout=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERr(1:LENGTH)
!
      HEADER_Pr=adjustL(HEADER_Pr)
      LENGTH=len_trim(HEADER_Pr)      
      RESTRTout_P=
     &     DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADER_Pr(1:LENGTH)
!
      HEADER_Pi=adjustL(HEADER_Pi)
      LENGTH=len_trim(HEADER_Pi)      
      INITin_P=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADER_Pi(1:LENGTH)
!

      HEADERt=adjustL(HEADERt)
      LENGTH=len_trim(HEADERt)      
      RESLTout=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERt(1:LENGTH)
!
      HEADERm=adjustL(HEADERm)
      LENGTH=len_trim(HEADERm)      
      RESLTgrid=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERm(1:LENGTH)
!
      HEADERf=adjustL(HEADERf)
      LENGTH=len_trim(HEADERf)   
      FORCEout=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERf(1:LENGTH)
!
      HEADERa=adjustL(HEADERa)
      LENGTH=len_trim(HEADERa)   
      animfile=DIRHED(1:LENGT)//'_'//SUBindex//'/'//HEADERa(1:LENGTH)
!
      GEOMin=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'geom_hpc'
!
      WALLin=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'wall_hpc'
!
      statis=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'statis'
!
      sld_wk=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'sldng_wk'
!
      suf_wk=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'eqrate_wk'

      HEADER_pro=adjustL(HEADER_pro)       !!!modify
      LENGTH=len_trim(HEADER_pro)          !!!modify
      PROBEout=DIRHED(1:LENGT)//'_'//SUBindex//'/'//
     &         HEADER_pro(1:LENGTH)//'_'//SUBindex//'.txt' !!!modify
!
!++++Modified by T. Unemura 030910++++++++++++++++++++++++++++++++++++++
      WLDISout=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'wall_dist'
!++++Modified by T. Unemura 030910++++++++++++++++++++++++++++++++++++++
      dflxTmpFname1=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'dflx_tmp1.dat'
      dflxTmpFname2=DIRHED(1:LENGT)//'_'//SUBindex//'/'//'dflx_tmp2.dat'
      return
!
 9999 continue
      ierror=1
      write(ifle,*) modnam,subnam
      end subroutine inputdata
!
!=================================================
      subroutine errmsg(vnam1,vnam2)
!=================================================
      character(*) :: vnam1,vnam2
      if( vnam1.ne.' ' ) return
      write(*,*) '### error : data error'
      write(*,*) 'lack of data : ',vnam2
      call FFRABORT(1,'module_hpc_input/errmsg')
      end subroutine errmsg
!
      end module module_hpc_input
!

!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      module module_hpcutil
!SS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! --- This module is real hpc module
!
        use module_hpc_input,
     &                 only : NPART,UCDFLAG,WALLFLAG,
     &                        GRIDFIL,BCFIL,INPFIL,UCDFIL,
     &                        WALFIL,DIRHED,
     &                        HEADERg, HEADERb, HEADERc, HEADERw,
     &                        HEADERs,HEADERi,HEADERr,HEADERt,HEADERf,
     &                        HEADERa,HEADERm,
     &                        HEADER_Pr,HEADER_Pi,
     &                        GRIDin,BCin,COMMin,
     &                        SRCin,INITin,RESLTout,RESTRTout,
     &                        INITin_P,RESTRTout_P,
     &                        FORCEout,GEOMin,WALLin,
     &                        statis,animfile,sld_wk,suf_wk
     &                        ,WLDISout,RESLTgrid
     &                        ,dflxTmpFname1,dflxTmpFname2
     &                        ,HEADER_pro,PROBEout   !!!modify
!
      include 'mpif.h'
!
      integer              :: PETOT, my_rank,SOLVER_COMM
      integer,public       :: errno
      integer,public       :: NPE,ROOT
!
      integer              :: NEIBPETOT, NEIBPETOTw
      integer              :: INTNODTOT, NODTOT, NODTOTw
      integer              :: NODtotDIFF
      integer              :: ELMtot
      integer              :: LENT, LENw
!
!
      integer ,allocatable :: NEIBPE(:)
      integer ,allocatable :: IMPORT_index(:), IMPORT_item(:),
     &                        EXPORT_index(:), EXPORT_item(:)
!
      integer ,allocatable :: IMPORT_indexIP(:),
     &                        EXPORT_indexIP(:)
!      integer, allocatable,save :: IPSND(:)
!
      integer ,allocatable :: IBC_P_NO(:,:,:,:)
      integer ,allocatable :: IBC_P_TOL(:,:,:)
!
      integer ,allocatable :: NEIBPEw(:)
      integer ,allocatable :: IMPORTw_index(:), IMPORTw_item(:),
     &                        EXPORTw_index(:), EXPORTw_item(:)
      
!
      
      integer ,allocatable :: NOD_IDg(:), ELM_IDg(:)       
      character(len=80)    :: LINE
!
!      integer,allocatable :: CPULST(:)
!
      end module module_hpcutil
!
