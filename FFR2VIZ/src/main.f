!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine main(argc,arrayArg,ierror)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      
      use FFRdata, only : nvrtx,cord,numVars,numBVars,wallDist,
     &                    nameVars,nameBVars,totdata,bnddata,
     &                    MPprocNum,outGGrpNum,cord_bak
      use FFRreaddata, only : icomp_avex,icomp_rmsx
      use ProgData, only : ypAverage,ypLogDiv,
     &                     divnumYPAV,wdisCalcType,
     &                     axis1,axis2,Radius,scaleFactor,
     &                     axisTag,wdisYPAVfn,
     &                     ypAvupfn,ypVal,gridDataMP,
     &                     ReTau,uTau,wdMax,wdMin,dyp

      use tovtk
!
      implicit none
!      
! --- 
!      
      integer,intent(in) :: argc
      character*56,intent(in) :: arrayArg(1:argc)
      integer,intent(out) :: ierror
!     
      character*80  :: hpcprefix
      character*255 :: hpcPath    ! extern function

      character*32 :: outFVprfx,outExt(2)
      character*80 :: gridData,Rsltfn,outFVfn
      character*84 :: gridDataSC
      character*255 :: filePath
      character*2 :: resultFormat
      character*4 :: gridFormat
!
      character*4 :: ffgFormatKey
      logical :: MPflag,procSelect
      logical,allocatable :: targProc(:)
! DEBUG
      logical :: gridOnly,procIDVals,vrtxIDVals
!      logical :: gridOnly,procIDVals
      real*8  :: gdScale
! for boundary y+
      logical :: boundaryYp
      character*80 :: bndYpWallfn
! for animation
      logical :: animeStat
      integer :: animeStep
      character*32 :: animeExt
      logical :: rsOpenStat,rsCloseStat
! for wall set t o0
      logical :: wallZero
! for relative velocity
      logical :: RelativeVelocity
      
      integer :: i,j,k,m,n,rtrnStat,ios,procSlctS,procSlctE
      integer :: len,l1,l2
      logical :: fileEx,inqChk,showUsage
      character*1 :: ans,arg1
      character(len=10) :: featureNameGF = ""
      character(len=10) :: featureNameRF = ""
      character(len=0)  :: versionGF = ""
      character(len=0)  :: versionRF = ""
! ffr2viz or ffrmovie
!      animeStat=.true.  ! == MOVIE
      animeStat=.false. ! ==VIZ
!
! --- Title 
!
      write(*,*) '                                              '
      write(*,*) '_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'
      write(*,*) '_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'
      if(.not.animeStat) then
         write(*,*) '      * FFR2VIZ *                             '
      else
         write(*,*) '      * FFRMOVIE *                             '
      end if
      write(*,*) '                 VERSION    = src_073.2 '
      write(*,*) '                 BUILD DATE = 2008-09-19-05:19:04 '
!      write(*,*) '           - (c) AdvanceSoft Corporation -    '
      write(*,*) '_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'
      write(*,*) '_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'
      write(*,*) '                                              '
!
! --- Default Value for Index
!
      gridData=' '
      outFVprfx=' '
      outFVfn=' '
      outExt='.fv'
      gridFormat='FV'
      ffgFormatKey='U'
      hpcprefix='hpc'
      Rsltfn=' '
      procSelect=.false.
      MPflag=.false.
      MPprocNum=1;
      resultFormat='FV'
      gridOnly=.false.
      procIDVals=.false.
! DEBUG
      vrtxIDVals=.false.
      gdScale=1.d0
      ypAverage=.false.
      wdisYPAVfn=' '
      ypAvupfn=' '
      divnumYPAV=1
      ypLogDiv=.false.
      ReTau=1.d0
      uTau=1.d0
      wdisCalcType=-1
      showUsage=.false.
!
      boundaryYp=.false.
      bndYpWallfn='wall_dist'
      wallZero=.false.
      RelativeVelocity=.false.

!
! --- Check default velaue
!
      do 1020 i=1,argc
        if(arrayArg(i)=='-m') then
          MPflag=.true.
          arg1=arrayArg(i+1)
          if(.not.((arg1>='0').and.(arg1<='9'))) then
            write(*,*)'You must assign integer number after -m.'
            showUsage=.true.
          end if
          read(arrayArg(i+1),*)MPprocNum
          write(*,*)'Data for MP calculation of ',MPprocNum,' CPUs.'
!
        else if(arrayArg(i)=='-s') then
          MPflag=.false.;     MPprocNum=1
!
        else if(arrayArg(i)=='-a') then
          procSelect=.false.
!
        else if(arrayArg(i)=='-p') then
          procSelect=.true.
          procSlctS=i+1;      procSlctE=i
          do 1021 j=i+1,argc
            arg1=arrayArg(j)
            if(.not.((arg1>='0').and.(arg1<='9'))) exit
            procSlctE=j
 1021     continue
          if(procSlctS>procSlctE) then
            write(*,*)'You must assign integer number after -p.'
            showUsage=.true.
          end if
!
        else if(arrayArg(i)=='-g') then
          arg1=arrayArg(i+1)
          if(arg1=='-') then
            write(*,*)'You must assign file name after -g.'
            showUsage=.true.
          end if
          gridData=trim(adjustl(arrayArg(i+1)))
!
        else if(arrayArg(i)=='-gs') then
          read(arrayArg(i+1),*)gdScale
          write(*,*)'Grid Scaling : (X,Y,Z) * ',gdScale
!
        else if(arrayArg(i)=='-o') then
          arg1=arrayArg(i+1)
          if(arg1=='-') then
            write(*,*)'You must assign file name after -o.'
            showUsage=.true.
          end if
          outFVprfx=trim(adjustl(arrayArg(i+1)))
!
        else if(arrayArg(i)=='-d') then
          arg1=arrayArg(i+1)
          if(arg1=='-') then
            write(*,*)'You must assign directory prefix after -d.'
            showUsage=.true.
          end if
          hpcprefix=trim(adjustl(arrayArg(i+1)))
!
        else if(arrayArg(i)=='-r') then
          arg1=arrayArg(i+1)
          if(arg1=='-') then
            write(*,*)'You must assign file name after -r.'
            showUsage=.true.
          end if
          Rsltfn=trim(adjustl(arrayArg(i+1)))
!
        else if(arrayArg(i)=='-rf') then
          if((arrayArg(i+1)=='FV').or.(arrayArg(i+1)=='fv')) then
            resultFormat='FV'
            outExt(1)='.fv'
! license
            featureNameRF = ""
            versionRF     = ""
          else if((arrayArg(i+1)=='vtk')
     &                           .or.(arrayArg(i+1)=='VTK')) then
            write(*,*)'*** output format is VTK'
            resultFormat='VT'
            outExt(1)='.vtk'
! license
            featureNameRF = "AFFR-VTK"//char(0)
            versionRF     = versionRF //char(0)
          else if((arrayArg(i+1)=='vtka')
     &                           .or.(arrayArg(i+1)=='VTKA')) then
            write(*,*)'*** output format is VTK (ascii)'
            resultFormat='VA'
            outExt(1)='.vtk'
! license
            featureNameRF = "AFFR-VTK"//char(0)
            versionRF     = versionRF //char(0)
!#ifdef FVASCII
          else if((arrayArg(i+1)=='FVA').or.(arrayArg(i+1)=='fva')) then
            resultFormat='FA'
            outExt(1)='.fv'
! license
            featureNameRF = ""
            versionRF     = ""
!#endif
          else if((arrayArg(i+1)=='FVS').or.(arrayArg(i+1)=='fvs')) then
            resultFormat='FS'
            outExt(1)='.gfv'
            outExt(2)='.rfv'
! license
            featureNameRF = ""
            versionRF     = ""
!          else if((arrayArg(i+1)=='GF').or.(arrayArg(i+1)=='gf')) then
!            resultFormat='GF'
!            outExt(1)='.gf'
!#ifdef AVS
          else if((arrayArg(i+1)=='AVS').or.(arrayArg(i+1)=='avs')) then
            resultFormat='AA'
            outExt(1)='.inp'
! license
            featureNameRF = "AFFR-AVS"//char(0)
            versionRF     = versionRF //char(0)
!#endif
!#ifdef FLUENT
          else if((arrayArg(i+1)=='FL') .or.(arrayArg(i+1)=='fl')) then
            resultFormat='FL'
            outExt(1)='.cas'
            outExt(2)='.dat'
! license
            featureNameRF = "AFFR-GB"//char(0)
            versionRF     = versionRF//char(0)
!#endif
!#ifdef SCRYU
          else if((arrayArg(i+1)=='SY')
     &                           .or.(arrayArg(i+1)=='sy')) then
            resultFormat='SY'
            outExt(1)='.fld'
! license
            featureNameRF = "AFFR-SY"//char(0)
            versionRF     = versionRF//char(0)
!#endif
!#ifdef STARCD
          else if((arrayArg(i+1)=='SC')
     &                           .or.(arrayArg(i+1)=='sc')) then
            resultFormat='SC'
            outExt(1)='.usr'
! license
            featureNameRF = "AFFR-SC"//char(0)
            versionRF     = versionRF//char(0)
!#endif
!
!#ifdef ENSIGHT
          else if((arrayArg(i+1)=='ES') .or.(arrayArg(i+1)=='es')) then
            resultFormat='ES'
            outExt(1)='.geo'
            outExt(2)='.case'
! license
            featureNameRF = "AFFR-ES"//char(0)
            versionRF     = versionRF//char(0)
!#endif
          else
            write(*,*)'You must assign file format key after -rf.'
            showUsage=.true.
          end if
!
        else if(arrayArg(i)=='-gf') then
          if((arrayArg(i+1)=='FV').or.(arrayArg(i+1)=='fv')) then
            gridFormat='FV'
!#ifdef STARCD
          else if((arrayArg(i+1)=='SC')
     &                       .or.(arrayArg(i+1)=='sc')) then
            gridFormat='SC'
! license
            featureNameGF = "AFFR-SC"//char(0)
            versionGF     = versionGF//char(0)
!#endif
!#ifdef GAMBIT
          else if((arrayArg(i+1)=='GB')
     &                       .or.(arrayArg(i+1)=='gb')) then
            gridFormat='GB'
! license
            featureNameGF = "AFFR-GB"//char(0)
            versionGF     = versionGF//char(0)
!          else if((arrayArg(i+1)=='GB-T')
!     &                       .or.(arrayArg(i+1)=='gb-t')) then
!            gridFormat='GB-T'
!          else if((arrayArg(i+1)=='GB-H')
!     &                       .or.(arrayArg(i+1)=='gb-h')) then
!            gridFormat='GB-H'
!#endif
!#ifdef SCRYU
          else if((arrayArg(i+1)=='SY')
     &                       .or.(arrayArg(i+1)=='sy')) then
            gridFormat='SY'
! license
            featureNameGF = "AFFR-SY"//char(0)
            versionGF     = versionGF//char(0)
!#endif
!
          else if((arrayArg(i+1)=='FF')
     &                       .or.(arrayArg(i+1)=='ff')) then
            gridFormat='FF'
! license
            featureNameGF = ""
            versionGF     = ""
          else if((arrayArg(i+1)=='GF')
     &                       .or.(arrayArg(i+1)=='gf')) then
            gridFormat='GF'
! license
            featureNameGF = ""
            versionGF     = ""
          else
            write(*,*)'You must assign grid format key after -gf.'
            showUsage=.true.
          end if
!
        else if(arrayArg(i)=='-nr') then
          gridOnly=.true.
!
        else if(arrayArg(i)=='-fm') then
          if((arrayArg(i+1)=='A').or.(arrayArg(i+1)=='a')) then
             ffgFormatKey='A'
          else
             ffgFormatKey='U'
          end if
!
        else if(arrayArg(i)=='-pv') then
          procIDVals=.true.
!
! DEBUG
        else if(arrayArg(i)=='-vn') then
          vrtxIDVals=.true.
!
! onishi under develop
c        else if(arrayArg(i)=='-byp') then
c          boundaryYp=.true.
c          arg1=arrayArg(i+1)
c          if(arg1=='-') then
c            write(*,*)'You must assign wall_dist file name after -byp.'
c            showUsage=.true.
c!            stop
c          end if
c          bndYpWallfn=trim(adjustl(arrayArg(i+1)))
!
        else if(arrayArg(i)=='-wz') then
          wallZero=.true.
!
        else if(arrayArg(i)=='-rv') then
          RelativeVelocity=.true.
!
        else if(arrayArg(i)=='-av') then
          read(arrayArg(i+1),*,iostat=ios)divnumYPAV
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
          arg1=arrayArg(i+2)
          if(arg1=='-') then
            write(*,*)'You must assign a file name after -av.'
            showUsage=.true.
          end if
          ypAvupfn=trim(adjustl(arrayArg(i+2)))
          ypAverage=.true.
!
        else if(arrayArg(i)=='-yp') then
          read(arrayArg(i+1),*,iostat=ios)ReTau
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
!
        else if(arrayArg(i)=='-up') then
          read(arrayArg(i+1),*,iostat=ios)uTau
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
!
        else if(arrayArg(i)=='-wdf') then
          arg1=arrayArg(i+1)
          if(arg1=='-') then
            write(*,*)'You must assign file name after -wdf.'
            showUsage=.true.
          end if
          wdisYPAVfn=trim(adjustl(arrayArg(i+1)))
          read(arrayArg(i+2),*,iostat=ios)scaleFactor
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
          wdisCalcType=0
!
        else if(arrayArg(i)=='-cyl') then
          arg1=arrayArg(i+1)
          if((arg1/='x').and.(arg1/='X').and.(arg1/='y').and.
     &            (arg1/='Y').and.(arg1/='z').and.(arg1/='Z')) then
            write(*,*)'You must assign x, y or z after -cyl.'
            showUsage=.true.
          end if
          axisTag=arg1
          read(arrayArg(i+2),*,iostat=ios)axis1
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
          read(arrayArg(i+3),*,iostat=ios)axis2
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
          read(arrayArg(i+4),*,iostat=ios)Radius
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
          read(arrayArg(i+5),*,iostat=ios)scaleFactor
          if(ios/=0) then
            write(*,*)'Option error.'
            showUsage=.true.
          end if
          wdisCalcType=2
!
        else if(arrayArg(i)=='-led') then
          ypLogDiv=.true.
!
        else if(arrayArg(i)=='-h') then
          showUsage=.true.
! DEBUG
        else
          if(arrayArg(i)(1:1)=='-') then
             write(*,*) '***'
             write(*,*) '*** Invalid Option:',arrayArg(i)
             write(*,*) '***'
             showUsage=.true.
          end if
        end if
 1020 continue
!
! --- HPC processor selection

!      if(gridOnly) then
!         if(MPflag .or. MPprocNum>1) then
!            write(*,*) ' *** Warning: -m option is ignored.'
!         end if
!         MPflag=.false.
!         MPprocNum=1
!      end if
!
      if((.not.(MPflag)).and.(procSelect)) procSelect=.false.
      allocate(targProc(0:MPprocNum-1));  targProc=.false.
      if(procSelect) then
        do 1200 i=procSlctS,procSlctE
          read(arrayArg(i),*,iostat=ios)m
          if(ios/=0) then
            write(*,*)'ABEND:CPU selection error.'
            showUsage=.true.
          end if
          if((m>=0).or.(m<MPprocNum)) then
            targProc(m)=.true.
          end if
 1200   continue
      else
        targProc(0:MPprocNum-1)=.true.
      end if
      outGGrpNum=0
      do 1210 i=0,MPprocNum-1
        if(targProc(i)) outGGrpNum=outGGrpNum+1
 1210 continue
      if(MPprocNum==1) then
        write(*,*)'Convert Single CPU data.'
      else
        write(*,*)'Merging ',outGGrpNum,' CPUs data.'
      end if
!
!
!
      if((gridData==' ').or.(outFVprfx==' ').or.(gridFormat==' ')) then
        write(*,*)'Lack of option arguments.'
        showUsage=.true.
      end if
      if(MPflag) then
        if((hpcprefix==' ')) then
          write(*,*)'Lack of option arguments.'
          showUsage=.true.
        end if
      else
        if((.not.(gridOnly)).and.(Rsltfn==' ')) then
          write(*,*)'Lack of option arguments.'
          showUsage=.true.
        end if
      end if
!
! --- Confirm file existence
!
      if(gridFormat=='FV') then
        inquire(file=gridData,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'Grid data file ',trim(gridData),' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
!#ifdef GAMBIT
      else if(gridFormat=='GB') then
        inquire(file=gridData,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'Grid data file ',trim(gridData),' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
!#endif
!#ifdef SCRYU
      else if(gridFormat=='SY') then
        inquire(file=gridData,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'SCRYU file',trim(gridData),' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
!#endif
      else if(gridFormat=='FF') then
        inquire(file=gridData,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'FrontFlowRed Grid file',trim(gridData),
     &              ' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
      else if(gridFormat=='GF') then
        inquire(file=gridData,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'FrontFlowRed Grid file',trim(gridData),
     &              ' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
!#ifdef STARCD
      else if(gridFormat=='SC') then
        write(gridDataSC,*)trim(adjustl(gridData)),'.inp'
        gridDataSC=adjustl(gridDataSC)
        inquire(file=gridDataSC,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'STAR-CD Input file ',trim(gridDataSC),
     &                                          ' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
        write(gridDataSC,*)trim(adjustl(gridData)),'.vrt'
        gridDataSC=adjustl(gridDataSC)
        inquire(file=gridDataSC,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'STAR-CD Vertex file ',trim(gridDataSC),
     &                                          ' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
        write(gridDataSC,*)trim(adjustl(gridData)),'.cel'
        gridDataSC=adjustl(gridDataSC)
        inquire(file=gridDataSC,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'STAR-CD Cell file ',trim(gridDataSC),
     &                                          ' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
        write(gridDataSC,*)trim(adjustl(gridData)),'.bnd'
        gridDataSC=adjustl(gridDataSC)
        inquire(file=gridDataSC,exist=fileEx)
        if(.not.(fileEx)) then
          write(*,*)'STAR-CD Boundary file ',trim(gridDataSC),
     &                                          ' cannot be found.'
          write(*,*)'Check file name and path.'
          showUsage=.true.
        end if
!#endif

!
      end if
!
!     Show Usage
!
      if(showUsage) then
         write(*,*) 'Usage:'
         if(.not.animeStat) then
            write(*,*) '  ffr2viz -r result.frontflow -g gridfile',
     &           ' -o outfile [options...]'
         else
            write(*,*) '  ffrmovie -r anim -g gridfile -o outfile',
     &           ' [options...]'
         end if
         write(*,*) 'options:'
         write(*,*) ' -m (n) : multi processor result, CPU number'
!         write(*,*) ' [-s] : single processor result'
!         write(*,*) ' [-a] : marge and convert all processor data'
!         write(*,*) ' -p (n1 n2 ..) : convert only specified',
!     &        'processor data'
         write(*,*) ' -g : grid data'
         write(*,*) ' -gs (scale) : grid file scale'
         write(*,*) ' -gf [FV,GF,FF] : grid data format'
         write(*,*) ' -fm [A,U] : grid file format for -gf [GF].',
     &        ' (A:Ascii, U:Unformatted)'
         write(*,*) ' -o (output file name) : output file name prefix'
         write(*,*) ' -d (directory prefix) : subdir prefix'
         write(*,*) ' -r (file name)  : result file'
         write(*,*) ' -rf [FV,GF] : result file format (default:FF)'
         write(*,*) ' -nr : No result (only grid)'
         write(*,*) ' -pv : output CPU ID as value'
! DEBUG
         write(*,*) ' -vn : output vertex ID as value'
!
         write(*,*) ' -wz : set wall velocity to 0'
         write(*,*) ' -rv : create RPM reference [rpm_velocity]',
     &        'you can convert Absolute -> Relative velocity with this.'
!
!         write(*,*) ' -byp (wall_dist file) : output Y+ value'
! under develop
!         write(*,*) ' -av (ndiv) (Ret) (outfile) : y+ averaging.'
!         write(*,*) ' -wdf (wall_dist file) : read wall distance',
!     &        'from file'
!         write(*,*) ' -cyl (x,y,z) (YZXaxis) (ZXYaxis) (R) (scale)',
!     &        ': calculate wall distance as cylinder'
!         write(*,*) ' -led : y+ log scale division'
         write(*,*) ''
         stop
       end if

!#######################################################################
! animation loop
      animeStep=-1
      animeExt=''
      rsOpenStat=.true.
      rsCloseStat=.true.
      do while(animeStep==-1 .or. animeStat)
         animeStep=animeStep+1
         if(animeStep>0) rsOpenStat =.false.
         if(animeStat)   rsCloseStat=.false.
         if(animeStat) then
!           ex) RES_0000.fv: animeExt='0000'
            write(animeExt,'(I4.4)') animeStep
         end if
!-------------------------
! ---    confirm overwrite
!-------------------------
         ans='n'
         inqChk=.true.
         outFVprfx=trim(outFVprfx)
         len=len_trim(outFVprfx)
         l1 =len_trim(outExt(1))
         l2 =len_trim(outExt(2))
!        outFVprfx:
!          ex) -rf FV -o RES.fv -> RES
         if(l1>0 .and. len>l1) then
            if(outFVprfx(len-l1+1:len)==outExt(1)(1:l1)) then
               outFVprfx(len-l1+1:len)=' '
            end if
         end if
         if(l2>0 .and. len>l2) then
            if(outFVprfx(len-l2+1:len)==outExt(2)(1:l2)) then
               outFVprfx(len-l2+1:len)=' '
            end if
         end if
!
         do while(inqChk)
            if (resultFormat == 'FL' .or. resultFormat == 'FS') then
               if(ans=='n') then
!                 '.cas','.gfv'
                  outFVfn=trim(outFVprfx)//trim(animeExt)//outExt(1)
                  inqChk=.true.
               else
!                 '.dat','.rfv'
                  outFVfn=trim(outFVprfx)//trim(animeExt)//outExt(2)
                  inqChk=.false.
               end if
               if(gridOnly) then
                  inqChk=.false.
               end if
            elseif (resultFormat == 'ES') then
               if(ans=='n') then
!                 'geo (geo0000)'
                  outFVfn=trim(outFVprfx)//outExt(1)//trim(animeExt)
                  inqChk=.true.
               else
!                 'case'
                  outFVfn=trim(outFVprfx)//outExt(2)
                  inqChk=.false.
               end if
               if(gridOnly) then
                  inqChk=.false.
               end if
            elseif (resultFormat == 'SC') then
               outFVfn=trim(outFVprfx)//trim(animeExt)
     &                                //'_pressu'//outExt(1)
               inqChk=.false.
            else
!              '.gf','.fv', '.inp', '.fld'
               outFVfn=trim(outFVprfx)//trim(animeExt)//outExt(1)
               inqChk=.false.
            end if

            if(animeStep<=0) then ! only first time
               inquire(file=outFVfn,exist=fileEx)
               if(fileEx) then
                  write(*,*)'File ',trim(outFVfn),
     &                          ' exists. Overwrite OK? [y/n]'
                  read(*,*)ans
                  if((ans/='y').and.(ans/='Y')) stop
! animation
                  if(animeStat .or.
     &               resultFormat == 'SC') then   ! anime==.true.
                     write(*,*)'All series of file xxxx',
     &                    trim(animeExt)//trim(outExt(1)),
     &                    ' will be overwritten. Are you sure OK? [y/n]'
                     read(*,*)ans
                     if((ans/='y').and.(ans/='Y')) stop
                  end if
               else
                  ans='y'
               end if
            else
               ans='y'
            end if
         end do
!        *** note!
         if (resultFormat == 'FL' .or. resultFormat == 'SC' .or.
     &       resultFormat == 'FS') then
            outFVfn=trim(outFVprfx)//trim(animeExt)
         elseif(resultFormat == 'ES') then
            outFVfn=trim(outFVprfx)
         end if

         if(animeStep==0) then ! #######################################
!---------------------
! ---       Read grid file
!---------------------
            call read_griddata(gridData,gridFormat,ffgFormatKey,
     &                         resultFormat,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'Error occured in grid data'
               stop
            end if

!           grid scaling
            cord(:,:)=gdScale*cord(:,:)
! DEBUG
            allocate(cord_bak(1:3,1:nvrtx))
            cord_bak(:,:)=cord(:,:)
!
!           read metis part file 'test.m.part.x' for HPC
            call read_partdata(MPprocNum,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'Error occured in xxx.m.part.x data'
               stop
            end if
!
            call setupGridGrobalInfo(gridFormat,resultFormat,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'Error occured in grid global data'
               stop
            end if
!
!----------------------------
! ---       Reading wall distance
!----------------------------
            if(ypAverage) then
               if(wdisCalcType<0) then
                  write(*,*)'ABEND:lack of option'
                  stop
               end if
               if(MPflag) then
                  write(*,*)'MP calculation is not compatible.'
                  stop
!---------------------------------------------------------
! ---       Reading wall distance for HPC and max. distance
!---------------------------------------------------------
               else
                  select case(wdisCalcType)
                  case(0)          !Reading from file
                     filePath=wdisYPAVfn
                     OPEN(21,FILE=filePath,FORM='UNFORMATTED',
     &                    action='read',iostat=ios)
                     read(21)i,j
                     if((i/=1).or.(j/=nvrtx)) then
                        write(*,*)'Wall dist data file error.'
                        stop
                     end if
                     allocate(wallDist(1:nvrtx))
                     read(21)wallDist(:)
                     close(21)
                  case(2)       ! Calculation as cylinder
                     allocate(wallDist(1:nvrtx))
                     if((axisTag=='x').or.(axisTag=='X')) then
                        do i=1,nvrtx
                           wallDist(i)=Radius-sqrt((cord(2,i)-axis1)**2+
     &                                             (cord(3,i)-axis2)**2)
                           if(wallDist(i)<0.0) wallDist(i)=0.0
                        end do
                     else if((axisTag=='y').or.(axisTag=='Y')) then
                        do i=1,nvrtx
                           wallDist(i)=Radius-sqrt((cord(3,i)-axis1)**2+
     &                                             (cord(1,i)-axis2)**2)
                           if(wallDist(i)<0.0) wallDist(i)=0.0
                        end do
                     else if((axisTag=='z').or.(axisTag=='Z')) then
                        do i=1,nvrtx
                           wallDist(i)=Radius-sqrt((cord(1,i)-axis1)**2+
     &                                             (cord(2,i)-axis2)**2)
                           if(wallDist(i)<0.0) wallDist(i)=0.0
                        end do
                     end if
                     wallDist(:)=wallDist(:)*scaleFactor
                  case default
                     write(*,*)'Optin error.'
                     stop
                  end select
                  wdMax=real(maxval(wallDist))
                  wdMin=real(minval(wallDist))
               end if
!
! --- 
!
               if(ypLogDiv) then ! log
                  dyp=(log10(wdMax*ReTau*scaleFactor)-
     &                 log10(wdMin*ReTau*scaleFactor))/real(divnumYPAV)
                  allocate(ypVal(0:divnumYPAV))
                  do i=0,divnumYPAV
                     ypVal(i)=10.0**(
     &                    real(i)*dyp+log10(wdMin*ReTau*scaleFactor))
                  end do
                  ypVal(:)=ypVal(:)/(ReTau*scaleFactor)
                  ypVal(divnumYPAV)=wdMax
               else             ! Linear 
                  dyp=(wdMax)/real(divnumYPAV)
                  allocate(ypVal(0:divnumYPAV))
                  do i=0,divnumYPAV
                     ypVal(i)=real(i)*dyp
                  end do
                  ypVal(divnumYPAV)=wdMax
               end if
            end if
!
         end if  ! if(animeStep==0) then ! #############################
!
         write(*,*)' -----------------------------------------------'
!------------------------------
! ---    Read FrontFlow Results
!------------------------------
!
         if(gridOnly) then
            call setNumVarsToZero(MPprocNum,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'Error occured in setNumvarsToZero'
               stop
            end if
            animeStat=.false.
         else

!
! ---    Read result.frontflow
!
! DEBUG
            call readAndStoreRslFFR(MPprocNum,Rsltfn,hpcprefix,
     &                              procIDVals,vrtxIDVals,
     &                              boundaryYp,bndYpWallfn,
     &                              animeStat,rsOpenStat,rsCloseStat,
     &                              wallZero,RelativeVelocity,rtrnStat)
!            call readAndStoreRslFFR(MPprocNum,Rsltfn,hpcprefix,
!     &                              procIDVals,boundaryYp,bndYpWallfn,
!     &                              animeStat,rsOpenStat,rsCloseStat,
!     &                              wallZero,RelativeVelocity,rtrnStat)
!
! animation
            if(rtrnStat==-1) then
!              exit animation loop
               write(*,*) ' End reading animation file.'
!#ifdef ENSIGHT
               if(resultFormat == 'ES') then
                  write(*,*) ' #### Add TIME to EnSight Case...'
                  call writeENSIGHTHeader_TIME(outFVfn,25,
     &                                  animeStat,animeStep,rtrnStat)
               end if
!#endif
               exit
            end if
            if(rtrnStat/=0) then
               write(*,*)'ABEND:Cannot read and store values.'
               stop
            end if

         end if                 ! from .not.(gridOnly)


         write(*,*) ' #### output file:',trim(outFVfn)

!------------------------------
! ---    Write Output Files
!------------------------------
!        --- FIELDVIEW --------------
         if (resultFormat == 'FV') then
            write(*,*)' #### Writing FV header...'
            call writeFVHeader(outFVfn,25,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'ABEND:Cannot create output file.'
               stop
            end if
            write(*,*)' #### writing FV Header completed.'
!
            write(*,*)' #### Writing FV variable data ...'
            call append_resultfv(outFVfn,25,MPprocNum,gridOnly,rtrnStat)
            write(*,*)' #### Appending data for FV completed.'
         else if (resultFormat == 'VT') then
           write(*,*) '*** Writing vtk data ...'
           call ffr_to_vtk(outFVfn, 25, MPprocNum, gridOnly, rtrnStat,
     &                     .false.)
           if (rtrnStat .ne. 0) then
             write(*,*) 'err> ffr_to_vtk'
             stop
           end if
         else if (resultFormat == 'VA') then
           write(*,*) '*** Writing vtk data (ascii) ...'
           call ffr_to_vtk(outFVfn, 25, MPprocNum, gridOnly, rtrnStat,
     &                     .true.)
           if (rtrnStat .ne. 0) then
             write(*,*) 'err> ffr_to_vtk'
             stop
           end if
!#ifdef FVASCII
!        --- FIELDVIEW --------------
         else if (resultFormat == 'FA') then
            write(*,*)' #### Writing FV(ASCII) header...'
            call writeFVHeader_ASCII(outFVfn,25,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'ABEND:Cannot create output file.'
               stop
            end if
            write(*,*)' #### writing FV(ASCII) Header completed.'
!
            write(*,*)' #### Writing FV(ASCII) variable data ...'
            call append_resultfv_ASCII(outFVfn,25,MPprocNum,
     &                                 gridOnly,rtrnStat)
            write(*,*)' #### Appending data for FV(ASCII) completed.'
!#endif
         else if (resultFormat == 'FS') then
            write(*,*)' #### Writing FV(SPLIT) header...'
            call writeFVHeader_SPLIT(outFVfn,25,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'ABEND:Cannot create output file.'
               stop
            end if
            write(*,*)' #### writing FV(SPLIT) Header completed.'
!
            write(*,*)' #### Writing FV(SPLIT) variable data ...'
            call append_resultfv_SPLIT(outFVfn,25,MPprocNum,
     &                                 gridOnly,rtrnStat)
            write(*,*)' #### Appending data for FV(SPLIT) completed.'
!#ifdef AVS
!        --- AVS --------------
         else if (resultFormat == 'AA') then
!
            write(*,*)' #### Writing AVS header...'
            call writeAVSHeader(outFVfn,25,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'ABEND:Cannot create output file.'
               stop
            end if
            write(*,*)' #### writing AVS Header completed.'
!
            write(*,*)' #### Writing AVS variable data ...'
            call append_resultavs(outFVfn,25,MPprocNum,
     &           targProc,gridOnly,gridFormat,rtrnStat)
            if(rtrnStat/=0) then
               write(*,*)'ABEND:Cannot append to output file.'
               stop
            end if
            write(*,*)' #### Appending data for AVS completed.'
!#endif
!
!#ifdef FLUENT
!        --- FLUENT --------------
         else if (resultFormat == 'FL') then

         else if (resultFormat == 'SY') then
         else if (resultFormat == 'SC') then
         else if (resultFormat == 'ES') then
         end if
!
! ------------------------------------------
!      
!        for next animation step
         if(animeStat) then ! anime==.true.
            deallocate(nameVars,totdata)
            deallocate(icomp_avex,icomp_rmsx)
! wall variables
!            deallocate(nameBVars,bnddata)
         end if

      end do ! do while(animeStep==-1 .or. animeStat)
!####################################################################

!  Y+/u+
      if(ypAverage) then
         call calculateYpUpAv(27,rtrnStat)
         if(rtrnStat/=0) then
            write(*,*)'ABEND:Cannot average uplus.'
            stop
         end if
      end if

      deallocate(targProc)
      deallocate(cord_bak)

      write(*,*) ' END: ffr2viz conversion is successfully completed.'
!
      ierror=0
      return

      end subroutine main

!      contains
!-------------------------------------------------------------------
      function hpcPath(procID,fileName,hpcprefix) result(path)
!-------------------------------------------------------------------
      implicit none
!
      character(*) :: fileName
      character(*) :: hpcprefix
      character(*) :: path
      integer :: procID

 5010 format(a,a,i4.4,a,a)
      write(path,5010) trim(adjustl(hpcprefix)),'_',procID,
     &     '/',trim(adjustl(fileName))
      path=adjustl(path)
      end function hpcPath
