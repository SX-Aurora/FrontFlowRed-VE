      SUBROUTINE OPEN_F(IUNIT,IMODE,FNAME,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM OPEN THE SPECIFIED FILE BY FORTRAN UNFORMATTED MODE  
C    
C    [IN]
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    CHAR+60 FNAME :FILE NAME 
C
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
C
      INTEGER IUNIT
      INTEGER IMODE
      CHARACTER*60 FNAME
      INTEGER IERROR
C
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.EQ.ASCII)THEN
        OPEN(IUNIT,FILE=FNAME,FORM='FORMATTED'
     &        ,ERR=100)
      ENDIF
      IF(IMODE.EQ.UNFMT)THEN
        OPEN(IUNIT,FILE=FNAME,FORM='UNFORMATTED'
     &       ,ERR=100)
      ENDIF
C
      RETURN
C
 100  CONTINUE
      IERROR=1    
      RETURN
      END
C
C
C
      SUBROUTINE CLOSE_F(IUNIT,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM CLOSE THE SPECIFIED FILE BY FORTRAN UNFORMATTED MODE  
C    
C    [IN]
C    INT*4   IUNIT :FILE BUFFER NUMBER
C
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
C
      INTEGER INUM
      INTEGER IERROR
C
      CLOSE(IUNIT,ERR=100)
      RETURN
C
 100  CONTINUE
      IERROR=1    
      RETURN
      END
C
C
C
      SUBROUTINE RW_C1_08(IRW,IMODE,IUNIT,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE CHARCTER*8 ARGUMENT BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX2  :SIZE OF ARRAY
C    INT*4   NUM2  :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    CHAR*8  VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      CHARACTER*8 VAR
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
          IF(IMODE.EQ.ASCII)READ(IUNIT,FMTAS,IOSTAT=IS,ERR=110)  VAR
          IF(IMODE.EQ.UNFMT)READ(IUNIT,      IOSTAT=IS,ERR=110)  VAR
 110  CONTINUE
      ENDIF
      IF(IRW.EQ.2)THEN
          IF(IMODE.EQ.ASCII)WRITE(IUNIT,FMTAS)  VAR
          IF(IMODE.EQ.UNFMT)WRITE(IUNIT      )  VAR
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_C1_30(IRW,IMODE,IUNIT,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE CHARCTER*30 ARGUMENT BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX2  :SIZE OF ARRAY
C    INT*4   NUM2  :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    CHAR*30 VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      CHARACTER*30 VAR
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
         IF(IMODE.EQ.ASCII) READ(IUNIT,FMTAN) VAR
         IF(IMODE.EQ.UNFMT) READ(IUNIT      ) VAR
      ENDIF
      IF(IRW.EQ.2)THEN
         IF(IMODE.EQ.ASCII) WRITE(IUNIT,FMTAN) VAR
         IF(IMODE.EQ.UNFMT) WRITE(IUNIT      ) VAR
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_C1_60(IRW,IMODE,IUNIT,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE CHARCTER*60 ARGUMENT BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX2  :SIZE OF ARRAY
C    INT*4   NUM2  :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    CHAR*60 VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      CHARACTER*60 VAR
      INTEGER     IERROR
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
         IF(IMODE.EQ.ASCII) READ(IUNIT,FMTAM) VAR
         IF(IMODE.EQ.UNFMT) READ(IUNIT      ) VAR
      ENDIF
      IF(IRW.EQ.2)THEN
         IF(IMODE.EQ.ASCII) WRITE(IUNIT,FMTAM) VAR
         IF(IMODE.EQ.UNFMT) WRITE(IUNIT      ) VAR
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_I1(IRW,IMODE,IUNIT,MAX,NUM,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE INTEGER*4 ARRAY BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    INT*4   VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      INTEGER     VAR(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) READ(IUNIT,FMTI) (VAR(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT) READ(IUNIT     ) (VAR(I),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) WRITE(IUNIT,FMTI) (VAR(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT) WRITE(IUNIT     ) (VAR(I),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_I2(IRW,IMODE,IUNIT,MAX,NUM,VAR1,VAR2,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE TWO INTEGER*4 ARRAYS BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    INT*4   VAR1   :ARRAY(1) TO READ-IN OR WRITTEN
C    INT*4   VAR2   :ARRAY(2) TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      INTEGER     VAR1(MAX)
      INTEGER     VAR2(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) 
     &  READ(IUNIT,FMTI)  (VAR1(I),VAR2(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT)  
     &  READ(IUNIT     )  (VAR1(I),VAR2(I),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) 
     &  WRITE(IUNIT,FMTI)  (VAR1(I),VAR2(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT)  
     &  WRITE(IUNIT     )  (VAR1(I),VAR2(I),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_I3(IRW,IMODE,IUNIT,MAX,NUM,VAR1,VAR2,VAR3,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE THREE INTEGER*4 ARRAYS BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    INT*4   VAR1   :ARRAY(1) TO READ-IN OR WRITTEN
C    INT*4   VAR2   :ARRAY(2) TO READ-IN OR WRITTEN
C    INT*4   VAR3   :ARRAY(3) TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      INTEGER     VAR1(MAX)
      INTEGER     VAR2(MAX)
      INTEGER     VAR3(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) 
     &  READ(IUNIT,FMTI) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)  
        IF(IMODE.EQ.UNFMT)  
     &  READ(IUNIT     ) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)   
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) 
     &  WRITE(IUNIT,FMTI) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)  
        IF(IMODE.EQ.UNFMT)  
     &  WRITE(IUNIT     ) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)   
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_IN(IRW,IMODE,IUNIT,MAX,MAX2,NUM,NUM2,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE INTEGER*4 2D-ARRAY BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   MAX2  :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C    INT*4   NUM2  :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    INT*4   VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX,MAX2
      INTEGER     NUM,NUM2
      INTEGER     VAR(MAX2,MAX),BUF
      INTEGER     IERROR
!
      INTEGER I,J
!
      IERROR=0
!
      IF(IRW.NE.0.AND.IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.NE.0.AND.NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.0)THEN
         IF(IMODE.EQ.ASCII) 
     &   READ(IUNIT,FMTI) ((BUF,J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   READ(IUNIT     ) ((BUF,J=1,NUM2),I=1,NUM)
      ENDIF
      IF(IRW.EQ.1)THEN
         IF(IMODE.EQ.ASCII) 
     &   READ(IUNIT,FMTI) ((VAR(J,I),J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   READ(IUNIT     ) ((VAR(J,I),J=1,NUM2),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
         IF(IMODE.EQ.ASCII) 
     &   WRITE(IUNIT,FMTI) ((VAR(J,I),J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   WRITE(IUNIT     ) ((VAR(J,I),J=1,NUM2),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_F1(IRW,IMODE,IUNIT,MAX,NUM,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE REAL*4 ARRAY BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*4  VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      REAL        VAR(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) READ(IUNIT,FMTF) (VAR(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT) READ(IUNIT     ) (VAR(I),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) WRITE(IUNIT,FMTF) (VAR(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT) WRITE(IUNIT     ) (VAR(I),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_F2(IRW,IMODE,IUNIT,MAX,NUM,VAR1,VAR2,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE TWO REAL*4 ARRAYS BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*4  VAR1   :ARRAY(1) TO READ-IN OR WRITTEN
C    REAL*4  VAR2   :ARRAY(2) TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      REAL        VAR1(MAX)
      REAL        VAR2(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) 
     &  READ(IUNIT,FMTF)  (VAR1(I),VAR2(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT)  
     &  READ(IUNIT     )  (VAR1(I),VAR2(I),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) 
     &  WRITE(IUNIT,FMTF)  (VAR1(I),VAR2(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT)  
     &  WRITE(IUNIT     )  (VAR1(I),VAR2(I),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_F3(IRW,IMODE,IUNIT,MAX,NUM,VAR1,VAR2,VAR3,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE THREE REAL*4 ARRAYS BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*4  VAR1   :ARRAY(1) TO READ-IN OR WRITTEN
C    REAL*4  VAR2   :ARRAY(2) TO READ-IN OR WRITTEN
C    REAL*4  VAR3   :ARRAY(3) TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      REAL        VAR1(MAX)
      REAL        VAR2(MAX)
      REAL        VAR3(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.NE.0.AND.NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) 
     &  READ(IUNIT,FMTF) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)  
        IF(IMODE.EQ.UNFMT)  
     &  READ(IUNIT     ) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)   
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) 
     &  WRITE(IUNIT,FMTF) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)  
        IF(IMODE.EQ.UNFMT)  
     &  WRITE(IUNIT     ) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)   
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_FN(IRW,IMODE,IUNIT,MAX,MAX2,NUM,NUM2,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE REAL*4 2D-ARRAY BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   MAX2  :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C    INT*4   NUM2  :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*4  VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX,MAX2
      INTEGER     NUM,NUM2
      REAL        VAR(MAX2,MAX),BUF
      INTEGER     IERROR
!
      INTEGER I,J
!
      IERROR=0
!
      IF(IRW.NE.0.AND.IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.0)THEN
         IF(IMODE.EQ.ASCII) 
     &   READ(IUNIT,FMTF) ((BUF,J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   READ(IUNIT     ) ((BUF,J=1,NUM2),I=1,NUM)
      ENDIF
      IF(IRW.EQ.1)THEN
         IF(IMODE.EQ.ASCII) 
     &   READ(IUNIT,FMTF) ((VAR(J,I),J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   READ(IUNIT     ) ((VAR(J,I),J=1,NUM2),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
         IF(IMODE.EQ.ASCII) 
     &   WRITE(IUNIT,FMTF) ((VAR(J,I),J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   WRITE(IUNIT     ) ((VAR(J,I),J=1,NUM2),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_D1(IRW,IMODE,IUNIT,MAX,NUM,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE REAL*8 ARRAY BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*8  VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      REAL*8      VAR(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) READ(IUNIT,FMTD) (VAR(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT) READ(IUNIT     ) (VAR(I),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) WRITE(IUNIT,FMTD) (VAR(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT) WRITE(IUNIT     ) (VAR(I),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_D2(IRW,IMODE,IUNIT,MAX,NUM,VAR1,VAR2,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE TWO REAL*8 ARRAYS BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*8  VAR1   :ARRAY(1) TO READ-IN OR WRITTEN
C    REAL*8  VAR2   :ARRAY(2) TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      REAL*8      VAR1(MAX)
      REAL*8      VAR2(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.NE.0.AND.NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) 
     &  READ(IUNIT,FMTD)  (VAR1(I),VAR2(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT)  
     &  READ(IUNIT     )  (VAR1(I),VAR2(I),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) 
     &  WRITE(IUNIT,FMTD)  (VAR1(I),VAR2(I),I=1,NUM)
        IF(IMODE.EQ.UNFMT)  
     &  WRITE(IUNIT     )  (VAR1(I),VAR2(I),I=1,NUM)
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_D3(IRW,IMODE,IUNIT,MAX,NUM,VAR1,VAR2,VAR3,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE THREE REAL*8 ARRAYS BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*8  VAR1   :ARRAY(1) TO READ-IN OR WRITTEN
C    REAL*8  VAR2   :ARRAY(2) TO READ-IN OR WRITTEN
C    REAL*8  VAR3   :ARRAY(3) TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX2
      INTEGER     NUM2
      REAL*8      VAR1(MAX)
      REAL*8      VAR2(MAX)
      REAL*8      VAR3(MAX)
      INTEGER     IERROR
!
      INTEGER I
!
      IERROR=0
!
      IF(IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.1)THEN
        IF(IMODE.EQ.ASCII) 
     &  READ(IUNIT,FMTD) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)  
        IF(IMODE.EQ.UNFMT)  
     &  READ(IUNIT     ) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)   
      ENDIF
      IF(IRW.EQ.2)THEN
        IF(IMODE.EQ.ASCII) 
     &  WRITE(IUNIT,FMTD) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)  
        IF(IMODE.EQ.UNFMT)  
     &  WRITE(IUNIT     ) (VAR1(I),VAR2(I),VAR3(I),I=1,NUM)   
      ENDIF
!
      RETURN
      END
C
C
C
      SUBROUTINE RW_DN(IRW,IMODE,IUNIT,MAX,MAX2,NUM,NUM2,VAR,IERROR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C    BASIC INPUT-OUTPUT  UTILITY FOR  GENERAL FILE VERSION 2.1              
C
C    THIS PROGRAM READ OR WRITE REAL*8 2D-ARRAY BY FORTARN 
C    UNFORMATTED  MODE.    
C    
C    [IN]
C    INT*4   IRW   :FLAR FOR ACTION(1:READ,2:WRITE)
C    INT*4   IMODE :TYPE OF INPUT OR OUTPUT MODE (ASCII OR UNFMT)
C    INT*4   IUNIT :FILE BUFFER NUMBER
C    INT*4   MAX   :SIZE OF ARRAY
C    INT*4   MAX2  :SIZE OF ARRAY
C    INT*4   NUM   :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C    INT*4   NUM2  :NUMBER OF ARRAY TO READ-IN OR WRITTEN
C
C    [INOUT]
C    REAL*8  VAR   :ARRAY TO READ-IN OR WRITTEN
C   
C    [OUT]    
C    INT*4   IERROR:ERROR CODE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      INCLUDE 'gf.h'
!
      INTEGER     IRW 
      INTEGER     IMODE
      INTEGER     IUNIT
      INTEGER     MAX,MAX2
      INTEGER     NUM,NUM2
      REAL*8      VAR(MAX2,MAX),BUF
      INTEGER     IERROR
!
      INTEGER I,J
!
      IERROR=0
!
      IF(IRW.NE.0.AND.IRW.NE.1.AND.IRW.NE.2)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(NUM.GT.MAX)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IMODE.NE.ASCII.AND.IMODE.NE.UNFMT)THEN
        IERROR=1
        RETURN 
      ENDIF 
!
      IF(IRW.EQ.0)THEN
         IF(IMODE.EQ.ASCII) 
     &   READ(IUNIT,FMTD) ((BUF,J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   READ(IUNIT     ) ((BUF,J=1,NUM2),I=1,NUM)
      ENDIF
      IF(IRW.EQ.1)THEN
         IF(IMODE.EQ.ASCII) 
     &   READ(IUNIT,FMTD) ((VAR(J,I),J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   READ(IUNIT     ) ((VAR(J,I),J=1,NUM2),I=1,NUM)
      ENDIF
      IF(IRW.EQ.2)THEN
         IF(IMODE.EQ.ASCII) 
     &   WRITE(IUNIT,FMTD) ((VAR(J,I),J=1,NUM2),I=1,NUM)
         IF(IMODE.EQ.UNFMT) 
     &   WRITE(IUNIT     ) ((VAR(J,I),J=1,NUM2),I=1,NUM)
      ENDIF
!
      RETURN
      END
!=======================================================================
      SUBROUTINE PRINTSTD(IOUT,MESSAGE)
!=======================================================================
      INTEGER IOUT
      CHARACTER*90  MESSAGE
 
      WRITE(IOUT,'(A90)')MESSAGE
      RETURN
      END
