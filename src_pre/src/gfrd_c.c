/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*    BASIC INPUT  UTILITY FOR GENERAL FILE VERSION 2.1                     */
/*                                                                          */
/*      Written by  Y. Yamade                                               */
/*      Since       9 Oct. 2002                                             */
/*      Reference   gfrd.f of LES3D                                         */
/*                                                                          */ 
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*    ARGUMENT LISTINGS                                                     */
/*      (1) INPUT                                                           */
/*INT*4    ITYPE       ; SPECIFIES THE TYPE OF ACTION TO BE TAKEN           */
/*           OPEN    --- OPEN THE SPECIFIED FILE, READ COMMENT STRINGS      */
/*                      OF THE FILE, RETURN READING MODE TO BE USED TO      */
/*                      READ THE FOLLOWING CONTENT OF THE FILE, THE         */
/*                      NUMBER OF THE COMMENT STRINGS AND THE COMMENT       */
/*                      STRINGS READ IN                                     */
/*           CLOSE   --- CLOSE THE FILE                                     */
/*           THE FOLLOWING THREE TYPES (IE. ANY, INQUIR, AND NEWSET)        */
/*          CAN BE SPECIFIED ONLY WHEN THE CURRENT FILE POSITION IS         */
/*          JUST BEFORE A HEADDER (IE. IFPOS=HEAD, SEE IFPOS ARGUMENT       */
/*          IN THE INPUT-OUTPUT SECTION FOR DETAIL OF THE FILE              */
/*          POSITION).                                                      */
/*           ANY     --- READ THE NEXT CONTENT OF THE FILE, REGARDLESS      */
/*                      OF ITS TYPE, RETURN THE TYPE OF THE CONTENT IN      */
/*                      ITSELF (ITYPE ARGUMENT).                            */
/*                       IF THE CONTENT IS A NEW SET, READ AND RETURN       */
/*                      THE NUMBER OF THE COMMENT STRINGS AND THE           */
/*                      COMMENT STRINGS READ IN. IF THE CONTENT IS AN       */
/*                      ARRAY DATA, READ AND RETURN THE KEYWORD, DATA       */
/*                      NAME AND DATA ITSELF. IF NO CONTENT EXISTS          */
/*                      BEYOND THE CURRENT FILE POSITION, ONLY ITYPE        */
/*                      ARGUMENT WILL BE RETURNED SET TO CLOSE.             */
/*           INQUIR  --- INQUIRE THE NEXT CONTENT OF THE FILE,              */
/*                      RETURN THE TYPE OF THE CONTENT IN ITSELF            */
/*                      (ITYPE ARGUMENT).                                   */
/*                       IF THE CONTENT IS A NEW SET, READ AND RETURN       */
/*                      THE NUMBER OF THE COMMENT STRINGS AND THE           */
/*                      COMMENT STRINGS READ IN. IF THE CONTENT IS AN       */
/*                      ARRAY DATA, READ AND RETURN THE KEYWORD, DATA       */
/*                      NAME.                                               */
/*                       IF NO CONTENT EXISTS BEYOND THE CURRENT FILE       */
/*                      POSITION, ONLY ITYPE ARGUMENT WILL BE RETURNED      */
/*                      SET TO CLOSE.                                       */
/*                       TO READ THE FOLLOWING DATA, CALL THIS ROUTINE      */
/*                      AGAIN WITH ITYPE ARGUMENT AS RETURNED.              */
/*                       TO SKIP THE FOLLOWING DATA, CALL THIS ROUTINE      */
/*                      AGAIN WITH ITYPE ARGUMENT WITH ITS SIGN CHANGED.    */
/*           NEWSET  --- SEEK THE FILE FROM THE CURRENT FILE POSITION       */
/*                      UNTIL A NEW SET OR THE END FILE IS ENCOUNTERED.     */
/*                       IF A NEW SET IS ENCOUNTERED, READ AND RETURN       */
/*                      THE NUMBER OF THE COMMENT STRINGS AND THE           */
/*                      COMMENT STRINGS READ IN.                            */
/*                       IF THE END FILE IS ENCOUNTERED, AN ERROR WILL      */
/*                      BE REPORTED AND THE SEQUENCE WILL BE RETURNED       */
/*                      WITH ITYPE ARGUMENT SET TO CLOSE.                   */
/*           THE FOLLOWING TYPES (IE. INTARY ,FLTARY AND DBLARY) CAN BE     */
/*          SPECIFIED TO READ AN ARRAY DATA OF THE TYPE DESCRIBED TO THE    */
/*          RIGHT. THE FILE POSITION CAN BE EITHER JUST BEFORE A HEADDER    */
/*          OR  DATA (IE. IFPOS=HEAD OR DATA, SEE IFPOS ARGUMENT IN THE     */
/*          INPUT-OUTPUT SECTION FOR DETAIL OF THE FILE POSITION).          */
/*          WITH THIS PARAMETER SPECIFICATION, THIS SUBROUTINE FUNCTIONS    */
/*          DIFFERENTLY DEPENDING ON THE CURRENT FILE POSITION SPECIFIED    */
/*          BY IFPOS ARGUMENT.                                              */
/*           IF THE CURRENT FILE POSITION IS JUST BEFORE A HEADDER,         */
/*          SEEK THE FILE FROM THE CURRENT FILE POSITION UNTIL AN ARRAY     */
/*          DATA OF THE SEPCIFIED TYPE WITH SPECIFIED KEYWORD, OR THE       */
/*          END FILE, IS ENCOUNTERED. IF THE SPECIFIED ARRAY DATA IS        */
/*          ENCOUNTERED, READ AND RETURN THE DATA. IF THE END FILE IS       */
/*          ENCOUNTERED, AN ERROR WILL BE REPORTED AND THE SEQUENCE         */
/*          WILL BE RETURNED WITH ITYPE ARGUMENT SET TO CLOSE.              */
/*           IF THE CURRENT FILE POSITION IS JUST BEFORE DATA, READ AND     */
/*          RETURN THE DATA AS DESCRIBED BY ITYPE ARGUMENT.                 */
/*           IF THESE TYPES (IE. ITYPE=INTARY ,FLTARY OR DBLARY) ARE        */
/*          SPECIFIED WITH ITS SIGN CHANGED, THIS SUBROUTINE FUNCTIONS THE  */
/*          SAME AS THE ORIGINAL SIGN, EXCEPT THAT FOR NEGATIVE SIGN THE    */
/*          DATA READING WILL BE SKIPPED.                                   */
/*           INTARY   --- READ  KEYWORD, DATA NAME, INTEGER ARRAY DATA      */
/*           FLTARY   --- READ  KEYWORD, DATA NAME, FLOAT   ARRAY DATA      */
/*           DBLARY   --- READ  KEYWORD, DATA NAME, DOUBLE  ARRAY DATA      */
/*    NOTES; IF EITHER NUM2 OR NUM ARGUMENTS (DESCRIBED LATER) IS LESS      */
/*          THAN 1, NO DATA WILL BE STORED TO ANY VARIABLE/ARRAYS EXCEPT    */
/*          THE DATA HEADDER, DATA NAME, AND THE ELEMENT NUMBERS.           */
/*INT*4    IDEST       ; SPECIFIES DESTINATION VARIABLE/ARRAYS TO STORE     */
/*                      THE ARRAY DATA READ IN FROM THE FILE WHEN           */
/*                      READING AN ARRAY DATA                               */
/*           ARY2D    --- THE DATA READ IN FROM THE FILE WILL BE STORED     */
/*                       TO 2-D ARRAY PASSED THROUGH IARY2D OR FARY2D       */
/*                       ARGUMENTS.                                         */
/*           ARY1D    --- THE DATA READ IN FROM THE FILE WILL BE STORED     */
/*                      TO 1-D ARRAYS PASSED THROUGH IARY1, IARY2, IARY3    */
/*                      OR FARY1, FARY2, FARY3 ARGUMENTS.                   */
/*           ARY0D    --- THE DATA READ IN FROM THE FILE WILL BE STORED     */
/*                      TO VARIABLE PASSED THROUGH IVAR OR FVAR             */
/*                      ARGUMENTS.                                          */
/*    NOTES; FOLLOWING CASES LISTED BELOW WILL BE REGARDED AS AN ERROR.     */
/*          IN SUCH CASES, NO DATA WILL BE STORED TO ANY VARIABLE/ARRAYS    */
/*          EXCEPT THE DATA HEADDER, DATA NAME, AND THE ELEMENT NUMBERS.    */
/*           (1) NUM2 > MAX2 OR NUM > MAX FOR ALL IDEST VALUES.             */
/*           (2) NUM2 > 3                 FOR IDEST = ARY1D                 */
/*           (3) NUM2 > 1    OR NUM > 1   FOR IDEST = ARY0D                 */
/*INT*4    IMODE       ; SPECIFIES THE TYPE OF INPUT  MODE                  */
/*           ASCII   --- ASCII               MODE                           */
/*           UNFMT   --- FORTRAN UNFORMATTED MODE                           */
/*           BINARY  --- BINARY              MODE ( NOT SUPPORTED )         */
/*    NOTES; THE CALLING PROGRAM UNIT DOES NOT HAVE TO SPECIFY IMODE        */
/*          ARGUMENT. THIS SUBROUTINE WILL AUTOMATICALLY DETERMINE, WHEN    */
/*          IT OPENS A FILE, THE INPUT MODE THAT SHOULD BE TAKEN TO READ    */
/*          THE FOLLOWING CONTENTS OF THE FILE. HENCE, THE CALLING          */
/*          PROGRAM UNIT SHOULD USE THE VALUE OF IMODE ARGUMENT THAT        */
/*          IS RETURNED WHEN OPENING A FILE, AND SHOULD NOT MODIFY THE      */
/*          RETURNED VALUE FOR THE SUBSEQUENT READ.                         */
/*INT*4    IUNIT       ; FILE BUFFER NUMBER TO BE USED TO ACCESS            */
/*                      THE SPECIFIED FILE                                  */
/*CHAR*(*) FILENM      ; FILE NAME TO BE OPENED                             */
/*    NOTES; FILENM ARGUMENT MUST BE SPECIFIED ONLY WHEN OPENING            */
/*          A FILE (IE. ITYPE=OPEN).                                        */
/*CHAR*8   KEYWRD     ; A CHARACTER STRING IDENTIFYING THE DATA TO BE       */
/*                     READ IN WHEN (ITYPE=INTARY ,FLTARY OR DBLARY) AND    */
/*                     IFPOS=HEAD.                                          */
/*    NOTES; KEYWRD ARGUMENT MUST BE SPCIFIED ONLY WHEN SEARCHING FOR AN    */
/*          ARRAY DATA (IE. ITYPE=INTARY OR FLTARY WITH IFPOS SET TO        */
/*          HEAD). IN ALL THE OTHER CASES, KEYWRD ARGUMENT WILL BE          */
/*          READ AND RETURNED IF THE READ IN DATA SET HAS AN ASSOCIATED     */
/*          KEYWORD.                                                        */
/*INT*4    MAX         ; THE FIRST  DIMENSION SIZE OF ARRAYS PASSED         */
/*INT*4    MAX2        ; THE SECOND DIMENSION SIZE OF ARRAYS PASSED,        */
/*                       THE FIRST  DIMENSION SIZE OF CHARACTER STRING      */
/*                      ARRAY COMMNT                                        */
/*         IUT0        ; FILE   NUMBER TO REPORT ERROR OCCURENCE            */
/*      (2) OUTPUT                                                          */
/*INT*4    NUM         ; FIRST  ELEMENT NUMBER OF ARRAY DATA READ IN        */
/*INT*4    NUM2        ; SECOND ELEMENT NUMBER OF ARRAY DATA READ IN        */
/*                      WHEN AN ARRAY DATA WAS READ IN                      */
/*                       NUMBER OF COMMENT STRINGS READ IN                  */
/*                      WHEN A FILE WAS OPENED OR A NEW SET WAS READ IN     */
/*    NOTES; NUM, AND NUM2 ARGUMENTS WILL BE RETURNED WHEN AN ARRAY DATA    */
/*          WAS READ IN (IE. ITYPE=INTARY ,FLTARY DBLARY,OR ITYPE=ANY       */
/*          AND THE NEXT CONTENT IS AN ARRAY DATA).                         */
/*           NUM2 ARGUMENT WILL ALSO BE RETURNED WHEN COMMENT STRINGS       */
/*          WERE READ IN (IE. ITYPE=OPEN, NEWSET, OR ITYPE=ANY, INQUIR      */
/*          AND THE NEXT CONTENT IS A NEW SET).                             */
/*CHAR*30  NAME        ; A CHARACTER STRING, INDICATING DATA NANE           */
/*    NOTES; NAME ARGUMENT WILL BE RETURNED WHEN AN ARRAY DATA SET WAS      */
/*          READ IN (IE. ITYPE=INTARY OR FLTARY OR DBLARY, OR ITYPE=ANY     */
/*           AND THE NEXT CONTENT IS AN ARRAY DATA).                        */
/*CHAR*60  COMMNT(MAX2); AN 1-D ARRAY OF CHARACTER STRINGS, INDICATING      */
/*                      THE GENERAL INFORMATION ON THE FILE OR A SET        */
/*    NOTES; COMMNT ARGUMENT WILL BE RETURNED ONLY WHEN OPENING A FILE      */
/*          (IE. ITYPE=OPEN), SEARCHING A NEW SET(IE. ITYPE=NEWSET),        */
/*          READING OR INQURING THE NEXT CONTENT OF THE FILE (IE.           */
/*          ITYPE=ANY OR INQUIR) AND THE CONTENT IS A NEW SET.              */
/*INT*4    IVAR        ; AN INTEGER VARIABLE  READ IN FROM THE FILE         */
/*                      WHEN ITYPE=INTARY, IDEST=ARY0D, AND NUM2=NUM=1      */
/*INT*4    IARY1(NUM)  ; AN 1-D INTEGER ARRAY READ IN FROM THE FILE         */
/*                      WHEN ITYPE=INTARY, IDEST=ARY1D, AND NUM2=(1,2,3)    */
/*INT*4    IARY2(NUM)  ; AN 1-D INTEGER ARRAY READ IN FROM THE FILE         */
/*                      WHEN ITYPE=INTARY, IDEST=ARY1D, AND NUM2=(  2,3)    */
/*INT*4    IARY3(NUM)  ; AN 1-D INTEGER ARRAY READ IN FROM THE FILE         */
/*                      WHEN ITYPE=INTARY, IDEST=ARY1D, AND NUM2=(    3)    */
/*INT*4    IARY2D(MAX2,; A  2-D INTEGER ARRAY READ IN FROM THE FILE         */
/*                 NUM) WHEN ITYPE=INTARY, IDEST=ARY2D                      */
/*REAL*4   FVAR        ; AN FLOAT(DOUBLE) VARIABLE  READ IN FROM THE FILE   */
/*                      WHEN ITYPE=FLTARY OR DBLARY,                        */
/*                      IDEST=ARY0D, AND NUM2=NUM=1                         */
/*REAL*4   FARY1(NUM)  ; AN 1-D FLOAT(DOUBLE)  ARRAY READ IN FROM THE FILE  */
/*                      WHEN ITYPE=FLTARY  OR DBLARY,                       */
/*                      IDEST=ARY1D, AND NUM2=(1,2,3)                       */
/*REAL*4   FARY2(NUM)  ; AN 1-D FLOAT(DOBLE) ARRAY READ IN FROM THE FILE    */
/*                      WHEN ITYPE=FLTARY  OR DBLARY,                       */
/*                      IDEST=ARY1D, AND NUM2=(2,3)                         */
/*REAL*4   FARY3(NUM)  ; AN 1-D FLOAT(DOUBLE) ARRAY READ IN FROM THE FILE   */
/*                      WHEN ITYPE=FLTARY  OR DBLARY,                       */
/*                      IDEST=ARY1D, AND NUM2=(3)                           */
/*REAL*4   FARY2D(MAX2,; A  2-D FLOAT(DOUBLE) ARRAY READ IN FROM THE FILE   */
/*                 NUM) WHEN ITYPE=FLTARY OD DBLARY,IDEST=ARY2D             */
/*INT*4    IERR        ; RETURN CODE WHOSE VALUE WILL BE EITHER             */
/*           0       --- INDICATING SUCCESSFUL TERMINATION                  */
/*       OR  1       --- INDICATING OCCURENCE OF SOME ERROR CONDITIONS      */
/*      (3) INPUT-OUPUT                                                     */
/*INT*4    IFPOS       ; SPECIFIES CURRENT FILE POSITION, TAKING ONE OF     */
/*                      THE FOLLOWING VALUES:                               */
/*           HEAD    --- JUST BEFORE A HEADDER                              */
/*           DATA    --- JUST BEFORE DATA                                   */
/*           END     --- JUST AFTER  THE END OF THE FILE                    */
/*           UNKNWN  --- UNKNOWN POSITION                                   */
/*    NOTES; AFTER EACH CALL TO THIS SUBROUTINE, IFPOS ARGUMENT IS          */
/*          SET TO AN APPROPRIATE VALUE. HENCE, THE CALLING PROGRAM         */
/*          UNIT DOES NOT HAVE TO SPECIFY THIS ARGUMENT, UNLESS ANOTHER     */
/*          SEQUENCE HAS CHANGED THE FILE POSITION AFTER THE LAST CALL      */
/*          TO THIS SUBROUTINE.                                             */
/*SUMMARY OF THE FUNCTIONS SUPPORTED IN THIS SUBROUTINE                     */
/*+===================================================================+     */
/*|   VALUE OF ITYPE    |    VALUE OF IFPOS   |                       |     */
/*+-------------------------------------------|      FUNCTIONS        |     */
/*| CALLING  | RETURNED | CALLING  | RETURNED |                       |     */
/*+===================================================================+     */
/*| OPEN     | OPEN     | UNDEFINED| HEAD     | OPEN FILE, RETURN MODE|     */
/*|          |          |          |          | READ COMMENNT STRINGS |     */
/*+-------------------------------------------------------------------+     */
/*| CLOSE    | CLOSE    | ANY      | UNCHANGED| CLOSE FILE            |     */
/*+-------------------------------------------------------------------+     */
/*| ANY      | CLOSE    | HEAD     | HEAD,END | READ NEXT DATA        |     */
/*|          | NEWSET   |          |OR UNKNWN |                       |     */
/*|          | INTARY OR|          |          |                       |     */
/*|          | FLTARY   |          |          |                       |     */
/*+-------------------------------------------------------------------+     */
/*| INQUIR   | CLOSE    | HEAD     | DATA,HEAD| INQUIRE NEXT DATA     |     */
/*|          | NEWSET   |          |END OR    |                       |     */
/*|          | INTARY OR|          |UNKNWN    |                       |     */
/*|          | FLTARY   |          |          |                       |     */
/*+-------------------------------------------------------------------+     */
/*| NEWSET   | NEWSET OR| HEAD     | HEAD,END | SEARCH FOR A NEXT NEW |     */
/*|          | CLOSE    |          |OR UNKNWN | SET, READ COMMENTS    |     */
/*+-------------------------------------------------------------------+     */
/*| INTARY   | UNCHANGED| HEAD     | HEAD,END | SEARCH FOR ARRAY DATA |     */
/*| FLTARY   |OR CLOSE  |          |OR UNKNWN | WITH SPECIFIED KEYWORD|     */
/*| DBLARY   |--------------------------------------------------------+     */
/*|          | UNCHANGED| DATA     | HEAD     | READ THE ARRAY DATA   |     */
/*|          |          |          |          | IMMEDIATELY BELOW     |     */
/*+===================================================================+     */

#include <stdio.h>
#include <string.h>
#include "gfc.h"

void gfrd_();

#if F2C001
#define open_f_   OPEN_F
#define close_f_  CLOSE_F
#define rw_c1_08_ RW_C1_08
#define rw_c1_30_ RW_C1_30
#define rw_c1_60_ RW_C1_60
#define rw_i1_    RW_I1
#define rw_i2_    RW_I2
#define rw_i3_    RW_I3
#define rw_in_    RW_IN
#define rw_f1_    RW_F1
#define rw_f2_    RW_F2
#define rw_f3_    RW_F3
#define rw_fn_    RW_FN
#define rw_d1_    RW_D1
#define rw_d2_    RW_D2
#define rw_d3_    RW_D3
#define rw_dn_    RW_DN
#endif

void gfrd_(int *ITYPE,int *IDEST,int *IMODE,int *IFPOS,int *ICAST,
           int *IUNIT,char *FILENMd,char *KEYWRD,char *NAME,char *COMMNT,
           int   *IVAR ,int   *IARY1 ,int   *IARY2,int   *IARY3,int   *IARY2D,
           float *FVAR ,float *FARY1 ,float *FARY2,float *FARY3,float *FARY2D,
           int *NUM,int *NUM2,int *maxd,int *max2d,int *IUT0,int *IERR)
{

  const char *ERMSGB = " ## SUBROUTINE GFRD  : FATAL      ERROR OCCURENCE; RETURNED";
  const char *EREXP1 = " AN ILLEGAL VALUE WAS SPECIFIED FOR A CONTROL PARAMETER    ";
  const char *EREXP2 = " THIS FILE IS NOT A GENERAL FILE VER. 2. FILE WAS CLOSED.  ";
  const char *EREXP3 = " INSUFFICIENT ARRAY SIZE. DATA HAS NOT BEEN READ IN.       ";
  const char *EREXP4 = " ALREADY AT THE END OF FILE. CANNOT CONTINUE READING.      ";
  const char *EREXP5 = " UNKNOWN FILE POSITION.      CANNOT CONTINUE READING.      ";
  const char *EREXP6 = " SPECIFIED DATA NOT FOUND TO THE END OF THE FILE.          ";
  const char *EREXP7 = " UNKNOWN HEADDER READ IN.    CANNOT CONTINUE READING.      ";
  const char *EREXP8 = " CANNOT OPEN THE SPECIFIED FILE.                           ";
  const char *EREXPA = " FILE POSITION DOES NOT MATCH TO SPECIFIED READING ACTION. ";
  const char *EREXPC = " FILE UNINT NUMBER IS ILLEAGAL                             ";

  const char *WRMSG0 = " ## SUBROUTINE GFRD  : WARNING                             ";
  const char *WRMSG1 = " ENDIAN WILl BE REVERSED WHEN READING IN BINARY MODE       ";
  const char *WRMSG2 = " THE DATA TYPE IS DIFFERENT BETWEEN DATA AND ARGUMENT      ";
  const char *WRMSG3 = " THE DATA WILL BE CASTED TO THE TYPE OF THE ARGUMENT       ";

  int   IRW0=0;
  int   IRW=1;

  int   MAX,MAX2;
  int   IDUM, ITYPEF, JRETRN;
  char  HEADRD[9] =" ";
  char  KEYF[9]   =" "; 
  char  NAMERD[31]=" ";
  char  FILENM[61]=" ";
  int   i,j,j1,j2,istore,iendian;
  int   ICAST0=0;
  float  FDUM,*FDUMR1,*FDUMR2,*FDUMR3;
  double DDUM,*DDUMR1,*DDUMR2,*DDUMR3;
  MAX =*maxd;
  MAX2=*max2d;  


  *IERR   = 0;
  JRETRN  = 0;
  strcpy(HEADRD,"NOTUNFMT");

  IDUM=IDUM;
  FDUM=FDUM;

  i=0;while( *(FILENMd+i) != ' ' && i<60){*(FILENM+i)=*(FILENMd+i);i++;}

  if(  *IUNIT >= FILEMAX ){ 
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXPC);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  /* [1] ENTRY FOR ITYPE = OPEN */
  if(*ITYPE == OPEN ){
      if( ( fp[*IUNIT] = fopen(FILENM,"r") ) == NULL){
          sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
          sprintf(MESSAGE,EREXP8);printstd(IUT0,MESSAGE);
          *IERR = 1;
          return;
  }

  fread ((char *)HEADRD,8,1,fp[*IUNIT] );
      if( strcmp_n(7,HEADRD,HEADER[BINARY]) ==0 ){
          *IMODE = BINARY ;
          fread ((char *)&iendian,4,1,fp[*IUNIT] );
          if(iendian==IENDIAN[0]){
              IREV[*IUNIT]=0;
          }
          else if(iendian==IENDIAN[1]){
              IREV[*IUNIT]=1;
              sprintf(MESSAGE,WRMSG0);printstd(IUT0,MESSAGE);
              sprintf(MESSAGE,WRMSG1);printstd(IUT0,MESSAGE);
          }
          else{
              sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
              sprintf(MESSAGE,EREXP2);printstd(IUT0,MESSAGE);
              fclose(fp[*IUNIT]);
              *IERR = 1 ;
              return ;
          }
      }
      else{
          fclose(fp[*IUNIT]);
          *IMODE=UNFMT;
          open_f_(IUNIT,IMODE,FILENMd,IERR);
          rw_c1_08_(&IRW,IMODE,IUNIT,HEADRD,IERR);
          if( strcmp_n(7,HEADRD,HEADER[UNFMT]) ==0 )*IMODE = UNFMT ;
          else{
              close_f_(IUNIT,IERR);
              *IMODE=ASCII;
              open_f_(IUNIT,IMODE,FILENMd,IERR);
              rw_c1_08_(&IRW,IMODE,IUNIT,HEADRD,IERR);
               if ( strcmp_n(7,HEADRD,HEADER[ASCII]) ==0 )*IMODE = ASCII ;
              else{
                  sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
                  sprintf(MESSAGE,EREXP2);printstd(IUT0,MESSAGE);
                  close_f_(IUNIT,IERR);
                  *IERR = 1 ;
                  return ;
              }/*  if ( strcmp(HEADRD,HEADER[ASCII]) ==0 ) */
          }/* if ( strcmp(HEADRD,HEADER[UNFMT]) ==0 )  */ 
      } /* if( ( fp[*IUNIT] = fopen(FILENM,"r") ) == NULL) */
    

      if( *IMODE == BINARY){
           fread ((char *)NUM2,4,1,fp[*IUNIT] );
           if(IREV[*IUNIT]==1) conv(NUM2,4);
      }  
      else{
           rw_i1_(&IRW,IMODE,IUNIT,&ONE,&ONE,NUM2,IERR); 
           if(*IERR!=0)return;
      }/* if( *IMODE == BINARY) */
      for(i=0;i<*NUM2;i++){

          istore=i;if(istore>MAX2)istore=MAX2;
          for(j=0;j<60;j++)*(COMMNT+istore*60+j)= ' ';

          if( *IMODE == BINARY){
              fread ((char *)COMMNT+istore*60,60,1,fp[*IUNIT] );
          }
          else{
              rw_c1_60_(&IRW,IMODE,IUNIT,COMMNT+istore*60,IERR);
          }/* if( *IMODE == BINARY) */

      }/* for(i=0;i<*NUM2;i++) */

      *IFPOS = HEAD;

      if(*NUM2 >= MAX2 ){
      sprintf(MESSAGE,"%d %d ",*NUM2,MAX2);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,ERMSGB);             printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP3);             printstd(IUT0,MESSAGE);
      *IERR = 1;}

      return;

  } /* if(*ITYPE == OPEN ) */

  /* [2] ENTRY FOR ITYPE = CLOSE */
  if(*ITYPE == CLOSE ){
      if(*IMODE==BINARY){
          fclose(fp[*IUNIT]);
      }
      else{
          close_f_(IUNIT,IERR);
      }
      return;
  } /* if(*ITYPE == CLOSE ) */

  /* [3] CHECK CONTROL PARAMETERS PASSED */

  if(*IFPOS == END){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP4);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if( (*IFPOS != HEAD) && (*IFPOS != DATA) ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP5);printstd(IUT0,MESSAGE);
      *IFPOS= UNKNWN;
      *IERR = 1;
      return;
  }

  if(    (*IFPOS == DATA) 
      && ( (*ITYPE == ANY) || (*ITYPE == INQUIR)  || (*ITYPE == NEWSET) ) ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXPA);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if( (*IMODE != ASCII) && (*IMODE != UNFMT) && (*IMODE != BINARY) ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP1);printstd(IUT0,MESSAGE);
      *IFPOS= UNKNWN;
      *IERR = 1;
      return;
  }

  if(   ( (*ITYPE==INTARY) || (*ITYPE==FLTARY) || (*ITYPE==DBLARY) || (*ITYPE==ANY  ) )
     && ( (*IDEST!=ARY2D ) && (*IDEST!=ARY1D ) && (*IDEST!=ARY0D) ) ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP1);printstd(IUT0,MESSAGE);
      *IFPOS= UNKNWN;
      *IERR = 1;
      return;
  }

  if( (*ITYPE==-INTARY) || (*ITYPE==-FLTARY) || (*ITYPE==-DBLARY) ){  
      JRETRN=-1;    
      *ITYPE=-*ITYPE;
  }

  if(    (*ITYPE!=ANY   ) && (*ITYPE!=INQUIR ) && (*ITYPE!=NEWSET) 
      && (*ITYPE!=INTARY) && (*ITYPE!=FLTARY ) && (*ITYPE!=DBLARY ) ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP1);printstd(IUT0,MESSAGE);
      *IFPOS= UNKNWN;
      *IERR = 1;
      return;
  }

  if( *IFPOS==DATA && JRETRN== 0 ) JRETRN=2;
  if( *IFPOS==DATA && JRETRN==-1 ) JRETRN=1;

  line1000: ;

  /* [4] READ NEXT HEADDER */
  if(  *IFPOS == HEAD){


     /* READ DATA HEADDER */
      if(*IMODE == BINARY){
          fread ((char *)HEADRD,8,1,fp[*IUNIT] );
      }
      else{
          rw_c1_08_(&IRW,IMODE,IUNIT,HEADRD,IERR);
      }/* (*IMODE == BINARY) */

       /* CHECK READ IN HEADDER */
      if(    (strcmp(HEADRD,HEADER[INTARY]) !=0 )
          && (strcmp(HEADRD,HEADER[FLTARY]) !=0 )
          && (strcmp(HEADRD,HEADER[DBLARY]) !=0 )
          && (strcmp(HEADRD,HEADER[NEWSET]) !=0 )
          && (strcmp(HEADRD,HEADER[CLOSE ]) !=0 )){
          sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
          sprintf(MESSAGE,EREXP7);printstd(IUT0,MESSAGE);
          *IFPOS= UNKNWN;
          *IERR = 1;
          return;
      }/* if( ...)  */


      /* CLOSE HEADDER */
      if( strcmp(HEADRD,HEADER[CLOSE ]) == 0){
          if(    *ITYPE != INQUIR && *ITYPE != ANY  ){
              sprintf(MESSAGE,EREXP6);printstd(IUT0,MESSAGE);
              *IERR = 1;
          } 
          *ITYPE= CLOSE;
          *IFPOS= END;
          return;
      }/* if( strcmp(HEADRD,HEADER[CLOSE ]) == 0) */

      /* NEWSET HEADDER */
      if( strcmp(HEADRD,HEADER[NEWSET] ) == 0 ){
           if( *IMODE == BINARY){
               fread ((char *)NUM2,4,1,fp[*IUNIT] );
               if(IREV[*IUNIT]==1) conv(NUM2,4);
      }
      else{
          rw_i1_(&IRW,IMODE,IUNIT,&ONE,&ONE,NUM2,IERR);
          if(*IERR!=0)return;
      }/* *IMODE == BINARY*/

          for(i=0;i<*NUM2;i++){

              istore=i;if(istore>MAX2)istore=MAX2;
              for(j=0;j<60;j++)*(COMMNT+istore*60+j)= ' ';

              if( *IMODE == BINARY){
                  fread ((char *)COMMNT+istore*60,60,1,fp[*IUNIT] );
              }
              else{
                  rw_c1_60_(&IRW,IMODE,IUNIT,COMMNT+istore*60,IERR);  
              }/* if( *IMODE == BINARY) */
  
          }/* for(i=0;i<*NUM2;i++) */ 

          /* Set IFPOS to HEAD */
          *IFPOS = HEAD;

          if( (*ITYPE==ANY) || (*ITYPE==INQUIR) || (*ITYPE==NEWSET) ){
              if(*NUM2 > MAX2){
                  sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
                  sprintf(MESSAGE,EREXP3);printstd(IUT0,MESSAGE);
                  *IERR = 1;
              } 
              *ITYPE= NEWSET;
              return;
          }/* (*ITYPE==ANY) || (*ITYPE==INQUIR) || (*ITYPE==NEWSET) */
      }/* if( strcmp(HEADRD,HEADER[NEWSET] ) == 0 ) */

      /* INTARY AND FLTARY HEADDER */
      if(   (strcmp(HEADRD,HEADER[INTARY])==0)  
        || (strcmp(HEADRD,HEADER[FLTARY])==0)  
        || (strcmp(HEADRD,HEADER[DBLARY])==0)  ) {

          if( strcmp(HEADRD,HEADER[INTARY])==0 )ITYPEF = INTARY ;
          if( strcmp(HEADRD,HEADER[FLTARY])==0 )ITYPEF = FLTARY ;
          if( strcmp(HEADRD,HEADER[DBLARY])==0 )ITYPEF = DBLARY ;

          for(j=0;j<30;j++)*(NAMERD+j)= ' ';

              if(*IMODE == BINARY){
                  fread ((char *)KEYF  , 8,1,fp[*IUNIT] );
                  fread ((char *)NAMERD,30,1,fp[*IUNIT] );
              }    
              else{
                  rw_c1_08_(&IRW,IMODE,IUNIT,KEYF,IERR);
                  rw_c1_30_(&IRW,IMODE,IUNIT,NAMERD,IERR);
              }/* if(*IMODE == BINARY) */ 

          /* Set IFPOS to DATA */
          *IFPOS = DATA;
         

          if( *ITYPE == INQUIR){
              *ITYPE = ITYPEF;
              strcpy(KEYWRD , KEYF  );
              for(j=0;j<30;j++)*(NAME+j)=*(NAMERD+j) ; 
              return;
          }

          if( ( *ITYPE == ANY) || 
              ( (ITYPEF == *ITYPE)  && (strcmp(KEYF,KEYWRD)==0 ) ) 
              ){
              if(JRETRN==-1){JRETRN = 1;}
              else          {JRETRN = 2;}
              *ITYPE = ITYPEF;
              strcpy(KEYWRD , KEYF  );
              for(j=0;j<30;j++)*(NAME+j)=*(NAMERD+j);
          }           

      }/* if(   (strcmp(HEADRD,HEADER[INTARY])==0)  ,,, )  */
  } /* if(  *IFPOS == HEAD) */

  /* [5] READ NEXT DATA */
  if( *IFPOS == DATA){

      /* Read Array Size */
      if( *IMODE == BINARY ){
          fread ((char *)NUM2,4,1,fp[*IUNIT] );
          fread ((char *)NUM ,4,1,fp[*IUNIT] );
          if(IREV[*IUNIT]==1){
              conv(NUM2,4);
              conv(NUM ,4);
          }
      }
      else{
          rw_i2_(&IRW,IMODE,IUNIT,&ONE,&ONE,NUM2,NUM,IERR);
          if(*IERR!=0)return;
      }/* if( *IMODE == BINARY) */

      if(*NUM2<=0 || *NUM<=0 )goto line150;

      /* CHECK TYPE OF THE DATA AND THE ARGUMENT */
      if(    ITYPE[0]==FLTARY && ITYPE[1]==DBLARY 
          || ITYPE[0]==DBLARY && ITYPE[1]==FLTARY  ){
          sprintf(MESSAGE,WRMSG0);printstd(IUT0,MESSAGE);
          sprintf(MESSAGE,WRMSG2);printstd(IUT0,MESSAGE);
          sprintf(MESSAGE," TYPE OF THE FILE DATA IS %s",HEADER[ITYPE[0]]);
          printstd(IUT0,MESSAGE);
          sprintf(MESSAGE," TYPE OF THE AURGEMENT IS %s",HEADER[ITYPE[1]]);
          printstd(IUT0,MESSAGE);
          if(*ICAST==0){
              sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
              *IERR=1;
              return;
          }
          else{
              sprintf(MESSAGE,WRMSG0);printstd(IUT0,MESSAGE);
              sprintf(MESSAGE,WRMSG3);printstd(IUT0,MESSAGE);
              ICAST0=1;
          }
      }

      if(    JRETRN <= 1 
          || ( *NUM  > MAX  )
          || ( *NUM2 > MAX2 ) 
          || ( *IDEST == ARY1D && *NUM2>3 )
          || ( *IDEST == ARY0D && (*NUM2>1 || *NUM>1 ) ) 
      ){
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&IDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&IDUM,4); }
              }
          }/* if(*ITYPE) */
          if( *ITYPE==INTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_in_(&IRW0,IMODE,IUNIT,&maxd,&max2d,NUM,NUM2,&IDUM,IERR);
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&FDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&FDUM,4); }
              }
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_fn_(&IRW0,IMODE,IUNIT,&maxd,&max2d,NUM,NUM2,&FDUM,IERR);
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&DDUM,8,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&DDUM,8); }
              }
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_dn_(&IRW0,IMODE,IUNIT,&maxd,&max2d,NUM,NUM2,&DDUM,IERR);
          }/* if(*ITYPE) */
	  goto line150;
      } /* if(    JRETRN <= 1  ... ) */

      if( *IDEST==ARY0D ){ 
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              fread ((char *)IVAR,4,1,fp[*IUNIT] );   
              if(IREV[*IUNIT]==1) { conv(IVAR,4); }
          }/* if(*ITYPE) */
          if( *ITYPE==INTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_i1_(&IRW,IMODE,IUNIT,&ONE,&ONE,IVAR,IERR);
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              fread ((char *)&FDUM,4,1,fp[*IUNIT] );   
              if(IREV[*IUNIT]==1) { conv(FVAR,4); }
              if(ICAST0==0){
                  *FVAR=FDUM;
              }
              else{
                  DDUM=(double)FDUM;
                  *(double *)FVAR=DDUM;
              }
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_f1_(&IRW,IMODE,IUNIT,&ONE,&ONE,&FDUM,IERR);
              if(ICAST0==0){
                  *FVAR=FDUM;    
              }
             else{
                  *FVAR=(double)FDUM;    
              } 
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              fread ((char *)&DDUM,8,1,fp[*IUNIT] );   
              if(IREV[*IUNIT]==1) { conv(FVAR,8); }
              if(ICAST0==0){
              *(double *)FVAR=DDUM;
              }
              else{
                  FDUM=(float)DDUM;
                  *FVAR=FDUM;
              } 
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_d1_(&IRW,IMODE,IUNIT,&ONE,&ONE,&DDUM,IERR);
              if(ICAST0==0){
                  *(double *)FVAR=DDUM;
              }
              else{
                  *FVAR=(float)DDUM;
              } 
          } /* if(*ITYPE) */

          if(*IERR!=0)return;
	  goto line150;
      } /* if( *IDEST==ARY0D ) */
 
      if( *NUM2==1 && *IDEST==ARY1D ){ 
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM;i++){
                  fread ((char *)&IDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&IDUM,4); }
                  *(IARY1+i)=IDUM;
              }
          }/* if(*ITYPE) */
          if( *ITYPE==INTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_i1_(&IRW,IMODE,IUNIT,maxd,NUM,IARY1,IERR);
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM;i++){
                  fread ((char *)&FDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&FDUM,4); }
                  if(ICAST0==0){
                      *(FARY1+i)=FDUM;
                  }
                  else{
                      DDUM=(double)FDUM;    
                      *(double *)(FARY1+i*2)=DDUM;
                  }
              }
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              if(ICAST0==0){
                  rw_f1_(&IRW,IMODE,IUNIT,maxd,NUM,FARY1,IERR);
              }
              else{
                  FDUMR1 = (float  *)malloc((*NUM)*4);
                  rw_f1_(&IRW,IMODE,IUNIT,NUM,NUM,FDUMR1,IERR);
                  for(i=0;i<*NUM;i++){
                      *(double *)(FARY1+i*2)=(double)FDUMR1[i];
                  }
                  free(FDUMR1);
              } 
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i< *NUM;i++){
                  fread ((char *)&DDUM,8,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&DDUM,8); }
                  if(ICAST0==0){
                      *(double *)(FARY1+i*2)=DDUM;
                  }
                  else{
                      FDUM=(float)DDUM;
                      *(FARY1+i)=FDUM;
                  }
              }
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              if(ICAST0==0){
                  rw_d1_(&IRW,IMODE,IUNIT,maxd,NUM,FARY1,IERR);
              }
              else{
                  DDUMR1 = (double *)malloc((*NUM)*8);
                  rw_d1_(&IRW,IMODE,IUNIT,NUM,NUM,DDUMR1,IERR);
                  for(i=0;i<*NUM;i++){
                      *(FARY1+i)=(float)DDUMR1[i];
                  }
                  free(DDUMR1);
              } 
 	  }/* if(*ITYPE) */

          if(*IERR!=0)return;
          goto line150;
      } /* if( *NUM2==1 && *IDEST==ARY1D ) */
 
      if( *NUM2==2 && *IDEST==ARY1D ){ 
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&IDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&IDUM,4); }
                  if((i%2)==0)*(IARY1+(i-0)/2)=IDUM;
                  if((i%2)==1)*(IARY2+(i-1)/2)=IDUM;
              }
          }/* if(*ITYPE) */
          if( *ITYPE==INTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              rw_i2_(&IRW,IMODE,IUNIT,maxd,NUM,IARY1,IARY2,IERR);
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && *IMODE==BINARY  ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&FDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&FDUM,4); }
                   if(ICAST0==0){
                      if((i%2)==0)*(FARY1+(i-0)/2)=FDUM;
                      if((i%2)==1)*(FARY2+(i-1)/2)=FDUM;
                  }  
                  else{
                      DDUM=(double)FDUM;
                      if((i%2)==0) *(double *)(FARY1+((i-0)/2)*2)=DDUM;
                      if((i%2)==1) *(double *)(FARY2+((i-1)/2)*2)=DDUM;
                  }
              }
	  }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              if(ICAST0==0){
                  rw_f2_(&IRW,IMODE,IUNIT,maxd,NUM,FARY1,FARY2,IERR);
              }
              else{
                  FDUMR1 = (float  *)malloc((*NUM)*4);
                  FDUMR2 = (float  *)malloc((*NUM)*4);
                  rw_f2_(&IRW,IMODE,IUNIT,NUM,NUM,FDUMR1,FDUMR2,IERR);
                  for(i=0;i<*NUM;i++){
                      *(double *)(FARY1+i*2)=(double)FDUMR1[i];
                      *(double *)(FARY2+i*2)=(double)FDUMR2[i];
                  }
                  free(FDUMR1);
                  free(FDUMR2);
              } 
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && *IMODE==BINARY  ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&DDUM,8,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&DDUM,8); }
                  if(ICAST0==0){
                      if((i%2)==0) *(double *)(FARY1+((i-0)/2)*2)=DDUM;
                      if((i%2)==1) *(double *)(FARY2+((i-1)/2)*2)=DDUM;
                  }
                  else{
                      FDUM=(float)DDUM;
                      if((i%2)==0)*(FARY1+(i-0)/2)=FDUM;
                      if((i%2)==1)*(FARY2+(i-1)/2)=FDUM;
                  }
              }
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              if(ICAST0==0){
                  rw_d2_(&IRW,IMODE,IUNIT,maxd,NUM,FARY1,FARY2,IERR);
              }
              else{
                  DDUMR1 = (double *)malloc((*NUM)*8);
                  DDUMR2 = (double *)malloc((*NUM)*8);
                  rw_d2_(&IRW,IMODE,IUNIT,NUM,NUM,DDUMR1,DDUMR2,IERR);
                  for(i=0;i<*NUM;i++){
                      *(FARY1+i)=(float)DDUMR1[i];
                      *(FARY2+i)=(float)DDUMR2[i];
                  }
                  free(DDUMR1);
                  free(DDUMR2);
              } 
	  }/* if(*ITYPE) */

          if(*IERR!=0)return;
          goto line150;
      } /* if( *NUM2==2 && *IDEST==ARY1D ) */
 
      if( *NUM2==3 && *IDEST==ARY1D ){ 
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&IDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&IDUM,4); }
                  if((i%3)==0)*(IARY1+(i-0)/3)=IDUM;
                  if((i%3)==1)*(IARY2+(i-1)/3)=IDUM;
                  if((i%3)==2)*(IARY3+(i-1)/3)=IDUM;
              }
          }/* if(*ITYPE) */
          if( *ITYPE==INTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
               rw_i3_(&IRW,IMODE,IUNIT,maxd,NUM,IARY1,IARY2,IARY3,IERR);
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&FDUM,4,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&FDUM,4); }
                  if(ICAST0==0){
                      if((i%3)==0)*(FARY1+(i-0)/3)=FDUM;
                      if((i%3)==1)*(FARY2+(i-1)/3)=FDUM;
                      if((i%3)==2)*(FARY3+(i-1)/3)=FDUM;
                  }
                  else{
                      DDUM=(double)FDUM;
                      if((i%3)==0) *(double *)(FARY1+((i-0)/3)*2)=DDUM;
                      if((i%3)==1) *(double *)(FARY2+((i-1)/3)*2)=DDUM;
                      if((i%3)==2) *(double *)(FARY3+((i-1)/3)*2)=DDUM;
                  }
              }
	  }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              if(ICAST0==0){
                  rw_f3_(&IRW,IMODE,IUNIT,maxd,NUM,FARY1,FARY2,FARY3,IERR);
              }
              else{
                  FDUMR1 = (float  *)malloc((*NUM)*4);
                  FDUMR2 = (float  *)malloc((*NUM)*4);
                  FDUMR3 = (float  *)malloc((*NUM)*4);
                  rw_f3_(&IRW,IMODE,IUNIT,NUM,NUM,FDUMR1,FDUMR2,FDUMR3,IERR);
                  for(i=0;i<*NUM;i++){
                      *(double *)(FARY1+i*2)=(double)FDUMR1[i];
                      *(double *)(FARY2+i*2)=(double)FDUMR2[i];
                      *(double *)(FARY3+i*2)=(double)FDUMR3[i];
                  }
                  free(FDUMR1);
                  free(FDUMR2);
                  free(FDUMR3);
              } 
          } /* if(*ITYPE) */  
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                  fread ((char *)&DDUM,8,1,fp[*IUNIT] ); 
                  if(IREV[*IUNIT]==1) { conv(&DDUM,8); }
                  if(ICAST0==0){
                      if((i%3)==0) *(double *)(FARY1+((i-0)/3)*2)=DDUM;
                      if((i%3)==1) *(double *)(FARY2+((i-1)/3)*2)=DDUM;
                      if((i%3)==2) *(double *)(FARY3+((i-1)/3)*2)=DDUM;
                  }
                  else{
                      FDUM=(float)DDUM;
                      if((i%3)==0)*(FARY1+(i-0)/3)=FDUM;
                      if((i%3)==1)*(FARY2+(i-1)/3)=FDUM;
                      if((i%3)==2)*(FARY3+(i-1)/3)=FDUM;
                  }
              }
	  }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              if(ICAST0==0){
                  rw_d3_(&IRW,IMODE,IUNIT,maxd,NUM,FARY1,FARY2,FARY3,IERR);
              }
              else{
                  DDUMR1 = (double *)malloc((*NUM)*8);
                  DDUMR2 = (double *)malloc((*NUM)*8);
                  DDUMR3 = (double *)malloc((*NUM)*8);
                  rw_d3_(&IRW,IMODE,IUNIT,NUM,NUM,DDUMR1,DDUMR2,DDUMR3,IERR);
                  for(i=0;i<*NUM;i++){
                      *(FARY1+i)=(float)DDUMR1[i];
                      *(FARY2+i)=(float)DDUMR2[i];
                      *(FARY3+i)=(float)DDUMR3[i];
                  }
                  free(DDUMR1);
                  free(DDUMR2);
                  free(DDUMR3);
              } 
          } /* if(*ITYPE) */  

          if(*IERR!=0)return;
          goto line150;
      } /* if( *NUM2==3 && *IDEST==ARY1D ) */

      j1=0;
      j2=0;                               
      if( *IDEST==ARY2D ){ 
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                   fread ((char *)&IDUM,4,1,fp[*IUNIT] ); 
                   if(IREV[*IUNIT]==1) { conv(&IDUM,4); }
                   *(IARY2D+j1*MAX2+j2)=IDUM;                
                   j2++;if(j2==*NUM2){j2=0;j1++;}
              }
	  }/* if(*ITYPE) */
          if( *ITYPE==INTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
               rw_in_(&IRW,IMODE,IUNIT,maxd,max2d,NUM,NUM2,IARY2D,IERR);
	  }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                   fread ((char *)&FDUM,4,1,fp[*IUNIT] ); 
                   if(IREV[*IUNIT]==1) { conv(&FDUM,4); }
                   if(ICAST0==0){
                       *(FARY2D+j1*MAX2+j2)=FDUM;                
                   }
                   else{
                       DDUM=(double)FDUM;
                       *(double *)(FARY2D+(j1*MAX2+j2)*2)=DDUM;                
                   }
                   j2++;if(j2==*NUM2){j2=0;j1++;}
              }
          }/* if(*ITYPE) */
          if( *ITYPE==FLTARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
               if(ICAST0==0){
                   rw_fn_(&IRW,IMODE,IUNIT,maxd,max2d,NUM,NUM2,FARY2D,IERR);
               }
               else{
                   FDUMR1 = (float  *)malloc((*NUM)*(*NUM2)*4);
                   rw_fn_(&IRW,IMODE,IUNIT,NUM,NUM2,NUM,NUM2,FDUMR1,IERR);
                   for(i=0;i<*NUM;i++){
                       for(j=0;j<*NUM2;j++){
                           *(double *)(FARY2D+(j+i*(*max2d))*2)=(double)FDUMR1[j+i*(*NUM2)];
                       }
                   }
                  free(FDUMR1);
               } 
	  }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i<(*NUM)*(*NUM2);i++){
                   fread ((char *)&DDUM,8,1,fp[*IUNIT] ); 
                   if(IREV[*IUNIT]==1) { conv(&DDUM,8); }
                   if(ICAST0==0){
                       *(double *)(FARY2D+(j1*MAX2+j2)*2)=DDUM;                
                    }
                    else{
                       FDUM=(float)DDUM;
                       *(FARY2D+j1*MAX2+j2)=FDUM;                
                    } 
                   j2++;if(j2==*NUM2){j2=0;j1++;}
              }
          }/* if(*ITYPE) */
          if( *ITYPE==DBLARY && (*IMODE==ASCII || *IMODE==UNFMT) ){
              if(ICAST0==0){
                  rw_dn_(&IRW,IMODE,IUNIT,maxd,max2d,NUM,NUM2,FARY2D,IERR);
              }
              else{
                  DDUMR1 = (double *)malloc((*NUM)*(*NUM2)*8);
                  rw_dn_(&IRW,IMODE,IUNIT,NUM,NUM2,NUM,NUM2,DDUMR1,IERR);
                  for(i=0;i<*NUM;i++){
                      for(j=0;j<*NUM2;j++){
                         *(FARY2D+j+i*(*max2d))=(float)DDUMR1[j+i*(*NUM2)];
                      }
                  }
                  free(DDUMR1);
              }/* if(ICAST) */  
          }/* if(*ITYPE) */

          if(*IERR!=0)return;
          goto line150;
      } /* if( *IDEST==ARY2D ) */

      line150: ;

      *IFPOS=HEAD;
  
      if(JRETRN>=1){
          if( (JRETRN==2) &&
           (    ( *NUM  > MAX  )
             || ( *NUM2 > MAX2 )
             || ( *IDEST == ARY1D && *NUM2>3 )
             || ( *IDEST == ARY0D && (*NUM2>1 || *NUM>1 ) ) 
           )
          ){
              sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
              sprintf(MESSAGE,EREXP3);printstd(IUT0,MESSAGE);
              *IERR = 1;} 
          return;
      }/* JRETRN>=1 */

  } /* if(  *IFPOS == DATA) */
  goto line1000;
}
