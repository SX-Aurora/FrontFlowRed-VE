/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  */
/*    BASIC OUTPUT UTILITY FOR GENERAL FILE VERSION 2.1                     */
/*                                                                          */
/*      Written by  Y. Yamade                                               */
/*      Since       9 Oct. 2002                                             */
/*      Reference   gfwrt.f of LES3D                                        */
/*                                                                          */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  */
/*                                                                          */
/*    ARGUMENT LISTINGS                                                     */
/*      (1) INPUT                                                           */
/*INT*4    ITYPE       ; SPECIFIES THE TYPE OF ACTION TO BE TAKEN           */
/*           OPEN     --- OPEN  FILE,     WRITE COMMENT STRINGS             */
/*           NEWSET   --- DEFINE NEW SET, WRITE COMMENT STRINGS             */
/*           CLOSE    --- CLOSE FILE                                        */
/*           INTARY   --- WRITE KEYWORD, DATA NAME, INTEGER ARRAY DATA      */
/*           FLTARY   --- WRITE KEYWORD, DATA NAME, FLOAT   ARRAY DATA      */
/*    NOTES; IF EITHER NUM2 OR NUM ARGUMENTS (DESCRIBED LATER) IS LESS      */
/*          THAN 1, NO DATA WILL BE WRITTEN TO THE FILE EXCEPT THE DATA     */
/*          HEADDER, DATA NAME, AND THE ELEMENT NUMBERS.                    */
/*                                                                          */
/*INT*4    ISRCE       ; SPECIFIES SOURCE VARIABLE/ARRAYS TO BE WRITTEN     */
/*                      TO THE FILE WHEN WRITING ARRAY DATA                 */
/*           ARY2D    --- 2-D ARRAY PASSED THROUGH IARY2D OR FARY2D         */
/*                       ARGUMENTS WILL BE WRITTEN TO THE FILE.             */
/*           ARY1D    --- 1-D ARRAYS PASSED THROUGH IARY1, IARY2, IARY3     */
/*                       OR FARY1, FARY2, FARY3 ARGUMENTS WILL BE           */
/*                       WRITTEN TO THE FILE.                               */
/*           ARY0D    --- VARIABLE PASSED THROUGH IVAR OR FVAR ARGUMENTS    */
/*                       WILL BE WRITTEN TO THE FILE.                       */
/*    NOTES; FOLLOWING CASES LISTED BELOW WILL BE REGARDED AS AN ERROR.     */
/*          IN SUCH CASES, NO DATA WILL WRITTEN TO THE FILE, EXCEPT THE     */
/*          DATA HEADDER, DATA NAME, AND THE ELEMENT NUMBERS.               */
/*           (1) NUM2 > MAX2              FOR ISRCE = ARY2D                 */
/*           (2) NUM2 > 3                 FOR ISRCE = ARY1D                 */
/*           (3) NUM2 > 1    OR NUM > 1   FOR ISRCE = ARY0D                 */
/*                                                                          */
/*INT*4    IMODE       ; SPECIFIES THE TYPE OF OUTPUT MODE                  */
/*           ASCII    --- ASCII               MODE                          */
/*           UNFMT    --- FORTRAN UNFORMATTED MODE                          */
/*           BINARY   --- BINARY              MODE ( NOT SUPPORTED )        */
/*INT*4    IUNIT       ; FILE BUFFER NUMBER TO BE USED TO ACCESS            */
/*                      THE SPECIFIED FILE                                  */
/*CHAR*(*) FILENM      ; FILE NAME TO BE OPENED                             */
/*    NOTES; FILENM ARGUMENT MUST BE SPECIFIED ONLY WHEN OPENING            */
/*          A FILE ( IE. ITYPE = OPEN ).                                    */
/*CHAR*8   KEYWRD      ; A CHARACTER STRING TO BE PUT ON EACH DATA,         */
/*                      SPECIFYING ITS KEY CHARACTERISTICS                  */
/*CHAR*30  NAME        ; A CHARACTER STRING TO BE PUT ON EACH DATA,         */
/*                      SPECIFYING ITS NAME.                                */
/*                                                                          */
/*CHAR*60  COMMNT(NUM2); AN 1-D ARRAY OF CHARACTER STRINGS TO BE PUT ON     */
/*                      EACH FILE AND EACH SET, WHEN NUM2 ARGUMENT IS       */
/*                      GREATER THAN 0.                                     */
/*    NOTES; KEYWRD ARGUMENT AND NAME ARGUMENT MUST BE SPECIFIED            */
/*          EXCEPT WHEN OPENING FILE, DEFINING A NEW SET, OR CLOSING        */
/*          FILE.                                                           */
/*           COMMNT ARGUMENT MUST BE SPECIFIED ONLY WHEN OPENING FILE       */
/*          OR DEFINING A NEW SET ( IE. ITYPE = OPEN AND NEWSET ).          */
/*INT*4    IVAR        ; AN INTEGER VARIABLE   TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=INTARY, ISRCE=ARY0D, AND NUM2=NUM=1      */
/*INT*4    IARY1(NUM)  ; AN 1-D INTEGER ARRAY  TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=INTARY, ISRCE=ARY1D, AND NUM2=(1,2,3)    */
/*INT*4    IARY2(NUM)  ; AN 1-D INTEGER ARRAY  TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=INTARY, ISRCE=ARY1D, AND NUM2=(  2,3)    */
/*INT*4    IARY3(NUM)  ; AN 1-D INTEGER ARRAY  TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=INTARY, ISRCE=ARY1D, AND NUM2=(    3)    */
/*INT*4    IARY2D(MAX2,; A  2-D INTEGER ARRAY  TO BE WRITTEN TO THE FILE    */
/*                 NUM) WHEN ITYPE=INTARY, ISRCE=ARY2D                      */
/*REAL*4   FVAR        ; AN FLOAT   VARIABLE   TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=FLTARY, ISRCE=ARY0D, AND NUM2=NUM=1      */
/*REAL*4   FARY1(NUM)  ; AN 1-D FLOAT   ARRAY  TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=FLTARY, ISRCE=ARY1D, AND NUM2=(1,2,3)    */
/*REAL*4   FARY2(NUM)  ; AN 1-D FLOAT   ARRAY  TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=FLTARY, ISRCE=ARY1D, AND NUM2=(  2,3)    */
/*REAL*4   FARY3(NUM)  ; AN 1-D FLOAT   ARRAY  TO BE WRITTEN TO THE FILE    */
/*                      WHEN ITYPE=FLTARY, ISRCE=ARY1D, AND NUM2=(    3)    */
/*REAL*4   FARY2D(MAX2,; A  2-D FLOAT   ARRAY  TO BE WRITTEN TO THE FILE    */
/*                 NUM) WHEN ITYPE=FLTARY, ISRCE=ARY2D                      */
/*INT*4    NUM         ; NUMBER OF FIRST  ARRAY ELEMENTS TO BE WRITTEN      */
/*                      TO THE FILE FOR WRITING ARRAY DATA                  */
/*INT*4    NUM2        ; NUMBER OF SECOND ARRAY ELEMENTS TO BE WRITTEN      */
/*                      TO THE FILE FOR WRITING ARRAY DATA,                 */
/*                       NUMBER OF COMMENT STRINGS       TO BE WRITTEN      */
/*                      TO THE FILE FOR OPENING FILE OR DEFINING A NEW      */
/*                      SET ( IE. ITYPE = OPEN AND NEWSET ).                */
/*INT*4    MAX2        ; THE SECOND DIMENSION SIZE OF 2-D ARRAYS PASSED     */
/*    NOTES; NUM, AND NUM2 ARGUMENTS MUST BE SPECIFIED WHEN WRITING AN      */
/*          ARRAY DATA (IE. ITYPE=INTARY OR FLTARY).                        */
/*           NUM2 ARGUMENT MUST ALSO BE SPECIFIED WHEN OPENING FILE OR      */
/*          DEFINING A NEW SET (IE. ITYPE = OPEN OR NEWSET).                */
/*           MAX2 ARGUMENT MUST BE SPECIFIED WHEN WRITING AN ARRAY DATA     */
/*          FROM 2-D SOURCE ARRAYS (IE. ITYPE=(INTARY OR FLTARY) AND        */
/*          (NUM2>3 OR ISRCE=ARY2D).                                        */
/*INT*4    IUT0        ; FILE   NUMBER TO REPORT ERROR OCCURENCE            */
/*      (2) OUPTUT                                                          */
/*INT*4    IERR        ; RETURN CODE WHOSE VALUE WILL BE EITHER             */
/*           0       --- INDICATING SUCCESSFUL TERMINATION                  */
/*       OR  1       --- INDICATING OCCURENCE OF SOME ERROR CONDITIONS      */
/****************************************************************************/
#include <stdio.h>
#include <string.h>
#include "gfc.h"

void gfwrt_();

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

void gfwrt_(int *ITYPE,int *ISRCE,int *IMODE,int *IREV,
           int *IUNIT,char *FILENMd,char *KEYWRD,char *NAME,char *COMMNT,
           int   *IVAR ,int   *IARY1 ,int   *IARY2,int   *IARY3,int   *IARY2D,
           float *FVAR ,float *FARY1 ,float *FARY2,float *FARY3,float *FARY2D,
           int *NUM,int *NUM2,int *max2d,int *IUT0,int *IERR)
{

  const char *ERMSGA = " ## SUBROUTINE GFWRT : RECOVERBLE ERROR OCCURENCE; CONTINUE";
  const char *ERMSGB = " ## SUBROUTINE GFWRT : FATAL      ERROR OCCURENCE; RETURNED";
  const char *ERRCV1 = " ## THE PARAMETER 'IREV'  WILL BE SET TO THE DEFAULT VALUE ";
  const char *EREXP1 = " AN ILLEGAL VALUE WAS SPECIFIED FOR A CONTROL PARAMETER    ";
  const char *EREXP2 = " SPECIFIED DATA SIZE (NUM2 AND/OR NUM) IS TO LARGE.        ";
  const char *EREXP3 = "  CANNOT OPEN THE SPECIFIED FILE.                          ";
  const char *EREXP4 = " SPECIFIED FILE NUMSER MUST BE LESS THAN 100.              ";

  const char *WRMSG0 = " ## SUBROUTINE GFWRT : WARNING                             ";
  const char *WRMSG1 = " ENDIAN WILl BE REVERSED WHEN WRITING IN BINARY MODE       ";

  int   IRW=2;

  int   MAX2;
  char  FILENM[61]=" ";
  int   i,j1,j2;
  int   NUMd,NUM2d;
  int   IDUM;
  float FDUM;
  double DDUM;
  *IERR = 0;

  MAX2=*max2d;  
  i=0;while( *(FILENMd+i) != ' ' && i<61){*(FILENM+i)=*(FILENMd+i);i++;}

  /*[1] CHECK PARAMTER */

  if(     *ITYPE!=OPEN   && *ITYPE!=CLOSE 
       && *ITYPE!=NEWSET && *ITYPE!=INTARY 
       && *ITYPE!=FLTARY && *ITYPE!=DBLARY   ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP1);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if( ( *ITYPE==INTARY && *ITYPE==FLTARY && *ITYPE==DBLARY) 
      &&
      ( *ISRCE!=ARY2D && *ISRCE!=ARY1D && *ISRCE!=ARY0D  )  ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP1);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if( *IMODE!=ASCII && *IMODE!=UNFMT && *IMODE!=BINARY){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP1);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if( ( *ITYPE==INTARY && *ITYPE==FLTARY && *ITYPE==DBLARY) 
      &&
         (    (*ISRCE==ARY2D && *NUM2 > MAX2         )
           || (*ISRCE==ARY1D && *NUM2 > 3            )  
           || (*ISRCE==ARY0D && (*NUM2>1 || *NUM>1)  )  
         )
    ){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP2);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if(  *IUNIT >= FILEMAX ){ 
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP4);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if( *IMODE==BINARY && ( *IREV!=0 && *IREV!=1 ) ){ 
      sprintf(MESSAGE,ERMSGA);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,ERRCV1);printstd(IUT0,MESSAGE);
      *IREV=0;
  }

  /*[2] OPEN */ 
  if(*ITYPE==OPEN){

      if(*IMODE==BINARY){
          if( (fp[*IUNIT] = fopen(FILENM,"w"))  == NULL ){
              sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
              sprintf(MESSAGE,EREXP3);printstd(IUT0,MESSAGE);
              *IERR = 1;
              return;
          }
      }
      else{
          open_f_(IUNIT,IMODE,FILENMd,IERR);
          if(*IERR!=0){
              sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
              sprintf(MESSAGE,EREXP3);printstd(IUT0,MESSAGE);
              *IERR = 1;
              return;
          }
      }

      j1=0;
      if(*IMODE==BINARY){
          fwrite ((char *)HEADER[BINARY],8,1,fp[*IUNIT] );
          if(*IREV==0)fwrite ((char *)(IENDIAN+0),4,1,fp[*IUNIT] );
          if(*IREV==1)fwrite ((char *)(IENDIAN+1),4,1,fp[*IUNIT] );
          NUM2d=*NUM2;
          if(*IREV==1){ conv(&NUM2d,4); } 
          fwrite ((char *)&NUM2d   ,4,1,fp[*IUNIT] );
          for(i=0;i<*NUM2;i++){ fwrite((char *)COMMNT+i*60,1,60,fp[*IUNIT] ); }
      }
      if(( *IMODE==ASCII || *IMODE==UNFMT )){
          rw_c1_08_(&IRW,IMODE,IUNIT,HEADER[*IMODE],IERR);
          rw_i1_(&IRW,IMODE,IUNIT,&ONE,&ONE,NUM2,IERR); 
          for(i=0;i<*NUM2;i++){ rw_c1_60_(&IRW,IMODE,IUNIT,COMMNT+i*60,IERR); }
          if(*IERR!=0)return;
      }

      if(  *IMODE==BINARY &&  *IREV==1  ){ 
          sprintf(MESSAGE,WRMSG0);printstd(IUT0,MESSAGE);
          sprintf(MESSAGE,WRMSG1);printstd(IUT0,MESSAGE);
      }
 
      return; 

  }/* if(*ITYPE==OPEN) */

  /*[3] NEWSET */
  if(*ITYPE==NEWSET){ 

      j1=0;
      if(*IMODE==BINARY){
          fwrite ((char *)HEADER[NEWSET],8,1,fp[*IUNIT] );
          NUM2d=*NUM2;
          if(*IREV==1){ conv(&NUM2d,4); } 
          fwrite ((char *)&NUM2d   ,4,1,fp[*IUNIT] );
          for(i=0;i<*NUM2;i++){ fwrite ((char *)COMMNT+i*60,1,60,fp[*IUNIT] ); }
      }
      if(( *IMODE==ASCII || *IMODE==UNFMT )){
          rw_c1_08_(&IRW,IMODE,IUNIT,HEADER[NEWSET],IERR);
          rw_i1_(&IRW,IMODE,IUNIT,&ONE,&ONE,NUM2,IERR); 
          for(i=0;i<*NUM2;i++){ rw_c1_60_(&IRW,IMODE,IUNIT,COMMNT+i*60,IERR); }
          if(*IERR!=0)return;
      }
      return; 
  }/* if(*ITYPE==NEWSET) */

  /*[4] CLOSE */
  if(*ITYPE==CLOSE){

      if(*IMODE==BINARY){
          fwrite ((char *)HEADER[CLOSE],1,8,fp[*IUNIT] );
          fclose(fp[*IUNIT]);
      }
      if(( *IMODE==ASCII || *IMODE==UNFMT )){
          rw_c1_08_(&IRW,IMODE,IUNIT,HEADER[CLOSE],IERR);
          close_f_(IUNIT,IERR);
          if(*IERR!=0)return;
      }
      return; 
  }/* if(*ITYPE==CLOSE) */

  /*[5] DATA (INTARY,FLTARY,DBLARY) */
  if(*ITYPE==INTARY || *ITYPE==FLTARY|| *ITYPE==DBLARY){

      if(*IMODE==BINARY){
          fwrite ((char *)HEADER[*ITYPE], 8,1,fp[*IUNIT] );
          fwrite ((char *)KEYWRD        , 8,1,fp[*IUNIT] );
          fwrite ((char *)NAME          ,30,1,fp[*IUNIT] );
          NUM2d=*NUM2;
           if(*IREV==1){ conv(&NUM2d,4); } 
          fwrite ((char *)&NUM2d   ,4,1,fp[*IUNIT] );
          NUMd =*NUM ;
          if(*IREV==1){ conv(&NUMd ,4); } 
          fwrite ((char *)&NUMd    ,4,1,fp[*IUNIT] );
      }
      if(( *IMODE==ASCII || *IMODE==UNFMT )){
          rw_c1_08_(&IRW,IMODE,IUNIT,HEADER[*ITYPE],IERR);
          rw_c1_08_(&IRW,IMODE,IUNIT,KEYWRD        ,IERR);
          rw_c1_30_(&IRW,IMODE,IUNIT,NAME          ,IERR);
          rw_i2_(&IRW,IMODE,IUNIT,&ONE,&ONE,NUM2,NUM,IERR); 
          if(*IERR!=0)return;
      } /* if(*IMODE==BINARY) */

      if( *NUM2<=0 || *NUM<=0 ) goto line150;

      if( *ISRCE == ARY0D){
          
          IDUM=*IVAR;
          FDUM=*FVAR;
          if(*ITYPE==DBLARY)DDUM=*(double *)FVAR; 

          if( *ITYPE==INTARY && *IMODE==BINARY ){
              if(*IREV==1){ conv(&IDUM,4); } 
              fwrite ((char *)&IDUM ,4,1,fp[*IUNIT] );
 	  }
          if( *ITYPE==INTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
              rw_i1_(&IRW,IMODE,IUNIT,&ONE,&ONE,&IDUM,IERR); 
 	  }
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              if(*IREV==1){ conv(&FDUM,4); } 
              fwrite ((char *)&FDUM ,4,1,fp[*IUNIT] );
	  }
          if( *ITYPE==FLTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
              rw_f1_(&IRW,IMODE,IUNIT,&ONE,&ONE,&FDUM,IERR); 
	  }
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              if(*IREV==1){ conv(&DDUM,8); } 
              fwrite ((char *)&DDUM ,8,1,fp[*IUNIT] );
	  }
          if( *ITYPE==DBLARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
              rw_d1_(&IRW,IMODE,IUNIT,&ONE,&ONE,&DDUM,IERR); 
	  } 

      if(*IERR!=0)return;    
      goto line150;
      } /* if( *ISRCE == ARY0D) */

      if( *NUM2 == 1 && *ISRCE == ARY1D ){

          j1=0; 
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  IDUM=*(IARY1+i);
                  if(*IREV==1){ conv(&IDUM,4); } 
                  fwrite ((char *)&IDUM ,4,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==INTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_i1_(&IRW,IMODE,IUNIT,NUM,NUM,IARY1,IERR); 
	  }
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  FDUM=*(FARY1+i);
                  if(*IREV==1){ conv(&FDUM,4); } 
                  fwrite ((char *)&FDUM ,4,1,fp[*IUNIT] );
              }
	  }
          if( *ITYPE==FLTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_f1_(&IRW,IMODE,IUNIT,NUM,NUM,FARY1,IERR); 
	  }
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  DDUM=*(double *)(FARY1+i*2);
                  if(*IREV==1){ conv(&DDUM,8); } 
                  fwrite ((char *)&DDUM ,8,1,fp[*IUNIT] );
              }
	  }
          if( *ITYPE==DBLARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_d1_(&IRW,IMODE,IUNIT,NUM,NUM,FARY1,IERR); 
	  }
          
      if(*IERR!=0)return;    
      goto line150;
      } /* if( *NUM2 == 1 && *ISRCE == ARY1D ) */

      if( *NUM2 == 2 && *ISRCE == ARY1D ){
          
          j1=0;
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  if((i%2)==0)IDUM=*(IARY1+(i-0)/2);
                  if((i%2)==1)IDUM=*(IARY2+(i-1)/2);
                  if(*IREV==1){ conv(&IDUM,4); } 
                  fwrite ((char *)&IDUM ,4,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==INTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_i2_(&IRW,IMODE,IUNIT,NUM,NUM,IARY1,IARY2,IERR);  
 	  }
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  if((i%2)==0)FDUM=*(FARY1+(i-0)/2);
                  if((i%2)==1)FDUM=*(FARY2+(i-1)/2);
                  if(*IREV==1){ conv(&FDUM,4); } 
                  fwrite ((char *)&FDUM ,4,1,fp[*IUNIT] );
              }
	  }
          if( *ITYPE==FLTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_f2_(&IRW,IMODE,IUNIT,NUM,NUM,FARY1,FARY2,IERR);  
          }
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  if((i%2)==0)DDUM=*(double *)(FARY1+((i-0)/2)*2);
                  if((i%2)==1)DDUM=*(double *)(FARY2+((i-1)/2)*2);
                  if(*IREV==1){ conv(&DDUM,8); } 
                  fwrite ((char *)&DDUM ,8,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==DBLARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_d2_(&IRW,IMODE,IUNIT,NUM,NUM,FARY1,FARY2,IERR);  
          }
         
      if(*IERR!=0)return;    
      goto line150;
      } /* if( *NUM2 == 2 && *ISRCE == ARY1D ) */

      if( *NUM2 == 3 && *ISRCE == ARY1D ){

          j1=0;
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  if((i%3)==0)IDUM=*(IARY1+(i-0)/3);
                  if((i%3)==1)IDUM=*(IARY2+(i-1)/3);
                  if((i%3)==2)IDUM=*(IARY3+(i-2)/3);
                  if(*IREV==1){ conv(&IDUM,4); } 
                  fwrite ((char *)&IDUM ,4,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==INTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_i3_(&IRW,IMODE,IUNIT,NUM,NUM,IARY1,IARY2,IARY3,IERR);  
          }
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  if((i%3)==0)FDUM=*(FARY1+(i-0)/3);
                  if((i%3)==1)FDUM=*(FARY2+(i-1)/3);
                  if((i%3)==2)FDUM=*(FARY3+(i-2)/3);
                  if(*IREV==1){ conv(&FDUM,4); } 
                  fwrite ((char *)&FDUM ,4,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==FLTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_f3_(&IRW,IMODE,IUNIT,NUM,NUM,FARY1,FARY2,FARY3,IERR);  
          }
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  if((i%3)==0)DDUM=*(double *)(FARY1+((i-0)/3)*2);
                  if((i%3)==1)DDUM=*(double *)(FARY2+((i-1)/3)*2);
                  if((i%3)==2)DDUM=*(double *)(FARY3+((i-2)/3)*2);
                  if(*IREV==1){ conv(&DDUM,8); } 
                  fwrite ((char *)&DDUM ,8,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==DBLARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_d3_(&IRW,IMODE,IUNIT,NUM,NUM,FARY1,FARY2,FARY3,IERR);  
          }
         
      if(*IERR!=0)return;    
      goto line150;
      } /* if( *NUM2 == 1 && *ISRCE == ARY1D ) */

      if( *ISRCE == ARY2D ){

          j1=0;
          j2=0;
          if( *ITYPE==INTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  IDUM=*(IARY2D+j1*MAX2+j2);                
                  j2++;if(j2==*NUM2){j2=0;j1++;}
                  if(*IREV==1){ conv(&IDUM,4); } 
                  fwrite ((char *)&IDUM ,4,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==INTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_in_(&IRW,IMODE,IUNIT,NUM,max2d,NUM,NUM2,IARY2D,IERR);  
          }
          if( *ITYPE==FLTARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  FDUM=*(FARY2D+j1*MAX2+j2);                
                  j2++;if(j2==*NUM2){j2=0;j1++;}
                  if(*IREV==1){ conv(&FDUM,4); } 
                  fwrite ((char *)&FDUM ,4,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==FLTARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_fn_(&IRW,IMODE,IUNIT,NUM,max2d,NUM,NUM2,FARY2D,IERR);  
          }
          if( *ITYPE==DBLARY && *IMODE==BINARY ){
              for(i=0;i< *NUM * (*NUM2);i++){
                  DDUM=*(double *)(FARY2D+(j1*MAX2+j2)*2);                
                  j2++;if(j2==*NUM2){j2=0;j1++;}
                  if(*IREV==1){ conv(&DDUM,8); } 
                  fwrite ((char *)&DDUM ,8,1,fp[*IUNIT] );
              }
          }
          if( *ITYPE==DBLARY && ( *IMODE==ASCII || *IMODE==UNFMT ) ){
               rw_dn_(&IRW,IMODE,IUNIT,NUM,max2d,NUM,NUM2,FARY2D,IERR);  
          }
          
      if(*IERR!=0)return;    
      goto line150;
      } /* if( *ISRCE == ARY2D ) */

      line150: ;
      return; 
 
  }/* if(*ITYPE==INTARY || *ITYPE==FLTARY) */

}
