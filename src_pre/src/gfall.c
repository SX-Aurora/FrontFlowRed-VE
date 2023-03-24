/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*   GENERAL FILE UTILITY FOR ACCESSING ANY TYPE OF DATA FILE               */
/*                                                                          */
/*      Written by  Y. Yamade                                               */
/*      Since       9 Oct. 2002                                             */
/*      Reference   gfmesh.f gfflow.f gfbou2.f gfbou3.f gfattr.f            */
/*                  gfover.f gfover.f gfsurf.f gfhist.f gfturb.f of LES3D   */
/*                                                                          */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                                          */
/*    ARGUMENT LISTINGS                                                     */
/*      (1) INPUT                                                           */
/*         IUT0        ; FILE NUMBER TO WRITE ERROR MESSAGE                 */
/*         IUT6        ; FILE NUMBER TO WRITE CALCULATION SEQUENCE          */
/*         IUNIT       ; FILE NUMBER TO ACCESS THE SPECIFIED FILE           */
/*         FILE        ; FILE NAME   TO BE OPENED                           */
/*         MCOMFL      ; MAX. NUMBER OF FILE COMMENT STRINGS                */
/*         MCOMST      ; MAX. NUMBER OF SET  COMMENT STRINGS                */
/*         IACT        ; SPECIFIES TYPE OF ACTION TO BE TAKEN               */
/*          (READ)  1 --- OPEN THE SPECIFIED FILE, READ FILE COMMENT        */
/*                       STRINGS, SEARCH AND READ A FLOW DATA SET, AND      */
/*                       CLOSE THE FILE                                     */
/*          (WRITE) 2 --- OPEN THE SPECIFIED FILE IN THE SPECIFIED MODE,    */
/*                       WRITE FILE COMMENT STRINGS, WRITE A FLOW DATA      */
/*                       SET, AND CLOSE THE FILE                            */
/*          (READ)  3 --- OPEN THE SPECIFIED FILE AND READ FILE COMMENT     */
/*                       STRINGS                                            */
/*          (WRITE) 4 --- OPEN THE SPECIFIED FILE IN THE SPECIFIED MODE,    */
/*                       AND WRITE FILE COMMENT STRINGS                     */
/*          (READ)  5 --- SEARCH AND READ THE NEXT FLOW DATA SET            */
/*          (WRITE) 6 --- APPEND-WRITE A FLOW DATA SET                      */
/*          (READ)  7 --- CLOSE THE FILE                                    */
/*          (WRITE) 8 --- CLOSE THE FILE                                    */
/*          NOTES ; GFALL  IS A SUBROUTINE WHICH READS OR WRITES A FLOW     */
/*                 DATA SET BY SEQUENTIALLY ACCESSING A SPECIFIED           */
/*                 GENERAL FILE. A FLOW DATA SET IS A SET OF DATA THAT      */
/*                 GFALL  DEFINES INTERNALLY.                               */
/*                  FOR IACT=1 AND 5, GFALL  WILL SEARCH AND READ A FLOW    */
/*                 DATA SET AND RETURN THE SEQUENCE TO THE CALLING          */
/*                 PROGRAM UNIT AS SOON AS ALL THE FLOW DATA HAVE BEEN      */
/*                 READ-IN ( AND THE FILE HAS BEEN CLOSED FOR IACT=1 ).     */
/*                                                                          */
/*                  A SET OF FLOW DATA MUST BE FOUND WITHIN A SINGLE        */
/*                 SET OF A GENERAL FILE. WHEN GFALL  ENCOUNTERS A          */
/*                 DEFINITION OF A NEW SET, IT WILL DISCARD ALL THE DATA    */
/*                 READ-IN SO FAR AND WILL RESUME THE SEARCH FROM THAT      */
/*                 POINT.                                                   */
/*                  WHEN GFALL  ENCOUNTERS A DATA THAT HAS ALREADY BEEN     */
/*                 READ-IN ( MULTIPLE DEFINITION ), IT WILL DISCARD THE     */
/*                 OLD DATA AND READ-IN THE NEWLY ENCOUNTERED DATA.         */
/*                  FOR ACTION TYPE READ, GFALL  INTERNALLY KEEPS TWO       */
/*                 VARIABLES ( 'IFPOS' AND 'IMODE' ) WHICH INDICATE THE     */
/*                 CURRENT FILE POSITION AND THE READING MODE. THUS,        */
/*                 NO MORE THAT ONE FILE CAN BE ACCESSED CONCURRENTLY       */
/*                 BY GFALL .                                               */
/*         IWRITE      ; SPECIFIES, FOR ACTION TYPE WRITE, WRITING MODE     */
/*                  1 --- ASCII                                             */
/*                  2 --- UNFORMAT                                          */
/*          NOTES ; FOR ACTION TYPE READ, ARGUMENT 'IWRITE' DOES NOT        */
/*                 HAVE TO BE SPECIFIED.                                    */
/*         INAME       ; SPECIFIES, FOR ACTION TYPE WRITE, WHETHER TO       */
/*                      USE THE INTERNALLY DEFINED DEFAULT NAMES OR         */
/*                      THE PASSED NAMES.                                   */
/*                  1 --- USE DEFAULT NAMES                                 */
/*                  2 --- USE PASSED  NAMES                                 */
/*          NOTES ; FOR ACTION TYPE READ, ARGUMENT 'INAME'  DOES NOT        */
/*                 HAVE TO BE SPECIFIED.                                    */
/*         ICAST   : SPECIFY THE ACTION WHEN THE DATA TYPE IS DIFFERENT     */
/*                   BETWEEN FILE-DATA AND ARGUMENT                         */
/*                  0 --- SET THE ERROR FLAG 'IERR'                         */
/*                  1 --- TEH READ-IN DATA WILL BE CASTED TO THE TYPE OF    */
/*                        THE ARGUMENT                                      */
/*         IDATA0  : SPECIFY THE ACTION WHEN THE DATA SIZE IS ZERO          */
/*                  0 --- NOT OUTPUT THE DATA WHOSE SIZE IS ZERO            */
/*                  1 ---     OUTPUT THE DATA WHOSE SIZE IS ZERO            */
/*          NOTES ; FOR ACTION TYPE READ, ARGUMENT 'IDATA0' DOES NOT        */
/*                 HAVE TO BE SPECIFIED.                                    */
/*         IALL    : SPECIFY THE ACTION WHEN ALL THE DATA CANNOT BE FOUND   */
/*                  0 --- SET IACT=7 AND RETURN                             */
/*                  1 --- SET THE ERROR FLAG 'IERR'                         */
/*          NOTES ; IF IALL IS NOT SET TO 1 ,GFALL  WILL NOT SET THE ERROR  */
/*                 FLAG 'IERR' EVEN IF IT COULD NOT READ ALL THE DATA       */
/*                 TO THE END OF THE FILE ( INCOMPLETE READ ), BECAUSE      */
/*                 THE CALLING PROGRAM UNIT IS NOT ALWAYS SURE THAT THERE   */
/*                 IS A COMPLETE SET OF DATA TO THE END OF THE FILE.        */
/*                 IN SUCH A CASE, GFALL  WILL SET THE VALUE OF 'IACT'      */
/*                 ARGUMENT TO 7 IN ORDER TO INFORM THE CALLING PROGRAM     */
/*                 UNIT OF THE INCOMPLETE READ. THE CALLING PROGRAM UNIT    */
/*                 CAN IDENTIFY THOSE DATA NOT READ-IN BY CHECKING 'NAME'   */
/*                 ARGUMENTS RETURNED UNSET.                                */
/*                 BUT THERE IS THE CASE THAT ALL THE DATA SHOULD BE READ   */  
/*                 , FOR EXAMPLE IT IS IMPOSSIBLE TO CONTINUE THE RUN       */
/*                 WITHOUT MESH DATA. WHEN READING SUCH CERTAINLY  NEEDED   */
/*                 DATA, IACT SHOULE BE SET TO 1                            */     
/*          NOTES ; FOR ACTION TYPE WRITE, ARGUMENT 'IALL' DOES NOT         */
/*                 HAVE TO BE SPECIFIED.                                    */
/*         ICHECK ; THIS ARGUMENT MUST BE SET AT THE END OF GFALL ARGUMENT  */
/*                 LIST , AND HAVE THE VALUE '999999'. IT IS THE FLAG FOR   */
/*                 CHECKING THAT NON-NAMED ARGUMENT LIST IS CORRELCTLLY     */
/*                 SET.                                                     */
/*         fmt     ; DATA KEYWORDS LIST.                                    */
/*                  NOTES ;                                                 */
/*                    KEWWOR IS DEFINED EXTERNALLY IN in gfc.h. IF          */
/*                   UNDEFINED KEYWORD IS SET TO fmt, THE EXCUTION WILL     */
/*                   OUTPUT ERROR MESSAGE AND TERMINATE. NUMBER OF          */
/*                   KEYWORDS MUET BE LESS THAN 100. THIS LIMITATION IS     */
/*                   DEFINED AS 'NDATA' IN THIS FILE.                       */  
/*                                                                          */
/*      (2) OUTPUT                                                          */
/*         IERR        ; RETURN CODE WHOSE VALUE WILL BE EITHER             */
/*                  0 --- INDICATING SUCCESSFUL TERMINATION                 */
/*               OR 1 --- INDICATING OCCURENCE OF SOME ERROR CONDITIONS     */
/*                                                                          */
/*      (3) INPUT-OUTPUT                                                    */
/*          NCOMFL      ; NUMBER OF FILE COMMENT STRINGS                    */
/*          COMFLE(ICOM); FILE COMMENT STRINGS                              */ 
/*          NCOMST      ; NUMBER OF SET  COMMENT STRINGS                    */
/*          COMSET(ICOM); SET  COMMENT STRINGS                              */
/*                                                                          */ 
/*         ...         ; NON-NAMED ARGUMENT LIST.                           */
/*                  NOTES ;                                                 */
/*                      THIS ARGUMEN HAVE TO BE SET DEPENDING ON THE DATA   */ 
/*                      KEYWORDS LIST (fmt). KEYWORDS IS CLASSIFIED INTO    */
/*                      FOLLOWING DATA TYPES:                               */
/*              TYPEI0: INTEGER  0-D VARIABLE  :IVAR(   1,  1)              */
/*              TYPEI1: INTEGER  1-D VARIABLE  :IVAR(   1,NUM)              */
/*              TYPEI2: INTEGER  1-D VARIABLE  :IVAR(   2,NUM)              */
/*              TYPEI3: INTEGER  1-D VARIABLE  :IVAR(   3,NUM)              */
/*              TYPEIN: INTEGER  2-D VARIABLE  :IVAR(NUM2,NUM)              */
/*              TYPEF0: FLOATING 0-D VARIABLE  :FVAR(   1,  1)              */
/*              TYPEF1: FLOATING 1-D VARIABLE  :FVAR(   1,NUM)              */
/*              TYPEF2: FLOATING 1-D VARIABLE  :FVAR(   2,NUM)              */
/*              TYPEF3: FLOATING 1-D VARIABLE  :FVAR(   3,NUM)              */
/*              TYPEFN: FLOATING 2-D VARIABLE  :FVAR(NUM2,NUM)              */
/*                                                                          */
/*                      FOLLOWING ARGUMENT IS NEED FOR EACH DATA TYPE       */
/*              TYPEI0: NAME,                   IVAR-0D                     */
/*              TYPEI1: NAME,MAX,    ,NUM ,     IVAR-1D(1    )              */
/*              TYPEI2: NAME,MAX,    ,NUM ,     IVAR-1D(1,2  )              */
/*              TYPEI3: NAME,MAX,    ,NUM ,     IVAR-1D(1,2,3)              */
/*              TYPEIN: NAME,MAX,MAX2,NUM2,NUM2,IVAR-2D                     */
/*              TYPEF0: NAME,                   FVAR-0D                     */
/*              TYPEF1: NAME,MAX,    ,NUM ,     FVAR-1D(1    )              */
/*              TYPEF2: NAME,MAX,    ,NUM ,     FVAR-1D(1,2  )              */
/*              TYPEF3: NAME,MAX,    ,NUM ,     FVAR-1D(1,2,3)              */
/*              TYPEFN: NAME,MAX,MAX2,NUM2,NUM2,FVAR-2D                     */
/*                                                                          */
/*               WHERE;                                                     */
/*               NAME IS NMAE OF DATA (CHARACTER*8)                         */  
/*               MAX,MAX2 ARE ARRAY SIZES OF ARGUMENTS                      */
/*               NUM,NUM2 ARE DATA SIZES OF ARGUMENTS                       */
/*               IVAR-0D,1D,2D ARE INTEGER  DATA ARGUMEN                    */
/*               FVAR-0D,1D,2D ARE FLOATING DATA ARGUMEN                    */
/*                                                                          */
/*                DATA TYPE IS DEFINED EXTERNALLY IN in gfc.h.              */
/*               THESE IMFOMATION OF NON-NAMED ARGUMENTS LIST WILL BE       */
/*               SET TO THE DATA TABLE 'table[i]' , WHISH IS THE STRUCTURE  */
/*               ARGUMENT DEFINED INTERNALLY IN GFALL .                     */    
/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "gfc.h"

#define NDATA  100

/* #define DEBUG0 */
/* #define DEBUG1 */
/* #define DEBUG000  */
/* #define DEBUG001  */ 

#if F2C001
void GFALL( 
#else
void gfall_( 
#endif
             int *, int *, int *, char *,
             int *, int *, char *,
             int *, int *, char *,
             int *, int *, int *, int *,int *, int *,int *,int *,
             char *, ... );


#if F2C001
void GFALL( 
#else
void gfall_( 
#endif
             int *IUT0, int *IUT6, int *IUNIT, char *FILEd
            ,int *MCOMFL, int *NCOMFL, char *COMFLE 
            ,int *MCOMST, int *NCOMST, char *COMSET 
            ,int *IACT, int *IWRITE, int *INAME,  int *IREV
            ,int *ICAST,int *IDATA0, int *IALL ,int *IERR
            ,char *fmt , ...)
{
  va_list ap         ;       /* ARGUMENT LIST                                 */
  TABLE table[NDATA] ;       /* DATE-TABLE                                    */
  int   key_list_ptr[NDATA]; /*                                               */
  int   *ICHECK      ;       /* (=999999) FLAG TO CHECK ARGUMENT SETTING      */

  char  cdefault=' ' ;       /* LOCAL CHAR  ARGUMENT TO INITIALIZE DATA-TABLE */  
  int   idefault=1   ;       /* LOCAL INT   ARGUMENT TO INITIALIZE DATA-TABLE */   
  int   idefault2=1  ;       /* LOCAL INT   ARGUMENT TO INITIALIZE DATA-TABLE */   
  int   idefault3=1  ;       /* LOCAL INT   ARGUMENT TO INITIALIZE DATA-TABLE */   
  float fdefault=0.0 ;       /* LOCAL FLOAT ARGUMENT TO INITIALIZE DATA-TABLE */  


  const char *ERMSGA = " ## SUBROUTINE GFALL : RECOVERBLE ERROR OCCURENCE; CONTINUE   ";
  const char *ERMSGB = " ## SUBROUTINE GFALL : FATAL      ERROR OCCURENCE; RETURNED   ";
  const char *ERRCV1 = " ## THE PARAMETER 'INAME' WILL BE SET TO THE DEFAULT VALUE    ";
  const char *EREXP1 = " ## KEYWORD IS NOT DEFINED IN gfc.h                           ";
  const char *EREXP2 = " ## ERROR OCCUERED WHEN MAKING DATA TABLE                     ";
  const char *EREXP3 = " ## AN ILLEGAL VALUE WAS SPECIFIED FOR A CONTROLLING PARAMETER";

  const char *WRMSG0 = " ## SUBROUTINE GFALL : WARNING                                ";
  const char *WRMSG1 = " ## GFALL: WARNING INCOMPLETE SET, DISCARD DATA               ";


  const char *OPNMSG[8] = {  " GFALL: LOOKING FOR A DATA SET                           "
                            ," GFALL: WRITING     A DATA SET                           "        
                            ," GFALL: OPENING     A DATA FILE                          "       
                            ," GFALL: OPENING     A DATA FILE                          "       
                            ," GFALL: LOOKING FOR THE NEXT DATA SET                    " 
                            ," GFALL: APPEND-WRITING     A DATA SET                    " 
                            ," GFALL: CLOSING     THE DATA FILE                        "     
                            ," GFALL: CLOSING     THE DATA FILE                        "     
                         } ;

 const char *ENDMSG[2] = {  "GFALL: SUCCESSFULLY RETURNING                             "
                           ,"GFALL: INCOMPLETE READ (END OF FILE)                      " 
                        };

 const char *ACTION[3] = {  "** PERSING ***"
                           ,"** READING OK ***"
                           ,"** WRITING OK ***" };

 const char *TYPEARG[3] = {  "INTEGER"
                            ,"SINGLE "
                            ,"DOUBLE " };

  const int  INAMED = 1;

  /*LOCAL ENTRY FOR BUFFER */ 
  char  *p                     ;
  char  headx[9]= " "          ;
  char  FILE[61]   = " "       ;
  char  CDUMMY[61]             ;
  char  NAMERD[31]             ;
  char  KEYWRD[9]= " "         ; 
  int   i,j,nbuf,idata         ;   
  int   ndata=0                ;
  int   keynum                 ; 
  int   TYPE0                  ;
  int   imatch=0               ;
  int   *NUM,*NUM2,NUMd,NUM2d  ;
  int   IERRDM                 ;
  int   *IVAR,ibuf             ;  
  float *FVAR,fbuf             ;

  /*LOCAL ENTRY FOR READ & WRITE MODE */ 
  int        ITYPE[2] ;       /*SPECIFIES THE TYPE OF ACTION TO BE TAKEN         */
  int        IDEST    ;       /*SPECIFIES DESTINATION VARIABLE/ARRAYS TO STORE   */
  int        ISRCE    ;       /*SPECIFIES SOURCE VARIABLE/ARRAYS TO BE WRITTEN   */
  static int IMODE[FILEMAX];  /*SPECIFIES THE TYPE OF INPUT  MODE                */
  static int IFPOS;           /*SPECIFIES CURRENT FILE POSITION                  */
  int JGOT[NDATA],JSTART,JALL;/*LOCAL ENTRY FOR CHECKING READ-IN DATA            */

/*********************************************************************************/
/* NOTE                                                                          */
/*  ITYPE[0] HAS A SAME SENSE AS ONE OF PREVIOUS VESION, AND ITYPE[1] IS NEEDED  */  
/* TO CHECK THE TYPE OF PASSED ARUGUMENT, WHICH IS INTEGER ,OR SINGLE FLOATING   */
/* POINT ,OR DOUBLE FLOATING POINT.                                              */
/*********************************************************************************/

#ifdef DEBUG0
for(i=1;i<=MAXKEY;i++){
printf("%2d %s  ",i,KEYWD0[i]);
switch( TYPE[i]){
case  TYPEI0 : printf("TYPEI0"); break;
case  TYPEI1 : printf("TYPEI1"); break;
case  TYPEI2 : printf("TYPEI2"); break;
case  TYPEI3 : printf("TYPEI3"); break;
case  TYPEIN : printf("TYPEIN"); break;
case  TYPEF0 : printf("TYPEF0"); break;
case  TYPEF1 : printf("TYPEF1"); break;
case  TYPEF2 : printf("TYPEF2"); break;
case  TYPEF3 : printf("TYPEF3"); break;
case  TYPEFN : printf("TYPEFN"); break;
default      : printf("NOT MATCH"); break;
}
printf("\n");
}
#endif

  /*GFALL: **********************************************************************/
  /*GFALL: [0    ] START OF EXCUTION OF GFALL */

  /*GFALL: **********************************************************************/
  /*GFALL: [1    ] PERSE FORMAT SENTENCE AND MAKE DATA-TABLE*/
  /*GFALL: [1-1  ] COUNT NUMBER OF KEYWORD */
  p=fmt;
  i=0;

  while(*p!='\0'&&*p!='!'){
      if(*p == '*'){
          key_list_ptr[ndata]=i;
          ndata++;
      } 
      p++;
      i++;
  }/*while(*p) */
  sprintf(MESSAGE,"** GFALL   *** :ndata :%d ",ndata); 
  printstd(IUT6,MESSAGE);

  /*GFALL: [1-2  ] INITIALIZE TABLE */
  for(i=0;i<ndata;i++){ 
      table[i].name   = &cdefault;
      table[i].keynum =  idefault;
      table[i].max2   = &idefault;
      table[i].max    = &idefault;
      table[i].num2   = &idefault2;
      table[i].num    = &idefault3;
      table[i].i0     = &idefault;
      table[i].i1     = &idefault;
      table[i].i2     = &idefault;
      table[i].i3     = &idefault;
      table[i].in     = &idefault;
      table[i].f0     = &fdefault;
      table[i].f1     = &fdefault;
      table[i].f2     = &fdefault;
      table[i].f3     = &fdefault;
      table[i].fn     = &fdefault;
  }/* for(i=0;i<ndata;i++)  */

  /*GFALL: [1-3  ] SEARCH KEYWORD AND CHECK KEYNUM */
  for(nbuf=0;nbuf<ndata;nbuf++){
      p=fmt+key_list_ptr[nbuf];
      headx[0]=*p; 
      for(i=0;i<7;i++){
          p++;
          headx[i+1]=*p;
      }
      imatch=0; 
      for(i=1;i<=MAXKEY;i++){
          if(strcmp(headx,KEYWD0[i])==0){
              table[nbuf].keynum=i;
              keynum=table[nbuf].keynum;  
              TYPE0=TYPE[keynum];
              if(   TYPE0==TYPEI0 || TYPE0==TYPEI1 || TYPE0==TYPEI2 
                 || TYPE0==TYPEI3 || TYPE0==TYPEIN){
                  table[nbuf].itype=INTARY;
              }
              else{
	          p++;CDUMMY[0]=*p;
	          p++;CDUMMY[1]=*p;
                  if(CDUMMY[0]=='%' && CDUMMY[1]=='D') {table[nbuf].itype=DBLARY;}
                  else                                 {table[nbuf].itype=FLTARY;}
              }
          imatch=1;
          }
      }
      if(imatch==0){
          sprintf(MESSAGE,ERMSGB);                        printstd(IUT0,MESSAGE);
          sprintf(MESSAGE," UNDEFINED KEYWORD :%s",headx);printstd(IUT0,MESSAGE);
          sprintf(MESSAGE,EREXP1);strcat(MESSAGE,headx);  printstd(IUT0,MESSAGE);
          *IERR=1; 
          return;
      }
      keynum=table[nbuf].keynum;
#ifdef DEBUG1
printf("%02d-th keynum=%2d:keyword=%8s:data-type=%s:default-name=%s \n"
,nbuf,keynum,KEYWD0[keynum],CTYPE[TYPE[keynum]],NAME0[keynum]);
#endif
  }

  /*GFALL: [1-4  ] MAKE DATA TABLE */
  va_start(ap,fmt);
  for(i=0;i<ndata;i++){
      keynum=table[i].keynum;  
      TYPE0 =table[i].itype-6;      
 
      switch(TYPE[keynum]){
          case TYPEI0: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEI0 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].max2 = &ONE;
              table[i].max  = &ONE;
              table[i].name = va_arg(ap, char *);
              table[i].i0   = va_arg(ap, int  *);
              break;

          case TYPEI1: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEI1 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char *);
              table[i].max  = va_arg(ap, int  *);
              table[i].max2 = &ONE;
              table[i].num  = va_arg(ap, int  *);
              table[i].num2 = &ONE;
              table[i].i1   = va_arg(ap, int  *);
              break;

          case TYPEI2: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEI2 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char *);
              table[i].max  = va_arg(ap, int  *);
              table[i].max2 = &TWO;
              table[i].num  = va_arg(ap, int  *);
              table[i].num2 = &TWO;
              table[i].i1   = va_arg(ap, int  *);
              table[i].i2   = va_arg(ap, int  *);
              break;

          case TYPEI3: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEI3 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char *);
              table[i].max  = va_arg(ap, int  *);
              table[i].max2 = &THREE;
              table[i].num  = va_arg(ap, int  *);
              table[i].num2 = &THREE;
              table[i].i1   = va_arg(ap, int  *);
              table[i].i2   = va_arg(ap, int  *);
              table[i].i3   = va_arg(ap, int  *);
              break;

          case TYPEIN: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEIN IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char *);
              table[i].max  = va_arg(ap, int  *);
              table[i].max2 = va_arg(ap, int  *);
              table[i].num  = va_arg(ap, int  *);
              table[i].num2 = va_arg(ap, int  *);
              table[i].in   = va_arg(ap, int  *);
              break;
 
          case TYPEF0: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEF0 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].max2 = &ONE;
              table[i].max  = &ONE;
              table[i].name = va_arg(ap, char  *);
              table[i].f0   = va_arg(ap, float *);
              break;

          case TYPEF1: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEF1 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char  *);
              table[i].max  = va_arg(ap, int   *);
              table[i].max2 = &ONE;
              table[i].num  = va_arg(ap, int   *);
              table[i].num2 = &ONE;
              table[i].f1   = va_arg(ap, float *);
              break;

          case TYPEF2: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEF2 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char *);
              table[i].max  = va_arg(ap, int   *);
              table[i].max2 = &TWO;
              table[i].num  = va_arg(ap, int   *);
              table[i].num2 = &TWO;
              table[i].f1   = va_arg(ap, float *);
              table[i].f2   = va_arg(ap, float *);
              break;

          case TYPEF3: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEF3 IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char *);
              table[i].max  = va_arg(ap, int   *);
              table[i].max2 = &THREE;
              table[i].num  = va_arg(ap, int   *);
              table[i].num2 = &THREE;
              table[i].f1   = va_arg(ap, float *);
              table[i].f2   = va_arg(ap, float *);
              table[i].f3   = va_arg(ap, float *);
              break;

          case TYPEFN: 
              sprintf(MESSAGE,"%s  :%8s: POINTER OF TYPEFN IS SET AS %8s"
                             ,ACTION[0],KEYWD0[keynum],TYPEARG[TYPE0]); 
              table[i].name = va_arg(ap, char  *);
              table[i].max  = va_arg(ap, int   *);
              table[i].max2 = va_arg(ap, int   *);
              table[i].num  = va_arg(ap, int   *);
              table[i].num2 = va_arg(ap, int   *);
              table[i].fn   = va_arg(ap, float *);
              break;
 
      }/* switch(TYPEkeynum]) */

      printstd(IUT6,MESSAGE);

  }/* for(i=0;i<ndata;i++) */

  /*GFALL: [1-5  ] CHECK :ARGUMENT IS CORRECTLY SET ? */
  ICHECK = va_arg(ap, int   *);
  if(*ICHECK != 999999){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP2);printstd(IUT0,MESSAGE);
      *IERR=1;
      return;
  }
  va_end(ap);


   /*GFALL: **********************************************************************/
  /*GFALL: [2    ] PREPARE FOR READ AND WRITE MODE */
  /*GFALL: [2-1  ] DUMMY SUBSTITUTION TO AVOID WARNING OF COMPILER */
  ibuf=0 ; IVAR=&ibuf;
  fbuf=0 ; FVAR=&fbuf;
  NUMd=0 ; NUM =&NUMd ; 
  NUM2d=0; NUM2=&NUM2d; 

  /*GFALL: [2-2  ] SET FILE NAME AS C CHARACTER CODE FORMAT */
  i=0;while( *(FILEd+i) != ' ' && i<61){*(FILE+i)=*(FILEd+i);i++;}

  /*GFALL: [2-3  ] SET INITIAL VALUES */
  *IERR = 0;
  *IREV = *IREV;
  JSTART=0;
  if(*IACT!=1 && *IACT!=5) JALL=2;

  if(*IACT%2==1) {
      for(idata=0;idata<ndata;idata++){
          *table[idata].num=0;
      }
  }


  /*GFALL: [2-4  ] CHECK PASSED CONTROL PARAMETERS */
  if(*IACT<1 || *IACT>8){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP3);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if(*IACT%2==0 && *IWRITE!=1&& *IWRITE!=2&& *IWRITE!=3){
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP3);printstd(IUT0,MESSAGE);
      *IERR = 1;
      return;
  }

  if(*IACT%2==0 && *INAME!=1&& *INAME!=2){
      sprintf(MESSAGE,ERMSGA);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,EREXP3);printstd(IUT0,MESSAGE);
      sprintf(MESSAGE,ERRCV1);printstd(IUT0,MESSAGE);
      *INAME=INAMED;
  }

  /*GFALL: **********************************************************************/
  sprintf(MESSAGE,"%s",OPNMSG[*IACT-1]);          printstd(IUT6,MESSAGE);
  sprintf(MESSAGE," FILENAME       =   %s",FILE); printstd(IUT6,MESSAGE);

  /*GFALL: [3    ] ACTION TYPE READ */
  if(*IACT%2==1){
    
  /*GFALL: [3-1  ] OPEN THE FILE */
      if(*IACT==1 || *IACT==3){
          ITYPE[0]=OPEN;
          gfrd_(ITYPE,&IDEST,&IMODE[*IUNIT],&IFPOS,ICAST,
                IUNIT,FILE,KEYWRD,NAMERD,COMFLE,
                IVAR,IVAR,IVAR,IVAR,IVAR,
                FVAR,FVAR,FVAR,FVAR,FVAR,
                NUM,NCOMFL,&ONE,MCOMFL,IUT0,IERR);
          if(ONE!=1)*IERR=1;
          if(*IERR!=0)goto line999;

          p=COMFLE;
          for(i=0;i<*NCOMFL;i++){
              for(j=0;j<59;j++)CDUMMY[j]=p[j];
              CDUMMY[59] = '\0';
              sprintf(MESSAGE,"FILE COMMENT %s ",CDUMMY);
              printstd(IUT6,MESSAGE); 
              p=p+60;
          }

      }/* if(*IACT==1 || *IACT==3) */

  /*GFALL: [3-2  ] SEARCH AND READ THE NEXT FLOW DATA SET */
      if(*IACT==1 || *IACT==5){
 
      line100: ;

  /*GFALL: [3-2-1] CHECK IF ALL DATA HAVE BEEN READ-IN */
      JALL=JSTART;
      for(i=0;i<ndata;i++){
          JALL=JALL*JGOT[i];
      }
      if(JALL==1) goto line200; 

  /*GFALL: [3-2-2] INQUIRE THE NEXT CONTENT IN THE FILE */ 
      ITYPE[0]=INQUIR;
      gfrd_(ITYPE,&IDEST,&IMODE[*IUNIT],&IFPOS,ICAST,
            IUNIT,FILE,KEYWRD,NAMERD,COMSET,
            IVAR,IVAR,IVAR,IVAR,IVAR,
            FVAR,FVAR,FVAR,FVAR,FVAR,
            NUM,NCOMST,&ONE,MCOMST,IUT0,IERR);
      if(ONE!=1)*IERR=1;
      if(*IERR!=0)goto line999;

  /*GFALL: [3-2-3] END OF THE FILE */
      if(ITYPE[0]==CLOSE){
          *IACT=7;
          goto line200;    
      } 

  /*GFALL: [3-2-3] A NEW SET */
      if(ITYPE[0]==NEWSET){
          if(JSTART==1){
              sprintf(MESSAGE,WRMSG0);printstd(IUT6,MESSAGE);
              sprintf(MESSAGE,WRMSG1);printstd(IUT6,MESSAGE);
          } 

          p=COMSET;
          for(i=0;i<*NCOMST;i++){
              for(j=0;j<59;j++)CDUMMY[j]=p[j];
              CDUMMY[59] = '\0';
              sprintf(MESSAGE,"SET  COMMENT %s ",CDUMMY);
              printstd(IUT6,MESSAGE); 
              p=p+60;
          }
          for(i=0;i<ndata;i++){
              JGOT[i]=0;
          }    
          JSTART=1;
          goto line100 ;
       }/* if(ITYPE[0]==NEWSET) */ 

  /*GFALL: [3-2-4] CHECK :KEYWORD IS MATCHED IN THE LIST OF gfc.h */
      imatch=0; 
      for(i=1;i<=MAXKEY;i++){
          if(strcmp(KEYWRD,KEYWD0[i])==0){
              imatch=1;
              keynum=i;
          }
      }
      if(imatch==0){
          sprintf(MESSAGE,ERMSGB);               printstd(IUT0,MESSAGE);
          sprintf(MESSAGE,"%s %s",EREXP1,KEYWRD);printstd(IUT0,MESSAGE);
          *IERR=1; 
          return;
      }

  /*GFALL: [3-2-5] CHECK :KEYWORD IS MATCHED IN THE LIST OF ARGUMENT */
      imatch=0; 
      for(i=0;i<ndata;i++){
          if(keynum==table[i].keynum&&JGOT[i]==0){
              imatch=1; 
              idata=i;
          }  
      }

  /*GFALL: [3-2-6] READ DATA */
      if(imatch==1){
          switch(TYPE[keynum]){
              case TYPEI0: IDEST= ARY0D;break;
              case TYPEI1: IDEST= ARY1D;break;
              case TYPEI2: IDEST= ARY1D;break;
              case TYPEI3: IDEST= ARY1D;break;
              case TYPEIN: IDEST= ARY2D;break;
              case TYPEF0: IDEST= ARY0D;break;
              case TYPEF1: IDEST= ARY1D;break;
              case TYPEF2: IDEST= ARY1D;break;
              case TYPEF3: IDEST= ARY1D;break;
              case TYPEFN: IDEST= ARY2D;break;
          }

          ITYPE[1]=table[idata].itype;
	  
#ifdef DEBUG000
printf("PASS 000 \n");
printf("PASS 000 :ITYPE :%d\n",ITYPE[0] );
printf("PASS 000 :KEYWRD:%s\n",KEYWRD);
printf("PASS 000 :keynum:%d\n",keynum);
printf("PASS 000 :idata :%d\n",idata );
printf("PASS 000 :imatch:%d\n",imatch);
printf("PASS 000 :IDEST :%d\n",IDEST );
printf("PASS 000 :MAX2  :%d\n",*table[idata].max2  );
printf("PASS 000 :MAX   :%d\n",*table[idata].max   );
#endif

          i=idata;
          gfrd_(ITYPE,&IDEST,&IMODE[*IUNIT],&IFPOS,ICAST,
                IUNIT,FILE,KEYWRD,NAMERD,COMSET,
                table[i].i0,table[i].i1,table[i].i2,table[i].i3,table[i].in,
                table[i].f0,table[i].f1,table[i].f2,table[i].f3,table[i].fn,
                table[i].num,table[i].num2,table[i].max,table[i].max2,IUT0,IERR);
          if(*IERR!=0)goto line999;
          JGOT[idata]=1;

          for(j=0;j<30;j++) *(table[idata].name+j)=*(NAMERD+j);
          j=30;             *(NAMERD+j)='\0';

          sprintf(MESSAGE,"%s  %s  (%8d,%8d)" 
                 ,ACTION[1],NAMERD,*table[idata].num2,*table[idata].num);
          if(TYPE[keynum]==TYPEI0){
              sprintf(CDUMMY,"%8d"   ,*table[idata].i0);
              strcat(MESSAGE,CDUMMY) ;
          }  
          if(TYPE[keynum]==TYPEF0){
              if(table[idata].itype==FLTARY){
                  sprintf(CDUMMY,"%14.6E",*table[idata].f0);
                  strcat(MESSAGE,CDUMMY) ;
              }
              if(table[idata].itype==DBLARY){
                  sprintf(CDUMMY,"%22.14E",*(double *)table[idata].f0);
                  strcat(MESSAGE,CDUMMY) ;
              }
          }
          printstd(IUT6,MESSAGE);

          goto line100;

      }/* if(imatch==1) */

  /*GFALL: [3-2-7] SKIP TO THE NEXT DATA */
      if(imatch==0){
          ITYPE[0]=-ITYPE[0];
          gfrd_(ITYPE,&IDEST,&IMODE[*IUNIT],&IFPOS,ICAST,
                IUNIT,FILE,KEYWRD,NAMERD,COMSET,
                IVAR,IVAR,IVAR,IVAR,IVAR,
                FVAR,FVAR,FVAR,FVAR,FVAR,
                NUM,NUM2,&ONE,&ONE,IUT0,IERR);
          if(ONE!=1)*IERR=1;
          if(*IERR!=0)goto line999;

          sprintf(MESSAGE,"%s  SKIP TO THE NEXT DATA           (%8s)" ,ACTION[1],KEYWRD);
          printstd(IUT6,MESSAGE);

          goto line100;
      }/* if(imatch==0) */
    
      line200: ;
  }/* if(*IACT==1 || *IACT==5) */

  line999: ;

  /*GFALL: [3-3  ] CLOSE THE FILE */
  if(*IACT==1 || *IACT==7 || *IERR!=0){

  ITYPE[0]=CLOSE;
  gfrd_(ITYPE,&IDEST,&IMODE[*IUNIT],&IFPOS,ICAST,
        IUNIT,FILE,KEYWRD,NAMERD,COMSET,
        IVAR,IVAR,IVAR,IVAR,IVAR,
        FVAR,FVAR,FVAR,FVAR,FVAR,
        NUM,NUM2,&ONE,&ONE,IUT0,&IERRDM);
  if(ONE!=1)*IERR=1;
  /* Note : IERRDM is used to avoid resetting IERR ARGUMENT */
  }/* if(*IACT==1 || *IACT==7 || *IERR!=0) */

  if(JALL==0&&*IALL==1){
      sprintf(MESSAGE,ENDMSG[1]);printstd(IUT6,MESSAGE);
      sprintf(MESSAGE,ERMSGB);printstd(IUT0,MESSAGE);
      *IERR=1;
      return;
  }

  if(*IERR==0){
      if(JALL==0){sprintf(MESSAGE,ENDMSG[1]);}
      else       {sprintf(MESSAGE,ENDMSG[0]);}
      printstd(IUT6,MESSAGE);

  }
  else{
      sprintf(MESSAGE,ERMSGB);
      printstd(IUT0,MESSAGE);
      
  }

  return;
  }/* if(*IACT%2==1) */

  /*GFALL: **********************************************************************/
  /*GFALL: [4    ] ACTION TYPE WRITE */
  if(*IACT%2==0){

  /*GFALL: [4-1  ] OPEN THE FILE IN THE SPECIFIED MODE */
  if(*IACT==2 || *IACT==4 ){
      ITYPE[0] = OPEN;
      gfwrt_(ITYPE,&ISRCE,IWRITE,IREV,
             IUNIT,FILE,KEYWRD,NAMERD,COMFLE,
             IVAR,IVAR,IVAR,IVAR,IVAR,
             FVAR,FVAR,FVAR,FVAR,FVAR,
             NUM,NCOMFL,&ONE,IUT0,IERR);
      if(ONE!=1)*IERR=1;
      if(*IERR!=0)goto line999;

      p=COMFLE;
      for(i=0;i<*NCOMFL;i++){
          for(j=0;j<59;j++)CDUMMY[j]=p[j];
          CDUMMY[59] = '\0';
          sprintf(MESSAGE,"FILE COMMENT %s ",CDUMMY);
          printstd(IUT6,MESSAGE); 
          p=p+60;
      }

  }/* if(*IACT==2 || *IACT==4 ) */

  /*GFALL: [4-2  ] APPEND-WRITE A FLOW DATA SET */   
  if(*IACT==2 || *IACT==6 ){

      ITYPE[0]=NEWSET;
      gfwrt_(ITYPE,&ISRCE,IWRITE,IREV,
             IUNIT,FILE,KEYWRD,NAMERD,COMSET,
             IVAR,IVAR,IVAR,IVAR,IVAR,
             FVAR,FVAR,FVAR,FVAR,FVAR,
             NUM,NCOMST,&ONE,IUT0,IERR);
      if(ONE!=1)*IERR=1;
      if(*IERR!=0)goto line1999;
     
      p=COMSET;
      for(i=0;i<*NCOMST;i++){
          for(j=0;j<59;j++)CDUMMY[j]=p[j];
          CDUMMY[59] = '\0';
          sprintf(MESSAGE,"SET  COMMENT %s ",CDUMMY);
          printstd(IUT6,MESSAGE); 
          p=p+60;
      }

  /*GFALL: [4-3  ] WRITE DATA */
      for(idata=0;idata<ndata;idata++){
   
          if(*table[idata].num>0 || *IDATA0==1 ){

 	
              keynum=table[idata].keynum;
              strcpy(KEYWRD,KEYWD0[table[idata].keynum ]);
              if(*INAME==1){
                  for(j=0;j<30;j++){
                     table[idata].name[j]=*(NAME0[keynum]+j);
                  }
              }

              switch(TYPE[keynum]){
                  case TYPEI0: ITYPE[0]=INTARY; ISRCE= ARY0D ; break;
                  case TYPEI1: ITYPE[0]=INTARY; ISRCE= ARY1D ; break;
                  case TYPEI2: ITYPE[0]=INTARY; ISRCE= ARY1D ; break;
                  case TYPEI3: ITYPE[0]=INTARY; ISRCE= ARY1D ; break;
                  case TYPEIN: ITYPE[0]=INTARY; ISRCE= ARY2D ; break;
                  case TYPEF0: ITYPE[0]=FLTARY; ISRCE= ARY0D ; break;
                  case TYPEF1: ITYPE[0]=FLTARY; ISRCE= ARY1D ; break;
                  case TYPEF2: ITYPE[0]=FLTARY; ISRCE= ARY1D ; break;
                  case TYPEF3: ITYPE[0]=FLTARY; ISRCE= ARY1D ; break;
                  case TYPEFN: ITYPE[0]=FLTARY; ISRCE= ARY2D ; break;
              }
              ITYPE[0]=table[idata].itype;

#ifdef DEBUG001          
printf("PASS 001 \n");
printf("PASS 001 :ITYPE :%d\n",ITYPE[0] );
printf("PASS 001 :KEYWRD:%s\n",KEYWRD);
printf("PASS 001 :keynum:%d\n",keynum);
printf("PASS 001 :idata :%d\n",idata );
printf("PASS 001 :imatch:%d\n",imatch);
printf("PASS 001 :ISRCE :%d\n",ISRCE );
printf("PASS 001 :MAX2  :%d\n",*table[idata].max2  );
printf("PASS 001 :MAX   :%d\n",*table[idata].max   );
printf("PASS 001 :IUT0  :%d\n",*IUT0);
#endif

              i=idata;
              gfwrt_(ITYPE,&ISRCE,IWRITE,IREV,
                     IUNIT,FILE,KEYWRD,table[idata].name,COMSET,
                     table[i].i0,table[i].i1,table[i].i2,table[i].i3,table[i].in,
                     table[i].f0,table[i].f1,table[i].f2,table[i].f3,table[i].fn,
                     table[i].num,table[i].num2,table[i].max2,IUT0,IERR);
              if(*IERR!=0)goto line1999;

              for(j=0;j<30;j++) *(NAMERD+j)=*(table[idata].name+j) ;
              j=30;             *(NAMERD+j)='\0';
             
              sprintf(MESSAGE,"%s  %s  (%8d,%8d)" 
                     ,ACTION[2],NAMERD,*table[idata].num2,*table[idata].num);
              if(TYPE[keynum]==TYPEI0){
                  sprintf(CDUMMY,"%8d"   ,*table[idata].i0);
                  strcat(MESSAGE,CDUMMY) ;
              }  
              if(TYPE[keynum]==TYPEF0){
                  if(table[idata].itype==FLTARY){
                      sprintf(CDUMMY,"%14.6E",*table[idata].f0);
                      strcat(MESSAGE,CDUMMY) ;
                  }
                  if(table[idata].itype==DBLARY){
                      sprintf(CDUMMY,"%14.6E",*(double *)table[idata].f0);
                      strcat(MESSAGE,CDUMMY) ;
                  }
              }
              printstd(IUT6,MESSAGE);

          }/* if(*table[idata].num>0) */

      }/* for(idata=0;idata<ndata;idata++) */

  }/* if(*IACT==2 || *IACT==6 ) */

  line1999: ;

  /*GFALL: [4-4  ] CLOSE THE FILE */
  if(*IACT==2 || *IACT==8 ){
      ITYPE[0] = CLOSE;
      gfwrt_(ITYPE,&ISRCE,IWRITE,IREV,
             IUNIT,FILE,KEYWRD,NAMERD,COMFLE,
             IVAR,IVAR,IVAR,IVAR,IVAR,
             FVAR,FVAR,FVAR,FVAR,FVAR,
             &ONE,&ONE,&ONE,IUT0,&IERRDM);
      if(ONE!=1)*IERR=1;
      /* Note : IERRDM is used to avoid resetting IERR ARGUMENT */
  }

  if(*IERR==0){
      if(JALL==0){sprintf(MESSAGE,ENDMSG[1]);}
      else       {sprintf(MESSAGE,ENDMSG[0]);}
      printstd(IUT6,MESSAGE);
  }
  else{
      sprintf(MESSAGE,ERMSGB);
      printstd(IUT0,MESSAGE);
  }/* if(*IERR==0) */

  return;

  }/* if(*IACT%2==0) */

}
