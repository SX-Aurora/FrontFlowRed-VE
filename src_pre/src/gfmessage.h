/* INCLUDE FILE FOR GENERAL FILE VER-2.1 */

 const char *ERMSGB = " ## SUBROUTINE GFALL : FATAL      ERROR OCCURENCE; RETURNED   ";
 const char *ERRCV1 = " ## THE PARAMETER WILL BE SET TO THE DEFAULT VALUE OF         ";
 const char *EREXP1 = " ## KEYWORD IS NOT DEFINED IN gfc.h                           ";
 const char *EREXP2 = " ## ERROR OCCUERED WHEN MAKING DATA TABLE                     ";
 const char *EREXP3 = " ## AN ILLEGAL VALUE WAS SPECIFIED FOR A CONTROLLING PARAMETER";
 const char *ERMSGA = " ## SUBROUTINE GFALL : RECOVERBLE ERROR OCCURENCE; CONTINUE   ";
 const char *ERMSGC = " ## SUBROUTINE GFALL : FATAL      ERROR REPORT   ; RETURNED   "; 


 const char *OPNMSG[8] = {  " GFALL: LOOKING FOR A DATA SET                              "
                           ," GFALL: WRITING     A DATA SET                              "        
                           ," GFALL: OPENING     A DATA FILE                             "       
                           ," GFALL: OPENING     A DATA FILE                             "       
                           ," GFALL: LOOKING FOR THE NEXT DATA SET                       " 
                           ," GFALL: APPEND-WRITING     A DATA SET                       " 
                           ," GFALL: CLOSING     THE DATA FILE                           "     
                           ," GFALL: CLOSING     THE DATA FILE                           "     
                         } ;

 const char *DIMMSG[2] = { "(2D)" , "(3D)"};

 const char *ENDMSG[2] = {  "GFALL: SUCCESSFULLY RETURNING                               "
                           ,"GFALL: INCOMPLETE READ (END OF FILE)                        " 
                        };

 const char *WRNMSG[2] = {  "GFALL: WARNING INCOMPLETE SET, DISCARD DATA                 "
                           ,"GFALL: WARNING MULTIPLE DATA DEFINITIONS                    "   
                        };

 const char *ACTION[3] = {  "** PERSING ***"
                           ,"** READING OK ***"
                           ,"** WRITING OK ***" };

 const char *TYPEARG[3] = {  "INTEGER"
                            ,"SINGLE "
                            ,"DOUBLE " };

/* END OF INCLUDE FILE FOR GENERAL FILE VER-1.2 */
