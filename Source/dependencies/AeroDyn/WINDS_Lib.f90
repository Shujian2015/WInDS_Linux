!==================================================================================================================================
! Some subroutines for debug
! Writen by Shujian Liu(shujian.liu@hotmail.com)
! Lase edited April 12, 2014
!................................................................................
!
! Notice that those data are stored in C:\WInDS_for_debug
! 
!==================================================================================================================================

MODULE WINDS_Library


   USE NWTC_Library
   USE AeroDyn_Types
   USE IFPORT   ! To use MAKEDIRQQ
   
   IMPLICIT        NONE


      ! ..... Public Subroutines ............
      PUBLIC :: SAVE_TO_TXT_number      
      PUBLIC :: SAVE_TO_TXT_1D
      PUBLIC :: SAVE_TO_TXT_2D

      
   
   
CONTAINS   
!==================================================================================================================================      
SUBROUTINE SAVE_TO_TXT_number(p, number_TO_PRINT, Number_NAME)   

   IMPLICIT                      NONE
      
      ! Passed Variables:
   REAL(DbKi),                    INTENT(IN   )  :: number_TO_PRINT 
   CHARACTER(LEN = *),            INTENT(IN   )  :: Number_NAME
   TYPE(AD_ParameterType),       INTENT(IN   )  :: P           ! Parameters
   !TYPE(AD_OtherStateType),      INTENT(IN   )  :: O!therState ! Initial other/optimization states    


   
   LOGICAL               :: FILE_EXISTS   
   CHARACTER(LEN = 1024) :: PATH
   
   LOGICAL               :: DIR_EXISTS
   LOGICAL(4)            :: RESULT_CREATE
   
            
   ! Create root DIR for all simulations   
   path   =  TRIM(p%FVM%WINDS_dir) // 'WInDS_for_debug/'   


   
   INQUIRE (DIRECTORY = TRIM(path), EXIST = DIR_EXISTS)
      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(path))       
   ENDIF     

   
   
   INQUIRE(FILE=TRIM(path)//Number_NAME//'.txt', EXIST=FILE_EXISTS)   ! file_exists will be TRUE if the file
                                                  ! exists and FALSE otherwise
   
   IF (file_exists) THEN
      OPEN (UNIT=2, FILE=TRIM(path)//Number_NAME//'.txt', ACTION="write", STATUS="OLD") 
   ELSE
      OPEN (UNIT=2, FILE=TRIM(path)//Number_NAME//'.txt', ACTION="write", STATUS="NEW")     
   END IF
   
   
   WRITE(2, '(1600F14.7)')( Number_TO_PRINT )

   CLOSE(2)



END SUBROUTINE SAVE_TO_TXT_number   

!==================================================================================================================================      
SUBROUTINE SAVE_TO_TXT_1D(ARRARY_TO_PRINT, ARRY_NAME)   

   IMPLICIT                      NONE
      
      ! Passed Variables:
   REAL(DbKi),DIMENSION(:),     INTENT(IN   )  :: ARRARY_TO_PRINT  
   CHARACTER(LEN = *),            INTENT(IN   )  :: ARRY_NAME
   !TYPE(AD_ParameterType),       INTENT(IN   )  :: P           ! Parameters
   !TYPE(AD_OtherStateType),      INTENT(IN   )  :: O!therState ! Initial other/optimization states    

      ! Local variables
   INTEGER(IntKi)        :: J_dim
   INTEGER(IntKi)        :: J   
   
   LOGICAL               :: FILE_EXISTS   
   CHARACTER(LEN = 1024) :: PATH
   LOGICAL               :: DIR_EXISTS
   LOGICAL(4)            :: RESULT_CREATE
   
   J_dim   =   SIZE(ARRARY_TO_PRINT,1)


   path='/WInDS_for_debug/'
   
   
   INQUIRE (DIRECTORY = TRIM(path), EXIST = DIR_EXISTS)
      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(path))       
   ENDIF     
   
   INQUIRE(FILE=TRIM(path)//ARRY_NAME//'.txt', EXIST=FILE_EXISTS)   ! file_exists will be TRUE if the file
                                                  ! exists and FALSE otherwise
   
   IF (file_exists) THEN
      OPEN (UNIT=2, FILE=TRIM(path)//ARRY_NAME//'.txt', ACTION="write", STATUS="OLD") 
   ELSE
      OPEN (UNIT=2, FILE=TRIM(path)//ARRY_NAME//'.txt', ACTION="write", STATUS="NEW")     
   END IF
   
   
   WRITE(2, '(1600F14.7)')( ARRARY_TO_PRINT(J) ,J=1,J_dim)

   CLOSE(2)



END SUBROUTINE SAVE_TO_TXT_1D
!==================================================================================================================================      
SUBROUTINE SAVE_TO_TXT_2D(p, ARRARY_TO_PRINT, ARRY_NAME)   

   IMPLICIT                      NONE
      
      ! Passed Variables:
   REAL(DbKi),DIMENSION(:,:),     INTENT(IN   )  :: ARRARY_TO_PRINT  
   CHARACTER(LEN = *),            INTENT(IN   )  :: ARRY_NAME
   TYPE(AD_ParameterType),       INTENT(IN   )  :: P           ! Parameters
   !TYPE(AD_OtherStateType),      INTENT(IN   )  :: O!therState ! Initial other/optimization states    

      ! Local variables
   INTEGER(IntKi)      :: I_dim
   INTEGER(IntKi)      :: J_dim
   INTEGER(IntKi)      :: I
   INTEGER(IntKi)      :: J   
   
   LOGICAL                 :: FILE_EXISTS   
   CHARACTER(LEN = 1024)   :: PATH
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE
   
   
   I_dim   =   SIZE(ARRARY_TO_PRINT,1)
   J_dim   =   SIZE(ARRARY_TO_PRINT,2)


      ! Create root DIR for all simulations  
   path   =  TRIM(p%FVM%WINDS_dir) //  'WInDS_outputs'   
   INQUIRE (DIRECTORY = TRIM(path), EXIST = DIR_EXISTS)      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(path))       
   ENDIF        
   
   path   =  TRIM(p%FVM%WINDS_dir) //  'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME)    
   INQUIRE (DIRECTORY = TRIM(path), EXIST = DIR_EXISTS)      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(path))       
   ENDIF     
   
   path   =  TRIM(p%FVM%WINDS_dir) //  'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME) //'/' // TRIM(p%FVM%CURRENT_TIME)//'_for_debug/'     
   INQUIRE (DIRECTORY = TRIM(path), EXIST = DIR_EXISTS)      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(path))       
   ENDIF       
   
   
   INQUIRE(FILE = TRIM(path)//ARRY_NAME//'.txt', EXIST=FILE_EXISTS)   ! file_exists will be TRUE if the file
                                                  ! exists and FALSE otherwise
   
   IF (file_exists) THEN
      OPEN (UNIT=120, FILE=TRIM(path)//ARRY_NAME//'.txt', ACTION="write", STATUS="OLD") 
   ELSE
      OPEN (UNIT=120, FILE=TRIM(path)//ARRY_NAME//'.txt', ACTION="write", STATUS="NEW")     
   END IF
   
   
   DO I=1,I_dim
      WRITE(120, '(1600F14.7)')( ARRARY_TO_PRINT(I,J) ,J=1,J_dim)
      !WRITE(120, '(1600d17.10)')( ARRARY_TO_PRINT(I,J) ,J=1,J_dim)
   END DO

   CLOSE(120)



END SUBROUTINE SAVE_TO_TXT_2D
!================================================================================   

END MODULE WINDS_Library
!**********************************************************************************************************************************