MODULE tree_library
    
   USE IFPORT   ! To use MAKEDIRQQ 
    
    
   IMPLICIT        NONE
        
! Define precision  (SingPrec.f90)          
   INTEGER, PARAMETER              :: B4Ki     = SELECTED_INT_KIND(  9 )           ! Kind for four-byte whole numbers
   INTEGER, PARAMETER              :: SiKi     = SELECTED_REAL_KIND(  6,  30 )     ! Kind for four-byte, floating-point numbers
   INTEGER, PARAMETER              :: R8Ki     = SELECTED_REAL_KIND( 14, 300 )     ! Kind for eight-byte floating-point numbers

   INTEGER, PARAMETER              :: IntKi    = B4Ki                              ! Default kind for integers 
   INTEGER, PARAMETER              :: ReKi     = SiKi                              ! Default kind for floating-point numbers
   INTEGER, PARAMETER              :: DbKi     = R8Ki                        ! Default kind for double floating-point numbers
                                                     
   REAL(ReKi)                      :: pi
    
   

      ! ..... Public Subroutines ............
   PUBLIC :: EqualRealNos     
   PUBLIC :: SAVE_TO_TXT_2D
   
   
CONTAINS
!================================================
   FUNCTION EqualRealNos ( ReNum1, ReNum2 )

      ! This function compares 2 real numbers and determines if they
      ! are "almost" equal, i.e. within some relative tolerance
      ! ("Safe Comparisons" suggestion from http://www.lahey.com/float.htm)

      ! passed variables

   REAL(DbKi), INTENT(IN )         :: ReNum1                            ! the first  real number to compare
   REAL(DbKi), INTENT(IN )         :: ReNum2                            ! the second real number to compare

   LOGICAL                         :: EqualRealNos                     ! the function definition -- returns .true. if the numbers are almost equal

      ! local variables
   REAL(DbKi), PARAMETER           :: Eps = EPSILON(ReNum1)             ! machine precision
   REAL(DbKi), PARAMETER           :: Tol = 100.0_DbKi*Eps / 2.0_DbKi   ! absolute tolerance (ignore the last 2 significant digits)

   REAL(DbKi)                      :: Fraction


      ! make sure we're never trying to get more precision than Tol

   Fraction = MAX( ABS(ReNum1+ReNum2), 1.0_DbKi )



      ! determine if ReNum1 and ReNum2 are approximately equal

   IF ( ABS(ReNum1 - ReNum2) <= Fraction*Tol ) THEN  ! the relative error
      EqualRealNos = .TRUE.
   ELSE
      EqualRealNos = .FALSE.
   ENDIF


   END FUNCTION EqualRealNos    
!==============================================================    
 SUBROUTINE SAVE_TO_TXT_2D(ARRARY_TO_PRINT, ARRY_NAME)   

   IMPLICIT                      NONE
      
      ! Passed Variables:
   REAL(ReKi),DIMENSION(:,:),     INTENT(IN   )  :: ARRARY_TO_PRINT  
   CHARACTER(LEN = *),            INTENT(IN   )  :: ARRY_NAME

      ! Local variables
   INTEGER(IntKi)      :: I_dim
   INTEGER(IntKi)      :: J_dim
   INTEGER(IntKi)      :: I
   INTEGER(IntKi)      :: J   
   INTEGER(IntKi)      :: err
   
   LOGICAL                 :: FILE_EXISTS   
   CHARACTER(LEN = 1024)   :: PATH
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE
   
   
   I_dim   =   SIZE(ARRARY_TO_PRINT,1)
   J_dim   =   SIZE(ARRARY_TO_PRINT,2)


   path='/Treecode_for_debug/'
   
   
   INQUIRE (DIRECTORY = TRIM(path), EXIST = DIR_EXISTS)
      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(path))   
      IF (.NOT. RESULT_CREATE) THEN
          WRITE(6,*) 'Fail to create path '
      END IF         
   ENDIF     
   
      
   
   INQUIRE(FILE = TRIM(path)//ARRY_NAME//'.txt', EXIST=FILE_EXISTS)   ! file_exists will be TRUE if the file
                                                                       !exists and FALSE otherwise
   
   IF (file_exists) THEN
      OPEN (UNIT=10, FILE=TRIM(path)//ARRY_NAME//'.txt', ACTION="write", STATUS="OLD") 
   ELSE
      OPEN (UNIT=10, FILE=TRIM(path)//ARRY_NAME//'.txt', ACTION="write", STATUS="NEW") 
   END IF
   
   
   DO I=1,I_dim
      WRITE(10, '(1600F14.7)')( ARRARY_TO_PRINT(I,J) ,J=1,J_dim)
      !WRITE(10, '(1600d17.10)')( ARRARY_TO_PRINT(I,J) ,J=1,J_dim)
   END DO

   CLOSE(10)

      
      

END SUBROUTINE SAVE_TO_TXT_2D
!================================================
       
END MODULE tree_library