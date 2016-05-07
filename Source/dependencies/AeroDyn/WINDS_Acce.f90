!**********************************************************************************************************************************    
! This module contains acceleration methods for WInDS, including GPU/OpenMP/Treecode/FMM....   
!
! "p%FVM%GPU_Parms%Accelerate" will lead to this part
!
! Last edited: 2014/05/24
!..................................................................................................
    

MODULE WINDS_Accelerate   

   use OMP_LIB      ! OpenMP
   USE NWTC_Library
   USE AeroDyn_Types
    
   IMPLICIT        NONE
 

      ! ..... Public Subroutines ............
   !PUBLIC :: BiotSavart_Acce
   
   PUBLIC :: BiotSavart_TEST      ! For test use ( to compare error and CPU time)
   PUBLIC :: BiotSavart_OpenMP    ! To accelerate Biot-Savart Law by OpenMP    
   

CONTAINS
!==================================================================================================================================
! SUBROUTINE BiotSavart_Acce()                        
    
!END SUBROUTINE BiotSavart_Acce     
!==================================================================================================================================
SUBROUTINE BiotSavart_OpenMP(F1_new, F2_new, Pts_new, GAMMA_new, RC_new, Uind_new,  & 
                           Parameters, OtherState, Discrete_states, ErrStat, ErrMess)
! Refer: http://people.sc.fsu.edu/~jburkardt/f_src/openmp/openmp.html
!        https://software.intel.com/en-us/node/432582
!....................................................................
 
   IMPLICIT                        NONE
   
      ! Passed variables
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: F1_new
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: F2_new 
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: Pts_new
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: GAMMA_new
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: RC_new
   REAL(DbKi), INTENT(INOUT), DIMENSION(:,:)        :: Uind_new   

   
   TYPE(AD_ParameterType),        INTENT(IN   )    :: Parameters       ! Parameters
   TYPE(AD_OtherStateType),       INTENT(IN   )    :: OtherState       ! Other/optimization states
   TYPE(AD_DiscreteStateType),    INTENT(IN   )    :: Discrete_states  ! Discrete states at t
   INTEGER(IntKi),                INTENT(  OUT)    :: ErrStat        ! The error status code
   CHARACTER(*),                  INTENT(  OUT)    :: ErrMess         ! The error message, if an error occurred

      ! Local variables
   INTEGER(IntKi)    :: I, J
   
   REAL(DbKi)     :: GMMA
   REAL(DbKi)     :: RC            
   REAL(DbKi)     :: PX
   REAL(DbKi)     :: PY
   REAL(DbKi)     :: PZ
   REAL(DbKi)     :: CO   
      
   REAL(DbKi)     :: X1
   REAL(DbKi)     :: Y1
   REAL(DbKi)     :: Z1
   REAL(DbKi)     :: X2
   REAL(DbKi)     :: Y2
   REAL(DbKi)     :: Z2
      
   REAL(DbKi)     :: X2X1
   REAL(DbKi)     :: Y2Y1
   REAL(DbKi)     :: Z2Z1
   REAL(DbKi)     :: L 
      
   REAL(DbKi)     :: PXX1
   REAL(DbKi)     :: pyy1
   REAL(DbKi)     :: pzz1
   REAL(DbKi)     :: pxx2
   REAL(DbKi)     :: pyy2
   REAL(DbKi)     :: pzz2
                                 
   REAL(DbKi)     :: R1
   REAL(DbKi)     :: R2
   REAL(DbKi)     :: R1DR2
   REAL(DbKi)     :: R1TR2
   REAL(DbKi)     :: LDR12
   REAL(DbKi)     :: CNU 
   REAL(DbKi)     :: UBAR 
   REAL(DbKi)     :: PRESUMX
   REAL(DbKi)     :: PRESUMY 
   REAL(DbKi)     :: PRESUMZ 
   REAL(DbKi)     :: DEN
   
   INTEGER(IntKi)    :: Len_fila   ! # of filament
   INTEGER(IntKi)    :: Len_pts    ! # of points      
   
   
     ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""  
   
    
   CO = Parameters%FVM%CO
      
   Len_fila   =   SIZE(F1_new , 1)
   Len_pts    =   SIZE(Pts_new, 1) 
   

   ! CALL OMP_SET_NUM_THREADS(Parameters%FVM%OpenMP_Parms%OpenMPCores)
   
   !$OMP PARALLEL DO &
   !$OMP   PRIVATE ( J, GMMA, RC, PX, PY, PZ,  X1, Y1, Z1, X2, Y2, Z2, X2X1, Y2Y1, Z2Z1, L, PXX1, pyy1, pzz1, &
   !$OMP            pxx2, pyy2, pzz2, R1, R2, R1DR2, R1TR2, LDR12, CNU, UBAR, PRESUMX, PRESUMY, PRESUMZ, DEN) 

   ! Loop the points of interest
   DO I = 1, Len_pts 
      
       PX = Pts_new(I, 1)
       PY = Pts_new(I, 2)
       PZ = Pts_new(I, 3)          
   
       PRESUMX = 0
       PRESUMY = 0
       PRESUMZ = 0

                      
       ! Loop the filaments  
       Do J = 1, Len_fila 
           
          X1   =  F1_new(J, 1)
          Y1   =  F1_new(J, 2)
          Z1   =  F1_new(J, 3)
   
          X2   =  F2_new(J, 1)
          Y2   =  F2_new(J, 2)
          Z2   =  F2_new(J, 3)
                             
                                                 
          GMMA =  GAMMA_new(J, 1)
          RC   =  RC_new(J, 1)
                             
          X2X1 =  X2 - X1
          Y2Y1 =  Y2 - Y1
          Z2Z1 =  Z2 - Z1
                             
          L    =  X2X1 * X2X1 + Y2Y1 * Y2Y1 + Z2Z1 * Z2Z1  ! Length of vortex filament (NOTE: L is L^2, as rc is rc^2)
   
          PXX1    =  PX - X1
          PYY1    =  PY - Y1
          PZZ1    =  PZ - Z1
          PXX2    =  PX - X2
          PYY2    =  PY - Y2
          PZZ2    =  PZ - Z2  
                             
          R1      =  SQRT( PXX1 * PXX1 + PYY1 * PYY1 + PZZ1 * PZZ1 )
          R2      =  SQRT( PXX2 * PXX2 + PYY2 * PYY2 + PZZ2 * PZZ2 )
          R1DR2   =  PXX1 * PXX2 + PYY1 * PYY2 + PZZ1 * PZZ2
          R1TR2   =  R1 * R2
          
          IF (Parameters%FVM%ViscFLAG) THEN                                 
              ! Vatistas core model (n=2) ...........................................
              LDR12   =  ( X2X1 * PXX1 + Y2Y1 * PYY1 + Z2Z1 * PZZ1) ** 2                              
              CNU     =  (( R1 * R1 ) - ( LDR12 / L ))                            
              CNU     =  CNU / SQRT(RC ** 2 + CNU  ** 2 ) 
              UBAR    =  CNU * GMMA / (TWOPI * 2) * (R1 + R2) / (R1TR2 * (R1TR2 + R1DR2))
                                 
              ! Check infinity and NAN
              IF ( EqualRealNos( L, 0.0_DbKi ) )                             UBAR = 0.0_DbKi 
              IF ( EqualRealNos( R1 * R1 - LDR12 / L, 0.0_DbKi ) )           UBAR = 0.0_DbKi 
              IF ( EqualRealNos( (R1TR2 * (R1TR2 + R1DR2)) , 0.0_DbKi ) )    UBAR = 0.0_DbKi 
              IF ( ISNAN(UBAR) )                                             UBAR = 0.0_DbKi  
                               
          ELSE ! Smoothing parameter .................................................
             DEN     =  R1TR2 * (R1TR2 + R1DR2) + (Parameters%FVM%DELTA * L)                           
             UBAR    =  (GMMA * (R1 + R2)) / (TWOPI * 2) 
             UBAR    =  UBAR / DEN 
                             
             ! Check infinity and NAN
             IF ( EqualRealNos( DEN, 0.0_DbKi ) )           UBAR = 0.0_DbKi 
             IF ( ISNAN(UBAR) )                             UBAR = 0.0_DbKi  
             
          END IF ! (p%FVM%ViscFLAG)  ...............................................      
          
                
          ! Influence is ignored beyond CO                             
          IF (R1 > CO)          UBAR = 0  
          IF (R2 > CO)          UBAR = 0 
                             
          PRESUMX = PRESUMX + UBAR * (PYY1 * PZZ2 - PZZ1 * PYY2)
          PRESUMY = PRESUMY + UBAR * (PZZ1 * PXX2 - PXX1 * PZZ2)
          PRESUMZ = PRESUMZ + UBAR * (PXX1 * PYY2 - PYY1 * PXX2)
                
       END DO !  J = 1, Len_fila 
       
       Uind_new(I , 1) =  PRESUMX      
       Uind_new(I , 2) =  PRESUMY            
       Uind_new(I , 3) =  PRESUMZ       

       
   END DO !  I = 1, Len_pts 
   
   !$OMP END PARALLEL DO
   
   
END SUBROUTINE BiotSavart_OpenMP
!==================================================================================================================================
SUBROUTINE BiotSavart_TEST(F1_new, F2_new, Pts_new, GAMMA_new, RC_new, Uind_new, Dim_fila, Dim_pts, & 
                           Parameters, OtherState, Discrete_states, ErrStat, ErrMess)
!....................................................................
 
   IMPLICIT                        NONE
   
      ! Passed variables
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: F1_new
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: F2_new 
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: Pts_new
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: GAMMA_new
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: RC_new
   REAL(DbKi), INTENT(INOUT), DIMENSION(:,:)        :: Uind_new   
   INTEGER(IntKi), INTENT(IN   ), DIMENSION(5)      :: Dim_fila   ! DIMENSION of filament
   INTEGER(IntKi), INTENT(IN   ), DIMENSION(5)      :: Dim_pts    ! DIMENSION of filament

   
   TYPE(AD_ParameterType),        INTENT(IN   )    :: Parameters       ! Parameters
   TYPE(AD_OtherStateType),       INTENT(IN   )    :: OtherState       ! Other/optimization states
   TYPE(AD_DiscreteStateType),    INTENT(IN   )    :: Discrete_states  ! Discrete states at t
   INTEGER(IntKi),                INTENT(  OUT)    :: ErrStat        ! The error status code
   CHARACTER(*),                  INTENT(  OUT)    :: ErrMess         ! The error message, if an error occurred

      ! Local variables
   INTEGER(IntKi)    :: I, J
   
   REAL(DbKi)     :: GMMA
   REAL(DbKi)     :: RC            
   REAL(DbKi)     :: PX
   REAL(DbKi)     :: PY
   REAL(DbKi)     :: PZ
   REAL(DbKi)     :: CO   
      
   REAL(DbKi)     :: X1
   REAL(DbKi)     :: Y1
   REAL(DbKi)     :: Z1
   REAL(DbKi)     :: X2
   REAL(DbKi)     :: Y2
   REAL(DbKi)     :: Z2
      
   REAL(DbKi)     :: X2X1
   REAL(DbKi)     :: Y2Y1
   REAL(DbKi)     :: Z2Z1
   REAL(DbKi)     :: L 
      
   REAL(DbKi)     :: PXX1
   REAL(DbKi)     :: pyy1
   REAL(DbKi)     :: pzz1
   REAL(DbKi)     :: pxx2
   REAL(DbKi)     :: pyy2
   REAL(DbKi)     :: pzz2
                                 
   REAL(DbKi)     :: R1
   REAL(DbKi)     :: R2
   REAL(DbKi)     :: R1DR2
   REAL(DbKi)     :: R1TR2
   REAL(DbKi)     :: LDR12
   REAL(DbKi)     :: CNU 
   REAL(DbKi)     :: UBAR 
   REAL(DbKi)     :: PRESUMX
   REAL(DbKi)     :: PRESUMY 
   REAL(DbKi)     :: PRESUMZ 
   REAL(DbKi)     :: DEN
   
   INTEGER(IntKi)    :: Len_fila   ! # of filament
   INTEGER(IntKi)    :: Len_pts    ! # of points      
   
   
    ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""    
      
   
   CO = Parameters%FVM%CO
      
   Len_fila   =   SIZE(F1_new , 1)
   Len_pts    =   SIZE(Pts_new, 1) 
   
   ! check:  Len_fila == Dim_fila(2)  * Dim_fila(4) * Dim_fila(5)
   !         Len_pts  == Dim_pts(2)  * Dim_pts(4) * Dim_pts(5) 
   
   ! Loop the points of interest
   DO I = 1, Len_pts 
      
       PX = Pts_new(I, 1)
       PY = Pts_new(I, 2)
       PZ = Pts_new(I, 3)          
   
       PRESUMX = 0
       PRESUMY = 0
       PRESUMZ = 0
                    
                      
       ! Loop the filaments  
       Do J = 1, Len_fila 
           
          X1   =  F1_new(J, 1)
          Y1   =  F1_new(J, 2)
          Z1   =  F1_new(J, 3)
   
          X2   =  F2_new(J, 1)
          Y2   =  F2_new(J, 2)
          Z2   =  F2_new(J, 3)
                             
                                                 
          GMMA =  GAMMA_new(J, 1)
          RC   =  RC_new(J, 1)
                             
          X2X1 =  X2 - X1
          Y2Y1 =  Y2 - Y1
          Z2Z1 =  Z2 - Z1
                             
          L    =  X2X1 * X2X1 + Y2Y1 * Y2Y1 + Z2Z1 * Z2Z1  ! Length of vortex filament (NOTE: L is L^2, as rc is rc^2)
   
          PXX1    =  PX - X1
          PYY1    =  PY - Y1
          PZZ1    =  PZ - Z1
          PXX2    =  PX - X2
          PYY2    =  PY - Y2
          PZZ2    =  PZ - Z2  
                             
          R1      =  SQRT( PXX1 * PXX1 + PYY1 * PYY1 + PZZ1 * PZZ1 )
          R2      =  SQRT( PXX2 * PXX2 + PYY2 * PYY2 + PZZ2 * PZZ2 )
          R1DR2   =  PXX1 * PXX2 + PYY1 * PYY2 + PZZ1 * PZZ2
          R1TR2   =  R1 * R2
          
          IF (Parameters%FVM%ViscFLAG) THEN                                 
              ! Vatistas core model (n=2) ...........................................
              LDR12   =  ( X2X1 * PXX1 + Y2Y1 * PYY1 + Z2Z1 * PZZ1) ** 2                              
              CNU     =  (( R1 * R1 ) - ( LDR12 / L ))                            
              CNU     =  CNU / SQRT(RC ** 2 + CNU  ** 2 ) 
              UBAR    =  CNU * GMMA / (TWOPI * 2) * (R1 + R2) / (R1TR2 * (R1TR2 + R1DR2))
                                 
              ! Check infinity and NAN
              IF ( EqualRealNos( L, 0.0_DbKi ) )                             UBAR = 0.0_DbKi 
              IF ( EqualRealNos( R1 * R1 - LDR12 / L, 0.0_DbKi ) )           UBAR = 0.0_DbKi 
              IF ( EqualRealNos( (R1TR2 * (R1TR2 + R1DR2)) , 0.0_DbKi ) )    UBAR = 0.0_DbKi 
              IF ( ISNAN(UBAR) )                                             UBAR = 0.0_DbKi  
                               
          ELSE ! Smoothing parameter .................................................
             DEN     =  R1TR2 * (R1TR2 + R1DR2) + (Parameters%FVM%DELTA * L)                           
             UBAR    =  (GMMA * (R1 + R2)) / (TWOPI * 2) 
             UBAR    =  UBAR / DEN 
                             
             ! Check infinity and NAN
             IF ( EqualRealNos( DEN, 0.0_DbKi ) )           UBAR = 0.0_DbKi
             IF ( ISNAN(UBAR) )                             UBAR = 0.0_DbKi 
             
          END IF ! (p%FVM%ViscFLAG)  ...............................................      
          
                
          ! Influence is ignored beyond CO                             
          IF (R1 > CO)          UBAR = 0  
          IF (R2 > CO)          UBAR = 0 
                             
          PRESUMX = PRESUMX + UBAR * (PYY1 * PZZ2 - PZZ1 * PYY2)
          PRESUMY = PRESUMY + UBAR * (PZZ1 * PXX2 - PXX1 * PZZ2)
          PRESUMZ = PRESUMZ + UBAR * (PXX1 * PYY2 - PYY1 * PXX2)
                
       END DO !  J = 1, Len_fila 

       Uind_new(I , 1) =  PRESUMX      
       Uind_new(I , 2) =  PRESUMY            
       Uind_new(I , 3) =  PRESUMZ     
       
   END DO !  I = 1, Len_pts 
   
   

END SUBROUTINE BiotSavart_TEST
                            
!=================================================================================================================================    
              
                           
                           
END MODULE WINDS_Accelerate
!**********************************************************************************************************************************    