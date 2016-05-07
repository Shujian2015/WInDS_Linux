!**********************************************************************************************************************************
! This module includes I/O subroutine for WInDS: read data from input file and output for Paraview.
! 
!==================================================================================================================================
MODULE WINDS_IO

   USE NWTC_Library
   USE AeroDyn_Types
   USE WINDS_Library
   USE IFPORT       ! To use MAKEDIRQQ
   USE OMP_LIB      ! OpenMP
   
   IMPLICIT        NONE

      ! ..... Public Subroutines ............
      
   PUBLIC :: WInDS_ReadInput            ! Read in the input file of WInDS(FVM)
   PUBLIC :: Initialize_paraview_files  ! Initialize paraview files 
   PUBLIC :: CreateVTUembedded          ! Write paraview files 
   PUBLIC :: Close_paraview_files       ! Close paraview files
   PUBLIC :: LB_load_AirfoilData        ! Load the dynamic stall airfoil data
   PUBLIC :: Write_DS_parameters        ! Write dynamic stall airfoil data
   PUBLIC :: WRITE_Treecode             ! Write the speedup/error of Biot-Savart Law
   PUBLIC :: WRITE_KJ                   ! Write the KJ iteration number and Max(D_gamma)
   PUBLIC :: WInDS_WriteSum             ! Write summay file for simulation
   PUBLIC :: WRITE_INTERNAL             ! Write internal varible to a file for Treecode. 
   PUBLIC :: Print_element              ! write the outputs file for elements 

   
CONTAINS   
!==================================================================================================================================    
SUBROUTINE WInDS_ReadInput(InitInp, P, xd, O, ErrStat, ErrMess )
! Read the input of WInDS. Currently it is named "AeroDyn_WInDS.dat" and 
! in the same folder with AeroDyn input file    
!...............................................................   
    
    IMPLICIT                      NONE
   
   ! Passed Variables:
   TYPE(AD_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(AD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   INTEGER,                      INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMess    
    

      ! Local Variables:

   INTEGER                    :: UnIn                 ! Logical unit for the input file.

   CHARACTER(1024)            :: LINE
   CHARACTER(1024)            :: FilePath             ! The path name of the input file (so files listed in it can be defined relative to the main input file location)
   CHARACTER(1024)            :: FilePathName         ! The whole path of WInDS input file + Name
   LOGICAL                    :: TEST1, TEST2
   INTEGER                    :: I
   
   CHARACTER(LEN = 1024)   :: TEMP1, TEMP2, TEMP3  ! debug
   
   
      ! Function definition      
   UnIn = 90    ! Some number for index of the file
   
   !-------------------------------------------------------------------------------------------------
   ! Open the WInDS input file  
   !-------------------------------------------------------------------------------------------------
   CALL GetPath( InitInp%ADFileName, FilePath )                    ! Count everything before (and including) the last "\" or "/". GetPath is defined in NWTC_IO.f90
   ! (sliu: currently the WInDS input file is in the same folder with AeroDyn input file)
   
   p%FVM%WINDS_dir = TRIM(FilePath)     ! record for future use
   
   FilePathName    = TRIM(FilePath)//'AeroDyn_WInDS.dat'
   
   CALL OpenFInpFile(UnIn, TRIM(FilePathName), ErrStat)   ! This routine opens a formatted input file. OpenFInpFile is defined in NWTC_IO.f90
   IF (ErrStat /= ErrID_None ) THEN
      ErrMess = ' Error (in WInDS): Cannot find/open the input file for WInDS under path:"' //TRIM(FilePathName)//'". Please check.'   
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Read the WInDS input file
   !-------------------------------------------------------------------------------------------------
      ! Read in the title
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Title', VarDescr='File title', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF   
   
      ! Read in the description  
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Desription', VarDescr='Desription', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF   
   
   
   !.........................................................................
      ! Read in the "GENERAL SETTING" part 
   !....................................................................    
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part1', VarDescr='"GENERAL SETTING" part', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
      ! Read in the total run time
   CALL ReadVar( UnIn, FilePathName, P%FVM%Total_Time, VarName='TMax', VarDescr='Total run time', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (TMax option). Please check.'
      RETURN   
   END IF 
   
   IF (P%FVM%Total_Time < 0) THEN
      ErrStat = ErrID_Fatal
      ErrMess = ' Error (in WInDS): Expecting postive value in TMax option.'
      RETURN
   END IF
      
   
      ! Read in the number of time steps stored
   CALL ReadVar( UnIn, FilePathName, p%FVM%NTP, VarName='NTP', VarDescr='Number of time steps stored', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%NTP < 3) THEN
      ErrStat = ErrID_Fatal
      ErrMess = ' Error (in WInDS): Expecting 3 or 4 in NTP option.'
      RETURN
   END IF
      
      ! Read in the DT_WINDS / DT_Timestep
   CALL ReadVar( UnIn, FilePathName, p%FVM%DT_RATIO , VarName='DtRatio', VarDescr='DT_WINDS / DT_Timestep', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (DtRatio option). Please check.'
      RETURN   
   END IF    
   
   IF (p%FVM%DT_RATIO < 0) THEN  ! Also check the value in subroutine "WINDS_SetParameters"
      ErrStat = ErrID_Fatal
      ErrMess = ' Error (in WInDS): Expecting postive value in DtRatio option.'
      RETURN
   END IF   
   
         
      ! Read in the whether steady inflow
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='SteadyFlag', VarDescr='Whether steady inflow', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (SteadyFlag option). Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%SteadyFlag =   .TRUE. 
         CASE ('FALSE')
            p%FVM%SteadyFlag =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in SteadyFlag option.'
            RETURN    
         END SELECT

         
   !....................................................................
      ! Read in the "BIOT-SAVART LAW" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part2', VarDescr='"BIOT-SAVART LAW" part ', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   
      ! Read in the Whether freeze the far wake
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='WakeFLAG', VarDescr='Whether to simplify the far wake', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (WakeFLAG option). Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%WakeFLAG =   .TRUE. 
         CASE ('FALSE')
            p%FVM%WakeFLAG =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in WakeFLAG option.'
            RETURN    
         END SELECT   
      
         
      ! Read in the Whether to apply induction to all wake nodes
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='RollFLAG', VarDescr='Whether to apply induction to all wake nodes', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (RollFLAG option). Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%roll =   .TRUE. 
         CASE ('FALSE')
            p%FVM%roll =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in Roll option.'
            RETURN    
         END SELECT   
   
         
       ! Read in RollDist
   CALL ReadVar( UnIn, FilePathName, p%FVM%RollDist, VarName='RollDist', VarDescr='Only work when RollFLGA is true.', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (RollDist option). Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%roll .AND.  p%FVM%RollDist < 1_IntKi) THEN
       ErrStat = ErrID_Fatal
       ErrMess =  ' Error (in WInDS): Expecting positive integer in RollDist option.'
       RETURN
   END IF   
   
       ! Read in WakeLength
   CALL ReadVar( UnIn, FilePathName, p%FVM%WakeDist, VarName='WakeLength', VarDescr='Number of timesteps beyond which the wake is cut off. -1 disables cutoff.', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (WakeLength option). Please check.'
      RETURN   
   END IF 
   
   !IF (p%FVM%roll .AND.  p%FVM%WakeLength < p%FVM%RollDist) THEN
   !    ErrStat = ErrID_Fatal
   !    ErrMess =  ' Error (in WInDS): Expecting WakeLength to be larger than RollDist in WakeLength option.'
   !    RETURN
   !END IF      
         
  
    
      ! Read in the Average wind speed
   CALL ReadVar( UnIn, FilePathName, p%FVM%AveSpeed , VarName='AveSpeed', VarDescr='Average wind speed', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (AveSpeed option). Please check.'
      RETURN   
   END IF    
   
   IF (p%FVM%AveSpeed< 0) THEN 
      ErrStat = ErrID_Fatal
      ErrMess = ' Error (in WInDS): Expecting postive value in AveSpeed option.'
      RETURN
   END IF   
   
   
         
       ! Read in the whether the "frozen" filaments will retain the induced velocity from the most recent timestep. 
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='UindPast', VarDescr='UindPast option', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (UindPast option). Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%UindPast =   .TRUE. 
         CASE ('FALSE')
            p%FVM%UindPast =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in UindPast option.'
            RETURN
      END SELECT      
     
   
   
   
   
       ! Read in the whether Vatistas viscous model
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='ViscFLAG', VarDescr='Vatistas viscous model', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (ViscFLAG option). Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%ViscFLAG =   .TRUE. 
         CASE ('FALSE')
            p%FVM%ViscFLAG =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in ViscFLAG option.'
            RETURN
      END SELECT   
   
   
       ! Read in the Index of Vatistas viscous model
   CALL ReadVar( UnIn, FilePathName, p%FVM%VISC, VarName='ViscModel', VarDescr='Index of Vatistas viscous model', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (ViscModel option). Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%VISC /= 2_IntKi) THEN
       ErrStat = ErrID_Fatal
       ErrMess =  ' Error (in WInDS): Expecting "2" in ViscModel option.'
       RETURN
   END IF   
   
       ! Read in the Smooth parameter
   CALL ReadVar( UnIn, FilePathName, p%FVM%DELTA, VarName='DELTA', VarDescr='Smooth parameter', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%DELTA < 0 .OR.  p%FVM%DELTA > 0.1) THEN
       ErrStat = ErrID_Fatal
       ErrMess =  ' Error (in WInDS): Expecting suitable value in DELTA option.'
       RETURN
   END IF    
   
   !....................................................................
      ! Read in the "NUMERICAL METHODS AND SOLUTION" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part3',                                 &
           VarDescr='"NUMERICAL METHODS AND SOLUTION" part ', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Numerical integration scheme
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='INTEG', VarDescr='Numerical integration scheme', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:3) )
   
   SELECT CASE ( TRIM(Line) )  
         CASE ( 'RK4' )
            p%FVM%INTEG = 'RK4'
         CASE ( 'RK2' )
            p%FVM%INTEG = 'RK2'    
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess =  ' Error (in WInDS): Expecting correct numerical integration scheme in INTEG option.(Only RK4 and RK2 is available now)'
            RETURN
         END SELECT
         
    
       ! Read in the Tolerance value for convergence of numerical methods
   CALL ReadVar( UnIn, FilePathName, p%FVM%TOL, VarName='Tolerance',                            &
            VarDescr='Tolerance value for convergence of numerical methods', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (Tolerance option). Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%TOL > 1.0 .OR. p%FVM%TOL < 0.0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMess = ' Error (in WInDS): Expecting appropriate value in Tolerance option(between 0 and 1.0).'
         RETURN
   END IF      
            
       ! Read in the Distance from wake nodes beyond which influence is negligible
   CALL ReadVar( UnIn, FilePathName, p%FVM%CO , VarName='CO',                          &
            VarDescr='Distance from wake nodes beyond which influence is negligible', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF    
   
       ! Read in the Relaxation value for fixed-point iteration
   CALL ReadVar( UnIn, FilePathName, p%FVM%RELAX , VarName='RELAX',                         &
              VarDescr='Relaxation value for fixed-point iteration', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Maximum number of iterations for Kutta-Joukowski theorem 
   CALL ReadVar( UnIn, FilePathName, p%FVM%MAXITER , VarName='MaxIter',                    &
             VarDescr='Maximum number of iterations for Kutta-Joukowski theorem ', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the whether to use quadratic extrapolation to guess next bound vorticity value
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='ExtrapWake', VarDescr='whether to use quadratic extrapolation to guess next bound vorticity value', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS (ExtrapWake option). Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%extrap_wake =   .TRUE. 
         CASE ('FALSE')
            p%FVM%extrap_wake =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in ExtrapWake option.'
            RETURN
      END SELECT      
   
   
   
   !....................................................................
      ! Read in the "BEM (First timestep in WInDS)" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part4', VarDescr='"BEM (First timestep in WInDS)" part', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the BEM Convergence tolerance
   CALL ReadVar( UnIn, FilePathName, p%FVM%BEM_Parms%TOL , VarName='BEMTOL', VarDescr='Convergence tolerance', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%BEM_Parms%TOL > 1.0 .OR. p%FVM%BEM_Parms%TOL < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMess = ' Error (in WInDS): Expecting appropriate value in BEMTOL option(between 0 and 1.0).'
      RETURN
   END IF         
      
       ! Read in the Maximum number of allowable iterations for BEM
   CALL ReadVar( UnIn, FilePathName,p%FVM%BEM_Parms%MAX_ITER, VarName='MAX_ITER',                   &
        VarDescr='Maximum number of allowable iterations for BEM', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Weighting factor on corrections to balance speed with stability
   CALL ReadVar( UnIn, FilePathName, p%FVM%BEM_Parms%WT , VarName='WEIGHT',                      &
              VarDescr='Weighting factor on corrections to balance speed with stability', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%BEM_Parms%WT > 1.0 .OR. p%FVM%BEM_Parms%WT < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMess =  ' Error (in WInDS): Expecting appropriate value in WEIGHT option(between 0 and 1.0).'
      RETURN
   END IF      
   
      
      
   !....................................................................
      ! Read in the "WIND SHEAR" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part5', VarDescr='Part5', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Wind shear (flag)
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='ShearFLAG', VarDescr='Wind shear (flag)', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%Shear_Parms%ShearFLAG =   .TRUE. 
         CASE ('FALSE')
            p%FVM%Shear_Parms%ShearFLAG =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in ShearFLAG option.'
            RETURN
      END SELECT      
   
       ! Read in the ShearType
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='ShearType', VarDescr='Shear Type', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
      
   CALL Conv2UC( LINE(1:1) )

   SELECT CASE ( TRIM(Line) )
         CASE ('1')
            p%FVM%Shear_Parms%model_type  =   1_IntKi 
         CASE ('2')
            p%FVM%Shear_Parms%model_type  =   2_IntKi
         CASE ('3')
            p%FVM%Shear_Parms%model_type  =   3_IntKi   
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "1", "2" or "3" in ShearType option.'
            RETURN
      END SELECT       
   
      
       ! Read in the Hub height
   CALL ReadVar( UnIn, FilePathName, p%FVM%Shear_Parms%z_ref , VarName='ZRef', VarDescr='Hub height', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   ! ****** sliu: check this

   
   !....................................................................
      ! Read in the "TOWER SHADOW" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part6', VarDescr='Part6', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Include induced velocity from tower effects
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='TWRFLAG',                   &
           VarDescr='Include induced velocity from tower effects', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%Twr_Parms%TWRFLAG =   .TRUE. 
         CASE ('FALSE')
            p%FVM%Twr_Parms%TWRFLAG =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in TWRFLAG option.'
            RETURN
         END SELECT      
         
       ! Read in the Tower UP or DOWN
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='TWRUPDOWN', VarDescr='Tower UP or DOWN', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:1) )

   SELECT CASE ( TRIM(Line) )
         CASE ('1')
            p%FVM%Twr_Parms%TWRUPDOWN  =   1_IntKi 
         CASE ('2')
            p%FVM%Twr_Parms%TWRUPDOWN  =   2_IntKi
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "1" or "2" in TWRUPDOWN option.'
            RETURN
      END SELECT       
   
       ! Read in the Tower METHOD
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='TWRMETHOD', VarDescr='Tower METHOD', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:1) )

   SELECT CASE ( TRIM(Line) )
         CASE ('1')
            p%FVM%Twr_Parms%TWRMETHOD  =   1_IntKi 
         CASE ('2')
            p%FVM%Twr_Parms%TWRMETHOD  =   2_IntKi
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "1" or "2" in TWRMETHOD option.'
            RETURN
         END SELECT        

         
         
   !....................................................................
      ! Read in the "GROUND EFFECTS" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part6', VarDescr='Part6', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Include vortex panels to model ground effects (flag)
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='GroundFLAG',                   &
           VarDescr='Include vortex panels to model ground effects (flag)', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%Ground_Parms%GroundFLAG =   .TRUE. 
         CASE ('FALSE')
            p%FVM%Ground_Parms%GroundFLAG =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in GroundFLAG option.'
            RETURN
         END SELECT     
         
       ! Read in the Method for calculating ground effects
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='GroundMethod',                   &
         VarDescr='Method for calculating ground effects', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('PANEL')
            p%FVM%Ground_Parms%Method =  'PANEL' 
         !CASE ('IMAGE')
         !   p%FVM%Ground_Parms%Method =  'IMAGE'
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            !ErrMess = ' Error (in WInDS): Expecting "PANEL" or "IMAGE" in GroundMethod option.'
            ErrMess = ' Error (in WInDS): Expecting "PANEL" in GroundMethod option. "IMAGE" is not available now'
            RETURN
         END SELECT        
      
       ! Read in the Square of panel quantity
   CALL ReadVar( UnIn, FilePathName, p%FVM%Ground_Parms%Sqrt_Panels ,                   &
         VarName='SqrtPanels', VarDescr='Square of panel quantity', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the GroundXmin
   CALL ReadVar( UnIn, FilePathName, p%FVM%Ground_Parms%Extent(1) , VarName='GroundXmin', VarDescr='GroundXmin', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the GroundXmax
   CALL ReadVar( UnIn, FilePathName, p%FVM%Ground_Parms%Extent(2) , VarName='GroundXmax', VarDescr='GroundXmax', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the GroundYmin
   CALL ReadVar( UnIn, FilePathName, p%FVM%Ground_Parms%Extent(3) , VarName='GroundYmin', VarDescr='GroundYmin', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the GroundYmax
   CALL ReadVar( UnIn, FilePathName, p%FVM%Ground_Parms%Extent(4) , VarName='GroundYmax', VarDescr='GroundYmax', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   !....................................................................
      ! Read in the "DYNAMIC STALL" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part7', VarDescr='DYNAMIC STALL', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
      
       ! Read in the Whether to use LB dynamic stall model
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='DS_Flag', VarDescr='Whether to use LB dynamic stall model', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF             
         
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%DS_Parms%DS_Flag =   .TRUE. 
         CASE ('FALSE')
            p%FVM%DS_Parms%DS_Flag =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in DS_Flag option.'
            RETURN
         END SELECT     
         
   IF (p%FVM%DS_Parms%DS_Flag == .TRUE.  .AND.  P%DStall == .FALSE.) THEN
      ErrStat = ErrID_Fatal
      ErrMess = ' Error (in WInDS): To use dynamic stall model in WInDS, user needs to set "Stall model" in AeroDyn input file to be "BEDDOES".'
      RETURN
   END IF         
   
       ! Read in the Whether to tune LB dynamic stall model
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='RelaxTune', VarDescr='Whether to tune LB dynamic stall model', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF             
         
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%DS_Parms%relax_tune =   .TRUE. 
         CASE ('FALSE')
            p%FVM%DS_Parms%relax_tune =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in RelaxTune option.'
            RETURN
         END SELECT     
            
      !  Read in the time to start ds   
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%start_t , VarName='StartTime', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF          
         
       ! Read in the Whether to load pre-calculated dynamic stall data
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='LoadData', VarDescr='Whether to load pre-calculated dynamic stall data', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF             
         
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )  ! Sliu: should remove this line
         CASE ('FALSE')
            p%FVM%DS_Parms%load_data =   .FALSE. 
         CASE ('TRUE')
            p%FVM%DS_Parms%load_data =   .TRUE.             
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "FALSE" or "TRUE" in LoadData option...'
            RETURN                        
            
         !CASE DEFAULT
         !   ErrStat = ErrID_Fatal
         !   ErrMess = ' Error (in WInDS): Expecting "FALSE" in LoadData option. Only to use pre-calculation here...'
         !   RETURN
         END SELECT     
    
   
      ! Selects File name for pre-calculated dynamic stall data
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%load_file , VarName='LoadData', VarDescr='File name for pre-calculated dynamic stall data', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF              
   

       
      !  Read in indicial coefs, [A1,A2,b1,b2]  
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%indicial(1) , VarName='Indicial1', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF          
   
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%indicial(2) , VarName='Indicial2', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF    
   
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%indicial(3) , VarName='Indicial3', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%indicial(4) , VarName='Indicial4', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF    
   
      !  Read in time constants, [Tp,Tf,Tv,Tvl]  
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%time_const(1) , VarName='TimeConst1', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF          
   
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%time_const(2) , VarName='TimeConst2', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF   
   
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%time_const(3) , VarName='TimeConst3', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF  
   
   CALL ReadVar( UnIn, FilePathName, p%FVM%DS_Parms%time_const(4) , VarName='TimeConst4', VarDescr='time to start ds', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF     
   
     ! Read in the Whether to write dynamic stall data to file
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='WriteData', VarDescr='Whether to write dynamic stall data to file', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF             
         
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%DS_Parms%write_data =   .TRUE. 
         CASE ('FALSE')
            p%FVM%DS_Parms%write_data =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in WriteData option.'
            RETURN
         END SELECT        
   

   !....................................................................
      ! Read in the "ANIMATION" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part8', VarDescr='ANIMATION', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Whether generate animation of wake evolution for Paraview
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='AnimFLAG',                   &
           VarDescr='Whether generate animation of wake evolution for Paraview', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%AnimFLAG  =   .TRUE. 
         CASE ('FALSE')
            p%FVM%AnimFLAG  =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in AnimFLAG option.'
            RETURN
         END SELECT   
         
 
   
   !....................................................................
      ! Read in the "PARALLEL COMPUTING" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part9', VarDescr='Part9', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
       ! Read in the Whether use Treecode/FMM/GPU/OpenMP
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='Accelerate', VarDescr='Whether use Treecode/FMM/GPU/OpenMP', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%OpenMP_Parms%Accelerate =   .TRUE. 
         CASE ('FALSE')
            p%FVM%OpenMP_Parms%Accelerate =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess =  ' Error (in WInDS): Expecting "TRUE" or "FALSE" in Accelerate option.'
            RETURN
         END SELECT     
         
      
         ! Read in number of cores used in OpenMP
   CALL ReadVar( UnIn, FilePathName, p%FVM%OpenMP_Parms%OpenMPCores , VarName='OpenMPCores', VarDescr='Number of OpenMP cores', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF        
        
   IF (p%FVM%OpenMP_Parms%OpenMPCores <=1 ) THEN
      ErrStat = ErrID_Fatal
      ErrMess =  ' Error (in WInDS): Expecting >1 integer in OpenMPCores option.'
      RETURN
   END IF      
         
   CALL OMP_SET_NUM_THREADS(p%FVM%OpenMP_Parms%OpenMPCores)
   
   !....................................................................
      ! Read in the "TREECODE ALGORITHM" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part10', VarDescr='Part10', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF          
         
         
       ! Read in the Whether use Treecode
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='TreeFlag', VarDescr='Whether use Treecode', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%Tree_Parms%TreeFlag =   .TRUE. 
         CASE ('FALSE')
            p%FVM%Tree_Parms%TreeFlag =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess =  ' Error (in WInDS): Expecting "TRUE" or "FALSE" in TreeFlag option.'
            RETURN
         END SELECT    
   
         
       ! Read in the Whether test and record the speedup (flag)
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='Speedup', VarDescr='Whether test and record the speedup', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )

   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%Tree_Parms%Speedup =   .TRUE.  
            p%FVM%Tree_Parms%Freq    =   'WHOL'
         CASE ('LAST')
            p%FVM%Tree_Parms%Speedup =   .TRUE.              
            p%FVM%Tree_Parms%Freq    =   'LAST'
         CASE ('TENS')
            p%FVM%Tree_Parms%Speedup =   .TRUE.              
            p%FVM%Tree_Parms%Freq    =   'TENS'            
         CASE ('HUND')
            p%FVM%Tree_Parms%Speedup =   .TRUE.              
            p%FVM%Tree_Parms%Freq    =   'HUND'
         CASE ('THOU')
            p%FVM%Tree_Parms%Speedup =   .TRUE.              
            p%FVM%Tree_Parms%Freq    =   'THOU' 
         CASE ('FALSE')
            p%FVM%Tree_Parms%Speedup =   .FALSE.            
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess = ' Error (in WInDS): Expecting "TRUE" or "FALSE" in Speedup option.' 
            RETURN
         END SELECT           
         
     
         
         
        ! Read in the parameter for parallel treecode  ! Sliu: should remove this line
   CALL ReadVar( UnIn, FilePathName, p%FVM%Tree_Parms%cores , VarName='Paralle', VarDescr='parallel treecode', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF    
   
   !IF (p%FVM%Tree_Parms%cores < 0 ) THEN
   !   ErrStat = ErrID_Fatal
   !   ErrMess =  ' Error (in WInDS): Expecting appropriate value in Parallel option(postive).'
   !   RETURN
   !END IF              
         
         
       ! Read in the theta (open box angle)
   CALL ReadVar( UnIn, FilePathName, p%FVM%Tree_Parms%theta , VarName='OpenAngle',                      &
              VarDescr='Theta (open box angle)', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%Tree_Parms%theta > 1.0 .OR. p%FVM%Tree_Parms%theta < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMess =  ' Error (in WInDS): Expecting appropriate value in OpenAngle option(between 0 and 1.0).'
      RETURN
   END IF         
   
   
       ! Read in the Taylor expansion order
   CALL ReadVar( UnIn, FilePathName, p%FVM%Tree_Parms%order , VarName='TaylorOrder', VarDescr='TaylorOrder', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF    
   
   IF (p%FVM%Tree_Parms%order > 10 .OR. p%FVM%Tree_Parms%order < 3 ) THEN
      ErrStat = ErrID_Fatal
      ErrMess =  ' Error (in WInDS): Expecting appropriate value in TaylorOrder option(between 3 and 10).'
      RETURN
   END IF     
   
      
   
       ! Read in the Maxparnode: Max particle quantity in one leaf
   CALL ReadVar( UnIn, FilePathName, p%FVM%Tree_Parms%maxparnode , VarName='Maxparnode', VarDescr='Max particle quantity in one leaf', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF    
   
   IF (p%FVM%Tree_Parms%maxparnode < 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMess =  ' Error (in WInDS): Expecting appropriate value in Maxparnode option(postive).'
      RETURN
   END IF        
   
       ! Read in the Dist_tol: Ignored beyond this distance
   CALL ReadVar( UnIn, FilePathName, p%FVM%Tree_Parms%dist_tol , VarName='Dist_tol',                      &
              VarDescr='Ignored beyond this distance', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF    
   
       ! Read in the Delta: Smoothing parameter
   CALL ReadVar( UnIn, FilePathName, p%FVM%Tree_Parms%delta , VarName='Delta',                      &
              VarDescr='Smoothing parameter', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   IF (p%FVM%Tree_Parms%delta > 0.1 .OR. p%FVM%Tree_Parms%delta < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMess =  ' Error (in WInDS): Expecting appropriate value in Delta option(between 0 and 0.1).'
      RETURN
   END IF       
   
         

   !....................................................................
      ! Read in the "OUTPUT" part 
   !....................................................................   
   CALL ReadVar( UnIn, FilePathName, LINE, VarName='Part11', VarDescr='Part11', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF          
   
   
       ! Read in the Whether to write Cl, Cd,  AOA and V_tot of each blade element
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='Element', VarDescr='Whether to write blade elment', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%element_output =   .TRUE. 
         CASE ('FALSE')
            p%FVM%element_output =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess =  ' Error (in WInDS): Expecting "TRUE" or "FALSE" in Element option.'
            RETURN
         END SELECT     
   
   
       ! Read in the Whether to write convergence of Kutta Joukowski subroutine
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='KJCoverg', VarDescr='Whether to write convergence of Kutta Joukowski subroutine', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%KJ_output =   .TRUE. 
         CASE ('FALSE')
            p%FVM%KJ_output =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess =  ' Error (in WInDS): Expecting "TRUE" or "FALSE" in KJCoverg option.'
            RETURN
         END SELECT   
         
         
       ! Read in the Whether to write summary
   CALL ReadVar( UnIn, FilePathName, LINE , VarName='SumPrint', VarDescr='Whether to write summary', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None )  THEN
      ErrMess = ' Error (in WInDS): Error occured when reading input file for WInDS. Please check.'
      RETURN   
   END IF 
   
   CALL Conv2UC( LINE(1:5) )
   
   SELECT CASE ( TRIM(Line) )
         CASE ('TRUE')
            p%FVM%WINDS_Sum =   .TRUE. 
         CASE ('FALSE')
            p%FVM%WINDS_Sum =   .FALSE.
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMess =  ' Error (in WInDS): Expecting "TRUE" or "FALSE" in SumPrint option.'
            RETURN
         END SELECT    
         

      
   !    ! Read in the 
   !CALL ReadVar( UnIn, FilePathName,  , VarName='', VarDescr='', ErrStat=ErrStat)
   !IF ( ErrStat /= ErrID_None ) RETURN    
   
  
   !-------------------------------------------------------------------------------------------------
   ! Close WInDS input file
   !-------------------------------------------------------------------------------------------------
   CLOSE(UnIn)   
           
END SUBROUTINE WInDS_ReadInput
!====================================================================================================
SUBROUTINE Initialize_paraview_files(P, xd,  O, ErrStat, ErrMess )
! Initialize Paraview files
!...................................................................................................    
    IMPLICIT                      NONE
   
   ! Passed Variables:
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Initial discrete states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   INTEGER(IntKi), INTENT(OUT)                  :: ErrStat
   CHARACTER(*), INTENT(OUT)                    :: ErrMess    
    

      ! Local Variables:
   INTEGER(IntKi)          :: IBlade
   CHARACTER(LEN = 1)      :: VisBlade
   CHARACTER(LEN = 1024)   :: PVDfilename 
   CHARACTER(LEN = 1024)   :: RootVisDIR
   CHARACTER(LEN = 1024)   :: WrapperVisDIR   
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE   
   LOGICAL                 :: FILE_EXISTS  
   CHARACTER(LEN = 1024)   :: LINE   
   CHARACTER(LEN = 1024)   :: TEMP
   CHARACTER(LEN = 1024)   :: filename   
   
   ! Create root DIR for all simulations
   RootVisDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs'  ! // p%FVM%AnimRootDIR         
   INQUIRE (DIRECTORY = TRIM(RootVisDIR), EXIST = DIR_EXISTS)
      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(RootVisDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
   ENDIF  
   
   RootVisDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME) ! // p%FVM%AnimRootDIR      
   INQUIRE (DIRECTORY = TRIM(RootVisDIR), EXIST = DIR_EXISTS)
      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(RootVisDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
   ENDIF     
   
   
   ! Create DIR for one case
   WrapperVisDIR  =   TRIM(p%FVM%CURRENT_TIME)//'_VISUAL'
   O%FVM_Other%FullDIR  =   TRIM(RootVisDIR)// '/' // TRIM(WrapperVisDIR)
   
   INQUIRE (DIRECTORY = TRIM(O%FVM_Other%FullDIR), EXIST = DIR_EXISTS)
      
   IF (.NOT. DIR_EXISTS) THEN
      RESULT_CREATE = MAKEDIRQQ (TRIM(O%FVM_Other%FullDIR))  
      
      IF ( .NOT. RESULT_CREATE )  THEN
         ErrStat = ErrID_Fatal 
         ErrMess = ' Error (in WInDS): Failed to create subdirectory for Paraview files. Please check.'
         RETURN   
      END IF       
   ENDIF   
   
   
   DO IBlade = 1, p%FVM%NB
      WRITE (VisBlade, "(I1)") (IBlade)
      PVDfilename    =   '/windstimeseries-blade'// TRIM(VisBlade) //'.pvd'

      INQUIRE(FILE = TRIM(O%FVM_Other%FullDIR)//TRIM(PVDfilename), EXIST=FILE_EXISTS)
      
      IF (.NOT. FILE_EXISTS) THEN
         OPEN (UNIT = IBlade+10, FILE = TRIM(O%FVM_Other%FullDIR)//TRIM(PVDfilename), ACTION="WRITE", STATUS="NEW")  
      ELSE
         OPEN (UNIT = IBlade+10, FILE = TRIM(O%FVM_Other%FullDIR)//TRIM(PVDfilename), ACTION="WRITE", STATUS="REPLACE")          
      END IF
   
      LINE = '<?xml version="1.0"?>'
      WRITE(IBlade+10, "(A)") (TRIM(LINE))
      
      LINE = '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      WRITE(IBlade+10, "(A)") (TRIM(LINE))      
      
      LINE = '<Collection>'
      WRITE(IBlade+10, "(A)") (TRIM(LINE)) 
            
   END DO ! IBlade = 1, p%FVM%NB
     

   


END SUBROUTINE Initialize_paraview_files
!====================================================================================================
SUBROUTINE CreateVTUembedded(IBlade, filename, P, xd, O, ErrStat, ErrMess )
! Write Paraview file at each WInDS timestep
!................................................................................................... 
    IMPLICIT                      NONE
   
   ! Passed Variables:
   INTEGER,                      INTENT(IN   )  :: IBlade
   CHARACTER(LEN = 1024),        INTENT(  OUT)  :: filename   
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Initial discrete states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess    
    

      ! Local Variables:
   CHARACTER(LEN = 1024)        :: TEMP1, TEMP2, TEMP3, TEMP4
   CHARACTER(LEN = 1024)        :: OutputName
   INTEGER(IntKi)               :: IElement
   INTEGER(IntKi)               :: ITimestep
   INTEGER(IntKi)               :: NUM
   CHARACTER(LEN = 1)           :: VisBlade
   CHARACTER(LEN = 1024)        :: timestring
   !INTEGER(IntKi)               :: N_plus
   INTEGER(IntKi)               :: NST
   CHARACTER(LEN = 1024)        :: DATATYPE
   LOGICAL                      :: FILE_EXISTS
   CHARACTER(LEN = 1024)        :: LINE
   INTEGER(IntKi)               :: NTW

   
         ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMess  = ""   
   
   
      
   WRITE (VisBlade, "(I1)") (IBlade)   
   WRITE (timestring, "(I6.6)") (O%WINDS_Timestep)    ! String of the integer With heading zeros
   
      
   filename =  'winds-data-blade-'  // TRIM(VisBlade) // '-time' // timestring(1:6) // '.vtu'
   
   DATATYPE = 'POINTDATA'   
      
   !N_plus   =  O%WINDS_Timestep + 1
   NST      =  p%FVM%NST
   NTW      = O%FVM_Other%NTW 
   
      
   INQUIRE(FILE = TRIM(O%FVM_Other%FullDIR)// '/' //TRIM(filename), EXIST=FILE_EXISTS)
      
   IF (.NOT. FILE_EXISTS) THEN
      OPEN (UNIT = 9, FILE = TRIM(O%FVM_Other%FullDIR)//'/'//TRIM(filename), ACTION="WRITE", STATUS="NEW")     
   ELSE
      OPEN (UNIT = 9, FILE = TRIM(O%FVM_Other%FullDIR)//'/'//TRIM(filename), ACTION="WRITE", STATUS="REPLACE")     
   END IF      
      
   ! HEADER
   LINE = '<?xml version="1.0"?>'
   WRITE(9, "(A)") (TRIM(LINE))        
      
   LINE = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
   WRITE(9, "(A)") (TRIM(LINE))   

   LINE = '  <UnstructuredGrid>'
   WRITE(9, "(A)") (TRIM(LINE))     
      
   TEMP1 = Num2LStr( (NTW+1) * NST ) 
   TEMP2 = Num2LStr( NTW * (NST-1)) 
      
   LINE = '    <Piece NumberOfPoints="'// TRIM(TEMP1) //'" NumberOfCells="'// TRIM(TEMP2) //'">'
   WRITE(9, "(A)") (TRIM(LINE))     
   
   
   ! POINTS
   LINE = '      <Points>'
   WRITE(9, "(A)") (TRIM(LINE)) 
   
   LINE = '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
   WRITE(9, "(A)") (TRIM(LINE))    
   
 
   DO ITimestep = 1, ntw+1  
      DO IElement = 1, NST
         WRITE(TEMP1, "(F13.8)") (O%FVM_Other%WAKE_DOMAIN(1, ITimestep, 1, IElement, IBlade))  
         WRITE(TEMP2, "(F13.8)") (O%FVM_Other%WAKE_DOMAIN(2, ITimestep, 1, IElement, IBlade)) 
         WRITE(TEMP3, "(F13.8)") (O%FVM_Other%WAKE_DOMAIN(3, ITimestep, 1, IElement, IBlade))  
   
         LINE = '          ' // TRIM(TEMP1) // ' ' // TRIM(TEMP2) // ' ' // TRIM(TEMP3) 
         WRITE(9, "(A)") (TRIM(LINE))      
      END DO
   END DO
   
   LINE = '        </DataArray>'
   WRITE(9, "(A)") (TRIM(LINE))       
   
   LINE = '      </Points>'
   WRITE(9, "(A)") (TRIM(LINE))   
   
   ! CELLS
   LINE = '      <Cells>'
   WRITE(9, "(A)") (TRIM(LINE))   
   
   LINE = '        <DataArray type="Int32" Name="connectivity" format="ascii">'
   WRITE(9, "(A)") (TRIM(LINE))     
   
   DO ITimestep = 0, ntw -1
      DO IElement = 0, NST -2
         NUM =  IElement + (ITimestep * NST)
          
         TEMP1= Num2LStr(NUM)  
         TEMP2= Num2LStr(NUM + 1) 
         TEMP3= Num2LStr(NUM + (NST + 1))  
         TEMP4= Num2LStr(NUM + NST)  
   
         LINE = '           ' // TRIM(TEMP1) // ' ' // TRIM(TEMP2) // ' ' // TRIM(TEMP3) // ' ' // TRIM(TEMP4) 
         WRITE(9, "(A)") (TRIM(LINE))      
      END DO
   END DO   
   
   LINE = '        </DataArray>'
   WRITE(9, "(A)") (TRIM(LINE))       
   
   LINE = '        <DataArray type="Int32" Name="offsets" format="ascii">'
   WRITE(9, "(A)") (TRIM(LINE))       

   DO ITimestep = 0, ntw-1
      DO IElement = 1, NST - 1
         NUM =  IElement + ITimestep * (NST-1)
          
         TEMP1 = Num2LStr(NUM * 4)  
   
         LINE = '           ' // TRIM(TEMP1)
         WRITE(9, "(A)") (TRIM(LINE))      
      END DO
   END DO      
   
   LINE = '        </DataArray>'
   WRITE(9, "(A)") (TRIM(LINE))    
      
   LINE = '        <DataArray type="UInt8" Name="types" format="ascii">'
   WRITE(9, "(A)") (TRIM(LINE))       

   DO ITimestep = 1, ntw
      DO IElement = 1, NST - 1         
         WRITE(TEMP1, "(I1)") ( 9 )     
         LINE = '           ' // TRIM(TEMP1)
         WRITE(9, "(A)") (TRIM(LINE))      
      END DO
   END DO      
      
   LINE = '        </DataArray>'
   WRITE(9, "(A)") (TRIM(LINE))       
   
   LINE = '      </Cells>'
   WRITE(9, "(A)") (TRIM(LINE))   
   
  
   ! POINT DATA
   IF (DATATYPE == 'POINTDATA') THEN
       
      OutputName = 'velocity'       ! def of std value
      LINE = '      <PointData Vectors="' // TRIM(OutputName) //'">'   
      WRITE(9, "(A)") (TRIM(LINE))    
      LINE = '        <DataArray type="Float64" Name="' // TRIM(OutputName) // '" NumberOfComponents="3" format="ascii">'   
      WRITE(9, "(A)") (TRIM(LINE))   
      
      DO ITimestep = 1, ntw+1      
         DO IElement = 1, NST             
         WRITE(TEMP1, "(E10.3)") (O%FVM_Other%VEL_DOMAIN(1, ITimestep, 1, IElement, IBlade))  
         WRITE(TEMP2, "(E10.3)") (O%FVM_Other%VEL_DOMAIN(2, ITimestep, 1, IElement, IBlade)) 
         WRITE(TEMP3, "(E10.3)") (O%FVM_Other%VEL_DOMAIN(3, ITimestep, 1, IElement, IBlade))  
   
         LINE = '          ' // TRIM(TEMP1) // ' ' // TRIM(TEMP2) // ' ' // TRIM(TEMP3) 
         WRITE(9, "(A)") (TRIM(LINE))      

         END DO
      END DO   
      
      LINE = '        </DataArray>'
      WRITE(9, "(A)") (TRIM(LINE))            
         
      LINE = '      </PointData>'
      WRITE(9, "(A)") (TRIM(LINE))  
            
   END IF ! (DATATYPE == 'POINTDATA')   
   
   ! FOOTER
   
   
   LINE = '    </Piece>'
   WRITE(9, "(A)") (TRIM(LINE))      
   
   LINE = '  </UnstructuredGrid>'
   WRITE(9, "(A)") (TRIM(LINE))      
   
   LINE = '</VTKFile>'
   WRITE(9, "(A)") (TRIM(LINE))    
   
   CLOSE ( 9 )
   
   

END SUBROUTINE CreateVTUembedded 
!====================================================================================================
SUBROUTINE Close_paraview_files(P, xd, O, ErrStat, ErrMess )
!...................................................................................................    
    IMPLICIT                      NONE
   
   ! Passed Variables:
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Initial discrete states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   INTEGER(IntKi), INTENT(OUT)                  :: ErrStat
   CHARACTER(*), INTENT(OUT)                    :: ErrMess    
    
   
      ! Local Variables:
   INTEGER(IntKi)          :: IBlade
   CHARACTER(LEN = 1024)   :: LINE     
   
         ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMess  = ""   
   
   
   DO IBlade = 1, p%FVM%NB
       
      LINE = '</Collection>'
      WRITE(IBlade+10, "(A)") (TRIM(LINE)) 
      
      LINE = '</VTKFile>'
      WRITE(IBlade+10, "(A)") (TRIM(LINE))  
            
      CLOSE(IBlade+10)      

   END DO
     
      

END SUBROUTINE Close_paraview_files
!====================================================================================================
SUBROUTINE LB_load_AirfoilData(p, O, xd, ErrStat, ErrMess)
!...................................................................................................  

  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
   TYPE(AD_DiscreteStateType),    INTENT(IN   )   :: xd             ! Discrete states at t   
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess
   
 
   ! Local Variables:
   INTEGER(IntKi)      :: NFOILID
   INTEGER(IntKi)      :: Check_name   
   INTEGER             :: UnIn                 ! Logical unit for the input file.
   CHARACTER(1024)     :: LINE
   CHARACTER(1024)     :: FilePathName         ! The whole path of WInDS input file + Name
   
   !LOGICAL             :: TEST1, TEST2
   !INTEGER             :: I
   
   
      ! Function definition      
   UnIn = 200    ! Some number for index of the file
   
   !-------------------------------------------------------------------------------------------------
   ! Open the WInDS input file  
   !-------------------------------------------------------------------------------------------------
   FilePathName    = TRIM(p%FVM%WINDS_dir)// TRIM(p%FVM%DS_Parms%load_file)
   
   CALL OpenFInpFile(UnIn, TRIM(FilePathName), ErrStat)   ! This routine opens a formatted input file. OpenFInpFile is defined in NWTC_IO.f90
   IF (ErrStat /= ErrID_None ) THEN
      ErrMess = ' Error (in WInDS): Cannot find/open the dynamic stall data file for WInDS under path:"' //TRIM(FilePathName)//   &
                 '". Please check. Otherwise, please change "LoadData" option to be True.'   
      RETURN
   END IF
   
   DO NFOILID = 1, P%AirFoil%NumFoil
        ! Read in the airfoil name
        CALL ReadVar( UnIn, FilePathName, LINE, VarName='Title', VarDescr='File title', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
           ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS. Please check.'
           RETURN   
        END IF          
      
        ! Read in A1
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%A1(NFOILID), VarName='A1', VarDescr='A1', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (A1). Please check.'
            RETURN   
        END IF       
      
        ! Read in A2
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%A2(NFOILID), VarName='A2', VarDescr='A2', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (A2). Please check.'
            RETURN   
        END IF            
      
        ! Read in b1
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%b1(NFOILID), VarName='b1', VarDescr='b1', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (b1). Please check.'
            RETURN   
        END IF      
        
        ! Read in  b2
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%b2(NFOILID), VarName='b2', VarDescr='b2', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (b2). Please check.'
            RETURN   
        END IF            
      
        ! Read in Tp
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%Tp(NFOILID), VarName='Tp', VarDescr='Tp', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (Tp). Please check.'
            RETURN   
        END IF    
      
        ! Read in Tf
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%Tf(NFOILID), VarName='Tf', VarDescr='Tf', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (Tf). Please check.'
            RETURN   
        END IF            
      
        ! Read in Tv
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%Tv(NFOILID), VarName='Tv', VarDescr='Tv', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (Tv). Please check.'
            RETURN   
        END IF      
        
        ! Read in Tvl
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%Tvl(NFOILID), VarName='Tvl', VarDescr='Tvl', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (Tvl). Please check.'
            RETURN   
        END IF            
      
        ! Read in Cn_alpha
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID), VarName='Cn_alpha', VarDescr='Cn_alpha', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (Cn_alpha). Please check.'
            RETURN   
        END IF      
                       
        ! Read in Alpha_0
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%Alpha_0(NFOILID), VarName='Alpha_0', VarDescr='Alpha_0', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (Alpha_0). Please check.'
            RETURN   
        END IF            
      
        ! Read in C_d_0
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%C_d_0(NFOILID), VarName='C_d_0', VarDescr='C_d_0', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (C_d_0). Please check.'
            RETURN   
        END IF      
        
        ! Read in C_m_0
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%C_m_0(NFOILID), VarName='C_m_0', VarDescr='C_m_0', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (C_m_0). Please check.'
            RETURN   
        END IF            
      
        ! Read in C_n_0
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%C_n_0(NFOILID), VarName='C_n_0', VarDescr='C_n_0', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (C_n_0). Please check.'
            RETURN   
        END IF    
      
        ! Read in alpha_1
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%alpha_1(NFOILID), VarName='alpha_1', VarDescr='alpha_1', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (alpha_1). Please check.'
            RETURN   
        END IF            
      
        ! Read in  alpha_2
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%alpha_2(NFOILID), VarName='alpha_2', VarDescr='alpha_2', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (alpha_2). Please check.'
            RETURN   
        END IF      
        
        ! Read in a(1)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%a(NFOILID, 1), VarName='a(1)', VarDescr='a(1)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (a(1)). Please check.'
            RETURN   
        END IF            
        
        ! Read in a(2)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%a(NFOILID, 2), VarName='a(2)', VarDescr='a(2)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (a(2)). Please check.'
            RETURN   
        END IF     
        
        ! Read in a(3)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%a(NFOILID, 3), VarName='a(3)', VarDescr='a(3)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (a(3)). Please check.'
            RETURN   
        END IF             

        ! Read in S(1)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%S(NFOILID, 1), VarName='S(1)', VarDescr='S(1)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (S(1)). Please check.'
            RETURN   
        END IF            
        
        ! Read in S(2)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%S(NFOILID, 2), VarName='S(2)', VarDescr='S(2)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (S(2)). Please check.'
            RETURN   
        END IF     
        
        ! Read in S(3)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%S(NFOILID, 3), VarName='S(3)', VarDescr='S(3)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (S(3)). Please check.'
            RETURN   
        END IF                     
        
        ! Read in c(1)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%c(NFOILID, 1), VarName='c(1)', VarDescr='c(1)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (c(1)). Please check.'
            RETURN   
        END IF            
        
        ! Read in c(2)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%c(NFOILID, 2), VarName='c(2)', VarDescr='c(2)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (c(2)). Please check.'
            RETURN   
        END IF     
        
        ! Read in c(3)
        CALL ReadVar( UnIn, FilePathName,  O%FVM_Other%airfoils_LB%c(NFOILID, 3), VarName='c(3)', VarDescr='c(3)', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (c(3)). Please check.'
            RETURN   
        END IF     
        
        ! Read in recovery_factor
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%recovery_factor(NFOILID), VarName='recovery_factor', VarDescr='recovery_factor', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (recovery_factor). Please check.'
            RETURN   
        END IF          
      
      
        ! Read in C_n_2
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%C_n_2(NFOILID), VarName='C_n_2', VarDescr='C_n_2', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (C_n_2). Please check.'
            RETURN   
        END IF      
        
        ! Read in C_n_1
        CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%C_n_1(NFOILID), VarName='C_n_1', VarDescr='C_n_1', ErrStat=ErrStat)
        IF ( ErrStat /= ErrID_None )  THEN
            ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (C_n_1). Please check.'
            RETURN   
        END IF            
        
        !! Read in ca_K1
        !CALL ReadVar( UnIn, FilePathName, O%FVM_Other%airfoils_LB%C_n_1(NFOILID), VarName='ca_K1', VarDescr='ca_K1', ErrStat=ErrStat)
        !IF ( ErrStat /= ErrID_None )  THEN
        !    ErrMess = ' Error (in WInDS): Error occured when reading Dynamic Stall data file for WInDS (ca_K1). Please check.'
        !    RETURN   
        !END IF      

   END DO !  I
   
   
   CLOSE(UnIn)
   

END SUBROUTINE LB_load_AirfoilData
!====================================================================================================
SUBROUTINE Write_DS_parameters(p, O, ErrStat, ErrMess)
!...................................................................................................    

  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(ad_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess
    
    ! Local Variables:   
   CHARACTER(LEN = 1024)   :: SumDIR
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE   
   CHARACTER(LEN = 1024)   :: Sumfilename
   LOGICAL                 :: FILE_EXISTS    
   CHARACTER(LEN = 1024)   :: LINE   
   CHARACTER(LEN = 1024)   :: TEMP1
   INTEGER(IntKi)          :: NFOILID
   INTEGER(IntKi)          :: Check_name
      

    ! Create root DIR for all simulations
   SumDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs'   
   
   INQUIRE (DIRECTORY = TRIM(SumDIR), EXIST = DIR_EXISTS)
      
   IF (.NOT. DIR_EXISTS) THEN
       RESULT_CREATE = MAKEDIRQQ (TRIM(SumDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          IF ( .NOT. RESULT_CREATE )  THEN       
              ErrStat = ErrID_Fatal          
              ErrMess = ' Error (in WInDS): Failed to create subdirectory for Write_DS_parameters. Please check.'         
              RETURN         
          END IF         
   ENDIF     
      
   Sumfilename    =  '/DynamicStall(output).dat'

   INQUIRE(FILE = TRIM(SumDIR)//TRIM(Sumfilename), EXIST=FILE_EXISTS)
      
   IF (.NOT. FILE_EXISTS) THEN
      OPEN (UNIT = 30, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="NEW")
   ELSE
      OPEN (UNIT = 30, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="REPLACE") 
   END IF       

       
   
   DO NFOILID = 1, P%AirFoil%NumFoil
       Check_name =  INDEX((p%AirFoil%FOILNM(NFOILID)), "AeroData") 
       IF (Check_name >=1) THEN  
          LINE =  'AirFoil Name: ' // TRIM(p%AirFoil%FOILNM(NFOILID)(Check_name+9:)) // '-------------------'  ! Just leave the airfoil file name, delete the path
       ELSE
          LINE =  'AirFoil Name: ' // TRIM(p%AirFoil%FOILNM(NFOILID)) // '-----------------------------------' ! Airfoil file name with the path
       END IF
       WRITE(30, "(A)") (TRIM(LINE))
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID) ) 
       LINE = TRIM(TEMP1)  // ' - Cn_alpha'      
       WRITE(30, "(A)") (TRIM(LINE))          

       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%C_d_0(NFOILID) ) 
       LINE = TRIM(TEMP1)  //  ' - C_d_0'      
       WRITE(30, "(A)") (TRIM(LINE))   
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%C_m_0(NFOILID) ) 
       LINE = TRIM(TEMP1)   //   ' - C_m_0'    
       WRITE(30, "(A)") (TRIM(LINE))          
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%C_n_0(NFOILID) ) 
       LINE = TRIM(TEMP1)  //  ' - C_n_0: '    
       WRITE(30, "(A)") (TRIM(LINE))        
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%Alpha_0(NFOILID) * 180/pi) 
       LINE = TRIM(TEMP1)  //  ' - Alpha_0 (deg)'     
       WRITE(30, "(A)") (TRIM(LINE))         
       
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%aoa_KHn(NFOILID) ) 
       LINE = TRIM(TEMP1)   //   ' - aoa_KHn'    
       WRITE(30, "(A)") (TRIM(LINE))          
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%aoa_KHp(NFOILID) ) 
       LINE = TRIM(TEMP1)  //  ' - aoa_KHp: '    
       WRITE(30, "(A)") (TRIM(LINE))    
       
       
       
       
       

       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%alpha_1(NFOILID) ) 
       LINE =  TRIM(TEMP1)  //  ' - alpha_1 (radius)'     
       WRITE(30, "(A)") (TRIM(LINE))          
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%alpha_2(NFOILID) ) 
       LINE = TRIM(TEMP1)  //   ' - alpha_2 (radius)'     
       WRITE(30, "(A)") (TRIM(LINE))  

       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%C_n_2(NFOILID) ) 
       LINE = TRIM(TEMP1) //  ' - C_n_2'       
       WRITE(30, "(A)") (TRIM(LINE))    
       
       WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%C_n_1(NFOILID) ) 
       LINE = TRIM(TEMP1)  //   ' - C_n_1'     
       WRITE(30, "(A)") (TRIM(LINE))           

       
   END DO !  NFOILID          
       
   CLOSE(30)  
       
   
END SUBROUTINE Write_DS_parameters
!====================================================================================================
SUBROUTINE WRITE_Treecode(Len_fila, Len_pts, time_original, time_accelerated, speedup_ratio, error, CALLER, p)
!...................................................................................................    
    IMPLICIT                      NONE
    
   ! Passed Variables:
   INTEGER(IntKi),       INTENT(IN   )  :: Len_fila         ! # of source filaments
   INTEGER(IntKi),       INTENT(IN   )  :: Len_pts          ! # of receiver points
   REAL(DbKi),           INTENT(IN   )  :: time_original    ! Time of Direct computation
   REAL(DbKi),           INTENT(IN   )  :: time_accelerated ! Time of high speed computation
   REAL(DbKi),           INTENT(IN   )  :: speedup_ratio    ! time_original / time_accelerated
   REAL(DbKi),           INTENT(IN   )  :: error            ! L1 norm error
   CHARACTER(*),         INTENT(IN   )  :: CALLER           ! Who calls this subroutine 
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
    
    ! Local Variables:   
   CHARACTER(LEN = 1024)   :: SumDIR
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE   
   CHARACTER(LEN = 1024)   :: Sumfilename
   LOGICAL                 :: FILE_EXISTS    
   CHARACTER(LEN = 1024)   :: LINE   
   CHARACTER(LEN = 1024)   :: TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6  
   
   
   
   
   SELECT CASE  (CALLER)    
       
   CASE ('START')
         ! Create root DIR for all simulations
       SumDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs'      
       INQUIRE (DIRECTORY = TRIM(SumDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(SumDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          !IF ( .NOT. RESULT_CREATE )  THEN        
          !    ErrStat = ErrID_Fatal          
          !    ErrMess = ' Error (in WInDS): Failed to create subdirectory for WRITE_Treecode. Please check.'         
          !    RETURN         
          !END IF         
       ENDIF            
       
       SumDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME)  !  'WInDS_Summary'      
       INQUIRE (DIRECTORY = TRIM(SumDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(SumDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
          !
          !IF ( .NOT. RESULT_CREATE )  THEN       
          !    ErrStat = ErrID_Fatal          
          !    ErrMess = ' Error (in WInDS): Failed to create subdirectory for WRITE_Treecode. Please check.'         
          !    RETURN         
          !END IF         
       ENDIF     
      
       Sumfilename    = '/' // TRIM(p%FVM%CURRENT_TIME) // '_Speedup.dat'

       INQUIRE(FILE = TRIM(SumDIR)//TRIM(Sumfilename), EXIST=FILE_EXISTS)
      
       IF (.NOT. FILE_EXISTS) THEN
          OPEN (UNIT = 20, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="NEW")   
       ELSE
          OPEN (UNIT = 20, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="REPLACE")               
       END IF       
       
       
       WRITE(20, "(A)") ('Settings for Treecode Algorithm')
       WRITE(20, "(A)") ('  Cores:       ' // TRIM( Num2LStr( p%FVM%OpenMP_Parms%OpenMPCores )))
       WRITE(20, "(A)") ('  OpenAngle:   ' // TRIM( Num2LStr( p%FVM%Tree_Parms%theta )))
       WRITE(20, "(A)") ('  TaylorOrder: ' // TRIM( Num2LStr( p%FVM%Tree_Parms%order )))
       WRITE(20, "(A)") ('  Maxparnode:  ' // TRIM( Num2LStr( p%FVM%Tree_Parms%maxparnode )))
       WRITE(20, "(A)") ('  Dist_tol:    ' // TRIM( Num2LStr( p%FVM%Tree_Parms%dist_tol )))
       WRITE(20, "(A)") ('  Smoothing:   ' // TRIM( Num2LStr( p%FVM%Tree_Parms%delta )))
       SELECT CASE (p%FVM%Tree_Parms%Freq) 
          CASE ('WHOL') 
              WRITE(20, "(A)") ('  Record:      ' // TRIM( 'All timesteps' ))
          CASE ('LAST') 
              WRITE(20, "(A)") ('  Record:      ' // TRIM( 'Only last timestep' ))
          CASE ('TENS')   
              WRITE(20, "(A)") ('  Record:      ' // TRIM( 'Every 10 timestep' ))
          CASE ('HUND') 
              WRITE(20, "(A)") ('  Record:      ' // TRIM( 'Every 100 timestep' ))
          CASE ('THOU') 
              WRITE(20, "(A)") ('  Record:      ' // TRIM( 'Every 100 timestep' ))
          END SELECT        
          
       WRITE(20, "(A)") ('.............................................................................')
       WRITE(20, "(A)") ('# of filaments.   # of Points.   Direct time.   '// &
                         'Treecode time.   Speedup ratio.    L2 norm error')
       
   CASE ('NONE') 
       
       TEMP1 = Num2LStr( Len_fila ) 
       TEMP2 = Num2LStr( Len_pts ) 
       WRITE(TEMP3, "(F13.8)") (time_original)  
       WRITE(TEMP4, "(F13.8)") (time_accelerated)  
       WRITE(TEMP5, "(F13.8)") (speedup_ratio) 
       WRITE(TEMP6, "(F13.8)") (error)  
   
       LINE = TRIM(TEMP1) // '              ' // TRIM(TEMP2) // '        '  // TRIM(TEMP3) // '  ' // & 
              TRIM(TEMP4) // '    ' //  TRIM(TEMP5)// '     ' // TRIM(TEMP6)
       
       WRITE(20, "(A)") (TRIM(LINE))         
       
              
   CASE ('END')
       
        CLOSE(20)  
       
   END SELECT    
   
END SUBROUTINE WRITE_Treecode
!====================================================================================================
SUBROUTINE WInDS_WriteSum(P, xd, O, ErrStat, ErrMess )
!
! Modified from SUBROUTINE RunTimes in FAST_Subs.f90
!...................................................................................................    
    IMPLICIT                      NONE
   
   ! Passed Variables:
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Initial discrete states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   !CHARACTER(*),                 INTENT(IN   )  :: CALLER       ! Who calls this subroutine 
   INTEGER(IntKi), INTENT(OUT)                  :: ErrStat
   CHARACTER(*), INTENT(OUT)                    :: ErrMess    
    
    ! Local Variables:   
   CHARACTER(LEN = 1024)   :: SumDIR
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE   
   CHARACTER(LEN = 1024)   :: Sumfilename
   LOGICAL                 :: FILE_EXISTS    
   CHARACTER(LEN = 1024)   :: LINE   
   CHARACTER(LEN = 1024)   :: TEMP1, TEMP2, TEMP3! , TEMP4, TEMP5, TEMP6  
   
   
   INTEGER(IntKi)               :: StrtTime (8)                                    ! Start time of simulation (including initialization)
   INTEGER(IntKi)               :: SimStrtTime (8)                                 ! Start time of simulation (after initialization)
   REAL(DbKi)                   :: UsrTime1                                        ! User CPU time for simulation initialization.
   REAL(DbKi)                   :: UsrTime2                                        ! User CPU time for simulation (without intialization)
   REAL(DbKi)                   :: ZTime                                           ! The final simulation time (not necessarially TMax)
   !REAL,OPTIONAL, INTENT(OUT)   :: UsrTime_out                                     ! User CPU time for entire run - optional value returned to calling routine

      ! Local variables

   REAL(DbKi)                   :: ClckTime                                        ! Elapsed clock time for the entire run.
   REAL(DbKi)                   :: ClckTimeSim                                     ! Elapsed clock time for the simulation phase of the run.
   REAL(DbKi)                   :: Factor                                          ! Ratio of seconds to a specified time period.
   REAL(DbKi)                   :: TRatio                                          ! Ratio of simulation time to elapsed clock time.
   REAL(DbKi), PARAMETER        :: SecPerDay = 24*60*60.0_ReKi                     ! Number of seconds per day

   REAL(DbKi)                   :: UsrTime                                         ! User CPU time for entire run.
   REAL(DbKi)                   :: UsrTimeSim                                      ! User CPU time for simulation (not including initialization).
   INTEGER(IntKi)               :: EndTimes (8)                                    ! An array holding the ending clock time of the simulation.

   CHARACTER( 8)                :: TimePer
   CHARACTER(MaxWrScrLen)       :: BlankLine  
             
 
   
       !....................................................................... 
         ! Create root DIR for all simulations
       SumDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs'   
       INQUIRE (DIRECTORY = TRIM(SumDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(SumDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          IF ( .NOT. RESULT_CREATE )  THEN        
              ErrStat = ErrID_Fatal          
              ErrMess = ' Error (in WInDS): Failed to create subdirectory for WInDS_WriteSum. Please check.'         
              RETURN         
          END IF         
       ENDIF     
              
       
       SumDIR    =  TRIM(p%FVM%WINDS_dir) //  'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME)       
       INQUIRE (DIRECTORY = TRIM(SumDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(SumDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          IF ( .NOT. RESULT_CREATE )  THEN       
              ErrStat = ErrID_Fatal          
              ErrMess = ' Error (in WInDS): Failed to create subdirectory for # WInDS_WriteSum #. Please check.'         
              RETURN         
          END IF         
       ENDIF     
      
       Sumfilename    = '/' // TRIM(p%FVM%CURRENT_TIME) // '_WInDS_sum.dat'
   
       INQUIRE(FILE = TRIM(SumDIR)//TRIM(Sumfilename), EXIST=FILE_EXISTS)
             
       IF (.NOT. FILE_EXISTS) THEN
          OPEN (UNIT = 40, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="NEW")     
       ELSE
          OPEN (UNIT = 40, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="REPLACE")           
       END IF      

       
       
   !.......................................................................                      
   WRITE(40, "(A)") ('This file is summary information for WInDS (free vortex aerodynamic method)')
   WRITE(40, "(A)") ('...........................................................................')

   
   StrtTime = O%FVM_Other%StrtTime
   SimStrtTime = O%FVM_Other%SimStrtTime
   UsrTime1 =  O%FVM_Other%UsrTime1             
   UsrTime2 =  O%FVM_Other%UsrTime2                           
   ZTime    =  O%FVM_Other%ZTime           
       
       
      ! Get the end times to compare with start times.

   CALL DATE_AND_TIME ( VALUES=EndTimes )
   CALL CPU_TIME ( UsrTime )
   UsrTime = MAX( 0.0_DbKi, UsrTime )  ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   

   ! Calculate the elapsed wall-clock time in seconds.

   ClckTime     = GetClockTime(StrtTime,      EndTimes)
  !ClckTimeInit = GetClockTime(StrtTime,   SimStrtTime)
   ClckTimeSim  = GetClockTime(SimStrtTime,   EndTimes)

      ! Calculate CPU times.

   UsrTime    = MAX( 0.0_DbKi, UsrTime - UsrTime1 )
   UsrTimeSim = MAX( 0.0_DbKi, UsrTime - UsrTime2 )


   IF ( .NOT. EqualRealNos( UsrTimeSim, 0.0_DbKi ) .AND. ZTime > 0.0_DbKi )  THEN

      TRatio = REAL(ZTime) / UsrTimeSim

      IF     ( UsrTime > SecPerDay )  THEN
         Factor = 1.0/SecPerDay
         TimePer = ' days'
      ELSEIF ( UsrTime >  3600.0 )  THEN
         Factor = 1.0/3600.0
         TimePer = ' hours'
      ELSEIF ( UsrTime >    60.0 )  THEN
         Factor = 1.0/60.0
         TimePer = ' minutes'
      ELSE
         Factor = 1.0
         TimePer = ' seconds'
      ENDIF

             
      WRITE(40, "(A)") ('Timing:')    
      WRITE(40, "(A)")  ('')
      WRITE(40, "(A)")  ( ' Total Real Time:       '//TRIM( Num2LStr( Factor*ClckTime      ) )//TRIM( TimePer ) )
      WRITE(40, "(A)")  ( ' Total CPU Time:        '//TRIM( Num2LStr( Factor*UsrTime       ) )//TRIM( TimePer ) )
!     WRITE(40, "(A)")  ( ' ')
!     WRITE(40, "(A)")  ( ' Simulation Real Time:  '//TRIM( Num2LStr( Factor*ClckTimeSim   ) )//TRIM( TimePer ) )
      WRITE(40, "(A)")  ( ' Simulation CPU Time:   '//TRIM( Num2LStr( Factor*UsrTimeSim    ) )//TRIM( TimePer ) )      
      WRITE(40, "(A)")  ( ' Simulated Time:        '//TRIM( Num2LStr( Factor*REAL( ZTime ) ) )//TRIM( TimePer ) )
      WRITE(40, "(A)")  ( ' Time Ratio (Sim/CPU):  '//TRIM( Num2LStr( TRatio ) ) )
      
   ENDIF

   WRITE(40, "(A)") ('...........................................................................')
   WRITE(40, "(A)") ('Parameters:')    
   WRITE(40, "(A)")  ('')   
   WRITE(40, "(A)") (' Dt of WInDS:        '//TRIM( Num2LStr( p%FVM%DT_WINDS))//' seconds')
   WRITE(40, "(A)") (' Frequency of WInDS: '//TRIM( Num2LStr( 1/p%FVM%DT_WINDS))//' Hz')
   WRITE(40, "(A)") (' Integration scheme: '//TRIM( p%FVM%INTEG ))
   
   IF (p%FVM%DS_Parms%DS_Flag) THEN
      WRITE(40, "(A)") (' Dynamic stall:      True')
   ELSE
      WRITE(40, "(A)") (' Dynamic stall:      False')
   END IF
   
   IF (p%FVM%Ground_Parms%GroundFLAG) THEN
      WRITE(40, "(A)") (' Ground effects:     True')
   ELSE
      WRITE(40, "(A)") (' Ground effects:     False')
   END IF   
   
   IF (p%FVM%Shear_Parms%ShearFLAG) THEN
      WRITE(40, "(A)") (' Shear model:        True')
   ELSE
      WRITE(40, "(A)") (' Shear model:        False')
   END IF   
        
   
   IF (p%FVM%Twr_Parms%TWRFLAG) THEN
      WRITE(40, "(A)") (' Tower shadow:       True')
   ELSE
      WRITE(40, "(A)") (' Tower shadow:       False')
   END IF      
   
   IF (p%FVM%AnimFLAG) THEN
      WRITE(40, "(A)") (' Animation:          True')
   ELSE
      WRITE(40, "(A)") (' Animation:          False')
   END IF         
   
   
   IF (p%FVM%WakeFLAG) THEN
      WRITE(40, "(A)") ('...........................................................................')
      WRITE(40, "(A)") ('Simplified wake model:')      
      WRITE(40, "(A)")  ('')   
      WRITE(40, "(A)") (' Cutoff after: '//TRIM( Num2LStr( p%FVM%WakeDist/ p%Blade%TipRadius/2  ))//' Diameter' )        
      WRITE(40, "(A)") (' Cutoff after: '//TRIM( Num2LStr( p%FVM%WakeNum))//' timesteps' ) 
      WRITE(40, "(A)") (' Freeze after: '//TRIM( Num2LStr( p%FVM%RollDist/ p%Blade%TipRadius/2 ))//' Diameter' )     
      WRITE(40, "(A)") (' Freeze after: '//TRIM( Num2LStr( O%FVM_Other%ntroll))//' timesteps' )

      
      IF (p%FVM%UindPast) THEN
         WRITE(40, "(A)") (' UindPast:     True') 
      ELSE
         WRITE(40, "(A)") (' UindPast:     False') 
      END IF      
   END IF
   
   WRITE(40, "(A)") ('...........................................................................')
   WRITE(40, "(A)") ('Acceleration:')    
   WRITE(40, "(A)") (' ')    
   WRITE(40, "(A)") (' CPU cores by user setting:   '//TRIM( Num2LStr(p%FVM%OpenMP_Parms%OpenMPCores))//' (Setting in AeroDyn_WInDS.dat)')
   WRITE(40, "(A)") (' Max CPU threads can be used: '//TRIM( Num2LStr(omp_get_num_threads()))//' (Setting in Window/Linux environment)')
   WRITE(40, "(A)") (' Max CPU cores can be used:   '//TRIM( Num2LStr(omp_get_num_procs()))//' (Setting in Window/Linux environment)')
   
   IF (p%FVM%Tree_Parms%TreeFlag) THEN
       WRITE(40, "(A)") (' Treecode:        True') 
       WRITE(40, "(A)") ('...........................................................................')
       WRITE(40, "(A)") ('Settings for Treecode Algorithm')
       WRITE(40, "(A)") (' ') 
       WRITE(40, "(A)") (' Cores:       ' // TRIM( Num2LStr( p%FVM%OpenMP_Parms%OpenMPCores )))
       WRITE(40, "(A)") (' OpenAngle:   ' // TRIM( Num2LStr( p%FVM%Tree_Parms%theta )))
       WRITE(40, "(A)") (' TaylorOrder: ' // TRIM( Num2LStr( p%FVM%Tree_Parms%order )))
       WRITE(40, "(A)") (' Maxparnode:  ' // TRIM( Num2LStr( p%FVM%Tree_Parms%maxparnode )))
       WRITE(40, "(A)") (' Dist_tol:    ' // TRIM( Num2LStr( p%FVM%Tree_Parms%dist_tol )))
       WRITE(40, "(A)") (' Smoothing:   ' // TRIM( Num2LStr( p%FVM%Tree_Parms%delta )))

       SELECT CASE (p%FVM%Tree_Parms%Freq) 
          CASE ('WHOL') 
              WRITE(40, "(A)") (' Record:      ' // TRIM( 'All timesteps' ))
          CASE ('LAST') 
              WRITE(40, "(A)") (' Record:      ' // TRIM( 'Only last timestep' ))
          CASE ('TENS')   
              WRITE(40, "(A)") (' Record:      ' // TRIM( 'Every 10 timestep' ))
          CASE ('HUND') 
              WRITE(40, "(A)") (' Record:      ' // TRIM( 'Every 100 timestep' ))
          CASE ('THOU') 
              WRITE(40, "(A)") (' Record:      ' // TRIM( 'Every 100 timestep' ))         
          END SELECT        
      
      
      
      
   ELSE
      WRITE(40, "(A)") (' Treecode:        False') 
   END IF
   
       
   CLOSE(40)  

   
   
CONTAINS
   ! .........................................................................
   FUNCTION GetClockTime(StartClockTime, EndClockTime)
   ! Copied from FAST_Subs.f90
   
   
   ! return the number of seconds between StartClockTime and EndClockTime
   
      REAL                         :: GetClockTime          ! Elapsed clock time for the simulation phase of the run.
      INTEGER   , INTENT(IN)       :: StartClockTime (8)                                 ! Start time of simulation (after initialization)
      INTEGER   , INTENT(IN)       :: EndClockTime (8)                                 ! Start time of simulation (after initialization)
   
   !bjj: This calculation will be wrong at certain times (e.g. if it's near midnight on the last day of the month), but to my knowledge, no one has complained...
      GetClockTime =       0.001*( EndClockTime(8) - StartClockTime(8) ) &  ! Is the milliseconds of the second (range 0 to 999) - local time
                     +           ( EndClockTime(7) - StartClockTime(7) ) &  ! Is the seconds of the minute (range 0 to 59) - local time
                     +      60.0*( EndClockTime(6) - StartClockTime(6) ) &  ! Is the minutes of the hour (range 0 to 59) - local time
                     +    3600.0*( EndClockTime(5) - StartClockTime(5) ) &  ! Is the hour of the day (range 0 to 23) - local time
                     + SecPerDay*( EndClockTime(3) - StartClockTime(3) )    ! Is the day of the month
   
   
   END FUNCTION     

END SUBROUTINE WInDS_WriteSum
!====================================================================================================
SUBROUTINE WRITE_KJ(Iteration, Timestep,  DG_MAX, CALLER, p)
!...................................................................................................    
    IMPLICIT                      NONE
    
   ! Passed Variables:
   INTEGER(IntKi),         INTENT(IN   )  :: Iteration         ! # of iteration when convergent
   INTEGER(IntKi),         INTENT(IN   )  :: Timestep          ! current timestep
   REAL(DbKi),             INTENT(IN   )  :: DG_MAX            ! max(D_Gamma)  
   CHARACTER(*),           INTENT(IN   )  :: CALLER            ! Who calls this subroutine 
   TYPE(AD_ParameterType), INTENT(IN   )  :: p                 ! Parameters

    ! Local Variables:   
   CHARACTER(LEN = 1024)   :: SumDIR
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE   
   CHARACTER(LEN = 1024)   :: Sumfilename
   LOGICAL                 :: FILE_EXISTS    
   CHARACTER(LEN = 1024)   :: LINE   
   CHARACTER(LEN = 1024)   :: TEMP1, TEMP2, TEMP3! , TEMP4, TEMP5, TEMP6  
   
   

   SELECT CASE  (TRIM(CALLER))    
       
   CASE ('START')
         ! Create root DIR for all simulations
       SumDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs'   
       INQUIRE (DIRECTORY = TRIM(SumDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(SumDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          !IF ( .NOT. RESULT_CREATE )  THEN        ! sliu 
          !    ErrStat = ErrID_Fatal          
          !    ErrMess = ' Error (in WInDS): Failed to create subdirectory for WInDS_Summary. Please check.'         
          !    RETURN         
          !END IF         
       ENDIF     
              
       
       SumDIR    =  TRIM(p%FVM%WINDS_dir) //  'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME)       
       INQUIRE (DIRECTORY = TRIM(SumDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(SumDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          !IF ( .NOT. RESULT_CREATE )  THEN        ! sliu 
          !    ErrStat = ErrID_Fatal          
          !    ErrMess = ' Error (in WInDS): Failed to create subdirectory for WInDS_Summary. Please check.'         
          !    RETURN         
          !END IF         
       ENDIF     
      
       Sumfilename    = '/' // TRIM(p%FVM%CURRENT_TIME) // '_KJ.dat'
   
       INQUIRE(FILE = TRIM(SumDIR)//TRIM(Sumfilename), EXIST=FILE_EXISTS)
      
       IF (.NOT. FILE_EXISTS) THEN
          OPEN (UNIT = 21, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="NEW")     
       ELSE
          OPEN (UNIT = 21, FILE = TRIM(SumDIR)//TRIM(Sumfilename), ACTION="WRITE", STATUS="REPLACE")     
          
       END IF       
              
       WRITE(21, "(A)") ('This file is used to output information for Kutta-Joukowski subroutine')
       WRITE(21, "(A)") ('.........................................................................')
       
       WRITE(TEMP1, "(F20.15)") ( p%FVM%TOL )
       LINE = 'Tolerance value for convergence of numerical methods:' // TEMP1
       WRITE(21, "(A)") (TRIM(LINE)) 
       TEMP1 = Num2LStr( p%FVM%MAXITER ) 
       LINE = 'The Maximum number of iterations for Kutta-Joukowski theorem:' // TEMP1
       WRITE(21, "(A)") (TRIM(LINE)) 
       
       WRITE(21, "(A)") ('........................................................................')
       WRITE(21, "(A)") ('Timestep.      Iteration.     Max(D_Gamma) at last iteration')
       
   CASE ('NONE') 
       
       TEMP1 = Num2LStr( Timestep ) 
       TEMP2 = Num2LStr( Iteration ) 
       WRITE(TEMP3, "(F20.15)") ( DG_MAX ) 
       LINE = TRIM(TEMP1) // ' ' // TRIM(TEMP2) // '     ' // TRIM(TEMP3) 
       
       WRITE(21, "(A)") (TRIM(LINE))         
       
              
   CASE ('END')
       
        CLOSE(21)  
       
   END SELECT    
   
END SUBROUTINE WRITE_KJ
!====================================================================================================
SUBROUTINE WRITE_INTERNAL(fila_x, fila_y, fila_z, points, GAMMA_array, RC_array, Uind_array, P, xd, O, ErrStat, ErrMess, Len_fila, Len_pts)
!...................................................................................................    
    IMPLICIT                      NONE
    
   ! Passed Variables:
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: fila_x, fila_y, fila_z 
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: points
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: GAMMA_array
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: RC_array
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: Uind_array   
   
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Initial discrete states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   INTEGER(IntKi), INTENT(OUT)                  :: ErrStat
   CHARACTER(*), INTENT(OUT)                    :: ErrMess    
   INTEGER(IntKi), INTENT(IN)                   :: Len_fila
   INTEGER(IntKi), INTENT(IN)                   :: Len_pts
   

    ! Local Variables:   
   CHARACTER(LEN = 1024)   :: TreeDIR
   LOGICAL                 :: DIR_EXISTS
   LOGICAL(4)              :: RESULT_CREATE   
   CHARACTER(LEN = 1024)   :: Treefilename
   LOGICAL                 :: FILE_EXISTS    
   CHARACTER(LEN = 1024)   :: LINE   
   CHARACTER(LEN = 1024)   :: TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6 
   CHARACTER(LEN=8)        :: DATE_TEMP
   CHARACTER(LEN=10)       :: TIME_TEMP
   
   
   INTEGER(IntKi)               :: StrtTime (8)                                    ! Start time of simulation (including initialization)
   INTEGER(IntKi)               :: SimStrtTime (8)                                 ! Start time of simulation (after initialization)
   REAL(DbKi)                   :: UsrTime1                                        ! User CPU time for simulation initialization.
   REAL(DbKi)                   :: UsrTime2                                        ! User CPU time for simulation (without intialization)
   REAL(DbKi)                   :: ZTime                                           ! The final simulation time (not necessarially TMax)
   !REAL,OPTIONAL, INTENT(OUT)   :: UsrTime_out                                     ! User CPU time for entire run - optional value returned to calling routine

      ! Local variables
   INTEGER(IntKi)               :: K, I 
   REAL(DbKi)                   :: ClckTime                                        ! Elapsed clock time for the entire run.
   REAL(DbKi)                   :: ClckTimeSim                                     ! Elapsed clock time for the simulation phase of the run.
   REAL(DbKi)                   :: Factor                                          ! Ratio of seconds to a specified time period.
   REAL(DbKi)                   :: TRatio                                          ! Ratio of simulation time to elapsed clock time.
   REAL(DbKi), PARAMETER        :: SecPerDay = 24*60*60.0_ReKi                     ! Number of seconds per day

   REAL(DbKi)                   :: UsrTime                                         ! User CPU time for entire run.
   REAL(DbKi)                   :: UsrTimeSim                                      ! User CPU time for simulation (not including initialization).
   INTEGER(IntKi)               :: EndTimes (8)                                    ! An array holding the ending clock time of the simulation.

   CHARACTER( 8)                :: TimePer
   CHARACTER(MaxWrScrLen)       :: BlankLine  
             
 

       !....................................................................... 
         ! Create root DIR for all simulations
       TreeDIR    =  TRIM(p%FVM%WINDS_dir) // 'WInDS_outputs'   
       INQUIRE (DIRECTORY = TRIM(TreeDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(TreeDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          IF ( .NOT. RESULT_CREATE )  THEN        
              ErrStat = ErrID_Fatal          
              ErrMess = ' Error (in WInDS): Failed to create subdirectory for WInDS_INTERNAL. Please check.'         
              RETURN         
          END IF         
       ENDIF     
 
          
       TreeDIR    =  TRIM(p%FVM%WINDS_dir) //  'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME)       
       INQUIRE (DIRECTORY = TRIM(TreeDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(TreeDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          IF ( .NOT. RESULT_CREATE )  THEN       
              ErrStat = ErrID_Fatal          
              ErrMess = ' Error (in WInDS): Failed to create subdirectory for # WInDS_INTERNAL #. Please check.'         
              RETURN         
          END IF         
       ENDIF     

    
       
       TreeDIR    =  TRIM(p%FVM%WINDS_dir) //  'WInDS_outputs/' //TRIM(p%FVM%CURRENT_TIME) //'/treecode'      
       INQUIRE (DIRECTORY = TRIM(TreeDIR), EXIST = DIR_EXISTS)      
       IF (.NOT. DIR_EXISTS) THEN
          RESULT_CREATE = MAKEDIRQQ (TRIM(TreeDIR))         !  refer: software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rfmkdir.htm
        
          IF ( .NOT. RESULT_CREATE )  THEN       
              ErrStat = ErrID_Fatal          
              ErrMess = ' Error (in WInDS): Failed to create subdirectory for ## WInDS_INTERNAL ##. Please check.'         
              RETURN         
          END IF         
       ENDIF       
       

       TEMP1 = Num2LStr( O%WINDS_Timestep )
       TEMP2 = Num2LStr( Len_fila )
       TEMP3 = Num2LStr( Len_pts )
       
       CALL DATE_AND_TIME(DATE = DATE_TEMP, TIME = TIME_TEMP)
       TEMP4   =  DATE_TEMP//'_'//TIME_TEMP(1:6)   ! hhmmss - CCYYMMDD      Refer: http://docs.oracle.com/cd/E19957-01/805-4942/6j4m3r8t2/index.html
                                    ! The current time, in the form hhmmss.sss, where hh is the hour, mm minutes, and ss.sss seconds and milliseconds.
                                    ! Date, in form CCYYMMDD, where CCYY is the four-digit year, MM the two-digit month, and DD the two-digit day of the month.
   
       
       Treefilename    = '/WInDS_internal_'// TRIM(TEMP1) // '_'// TRIM(TEMP2) // '_' //TRIM(TEMP3) // '_' //TRIM(TEMP4)  // '.dat'
   
       INQUIRE(FILE = TRIM(TreeDIR)//TRIM(Treefilename), EXIST=FILE_EXISTS)
             
       IF (.NOT. FILE_EXISTS) THEN
          OPEN (UNIT = 50, FILE = TRIM(TreeDIR)//TRIM(Treefilename), ACTION="WRITE", STATUS="NEW")     
       ELSE
          OPEN (UNIT = 50, FILE = TRIM(TreeDIR)//TRIM(Treefilename), ACTION="WRITE", STATUS="REPLACE")           
       END IF      
       
       
       WRITE(50, "(A)") ('Timestep number:              ' // TRIM(Num2LStr(O%WINDS_Timestep)))
       WRITE(50, "(A)") ('Real time (min):              ' // TRIM(Num2LStr(O%WINDS_Timestep*p%FVM%DT_WINDS/60)))
       WRITE(50, "(A)") ('Number of filaments (source): ' // TRIM(TEMP2))      
       WRITE(50, "(A)") ('Number of points (receiver):  ' // TRIM(TEMP3))            
       WRITE(50, "(A)") ('.......................................................')
       WRITE(50, "(A)") ('fila_x:')

       DO k=1, Len_fila
          WRITE(50, "(A)") (TRIM(Num2LStr(fila_x(K,1)))// '   '// TRIM(Num2LStr(fila_x(K,2)))// '   '// TRIM(Num2LStr(fila_x(K,3))) ) 
       END DO
       
       WRITE(50, "(A)") ('.......................................................')
       WRITE(50, "(A)") ('fila_y:')

       DO k=1, Len_fila
          WRITE(50, "(A)") (TRIM(Num2LStr(fila_y(K,1)))// '   '// TRIM(Num2LStr(fila_y(K,2)))// '   '// TRIM(Num2LStr(fila_y(K,3))) ) 
       END DO       
       
       WRITE(50, "(A)") ('.......................................................')
       WRITE(50, "(A)") ('fila_z:')

       DO k=1, Len_fila
          WRITE(50, "(A)") (TRIM(Num2LStr(fila_z(K,1)))// '   '// TRIM(Num2LStr(fila_z(K,2)))// '   '// TRIM(Num2LStr(fila_z(K,3))) ) 
       END DO       
              
       WRITE(50, "(A)") ('.......................................................')
       WRITE(50, "(A)") ('GAMMA:')

       DO k=1, Len_fila
          WRITE(50, "(A)") (TRIM(Num2LStr(GAMMA_array(K,1))) ) 
       END DO              
       
       WRITE(50, "(A)") ('.......................................................')
       WRITE(50, "(A)") ('RC:')

       DO k=1, Len_fila
          WRITE(50, "(A)") (TRIM(Num2LStr(RC_array(K,1))) ) 
       END DO              
          
       WRITE(50, "(A)") ('.......................................................')
       WRITE(50, "(A)") ('Points:')

       DO k=1, Len_pts
          WRITE(50, "(A)") (TRIM(Num2LStr(points(K,1)))// '   '// TRIM(Num2LStr(points(K,2)))// '   '// TRIM(Num2LStr(points(K,3))) ) 
       END DO       
   
       
      CLOSE(50)  
   
    
    
    
END SUBROUTINE WRITE_INTERNAL
!====================================================================================================
SUBROUTINE Print_element(P, xd, O, ErrStat, ErrMess)
!
!...................................................................................................    
    IMPLICIT                      NONE
   
   ! Passed Variables:
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Initial discrete states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   INTEGER(IntKi), INTENT(OUT)                  :: ErrStat
   CHARACTER(*), INTENT(OUT)                    :: ErrMess    
    
   

      ! Local variables
   REAL(DbKi), DIMENSION(p%Element%NElm, p%NumBl)               :: PRINT_NAME1
   REAL(DbKi), DIMENSION(p%Element%NElm + 1, p%NumBl)           :: PRINT_NAME2
   REAL(DbKi), DIMENSION(p%Element%NElm, 3)                     :: PRINT_NAME3   
   
   INTEGER(IntKi)      :: IDim        ! For dimensions
   INTEGER(IntKi)      :: ITimestep   ! For timesteps
   
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB 
   
   CHARACTER(LEN=1024)                  :: temp_number   
   CHARACTER(LEN=2), DIMENSION(20)      :: Rank_num   
   CHARACTER(LEN=1)                     :: temp_1 
   CHARACTER(LEN=2)                     :: temp_2 
   CHARACTER(LEN=3)                     :: temp_3 
 
   INTEGER                    :: IBlade
   INTEGER                    :: IElement
   INTEGER                    :: Node              ! Node index.
  
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl  

   
   WRITE (temp_number , "(I6.6)") ( O%WINDS_Timestep )   
                      
   ! CL
   PRINT_NAME1 = 0.0
   DO IElement = 1, NS
       DO IBlade =  1,NB
           PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_CL(1, O%WINDS_Timestep, 1, IElement, IBlade) 
       END DO  
   END DO  
   CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_cl_'//TRIM(temp_number))    ! Write cl to text
                       
   ! CD
   PRINT_NAME1 = 0.0
   DO IElement = 1, NS
       DO IBlade =  1,NB
           PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_CD(1, O%WINDS_Timestep, 1, IElement, IBlade) 
       END DO  
   END DO  
   CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_cd_'//TRIM(temp_number))    ! Write cd to text    
                       
                       
   ! AOA
   PRINT_NAME1 = 0.0
   DO IElement = 1, NS
       DO IBlade =  1,NB
          PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_AOA(1, O%WINDS_Timestep, 1, IElement, IBlade) 
       END DO    
   END DO    
   CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_aoa_'// TRIM(temp_number))  
                           

   !V_tot                           
   DO IBlade =  1,NB
       WRITE (temp_1, "(I1)") (IBlade)    ! String of the integer With heading zeros 
       PRINT_NAME3 = 0.0
       DO IElement = 1, NS
          DO IDim =  1,3
              PRINT_NAME3(IElement, IDIM) = O%FVM_Other%KJ%VEL_TOT(IDIM, O%WINDS_Timestep, 1, IElement, IBlade) 
          END DO    
       END DO    
       CALL SAVE_TO_TXT_2D(p, PRINT_NAME3 , 'WINDS_vtot_'//temp_1//'_blade_'// TRIM(temp_number)) 
    END DO               
             

END SUBROUTINE print_element
!====================================================================================================







END MODULE WINDS_IO
!**********************************************************************************************************************************   