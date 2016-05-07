!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
!
!    Lidar module, a submodule of InflowWind
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! File last committed: $Date: 2015-02-09 11:06:27 -0700 (Mon, 09 Feb 2015) $
! (File) Revision #: $Rev: 140 $
! URL: $HeadURL: https://windsvn.nrel.gov/InflowWind/branches/modularization/Source/Lidar.f90 $
!**********************************************************************************************************************************
MODULE Lidar

   USE InflowWind_Subs
   USE Lidar_Types

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: Lidar_Ver = ProgDesc( 'Lidar', 'v1.00.00a-bjj', '30-Jan-2015' )
   CHARACTER(*),   PARAMETER            :: Lidar_Nickname = 'Lidar'
   
   
   REAL(ReKi),     PARAMETER            :: BeamRad      =  0.028                            
   REAL(ReKi),     PARAMETER            :: LsrWavLen    =  0.000001565                      ! Laser wavelength

   
! ==================================================================================================="

      
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Lidar_Init                           ! Initialization routine
   PUBLIC :: Lidar_End                            ! Ending routine (includes clean up)
   PUBLIC :: Lidar_CalcOutput                     ! Routine for computing outputs
 
   
!bjj: to do:
! + add mesh to map nacelle rotor apex position in ElastoDyn to lidar location (possibly an array of lidars) in InflowWind
! + add input file (part of InflowWind input file)
!    - number of lidars, type, location, number of pulse range gates, etc
!    - initial measurement position(s)
!    - scan pattern & associated values [remove this functionality from Matlab]
! + add subroutine for scanning patterns
! future work:
! + do we want to know if the blade is in front of the lidar so we can return garbage to simulate that scenario, too?
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Lidar_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
! note that we're calling this with the InflowWind data types, so that data is INOUT instead of OUT
!..................................................................................................................................

   TYPE(IfW_InitInputType),         INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(IfW_InputType),             INTENT(INOUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(IfW_ParameterType),         INTENT(INOUT)  :: p           ! Parameters
   TYPE(IfW_ContinuousStateType),   INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(IfW_DiscreteStateType),     INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(IfW_ConstraintStateType),   INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(IfW_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Initial other/optimization states
   TYPE(IfW_OutputType),            INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                  !   only the output mesh is initialized)
   REAL(DbKi),                      INTENT(IN )    :: Interval    ! Coupling interval in seconds: the rate that
                                                                  !   (1) Lidar_UpdateStates() is called in loose coupling &
                                                                  !   (2) Lidar_UpdateDiscState() is called in tight coupling.
                                                                  !   Input is the suggested time from the glue code;
                                                                  !   Output is the actual coupling interval that will be used
                                                                  !   by the glue code.
   TYPE(IfW_InitOutputType),        INTENT(INOUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

   INTEGER(IntKi)                                  :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(LEN(ErrMsg))                          :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'Lidar_Init'
   
   REAL(ReKi)                                      :: TempWindSpeed(3)

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

      
   !   ! Initialize the NWTC Subroutine Library
   !
   !CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   !CALL DispNVD( Lidar_Ver )

      !............................................................................................      
      ! Read the input file and validate the data
      !............................................................................................      
   !p%RootName = TRIM(InitInp%RootName)//'.'//Lidar_Nickname ! all of the output file names from this module will end with '.ModName'
      
      
      !............................................................................................
      ! Define parameters here:
      !............................................................................................
      
   p%lidar%RotorApexOffsetPos = InitInp%lidar%RotorApexOffsetPos
      
   p%lidar%SensorType = InitInp%lidar%SensorType      
   IF (p%lidar%SensorType == SensorType_None) THEN
      p%lidar%NumPulseGate = 0
   ELSEIF (p%lidar%SensorType == SensorType_SinglePoint) THEN
      p%lidar%NumPulseGate = 1
   ELSE
      
         ! variables for both pulsed and continuous-wave lidars
      
      CALL InflowWind_GetMean(0.0_DbKi, InitInp%lidar%TMax, Interval, InitInp%lidar%HubPosition, TempWindSpeed, &
                              p, x, xd, z, OtherState, ErrStat2, ErrMsg2 )
     
      p%lidar%SpatialRes     =  0.5_ReKi*TempWindSpeed(1)*Interval      
      p%lidar%RayRangeSq     =  (Pi*(BeamRad**2)/LsrWavLen)**2
   
      p%lidar%LidRadialVel   = InitInp%lidar%LidRadialVel  !.FALSE.
   

      IF (p%lidar%SensorType == SensorType_ContinuousLidar) THEN
      
         p%lidar%WtFnTrunc    = 0.02_ReKi   
         p%lidar%NumPulseGate = 1
   
      ELSEIF (p%lidar%SensorType == SensorType_PulsedLidar) THEN
      
         p%lidar%WtFnTrunc    = 0.01_ReKi      
         p%lidar%NumPulseGate = InitInp%lidar%NumPulseGate
            
            ! values for the WindCube
         p%lidar%DeltaP        = 30.0_ReKi
         p%lidar%DeltaR        = 30.0_ReKi
         p%lidar%PulseRangeOne = 35.0_ReKi
            
         p%lidar%r_p           = p%lidar%DeltaR/(2.0_ReKi*SQRT(LOG(2.0_ReKi)))
         
      ELSE
         
         CALL SetErrStat(ErrID_Fatal, "Invalid sensor type.", ErrStat, ErrMsg, RoutineName)
         RETURN
         
      END IF
      
   END IF  
   
  
      !............................................................................................
      ! Define initial system states here:
      !............................................................................................

   !x%lidar%DummyContState           = 0.0_ReKi
   !xdlidar%%DummyDiscState          = 0.0_ReKi
   !z%lidar%DummyConstrState         = 0.0_ReKi
   !OtherState%lidar%DummyOtherState = 0.0_ReKi
   !

      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................

   u%lidar%LidPosition = InitInp%lidar%HubPosition
   u%lidar%MsrPosition = InitInp%lidar%HubPosition + (/ 50.0, 0.0, 0.0 /) !bjj: todo FIXME  with initial guess of lidar focus.
   u%lidar%PulseLidEl  = 0.0_ReKi
   u%lidar%PulseLidAz  = 0.0_ReKi
   
   
      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................
   !CALL AllocAry( y%WriteOutput, p%NumOuts, 'WriteOutput', ErrStat2, ErrMsg2 )
   !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   !   IF (ErrStat >= AbortErrLev) RETURN
   !y%WriteOutput = 0
   
   
   CALL AllocAry( y%lidar%LidSpeed, p%lidar%NumPulseGate, 'y%lidar%LidSpeed', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );

   !CALL AllocAry( y%LidErr, p%NumPulseGate, 'y%LidErr', ErrStat2, ErrMsg2 )
   !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
   CALL AllocAry( y%lidar%WtTrunc, p%lidar%NumPulseGate, 'y%lidar%WtTrunc', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
   
      
   IF (ErrStat >= AbortErrLev) RETURN
   y%lidar%LidSpeed = 0.0
   y%lidar%WtTrunc  = 0.0
               
      !............................................................................................
      ! Define initialization-routine output here:
      !............................................................................................
   !CALL AllocAry( InitOut%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', ErrStat2, ErrMsg2 )
   !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   !CALL AllocAry( InitOut%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', ErrStat2, ErrMsg2 )
   !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   !   
   !IF (ErrStat >= AbortErrLev) RETURN
   
   !InitOut%WriteOutputHdr = p%OutParam(1:p%NumOuts)%Name
   !InitOut%WriteOutputUnt = p%OutParam(1:p%NumOuts)%Units     
   !InitOut%Ver = Lidar_Ver            
   
                  
   RETURN
   
END SUBROUTINE Lidar_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Lidar_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(IfW_InputType),             INTENT(INOUT)  :: u           ! System inputs
      TYPE(IfW_ParameterType),         INTENT(INOUT)  :: p           ! Parameters
      TYPE(IfW_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states
      TYPE(IfW_DiscreteStateType),     INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(IfW_ConstraintStateType),   INTENT(INOUT)  :: z           ! Constraint states
      TYPE(IfW_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(IfW_OutputType),            INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""
         
         ! Close files here:
         

         ! Destroy the input data:

      CALL Lidar_DestroyInput( u%lidar, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL Lidar_DestroyParam( p%lidar, ErrStat, ErrMsg )


         ! Destroy the state data:

      !CALL Lidar_DestroyContState(   x,           ErrStat, ErrMsg )
      !CALL Lidar_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      !CALL Lidar_DestroyConstrState( z,           ErrStat, ErrMsg )
      !CALL Lidar_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )


         ! Destroy the output data:

      CALL Lidar_DestroyOutput( y%lidar, ErrStat, ErrMsg )




END SUBROUTINE Lidar_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Lidar_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! Routine for computing outputs, used in both loose and tight coupling.
! Note, this breaks the framework because we're passing the IfW types instead of the Lidar types... this is necessary to get
!    appropriate wind speeds for the lidar measurements. 
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   )  :: t                    ! Current simulation time in seconds
   TYPE(IfW_InputType),             INTENT(IN   )  :: u                    ! Inputs at t
   TYPE(IfW_ParameterType),         INTENT(IN   )  :: p                    ! Parameters
   TYPE(IfW_ContinuousStateType),   INTENT(IN   )  :: x                    ! Continuous states at t
   TYPE(IfW_DiscreteStateType),     INTENT(IN   )  :: xd                   ! Discrete states at t
   TYPE(IfW_ConstraintStateType),   INTENT(IN   )  :: z                    ! Constraint states at t
   TYPE(IfW_OtherStateType),        INTENT(INOUT)  :: OtherState           ! Other/optimization states
   TYPE(IfW_OutputType),            INTENT(INOUT)  :: y                    ! Outputs computed at t (Input only so that mesh con-
                                                                           !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat              ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
            
   REAL(ReKi)                                      :: FocDist              ! Focus Distance of the laser
   REAL(ReKi)                                      :: FocDistMod           ! Focus Distance of the laser?
   REAL(ReKi)                                      :: LidPhi               ! angle used with LidTheta to describe the direction the lidar is pointed
   REAL(ReKi)                                      :: LidTheta             ! angle used with LidPhi to describe the direction the lidar is pointed
   REAL(ReKi)                                      :: LidRange             ! lidar range
                                                                           
   REAL(ReKi)                                      :: WtFuncSum            ! sum of weight function, used to normalize summation to 1
   REAL(ReKi)                                      :: LidWt                ! The weighting function value
   REAL(ReKi)                                      :: LidWtMax             ! maximum weighting function value?
   REAL(ReKi)                                      :: LidWtRatio           ! LidWt/LidWtMax
   REAL(ReKi)                                      :: LidDirUnVec(3)       ! lidar look direction unit vector
            
   REAL(ReKi)                                      :: Distance(3)          ! distance vector between input measurement and lidar positions
   
   TYPE(IfW_InputType)                             :: Input                ! position where wind speed should be returned
   TYPE(IfW_OutputType)                            :: Output               ! velocity at Input%Position
   
   REAL(ReKi)                                      :: OutputVelocity(3)
   
      
   INTEGER(IntKi)                                  :: IRangeGt
   INTEGER(IntKi)                                  :: ErrStat2
   CHARACTER(LEN(ErrMsg))                          :: ErrMsg2
         
   CHARACTER(*), PARAMETER                         :: RoutineName = 'Lidar_CalcOutput'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (p%lidar%SensorType == SensorType_None) RETURN
   
      
      ! allocate arrays to compute outputs
   CALL AllocAry(Input%Position, 3,1, 'Input%Position',ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL AllocAry(Output%Velocity, 3,1, 'Output%Position',ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN      
   END IF
   
   !...............................................................................................................................   
   ! Compute the outputs
   !...............................................................................................................................   


   IF (p%lidar%SensorType == SensorType_SinglePoint) THEN
      
      !get lidar speed at the focal point to see if it is out of bounds   
      Input%Position(:,1) = u%lidar%MsrPosition      
      CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                
      
      y%lidar%LidSpeed = Output%Velocity(:,1)
      y%lidar%WtTrunc  = 1.0_ReKi
         
   ELSEIF (p%lidar%SensorType == SensorType_ContinuousLidar) THEN
      !calculate the focal distance of the lidar as well as the modified focal distance so that the peak of the weighting func
      !is at the intended focal distance
   
      Distance   = u%lidar%MsrPosition - u%lidar%LidPosition      
      FocDist    = SQRT( DOT_PRODUCT( Distance, Distance ) ) !TwoNorm
      
      IF(EqualRealNos(FocDist,0.0_ReKi)) THEN ! Avoid division-by-zero
         y%lidar%LidSpeed = -99.0
         y%lidar%WtTrunc  = 0.0         
         CALL SetErrStat(ErrID_Fatal,"Measurement position cannot be the same as the lidar position.", ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
      FocDistMod = (p%lidar%RayRangeSq - SQRT(p%lidar%RayRangeSq**2 - 4*p%lidar%RayRangeSq*(FocDist**2)))/(2*FocDist);
   
   
      !Find angles that the lidar is pointed at
      LidPhi = ATAN(Distance(3)/SQRT(Distance(1)**2 + Distance(2)**2))
      !LidTheta = ATAN(Distance(2)/ABS(Distance(1)))
      LidTheta = ATAN2(Distance(2),-Distance(1))
   
   
      !calculate the unit vector of the lidar look direction
      LidDirUnVec = Distance / FocDist
      LidWt = 1.0/(FocDist**2 + ((1 - FocDist/FocDistMod)**2)*p%lidar%RayRangeSq)
      LidWtMax = LidWt
      LidWtRatio = 1.0_ReKi !LidWt/LidWtMax
   
      !get lidar speed at the focal point to see if it is out of bounds   
      Input%Position(:,1) = u%lidar%MsrPosition      
      CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                
    
         !if out of bounds
      IF (ErrStat >= AbortErrLev) THEN
         !y%lidar%LidErr = 1
         y%lidar%LidSpeed = -99.0
         CALL Cleanup()
         RETURN !escape function
      ENDIF
    
      y%lidar%LidSpeed = LidWt*DOT_PRODUCT(-1*LidDirUnVec,Output%Velocity(:,1))
    
      WtFuncSum = LidWt
      y%lidar%WtTrunc = p%lidar%WtFnTrunc
    
      !initialize lidar range
      LidRange = 0.
   
   
   
      !calculate the weighted lidar returns
      DO 
      
         !escape loop if desired truncation point has been reached 
         IF (LidWtRatio < p%lidar%WtFnTrunc) THEN
            EXIT
         ENDIF
       
       
         !calculate range of current beam point
         LidRange = LidRange + p%lidar%SpatialRes
          
         LidWt = 1.0/((FocDist + LidRange)**2 + ((1 - (FocDist + LidRange)/FocDistMod)**2)*p%lidar%RayRangeSq)
         LidWtRatio = LidWt/LidWtMax

       
         !trunc point is behind lidar
         IF (LidRange > FocDist) THEN
            IF (NWTC_VerboseLevel == NWTC_Verbose) &
               CALL SetErrStat( ErrID_Info, "Lidar truncation point is behind the lidar. Truncation ratio is "//trim(num2lstr(LidWtRatio))//'.', ErrStat, ErrMsg, RoutineName)  ! set informational message about point being behind lidar
            !y%LidErr = 3
            y%lidar%WtTrunc = LidWtRatio
            EXIT
         ENDIF    
          
          
         !calculate points to scan for current beam point
         Input%Position(3,1) = u%lidar%LidPosition(3) + SIN(LidPhi)*(LidRange + FocDist)
         Input%Position(1,1) = u%lidar%LidPosition(1) - COS(LidTheta)*COS(LidPhi)*(LidRange + FocDist)
         Input%Position(2,1) = u%lidar%LidPosition(2) + SIN(LidTheta)*COS(LidPhi)*(LidRange + FocDist)                     
          
         CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, ErrStat2, ErrMsg2 )    
            IF (ErrStat2 >= AbortErrLev ) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
            !y%LidErr = 2
               y%lidar%WtTrunc = LidWtRatio
               EXIT
            ENDIF
                
         OutputVelocity = Output%Velocity(:,1)
          
       
         !calculate points to scan for current beam point
         Input%Position(3,1) = u%lidar%LidPosition(3) + SIN(LidPhi)*(FocDist - LidRange)
         Input%Position(1,1) = u%lidar%LidPosition(1) - COS(LidTheta)*COS(LidPhi)*(FocDist - LidRange)
         Input%Position(2,1) = u%lidar%LidPosition(2) + SIN(LidTheta)*COS(LidPhi)*(FocDist - LidRange)
       
         CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, ErrStat2, ErrMsg2 )      
            IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
            !y%lidar%LidErr = 2
               y%lidar%WtTrunc = LidWtRatio
               EXIT
            ENDIF
       
      
         y%lidar%LidSpeed = y%lidar%LidSpeed + LidWt*DOT_PRODUCT(-1*LidDirUnVec, OutputVelocity + Output%Velocity(:,1))      
         WtFuncSum = WtFuncSum + 2*LidWt
           
       
      END DO
   
      !Normalize the weighting function summation to 1
   
      IF ( p%lidar%LidRadialVel )  THEN
            !This detects the radial component
            y%lidar%LidSpeed = y%lidar%LidSpeed/WtFuncSum
      ELSE
            !This returns the 'x' component estimate
            y%lidar%LidSpeed = -1*y%lidar%LidSpeed/(WtFuncSum*LidDirUnVec(1))
      ENDIF   
   
   ELSE !p%SensorType == SensorType_PulsedLidar
      
      
      
      LidDirUnVec(1) = -1*COS(u%lidar%PulseLidEl)
      LidDirUnVec(2) = SIN(u%lidar%PulseLidEl)*SIN(u%lidar%PulseLidAz)
      LidDirUnVec(3) = SIN(u%lidar%PulseLidEl)*COS(u%lidar%PulseLidAz)
   
      
      DO IRangeGt = 1,p%lidar%NumPulseGate
   
         !y%lidar%LidErr(IRangeGt) = 0
   
         !get lidar speed at the focal point to see if it is out of bounds
         Input%Position(:,1) = u%lidar%LidPosition + LidDirUnVec*(p%lidar%PulseRangeOne + (IRangeGt-1)*p%lidar%DeltaP)                 
         CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                
         LidWt = NWTC_ERF((p%lidar%DeltaP/2)/p%lidar%r_p)/p%lidar%DeltaP        
         LidWtMax = LidWt
         LidWtRatio = 1.0_ReKi !LidWt/LidWtMax
       
        
         !if out of bounds
         IF (AbortErrLev >= AbortErrLev) THEN
            !y%LidErr(IRangeGt) = 1
            y%lidar%LidSpeed(IRangeGt) = -99
            CALL Cleanup()
            RETURN !escape function
         ENDIF
        
         y%lidar%LidSpeed(IRangeGt) = LidWt*DOT_PRODUCT(-1*LidDirUnVec,Output%Velocity(:,1))
        
         WtFuncSum = LidWt
         y%lidar%WtTrunc(IRangeGt) = p%lidar%WtFnTrunc
        
         !initialize lidar range
         LidRange = 0.
   
         DO
            !escape loop if desired truncation point has been reached 
            IF (LidWtRatio < p%lidar%WtFnTrunc) THEN
                  EXIT
            ENDIF
           
            !calculate range of current beam point
            LidRange = LidRange + p%lidar%SpatialRes
              
            LidWt = (NWTC_ERF((LidRange + p%lidar%DeltaP/2.)/p%lidar%r_p) - NWTC_ERF((LidRange - p%lidar%DeltaP/2.)/p%lidar%r_p))/(2.*p%lidar%DeltaP)
            LidWtRatio = LidWt/LidWtMax
           
            
               !trunc point is behind lidar
            IF (LidRange > (p%lidar%PulseRangeOne + (IRangeGt-1)*p%lidar%DeltaP)) THEN
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Info, "Lidar truncation point at gate "//trim(num2lstr(IRangeGt))//" is behind the lidar. Truncation ratio is "&
                                 //trim(num2lstr(LidWtRatio))//'.', ErrStat, ErrMsg, RoutineName)  ! set informational message about point being behind lidar
               !y%LidErr(IRangeGt) = 3
               y%lidar%WtTrunc(IRangeGt) = LidWtRatio
               EXIT
            ENDIF
                                     
           
            !calculate points to scan for current beam point
            Input%Position(:,1) = u%lidar%LidPosition + LidDirUnVec*(p%lidar%PulseRangeOne + (IRangeGt-1)*p%lidar%DeltaP + LidRange)              
            CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, ErrStat2, ErrMsg2 )   
               IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed at gate "//trim(num2lstr(IRangeGt))//" truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
                  !y%LidErr(IRangeGt) = 2
                  y%lidar%WtTrunc(IRangeGt) = LidWtRatio
                  EXIT
               ENDIF
            OutputVelocity = Output%Velocity(:,1)
                                                                               
            !calculate points to scan for current beam point
            Input%Position(:,1) = u%lidar%LidPosition + LidDirUnVec*(p%lidar%PulseRangeOne + (IRangeGt-1)*p%lidar%DeltaP - LidRange)    
            CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, ErrStat2, ErrMsg2 )      
               IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed at gate "//trim(num2lstr(IRangeGt))//" truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
                  !y%lidar%LidErr(IRangeGt) = 2
                  y%lidar%WtTrunc(IRangeGt) = LidWtRatio
                  EXIT
               ENDIF
           
           
            y%lidar%LidSpeed(IRangeGt) = y%lidar%LidSpeed(IRangeGt) + LidWt*DOT_PRODUCT(-1*LidDirUnVec,Output%Velocity(:,1) + OutputVelocity)           
            WtFuncSum = WtFuncSum + 2*LidWt
           
         END DO

       
         IF ( p%lidar%LidRadialVel )  THEN
            !This detects the radial component
            y%lidar%LidSpeed(IRangeGt) = y%lidar%LidSpeed(IRangeGt)/WtFuncSum
         ELSE
            !This returns the 'x' component estimate
            y%lidar%LidSpeed(IRangeGt) = -1*y%lidar%LidSpeed(IRangeGt)/(LidDirUnVec(1)*WtFuncSum)
         ENDIF
   
      END DO      
      
      
   END IF !type of lidar measurements
         
   
   CALL Cleanup()
         
   RETURN
CONTAINS
   SUBROUTINE Cleanup()
      
      IF (ALLOCATED(Input%Position)) DEALLOCATE(Input%Position)
      IF (ALLOCATED(Output%Velocity)) DEALLOCATE(Output%Velocity)
   
   END SUBROUTINE Cleanup
   
END SUBROUTINE Lidar_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Lidar
!**********************************************************************************************************************************
