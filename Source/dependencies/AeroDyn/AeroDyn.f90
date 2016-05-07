!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
! File last committed: $Date: 2015-03-13 12:04:56 -0600 (Fri, 13 Mar 2015) $
! (File) Revision #: $Rev: 177 $
! URL: $HeadURL: https://windsvn.nrel.gov/AeroDyn/trunk/Source/AeroDyn.f90 $
!**********************************************************************************************************************************
!
!!!!!!! Things to care for WInDS:  Transfer max time step   ! sliu    !!!!!!!!!

    
    
MODULE AeroDyn

   USE AeroDyn_Types
   USE AeroSubs
   USE NWTC_Library

   
   !-----Umass WINDS-----
   USE WINDS          ! Wake Induced Dynamics Simulator
   USE WINDS_IO       ! Deal with input and ouput files
   USE WINDS_DS       ! LB Dynamic stall
   USE WINDS_Library  ! For debug use, write internal variables to txt file...(Feel free to ask to auther for this)
   !-----Umass WINDS-----   
   

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: AD_Ver = ProgDesc( 'AeroDyn', 'v14.03.01a-bjj', '13-Mar-2015' )

      ! ..... Public Subroutines ............

   PUBLIC :: AD_Init                           ! Initialization routine
   PUBLIC :: AD_End                            ! Ending routine (includes clean up)

   PUBLIC :: AD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   PUBLIC :: AD_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: AD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: AD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: AD_UpdateDiscState                ! Tight coupling routine for updating discrete states


! Note that the following routines will be updated with new definitions of arrays returned (no longer one-byte arrays)
!   PUBLIC :: AD_Pack                           ! Routine to pack (save) data into one array of bytes
!   PUBLIC :: AD_Unpack                         ! Routine to unpack an array of bytes into data structures usable by the module

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_Init( InitInp, u, p, x, xd, z, O, y, Interval, InitOut, ErrStat, ErrMess )
!..................................................................................................................................
   USE               AeroGenSubs,   ONLY: ElemOpen
   USE DWM
   IMPLICIT NONE

   TYPE(AD_InitInputType),       INTENT(INOUT)  :: InitInp     ! Input data for initialization routine
   TYPE(AD_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(AD_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(AD_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(AD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(AD_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD_OtherStateType),      INTENT(  OUT)  :: O !therState  Initial other/optimization states
   TYPE(AD_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                               !   only the output mesh is initialized)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                               !   (1) AD_UpdateStates() is called in loose coupling &
                                                               !   (2) AD_UpdateDiscState() is called in tight coupling.
                                                               !   Input is the suggested time from the glue code;
                                                               !   Output is the actual coupling interval that will be used
                                                               !   by the glue code.
   TYPE(AD_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


      ! Internal variables
   REAL(ReKi)                       :: CosPrecone
   REAL(ReKi)                       :: DTip, ElemRad, Dhub, Rhub     ! variables for calculating hub- and tip-loss constants
   REAL(ReKi)                       :: HubRadius
   REAL(ReKi)                       :: MeanWind
   REAL(ReKi)                       :: TipRadius
   REAL(ReKi)                       :: TmpVar
   REAL(ReKi)                       :: TmpPos(3)
   REAL(ReKi)                                :: TwrNodeHt                     ! The height of the current tower node.

   INTEGER                          :: IB, IE 
   INTEGER                          :: IELM

   CHARACTER(1024)                  :: Title

   INTEGER                                   :: Elem                          ! Index for mesh element.
   INTEGER                                   :: InterpIndx                 ! Index telling the interpolation routine where to start in the array.
   INTEGER                                   :: Node                          ! Index used to pull points out of the array of values at given node location.
   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.

   CHARACTER(LEN(ErrMess))                   :: ErrMessLcl          ! Error message returned by called routines.
   CHARACTER(*), PARAMETER                   :: RoutineName = 'AD_Init'

         ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMess  = ""
   InterpIndx = 1
   !-------------------------------------------------------------------------------------------------
   ! Check that the module hasn't already been initialized.
   !-------------------------------------------------------------------------------------------------

   IF ( p%Initialized ) THEN
      CALL SetErrStat( ErrID_Warn,'AeroDyn has already been initialized.',ErrStat,ErrMess,RoutineName)
      RETURN
   ELSE
      p%Initialized = .TRUE.
      CALL NWTC_Init( )
   END IF

   
         ! Display the module information

   CALL DispNVD( AD_Ver )
   
   InitOut%Ver = AD_Ver
   O%FirstWarn = .TRUE.
   !-------------------------------------------------------------------------------------------------
   ! Set up AD variables
   !-------------------------------------------------------------------------------------------------

   p%LinearizeFlag     = .FALSE.             ! InitInp%LinearizeFlag
   p%Blade%BladeLength = InitInp%TurbineComponents%BladeLength
   p%DtAero            = Interval            ! set the default DT here; may be overwritten later, when we read the input file in AD_GetInput()
   p%UseDWM            = InitInp%UseDWM

         ! Define parameters here:
   
   !.................................................................
   ! Umass WInDS.....................................................
   ! Get the current time
   CALL DATE_AND_TIME ( Values=O%FVM_Other%StrtTime )                        ! Let's time the whole simulation
   CALL CPU_TIME ( O%FVM_Other%UsrTime1 )                                    ! Initial time (this zeros the start time when used as a MATLAB function)
   O%FVM_Other%UsrTime1 = MAX( 0.0_DbKi, O%FVM_Other%UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   ! Umass WInDS....................................................
   !.................................................................
   
   
   p%WrOptFile   = InitInp%WrSumFile

   p%NumBl   = SIZE( InitInp%TurbineComponents%Blade )
   IF ( p%NumBl < 1 ) THEN
      CALL SetErrStat( ErrID_Fatal,'AeroDyn cannot run without blades in the model.',ErrStat,ErrMess,RoutineName)
      RETURN
   END IF
!bjj: what's the difference between p%NumBl, p%Blade%NB, and InitInp%NumBl?
!MLB: Heck if I know!

         ! Define initial system states here:
   !-------------------------------------------------------------------------------------------------
   ! Read the AeroDyn input file and open the output file if requested
   ! bjj: these should perhaps be combined
   !-------------------------------------------------------------------------------------------------
   CALL AD_GetInput(InitInp, P, x, xd, z, O, y, ErrStatLcl, ErrMessLcl )
      CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
      IF (ErrStat >= AbortErrLev ) RETURN 

   p%WindFileName      = InitInp%WindFileName ! InitInp%WindFileName  gets set in AD_GetInput


      ! allocate variables for aerodyn forces
   p%LinearizeFlag     = .FALSE.

   Interval = p%DtAero
   
   !   IF ( ABS( InitInp%IfW_InitInputs%lidar,HubPosition(3) - p%Rotor%HH ) > 0.1*( InitInp%IfW_InitInputs%lidar,HubPosition(3) ) )  THEN  
   !
   !      CALL ProgWarn( ' The ElastoDyn hub height ('//TRIM(Num2LStr( InitInp%IfW_InitInputs%lidar,HubPosition(3) ))//') and AeroDyn input'// &
   !                    ' reference hub height ('//TRIM(Num2LStr(p%Rotor%HH))//') differ by more than 10%.' )
   !   ENDIF
   
   IF (.NOT. EqualRealNos( InitInp%IfW_InitInputs%lidar%HubPosition(3), p%Rotor%HH ) ) THEN
      call SetErrStat( ErrID_Warn,'ElastoDyn and AeroDyn have different hub heights',ErrStat,ErrMess,RoutineName)      
   END IF
   

   IF ( .NOT. ALLOCATED( o%StoredForces  )) THEN
      CALL AllocAry(O%StoredForces, 3,p%Element%NELM,p%NumBl,'O%StoredForces',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
   IF ( .NOT. ALLOCATED( o%StoredMoments ))  THEN
      CALL AllocAry(O%StoredMoments, 3,p%Element%NELM,p%NumBl,'O%StoredForces',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
     
   IF (.NOT. ALLOCATED(O%Element%W2) ) THEN
      CALL AllocAry(O%Element%W2, p%Element%NELM, p%NumBl,'O%Element%W2',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF

   IF (.NOT. ALLOCATED(O%Element%Alpha) ) THEN
      CALL AllocAry(O%Element%Alpha, p%Element%NELM, p%NumBl,'O%Element%Alpha',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
   IF (ErrStat >= AbortErrLev ) RETURN 

   
   P%UnWndOut = -1
   P%UnElem = -1   
   IF ( p%ElemPrn )  THEN
      CALL ElemOpen ( TRIM( InitInp%OutRootName )//'.AD.out', P, O, ErrStat, ErrMess, AD_Ver )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
         IF (ErrStat >= AbortErrLev ) RETURN 
   END IF
      

   !-------------------------------------------------------------------------------------------------
   ! Calculate the rotor and hub radaii from the input values
   !-------------------------------------------------------------------------------------------------
   HubRadius = DOT_PRODUCT( InitInp%TurbineComponents%Blade(1)%Position(:)        &
                          - InitInp%TurbineComponents%Hub%Position(:),            &
                            InitInp%TurbineComponents%Blade(1)%Orientation(3,:) )

   DO IB = 2,p%NumBl
      TmpVar    = DOT_PRODUCT( InitInp%TurbineComponents%Blade(IB)%Position(:)    &
                             - InitInp%TurbineComponents%Hub%Position(:),         &
                               InitInp%TurbineComponents%Blade(IB)%Orientation(3,:) )
      IF ( ABS( TmpVar - HubRadius ) > 0.001 ) THEN ! within 1 mm
         CALL ProgWarn( ' AeroDyn\AD_Init() calculated HubRadius is not the same for all '// &
                           'blades. Using value from blade 1.' )
         EXIT
      END IF
   END DO !IB

   TipRadius = InitInp%TurbineComponents%BladeLength + HubRadius

   CosPrecone = ASIN( DOT_PRODUCT( InitInp%TurbineComponents%Blade(1)%Orientation(3,:), &
                                   InitInp%TurbineComponents%Hub%Orientation(1,:) ) )  ! precone angle -- do COS later

   DO IB = 2,p%NumBl
      TmpVar  = ASIN( DOT_PRODUCT( InitInp%TurbineComponents%Blade(IB)%Orientation(3,:), &
                                   InitInp%TurbineComponents%Hub%Orientation(1,:) ) )
      IF ( ABS( TmpVar - CosPrecone ) > 0.009 ) THEN     ! within ~ 1/2 degree
         CALL ProgWarn( ' AeroDyn\AD_Init() calculated precone angle is not the same for all'// &
                           ' blades. Using value from blade 1.' )
         EXIT
      END IF
   END DO !IBld

   CosPrecone = COS( CosPrecone )

   p%Blade%R = TipRadius * CosPrecone
   RHub = HubRadius * CosPrecone
   p%HubRad = RHub

      ! Check that the AeroDyn input DR and RElm match (use the HubRadius and TipRadius to verify)
      ! before using them to calculate the tip- and hub-loss constants
   CALL CheckRComp( P, x, xd, z, O, y, ErrStat, ErrMess, &
                    InitInp%ADFileName, HubRadius, TipRadius )

   IF ( ErrStat /= ErrID_None ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Calculate tip-loss constants
   !-------------------------------------------------------------------------------------------------
   DO IElm = 1,p%Element%NElm  ! Loop through all blade elements

      ElemRad = p%Element%RELM(IElm)*CosPrecone

      IF( ElemRad == 0.0 )  THEN  !BJJ: should this be 0.001 (or another small number) instead of exactly 0.0?
         CALL SetErrStat( ErrID_Fatal,'Error calculating tip loss constant for element '//TRIM(Int2LStr(IElm))//&
                          '. Division by zero.',ErrStat,ErrMess,RoutineName)
         
         RETURN
      ELSE
         DTip         = p%Blade%R - ElemRad
         p%Element%TLCNST(IElm) = 0.5 * p%NumBl * DTip / ElemRad
      ENDIF

   ENDDO             ! IElm - all blade elements


   !-------------------------------------------------------------------------------------------------
   ! Calculate hub-loss constants
   !-------------------------------------------------------------------------------------------------
   IF ( RHub > 0.001 )  THEN

      DO Ielm = 1,p%Element%NELM  ! Loop through all blade elements

         ElemRad = p%Element%RELM(Ielm)*CosPrecone  ! Use only the precone angle of blade 1 (assumed very similar to other blades)

         DHub         = ElemRad - RHub
         p%Element%HLCNST(Ielm) = 0.5 * p%NumBl * DHub / RHub

      ENDDO             ! IELM - all blade elements

   ELSE

      p%Element%HLCNST(:) = 0.0

   ENDIF



      !-------------------------------------------------------------------------------------------------
      ! Interpolate the tower diameter at ElastoDyn's tower nodes if we will be computing tower aerodynamics.
      !-------------------------------------------------------------------------------------------------

   IF ( p%TwrProps%CalcTwrAero )  THEN

         !-------------------------------------------------------------------------------------------------
         ! IMPORTANT NOTES:
         !     o  Supposedly, the glue code will not try to do anything with the tower-aero mesh if is is
         !        not created, so the creation is inside the test for CalcTwrAero.
         !     o  The tower properties from AeroDyn's tower file are for heights from the origin (ground or
         !        MSL) to the hub height--not the top of the tower.
         !     o  For now, we are allowing only one set of Cd for the entire tower.
         !     o  InterpIndx is initialize to 1 at compile time.
         !-------------------------------------------------------------------------------------------------


         ! Create the mesh for the tower aerodynamics.

      CALL MeshCreate ( BlankMesh       = u%Twr_InputMarkers      &
                      , IOS             = COMPONENT_INPUT         &
                      , NNodes          = InitInp%NumTwrNodes     &
                      , Orientation     = .TRUE.                  &
                      , TranslationDisp = .TRUE.                  &
                      , TranslationVel  = .TRUE.                  &
                      , ErrStat         = ErrStatLcl              &
                      , ErrMess         = ErrMessLcl              )

      CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
      IF ( ErrStat >= AbortErrLev )  RETURN


         ! Set the positions of the nodes.  MeshCreate() allocated the Position array.

      DO Node = 1,u%Twr_InputMarkers%Nnodes
         CALL MeshPositionNode ( Mesh  = u%Twr_InputMarkers          &
                                ,INode = Node                        &
                                ,Pos   = InitInp%TwrNodeLocs(:,Node) &  
                                ,ErrStat   = ErrStatLcl              &
                                ,ErrMess   = ErrMessLcl              )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

      END DO         


         ! Construct the tower with Line-2 elements.

      DO Elem=1,u%Twr_InputMarkers%Nnodes-1

         CALL MeshConstructElement ( Mesh     = u%Twr_InputMarkers &
                                   , Xelement = ELEMENT_LINE2      &
                                   , P1       = Elem               &
                                   , P2       = Elem+1             &
                                   , ErrStat  = ErrStatLcl         &
                                   , ErrMess  = ErrMessLcl         )

         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

      ENDDO


         ! Commit the mesh to the funny farm.

      CALL MeshCommit ( u%Twr_InputMarkers, ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN


         ! Copy the input mesh to create the output mesh.  Does

      CALL MeshCopy ( SrcMesh  = u%Twr_InputMarkers &
                    , DestMesh = y%Twr_OutputLoads  &
                    , CtrlCode = MESH_SIBLING       &
                    , Force    = .TRUE.             &
                    , ErrStat  = ErrStatLcl         &
                    , ErrMess  = ErrMessLcl         )

         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN


         ! Check to ensure that the user did not specify more than one set of Cd(Re) tables.  Temporary restriction.

      IF ( p%TwrProps%NTwrCD /= 1 )  THEN
         CALL SetErrStat(ErrID_Fatal,'You must have one and only one set of drag coefficients for the AeroDyn tower file.',ErrStat,ErrMess,RoutineName )
         RETURN
      END IF


         ! Build the TwrNodeWidth array.

      p%TwrProps%NumTwrNodes = InitInp%NumTwrNodes

      IF (.NOT. ALLOCATED( p%TwrProps%TwrNodeWidth ) ) THEN
         CALL AllocAry( p%TwrProps%TwrNodeWidth, p%TwrProps%NumTwrNodes, "array for tower widths at ED node locations", ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN
      END IF

      DO Node=1,p%TwrProps%NumTwrNodes

         TwrNodeHt = InitInp%TwrNodeLocs(3,Node)/p%Rotor%HH

         p%TwrProps%TwrNodeWidth(Node) = InterpStp( TwrNodeHt, p%TwrProps%TwrHtFr, p%TwrProps%TwrWid, InterpIndx, p%TwrProps%NTwrHT )

      END DO ! Node

   END IF ! ( p%TwrProps%CalcTwrAero )


   !-------------------------------------------------------------------------------------------------
   ! Write the summary (opt) file, then close it
   !-------------------------------------------------------------------------------------------------

   IF (p%WrOptFile) THEN

      CALL ADOut(InitInp, P, O, AD_Ver, TRIM(InitInp%OutRootName)//'.AD.sum', ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Initialize the wind inflow module
   !-------------------------------------------------------------------------------------------------

   InitInp%IfW_InitInputs%WindFileName = p%WindFileName
   InitInp%IfW_InitInputs%ReferenceHeight = p%Rotor%HH
   InitInp%IfW_InitInputs%Width = 2 * p%Blade%R
   InitInp%IfW_InitInputs%WindFileType = DEFAULT_WindNumber

!bjj: make sure p%Element%NElm, p%NumBl, and u%Twr_InputMarkers%Nnodes are set first

   CALL IfW_Init( InitInp%IfW_InitInputs,   o%IfW_Inputs,    p%IfW_Params,                          &
                     x%IfW_ContStates, xd%IfW_DiscStates,   z%IfW_ConstrStates,    O%IfW_OtherStates,   &
                     y%IfW_Outputs,    Interval,  InitOut%IfW_InitOutput,   ErrStatLcl,    ErrMessLcl )

   CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
   IF (ErrStat >= AbortErrLev) RETURN
   
   
         ! Allocate the InflowWind derived types.  This is OK for now because InflowWind is being used as a submodule of AeroDyn.
   IF (.NOT. ALLOCATED(o%IfW_Inputs%Position) ) THEN      
      CALL AllocAry( o%IfW_Inputs%Position, 3, p%Element%NElm*p%NumBl + u%Twr_InputMarkers%Nnodes + 1, "position vectors to be sent to IfW_CalcOutput", ErrStatLcl, ErrMessLcl )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
   END IF     
   
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Calling the DWM
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   IF ( p%UseDWM ) THEN   
      InitInp%DWM_InitInputs%IfW_InitInputs%WindFileName    = p%WindFileName
      InitInp%DWM_InitInputs%IfW_InitInputs%ReferenceHeight = p%Rotor%HH
      InitInp%DWM_InitInputs%IfW_InitInputs%Width           = 2 * p%Blade%R
      InitInp%DWM_InitInputs%IfW_InitInputs%WindFileType    = DEFAULT_WindNumber !bjj: I also set this in DWM for future when we don't use InflowWind in AeroDyn
   
         ! bjj: all this stuff should be put in DWM_Init.....>
      p%DWM_Params%RR              = p%Blade%R
      p%DWM_Params%BNum            = p%NumBl
      p%DWM_Params%ElementNum      = O%ElOut%NumElOut  !bjj: NumElOut is the number of elements to be printed in an output file. I really think you want the number of blade elements. I guess we should check that NumElOut is the same as p%Element%NElm
      p%DWM_Params%air_density     = p%Wind%Rho
   
      IF (.NOT. ALLOCATED(o%DWM_otherstates%Nforce  )) ALLOCATE ( o%DWM_otherstates%Nforce(  p%Element%NElm,p%NumBl),STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM Nforce array', ErrStat,ErrMess,RoutineName )
      IF (.NOT. ALLOCATED(o%DWM_otherstates%blade_dr)) ALLOCATE ( o%DWM_otherstates%blade_dr(p%Element%NElm),        STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM blade_dr array', ErrStat,ErrMess,RoutineName )
      IF (.NOT. ALLOCATED(p%DWM_Params%ElementRad   )) ALLOCATE ( p%DWM_Params%ElementRad(   p%Element%NElm),        STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM ElementRad array', ErrStat,ErrMess,RoutineName )
     
      o%DWM_otherstates%blade_dr(:) = p%Blade%DR(:)
      p%DWM_Params%ElementRad(:)    = p%Element%RELM(:)   
   
      CALL DWM_Init( InitInp%DWM_InitInputs, O%DWM_Inputs, p%DWM_Params, x%DWM_ContStates, xd%DWM_DiscStates, z%DWM_ConstrStates, & 
                     O%DWM_OtherStates, y%DWM_Outputs, Interval, InitOut%DWM_InitOutput, ErrStatLcl, ErrMessLcl )
   

      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
   END IF !UseDWM
   
   !-------------------------------------------------------------------------------------------------
   ! Turn off dynamic inflow for wind less than 8 m/s (per DJL: 8 m/s is really just an empirical guess)
   ! DJL: Comment out this code when using new proposed GDW check in ELEMFRC
   ! BJJ: FIX THIS!!!!
   !-------------------------------------------------------------------------------------------------

   IF (p%DynInfl) THEN

      IF ( p%IfW_Params%WindFileType /= FF_WINDNumber ) THEN
         CALL SetErrStat(ErrID_Info,'Dynamic inflow will not check for low mean wind speed.',ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN
      ELSE IF ( O%IfW_OtherStates%FFWind%MeanFFWS  < 8.0 ) THEN
         p%DynInfl = .FALSE.
         CALL SetErrStat(ErrID_Info,'Estimated average wind speed in FF wind file is less than 8 m/s. Dynamic Inflow will be turned off.',ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN
      END IF

   ENDIF

   !-------------------------------------------------------------------------------------------------
   ! Set initial guesses for inputs:
   !-------------------------------------------------------------------------------------------------
   
   !..........
   ! u%TurbineComponents
   !..........

   CALL AD_CopyAeroConfig( InitInp%TurbineComponents, u%TurbineComponents, MESH_NEWCOPY, ErrStatLcl, ErrMessLcl )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
   
   !..........
   ! u%InputMarkers (blade meshes):
   !..........
   
   ALLOCATE( u%InputMarkers(p%NumBl), STAT=ErrStatLcl )
   IF (ErrStatLcl /= 0 ) THEN
      CALL SetErrStat ( ErrID_Fatal, 'Could not allocate u%InputMarkers (meshes)', ErrStat,ErrMess,RoutineName )
      RETURN
   END IF


   DO IB = 1, p%NumBl
      CALL MeshCreate( BlankMesh      = u%InputMarkers(IB)    &
                     ,IOS            = COMPONENT_INPUT        &
                     ,NNodes         = p%Element%NELM         &
                     ,Orientation    = .TRUE.                 &
                     ,TranslationVel = .TRUE.                 &
                     ,TranslationAcc = .TRUE.                 &  !bjj: added for MHK turbines
                     ,RotationVel    = .TRUE.                 &
                     ,nScalars       = 2                      &  ! scalar 1 is W, scalar 2 is Alpha
                     ,ErrStat        = ErrStatLcl             &
                     ,ErrMess        = ErrMessLcl             )

      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
      ! create the elements
      DO IE = 1, p%Element%NELM-1 ! construct the blades into Line2 elements
         CALL MeshConstructElement ( Mesh = u%InputMarkers(IB)    &
                                  ,Xelement = ELEMENT_LINE2       &
                                  ,P1       = IE                  &
                                  ,P2       = IE+1                &
                                  ,ErrStat  = ErrStatLcl          &
                                  ,ErrMess  = ErrMessLcl          )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         
      ENDDO
     
      ! position/orient the nodes
      DO IE = 1, p%Element%NELM
         TmpPos(1) = 0.
         TmpPos(2) = 0.
         TmpPos(3) = p%Element%Relm(IE) - HubRadius
         CALL MeshPositionNode ( Mesh = u%InputMarkers(IB)              &
                                 ,INode = IE                            &
                                 ,Pos= TmpPos                           &  ! this info comes from FAST (not yet)
                                 ,ErrStat   = ErrStatLcl                &
                                 ,ErrMess   = ErrMessLcl                )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

         ! RELATIVE ORIENTATION OF BLADE ELEMENTS
         u%InputMarkers(IB)%Orientation(1,1,IE) = COS( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(2,1,IE) = SIN( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(3,1,IE) = SIN( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(1,2,IE) = -1. * u%InputMarkers(IB)%Orientation(2,1,IE)
         u%InputMarkers(IB)%Orientation(2,2,IE) =       u%InputMarkers(IB)%Orientation(1,1,IE)
         u%InputMarkers(IB)%Orientation(3,2,IE) = 0.0
         u%InputMarkers(IB)%Orientation(1,3,IE) = 0.0
         u%InputMarkers(IB)%Orientation(2,3,IE) = 0.0
         u%InputMarkers(IB)%Orientation(3,3,IE) = 1.0
      ENDDO
     
       CALL MeshCommit ( Mesh = u%InputMarkers(IB)    &
                        ,ErrStat  = ErrStatLcl        &
                        ,ErrMess  = ErrMessLcl        )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
          
   ENDDO
         

   !..........
   ! u%Twr_InputMarkers (tower meshes):
   !..........
     
   !bjj: done above in section for IF Tower Loads is on
   

   !..........
   ! u%MulTabLoc:
   !..........
   
   IF (.NOT. ALLOCATED(u%MulTabLoc)) THEN
      ALLOCATE( u%MulTabLoc(p%Element%NELM, p%NumBl), STAT = ErrStatLcl )
      IF (ErrStatLcl /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Could not allocate u%MulTabLoc', ErrStat,ErrMess,RoutineName )
         RETURN
      END IF
   END IF

   u%MulTabLoc(:,:) = 0.0
   
   
   !-------------------------------------------------------------------------------------------------
   ! Allocate space for outputs and set up output meshes:
   !-------------------------------------------------------------------------------------------------
   
   !..........
   ! y%OutputLoads (blade meshes):
   !..........
   
   
   ALLOCATE( y%OutputLoads(p%NumBl), STAT = ErrStatLcl )
      IF (ErrStatLcl /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Could not allocate y%OutputLoads (meshes)', ErrStat,ErrMess,RoutineName )
         RETURN
      END IF
   
   DO IB = 1, p%NumBl

       CALL MeshCopy ( SrcMesh  = u%InputMarkers(IB)  &
                      ,DestMesh = y%OutputLoads(IB)   &
                      ,CtrlCode = MESH_SIBLING        &
                      ,Force    = .TRUE.              &
                      ,Moment   = .TRUE.              &
                      ,ErrStat  = ErrStatLcl          &
                      ,ErrMess  = ErrMessLcl          )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
   ENDDO
   
   
   !..........
   ! y%Twr_OutputLoads (tower meshes):
   !..........

   !bjj: done above in section for IF Tower Loads is on
   
     
   !-------------------------------------------------------------------------------------------------
   ! Initialize AeroDyn variables not initialized elsewhere (except in module initialization)
   ! and return
   !-------------------------------------------------------------------------------------------------
   o%InducedVel%SumInfl  = 0.0_ReKi
   o%Rotor%AvgInfl       = 0.0_ReKi
   o%OldTime             = 0.0_DbKi
   O%SuperSonic          = .FALSE.   
   o%NoLoadsCalculated   = .TRUE.

   p%TwoPiNB     = TwoPi / REAL( p%NumBl, ReKi )

   p%Initialized = .TRUE.
   
        
   
   DO ie = 1, maxInfl
      p%DynInflow%xMinv(ie) = PIBY2 / hfunc(MRvector(ie), NJvector(ie))   !bjj: this is really just a Fortran parameter, too.
   END DO !ie   

    
   !************************************************************************************************************   
   !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...   
   !............................................................................................................
   ! Umass WInDS starts
   !------------------------------------
   
   p%Blade%TipRadius      =    TipRadius         ! sliu: maybe do this: TipRadius = InitInp%TurbineComponents%BladeLength + HubRadius
   p%Blade%HubRadius      =    HubRadius         ! = DOT_PRODUCT( InitInp%TurbineComponents%Blade(1)%Position(:) - InitInp%TurbineComponents%Hub%Position(:),  InitInp%TurbineComponents%Blade(1)%Orientation(3,:) ) 

   p%FVM%UseWINDS         =    .TRUE.            ! whether to use WINDS
   O%Aerodyn_Timestep     =    -1                ! The timestep used in FAST and AeroDyn (Same as n_t_global in FAST_Prog.f90)
                                                 ! First call at FAST_Solution0. Second call at FAST_Solution and n_t_global = 0                                                 
   
   O%WINDS_Timestep       =    1                 ! The timestep in  WINDS. Start from 1.  (O%WINDS_Timestep - 1) * Dt_Ratio = O%Aerodyn_Timestep

  
   O%FVM_Other%TIME%Time_Total  = 0.0            ! Total time of WINDS
   O%FVM_Other%TIME%Time_biotsavart  = 0.0       ! Total time of Inducevelocity subroutine
   O%FVM_Other%TIME%Time_biotsavart_Acce  = 0.0 ! Total time of biotsavart subroutine
   


   ! Read input file
   CALL WInDS_ReadInput(InitInp, P, xd, O, ErrStat, ErrMess )
   IF (ErrStat /= 0 ) RETURN

      ! Set parameters
   CALL WINDS_SetParameters( p, O, ErrStat, ErrMess )
   IF (ErrStat /= 0 ) RETURN
   
     ! Allocate valuables
   CALL WINDS_Allocate( p, O, xd, ErrStat, ErrMess )
   IF (ErrStat /= 0 ) RETURN
   
     ! Wind shear 
   IF (p%FVM%Shear_Parms%ShearFLAG ) THEN
      CALL WINDS_Shear_Model(p, O, ErrStat, ErrMess)
      IF (ErrStat /= 0 ) RETURN
   END IF
         
     ! Ground effects
   IF (p%FVM%Ground_Parms%GroundFLAG  .AND.  p%FVM%Ground_Parms%METHOD == 'PANEL') THEN
      CALL WINDS_Ground_model(p, O, ErrStat, ErrMess)
      IF (ErrStat /= 0 ) RETURN
   END IF   
      
   ! Calculate varibles for LB dynamic stall
   IF (p%FVM%DS_Parms%DS_Flag) THEN
       
      ! sliu: do not load data any more .. 
      IF (p%FVM%DS_Parms%load_data) THEN 
         CALL LB_load_AirfoilData(p, O, xd, ErrStat, ErrMess) 
      ELSE
         CALL LB_Initialize_AirfoilData(p, O, xd, ErrStat, ErrMess)
      END IF      
      p%FVM%DS_Parms%start_n = FLOOR( p%FVM%DS_Parms%start_t / p%FVM%DT_WINDS + 1 )
      ! Write these variables into txt file for debug purpose
      IF (p%FVM%DS_Parms%write_data) THEN
         CALL Write_DS_parameters(p, O, ErrStat, ErrMess)
      END IF
      IF ( ErrStat /= 0 )  RETURN   
   ELSE 
      p%FVM%DS_Parms%start_n = p%FVM%NT + 1
               
   END IF      
   
      ! Initialize paraview animation files (.pvd) 
   IF (p%FVM%AnimFLAG) THEN
      CALL Initialize_paraview_files(P, xd, O, ErrStat, ErrMess )
      IF (ErrStat /= 0 ) RETURN
   END IF
   
       ! Record the speedup and error of N-body algorithm or parallel computation
   IF (p%FVM%Tree_Parms%Speedup) THEN   
       CALL WRITE_Treecode(0_IntKi, 0_IntKi, 0.0_DbKi,  0.0_DbKi, 0.0_DbKi, 0.0_DbKi, 'START', p)
   END IF
   
      ! Record the iteration number of KJ
   IF (p%FVM%KJ_output) THEN
      CALL WRITE_KJ(0_IntKi, 0_IntKi, 0.0_DbKi,  'START', p)
   END IF   
   

       
   
   !------------------------------------
   ! Umass WInDS ends
   !............................................................................................................
   !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...   
   !************************************************************************************************************     
   
  
   
   
   
   
   RETURN

END SUBROUTINE AD_Init

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMess )
! This routine is called at the end of the simulation.
!..................................................................................................................................
      USE DWM_Types
      USE DWM

      TYPE(AD_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(AD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(AD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(AD_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""



      !************************************************************************************************************
      !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...
      !............................................................................................................
      ! Umass WInDS 
      
         ! Close the Paraview input file
      IF (p%FVM%AnimFLAG) THEN
         CALL Close_paraview_files(P, xd, OtherState, ErrStat, ErrMess )
         IF (ErrStat /= 0 ) RETURN
      END IF
      
         ! Record the speedup and error of N-body algorithm or parallel computation
      IF (p%FVM%Tree_Parms%Speedup) THEN   
          CALL WRITE_Treecode(0_IntKi, 0_IntKi, 0.0_DbKi,  0.0_DbKi, 0.0_DbKi, 0.0_DbKi, 'END', p)
      END IF
      
      ! Record the iteration number of KJ
      IF (p%FVM%KJ_output) THEN
         CALL WRITE_KJ(0_IntKi, 0_IntKi, 0.0_DbKi, 'END', p)
      END IF 
      
      ! Summary file
      IF (p%FVM%WINDS_Sum) THEN      
         CALL WInDS_WriteSum(P, xd, OtherState, ErrStat, ErrMess )
      END IF 
      !............................................................................................................
      !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...
      !************************************************************************************************************
      
      
      
      
         ! Place any last minute operations or calculations here:
         
      IF (p%UseDWM ) THEN
         !----- Call the DWM ------- 
      
         CALL DWM_End( OtherState%DWM_Inputs, p%DWM_Params, x%DWM_ContStates, xd%DWM_DiscStates, z%DWM_ConstrStates, &
                                           OtherState%DWM_OtherStates, y%DWM_Outputs, ErrStat, ErrMess )
      END IF ! UseDWM
      
      !--------------------------

      CALL IfW_End(  OtherState%IfW_Inputs, p%IfW_Params, x%IfW_ContStates, xd%IfW_DiscStates, z%IfW_ConstrStates, &
                     OtherState%IfW_OtherStates, y%IfW_Outputs, ErrStat, ErrMess )

         ! Close files here:

      ! AD_IOParams
   IF (P%UnEc > 0)    CLOSE(P%UnEc) ! not currently used

   IF (P%UnWndOut > 0) CLOSE(P%UnWndOut)
   IF (P%UnElem   > 0) CLOSE(P%UnElem)

         ! Destroy the input data:

      CALL AD_DestroyInput( u, ErrStat, ErrMess )


         ! Destroy the parameter data:

      CALL AD_DestroyParam( p, ErrStat, ErrMess )


         ! Destroy the state data:

      CALL AD_DestroyContState(   x,           ErrStat, ErrMess )
      CALL AD_DestroyDiscState(   xd,          ErrStat, ErrMess )
      CALL AD_DestroyConstrState( z,           ErrStat, ErrMess )
      CALL AD_DestroyOtherState(  OtherState,  ErrStat, ErrMess )


         ! Destroy the output data:

      CALL AD_DestroyOutput( y, ErrStat, ErrMess )




END SUBROUTINE AD_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMess )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n           ! Current simulation time step n = 0,1,...
      TYPE(AD_InputType),                 INTENT(INOUT) :: u(:)        ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)   ! Times associated with u(:), in seconds
      TYPE(AD_ParameterType),             INTENT(IN   ) :: p           ! Parameters
      TYPE(AD_ContinuousStateType),       INTENT(INOUT) :: x           ! Input: Continuous states at t;
                                                                       !   Output: Continuous states at t + Interval
      TYPE(AD_DiscreteStateType),         INTENT(INOUT) :: xd          ! Input: Discrete states at t;
                                                                       !   Output: Discrete states at t  + Interval
      TYPE(AD_ConstraintStateType),       INTENT(INOUT) :: z           ! Input: Initial guess of constraint states at t;
                                                                       !   Output: Constraint states at t
      TYPE(AD_OtherStateType),            INTENT(INOUT) :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat     ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMess     ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(AD_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at Time
      TYPE(AD_ConstraintStateType)                 :: z_Residual  ! Residual of the constraint state equations (Z)

      INTEGER(IntKi)                                    :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMess))                            :: ErrMess2     ! Error message if ErrStat2 /= ErrID_None

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""





END SUBROUTINE AD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_CalcOutput( Time, u, p, x, xd, z, O, y, ErrStat, ErrMess )
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
      USE               AeroGenSubs,   ONLY: ElemOut
      USE               DWM_Types
      USE               DWM

      REAL(DbKi),                   INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O!therState ! Other/optimization states
      TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


      ! Local variables
   REAL(DbKi), PARAMETER      :: OnePlusEpsilon = 1 + EPSILON(Time)

   REAL(ReKi)                 :: VNElement
   REAL(ReKi)                 :: VelNormalToRotor2
   REAL(ReKi)                 :: VNWind
   REAL(ReKi)                 :: VTTotal
   REAL(ReKi)                 :: DFN
   REAL(ReKi)                 :: DFT
   REAL(ReKi)                 :: PMA
   REAL(ReKi)                 :: SPitch                     ! sine of PitNow
   REAL(ReKi)                 :: CPitch                     ! cosine of PitNow

   REAL(ReKi)                 :: AvgVelNacelleRotorFurlYaw
   REAL(ReKi)                 :: AvgVelTowerBaseNacelleYaw
   REAL(ReKi)                 :: AvgVelTowerBaseYaw
   REAL(ReKi)                 :: AzimuthAngle
   REAL(ReKi)                 :: rNacelleHub   (2)
   REAL(ReKi)                 :: rLocal
   REAL(ReKi)                 :: rRotorFurlHub (2)
   REAL(ReKi)                 :: rTowerBaseHub (2)

   REAL(ReKi)                 :: tmpVector     (3)
   REAL(ReKi)                 :: VelocityVec   (3)

   INTEGER                    :: ErrStatLcL        ! Error status returned by called routines.
   INTEGER                    :: IBlade
   INTEGER                    :: IElement
   INTEGER                    :: Node              ! Node index.

   INTEGER                    :: I
   CHARACTER(LEN(ErrMess))    :: ErrMessLcl          ! Error message returned by called routines.

   !!.....................................................
   !! Umass WINDS
   !REAL(DbKi)          :: Time_1  ! To record CPU time(Start time) 
   !REAL(DbKi)          :: Time_2  ! To record CPU time(End time)  
   !! Umass WINDS    
   !!.....................................................   
   
   !..................................... ...............     
   ! umass debug, to be deleted
   REAL(DbKi), DIMENSION(p%Element%NElm, p%NumBl)               :: PRINT_NAME1
   REAL(DbKi), DIMENSION(p%Element%NElm+1, p%NumBl)             :: PRINT_NAME2
   REAL(DbKi), DIMENSION(p%Element%NElm, 3)                     :: PRINT_NAME3   
   
   INTEGER(IntKi)      :: IDim        ! For dimensions
   INTEGER(IntKi)      :: ITimestep   ! For timesteps
   
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB 
   
   
   CHARACTER(LEN=2), DIMENSION(20)      :: Rank_num
   CHARACTER(LEN=1024)                  :: temp_number
   CHARACTER(LEN=1)                     :: temp_1 
   CHARACTER(LEN=2)                     :: temp_2 
   CHARACTER(LEN=3)                     :: temp_3 
      
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl     
   ! umass debug, to be deleted
   !....................................................           
   
   


   ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMess  = ""

   !-------------------------------------------------------------------------------------------------
   ! Check that the module has been initialized.
   !-------------------------------------------------------------------------------------------------
   IF ( .NOT. p%Initialized ) THEN
      ErrStat = ErrID_Fatal
      ErrMess = 'AeroDyn must be initialized before trying to calculate aerodynamic loads.'
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Determine if loads should be recalculated or just returned
   !-------------------------------------------------------------------------------------------------
      ! NOTE: Time is scaled by OnePlusEps to ensure that loads are calculated at every
   !       time step when DTAero = DT, even in the presence of numerical precision errors.

   IF ( o%NoLoadsCalculated .OR. ( Time*OnePlusEpsilon - o%OldTime ) >= p%DTAERO )  THEN
         ! It's time to update the aero forces

         ! First we reset the DTAERO parameters for next time
      o%DT      = Time - o%OldTime     !bjj: DT = 0 on first step,
                                       !but the subroutines that use DT check for NoLoadsCalculated (or time > 0)
      o%OldTime = Time

   ELSE IF ( .NOT. p%LinearizeFlag ) THEN

         ! Return the previously-calculated loads

!      CurrentOutputs = ADCurrentLoads

      DO IBlade=1,p%NumBl
       DO IElement=1,p%Element%Nelm
         y%OutputLoads(IBlade)%Force(:,IElement)  = o%StoredForces(:,IElement,IBlade)
         y%OutputLoads(IBlade)%Moment(:,IElement) = o%StoredMoments(:,IElement,IBlade)
       ENDDO
      ENDDO

      IF ( O%FirstWarn ) THEN
         CALL SetErrStat ( ErrID_Warn, 'AeroDyn was designed for an explicit-loose coupling scheme. '//&
            'Using last calculated values from AeroDyn on all subsequent calls until time is advanced. '//&
            'Warning will not be displayed again.', ErrStat,ErrMess,'AD_CalcOutput' )
         O%FirstWarn = .FALSE.       
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
      END IF
            
      RETURN

   ENDIF

   
      ! Fill input array for InflowWind
      
   Node = 0
   DO IBlade = 1,p%NumBl
      DO IElement = 1,p%Element%NElm      
         Node = Node + 1
         o%IfW_Inputs%Position(:,Node) = u%InputMarkers(IBlade)%Position(:,IElement)
      END DO
   END DO
   
   DO Node=1,u%Twr_InputMarkers%Nnodes
      o%IfW_Inputs%Position(:,p%Element%NElm*p%NumBl + Node) = u%Twr_InputMarkers%TranslationDisp(:,Node) + u%Twr_InputMarkers%Position(:,Node)
   END DO      
   o%IfW_Inputs%Position(:,p%Element%NElm*p%NumBl+u%Twr_InputMarkers%Nnodes+1) = (/0.0_ReKi, 0.0_ReKi, p%Rotor%HH /)

      ! Calculate the wind-speed inputs
  
   CALL IfW_CalcOutput( Time, o%IfW_Inputs, p%IfW_Params, &
                           x%IfW_ContStates, xd%IfW_DiscStates, z%IfW_ConstrStates, O%IfW_OtherStates, &   ! States -- none in this case
                           y%IfW_Outputs, ErrStatLcl, ErrMessLcl )

   CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   END IF

   
      ! set outputs from InflowWind
      
   if (.not. allocated( y%IfW_Outputs%WriteOutput ) ) then
      CALL AllocAry( y%IfW_Outputs%WriteOutput,3 + min(5,p%IfW_Params%lidar%NumPulseGate), "y%IfW_Outputs%WriteOutput", ErrStatLcl, ErrMessLcl )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF   
   end if

   y%IfW_Outputs%WriteOutput(1:3) = y%IfW_Outputs%Velocity(:,p%Element%NElm*p%NumBl+u%Twr_InputMarkers%Nnodes+1)
   do i=1,min(5,p%IfW_Params%lidar%NumPulseGate)
      y%IfW_Outputs%WriteOutput(3+i) = y%IfW_Outputs%lidar%lidSpeed(i)
   end do
   

   !-------------------------------------------------------------------------------------------------
   ! Calculate the forces and moments for the blade: SUBROUTINE AeroFrcIntrface( FirstLoop, JElemt, DFN, DFT, PMA )
   !-------------------------------------------------------------------------------------------------

      ! calculate rotor speed
      ! note: Subtracting the RotorFurl rotational velocity for REVS is needed to get the
      ! same answers as before v13.00.00. RotorFurl shouldn't be needed.

   o%Rotor%REVS = ABS( DOT_PRODUCT( u%TurbineComponents%Hub%RotationVel(:) - u%TurbineComponents%RotorFurl%RotationVel(:), &
                                    u%TurbineComponents%Hub%Orientation(1,:) ) )


      ! calculate yaw angle
      ! note: YawAng should use the Hub instead of the RotorFurl, but it is calculated this way to
      ! get the same answers as previous version.
   o%Rotor%YawAng = ATAN2( -1.*u%TurbineComponents%RotorFurl%Orientation(1,2), u%TurbineComponents%RotorFurl%Orientation(1,1) )
   o%Rotor%SYaw   = SIN( o%Rotor%YawAng )
   o%Rotor%CYaw   = COS( o%Rotor%YawAng )

      ! tilt angle
      ! note: tilt angle should use the Hub instead of RotorFurl, but it needs hub to get the same
      ! answers as the version before v13.00.00

   o%Rotor%Tilt = ATAN2( u%TurbineComponents%RotorFurl%Orientation(1,3), &
                         SQRT( u%TurbineComponents%RotorFurl%Orientation(1,1)**2 + &
                         u%TurbineComponents%RotorFurl%Orientation(1,2)**2 ) )

   o%Rotor%CTilt     = COS( o%Rotor%Tilt )
   o%Rotor%STilt     = SIN( o%Rotor%Tilt )


      ! HubVDue2Yaw - yaw velocity due solely to yaw

  AvgVelNacelleRotorFurlYaw = u%TurbineComponents%RotorFurl%RotationVel(3) - u%TurbineComponents%Nacelle%RotationVel(3)
  AvgVelTowerBaseNacelleYaw = u%TurbineComponents%Nacelle%RotationVel(3)   - u%TurbineComponents%Tower%RotationVel(3)
  AvgVelTowerBaseYaw        = u%TurbineComponents%Tower%RotationVel(3)

  rRotorFurlHub(1:2)        = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%RotorFurl%Position(1:2)
  rNacelleHub(1:2)          = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%Nacelle%Position(1:2)
  rTowerBaseHub(1:2)        = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%Tower%Position(1:2)

  o%Rotor%YawVel =   ( AvgVelNacelleRotorFurlYaw * rRotorFurlHub(2) + AvgVelTowerBaseNacelleYaw * rNacelleHub(2) &
                         + AvgVelTowerBaseYaw * rTowerBaseHub(2) ) * o%Rotor%SYaw &
                  - ( AvgVelNacelleRotorFurlYaw * rRotorFurlHub(1) + AvgVelTowerBaseNacelleYaw * rNacelleHub(1) &
                         + AvgVelTowerBaseYaw * rTowerBaseHub(1) ) * o%Rotor%CYaw


   !.................................................................................................
   ! start of NewTime routine
   !.................................................................................................

   o%Rotor%AvgInfl = o%InducedVel%SumInfl * 2.0 / (p%Blade%R*p%Blade%R*p%NumBl)  ! Average inflow from the previous time step
   o%InducedVel%SumInfl = 0.0   ! reset to sum for the current time step

   CALL DiskVel(Time, P, O, ErrStatLcl, ErrMessLcl)  ! Get a sort of "Average velocity" - sets a bunch of stored variables...
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput/DiskVel' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
         
   IF ( P%DStall ) CALL BedUpdate( O )   ! that's an 'O' as in 'OtherState'

   ! Enter the dynamic inflow routines here

   IF ( p%Wake )  THEN
      CALL Inflow(Time, P, O, ErrStatLcl, ErrMessLcl)
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput/Inflow' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF   
   END IF      
   
       !bjj: perhaps we should send NoLoadsCalculated to initialize dynamic inflow [subroutine Infinit()]
       !bjj: instead of the check that time > 0...?

   !.................................................................................................
   ! end of NewTime routine
   !.................................................................................................

   Node = 0
   DO IBlade = 1,p%NumBl

         ! calculate the azimuth angle ( we add pi because AeroDyn defines 0 as pointing downward)
         ! note: the equation below should use TurbineComponents%Blade markers, but this is used to get the
         ! same answers as the previous version (before v13.00.00)

      AzimuthAngle = ATAN2( -1.*DOT_PRODUCT( u%TurbineComponents%Hub%Orientation(3,:),         &
                                             u%TurbineComponents%RotorFurl%Orientation(2,:) ), &
                                DOT_PRODUCT( u%TurbineComponents%Hub%Orientation(3,:),         &
                                             u%TurbineComponents%RotorFurl%Orientation(3,:) )  ) + pi + (IBlade - 1)*p%TwoPiNB



      DO IElement = 1,p%Element%NElm

            ! calculate element pitch

         o%Element%PitNow    = -1.*ATAN2( -1.*DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(2,:,IElement) ) , &
                                              DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(1,:,IElement) )   )

         SPitch    = SIN( o%Element%PitNow )
         CPitch    = COS( o%Element%PitNow )


            ! calculate distance between hub and element

         tmpVector = u%InputMarkers(IBlade)%Position(:,IElement) - u%TurbineComponents%Hub%Position(:)
         rLocal = SQRT(   DOT_PRODUCT( tmpVector, u%TurbineComponents%Hub%Orientation(2,:) )**2  &
                        + DOT_PRODUCT( tmpVector, u%TurbineComponents%Hub%Orientation(3,:) )**2  )

            ! determine if MulTabLoc should be set.  
     
         O%AirFoil%MulTabLoc = u%MulTabLoc(IElement,IBlade)
         
         !-------------------------------------------------------------------------------------------
         ! Get wind velocity components; calculate velocity normal to the rotor squared
         ! Save variables for printing in a file later;
         !-------------------------------------------------------------------------------------------
         Node = Node + 1
         VelocityVec(:)    = AD_WindVelocityWithDisturbance( Time, u, p, x, xd, z, O, y, ErrStatLcl, ErrMessLcl, &
                                                             u%InputMarkers(IBlade)%Position(:,IElement), &
                                                             y%IfW_Outputs%Velocity(:,Node) )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF
         
         !-------------------------------------------------------------------------------------------
         ! DWM wind input update phase 1
         !-------------------------------------------------------------------------------------------
         IF (p%UseDWM) THEN
            !bjj: FIX THIS!!!!         
            !bjj: where do p%DWM_Params%RTPD%SimulationOrder_index and p%DWM_Params%RTPD%upwindturbine_number get set?
         
            IF ( p%DWM_Params%RTPD%SimulationOrder_index > 1) THEN
               IF(  p%DWM_Params%RTPD%upwindturbine_number /= 0 ) THEN
                       
                  o%DWM_otherstates%position_y = u%InputMarkers(IBlade)%Position(2,IElement)
                                         
                  o%DWM_otherstates%position_z = u%InputMarkers(IBlade)%Position(3,IElement)
              
                  o%DWM_otherstates%velocity_wake_mean = 1
              
                  DO I = 1,p%DWM_Params%RTPD%upwindturbine_number
                     o%DWM_otherstates%DWM_tb%Aerodyn_turbine_num = I
                 
                     CALL   DWM_phase1( Time, O%DWM_Inputs, p%DWM_Params, x%DWM_contstates, xd%DWM_discstates, z%DWM_constrstates, &
                                           o%DWM_otherstates, y%DWM_outputs, ErrStatLcl, ErrMessLcl )
                 
                     CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput/DWM_phase1' )
                     IF (ErrStat >= AbortErrLev) THEN
                        CALL CleanUp()
                        RETURN
                     END IF   
                                                           
                     o%DWM_otherstates%velocity_wake_mean = (1-((1-o%DWM_otherstates%velocity_wake_mean)**2 + (1-o%DWM_otherstates%shifted_velocity_aerodyn)**2)**0.5)
                  END DO
              
                  o%DWM_otherstates%velocity_wake_mean    = o%DWM_otherstates%velocity_wake_mean * p%DWM_Params%Wind_file_Mean_u
              
                  VelocityVec(1) = (VelocityVec(1) - p%DWM_Params%Wind_file_Mean_u)*(O%DWM_Inputs%Upwind_result%upwind_small_TI(I)/p%DWM_Params%TI_amb) &
                                  + o%DWM_otherstates%velocity_wake_mean
              
               END IF
            END IF
                                     
           !------------------------DWM PHASE 2-----------------------------------------------
            IF (Time > 50.00 ) THEN
               o%DWM_otherstates%U_velocity           = VelocityVec(1)
               o%DWM_otherstates%V_velocity           = VelocityVec(2)
               o%DWM_otherstates%NacYaw               = o%Rotor%YawAng 
               o%DWM_otherstates%DWM_tb%Blade_index   = IBlade
               o%DWM_otherstates%DWM_tb%Element_index = IElement    

               CALL   DWM_phase2( Time, O%DWM_Inputs, p%DWM_Params, x%DWM_contstates, xd%DWM_discstates, z%DWM_constrstates, &
                                           o%DWM_otherstates, y%DWM_outputs, ErrStatLcl, ErrMessLcl )
            
                     CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput/DWM_phase1' )
                     IF (ErrStat >= AbortErrLev) THEN
                        CALL CleanUp()
                        RETURN
                     END IF   
            
         
               !CALL CalVelScale(VelocityVec(1),VelocityVec(2),O%o%DWM_otherstatesutputType,DWM_ConstraintStateType)
         
               !CALL turbine_average_velocity( VelocityVec(1), IBlade, IElement, O%o%DWM_otherstatesutputType,AD_ParameterType,DWM_ConstraintStateType)
            END IF
         END IF ! UseDWM
            
         !-----------------------------------------------------------------------------------------------------------------------
         
         
         
         
         VelNormalToRotor2 = ( VelocityVec(3) * o%Rotor%STilt + (VelocityVec(1) * o%Rotor%CYaw               &
                             - VelocityVec(2) * o%Rotor%SYaw) * o%Rotor%CTilt )**2

         !-------------------------------------------------------------------------------------------
         ! reproduce GetVNVT routine:
         !-------------------------------------------------------------------------------------------
         tmpVector =  -1.*SPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                        + CPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VTTotal   =     DOT_PRODUCT( tmpVector, VelocityVec - u%InputMarkers(IBlade)%TranslationVel(:,IElement)  )

         tmpVector =     CPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                       + SPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VNWind    =     DOT_PRODUCT( tmpVector, VelocityVec )
         VNElement = -1.*DOT_PRODUCT( tmpVector, u%InputMarkers(IBlade)%TranslationVel(:,IElement ) )

         !-------------------------------------------------------------------------------------------
         ! Get blade element forces and induced velocity
         !-------------------------------------------------------------------------------------------
         CALL ELEMFRC( p, O, ErrStatLcl, ErrMessLcl,                             &
                       AzimuthAngle, rLocal, IElement, IBlade, VelNormalToRotor2, VTTotal, VNWind, &
                       VNElement, DFN, DFT, PMA, o%NoLoadsCalculated )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF
            
         !-------------------------------------------------------------------------------------------
         ! Set up dynamic inflow parameters
         !-------------------------------------------------------------------------------------------
         IF ( p%DynInfl .OR. O%DynInit ) THEN
            CALL GetRM (P, O, ErrStatLcl, ErrMessLcl, &
                        rLocal, DFN, DFT, AzimuthAngle, IElement, IBlade)
               CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput' )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL CleanUp()
                  RETURN
               END IF
            
         ENDIF
         
         IF (p%UseDWM) THEN
             o%DWM_otherstates%Nforce(IElement,IBlade) = DFN              ! 12.4.2014 add by yh
         END IF ! UseDWM

         o%StoredForces(1,IElement,IBlade)  = ( DFN*CPitch + DFT*SPitch ) / p%Blade%DR(IElement)
         o%StoredForces(2,IElement,IBlade)  = ( DFN*SPitch - DFT*CPitch ) / p%Blade%DR(IElement)
         o%StoredForces(3,IElement,IBlade)  = 0.0

         o%StoredMoments(1,IElement,IBlade)  = 0.0
         o%StoredMoments(2,IElement,IBlade)  = 0.0
         o%StoredMoments(3,IElement,IBlade)  = PMA / p%Blade%DR(IElement)

!      DO IBlade=1,p%NumBl
!       DO IElement=1,p%Element%Nelm
!         y%OutputLoads(IBlade)%Force(:,IElement)  = o%StoredForces(:,IElement,IBlade)
!         y%OutputLoads(IBlade)%Moment(:,IElement) = o%StoredMoments(:,IElement,IBlade)
!       ENDDO
!!      ENDDO

            ! save velocities for output, if requested

         IF ( O%ElOut%WndElPrList(IElement) > 0 ) THEN
            O%ElOut%SaveVX( O%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(1)
            O%ElOut%SaveVY( O%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(2)
            O%ElOut%SaveVZ( O%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(3)
         ENDIF

         
         
         !************************************************************************************************************
         !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...
         !............................................................................................................
         ! Umass WInDS (sliu: I tried to call AD_WindVelocityWithDisturbance in my code, but the output wind speed in .outb get zero)
   
         ! Get the inflow wind speed
         IF (O%Aerodyn_Timestep == -1) THEN
            O%FVM_Other%WIND_INFTY(1, 1, 1, IElement, IBlade) = VelocityVec(1)
            O%FVM_Other%WIND_INFTY(2, 1, 1, IElement, IBlade) = VelocityVec(2)
            O%FVM_Other%WIND_INFTY(3, 1, 1, IElement, IBlade) = VelocityVec(3)       
         
            O%FVM_Other%WIND_INFTYM(1, 1, 1, IElement, IBlade) = SQRT(VelocityVec(1)**2 + VelocityVec(2)**2 + VelocityVec(3)**2 )    
         END IF 
         !............................................................................................................          
         !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...
         !************************************************************************************************************         
         
         
         

      END DO !IElement

      IF ( IBlade == 1 .AND. p%ElemPrn ) THEN
         O%ElOut%VXSAV  = VelocityVec(1)
         O%ElOut%VYSAV  = VelocityVec(2)
         O%ElOut%VZSAV  = VelocityVec(3)
      ENDIF


   END DO !IBlade

   O%NoLoadsCalculated = .FALSE.


   DO IBlade=1,p%NumBl
     DO IElement=1,p%Element%Nelm
       y%OutputLoads(IBlade)%Force(:,IElement)  = o%StoredForces(:,IElement,IBlade)
       y%OutputLoads(IBlade)%Moment(:,IElement) = o%StoredMoments(:,IElement,IBlade)
     ENDDO
   ENDDO
   
   
   !------------------------DWM PHASE 3-----------------------------------------------
   IF (p%UseDWM) THEN
   
      IF (Time > 50.00 ) THEN !BJJ: why is 50 hard-coded here and above???
            
         !o%DWM_otherstates%Nforce(:,:)    = o%DWM_otherstates%DFN_DWM(:,:) 
         CALL   DWM_phase3( Time, O%DWM_Inputs, p%DWM_Params, x%DWM_contstates, xd%DWM_discstates, z%DWM_constrstates, &
                                o%DWM_otherstates, y%DWM_outputs, ErrStatLcl, ErrMessLcl )
    
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD_CalcOutput/DWM_phase3' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF   
      
         !CALL filter_average_induction_factor( AD_ParameterType, DWM_ConstraintStateType, O%o%DWM_otherstatesutputType )
      END IF
   END IF !UseDWM
   
   !-----------------------------------------------------------------------------------

        

      ! Loop through all the tower nodes to calculate the aerodynamic loads on the tower if aerodynamics were requested.

   IF ( p%TwrProps%CalcTwrAero )  THEN

      DO Node=1,u%Twr_InputMarkers%Nnodes


            ! Calculate the aerodynamic load on this tower node: TwrAeroLoads ( p, Node, NodeDCMGbl, NodeVelGbl, NodeWindVelGbl, NodeFrcGbl )

         CALL TwrAeroLoads ( p, Node, u%Twr_InputMarkers%Orientation(:,:,Node), u%Twr_InputMarkers%TranslationVel(:,Node) &
                           , y%IfW_Outputs%Velocity(:,Node+p%NumBl*p%Element%NElm), y%Twr_OutputLoads%Force(:,Node) )

      END DO ! Node

   END IF ! ( p%TwrProps%CalcTwrAero )

   !................................................................................................
      
   !************************************************************************************************************
   !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...
   !............................................................................................................
   ! Umass WInDS starts
   !
   ! The WINDS main part:
   !------------------------------------     
      
   IF ( p%FVM%UseWINDS ) THEN  
      O%FVM_Other%ZTime = time  ! The current simulation time (actual or time of prediction)
      !........................................................
      ! The first timestep
         
      IF ( O%Aerodyn_Timestep == -1 )  THEN

             ! Start time
            CALL DATE_AND_TIME ( Values=O%FVM_Other%SimStrtTime )
            CALL CPU_TIME ( O%FVM_Other%UsrTime2 )                                                    ! Initial CPU time   
            O%FVM_Other%UsrTime2 = MAX( 0.0_DbKi, O%FVM_Other%UsrTime2 )  ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
            
               ! Calculate positions and velocity of blades   
            CALL WINDS_Kinematics(u, p, O, xd, O%WINDS_Timestep, ErrStat, ErrMess)
            IF (ErrStat /= 0 ) RETURN
            
            CALL WINDS_Velocity(u, p, O, xd, O%WINDS_Timestep, ErrStat, ErrMess)    
            IF (ErrStat /= 0 ) RETURN
            
            CALL WINDS_FVMInitial(u, p, O, xd, ErrStat, ErrMess, x, z, y)   !  "initials" in WInDS            
            IF (ErrStat /= 0 ) RETURN
            
            DO IBlade=1,p%NumBl
               DO IElement=1,p%Element%Nelm
                  y%OutputLoads(IBlade)%Force(1,IElement)  = O%FVM_Other%StoredForces(1, 1, 1, IElement, IBlade)   
                  y%OutputLoads(IBlade)%Force(2,IElement)  = O%FVM_Other%StoredForces(2, 1, 1, IElement, IBlade)   
                  y%OutputLoads(IBlade)%Force(3,IElement)  = O%FVM_Other%StoredForces(3, 1, 1, IElement, IBlade)   
                  y%OutputLoads(IBlade)%Moment(1,IElement) = O%FVM_Other%StoredMoments(1, 1, 1, IElement, IBlade)
                  y%OutputLoads(IBlade)%Moment(2,IElement) = O%FVM_Other%StoredMoments(2, 1, 1, IElement, IBlade)
                  y%OutputLoads(IBlade)%Moment(3,IElement) = O%FVM_Other%StoredMoments(3, 1, 1, IElement, IBlade)                      
                   
                   
                  O%FVM_Other%PreviousForces(1, 1, 1, IElement, IBlade)  = o%StoredForces(1,IElement,IBlade)  
                  O%FVM_Other%PreviousForces(2, 1, 1, IElement, IBlade)  = o%StoredForces(2,IElement,IBlade) 
                  O%FVM_Other%PreviousForces(3, 1, 1, IElement, IBlade)  = o%StoredForces(3,IElement,IBlade)   
                  O%FVM_Other%PreviousMoments(1, 1, 1, IElement, IBlade) = o%StoredMoments(1,IElement,IBlade)
                  O%FVM_Other%PreviousMoments(2, 1, 1, IElement, IBlade) = o%StoredMoments(2,IElement,IBlade)
                  O%FVM_Other%PreviousMoments(3, 1, 1, IElement, IBlade) = o%StoredMoments(3,IElement,IBlade)                               
               ENDDO
            ENDDO              
                       !................................................. 
                       ! Debug option to ouput blade element data.
                       IF (p%FVM%element_output) THEN                             
                           CALL print_element(P, xd, O, ErrStat, ErrMess)                           
                          ! ! umass debug.......................................
                          ! PRINT_NAME1 = 0.0
                          ! WRITE (temp_number, "(I6.6)") (1)    ! String of the integer With heading zeros
                          ! 
                          ! ! CL
                          ! DO IElement = 1, NS
                          !     DO IBlade =  1,NB
                          !        PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_CL(1, 1, 1, IElement, IBlade) 
                          !     END DO    
                          ! END DO    
                          ! CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_cl_'// TRIM(temp_number))
                          ! 
                          ! ! CD
                          ! PRINT_NAME1 = 0.0
                          ! DO IElement = 1, NS
                          !     DO IBlade =  1,NB
                          !        PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_CD(1, 1, 1, IElement, IBlade) 
                          !     END DO    
                          ! END DO    
                          ! CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_cd_'// TRIM(temp_number))   
                          ! 
                          ! ! AOA
                          ! PRINT_NAME1 = 0.0
                          ! DO IElement = 1, NS
                          !     DO IBlade =  1,NB
                          !        PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_AOA(1, 1, 1, IElement, IBlade) 
                          !     END DO    
                          ! END DO    
                          ! CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_aoa_'// TRIM(temp_number))  
                          ! 
                          ! !V_tot                           
                          ! DO IBlade =  1,NB
                          !    WRITE (temp_1, "(I1)") (IBlade)    ! String of the integer With heading zeros 
                          !    PRINT_NAME3 = 0.0
                          !    DO IElement = 1, NS
                          !        DO IDim =  1,3
                          !           PRINT_NAME3(IElement, IDIM) = O%FVM_Other%KJ%VEL_TOT(IDIM, 1, 1, IElement, IBlade) 
                          !        END DO    
                          !    END DO    
                          !    CALL SAVE_TO_TXT_2D(p, PRINT_NAME3 , 'WINDS_vtot_'//temp_1//'_blade_'// TRIM(temp_number)) 
                          ! END DO                            
                          !! umass debug.......................................
                       END IF  
                       !.................................................

            ! ! End time
            !CALL CPU_TIME(Time_2) 
            !O%FVM_Other%TIME%Time_Total = O%FVM_Other%TIME%Time_Total + Time_2 - Time_1  
              
            O%WINDS_Timestep =O%WINDS_Timestep + 1      ! WInDS internal timestep, which is 1 ,2, 3...

      !........................................................
      ! The global timesteps used by WINDS
           
      ELSE IF (MOD(O%Aerodyn_Timestep , p%FVM%DT_RATIO) == 0 .AND. O%Aerodyn_Timestep /=0) THEN
            
                ! IF ( p%FVM%SteadyFlag ) THEN ! Make it steady flow
               DO IBlade = 1,p%NumBl
                  DO IElement=1,p%Element%Nelm   
                     O%FVM_Other%WIND_INFTY(1, O%WINDS_Timestep, 1, IElement, IBlade) = O%FVM_Other%WIND_INFTY(1, 1, 1, IElement, IBlade)
                     O%FVM_Other%WIND_INFTY(2, O%WINDS_Timestep, 1, IElement, IBlade) = O%FVM_Other%WIND_INFTY(2, 1, 1, IElement, IBlade)
                     O%FVM_Other%WIND_INFTY(3, O%WINDS_Timestep, 1, IElement, IBlade) = O%FVM_Other%WIND_INFTY(3, 1, 1, IElement, IBlade) 
                     
                     O%FVM_Other%WIND_INFTYM(1, O%WINDS_Timestep, 1, IElement, IBlade) = O%FVM_Other%WIND_INFTYM(1,1,1,IElement,IBlade)                
                  END DO                  
               END DO             
               !END IF  ! p%FVM%SteadyFlag            
          
                 ! Check if Wake should be cut off
               CALL WINDS_check_cutoff(u, p, O, xd, O%WINDS_Timestep, ErrStat, ErrMess)               
               IF (ErrStat /= 0 ) RETURN
               
                 ! Calculate positions and velocity of blades   
               CALL WINDS_Kinematics(u, p, O, xd, O%WINDS_Timestep, ErrStat, ErrMess)
               IF (ErrStat /= 0 ) RETURN
               
               CALL WINDS_Velocity(u, p, O, xd, O%WINDS_Timestep, ErrStat, ErrMess)    
               IF (ErrStat /= 0 ) RETURN
               
               CALL WINDS_FVM(u, p, O, xd, O%WINDS_Timestep, ErrStat, ErrMess)   ! Part of the WInDS main driver       
               IF (ErrStat /= 0 ) RETURN
               
               DO IBlade=1,p%NumBl
                  DO IElement=1,p%Element%Nelm
                     y%OutputLoads(IBlade)%Force(1,IElement)  = O%FVM_Other%StoredForces(1, O%WINDS_Timestep, 1, IElement, IBlade)   
                     y%OutputLoads(IBlade)%Force(2,IElement)  = O%FVM_Other%StoredForces(2, O%WINDS_Timestep, 1, IElement, IBlade)   
                     y%OutputLoads(IBlade)%Force(3,IElement)  = O%FVM_Other%StoredForces(3, O%WINDS_Timestep, 1, IElement, IBlade)   
                     y%OutputLoads(IBlade)%Moment(1,IElement) = O%FVM_Other%StoredMoments(1, O%WINDS_Timestep, 1, IElement, IBlade)
                     y%OutputLoads(IBlade)%Moment(2,IElement) = O%FVM_Other%StoredMoments(2, O%WINDS_Timestep, 1, IElement, IBlade)
                     y%OutputLoads(IBlade)%Moment(3,IElement) = O%FVM_Other%StoredMoments(3, O%WINDS_Timestep, 1, IElement, IBlade)                 
                 
                     O%FVM_Other%PreviousForces(1, 1, 1, IElement, IBlade)  = y%OutputLoads(IBlade)%Force(1,IElement)  
                     O%FVM_Other%PreviousForces(2, 1, 1, IElement, IBlade)  = y%OutputLoads(IBlade)%Force(2,IElement) 
                     O%FVM_Other%PreviousForces(3, 1, 1, IElement, IBlade)  = y%OutputLoads(IBlade)%Force(3,IElement)    
                     O%FVM_Other%PreviousMoments(1, 1, 1, IElement, IBlade) = y%OutputLoads(IBlade)%Moment(1,IElement)
                     O%FVM_Other%PreviousMoments(2, 1, 1, IElement, IBlade) = y%OutputLoads(IBlade)%Moment(2,IElement) 
                     O%FVM_Other%PreviousMoments(3, 1, 1, IElement, IBlade) = y%OutputLoads(IBlade)%Moment(3,IElement)                       
                  ENDDO
               ENDDO  
                
                       ! Debug option to ouput blade element data.
                       IF (p%FVM%element_output) THEN
                           CALL print_element(P, xd, O, ErrStat, ErrMess) 
                          !!!..........................................................
                          !!!Umass debug              
                          ! WRITE (temp_number , "(I6.6)") ( O%WINDS_Timestep )   
                          !
                          ! ! CL
                          ! PRINT_NAME1 = 0.0
                          ! DO IElement = 1, NS
                          !     DO IBlade =  1,NB
                          !         PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_CL(1, O%WINDS_Timestep, 1, IElement, IBlade) 
                          !     END DO  
                          ! END DO  
                          ! CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_cl_'//TRIM(temp_number))    ! Write cl to text
                          !
                          ! ! CD
                          ! PRINT_NAME1 = 0.0
                          ! DO IElement = 1, NS
                          !     DO IBlade =  1,NB
                          !         PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_CD(1, O%WINDS_Timestep, 1, IElement, IBlade) 
                          !     END DO  
                          ! END DO  
                          ! CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_cd_'//TRIM(temp_number))    ! Write cd to text    
                          !
                          !
                          ! ! AOA
                          ! PRINT_NAME1 = 0.0
                          ! DO IElement = 1, NS
                          !     DO IBlade =  1,NB
                          !        PRINT_NAME1(IElement, IBlade) = O%FVM_Other%PERF_AOA(1, O%WINDS_Timestep, 1, IElement, IBlade) 
                          !     END DO    
                          ! END DO    
                          ! CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_aoa_'// TRIM(temp_number))  
                          ! 
                          !
                          !  !V_tot                           
                          ! DO IBlade =  1,NB
                          !    WRITE (temp_1, "(I1)") (IBlade)    ! String of the integer With heading zeros 
                          !    PRINT_NAME3 = 0.0
                          !    DO IElement = 1, NS
                          !        DO IDim =  1,3
                          !           PRINT_NAME3(IElement, IDIM) = O%FVM_Other%KJ%VEL_TOT(IDIM, O%WINDS_Timestep, 1, IElement, IBlade) 
                          !        END DO    
                          !    END DO    
                          !    CALL SAVE_TO_TXT_2D(p, PRINT_NAME3 , 'WINDS_vtot_'//temp_1//'_blade_'// TRIM(temp_number)) 
                          ! END DO                            
                          ! !!Umass debug         
                          ! !!..........................................................
                       END IF   
                       
                    
              O%WINDS_Timestep = O%WINDS_Timestep + 1      ! WInDS internal timestep, which is 1 ,2, 3...
              
              ! sliu: Curently, I don't access to get the FAST simulation time. User needs to set the simulation time seperately....
               IF (O%WINDS_Timestep > p%FVM%NT ) THEN
                   ErrStat = ErrID_Fatal
                   ErrMESS = ' The user setting simulation time of WInDS is shorter than the FAST simulation time. '//&
                              ' Please check the setting. '
                   RETURN
               END IF
               
         !........................................................
         ! The global timesteps ignored by WINDS 
                
      ELSE
              ! Copy the load from the previous timestep(refer to O%Aerodyn_Timestep)
            DO IBlade=1,p%NumBl
               DO IElement=1,p%Element%Nelm
                  y%OutputLoads(IBlade)%Force(1,IElement)  = O%FVM_Other%PreviousForces(1, 1, 1, IElement, IBlade)   
                  y%OutputLoads(IBlade)%Force(2,IElement)  = O%FVM_Other%PreviousForces(2, 1, 1, IElement, IBlade)
                  y%OutputLoads(IBlade)%Force(3,IElement)  = O%FVM_Other%PreviousForces(3, 1, 1, IElement, IBlade)  
                  y%OutputLoads(IBlade)%Moment(1,IElement) = O%FVM_Other%PreviousMoments(1, 1, 1, IElement, IBlade)
                  y%OutputLoads(IBlade)%Moment(2,IElement) = O%FVM_Other%PreviousMoments(2, 1, 1, IElement, IBlade)
                  y%OutputLoads(IBlade)%Moment(3,IElement) = O%FVM_Other%PreviousMoments(3, 1, 1, IElement, IBlade)
               ENDDO
            ENDDO         
      END IF ! ( O%Aerodyn_Timestep == -1 ) 
         
      O%Aerodyn_Timestep = O%Aerodyn_Timestep  + 1   ! Update the global timestep


   ENDIF  ! p%FVM%UseWINDS
   
   !------------------------------------
   ! Umass WInDS ends
   !............................................................................................................
   !....WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...WInDS...   
   !************************************************************************************************************

   
   
     
   
   
   CALL ElemOut(time, P, O )
   
   CALL CleanUp (  )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
   SUBROUTINE CleanUp ( )


      ! This subroutine cleans up the parent routine before exiting.


      !   ! Deallocate the IfW_Inputs%Position array if it had been allocated.
      !
      !CALL IfW_DestroyInput( IfW_Inputs, ErrStatLcl, ErrMessLcl )


      RETURN

   END SUBROUTINE CleanUp 
   

END SUBROUTINE AD_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMess )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Compute the first time derivatives of the continuous states here:

!     dxdt%DummyDiscState = 0.


END SUBROUTINE AD_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMess )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                  !   Output: Discrete states at Time + Interval
      TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE AD_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMess )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using
                                                                  !     the input values described above
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Solve for the constraint states here:

!      z_residual%DummyConstrState = 0.

END SUBROUTINE AD_CalcConstrStateResidual



!====================================================================================================
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------------------------------------------------------------------------------------------------------------------------

END MODULE AeroDyn
!**********************************************************************************************************************************

