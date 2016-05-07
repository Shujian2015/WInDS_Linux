MODULE IfW_FFWind
!
!  This module uses full-field binary wind files to determine the wind inflow.
!  This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
!  and that all units are specified in the metric system (using meters and seconds).
!  Data is shifted by half the grid width to account for turbine yaw (so that data in the X
!  direction actually starts at -1*OtherStates%FFYHWid meters).
!
!  Created 25-Sep-2009 by B. Jonkman, National Renewable Energy Laboratory
!     using subroutines and modules from AeroDyn v12.58
!
!----------------------------------------------------------------------------------------------------
!  Feb 2013    v2.00.00          A. Platt
!     -- updated to the new framework
!     -- Modified to use NWTC_Library v. 2.0
!     -- Note:  Jacobians are not included in this version.
!
!----------------------------------------------------------------------------------------------------
! File last committed: $Date: 2015-02-09 12:53:34 -0700 (Mon, 09 Feb 2015) $
! (File) Revision #: $Rev: 141 $
! URL: $HeadURL: https://windsvn.nrel.gov/InflowWind/branches/modularization/Source/IfW_FFWind.f90 $
!----------------------------------------------------------------------------------------------------
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
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
!----------------------------------------------------------------------------------------------------

   USE                                          NWTC_Library
   USE                                          IfW_FFWind_Types

   IMPLICIT                                     NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_FFWind_Ver = ProgDesc( 'IfW_FFWind', 'v2.00.00', '17-Sep-2013' )

   PUBLIC                                    :: IfW_FFWind_Init
   PUBLIC                                    :: IfW_FFWind_End
   PUBLIC                                    :: IfW_FFWind_CalcOutput


      ! The following do not contain anything since there are no states.
   PUBLIC                                    :: IfW_FFWind_UpdateStates
   PUBLIC                                    :: IfW_FFWind_CalcContStateDeriv
   PUBLIC                                    :: IfW_FFWind_UpdateDiscState
   PUBLIC                                    :: IfW_FFWind_CalcConstrStateResidual



CONTAINS
!====================================================================================================
SUBROUTINE IfW_FFWind_Init(InitData,   InputGuess, ParamData,                       &
                           ContStates, DiscStates, ConstrStates,     OtherStates,   &
                           OutData,    Interval,   InitOutData,      ErrStat,       ErrMsg)
   !-------------------------------------------------------------------------------------------------
   !  This routine is used read the full-field turbulence data.
   !  09/25/1997  - Created by M. Buhl from GETFILES in ViewWind.
   !  09/23/2009  - modified by B. Jonkman: this subroutine was split into several subroutines (was ReadFF)
   !  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !-------------------------------------------------------------------------------------------------

   IMPLICIT                       NONE


      ! Passed Variables
   TYPE(IfW_FFWind_InitInputType),        INTENT(IN   )  :: InitData       ! Initialization data passed to the module
   TYPE(IfW_FFWind_InputType),            INTENT(  OUT)  :: InputGuess     ! Initialized input data variable
   TYPE(IfW_FFWind_ParameterType),        INTENT(  OUT)  :: ParamData      ! Parameters
   TYPE(IfW_FFWind_ContinuousStateType),  INTENT(  OUT)  :: ContStates     ! Continuous States  (unused)
   TYPE(IfW_FFWind_DiscreteStateType),    INTENT(  OUT)  :: DiscStates     ! Discrete States    (unused)
   TYPE(IfW_FFWind_ConstraintStateType),  INTENT(  OUT)  :: ConstrStates   ! Constraint States  (unused)
   TYPE(IfW_FFWind_OtherStateType),       INTENT(  OUT)  :: OtherStates    ! Other State data   (storage for the main data)
   TYPE(IfW_FFWind_OutputType),           INTENT(  OUT)  :: OutData        ! Initial output
   TYPE(IfW_FFWind_InitOutputType),       INTENT(  OUT)  :: InitOutData    ! Initial output

   REAL(DbKi),                            INTENT(IN   )  :: Interval       ! Time Interval to use (passed through here)


      ! Error Handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat        ! determines if an error has been encountered
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg         ! Message about errors


      ! Temporary variables for error handling
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg      ! temporary error message


      ! Local Variables:

   REAL(ReKi)                                            :: TI      (3)    ! turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   REAL(ReKi)                                            :: BinTI   (3)    ! turbulence intensities of the wind components as defined in the FF binary file, not necessarially the actual TI
   REAL(ReKi)                                            :: UBar
   REAL(ReKi)                                            :: ZCenter

   INTEGER(B2Ki)                                         :: Dum_Int2
   INTEGER(IntKi)                                        :: DumInt
   INTEGER(IntKi)                                        :: I
   LOGICAL                                               :: CWise
   LOGICAL                                               :: Exists
   CHARACTER( 1028 )                                     :: SumFile        ! length is LEN(ParamData%WindFileName) + the 4-character extension.
   CHARACTER( 1028 )                                     :: TwrFile        ! length is LEN(ParamData%WindFileName) + the 4-character extension.



      !-------------------------------------------------------------------------------------------------
      ! Initialize temporary variables
      !-------------------------------------------------------------------------------------------------

   ErrMsg      = ''
   ErrStat     = ErrID_None

   TmpErrMsg   = ''
   TmpErrStat  = ErrID_None


      !-------------------------------------------------------------------------------------------------
      ! Set values for unused output types.
      !     This is just to keep the compiler from complaining (not really necessary).
      !-------------------------------------------------------------------------------------------------

      ! Allocate the empty position array.
   CALL AllocAry( InputGuess%Position, 3, 1, &
                  'Empty position array in initialization.', TmpErrStat, TmpErrMsg )
   InputGuess%Position(:,1)      = 0.0
   ErrStat  = MAX(TmpErrStat, ErrStat)
   ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Allocate the empty velocity array.
   CALL AllocAry( OutData%Velocity, 3, 1, &
                  'Empty velocity array in initialization.', TmpErrStat, TmpErrMsg )
   OutData%Velocity(:,1)         = 0.0
   ErrStat  = MAX(TmpErrStat, ErrStat)
   ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   IF ( ErrStat >= AbortErrLev ) RETURN

   ContStates%DummyContState     = 0.0
   DiscStates%DummyDiscState     = 0.0
   ConstrStates%DummyConstrState = 0.0


      !-------------------------------------------------------------------------------------------------
      ! Check that it's not already initialized
      !-------------------------------------------------------------------------------------------------

   IF ( ParamData%Initialized ) THEN
      ErrMsg   = ' FFWind has already been initialized.'
      ErrStat  = ErrId_Warn      ! no reason to stop the program over this
      RETURN
   ELSE
      ErrStat = ErrId_None
   ENDIF


      ! Get a unit number to use

   CALL GetNewUnit(OtherStates%UnitWind, TmpErrStat, TmpErrMsg)
   IF (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(TmpErrStat, ErrStat)      
      ErrMsg   = TRIM(ErrMsg)//"IfW_FFWind: "//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData
      !-------------------------------------------------------------------------------------------------

   ParamData%ReferenceHeight  =  InitData%ReferenceHeight       ! Height of the turbine
   ParamData%Width            =  InitData%Width                 ! width of the wind profile box
   ParamData%WindFileName     =  InitData%WindFileName          ! Filename of the FF wind file


      !----------------------------------------------------------------------------------------------
      ! Open the binary file, read its "header" (first 2-byte integer) to determine what format
      ! binary file it is, and close it.
      !----------------------------------------------------------------------------------------------

   CALL OpenBInpFile (OtherStates%UnitWind, TRIM(ParamData%WindFileName), TmpErrStat, TmpErrMsg)
   ErrStat  = MAX(TmpErrStat, ErrStat)
   ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Read the first binary integer from the file to get info on the type.
      ! Cannot use library read routines since this is a 2-byte integer.
   READ ( OtherStates%UnitWind, IOSTAT=TmpErrStat )  Dum_Int2
   CLOSE( OtherStates%UnitWind )

   IF (TmpErrStat /= 0) THEN
      ErrMsg   = ' Error reading first binary integer from file "'//TRIM(ParamData%WindFileName)//'."'
      ErrStat  = ErrID_Fatal
      RETURN
   ENDIF


      !----------------------------------------------------------------------------------------------
      ! Read the files to get the required FF data.
      !----------------------------------------------------------------------------------------------
   DumInt = Dum_Int2  ! change to default INTEGER, instead of INT(2) to compare in SELECT below

   SELECT CASE (DumInt)

      CASE ( 7, 8 )                                                    ! TurbSim binary format

         CALL Read_TurbSim_FF(OtherStates%UnitWind, TmpErrStat, TmpErrMsg)
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
            IF ( ErrStat >= AbortErrLev ) THEN
               CLOSE ( OtherStates%UnitWind )
               RETURN
            END IF
         
      CASE ( -1, -2, -3, -99 )                                         ! Bladed-style binary format

         !...........................................................................................
         ! Create full-field summary file name from binary file root name.  Also get tower file
         ! name.
         !...........................................................................................

            CALL GetRoot(ParamData%WindFileName, SumFile)

            TwrFile = TRIM(SumFile)//'.twr'
            SumFile = TRIM(SumFile)//'.sum'


         !...........................................................................................
         ! Read the summary file to get necessary scaling information
         !...........................................................................................

            CALL Read_Summary_FF (OtherStates%UnitWind, TRIM(SumFile), CWise, ZCenter, TI, TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( OtherStates%UnitWind )
                  RETURN
               END IF
            
            UBar = OtherStates%MeanFFWS      ! temporary storage .... this is our only check to see if the summary and binary files "match"


         !...........................................................................................
         ! Open the binary file and read its header
         !...........................................................................................

            CALL OpenBInpFile (OtherStates%UnitWind, TRIM(ParamData%WindFileName), TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( OtherStates%UnitWind )
                  RETURN
               END IF
            
            IF ( Dum_Int2 == -99 ) THEN                                                      ! Newer-style BLADED format
               CALL Read_Bladed_FF_Header1 (OtherStates%UnitWind, BinTI, TmpErrStat, TmpErrMsg)
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( OtherStates%UnitWind )
                     RETURN
                  END IF

                  ! If the TIs are also in the binary file (BinTI > 0),
                  ! use those numbers instead of ones from the summary file

               DO I =1,OtherStates%NFFComp
                  IF ( BinTI(I) > 0 ) TI(I) = BinTI(I)
               ENDDO

            ELSE
               CALL Read_Bladed_FF_Header0 (OtherStates%UnitWind, TmpErrStat, TmpErrMsg)     ! Older-style BLADED format
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( OtherStates%UnitWind )
                     RETURN
                  END IF

            ENDIF



         !...........................................................................................
         ! Let's see if the summary and binary FF wind files go together before continuing.
         !...........................................................................................

            IF ( ABS( UBar - OtherStates%MeanFFWS ) > 0.1 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error: Incompatible mean hub-height wind speeds in FF wind files. '//&
                           '(Check that the .sum and .wnd files were generated together.)', ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
               RETURN
            ENDIF


         !...........................................................................................
         ! Calculate the height of the bottom of the grid
         !...........................................................................................

            OtherStates%GridBase = ZCenter - OtherStates%FFZHWid         ! the location, in meters, of the bottom of the grid


         !...........................................................................................
         ! Read the binary grids (converted to m/s) and close the file
         !...........................................................................................

            CALL Read_Bladed_Grids( OtherStates, CWise, TI, TmpErrStat, TmpErrMsg)
            CLOSE ( OtherStates%UnitWind )

               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
               IF ( ErrStat >= AbortErrLev ) RETURN
            
         !...........................................................................................
         ! Read the tower points file
         !...........................................................................................

            INQUIRE ( FILE=TRIM(TwrFile) , EXIST=Exists )

            IF (  Exists )  THEN
               CALL Read_FF_Tower( OtherStates, TmpErrStat, TmpErrMsg )
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( OtherStates%UnitWind )
                     RETURN
                  END IF
            ELSE
               OtherStates%NTGrids = 0
            ENDIF


      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, ' Error: Unrecognized binary wind file type.', ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
         RETURN

   END SELECT


   IF (ParamData%Periodic) THEN
      OtherStates%InitXPosition = 0                ! start at the hub
      OtherStates%TotalTime     = OtherStates%NFFSteps*OtherStates%FFDTime
   ELSE
      OtherStates%InitXPosition = OtherStates%FFYHWid          ! start half the grid with ahead of the turbine
      OtherStates%TotalTime     = (OtherStates%NFFSteps-1)*OtherStates%FFDTime
   ENDIF

   ParamData%Initialized = .TRUE.


      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information
      !-------------------------------------------------------------------------------------------------

   InitOutData%HubHeight   = ParamData%ReferenceHeight
   InitOutdata%Ver         = IfW_FFWind_Ver


      ! Allocate and populate the OutputHdr array (contains names of outputable values)

   CALL AllocAry( InitOutData%WriteOutputHdr, 3, 'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
      IF ( ErrStat >= AbortErrLev ) RETURN
   
   InitOutData%WriteOutputHdr(1) = 'WindVxi'
   InitOutData%WriteOutputHdr(2) = 'WindVyi'
   InitOutData%WriteOutputHdr(3) = 'WindVzi'


      ! Allocate and populate the OutputUnt array (contains units of outputable values)

   CALL AllocAry( InitOutData%WriteOutputUnt, 3, 'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_FFWind_Init' )                  
      IF ( ErrStat >= AbortErrLev ) RETURN

   InitOutData%WriteOutputUnt(1) = '(m/s)'
   InitOutData%WriteOutputUnt(2) = '(m/s)'
   InitOutData%WriteOutputUnt(3) = '(m/s)'


   RETURN


   CONTAINS

   !====================================================================================================
   SUBROUTINE Read_Summary_FF ( UnitWind, FileName, CWise, ZCenter, TI, ErrStat, ErrMsg )
   ! This subroutine reads the text summary file to get normalizing parameters, the location of the
   ! grid, and the direction the grid was written to the binary file
   !
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------

         ! Passed variables
      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind       ! unit number for the file to open
      CHARACTER(*),                       INTENT(IN   )  :: FileName       ! name of the summary file
      LOGICAL,                            INTENT(  OUT)  :: CWise          ! rotation (for reading the order of the binary data)
      REAL(ReKi),                         INTENT(  OUT)  :: ZCenter        ! the height at the center of the grid
      REAL(ReKi),                         INTENT(  OUT)  :: TI      (3)    ! turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        ! returns 0 if no error encountered in the subroutine
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         ! holds the error messages

        ! Local variables
      REAL(ReKi)                                         :: ZGOffset       ! The vertical offset of the turbine on rectangular grid (allows turbulence not centered on turbine hub)

      INTEGER, PARAMETER                                 :: NumStrings = 6 ! number of strings to be looking for in the file

      INTEGER(IntKi)                                     :: FirstIndx      ! The first character of a line where data is located
      INTEGER(IntKi)                                     :: I              ! A loop counter
      INTEGER(IntKi)                                     :: LastIndx       ! The last  character of a line where data is located
      INTEGER(IntKi)                                     :: LineCount      ! Number of lines that have been read in the file

      LOGICAL                                            :: StrNeeded(NumStrings)   ! if the string has been found

      CHARACTER(1024)                                    :: LINE           ! temporary storage for reading a line from the file

         ! Temporary variables for error handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! temporary error status
      CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg      ! temporary error message

         !----------------------------------------------------------------------------------------------
         ! Initialize some variables
         !----------------------------------------------------------------------------------------------

      ErrStat              = ErrID_None
      ErrMsg               = ''

      LineCount            = 0
      StrNeeded(:)         = .TRUE.
      ZGOffset             = 0.0
      OtherStates%RefHt                = 0.0
      ParamData%Periodic   = .FALSE.

         !----------------------------------------------------------------------------------------------
         ! Open summary file.
         !----------------------------------------------------------------------------------------------

      CALL OpenFInpFile ( OtherStates%UnitWind, TRIM( FileName ), TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Read_Summary_FF' )                  
      IF ( ErrStat >= AbortErrLev ) RETURN


         !----------------------------------------------------------------------------------------------
         ! Read the summary file.
         !----------------------------------------------------------------------------------------------

      ! Here are the strings we're looking for, in this order:
      ! 1) 'CLOCKWISE'
      ! 2) 'HUB HEIGHT'
      ! 3)     (unused; decided we didn't need to read data also stored in the binary file)
      ! 4) 'UBAR'
      ! 5) 'HEIGHT OFFSET' (optional)
      ! 6) 'PERIODIC' (optional)
         
         
      DO WHILE ( ( ErrStat == ErrID_None ) .AND. StrNeeded(NumStrings) )

         LineCount = LineCount + 1

         READ ( OtherStates%UnitWind, '(A)', IOSTAT=TmpErrStat ) LINE
         IF ( TmpErrStat /= 0 ) THEN

            IF ( StrNeeded(1) .OR. StrNeeded(2) .OR. StrNeeded(4)  ) THEN  ! the "HEIGHT OFFSET" and "PERIODIC" parameters are not necessary.  We'll assume they are zero/false if we didn't find it.
               CALL SetErrStat( ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'// &
                           TRIM(FileName)//'". Could not find all of the required parameters.', ErrStat, ErrMsg, 'Read_Summary_FF' )                  
               RETURN
            ELSE
               EXIT
            ENDIF

         ENDIF

         CALL Conv2UC ( LINE )


         IF ( StrNeeded(1) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #1: Get the rotation direction, using the string "CLOCKWISE"
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'CLOCKWISE' ) > 0 ) THEN

               READ (LINE, *, IOSTAT = TmpErrStat)  CWise          ! Look for True/False values

               IF ( TmpErrStat /= 0 ) THEN                         ! Look for Yes/No values instead

                  LINE = ADJUSTL ( LINE )                      ! Remove leading spaces from input line

                  SELECT CASE (LINE(1:1) )
                     CASE ('Y')
                        CWise = .TRUE.
                     CASE ('N')
                        CWise = .FALSE.
                     CASE DEFAULT
                        CALL SetErrStat( ErrID_Fatal, ' Error reading rotation direction (CLOCKWISE) from FF summary file.', ErrStat, ErrMsg, 'Read_Summary_FF' )                  
                        RETURN
                  END SELECT

               ENDIF ! TmpErrStat /= 0
               StrNeeded(1) = .FALSE.

            ENDIF   ! INDEX for "CLOCKWISE"

         ELSEIF ( StrNeeded(2) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #2: Get the hub height, using the strings "HUB HEIGHT" or "ZHUB"
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'HUB HEIGHT' ) > 0 .OR. INDEX( LINE, 'ZHUB' ) > 0 ) THEN

               READ (LINE, *, IOSTAT = TmpErrStat) OtherStates%RefHt

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading hub height from FF summary file.', ErrStat, ErrMsg, 'Read_Summary_FF' )                  
                  RETURN
               ENDIF
               StrNeeded(2) = .FALSE.

            ENDIF !INDEX for "HUB HEIGHT" or "ZHUB"

   !      ELSEIF ( StrNeeded(3) ) THEN
   !
   !         !-------------------------------------------------------------------------------------------
   !         ! #3: Get the grid width (& height, if available), using the strings "GRID WIDTH" or "RDIAM"
   !         !    If GRID HEIGHT is specified, use it, too. -- THIS IS UNNECESSARY AS IT'S STORED IN THE BINARY FILE
   !         !-------------------------------------------------------------------------------------------

         ELSEIF ( StrNeeded(4) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #4: Get the mean wind speed "UBAR" and turbulence intensities from following lines for
            !     scaling Bladed-style FF binary files
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'UBAR') > 0 ) THEN

               FirstIndx = INDEX( LINE, '=' ) + 1        ! Look for the equal siqn to find the number we're looking for

               READ ( LINE( FirstIndx:LEN(LINE) ), *, IOSTAT=TmpErrStat ) OtherStates%MeanFFWS

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading UBar binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, 'Read_Summary_FF' )                  
                  RETURN
               ENDIF

               DO I = 1,3

                  LineCount = LineCount + 1

                  READ ( OtherStates%UnitWind, '(A)', IOSTAT=TmpErrStat ) LINE
                  IF ( TmpErrStat /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'//TRIM(FileName)//&
                                          '". Could not find all of the required parameters.', ErrStat, ErrMsg, 'Read_Summary_FF' )                  
                     RETURN
                  ENDIF

                  FirstIndx = INDEX( LINE, '=' ) + 1     ! Read the number between the = and % signs
                  LastIndx  = INDEX( LINE, '%' ) - 1

                  IF ( LastIndx <= FirstIndx ) LastIndx = LEN( LINE )   ! If there's no % sign, read to the end of the line

                  READ ( LINE( FirstIndx:LastIndx ), *, IOSTAT=TmpErrStat ) TI(I)
                  IF ( TmpErrStat /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(I))// &
                                 ') binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, 'Read_Summary_FF' )                  
                     RETURN
                  ENDIF

               ENDDO !I

               StrNeeded(4) = .FALSE.

             ENDIF

         ELSEIF ( StrNeeded(5) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #5: Get the grid "HEIGHT OFFSET", if it exists (in TurbSim). Otherwise, assume it's zero
            !           ZGOffset = HH - GridBase - OtherStates%FFZHWid
            !-------------------------------------------------------------------------------------------
            IF ( INDEX( LINE, 'HEIGHT OFFSET' ) > 0  ) THEN

               FirstIndx = INDEX ( LINE, '=' ) + 1

               READ ( LINE( FirstIndx:LEN(LINE) ), *, IOSTAT=TmpErrStat ) ZGOffset

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading height offset from FF summary file.', ErrStat, ErrMsg, 'Read_Summary_FF' )                  
                  RETURN
               ENDIF

               StrNeeded(5) = .FALSE.

            ENDIF !INDEX for "HEIGHT OFFSET"

         ELSEIF ( StrNeeded(6) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #5: Get the grid "PERIODIC", if it exists (in TurbSim). Otherwise, assume it's
            !        not a periodic file
            !-------------------------------------------------------------------------------------------
            IF ( INDEX( LINE, 'PERIODIC' ) > 0  ) THEN

               ParamData%Periodic   = .TRUE.
               StrNeeded(6)         = .FALSE.

            ENDIF !INDEX for "PERIODIC"

         ENDIF ! StrNeeded


      ENDDO !WHILE


      !-------------------------------------------------------------------------------------------------
      ! Close the summary file
      !-------------------------------------------------------------------------------------------------

      CLOSE ( OtherStates%UnitWind )


      !-------------------------------------------------------------------------------------------------
      ! Calculate the height of the grid center
      !-------------------------------------------------------------------------------------------------

       ZCenter  = OtherStates%RefHt - ZGOffset


   END SUBROUTINE Read_Summary_FF
   !====================================================================================================
   SUBROUTINE Read_TurbSim_FF(UnitWind, ErrStat, ErrMsg)
   ! This subroutine reads the binary TurbSim-format FF file (.bts).  It fills the FFData array with
   ! velocity data for the grids and fills the FFTower array with velocities at points on the tower
   ! (if data exists).
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------

         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind       ! unit number for the wind file
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        ! error status return value (0=no error; non-zero is error)
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         ! message about the error encountered

         ! Local Variables:

      REAL(SiKi)                                         :: Dum_Real4      ! dummy 4-byte real number
      INTEGER(B1Ki)                                      :: Dum_Int1       ! dummy 1-byte integer
      INTEGER(B2Ki)                                      :: Dum_Int2       ! dummy 2-byte integer
      INTEGER(B4Ki)                                      :: Dum_Int4       ! dummy 4-byte integer

      INTEGER(IntKi)                                     :: IC             ! loop counter for wind components
      INTEGER(IntKi)                                     :: IT             ! loop counter for time
      INTEGER(IntKi)                                     :: IY             ! loop counter for y
      INTEGER(IntKi)                                     :: IZ             ! loop counter for z
      INTEGER(IntKi)                                     :: NChar          ! number of characters in the description string

      REAL(SiKi)                                         :: Vslope(3)      ! slope  for "un-normalizing" data
      REAL(SiKi)                                         :: Voffset(3)     ! offset for "un-normalizing" data

      LOGICAL                                            :: FirstWarn      ! we don't need to print warning for each character that exceeds the 
      CHARACTER(1024)                                    :: DescStr        ! description string contained in the file

      
         ! Temporary variables for error handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! temporary error status
      CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg      ! temporary error message


      OtherStates%NFFComp = 3                                              ! this file contains 3 wind components
      ErrStat = ErrID_None
      ErrMsg  = ""
      
   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------

      CALL OpenBInpFile (OtherStates%UnitWind, TRIM(ParamData%WindFileName), TmpErrStat, TmpErrMsg)
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
      IF ( ErrStat >= AbortErrLev ) RETURN

      !-------------------------------------------------------------------------------------------------
      ! Read the header information
      !-------------------------------------------------------------------------------------------------
            ! Read in the 2-byte integer. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int2             ! the file identifier, INT(2)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the file identifier in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            ParamData%Periodic = Dum_Int2 == INT( 8, B2Ki) ! the number 7 is used for non-periodic wind files; 8 is periodic wind


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of grid points vertically, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of z grid points in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%NZGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of grid points laterally, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of y grid points in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%NYGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of tower points, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of tower points in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%NTGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of time steps, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of time steps in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%NFFSteps = Dum_Int4


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in vertical direction (dz), REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dz in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%InvFFZD = 1.0/Dum_Real4                            ! 1/dz
            OtherStates%FFZHWid = 0.5*(OtherStates%NZGrids-1)*Dum_Real4                ! half the grid height


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in lateral direction (dy), REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dy in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%InvFFYD = 1.0 / Dum_Real4                          ! 1/dy
            OtherStates%FFYHWid = 0.5*(OtherStates%NYGrids-1)*Dum_Real4                ! half grid grid width


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in time (dt), REAL(4), in m/s
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dt in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%FFDTime = Dum_Real4
            OtherStates%FFRate  = 1.0/OtherStates%FFDTime


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! the mean wind speed at hub height, REAL(4), in m/s
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading mean wind speed in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%MeanFFWS = Dum_Real4
            OtherStates%InvMFFWS = 1.0 / OtherStates%MeanFFWS


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! height of the hub, REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading zHub in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%RefHt = Dum_Real4


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! height of the bottom of the grid, REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading GridBase in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
               RETURN
            ENDIF
            OtherStates%GridBase = Dum_Real4

    !        ZGOffset = OtherStates%RefHt - OtherStates%GridBase  - OtherStates%FFZHWid


         !----------------------------------------------------------------------------------------------
         ! Read the binary scaling factors
         !----------------------------------------------------------------------------------------------

            DO IC = 1,OtherStates%NFFComp
                  ! Read in the 4-byte real. Can't use library read routines for this.
               READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Vslope(IC)     ! the IC-component slope for scaling, REAL(4)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading Vslope('//Num2LStr(IC)//') in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
                     RETURN
                  ENDIF


                  ! Read in the 4-byte real. Can't use library read routines for this.
               READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Voffset(IC)    ! the IC-component offset for scaling, REAL(4)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading Voffset('//Num2LStr(IC)//') in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
                     RETURN
                  ENDIF

            ENDDO !IC


         !----------------------------------------------------------------------------------------------
         ! Read the description string: "Generated by TurbSim (vx.xx, dd-mmm-yyyy) on dd-mmm-yyyy at hh:mm:ss."
         !----------------------------------------------------------------------------------------------

            ! Read in the 4-byte integer. Can't use library read routines for this.
            READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int4          ! the number of characters in the description string, max 200, INT(4)
               IF ( TmpErrStat /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading NCHAR in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' )                  
                  RETURN
               ENDIF
               nchar = Dum_Int4

            DescStr = ''                                       ! Initialize the description string
            FirstWarn = .true.
            
            DO IC=1,nchar

                  ! Read in the 1-byte integer. Can't use library read routines for this.
               READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int1       ! the ASCII integer representation of the character, INT(1)
               IF ( TmpErrStat /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading description line in the FF binary file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' ) 
                  RETURN
               ENDIF

               IF ( LEN(DescStr) >= IC ) THEN
                  DescStr(IC:IC) = ACHAR( Dum_Int1 )              ! converted ASCII characters
               ELSEIF ( FirstWarn ) THEN
                  FirstWarn = .FALSE.
                  CALL SetErrStat( ErrID_Info, ' Description string was too long for variable.'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' ) 
               ENDIF

            ENDDO !IC


      !-------------------------------------------------------------------------------------------------
      ! Get the grid and tower velocities
      !-------------------------------------------------------------------------------------------------

         ! this could take a while, so we'll write a message indicating what's going on:
         
         CALL WrScr( NewLine//'   Reading a '//TRIM( Num2LStr(OtherStates%NYGrids) )//'x'//TRIM( Num2LStr(OtherStates%NZGrids) )//  &
                    ' grid ('//TRIM( Num2LStr(OtherStates%FFYHWid*2) )//' m wide, '// &
                    TRIM( Num2LStr(OtherStates%GridBase) )//' m to '// &
                    TRIM( Num2LStr(OtherStates%GridBase+OtherStates%FFZHWid*2) )//&
                    ' m above ground) with a characteristic wind speed of '// &
                    TRIM( Num2LStr(OtherStates%MeanFFWS) )//' m/s. '//TRIM(DescStr) )


      !----------------------------------------------------------------------------------------------
      ! Allocate arrays for the FF grid as well as the tower points, if they exist
      !----------------------------------------------------------------------------------------------

         IF ( .NOT. ALLOCATED( OtherStates%FFData ) ) THEN
            CALL AllocAry( OtherStates%FFData, OtherStates%NZGrids, OtherStates%NYGrids, OtherStates%NFFComp, OtherStates%NFFSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'READ_TurbSim_FF' ) 
            IF ( ErrStat >= AbortErrLev ) RETURN
         ENDIF


         IF ( OtherStates%NTGrids > 0 ) THEN

            IF ( .NOT. ALLOCATED( OtherStates%FFTower ) ) THEN
               CALL AllocAry( OtherStates%FFTower, OtherStates%NFFComp, OtherStates%NTGrids, OtherStates%NFFSteps, &
                     'Tower wind file data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'READ_TurbSim_FF' ) 
               IF ( ErrStat >= AbortErrLev ) RETURN
            ENDIF

         ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Read the 16-bit data and scale it to 32-bit reals
      !-------------------------------------------------------------------------------------------------

         ! Loop through time.

         DO IT=1,OtherStates%NFFSteps

            !...........................................................................................
            ! Read grid data at this time step.
            !...........................................................................................

            DO IZ=1,OtherStates%NZGrids
               ! Zgrid(IZ) = Z1 + (IZ-1)*dz                 ! Vertical location of grid data point, in m relative to ground

               DO IY=1,OtherStates%NYGrids
                  ! Ygrid(IY) = -0.5*(ny-1)*dy + (IY-1)*dy  ! Horizontal location of grid data point, in m relative to tower centerline

                  DO IC=1,OtherStates%NFFComp                           ! number of wind components (U, V, W)

                        ! Read in the 2-byte integer. Can't use library read routines for this.
                     READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                     IF ( TmpErrStat /= 0 )  THEN                        
                        CALL SetErrStat( ErrID_Fatal, ' Error reading grid wind components in the FF binary file "'// &
                                    TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'READ_TurbSim_FF' ) 
                        RETURN
                     ENDIF

                     OtherStates%FFData(IZ,IY,IC,IT) = ( Dum_Int2 - Voffset(IC) ) / VSlope(IC)

                  ENDDO !IC

               ENDDO !IY

            ENDDO ! IZ


            !...........................................................................................
            ! Read the tower data at this time step.
            !...........................................................................................

            DO IZ=1,OtherStates%NTGrids         ! If NTGrids<1, there are no tower points & FFTower is not allocated

               ! Ytower     = 0               ! Lateral location of the tower data point, in m relative to tower centerline
               ! Ztower(IZ) = Z1 - (IZ-1)*dz  ! Vertical location of tower data point, in m relative to ground

               DO IC=1,OtherStates%NFFComp   ! number of wind components

                     ! Read in a 2-byte integer. Can't use library routines for this.
                  READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading tower wind components in the FF binary file "'//TRIM(ParamData%WindFileName)//'."'&
                                      , ErrStat, ErrMsg, 'READ_TurbSim_FF' )                    
                     RETURN
                  ENDIF

                  OtherStates%FFTower(IC,IZ,IT) = ( Dum_Int2 - Voffset(IC) ) / VSlope(IC)  ! wind-component scaled to m/s

               ENDDO !IC

            ENDDO ! IZ


         ENDDO ! IT

      !-------------------------------------------------------------------------------------------------
      ! close the file and return
      !-------------------------------------------------------------------------------------------------

      CLOSE ( OtherStates%UnitWind )

      IF ( ParamData%Periodic ) THEN
         TmpErrMsg   = '   Processed '//TRIM( Num2LStr( OtherStates%NFFSteps ) )//' time steps of '// &
                        TRIM( Num2LStr ( OtherStates%FFRate ) )//'-Hz full-field data (period of '// &
                        TRIM( Num2LStr( OtherStates%FFDTime*( OtherStates%NFFSteps ) ) )//' seconds).'
      ELSE
         TmpErrMsg   = '   Processed '//TRIM( Num2LStr( OtherStates%NFFSteps ) )//' time steps of '// &
                        TRIM( Num2LStr ( OtherStates%FFRate ) )//'-Hz full-field data ('// &
                        TRIM( Num2LStr( OtherStates%FFDTime*( OtherStates%NFFSteps - 1 ) ) )//' seconds).'
      ENDIF
  
      CALL WrScr( NewLine//TRIM(TmpErrMsg) ) 
      
!     CALL SetErrStat( ErrID_Info, TmpErrMsg, ErrStat, ErrMsg, 'READ_TurbSim_FF' )                    
      
      RETURN

   END SUBROUTINE READ_TurbSim_FF

   !====================================================================================================
   SUBROUTINE Read_Bladed_FF_Header0 (UnitWind, ErrStat, ErrMsg)
   !   Reads the binary headers from the turbulence files of the old Bladed variety.  Note that
   !   because of the normalization, neither OtherStates%NZGrids or OtherStates%NYGrids are larger than 32 points.
   !   21-Sep-2009 - B. Jonkman, NREL/NWTC.
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------


      IMPLICIT                                                 NONE


         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


         ! Local Variables:
      REAL(ReKi)                                         :: FFXDelt
      REAL(ReKi)                                         :: FFYDelt
      REAL(ReKi)                                         :: FFZDelt

      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(IntKi)                                     :: I


         ! Temporary Error Handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! for checking the IOSTAT from a READ or Open statement
      CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg      ! Temporary ErrMsg


      !-------------------------------------------------------------------------------------------------
      ! Initializations
      !-------------------------------------------------------------------------------------------------

      ErrStat  = ErrID_None
      ErrMsg   = ''


      !-------------------------------------------------------------------------------------------------
      ! Read the header (file has just been opened)
      !-------------------------------------------------------------------------------------------------

         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! -NFFC (file ID)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of wind components from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         OtherStates%NFFComp = -1*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta z (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dz from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         FFZDelt = 0.001*Dum_Int2
         OtherStates%InvFFZD = 1.0/FFZDelt


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta y (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         FFYDelt = 0.001*Dum_Int2
         OtherStates%InvFFYD = 1.0/FFYDelt


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta x (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dx from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         FFXDelt = 0.001*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! half the number of time steps

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of time steps from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         OtherStates%NFFSteps = 2*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 10 times the mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         OtherStates%MeanFFWS = 0.1*Dum_Int2
         OtherStates%InvMFFWS = 1.0/OtherStates%MeanFFWS
         OtherStates%FFDTime  = FFXDelt/OtherStates%MeanFFWS
         OtherStates%FFRate   = 1.0/OtherStates%FFDTime


      DO I = 1,5

         ! Read 2-byte integer. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                              ! unused variables: zLu, yLu, xLu, dummy, random seed

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 2-byte integers from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
               RETURN
            ENDIF

      END DO


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 1000*nz

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading nz from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         OtherStates%NZGrids  = Dum_Int2/1000
         OtherStates%FFZHWid  = 0.5*FFZDelt*( OtherStates%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 1000*ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
            RETURN
         ENDIF
         OtherStates%NYGrids  = Dum_Int2/1000
         OtherStates%FFYHWid  = 0.5*FFYDelt*( OtherStates%NYGrids - 1 )


      IF (OtherStates%NFFComp == 3) THEN

         DO I=1,6

               ! Read 2-byte integer. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                           ! unused variables: zLv, yLv, xLv, zLw, yLw, xLw

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 2-byte length scales from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header0')
                  RETURN
               ENDIF

         ENDDO !I

      ENDIF !NFFComp


      RETURN

   END SUBROUTINE Read_Bladed_FF_Header0
   !====================================================================================================
   SUBROUTINE Read_Bladed_FF_Header1 (UnitWind, TI, ErrStat, ErrMsg)
   !   Reads the binary headers from the turbulence files of the new Bladed variety.
   !   16-May-2002 - Windward Engineering.
   !   21-Sep-2009 - B. Jonkman, NREL.  updated to trap errors and add extra parameters for MANN model
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------


      IMPLICIT                                              NONE


         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind
      REAL(ReKi),                         INTENT(  OUT)  :: TI(3)
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


         ! Local Variables:

      REAL(ReKi)                                         :: FFXDelt
      REAL(ReKi)                                         :: FFYDelt
      REAL(ReKi)                                         :: FFZDelt

      REAL(SiKi)                                         :: Dum_Real4
      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(B4Ki)                                      :: Dum_Int4

      INTEGER(IntKi)                                     :: I
      INTEGER(IntKi)                                     :: TurbType


         ! Temporary Error Handling
      INTEGER(IntKi)                                     :: TmpErrStat
      CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg


      !-------------------------------------------------------------------------------------------------
      ! Initializations
      !-------------------------------------------------------------------------------------------------

      ErrStat  = ErrID_None
      ErrMsg   = ''

      TI(:) = -1                                                                                !Initialize to -1 (not all models contain TI)

      !-------------------------------------------------------------------------------------------------
      ! File reading
      !-------------------------------------------------------------------------------------------------

         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! -99 (file ID)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading integer from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! turbulence type

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading turbulence type from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         TurbType = Dum_Int2


      SELECT CASE (TurbType)
         CASE(1, 2)
            !----------------------------------------
            !1-component Von Karman (1) or Kaimal (2)
            !----------------------------------------
               OtherStates%NFFComp = 1

         CASE(3, 5)
            !----------------------------------------
            !3-component Von Karman (3) or IEC-2
            ! Kaimal (5)
            !----------------------------------------
               OtherStates%NFFComp = 3

         CASE(4)
            !----------------------------------------
            !improved Von Karman
            !----------------------------------------

                  ! Read 2-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of components (should be 3)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of components from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                     RETURN
                  ENDIF
                  OtherStates%NFFComp = Dum_Int4

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Latitude (deg)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading latitude from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                     RETURN
                  ENDIF

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Roughness length (m)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading roughness length from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                     RETURN
                  ENDIF

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Reference height (m) = Z(1) + GridHeight / 2.0

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading reference height from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                     RETURN
                  ENDIF


               DO I = 1,3
                     ! Read 4-byte real. Can't use library routines for this.
                  READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                    ! TI(u, v, w) (%)

                     IF (TmpErrStat /= 0) THEN
                        CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//'TRIM(Num2LStr(I))'//') from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                        RETURN
                     ENDIF
                     TI(I) = Dum_Real4                                                          ! This overwrites the TI read in the summary file

               END DO !I


         CASE (7, 8)
            !----------------------------------------
            ! General Kaimal (7) or  Mann model (8)
            !----------------------------------------

                  ! Read 4-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of bytes in header

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of header records from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                     RETURN
                  ENDIF

                  ! Read 4-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of components

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of data from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                     RETURN
                  ENDIF
                  OtherStates%NFFComp = Dum_Int4


         CASE DEFAULT

            CALL SetErrStat( ErrID_Warn, ' InflowWind does not recognize the full-field turbulence file type ='// &
                        TRIM(Num2LStr(TurbType))//'.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            IF (ErrStat >= AbortErrLev) RETURN

      END SELECT !TurbType


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                                ! delta z (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dz from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         FFZDelt = Dum_Real4
         OtherStates%InvFFZD = 1.0/FFZDelt


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! delta y (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         FFYDelt = Dum_Real4
         OtherStates%InvFFYD = 1.0/FFYDelt

         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! delta x (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dx from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         FFXDelt = Dum_Real4


         ! Read 4-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! half the number of time steps

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of time steps from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         OtherStates%NFFSteps = 2*Dum_Int4


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         OtherStates%MeanFFWS = Dum_Real4
         OtherStates%InvMFFWS = 1.0/OtherStates%MeanFFWS
         OtherStates%FFDTime  = FFXDelt/OtherStates%MeanFFWS
         OtherStates%FFRate   = 1.0/OtherStates%FFDTime


      DO I = 1,3

            ! Read 4-byte real. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                            ! unused variables: zLu, yLu, xLu

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte length scales from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
               RETURN
            ENDIF

      END DO


      DO I = 1,2

         ! Read 4-byte integer. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                             ! unused variables: dummy, random seed

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte integers from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
               RETURN
            ENDIF

      END DO


         ! Read 4-integer real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! nz

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading nz from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         OtherStates%NZGrids  = Dum_Int4
         OtherStates%FFZHWid  = 0.5*FFZDelt*( OtherStates%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 4-integer real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
            RETURN
         ENDIF
         OtherStates%NYGrids  = Dum_Int4
         OtherStates%FFYHWid  = 0.5*FFYDelt*( OtherStates%NYGrids - 1 )


      IF (OtherStates%NFFComp == 3) THEN

         DO I=1,6

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: zLv, yLv, xLv, zLw, yLw, xLw

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte length scales from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

         ENDDO !I

      ENDIF !NFFComp



      IF ( TurbType == 7 ) THEN     ! General Kaimal model

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variable: coherence decay constant

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading coherence decay constant from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: coherence scale parameter in m

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading coherence scale parameter from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

      ELSE IF ( TurbType == 8 ) THEN     ! Mann model

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: shear parameter (gamma), scale length

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,4

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,3

               ! Read 4-integer real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                          ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,3

               ! Read 4-integer real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                          ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, 'Read_Bladed_FF_Header1')
                  RETURN
               ENDIF

         ENDDO !I


      ENDIF !TurbType


      RETURN

   END SUBROUTINE Read_Bladed_FF_Header1
   !====================================================================================================
   SUBROUTINE Read_Bladed_Grids ( OtherStates, CWise, TI, ErrStat, ErrMsg )
   ! This subroutine continues reading OtherStates%UnitWind, starting after the headers have been read.
   ! It reads the Grids and converts the data to un-normalized wind speeds in m/s.
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------

      IMPLICIT                                              NONE

         ! Passed variables

      TYPE(IfW_FFWind_OtherStateType),    INTENT(INOUT)  :: OtherStates
      LOGICAL,                            INTENT(IN   )  :: CWise
      REAL(ReKi),                         INTENT(IN   )  :: TI      (3)    ! turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg

      REAL(ReKi),    PARAMETER                           :: FF_Offset(3) = (/ 1.0, 0.0, 0.0 /)  ! used for "un-normalizing" the data

      INTEGER(IntKi)                                     :: CFirst
      INTEGER(IntKi)                                     :: CLast
      INTEGER(IntKi)                                     :: CStep
      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(IntKi)                                     :: I
      INTEGER(IntKi)                                     :: IC
      INTEGER(IntKi)                                     :: IR
      INTEGER(IntKi)                                     :: IT

      INTEGER(IntKi)                                     :: TmpNumSteps

         ! Temporary variables for error handling

      INTEGER(IntKi)                                     :: TmpErrStat     ! for checking the result of IOSTAT on READ or Open statements
      CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg


      !-------------------------------------------------------------------------------------------------
      ! Generate an informative message. Initialize the ErrStat.
      !-------------------------------------------------------------------------------------------------
         ! This could take a while, so we'll write a message to tell users what's going on:
         
      CALL WrScr( NewLine//'   Reading a '//TRIM( Num2LStr(OtherStates%NYGrids) )//'x'//TRIM( Num2LStr(OtherStates%NZGrids) )//  &
                  ' grid ('//TRIM( Num2LStr(OtherStates%FFYHWid*2) )//' m wide, '// &
                  TRIM( Num2LStr(OtherStates%GridBase) )//' m to '// &
                  TRIM( Num2LStr(OtherStates%GridBase+OtherStates%FFZHWid*2) )//&
                  ' m above ground) with a characteristic wind speed of '//TRIM( Num2LStr(OtherStates%MeanFFWS) )//' m/s. ' )
      ErrMsg   = ""
      ErrStat  =  ErrID_None


      !-------------------------------------------------------------------------------------------------
      ! Allocate space for the FF array
      !-------------------------------------------------------------------------------------------------

      TmpNumSteps = OtherStates%NFFSteps + 1       ! add another step, just in case there is an odd number of steps.

   !bjj: should we reorganize this FFData array so we access the data faster?

      IF ( .NOT. ALLOCATED( OtherStates%FFData ) ) THEN
         CALL AllocAry( OtherStates%FFData, OtherStates%NZGrids,OtherStates%NYGrids,OtherStates%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Read_Bladed_Grids')
         IF ( ErrStat >= AbortErrLev ) RETURN

      ELSE
         IF (SIZE(OtherStates%FFData,1) /= OtherStates%NZGrids .OR. SIZE(OtherStates%FFData,2) /= OtherStates%NYGrids .OR. &
             SIZE(OtherStates%FFData,3) /= OtherStates%NFFComp .OR. SIZE(OtherStates%FFData,3) /= TmpNumSteps ) THEN

               ! Let's make the array the correct size (we should never get here, but you never know)

            DEALLOCATE( OtherStates%FFData )

            CALL AllocAry( OtherStates%FFData, OtherStates%NZGrids,OtherStates%NYGrids,OtherStates%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Read_Bladed_Grids')
               IF ( ErrStat >= AbortErrLev ) RETURN            

         ENDIF !Incorrect size
      ENDIF ! allocated

      !-------------------------------------------------------------------------------------------------
      ! Initialize the data and set column indexing to account for direction of turbine rotation (CWise)
      !-------------------------------------------------------------------------------------------------

      OtherStates%FFData(:,:,:,:) = 0.0                        ! we may have only one component

      IF ( CWise )  THEN
         CFirst    = OtherStates%NYGrids
         CLast     = 1
         CStep     = -1
      ELSE
         CFirst    = 1
         CLast     = OtherStates%NYGrids
         CStep     = 1
      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Loop through all the time steps, reading the data and converting to m/s
      !-------------------------------------------------------------------------------------------------
   !bjj: should we reorganize this FFData array so we access the data faster?

      OtherStates%NFFSteps = TmpNumSteps

   TIME_LOOP:  DO IT=1,TmpNumSteps     ! time (add 1 to see if there is an odd number of grids)

         DO IR=1,OtherStates%NZGrids               ! the rows (vertical)

            DO IC=CFirst,CLast,CStep   ! the columns (lateral)

               DO I=1,OtherStates%NFFComp          ! wind components (U, V, W)

                     ! Get the next integer from the file.
                     ! This is a 2-byte integer, so we can't use the library read routines.
                  READ (OtherStates%UnitWind,IOStat=TmpErrStat)  Dum_Int2
                  IF (TmpErrStat /= 0) THEN
                     IF ( IT == TmpNumSteps ) THEN ! There really were an even number of steps
                        OtherStates%NFFSteps = TmpNumSteps - 1
                        ErrStat  = 0
                        EXIT TIME_LOOP
                     ELSE
                        CALL SetErrStat( ErrID_Fatal, ' Error reading binary data file. '// &
                                    'ic = '//TRIM(Num2LStr(ic))// &
                                    ', ir = '//TRIM(Num2LStr(ir))// &
                                    ', it = '//TRIM(Num2LStr(it))// &
                                    ', nffsteps = '//TRIM(Num2LStr(OtherStates%NFFSteps)), ErrStat, ErrMsg, 'Read_Bladed_Grids')
                        RETURN
                     ENDIF
                  ELSE
                     OtherStates%FFData(IR,IC,I,IT) = OtherStates%MeanFFWS*(FF_Offset(I)+0.00001*TI(I)*Dum_Int2)
                  ENDIF

               END DO !I

            END DO !IC

         END DO !IR

      END DO TIME_LOOP !IT

      IF ( ParamData%Periodic ) THEN
         TmpErrMsg = '   Processed '//TRIM( Num2LStr( OtherStates%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( OtherStates%FFRate ) )//'-Hz full-field data (period of '// &
                    TRIM( Num2LStr( OtherStates%FFDTime*OtherStates%NFFSteps ) )//' seconds).'
         
      ELSE
         TmpErrMsg= '   Processed '//TRIM( Num2LStr( OtherStates%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( OtherStates%FFRate ) )//'-Hz full-field data ('// &
                    TRIM( Num2LStr( OtherStates%FFDTime*( OtherStates%NFFSteps - 1 ) ) )//' seconds).'
      ENDIF
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )
      !CALL SetErrStat( ErrID_Info, TmpErrMsg, ErrStat, ErrMsg, 'Read_Bladed_Grids' )
      
      

   END SUBROUTINE Read_Bladed_Grids
   !====================================================================================================
   SUBROUTINE Read_FF_Tower( OtherStates, ErrStat, ErrMsg )
   ! This subroutine reads the binary tower file that corresponds with the Bladed-style FF binary file.
   ! The FF grid must be read before this subroutine is called! (many checks are made to ensure the
   ! files belong together)
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------

         ! Passed Variables:

      TYPE(IfW_FFWind_OtherStateType),    INTENT(INOUT)  :: OtherStates    ! unit number for the wind file
!      TYPE(IfW_FFWind_ParameterType),     INTENT(IN   )  :: ParamData      ! Parameters
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        ! error status return value (0=no error; non-zero is error)
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         ! a message for errors that occur

         ! Local Variables:

      REAL(SiKi)                                         :: Dum_Real4      ! dummy 4-byte real number
      INTEGER(B2Ki)                                      :: Dum_Int2       ! dummy 2-byte integer
      INTEGER(B4Ki)                                      :: Dum_Int4       ! dummy 4-byte integer

      INTEGER(IntKi)                                     :: IC             ! loop counter for wind components
      INTEGER(IntKi)                                     :: IT             ! loop counter for time
      INTEGER(IntKi)                                     :: IZ             ! loop counter for z

      REAL(ReKi),    PARAMETER                           :: TOL = 1E-4     ! tolerence for wind file comparisons

      REAL(ReKi),    PARAMETER                           :: FF_Offset(3) = (/ 1.0, 0.0, 0.0 /)  ! used for "un-normalizing" the data
      REAL(SiKi)                                         :: TI       (3)   ! scaling values for "un-normalizing the data" [approx. turbulence intensities of the wind components]


         ! Temporary Error Handling

      INTEGER(IntKi)                                     :: TmpErrStat     ! IOSTAT value.
      CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg

      !-------------------------------------------------------------------------------------------------
      ! Initialization
      !-------------------------------------------------------------------------------------------------

      ErrMsg   = ''
      ErrStat  = ErrID_None

      OtherStates%NTGrids = 0

      IF ( OtherStates%NFFComp /= 3 ) THEN
         ErrMsg   = ' Error: Tower binary files require 3 wind components.'
         ErrStat  = ErrID_Fatal  
         RETURN
      ELSE
         ErrStat = ErrID_None
         ErrMsg  = ""
      ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Open the file
      !-------------------------------------------------------------------------------------------------

      CALL OpenBInpFile (OtherStates%UnitWind, TRIM(ParamData%WindFileName), TmpErrStat, TmpErrMsg)
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Read_FF_Tower' )
         IF (ErrStat >= AbortErrLev) RETURN

      !-------------------------------------------------------------------------------------------------
      ! Read the header information and check that it's compatible with the FF Bladed-style binary
      ! parameters already read.
      !-------------------------------------------------------------------------------------------------
            ! This is a 4-byte real, so we can't use the library read routines.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! dz, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dz in the binary tower file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*OtherStates%InvFFZD-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Resolution in the FF binary file does not match the tower file.', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! dx, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dx in the binary tower file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*OtherStates%InvMFFWS/OtherStates%FFDTime-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Time resolution in the FF binary file does not match the tower file.', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! Zmax, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading Zmax in the binary tower file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4/OtherStates%GridBase-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Height in the FF binary file does not match the tower file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF


            ! This is a 4-byte integer, so we can't use the library read routines.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                ! NumOutSteps [4-byte INTEGER]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading NumOutSteps in the binary tower file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF

            IF ( Dum_Int4 /= OtherStates%NFFSteps ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Number of time steps in the FF binary file does not match the tower file.', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF


            ! This is a 4-byte integer, so we can't use the library read routines.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                ! NumZ      [4-byte INTEGER]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading NumZ in the binary tower file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF
            OtherStates%NTGrids = Dum_Int4


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! UHub      [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading UHub in the binary tower file "'//TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*OtherStates%InvMFFWS - 1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Mean wind speed in the FF binary file does not match the tower file.', ErrStat, ErrMsg, 'Read_FF_Tower' )
               OtherStates%NTGrids  = 0
               RETURN
            ENDIF


         DO IC=1,3
               ! Read the TI values fromthe tower file: 4-byte reals.
               
               !bjj: not sure you can call this routine to read from a binary file...
            !CALL ReadVar( OtherStates%UnitWind, TRIM(ParamData%WindFileName), TI(IC), 'TI('//TRIM(Num2LStr(IC))//')', 'TI value for u,v, or w', TmpErrStat, TmpErrMsg )
            !IF (TmpErrStat /= ErrID_None) THEN
            !   OtherStates%NTGrids  = 0
            !   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Read_FF_Tower' )
            !   IF (ErrStat >= AbortErrLev) RETURN
            !ENDIF
            !
            READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   TI(IC)               ! TI(u), TI(v), TI(w)  [4-byte REAL]
            
            IF (TmpErrStat /= 0) THEN
               OtherStates%NTGrids  = 0
               CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(IC))//') in the binary tower file "' &
                               //TRIM( ParamData%WindFileName )//'."', ErrStat, ErrMsg, 'Read_FF_Tower' )
               RETURN
            ENDIF
         
         END DO

      !----------------------------------------------------------------------------------------------
      ! Allocate arrays for the tower points
      !----------------------------------------------------------------------------------------------

         IF ( OtherStates%NTGrids > 0 ) THEN

            IF ( .NOT. ALLOCATED( OtherStates%FFTower ) ) THEN
               CALL AllocAry( OtherStates%FFTower, OtherStates%NFFComp, OtherStates%NTGrids, OtherStates%NFFSteps, &
                  'Tower wind data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Read_FF_Tower' )
               IF (ErrStat >= AbortErrLev) RETURN
               
            ELSE
               ! Check sizes here!
            ENDIF

         ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Read the 16-bit time-series data and scale it to 32-bit reals
      !-------------------------------------------------------------------------------------------------

         ! Loop through time.

         DO IT=1,OtherStates%NFFSteps

            DO IZ=1,OtherStates%NTGrids         ! If NTGrids<1, there are no tower points & FFTower is not allocated

               ! Ytower     = 0               ! Lateral location of the tower data point, in m relative to tower centerline
               ! Ztower(IZ) = Z1 - (IZ-1)*dz  ! Vertical location of tower data point, in m relative to ground

               DO IC=1,OtherStates%NFFComp   ! number of wind components

                     ! Read in the 2-byte integer. Can't use library read routines for this.
                  READ (OtherStates%UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading binary tower data file. it = '//TRIM(Num2LStr(it))// &
                                    ', nffsteps = '//TRIM(Num2LStr(OtherStates%NFFSteps)), ErrStat, ErrMsg, 'Read_FF_Tower' )
                     OtherStates%NTGrids  = 0
                     RETURN
                  ENDIF

                  OtherStates%FFTower(IC,IZ,IT) = OtherStates%MeanFFWS*(FF_Offset(IC)+0.00001*TI(IC)*Dum_Int2)   ! wind-component scaled to m/s

               ENDDO !IC

            ENDDO ! IZ


         ENDDO ! IT

      !-------------------------------------------------------------------------------------------------
      ! Close the file
      !-------------------------------------------------------------------------------------------------
      CLOSE ( OtherStates%UnitWind )

      TmpErrMsg = '   Processed '//TRIM( Num2LStr(OtherStates%NFFSteps) )//' time steps of '// &
            TRIM( Num2LStr(OtherStates%NTGrids) )//'x1 tower data grids.'
      
      !CALL SetErrStat( ErrID_Info, ErrMsgLcl, ErrStat, ErrMsg, 'Read_FF_Tower' )
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )
      
      RETURN

   END SUBROUTINE Read_FF_Tower
END SUBROUTINE IfW_FFWind_Init
!====================================================================================================


!====================================================================================================
SUBROUTINE IfW_FFWind_CalcOutput(Time,    InData,        ParamData,                       &
                           ContStates,    DiscStates,    ConstrStates,     OtherStates,   &
                           OutData,       ErrStat,       ErrMsg)
   !-------------------------------------------------------------------------------------------------
   ! This routine acts as a wrapper for the GetWindSpeed routine. It steps through the array of input
   ! positions and calls the GetWindSpeed routine to calculate the velocities at each point.
   !
   ! There are inefficiencies in how this set of routines is coded, but that is a problem for another
   ! day. For now, it merely needs to be functional. It can be fixed up and made all pretty later.
   !
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !-------------------------------------------------------------------------------------------------

   IMPLICIT                                                 NONE

      ! Passed Variables
   REAL(DbKi),                            INTENT(IN   )  :: Time           ! time from the start of the simulation
   TYPE(IfW_FFWind_InputType),            INTENT(IN   )  :: InData         ! Input Data
   TYPE(IfW_FFWind_ParameterType),        INTENT(IN   )  :: ParamData      ! Parameters
   TYPE(IfW_FFWind_ContinuousStateType),  INTENT(IN   )  :: ContStates     ! Continuous States  (unused)
   TYPE(IfW_FFWind_DiscreteStateType),    INTENT(IN   )  :: DiscStates     ! Discrete States    (unused)
   TYPE(IfW_FFWind_ConstraintStateType),  INTENT(IN   )  :: ConstrStates   ! Constraint States  (unused)
   TYPE(IfW_FFWind_OtherStateType),       INTENT(INOUT)  :: OtherStates    ! Other State data   (storage for the main data)
   TYPE(IfW_FFWind_OutputType),           INTENT(  OUT)  :: OutData        ! Initial output

      ! Error handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat        ! error status
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg         ! The error message

      ! local variables
   INTEGER(IntKi)                                        :: NumPoints      ! Number of points specified by the InData%Position array

      ! local counters
   INTEGER(IntKi)                                        :: PointNum       ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg      ! temporary error message


      !-------------------------------------------------------------------------------------------------
      ! Check that the module has been initialized.
      !-------------------------------------------------------------------------------------------------

      IF ( .NOT. ParamData%Initialized ) THEN
         ErrMsg   = ' Initialialize the FFWind module before calling its subroutines.'
         ErrStat  = ErrID_Fatal
         RETURN
      ELSE
         ErrStat = ErrID_None
         ErrMsg   = ''
      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------

   TmpErrStat  = ErrID_None

      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   TmpErrMsg   = ""
   NumPoints   =  SIZE(InData%Position,2)
      ! Allocate Velocity output array
   CALL AllocAry( OutData%Velocity, 3, NumPoints, "Velocity matrix at timestep", TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'CalcOutput' )
      IF (ErrStat >= AbortErrLev) RETURN


      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      OutData%Velocity(:,PointNum) = FF_Interp(Time,InData%Position(:,PointNum),ParamData,OtherStates,TmpErrStat,TmpErrMsg)

         ! Error handling
      IF (TmpErrStat /= ErrID_None) THEN  !  adding this so we don't have to convert numbers to strings every time
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, "IfW_FFWind:CalcOutput [position=("//   &
                                                      TRIM(Num2LStr(InData%Position(1,PointNum)))//", "// &
                                                      TRIM(Num2LStr(InData%Position(2,PointNum)))//", "// &
                                                      TRIM(Num2LStr(InData%Position(3,PointNum)))//")]" )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF

   ENDDO


   RETURN

CONTAINS
   !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
   FUNCTION FF_Interp(Time, Position, ParamData, OtherStates, ErrStat, ErrMsg)
   !    This function is used to interpolate into the full-field wind array or tower array if it has
   !    been defined and is necessary for the given inputs.  It receives X, Y, Z and
   !    TIME from the calling routine.  It then computes a time shift due to a nonzero X based upon
   !    the average windspeed.  The modified time is used to decide which pair of time slices to interpolate
   !    within and between.  After finding the two time slices, it decides which four grid points bound the
   !    (Y,Z) pair.  It does a bilinear interpolation for each time slice. Linear interpolation is then used
   !    to interpolate between time slices.  This routine assumes that X is downwind, Y is to the left when
   !    looking downwind and Z is up.  It also assumes that no extrapolation will be needed.
   !
   !    If tower points are used, it assumes the velocity at the ground is 0.  It interpolates between
   !    heights and between time slices, but ignores the Y input.
   !
   !    11/07/1994 - Created by M. Buhl from the original TURBINT.
   !    09/25/1997 - Modified by M. Buhl to use f90 constructs and new variable names.  Renamed to FF_Interp.
   !    09/23/2009 - Modified by B. Jonkman to use arguments instead of modules to determine time and position.
   !                 Height is now relative to the ground
   !   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !
   !----------------------------------------------------------------------------------------------------

      IMPLICIT                                              NONE

      REAL(DbKi),                         INTENT(IN   )  :: Time
      REAL(ReKi),                         INTENT(IN   )  :: Position(3)    ! takes the place of XGrnd, YGrnd, ZGrnd
      TYPE(IfW_FFWind_ParameterType),     INTENT(IN   )  :: ParamData      ! Parameters
      TYPE(IfW_FFWind_OtherStateType),    INTENT(INOUT)  :: OtherStates    ! Other State data   (storage for the main data)
      REAL(ReKi)                                         :: FF_Interp(3)   ! The U, V, W velocities

      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg

         ! Local Variables:

      REAL(ReKi)                                         :: TimeShifted
      REAL(ReKi),PARAMETER                               :: Tol = 1.0E-3   ! a tolerance for determining if two reals are the same (for extrapolation)
      REAL(ReKi)                                         :: W_YH_Z
      REAL(ReKi)                                         :: W_YH_ZH
      REAL(ReKi)                                         :: W_YH_ZL
      REAL(ReKi)                                         :: W_YL_Z
      REAL(ReKi)                                         :: W_YL_ZH
      REAL(ReKi)                                         :: W_YL_ZL
      REAL(ReKi)                                         :: Wnd(2)
      REAL(ReKi)                                         :: T
      REAL(ReKi)                                         :: TGRID
      REAL(ReKi)                                         :: Y
      REAL(ReKi)                                         :: YGRID
      REAL(ReKi)                                         :: Z
      REAL(ReKi)                                         :: ZGRID

      INTEGER(IntKi)                                     :: IDIM
      INTEGER(IntKi)                                     :: IG
      INTEGER(IntKi)                                     :: IT
      INTEGER(IntKi)                                     :: ITHI
      INTEGER(IntKi)                                     :: ITLO
      INTEGER(IntKi)                                     :: IYHI
      INTEGER(IntKi)                                     :: IYLO
      INTEGER(IntKi)                                     :: IZHI
      INTEGER(IntKi)                                     :: IZLO

      LOGICAL                                            :: OnGrid

      !-------------------------------------------------------------------------------------------------
      ! Initialize variables
      !-------------------------------------------------------------------------------------------------

      FF_Interp(:)          = 0.0                         ! the output velocities (in case OtherStates%NFFComp /= 3)
      Wnd(:)                = 0.0                         ! just in case we're on an end point

      ErrStat               = ErrID_None
      ErrMsg               = ""
      
      !-------------------------------------------------------------------------------------------------
      ! Find the bounding time slices.
      !-------------------------------------------------------------------------------------------------

      ! Perform the time shift.  At time=0, a point half the grid width downstream (OtherStates%FFYHWid) will index into the zero time slice.
      ! If we did not do this, any point downstream of the tower at the beginning of the run would index outside of the array.
      ! This all assumes the grid width is at least as large as the rotor.  If it isn't, then the interpolation will not work.


      TimeShifted = TIME + ( OtherStates%InitXPosition - Position(1) )*OtherStates%InvMFFWS    ! in distance, X: InputInfo%Position(1) - OtherStates%InitXPosition - TIME*OtherStates%MeanFFWS


      IF ( ParamData%Periodic ) THEN ! translate TimeShifted to ( 0 <= TimeShifted < OtherStates%TotalTime )

         TimeShifted = MODULO( TimeShifted, OtherStates%TotalTime )
             ! If TimeShifted is a very small negative number, modulo returns the incorrect value due to internal rounding errors.
             ! See bug report #471
         IF (TimeShifted == OtherStates%TotalTime) TimeShifted = 0.0_ReKi

         TGRID = TimeShifted*OtherStates%FFRate
         ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
         T     = TGRID - ITLO             ! a value between 0 and 1 that indicates a relative location between ITLO and ITHI

         ITLO = ITLO + 1
         IF ( ITLO == OtherStates%NFFSteps ) THEN
            ITHI = 1
         ELSE
            ITHI = ITLO + 1
         ENDIF


      ELSE

         TGRID = TimeShifted*OtherStates%FFRate
         ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
         T     = TGRID - ITLO             ! a value between 0 and 1 that indicates a relative location between ITLO and ITHI

         ITLO = ITLO + 1                  ! add one since our grids start at 1, not 0
         ITHI = ITLO + 1

         IF ( ITLO >= OtherStates%NFFSteps .OR. ITLO < 1 ) THEN
            IF ( ITLO == OtherStates%NFFSteps  ) THEN
               ITHI = ITLO
               IF ( T <= TOL ) THEN ! we're on the last point
                  T = 0.0
               ELSE  ! We'll extrapolate one dt past the last value in the file
                  ITLO = ITHI - 1
               ENDIF
            ELSE
               ErrMsg   = ' Error: FF wind array was exhausted at '//TRIM( Num2LStr( REAL( TIME,   ReKi ) ) )// &
                          ' seconds (trying to access data at '//TRIM( Num2LStr( REAL( TimeShifted, ReKi ) ) )//' seconds).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ENDIF

      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Find the bounding rows for the Z position. [The lower-left corner is (1,1) when looking upwind.]
      !-------------------------------------------------------------------------------------------------

      ZGRID = ( Position(3) - OtherStates%GridBase )*OtherStates%InvFFZD

      IF (ZGRID > -1*TOL) THEN
         OnGrid = .TRUE.

         IZLO = INT( ZGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
         IZHI = IZLO + 1

         Z = ZGRID - ( IZLO - 1 )            ! a value between 0 and 1 that indicates a relative location between IZLO and IZHI

         IF ( IZLO < 1 ) THEN
            IF ( IZLO == 0 .AND. Z >= 1.0-TOL ) THEN
               Z    = 0.0
               IZLO = 1
            ELSE
               ErrMsg   = ' Error: FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                          TRIM(Num2LStr(Position(3)))//' m is below the grid).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ELSEIF ( IZLO >= OtherStates%NZGrids ) THEN
            IF ( IZLO == OtherStates%NZGrids .AND. Z <= TOL ) THEN
               Z    = 0.0
               IZHI = IZLO                   ! We're right on the last point, which is still okay
            ELSE
               ErrMsg   = ' Error: FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                          TRIM(Num2LStr(Position(3)))//' m is above the grid).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ENDIF

      ELSE

         OnGrid = .FALSE.  ! this is on the tower

         IF ( OtherStates%NTGrids < 1 ) THEN
            ErrMsg   = ' Error: FF wind array boundaries violated. Grid too small in Z direction '// &
                       '(height (Z='//TRIM(Num2LStr(Position(3)))//' m) is below the grid and no tower points are defined).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF

         IZLO = INT( -1.0*ZGRID ) + 1            ! convert REAL to INTEGER, then add one since our grids start at 1, not 0


         IF ( IZLO >= OtherStates%NTGrids ) THEN  !our dz is the difference between the bottom tower point and the ground
            IZLO = OtherStates%NTGrids

            Z    = 1.0 - Position(3) / (OtherStates%GridBase - (IZLO-1)/OtherStates%InvFFZD) !check that this isn't 0
         ELSE
            Z    = ABS(ZGRID) - (IZLO - 1)
         ENDIF
         IZHI = IZLO + 1

      ENDIF


      IF ( OnGrid ) THEN      ! The tower points don't use this

         !-------------------------------------------------------------------------------------------------
         ! Find the bounding columns for the Y position. [The lower-left corner is (1,1) when looking upwind.]
         !-------------------------------------------------------------------------------------------------

            YGRID = ( Position(2) + OtherStates%FFYHWid )*OtherStates%InvFFYD    ! really, it's (Position(2) - -1.0*OtherStates%FFYHWid)

            IYLO = INT( YGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
            IYHI = IYLO + 1

            Y    = YGRID - ( IYLO - 1 )         ! a value between 0 and 1 that indicates a relative location between IYLO and IYHI

            IF ( IYLO >= OtherStates%NYGrids .OR. IYLO < 1 ) THEN
               IF ( IYLO == 0 .AND. Y >= 1.0-TOL ) THEN
                  Y    = 0.0
                  IYLO = 1
               ELSE IF ( IYLO == OtherStates%NYGrids .AND. Y <= TOL ) THEN
                  Y    = 0.0
                  IYHI = IYLO                   ! We're right on the last point, which is still okay
               ELSE
                  ErrMsg   = ' Error FF wind array boundaries violated: Grid too small in Y direction. Y='// &
                             TRIM(Num2LStr(Position(2)))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*OtherStates%FFYHWid))// &
                             ', '//TRIM(Num2LStr(OtherStates%FFYHWid))//']'
                  ErrStat = ErrID_Fatal         ! we don't return anything
                  RETURN
               ENDIF
            ENDIF

         !-------------------------------------------------------------------------------------------------
         ! Interpolate on the grid
         !-------------------------------------------------------------------------------------------------

         DO IDIM=1,OtherStates%NFFComp       ! all the components

            IT = ITLO            ! Start using the ITLO slice

            DO IG=1,2            ! repeat for 2 time slices (by changing the value of IT. note that we can't loop from IXLO to IXHI because they could be OtherStates%NFFSteps and 1 respectively)

               !-------------------------------------------------------------------------------------------
               ! Get the wind velocity values for the four corners of the grid for this time.
               !-------------------------------------------------------------------------------------------

               W_YL_ZL = OtherStates%FFData( IZLO, IYLO, IDIM, IT )
               W_YL_ZH = OtherStates%FFData( IZHI, IYLO, IDIM, IT )
               W_YH_ZL = OtherStates%FFData( IZLO, IYHI, IDIM, IT )
               W_YH_ZH = OtherStates%FFData( IZHI, IYHI, IDIM, IT )


               !-------------------------------------------------------------------------------------------
               ! Interpolate within the grid for this time.
               !-------------------------------------------------------------------------------------------

               W_YL_Z  = ( W_YL_ZH - W_YL_ZL )*Z + W_YL_ZL
               W_YH_Z  = ( W_YH_ZH - W_YH_ZL )*Z + W_YH_ZL
               Wnd(IG) = ( W_YH_Z  - W_YL_Z  )*Y + W_YL_Z

               IT = ITHI            ! repeat for the using the ITHI slice

            END DO !IG

            !----------------------------------------------------------------------------------------------
            ! Interpolate between the two times.
            !----------------------------------------------------------------------------------------------

            FF_Interp(IDIM) = ( Wnd(2) - Wnd(1) ) * T + Wnd(1)    ! interpolated velocity

         END DO !IDIM

      ELSE

      !-------------------------------------------------------------------------------------------------
      ! Interpolate on the tower array
      !-------------------------------------------------------------------------------------------------

         DO IDIM=1,OtherStates%NFFComp    ! all the components

            IT = ITLO            ! Start using the ITLO slice

            DO IG=1,2            ! repeat for 2 time slices (by changing the value of IT. note that we can't loop from IXLO to IXHI because they could be OtherStates%NFFSteps and 1 respectively)

               !-------------------------------------------------------------------------------------------
               ! Get the wind velocity values for the two corners of the grid for this time.
               !-------------------------------------------------------------------------------------------

               W_YH_ZL = OtherStates%FFTower( IDIM, IZLO, IT )

               IF ( IZHI > OtherStates%NTGrids ) THEN
                  W_YH_ZH = 0.0
               ELSE
                  W_YH_ZH = OtherStates%FFTower( IDIM, IZHI, IT )
               ENDIF


               !-------------------------------------------------------------------------------------------
               ! Interpolate within the grid for this time.
               !-------------------------------------------------------------------------------------------

               Wnd(IG) = ( W_YH_ZH - W_YH_ZL )*Z + W_YH_ZL

               IT = ITHI            ! repeat for the using the ITHI slice

            END DO !IG

            !----------------------------------------------------------------------------------------------
            ! Interpolate between the two times.
            !----------------------------------------------------------------------------------------------

            FF_Interp(IDIM) = ( Wnd(2) - Wnd(1) ) * T + Wnd(1)    ! interpolated velocity

         END DO !IDIM

      ENDIF ! OnGrid

      RETURN

   END FUNCTION FF_Interp
END SUBROUTINE IfW_FFWind_CalcOutput


!====================================================================================================
SUBROUTINE IfW_FFWind_End( InData,     ParamData,                                &
                           ContStates, DiscStates, ConstrStates,  OtherStates,   &
                           OutData,                                              &
                           ErrStat,    ErrMsg)
   !-------------------------------------------------------------------------------------------------
   !  This subroutine cleans up any data that is still allocated.  The (possibly) open files are
   !  closed in InflowWindMod.
   !
   !  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !-------------------------------------------------------------------------------------------------

      ! Passed Variables
   TYPE(IfW_FFWind_InputType),            INTENT(INOUT)  :: InData         ! Initialized input data variable
   TYPE(IfW_FFWind_ParameterType),        INTENT(INOUT)  :: ParamData      ! Parameters
   TYPE(IfW_FFWind_ContinuousStateType),  INTENT(INOUT)  :: ContStates     ! Continuous States  (unused)
   TYPE(IfW_FFWind_DiscreteStateType),    INTENT(INOUT)  :: DiscStates     ! Discrete States    (unused)
   TYPE(IfW_FFWind_ConstraintStateType),  INTENT(INOUT)  :: ConstrStates   ! Constraint States  (unused)
   TYPE(IfW_FFWind_OtherStateType),       INTENT(INOUT)  :: OtherStates    ! Other State data   (storage for the main data)
   TYPE(IfW_FFWind_OutputType),           INTENT(INOUT)  :: OutData        ! Initial output


      ! Error Handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat        ! determines if an error has been encountered
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg         ! Message about errors


      ! Local Variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg      ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None


      ! Destroy the input data

   CALL IfW_FFWind_DestroyInput(       InData,        TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! Destroy parameter data

   CALL IfW_FFWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! Destroy the state data

   CALL IfW_FFWind_DestroyContState(   ContStates,    TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF

   CALL IfW_FFWind_DestroyDiscState(   DiscStates,    TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF

   CALL IfW_FFWind_DestroyConstrState( ConstrStates,  TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF

   CALL IfW_FFWind_DestroyOtherState(  OtherStates,   TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! Destroy the output data

   CALL IfW_FFWind_DestroyOutput(      OutData,       TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! flag as uninitialized
   ParamData%Initialized = .FALSE.


END SUBROUTINE IfW_FFWind_End

!====================================================================================================
!====================================================================================================
!====================================================================================================



!====================================================================================================
! The following are generic routines required by the framework.
!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!====================================================================================================
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_FFWind_UpdateStates( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_FFWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_FFWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_FFWind_ContinuousStateType),  INTENT(INOUT)  :: x           ! Input: Continuous states at Time;
                                                                           ! Output: Continuous states at Time + Interval
      TYPE(IfW_FFWind_DiscreteStateType),    INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                           ! Output: Discrete states at Time  + Interval
      TYPE(IfW_FFWind_ConstraintStateType),  INTENT(INOUT)  :: z           ! Input: Initial guess of constraint states at Time;
                                                                           ! Output: Constraint states at Time
      TYPE(IfW_FFWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(IfW_FFWind_ContinuousStateType)                  :: dxdt        ! Continuous state derivatives at Time
      TYPE(IfW_FFWind_ConstraintStateType)                  :: z_Residual  ! Residual of the constraint state equations (Z)

      INTEGER(IntKi)                                        :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                                :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""



         ! Solve for the constraint states (z) here:

         ! Check if the z guess is correct and update z with a new guess.
         ! Iterate until the value is within a given tolerance.

      CALL IfW_FFWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_FFWind_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      ENDIF

      ! DO WHILE ( z_Residual% > tolerance )
      !
      !  z =
      !
      !  CALL IfW_FFWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      !  IF ( ErrStat >= AbortErrLev ) THEN
      !     CALL IfW_FFWind_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
      !     ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
      !     RETURN
      !  ENDIF
      !
      ! END DO


         ! Destroy z_Residual because it is not necessary for the rest of the subroutine:

      CALL IfW_FFWind_DestroyConstrState( z_Residual, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



         ! Get first time derivatives of continuous states (dxdt):

      CALL IfW_FFWind_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_FFWind_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      ENDIF


         ! Update discrete states:
         !   Note that xd [discrete state] is changed in IfW_FFWind_UpdateDiscState(), so IfW_FFWind_CalcOutput(),
         !   IfW_FFWind_CalcContStateDeriv(), and IfW_FFWind_CalcConstrStates() must be called first (see above).

      CALL IfW_FFWind_UpdateDiscState(Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_FFWind_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      ENDIF


         ! Integrate (update) continuous states (x) here:

      !x = function of dxdt and x


         ! Destroy dxdt because it is not necessary for the rest of the subroutine

      CALL IfW_FFWind_DestroyContState( dxdt, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



END SUBROUTINE IfW_FFWind_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_FFWind_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_FFWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_FFWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_FFWind_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_FFWind_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(IfW_FFWind_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(IfW_FFWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(IfW_FFWind_ContinuousStateType),  INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Compute the first time derivatives of the continuous states here:

      dxdt%DummyContState = 0.0_ReKi


END SUBROUTINE IfW_FFWind_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_FFWind_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_FFWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_FFWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_FFWind_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_FFWind_DiscreteStateType),    INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                           !   Output: Discrete states at Time + Interval
      TYPE(IfW_FFWind_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(IfW_FFWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE IfW_FFWind_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_FFWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_FFWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_FFWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_FFWind_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_FFWind_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(IfW_FFWind_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(IfW_FFWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(IfW_FFWind_ConstraintStateType),  INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using
                                                                           ! the input values described above
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0.0_ReKi

END SUBROUTINE IfW_FFWind_CalcConstrStateResidual



!====================================================================================================
END MODULE IfW_FFWind