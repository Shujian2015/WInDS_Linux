!**********************************************************************************************************************************
! FAST_Prog.f90, FAST_IO.f90, FAST_Types.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
!
!    This file is part of FAST.
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
! File last committed: $Date: 2015-03-31 14:15:46 -0600 (Tue, 31 Mar 2015) $
! (File) Revision #: $Rev: 954 $
! URL: $HeadURL: https://windsvn.nrel.gov/FAST/branches/BJonkman/Source/FAST_Mods.f90 $
!**********************************************************************************************************************************
MODULE FAST_ModTypes

   USE NWTC_Library
   USE FAST_Types

   TYPE(ProgDesc), PARAMETER :: FAST_Ver    = &
                                ProgDesc( 'FAST', 'v8.10.00a-bjj', '31-Mar-2015' ) ! The version number of this module
         
   !..................................................................
   
   INTEGER(IntKi), PARAMETER :: Type_LandBased          = 1
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Fixed     = 2
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Floating  = 3
   
   INTEGER(IntKi), PARAMETER :: STATE_CURR              = 1
   INTEGER(IntKi), PARAMETER :: STATE_PRED              = 2
   
         
   INTEGER(IntKi), PARAMETER :: SizeJac_ED_HD  = 12
   
#if defined COMPILE_SIMULINK || defined COMPILE_LABVIEW
   !bjj: 2015-03-03: not sure this is still necessary...
   INTEGER(B2Ki),  PARAMETER :: OutputFileFmtID = FileFmtID_WithTime            ! We cannot guarantee the output time step is constant in binary files
#else
   INTEGER(B2Ki),  PARAMETER :: OutputFileFmtID = FileFmtID_WithoutTime         ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)
#endif  

   LOGICAL,        PARAMETER :: GenerateAdamsModel = .FALSE.



END MODULE FAST_ModTypes
!=======================================================================

