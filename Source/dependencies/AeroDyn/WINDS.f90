!**********************************************************************************************************************************
! Wake Induced Dynamics Simulator (WInDS)
!==================================================================================================================================
!
! The work is based on the WInDS (Matlab code) by Dr.Thomas Sebastian, updated by Dr.Matthew Lackner, Nathaniel DeVelder, Evan Gaertner, David Lenz and Andrew Sciotti.
!
! REFERENCES:
! 1. Sebastian, Thomas, "The Aerodynamics and Near Wake of an Offshore Floating Horizontal Axis Wind Turbine" (2012). Dissertations. Paper 516.
! 2. deVelder, Nathaniel B., "Free Wake Potential Flow Vortex Wind Turbine Modeling: Advances in Parallel Processing and Integration of Ground Effects". 
!    Masters Theses 1896 - February 2014. Paper 1176.
! 3. Gaertner, Evan M., "Modeling Dynamic Stall for a Free Vortex Wake Model of a Floating Offshore Wind Turbine". Masters Theses - September 2014.
! 4. Sebastian, T., and M. A. Lackner. "Development of a free vortex wake method code for offshore floating wind turbines." 
!    Renewable Energy 46 (2012): 269-275.
! 5. Sebastian, T., and M. A. Lackner. "Characterization of the unsteady aerodynamics of offshore floating wind turbines." 
!    Wind Energy 16.3 (2013): 339-352.
! 6. Sebastian, Thomas, and Matthew Lackner. "Analysis of the Induction and Wake Evolution of an Offshore Floating Wind Turbine." 
!    Energies (19961073) 5.4 (2012).
! 7. Lackner, Matthew A., Nathaniel deVelder, and Thomas Sebastian. "On 2D and 3D potential flow models of upwind wind turbine tower interference." 
!    Computers & Fluids 71 (2013): 375-379.
! 8. David Lenz. Implementation Cutoff & Freeze.
! 9. Jonkman, B., et al. "NWTC Programmer’s Handbook: A Guide for Software Development Within the FAST Computer-Aided Engineering Tool." 
!    NREL/TP-Draft Version, Golden, CO: National Renewable Energy Laboratory (2012).
! 10. Algorithmic Outline of Unsteady Aerodynamics (AERODYN) Modules. Project WE-201103. FINAL REPORT. September 2, 2011.
!
!.................................................................................................................................
!
! Fortran code is written by Shujian Liu (shujian.liu@hotmail.com)
!
! Last edited Aug 6, 2015
!==================================================================================================================================
! A few notes:
!   1) Data type: 5D array -> (Dimension index(1~3), Time index(1~nt), Timestep stored(1~nt_p),  Radial index(1~ns/nst), Blade index(1~nb))
!   2) Angles are in radian, not degree.
!==================================================================================================================================

MODULE WINDS

   USE NWTC_Library 
   USE AeroDyn_Types
   USE WINDS_IO          ! Handle input and output
   USE WINDS_Accelerate  ! Acceleration of Biot-Savart Law 
   USE WINDS_DS          ! LB Dynamic stall
   USE WINDS_Library     ! Some subroutines for debug
   USE WINDS_Treecode    ! Treecode algorithm to speedup 
   

   
   IMPLICIT        NONE

      ! ..... Public Subroutines ............
      
   PUBLIC :: WINDS_SetParameters   ! Set the parameters of WInDS(FVM) module
   PUBLIC :: WINDS_Allocate        ! Allocate necessary memory for the WInDS(FVM) module
   PUBLIC :: WINDS_Kinematics      ! Computes the station locations of each blade in the inertial coordinate system
   PUBLIC :: WINDS_Velocity        ! Computes the velocity contributions due to turbine and platform motions and freestream flow in the inertial and blade coordinate systems.
   PUBLIC :: WINDS_FVMInitial      ! Calculate aerodynamic loads of the first timestep by BEM
   PUBLIC :: WINDS_FVM             ! Use free vortex wake method to calculate the aerodynamic loads after first timestep
   


CONTAINS
!==================================================================================================================================
SUBROUTINE WINDS_SetParameters( p, O, ErrStat, ErrMess )
! This subroutine sets the parameters, based on the data stored in InputFileData
! Called from: AD_Init (in AeroDyn.f90)
! (~ part of WInDS.m, constants.m, NRELrotor.m)
!
!
! ************ Pseudo code *********
!  Basic parameters
!  Constants with Ramasamy-Leishman vortex model
!  Current date and time
!  Parameters for blade element
!  Cutoff parameters
!....................................................................

   IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_ParameterType),   INTENT(INOUT)    :: p              ! The module's parameter data
   TYPE(AD_OtherStateType),  INTENT(INOUT)    :: O              ! Other/optimization states   
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        ! The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMess         ! The error message, if an error occurred

      ! Local variables
   INTEGER(IntKi)      :: IBlade     ! Index for blade
   INTEGER(IntKi)      :: IElement   ! Index for blade station
   INTEGER(IntKi)      :: IElement2  ! Index for blade trailing nodes

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS 
      
   CHARACTER(LEN=8)    :: DATE_WINDS
   CHARACTER(LEN=10)   :: TIME_WINDS

   
   ! ....................................................................................................
      ! Basic parameters for simulation
   p%FVM%NB         =   p%NumBl           ! Number of blades
   p%FVM%NST        =   p%Element%NElm + 1   ! Number of trailing nodes (number of station +1) 
   p%FVM%NS         =   p%Element%NElm       ! Number of shed nodes (stations)      
   
   NST              =   p%FVM%NST   ! just for convenience
   NS               =   p%FVM%NS    ! just for convenience 
   

   p%FVM%DT_WINDS   =   p%FVM%DT_RATIO  * p%DtAero             ! Timestep duration in WInDS
   
   IF ( p%FVM%DT_WINDS >= 0.2 .OR. p%FVM%DT_WINDS <= 0.01)  THEN
      ErrStat = ErrID_Fatal 
      ErrMess = ' Error (in WInDS): The timestep duration is recommanded to be 0.01-0.2 sec. Please check DT_RATIO option.'
      RETURN   
   END IF    
   
   p%FVM%NT         =   P%FVM%Total_Time / p%FVM%DT_WINDS + 1  ! sliu: would be better if obtained from FAST. Notice "+1", plus one for 0 time.    
     
   !-----------------------------------------------------------------------------------------------------
    ! Constants associated with Ramasamy-Leishman vortex model
   p%FVM%RL_Model%ALPHA  =  1.25643_DbKi
   !p%FVM%RL_Model%NU     =  1e-5_DbKi
   p%FVM%RL_Model%NU     =  1.4118e-5_DbKi   
   p%FVM%RL_Model%DELTA  =  100_DbKi
   p%FVM%RL_Model%A1     =  6.5e-5_DbKi   
   
      ! Atmospheric properties   
   p%FVM%RHO      = 1.23_DbKi    ! or p%FVM%RHO = p%Wind%Rho  
      
      ! Gravity  
   p%FVM%Gravity  =  9.81_DbKi   
   
   
   !-----------------------------------------------------------------------------------------------------
     ! Current date and time
   CALL DATE_AND_TIME(DATE = DATE_WINDS, TIME = TIME_WINDS)
   p%FVM%CURRENT_TIME   =  DATE_WINDS//'_'//TIME_WINDS(1:6)   ! hhmmss - CCYYMMDD      Refer: http://docs.oracle.com/cd/E19957-01/805-4942/6j4m3r8t2/index.html
                                    ! The current time, in the form hhmmss.sss, where hh is the hour, mm minutes, and ss.sss seconds and milliseconds.
                                    ! Date, in form CCYYMMDD, where CCYY is the four-digit year, MM the two-digit month, and DD the two-digit day of the month.
   
   !-----------------------------------------------------------------------------------------------------
      ! Parameters for blade element
   IF (.NOT. ALLOCATED( p%BLADE%RTrail ) )      ALLOCATE( p%BLADE%RTrail(NST)) 
   IF (.NOT. ALLOCATED( p%BLADE%AeroCen ) )     ALLOCATE( p%BLADE%AeroCen(NS))    
   IF (.NOT. ALLOCATED( p%BLADE%RNodes ) )      ALLOCATE( p%BLADE%RNodes(NS)) 

   p%BLADE%RTrail     =  0.0_DbKi
   p%BLADE%AeroCen    =  0.25_DbKi   ! sliu: thin airfoil assumption
   p%BLADE%RNodes     =  0.0_DbKi

   
   ! In the Matlab WInDS, these data is stored in NRELrotor.m
   p%BLADE%RNodes(1) =  p%HubRad + p%Blade%DR(1)/2
   DO IElement = 2, p%FVM%NS
      p%BLADE%RNodes(IElement) = p%BLADE%RNodes(IElement-1) +(p%Blade%DR(IElement-1) + p%Blade%DR(IElement))/2
   END DO

   p%BLADE%RTrail(1) = p%HubRad
   DO IElement2 = 2, p%FVM%NST
      p%BLADE%RTrail (IElement2) = p%BLADE%RTrail(IElement2 -1) + p%Blade%DR(IElement2 -1)
   END DO
       

   
   !-----------------------------------------------------------------------------------------------------
   ! Convert Cutoff parameter from [m] to number of timesteps and process information for initialization Timestep
   ! Full wake (the original model)
   IF (.NOT. p%FVM%WakeFLAG) THEN          
      p%FVM%WakeDist = -1
      p%FVM%UindPast = .FALSE. 
   END IF
   
       
   p%FVM%WakeNum = ceiling( p%FVM%WakeDist/ p%FVM%DT_WINDS / p%FVM%AveSpeed ) ! Left p%FVM%WakeNum is in timestep and right p%FVM%WakeDist is in meter
   
   IF (p%FVM%WakeNum <= 5_IntKi .OR. p%FVM%WakeNum >= p%FVM%NT ) THEN
      p%FVM%NTW_total = p%FVM%NT             ! Store all wake nodes if cutoff is disabled or set higher than total node count  
      O%FVM_Other%NTW = 1_IntKi  ! Initialize user.ntw with 1 for next functions
   ELSE
      p%FVM%NTW_total = p%FVM%WakeNum   ! Number of timesteps to be stored
      O%FVM_Other%NTW = 1_IntKi  ! Initialize user.ntw with 1 for next functions 
   END IF    
  
   
   
    ! Convert and process the roll parameter for the initialisation step
   IF (p%FVM%WakeFLAG) THEN
      O%FVM_Other%ntroll = ceiling( p%FVM%RollDist / p%FVM%DT_WINDS / p%FVM%AveSpeed ) ! O%FVM_Other%ntroll is in timestep and p%FVM%RollDist is in meter
   ELSE     
      IF (p%FVM%Roll) THEN
         O%FVM_Other%ntroll = p%FVM%NTW_total
      ELSE
         O%FVM_Other%ntroll = 0_IntKi       
      END IF
   END IF
  
   
   
END SUBROUTINE WINDS_SetParameters

!==================================================================================================================================
SUBROUTINE WINDS_Allocate( p, O, xd, ErrStat, ErrMess )
! This subroutine sets the parameters, based on the data stored in InputFileData
! Called from: AD_Init (in AeroDyn.f90)
! (~ initialize_pos_vel_perf.m / initials.m / some other functions)
!....................................................................

   IMPLICIT                        NONE

      ! Passed variables

   TYPE(AD_ParameterType),      INTENT(IN   )   :: p              ! The module's parameter data
   TYPE(AD_OtherStateType),     INTENT(INOUT)   :: O              ! Other/optimization states   
   TYPE(AD_DiscreteStateType),  INTENT(IN   )   :: xd          ! Discrete states at t   
   INTEGER(IntKi),              INTENT(  OUT)   :: ErrStat        ! The error status code
   CHARACTER(*),                INTENT(  OUT)   :: ErrMess        ! The error message, if an error occurred


      ! Local variables
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS  
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB 
   INTEGER(IntKi)      :: NTP
   INTEGER(IntKi)      :: NTW
      
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT 
   NT  = NT + 3    ! To make sure when using RK4, the array is big enough ...  

   NB  = p%NumBl
   NTP = p%FVM%NTP
   NTW = p%FVM%NTW_total

      !-----------------------------------------------------------------
      ! position variables  (initialize_pos_vel_perf.m) 
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_AEROCENT ) ) THEN
      ALLOCATE( O%FVM_Other%POS_AEROCENT(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_AEROCENT.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%POS_AEROCENT         =   0.0_DbKi         
      
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_LEAD ) ) THEN
      ALLOCATE( O%FVM_Other%POS_LEAD(3, NT, 1, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_LEAD.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF
   END IF   
   O%FVM_Other%POS_LEAD          =   0.0_DbKi   ! not used, but exists in Matlab code
      
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_BOUND ) ) THEN
      ALLOCATE( O%FVM_Other%POS_BOUND(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_BOUND.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%POS_BOUND         =   0.0_DbKi
   
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_COLLOC ) ) THEN
      ALLOCATE( O%FVM_Other%POS_COLLOC(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_COLLOC.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF   
   O%FVM_Other%POS_COLLOC        =   0.0_DbKi   ! not used, but exists in Matlab code
   
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_QUARTER ) ) THEN
      ALLOCATE( O%FVM_Other%POS_QUARTER(3, NT, 1, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_QUARTER.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF   
   O%FVM_Other%POS_QUARTER       =   0.0_DbKi
      
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_TRAIL ) ) THEN
      ALLOCATE( O%FVM_Other%POS_TRAIL(3, NT, 1, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_TRAIL .'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%POS_TRAIL         =   0.0_DbKi
              
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_END ) ) THEN
      ALLOCATE( O%FVM_Other%POS_END(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_END.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%POS_END           =   0.0_DbKi 
   
    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_NODES_BXN ) ) THEN
      ALLOCATE( O%FVM_Other%POS_NODES_BXN(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_NODES_BXN.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%POS_NODES_BXN     =   0.0_DbKi 
   
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_NODES_BYN ) ) THEN
      ALLOCATE( O%FVM_Other%POS_NODES_BYN(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_NODES_BYN.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%POS_NODES_BYN     =   0.0_DbKi 
   
   IF (.NOT. ALLOCATED( O%FVM_Other%POS_NODES_BZN ) ) THEN
      ALLOCATE( O%FVM_Other%POS_NODES_BZN(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%POS_NODES_BZN.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF           
   O%FVM_Other%POS_NODES_BZN     =   0.0_DbKi
   

     !------------------------------------------------------------------
     ! velocity variables  ( ~initialize_pos_vel_perf)
     !....................................
  IF (.NOT. ALLOCATED( O%FVM_Other%VEL_AEROCENT ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_AEROCENT(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_AEROCENT.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%VEL_AEROCENT         =   0.0_DbKi      
      
      
   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_BOUND ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_BOUND(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_BOUND.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF   
   O%FVM_Other%vel_BOUND         =   0.0_DbKi       
      
   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_BLADE ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_BLADE(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_BLADE.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF         
   O%FVM_Other%vel_BLADE         =   0.0_DbKi
      
   
   
     !------------------------------------------------------------------
     ! Performance variables
     !....................................     
   IF (.NOT. ALLOCATED( O%FVM_Other%PERF_CL ) ) THEN
      ALLOCATE( O%FVM_Other%PERF_CL(1, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%PERF_CL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF             
   O%FVM_Other%PERF_CL           =   0.0_DbKi
    
   IF (.NOT. ALLOCATED( O%FVM_Other%PERF_CD ) ) THEN
      ALLOCATE( O%FVM_Other%PERF_CD(1, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%PERF_CD.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%PERF_CD           =   0.0_DbKi
                    
   IF (.NOT. ALLOCATED( O%FVM_Other%PERF_AOA ) ) THEN
      ALLOCATE( O%FVM_Other%PERF_AOA(1, NT, 1, NS, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%PERF_AOA.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF   
   O%FVM_Other%PERF_AOA          =   0.0_DbKi

   
   IF (.NOT. ALLOCATED( O%FVM_Other%PERF_A ) ) THEN
      ALLOCATE( O%FVM_Other%PERF_A(1, NT, 1, NS, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%PERF_A.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF        
   O%FVM_Other%PERF_A            =   0.0_DbKi

   IF (.NOT. ALLOCATED( O%FVM_Other%PERF_CM ) ) THEN
      ALLOCATE( O%FVM_Other%PERF_CM(1, NT, 1, NS, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%PERF_CM.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF        
   O%FVM_Other%PERF_CM           =   0.0_DbKi 


     !------------------------------------------------------------------
     ! Wake and velocity ( ~initials.m)
     !....................................    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_DOMAIN ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_DOMAIN(3, NTW + 1, NTP, NST, NB ), STAT = ErrStat )  
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_DOMAIN.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF        
   O%FVM_Other%WAKE_DOMAIN           =   0.0_DbKi 
 

   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_DOMAIN ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_DOMAIN(3, NTW + 1, NTP, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_DOMAIN.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_DOMAIN            =   0.0_DbKi 
          
 
   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_DOMAIN_RK ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_DOMAIN_RK(3, NTW + 1, 3, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_DOMAIN_RK.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_DOMAIN_RK         =   0.0_DbKi    
   

    IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UIND ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_UIND(3, NTW + 1, NTP, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UIND.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF        
   O%FVM_Other%VEL_UIND              =   0.0_DbKi

   
              
   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UINDB ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_UINDB(3, NT +1, NTP, NS, NB ), STAT = ErrStat )    
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UINDB.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_UINDB                     =   0.0_DbKi
   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UIND_RK ) ) THEN
      ALLOCATE( O%FVM_Other%VEL_UIND_RK(3, NTW + 1, 3, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UIND_RK.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF        
   O%FVM_Other%VEL_UIND_RK              =   0.0_DbKi
   

   ! Ground effects....
   IF (p%FVM%Ground_Parms%GroundFLAG) THEN
       
       IF (p%FVM%Ground_Parms%Method  == 'PANEL' ) THEN
           
           IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UIND_GROUND ) ) THEN
              ALLOCATE( O%FVM_Other%VEL_UIND_GROUND(3, NTW + 1, NTP, NST, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%VEL_UIND_GROUND.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%VEL_UIND_GROUND       =   0.0_DbKi
         
           
           IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UINDB_GROUND ) ) THEN
              ALLOCATE( O%FVM_Other%VEL_UINDB_GROUND(3, NTW + 1, NTP, NS, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%VEL_UINDB_GROUND.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%VEL_UINDB_GROUND       =   0.0_DbKi
       
       ELSE IF (p%FVM%Ground_Parms%Method  == 'IMAGE' ) THEN

           IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UIND_MIRROR ) ) THEN
              ALLOCATE( O%FVM_Other%VEL_UIND_MIRROR(3, NTW + 1, NTP, NST, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%VEL_UIND_MIRROR.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%VEL_UIND_MIRROR       =   0.0_DbKi
  
           
           IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UINDB_MIRROR ) ) THEN
              ALLOCATE( O%FVM_Other%VEL_UINDB_MIRROR(3, NTW + 1, NTP, NS, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%VEL_UINDB_MIRROR.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%VEL_UINDB_MIRROR       =   0.0_DbKi
           
           IF (.NOT. ALLOCATED( O%FVM_Other%GROUND%Vind_mirtmp_SHED ) ) THEN
              ALLOCATE( O%FVM_Other%GROUND%Vind_mirtmp_SHED(3, NTW + 1, 1, NST, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%GROUND%Vind_mirtmp_SHED.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%GROUND%Vind_mirtmp_SHED  =   0.0_DbKi           
           
           IF (.NOT. ALLOCATED( O%FVM_Other%GROUND%Vind_mirtmp_TRAIL ) ) THEN
              ALLOCATE( O%FVM_Other%GROUND%Vind_mirtmp_TRAIL(3, NTW + 1, 1, NST, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%GROUND%Vind_mirtmp_TRAIL.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%GROUND%Vind_mirtmp_TRAIL  =   0.0_DbKi             
           
           IF (.NOT. ALLOCATED( O%FVM_Other%GROUND%mirror ) ) THEN
              ALLOCATE( O%FVM_Other%GROUND%mirror(3, NTW, 1, NST, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%GROUND%mirror.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%GROUND%mirror  =   0.0_DbKi                   
           
           IF (.NOT. ALLOCATED( O%FVM_Other%GROUND%wake_mirror ) ) THEN
              ALLOCATE( O%FVM_Other%GROUND%wake_mirror(3, NTW, 1, NST, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%GROUND%wake_mirror.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%GROUND%wake_mirror  =   0.0_DbKi            
           
            IF (.NOT. ALLOCATED( O%FVM_Other%GROUND%gamma_shed_mirror ) ) THEN
              ALLOCATE( O%FVM_Other%GROUND%gamma_shed_mirror(1, NTW, 1, NS, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%GROUND%gamma_shed_mirror.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%GROUND%gamma_shed_mirror  =   0.0_DbKi      
           
           IF (.NOT. ALLOCATED( O%FVM_Other%GROUND%gamma_trail_mirror ) ) THEN
              ALLOCATE( O%FVM_Other%GROUND%gamma_trail_mirror(1, NTW -1, 1, NST, NB ), STAT = ErrStat )   
              IF ( ErrStat /= 0 ) THEN
                 ErrMess = ' Error allocating space for O%FVM_Other%GROUND%gamma_trail_mirror.'
                 ErrStat = ErrID_Fatal
                 RETURN         
              END IF   
           END IF        
           O%FVM_Other%GROUND%gamma_trail_mirror  =   0.0_DbKi  
           
   
       END IF !  (p%FVM%Ground_Parms%Method == 'PANEL' )   
       
   END IF !  (p%FVM%Ground_Parms%GroundFLAG)   
   
   
      
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_RE_SHED ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_RE_SHED(1, NTW + 1, NTP, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_RE_SHED.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%WAKE_RE_SHED          =   0.0_DbKi 
      
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_RE_TRAIL ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_RE_TRAIL(1, NTW + 1, NTP, NST, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_RE_TRAIL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%WAKE_RE_TRAIL         =   0.0_DbKi

   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_RC_SHED ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_RC_SHED(1, NTW + 1, NTP, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_RC_SHED.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF     
   O%FVM_Other%WAKE_RC_SHED          =   0.0_DbKi 
   
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_RC_TRAIL ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_RC_TRAIL(1, NTW+1, NTP, NST, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_RC_TRAIL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%WAKE_RC_TRAIL         =   0.0_DbKi
   

   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_LENGTH_SHED ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_LENGTH_SHED(1, NTW+1, NTP, NS, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_LENGTH_SHED.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%WAKE_LENGTH_SHED     =   0.0_DbKi
   
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_LENGTH_TRAIL ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_LENGTH_TRAIL(1, NTW+1, NTP, NST, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_LENGTH_TRAIL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%WAKE_LENGTH_TRAIL    =   0.0_DbKi
   
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_RC_EFF_SHED ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_RC_EFF_SHED(1, NTW+1, NTP, NS, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_RC_EFF_SHED.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF    
   O%FVM_Other%WAKE_RC_EFF_SHED    =   0.0_DbKi 
     
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_RC_EFF_TRAIL ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_RC_EFF_TRAIL(1, NTW+1, NTP, NST, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_RC_EFF_TRAIL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%WAKE_RC_EFF_TRAIL    =   0.0_DbKi
      
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_GAMMA_SHED ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_GAMMA_SHED(1, NTW+1, NTP, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_GAMMA_SHED.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF      
   O%FVM_Other%WAKE_GAMMA_SHED      =   0.0_DbKi
      
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_GAMMA_TRAIL ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_GAMMA_TRAIL(1, NTW+1, NTP, NST, NB), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_GAMMA_TRAIL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%WAKE_GAMMA_TRAIL     =   0.0_DbKi
    
   
   
   
     !------------------------------------------------------------------
     ! initials.m 
     !....................................
   IF (.NOT. ALLOCATED( O%FVM_Other%WIND_INFTY ) ) THEN
      ALLOCATE( O%FVM_Other%WIND_INFTY(3, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WIND_INFTY.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF         
   O%FVM_Other%WIND_INFTY      =   0.0_DbKi

   IF (.NOT. ALLOCATED( O%FVM_Other%WIND_INFTYM ) ) THEN
      ALLOCATE( O%FVM_Other%WIND_INFTYM(1, NT, 1, NS, NB ), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WIND_INFTYM.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF         
   O%FVM_Other%WIND_INFTYM     =   0.0_DbKi 


     !------------------------------------------------------------------
     !  Initials.m  Line 240
     !....................................
   IF (.NOT. ALLOCATED( O%FVM_Other%WAKE_R0 ) ) THEN
      ALLOCATE( O%FVM_Other%WAKE_R0(NT), STAT = ErrStat )   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%WAKE_R0.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF     
   O%FVM_Other%WAKE_R0        =   0.0_DbKi 
  
 
   IF ( .NOT. ALLOCATED( o%FVM_Other%StoredForces ))  THEN
      ALLOCATE(o%FVM_Other%StoredForces(3, NT, 1, NS, NB), STAT = ErrStat)
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%StoredForces.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%StoredForces    =   0.0_DbKi
   
     IF ( .NOT. ALLOCATED( o%FVM_Other%StoredMoments ))  THEN
      ALLOCATE(o%FVM_Other%StoredMoments(3, NT, 1, NS, NB), STAT = ErrStat)
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for o%FVM_Other%StoredMoments.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%StoredMoments   =   0.0_DbKi 
   
   
   IF ( .NOT. ALLOCATED( o%FVM_Other%PreviousForces ))  THEN
      ALLOCATE(o%FVM_Other%PreviousForces(3, 1, 1, NS, NB), STAT = ErrStat)
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%PreviousForces.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%PreviousForces    =   0.0_DbKi
   
     IF ( .NOT. ALLOCATED( o%FVM_Other%PreviousMoments ))  THEN
      ALLOCATE(o%FVM_Other%PreviousMoments(3, 1, 1, NS, NB), STAT = ErrStat)
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for o%FVM_Other%PreviousMoments.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%PreviousMoments   =   0.0_DbKi    
   
   IF ( .NOT. ALLOCATED( O%FVM_Other%VEL_UIND_SHED ))  THEN
      ALLOCATE(O%FVM_Other%VEL_UIND_SHED(3, NT, 1, NST, NB), STAT = ErrStat)
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UIND_SHED.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_UIND_SHED   =   0.0_DbKi    

   IF ( .NOT. ALLOCATED( O%FVM_Other%VEL_UIND_TRAIL ))  THEN
      ALLOCATE(O%FVM_Other%VEL_UIND_TRAIL(3, NT, 1, NST, NB), STAT = ErrStat)
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UIND_TRAIL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_UIND_TRAIL   =   0.0_DbKi    

   
   IF ( .NOT. ALLOCATED( O%FVM_Other%VEL_UINDB_SHED ))  THEN
      ALLOCATE(O%FVM_Other%VEL_UINDB_SHED(3, 1, 1, NS, NB), STAT = ErrStat)   
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UINDB_SHED.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_UINDB_SHED   =   0.0_DbKi  
   
   IF ( .NOT. ALLOCATED( O%FVM_Other%VEL_UINDB_TRAIL ))  THEN
      ALLOCATE(O%FVM_Other%VEL_UINDB_TRAIL(3, 1, 1, NS, NB), STAT = ErrStat)    
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UINDB_TRAIL.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_UINDB_TRAIL   =   0.0_DbKi  

   IF ( .NOT. ALLOCATED( O%FVM_Other%VEL_UINDB_SHED_pre ))  THEN
      ALLOCATE(O%FVM_Other%VEL_UINDB_SHED_pre(3, 1, 1, NS, NB), STAT = ErrStat)  
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UINDB_SHED_pre.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_UINDB_SHED_pre   =   0.0_DbKi  
   
   IF ( .NOT. ALLOCATED( O%FVM_Other%VEL_UINDB_TRAIL_pre ))  THEN
      ALLOCATE(O%FVM_Other%VEL_UINDB_TRAIL_pre(3, 1, 1, NS, NB), STAT = ErrStat) 
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%VEL_UINDB_TRAIL_pre.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF   
   END IF 
   O%FVM_Other%VEL_UINDB_TRAIL_pre   =   0.0_DbKi  
      
   
     !------------------------------------------------------------------
     ! These are used for KJ interations:
     !....................................
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%DG ) )           ALLOCATE( O%FVM_Other%KJ%DG(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%GAMMA1 ) )       ALLOCATE( O%FVM_Other%KJ%GAMMA1(1, 1, 1, NS,  NB))      

   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%CM ) )           ALLOCATE( O%FVM_Other%KJ%CM(1, 1, 1, NS,  NB))  
   
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%CL ) )           ALLOCATE( O%FVM_Other%KJ%CL(1, 1, 1, NS,  NB))   
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%CD ) )           ALLOCATE( O%FVM_Other%KJ%CD(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%VEL_ROT ) )      ALLOCATE( O%FVM_Other%KJ%VEL_ROT(3, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%VEL_TOT ) )      ALLOCATE( O%FVM_Other%KJ%VEL_TOT(3, NT, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%U ) )            ALLOCATE( O%FVM_Other%KJ%U(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%V) )             ALLOCATE( O%FVM_Other%KJ%V(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%W ) )            ALLOCATE( O%FVM_Other%KJ%W(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%Vinf ) )         ALLOCATE( O%FVM_Other%KJ%Vinf(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%Vtot ) )         ALLOCATE( O%FVM_Other%KJ%Vtot(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%AOA ) )          ALLOCATE( O%FVM_Other%KJ%AOA(1, 1, 1, NS,  NB)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%KJ%CM ) )           ALLOCATE( O%FVM_Other%KJ%CM(1, 1, 1, NS,  NB))     
   
   
   ! Space for LB dynamic stall
   IF (p%FVM%DS_Parms%DS_Flag) THEN
      CALL LB_Initialize_Variables(p, O, xd, ErrStat, ErrMess)
      IF ( ErrStat /= 0 ) THEN
         ErrMess = ' Error allocating space for O%FVM_Other%DS.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF 
   END IF
   
  
   
   
END SUBROUTINE WINDS_Allocate
!==================================================================================================================================
SUBROUTINE WINDS_check_cutoff(u, p, O, xd, N, ErrStat, ErrMess)
! Called from: AD_CalcOutput  (in AeroDyn.f90)
! (~ check_cutoff.m)
!
! Original written by David Lenz!
!....................................................................
  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_InputType),            INTENT(IN   )  :: u           ! Inputs at t
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep, counts from 1    
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess

   
   ! Process the Cutoff parameter   
   IF (p%FVM%WakeNum <= N .AND. p%FVM%WakeNum > 5_IntKi ) THEN
      O%FVM_Other%NTW = p%FVM%WakeNum
   ELSE
      O%FVM_Other%NTW = N
   END IF    
         
   ! Process the roll parameter to cut induction
   IF (p%FVM%WakeFLAG) THEN
      IF (O%FVM_Other%NTW < p%FVM%RollDist) THEN
         O%FVM_Other%ntroll = O%FVM_Other%NTW 
      ELSE  
         O%FVM_Other%ntroll = p%FVM%RollDist
      END IF      
   ELSE     
      IF (p%FVM%roll) THEN 
         O%FVM_Other%ntroll = O%FVM_Other%NTW 
      ELSE 
         O%FVM_Other%ntroll = 0_IntKi
      END IF      
   END IF   
   
   
   
END SUBROUTINE WINDS_check_cutoff

!==================================================================================================================================
SUBROUTINE WINDS_Ground_Model(p, O, ErrStat, ErrMess)
! ~ ground_model.m
!....................................................................
  IMPLICIT                        NONE

      ! Passed variables
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState 
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess   
   
   
   INTEGER(IntKi)      :: i
   INTEGER(IntKi)      :: j
   INTEGER(IntKi)      :: k
   INTEGER(IntKi)      :: NP  ! sqrt of panels # 

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS  
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB 
   INTEGER(IntKi)      :: NTP
      
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   NTP = p%FVM%NTP        
   
   
   !O%FVM_Other%Ground%Np     = p%FVM%Ground_Parms%Sqrt_Panels  ! square root of number of panels
   NP = p%FVM%Ground_Parms%Sqrt_Panels

 ! As different with the Matlab WInDS, the ground geometry is not saved.

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%xp ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%xp(NP + 1 , 1), STAT = ErrStat )   
   END IF   

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%yp ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%yp(NP + 1 , 1), STAT = ErrStat )   
   END IF   

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%zp ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%zp(NP + 1 , 1), STAT = ErrStat )   
   END IF  

   ! Mesh
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%xp_mesh ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%xp_mesh(NP+1, NP+1, NP+1), STAT = ErrStat )   
   END IF     
   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%yp_mesh ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%yp_mesh(NP+1, NP+1, NP+1), STAT = ErrStat )   
   END IF        
   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%zp_mesh ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%zp_mesh(NP+1, NP+1, NP+1), STAT = ErrStat )   
   END IF        
   
     ! panel strengths
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%Gamma ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%Gamma(NP **2 , p%FVM%NT), STAT = ErrStat )   
   END IF     
   O%FVM_Other%Ground%Gamma  =  0.0_DbKi

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%p ) ) THEN
      ALLOCATE( O%FVM_Other%Ground%p(NP **2 , 3), STAT = ErrStat )   
   END IF     
   
   ! For inducedvelocity subroutine
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%P_source ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%P_source(3, 1, 1, NP **2, 1), STAT = ErrStat )   
   END IF        
   
   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UIND_SHED_GROUND_P ) ) THEN             
      ALLOCATE( O%FVM_Other%VEL_UIND_SHED_GROUND_P(3, 1, 1, NP **2, 1), STAT = ErrStat )   
   END IF        
      
   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UIND_TRAIL_GROUND_P ) ) THEN             
      ALLOCATE( O%FVM_Other%VEL_UIND_TRAIL_GROUND_P(3, 1, 1, NP **2, 1), STAT = ErrStat )   
   END IF         

   IF (.NOT. ALLOCATED( O%FVM_Other%VEL_UIND_GROUND_P ) ) THEN             
      ALLOCATE( O%FVM_Other%VEL_UIND_GROUND_P(3, 1, 1, NP **2, 1), STAT = ErrStat )   
   END IF   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_n ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_n(NP **2, 1), STAT = ErrStat )   
   END IF     
   
   ! For Induced_Velocity_Ground_Panels subroutine   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%gamma_grid ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%gamma_grid(1, NP, 1, NP, 1), STAT = ErrStat )   
   END IF    

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%rc_grid ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%rc_grid(1, NP, 1, NP, 1), STAT = ErrStat )   
   END IF  
   
   !...........      
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%grid_lefttop ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%grid_lefttop(3, NP, 1, NP, 1), STAT = ErrStat )   
   END IF       

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%grid_righttop ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%grid_righttop(3, NP, 1, NP, 1), STAT = ErrStat )   
   END IF    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%grid_rightbottom ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%grid_rightbottom(3, NP, 1, NP, 1), STAT = ErrStat )   
   END IF    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%grid_leftbottom ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%grid_leftbottom(3, NP, 1, NP, 1), STAT = ErrStat )   
   END IF       
   
   !.........
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_top_main ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_top_main(3, NT, 1, NST, NB), STAT = ErrStat )   
   END IF         
   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_right_main ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_right_main(3, NT, 1, NST, NB), STAT = ErrStat )   
   END IF            

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_bottom_main ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_bottom_main(3, NT, 1, NST, NB), STAT = ErrStat )   
   END IF

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_left_main ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_left_main(3, NT, 1, NST, NB), STAT = ErrStat )   
   END IF
   
   !.........
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_top_KJ ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_top_KJ(3, 1, 1, NST, NB), STAT = ErrStat )   
   END IF         
   
   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_right_KJ ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_right_KJ(3, 1, 1, NST, NB), STAT = ErrStat )   
   END IF            

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_bottom_KJ ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_bottom_KJ(3, 1, 1, NST, NB), STAT = ErrStat )   
   END IF

   IF (.NOT. ALLOCATED( O%FVM_Other%Ground%V_ind_left_KJ ) ) THEN             
      ALLOCATE( O%FVM_Other%Ground%V_ind_left_KJ(3, 1, 1, NST, NB), STAT = ErrStat )   
   END IF   
   
   
      
    ! Define geometry of panels. note this would happen in an initialization and not need to be called at each time step 
   O%FVM_Other%Ground%x_min  = p%FVM%Ground_Parms%Extent(1)
   O%FVM_Other%Ground%x_max  = p%FVM%Ground_Parms%Extent(2)
   O%FVM_Other%Ground%dx     = (O%FVM_Other%Ground%x_max -  O%FVM_Other%Ground%x_min) / NP    


   DO i = 1, NP+1
      O%FVM_Other%Ground%xp(i , 1) = O%FVM_Other%Ground%x_min + O%FVM_Other%Ground%dx * (i-1)
   END DO   

   O%FVM_Other%Ground%y_min  = p%FVM%Ground_Parms%Extent(3)
   O%FVM_Other%Ground%y_max  = p%FVM%Ground_Parms%Extent(4)
   O%FVM_Other%Ground%dy     = (O%FVM_Other%Ground%y_max -  O%FVM_Other%Ground%y_min) / NP   


   DO i = 1, NP+1
      O%FVM_Other%Ground%yp(i , 1) = O%FVM_Other%Ground%y_min + O%FVM_Other%Ground%dy * (i-1)
   END DO   
   
   O%FVM_Other%Ground%zp = 0.0_DbKi
   
    ! define the nodes of the vortex panel grid  
   DO i = 1, NP +1
      DO j = 1, NP +1
         DO k = 1, NP +1
            O%FVM_Other%Ground%xp_mesh(i, j ,k)  = O%FVM_Other%Ground%xp(j , 1)
            O%FVM_Other%Ground%yp_mesh(i, j ,k)  = O%FVM_Other%Ground%yp(i , 1) 
            O%FVM_Other%Ground%zp_mesh(i, j ,k)  = O%FVM_Other%Ground%zp(k , 1)             
         END DO
      END DO
   END DO
   
   
     ! cut off distance fraction for vortex segments 
   O%FVM_Other%Ground%co     =  0.05_DbKi
   

   
   i = 1_IntKi
   DO j = 1, NP
      DO k = 1, NP
         O%FVM_Other%Ground%p(i,1)  =  O%FVM_Other%Ground%xp_mesh(j,k,1) + 0.5 * (O%FVM_Other%Ground%xp_mesh(j,k+1,1) - &
                                                                                    O%FVM_Other%Ground%xp_mesh(j,k,1))
                                                                                    
         O%FVM_Other%Ground%p(i,2)  =  O%FVM_Other%Ground%yp_mesh(j,k,1) + 0.5 * (O%FVM_Other%Ground%yp_mesh(j+1,1,1) - &
                                                                                    O%FVM_Other%Ground%yp_mesh(j,1,1))    
                                                                                      
         O%FVM_Other%Ground%p(i,3)  =  O%FVM_Other%Ground%zp_mesh(1,1,k) + 0.5 * (O%FVM_Other%Ground%xp_mesh(1,1,k+1) - &
                                                                                    O%FVM_Other%Ground%xp_mesh(1,1,k))      
         i = i + 1                                                                           
                                                                                    
      END DO
   END DO    


   CALL ground_influence_coefficients(p, O, ErrStat, ErrMess)
   IF (ErrStat /= 0 ) RETURN
   

   DO j = 1, NP * NP
      O%FVM_Other%Ground%p_source(1, 1, 1, j, 1)  =  O%FVM_Other%Ground%p(j, 1)                               
      O%FVM_Other%Ground%p_source(2, 1, 1, j, 1)  =  O%FVM_Other%Ground%p(j, 2)                                 
      O%FVM_Other%Ground%p_source(3, 1, 1, j, 1)  =  O%FVM_Other%Ground%p(j, 3)                                                                         
   END DO        
   
   
   

CONTAINS
   !...............................................................................................................................
   SUBROUTINE ground_influence_coefficients(p, O, ErrStat, ErrMess)
   
      IMPLICIT                        NONE

        ! Passed variables
      TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState 
      INTEGER,                       INTENT(INOUT)  :: ErrStat
      CHARACTER(*),                  INTENT(INOUT)  :: ErrMess   
   
      REAL(DbKi)          :: summation
      
      INTEGER(IntKi)      :: i       
      INTEGER(IntKi)      :: j
      INTEGER(IntKi)      :: m
      INTEGER(IntKi)      :: n
      INTEGER(IntKi)      :: k
 
      
      INTEGER(IntKi)      :: Np     ! Square root of number of panels
      INTEGER(IntKi)      :: NN     ! The amount of panels = Np * Np
      REAL(DbKi)          :: dx
      REAL(DbKi)          :: dy
      

      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: r1
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: r2      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: r3
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: r4
      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: n1
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: n2      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: n3
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: n4      
      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: t1
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: t2      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: t3
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: t4
      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: c1
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: c2      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: c3
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: c4  
      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: a1
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: a2      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: a3
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: a4   
      
!       REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: n_vec
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: n_vec_transpose
      
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: A 

      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: x 
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: y 
      REAL(DbKi) , DIMENSION(:,:), ALLOCATABLE  :: z 


      Np  =  p%FVM%Ground_Parms%Sqrt_Panels   ! or 
      NN  =  Np ** 2
      dx  =  O%FVM_Other%Ground%dx
      dy  =  O%FVM_Other%Ground%dy      
      
      
      IF (.NOT. ALLOCATED( r1 ) ) THEN
         ALLOCATE( r1(NN , 3), STAT = ErrStat )   
      END IF              
      IF (.NOT. ALLOCATED( r2 ) ) THEN
         ALLOCATE( r2(NN , 3), STAT = ErrStat )   
      END IF          
      IF (.NOT. ALLOCATED( r3 ) ) THEN
         ALLOCATE( r3(NN , 3), STAT = ErrStat )   
      END IF                      
      IF (.NOT. ALLOCATED( r4 ) ) THEN
         ALLOCATE( r4(NN , 3), STAT = ErrStat )   
      END IF          
      r1 = 0.0_DbKi
      r2 = 0.0_DbKi
      r3 = 0.0_DbKi
      r4 = 0.0_DbKi           
      
      IF (.NOT. ALLOCATED( n1 ) ) THEN
         ALLOCATE( n1(NN , 1), STAT = ErrStat )   
      END IF              
      IF (.NOT. ALLOCATED( n2 ) ) THEN
         ALLOCATE( n2(NN , 1), STAT = ErrStat )   
      END IF          
      IF (.NOT. ALLOCATED( n3 ) ) THEN
         ALLOCATE( n3(NN , 1), STAT = ErrStat )   
      END IF                      
      IF (.NOT. ALLOCATED( n4 ) ) THEN
         ALLOCATE( n4(NN , 1), STAT = ErrStat )   
      END IF        
      n1 = 0.0_DbKi
      n2 = 0.0_DbKi
      n3 = 0.0_DbKi
      n4 = 0.0_DbKi          
      
      IF (.NOT. ALLOCATED( t1 ) ) THEN
         ALLOCATE( t1(NN , 1), STAT = ErrStat )   
      END IF              
      IF (.NOT. ALLOCATED( t2 ) ) THEN
         ALLOCATE( t2(NN , 1), STAT = ErrStat )   
      END IF          
      IF (.NOT. ALLOCATED( t3 ) ) THEN
         ALLOCATE( t3(NN , 1), STAT = ErrStat )   
      END IF                      
      IF (.NOT. ALLOCATED( t4 ) ) THEN
         ALLOCATE( t4(NN , 1), STAT = ErrStat )   
      END IF          
      t1 = 0.0_DbKi
      t2 = 0.0_DbKi
      t3 = 0.0_DbKi
      t4 = 0.0_DbKi         
      
      IF (.NOT. ALLOCATED( c1 ) ) THEN
         ALLOCATE( c1(NN , 3), STAT = ErrStat )   
      END IF              
      IF (.NOT. ALLOCATED( c2 ) ) THEN
         ALLOCATE( c2(NN , 3), STAT = ErrStat )   
      END IF          
      IF (.NOT. ALLOCATED( c3 ) ) THEN
         ALLOCATE( c3(NN , 3), STAT = ErrStat )   
      END IF                      
      IF (.NOT. ALLOCATED( c4 ) ) THEN
         ALLOCATE( c4(NN , 3), STAT = ErrStat )   
      END IF          
      c1 = 0.0_DbKi
      c2 = 0.0_DbKi
      c3 = 0.0_DbKi
      c4 = 0.0_DbKi           
 
      
      IF (.NOT. ALLOCATED( a1 ) ) THEN
         ALLOCATE( a1(NN , 3), STAT = ErrStat )   
      END IF              
      IF (.NOT. ALLOCATED( a2 ) ) THEN
         ALLOCATE( a2(NN , 3), STAT = ErrStat )   
      END IF          
      IF (.NOT. ALLOCATED( a3 ) ) THEN
         ALLOCATE( a3(NN , 3), STAT = ErrStat )   
      END IF                      
      IF (.NOT. ALLOCATED( a4 ) ) THEN
         ALLOCATE( a4(NN , 3), STAT = ErrStat )   
      END IF          
      a1 = 0.0_DbKi
      a2 = 0.0_DbKi
      a3 = 0.0_DbKi
      a4 = 0.0_DbKi     

      
      
      IF (.NOT. ALLOCATED( O%FVM_Other%Ground%n_vec ) ) THEN
         ALLOCATE( O%FVM_Other%Ground%n_vec(NN , 3), STAT = ErrStat )   
      END IF          
      O%FVM_Other%Ground%n_vec = 0.0_DbKi      
      
      IF (.NOT. ALLOCATED( n_vec_transpose ) ) THEN
         ALLOCATE( n_vec_transpose(3, 1), STAT = ErrStat )   
      END IF          
      n_vec_transpose = 0.0_DbKi      
      
      IF (.NOT. ALLOCATED( A ) ) THEN
         ALLOCATE( a(NN , NN), STAT = ErrStat )   
      END IF          
      a = 0.0_DbKi

            

      IF (.NOT. ALLOCATED( O%FVM_Other%Ground%b ) ) THEN
         ALLOCATE( O%FVM_Other%Ground%b(NN , NN), STAT = ErrStat )   
      END IF          
      O%FVM_Other%Ground%b = 0.0_DbKi      


      IF (.NOT. ALLOCATED( x ) ) THEN
         ALLOCATE( x(NN , 1), STAT = ErrStat )   
      END IF          
      IF (.NOT. ALLOCATED( y ) ) THEN
         ALLOCATE( y(NN , 1), STAT = ErrStat )   
      END IF   
      IF (.NOT. ALLOCATED( z ) ) THEN
         ALLOCATE( z(NN , 1), STAT = ErrStat )   
      END IF   

      
      x(1:NN , 1) = O%FVM_Other%Ground%p(1:NN , 1)
      y(1:NN , 1) = O%FVM_Other%Ground%p(1:NN , 2)
      z(1:NN , 1) = O%FVM_Other%Ground%p(1:NN , 3)
      
      O%FVM_Other%Ground%co = O%FVM_Other%Ground%co * (dx+dy)/2
      
      j = 1_IntKi
      
      DO n = 1, Np      ! count through y direction
         DO m = 1, Np   ! count through x direction
      
               r1(:,1) = x(:,1) - O%FVM_Other%Ground%xp_mesh(n,m,1)
               r1(:,2) = y(:,1) - O%FVM_Other%Ground%yp_mesh(n,m,1)
               r1(:,3) = z(:,1) - O%FVM_Other%Ground%zp_mesh(1,1,m)


               r2(:,1) = x(:,1) - O%FVM_Other%Ground%xp_mesh(n,m+1,1)
               r2(:,2) = y(:,1) - O%FVM_Other%Ground%yp_mesh(n,m+1,1)
               r2(:,3) = z(:,1) - O%FVM_Other%Ground%zp_mesh(1,1,m+1)
      
               r3(:,1) = x(:,1) - O%FVM_Other%Ground%xp_mesh(n+1,m+1,1)
               r3(:,2) = y(:,1) - O%FVM_Other%Ground%yp_mesh(n+1,m+1,1)
               r3(:,3) = z(:,1) - O%FVM_Other%Ground%zp_mesh(1,1,m+1)      
      
               r4(:,1) = x(:,1) - O%FVM_Other%Ground%xp_mesh(n+1,m,1)
               r4(:,2) = y(:,1) - O%FVM_Other%Ground%yp_mesh(n+1,m,1)
               r4(:,3) = z(:,1) - O%FVM_Other%Ground%zp_mesh(1,1,m)      
      
               n1(:,1) = (sum((r1 ** 2),2)) ** 0.5
               n2(:,1) = (sum((r2 ** 2),2)) ** 0.5
               n3(:,1) = (sum((r3 ** 2),2)) ** 0.5
               n4(:,1) = (sum((r4 ** 2),2)) ** 0.5     
      
            ! ------------------------------------
               DO i = 1, NN
                  summation = r1(i, 1) * r2(i, 1) + r1(i, 2) * r2(i, 2) + r1(i, 3) * r2(i, 3) 
                  t1(i,1) = (n1(i,1) + n2(i,1))/(n1(i,1) * n2(i,1) * (n1(i,1) * n2(i,1) + summation + O%FVM_Other%Ground%co ** 2))
               END DO    !t1(:,1) = (n1 + n2) / (n1 * n2 * (n1 * n2 + sum(r1 * r2 , 2) + O%FVM_Other%Ground%co ** 2))
               
               DO i = 1, NN
                  c1(i, 1) =  r1(i, 2) * r2(i, 3) - r1(i, 3) * r2(i, 2) 
                  c1(i, 2) =  r1(i, 3) * r2(i, 1) - r1(i, 1) * r2(i, 3) 
                  c1(i, 3) =  r1(i, 1) * r2(i, 2) - r1(i, 2) * r2(i, 1)     
               END DO  

               DO i = 1, NN
                  A1(i, 1) = c1(i, 1) * t1(i, 1)
                  A1(i, 2) = c1(i, 2) * t1(i, 1)
                  A1(i, 3) = c1(i, 3) * t1(i, 1)
               END DO 
               
            ! ------------------------------------   
               DO i = 1, NN
                  summation = r2(i, 1) * r3(i, 1) + r2(i, 2) * r3(i, 2) + r2(i, 3) * r3(i, 3) 
                  t2(i,1) = (n2(i,1) + n3(i,1))/(n2(i,1) * n3(i,1) * (n2(i,1) * n3(i,1) + summation + O%FVM_Other%Ground%co ** 2))
               END DO    !t2(:,1) = (n2 + n3) / (n2 * n3 * (n2 * n3 + sum(r2 * r3 , 2) + O%FVM_Other%Ground%co ** 2))
               
               DO i = 1, NN
                  c2(i, 1) =  r2(i, 2) * r3(i, 3) - r2(i, 3) * r3(i, 2) 
                  c2(i, 2) =  r2(i, 3) * r3(i, 1) - r2(i, 1) * r3(i, 3) 
                  c2(i, 3) =  r2(i, 1) * r3(i, 2) - r2(i, 2) * r3(i, 1)     
               END DO  

               DO i = 1, NN
                  A2(i, 1) = c2(i, 1) * t2(i, 1)
                  A2(i, 2) = c2(i, 2) * t2(i, 1)
                  A2(i, 3) = c2(i, 3) * t2(i, 1)
               END DO                
               
               
            ! ------------------------------------   
                DO i = 1, NN
                  summation = r3(i, 1) * r4(i, 1) + r3(i, 2) * r4(i, 2) + r3(i, 3) * r4(i, 3) 
                  t3(i,1) = (n3(i,1) + n4(i,1))/(n3(i,1) * n4(i,1) * (n3(i,1) * n4(i,1) + summation + O%FVM_Other%Ground%co ** 2))
                END DO    ! t3(:,1) = (n3 + n4) / (n3 * n4 * (n3 * n4 + sum(r3 * r4 , 2) + O%FVM_Other%Ground%co ** 2))
                
               DO i = 1, NN
                  c3(i, 1) =  r3(i, 2) * r4(i, 3) - r3(i, 3) * r4(i, 2) 
                  c3(i, 2) =  r3(i, 3) * r4(i, 1) - r3(i, 1) * r4(i, 3) 
                  c3(i, 3) =  r3(i, 1) * r4(i, 2) - r3(i, 2) * r4(i, 1)     
               END DO  

               DO i = 1, NN
                  A3(i, 1) = c3(i, 1) * t3(i, 1)
                  A3(i, 2) = c3(i, 2) * t3(i, 1)
                  A3(i, 3) = c3(i, 3) * t3(i, 1)
               END DO                 

               
            ! ------------------------------------   
               DO i = 1, NN
                  summation = r4(i, 1) * r1(i, 1) + r4(i, 2) * r1(i, 2) + r4(i, 3) * r1(i, 3) 
                  t4(i,1) = (n4(i,1) + n1(i,1))/(n4(i,1) * n1(i,1) * (n4(i,1) * n1(i,1) + summation + O%FVM_Other%Ground%co ** 2))
               END DO    ! t4(:,1) = (n4 + n1) / (n4 * n1 * (n4 * n1 + sum(r4 * r1 , 2) + O%FVM_Other%Ground%co ** 2))
               
               DO i = 1, NN
                  c4(i, 1) =  r4(i, 2) * r1(i, 3) - r4(i, 3) * r1(i, 2) 
                  c4(i, 2) =  r4(i, 3) * r1(i, 1) - r4(i, 1) * r1(i, 3) 
                  c4(i, 3) =  r4(i, 1) * r1(i, 2) - r4(i, 2) * r1(i, 1)     
               END DO  

               DO i = 1, NN
                  A4(i, 1) = c4(i, 1) * t4(i, 1)
                  A4(i, 2) = c4(i, 2) * t4(i, 1)
                  A4(i, 3) = c4(i, 3) * t4(i, 1)
               END DO                 

             ! ------
               O%FVM_Other%Ground%n_vec(J,:) = c2(1, :) / SQRT(c2(1,1)**2 + c2(1,2)**2 + c2(1,3)**2 )
               
               n_vec_transpose(1,1) = O%FVM_Other%Ground%n_vec(J,1)
               n_vec_transpose(2,1) = O%FVM_Other%Ground%n_vec(J,2)
               n_vec_transpose(3,1) = O%FVM_Other%Ground%n_vec(J,3) 
               
               A(:, j) = MATMUL ( (A1 + A2 + A3 + A4) / 4 / PI ,  n_vec_transpose(1:3,1) )
               
               
               j = j + 1
            
         END DO
      END DO
      
      CALL Numerical_inverse(A, O%FVM_Other%Ground%b, Np*NP)

               
               
               
     ! Deallocate
      IF (ALLOCATED ( r1 )) DEALLOCATE ( r1 )
      IF (ALLOCATED ( r2 )) DEALLOCATE ( r2 )     
      IF (ALLOCATED ( r3 )) DEALLOCATE ( r3 )
      IF (ALLOCATED ( r4 )) DEALLOCATE ( r4 )      
      IF (ALLOCATED ( n1 )) DEALLOCATE ( n1 )      
      IF (ALLOCATED ( n2 )) DEALLOCATE ( n2 )      
      IF (ALLOCATED ( n3 )) DEALLOCATE ( n3 )      
      IF (ALLOCATED ( n4 )) DEALLOCATE ( n4 )      
      IF (ALLOCATED ( t1 )) DEALLOCATE ( t1 )      
      IF (ALLOCATED ( t2 )) DEALLOCATE ( t2 )      
      IF (ALLOCATED ( t3 )) DEALLOCATE ( t3 )      
      IF (ALLOCATED ( t4 )) DEALLOCATE ( t4 )      
      IF (ALLOCATED ( c1 )) DEALLOCATE ( c1 )
      IF (ALLOCATED ( c2 )) DEALLOCATE ( c2 )      
      IF (ALLOCATED ( c3 )) DEALLOCATE ( c3 )      
      IF (ALLOCATED ( c4 )) DEALLOCATE ( c4 )      
      IF (ALLOCATED ( a1 )) DEALLOCATE ( a1 )      
      IF (ALLOCATED ( a2 )) DEALLOCATE ( a2 )   
      IF (ALLOCATED ( a3 )) DEALLOCATE ( a3 )   
      IF (ALLOCATED ( a4 )) DEALLOCATE ( a4 )   
      IF (ALLOCATED ( n_vec_transpose )) DEALLOCATE ( n_vec_transpose )         
      IF (ALLOCATED ( a )) DEALLOCATE ( a )  
      IF (ALLOCATED ( x )) DEALLOCATE ( x )  
      IF (ALLOCATED ( y )) DEALLOCATE ( y )  
      IF (ALLOCATED ( z )) DEALLOCATE ( z )        
      
   END SUBROUTINE ground_influence_coefficients
   
   !==================================================================================================================================
   SUBROUTINE Numerical_inverse(a, Inv_a, n)
   ! Computing Inverse matrix. Based on the Doolittle LU method
   ! Reference: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
   !            http://en.wikipedia.org/wiki/LU_decomposition#Doolittle_algorithm
   !....................................................................

      IMPLICIT NONE 

         ! Passed variables
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE,   INTENT(INOUT)  :: a           ! Square matrix
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE,   INTENT(INOUT)  :: Inv_a       ! Inverse of a
      INTEGER(IntKi),                            INTENT(IN   )  :: n           ! The matrix a has size (N,N)

       ! Internal variables
      REAL(DbKi)                :: L(n,n), U(n,n), b(n), d(n), x(n)
      REAL(DbKi)                :: coeff
      INTEGER(IntKi)            :: i, j, k

      !.........................................
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L = 0.0
      U = 0.0
      b = 0.0

      !.........................................
      ! step 1: forward elimination
      DO k = 1, n-1
         DO i = k+1,n
            coeff = a(i,k) / a(k,k)
            L(i,k) = coeff
            DO j = k+1,n
               a(i,j) = a(i,j) - coeff * a(k,j)
            END DO
         END DO
      END DO

      !.........................................
       ! Step 2: prepare L and U matrices 
       ! L matrix is a matrix of the elimination coefficient
       ! + the diagonal elements are 1.0
      DO i=1,n
        L(i,i) = 1.0
      END DO
       ! U matrix is the upper triangular part of A
      DO j = 1,n
        DO i = 1,j
          U(i,j) = a(i,j)
        END DO
      END DO

       !.........................................
       ! Step 3: compute columns of the inverse matrix Inv_a
      DO k = 1,n
         b(k) = 1.0
         d(1) = b(1)
     
        ! Step 3a: Solve Ld=b using the forward substitution
         DO i = 2,n
            d(i) = b(i)
            DO j = 1,i-1
               d(i) = d(i) - L(i,j) * d(j)
            END DO
         END DO
     
       ! Step 3b: Solve Ux=d using the back substitution
         x(n) = d(n) / U(n,n)
         DO i = n-1,1,-1
            x(i) = d(i)
            DO j = n,i+1,-1
              x(i) = x(i) - U(i,j) * x(j)
            END DO
            x(i) = x(i) / u(i,i)
         END DO
     
        ! Step 3c: fill the solutions x(n) into column k of Inv_a
         DO i=1,n
            Inv_a(i,k) = x(i)
         END DO
         b(k)=0.0
      END DO


   END SUBROUTINE Numerical_inverse   
   !======================================================================
   
   
END SUBROUTINE WINDS_Ground_Model
!==================================================================================================================================
SUBROUTINE WINDS_Kinematics(u, p, O, xd, N, ErrStat, ErrMess)
! Computes the station locations of each blade in the inertial coordinate system
! Called from: AD_CalcOutput  (in AeroDyn.f90)
! (~ kinematics.m but much simpler)
!....................................................................
  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_InputType),            INTENT(IN   )  :: u           ! Inputs at t
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep, counts from 1    
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess


      ! Local variables
   INTEGER(IntKi)      :: IDim         ! Index for dimensions
   INTEGER(IntKi)      :: IBlade       ! Index for blade
   INTEGER(IntKi)      :: IElement     ! Index for blade stations
   INTEGER(IntKi)      :: IElement2    ! Index for blade trailing nodes
   INTEGER(IntKi)      :: ITimestep    ! Index for timesteps         
      
   REAL(DbKi), DIMENSION(p%Element%NElm+1, p%NumBl)      :: BLADE_QUARTER
   REAL(DbKi), DIMENSION(p%Element%NElm+1, p%NumBl)      :: BLADE_TRAIL
   REAL(DbKi), DIMENSION(p%Element%NElm, p%NumBl)        :: BLADE_BOUND
   REAL(DbKi), DIMENSION(p%Element%NElm, p%NumBl)        :: BLADE_END

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS 
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB    


   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl

   
   
   
      ! Aerodynamic center of every blade elememt
   DO IBlade = 1, NB
      DO IElement = 1, NS 
         O%FVM_Other%POS_AEROCENT(1, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%Position(1, IElement)
         O%FVM_Other%POS_AEROCENT(2, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%Position(2, IElement)
         O%FVM_Other%POS_AEROCENT(3, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%Position(3, IElement)
      END DO
   END DO

   
   ! Position of spanwise stations and nodes in inertial coordinate system
  
   DO IBlade = 1, NB
      DO IElement = 1, NS 
         O%FVM_Other%POS_END(1, N, 1, IElement, IBlade) = O%FVM_Other%POS_AEROCENT(1, N, 1, IElement, IBlade) +     &  
                                                          (1 - p%Blade%AEROCEN(IElement)) * p%blade%C(IElement) *   &
                                                          u%InputMarkers(IBlade)%Orientation(2,1,IElement)  
         O%FVM_Other%POS_END(2, N, 1, IElement, IBlade) = O%FVM_Other%POS_AEROCENT(2, N, 1, IElement, IBlade) +     &
                                                          (1 - p%Blade%AEROCEN(IElement)) * p%blade%C(IElement)     &  
                                                           * u%InputMarkers(IBlade)%Orientation(2,2,IElement) 
         O%FVM_Other%POS_END(3, N, 1, IElement, IBlade) = O%FVM_Other%POS_AEROCENT(3, N, 1, IElement, IBlade) +    &
                                                          (1 - p%Blade%AEROCEN(IElement)) * p%blade%C(IElement)    &
                                                           * u%InputMarkers(IBlade)%Orientation(2,3,IElement) 
      END DO
   END DO 
  
   DO IBlade = 1, NB
      DO IElement = 1, NS 
         O%FVM_Other%POS_BOUND(1, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%Position(1,IElement)
         O%FVM_Other%POS_BOUND(2, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%Position(2,IElement)
         O%FVM_Other%POS_BOUND(3, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%Position(3,IElement)
         
      END DO
   END DO   

   

       ! 2nd order interpolation
   DO IBlade = 1, NB

      DO IDim = 1, 3
         O%FVM_Other%POS_QUARTER(IDim, N, 1, 1, IBlade) =                                                           & 
                                                    QUADRATIC_INTERP(O%FVM_Other%POS_BOUND(IDim, N, 1, 1, IBlade),  &  
                                                                     O%FVM_Other%POS_BOUND(IDim, N, 1, 2, IBlade),  &
                                                                     O%FVM_Other%POS_BOUND(IDim, N, 1, 3, IBlade),  &
                                                                     p%BLADE%RNodes(1), p%BLADE%RNodes(2), p%BLADE%RNodes(3), &
                                                                     p%BLADE%RTrail(1))
         

         DO IElement2 = 2, NST - 2
            O%FVM_Other%POS_QUARTER(IDim, N, 1, IElement2, IBlade) =                                                              &  
                                                       QUADRATIC_INTERP(O%FVM_Other%POS_BOUND(IDim, N, 1, IElement2 -1, IBlade),  &    
                                                                        O%FVM_Other%POS_BOUND(IDim, N, 1, IElement2, IBlade),     &
                                                                        O%FVM_Other%POS_BOUND(IDim, N, 1, IElement2 +1, IBlade),  &
                                                                        p%BLADE%RNodes(IElement2 -1),                             &
                                                                        p%BLADE%RNodes(IElement2),                                &
                                                                        p%BLADE%RNodes(IElement2 +1),                             &
                                                                        p%BLADE%RTrail(IElement2)  )
         END DO             
             
         O%FVM_Other%POS_QUARTER(IDim, N, 1, NST -1, IBlade) = QUADRATIC_INTERP(O%FVM_Other%POS_BOUND(IDim, N, 1, NST -3, IBlade),  &  
                                                                          O%FVM_Other%POS_BOUND(IDim, N, 1, NST -2, IBlade),        &
                                                                          O%FVM_Other%POS_BOUND(IDim, N, 1, NST -1, IBlade),        &
                                                                          p%BLADE%RNodes(NST - 3), p%BLADE%RNodes(NST - 2),         &
                                                                          p%BLADE%RNodes(NST - 1), p%BLADE%RTrail(NST -1) )
         
         O%FVM_Other%POS_QUARTER(IDim, N, 1, NST, IBlade) = QUADRATIC_INTERP(O%FVM_Other%POS_BOUND(IDim, N, 1, NST -3, IBlade),      &  
                                                                          O%FVM_Other%POS_BOUND(IDim, N, 1, NST -2, IBlade),         &
                                                                          O%FVM_Other%POS_BOUND(IDim, N, 1, NST -1, IBlade),         &
                                                                          p%BLADE%RNodes(NST - 3), p%BLADE%RNodes(NST - 2),          &
                                                                          p%BLADE%RNodes(NST - 1), p%BLADE%RTrail(NST) )                  
      END DO
   END DO         

   
   DO IBlade = 1, NB
      DO IDim = 1, 3
         O%FVM_Other%POS_TRAIL(IDim, N, 1, 1, IBlade) = QUADRATIC_INTERP(O%FVM_Other%POS_END(IDim, N, 1, 1, IBlade),          &  
                                                                     O%FVM_Other%POS_END(IDim, N, 1, 2, IBlade),              &
                                                                     O%FVM_Other%POS_END(IDim, N, 1, 3, IBlade),              &
                                                                     p%BLADE%RNodes(1), p%BLADE%RNodes(2), p%BLADE%RNodes(3), &
                                                                     p%BLADE%RTrail(1))
         

         DO IElement2 = 2, NST - 2
            O%FVM_Other%POS_TRAIL(IDim, N, 1, IElement2, IBlade) =                                                                 & 
                                                          QUADRATIC_INTERP(O%FVM_Other%POS_END(IDim, N, 1, IElement2 -1, IBlade),  &    
                                                                           O%FVM_Other%POS_END(IDim, N, 1, IElement2, IBlade),     &
                                                                           O%FVM_Other%POS_END(IDim, N, 1, IElement2 +1, IBlade),  &
                                                                           p%BLADE%RNodes(IElement2 -1),                           &
                                                                           p%BLADE%RNodes(IElement2),                              &
                                                                           p%BLADE%RNodes(IElement2 +1),                           &
                                                                           p%BLADE%RTrail(IElement2) )
         END DO             
             
         O%FVM_Other%POS_TRAIL(IDim, N, 1, NST -1, IBlade) = QUADRATIC_INTERP(O%FVM_Other%POS_END(IDim, N, 1, NST -3, IBlade),  &  
                                                                          O%FVM_Other%POS_END(IDim, N, 1, NST -2, IBlade),      &
                                                                          O%FVM_Other%POS_END(IDim, N, 1, NST -1, IBlade),      &
                                                                          p%BLADE%RNodes(NST - 3), p%BLADE%RNodes(NST - 2),     &
                                                                          p%BLADE%RNodes(NST - 1), p%BLADE%RTrail(NST -1) )
         
         O%FVM_Other%POS_TRAIL(IDim, N, 1, NST, IBlade) = QUADRATIC_INTERP(O%FVM_Other%POS_END(IDim, N, 1, NST -3, IBlade),     &  
                                                                          O%FVM_Other%POS_END(IDim, N, 1, NST -2, IBlade),      &
                                                                          O%FVM_Other%POS_END(IDim, N, 1, NST -1, IBlade),      &
                                                                          p%BLADE%RNodes(NST - 3), p%BLADE%RNodes(NST - 2),     &
                                                                          p%BLADE%RNodes(NST - 1), p%BLADE%RTrail(NST) )                  
      END DO
   END DO     
   

CONTAINS
   ! ===============================================================
   FUNCTION QUADRATIC_INTERP(X0, X1, X2, L0, L1, L2, Lx)
   ! Numerical Methods for Engineers, 6th edition. pp: 491
      REAL(DbKi)               :: QUADRATIC_INTERP
      REAL(DbKi), INTENT(IN)   :: X0, X1, X2, L0, L1, L2, Lx
      
      ! Local
      REAL(DbKi) :: B0, B1, B2
      
      B0 = X0
      B1 = (X1 - X0)/(L1 - L0)
      B2 = ((X2 - X1)/(L2 - L1) - (X1 - X0)/(L1 - L0)) / (L2 - L0)
      QUADRATIC_INTERP = B0 + B1 * (Lx - L0) + B2 * (Lx - L0) * (Lx - L1)
      

   END FUNCTION QUADRATIC_INTERP
   ! ===============================================================
   
END SUBROUTINE WINDS_Kinematics

!==================================================================================================================================
SUBROUTINE WINDS_Velocity(u, p, O, xd, N, ErrStat, ErrMess)
! Computes the velocity contributions due to turbine and platform motions and freestream flow in the inertial and blade coordinate systems.
! Called from: AD_CalcOutput  (in AeroDyn.f90)
! (~ velocity.m but much simpler)
!....................................................................   

  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_InputType),            INTENT(IN   )  :: u           ! Inputs at t
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(ad_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep, counts from 1    
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess

      ! Local variables
   INTEGER(IntKi)      :: IDim    ! For dimensions
   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: IElement2
   INTEGER(IntKi)      :: ITimestep
   

   
   REAL(DbKi), DIMENSION(1:3)          :: TEMP1 
   REAL(DbKi), DIMENSION(1:3)          :: TEMP2
   REAL(DbKi), DIMENSION(1:3)          :: TEMP3
   REAL(DbKi), DIMENSION(1:3)          :: TEMP4 
   REAL(DbKi), DIMENSION(1:3)          :: CROSS_TEMP
   REAL(DbKi), DIMENSION(1:3)          :: CROSS_TEMP2
   
   INTEGER(IntKi)      :: NST
   INTEGER(IntKi)      :: NS 
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB  
   
   REAL(DbKi)          :: temp
   
   REAL(DbKi)          :: point(3), U_inf(3), temp_vel(3)      

  
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   

   ! Transfer the velocity of aerodynamic center (in the inertial coordinate system)     ! sliu: it is useless but exists in Matlab code
   DO IBlade = 1, NB           
      DO IElement = 1, NS 
         O%FVM_Other%VEL_AEROCENT(1, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%TranslationVel(1,IElement)  ! TranslationVel is in the inertial coordinate system (Aerodynamic center)
         O%FVM_Other%VEL_AEROCENT(2, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%TranslationVel(2,IElement)
         O%FVM_Other%VEL_AEROCENT(3, N, 1, IElement, IBlade) = u%InputMarkers(IBlade)%TranslationVel(3,IElement)
      END DO
   END DO
   
   
     ! Compute kinematically-derived inertial velocities using central differencing
   CALL Bwdiff_bound(p, O )    ! <---- get : O%FVM_Other%VEL_BOUND

   DO IBlade = 1, NB
      DO IElement = 1, NS 
         IF (p%FVM%Shear_Parms%ShearFLAG) THEN             
            point(1) = O%FVM_Other%POS_BOUND(1, N, 1, IElement, IBlade)
            point(2) = O%FVM_Other%POS_BOUND(2, N, 1, IElement, IBlade)
            point(3) = O%FVM_Other%POS_BOUND(3, N, 1, IElement, IBlade)
                      
            U_inf(1) = O%FVM_Other%WIND_INFTY(1, N, 1, 1, IBlade)
            U_inf(2) = O%FVM_Other%WIND_INFTY(2, N, 1, 1, IBlade)
            U_inf(3) = O%FVM_Other%WIND_INFTY(3, N, 1, 1, IBlade)
            
            temp_vel = 0.0_DbKi
          
            CALL Shear_calc(point, U_inf, temp_vel, p, O, N, ErrStat, ErrMess) 
            IF (ErrStat /= 0 ) RETURN
            
            O%FVM_Other%VEL_BOUND(1, N, 1, IElement, IBlade) = - O%FVM_Other%VEL_BOUND(1, N, 1, IElement, IBlade) + temp_vel(1)
            O%FVM_Other%VEL_BOUND(2, N, 1, IElement, IBlade) = - O%FVM_Other%VEL_BOUND(2, N, 1, IElement, IBlade) + temp_vel(2)
            O%FVM_Other%VEL_BOUND(3, N, 1, IElement, IBlade) = - O%FVM_Other%VEL_BOUND(3, N, 1, IElement, IBlade) + temp_vel(3)   
             
         ELSE 
            O%FVM_Other%VEL_BOUND(1, N, 1, IElement, IBlade) = - O%FVM_Other%VEL_BOUND(1, N, 1, IElement, IBlade)   &
                                                               + O%FVM_Other%WIND_INFTY(1, N, 1, IElement, IBlade)
            O%FVM_Other%VEL_BOUND(2, N, 1, IElement, IBlade) = - O%FVM_Other%VEL_BOUND(2, N, 1, IElement, IBlade)   & 
                                                               + O%FVM_Other%WIND_INFTY(2, N, 1, IElement, IBlade)
            O%FVM_Other%VEL_BOUND(3, N, 1, IElement, IBlade) = - O%FVM_Other%VEL_BOUND(3, N, 1, IElement, IBlade)   & 
                                                               + O%FVM_Other%WIND_INFTY(3, N, 1, IElement, IBlade)
         END IF !(p%FVM%Shear_Parms%ShearFLAG)
      END DO
   END DO
   


   ! These valuables are used in subroutine KuttaJoukowski
   O%FVM_Other%POS_NODES_BXN = 0.0_DbKi
   O%FVM_Other%POS_NODES_BYN = 0.0_DbKi
   O%FVM_Other%POS_NODES_BZN = 0.0_DbKi
      
 
   DO IBlade = 1, NB
      DO IElement = 1, NS
         DO IDim = 1, 3
            O%FVM_Other%POS_NODES_BXN(IDim,N,1,IElement,IBlade)  =  u%InputMarkers(IBlade)%Orientation(2,IDim,IElement)    ! CoordS%te2 in ElastoDyn
         END DO             
      END DO
   END DO      
   
   
   DO IBlade = 1, NB
      DO IElement = 1, NS
         temp = 0.0_DbKi
         DO IDim = 1, 3
            temp = temp + (O%FVM_Other%POS_NODES_BXN(IDim,N,1,IElement,IBlade)) ** 2
         END DO
         temp = SQRT(temp)
         DO IDim = 1, 3 
            O%FVM_Other%POS_NODES_BXN(IDim,N,1,IElement,IBlade)  =  O%FVM_Other%POS_NODES_BXN(IDim,N,1,IElement,IBlade) / temp
         END DO
      END DO
   END DO    
   
   

   DO IElement = 1, NS-1
      O%FVM_Other%POS_NODES_BZN(:, N, 1, IElement,:) = O%FVM_Other%POS_BOUND(:, N, 1, IElement +1,:) -       &
                                                      O%FVM_Other%POS_BOUND(:, N, 1, IElement,:) 
   END DO
   

   DO IBlade = 1, NB
      DO IElement = 1, NS -1
         temp = 0.0_DbKi
         DO IDim = 1, 3
            temp = temp + (O%FVM_Other%POS_NODES_BZN(IDim,N,1,IElement,IBlade)) ** 2
         END DO
         temp = SQRT(temp)
         DO IDim = 1, 3 
            O%FVM_Other%POS_NODES_BZN(IDim,N,1,IElement,IBlade)  =  O%FVM_Other%POS_NODES_BZN(IDim,N,1,IElement,IBlade) / temp
         END DO
      END DO
   END DO    
   
   
   
   DO IElement = NS-1, 1, -1 
      O%FVM_Other%POS_NODES_BZN(:,N,1,IElement +1,:) =  O%FVM_Other%POS_NODES_BZN(:,N,1,IElement,:)
   END DO
   
   DO IElement = 1, NS      !  another method: O%FVM_Other%POS_NODES_BYN is same direction with u%InputMarkers(IBlade)%Orientation(1,:,:)  and CoordS%te1 in ElastoDyn
      DO IBlade = 1, NB
         temp3(1) = O%FVM_Other%POS_NODES_BZN(1,N,1,IElement,IBlade)
         temp3(2) = O%FVM_Other%POS_NODES_BZN(2,N,1,IElement,IBlade)
         temp3(3) = O%FVM_Other%POS_NODES_BZN(3,N,1,IElement,IBlade)
         temp4(1) = O%FVM_Other%POS_NODES_BXN(1,N,1,IElement,IBlade)
         temp4(2) = O%FVM_Other%POS_NODES_BXN(2,N,1,IElement,IBlade)
         temp4(3) = O%FVM_Other%POS_NODES_BXN(3,N,1,IElement,IBlade)         
         
         CROSS_TEMP2 = CROSS_PRODUCT_Db(temp3, temp4)
         O%FVM_Other%POS_NODES_BYN(1,N,1,IElement,IBlade) = CROSS_TEMP2(1)
         O%FVM_Other%POS_NODES_BYN(2,N,1,IElement,IBlade) = CROSS_TEMP2(2)
         O%FVM_Other%POS_NODES_BYN(3,N,1,IElement,IBlade) = CROSS_TEMP2(3)
      END DO
   END DO
   
         ! VEL_BLADE: the velocites along the blade quarter-chord (lifting-line) in the blade coordinate system
   O%FVM_Other%VEL_BLADE(1,N,1,:,:)  =  O%FVM_Other%POS_NODES_BXN(1,N,1,:,:) *  O%FVM_Other%VEL_BOUND(1,N,1,:,:) + &   
                                        O%FVM_Other%POS_NODES_BXN(2,N,1,:,:) *  O%FVM_Other%VEL_BOUND(2,N,1,:,:) + &
                                        O%FVM_Other%POS_NODES_BXN(3,N,1,:,:) *  O%FVM_Other%VEL_BOUND(3,N,1,:,:)       
   
   O%FVM_Other%VEL_BLADE(2,N,1,:,:)  =  O%FVM_Other%POS_NODES_BYN(1,N,1,:,:) *  O%FVM_Other%VEL_BOUND(1,N,1,:,:) + &
                                        O%FVM_Other%POS_NODES_BYN(2,N,1,:,:) *  O%FVM_Other%VEL_BOUND(2,N,1,:,:) + &
                                        O%FVM_Other%POS_NODES_BYN(3,N,1,:,:) *  O%FVM_Other%VEL_BOUND(3,N,1,:,:)      
   
   O%FVM_Other%VEL_BLADE(3,N,1,:,:)  =  O%FVM_Other%POS_NODES_BZN(1,N,1,:,:) *  O%FVM_Other%VEL_BOUND(1,N,1,:,:) + &
                                        O%FVM_Other%POS_NODES_BZN(2,N,1,:,:) *  O%FVM_Other%VEL_BOUND(2,N,1,:,:) + &
                                        O%FVM_Other%POS_NODES_BZN(3,N,1,:,:) *  O%FVM_Other%VEL_BOUND(3,N,1,:,:)      

   
CONTAINS   
   !....................................................................
   SUBROUTINE Bwdiff_bound(p, O) !, y, dy)
   ! Computes the derivative of the function 'y' wrt 'x' numerically via centered difference method
   ! (~ bwdiff.m )
   !....................................................................

  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState ! Other/optimization states
   !REAL(DbKi),DIMENSION(:,:,:,:), INTENT(IN   )  :: y           ! Array of dependent values such that f(x)=y
   !REAL(DbKi),DIMENSION(:,:,:,:), INTENT(INOUT)  :: dy          ! Array of derivative values of y wrt x
   
   

      ! Local variables
   INTEGER(IntKi)      :: N  
   
   INTEGER(IntKi)      :: NST
   INTEGER(IntKi)      :: NS 
   INTEGER(IntKi)      :: NT  
   INTEGER(IntKi)      :: NB   

   INTEGER(IntKi)      :: IDim    ! Index for dimensions
   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: IElement2
   INTEGER(IntKi)      :: ITimestep
   REAL(DbKi)          :: DT      ! Duration of between timesteps

   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl 
     
   DT = p%FVM%DT_WINDS
   N  = O%WINDS_Timestep    
   
                              
   IF (N == 1) O%FVM_Other%VEL_BOUND(:,N,1,:,:) = 0.0_DbKi
   
   IF (N == 2) O%FVM_Other%VEL_BOUND(:,N,1,:,:) = (O%FVM_Other%POS_BOUND(:,N,1,:,:) - O%FVM_Other%POS_BOUND(:,N-1,1,:,:))/DT
      
   IF (N == 3) O%FVM_Other%VEL_BOUND(:,N,1,:,:) = (3* O%FVM_Other%POS_BOUND(:,N,1,:,:) - 4* O%FVM_Other%POS_BOUND(:,N-1,1,:,:) +        &
                                                  O%FVM_Other%POS_BOUND(:,N-2,1,:,:))/ 2 / DT

   IF (N == 4) O%FVM_Other%VEL_BOUND(:,N,1,:,:) = (11* O%FVM_Other%POS_BOUND(:,N,1,:,:) - 18* O%FVM_Other%POS_BOUND(:,N-1,1,:,:) +       & 
                                            9* O%FVM_Other%POS_BOUND(:,N-2,1,:,:) - 2* O%FVM_Other%POS_BOUND(:,N-3,1,:,:))/ 6 / DT

   IF (N == 5) O%FVM_Other%VEL_BOUND(:,N,1,:,:) = (25* O%FVM_Other%POS_BOUND(:,N,1,:,:) - 48* O%FVM_Other%POS_BOUND(:,N-1,1,:,:) +       & 
                                            36* O%FVM_Other%POS_BOUND(:,N-2,1,:,:) - 16* O%FVM_Other%POS_BOUND(:,N-3,1,:,:) +      &
                                            3* O%FVM_Other%POS_BOUND(:,N-4,1,:,:))/ 12 / DT

   IF (N >=6 ) O%FVM_Other%VEL_BOUND(:,N,1,:,:) = (137* O%FVM_Other%POS_BOUND(:,N,1,:,:) - 300* O%FVM_Other%POS_BOUND(:,N-1,1,:,:)+     &
                                       300* O%FVM_Other%POS_BOUND(:,N-2,1,:,:) - 200* O%FVM_Other%POS_BOUND(:,N-3,1,:,:) +         &
                                       75* O%FVM_Other%POS_BOUND(:,N-4,1,:,:) - 12 * O%FVM_Other%POS_BOUND(:,N-5,1,:,:))/ 60 / DT                              
                              
                              
   

   END SUBROUTINE Bwdiff_bound
   !....................................................................
   FUNCTION Cross_Product_Db(Vector1, Vector2)

      ! This function computes the cross product of two 3-element arrays:
      ! Cross_Product = Vector1 X Vector2 (resulting in a vector)


      ! Argument declarations.

   REAL(DbKi), INTENT(IN )         :: Vector1       (3)
   REAL(DbKi), INTENT(IN )         :: Vector2       (3)

      ! Function definition
   REAL(DbKi)                      :: Cross_Product_Db (3)        ! = Vector1 X Vector2 (resulting in a vector)


   Cross_Product_Db(1) = Vector1(2)*Vector2(3) - Vector1(3)*Vector2(2)
   Cross_Product_Db(2) = Vector1(3)*Vector2(1) - Vector1(1)*Vector2(3)
   Cross_Product_Db(3) = Vector1(1)*Vector2(2) - Vector1(2)*Vector2(1)


   
   END FUNCTION Cross_Product_Db
   !....................................................................
   
END SUBROUTINE WINDS_Velocity

!==================================================================================================================================
SUBROUTINE WINDS_FVMInitial(u, p, O, xd, ErrStat, ErrMess, x, z, y)
! Use blade element and momentum method to calculate the induced velocity at first timestep
! Similar with "initials.m" in WInDS by Thomas Sebastian
! Called from: AD_CalcOutput  (in AeroDyn.f90)
! (~ initials.m / BEM.m)
!....................................................................
   IMPLICIT                        NONE


      ! Passed variables  
   TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O           ! Other/optimization states
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None

   TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
   TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
   TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
                                                                       
      ! Internal variables
   INTEGER(IntKi)             :: IDim       ! Index for dimensions
   INTEGER(IntKi)             :: IBlade     ! Index for blade
   INTEGER(IntKi)             :: IElement   ! Index for blade stations
   INTEGER(IntKi)             :: IElement2  ! Index for blade trailing nodes
   INTEGER(IntKi)             :: ITimestep            

   INTEGER(IntKi)             :: N   !=  O%WINDS_Timestep 
   INTEGER(IntKi)             :: NST != p%FVM%NST
   INTEGER(IntKi)             :: NS  != p%FVM%NS
   INTEGER(IntKi)             :: NT  != p%FVM%NT
   INTEGER(IntKi)             :: NB  != p%FVM%NB

   REAL(DbKi)                 :: Temp1, Temp2
   
   CHARACTER(LEN = 1024)      :: LINE   
   CHARACTER(LEN = 1024)      :: TEMP
   CHARACTER(LEN = 1024)      :: filename 
   
   
   N   =  O%WINDS_Timestep       
   NST = p%Element%NElm + 1
   NB  = p%NumBl
   NS  = p%Element%NElm
   
   
   ! IF ( p%FVM%SteadyFlag ) THEN ! Make it steady flow         
   O%FVM_Other%Wind_Mean = 0.0_DbKi
      
   DO IBlade = 1,p%NumBl
       DO IElement = 1,p%Element%Nelm
           O%FVM_Other%Wind_Mean(1) = O%FVM_Other%Wind_Mean(1) + O%FVM_Other%WIND_INFTY(1, 1, 1, IElement, IBlade)
           O%FVM_Other%Wind_Mean(2) = O%FVM_Other%Wind_Mean(2) + O%FVM_Other%WIND_INFTY(2, 1, 1, IElement, IBlade)
           O%FVM_Other%Wind_Mean(3) = O%FVM_Other%Wind_Mean(3) + O%FVM_Other%WIND_INFTY(3, 1, 1, IElement, IBlade)
       END DO
   END DO
            
   O%FVM_Other%Wind_Mean(:) = O%FVM_Other%Wind_Mean(:) / p%NumBl / p%Element%Nelm

   DO IBlade = 1,p%NumBl
       DO IElement=1,p%Element%Nelm   
           O%FVM_Other%WIND_INFTY(1, 1, 1, IElement, IBlade) = O%FVM_Other%Wind_Mean(1)
           O%FVM_Other%WIND_INFTY(2, 1, 1, IElement, IBlade) = O%FVM_Other%Wind_Mean(2)
           O%FVM_Other%WIND_INFTY(3, 1, 1, IElement, IBlade) = O%FVM_Other%Wind_Mean(3)  
           
           O%FVM_Other%WIND_INFTYM(1, 1, 1, IElement, IBlade)  =     &
                                  SQRT(O%FVM_Other%Wind_Mean(1)**2 + O%FVM_Other%Wind_Mean(2)**2 + O%FVM_Other%Wind_Mean(3)**2 )                   
       END DO
   END DO             
   !END IF  ! p%FVM%SteadyFlag  
   
     

     ! Substitute in initial values and truncate size of variables by timestep
   O%FVM_Other%WAKE_DOMAIN(:,1,1,:,:)  =  O%FVM_Other%POS_QUARTER(:,1,1,:,:)  
   O%FVM_Other%WAKE_DOMAIN(:,2,1,:,:)  =  O%FVM_Other%POS_TRAIL(:,1,1,:,:)   
   
   
   CALL shear(O, p, 'MAIN', 1 ,N, ErrStat, ErrMess)   
   IF (ErrStat /= 0 ) RETURN
   
   O%FVM_Other%VEL_DOMAIN_RK(:,:,1,:,:)  =  O%FVM_Other%VEL_DOMAIN(:,:,1,:,:)
   O%FVM_Other%VEL_DOMAIN_RK(:,:,2,:,:)  =  O%FVM_Other%VEL_DOMAIN(:,:,1,:,:)
   O%FVM_Other%VEL_DOMAIN_RK(:,:,3,:,:)  =  O%FVM_Other%VEL_DOMAIN(:,:,1,:,:) 
   
            
      ! Define initial induced velocities via 1st-order methods         
   CALL BEM(u, p, xd, O, ErrStat, ErrMess)  
   IF (ErrStat /= 0 ) RETURN
   
      ! Define initial vortex strength      
            ! Use Kutta-Joukowski theorem to define bound circulation strength         
   DO IBlade = 1, NB     
      DO IElement = 1, NS   
         O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement, IBlade) = 0.5 * O%FVM_Other%WIND_INFTYM(1,1,1,1,1) * p%Blade%C(IElement) *     &
                                                                    O%FVM_Other%PERF_CL(1, 1, 1, IElement, 1)     
     
      END DO !IElement
   END DO !IBlade      



            ! Compute spanwise change in bound filament to compute first set of trailing filaments            
   DO IBlade = 1,NB    
      DO IElement2 = 1, NST     
         IF (IElement2 == 1) THEN ! First element
            O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1, 1, IElement2, IBlade) = O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2, IBlade)
            
         ELSE IF (IElement2 == NST) THEN  ! Last element
            O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1, 1, IElement2, IBlade) = -O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2-1, IBlade) 
            
         ELSE ! Not first or last
            O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1, 1, IElement2, IBlade) = O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2, IBlade) &
                                                                     - O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2 -1, IBlade) 
         ENDIF  
      END DO !IElement2        
   END DO !IBlade      



            ! Shed filaments computed via spanwise summation of trailing filaments (ensure Kelvin's theorem is satisfied)
     Do IBlade = 1, NB  
         Temp1 = 0.0 
         
         Do IElement = 1, NS 
            Temp1 = Temp1 + O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1, 1, IElement, IBlade)
            O%FVM_Other%WAKE_GAMMA_SHED(1, 2, 1, IElement, IBlade) =  - Temp1
         END DO
     END DO            
            
  
   
      ! Modify core size using Ramasamy-Leishman model
   CALL VCORE(p, O, xd, 1, ErrStat, ErrMess)
   IF (ErrStat /= 0 ) RETURN

   
   ! Output First Timestep to Paraview file    
   IF (p%FVM%AnimFLAG) THEN
      DO IBlade = 1, p%FVM%NB  
         CALL CreateVTUembedded(IBlade, filename, P, xd, O, ErrStat, ErrMess )
         IF (ErrStat /= 0 ) RETURN
            
         WRITE(TEMP, "(F0.3)") ( (O%WINDS_Timestep - 1) * p%FVM%DT_WINDS )  
      
         LINE = '<DataSet timestep="' // TEMP(1:9) // '" group="" part="0" file="'//  TRIM(filename) //'"/>'
         WRITE(IBlade+10, "(A)") (TRIM(LINE))         
      END DO   
   END IF
   
   
END SUBROUTINE WINDS_FVMInitial

!==================================================================================================================================
SUBROUTINE WINDS_FVM(u, p, O, xd, N, ErrStat, ErrMess)
! Use free vortex wake method to calculate the aerodynamic loads after first timestep
! Called from: AD_CalcOutput (in AeroDyn.f90)
! (~ WInDS.m  Line.185-212 )
!....................................................................

  IMPLICIT                        NONE

      ! Passed variables          
   TYPE(AD_InputType),            INTENT(IN   )  :: u           ! Inputs at t
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep   
   TYPE(aD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states

   INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   INTEGER(IntKi)      :: IDim    ! For dimensions
   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: IElement2
   INTEGER(IntKi)      :: ITimestep

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB    
   
   CHARACTER(LEN = 1024)   :: LINE   
   CHARACTER(LEN = 1024)   :: TEMP
   CHARACTER(LEN = 1024)   :: filename 
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   
   
      ! Store previous time step data   
    CALL CUR_2_PREV(O, p, xd, N, ErrStat, ErrMess)      
   IF (ErrStat /= 0 ) RETURN

      ! Update velocity domain based on last wake domain
    CALL Shear(O, p, 'MAIN', 1 ,N, ErrStat, ErrMess)  
    IF (ErrStat /= 0 ) RETURN

      ! Numerically convect wake nodes to time+1   
    IF  ( N /= p%FVM%NT ) THEN 
       CALL NUM_ADVECT(O, p, xd, N, ErrStat, ErrMess) 
       IF (ErrStat /= 0 ) RETURN
    END IF
              
      ! Modify vortex core size via Ramasamy-Leishman model and include effect of filament stretching from previous timestep
   CALL VCORE(p, O, xd, N, ErrStat, ErrMess)   
   IF (ErrStat /= 0 ) RETURN
   
      ! Compute strength of new bound vortex via Kutta-Joukowski theorem
   CALL KuttaJoukowski(O, p, xd, N, ErrStat, ErrMess)
   IF (ErrStat /= 0 ) RETURN

      ! Compute performance(loads)
   CALL Aero_Loads(u, p, xd, O, ErrStat, ErrMess)  
   IF (ErrStat /= 0 ) RETURN
   
   IF (p%FVM%AnimFLAG) THEN
      DO IBlade = 1, p%FVM%NB  
         CALL CreateVTUembedded(IBlade, filename, P, xd, O, ErrStat, ErrMess )
         IF (ErrStat /= 0 ) RETURN
            
         WRITE(TEMP, "(F7.3)") ( (O%WINDS_Timestep - 1) * p%FVM%DT_WINDS )  
         LINE = '<DataSet timestep="' // TRIM(TEMP) // '" group="" part="0" file="' // TRIM(filename) //'"/>'
         WRITE(IBlade+10, "(A)") (TRIM(LINE))         
      END DO   
   END IF


END SUBROUTINE WINDS_FVM

!==================================================================================================================================
SUBROUTINE WINDS_Shear_Model(p, O, ErrStat, ErrMess)
! ~ shear_model.m
!  Defines parameters used in shear model
!....................................................................

  IMPLICIT                        NONE

      ! Passed variables
   TYPE(AD_ParameterType),        INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState 
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess   


   IF (p%FVM%Shear_Parms%model_type == 1) THEN       
      p%FVM%Shear_Parms%Alpha = 0.144         ! Power law exponent
      
   ELSEIF (p%FVM%Shear_Parms%model_type == 2) THEN       
      p%FVM%Shear_Parms%z0 = 0.01             ! Roughness length
      
   END IF

END SUBROUTINE WINDS_Shear_Model

!==================================================================================================================================
SUBROUTINE Aero_Loads(u, p, xd, O, ErrStat, ErrMess)
! Calculate the aerodynamic loads on the blades (as "perform" function in Matlab WInDS)
! Called from: WINDS_FVM (in WINDS.f90)
! (~ perform.m)
!....................................................................

  IMPLICIT                        NONE


      ! Passed variables

   TYPE(AD_InputType),            INTENT(IN   )  :: u           ! Inputs at Time      
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O           ! Other/optimization states
   INTEGER(IntKi),                INTENT(INOUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi)      :: IDim    ! For dimensions
   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: ITimestep
   
   REAL(DbKi), DIMENSION(SIZE(O%FVM_Other%VEL_BLADE(:,:,1,1,:)))      :: velX
   
      ! From AeroDyn
   REAL(DbKi)                 :: Elem_pitch ! element pitch     
   REAL(DbKi)                 :: PHI     ! Inflow angle  
   REAL(DbKi)                 :: QA      ! Force base factor
   REAL(DbKi)                 :: CPHI    ! Cos(PHI)
   REAL(DbKi)                 :: SPHI    ! Sin(PHI)
   REAL(DbKi)                 :: SPitch  ! Sin(Pitch)
   REAL(DbKi)                 :: CPitch  ! Cos(Pitch)
   REAL(DbKi)                 :: W2      ! U^2 + V^2

   REAL(DbKi)                 :: CLA
   REAL(DbKi)                 :: CDA
   REAL(DbKi)                 :: PMA
   REAL(DbKi)                 :: CMA  
   REAL(DbKi)                 :: DFN
   REAL(DbKi)                 :: DFT
   
   
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS 
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB     
   INTEGER(IntKi)      :: N   


   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   
   N = O%WINDS_Timestep  
   

      ! From AeroDyn(Subroutine ELEMFRC)   
   DO IElement = 1, NS
      DO IBlade = 1, NB  
         Elem_pitch = -1.*ATAN2( -1.*DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(2,:,IElement) ) , &
                                              DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(1,:,IElement) )   ) 
                    
         PHI      = O%FVM_Other%PERF_AOA(1, N, 1, IElement, IBLADE) + Elem_pitch  ! Elem_pitch is blade element pitch, so no Twist angle needed  .... it has opposite sign with Matlab WInDS, but same with BEM in AeroDyn
         CLA      = O%FVM_Other%PERF_CL(1, N, 1, IElement, IBLADE) 
         CDA      = O%FVM_Other%PERF_CD(1, N, 1, IElement, IBLADE) 
         CMA      = O%FVM_Other%PERF_CM(1, N, 1, IElement, IBLADE)      
         
         W2 = O%FVM_Other%KJ%VEL_TOT(1, N, 1, IElement, IBLADE) **2  + O%FVM_Other%KJ%VEL_TOT(2, N, 1, IElement, IBLADE) **2
                  
         QA       = 0.5 * P%Wind%RHO * W2 * P%Blade%C(IElement) ! sliu: do not mutiply P%Blade%DR(IElement) here and do not divide it in Forces (protential bug in AeroDyn)...

         CPHI     = COS( PHI )
         SPHI     = SIN( PHI )
         DFN      = ( CLA * CPHI + CDA * SPHI ) * QA
         DFT      = ( CLA * SPHI - CDA * CPHI ) * QA

         IF ( P%PMOMENT ) THEN
            PMA  = CMA * QA * P%Blade%C(IElement)
         ELSE
            PMA  = 0.0_DbKi
            CMA  = 0.0_DbKi
         ENDIF
         
         SPitch    = SIN( Elem_pitch )
         CPitch    = COS( Elem_pitch )
         
         
         O%FVM_Other%StoredForces(1, N, 1, IElement, IBLADE)   = ( DFN*CPitch + DFT*SPitch ) 
         O%FVM_Other%StoredForces(2, N, 1, IElement, IBLADE)   = ( DFN*SPitch - DFT*CPitch ) 
         O%FVM_Other%StoredForces(3, N, 1, IElement, IBLADE)   = 0.0_DbKi

         O%FVM_Other%StoredMoments(1, N, 1, IElement, IBLADE)  = 0.0_DbKi
         O%FVM_Other%StoredMoments(2, N, 1, IElement, IBLADE)  = 0.0_DbKi
         O%FVM_Other%StoredMoments(3, N, 1, IElement, IBLADE)  = PMA
         
      END DO ! IBlade 
   END DO ! IElement


END SUBROUTINE Aero_Loads

!==================================================================================================================================
SUBROUTINE BiotSavart(F1_all, F2_all, Pts_all, GAMMA_all, RC_all, U_ind_all, P, O, xd, ErrStat, ErrMess)   
! Biot-Savart Law to calculate the induced velocity
! Called from: InducedVelocity (in WINDS.f90)
! (~ BiotSavart.m)
!
! ********Input(s)********
! F1                    Array containing first point of each vortex filament
! F2                    Array containing second point of each vortex filament
! Pts                   Array containing points of interest (where induction is computed)
! gamma                 Array of vortex filament circulation strengths
! rc                    Vortex core sizes (actually radius squared for code speed-up)
! p, O, xd     Contains parameters and states....
!
! ********Output(s)********
! uind      Array of induced velocity at each of the points P due to contributions from filaments defined by F1 and F2
!
!
!
! **************** Pseudo code ****************
!
!IF (p%FVM%OpenMP_Parms%Accelerate) THEN    
!     Reshape the arrays
!     IF (p%FVM%Tree_Parms%TreeFlag .AND. Len_pts > 5000) THEN 
!          Run parallelized Treecode
!          IF (p%FVM%Tree_Parms%Speedup) THEN              
!              Record the speedup and error
!          END IF
!     ELSE                         
!          Use direct parallel computation                    
!     END IF
!ELSE  
!     Run on single core.
!END IF
!....................................................................
   
     IMPLICIT                        NONE
   
         
         ! Passing variables
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:,:,:,:) :: F1_all
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:,:,:,:) :: F2_all
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:,:,:,:) :: Pts_all   
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:,:,:,:) :: GAMMA_all  
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:,:,:,:) :: RC_all      
      REAL(DbKi), INTENT(INOUT), DIMENSION(:,:,:,:,:) :: U_ind_all  
      
      TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
      TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None   
   
   
      ! Internal variables
      INTEGER(IntKi)    :: I_1
      INTEGER(IntKi)    :: I_2  ! Counter of point age index      (points of interest)
      INTEGER(IntKi)    :: I_4  ! Counter of point element index  (points of interest)
      INTEGER(IntKi)    :: I_5  ! Counter of point blade index    (points of interest)
      
      INTEGER(IntKi)    :: J_2  ! Counter of filament age index
      INTEGER(IntKi)    :: J_4  ! Counter of filament element index
      INTEGER(IntKi)    :: J_5  ! Counter of filament blade index   
      
      REAL(DbKi)        :: Time_1, Time_2  
   
      ! For non-acceleration:    
      REAL(DbKi)     :: GMMA, RC, PX, PY, PZ, CO         
      REAL(DbKi)     :: X1, Y1, Z1, X2, Y2, Z2      
      REAL(DbKi)     :: X2X1, Y2Y1, Z2Z1, L       
      REAL(DbKi)     :: PXX1, pyy1, pzz1, pxx2, pyy2, pzz2                                 
      REAL(DbKi)     :: R1, R2, R1DR2, R1TR2, LDR12, CNU, UBAR 
      REAL(DbKi)     :: PRESUMX, PRESUMY, PRESUMZ, DEN

      
      ! For acceleration:       
      INTEGER(IntKi)      :: K, I 
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: F1_new
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: F2_new   
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: Pts_new
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: GAMMA_new
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: RC_new
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: Uind_new   
      INTEGER(IntKi), DIMENSION(5)                   :: Dim_fila   ! DIMENSION of filament
      INTEGER(IntKi), DIMENSION(5)                   :: Dim_pts    ! DIMENSION of filament
      INTEGER(IntKi)                                 :: Len_fila   ! # of filament
      INTEGER(IntKi)                                 :: Len_pts    ! # of points   
   
      
      REAL(DbKi)                                     :: sum_error 
      REAL(DbKi)                                     :: time_original
      REAL(DbKi)                                     :: time_accelerated
      REAL(DbKi)                                     :: speedup_ratio
      REAL(DbKi)                                     :: Error
      
      ! For Treecode
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: fila_x  ! (:,1) Start point; (:,2) End point; (:,3) mid-point;  
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: fila_y  ! (:,1) Start point; (:,2) End point; (:,3) mid-point;
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: fila_z  ! (:,1) Start point; (:,2) End point; (:,3) mid-point;     
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: GAMMA_new_copy
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: RC_new_copy

      
      IF (p%FVM%OpenMP_Parms%Accelerate) THEN ! .....................................................................  

          Dim_fila(1) = 1
          Dim_fila(2) = SIZE(F1_all , 2)
          Dim_fila(3) = 1
          Dim_fila(4) = SIZE(F1_all , 4)
          Dim_fila(5) = SIZE(F1_all , 5)
          Len_fila    = Dim_fila(2)  * Dim_fila(4) * Dim_fila(5)
           
          Dim_pts(1) = 1
          Dim_pts(2) = SIZE(Pts_all , 2)
          Dim_pts(3) = 1
          Dim_pts(4) = SIZE(Pts_all , 4)
          Dim_pts(5) = SIZE(Pts_all , 5)
          Len_pts    = Dim_pts(2)  * Dim_pts(4) * Dim_pts(5)   
          
          
          IF (.NOT. ALLOCATED( F1_new ) )        ALLOCATE( F1_new    (Len_fila, 3))
          IF (.NOT. ALLOCATED( F2_new ) )        ALLOCATE( F2_new    (Len_fila, 3))
          IF (.NOT. ALLOCATED( GAMMA_new ) )     ALLOCATE( GAMMA_new (Len_fila, 1))
          IF (.NOT. ALLOCATED( RC_new ) )        ALLOCATE( RC_new    (Len_fila, 1)) 
          IF (.NOT. ALLOCATED( Pts_new ) )       ALLOCATE( Pts_new   (Len_pts, 3))  
          IF (.NOT. ALLOCATED( Uind_new ) )      ALLOCATE( Uind_new  (Len_pts, 3)) 
           
          
          IF (.NOT. ALLOCATED( fila_x ) ) ALLOCATE( fila_x(Len_fila,3)) 
          IF (.NOT. ALLOCATED( fila_y ) ) ALLOCATE( fila_y(Len_fila,3)) 
          IF (.NOT. ALLOCATED( fila_z ) ) ALLOCATE( fila_z(Len_fila,3)) 
          

          K = 1
          ! Loop the points of interest
          DO I_5 = 1, SIZE(Pts_all , 5)       ! Counter of point blade index     
             DO I_4 = 1, SIZE(Pts_all , 4)      ! Counter of point element index
               DO I_2 = 1, SIZE(Pts_all , 2)      ! Counter of point age index                   
                  Pts_new(K, 1) = Pts_all(1, I_2, 1, I_4, I_5)
                  Pts_new(K, 2) = Pts_all(2, I_2, 1, I_4, I_5)
                  Pts_new(K, 3) = Pts_all(3, I_2, 1, I_4, I_5)      
                  K = K + 1
               END DO
            END DO
         END DO
         
         K = 1
         ! Loop the filaments  
         Do J_5 = 1, SIZE(F1_all , 5)     ! Counter of filament blade index  
            DO J_4 = 1, SIZE(F1_all , 4)      ! Counter of filament element index
               DO J_2 = 1, SIZE(F1_all , 2)     ! Counter of filament age index                         
                  F1_new(K, 1)    =  F1_all   (1, J_2, 1, J_4, J_5)
                  F1_new(K, 2)    =  F1_all   (2, J_2, 1, J_4, J_5)
                  F1_new(K, 3)    =  F1_all   (3, J_2, 1, J_4, J_5)
                  F2_new(K, 1)    =  F2_all   (1, J_2, 1, J_4, J_5)
                  F2_new(K, 2)    =  F2_all   (2, J_2, 1, J_4, J_5)
                  F2_new(K, 3)    =  F2_all   (3, J_2, 1, J_4, J_5)            
                  GAMMA_new(K, 1) =  GAMMA_all(1, J_2, 1, J_4, J_5)
                  RC_new(K, 1)    =  RC_all   (1, J_2, 1, J_4, J_5)
                  
                  
                  ! Rearrange the array for Treecode. This way is easier for Treecode to build the tree structure.
                  fila_x(K,1)        = F1_new(K, 1)   ! Start point, x
                  fila_y(K,1)        = F1_new(K, 2)   ! Start point, y
                  fila_z(K,1)        = F1_new(K, 3)   ! Start point, z
                  fila_x(K,2)        = F2_new(K, 1)   ! End point, x
                  fila_y(K,2)        = F2_new(K, 2)   ! End point, y
                  fila_z(K,2)        = F2_new(K, 3)   ! End point, z    
                  fila_x(K,3)        = F1_new(K, 1)/2 + F2_new(K, 1)/2    ! Mid-point, x
                  fila_y(K,3)        = F1_new(K, 2)/2 + F2_new(K, 2)/2    ! Mid-point, y
                  fila_z(K,3)        = F1_new(K, 3)/2 + F2_new(K, 3)/2    ! Mid-point, z
                  
                  K = K + 1 
               END DO
            END DO
         END DO
           
         
         !!.........................................................
         !! Write internal varible to a file for Treecode.  Sliu: Comment out this part
         !IF (MOD(O%WINDS_Timestep,INT(60/p%FVM%DT_WINDS)) == 0) THEN
         !  CALL WRITE_INTERNAL(fila_x, fila_y, fila_z, Pts_new, GAMMA_new, RC_new, Uind_new, P, xd, O, ErrStat, ErrMess, Len_fila, Len_pts)
         !END IF   
         !.........................................................
         
         
         
            !.......................................................
            ! Use Parallelized Treecode algorithm
         IF (p%FVM%Tree_Parms%TreeFlag .AND. Len_pts > 5000) THEN 
             IF (.NOT. ALLOCATED( GAMMA_new_copy ) )     ALLOCATE( GAMMA_new_copy (Len_fila, 1))
             IF (.NOT. ALLOCATED( RC_new_copy ) )        ALLOCATE( RC_new_copy    (Len_fila, 1)) 
             
             GAMMA_new_copy = GAMMA_new
             RC_new_copy = RC_new
             
             CALL CPU_TIME(Time_1) ! Start time
             CALL BiotSavart_Treecode(fila_x, fila_y, fila_z, Pts_new, GAMMA_new_copy, RC_new_copy, Uind_new, P, O, xd, ErrStat, ErrMess)
             
             CALL CPU_TIME(Time_2) ! End time
             time_accelerated                       =  Time_2 - Time_1 
             O%FVM_Other%TIME%Time_biotsavart_acce  =  O%FVM_Other%TIME%Time_biotsavart_acce + time_accelerated         
             
             IF ( ALLOCATED( GAMMA_new_copy ) )     DEALLOCATE( GAMMA_new_copy )
             IF ( ALLOCATED( RC_new_copy ) )        DEALLOCATE( RC_new_copy  )         
             
             
             !...............................................................
             ! IF user sets to record the speedup and error.
             IF (p%FVM%Tree_Parms%Speedup) THEN
                SELECT CASE (p%FVM%Tree_Parms%Freq) 
                   CASE ('WHOL') 
                        CALL Tree_Record(Uind_new, F1_new, F2_new, Pts_new, GAMMA_new, RC_new, P, O, xd, ErrStat, ErrMess) 
                   CASE ('LAST') 
                       IF (O%WINDS_Timestep == p%FVM%NT-1) THEN
                          ! Just for the last timestep
                           CALL Tree_Record(Uind_new, F1_new, F2_new, Pts_new, GAMMA_new, RC_new, P, O, xd, ErrStat, ErrMess) 
                       END IF  
                   CASE ('TENS')   
                       IF (MOD(O%WINDS_Timestep, 10_IntKi) == 0 ) THEN
                          !  Every 10 timestep
                          CALL Tree_Record(Uind_new, F1_new, F2_new, Pts_new, GAMMA_new, RC_new, P, O, xd, ErrStat, ErrMess)    
                       END IF                         
                   CASE ('HUND') 
                       IF (MOD(O%WINDS_Timestep, 100_IntKi) == 0 ) THEN
                          !  Every 100 timestep
                         CALL Tree_Record(Uind_new, F1_new, F2_new, Pts_new, GAMMA_new, RC_new, P, O, xd, ErrStat, ErrMess) 
                       END IF  
                   CASE ('THOU') 
                       IF (MOD(O%WINDS_Timestep, 1000_IntKi) == 0 ) THEN
                           ! Every 1000 timestep
                          CALL Tree_Record(Uind_new, F1_new, F2_new, Pts_new, GAMMA_new, RC_new, P, O, xd, ErrStat, ErrMess)           
                       END IF                      
                   END SELECT 
                   IF (ErrStat /= 0 ) RETURN
             END IF ! (p%FVM%Tree_Parms%Speedup)                
             
         ELSE              
            !.......................................................
            ! Use direct parallel computation           
            CALL BiotSavart_OpenMP(F1_new, F2_new, Pts_new, GAMMA_new, RC_new, Uind_new, P, O, xd, ErrStat,ErrMess)
         END IF   
         IF (ErrStat /= 0 ) RETURN
         
         
         K = 1
         ! Loop the points of interest. Reshape to get the output for this subroutine
         DO I_5 = 1, SIZE(Pts_all , 5)       ! Counter of point blade index     
            DO I_4 = 1, SIZE(Pts_all , 4)      ! Counter of point element index
               DO I_2 = 1, SIZE(Pts_all , 2)      ! Counter of point age index         
                  U_ind_all(1, I_2, 1, I_4, I_5) = Uind_new(K,1)
                  U_ind_all(2, I_2, 1, I_4, I_5) = Uind_new(K,2)
                  U_ind_all(3, I_2, 1, I_4, I_5) = Uind_new(K,3)
                  K = K + 1
               END DO
            END DO
         END DO                      

          
          
      ELSE  ! Non-acceleration ...............................................................................    
          
         U_ind_all  = 0.0_DbKi            
         CO = p%FVM%CO
      
          ! Loop the points of interest
         DO I_5 = 1, SIZE(Pts_all , 5)       ! Counter of point blade index     
            DO I_4 = 1, SIZE(Pts_all , 4)      ! Counter of point element index
               DO I_2 = 1, SIZE(Pts_all , 2)      ! Counter of point age index
                   
                  PX = Pts_all(1, I_2, 1, I_4, I_5)
                  PY = Pts_all(2, I_2, 1, I_4, I_5)
                  PZ = Pts_all(3, I_2, 1, I_4, I_5)          
   
                  PRESUMX = 0.0
                  PRESUMY = 0.0
                  PRESUMZ = 0.0
                    
                      
                  ! Loop the filaments  
                  Do J_5 = 1, SIZE(F1_all , 5)     ! Counter of filament blade index  
                      DO J_4 = 1, SIZE(F1_all , 4)      ! Counter of filament element index
                          DO J_2 = 1, SIZE(F1_all , 2)     ! Counter of filament age index                     
                              X1      =  F1_all(1, J_2, 1, J_4, J_5)
                              Y1      =  F1_all(2, J_2, 1, J_4, J_5)
                              Z1      =  F1_all(3, J_2, 1, J_4, J_5)
   
                              X2      =  F2_all(1, J_2, 1, J_4, J_5)
                              Y2      =  F2_all(2, J_2, 1, J_4, J_5)
                              Z2      =  F2_all(3, J_2, 1, J_4, J_5)
                             
                                                 
                              GMMA    =  GAMMA_all(1, J_2, 1, J_4, J_5)
                              RC      =  RC_all(1, J_2, 1, J_4, J_5)
                             
                             
                              X2X1    =  X2 - X1
                              Y2Y1    =  Y2 - Y1
                              Z2Z1    =  Z2 - Z1
                             
                              L       =  X2X1 * X2X1 + Y2Y1 * Y2Y1 + Z2Z1 * Z2Z1  ! Length of vortex filament (NOTE: L is L^2, as rc is rc^2)
   
                              PXX1    =  px - X1
                              pyy1    =  py - y1
                              pzz1    =  pz - z1
                              pxx2    =  px - x2
                              pyy2    =  py - y2
                              pzz2    =  pz - z2  
                             
                              R1      =  SQRT( PXX1 * PXX1 + PYY1 * PYY1 + PZZ1 * PZZ1 )
                              R2      =  SQRT( PXX2 * PXX2 + PYY2 * PYY2 + PZZ2 * PZZ2 )
                              R1DR2   =  PXX1 * PXX2 + PYY1 * PYY2 + PZZ1 * PZZ2
                              R1TR2   =  R1 * R2
                
                              IF (p%FVM%ViscFLAG) THEN
                                 ! Vatistas core model (n=2) .................................................
                                 LDR12   =  ( X2X1 * PXX1 + Y2Y1 * PYY1 + Z2Z1 * PZZ1) ** 2                              
                                 CNU     =  (( R1 * R1 ) - ( LDR12 / L ))                            
                                 CNU     =  CNU / SQRT(RC ** 2 + CNU  ** 2 )   !sliu: CNU = CNU * (RC ** 2 + CNU  ** 2 ) ** (-1/2) ! gives wrong result
                                 UBAR    =  CNU * GMMA / (TWOPI * 2) * (R1 + R2) / (R1TR2 * (R1TR2 + R1DR2))
                             
                                 ! Check infinity and NAN
                                 IF ( EqualRealNos( L, 0.0_DbKi ) )                             UBAR = 0.0_DbKi 
                                 IF ( EqualRealNos( R1 * R1 - LDR12 / L, 0.0_DbKi ) )           UBAR = 0.0_DbKi 
                                 IF ( EqualRealNos( (R1TR2 * (R1TR2 + R1DR2)) , 0.0_DbKi ) )    UBAR = 0.0_DbKi 
                                 IF ( ISNAN(UBAR) )                                             UBAR = 0.0_DbKi  
                                                              
                              ELSE ! Smoothing parameter .................................................
                                 DEN     =  R1TR2 * (R1TR2 + R1DR2) + (p%FVM%DELTA * L)                           
                                 UBAR    =  (GMMA * (R1 + R2)) / (TWOPI * 2) 
                                 UBAR    =  UBAR / DEN 
                             
                                 ! Check infinity and NAN
                                 IF ( EqualRealNos( DEN, 0.0_DbKi ) )         UBAR = 0.0_DbKi 
                                 IF ( ISNAN(UBAR) )                           UBAR = 0.0_DbKi  
                              END IF ! (p%FVM%ViscFLAG)
                          
                                                       
                              ! Influence is ignored beyond CO                             
                              IF (R1>CO)          UBAR = 0  
                              IF (R2>CO)          UBAR = 0 
                             
                              PRESUMX = PRESUMX + UBAR * (PYY1 * PZZ2 - PZZ1 * PYY2)
                              PRESUMY = PRESUMY + UBAR * (PZZ1 * PXX2 - PXX1 * PZZ2)
                              PRESUMZ = PRESUMZ + UBAR * (PXX1 * PYY2 - PYY1 * PXX2)
                             
   
                          END DO  
                      END DO 
                  END DO 
                          
                  U_ind_all(1, I_2, 1, I_4, I_5) =  PRESUMX      
                  U_ind_all(2, I_2, 1, I_4, I_5) =  PRESUMY            
                  U_ind_all(3, I_2, 1, I_4, I_5) =  PRESUMZ            
                          
               END DO  ! I_2 = 1, SIZE(Pts,2)      ! Counter of point age index
            END DO ! I_4 = 1, SIZE(Pts,4)      ! Counter of point element index               
         END DO ! I_5 = 1, SIZE(Pts,5)       ! Counter of point blade index  
         
      END IF ! (p%FVM%OpenMP_Parms%Accelerate)
      
CONTAINS
   ! ............................................................................................................................
   SUBROUTINE  Tree_Record(Uind_new, F1_new, F2_new, Pts_new, GAMMA_new, RC_new, P, O, xd, ErrStat, ErrMess)   
   
     IMPLICIT                        NONE
   
         
         ! Passing variables
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: Uind_new  ! Result by Treecode
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: F1_new
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: F2_new 
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: Pts_new
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: GAMMA_new
      REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: RC_new

      
      TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
      TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None         
   
  

      ! Local variables      
      INTEGER(IntKi)    :: Len_fila   ! # of filament
      INTEGER(IntKi)    :: Len_pts    ! # of points   
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE        :: Uind_new_test  ! For test use (compare error and CPU time)   
      INTEGER(IntKi)    :: I_1, k
      REAL(DbKi)        :: Time_1, Time_2     
      REAL(DbKi)        :: sum1, sum2      

      !! For debug
      !CHARACTER(LEN=1024)                  :: temp_number
      !REAL(DbKi), DIMENSION(100,2)         :: PRINT_NAME1
      !! For debug
            
      
      
      Len_fila   =   SIZE(F1_new , 1)
      Len_pts    =   SIZE(Pts_new, 1)    
             
      IF (.NOT. ALLOCATED( Uind_new_test ) ) ALLOCATE( Uind_new_test(Len_pts, 3))           
             
      CALL CPU_TIME(Time_1) ! Start time          
             
      CALL BiotSavart_OpenMP(F1_new, F2_new, Pts_new, GAMMA_new, RC_new, Uind_new_test, P, O, xd, ErrStat,ErrMess)
      IF (ErrStat /= 0 ) RETURN
             
      CALL CPU_TIME(Time_2) ! End time
      
      time_original                    =  Time_2 - Time_1       
      !O%FVM_Other%TIME%Time_biotsavart =  O%FVM_Other%TIME%Time_biotsavart + time_original  
      speedup_ratio                    =  time_original / time_accelerated 
                
      !K = 1
      !sum_error = 0.0_DbKi

      ! Loop the points of interest
      !DO k = 1, Len_pts
      !   DO I_1 = 1, 3 
      !      IF (.NOT. EqualRealNos( Uind_new_test(k, I_1), 0.0_DbKi ) ) THEN
      !         sum_error = sum_error + ABS(Uind_new(k, I_1) / Uind_new_test(k, I_1) - 1)       
      !      END IF
      !   END DO
      !END DO
      !Error = sum_error / Len_pts / 3    ! record the L1 norm error  
      
      K = 1
      sum1 = 0.0_DbKi
      sum2 = 0.0_DbKi

      ! Loop the points of interest to get L2 norm error. 
      DO k = 1, Len_pts
         DO I_1 = 1, 3             
            sum1 = sum1+ (Uind_new_test(k, I_1)-Uind_new(k, I_1))**2
            sum2 = sum2+ (Uind_new_test(k, I_1))**2

            !IF (.NOT. EqualRealNos( Uind_new_test(k, I_1), 0.0_DbKi ) ) THEN
               !sum_error = sum_error + ABS(Uind_new(k, I_1) / Uind_new_test(k, I_1) - 1)       
            !END IF
         END DO
      END DO
      Error = SQRT(sum1/sum2)
       
             
      CALL WRITE_Treecode(Len_fila, Len_pts, time_original,  time_accelerated, speedup_ratio,  Error, 'NONE', p)
              
             
             !! Record the error for each points at each timestep
             !IF (Len_pts > 1000) THEN
             !   PRINT_NAME1 = 0.0_DbKi
             !   WRITE (temp_number, "(I6.6)") (O%WINDS_Timestep)    ! String of the integer With heading zeros
             !          
             !   DO k=1,100
             !      PRINT_NAME1(k, 1) = Uind_new(k, 1) 
             !      PRINT_NAME1(k, 2) = Uind_new_test(k, 1) 
             !   END DO
             !   CALL SAVE_TO_TXT_2D(p, PRINT_NAME1 , 'WINDS_treecode_'// TRIM(temp_number))
             !END IF
             
      IF ( ALLOCATED( Uind_new_test ) ) DEALLOCATE( Uind_new_test)   ! Just use to check error....   
   
   
   
   END SUBROUTINE Tree_Record
   ! ............................................................................................................................
END SUBROUTINE BiotSavart     

!==================================================================================================================================
SUBROUTINE CUR_2_PREV(O, p, xd, N, ErrStat, ErrMess)
! Copy wake and vel from "current" timestep to "previous" 
! Called from: WINDS_FVM (in WINDS.f90)
! (~ cur_2_prev.m)
!....................................................................

  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState ! Other/optimization states
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep, counts from 1    
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess
      
   INTEGER(IntKi)      :: NTP
   
   NTP = p%FVM%NTP       
   
   O%FVM_Other%WAKE_DOMAIN(:,:,2:NTP,:,:)   =   O%FVM_Other%WAKE_DOMAIN(:,:,1:NTP-1,:,:)

   O%FVM_Other%VEL_DOMAIN(:,:,2:NTP,:,:)    =   O%FVM_Other%VEL_DOMAIN(:,:,1:NTP-1,:,:)
   O%FVM_Other%VEL_UIND(:,:,2:NTP,:,:)      =   O%FVM_Other%VEL_UIND(:,:,1:NTP-1,:,:)
   O%FVM_Other%VEL_UINDB(:,:,2:NTP,:,:)     =   O%FVM_Other%VEL_UINDB(:,:,1:NTP-1,:,:)      
   
 
   IF (p%FVM%Ground_Parms%GroundFLAG) THEN 
       
      SELECT CASE (p%FVM%Ground_Parms%Method)
         CASE ( 'PANEL') 
            O%FVM_Other%vel_uind_ground(:,:,2:ntp,:,:)  = O%FVM_Other%vel_uind_ground(:,:,1:ntp-1,:,:)
            O%FVM_Other%vel_uindb_ground(:,:,2:ntp,:,:) = O%FVM_Other%vel_uindb_ground(:,:,1:ntp-1,:,:)
         
         CASE ('IMAGE')
            O%FVM_Other%vel_uind_mirror(:,:,2:ntp,:,:)  = O%FVM_Other%vel_uind_mirror(:,:,1:ntp-1,:,:)
            O%FVM_Other%vel_uindb_mirror(:,:,2:ntp,:,:) = O%FVM_Other%vel_uindb_mirror(:,:,1:ntp-1,:,:)
         END SELECT
         
   END IF  
   
   

   O%FVM_Other%WAKE_RE_SHED(:,:,2:NTP,:,:)        =  O%FVM_Other%WAKE_RE_SHED(:,:,1:NTP-1,:,:)
   O%FVM_Other%WAKE_RE_TRAIL(:,:,2:NTP,:,:)       =  O%FVM_Other%WAKE_RE_TRAIL(:,:,1:NTP-1,:,:)

   O%FVM_Other%WAKE_RC_SHED(:,:,2:NTP,:,:)        =  O%FVM_Other%WAKE_RC_SHED(:,:,1:NTP-1,:,:)
   O%FVM_Other%WAKE_RC_TRAIL(:,:,2:NTP,:,:)       =  O%FVM_Other%WAKE_RC_TRAIL(:,:,1:NTP-1,:,:)

   O%FVM_Other%WAKE_LENGTH_SHED(:,:,2:NTP,:,:)    =  O%FVM_Other%WAKE_LENGTH_SHED(:,:,1:NTP-1,:,:)
   O%FVM_Other%WAKE_LENGTH_TRAIL(:,:,2:NTP,:,:)   =  O%FVM_Other%WAKE_LENGTH_TRAIL(:,:,1:NTP-1,:,:)

   O%FVM_Other%WAKE_GAMMA_SHED(:,:,2:NTP,:,:)     =  O%FVM_Other%WAKE_GAMMA_SHED(:,:,1:NTP-1,:,:)
   O%FVM_Other%WAKE_GAMMA_TRAIL(:,:,2:NTP,:,:)    =  O%FVM_Other%WAKE_GAMMA_TRAIL(:,:,1:NTP-1,:,:)
   
   O%FVM_Other%WAKE_RC_EFF_SHED(:,:,2:NTP,:,:)    =  O%FVM_Other%WAKE_RC_EFF_SHED(:,:,1:NTP-1,:,:)
   O%FVM_Other%WAKE_RC_EFF_TRAIL(:,:,2:NTP,:,:)   =  O%FVM_Other%WAKE_RC_EFF_TRAIL(:,:,1:NTP-1,:,:)




END SUBROUTINE CUR_2_PREV

!==================================================================================================================================
SUBROUTINE InducedVelocity(O, p, xd, N, CALLER, ErrStat, ErrMess)
! Calculate the induced velocity via Biot-Savart Law in different situations
! (~ InducedVelocity.m)
!....................................................................

  IMPLICIT                        NONE

      ! Passed variables      
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states   
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep 
   CHARACTER(*),                  INTENT(IN   )  :: CALLER      ! Who calls this subroutine 
   
   INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None


      ! Internal variables
   INTEGER(IntKi)      :: IDim    ! For dimensions
   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: IElement2
   INTEGER(IntKi)      :: ITimestep
   
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB  
   INTEGER(IntKi)      :: NP
   INTEGER(IntKi)      :: NTW
   INTEGER(IntKi)      :: ntroll
   
   INTEGER(IntKi)      :: J, c1, c2, K   ! counter 
   
   INTEGER(IntKi), DIMENSION(p%Element%NElm-1) :: ns_ind_1
   INTEGER(IntKi), DIMENSION(p%Element%NElm-1) :: ns_ind_2
   INTEGER(IntKi), DIMENSION(p%Element%NElm-1) :: ns_ind_3
   
   
   CHARACTER(LEN = 1024)   :: TEMP10  ! debug
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   
   NP  = p%FVM%Ground_Parms%Sqrt_Panels
      
   NTW = O%FVM_Other%NTW
   ntroll = O%FVM_Other%ntroll

   SELECT CASE  (CALLER)
   
      !..................................................................................................................... 
      CASE ('MAIN') ! Called from Main WIndS Time Series
         IF (ntroll > 0) THEN
             ! The velocity induced by shed vorticity    
            CALL BiotSavart(O%FVM_Other%WAKE_DOMAIN    (:,  1:NTW,        1:1,  1:(NST-1), :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                           O%FVM_Other%WAKE_DOMAIN     (:,  1:NTW,        1:1,  2:NST,     :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                           O%FVM_Other%WAKE_DOMAIN     (:,  2:ntroll,     1:1,  :,         :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                           O%FVM_Other%WAKE_GAMMA_SHED (:,  1:NTW,        1:1,  :,         :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                           O%FVM_Other%WAKE_RC_EFF_SHED(:,  1:NTW,        1:1,  :,         :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                           O%FVM_Other%VEL_UIND_SHED   (:,  1:(ntroll-1), 1:1,  :,         :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                           P, O, xd, ErrStat, ErrMess)
            
               ! The velocity induced by trailing vorticity 
            CALL BiotSavart(O%FVM_Other%WAKE_DOMAIN     (:,  2:NTW,        1:1,  :,   :),     &  ! F1:     Filament start          (1:3, 2:N,     1,  1:NST,  1:NB)
                            O%FVM_Other%WAKE_DOMAIN     (:,  1:(NTW-1),    1:1,  :,   :),     &  ! F2:     Filament end            (1:3, 1:N-1,   1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_DOMAIN      (:,  2:ntroll,     1:1,  :,   :),     &  ! Pts:    Points of interest      (1:3, 2:N,     1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_GAMMA_TRAIL (:,  1:(NTW-1),    1:1,  :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N-1,   1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_RC_EFF_TRAIL(:,  1:(NTW-1),    1:1,  :,   :),     &  ! RC:     Radius of filament      (1,   1:N-1,   1,  1:NST,  1:NB)
                           O%FVM_Other%VEL_UIND_TRAIL   (:,  1:(ntroll-1), 1:1,  :,   :),     &  ! Output: Induced velocity U      (1:3, 1:N-1,   1,  1:NST,  1:NB)
                           P, O, xd, ErrStat, ErrMess)
            IF (ErrStat /= 0 ) RETURN
                                
               ! Sum the induced velocity contributions due to shed and trailing filaments            
            IF (p%FVM%UindPast) THEN            
               O%FVM_Other%VEL_UIND(:, (ntroll+1):NTW, 1, :, :) = O%FVM_Other%VEL_UIND(:, ntroll:(NTW-1), 1, :, :) 
            END IF                                
            O%FVM_Other%VEL_UIND(:, 2:ntroll, 1, :, :) =  O%FVM_Other%VEL_UIND_SHED(:, 1:ntroll-1, 1, :, :) +       &
                                                         O%FVM_Other%VEL_UIND_TRAIL(:, 1:ntroll-1, 1, :, :)       
  
   
         END IF ! (ntroll > 0)
         
            ! Add the total induced velocity in the wake to the freestream velocity
         O%FVM_Other%VEL_DOMAIN(:, 2:NTW, 1, :, :) = O%FVM_Other%VEL_DOMAIN(:, 2:NTW, 1, :, :) +      &
                                                       O%FVM_Other%VEL_UIND(:, 2:NTW, 1, :, :)    
         
          
          
          ! Compute effects of ground panels.
   
         IF (p%FVM%Ground_Parms%GroundFLAG) THEN              
                         
              IF (p%FVM%Ground_Parms%Method  ==  'PANEL' ) THEN
                  
                  ! Compute induced velocity at all ground panel centers points due to wake and then calculate panel strengths
                  CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN           (:,  1:NTW,   1:1,  1:(NST-1), :),     &  ! F1:     Filament start
                                   O%FVM_Other%WAKE_DOMAIN           (:,  1:NTW,   1:1,  2:NST,     :),     &  ! F2:     Filament end
                                   O%FVM_Other%Ground%P_source       (:,  :,       1:1,  :,         :),     &  ! Pts:    Points of interest      (1:3,  1,  1,  1:(NP**2) , 1)
                                   O%FVM_Other%WAKE_GAMMA_SHED       (:,  1:NTW,   1:1,  :,         :),     &  ! GAMMA:  Circulation of filament
                                   O%FVM_Other%WAKE_RC_EFF_SHED      (:,  1:NTW,   1:1,  :,         :),     &  ! RC:     Radius of filament 
                                   O%FVM_Other%VEL_UIND_SHED_GROUND_P(:,  :,       1:1,  :,         :),     &  ! Output: Induced velocity U      (1:3,  1,  1,  1:(NP**2) , 1)
                                   P, O, xd, ErrStat, ErrMess)

                  CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN            (:,  2:NTW,   1:1,  :,    :),     &  ! F1:     Filament start
                                   O%FVM_Other%WAKE_DOMAIN            (:,  1:NTW-1, 1:1,  :,    :),     &  ! F2:     Filament end
                                   O%FVM_Other%Ground%P_source        (:,  :,       1:1,  :,    :),     &  ! Pts:    Points of interest      (1:3,  1,  1,  1:(NP**2) , 1)
                                   O%FVM_Other%WAKE_GAMMA_TRAIL       (:,  1:NTW-1, 1:1,  :,    :),     &  ! GAMMA:  Circulation of filament
                                   O%FVM_Other%WAKE_RC_EFF_TRAIL      (:,  1:NTW-1, 1:1,  :,    :),     &  ! RC:     Radius of filament     
                                   O%FVM_Other%VEL_UIND_TRAIL_GROUND_P(:,  :,       1:1,  :,    :),     &  ! Output: Induced velocity U  (1:3,  1,  1,  1:(NP**2) , 1)
                                   P, O, xd, ErrStat, ErrMess)                                   
                  IF (ErrStat /= 0 ) RETURN 
                                   
                   O%FVM_Other%VEL_UIND_GROUND_P(1,:,:,:,:)  = O%FVM_Other%VEL_UIND_SHED_GROUND_P(1,:,:,:,:)  +  &
                                           O%FVM_Other%VEL_UIND_TRAIL_GROUND_P(1,:,:,:,:)  + O%FVM_Other%WIND_INFTY(1, N, 1, 1, 1)            
                   O%FVM_Other%VEL_UIND_GROUND_P(2,:,:,:,:)  = O%FVM_Other%VEL_UIND_SHED_GROUND_P(2,:,:,:,:)  +   &
                                           O%FVM_Other%VEL_UIND_TRAIL_GROUND_P(2,:,:,:,:)  +  O%FVM_Other%WIND_INFTY(2, N, 1, 1, 1) 
                   O%FVM_Other%VEL_UIND_GROUND_P(3,:,:,:,:)  = O%FVM_Other%VEL_UIND_SHED_GROUND_P(3,:,:,:,:)  +  &
                                            O%FVM_Other%VEL_UIND_TRAIL_GROUND_P(3,:,:,:,:)  + O%FVM_Other%WIND_INFTY(3, N, 1, 1, 1)                    
                   
                   
                   DO J = 1 , NP **2
                       O%FVM_Other%Ground%V_n(J, 1) =                                                                &
                               - (O%FVM_Other%VEL_UIND_GROUND_P(1, 1 , 1, J, 1) * O%FVM_Other%Ground%n_vec(J , 1) +  &
                                  O%FVM_Other%VEL_UIND_GROUND_P(2, 1 , 1, J, 1) * O%FVM_Other%Ground%n_vec(J , 2) +  &
                                  O%FVM_Other%VEL_UIND_GROUND_P(3, 1 , 1, J, 1) * O%FVM_Other%Ground%n_vec(J , 3) ) 
                   END DO
                   
                                   
                   O%FVM_Other%Ground%Gamma(: , N) =   MATMUL(  O%FVM_Other%Ground%B ,   O%FVM_Other%Ground%V_n(:, 1) )       
                                   
                                   
                      
                   !  Add the total induced velocity in the wake to the freestream velocity
                   CALL Induced_Velocity_Ground_Panels(O, p, xd, N, 'MAIN', ErrStat, ErrMess)
                   IF (ErrStat /= 0 ) RETURN
                   O%FVM_Other%VEL_DOMAIN(:, 2:NTW, 1, :, :) = O%FVM_Other%VEL_DOMAIN(:, 2:NTW, 1, :, :) +   &
                                                               O%FVM_Other%vel_uind_ground(:,2:NTW,1,:,:)                
                                  
                                   
              ELSE IF  (p%FVM%Ground_Parms%Method  ==  'IMAGE' ) THEN
                  ! Compute induced velocity at all wake points using method of images   
                   CALL Induced_Velocity_Ground_Mirror(O, p, xd, N, 'MAIN', ErrStat, ErrMess)
                   IF (ErrStat /= 0 ) RETURN
                   O%FVM_Other%VEL_DOMAIN(:, 2:NTW, 1, :, :) = O%FVM_Other%VEL_DOMAIN(:, 2:NTW, 1, :, :) + &
                                                             O%FVM_Other%vel_uind_mirror(:,2:NTW,1,:,:)                     
                  
                  
              END IF !  (p%FVM%Ground_Parms%Method  ==  'PANEL' )
              
          END IF !  (p%FVM%Ground_Parms%GroundFLAG) 
          
          
             
     !..................................................................................................................... 
      CASE ('KJ_PRE')  ! Called From Kutta-Joukowski: Pre-iteration
         
         IF  (O%WINDS_Timestep >= p%FVM%DS_Parms%start_n .AND. p%FVM%DS_Parms%DS_Flag) THEN
             ! Neglect effects of a shed vortices on own bound vortex when using dynamic stall
             DO c1 = 1, NB
                DO c2 = 1, NS
                   DO K = 1, c2-1
                      ns_ind_1(K) = K 
                      ns_ind_2(K) = K+1
                      ns_ind_3(K) = K
                   END DO
                   DO K = c2, NS-1
                      ns_ind_1(K) = K+1
                      ns_ind_2(K) = K+2
                      ns_ind_3(K) = K+1
                   END DO
                   
                   CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN       (:,  3:NTW+1,   1:1,  ns_ind_1, c1:c1),     &  ! F1:     Filament start  
                                    O%FVM_Other%WAKE_DOMAIN       (:,  3:NTW+1,   1:1,  ns_ind_2, c1:c1),     &  ! F2:     Filament end       
                                    O%FVM_Other%POS_BOUND         (:,  N:N,       1:1,  c2:c2,    c1:c1),     &  ! Pts:    Points of interest     
                                    O%FVM_Other%WAKE_GAMMA_SHED   (:,  3:NTW+1,   1:1,  ns_ind_3, c1:c1),     &  ! GAMMA:  Circulation of filament 
                                    O%FVM_Other%WAKE_RC_EFF_SHED  (:,  3:NTW+1,   1:1,  ns_ind_3, c1:c1),     &  ! RC:     Radius of filament     
                                    O%FVM_Other%VEL_UINDB_SHED_PRE(:,  1:1,       1:1,  c2:c2,    c1:c1),     &  ! Output: Induced velocity U     
                                    P, O, xd, ErrStat, ErrMess)
                END DO
             END DO

         ELSE ! No dynamic stall, induced velocity calc includes all shed vortices
               ! The velocity induced by shed vorticity   
            CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN       (:,  3:NTW+1,   1:1,   1:(NST-1), :),     &  ! F1:     Filament start           (1:3,  3:N+1,   1,  1:NS,  1:NB)
                             O%FVM_Other%WAKE_DOMAIN       (:,  3:NTW+1,   1:1,   2:NST,     :),     &  ! F2:     Filament end             (1:3,  3:N+1,   1,  2:NST, 1:NB)
                             O%FVM_Other%POS_BOUND         (:,  N:N,       1:1,   :,         :),     &  ! Pts:    Points of interest       (1:3,  N:N,     1,  1:NS,  1:NB)
                             O%FVM_Other%WAKE_GAMMA_SHED   (:,  3:NTW+1,   1:1,   :,         :),     &  ! GAMMA:  Circulation of filament  (1,    3:N+1,   1,  1:NS,  1:NB)
                             O%FVM_Other%WAKE_RC_EFF_SHED  (:,  3:NTW+1,   1:1,   :,         :),     &  ! RC:     Radius of filament       (1,    3:N+1,   1,  1:NS,  1:NB)
                             O%FVM_Other%VEL_UINDB_SHED_PRE(:,  1:1,       1:1,   :,         :),     &  ! Output: Induced velocity U       (1:3,  1:1,     1,  1:NS,  1:NB)
                             P, O, xd, ErrStat, ErrMess)
                             
         END IF ! (O%WINDS_Timestep >= p%FVM%DS_Parms%start_n .AND. p%FVM%DS_Parms%DS_Flag)     
    
            ! The velocity induced by trailing vorticity 
         CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN        (:,  3:NTW+1,  1:1,    :,     :),     &  ! F1:     Filament start              (1:3,  3:N+1,   1,  1:NST,  1:NB)
                          O%FVM_Other%WAKE_DOMAIN        (:,  2:NTW,    1:1,    :,     :),     &  ! F2:     Filament end                (1:3,  2:N,     1,  1:NST,  1:NB)
                          O%FVM_Other%POS_BOUND          (:,  N:N,      1:1,    :,     :),     &  ! Pts:    Points of interest          (1:3,  N:N,     1,  1:NS,   1:NB)
                          O%FVM_Other%WAKE_GAMMA_TRAIL   (:,  2:NTW,    1:1,    :,     :),     &  ! GAMMA:  Circulation of filament     (1,    2:N,     1,  1:NST,  1:NB)
                          O%FVM_Other%WAKE_RC_EFF_TRAIL  (:,  2:NTW,    1:1,    :,     :),     &  ! RC:     Radius of filament          (1,    2:N,     1,  1:NST,  1:NB)
                          O%FVM_Other%VEL_UINDB_TRAIL_PRE(:,  1:1,      1:1,    :,     :),     &  ! Output: Induced velocity U          (1:3,  1:1,     1,  1:NS,   1:NB)
                          P, O, xd, ErrStat, ErrMess)
         !IF (ErrStat /= 0 ) RETURN 
                           
            ! Compute tower and ground effects outside while loop as they do not depend on updated circulation                
         IF (p%FVM%Twr_Parms%TWRFLAG) THEN
             CALL TOWER_SHADOW(O, p, xd, N, ErrStat, ErrMess)    
             !IF (ErrStat /= 0 ) RETURN
         END IF
           
         IF (p%FVM%Ground_Parms%GroundFLAG) THEN              
               
             IF (p%FVM%Ground_Parms%Method  ==  'PANEL' ) THEN  
                  ! Compute induced velocity at all wake points due to ground panels        
                CALL Induced_Velocity_Ground_Panels(O, p, xd, N, 'KJ', ErrStat, ErrMess)
                IF (ErrStat /= 0 ) RETURN
              
             ELSE IF  (p%FVM%Ground_Parms%Method  ==  'IMAGE' ) THEN
                  !Compute induced velocity at all wake points using method of images
                CALL Induced_Velocity_Ground_Mirror(O, p, xd, N, 'KJ', ErrStat, ErrMess)
                IF (ErrStat /= 0 ) RETURN
               
             END IF !  (p%FVM%Ground_Parms%Method  ==  'PANEL' )
         
         END IF !  (p%FVM%Ground_Parms%GroundFLAG) 
                        
              
              
   
     !..................................................................................................................... 
      CASE ('KJ_ITER')   ! Called from Kutta-Joukowski: Main Iteration  
          
         IF  (O%WINDS_Timestep >= p%FVM%DS_Parms%start_n .AND. p%FVM%DS_Parms%DS_Flag) THEN
             ! Neglect effects of a shed vortices on own bound vortex when using dynamic stall
             DO c1 = 1, NB
                DO c2 = 1, NS
                   DO K = 1, c2-1
                      ns_ind_1(K) = K 
                      ns_ind_2(K) = K+1
                      ns_ind_3(K) = K
                   END DO
                   DO K = c2, NS-1
                      ns_ind_1(K) = K+1
                      ns_ind_2(K) = K+2
                      ns_ind_3(K) = K+1
                   END DO
                   
                   CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN       (:,  1:2,   1:1,  ns_ind_1, c1:c1),     &  ! F1:     Filament start          
                                    O%FVM_Other%WAKE_DOMAIN       (:,  1:2,   1:1,  ns_ind_2, c1:c1),     &  ! F2:     Filament end             
                                    O%FVM_Other%POS_BOUND         (:,  N:N,   1:1,  c2:c2,    c1:c1),     &  ! Pts:    Points of interest        
                                    O%FVM_Other%WAKE_GAMMA_SHED   (:,  1:2,   1:1,  ns_ind_3, c1:c1),     &  ! GAMMA:  Circulation of filament  
                                    O%FVM_Other%WAKE_RC_EFF_SHED  (:,  1:2,   1:1,  ns_ind_3, c1:c1),     &  ! RC:     Radius of filament       
                                    O%FVM_Other%VEL_UINDB_SHED    (:,  1:1,   1:1,  c2:c2,    c1:c1),     &  ! Output: Induced velocity U     
                                    P, O, xd, ErrStat, ErrMess)
                END DO
             END DO          
          
          ELSE ! No dynamic stall, induced velocity calc includes all shed vortices
                 ! The velocity induced by shed vorticity   
             CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN       (:,  1:2,   1:1,  1:NST-1, :),     &  ! F1:     Filament start            (1:3,  1:2,   1,  1:NS,  1:NB)
                              O%FVM_Other%WAKE_DOMAIN       (:,  1:2,   1:1,  2:NST,   :),     &  ! F2:     Filament end              (1:3,  1:2,   1,  2:NST, 1:NB)
                              O%FVM_Other%POS_BOUND         (:,  N:N,   1:1,  :,       :),     &  ! Pts:    Points of interest        (1:3,  N:N,   1,  1:NST, 1:NB)
                              O%FVM_Other%WAKE_GAMMA_SHED   (:,  1:2,   1:1,  :,       :),     &  ! GAMMA:  Circulation of filament   (1,    1:2,   1,  1:NS,  1:NB)
                              O%FVM_Other%WAKE_RC_EFF_SHED  (:,  1:2,   1:1,  :,       :),     &  ! RC:     Radius of filament        (1,    1:2,   1,  1:NS,  1:NB)
                              O%FVM_Other%VEL_UINDB_SHED    (:,  1:1,   1:1,  :,       :),     &  ! Output: Induced velocity U        (1:3,  1:1,   1,  1:NST, 1:NB)
                              P, O, xd, ErrStat, ErrMess)
          END IF ! (O%WINDS_Timestep >= p%FVM%DS_Parms%start_n .AND. p%FVM%DS_Parms%DS_Flag)     
                           
               ! The velocity induced by trailing vorticity 
          CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN        (:,  2:2,  1:1,   :,   :),     &  ! F1:     Filament start           (1:3,  2:2,   1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_DOMAIN        (:,  1:1,  1:1,   :,   :),     &  ! F2:     Filament end             (1:3,  1:1,   1,  1:NST,  1:NB)
                           O%FVM_Other%POS_BOUND          (:,  N:N,  1:1,   :,   :),     &  ! Pts:    Points of interest       (1:3,  N:N,   1,  1:NS,   1:NB)
                           O%FVM_Other%WAKE_GAMMA_TRAIL   (:,  1:1,  1:1,   :,   :),     &  ! GAMMA:  Circulation of filament  (1,    1:2,   1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_RC_EFF_TRAIL  (:,  1:1,  1:1,   :,   :),     &  ! RC:     Radius of filament       (1,    1:2,   1,  1:NST,  1:NB)
                           O%FVM_Other%VEL_UINDB_TRAIL    (:,  1:1,  1:1,   :,   :),     &  ! Output: Induced velocity U       (1:3,  1:1,   1,  1:NS,   1:NB)
                           P, O, xd, ErrStat, ErrMess)    
          IF (ErrStat /= 0 ) RETURN
          
          O%FVM_Other%VEL_UINDB(:,N,1,:,:) = (O%FVM_Other%VEL_UINDB_SHED_PRE(:,1,1,:,:) + O%FVM_Other%VEL_UINDB_SHED(:,1,1,:,:)) &
                                            + (O%FVM_Other%VEL_UINDB_TRAIL_PRE(:,1,1,:,:) + O%FVM_Other%VEL_UINDB_TRAIL(:,1,1,:,:)) 
                    
          
          IF (p%FVM%Twr_Parms%TWRFLAG) THEN   
              ! Add in tower shadow uind
             O%FVM_Other%VEL_UINDB(:,N,1,:,:) = O%FVM_Other%VEL_UINDB(:,N,1,:,:) + O%FVM_Other%VEL_UIND_TOWER(:,N,1,:,:)
              
          END IF ! (p%FVM%Twr_Parms%TWRFLAG)
          
          
          IF (p%FVM%Ground_Parms%GroundFLAG) THEN              
                         
              IF (p%FVM%Ground_Parms%Method  ==  'PANEL' ) THEN  
                  O%FVM_Other%VEL_UINDB(:,N,1,:,:)  =  O%FVM_Other%VEL_UINDB(:,N,1,:,:)  +  O%FVM_Other%VEL_UINDB_GROUND(:,N,1,:,:)
                  
              !ELSE IF  (p%FVM%Ground_Parms%Method  ==  'IMAGE' )  ! sliu: to be written later
                  
              END IF !  (p%FVM%Ground_Parms%Method  ==  'PANEL' )
              
          END IF !  (p%FVM%Ground_Parms%GroundFLAG)           
          
          
          
      
     !..................................................................................................................... 
      CASE ('RK2','RK4')  
         IF (ntroll > 0) THEN
             ! The velocity induced by shed vorticity    
             CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN  (:,  1:NTW+1,    1:1,  1:NST-1,  :),     &  ! F1:     Filament start         (1:3, 1:N+1,     1,  1:(NST-1), 1:NB)
                           O%FVM_Other%WAKE_DOMAIN     (:,  1:NTW+1,    1:1,  2:NST,    :),     &  ! F2:     Filament end           (1:3, 1:N+1,     1,  2:NST,     1:NB)
                           O%FVM_Other%WAKE_DOMAIN     (:,  3:ntroll+1, 1:1,  :,        :),     &  ! Pts:    Points of interest     (1:3, 3:N+1,     1,  1:NST,     1:NB)
                           O%FVM_Other%WAKE_GAMMA_SHED (:,  1:NTW+1,    1:1,  :,        :),     &  ! GAMMA:  Circulation of filament(1,   1:N+1,     1,  1:(NST-1), 1:NB)
                           O%FVM_Other%WAKE_RC_EFF_SHED(:,  1:NTW+1,    1:1,  :,        :),     &  ! RC:     Radius of filament     (1,   1:N+1,     1,  1:(NST-1), 1:NB)
                           O%FVM_Other%VEL_UIND_SHED   (:,  1:ntroll-1, 1:1,  :,        :),     &  ! Output: Induced velocity U     (1:3, 1:N-1,     1,  1:NST,     1:NB)
                           P, O, xd, ErrStat, ErrMess)
            
                 ! The velocity induced by trailing vorticity 
             CALL  BiotSavart(O%FVM_Other%WAKE_DOMAIN   (:,  2:NTW+1,      1:1,   :,   :),     &  ! F1:     Filament start         (1:3, 2:N+1,  1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_DOMAIN      (:,  1:NTW,        1:1,   :,   :),     &  ! F2:     Filament end           (1:3, 1:N,    1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_DOMAIN      (:,  3:ntroll+1,   1:1,   :,   :),     &  ! Pts:    Points of interest     (1:3, 3:N+1,  1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_GAMMA_TRAIL (:,  1:NTW,        1:1,   :,   :),     &  ! GAMMA:  Circulation of filament(1,   1:N,    1,  1:NST,  1:NB)
                           O%FVM_Other%WAKE_RC_EFF_TRAIL(:,  1:NTW,        1:1,   :,   :),     &  ! RC:     Radius of filament     (1,   1:N,    1,  1:NST,  1:NB)
                           O%FVM_Other%VEL_UIND_TRAIL   (:,  1:(ntroll-1), 1:1,   :,   :),     &  ! Output: Induced velocity U     (1:3, 1:N-1,  1,  1:NST,  1:NB)
                           P, O, xd, ErrStat, ErrMess)
            ! IF (ErrStat /= 0 ) RETURN      
         
               ! Sum the induced velocity contributions due to shed and trailing filaments
            IF (p%FVM%UindPast) THEN            
               O%FVM_Other%VEL_UIND_RK(:, (ntroll+2):NTW+1, 1, :, :) = O%FVM_Other%VEL_UIND_RK(:, (ntroll+1):NTW, 1, :, :) 
            END IF                             
             O%FVM_Other%VEL_UIND_RK(:, 3:(ntroll+1), O%FVM_Other%RK_counter, :, :) =                                  &
                                                            O%FVM_Other%VEL_UIND_SHED(:, 1:(ntroll-1), 1, :, :) +     & 
                                                            O%FVM_Other%VEL_UIND_TRAIL(:, 1:(ntroll-1), 1, :, :)     
            
         END IF    
         
         O%FVM_Other%VEL_DOMAIN_RK(:, 3:(NTW+1), O%FVM_Other%RK_counter, :, :) =                                     &
                                              O%FVM_Other%VEL_DOMAIN_RK(:, 3:(NTW+1), O%FVM_Other%RK_counter, :, :) +   &
                                              O%FVM_Other%VEL_UIND_RK(:, 3:(NTW+1),  O%FVM_Other%RK_counter, :, :)
            
            
      END SELECT
 
    

END SUBROUTINE InducedVelocity

!==================================================================================================================================
SUBROUTINE Induced_Velocity_Ground_Mirror(O, p, xd, N, CALLER, ErrStat, ErrMess)
! calculates the induced velocity of ground effects (method of images) via Biot-Savart Law.
!..........................................................................

  IMPLICIT                        NONE

      ! Passed variables      
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states   
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep 
   CHARACTER(*),                  INTENT(IN   )  :: CALLER      ! Who calls this subroutine 
   
   INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None


      ! Internal variables
   INTEGER(IntKi)      :: IDim    ! For dimensions
   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: IElement2
   INTEGER(IntKi)      :: ITimestep
   
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB 
   
   INTEGER(IntKi)      :: xdim
   INTEGER(IntKi)      :: ydim
   INTEGER(IntKi)      :: zdim
   INTEGER(IntKi)      :: I
   INTEGER(IntKi)      :: J
   
   INTEGER(IntKi)      :: Np
   
   INTEGER(IntKi)      :: steps   ! stored steps #?
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   
   Np = p%FVM%Ground_Parms%Sqrt_Panels
   
   steps = 1
   xdim = Np
   ydim = Np   
   zdim = 1
   
   
   
   O%FVM_Other%Ground%MIRROR = 1.0
   
   ! Set the Z component negative to produce mirrored wake
   O%FVM_Other%Ground%MIRROR(3,:,:,:,:)  = O%FVM_Other%Ground%MIRROR(3,:,:,:,:) ** (-1.0)
   
   DO IBlade = 1, NB
      DO IElement2 = 1, NST
         DO ITimestep = 1, N
            DO IDim = 1, 3
               O%FVM_Other%Ground%WAKE_MIRROR(1, ITimestep, 1, IElement, IBlade) =                        &
                                           O%FVM_Other%WAKE_DOMAIN(1, ITimestep, 1, IElement, IBlade) *   &
                                           O%FVM_Other%Ground%MIRROR(1, ITimestep, 1, IElement, IBlade)
            END DO
         END DO
      END DO
   END DO
   
   O%FVM_Other%Ground%GAMMA_SHED_MIRROR(:,1:N,1,:,:)  = O%FVM_Other%WAKE_GAMMA_SHED(:,1:N,1,:,:) ** (-1.0)
   O%FVM_Other%Ground%GAMMA_TRAIL_MIRROR(:,1:N,1,:,:)  = O%FVM_Other%WAKE_GAMMA_TRAIL(:,1:N-1,1,:,:) ** (-1.0)
   
   
   
   
   
   IF (CALLER == 'MAIN') THEN
       CALL BiotSavart(O%FVM_Other%Ground%WAKE_MIRROR      (:,   1:N,    1:1,  1:NS,   :),     &  ! F1:     Filament start        
                       O%FVM_Other%Ground%WAKE_MIRROR      (:,   1:N,    1:1,  2:NST,  :),     &  ! F2:     Filament end            
                       O%FVM_Other%WAKE_DOMAIN             (:,   2:N,    1:1,  :,     :),     &  ! Pts:    Points of interest     
                       O%FVM_Other%Ground%GAMMA_SHED_MIRROR(:,   :,      :,    :,     :),     &  ! GAMMA:  Circulation of filament 
                       O%FVM_Other%WAKE_RC_EFF_SHED        (:,   :,      :,    :,     :),     &  ! RC:     Radius of filament      
                       O%FVM_Other%Ground%Vind_mirtmp_SHED (:,   1:N-1,  :,    :,     :),     &  ! Output: Induced velocity U     
                       P, O, xd, ErrStat, ErrMess)
                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%WAKE_MIRROR       (:,   2:N,    1:1,  :,   :),     &  ! F1:     Filament start        
                       O%FVM_Other%Ground%WAKE_MIRROR       (:,   1:N-1,  1:1,  :,  :),     &  ! F2:     Filament end            
                       O%FVM_Other%WAKE_DOMAIN              (:,   2:N,    1:1,  :,     :),     &  ! Pts:    Points of interest     
                       O%FVM_Other%Ground%GAMMA_TRAIL_MIRROR(:,   :,      :,    :,     :),     &  ! GAMMA:  Circulation of filament 
                       O%FVM_Other%WAKE_RC_EFF_TRAIL        (:,   :,      :,    :,     :),     &  ! RC:     Radius of filament      
                       O%FVM_Other%Ground%Vind_mirtmp_TRAIL (:,   1:N-1,  :,    :,     :),     &  ! Output: Induced velocity U     
                       P, O, xd, ErrStat, ErrMess)
       IF (ErrStat /= 0 ) RETURN
                       
       ! Sum the induced velocity contributions due to shed and trailing filaments
       O%FVM_Other%VEL_UIND_MIRROR(1:3, 2:N, 1, 1:NST, 1:3) = O%FVM_Other%Ground%Vind_mirtmp_SHED(1:3, 1:N-1, 1, 1:NST, 1:3)  +   &
                                                              O%FVM_Other%Ground%Vind_mirtmp_TRAIL(1:3, 1:N-1, 1, 1:NST, 1:3)
   
   
              
                       
   ELSE  IF (CALLER == 'KJ') THEN
       CALL BiotSavart(O%FVM_Other%Ground%WAKE_MIRROR      (:,   1:N,    1:1,    1:NS,   :),     &  ! F1:     Filament start        
                       O%FVM_Other%Ground%WAKE_MIRROR      (:,   1:N,    1:1,    2:NST,  :),     &  ! F2:     Filament end            
                       O%FVM_Other%POS_BOUND               (:,   N:N,    1:1,  :,     :),     &  ! Pts:    Points of interest     
                       O%FVM_Other%Ground%GAMMA_SHED_MIRROR(:,   :,      :,    :,     :),     &  ! GAMMA:  Circulation of filament 
                       O%FVM_Other%WAKE_RC_EFF_SHED        (:,   :,      :,    :,     :),     &  ! RC:     Radius of filament      
                       O%FVM_Other%Ground%VindB_mirtmp_SHED(:,   1:1,    :,    :,     :),     &  ! Output: Induced velocity U     
                       P, O, xd, ErrStat, ErrMess)
                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%WAKE_MIRROR       (:,   2:N,    1:1,    :,     :),     &  ! F1:     Filament start        
                       O%FVM_Other%Ground%WAKE_MIRROR       (:,   1:N-1,  1:1,    :,     :),     &  ! F2:     Filament end            
                       O%FVM_Other%POS_BOUND                (:,   N:N,    1:1,    :,     :),     &  ! Pts:    Points of interest     
                       O%FVM_Other%Ground%GAMMA_TRAIL_MIRROR(:,   :,      :,      :,     :),     &  ! GAMMA:  Circulation of filament 
                       O%FVM_Other%WAKE_RC_EFF_TRAIL        (:,   :,      :,      :,     :),     &  ! RC:     Radius of filament      
                       O%FVM_Other%Ground%VindB_mirtmp_TRAIL(:,   1:1,    :,      :,     :),     &  ! Output: Induced velocity U     
                       P, O, xd, ErrStat, ErrMess)
       IF (ErrStat /= 0 ) RETURN
                       
       ! Sum the induced velocity contributions due to shed and trailing filaments
       O%FVM_Other%VEL_UINDB_MIRROR(1:3, 1, 1, 1:NST, 1:3) = O%FVM_Other%Ground%VindB_mirtmp_SHED(1:3, 1, 1, 1:NST, 1:3) +  &
                                                               O%FVM_Other%Ground%VindB_mirtmp_TRAIL(1:3, 1, 1, 1:NST, 1:3)         
   
                      
   END IF ! (CALLER == 'MAIN')    
   
             

END SUBROUTINE Induced_Velocity_Ground_Mirror


!==================================================================================================================================
SUBROUTINE Induced_Velocity_Ground_Panels(O, p, xd, N, CALLER, ErrStat, ErrMess)
! calculates the induced velocity of ground effects (panels method) via Biot-Savart Law.
!..........................................................................

  IMPLICIT                        NONE

      ! Passed variables      
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states   
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep 
   CHARACTER(*),                  INTENT(IN   )  :: CALLER      ! Who calls this subroutine 
   
   INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None


      ! Internal variables
   INTEGER(IntKi)      :: IDim    ! For dimensions
   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: IElement2
   INTEGER(IntKi)      :: ITimestep
   
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB 
   
   INTEGER(IntKi)      :: xdim
   INTEGER(IntKi)      :: ydim
   INTEGER(IntKi)      :: zdim
   INTEGER(IntKi)      :: I
   INTEGER(IntKi)      :: J
   
   INTEGER(IntKi)      :: Np
   
   INTEGER(IntKi)      :: steps   ! stored steps #?
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   
   Np = p%FVM%Ground_Parms%Sqrt_Panels
   
   steps = 1
   xdim = Np
   ydim = Np   
   zdim = 1
   
   DO I = 1, ydim
      DO J = 1, xdim
         O%FVM_Other%Ground%grid_lefttop(1,i,1,j,1)     =  O%FVM_Other%Ground%xp_mesh(i,j,1)
         O%FVM_Other%Ground%grid_lefttop(2,i,1,j,1)     =  O%FVM_Other%Ground%yp_mesh(i,j,1)   
         O%FVM_Other%Ground%grid_lefttop(3,i,1,j,1)     =  O%FVM_Other%Ground%zp_mesh(i,j,1)   
         
         O%FVM_Other%Ground%grid_righttop(1,i,1,j,1)    =  O%FVM_Other%Ground%xp_mesh(i+1,j,1)
         O%FVM_Other%Ground%grid_righttop(2,i,1,j,1)    =  O%FVM_Other%Ground%yp_mesh(i+1,j,1)   
         O%FVM_Other%Ground%grid_righttop(3,i,1,j,1)    =  O%FVM_Other%Ground%zp_mesh(i+1,j,1)            
         
         O%FVM_Other%Ground%grid_rightbottom(1,i,1,j,1) =  O%FVM_Other%Ground%xp_mesh(i+1,j+1,1)
         O%FVM_Other%Ground%grid_rightbottom(2,i,1,j,1) =  O%FVM_Other%Ground%yp_mesh(i+1,j+1,1)   
         O%FVM_Other%Ground%grid_rightbottom(3,i,1,j,1) =  O%FVM_Other%Ground%zp_mesh(i+1,j+1,1)      
         
         O%FVM_Other%Ground%grid_leftbottom(1,i,1,j,1)  =  O%FVM_Other%Ground%xp_mesh(i,j+1,1)
         O%FVM_Other%Ground%grid_leftbottom(2,i,1,j,1)  =  O%FVM_Other%Ground%yp_mesh(i,j+1,1)   
         O%FVM_Other%Ground%grid_leftbottom(3,i,1,j,1)  =  O%FVM_Other%Ground%zp_mesh(i,j+1,1)             
         
         O%FVM_Other%Ground%gamma_grid(1,i,1,j,1)       =  O%FVM_Other%Ground%Gamma(j + (i-1) * xdim, N)
         O%FVM_Other%Ground%rc_grid   (1,i,1,j,1)       =  0.00
         
      END DO
   END DO
   
   IF (CALLER == 'MAIN') THEN
       CALL BiotSavart(O%FVM_Other%Ground%grid_righttop (:,   :,      :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_lefttop  (:,   :,      :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%WAKE_DOMAIN          (:,   2:N,    1:1,  :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid    (:,   :,      :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid       (:,   :,      :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_top_main(:,   1:N-1,  :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)
                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%grid_rightbottom(:,   :,     :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_righttop   (:,   :,     :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%WAKE_DOMAIN            (:,   2:N,   1:1,  :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid      (:,   :,     :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid         (:,   :,     :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_right_main(:,   1:N-1, :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%grid_leftbottom (:,   :,     :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_rightbottom(:,   :,     :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%WAKE_DOMAIN            (:,   2:N,   1:1,  :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid      (:,   :,     :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid         (:,   :,     :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_bottom_main(:,  1:N-1, :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)                               
                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%grid_lefttop   (:,   :,     :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_leftbottom(:,   :,     :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%WAKE_DOMAIN           (:,   2:N,   1:1,  :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid     (:,   :,     :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid        (:,   :,     :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_left_main(:,   1:N-1, :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)                                
       IF (ErrStat /= 0 ) RETURN
       
       ! Sum the induced velocities
       O%FVM_Other%vel_uind_ground(1:3, 2:N, 1, 1:NST, 1:3) = O%FVM_Other%Ground%V_ind_top_main   (1:3, 1:N-1, 1, 1:NST, 1:3) +   &
                                                              O%FVM_Other%Ground%V_ind_bottom_main(1:3, 1:N-1, 1, 1:NST, 1:3) +   &
                                                              O%FVM_Other%Ground%V_ind_left_main  (1:3, 1:N-1, 1, 1:NST, 1:3) +   &
                                                              O%FVM_Other%Ground%V_ind_right_main (1:3, 1:N-1, 1, 1:NST, 1:3)
   
              
                       
   ELSE IF (CALLER == 'KJ')  THEN
       
       CALL BiotSavart(O%FVM_Other%Ground%grid_righttop(:,   :,    :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_lefttop (:,   :,    :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%POS_BOUND           (:,   N:N,  :,    :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid   (:,   :,    :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid      (:,   :,    :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_top_KJ (:,   :,    :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)
                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%grid_rightbottom(:,   :,    :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_righttop   (:,   :,    :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%POS_BOUND              (:,   N:N,  :,    :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid      (:,   :,    :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid         (:,   :,    :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_right_KJ  (:,   :,    :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%grid_leftbottom (:,   :,    :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_rightbottom(:,   :,    :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%POS_BOUND              (:,   N:N,  :,    :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid      (:,   :,    :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid         (:,   :,    :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_bottom_KJ (:,   :,    :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)                               
                       
                       
       CALL BiotSavart(O%FVM_Other%Ground%grid_lefttop   (:,   :,    :,    :,   :),     &  ! F1:     Filament start          (1:3, 1:N,     1,  1:(NST-1), 1:NB)
                       O%FVM_Other%Ground%grid_leftbottom(:,   :,    :,    :,   :),     &  ! F2:     Filament end            (1:3, 1:N,     1,  2:NST,     1:NB)
                       O%FVM_Other%POS_BOUND             (:,   N:N,  :,    :,   :),     &  ! Pts:    Points of interest      (1,   2:N,     1,  1:NST,     1:NB)
                       O%FVM_Other%Ground%gamma_grid     (:,   :,    :,    :,   :),     &  ! GAMMA:  Circulation of filament (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%rc_grid        (:,   :,    :,    :,   :),     &  ! RC:     Radius of filament      (1,   1:N,     1,  1:NST-1,   1:NB)
                       O%FVM_Other%Ground%V_ind_left_KJ  (:,   :,    :,    :,   :),     &  ! Output: Induced velocity U      (1,   1:(N-1), 1,  1:NST,     1:NB)
                       P, O, xd, ErrStat, ErrMess)     
       IF (ErrStat /= 0 ) RETURN                
                       
       ! Sum the induced velocities
       O%FVM_Other%vel_uindb_ground(1:3, 1, 1, 1:NS, 1:3) =  O%FVM_Other%Ground%V_ind_top_KJ   (1:3, 1, 1, 1:NS, 1:3)   +  &
                                                             O%FVM_Other%Ground%V_ind_bottom_KJ(1:3, 1, 1, 1:NS, 1:3) +  &
                                                             O%FVM_Other%Ground%V_ind_left_KJ  (1:3, 1, 1, 1:NS, 1:3)   +  &
                                                             O%FVM_Other%Ground%V_ind_right_KJ (1:3, 1, 1, 1:NS, 1:3)
   
                      
   END IF ! (CALLER == 'MAIN')    
   
       
END SUBROUTINE Induced_Velocity_Ground_Panels

!==================================================================================================================================
SUBROUTINE KuttaJoukowski(O, p, xd, N, ErrStat, ErrMess)
! Computes the bound vortex filament strength via Kutta-Joukowski theorem, solving via fixed-point iteration.
! (~ KuttaJoukowski.m)
!....................................................................

      ! Passed variables      
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep 
   INTEGER,                       INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess
   
      
      ! Internal variables
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS  
   INTEGER(IntKi)      :: NT
   INTEGER(IntKi)      :: NB 
   
   REAL(DbKi)          :: DG_MAX
   INTEGER(IntKi)      :: Iiter     ! The counter of iteration
   REAL(DbKi),DIMENSION(p%FVM%MAXITER) :: converge
   REAL(DbKi)          :: AVE_Conv
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
      

   Iiter   =  1_IntKi
   DG_MAX  =  1.0_DbKi
   O%FVM_Other%KJ%DG     =  1.0_DbKi
   O%FVM_Other%KJ%Relax  =  p%FVM%Relax      ! If not copy p%FVM%Relax into O%FVM_Other%KJ%Relax, it should be "INTENT(INOUT)  :: p", which conflicts with SUBROUTINE AD_CalcOutput.
   
      ! .............................................................................
      ! Compute velocity due to everything but newest shed and trailed vortices before while loop
   CALL InducedVelocity(O, p, xd, N, 'KJ_PRE', ErrStat, ErrMess)
   IF (ErrStat /= 0 ) RETURN
   
      ! .............................................................................
      !Update dynamic stall variables for new time step
   IF (p%FVM%DS_Parms%DS_Flag) THEN
      IF (O%WINDS_Timestep > p%FVM%DS_Parms%start_n ) THEN
          O%FVM_Other%DS%global_previous_alpha  = O%FVM_Other%DS%global_current_alpha 
          O%FVM_Other%DS%global_previous_sigma1 = O%FVM_Other%DS%global_current_sigma1
          O%FVM_Other%DS%global_previous_sigma3 = O%FVM_Other%DS%global_current_sigma3 
          
          O%FVM_Other%DS%attached_previous_q    = O%FVM_Other%DS%attached_current_q   
          O%FVM_Other%DS%attached_previous_X1   = O%FVM_Other%DS%attached_current_X1   
          O%FVM_Other%DS%attached_previous_X2   = O%FVM_Other%DS%attached_current_X2  
          O%FVM_Other%DS%attached_previous_K_alpha   = O%FVM_Other%DS%attached_current_K_alpha  
          O%FVM_Other%DS%attached_previous_dK_alpha  = O%FVM_Other%DS%attached_current_dK_alpha 
          O%FVM_Other%DS%attached_previous_K_q   =  O%FVM_Other%DS%attached_current_K_q  
          O%FVM_Other%DS%attached_previous_dK_q  =  O%FVM_Other%DS%attached_current_dK_q  
          O%FVM_Other%DS%attached_previous_ddK_q  =  O%FVM_Other%DS%attached_current_ddK_q           
          O%FVM_Other%DS%attached_previous_C_p_n =  O%FVM_Other%DS%attached_current_C_p_n 
          
          O%FVM_Other%DS%TEsep_previous_D_p      =  O%FVM_Other%DS%TEsep_current_D_p  
          O%FVM_Other%DS%TEsep_previous_D_f      =  O%FVM_Other%DS%TEsep_current_D_f    
          O%FVM_Other%DS%TEsep_previous_f_prime  =  O%FVM_Other%DS%TEsep_current_f_prime  
          O%FVM_Other%DS%TEsep_previous_f_2prime =  O%FVM_Other%DS%TEsep_current_f_2prime 
          O%FVM_Other%DS%LEsep_previous_tau_v    =  O%FVM_Other%DS%LEsep_current_tau_v   
          O%FVM_Other%DS%LEsep_previous_C_v_n    =  O%FVM_Other%DS%LEsep_current_C_v_n  
          O%FVM_Other%DS%LEsep_previous_C_v      =  O%FVM_Other%DS%LEsep_current_C_v   
          
          O%FVM_Other%DS%global_ca_mdl_change      = 0_IntKi
          O%FVM_Other%DS%global_ca_mdl             = 0_IntKi
          O%FVM_Other%DS%global_ca_mdl_over_write  = 0_IntKi
      END IF
   END IF
   
   
   ! .............................................................................
   ! Main Kutta-Joukowski iteration     
   IF ( p%FVM%DS_Parms%relax_tune  .AND. O%WINDS_Timestep <=2 ) THEN
       
       converge = 0.0_DbKi
       converge(1) = 1.0_DbKi
       
       DO WHILE (  converge(Iiter) > p%FVM%TOL  .AND.  (Iiter < p%FVM%MAXITER)  )
          O%FVM_Other%KJ%GAMMA1(1, 1, 1, 1:NS, 1:NB) = O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, 1:NS, 1:NB)
          CALL KJ(O, p, DG_MAX, Iiter, ErrStat, ErrMess)          
          converge(Iiter+1) = DG_MAX
          
          IF (Iiter>=5) THEN
              AVE_Conv = SUM(converge(Iiter-2: Iiter))/3
              IF (converge(Iiter+1) > AVE_Conv  .AND.  O%FVM_Other%KJ%Relax > 0.01) THEN              
                 O%FVM_Other%KJ%Relax = O%FVM_Other%KJ%Relax * 0.98 
              END IF
          END IF
          
          IF (ErrStat /= 0 ) RETURN
  
          Iiter = Iiter + 1 
       END DO
       
   ELSE
 
      DO WHILE (  DG_MAX > p%FVM%TOL  .AND.  (Iiter < p%FVM%MAXITER)  )
       
         O%FVM_Other%KJ%GAMMA1(1, 1, 1, 1:NS, 1:NB) = O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, 1:NS, 1:NB)

         CALL KJ(O, p, DG_MAX, Iiter, ErrStat, ErrMess)
         IF (ErrStat /= 0 ) RETURN

         Iiter = Iiter + 1
      
         IF (Iiter > p%FVM%MAXITER) THEN 
             ErrStat = ErrID_Fatal
             ErrMess = ' Error (in WInDS): Kutta-Joukowski subroutine cannot converge.'
             RETURN         
         END IF
          
      END DO     
      
      
        ! Record the iteration number of KJ
      IF (p%FVM%KJ_output) THEN
         CALL WRITE_KJ(Iiter, N, DG_MAX, 'NONE', p)
      END IF         
      
      
   END IF ! ( p%FVM%DS_Parms%relax_tune  .AND. O%WINDS_Timestep <=2 )
   
   
   
   

      
  
CONTAINS
   !...............................................................................................................................
   SUBROUTINE KJ(O, p, DG_MAX, Iiter, ErrStat, ErrMess)

         ! Passed variables      
      TYPE(AD_ParameterType),          INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_OtherStateType),         INTENT(INOUT)  :: O!therState ! Other/optimization states
      REAL(DbKi),                      INTENT(INOUT)  :: DG_MAX
      INTEGER(IntKi),                  INTENT(IN   )  :: Iiter     ! The counter of iteration
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMess

         ! Internal variables
      INTEGER(IntKi)             :: N 

      INTEGER(IntKi)             :: NST
      INTEGER(IntKi)             :: NS  
      INTEGER(IntKi)             :: NT  
      INTEGER(IntKi)             :: NB   
    
      INTEGER(IntKi)             :: IBlade     
      INTEGER(IntKi)             :: ITimestep
      INTEGER(IntKi)             :: IElement  
      INTEGER(IntKi)             :: IElement2
      
      REAL(DbKi)                 :: Temp1, Temp2         

      REAL(DbKi)                 :: CLA, CDA, CMA  ! For the CLCD from AeroDyn  
      INTEGER(IntKi)             :: Check_name
      
      
      CHARACTER(LEN = 1024)   :: TEMP10  ! debug
      
      NST = p%Element%NElm + 1
      NS  = p%Element%NElm
      NT  = p%FVM%NT
      NB  = p%NumBl
      
      N   =  O%WINDS_Timestep 
              
      O%FVM_Other%KJ%CL  =  O%FVM_Other%PERF_CL(1, N-1, 1, NS, NB)
      O%FVM_Other%KJ%CD  =  O%FVM_Other%PERF_CD(1, N-1, 1, NS, NB)
      O%FVM_Other%KJ%CM  =  O%FVM_Other%PERF_CM(1, N-1, 1, NS, NB)   
      O%FVM_Other%KJ%VEL_ROT = 0.0_DbKi
      O%FVM_Other%KJ%VEL_TOT = 0.0_DbKi
      O%FVM_Other%KJ%U       = 0.0_DbKi
      O%FVM_Other%KJ%V       = 0.0_DbKi
      O%FVM_Other%KJ%W       = 0.0_DbKi
      O%FVM_Other%KJ%Vinf    = 0.0_DbKi
      O%FVM_Other%KJ%Vtot    = 0.0_DbKi
      O%FVM_Other%KJ%AOA     = 0.0_DbKi  
      
      O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, 1:NS, 1:NB) = O%FVM_Other%KJ%GAMMA1(1, 1, 1, 1:NS, 1:NB)


         ! Compute induced velocity on lifting line due to shed and trailing filament induction
      CALL InducedVelocity(O, p, xd, N, 'KJ_ITER', ErrStat, ErrMess)
      !IF (ErrStat /= 0 ) RETURN

         ! Perform coordinate transformation on induced velocity (inertial to blade)    ! see SUBROUTINE WINDS_Velocity
      O%FVM_Other%KJ%VEL_ROT(1,1,1,1:NS,1:NB) = O%FVM_Other%POS_NODES_BXN(1,N,1,1:NS,1:NB) *                                &
                                                                                  O%FVM_Other%VEL_UINDB(1,N,1,1:NS,1:NB) +  &
                                     O%FVM_Other%POS_NODES_BXN(2,N,1,1:NS,1:NB) * O%FVM_Other%VEL_UINDB(2,N,1,1:NS,1:NB) +  &
                                     O%FVM_Other%POS_NODES_BXN(3,N,1,1:NS,1:NB) * O%FVM_Other%VEL_UINDB(3,N,1,1:NS,1:NB)
                       
      O%FVM_Other%KJ%VEL_ROT(2,1,1,1:NS,1:NB) = O%FVM_Other%POS_NODES_BYN(1,N,1,1:NS,1:NB) *                                 &
                                                                                   O%FVM_Other%VEL_UINDB(1,N,1,1:NS,1:NB) +  &  
                                      O%FVM_Other%POS_NODES_BYN(2,N,1,1:NS,1:NB) * O%FVM_Other%VEL_UINDB(2,N,1,1:NS,1:NB) +  & 
                                      O%FVM_Other%POS_NODES_BYN(3,N,1,1:NS,1:NB) * O%FVM_Other%VEL_UINDB(3,N,1,1:NS,1:NB)
    
      O%FVM_Other%KJ%VEL_ROT(3,1,1,1:NS,1:NB) = O%FVM_Other%POS_NODES_BZN(1,N,1,1:NS,1:NB) *                                 &
                                                                                   O%FVM_Other%VEL_UINDB(1,N,1,1:NS,1:NB) +  & 
                                      O%FVM_Other%POS_NODES_BZN(2,N,1,1:NS,1:NB) * O%FVM_Other%VEL_UINDB(2,N,1,1:NS,1:NB) +  &  
                                      O%FVM_Other%POS_NODES_BZN(3,N,1,1:NS,1:NB) * O%FVM_Other%VEL_UINDB(3,N,1,1:NS,1:NB)

         ! Compute effective wind in blade coordinate system
      O%FVM_Other%KJ%VEL_TOT(1:3, N, 1, 1:NS, 1:NB)  =  O%FVM_Other%VEL_BLADE(1:3, N, 1, 1:NS, 1:NB) +                       &
                                                                               O%FVM_Other%KJ%VEL_ROT(1:3, 1, 1, 1:NS, 1:NB) 
      O%FVM_Other%KJ%U(1, 1, 1, 1:NS, 1:NB)          =  O%FVM_Other%KJ%VEL_TOT(1, N, 1, 1:NS, 1:NB)
      O%FVM_Other%KJ%V(1, 1, 1, 1:NS, 1:NB)          =  O%FVM_Other%KJ%VEL_TOT(2, N, 1, 1:NS, 1:NB)
      O%FVM_Other%KJ%W(1, 1, 1, 1:NS, 1:NB)          =  O%FVM_Other%KJ%VEL_TOT(3, N, 1, 1:NS, 1:NB)      
      
      O%FVM_Other%KJ%Vinf(1, 1, 1, 1:NS, 1:NB) = SQRT(SUM(O%FVM_Other%VEL_BLADE(1:3, N, 1, 1:NS, 1:NB) ** 2, 1))
      O%FVM_Other%KJ%Vtot(1, 1, 1, 1:NS, 1:NB) = SQRT(SUM(O%FVM_Other%KJ%VEL_TOT(1:3, N, 1, 1:NS, 1:NB) ** 2, 1))         
         
         ! Compute angle of attack and sideslip angle
      O%FVM_Other%KJ%AOA   =  ATAN2(- O%FVM_Other%KJ%V, O%FVM_Other%KJ%U) ! angle in radian, not degrees...

                  
        ! Find coef. of lift and drag 
      IF (O%WINDS_Timestep > p%FVM%DS_Parms%start_n  .AND.  p%FVM%DS_Parms%DS_Flag ) THEN
          ! Call Leishman-Beddoes dynamic stall model
         DO IBlade = 1, NB   
            DO IElement = 1, NS
               Check_name =  INDEX((p%AirFoil%FOILNM(P%AirFoil%NFOIL(IElement))), "Cylinder")    
               IF (Check_name == 0) THEN ! Not "Cylinder"
                  CALL LB_DynStall(p, O, xd, Iiter, ErrStat, ErrMess, IElement, IBlade, CLA, CDA, CMA, P%AirFoil%NFOIL(IElement))                         
               ELSE  ! "Cylinder"
                   CALL CLCD_FVM(p, O, xd, ErrStat, ErrMess, O%FVM_Other%KJ%AOA(1, 1, 1, IElement, IBlade), CLA, CDA, CMA, & 
                                 P%AirFoil%NFOIL(IElement))    
               END IF   
               
               IF ( ErrStat /= ErrID_None )  THEN      
                  ErrMess = ' Error (in WInDS): Error occured when finding Cl in Kutta-Joukowski/LB_DynStall subroutine. Please check.'      
                  RETURN      
               END IF 
               
               IF ( ISNAN(CLA) )  THEN
                  ErrStat =  ErrID_Fatal
                  ErrMess = ' Error (in WInDS): Cl is Nan in Kutta-Joukowski/LB_DynStall subroutine. Please check.'   
                  RETURN      
               END IF 
               
               IF ( ISNAN(CDA) )  THEN
                  ErrStat =  ErrID_Fatal
                  ErrMess = ' Error (in WInDS): Cd is Nan in Kutta-Joukowski/LB_DynStall subroutine. Please check.'      
                  RETURN      
               END IF

               O%FVM_Other%KJ%CL(1,1,1,IElement,IBlade) = CLA
               O%FVM_Other%KJ%CD(1,1,1,IElement,IBlade) = CDA    
               O%FVM_Other%KJ%CM(1,1,1,IElement,IBlade) = CMA 
            END DO
         END DO    
          
      ELSE
            ! Interpolate over airfoil data tables    ! ~ SUBROUTINE CLCD in AeroSubs.f90         
         DO IBlade = 1, NB   
            DO IElement = 1, NS
                 ! The routine used in original AeroDyn
               CALL CLCD_FVM(p, O, xd, ErrStat, ErrMess, O%FVM_Other%KJ%AOA(1, 1, 1, IElement, IBlade), CLA, CDA, CMA, & 
                             P%AirFoil%NFOIL(IElement))                         
   
               IF ( ErrStat /= ErrID_None )  THEN      
                  ErrMess = ' Error (in WInDS): Error occured when finding Cl in Kutta-Joukowski subroutine. Please check.'      
                  RETURN      
               END IF 

               O%FVM_Other%KJ%CL(1,1,1,IElement,IBlade) = CLA
               O%FVM_Other%KJ%CD(1,1,1,IElement,IBlade) = CDA    
               O%FVM_Other%KJ%CM(1,1,1,IElement,IBlade) = CMA 
            END DO
         END DO
      END IF ! (O%WINDS_Timestep > p%FVM%DS_Parms%start_n  .AND.  p%FVM%DS_Parms%DS_Flag )   
         
 
        ! Compute bound vorticity via Kutta-Joukowski theorem
      O%FVM_Other%KJ%DG = 0.0
      
      DO IBlade = 1, NB           
         DO IElement = 1, NS
            O%FVM_Other%KJ%GAMMA1(1, 1, 1, IElement, IBlade) = 0.5 * O%FVM_Other%KJ%Vinf(1, 1, 1, IElement, IBlade) *    &
                                                         p%Blade%C(IElement) *  O%FVM_Other%KJ%CL(1, 1, 1, IElement, IBlade)
         END DO
      END DO         
        
        
        ! Compute performance variables and coefficients        
      O%FVM_Other%KJ%DG(1,1,1,1:NS,1:NB) = O%FVM_Other%KJ%GAMMA1(1,1,1,1:NS,1:NB) -                   &
                                                              O%FVM_Other%wake_gamma_shed(1,1,1,1:NS,1:NB)  !Change in bound vorticity between iterations

      O%FVM_Other%wake_gamma_shed(1,1,1,1:NS,1:NB)  =  O%FVM_Other%wake_gamma_shed(1,1,1,1:NS,1:NB) +     &
                                                           O%FVM_Other%KJ%Relax * O%FVM_Other%KJ%DG(1,1,1,1:NS,1:NB)


      Do IBlade = 1,NB    
         Do IElement2 = 1, NST     
            IF (IElement2 == 1) THEN ! First element
               O%FVM_Other%WAKE_gamma_trail(1,1,1,IElement2,IBlade) =      &                                  
                                                        O%FVM_Other%WAKE_gamma_shed(1,1,1,IElement2,IBlade)
            ELSE IF (IElement2 == NST) THEN  ! Last element
               O%FVM_Other%WAKE_gamma_trail(1,1,1,IElement2,IBlade) =      &   
                                                        - O%FVM_Other%WAKE_gamma_shed(1,1,1,IElement2-1,IBlade) 
            ELSE ! Not first or last
               O%FVM_Other%WAKE_gamma_trail(1,1,1,IElement2,IBlade) =      &
                                                           O%FVM_Other%WAKE_gamma_shed(1,1,1,IElement2,IBlade)    &
                                                      - O%FVM_Other%WAKE_gamma_shed(1,1,1,IElement2 -1,IBlade) 
            ENDIF  
         END DO !IElement2        
      END DO !IBlade      

     
      
      Do IBlade = 1, NB  
         Temp1 = 0.0 
         Temp2 = 0.0        
         
         Do IElement = 1, NS 
            Temp1 = Temp1 + O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1 ,1, IElement,  IBlade)
            Temp2 = Temp2 + O%FVM_Other%WAKE_GAMMA_TRAIL(1, 2 ,1, IElement,  IBlade)
            O%FVM_Other%WAKE_GAMMA_SHED(1,2,1,IElement,IBlade) = Temp2 - Temp1
         END DO
      END DO            

 
      
      O%FVM_Other%KJ%DG = O%FVM_Other%KJ%DG / ( ABS(O%FVM_Other%KJ%GAMMA1) + 1 )

      
      
      DG_MAX = 0.0
      
      Do IBlade = 1, NB  
         Do IElement = 1, NS       
            IF (DG_MAX  <= ABS(O%FVM_Other%KJ%DG(1, 1, 1, IElement, IBlade))  )  THEN
                DG_MAX = ABS(O%FVM_Other%KJ%DG(1, 1, 1, IElement, IBlade))
            ENDIF            
         END DO
      END DO          
      
      
         ! Compute performance variables and coefficients
      O%FVM_Other%PERF_CL(1,N,1,1:NS,1:NB)    =  O%FVM_Other%KJ%CL(1,1,1,1:NS,1:NB)
      O%FVM_Other%PERF_CD(1,N,1,1:NS,1:NB)    =  O%FVM_Other%KJ%CD(1,1,1,1:NS,1:NB)
      O%FVM_Other%PERF_AOA(1,N,1,1:NS,1:NB)   =  O%FVM_Other%KJ%AOA(1,1,1,1:NS,1:NB)
      
      O%FVM_Other%PERF_CM(1,N,1,1:NS,1:NB)    =  O%FVM_Other%KJ%CM(1,1,1,1:NS,1:NB) ! Matlab WInDS does not have this, but FAST needs.

      
   END SUBROUTINE KJ
   !............................................................
   
END SUBROUTINE KuttaJoukowski

!==================================================================================================================================
SUBROUTINE NUM_ADVECT(O, p, xd, N, ErrStat, ErrMess) 
! numerically convects wake nodes to next timestep.
!...........................................................

      ! Passed variables   
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states      
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep 
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess 

   
      ! Compute induced velocity at all wake points due to wake, store previous
   CALL InducedVelocity(O, p, xd, N, 'MAIN', ErrStat, ErrMess)
   IF (ErrStat /= 0 ) RETURN
      
      ! Numerically convect wake nodes to time+1   
   SELECT CASE  ( p%FVM%INTEG )  
          
         CASE ('RK4')         ! Runge-Kutta 4
            CALL Numerical_RK4(O, p, xd, N, ErrStat, ErrMess) 
            IF (ErrStat /= 0 ) RETURN
            
         CASE ('RK2')         ! Runge-Kutta 2
            CALL Numerical_RK2(O, p, xd, N, ErrStat, ErrMess)  
            IF (ErrStat /= 0 ) RETURN
  
         !CASE ('FE0')        ! Foward euler
         !   CALL Numerical_FE(O, p, xd, N, ErrStat, ErrMess)     
         !   IF (ErrStat /= 0 ) RETURN
    
  END SELECT ! p%FVM%integ    


   
CONTAINS
   !============================================================================
   SUBROUTINE Numerical_FE(O, p, xd, N, ErrStat, ErrMess) 
   ! Numerically convects wake nodes to time+1 via forward Euler numerical integration
   ! (~ numerical\fe.m)
   !....................................................................

  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! Timestep , n counts from 1
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess
   
   
      ! Internal variables
   REAL(DbKi)        :: DT         ! Duration of between timesteps

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS  
   INTEGER(IntKi)      :: NT  
   INTEGER(IntKi)      :: NB    

   
   INTEGER(IntKi)    :: IBlade
   INTEGER(IntKi)    :: ITimestep
   INTEGER(IntKi)    :: IElement
   INTEGER(IntKi)    :: IElement2

   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl

   DT   =   p%FVM%DT_WINDS

                   
      !.................................................................      
      ! Integrate to time+1 via forward Euler
     
   Do ITimestep = 2, N   
      O%FVM_Other%WAKE_DOMAIN(:, ITimestep +1 ,1 ,: , :)  =  O%FVM_Other%WAKE_DOMAIN(:,ITimestep,2,:,:) +    &
                                                              DT * O%FVM_Other%VEL_DOMAIN(:,ITimestep,1,:,:)
   END DO !ITimestep 
    
   O%FVM_Other%WAKE_DOMAIN(:,1,1,:,:)   =   O%FVM_Other%POS_QUARTER(:,N,1,:,:)  
   O%FVM_Other%WAKE_DOMAIN(:,2,1,:,:)   =   O%FVM_Other%POS_TRAIL(:,N,1,:,:)   

   CALL UPDATE_WAKE(O, p, xd, N+1, ErrStat, ErrMess,'FE') 


   END SUBROUTINE Numerical_FE

   !============================================================================
   SUBROUTINE Numerical_RK2(O, p, xd, N, ErrStat, ErrMess) 
   ! Numerically convects wake nodes to time+1 via 2nd-order Runge-Kutta numerical integration
   ! (~ numerical\rk2.m)
   !....................................................................
  
  IMPLICIT                        NONE
  
  
      ! Passed variables
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! Timestep , n counts from 1
   INTEGER,                       INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess
   
   
      ! Internal variables
   REAL(DbKi)        :: DT         ! Duration of between timesteps
  
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB    
   INTEGER(IntKi)      :: NTW
   
   INTEGER(IntKi)    :: IBlade
   INTEGER(IntKi)    :: ITimestep
   INTEGER(IntKi)    :: IElement
   INTEGER(IntKi)    :: IElement2
   REAL(DbKi)        :: Temp1, Temp2 
  
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NT  = p%FVM%NT
   NB  = p%NumBl
   NTW = O%FVM_Other%NTW
   
   DT  = p%FVM%DT_WINDS
  
       
      !-------------------------------------------------   
      ! % first step in rk2 (forward Euler)
      !-------------------------------------------------      
   Do ITimestep = 2, NTW   
      O%FVM_Other%WAKE_DOMAIN(:, ITimestep +1 ,1 ,: , :)  =  O%FVM_Other%WAKE_DOMAIN(:,ITimestep,2,:,:) +    &
                                                             DT * O%FVM_Other%VEL_DOMAIN(:,ITimestep,1,:,:)
   END DO !ITimestep 
    
   O%FVM_Other%WAKE_DOMAIN(:,1,1,:,:)   =   O%FVM_Other%POS_QUARTER(:,N,1,:,:)  
   O%FVM_Other%WAKE_DOMAIN(:,2,1,:,:)   =   O%FVM_Other%POS_TRAIL(:,N,1,:,:)      
  
       !-------------------------------------------------  
       ! Update shed and trailing filament strength
       !-------------------------------------------------
   CALL UPDATE_WAKE(O, p, xd, N+1, ErrStat, ErrMess, 'RK2')     
   IF (ErrStat /= 0 ) RETURN 
       
   
       !-------------------------------------------------  
       ! Update velocity domain based on new wake domain
       !-------------------------------------------------
    O%FVM_Other%RK_counter = 1
    CALL shear(O, p, 'RK', 1 ,N, ErrStat, ErrMess)
    IF (ErrStat /= 0 ) RETURN 
    
       !-------------------------------------------------------------------------------------------------------------------   
       ! Modify vortex core size via Ramasamy-Leishman model and include effect of filament stretching from previous timestep
       !------------------------------------------------- 
   CALL VCORE(p, O, xd, N + 1, ErrStat, ErrMess)   
   IF (ErrStat /= 0 ) RETURN 
      
       !------------------------------------------------- 
       ! Compute strength of new bound vortex via Kutta-Joukowski theorem
       !------------------------------------------------- 
   CALL KuttaJoukowski(O, p, xd, N, ErrStat, ErrMess)
   IF (ErrStat /= 0 ) RETURN 
      
        ! Compute induced velocity at all points
   O%FVM_Other%RK_counter  = 1     
   CALL InducedVelocity(O, p, xd, N, 'RK2', ErrStat, ErrMess)       
   IF (ErrStat /= 0 ) RETURN 
   
  
       ! Correct integration estimate 
   O%FVM_Other%WAKE_DOMAIN(:,3:(NTW+1),1,:,:) =  O%FVM_Other%WAKE_DOMAIN(:,2:NTW,2,:,:)  +       &
                                               DT/2 * (  O%FVM_Other%VEL_UIND(:,2:NTW,1,:,:)     &
                                               + O%FVM_Other%VEL_DOMAIN_rk(:,3:(NTW+1),1,:,:) )
  
  
   END SUBROUTINE Numerical_RK2

   !============================================================================
   SUBROUTINE Numerical_RK4(O, p, xd, N, ErrStat, ErrMess) 
   ! Numerically convects wake nodes to time+1 via 4th-order Runge-Kutta numerical integration
   ! (~ numerical\rk4.m)
   !....................................................................

  IMPLICIT                        NONE


      ! Passed variables
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! Timestep , n counts from 1
   INTEGER,                       INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess
   
   
      ! Internal variables
   REAL(DbKi)          :: DT         ! Duration of between timesteps

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS
   INTEGER(IntKi)      :: NT 
   INTEGER(IntKi)      :: NB    
   INTEGER(IntKi)      :: NTW

   INTEGER(IntKi)      :: IBlade
   INTEGER(IntKi)      :: ITimestep
   INTEGER(IntKi)      :: IElement
   INTEGER(IntKi)      :: IElement2
   REAL(DbKi)          :: Temp1, Temp2 

   CHARACTER(LEN = 1024)   :: TEMP10  ! debug
   
   
   NST =  p%Element%NElm + 1
   NS  =  p%Element%NElm
   NT  =  p%FVM%NT
   NB  =  p%NumBl
   NTW =  O%FVM_Other%NTW
   
   DT  =  p%FVM%DT_WINDS

          
      !==========================================================================================================================================   
      ! First step in rk4 (forward Euler)
      !=================================================      
   O%FVM_Other%WAKE_DOMAIN(:, 3:NTW+1 ,1 ,: , :)  =  O%FVM_Other%WAKE_DOMAIN(:, 2:NTW, 2, :, :) +    &
                                                     DT / 2 * O%FVM_Other%VEL_DOMAIN(:, 2:NTW, 1, :, :)    
   O%FVM_Other%WAKE_DOMAIN(:,1,1,:,:)   =   O%FVM_Other%POS_QUARTER(:,N,1,:,:)  
   O%FVM_Other%WAKE_DOMAIN(:,2,1,:,:)   =   O%FVM_Other%POS_TRAIL(:,N,1,:,:) 

    !-------------------------------------------------  
    ! Update shed and trailing filament strength
   CALL UPDATE_WAKE(O, p, xd, N+1, ErrStat, ErrMess,'RK4')     
   !IF (ErrStat /= 0 ) RETURN 
     
    !-------------------------------------------------  
    ! Update velocity domain based on new wake domain
   CALL shear(O, p, 'RK', 1 ,N, ErrStat, ErrMess)
   !IF (ErrStat /= 0 ) RETURN 
   
    !-------------------------------------------------------------------------------------------------------------------   
    ! Modify vortex core size via Ramasamy-Leishman model and include effect of filament stretching from previous timestep
   CALL VCORE(p, O, xd, N + 1, ErrStat, ErrMess)   
   !IF (ErrStat /= 0 ) RETURN 
   
    !------------------------------------------------- 
    ! Compute strength of new bound vortex via Kutta-Joukowski theorem
   CALL KuttaJoukowski(O, p, xd, N, ErrStat, ErrMess)
   !IF (ErrStat /= 0 ) RETURN 
   
    !============================================================================================================================================   
    ! Second step in rk4 
    !=================================================      
        ! Compute induced velocity at all points
    O%FVM_Other%RK_counter  = 1     
    CALL InducedVelocity(O, p, xd, N, 'RK4', ErrStat, ErrMess)       
    !IF (ErrStat /= 0 ) RETURN 

    
       !Correct integration estimate 
    O%FVM_Other%WAKE_DOMAIN(:,3:(NTW+1),1,:,:) =  O%FVM_Other%WAKE_DOMAIN(:,2:NTW,2,:,:)  +       &
                                                DT/2 *  O%FVM_Other%VEL_DOMAIN_rk(:,3:(NTW+1),1,:,:) 

   !  !-------------------------------------------------  
   !  ! Update shed and trailing filament strength
   !CALL UPDATE_WAKE(O, p, xd, N+1, ErrStat, ErrMess,'RK4')   
   !IF (ErrStat /= 0 ) RETURN 
   !   
     !-------------------------------------------------  
     ! Update velocity domain based on new wake domain
    CALL shear(O, p, 'RK', 2 ,N, ErrStat, ErrMess)
    IF (ErrStat /= 0 ) RETURN 
    
      !-------------------------------------------------------------------------------------------------------------------   
      ! Modify vortex core size via Ramasamy-Leishman model and include effect of filament stretching from previous timestep
    CALL VCORE(p, O, xd, N + 1, ErrStat, ErrMess)   
    IF (ErrStat /= 0 ) RETURN 
    
       !------------------------------------------------- 
       ! Compute strength of new bound vortex via Kutta-Joukowski theorem
    CALL KuttaJoukowski(O, p, xd, N, ErrStat, ErrMess)      
    IF (ErrStat /= 0 ) RETURN   
      
      !=====================================================================================================================================          
      !third step in rk4 
      !=================================================  
        ! Compute induced velocity at all points
    O%FVM_Other%RK_counter  = 2         
    CALL InducedVelocity(O, p, xd, N, 'RK4', ErrStat, ErrMess)       
    IF (ErrStat /= 0 ) RETURN 

       !Correct integration estimate 
    O%FVM_Other%WAKE_DOMAIN(:,3:(NTW+1),1,:,:) =  O%FVM_Other%WAKE_DOMAIN(:,2:NTW,2,:,:)  +       &
                                               DT *  O%FVM_Other%VEL_DOMAIN_rk(:,3:(NTW+1),2,:,:) 

   !  !-------------------------------------------------  
   !  ! Update shed and trailing filament strength
   !CALL UPDATE_WAKE(O, p, xd, N+1, ErrStat, ErrMess,'RK4') 
   !IF (ErrStat /= 0 ) RETURN 
   
     !-------------------------------------------------  
     ! Update velocity domain based on new wake domain
   CALL shear(O, p, 'RK', 3 ,N, ErrStat, ErrMess)
   IF (ErrStat /= 0 ) RETURN 
   
     !-------------------------------------------------------------------------------------------------------------------   
     ! Modify vortex core size via Ramasamy-Leishman model and include effect of filament stretching from previous timestep
   CALL VCORE(p, O, xd, N + 1, ErrStat, ErrMess)   
   IF (ErrStat /= 0 ) RETURN 
   
     !------------------------------------------------- 
     ! Compute strength of new bound vortex via Kutta-Joukowski theorem
   CALL KuttaJoukowski(O, p, xd, N, ErrStat, ErrMess)      
   IF (ErrStat /= 0 ) RETURN         
      
      
      !=====================================================================================================================================          
      ! fourth step in rk4 
      !=================================================        
        ! Compute induced velocity at all points
    O%FVM_Other%RK_counter  = 3        
    CALL InducedVelocity(O, p, xd, N, 'RK4', ErrStat, ErrMess)     
    IF (ErrStat /= 0 ) RETURN 

       !Correct integration estimate 
    O%FVM_Other%WAKE_DOMAIN(:,3:(NTW+1),1,:,:) =  O%FVM_Other%WAKE_DOMAIN(:,2:NTW,2,:,:)  +           &
                                                  DT/6 *( O%FVM_Other%VEL_DOMAIN(:,2:NTW,1,:,:)   +  &      
                                                         2 * O%FVM_Other%VEL_DOMAIN_rk(:,3:(NTW+1),1,:,:) +  &        
                                                         2 * O%FVM_Other%VEL_DOMAIN_rk(:,3:(NTW+1),2,:,:) +  &
                                                             O%FVM_Other%VEL_DOMAIN_rk(:,3:(NTW+1),3,:,:) )
      
      
   END SUBROUTINE Numerical_RK4

   !============================================================================
        
END SUBROUTINE NUM_ADVECT

!==================================================================================================================================
SUBROUTINE Shear(O, p, CALLER, rk_count, N, ErrStat, ErrMess)
! (~ shear.m)
! calculates the inflow velocity of all points based on selected shear model
!....................................................................

      ! Passed variables   
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState  ! Other/optimization states      
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   CHARACTER(*),                  INTENT(IN   )  :: CALLER      ! Who calls this subroutine 
   INTEGER(IntKi),                INTENT(IN   )  :: rk_count     
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! The number of timestep 
   INTEGER,                       INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess
   
       ! Internal variables
   INTEGER(IntKi)             :: NTW
   INTEGER(IntKi)             :: IDim    
   INTEGER(IntKi)             :: IBlade
   INTEGER(IntKi)             :: IElement2
   INTEGER(IntKi)             :: ITimestep       
 
   REAL(DbKi)                 :: point(3), U_inf(3), temp_vel(3)   
   NTW = O%FVM_Other%NTW
   
   SELECT CASE  ( CALLER )   
      CASE ('MAIN')       
         IF (p%FVM%Shear_Parms%ShearFLAG) THEN
             DO ITimestep = 1, NTW
                DO IBlade = 1, p%NumBl
                   Do IElement2 = 1, p%Element%NElm +1
                      point(1) = O%FVM_Other%WAKE_DOMAIN(1, ITimestep, 1, IElement2, IBlade)
                      point(2) = O%FVM_Other%WAKE_DOMAIN(2, ITimestep, 1, IElement2, IBlade)
                      point(3) = O%FVM_Other%WAKE_DOMAIN(3, ITimestep, 1, IElement2, IBlade)
                      
                      U_inf(1) = O%FVM_Other%WIND_INFTY(1, ITimestep, 1, 1, IBlade)
                      U_inf(2) = O%FVM_Other%WIND_INFTY(2, ITimestep, 1, 1, IBlade)
                      U_inf(3) = O%FVM_Other%WIND_INFTY(3, ITimestep, 1, 1, IBlade)
                      
                      CALL Shear_calc(point, U_inf, temp_vel, p, O, N, ErrStat, ErrMess)
                      IF (ErrStat /= 0 ) RETURN
                      
                      O%FVM_Other%VEL_DOMAIN(1, ITimestep, 1, IElement2, IBlade) = temp_vel(1)
                      O%FVM_Other%VEL_DOMAIN(2, ITimestep, 1, IElement2, IBlade) = temp_vel(2)                      
                      O%FVM_Other%VEL_DOMAIN(3, ITimestep, 1, IElement2, IBlade) = temp_vel(3)                   
                   END DO  ! IElement2
                END DO
             END DO
             
         ELSE
            Do IBlade = 1, p%NumBl
               Do ITimestep = 1, NTW   
                  Do IDim = 1, 3  
                     Do IElement2 = 1, p%Element%NElm + 1
                        O%FVM_Other%VEL_DOMAIN(IDim, ITimestep, 1, IElement2, IBlade ) =                   &
                                                             O%FVM_Other%WIND_INFTY(IDim, 1, 1, 1, IBlade)   ! sliu: what if not steady inflow
                     END DO  ! IElement2
                  END DO  ! IDim 
               END DO  ! ITimestep
            END DO  ! IBlade                                             
         END IF
     
      CASE ('RK')   
          
         IF (p%FVM%Shear_Parms%ShearFLAG) THEN
             DO ITimestep = 1, NTW +1
                DO IBlade = 1, p%NumBl
                   Do IElement2 = 1, p%Element%NElm +1
                      point(1) = O%FVM_Other%WAKE_DOMAIN(1, ITimestep, 1, IElement2, IBlade)
                      point(2) = O%FVM_Other%WAKE_DOMAIN(2, ITimestep, 1, IElement2, IBlade)
                      point(3) = O%FVM_Other%WAKE_DOMAIN(3, ITimestep, 1, IElement2, IBlade)
                      
                      !U_inf(1) = O%FVM_Other%WIND_INFTY(1, ITimestep, 1, 1, IBlade)
                      !U_inf(2) = O%FVM_Other%WIND_INFTY(2, ITimestep, 1, 1, IBlade)
                      !U_inf(3) = O%FVM_Other%WIND_INFTY(3, ITimestep, 1, 1, IBlade) 
                      
                      U_inf(1) = O%FVM_Other%WIND_INFTY(1, 1, 1, 1, IBlade)
                      U_inf(2) = O%FVM_Other%WIND_INFTY(2, 1, 1, 1, IBlade)
                      U_inf(3) = O%FVM_Other%WIND_INFTY(3, 1, 1, 1, IBlade)                      
                      
                      CALL Shear_calc(point, U_inf, temp_vel, p, O, N+1, ErrStat, ErrMess )
                      IF (ErrStat /= 0 ) RETURN
                      
                      O%FVM_Other%VEL_DOMAIN_rk(1, ITimestep, rk_count, IElement2, IBlade) = temp_vel(1)
                      O%FVM_Other%VEL_DOMAIN_rk(2, ITimestep, rk_count, IElement2, IBlade) = temp_vel(2)                      
                      O%FVM_Other%VEL_DOMAIN_rk(3, ITimestep, rk_count, IElement2, IBlade) = temp_vel(3)                  
                    
                   END DO  ! IElement2
                END DO
             END DO
             
         ELSE
            Do IBlade = 1, p%NumBl
               Do ITimestep = 1, NTW + 1  
                  Do IDim = 1, 3  
                     Do IElement2 = 1, p%Element%NElm + 1
                        O%FVM_Other%VEL_DOMAIN_RK(IDim, ITimestep, rk_count, IElement2, IBlade ) =     &
                                                                  O%FVM_Other%WIND_INFTY(IDim, 1, 1, 1, IBlade)
                     END DO ! IElement2
                  END DO  ! IDim 
               END DO  ! ITimestep
            END DO  ! IBlade                                             
         END IF
         
   END SELECT        
 
   
END SUBROUTINE Shear

!==================================================================================================================================
SUBROUTINE Shear_calc(point, U_inf, temp_vel, p, O, N, ErrStat, ErrMess)
! calculates the inflow velocity of one point based on selected shear model
! Called from 'WINDS_Velocity' or 'Shear'.
!...............................................................................

        ! Passed variables   
      REAL(DbKi),DIMENSION(3),         INTENT(IN   )  :: point
      REAL(DbKi),DIMENSION(3),         INTENT(IN   )  :: U_inf
      REAL(DbKi),DIMENSION(3),         INTENT(  OUT)  :: temp_vel    
      TYPE(AD_OtherStateType),         INTENT(INOUT)  :: O!therState  ! Other/optimization states      
      TYPE(AD_ParameterType),          INTENT(IN   )  :: p           ! Parameters
      INTEGER(IntKi),                  INTENT(IN   )  :: N           ! The number of timestep 
      INTEGER,                       INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess
      
      
        ! Internal variables      
      REAL(DbKi)                 :: z
      REAL(DbKi)                 :: z_ref
      REAL(DbKi)                 :: Alpha
      REAL(DbKi)                 :: z0
      
    
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""
   
        

      z      =  point(3)
      z_ref  =  p%FVM%Shear_Parms%z_ref
      
      IF (p%FVM%Shear_Parms%model_type == 1) THEN       
         Alpha       =  p%FVM%Shear_Parms%Alpha
         temp_vel(1) =  U_inf(1) * (z/z_ref) ** alpha
         temp_vel(2) =  U_inf(2) * (z/z_ref) ** alpha
         temp_vel(3) =  U_inf(3)      
         
      ELSEIF (p%FVM%Shear_Parms%model_type == 2) THEN 
         z0        =    p%FVM%Shear_Parms%z0
         temp_vel(1) =  U_inf(1) * log(z/z0) / log(z_ref/z0)
         temp_vel(2) =  U_inf(2) * log(z/z0) / log(z_ref/z0)
         temp_vel(3) =  U_inf(3)            
         
      END IF   
   
   
END SUBROUTINE Shear_calc
   
!====================================================================================================
SUBROUTINE TOWER_SHADOW(O, p, xd, N, ErrStat, ErrMess) 
! WInDS Tower Shadow Velocity Deficit Code
!..........................................................

   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState 
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! Timestep , n counts from 1
   INTEGER,                       INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess

      ! Internal variables
   INTEGER(IntKi)             :: IDim       ! Index for dimensions
   INTEGER(IntKi)             :: IBlade     ! Index for blade
   INTEGER(IntKi)             :: IElement   ! Index for blade stations
   INTEGER(IntKi)             :: IElement2  ! Index for blade trailing nodes
   INTEGER(IntKi)             :: ITime 

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS 
   INTEGER(IntKi)      :: NB    
    
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""
   
     

   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NB  = p%NumBl
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   O%FVM_Other%VEL_UIND_TOWER = 0.0
   ! O%FVM_Other%VEL_UIND_TOWER(:,N,1,:,:) 
   

END SUBROUTINE TOWER_SHADOW

!====================================================================================================
SUBROUTINE UPDATE_WAKE(O, p, xd, N_plus, ErrStat, ErrMess, CALLER) 
! updates bound vortices and latest trailed and shed vortices
!....................................................................
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState 
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N_plus       ! Timestep , n counts from 1
   INTEGER,                       INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess
   CHARACTER(*),                  INTENT(IN   )  :: CALLER      ! Who calls this subroutine 

   
      ! Internal variables
   INTEGER(IntKi)             :: IDim       ! Index for dimensions
   INTEGER(IntKi)             :: IBlade     ! Index for blade
   INTEGER(IntKi)             :: IElement   ! Index for blade stations
   INTEGER(IntKi)             :: IElement2  ! Index for blade trailing nodes
   INTEGER(IntKi)             :: ITime 

   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS 
   INTEGER(IntKi)      :: NB 
   REAL(DbKi)          :: DT
   REAL(DbKi), DIMENSION(1,1,1,p%Element%NElm,p%NumBl)    :: b0, b1, b2
   REAL(DbKi)          :: Temp1, Temp2  
   INTEGER(IntKi)      :: NTW
   
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""
   
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NB  = p%NumBl
   DT  = p%FVM%DT_WINDS
   NTW = O%FVM_Other%NTW +1  
   
   IF (O%WINDS_Timestep <= 8  .OR. (.NOT.p%FVM%extrap_wake)  .OR. CALLER=='AB') THEN   
      O%FVM_Other%WAKE_GAMMA_SHED(:,1,1,:,:) = O%FVM_Other%WAKE_GAMMA_SHED(:,1,2,:,:)
   ELSE
      b0(1,1,1,:,:) = O%FVM_Other%WAKE_GAMMA_SHED(1,1,4,:,:)
      b1(1,1,1,:,:) = (O%FVM_Other%WAKE_GAMMA_SHED(1,1,3,:,:) - O%FVM_Other%WAKE_GAMMA_SHED(1,1,4,:,:))/ DT
      b2(1,1,1,:,:) = (O%FVM_Other%WAKE_GAMMA_SHED(1,1,2,:,:) - 2 * O%FVM_Other%WAKE_GAMMA_SHED(1,1,3,:,:) + O%FVM_Other%WAKE_GAMMA_SHED(1,1,4,:,:))/ (2*DT*DT)
      
      SELECT CASE  (CALLER)
         CASE ('FE','RK2') 
              O%FVM_Other%WAKE_GAMMA_SHED(1,1,1,:,:) = b0(1,1,1,:,:) + 3*b1(1,1,1,:,:)*DT + 6*b2(1,1,1,:,:)*DT*DT
         CASE ('RK4')  
              O%FVM_Other%WAKE_GAMMA_SHED(1,1,1,:,:) = b0(1,1,1,:,:) + (5/2)*b1(1,1,1,:,:)*DT + (15/4)*b2(1,1,1,:,:)*DT*DT
      END SELECT
      
   END IF
      
      
      !..........................................................................................
      ! Compute spanwise change in bound filament to compute first set of trailing filaments            
   DO IBlade = 1,NB    
      DO IElement2 = 1, NST     
         IF (IElement2 == 1) THEN ! First element
            O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1, 1, IElement2, IBlade) = O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2, IBlade)
            
         ELSE IF (IElement2 == NST) THEN  ! Last element
            O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1, 1, IElement2, IBlade) = - O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2-1, IBlade) 
            
         ELSE ! Not first or last
            O%FVM_Other%WAKE_GAMMA_TRAIL(1, 1, 1, IElement2, IBlade) = O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2, IBlade) &
                                                                    - O%FVM_Other%WAKE_GAMMA_SHED(1, 1, 1, IElement2 -1, IBlade) 
         ENDIF  
      END DO !IElement2        
   END DO !IBlade      


      !..........................................................................................   
      !Previous set of trailing filaments becomes new set of trailing filaments
   O%FVM_Other%WAKE_GAMMA_TRAIL(:, 2:NTW, 1, :, :)  =   O%FVM_Other%WAKE_GAMMA_TRAIL(:, 1:NTW-1, 2, :, :);   
   
   
      !..........................................................................................
      ! Shed filaments computed via spanwise summation of trailing filaments (ensure Kelvin's theorem is satisfied)
   Do IBlade = 1, NB        
      DO ITime  = 1, NTW -1 
         Temp1 = 0.0 
         Temp2 = 0.0 
         Do IElement = 1, NS    
            Temp1 = Temp1 + O%FVM_Other%WAKE_GAMMA_TRAIL(1, ITime, 1, IElement, IBlade)
            Temp2 = Temp2 + O%FVM_Other%WAKE_GAMMA_TRAIL(1, ITime +1 ,1, IElement,  IBlade)
            O%FVM_Other%WAKE_GAMMA_SHED(1, ITime+1, 1, IElement, IBlade) = Temp2 - Temp1 
         END DO   
      END DO
      
      ! Last timestep:
      Temp1 = 0.0 
      Do IElement = 1, NS    
         Temp1 = Temp1 + O%FVM_Other%WAKE_GAMMA_TRAIL(1, NTW, 1, IElement, IBlade)
         O%FVM_Other%WAKE_GAMMA_SHED(1, NTW +1, 1, IElement, IBlade) =  - Temp1 
      END DO        
   END DO            
   
     
   
END SUBROUTINE UPDATE_WAKE

!==================================================================================================================================
SUBROUTINE VCORE(p, O, xd, N, ErrStat, ErrMess)
! Computes the effective vortex filament core size using the Ramasamy-Leishman model and filament stretching.
! (~ vcore.m)
!....................................................................

  IMPLICIT                        NONE

      ! Passed variables
   TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState 
   TYPE(AD_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER(IntKi),                INTENT(IN   )  :: N           ! Timestep , n counts from 1
   INTEGER,                       INTENT(INOUT)  :: ErrStat
   CHARACTER(*),                  INTENT(INOUT)  :: ErrMess   
   
   
       ! Internal variables
   INTEGER(IntKi)      :: NST 
   INTEGER(IntKi)      :: NS 
   INTEGER(IntKi)      :: NB  
   INTEGER(IntKi)      :: NTW

   INTEGER(IntKi)             :: IDim    
   INTEGER(IntKi)             :: IBlade
   INTEGER(IntKi)             :: IElement
   INTEGER(IntKi)             :: IElement2
   INTEGER(IntKi)             :: ITimestep       
   REAL(DbKi)                 :: T0         ! T0 define initial vortex core size using Ramasamy-Leishman model
   
   CHARACTER(LEN = 1024)   :: TEMP10  ! debug
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NB  = p%NumBl
   NTW = O%FVM_Other%NTW +1  

      

      ! Define initial vortex core size       (sliu: this part should be moved to initials.m)
   O%FVM_Other%TipSpdRat = p%Blade%R * O%Rotor%REVS / O%FVM_Other%WIND_INFTYM(1,1,1,1,1)    ! REVS: Rotor rotational speed (i.e. RPM in rad/sec)
   T0 = TwoPi * p%Blade%TipRadius /(12 * O%FVM_Other%TipSpdRat * O%FVM_Other%WIND_INFTYM(1,1,1,1,1)) 
   
   O%FVM_Other%WAKE_R0(1) = SQRT(4 * p%FVM%RL_Model%ALPHA  * p%FVM%RL_Model%NU * p%FVM%RL_Model%DELTA * T0)   

   
         ! Compute vortex Re #
   Do IElement = 1, NS
      DO ITimestep = 1,(NTW+1)
         DO IBlade = 1, NB
            O%FVM_Other%WAKE_RE_SHED(:, ITimestep, 1, IElement, IBlade) =              &
                                    ABS(O%FVM_Other%WAKE_GAMMA_SHED(:, ITimestep, 1, IElement, IBlade)/ p%FVM%RL_Model%NU)
         END DO
      END DO  
   END DO
   

   
   Do IElement2 = 1, NST
      DO ITimestep = 1, NTW
         DO IBlade = 1, NB     
            O%FVM_Other%WAKE_RE_TRAIL(1, ITimestep, 1, IElement2, IBlade) =              &
                                         ABS(O%FVM_Other%WAKE_GAMMA_TRAIL(1, ITimestep, 1, IElement2, IBlade)/ p%FVM%RL_Model%NU)
         END DO
      END DO  
   END DO
    
    
      ! Modify coresize using Ramasamy-Leishman model  
   Do IElement = 1, NS
      DO ITimestep = 1,(NTW+1)
         DO IBlade = 1, NB    
            O%FVM_Other%WAKE_RC_SHED(1, ITimestep, 1, IElement, IBlade) = (O%FVM_Other%WAKE_R0(1)**2 +      &
                                4 * p%FVM%RL_Model%ALPHA * p%FVM%RL_Model%NU *( 1 + p%FVM%RL_Model%A1 *          &
                                O%FVM_Other%WAKE_RE_SHED(1, ITimestep, 1, IElement, IBlade)) * (p%FVM%DT_WINDS) * (N-1))
         END DO
      END DO  
   END DO        
        
   Do IElement2 = 1, NST
      DO ITimestep = 1, NTW
         DO IBlade = 1, NB    
            O%FVM_Other%WAKE_RC_TRAIL(1, ITimestep, 1, IElement2, IBlade) = (O%FVM_Other%WAKE_R0(1) **2 +       &
                            4 * p%FVM%RL_Model%ALPHA * p%FVM%RL_Model%NU *( 1 + p%FVM%RL_Model%A1 *                  &
                            O%FVM_Other%WAKE_RE_TRAIL(1, ITimestep, 1, IElement2, IBlade)) * (p%FVM%DT_WINDS) * (N-1))
         END DO
      END DO  
   END DO         
             
   Do IElement = 1, NS
      DO ITimestep = 1,(NTW+1)
         DO IBlade = 1, NB                
            O%FVM_Other%WAKE_RC_EFF_SHED(1, ITimestep, 1, IElement, IBlade) =                 &
                                                       O%FVM_Other%WAKE_RC_SHED(1, ITimestep, 1, IElement, IBlade)
         END DO
      END DO  
   END DO        
   
   Do IElement2 = 1, NST
      DO ITimestep = 1, NTW
         DO IBlade = 1, NB   
            O%FVM_Other%WAKE_RC_EFF_TRAIL(1, ITimestep, 1, IElement2, IBlade)  =               & 
                                                      O%FVM_Other%WAKE_RC_TRAIL(1, ITimestep, 1, IElement2, IBlade)
         END DO
      END DO  
   END DO                   
    
    
      ! Determine filament lengths, then apply filament stretching
   IF (O%FVM_Other%ntroll > 0) THEN
      Do IElement2 = 1, NST-1
         DO ITimestep = 1,(NTW+1)
            DO IBlade = 1, NB  
               O%FVM_Other%WAKE_LENGTH_SHED(1, ITimestep, 1, IElement2, IBlade) =                                                &
                                  ( O%FVM_Other%WAKE_DOMAIN(1, ITimestep, 1, IElement2, IBlade) -                                &
                                                         O%FVM_Other%WAKE_DOMAIN(1, ITimestep, 1, IElement2 +1, IBlade) )** 2    &     
                                +  ( O%FVM_Other%WAKE_DOMAIN(2, ITimestep, 1, IElement2, IBlade) -                               &
                                                         O%FVM_Other%WAKE_DOMAIN(2, ITimestep, 1, IElement2 +1, IBlade) )** 2    &
                                +   ( O%FVM_Other%WAKE_DOMAIN(3, ITimestep, 1, IElement2, IBlade) -                           &
                                                          O%FVM_Other%WAKE_DOMAIN(3, ITimestep, 1, IElement2 +1, IBlade) )** 2     ! Length ^2
            END DO
         END DO  
      END DO  
      
      Do IElement2 = 1, NST
         DO ITimestep = 1, NTW
            DO IBlade = 1, NB  
               O%FVM_Other%WAKE_LENGTH_TRAIL(1, ITimestep, 1, IElement2, IBlade) =                                             &
                               (  O%FVM_Other%WAKE_DOMAIN(1, ITimestep +1, 1, IElement2, IBlade) -                             &
                                                        O%FVM_Other%WAKE_DOMAIN(1, ITimestep, 1, IElement2, IBlade) )** 2      &
                            +  (  O%FVM_Other%WAKE_DOMAIN(2, ITimestep +1, 1, IElement2, IBlade) -                             &
                                                        O%FVM_Other%WAKE_DOMAIN(2, ITimestep, 1, IElement2, IBlade) )** 2      &                
                            +  (  O%FVM_Other%WAKE_DOMAIN(3, ITimestep +1, 1, IElement2, IBlade) -                             &
                                                        O%FVM_Other%WAKE_DOMAIN(3, ITimestep, 1, IElement2, IBlade) )** 2      ! Length ^2
            END DO
         END DO  
      END DO      
   
       !Effective vortex filament core size due to filament stretching between current time and time-1  
      CALL FilamentMod(p, O, N)
   END IF  ! (p%FVM%roll)   
    
CONTAINS
   !...............................................................................................................................
   SUBROUTINE Filamentmod(p, O, N)   
   ! Computes the effective vortex filament core size due to filament stretching between timesteps.
   !....................................................................... 

      TYPE(AD_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_OtherStateType),       INTENT(INOUT)  :: O!therState 
      INTEGER(IntKi),                INTENT(IN   )  :: N           ! Timestep , n counts from 1    
    
          ! Internal variables
      INTEGER(IntKi)      :: NST
      INTEGER(IntKi)      :: NS 
      INTEGER(IntKi)      :: NB     
      INTEGER(IntKi)      :: NTW
      
      INTEGER(IntKi)              :: IDim    
      INTEGER(IntKi)              :: IBlade
      INTEGER(IntKi)              :: IElement
      INTEGER(IntKi)              :: IElement2
      INTEGER(IntKi)              :: ITimestep       
       
      REAL(DbKi), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: TRAILNEW
      REAL(DbKi), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: TRAILOLD      
      REAL(DbKi), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: SHEDNEW      
      REAL(DbKi), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: SHEDOLD     
      
      REAL(DbKi), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: STRAIN_TRAIL 
      REAL(DbKi), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: STRAIN_SHED   
   
      NST = p%Element%NElm + 1
      NS  = p%Element%NElm
      NB  = p%NumBl
      NTW = O%FVM_Other%NTW +1  
      
      IF (.NOT. ALLOCATED( TRAILNEW ) )      ALLOCATE(TRAILNEW(1, NTW, 1, NST, NB))    
      IF (.NOT. ALLOCATED( TRAILOLD ) )      ALLOCATE(TRAILOLD(1, NTW, 1, NST, NB)) 
      IF (.NOT. ALLOCATED( SHEDNEW ) )       ALLOCATE(SHEDNEW(1, NTW, 1, NS,NB)) 
      IF (.NOT. ALLOCATED( SHEDOLD ) )       ALLOCATE(SHEDOLD(1, NTW, 1, NS,NB))      
      IF (.NOT. ALLOCATED( STRAIN_TRAIL ) )  ALLOCATE(STRAIN_TRAIL(1, NTW, 1, NST, NB)) 
      IF (.NOT. ALLOCATED( STRAIN_SHED ) )   ALLOCATE(STRAIN_SHED(1, NTW, 1, NS,NB))                
    
      IF (N>3) THEN
         TRAILNEW(1, 1:NTW, 1, 1:NST, 1:NB) = SQRT(O%FVM_Other%WAKE_LENGTH_TRAIL(1, 1:NTW, 1, 1:NST, 1:NB))       ! Note that wake.length is really L^2
         TRAILOLD(1, 1:NTW, 1, 1:NST, 1:NB) = SQRT(O%FVM_Other%WAKE_LENGTH_TRAIL(1, 1:NTW, 2, 1:NST, 1:NB))
         SHEDNEW(1, 1:NTW, 1, 1:NS, 1:NB)   = SQRT(O%FVM_Other%WAKE_LENGTH_SHED(1, 1:NTW, 1, 1:NS, 1:NB))
         SHEDOLD(1, 1:NTW, 1, 1:NS, 1:NB)   = SQRT(O%FVM_Other%WAKE_LENGTH_SHED(1, 1:NTW, 2, 1:NS, 1:NB))
         
            ! Compute strain of trailing and shed filaments
         STRAIN_TRAIL = (TRAILNEW - TRAILOLD) / TRAILOLD
         STRAIN_SHED  = (SHEDNEW  - SHEDOLD) / SHEDOLD

            ! Equations modified as rc and re_eff are squared    
         O%FVM_Other%WAKE_RC_EFF_TRAIL(1, 1:NTW, 1, 1:NST, 1:NB) = O%FVM_Other%WAKE_RC_TRAIL(1, 1:NTW, 1, 1:NST, 1:NB) *   &
                                                                   (1 / (1 + STRAIN_TRAIL(1, 1:NTW, 1, 1:NST, 1:NB)))
         O%FVM_Other%WAKE_RC_EFF_SHED(1, 1:NTW, 1, 1:NS, 1:NB)  = O%FVM_Other%WAKE_RC_SHED(1, 1:NTW, 1, 1:NS, 1:NB) *      &
                                                                   (1 / (1 + STRAIN_SHED(1, 1:NTW, 1, 1:NS, 1:NB)))
      END IF ! N>3 


      IF (ALLOCATED( TRAILNEW ) )      DEALLOCATE(TRAILNEW) 
      IF (ALLOCATED( TRAILOLD ) )      DEALLOCATE(TRAILOLD) 
      IF (ALLOCATED( SHEDNEW ) )       DEALLOCATE(SHEDNEW) 
      IF (ALLOCATED( SHEDOLD ) )       DEALLOCATE(SHEDOLD)      
      IF (ALLOCATED( STRAIN_TRAIL ) )  DEALLOCATE(STRAIN_TRAIL) 
      IF (ALLOCATED( STRAIN_SHED ) )   DEALLOCATE(STRAIN_SHED)    

    
    
   END SUBROUTINE Filamentmod
   !...............................................................................................................................

   
END SUBROUTINE VCORE


!==================================================================================================================================
!==================================================================================================================================
!========================                                                                             =============================
!========================                                                                             =============================
!========================                                                                             =============================
!========================           The lines below are modified from AeroDyn v14.02.00c-mlb          =============================
!========================                                                                             =============================
!========================                                                                             =============================
!========================                                                                             =============================
!=================================================================================================================================
!==================================================================================================================================
SUBROUTINE BEM(u, p, xd, O, ErrStat, ErrMess)
! This part is modified from AeroDyn v14.02.00c-mlb and Matlab WInDS
! Some comments come from:  Algorithmic Outline of Unsteady Aerodynamics (AERODYN) Modules. Project WE-201103. Rick Damiani, Ph.D., P.E.
!============================      

      TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O!therState ! Other/optimization states

      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None

      ! Local variables
      INTEGER(IntKi)      :: NST != p%FVM%NST
      INTEGER(IntKi)      :: NS  != p%FVM%NS
      INTEGER(IntKi)      :: NT  != p%FVM%NT
      INTEGER(IntKi)      :: NB  != p%FVM%NB     

      INTEGER(IntKi)      :: IDim    ! For dimensions
      INTEGER(IntKi)      :: IBlade
      INTEGER(IntKi)      :: IElement


      REAL(DbKi)                 :: Diff_A   ! axial induction factor residual
      REAL(DbKi)                 :: Diff_AP  ! tangential induction factor residual
      REAL(DbKi)                 :: SPitch                     ! sine of PitNow
      REAL(DbKi)                 :: CPitch                     ! cosine of PitNow   
      REAL(DbKi)                 :: tmpVector     (3)    
      REAL(DbKi)                 :: rLocal                     ! Local radial distance of element in the plane of rotation
      REAL(DbKi)                 :: VelocityVec   (3)     
      REAL(DbKi)                 :: VelNormalToRotor2   
      REAL(DbKi)                 :: VTTotal   
      REAL(DbKi)                 :: VNElement
      REAL(DbKi)                 :: VNWind

      REAL(DbKi)                 :: LAMBDAR ! Local speed ratio
      REAL(DbKi)                 :: SIGMAP  ! Local solidity  ( = SOLFACT * VNROTOR2, in AeroDyn)
      REAL(DbKi)                 :: TWST   
      REAL(DbKi)                 :: PTCH  
   
      REAL(DbKi)                 :: A
      REAL(DbKi)                 :: AP   
      LOGICAL                    :: Iter_flag
   
      REAL(DbKi)                 :: QA
      REAL(DbKi)                 :: CLA
      REAL(DbKi)                 :: CDA
      REAL(DbKi)                 :: CMA
   
      REAL(DbKi)                 :: DFN
      REAL(DbKi)                 :: DFT
      REAL(DbKi)                 :: PMA

      REAL(DbKi)                 :: PHI
      REAL(DbKi)                 :: ALPHA

      REAL(DbKi)                 :: CPHI, SPHI
      REAL(DbKi)                 :: W2
    
      INTEGER(IntKi)             :: Itera_counter      
       
      
      NST = p%Element%NElm + 1
      NS  = p%Element%NElm
      
      NT  = p%FVM%NT
      NB  = p%NumBl

      
      DO IElement = 1, NS  
          
         IBlade = 1     ! Assume all the blades have some Cl and Cd 
          
         ! element pitch angle
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
         !-------------------------------------------------------------------------------------------         
         VelocityVec(1) = O%FVM_Other%WIND_INFTY(1, 1, 1, IElement, IBlade)
         VelocityVec(2) = O%FVM_Other%WIND_INFTY(2, 1, 1, IElement, IBlade) 
         VelocityVec(3) = O%FVM_Other%WIND_INFTY(3, 1, 1, IElement, IBlade)
         
         VelNormalToRotor2 = ( VelocityVec(3) * o%Rotor%STilt + (VelocityVec(1) * o%Rotor%CYaw               &
                             - VelocityVec(2) * o%Rotor%SYaw) * o%Rotor%CTilt )**2                                    ! Square of the wind velocity component normal to the rotor

         tmpVector =  -1.*SPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                        + CPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VTTotal   =     DOT_PRODUCT( tmpVector, VelocityVec - u%InputMarkers(IBlade)%TranslationVel(:,IElement)  )   ! Component in the plane of rotation of (wind ?Deflection translational velocity )

         tmpVector =     CPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                       + SPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VNWind    =     DOT_PRODUCT( tmpVector, VelocityVec )                                                        ! Component normal to the plane of rotation of wind velocity alone
         VNElement = -1.*DOT_PRODUCT( tmpVector, u%InputMarkers(IBlade)%TranslationVel(:,IElement ) )                 ! Component normal to the plane of rotation of deflection translational velocity alone      
          
         LAMBDAR =  p%BLADE%RNodes(IElement) * O%Rotor%REVS / O%FVM_Other%WIND_INFTYM(1, 1, 1, 1, 1)                  ! Local speed ratio
         SIGMAP  =  p%NumBl * P%Blade%C(IElement) / ( TWOPI * p%BLADE%RNodes(IElement))                            ! Local solidity  
         TWST    =  - P%Element%TWIST(IElement) 
         PTCH    =  o%Element%PitNow 

         !-------------------------------------------------------------------------------------------
         ! Initial values for axial and tangential induction factors
         !-------------------------------------------------------------------------------------------
         A  = REAL( 0.25 * ( 2 + Pi * lambdar * sigmap -         &
                       SQRT(4 - 4 * Pi * lambdar * sigmap + Pi * lambdar ** 2 * sigmap * (8*(twst + ptch) + Pi * sigmap))))
         AP = 0.0
         
         Iter_flag = .TRUE.
         Itera_counter = 0
         
         DO WHILE (Iter_flag) 
             CALL BEM_iteration(Itera_counter, Iter_flag, PHI, A, AP, CLA, CDA, CMA, ALPHA,   &
                              VNWind, VNElement, VTTotal, W2, SIGMAP, DIFF_A, DIFF_AP, IElement, O, p, ErrStat, ErrMess)
             IF (ErrStat /= 0 ) RETURN
         END DO
         
         
         ! Copy parameters and calculate force
         DO IBlade = 1, NB
            O%FVM_Other%PERF_CL(1, 1, 1, IElement, IBlade)    =  CLA
            O%FVM_Other%PERF_CD(1, 1, 1, IElement, IBlade)    =  CDA
            O%FVM_Other%PERF_CM(1, 1, 1, IElement, IBlade)    =  CMA
            O%FVM_Other%PERF_AOA(1, 1, 1, IElement, IBlade)   =  ALPHA            
            
            
            QA       = 0.5 * P%Wind%RHO * W2 * P%Blade%DR(IElement) * P%Blade%C(IElement)
            CPHI     = COS( PHI )
            SPHI     = SIN( PHI )
            DFN      = ( CLA * CPHI + CDA * SPHI ) * QA
            DFT      = ( CLA * SPHI - CDA * CPHI ) * QA

            PMA  = CMA * QA * P%Blade%C(IElement)    
         
            O%FVM_Other%StoredForces(1, 1, 1, IElement, IBlade)   = ( DFN*CPitch + DFT*SPitch ) / p%Blade%DR(IElement)
            O%FVM_Other%StoredForces(2, 1, 1, IElement, IBlade)   = ( DFN*SPitch - DFT*CPitch ) / p%Blade%DR(IElement)
            O%FVM_Other%StoredForces(3, 1, 1, IElement, IBlade)   = 0.0

            O%FVM_Other%StoredMoments(1, 1, 1, IElement, IBlade)  = 0.0
            O%FVM_Other%StoredMoments(2, 1, 1, IElement, IBlade)  = 0.0
            O%FVM_Other%StoredMoments(3, 1, 1, IElement, IBlade)  = PMA / p%Blade%DR(IElement)                 
         END DO !IBlade         
         
      END DO !IElement 
      
      
CONTAINS
!=============================================   
   SUBROUTINE BEM_iteration(Itera_counter, Iter_flag, PHI, A, AP, CLA, CDA, CMA, ALPHA,  &
     VNWind, VNElement, VTTotal, W2, SIGMAP, DIFF_A, DIFF_AP, J, O, p, ErrStat, ErrMess)
!..............................................   

      INTEGER(IntKi),            INTENT(INOUT)  :: Itera_counter  ! Iteration counter
      LOGICAL,                   INTENT(  OUT)  :: Iter_flag      ! Whether continue iteration
      REAL(DbKi),                INTENT(  OUT)  :: PHI            ! Spanwise inflow Angle
      REAL(DbKi),                INTENT(INOUT)  :: A              ! Spanwise axial induction factor
      REAL(DbKi),                INTENT(INOUT)  :: AP             ! Spanwise tangential induction factor
      REAL(DbKi),                INTENT(INOUT)  :: CLA            ! Spanwise lift coefficient
      REAL(DbKi),                INTENT(INOUT)  :: CDA            ! Spanwise drag coefficient 
      REAL(DbKi),                INTENT(INOUT)  :: CMA            ! Moment coefficient 
      REAL(DbKi),                INTENT(INOUT)  :: ALPHA          ! Spanwise angle of attack
      REAL(DbKi),                INTENT(IN   )  :: VNWind         ! Component of wind velocity normal to rotational plane (rotor)
      REAL(DbKi),                INTENT(IN   )  :: VNElement      ! Component of the translational velocity normal to plane of rotor 
      REAL(DbKi),                INTENT(IN   )  :: VTTotal        ! Tangential velocity accounting for wind and rotation/deflections of blade
      REAL(DbKi),                INTENT(INOUT)  :: W2             ! Relative velocity squared over blade element
      REAL(DbKi),                INTENT(IN   )  :: SIGMAP         ! Local solidity 
      REAL(DbKi),                INTENT(INOUT)  :: DIFF_A         ! Axial induction factor residual
      REAL(DbKi),                INTENT(INOUT)  :: DIFF_AP        ! Tangential induction factor residual
      INTEGER(IntKi),            INTENT(IN   )  :: J              ! Element index (IElement)   
      TYPE(AD_OtherStateType),   INTENT(INOUT)  :: O!therState          ! Other/optimization states   
      TYPE(AD_ParameterType),    INTENT(IN   )  :: p                    ! Parameters
   
      INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat              ! Error status of the operation
      CHARACTER(*),              INTENT(  OUT)  :: ErrMess              ! Error message if ErrStat /= ErrID_None
      
   
      ! Local variables
      REAL(DbKi)                 :: A0    ! Old value for A: induction factor
      REAL(DbKi)                 :: AP0   ! Old value for AP: induction factor
      REAL(DbKi)                 :: VNA   ! Effective normal to plane of rotation velocity
      REAL(DbKi)                 :: VTA   ! Effective tangent to plane of rotation velocity
      REAL(DbKi)                 :: CPHI  ! cos(PHI)
      REAL(DbKi)                 :: SPHI  ! sin(PHI)
      REAL(DbKi)                 :: CT    ! elemental thrust coefficient
      REAL(DbKi)                 :: TIPLOSS, HUBLOSS, LOSS     
      REAL(DbKi)                 :: PSI   ! Total azimuthal angle of skew
      REAL(DbKi)                 :: SANG  ! Sin(ANGFLW) (Wind)
      REAL(DbKi)                 :: BB    ! 15p/64* SQRT( (1. - SANG )/(1. + SANG)) 
      
      
      Itera_counter = Itera_counter + 1
            
      IF ( Itera_counter > p%FVM%BEM_Parms%MAX_ITER) THEN
         ErrStat = ErrID_Fatal
         ErrMess = ' Error (in WInDS): BEM for the 1st timestep cannot converge.'
         RETURN
      END IF
      
      ! Save previous values of axial and tangential induction factors
      A0  = A   ! Old value for A: induction factor
      AP0 = AP  ! Old value for AP: induction factor
      
      ! Compute inflow angle and angle of attack
      VNA    = VNWind * ( 1. - A0 ) + VNElement       ! Effective normal to plane of rotation velocity
      VTA    = VTTotal  * ( 1. + AP0 )                ! Effective tangent to plane of rotation velocity

      PHI    = ATAN2( VNA, VTA )                      ! Spanwise inflow Angle
      ALPHA  = PHI - o%Element%PITNOW                 ! AOA

      ALPHA = MODULO( ALPHA, TwoPi )
      IF ( ALPHA > Pi )   ALPHA = ALPHA - TwoPi       ! To ensures that Angle lies between -pi and pi.

      ! Compute lift and drag coefficients
      CALL CLCD_FVM ( P,  O, xd, ErrStat, ErrMess, ALPHA, CLA, CDA, CMA, P%AirFoil%NFoil(J) )
      IF (ErrStat /= 0 ) RETURN

      ! Compute elemental thrust coefficient
      CPHI     = COS( PHI )
      SPHI     = SIN( PHI )
      CT = SIGMAP * (1 - a)**2 * (CLA * CPHI + CDA * SPHI) / SPHI ** 2
            
            
      ! Compute loss correction factor due to tip and hub losses       
      TIPLOSS = 2 / Pi * ACOS(EXP(-(p%NumBl * (p%Blade%TipRadius - p%BLADE%RNodes(J)) / (2 * p%BLADE%RNodes(J) * SPHI))))      
      HUBLOSS = 2 / Pi * ACOS(EXP(-(p%NumBl * (p%BLADE%RNodes(J) - p%Blade%HubRadius) / (2 * p%Blade%HubRadius * SPHI))))           
      LOSS    = TIPLOSS * HUBLOSS       
            
            
      ! Compute axial induction factor using conventional BEM theory
      W2   = VNA * VNA + VTA * VTA          ! W2: Relative velocity squared over blade element
      SPHI = VNA/SQRT( W2 )                 ! SIN(phi)    
      CPhi = COS( Phi )                     ! Cos(PHI)     
      A = REAL ((1 + 4 * LOSS * SPHI **2 / (sigmap * ( CLA * CPhi + CDA * SPHI))) **(-1)) 
      
  
      ! Compute axial induction factor using modified Glauert correction on highly loaded gridpoints
      IF ( CT > 0.96 * LOSS ) THEN
          IF (DIFF_A > p%FVM%BEM_Parms%TOL .OR. DIFF_AP > p%FVM%BEM_Parms%TOL) THEN              
             A = REAL((18 * LOSS - 20 - 3 * SQRT (CT * (50 - 36 * LOSS) + 12 * LOSS * (3 * LOSS - 4))) / (36 * LOSS - 50))
          END IF
      END IF
      
      ! Compute tangential induction factor
       AP = (4 * LOSS * CPHI * SPHI / (SIGMAP * (CLA * SPHI - CDA * CPHI)) - 1) ** (-1)
      
      ! Compute residuals
      DIFF_A  = ABS(A0 - A)
      DIFF_AP = ABS(AP0 - AP)
    
      ! Apply corrective weighting for convergence stability
      IF (p%FVM%BEM_Parms%WT > 0) THEN
        A  = A0 + p%FVM%BEM_Parms%WT * (A - A0)
        AP = AP0 + p%FVM%BEM_Parms%WT * (AP - AP0)
      END IF
      
      ! Decide whether to continue iteration.
      IF ( DIFF_A < p%FVM%BEM_Parms%TOL) THEN
         IF ( DIFF_AP < p%FVM%BEM_Parms%TOL) THEN
             Iter_flag = .FALSE.
         END IF
      END IF
      
   END SUBROUTINE BEM_iteration
   !========================================================   
    
      
END SUBROUTINE BEM

!==================================================================================================================================
SUBROUTINE CLCD_FVM( P,  O, xd,  ErrStat, ErrMess, ALPHA, CLA, CDA, CMA, I)    
! (From AeroSubs.f90 of AeroDyn module)
!...........................................................................
!   This subroutine interpolates airfoil coefficients from a table of airfoil data.  The table must consist
!   of ALPHA, CL and CD over the entire range of angles that will be encountered.
  !
 ! VARIABLES:
 !    CLA      = Returned value of lift coefficient
 !    CDA      = Returned value of drag coeff
 !    CMA      = Returned value of pitching moment coeff
 !    ALPHA    = Angle of attack (radians)
 !    AL       = Array containing the angle of attack
 !    CL       = Array containing the lift coeffs. at AL(I)
 !    CD       = Array containing the drag coeffs. at AL(I)
 !    CM       = Array containing the moment coeffs. at AL(I)
 !    I        = Airfoil ID for this element, equal to NFoil(J), where J is the index identifying the blade element
 !    MulTabLoc= Multiple airfoil table location for this element
 !    MulTabMet= Array containing the multiple airfoil table metric
 ! ******************************************************
!USE                           Airfoil
   IMPLICIT                      NONE
   
   
      ! Passed Variables:
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O!therState ! Initial other/optimization states
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t   
   INTEGER,                      INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMess

   ! Passed Variables:
   REAL(DbKi),INTENT(INOUT)   :: ALPHA
   REAL(DbKi),INTENT(OUT)     :: CDA
   REAL(DbKi),INTENT(OUT)     :: CLA
   REAL(DbKi),INTENT(OUT)     :: CMA

   INTEGER   ,INTENT(IN)      :: I      ! NFOIL(J)

   ! Local Variables:

   REAL(DbKi)                 :: CDA1
   REAL(DbKi)                 :: CDA2
   REAL(DbKi)                 :: CLA1
   REAL(DbKi)                 :: CLA2
   REAL(DbKi)                 :: CMA1
   REAL(DbKi)                 :: CMA2
   REAL(DbKi)                 :: P1
   REAL(DbKi)                 :: P2

   INTEGER                    :: N1
   INTEGER                    :: N1P1
   INTEGER                    :: N2
   INTEGER                    :: N2P1
   INTEGER                    :: NTAB

   ErrStat = ErrID_None
   ErrMess = ""

   IF (.NOT. ALLOCATED(P%AirFoil%NFoil) ) THEN
      CDA = 0
      CLA = 0
      CMA = 0
      ErrStat = ErrID_Fatal
      RETURN
   ELSE
      ErrStat = ErrID_None
   END IF

   NTAB = P%AirFoil%NLIFT(I)

   IF ( ( ALPHA < O%AirFoil%AL(I,1) ) .OR. ( ALPHA > O%AirFoil%AL(I,NTAB) ) )   THEN
   !bjj: This error message isn't necessarially accurate:
      CALL ProgAbort( ' Angle of attack = '//TRIM(Num2LStr(ALPHA*R2D))// &
                      ' deg is outside data table range. '// & !Blade #'//TRIM(Int2LStr(IBLADE))//&
                      ' Airfoil '//TRIM(Int2LStr(I))//'.' )
   !                   ' element '//TRIM(Int2LStr(J))//'.' )

      ErrStat = ErrID_Fatal
      RETURN
   ENDIF

   ALPHA = MIN( MAX( ALPHA, O%AirFoil%AL(I,1) ), O%AirFoil%AL(I,NTAB) )
   CALL LocateBin_Db (ALPHA, O%AirFoil%AL(I,1:NTAB), N1, NTAB )

   IF (N1 == 0) THEN
      N1   = 1
      N1P1 = 2
      P1   = 0.0
   ELSEIF(N1 == NTAB) THEN
      N1P1 = N1
      N1   = N1 - 1
      P1   = 1.0
   ELSE
      N1P1 = N1 + 1
      P1   = ( ALPHA - O%AirFoil%AL(I, N1) )/( O%AirFoil%AL(I, N1P1) - O%AirFoil%AL(I, N1) )
   END IF




    ! If the element has multiple airfoil tables, do a 2-D linear interpolation
    !  for Cl and CD

   IF (P%AirFoil%NTables(I) > 1) THEN

      O%AirFoil%MulTabLoc = MIN( MAX( O%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1) ), P%AirFoil%MulTabMet(I,P%AirFoil%NTables(I)))
      CALL LocateBin (O%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1:P%AirFoil%NTables(I)),N2,P%AirFoil%NTables(I))

      IF (N2 == 0) THEN
         N2   = 1
         N2P1 = 2
         P2   = 0.0
      ELSE IF ( N2 == P%AirFoil%NTables(I) ) THEN
         N2P1 = N2
         N2   = N2 - 1
         P2   = 1.0
      ELSE
         N2P1 = N2 + 1
         P2   = (O%AirFoil%MulTabLoc - P%AirFoil%MulTabMet(I,N2))/(P%AirFoil%MulTabMet(I,N2P1)-P%AirFoil%MulTabMet(I,N2))
      END IF

      CLA1 = O%AirFoil%CL(I,N1,N2) + P1 * ( O%AirFoil%CL(I,N1P1,N2) - O%AirFoil%CL(I,N1,N2) )
      CDA1 = O%AirFoil%CD(I,N1,N2) + P1 * ( O%AirFoil%CD(I,N1P1,N2) - O%AirFoil%CD(I,N1,N2) )
      CMA1 = O%AirFoil%CM(I,N1,N2) + P1 * ( O%AirFoil%CM(I,N1P1,N2) - O%AirFoil%CM(I,N1,N2) )

      CLA2 = O%AirFoil%CL(I,N1,N2P1) + P1 * ( O%AirFoil%CL(I,N1P1,N2P1) - O%AirFoil%CL(I,N1,N2P1) )
      CDA2 = O%AirFoil%CD(I,N1,N2P1) + P1 * ( O%AirFoil%CD(I,N1P1,N2P1) - O%AirFoil%CD(I,N1,N2P1) )
      CMA2 = O%AirFoil%CM(I,N1,N2P1) + P1 * ( O%AirFoil%CM(I,N1P1,N2P1) - O%AirFoil%CM(I,N1,N2P1) )

      CLA = CLA1 + P2 * ( CLA2 - CLA1 )
      CDA = CDA1 + P2 * ( CDA2 - CDA1 )
      CMA = CMA1 + P2 * ( CMA2 - CMA1 )

   ELSE

      CLA  = O%AirFoil%CL(I,N1,1) + P1 * ( O%AirFoil%CL(I,N1P1,1) - O%AirFoil%CL(I,N1,1) )
      CDA  = O%AirFoil%CD(I,N1,1) + P1 * ( O%AirFoil%CD(I,N1P1,1) - O%AirFoil%CD(I,N1,1) )
      CMA  = O%AirFoil%CM(I,N1,1) + P1 * ( O%AirFoil%CM(I,N1P1,1) - O%AirFoil%CM(I,N1,1) )

   ENDIF ! (P%AirFoil%NTables(I) > 1)
   
CONTAINS
   ! ====================================================================================================
   SUBROUTINE LocateBin_Db( XVal, XAry, Ind, AryLen )

      ! This subroutine finds the lower-bound index of an input x-value located in an array.
      ! On return, Ind has a value such that
      !           XAry(Ind) <= XVal < XAry(Ind+1), with the exceptions that
      !             Ind = 0 when XVal < XAry(1), and
      !          Ind = AryLen when XAry(AryLen) <= XVal.
      !
      ! It uses a binary interpolation scheme that takes about log(AryLen)/log(2) steps to converge.
      ! If the index doesn't change much between calls, LocateStp() may be a better option.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(OUT)         :: Ind                                             ! Final (low) index into the array.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(DbKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.


      ! Local declarations.

   INTEGER                      :: IHi                                             ! The high index into the arrays.
   INTEGER                      :: IMid                                            ! The mid-point index between IHi and Ind.



      ! Let's check the limits first.

   IF ( XVal < XAry(1) )  THEN
      Ind = 0
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      Ind = AryLen
   ELSE
         ! Let's interpolate!

      Ind  = 1
      IHi  = AryLen

      DO WHILE ( IHi-Ind > 1 )

         IMid = ( IHi + Ind )/2

         IF ( XVal >= XAry(IMid) ) THEN
            Ind = IMid
         ELSE
            IHi = IMid
         END IF

      END DO

   END IF

   RETURN
   END SUBROUTINE LocateBin_Db
   ! ====================================================================================================

END SUBROUTINE CLCD_FVM

!====================================================================================================

    
END MODULE WINDS
!**********************************************************************************************************************************