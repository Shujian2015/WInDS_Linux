MODULE WINDS_DS
   ! Reference: 
   !     Gaertner, Evan M., "Modeling Dynamic Stall for a Free Vortex Wake Model of a Floating Offshore Wind Turbine" (2014).
   !     Masters Theses May 2014-current. Paper 85.
   !     http://scholarworks.umass.edu/masters_theses_2/85/
    
   USE NWTC_Library 
   USE AeroDyn_Types
   USE WINDS_IO          ! Handle input and output
   !USE WINDS_Library    ! Some subroutines for debug

   
   ! Sliu: problem for NTables
   ! All angles are in radius while in Matlab code angles are in degrees.
   
   
   
   IMPLICIT        NONE

      ! ..... Public Subroutines ............
      
   PUBLIC :: LB_Initialize_Variables
   PUBLIC :: LB_Initialize_AirfoilData
   PUBLIC :: LB_AfData_spline_fit ! New
   PUBLIC :: LB_DynStall   


CONTAINS
!==================================================================================================================================
SUBROUTINE LB_Initialize_Variables(p, O, xd, ErrStat, ErrMess)
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
   INTEGER(IntKi)      :: NB 
   INTEGER(IntKi)      :: NF
   INTEGER(IntKi)      :: NumCl ! the maximum number of lines (for only the data such as AOA or Cl) across all files
   
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""
   
      
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NB  = p%NumBl
   NF  = P%AirFoil%NumFoil
   NumCl = P%AirFoil%NumCl ! the maximum number of lines (for only the data such as AOA or Cl) across all files

   
   
   
   
     !------------------------------------------------------------------
     ! These are used for Dynamic Stall:
     !....................................  
   
   ! Recursive Variables
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_previous_alpha ) )         ALLOCATE( O%FVM_Other%DS%global_previous_alpha(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_current_alpha ) )          ALLOCATE( O%FVM_Other%DS%global_current_alpha(NS,  NB))    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_previous_sigma1 ) )        ALLOCATE( O%FVM_Other%DS%global_previous_sigma1(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_current_sigma1 ) )         ALLOCATE( O%FVM_Other%DS%global_current_sigma1(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_previous_sigma3 ) )        ALLOCATE( O%FVM_Other%DS%global_previous_sigma3(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_current_sigma3 ) )         ALLOCATE( O%FVM_Other%DS%global_current_sigma3(NS,  NB))    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_q ) )           ALLOCATE( O%FVM_Other%DS%attached_previous_q(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_q ) )            ALLOCATE( O%FVM_Other%DS%attached_current_q(NS,  NB))    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_X1 ) )          ALLOCATE( O%FVM_Other%DS%attached_previous_X1(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_X1 ) )           ALLOCATE( O%FVM_Other%DS%attached_current_X1(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_X2 ) )          ALLOCATE( O%FVM_Other%DS%attached_previous_X2(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_X2 ) )           ALLOCATE( O%FVM_Other%DS%attached_current_X2(NS,  NB))   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_K_alpha ) )     ALLOCATE( O%FVM_Other%DS%attached_previous_K_alpha(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_K_alpha ) )      ALLOCATE( O%FVM_Other%DS%attached_current_K_alpha(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_dK_alpha ) )    ALLOCATE( O%FVM_Other%DS%attached_previous_dK_alpha(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_dK_alpha ) )     ALLOCATE( O%FVM_Other%DS%attached_current_dK_alpha(NS,  NB))   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_K_q ) )         ALLOCATE( O%FVM_Other%DS%attached_previous_K_q(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_K_q ) )          ALLOCATE( O%FVM_Other%DS%attached_current_K_q(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_dK_q ) )        ALLOCATE( O%FVM_Other%DS%attached_previous_dK_q(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_dK_q ) )         ALLOCATE( O%FVM_Other%DS%attached_current_dK_q(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_ddK_q ) )        ALLOCATE( O%FVM_Other%DS%attached_previous_ddK_q(NS,  NB))    
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_ddK_q ) )         ALLOCATE( O%FVM_Other%DS%attached_current_ddK_q(NS,  NB))
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_previous_C_p_n ) )       ALLOCATE( O%FVM_Other%DS%attached_previous_C_p_n(NS,  NB))        
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%attached_current_C_p_n ) )        ALLOCATE( O%FVM_Other%DS%attached_current_C_p_n(NS,  NB))    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_previous_D_p ) )            ALLOCATE( O%FVM_Other%DS%TEsep_previous_D_p(NS,  NB))     
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_current_D_p ) )             ALLOCATE( O%FVM_Other%DS%TEsep_current_D_p(NS,  NB))       
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_previous_D_f ) )            ALLOCATE( O%FVM_Other%DS%TEsep_previous_D_f(NS,  NB))        
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_current_D_f ) )             ALLOCATE( O%FVM_Other%DS%TEsep_current_D_f(NS,  NB))  
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_previous_f_prime ) )        ALLOCATE( O%FVM_Other%DS%TEsep_previous_f_prime(NS,  NB))        
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_current_f_prime ) )         ALLOCATE( O%FVM_Other%DS%TEsep_current_f_prime(NS,  NB))        
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_previous_f_2prime ) )       ALLOCATE( O%FVM_Other%DS%TEsep_previous_f_2prime(NS,  NB))     
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%TEsep_current_f_2prime ) )        ALLOCATE( O%FVM_Other%DS%TEsep_current_f_2prime(NS,  NB))  
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_previous_tau_v ) )          ALLOCATE( O%FVM_Other%DS%LEsep_previous_tau_v(NS,  NB))      
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_current_tau_v ) )           ALLOCATE( O%FVM_Other%DS%LEsep_current_tau_v(NS,  NB))     
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_previous_C_v_n ) )          ALLOCATE( O%FVM_Other%DS%LEsep_previous_C_v_n(NS,  NB))        
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_current_C_v_n ) )           ALLOCATE( O%FVM_Other%DS%LEsep_current_C_v_n(NS,  NB))    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_previous_C_v ) )            ALLOCATE( O%FVM_Other%DS%LEsep_previous_C_v(NS,  NB))     
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_current_C_v ) )             ALLOCATE( O%FVM_Other%DS%LEsep_current_C_v(NS,  NB))      
   

   O%FVM_Other%DS%global_previous_alpha  =   0.0_DbKi 
   O%FVM_Other%DS%global_current_alpha   =   0.0_DbKi 
   
   O%FVM_Other%DS%global_previous_sigma1 =   1.0_DbKi 
   O%FVM_Other%DS%global_current_sigma1  =   1.0_DbKi 
   O%FVM_Other%DS%global_previous_sigma3 =   1.0_DbKi     
   O%FVM_Other%DS%global_current_sigma3  =   1.0_DbKi  
   
   O%FVM_Other%DS%attached_previous_q    =   0.0_DbKi 
   O%FVM_Other%DS%attached_current_q     =   0.0_DbKi 
   
   O%FVM_Other%DS%attached_previous_X1   =   0.0_DbKi 
   O%FVM_Other%DS%attached_current_X1    =   0.0_DbKi 
   O%FVM_Other%DS%attached_previous_X2   =   0.0_DbKi 
   O%FVM_Other%DS%attached_current_X2    =   0.0_DbKi 
   
   O%FVM_Other%DS%attached_previous_K_alpha   =  0.0_DbKi   
   O%FVM_Other%DS%attached_current_K_alpha    =  0.0_DbKi 
   O%FVM_Other%DS%attached_previous_dK_alpha  =  0.0_DbKi   
   O%FVM_Other%DS%attached_current_dK_alpha   =  0.0_DbKi 
   
   O%FVM_Other%DS%attached_previous_K_q   =   0.0_DbKi 
   O%FVM_Other%DS%attached_current_K_q    =   0.0_DbKi 
   O%FVM_Other%DS%attached_previous_dK_q  =   0.0_DbKi  
   O%FVM_Other%DS%attached_current_dK_q   =   0.0_DbKi 
   O%FVM_Other%DS%attached_previous_ddK_q  =   0.0_DbKi  
   O%FVM_Other%DS%attached_current_ddK_q   =   0.0_DbKi   
   
   O%FVM_Other%DS%attached_previous_C_p_n =   0.0_DbKi    
   O%FVM_Other%DS%attached_current_C_p_n  =   0.0_DbKi 
   
   O%FVM_Other%DS%TEsep_previous_D_p      =   0.0_DbKi 
   O%FVM_Other%DS%TEsep_current_D_p       =   0.0_DbKi 
   
   O%FVM_Other%DS%TEsep_previous_D_f      =   0.0_DbKi    
   O%FVM_Other%DS%TEsep_current_D_f       =   0.0_DbKi 
   
   O%FVM_Other%DS%TEsep_previous_f_prime  =   0.0_DbKi    
   O%FVM_Other%DS%TEsep_current_f_prime   =   0.0_DbKi  
   O%FVM_Other%DS%TEsep_previous_f_2prime =   0.0_DbKi 
   O%FVM_Other%DS%TEsep_current_f_2prime  =   0.0_DbKi 
   
   O%FVM_Other%DS%LEsep_previous_tau_v   =    0.0_DbKi 
   O%FVM_Other%DS%LEsep_current_tau_v    =    0.0_DbKi 
   
   O%FVM_Other%DS%LEsep_previous_C_v_n   =    0.0_DbKi     
   O%FVM_Other%DS%LEsep_current_C_v_n    =    0.0_DbKi 
   
   O%FVM_Other%DS%LEsep_previous_C_v     =    0.0_DbKi  
   O%FVM_Other%DS%LEsep_current_C_v      =    0.0_DbKi 
   
   
   ! Non-periodically updated variables 
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_sigma2 ) )                 ALLOCATE( O%FVM_Other%DS%global_sigma2(NS,  NB))      
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_vortex_shed_start ) )       ALLOCATE( O%FVM_Other%DS%LEsep_vortex_shed_start(NS,  NB))    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%LEsep_seperating ) )              ALLOCATE( O%FVM_Other%DS%LEsep_seperating(NS,  NB))     
   
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_ca_mdl_change ) )          ALLOCATE( O%FVM_Other%DS%global_ca_mdl_change(NS,  NB))        
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_ca_mdl ) )                 ALLOCATE( O%FVM_Other%DS%global_ca_mdl(NS,  NB, p%FVM%MAXITER))      
   IF (.NOT. ALLOCATED( O%FVM_Other%DS%global_ca_mdl_over_write ) )      ALLOCATE( O%FVM_Other%DS%global_ca_mdl_over_write(NS,  NB))      

   O%FVM_Other%DS%global_sigma2             = 5.0_DbKi 
   O%FVM_Other%DS%LEsep_vortex_shed_start   = 0_IntKi 
   O%FVM_Other%DS%LEsep_seperating          = 0_IntKi
   O%FVM_Other%DS%global_ca_mdl_change      = 0_IntKi
   O%FVM_Other%DS%global_ca_mdl             = 0_IntKi
   O%FVM_Other%DS%global_ca_mdl_over_write  = 0_IntKi
   

     !------------------------------------------------------------------
     ! These are used for LB airfoils:
     !....................................  
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%A1 ) )         ALLOCATE( O%FVM_Other%airfoils_LB%A1(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%A2 ) )         ALLOCATE( O%FVM_Other%airfoils_LB%A2(NF)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%b1 ) )         ALLOCATE( O%FVM_Other%airfoils_LB%b1(NF)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%b2 ) )         ALLOCATE( O%FVM_Other%airfoils_LB%b2(NF)) 
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%Tp ) )         ALLOCATE( O%FVM_Other%airfoils_LB%Tp(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%Tf ) )         ALLOCATE( O%FVM_Other%airfoils_LB%Tf(NF)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%Tv ) )         ALLOCATE( O%FVM_Other%airfoils_LB%Tv(NF)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%Tvl ) )        ALLOCATE( O%FVM_Other%airfoils_LB%Tvl(NF)) 
   
   O%FVM_Other%airfoils_LB%A1  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%A2  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%b1  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%b2  = 0.0_DbKi 
   
   O%FVM_Other%airfoils_LB%Tp  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%Tf  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%Tv  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%Tvl = 0.0_DbKi    
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%cn ) )         ALLOCATE( O%FVM_Other%airfoils_LB%cn(NF, p%AirFoil%NumCL))    ! NumCL: the maximum number of lines (for only the data such as AOA or Cl) across all files
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%ca ) )         ALLOCATE( O%FVM_Other%airfoils_LB%ca(NF, p%AirFoil%NumCL)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%Cn_alpha ) )   ALLOCATE( O%FVM_Other%airfoils_LB%Cn_alpha(NF)) 
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%Alpha_0 ) )    ALLOCATE( O%FVM_Other%airfoils_LB%Alpha_0(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%C_d_0 ) )      ALLOCATE( O%FVM_Other%airfoils_LB%C_d_0(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%C_m_0 ) )      ALLOCATE( O%FVM_Other%airfoils_LB%C_m_0(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%C_n_0 ) )      ALLOCATE( O%FVM_Other%airfoils_LB%C_n_0(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%alpha_1 ) )    ALLOCATE( O%FVM_Other%airfoils_LB%alpha_1(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%alpha_2 ) )    ALLOCATE( O%FVM_Other%airfoils_LB%alpha_2(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%a ) )          ALLOCATE( O%FVM_Other%airfoils_LB%a(NF,3))  
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%S ) )          ALLOCATE( O%FVM_Other%airfoils_LB%S(NF,3))     
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%c ) )          ALLOCATE( O%FVM_Other%airfoils_LB%c(NF,3))     
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%recovery_factor ) )      ALLOCATE( O%FVM_Other%airfoils_LB%recovery_factor(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%C_n_2 ) )      ALLOCATE( O%FVM_Other%airfoils_LB%C_n_2(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%C_n_1 ) )      ALLOCATE( O%FVM_Other%airfoils_LB%C_n_1(NF))    
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%ca_K1 ) )      ALLOCATE( O%FVM_Other%airfoils_LB%ca_K1(NF))    
   
   
   O%FVM_Other%airfoils_LB%cn       = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%ca       = 0.0_DbKi    
   O%FVM_Other%airfoils_LB%Cn_alpha  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%Alpha_0  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%C_d_0    = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%C_m_0    = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%C_n_0    = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%alpha_1  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%alpha_2  = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%a        = 0.0_DbKi    
   O%FVM_Other%airfoils_LB%S        = 0.0_DbKi    
   O%FVM_Other%airfoils_LB%c        = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%recovery_factor = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%C_n_2    = 0.0_DbKi 
   O%FVM_Other%airfoils_LB%C_n_1    = 0.0_DbKi
   O%FVM_Other%airfoils_LB%ca_K1    = 0.0_DbKi
   
   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%aoa_KHn ) )      ALLOCATE( O%FVM_Other%airfoils_LB%aoa_KHn(NF))
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%aoa_KHp ) )      ALLOCATE( O%FVM_Other%airfoils_LB%aoa_KHp(NF))
   
   O%FVM_Other%airfoils_LB%aoa_KHn  = 0.0_DbKi
   O%FVM_Other%airfoils_LB%aoa_KHp  = 0.0_DbKi
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_size ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_size(NF, 8))

   O%FVM_Other%airfoils_LB%pp_size = 0_IntKi
   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_f_breaks ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_f_breaks(NF, NumCl))
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_f_coefs ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_f_coefs(NF, NumCl, 4))
   
   O%FVM_Other%airfoils_LB%pp_f_breaks  = 0.0_DbKi
   O%FVM_Other%airfoils_LB%pp_f_coefs  = 0.0_DbKi
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_K_Ca_breaks ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_K_Ca_breaks(NF, NumCl))
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_K_Ca_coefs ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_K_Ca_coefs(NF, NumCl, 4))
   
   O%FVM_Other%airfoils_LB%pp_K_Ca_breaks  = 0.0_DbKi
   O%FVM_Other%airfoils_LB%pp_K_Ca_coefs  = 0.0_DbKi   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_K_Cn_breaks ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_K_Cn_breaks(NF, NumCl))
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_K_Cn_coefs ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_K_Cn_coefs(NF, NumCl, 4))
   
   O%FVM_Other%airfoils_LB%pp_K_Cn_breaks  = 0.0_DbKi
   O%FVM_Other%airfoils_LB%pp_K_Cn_coefs  = 0.0_DbKi   
   
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_K_Cm_breaks ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_K_Cm_breaks(NF, NumCl))
   IF (.NOT. ALLOCATED( O%FVM_Other%airfoils_LB%pp_K_Cm_coefs ) )      ALLOCATE( O%FVM_Other%airfoils_LB%pp_K_Cm_coefs(NF, NumCl, 4))
   
   O%FVM_Other%airfoils_LB%pp_K_Cm_breaks  = 0.0_DbKi
   O%FVM_Other%airfoils_LB%pp_K_Cm_coefs  = 0.0_DbKi   
   
   
   

END SUBROUTINE LB_Initialize_Variables     
!==================================================================================================================================     
SUBROUTINE LB_Initialize_AirfoilData(p, O, xd, ErrStat, ErrMess)
! Compute additional airfoil data needed for the dynamic stall                                                                                        -
! model on initialization.  Also call a subroutine to perform curve                                                                                               -
! fitting to the 2D static airfoil data
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
   INTEGER(IntKi)      :: NB 
   
   INTEGER(IntKi)      :: naf
   INTEGER(IntKi)      :: NFOILID, J, K, M
   INTEGER(IntKi)      :: nAoA 
   REAL(DbKi), dimension(P%AirFoil%NumCL)   :: Cn_temp
   INTEGER(IntKi)      :: Check_name
   
   !REAL(DbKi)          :: StallAoA
      
   REAL(DbKi)          :: aoa_Cn2
   REAL(DbKi)          :: C_n_1_v1 
   REAL(DbKi)          :: C_n_1_v2 
   REAL(DbKi)          :: temp_angle
   
   CHARACTER(LEN = 1024)   :: TEMP1  ! debug
      
      
      
   !StallAoA = pi/12  ! Sliu: Not correct
   
   NST = p%Element%NElm + 1
   NS  = p%Element%NElm
   NB  = p%NumBl
   

   DO NFOILID = 1, P%AirFoil%NumFoil
      IF (P%AirFoil%NTables(NFOILID) > 1) THEN 
         ErrMess = ' WInDS error: Currently WInDS does not support multi-column airfoil file.'
         ErrStat = ErrID_Fatal
         RETURN   
      ELSE  
          

         Check_name =  INDEX((p%AirFoil%FOILNM(NFOILID)), "Cylinder") 

         IF (Check_name == 0) THEN ! Not "Cylinder"
            O%FVM_Other%airfoils_LB%A1(NFOILID)  = p%FVM%DS_Parms%indicial(1)
            O%FVM_Other%airfoils_LB%A2(NFOILID)  = p%FVM%DS_Parms%indicial(2) 
            O%FVM_Other%airfoils_LB%b1(NFOILID)  = p%FVM%DS_Parms%indicial(3) 
            O%FVM_Other%airfoils_LB%b2(NFOILID)  = p%FVM%DS_Parms%indicial(4) 
   
            O%FVM_Other%airfoils_LB%Tp(NFOILID)  = p%FVM%DS_Parms%time_const(1)
            O%FVM_Other%airfoils_LB%Tf(NFOILID)  = p%FVM%DS_Parms%time_const(2)
            O%FVM_Other%airfoils_LB%Tv(NFOILID)  = p%FVM%DS_Parms%time_const(3) 
            O%FVM_Other%airfoils_LB%Tvl(NFOILID) = p%FVM%DS_Parms%time_const(4)
         
           
            ! Calc Static Profiles for Cn and Ca from Cl and Cd
            nAoA = P%AirFoil%NLIFT(NFOILID)  ! Number of AoA lines stored in each airfoil file
            
            DO J = 1, nAoA-1
               IF ( EqualRealNos( O%AirFoil%AL(NFOILID,J), O%AirFoil%AL(NFOILID,J+1) ) ) THEN
                   ErrMess = ' WInDS error: There are two lines in one airfoil file with same angle. Please delete one.' 
                   ErrStat = ErrID_Fatal
                   RETURN       
               END IF               
            END DO ! j       
            
                     
            DO J = 1, nAoA
               O%FVM_Other%airfoils_LB%cn(NFOILID, J) = O%AirFoil%CL(NFOILID,J,1) * COS(O%AirFoil%AL(NFOILID,J)) & 
                                                        + O%AirFoil%CD(NFOILID,J,1) * SIN(O%AirFoil%AL(NFOILID,J))  ! O%AirFoil%AL in radians
               O%FVM_Other%airfoils_LB%ca(NFOILID, J) = - O%AirFoil%CD(NFOILID,J,1) * COS(O%AirFoil%AL(NFOILID,J)) & 
                                                        + O%AirFoil%CL(NFOILID,J,1) * SIN(O%AirFoil%AL(NFOILID,J))
            END DO ! j
            
            
            ! Cn-alpha slope
            Cn_temp = 0.0_DbKi
            
            DO J = 2, nAoA
               IF (O%AirFoil%AL(NFOILID,J) >= O%Beddoes%AOL( NFOILID, 1 ) .AND. O%AirFoil%AL(NFOILID,J) < O%Beddoes%StallAngle( NFOILID, 1 )) THEN ! StallAoA) THEN 
                  Cn_temp(J) = (O%FVM_Other%airfoils_LB%cn(NFOILID, J) - O%FVM_Other%airfoils_LB%cn(NFOILID, J-1)) / &
                               ( O%AirFoil%AL(NFOILID,J) -  O%AirFoil%AL(NFOILID,J-1) )                       
               END IF
            END DO ! j     
           
  
            O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID) = MAXVAL(Cn_temp(1:nAoA))

    
             ! Find zero crossing points  
            DO J = 2, nAoA
               IF (O%AirFoil%AL(NFOILID,J) >= -(pi/12) .AND. O%AirFoil%AL(NFOILID,J) <= O%Beddoes%StallAngle( NFOILID, 1 )) THEN !StallAoA ) THEN 
                  IF (O%FVM_Other%airfoils_LB%cn(NFOILID, J-1) < 0.0 .AND. O%FVM_Other%airfoils_LB%cn(NFOILID, J) >= 0.0) THEN
                     O%FVM_Other%airfoils_LB%C_d_0(NFOILID) = interp1D_1(O%FVM_Other%airfoils_LB%cn(NFOILID, J-1:J), O%AirFoil%CD(NFOILID, J-1:J, 1), 0.0_DbKi)
                     O%FVM_Other%airfoils_LB%C_m_0(NFOILID) = interp1D_1(O%FVM_Other%airfoils_LB%cn(NFOILID, J-1:J), O%AirFoil%CM(NFOILID, J-1:J, 1), 0.0_DbKi) 
                  END IF
                  IF (O%AirFoil%AL(NFOILID, J-1) < 0.0 .AND. O%AirFoil%AL(NFOILID, J) >= 0.0 )  THEN
                     O%FVM_Other%airfoils_LB%C_n_0(NFOILID) = interp1D_2(O%AirFoil%AL(NFOILID, J-1:J), O%FVM_Other%airfoils_LB%cn(NFOILID, J-1:J), 0.0_DbKi)                    
                  END IF
               END IF            
            END DO ! j     
            
            O%FVM_Other%airfoils_LB%Alpha_0(NFOILID) = O%Beddoes%AOL( NFOILID, 1 ) ! Both in radius
             
            ! Fits to airfoil data and TE seperation point
            CALL LB_AfData_spline_fit(NFOILID, p, O, xd, ErrStat, ErrMess)
           
            ! Critical normal force coefficient
            aoa_Cn2 = O%FVM_Other%airfoils_LB%alpha_2(NFOILID) 
            O%FVM_Other%airfoils_LB%C_n_2(NFOILID) =  O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID) * aoa_Cn2 + O%FVM_Other%airfoils_LB%C_n_0(NFOILID)
            
            C_n_1_v1 = 0.0_DbKi
            
            DO J = 2, nAoA
               IF (O%AirFoil%AL(NFOILID,J) >= -(pi/12) .AND. O%AirFoil%AL(NFOILID,J) <= O%Beddoes%StallAngle( NFOILID, 1 )) THEN !StallAoA ) THEN 
                  IF (O%FVM_Other%airfoils_LB%cn(NFOILID, J) > C_n_1_v1) THEN
                     C_n_1_v1 = O%FVM_Other%airfoils_LB%cn(NFOILID, J)
                  END IF
               END IF            
            END DO ! j     
                       
            C_n_1_v1 = C_n_1_v1 + O%FVM_Other%airfoils_LB%C_n_0(NFOILID)
            
            temp_angle = O%FVM_Other%airfoils_LB%alpha_1(NFOILID)
            DO J=1, nAoA-1
               IF (temp_angle >= O%AirFoil%AL(NFOILID,J) .AND. temp_angle <= O%AirFoil%AL(NFOILID,J+1)) THEN
                  K=J
                  EXIT
               END IF   
            END DO
            
            C_n_1_v2 = interp1D_2(O%AirFoil%AL(NFOILID, K:K+1), O%FVM_Other%airfoils_LB%cn(NFOILID, K:K+1), temp_angle)
            C_n_1_v2 = C_n_1_v2 + O%FVM_Other%airfoils_LB%C_n_0(NFOILID)
            
            O%FVM_Other%airfoils_LB%C_n_1(NFOILID) = MAX(C_n_1_v1, C_n_1_v2)
                      
         END IF ! NOT "Cylinder"   
      END IF ! (P%AirFoil%NTables(NFOILID) > 1) 
   END DO !  
   
   

   

END SUBROUTINE LB_Initialize_AirfoilData     
!==================================================================================================================================  
SUBROUTINE LB_AfData_spline_fit(NFOILID, p, O, xd, ErrStat, ErrMess)
 

   IMPLICIT                        NONE

      ! Passed variables

   INTEGER(IntKi),              INTENT(IN   )   :: NFOILID        ! Counter for airfoil file 
   TYPE(AD_ParameterType),      INTENT(IN   )   :: p              ! The module's parameter data
   TYPE(AD_OtherStateType),     INTENT(INOUT)   :: O              ! Other/optimization states   
   TYPE(AD_DiscreteStateType),  INTENT(IN   )   :: xd             ! Discrete states at t   
   INTEGER(IntKi),              INTENT(  OUT)   :: ErrStat        ! The error status code
   CHARACTER(*),                INTENT(  OUT)   :: ErrMess        ! The error message, if an error occurred

   
    
   

   !   ! Local variables
   REAL(DbKi)      :: aoa_start 
   REAL(DbKi)      :: aoa_end
   INTEGER(IntKi)  :: n_rg 
   INTEGER(IntKi)  :: I , J, K   
   REAL(DbKi)      :: cn_temp 
   REAL(DbKi)      :: Cn_alpha_temp 
   REAL(DbKi)      :: aoa_temp
   REAL(DbKi)      :: Alpha_0_temp
   REAL(DbKi)      :: temp    
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: AOA
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cn
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: ca
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cl   
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cd
   REAL(DbKi)  :: Cn_alpha
   REAL(DbKi)  :: Alpha_0
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cm
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: f
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: Cn_f
   INTEGER(IntKi)  :: idx_DS_p   
   INTEGER(IntKi)  :: idx_DS_n
   INTEGER(IntKi), DIMENSION(:), ALLOCATABLE :: rg_KH
   INTEGER(IntKi), DIMENSION(:), ALLOCATABLE :: rg_DS_p
   INTEGER(IntKi), DIMENSION(:), ALLOCATABLE :: rg_DS_n
   
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: f1 
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: aoa1
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: ca1 
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cd1 
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cn1 
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cm1 

   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: aoa2
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: ca2
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cn2
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cm2

   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: aoa3 
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: ca3
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cn3
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: cm3
   
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: aoa_f
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: f_temp
   
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: temp_coefs
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: f_mdl_aoa
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Ca_mdl1_aoa
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Ca_mdl1_err
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_KH_Ca
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Ca_mdl2p_aoa
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Ca_mdl2p_err
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_DSp_Ca
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Ca_mdl2n_aoa
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Ca_mdl2n_err
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_DSn_Ca
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cn_mdl2p_aoa
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cn_mdl2p_err
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_DSp_Cn
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cn_mdl2n_aoa
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cn_mdl2n_err
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_DSn_Cn  
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cm_mdl1_aoa 
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cm_mdl1_err 
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: xcp_mdl1_aoa
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_KH_Cm
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cm_mdl2p_err
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_DSp_Cm    
   REAL(DbKi), DIMENSION(:),   ALLOCATABLE :: Cm_mdl2n_err
   REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: pp_K_DSn_Cm  
   
   INTEGER(IntKi)  :: curve_fit_limit1
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: dca
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: aoa_m   
   REAL(DbKi)  :: ddca_min
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: ddca
   INTEGER(IntKi)  :: curve_fit_limit2
   !
   !CHARACTER(LEN = 1024)   :: TEMP1  ! debug
   
   
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""
   
      
   
   ! Set bounds on airfoil data
   aoa_start = - pi
   aoa_end = pi   
   n_rg = 0_IntKi
   
   DO I = 1, P%AirFoil%NLIFT(NFOILID)
      IF (O%AirFoil%AL(NFOILID,I) >= aoa_start .AND. O%AirFoil%AL(NFOILID,I) <= aoa_end) THEN
         cn_temp = O%FVM_Other%airfoils_LB%cn(NFOILID, I)
         Cn_alpha_temp = O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID)  ! in radians
         aoa_temp = O%AirFoil%AL(NFOILID,I)                         ! in radians
         Alpha_0_temp = O%FVM_Other%airfoils_LB%Alpha_0(NFOILID)    ! in radians
         
         temp = cn_temp / (Cn_alpha_temp *( aoa_temp - Alpha_0_temp ))
         
         IF ( temp > 0.0_DbKi ) THEN
            n_rg = n_rg +1  
         END IF
         
      END IF
   END DO

    
   IF (.NOT. ALLOCATED( AOA ) )      ALLOCATE( AOA(n_rg))    
   IF (.NOT. ALLOCATED( cn  ) )      ALLOCATE( cn(n_rg) )    
   IF (.NOT. ALLOCATED( ca ) )       ALLOCATE( ca(n_rg) )    
   IF (.NOT. ALLOCATED( cl ) )       ALLOCATE( cl(n_rg) )    
   IF (.NOT. ALLOCATED( cd ) )       ALLOCATE( cd(n_rg) ) 
   IF (.NOT. ALLOCATED( cm ) )       ALLOCATE( cm(n_rg) )
   IF (.NOT. ALLOCATED( f ) )        ALLOCATE( f(n_rg)  )    
   IF (.NOT. ALLOCATED( Cn_f ) )     ALLOCATE( Cn_f(n_rg)  ) 
   
   K = 1_IntKi
   
   DO I = 1, P%AirFoil%NLIFT(NFOILID)
      IF (O%AirFoil%AL(NFOILID,I) >= aoa_start .AND. O%AirFoil%AL(NFOILID,I) <= aoa_end) THEN
         cn_temp = O%FVM_Other%airfoils_LB%cn(NFOILID, I)
         Cn_alpha_temp = O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID)
         aoa_temp = O%AirFoil%AL(NFOILID,I)
         Alpha_0_temp = O%FVM_Other%airfoils_LB%Alpha_0(NFOILID)
         
         temp = cn_temp / (Cn_alpha_temp *( aoa_temp - Alpha_0_temp ))
         
         IF ( temp > 0.0_DbKi ) THEN             
            aoa(K) = O%AirFoil%AL(NFOILID,I)
            cn(K)  = O%FVM_Other%airfoils_LB%cn(NFOILID, I)
            ca(K)  = O%FVM_Other%airfoils_LB%ca(NFOILID, I)
            cl(K) = O%AirFoil%CL(NFOILID,I,1)
            cd(K) = O%AirFoil%CD(NFOILID,I,1)
            cm(K) = O%AirFoil%CM(NFOILID,I,1)               
            K = K +1  
         END IF         
      END IF
   END DO

      
   Cn_alpha = O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID)  ! in radians
   Alpha_0  = O%FVM_Other%airfoils_LB%Alpha_0(NFOILID)   ! in radians
   
   
    ! Find deep stall break point, where trailing edge separation point equation diverges from the data  
   DO I = 1, n_rg 
      f(I) = (2*SQRT(cn(I)/( Cn_alpha*(aoa(I)-Alpha_0)))-1)**2   
      Cn_f(I) = Cn_alpha * (aoa(I)-Alpha_0) * ((1+SQRT(f(I)))/2)**2      
   END DO   
   
   DO I = 1, n_rg
      IF ( (Cn_f(I)-cn(I))> 1e-10 .AND. aoa(I) >0.0) THEN
          idx_DS_p = I -1
          EXIT          
      END IF   
   END DO   
   
   DO I = 1, n_rg 
      IF ( ABS(Cn_f(I)-cn(I)) < 1e-10 ) THEN
          idx_DS_n = I
          EXIT          
      END IF   
   END DO      
   
   O%FVM_Other%airfoils_LB%aoa_KHn(NFOILID) = aoa(idx_DS_n)
   O%FVM_Other%airfoils_LB%aoa_KHp(NFOILID) = aoa(idx_DS_p)
   
    ! Separate data for fits in different regions
   IF (.NOT. ALLOCATED( rg_KH ) )      ALLOCATE( rg_KH(idx_DS_p - idx_DS_n + 1))
   IF (.NOT. ALLOCATED( rg_DS_p ) )    ALLOCATE( rg_DS_p(n_rg - idx_DS_p + 1))
   IF (.NOT. ALLOCATED( rg_DS_n ) )    ALLOCATE( rg_DS_n(idx_DS_n))
   
   IF (.NOT. ALLOCATED( f1 ) )      ALLOCATE( f1(idx_DS_p - idx_DS_n + 1))
   IF (.NOT. ALLOCATED( aoa1 ) )    ALLOCATE( aoa1(idx_DS_p - idx_DS_n + 1))
   IF (.NOT. ALLOCATED( ca1 ) )     ALLOCATE( ca1(idx_DS_p - idx_DS_n + 1))
   IF (.NOT. ALLOCATED( cd1 ) )     ALLOCATE( cd1(idx_DS_p - idx_DS_n + 1))
   IF (.NOT. ALLOCATED( cn1 ) )     ALLOCATE( cn1(idx_DS_p - idx_DS_n + 1))
   IF (.NOT. ALLOCATED( cm1 ) )     ALLOCATE( cm1(idx_DS_p - idx_DS_n + 1))
   
   IF (.NOT. ALLOCATED( aoa2 ) )    ALLOCATE( aoa2(n_rg - idx_DS_p + 1))   
   IF (.NOT. ALLOCATED( ca2 ) )    ALLOCATE( ca2(n_rg - idx_DS_p + 1))   
   IF (.NOT. ALLOCATED( cn2 ) )    ALLOCATE( cn2(n_rg - idx_DS_p + 1))   
   IF (.NOT. ALLOCATED( cm2 ) )    ALLOCATE( cm2(n_rg - idx_DS_p + 1))   

   IF (.NOT. ALLOCATED( aoa3 ) )    ALLOCATE( aoa3(idx_DS_n))   
   IF (.NOT. ALLOCATED( ca3 ) )    ALLOCATE( ca3(idx_DS_n))   
   IF (.NOT. ALLOCATED( cn3 ) )    ALLOCATE( cn3(idx_DS_n))   
   IF (.NOT. ALLOCATED( cm3 ) )    ALLOCATE( cm3(idx_DS_n))        
   
   DO I = 1, idx_DS_p - idx_DS_n + 1
      rg_KH(I) = idx_DS_n + I - 1   ! Kirchhoff-Helmholtz flate plate model
      
      f1(I)   = f(rg_KH(I))
      aoa1(I) = aoa(rg_KH(I))
      ca1(I)  = ca(rg_KH(I))
      cd1(I)  = cd(rg_KH(I))
      cn1(I)  = cn(rg_KH(I))
      cm1(I)  = cm(rg_KH(I))            
   END DO   
   
   DO I = 1, n_rg - idx_DS_p + 1
      rg_DS_p(I) = idx_DS_p + I - 1 ! deep stall, positive AoA
      
      aoa2(I) = aoa(rg_DS_p(I))
      ca2(I)  = ca(rg_DS_p(I))
      cn2(I)  = cn(rg_DS_p(I))
      cm2(I)  = cm(rg_DS_p(I))   
   END DO
   
   DO I = 1, idx_DS_n
      rg_DS_n(I) = I ! deep stall, negative AoA
      
      aoa3(I) = aoa(rg_DS_n(I))
      ca3(I)  = ca(rg_DS_n(I))
      cn3(I)  = cn(rg_DS_n(I))
      cm3(I)  = cm(rg_DS_n(I))         
   END DO  
   
   ! f Kirchhoff-Helmholtz spline fit
   IF (.NOT. ALLOCATED( aoa_f ) )      ALLOCATE( aoa_f(n_rg))    
   IF (.NOT. ALLOCATED( f_temp ) )     ALLOCATE( f_temp(n_rg) )   
   
   K=1
   DO I = 1, idx_DS_n -1
      aoa_f(K) = aoa3(I)
      f_temp(K)  = 0.0_DbKi  
      K=K+1
   END DO 
   
   DO I = 1, idx_DS_p - idx_DS_n + 1
      aoa_f(K)  = aoa1(I)
      f_temp(K) = f1(I) 
      K=K+1
   END DO        
   
   DO I = 2, n_rg - idx_DS_p + 1
      aoa_f(K)  = aoa2(I)
      f_temp(K) = 0.0_DbKi  
      K=K+1      
   END DO            

   
   IF (.NOT. ALLOCATED( temp_coefs ) )      ALLOCATE( temp_coefs(n_rg-1,4))   
   
   CALL spline(aoa_f, f_temp, temp_coefs)

   
   DO I=1, n_rg   
      O%FVM_Other%airfoils_LB%pp_f_breaks(NFOILID, I)=aoa_f(I)
   END DO
   DO I=1, n_rg-1      
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_f_coefs(NFOILID, I, J) =  temp_coefs(I,J)
      END DO     
   END DO
   
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 1) = size(aoa_f)   ! Acually size of O%FVM_Other%airfoils_LB%pp_K_Ca_breaks(K) 
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 2) = size(aoa_f) -1   ! Acually size of O%FVM_Other%airfoils_LB%pp_K_Ca_coefs(K)    
      
     
   IF (.NOT. ALLOCATED( f_mdl_aoa ) )      ALLOCATE( f_mdl_aoa(size(aoa1)))   
   CALL ppval(aoa_f, temp_coefs, aoa1, f_mdl_aoa)
   
   DO I=1 , size(aoa1)
      IF (f_mdl_aoa(I) < 0.0_DbKi)  f_mdl_aoa = 0.0_DbKi
   END DO
   
    ! Ca Kirchhoff-Helmholtz error spline fit
   IF (.NOT. ALLOCATED( Ca_mdl1_aoa ) )      ALLOCATE( Ca_mdl1_aoa(size(aoa1)))  
   IF (.NOT. ALLOCATED( Ca_mdl1_err ) )      ALLOCATE( Ca_mdl1_err(size(aoa1))) 
   DO I=1 , size(aoa1)
      Ca_mdl1_aoa(I) = (Cn_alpha*((aoa1(I)-Alpha_0)**2)*sqrt(f_mdl_aoa(I)))
      Ca_mdl1_err(I) = ca1(I) - Ca_mdl1_aoa(I)
   END DO   
   
   IF (.NOT. ALLOCATED( pp_K_KH_Ca ) )      ALLOCATE( pp_K_KH_Ca(size(aoa1)-1, 4)) 
   CALL spline(aoa1, Ca_mdl1_err, pp_K_KH_Ca)
   
   
   ! Ca deep stall error spline fit, positive AoA
   IF (.NOT. ALLOCATED( Ca_mdl2p_aoa ) )      ALLOCATE( Ca_mdl2p_aoa(size(aoa2)))  
   IF (.NOT. ALLOCATED( Ca_mdl2p_err ) )      ALLOCATE( Ca_mdl2p_err(size(aoa2))) 
   
   Ca_mdl2p_aoa = 0.0_DbKi
   DO I=1 , size(aoa2)
      Ca_mdl2p_err(I) = ca2(I) - Ca_mdl2p_aoa(I)
   END DO   
   
   IF (.NOT. ALLOCATED( pp_K_DSp_Ca ) )      ALLOCATE( pp_K_DSp_Ca(size(aoa2)-1, 4)) 
   CALL spline(aoa2, Ca_mdl2p_err, pp_K_DSp_Ca)   
   
   ! Ca deep stall error spline fit, negative AoA
   IF (.NOT. ALLOCATED( Ca_mdl2n_aoa ) )      ALLOCATE( Ca_mdl2n_aoa(size(aoa3)))  
   IF (.NOT. ALLOCATED( Ca_mdl2n_err ) )      ALLOCATE( Ca_mdl2n_err(size(aoa3))) 
   
   Ca_mdl2n_aoa = 0.0_DbKi
   DO I=1 , size(aoa3)
      Ca_mdl2n_err(I) = ca3(I) - Ca_mdl2n_aoa(I)
   END DO   
   
   IF (.NOT. ALLOCATED( pp_K_DSn_Ca ) )      ALLOCATE( pp_K_DSn_Ca(size(aoa3)-1, 4)) 
   CALL spline(aoa3, Ca_mdl2n_err, pp_K_DSn_Ca)      
   
   !  Combine Ca spline fits
   k=1
   DO I= 1, size(aoa3)-1
      O%FVM_Other%airfoils_LB%pp_K_Ca_breaks(NFOILID, K) = aoa3(I)
      K=K+1
   END DO
   DO I= 1, size(aoa1)
      O%FVM_Other%airfoils_LB%pp_K_Ca_breaks(NFOILID, K) = aoa1(I)
      K=K+1
   END DO   
   DO I= 2, size(aoa2)
      O%FVM_Other%airfoils_LB%pp_K_Ca_breaks(NFOILID, K) = aoa2(I)
      K=K+1
   END DO      
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 3) = K-1   ! Acually size (non-zeros) of O%FVM_Other%airfoils_LB%pp_K_Ca_breaks(K) 
   
   K=1
   DO I= 1, size(aoa3)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Ca_coefs(NFOILID, K,J) = pp_K_DSn_Ca(I,J)
      END DO
      K=K+1
   END DO
   DO I= 1, size(aoa1)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Ca_coefs(NFOILID, K,J) = pp_K_KH_Ca(I,J)
      END DO
      K=K+1
   END DO   
   DO I= 1, size(aoa2)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Ca_coefs(NFOILID, K,J) = pp_K_DSp_Ca(I,J)
      END DO
      K=K+1
   END DO        
   
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 4) = K-1   ! Acually size (non-zeros) of O%FVM_Other%airfoils_LB%pp_K_Ca_coefs(K) 

   !.......................................................
   ! Cn deep stall error spline fit, positive AoA
   IF (.NOT. ALLOCATED( Cn_mdl2p_aoa ) )      ALLOCATE( Cn_mdl2p_aoa(size(aoa2)))  
   IF (.NOT. ALLOCATED( Cn_mdl2p_err ) )      ALLOCATE( Cn_mdl2p_err(size(aoa2))) 
   
   DO I=1 , size(aoa2)
      Cn_mdl2p_aoa(I) = 2*sin(aoa2(I))
      Cn_mdl2p_err(I) = Cn2(I) - Cn_mdl2p_aoa(I)
   END DO   
   
   IF (.NOT. ALLOCATED( pp_K_DSp_Cn ) )      ALLOCATE( pp_K_DSp_Cn(size(aoa2)-1, 4)) 
   CALL spline(aoa2, Cn_mdl2p_err, pp_K_DSp_Cn)   
   
   ! Cn deep stall error spline fit, negative AoA
   IF (.NOT. ALLOCATED( Cn_mdl2n_aoa ) )      ALLOCATE( Cn_mdl2n_aoa(size(aoa3)))  
   IF (.NOT. ALLOCATED( Cn_mdl2n_err ) )      ALLOCATE( Cn_mdl2n_err(size(aoa3))) 
   
   DO I=1 , size(aoa3)
      Cn_mdl2n_aoa(I) = 2*sin(aoa3(I)) 
      Cn_mdl2n_err(I) = Cn3(I) - Cn_mdl2n_aoa(I)
   END DO   
   
   IF (.NOT. ALLOCATED( pp_K_DSn_Cn ) )      ALLOCATE( pp_K_DSn_Cn(size(aoa3)-1, 4)) 
   CALL spline(aoa3, Cn_mdl2n_err, pp_K_DSn_Cn)      
   
   !  Combine Cn spline fits
   k=1
   DO I= 1, size(aoa3)
      O%FVM_Other%airfoils_LB%pp_K_Cn_breaks(NFOILID, K) = aoa3(I)
      K=K+1
   END DO
   DO I= 1, size(aoa2)
      O%FVM_Other%airfoils_LB%pp_K_Cn_breaks(NFOILID, K) = aoa2(I)
      K=K+1
   END DO      
   
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 5) = K-1   ! Acually size (non-zeros) of O%FVM_Other%airfoils_LB%pp_K_Cn_breaks(K) 
   
   K=1
   DO I= 1, size(aoa3)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Cn_coefs(NFOILID, K,J) = pp_K_DSn_Cn(I,J)
      END DO
      K=K+1
   END DO   
   
   DO J=1,4
      O%FVM_Other%airfoils_LB%pp_K_Cn_coefs(NFOILID, K,J) = 0.0_DbKi
   END DO
   K=K+1
 
   DO I= 1, size(aoa2)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Cn_coefs(NFOILID, K,J) = pp_K_DSp_Cn(I,J)
      END DO
      K=K+1
   END DO        
   
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 6) = K-1   ! Acually size (non-zeros) of O%FVM_Other%airfoils_LB%pp_K_Cn_coefs
     
   
   !.......................................................
    ! Cm Kirchhoff-Helmholtz error spline fit
   IF (.NOT. ALLOCATED( Cm_mdl1_aoa ) )      ALLOCATE( Cm_mdl1_aoa(size(aoa1)))  
   IF (.NOT. ALLOCATED( Cm_mdl1_err ) )      ALLOCATE( Cm_mdl1_err(size(aoa1))) 
   IF (.NOT. ALLOCATED( xcp_mdl1_aoa ) )     ALLOCATE( xcp_mdl1_aoa(size(aoa1))) 

   DO I=1 , size(aoa1)
      xcp_mdl1_aoa(I) = (5*(1-sqrt(f_mdl_aoa(I)))**2+4*(sqrt(f_mdl_aoa(I))-1))/16
      Cm_mdl1_aoa(I) = cn1(I)*xcp_mdl1_aoa(I) + O%FVM_Other%airfoils_LB%C_m_0(NFOILID)
      Cm_mdl1_err(I) = Cm1(I) - Cm_mdl1_aoa(I)
   END DO   
   
   IF (.NOT. ALLOCATED( pp_K_KH_Cm ) )      ALLOCATE( pp_K_KH_Cm(size(aoa1)-1, 4)) 
   CALL spline(aoa1, Cm_mdl1_err, pp_K_KH_Cm)
   
   
   ! Cm deep stall error spline fit, positive AoA
   IF (.NOT. ALLOCATED( Cm_mdl2p_err ) )      ALLOCATE( Cm_mdl2p_err(size(aoa2))) 
   
   Cm_mdl2p_err(:) = Cm2(:) 
   
   IF (.NOT. ALLOCATED( pp_K_DSp_Cm ) )      ALLOCATE( pp_K_DSp_Cm(size(aoa2)-1, 4)) 
   CALL spline(aoa2, Cm_mdl2p_err, pp_K_DSp_Cm)   
   
   ! Cm deep stall error spline fit, negative AoA
   IF (.NOT. ALLOCATED( Cm_mdl2n_err ) )      ALLOCATE( Cm_mdl2n_err(size(aoa3))) 
   
   Cm_mdl2n_err(:)  = Cm3(:)   
   
   IF (.NOT. ALLOCATED( pp_K_DSn_Cm ) )      ALLOCATE( pp_K_DSn_Cm(size(aoa3)-1, 4)) 
   CALL spline(aoa3, Cm_mdl2n_err, pp_K_DSn_Cm)      
   
   !  Combine Cm spline fits
   k=1
   DO I= 1, size(aoa3)-1
      O%FVM_Other%airfoils_LB%pp_K_Cm_breaks(NFOILID, K) = aoa3(I)
      K=K+1
   END DO
   DO I= 1, size(aoa1)
      O%FVM_Other%airfoils_LB%pp_K_Cm_breaks(NFOILID, K) = aoa1(I)
      K=K+1
   END DO   
   DO I= 2, size(aoa2)
      O%FVM_Other%airfoils_LB%pp_K_Cm_breaks(NFOILID, K) = aoa2(I)
      K=K+1
   END DO      
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 7) = K-1   ! Acually size (non-zeros) of O%FVM_Other%airfoils_LB%pp_K_Cm_breaks(K) 
   
   K=1
   DO I= 1, size(aoa3)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Cm_coefs(NFOILID, K,J) = pp_K_DSn_Cm(I,J)
      END DO
      K=K+1
   END DO
   DO I= 1, size(aoa1)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Cm_coefs(NFOILID, K,J) = pp_K_KH_Cm(I,J)
      END DO
      K=K+1
   END DO   
   DO I= 1, size(aoa2)-1
      DO J=1,4
         O%FVM_Other%airfoils_LB%pp_K_Cm_coefs(NFOILID, K,J) = pp_K_DSp_Cm(I,J)
      END DO
      K=K+1
   END DO        
   
   O%FVM_Other%airfoils_LB%pp_size(NFOILID, 8) = K-1   ! Acually size (non-zeros) of O%FVM_Other%airfoils_LB%pp_K_Cm_coefs(K)    
   
   
   !..........................................
   ! Find first and second break points
   DO J = 1, n_rg 
      IF ( f(J)>0.7 ) THEN
         curve_fit_limit1 = J      
      END IF       
   END DO    
    
   O%FVM_Other%airfoils_LB%alpha_1(NFOILID)  = interp1D_3(f(curve_fit_limit1:curve_fit_limit1+1),  aoa(curve_fit_limit1:curve_fit_limit1+1), 0.7_DbKi)  ! in radius
    
   ! Find alpha_2, the second break point 
   IF (.NOT. ALLOCATED( dca ) )       ALLOCATE( dca(n_rg -1))    
   IF (.NOT. ALLOCATED( aoa_m ) )     ALLOCATE( aoa_m(n_rg -1))    
   DO J = 2, n_rg 
      dca(J-1)   = (ca(j)-ca(J-1))/(aoa(J)-aoa(J-1))
      aoa_m(J-1) = (aoa(J)+aoa(J-1))/2 
   END DO      
    
   IF (.NOT. ALLOCATED( ddca ) )     ALLOCATE( ddca(n_rg -2)) 
   ddca_min = 0.0_DbKi   
   DO J = 2, n_rg -1 
      IF (aoa_m(J) > O%FVM_Other%airfoils_LB%alpha_1(NFOILID) +2*pi/180) THEN ! sliu: 2*pi/180 is very small angle
         ddca(J-1) =  dca(J)-dca(J-1)   
         IF (ddca(j-1) < ddca_min) THEN
            ddca_min = ddca(J-1)
            O%FVM_Other%airfoils_LB%alpha_2(NFOILID) = aoa(J)
            curve_fit_limit2 = J 
         END IF
      END IF 
   END DO  

      
   
   
   !!............debug..................
   !!IF (NFOILID==3) THEN
   !!   DO I = 1, n2 !P%AirFoil%NLIFT(NFOILID)
   !!      WRITE(TEMP1, "(F13.8)") ( f(I) )     
   !!      WRITE(6, "(A)") (TRIM("f: "//TEMP1))  
   !!   END DO         
   !!END IF    
   !!
   !!
   !!IF (NFOILID==3) THEN
   !   !WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%airfoils_LB%Cn_alpha(NFOILID) )     
   !   !WRITE(6, "(A)") (TRIM("Cn_alpha(3): "//TEMP1))  
   !   ! 
   !   !WRITE(TEMP1, "(I)") ( curve_fit_limit1 )     
   !   !WRITE(6, "(A)") (TRIM("curve_fit_limit1: "//TEMP1)) 
   !   !
   !   !WRITE(TEMP1, "(I)") ( curve_fit_limit2 )     
   !   !WRITE(6, "(A)") (TRIM("curve_fit_limit2: "//TEMP1)) 
   !!END IF
   !!!............debug.....................

   
   
   
   IF (ALLOCATED (AOA))      DEALLOCATE (AOA)
   IF (ALLOCATED (cn))       DEALLOCATE (cn)
   IF (ALLOCATED (ca))       DEALLOCATE (ca)
   IF (ALLOCATED (cl))       DEALLOCATE (cl )
   IF (ALLOCATED (cd))       DEALLOCATE (cd)
   IF (ALLOCATED (cm))       DEALLOCATE (cm)
   IF (ALLOCATED (f))        DEALLOCATE (f)
   IF (ALLOCATED (Cn_f))     DEALLOCATE (Cn_f)
   IF (ALLOCATED (rg_KH))    DEALLOCATE (rg_KH)
   IF (ALLOCATED (rg_DS_p))  DEALLOCATE (rg_DS_p)
   IF (ALLOCATED (rg_DS_n))  DEALLOCATE (rg_DS_n)
   
   IF (ALLOCATED (f1))       DEALLOCATE (f1)
   IF (ALLOCATED (aoa1))     DEALLOCATE (aoa1)
   IF (ALLOCATED (ca1))      DEALLOCATE (ca1)
   IF (ALLOCATED (cd1))      DEALLOCATE (cd1)
   IF (ALLOCATED (cn1))      DEALLOCATE (cn1) 
   IF (ALLOCATED (cm1))      DEALLOCATE (cm1)

   IF (ALLOCATED (aoa2))     DEALLOCATE (aoa2)
   IF (ALLOCATED (ca2))      DEALLOCATE (ca2)
   IF (ALLOCATED (cn2))      DEALLOCATE (cn2)
   IF (ALLOCATED (cm2))      DEALLOCATE (cm2)

   IF (ALLOCATED (aoa3))     DEALLOCATE (aoa3)
   IF (ALLOCATED (ca3))      DEALLOCATE (ca3)
   IF (ALLOCATED (cn3))      DEALLOCATE (cn3)
   IF (ALLOCATED (cm3))      DEALLOCATE (cm3)
   
   IF (ALLOCATED (aoa_f))    DEALLOCATE (aoa_f)
   IF (ALLOCATED (f_temp))   DEALLOCATE (f_temp)
   
   IF (ALLOCATED (temp_coefs))      DEALLOCATE (temp_coefs)
   IF (ALLOCATED (f_mdl_aoa))       DEALLOCATE (f_mdl_aoa)
   IF (ALLOCATED (Ca_mdl1_aoa))     DEALLOCATE (Ca_mdl1_aoa)
   IF (ALLOCATED (Ca_mdl1_err))     DEALLOCATE (Ca_mdl1_err)
   IF (ALLOCATED (pp_K_KH_Ca))      DEALLOCATE (pp_K_KH_Ca)
   IF (ALLOCATED (Ca_mdl2p_aoa))    DEALLOCATE (Ca_mdl2p_aoa)
   IF (ALLOCATED (Ca_mdl2p_err))    DEALLOCATE (Ca_mdl2p_err)
   IF (ALLOCATED (pp_K_DSp_Ca))     DEALLOCATE ( pp_K_DSp_Ca)
   IF (ALLOCATED (Ca_mdl2n_aoa))    DEALLOCATE (Ca_mdl2n_aoa)
   IF (ALLOCATED (Ca_mdl2n_err))    DEALLOCATE (Ca_mdl2n_err)
   IF (ALLOCATED (pp_K_DSn_Ca))     DEALLOCATE (pp_K_DSn_Ca)
   IF (ALLOCATED (Cn_mdl2p_aoa))    DEALLOCATE (Cn_mdl2p_aoa)
   IF (ALLOCATED (Cn_mdl2p_err))    DEALLOCATE (Cn_mdl2p_err)
   IF (ALLOCATED (pp_K_DSp_Cn))     DEALLOCATE (pp_K_DSp_Cn)
   IF (ALLOCATED (Cn_mdl2n_aoa))    DEALLOCATE (Cn_mdl2n_aoa)
   IF (ALLOCATED (Cn_mdl2n_err))    DEALLOCATE (Cn_mdl2n_err)
   IF (ALLOCATED (pp_K_DSn_Cn))     DEALLOCATE (pp_K_DSn_Cn)
   IF (ALLOCATED (Cm_mdl1_aoa))     DEALLOCATE (Cm_mdl1_aoa)
   IF (ALLOCATED (Cm_mdl1_err))     DEALLOCATE (Cm_mdl1_err)
   IF (ALLOCATED (xcp_mdl1_aoa))    DEALLOCATE (xcp_mdl1_aoa)
   IF (ALLOCATED (pp_K_KH_Cm))      DEALLOCATE (pp_K_KH_Cm)
   IF (ALLOCATED (Cm_mdl2p_err))    DEALLOCATE (Cm_mdl2p_err)
   IF (ALLOCATED (pp_K_DSp_Cm))     DEALLOCATE (pp_K_DSp_Cm)
   IF (ALLOCATED (Cm_mdl2n_err))    DEALLOCATE (Cm_mdl2n_err)
   IF (ALLOCATED (pp_K_DSn_Cm))     DEALLOCATE (pp_K_DSn_Cm)
   
   IF (ALLOCATED (dca))       DEALLOCATE (dca)
   IF (ALLOCATED (aoa_m))     DEALLOCATE (aoa_m)
   IF (ALLOCATED (ddca))      DEALLOCATE (ddca)
   
   
   

CONTAINS
   !...............................................................................................................................
   SUBROUTINE spline(x, y, s)      
   ! This function determines the piece-wise coefficients for a cubic spline
   ! fit to x and y data. The theory is outlined in
   ! Mathews, J.H., Fink, K.K. (2004). "Numerical Methods Using Matlab, 4th Ed."
   !     Prentice-Hall Inc, Upper Saddle River, NJ. ISBN: 0-13-065248-2.  
   
      IMPLICIT                        NONE
      
         ! Passed variables 
      REAL(DbKi),DIMENSION(:),             INTENT(IN   )  :: x          
      REAL(DbKi),DIMENSION(:),             INTENT(IN   )  :: y          
      REAL(DbKi),DIMENSION(size(x)-1,4),   INTENT(  OUT)  :: s  !coefs    
   
      INTEGER(IntKi)   :: i, j
      INTEGER(IntKi)   :: n ! = size(x)-1      
      REAL(DbKi), DIMENSION(size(x)-1)           :: h, d, a, b, c
      REAL(DbKi), DIMENSION(size(x)-2)           :: V, m_interior
      REAL(DbKi), DIMENSION(size(x)-2,size(x)-2) :: HH, inv_HH
      REAL(DbKi), DIMENSION(size(x))             :: m
      
      n = size(x)-1   
      
      h = 0.0_DbKi   
      d = 0.0_DbKi  
      a = 0.0_DbKi  
      b = 0.0_DbKi  
      c = 0.0_DbKi  
      V = 0.0_DbKi  
      HH = 0.0_DbKi  
      s = 0.0_DbKi  
      
       ! First derivatives (slope between points)
      do i=1,n
        h(i) = x(i+1) - x(i)
        d(i) = (y(i+1) - y(i)) / h(i)
      end do   
      
      
       !Solve for second derivative coefficients
      do i=2,n
          a(i) = h(i-1)
          b(i) = 2*(h(i-1) + h(i))
          c(i) = h(i)

          V(i-1) = 6*(d(i) - d(i-1))
         
          do j = 2,n
              if (j == i-1) THEN
                  HH(i-1,j-1) = a(i)
              else if (j == i) THEN
                  HH(i-1,j-1) = b(i)
              else if (j == i+1) THEN
                  HH(i-1,j-1) = c(i)
              end if
          end do   
      end do         
      
       ! Natural spline end point constraints
      HH(1,1) = 2*(h(1)+h(2))
      HH(1,2) = h(2)
      HH(n-1,n-1) = 2*(h(n-1)+h(n))
      HH(n-1,n-2) = h(n-1)     
      
       ! Solve for interior points of m by matrix inversion (in Matlab: m_interior = H\V)
      call matrix_inverse(n-1, HH, inv_HH )

      
      do i = 1, n-1
         m_interior(i) = 0.0        
         do j = 1, n-1
            m_interior(i) = m_interior(i) + inv_HH(i, j)*V(j) 
         end do   
      end do   
           
      
       ! Natural spline end point constraints
      m = 0.0_DbKi  
      do i= 2,n
         m(i) = m_interior(i-1)
      end do   
      
       ! Find piecewise coefficients
       ! Solution form:
       !   f_k(w) = s(k,1)*w^3 + s(k,2)*w^2 + s(k,3)*w + s(k,4);
      do i = 1,n
          s(i,1) = (m(i+1)-m(i))/(6*h(i))
          s(i,2) = m(i)/2
          s(i,3) = d(i)-(h(i)*(2*m(i)+m(i+1)))/6
          s(i,4) = y(i)
      end do         
            
   
    END SUBROUTINE spline
   !...............................................................................................................................
   SUBROUTINE matrix_inverse(n, a, Inv_a)
   ! Computing Inverse matrix. Based on the Doolittle LU method
   ! Reference: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
   !            http://en.wikipedia.org/wiki/LU_decomposition#Doolittle_algorithm
   !....................................................................

      IMPLICIT NONE 

         ! Passed variables
      INTEGER(IntKi),               INTENT(IN   )  :: n           ! The matrix a has size (N,N)
      REAL(DbKi), DIMENSION(n,n),   INTENT(INOUT)  :: a           ! Square matrix
      REAL(DbKi), DIMENSION(n,n),   INTENT(  OUT)  :: Inv_a       ! Inverse of a      

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


   END SUBROUTINE matrix_inverse   
   !...............................................................................................................................
   SUBROUTINE ppval(x, coefs, xx, yy)   
   ! Evaluate piecewise polynomial
   ! Refer to Matlab function ppval.m
   
      IMPLICIT                        NONE
      
         ! Passed variables 
      REAL(DbKi),DIMENSION(:),         INTENT(IN   )  :: x      
      REAL(DbKi),DIMENSION(:,:),       INTENT(IN   )  :: coefs    
      REAL(DbKi),DIMENSION(:),         INTENT(IN   )  :: xx   
      REAL(DbKi),DIMENSION(size(xx)),  INTENT(  OUT)  :: yy        
  
      
   
      INTEGER(IntKi)   :: n1, n2
      INTEGER(IntKi)   :: I, J, k, M

      n1 = size(x)
      n2 = size(xx)
      
      
      DO I=1, n2
         DO J =1,n1-1
            IF (xx(I) > x(J)  .AND. xx(I) <= x(J+1)) THEN
               K = J
               EXIT
            END IF
         END DO
         
         yy(I)=coefs(K,1)
         DO M=2,4
            yy(I)=(xx(I) - x(K)) *yy(I)+ coefs(K,M)
         END DO              

      END DO
          
   
    END SUBROUTINE ppval
   !======================================================================



END SUBROUTINE LB_AfData_spline_fit 
!==================================================================================================================================  
SUBROUTINE LB_DynStall(p, O, xd, Iiter, ErrStat, ErrMess, IElement, IBlade, C_l, C_d, C_m, I)
!  
!........................................................................

   IMPLICIT                        NONE

      ! Passed variables

   TYPE(AD_ParameterType),      INTENT(IN   )   :: p              ! The module's parameter data
   TYPE(AD_OtherStateType),     INTENT(INOUT)   :: O              ! Other/optimization states   
   TYPE(AD_DiscreteStateType),  INTENT(IN   )   :: xd             ! Discrete states at t   
   INTEGER(IntKi),              INTENT(  OUT)   :: ErrStat        ! The error status code
   CHARACTER(*),                INTENT(  OUT)   :: ErrMess        ! The error message, if an error occurred
   INTEGER(IntKi),INTENT(IN)  :: IElement
   INTEGER(IntKi),INTENT(IN)  :: IBlade
   REAL(DbKi),INTENT(OUT)     :: C_d
   REAL(DbKi),INTENT(OUT)     :: C_l
   REAL(DbKi),INTENT(OUT)     :: C_m
   INTEGER(IntKi),INTENT(IN)  :: I         ! NFOIL(J)
   INTEGER(IntKi), INTENT(IN) :: Iiter     ! The counter of iteration
   

      ! Local variables
   REAL(DbKi)   :: aoa
   REAL(DbKi)   :: Vtot
   REAL(DbKi)   :: dalpha
   REAL(DbKi)   :: A1, A2, b1, b2
   REAL(DbKi)   :: C_n_alpha
   REAL(DbKi)   :: DT          ! timestep duration of WInDS
   REAL(DbKi)   :: DS
   REAL(DbKi)   :: a           ! Sound speed
   REAL(DbKi)   :: M           ! Mach number
   REAL(DbKi)   :: Beta 
   REAL(DbKi)   :: T_I           
   REAL(DbKi)   :: alpha_e      ! Eq(3.10)
   REAL(DbKi)   :: C_n_c        ! Eq(3.9)
   REAL(DbKi)   :: k_alpha      ! Eq(3.15)
   REAL(DbKi)   :: C_n_nc_alpha ! Eq(3.18)
   REAL(DbKi)   :: k_q          ! Eq(3.16)
   REAL(DbKi)   :: C_n_nc_q     ! Eq(3.17)
   REAL(DbKi)   :: C_n_nc
   REAL(DbKi)   :: C_n_prime ! Eq(3.28)
   REAL(DbKi)   :: alpha_f   ! Eq(3.30) 
   REAL(DbKi)   :: f
   REAL(DbKi)   :: E_f       ! Eq(3.45)
   REAL(DbKi)   :: E_l_f     ! Eq(3.46)
   REAL(DbKi)   :: C_f_n 
   REAL(DbKi)   :: K_n       ! Eq(3.38)
   REAL(DbKi)   :: E_v       ! Eq(3.48)
   REAL(DbKi)   :: E_v_l     ! Eq(3.49)
   REAL(DbKi)   :: C_n       ! Eq(3.41)
   REAL(DbKi)   :: D, E      ! Eq(3.36)
   REAL(DbKi)   :: phi       ! Eq(3.36)
   REAL(DbKi)   :: C_a       ! Eq(3.35)
   REAL(DbKi)   :: C_d_p     ! last two parts of Eq(3.43)
   INTEGER(IntKi) :: check
   INTEGER(IntKi) :: n1,n2 ! actual size
   REAL(DbKi)   :: K_Cn, K_Ca, K_Cm
   REAL(DbKi)   :: C_m_f, X_cp_f, C_m_nc_q, C_m_v, X_cp_v
   
   
   CHARACTER(LEN = 1024)   :: TEMP1, TEMP2, TEMP3  ! debug
   
   
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""
   
   
   
   C_l = 0.0
   C_d = 0.0
   
   aoa   = O%FVM_Other%KJ%AOA(1, 1, 1, IElement, IBlade)
   Vtot  = O%FVM_Other%KJ%Vtot(1, 1, 1, IElement, IBlade)
   O%FVM_Other%DS%global_current_alpha(IElement, IBlade) = aoa
   dalpha = aoa - O%FVM_Other%DS%global_previous_alpha(IElement, IBlade)
   O%FVM_Other%DS%attached_current_q(IElement, IBlade) =  dalpha * p%Blade%C(IElement) / Vtot
   A1 = O%FVM_Other%airfoils_LB%A1(I)
   A2 = O%FVM_Other%airfoils_LB%A2(I)
   b1 = O%FVM_Other%airfoils_LB%b1(I)   
   b2 = O%FVM_Other%airfoils_LB%b2(I)
   C_n_alpha = O%FVM_Other%airfoils_LB%Cn_alpha(I)
   
   DT = p%FVM%DT_WINDS
   DS = (2/p%Blade%C(IElement)) * Vtot * DT

   !....................................................................................
   ! Unsteady Attached Flow
   !...........................................................   
   a = 343.2    ! m/s speed of sound at 20C
   M = Vtot / a
   Beta = sqrt(1-M**2)
   T_I  = p%Blade%C(IElement) / a   
   
   !IF ( M**2 > 1.0)  THEN
   !   ErrStat = ErrID_Fatal 
   !   ErrMess = ' Error (in WInDS): Vtot is too large in LB_DynStall subroutine.'
   !   RETURN   
   !END IF          
   !.............................................
   ! Dynamic Stall Model
    ! Effective angle of Attack
   O%FVM_Other%DS%attached_current_X1(IElement, IBlade) = O%FVM_Other%DS%attached_previous_X1(IElement, IBlade) * exp(-b1*DS)   &
                                                             +(A1/6)*dalpha*(1+4*exp(-b1*(DS/2))+exp(-b1*DS))   ! Eq(3.11)
   O%FVM_Other%DS%attached_current_X2(IElement, IBlade) = O%FVM_Other%DS%attached_previous_X2(IElement, IBlade) * exp(-b2*ds)   &
                                                             +(A2/6)*dalpha*(1+4*exp(-b2*(DS/2))+exp(-b2*DS))   ! Eq(3.12)
   alpha_e = aoa - O%FVM_Other%DS%attached_current_X1(IElement, IBlade) - O%FVM_Other%DS%attached_current_X2(IElement, IBlade)  ! Eq(3.10)
   
    ! Circulatory lift for changes in angle of attack
   C_n_c = C_n_alpha * alpha_e         ! Eq(3.9)
   
    ! Noncirculatory normal force for changes in angle of attack
   k_alpha = 0.75 / ((1-M)+pi*(Beta**2)*(M**2)*(A1*b1+A2*b2))    ! Eq(3.15)
   O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) = dalpha / DT
   O%FVM_Other%DS%attached_current_dK_alpha(IElement, IBlade) = &
                O%FVM_Other%DS%attached_previous_dK_alpha(IElement, IBlade) * exp(-DT/(2*k_alpha*T_I)) +  & 
                (O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) - O%FVM_Other%DS%attached_previous_K_alpha(IElement, IBlade))  &
                 * exp(-DT/(2*k_alpha*T_I)) ! Eq(3.21)
   C_n_nc_alpha  = ((4*k_alpha*T_I)/M) * ( O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) -  &
                    O%FVM_Other%DS%attached_current_dK_alpha(IElement, IBlade) ) ! Eq(3.18)
   
    ! Noncirculatory normal force for changes in pitch rate 
   k_q = 0.75 / ((1-M)+2*pi*(Beta**2)*(M**2)*(A1*b1+A2*b2))     ! Eq(3.16)
   O%FVM_Other%DS%attached_current_K_q(IElement, IBlade) =    &
                ( O%FVM_Other%DS%attached_current_q(IElement, IBlade) - O%FVM_Other%DS%attached_previous_q(IElement, IBlade) )/ DT  ! Eq(3.20)
   O%FVM_Other%DS%attached_current_dK_q(IElement, IBlade) = &
                O%FVM_Other%DS%attached_previous_dK_q(IElement, IBlade) * exp(-DT/(2*k_q*T_I)) +  & 
                (O%FVM_Other%DS%attached_current_K_q(IElement, IBlade) - O%FVM_Other%DS%attached_previous_K_q(IElement, IBlade)) * exp(-DT/(2*k_q*T_I))  ! Eq(3.22)
   C_n_nc_q  = ((k_q*T_I)/M) * ( O%FVM_Other%DS%attached_current_K_q(IElement, IBlade) - O%FVM_Other%DS%attached_current_dK_q(IElement, IBlade) )        ! Eq(3.17)

     ! Total normal force under attached flow
   C_n_nc = C_n_nc_alpha + C_n_nc_q
   O%FVM_Other%DS%attached_current_C_p_n(IElement, IBlade) = C_n_c + C_n_nc     ! Eq(3.23)
  
   
   !....................................................................................
   ! Trailing Edge Flow Separation
   !...........................................................
     ! Ersatz lift parameter for leading edge seperation onset
   O%FVM_Other%DS%TEsep_current_D_p(IElement, IBlade) = O%FVM_Other%DS%TEsep_previous_D_p(IElement, IBlade) * exp(-DS/O%FVM_Other%airfoils_LB%Tp(I)) +  & 
                         ( O%FVM_Other%DS%attached_current_C_p_n(IElement, IBlade) - O%FVM_Other%DS%attached_previous_C_p_n(IElement, IBlade))    &
                          * exp(-DS/( 2*O%FVM_Other%airfoils_LB%Tp(I)))  ! Eq(3.29)
   C_n_prime = O%FVM_Other%DS%attached_current_C_p_n(IElement, IBlade) - O%FVM_Other%DS%TEsep_current_D_p(IElement, IBlade)    ! Eq(3.28)
   
      ! Leading Edge Vortex Position
   IF (C_n_prime > O%FVM_Other%airfoils_LB%C_n_1(I) .AND. dalpha > 0) THEN
       O%FVM_Other%DS%LEsep_vortex_shed_start(IElement, IBlade) = 1_IntKi
   END IF
   
   IF (O%FVM_Other%DS%LEsep_vortex_shed_start(IElement, IBlade) == 1_IntKi) THEN
       O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) = O%FVM_Other%DS%LEsep_previous_tau_v(IElement, IBlade) + DS * 0.45_DbKi
   ELSE 
       O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) = 0.0_DbKi
   END IF
   
   IF (O%FVM_Other%DS%LEsep_vortex_shed_start(IElement, IBlade) == 1_IntKi) THEN
       IF ( C_n_prime < O%FVM_Other%airfoils_LB%C_n_1(I) ) THEN
           IF ( dalpha<0 ) THEN
               IF (O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) > O%FVM_Other%airfoils_LB%Tvl(I))  THEN
                   O%FVM_Other%DS%LEsep_vortex_shed_start(IElement, IBlade) = 0_IntKi
               END IF
           END IF
       END IF
   END IF
   
    !Effective trailing edge seperation point
   alpha_f = (C_n_prime) / C_n_alpha         ! Eq(3.30)    
      
   
   
   if (isnan(alpha_f)) THEN
      print *, "alpha_f : ", Vtot, Beta, O%FVM_Other%DS%attached_current_C_p_n(IElement, IBlade), O%FVM_Other%DS%TEsep_current_D_p(IElement, IBlade) 
   end if      
   
   
    ! Unsteady trailing edge flow separation using spline fit, f'  
   n1=O%FVM_Other%airfoils_LB%pp_size(I,1)
   n2=O%FVM_Other%airfoils_LB%pp_size(I,2)
   CALL ppval(O%FVM_Other%airfoils_LB%pp_f_breaks(I,1:n1), O%FVM_Other%airfoils_LB%pp_f_coefs(I,1:n2,1:4), n1, alpha_f, O%FVM_Other%DS%TEsep_current_f_prime(IElement, IBlade)) 
 
   
   IF (O%FVM_Other%DS%TEsep_current_f_prime(IElement, IBlade) < 0.0_DbKi ) O%FVM_Other%DS%TEsep_current_f_prime(IElement, IBlade) = 0.0_DbKi
   

   ! Time lagged unsteady trailing edge flow seperation point 
   E_f = EXP(-DS/O%FVM_Other%airfoils_LB%Tf(I))  ! Eq(3.45)
   E_l_f = E_f ** O%FVM_Other%DS%global_previous_sigma1(IElement, IBlade)  ! Eq(3.46)
   O%FVM_Other%DS%TEsep_current_D_f(IElement, IBlade) = O%FVM_Other%DS%TEsep_previous_D_f(IElement, IBlade)*E_f + &
                                                        (O%FVM_Other%DS%TEsep_current_f_prime(IElement, IBlade) - &
                                                        O%FVM_Other%DS%TEsep_previous_f_prime(IElement, IBlade))*SQRT(E_l_f) ! Eq(3.44)
   O%FVM_Other%DS%TEsep_current_f_2prime(IElement, IBlade) = O%FVM_Other%DS%TEsep_current_f_prime(IElement, IBlade) - &
                                                             O%FVM_Other%DS%TEsep_current_D_f(IElement, IBlade) ! Eq(3.31)
   
    ! Normal force with unsteady trailing edge flow seperation 
   n1=O%FVM_Other%airfoils_LB%pp_size(I,5)
   n2=O%FVM_Other%airfoils_LB%pp_size(I,6)
   IF (alpha_f <= O%FVM_Other%airfoils_LB%aoa_KHn(I) .OR. alpha_f >= O%FVM_Other%airfoils_LB%aoa_KHp(I)) THEN
      CALL ppval(O%FVM_Other%airfoils_LB%pp_K_Cn_breaks(I,1:n1), O%FVM_Other%airfoils_LB%pp_K_Cn_coefs(I,1:n2,1:4), n1, alpha_f, K_Cn ) 
      C_f_n = 2*sin(alpha_f) + K_Cn
   ELSE
      C_f_n = C_n_alpha * (alpha_f - O%FVM_Other%airfoils_LB%Alpha_0(I)) * ((1+sqrt(O%FVM_Other%DS%TEsep_current_f_2prime(IElement, IBlade) ))/2)**2
   END IF
   
   
   !....................................................................................
   ! Leading Edge Flow Separation
   !...........................................................   
    ! Normal force from vortex shedding
   K_n = ((1+SQRT(O%FVM_Other%DS%TEsep_current_f_2prime(IElement, IBlade)))**2) / 4   ! Eq(3.38)
   O%FVM_Other%DS%LEsep_current_C_v(IElement, IBlade) =  C_n_c*(1-K_n)      ! Eq(3.37)
   E_v = EXP(-DS/O%FVM_Other%airfoils_LB%Tv(I))                  ! Eq(3.48)
   E_v_l = E_v ** O%FVM_Other%DS%global_previous_sigma3(IElement, IBlade)   ! Eq(3.49)
   O%FVM_Other%DS%LEsep_current_C_v_n(IElement, IBlade) = O%FVM_Other%DS%LEsep_previous_C_v_n(IElement, IBlade)*E_v_l +     &
                         (O%FVM_Other%DS%LEsep_current_C_v(IElement, IBlade) - O%FVM_Other%DS%LEsep_previous_C_v(IElement, IBlade)) *SQRT(E_v_l)  ! Eq(3.39)
   IF (O%FVM_Other%DS%LEsep_current_C_v_n(IElement, IBlade) < 0 .OR.  (O%FVM_Other%DS%LEsep_vortex_shed_start(IElement, IBlade) ==0  .AND.  dalpha<0)) THEN
      O%FVM_Other%DS%LEsep_current_C_v_n(IElement, IBlade) = 0.0_DbKi
   END IF
   
    ! Total Unsteady Normal Force   
   C_n = C_n_nc + C_f_n + O%FVM_Other%DS%LEsep_current_C_v_n(IElement, IBlade) ! old: + O%FVM_Other%airfoils_LB%C_n_0(I) ! Eq(3.41)
  
   
    ! Unsteady Axial Force Ca and Pitching Moment Cm
   n1=O%FVM_Other%airfoils_LB%pp_size(I,3)
   n2=O%FVM_Other%airfoils_LB%pp_size(I,4)   
   CALL ppval(O%FVM_Other%airfoils_LB%pp_K_Ca_breaks(I,1:n1), O%FVM_Other%airfoils_LB%pp_K_Ca_coefs(I,1:n2,1:4), n1, aoa, K_Ca ) 
   
   
   n1=O%FVM_Other%airfoils_LB%pp_size(I,7)
   n2=O%FVM_Other%airfoils_LB%pp_size(I,8)      
   CALL ppval(O%FVM_Other%airfoils_LB%pp_K_Cm_breaks(I,1:n1), O%FVM_Other%airfoils_LB%pp_K_Cm_coefs(I,1:n2,1:4), n1, aoa, K_Cm )       
   
   IF (aoa < O%FVM_Other%airfoils_LB%aoa_KHn(I) .OR. aoa > O%FVM_Other%airfoils_LB%aoa_KHp(I)) THEN
      C_a = K_Ca
      C_m_f = K_Cm
   ELSE
      C_a = K_Ca + C_n_alpha*(aoa-O%FVM_Other%airfoils_LB%Alpha_0(I))**2 * sqrt(O%FVM_Other%DS%TEsep_current_f_2prime(IElement, IBlade))
      
      x_cp_f = (5*(1-sqrt(O%FVM_Other%DS%TEsep_current_f_2prime(IElement, IBlade)))**2+4*(sqrt(O%FVM_Other%DS%TEsep_current_f_2prime(IElement, IBlade))-1))/16
      C_m_f = x_cp_f*C_n_c + O%FVM_Other%airfoils_LB%C_m_0(I) + K_Cm  
   END IF

   ! Noncirculatory pitching moment coef for changes in pitch rate
   O%FVM_Other%DS%attached_current_ddK_q(IElement, IBlade) = &
                            O%FVM_Other%DS%attached_previous_ddK_q(IElement, IBlade) *  exp(-DT/(k_q**2 * T_I)) + &
                            ( O%FVM_Other%DS%attached_current_K_q(IElement, IBlade) - &
                             O%FVM_Other%DS%attached_previous_K_q(IElement, IBlade)  ) * SQRT(exp(-DT/(k_q**2 * T_I)))  
   C_m_nc_q= C_n_nc_q/4 -((k_q**2*T_I)/(3*M))*(O%FVM_Other%DS%attached_current_K_q(IElement, IBlade) -  &
             O%FVM_Other%DS%attached_current_ddK_q(IElement, IBlade))
   
   
   ! Pitching moment coef contribution from vortex shedding
   x_cp_v =0.2*(1-cos(pi*O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) / O%FVM_Other%airfoils_LB%Tvl(I)))
   C_m_v = -x_cp_v * O%FVM_Other%DS%LEsep_current_C_v_n(IElement, IBlade)
   
   
    ! Total unsteady pitching moment coef
   C_m = C_m_nc_q + C_m_f + C_m_v      
   
    ! Unsteady Lift and Drag from Cn and Ca
   C_d = C_n*sin(aoa) - C_a*cos(aoa) 
   C_l = C_n*cos(aoa) + C_a*sin(aoa)    ! Eq(3.42)
   

   !! ................ debug
   !IF (Iiter == 1 .AND. IElement == 4 .AND. IBlade==1 .AND. O%WINDS_Timestep == 2) THEN
   !   WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%KJ%VEL_TOT(1, 2, 1, 4, 1) ) 
   !   WRITE(TEMP2, "(F13.8)") ( O%FVM_Other%KJ%VEL_TOT(2, 2, 1, 4, 1) )
   !   WRITE(TEMP3, "(F13.8)") ( O%FVM_Other%KJ%VEL_TOT(3, 2, 1, 4, 1) )      
   !   WRITE(6, "(A)") (TRIM("...%KJ%VEL_TOT: "//TEMP1) ) 
   !   WRITE(6, "(A)") (TRIM(TEMP2)) 
   !   WRITE(6, "(A)") (TRIM(TEMP3)) 
   !   
   !END IF
   !!...................debug

   
   ! ................ debug
   !IF (Iiter == 20 .AND. IElement == 4 .AND. IBlade==1 .AND. O%WINDS_Timestep == 2) THEN
   !    
   !    
   !   WRITE(TEMP1, "(F13.8)") ( aoa ) 
   !   WRITE(6, "(A)") (TRIM("aoa: "//TEMP1))  
   !   !WRITE(TEMP1, "(F13.8)") ( O%FVM_Other%DS%global_previous_alpha(IElement, IBlade) ) 
   !   !WRITE(6, "(A)") (TRIM("...%global_previous_alpha: "//TEMP1))       
   !   !
   !  
   !   WRITE(TEMP1, "(F13.8)") ( C_n ) 
   !   WRITE(6, "(A)") (TRIM("C_n: "//TEMP1))       
   !   
   !    WRITE(TEMP1, "(F13.8)") ( C_a ) 
   !   WRITE(6, "(A)") (TRIM("C_a: "//TEMP1))       
   !   
   !   WRITE(TEMP1, "(F13.8)") ( C_l ) 
   !   WRITE(6, "(A)") (TRIM("C_l: "//TEMP1)) 
   !   
   !   WRITE(6, "(A)") ("........................... ")  
   !  
   !END IF
   !...................debug
      
     
   
   
    ! Time Constant Modification, sigma1
   IF (O%FVM_Other%DS%TEsep_current_f_2prime(IElement, IBlade) < O%FVM_Other%DS%TEsep_previous_f_2prime(IElement, IBlade)) THEN
       ! Flow is seperating   Table 3.1
       O%FVM_Other%DS%LEsep_seperating(IElement, IBlade) = 1
       IF (C_n_prime <= O%FVM_Other%airfoils_LB%C_n_1(I)) THEN
          O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 1.0_DbKi
       ELSE
          O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 1.75_DbKi
       END IF
       IF (O%FVM_Other%DS%TEsep_previous_f_2prime(IElement, IBlade) <= 0.7) THEN
          O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 2.0_DbKi
       END IF
       IF (O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) < 0) THEN
          O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 2.0_DbKi
       END IF
   ELSE        
      ! Flow is reattaching    Table 3.2
       O%FVM_Other%DS%LEsep_seperating(IElement, IBlade) = 0
       O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 0.5_DbKi
       
       IF (C_n_prime <= O%FVM_Other%airfoils_LB%C_n_1(I)) THEN
          O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 0.5_DbKi
       END IF
       
       IF (O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) >= 0 ) THEN
          IF (O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) <= O%FVM_Other%airfoils_LB%Tvl(I))  THEN
             O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 0.25_DbKi
          END IF
       END IF
       
       IF (O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) > 0) THEN
          O%FVM_Other%DS%global_current_sigma1(IElement, IBlade) = 0.75
       END IF
       
   END IF
   
    ! Time Constant Modification, sigma3   Table 3.3
   check = 0_IntKi
   O%FVM_Other%DS%global_current_sigma3(IElement, IBlade)  = 1
   IF ( O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) >= O%FVM_Other%airfoils_LB%Tvl(I)) THEN
      IF (O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade)<= O%FVM_Other%airfoils_LB%Tvl(I)*2 ) THEN
         O%FVM_Other%DS%global_current_sigma3(IElement, IBlade)  = 3.0 
         check = 1_IntKi
      END IF
   END IF
   
   IF (O%FVM_Other%DS%LEsep_seperating(IElement, IBlade) == 0_IntKi) THEN
      O%FVM_Other%DS%global_current_sigma3(IElement, IBlade)  = 4.0
      check = 1_IntKi
   END IF

   IF ( O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) >= 0) THEN
      IF (O%FVM_Other%DS%LEsep_current_tau_v(IElement, IBlade) <= O%FVM_Other%airfoils_LB%Tvl(I) ) THEN
         O%FVM_Other%DS%global_current_sigma3(IElement, IBlade)  = 1
         IF (O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) < 0 ) THEN
            O%FVM_Other%DS%global_current_sigma3(IElement, IBlade)  = 2.0
         END IF
         check = 1_IntKi
      END IF
   END IF
   
   IF (check == 0) THEN
      IF (O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) < 0 ) THEN
         O%FVM_Other%DS%global_current_sigma3(IElement, IBlade)  = 4.0
      END IF      
   END IF
   
   IF (O%FVM_Other%DS%LEsep_seperating(IElement, IBlade) == 0) THEN
       IF (O%FVM_Other%DS%attached_current_K_alpha(IElement, IBlade) < 0) THEN
          O%FVM_Other%DS%global_current_sigma3(IElement, IBlade)  = 1.0 
       END IF
   END IF

   
   
CONTAINS
     !...............................................................................................................................
   SUBROUTINE ppval(x, coefs, n, xx, yy)   
   ! Evaluate piecewise polynomial
   ! Refer to Matlab function ppval.m
   
      IMPLICIT                        NONE
      
         ! Passed variables 
      REAL(DbKi),DIMENSION(:),         INTENT(IN   )  :: x      
      REAL(DbKi),DIMENSION(:,:),       INTENT(IN   )  :: coefs 
      INTEGER(IntKi),                  INTENT(IN   )  :: n
      REAL(DbKi),                      INTENT(IN   )  :: xx   
      REAL(DbKi),                      INTENT(  OUT)  :: yy        
  
   
      INTEGER(IntKi)   :: J, k, M

      DO J =1,n-1
         IF (xx >= x(J)  .AND. xx <= x(J+1)) THEN
            K = J
            EXIT
         END IF
      END DO
      
      yy=coefs(K,1)
      DO M=2,4
         yy=(xx - x(k)) * yy + coefs(K,M)
      END DO             
   
    END SUBROUTINE ppval
   !======================================================================
   
   
   
END SUBROUTINE LB_DynStall     
!==================================================================================================================================  
FUNCTION interp1D_1(X,Y,vx)
! ID interpolation to find y when x=vx. Y=f(X).
!
! ! X is DbKi and Y is ReKi
!...................................................

   REAL(DbKi)               :: interp1D_1
   REAL(DbKi),Dimension(1,2), INTENT(IN)   :: X
   REAL(ReKi),Dimension(1,2), INTENT(IN)   :: Y   
   REAL(DbKi)               :: vx

   interp1D_1 = Y(1,1) + (Y(1,2) - Y(1,1)) * (vx - X(1,1)) / (X(1,2) - X(1,1))  
     
END FUNCTION  interp1D_1
!==================================================================================================================================  
FUNCTION interp1D_2(X,Y,vx)
! ID interpolation to find y when x=vx. Y=f(X).
!
! X is ReKi and Y is DbKi
!...................................................

   REAL(DbKi)               :: interp1D_2
   REAL(ReKi),Dimension(1,2), INTENT(IN)   :: X
   REAL(DbKi),Dimension(1,2), INTENT(IN)   :: Y   
   REAL(DbKi)               :: vx

   interp1D_2 = Y(1,1) + (Y(1,2) - Y(1,1)) * (vx - X(1,1)) / (X(1,2) - X(1,1))  
     
END FUNCTION  interp1D_2
!==================================================================================================================================  
FUNCTION interp1D_3(X,Y,vx)
! ID interpolation to find y when x=vx. Y=f(X).
! 
! X and Y are DbKi
!...................................................

   REAL(DbKi)               :: interp1D_3
   REAL(DbKi),Dimension(2), INTENT(IN)   :: X, Y
   REAL(DbKi)               :: vx

   interp1D_3 = Y(1) + (Y(2) - Y(1)) * (vx - X(1)) / (X(2) - X(1))  
     
END FUNCTION  interp1D_3
!==================================================================================================================================    
    
END MODULE WINDS_DS    