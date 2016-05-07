   ! NOTE: This source file contains dummy placeholders for ALL of the
   !       user-specified routines available in FAST.  These routines
   !       are as follows:
   !          Routine       Description
   !          ------------  ---------------------------------------------------
   !          PitchCntrl()  User-specified blade pitch control (either
   !                        independent or rotor-collective) model.
   !          UserGen()     User-specified generator torque and power model.
   !          UserHSSBr()   User-specified high-speed shaft brake model.
   !          UserPtfmLd()  User-specified platform loading model.
   !          UserRFrl()    User-specified rotor-furl spring/damper model.
   !          UserTeet()    User-specified rotor-teeter spring/damper model.
   !          UserTFin()    User-specified tail fin aerodynamics model.
   !          UserTFrl()    User-specified tail-furl spring/damper model.
   !          UserVSCont()  User-specified variable-speed torque and power
   !                        control model.
   !          UserYawCont() User-specified nacelle-yaw control model.
   !       In order to interface FAST with your own user-specified routines,
   !       you can develop your own logic within these dummy placeholders and
   !       recompile FAST; OR comment out the appropriate dummy placeholders,
   !       create your own routines in their own source files, and recompile
   !       FAST while linking in these additional source files.  For example,
   !       the executable version of FAST that is distributed with the FAST
   !       archive is linked with the example PitchCntrl() routine contained in
   !       source file PitchCntrl_ACH.f90 and the example UserGen() and
   !       UserVSCont() routines contained in source file UserVSCont_KP.f90;
   !       thus, the dummy placeholders for routines PitchCntrl(), UserGen(),
   !       and UserVSCont() are commented out within this source file.  The
   !       example pitch controller was written by Craig Hansen (ACH) and the
   !       example generator and variable speed controllers were written by
   !       Kirk Pierce (KP).  Please see the aforementioned source files for
   !       additional information on these example user-specified routines.

   ! NOTE: If you (the user) wants to access the current value of ANY of the
   !       output parameters available as outputs from FAST from your
   !       user-defined routines, then do the following:
   !          (1) USE MODULE Output() in your routine.
   !          (2) Access the output parameter by typing "AllOuts(OutName)",
   !              where OutName is the PRIMARY name of the output parameter.
   !              For example, to access the current value of the in-plane
   !              bending moment at the root of blade 1 (in kN�m), type in
   !              "AllOuts(RootMxc1)", since RootMxc1 is the primary name of
   !              this output parameter--RootMIP1 will not work in place of
   !              RootMxc1, since it is a SECONDARY name.  Also, you CANNOT use
   !              the prefixes ("-", "_", "m", or "M") in front of OutName to
   !              reverse the sign of the selected output channel.
   !       Note that OutName DOES NOT have to be one of the output parameters
   !       you listed in OutList from your primary input file.  Also note that
   !       this technique WILL also work for user-defined routines written for
   !       ADAMS datasets extracted using the FAST-to-ADAMS preprocessor.

!=======================================================================
!SUBROUTINE PitchCntrl ( BlPitch, ElecPwr, LSS_Spd, TwrAccel, NumBl, ZTime, DT, DirRoot, BlPitchCom )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! blade pitch control model (either independent or rotor-collective).
!   ! Modify this code to create your own model.
!
!
!USE                             Precision
!
!
!IMPLICIT                        NONE
!
!
!   ! Passed variables:
!
!INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
!
!REAL(ReKi), INTENT(IN )      :: BlPitch   (NumBl)                               ! Current values of the blade pitch angles, rad.
!REAL(DbKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
!REAL(ReKi), INTENT(IN )      :: ElecPwr                                         ! Electrical power, watts.
!REAL(ReKi), INTENT(IN )      :: LSS_Spd                                         ! LSS speed (rad/s)
!REAL(ReKi), INTENT(OUT)      :: BlPitchCom(NumBl)                               ! Commanded blade pitch angles (demand pitch angles), rad.
!REAL(ReKi), INTENT(IN )      :: TwrAccel                                        ! Tower Acceleration, m/s^2.
!REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!
!BlPitchCom = 0.0
!
!
!
!RETURN
!END SUBROUTINE PitchCntrl
!=======================================================================
!SUBROUTINE UserGen ( HSS_Spd, LSS_Spd, NumBl, ZTime, DT, GenEff, DelGenTrq, DirRoot, GenTrq, ElecPwr )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! generator torque and power model.  Modify this code to create your
!   ! own model.
!
!   ! NOTE: If you (the user) wants to switch on-or-off the generator DOF at
!   !       runtime from this user-defined routine, then do the following:
!   !          (1) USE MODULE DOFs().
!   !          (2) Type in "DOF_Flag(DOF_GeAz) = VALUE" where VALUE = .TRUE. or
!   !              .FALSE. depending on whether you want to turn-on or turn-off
!   !              the DOF, respectively.  Turning off the DOF forces the
!   !              current RATE to remain fixed.  If the rate is currently zero,
!   !              the current POSITION will remain fixed as well.
!   !       Note that this technique WILL NOT work for user-defined routines
!   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
!   !       preprocessor.
!
!
!USE                             Precision
!
!
!IMPLICIT                        NONE
!
!
!   ! Passed Variables:
!
!INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
!
!REAL(ReKi), INTENT(IN )      :: DelGenTrq                                       ! Pertubation in generator torque used during FAST linearization (zero otherwise), N-m.
!REAL(DbKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
!REAL(ReKi), INTENT(OUT)      :: ElecPwr                                         ! Electrical power (account for losses), watts.
!REAL(ReKi), INTENT(IN )      :: GenEff                                          ! Generator efficiency, (-).
!REAL(ReKi), INTENT(OUT)      :: GenTrq                                          ! Electrical generator torque, N-m.
!REAL(ReKi), INTENT(IN )      :: LSS_Spd                                         ! LSS speed, rad/s.
!REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
!REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!
!GenTrq  = 0.0 + DelGenTrq  ! Make sure to add the pertubation on generator torque, DelGenTrq.  This is used only for FAST linearization (it is zero otherwise).
!
!
!   ! The generator efficiency is either additive for motoring,
!   !   or subtractive for generating power.
!
!IF ( GenTrq > 0.0 )  THEN
!   ElecPwr = GenTrq*HSS_Spd*GenEff
!ELSE
!   ElecPwr = GenTrq*HSS_Spd/GenEff
!ENDIF
!
!
!
!RETURN
!END SUBROUTINE UserGen
!=======================================================================
SUBROUTINE UserHSSBr ( GenTrq, ElecPwr, HSS_Spd, NumBl, ZTime, DT, DirRoot, HSSBrFrac )


   ! This is a dummy routine for holding the place of a user-specified
   ! HSS brake model.  This routine must specify the fraction
   ! (HSSBrFrac) of full torque to be applied to the HSS by the HSS
   ! brake.  The magnitude of the full torque (HSSBrFrac = 1.0) equals
   ! HSSBrTqF from the primary input file.  Modify this code to create
   ! your own model.

   ! NOTE: If you (the user) wants to switch on-or-off the generator DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_GeAz) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).

REAL(DbKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
REAL(ReKi), INTENT(IN )      :: ElecPwr                                         ! Electrical power (account for losses), watts.
REAL(ReKi), INTENT(IN )      :: GenTrq                                          ! Electrical generator torque, N-m.
REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
REAL(ReKi), INTENT(OUT)      :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-).
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



HSSBrFrac = 0.0   ! NOTE: This must be specified as a real number between 0.0 (off - no brake torque) and 1.0 (full - max brake torque = HSSBrTqF); FAST/ADAMS will Abort otherwise.



RETURN
END SUBROUTINE UserHSSBr
!=======================================================================
SUBROUTINE UserPtfmLd ( X, XD, ZTime, DirRoot, PtfmAM, PtfmFt )


   ! This is a dummy routine for holding the place of a user-specified
   ! platform loading model.  Modify this code to create your own model.
   ! The local variables and associated calculations below provide a
   ! template for making this user-specified platform loading model
   ! include linear 6x6 damping and stiffness matrices.  These are
   ! provided as an example only and can be modified or deleted as
   ! desired by the user without detriment to the interface (i.e., they
   ! are not necessary for the interface).

   ! The platform loads returned by this routine should contain contributions
   !   from any external load acting on the platform other than loads
   !   transmitted from the wind turbine.  For example, these loads should
   !   contain contributions from foundation stiffness and damping [not
   !   floating] or mooring line restoring and damping [floating], as well as
   !   hydrostatic and hydrodynamic contributions [offshore].  The platform
   !   loads will be applied on the platform at the instantaneous platform
   !   reference position within FAST and ADAMS.

   ! This routine assumes that the platform loads are transmitted through a
   !   medium like soil [foundation] and/or water [offshore], so that added
   !   mass effects are important.  Consequently, the routine assumes that the
   !   total platform load can be written as:
   !
   ! PtfmF(i) = SUM( -PtfmAM(i,j)*XDD(j), j=1,2,..,6) + PtfmFt(i) for i=1,2,...,6
   !
   ! where,
   !   PtfmF(i)    = the i'th component of the total load applied on the
   !                 platform; positive in the direction of positive motion of
   !                 the i'th DOF of the platform
   !   PtfmAM(i,j) = the (i,j) component of the platform added mass matrix
   !                 (output by this routine)
   !   XDD(j)      = the j'th component of the platform acceleration vector
   !   PtfmFt(i)   = the i'th component of the portion of the platform load
   !                 associated with everything but the added mass effects;
   !                 positive in the direction of positive motion of the i'th
   !                 DOF of the platform (output by this routine)

   ! The order of indices in all arrays passed to and from this routine is as
   !   follows:
   !      1 = Platform surge / xi-component of platform translation (internal DOF index = DOF_Sg)
   !      3 = Platform sway  / yi-component of platform translation (internal DOF index = DOF_Sw)
   !      3 = Platform heave / zi-component of platform translation (internal DOF index = DOF_Hv)
   !      4 = Platform roll  / xi-component of platform rotation    (internal DOF index = DOF_R )
   !      5 = Platform pitch / yi-component of platform rotation    (internal DOF index = DOF_P )
   !      6 = Platform yaw   / zi-component of platform rotation    (internal DOF index = DOF_Y )

   ! NOTE: The added mass matrix returned by this routine, PtfmAM, must be
   !       symmetric.  FAST and ADAMS will abort otherwise.
   !
   !       Please also note that the hydrostatic restoring contribution to the
   !       hydrodynamic force returned by this routine should not contain the
   !       effects of body weight, as is often done in classical marine
   !       hydrodynamics.  The effects of body weight are included within FAST
   !       and ADAMS.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(OUT)      :: PtfmAM (6,6)                                    ! Platform added mass matrix, kg, kg-m, kg-m^2.
REAL(ReKi), INTENT(OUT)      :: PtfmFt   (6)                                    ! The 3 components of the portion of the platform force (in N  ) acting at the platform reference and the 3 components of the portion of the platform moment (in N-m  ) acting at the platform reference associated with everything but the added-mass effects; positive forces are in the direction of motion.
REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement    (in m  )        of the platform reference and the 3 components of the rotational displacement        (in rad  )        of the platform relative to the inertial frame.
REAL(ReKi), INTENT(IN )      :: XD       (6)                                    ! The 3 components of the translational velocity        (in m/s)        of the platform reference and the 3 components of the rotational (angular) velocity  (in rad/s)        of the platform relative to the inertial frame.
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.


   ! Local Variables:

REAL(ReKi)                   :: Damp   (6,6)                                    ! Damping matrix.
REAL(ReKi)                   :: Stff   (6,6)                                    ! Stiffness/restoring matrix.

INTEGER(4)                   :: I                                               ! Generic index.
INTEGER(4)                   :: J                                               ! Generic index.



Damp  (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

Stff  (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

PtfmAM(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

PtfmFt(1)   = 0.0
PtfmFt(2)   = 0.0
PtfmFt(3)   = 0.0
PtfmFt(4)   = 0.0
PtfmFt(5)   = 0.0
PtfmFt(6)   = 0.0

DO J = 1,6
   DO I = 1,6
      PtfmFt(I) = PtfmFt(I) - Damp(I,J)*XD(J) - Stff(I,J)*X(J)
   ENDDO
ENDDO



RETURN
END SUBROUTINE UserPtfmLd
!=======================================================================
SUBROUTINE UserRFrl ( RFrlDef, RFrlRate, ZTime, DirRoot, RFrlMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! rotor-furl spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the rotor-furl DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_RFrl) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the rotor-furl hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: RFrlDef                                         ! Rotor-furl angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: RFrlMom                                         ! Rotor-furl restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: RFrlRate                                        ! Rotor-furl angular rate, rad/s
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(*), INTENT(IN )    :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



RFrlMom = 0.0



RETURN
END SUBROUTINE UserRFrl
!=======================================================================
SUBROUTINE UserTeet ( TeetDef, TeetRate, ZTime, DirRoot, TeetMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! teeter spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the teeter DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_Teet) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the teeter hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: TeetDef                                         ! Rotor-teeter angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: TeetMom                                         ! Rotor-teeter restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: TeetRate                                        ! Rotor-teeter angular rate, rad/s
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TeetMom = 0.0



RETURN
END SUBROUTINE UserTeet
!=======================================================================
SUBROUTINE UserTFin ( TFrlDef , TFrlRate, ZTime   , DirRoot, &
                      TFinCPxi, TFinCPyi, TFinCPzi,          &
                      TFinCPVx, TFinCPVy, TFinCPVz,          &
                      TFinAOA , TFinQ   ,                    &
                      TFinCL  , TFinCD  ,                    &
                      TFinKFx , TFinKFy                        )


   ! This is a dummy routine for holding the place of user-specified
   ! computations for tail fin aerodynamic loads.  Modify this code to
   ! create your own logic.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(OUT)      :: TFinAOA                                         ! Angle-of-attack between the relative wind velocity and tail fin chordline, rad.
REAL(ReKi), INTENT(OUT)      :: TFinCD                                          ! Tail fin drag            coefficient resulting from current TFinAOA, (-).
REAL(ReKi), INTENT(OUT)      :: TFinCL                                          ! Tail fin lift            coefficient resulting from current TFinAOA, (-).
REAL(ReKi), INTENT(IN )      :: TFinCPVx                                        ! Absolute Velocity of the tail center-of-pressure along tail fin chordline pointing toward tail fin trailing edge, m/s.
REAL(ReKi), INTENT(IN )      :: TFinCPVy                                        ! Absolute Velocity of the tail center-of-pressure normal to plane of tail fin pointing towards suction surface   , m/s.
REAL(ReKi), INTENT(IN )      :: TFinCPVz                                        ! Absolute Velocity of the tail center-of-pressure in plane of tail fin normal to chordline and nominally upward  , m/s.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Improve the description of input arguments TFinCPxi, TFinCPyi, and
!jmj   TFinCPzi:
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPxi                                        ! Downwind distance from the inertial frame origin to the tail fin center-of-pressure, m.
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPyi                                        ! Lateral  distance from the inertial frame origin to the tail fin center-of-pressure, m.
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPzi                                        ! Vertical distance from the inertial frame origin to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPxi                                        ! Downwind distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPyi                                        ! Lateral  distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPzi                                        ! Vertical distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi), INTENT(OUT)      :: TFinKFx                                         ! Aerodynamic force  at the tail fin center-of-pressure (point K) along tail fin chordline pointing toward tail fin trailing edge, N.
REAL(ReKi), INTENT(OUT)      :: TFinKFy                                         ! Aerodynamic force  at the tail fin center-of-pressure (point K) normal to plane of tail fin pointing towards suction surface   , N.
REAL(ReKi), INTENT(OUT)      :: TFinQ                                           ! Dynamic pressure of the relative wind velocity, Pa.
REAL(ReKi), INTENT(IN )      :: TFrlDef                                         ! Tail-furl angular deflection, rad.
REAL(ReKi), INTENT(IN )      :: TFrlRate                                        ! Tail-furl angular rate, rad/s
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TFinAOA = 0.0
TFinCL  = 0.0
TFinCD  = 0.0
TFinQ   = 0.0
TFinKFx = 0.0
TFinKFy = 0.0



RETURN
END SUBROUTINE UserTFin
!=======================================================================
SUBROUTINE UserTFrl ( TFrlDef, TFrlRate, ZTime, DirRoot, TFrlMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! tail-furl spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the tail-furl DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_TFrl) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the tail-furl hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: TFrlDef                                         ! Tail-furl angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: TFrlMom                                         ! Tail-furl restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: TFrlRate                                        ! Tail-furl angular rate, rad/s
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TFrlMom = 0.0



RETURN
END SUBROUTINE UserTFrl
!=======================================================================
SUBROUTINE UserTwrLd ( JNode, X, XD, ZTime, DirRoot, TwrAM, TwrFt )

   ! This is a dummy routine for holding the place of a user-specified
   ! tower loading model.  Modify this code to create your own model.
   ! The local variables and associated calculations below provide a
   ! template for making this user-specified tower loading model
   ! include linear 6x6 damping and stiffness matrices.  These are
   ! provided as an example only and can be modified or deleted as
   ! desired by the user without detriment to the interface (i.e., they
   ! are not necessary for the interface).

   ! The tower loads returned by this routine should contain contributions from
   !   any external load acting on the current tower element (indicated by
   !   JNode) other than loads transmitted from tower aerodynamics.  For
   !   example, these tower forces should contain contributions from foundation
   !   stiffness and damping [not floating] or mooring line/guy wire restoring
   !   and damping, as well as hydrostatic and hydrodynamic contributions
   !   [offshore].

   ! This routine assumes that the tower loads are transmitted through a medium
   !   like soil [foundation] and/or water [offshore], so that added mass
   !   effects are important.  Consequently, the routine assumes that the total
   !   load per unit length on the current tower element can be written as:
   !
   ! TwrF(i) = SUM( -TwrAM(i,j)*XDD(j), j=1,2,..,6) + TwrFt(i) for i=1,2,...,6
   !
   ! where,
   !   TwrF(i)    = the i'th component of the total load per unit length
   !                applied on the current tower element; positive in the
   !                direction of positive motion of the i'th DOF of the current
   !                tower element
   !   TwrAM(i,j) = the (i,j) component of the tower added mass matrix per unit
   !                length (output by this routine)
   !   XDD(j)     = the j'th component of the current tower element
   !                acceleration vector
   !   TwrFt(i)   = the i'th component of the portion of the current tower
   !                element load per unit length associated with everything but
   !                the added mass effects; positive in the direction of
   !                positive motion of the i'th DOF of the current tower
   !                element (output by this routine)

   ! The order of indices in all arrays passed to and from this routine is as
   !   follows:
   !      1 = Current tower element surge / xi-component of translation
   !      3 = Current tower element sway  / yi-component of translation
   !      3 = Current tower element heave / zi-component of translation
   !      4 = Current tower element roll  / xi-component of rotation
   !      5 = Current tower element pitch / yi-component of rotation
   !      6 = Current tower element yaw   / zi-component of rotation

   ! NOTE: The added mass matrix returned by this routine, TwrAM, must be
   !       symmetric.  FAST and ADAMS will abort otherwise.
   !
   !       Please also note that the hydrostatic restoring contribution to the
   !       hydrodynamic force returned by this routine should not contain the
   !       effects of body weight, as is often done in classical marine
   !       hydrodynamics.  The effects of body weight are included within FAST
   !       and ADAMS.

USE                             NWTC_Library

IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(OUT)      :: TwrAM  (6,6)                                    ! Added mass matrix per unit length of current tower element, kg/m, kg-m/m, kg-m^2/m.
REAL(ReKi), INTENT(OUT)      :: TwrFt    (6)                                    ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force per unit length (in N/m) at the current tower element and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the tower moment per unit length (in N-m/m) acting at the current tower element associated with everything but the added-mass effects; positive forces are in the direction of motion.
REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement (in m  ) of the current tower node and the 3 components of the rotational displacement       (in rad  ) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
REAL(ReKi), INTENT(IN )      :: XD       (6)                                    ! The 3 components of the translational velocity     (in m/s) of the current tower node and the 3 components of the rotational (angular) velocity (in rad/s) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

INTEGER(4), INTENT(IN )      :: JNode                                           ! The number of the current tower node / element, (-). [1 to TwrNodes]

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.


   ! Local Variables:

REAL(ReKi)                   :: Damp   (6,6)                                    ! Damping matrix.
REAL(ReKi)                   :: Stff   (6,6)                                    ! Stiffness/restoring matrix.

INTEGER(4)                   :: I                                               ! Generic index.
INTEGER(4)                   :: J                                               ! Generic index.



Damp (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

Stff (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

TwrAM(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

TwrFt(1)   = 0.0
TwrFt(2)   = 0.0
TwrFt(3)   = 0.0
TwrFt(4)   = 0.0
TwrFt(5)   = 0.0
TwrFt(6)   = 0.0

DO J = 1,6
   DO I = 1,6
      TwrFt(I) = TwrFt(I) - Damp(I,J)*XD(J) - Stff(I,J)*X(J)
   ENDDO
ENDDO



RETURN
END SUBROUTINE UserTwrLd

!=======================================================================
!SUBROUTINE UserVSCont ( HSS_Spd, LSS_Spd, NumBl, ZTime, DT, GenEff, DelGenTrq, DirRoot, GenTrq, ElecPwr )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! variable-speed torque and power control model.  Modify this code to
!   ! create your own model.
!
!   ! NOTE: If you (the user) wants to switch on-or-off the generator DOF at
!   !       runtime from this user-defined routine, then do the following:
!   !          (1) USE MODULE DOFs().
!   !          (2) Type in "DOF_Flag(DOF_GeAz) = VALUE" where VALUE = .TRUE. or
!   !              .FALSE. depending on whether you want to turn-on or turn-off
!   !              the DOF, respectively.  Turning off the DOF forces the
!   !              current RATE to remain fixed.  If the rate is currently zero,
!   !              the current POSITION will remain fixed as well.
!   !       Note that this technique WILL NOT work for user-defined routines
!   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
!   !       preprocessor.
!
!
!USE                             Precision
!
!
!IMPLICIT                        NONE
!
!
!   ! Passed Variables:
!
!INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
!
!REAL(ReKi), INTENT(IN )      :: DelGenTrq                                       ! Pertubation in generator torque used during FAST linearization (zero otherwise), N-m.
!REAL(DbKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
!REAL(ReKi), INTENT(OUT)      :: ElecPwr                                         ! Electrical power (account for losses), watts.
!REAL(ReKi), INTENT(IN )      :: LSS_Spd                                         ! LSS speed, rad/s.
!REAL(ReKi), INTENT(IN )      :: GenEff                                          ! Generator efficiency, (-).
!REAL(ReKi), INTENT(OUT)      :: GenTrq                                          ! Electrical generator torque, N-m.
!REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
!REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!
!GenTrq  = 0.0 + DelGenTrq  ! Make sure to add the pertubation on generator torque, DelGenTrq.  This is used only for FAST linearization (it is zero otherwise).
!
!
!   ! The generator efficiency is either additive for motoring,
!   !   or subtractive for generating power.
!
!IF ( GenTrq > 0.0 )  THEN
!   ElecPwr = GenTrq*HSS_Spd*GenEff
!ELSE
!   ElecPwr = GenTrq*HSS_Spd/GenEff
!ENDIF
!
!
!
!RETURN
!END SUBROUTINE UserVSCont
!=======================================================================
SUBROUTINE UserYawCont ( YawPos, YawRate, WindDir, YawError, NumBl, ZTime, DT, DirRoot, YawPosCom, YawRateCom )


   ! This is a dummy routine for holding the place of a user-specified
   ! nacelle-yaw controller.  Modify this code to create your own device.


   ! As indicated, the yaw controller must always specify a command (demand)
   !   yaw angle, YawPosCom, AND command (demand) yaw rate, YawRateCom.
   !   Normally, you should correlate these commands so that the commanded yaw
   !   angle is the integral of the commanded yaw rate, or likewise, the
   !   commanded yaw rate is the derivative of the commanded yaw angle.  FAST
   !   WILL NOT compute these correlations for you and DOES NOT check to
   !   ensure that they are correlated.  In some situations, it is desirable to
   !   set one of the commands (either yaw angle OR yaw rate) to ZERO depending
   !   on the desired transfer function of FAST's built-in actuator model (see
   !   below for a discussion of FAST's built-in actuator model).  In general,
   !   the commanded yaw angle and rate SHOULD NEVER be defined independent of
   !   each other with BOTH commands NONZERO.


   ! The yaw controller's effect on the FAST model depends on whether or not
   !   the yaw DOF is enabled as follows:
   !
   ! YawDOF = False - If the yaw DOF is disabled, then the commanded yaw angle
   !                  and rate will be the ACTUAL yaw angle and yaw rate used
   !                  internally by FAST (in general, you should ensure these
   !                  are correlated).  In this case, any desired actuator
   !                  effects should be built within this controller.  Also in
   !                  this case, FAST WILL NOT compute the correlated yaw
   !                  acceleration, but assume that it is ZERO.  If the
   !                  commanded yaw rate is zero while the commanded yaw angle
   !                  is changing in time, then the yaw controller's effect
   !                  on yaw angle is the identical to routine PitchCntrl()'s
   !                  effect on pitch angle (i.e., routine PitchCntrl()
   !                  commands changes in pitch angle with no associated
   !                  changes in pitch rate or pitch acceleration).  For yaw
   !                  control, this situation should be avoided however, since
   !                  yaw-induced gyroscopic pitching loads on the turbine
   !                  brought about by the yaw rate may be significant.
   !
   ! YawDOF = True  - If the yaw DOF is enabled, then the commanded yaw angle
   !                  and rate, YawPosCom and YawRateCom, become the neutral
   !                  yaw angle, YawNeut, and neutral yaw rate, YawRateNeut, in
   !                  FAST's built-in second-order actuator model defined by
   !                  inputs YawSpr and YawDamp.


   ! Description of FAST's built-in actuator model:
   !
   ! In the time-domain, FAST's built-in actuator model is defined as follows:
   !
   ! YawIner*YawAccel + YawDamp*YawRate + YawSpr*YawPos
   !                             = YawDamp*YawRateNeut + YawSpr*YawNeut + YawTq
   !
   ! so that the transmitted torque is:
   !
   ! YawMom = YawSpr*( YawPos - YawNeut ) + YawDamp*( YawRate - YawRateNeut )
   !
   ! where,
   !   YawSpr      = nacelle-yaw spring constant (defined in FAST's primary
   !                 input file)
   !   YawDamp     = nacelle-yaw damping constant (defined in FAST's primary
   !                 input file)
   !   YawIner     = instantaneous inertia of the nacelle and rotor about the
   !                 yaw axis
   !   YawNeut     = the commanded (neutral) yaw angle = YawPosCom
   !   YawRateNeut = the commanded (neutral) yaw rate  = YawRateCom
   !   YawPos      = yaw angle (position)
   !   YawRate     = yaw rate
   !   YawAccel    = yaw acceleration
   !   YawTq       = torque about the yaw axis applied by external forces above
   !                 the yaw bearing, such as wind loading
   !   YawMom      = torque transmitted through the yaw bearing
   !
   ! If the commanded yaw angle and rate are correlated (so that the commanded
   !   yaw angle is the integral of the commanded yaw rate, or likewise, the
   !   commanded yaw rate is the derivative of the commanded yaw angle), then
   !   FAST's built-in second-order actuator model will have the following
   !   characteristic transfer function:
   !
   !               YawDamp*s + YawSpr             2*Zeta*OmegaN*s + OmegaN^2
   ! T(s) = -------------------------------- = --------------------------------
   !        YawIner*s^2 + YawDamp*s + YawSpr   s^2 + 2*Zeta*OmegaN*s + OmegaN^2
   !
   ! where,
   !   T(s)    = the transfer function of FAST's built-in 2nd order actuator
   !             model
   !   OmegaN  = SQRT(YawSpr/YawIner) = yaw actuator natural frequency
   !   Zeta    = YawDamp/(2*SQRT(YawSpr*YawIner)) = yaw actuator damping ratio
   !             in fraction of critical
   !
   ! If only the yaw angle is commanded, and YawRateCom is zeroed, then the
   !   charecteristic transfer function of FAST's built-in second-order
   !   actuator model simplifies to:
   !
   !                     YawSpr                            OmegaN^2
   ! T(s) = -------------------------------- = --------------------------------
   !        YawIner*s^2 + YawDamp*s + YawSpr   s^2 + 2*Zeta*OmegaN*s + OmegaN^2
   !
   ! If only the yaw rate is commanded, and YawPosCom is zeroed, then the
   !   charecteristic transfer function of FAST's built-in second-order
   !   actuator model simplifies to:
   !
   !                    YawDamp                         2*Zeta*OmegaN
   ! T(s) = -------------------------------- = --------------------------------
   !        YawIner*s^2 + YawDamp*s + YawSpr   s^2 + 2*Zeta*OmegaN*s + OmegaN^2


   ! NOTE: If you (the user) wants to switch on-or-off the yaw DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_Yaw) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF acts is like
   !              setting YawDOF to False.
   !       This technique is useful, for example, if the yaw bearing has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).

REAL(DbKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
REAL(ReKi), INTENT(IN )      :: WindDir                                         ! Current horizontal hub-height wind direction (positive about the zi-axis), rad.
REAL(ReKi), INTENT(IN )      :: YawError                                        ! Current nacelle-yaw error estimate (positve about the zi-axis), rad.
REAL(ReKi), INTENT(IN )      :: YawPos                                          ! Current nacelle-yaw angular position, rad.
REAL(ReKi), INTENT(OUT)      :: YawPosCom                                       ! Commanded nacelle-yaw angular position (demand yaw angle), rad.
REAL(ReKi), INTENT(IN )      :: YawRate                                         ! Current nacelle-yaw angular rate, rad/s.
REAL(ReKi), INTENT(OUT)      :: YawRateCom                                      ! Commanded nacelle-yaw angular rate (demand yaw rate), rad/s.
REAL(DbKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



YawPosCom  = 0.0
YawRateCom = 0.0



RETURN
END SUBROUTINE UserYawCont
!=======================================================================
