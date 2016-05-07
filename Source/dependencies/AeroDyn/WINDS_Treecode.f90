!**********************************************************************************************************************************    
!   Primary author:  Hans Johnston  (johnston@math.umass.edu) 
!   Department of Mathematics
!   University of Massachusetts, Amherst
!   February 2011 
!   Reference: http://www.math.umass.edu/~johnston/newtreecode.html
!    
!
!   Modified by Shujian Liu (shujian@umass.edu)
!   Last edited: 2015/08/06    
!    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE WINDS_Treecode   

   USE OMP_LIB      ! OpenMP
   USE NWTC_Library
   USE AeroDyn_Types
    
   IMPLICIT        NONE
 

      ! ..... Public Subroutines ............
   PUBLIC :: BiotSavart_Treecode  ! To accelerate Biot-Savart Law by Treecode Algorithn....
   
   
! global variables to track tree levels  
   INTEGER(IntKi) :: minlevel,maxlevel
   
! node pointer and node type declarations
   TYPE tnode_pointer
        TYPE(tnode), POINTER :: p_to_tnode
   END TYPE tnode_pointer
   TYPE tnode
        INTEGER(IntKi)   :: numpar,ibeg,iend
        REAL(DbKi)       :: x_min,y_min,z_min
        REAL(DbKi)       :: x_max,y_max,z_max
        REAL(DbKi)       :: x_mid,y_mid,z_mid
        REAL(DbKi)       :: radius,sqradius,aspect
        INTEGER(IntKi)   :: level,num_children,exist_ms
        REAL(DbKi),DIMENSION(:,:,:,:),POINTER :: ms
        TYPE(tnode_pointer), DIMENSION(8) :: child
   END TYPE tnode   
   

CONTAINS
!=================================================================================================================================   
SUBROUTINE BiotSavart_Treecode(fila_x, fila_y, fila_z, points, GAMMA_array, RC_array, Uind_array,  &
                               Parameters, OtherState, Discrete_states, ErrStat, ErrMess)

! INPUTS:
!       fila_x: x coordinates of source vortex filaments
!       fila_y: y coordinates of source vortex filaments
!       fila_z: z coordinates of source vortex filaments
!       points: Target/receiving particles
!       GAMMA_array: Strength of source vortex filaments
!       RC_array: Core radius of source vortex filaments
!       Parameters: parameters for Treecode Algorithm
!       OtherState: not used here
!       Discrete_states: not used here
!
! OOTPUTS:
!       Uind_array: induced velocity for all the target/receiving particles
!       ErrStat, ErrMess: error information
!
!....................................................................
 
   IMPLICIT                        NONE
   
      ! Passed variables
   REAL(DbKi), INTENT(INOUT), DIMENSION(:,:)        :: fila_x, fila_y, fila_z 
   REAL(DbKi), INTENT(IN   ), DIMENSION(:,:)        :: points
   REAL(DbKi), INTENT(INOUT), DIMENSION(:,:)        :: GAMMA_array
   REAL(DbKi), INTENT(INOUT), DIMENSION(:,:)        :: RC_array
   REAL(DbKi), INTENT(INOUT), DIMENSION(:,:)        :: Uind_array   


   
   TYPE(AD_ParameterType),        INTENT(IN   )    :: Parameters       ! Parameters
   TYPE(AD_OtherStateType),       INTENT(IN   )    :: OtherState       ! Other/optimization states
   TYPE(AD_DiscreteStateType),    INTENT(IN   )    :: Discrete_states  ! Discrete states at t
   INTEGER(IntKi),                INTENT(  OUT)    :: ErrStat        ! The error status code
   CHARACTER(*),                  INTENT(  OUT)    :: ErrMess         ! The error message, if an error occurred


   ! Parameters for Treecode
      REAL(DbKi)     :: theta
      INTEGER(IntKi) :: order
      INTEGER(IntKi) :: maxparnode
      REAL(DbKi)     :: dist_tol
      REAL(DbKi)     :: delta        

      
   ! local variables
      REAL(DbKi),ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2
      
! (formerly global) variables for taylor expansions

      INTEGER(IntKi) :: torder,torderlim,orderoffset

! (formerly global) variables used when computing potential/force

      REAL(DbKi) :: thetasq,gdist_tolsq 

! OpenMP variables

      INTEGER(IntKi) :: nthreads,myid,chunk_size,ibeg,iend       
            
! local variables

      TYPE(tnode),POINTER :: troot
      INTEGER(IntKi) :: i  ! Counter of OpenMP cores
      INTEGER(IntKi) :: level,err
      REAL(DbKi), DIMENSION(6) :: xyzminmax ! x,y,z min and max
      REAL(DbKi) :: t1 ! temporary variable for Taylor coeffs
      INTEGER(IntKi) :: j,k ! shujian
      
   ! runtime parameters

      INTEGER(IntKi) :: numfila, numtars
         
   ! arrays for coordinates, charge, force calculations and permutation info
            
      INTEGER(IntKi),ALLOCATABLE,DIMENSION(:) :: perm
            
    ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMess  = ""      
                  
      
   ! Set parameters    
      theta      = Parameters%FVM%Tree_Parms%theta  
      order      = Parameters%FVM%Tree_Parms%order   
      maxparnode = Parameters%FVM%Tree_Parms%maxparnode 
      dist_tol   = Parameters%FVM%Tree_Parms%dist_tol 
      delta      = Parameters%FVM%Tree_Parms%delta         
      
      numtars  = SIZE(points, 1)
      numfila  = SIZE(fila_x, 1)
            
      
   ! (formerly global) integers & reals: TORDER, TORDERLIM and THETASQ

      torder = order
      torderlim = torder+1
      
      thetasq = theta*theta
      gdist_tolsq = dist_tol*dist_tol
      
   ! allocate global Taylor expansion variables

      ALLOCATE(cf(0:torder),cf1(torder+1),cf2(torder+1),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error (in WInDS): Error allocating Taylor variables! (In Treecode Algorithm)'
         STOP
      END IF

   ! initialize arrays for Taylor sums and coeffs

      DO i=0,torder
         cf(i)=-(REAL(i)+1.0_DbKi)
      END DO
      DO i=1,torderlim
         t1=1.0_DbKi/REAL(i)
         cf1(i)=1.0_DbKi-0.5_DbKi*t1
         cf2(i)=1.0_DbKi-t1
      END DO     

      xyzminmax(1)=MINVAL(fila_x(1:numfila,3))
      xyzminmax(2)=MAXVAL(fila_x(1:numfila,3))
      xyzminmax(3)=MINVAL(fila_y(1:numfila,3)) 
      xyzminmax(4)=MAXVAL(fila_y(1:numfila,3)) 
      xyzminmax(5)=MINVAL(fila_z(1:numfila,3))
      xyzminmax(6)=MAXVAL(fila_z(1:numfila,3))

   ! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(troot)  
      
      
   ! set global variables to track tree levels during construction

      level=0
      minlevel=50000
      maxlevel=0
      

     
      !..............................................................................................................................
      ! Create the tree structure
      ALLOCATE(perm(numfila))
      CALL CREATE_TREE(troot, 1, numfila, fila_x, fila_y, fila_z, GAMMA_array, RC_array, perm, torder, maxparnode, xyzminmax, level, numfila)
    
        
      !..............................................................................................................................
      ! Compute the induced velocity      
      nthreads = omp_get_max_threads()         
      IF (Parameters%FVM%Tree_Parms%cores > 0 .AND. Parameters%FVM%Tree_Parms%cores < nthreads) THEN
          nthreads = Parameters%FVM%Tree_Parms%cores
      END IF                             
      !nthreads = 4   ! sliu: for debug   
      chunk_size = (numtars+nthreads-1)/nthreads          
      
      !$OMP PARALLEL DO PRIVATE(myid,ibeg,iend,i)   &
      !$OMP             SHARED(troot,fila_x,fila_y,fila_z,GAMMA_array,RC_array,Uind_array,numtars,  &
      !$OMP                 cf,cf1,cf2,torder,torderlim,thetasq,gdist_tolsq, &
      !$OMP                 points,numfila,chunk_size, nthreads)   
         do i = 0,nthreads-1
            myid = omp_get_thread_num()   
            ! myid = Parameters%FVM%OpenMP_Parms%OpenMPCores  ! sliu: bug
            ibeg = myid*chunk_size+1
            iend = MIN0((myid+1)*chunk_size,numtars)
          
            CALL TREE_COMPVel(troot, points, numtars,ibeg,iend, &
                              cf,cf1,cf2,torder,torderlim,       &
                              thetasq,gdist_tolsq,               &
                              fila_x, fila_y, fila_z, GAMMA_array, RC_array, Uind_array, numfila)
         end do             
      !$OMP END PARALLEL DO       
        
      
         
         
      !..............................................................................................................................
      ! Call CLEANUP to deallocate global variables and tree structure.
      CALL CLEAN_UP(troot)      
      IF ( ALLOCATED( perm ) )         DEALLOCATE( perm )

   
   
END SUBROUTINE BiotSavart_Treecode   
    
!=================================================================================================================================                               
RECURSIVE SUBROUTINE CREATE_TREE(p, ibeg, iend, fila_x, fila_y, fila_z, g, rc, perm, torder, maxparnode, xyzmm, level, numfila)

! CREATE_TREE recursively create the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box.   
!............................................................................................................................

      IMPLICIT NONE

      TYPE(tnode), POINTER,             INTENT(INOUT) :: p
      INTEGER(IntKi),                   INTENT(IN   ) :: ibeg, iend, maxparnode, level, numfila, torder
      INTEGER(IntKi),DIMENSION(numfila),INTENT(INOUT) :: perm      
      REAL(DbKi), DIMENSION(:,:),       INTENT(INOUT) :: fila_x
      REAL(DbKi), DIMENSION(:,:),       INTENT(INOUT) :: fila_y 
      REAL(DbKi), DIMENSION(:,:),       INTENT(INOUT) :: fila_z 
      REAL(DbKi), DIMENSION(:,:),       INTENT(INOUT) :: g, rc
      REAL(DbKi), DIMENSION(6),         INTENT(IN   ) :: xyzmm

! local variables

      REAL(DbKi)                     :: x_mid,y_mid,z_mid,xl,yl,zl,lmax,t1,t2,t3
      INTEGER(IntKi), DIMENSION(8,2) :: ind
      REAL(DbKi),     DIMENSION(6,8) :: xyzmms
      INTEGER(IntKi)                 :: i,j,limin,limax,err,loclev,numposchild
      REAL(DbKi),     DIMENSION(6)   :: lxyzmm
     
! allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*)   ' Error (in WInDS): Error allocating pointer! (In Treecode Algorithm)'
         RETURN              
      END IF

! set node fields: number of particles, exist_ms
! and xyz bounds 

      p%numpar=iend-ibeg+1
      p%exist_ms=0

      ! shrink the box
      p%x_min=MINVAL(fila_x(ibeg:iend,3)) ! fila_x(ibeg:iend,3) middle points of the filaments...
      p%x_max=MAXVAL(fila_x(ibeg:iend,3)) 
      p%y_min=MINVAL(fila_y(ibeg:iend,3)) 
      p%y_max=MAXVAL(fila_y(ibeg:iend,3)) 
      p%z_min=MINVAL(fila_z(ibeg:iend,3)) 
      p%z_max=MAXVAL(fila_z(ibeg:iend,3)) 


! compute aspect ratio

      xl=p%x_max-p%x_min
      yl=p%y_max-p%y_min
      zl=p%z_max-p%z_min

      lmax=MAX(xl,yl,zl)
      t1=lmax
      t2=MIN(xl,yl,zl)
      IF (t2 .NE. 0.0_DbKi) THEN
         p%aspect=t1/t2
      ELSE
         p%aspect=0.0_DbKi
      END IF

! midpoint coordinates , RADIUS and SQRADIUS 

      p%x_mid=(p%x_max+p%x_min)/2.0_DbKi
      p%y_mid=(p%y_max+p%y_min)/2.0_DbKi
      p%z_mid=(p%z_max+p%z_min)/2.0_DbKi
      t1=p%x_max-p%x_mid
      t2=p%y_max-p%y_mid
      t3=p%z_max-p%z_mid
      p%sqradius=t1*t1+t2*t2+t3*t3
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%ibeg=ibeg
      p%iend=iend
      p%level=level
      IF (maxlevel .LT. level) THEN
         maxlevel=level
      END IF
      p%num_children=0
      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO
      
! compute moments!  all threads will share tree

      IF (p%exist_ms .EQ. 0) THEN
         ALLOCATE(p%ms(0:torder,0:torder,0:torder,3),STAT=err)
          IF (err .NE. 0) THEN
              WRITE(6,*) 'Error (in WInDS): Error allocating node moments! (In Treecode Algorithm)'
              STOP
          END IF
          CALL COMP_MS(p,fila_x,fila_y,fila_z,g,torder,numfila)
          p%exist_ms=1
       END IF   

      IF (p%numpar .GT. maxparnode) THEN

! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1

         xyzmms(1,1)=p%x_min
         xyzmms(2,1)=p%x_max
         xyzmms(3,1)=p%y_min
         xyzmms(4,1)=p%y_max
         xyzmms(5,1)=p%z_min
         xyzmms(6,1)=p%z_max
         ind(1,1)=ibeg
         ind(1,2)=iend
         x_mid=p%x_mid
         y_mid=p%y_mid
         z_mid=p%z_mid

         CALL PARTITION_8(fila_x,fila_y,fila_z,g,rc,perm,xyzmms,xl,yl,zl,lmax,    &
                          numposchild,x_mid,y_mid,z_mid,ind,numfila)

! create children if indicated and store info in parent

         loclev=level+1
         DO i=1,numposchild
            IF (ind(i,1) .LE. ind(i,2)) THEN
               p%num_children=p%num_children+1
               lxyzmm=xyzmms(:,i)
               CALL CREATE_TREE(p%child(p%num_children)%p_to_tnode,           & 
                              ind(i,1),ind(i,2),fila_x,fila_y,fila_z,g,rc,perm,  &
                              torder, maxparnode,lxyzmm,loclev,numfila)
            END IF
         END DO
      ELSE
         IF (level .LT. minlevel) THEN
            minlevel=level
         END IF
      END IF   

      
CONTAINS
   !........................................................................................................................
    SUBROUTINE PARTITION_8(fila_x,fila_y,fila_z,GAMMA_array,RC_array,perm,xyzmms,xl,yl,zl,lmax,    &
                           numposchild,x_mid,y_mid,z_mid,ind,numfila)

! PARTITION_8 determines the particle indices of the eight sub boxes
! containing the particles after the box defined by particles I_BEG
! to I_END is divided by its midpoints in each coordinate direction.
! The determination of the indices is accomplished by the subroutine
! PARTITION. A box is divided in a coordinate direction as long as the
! resulting aspect ratio is not too large. This avoids the creation of
! "narrow" boxes in which Talyor expansions may become inefficient.
! On exit the INTEGER array IND (dimension 8 x 2) contains
! the indice limits of each new box (node) and NUMPOSCHILD the number 
! of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
! that box J is empty.


    IMPLICIT NONE

      INTEGER(IntKi),                   INTENT(IN   ) :: numfila
      REAL(DbKi),DIMENSION(numfila,3),  INTENT(INOUT) :: fila_x,fila_y,fila_z
      REAL(DbKi),DIMENSION(numfila,1),  INTENT(INOUT) :: GAMMA_array, RC_array
      INTEGER(IntKi),DIMENSION(numfila),INTENT(INOUT) :: perm   
      INTEGER(IntKi), DIMENSION(8,2),   INTENT(INOUT) :: ind
      REAL(DbKi),DIMENSION(6,8),        INTENT(INOUT) :: xyzmms
      REAL(DbKi),                       INTENT(IN   ) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax
      INTEGER(IntKi),                   INTENT(INOUT) :: numposchild

! local variables

      INTEGER(IntKi) :: temp_ind,i
      REAL(DbKi) :: critlen

      numposchild=1
      critlen=lmax/sqrt(2.0_DbKi)

      IF (xl .GE. critlen) THEN
         CALL PARTITION(fila_x,fila_y,fila_z,GAMMA_array,RC_array,perm,ind(1,1),ind(1,2),      &
                        x_mid,temp_ind,numfila)

         ind(2,1)=temp_ind+1
         ind(2,2)=ind(1,2)
         ind(1,2)=temp_ind
         xyzmms(:,2)=xyzmms(:,1)
         xyzmms(2,1)=x_mid
         xyzmms(1,2)=x_mid
         numposchild=2*numposchild
      END IF 
 
      IF (yl .GE. critlen) THEN
         DO i=1,numposchild  
            CALL PARTITION(fila_y,fila_x,fila_z,GAMMA_array,RC_array,perm,ind(i,1),ind(i,2),     & 
                           y_mid,temp_ind,numfila)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(4,i)=y_mid
            xyzmms(3,numposchild+i)=y_mid
         END DO
         numposchild=2*numposchild
      END IF

      IF (zl .GE. critlen) THEN
         DO i=1,numposchild
            CALL PARTITION(fila_z,fila_x,fila_y,GAMMA_array,RC_array,perm,ind(i,1),ind(i,2),     & 
                           z_mid,temp_ind,numfila)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(6,i)=z_mid
            xyzmms(5,numposchild+i)=z_mid
         END DO
         numposchild=2*numposchild
      END IF

      RETURN 
    END SUBROUTINE PARTITION_8
   !........................................................................................................................  
   SUBROUTINE PARTITION(a,b,c,g,rc,perm,ibeg,iend,val,midind,numfila)

   ! PARTITION determines the index MIDIND, after partitioning
   ! in place the  arrays A,B,C and G,  such that 
   ! A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL. 
   ! If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
   ! is returned as IBEG-1. 

      IMPLICIT NONE
 
      INTEGER(IntKi),                   INTENT(IN   ) :: numfila,ibeg,iend
      INTEGER(IntKi),DIMENSION(numfila),INTENT(INOUT) :: perm
      REAL(DbKi),DIMENSION(numfila,3),  INTENT(INOUT) :: a,b,c
      REAL(DbKi),DIMENSION(numfila,1),  INTENT(INOUT) :: g,rc
      INTEGER(IntKi),                   INTENT(INOUT) :: midind   
      REAL(DbKi),                       INTENT(IN   ) :: val

! local variables

      REAL(DbKi),DIMENSION(1,3)  :: temp_a, temp_b, temp_c
      REAL(DbKi)                 :: temp_g, temp_rc
      INTEGER(IntKi)             :: lower, upper, tperm

      IF (ibeg .LT. iend) THEN

! temporarily store IBEG entries and set A(IBEG)=VAL for 
! the partitoning algorithm.  

         temp_a(1,1:3)=a(ibeg,1:3)
         temp_b(1,1:3)=b(ibeg,1:3)
         temp_c(1,1:3)=c(ibeg,1:3)
         temp_g = g(ibeg,1)
         temp_rc = rc(ibeg,1)
         tperm = perm(ibeg)
         a(ibeg,3)=val 
         upper=ibeg
         lower=iend

         DO WHILE (upper .NE. lower)
            DO WHILE ((upper .LT. lower) .AND. (val .LT. a(lower,3)))
                  lower=lower-1
            END DO
            IF (upper .NE. lower) THEN
               a(upper,1:3) = a(lower,1:3)
               b(upper,1:3) = b(lower,1:3)
               c(upper,1:3) = c(lower,1:3)
               g(upper,1)   = g(lower,1)
               rc(upper,1)  = rc(lower,1)
               perm(upper)  = perm(lower)
            END IF
            DO WHILE ((upper .LT. lower) .AND. (val .GE. a(upper,3)))
                  upper=upper+1
            END DO
            IF (upper .NE. lower) THEN
               a(lower,1:3)= a(upper,1:3)
               b(lower,1:3)= b(upper,1:3)
               c(lower,1:3)= c(upper,1:3)
               g(lower,1)  = g(upper,1)
               rc(lower,1) = rc(upper,1)
               perm(lower) = perm(upper)
            END IF
         END DO
         midind=upper

! replace temp_a in position UPPER and change MIDIND if temp_a > VAL 

         IF (temp_a(1,3) .GT. val) THEN
            midind=upper-1
         END IF
         a(upper,1:3)= temp_a(1,1:3)
         b(upper,1:3)= temp_b(1,1:3)
         c(upper,1:3)= temp_c(1,1:3)
         g(upper,1)  = temp_g
         rc(upper,1) = temp_rc
         perm(upper) = tperm

      ELSEIF (ibeg .EQ. iend) THEN
         IF (a(ibeg,3) .LE. val) THEN
            midind=ibeg
         ELSE
            midind=ibeg-1
         END IF
      ELSE
         midind=ibeg-1
      END IF

   !RETURN
   END SUBROUTINE PARTITION
   !======================================================================================
   SUBROUTINE COMP_MS(p, fila_x, fila_y, fila_z, g, torder, numfila)
   IMPLICIT NONE

   ! COMP_MS computes the moments for node P needed in the Taylor approximation

      INTEGER(IntKi),                 INTENT(IN   ) :: torder,numfila
      TYPE(tnode),POINTER,            INTENT(INOUT) :: p 
      REAL(DbKi),DIMENSION(numfila,3),INTENT(IN   ) :: fila_x, fila_y, fila_z
      REAL(DbKi),DIMENSION(numfila,1),INTENT(IN   ) :: g

   ! local variables

      INTEGER(IntKi)          :: i,k1,k2,k3
      INTEGER(IntKi)          :: j, j1
      REAL(DbKi)              :: dx,dy,dz,tx,ty,tz
      INTEGER(IntKi)          :: simpsons
      REAL(DbKi)              :: wx,wy,wz
      REAL(DbKi),DIMENSION(3) :: dx_sim,dy_sim,dz_sim,wx_sim,wy_sim,wz_sim
     
      simpsons = 2  ! Need to be tuned, can be: 0, 2, 4, 8...
      p%ms = 0.0_DbKi
      
      
      IF (simpsons==0) THEN      
         DO i=p%ibeg,p%iend
            dx=p%x_mid - fila_x(i,3)
            dy=p%y_mid - fila_y(i,3)
            dz=p%z_mid - fila_z(i,3)
            
            wx = g(i,1)*(fila_x(i,2)-fila_x(i,1))/4/pi
            wy = g(i,1)*(fila_y(i,2)-fila_y(i,1))/4/pi
            wz = g(i,1)*(fila_z(i,2)-fila_z(i,1))/4/pi           
            
            tz=1.0_DbKi
            DO k3=0,torder
               ty=1.0_DbKi
               DO k2=0,torder-k3
                  tx=1.0_DbKi 
                  DO k1=0,torder-k3-k2
                     p%ms(k1,k2,k3,1)=p%ms(k1,k2,k3,1)+tx*ty*tz*wx
                     p%ms(k1,k2,k3,2)=p%ms(k1,k2,k3,2)+tx*ty*tz*wy
                     p%ms(k1,k2,k3,3)=p%ms(k1,k2,k3,3)+tx*ty*tz*wz
                     tx=tx*dx
                  END DO
                  ty=ty*dy
               END DO
               tz=tz*dz
            END DO
         END DO
         
      ELSE  ! Use Simpsons' rule
         DO i=p%ibeg,p%iend
            DO j=0, simpsons-1
               dx_sim(1)=p%x_mid - (fila_x(i,1) + (fila_x(i,2)-fila_x(i,1))*j/simpsons)
               dy_sim(1)=p%y_mid - (fila_y(i,1) + (fila_y(i,2)-fila_y(i,1))*j/simpsons)
               dz_sim(1)=p%z_mid - (fila_z(i,1) + (fila_z(i,2)-fila_z(i,1))*j/simpsons)
               
               dx_sim(2)=p%x_mid - (fila_x(i,1) + (fila_x(i,2)-fila_x(i,1))*(j+0.5)/simpsons)
               dy_sim(2)=p%y_mid - (fila_y(i,1) + (fila_y(i,2)-fila_y(i,1))*(j+0.5)/simpsons)
               dz_sim(2)=p%z_mid - (fila_z(i,1) + (fila_z(i,2)-fila_z(i,1))*(j+0.5)/simpsons)
               
               dx_sim(3)=p%x_mid - (fila_x(i,1) + (fila_x(i,2)-fila_x(i,1))*(j+1)/simpsons)
               dy_sim(3)=p%y_mid - (fila_y(i,1) + (fila_y(i,2)-fila_y(i,1))*(j+1)/simpsons)
               dz_sim(3)=p%z_mid - (fila_z(i,1) + (fila_z(i,2)-fila_z(i,1))*(j+1)/simpsons)
               
               wx_sim(1) = g(i,1)*(fila_x(i,2)-fila_x(i,1))/(4*pi)/(6*simpsons)
               wy_sim(1) = g(i,1)*(fila_y(i,2)-fila_y(i,1))/(4*pi)/(6*simpsons)
               wz_sim(1) = g(i,1)*(fila_z(i,2)-fila_z(i,1))/(4*pi)/(6*simpsons)                  
          
               wx_sim(2) = wx_sim(1) * 4
               wy_sim(2) = wy_sim(1) * 4
               wz_sim(2) = wz_sim(1) * 4                
          
               wx_sim(3) = wx_sim(1)
               wy_sim(3) = wy_sim(1)
               wz_sim(3) = wz_sim(1)            

               DO j1 = 1,3
                  tz=1.0_DbKi
                  DO k3=0,torder
                     ty=1.0_DbKi
                     DO k2=0,torder-k3
                        tx=1.0_DbKi 
                        DO k1=0,torder-k3-k2
                           p%ms(k1,k2,k3,1)=p%ms(k1,k2,k3,1)+tx*ty*tz*wx_sim(j1)
                           p%ms(k1,k2,k3,2)=p%ms(k1,k2,k3,2)+tx*ty*tz*wy_sim(j1)
                           p%ms(k1,k2,k3,3)=p%ms(k1,k2,k3,3)+tx*ty*tz*wz_sim(j1)
                           tx=tx*dx_sim(j1)
                        END DO
                        ty=ty*dy_sim(j1)
                     END DO
                     tz=tz*dz_sim(j1)
                  END DO
               END DO
            END DO    
         END DO          
      END IF
      
         
   !RETURN
   END SUBROUTINE COMP_MS
   !...................................................................................... 

   
END SUBROUTINE CREATE_TREE  
!=================================================================================================================================                               
 SUBROUTINE TREE_COMPVel(p, points, numtars,                &
                         pnum_beg,pnum_end,     &
                         cf,cf1,cf2,torder,torderlim,    &
                         thetasq,gdist_tolsq,            &
                         fila_x, fila_y, fila_z, g, rc, vel, numfila)

! TREE_COMPF is the driver routine which calls COMPF_TREE for each target
! particle, setting the global variables TARPOS before the call. 
! P is the root node of the tree. 
      IMPLICIT NONE

      
      INTEGER,INTENT(IN) :: numfila,numtars, &
                            pnum_beg,pnum_end,torder,torderlim
      TYPE(tnode),POINTER,INTENT(IN   )             :: p  
      REAL(DbKi),DIMENSION(numtars,3),INTENT(IN   ) :: points
      REAL(DbKi),DIMENSION(numfila,3),INTENT(IN   ) :: fila_x, fila_y, fila_z
      REAL(DbKi),DIMENSION(numfila,1),INTENT(IN   ) :: g, rc
      REAL(DbKi),DIMENSION(numtars,3),INTENT(INOUT) :: vel
      REAL(DbKi),INTENT(IN)  :: cf(0:torder),cf1(torder+1),      &
                                cf2(torder+1),thetasq,gdist_tolsq
     
! local variables
      REAL(DbKi),ALLOCATABLE,DIMENSION(:,:,:) :: b1
      INTEGER(IntKi)          :: i,j,err
      REAL(DbKi),DIMENSION(3) :: vel_local, tarpos
      
      ALLOCATE(b1(0:torderlim,0:torderlim,0:torderlim), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error (in WInDS): Error allocating b1 in TREE_COMPFP! (In Treecode Algorithm)'
         STOP
      END IF      
      
      

      IF (p%num_children .EQ. 0) THEN
          DO i = pnum_beg,pnum_end
             tarpos(1)=points(i,1)
             tarpos(2)=points(i,2)
             tarpos(3)=points(i,3)  
             vel_local  = 0.0_DbKi
             CALL COMVel_TREE(p, tarpos, vel_local, fila_x, fila_y, fila_z, g, rc, cf,cf1,cf2,b1,     &
                             thetasq,gdist_tolsq, torder, torderlim, numfila) 
             vel(i,1:3) = vel_local
          END DO   
      ELSE
         DO i = pnum_beg,pnum_end 
            tarpos(1)=points(i,1)
            tarpos(2)=points(i,2)
            tarpos(3)=points(i,3)    
            vel(i,1:3) = 0.0_DbKi
            vel_local  = 0.0_DbKi
            DO j=1,p%num_children
               CALL COMVel_TREE(p%child(j)%p_to_tnode, tarpos, vel_local, fila_x, fila_y, fila_z, g, rc,      &
                               cf,cf1,cf2,b1,           &
                               thetasq,gdist_tolsq,     &
                               torder,torderlim,numfila)
            
               vel(i,1:3) = vel(i,1:3) + vel_local            
            END DO
         END DO   
      END IF


      
CONTAINS
   !======================================================================================
   RECURSIVE SUBROUTINE COMVel_TREE(p, tarpos,vel,fila_x, fila_y, fila_z, g, rc,  &
                                    cf,cf1,cf2,              &
                                    b1,thetasq,gdist_tolsq,  &
                                    torder,torderlim,numfila)
   IMPLICIT NONE

      INTEGER,                        INTENT(IN   ) :: torder,torderlim
      INTEGER(IntKi),                 INTENT(IN   ) :: numfila
      TYPE(tnode),POINTER,            INTENT(IN   ) :: p  
      REAL(DbKi),DIMENSION(3),        INTENT(IN   ) :: tarpos 
      REAL(DbKi),DIMENSION(3),        INTENT(  OUT) :: vel
      REAL(DbKi),DIMENSION(numfila,3),INTENT(IN   ) :: fila_x, fila_y, fila_z
      REAL(DbKi),DIMENSION(numfila,1),INTENT(IN   ) :: g, rc
      REAL(DbKi),INTENT(IN) :: cf(0:torder),cf1(torder+1),  &
                               cf2(torder+1),thetasq,gdist_tolsq
      REAL(DbKi),DIMENSION(0:torderlim,0:torderlim,0:torderlim), INTENT(INOUT) :: b1
     
! local variables

      REAL(DbKi),DIMENSION(3) :: vel_local,vect_a
      REAL(DbKi) :: tx,ty,tz,distsq,t1,t2,cfj,cfk
      INTEGER(IntKi) :: i,j,k,err

! determine DISTSQ for MAC test

      tx=tarpos(1)-p%x_mid
      ty=tarpos(2)-p%y_mid
      tz=tarpos(3)-p%z_mid
      distsq=tx*tx+ty*ty+tz*tz

! intialize induced velocity

      vel=0.0_DbKi


! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.

      IF ((p%sqradius .LT. distsq*thetasq) .AND.  (p%sqradius .NE. 0.0_DbKi)) THEN
         vel_local=0.0_DbKi   
         
         CALL COMP_TCOEFF(p,tarpos,cf,cf1,cf2,b1,torder,torderlim)
         DO k=0,torder
            !cfk=cf(k)
            DO j=0,torder-k
               !cfj=cf(j)
               DO i=0,torder-k-j
                  vect_a(1) = (i+1)* b1(i+1,j,k)   ! i+1 = cf(i)
                  vect_a(2) = (j+1)* b1(i,j+1,k)
                  vect_a(3) = (k+1)* b1(i,j,k+1)
                   
                  vel_local(1) = vel_local(1) + vect_a(2)*p%ms(i,j,k,3) - vect_a(3)*p%ms(i,j,k,2) 
                  vel_local(2) = vel_local(2) + vect_a(3)*p%ms(i,j,k,1) - vect_a(1)*p%ms(i,j,k,3) 
                  vel_local(3) = vel_local(3) + vect_a(1)*p%ms(i,j,k,2) - vect_a(2)*p%ms(i,j,k,1) 
               END DO
            END DO
         END DO
         vel = vel + vel_local
      
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.

         IF (p%num_children .EQ. 0) THEN
            vel_local=0.0_DbKi 
            CALL COMVel_DIRECT(tarpos,vel_local,p%ibeg,p%iend, fila_x, fila_y, fila_z, g, rc, numfila)
            vel = vel_local
         ELSE
            DO i=1,p%num_children
               vel_local = 0.0_DbKi 
               CALL COMVel_TREE(p%child(i)%p_to_tnode, tarpos, vel_local, fila_x, fila_y, fila_z, g, rc, &
                                cf,cf1,cf2, b1, thetasq, gdist_tolsq, torder, torderlim, numfila)               
               vel = vel + vel_local
            END DO  
         END IF 
      END IF

      !RETURN
   END SUBROUTINE COMVel_TREE    
   !======================================================================================
   SUBROUTINE COMP_TCOEFF(p,tarpos,cf,cf1,cf2,b1,torder,torderlim)
      IMPLICIT NONE

   ! COMP_TCOEFF computes the Taylor coefficients of the potential
   ! using a recurrence formula.  The center of the expansion is the
   ! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.  

      TYPE(tnode),POINTER,    INTENT(IN   ) :: p      
      REAL(DbKi),DIMENSION(3),INTENT(IN   ) :: tarpos
      INTEGER(IntKi),         INTENT(IN   ) :: torder,torderlim     
      REAL(DbKi),INTENT(IN) :: cf(0:torder),cf1(torder+1), cf2(torder+1)
      REAL(DbKi),DIMENSION(0:torderlim,0:torderlim,0:torderlim), INTENT(INOUT) :: b1
     
   ! local varaibles

      REAL(DbKi)     :: dx,dy,dz,tdx,tdy,tdz,fac,sqfac,t1
      INTEGER(IntKi) :: i,j,k,tp1
      
      !REAL(DbKi)     ::  deltasq ! sliu
   
    ! setup variables

      dx=tarpos(1)-p%x_mid
      dy=tarpos(2)-p%y_mid
      dz=tarpos(3)-p%z_mid
      tdx=2.0_DbKi*dx
      tdy=2.0_DbKi*dy
      tdz=2.0_DbKi*dz
      !deltasq=0.0001_DbKi**2
      !fac=-1.0_DbKi/(dx*dx+dy*dy+dz*dz+deltasq)   
      fac=-1.0_DbKi/(dx*dx+dy*dy+dz*dz)   ! sliu: + deltasq
      sqfac=SQRT(-fac)

   ! 0th coeff or function val 

      b1(0,0,0)=sqfac

   ! 2 indices are 0
   
      b1(0,0,1)=fac*dz*sqfac
      b1(0,1,0)=fac*dy*sqfac
      b1(1,0,0)=fac*dx*sqfac
  
      DO i=2,torderlim
         b1(0,0,i)=fac*(tdz*cf1(i)*b1(0,0,i-1)+cf2(i)*b1(0,0,i-2))
         b1(0,i,0)=fac*(tdy*cf1(i)*b1(0,i-1,0)+cf2(i)*b1(0,i-2,0))
         b1(i,0,0)=fac*(tdx*cf1(i)*b1(i-1,0,0)+cf2(i)*b1(i-2,0,0))
      END DO

   ! 1 index 0, 1 index 1, other >=1
   
      b1(0,1,1)=fac*(dz*b1(0,1,0)) *3
      b1(1,0,1)=fac*(dz*b1(1,0,0)) *3
      b1(1,1,0)=fac*(dy*b1(1,0,0)) *3

      DO i=2,torderlim-1
         b1(0,1,i)=fac*(dy*b1(0,0,i)+tdz*b1(0,1,i-1)+b1(0,1,i-2))
         b1(1,0,i)=fac*(dx*b1(0,0,i)+tdz*b1(1,0,i-1)+b1(1,0,i-2))
         b1(0,i,1)=fac*(dz*b1(0,i,0)+tdy*b1(0,i-1,1)+b1(0,i-2,1))
         b1(1,i,0)=fac*(dx*b1(0,i,0)+tdy*b1(1,i-1,0)+b1(1,i-2,0))
         b1(i,1,0)=fac*(dy*b1(i,0,0)+tdx*b1(i-1,1,0)+b1(i-2,1,0))
         b1(i,0,1)=fac*(dz*b1(i,0,0)+tdx*b1(i-1,0,1)+b1(i-2,0,1))
      END DO
  
   ! 1 index 0, others >= 2
      
      !DO i=2,torderlim-2   ! sliu
      !   DO j=2,torderlim-i
      
       DO j=2,torderlim-2   
          DO i=2,torderlim-j 
            b1(0,i,j)=fac*(tdy*cf1(i)*b1(0,i-1,j)+tdz*b1(0,i,j-1)+cf2(i)*b1(0,i-2,j)+b1(0,i,j-2))   
            b1(i,0,j)=fac*(tdx*cf1(i)*b1(i-1,0,j)+tdz*b1(i,0,j-1)+cf2(i)*b1(i-2,0,j)+b1(i,0,j-2))
            b1(i,j,0)=fac*(tdx*cf1(i)*b1(i-1,j,0)+tdy*b1(i,j-1,0)+cf2(i)*b1(i-2,j,0)+b1(i,j-2,0))
         END DO
      END DO

 
   ! 2 indices 1, other >= 1
   ! b(1,1,1) is correct, but a little tricky!   
   ! b1(1,1,1)=5.0*dz*fac*b1(1,1,0)

      b1(1,1,1)=fac*(dx*b1(0,1,1)+tdy*b1(1,0,1)+tdz*b1(1,1,0))
      DO i=2,torderlim-2
         b1(1,1,i)=fac*(dx*b1(0,1,i)+tdy*b1(1,0,i)+tdz*b1(1,1,i-1) +b1(1,1,i-2))
         b1(1,i,1)=fac*(dx*b1(0,i,1)+tdy*b1(1,i-1,1)+tdz*b1(1,i,0) +b1(1,i-2,1))
         b1(i,1,1)=fac*(dy*b1(i,0,1)+tdx*b1(i-1,1,1)+tdz*b1(i,1,0) +b1(i-2,1,1))
      END DO

   ! 1 index 1, others >=2
      
   !   DO i=2,torderlim-3
   !      DO j=2,torderlim-i
    
       DO j=2,torderlim-3    ! sliu
          DO i=2,torderlim-j
            b1(1,i,j)=fac*(dx*b1(0,i,j)+tdy*b1(1,i-1,j)+tdz*b1(1,i,j-1) +b1(1,i-2,j)+b1(1,i,j-2))
            b1(i,1,j)=fac*(dy*b1(i,0,j)+tdx*b1(i-1,1,j)+tdz*b1(i,1,j-1) +b1(i-2,1,j)+b1(i,1,j-2))
            b1(i,j,1)=fac*(dz*b1(i,j,0)+tdx*b1(i-1,j,1)+tdy*b1(i,j-1,1) +b1(i-2,j,1)+b1(i,j-2,1))
         END DO
      END DO

   ! all indices >=2
      DO k=2,torderlim-4
         DO j=2,torderlim-2-k
            DO i=2,torderlim-k-j
               b1(i,j,k)=fac*(tdx*cf1(i)*b1(i-1,j,k)+tdy*b1(i,j-1,k)   &
                             +tdz*b1(i,j,k-1)+cf2(i)*b1(i-2,j,k)       &
                             +b1(i,j-2,k)+b1(i,j,k-2)) 
            END DO
         END DO
      END DO  
   
   
   

      !RETURN
   END SUBROUTINE COMP_TCOEFF    
   !======================================================================================
   SUBROUTINE COMVel_DIRECT(tarpos, vel, ibeg, iend, fila_x, fila_y, fila_z, g, r, numfila)
      IMPLICIT NONE

! COMPF_DIRECT directly computes the force on the current target
! particle determined by the global variable TARPOS. 
      REAL(DbKi),DIMENSION(3),        INTENT(IN   ) :: tarpos
      INTEGER(IntKi),                 INTENT(IN   ) :: ibeg,iend,numfila
      REAL(DbKi),DIMENSION(3),        INTENT(  OUT) :: vel
      REAL(DbKi),DIMENSION(numfila,3),INTENT(IN   ) :: fila_x, fila_y, fila_z
      REAL(DbKi),DIMENSION(numfila,1),INTENT(IN   ) :: g,r

! local variables

      INTEGER(IntKi) :: i
      !REAL(DbKi)     :: t1,t2,tx,ty,tz
      
      REAL(DbKi)     :: GMMA, X1, Y1, Z1, X2, Y2, Z2, X2X1, Y2Y1, Z2Z1, L, PXX1, PYY1, PZZ1, PXX2, PYY2, PZZ2
      REAL(DbKi)     :: R1, R2, R1DR2, R1TR2, DEN, UBAR  
      REAL(DbKi)     :: LDR12, CNU, CO, RC 

      !CO = p%FVM%CO ! sliu
      vel = 0.0_DbKi
      

      DO i=ibeg,iend
          X1   =  fila_x(I,1)
          Y1   =  fila_y(I,1)
          Z1   =  fila_z(I,1)
   
          X2   =  fila_x(I,2)
          Y2   =  fila_y(I,2)
          Z2   =  fila_z(I,2)          
          
          GMMA =  g(i,1)
          RC   =  r(i,1)
          
          
          X2X1 =  X2 - X1
          Y2Y1 =  Y2 - Y1
          Z2Z1 =  Z2 - Z1
                             
          L    =  X2X1 * X2X1 + Y2Y1 * Y2Y1 + Z2Z1 * Z2Z1  ! Length of vortex filament (NOTE: L is L^2, as rc is rc^2)
   
          PXX1    =  tarpos(1) - X1
          PYY1    =  tarpos(2) - Y1
          PZZ1    =  tarpos(3) - Z1
          PXX2    =  tarpos(1) - X2
          PYY2    =  tarpos(2) - Y2
          PZZ2    =  tarpos(3) - Z2  
                             
          R1      =  SQRT( PXX1 * PXX1 + PYY1 * PYY1 + PZZ1 * PZZ1 )
          R2      =  SQRT( PXX2 * PXX2 + PYY2 * PYY2 + PZZ2 * PZZ2 )
          R1DR2   =  PXX1 * PXX2 + PYY1 * PYY2 + PZZ1 * PZZ2
          R1TR2   =  R1 * R2          
          
          !IF (p%FVM%ViscFLAG) THEN
             ! Vatistas core model (n=2) .................................................
             LDR12   =  ( X2X1 * PXX1 + Y2Y1 * PYY1 + Z2Z1 * PZZ1) ** 2 
             CNU     =  (( R1 * R1 ) - ( LDR12 / L ))    
             CNU     =  CNU / SQRT(RC ** 2 + CNU  ** 2 )   !sliu: CNU = CNU * (RC ** 2 + CNU  ** 2 ) ** (-1/2) ! gives wrong result
             UBAR    =  CNU * GMMA / (TWOPI * 2) * (R1 + R2) / (R1TR2 * (R1TR2 + R1DR2))
                             
             ! Check infinity and NAN
             IF ( EqualRealNos( L, 0.0_DbKi ) )                             UBAR = 0 
             IF ( EqualRealNos( R1 * R1 - LDR12 / L, 0.0_DbKi ) )           UBAR = 0 
             IF ( EqualRealNos( (R1TR2 * (R1TR2 + R1DR2)) , 0.0_DbKi ) )    UBAR = 0 
             IF ( ISNAN(UBAR) )                                             UBAR = 0  
                                                              
          !ELSE ! Smoothing parameter .................................................
          !   DEN     =  R1TR2 * (R1TR2 + R1DR2) + (p%FVM%DELTA * L)                           
          !   UBAR    =  (GMMA * (R1 + R2)) / (TWOPI * 2) 
          !   UBAR    =  UBAR / DEN 
          !                   
          !   ! Check infinity and NAN
          !   IF ( EqualRealNos( DEN, 0.0 ) )                UBAR = 0 
          !   IF ( ISNAN(UBAR) )                             UBAR = 0  
          !END IF ! (p%FVM%ViscFLAG)
                             
         
          ! Influence is ignored beyond CO                             
          !IF (R1>CO)          UBAR = 0  
          !IF (R2>CO)          UBAR = 0           
          
          
          vel(1) = vel(1) + UBAR * (PYY1 * PZZ2 - PZZ1 * PYY2)
          vel(2) = vel(2) + UBAR * (PZZ1 * PXX2 - PXX1 * PZZ2)
          vel(3) = vel(3) + UBAR * (PXX1 * PYY2 - PYY1 * PXX2)

      END DO   

      !RETURN
      END SUBROUTINE COMVel_DIRECT 
   !======================================================================================
                                      
                         
!RETURN
END SUBROUTINE TREE_COMPVel
!=================================================================================================================================   
SUBROUTINE CLEAN_UP(p)
! CLEANUP deallocates allocated global variables and then
! calls recursive routine REMOVE_NODE to delete the tree.

      IMPLICIT NONE


      TYPE(tnode),POINTER :: p      

! local variables
  
      INTEGER :: err


      CALL REMOVE_NODE(p)
      DEALLOCATE(p, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error (in WInDS): Error deallocating root node! (In Treecode Algorithm)'
         STOP
      END IF 
      NULLIFY(p)
      
CONTAINS
   !==========================================================================================
   RECURSIVE SUBROUTINE REMOVE_NODE(p)
   ! REMOVE_NODE recursively removes each node from the
   ! tree and deallocates its memory for MS array if it
   ! exits.

      IMPLICIT NONE


      TYPE(tnode), POINTER :: p 

   ! local variables

      INTEGER :: i,err

      IF (p%exist_ms .EQ. 1) THEN
         DEALLOCATE(p%ms,STAT=err)
         IF (err .NE. 0) THEN
            WRITE(6,*) 'Error (in WInDS): Error deallocating node MS! (In Treecode Algorithm)'
            STOP
         END IF               
      END IF

      IF (p%num_children .GT. 0) THEN
          DO i=1,p%num_children
            CALL REMOVE_NODE(p%child(i)%p_to_tnode)
            DEALLOCATE(p%child(i)%p_to_tnode,STAT=err)
            IF (err .NE. 0) THEN
               WRITE(6,*) 'Error deallocating node child! '
               STOP
            END IF                           
          END DO
      END IF 

      !RETURN                
      END SUBROUTINE REMOVE_NODE      
      !=============================================================================== 

!RETURN
END SUBROUTINE CLEAN_UP

      
END MODULE WINDS_Treecode
!**********************************************************************************************************************************    