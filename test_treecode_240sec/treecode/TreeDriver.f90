PROGRAM TREEDRIVER
   USE OMP_LIB      
   USE WINDS_Treecode
   USE tree_library 
   USE direct_comp 
   
   IMPLICIT NONE
   
   
   ! arrarys for WInDS...
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: F1_new
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: F2_new
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: Pts_new  
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: GAMMA_new, GAMMA_new_copy 
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: RC_new, RC_new_copy      
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: Uind_new, Uind_new_test  
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: fila_x, fila_x_copy   ! (:,1) Start point; (:,2) End point; (:,3) mid-point;  
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: fila_y, fila_y_copy   ! (:,1) Start point; (:,2) End point; (:,3) mid-point;
      REAL(DbKi), DIMENSION(:,:), ALLOCATABLE :: fila_z, fila_z_copy   ! (:,1) Start point; (:,2) End point; (:,3) mid-point;
      
      
      INTEGER(IntKi)                          :: Len_fila   ! # of filament
      INTEGER(IntKi)                          :: Len_pts    ! # of points   
      INTEGER(IntKi)                          :: i, k, j1, j2, j3 
      CHARACTER(LEN=1024)                     :: input_file
      CHARACTER(LEN=1024)                     :: line  
      
      REAL(DbKi), DIMENSION(5) :: angle
      INTEGER(IntKi), DIMENSION(3) :: order 
      INTEGER(IntKi), DIMENSION(4) :: maxpar 
      REAL(DbKi), DIMENSION(3, 5, 4) :: y, z


      
      REAL(DbKi) :: start, finish, direct_time, tree_time, l2_err, sum1, sum2, dist_tol, delta
      INTEGER(IntKi)  :: cores 
           
! Set parameters

      Len_fila = 129600
      Len_pts =  129546
      input_file = "WInDS_internal_2400_129600_129546_20160304_154329.dat"

      angle   = (/ 0.1 , 0.2 , 0.3,   0.4,   0.5 /)
      order   = (/ 5,     6,     7  /)
      maxpar  = (/ 250, 500,  750,   1000 /)       
      dist_tol= 10000
      delta   = 0.001
      cores   = 8    
      
      call omp_set_num_threads(cores)
      
      
      
      IF (.NOT. ALLOCATED( F1_new ) )        ALLOCATE( F1_new    (Len_fila, 3))
      IF (.NOT. ALLOCATED( F2_new ) )        ALLOCATE( F2_new    (Len_fila, 3))
      IF (.NOT. ALLOCATED( GAMMA_new ) )     ALLOCATE( GAMMA_new (Len_fila, 1))
      IF (.NOT. ALLOCATED( RC_new ) )        ALLOCATE( RC_new    (Len_fila, 1)) 
      IF (.NOT. ALLOCATED( Pts_new ) )       ALLOCATE( Pts_new   (Len_pts, 3))  
      IF (.NOT. ALLOCATED( Uind_new ) )      ALLOCATE( Uind_new  (Len_pts, 3))  
      IF (.NOT. ALLOCATED( Uind_new_test ) ) ALLOCATE( Uind_new_test (Len_pts, 3))  
      
      IF (.NOT. ALLOCATED( fila_x ) ) ALLOCATE( fila_x(Len_fila,3)) 
      IF (.NOT. ALLOCATED( fila_y ) ) ALLOCATE( fila_y(Len_fila,3)) 
      IF (.NOT. ALLOCATED( fila_z ) ) ALLOCATE( fila_z(Len_fila,3)) 
          
          
      
! Read in data
      open(12, file=TRIM(input_file))

      do i=1,6
         read (12,*) line
      enddo      
      
      do i=1,Len_fila
         read (12,*) fila_x(i,1), fila_x(i,2), fila_x(i,3)
      enddo
      
      read (12,*) line
      read (12,*) line
      
      do i=1,Len_fila
         read (12,*) fila_y(i,1), fila_y(i,2), fila_y(i,3)
      enddo     
      
      read (12,*) line
      read (12,*) line
      
      do i=1,Len_fila
         read (12,*) fila_z(i,1), fila_z(i,2), fila_z(i,3)
      enddo    
      
      
      read (12,*) line
      read (12,*) line
      
      do i=1,Len_fila
         read (12,*) GAMMA_new(i,1)
      enddo           
      
      read (12,*) line
      read (12,*) line
      
      do i=1,Len_fila
         read (12,*) RC_new(i,1)
      enddo  
      
      read (12,*) line
      read (12,*) line
      
      do i=1,Len_pts
         read (12,*) Pts_new(i,1), Pts_new(i,2), Pts_new(i,3)
      enddo                 
 
      close(12)
      
! Call direct calculation
      
      do i=1, Len_fila
         F1_new(i, 1)    =     fila_x(i,1)     ! Start point, x
         F1_new(i, 2)    =     fila_y(i,1)     ! Start point, y
         F1_new(i, 3)    =     fila_z(i,1)     ! Start point, z
         F2_new(i, 1)    =     fila_x(i,2)     ! End point, x
         F2_new(i, 2)    =     fila_y(i,2)     ! End point, y
         F2_new(i, 3)    =     fila_z(i,2)     ! End point, z              
      enddo    
      
      call cpu_time(start)
      CALL BiotSavart_direct(F1_new, F2_new, Pts_new, GAMMA_new, RC_new, Uind_new_test, Len_fila, Len_pts, 0.001_DbKi)
      call cpu_time(finish)
      direct_time=finish-start
      
      print *, 'Direct time:', direct_time
      !print *, Uind_new_test(1:10,:)
      
! Call treecode
      
      IF (.NOT. ALLOCATED( fila_x_copy ) )        ALLOCATE( fila_x_copy(Len_fila,3)) 
      IF (.NOT. ALLOCATED( fila_y_copy ) )        ALLOCATE( fila_y_copy(Len_fila,3)) 
      IF (.NOT. ALLOCATED( fila_z_copy ) )        ALLOCATE( fila_z_copy(Len_fila,3))       
      IF (.NOT. ALLOCATED( GAMMA_new_copy ) )     ALLOCATE( GAMMA_new_copy (Len_fila, 1))
      IF (.NOT. ALLOCATED( RC_new_copy ) )        ALLOCATE( RC_new_copy    (Len_fila, 1))      
     

      DO j1=1,5
          DO j2=1,3 
              DO j3=1,4     
                  fila_x_copy = fila_x
                  fila_y_copy = fila_y  
                  fila_z_copy = fila_z
                  GAMMA_new_copy = GAMMA_new
                  RC_new_copy = RC_new
             
                  CALL CPU_TIME(start) ! Start time
                  CALL BiotSavart_Treecode(fila_x_copy, fila_y_copy, fila_z_copy, Pts_new, GAMMA_new_copy, RC_new_copy, Uind_new, angle(j1), order(j2), maxpar(j3), dist_tol, delta, cores)
                  CALL CPU_TIME(finish) ! End time 
                  tree_time=finish-start
      
                  sum1=0.0
                  sum2=0.0
                  ! Loop the points of interest to get L2 norm error. 
                  DO i = 1, Len_pts
                     DO k = 1, 3             
                        sum1 = sum1+ (Uind_new_test(i, k)-Uind_new(i, k))**2
                        sum2 = sum2+ (Uind_new_test(i, k))**2
                     END DO
                  END DO
                  l2_err = SQRT(sum1/sum2)

                  y(j2, j1, j3) = l2_err
                  z(j2, j1, j3) = direct_time/tree_time
                  
                  print *, '........................................'
                  print *, 'angle:', angle(j1), '--- order:', order(j2), '--- leaf:', maxpar(j3)
                  print *, 'Tree time:', tree_time
                  print *, 'Tree l2 error:', l2_err
                  print *, 'Speedup:', direct_time/tree_time


                  ! print *, Uind_new(1:10,:)
      
              enddo 
          enddo 
      enddo 

      print *, '........................................'
      print *, 'angle:', angle
      print *, 'order:', order
      print *, 'maxpar:', maxpar 

      print *, '........................................'

      DO j2=1,3
          print *, y(j2, :, :)
          print *, ' '
      enddo 

      print *, '........................................' 

      DO j2=1,3
          print *, z(j2, :, :)
          print *, ' '
      enddo 

                  
      !write(6,*) 'enter any number to finish: '
      !read(5,*) k

END PROGRAM TREEDRIVER