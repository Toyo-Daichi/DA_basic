! Created on 2020.4.1~
! @author: Toyo_Daichi

program lorenz63

  use kinddef
  use lorenz63_prm

  implicit none

  ! --- setting Parameters
  integer :: nt_asm       ! Period of data assimilation
  integer :: nt_prd       ! Period of prediction
  integer :: obs_interval ! Interval of observation
  
  real(r_size), parameter  :: size_noise_obs = 1.0d0

  character(8)  :: da_method
  character(8)  :: enkf_method
  character(12) :: intg_method
  
  ! --- Physical variable
  real(r_size), allocatable :: x_true(:), y_true(:), z_true(:)
  real(r_size), allocatable :: x_sim(:), y_sim(:), z_sim(:)
  real(r_size), allocatable :: x_anl(:), y_anl(:), z_anl(:)
  real(r_size), allocatable :: x_obs(:), y_obs(:), z_obs(:)

  real(r_size), allocatable :: x_anl_m(:,:), y_anl_m(:,:), z_anl_m(:,:)
  real(r_size), allocatable :: x_prtb(:), y_prtb(:), z_prtb(:)
  
  real(r_size), allocatable :: Ef(:,:), Ea(:,:) 
  real(r_size), allocatable :: obs_prtbmtx(:,:)
  real(r_size), allocatable :: x_innov(:), y_innov(:), z_innov(:)

  ! --- Temporal state vector
  real(r_size), allocatable :: xt_vec(:,:)
  real(r_size), allocatable :: yt_vec(:,:)
  
  ! --- Matrix(element 1:x, 2:y, 3:z)
  ! +++ default setting
  ! Pf = ( Pxx: 1.0  Pxy: 0.0 Pxz: 0.0
  !        Pyx: 0.0  Pyy: 1.0 Pyz: 0.0 
  !        Pzx: 0.0  Pzy: 0.0 Pzz: 1.0 )
  
  real(r_size), allocatable :: JM(:,:)  ! state transient matrix (=Jacobean　matrix)
  real(r_size), allocatable :: Pf(:,:)  ! Forecast error convariance matrix (in KF, EnKF)
  real(r_size), allocatable :: Pa(:,:)  ! Analysis error convariance matrix
  real(r_size), allocatable ::  I(:,:)  ! eye_matrix
  real(r_size), allocatable ::  R(:,:)  ! Observation error convariance matrix
  real(r_size), allocatable :: Kg(:,:)  ! Kalman gain
  real(r_size), allocatable :: Kh(:,:)  ! Kalman gain for pertuvation
  real(r_size), allocatable ::  H(:,:)  ! Observation operator
  real(r_size), allocatable :: Pr(:,:)

  real(r_size), allocatable :: obs_mtx(:,:)
  real(r_size), allocatable :: obs_inv(:,:)
  
  ! --- for EnSRF
  real(r_size), allocatable :: obs_srf_mtx_v1(:,:), obs_srf_mtx_v2(:,:)
  real(r_size), allocatable :: obs_srf_inv_v1(:,:), obs_srf_inv_v2(:,:)

  real(r_size), allocatable :: hx(:,:)
  real(r_size), allocatable :: hdxf(:,:)
  real(r_size), allocatable :: work(:)
  
  ! --- Output control
  logical                   :: opt_veach = .false.
  real(r_size), allocatable :: yt_vec4out(:, :) ! for output
  integer, parameter        :: output_interval = 1
  character(256)            :: output_file
  character(256)            :: output_file_error_covariance
  character(1096)           :: linebuf
  
  ! --- Working variable
  !real(r_size) :: x_trend, y_trend, z_trend ! for tl_model
  integer :: it, il
  integer :: mems, imem
  integer :: ierr
  integer :: lda
  real(r_size) :: dt_obs
  real(r_dble) :: Gnoise ! Gaussian noise
  real(r_size) :: alpha  ! inflation
  real(r_size), parameter:: undef = -999.e0
  
  ! --- for matrix calculation consistent
  integer      :: ipiv, lwork
  
  !======================================================================
  ! Data assimilation
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------
  real(r_size) :: x_tinit, y_tinit, z_tinit
  real(r_size) :: x_sinit, y_sinit, z_sinit
  
  namelist /set_dim/ nx, ny
  namelist /set_parm/ nt_asm, nt_prd, obs_interval
  namelist /da_setting/ da_method, alpha
  namelist /intg_setting/ intg_method
  namelist /enkf_setting/ mems, enkf_method
  namelist /initial_score/ x_tinit, y_tinit, z_tinit, x_sinit, y_sinit, z_sinit
  namelist /output/ output_file, output_file_error_covariance, opt_veach
  
  read(5, nml=set_dim, iostat = ierr)
  read(5, nml=set_parm, iostat = ierr)
  read(5, nml=da_setting, iostat = ierr)
  read(5, nml=intg_setting, iostat = ierr)
  read(5, nml=enkf_setting, iostat = ierr)
  read(5, nml=initial_score, iostat = ierr)
  read(5, nml=output, iostat = ierr)
  
  ! +++ name list state num(io) check
  if (ierr < 0 ) then
    write(6,*) '   Msg : Main[ .sh /  @namelist ] '
    write(6,*) '   Not found namelist.        '
    write(6,*) '   Use default values.        '
  else if (ierr > 0) then
    write(6,*) '   Msg : Main[ .sh /  @namelist ] '
    write(6,*) '   *** Warning : Not appropriate names in namelist !! Check !!'
    write(6,*) '   Stop : lorenz63_main.f90              '
    stop
  end if
  
  ! +++ display namelist
  write(6,*) 'Data assimlation method :: ', da_method
  write(6,*) 'Integral method         :: ', intg_method
  
  allocate(x_true(0:nt_asm+nt_prd), y_true(0:nt_asm+nt_prd), z_true(0:nt_asm+nt_prd))
  allocate(x_anl(0:nt_asm+nt_prd), y_anl(0:nt_asm+nt_prd), z_anl(0:nt_asm+nt_prd))
  allocate(x_sim(0:nt_asm+nt_prd), y_sim(0:nt_asm+nt_prd), z_sim(0:nt_asm+nt_prd))
  allocate(x_anl_m(0:nt_asm, mems), y_anl_m(0:nt_asm, mems), z_anl_m(0:nt_asm, mems))
  allocate(x_prtb(mems), y_prtb(mems), z_prtb(mems))
  allocate(x_obs(0:nt_asm/obs_interval), y_obs(0:nt_asm/obs_interval), z_obs(0:nt_asm/obs_interval))
  
  allocate(xt_vec(nx,1), yt_vec(ny,1))
  allocate(Ef(nx,mems), Ea(nx,mems), obs_prtbmtx(ny, mems))
  
  allocate(x_innov(mems), y_innov(mems), z_innov(mems))
  allocate(yt_vec4out(ny, 0:nt_asm+nt_prd))
  
  ! covariance matrix set.
  allocate(Pf(1:nx, 1:nx))
  allocate(Pa(1:nx, 1:nx))
  allocate(JM(1:nx, 1:nx))
  allocate( I(1:nx, 1:nx))
  allocate(Kg(1:nx, 1:ny))
  allocate(Kh(1:nx, 1:ny))
  allocate( H(1:ny, 1:nx))
  allocate( R(1:ny, 1:ny))
  allocate(Pr(1:ny, 1:ny))

  allocate(obs_mtx(1:ny, 1:ny), obs_inv(1:ny, 1:ny))
  allocate(obs_srf_mtx_v1(1:ny, 1:ny), obs_srf_mtx_v2(1:ny, 1:ny))
  allocate(obs_srf_inv_v1(1:ny, 1:ny), obs_srf_inv_v2(1:ny, 1:ny))
  
  allocate(hx(ny, 1))
  allocate(hdxf(ny, 1))
  
  lda = ny; lwork = ny
  allocate(work(1:ny))
  
  !----------------------------------------------------------------------
  ! +++ initial setting
  !----------------------------------------------------------------------
  
  x_true(0) = x_tinit; y_true(0) = y_tinit; z_true(0) = z_tinit
  x_sim(0)  = x_sinit; y_sim(0)  = y_sinit; z_sim(0)  = z_sinit
  
  Pf = 0.d0; Pa = 0.d0; I = 0.d0
  H  = 0.d0; R  = 0.d0
  Kg = 0.d0; Kh = 0.0d0
  
  forall ( il=1:nx ) Pf(il, il) = 1.0d0
  forall ( il=1:nx ) Pa(il, il) = 1.0d0
  forall ( il=1:ny ) R(il, il)  = size_noise_obs
  forall ( il=1:ny ) H(il, il)  = 1.0d0
  forall ( il=1:nx ) I(il, il)  = 1.0d0
  
  !----------------------------------------------------------------------
  ! +++ if 3D var method (Not tested ... )
  !Pf(1,1) = 0.01; Pf(1,2) = 0.01; Pf(1,3) = 0.01
  !Pf(2,1) = 0.01; Pf(2,2) = 0.01; Pf(2,3) = 0.01
  !Pf(3,1) = 0.01; Pf(3,2) = 0.01; Pf(3,3) = 0.01
  !----------------------------------------------------------------------
  
  ! --- Initialization of random number generator
  call random_seed()
  
  ! --- Sec2. True field and observations
  do it = 1, nt_asm+nt_prd
    ! forward time step
    
    call cal_Lorenz(                           &
    x_true(it-1), y_true(it-1), z_true(it-1),  & ! IN
    x_k(1), y_k(1), z_k(1)                     & ! OUT
    )
    
    !------------------------------------------------------- 
    ! +++ Euler method
    !------------------------------------------------------- 
    if ( trim(intg_method) == 'Euler' ) then
      x_true(it) = x_true(it-1) + dt*x_k(1)
      y_true(it) = y_true(it-1) + dt*y_k(1)
      z_true(it) = z_true(it-1) + dt*z_k(1)
      
      !------------------------------------------------------- 
      ! +++ Runge-Kutta method
      !------------------------------------------------------- 
    else if ( trim(intg_method) == 'Runge-Kutta' ) then 
      
      call Lorenz63_Runge_Kutta(                  &
      x_true(it-1), y_true(it-1), z_true(it-1), & ! IN
      x_true(it), y_true(it), z_true(it)        & ! OUT
      )
      
    end if
    
    ! making observations
    if ((mod(it, obs_interval) == 0) .and. (it <= nt_asm)) then
      
      call gaussian_noise(sqrt(R(1,1)), Gnoise)
      ! Generate observation by adding Gaussian noise to true value
      x_obs(it/obs_interval) = x_true(it) + Gnoise
      
      call gaussian_noise(sqrt(R(2,2)), Gnoise)
      y_obs(it/obs_interval) = y_true(it) + Gnoise
      
      call gaussian_noise(sqrt(R(3,3)), Gnoise)
      z_obs(it/obs_interval) = z_true(it) + Gnoise
      
      ! +++ for debug
      !write(6,*) 'time_step, x_obs, y_obs, z_obs', &
      !it, x_obs(it/obs_interval), y_obs(it/obs_interval), z_obs(it/obs_interval)
    end if
  end do
  
  ! --- Sec3. Simulation run without DA
  do it = 1, nt_asm+nt_prd
    ! forward time step
    
    call cal_Lorenz(                           &
    x_sim(it-1), y_sim(it-1), z_sim(it-1),     & ! IN
    x_k(1), y_k(1), z_k(1)                     & ! OUT
    )
    
    !------------------------------------------------------- 
    ! +++ Euler method
    !------------------------------------------------------- 
    if ( trim(intg_method) == 'Euler' ) then
      x_sim(it) = x_sim(it-1) + dt * x_k(1)
      y_sim(it) = y_sim(it-1) + dt * y_k(1)
      z_sim(it) = z_sim(it-1) + dt * z_k(1)
      
      !------------------------------------------------------- 
      ! +++ Runge-Kutta method
      !------------------------------------------------------- 
    else if ( trim(intg_method) == 'Runge-Kutta' ) then 
      
      call Lorenz63_Runge_Kutta(               &
      x_sim(it-1), y_sim(it-1), z_sim(it-1), & ! IN
      x_sim(it), y_sim(it), z_sim(it)        & ! OUT
      )
      
    end if
  end do
  
  !====================================================================== 
  ! --- Sec4. Data assimilation

  if ( da_method == 'KF' ) then
    x_anl(0) = x_sim(0)
    y_anl(0) = y_sim(0)
    z_anl(0) = z_sim(0)
    
    do it = 1, nt_asm
      if ( opt_veach ) call output_Pf(it, nt_asm, Pa)
      
      ! 4.1: Time integration
      if ( trim(intg_method) == 'Euler' ) then
        !------------------------------------------------------- 
        ! +++ Euler method
        !------------------------------------------------------- 
        
        call cal_Lorenz(                         &
        x_anl(it-1), y_anl(it-1), z_anl(it-1),   & ! IN
        x_k(1), y_k(1), z_k(1)                   & ! OUT
        )
        
        x_anl(it) = x_anl(it-1) + dt * x_k(1)
        y_anl(it) = y_anl(it-1) + dt * y_k(1)
        z_anl(it) = z_anl(it-1) + dt * z_k(1)
          
      else if ( trim(intg_method) == 'Runge-Kutta' ) then 
          !------------------------------------------------------- 
          ! +++ Runge-Kutta method
          !------------------------------------------------------- 
          
        call Lorenz63_Runge_Kutta(                &
        x_anl(it-1), y_anl(it-1), z_anl(it-1),    & ! IN
        x_anl(it), y_anl(it), z_anl(it)           & ! OUT
        )
      end if
        
      if (mod(it, obs_interval) == 0) then
        write(6,*) 'Data assim. time == ', it*dt, 'sec.'
        write(6,*) ''
        write(6,*) '  TRUTH   = ', x_true(it), y_true(it), z_true(it)
        write(6,*) '  PREDICT = ', x_sim(it), y_sim(it), z_sim(it)
        write(6,*) '  OBSERVE = ', x_obs(it/obs_interval), y_obs(it/obs_interval), z_obs(it/obs_interval)
        write(6,*) ''

        !--------------------------------------------------------------------------------
        ! 4.2: Kalman fileter
        ! +++ 4.2.1 State Transient Matrix
        ! >> solve from TL function
        !call TL_Lorez63_Runge_Kutta(             &
        !  x_anl(it-1), y_anl(it-1), z_sim(it-1), & ! IN:  state
        !  Pa(1,1), Pa(1,2), Pa(1,3),             & ! IN:  dX 1st. step
        !  Pa(2,1), Pa(2,2), Pa(2,3),             & ! IN:  dY 2nd. step
        !  Pa(3,1), Pa(3,2), Pa(3,3),             & ! IN:  dZ 3rd. step
        !  x_trend, y_trend, z_trend              & ! OUT: trend
        !)
        ! +++ for check
        ! write(6,*) x_trend, y_trend, z_trend
        
        ! >> Jacobian　matrix (*** origin form)
        !JM(1,1) = -sig              ;JM(1,2) = sig             ;JM(1,3) = 0.0d0
        !JM(2,1) = gamm -z_anl(it-1) ;JM(2,2) = -1.d0           ;JM(2,3) = -x_anl(it-1)
        !JM(3,1) = y_anl(it-1)       ;JM(3,2) = x_anl(it-1)     ;JM(3,3) = -b
        !--------------------------------------------------------------------------------
        
        ! >> Jacobian　matrix with time evolution
        dt_obs = dt*obs_interval
        JM(1,1) = 1.0d0 -dt_obs*sig          ;JM(1,2) = dt_obs*sig          ;JM(1,3) = 0.0d0
        JM(2,1) = dt_obs*(gamm -z_anl(it-1)) ;JM(2,2) = 1.d0 - dt_obs       ;JM(2,3) = -dt_obs*x_anl(it-1)
        JM(3,1) = dt_obs*y_anl(it-1)         ;JM(3,2) = dt_obs*x_anl(it-1) ;JM(3,3) = 1.0d0 -dt_obs*b
        
        write(6,*) '  ANALYSIS ERROR COVARIANCE before one step.'
        call confirm_matrix(Pa, nx, nx)
        
        !------------------------------------------------------- 
        ! +++ making Pf
        Pf = matmul(JM, matmul(Pa, transpose(JM)))
        Pf = Pf*(1.0d0 + alpha)
        !------------------------------------------------------- 
        
        write(6,*) '  PREDICTION ERROR COVARIANCE on present step'
        call confirm_matrix(Pf, nx, nx)
        write(6,*) ''
        
        !------------------------------------------------------- 
        ! >> 4.2.3 Kalman gain: Weighting of model result and obs.
        ! +++ making inverse matrix
        !------------------------------------------------------- 
        
        obs_mtx = matmul(H, matmul(Pf, transpose(H))) + R
        obs_inv = obs_mtx
        
        call dgetrf(ny, ny, obs_inv, lda, ipiv, ierr)
        call dgetri(ny, obs_inv, lda, ipiv, work, lwork, ierr)
        ! It will be changed for some reason, so support.
        intg_method = 'Runge-Kutta'
        
        !------------------------------------------------------- 
        ! +++ inverse matrix calculate for 2x2 on formula
        ! >> call inverse_matrix_for2x2(obs_mtx, obs_inv)
        !------------------------------------------------------- 
        
        Kg = matmul(matmul(Pf, transpose(H)), obs_inv)
        
        !-------------------------------------------------------
        ! >> 4.2.4 calculate innovation and correlation
        ! +++ Kalman Filter Main equation
        !-------------------------------------------------------
        xt_vec(1,1) = x_anl(it)
        xt_vec(2,1) = y_anl(it)
        xt_vec(3,1) = z_anl(it)
        
        yt_vec(1,1) = x_obs(it/obs_interval)
        yt_vec(2,1) = y_obs(it/obs_interval)
        yt_vec(3,1) = z_obs(it/obs_interval)
        
        hx = matmul(H, xt_vec)
        hdxf = yt_vec - hx
        
        xt_vec = xt_vec + matmul(Kg, hdxf)
        x_anl(it) = xt_vec(1,1); y_anl(it) = xt_vec(2,1); z_anl(it) = xt_vec(3,1)
        
        ! >> 4.2.5 analysis error covariance matrix
        Pa = matmul(I - matmul(Kg, H), Pf)
      end if
    end do

  else if ( da_method == 'EnKF' ) then
    ! +++ initial setting
    do imem = 1, mems
      call gaussian_noise(sqrt(Pa(1,1)), Gnoise)
      x_anl_m(0, imem) = x_sim(0) + Gnoise
      call gaussian_noise(sqrt(Pa(2,2)), Gnoise)
      y_anl_m(0, imem) = y_sim(0) + Gnoise
      call gaussian_noise(sqrt(Pa(3,3)), Gnoise)
      z_anl_m(0, imem) = z_sim(0) + Gnoise
    end do
      
    x_anl(0) = sum(x_anl_m(0, 1:mems))/mems 
    y_anl(0) = sum(y_anl_m(0, 1:mems))/mems 
    z_anl(0) = sum(z_anl_m(0, 1:mems))/mems
    
    do it = 1, nt_asm
      if ( opt_veach ) call output_Pf(it, nt_asm, Pf)
        
      ! 4.1: Time integration
      do imem = 1, mems
        ! 4.1: Time integration
        call cal_Lorenz(                                                   &
          x_anl_m(it-1, imem), y_anl_m(it-1, imem), z_anl_m(it-1, imem),   & ! IN
          x_k(1), y_k(1), z_k(1)                                           & ! OUT
        )
        
        !------------------------------------------------------- 
        ! +++ Euler method
        if ( trim(intg_method) == 'Euler' ) then
          x_anl_m(it, imem) = x_anl_m(it-1, imem) + dt * x_k(1)
          y_anl_m(it, imem) = y_anl_m(it-1, imem) + dt * y_k(1)
          z_anl_m(it, imem) = z_anl_m(it-1, imem) + dt * z_k(1)
          
          !------------------------------------------------------- 
          ! +++ Runge-Kutta method
        else if ( trim(intg_method) == 'Runge-Kutta' ) then
          
          call Lorenz63_Runge_Kutta(                                          &
            x_anl_m(it-1, imem), y_anl_m(it-1, imem), z_anl_m(it-1, imem),    & ! IN
            x_anl_m(it, imem), y_anl_m(it, imem), z_anl_m(it, imem)           & ! OUT
          )
        end if
      end do
        
      if(mod(it, obs_interval) == 0) then
        write(6,*) 'Data assim. time == ', it*dt, 'sec.'
        write(6,*) ''
        write(6,*) '  TRUTH   = ', x_true(it), y_true(it), z_true(it)
        write(6,*) '  PREDICT = ', x_sim(it), y_sim(it), z_sim(it)
        write(6,*) '  OBSERVE = ', x_obs(it/obs_interval), y_obs(it/obs_interval), z_obs(it/obs_interval)
        write(6,*) ''
        x_anl(it) = sum(x_anl_m(it, 1:mems))/mems
        y_anl(it) = sum(y_anl_m(it, 1:mems))/mems
        z_anl(it) = sum(z_anl_m(it, 1:mems))/mems
        Pf = 0.0d0
        do imem = 1, mems
          x_prtb(imem) = x_anl_m(it, imem) - x_anl(it)
          y_prtb(imem) = y_anl_m(it, imem) - y_anl(it)
          z_prtb(imem) = z_anl_m(it, imem) - z_anl(it)
          
          !------------------------------------------------------- 
          ! making Pf (one shot)
          ! +++ Dispersion
          Pf(1,1) = Pf(1,1) + x_prtb(imem)**2/(mems-1)
          Pf(2,2) = Pf(2,2) + y_prtb(imem)**2/(mems-1)
          Pf(3,3) = Pf(3,3) + z_prtb(imem)**2/(mems-1)
          
          ! +++　Covariance(x,y; x,z; y,z)
          Pf(1,2) = Pf(1,2) + x_prtb(imem)*y_prtb(imem)/(mems-1)
          Pf(2,1) = Pf(1,2)
          Pf(1,3) = Pf(1,3) + x_prtb(imem)*z_prtb(imem)/(mems-1)
          Pf(3,1) = Pf(1,3)
          Pf(2,3) = Pf(2,3) + y_prtb(imem)*z_prtb(imem)/(mems-1)
          Pf(3,2) = Pf(1,3)
          !------------------------------------------------------- 
        end do
       Pf = Pf*(1.0d0 + alpha)
        
        write(6,*) ' PREDICTION ERROR COVARIANCE on present step'
        call confirm_matrix(Pf, nx, nx)
        write(6,*) ''
        
        obs_mtx = matmul(H, matmul(Pf, transpose(H))) + R
        obs_inv = obs_mtx
        
        call dgetrf(ny, ny, obs_inv, lda, ipiv, ierr)
        call dgetri(ny, obs_inv, lda, ipiv, work, lwork, ierr)
        intg_method = 'Runge-Kutta'
        
        write(6,*) ' KALMAN GAIN WEIGHTING MATRIX '
        Kg = matmul(matmul(Pf, transpose(H)), obs_inv)
          
        call confirm_matrix(Kg, nx, ny)
        write(6,*) ''
        
        if ( trim(enkf_method) == 'PO' ) then
          !------------------------------------------------------- 
          ! >> EnKF
          ! +++ Pertuturbed observation method (PO)
          ! (Phase. plus OBS pertubation)
          !------------------------------------------------------- 
          do imem = 1, mems
            call gaussian_noise(sqrt(R(1,1)), Gnoise)
            x_innov(imem) = x_obs(it/obs_interval) + Gnoise
            obs_prtbmtx(1,imem) = x_innov(imem)
            call gaussian_noise(sqrt(R(2,2)), Gnoise)
            y_innov(imem) = y_obs(it/obs_interval) + Gnoise
            obs_prtbmtx(2,imem) = y_innov(imem)
            call gaussian_noise(sqrt(R(3,3)), Gnoise)
            z_innov(imem) = z_obs(it/obs_interval) + Gnoise
            obs_prtbmtx(3,imem) = z_innov(imem)
            
            xt_vec(1,1) = x_anl_m(it, imem)
            xt_vec(2,1) = y_anl_m(it, imem)
            xt_vec(3,1) = z_anl_m(it, imem)
            
            yt_vec(1,1) = obs_prtbmtx(1,imem)
            yt_vec(2,1) = obs_prtbmtx(2,imem)
            yt_vec(3,1) = obs_prtbmtx(3,imem)
            
            hx = matmul(H, xt_vec)
            hdxf = yt_vec - hx
            
            xt_vec = xt_vec + matmul(Kg, hdxf)
            x_anl_m(it, imem) = xt_vec(1,1); y_anl_m(it, imem) = xt_vec(2,1); z_anl_m(it, imem) = xt_vec(3,1)
          end do
          
          x_anl(it) = sum(x_anl_m(it, 1:mems))/mems
          y_anl(it) = sum(y_anl_m(it, 1:mems))/mems
          z_anl(it) = sum(z_anl_m(it, 1:mems))/mems
          
        else if ( trim(enkf_method) == 'SRF' ) then
          !----------------------------------------------------------- 
          ! >> EnKF (-> Serial EnSRF)
          ! +++ Squared root fileter method (SRF)
          ! Ea = (I - K^H)*Ef
          ! K^ = 
          ! Pf H^T [(H Pf H^T + R)^-1/2]^T [(H Pf H^T + R)^-1/2]^-1
          !----------------------------------------------------------- 
          
          xt_vec(1,1) = x_anl(it); xt_vec(2,1) = y_anl(it); xt_vec(3,1) = z_anl(it)
          
          Ef(1, :) = x_prtb(:)
          Ef(2, :) = y_prtb(:)
          Ef(3, :) = z_prtb(:)
          
          obs_srf_mtx_v1 = matmul(H, matmul(Pf, transpose(H))) + R
          call f01epf(ny, obs_srf_mtx_v1, lda, ierr)
          obs_srf_inv_v1 = obs_srf_mtx_v1

          write(6,*) ' OBS INVERSE MATRIX '
          call confirm_matrix(obs_srf_mtx_v1, ny, ny)
          write(6,*) ''
          call confirm_matrix(obs_srf_inv_v1, ny, ny)
          write(6,*) ''
          
          call dgetrf(ny, ny, obs_srf_inv_v1, lda, ipiv, ierr)
          call dgetri(ny, obs_srf_inv_v1, lda, ipiv, work, lwork, ierr)
          
          obs_srf_mtx_v2 = sqrt(obs_srf_mtx_v1) + sqrt(R)
          obs_srf_inv_v2 = obs_srf_mtx_v2
          
          call dgetrf(ny, ny, obs_srf_inv_v2, lda, ipiv, ierr)
          call dgetri(ny, obs_srf_inv_v2, lda, ipiv, work, lwork, ierr)
          
          Kh = matmul(matmul(Pf, transpose(H)), transpose(obs_srf_inv_v1))
          Kh = matmul(Kh, obs_srf_inv_v2)
          
          write(6,*) ' KALMAN GAIN HAT WEIGHTING MATRIX '
          call confirm_matrix(Kh, nx, ny)
          write(6,*) ''
          Ea = matmul(I - matmul(Kh, H), Ef)

          xt_vec(1,1) = xt_vec(1,1) + sum(Ea(1,:))/mems
          xt_vec(2,1) = xt_vec(2,1) + sum(Ea(2,:))/mems
          xt_vec(3,1) = xt_vec(3,1) + sum(Ea(3,:))/mems

          x_anl(it) = xt_vec(1,1)
          y_anl(it) = xt_vec(2,1)
          z_anl(it) = xt_vec(3,1)

        end if
      end if
    end do
  end if
  
  ! --- Sec5. Prediction after Data assimilation
    do it = nt_asm+1, nt_asm+nt_prd
      write(6,*) 'Data assim. time step: ', it
      ! forward time step
      
      call cal_Lorenz(                           &
      x_anl(it-1), y_anl(it-1), z_anl(it-1),     & ! IN
      x_k(1), y_k(1), z_k(1)                     & ! OUT
      )
      
      !------------------------------------------------------- 
      ! +++ Euler method
      if ( trim(intg_method) == 'Euler' ) then
        x_anl(it) = x_anl(it-1) + dt * x_k(1)
        y_anl(it) = y_anl(it-1) + dt * y_k(1)
        z_anl(it) = z_anl(it-1) + dt * z_k(1)
        
        !------------------------------------------------------- 
        ! +++ Runge-Kutta method
      else if ( trim(intg_method) == 'Runge-Kutta' ) then 
        
        call Lorenz63_Runge_Kutta(                &
        x_anl(it-1), y_anl(it-1), z_anl(it-1),    & ! IN
        x_anl(it), y_anl(it), z_anl(it)           & ! OUT
        )
        
      end if
    end do
    
    ! --- Sec6. Writing OUTPUT
    if ( opt_veach ) then
      yt_vec4out(:, 0:nt_asm+nt_prd) = undef
      do it = 1, nt_asm
        if (mod(it, obs_interval) == 0) then
          yt_vec4out(1, it) = x_obs(it/obs_interval)
          yt_vec4out(2, it) = y_obs(it/obs_interval)
          yt_vec4out(3, it) = z_obs(it/obs_interval)
        end if
      end do
      
      open (1, file=trim(output_file), status='replace')
      write(1,*) 'timestep, x_true, y_true, z_true, x_sim, y_sim, z_sim, x_anl, y_anl, z_anl, x_obs, y_obs, z_obs'
      do it = 0, nt_asm+nt_prd
        if (mod(it, output_interval) == 0) then
          write(linebuf, '(f5.2, ",", 9(f12.7, ","), 2(F7.2, ","), F7.2)')  & 
            dt*it,                                      &
            x_true(it), y_true(it), z_true(it),         &
            x_sim(it), y_sim(it), z_sim(it),            &
            x_anl(it), y_anl(it), z_anl(it),            &
            yt_vec4out(1, it), yt_vec4out(2, it), yt_vec4out(3,it)
          call del_spaces(linebuf)
          write(1, '(a)') trim(linebuf)
        end if
      end do
      close(1)
      write(6,*) '-------------------------------------------------------'
      write(6,*) '+++ Check Writing output system,  '
      write(6,*) ' && Successfuly output !!!        '  
   
    else if ( .not. opt_veach ) then
      write(6,*) '-------------------------------------------------------'
      write(6,*) '+++ Check calculation system,  '
      write(6,*) ' && Successfuly calculate !!!  '  

    end if

contains

  subroutine inverse_matrix_for2x2(         &
    matrix,                                 & ! IN:  input matrix
    inv_matrix                              & ! OUT: inverse matrix
    )
    
    implicit none
    real(r_size)               :: inv_prm
    real(r_size), intent(in)   :: matrix(2,2) 
    real(r_size), intent(out)  :: inv_matrix(2,2)

    write(6,*) ' +++ calculate on 2x2 inverse formula'
    inv_prm = 1.0d0 / ( matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1) )
    inv_matrix(1,1) =  inv_prm*matrix(2,2); inv_matrix(1,2) = -inv_prm*matrix(1,2)
    inv_matrix(2,1) = -inv_prm*matrix(2,1); inv_matrix(2,2) =  inv_prm*matrix(1,1)

    return
  end subroutine

  subroutine output_Pf(it, last_step, error_covariance_matrix)
    implicit none
    integer, intent(in)       :: it, last_step
    real(r_size), intent(in)  :: error_covariance_matrix(3,3)

    if ( it == 1 ) then
      open(2, file=trim(output_file_error_covariance), status='replace')
      write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        ! +++ for debug
        !write(6,*) '+++ err covariance matrix 1st. step'
        !write(6,*) error_covariance_matrix(:,:)
        
      else if ( it /= 1 .and. it /= last_step) then
        write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        ! +++ for debug
        !write(6,*) '+++ err covariance matrix 2nd. step ~'
        !write(6,*) error_covariance_matrix(:,:)
        
      else if ( it == last_step ) then
        write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        ! +++ for debug
        !write(6,*) '+++ err covariance matrix last step '
        !write(6,*) error_covariance_matrix(:,:)
      close(2)
    end if
      
  end subroutine

  subroutine gaussian_noise(    &
    size_gnoise,                & ! IN: One case is observation error diag.
    gnoise                      & ! OUT
    )
    ! Generate Gaussian Noise (Gnoise) from uniform random number
    ! based on Box-Muller method

    implicit none
    real(kind=r_dble),intent(in)  :: size_gnoise
    real(kind=r_dble),intent(out) :: gnoise
  
    ! constant
    real(kind=r_dble),parameter :: pi=3.14159265358979
  
    ! working variable
    real(kind=r_dble) :: noise1,noise2
  
    call random_number(noise1)
    call random_number(noise2)
    gnoise=size_gnoise*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
  
  end subroutine gaussian_noise
  
  subroutine del_spaces(space)
    implicit none
    character(*), intent(inout) :: space
    character(len=len(space))   :: tmp
    integer ::  i, j

    j = 1
    do i = 1, len(space)
      if (space(i:i)==' ') cycle
      tmp(j:j) = space(i:i)
      j = j + 1
    end do
    space = tmp(1:j-1)
  end subroutine del_spaces

  subroutine confirm_matrix(X,N,M)
    implicit none

    integer, intent(in)          :: N,M
    double precision, intent(in) :: X(N,M)
    integer                      :: i, j
    
    do i=1,n
      do j=1,M
        write(*,fmt='(f15.8)',advance='no') X(i,j)
      end do
      write(*,*)
    end do
    !print *, "==============================="
  end subroutine


end program lorenz63
