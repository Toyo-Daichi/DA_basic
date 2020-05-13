! Created on 2020.4.1
! @author: Toyo_Daichi

program lorenz63

  use kinddef
  use lorenz63_prm

  implicit none

  ! --- setting Parameters
  integer :: nt_asm       ! Period of data assimilation
  integer :: nt_prd       ! Period of prediction
  integer :: obs_interval ! Interval of observation

  real(r_size), parameter  :: dt = 1.0d-2 ! Time step

  character(8)  :: da_method
  character(12) :: intg_method
  
  ! --- Physical variable
  real(r_size), allocatable :: x_true(:), y_true(:), z_true(:)
  real(r_size), allocatable :: x_sim(:), y_sim(:), z_sim(:)
  
  real(r_size), allocatable :: x_da(:), y_da(:), z_da(:)
  real(r_size), allocatable :: x_obs(:), y_obs(:)

  real(r_size), allocatable :: x_da_m(:, :), y_da_m(:, :), z_da_m(:, :)
  real(r_size), allocatable :: x_prtb(:), y_prtb(:), z_prtb(:)
  
  real(r_size), allocatable :: obs_ens(:,:)
  real(r_size), allocatable :: x_innov(:), y_innov(:)

  ! --- Matrix(element 1:x, 2:y, 3:z)
  ! +++ default setting
  ! Pf = (  Pxx: 1.0  Pxy: 0.0 Pxz: 0.0
  !         Pyx: 0.0  Pyy: 1.0 Pyz: 0.0 
  !         Pzx: 0.0  Pzy: 0.0 Pzz: 1.0 )

  real(r_size) :: x_a(nx,1)
  real(r_size) :: yt_obs(ny,1)
  
  real(r_size) :: JM(nx,nx)  ! state transient matrix
  real(r_size) :: Pf(nx,nx)  ! Forecast error convariance matrix (in KF, EnKF)
  real(r_size) :: Pa(nx,nx)  ! Analysis error convariance matrix
  real(r_size) ::  I(nx,nx)  ! eye_matrix
  real(r_size) ::  R(ny,ny)  ! Observation error convariance matrix
  real(r_size) :: Kg(nx,ny)  ! Kalman gain
  real(r_size) ::  H(ny,nx)  ! Observation operator
  
  ! --- Output control
  real(r_size), allocatable :: obs_chr(:, :)
  integer, parameter        :: output_interval = 1
  logical                   :: opt_veach = .true.
  character(256)            :: output_file
  character(256)            :: output_file_error_covariance
  character(1096)           :: linebuf

  ! --- Working variable
  integer :: it, il
  integer :: mems, imem
  integer :: ierr
  integer :: iflag
  integer :: iter
  real(r_dble) :: Gnoise ! Gaussian noise
  real(r_dble) :: delta ! Gaussian noise
  real(r_size) :: alpha  ! inflation
  real(r_size), parameter:: undef = -999.e0

  ! --- matrix calculation
  ! for inverse
  integer      :: ipiv(1:nx), lwork
  real(r_size) :: lwork0
  real(r_size) :: hx(ny,1)
  real(r_size) :: hdxf(ny,1)
  
  real(r_size), allocatable :: work_on(:)
  
  ! --- Inverse matrix formula for 2x2
  real(r_size) :: inv_matrix(ny,ny)
  real(r_size) :: R_HPHt(ny,ny), HPHt(ny,ny)
  real(r_size) :: eye_matrix(ny,ny)

  !======================================================================
  ! Data assimilation
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------
  real(r_size) :: x_tinit, y_tinit, z_tinit
  real(r_size) :: x_sinit, y_sinit, z_sinit
  real(r_size) :: x_e, y_e, z_e
  real(r_size) :: Pf_init(9), B_init(9)
  real(r_size) :: R_init(4)
  real(r_size) :: Kg_init(6)
  real(r_size) :: H_init(6)
  
  namelist /set_parm/ nt_asm, nt_prd, obs_interval
  namelist /da_setting/ da_method, alpha
  namelist /intg_setting/ intg_method
  namelist /ensemble_size/ mems
  namelist /initial_score/ x_tinit, y_tinit, z_tinit, x_sinit, y_sinit, z_sinit
  namelist /initial_matrix/ Pf_init, B_init, R_init, Kg_init, H_init
  namelist /output/ output_file, output_file_error_covariance, opt_veach

  read(5, nml=set_parm, iostat = ierr)
  read(5, nml=da_setting, iostat = ierr)
  read(5, nml=intg_setting, iostat = ierr)
  if ( trim(da_method) == 'EnKF' ) then
    read(5, nml=ensemble_size, iostat = ierr)
  end if 
  read(5, nml=initial_score, iostat = ierr)
  read(5, nml=initial_matrix, iostat = ierr)
  read(5, nml=output, iostat = ierr)
  ! name list io check
  if (ierr < 0 ) then
    write(6,*) '   Msg : Main[ .sh /  @namelist ] '
    write(6,*) '   Not found namelist.        '
    write(6,*) '   Use default values.        '
  else if (ierr > 0) then
    write(6,*) trim(output_file), opt_veach
    write(6,*) '   Msg : Main[ .sh /  @namelist ] '
    write(6,*) '   *** Warning : Not appropriate names in namelist !! Check !!'
    write(6,*) '   Stop : lorenz63_main.f90              '
    stop
  end if

  ! +++ display namelist
  write(6,*) 'Data assimlation method :: ', da_method
  write(6,*) 'Integral method         :: ', intg_method

  allocate(x_true(0:nt_asm+nt_prd), y_true(0:nt_asm+nt_prd), z_true(0:nt_asm+nt_prd))
  allocate(x_sim(0:nt_asm+nt_prd), y_sim(0:nt_asm+nt_prd), z_sim(0:nt_asm+nt_prd))
  allocate(x_da(0:nt_asm+nt_prd), y_da(0:nt_asm+nt_prd), z_da(0:nt_asm+nt_prd))
  allocate(x_da_m(0:nt_asm, mems), y_da_m(0:nt_asm, mems), z_da_m(0:nt_asm, mems))
  allocate(x_prtb(mems), y_prtb(mems), z_prtb(mems))

  allocate(x_obs(0:nt_asm/obs_interval), y_obs(0:nt_asm/obs_interval))
  allocate(obs_ens(2, mems))
  allocate(x_innov(mems), y_innov(mems))
  allocate(obs_chr(2, 0:nt_asm+nt_prd))

  !----------------------------------------------------------------------
  ! +++ initial setting
  
  x_true(0) = x_tinit; y_true(0) = y_tinit; z_true(0) = z_tinit
  x_sim(0)  = x_sinit; y_sim(0)  = y_sinit; z_sim(0)  = z_sinit
  
  Pf(1,1) = Pf_init(1); Pf(1,2) = Pf_init(2); Pf(1,3) = Pf_init(3)
  Pf(2,1) = Pf_init(4); Pf(2,2) = Pf_init(5); Pf(2,3) = Pf_init(6)
  Pf(3,1) = Pf_init(7); Pf(3,2) = Pf_init(8); Pf(3,3) = Pf_init(9)
  
  Pa = Pf
  
  R(1,1) = R_init(1); R(1,2) = R_init(2)
  R(2,1) = R_init(3); R(2,2) = R_init(4)
  
  Kg(1,1) = Kg_init(1); Kg(1,2) = Kg_init(2)
  Kg(2,1) = Kg_init(3); Kg(2,2) = Kg_init(4)
  Kg(3,1) = Kg_init(5); Kg(3,2) = Kg_init(6)

  H(1,1) = H_init(1); H(1,2) = H_init(2); H(1,3) = H_init(3)
  H(2,1) = H_init(4); H(2,2) = H_init(5); H(2,3) = H_init(6)

  I = 0.0d0
  forall ( il=1:nx )  I(il, il) = 1.0d0
  
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
    if ( trim(intg_method) == 'Euler' ) then
      x_true(it) = x_true(it-1) + dt * x_k(1)
      y_true(it) = y_true(it-1) + dt * y_k(1)
      z_true(it) = z_true(it-1) + dt * z_k(1)
      
      !------------------------------------------------------- 
      ! +++ Runge-Kutta method
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
      
      write(6,*) 'time_step, x_obs, y_obs', it, x_obs(it/obs_interval), y_obs(it/obs_interval)
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
    if ( trim(intg_method) == 'Euler' ) then
      x_sim(it) = x_sim(it-1) + dt * x_k(1)
      y_sim(it) = y_sim(it-1) + dt * y_k(1)
      z_sim(it) = z_sim(it-1) + dt * z_k(1)
      
      !------------------------------------------------------- 
      ! +++ Runge-Kutta method
    else if ( trim(intg_method) == 'Runge-Kutta' ) then 
      
      call Lorenz63_Runge_Kutta(               &
        x_sim(it-1), y_sim(it-1), z_sim(it-1), & ! IN
        x_sim(it), y_sim(it), z_sim(it)        & ! OUT
        )
        
      end if
    end do
    
    ! --- Sec4. Data assimilation
    if ( da_method == 'KF' ) then
      x_da(0) = x_sim(0)
      y_da(0) = y_sim(0)
      z_da(0) = z_sim(0)
      
      do it = 1, nt_asm
        if ( opt_veach ) then
          call write_error_covariance_matrix(it, nt_asm, Pa)
        end if  

        write(6,*) 'Data assim. time step: ', it
        
        ! 4.1: Time integration
        call cal_Lorenz(                         &
        x_da(it-1), y_da(it-1), z_da(it-1),      & ! IN
        x_k(1), y_k(1), z_k(1)                   & ! OUT
        )
        
        !------------------------------------------------------- 
        ! +++ Euler method
        if ( trim(intg_method) == 'Euler' ) then
          x_da(it) = x_da(it-1) + dt * x_k(1)
          y_da(it) = y_da(it-1) + dt * y_k(1)
          z_da(it) = z_da(it-1) + dt * z_k(1)

          !------------------------------------------------------- 
          ! +++ Runge-Kutta method
        else if ( trim(intg_method) == 'Runge-Kutta' ) then 

          call Lorenz63_Runge_Kutta(             &
          x_da(it-1), y_da(it-1), z_da(it-1),    & ! IN
          x_da(it), y_da(it), z_da(it)           & ! OUT
          )
        end if

        if (mod(it, obs_interval) == 0) then
          !------------------------------------------------------- 
          ! 4.2: Kalman fileter
          !------------------------------------------------------- 
          ! +++ 4.2.1 State Transient Matrix
          ! >> original form
          ! JM(1,1) = 1.0d0 - dt*sig;         JM(1,2) = dt*sig;         JM(1,3) = 0.0d0
          ! JM(2,1) = dt*(gamm - z_da(it-1)); JM(2,2) = 1.0d0 - dt;     JM(2,3) = -dt*x_da(it-1)
          ! JM(3,1) = dt*y_da(it-1);          JM(3,2) = dt*x_da(it-1);  JM(3,3) = 1.0d0 - dt*b
        
          do il = 1, nx
            delta = 1.0d0-3
            call Lorenz63_Runge_Kutta( &
              x_da(it-1)+delta, &
              y_da(it-1)+delta, &
              z_da(it-1)+delta, &
              x_e, y_e, z_e     &
            )
            JM(1,il) = ( x_e - x_da(it) )/ delta
            JM(2,il) = ( y_e - y_da(it) )/ delta
            JM(3,il) = ( z_e - z_da(it) )/ delta
          end do

          call confirm_matrix(JM,nx,nx)

          Pf = matmul(JM, matmul(Pa, transpose(JM)))*(1.0d0 + alpha)
          !------------------------------------------------------- 
          ! >> 4.2.3 Kalman gain: Weighting of model result and obs.
          ! (Note) Observation in x,y ----> component 2 (x,y)
          ! calculate inverse matrix @R_HPHt
          ! *** ref:
          ! http://www.rcs.arch.t.u-tokyo.ac.jp/kusuhara/tips/linux/fortran.html
          
          HPHt = matmul(H, matmul(Pf, transpose(H)))
          R_HPHt = R + HPHt
          
          !------------------------------------------------------- 
          ! +++ inverse matrix calculate for 2x2 on formula
          call inverse_matrix_for2x2(R_HPHt, inv_matrix)
          
          eye_matrix = matmul(R_HPHt, inv_matrix)
          !write(6,*) ' +++ confirm eye matrix calculate ... '
          !write(6,*) eye_matrix
          Kg = matmul(matmul(Pf, transpose(H)), inv_matrix)
          
          !------------------------------------------------------- q
          ! >> 4.2.4 calculate innovation and correlation
          ! +++ Kalman Filter Main equation
          x_a(1,1) = x_da(it);  x_a(2,1) = y_da(it); x_a(3,1) = z_da(it)
          yt_obs(1,1) = x_obs(it/obs_interval); yt_obs(2,1) = y_obs(it/obs_interval)

          hx = matmul(H, x_a)
          hdxf = yt_obs - hx
          x_a = x_a + matmul(Kg, hdxf)
          call confirm_matrix(Kg, nx, ny)
          x_da(it) = x_a(1,1); y_da(it) = x_a(2,1); z_da(it) = x_a(3,1)
          
          ! >> 4.2.5 analysis error covariance matrix
          Pa = matmul((I - matmul(Kg, H)), Pf)
          call confirm_matrix(Pa, nx, nx)
        end if
      end do
      
    else if ( da_method == 'EnKF' ) then
      ! +++ initial setting
      do imem = 1, mems
        call gaussian_noise(sqrt(Pa(1,1)), Gnoise)
        x_da_m(0, imem) = x_sim(0) + Gnoise
        call gaussian_noise(sqrt(Pa(2,2)), Gnoise)
        y_da_m(0, imem) = y_sim(0) + Gnoise
        call gaussian_noise(sqrt(Pa(3,3)), Gnoise)
        z_da_m(0, imem) = z_sim(0) + Gnoise
      end do
      
      x_da(0) = sum(x_da_m(0, 1:mems))/mems 
      y_da(0) = sum(x_da_m(0, 1:mems))/mems 
      z_da(0) = sum(x_da_m(0, 1:mems))/mems
      
      do it = 1, nt_asm
        if ( opt_veach ) then
          call write_error_covariance_matrix(it, nt_asm, Pf)
        end if 
        ! 4.1: Time integration
        do imem = 1, mems
          ! 4.1: Time integration
          call cal_Lorenz(                                              &
          x_da_m(it-1, imem), y_da_m(it-1, imem), z_da_m(it-1, imem),   & ! IN
          x_k(1), y_k(1), z_k(1)                                        & ! OUT
          )
          
          !------------------------------------------------------- 
          ! +++ Euler method
          if ( trim(intg_method) == 'Euler' ) then
            x_da_m(it, imem) = x_da_m(it-1, imem) + dt * x_k(1)
            y_da_m(it, imem) = y_da_m(it-1, imem) + dt * y_k(1)
            z_da_m(it, imem) = z_da_m(it-1, imem) + dt * z_k(1)
            
            !------------------------------------------------------- 
            ! +++ Runge-Kutta method
          else if ( trim(intg_method) == 'Runge-Kutta' ) then 
            
            call Lorenz63_Runge_Kutta(                                     &
            x_da_m(it-1, imem), y_da_m(it-1, imem), z_da_m(it-1, imem),    & ! IN
            x_da_m(it, imem), y_da_m(it, imem), z_da_m(it, imem)           & ! OUT
            )
          end if
        end do
        
        if(mod(it, obs_interval) == 0) then
          x_da(it) = sum(x_da_m(it, 1:mems))/mems
          y_da(it) = sum(y_da_m(it, 1:mems))/mems
          z_da(it) = sum(z_da_m(it, 1:mems))/mems
          Pf = 0.0d0
          do imem = 1, mems
            x_prtb(imem) = x_da_m(it, imem) - x_da(it)
            y_prtb(imem) = y_da_m(it, imem) - y_da(it)
            z_prtb(imem) = z_da_m(it, imem) - z_da(it)
            
            !------------------------------------------------------- 
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
          end do
          
          HPHt = matmul(H, matmul(Pf, transpose(H)))
          R_HPHt = R + HPHt
          
          !------------------------------------------------------- 
          ! +++ inverse matrix calculate for 2x2 on formula
          call inverse_matrix_for2x2(R_HPHt, inv_matrix)
          
          eye_matrix = matmul(R_HPHt, inv_matrix)
          !write(6,*) ' +++ confirm eye matrix calculate ... '
          !write(6,*) eye_matrix
          Kg = matmul(matmul(Pf, transpose(H)), inv_matrix)
          call confirm_matrix(Kg, nx, ny)
          
          do imem = 1, mems
            !------------------------------------------------------- 
            ! +++ Pertuturbed observation method (PO)
            call gaussian_noise(sqrt(R(1,1)), Gnoise)
            x_innov(imem) = x_obs(it/obs_interval) + Gnoise
            obs_ens(1,imem) = x_innov(imem)
            call gaussian_noise(sqrt(R(2,2)), Gnoise)
            y_innov(imem) = y_obs(it/obs_interval) + Gnoise
            obs_ens(2,imem) = y_innov(imem)
            
            x_a(1,1) = x_da_m(it, imem); x_a(2,1) = y_da_m(it, imem); x_a(3,1) = z_da_m(it, imem)
            yt_obs(1,1) = obs_ens(1,imem); yt_obs(2,1) = obs_ens(2,imem)
            hx = matmul(H, x_a)
            hdxf = yt_obs - hx

            x_a = x_a + matmul(Kg, hdxf)
            x_da_m(it, imem) = x_a(1,1); y_da_m(it, imem) = x_a(2,1); z_da_m(it, imem) = x_a(3,1)
          end do

        end if
        
        x_da(it) = sum(x_da_m(it, 1:mems))/mems
        y_da(it) = sum(y_da_m(it, 1:mems))/mems
        z_da(it) = sum(z_da_m(it, 1:mems))/mems

        ! >> 4.2.5 analysis error covariance matrix
        Pa = Pf - matmul(matmul(Kg, H), Pf)
        Pf = Pa
      end do
 
    end if
  
  ! --- Sec5. Prediction after Data assimilation
    do it = nt_asm+1, nt_asm+nt_prd
      write(6,*) 'Data assim. time step: ', it
      ! forward time step
      
      call cal_Lorenz(                           &
      x_da(it-1), y_da(it-1), z_da(it-1),        & ! IN
      x_k(1), y_k(1), z_k(1)                     & ! OUT
      )
      
      !------------------------------------------------------- 
      ! +++ Euler method
      if ( trim(intg_method) == 'Euler' ) then
        x_da(it) = x_da(it-1) + dt * x_k(1)
        y_da(it) = y_da(it-1) + dt * y_k(1)
        z_da(it) = z_da(it-1) + dt * z_k(1)
        
        !------------------------------------------------------- 
        ! +++ Runge-Kutta method
      else if ( trim(intg_method) == 'Runge-Kutta' ) then 
        
        call Lorenz63_Runge_Kutta(               &
        x_da(it-1), y_da(it-1), z_da(it-1),    & ! IN
        x_da(it), y_da(it), z_da(it)           & ! OUT
        )
        
      end if
    end do
    
    ! --- Sec6. Writing OUTPUT
    if ( opt_veach ) then
      obs_chr(:, 0:nt_asm+nt_prd) = undef
      do it = 1, nt_asm
        if (mod(it, obs_interval) == 0) then
          obs_chr(1, it) = x_obs(it/obs_interval)
          obs_chr(2, it) = y_obs(it/obs_interval)
        end if
      end do
      
      open (1, file=trim(output_file), status='replace')
      write(1,*) 'timestep, x_true, y_true, z_true, x_sim, y_sim, z_sim, x_da, y_da, z_da, x_obs, y_obs'
      do it = 0, nt_asm+nt_prd
        if (mod(it, output_interval) == 0) then
          write(linebuf, '(f5.2, ",", 9(f12.7, ","), F7.2, ",", F7.2)')  & 
            dt*it,                                      &
            x_true(it), y_true(it), z_true(it),         &
            x_sim(it), y_sim(it), z_sim(it),            &
            x_da(it), y_da(it), z_da(it),               &
            obs_chr(1, it), obs_chr(2, it)
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

  subroutine Lorenz63_Runge_Kutta(  &
    x_in, y_in, z_in,               & ! IN: previous step score 
    x_out, y_out, z_out             & ! OUT: Runge-Kutta method score 
  )
    
    use lorenz63_prm
    implicit none
    
    real(r_size), intent(in)    :: x_in, y_in, z_in
    real(r_size), intent(out)   :: x_out, y_out, z_out
    
    x_cal(1) = x_in + 0.5*x_k(1)*dt
    y_cal(1) = y_in + 0.5*y_k(1)*dt
    z_cal(1) = z_in + 0.5*z_k(1)*dt
    
    call cal_Lorenz(                         &
    x_cal(1), y_cal(1), z_cal(1),            & ! IN
    x_k(2), y_k(2), z_k(2)                   & ! OUT
    )
    
    x_cal(2) = x_in + 0.5*x_k(2)*dt 
    y_cal(2) = y_in + 0.5*y_k(2)*dt 
    z_cal(2) = z_in + 0.5*z_k(2)*dt
    
    call cal_Lorenz(                         &
    x_cal(2), y_cal(2), z_cal(2),            & ! IN
    x_k(3), y_k(3), z_k(3)                   & ! OUT
    )
    
    x_cal(3) = x_in + x_k(3)*dt
    y_cal(3) = y_in + y_k(3)*dt
    z_cal(3) = z_in + z_k(3)*dt
    
    call cal_Lorenz(                         &
    x_cal(3), y_cal(3), z_cal(3),            & ! IN
    x_k(4), y_k(4), z_k(4)                   & ! OUT
    )
    
    x_out = x_in + dt * (x_k(1) + 2*x_k(2) + 2*x_k(3) + x_k(4)) / 6.0d0
    y_out = y_in + dt * (y_k(1) + 2*y_k(2) + 2*y_k(3) + y_k(4)) / 6.0d0
    z_out = z_in + dt * (z_k(1) + 2*z_k(2) + 2*z_k(3) + z_k(4)) / 6.0d0
    
    return
  end subroutine Lorenz63_Runge_Kutta
  
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

  subroutine write_error_covariance_matrix(it, last_step, error_covariance_matrix)
    implicit none
    integer, intent(in)       :: it, last_step
    real(r_size), intent(in)  :: error_covariance_matrix(3,3)

    if ( it == 1 ) then
      open(2, file=trim(output_file_error_covariance), status='replace')
      write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        write(6,*) '+++ err covariance matrix 1st. step'
        !write(6,*) error_covariance_matrix(:,:)
        
      else if ( it /= 1 .and. it /= last_step) then
        write(6,*) '+++ err covariance matrix 2nd. step ~'
        write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        !write(6,*) error_covariance_matrix(:,:)
        
      else if ( it == last_step ) then
        write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        write(6,*) '+++ err covariance matrix last step '
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
    print *, "==============================="
  end subroutine


end program lorenz63
