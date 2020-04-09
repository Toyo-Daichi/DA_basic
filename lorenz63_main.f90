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

  real(r_size), parameter  :: dt   = 1.0d-2 ! Time step

  character(8) :: da_method
  
  ! --- Physical variable
  real(r_size), allocatable :: x_true(:), y_true(:), z_true(:)
  real(r_size), allocatable :: x_sim(:), y_sim(:), z_sim(:)
  
  real(r_size), allocatable :: x_da(:), y_da(:), z_da(:)
  real(r_size), allocatable :: x_obs(:)

  real(r_size), allocatable :: x_da_m(:, :), y_da_m(:, :), z_da_m(:, :)
  real(r_size), allocatable :: x_prtb(:), y_prtb(:), z_prtb(:)

  ! --- For calculation, Runge-Kutta method
  real(r_size) :: x_k(4), y_k(4), z_k(4)

  ! --- Matrix(element 1:x, 2:y, 3:z)
  ! +++ default setting
  ! Pf = (  Pxx: 1.0  Pxy: 0.0 Pxz: 0.0
  !         Pyx: 0.0  Pyy: 1.0 Pyz: 0.0 
  !         Pzx: 0.0  Pzy: 0.0 Pzz: 1.0 )
  
  real(r_size) :: M(3,3)   ! state transient matrix
  real(r_size) :: Pf(3,3)  ! Forecast error convariance matrix (in KF, EnKF)
  real(r_size) :: Pa(3,3)  ! Analysis error convariance matrix
  real(r_size) :: R(2,2)   ! Observation error convariance matrix
  real(r_size) :: Kg(3,2)  ! Kalman gain
  real(r_size) :: H(2,3)   ! Observation operator
  
  !  >> -----------------------------------------------------------------
  ! +++ 4Dvar_Ajoint model : Cost function and ajoint variable
  real(r_size) :: B(3,3)   ! Background error convariance matrix
  real(r_size) :: J        ! Cost function
  real(r_size) :: Jold     ! Cost function in a previous iteration
  real(r_size) :: x_b, v_b ! First guess of initial value

  real(r_size), allocatable :: adx(:), ady(:), adz(:)
  real(r_size), allocatable :: x_save(:), y_save(:), z_save(:)
  real(r_size), allocatable :: x_tmp(:), y_tmp(:), z_tmp(:)
  real(r_size), parameter   :: alpx = 0.02         ! Coefficient for minization
  real(r_size), parameter   :: alpy = 0.02
  real(r_size), parameter   :: alpz = 0.02
  integer, parameter        :: iter_max = 500      ! maxmum number of iteration
  real(r_size), parameter   :: cond_iter = 1.0d-4  ! condition for iteration end ***(Jold - J)/Jold

  ! --- Output control
  character(7),allocatable :: obs_chr(:)
  integer, parameter       :: output_interval = 20
  !logical                  :: opt_beach = .false.
  character(256)           :: linebuf
  character(256)           :: output_file

  ! --- Working variable
  integer :: it
  integer :: mems, imem
  integer :: ierr
  integer :: iflag
  integer :: iter
  real(r_size) :: x_innov
  real(r_size) :: Ptmp(3,3)
  real(r_size) :: noise1, noise2, Gnoise ! Gaussian noise

  !======================================================================
  ! Data assimilation
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------
  real(r_size) :: x_tinit, y_tinit, z_tinit
  real(r_size) :: x_sinit, y_sinit, z_sinit
  real(r_size) :: Pf_init(9), B_init(9)
  real(r_size) :: R_init(4)
  real(r_size) :: Kg_init(6)
  real(r_size) :: H_init(6)
  
  namelist /set_parm/ nt_asm, nt_prd, obs_interval
  namelist /da_setting/ da_method
  namelist /ensemble_size/ mems
  namelist /initial_osc/ x_tinit, y_tinit, z_tinit,
                         x_sinit, y_sinit, z_sinit
  namelist /initial_matrix/ Pf_init, B_init, R_init, Kg_init, H_init
  namelist /output/ output_file ! opt_beach
 
  read(5, nml=set_parm, iostat = ierr)
  read(5, nml=da_setting, iostat = ierr)
  if ( trim(da_method) == 'EnKF' ) then
    read(5, nml=ensemble_size, iostat = ierr)
  end if 
  read(5, nml=initial_osc, iostat = ierr)
  read(5, nml=initial_que, iostat = ierr)
  read(5, nml=output, iostat = ierr)
  ! name list io check
  if (ierr < 0 ) then
    write(6,*) '   Msg : Main[ .sh /  @namelist ] '
    write(6,*) '   Not found namelist.        '
    write(6,*) '   Use default values.        '
  else if (ierr > 0) then
    write(6,*) '   Msg : Main[ .sh /  @namelist ] '
    write(6,*) '   *** Warning : Not appropriate names in namelist !! Check !!'
    write(6,*) '   Stop : oscillation.f90              '
    stop
  end if

  allocate(x_true(0:nt_asm+nt_prd), y_true(0:nt_asm+nt_prd), z_true(0:nt_asm+nt_prd))
  allocate(x_sim(0:nt_asm+nt_prd), y_sim(0:nt_asm+nt_prd), z_sim(0:nt_asm+nt_prd))
  allocate(x_da(0:nt_asm+nt_prd), y_da(0:nt_asm+nt_prd), z_da(0:nt_asm+nt_prd))

  allocate(x_da_m(0:nt_asm, mems), y_da_m(0:nt_asm, mems), z_da_m(0:nt_asm, mems))
  allocate(x_tmp(0:nt_asm), y_tmp(0:nt_asm), z_tmp(0:nt_asm))
  allocate(x_prtb(mems), y_prtb(mems), z_prtb(mems))
  allocate(adx(0:nt_asm), ady(0:nt_asm), adz(0:nt_asm))
  allocate(x_save(0:nt_asm), y_save(0:nt_asm), z_save(0:nt_asm))
  allocate(x_obs(0:nt_asm/obs_interval))
  allocate(obs_chr(0:nt_asm+nt_prd))

  !----------------------------------------------------------------------
  ! +++ initial setting
  
  x_true(0) = x_tinit; y_true(0) = y_tinit; z_true(0) = z_tinit
  x_sim(0)  = x_sinit; y_sim(0)  = y_sinit; z_sim(0)  = z_sinit
  
  Pf(1,1) = Pf_init(1); Pf(1,2)=Pf_init(2); Pf(1,3)=Pf_init(3)
  Pf(2,1) = Pf_init(4); Pf(2,2)=Pf_init(5); Pf(2,3)=Pf_init(6)
  Pf(3,1) = Pf_init(7); Pf(3,2)=Pf_init(8); Pf(2,3)=Pf_init(9)
  
  Pa = Pf
  
  R(1,1) = R_init(1); R(1,2) = R_init(2)
  R(2,1) = R_init(3); R(2,2) = R_init(4)
  
  Kg(1,1) = Kg_init(1); Kg(1,2) = Kg_init(2)
  Kg(2,1) = Kg_init(1); Kg(2,2) = Kg_init(2)
  Kg(3,1) = Kg_init(1); Kg(3,2) = Kg_init(2)

  H(1,1) = H_init(1); H(1,2) = H_init(2); H(1,3) = H_init(3)
  H(2,1) = H_init(4); H(2,2) = H_init(5); H(2,3) = H_init(6)
  
  ! --- Initialization of random number generator
  call random_seed()
  
  ! --- Sec2. True field and observations(Runge-Kutta method)
  do it = 1, nt_asm+nt_prd
    ! forward time step

    call cal_Lorenz(                           &
      x_true(it-1), y_true(it-1), z_true(it-1),   & ! IN
      x_k(1), y_k(1), z_k(1)                   & ! OUT
    )

    x_cal(1) = x_true(it-1) + 0.5*x_k(1)*dt
    y_cal(1) = y_true(it-1) + 0.5*y_k(1)*dt
    z_cal(1) = z_true(it-1) + 0.5*z_k(1)*dt

    call cal_Lorenz(                           &
    x_cal(1), y_cal(1), z_cal(1),              & ! IN
    x_k(2), y_k(2), z_k(2)                     & ! OUT
    )

    x_cal(2) = x_true(it-1) + 0.5*x_k(2)*dt 
    y_cal(2) = y_true(it-1) + 0.5*y_k(2)*dt 
    z_cal(2) = z_true(it-1) + 0.5*z_k(2)*dt

    call cal_Lorenz(                           &
    x_cal(2), y_cal(2), z_cal(2),              & ! IN
    x_k(3), y_k(3), z_k(3)                     & ! OUT
    )

    x_cal(3) = x_true(it-1) + x_k(3)*dt
    y_cal(3) = y_true(it-1) + y_k(3)*dt
    y_cal(3) = z_true(it-1) + z_k(3)*dt

    call cal_Lorenz(                           &
    x_cal(3), y_cal(3), z_cal(3),              & ! IN
    x_k(4), y_k(4), z_k(4)                     & ! OUT
    )

    x_true(it) = x_true(it-1) + dt / 6.0d0 * ( x_k(1) + 2*x_k(2) + 2*x_k(3) + x_k(4))
    y_true(it) = y_true(it-1) + dt / 6.0d0 * ( y_k(1) + 2*y_k(2) + 2*y_k(3) + y_k(4))
    z_true(it) = z_true(it-1) + dt / 6.0d0 * ( z_k(1) + 2*z_k(2) + 2*z_k(3) + z_k(4))

    ! making observations
    if ((mod(it, obs_interval) == 0) .and. (it <= nt_asm)) then
      ! Generate Gaussian Noise (Gnoise) from uniform random number
      ! based on Box-Muller method
      call random_number(noise1)
      call random_number(noise2)
      Gnoise = sqrt(R(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
      ! Generate observation by adding Gaussian noise to true value
      x_obs(it/obs_interval) = x_true(it) + Gnoise
    end if
  end do
  
  
  ! --- Sec5. Prediction after Data assimilation
  do it = nt_asm+1, nt_asm + nt_prd
    x_da(it) = x_da(it-1) + dt*v_da(it-1)
    v_da(it) = -(k * dt / mass) * x_da(it-1) + (1.0d0 - dump * dt / mass ) * v_da(it-1)
  end do

  !----------------------------------------------------------------------
  ! +++ OUTPUT Main
  !----------------------------------------------------------------------
  ! Identical Twin Experiment
  
  obs_chr(0:nt_asm+nt_prd) = 'None'
  do it = 1, nt_asm
    if (mod(it, obs_interval) == 0) then
      write(obs_chr(it), '(F7.3)') x_obs(it/obs_interval)
    end if
  end do
  
  write(6,*) '  -------- Identical Twin Experiment --------  '
  write(6,*) '  ***  Method: ', da_method
  write(6,*)
  write(6,'(A,F7.2,A,F7.2)') 'Assimilation Period: t= ', 0.0, '-', dt*nt_asm
  write(6,'(A,F7.2,A,F7.2)') 'Prediction   Period: t= ', 0.0, '-', dt*(nt_asm+nt_prd)
  write(6,*)
  
  write(6,*) '  >>> Assimilation Period: x '
  write(6,*) ' [Time]   [True]   [No assim] [Assim] [Observation] '
  do it = 0, nt_asm
    if (mod(it, output_interval) == 0) then
      write(6, '(F7.2, 3F10.3, 4X, A)') dt*it, x_true(it), x_sim(it), x_da(it), obs_chr(it)
    end if
  end do

  write(6,*) '  >>> Prediction Period: x '
  do it = nt_asm+1, nt_asm+nt_prd
    if (mod(it, output_interval) == 0) then
      write(6, '(F7.2, 3F10.3)') dt*it, x_true(it), x_sim(it), x_da(it)
    end if
  end do

  write(6,*) '  >>> Assimilation Period: v '
  write(6,*) ' [Time]   [True]   [No assim] [Assim] '
  do it = 0, nt_asm
    if (mod(it, output_interval) == 0) then
      write(6, '(F7.2, 3F10.3)') dt*it, x_true(it), x_sim(it), x_da(it)
    end if
  end do

  write(6,*) '  >>> Prediction Period: v '
  do it = nt_asm+1, nt_asm+nt_prd
    if (mod(it, output_interval) == 0) then
      write(6, '(F7.2, 3F10.3)') dt*it, v_true(it), v_sim(it), v_da(it)
    end if
  end do

  open (1, file=trim(output_file), status='replace')
    write(1,*) 'timestep, x_true, x_sim, x_da, v_true, v_sim, v_da, obs_data'
    do it = 0, nt_asm+nt_prd
      if (mod(it, output_interval) == 0) then
        write(linebuf, *) dt*it, ',', x_true(it), ',', x_sim(it), ',', x_da(it), ',', &
                          v_true(it), ',', v_sim(it), ',', v_da(it), ',', obs_chr(it)
        call del_spaces(linebuf)
        write (1, '(a)') trim(linebuf)
      end if
    end do
  close(1)

  contains

  subroutine del_spaces(space)
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

end program lorenz63
