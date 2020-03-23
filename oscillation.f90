! Created on 2020.3.21
! @author: Toyo_Daichi
! (ref: http://www.itonwp.sci.u-ryukyu.ac.jp/itokosk.html) 

program oscillation
  implicit none

  ! --- setting Parameters
  integer, parameter  :: nt_asm       = 400 ! Period of data assimilation
  integer, parameter  :: nt_prd       = 400 ! Period of prediction
  integer, parameter  :: obs_interval = 40  ! Interval of observation
  integer, parameter  :: r_size       = 8   ! Byte 

  real(r_size), parameter  :: mass = 1.0d0
  real(r_size), parameter  :: k    = 0.5d0
  real(r_size), parameter  :: dump = 0.3d0  ! Damping coefficinet
  real(r_size), parameter  :: dt   = 1.0d-2 ! Time step
  real(r_size), parameter  :: pi   = 3.14159265358979d0
  
  ! --- Physical variable
  real(r_size) :: x_true(0:nt_asm+nt_prd), v_true(0:nt_asm+nt_prd)
  real(r_size) :: x_sim(0:nt_asm+nt_prd), v_sim(0:nt_asm+nt_prd)
  
  real(r_size) :: x_da(0:nt_asm+nt_prd), v_da(0:nt_asm+nt_prd)
  real(r_size) :: x_obs(0:nt_asm/obs_interval)
  
  ! --- Matrix(element 1:x, 2:v)
  real(r_size) :: M(2,2)   ! state transient matrix
  real(r_size) :: Pf(2,2)  ! Forecast error convariance matrix
  real(r_size) :: Pa(2,2)  ! Analysis error convariance matrix
  real(r_size) :: R(1,1)   ! Observation error convariance matrix
  real(r_size) :: Kg(2,1)  ! Kalman gain
  real(r_size) :: H(1,2)   ! Observation operator
  
  ! --- Output control
  integer, parameter  :: output_interval = 20
  character(7)        :: obs_chr(0:nt_asm+nt_prd)
  character(256)      :: linebuf
  
  ! --- Working variable
  integer :: it
  integer :: irec
  integer :: ierr
  integer :: iflag
  real(r_size) :: x_innov
  real(r_size) :: Ptmp(2,2)
  real(r_size) :: noise1, noise2, Gnoise ! Gaussian noise
  !======================================================================
  ! Data assimilation
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist
  !----------------------------------------------------------------------
  real(r_size) :: x_tinit, v_tinit
  real(r_size) :: x_sinit, v_sinit
  real(r_size) :: Pf_init(4)
  real(r_size) :: R_init
  real(r_size) :: Kg_init(2)
  real(r_size) :: H_init(2)
  
  namelist /initial_osc/ &
  x_tinit, v_tinit, &
  x_sinit, v_sinit, &
  Pf_init, &
  R_init,  &
  Kg_init, &
  H_init
  
  open(1, file='initial_Osc.cnf')
    read(1, nml=initial_osc, iostat=ierr)
  close(1)
  
  ! namelist check for debug
  if (ierr < 0) then
    write(6,*) '   Msg : Main[ main /  initial_osc ] '
    write(6,*) '   Not found namelist.        '
    write(6,*) '   Use default values.        '
  else if (ierr > 0) then
    write(6,*) '   Msg : Main[ main /  initial_osc ] '
    write(6,*) '   *** Warning : Not appropriate names in namelist !! Check !!'
    write(6,*) '   Stop : oscillation.f90              '
    stop
  end if

  x_true(0) = x_tinit; v_true(0) = v_tinit
  x_sim(0)  = x_sinit; v_sim(0)  = v_sinit
  
  Pf(1,1) = Pf_init(1); Pf(1,2)=Pf_init(2); Pf(2,1)=Pf_init(3); Pf(2,2)=Pf_init(4)
  Pa = Pf
  
  R(1,1) = R_init
  
  Kg(1:2,1:1) = Kg_init(1)
  H(1,1) = H_init(1); H(1,2) = H_init(2)
  
  ! --- Initialization of random number generator
  call random_seed()
  
  ! --- Sec2. True field and observations
  do it = 1, nt_asm+nt_prd
    ! forward time step
    x_true(it) = x_true(it-1) + dt * v_true(it-1)
    v_true(it) = -(k * dt / mass) * x_true(it-1) + (1.0d0 - dump * dt / mass ) * v_true(it-1)
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
  
  ! --- Sec3. Simulation run without DA
  do it = 1, nt_asm + nt_prd
    x_sim(it) = x_sim(it-1) + dt * v_sim(it-1)  
    v_sim(it) = -(k * dt / mass) * x_sim(it-1) + (1.0d0 - dump * dt / mass ) * v_sim(it-1)
  end do

  ! --- Sec4. Data assimilation
  x_da(0) = x_sim(0) 
  v_da(0) = v_sim(0)
  
  do it = 1, nt_asm
    ! 4.1: Time integration
    x_da(it) = x_da(it-1) + dt*v_da(it-1)
    v_da(it) = -(k * dt / mass) * x_da(it-1) + (1.0d0 - dump * dt / mass ) * v_da(it-1)
    ! 4.2: Kalman filter
    ! >> 4.2.1 State Transient Matix
    M(1,1) = 1.0d0
    M(1,2) = dt
    M(2,1) = -k * dt / mass
    M(2,2) = 1.0d0 - dump * dt / mass

    ! >> 4.2.2 Lyapunov equation: Obtain Pf
    Ptmp = transpose(M)
    Ptmp = matmul(Pf, Ptmp)
    Pf   = matmul(M, Ptmp)
    if (mod(it, obs_interval) == 0) then
      ! >> 4.2.3 Kalman gain: Weighting of model result and obs.
      ! (Note) Observation only in x ----> component 1(x) only 
      !        In this case, Inverse matrix ----> scalar matrix
      Kg(1,1) = Pf(1,1) / (R(1,1) + Pf(1,1)) 
      Kg(2,1) = Pf(2,1) / (R(1,1) + Pf(1,1)) 
      ! >> 4.2.4 calculate innovation and correction
      x_innov  = x_obs(it / obs_interval) - x_da(it)
      x_da(it) = x_da(it) + Kg(1,1) * x_innov  
      v_da(it) = v_da(it) + Kg(2,1) * x_innov
      ! >> 4.2.5 analysis error covariance matrix
      Pa = Pf - matmul(matmul(Kg, H), Pf)
      Pf = Pa
      write(6,'(A, F8.2, A, 5F10.3)') 'time = ', dt*it, ', Pa =,', Pa
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
  write(6,*) '  -------- Method: Kalman Fileter    --------  '
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

  open (1, file='./output/oscillation_KF.csv', status='replace')
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

  subroutine del_spaces(s)
    character (*), intent (inout) :: s
    character (len=len(s)) tmp
    integer i, j
    j = 1
    do i = 1, len(s)
      if (s(i:i)==' ') cycle
      tmp(j:j) = s(i:i)
      j = j + 1
    end do
    s = tmp(1:j-1)
  end subroutine del_spaces

end program oscillation