! Created on 2020.3.21
! @author: Toyo_Daichi
! (ref: http://www.itonwp.sci.u-ryukyu.ac.jp/itokosk.html) 

program oscillation

  use kinddef

  implicit none

  ! --- setting Parameters
  integer :: nt_asm       ! Period of data assimilation
  integer :: nt_prd       ! Period of prediction
  integer :: obs_interval ! Interval of observation

  real(r_size), parameter  :: mass = 1.0d0
  real(r_size), parameter  :: k    = 0.5d0
  real(r_size), parameter  :: dump = 0.3d0  ! Damping coefficinet
  real(r_size), parameter  :: dt   = 1.0d-2 ! Time step
  real(r_size), parameter  :: pi   = 3.14159265358979d0

  character(8) :: da_method
  
  ! --- Physical variable
  real(r_size), allocatable :: x_true(:), v_true(:)
  real(r_size), allocatable :: x_sim(:), v_sim(:)
  
  real(r_size), allocatable :: x_da_m(:, :), v_da_m(:, :)
  real(r_size), allocatable :: x_prtb(:), v_prtb(:)

  real(r_size), allocatable :: x_da(:), v_da(:) ! if EnKF: ensemble mean
  real(r_size), allocatable :: x_obs(:)
  
  ! --- Matrix(element 1:x, 2:v)
  ! +++ default setting
  ! Pf = (  Pxx: 1.0  Pxv: 0.0
  !         Pyx: 0.0  Pyy: 1.0  )
  
  real(r_size) :: M(2,2)   ! state transient matrix
  real(r_size) :: Pf(2,2)  ! Forecast error convariance matrix (in KF, EnKF)
  real(r_size) :: Pa(2,2)  ! Analysis error convariance matrix
  real(r_size) :: R(1,1)   ! Observation error convariance matrix
  real(r_size) :: Kg(2,1)  ! Kalman gain
  real(r_size) :: H(1,2)   ! Observation operator
  
  !  >> -----------------------------------------------------------------
  ! +++ 4Dvar_Ajoint model : Cost function and ajoint variable
  real(r_size) :: B(2,2)   ! Background error convariance matrix
  real(r_size) :: J        ! Cost function
  real(r_size) :: Jold     ! Cost function in a previous iteration
  real(r_size) :: x_b, v_b ! First guess of initial value

  real(r_size), allocatable :: adx(:), adv(:)
  real(r_size), allocatable :: x_save(:), v_save(:)
  real(r_size), allocatable :: x_tmp(:), v_tmp(:)
  real(r_size), parameter   :: alpx = 0.02         ! Coefficient for minization
  real(r_size), parameter   :: alpv = 0.02
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
  real(r_size) :: Ptmp(2,2)
  real(r_size) :: noise1, noise2, Gnoise ! Gaussian noise

  !======================================================================
  ! Data assimilation
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------
  real(r_size) :: x_tinit, v_tinit
  real(r_size) :: x_sinit, v_sinit
  real(r_size) :: Pf_init(4), B_init(4)
  real(r_size) :: R_init
  real(r_size) :: Kg_init(2)
  real(r_size) :: H_init(2)
  
  namelist /set_parm/ nt_asm, nt_prd, obs_interval
  namelist /da_setting/ da_method
  namelist /ensemble_size/ mems
  namelist /initial_osc/ x_tinit, v_tinit, x_sinit, v_sinit
  namelist /initial_que/ Pf_init, B_init, R_init, Kg_init, H_init
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

  allocate(x_true(0:nt_asm+nt_prd), v_true(0:nt_asm+nt_prd))
  allocate(x_sim(0:nt_asm+nt_prd), v_sim(0:nt_asm+nt_prd))
  allocate(x_da_m(0:nt_asm, mems), v_da_m(0:nt_asm, mems))
  allocate(x_prtb(mems), v_prtb(mems))
  allocate(x_da(0:nt_asm+nt_prd), v_da(0:nt_asm+nt_prd))
  allocate(x_tmp(0:nt_asm), v_tmp(0:nt_asm))
  allocate(adx(0:nt_asm), adv(0:nt_asm))
  allocate(x_save(0:nt_asm), v_save(0:nt_asm))
  allocate(x_obs(0:nt_asm/obs_interval))
  allocate(obs_chr(0:nt_asm+nt_prd))

  !----------------------------------------------------------------------
  ! +++ initial setting
  
  x_true(0) = x_tinit; v_true(0) = v_tinit
  x_sim(0)  = x_sinit; v_sim(0)  = v_sinit
  
  Pf(1,1) = Pf_init(1); Pf(1,2)=Pf_init(2)
  Pf(2,1) = Pf_init(3); Pf(2,2)=Pf_init(4)
  B(1,1)  = B_init(1);  B(1,2)=B_init(2)
  B(2,1)  = B_init(3);  B(2,2)=B_init(4)
  
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
  if ( da_method == 'KF' ) then
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
        write(6,*) 'Pf(1,1) / (R(1,1) + Pf(1,1))'
        write(6,*) '=', Pf(1,1), '/', R(1,1) , '+', Pf(1,1)
        write(6,*) 'Kg(1,1) = ', Kg(1,1)
        Kg(1,1) = Pf(1,1) / (R(1,1) + Pf(1,1))
        write(6,*) 'Pf(2,1) / (R(1,1) + Pf(1,1))'
        write(6,*) '=', Pf(2,1), '/', R(1,1) , '+', Pf(1,1)
        write(6,*) 'Kg(2,1) = ', Kg(1,1)
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
    
  else if ( da_method == 'EnKF' ) then
    ! making ensemble intial score
    do imem = 1, mems
      ! Generate Gaussian Noise (Gnoise) from uniform random number
      ! based on Box-Muller method
      call random_number(noise1)
      call random_number(noise2)
      Gnoise = sqrt(Pa(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
      x_da_m(0,imem) = x_sim(0) + Gnoise ! perturbation

      call random_number(noise1)
      call random_number(noise2)
      Gnoise = sqrt(Pa(2,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
      v_da_m(0,imem) = v_sim(0) + Gnoise ! perturbation
    end do

    x_da(0) = sum(x_da_m(0,1:mems))/mems
    v_da(0) = sum(v_da_m(0,1:mems))/mems

    do it = 1, nt_asm
      ! 4.1: Time integration
      do imem = 1, mems
        x_da_m(it, imem) = x_da_m(it-1, imem) + dt*v_da_m(it-1, imem)
        v_da_m(it, imem) = -(k * dt / mass) * x_da_m(it-1, imem) + (1.0d0 - dump * dt / mass ) * v_da_m(it-1, imem)
      end do

      if(mod(it,obs_interval) == 0) then
        x_da(it) = sum(x_da_m(it,1:mems))/mems
        v_da(it) = sum(v_da_m(it,1:mems))/mems
        Pf = 0.0d0
        do imem = 1, mems
          x_prtb(imem) = x_da_m(it, imem) - x_da(it)
          v_prtb(imem) = v_da_m(it, imem) - v_da(it)
          Pf(1,1) = Pf(1,1) + x_prtb(imem)**2/(mems-1)
          Pf(1,2) = Pf(1,2) + x_prtb(imem)*v_prtb(imem)/(mems-1)
          Pf(2,1) = Pf(2,1) + v_prtb(imem)*x_prtb(imem)/(mems-1)
          Pf(2,2) = Pf(2,2) + v_prtb(imem)**2/(mems-1)
        end do
        ! Section 4-2-2: Kalman gain: Weighting of model result and obs.
        ! (Note) Obsevation only in x --> component 1(x) only
        !        In this case, inverse matrix --> scalar inverse
        Kg(1,1) = Pf(1,1)/(R(1,1)+Pf(1,1))
        Kg(2,1) = Pf(2,1)/(R(1,1)+Pf(1,1))
        ! Section 4-2-3: calculate innovation and correction
        do imem=1, mems
          ! Generate Gaussian Noise (gnoise) from uniform random number 
          call random_number(noise1)
          call random_number(noise2)
          Gnoise = sqrt(R(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
          
          x_innov = x_obs(it/obs_interval)+Gnoise-x_da_m(it,imem)
          x_da_m(it,imem) = x_da_m(it,imem)+Kg(1,1)*x_innov
          v_da_m(it,imem) = v_da_m(it,imem)+Kg(2,1)*x_innov
        end do
        ! Section 4-2-4: analysis error covariance matrix
        Pa=Pf-matmul(matmul(Kg,H),Pf)
        Pf=Pa
        write(*,'(A,F8.2,A,5F10.3)') "time=",dt*it, ", Pa=", Pa
      end if
      x_da(it)=sum(x_da_m(it,1:mems))/mems
      v_da(it)=sum(v_da_m(it,1:mems))/mems
    end do

  else if ( da_method == 'Ajoint' ) then
    x_tmp(0) = x_sim(0)
    v_tmp(0) = v_sim(0)
    x_b = x_tmp(0) ! set First guess value
    v_b = v_tmp(0)

    do iter = 1, iter_max
      ! Section 4-1: initiallization for ajoint model run
      J = 0.0d0
      adx(0:nt_asm) = 0.0d0
      adv(0:nt_asm) = 0.0d0
      ! Section 4-2: Forward model run
      do it = 1, nt_asm
        x_tmp(it) = x_tmp(it-1) + dt*v_tmp(it-1)
        v_tmp(it) = -(k * dt / mass) * x_tmp(it-1) + (1.0d0 - dump * dt / mass ) * v_tmp(it-1)
        if (mod(it, obs_interval) == 0) then
          ! Calculate Cost function
          J = J + 0.5 * (x_obs(it / obs_interval) - x_tmp(it))**2 / R(1, 1)
        end if
      end do
      ! Section 4-3: ajoint model run
      do it = nt_asm, 1, -1
        if (mod(it, obs_interval) ==0) then
          ! Calculate misfit and chanfge ajoint variable
          adx(it) = adx(it) + (x_tmp(it) - x_obs(it/obs_interval))/R(1,1)
        end if
        adx(it-1) = adx(it) -k*dt / mass*adv(it)
        adv(it-1) = dt * adx(it) + (1.0d0 - dump*dt / mass) * adv(it)
      end do

      ! Section 4-4: Consider background covariance
      J = J + 0.5d0 * (x_tmp(0) - x_b)**2/ B(1,1) &
      &     + 0.5d0 * (v_tmp(0) - v_b)**2/ B(2,2)
      adx(0) = adx(0) + (x_tmp(0) -x_b) / B(1,1)
      adv(0) = adv(0) + (v_tmp(0) -v_b) / B(1,1)
      ! Section 4-5: Check the end of iteration
      if ((iter > 1) .and. (Jold < J)) then
        x_da(0:nt_asm) = x_tmp(0:nt_asm)
        v_da(0:nt_asm) = v_tmp(0:nt_asm)

        write(6,*)
        write(6,*) 'Cost function increases from the previous iteration.'
        write(6,*) 'Replace DA results with those in previous iteration and exit DA.'
        write(6,'(A,1X,i3,1X,A,3(F9.3,1X))') &
        & 'iteration =', iter, '; ajoint x,v =', adx(0), adv(0)
        write(6,'(2X,A,F9.3,1X,A,1X,3(F7.3,2X))') &
        & 'J =', J, '; x(0),v(0) =>', x_da(0), v_da(0)
        write(6,*)
        exit
      
      else if ((iter > 1) .and. (cond_iter > (Jold-J)/ Jold)) then
        x_da(0:nt_asm) = x_tmp(0:nt_asm)
        v_da(0:nt_asm) = v_tmp(0:nt_asm)

        write(6,*)
        write(6,*) 'Differences between J and Jold become small => exit DA.'
        write(6,'(A,1X,i3,1X,A,3(F9.3,1X))') &
        & 'iteration =', iter, '; ajoint x,v =', adx(0), adv(0)
        write(6,'(2X,A,F9.3,1X,A,1X,3(F7.3,2X))') &
        & 'J =', J, '; x(0),v(0) =>', x_da(0), v_da(0)
        write(6,*)
        exit
      
      else if (iter == iter_max) then
        x_da(0:nt_asm) = x_tmp(0:nt_asm)
        v_da(0:nt_asm) = v_tmp(0:nt_asm)

        write(6,*)
        write(6,*) 'Maximum number of iteration reached'
        write(6,'(A,1X,i3,1X,A,3(F9.3,1X))') &
        & 'iteration =', iter, '; ajoint x,v =', adx(0), adv(0)
        write(6,'(2X,A,F9.3,1X,A,1X,3(F7.3,2X))') &
        & 'J =', J, '; x(0),v(0) =>', x_da(0), v_da(0)
        write(6,*)
        exit
      end if

      ! Section 4-6: save values
      Jold = J

      ! Section 4-7: Update of x and v
      x_tmp(0) = x_tmp(0) - alpx*adx(0)
      v_tmp(0) = v_tmp(0) - alpv*adv(0)

      ! OUTPUT
      write(6,'(A,1X,i3,1X,A,1X,2F9.3,1X,A,1X,F9.3,1X,A,1X,3(F7.3,2X))') &
      & 'iteration =  ',   iter,                  &
      & '; x(0), v(0) = ', x_tmp(0), v_tmp(0),    &
      & '; J = ', J,                              &
      & '; ajoint x,v = ', adx(0), adv(0)

    end do
  end if
  
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

end program oscillation
