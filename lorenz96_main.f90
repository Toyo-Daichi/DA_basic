! Created on 2020.5.2 ~
! @author: Toyo_Daichi

program lorenz96_main

  use common
  use common_mtx
  use lorenz96_prm
  use lorenz96_cal
  
  implicit none

  ! --- settign exp.
  character(8)  :: tool, ts_check
  character(8)  :: da_method
  character(8)  :: enkf_method
  character(12) :: intg_method

  ! --- setting variance dimension
  ! *** x(timescale, nx)
  real(r_size), allocatable :: x_true(:,:)
  real(r_size), allocatable :: x_anl(:,:)
  real(r_size), allocatable :: x_sim(:,:)
  real(r_size), allocatable :: x_obs(:,:)
  ! *** x_mem(timescale, nx, member), x_prtb(nx, member)
  real(r_size), allocatable :: x_anl_m(:,:,:)
  real(r_size), allocatable :: x_prtb(:,:)
  
  ! --- For spinup
  ! *** x(timescale, nx)
  real(r_size), allocatable :: x_out(:,:)
  real(r_size), allocatable :: x_init(:)

  ! --- Temporal state vector
  real(r_size), allocatable :: xt_vec(:,:)
  real(r_size), allocatable :: yt_vec(:,:)
  
  ! *** Various parameters
  real(r_size), parameter   :: size_noise_obs = 1.0d0
  real(r_size)              :: gnoise, alpha
  real(r_dble)              :: delta
  integer                   :: mems
  
  ! --- Data assimilation exp.
  ! +++ Matrix setting
  real(r_size), allocatable :: JM(:,:)  ! state transient matrix (=Jacobean　matrix)
  real(r_size), allocatable :: Pf(:,:)  ! Forecast error convariance matrix (in KF, EnKF)
  real(r_size), allocatable :: Pa(:,:)  ! Analysis error convariance matrix
  real(r_size), allocatable ::  I(:,:)  ! eye_matrix
  real(r_size), allocatable ::  R(:,:)  ! Observation error convariance matrix
  real(r_size), allocatable :: Kg(:,:)  ! Kalman gain
  real(r_size), allocatable :: Kh(:,:)  ! Kalman gain for pertuvation
  real(r_size), allocatable ::  H(:,:)  ! Observation operator
  
  real(r_size), allocatable :: Ef(:,:), Ea(:,:) ! Ensemble pertubation
  real(r_size), allocatable :: obs_prtbmtx(:,:) ! Ensemble pertubation

  real(r_size), allocatable :: hx(:,:)
  real(r_size), allocatable :: hdxf(:,:)
  real(r_size), allocatable :: work(:)

  real(r_size), allocatable :: obs_mtx(:,:)
  real(r_size), allocatable :: obs_inv(:,:)

  ! --- for EnSRF
  real(r_size), allocatable :: obs_srf_mtx_v1(:,:), obs_srf_mtx_v2(:,:)
  real(r_size), allocatable :: obs_srf_inv_v1(:,:), obs_srf_inv_v2(:,:)
  real(r_size) :: sqrt_mems

  ! --- Output control
  logical, save         :: opt_veach = .false.
  logical, save         :: da_veach  = .false.
  character(256)        :: initial_true_file
  character(256)        :: initial_sim_file
  character(256)        :: output_true_file
  character(256)        :: output_anl_file
  character(256)        :: output_sim_file
  character(256)        :: output_obs_file
  character(256)        :: output_errcov_file
  
  ! --- Working variable
  character(1086) :: linebuf
  character(36)   :: cfmt, cfmt_obs
  character(4)    :: cfmt_num, cfmt_obsnum
  integer         :: spinup_period, normal_period
  integer         :: kt_oneday
  integer         :: ix, iy, it, il, imem, ierr, lda, ipiv, lwork

  integer, parameter :: one_loop = 1

  !======================================================================
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------
  
  namelist /set_parm/ nx, dt, force, oneday
  namelist /set_exp/ tool, ts_check, intg_method
  namelist /set_da_exp/ da_veach, mems, da_method, alpha
  namelist /enkf_setting/ mems, enkf_method
  namelist /set_period/ spinup_period, normal_period
  namelist /set_mobs/ obs_xintv, obs_tintv
  namelist /output/ initial_true_file, initial_sim_file, &
    output_true_file, output_anl_file, output_sim_file, output_obs_file, output_errcov_file, opt_veach
  
  read(5, nml=set_parm, iostat=ierr)
  read(5, nml=set_exp, iostat=ierr)
  read(5, nml=set_da_exp, iostat=ierr)
  read(5, nml=enkf_setting, iostat = ierr)
  read(5, nml=set_period, iostat=ierr)
  read(5, nml=set_mobs, iostat=ierr)
  read(5, nml=output, iostat=ierr)
  ! name list io check
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
  
  kt_oneday = int(oneday/dt) ! Unit change for 1day
  if ( trim(tool) == 'spinup' ) then 
    allocate(x_true(1,nx), x_out(1, nx))
    allocate(x_init(1:nx))
    da_veach = .false.
  else if ( trim(tool) == 'normal' ) then
    allocate(x_true(0:kt_oneday*normal_period, nx))
    allocate(x_init(nx))
    
    ! --- OBS setting
    ny = int(nx/obs_xintv)
    obs_time = int((kt_oneday*normal_period)/obs_tintv)
    lda = ny
    lwork = ny
    allocate(x_obs(obs_time,ny))
    allocate(work(ny))
    
    ! --- Temporal state vector
    allocate(xt_vec(nx,1))
    allocate(yt_vec(ny,1))
    ! --- EnSRF setting
    sqrt_mems = sqrt(float(mems-1))
    
    !----------------------------------------------------------------------
    ! +++ Data assim. set 
    !----------------------------------------------------------------------
    if ( da_veach ) then
      allocate(x_anl(0:kt_oneday*normal_period, nx))
      allocate(x_sim(0:kt_oneday*normal_period, nx))
      if ( da_method == 'EnKF' ) then
        allocate(x_anl_m(0:kt_oneday*normal_period, nx, mems))
        allocate(x_prtb(nx, mems))
      end if

      ! covariance matrix set.
      allocate(Pf(1:nx, 1:nx))
      allocate(Pa(1:nx, 1:nx))
      allocate(JM(1:nx, 1:nx))
      allocate( I(1:nx, 1:nx))
      allocate(Kg(1:nx, 1:ny))
      allocate( H(1:ny, 1:nx))
      allocate( R(1:ny, 1:ny))

      allocate(Ef(nx,mems), Ea(nx,mems), obs_prtbmtx(ny, mems))

      ! for kalmangain inverse calculate.
      allocate(obs_mtx(1:ny, 1:ny), obs_inv(1:ny, 1:ny))
      allocate(obs_srf_mtx_v1(1:ny, 1:ny), obs_srf_mtx_v2(1:ny, 1:ny))
      allocate(obs_srf_inv_v1(1:ny, 1:ny), obs_srf_inv_v2(1:ny, 1:ny))
      
      allocate(hx(ny,1))
      allocate(hdxf(ny,1))
      
    end if
  end if

  !======================================================================
  !
  ! --- Sec.2 lorenz96 calculation
  !----------------------------------------------------------------------
  ! +++ making true run.
  !----------------------------------------------------------------------

  write(6,*) 'Exp. setting            :: ', tool
  write(6,*) 'Integral method         :: ', intg_method
  
  if ( trim(tool) == 'spinup' ) then
    call com_randn(nx, x_true(1,:))
    x_true(1,:) = x_true(1,:)*5.0d0
    call ting_rk4(kt_oneday*spinup_period, x_true(1,:), x_out(1,:))
    
  else if ( trim(tool) == 'normal' ) then
    open(2, file=trim(initial_true_file), form='formatted', status='old')
      read(2,*) x_init
    close(2)
    x_true(0, :) = x_init
    do it = 1, kt_oneday*normal_period
      call ting_rk4(one_loop, x_true(it-1,:), x_true(it,:))
      
      !-------------------------------------------------------------------
      ! +++ making obs score
      !-------------------------------------------------------------------
      if ( mod (it,obs_tintv)==0 ) then
        do ix=1,nx
          if( mod(ix,obs_xintv)==0 ) then
            call gaussian_noise(size_noise_obs, gnoise)
            x_obs(it/obs_tintv, ix/obs_xintv) = x_true(it, ix) + gnoise
          endif
        enddo
      endif
    end do
    close(2)
  end if
  
  !======================================================================
  !
  ! --- Sec.3 Data Assimlation
  DataAssim_exp_set :&
  if ( da_veach ) then
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Data assimilation exp. start '
    write(6,*) ' >> Data assimilation method  :: ', da_method
    
    ! making identity matrix
    Pf = 0.d0; Pa = 0.d0; I = 0.d0
    Kg = 0.d0;  H = 0.d0; R = 0.d0

    forall ( il=1:nx ) Pf(il, il) = 1.0d0 
    forall ( il=1:nx ) Pa(il, il) = 1.0d0
    forall ( il=1:nx )  I(il, il) = 1.0d0
    forall ( il=1:ny )  R(il, il) = size_noise_obs
    forall ( il=1:ny )  H(il, il) = 1.0d0 ! not regular matrix
    
    !-------------------------------------------------------------------
    ! +++ making sim score
    !-------------------------------------------------------------------
    open(2, file=trim(initial_sim_file), form='formatted', status='old')
    read(2,*) x_init
    close(2)
    x_sim(0,:) = x_init
    do it = 1, kt_oneday*normal_period
      call ting_rk4(one_loop, x_sim(it-1,:), x_sim(it,:))
    end do
    
    !-------------------------------------------------------------------
    ! +++ Data assimilation
    !
    ! 1st. Extend Kalman Filter
    ! *** Key: Pf = JM*Pa*JM^T
    !-------------------------------------------------------------------
    
    DataAssim_exp_kind :&
    if ( trim(da_method) == 'KF' ) then
      ! --- initialize
      x_anl(0,:) = x_sim(0,:)
      
      data_assim_KF_loop :&
      do it = 1, kt_oneday*normal_period
        write(6,*) 'Data assim. time step: ', it
        write(6,*) ''
        call ting_rk4(one_loop, x_anl(it-1,:), x_anl(it,:))
        
        if ( mod(it, obs_tintv)==0 ) then
          write(6,*) '  TRUTH    = ', x_true(it,1:5), '...'
          write(6,*) '  PREDICT  = ', x_sim(it,1:5), '...'
          write(6,*) '  OBSERVE  = ', x_obs(it,1:5), '...'
          write(6,*) '  ANALYSIS (BEFORE) = ', x_anl(it,1:5), '...'
          write(6,*) ''
          
          !-------------------------------------------------------------------
          ! +++ Making Jacobian matrix.
          ! >> Please note that time changes are included...
          !-------------------------------------------------------------------
          delta = dt*obs_tintv
          JM = 0.0d0
          JM(1,1)     = 1.0d0 - delta
          JM(1,2)     = x_anl(it-1,nx)*delta
          JM(1,nx-1)  = -x_anl(it-1,nx)*delta
          JM(1,nx)    = (x_anl(it-1,2)-x_anl(it-1,nx-1))*delta
          
          JM(2,1)     = (x_anl(it-1,3)-x_anl(it-1,nx))*delta
          JM(2,2)     = 1.0d0 - delta
          JM(2,3)     = x_anl(it-1,1)*delta
          JM(2,nx)    = -x_anl(it-1,1)*delta 
          
          do il = 3, nx-1
            JM(il,il-2) = -x_anl(it-1,il-1)*delta
            JM(il,il-1) = (x_anl(it-1,il+1)-x_anl(it-1,il-2))*delta
            JM(il,il)   = 1.0d0 - delta
            JM(il,il+1) = x_anl(it-1,il-1)*delta
          end do
          
          JM(nx,1) = x_anl(it-1,nx-1)*delta
          JM(nx,nx-2) = -x_anl(it-1,nx-1)*delta
          JM(nx,nx-1) = (x_anl(it-1,1)-x_anl(it-1,nx-2))*delta
          JM(nx,nx) = 1.0d0 - delta
          
          Pf = matmul(matmul(JM, Pa), transpose(JM))*(1.0d0 + alpha)
          write(6,*) '  PREDICTION ERROR COVARIANCE on present step'
          call confirm_matrix(Pf, nx, nx)
          write(6,*) ''
          
          !-----------------------------------------------------------
          ! +++ making inverse matrix
          !-----------------------------------------------------------
          obs_mtx = matmul(matmul(H, Pf), transpose(H)) + R
          obs_inv = obs_mtx
          call dgetrf(ny, ny, obs_inv, lda, ipiv, ierr)
          call dgetri(ny, obs_inv, lda, ipiv, work, lwork, ierr)
          
          !-----------------------------------------------------------
          ! +++ inflation mode
          !-----------------------------------------------------------
          write(6,*) ' KALMAN GAIN WEIGHTING MATRIX '
          Kg = matmul(matmul(Pf, transpose(H)), obs_inv)
          call confirm_matrix(Kg, nx, ny)
          write(6,*) ''

          x_anl(it,:) = x_anl(it,:) + matmul(Kg, (x_obs(it/obs_tintv,:) - matmul(H, x_anl(it,:))))
          write(6,*) '  ANALYSIS (AFTER) = ', x_anl(it,1:5), '...'
          
          Pa = matmul((I - matmul(Kg, H)), Pf)
          write(6,*) '  ANALYSIS ERROR COVARIANCE on present step.'
          call confirm_matrix(Pa, nx, nx)
          write(6,*) ''

          call output_errcov(nx, it/obs_tintv, kt_oneday*normal_period/obs_tintv, Pa)
        end if
        
      end do &
      data_assim_KF_loop
      
    else if ( trim(da_method) == 'EnKF' ) then
      !-----------------------------------------------------------
      ! +++ making ensemble initial score
      !-----------------------------------------------------------
      do ix = 1, nx
        do imem = 1, mems
          call gaussian_noise(size_noise_obs, gnoise)
          x_anl_m(0, ix, imem) = x_sim(0, ix) + gnoise
        end do
      end do

      do ix = 1, nx
        x_anl(0, ix) = sum(x_anl_m(0,ix,1:mems))/mems
      end do
      
      !-----------------------------------------------------------
      ! +++ assimilation loop start
      !-----------------------------------------------------------
      data_assim_EnKF_loop :&
      do it = 1, kt_oneday*normal_period
        write(6,*) 'Data assim. time step: ', it

        write(6,*) ' CHECK ENSEMBLE PART FOR UPDATE (BEFORE) on ', it-1
        write(6,*) x_anl_m(it-1,1,1:5)
        
        do imem = 1, mems
          call ting_rk4(one_loop, x_anl_m(it-1,:,imem), x_anl_m(it,:,imem))
        end do

        
        if ( mod(it, obs_tintv) == 0 .and. it /= 0) then
          write(6,*) '  TRUTH    = ', x_true(it,1:5), '...'
          write(6,*) '  PREDICT  = ', x_sim(it,1:5), '...'
          write(6,*) '  OBSERVE  = ', x_obs(it,1:5), '...'
          do ix = 1, nx
            x_anl(it, ix) = sum(x_anl_m(it,ix,1:mems))/mems
          end do
          write(6,*) '  ANALYSIS (BEFORE) = ', x_anl(it,1:5), '...'
          write(6,*) ''
          Pf = 0.0d0
          
          !=======================================================
          !
          ! Pf = Ef*Ef^T / (mems-1) 
          !------------------------------------------------------- 
          ! +++ making error covariance matix by pertubation
          ! >> this step is not include Ef sqrt(mems-1)
          !    Check Pf = Pf/(mems-1)
          !------------------------------------------------------- 
          do imem = 1, mems
            x_prtb(:,imem) = x_anl_m(it,:,imem) - x_anl(it,:)
            Ef(:,imem) = x_prtb(:,imem)
          end do

          write(6,*) ' PREDICTION ENSEMBLE VECTOR '
          call confirm_matrix(Ef, nx, mems)
          write(6,*) ''

          write(6,*) '  PREDICTION ERROR COVARIANCE on present step'
          Pf = matmul(Ef, transpose(Ef))
          Pf = Pf/(mems-1)
          call confirm_matrix(Pf, nx, nx)
          write(6,*) ''

          !-----------------------------------------------------------
          ! +++ making inverse matrix
          !-----------------------------------------------------------
          obs_mtx = matmul(matmul(H, Pf), transpose(H)) + R
          obs_inv = obs_mtx
          call dgetrf(ny, ny, obs_inv, lda, ipiv, ierr)
          call dgetri(ny, obs_inv, lda, ipiv, work, lwork, ierr)
          
          !-----------------------------------------------------------
          ! +++ adaptive inflation mode
          !-----------------------------------------------------------
          Pf = Pf*(1.0d0 + alpha)
          write(6,*) ' KALMAN GAIN WEIGHTING MATRIX '
          Kg = matmul(matmul(Pf, transpose(H)), obs_inv)
          call confirm_matrix(Kg, nx, ny)
          write(6,*) ''
          
          da_enkf_method :&
          if ( trim(enkf_method) == 'PO' ) then
            !------------------------------------------------------- 
            ! >> EnKF
            ! +++ Pertuturbed observation method (PO)
            ! (Phase. plus OBS pertubation)
            ! x_anl = x_sim + Kg(y +"e" -H*x_sim)
            !------------------------------------------------------- 
            
            ! +++ making obs pertubation err.
            do iy = 1, ny
              do imem = 1, mems
                call gaussian_noise(size_noise_obs, gnoise)
                obs_prtbmtx(iy, imem) = x_obs(it/obs_tintv, iy) + gnoise 
              end do
            end do

            do imem = 1, mems
              ! for matrix calculation
              xt_vec(:,1) = x_anl_m(it,:,imem)
              yt_vec(:,1) = obs_prtbmtx(:,imem)

              ! simple KF formula
              hx = matmul(H, xt_vec)
              hdxf = yt_vec - hx
              xt_vec = xt_vec + matmul(Kg, hdxf)

              ! fix matrix -> each score
              x_anl_m(it,:,imem) =  xt_vec(:,1)
            end do

            do ix = 1, nx
              x_anl(it, ix) = sum(x_anl_m(it, ix, :))/mems
            end do
            write(6,*) '  ANALYSIS (AFTER) = ', x_anl(it,1:5), '...'
            write(6,*) ''
            write(6,*) ' CHECK ENSEMBLE PART FOR UPDATE (AFTER) on ', it 
            write(6,*) x_anl_m(it,1,1:5)
            write(6,*) ''
          
          else if ( trim(enkf_method) == 'SRF' ) then
            !----------------------------------------------------------- 
            ! >> EnKF (-> Serial EnSRF)
            ! +++ Squared root fileter method (SRF)
            ! Ea = (I - K^H)*Ef
            ! K^ = 
            ! Pf H^T [(H Pf H^T + R)^-1/2]^T [(H Pf H^T + R)^-1/2]^-1
            !----------------------------------------------------------- 

            ! +++ (1) Average step 
            x_anl(it,:) = x_anl(it,:) + matmul(Kg, (x_obs(it/obs_tintv,:) - matmul(H, x_anl(it,:))))
            write(6,*) '  ANALYSIS (AVERAGE STEP) = ', x_anl(it,1:5), '...'
            write(6,*) ''

            ! +++ (2) Pertubation step
            obs_srf_mtx_v1 = matmul(H, matmul(Pf, transpose(H))) + R
            call mtx_sqrt(ny, obs_srf_mtx_v1, obs_srf_inv_v1)
            
            call dgetrf(ny, ny, obs_srf_inv_v1, lda, ipiv, ierr)
            call dgetri(ny, obs_srf_inv_v1, lda, ipiv, work, lwork, ierr)
            
            obs_srf_mtx_v2 = obs_srf_mtx_v1 + R
            call mtx_sqrt(ny, obs_srf_mtx_v2, obs_srf_inv_v2)
            obs_srf_inv_v2 = obs_srf_mtx_v2
            
            call dgetrf(ny, ny, obs_srf_inv_v2, lda, ipiv, ierr)
            call dgetri(ny, obs_srf_inv_v2, lda, ipiv, work, lwork, ierr)

            Kh = matmul(matmul(Pf, transpose(H)), transpose(obs_srf_inv_v1))
            Kh = matmul(Kh, obs_srf_inv_v2)
          
            write(6,*) ' KALMAN GAIN HAT WEIGHTING MATRIX '
            call confirm_matrix(Kh, nx, ny)
            write(6,*) ''
          
            write(6,*) ' ANALYSIS ENSEMBLE VECTOR '
            Ef = Ef/sqrt_mems
            Ea = matmul(I - matmul(Kh,H), Ef)
            call confirm_matrix(Ea, nx, mems)

            write(6,*) ' && AVERAGE PERTUBATION   '
            write(6,*) sum(Ea(1,:))/mems, sum(Ea(2,:))/mems, sum(Ea(3,:))/mems, &
                       sum(Ea(4,:))/mems, sum(Ea(5,:))/mems, '...'
            write(6,*) ''

            ! +++ for next step, initiallize
            do imem = 1, mems
              x_anl_m(it, :, imem) = x_anl(it, :) + Ea(:,imem)
            end do

            obs_srf_mtx_v1 = 0.0d0; obs_srf_mtx_v2 = 0.0d0
            obs_srf_inv_v1 = 0.0d0; obs_srf_inv_v2 = 0.0d0

            ! +++ Union step (ave. + pertb)
            do ix = 1, nx
              x_anl(it, ix) = x_anl(it, ix) + sum(Ea(ix, :))/mems
            end do
            write(6,*) '  ANALYSIS (PRTB STEP) = ', x_anl(it,1:5), '...'
            write(6,*) ''

          end if &
          da_enkf_method
        end if

      end do &
      data_assim_EnKF_loop
      
    end if &
    DataAssim_exp_kind

  end if &
  DataAssim_exp_set
  
  !======================================================================
  !
  ! --- Sec.* Writing OUTPUT
  if ( opt_veach ) then
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Check Writing output system,  '
    if (nx <= 100 ) then 
      cfmt = '(xx(F12.7, ","), F12.7)'
      write(cfmt_num,"(I2)") nx-1
      cfmt(2:3) = cfmt_num
    else if (nx >= 100 ) then
      cfmt = 'xxx(F12.7, ","), F12.7)'
      write(cfmt_num, "(I3)") nx-1
      cfmt(2:4) = cfmt_num
    end if
    
    ! select open file
    if ( trim(tool) == 'spinup' ) then
      if ( trim(ts_check) == 'true' ) then
        open(2, file=trim(initial_true_file), form='formatted', status='replace')
      else if (trim(ts_check) == 'sim' ) then
        open(2, file=trim(initial_sim_file), form='formatted', status='replace')
      end if
      
      write(linebuf, trim(cfmt)) x_out(0,:)
      call del_spaces(linebuf)
      write(2,'(a)') linebuf
      write(6,*) ' && Successfuly output !!!        ' 
      
    else if ( trim(tool) == 'normal' ) then
      open(21, file=trim(output_true_file), form='formatted', status='replace')
      open(22, file=trim(output_anl_file), form='formatted', status='replace')
      open(23, file=trim(output_sim_file), form='formatted', status='replace')
      do it = 0, kt_oneday*normal_period
        
        write(linebuf, trim(cfmt)) x_true(it,:)
        call del_spaces(linebuf)
        write(21,'(a)') linebuf
        write(linebuf, trim(cfmt)) x_anl(it,:)
        call del_spaces(linebuf)
        write(22,'(a)') linebuf
        write(linebuf, trim(cfmt)) x_sim(it,:)
        call del_spaces(linebuf)
        write(23,'(a)') linebuf
      end do
      close(21)
      close(22)
      close(23)
      
      open(24, file=trim(output_obs_file), form='formatted', status='replace')
        do it = 1, kt_oneday*normal_period/obs_tintv
          cfmt_obs = '(xx(F15.7, ","), F15.7)'
          write(cfmt_obsnum,"(I2)") ny-1
          cfmt_obs(2:3) = cfmt_obsnum
          write(linebuf, trim(cfmt_obs)) x_obs(it, :)
          call del_spaces(linebuf)
          write(24,'(a)') linebuf
        end do
      close(24)
      write(6,*) ' && Successfuly output !!!        '  
    endif
    
  else if ( .not. opt_veach ) then
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Check calculation system,  '
    write(6,*) ' && Successfuly calculate !!!  '  
    
  end if
  
  !======================================================================
  !
  ! +++ Tidy up
  ! deallocate(x_true, x_init)
  
contains
  
  subroutine gaussian_noise(    &
    size_gnoise,                & ! IN: observation error diag.
    gnoise                      & ! OUT
    )
    ! Generate Gaussian Noise (Gnoise) from uniform random number
    ! based on Box-Muller method

    implicit none
    real(kind=r_size),intent(in)  :: size_gnoise
    real(kind=r_size),intent(out) :: gnoise
  
    ! constant
    real(kind=r_dble),parameter :: pi=3.14159265358979
  
    ! working variable
    real(kind=r_dble) :: noise1,noise2
  
    call random_number(noise1)
    call random_number(noise2)
    gnoise=size_gnoise*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
  
  end subroutine gaussian_noise

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

  subroutine output_errcov(nx, it, last_step, errcov_mtx)
    implicit none
    integer, intent(in)       :: nx, it, last_step
    real(r_size), intent(in)  :: errcov_mtx(nx,nx)

    character(100000) :: linebuf
    character(36)     :: cfmt
    character(4)      :: cfmt_num

    logical, save   :: opt_veach = .false.

    cfmt = '(xxxx(F12.7, ","), F12.7)'
    write(cfmt_num,"(I4)") nx*nx -1
    cfmt(2:5) = cfmt_num

    if ( it == 1 ) then
      open(2, file=trim(output_errcov_file), status='replace')
        write(linebuf, trim(cfmt)) errcov_mtx(:,:)
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        write(6,*) '+++ err covariance matrix 1st. step'
        
      else if ( it /= 0 .and. it /= last_step) then
        write(6,*) '+++ err covariance matrix 2nd. step ~'
        write(linebuf, trim(cfmt)) errcov_mtx(:,:)
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        
      else if ( it == last_step ) then
        write(linebuf, trim(cfmt)) errcov_mtx(:,:)
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        write(6,*) '+++ err covariance matrix last step '
      close(2)
    end if
      
  end subroutine


end program lorenz96_main