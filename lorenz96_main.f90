! Created on 2020.5.2 ~
! @author: Toyo_Daichi

program lorenz96_main
  use common
  use lorenz96_prm
  use lorenz96_cal
  
  implicit none

  ! --- setting parameter
  ! *** x(timescale, nx)
  real(r_size), allocatable :: x_true(:,:)
  real(r_size), allocatable :: x_anl(:,:)
  real(r_size), allocatable :: x_sim(:,:)
  ! *** x_mem(timescale, nx, member), x_prtb(nx, member)
  real(r_size), allocatable :: x_anl_mem(:,:,:)
  real(r_size), allocatable :: x_prtb(:,:)

  ! --- For spinup && OBS info
  real(r_size), allocatable :: x_out(:,:)
  real(r_size), allocatable :: x_tmp(:)
  real(r_size), allocatable :: x_e(:)
  real(r_size), allocatable :: yt_obs(:,:)
  real(r_size), allocatable :: yt_obs_ens(:,:)
  ! ***
  real(r_size), parameter   :: size_noise_obs = 1.0d0
  integer                   :: ny, nt
  integer                   :: obs_xintv, obs_tintv
  integer                   :: mems

  ! --- settign exp.
  character(8)  :: tool, ts_check, da_method
  character(12) :: intg_method
  
  ! --- Data assimilation exp.
  ! *** Extend Kalamn Filter
  real(r_size), allocatable :: Pf(:,:)
  real(r_size), allocatable :: Pa(:,:)
  real(r_size), allocatable :: JM(:,:)
  real(r_size), allocatable ::  I(:,:)
  real(r_size), allocatable :: Kg(:,:)
  real(r_size), allocatable ::  H(:,:)
  real(r_size), allocatable ::  R(:,:)

  real(r_size), allocatable :: obs_matrix(:,:)
  real(r_size), allocatable :: obs_inv_matrix(:,:)
  real(r_size), allocatable :: eye_matrix(:,:)
  real(r_size), allocatable :: work(:)

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
  real(r_size)    :: gnoise, alpha
  real(r_dble)    :: delta
  integer         :: spinup_period, normal_period
  integer         :: kt_oneday
  integer         :: ix, it, il, imem, ierr, lda, ipiv, lwork

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
  namelist /set_period/ spinup_period, normal_period
  namelist /set_mobs/ obs_xintv, obs_tintv
  namelist /output/ initial_true_file, initial_sim_file, &
    output_true_file, output_anl_file, output_sim_file, output_obs_file, output_errcov_file, opt_veach
  
  read(5, nml=set_parm, iostat=ierr)
  read(5, nml=set_exp, iostat=ierr)
  read(5, nml=set_da_exp, iostat=ierr)
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
    allocate(x_tmp(nx))
    da_veach = .false.
  else if ( trim(tool) == 'normal' ) then
    allocate(x_true(0:kt_oneday*normal_period, nx))
    allocate(x_tmp(nx))
    
    ! Obs set.
    ny = int(nx/obs_xintv)
    nt = int((kt_oneday*normal_period)/obs_tintv)
    allocate(yt_obs(nt,ny))
    lda = ny; lwork = ny
    allocate(work(ny))
    
    !----------------------------------------------------------------------
    ! +++ Data assim. set 
    !----------------------------------------------------------------------
    if ( da_veach ) then
      allocate(x_anl(0:kt_oneday*normal_period, nx))
      allocate(x_sim(0:kt_oneday*normal_period, nx))
      if ( da_method == 'EnKF' ) then
        allocate(x_anl_mem(0:kt_oneday*normal_period, nx, mems))
        allocate(x_prtb(nx, mems))
        allocate(yt_obs_ens(ny, mems))
      end if
      ! covariance matrix set.
      allocate(x_e(1:nx))
      allocate(Pf(1:nx, 1:nx))
      allocate(Pa(1:nx, 1:nx))
      allocate(JM(1:nx, 1:nx))
      allocate( I(1:nx, 1:nx))
      allocate(Kg(1:nx, 1:ny))
      allocate( H(1:ny, 1:nx))
      ! for kalmangain inverse calculate.
      allocate(R(1:ny, 1:ny))
      allocate(eye_matrix(1:ny, 1:ny))
      allocate(obs_matrix(1:ny, 1:ny))
      allocate(obs_inv_matrix(1:ny, 1:ny))
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
      read(2,*) x_tmp
    close(2)
    x_true(0, :) = x_tmp
    do it = 1, kt_oneday*normal_period
      call ting_rk4(one_loop, x_true(it-1,:), x_true(it,:))
      
      !-------------------------------------------------------------------
      ! +++ making obs score
      !-------------------------------------------------------------------
      if ( mod (it,obs_tintv)==0 ) then
        do ix=1,nx
          if( mod(ix,obs_xintv)==0 ) then
            call gaussian_noise(size_noise_obs, gnoise)
            yt_obs(it/obs_tintv, ix/obs_xintv) = x_true(it, ix) + gnoise
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
    ! +++ making NoDA score
    !-------------------------------------------------------------------
    open(2, file=trim(initial_sim_file), form='formatted', status='old')
    read(2,*) x_tmp
    close(2)
    x_sim(0,:) = x_tmp
    do it = 1, kt_oneday*normal_period
      call ting_rk4(one_loop, x_sim(it-1,:), x_sim(it,:))
    end do
    
    !-------------------------------------------------------------------
    ! +++ Data assimilation
    DataAssim_exp_kind :&
    if ( trim(da_method) == 'KF' ) then
      ! --- initialize
      x_anl(0,:) = x_sim(0,:)
      
      data_assim_KF_loop :&
      do it = 1, kt_oneday*normal_period
        write(6,*) 'Data assim. time step: ', it
        call ting_rk4(one_loop, x_anl(it-1,:), x_anl(it,:))
        
        if ( mod(it, obs_tintv)==0 ) then
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

          call confirm_matrix(JM, nx, nx)
          Pf = matmul(matmul(JM, Pa), transpose(JM))*(1.0d0 + alpha)
          !-----------------------------------------------------------
          ! +++ making inverse matrix
          !-----------------------------------------------------------
          obs_matrix = matmul(matmul(H, Pf), transpose(H)) + R
          obs_inv_matrix = obs_matrix
          call dgetrf(ny, ny, obs_inv_matrix, lda, ipiv, ierr)
          call dgetri(ny, obs_inv_matrix, lda, ipiv, work, lwork, ierr)
          eye_matrix = matmul(obs_matrix, obs_inv_matrix)
          
          !-----------------------------------------------------------
          ! +++ adaptive inflation mode
          !-----------------------------------------------------------
          Kg = matmul(matmul(Pf, transpose(H)), obs_inv_matrix)
          x_anl(it,:) = x_anl(it,:) + matmul(Kg, (yt_obs(it/obs_tintv,:) - matmul(H, x_anl(it,:))))
          Pa = matmul((I - matmul(Kg, H)), Pf)
          call confirm_matrix(Pa, nx, nx)
          call write_errcov(nx, it/obs_tintv, kt_oneday*normal_period/obs_tintv, Pa)
        end if
      end do &
      data_assim_KF_loop
      
    else if ( trim(da_method) == 'EnKF' ) then
      !-----------------------------------------------------------
      ! +++ making ensemble initial score
      !-----------------------------------------------------------
      do imem = 1, mems
        call gaussian_noise(size_noise_obs, gnoise)
        x_anl_mem(0, :, imem) = x_sim(0, :) + gnoise
      end do
      do ix = 1, nx
        x_anl(0, ix) = sum(x_anl_mem(0,ix,1:mems))/mems
      end do
      
      !-----------------------------------------------------------
      ! +++ assimilation loop start
      data_assim_EnKF_loop :&
      do it = 1, kt_oneday*normal_period
        write(6,*) 'Data assim. time step: ', it
        do imem = 1, mems
          call ting_rk4(one_loop, x_anl_mem(it-1,:, imem), x_anl_mem(it,:, imem))
        end do
        
        if ( mod(it, obs_tintv)==0 .and. it /= 0) then
          Pf = 0.0d0
          do ix = 1, nx
            x_anl(it, ix) = sum(x_anl_mem(it,ix,1:mems))/mems
          end do
          do imem = 1, mems
            x_prtb(:, imem) = x_anl_mem(it, :, imem) - x_anl(it, :)
            !------------------------------------------------------- 
            ! +++ making ptbmtx
            !------------------------------------------------------- 
            forall ( il=1:nx ) Pf(il,il) = Pf(il,il) + x_prtb(il,imem)**2/(mems-1)
            do il = 1, nx-1
              Pf(il,il+1) = Pf(il,il+1) + x_prtb(il, imem)*x_prtb(il+1,imem)/(mems-1)
              Pf(il+1,il) = Pf(il,il+1)
            end do
            Pf(nx,1) = Pf(nx,1) + x_prtb(nx, imem)*x_prtb(1,imem)/(mems-1)
          end do
          
          !-----------------------------------------------------------
          ! +++ making inverse matrix
          !-----------------------------------------------------------
          obs_matrix = matmul(matmul(H, Pf), transpose(H)) + R
          obs_inv_matrix = obs_matrix
          call dgetrf(ny, ny, obs_inv_matrix, lda, ipiv, ierr)
          call dgetri(ny, obs_inv_matrix, lda, ipiv, work, lwork, ierr)
          eye_matrix = matmul(obs_matrix, obs_inv_matrix)
          !call confirm_matrix(Pf, nx, nx)
          
          !-----------------------------------------------------------
          ! +++ adaptive inflation mode
          !-----------------------------------------------------------
          Pf = Pf*(1.0d0 + 0.1d0)
          Kg = matmul(matmul(Pf, transpose(H)), obs_inv_matrix)
          
          !-----------------------------------------------------------
          ! +++ Pertuturbed observation method (PO)
          !-----------------------------------------------------------
          !do imem = 1, mems
          !  do iy = 1, ny
          !  call gaussian_noise(size_noise_obs, gnoise)
          !  yt_obs_ens(iy, imem) = yt_obs(it, iy) + gnoise
          !
          !nd do

          x_anl(it,:) = x_anl(it,:) + matmul(Kg, (yt_obs(it/obs_tintv,:) - matmul(H, x_anl(it,:))))
          Pa = Pf - matmul(matmul(Kg, H), Pf)
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
          write(linebuf, trim(cfmt_obs)) yt_obs(it, :)
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
  ! +++ tidy up
  
  deallocate(x_true, x_tmp)
  
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

  subroutine write_errcov(nx, it, last_step, error_covariance_matrix)
    implicit none
    integer, intent(in)       :: nx, it, last_step
    real(r_size), intent(in)  :: error_covariance_matrix(nx,nx)

    character(100000) :: linebuf
    character(36)     :: cfmt
    character(4)      :: cfmt_num

    logical, save   :: opt_veach = .false.

    cfmt = '(xxxx(F12.7, ","), F12.7)'
    write(cfmt_num,"(I4)") nx*nx -1
    cfmt(2:5) = cfmt_num

    if ( it == 1 ) then
      open(2, file=trim(output_errcov_file), status='replace')
      write(6,*) cfmt
        write(linebuf, trim(cfmt)) error_covariance_matrix(:,:)
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        write(6,*) '+++ err covariance matrix 1st. step'
        
      else if ( it /= 0 .and. it /= last_step) then
        write(6,*) '+++ err covariance matrix 2nd. step ~'
        write(linebuf, trim(cfmt)) error_covariance_matrix(:,:)
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        
      else if ( it == last_step ) then
        write(linebuf, trim(cfmt)) error_covariance_matrix(:,:)
        call del_spaces(linebuf)
        write(2, '(a)') trim(linebuf)
        write(6,*) '+++ err covariance matrix last step '
      close(2)
    end if
      
  end subroutine


end program lorenz96_main