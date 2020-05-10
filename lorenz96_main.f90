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
  real(r_size), allocatable :: x_NoDA(:,:)
  real(r_size), allocatable :: x_DA(:,:)
  
  ! --- For spinup && OBS info
  real(r_size), allocatable :: x_out(:,:)
  real(r_size), allocatable :: x_tmp(:)
  real(r_size), allocatable :: yt_obs(:,:)
  ! ***
  real(r_size), parameter   :: size_noise_obs = 1.0d-2
  integer                   :: ny, nt
  integer                   :: obs_xintv, obs_tintv

  ! --- settign exp.
  character(8)  :: tool, da_method
  character(12) :: intg_method
  
  ! --- Data assimilation exp.
  ! *** Extend Kalamn Filter
  real(r_size), allocatable :: Pf(:,:)
  real(r_size), allocatable :: Pa(:,:)
  real(r_size), allocatable :: Kg(:,:)
  real(r_size), allocatable ::  H(:,:)
  real(r_size), allocatable ::  R(:,:)

  real(r_size), allocatable :: obs_matrix(:,:)
  real(r_size), allocatable :: obs_inv_matrix(:,:)
  real(r_size), allocatable :: work(:)

  ! --- Output control
  logical, save         :: opt_veach = .false.
  logical, save         :: da_veach  = .false.
  character(256)        :: initial_true_file
  character(256)        :: initial_sim_file
  character(256)        :: output_file
  
  ! --- Working variable
  character(512)  :: linebuf
  character(24)   :: cfmt
  character(4)    :: cfmt_num
  real(r_size)    :: gnoise
  integer         :: spinup_period, normal_period
  integer         :: kt_oneday
  integer         :: ix, it, il, ir, ierr, lda, ipiv, lwork

  real(r_size), parameter :: alpha = 0.05d0
  integer, parameter      :: one_loop=1

  !======================================================================
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------
  
  namelist /set_parm/ nx, dt, force, oneday
  namelist /set_exp/ tool, intg_method
  namelist /set_da_exp/ da_veach, da_method
  namelist /set_period/ spinup_period, normal_period
  namelist /set_mobs/ obs_xintv, obs_tintv
  namelist /output/ initial_true_file, initial_sim_file, output_file, opt_veach
  
  read(5, nml=set_parm, iostat=ierr)
  read(5, nml=set_exp, iostat=ierr)
  read(5, nml=set_da_exp, iostat=ierr)
  read(5, nml=set_period, iostat=ierr)
  read(5, nml=set_mobs, iostat=ierr)
  
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
  read(5, nml=output, iostat=ierr)
  
  kt_oneday = int(oneday/dt) ! Unit change for 1day
  if ( trim(tool) == 'spinup' ) then 
    allocate(x_true(1,nx), x_out(1, nx))
    allocate(x_tmp(nx))
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
      allocate(x_NoDA(0:kt_oneday*normal_period, nx))
      allocate(x_DA(0:kt_oneday*normal_period, nx))
      ! covariance matrix set.
      allocate(Pf(1:nx, 1:nx))
      allocate(Pa(1:nx, 1:nx))
      allocate(Kg(1:nx, 1:ny))
      allocate( H(1:ny, 1:nx))
      ! for kalmangain inverse calculate.
      allocate(R(1:ny, 1:ny))
      allocate(obs_matrix(1:ny, 1:ny))
      allocate(obs_inv_matrix(1:ny, 1:ny))
    end if
  end if

  !======================================================================
  !
  ! --- Sec.2 lorenz96 calculation
  !----------------------------------------------------------------------
  ! +++ display namelist
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
      call ting_rk4(it, x_true(it-1,:), x_true(it,:))
      
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
  if ( da_veach ) then
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Data assimilation exp. start '
    write(6,*) ' >> Data assimilation method  :: ', da_method
    
    !-------------------------------------------------------------------
    ! +++ making NoDA score
    !-------------------------------------------------------------------
    open(2, file=trim(initial_sim_file), form='formatted', status='old')
      read(2,*) x_tmp
    close(2)
    x_NoDA(0,:) = x_tmp
    do it = 1, kt_oneday*normal_period
      call ting_rk4(kt_oneday*normal_period, x_NoDA(it-1,:), x_NoDA(it,:))
    end do

    !-------------------------------------------------------------------
    ! +++ Data assimilation
    if ( trim(da_method) == 'KF' ) then
      ! --- initialize
      x_DA(0,:) = x_NoDA(0,:)
      Pf = 0.d0; Pa = 0.d0
      Kg = 0.d0;  H = 0.d0; R = 0.d0
      
      ! making identity matrix
      forall ( il=1:nx, ir=1:nx ) 
        Pf(il, ir) = 1.0d0 
        Pa(il, ir) = 1.0d0
      end forall
      forall ( il=1:nx, ir=1:ny ) Kg(il, ir) = 1.0d0
      forall ( il=1:ny, ir=1:nx )  H(il, ir) = 1.0d0
      forall ( il=1:ny, ir=1:ny )  R(il, ir) = size_noise_obs

      data_assim_loop :&
      do it = 1, kt_oneday*normal_period
        write(6,*) 'Data assim. time step: ', it
        call ting_rk4(kt_oneday*normal_period, x_DA(it-1,:), x_DA(it,:))

        if ( mod(it, obs_tintv)==0 ) then
          call tinteg_rk4_ptbmtx( alpha, one_loop, nx, x_DA(it,:), Pa, Pf)
          !-----------------------------------------------------------
          ! +++ making inverse matrix
          !-----------------------------------------------------------
          obs_matrix = matmul(matmul(H, Pf), transpose(H)) + R
          obs_inv_matrix = obs_matrix

          call dgetrf(ny, ny, obs_inv_matrix, lda, ipiv, ierr)
          call dgetri(ny, obs_inv_matrix, lda, ipiv, work, lwork, ierr)

          Kg = matmul(matmul(Pf, transpose(H)), obs_inv_matrix)
          x_DA(it,:) = x_DA(it,:) + matmul(Kg, (yt_obs(it/obs_tintv,:) - matmul(H, x_DA(it,:))))
          Pa = Pf - matmul(matmul(Kg, H), Pf)
        end if
        write(6,*) matmul(obs_matrix, obs_inv_matrix)
      end do &
      data_assim_loop

    else if ( trim(da_method) == 'EnKF' ) then
    end if
  end if

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
      cfmt = '(xxx(F12.7, ","), F12.7)'
      write(cfmt_num, "(I3)") nx-1
      cfmt(2:4) = cfmt_num
    end if
    
    ! select open file
    if ( trim(tool) == 'spinup' ) then
      open(2, file=trim(initial_true_file), form='formatted', status='replace')
      write(linebuf, cfmt) x_out(0,:)
      call del_spaces(linebuf)
      write(2,'(a)') linebuf
      
    else if ( trim(tool) == 'normal' ) then
      open(2, file=trim(output_file), form='formatted', status='replace')
      do it = 0, kt_oneday*normal_period, kt_oneday/4
        print *, it
        write(linebuf, cfmt) x_true(it,:)
        call del_spaces(linebuf)
        write(2,'(a)') linebuf
      end do
    endif
    close(2)
    write(6,*) ' && Successfuly output !!!        '  
    
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

end program lorenz96_main