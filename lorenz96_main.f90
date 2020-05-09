! Created on 2020.5.2 ~
! @author: Toyo_Daichi

program lorenz96_main
  use common
  use lorenz96_prm
  use lorenz96_cal
  
  implicit none

  ! --- setting parameter
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
  
  ! --- Output control
  logical, save         :: opt_veach = .false.
  character(256)        :: initial_file
  character(256)        :: output_file
  
  ! --- Working variable
  character(512)  :: linebuf
  character(24)   :: cfmt
  character(4)    :: cfmt_num
  real(r_size)    :: gnoise
  integer         :: spinup_period, normal_period
  integer         :: kt_oneday
  integer         :: ix, it, ierr

  !======================================================================
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------
  
  namelist /set_parm/ nx, dt, force, oneday
  namelist /set_exp/ tool, da_method, intg_method
  namelist /set_period/ spinup_period, normal_period
  namelist /set_mobs/ obs_xintv, obs_tintv
  namelist /output/ initial_file, output_file, opt_veach
  
  read(5, nml=set_parm, iostat=ierr)
  read(5, nml=set_exp, iostat=ierr)
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
    ny = int(nx/obs_xintv)
    nt = int((kt_oneday*normal_period)/obs_tintv)
    allocate(yt_obs(nt,ny))
  end if
  
  !======================================================================
  !
  ! --- Sec.2 lorenz96 calculation
  !----------------------------------------------------------------------
  ! +++ display namelist
  !----------------------------------------------------------------------

  write(6,*) 'Exp. setting            :: ', tool
  write(6,*) 'Data assimlation method :: ', da_method
  write(6,*) 'Integral method         :: ', intg_method
  
  if ( trim(tool) == 'spinup' ) then
    call com_randn(nx, x_true(1,:))
    x_true(1,:) = x_true(1,:)*5.0d0
    call ting_rk4(kt_oneday*spinup_period, x_true(1,:), x_out(1,:))
    
  else if ( trim(tool) == 'normal' ) then
    open(2, file=trim(initial_file), form='formatted', status='old')
    read(2,*) x_tmp
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
      open(2, file=trim(initial_file), form='formatted', status='replace')
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