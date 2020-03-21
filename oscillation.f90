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
  real(r_size), parameter  :: pi   = 3.1415926535r_size979d0
  
  ! --- Physical variable
  real(r_size) :: x_t(0:nt_asm+nt_prd), v_t(0:nt_asm+nt_prd)
  real(r_size) :: x_s(0:nt_asm+nt_prd), v_s(0:nt_asm+nt_prd)

  real(r_size) :: x_da(0:nt_asm+nt_prd), v_da(0:nt_asm+nt_prd)
  real(r_size) :: x_obs(0:nt_asm/obs_interval)
  
  ! --- Matrix(element 1:x, 2:v)
  real(r_size) :: M(2,2)   ! state transient matrix
  real(r_size) :: Pf(2,2)  ! Forecast error convariance matrix
  real(r_size) :: Pa(2,2)  ! Analysis error convariance matrix
  real(r_size) :: R(1,1)   ! Observation error convariance matrix
  real(r_size) :: Kg(2,1)  ! Kalman gain
  real(r_size) :: H(1,1)   ! Observation operator
  
  ! --- Output control
  integer, parameter  :: output_interval = 20
  character(len=7)    :: obs_chr(0:nt_asm)
  
  ! --- Working variable
  integer :: it
  integer :: irec
  integer :: ierr
  integer :: iflag
  real(r_size) :: x_innov
  real(r_size) :: Ptmp(2,2)
  real(r_size) :: noise1, noise2, Gnoise ! Gaussian noise
  
  ! --- Initialization of random number generator
  call random_seed()
  
  !======================================================================
  !
  !----------------------------------------------------------------------
  ! +++ open namelist
  !----------------------------------------------------------------------
  ! --- Input control
  namelist /initial_osc/ &
    x_t(0), v_t(0), x_s(0), v_s(0),     &
    Pf(1,1), Pf(1,2), Pf(2,1), Pf(2,2), &
    R(1,1), Kg(1:2, 1:1),               &
    H(1,1), H(1,2)
  
  open(1, file='mainOscillation.cnf')
    read(1, nml=initial, iostat=ierr)
  close(1)
  
  ! namelist check
  if (ierr < 0 ) then
    write(6,*) '   Msg : Main[ main /  initial_osc ] '
    write(6,*) '   Not found namelist.        '
    write(6,*) '   Use default values.        '
  else if (ierr>0) then
    write(6,*) '   Msg : Main[ main /  initial_osc ] '
    write(6,*) '   *** Warning : Not appropriate names in namelist !! Check !!'
    write(6,*) '   Stop : oscillation.f90              '
    stop
  end if

  