! Created on 2020.5.2
! @author: Toyo_Daichi

program lorenz96_spinup
  use common
  use lorenz96_prm
  use lorenz96_main
  
  implicit none

  ! --- setting parameter
  real(r_size), allocatable :: x(:)

  character(8)  :: da_method
  character(12) :: intg_method
  
  ! --- Output control
  logical, save, public :: opt_veach = .true.
  character(256)        :: output_file
  
  ! --- Working variable
  integer         :: kt_oneday
  integer         :: it, ierr

  !======================================================================
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------

  namelist /set_parm/ nx, dt, force, oneday
  namelist /set_exp/ da_method, intg_method
  namelist /output/ output_file, opt_veach
  
  read(5, nml=set_parm, iostat=ierr)
  read(5, nml=set_exp, iostat=ierr)
  read(5, nml=output, iostat=ierr)

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

  allocate(x(nx))

  ! +++ display namelist
  write(6,*) 'Data assimlation method :: ', da_method
  write(6,*) 'Integral method         :: ', intg_method

  kt_oneday = int(oneday/dt)

  call com_randn(nx, x)
  x = x*5.0d0
  call ting_rk4(kt_oneday*360*100, x, x)

  ! --- Sec. Writing OUTPUT
  if ( opt_veach ) then
    
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Check Writing output system,  '
    write(6,*) ' && Successfuly output !!!        '  
 
  else if ( .not. opt_veach ) then
    write(6,*) x
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Check calculation system,  '
    write(6,*) ' && Successfuly calculate !!!  '  

  end if

  stop
end program lorenz96_spinup