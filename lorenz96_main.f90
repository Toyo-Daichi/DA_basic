! Created on 2020.5.2 ~
! @author: Toyo_Daichi

program lorenz96_main
  use common
  use lorenz96_prm
  use lorenz96_cal
  
  implicit none

  ! --- setting parameter
  real(r_size), allocatable :: x_in(:)
  real(r_size), allocatable :: x_out(:)

  character(8)  :: tool, da_method
  character(12) :: intg_method
  
  ! --- Output control
  logical, save         :: opt_veach = .false.
  character(256)        :: initial_file
  character(256)        :: output_file
  
  ! --- Working variable
  character(1096) :: linebuf
  character(24)   :: cfmt
  character(4)    :: cfmt_num
  integer         :: spinup_period, normal_period
  integer         :: kt_oneday
  integer         :: it, ierr

  !======================================================================
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------

  namelist /set_parm/ nx, dt, force, oneday
  namelist /set_exp/ tool, da_method, intg_method
  namelist /set_period/ spinup_period, normal_period
  namelist /output/ initial_file, output_file, opt_veach
  
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

  allocate(x_in(nx), x_out(nx))

  !======================================================================
  !
  ! --- Sec.2 lorenz96 calculation
  ! +++ display namelist
  write(6,*) 'Exp. setting            :: ', tool
  write(6,*) 'Data assimlation method :: ', da_method
  write(6,*) 'Integral method         :: ', intg_method

  kt_oneday = int(oneday/dt) ! Unit change for 1day

  if ( trim(tool) == 'spinup' ) then  
    call com_randn(nx, x_in)
    x_in = x_in*5.0d0
    call ting_rk4(kt_oneday*spinup_period, x_in, x_out)
  else if ( trim(tool) == 'normal' ) then
    open(2, file=trim(initial_file), form='formatted', status='old')
      read(2,*) x_in
    close(2)
    call ting_rk4(kt_oneday*normal_period, x_in, x_out)
  end if
  
  !======================================================================
  !
  ! --- Sec.* Writing OUTPUT
  if ( opt_veach ) then
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Check Writing output system,  '
    if (nx .lt. 100 ) then 
      cfmt = '(xx(F12.7, ","), F12.7)'
      write(cfmt_num,"(I2)") nx-1
      cfmt(2:3) = cfmt_num
    else if (nx .ge. 100 ) then
      cfmt = '(xxx(F12.7, ","), F12.7)'
      write(cfmt_num, "(I3)") nx-1
      cfmt(2:4) = cfmt_num
    end if
    open(2, file=trim(initial_file), form='formatted', status='replace')
      write(linebuf, cfmt) x_out
      call del_spaces(linebuf)
      write(2,'(a)') linebuf
    close(2)
    write(6,*) ' && Successfuly output !!!        '  
 
  else if ( .not. opt_veach ) then
    write(6,*) '-------------------------------------------------------'
    write(6,*) '+++ Check calculation system,  '
    write(6,*) ' && Successfuly calculate !!!  '  

  end if

  stop
end program lorenz96_main