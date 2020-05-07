! Created on 2020.5.2 ~
! @author: Toyo_Daichi

program lorenz96_main
  use common
  use lorenz96_prm
  use lorenz96_cal
  
  implicit none

  ! --- setting parameter
  real(r_size), allocatable :: x_in(:,:)
  real(r_size), allocatable :: x_out(:,:)
  real(r_size), allocatable :: x_tmp(:)

  character(8)  :: tool, da_method
  character(12) :: intg_method
  
  ! --- Output control
  logical, save         :: opt_veach = .false.
  character(256)        :: initial_file
  character(256)        :: output_file
  
  ! --- Working variable
  character(512) :: linebuf
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
  read(5, nml=set_period, iostat=ierr)
  
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
    allocate(x_in(1, nx), x_out(1, nx))
    allocate(x_tmp(nx))
  else if ( trim(tool) == 'normal' ) then
    allocate(x_in(0:kt_oneday*normal_period, nx), x_out(0:kt_oneday*normal_period, nx))
    allocate(x_tmp(nx))
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
    call com_randn(nx, x_in(1,:))
    x_in(1,:) = x_in(1,:)*5.0d0
    call ting_rk4(kt_oneday*spinup_period, x_in(1,:), x_out(1,:))
    
  else if ( trim(tool) == 'normal' ) then
    open(2, file=trim(initial_file), form='formatted', status='old')
    read(2,*) x_tmp
    x_in(0, :) = x_tmp; x_out(0, :) = x_tmp
    do it = 1, kt_oneday*normal_period
      call ting_rk4(it, x_in(it-1,:), x_out(it,:))
      x_in(it, :) = x_out(it, :)
    end do
    close(2)
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
    
    ! select open file
    if ( trim(tool) == 'spinup' ) then
      open(2, file=trim(initial_file), form='formatted', status='replace')
        write(linebuf, cfmt) x_out(0,:)
        call del_spaces(linebuf)
        write(2,'(a)') linebuf

    else if ( trim(tool) == 'normal' ) then
      open(2, file=trim(output_file), form='formatted', status='replace')
      do it = 0, kt_oneday*normal_period
        write(linebuf, cfmt) x_out(it,:)
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
  deallocate(x_in, x_out, x_tmp)
  stop
end program lorenz96_main