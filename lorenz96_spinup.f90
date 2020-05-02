! Created on 2020.5.2
! @author: Toyo_Daichi

program lorenz96_spinup
  use kinddef
  use lorenz96_prm
  use lorenz96_main
  
  implicit none

  ! --- setting parameter
  real(r_size) :: x(nx)
  integer      :: kt_oneday

  ! --- Output control
  logical         :: opt_veach = .true.
  character(256)  :: output_file

  ! --- Working variable
  integer         :: it, ierr

  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------

  namelist /output/ output_file, opt_veach

  read(5, nml=output, iostat=ierr)

  kt_oneday = int(oneday/dt)

  call com_randn(nx, x)
  x = x*5.0d0
  call titeg_tk4(kt_oneday*360*100, x, x)

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