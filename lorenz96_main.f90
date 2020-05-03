! Created on 2020.5.2
! @author: Toyo_Daichi

module lorenz96_main

  use common
  use lorenz96_prm

  public   :: ting_rk4, Lorenz96_core
  private  :: write_Lorenz96_output ,del_spaces

  contains
 
  subroutine ting_rk4(kt, x_in, x_out)
    implicit none

    integer, intent(in)  :: kt

    real(r_size), intent(in)  :: x_in(1:nx) 
    real(r_size), intent(out) :: x_out(1:nx)

    real(r_size), allocatable :: x(:), xtmp(:)
    real(r_size), allocatable :: q1(:), q2(:), q3(:), q4(:)

    ! --- Working variable
    integer :: ik
    
    allocate(x(1:nx), xtmp(1:nx))
    allocate(q1(1:nx), q2(1:nx), q3(1:nx), q4(1:nx))
    
    x(:) = x_in(:)
    
    ! --- time integration start
    do ik = 1, kt
      xtmp(:) = x(:)
      call Lorenz96_core(xtmp, q1)
      xtmp(:) = x(:) + 0.5d0 * q1(:)
      call Lorenz96_core(xtmp, q2)
      xtmp(:) = x(:) + 0.5d0 * q2(:)
      call Lorenz96_core(xtmp, q3)
      xtmp(:) = x(:) + q3(:)
      call Lorenz96_core(xtmp, q4)
      x(:) = x(:) + (q1(:) + 2.0d0*q2(:) + 2.0d0*q3(:) + q4(:))/6.0d0
      if ( opt_veach ) then
        call write_Lorenz96_output(ik, kt, x)
      end if
    end do
    x_out(:) = x(:)

    ! --- tidy up
    deallocate(xtmp, q1, q2, q3, q4)
    return
  end subroutine ting_rk4

  subroutine Lorenz96_core(x_in, x_out)
    implicit none
    real(r_size), intent(in)  :: x_in(1:nx)
    real(r_size), intent(out) :: x_out(1:nx)

    ! --- Working variable
    integer :: i
    
    x_out(1) = x_in(nx)*(x_in(2) - x_in(nx-1)) - x_in(1) + force
    x_out(2) = x_in(1)*(x_in(3) - x_in(nx)) - x_in(2) + force
    
    do i = 3, nx-1
      x_out(i) = x_in(i-1)*(x_in(i+1) - x_in(i-2)) - x_in(i) + force
    end do
    x_out(nx) = x_in(nx-1)*(x_in(1) - x_in(nx-2)) - x_in(nx) + force
    x_out = dt* x_out(:)
    
    return
  end subroutine Lorenz96_core

  subroutine write_Lorenz96_output(it, last_step, x_in)
    implicit none
    integer, intent(in)       :: it, last_step
    real(r_size), intent(in)  :: x_in(1:nx)
    character(1096)           :: linebuf
  

    if ( it == 1 ) then
      open(2, trim(output_file), status='replace')
      write(linebuf, )
      call del_spaces(linebuf)
      write(2, '(a)') trim(linebuf)
      write(6,*) '+++ err covariance matrix 1st. step'
      write(6,*) error_covariance_matrix(:,:)
        
    else if ( it /= 1 .and. it /= last_step) then
      write(6,*) '+++ err covariance matrix 2nd. step ~'
      write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
      call del_spaces(linebuf)
      write(2, '(a)') trim(linebuf)
      write(6,*) error_covariance_matrix(:,:)
        
    else if ( it == last_step ) then
      write(linebuf, '(8(f12.5, ","), f12.5)') error_covariance_matrix
      call del_spaces(linebuf)
      write(2, '(a)') trim(linebuf)
      write(6,*) '+++ err covariance matrix last step '
      write(6,*) error_covariance_matrix(:,:)
      close(2)
    end if
      
  end subroutine write_Lorenz96_output
  
  subroutine del_spaces(space)
    implicit none
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

    return
  end subroutine del_spaces

end module lorenz96_main
