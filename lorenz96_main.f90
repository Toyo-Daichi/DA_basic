! Created on 2020.5.2
! @author: Toyo_Daichi

module lorenz96_main

  use kinddef
  use lorenz96_prm

  private  :: none ! None
  public   :: ting_rk4, Lorenz96_core, com_randn, del_spaces

  contains

  !======================================================================
  ! +++ Methods of Lorenz96
  !
  ! --- Sec.1 Input control
  !----------------------------------------------------------------------
  ! +++ open namelist, allocation
  !----------------------------------------------------------------------

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
    do ik = 1, ik
      xtmp(:) = x(:)
      call Lorenz96_core(xtmp, q1)
      xtmp(:) = x(:) + 0.5d0 * q1(:)
      call Lorenz96_core(xtmp, q2)
      xtmp(:) = x(:) + 0.5d0 * q2(:)
      call Lorenz96_core(xtmp, q3)
      xtmp(:) = x(:) + q3(:)
      call Lorenz96_core(xtmp, q4)
      x(:) = x(:) + (q1(:) + 2.0d0*q2(:) + 2.0d0*q3(:) + q4(:))/6.0d0
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
  
  subroutine com_randn(ndim, var)
    implicit none
    
    integer, intent(in)       :: ndim
    real(r_size), intent(out) :: var(1:ndim)
    real(r_dble), parameter   :: pi=3.14159265358979
    
    real(r_size)  :: rnd(2)
    real(r_dble)  :: genrand_res53
    logical, save :: first = .true.
    
    ! --- Working variable
    integer :: idate(8)
    integer :: i, iseed

    if (first) then
      call date_and_time(values=idate)
      iseed = idate(8) + idate(7)*1000
      call init_gen_rand(iseed)
      first = .false.
    end if

    if (mod(ndim,2) == 0) then
      do i = 1, ndim/2
        rnd(1) = genrand_res53()
        rnd(2) = genrand_res53()
        var(i*2-1) = sqrt(-2.0d0*log(rnd(1))*sin(2.0d0*pi*rnd(2)))
        var(i*2) = sqrt(-2.0d0*log(rnd(1))*cos(2.0d0*pi*rnd(2)))
      end do
    else
      do i = 1, (ndim-1)/2
        rnd(1) = genrand_res53()
        rnd(2) = genrand_res53()
        var(i*2-1) = sqrt(-2.0d0*log(rnd(1))*sin(2.0d0*pi*rnd(2)))
        var(i*2) = sqrt(-2.0d0*log(rnd(1))*cos(2.0d0*pi*rnd(2)))
      end do
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(ndim) = sqrt(-2.0d0*log(rnd(1))*cos(2.0d0*pi*rnd(2)))
    end if
    return
  end subroutine com_randn
  
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
