! Created on 2020.5.2
! @author: Toyo_Daichi
!
! +++ reference
!  >> https://github.com/takemasa-miyoshi/letkf

module lorenz96_cal

  use common
  use lorenz96_prm

  public   :: ting_rk4, Lorenz96_core, tinteg_rk4_ptbmtx, del_spaces
  private  :: none

contains

  !======================================================================
  !
  ! --- Sec.1 lorenz96 calculation
  !----------------------------------------------------------------------
! +++ Lorenz96 equation
!----------------------------------------------------------------------
  
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

  !----------------------------------------------------------------------
  ! +++ Runge Kutta method
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
    end do
    x_out(:) = x(:)

    ! --- tidy up
    deallocate(xtmp, q1, q2, q3, q4)
    return
  end subroutine ting_rk4

  !======================================================================
  !
  ! --- Sec.2  Time integration of Perturbation Matrix
  !----------------------------------------------------------------------
  ! +++ M P M^T
  !----------------------------------------------------------------------
  subroutine tinteg_rk4_ptbmtx(   &
    alpha, kt, nx,                & ! IN:  loop num
    x_in,                         & ! IN:  input score
    Pa,                           & ! IN:  analysis matrix
    Pf                            & ! OUT: forecast matrix
  )

    implicit none
    real(r_size), intent(in)  :: alpha ! NL(x+alpha*dx) = NL(x) + alpha*dxf
    integer, intent(in)       :: kt, nx 
    real(r_size), intent(in)  :: x_in(1:nx)
    real(r_size), intent(in)  :: Pa(1:nx, 1:nx)
    real(r_size), intent(out) :: Pf(1:nx, 1:nx)

    ! --- Working variable
    real(r_size), allocatable :: work1(:), work2(:)
    integer :: i

    allocate(work1(1:nx), work2(1:nx))
    call ting_rk4(kt, x_in, work1)
    do i = 1, nx
      work2(:) = Pa(:,i) * alpha + x_in(:)
      call ting_rk4(kt, work2, work2)
      Pf(:,i) = (work2 - work1) / alpha
    end do

    ! --- tidy up
    deallocate(work1, work2)

    return 
  end subroutine tinteg_rk4_ptbmtx


  !======================================================================
  ! +++ Useful tools
  !----------------------------------------------------------------------
  ! >> Localization
  !    Schur Product
  !-----------------------------------------------------------------------
  SUBROUTINE enkf_schur(scale,dist,factor)
    IMPLICIT NONE
    REAL(r_size),INTENT(IN)  :: scale
    REAL(r_size),INTENT(IN)  :: dist
    REAL(r_size),INTENT(OUT) :: factor
    REAL(r_size) :: a,b
  
    a = scale * SQRT(10.0d0/3.0d0)
    b = dist / a

    IF( dist <= a ) THEN
      factor = 1.0d0 -0.25d0*b**5 + 0.5d0*b**4 + 5.0d0/8.0d0*b**3 &
        & - 5.0d0/3.0d0*b**2
    ELSE IF( dist <= 2*a ) THEN
      factor = 1.0d0/12.0d0*b**5 - 0.5d0*b**4 + 5.0d0/8.0d0*b**3 &
        & + 5.0d0/3.0d0*b**2 - 5.0d0*b + 4.0d0 - 2.0d0/3.0d0/b
    !ELSE IF( a == 0.0d0 ) THEN
    !  factor = 1.0d0
    ELSE
      factor = 0.0d0
    END IF
  
    RETURN
  END SUBROUTINE enkf_schur

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

end module lorenz96_cal
