!
! Lorenz (1963) model
!

module lorenz63_prm
  
  use kinddef

  ! --- matrix size
  integer, public   :: nx, ny

  ! --- Lorenz parameter
  real(r_size), parameter :: sig  = 10.0d0
  real(r_size), parameter :: gamm = 28.0d0
  real(r_size), parameter :: b    = 2.666666666666666667d0 !(= 8/3)
  real(r_size), parameter :: dt   = 1.0d-2 ! Time step

  ! --- For calculation, Runge-Kutta method
  real(r_size), public :: x_cal(3), y_cal(3), z_cal(3)
  real(r_size), public :: x_k(4), y_k(4), z_k(4)

end module lorenz63_prm
