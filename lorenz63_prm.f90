!
! Lorenz (1963) model
!
!  dx/dt = -sig*(x+y)
!  dy/dt = -x*z+gamm*x-y
!  dz/dt = x*y-b*z

module lorenz63_prm
  
  use kinddef

  ! --- matrix size
  integer, parameter      :: Nx   = 3
  integer, parameter      :: Nobs = 2

  ! --- Lorenz parameter
  real(r_size), parameter :: sig  = 10.0d0
  real(r_size), parameter :: gamm = 28.0d0
  real(r_size), parameter :: b    = 2.666666666666666667d0 !(= 8/3)

  ! --- For calculation, Runge-Kutta method
  real(r_size) :: x_cal(3), y_cal(3), z_cal(3)
  real(r_size) :: x_k(4), y_k(4), z_k(4)

end module lorenz63_prm
