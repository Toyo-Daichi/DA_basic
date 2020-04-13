!
! Lorenz (1963) Parameter
!

module lorenz63_prm
  
  use kinddef

  real(r_size), parameter :: sig  = 10.0d0
  real(r_size), parameter :: gamm = 28.0d0
  real(r_size), parameter :: b    = 2.666666666666666667d0 !(= 8/3)

  ! --- For calculation, Runge-Kutta method
  real(r_size) :: x_cal(3), y_cal(3), z_cal(3)
  real(r_size) :: x_k(4), y_k(4), z_k(4)

end module lorenz63_prm
