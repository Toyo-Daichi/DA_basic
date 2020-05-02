!
! Lorenz (1996) model
!

module lorenz96_prm
  
  use kinddef

  private

  ! --- Lorenz96 parameter
  integer, parameter, public :: nx = 40 ! number of grid point

  real(r_size), save, public :: dt = 0.005d0
  real(r_size), save, public :: force = 8.0d0
  real(r_size), save, public :: oneday = 0.2d0

end module lorenz96_prm
