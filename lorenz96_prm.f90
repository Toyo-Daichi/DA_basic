!
! Lorenz (1996) model parameter
!

module lorenz96_prm
  use common
  private :: none 
  
  integer, public            :: nx  ! number of grid point
  real(r_size), save, public :: dt
  real(r_size), save, public :: force
  real(r_size), save, public :: oneday

end module lorenz96_prm
