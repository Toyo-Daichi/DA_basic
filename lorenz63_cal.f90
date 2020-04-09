!
! Lorenz (1963) equation
!

subroutine cal_Lorenz( &
    x, y, z,           & ! IN : previous step score 
    r_x, r_y, r_z      & ! OUT: Lorenz63 score  
  )
  
  use kinddef
  use lorenz63_prm
  
  implicit none
  real(r_size), intent(in)  :: x, y, z
  real(r_size), intent(out) :: r_x, r_y, r_z

  r_x = -sig*x + sig*y
  r_y = gamm*x - y -x*z
  r_z = x*y - b*z

  return
end subroutine cal_Lorenz