!
! Lorenz (1963) equation
!
!  dx/dt = -sig*(x+y)
!  dy/dt = -x*z+gamm*x-y
!  dz/dt = x*y-b*z

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
  r_y = gamm*x -y -x*z
  r_z = x*y - b*z
  
  return
end subroutine cal_Lorenz

subroutine Lorenz63_Runge_Kutta(  &
  x_in, y_in, z_in,               & ! IN: previous step score 
  x_out, y_out, z_out             & ! OUT: Runge-Kutta method score 
)
  
  use lorenz63_prm
  implicit none
  
  real(r_size), intent(in)    :: x_in, y_in, z_in
  real(r_size), intent(out)   :: x_out, y_out, z_out
  
  call cal_Lorenz(                           &
    x_in, y_in, z_in,                        & ! IN
    x_k(1), y_k(1), z_k(1)                   & ! OUT
  )
  
  x_cal(1) = x_in + 0.5*x_k(1)*dt
  y_cal(1) = y_in + 0.5*y_k(1)*dt
  z_cal(1) = z_in + 0.5*z_k(1)*dt
  
  call cal_Lorenz(                         &
  x_cal(1), y_cal(1), z_cal(1),            & ! IN
  x_k(2), y_k(2), z_k(2)                   & ! OUT
  )
  
  x_cal(2) = x_in + 0.5*x_k(2)*dt 
  y_cal(2) = y_in + 0.5*y_k(2)*dt 
  z_cal(2) = z_in + 0.5*z_k(2)*dt
  
  call cal_Lorenz(                         &
  x_cal(2), y_cal(2), z_cal(2),            & ! IN
  x_k(3), y_k(3), z_k(3)                   & ! OUT
  )
  
  x_cal(3) = x_in + x_k(3)*dt
  y_cal(3) = y_in + y_k(3)*dt
  z_cal(3) = z_in + z_k(3)*dt
  
  call cal_Lorenz(                         &
  x_cal(3), y_cal(3), z_cal(3),            & ! IN
  x_k(4), y_k(4), z_k(4)                   & ! OUT
  )
  
  x_out = x_in + dt * (x_k(1) + 2*x_k(2) + 2*x_k(3) + x_k(4)) / 6.0d0
  y_out = y_in + dt * (y_k(1) + 2*y_k(2) + 2*y_k(3) + y_k(4)) / 6.0d0
  z_out = z_in + dt * (z_k(1) + 2*z_k(2) + 2*z_k(3) + z_k(4)) / 6.0d0
  
  return
end subroutine Lorenz63_Runge_Kutta

subroutine cal_TL_Lorenz(   &
  dx, dy, dz, xb, yb, zb,   & ! IN
  fx, fy, fz                & ! OUT
  )
  
  use kinddef
  use lorenz63_prm
  
  implicit none
  real(r_size), intent(in)  :: dx, dy, dz
  real(r_size), intent(in)  :: xb, yb, zb
  real(r_size), intent(out) :: fx, fy, fz

  fx =      -sig*dx   +sig*dy 
  fy = (gamm-zb)*dx -1.0d0*dy  -xb*dz
  fz =        yb*dx    +xb*dy   -b*dz

  return
end subroutine cal_TL_Lorenz

subroutine TL_Lorez63_Runge_Kutta( &
  x_state, y_state, z_state,       & ! IN
  dx, dy, dz,                      & ! IN
  x_trend, y_trend, z_trend        & ! OUT
)

  use kinddef
  use lorenz63_prm

  implicit none
  real(r_size), intent(in)  :: x_state, y_state, z_state
  real(r_size), intent(in)  :: dx, dy, dz
  real(r_size), intent(out) :: x_trend, y_trend, z_trend

  call cal_TL_Lorenz(          &
    x_state, y_state, z_state, &
    dx, dy, dz,                &
    x_k(1), y_k(1), z_k(1)       &
    )
    
    x_cal(1) = x_state + 0.5*x_k(1)*dt
    y_cal(1) = y_state + 0.5*y_k(1)*dt
    z_cal(1) = z_state + 0.5*z_k(1)*dt
    
    call cal_TL_Lorenz(                      &
    x_cal(1), y_cal(1), z_cal(1),            & ! IN
    dx, dy, dz,                              & ! IN
    x_k(2), y_k(2), z_k(2)                   & ! OUT
    )
    
    x_cal(2) = x_state + 0.5*x_k(2)*dt 
    y_cal(2) = y_state + 0.5*y_k(2)*dt 
    z_cal(2) = z_state + 0.5*z_k(2)*dt
    
    call cal_TL_Lorenz(                      &
    x_cal(2), y_cal(2), z_cal(2),            & ! IN
    dx, dy, dz,                              & ! IN
    x_k(3), y_k(3), z_k(3)                   & ! OUT
    )
    
    x_cal(3) = x_state + x_k(3)*dt
    y_cal(3) = y_state + y_k(3)*dt
    z_cal(3) = z_state + z_k(3)*dt
    
    call cal_TL_Lorenz(                      &
    x_cal(3), y_cal(3), z_cal(3),            & ! IN
    dx, dy, dz,                              & ! IN
    x_k(4), y_k(4), z_k(4)                   & ! OUT
    )
    
  x_trend = x_state + dt * (x_k(1) + 2*x_k(2) + 2*x_k(3) + x_k(4)) / 6.0d0
  y_trend = y_state + dt * (y_k(1) + 2*y_k(2) + 2*y_k(3) + y_k(4)) / 6.0d0
  z_trend = z_state + dt * (z_k(1) + 2*z_k(2) + 2*z_k(3) + z_k(4)) / 6.0d0

  return
end subroutine TL_Lorez63_Runge_Kutta