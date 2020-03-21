! Created on 2020.3.21
! @author: Toyo_Daichi
! (ref: http://www.itonwp.sci.u-ryukyu.ac.jp/itokosk.html) 

program oscillation
  implicit none

  ! --- setting Parameters
  integer, parameter  :: nt_asm       = 400 ! Period of data assimilation
  integer, parameter  :: nt_prd       = 400 ! Period of prediction
  integer, parameter  :: obs_interval = 40  ! Interval of observation

  real(8), parameter  :: mass = 1.0d0
  real(8), parameter  :: k    = 0.5d0
  real(8), parameter  :: dump = 0.3d0  ! Damping coefficinet
  real(8), parameter  :: dt   = 1.0d-2 ! Time step
  real(8), parameter  :: pi   = 3.14159265358979d0
  
  ! --- Physical variable
  real(8) :: x_t(0:nt_asm+nt_prd), v_t(0:nt_asm+nt_prd)
  real(8) :: x_s(0:nt_asm+nt_prd), v_s(0:nt_asm+nt_prd)

  real(8) :: x_da(0:nt_asm+nt_prd), v_da(0:nt_asm+nt_prd)
  real(8) :: x_obs(0:nt_asm/obs_interval)

  


  
  