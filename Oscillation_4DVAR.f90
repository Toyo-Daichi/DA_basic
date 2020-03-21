!***********************
! Oscillation (4DVAR)
! 2012. 8.29 K. Ito
!***********************

  program Oscillation_4DVAR
  implicit none
  
  !*****
  !---Parameters
  integer, parameter :: nt_asm=100 ! Period of data assimilation
  integer, parameter :: nt_prd=100 ! Period of prediction
  integer, parameter :: obs_interval=20 ! Inteval of obsevation
  real(8), parameter :: mass=1.0d0
  real(8), parameter :: k=0.5d0
  real(8), parameter :: dump=0.3d0 ! Damping coefficient
  real(8), parameter :: dt=1.0d-2 ! Time step
  real(8), parameter :: pi=3.14159265358979d0
  !---Physical variable
  real(8) :: x_t(0:nt_asm+nt_prd),  v_t(0:nt_asm+nt_prd)
  real(8) :: x_s(0:nt_asm+nt_prd),  v_s(0:nt_asm+nt_prd)
  real(8) :: x_da(0:nt_asm+nt_prd), v_da(0:nt_asm+nt_prd)
  real(8) :: x_tmp(0:nt_asm), v_tmp(0:nt_asm)
  real(8) :: x_b, v_b ! First guess of initial value
  real(8) :: x_obs(0:nt_asm/obs_interval)
  !---Cost function and adjoint variable
  real(8) :: J ! Costfunction
  real(8) :: Jold ! Costfunction in a previous iteration
  real(8) :: adx(0:nt_asm), adv(0:nt_asm)
  real(8) :: x_save(0:nt_asm), v_save(0:nt_asm)
  real(8) :: alpx = 0.02, alpv=0.02 ! Coefficient for minimization
  integer :: iter_max = 500 ! maximum number of iteration
  real(8) :: cond_iter = 1.0d-4 ! condition for iteration end ((Jold-J)/Jold)
  !---Output control
  integer, parameter :: output_interval=20
  character(len=7) :: obs_chr(0:nt_asm)
  !---Matrix(element 1:x, element 2:v)
  real(8) :: B(2,2)  ! Background error covariance matrix
  real(8) :: R(1,1)  ! Observation error covariance matrix
  !---Working var iable
  integer :: it
  real(8) :: x_innov
  integer :: iter ! number of iteration
  real(8) :: noise1,noise2,gnoise ! Gaussian noize

  !*****Initialization of random number generator
  call random_seed()

  !*****Section 1: Initial value setting
  x_t(0)=5.0d0
  v_t(0)=0.0d0
  x_s(0)=4.0d0
  v_s(0)=1.0d0
  B(1,1)=1.0d0 ; B(1,2)=0.0d0
  B(2,1)=0.0d0 ; B(2,2)=1.0d0
  R(1,1)=0.1d0
  
  !*****Section 2: True Field and Observations
  do it=1,nt_asm+nt_prd
    !---forward time step
    x_t(it) = x_t(it-1) + dt * v_t(it-1)
    v_t(it) = -(k * dt / mass) * x_t(it-1)             &
    &       + (1.0d0 - dump * dt / mass) * v_t(it-1)
    !---Observation
    if ((mod(it,obs_interval) == 0) .and. (it <= nt_asm)) then
      ! Generate Gaussian Noise (gnoise) from uniform random number 
      ! (based on Box-Muller Method)
      call random_number(noise1)
      call random_number(noise2)
      gnoise=sqrt(R(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
      ! Generate observation by adding Gaussian noise to true value
      x_obs(it/obs_interval)=x_t(it) + gnoise
    end if
  end do

  !*****Section 3: Simulation without DA
  do it=1,nt_asm+nt_prd
    x_s(it) = x_s(it-1) + dt * v_s(it-1)
    v_s(it) = -(k * dt / mass) * x_s(it-1)          &
    &       + (1.0d0 - dump * dt / mass) * v_s(it-1)   
  end do

  !*****Section 4: Data assimilation: 4DVAR
  x_tmp(0)=x_s(0)
  v_tmp(0)=v_s(0)
  x_b = x_tmp(0)   !set First guess value
  v_b = v_tmp(0)   !set First guess value
  do iter=1,iter_max
    !---Section 4-1: initialization for adjoint model run
    J=0.0d0
    adx(0:nt_asm)=0.0d0
    adv(0:nt_asm)=0.0d0
    !---Section 4-2: Forward model run
    do it=1,nt_asm
      x_tmp(it) = x_tmp(it-1) + dt * v_tmp(it-1)
      v_tmp(it) = -(k * dt / mass) * x_tmp(it-1)         &
      &       + (1.0d0 - dump * dt / mass) * v_tmp(it-1)  
      if(mod(it,obs_interval) == 0) then
        !---Calculate costfunction
        J = J + 0.5 * (x_obs(it / obs_interval) - x_tmp(it)) ** 2 / R(1,1)
      end if
    end do
    !---Section 4-3: adjoint model run
    do it=nt_asm,1,-1
      if(mod(it,obs_interval) == 0) then
        !---Calculate misfit and change adjoint variable
        adx(it) = adx(it) + (x_tmp(it)-x_obs(it/obs_interval))/R(1,1)
      end if
      adx(it-1) = adx(it) - k * dt / mass * adv(it) 
      adv(it-1) = dt * adx(it) + (1.0d0 - dump * dt / mass) * adv(it)
    end do
    !---Section 4-4: Consider background covariance
    J = J + 0.5d0 * (x_tmp(0)-x_b)**2 / B(1,1) &
    &     + 0.5d0 * (v_tmp(0)-v_b)**2 / B(2,2)
    adx(0)=adx(0)+(x_tmp(0)-x_b)/B(1,1)
    adv(0)=adv(0)+(v_tmp(0)-v_b)/B(2,2)
    !---Section 4-5: Check the end of iteration
    if((iter > 1) .and. (Jold < J)) then
      x_da(0:nt_asm)=x_tmp(0:nt_asm)
      v_da(0:nt_asm)=v_tmp(0:nt_asm)
      write(*,*)
      write(*,*) "Cost function increases from the previous iteration."
      write(*,*) "Replace DA results with those in previous iteration and exit DA."
      write(*,'(A,1X,i3,1X,A,3(F9.3,1X))') &
      &   "iteration=",iter, "; adjoint x,v=", adx(0),adv(0)
      write(*,'(10X,A,F9.3,1X,A,1X,3(F7.3,2X))') "J=",J, &
      &   "; x(0),v(0)=>",x_da(0),v_da(0)
      write(*,*)
      exit
    endif
    if((iter > 1) .and. (cond_iter > (Jold-J) / Jold)) then
      x_da(0:nt_asm)=x_tmp(0:nt_asm)
      v_da(0:nt_asm)=v_tmp(0:nt_asm)
      write(*,*)
      write(*,*) "Differences between J and Jold become small => exit DA."
      write(*,'(A,1X,i3,1X,A,3(F9.3,1X))') &
      &   "iteration=",iter, "; adjoint x,v=", adx(0),adv(0)
      write(*,'(10X,A,F9.3,1X,A,1X,3(F7.3,2X))') "J=",J, &
      &   "; x(0),v(0)=>",x_da(0),v_da(0)
      write(*,*)
      exit
    endif
    if(iter == iter_max) then
      x_da(0:nt_asm)=x_tmp(0:nt_asm)
      v_da(0:nt_asm)=v_tmp(0:nt_asm)
      write(*,*)
      write(*,*) "Maximum number of iteration reached"
      write(*,'(A,1X,i3,1X,A,3(F9.3,1X))') &
      &   "iteration=",iter, "; adjoint x,v=", adx(0),adv(0)
      write(*,'(10X,A,F9.3,1X,A,1X,3(F7.3,2X))') "J=",J, &
      &   "; x(0),v(0)=>",x_da(0),v_da(0)
      write(*,*)
      exit
    endif
    !---Section 4-6: save values
    Jold=J
    !---Section 4-7: Update of x and v
    x_tmp(0)=x_tmp(0)-alpx*adx(0)
    v_tmp(0)=v_tmp(0)-alpv*adv(0)
    !---Print the result of iteration
    write(*,'(A,1X,i3,1X,A,1X,2(F7.3,2X,A),F9.3)') &
    &    "iteration=",iter, "; x(0)=",x_tmp(0),"; v(0)=",v_tmp(0), "; J=",J
    write(*,'(6X,2(A,1X,F8.3,1X))') "adjoint x=", adx(0),"adjoint v=", adv(0)
  end do
  
  !*****Section 5: Prediction after Data Assimilation
  do it=nt_asm+1,nt_asm+nt_prd
    x_da(it) = x_da(it-1) + dt * v_da(it-1)
    v_da(it) = -(k * dt / mass) * x_da(it-1)          &
    &       + (1.0d0 - dump * dt / mass) * v_da(it-1)    
  end do

  !*****Output on screen
  !---Preparation output: conversion real --> character
  obs_chr(0:nt_asm)="No obs"
  do it=1,nt_asm
    if (mod(it,obs_interval) == 0) then
      write(obs_chr(it), '(F7.3)') x_obs(it/obs_interval)
    end if
  end do
  !---Output Main
  write(*,*) "#######################################"
  write(*,*) "###    Identical Twin Experiment    ###"
  write(*,*) "###    Adjoint method (page 137)    ###"
  write(*,*) "#######################################"
  write(*,*) 
  write(*,'(A,F7.2,A,F7.2)') "Assmilation Period: t=",        &
  &                          0.0,"-",dt*nt_asm
  write(*,'(A,F7.2,A,F7.2)') "Prediction Period:  t=",         &
  &                          dt*nt_asm,"-", dt*(nt_asm+nt_prd)
  write(*,*) 
  write(*,*) "*******Assimilation Period: x ********"
  write(*,*) " [Time]   [True]  [No Assim]  [Assim]  [Observation]"
  do it=0,nt_asm
    if (mod(it,output_interval) == 0) then
      write(*,'(F7.2,3F10.3,4X,A)') dt*it, x_t(it), x_s(it),x_da(it),obs_chr(it)
    endif
  end do
  write(*,*) "*******Prediction Period: x ********"
  do it=nt_asm+1,nt_asm+nt_prd
    if (mod(it,output_interval) == 0) then
      write(*,'(F7.2,3F10.3)') dt*it, x_t(it), x_s(it),x_da(it)
    endif
  end do
  write(*,*) 
  write(*,*) "*******Assimilation Period: v ********"
  write(*,*) " [Time]   [True]  [No Assim]  [Assim]"
  do it=0,nt_asm
    if (mod(it,output_interval) == 0) then
      write(*,'(F7.2,3F10.3)') dt*it, v_t(it), v_s(it), v_da(it)
    endif
  end do
  write(*,*) "*******Prediction Period: v ********"
  do it=nt_asm+1,nt_asm+nt_prd
    if (mod(it,output_interval) == 0) then
      write(*,'(F7.2,3F10.3)') dt*it, v_t(it), v_s(it), v_da(it)
    endif
  end do
  
  end program Oscillation_4DVAR
