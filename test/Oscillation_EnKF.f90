!***********************
! Oscillation (EnKF)
! 2012. 8.29 K. Ito
!***********************

  program Oscillation_KF
  implicit none
  
  !*****
  !---Parameters
  integer, parameter :: nt_asm=100 ! Period of data assimilation
  integer, parameter :: nt_prd=100 ! Period of prediction
  integer, parameter :: obs_interval=40 ! Inteval of obsevation
  real(8), parameter :: mass=1.0d0
  real(8), parameter :: k=0.5d0
  real(8), parameter :: dump=0.3d0 ! Damping coefficient
  real(8), parameter :: dt=1.0d-2 ! Time step
  real(8), parameter :: pi=3.14159265358979d0
  integer, parameter :: nens=1000
  !---Physical variable
  real(8) :: x_t(0:nt_asm+nt_prd),  v_t(0:nt_asm+nt_prd)
  real(8) :: x_s(0:nt_asm+nt_prd),  v_s(0:nt_asm+nt_prd)
  real(8) :: x_da_m(0:nt_asm+nt_prd,1:nens), v_da_m(0:nt_asm+nt_prd,1:nens)
  real(8) :: x_da(0:nt_asm+nt_prd), v_da(0:nt_asm+nt_prd) !ensemble mean
  real(8) :: x_obs(0:nt_asm/obs_interval)
  !---Output control
  integer, parameter :: output_interval=20 ! Interval of output monitoring
  character(len=7) :: obs_chr(0:nt_asm)
  !---Matrix(element 1=>x, element 2=>v)
  real(8) :: M(2,2)  ! state transient matrix
  real(8) :: Pf(2,2) ! Forecast error covariance matrix
  real(8) :: Pa(2,2) ! Analysis error covariance matrix
  real(8) :: R(1,1)  ! Observation error covariance matrix
  real(8) :: Kg(2,1) ! Kalman gain
  real(8) :: H(1,2)  ! Observation operator
  !---Working variable
  integer :: it,iens
  real(8) :: x_innov
  real(8) :: x_prtb(nens),v_prtb(nens)
  real(8) :: Ptmp(2,2)
  real(8) :: noise1,noise2,gnoise ! Gaussian noise

  !*****Initialization of random number generator
  call random_seed()

  !*****Section 1: Initial value setting
  x_t(0)=5.0d0
  v_t(0)=0.0d0
  x_s(0)=4.0d0
  v_s(0)=1.0d0
  Pf(1,1)=1.0d0; Pf(1,2)=0.0d0; Pf(2,1)=0.0d0; Pf(2,2)=1.0d0
  Pa=Pf
  R(1,1)=0.1d0
  Kg(1:2,1)=0.0d0
  H(1,1)=1.0d0
  H(1,2)=0.0d0

  !*****Section 2: True Field and Observations
  do it=1,nt_asm+nt_prd
    !---forward time step
    x_t(it) = x_t(it-1) + dt * v_t(it-1)
    v_t(it) = -(k * dt / mass) * x_t(it-1)          &
    &       + (1.0d0 - dump * dt / mass) * v_t(it-1) 
    !---Observation
    if ((mod(it,obs_interval) == 0) .and. (it <= nt_asm)) then
      ! Generate Gaussian Noise (gnoise) from uniform random number 
      ! (based on Box-Muller Method)
      call random_number(noise1)
      call random_number(noise2)
      gnoise=sqrt(R(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
      ! Generate observation by adding Gaussian noise to true value
      x_obs(it/obs_interval)=x_t(it)+gnoise
    end if
  end do

  !*****Section 3: Simulation run without DA
  do it=1,nt_asm+nt_prd
    x_s(it) = x_s(it-1) + dt * v_s(it-1)
    v_s(it) = -(k * dt / mass) * x_s(it-1)          &
    &       + (1.0d0 - dump * dt / mass) * v_s(it-1) 
  end do

  !*****Section 4: Data Assimilation
  do iens=1,nens
    ! Generate Gaussian Noise (gnoise) from uniform random number 
    call random_number(noise1)
    call random_number(noise2)
    gnoise=sqrt(Pa(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
    x_da_m(0,iens)=x_s(0)+gnoise
    ! Generate Gaussian Noise (gnoise) from uniform random number 
    call random_number(noise1)
    call random_number(noise2)
    gnoise=sqrt(Pa(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
    v_da_m(0,iens)=v_s(0)+gnoise
  end do
  x_da(0)=sum(x_da_m(0,1:nens))/nens
  v_da(0)=sum(v_da_m(0,1:nens))/nens

  do it=1,nt_asm
    !---Section 4-1: time integration
    do iens=1,nens
      x_da_m(it,iens) = x_da_m(it-1,iens) + dt * v_da_m(it-1,iens)
      v_da_m(it,iens) = -(k * dt / mass) * x_da_m(it-1,iens)   &
      &       + (1.0d0 - dump * dt / mass) * v_da_m(it-1,iens) 
    end do
    !---Section 4-2: Kalman Filter
    ! Section 4-2-1: Obtain Pf
    if(mod(it,obs_interval) == 0) then
      x_da(it)=sum(x_da_m(it,1:nens))/nens
      v_da(it)=sum(v_da_m(it,1:nens))/nens
      Pf=0.0d0
      do iens=1,nens
        x_prtb(iens)=x_da_m(it,iens)-x_da(it)
        v_prtb(iens)=v_da_m(it,iens)-v_da(it)
        Pf(1,1)=Pf(1,1)+x_prtb(iens)**2/(nens-1)
        Pf(1,2)=Pf(1,2)+x_prtb(iens)*v_prtb(iens)/(nens-1)
        Pf(2,1)=Pf(2,1)+v_prtb(iens)*x_prtb(iens)/(nens-1)
        Pf(2,2)=Pf(2,2)+v_prtb(iens)**2/(nens-1)
      end do
      ! Section 4-2-2: Kalman gain: Weighting of model result and obs.
      ! (Note) Obsevation only in x --> component 1(x) only
      !        In this case, inverse matrix --> scalar inverse
      Kg(1,1)=Pf(1,1)/(R(1,1)+Pf(1,1))
      Kg(2,1)=Pf(2,1)/(R(1,1)+Pf(1,1))
      ! Section 4-2-3: calculate innovation and correction
      do iens=1,nens
        ! Generate Gaussian Noise (gnoise) from uniform random number 
        call random_number(noise1)
        call random_number(noise2)
        gnoise=sqrt(R(1,1))*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)
        x_innov=x_obs(it/obs_interval)+gnoise-x_da_m(it,iens)
        x_da_m(it,iens)=x_da_m(it,iens)+Kg(1,1)*x_innov
        v_da_m(it,iens)=v_da_m(it,iens)+Kg(2,1)*x_innov
      end do
      ! Section 4-2-4: analysis error covariance matrix
      Pa=Pf-matmul(matmul(Kg,H),Pf)
      Pf=Pa
      write(*,'(A,F8.2,A,5F10.3)') "time=",dt*it, ", Pa=", Pa
    end if
    x_da(it)=sum(x_da_m(it,1:nens))/nens
    v_da(it)=sum(v_da_m(it,1:nens))/nens
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
  write(*,*) "###    Ensemble Kalman Filter       ###"
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
  
  end program Oscillation_KF