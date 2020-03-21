!----------------------------------------
! Sample program: background error covariances
! Based on Lorenz (1996) model
!
! 2014/08/02 Kosuke Ito
!----------------------------------------
program Lorenz96_4DVAR

  use kinddef
  use minimizer

  implicit none

  ! numerical setting
  real(kind=r_size) :: init_time=0.0
  real(kind=r_size) :: dt=1.0d-4 ! time interval
  integer(4),parameter :: nx=40 ! number of grid points
  integer(4),parameter :: nv=nx
  integer(4),parameter :: nstep=20000 ! num of total step
  integer(4),parameter :: nstep_cycle=100 ! num of time steps in one assimilation window

  ! observations
  integer(4),parameter :: obs_tintv=10
  integer(4),parameter :: obs_xintv=20 
  integer(4),parameter :: ny=nx/obs_xintv
  real(kind=r_size) :: size_noise_obs=1.0d-2
  real(kind=r_size),allocatable :: yt_obs(:,:)

  ! true, background, and analysis state
  real(kind=r_size) :: x_t(nx),x_NoDA(nx),x_DA(nx)

  ! innovation
  real(kind=r_size) :: y_innov(ny)

  ! incremant
  real(kind=r_size) :: x_increment(nx)

  ! physical variable
  real(kind=r_size) :: x_DA_b(nx) ! 4dvar background
  real(kind=r_size) :: x_tmp(nx)

  ! forcing term
  real(kind=r_size) :: F=8.0d0

  ! covariance
  real(kind=r_size) :: vx_Bsqrt(nv,nx) ! squared root of background error covariances

  ! errors
  real(kind=r_size) :: total_error_NoDA,total_error_DA
  real(kind=r_size) :: error_NoDA,error_DA

  ! minimization
  ! nx: 3DVAR w/climatological B
  integer(4) :: iter,iflag
  integer(4) :: ivar_max=30

  ! cycles
  integer(4) :: ncycle

  ! adjoint variables
  real(kind=r_size) :: xt_Basic(nx,nstep_cycle)
  real(kind=r_size) :: x_ad(nx)

  ! control variable
  real(kind=r_size) :: costf
  real(kind=r_size) :: costf_b ! J_b
  real(kind=r_size) :: costf_o ! J_obs
  real(kind=r_size),allocatable :: v_dJodv(:) ! dJo/dv
  real(kind=r_size),allocatable :: v_dJbdv(:) ! dJb/dv
  real(kind=r_size),allocatable :: v_dJdv(:) ! dJ/dv
  real(kind=r_size),allocatable :: v_ctrl(:) ! v

  ! working variable
  real(kind=r_size) time,gnoise
  real(kind=r_size) :: x_d(nx)
  integer(4) :: it,ix,ix1,ix2,ivar
  integer(4) :: icycle ! 4dvar cycles

  ! --- set control variable
  allocate(v_dJodv(nv)     )
  allocate(v_dJbdv(nv)     )
  allocate(v_dJdv(nv)     )
  allocate(v_ctrl(nv)      )
  allocate(yt_obs(ny,nstep_cycle))

  ! --- initial check
  if(nstep_cycle < obs_tintv ) then
    write(*,*) "(nstep_cycle) cannot be smaller than (obs_tintv): stop"
    stop
  endif
  
  ! --- ncycle
  ncycle=nstep/nstep_cycle
  if(mod(nstep, nstep_cycle) .ne. 0) then
    write(*,*) "(nstep) must be (integer)x(nstep_cycle): stop"
    write(*,*) "nstep=",nstep," nstep_cycle=",nstep_cycle
    stop
  endif

  ! --- prepare to read true init
  x_t =(/  5.991663, -2.762711,  0.6026607, -0.3006261, 0.9643651, 2.501237, 4.998241,   &
  &        2.892934, -3.621429,  2.461397,   1.8771800, 6.2233970, 6.522135, -0.8861141, &
  &        4.818257,  7.467416, -3.218989,   2.5829040, 3.2673140, 6.761271, 3.266266,   &
  &       -2.644757,  5.393810,  5.37906,    1.7756600, 3.5628790, 6.026888, -2.026147,  &
  &        2.392012, -2.316495, -0.08800707, 8.2444690, 3.4294290, 4.909702, 5.242771,   &
  &       -2.830564,  6.541541,  4.651697,   1.7250000, 6.2584440  /)

  ! --- initial value setting for NoDA & DA
  time=init_time
  total_error_DA=0.0
  total_error_NoDA=0.0  
  x_NoDA = (/ &
  &  6.716308, 3.591086,  0.7971129, 3.654189,  3.710051,  5.701009,  6.782502, &
  & -2.996980, 2.160929,  2.218730,  4.203613,  5.745687, -0.4557182,-1.004548, &
  &  2.827042, 8.029896,  1.865615, -0.8915706,-1.471266,  3.359117,  6.136018, &
  &  6.040505, 6.758795, -2.895099, -1.628105,  3.872653,  1.848279,  2.076902, &
  &  8.664373, 0.838648, -3.903258,  3.6485,    0.2785346, 1.752599,  4.871539, &
  &  6.265904, 2.599438, -2.982158, -1.805941, -2.822848 /)
  x_DA=x_NoDA

  ! --- Setting of squared Background error covariance
  vx_Bsqrt=0

!   ! --- diagonal B^{1/2} matrix
!    do ix1=1,nx
!      vx_Bsqrt(ix1,  ix1)= 5.0
!    enddo

  ! --- prescribed B^{1/2} matrix
  vx_Bsqrt(1,1)=5.0; vx_Bsqrt(2,1)=0.2; vx_Bsqrt(3,1)=-0.5; vx_Bsqrt(nx-1,1)=-0.5; vx_Bsqrt(nx,1)= 0.2
  vx_Bsqrt(1,2)=0.2; vx_Bsqrt(2,2)=5.0; vx_Bsqrt(3,2)= 0.2; vx_Bsqrt(4   ,2)=-0.5; vx_Bsqrt(nx,2)=-0.5
  do ix1=3,nx-2
    vx_Bsqrt(ix1-2,ix1)=-0.5
    vx_Bsqrt(ix1-1,ix1)= 0.2
    vx_Bsqrt(ix1,  ix1)= 5.0
    vx_Bsqrt(ix1+1,ix1)= 0.2
    vx_Bsqrt(ix1+2,ix1)=-0.5
  enddo
  vx_Bsqrt(1,nx-1)=-0.5; vx_Bsqrt(nx-3,nx-1)=-0.5; vx_Bsqrt(nx-2,nx-1)= 0.2; vx_Bsqrt(nx-1,nx-1)= 5.0; vx_Bsqrt(nx,nx-1)=0.2
  vx_Bsqrt(1,nx  )= 0.2; vx_Bsqrt(2,nx     )=-0.5; vx_Bsqrt(nx-2,nx  )=-0.5; vx_Bsqrt(nx-1,nx  )= 0.2; vx_Bsqrt(nx,nx  )=5.0

  !*** cycles
  do icycle=1,ncycle

    ! --- true and NoDA calculation
    do it=1,nstep_cycle
      ! time
      time=time+dt
      ! true time integration
      call x_DvDtEuler_x(nx,dt,x_t,F)
      ! NoDA time integration
      call x_DvDtEuler_x(nx,dt,x_NoDA,F)
      ! generate obs
      if (mod(it,obs_tintv)==0) then    
        do ix=1,nx
          if(mod(ix,obs_xintv)==0) then
            call gaussian_noise(size_noise_obs,gnoise)
            yt_obs(ix/obs_xintv,it)=x_t(ix)+gnoise
          endif
        enddo
      endif  
    enddo

    ! *** 4DVAR
    ! --- save the first guess
    x_DA_b=x_DA
    x_tmp=x_DA_b

    ! --- iteration
    call initialize_minimizer(nv,iter,iflag)

    x_increment(1:nx)=0.0
    y_innov=0.0
    v_ctrl(1:nv)=0.0

    do ivar=1,ivar_max

      costf=0.0
      costf_o=0.0

      ! --- step 1-1: forward model run to calculate Jo term
      do it=1,nstep_cycle
        call x_DvDtEuler_x(nx,dt,x_tmp,F)
        xt_Basic(1:nx,it)=x_tmp(1:nx)
        if (mod(it,obs_tintv)==0) then
          do ix=1,nx
            if(mod(ix,obs_xintv)==0) then
              if (ivar==1) then
                y_innov(ix/obs_xintv) = yt_obs(ix/obs_xintv,it) - xt_Basic(ix,it)
              endif
              costf_o = costf_o &
              & + 0.5 * (xt_Basic(ix,it)-yt_obs(ix/obs_xintv,it))**2 / (size_noise_obs**2)
            endif
          enddo
        endif
      enddo

      ! --- step 1-2: backward model run to derive dJo/dx (=x_ad)
      x_ad(1:nx)=0.0
      do it=nstep_cycle,1,-1
        ! calculate cost_function
        if (mod(it,obs_tintv)==0) then
          do ix=1,nx
            if(mod(ix,obs_xintv)==0) then
              x_ad(ix) = x_ad(ix) &
              & + (xt_Basic(ix,it)-yt_obs(ix/obs_xintv,it)) / (size_noise_obs**2)
            endif
          enddo
        endif    
        ! backward time calculation
        call x_DvDtEulerAd_x(nx,dt,xt_Basic(1:nx,it),x_ad,F)
      enddo 

      ! --- step 2: conversion dJo/dx --> dJo/dv
      v_dJodv(1:nv)=0.0
      v_dJodv=matmul(vx_Bsqrt,x_ad)

      ! --- step 3: calculation Jb term and dJb/dv
      costf_b=0.0
      do ix=1,nx
        costf_b = costf_b + 0.5 * v_ctrl(ix) ** 2 
      enddo
      v_dJbdv=v_ctrl

      ! --- step 4: minimization
      v_dJdv = v_dJbdv + v_dJodv
      costf = costf_b + costf_o
      call minimize(nv,v_ctrl,costf,v_dJdv,iter,iflag)

      ! --- step 5: conversion delta v0 --> delta x0 and update
      x_increment(1:nx)=0.0
      x_increment=matmul(vx_Bsqrt,v_ctrl)
      x_tmp=x_DA_b+x_increment

      ! --- check for convergence of minimization
      if ((iflag==0).or.(iflag<0)) then
        exit
      endif

    end do

    ! end of minimizer
    call Terminate_minimizer

    ! --- update the initial condition and trajectory
    do ix=1,nx
      x_DA(ix)=x_DA_b(ix)+x_increment(ix)
    enddo
    do it=1,nstep_cycle
      call x_DvDtEuler_x(nx,dt,x_DA,F)
    enddo

    ! calculate errors
    error_NoDA=0
    do ix=1,nx
      error_NoDA=error_NoDA+(x_NoDA(ix)-x_t(ix))**2
    enddo
    error_NoDA=sqrt(error_NoDA/nx)
    total_error_NoDA=total_error_NoDA+error_NoDA

    error_DA=0
    do ix=1,nx
      error_DA=error_DA+(x_DA(ix)-x_t(ix))**2
    enddo
    error_DA=sqrt(error_DA/nx)
    total_error_DA=total_error_DA+error_DA

  write(*,*) "cycle: ", icycle, ", error in NoDA and DA: ",error_NoDA," ",error_DA

  enddo

  ! output total error
  write(*,*) "-------------"
  write(*,*) "total error ",total_error_NoDA," ",total_error_DA
       

end program Lorenz96_4DVAR

!*****************************************

subroutine x_DvDtEuler_x(nx,dt,x_val,F)

  use kinddef
  implicit none

  ! dimension
  integer,intent(in) :: nx  

  ! physical variable
  real(kind=r_size),intent(in)  :: dt,F
  real(kind=r_size),intent(inout)  :: x_val(nx)
  
  ! working variable
  real(kind=r_size) :: x_dval(nx)
  integer(4) :: it,ix

  do ix=3,nx-1
    x_dval(ix)=x_val(ix-1)*(x_val(ix+1)-x_val(ix-2))-x_val(ix)+F
  enddo
  x_dval(1 )=x_val(nx  )*(x_val(2)-x_val(nx-1))-x_val(1 )+F
  x_dval(2 )=x_val(1   )*(x_val(3)-x_val(nx  ))-x_val(2 )+F
  x_dval(nx)=x_val(nx-1)*(x_val(1)-x_val(nx-2))-x_val(nx)+F
  x_val=x_val+x_dval*dt
  
end subroutine x_DvDtEuler_x

!*****************************

subroutine x_DvDtEulerAd_x(nx,dt,x_Basic,x_ad,F)

  use kinddef
  implicit none

  integer,intent(in) :: nx  
  real(kind=r_size),intent(in) :: dt,F
  real(kind=r_size),intent(in) :: x_Basic(nx)
  real(kind=r_size),intent(inout) :: x_ad(nx)
  real(kind=r_size) :: x_dad(nx)

  integer :: ix

  x_dad(1:nx)=0.0

  do ix=3,nx-1
    x_dad(ix+1)=x_dad(ix+1)+x_basic(ix-1)*dt*x_ad(ix)
    x_dad(ix  )=x_dad(ix  )+(1-dt)*x_ad(ix)
    x_dad(ix-1)=x_dad(ix-1)+(x_basic(ix+1)-x_basic(ix-2))*dt*x_ad(ix)
    x_dad(ix-2)=x_dad(ix-2)-x_basic(ix-1)*dt*x_ad(ix)
  enddo
  x_dad(2   )=x_dad(2   )+x_basic(nx-1)*dt*x_ad(1)
  x_dad(1   )=x_dad(1   )+(1-dt)*x_ad(1)
  x_dad(nx  )=x_dad(nx  )+(x_basic(2)-x_basic(nx-1))*dt*x_ad(1)
  x_dad(nx-1)=x_dad(nx-1)-x_basic(nx)*dt*x_ad(1)

  x_dad(3   )=x_dad(3   )+x_basic(1)*dt*x_ad(2)
  x_dad(2   )=x_dad(2   )+(1-dt)*x_ad(2)
  x_dad(1   )=x_dad(1   )+(x_basic(3)-x_basic(nx-1))*dt*x_ad(2)
  x_dad(nx  )=x_dad(nx  )-x_basic(1)*dt*x_ad(2)

  x_dad(1   )=x_dad(1   )+x_basic(nx-1)*dt*x_ad(nx)
  x_dad(nx  )=x_dad(nx  )+(1-dt)*x_ad(nx)
  x_dad(nx-1)=x_dad(nx-1)+(x_basic(1   )-x_basic(nx-2))*dt*x_ad(nx)
  x_dad(nx-2)=x_dad(nx-2)-x_basic(nx-1)*dt*x_ad(nx)

  x_ad(1:nx)=x_dad(1:nx)

end subroutine x_DvDtEulerAd_x

!*****************************

subroutine gaussian_noise(size_gnoise,gnoise)

  use kinddef
  implicit none

  real(kind=r_size),intent(in) :: size_gnoise
  real(kind=r_size),intent(out) :: gnoise

  ! constant
  real(kind=r_dble),parameter :: pi=3.14159265358979

  ! working variable
  real(kind=r_dble) :: noise1,noise2

  call random_number(noise1)
  call random_number(noise2)
  gnoise=size_gnoise*sqrt(-2.0d0*log(1.0d0-noise1))*cos(2.0d0*pi*noise2)

end subroutine gaussian_noise
