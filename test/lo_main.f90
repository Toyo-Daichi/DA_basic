!
!  ****  Lorenz Chaos model  ****
!
   program lorenzrk

      use mod_para_l

      implicit none

      character(60):: apara
      real(8):: dt,   t
      integer:: nt,   it,   nft,  i

      real(8),allocatable:: X(:), Y(:), Z(:)

      real(8):: X_o(3), X_m(3), X_f(3)
      real(8):: k1(3),  k2(3),  k3(3),  k4(3)

      open ( 10, file='lo_para.dat', status='old' )

!  ****  read parameters  ****

      read ( 10,* )
      read ( 10,* ) apara
      read ( 10,* ) sig
      write ( 6,* ) apara
      write ( 6,* ) sig

      read ( 10,* ) apara
      read ( 10,* ) gamm
      write ( 6,* ) apara 
      write ( 6,* ) gamm

      read ( 10,* ) apara
      read ( 10,* ) b
      write ( 6,* ) apara
      write ( 6,* ) b

      read ( 10,* ) apara
      read ( 10,* ) dt
      write ( 6,* ) apara
      write ( 6,* ) dt
 
      read ( 10,* ) apara
      read ( 10,* ) nt
      write ( 6,* ) apara
      write ( 6,* ) nt

      read ( 10,* ) apara
      read ( 10,* ) nft
      write ( 6,* ) apara
      write ( 6,* ) nft

!  ****  Set dimension of X, Y, Z  ****

      allocate ( X(0:nt),Y(0:nt),Z(0:nt) )

!  ****  Inital condition  ****

      it = 0
      t  = 0.0d0

      read ( 10,* ) apara
      read ( 10,* ) X(0)
      write ( 6,* ) apara
      write ( 6,* ) X(0)

      read ( 10,* ) apara
      read ( 10,* ) Y(0)
      write ( 6,* ) apara
      write ( 6,* ) Y(0)

      read ( 10,* ) apara
      read ( 10,* ) Z(0)
      write ( 6,* ) apara
      write ( 6,* ) Z(0)

!  ****  open output file and write inital condition  ****

      open ( 20, file='lorenz1.dat' )

      write ( 20,* ) it, t, X(0), Y(0), Z(0)

!  ****  Time marching loop  ***********************************************************

      do it = 1,nt
        i = it - 1
        t = dfloat(it) * dt

!   Runge-Kutta method

        X_o(1) = X(i)
        X_o(2) = Y(i)
        X_o(3) = Z(i)

        call cal_RHS( X_o, k1 )

        X_m(:) = X_o(:) + 0.5d0 * dt * k1(:)

        call  cal_RHS( X_m, k2 )

        X_m(:) = X_o(:) + 0.5d0 * dt * k2(:)

        call  cal_RHS( X_m, k3 )
 
        X_f(:) = X_o(:) + dt * k3(:)

        call  cal_RHS( X_f, k4 )

        X_f(:) = X_o(:) + dt / 6.0d0 * ( k1(:) + 2.0d0*k2(:) + 2.0d0*k3(:) + k4(:) )

        X(i+1) = X_f(1)
        Y(i+1) = X_f(2)
        Z(i+1) = X_f(3)

!  ****  write results in each nft steps  ****

        if ( mod(it,nft) == 0 ) write (20,*) it, t, X(it), Y(it), Z(it)

      end do

!  ****  End of Time marching loop  ****************************************************

      close (10)
      close (20)

      stop

    end program lorenzrk
