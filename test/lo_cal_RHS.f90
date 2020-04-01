!
!    calculate RHS of Lorenz Equation
!
   subroutine cal_RHS( X, k )

      use mod_para_l

      implicit none

      real(8):: X(3), k(3)

!  Calculate RHS

      k(1) = -sig * X(1) + sig * X(2)
      k(2) = -X(1) * X(3) + gamm * X(1) - X(2) 
      k(3) = X(1) * X(2) - b * X(3)

      return

   end subroutine cal_RHS
