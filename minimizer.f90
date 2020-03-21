! =============================================================================
!
!  MINIMIZE: Minimizer for JNoVA
!            module to use lbfgs.f (Original is includeed in MSM4DVAR)
!
! =============================================================================
module minimizer
!
! -----------------------------------------------------------------------------
!
!  HISTORY:
!  2002. 7.31 Y.Honda / Modify for JNoVA
!
! -----------------------------------------------------------------------------
!
!  COMMENT: Require 'lbfgs.f'
!
! -----------------------------------------------------------------------------
!  use vardef
!
! -----------------------------------------------------------------------------
!   Declaration of Variables
! -----------------------------------------------------------------------------
!
  implicit none

  ! ====================================================================
  !   >>> Variables required lbfgs.f(LBFGS Minimization)
  ! ====================================================================
  integer(4),          parameter :: nsave = 10
  real(8), parameter :: epsln = 1.0d-4

  ! integer(4),          save, private              :: iter
  integer(4),          save, private              :: iprint(2)
  integer(4),          save, private              :: point
  logical,             save, private              :: diagco
  real(8), save, private, allocatable :: zdiag(:)
  real(8), save, private, allocatable :: zs(:), zy(:), zw(:)

contains
! =============================================================================
! 
!   INITIALIZE_MINIMIZER: Initialize Minimizer parameters
!
! =============================================================================
  subroutine Initialize_Minimizer(vsiz, iter, iflag)
!
! -----------------------------------------------------------------------------
!
!   HISTORY: 2002. 7.31 Y.Honda - First Code
!
! -----------------------------------------------------------------------------
    implicit none

    integer(4), intent(in) :: vsiz
    integer(4), intent(out) :: iter, iflag

    iter = 0

    iflag = 0

    iprint(1) = 1; iprint(2) = 0
    diagco = .false.

    allocate( zdiag(1:vsiz) )
    allocate( zs(1:vsiz * nsave) )
    allocate( zy(1:vsiz * nsave) )
    allocate( zw(1:vsiz + 2 * nsave) )

    zdiag = 0.0
    zs    = 0.0
    zy    = 0.0
    zw    = 0.0

    return
  end subroutine Initialize_Minimizer

! =============================================================================
! 
!   MINIMIZE: Minimize by VA15AD (LBFGS)
!
! =============================================================================
  subroutine Minimize(vsiz, xctl, costf, costg, iter, iflag)
!
! -----------------------------------------------------------------------------
!
!   HISTORY: 2002. 7.31 Y.Honda - First Code
!            2002. 9.12 Y.Honda 
!              - Automatic Justification of the precision
!                of variables, xctl, costf, costg
!                "R_DBLE" is always used in minimization process.
!
! -----------------------------------------------------------------------------
!
!   COMMENTS:
!
!    Return Value
!      iflag = 1 or 2 : continue to minimize
!      iflag =   0    : succeed to find the minimum point
!            <   0    : abnormal termination by various reason.
!                       see the output file, 'lbfgs.txt'.
!
!    Check the precision of variable using 'costf'.
!
! -----------------------------------------------------------------------------
!
    implicit none
    integer(4),          intent(in)    :: vsiz
    integer(4),          intent(inout) :: iter
    integer(4),          intent(out)   :: iflag
    real(8), intent(inout) :: xctl(vsiz)
    real(8), intent(inout) :: costf
    real(8), intent(inout) :: costg(vsiz)

    integer             :: i
    real(8) :: axctl(vsiz)
    real(8) :: acostf
    real(8) :: acostg(vsiz)

    call va15ad(vsiz, nsave, xctl, costf, costg, diagco, zdiag, iprint,     &
      &         epsln, zs, zy, point, zw, iflag, iter)

    return
  end subroutine Minimize

! =============================================================================
! 
!   TERMINATE_MINIMIZER: Terminate Minimizer
!
! =============================================================================
  subroutine Terminate_Minimizer
!
! -----------------------------------------------------------------------------
!
!   HISTORY: 2002. 8. 1 Y.Honda - First Code
!
! -----------------------------------------------------------------------------
    implicit none

    deallocate( zdiag )
    deallocate( zs )
    deallocate( zy )
    deallocate( zw )

    return
  end subroutine Terminate_Minimizer

end module minimizer
