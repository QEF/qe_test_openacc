!----------------------------------------------------------------------
MODULE exch_lda_acc
!----------------------------------------------------------------------
!! LDA exchange - first, simple experiment with openacc.
!
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE slater_d( rs, ex, vx )
!$acc routine (slater_d) seq
  !
  IMPLICIT NONE
  !!
  REAL(8), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(8), INTENT(OUT) :: ex
  !! Exchange energy (per unit volume)
  REAL(8), INTENT(OUT) :: vx
  !! Exchange potential
  !
  ! ... local variables
  !
  ex =  -0.687247939924714d0* 2.0d0/3.0d0/ rs
  vx = 4.d0 / 3.d0 * (-0.687247939924714d0) *2.0d0/3.0d0 / rs
  !
  RETURN
  !
END SUBROUTINE slater_d
!
END MODULE
