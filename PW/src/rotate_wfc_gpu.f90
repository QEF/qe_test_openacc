!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gpu &
            ( npwx, npw, nstart, gstart, nbnd, psi_d, npol, overlap, evc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine (maybe it should be an interface) for
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi ( atomic or random wavefunctions ).
  ! ... Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  ! ... Calls h_psi, s_psi to calculate H|psi> ans S|psi>
  ! ... It only uses an auxiliary array of the same size as psi.
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : use_para_diag, gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, gstart, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
    ! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi_d(npwx*npol,nstart), evc_d(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP), INTENT(OUT) :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE)       :: psi_d, evc_d, e_d
#endif
  COMPLEX(DP), ALLOCATABLE :: psi_h(:,:), evc_h(:,:)
  REAL(DP), ALLOCATABLE    :: e_h(:)
  EXTERNAL h_psi_gpu, s_psi_gpu
  !
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
  !
  CALL start_clock_gpu( 'wfcrot' )

  IF( use_para_diag ) THEN
     !
     ! use data distributed subroutine
     !
     IF ( gamma_only ) THEN
        !
        CALL protate_wfc_gamma_gpu ( h_psi_gpu, s_psi_gpu, overlap, &
                                     npwx, npw, nstart, nbnd, psi_d, evc_d, e_d )
        !
     ELSE
        !
        CALL protate_wfc_k_gpu ( h_psi_gpu, s_psi_gpu, overlap, &
                                 npwx, npw, nstart, nbnd, npol, psi_d, evc_d, e_d )
        !
     END IF
     !
  ELSE
     !
     ! use serial subroutines
     !
     IF ( gamma_only ) THEN
        !
        CALL rotate_wfc_gamma_gpu ( h_psi_gpu, s_psi_gpu, overlap, &
                                    npwx, npw, nstart, nbnd, psi_d, evc_d, e_d )
        !
     ELSE
        !
        CALL rotate_wfc_k_gpu ( h_psi_gpu, s_psi_gpu, overlap, &
                                npwx, npw, nstart, nbnd, npol, psi_d, evc_d, e_d )
        !
     END IF
     !
  END IF
  !
  CALL stop_clock_gpu( 'wfcrot' )
  !
END SUBROUTINE rotate_wfc_gpu
