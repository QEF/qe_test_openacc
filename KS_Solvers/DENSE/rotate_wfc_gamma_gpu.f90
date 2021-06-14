!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gamma_gpu( h_psi_gpu, s_psi_gpu, overlap, &
                             npwx, npw, nstart, nbnd, psi_d, evc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi_d(-G)=psi_d*(G), except G=0
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#else
#define cublasDGEMM dgemm
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, &
          nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp_bands_util, ONLY : gstart ! index of the first nonzero G 
  USE mp,            ONLY : mp_sum 
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi_d, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  LOGICAL :: overlap
    ! if .FALSE. : S|psi_d> not needed
  COMPLEX(DP) :: psi_d(npwx,nstart), evc_d(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, evc_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER                  :: npw2, npwx2
  COMPLEX(DP), ALLOCATABLE :: aux_d(:,:)
  REAL(DP),    ALLOCATABLE :: hr_d(:,:), sr_d(:,:), vr_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: aux_d, hr_d, sr_d, vr_d, en_d
#endif
  INTEGER :: n_start, n_end, my_n, i, j
  !
  EXTERNAL  h_psi_gpu,    s_psi_gpu
    ! h_psi(npwx,npw,nvec,psi_d,hpsi)
    !     calculates H|psi_d>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi_d> (if needed)
    !     Vectors psi_d,hpsi,spsi are dimensioned (npwx,npol,nvec)

  npw2  = 2 * npw
  npwx2 = 2 * npwx

  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )

  !
  ALLOCATE( aux_d(  npwx, nstart ) )
  ALLOCATE( hr_d( nstart, nstart ) )
  ALLOCATE( sr_d( nstart, nstart ) )
  ALLOCATE( vr_d( nstart, nstart ) )
  ALLOCATE( en_d( nstart ) )
  call start_clock('rotwfcg'); !write(*,*) 'start rotwfcg' ; FLUSH(6)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi_d(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) THEN
     !$cuf kernel do(1)
     DO i=1,nstart
        psi_d(1,i) = CMPLX( DBLE( psi_d(1,i) ), 0.D0,kind=DP)
     END DO
  END IF
  !
  call start_clock('rotwfcg:hpsi'); !write(*,*) 'start rotwfcg:hpsi' ; FLUSH(6)
  CALL h_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
  call stop_clock('rotwfcg:hpsi'); !write(*,*) 'stop rotwfcg:hpsi' ; FLUSH(6)
  !
  call start_clock('rotwfcg:hc'); !write(*,*) 'start rotwfcg:hc' ; FLUSH(6)
  hr_d=0.D0
  CALL divide(inter_bgrp_comm,nstart,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
  CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi_d, &
                    npwx2, aux_d(1,n_start), npwx2, 0.D0, hr_d(1,n_start), nstart )
  IF ( gstart == 2 ) call CGcudaDGER( nstart, my_n, -1.D0, psi_d, &
                                      npwx2, aux_d(1,n_start), npwx2, hr_d(1,n_start), nstart )
  CALL mp_sum( hr_d, inter_bgrp_comm )
  !
  CALL mp_sum( hr_d, intra_bgrp_comm )
  !
  sr_d=0.D0
  IF ( overlap ) THEN
     !
     CALL s_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
     !
     if (n_start .le. n_end) &
     CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi_d, &
                       npwx2, aux_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     IF ( gstart == 2 ) CALL CGcudaDGER( nstart, my_n, -1.D0, psi_d, &
                                         npwx2, aux_d(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi_d, &
                       npwx2, psi_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     IF ( gstart == 2 ) CALL CGcudaDGER( nstart, my_n, -1.D0, psi_d, &
                                         npwx2, psi_d(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !
  END IF
  CALL mp_sum( sr_d, inter_bgrp_comm )
  !
  CALL mp_sum( sr_d, intra_bgrp_comm )
  call stop_clock('rotwfcg:hc'); !write(*,*) 'stop rotwfcg:hc' ; FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotwfcg:diag'); !write(*,*) 'start rotwfcg:diag' ; FLUSH(6)
  CALL diaghg( nstart, nbnd, hr_d, sr_d, nstart, en_d, vr_d, me_bgrp, root_bgrp, intra_bgrp_comm )
  call stop_clock('rotwfcg:diag'); !write(*,*) 'stop rotwfcg:diag' ; FLUSH(6)
  call start_clock('rotwfcg:evc_d'); !write(*,*) 'start rotwfcg:evc_d' ; FLUSH(6)
  !
  !$cuf kernel do(1)
  DO i=1, nbnd
     e_d(i) = en_d(i)
  END DO
  !
  ! ... update the basis set
  !
  aux_d=(0.D0,0.D0)
  if (n_start .le. n_end) &
  CALL cublasDGEMM( 'N','N', npw2, nbnd, my_n, 1.D0, psi_d(1,n_start), &
                     npwx2, vr_d(n_start,1), nstart, 0.D0, aux_d, npwx2 )
  CALL mp_sum( aux_d, inter_bgrp_comm )
  !
  !$cuf kernel do(2)
  DO i=1, nbnd
     DO j=1, npwx
        evc_d(j,i) = aux_d(j,i)
     END DO
  END DO
  call stop_clock('rotwfcg:evc_d'); !write(*,*) 'stop rotwfcg:evc_d' ; FLUSH(6)
  !
  DEALLOCATE( en_d )
  DEALLOCATE( vr_d )
  DEALLOCATE( sr_d )
  DEALLOCATE( hr_d )
  DEALLOCATE( aux_d )
  call stop_clock('rotwfcg'); !write(*,*) 'stop rotwfcg' ; FLUSH(6)
  !call print_clock('rotwfcg')
  !call print_clock('rotwfcg:hpsi')
  !call print_clock('rotwfcg:hc')
  !call print_clock('rotwfcg:diag')
  !call print_clock('rotwfcg:evc_d')
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_gamma_gpu

! In principle this can go away .......
SUBROUTINE CGcudaDGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
!     .. Scalar Arguments ..
    DOUBLE PRECISION ::  ALPHA
    INTEGER          ::   INCX, INCY, LDA, M, N
!     .. Array Arguments ..
    DOUBLE PRECISION :: A( LDA, * ), X( * ), Y( * )
#if defined(__CUDA)
    attributes(device) :: A, X, Y
#endif
    CALL DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

END SUBROUTINE CGcudaDGER
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_wfc_gamma_gpu( h_psi_gpu, s_psi_gpu, overlap, &
                                  npwx, npw, nstart, nbnd, psi_d, evc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_wfc for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
  USE util_param,       ONLY : DP
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp_bands_util,    ONLY : gstart ! index of the first nonzero G
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi_d(npwx,nstart), evc_d(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, evc_d, e_d
#endif
  !
  ! ... local variables:
  !
  INTEGER             :: npw2, npwx2
  COMPLEX(DP), ALLOCATABLE :: aux_d(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: aux_d
#endif
  !
  INTEGER :: idesc(LAX_DESC_SIZE)
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  LOGICAL :: do_distr_diag_inside_bgrp
  INTEGER :: ortho_parent_comm
  INTEGER, ALLOCATABLE :: idesc_ip( :, :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  !
  INTEGER :: i
  !
  EXTERNAL  h_psi_gpu,    s_psi_gpu
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  call start_clock('protwfcg')
  !
  CALL laxlib_getval( do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp, &
       ortho_parent_comm = ortho_parent_comm )
  CALL desc_init( nstart, nx, la_proc, idesc, rank_ip, idesc_ip )
  !
  npw2  = 2 * npw
  npwx2 = 2 * npwx

  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )
  !
  ALLOCATE( aux_d(  npwx, nstart ) )
  ALLOCATE( hr( nx, nx ) )
  ALLOCATE( sr( nx, nx ) )
  ALLOCATE( vr( nx, nx ) )
  ALLOCATE( en( nstart ) )

  aux_d=(0.0_DP,0.0_DP)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  IF ( gstart == 2 ) THEN
     !$cuf kernel do(1)
     DO i=1,nstart
        psi_d(1,i) = CMPLX( DBLE( psi_d(1,i) ), 0.D0, kind=DP )
     END DO
  ENDIF
  !
  call start_clock('protwfcg:hpsi')
  CALL h_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
  call stop_clock('protwfcg:hpsi')
  !
  call start_clock('protwfcg:hc')
  CALL compute_distmat_gpu( hr, psi_d, aux_d )
  !
  IF ( overlap ) THEN
     !
     CALL s_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
     CALL compute_distmat_gpu( sr, psi_d, aux_d )
     !
  ELSE
     !
     CALL compute_distmat_gpu( sr, psi_d, psi_d )
     !
  END IF
  call stop_clock('protwfcg:hc')
  !
  ! ... Diagonalize
  !
  call start_clock('protwfcg:diag')
  IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg en and vr are the same across ortho_parent_comm
     ! only the first bgrp performs the diagonalization
     IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nstart, hr, sr, nx, en, vr, idesc )
     IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
       CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
       CALL mp_bcast( en, root_bgrp_id, inter_bgrp_comm )
     ENDIF
  ELSE
     CALL pdiaghg( nstart, hr, sr, nx, en, vr, idesc )
  END IF
  call stop_clock('protwfcg:diag')
  !
  e_d(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  call start_clock('protwfcg:evc')
  CALL refresh_evc_gpu( )
  !
  evc_d(:,:) = aux_d(:,1:nbnd)
  call stop_clock('protwfcg:evc')
  !
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( aux_d )
  !
  DEALLOCATE( idesc_ip )
  DEALLOCATE( rank_ip )
  call stop_clock('protwfcg')
  !
  RETURN
  !
CONTAINS
  !
  !
  SUBROUTINE compute_distmat_gpu( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
#if defined(__CUDA)
     use cublas
#endif
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     COMPLEX(DP), ALLOCATABLE :: work_d(:,:)
#if defined(__CUDA)
     attributes(DEVICE) :: v, w, work_d
#endif
     !
     ALLOCATE( work_d( nx, nx ) )
     !
     work_d = ZERO
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = 0.0d0
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = idesc_ip(LAX_DESC_NC, 1, ipc )
        ic = idesc_ip(LAX_DESC_IC, 1, ipc )
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = idesc_ip(LAX_DESC_NR, ipr, ipc )
           ir = idesc_ip(LAX_DESC_IR, ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0, v(1,ir), npwx2, w(1,ic), npwx2, 0.D0, work_d, nx )

           IF ( gstart == 2 ) &
              CALL CGcudaDGER( nr, nc, -1.D0, v(1,ir), npwx2, w(1,ic), npwx2, work_d, nx )

           work = work_d

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO

     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL laxlib_dsqmsym( nstart, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat_gpu
  !
  !
  SUBROUTINE refresh_evc_gpu( )
     !
#if defined(__CUDA)
     use cublas
#endif
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta
     !
     COMPLEX(DP), ALLOCATABLE :: work_d(:,:)
#if defined(__CUDA)
     attributes(DEVICE) :: work_d
#endif
     !
     ALLOCATE( work_d( nx, nx ) )
     !
     work_d = ZERO
     !
     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = idesc_ip(LAX_DESC_NC, 1, ipc )
        ic = idesc_ip(LAX_DESC_IC, 1, ipc )
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = 0.0d0

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = idesc_ip(LAX_DESC_NR, ipr, ipc )
              ir = idesc_ip(LAX_DESC_IR, ipr, ipc )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vr(:,1:nc), root, ortho_parent_comm )
                 !
                 work_d = vr
                 !
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, psi_d(1,ir), npwx2, work_d, nx, beta, aux_d(1,ic), npwx2 )
                 !
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 !
                 work_d = vtmp
                 !
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, psi_d(1,ir), npwx2, work_d, nx, beta, aux_d(1,ic), npwx2 )
                 !
              END IF
              !
              beta = 1.0d0

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )
     !
     DEALLOCATE( work_d )
     !
     RETURN
  END SUBROUTINE refresh_evc_gpu
  !
END SUBROUTINE protate_wfc_gamma_gpu
