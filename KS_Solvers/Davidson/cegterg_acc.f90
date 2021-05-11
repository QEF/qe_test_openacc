!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cegterg_acc( h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &
                    npw, npwx, nvec, nvecx, npol, evc, ethr, &
                    e, btype, notcnv, lrot, dav_iter, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
#if defined(__CUDA)
  use cublas
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
                            nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum, mp_gather, mp_bcast, mp_size,&
                            mp_type_create_column_section, mp_type_free
  USE device_fbuff_m,      ONLY : gbuf => pin_buf
  USE device_memcpy_m, ONLY : dev_memcpy, dev_memset, dev_memcpy
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! umber of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol,nvec)
    !  evc contains the  refined estimates of the eigenvectors  
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! total number of indivitual hpsi
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters
  INTEGER :: n_start, n_end, my_n
  INTEGER :: column_section_type
    ! defines a column section for communication
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
    ! receive counts and memory offsets
  COMPLEX(DP), POINTER  :: pinned_buffer(:,:)
    ! auxiliary variable for performing MPI operation and overcome CUDAFortran limitations
  REAL(DP), ALLOCATABLE :: ew_host(:)
  REAL(DP), ALLOCATABLE :: e_host(:)
    ! auxiliary variables for performing dot product
  INTEGER :: i,j,k, ipol
    !
  !
  REAL(DP), EXTERNAL :: KSddot
  !
  EXTERNAL  h_psi_gpu,    s_psi_gpu,    g_psi_gpu
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
!civn 
  write(*,*) 'using cegterg_acc'
!
  nhpsi = 0
  CALL start_clock( 'cegterg' ); !write(*,*) 'start cegterg' ; FLUSH(6)
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
!civn 
!$acc data deviceptr( evc(npwx*npol,nvec), e(nvec) )   
! 
  !
  ALLOCATE(  psi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx*npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ALLOCATE( sc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc ', ABS(ierr) )
  ALLOCATE( hc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc ', ABS(ierr) )
  ALLOCATE( vc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc ', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew ', ABS(ierr) )
!civn 
!$acc data create(  psi( npwx*npol, nvecx ), hpsi( npwx*npol, nvecx ), sc( nvecx, nvecx ), hc( nvecx, nvecx ),vc( nvecx, nvecx ),  ew( nvecx ) ) 
!$acc host_data use_device( psi, hpsi, sc, hc, vc, ew ) 
!$acc data if(uspp) create( spsi( npwx*npol, nvecx ) )  
!$acc host_data if(uspp) use_device( spsi )  
!
  ALLOCATE( ew_host( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew_host ', ABS(ierr) )
  ALLOCATE( e_host( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate e_host ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate conv ', ABS(ierr) )
  ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
  !
  ! This buffer is used to perform MPI calls with non-contiguous slices.
  ! In order to limit the number of allocated buffers, a rather large,
  ! but hopefully 'repetitive' size is selected (as of today buffers are
  ! selected according to the leading dimension(s) )
  !
  CALL gbuf%lock_buffer(pinned_buffer, (/nvecx, nvecx/), ierr)
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  CALL dev_memcpy(psi, evc, (/ 1 , npwx*npol /), 1, &
                                (/ 1 , nvec /), 1)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_gpu( npwx, npw, nvec, psi, hpsi ) ; nhpsi = nhpsi + nvec
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi, spsi )
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
  CALL start_clock( 'cegterg:init' )
  !
  CALL divide_all(inter_bgrp_comm,nbase,n_start,n_end,recv_counts,displs)
  CALL mp_type_create_column_section(sc(1,1), 0, nbase, nvecx, column_section_type)
  my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
  !
  if (n_start .le. n_end) &
  CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi, kdmx, hpsi(1,n_start), &
              kdmx, ZERO, hc(1,n_start), nvecx )
  !
  if (n_start .le. n_end) then
     !
     !pinned_buffer(1:nbase, n_start:n_end) = hc( 1:nbase, n_start:n_end )
     !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, hc( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
     CALL dev_memcpy( pinned_buffer, hc, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
     CALL mp_sum( pinned_buffer(1:nbase, n_start:n_end), intra_bgrp_comm )
     !hc( 1:nbase, n_start:n_end ) = pinned_buffer(1:nbase, n_start:n_end)
     !ierr = cudaMemcpy2D( hc(1, n_start) , nvecx, pinned_buffer( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
     CALL dev_memcpy( hc, pinned_buffer, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
     !
  end if
  CALL mp_gather( hc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
  !
  IF ( uspp ) THEN
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi, kdmx, spsi(1,n_start), kdmx, &
                 ZERO, sc(1,n_start), nvecx )
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi, kdmx, psi(1,n_start), kdmx, &
                 ZERO, sc(1,n_start), nvecx )
     !
  END IF
  !
  if ((n_start .le. n_end) .and. (mp_size(intra_bgrp_comm) > 1 )) then
     !pinned_buffer(1:nbase, n_start:n_end) = sc( 1:nbase, n_start:n_end )
     !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
     CALL dev_memcpy( pinned_buffer, sc, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
     CALL mp_sum( pinned_buffer( 1:nbase, n_start:n_end ), intra_bgrp_comm )
     !sc( 1:nbase, n_start:n_end ) = pinned_buffer(1:nbase, n_start:n_end)
     !ierr = cudaMemcpy2D( sc(1, n_start) , nvecx, pinned_buffer( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
     CALL dev_memcpy( sc, pinned_buffer, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
  end if
  CALL mp_gather( sc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
  !
  CALL mp_type_free( column_section_type )
  !
!$acc parallel loop
  DO n = 1, nbase
     !
     ! ... the diagonal of hc and sc must be strictly real
     !
     hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
     sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
     !
     DO m = n + 1, nbase
        !
        hc(n,m) = CONJG( hc(m,n) )
        sc(n,m) = CONJG( sc(m,n) )
        !
     END DO
     !
  END DO
  !
  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
     CALL dev_memset(vc, ZERO, (/1, nbase/), 1, (/1, nbase/), 1)
     !
!$acc parallel loop
     DO n = 1, nbase
        !
        e(n) = REAL( hc(n,n) )
        !
        vc(n,n) = ONE
        !
     END DO
     !
     CALL mp_bcast( e, root_bgrp_id, inter_bgrp_comm )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hc, sc, nvecx, ew, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
     CALL dev_memcpy (e, ew, (/ 1, nvec /), 1 )
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'cegterg:update' )
     !
     !  ======== FROM HERE =====
     !np = 0
     !
     !DO n = 1, nvec
     !   !
     !   IF ( .NOT. conv(n) ) THEN
     !      !
     !      ! ... this root not yet converged ... 
     !      !
     !      np = np + 1
     !      !
     !      ! ... reorder eigenvectors so that coefficients for unconverged
     !      ! ... roots come first. This allows to use quick matrix-matrix 
     !      ! ... multiplications to set a new basis vector (see below)
     !      !
     !      IF ( np /= n ) vc(:,np) = vc(:,n)
     !      !
     !      ! ... for use in g_psi
     !      !
     !      ew(nbase+np) = e(n)
     !      !
     !   END IF
     !   !
     !END DO
     ! ========= TO HERE, REPLACED BY =======

     CALL reorder_evals_cevecs(nbase, nvec, nvecx, conv, e, ew, vc)
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
     IF ( uspp ) THEN
        !
        if (n_start .le. n_end) &
        CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, spsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, psi(1,nb1), kdmx )
        !     
     ELSE
        !
        if (n_start .le. n_end) &
        CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, psi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, psi(1,nb1), kdmx )
        !
     END IF
! NB: must not call mp_sum over inter_bgrp_comm here because it is done later to the full correction
     !
!$acc parallel loop
     DO np=1,notcnv
!$acc loop
        DO ipol = 1, npol
!$acc loop
           DO k=1,npwx
             psi(k + (ipol-1)*npwx,nbase+np) = - ew(nbase+np)*psi(k + (ipol-1)*npwx,nbase+np)
           END DO
        END DO
     END DO
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, hpsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                 ONE, psi(1,nb1), kdmx )
     CALL mp_sum( psi(:,nb1:nbase+notcnv), inter_bgrp_comm )
     !
     ! clean up garbage if there is any
     IF (npw < npwx) CALL dev_memset(psi, ZERO, [npw+1,npwx], 1, [nb1, nbase+notcnv])
     IF (npol == 2)  CALL dev_memset(psi, ZERO, [npwx+npw+1,2*npwx], 1, [nb1, nbase+notcnv])
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi(1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     !!! == OPTIMIZE HERE ==
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew_host(n) = KSDdot( 2*npw, psi(1,nbn), 1, psi(1,nbn), 1 )
           !
        ELSE
           !
           ew_host(n) = KSDdot( 2*npw, psi(1,nbn), 1, psi(1,nbn), 1 ) + &
                        KSDdot( 2*npw, psi(npwx+1,nbn), 1, psi(npwx+1,nbn), 1 )
           !
        END IF
        !
     END DO
     !
     CALL mp_sum( ew_host( 1:notcnv ), intra_bgrp_comm )
!civn 
!     ew(1:notcnv) = ew_host(1:notcnv)
!$acc data copyin( ew_host( nvecx ) ) 
!$acc parallel loop
     DO i = 1, notcnv
       ew(i) = ew_host(i)
     END DO  
!$acc end data
     !
!$acc parallel loop
     DO i = 1,notcnv
!$acc loop
        DO ipol = 1,npol
!$acc loop
           DO k=1,npw
             psi(k + (ipol-1)*npwx,nbase+i) = psi(k+(ipol-1)*npwx,nbase+i)/SQRT( ew(i) )
           END DO
        END DO
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi_gpu( npwx, npw, notcnv, psi(:,nb1), hpsi(:,nb1) ) ; nhpsi = nhpsi + notcnv
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi(1,nb1), spsi(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
     CALL divide_all(inter_bgrp_comm,nbase+notcnv,n_start,n_end,recv_counts,displs)
     CALL mp_type_create_column_section(sc(1,1), nbase, notcnv, nvecx, column_section_type)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     !
     CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, hpsi(1,nb1), kdmx, psi(1,n_start), kdmx, &
                 ZERO, hc(nb1,n_start), nvecx )
     !
     if ((n_start .le. n_end) .and. (mp_size(intra_bgrp_comm) > 1 )) then
        !pinned_buffer(nb1:nbase+notcnv, n_start:n_end) = hc( nb1:nbase+notcnv, n_start:n_end )
        !ierr = cudaMemcpy2D( pinned_buffer(nb1, n_start) , nvecx, hc( nb1, n_start ), nvecx, notcnv, n_end-n_start+1 )
        CALL dev_memcpy( pinned_buffer, hc, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
        CALL mp_sum( pinned_buffer( nb1:nbase+notcnv, n_start:n_end ), intra_bgrp_comm )
        !hc( nb1:nbase+notcnv, n_start:n_end ) = pinned_buffer(nb1:nbase+notcnv, n_start:)
        !ierr = cudaMemcpy2D(  hc( nb1, n_start ), nvecx, pinned_buffer(nb1,n_start), nvecx, notcnv, n_end-n_start+1 )
        CALL dev_memcpy( hc, pinned_buffer, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
     end if
     CALL mp_gather( hc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
     !
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     IF ( uspp ) THEN
        !
        CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, spsi(1,nb1), kdmx, psi(1,n_start), kdmx, &
                    ZERO, sc(nb1,n_start), nvecx )
        !     
     ELSE
        !
        CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, psi(1,nb1), kdmx, psi(1,n_start), kdmx, &
                    ZERO, sc(nb1,n_start), nvecx )
        !
     END IF
     !
     if ( (n_start .le. n_end) .and. (mp_size(intra_bgrp_comm) > 1 ) ) then
        !pinned_buffer( nb1:nbase+notcnv, n_start:n_end ) = sc( nb1:nbase+notcnv, n_start:n_end )
        !ierr = cudaMemcpy2D( pinned_buffer(nb1, n_start) , nvecx, sc( nb1, n_start ), nvecx, notcnv, n_end-n_start+1 )
        CALL dev_memcpy( pinned_buffer, sc, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
        CALL mp_sum( pinned_buffer( nb1:nbase+notcnv, n_start:n_end ), intra_bgrp_comm )
        !sc( nb1:nbase+notcnv, n_start:n_end ) = pinned_buffer( nb1:nbase+notcnv, n_start:n_end )
        !ierr = cudaMemcpy2D(  sc( nb1, n_start ), nvecx, pinned_buffer(nb1,n_start), nvecx, notcnv, n_end-n_start+1 )
        CALL dev_memcpy( sc, pinned_buffer, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
     end if
     CALL mp_gather( sc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
     !
     CALL mp_type_free( column_section_type )
     !
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
!$acc parallel loop
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real
        !
        IF( n>=nb1 ) THEN
           hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
           sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
        ENDIF
        !
        DO m = MAX(n+1,nb1), nbase
           !
           hc(n,m) = CONJG( hc(m,n) )
           sc(n,m) = CONJG( sc(m,n) )
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hc, sc, nvecx, ew, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
     ! ... test for convergence (on the CPU)
     !
!civn 
!    ew_host(1:nvec) = ew(1:nvec)
!    e_host(1:nvec) = e(1:nvec)
!$acc data copy( ew_host(nvecx), e_host(nvecx) )
!$acc parallel loop
     DO i = 1, nvec
       ew_host(i) = ew(i)
       e_host(i) = e(i)
     END DO 
!$acc end data
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew_host(1:nvec) - e_host(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew_host(1:nvec) - e_host(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     CALL dev_memcpy (e, ew, (/ 1, nvec /) )
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. &
          nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'cegterg:last' )
        !
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, psi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, evc, kdmx )
        CALL mp_sum( evc, inter_bgrp_comm )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        CALL dev_memcpy(psi, evc, (/ 1, npwx*npol /), 1, &
                                      (/ 1, nvec /), 1)
        !
        IF ( uspp ) THEN
           !
           CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, spsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                       ZERO, psi(1,nvec+1), kdmx)
           CALL dev_memcpy(spsi, psi(:,nvec+1:), &
                                        (/1, npwx*npol/), 1, &
                                        (/1, nvec/), 1)
           CALL mp_sum( spsi(:,1:nvec), inter_bgrp_comm )
           !
        END IF
        !
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, hpsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, psi(1,nvec+1), kdmx )
       CALL dev_memcpy(hpsi, psi(:,nvec+1:), &
                                       (/1, npwx*npol/), 1, &
                                       (/1, nvec/), 1)
        CALL mp_sum( hpsi(:,1:nvec), inter_bgrp_comm )
        !
        ! ... refresh the reduced hamiltonian 
        !
        nbase = nvec
        !
        ! These variables are set to ZERO in the CUF Kernel below
        !hc(1:nbase,1:nbase) = ZERO
        !sc(1:nbase,1:nbase) = ZERO
        !vc(1:nbase,1:nbase) = ZERO
        !
!$acc parallel loop
        DO n = 1, nbase
           DO j = 1, nbase
              !
              IF ( j == n ) THEN
                 hc(j,n) = CMPLX( e(n), 0.0_DP ,kind=DP)
                 !
                 sc(j,n) = ONE
                 vc(j,n) = ONE
              ELSE
                 hc(j,n) = ZERO; sc(j,n) = ZERO; vc(j,n) = ZERO
              END IF
           END DO
           !
        END DO
        !
        CALL stop_clock( 'cegterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  CALL gbuf%release_buffer(pinned_buffer, ierr)
!civn 
!$acc end host_data
!$acc end data
!$acc end host_data
!$acc end data
  DEALLOCATE( recv_counts )
  DEALLOCATE( displs )
  DEALLOCATE( conv )
  DEALLOCATE( e_host, ew_host, ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
!civn 
!$acc end data
  !
  CALL stop_clock( 'cegterg' ); !write(*,*) 'stop cegterg' ; FLUSH(6)
  !call print_clock( 'cegterg' )
  !call print_clock( 'cegterg:init' )
  !call print_clock( 'cegterg:diag' )
  !call print_clock( 'cegterg:update' )
  !call print_clock( 'cegterg:overlap' )
  !call print_clock( 'cegterg:last' )
  !
  RETURN
  !
END SUBROUTINE cegterg_acc
