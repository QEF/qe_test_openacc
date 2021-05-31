!
SUBROUTINE ppcg_k_idx( h_psi, s_psi, overlap, precondition, &
                 npwx, npw, nbnd, npol, psi, e, btype, &
                 ethr, maxter, notconv, avg_iter, sbsize, rr_step, scf_iter)
  !
  !----------------------------------------------------------------------------
  !
  ! E.V. Ignore btype, use ethr as threshold on subspace residual subspace
  ! SdG  restore btype use in the eigenvalue locking procedure
  ! TO DO:
  !        - use divide_all/mp_allgather instead of divie/mp_sum 
  !          (mp_type_create_column_section to get column_type) cfr ParO
  !        - allocate p, hp, sp as (sbsize) nbnd/ngrp 
  !
  USE util_param,         ONLY : DP, stdout
  USE mp,                 ONLY : mp_bcast, mp_root_sum, mp_sum
  USE mp_bands_util,      ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  !
  IMPLICIT NONE
  include 'laxlib.fh'
  COMPLEX (DP), PARAMETER :: C_ONE = (1.D0,0.D0), C_ZERO = (0.D0,0.D0)
  !
  ! ... I/O variables
  !
  LOGICAL,      INTENT(IN)    :: overlap ! whether the eigenvalue problem is a generalized one or not
  INTEGER,      INTENT(IN)    :: npwx, npw, nbnd, npol, maxter
  ! maximum number of PW for wavefunctions
  ! the number of plane waves
  ! number of bands
  ! number of independent spin components of the wfc.
  ! maximum number of iterations
  COMPLEX (DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  INTEGER,      INTENT(IN)    :: btype(nbnd) ! one if the corresponding state has to be
                                             ! ...converged to full accuracy, zero otherwise (src/pwcom.f90)
  REAL (DP),    INTENT(INOUT) :: e(nbnd)
  REAL (DP),    INTENT(IN)    :: precondition(npw), ethr
  ! the diagonal preconditioner
  ! the convergence threshold for eigenvalues
  INTEGER,      INTENT(OUT)   :: notconv
  REAL(DP),     INTENT(OUT)   :: avg_iter
  ! number of notconverged elements
  ! average number of iterations in PPCG
  INTEGER,      INTENT(IN)   ::  sbsize, rr_step      ! sub-block size (num. of vectors in sub-blocks)
                                                      ! ...to be used in PPCG block splitting. by default, sbsize=1
                                                      ! run the Rayleigh Ritz procedure every rr_step
  INTEGER, INTENT(IN)         :: scf_iter             ! this variable should be removed in future version, used for timing purpose
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE ::  hpsi(:,:), spsi(:,:), w(:,:), hw(:,:), sw(:,:), p(:,:), hp(:,:), sp(:,:)
  COMPLEX(DP), ALLOCATABLE ::  buffer1(:, :) 
  COMPLEX(DP), ALLOCATABLE ::  K(:,:), K_store(:,:), M(:,:), M_store(:,:), cwork(:)
  REAL(DP), ALLOCATABLE    ::  rwork(:)
  INTEGER,  ALLOCATABLE    ::  iwork(:)
  REAL (DP)                ::  trG, trG1, trdif, trtol, lock_tol
  COMPLEX(DP)              ::  coord_psi(sbsize,sbsize), coord_w(sbsize,sbsize), coord_p(sbsize,sbsize),  &
                               G(nbnd,nbnd), G1(nbnd, nbnd)
  REAL (DP)                ::  D(3*sbsize)
  INTEGER                  ::  nsb, sbsize_last, sbsize3, dimp,                        &
                               l, i, j, iter, total_iter, ierr = 0, info = 0,                       &
                               col_idx(sbsize), lcwork = -1, lrwork = -1, liwork = -1
  INTEGER                  ::  nact, nact_old, act_idx(nbnd), idx(nbnd)
                               ! number of active columns
                               ! number of active columns on previous iteration
                               ! indices of active columns
                               ! auxiliary indices
  INTEGER                  ::  n_start, n_end, my_n ! auxiliary indices for band group parallelization
  INTEGER                  ::  print_info     ! If > 0 then iteration information is printed
  INTEGER                  ::  kdim, kdimx, ipol, ibnd
  LOGICAL                  ::  clean
  !
  EXTERNAL h_psi, s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
  REAL(DP) :: rnrm 
  !
  COMPLEX (DP), ALLOCATABLE    ::  Gl(:,:)
  !
  INTEGER :: idesc(LAX_DESC_SIZE)
  ! descriptor of the current distributed Gram matrix
  LOGICAL :: la_proc
  ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  ! matrix distribution descriptors
  LOGICAL :: force_repmat    ! = .TRUE. to force replication of the Gram matrices for Cholesky
                             ! Needed if the sizes of these  matrices become too small after locking
  REAL                          :: res_array(maxter)

  INTEGER, PARAMETER   :: blocksz = 256 ! used to optimize some omp parallel do loops
  INTEGER   :: nblock
  !
  INTEGER :: np_ortho(2), ortho_parent_comm, ortho_cntx
  LOGICAL :: do_distr_diag_inside_bgrp
  !
  !
  res_array     = 0.0
  !
  CALL start_clock( 'ppcg_k' )
!civn 
  write(stdout,*) 'using ppcg_k_idx'
!
  !
  !  ... Initialization and validation
  CALL laxlib_getval( np_ortho = np_ortho, ortho_cntx = ortho_cntx, &
       ortho_parent_comm = ortho_parent_comm, &
       do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
  !
!civn 
  print_info =  3 
!
  sbsize3 = sbsize*3
  kdim  = npwx*(npol-1) + npw
  kdimx = npwx*npol
  clean = (npw < npwx) .AND. ( npol == 2 )
  if (npol> 2) CALL errore( 'ppcg ',' wrong npol value: npol > 2 ', npol )
  if (npol<=0) CALL errore( 'ppcg ',' non positive npol value: errcode = 1+abs(npol) ', 1+abs(npol) )
  !
  nact      =  nbnd
  nact_old  =  nbnd
  lock_tol  =  SQRT( ethr )  ! 1.D-4*SQRT( ethr )  ! 1.0D-7
  trdif     =  -1.D0
  ! trdif is used for stopping. It equals either
  ! the difference of traces of psi'*hpsi on two consecutive
  ! itrations or -1.D0 when this difference is undefined
  ! (e.g., initial iteration or the iteration following RR)
  trG    =  0.D0
  trG1   =  0.D0
  ! trG and trG1 are use to evaluate the difference of traces
  ! of psi'*hpsi on two consecutive iterations
  iter      =  1
  force_repmat = .FALSE.
  !
  CALL allocate_all
  !
  ! ... Compute block residual w = hpsi - psi*(psi'hpsi) (psi is orthonormal on input)
  !
  call start_clock('ppcg:hpsi')
  if (clean)  psi(npw+1:npwx,:) = C_ZERO
  CALL h_psi( npwx, npw, nbnd, psi, hpsi ) 
  if (clean) hpsi(npw+1:npwx,:) = C_ZERO
  if (overlap) then 
    CALL s_psi( npwx, npw, nbnd, psi, spsi) 
    if (clean) spsi(npw+1:npwx,:) = C_ZERO
  end if 
  avg_iter = 1.d0
  call stop_clock('ppcg:hpsi')
  !
  !     G = psi'hpsi
  call start_clock('ppcg:zgemm')
  G = C_ZERO
  CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
  if (n_start .le. n_end) &
  CALL ZGEMM('C','N', nbnd, my_n, kdim, C_ONE, psi, kdimx, hpsi(1,n_start), kdimx, C_ZERO, G(1,n_start), nbnd)
  CALL mp_sum( G, inter_bgrp_comm )
  CALL mp_sum( G, intra_bgrp_comm )
  call stop_clock('ppcg:zgemm')
  !
  !    w = hpsi - spsi*G
  call start_clock('ppcg:zgemm')
  IF ( my_bgrp_id /= root_bgrp_id ) call threaded_memset( w, 0.d0, 2*kdimx*nact )
  w(1:kdimx, 1:nact ) = hpsi(1:kdimx, 1:nact )
  CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
  if (overlap) then
     if (n_start .le. n_end) &
     CALL ZGEMM('N','N',kdim, nbnd, my_n, -C_ONE,spsi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, w, kdimx)
  else
     if (n_start .le. n_end) &
     CALL ZGEMM('N','N',kdim, nbnd, my_n, -C_ONE, psi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, w, kdimx)
  end if
  CALL mp_sum( w, inter_bgrp_comm )
  call stop_clock('ppcg:zgemm')
  !
  !
  ! ... Lock converged eigenpairs (set up act_idx and nact and store current nact in nact_old)
  call start_clock('ppcg:lock')
  nact_old = nact;
  CALL lock_epairs(kdim, nbnd, btype, w, kdimx, lock_tol, nact, act_idx)
  call stop_clock('ppcg:lock')
  !
  ! pack the active bands to the first columns (act_idx)
  CALL reshape_array(psi, kdimx, nact, nbnd, act_idx, buffer1, 1) 
  CALL reshape_array(w, kdimx, nact, nbnd, act_idx, buffer1, 1) 
  CALL reshape_array(hpsi, kdimx, nact, nbnd, act_idx, buffer1, 1) 
  if(overlap) CALL reshape_array(spsi, kdimx, nact, nbnd, act_idx, buffer1, 1) 
  !
  ! ... Set up iteration parameters after locking
  CALL setup_param
  !
  G1(1:nact,1:nact) = G(act_idx(1:nact), act_idx(1:nact))
  trG = get_trace( G1, nbnd, nact )
  !
  ! Print initial info ...
  IF (print_info >= 1)  THEN
     WRITE(stdout, '("Ethr: ",1pD9.2,", npw: ", I10, ", nbnd: ", I10, " , ",  &
              & "maxter: ",I5, ", sbsize:  ", I10,", nsb: ", I10 ,", nact: ", &
              & I10, ", trtol: ", 1pD9.2 )')  ethr, npw, nbnd, maxter, sbsize, nsb, nact, trtol
     IF (print_info == 3) THEN
        res_array(iter) = print_rnrm( w, kdim, nbnd, kdimx)
        WRITE(stdout,'("Res. norm:  ", 1pD9.2)') res_array(iter)
     END IF
     FLUSH( stdout )
  END IF
  !
  !---Begin the main loop
  !
  DO WHILE ( ((trdif > trtol) .OR. (trdif == -1.D0))  .AND. (iter <= maxter) .AND. (nact > 0) )
     !
     ! ... apply the diagonal preconditioner
     !
     nblock = (npw-1) / blocksz +1         ! used to optimize some omp parallel do loops
     !$omp parallel do collapse(3)
     DO j = 1, nact ; DO ipol=0,npol-1 ; DO i=1,nblock
        w(1+(i-1)*blocksz+npwx*ipol:MIN(i*blocksz,npw)+npwx*ipol, j) =    &
            w(1+(i-1)*blocksz+npwx*ipol:MIN(i*blocksz,npw)+npwx*ipol, j) / &
                precondition(1+(i-1)*blocksz:MIN(i*blocksz,npw))
     END DO ; END DO ; END DO
     !$omp end parallel do
     !
     call start_clock('ppcg:zgemm')
     G(1:nbnd,1:nact) = C_ZERO  
     CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
     if (overlap) then
        if (n_start .le. n_end) &
        CALL ZGEMM( 'C','N', my_n, nact, kdim, C_ONE, spsi(1,n_start), kdimx, w, kdimx, C_ZERO, G(n_start,1), nbnd )
     else
        if (n_start .le. n_end) &
        CALL ZGEMM( 'C','N', my_n, nact, kdim, C_ONE,psi(1,n_start), kdimx, w, kdimx, C_ZERO, G(n_start,1), nbnd )
     end if
     CALL mp_sum( G(1:nbnd,1:nact), inter_bgrp_comm )
     CALL mp_sum( G(1:nbnd,1:nact), intra_bgrp_comm )
     call stop_clock('ppcg:zgemm')
     !
     !     w = w - psi*G
     call start_clock('ppcg:zgemm')
     IF ( my_bgrp_id /= root_bgrp_id ) call threaded_memset( w, 0.d0, 2*kdimx*nact )  
     if (n_start .le. n_end) &
     CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, psi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, w, kdimx)
     CALL mp_sum( w(:,1:nact), inter_bgrp_comm )
     call stop_clock('ppcg:zgemm')
     !
     ! ... Compute h*w
     call start_clock('ppcg:hpsi')
     CALL h_psi( npwx, npw, nact, w, hw )      
     if (clean) hw (npw+1:npwx,1:nact) = C_ZERO
     if (overlap) then ! ... Compute s*w
        CALL s_psi( npwx, npw, nact, w, sw )   
        if (clean) sw(npw+1:npwx,1:nact) = C_ZERO
     end if
     avg_iter = avg_iter + nact/dble(nbnd)
     call stop_clock('ppcg:hpsi')
     !
     ! ... orthogonalize p against psi and w
!ev     IF ( MOD(iter, rr_step) /= 1 ) THEN    ! In this case, P is skipped after each RR
     IF ( iter  /=  1 ) THEN
        !  G = spsi'p
        call start_clock('ppcg:zgemm')
        G(1:nact,1:nact) = C_ZERO 
        CALL divide(inter_bgrp_comm,nact,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nact,n_start,n_end
        if (overlap) then
           if (n_start .le. n_end) &
           CALL ZGEMM('C','N', my_n, nact, kdim, C_ONE, spsi(1,n_start), kdimx, p, kdimx, C_ZERO, G(n_start,1), nbnd)
        else
           if (n_start .le. n_end) &
           CALL ZGEMM('C','N', my_n, nact, kdim, C_ONE, psi(1,n_start), kdimx, p, kdimx, C_ZERO, G(n_start,1), nbnd)
        end if
        CALL mp_sum( G(1:nact,1:nact), inter_bgrp_comm )
        CALL mp_sum( G(1:nact,1:nact), intra_bgrp_comm )
        call stop_clock('ppcg:zgemm')
        !
        ! p = p - psi*G, hp = hp - hpsi*G, sp = sp - spsi*G
        call start_clock('ppcg:zgemm')
        IF ( my_bgrp_id /= root_bgrp_id ) call threaded_memset( p, 0.d0, 2*kdimx*nact )  
        if (n_start .le. n_end) & ! could be done differently
        CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, psi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, p, kdimx)
        CALL mp_sum( p(:,1:nact), inter_bgrp_comm )
        call stop_clock('ppcg:zgemm')
        !
        call start_clock('ppcg:zgemm')
        IF ( my_bgrp_id /= root_bgrp_id ) call threaded_memset( hp, 0.d0, 2*kdimx*nact )
        if (n_start .le. n_end) &
        CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, hpsi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, hp, kdimx)
        CALL mp_sum( hp(:,1:nact), inter_bgrp_comm )
        call stop_clock('ppcg:zgemm')
        !
        if (overlap) then
           call start_clock('ppcg:zgemm')
           IF ( my_bgrp_id /= root_bgrp_id ) call threaded_memset( sp, 0.d0, 2*kdimx*nact )
           if (n_start .le. n_end) &
           CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, spsi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, sp, kdimx)
           CALL mp_sum( sp(:,1:nact), inter_bgrp_comm )
           call stop_clock('ppcg:zgemm')
        end if
     END IF
     !
     !  ... for each sub-block construct the small projected matrices K and M
     !      and store in K_store and M_store
     !
     K_store = C_ZERO
     M_store = C_ZERO
     !
     CALL divide(inter_bgrp_comm,nsb,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nsb,n_start,n_end
     DO j = n_start, n_end
        !
        ! Get size of the sub-block and define indices of the corresponding columns
        IF ( j < nsb )  THEN
           l = sbsize
        ELSE
           l = sbsize_last
        END IF
        col_idx(1:l) = act_idx(  (/ (i, i = (j-1)*sbsize + 1, (j-1)*sbsize + l) /) ) !civn try 2remove
        !
        ! ... form the local Gramm matrices (K,M)
        K = C_ZERO
        M = C_ZERO
        !
!civn sbsize --> l  ? 
        call start_clock('ppcg:zgemm')
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, hpsi(1,(j-1)*sbsize + 1), &
                                                                                        kdimx, C_ZERO, K, sbsize3)
        if (overlap) then
           if (clean) spsi(npw+1:npwx,(j-1)*sbsize+1:(j-1)*sbsize+l) = C_ZERO 
           CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, spsi(1,(j-1)*sbsize + 1), &
                                                                                       kdimx, C_ZERO, M, sbsize3)
        else
           CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, psi(1,(j-1)*sbsize + 1), &
                                                                                       kdimx, C_ZERO, M, sbsize3)
        end if
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, hw(1,(j-1)*sbsize + 1), &
                                                                           kdimx, C_ZERO, K(l+1, l+1), sbsize3)
        if (overlap) then
           CALL ZGEMM('C','N', l, l, kdim, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, sw(1,(j-1)*sbsize + 1), kdimx, &
                                                                                    C_ZERO, M(l+1, l+1 ), sbsize3)
        else
           CALL ZGEMM('C','N', l, l, kdim, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, w(1,(j-1)*sbsize + 1), kdimx, &
                                                                                     C_ZERO, M(l+1, l+1 ), sbsize3)
        end if
        if (clean) hw(npw+1:npwx,(j-1)*sbsize+1:(j-1)*sbsize+l) = C_ZERO
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, hw(1,(j-1)*sbsize + 1), kdimx, &
                                                                                        C_ZERO, K(1, l+1), sbsize3)
        if (overlap) then
           CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, sw(1,(j-1)*sbsize + 1), kdimx, &
                                                                                          C_ZERO, M(1, l+1), sbsize3)
        else
           CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, w(1,(j-1)*sbsize + 1), kdimx, &
                                                                                          C_ZERO, M(1, l+1), sbsize3)
        end if
        call stop_clock('ppcg:zgemm')
        !
        ! ---
        !
!ev        IF ( MOD(iter,rr_step) /= 1 ) THEN   ! In this case, P is skipped after each RR
        IF ( iter  /= 1 ) THEN
          call start_clock('ppcg:zgemm')
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, p(1,(j-1)*sbsize + 1), kdimx, hp(1,(j-1)*sbsize + 1), &
                                                                       kdimx, C_ZERO, K(2*l + 1, 2*l+1), sbsize3)
          if (overlap) then
             CALL ZGEMM('C','N', l, l, kdim, C_ONE, p(1,(j-1)*sbsize + 1), kdimx, sp(1,(j-1)*sbsize + 1), &
                                                                         kdimx, C_ZERO, M(2*l + 1, 2*l+1), sbsize3)
          else
             CALL ZGEMM('C','N', l, l, kdim, C_ONE, p(1,(j-1)*sbsize + 1), kdimx, p(1,(j-1)*sbsize + 1), &
                                                                         kdimx, C_ZERO, M(2*l + 1, 2*l+1), sbsize3)
          end if
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, hp(1,(j-1)*sbsize + 1), kdimx, &
                                                                                         C_ZERO, K(1, 2*l+1), sbsize3)
          if (overlap) then
             CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, sp(1,(j-1)*sbsize + 1), kdimx, &
                                                                                         C_ZERO, M(1, 2*l+1), sbsize3)
          else
             CALL ZGEMM('C','N', l, l, kdim, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, p(1,(j-1)*sbsize + 1), kdimx, &
                                                                                       C_ZERO, M(1, 2*l+1), sbsize3)
          end if
          call stop_clock('ppcg:zgemm')
          !
          ! ---
          !
          call start_clock('ppcg:zgemm')
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, hp(1,(j-1)*sbsize + 1), kdimx, &
                                                                                   C_ZERO, K(l+1, 2*l+1), sbsize3)
          if (overlap) then
             CALL ZGEMM('C','N', l, l, kdim, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, sp(1,(j-1)*sbsize + 1), kdimx, &
                                                                            C_ZERO, M(l+1, 2*l+1), sbsize3)
          else
             CALL ZGEMM('C','N', l, l, kdim, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, p(1,(j-1)*sbsize + 1), kdimx, &
                                                                              C_ZERO, M(l+1, 2*l+1), sbsize3)
          end if
          call stop_clock('ppcg:zgemm')
          !
        END IF
        !
        ! ... store the projected matrices
        K_store(:, (j-1)*sbsize3 + 1 : j*sbsize3 ) = K
        M_store(:, (j-1)*sbsize3 + 1 : j*sbsize3 ) = M
        !
     END DO
     CALL mp_sum(K_store,inter_bgrp_comm)
     CALL mp_sum(M_store,inter_bgrp_comm)
     CALL mp_sum(K_store,intra_bgrp_comm)
     CALL mp_sum(M_store,intra_bgrp_comm)
     !
     ! ... perform nsb 'separate RQ minimizations' and update approximate subspace

     idx(:) = 0 ! find the inactive columns to be kept by root_bgrp_id only
     idx(( n_start-1) * sbsize + 1 : min( n_end * sbsize, nact) ) = 1 
     if (my_bgrp_id == root_bgrp_id) then
        idx(nact + 1:nbnd) = 1 
     end if
     !
     DO j = n_start, n_end
       !
       ! Get size of the sub-block and define indices of the corresponding columns
       IF ( j < nsb )  THEN
          l = sbsize
       ELSE
          l = sbsize_last
       END IF

       col_idx(1:l) = act_idx( (/ (i, i = (j-1)*sbsize + 1, (j-1)*sbsize + l) /)  )
       !
       K = K_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
       M = M_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
       !
       lcwork =       6*sbsize +  9*sbsize**2
       lrwork =  1 + 15*sbsize + 18*sbsize**2
       liwork =  3 + 15*sbsize
!ev       IF ( MOD(iter, rr_step) /= 1) THEN    ! set the dimension of the separate projected eigenproblem
       IF ( iter /= 1 ) THEN    ! set the dimension of the separate projected eigenproblem
          dimp = 3*l
       ELSE
          dimp = 2*l
       END IF
       !
       CALL ZHEGVD(1, 'V','U', dimp, K, sbsize3, M, sbsize3, D, cwork, lcwork, rwork, lrwork, iwork, liwork, info)
       IF (info /= 0) THEN
       ! reset the matrix and try again with psi and w only
          K = K_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
          M = M_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
          dimp = 2*l
          CALL ZHEGVD(1, 'V','U', dimp, K, sbsize3, M, sbsize3, D, cwork, lcwork, rwork, lrwork, iwork, liwork, info)
          IF (info /= 0) THEN
             CALL errore( 'ppcg ',' zhegvd failed ', info )
             STOP
          END IF
       END IF
       !
       coord_psi(1 : l, 1 : l) = K(1 : l, 1 : l)
       coord_w(1 : l, 1 : l) = K(l+1 : 2*l, 1 : l)
       !
       ! ... update the sub-block of P and AP
!ev       IF ( MOD(iter, rr_step) /= 1 ) THEN
!sdg      IF ( iter /= 1 ) THEN
       IF ( dimp == 3*l ) THEN
          !
          coord_p(1 : l, 1 : l) = K(2*l+1 : 3*l, 1 : l)
          !
          call start_clock('ppcg:zgemm')
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, p(1,(j-1)*sbsize + 1), kdimx, coord_p, sbsize, C_ZERO, buffer1, kdimx)
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, coord_w, sbsize, C_ONE, buffer1, kdimx)
          p(:,(j-1)*sbsize + 1: (j-1)*sbsize + l)  = buffer1(:,1:l)
          call stop_clock('ppcg:zgemm')
          !
          call start_clock('ppcg:zgemm')
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, hp(1,(j-1)*sbsize + 1), kdimx, coord_p, sbsize, C_ZERO, buffer1, kdimx)
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, hw(1,(j-1)*sbsize + 1), kdimx, coord_w, sbsize, C_ONE, buffer1, kdimx)
          hp(:,(j-1)*sbsize + 1: (j-1)*sbsize + l)  = buffer1(:,1:l)
          call stop_clock('ppcg:zgemm')
          !
          if (overlap) then
             call start_clock('ppcg:zgemm')
             CALL ZGEMM('N','N', kdim, l, l, C_ONE, sp(1,(j-1)*sbsize + 1), kdimx, coord_p, sbsize, C_ZERO, buffer1, kdimx)
             CALL ZGEMM('N','N', kdim, l, l, C_ONE, sw(1,(j-1)*sbsize + 1), kdimx, coord_w, sbsize, C_ONE, buffer1, kdimx)
             sp(:,(j-1)*sbsize + 1: (j-1)*sbsize + l)  = buffer1(:,1:l)
             call stop_clock('ppcg:zgemm')
          end if
          !
       ELSE
          !
          call start_clock('ppcg:zgemm')
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, w(1,(j-1)*sbsize + 1), kdimx, coord_w, sbsize, C_ZERO, buffer1, kdimx)
          p(:,(j-1)*sbsize + 1: (j-1)*sbsize + l) = buffer1(:, 1:l)
          call stop_clock('ppcg:zgemm')
          !
          call start_clock('ppcg:zgemm')
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, hw(1,(j-1)*sbsize + 1), kdimx, coord_w, sbsize, C_ZERO, buffer1, kdimx)
          hp(:,(j-1)*sbsize + 1: (j-1)*sbsize + l) = buffer1(:, 1:l)
          call stop_clock('ppcg:zgemm')
          !
          if (overlap) then
             call start_clock('ppcg:zgemm')
             CALL ZGEMM('N','N', kdim, l, l, C_ONE, sw(1,(j-1)*sbsize + 1), kdimx, coord_w, sbsize, c_ZERO, buffer1, kdimx)
             sp(:,(j-1)*sbsize + 1: (j-1)*sbsize + l) = buffer1(:, 1:l)
             call stop_clock('ppcg:zgemm')
          end if
          !
       END IF
       !
       ! Update the sub-blocks of psi and hpsi (and spsi)
       call start_clock('ppcg:zgemm')
       CALL ZGEMM('N','N', kdim, l, l, C_ONE, psi(1,(j-1)*sbsize + 1), kdimx, coord_psi, sbsize, C_ZERO, buffer1, kdimx)
       psi(:, (j-1)*sbsize + 1: (j-1)*sbsize + l)  = buffer1(:,1:l) + p(:,(j-1)*sbsize + 1: (j-1)*sbsize + l)
       call stop_clock('ppcg:zgemm')
       !
       call start_clock('ppcg:zgemm')
       CALL ZGEMM('N','N', kdim, l, l, C_ONE, hpsi(1,(j-1)*sbsize + 1), kdimx, coord_psi, sbsize, C_ZERO, buffer1, kdimx)
       hpsi(:, (j-1)*sbsize + 1: (j-1)*sbsize + l)  = buffer1(:,1:l) + hp(:,(j-1)*sbsize + 1: (j-1)*sbsize + l)
       call stop_clock('ppcg:zgemm')
       !
       if (overlap) then
          call start_clock('ppcg:zgemm')
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, spsi(1,(j-1)*sbsize + 1), kdimx, coord_psi, sbsize, C_ZERO, buffer1, kdimx)
          spsi(:, (j-1)*sbsize + 1: (j-1)*sbsize + l) = buffer1(:,1:l) + sp(:,(j-1)*sbsize + 1: (j-1)*sbsize + l)
          call stop_clock('ppcg:zgemm')
       end if
       !
       idx(col_idx(1:l)) = 1 ! keep track of which columns this bgrp has acted on  
       !
     END DO  ! end 'separate RQ minimizations'
     ! set to zero the columns not assigned to this bgrp, inactive colums are assigned to root_bgrp
     do j=1,nbnd
        if (idx(j)==0) then
          psi (:,j) = C_ZERO 
          hpsi (:,j) = C_ZERO 
          p(:,j) = C_ZERO 
          hp(:,j) = C_ZERO
          if(overlap) then 
            spsi (:,j) = C_ZERO 
            sp(:,j) = C_ZERO
          end if 
        end if
     end do
     CALL mp_sum(psi ,inter_bgrp_comm)
     CALL mp_sum(hpsi,inter_bgrp_comm)
     CALL mp_sum(p ,inter_bgrp_comm)
     CALL mp_sum(hp,inter_bgrp_comm)
     if(overlap) then 
       CALL mp_sum(spsi,inter_bgrp_comm)
       CALL mp_sum(sp,inter_bgrp_comm)
     end if
    !
    !
    ! ... Perform the RR procedure every rr_step
!    IF ( (MOD(iter, rr_step) == 0) .AND. (iter /= maxter) ) THEN
    IF ( MOD(iter, rr_step) == 0 ) THEN
       !
       !
       ! unpack the active bands to the sparse order to compute eigenvalues 
       CALL reshape_array(psi, kdimx, nact, nbnd, act_idx, buffer1, -1) 
       CALL reshape_array(hpsi, kdimx, nact, nbnd, act_idx, buffer1, -1) 
       CALL reshape_array(w, kdimx, nact, nbnd, act_idx, buffer1, -1) 
       CALL reshape_array(hw, kdimx, nact, nbnd, act_idx, buffer1, -1) 
       CALL reshape_array(p, kdimx, nact, nbnd, act_idx, buffer1, -1) 
       CALL reshape_array(hp, kdimx, nact, nbnd, act_idx, buffer1, -1) 
       if(overlap) then 
         CALL reshape_array(spsi, kdimx, nact, nbnd, act_idx, buffer1, -1) 
         CALL reshape_array(sw, kdimx, nact, nbnd, act_idx, buffer1, -1) 
         CALL reshape_array(sp, kdimx, nact, nbnd, act_idx, buffer1, -1) 
       end if
       !
       call start_clock('ppcg:RR')
       if(overlap) then 
         CALL extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi, spsi )
       else
         CALL extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi )
       end if 
       call stop_clock('ppcg:RR')
       !
       IF (print_info >= 2) WRITE(stdout, *) 'RR has been invoked.' ; !CALL flush( stdout )
       !
       ! ... Compute the new residual vector block by evaluating
       !     residuals for individual eigenpairs in psi and e
       nblock = (kdim-1) / blocksz + 1         ! used to optimize some omp parallel do loops
       if (overlap) then
          !$omp parallel do collapse(2)
          DO j = 1, nbnd ; DO i=1,nblock
             w( 1+(i-1)*blocksz:MIN(i*blocksz,kdim) ,j ) = hpsi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j ) &
                                                  - spsi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j )*e( j )
          END DO ; END DO
          !$omp end parallel do
       else
          !$omp parallel do collapse(2)
          DO j = 1, nbnd ; DO i=1,nblock
             w( 1+(i-1)*blocksz:MIN(i*blocksz,kdim) ,j ) = hpsi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j ) &
                                                          -  psi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j )*e( j )
          END DO ; END DO
          !$omp end parallel do
       end if
       !
       ! ... Lock converged eigenpairs (set up act_idx and nact)
       call start_clock('ppcg:lock')
       nact_old = nact;
       CALL lock_epairs(kdim, nbnd, btype, w, kdimx, lock_tol, nact, act_idx)
       call stop_clock('ppcg:lock')
       !
       ! pack the active bands to the first columns (act_idx)
       CALL reshape_array(psi, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
       CALL reshape_array(hpsi, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
       CALL reshape_array(w, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
       CALL reshape_array(hw, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
       CALL reshape_array(p, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
       CALL reshape_array(hp, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
       if(overlap) then 
         CALL reshape_array(spsi, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
         CALL reshape_array(sw, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
         CALL reshape_array(sp, kdimx, nact, nbnd, act_idx, buffer1(1,1), 1) 
       end if 
       !
       ! ... Set up iteration parameters after locking
       CALL setup_param
       !
       !
       trG1  =  0.D0
       trdif = -1.D0
       trG   = SUM( e(act_idx(1:nact)) )
       !
    ELSE
       !
!! EV begin
      ! ... orthogonalize psi and update hpsi accordingly
      IF ( .NOT. force_repmat ) THEN
         !
         call start_clock('ppcg:cholQR')
         if (overlap) then
            CALL cholQR_dmat(kdim, nact, psi, spsi, kdimx, Gl, idesc)
         else
            CALL cholQR_dmat(kdim, nact, psi, psi, kdimx, Gl, idesc)
         end if
         call stop_clock('ppcg:cholQR')
         !
         call start_clock('ppcg:ZTRSM')
         CALL zgemm_dmat( kdim, nact, kdimx, idesc, C_ONE, hpsi, Gl, C_ZERO, hpsi )
         call stop_clock('ppcg:ZTRSM')
         !
         if (overlap) then
            call start_clock('ppcg:ZTRSM')
            CALL zgemm_dmat( kdim, nact, kdimx, idesc, C_ONE, spsi, Gl, C_ZERO, spsi )
            call stop_clock('ppcg:ZTRSM')
         end if
         !
      ELSE
         !
         call start_clock('ppcg:cholQR')
         if (overlap) then
            CALL cholQR(kdim, nact, psi, spsi, kdimx, G, nbnd)
         else
            CALL cholQR(kdim, nact, psi, psi, kdimx, G, nbnd)
         end if
         call stop_clock('ppcg:cholQR')
         !
         call start_clock('ppcg:ZTRSM')
         CALL ZTRSM('R', 'U', 'N', 'N', kdim, nact, C_ONE, G, nbnd, hpsi, kdimx)
         call stop_clock('ppcg:ZTRSM')
         !
         if (overlap) then
            call start_clock('ppcg:ZTRSM')
            CALL ZTRSM('R', 'U', 'N', 'N', kdim, nact, C_ONE, G, nbnd, spsi, kdimx)
            call stop_clock('ppcg:ZTRSM')
         end if
         !
      END IF
!! EV end
       !
       ! ... Compute the new subspace residual for active columns
       !
       !  G = psi'hpsi
       call start_clock('ppcg:zgemm')
       G = C_ZERO
       CALL divide(inter_bgrp_comm,nact,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nact,n_start,n_end
       if (n_start .le. n_end) &
       CALL ZGEMM('C','N', nact, my_n, kdim, C_ONE, psi, kdimx, hpsi(1,n_start), kdimx, C_ZERO, G(1,n_start), nbnd)
       CALL mp_sum(G(1:nact,1:nact), inter_bgrp_comm)
       CALL mp_sum(G(1:nact,1:nact), intra_bgrp_comm)
       call stop_clock('ppcg:zgemm')
       !
       ! w = hpsi - spsi*G
       call start_clock('ppcg:zgemm')
       IF ( my_bgrp_id /= root_bgrp_id ) then 
         call threaded_memset( w, 0.d0, 2*kdimx*nact )
       ELSE
         w(:,1:nact) = hpsi(:,1:nact)
       END IF 
       if (overlap) then
          if (n_start .le. n_end) &
          CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, spsi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, w, kdimx)
       else
          if (n_start .le. n_end) &
          CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, psi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, w, kdimx)
       end if
       CALL mp_sum( w(:,1:nact), inter_bgrp_comm )
       call stop_clock('ppcg:zgemm')
       !
       ! ... Compute trace of the projected matrix on current iteration
       ! trG1  = get_trace( G(1:nact, 1:nact), nact )
       trG1  = get_trace( G, nbnd, nact )
       trdif = ABS(trG1 - trG)
       trG   = trG1
       !
    END IF
    !
    ! Print iteration info ...
    IF (print_info >= 1) THEN
       WRITE(stdout, '("iter: ", I5, " nact = ", I5, ", trdif = ", 1pD9.2, ", trtol = ", 1pD9.2 )') iter, nact, trdif, trtol
       IF (print_info == 3) THEN
          res_array(iter) = print_rnrm( w, kdim, nbnd, kdimx)
          WRITE(stdout,'("Res. norm:  ", 1pD9.2)') res_array(iter)
       END IF
       FLUSH( stdout )
    END IF
    !
    total_iter = iter
    iter   = iter + 1
    !
    !
 END DO   !---End the main loop
 !
 ! unpack the (no-more-)active bands to the sparse order to compute eigenvalues and return the output psi
 CALL reshape_array(psi, kdimx, nact, nbnd, act_idx, buffer1, -1) 
!
! IF (nact > 0) THEN
 IF ( MOD(iter-1, rr_step) /= 0 ) THEN        ! if RR has not just been performed
 ! if nact==0 then the RR has just been performed
 ! in the main loop
    call start_clock('ppcg:RR')
    !
    ! unpack the (no-more-)active bands to the sparse order to compute eigenvalues and return the output 
    CALL reshape_array(hpsi, kdimx, nact, nbnd, act_idx, buffer1, -1) 
    CALL reshape_array(w, kdimx, nact, nbnd, act_idx, buffer1, -1) 
    if(overlap) then 
      CALL reshape_array(spsi, kdimx, nact, nbnd, act_idx, buffer1, -1) 
      CALL extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi, spsi )
    else
      CALL extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi )
    end if 
    call stop_clock('ppcg:RR')
    !
    ! ... Compute residuals
    if (overlap) then
       !$omp parallel do collapse(2)
       DO j = 1, nbnd ; DO i=1,nblock
          w( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j ) = hpsi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j ) &
                                               - spsi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j )*e( j )
       END DO ; END DO
       !$omp end parallel do
    else
       !$omp parallel do collapse(2)
       DO j = 1, nbnd ; DO i=1,nblock
          w( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j ) = hpsi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j ) &
                                               -  psi( 1+(i-1)*blocksz:MIN(i*blocksz,kdim), j )*e( j )
       END DO ; END DO
       !$omp end parallel do
    end if
    !
    ! ... Get the number of converged eigenpairs and their indices
    !     Note: The tolerance is 10*lock_tol, i.e., weaker tham lock_tol
! E.V. notconv issue should be addressed
    call start_clock('ppcg:lock')
    CALL lock_epairs(kdim, nbnd, btype, w, kdimx, 10*lock_tol, nact, act_idx)
    call stop_clock('ppcg:lock')
    !
 END IF
 !
! E.V. notconv issue comment
 notconv = 0 ! nact
 !
 IF (print_info >= 1) THEN
    WRITE(stdout, *) '-----------PPCG result summary ...  ----------------'
    WRITE(stdout, '("avg_iter: ", f6.2,  ", notconv: ", I5)') avg_iter, notconv
    FLUSH( stdout )
 END IF
 !
 CALL deallocate_all
 !
 CALL stop_clock( 'ppcg_k' )
 !
!!!EV-BANDS
if (print_info == 3) then
    write (stdout,'(1pD9.2)') ( res_array(j), j=1,maxter )
end if
 !
 RETURN
   !
CONTAINS
  !
  SUBROUTINE allocate_all
    !
    ! This subroutine allocates memory for the eigensolver
    !
    INTEGER :: nx
    ! maximum local block dimension
    !
    ALLOCATE ( buffer1(kdimx,sbsize), stat = ierr ) 
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate buffer1 ', ABS(ierr) )
    ALLOCATE ( hpsi(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate hpsi ', ABS(ierr) )
    if (overlap) ALLOCATE ( spsi(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate spsi ', ABS(ierr) )
    ALLOCATE ( w(kdimx,nbnd), hw(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate w and hw ', ABS(ierr) )
    if (overlap) ALLOCATE ( sw(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate sw ', ABS(ierr) )
    ALLOCATE ( p(kdimx,nbnd), hp(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate p and hp ', ABS(ierr) )
    if (overlap) ALLOCATE ( sp(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate sp ', ABS(ierr) )
    ALLOCATE ( K(sbsize3, sbsize3), M(sbsize3,sbsize3), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate K and M ', ABS(ierr) )
    ALLOCATE ( cwork( 1 + 18*sbsize + 18*sbsize**2 ), rwork( 1 + 18*sbsize + 18*sbsize**2 ), iwork(3 + 15*sbsize), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate lapack work arrays ', ABS(ierr) )
    !
    CALL desc_init( nbnd, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip  )
    !
    IF ( la_proc ) THEN
       ALLOCATE( Gl( nx, nx ), STAT=ierr )
    ELSE
       ALLOCATE( Gl( 1, 1 ), STAT=ierr )
    END IF
    IF( ierr /= 0 ) CALL errore( 'ppcg ',' cannot allocate Gl ', ABS(ierr) )

  END SUBROUTINE allocate_all
  !
  !
  !
  SUBROUTINE cholQR(n, k, X, SX, ldx,  R, ldr)
    !
    !  This subroutine orthogonalizes X using the Choleski decomposition.
    !  The matrix X'X and its Choleski decomposition is replicated on all processors.
    !  TBD: If the Choleslki orthogonalization fails, a slower but more accuarte
    !  Householder QR is performed.
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    INTEGER,     INTENT (IN) :: n, k, ldx, ldr
    COMPLEX(DP), INTENT (INOUT) :: X(ldx,*), SX(ldx,*)
    COMPLEX(DP), INTENT(OUT) :: R(ldr,*)
    !
    ! ... local variables
    !
    COMPLEX(DP):: XTX(k,k)
    INTEGER    :: ierr = 0, info = 0
    !
    ! ... do Cholesky QR unless X is rank deficient
    !
    CALL ZGEMM('C','N', k, k, n, C_ONE, X, ldx, SX, ldx, C_ZERO, XTX, k)
    !
    CALL mp_sum( XTX, intra_bgrp_comm )
    !
    CALL ZPOTRF('U', k, XTX, k, info)
    IF ( info == 0 ) THEN
       !
       ! ... triangualar solve
       CALL ZTRSM('R', 'U', 'N', 'N', n, k, C_ONE, XTX, k, X, ldx)
       !
    ELSE
       ! TBD: QR
       WRITE(stdout,*) '[Q, R] = qr(X, 0) failed'
       STOP
!           ![X, R] = qr(X, 0)
!           ! QR without forming Q
!           !
!           call zgeqrf(nn,kk,X,nn,tau,wqr,lwqr,info)
!           !
!           if (info<0) then
!              write(*,*), '[Q, R] = qr(X, 0) failed'
!              stop
!           endif
!           ! QR: constructing Q
!           !
!           call zungqr(nn,kk,kk,X,nn,tau,wqr,lwqr,info)
!           !
    END IF
    !
    ! ... also return R factor if needed
    !
    CALL ZLACPY('U', k, k,  XTX, k,  R, ldr)
    !
    RETURN
    !
  END SUBROUTINE cholQR
  !
  !
  !
  SUBROUTINE cholQR_dmat(kdim, k, X, SX, kdimx,  Rl, idesc)
    !
    ! Distributed version of cholQR
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    INTEGER,     INTENT (IN) :: kdim, k, kdimx
    COMPLEX(DP), INTENT (INOUT) :: X(kdimx,k)
    COMPLEX(DP), INTENT (IN) :: SX(kdimx,k)
    INTEGER, INTENT (IN)  :: idesc(:)
    COMPLEX(DP), INTENT(OUT) :: Rl(:, :)
    ! inverse of the upper triangular Cholesky factor
    !
    ! ... local variables
    !
    COMPLEX(DP)  :: buffer(kdimx,k)
    COMPLEX(DP), ALLOCATABLE   :: XTXl(:,:)
    INTEGER    :: nx, ierr = 0
#ifdef __SCALAPACK
    INTEGER     :: desc_sca( 16 ), info
#endif
    !
    nx = idesc(LAX_DESC_NRCX)
    !
    IF ( la_proc ) THEN
       !
       ALLOCATE( XTXl( nx, nx ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate XTXl ', ABS(ierr) )
       !
    ELSE
       !
       ALLOCATE( XTXl( 1, 1 ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate XTXl ', ABS(ierr) )
       !
    END IF
    !
    ! ... Perform Cholesky of X'X
    !
    CALL compute_distmat(XTXl, idesc, X, SX, k)
    !
    IF ( la_proc ) THEN
       !
#ifdef __SCALAPACK
       CALL descinit( desc_sca, k, k, nx, nx, 0, 0, ortho_cntx, SIZE( XTXl, 1 ), info )
       IF( info /= 0 ) CALL errore( ' ppcg ', ' descinit ', ABS( info ) )
       !
       CALL PZPOTRF( 'U', k, XTXl, 1, 1, desc_sca, info )
!      IF( info /= 0 ) CALL errore( ' ppcg ', ' problems computing cholesky ', ABS( info ) )
       !
       IF ( info == 0 ) THEN
          !
!          ! set the lower triangular part to zero
!          CALL sqr_zsetmat( 'L', k, C_ZERO, XTXl, size(XTXl,1), desc )
          !
          ! find inverse of the upper triangular Cholesky factor R
          CALL PZTRTRI( 'U', 'N', k, XTXl, 1, 1, desc_sca, info )
          IF( info /= 0 ) CALL errore( ' ppcg ', ' problems computing inverse ', ABS( info ) )
          !
          ! set the lower triangular part to zero
          CALL sqr_setmat( 'L', k, C_ZERO, XTXl, size(XTXl,1), idesc )
       !
       ELSE
       ! TBD: QR
          WRITE(stdout,*) '[Q, R] = qr(X, 0) failed'
          STOP
!           ![X, R] = qr(X, 0)
!           ! QR without forming Q
!           !
!           call zgeqrf(nn,kk,X,nn,tau,wqr,lwqr,info)
!           !
!           if (info<0) then
!              write(*,*), '[Q, R] = qr(X, 0) failed'
!              stop
!           endif
!           ! QR: constructing Q
!           !
!           call zungqr(nn,kk,kk,X,nn,tau,wqr,lwqr,info)
!           !
       END IF
#else
       CALL laxlib_pzpotrf( XTXl, nx, k, idesc )
       !
       CALL laxlib_pztrtri ( XTXl, nx, k, idesc )
#endif
    !
    !
    END IF
    !
    CALL zgemm_dmat( kdim, k, kdimx, idesc, C_ONE, X, XTXl, C_ZERO, buffer )
    !
    X = buffer
    ! ... also return R factor
    Rl = XTXl
    !
    DEALLOCATE(XTXl)
    !
    RETURN
    !
  END SUBROUTINE cholQR_dmat
  !
  !
  !
  SUBROUTINE lock_epairs(kdim, nbnd, btype, w, kdimx, tol, nact, act_idx)
     !
     ! Lock converged eigenpairs: detect "active" columns of w
     ! by checking if each column has norm greater than tol.
     ! Returns the number of active columns and their indices.
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER,     INTENT (IN) :: kdim, kdimx, nbnd, btype(nbnd)
     COMPLEX(DP), INTENT (IN) :: w(kdimx,nbnd)
     REAL(DP),    INTENT (IN) :: tol
     INTEGER,     INTENT(OUT) :: nact, act_idx(nbnd)
     !
     ! ... local variables
     !
     INTEGER         :: j
     REAL(DP)        :: rnrm_store(nbnd), band_tollerance 
     REAL(DP), EXTERNAL :: DDOT

     !
     nact = 0
     ! ... Compute norms of each column of psi
     rnrm_store = 0.D0
     CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
     DO j = n_start, n_end
        !
        rnrm_store(j)   =  DDOT(2*kdim, w(:,j), 1, w(:,j), 1)
        !
     END DO
     CALL mp_sum( rnrm_store, inter_bgrp_comm )
     !
     CALL mp_sum( rnrm_store, intra_bgrp_comm )
     !
     DO j = 1, nbnd
        !
        if ( btype(j) == 0 ) then
             band_tollerance = max(2.5*tol,1.d-3)
        else
             band_tollerance = tol
        end if
        !
        rnrm_store(j) = SQRT( rnrm_store(j) )
        !
        IF ( (print_info >= 2) .AND. (iter > 1) )  THEN
          write(stdout, '( "Eigenvalue ", I5, " = ", 1pe12.4, ". Residual norm = ",  1pe9.2)') &
                      j, e(j), rnrm_store(j)
        END IF
        !
        IF ( rnrm_store(j) > band_tollerance ) THEN
           nact = nact + 1
           act_idx(nact) = j
        END IF
        !
     END DO
     !
     !
  END SUBROUTINE lock_epairs
  !
  !
  SUBROUTINE setup_param
     !
     ! Based on the information about active columns of psi,
     ! set up iteration parameters, such as:
     !   - number of subblock
     !   - size of the last sub-block
     !   - tolerance level for the trace difference of reduced hamiltonians
     !     on consecutive iterations
     !   - replicate or re-distribute the Gram matrix G
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
!     INTEGER,     INTENT (IN) :: nbnd, nact, act_idx(nbnd), sbsize
!     INTEGER,     INTENT(OUT) :: nsb, sbsize_last
!     REAL(DP),    INTENT(OUT) :: tol
     !
     ! ... local variables
     !
     INTEGER         :: ierr
     INTEGER :: nx, act_thresh
     ! maximum local block dimension
     ! Threshold on the number of active vectors that determines if
     ! Choleski is distributed or replicated
     !
     !
     ! ... Compute the number of sub-blocks and the size of the last sub-block
     sbsize_last = sbsize
     nsb   =  FLOOR( DBLE(nact)/DBLE(sbsize) )
     IF ( MOD( nact, sbsize ) /= 0) THEN
        sbsize_last =  nact - sbsize*nsb
        nsb         =  nsb + 1
     END IF
     !
     ! ... Compute the current tolerance level for nact pairs
     trtol   = ethr*SQRT(DBLE(nact)) ! MIN( ethr*SQRT(DBLE(nact)), 1.0D-2 )
     !
     ! ... Redistribute matrices for Cholesky because number of active vectors has changed
     !
     ! Don't want to run distributed Cholesky for matrices of size <= 100
     act_thresh = MAX( 100, np_ortho(1) )
     !
     IF ( ( nact > act_thresh ) .AND. (nact /= nact_old) ) THEN
        !
        IF ( ALLOCATED(Gl) ) DEALLOCATE(Gl)
        !
        CALL desc_init( nact, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip  )
        !
        IF ( la_proc ) THEN
           !
           ALLOCATE( Gl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'ppcg ',' cannot allocate Gl ', ABS(ierr) )
           !
        ELSE
           !
           ALLOCATE( Gl( 1, 1 ), STAT=ierr )
           IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate Gl ', ABS(ierr) )
           !
        END IF
        !
        force_repmat = .FALSE.
        !
     ELSE
       ! replicate Cholesky because number of active vectors is small
       !
       IF (nact <= act_thresh) THEN
          !
          force_repmat = .TRUE.
          !
          IF ( ALLOCATED(Gl) ) DEALLOCATE(Gl)
          !
        ELSE
          !
          force_repmat = .FALSE.
          !
        END IF
           !
     END IF
     !
     !
     ! ... Allocate storage for the small matrices in separate minimization loop
     IF ( ALLOCATED(K_store) ) DEALLOCATE(K_store)
     IF ( ALLOCATED(M_store) ) DEALLOCATE(M_store)
     ALLOCATE ( K_store(sbsize3, sbsize3*nsb), M_store(sbsize3,sbsize3*nsb), stat = ierr )
     IF (ierr /= 0) &
        CALL errore( 'ppcg ',' cannot allocate K_store and M_store ', ABS(ierr) )
     !
     RETURN
     !
  END SUBROUTINE setup_param
  !
  !
  !
  !
  SUBROUTINE extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi, spsi)
     !
     ! Perform the Rayleigh-Ritz to "rotate" psi to eigenvectors
     ! The Gram matrices are distributed over la processor group
     !
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER,     INTENT (IN)  :: kdim, kdimx, nbnd
     REAL(DP),    INTENT (OUT) :: e(nbnd)
     COMPLEX(DP), INTENT (INOUT)  :: psi(kdimx,nbnd), hpsi(kdimx,nbnd)
     COMPLEX(DP), INTENT (INOUT), OPTIONAL  :: spsi(kdimx,nbnd)
!     LOGICAL,     INTENT(IN), OPTIONAL :: ortho
     ! ortho = .true.(default) then orthogonalization of psi prior to the RR is enforced
     !
     ! ... local variables
     !
     COMPLEX (DP), ALLOCATABLE    :: Hl(:,:), Sl(:,:)
     ! local part of projected Hamiltonian and of the overlap matrix
     INTEGER :: idesc(LAX_DESC_SIZE)
     !
     ! Matrix distribution descriptors to temporary store the "global" current descriptor
     LOGICAL :: la_proc_store
     ! flag to distinguish procs involved in linear algebra
     INTEGER, ALLOCATABLE :: irc_ip_store( : )
     INTEGER, ALLOCATABLE :: nrc_ip_store( : )
     INTEGER, ALLOCATABLE :: rank_ip_store( :, : )
     !
     COMPLEX(DP), ALLOCATABLE  :: vl(:,:)
     COMPLEX(DP), ALLOCATABLE  :: buffer2(:,:)
!     REAL(DP)                  :: R(nbnd, nbnd)
     INTEGER                   :: nx
!     LOGICAL                   :: do_orth
     !
     if (.not. overlap .and. present(spsi) .or. &
         overlap .and. .not. present(spsi) ) Call errore('ppcg : extract_epairs_dmat', 'wrong overlap', 1) 
     !
     ALLOCATE ( buffer2(kdimx,nbnd), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate buffer2 ', ABS(ierr) )
     !
     ! Store current distributed matrix descriptor information
     ALLOCATE ( irc_ip_store( np_ortho(1)  ), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate irc_ip_store ', ABS(ierr) )
     !
     ALLOCATE ( nrc_ip_store( np_ortho(1) ), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate nrc_ip_store ', ABS(ierr) )
     !
     ALLOCATE ( rank_ip_store( np_ortho(1), np_ortho(2) ), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate rank_ip_store ', ABS(ierr) )
     !
     irc_ip_store  = irc_ip
     nrc_ip_store  = nrc_ip
     rank_ip_store = rank_ip
     !
     CALL desc_init( nbnd, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip  )
     !
     IF ( la_proc ) THEN
        !
        ALLOCATE( vl( nx, nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate vl ', ABS(ierr) )
        !
        ALLOCATE( Sl( nx, nx ), STAT=ierr )
        IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate Sl ', ABS(ierr) )
        !
        ALLOCATE( Hl( nx, nx ), STAT=ierr )
        IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate Hl ', ABS(ierr) )
        !
     ELSE
        !
        ALLOCATE( vl( 1, 1 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
        !
        ALLOCATE( Sl( 1, 1 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate Sl ', ABS(ierr) )
        !
        ALLOCATE( Hl( 1, 1 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate Hl ', ABS(ierr) )
        !
     END IF
     !
     !
     !  G = psi'*hpsi
     CALL compute_distmat(Hl, idesc, psi, hpsi, nbnd)
     if (overlap) then
        CALL compute_distmat(Sl, idesc, psi, spsi, nbnd)
     else
        CALL compute_distmat(Sl, idesc, psi,  psi, nbnd)
     end if
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg e and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id) CALL pdiaghg( nbnd, Hl, Sl, nx, e, vl, idesc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
          CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
          CALL mp_bcast( e,  root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pdiaghg( nbnd, Hl, Sl, nx, e, vl, idesc )
     END IF
     !
     ! "Rotate" psi to eigenvectors
     !
     CALL zgemm_dmat( kdim, nbnd, kdimx, idesc, C_ONE,  psi, vl, C_ZERO,  buffer2 )
     psi = buffer2
     CALL zgemm_dmat( kdim, nbnd, kdimx, idesc, C_ONE, hpsi, vl, C_ZERO, buffer2 )
     hpsi = buffer2
     if (overlap) then 
       CALL zgemm_dmat( kdim, nbnd, kdimx, idesc, C_ONE, spsi, vl, C_ZERO, buffer2 )
       spsi = buffer2
     end if 
     !
     ! Restore current "global" distributed  matrix descriptor
     !
     irc_ip  = irc_ip_store
     nrc_ip  = nrc_ip_store
     rank_ip = rank_ip_store
     !
     DEALLOCATE ( irc_ip_store, nrc_ip_store, rank_ip_store )
     DEALLOCATE( buffer2 ) 
     DEALLOCATE(Hl, Sl)
     DEALLOCATE(vl)
     !
  END SUBROUTINE extract_epairs_dmat
  !
  !
  !
  REAL(DP) FUNCTION get_trace(G, ld, k)
     !
     !  This function returns trace of a k-by-k matrix G
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER,     INTENT (IN) :: ld, k
     COMPLEX(DP), INTENT (IN) :: G(ld,*)
     !
     ! ... local variables
     !
     INTEGER    :: j
     !
     get_trace = 0.D0
     DO j = 1, k
        get_trace = get_trace + DBLE(G(j, j))
     END DO
     !
     RETURN
     !
  END FUNCTION get_trace
  !
  !
  SUBROUTINE deallocate_all
    !
    ! This subroutine releases the allocated memory
    !
    IF ( ALLOCATED(buffer1) )    DEALLOCATE ( buffer1 )
    IF ( ALLOCATED(hpsi) )    DEALLOCATE ( hpsi )
    IF ( ALLOCATED(spsi) )    DEALLOCATE ( spsi )
    IF ( ALLOCATED(w) )       DEALLOCATE ( w )
    IF ( ALLOCATED(hw) )      DEALLOCATE ( hw )
    IF ( ALLOCATED(sw) )      DEALLOCATE ( sw )
    IF ( ALLOCATED(p) )       DEALLOCATE ( p )
    IF ( ALLOCATED(hp) )      DEALLOCATE ( hp )
    IF ( ALLOCATED(sp) )      DEALLOCATE ( sp )
    IF ( ALLOCATED(K) )       DEALLOCATE ( K )
    IF ( ALLOCATED(M) )       DEALLOCATE ( M )
    IF ( ALLOCATED(K_store) ) DEALLOCATE ( K_store )
    IF ( ALLOCATED(M_store) ) DEALLOCATE ( M_store )
    IF ( ALLOCATED(cwork) )   DEALLOCATE ( cwork )
    IF ( ALLOCATED(rwork) )   DEALLOCATE ( rwork )
    IF ( ALLOCATED(iwork) )   DEALLOCATE ( iwork )
    IF ( ALLOCATED(irc_ip) )  DEALLOCATE ( irc_ip )
    IF ( ALLOCATED(nrc_ip) )  DEALLOCATE ( nrc_ip )
    IF ( ALLOCATED(rank_ip) ) DEALLOCATE ( rank_ip )
    IF ( ALLOCATED(Gl) )      DEALLOCATE ( Gl )
!    IF ( ALLOCATED(G) )       DEALLOCATE( G )
!    IF ( ALLOCATED(Sl) )      DEALLOCATE( Sl )
    !
  END SUBROUTINE deallocate_all
  !
  !
  SUBROUTINE compute_distmat( dm, idesc, v, w, k)
!    Copy-paste from pcegterg and desc added as a parameter
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     COMPLEX(DP), INTENT(OUT)   :: dm( :, : )
     INTEGER, INTENT(IN) :: idesc(:)
     COMPLEX(DP), INTENT(IN) :: v(:,:), w(:,:)
     INTEGER, INTENT(IN)     :: k
     ! global size of dm = number of vectors in v and w blocks
     !
     ! ... local variables
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     INTEGER :: nx
     ! maximum local block dimension
     !
     nx = idesc(LAX_DESC_NRCX)
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = C_ZERO
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C','N', nr, nc, kdim, C_ONE, v(1,ir), kdimx, w(1,ic), kdimx, C_ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     !
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL laxlib_zsqmher( k, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE zgemm_dmat( n, k, ld, idesc, alpha, X, Gl, beta, Y  )
!    Copy-paste from refresh_evc in pregterg with some modifications
     !
     ! Compute Y = alpha*(X*G) + beta*Y, where G is distributed across la processor group
     ! and is given by local block Gl.
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER, INTENT(IN) :: n, k, ld
     ! number of rows of X and Y
     ! number of columns of X,Y  and size/leading dimension of (global) G
     ! leading dimension of X and Y
     INTEGER, INTENT(IN) :: idesc(:)
     ! descriptor of G
     COMPLEX(DP), INTENT(IN)      ::  alpha, beta
     COMPLEX(DP), INTENT (IN)     ::  X(ld, k)
     COMPLEX(DP), INTENT (INOUT)  ::  Y(ld, k)
     COMPLEX(DP), INTENT(IN)      ::  Gl( :, :)
     !
     ! ... local variables
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: Gltmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: Xtmp( :, : )
     COMPLEX(DP) :: gamm
     INTEGER :: nx
     !
     nx = idesc(LAX_DESC_NRCX)
     !
     ALLOCATE( Gltmp( nx, nx ) )
     ALLOCATE( Xtmp( ld, k ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= k ) THEN
           !
           nc = min( nc, k - ic + 1 )
           !
           gamm = C_ZERO

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )
              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( Gl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N','N', n, nc, nr, C_ONE, X(1,ir), ld, Gl, nx, gamm, Xtmp(1,ic), ld )
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( Gltmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N','N', n, nc, nr, C_ONE, X(1,ir), ld, Gltmp, nx, gamm, Xtmp(1,ic), ld )
              END IF
              !
              gamm = C_ONE
              !
           END DO
           !
        END IF
        !
     END DO
     !
     IF (beta /= 0.D0) THEN
        Y = alpha*Xtmp + beta*Y
     ELSE
       Y = alpha*Xtmp
     END IF
     !
     DEALLOCATE( Gltmp )
     DEALLOCATE( Xtmp )
     !
     RETURN
     !
  END SUBROUTINE zgemm_dmat
  !
  REAL(DP) FUNCTION print_rnrm(w, n, m, nx)
     !
     !  Compute the Frobenius norm of w 
     !
     INTEGER :: n, m, nx
     COMPLEX(DP), INTENT(IN) :: w(nx, m) 
     REAL(DP) :: rnrm   
     REAL(DP), EXTERNAL       :: ZLANGE
     ! 
     rnrm = 0.0d0
     rnrm = ZLANGE('F', n, m, w, nx, w)  ! work is not referenced for Frobenius norm
     rnrm = rnrm**2
     CALL mp_sum( rnrm, intra_bgrp_comm )
     rnrm = SQRT(ABS(rnrm))
     print_rnrm = rnrm
     !   
     RETURN
     !
  END FUNCTION print_rnrm
  !
  SUBROUTINE reshape_array(psi, kdimx, nact, nbnd, act_idx, buffer, isgn)
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: kdimx, nact, nbnd, act_idx(nbnd), isgn
    COMPLEX(DP), INTENT(INOUT) :: psi(kdimx, nbnd)  
    COMPLEX(DP), INTENT(INOUT) :: buffer(kdimx)  
    INTEGER :: ibnd
    !
    if(isgn.gt.0) then 
      do ibnd = 1, nact
        if(ibnd.ne.act_idx(ibnd)) then 
          buffer(1:kdimx) = psi(1:kdimx,ibnd) 
          psi(1:kdimx,ibnd) = psi(1:kdimx,act_idx(ibnd))
          psi(1:kdimx,act_idx(ibnd)) = buffer(1:kdimx) 
        end if 
      end do 
    elseif(isgn.lt.0) then 
      do ibnd = 1, nact 
        if(ibnd.ne.act_idx(ibnd)) then 
          buffer(1:kdimx) = psi(1:kdimx,act_idx(ibnd)) 
          psi(1:kdimx,act_idx(ibnd)) = psi(1:kdimx,ibnd)
          psi(1:kdimx,ibnd) = buffer(1:kdimx) 
        end if 
      end do 
    else
      Call errore("reshape_array", "wrong isgn", 1)
    endif 
    !
    RETURN  
    !
  END SUBROUTINE reshape_array
  !
  REAL(DP) FUNCTION get_trace2 (x, n, m) 
    IMPLICIT NONE 
    INTEGER :: i, j, n, m 
    COMPLEX(DP) :: x(n,m)
    get_trace2 = 0.0d0
    DO i = 1, n 
      DO j = 1, m
        get_trace2 = get_trace2 + x(i,j) * conjg( x(i,j) ) 
      END DO 
    END DO 
    RETURN
  END FUNCTION get_trace2  
  ! 
END SUBROUTINE ppcg_k_idx
