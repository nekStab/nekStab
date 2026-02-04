# Modal Analysis Implementation Plan for nekStab

## Overview

Integrate POD, DMD, and SPOD natively into nekStab using Nek5000 data structures and MPI patterns.

| Method | Output | Best For |
|--------|--------|----------|
| **POD** | Spatial modes ranked by energy | Dominant structures, ROMs |
| **DMD** | Modes with growth rates + frequencies | Stability, transient dynamics |
| **SPOD** | Per-frequency energy-ranked modes | Stationary turbulence |

**Design Principles:**
- Use `krylov_vector` type (separate vx, vy, vz, pr, t arrays)
- Use existing `k_dot`, `k_norm`, `k_copy`, `k_add2s2` operations
- Use `glsc3()` for weighted inner products (MPI-safe)
- Use `load_fld()`, `outpost2()` for I/O
- Match nekStab coding style (fixed-form compatible)

---

## nekStab Data Structures

### krylov_vector Type (existing)
```fortran
type :: krylov_vector
   real, dimension(lv) :: vx, vy, vz     ! Velocity (lx1×ly1×lz1×lelv)
   real, dimension(lp) :: pr             ! Pressure
   real, dimension(lv, ldimt) :: t       ! Temperature/scalars
   real :: time                          ! Time stamp
end type krylov_vector
```

### Snapshot Storage
```fortran
! Array of krylov_vectors (like in krylov_schur)
type(krylov_vector), allocatable :: snaps(:)  ! (nsnap)
allocate(snaps(nsnap))
```

### Inner Product (existing k_dot)
```fortran
! From krylov_subspace.f90 - uses glsc3 with bm1s mass matrix
call k_dot(alpha, p, q)  ! alpha = <p,q>_E (energy norm)

! Internally does:
alpha = glsc3(p%vx, q%vx, bm1s, nv) + glsc3(p%vy, q%vy, bm1s, nv)
if (if3d) alpha = alpha + glsc3(p%vz, q%vz, bm1s, nv)
if (ifto) alpha = alpha + glsc3(p%t(:,1), q%t(:,1), bm1s, nt)
```

---

## Mathematical Foundation

### POD
```
C(i,j) = <snaps(i), snaps(j)>_E / nsnap    (nsnap × nsnap correlation)
C v = λ v                                   (eigenvalue problem)
Φ_k = Σ_i v_k(i) snaps(i) / √(λ_k nsnap)   (spatial modes)
```

### DMD (Projected)
```
Q_1 = [snaps(1), ..., snaps(n-1)]
Q_2 = [snaps(2), ..., snaps(n)]
Q_1 ≈ U Σ V^T                    (truncated SVD via eigendecomposition)
Ã = U^H Q_2 V Σ^{-1}             (r × r projected operator)
Ã w = μ w                        (eigenvalues)
Φ = U w                          (projected DMD modes)
ω = log(μ)/dt                    (continuous eigenvalue)
```

### SPOD
```
n_blk = (nsnap - n_fft) / (n_fft - n_overlap) + 1

For each frequency k:
  Q̂_k = FFT of windowed blocks at frequency k    (complex, n_blk vectors)
  CSD(i,j) = <Q̂_k(i), Q̂_k(j)>_E / n_blk         (Hermitian)
  CSD v = λ v                                     (eigenvalue problem)
  Ψ = Σ_i v(i) Q̂_k(i) / √(λ n_blk)              (SPOD modes)
```

---

## User Parameters

### NEKSTAB.inc Declarations (to add)
```fortran
!     Modal analysis parameters (Mode 6)
      character(len=80) modal_prefix
      integer modal_nsnap, modal_nsave, dmd_rank, spod_nfft, spod_noverlap
      real modal_dt
      logical ifpod, ifdmd, ifspod

      common /modal_c/ modal_prefix
      common /modal_i/ modal_nsnap, modal_nsave, dmd_rank, spod_nfft, spod_noverlap
      common /modal_r/ modal_dt
      common /modal_l/ ifpod, ifdmd, ifspod
```

### Set in `nekStab_usrchk` (`.usr` file):

```fortran
subroutine nekStab_usrchk
   implicit none
   include 'SIZE'
   include 'TOTAL'

!  ═══════════════════════════════════════════════════════════════════
!  MODAL ANALYSIS (Mode 6)
!  ═══════════════════════════════════════════════════════════════════

!  Snapshot source
   modal_prefix = 'DNS'     ! File prefix (e.g., DNS_cyl0.f00001)
   modal_nsnap  = 500       ! Number of snapshots
   modal_dt     = 0.1       ! Time between snapshots
   modal_nsave  = 10        ! Modes to save

!  Which methods to run
   ifpod  = .true.          ! Compute POD
   ifdmd  = .true.          ! Compute DMD
   ifspod = .true.          ! Compute SPOD

!  DMD options
   dmd_rank = 0             ! SVD rank (0 = auto 99% energy)

!  SPOD options
   spod_nfft = 128          ! FFT block length
   spod_noverlap = 64       ! Overlap (typically nfft/2)

   return
end subroutine
```

---

## Algorithm (Mode 6)

```fortran
subroutine modal_analysis
   use krylov_subspace
   implicit none
   include 'SIZE'
   include 'TOTAL'

   type(krylov_vector), allocatable :: snaps(:)
   type(krylov_vector) :: mean_snap
   real, allocatable :: C(:,:), eigvals(:), eigvecs(:,:)
   integer :: i, j

!  ─────────────────────────────────────────────────────────────────
!  PHASE 1: Load snapshots
!  ─────────────────────────────────────────────────────────────────
   allocate(snaps(modal_nsnap))
   call load_files(snaps, modal_nsnap, modal_nsnap, modal_prefix)

!  Compute and subtract mean
   call compute_mean(snaps, modal_nsnap, mean_snap)
   call subtract_mean(snaps, modal_nsnap, mean_snap)

!  ─────────────────────────────────────────────────────────────────
!  PHASE 2: POD
!  ─────────────────────────────────────────────────────────────────
   if (ifpod) then
      allocate(C(modal_nsnap, modal_nsnap))
      allocate(eigvals(modal_nsnap), eigvecs(modal_nsnap, modal_nsnap))

      ! Form correlation matrix using k_dot
      do j = 1, modal_nsnap
         do i = 1, j
            call k_dot(C(i,j), snaps(i), snaps(j))
            C(j,i) = C(i,j)  ! Symmetric
         end do
      end do
      C = C / modal_nsnap

      ! Eigensolve (DSYEV)
      call eig_symmetric(C, eigvals, eigvecs, modal_nsnap)

      ! Reconstruct and save modes
      call pod_reconstruct_modes(snaps, eigvecs, eigvals, modal_nsnap)

      deallocate(C, eigvals, eigvecs)
   end if

!  ─────────────────────────────────────────────────────────────────
!  PHASE 3: DMD
!  ─────────────────────────────────────────────────────────────────
   if (ifdmd) then
      call dmd_compute(snaps, modal_nsnap, modal_dt, dmd_rank)
   end if

!  ─────────────────────────────────────────────────────────────────
!  PHASE 4: SPOD
!  ─────────────────────────────────────────────────────────────────
   if (ifspod) then
      call spod_compute(snaps, modal_nsnap, modal_dt, spod_nfft, spod_noverlap)
   end if

   deallocate(snaps)
end subroutine
```

---

## Key Subroutines

### Mean Computation
```fortran
subroutine compute_mean(snaps, nsnap, mean_snap)
   use krylov_subspace
   implicit none
   include 'SIZE'
   include 'TOTAL'

   integer, intent(in) :: nsnap
   type(krylov_vector), intent(in) :: snaps(nsnap)
   type(krylov_vector), intent(out) :: mean_snap
   integer :: i
   real :: scale

   call k_zero(mean_snap)
   do i = 1, nsnap
      call k_add2(mean_snap, snaps(i))
   end do
   scale = 1.0d0 / nsnap
   call k_cmult(mean_snap, scale)
end subroutine
```

### Mean Subtraction
```fortran
subroutine subtract_mean(snaps, nsnap, mean_snap)
   use krylov_subspace
   implicit none
   include 'SIZE'
   include 'TOTAL'

   integer, intent(in) :: nsnap
   type(krylov_vector), intent(inout) :: snaps(nsnap)
   type(krylov_vector), intent(in) :: mean_snap
   integer :: i

   do i = 1, nsnap
      call k_sub2(snaps(i), mean_snap)  ! snaps(i) = snaps(i) - mean
   end do
end subroutine
```

### POD Mode Reconstruction
```fortran
subroutine pod_reconstruct_modes(snaps, eigvecs, eigvals, nsnap)
   use krylov_subspace
   implicit none
   include 'SIZE'
   include 'TOTAL'

   integer, intent(in) :: nsnap
   type(krylov_vector), intent(in) :: snaps(nsnap)
   real, intent(in) :: eigvecs(nsnap, nsnap), eigvals(nsnap)

   type(krylov_vector) :: mode
   real :: scale
   integer :: m, i
   character(len=3) :: prefix

   prefix = 'pod'

!  Modes are sorted descending by eigenvalue (done in eig_symmetric)
   do m = 1, min(modal_nsave, nsnap)
      call k_zero(mode)

      ! Φ_m = Σ_i v_m(i) * snaps(i) / sqrt(λ_m * nsnap)
      do i = 1, nsnap
         call k_add2s2(mode, snaps(i), eigvecs(i, m))
      end do

      if (eigvals(m) > 0) then
         scale = 1.0d0 / sqrt(eigvals(m) * nsnap)
         call k_cmult(mode, scale)
      end if

      ! Output mode
      call nopcopy(vx, vy, vz, pr, t, mode%vx, mode%vy, mode%vz, mode%pr, mode%t)
      call outpost2(vx, vy, vz, pr, t, 0, prefix)
   end do

!  Write eigenvalue spectrum
   if (nid == 0) then
      open(unit=77, file='pod_energy.dat', status='replace')
      write(77, '(A)') '# mode  eigenvalue  cumulative_%'
      do m = 1, nsnap
         write(77, '(I6, 2E16.8)') m, eigvals(m), &
              100.0 * sum(eigvals(1:m)) / sum(eigvals)
      end do
      close(77)
   end if
end subroutine
```

### Eigenvalue Solver Wrapper
```fortran
subroutine eig_symmetric(A, eigvals, eigvecs, n)
!  Wrapper for LAPACK DSYEV
!  Returns eigenvalues in DESCENDING order (largest first)

   implicit none
   integer, intent(in) :: n
   real, intent(in) :: A(n, n)
   real, intent(out) :: eigvals(n), eigvecs(n, n)

   real, allocatable :: work(:), Acopy(:,:)
   real :: tmp_val, tmp_vec(n)
   integer :: lwork, info, i, j

   allocate(Acopy(n,n))
   Acopy = A

!  Query optimal workspace
   allocate(work(1))
   call dsyev('V', 'U', n, Acopy, n, eigvals, work, -1, info)
   lwork = int(work(1))
   deallocate(work)
   allocate(work(lwork))

!  Compute eigenvalues/vectors (ascending order)
   call dsyev('V', 'U', n, Acopy, n, eigvals, work, lwork, info)
   eigvecs = Acopy

!  Reverse to descending order
   do i = 1, n/2
      j = n - i + 1
      tmp_val = eigvals(i)
      eigvals(i) = eigvals(j)
      eigvals(j) = tmp_val
      tmp_vec = eigvecs(:, i)
      eigvecs(:, i) = eigvecs(:, j)
      eigvecs(:, j) = tmp_vec
   end do

   deallocate(work, Acopy)
end subroutine
```

---

## Complex Fields for DMD/SPOD

**Problem**: DMD and SPOD modes are complex, but Nek5000 fields are real.

**Solution**: Store real and imaginary parts separately:

```fortran
! For complex mode = mode_re + i * mode_im
type(krylov_vector) :: mode_re, mode_im

! Output as separate files
call nopcopy(vx,vy,vz,pr,t, mode_re%vx,mode_re%vy,mode_re%vz,mode_re%pr,mode_re%t)
call outpost2(vx, vy, vz, pr, t, 0, 'dRe')  ! Real part

call nopcopy(vx,vy,vz,pr,t, mode_im%vx,mode_im%vy,mode_im%vz,mode_im%pr,mode_im%t)
call outpost2(vx, vy, vz, pr, t, 0, 'dIm')  ! Imaginary part
```

**Output naming**:
- POD: `pod_<SESSION>0.f*` (real modes)
- DMD: `dRe_<SESSION>0.f*`, `dIm_<SESSION>0.f*`
- SPOD: `sRe_<SESSION>0.f*`, `sIm_<SESSION>0.f*`

---

## Complex Inner Product for SPOD

For SPOD, we need Hermitian inner product of complex fields:

```fortran
subroutine k_dot_complex(alpha, p_re, p_im, q_re, q_im)
!  Compute <p, q> = p^H * W * q for complex krylov_vectors
!  p = p_re + i*p_im,  q = q_re + i*q_im
!  <p,q> = <p_re,q_re> + <p_im,q_im> + i*(<p_re,q_im> - <p_im,q_re>)

   use krylov_subspace
   implicit none
   include 'SIZE'
   include 'TOTAL'

   type(krylov_vector), intent(in) :: p_re, p_im, q_re, q_im
   complex(kind=kind(0.0d0)), intent(out) :: alpha

   real :: rr, ii, ri, ir

   call k_dot(rr, p_re, q_re)  ! <p_re, q_re>
   call k_dot(ii, p_im, q_im)  ! <p_im, q_im>
   call k_dot(ri, p_re, q_im)  ! <p_re, q_im>
   call k_dot(ir, p_im, q_re)  ! <p_im, q_re>

   alpha = cmplx(rr + ii, ri - ir, kind=kind(0.0d0))
end subroutine
```

---

## SPOD Implementation Sketch

```fortran
subroutine spod_compute(snaps, nsnap, dt, nfft, noverlap)
   use krylov_subspace
   use fourier_fftw
   implicit none
   include 'SIZE'
   include 'TOTAL'

   integer, intent(in) :: nsnap, nfft, noverlap
   real, intent(in) :: dt
   type(krylov_vector), intent(in) :: snaps(nsnap)

   integer :: nblk, nfreq, iblk, ifreq, i, j, blk_start
   real, allocatable :: window(:)
   real :: win_weight
   type(krylov_vector), allocatable :: blk_re(:), blk_im(:)  ! FFT results per block
   complex(kind=kind(0.0d0)), allocatable :: CSD(:,:)

!  Number of blocks
   nblk = (nsnap - nfft) / (nfft - noverlap) + 1
   nfreq = nfft / 2 + 1

!  Allocate arrays
   allocate(window(nfft))
   allocate(blk_re(nblk), blk_im(nblk))
   allocate(CSD(nblk, nblk))

!  Hamming window
   call hamming_window(nfft, window, win_weight)

!  For each frequency
   do ifreq = 1, nfreq

!     FFT each block at this frequency
      do iblk = 1, nblk
         blk_start = (iblk - 1) * (nfft - noverlap) + 1

!        Apply window and FFT (per spatial point)
         call fft_block_at_freq(snaps(blk_start:blk_start+nfft-1), &
                                window, win_weight, dt, ifreq, nfft, &
                                blk_re(iblk), blk_im(iblk))
      end do

!     Form CSD matrix (Hermitian: only compute upper triangle)
      do j = 1, nblk
         do i = 1, j
            call k_dot_complex(CSD(i,j), blk_re(i), blk_im(i), &
                                          blk_re(j), blk_im(j))
            CSD(j,i) = conjg(CSD(i,j))
         end do
      end do
      CSD = CSD / nblk

!     Eigensolve and reconstruct modes
      call spod_eigen_and_save(CSD, blk_re, blk_im, nblk, ifreq, nfreq, dt)

   end do

   deallocate(window, blk_re, blk_im, CSD)
end subroutine
```

---

## Memory Considerations

### Snapshot Storage
```
Memory = nsnap × (3 × lv + lp + ldimt × lv) × 8 bytes
       ≈ nsnap × 4 × lv × 8 bytes  (for velocity-only)

Example: lv = 1M points, nsnap = 500
Memory ≈ 500 × 4 × 1M × 8 = 16 GB
```

### For Large 3D Cases
- Load snapshots in chunks for SPOD (block-by-block)
- Or use streaming algorithm with incremental SVD
- Online mode (0.2) computes during DNS without storing all snapshots

---

## Helper Routines Needed

### FFT Per Block (for SPOD)
```fortran
subroutine fft_block_at_freq(block_snaps, window, win_weight, dt, ifreq, nfft, out_re, out_im)
!  FFT a block of snapshots and extract frequency component ifreq
!
!  For each spatial point (i.e., each entry in vx, vy, vz, t):
!    1. Apply window: x_windowed(j) = window(j) * block_snaps(j)%field(pt)
!    2. FFT the time series
!    3. Extract real/imag at frequency ifreq
!    4. Apply normalization: 2/(nfft * win_weight * dt) for one-sided

   use krylov_subspace
   use fourier_fftw
   implicit none
   include 'SIZE'
   include 'TOTAL'

   integer, intent(in) :: nfft, ifreq
   real, intent(in) :: dt, win_weight
   real, intent(in) :: window(nfft)
   type(krylov_vector), intent(in) :: block_snaps(nfft)
   type(krylov_vector), intent(out) :: out_re, out_im

   real :: time_series(nfft), fft_re(nfft), fft_im(nfft)
   real :: norm_factor
   integer :: ipt, j

!  Normalization: 2 for one-sided, divided by (nfft * window_energy * dt)
!  Factor of 2 accounts for discarding negative frequencies
   norm_factor = 2.0d0 / (real(nfft) * win_weight * dt)
   if (ifreq == 1 .or. ifreq == nfft/2+1) norm_factor = norm_factor / 2.0d0  ! DC and Nyquist

   call k_zero(out_re)
   call k_zero(out_im)

!  Process velocity vx
   do ipt = 1, nv
      do j = 1, nfft
         time_series(j) = window(j) * block_snaps(j)%vx(ipt)
      end do
      call fft_r2c_1d(time_series, fft_re, fft_im, nfft)
      out_re%vx(ipt) = fft_re(ifreq) * norm_factor
      out_im%vx(ipt) = fft_im(ifreq) * norm_factor
   end do

!  Similarly for vy, vz (if 3D), t (if thermal)
!  ... (same pattern for each field)

end subroutine
```

### Hamming Window
```fortran
subroutine hamming_window(n, window, win_weight)
   implicit none
   integer, intent(in) :: n
   real, intent(out) :: window(n), win_weight
   real, parameter :: PI = 3.14159265358979323846d0
   integer :: i

   do i = 1, n
      window(i) = 0.54d0 - 0.46d0 * cos(2.0d0 * PI * (i-1) / (n-1))
   end do

!  Window energy (for normalization)
   win_weight = sum(window**2) / n
end subroutine
```

## Implementation Phases

### Phase 1: Infrastructure
- [ ] Add modal parameters to NEKSTAB.inc
- [ ] Add modal_analysis dispatcher to main.f90 (Mode 6)
- [ ] Implement compute_mean, subtract_mean
- [ ] Implement eig_symmetric wrapper (DSYEV)
- [x] k_add2s2 added to krylov_subspace.f90
- [ ] Test k_dot consistency

### Phase 2: POD
- [ ] Implement correlation matrix formation using k_dot
- [ ] Implement pod_reconstruct_modes
- [ ] Output pod_energy.dat
- [ ] Validate on cylinder (St ≈ 0.167)

### Phase 3: DMD
- [ ] Implement SVD via eigendecomposition of Q^T Q
- [ ] Implement projected DMD operator
- [ ] Complex mode output (dRe, dIm)
- [ ] Output dmd_spectrum.dat (freq, growth, amplitude)

### Phase 4: SPOD
- [ ] Implement k_dot_complex
- [ ] Implement fft_block_at_freq using fourier_fftw
- [ ] Implement hamming_window
- [ ] Implement CSD formation
- [ ] Implement eig_hermitian wrapper (ZHEEV) - see below
- [ ] Complex mode output (sRe, sIm)
- [ ] Output spod_spectrum.dat

### ZHEEV Wrapper (for SPOD)
```fortran
subroutine eig_hermitian(A, eigvals, eigvecs, n)
!  Wrapper for LAPACK ZHEEV (complex Hermitian eigenvalue problem)
!  Returns eigenvalues in DESCENDING order

   implicit none
   integer, intent(in) :: n
   complex(kind=kind(0.0d0)), intent(in) :: A(n, n)
   real, intent(out) :: eigvals(n)
   complex(kind=kind(0.0d0)), intent(out) :: eigvecs(n, n)

   complex(kind=kind(0.0d0)), allocatable :: work(:), Acopy(:,:)
   real, allocatable :: rwork(:)
   complex(kind=kind(0.0d0)) :: tmp_vec(n)
   real :: tmp_val
   integer :: lwork, info, i, j

   allocate(Acopy(n,n), rwork(3*n-2))
   Acopy = A

!  Query optimal workspace
   allocate(work(1))
   call zheev('V', 'U', n, Acopy, n, eigvals, work, -1, rwork, info)
   lwork = int(real(work(1)))
   deallocate(work)
   allocate(work(lwork))

!  Compute eigenvalues/vectors (ascending order)
   call zheev('V', 'U', n, Acopy, n, eigvals, work, lwork, rwork, info)
   eigvecs = Acopy

!  Reverse to descending order
   do i = 1, n/2
      j = n - i + 1
      tmp_val = eigvals(i)
      eigvals(i) = eigvals(j)
      eigvals(j) = tmp_val
      tmp_vec = eigvecs(:, i)
      eigvecs(:, i) = eigvecs(:, j)
      eigvecs(:, j) = tmp_vec
   end do

   deallocate(work, rwork, Acopy)
end subroutine
```

### Phase 5: Polish
- [ ] Add checkpointing
- [ ] Validate against PySPOD
- [ ] Document in README

---

## Output Files

| Method | Modes | Spectrum |
|--------|-------|----------|
| POD | `pod_<SESSION>0.f*` | `pod_energy.dat` |
| DMD | `dRe_<SESSION>0.f*`, `dIm_<SESSION>0.f*` | `dmd_spectrum.dat` |
| SPOD | `sRe_<SESSION>0.f*`, `sIm_<SESSION>0.f*` | `spod_spectrum.dat` |

---

## Main Dispatcher Integration

```fortran
! In main.f90, add case 6:

case (6)  ! Modal analysis
   if (nid == 0) write(6,*) 'Modal analysis (POD/DMD/SPOD)...'
   call modal_analysis
   call nek_end
```

---

## Summary

- Uses existing nekStab `krylov_vector` type and operations
- Uses `k_dot` / `glsc3` for MPI-safe inner products
- Uses `k_add2s2` for efficient mode reconstruction (p = p + c*q)
- Uses `load_fld` / `outpost2` for I/O
- Complex modes stored as separate real/imag files
- Matches nekStab coding conventions
