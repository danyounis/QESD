! Library of Mathematical Functions (math.f08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 3/29/2020
! Revised: 11/29/2022
!
! References:
! [1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
! Numerical Recipes in Fortran 90 (Cambridge University Press, Cambridge, 2001).

module math

use prec
use omp_lib

implicit none

real(8), parameter, public :: pi = 4.0d0*datan(1.0d0)
complex, parameter, public :: i = (0,1)

interface linspace
    module procedure linspace_n
    module procedure linspace_d
end interface linspace

interface rotM90
    module procedure rotM90R
    module procedure rotM90C
end interface

interface fliplr
    module procedure fliplr_ra
    module procedure fliplr_ca
    module procedure fliplr_rm
    module procedure fliplr_cm
end interface fliplr

interface diag
    module procedure diag_r
    module procedure diag_c
end interface diag

interface inverse
    module procedure inverse_r
    module procedure inverse_c
end interface inverse

interface Kproduct
    module procedure Kproduct_r
    module procedure Kproduct_c
end interface Kproduct

interface tridiag_matmul_dmat
    module procedure tridiag_matmul_dmat_sp_rrr
    module procedure tridiag_matmul_dmat_sp_crc
    module procedure tridiag_matmul_dmat_sp_ccr
    module procedure tridiag_matmul_dmat_sp_ccc
end interface tridiag_matmul_dmat

interface tridiag_matmul_cvec
    module procedure tridiag_matmul_cvec_sp_rrr
    module procedure tridiag_matmul_cvec_sp_crc
    module procedure tridiag_matmul_cvec_sp_ccr
    module procedure tridiag_matmul_cvec_sp_ccc
end interface tridiag_matmul_cvec

interface tridiag_fbwd_subs
    module procedure tridiag_fbwd_subs_sp_rrr
    module procedure tridiag_fbwd_subs_sp_crc
    module procedure tridiag_fbwd_subs_sp_rcc
    module procedure tridiag_fbwd_subs_sp_ccc
end interface tridiag_fbwd_subs

interface fft_shift
    module procedure fft_shift_1d
    module procedure fft_shift_2d
end interface fft_shift

interface fft
    module procedure fft_1d
    module procedure fft_2d
end interface fft

interface ifft
    module procedure ifft_1d
    module procedure ifft_2d
end interface ifft

interface real_stagger_complex
    module procedure real_stagger_complex_1d
    module procedure real_stagger_complex_2d
end interface real_stagger_complex

interface complex_stagger_real
    module procedure complex_stagger_real_1d
    module procedure complex_stagger_real_2d
end interface complex_stagger_real

interface zpad_signal
    module procedure zero_pad_signal_r1d
    module procedure zero_pad_signal_c1d
end interface zpad_signal

interface deriv
    module procedure deriv_r1d
    module procedure deriv_c1d
    module procedure deriv_r2d
    module procedure deriv_c2d
end interface deriv

interface grad
    module procedure grad_r1d
    module procedure grad_c1d
    module procedure grad_r2d
    module procedure grad_c2d
end interface grad

interface trapz
    module procedure trapz_r1d
    module procedure trapz_c1d
    module procedure trapz_r2d_part
    module procedure trapz_r2d_full
    module procedure trapz_c2d_part
    module procedure trapz_c2d_full
end interface trapz

interface phase_unwrap
    module procedure phase_unwrap_1d
    module procedure phase_unwrap_2d
end interface phase_unwrap

interface bcuint
    module procedure bcuint_r
    module procedure bcuint_c
end interface bcuint

contains

function linspace_n(a,b,n) result(vec)
    implicit none
    real(num), intent(in)   :: a, b
    integer,   intent(in)   :: n
    real(num), dimension(n) :: vec

    real(num) :: dx
    integer :: j

    dx = (b-a)/(n-1)
    ForAll(j=1:n) vec(j) = a + (j-1)*dx

    return
end function linspace_n

function linspace_d(a,b,dx) result(vec)
    implicit none
    real(num), intent(in) :: a, b
    real(num), intent(in) :: dx
    real(num), dimension(nint((b-a)/dx)+1) :: vec
    integer :: n, j

    n = nint((b-a)/dx)+1
    ForAll(j=1:n) vec(j) = a + (j-1)*dx

    return
end function linspace_d

function zeros(n) result(vec)
    implicit none
    integer, intent(in) :: n
    real(num), dimension(n) :: vec
    vec = 0.0_num
    return
end function zeros

function ones(n) result(vec)
    implicit none
    integer, intent(in) :: n
    real(num), dimension(n) :: vec
    vec = 1.0_num
    return
end function ones

function Kdelta(n,m) result(r)
    implicit none
    integer, intent(in) :: n, m
    integer :: r

    select case (n == m)
    case (.true.)
        r = 1
    case (.false.)
        r = 0
    end select

    return
end function Kdelta

function rotM90R(M) result(V)
    implicit none
    real(num), intent(in), dimension(:,:) :: M
    real(num), dimension(size(M,1),size(M,2)) :: V
    integer :: j, k

    j = 0
    do k=ubound(M,2),lbound(M,2),-1
        j = j + 1
        V(j,:) = M(:,k)
    end do

    return
end function rotM90R

function rotM90C(M) result(V)
    implicit none
    complex(num), intent(in), dimension(:,:) :: M
    complex(num), dimension(size(M,1),size(M,2)) :: V
    integer :: j, k

    j = 0
    do k=ubound(M,2),lbound(M,2),-1
        j = j + 1
        V(j,:) = M(:,k)
    end do

    return
end function rotM90C

function fliplr_ra(arr) result(f_arr)
    implicit none
    real(num), intent(in), dimension(:) :: arr
    real(num), dimension(size(arr)) :: f_arr
    integer :: n, j
    n = size(arr,1)
    do j=n,1,-1
        f_arr(n-j+1) = arr(j)
    end do
    return
end function fliplr_ra

function fliplr_ca(arr) result(f_arr)
    implicit none
    complex(num), intent(in), dimension(:) :: arr
    complex(num), dimension(size(arr)) :: f_arr
    f_arr = fliplr(real(arr)) + i*fliplr(aimag(arr))
    return
end function fliplr_ca

function fliplr_rm(M,ax) result(f_M)
    implicit none
    real(num), intent(in) :: M(:,:)
    integer, intent(in) :: ax

    real(num), dimension(size(M,1),size(M,2)) :: f_M
    integer :: nn(2), j1, j2

    nn = (/size(M,1),size(M,2)/)

    select case (ax)
    case (0)
        do j1=1,nn(1)
            f_M(j1,:) = fliplr(M(j1,:))
        end do
        do j2=1,nn(2)
            f_M(:,j2) = fliplr(f_M(:,j2))
        end do
    case (1)
        do j2=1,nn(2)
            f_M(:,j2) = fliplr(M(:,j2))
        end do
    case (2)
        do j1=1,nn(1)
            f_M(j1,:) = fliplr(M(j1,:))
        end do
    end select

    return
end function fliplr_rm

function fliplr_cm(M,ax) result(f_M)
    implicit none
    complex(num), intent(in) :: M(:,:)
    integer, intent(in) :: ax
    complex(num), dimension(size(M,1),size(M,2)) :: f_M
    f_M = fliplr(real(M),ax) + i*fliplr(aimag(M),ax)
    return
end function fliplr_cm

function identity(n) result(A)
    implicit none
    integer,   intent(in) :: n
    real(num), dimension(n,n) :: A
    integer :: j
    A = 0.0_num
    ForAll(j=1:n) A(j,j) = 1.0_num
    return
end function identity

function diag_r(arr) result(M)
    implicit none
    real(num), intent(in), dimension(:) :: arr
    real(num), dimension(size(arr),size(arr)) :: M
    integer :: j
    M = 0.0_num
    ForAll(j=1:size(arr,1)) M(j,j) = arr(j)
    return
end function diag_r

function diag_c(arr) result(M)
    implicit none
    complex(num), intent(in), dimension(:) :: arr
    complex(num), dimension(size(arr),size(arr)) :: M
    integer :: j
    M = 0.0_num
    ForAll(j=1:size(arr,1)) M(j,j) = arr(j)
    return
end function diag_c

function inverse_r(M) result(C)
    implicit none
    real(num), intent(in) :: M(:,:)

    real(num), dimension(size(M,1)) :: b, d, x
    real(num), dimension(size(M,1),size(M,1)) :: C, L, U, a

    real(num) :: coeff
    integer :: n, o, j, k

    n = size(M,1); a = M;
    L = 0.0_num; U = 0.0_num; b = 0.0_num;

    ! forward elimination
    do k=1,n-1
        do o=k+1,n
            coeff = a(o,k)/a(k,k)
            L(o,k) = coeff
            do j=k+1,n
                a(o,j) = a(o,j) - coeff*a(k,j)
            end do
        end do
    end do

    ! construct L and U matrices
    ForAll(o=1:n) L(o,o) = 1.0_num
    do j=1,n
        do o=1,j
            U(o,j) = a(o,j)
        end do
    end do

    ! compute columns of the inverse matrix C
    do k=1,n
        b(k) = 1.0_num
        d(1) = b(1)
        ! solve Ld=b using forward substitution
        do o=2,n
            d(o) = b(o)
            do j=1,o-1
                d(o) = d(o) - L(o,j)*d(j)
            end do
        end do
        ! solve Ux=d using backward substitution
        x(n) = d(n)/U(n,n)
        do o=n-1,1,-1
            x(o) = d(o)
            do j=n,o+1,-1
                x(o) = x(o) - U(o,j)*x(j)
            end do
            x(o) = x(o)/U(o,o)
        end do
        ! fill the solutions x(n) into column k of C
        do o=1,n
            C(o,k) = x(o)
        end do
        b(k) = 0.0_num
    end do

    return
end function inverse_r

function inverse_c(M) result(C)
    implicit none
    complex(num), intent(in) :: M(:,:)

    complex(num), dimension(size(M,1)) :: b, d, x
    complex(num), dimension(size(M,1),size(M,1)) :: C, L, U, a

    complex(num) :: coeff
    integer :: n, o, j, k

    n = size(M,1); a = M;
    L = (0.0_num,0.0_num); U = (0.0_num,0.0_num); b = (0.0_num,0.0_num);

    ! forward elimination
    do k=1,n-1
        do o=k+1,n
            coeff = a(o,k)/a(k,k)
            L(o,k) = coeff
            do j=k+1,n
                a(o,j) = a(o,j) - coeff*a(k,j)
            end do
        end do
    end do

    ! construct L and U matrices
    ForAll(o=1:n) L(o,o) = (1.0_num,1.0_num)
    do j=1,n
        do o=1,j
            U(o,j) = a(o,j)
        end do
    end do

    ! compute columns of the inverse matrix C
    do k=1,n
        b(k) = (1.0_num,1.0_num)
        d(1) = b(1)
        ! solve Ld=b using forward substitution
        do o=2,n
            d(o) = b(o)
            do j=1,o-1
                d(o) = d(o) - L(o,j)*d(j)
            end do
        end do
        ! solve Ux=d using backward substitution
        x(n) = d(n)/U(n,n)
        do o=n-1,1,-1
            x(o) = d(o)
            do j=n,o+1,-1
                x(o) = x(o) - U(o,j)*x(j)
            end do
            x(o) = x(o)/U(o,o)
        end do
        ! fill the solutions x(n) into column k of C
        do o=1,n
            C(o,k) = x(o)
        end do
        b(k) = (0.0_num,0.0_num)
    end do

    return
end function inverse_c

function Kproduct_r(A,B) result(AB)
    implicit none
    real(num), intent(in), dimension(:,:) :: A, B
    complex(num) :: AB(size(A,1)*size(B,1),size(A,2)*size(B,2))

    integer :: R,RA,RB,C,CA,CB
    integer :: J,K

    RA = size(A,1)
    CA = size(A,2)
    RB = size(B,1)
    CB = size(B,2)

    R = 0
    do J=1,RA
        C = 0
        do K=1,CA
            AB(R+1:R+RB,C+1:C+CB) = A(J,K)*B
            C = C + CB
        end do
        R = R + RB
    end do
    return
end function Kproduct_r

function Kproduct_c(A,B) result(AB)
    implicit none
    complex(num), intent(in), dimension(:,:) :: A, B
    complex(num) :: AB(size(A,1)*size(B,1),size(A,2)*size(B,2))

    integer :: R,RA,RB,C,CA,CB
    integer :: J,K

    RA = size(A,1)
    CA = size(A,2)
    RB = size(B,1)
    CB = size(B,2)

    R = 0
    do J=1,RA
        C = 0
        do K=1,CA
            AB(R+1:R+RB,C+1:C+CB) = A(J,K)*B
            C = C + CB
        end do
        R = R + RB
    end do
    return
end function Kproduct_c

subroutine tridiag_matmul_dmat_sp_rrr(y,A,V)
    implicit none
    real(num), intent(out) :: y(:)
    real(num), intent(in)  :: A(:)
    real(num), intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1)
    y(2) = A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=3, 3*n-6, 3
        y(j) = A(j)*V(j/3)
        y(j+1) = A(j+1)*V(j/3 + 1)
        y(j+2) = A(j+2)*V(j/3 + 2)
    end do
    !$OMP END PARALLEL DO
    y(3*n-3) = A(3*n-3)*V(n-1)
    y(3*n-2) = A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_dmat_sp_rrr

subroutine tridiag_matmul_dmat_sp_crc(y,A,V)
    implicit none
    complex(num), intent(out) :: y(:)
    real(num),    intent(in)  :: A(:)
    complex(num), intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1)
    y(2) = A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=3, 3*n-6, 3
        y(j) = A(j)*V(j/3)
        y(j+1) = A(j+1)*V(j/3 + 1)
        y(j+2) = A(j+2)*V(j/3 + 2)
    end do
    !$OMP END PARALLEL DO
    y(3*n-3) = A(3*n-3)*V(n-1)
    y(3*n-2) = A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_dmat_sp_crc

subroutine tridiag_matmul_dmat_sp_ccr(y,A,V)
    implicit none
    complex(num), intent(out) :: y(:)
    complex(num), intent(in)  :: A(:)
    real(num),    intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1)
    y(2) = A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=3, 3*n-6, 3
        y(j) = A(j)*V(j/3)
        y(j+1) = A(j+1)*V(j/3 + 1)
        y(j+2) = A(j+2)*V(j/3 + 2)
    end do
    !$OMP END PARALLEL DO
    y(3*n-3) = A(3*n-3)*V(n-1)
    y(3*n-2) = A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_dmat_sp_ccr

subroutine tridiag_matmul_dmat_sp_ccc(y,A,V)
    implicit none
    complex(num), intent(out) :: y(:)
    complex(num), intent(in)  :: A(:)
    complex(num), intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1)
    y(2) = A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=3, 3*n-6, 3
        y(j) = A(j)*V(j/3)
        y(j+1) = A(j+1)*V(j/3 + 1)
        y(j+2) = A(j+2)*V(j/3 + 2)
    end do
    !$OMP END PARALLEL DO
    y(3*n-3) = A(3*n-3)*V(n-1)
    y(3*n-2) = A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_dmat_sp_ccc

subroutine tridiag_matmul_cvec_sp_rrr(y,A,V)
    implicit none
    real(num), intent(out) :: y(:)
    real(num), intent(in)  :: A(:)
    real(num), intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1) + A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=2,n-1
        y(j) = A(3*j-3)*V(j-1) + A(3*j-2)*V(j) + A(3*j-1)*V(j+1)
    end do
    !$OMP END PARALLEL DO
    y(n) = A(3*n-3)*V(n-1) + A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_cvec_sp_rrr

subroutine tridiag_matmul_cvec_sp_crc(y,A,V)
    implicit none
    complex(num), intent(out) :: y(:)
    real(num),    intent(in)  :: A(:)
    complex(num), intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1) + A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=2,n-1
        y(j) = A(3*j-3)*V(j-1) + A(3*j-2)*V(j) + A(3*j-1)*V(j+1)
    end do
    !$OMP END PARALLEL DO
    y(n) = A(3*n-3)*V(n-1) + A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_cvec_sp_crc

subroutine tridiag_matmul_cvec_sp_ccr(y,A,V)
    implicit none
    complex(num), intent(out) :: y(:)
    complex(num), intent(in)  :: A(:)
    real(num),    intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1) + A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=2,n-1
        y(j) = A(3*j-3)*V(j-1) + A(3*j-2)*V(j) + A(3*j-1)*V(j+1)
    end do
    !$OMP END PARALLEL DO
    y(n) = A(3*n-3)*V(n-1) + A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_cvec_sp_ccr

subroutine tridiag_matmul_cvec_sp_ccc(y,A,V)
    implicit none
    complex(num), intent(out) :: y(:)
    complex(num), intent(in)  :: A(:)
    complex(num), intent(in)  :: V(:)

    integer :: n, j
    n = size(V)

    y(1) = A(1)*V(1) + A(2)*V(2)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=2,n-1
        y(j) = A(3*j-3)*V(j-1) + A(3*j-2)*V(j) + A(3*j-1)*V(j+1)
    end do
    !$OMP END PARALLEL DO
    y(n) = A(3*n-3)*V(n-1) + A(3*n-2)*V(n)

    return
end subroutine tridiag_matmul_cvec_sp_ccc

subroutine tridiag_fbwd_subs_sp_rrr(A,x,b)
    implicit none
    real(num), intent(in)    :: A(:)
    real(num), intent(out)   :: x(:)
    real(num), intent(inout) :: b(:)

    integer :: n, j

    ! temp vars
    real(num) :: Ap(size(b))
    n = size(b)

    ! fwd. subs. Ap.x=b
    Ap(1) = A(1)
    do j=1,n-1
        Ap(j+1) = A(4+3*(j-1)) - (A(3+3*(j-1))/Ap(j))*A(2+3*(j-1))
        b(j+1) = b(j+1) - (A(3+3*(j-1))/Ap(j))*b(j)
    end do

    ! bwd. subs. x=Ap\b
    x(n) = b(n)/Ap(n)
    do j=n-1, 1, -1
        x(j) = (b(j) - A(3*j-1)*x(j+1))/Ap(j)
    end do

    return
end subroutine tridiag_fbwd_subs_sp_rrr

subroutine tridiag_fbwd_subs_sp_crc(A,x,b)
    implicit none
    complex(num), intent(in)    :: A(:)
    real(num),    intent(out)   :: x(:)
    complex(num), intent(inout) :: b(:)

    integer :: n, j

    ! temp vars
    complex(num) :: Ap(size(b))
    n = size(b)

    ! fwd. subs. Ap.x=b
    Ap(1) = A(1)
    do j=1,n-1
        Ap(j+1) = A(4+3*(j-1)) - (A(3+3*(j-1))/Ap(j))*A(2+3*(j-1))
        b(j+1) = b(j+1) - (A(3+3*(j-1))/Ap(j))*b(j)
    end do

    ! bwd. subs. x=Ap\b
    x(n) = b(n)/Ap(n)
    do j=n-1, 1, -1
        x(j) = (b(j) - A(3*j-1)*x(j+1))/Ap(j)
    end do

    return
end subroutine tridiag_fbwd_subs_sp_crc

subroutine tridiag_fbwd_subs_sp_rcc(A,x,b)
    implicit none
    real(num),    intent(in)    :: A(:)
    complex(num), intent(out)   :: x(:)
    complex(num), intent(inout) :: b(:)

    integer :: n, j

    ! temp vars
    real(num) :: Ap(size(b))
    n = size(b)

    ! fwd. subs. Ap.x=b
    Ap(1) = A(1)
    do j=1,n-1
        Ap(j+1) = A(4+3*(j-1)) - (A(3+3*(j-1))/Ap(j))*A(2+3*(j-1))
        b(j+1) = b(j+1) - (A(3+3*(j-1))/Ap(j))*b(j)
    end do

    ! bwd. subs. x=Ap\b
    x(n) = b(n)/Ap(n)
    do j=n-1, 1, -1
        x(j) = (b(j) - A(3*j-1)*x(j+1))/Ap(j)
    end do

    return
end subroutine tridiag_fbwd_subs_sp_rcc

subroutine tridiag_fbwd_subs_sp_ccc(A,x,b)
    implicit none
    complex(num), intent(in)    :: A(:)
    complex(num), intent(out)   :: x(:)
    complex(num), intent(inout) :: b(:)

    integer :: n, j

    ! temp vars
    complex(num) :: Ap(size(b))
    n = size(b)

    ! fwd. subs. Ap.x=b
    Ap(1) = A(1)
    do j=1,n-1
        Ap(j+1) = A(4+3*(j-1)) - (A(3+3*(j-1))/Ap(j))*A(2+3*(j-1))
        b(j+1) = b(j+1) - (A(3+3*(j-1))/Ap(j))*b(j)
    end do

    ! bwd. subs. x=Ap\b
    x(n) = b(n)/Ap(n)
    do j=n-1, 1, -1
        x(j) = (b(j) - A(3*j-1)*x(j+1))/Ap(j)
    end do

    return
end subroutine tridiag_fbwd_subs_sp_ccc

subroutine svdcmp(a,m,n,w,v)
    implicit none
    integer, intent(in) :: m, n
    real(num), intent(inout) :: a(m,n)
    real(num), intent(out) :: v(n,n), w(n)
    integer, parameter :: NMAX = 500 ! max anticipated value of n

    real(num) :: anorm, c, f, g, h, s, scale, x, y, z, rv1(NMAX)
    integer :: i, its, j, jj, k, l, nm

    g = 0.0
    scale = 0.0
    anorm = 0.0

    do i=1,n
        l = i+1
        rv1(i) = scale*g
        g = 0.0
        s = 0.0
        scale = 0.0
        if (i <= m) then
            do k=i,m
                scale = scale + abs(a(k,i))
            end do
            if (scale /= 0.0) then
                do k=i,m
                    a(k,i) = a(k,i)/scale
                    s = s + a(k,i)*a(k,i)
                end do
                f = a(i,i)
                g = -sign(sqrt(s),f)
                h = f*g - s
                a(i,i) = f - g
                do j=l,n
                    s = 0.0
                    do k=i,m
                        s = s + a(k,i)*a(k,j)
                    end do
                    f = s/h
                    do k=i,m
                        a(k,j) = a(k,j) + f*a(k,i)
                    end do
                end do
                do k=i,m
                    a(k,i) = scale*a(k,i)
                end do
            end if
        end if
        w(i) = scale*g
        g = 0.0
        s = 0.0
        scale = 0.0
        if ((i <= m).and.(i /= n)) then
            do k=l,n
                scale = scale + abs(a(i,k))
            end do
            if (scale /= 0.0) then
                do k=l,n
                    a(i,k) = a(i,k)/scale
                    s = s + a(i,k)*a(i,k)
                end do
                f = a(i,l)
                g = -sign(sqrt(s),f)
                h = f*g - s
                a(i,l) = f - g
                do k=l,n
                    rv1(k) = a(i,k)/h
                end do
                do j=l,m
                    s = 0.0
                    do k=l,n
                        s = s + a(j,k)*a(i,k)
                    end do
                    do k=l,n
                        a(j,k) = a(j,k) + s*rv1(k)
                    end do
                end do
                do k=l,n
                    a(i,k) = scale*a(i,k)
                end do
            end if
        end if
        anorm = max( anorm, (abs(w(i))+abs(rv1(i))) )
    end do
    do i=n,1,-1
        if (i < n) then
            if (g /= 0.0) then
                do j=l,n
                    v(j,i) = (a(i,j)/a(i,l))/g
                end do
                do j=l,n
                    s = 0.0
                    do k=l,n
                        s = s + a(i,k)*v(k,j)
                    end do
                    do k=l,n
                        v(k,j) = v(k,j) + s*v(k,i)
                    end do
                end do
            end if
            do j=l,n
                v(i,j) = 0.0
                v(j,i) = 0.0
            end do
        end if
        v(i,i) = 1.0
        g = rv1(i)
        l = i
    end do
    do i=min(m,n),1,-1
        l = i+1
        g = w(i)
        do j=l,n
            a(i,j) = 0.0
        end do
        if (g /= 0.0) then
            g = 1.0/g
            do j=l,n
                s = 0.0
                do k=l,m
                    s = s + a(k,i)*a(k,j)
                end do
                f = (s/a(i,i))*g
                do k=i,m
                    a(k,j) = a(k,j) + f*a(k,i)
                end do
            end do
            do j=i,m
                a(j,i) = a(j,i)*g
            end do
        else
            do j=i,m
                a(j,i) = 0.0
            end do
        end if
        a(i,i) = a(i,i) + 1.0
    end do
    do k=n,1,-1
        do its=1,30
            do l=k,1,-1
                nm = l-1
                if ((abs(rv1(l))+anorm) == anorm) goto 2
                if ((abs(w(nm))+anorm) == anorm) goto 1
            end do
1           c = 0.0
            s = 1.0
            do i=l,k
                f = s*rv1(i)
                rv1(i) = c*rv1(i)
                if ((abs(f)+anorm) == anorm) goto 2
                g = w(i)
                h = pythag(f,g)
                w(i) = h
                h = 1.0/h
                c = g*h
                s = -f*h
                do j=1,m
                    y = a(j,nm)
                    z = a(j,i)
                    a(j,nm) = y*c + z*s
                    a(j,i) = -y*s + z*c
                end do
            end do
2           z = w(k)
            if (l == k) then
                if (z < 0.0) then
                    w(k) = -z
                    do j=1,n
                        v(j,k) = -v(j,k)
                    end do
                end if
                goto 3
            end if
            ! if (its == 30) pause 'no convergence in svdcmp'
            x = w(l)
            nm = k-1
            y = w(nm)
            g = rv1(nm)
            h = rv1(k)
            f = ((y-z)*(y+z) + (g-h)*(g+h))/(2.0*h*y)
            g = pythag(f,1.0_num)
            f = ((x-z)*(x+z)+ h*((y/(f+sign(g,f)))-h))/x
            c = 1.0
            s = 1.0
            do j=l,nm
                i = j+1
                g = rv1(i)
                y = w(i)
                h = s*g
                g = c*g
                z = pythag(f,h)
                rv1(j) = z
                c = f/z
                s = h/z
                f = (x*c) + (g*s)
                g = -(x*s) + (g*c)
                h = y*s
                y = y*c
                do jj=1,n
                    x = v(jj,j)
                    z = v(jj,i)
                    v(jj,j) = (x*c) + (z*s)
                    v(jj,i) = -(x*s) + (z*c)
                end do
                z = pythag(f,h)
                w(j) = z
                if (z /= 0.0) then
                    z = 1.0/z
                    c = f*z
                    s = h*z
                end if
                f = (c*g) + (s*y)
                x = -(s*g) + (c*y)
                do jj=1,m
                    y = a(jj,j)
                    z = a(jj,i)
                    a(jj,j) = (y*c) + (z*s)
                    a(jj,i) = -(y*s) + (z*c)
                end do
            end do
            rv1(l) = 0.0
            rv1(k) = f
            w(k) = x
        end do
3       continue
    end do
    return
end subroutine svdcmp

real(num) function pythag(a,b)
    implicit none
    real(num), intent(in) :: a, b
    real(num) :: absa, absb

    absa = abs(a)
    absb = abs(b)
    if (absa > absb) then
        pythag = absa*sqrt(1.0 + (absb/absa)**2)
    else
        if (absb == 0.0) then
            pythag = 0.0
        else
            pythag = absb*sqrt(1.0 + (absa/absb)**2)
        end if
    end if

    return
end function pythag

function fft_freq(n,delta,shift) result(freq)
    implicit none
    integer, intent(in) :: n
    real(num), intent(in) :: delta
    logical, intent(in) :: shift
    real(num), dimension(n) :: freq
    integer :: j

    freq = 0.0_num
    do j=1,n
        if ((j >= 1).and.(j <= n/2)) then
            freq(j) = real(j-1,num)
        else if ((j > n/2).and.(j <= n)) then
            freq(j) = real(j-1-n,num)
        end if
    end do

    freq = freq/real(n*delta,num)

    if (shift) freq = cshift(freq,n/2)
    return
end function fft_freq

subroutine fft_shift_1d(func)
    implicit none
    complex(num), intent(inout) :: func(:)
    integer :: n
    n = size(func)
    func = cshift(func,n/2)
    return
end subroutine fft_shift_1d

subroutine fft_shift_2d(func)
    implicit none
    complex(num), intent(inout) :: func(:,:)
    integer :: nn(2), j
    nn = (/size(func,1),size(func,2)/)
    do j=1,nn(1)
        call fft_shift(func(j,:))
    end do
    do j=1,nn(2)
        call fft_shift(func(:,j))
    end do
    return
end subroutine fft_shift_2d

subroutine four1(data,nn,isign)
    implicit none
    integer,   intent(in) :: nn, isign
    real(num), intent(inout), dimension(2*nn) :: data

    integer   :: i, istep, j, m, mmax, n
    real(num) :: tempi, tempr
    real(num) :: theta, wi, wpi, wpr, wr, wtemp

    n = 2*nn
    j = 1
    do i=1,n,2
        if (j > i) then
            tempr = data(j)
            tempi = data(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
        end if
        m = n/2
1       if ((m >= 2).and.(j > m)) then
            j = j - m
            m = m/2
            goto 1
        end if
        j = j + m
    end do

    mmax = 2

2   if (n > mmax) then
        istep = 2*mmax
        theta = 2.0*pi/real(isign*mmax,num)
        wpr = -2.0*sin(0.5*theta)**2
        wpi = sin(theta)
        wr = 1.0
        wi = 0.0
        do m = 1,mmax,2
            do i = m,n,istep
                j = i + mmax
                tempr = wr*data(j) - wi*data(j+1)
                tempi = wr*data(j+1) + wi*data(j)
                data(j) = data(i) - tempr
                data(j+1) = data(i+1) - tempi
                data(i) = data(i) + tempr
                data(i+1) = data(i+1) + tempi
            end do
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
        end do
        mmax = istep
        goto 2
    end if
    return
end subroutine four1

subroutine four2(data,nn,isign)
    implicit none
    integer,   intent(in)  :: nn(2), isign
    real(num), intent(inout), dimension(2*nn(1)*nn(2)) :: data

    integer, parameter :: ndim = 2

    integer :: i1, i2, i2rev, i3, i3rev, ibit, idim, &
    ifp1, ifp2, ip1, ip2, ip3, k1, k2, n, nprev, nrem, ntot

    real(num) :: tempi, tempr
    real(num) :: theta, wi, wpi, wpr, wr, wtemp

    ntot = 1
    do idim=1,ndim
        ntot = ntot*nn(idim)
    end do
    nprev = 1
    do idim=1,ndim
        n = nn(idim)
        nrem = ntot/(n*nprev)
        ip1 = 2*nprev
        ip2 = ip1*n
        ip3 = ip2*nrem
        i2rev = 1
        do i2=1,ip2,ip1
            if (i2 < i2rev) then
                do i1=i2,i2+ip1-2,2
                    do i3=i1,ip3,ip2
                        i3rev = i2rev + i3 - i2
                        tempr = data(i3)
                        tempi = data(i3+1)
                        data(i3) = data(i3rev)
                        data(i3+1) = data(i3rev+1)
                        data(i3rev) = tempr
                        data(i3rev+1) = tempi
                    end do
                end do
            end if
            ibit = ip2/2
1           if ((ibit >= ip1).and.(i2rev > ibit)) then
                i2rev = i2rev - ibit
                ibit = ibit/2
                goto 1
            end if
            i2rev = i2rev + ibit
        end do
        ifp1 = ip1
2       if (ifp1 < ip2) then
            ifp2 = 2*ifp1
            theta = isign*2.0*pi/(ifp2/ip1)
            wpr = -2.0*sin(0.5*theta)**2
            wpi = sin(theta)
            wr = 1.0
            wi = 0.0
            do i3=1,ifp1,ip1
                do i1=i3,i3+ip1-2,2
                    do i2=i1,ip3,ifp2
                        k1 = i2
                        k2 = k1+ifp1
                        tempr = wr*data(k2) - wi*data(k2+1)
                        tempi = wr*data(k2+1) + wi*data(k2)
                        data(k2) = data(k1) - tempr
                        data(k2+1) = data(k1+1) - tempi
                        data(k1) = data(k1) + tempr
                        data(k1+1) = data(k1+1) + tempi
                    end do
                end do
                wtemp = wr
                wr = wr*wpr - wi*wpi + wr
                wi = wi*wpr + wtemp*wpi + wi
            end do
            ifp1 = ifp2
            goto 2
        end if
        nprev = n*nprev
    end do
    return
end subroutine four2

function fft_1d(func) result(ft_func)
    implicit none
    complex(num), intent(in) :: func(:)

    complex(num), dimension(size(func)) :: ft_func
    real(num),    dimension(2*size(func)) :: ft_temp

    integer :: n
    n = size(func)

    ft_temp = real_stagger_complex(func,n)
    call four1(ft_temp,n,1)
    ft_func = complex_stagger_real(ft_temp,n)
    return
end function fft_1d

function fft_2d(func) result(ft_func)
    implicit none
    complex(num), intent(in) :: func(:,:)

    complex(num), dimension(size(func,1),size(func,2)) :: ft_func
    real(num),    dimension(2*size(func,1)*size(func,2)) :: ft_temp

    integer :: nn(2)
    nn = (/size(func,1),size(func,2)/)

    ft_temp = real_stagger_complex(func,nn)
    call four2(ft_temp,nn,1)
    ft_func = complex_stagger_real(ft_temp,nn)
    return
end function fft_2d

function ifft_1d(func) result(ift_func)
    implicit none
    complex(num), intent(in) :: func(:)

    complex(num), dimension(size(func)) :: ift_func
    real(num),    dimension(2*size(func)) :: ift_temp

    integer :: n
    n = size(func)

    ift_temp = real_stagger_complex(func,n)
    call four1(ift_temp,n,-1)
    ift_func = complex_stagger_real(ift_temp,n)/real(n,num)
    return
end function ifft_1d

function ifft_2d(func) result(ift_func)
    implicit none
    complex(num), intent(in) :: func(:,:)

    complex(num), dimension(size(func,1),size(func,2)) :: ift_func
    real(num),    dimension(2*size(func,1)*size(func,2)) :: ift_temp

    integer :: nn(2)
    nn = (/size(func,1),size(func,2)/)

    ift_temp = real_stagger_complex(func,nn)
    call four2(ift_temp,nn,-1)
    ift_func = complex_stagger_real(ift_temp,nn)/real(nn(1)*nn(2),num)
    return
end function ifft_2d

function real_stagger_complex_1d(g,n) result(f)
    implicit none
    integer,      intent(in)               :: n
    complex(num), intent(in), dimension(n) :: g

    real(num), dimension(2*n) :: f
    integer :: j

    do j=1,n
        f(2*j-1) = real(g(j))
        f(2*j) = aimag(g(j))
    end do
    return
end function real_stagger_complex_1d

function real_stagger_complex_2d(g,nn) result(f)
    implicit none
    integer,      intent(in), dimension(2) :: nn
    complex(num), intent(in), dimension(nn(1),nn(2)) :: g

    real(num), dimension(2*nn(1)*nn(2)) :: f
    integer :: j1, j2, k

    f = 0.0_num; k = 1;
    do j2=1,nn(2)
        do j1=1,nn(1)
            f(k) = real(g(j1,j2))
            f(k+1) = aimag(g(j1,j2))
            k = k + 2
        end do
    end do
    return
end function real_stagger_complex_2d

function complex_stagger_real_1d(f,n) result(g)
    implicit none
    integer,      intent(in)                 :: n
    real(num),    intent(in), dimension(2*n) :: f

    complex(num), dimension(n) :: g
    integer :: j

    do j=1,n
        g(j) = f(2*j-1) + i*f(2*j)
    end do
    return
end function complex_stagger_real_1d

function complex_stagger_real_2d(f,nn) result(g)
    implicit none
    integer,   intent(in), dimension(2) :: nn
    real(num), intent(in), dimension(2*nn(1)*nn(2)) :: f

    complex(num), dimension(nn(1),nn(2)) :: g
    integer :: j1, j2, k

    g = 0.0_num; k = 1;
    do j2=1,nn(2)
        do j1=1,nn(1)
            g(j1,j2) = f(k) + i*f(k+1)
            k = k + 2
        end do
    end do
    return
end function complex_stagger_real_2d

subroutine zero_pad_signal_r1d(x,f,n,dx)
    implicit none
    real(num), allocatable, intent(inout) :: x(:), f(:)
    real(num), intent(out) :: dx
    integer, intent(in) :: n

    real(num), allocatable :: y(:), g(:)
    allocate(y(n), g(n))

    ! extend time domain, preserving step-size
    y = linspace(x(1), x(size(x)) + (n-size(x))*(x(2)-x(1)), n)
    call move_alloc(y,x)
    dx = x(2) - x(1)

    ! zero-pad signal
    g = 0.0_num
    g(1:size(f)) = f
    call move_alloc(g,f)

    return
end subroutine zero_pad_signal_r1d

subroutine zero_pad_signal_c1d(x,f,n,dx)
    implicit none
    real(num), allocatable, intent(inout) :: x(:)
    complex(num), allocatable, intent(inout) :: f(:)
    real(num), intent(out) :: dx
    integer, intent(in) :: n

    real(num), allocatable :: y(:)
    complex(num), allocatable :: g(:)
    allocate(y(n), g(n))

    ! extend time domain, preserving step-size
    y = linspace(x(1), x(size(x)) + (n-size(x))*(x(2)-x(1)), n)
    call move_alloc(y,x)
    dx = x(2) - x(1)

    ! zero-pad signal
    g = 0.0_num
    g(1:size(f)) = f
    call move_alloc(g,f)

    return
end subroutine zero_pad_signal_c1d

logical function is_pow2(n)
    implicit none
    integer, intent(in) :: n

    is_pow2 = .false.

    if (floor(log(real(n,num))/log(2.0_num)) == &
      ceiling(log(real(n,num))/log(2.0_num))) is_pow2 = .true.

    return
end function is_pow2

integer function next_pow2(n)
    implicit none
    integer, intent(in) :: n
    next_pow2 = 2**(ceiling(log(real(n,num))/log(2.0_num)))
    return
end function next_pow2

subroutine deriv_r1d(f,dx,df)
    implicit none
    real(num), intent(in) :: f(:), dx
    real(num), intent(out) :: df(:)

    integer :: n, j
    n = size(f)

    do j=1,2
        df(j) = -25*f(j) + 48*f(j+1) - 36*f(j+2) + 16*f(j+3) - 3*f(j+4)
    end do
    do j=3,n-2
        df(j) = -f(j+2) + 8*f(j+1) - 8*f(j-1) + f(j-2)
    end do
    do j=n-1,n
        df(j) = 25*f(j) - 48*f(j-1) + 36*f(j-2) - 16*f(j-3) + 3*f(j-4)
    end do

    df = df/(12.0*dx)

    return
end subroutine deriv_r1d

subroutine deriv_c1d(f,dx,df)
    implicit none
    complex(num), intent(in) :: f(:)
    real(num), intent(in) :: dx
    complex(num), intent(out) :: df(:)

    integer :: n, j
    n = size(f)

    do j=1,2
        df(j) = -25*f(j) + 48*f(j+1) - 36*f(j+2) + 16*f(j+3) - 3*f(j+4)
    end do
    do j=3,n-2
        df(j) = -f(j+2) + 8*f(j+1) - 8*f(j-1) + f(j-2)
    end do
    do j=n-1,n
        df(j) = 25*f(j) - 48*f(j-1) + 36*f(j-2) - 16*f(j-3) + 3*f(j-4)
    end do

    df = df/(12.0*dx)

    return
end subroutine deriv_c1d

subroutine deriv_r2d(f,dr,ax,df)
    implicit none
    real(num), intent(in) :: f(:,:), dr
    real(num), intent(out) :: df(:,:)
    integer, intent(in) :: ax

    integer :: nn(2), j
    nn = (/size(f,1),size(f,2)/)

    select case (ax)
    case (1)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(2)
            call deriv(f(:,j),dr,df(:,j))
        end do
        !$OMP END PARALLEL DO
    case (2)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(1)
            call deriv(f(j,:),dr,df(j,:))
        end do
        !$OMP END PARALLEL DO
    end select

    return
end subroutine deriv_r2d

subroutine deriv_c2d(f,dr,ax,df)
    implicit none
    complex(num), intent(in) :: f(:,:)
    real(num), intent(in) :: dr
    complex(num), intent(out) :: df(:,:)
    integer, intent(in) :: ax

    integer :: nn(2), j
    nn = (/size(f,1),size(f,2)/)

    select case (ax)
    case (1)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(2)
            call deriv(f(:,j),dr,df(:,j))
        end do
        !$OMP END PARALLEL DO
    case (2)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(1)
            call deriv(f(j,:),dr,df(j,:))
        end do
        !$OMP END PARALLEL DO
    end select

    return
end subroutine deriv_c2d

function grad_r1d(f,del) result(gradf)
    implicit none
    real(num), intent(in) :: f(:), del
    real(num), dimension(size(f)) :: gradf

    integer :: n, j
    n = size(f)

    gradf(1) = -3*f(1) + 4*f(2) - f(3)
    do j=2,n-1
        gradf(j) = f(j+1) - f(j-1)
    end do
    gradf(n) = 3*f(n) - 4*f(n-1) + f(n-2)

    gradf = gradf/(2.0*del)
    return
end function grad_r1d

function grad_c1d(f,del) result(gradf)
    implicit none
    complex(num), intent(in) :: f(:)
    real(num),    intent(in) :: del
    complex(num), dimension(size(f)) :: gradf
    gradf = grad(real(f),del) + i*grad(aimag(f),del)
    return
end function grad_c1d

function grad_r2d(f,del,ax) result(gradf)
    implicit none
    real(num), intent(in) :: f(:,:), del
    integer,   intent(in) :: ax

    real(num), dimension(size(f,1),size(f,2)) :: gradf
    integer :: nn(2), j

    nn = (/size(f,1),size(f,2)/)
    gradf = 0.0_num

    select case (ax)
    case (1)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(2)
            gradf(:,j) = grad(f(:,j),del)
        end do
        !$OMP END PARALLEL DO
    case (2)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(1)
            gradf(j,:) = grad(f(j,:),del)
        end do
        !$OMP END PARALLEL DO
    end select

    return
end function grad_r2d

function grad_c2d(f,del,ax) result(gradf)
    implicit none
    complex(num), intent(in) :: f(:,:)
    real(num),    intent(in) :: del
    integer,      intent(in) :: ax

    complex(num), dimension(size(f,1),size(f,2)) :: gradf
    integer :: nn(2), j

    nn = (/size(f,1),size(f,2)/)
    gradf = 0.0_num

    select case (ax)
    case (1)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(2)
            gradf(:,j) = grad(f(:,j),del)
        end do
        !$OMP END PARALLEL DO
    case (2)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(1)
            gradf(j,:) = grad(f(j,:),del)
        end do
        !$OMP END PARALLEL DO
    end select

    return
end function grad_c2d

function diff(x) result(dx)
    implicit none
    real(num), intent(in) :: x(:)
    real(num), dimension(size(x)-1) :: dx
    integer :: n, j
    n = size(x)
    ForAll(j=1:n-1) dx(j) = x(j+1) - x(j)
    return
end function diff

function simint(y,y0,dx) result(inty)
    implicit none
    real(num), intent(in) :: y(:), y0, dx
    real(num), dimension(size(y)) :: inty

    real(num) :: fac, term(5)
    integer   :: n, ip, j, k

    n = size(y)
    inty(1) = y0
    fac = dx/3.0

    term(1) = 1.25*y(1)
    term(2) = -0.25*y(1)
    k = 1
    do ip=3,n,2
        j = ip-1
        term(3) = 2.0*y(j)
        term(4) = 1.25*y(ip)
        term(5) = -0.25*y(ip)
        inty(j) = inty(k) + fac*(term(1)+term(3)+term(5))
        inty(ip) = inty(j) + fac*(term(2)+term(3)+term(4))
        term(1) = term(4)
        term(2) = term(5)
        k = ip
        inty(n) = inty(k) + fac*(-0.25*y(j) + 2.0*y(k) + 1.25*y(n))
    end do
    return
end function simint

function trapz_r1d(f) result(s)
    implicit none
    real(num), intent(in) :: f(:)

    real(num) :: s
    integer :: j

    s = 0.0_num
    do j=2,size(f)
        s = s + 0.5*(f(j-1) + f(j))
    end do
    return
end function trapz_r1d

function trapz_c1d(f) result(s)
    implicit none
    complex(num), intent(in) :: f(:)
    complex(num) :: s
    s = trapz(real(f)) + i*trapz(aimag(f))
    return
end function trapz_c1d

function trapz_r2d_part(f,ax) result(s)
    implicit none
    real(num), intent(in) :: f(:,:)
    integer, intent(in) :: ax

    real(num), dimension(size(f,ax)) :: s

    integer :: nn(2), j
    nn = (/size(f,1),size(f,2)/)

    select case (ax)
    case (1)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(1)
            s(j) = trapz(f(j,:))
        end do
        !$OMP END PARALLEL DO
    case (2)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nn(2)
            s(j) = trapz(f(:,j))
        end do
        !$OMP END PARALLEL DO
    end select

    return
end function trapz_r2d_part

function trapz_c2d_part(f,ax) result(s)
    implicit none
    complex(num), intent(in) :: f(:,:)
    integer, intent(in) :: ax
    complex(num), dimension(size(f,ax)) :: s
    s = trapz(real(f),ax) + i*trapz(aimag(f),ax)
    return
end function trapz_c2d_part

function trapz_r2d_full(f) result(s)
    implicit none
    real(num), intent(in) :: f(:,:)
    real(num) :: s

    real(num), dimension(size(f,1)) :: s1
    integer :: j

    do j=1,size(f,1)
        s1(j) = trapz(f(j,:))
    end do

    s = trapz(s1)
    return
end function trapz_r2d_full

function trapz_c2d_full(f) result(s)
    implicit none
    complex(num), intent(in) :: f(:,:)
    complex(num) :: s
    s = trapz(real(f)) + i*trapz(aimag(f))
    return
end function trapz_c2d_full

subroutine phase_unwrap_1d(f)
    implicit none
    real(num), intent(inout) :: f(:)

    real(num) :: delta
    integer :: n, j
    n = size(f)

    do j=2,n
        delta = f(j) - f(j-1)
        if (delta > pi) then
            f(j:) = f(j:) - 2.*pi
        else if (delta < -pi) then
            f(j:) = f(j:) + 2.*pi
        end if
    end do
    return
end subroutine phase_unwrap_1d

subroutine phase_unwrap_2d(f)
    implicit none
    real(num), intent(inout) :: f(:,:)

    integer :: nn(2), j
    nn = (/size(f,1),size(f,2)/)

    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=1,nn(1)
        call phase_unwrap(f(j,:))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(SHARED)
    do j=1,nn(2)
        call phase_unwrap(f(:,j))
    end do
    !$OMP END PARALLEL DO

    return
end subroutine phase_unwrap_2d

subroutine bcuint_r(y,y1,y2,y12,x1l,x2l,x1,x2,dx,ansy,ansy1,ansy2)
    implicit none
    real(num), intent(in) :: y(4), y1(4), y2(4), y12(4)
    real(num), intent(in) :: x1l, x2l, x1, x2, dx(2)
    real(num), intent(out) :: ansy, ansy1, ansy2

    real(num) :: x1p, x2p, v(16), a(4,4)
    integer :: n

    x1p = (x1-x1l)/dx(1)
    x2p = (x2-x2l)/dx(2)

    ! pack vector of data on surrounding 4 corners
    do n=1,4
        v(n) = y(n)
        v(n+4) = y1(n)*dx(1)
        v(n+8) = y2(n)*dx(2)
        v(n+12) = y12(n)*product(dx)
    end do

    ! matrix elements of interpolating polynomial coefficients
    a(1,1) = v(1)
    a(2,1) = v(5)
    a(3,1) = -3*v(1) + 3*v(2) - 2*v(5) - v(6)
    a(4,1) =  2*v(1) - 2*v(2) + v(5) + v(6)

    a(1,2) = v(9)
    a(2,2) = v(13)
    a(3,2) = -3*v(9) + 3*v(10) - 2*v(13) - v(14)
    a(4,2) =  2*v(9) - 2*v(10) + v(13) + v(14)

    a(1,3) = -3*v(1) + 3*v(4) - 2*v(9) - v(12)
    a(2,3) = -3*v(5) + 3*v(8) - 2*v(13) - v(16)
    a(3,3) = 9*v(1) - 9*v(2) + 9*v(3) - 9*v(4) + 6*v(5) + 3*v(6) - 3*v(7) - 6*v(8) &
        + 6*v(9) - 6*v(10) - 3*v(11) + 3*v(12) + 4*v(13) + 2*v(14) + v(15) + 2*v(16)
    a(4,3) = -6*v(1) + 6*v(2) - 6*v(3) + 6*v(4) - 3*v(5) - 3*v(6) + 3*v(7) + 3*v(8) &
        - 4*v(9) + 4*v(10) + 2*v(11) - 2*v(12) - 2*v(13) - 2*v(14) - v(15) - v(16)

    a(1,4) = 2*v(1) - 2*v(4) + v(9) + v(12)
    a(2,4) = 2*v(5) - 2*v(8) + v(13) + v(16)
    a(3,4) = -6*v(1) + 6*v(2) - 6*v(3) + 6*v(4) - 4*v(5) - 2*v(6) + 2*v(7) + 4*v(8) &
        - 3*v(9) + 3*v(10) + 3*v(11) - 3*v(12) - 2*v(13) - v(14) - v(15) - 2*v(16)
    a(4,4) = 4*v(1) - 4*v(2) + 4*v(3) - 4*v(4) + 2*v(5) + 2*v(6) - 2*v(7) - 2*v(8) &
        + 2*v(9) - 2*v(10) - 2*v(11) + 2*v(12) + v(13) + v(14) + v(15) + v(16)

    ! compute interpolated values
    ansy = dot_product((/1.0_num,x1p,x1p**2,x1p**3/), matmul(a, (/1.0_num,x2p,x2p**2,x2p**3/)))
    ansy1 = dot_product((/0.0_num,1.0_num,2.*x1p,3.*(x1p**2)/), matmul(a, (/1.0_num,x2p,x2p**2,x2p**3/)))/dx(1)
    ansy2 = dot_product((/1.0_num,x1p,x1p**2,x1p**3/), matmul(a, (/0.0_num,1.0_num,2.*x2p,3.*(x2p**2)/)))/dx(2)

    return
end subroutine bcuint_r

subroutine bcuint_r_old(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
    implicit none
    real(num), intent(in) :: y(4), y1(4), y2(4), y12(4)
    real(num), intent(in) :: x1l, x1u, x2l, x2u, x1, x2
    real(num), intent(out) :: ansy, ansy1, ansy2

    integer :: j
    real(num) :: t, u, c(4,4)

    call bcucof(y, y1, y2, y12, x1u-x1l, x2u-x2l, c)

    t = (x1-x1l)/(x1u-x1l)
    u = (x2-x2l)/(x2u-x2l)

    ansy  = 0.0_num
    ansy1 = 0.0_num
    ansy2 = 0.0_num

    do j=4,1,-1
        ansy  = t*ansy + ((c(j,4)*u + c(j,3))*u + c(j,2))*u + c(j,1)
        ansy1 = u*ansy1 + (3.*c(4,j)*t + 2.*c(3,j))*t + c(2,j)
        ansy2 = t*ansy2 + (3.*c(j,4)*u + 2.*c(j,3))*u + c(j,2)
    end do
    ansy1 = ansy1/(x1u-x1l)
    ansy2 = ansy2/(x2u-x2l)

    return
end subroutine bcuint_r_old

subroutine bcuint_c(y,y1,y2,y12,x1l,x2l,x1,x2,dx,ansy,ansy1,ansy2)
    implicit none
    complex(num), intent(in) :: y(4), y1(4), y2(4), y12(4)
    real(num),    intent(in) :: x1l, x2l, x1, x2, dx(2)
    complex(num), intent(out) :: ansy, ansy1, ansy2

    real(num) :: re_ansy,  im_ansy
    real(num) :: re_ansy1, im_ansy1
    real(num) :: re_ansy2, im_ansy2

    call bcuint(real(y,num), real(y1,num), real(y2,num), real(y12,num), &
        x1l, x2l, x1, x2, dx, re_ansy, re_ansy1, re_ansy2)
    call bcuint(aimag(y), aimag(y1), aimag(y2), aimag(y12), &
        x1l, x2l, x1, x2, dx, im_ansy, im_ansy1, im_ansy2)

    ansy  = re_ansy  + i*im_ansy
    ansy1 = re_ansy1 + i*im_ansy1
    ansy2 = re_ansy2 + i*im_ansy2

    return
end subroutine bcuint_c

subroutine bcucof(y,y1,y2,y12,d1,d2,c)
    implicit none
    real(num), intent(in) :: y(4), y1(4), y2(4), y12(4), d1, d2
    real(num), intent(out) :: c(4,4)

    real :: d1d2, xx, cl(16), wt(16,16), x(16)
    integer :: n, m, k, l

    SAVE wt
    DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,&
            10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,&
            4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,&
            10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,&
            0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,&
            10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,&
            5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,&
            10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

    d1d2 = d1*d2

    ! pack a temporary vector x
    do n=1,4
        x(n) = y(n)
        x(n+4) = y1(n)*d1
        x(n+8) = y2(n)*d2
        x(n+12) = y12(n)*d1d2
    end do

    ! matrix multiply by the stored table
    do n=1,16
        xx = 0.0
        do k=1,16
            xx = xx + wt(n,k)*x(k)
        end do
        cl(n) = xx
    end do

    ! unpack the result into the output table
    l = 0
    do n=1,4
        do m=1,4
            l = l+1
            c(n,m) = cl(l)
        end do
    end do

    return
end subroutine bcucof

recursive function LegendrePoly(n,x) result(r)
    implicit none
    integer, intent(in) :: n
    real(num), intent(in) :: x
    real(num) :: r

    select case (n)
    case (0)
        r = 1.0_num
    case (1)
        r = x
    case default
        r = (2*n-1)*x*LegendrePoly(n-1,x) - (n-1)*LegendrePoly(n-2,x)
        r = r/real(n,num)
    end select

    return
end function LegendrePoly

function LegendrePolySeq(n,x) result(s)
    implicit none
    integer, intent(in) :: n
    real(num), intent(in) :: x
    real(num), dimension(n) :: s
    integer :: j

    s(1) = 1.0_num
    s(2) = x

    do j=2,n-1
        s(j+1) = (2*j-1)*x*s(j) - (j-1)*s(j-1)
        s(j+1) = s(j+1)/real(j,num)
    end do

    return
end function LegendrePolySeq

pure recursive function factorial(n) result(r)
    implicit none
    integer, intent(in) :: n
    integer(kind=16) :: r
    select case (n)
    case (0:1)
        r = 1
    case default
        r = n*factorial(n-1)
    end select
    return
end function factorial

subroutine cache_factorial(r)
    implicit none
    integer(kind=16), intent(out) :: r(33)
    integer :: j
    r(1) = 1
    do j=2,33
        r(j) = j*r(j-1)
    end do
    return
end subroutine cache_factorial

real(num) function winHann(t, tau)
    implicit none
    real(num), intent(in) :: t, tau
    winHann = 0.5/tau * (1.0 - cos(2.*pi*t/tau))
    return
end function winHann

subroutine init_RNG()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: dt(8), pid, n, iu, iostat
    integer(int64) :: t

    call random_seed(size=n)
    allocate(seed(n))

    ! First query OS for random seed
    open(newunit=iu, file='/dev/urandom', access='stream', &
        form='unformatted', action='read', status='old', iostat=iostat)
    if (iostat == 0) then
        read(iu) seed
        close(iu)
    ! Fallback to XOR:ing the current time and pid
    ! useful for launching multiple instances in parallel
    else
        call system_clock(t)
        if (t == 0) then
            call date_and_time(values=dt)
            t = 365_int64*24*60*60*1000*(dt(1)-1970) &
                + 31_int64*24*60*60*1000*dt(2) &
                + 24_int64*60*60*1000*dt(3) &
                + 60*60*1000*dt(5) + 60*1000*dt(6) &
                + 1000*dt(7) + dt(8)
        end if
        pid = getpid()
        t = ieor(t, int(pid,kind(t)))
        seed(:) = lcg(t)
    end if
    call random_seed(put=seed)
    return

contains

    integer function lcg(s)
        use iso_fortran_env, only: int64
        implicit none
        integer(int64), intent(inout) :: s
        if (s == 0) then
            s = 104729
        else
            s = mod(s, 4294967296_int64)
        end if
        s = mod(s*279470273_int64, 4294967291_int64)
        lcg = int(mod(s, int(huge(0),int64)), kind(0))
        return
    end function lcg

end subroutine init_RNG

end module math
