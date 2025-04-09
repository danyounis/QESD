program task2
    use prec
    use math
    use optimize
    use omp_lib

    implicit none

    ! parameters
    ! adim: 14-angle parametrization, sdim: 3-qubit state dimensionality
    integer, parameter :: adim = 14, sdim = 8
    real(num), parameter :: dsqrt2 = sqrt(2.0_num)
    complex(num), parameter :: ReNum = (1.0_num, 0.0_num)

    ! general variables
    real(num) :: time, tic, toc
    ! r: rho 8-tuple, v: witness operator 8-tuple
    real(num) :: r(sdim), v(sdim)

    ! extra variables
    real(num) :: norm, grad_fin(sdim+1,sdim)

    type(OptimizeND_NelderMead) :: G
    type(OptimizeND_NelderMead) :: H

    integer :: n, iu

    call OMP_set_num_threads(OMP_get_max_threads())
    call OMP_set_nested(.false.)
    call init_RNG()

    open(newunit=iu, file='time.param', status='old', action='read')
    read(iu,*) time
    close(iu)

    ! set $\rho$ 8-tuple
    r = RHO(time)

    call G%create(func=Gmax, ndim=sdim, ftol=real(1.d-6,num), itmax=-1, warn=.false.)

    ! G%aux1: optimal 14-angle for ea. Nelder-Mead vertex
    ! G%aux2: assists G%aux1 in assignment to the correct location
    allocate(G%aux1(sdim+1,adim), G%aux2(1,adim))

    call H%create(func=Hpsi, ndim=adim, ftol=real(1.d-6,num), itmax=-1, warn=.false.)

    call init_exterior(G)
    tic = OMP_get_wtime()
    call G%minimize()
    toc = OMP_get_wtime()

    ! OUTPUT DATA !

    open(newunit=iu, file='feval.dat', status='replace')
    do n=1,sdim+1
        write(iu,*) -G%y(n)
    end do
    close(iu)

    open(newunit=iu, file='vstar.dat', status='replace')
    do n=1,sdim+1
        write(iu,*) G%p(n,:)
    end do
    close(iu)

    open(newunit=iu, file='astar.dat', status='replace')
    do n=1,sdim+1
        write(iu,*) G%aux1(n,:)
    end do
    close(iu)

    open(newunit=iu, file='grad.dat', status='replace')
    do n=1,sdim+1
        grad_fin(n,:) = r - TrPmatPSI(PSI(G%aux1(n,1:3),G%aux1(n,4:9),G%aux1(n,10:13),G%aux1(n,14)))
        write(iu,*) grad_fin(n,:)
    end do
    close(iu)

    open(newunit=iu, file='norm.dat', status='replace')
    do n=1,sdim+1
        write(iu,*) norm2(grad_fin(n,:))
    end do
    close(iu)

    open(newunit=iu, file='clock.dat', status='replace')
    write(iu,*) toc-tic
    close(iu)

    call exit()

contains

! OBJECTIVE FUNCTIONS ------------------------------------------------------------------------------

! $H(\psi|\mathbf{v}) = \mathbf{v}\cdot(\mathbf{r}-\mathbf{p}) + E(\psi)$; input $\Psi$ angles (x)
real(num) function Hpsi(x)
    implicit none
    real(num), intent(in) :: x(:)
    real(num) :: p(sdim)

    p = TrPmatPSI(PSI(x(1:3),x(4:9),x(10:13),x(14)))
    Hpsi = sum(v*(r-p)) + F123(x(10:13),x(14))

    return
end function Hpsi

function grad_Hpsi(x)
    implicit none
    real(num), intent(in) :: x(:)
    real(num) :: grad_Hpsi(size(x))
    real(num) :: err
    integer :: j

    do j=1,adim
        grad_Hpsi(j) = dfridr(Hpsi,j,x,real(1.d-2,num),err)
    end do

    return
end function grad_Hpsi

! $G(\mathbf{w}) = \min_\psi H(\psi|\mathbf{w})$; input witness operator 8-tuple

real(num) function Gmax(w)
    implicit none
    real(num), intent(in) :: w(:)
    integer, parameter :: nmins = 500
    real(num) :: My(nmins), Mp(nmins,adim) ! thread master variables
    integer :: n, k

    ! set the current witness operator
    v = w

    ! parallel interior minimizations to thoroughly explore the space
    !$OMP PARALLEL DO SHARED(My,Mp) PRIVATE(k) FIRSTPRIVATE(H)
    do n=1,nmins
        call init_interior(H)
        call H%minimize()

        k = minloc(H%y,1)
        My(n) = H%y(k)
        Mp(n,:) = H%p(k,:)
    end do
    !$OMP END PARALLEL DO

    k = minloc(My,1)
    Gmax = -My(k)
    G%aux2(1,:) = Mp(k,:)

    return
end function Gmax

! SUPPORTING ROUTINES ------------------------------------------------------------------------------

! initialize the ($\theta$,$\phi$,$\alpha$,$\beta$) angles randomly

subroutine init_interior(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    integer :: j

    this%iter = 0

    do j=1,adim+1
        call random_number(this%p(j,:))
        this%p(j,1:3) = (pi/2.)*this%p(j,1:3) ! $\theta$
        this%p(j,4:9) = (2.*pi)*this%p(j,4:9) ! $\phi$
        this%p(j,10:13) = (pi/2.)*this%p(j,10:13) ! $\alpha$
        this%p(j,14) = pi*this%p(j,14) ! $\beta$
        this%y(j) = this%func(this%p(j,:)) ! function value
    end do

    return
end subroutine init_interior

! initialize the $\vec{v}$ components randomly

subroutine init_exterior(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    integer :: j

    this%iter = 0

    do j=1,sdim+1
        call random_number(this%p(j,:))
        this%p(j,:) = 2.d+3*this%p(j,:) - 1.d+3 ! $\vec{v}$
        this%y(j) = this%func(this%p(j,:)) ! function value
    end do

    return
end subroutine init_exterior

! SUB-CALCULATION FUNCTIONS ------------------------------------------------------------------------

! in 3-qubit-space, where (0,1)=(ground,excited), use 8-component state ordering
! Ord = (|000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>)

! RHO(t) = t*|GHZ><GHZ| + (1-t)*|W><W|
! expanded as an 8-tuple in terms of the basis P-matrices
function RHO(t)
    implicit none
    real(num) :: RHO(sdim)
    real(num), intent(in) :: t

    RHO = 0.0_num

    RHO(1) = t/dsqrt2
    RHO(2) = t/dsqrt2
    RHO(5) = 1.-t

    return
end function RHO

! 3-qubit state parametrized by 3 (t)heta's, 6 (p)hi's, 4 (a)lpha's, and 1 (b)eta angle
function PSI(t,p,a,b)
    implicit none
    complex(num) :: PSI(sdim)
    real(num), intent(in) :: t(3), p(6), a(4), b

    complex(num), dimension(2,2) :: M1, M2, M3
    complex(num), dimension(8,8) :: U

    M1(1,:) = (/ exp(+i*p(1))*cos(t(1)), exp(+i*p(2))*sin(t(1))/)
    M1(2,:) = (/-exp(-i*p(2))*sin(t(1)), exp(-i*p(1))*cos(t(1))/)

    M2(1,:) = (/ exp(+i*p(3))*cos(t(2)), exp(+i*p(4))*sin(t(2))/)
    M2(2,:) = (/-exp(-i*p(4))*sin(t(2)), exp(-i*p(3))*cos(t(2))/)

    M3(1,:) = (/ exp(+i*p(5))*cos(t(3)), exp(+i*p(6))*sin(t(3))/)
    M3(2,:) = (/-exp(-i*p(6))*sin(t(3)), exp(-i*p(5))*cos(t(3))/)

    U = Kproduct(Kproduct(M1,M2),M3)

    PSI = matmul(U,(/cos(a(1))*ReNum, 0.*ReNum, 0.*ReNum, 0.*ReNum, &
        sin(a(1))*cos(a(2))*exp(i*b), &
        sin(a(1))*sin(a(2))*cos(a(3))*ReNum, &
        sin(a(1))*sin(a(2))*sin(a(3))*cos(a(4))*ReNum, &
        sin(a(1))*sin(a(2))*sin(a(3))*sin(a(4))*ReNum/))

    return
end function PSI

! returns p(i) = trace over |PSI><PSI| times the i-th basis P-matrix; input generic 3-qubit state (S)
function TrPmatPSI(S) result(p)
    implicit none
    real(num) :: p(sdim)
    complex(num), intent(in) :: S(sdim)

    ! Im ~ 0
    p(1) = real((abs(S(1))**2 + abs(S(8))**2)/dsqrt2, num)
    p(2) = real((conjg(S(8))*S(1) + conjg(S(1))*S(8))/dsqrt2, num)
    p(3) = real(i*(conjg(S(8))*S(1) - conjg(S(1))*S(8))/dsqrt2, num)
    p(4) = real((abs(S(1))**2 - abs(S(8))**2)/dsqrt2, num)
    p(5) = real((conjg(S(2)) + conjg(S(3)) + conjg(S(5)))*(S(2) + S(3) + S(5))/3., num)
    p(6) = real((conjg(S(5))*(2.*S(5) - S(3) - S(2)) + &
                  conjg(S(2))*(2.*S(2) - S(5) - S(3)) + &
                  conjg(S(3))*(2.*S(3) - S(5) - S(2)))/3./dsqrt2, num)
    p(7) = real((conjg(S(4)) + conjg(S(6)) + conjg(S(7)))*(S(4) + S(6) + S(7))/3., num)
    p(8) = real((conjg(S(7))*(2.*S(7) - S(6) - S(4)) + &
                  conjg(S(4))*(2.*S(4) - S(7) - S(6)) + &
                  conjg(S(6))*(2.*S(6) - S(7) - S(4)))/3./dsqrt2, num)

    return
end function TrPmatPSI

! concurrence fill; input 4 (a)lpha and 1 (b)eta angle
real(num) function F123(a,b)
    implicit none
    real(num), intent(in) :: a(4), b
    real(num) :: x, y, z
    real(num) :: c(4), s(4), c2(4), s2(4)

    c = (/cos(a(1)), cos(a(2)), cos(a(3)), cos(a(4))/)
    s = (/sin(a(1)), sin(a(2)), sin(a(3)), sin(a(4))/)
    c2 = (/cos(a(1))**2, cos(a(2))**2, cos(a(3))**2, cos(a(4))**2/)
    s2 = (/sin(a(1))**2, sin(a(2))**2, sin(a(3))**2, sin(a(4))**2/)

    x = (sin(2*a(1))**2)*s2(2)
    y = 4*s2(1)*s2(2)*s2(3)*&
        (c2(1) + s2(1)*(c2(3)*c2(4)*s2(2) - 2*c(2)*c(3)*c(4)*s(2)*s(4)*cos(b) + c2(2)*s2(4)))
    z = 4*s2(1)*s2(2)*(s2(1)*s2(3)*&
        (c2(3)*c2(4)*s2(2) - 2*c(2)*c(3)*c(4)*s(2)*s(4)*cos(b) + c2(2)*s2(4)) + c2(1)*(c2(3) + s2(3)*s2(4)))

    F123 = ((1./3._num)*(x+y+z)*(-x+y+z)*(x-y+z)*(x+y-z))**(1./4._num)

    return
end function F123

end program task2
