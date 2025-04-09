program task5m
    use prec
    use math
    use optimize
    use omp_lib

    implicit none

    ! parameters
    integer, parameter :: adim = 14 ! angle parametrization
    integer, parameter :: sdim = 8 ! 3-qubit state dimensionality
    integer, parameter :: vdim = 6 ! witness operator dimensionality
    real(num), parameter :: dsqrt2 = sqrt(2.0_num)
    real(num), parameter :: dsqrt3 = sqrt(3.0_num)
    real(num), parameter :: dsqrt6 = sqrt(6.0_num)
    complex(num), parameter :: ReNum = (1.0_num, 0.0_num)

    ! input parameters for (e)xterior & (i)nterior optimization
    real(num) :: ftol_e, ftol_i
    real(num) :: vrange(2), wigperc, alpha
    integer :: itmax_e, itmax_i, nmaxs, nmins
    logical :: srestart
    character(len=2) :: imethod

    ! general variables
    real(num) :: time, tic, toc
    integer :: l, lm=1
    ! r: rho 6-tuple, v: witness operator 6-tuple
    real(num) :: r(vdim), v(vdim)

    ! interior-step thread master variables
    real(num), allocatable :: My(:), Mp(:,:)

    type(OptimizeND_NelderMead), allocatable :: G(:)
    type(OptimizeND_NelderMead) :: Hnm
    type(OptimizeND_ConjugateGradient) :: Hcg

    call OMP_set_num_threads(OMP_get_max_threads())
    call OMP_set_nested(.false.)
    call init_RNG()
    call read_deck

    allocate(G(nmaxs), My(nmins), Mp(nmins,adim))

    ! set $\rho$ 8-tuple
    r = RHO(time)

    do l=1,nmaxs
        call G(l)%create(func=Gmax, ndim=vdim, ftol=ftol_e, itmax=itmax_e, warn=.false.)

        ! G(l)%aux1: optimal angles for ea. Nelder-Mead vertex
        ! G(l)%aux2: assists G%aux1 in assignment to the correct location
        allocate(G(l)%aux1(vdim+1,adim), G(l)%aux2(1,adim))
    end do

    select case (imethod)
    case ('NM')
        call Hnm%create(func=Hpsi, ndim=adim, ftol=ftol_i, itmax=itmax_i, warn=.false.)
    case ('CG')
        call Hcg%create(func=Hpsi, gfunc=grad_Hpsi, alpha=alpha, &
            ndim=adim, ftol=ftol_i, itmax=itmax_i, warn=.false.)
    end select

    tic = OMP_get_wtime()
    do l=1,nmaxs
        call init_exterior(G(l))
        call G(l)%minimize()

        ! output this iteration's entanglement measure
        ! if it exceeds the previously-known maximum
        if ( (l==1) .or. ((l>1).and.(maxval(-G(l)%y) > maxval(-G(lm)%y))) ) then
            lm = l
            call write_data
        end if
    end do
    toc = OMP_get_wtime()

    call exit()

contains

! OBJECTIVE FUNCTIONS ------------------------------------------------------------------------------

! $H(\psi|\mathbf{v}) = \mathbf{v}\cdot(\mathbf{r}-\mathbf{p}) + E(\psi)$; input $\Psi$ angles (x)
real(num) function Hpsi(x)
    implicit none
    real(num), intent(in) :: x(:)
    real(num) :: p(vdim)

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

! $G(\mathbf{v}) = \min_\psi H(\psi|\mathbf{v})$; input witness operator 8-tuple

real(num) function Gmax(x)
    implicit none
    real(num), intent(in) :: x(:)
    integer :: n, k

    ! set the current witness operator
    v = x

    ! parallel interior minimizations to thoroughly explore the space
    select case (imethod)
    case ('NM')

        !$OMP PARALLEL DO SHARED(My,Mp) PRIVATE(k) FIRSTPRIVATE(Hnm)
        do n=1,nmins
            call init_interior_nm(Hnm)
            call Hnm%minimize()

            k = minloc(Hnm%y,1)
            My(n) = Hnm%y(k)
            Mp(n,:) = Hnm%p(k,:)
        end do
        !$OMP END PARALLEL DO

    case ('CG')

        !$OMP PARALLEL DO SHARED(My,Mp) FIRSTPRIVATE(Hcg)
        do n=1,nmins
            call init_interior_cg(Hcg)
            call Hcg%minimize()

            My(n) = Hcg%y
            Mp(n,:) = Hcg%p
        end do
        !$OMP END PARALLEL DO

    end select

    k = minloc(My,1)
    Gmax = -My(k)
    G(l)%aux2(1,:) = Mp(k,:)

    return
end function Gmax

! SUPPORTING ROUTINES ------------------------------------------------------------------------------

! initialize the ($\theta$,$\phi$,$\alpha$,$\beta$) angles randomly

subroutine init_interior_nm(this)
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
end subroutine init_interior_nm

subroutine init_interior_cg(this)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(inout) :: this

    this%iter = 0

    call random_number(this%p)
    this%p(1:3) = (pi/2.)*this%p(1:3) ! $\theta$
    this%p(4:9) = (2.*pi)*this%p(4:9) ! $\phi$
    this%p(10:13) = (pi/2.)*this%p(10:13) ! $\alpha$
    this%p(14) = pi*this%p(14) ! $\beta$
    this%y = this%func(this%p) ! function value

    return
end subroutine init_interior_cg

! initialize the $\vec{v}$ components randomly (but within the user-defined range)

subroutine init_exterior(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    real(num) :: oldp(vdim)
    integer :: j, iu
    logical :: exist

    this%iter = 0

    if (srestart) then

        inquire(file='vstar.dat', exist=exist)
        if (.not.exist) then
            print '(A)', 'Error: To soft restart, I need a vstar.dat file.'
            print '(A)', 'Terminating execution.'
            call exit(1)
        end if

        open(newunit=iu, file='vstar.dat', status='old', action='read')
        read(iu,*) oldp(:)
        close(iu)

        do j=1,vdim+1
            call random_number(this%p(j,:))
            this%p(j,:) = oldp + wigperc*(2.*this%p(j,:)-1.)/100. ! $\vec{v}$
            this%y(j) = this%func(this%p(j,:)) ! function value
        end do

    else

        do j=1,vdim+1
            call random_number(this%p(j,:))
            this%p(j,:) = (vrange(2)-vrange(1))*this%p(j,:) + vrange(1) ! $\vec{v}$
            this%y(j) = this%func(this%p(j,:)) ! function value
        end do

    end if

    return
end subroutine init_exterior

! SUB-CALCULATION FUNCTIONS ------------------------------------------------------------------------

! in 3-qubit-space, where (0,1)=(ground,excited), use 8-component state ordering
! Ord = (|000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>)

! see Task 5 notes for the density matrix description
function RHO(t)
    implicit none
    real(num) :: RHO(vdim)
    real(num), intent(in) :: t
    real(num) :: p, q

    p = sqrt(exp(-t))
    q = sqrt(1.-exp(-t))

    RHO(1) = q**4
    RHO(2) = 2.*(p**2)*(q**2)/dsqrt3
    RHO(3) = (p**4)/dsqrt3
    RHO(4) = 0.
    RHO(5) = dsqrt2*(p**2)*(q**2)/dsqrt3
    RHO(6) = dsqrt2*(p**4)/dsqrt3

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
    real(num) :: p(vdim)
    complex(num), intent(in) :: S(sdim)

    ! Im ~ 0
    p(1) = real(abs(S(1))**2, num)
    p(2) = real((abs(S(2))**2 + abs(S(3))**2 + abs(S(5))**2)/dsqrt3, num)
    p(3) = real((abs(S(4))**2 + abs(S(6))**2 + abs(S(7))**2)/dsqrt3, num)
    p(4) = real(abs(S(8))**2, num)
    p(5) = real((conjg(S(2))*(S(3)+S(5)) + conjg(S(3))*(S(2)+S(5)) + conjg(S(5))*(S(2)+S(3)))/dsqrt6, num)
    p(6) = real((conjg(S(4))*(S(6)+S(7)) + conjg(S(6))*(S(4)+S(7)) + conjg(S(7))*(S(4)+S(6)))/dsqrt6, num)

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

! DATA I/O ROUTINES --------------------------------------------------------------------------------

subroutine read_deck
    implicit none
    character(len=128) :: label
    character(len=3) :: unit
    integer :: iu
    open(newunit=iu, file='input-T5m.deck', status='old', action='read')
    read(iu,*) ! density matrix
    read(iu,*) label, time
    read(iu,*)
    read(iu,*) ! exterior optimization (Nelder-Mead)
    read(iu,*) label, nmaxs
    read(iu,*) label, ftol_e
    read(iu,*) label, itmax_e
    read(iu,*) label, vrange(1), vrange(2) ! range to initialize $\vec{v}$
    read(iu,*)
    read(iu,*) ! soft restart?
    read(iu,*) label, srestart
    read(iu,*) label, wigperc ! percent ($\pm$) to wiggle $\vec{v}*$
    read(iu,*)
    read(iu,*) ! interior optimization
    read(iu,*) label, imethod ! NM/Nelder-Mead or CG/Conjugate-Gradient
    read(iu,*) label, nmins
    read(iu,*) label, ftol_i
    read(iu,*) label, itmax_i
    read(iu,*) label, alpha
    close(iu)
    return
end subroutine read_deck

subroutine write_data
    implicit none
    real(num) :: grad(vdim+1,vdim)
    integer :: n, iu

    open(newunit=iu, file='feval.dat', status='replace')
    do n=1,vdim+1
        write(iu,*) -G(lm)%y(n)
    end do
    close(iu)

    open(newunit=iu, file='vstar.dat', status='replace')
    do n=1,vdim+1
        write(iu,*) G(lm)%p(n,:)
    end do
    close(iu)

    open(newunit=iu, file='astar.dat', status='replace')
    do n=1,vdim+1
        write(iu,*) G(lm)%aux1(n,:)
    end do
    close(iu)

    open(newunit=iu, file='grad.dat', status='replace')
    do n=1,vdim+1
        grad(n,:) = r - TrPmatPSI(PSI(G(lm)%aux1(n,1:3),G(lm)%aux1(n,4:9),G(lm)%aux1(n,10:13),G(lm)%aux1(n,14)))
        write(iu,*) grad(n,:)
    end do
    close(iu)

    open(newunit=iu, file='norm.dat', status='replace')
    do n=1,vdim+1
        write(iu,*) norm2(grad(n,:))
    end do
    close(iu)

    return
end subroutine write_data

end program task5m
