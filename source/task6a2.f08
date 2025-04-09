program task6a2
    use prec
    use math
    use optimize
    use omp_lib

    implicit none

    ! parameters
    integer, parameter :: sdim = 8 ! 3-qubit state dimensionality
    integer, parameter :: r = 8 ! rank / no. chi angles
    real(num), parameter :: dsqrt3 = sqrt(3.0_num)
    ! user-set parameters
    integer :: k ! order+rank (min: 11)
    integer :: na ! no. theta & phi angles (=r*(2*k-r), min: 112)
    ! angle bounds
    real(num), parameter :: rng_t(2) = (/0.0_num, real(pi,num)/2./) ! $\theta$
    real(num), parameter :: rng_p(2) = (/-real(pi,num), real(pi,num)/) ! $\phi$
    real(num), parameter :: rng_c(2) = (/-real(pi,num), real(pi,num)/) ! $\chi$

    ! input parameters for optimization
    real(num) :: ftol, dt, dp, dc
    integer :: itmax, nmins

    ! global calculation variables
    complex(num) :: eigen_vals(r), eigen_vecs(r,sdim)
    complex(num), allocatable :: PSI(:,:), PSI_TILDE(:,:), A(:,:)
    real(num), allocatable :: prob(:)

    ! general variables
    type(OptimizeND_NelderMead) :: H
    real(num), allocatable :: My(:), Mp(:,:)
    real(num) :: time, Theta, cminval, tic, toc
    integer :: l

    call OMP_set_num_threads(OMP_get_max_threads())
    call OMP_set_nested(.false.)
    call init_RNG()
    call read_deck()

    if (k<11) then
        print '(A)', 'Recommended minimum k value (order+rank) is 11.'
        print '(A)', 'Please correct the input deck.'
        call exit(1)
    end if

    allocate(PSI(k,sdim), PSI_TILDE(k,sdim), A(k,r), prob(k))
    na = int(r*(2*k-r))

    ! cache eigenvalues/vectors
    eigen_vals = lambda()
    eigen_vecs = vex()

    call H%create(func=F, ndim=int(2*na+r), ftol=ftol, itmax=itmax, warn=.false.)

    allocate(My(H%ndim+1), Mp(H%ndim+1,H%ndim))
    cminval = huge(0.0_num)

    do l=1,nmins
        call H%reset()
        call init_opt(H)
        call H%minimize()
        if (minval(H%y) < cminval) then
            cminval = minval(H%y)
            My = H%y
            Mp = H%p
        end if
    end do

    ! tic = OMP_get_wtime()
    ! toc = OMP_get_wtime()
    ! print *, toc-tic, ' secs'

    call write_data()
    call exit()

contains

! OBJECTIVE FUNCTION -------------------------------------------------------------------------------

! $F(\rho)=\sum_i p_i F(\psi_i)$, input angles $(\vec{\theta}_{\rm[na]},\,\vec{\phi}_{\rm[na]},\,\vec{\chi}_{[r=8]})$
real(num) function F(x)
    implicit none
    real(num), intent(in) :: x(:)
    integer :: n, j

    F = 0.0_num
    PSI_TILDE = (0.0_num, 0.0_num)

    A = Amat(x(1:na),x(na+1:2*na),x(2*na+1:2*na+r))

    do n=1,k
    do j=1,r
        PSI_TILDE(n,:) = PSI_TILDE(n,:) + A(n,j)*sqrt(eigen_vals(j))*eigen_vecs(j,:)
    end do
    end do

    do n=1,k
        prob(n) = sum(real(conjg(PSI_TILDE(n,:))*PSI_TILDE(n,:),num))
        PSI(n,:) = PSI_TILDE(n,:)/sqrt(prob(n))
        F = F + prob(n)*FILL(PSI(n,:))
    end do

    return
end function F

! SUPPORTING ROUTINES ------------------------------------------------------------------------------

! initialize the $(\vec{\theta}_{\rm[na]},\,\vec{\phi}_{\rm[na]},\,\vec{\chi}_{[r=8]})$ angles randomly

subroutine init_opt(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    integer :: n

    this%iter = 0

    call random_number(this%p(1,:))
    ! 1st vertex, random guess in angle bounds
    this%p(1,1:na) = (rng_t(2)-rng_t(1))*this%p(1,1:na) + rng_t(1) ! $\vec{\theta}$
    this%p(1,na+1:2*na) = (rng_p(2)-rng_p(1))*this%p(1,na+1:2*na) + rng_p(1) ! $\vec{\phi}$
    this%p(1,2*na+1:2*na+r) = (rng_c(2)-rng_c(1))*this%p(1,2*na+1:2*na+r) + rng_c(1) ! $\vec{\chi}$
    ! function value
    this%y(1) = this%func(this%p(1,:))

    ! other vertices obtained by taking unit steps in each direction
    do n=2,this%ndim+1
        this%p(n,:) = this%p(1,:)

        if (n <= na+1) then
            this%p(n,n-1) = this%p(n,n-1) + dt
        else if ((n > na+1).and.(n <= 2*na+1)) then
            this%p(n,n-1) = this%p(n,n-1) + dp
        else
            this%p(n,n-1) = this%p(n,n-1) + dc
        end if

        ! function value
        this%y(n) = this%func(this%p(n,:))
    end do

    return
end subroutine init_opt

! SUB-CALCULATION FUNCTIONS ------------------------------------------------------------------------

! in 3-qubit-space, where (0,1)=(ground,excited), use 8-component state ordering
! Ord = (|000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>)

! see Task 6 notes for the density matrix description

! density matrix eigenvalues
function lambda()
    implicit none
    complex(num) :: lambda(r)
    real(num) :: p, q

    p = sqrt(exp(-time))
    q = sqrt(1.-exp(-time))

    lambda(1) = (p**4)*(q**2)*(cos(Theta)**2)
    lambda(2) = lambda(1)
    lambda(3) = (p**2)*(q**4)*(cos(Theta)**2)
    lambda(4) = lambda(3)

    lambda(5) = 0.5*(p**6)*(cos(Theta)**2) + 0.25*(p**2)*(1.+(q**4)+((q**4)-1.)*cos(2.*Theta)) &
        - 0.25*(p**2)*sqrt(-16.*(p**4)*(q**4)*(cos(Theta)**4) + (1.+(p**4)+(q**4)+((p**4)+(q**4)-1.)*cos(2.*Theta))**2)

    lambda(6) = 0.5*(p**6)*(cos(Theta)**2) + 0.25*(p**2)*(1.+(q**4)+((q**4)-1.)*cos(2.*Theta)) &
        + 0.25*(p**2)*sqrt(-16.*(p**4)*(q**4)*(cos(Theta)**4) + (1.+(p**4)+(q**4)+((p**4)+(q**4)-1.)*cos(2.*Theta))**2)

    lambda(7) = 0.5*(q**6)*(cos(Theta)**2) + 0.25*(q**2)*(1.+(p**4)+((p**4)-1.)*cos(2.*Theta)) &
        - 0.25*(q**2)*sqrt(-16.*(p**4)*(q**4)*(cos(Theta)**4) + (1.+(p**4)+(q**4)+((p**4)+(q**4)-1.)*cos(2.*Theta))**2)

    lambda(8) = 0.5*(q**6)*(cos(Theta)**2) + 0.25*(q**2)*(1.+(p**4)+((p**4)-1.)*cos(2.*Theta)) &
        + 0.25*(q**2)*sqrt(-16.*(p**4)*(q**4)*(cos(Theta)**4) + (1.+(p**4)+(q**4)+((p**4)+(q**4)-1.)*cos(2.*Theta))**2)

    return
end function lambda

! density matrix eigenvectors
function vex()
    implicit none
    complex(num) :: vex(r,sdim)
    real(num) :: p, q
    integer :: n

    p = sqrt(exp(-time))
    q = sqrt(1.-exp(-time))

    vex = (0.0_num, 0.0_num)

    vex(1,4) = -1.0_num
    vex(1,7) = 1.0_num

    vex(2,4) = -1.0_num
    vex(2,6) = 1.0_num

    vex(3,2) = -1.0_num
    vex(3,5) = 1.0_num

    vex(4,2) = -1.0_num
    vex(4,3) = 1.0_num

    vex(5,2) = (2.*p**2*(q**4-p**4)*cos(Theta)/sin(Theta) &
        - sqrt(cmplx(p**4*(-16.*p**4*q**4*cos(Theta)**4 &
        + (1.+p**4+q**4+(-1.+p**4+q**4)*cos(2.*Theta))**2)))/cos(Theta)/sin(Theta) &
        + 2.*p**2*tan(Theta))/(4.*dsqrt3*p**4)
    vex(5,3) = vex(5,2)
    vex(5,5) = vex(5,2)
    vex(5,8) = 1.0_num

    vex(6,2) = (2.*p**2*(q**4-p**4)*cos(Theta)/sin(Theta) &
        + sqrt(cmplx(p**4*(-16.*p**4*q**4*cos(Theta)**4 &
        + (1.+p**4+q**4+(-1.+p**4+q**4)*cos(2.*Theta))**2)))/cos(Theta)/sin(Theta) &
        + 2.*p**2*tan(Theta))/(4.*dsqrt3*p**4)
    vex(6,3) = vex(6,2)
    vex(6,5) = vex(6,2)
    vex(6,8) = 1.0_num

    vex(7,1) = -(dsqrt3*(-2.*q**6*cos(Theta)**2 &
        + q**2*(-1.+p**4+(1.+p**4)*cos(2.*Theta)) &
        + sqrt(cmplx(q**4*(-16.*p**4*q**4*cos(Theta)**4 &
        + (1.+p**4+q**4+(-1.+p**4+q**4)*cos(2.*Theta))**2))))/cos(Theta)/sin(Theta))/(4.*p**2*q**2 + epsilon(0.0_num))
    vex(7,4) = 1.0_num
    vex(7,6) = 1.0_num
    vex(7,7) = 1.0_num

    vex(8,1) = (dsqrt3*(2.*q**6*cos(Theta)**2 &
        - q**2*(-1.+p**4+(1.+p**4)*cos(2.*Theta)) &
        + sqrt(cmplx(q**4*(-16.*p**4*q**4*cos(Theta)**4 &
        + (1.+p**4+q**4+(-1.+p**4+q**4)*cos(2.*Theta))**2))))/cos(Theta)/sin(Theta))/(4.*p**2*q**2 + epsilon(0.0_num))
    vex(8,4) = 1.0_num
    vex(8,6) = 1.0_num
    vex(8,7) = 1.0_num

    ! normalization is critical !
    ForAll(n=1:r) &
        vex(n,:) = vex(n,:)/sqrt(sum(real(conjg(vex(n,:))*vex(n,:),num)))

    return
end function vex

! inverse of rotation matrices, input (n)umber (1 thru k-1) and (t)heta_i, (p)hi_i
function g(n,t,p)
    implicit none
    complex(num) :: g(k,k)
    integer, intent(in) :: n
    real(num), intent(in) :: t, p
    integer :: j

    g = (0.0_num, 0.0_num)
    ForAll(j=1:k) g(j,j) = 1.0_num

    g(n,n) = exp(-i*p)*cos(t)
    g(n,n+1) = -exp(-i*p)*sin(t)
    g(n+1,n) = exp(i*p)*sin(t)
    g(n+1,n+1) = exp(i*p)*cos(t)

    return
end function g

! k-by-r matrix of (c)hi's
function Rmat(c)
    implicit none
    complex(num) :: Rmat(k,r)
    real(num), intent(in) :: c(r)
    integer :: j

    Rmat = (0.0_num, 0.0_num)
    ForAll(j=1:r) Rmat(j,j) = exp(i*c(j))

    return
end function Rmat

! A matrix, see "notes/Task 5 - 65 parameters.pdf", Eq. (6).
! input vectors of (t)heta, (p)hi, and (c)hi angles
function Amat(t,p,c)
    implicit none
    complex(num) :: Amat(k,r)
    real(num), intent(in) :: t(na), p(na), c(r)
    integer :: j, n, m

    j = 0
    Amat = Rmat(c)

    do n=r,1,-1
    do m=n,k-1
        j = j + 1
        Amat = matmul(g(m,t(j),p(j)), Amat)
    end do
    end do

    return
end function Amat

! concurrence fill; input generic 3-qubit state (S)
real(num) function FILL(S)
    implicit none
    complex(num), intent(in) :: S(sdim)
    real(num) :: C(3), Qn

    C(1) = 4.*DET_RHO_1(S) ! C^2_{1(23)}
    C(2) = 4.*DET_RHO_2(S) ! C^2_{2(31)}
    C(3) = 4.*DET_RHO_3(S) ! C^2_{3(12)}

    Qn = 0.5*sum(C)

    FILL = ((16.*Qn/3.)*(Qn-C(1))*(Qn-C(2))*(Qn-C(3)))**(1./4.)

    return
end function FILL

! DET_RHO_x: partial traces for a generic 3-qubit state (S)

real(num) function DET_RHO_1(S)
    implicit none
    complex(num), intent(in) :: S(sdim)
    complex(num) :: temp

    temp = (abs(S(2))**2)*(abs(S(5))**2) - S(2)*S(5)*conjg(S(1))*conjg(S(6)) + &
           (abs(S(3))**2)*(abs(S(5))**2) - S(3)*S(5)*conjg(S(1))*conjg(S(7)) + &
           (abs(S(4))**2)*(abs(S(5))**2) - S(4)*S(5)*conjg(S(1))*conjg(S(8)) - &
           S(1)*S(6)*conjg(S(2))*conjg(S(5)) + (abs(S(1))**2)*(abs(S(6))**2) + &
           (abs(S(3))**2)*(abs(S(6))**2) - S(3)*S(6)*conjg(S(2))*conjg(S(7)) + &
           (abs(S(4))**2)*(abs(S(6))**2) - S(4)*S(6)*conjg(S(2))*conjg(S(8)) - &
           S(1)*S(7)*conjg(S(3))*conjg(S(5)) + (abs(S(1))**2)*(abs(S(7))**2) - &
           S(2)*S(7)*conjg(S(3))*conjg(S(6)) + (abs(S(2))**2)*(abs(S(7))**2) + &
           (abs(S(4))**2)*(abs(S(7))**2) - S(4)*S(7)*conjg(S(3))*conjg(S(8)) - &
           S(1)*S(8)*conjg(S(4))*conjg(S(5)) + (abs(S(1))**2)*(abs(S(8))**2) - &
           S(2)*S(8)*conjg(S(4))*conjg(S(6)) + (abs(S(2))**2)*(abs(S(8))**2) - &
           S(3)*S(8)*conjg(S(4))*conjg(S(7)) + (abs(S(3))**2)*(abs(S(8))**2);

    DET_RHO_1 = real(temp,num) ! Im ~ 0

    return
end function DET_RHO_1

real(num) function DET_RHO_2(S)
    implicit none
    complex(num), intent(in) :: S(sdim)
    complex(num) :: temp

    temp = (abs(S(2))**2)*(abs(S(3))**2) - S(2)*S(3)*conjg(S(1))*conjg(S(4)) - &
           S(1)*S(4)*conjg(S(2))*conjg(S(3)) + (abs(S(1))**2)*(abs(S(4))**2) + &
           (abs(S(3))**2)*(abs(S(5))**2) - S(3)*S(5)*conjg(S(1))*conjg(S(7)) + &
           (abs(S(4))**2)*(abs(S(5))**2) - S(4)*S(5)*conjg(S(2))*conjg(S(7)) + &
           (abs(S(3))**2)*(abs(S(6))**2) - S(3)*S(6)*conjg(S(1))*conjg(S(8)) + &
           (abs(S(4))**2)*(abs(S(6))**2) - S(4)*S(6)*conjg(S(2))*conjg(S(8)) - &
           S(1)*S(7)*conjg(S(3))*conjg(S(5)) + (abs(S(1))**2)*(abs(S(7))**2) - &
           S(2)*S(7)*conjg(S(4))*conjg(S(5)) + (abs(S(2))**2)*(abs(S(7))**2) + &
           (abs(S(6))**2)*(abs(S(7))**2) - S(6)*S(7)*conjg(S(5))*conjg(S(8)) - &
           S(1)*S(8)*conjg(S(3))*conjg(S(6)) + (abs(S(1))**2)*(abs(S(8))**2) - &
           S(2)*S(8)*conjg(S(4))*conjg(S(6)) + (abs(S(2))**2)*(abs(S(8))**2) - &
           S(5)*S(8)*conjg(S(6))*conjg(S(7)) + (abs(S(5))**2)*(abs(S(8))**2);

    DET_RHO_2 = real(temp,num) ! Im ~ 0

    return
end function DET_RHO_2

real(num) function DET_RHO_3(S)
    implicit none
    complex(num), intent(in) :: S(sdim)
    complex(num) :: temp

    temp = (abs(S(2))**2)*(abs(S(3))**2) - S(2)*S(3)*conjg(S(1))*conjg(S(4)) - &
           S(1)*S(4)*conjg(S(2))*conjg(S(3)) + (abs(S(1))**2)*(abs(S(4))**2) + &
           (abs(S(2))**2)*(abs(S(5))**2) - S(2)*S(5)*conjg(S(1))*conjg(S(6)) + &
           (abs(S(4))**2)*(abs(S(5))**2) - S(4)*S(5)*conjg(S(3))*conjg(S(6)) - &
           S(1)*S(6)*conjg(S(2))*conjg(S(5)) + (abs(S(1))**2)*(abs(S(6))**2) - &
           S(3)*S(6)*conjg(S(4))*conjg(S(5)) + (abs(S(3))**2)*(abs(S(6))**2) + &
           (abs(S(2))**2)*(abs(S(7))**2) - S(2)*S(7)*conjg(S(1))*conjg(S(8)) + &
           (abs(S(4))**2)*(abs(S(7))**2) - S(4)*S(7)*conjg(S(3))*conjg(S(8)) + &
           (abs(S(6))**2)*(abs(S(7))**2) - S(6)*S(7)*conjg(S(5))*conjg(S(8)) - &
           S(1)*S(8)*conjg(S(2))*conjg(S(7)) + (abs(S(1))**2)*(abs(S(8))**2) - &
           S(3)*S(8)*conjg(S(4))*conjg(S(7)) + (abs(S(3))**2)*(abs(S(8))**2) - &
           S(5)*S(8)*conjg(S(6))*conjg(S(7)) + (abs(S(5))**2)*(abs(S(8))**2);

    DET_RHO_3 = real(temp,num) ! Im ~ 0

    return
end function DET_RHO_3

! DATA I/O ROUTINES --------------------------------------------------------------------------------

subroutine read_deck
    implicit none
    character(len=128) :: label
    character(len=3) :: unit
    integer :: iu, n, m
    open(newunit=iu, file='input-T6.deck', status='old', action='read')
    read(iu,*) ! density matrix
    read(iu,*) label, time
    read(iu,*) label, Theta, unit ! unit: 'pi' or 'rad'
    if (trim(adjustl(unit)) == 'pi') Theta = pi*Theta
    read(iu,*) label, k ! order+rank
    read(iu,*)
    read(iu,*) ! initialization parameters
    read(iu,*) label, nmins
    read(iu,*) label, dt; dt = pi*dt;
    read(iu,*) label, dp; dp = pi*dp;
    read(iu,*) label, dc; dc = pi*dc;
    read(iu,*)
    read(iu,*) ! termination parameters
    read(iu,*) label, ftol
    read(iu,*) label, itmax
    close(iu)
    return
end subroutine read_deck

subroutine write_data
    implicit none
    integer :: n, iu

    open(newunit=iu, file='yvals.dat', status='replace')
    do n=1,H%ndim+1
        write(iu,*) My(n)
    end do
    close(iu)

    open(newunit=iu, file='pvecs.dat', status='replace')
    do n=1,H%ndim+1
        write(iu,*) Mp(n,:)
    end do
    close(iu)

    return
end subroutine write_data

end program task6a2
