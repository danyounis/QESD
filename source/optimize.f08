! Optimization Module (optimize.f08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 12/7/2021
! Revised: 11/29/2022
!
! References:
! [1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
! Numerical Recipes in Fortran 90 (Cambridge University Press, Cambridge, 2001).

module optimize

use prec

implicit none

type, public :: OptimizeND_NelderMead ! (NM)
    procedure(vsfunc), pointer, nopass :: func => null()
    real(num), allocatable :: y(:), p(:,:)
    real(num), allocatable :: aux1(:,:), aux2(:,:)
    real(num) :: ftol
    integer :: ndim, iter, itmax
    logical :: warn
contains
! core
    procedure :: create => create_NM
    procedure :: minimize => amoeba
    procedure :: destroy => destroy_NM
    procedure :: reset => reset_NM
! support
    procedure :: amotry => amotry_NM
end type OptimizeND_NelderMead

type, public :: OptimizeND_ConjugateGradient ! (CG)
    procedure(vsfunc), pointer, nopass :: func => null()
    procedure(vvfunc), pointer, nopass :: gfunc => null()
    real(num), allocatable :: xi(:), p(:)
    real(num) :: y, alpha, ftol
    integer :: ndim, iter, itmax
    logical :: warn
contains
! core
    procedure :: create => create_CG
    procedure :: minimize => frprmn
    procedure :: destroy => destroy_CG
    procedure :: reset => reset_CG
! support
    procedure :: linmin => linmin_CG
    procedure :: mnbrak => mnbrak_CG
    procedure :: dbrent => dbrent_CG
end type OptimizeND_ConjugateGradient

type, public :: OptimizeND_ParticleSwarm ! (PS)
    procedure(vsfunc), pointer, nopass :: func => null()
    real(num), allocatable :: x(:,:), v(:,:)
    real(num), allocatable :: xbest(:,:), func_xbest(:)
    real(num), allocatable :: p(:)
    real(num) :: y, ftol
    real(num) :: w, c(2)
    integer :: ndim, npart, iter, itmax
    logical :: warn
contains
! core
    procedure :: create => create_PS
    procedure :: minimize => nemo
    procedure :: destroy => destroy_PS
    procedure :: reset => reset_PS
! support
    procedure :: span => flockspan_PS
end type OptimizeND_ParticleSwarm

type, public :: OptimizeND_GuidedMonteCarlo ! (GMC)
    procedure(vsfunc), pointer, nopass :: func => null()
    real(num), allocatable, dimension(:) :: p, mag_step, mag_delta
    real(num) :: y, step(0:20)
    integer :: ndim, iter, itmax, iroc
    logical :: warn
contains
! core
    procedure :: create => create_GMC
    procedure :: minimize => gmcmn
    procedure :: destroy => destroy_GMC
    procedure :: reset => reset_GMC
! support
end type OptimizeND_GuidedMonteCarlo

abstract interface
    ! scalar->scalar function
    function ssfunc(x)
        use prec, only: num
        implicit none
        real(num), intent(in) :: x
        real(num) :: ssfunc
    end function ssfunc

    ! vector->scalar function
    function vsfunc(x)
        use prec, only: num
        implicit none
        real(num), intent(in) :: x(:)
        real(num) :: vsfunc
    end function vsfunc

    ! vector->vector function
    function vvfunc(x)
        use prec, only: num
        implicit none
        real(num), intent(in) :: x(:)
        real(num) :: vvfunc(size(x))
    end function vvfunc
end interface

interface dfridr
    module procedure dfridr_ss
    module procedure dfridr_vs
end interface dfridr

contains

! create_ subroutines ----------------------------------------------------------

subroutine create_NM(this,func,ndim,ftol,itmax,warn)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    procedure(vsfunc) :: func
    real(num), intent(in) :: ftol
    integer, intent(in) :: ndim, itmax
    logical, intent(in) :: warn

    this%func => func
    this%ndim = ndim
    this%ftol = ftol
    if (this%ftol == 0.0_num) this%ftol = real(1.d-6,num)
    this%itmax = itmax
    if (this%itmax == -1) this%itmax = 5000
    this%warn = warn
    this%iter = 0

    allocate(this%y(this%ndim+1), this%p(this%ndim+1,this%ndim))

    return
end subroutine create_NM

subroutine create_CG(this,func,gfunc,alpha,ndim,ftol,itmax,warn)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(inout) :: this
    procedure(vsfunc) :: func
    procedure(vvfunc) :: gfunc
    real(num), intent(in) :: alpha, ftol
    integer, intent(in) :: ndim, itmax
    logical, intent(in) :: warn

    this%func => func
    this%gfunc => gfunc
    this%ndim = ndim
    this%ftol = ftol
    if (this%ftol == 0.0_num) this%ftol = real(1.d-6,num)
    this%itmax = itmax
    if (this%itmax == -1) this%itmax = 200
    this%warn = warn
    this%iter = 0
    this%alpha = alpha

    allocate(this%xi(this%ndim), this%p(this%ndim))

    return
end subroutine create_CG

subroutine create_PS(this,func,ndim,npart,ftol,itmax,warn,w,c)
    implicit none
    class(OptimizeND_ParticleSwarm), intent(inout) :: this
    procedure(vsfunc) :: func
    real(num), intent(in) :: ftol
    integer, intent(in) :: ndim, npart, itmax
    logical, intent(in) :: warn
    real(num), intent(in), optional :: w, c(2)

    this%func => func
    this%ndim = ndim
    this%npart = npart
    this%ftol = ftol
    if (this%ftol == 0.0_num) this%ftol = real(1.d-6,num)
    this%itmax = itmax
    if (this%itmax == -1) this%itmax = 5000
    this%warn = warn
    this%iter = 0

    ! default values
    this%w = 0.8_num
    this%c = (/0.1_num, 0.1_num/)

    if (present(w)) then
    if ((w < 0.0_num).or.(abs(w) > 1.0_num)) then
        print '(A)', 'The inertia weight constant must be non-negative and less than 1.'
        print '(A)', 'Terminating execution.'
        call exit(1)
    else
        this%w = w
    end if
    end if

    if (present(c)) then
    if (any(c < 0.0_num).or.any(abs(c) > 1.0_num)) then
        print '(A)', 'The cognitive & social coefficients must be non-negative and less than 1.'
        print '(A)', 'Terminating execution.'
        call exit(1)
    else
        this%c = c
    end if
    end if

    allocate(this%x(this%npart,this%ndim), this%v(this%npart,this%ndim))
    allocate(this%xbest(this%npart,this%ndim), this%func_xbest(this%npart))
    allocate(this%p(this%ndim))

    return
end subroutine create_PS

subroutine create_GMC(this,func,mag_step,mag_delta,ndim,itmax,iroc,warn)
    implicit none
    class(OptimizeND_GuidedMonteCarlo), intent(inout) :: this
    procedure(vsfunc) :: func
    real(num), intent(in) :: mag_step(:), mag_delta(:)
    integer, intent(in) :: ndim, itmax, iroc
    logical, intent(in) :: warn

    this%func => func
    this%ndim = ndim
    this%itmax = itmax
    if (this%itmax == -1) this%itmax = int(1e+4)
    this%iroc = iroc
    if (this%iroc == -1) this%iroc = 3
    this%warn = warn
    this%iter = 0

    allocate(this%p(this%ndim), this%mag_delta(this%ndim), this%mag_step(this%ndim))

    ! magnitudes of step guides in ea. direction
    this%mag_step = mag_step

    ! magnitudes of small perturbations in ea. direction
    this%mag_delta = mag_delta

    ! Delgoda and Pulfer's guides (arbitrary units) based on frustration (0:20)
    this%step = real((/90.,90.,90.,90.,45.,22.5,11.25,5.625,2.8125,1.40625,&
        1.40625,1.40625,0.703125,0.3515625,0.17578125,0.087890625,0.0439453125,&
        0.02197265625,0.002197265625,0.0002197265625,0.00002197265625/),num)

    ! normalize (later to be multiplied by mag_step)
    this%step = this%step/maxval(this%step)

    return
end subroutine create_GMC

! destroy_ subroutines ---------------------------------------------------------

subroutine destroy_NM(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    this%func => null()
    if (allocated(this%y)) deallocate(this%y)
    if (allocated(this%p)) deallocate(this%p)
    return
end subroutine destroy_NM

subroutine destroy_CG(this)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(inout) :: this
    this%func => null()
    this%gfunc => null()
    if (allocated(this%xi)) deallocate(this%xi)
    if (allocated(this%p)) deallocate(this%p)
    return
end subroutine destroy_CG

subroutine destroy_PS(this)
    implicit none
    class(OptimizeND_ParticleSwarm), intent(inout) :: this
    this%func => null()
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%v)) deallocate(this%v)
    if (allocated(this%xbest)) deallocate(this%xbest)
    if (allocated(this%func_xbest)) deallocate(this%func_xbest)
    return
end subroutine destroy_PS

subroutine destroy_GMC(this)
    implicit none
    class(OptimizeND_GuidedMonteCarlo), intent(inout) :: this
    this%func => null()
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%mag_step)) deallocate(this%mag_step)
    if (allocated(this%mag_delta)) deallocate(this%mag_delta)
    return
end subroutine destroy_GMC

! reset_ subroutines -----------------------------------------------------------

subroutine reset_NM(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    this%iter = 0
    this%y = 0.0_num
    this%p = 0.0_num
    return
end subroutine reset_NM

subroutine reset_CG(this)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(inout) :: this
    this%iter = 0
    this%y = 0.0_num
    this%p = 0.0_num
    this%xi = 0.0_num
    return
end subroutine reset_CG

subroutine reset_PS(this)
    implicit none
    class(OptimizeND_ParticleSwarm), intent(inout) :: this
    this%iter = 0
    this%y = 0.0_num
    this%p = 0.0_num
    this%x = 0.0_num
    this%v = 0.0_num
    this%xbest = 0.0_num
    this%func_xbest = 0.0_num
    return
end subroutine reset_PS

subroutine reset_GMC(this)
    implicit none
    class(OptimizeND_GuidedMonteCarlo), intent(inout) :: this
    this%iter = 0
    this%y = 0.0_num
    this%p = 0.0_num
    return
end subroutine reset_GMC

! Nelder-Mead procedures -------------------------------------------------------

subroutine amoeba(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    real(num), parameter :: TINY=1.d-10

    integer :: i, n, ilo, ihi, inhi
    real(num) :: swap, rtol, ytry, ysave, psum(this%ndim)

1   do n=1,this%ndim
        psum(n) = sum(this%p(:,n))
    end do

2   ilo = 1

    if (this%y(1) > this%y(2)) then
        ihi = 1
        inhi = 2
    else
        ihi = 2
        inhi = 1
    end if

    do i=1,this%ndim+1
        if (this%y(i) <= this%y(ilo)) ilo = i
        if (this%y(i) > this%y(ihi)) then
            inhi = ihi
            ihi = i
        else if (this%y(i) > this%y(inhi)) then
            if (i /= ihi) inhi = i
        end if
    end do

    rtol = 2.*abs(this%y(ihi)-this%y(ilo))/(abs(this%y(ihi))+abs(this%y(ilo))+TINY)

    if (rtol < this%ftol) then
        swap = this%y(1)
        this%y(1) = this%y(ilo)
        this%y(ilo) = swap
        do n=1,this%ndim
            swap = this%p(1,n)
            this%p(1,n) = this%p(ilo,n)
            this%p(ilo,n) = swap
        end do
        return
    end if

    if (this%iter >= this%itmax) then
        if (this%warn) print '(A)', 'Warning: amoeba exceeded ITMAX'
        return
    end if

    this%iter = this%iter + 2
    ytry = this%amotry(psum,ihi,-1.0)

    if (ytry <= this%y(ilo)) then
        ytry = this%amotry(psum,ihi,2.0)
    else if (ytry >= this%y(inhi)) then
        ysave = this%y(ihi)
        ytry = this%amotry(psum,ihi,0.5)
        if (ytry >= ysave) then
            do i=1,this%ndim+1
                if (i /= ilo) then
                    psum = 0.5*(this%p(i,:) + this%p(ilo,:))
                    this%p(i,:) = psum
                    this%y(i) = this%func(psum)
                end if
            end do
            this%iter = this%iter + this%ndim
            goto 1
        end if
    else
        this%iter = this%iter - 1
    end if
    goto 2

    return
end subroutine amoeba

real(num) function amotry_NM(this,psum,ihi,fac)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    real(num), intent(inout) :: psum(this%ndim)
    integer, intent(in) :: ihi
    real, intent(in) :: fac

    real(num) :: fac1, fac2, ytry, ptry(this%ndim)

    fac1 = real(1.0_num-fac,num)/real(this%ndim,num)
    fac2 = fac1 - real(fac,num)

    ptry = fac1*psum - fac2*this%p(ihi,:)
    ytry = this%func(ptry)

    if (ytry < this%y(ihi)) then
        this%y(ihi) = ytry
        psum = psum - this%p(ihi,:) + ptry
        this%p(ihi,:) = ptry
    end if

    amotry_NM = ytry

    return
end function amotry_NM

! Conjugate-Gradient procedures ------------------------------------------------

subroutine frprmn(this)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(inout) :: this
    real(num), parameter :: EPS=1.d-10

    integer :: its
    real(num) :: dgg, fp, gam, gg
    real(num), dimension(this%ndim) :: g, h

    fp = this%func(this%p)
    this%xi = this%gfunc(this%p)
    g = -this%xi
    h = g
    this%xi = h

    do its=1,this%itmax
        this%iter = its
        call this%linmin()
        if (2.*abs(this%y-fp) <= this%ftol*(abs(this%y)+abs(fp)+EPS)) return
        fp = this%y
        this%xi = this%gfunc(this%p)
        gg = sum(g**2)
        dgg = sum((this%xi + g)*this%xi) ! Polak-Ribiere (improved)
        ! dgg = sum(this%xi**2) ! Fletcher-Reeves
        if (gg == 0.0_num) return
        gam = dgg/gg
        g = -this%xi
        h = g + gam*h
        this%xi = h
    end do

    if (this%warn) print '(A)', 'Warning: frprmn exceeded ITMAX'

    return
end subroutine frprmn

subroutine linmin_CG(this)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(inout) :: this

    real(num), parameter :: tol=1.d-4 ! dbrent tolerance
    real(num) :: ax, bx, fa, fb, fx, xmin, xx

    ! initial guess for brackets
    ax = 0.0_num
    xx = this%alpha

    call this%mnbrak(ax,xx,bx,fa,fx,fb)
    this%y = this%dbrent(ax,xx,bx,tol,xmin)

    this%xi = xmin*this%xi
    this%p = this%p + this%xi

    return
end subroutine linmin_CG

subroutine mnbrak_CG(this,ax,bx,cx,fa,fb,fc)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(in) :: this
    real(num), intent(inout) :: ax, bx, cx, fa, fb, fc
    real(num), parameter :: GOLD=(1.0_num+sqrt(5.0_num))/2.0_num, GLIMIT=100.0_num, TINY=1.d-20
    real(num) :: fu, q, r, u, ulim

    fa = this%func(this%p + ax*this%xi)
    fb = this%func(this%p + bx*this%xi)

    ! swap roles to search downhill
    if (fb > fa) then
        ax = ax + bx
        bx = ax - bx
        ax = ax - bx

        fa = fa + fb
        fb = fa - fb
        fa = fa - fb
    end if

    cx = bx + GOLD*(bx-ax)
    fc = this%func(this%p + cx*this%xi)

1   if (fb >= fc) then
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx - ((bx-cx)*q - (bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim = bx + GLIMIT*(cx-bx)
        if ((bx-u)*(u-cx) > 0.) then
            fu = this%func(this%p + u*this%xi)

            if (fu < fc) then
                ax = bx
                fa = fb
                bx = u
                fb = fu
                return
            else if (fu > fb) then
                cx = u
                fc = fu
                return
            end if
            u = cx + GOLD*(cx-bx)
            fu = this%func(this%p + u*this%xi)
        else if ((cx-u)*(u-ulim) > 0.) then
            fu = this%func(this%p + u*this%xi)
            if (fu < fc) then
                bx = cx
                cx = u
                u = cx + GOLD*(cx-bx)
                fb = fc
                fc = fu
                fu = this%func(this%p + u*this%xi)
            end if
        else if ((u-ulim)*(ulim-cx) >= 0.) then
            u = ulim
            fu = this%func(this%p + u*this%xi)
        else
            u = cx + GOLD*(cx-bx)
            fu = this%func(this%p + u*this%xi)
        end if
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        goto 1
    end if

    return
end subroutine mnbrak_CG

real(num) function dbrent_CG(this,ax,bx,cx,ftol,xmin)
    implicit none
    class(OptimizeND_ConjugateGradient), intent(in) :: this
    real(num), intent(in) :: ax, bx, cx, ftol
    real(num), intent(out) :: xmin
    real(num), parameter :: ZEPS=1.d-10

    integer :: iter
    logical :: ok1,ok2
    real(num) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,&
        u,u1,u2,v,w,x,xm

    a = min(ax,cx)
    b = max(ax,cx)
    v = bx
    w = v
    x = v
    e = 0.0_num
    fx = this%func(this%p + x*this%xi)
    fv = fx
    fw = fx
    dx = sum(this%xi*this%gfunc(this%p + x*this%xi))
    dv = dx
    dw = dx

    do iter=1,this%itmax
        xm = 0.5*(a+b)
        tol1 = ftol*abs(x) + ZEPS
        tol2 = 2.*tol1
        if (abs(x-xm) <= (tol2-0.5*(b-a))) goto 3
        if (abs(e) > tol1) then
            d1 = 2.*(b-a)
            d2 = d1
            if (dw /= dx) d1 = (w-x)*dx/(dx-dw)
            if (dv /= dx) d2 = (v-x)*dx/(dx-dv)
            u1 = x+d1
            u2 = x+d2
            ok1 = ((a-u1)*(u1-b)>0.).and.(dx*d1<=0.)
            ok2 = ((a-u2)*(u2-b)>0.).and.(dx*d2<=0.)
            olde = e
            e = d
            if (.not.(ok1.or.ok2)) then
                goto 1
            else if (ok1.and.ok2) then
                if (abs(d1) < abs(d2)) then
                    d = d1
                else
                    d = d2
                end if
            else if (ok1) then
                d = d1
            else
                d = d2
            end if
            if (abs(d) > abs(0.5*olde)) goto 1
            u = x+d
            if (((u-a) < tol2).or.((b-u) < tol2)) d = sign(tol1,xm-x)
            goto 2
        end if
1       if (dx >= 0.) then
            e = a-x
        else
            e = b-x
        end if
        d = 0.5*e
2       if (abs(d) >= tol1) then
            u = x+d
            fu = this%func(this%p + u*this%xi)
        else
            u = x + sign(tol1,d)
            fu = this%func(this%p + u*this%xi)
            if (fu > fx) goto 3
        end if
        du = sum(this%xi*this%gfunc(this%p + u*this%xi))
        if (fu <= fx) then
            if (u >= x) then
                a = x
            else
                b = x
            end if
            v = w
            fv = fw
            dv = dw
            w = x
            fw = fx
            dw = dx
            x = u
            fx = fu
            dx = du
        else
            if (u < x) then
                a = u
            else
                b = u
            end if
            if ((fu <= fw).or.(w == x)) then
                v = w
                fv = fw
                dv = dw
                w = u
                fw = fu
                dw = du
            else if ((fu <= fv).or.(v == x).or.(v == w)) then
                v = u
                fv = fu
                dv = du
            end if
        end if
    end do

    if (this%warn) print '(A)', 'Warning: dbrent_CG exceeded ITMAX'
3   xmin = x
    dbrent_CG = fx

    return
end function dbrent_CG

! Particle Swarm procedures ----------------------------------------------------

subroutine nemo(this)
    implicit none
    class(OptimizeND_ParticleSwarm), intent(inout) :: this
    real(num) :: r(2), temp
    integer :: i, n

    this%xbest = this%x

    do n=1,this%npart
        this%func_xbest(n) = this%func(this%xbest(n,:))
    end do

    this%p = this%xbest(minloc(this%func_xbest,1),:)
    this%y = minval(this%func_xbest)

    do i=1,this%itmax

        call random_number(r) ! direction vector

        ! parallel-advance all particles
        !$OMP PARALLEL DO SHARED(this) PRIVATE(temp)
        do n=1,this%npart
            ! update velocity
            this%v(n,:) = this%w*this%v(n,:) &
                + this%c(1)*r(1)*(this%xbest(n,:) - this%x(n,:)) &
                + this%c(2)*r(2)*(this%p - this%x(n,:))
            ! update position
            this%x(n,:) = this%x(n,:) + this%v(n,:)

            ! re-evaluate objective function
            temp = this%func(this%x(n,:))
            ! update master list if necessary
            if (temp < this%func_xbest(n)) then
                this%xbest(n,:) = this%x(n,:)
                this%func_xbest(n) = temp
            end if
        end do
        !$OMP END PARALLEL DO

        ! update the swarm's best point
        this%p = this%xbest(minloc(this%func_xbest,1),:)
        this%y = minval(this%func_xbest)

        this%iter = this%iter + 1
        if (this%span() < this%ftol) return

    end do

    if (this%warn) print '(A)', 'Warning: nemo exceeded ITMAX'

    return
end subroutine nemo

real(num) function flockspan_PS(this)
    implicit none
    class(OptimizeND_ParticleSwarm), intent(in) :: this
    real(num) :: dist
    integer :: i, j

    flockspan_PS = norm2(this%x(1,:)-this%x(2,:))

    do i=1,this%npart
        do j=i+1,this%npart
            dist = norm2(this%x(i,:)-this%x(j,:))
            if (dist < flockspan_PS) flockspan_PS = dist
        end do
    end do

    return
end function flockspan_PS

! Guided Monte Carlo procedures ------------------------------------------------

subroutine gmcmn(this)
    implicit none
    class(OptimizeND_GuidedMonteCarlo), intent(inout) :: this
    real(num), dimension(this%ndim) :: b, delta, rx, cx
    real(num) :: fb, fx, rxi, cxi
    integer :: fcount, istep, iroc
    fcount=0; istep=0; iroc=0;

    ! call init_RNG()
    call random_number(delta)
    b = this%p + this%mag_delta*delta
    fb = huge(0.0_num)

2   fx = this%func(this%p)
    this%iter = this%iter + 1

    if ((this%itmax /= -2).and.(this%iter >= this%itmax)) then
        ! set the solution as the best point and terminate
        this%p = b
        this%y = this%func(this%p)
        if (this%warn) &
            print '(A)', 'Warning: Optimization terminated prematurely (ITMAX reached)'
        return
    end if

    ! frustration check
    if (fx <= fb) then
        b = this%p
        fb = fx
        fcount = 0
        goto 5
    else
        fcount = fcount + 1
    end if

    if ((fcount>20).and.(iroc>this%iroc)) then
        ! normal termination; minimum identified
        this%p = b; this%y = this%func(this%p);
        return
    else if (fcount>20) then
        istep = istep + 1
        fcount = 0
    end if

    if (istep>this%ndim) then
        iroc = iroc + 1
        istep = 0
        fcount = 0
    end if

5   if (istep == 0) then
        ! search all variables simultaneously
        call random_number(rx); rx = rx*this%step(fcount)*this%mag_step; ! step-sizes
        call random_number(cx); cx = 2.*cx - 1.; ! direction
        this%p = b + rx*cx
    else ! istep > 0
        ! search only p(istep), freezing all others in their best
        call random_number(rxi); rxi = rxi*this%step(fcount)*this%mag_step(istep); ! step-size
        call random_number(cxi); cxi = 2.*cxi - 1.; ! direction
        this%p(istep) = b(istep) + rxi*cxi
    end if

    goto 2

    return
end subroutine gmcmn

! 1D optimization procedures ---------------------------------------------------

subroutine golden(func,ax,bx,cx,ftol,xmin,fmin)
    implicit none
    procedure(ssfunc) :: func
    real(num), intent(in) :: ax, bx, cx, ftol
    real(num), intent(out) :: xmin, fmin
    real(num), parameter :: R=(sqrt(5.0_num)-1.0_num)/2.0_num, C=1.0_num-R
    real(num) :: f1, f2, x0, x1, x2, x3

    x0 = ax
    x3 = cx

    if (abs(cx-bx) > abs(bx-ax)) then
        x1 = bx
        x2 = bx + C*(cx-bx)
    else
        x2 = bx
        x1 = bx - C*(bx-ax)
    end if

    f1 = func(x1)
    f2 = func(x2)

1   if (abs(x3-x0) > ftol*(abs(x1)+abs(x2))) then
        if (f2 < f1) then
            x0 = x1
            x1 = x2
            x2 = R*x1 + C*x3
            f1 = f2
            f2 = func(x2)
        else
            x3 = x2
            x2 = x1
            x1 = R*x2 + C*x0
            f2 = f1
            f1 = func(x1)
        end if
        goto 1
    end if

    if (f1 < f2) then
        fmin = f1
        xmin = x1
    else
        fmin = f2
        xmin = x2
    end if

    return
end subroutine golden

subroutine brent(func,ax,bx,cx,ftol,xmin,fmin)
    implicit none
    procedure(ssfunc) :: func
    real(num), intent(in) :: ax, bx, cx, ftol
    real(num), intent(out) :: xmin, fmin
    integer, parameter :: ITMAX=100
    real(num), parameter :: CGOLD=(3.0_num-sqrt(5.0_num))/2.0_num, ZEPS=1.d-10

    real(num) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    integer :: iter

    a = min(ax,cx)
    b = max(ax,cx)
    v = bx
    w = v
    x = v
    e = 0.0_num
    fx = func(x)
    fv = fx
    fw = fx

    do iter=1,ITMAX
        xm = 0.5*(a+b)
        tol1 = ftol*abs(x) + ZEPS
        tol2 = 2.*tol1
        if (abs(x-xm) <= (tol2-0.5*(b-a))) goto 3
        if (abs(e) > tol1) then
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q - (x-w)*r
            q = 2.*(q-r)
            if (q > 0.0_num) p = -p
            q = abs(q)
            etemp = e
            e = d
            if ((abs(p) >= abs(0.5*q*etemp)).or.(p <= q*(a-x)).or.(p >= q*(b-x))) goto 1
            d = p/q
            u = x + d
            if ((u-a < tol2).or.(b-u < tol2)) d = sign(tol1,xm-x)
            goto 2
        end if
1       if (x >= xm) then
            e = a-x
        else
            e = b-x
        end if
        d = CGOLD*e
2       if (abs(d) >= tol1) then
            u = x + d
        else
            u = x + sign(tol1,d)
        end if
        fu = func(u)
        if (fu <= fx) then
            if (u >= x) then
                a = x
            else
                b = x
            end if
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
        else
            if (u < x) then
                a = u
            else
                b = u
            end if
            if ((fu <= fw).or.(w == x)) then
                v = w
                fv = fw
                w = u
                fw = fu
            else if ((fu <= fv).or.(v == x).or.(v == w)) then
                v = u
                fv = fu
            end if
        end if
    end do

    print '(A)', 'Warning: brent exceeded ITMAX'
3   xmin = x
    fmin = fx

    return
end subroutine brent

real(num) function dbrent(func,dfunc,ax,bx,cx,ftol,xmin)
    implicit none
    procedure(ssfunc) :: func
    procedure(ssfunc) :: dfunc
    real(num), intent(in) :: ax, bx, cx, ftol
    real(num), intent(out) :: xmin
    integer, parameter :: ITMAX=100
    real(num), parameter :: ZEPS=1.d-10

    integer :: iter
    logical :: ok1,ok2
    real(num) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,&
        u,u1,u2,v,w,x,xm

    a = min(ax,cx)
    b = max(ax,cx)
    v = bx
    w = v
    x = v
    e = 0.0_num
    fx = func(x)
    fv = fx
    fw = fx
    dx = dfunc(x)
    dv = dx
    dw = dx

    do iter=1,ITMAX
        xm = 0.5*(a+b)
        tol1 = ftol*abs(x) + ZEPS
        tol2 = 2.*tol1
        if (abs(x-xm) <= (tol2-0.5*(b-a))) goto 3
        if (abs(e) > tol1) then
            d1 = 2.*(b-a)
            d2 = d1
            if (dw /= dx) d1 = (w-x)*dx/(dx-dw)
            if (dv /= dx) d2 = (v-x)*dx/(dx-dv)
            u1 = x+d1
            u2 = x+d2
            ok1 = ((a-u1)*(u1-b)>0.).and.(dx*d1<=0.)
            ok2 = ((a-u2)*(u2-b)>0.).and.(dx*d2<=0.)
            olde = e
            e = d
            if (.not.(ok1.or.ok2)) then
                goto 1
            else if (ok1.and.ok2) then
                if (abs(d1) < abs(d2)) then
                    d = d1
                else
                    d = d2
                end if
            else if (ok1) then
                d = d1
            else
                d = d2
            end if
            if (abs(d) > abs(0.5*olde)) goto 1
            u = x+d
            if (((u-a) < tol2).or.((b-u) < tol2)) d = sign(tol1,xm-x)
            goto 2
        end if
1       if (dx >= 0.) then
            e = a-x
        else
            e = b-x
        end if
        d = 0.5*e
2       if (abs(d) >= tol1) then
            u = x+d
            fu = func(u)
        else
            u = x + sign(tol1,d)
            fu = func(u)
            if (fu > fx) goto 3
        end if
        du = dfunc(u)
        if (fu <= fx) then
            if (u >= x) then
                a = x
            else
                b = x
            end if
            v = w
            fv = fw
            dv = dw
            w = x
            fw = fx
            dw = dx
            x = u
            fx = fu
            dx = du
        else
            if (u < x) then
                a = u
            else
                b = u
            end if
            if ((fu <= fw).or.(w == x)) then
                v = w
                fv = fw
                dv = dw
                w = u
                fw = fu
                dw = du
            else if ((fu <= fv).or.(v == x).or.(v == w)) then
                v = u
                fv = fu
                dv = du
            end if
        end if
    end do

    print '(A)', 'Warning: dbrent exceeded ITMAX'
3   xmin = x
    dbrent = fx

    return
end function dbrent

subroutine mnbrak(func,ax,bx,cx,fa,fb,fc)
    implicit none
    procedure(ssfunc) :: func
    real(num), intent(inout) :: ax, bx, cx, fa, fb, fc
    real(num), parameter :: GOLD=(1.0_num+sqrt(5.0_num))/2.0_num, GLIMIT=100.0_num, TINY=1.d-20
    real(num) :: fu, q, r, u, ulim

    fa = func(ax)
    fb = func(bx)

    ! swap roles to search downhill
    if (fb > fa) then
        ax = ax + bx
        bx = ax - bx
        ax = ax - bx

        fa = fa + fb
        fb = fa - fb
        fa = fa - fb
    end if

    cx = bx + GOLD*(bx-ax)
    fc = func(cx)

1   if (fb >= fc) then
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx - ((bx-cx)*q - (bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim = bx + GLIMIT*(cx-bx)
        if ((bx-u)*(u-cx) > 0.) then
            fu = func(u)

            if (fu < fc) then
                ax = bx
                fa = fb
                bx = u
                fb = fu
                return
            else if (fu > fb) then
                cx = u
                fc = fu
                return
            end if
            u = cx + GOLD*(cx-bx)
            fu = func(u)
        else if ((cx-u)*(u-ulim) > 0.) then
            fu = func(u)
            if (fu < fc) then
                bx = cx
                cx = u
                u = cx + GOLD*(cx-bx)
                fb = fc
                fc = fu
                fu = func(u)
            end if
        else if ((u-ulim)*(ulim-cx) >= 0.) then
            u = ulim
            fu = func(u)
        else
            u = cx + GOLD*(cx-bx)
            fu = func(u)
        end if
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        goto 1
    end if

    return
end subroutine mnbrak

! Miscellaneous procedures -----------------------------------------------------

real(num) function dfridr_ss(func,x,h,err)
    implicit none
    procedure(ssfunc) :: func
    real(num), intent(in) :: x, h
    real(num), intent(out) :: err

    integer, parameter :: NTAB=10
    real(num), parameter :: CON=1.4_num, CON2=CON*CON, BIG=1.d+30, SAFE=2.0_num

    integer :: j, k
    real(num) :: errt, fac, hh, a(NTAB,NTAB)

    if (h == 0.0_num) then
        dfridr_ss = 0.0_num
        print '(A)', 'Error: step-size (h) must be non-zero in dfridr'
        return
    end if

    hh = h
    a(1,1) = (func(x+hh)-func(x-hh))/(2.*hh)
    err = BIG

    do k=2,NTAB
        hh = hh/CON
        a(1,k) = (func(x+hh)-func(x-hh))/(2.*hh)
        fac = CON2
        do j=2,k
            a(j,k) = (fac*a(j-1,k)-a(j-1,k-1))/(fac-1.)
            fac = CON2*fac
            errt = max(abs(a(j,k)-a(j-1,k)),abs(a(j,k)-a(j-1,k-1)))
            if (errt <= err) then
                err = errt
                dfridr_ss = a(j,k)
            end if
        end do
        if (abs(a(k,k)-a(k-1,k-1)) >= SAFE*err) return
    end do

    return
end function dfridr_ss

real(num) function dfridr_vs(func,n,x,h,err)
    implicit none
    procedure(vsfunc) :: func
    integer, intent(in) :: n
    real(num), intent(in) :: x(:), h
    real(num), intent(out) :: err

    integer, parameter :: NTAB=10
    real(num), parameter :: CON=1.4_num, CON2=CON*CON, BIG=1.d+30, SAFE=2.0_num

    integer :: j, k
    real(num) :: ei(size(x)), errt, fac, hh, a(NTAB,NTAB)

    if (h == 0.0_num) then
        dfridr_vs = 0.0_num
        print '(A)', 'Error: step-size (h) must be non-zero in dfridr'
        return
    end if

    ei = 0.
    ei(n) = 1.

    hh = h
    a(1,1) = (func(x+hh*ei)-func(x-hh*ei))/(2.*hh)
    err = BIG

    do k=2,NTAB
        hh = hh/CON
        a(1,k) = (func(x+hh*ei)-func(x-hh*ei))/(2.*hh)
        fac = CON2
        do j=2,k
            a(j,k) = (fac*a(j-1,k)-a(j-1,k-1))/(fac-1.)
            fac = CON2*fac
            errt = max(abs(a(j,k)-a(j-1,k)),abs(a(j,k)-a(j-1,k-1)))
            if (errt <= err) then
                err = errt
                dfridr_vs = a(j,k)
            end if
        end do
        if (abs(a(k,k)-a(k-1,k-1)) >= SAFE*err) return
    end do

    return
end function dfridr_vs

end module optimize
