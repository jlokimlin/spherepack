!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK version 3.2                   *
!     *                                                               *
!     *       A Package of Fortran Subroutines and Programs           *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!
!     This version of compute_gaussian_latitudes_and_weights implements
!     the method presented in:
!     P. N. Swarztrauber, Computing the points and weights for
!     Gauss-Legendre quadrature, SIAM J. Sci. Comput.,
!     24(2002) pp. 945-954.
!
!     The w and lwork arrays are dummy and included only to
!     permit a simple pluggable exchange with the
!     old compute_gaussian_latitudes_and_weights in previous versions of SPHEREPACK
!
!
!     gauss points and weights are computed using the fourier-newton
!     described in "on computing the points and weights for 
!     gauss-legendre quadrature", paul n. swarztrauber, siam journal 
!     on scientific computing that has been accepted for publication.
!     This routine is faster and more accurate than older program
!     with the same name.
!
!     subroutine compute_gaussian_latitudes_and_weights computes the nlat gaussian colatitudes and weights
!     in real. the colatitudes are in radians and lie in the
!     in the interval (0, pi).
!
!     input parameters
!
!     nlat    the number of gaussian colatitudes in the interval (0, pi)
!             (between the two poles).  nlat must be greater than zero.
!
!     w       unused real variable that permits a simple
!             exchange with the old routine with the same name
!             in spherepack.
!
!     lwork   unused variable that permits a simple exchange with the
!             old routine with the same name in spherepack.
!
!     output parameters
!
!     theta   a real array with length nlat
!             containing the gaussian colatitudes in
!             increasing radians on the interval (0, pi).
!
!     wts     a real array with lenght nlat
!             containing the gaussian weights.
!
!     ierror = 0 no errors
!            = 1 if nlat <= 0
!
module gaussian_latitudes_and_weights_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: compute_gaussian_latitudes_and_weights

    !------------------------------------------------------------------
    ! Dictionary: Variables confined to the module
    !------------------------------------------------------------------
    real (wp), parameter :: ZERO = 0.0_wp
    real (wp), parameter :: HALF = 0.5_wp
    real (wp), parameter :: ONE = 1.0_wp
    real (wp), parameter :: TWO = 2.0_wp
    real (wp), parameter :: THREE = 3.0_wp
    !------------------------------------------------------------------

contains

    subroutine compute_gaussian_latitudes_and_weights(nlat, theta, wts, w, lwork, ierror)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        real (wp),    intent (out) :: theta(nlat)
        real (wp),    intent (out) :: wts(nlat)
        real (wp),    intent (in)  :: w
        integer (ip), intent (in)  :: lwork
        integer (ip), intent (out) :: ierror
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)         :: i, idx, it, nix
        integer (ip)         :: nhalf, half_nlat
        real (wp), parameter :: HALF_PI = PI/2
        real (wp), parameter :: SQRT3 = sqrt(THREE)
        real (wp)            :: dt, half_dt
        real (wp)            :: zprev, zlast, eff_zero
        real (wp)            :: zhold, pb, dpb, dcor, cz
        real (wp)            :: eps, sgnd
        !----------------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (nlat <= 0) then
            ierror = 1
            return
        else
            ierror = 0
        end if

        !
        !==> compute weights and points analytically when nlat=1, 2
        !
        select case (nlat)
            case(1)
                theta = HALF_PI
                wts = TWO
            case(2)
                theta(1) = acos(SQRT3/3)
                theta(2) = acos(-SQRT3/3)
                wts = ONE
            case default
                eps = sqrt(epsilon(ONE))
                eps = eps * sqrt(eps)
                half_nlat = nlat/2
                nhalf = (nlat+1)/2
                idx = half_nlat+2

                call compute_fourier_coefficients(nlat, cz, theta(half_nlat+1), wts(half_nlat+1))

                dt = HALF_PI/nhalf
                half_dt = dt/2
                !
                !==> Estimate first point next to theta = pi/2
                !
                select case (mod(nlat,2))
                    case (0)
                        !
                        !==> nlat even
                        !
                        eff_zero = HALF_PI-half_dt
                        nix = nhalf
                    case default
                        !
                        !==> nlat odd
                        !
                        eff_zero = HALF_PI-dt
                        zprev = HALF_PI
                        nix = nhalf-1
                end select

                start_iteration: do
                    it = 0
                    newton_iteration: do
                        it = it+1
                        zlast = eff_zero
                        !
                        !==> Newton iterations
                        !
                        call compute_legendre_poly_and_deriv(nlat, eff_zero, cz, theta(half_nlat+1), wts(half_nlat+1), pb, dpb)

                        dcor = pb/dpb
                        sgnd = ONE

                        if (dcor /= ZERO) sgnd = dcor/abs(dcor)

                        dcor = sgnd * min(abs(dcor), 0.2_wp * dt)
                        eff_zero = eff_zero-dcor

                        !
                        !==> Repeat iteration
                        !
                        if (abs(eff_zero-zlast) - eps*abs(eff_zero) > ZERO) cycle newton_iteration

                        theta(nix) = eff_zero
                        zhold = eff_zero
                        !
                        ! ==> yakimiw's formula permits using old pb and dpb
                        !
                        wts(nix) = (2*nlat+1)/(dpb+pb*cos(zlast)/sin(zlast))**2
                        nix = nix-1

                        if (nix /= 0) then

                            if (nix == nhalf-1)  then
                                eff_zero = THREE * eff_zero - PI
                            else if (nix < nhalf-1)  then
                                eff_zero = TWO * eff_zero-zprev
                            end if

                            zprev = zhold

                            !
                            !==> Re-initialize loop
                            !
                            cycle start_iteration
                        end if
                        exit newton_iteration
                    end do newton_iteration
                    exit start_iteration
                end do start_iteration
                !
                !==> Extend points and weights via symmetries
                !
                if (mod(nlat,2) /= 0) then

                    theta(nhalf) = HALF_PI

                    call compute_legendre_poly_and_deriv(nlat, HALF_PI, cz, theta(half_nlat+1), wts(half_nlat+1), pb, dpb)

                    wts(nhalf) = real(2*nlat+1, kind=wp)/(dpb**2)

                end if

                do i=1, half_nlat
                    wts(nlat-i+1) = wts(i)
                    theta(nlat-i+1) = PI-theta(i)
                end do
                !
                !==> Set weights
                !
                wts = TWO * wts/sum(wts)
        end select


    end subroutine compute_gaussian_latitudes_and_weights


    pure subroutine compute_fourier_coefficients( &
        n, cz, legendre_poly_coeff, legendre_deriv_coeff)
        !
        ! Purpose:
        !
        ! Computes the fourier coefficients of the legendre
        ! polynomial p_n^0 and its derivative.
        ! n is the degree and n/2 or (n+1)/2
        ! coefficients are returned in cp depending on whether
        ! n is even or odd. The same number of coefficients
        ! are returned in dcp. For n even the constant
        ! coefficient is returned in cz.
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: cz
        real (wp),    intent (out) :: legendre_poly_coeff(n/2+1)
        real (wp),    intent (out) :: legendre_deriv_coeff(n/2+1)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: j
        real (wp)    :: t1, t2, t3, t4
        !----------------------------------------------------------------------

        associate( &
            ncp => (n+1)/2, &
            cp => legendre_poly_coeff, &
            dcp => legendre_deriv_coeff &
            )

            t1 = -ONE
            t2 = real(n + 1, kind=wp)
            t3 = ZERO
            t4 = real(2*n + 1, kind=wp)

            select case (mod(n,2))
                case (0)
                    !
                    !==> n even
                    !
                    cp(ncp) = ONE

                    do j = ncp, 2, -1
                        t1 = t1+TWO
                        t2 = t2-ONE
                        t3 = t3+ONE
                        t4 = t4-TWO
                        cp(j-1) = (t1*t2)/(t3*t4)*cp(j)
                    end do

                    t1 = t1+TWO
                    t2 = t2-ONE
                    t3 = t3+ONE
                    t4 = t4-TWO
                    cz = (t1*t2)/(t3*t4)*cp(1)

                    do j=1, ncp
                        dcp(j) = real(2*j, kind=wp)*cp(j)
                    end do

                case default
                    !
                    !==> odd
                    !
                    cp(ncp) = ONE
                    do j = ncp-1, 1, -1
                        t1 = t1+TWO
                        t2 = t2-ONE
                        t3 = t3+ONE
                        t4 = t4-TWO
                        cp(j) = (t1*t2)/(t3*t4)*cp(j+1)
                    end do

                    do j=1, ncp
                        dcp(j) = real(2*j-1, kind=wp)*cp(j)
                    end do

            end select
        end associate

    end subroutine compute_fourier_coefficients


    pure subroutine compute_legendre_poly_and_deriv( &
        n, theta, cz, legendre_poly_coeff, legendre_deriv_coeff, &
        legendre_poly, legendre_deriv)
        !
        ! Purpose:
        !
        ! Computes pn(theta) and its derivative dpb(theta) with
        ! respect to theta
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (in)  :: theta
        real (wp),    intent (in)  :: cz
        real (wp),    intent (in)  :: legendre_poly_coeff(n/2+1)
        real (wp),    intent (in)  :: legendre_deriv_coeff(n/2+1)
        real (wp),    intent (out) :: legendre_poly
        real (wp),    intent (out) :: legendre_deriv
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: k, kdo
        real (wp)    :: cost, sint, temp
        !----------------------------------------------------------------------

        associate( &
            cos2t => cos(TWO * theta), &
            sin2t => sin(TWO * theta), &
            cp => legendre_poly_coeff, &
            dcp => legendre_deriv_coeff, &
            pb => legendre_poly, &
            dpb => legendre_deriv &
            )

            select case (mod(n,2))
                case (0)
                    !
                    !==> n even
                    !
                    kdo = n/2
                    pb = HALF * cz
                    dpb = ZERO
                    if (n > 0) then
                        cost = cos2t
                        sint = sin2t
                        do k=1, kdo
                            pb = pb+cp(k)*cost
                            dpb = dpb-dcp(k)*sint
                            temp = cos2t*cost-sin2t*sint
                            sint = sin2t*cost+cos2t*sint
                            cost = temp
                        end do
                    end if
                case default
                    !
                    !==> n odd
                    !
                    kdo = (n + 1)/2
                    pb = ZERO
                    dpb = ZERO
                    cost = cos(theta)
                    sint = sin(theta)
                    do k=1, kdo
                        pb = pb+cp(k)*cost
                        dpb = dpb-dcp(k)*sint
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
            end select
        end associate

    end subroutine compute_legendre_poly_and_deriv

end module gaussian_latitudes_and_weights_routines
