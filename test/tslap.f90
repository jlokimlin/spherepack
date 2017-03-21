!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                          Spherepack                           *
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
!     3/97
!
!     a program for testing slap, islap (ec, es, gc, gs)
!
!     (1) set a scalar field s as poly in x, y, z restricted to sphere
!
!     (2) compute scalar laplacian in array sclp using slap(ec, es, gc, gs)
!
!     (3) compare (2) with analytic scalar laplacian in sclpe
!
!     (4) compute the inverse  of (2) and compare with (1)
!
program tslap

    use spherepack

    ! Explicit typing only
    implicit none

    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: d2xdp2
    real(wp) :: d2xdt2
    real(wp) :: d2ydp2
    real(wp) :: d2ydt2
    real(wp) :: d2zdp2
    real(wp) :: d2zdt2
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: dxdp
    real(wp) :: dxdt
    real(wp) :: dydp
    real(wp) :: dydt
    real(wp) :: dzdp
    real(wp) :: dzdt
    real(wp) :: err2
    real(wp) :: err2s
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: error_flag
    integer(ip) :: j
    integer(ip) :: k
    real(wp) :: phi
    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    integer(ip), parameter :: nlat=15, nlon= 22, nt = 3
    integer(ip), parameter :: mdab = (nlon+2)/2
    integer(ip), parameter :: isym = 0
    real(wp), allocatable  :: wavetable(:)
    real(wp), dimension(mdab, nlat, nt) :: a, b
    real(wp), dimension(nlat, nlon, nt) :: s, se, sclp, sclpe
    real(wp), dimension(nt)             :: perturbation, xlm
    real(wp), dimension(nlat)           :: gaussian_latitudes, gaussian_weights
    real(wp), parameter                 :: ZERO = 0.0_wp, TWO = 2.0_wp

    ! Print dimension variables
    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")

    ! Set equally spaced colatitude and longitude increments
    dphi = TWO_PI/nlon
    dlat = PI/(nlat-1)

    ! Compute nlat-many gaussian points
    call compute_gaussian_latitudes_and_weights(nlat, &
        gaussian_latitudes, gaussian_weights, error_flag)

    call name("compute_gaussian_latitudes_and_weights")
    call iout(error_flag, " error_flag = ")
    call vecout(gaussian_latitudes, "gaussian_latitudes", nlat)

    ! Set helmholtz constant zero for laplacian inversion
    xlm = ZERO

    ! Test all laplacian and inversion routines
    do icase=1, 4
        call name("****")
        call name("****")
        call iout(icase, "icas")

        ! Set scalar field as poly in x, y, z restricted to the sphere
        do k=1, nt
            do j=1, nlon
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    select case (icase)
                        case(0:2)
                            theta = real(i - 1, kind=wp) * dlat
                        case default
                            theta = gaussian_latitudes(i)
                    end select
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    d2xdt2 = -x
                    dxdp = -cost*sinp
                    d2xdp2 = -x
                    dydt = cost*sinp
                    d2ydt2 = -y
                    dydp = cost*cosp
                    d2ydp2 = -y
                    dzdt = -sint
                    d2zdt2 = -z
                    dzdp = ZERO
                    d2zdp2 = ZERO
                    select case (k)
                        case (1)
                            s(i, j, k) = x+y
                            sclpe(i, j, k) = -TWO * (x+y)
                        case (2)
                            s(i, j, k) = x+z
                            sclpe(i, j, k) = -TWO * (x+z)
                        case (3)
                            s(i, j, k) = y+z
                            sclpe(i, j, k) = -TWO * (y+z)
                    end select
                end do
            end do
        end do

        ! Store exact solution
        se = s

        select case (icase)
            case (1)
        		
                call name("**ec")
        		
                call initialize_shaec(nlat, nlon, wavetable, error_flag)
                call name("shai")
                call iout(error_flag, "error_flag = ")
        		
                call shaec(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wavetable, error_flag)
                call name("sha ")
                call iout(error_flag, "error_flag = ")
        		
                call initialize_shsec(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call slapec(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, error_flag)
                call name("slap")
                call iout(error_flag, "error_flag = ")
            case (2)
        		
                call name("**es")
        		
                call initialize_shaes(nlat, nlon, wavetable, error_flag)
                call name("shai")
                call iout(error_flag, "error_flag = ")
        		
                call shaes(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wavetable, error_flag)
                call name("sha ")
                call iout(error_flag, "error_flag = ")
        		
                call initialize_shses(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call slapes(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, error_flag)
                call name("slap")
                call iout(error_flag, "error_flag = ")
            case (3)
        		
                call name("**gc")
        		
                call initialize_shagc(nlat, nlon, wavetable, error_flag)
                call name("shai")
                call iout(error_flag, "error_flag = ")
        		
                call shagc(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wavetable, error_flag)
                call name("sha ")
                call iout(error_flag, "error_flag = ")
        		
                call initialize_shsgc(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call slapgc(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, error_flag)
                call name("slap")
                call iout(error_flag, "error_flag = ")
            case (4)
        		
                call name("**gs")
        		
                call initialize_shags(nlat, nlon, wavetable, error_flag)
                call name("shai")
                call iout(error_flag, "error_flag = ")
        		
                call shags(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wavetable, error_flag)
                call name("sha ")
                call iout(error_flag, "error_flag = ")
        		
                call initialize_shsgs(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call slapgs(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, error_flag)
                call name("slap")
                call iout(error_flag, "error_flag = ")
        end select

        ! Compute discretization error" in sclp
        err2 = norm2(sclpe - sclp)
        call vout(err2, "err2")

        ! invert sclp
        s = ZERO
        select case (icase)
            case (1)
        		
                call initialize_shaec(nlat, nlon, wavetable, error_flag)
                call shaec(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, error_flag)
        		
                call initialize_shsec(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call islapec(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, perturbation, error_flag)
                call name("isla")
                call iout(error_flag, "error_flag = ")
                call vecout(perturbation, "perturbation = ", nt)
            case (2)
        		
                call initialize_shaes(nlat, nlon, wavetable, error_flag)
                call shaes(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, error_flag)
        		
                call initialize_shses(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call islapes(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, perturbation, error_flag)
                call name("isla")
                call iout(error_flag, "error_flag = ")
                call vecout(perturbation, "perturbation = ", nt)
            case (3)
        		
                call initialize_shagc(nlat, nlon, wavetable, error_flag)
                call shagc(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, wavetable, error_flag)
        		
                call initialize_shsgc(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call islapgc(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, perturbation, error_flag)
                call name("isla")
                call iout(error_flag, "error_flag = ")
                call vecout(perturbation, "perturbation = ", nt)
            case (4)
        		
                call name("**gs")
        		
                call initialize_shags(nlat, nlon, wavetable, error_flag)
                call shags(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, error_flag)
        		
                call initialize_shsgs(nlat, nlon, wavetable, error_flag)
                call name("shsi")
                call iout(error_flag, "error_flag = ")
        		
                call islapgs(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
                    wavetable, perturbation, error_flag)
                call name("isla")
                call iout(error_flag, "error_flag = ")
                call vecout(perturbation, "perturbation = ", nt)
        end select

        ! Compare s with original
        err2s = norm2(se - s)
        call vout(err2s, "errs")
    end do

    ! Release memory
    deallocate (wavetable)

end program tslap
