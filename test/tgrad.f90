!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Spherepack                            *
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
!     12/96
!
!     a program for testing all gradient and inverse gradient routines
!
!     (1) first a scalar field is set in st by restricting a poly in x, y, z
!         to the sphere surface
!
!     (2) a scalar analysis is used to compute the coefs a, b of st
!
!     (3) a, b are input to the various gradient routines to compute a vector field
!         (v, w)
!
!     (4) the vector field (v, w) is compared with the gradient of st obtained
!         analytically
!
!     (5) the inverse gradient of (v, w) is computed and compared with (1)
!
program tgrad

    use spherepack

    ! Explicit typing only
    implicit none

    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: dsfdp
    real(wp) :: dsfdt
    real(wp) :: dxdp
    real(wp) :: dxdt
    real(wp) :: dydp
    real(wp) :: dydt
    real(wp) :: dzdp
    real(wp) :: dzdt
    real(wp) :: err2s
    real(wp) :: err2v
    real(wp) :: err2w
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: ierror
    integer(ip) :: j
    integer(ip) :: k
    real(wp) :: phi
    real(wp) :: se
    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: ve
    real(wp) :: we
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    !
    !     set dimensions with parameter statements
    !
    integer(ip), parameter              :: nlat = 33, nlon = 18, nt = 4
    integer(ip), parameter              :: isym = 0, ityp = 0
    integer(ip), parameter              :: mdb = (nlon + 1)/2, mdab = (nlon+2)/2
    real(wp), dimension(mdb, nlat, nt)  :: br, bi, cr, ci
    real(wp), dimension(mdab, nlat, nt) :: a, b
    real(wp), dimension(nlat, nlon, nt) :: sf, v, w
    real(wp), dimension(nlat)           :: gaussian_latitudes, gaussian_weights
    real(wp), allocatable               :: wavetable(:)

    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")

    ! Set equally spaced colatitude and longitude increments
    dphi = TWO_PI/nlon
    dlat = PI/(nlat-1)

    ! Compute nlat-many gaussian latitudinal points
    call compute_gaussian_latitudes_and_weights(nlat, gaussian_latitudes, gaussian_weights, ierror)

    call name("gaqd")
    call iout(ierror, " error_flag")
    call vecout(gaussian_latitudes, "gaussian_latitudes", nlat)
    !
    !     test all analysis and synthesis subroutines
    !
    do icase=1, 4
        !
        !     icase=1 test gradec, igradec
        !     icase=2 test grades, igrades
        !     icase=3 test gradgc, igradgc
        !     icase=4 test gradgs, igradgs
        !
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !
        !     set scalar field as polys in x, y, z and set vector field (v, w) using
        !     v = dsfdt, w = 1/sint*dsfdp
        !
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    dxdp = -sint*sinp
                    dydt = cost*sinp
                    dydp = sint*cosp
                    dzdt = -sint
                    dzdp = 0.0
                    select case (k)
                        case (1)
                            sf(i, j, k) = x*y
                            dsfdt = x*dydt+y*dxdt
                            dsfdp = x*dydp+y*dxdp
                            v(i, j, k) = dsfdt
                            w(i, j, k) = (cosp*dydp+sinp*dxdp)
                        case (2)
                            sf(i, j, k) = x*z
                            dsfdp = x*dzdp+z*dxdp
                            dsfdt = x*dzdt+z*dxdt
                            v(i, j, k) = dsfdt
                            w(i, j, k) = cosp*dzdp-z*sinp
                        case (3)
                            sf(i, j, k) = y*z
                            dsfdt = y*dzdt + z*dydt
                            dsfdp = y*dzdp+ z*dydp
                            v(i, j, k) = dsfdt
                            w(i, j, k) = sinp*dzdp + z*cosp
                        case (4)
                            sf(i, j, k) = x*y*z
                            dsfdt = x*y*dzdt + x*z*dydt + y*z*dxdt
                            dsfdp = x*y*dzdp + x*z*dydp + y*z*dxdp
                            v(i, j, k) = dsfdt
                            w(i, j, k) = cosp*y*dzdp+cosp*z*dydp+sinp*z*dxdp
                    end select
                end do
            end do
        end do

        select case(icase)
            case(1)

                call name("**ec")
                !
                !     analyze scalar field st
                !
                call initialize_shaec(nlat, nlon, wavetable, ierror)
                call name("shai")
                call iout(ierror, "error_flag = ")

                call shaec(nlat, nlon, isym, nt, sf, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("sha ")
                call iout(ierror, "error_flag = ")

                !     compute gradient of st in v, w
                call initialize_vhsec(nlat, nlon, wavetable, ierror)
                call name("vhci")
                call iout(ierror, "error_flag = ")

                call gradec(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("grad")
                call iout(ierror, "error_flag = ")

            case(2)

                call name("**es")
                !
                !     analyze scalar field st
                !
                call initialize_shaes(nlat, nlon, wavetable, ierror)
                call name("shai")
                call iout(ierror, "error_flag = ")

                call shaes(nlat, nlon, isym, nt, sf, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("sha ")
                call iout(ierror, "error_flag = ")

                !     compute gradient of st in v, w
                call initialize_vhses(nlat, nlon, wavetable, ierror)
                call name("vhsi")
                call iout(ierror, "error_flag = ")

                call grades(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("grad")
                call iout(ierror, "error_flag = ")

            case(3)


                call name("**gc")
                call initialize_shagc(nlat, nlon, wavetable, ierror)
                call name("shai")
                call iout(ierror, "error_flag = ")

                call shagc(nlat, nlon, isym, nt, sf, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("sha ")
                call iout(ierror, "error_flag = ")


                !     compute gradient of st in v, w
                call initialize_vhsgc(nlat, nlon, wavetable, ierror)
                call name("vhgc")
                call iout(ierror, "error_flag = ")

                call gradgc(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("grad")
                call iout(ierror, "error_flag = ")


            case(4)

                call name("**gs")

                call initialize_shags(nlat, nlon, wavetable, ierror)
                call name("shai")
                call iout(ierror, "error_flag = ")

                call shags(nlat, nlon, isym, nt, sf, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("sha ")
                call iout(ierror, "error_flag = ")

                !     compute gradient of st in v, w
                call initialize_vhsgs(nlat, nlon, wavetable, ierror)
                call name("vhgs")
                call iout(ierror, "error_flag = ")

                call gradgs(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, mdab, nlat, wavetable, ierror)
                call name("grad")
                call iout(ierror, "error_flag = ")
        end select

        !     compute "error" in v, w
        err2v = 0.0
        err2w = 0.0
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    dxdp = -sint*sinp
                    dydt = cost*sinp
                    dydp = sint*cosp
                    dzdt = -sint
                    dzdp = 0.0
                    select case (k)
                        case (1)
                            dsfdt = x*dydt+y*dxdt
                            dsfdp = x*dydp+y*dxdp
                            ve = dsfdt
                            we = (cosp*dydp+sinp*dxdp)
                        case (2)
                            dsfdp = x*dzdp+z*dxdp
                            dsfdt = x*dzdt+z*dxdt
                            ve = dsfdt
                            we = cosp*dzdp-z*sinp
                        case (3)
                            dsfdt = y*dzdt + z*dydt
                            dsfdp = y*dzdp+ z*dydp
                            ve = dsfdt
                            we = sinp*dzdp + z*cosp
                        case (4)
                            dsfdt = x*y*dzdt + x*z*dydt + y*z*dxdt
                            dsfdp = x*y*dzdp + x*z*dydp + y*z*dxdp
                            ve = dsfdt
                            we = cosp*y*dzdp+cosp*z*dydp+sinp*z*dxdp
                    end select
                    err2v = err2v + (v(i, j, k)-ve)**2
                    err2w = err2w + (w(i, j, k)-we)**2
                end do
            end do
        end do
        !
        !     set and print least squares error in v, w
        !
        err2v = sqrt(err2v/(nt*nlat*nlon))
        err2w = sqrt(err2w/(nt*nlat*nlon))
        call vout(err2v, "errv")
        call vout(err2w, "errw")
        !
        !     now recompute sf by inverting (v, w) using igrad(ec, es, gc, gs)
        !
        select case (icase)
            case(1)

                call initialize_vhaec(nlat, nlon, wavetable, ierror)
                call name("vhai")
                call iout(ierror, "error_flag = ")

                call vhaec(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, &
                    mdb, nlat, wavetable, ierror)
                call name("vha ")
                call iout(ierror, "error_flag = ")

                call initialize_shsec(nlat, nlon, wavetable, ierror)
                call name("shec")
                call iout(ierror, "error_flag = ")

                call igradec(nlat, nlon, isym, nt, sf, nlat, nlon, br, bi, &
                    mdb, nlat, wavetable, ierror)
                call name("igra")
                call iout(ierror, "error_flag = ")

            case(2)
                !
                !     analyze vector field (v, w)
                !
                call initialize_vhaes(nlat, nlon, wavetable, ierror)
                call name("vhai")
                call iout(ierror, "error_flag = ")

                call vhaes(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, &
                    mdb, nlat, wavetable, ierror)
                call name("vha ")
                call iout(ierror, "error_flag = ")

                call initialize_shses(nlat, nlon, wavetable, ierror)
                call name("shes")
                call iout(ierror, "error_flag = ")

                call igrades(nlat, nlon, isym, nt, sf, nlat, nlon, br, bi, &
                    mdb, nlat, wavetable, ierror)
                call name("igra")
                call iout(ierror, "error_flag = ")

            case(3)

                call initialize_vhagc(nlat, nlon, wavetable, ierror)
                call name("vhai")
                call iout(ierror, "error_flag = ")

                call vhagc(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, &
                    mdb, nlat, wavetable, ierror)
                call name("vha ")
                call iout(ierror, "error_flag = ")

                call initialize_shsgc(nlat, nlon, wavetable, ierror)
                call name("shgc")
                call iout(ierror, "error_flag = ")

                call igradgc(nlat, nlon, isym, nt, sf, nlat, nlon, br, bi, &
                    mdb, nlat, wavetable, ierror)
                call name("igra")
                call iout(ierror, "error_flag = ")

            case(4)
                call initialize_vhags(nlat, nlon, wavetable, ierror)
                call name("vhai")
                call iout(ierror, "error_flag = ")

                call vhags(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, &
                    mdb, nlat, wavetable, ierror)
                call name("vha ")
                call iout(ierror, "error_flag = ")

                call initialize_shsgs(nlat, nlon, wavetable, ierror)
                call name("shgs")
                call iout(ierror, "error_flag = ")

                call igradgs(nlat, nlon, isym, nt, sf, nlat, nlon, br, bi, &
                    mdb, nlat, wavetable, ierror)
                call name("igra")
                call iout(ierror, "error_flag = ")

        end select

        !     compare this sf with original
        !
        err2s = 0.0
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase > 2) theta=gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    select case (k)
                        case (1)
                            se = x*y
                        case (2)
                            se = x*z
                        case (3)
                            se = y*z
                        case (4)
                            se = x*y*z
                    end select
                    err2s = err2s + (sf(i, j, k)-se)**2
                end do
            end do
        end do
        err2s = sqrt(err2s/(nlat*nlon*nt))
        call vout(err2s, "errs")
    end do

    ! Release memory
    deallocate (wavetable)

end program tgrad

