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
!
!
!
!     12/96
!
!     a program for testing all divergence and inverse divergence routines
!
!     (1) first set an irrotational vector field v, w
!
!     (2) compute the coefficients br, bi, cr, ci of (v, w) using vector analysis
!
!     (3) compute the divergence of (v, w) using divec, dives, divgc, divgs
!
!     (4) analystically compute the divergence and compare with (3)
!
!     (5) invert the divergence and compare with the irrotational (v, w)
!
program tdiv

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
    real(wp) :: dsfdp
    real(wp) :: dsfdt
    real(wp) :: dvdt
    real(wp) :: dve
    
    real(wp) :: dx2dt2
    real(wp) :: dxdp
    real(wp) :: dxdt
    real(wp) :: dydp
    real(wp) :: dydt
    real(wp) :: dzdp
    real(wp) :: dzdt
    real(wp) :: err2
    real(wp) :: err2v
    real(wp) :: err2w
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: error_flag
    integer(ip) :: j
    integer(ip) :: k
    real(wp) :: phi
    real(wp) :: sf
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
    integer(ip), parameter              :: nlat = 15, nlon =  9, nt = 2
    integer(ip), parameter              :: isym = 0
    integer(ip), parameter              :: mdab = (nlon+2)/2, mdb = (nlon + 1)/2
    real(wp), dimension(mdb, nlat, nt)  :: br, bi, cr, ci
    real(wp), dimension(mdab, nlat, nt) :: a, b
    real(wp), dimension(nlat, nlon, nt) :: dv, v, w
    real(wp), dimension(nlat)           :: gaussian_latitudes, gaussian_weights
    real(wp)                            :: pertrb(nt)
    real(wp), allocatable               :: wavetable(:)
    real(wp), parameter                 :: ZERO = 0.0_wp, TWO = 2.0_wp, SIX = 6.0_wp

    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")

    ! Set equally spaced colatitude and longitude increments
    dphi = TWO_PI/nlon
    dlat = PI/(nlat-1)

    ! Compute nlat-many gaussian latitudinal points
    call compute_gaussian_latitudes_and_weights(nlat, gaussian_latitudes, gaussian_weights, error_flag)

    call name("compute_gaussian_latitudes_and_weights")
    call iout(error_flag, " ier")
    call vecout(gaussian_latitudes, "gaussian_latitudes", nlat)
    !
    !     test all divergence and inverse divergence subroutines
    !
    do icase=1, 4
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !     set vector field v, w
        !
        do k=1, nt
            do j=1, nlon
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    select case (icase)
                        case (0:2)
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
                    d2xdt2 = -sint*cosp
                    dxdp = -sint*sinp
                    d2xdp2 = -sint*cosp
                    dydt = cost*sinp
                    d2ydt2 = -sint*sinp
                    dydp = sint*cosp
                    d2ydp2 = -sint*sinp
                    dzdt = -sint
                    d2zdt2 = -cost
                    dzdp = ZERO
                    d2zdp2 = ZERO
                    if (k==1) then
                        sf = x*y
                        dsfdt = x*dydt+y*dxdt
                        dsfdp = x*dydp+y*dxdp
                        v(i, j, k) = dsfdt
                        w(i, j, k) = cosp*dydp+sinp*dxdp
                        !              dv = 1/sint*(d(sint*v)/dt + dw/dp)
                        dvdt = x*d2ydt2 + TWO * dxdt*dydt + y*d2xdt2
                        !              1/sint*dwdp = 1/sint*(cosp*d2ydp2-sinp*dydp+sinp*d2xdp2+cosp*dxdp)
                        !                          = -4.*sinp*cosp
                        dv(i, j, k) = dvdt + cost*(cosp*dydt+sinp*dxdt) -4.*cosp*sinp
                    else if (k==2) then
                        sf = x*z
                        dsfdt = x*dzdt+z*dxdt
                        dsfdp = x*dzdp+z*dxdp
                        v(i, j, k) = dsfdt
                        w(i, j, k) = -cost*sinp
                        dvdt = x*d2zdt2+TWO * dzdt*dxdt + z*dx2dt2
                        !              dv = 1/sint*(d(sint*v)/dt + dw/dp)
                        dv(i, j, k) = -SIX * cost*sint*cosp
                    else if (k==3) then
                        sf = y*z
                        dsfdt = y*dzdt+z*dydt
                        dsfdp = y*dzdp+z*dydp
                        v(i, j, k) = dsfdt
                        w(i, j, k) = z*cosp
                        dv(i, j, k) = -SIX * cost*sint*sinp
                    end if
                end do
            end do
        end do

        if (icase==1) then

            call name("**ec")
            !
            !     analyze vector field
            !
            call initialize_vhaec(nlat, nlon, wavetable, error_flag)
            call name("shs ")
            call iout(error_flag, "error_flag = ")

            call vhaec(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdb, &
                nlat, wavetable, error_flag)
            call name("vha ")
            call iout(error_flag, "error_flag = ")

            !     compute divergence of (v, w) in dv
            call initialize_shsec(nlat, nlon, wavetable, error_flag)
            call name("shsi")
            call iout(error_flag, "error_flag = ")

            call divec(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wavetable, error_flag)
            call name("div ")
            call iout(error_flag, "error_flag = ")
            call iout(nlat, "nlat")
            call iout(nlon, "nlon")


        else if (icase==2) then

            call name("**es")
            call initialize_shses(nlat, nlon, wavetable, error_flag)
            call name("shsi")
            call iout(error_flag, "error_flag = ")

            call dives(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wavetable, error_flag)

            call name("div ")
            call iout(error_flag, "error_flag = ")

        else if (icase == 3) then

            call name("**gc")

            call initialize_shsgc(nlat, nlon, wavetable, error_flag)
            call name("shsi")
            call iout(error_flag, "error_flag = ")

            call divgc(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wavetable, error_flag)
            call name("div ")
            call iout(error_flag, "error_flag = ")

        else if (icase == 4) then

            call name("**gs")

            call initialize_shsgs(nlat, nlon, wavetable, error_flag)
            call name("shsi")
            call iout(error_flag, "error_flag = ")

            call divgs(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wavetable, error_flag)
            call name("div ")
            call iout(error_flag, "error_flag = ")

        end if

        err2 = ZERO
        do k=1, nt
            do j=1, nlon
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    select case (icase)
                        case (0:2)
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
                    d2xdt2 = -sint*cosp
                    dxdp = -sint*sinp
                    d2xdp2 = -sint*cosp
                    dydt = cost*sinp
                    d2ydt2 = -sint*sinp
                    dydp = sint*cosp
                    d2ydp2 = -sint*sinp
                    dzdt = -sint
                    d2zdt2 = -cost
                    dzdp = ZERO
                    d2zdp2 = ZERO
                    if (k==1) then
                        sf = x*y
                        dvdt = x*d2ydt2 + TWO * dxdt*dydt + y*d2xdt2
                        dve = dvdt + cost*(cosp*dydt+sinp*dxdt) - 4.*cosp*sinp
                    else if (k==2) then
                        !              sf = x*z
                        dve = -SIX * sint*cost*cosp
                    else if (k==3) then
                        !              sf = y*z
                        dve = -SIX * cost*sint*sinp
                    end if
                    err2 = err2 + (dv(i, j, k)-dve)**2
                end do
            end do
        end do
        !
        !     set and print least squares error in dv
        !
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2, "err2")

        !     now recompute (v, w) inverting dv using idiv(ec, es, gc, gs)
        v = ZERO
        w = ZERO


        if (icase==1) then

            call name("**ec")


            call initialize_shaec(nlat, nlon, wavetable, error_flag)
            call name("shai")
            call iout(error_flag, "error_flag = ")

            call shaec(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wavetable, error_flag)
            call name("sha ")
            call iout(error_flag, "error_flag = ")

            call initialize_vhsec(nlat, nlon, wavetable, error_flag)
            call name("idvi")
            call iout(error_flag, "error_flag = ")

            call idivec(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wavetable, pertrb, error_flag)
            call name("idiv")
            call iout(error_flag, "error_flag = ")
            call vecout(pertrb, "prtb", nt)

        else if (icase==2) then

            call name("**es")


            call initialize_shaes(nlat, nlon, wavetable, error_flag)
            call name("shai")
            call iout(error_flag, "error_flag = ")

            call shaes(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wavetable, error_flag)
            call name("sha ")
            call iout(error_flag, "error_flag = ")

            call initialize_vhses(nlat, nlon, wavetable, error_flag)
            call name("idvi")
            call iout(error_flag, "error_flag = ")

            call idives(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wavetable, pertrb, error_flag)
            call name("idiv")
            call iout(error_flag, "error_flag = ")
            call vecout(pertrb, "prtb", nt)

        else if (icase==3) then

            call name("**gc")


            call initialize_shagc(nlat, nlon, wavetable, error_flag)
            call name("shai")
            call iout(error_flag, "error_flag = ")

            call shagc(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wavetable, error_flag)
            call name("sha ")
            call iout(error_flag, "error_flag = ")

            call initialize_vhsgc(nlat, nlon, wavetable, error_flag)
            call name("idvi")
            call iout(error_flag, "error_flag = ")

            call idivgc(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wavetable, pertrb, error_flag)
            call name("idiv")
            call iout(error_flag, "error_flag = ")
            call vecout(pertrb, "prtb", nt)

        else if (icase==4) then

            call name("**gs")


            call initialize_shags(nlat, nlon, wavetable, error_flag)
            call name("shai")
            call iout(error_flag, "error_flag = ")

            call shags(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wavetable, error_flag)
            call name("sha ")
            call iout(error_flag, "error_flag = ")

            call initialize_vhsgs(nlat, nlon, wavetable, error_flag)
            call name("idvi")
            call iout(error_flag, "error_flag = ")

            call idivgs(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wavetable, pertrb, error_flag)
            call name("idiv")
            call iout(error_flag, "error_flag = ")
            call vecout(pertrb, "prtb", nt)

        end if

        !     compare this v, w with original
        err2v = ZERO
        err2w = ZERO
        do k=1, nt
            do j=1, nlon
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    select case (icase)
                        case (0:2)
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
                    dxdp = -sint*sinp
                    dydt = cost*sinp
                    dydp = sint*cosp
                    dzdt = -sint
                    dzdp = ZERO
                    if (k==1) then
                        sf = x*y
                        dsfdt = x*dydt+y*dxdt
                        dsfdp = x*dydp+y*dxdp
                        ve = dsfdt
                        we = cosp*dydp+sinp*dxdp
                    else if (k==2) then
                        sf = x*z
                        dsfdt = x*dzdt+z*dxdt
                        dsfdp = x*dzdp+z*dxdp
                        ve = dsfdt
                        we = -cost*sinp
                    else if (k==3) then
                        sf = y*z
                        dsfdt = y*dzdt+z*dydt
                        dsfdp = y*dzdp+z*dydp
                        ve = dsfdt
                        we = z*cosp
                    else if (k==4) then
                    end if
                    err2v = err2v + (v(i, j, k)-ve)**2
                    err2w = err2w + (w(i, j, k)-we)**2
                end do
            end do
        end do
        err2v = sqrt(err2v/(nlat*nlon*nt))
        err2w = sqrt(err2w/(nlat*nlon*nt))
        call vout(err2v, "errv")
        call vout(err2w, "errw")
    end do

    ! Release memory
    deallocate (wavetable)

end program tdiv

