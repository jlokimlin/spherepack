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
!     6/98
!
!     a program for testing all vorticity and ivnerse vorticity routines
!
!
!     (1) first set a stream function and velocity potential scalar fields as
!         polys in x, y, z restricted to the sphere
!
!     (2) derive a vector field (v, w) from (1)
!
!     (3) compute the vorticity vt of (2) and compare with the vorticity
!         computed analytically
!
!     (4) compute vector field (ve, we) using br, bi, cr, ci from (v, w) with
!         br=bi=0.0
!
!     (5) invert the vorticity in (3) and compare with (4)
!
program tvrt

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
    real(wp) :: err2v
    real(wp) :: err2w
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: ierror
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: nmax
    real(wp) :: phi

    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: vte
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    !
    !     set dimensions with parameter statements
    integer(ip), parameter :: nlat= 24, nlon= 14, nt = 3
    integer(ip), parameter :: isym = 0, ityp = 0
    integer(ip), parameter :: mdab = (nlon+2)/2, mdc = (nlon+1)/2
    real(wp), dimension(mdc, nlat, nt)  :: br, bi, cr, ci
    real(wp), dimension(mdab, nlat, nt) :: a, b
    real(wp), dimension(nlat, nlon, nt) :: vt, v, w, ve, we
    real(wp), dimension(nlat)           :: thetag, dwts
    real(wp)                            :: pertrb(nt)
    real(wp), allocatable               :: wavetable(:)
    real(wp), parameter                 :: ZERO = 0.0_wp, TWO = 2.0_wp

    ! Set dimension variables
    nmax = max(nlat, nlon)
    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")
    !
    !     set equally spaced colatitude and longitude increments
    !
    dphi = TWO_PI/nlon
    dlat = PI/(nlat-1)
    !
    !     compute nlat gaussian points in thetag
    !
    call compute_gaussian_latitudes_and_weights(nlat, thetag, dwts, ierror)

    call name("gaqd")
    call iout(ierror, " error_flag = ")
    call vecout(thetag, "gaussian_latitudes", nlat)
    !
    !     test all vorticity subroutines
    !
    do icase=1, 4
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !
        !     set scalar stream and velocity potential fields as polys in x, y, z
        !     and then set v, w from st, sv scalar fields
        !
        do k=1, nt
            do j=1, nlon
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    select case (icase)
                        case (0:2)
                            theta = real(i -1, kind=wp) * dlat
                        case default
                            theta=thetag(i)
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
                    select case (k)
                        case (1)
                            !              st(i, j, k) = x
                            !              sv(i, j, k) = y
                            !
                            !          v = -1/sin(theta)*dstdp + dsvdt
                            !
                            !          w =  1/sin(theta)*dsvdp + dstdt
                            !
                            v(i, j, k) = sinp + cost * sinp
                            w(i, j, k) = cosp + cost * cosp
                            vt(i, j, k) = -TWO * sint * cosp
                        case (2)
                            !              st = y
                            !              sv = z
                            v(i, j, k) = -cosp-sint
                            w(i, j, k) = cost * sinp
                            !         sint*vt = -dvdp + sint*dwdt + cost*w
                            !                 = sinp + sint*(-sint*sinp)+cost*cost*sinp
                            !                 = sinp + (cost**2-sint**2)*sinp
                            vt(i, j, k) = -TWO * sint * sinp
                        case (3)
                            !           st = x
                            !           sv = z
                            v(i, j, k) = sinp - sint
                            w(i, j, k) = cost * cosp
                            !     sint*vt = -cosp-sint*sint*sinp+cost*cost*cosp
                            !             = -cosp + (1-2.*sint**2)*cosp =
                            vt(i, j, k) = -TWO * sint * cosp
                    end select
                end do
            end do
        end do

        if (icase==1) then

            call name("**ec")
            !
            !     analyze vector field
            !
            call initialize_vhaec(nlat, nlon, wavetable, ierror)
            call name("vhai")
            call iout(ierror, "error_flag = ")

            call vhaec(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdc, &
                nlat, wavetable, ierror)
            call name("vha ")
            call iout(ierror, "error_flag = ")

            !     compute vorticity of (v, w) in vt
            !

            call initialize_shsec(nlat, nlon, wavetable, ierror)

            call name("vrti")
            call iout(ierror, "error_flag = ")

            call vrtec(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wavetable, ierror)

            call name("vrt ")
            call iout(ierror, "error_flag = ")
            call iout(nlat, "nlat")
            call iout(nlon, "nlon")

        else if (icase==2) then

            call name("**es")
            call initialize_shses(nlat, nlon, wavetable, ierror)

            call name("vrti")
            call iout(ierror, "error_flag = ")

            call vrtes(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wavetable, ierror)

            call name("vrt ")
            call iout(ierror, "error_flag = ")

        else if (icase == 3) then

            call name("**gc")

            call initialize_shsgc(nlat, nlon, wavetable, ierror)

            call name("vrti")
            call iout(ierror, "error_flag = ")

            call vrtgc(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wavetable, ierror)

            call name("vrt ")
            call iout(ierror, "error_flag = ")

        else if (icase == 4) then

            call name("**gs")

            call initialize_shsgs(nlat, nlon, wavetable, ierror)

            call name("vrti")
            call iout(ierror, "error_flag = ")

            call vrtgs(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wavetable, ierror)

            call name("vrt ")
            call iout(ierror, "error_flag = ")
        end if

        !     compute "error" in vt
        !
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
                            theta=thetag(i)
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
                    select case (k)
                        case (1, 3)
                            vte = -TWO * sint * cosp
                        case (2)
                            vte = -TWO * sint * sinp
                    end select
                    err2 = err2 + (vt(i, j, k)-vte)**2
                end do
            end do
        end do
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2, "err2")

        !     now recompute (v, w) inverting vt using ivrt(ec, es, gc, gs)
        !     and compare with (ve, we) generated by synthesizing br, bi, cr, ci
        !     with br=bi=0.0
        br = ZERO
        bi = ZERO

        select case (icase)
            case (1)
        		
                call name("**ec")
        		
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call initialize_vhsec(nlat, nlon, wavetable, ierror)
                call name("initialize_vhsec")
                call iout(ierror, "error_flag = ")
        		
                call vhsec(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wavetable, ierror)
        		
                call name("vhsec")
                call iout(ierror, "error_flag = ")
        		
                call initialize_shaec(nlat, nlon, wavetable, ierror)
                call name("initialize_shaec")
                call iout(ierror, "error_flag = ")
        		
                call shaec(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, ierror)
                call name("shaec")
                call iout(ierror, "error_flag = ")


        		
        		
                call initialize_vhsec(nlat, nlon, wavetable, ierror)
                call name("initialize_vhsec")
                call iout(ierror, "error_flag = ")
        		
                call ivrtec(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, pertrb, ierror)
                call name("ivrtec")
                call iout(ierror, "error_flag = ")
                call vout(pertrb, "prtb")
            case (2)
        		
                call name("**es")
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call initialize_vhses(nlat, nlon, wavetable, ierror)
                call name("vhsi")
                call iout(ierror, "error_flag = ")
        		
                call vhses(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wavetable, ierror)
                call name("vhs ")
                call iout(ierror, "error_flag = ")
        		
        		
                call initialize_shaes(nlat, nlon, wavetable, ierror)
                call name("shai")
                call iout(ierror, "error_flag = ")
        		
                call shaes(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, ierror)
                call name("sha ")
                call iout(ierror, "error_flag = ")


        		
                call initialize_vhses(nlat, nlon, wavetable, ierror)
                call name("ivti")
                call iout(ierror, "error_flag = ")
        		
                call ivrtes(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, pertrb, ierror)
                call name("ivrt")
                call iout(ierror, "error_flag = ")
                call vout(pertrb, "prtb")
            case (3)
        		
                call name("**gc")
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call initialize_vhsgc(nlat, nlon, wavetable, ierror)
                call name("initialize_vhsgc")
                call iout(ierror, "error_flag = ")
        		
                call vhsgc(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wavetable, ierror)
        		
                call name("vhsgc")
                call iout(ierror, "error_flag = ")
        		
                call initialize_shagc(nlat, nlon, wavetable, ierror)
                call name("initialize_shagc")
                call iout(ierror, "error_flag = ")
        		
                call shagc(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, ierror)
                call name("shagc")
                call iout(ierror, "error_flag = ")


        		
                call initialize_vhsgc(nlat, nlon, wavetable, ierror)
                call name("initialize_vhsgc")
                call iout(ierror, "error_flag = ")
        		
                call ivrtgc(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, pertrb, ierror)
        		
                call name("ivrtgc")
                call iout(ierror, "error_flag = ")
                call vout(pertrb, "prtb")
            case (4)
        		
                call name("**gs")
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call initialize_vhsgs(nlat, nlon, wavetable, ierror)
                call name("vhsi")
                call iout(ierror, "error_flag = ")
        		
                call vhsgs(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wavetable, ierror)
                call name("vhs ")
                call iout(ierror, "error_flag = ")
        		
                call initialize_shags(nlat, nlon, wavetable, ierror)
                call name("shai")
                call iout(ierror, "error_flag = ")
        		
                call shags(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, ierror)
                call name("sha ")
                call iout(ierror, "error_flag = ")


        		
                call initialize_vhsgs(nlat, nlon, wavetable, ierror)
                call name("ivti")
                call iout(ierror, "error_flag = ")
        		
                call ivrtgs(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wavetable, pertrb, ierror)
                call name("ivrt")
                call iout(ierror, "error_flag = ")
                call vout(pertrb, "prtb")
        end select

        !
        !     compare this v, w with ve, we
        !
        err2v = ZERO
        err2w = ZERO
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    err2v = err2v + (v(i, j, k)-ve(i, j, k))**2
                    err2w = err2w + (w(i, j, k)-we(i, j, k))**2
                end do
            end do
        end do
        err2v = sqrt(err2v/(nlat*nlon*nt))
        err2w = sqrt(err2w/(nlat*nlon*nt))
        call vout(err2v, "errv")
        call vout(err2w, "errw")
    !
    !     end of icase loop
    !
    end do

    ! Release memory
    deallocate (wavetable)

end program tvrt
