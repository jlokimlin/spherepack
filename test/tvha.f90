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
!     *       A Package of Fortran77 Subroutines and Programs         *
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
!     11/96
!
!     a program for testing all vector analysis and synthesis subroutines
!
!     (1) first a scalar stream function and a scalar velocity potential function
!         are set in st, sv by restricting polys in x, y, z to the sphere surface
!
!     (2) the vector vield (v, w) is set by analytically differenting the scalar fields in (1)
!         using the standard formula relating a vector field and the stream and velocity
!         potential scalar fields in colatitude X longitude spherical coordinates
!
!          v = -1/sin(theta)*d(st)/dphi + d(sv)/dtheta
!
!          w =  1/sin(theta)*d(sv)/dphi + d(st)/dtheta
!
!     (3) a vector analysis is performed on (v, w)
!
!     (4) a vector synthesis is performed using coeffs from (3)
!
!     (5) the synthesized vector field from (4) is compared with the vector field from (2)
!
!     note:  vhaec, vhaes, vhagc, vhags, vhsec, vhses, vhsgc, vhsgs are all tested!
!
program test_all_vector_analysis_and_synthesis_routines

    use spherepack

    ! Explicit typing only
    implicit none
    
    ! Dictionary
    integer(ip), parameter              :: NLAT= 25, NLON= 19, NT = 2
    integer(ip), parameter              :: ITYP = 0
    real(wp), dimension(NLAT, NLAT, NT) :: br, bi, cr, ci
    real(wp), dimension(NLAT, NLON, NT) :: st, sv, v, w
    real(wp), dimension(NLAT)           :: gaussian_latitudes, gaussian_weights
    real(wp), allocatable               :: wavetable(:)
    real(wp)                            :: cosp, cost, dlat, dphi
    real(wp)                            :: dstdp, dstdt, dsvdp, dsvdt
    real(wp)                            :: dxdp, dxdt, dydp, dydt, dzdp, dzdt
    real(wp)                            :: err2v, err2w
    integer(ip)                         :: i, j, k, icase
    integer(ip)                         :: error_flag
    real(wp)                            :: phi, sinp, sint, theta
    real(wp)                            :: ve, we, x, y, z

    call name("Testing all vector analysis and synthesis procedures")

    ! Set dimension variables
    call iout(NLAT, "nlat")
    call iout(NLON, "nlon")
    call iout(NT, "  nt")

    ! Set equally spaced colatitude and longitude increments
    dphi = TWO_PI/NLON
    dlat = PI/(NLAT-1)

    ! Compute nlat-many gaussian latitudinal points
    call compute_gaussian_latitudes_and_weights(NLAT, gaussian_latitudes, gaussian_weights, error_flag)

    call name("gaqd")
    call iout(error_flag, " ier")
    call vecout(gaussian_latitudes, "thtg", NLAT)

    ! Test all analysis and synthesis subroutines
    do icase=1, 4

        call name("*******************************")

        ! Set scalar stream and velocity potential fields as polys in x, y, z
        ! and then set v, w from st, sv scalar fields
        do k=1, NT
            do j=1, NLON
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, NLAT
                    select case (icase)
                        case (0:2)
                            theta = real(i-1, kind=wp) * dlat
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
                    dzdp = 0.0_wp
                    select case (k)
                        case (1)
                            st(i, j, k) = x*y
                            sv(i, j, k) = y*z
                            dstdt = x*dydt+y*dxdt
                            dstdp = x*dydp+y*dxdp
                            dsvdt = y*dzdt+z*dydt
                            dsvdp = y*dzdp+z*dydp
                            v(i, j, k) = -(cosp*dydp+sinp*dxdp) + dsvdt
                            w(i, j, k) = sinp*dzdp + cost*cosp + dstdt
                        case (2)
                            st(i, j, k) = x*z
                            sv(i, j, k) = x*y
                            dstdp = x*dzdp+z*dxdp
                            dstdt = x*dzdt+z*dxdt
                            dsvdp = x*dydp+y*dxdp
                            dsvdt = x*dydt+y*dxdt

                            ! v = -1/sin(theta)*d(st)/dphi + d(sv)/dtheta
                            v(i, j, k) = z*sinp + dsvdt

                            ! w =  1/sin(theta)*d(sv)/dphi + d(st)/dtheta
                            w(i, j, k) = cosp*dydp+ sinp*dxdp + dstdt
                    end select
                end do
            end do
        end do

        select case (icase)
            case (1)
        		
                call name("testing vhaec and vhsec")
                call initialize_vhaec(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhaec")
                call iout(error_flag, "error_flag = ")
                call vhaec(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhaec")
                call iout(error_flag, "error_flag = ")
        		
                !  Now synthesize v, w from br, bi, cr, ci and compare with original
                call initialize_vhsec(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhsec")
                call iout(error_flag, "error_flag = ")
                call vhsec(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhsec")
                call iout(error_flag, "error_flag = ")
            case (2)
        		
                call name("testing vhaes and vhses")
                call initialize_vhaes(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhaes")
                call iout(error_flag, "error_flag = ")
                call vhaes(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhaes")
                call iout(error_flag, "error_flag = ")
        		
                ! Now synthesize v, w from br, bi, cr, ci and compare with original
                call initialize_vhses(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhses")
                call iout(error_flag, "error_flag = ")
                call vhses(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhses")
                call iout(error_flag, "error_flag = ")
            case (3)

                call name("testing vhagc and vhsgc")
                call iout(NLAT, "nlat")
                call initialize_vhagc(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhagc")
                call iout(error_flag, "error_flag = ")
                call vhagc(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhagc")
                call iout(error_flag, "error_flag = ")

                ! Now synthesize v, w from br, bi, cr, ci and compare with original
                call initialize_vhsgc(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhsgc")
                call iout(error_flag, "error_flag = ")
        		
                call vhsgc(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhsgc")
                call iout(error_flag, "error_flag = ")
            case (4)
        		
                call name("testing vhags and vhsgs")
                call initialize_vhags(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhags")
                call iout(error_flag, "error_flag = ")
                call vhags(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhags")
                call iout(error_flag, "error_flag = ")

                ! Now synthesize v, w from br, bi, cr, ci and compare with original
                call initialize_vhsgs(NLAT, NLON, wavetable, error_flag)
                call name("initialize_vhsgs")
                call iout(error_flag, "error_flag = ")
                call vhsgs(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wavetable, error_flag)
                call name("vhsgs")
                call iout(error_flag, "error_flag = ")
        end select

        err2v = 0.0_wp
        err2w = 0.0_wp
        do k=1, NT
            do j=1, NLON
                phi = real(j - 1, kind=wp)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, NLAT
                    select case (icase)
                        case (0:2)
                            theta = real(i-1, kind=wp) * dlat
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
                    dzdp = 0.0
                    select case (k)
                        case (1)
                            st(i, j, k) = x*y
                            sv(i, j, k) = y*z
                            dstdt = x*dydt+y*dxdt
                            dstdp = x*dydp+y*dxdp
                            dsvdt = y*dzdt+z*dydt
                            dsvdp = y*dzdp+z*dydp
                            ve = -(cosp*dydp+sinp*dxdp) + dsvdt
                            we = sinp*dzdp + cost*cosp + dstdt
                        case (2)
                            st(i, j, k) = x*z
                            sv(i, j, k) = x*y
                            dstdp = x*dzdp+z*dxdp
                            dstdt = x*dzdt+z*dxdt
                            dsvdp =  x*dydp+y*dxdp
                            dsvdt = x*dydt+y*dxdt
                            ve = z*sinp + dsvdt
                            we = cosp*dydp+ sinp*dxdp + dstdt
                    end select
                    err2v = err2v + (v(i, j, k) - ve)**2
                    err2w = err2w + (w(i, j, k) - we)**2
                end do
            end do
        end do

        ! Set and print least squares error in v, w
        err2v = sqrt(err2v/(NT*NLAT*NLON))
        err2w = sqrt(err2w/(NT*NLAT*NLON))
        call vout(err2v, "errv")
        call vout(err2w, "errw")
    end do

    ! Release memory
    deallocate (wavetable)

end program test_all_vector_analysis_and_synthesis_routines
