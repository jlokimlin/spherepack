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
!     *       A Package of Fortran subroutines and programs           *
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
!     a program for testing all scalar analysis and synthesis subroutines
!
program test_all_scalar_analysis_and_synthesis_routines

    use spherepack

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter            :: NLAT= 15, NLON= 18, NT = 3
    integer(ip), parameter            :: ISYM = 0
    real(wp)                          :: cosp, cost, dlat, dphi
    integer(ip)                       :: i, j, k, icase, error_flag
    real(wp)                          :: phi, sinp, sint
    real(wp)                          :: theta, xyzk, err2
    real(wp)                          :: scalar_function(NLAT,NLON,NT)
    real(wp), allocatable             :: wavetable(:)
    real(wp), dimension(NLAT,NLAT,NT) :: a, b
    real(wp), dimension(NLAT)         :: gaussian_latitudes, gaussian_weights

    call name("Testing all scalar analysis and synthesis procedures")

    ! Print dimension variables
    call iout(NLAT, "nlat")
    call iout(NLON, "nlon")
    call iout(NT, "  nt")

    ! Set equally spaced colatitude and longitude increments
    dphi = TWO_PI/NLON
    dlat = pi/(NLAT-1)

    ! Compute nlat-many gaussian points in thetag
    call compute_gaussian_latitudes_and_weights(NLAT, gaussian_latitudes, gaussian_weights, error_flag)

    call name("gaqd")
    call iout(error_flag, " error flag")
    call vecout(gaussian_latitudes, "gaussian_latitudes", NLAT)

    ! Test all analysis and synthesis subroutines
    do icase=1, 4

        call name("*****************************************")

        ! Set scalar field as (x*y*z)**k) restricted to the sphere
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
                    xyzk = (sint*(sint*cost*sinp*cosp))**k
                    ! scalar_function(i, j, k) = exp(xyzk)
                    scalar_function(i, j, k) = xyzk
                end do
            end do
        end do

        select case (icase)
            case (1)

                call name("testing shaec and shsec")
                call initialize_shaec(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shaec")
                call iout(error_flag, "error_flag")
                call shaec(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shaec")
                call iout(error_flag, "error_flag")

                call initialize_shsec(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shsec")
                call iout(error_flag, "error_flag")
                call shsec(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shsec")
                call iout(error_flag, "error_flag")

            case (2)

                call name("testing shaes and shses")
                call initialize_shaes(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shaes")
                call iout(error_flag, "error_flag")
                call shaes(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shaes")
                call iout(error_flag, "error_flag")

                call initialize_shses(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shses")
                call iout(error_flag, "error_flag")
                call shses(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shses")
                call iout(error_flag, "error_flag")

            case (3)
        		
                call name("testing shagc and shsgc")
                call initialize_shagc(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shagc")
                call iout(error_flag, "error_flag")
                call shagc(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shagc")
                call iout(error_flag, "error_flag")
        		
                call initialize_shsgc(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shsgc")
                call iout(error_flag, "error_flag")
                call shsgc(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shsgc")
                call iout(error_flag, "error_flag")

            case (4)
        		
                call name("testing shags and shsgs")
                call initialize_shags(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shags")
                call iout(error_flag, "error_flag")
                call shags(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shags")
                call iout(error_flag, "error_flag")
        		
                call initialize_shsgs(NLAT, NLON, wavetable, error_flag)
                call name("initialize_shsgs")
                call iout(error_flag, "error_flag")
                call shsgs(NLAT, NLON, ISYM, NT, scalar_function, NLAT, NLON, a, b, NLAT, NLAT, wavetable, error_flag)
                call name("shsgs")
                call iout(error_flag, "error_flag")
        end select

        ! Compute "error" in s
        err2 = 0.0_wp
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
                    xyzk = (sint*(sint*cost*sinp*cosp))**k
                    ! err2 = err2 +(exp(xyzk)-s(i, j, k))**2
                    err2 = err2 + (xyzk-scalar_function(i, j, k))**2
                end do
            end do
        end do
        err2 = sqrt(err2/(NT*NLAT*NLON))
        call vout(err2, "err2")
    end do

    ! Release memory
    deallocate (wavetable)

end program test_all_scalar_analysis_and_synthesis_routines
