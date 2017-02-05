!
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!  .                                                             .
!  .                  copyright (c) 1998 by UCAR                 .
!  .                                                             .
!  .       University Corporation for Atmospheric Research       .
!  .                                                             .
!  .                      all rights reserved                    .
!  .                                                             .
!  .                                                             .
!  .                          spherepack                         .
!  .                                                             .
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
!
!
! ... file advec.f90
!
!     program advec solves the time-dependent linear advection 
!     equation for geopotential phi using the spherepack software
!
!          d(phi)/dt = -(u, v) .dot. gradient(phi)
!
!                    = -(u*grad_phi_lon + v*grad_phi_lat)
!
! ... required files
!
!     gradgc.f90, shagc.f90, shsgc.f90, vhsgc.f90, type_SpherepackAux.f90 type_RealPeriodicFastFourierTransform.f90, compute_gaussian_latitudes_and_weights.f90
!
!
! definitions:
!
!
!     nlat          number of gaussian latitudes excluding poles
!     nlon          number of distinct longitudes
!     omega         rotation rate of earth in radians per second
!     alpha         angle between axis of rotation and the coordinate
!                   axis
!     beta          latitude of the cosine bell
!     aa            radius of earth in meters
!     ncycle        exit number
!     time          model time in seconds
!     dt            time step
!     lambda        longitude
!     theta         latitude
!
!   the first dimension of the following two dimensional arrays
!   corresponds to the latitude index with values i=1, ..., nlat
!   where i=1 is the northern most gaussian point thetag(i)
!   and i=nlat is the southern most gaussian point thetag(nlat).
!   the second dimension is longitude with values j=1, ..., nlon
!   where j=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!
!     thetag(i)           vector of gaussian points on the full sphere which
!                         have north to south orientation as i=1, ..., nlat
!
!     u(i, j)               east longitudinal velocity component
!
!     v(i, j)               latitudinal velocity component
!
!     phi(i, j)             the geopotential at t = time
!
!     phi_new(i, j)         the geopotential at t=time+dt
!
!     phi_old(i, j)         the geopotential at t=time-dt
!
!     grad_phi_lon(i, j)    the longitudinal derivative component of
!                          the gradient of phi
!
!                          grad_phi_lon = 1/(cos(theta))*d(phi)/dlambda
!
!
!     grad_phi_lat(i, j)   the latitudinal derivative component of
!                         the gradient of phi
!
!                         grad_phi_lat = d(phi)/dtheta
!
!
!   the following two dimensional arrays are nonzero in the triangle
!   n=1, ..., nlat and m less than or equal to n.
!
!     ar(m, n), br(m, n)    spectral coefficients of phi
!
program advec

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        wp, & ! Working precision
        ip    ! Integer precision

    use type_AdvectionSolver, only: &
        AdvectionSolver, &
        TIME_TO_CIRCUMVENT_THE_EARTH

    ! Explicit typing only
    implicit none

    ! Dictionary
    type(AdvectionSolver)             :: solver
    integer(ip), parameter            :: NLONS = 45
    integer(ip), parameter            :: NLATS = 23
    integer(ip)                       :: i, j ! Counters
    integer(ip)                       :: mprint, ncycle, ntime
    real(wp), dimension(NLATS, NLONS) :: u, v, phi, phi_old, phi_new
    real(wp), dimension(NLATS, NLONS) :: exact_phi, grad_phi, grad_phi_lon, grad_phi_lat
    real(wp)                          :: p0_l2, p0_max, time, htime
    real(wp),         parameter       :: ZERO = 0.0_wp, HALF = 0.5_wp
    character(len=*), parameter       :: WRITE_FMT = &
        '(a, i10, a, f10.2/, a, f10.0, a, i10/, a, i10, '&
        //'a, 1pe15.6/, a, 1pe15.6, a, 0pf10.2/a, 1pe15.6, a, 1pe15.6)'

    write( stdout, '(/a)') '     advec *** TEST RUN *** '

    ! Allocate memory
    call solver%create(nlat=NLATS, nlon=NLONS)

    ! Set vector velocities and cosine bell in geopotential
    call solver%get_vector_velocities(u, v)

    ! Compute geopotential at t=-dt in phi_old and at t=0.0 in phi
    ! to start up leapfrog scheme
    associate( dt => solver%TIME_STEP )
        call solver%get_geopotential(-dt, phi_old)
        call solver%get_geopotential(ZERO, phi)
    end associate

    ! Smooth geopotential at t=-dt and t=0. by synthesizing after analysis
    call solver%perform_scalar_analysis(phi_old)
    call solver%perform_scalar_synthesis(phi_old)
    call solver%perform_scalar_analysis(phi)
    call solver%perform_scalar_synthesis(phi)

    ! Compute l2 and max norms of geopotential at t=0.
    p0_l2 = norm2(phi)
    p0_max = maxval(abs(phi))

    ! Set number of time steps for 12 days (time to circumvent the earth)
    associate( dt => solver%TIME_STEP )
        ntime = int(real(TIME_TO_CIRCUMVENT_THE_EARTH, kind=wp)/dt + HALF, kind=ip)
        mprint = ntime/12
        time = ZERO
        ncycle = 0
    end associate

    !  Start time loop
    time_loop: do while (ncycle <= ntime)

        !  Compute gradient of phi at current time
        call solver%get_gradient(phi, grad_phi_lat, grad_phi_lon)

        !  Compute the time derivative of phi, note that the sign
        !  of the last term is positive because the gradient is
        !  computed with respect to colatitude rather than latitude.
        grad_phi = -u * grad_phi_lon + v * grad_phi_lat

        ! Write variables to standard output
        if (mod(ncycle, mprint) == 0) then

            ! Compute exact solution
            call solver%get_geopotential(time, exact_phi)

            !  Compute errors
            associate( &
                htime => time/3600, &
                dt => solver%TIME_STEP, &
                nlats => solver%NUMBER_OF_LATITUDES, &
                nlons => solver%NUMBER_OF_LONGITUDES, &
                omega => solver%ROTATION_RATE_OF_EARTH, &
                hzero => solver%MAXIMUM_VALUE_OF_COSINE_BELL, &
                alphad => solver%TILT_ANGLE_IN_DEGREES, &
                errm => maxval(abs(exact_phi-phi))/p0_max, &
                err2 => norm2(exact_phi-phi)/p0_l2 &
                )

                write( stdout, '(/a)' ) ' advecting cosine bell, test case 2'
                write( stdout, fmt = WRITE_FMT ) &
                    ' exit number              '  , ncycle, &
                    '  model time in  hours      ', htime, &
                    ' time step in seconds      ' , dt, &
                    ' number of latitudes       ' , nlats, &
                    ' number of longitudes      ' , nlons, &
                    ' rotation rate        '      , omega, &
                    ' mean height          '      , hzero, &
                    ' tilt angle                ' , alphad, &
                    ' max geopot. error    '      , errm, &
                    ' RMS geopot. error    '      , err2
            end associate
        end if

        !  Update various quantities
        associate( dt => solver%TIME_STEP )

            ! Increment
            time = time + dt
            ncycle = ncycle + 1

            ! Update phi_old, phi for next time step
            phi_new = phi_old + 2.0_wp * dt * grad_phi
            phi_old = phi
            phi = phi_new
        end associate
    end do time_loop

    !  Release memory
    call solver%destroy()

end program advec
