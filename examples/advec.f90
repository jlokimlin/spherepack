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
!          d(phi)/dt = -(u,v) .dot. gradient(phi)
!
!                    = -(u*grad_phi_lon + v*grad_phi_lat)
!
! ... required files
!
!     gradgc.f90,shagc.f90,shsgc.f90,vhsgc.f90,type_SpherepackAux.f90 type_HFFTpack.f90,gaqd.f90
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
!   corresponds to the latitude index with values i=1,...,nlat
!   where i=1 is the northern most gaussian point thetag(i)
!   and i=nlat is the southern most gaussian point thetag(nlat).
!   the second dimension is longitude with values j=1,...,nlon
!   where j=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!
!     thetag(i)           vector of gaussian points on the full sphere which
!                         have north to south orientation as i=1,...,nlat
!
!     u(i,j)               east longitudinal velocity component
!
!     v(i,j)               latitudinal velocity component
!
!     phi(i,j)             the geopotential at t = time
!
!     phi_new(i,j)         the geopotential at t=time+dt
!
!     phi_old(i,j)         the geopotential at t=time-dt
!
!     grad_phi_lon(i,j)    the longitudinal derivative component of
!                          the gradient of phi
!
!                          grad_phi_lon = 1/(cos(theta))*d(phi)/dlambda
!
!
!     grad_phi_lat(i,j)   the latitudinal derivative component of
!                         the gradient of phi
!
!                         grad_phi_lat = d(phi)/dtheta
!
!
!   the following two dimensional arrays are nonzero in the triangle
!   n=1,...,nlat and m less than or equal to n.
!
!     ar(m,n),br(m,n)    spectral coefficients of phi
!
module type_AdvectionSolver

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        GaussianSphere

    ! Explicit typing only
    implicit none

    !----------------------------------------------------------------------
    ! Dictionary: global variables
    !----------------------------------------------------------------------
    real (wp),    parameter :: PI = acos(-1.0_wp)
    real (wp),    parameter :: DEGREES_TO_RADIANS = PI/180
    integer (ip), parameter :: TIME_TO_CIRCUMVENT_THE_EARTH = 12*24*3600
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, public, extends (GaussianSphere) :: AdvectionSolver
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        real (wp) :: TIME_STEP = 6.0e+2_wp
        real (wp) :: ROTATION_RATE_OF_EARTH = 2.0_wp * PI/TIME_TO_CIRCUMVENT_THE_EARTH
        real (wp) :: LATITUDE_OF_COSINE_BELL = PI/6
        real (wp) :: RADIUS_OF_EARTH_IN_METERS = 1.0_wp/3
        real (wp) :: TILT_ANGLE_IN_DEGREES = 60.0_wp
        real (wp) :: TILT_ANGLE = 60.0_wp * DEGREES_TO_RADIANS
        real (wp) :: MAXIMUM_VALUE_OF_COSINE_BELL = 1.0e+3_wp
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public :: get_geopotential
        procedure, public :: get_vector_velocities
        procedure, nopass :: atanxy
        !----------------------------------------------------------------------
    end type AdvectionSolver


contains


    subroutine get_geopotential(this, t, geopot)
        !
        !     computes advecting cosine bell on a tilted grid a time t.
        !
        ! input parameters
        !
        !     t      time in seconds
        !
        !     alpha  tilt angle in radians
        !
        !     beta   colatitude of cosine bell in untilted coordinate
        !            system in radians
        !
        !     omega  angular velocity in radians per second
        !
        !     hzero  maximum value of cosine bell
        !
        !     re     radius of support for cosine bell in radians
        !
        !     nlat   number of latitudes including the poles
        !
        !     nlon   number of distinct longitude lines
        !
        !     idim   first dimension of output array h
        !
        !     colat  vector of Gauss colatitude grid points
        !
        ! output parameter
        !
        !     geopot     an nlat by nlon array containing the geopotential
        !
        !             on a tilted grid
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (AdvectionSolver), intent (in out) :: this
        real (wp),               intent (in)     :: t
        real (wp),               intent (out)    :: geopot(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: i, j
        real (wp)    :: xc, yc, zc, x1, y1, z1
        real (wp)    :: lambda, theta, cost, sint, sth ,cthclh
        real (wp)    :: cthslh,lhat,cosl,sinl,cth ,that, r
        !----------------------------------------------------------------------

        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES, &
            omega => this%ROTATION_RATE_OF_EARTH, &
            beta => this%LATITUDE_OF_COSINE_BELL, &
            re => this%RADIUS_OF_EARTH_IN_METERS, &
            alpha => this%TILT_ANGLE, &
            hzero => this%MAXIMUM_VALUE_OF_COSINE_BELL,&
            colat => this%grid%latitudes &
            )

            associate( lambdc => omega*t )

                call sph2cart(1.0_wp, beta, lambdc, xc, yc, zc)

            end associate


            associate( &
                cosa => cos(alpha), &
                sina => sin(alpha) &
                )
                do j=1,nlon
                    lambda = this%grid%longitudes(j)
                    associate( &
                        cosp => cos(lambda), &
                        sinp => sin(lambda) &
                        )
                        do i=1,nlat
                            theta = colat(i)
                            cost = cos(theta)
                            sint = sin(theta)
                            sth = cosa*cost+sina*sint*cosp
                            cthclh = cosa*sint*cosp-sina*cost
                            cthslh = sint*sinp
                            lhat = atanxy(cthclh,cthslh)
                            cosl = cos(lhat)
                            sinl = sin(lhat)
                            cth = cosl*cthclh+sinl*cthslh
                            that = atanxy(sth,cth)
                            !
                            !== Compute scaled radial unit vector
                            !
                            call sph2cart(1.0_wp, that, lhat, x1, y1, z1)
                            !
                            !==> compute distance
                            !
                            associate( dist => norm2([x1-xc, y1-yc, z1-zc]) )
                                !
                                !==> Initialize geopotential
                                !
                                geopot(i,j) = 0.0_wp

                                if (dist >= re) cycle

                                r = 2.0_wp * asin(dist/2)

                            end associate

                            if (r >= re) cycle

                            !
                            !==> Set geopotential
                            geopot(i,j) = hzero * (cos(r*PI/re)+1.0_wp)/2
                        end do
                    end associate
                end do
            end associate
        end associate

    end subroutine get_geopotential


    subroutine get_vector_velocities(this, u, v)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (AdvectionSolver), intent (in out) :: this
        real (wp),               intent (out)    :: u(:,:)
        real (wp),               intent (out)    :: v(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: i, j
        real (wp)    :: sinp, cosp, sint, cost, sinl, cosl
        real (wp)    :: cth, xlhat, uhat
        !----------------------------------------------------------------------

        associate( &
            nlats => this%NUMBER_OF_LATITUDES, &
            nlons => this%NUMBER_OF_LONGITUDES, &
            colat => this%grid%latitudes, &
            omega => this%ROTATION_RATE_OF_EARTH, &
            alpha => this%TILT_ANGLE &
            )

            do j=1,nlons
                associate( xlm => this%grid%longitudes(j) )
                    sinp = sin(xlm)
                    cosp = cos(xlm)
                    do i=1,nlats
                        cost = cos(colat(i))
                        sint = sin(colat(i))
                        associate( &
                            sth => cos(alpha)*cost+sin(alpha)*sint*cosp, &
                            cthclh => cos(alpha)*sint*cosp-sin(alpha)*cost, &
                            cthslh => sint*sinp &
                            )

                            xlhat = this%atanxy(cthclh,cthslh)
                            cosl = cos(xlhat)
                            sinl = sin(xlhat)
                            cth = cosl*cthclh+sinl*cthslh
                            uhat = omega*cth
                            u(i,j) = (cos(alpha)*sinp*sinl+cosp*cosl)*uhat
                            v(i,j) = (cos(alpha)*cost*cosp*sinl &
                                -cost*sinp*cosl+sin(alpha)*sint*sinl)*uhat

                        end associate
                    end do
                end associate
            end do
        end associate

    end subroutine get_vector_velocities


    pure function atanxy(x,y) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp), intent (in) :: x
        real (wp), intent (in) :: y
        real                   :: return_value
        !--------------------------------------------------------------------------------

        return_value = 0.0_wp

        if (x == 0.0_wp .and. y == 0.0_wp) return

        return_value = atan2(y,x)

    end function atanxy


    pure subroutine sph2cart(r, theta, phi, x, y, z)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp), intent (in)  :: r
        real (wp), intent (in)  :: theta
        real (wp), intent (in)  :: phi
        real (wp), intent (out) :: x, y, z
        !----------------------------------------------------------------------

        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)

    end subroutine sph2cart


end module type_AdvectionSolver

program advec

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use type_AdvectionSolver

    ! Explicit typing only
    implicit none

    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    type (AdvectionSolver)        :: solver
    integer (ip), parameter       :: NLONS = 45
    integer (ip), parameter       :: NLATS = 23
    integer (ip)                  :: i, j !! Counters
    integer (ip)                  :: mprint, ncycle, ntime
    real (wp)                     :: u(NLATS,NLONS),v(NLATS,NLONS)
    real (wp)                     :: phi_old(NLATS,NLONS),phi_new(NLATS,NLONS)
    real (wp)                     :: phi(NLATS,NLONS)
    real (wp)                     :: exact_phi(NLATS,NLONS)
    real (wp)                     :: grad_phi(NLATS,NLONS)
    real (wp)                     :: grad_phi_lon(NLATS,NLONS),grad_phi_lat(NLATS,NLONS)
    real (wp)                     :: p0_l2, p0_max, time, htime
    character (len=*), parameter  :: write_format = &
        '( a,i10 ,a,f10.2/ ,a,f10.0 ,a,i10/ ,a,i10 '&
        //',a,1pe15.6/,a,1pe15.6 ,a,0pf10.2/ a,1pe15.6 ,a,1pe15.6)'
    !----------------------------------------------------------------------

    write( stdout, '(/A)') '     advec *** TEST RUN *** '

    !
    !==> Allocate memory
    !
    call solver%create(nlat=NLATS, nlon=NLONS)

    !
    !==> set vector velocities and cosine bell in geopotential
    !
    call solver%get_vector_velocities(u, v)



    associate( dt => solver%TIME_STEP )
        !
        !==> Compute geopotential at t=-dt in phi_old and at t=0.0 in phi
        !    to start up leapfrog scheme
        !
        call solver%get_geopotential(-dt, phi_old)
        call solver%get_geopotential(0.0_wp, phi)

    end associate


    !
    !==> smooth geopotential at t=-dt and t=0. by synthesizing after analysis
    !
    call solver%perform_scalar_analysis(phi_old)
    call solver%perform_scalar_synthesis(phi_old)
    call solver%perform_scalar_analysis(phi)
    call solver%perform_scalar_synthesis(phi)


    !
    !==> compute l2 and max norms of geopotential at t=0.
    !
    p0_l2 = norm2(phi)
    p0_max = maxval(abs(phi))


    associate( dt => solver%TIME_STEP )
        !
        !==> set number of time steps for 12 days (time to circumvent the earth)
        !
        ntime = int(real(TIME_TO_CIRCUMVENT_THE_EARTH, kind=wp)/dt + 0.5_wp, kind=ip)
        mprint = ntime/12
        time = 0.0_wp
        ncycle = 0

    end associate


    !
    !==> Start time loop
    !
    time_loop: do while (ncycle <= ntime)
        !
        !==> Compute gradient of phi at current time
        !
        call solver%get_gradient(phi, grad_phi_lat, grad_phi_lon)
        !
        !==> Compute the time derivative of phi, note that the sign
        !    of the last term is positive because the gradient is
        !    computed with respect to colatitude rather than latitude.
        !
        grad_phi = -u * grad_phi_lon + v * grad_phi_lat
        !
        !==> write variables to standard output
        !
        if (mod(ncycle,mprint) == 0) then
            !
            !==> Compute exact solution
            !
            call solver%get_geopotential(time, exact_phi)
            !
            !==> Compute errors
            !
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

                write( stdout, '(/A)' ) ' advecting cosine bell, test case 2'
                write( stdout, fmt = write_format ) &
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

        !
        !==> Update various quantities
        !
        associate( dt => solver%TIME_STEP )
            !
            !==> increment
            !
            time = time + dt
            ncycle = ncycle + 1
            !
            !==> update phi_old,phi for next time step
            !
            phi_new = phi_old + 2.0_wp * dt * grad_phi
            phi_old = phi
            phi = phi_new

        end associate
        !
        !==> end of time loop
        !
    end do time_loop

    !
    !==> Release memory
    !
    call solver%destroy()

end program advec
