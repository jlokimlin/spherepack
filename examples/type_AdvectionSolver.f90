module type_AdvectionSolver

    use, intrinsic :: ISO_Fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        GaussianSphere, &
        PI

    ! Explicit typing only
    implicit none

    private
    public :: DEGREES_TO_RADIANS
    public :: TIME_TO_CIRCUMVENT_THE_EARTH

    !----------------------------------------------------------------------
    ! Dictionary: global variables
    !----------------------------------------------------------------------
    real(wp),    parameter :: DEGREES_TO_RADIANS = PI/180
    integer(ip), parameter :: TIME_TO_CIRCUMVENT_THE_EARTH = 12*24*3600
    !----------------------------------------------------------------------

    
    type, public, extends(GaussianSphere) :: AdvectionSolver
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
        real(wp) :: TIME_STEP = 6.0e+2_wp
        real(wp) :: ROTATION_RATE_OF_EARTH = 2.0_wp * PI/TIME_TO_CIRCUMVENT_THE_EARTH
        real(wp) :: LATITUDE_OF_COSINE_BELL = PI/6
        real(wp) :: RADIUS_OF_EARTH_IN_METERS = 1.0_wp/3
        real(wp) :: TILT_ANGLE_IN_DEGREES = 60.0_wp
        real(wp) :: TILT_ANGLE = 60.0_wp * DEGREES_TO_RADIANS
        real(wp) :: MAXIMUM_VALUE_OF_COSINE_BELL = 1.0e+3_wp
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
        !----------------------------------------------------------------------
        procedure, public :: get_geopotential
        procedure, public :: get_vector_velocities
        procedure, nopass :: atanxy
        !----------------------------------------------------------------------
    end type AdvectionSolver


contains


    subroutine get_geopotential(self, t, geopot)
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
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(AdvectionSolver), intent(inout)  :: self
        real(wp),               intent(in)     :: t
        real(wp),               intent(out)    :: geopot(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: i, j
        real(wp)    :: xc, yc, zc, x1, y1, z1
        real(wp)    :: lambda, theta, cost, sint, sth ,cthclh
        real(wp)    :: cthslh,lhat,cosl,sinl,cth ,that, r
        !----------------------------------------------------------------------

        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            omega => self%ROTATION_RATE_OF_EARTH, &
            beta => self%LATITUDE_OF_COSINE_BELL, &
            re => self%RADIUS_OF_EARTH_IN_METERS, &
            alpha => self%TILT_ANGLE, &
            hzero => self%MAXIMUM_VALUE_OF_COSINE_BELL,&
            colat => self%grid%latitudes &
            )

            associate( lambdc => omega*t )

                call sph2cart(1.0_wp, beta, lambdc, xc, yc, zc)

            end associate


            associate( &
                cosa => cos(alpha), &
                sina => sin(alpha) &
                )
                do j=1,nlon
                    lambda = self%grid%longitudes(j)
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
                            !  compute distance
                            !
                            associate( dist => norm2([x1-xc, y1-yc, z1-zc]) )
                                !
                                !  Initialize geopotential
                                !
                                geopot(i,j) = 0.0_wp

                                if (dist >= re) cycle

                                r = 2.0_wp * asin(dist/2)

                            end associate

                            if (r >= re) cycle

                            !
                            !  Set geopotential
                            geopot(i,j) = hzero * (cos(r*PI/re)+1.0_wp)/2
                        end do
                    end associate
                end do
            end associate
        end associate

    end subroutine get_geopotential


    subroutine get_vector_velocities(self, u, v)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(AdvectionSolver), intent(inout)  :: self
        real(wp),               intent(out)    :: u(:,:)
        real(wp),               intent(out)    :: v(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: i, j
        real(wp)    :: sinp, cosp, sint, cost, sinl, cosl
        real(wp)    :: cth, xlhat, uhat
        !----------------------------------------------------------------------

        associate( &
            nlats => self%NUMBER_OF_LATITUDES, &
            nlons => self%NUMBER_OF_LONGITUDES, &
            colat => self%grid%latitudes, &
            omega => self%ROTATION_RATE_OF_EARTH, &
            alpha => self%TILT_ANGLE &
            )

            do j=1,nlons
                associate( xlm => self%grid%longitudes(j) )
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

                            xlhat = self%atanxy(cthclh,cthslh)
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
        ! Dummy arguments
        !----------------------------------------------------------------------
        real(wp), intent(in) :: x
        real(wp), intent(in) :: y
        real                   :: return_value
        !--------------------------------------------------------------------------------

        return_value = 0.0_wp

        if (x == 0.0_wp .and. y == 0.0_wp) return

        return_value = atan2(y,x)

    end function atanxy


    pure subroutine sph2cart(r, theta, phi, x, y, z)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        real(wp), intent(in)  :: r
        real(wp), intent(in)  :: theta
        real(wp), intent(in)  :: phi
        real(wp), intent(out) :: x, y, z
        !----------------------------------------------------------------------

        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)

    end subroutine sph2cart


end module type_AdvectionSolver
