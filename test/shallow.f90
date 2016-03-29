!
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
!     Non-linear steady-state geostropic flow in a shallow water model.
!
!     Errors should be O(10E-5) or less in 32-bit precision, O(10E-7) or less
!     in 64-bit precision.
!
!
!     The nonlinear shallow-water equations on the sphere are
!     solved using a spectral method based on the spherical harmonics.
!     the method is described in the paper:
!
! [1] p. n. swarztrauber, spectral transform methods for solving
!     the shallow-water equations on the sphere, p.n. swarztrauber, 
!     monthly weather review, vol. 124, no. 4, april 1996, pp. 730-744.
!
!     this program implements test case 3 (steady nonlinear rotated flow)
!     in the paper:
!
! [2] d.l. williamson, j.b. drake, j.j. hack, r. jakob, and
!     p.n. swarztrauber, j. comp. phys., a standard test set
!     for numerical approximations to the shallow-water
!     equations in spherical geometry, j. comp. phys., 
!     vol. 102, no. 1, sept. 1992, pp. 211-224.
!
! definitions:
!
!
!     nlat          number of latitudes
!     nlon          number of distinct longitudes
!     ntrunc        max wave number
!     omega         rotation rate of earth in radians per second
!     aa            radius of earth in meters
!     pzero         mean height of geopotential
!     uzero         maximum velocity
!     alpha         tilt angle of the rotated grid
!     ncycle        cycle number
!     time          model time in seconds
!     dt            time step
!     lambda        longitude
!     theta         colatitude
!
!   the second dimension of the following two dimensional arrays
!   corresponds to the latitude index with values j=1, ..., nlat
!   going from north to south.
!   the second dimension is longitude with values i=1, ..., nlon
!   where i=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!     u(i, j)       east longitudinal velocity component at t=time
!     v(i, j)       latitudinal velocity component at t=time
!     p(i, j)       +pzero = geopotential at t=time
!
!     divg(i, j)    divergence (d/dtheta (cos(theta) v)
!                                          + du/dlambda)/cos(theta)
!     vrtg(i, j)    vorticity  (d/dtheta (cos(theta) u)
!                                          - dv/dlambda)/cos(theta)
!
!     uxact(i, j)   the "exact" longitudinal velocity component
!     vxact(i, j)   the "exact" latitudinal  velocity component
!     pxact(i, j)   the "exact" geopotential
!
program shallow

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT, &
        compiler_version, &
        compiler_options

    use type_RegularSphere, only: &
        RegularSphere
    use type_GaussianSphere, only: &
        GaussianSphere

    ! Explicit typing only
    implicit none

    !--------------------------------------------------------------------------------
    ! Dictionary
    !--------------------------------------------------------------------------------
    integer (ip), parameter :: NLONS=128
    integer (ip), parameter :: NLATS=NLONS/2! + 1
    integer (ip), parameter :: NTRUNC=NLATS-1
    integer (ip), parameter :: NL = 90
    integer (ip), parameter :: NMDIM = (NTRUNC+1)*(NTRUNC+2)/2
    integer (ip)            :: MAXIMUM_NUMBER_OF_TIME_ITERATIONS
    integer (ip)            :: MPRINT
    integer (ip)            :: i, j !! Counters
    integer (ip)            :: cycle_number
    integer (ip)            :: n_save1, n_save2,  n_old, n_now, n_new !! Iteration counters
    real (wp)               :: phlt(361)
    real (wp), dimension(NLONS, NLATS) :: uxact
    real (wp), dimension(NLONS, NLATS) :: vxact
    real (wp), dimension(NLONS, NLATS) :: pxact
    real (wp), dimension(NLONS, NLATS) :: u
    real (wp), dimension(NLONS, NLATS) :: v
    real (wp), dimension(NLONS, NLATS) :: p
    real (wp), dimension(NLONS, NLATS) :: f
    real (wp), dimension(NLONS, NLATS) :: coslat
    real (wp), dimension(NLONS, NLATS) :: ug
    real (wp), dimension(NLONS, NLATS) :: vg
    real (wp), dimension(NLONS, NLATS) :: pg
    real (wp), dimension(NLONS, NLATS) :: vrtg
    real (wp), dimension(NLONS, NLATS) :: divg
    real (wp), dimension(NLONS, NLATS) :: scrg1
    real (wp), dimension(NLONS, NLATS) :: scrg2
    real (wp), dimension(NLATS, NLONS) :: temp, temp1, temp2
    complex (wp), dimension(NMDIM) :: vort_spec
    complex (wp), dimension(NMDIM) :: div_spec
    complex (wp), dimension(NMDIM) :: p_spec
    complex (wp), dimension(NMDIM) :: scr_spec
    complex (wp), dimension(NMDIM, 3) :: dvrtdt_spec
    complex (wp), dimension(NMDIM, 3) :: ddivdt_spec
    complex (wp), dimension(NMDIM, 3) :: dpdt_spec
    real (wp), parameter :: RADIUS_OF_EARTH_IN_METERS = 6.37122e+6_wp
    real (wp), parameter :: PI = acos( -1.0_wp )
    real (wp), parameter :: HALF_PI = 0.5_wp * PI
    real (wp), parameter :: RADIAN_UNIT = PI/180.0_wp
    real (wp), parameter :: ROTATION_RATE_OF_EARTH = 7.292e-5_wp
    real (wp), parameter :: DT = 600.0_wp
    real (wp), parameter :: TILT_ANGLE = 60.0_wp
    real (wp), parameter :: ALPHA = RADIAN_UNIT * TILT_ANGLE
    
    
    real (wp), parameter :: U_ZERO = 40.0_wp
    real (wp), parameter :: P_ZERO = 2.94e+4_wp
    real (wp), parameter :: F_ZERO = 2.0_wp * ROTATION_RATE_OF_EARTH
    real (wp) :: time
    real (wp) :: lambda_hat
    real (wp) :: u_hat
    real (wp) :: theta_hat
    real (wp) :: evmax
    real (wp) :: epmax
    real (wp) :: dvmax
    real (wp) :: dpmax
    real (wp) :: model_time_in_hours
    real (wp) :: dvgm
    real (wp) :: v2max
    real (wp) :: p2max
    real (wp) :: vmax
    real (wp)                       :: pmax
    character (len=:), allocatable  :: model_write_format, error_write_format
    type (GaussianSphere)            :: sphere
    !--------------------------------------------------------------------------------

    write( stdout, '(A)' ) ''
    write( stdout, '(A)') ' *** Test program for TYPE(GaussianSphericalHarmonic) ***'
    write( stdout, '(A)') ''
    write( stdout, '(A)') 'Non-linear steady-state geostropic flow in a shallow water model'
    write( stdout, '(A)') ''
    write( stdout, '(A, I11)') 'Triangular trunction number  = ', NTRUNC
    write( stdout, '(A, I11)') 'Number of gaussian latitudes = ', NLATS
    write( stdout, '(A)') ''

    ! Set constants
    MAXIMUM_NUMBER_OF_TIME_ITERATIONS = nint( 864.0e+2_wp * 5.0_wp/DT, kind=ip)
    MPRINT = MAXIMUM_NUMBER_OF_TIME_ITERATIONS/10

    !  ==> initialize regular (equally-spaced) spherical harmonic
    call sphere%create( &
        nlat=NLATS, nlon=NLONS, ntrunc=NTRUNC, rsphere=RADIUS_OF_EARTH_IN_METERS)

    ! ==> compute the derivative of the unrotated geopotential
    !     p as a function of latitude
    associate( rsphere => sphere%RADIUS_OF_SPHERE )
        do i = 1, NL - 2
            associate( theta => real(i, kind=wp) * (PI/(NL-1)) )
                u_hat = &
                    initialize_unrotated_longitudinal_velocity( U_ZERO, HALF_PI - theta )

                phlt(i) = &
                    cos(theta) * u_hat * ( u_hat/sin(theta) + rsphere * F_ZERO )/(NL-1)
            end associate
        end do
    end associate

    !     compute sine transform of the derivative of the geopotential
    !     for the purpose of computing the geopotential by integration
    !     see equation (3.9) in reference [1] above

    call compute_sine_transform(phlt(1:NL-2))

    !     compute the cosine coefficients of the unrotated geopotential
    !     by the formal integration of the sine series representation

    do i = 1, NL - 2
        phlt(i) = -phlt(i)/i
    end do

    !     phlt(i) contains the coefficients in the cosine series
    !     representation of the unrotated geopotential that are used
    !     below to compute the geopotential on the rotated grid.
    !
    !     compute the initial values of  east longitudinal
    !     and latitudinal velocities u and v as well as the
    !     geopotential p and coriolis f on the rotated grid.
    !


    do j = 1, NLONS
        associate( &
            latitudes => sphere%grid%latitudes, &
            lambda => sphere%grid%longitudes(j) &
            )
            associate( &
                cos_a => cos(ALPHA), &
                sin_a => sin(ALPHA), &
                cos_l => cos(lambda), &
                sin_l => sin(lambda) &
                )
                do i = 1, NLATS

                    !     lambda is longitude, theta is colatitude, and pi/2-theta is
                    !     latitude on the rotated grid. lhat and that are longitude
                    !     and colatitude on the unrotated grid. see text starting at
                    !     equation (3.10)
                    !
                    associate( theta => HALF_PI-asin(latitudes(i)) )
                        associate( &
                            cos_t => cos(theta), &
                            sin_t => sin(theta) &
                            )
                            associate( &
                                sint   => cos_a*cos_t+sin_a*sin_t*cos_l, &
                                cthclh => cos_a*sin_t*cos_l-sin_a*cos_t, &
                                cthslh => sin_t*sin_l &
                                )
                                lambda_hat = atanxy(cthclh, cthslh)
                                associate( &
                                    cos_lh => cos(lambda_hat), &
                                    sin_lh => sin(lambda_hat) &
                                    )
                                    associate( cost => cos_lh*cthclh+sin_lh*cthslh )
                                        theta_hat = atanxy(sint, cost)
                                        u_hat = initialize_unrotated_longitudinal_velocity(U_ZERO, HALF_PI-theta_hat)
                                        pxact(j, i) = compute_cosine_transform(theta_hat, phlt)
                                        uxact(j, i) = u_hat*(cos_a*sin_l*sin_lh+cos_l*cos_lh)
                                        vxact(j, i) = u_hat*(cos_a*cos_l*sin_lh*cos_t-cos_lh*sin_l*cos_t+sin_a*sin_lh*sin_t)
                                        f(j, i) = F_ZERO * sint
                                        coslat(j, i) = sqrt(1.0_wp - latitudes(i)**2)
                                    end associate
                                end associate
                            end associate
                        end associate
                    end associate
                end do
            end associate
        end associate
    end do

    vmax = 0.0_wp
    pmax = 0.0_wp
    v2max = 0.0_wp
    p2max = 0.0_wp
    do j = 1, NLATS
        do i = 1, NLONS
            v2max = v2max + uxact(i, j)**2 + vxact(i, j)**2
            p2max = p2max + pxact(i, j)**2
            vmax = max(abs(uxact(i, j)), abs(vxact(i, j)), vmax)
            pmax = max(abs(pxact(i, j)), pmax)
        end do
    end do
    !
    !==> initialize first time step
    !
    u = uxact
    v = vxact
    p = pxact
    ug = u*coslat
    vg = v*coslat
    pg = p

    !
    !==>  compute spectral coeffs of initial vrt, div, p.
    !
    ! transpose data.
    ! minus sign to account for difference between
    ! mathematical and geophysical spherical coords.

    temp1 = transpose(ug)
    temp2 = -transpose(vg)
    call sphere%get_vorticity_and_divergence_coefficients_from_velocities( &
        temp1, temp2, vort_spec, div_spec)

    temp = transpose(p)
    call sphere%analyze_into_complex_spectral_coefficients( temp, p_spec)
    !
    !==> time step loop
    !
    n_new = 1
    n_now = 2
    n_old = 3

    !
    !==> Allocate memory for formatted model write statements
    !
    allocate( model_write_format, source=&
        '(A, i10, A, f10.2/, A, f10.0, A, i10/, A, i10, '&
        //'A, i10/, A, 1pe15.6, A, 1pe15.6, /A, 1pe15.6, A, 1pe15.6)' &
        )

    !
    !==> Allocate memory for formatted error write statements
    !
    allocate( error_write_format, source='(2(A, 1pe15.6)/, A, 1pe15.6)' )

    do cycle_number = 0, MAXIMUM_NUMBER_OF_TIME_ITERATIONS

        time = real(cycle_number, kind=wp)*DT

        !==> Inverse transform to get vort and phig on grid
        !
        call sphere%synthesize_from_complex_spectral_coefficients(vort_spec, temp )
        vrtg = transpose(temp)

        call sphere%synthesize_from_complex_spectral_coefficients(p_spec, temp)
        pg = transpose(temp)

        !
        !==> synthesize velocites u and v from spectral coeffs. of vort and div
        !
        call sphere%get_velocities_from_vorticity_and_divergence_coefficients( &
            vort_spec, div_spec, temp1, temp2 )
        ug = transpose(temp1)
        vg = -transpose(temp2)
        !
        !==> compute error statistics
        !
        if (mod(cycle_number, MPRINT) == 0) then
            ! ==> synthesize divergence from spectral coefficients
            call sphere%synthesize_from_complex_spectral_coefficients(div_spec, temp)
            divg = transpose(temp)
            u = ug
            v = vg
            p = pg
            model_time_in_hours = time/3600.0_wp

            !
            !==> Write model to console
            !
            write( stdout, '(A)' ) ''
            write( stdout, '(A)' ) ' steady nonlinear rotated flow:'
            write( stdout, fmt = model_write_format ) &
                ' cycle number              ', cycle_number, &
                ' model time in  hours      ', model_time_in_hours, &
                ' time step in seconds      ', DT, &
                ' number of latitudes       ', NLATS,    &
                ' number of longitudes      ', NLONS,    &
                ' max wave number           ', NTRUNC,    &
                ' rotation rate        ', ROTATION_RATE_OF_EARTH,   &
                ' mean height          ', P_ZERO,     &
                ' maximum velocity     ', U_ZERO,      &
                ' tilt angle           ', TILT_ANGLE

            !
            !==> Initialize equantities
            !
            dvgm = 0.0_wp
            dvmax = 0.0_wp
            dpmax = 0.0_wp
            evmax = 0.0_wp
            epmax = 0.0_wp

            do j=1, NLATS
                do i=1, NLONS
                    dvgm = &
                        max(dvgm, abs(divg(i, j)))
                    dvmax = &
                        dvmax+(u(i, j)-uxact(i, j))**2+(v(i, j)-vxact(i, j))**2
                    dpmax = &
                        dpmax+(p(i, j)-pxact(i, j))**2
                    evmax = &
                        max(evmax, abs(v(i, j)-vxact(i, j)), abs(u(i, j)-uxact(i, j)))
                    epmax = &
                        max(epmax, abs(p(i, j)-pxact(i, j)))
                end do
            end do

            dvmax = sqrt(dvmax/v2max)
            dpmax = sqrt(dpmax/p2max)
            evmax = evmax/vmax
            epmax = epmax/pmax

            !
            !===> Write errors to console
            !
            write( stdout, fmt = error_write_format ) &
                ' max error in velocity', evmax, &
                ' max error in geopot. ', epmax, &
                ' l2 error in velocity ', dvmax, &
                ' l2 error in geopot.  ', dpmax, &
                ' maximum divergence   ', dvgm
        end if

        !==> Compute right-hand sides of prognostic eqns.

        scrg1 = ug * ( vrtg + f )
        scrg2 = vg * ( vrtg + f )

        temp1 = transpose(scrg1)
        temp2 = -transpose(scrg2 )
        call sphere%get_vorticity_and_divergence_coefficients_from_velocities( &
            temp1, temp2, ddivdt_spec(:, n_new), dvrtdt_spec(:, n_new) )

        dvrtdt_spec(:, n_new) = -dvrtdt_spec(:, n_new)

        scrg1 = ug * ( pg + P_ZERO )
        scrg2 = vg * ( pg + P_ZERO )

        temp1 = transpose(scrg1)
        temp2 = -transpose(scrg2)
        call sphere%get_vorticity_and_divergence_coefficients_from_velocities( &
            temp1, temp2, scr_spec, dpdt_spec(:, n_new))

        dpdt_spec(:, n_new) = -dpdt_spec(:, n_new)

        scrg1 = pg + 0.5_wp * ( (ug**2 + vg**2) )!/ coslat**2 )

        temp1 = transpose(scrg1 )
        call sphere%analyze_into_complex_spectral_coefficients( temp1, scr_spec)

        associate( lap => sphere%laplacian_coefficients )
            ddivdt_spec(:, n_new) = ddivdt_spec(:, n_new) - lap * scr_spec
        end associate

        !==> update vrt and div with third-order adams-bashforth.

        !==> forward euler, then 2nd-order adams-bashforth time steps to start.

        select case (cycle_number)
            case (0)
                dvrtdt_spec(:, n_now) = dvrtdt_spec(:, n_new)
                dvrtdt_spec(:, n_old) = dvrtdt_spec(:, n_new)
                ddivdt_spec(:, n_now) = ddivdt_spec(:, n_new)
                ddivdt_spec(:, n_old) = ddivdt_spec(:, n_new)
                dpdt_spec(:, n_now) = dpdt_spec(:, n_new)
                dpdt_spec(:, n_old) = dpdt_spec(:, n_new)
            case (1)
                dvrtdt_spec(:, n_old) = dvrtdt_spec(:, n_new)
                ddivdt_spec(:, n_old) = ddivdt_spec(:, n_new)
                dpdt_spec(:, n_old) = dpdt_spec(:, n_new)
        end select

        vort_spec = &
            vort_spec + DT * (&
            (23.0_wp/12.0_wp) * dvrtdt_spec(:, n_new) &
            - (16.0_wp/12.0_wp) * dvrtdt_spec(:, n_now) &
            + (5.0_wp/12.0_wp) * dvrtdt_spec(:, n_old) )

        div_spec = &
            div_spec + DT *( &
            (23.0_wp/12.0_wp) * ddivdt_spec(:, n_new) &
            - (16.0_wp/12.0_wp) * ddivdt_spec(:, n_now) &
            + (5.0_wp/12.0_wp) * ddivdt_spec(:, n_old) )

        p_spec = &
            p_spec + DT * (&
            (23.0_wp/12.0_wp) * dpdt_spec(:, n_new) &
            - (16.0_wp/12.0_wp) * dpdt_spec(:, n_now) &
            + (5.0_wp/12.0_wp) * dpdt_spec(:, n_old) )

        !==> switch indices

        n_save1 = n_new
        n_save2 = n_now
        n_new = n_old
        n_now = n_save1
        n_old = n_save2

    !==> end time step loop
    end do

    !
    !==> Print compiler info
    !
    write( stdout, '(A)' ) ''
    write( stdout, '(4A)' ) 'This file was compiled by ', &
        compiler_version(), ' using the options ', &
        compiler_options()
    write( stdout, '(A)' ) ''

    !
    !==>  Release memory
    !
    call sphere%destroy()
    deallocate( model_write_format )
    deallocate( error_write_format )

contains

    pure function initialize_unrotated_longitudinal_velocity( &
        amp, thetad ) result (return_value)
        !
        !     computes the initial unrotated longitudinal velocity
        !     see section 3.3.
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in) :: amp
        real (wp), intent (in) :: thetad
        real (wp)              :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (wp), parameter :: ZERO = nearest(1.0_wp, 1.0_wp)-nearest(1.0_wp, -1.0_wp)
        real (wp), parameter :: PI = acos( -1.0_wp )
        real (wp)            :: x
        !--------------------------------------------------------------------------------

        associate( &
            thetab => -PI/6.0_wp, &
            thetae => PI/2.0_wp, &
            xe => 3.0e-1_wp &
            )

            x =xe*(thetad-thetab)/(thetae-thetab)

            return_value = 0.0_wp

            if(x <= ZERO .or. x >= xe) return

            associate( arg => (-1.0_wp/x) - (1.0_wp/(xe-x)) + (4.0_wp/xe) )
                return_value = amp * exp( arg )
            end associate
        end associate

    end function initialize_unrotated_longitudinal_velocity


    pure function atanxy( x, y ) result (return_value)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in) :: x
        real (wp), intent (in) :: y
        real (wp)              :: return_value
        !--------------------------------------------------------------------------------
        real (wp), parameter :: ZERO = nearest(1.0_wp, 1.0_wp)-nearest(1.0_wp, -1.0_wp)
        !--------------------------------------------------------------------------------

        ! Initialize return value
        return_value = 0.0_wp

        if ( x == ZERO .and. y == ZERO ) return

        return_value = atan2( y, x )

    end function atanxy


    subroutine compute_sine_transform( x )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in out) :: x(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)           :: i, j !! Counters
        real (wp), allocatable ::  w(:)
        !--------------------------------------------------------------------------------

        associate( n => size(x) )
            ! Allocate memory
            allocate( w(n) )
            ! Associate various quantities
            associate( arg => acos(-1.0_wp)/(n+1) )
                do j = 1, n
                    w(j) = 0.0_wp
                    do i = 1, n
                        associate( sin_arg => real(i*j, kind=wp)*arg )
                            w(j) = w(j)+x(i)*sin(sin_arg)
                        end associate
                    end do
                end do
            end associate
        end associate

        x = 2.0_wp * w

        ! Release memory
        deallocate(w)

    end subroutine compute_sine_transform


    pure function compute_cosine_transform(theta, cf) result (return_value)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent(in) :: theta
        real (wp), intent(in) :: cf(:)
        real (wp)             :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)          :: i !! Counter
        !--------------------------------------------------------------------------------

        return_value = 0.0_wp

        associate( n => size(cf) )
            do i=1, n
                return_value = return_value + cf(i)*cos(i*theta)
            end do
        end associate

    end function compute_cosine_transform


end program shallow
