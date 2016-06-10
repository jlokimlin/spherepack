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
! ... file shallow.f90
!
!     program shallow solves the nonlinear shallow-water equations
!     on the sphere using spherepack software.
!
! ... required spherepack files
!
!     vtses.f90, dives.f90, vrtes.f90, grades.f90, type_SpherepackAux.f90, type_HFFTpack.f90,
!     vhaes.f90,vhses.f90,shaes.f90,shses.f90
!
!
!     the nonlinear shallow-water equations on the sphere are
!     solved using a spectral method based on the spherical
!     vector harmonics. the method is described in the paper:
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
!     nlat          number of latitudes including poles
!     nlon          number of distinct longitudes
!     mmode         max wave number
!     omega         rotation rate of earth in radians per second
!     aa            radius of earth in meters
!     pzero         mean height of geopotential
!     uzero         maximum velocity
!     alpha         tilt angle of the rotated grid
!     ncycle        exit number
!     time          model time in seconds
!     dt            time step
!     lambda        longitude
!     theta         colatitude
!
!   the first dimension of the following two dimensional arrays
!   corresponds to the latitude index with values i=1,...,nlat
!   where i=1 is the north pole and i=nlat is the south pole.
!   the second dimension is longitude with values j=1,...,nlon
!   where j=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!     u(i,j)       east longitudinal velocity component at t=time
!     v(i,j)       latitudinal velocity component at t=time
!     p(i,j)       +pzero = geopotential at t=time
!
!     unew(i,j)    east longitudinal velocity component at t=time+dt
!     vnew(i,j)    latitudinal velocity component at t=time+dt
!     pnew(i,j)    +pzero = geopotential at t=time+dt
!
!     uold(i,j)    east longitudinal velocity component at t=time-dt
!     vold(i,j)    latitudinal velocity component at t=time-dt
!     pold(i,j)    +pzero = geopotential at t=time-dt
!
!     divg(i,j)    divergence (d/dtheta (cos(theta) v)
!                                          + du/dlambda)/cos(theta)
!     vort(i,j)    vorticity  (d/dtheta (cos(theta) u)
!                                          - dv/dlambda)/cos(theta)
!
!     ut(i,j)      latitudinal derivative of longitudinal
!                  velocity component
!     vt(i,j)      latitudinal derivative of latitudinal
!                  velocity component
!
!     dudt(i,j)    time derivative of longitudinal velocity component
!     dvdt(i,j)    time derivative of latitudinal  velocity component
!     dpdt(i,j)    time derivative of geopotential
!
!     gpdl(i,j)    first component of the gradient of p(i,j)
!                  the longitudinal derivative of the geopotential
!                  divided by the cosine of the latitude
!
!     gpdt(i,j)    second component of the gradient of p(i,j)
!                  the latitudinal derivative of the geopotential
!
!     uxact(i,j)   the "exact" longitudinal veloctiy component
!     vxact(i,j)   the "exact" latitudinal  veloctiy component
!     uxact(i,j)   the "exact" geopotential
!
!     f(i,j)       the coriolis force on rotated grid
!
!   the following two dimensional arrays are nonzero in the triangle
!   n=1,...,nlat and m less than or equal to n.
!
!     a(m,n),b(m,n)    spectral coefficients of the geopotential
!
!     br(m,n),bi(m,n)  spectral coefficients of the velocity
!     cr(m,n),ci(m,n)  vector [u(i,j),v(i,j)]
!
!
!     phlt(i)      the coefficients in the cosine series
!                  representation of the unrotated geopotential
!


module type_ShallowWaterSolver

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        Regularsphere

    ! Explicit typing only
    implicit none

    ! Declare derived data type
    type, extends (RegularSphere) :: ShallowWaterSolver
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure         :: spec_to_grid
        procedure         :: grid_to_spec
        procedure         :: get_vrtdivspec
        procedure         :: get_uv
        procedure, nopass :: atanxy
        procedure, nopass :: get_initial_velocity
        procedure, nopass :: sine_transform
        procedure, nopass :: cosine_transform
        !----------------------------------------------------------------------
    end type ShallowWaterSolver

contains


    function get_initial_velocity(amp, thetad) result (return_value)
        implicit none
        !
        !     computes the initial unrotated longitudinal velocity
        !     see section 3.03.0
        !
        real (wp) :: amp
        real (wp) :: thetad
        real (wp) :: return_value
        !------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------
        real (wp), parameter :: pi = acos(-1.0_wp)
        real (wp)            :: thetab, thetae, xe, x
        !------------------------------------------------------

        thetab = -pi/6
        thetae = pi/2
        xe = 3.0e-1_wp

        x = xe*(thetad-thetab)/(thetae-thetab)

        return_value = 0.0_wp

        if (x <= 0.0_wp .or. x >= xe) return

        return_value = amp*exp(-1.0_wp /x-1.0_wp/(xe-x)+4.0_wp/xe)

    end function get_initial_velocity

    pure function atanxy(x, y) result (return_value)

        real (wp), intent (in) :: x
        real (wp), intent (in) :: y
        real (wp)              :: return_value

        return_value = 0.0_wp

        if (x == 0.0_wp .and. y == 0.0_wp) return

        return_value = atan2(y, x)

    end function atanxy


    subroutine sine_transform(n, x)

        !     computes the sine transform
        integer (ip), intent (in)     :: n
        real (wp),    intent (in out) :: x(n)

        integer (ip) ::  i, j
        real (wp)    :: w(n), arg

        arg = acos(-1.0_wp)/(n+1)

        do j=1, n
            w(j) = 0.0_wp
            do i=1, n
                w(j) = w(j)+x(i)*sin(real(i*j, kind=wp)*arg)
            end do
        end do

        do i=1, n
            x(i) = 2.0_wp *w(i)
        end do

    end subroutine sine_transform


    pure function cosine_transform(theta, n, cf) result (return_value)

        !     computes the cosine transform
        real (wp),    intent (in) :: theta
        integer (ip), intent (in) :: n
        real (wp),    intent (in) :: cf(n)
        real (wp)                 :: return_value

        integer (ip) ::  i

        return_value = 0.0_wp

        do i=1, n
            return_value = return_value+cf(i)*cos(real(i, kind=wp)*theta)
        end do

    end function cosine_transform



    subroutine grid_to_spec(this, datagrid, dataspec)

        ! converts complex spectral coefficients (dataspec) to
        ! gridded data array (datagrid).

        class (ShallowWaterSolver), intent (in out) :: this
        real (wp),                  intent (in)     :: datagrid(:,:)
        complex (wp),               intent (out)    :: dataspec(:)

        integer (ip)           :: n, m
        real (wp), allocatable :: temp(:,:)


        associate( &
            nlat => size(datagrid, dim=2),&
            nlon => size(datagrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( temp(nlat,nlon) )

        end associate

        !
        !==> Transpose data
        !
        temp = transpose(datagrid)

        !
        !==> spherical harmonic analysis
        !
        call this%perform_scalar_analysis(temp)

        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients &
            )

            !
            !==> Fill complex array dataspec with result.
            !
            dataspec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp)
            !
            !==> Reset constants
            !
            a = 0.0_wp
            b = 0.0_wp
        end associate
        !
        !==> Release memory
        !
        deallocate( temp )

    end subroutine grid_to_spec


    subroutine spec_to_grid(this, dataspec, datagrid)

        ! converts complex spectral coefficients (dataspec) to
        ! gridded data array (datagrid).

        class (ShallowWaterSolver), intent (in out) :: this
        real (wp),                  intent (in out) :: datagrid(:,:)
        complex (wp),               intent (in)     :: dataspec(:)

        integer (ip)           :: n, m, nm, i
        real (wp), allocatable :: temp(:,:)


        associate( &
            nlat => size(datagrid, dim=2),&
            nlon => size(datagrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( temp(nlat,nlon) )

        end associate


        !
        !==> fill two real arrays (a, b) with contents of dataspec.
        !
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients &
            )

            a = 0.0_wp
            b = 0.0_wp
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = 2.0_wp * real(dataspec(nm))
                    b(m, n) = 2.0_wp * aimag(dataspec(nm))
                end do
            end do

            !
            !==> Perform spherical harmonic synthesis
            !
            call this%perform_scalar_synthesis(temp)

            !
            !==> Reset coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
        end associate
        !
        !==> Transpose data
        !
        datagrid = transpose(temp)

        !
        !==> Release memory
        !
        deallocate( temp )

    end subroutine spec_to_grid



    subroutine get_vrtdivspec(this, ugrid, vgrid, vrtspec, divspec)
        !
        ! Purpose:
        !
        ! Calculate spectral coefficients of vorticity and divergence
        ! (vrtspec, divspec) given input gridded winds (ugrid, vgrid).
        !
        class (ShallowWaterSolver), intent (in out) :: this
        real (wp),                  intent (in)     :: ugrid(:,:)
        real (wp),                  intent (in)     :: vgrid(:,:)
        complex (wp),               intent (out)    :: vrtspec(:)
        complex (wp),               intent (out)    :: divspec(:)

        real (wp)              :: fn
        real (wp), allocatable :: v(:,:), w(:,:), sqnn(:)
        integer (ip)           :: n, m !! Counters

        associate( &
            nlat => size(ugrid, dim=2),&
            nlon => size(ugrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( v(nlat,nlon) )
            allocate( w(nlat,nlon) )
            allocate( sqnn(nlat) )

        end associate

        !
        !==> Transpose data.
        !    minus sign to account for difference between
        !    mathematical and geophysical spherical coords
        !
        v = -transpose(vgrid)
        w = transpose(ugrid)

        !
        !==> Calculate vector spherical harmonic analysis.
        !
        call this%vector_analysis_from_spherical_components(v, w)

        !
        !==> Multiply vector harmonic coefficients of winds by
        !    appropriate factors to convert into vorticity and
        !    divergence coefficients.
        !
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            nlat => this%NUMBER_OF_LATITUDES, &
            rsphere => this%RADIUS_OF_SPHERE, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )

            do n=1, nlat
                fn = real(n - 1, kind=wp)
                sqnn(n) = sqrt(fn * (fn + 1.0_wp))
            end do

            a = 0.0_wp
            b = 0.0_wp
            do n=1, nlat
                a(:,n) = -(sqnn(n)/rsphere)*br(:,n)
                b(:,n) = -(sqnn(n)/rsphere)*bi(:,n)
            end do

            divspec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp)

            a = 0.0_wp
            b = 0.0_wp
            do n=1, nlat
                a(:,n) = (sqnn(n)/rsphere)*cr(:,n)
                b(:,n) = (sqnn(n)/rsphere)*ci(:,n)
            end do

            vrtspec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp)

            !
            !==> Reset coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
            br = 0.0_wp
            bi = 0.0_wp
            cr = 0.0_wp
            ci = 0.0_wp
        end associate

        !
        !==> Release memory
        !
        deallocate( v )
        deallocate( w )
        deallocate( sqnn )

    end subroutine get_vrtdivspec


    subroutine get_uv(this, vrtspec, divspec, ugrid, vgrid)
        !
        ! Purpose:
        !
        ! Given spectral coefficients of vorticity and divergence
        ! (vrtspec, divspec) compute gridded winds (ugrid, vgrid)
        !

        class (ShallowWaterSolver), intent (in out) :: this
        complex (wp),               intent (in)     :: vrtspec(:)
        complex (wp),               intent (in)     :: divspec(:)
        real (wp),                  intent (out)    :: ugrid(:,:)
        real (wp),                  intent (out)    :: vgrid(:,:)

        real (wp)              :: fn
        real (wp), allocatable :: v(:,:), w(:,:), isqnn(:)
        integer (ip)           :: n, m, nm, i !! Counters

        associate( &
            nlat => size(ugrid, dim=2),&
            nlon => size(ugrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( v(nlat,nlon) )
            allocate( w(nlat,nlon) )
            allocate( isqnn(nlat) )

        end associate

        !
        !==> multiply spectral coefficients of vorticity and divergence
        !    by appropriate factors to convert them into vector harmonic
        !    coefficients of winds.
        !
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            nlat => this%NUMBER_OF_LATITUDES, &
            rsphere => this%RADIUS_OF_SPHERE, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )

            isqnn(1) = 0.0_wp
            do n=2, nlat
                fn = real(n - 1, wp)
                isqnn(n) = rsphere/sqrt(fn*(fn+1.0_wp))
            end do

            a = 0.0_wp
            b = 0.0_wp
            br = 0.0_wp
            bi = 0.0_wp

            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = -2.0_wp * real(divspec(nm))
                    b(m, n) = -2.0_wp * aimag(divspec(nm))
                end do
            end do

            do n=1, nlat
                br(:,n) = isqnn(n)*a(:,n)
                bi(:,n) = isqnn(n)*b(:,n)
            end do

            cr = 0.0_wp
            ci = 0.0_wp
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = 2.0_wp * real(vrtspec(nm))
                    b(m, n) = 2.0_wp * aimag(vrtspec(nm))
                end do
            end do

            do n=1, nlat
                cr(:,n) = isqnn(n)*a(:,n)
                ci(:,n) = isqnn(n)*b(:,n)
            end do



            !
            !==> compute vector harmonic synthesis to get winds on grid.
            !
            call this%perform_vector_synthesis(v, w)

            !
            !==> Reset coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
            br = 0.0_wp
            bi = 0.0_wp
            cr = 0.0_wp
            ci = 0.0_wp

        end associate
        !
        !==> Transpose data
        !    minus sign to account for differences
        !    between mathematical and geophysical spherical coords.
        !
        vgrid = -transpose(v)
        ugrid = transpose(w)

        !
        !==> Release memory
        !
        deallocate( v )
        deallocate( w )
        deallocate( isqnn )

    end subroutine get_uv

end module type_ShallowWaterSolver


program shallow

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use type_ShallowWaterSolver

    ! Explicit typing only
    implicit none

    !--------------------------------------------------------------------------------
    ! Dictionary
    !--------------------------------------------------------------------------------
    integer (ip), parameter           :: nlon = 128
    integer (ip), parameter           :: nlat = nlon/2+1
    integer (ip), parameter           :: ntrunc = 42
    integer (ip), parameter           :: nmdim = (ntrunc+1)*(ntrunc+2)/2
    integer (ip), parameter           :: nl = 91
    integer (ip), parameter           :: nlm1 = nl-1
    integer (ip), parameter           :: nlm2 = nl-2
    integer (ip)                      :: itmax, mprint, i, j, ncycle
    integer (ip)                      :: temp_save_new, temp_save_now,old,now,new
    real (wp), dimension(nlon, nlat)  :: uxact, vxact, pxact, u, v, p, f
    real (wp), dimension(nlon, nlat)  :: ug, vg, pg, vrtg, divg, scrg1, scrg2
    real (wp)                         :: phlt(nlm2)
    real (wp)                         :: theta_mesh, lambda, lhat, uhat, aa, uzero, pzero, pi, half_pi, radian_unit, omega, alphad
    real (wp)                         :: alpha, fzero, dt, cfn, dlath, theta, sth
    real (wp)                         :: cth, cosa, sina, lambda_mesh, st, ct, cthclh, cthslh
    real (wp)                         :: clh, time, that, sl, slh, evmax, epmax, dvmax, dpmax, htime, dvgm, cl
    real (wp)                         :: v2max, p2max, vmax, pmax
    complex (wp), dimension(nmdim)    :: vrtnm, divnm, pnm, scrnm
    complex (wp), dimension(nmdim, 3) :: dvrtdtnm, ddivdtnm, dpdtnm
    type (ShallowWaterSolver)         :: solver
    character (len=:), allocatable    :: write_fmt
    !--------------------------------------------------------------------------------

    write( stdout, '(/A)') '     shallow *** TEST RUN *** '

    !
    !==> Initialize constants
    !
    pi = acos(-1.0_wp)
    half_pi = pi/2
    radian_unit = pi/180.0_wp
    aa = 6.37122e+6_wp
    omega = 7.292e-5_wp
    fzero = 2.0_wp * omega
    uzero = 40.0_wp
    pzero = 2.94e+4_wp
    alphad = 60.0_wp
    alpha = radian_unit * alphad
    dt = 300.0_wp
    itmax = nint(86400.0_wp * 5.0_wp/dt, kind=ip)
    mprint = itmax/10

    !
    !==> Allocate memory
    !
    call solver%create(nlat=nlat, nlon=nlon, ntrunc=ntrunc, rsphere=aa)

    allocate( write_fmt, source=&
        '(A, i10, A, f10.2/, A, f10.0, A, i10/, A, i10, '&
        //'A, i10/, A, 1pe15.6, A, 1pe15.6, /A, 1pe15.6, A, 1pe15.6)' )

    !
    !==> compute the derivative of the unrotated geopotential
    !    p as a function of latitude
    !
    cfn = 1.0_wp / nlm1
    dlath = pi / nlm1
    do i=1, nlm2
        theta = i*dlath
        sth = sin(theta)
        cth = cos(theta)
        uhat = solver%get_initial_velocity(uzero, half_pi-theta)
        phlt(i) = cfn*cth*uhat*(uhat/sth+aa*fzero)
    end do
    !
    !     compute sine transform of the derivative of the geopotential
    !     for the purpose of computing the geopotential by integration
    !     see equation (3.9) in reference [1] above
    !
    call solver%sine_transform(nlm2, phlt)
       !
       !     compute the cosine coefficients of the unrotated geopotential
       !     by the formal integration of the sine series representation
       !
    do i=1, nlm2
        phlt(i) = -phlt(i)/i
    end do
    !
    !     phlt(i) contains the coefficients in the cosine series
    !     representation of the unrotated geopotential that are used
    !     below to compute the geopotential on the rotated grid.0
    !
    !     compute the initial values of  east longitudinal
    !     and latitudinal velocities u and v as well as the
    !     geopotential p and coriolis f on the rotated grid.0
    !
    cosa = cos(alpha)
    sina = sin(alpha)
    lambda_mesh = (2.0_wp*pi)/nlon
    theta_mesh = pi/(nlat-1)

    do j=1, nlon
        lambda = real(j - 1, wp) * lambda_mesh
        cl = cos(lambda)
        sl = sin(lambda)
        do i=1, nlat
            !
            !     lambda is longitude, theta is colatitude, and pi/2-theta is
            !     latitude on the rotated grid.0 lhat and that are longitude
            !     and colatitude on the unrotated grid.0 see text starting at
            !     equation (3.010)
            !
            theta = real(i - 1, wp)*theta_mesh
            st = cos(theta)
            ct = sin(theta)
            sth = cosa*st+sina*ct*cl
            cthclh = cosa*ct*cl-sina*st
            cthslh = ct*sl
            lhat = atanxy(cthclh, cthslh)
            clh = cos(lhat)
            slh = sin(lhat)
            cth = clh*cthclh+slh*cthslh
            that = solver%atanxy(sth, cth)
            uhat = solver%get_initial_velocity(uzero, half_pi-that)
            pxact(j, i) = solver%cosine_transform(that, nlm2, phlt)
            uxact(j, i) = uhat*(cosa*sl*slh+cl*clh)
            vxact(j, i) = uhat*(cosa*cl*slh*st-clh*sl*st+sina*slh*ct)
            f(j, i) = fzero*sth
        end do
    end do

    vmax = 0.0_wp
    pmax = 0.0_wp
    v2max = 0.0_wp
    p2max = 0.0_wp

    do j=1, nlat
        do i=1, nlon
            v2max = v2max+uxact(i, j)**2+vxact(i, j)**2
            p2max = p2max+pxact(i, j)**2
            vmax = max(abs(uxact(i, j)), abs(vxact(i, j)), vmax)
            pmax = max(abs(pxact(i, j)), pmax)
        end do
    end do
    !
    !     initialize first time step
    !
    u = uxact
    v = vxact
    p = pxact
    ug = u
    vg = v
    pg = p

    !
    !==> Compute spectral coeffs of initial vrt, div, p
    !
    call solver%get_vrtdivspec(ug, vg, vrtnm, divnm)
    call solver%grid_to_spec(pg, pnm)

    !
    !==> time step loop
    !
    new = 1
    now = 2
    old = 3
    do ncycle = 0, itmax

        time = real(ncycle, kind=wp) * dt

        !
        !==> Inverse transform to get vort and phig on grid
        !
        call solver%spec_to_grid(vrtnm, vrtg)
        call solver%spec_to_grid(pnm, pg)

        !
        !==> compute u and v on grid from spectral coeffs of vort and div.
        !
        call solver%get_uv(vrtnm, divnm, ug, vg)

        !==> compute error statistics

        if (mod(ncycle, mprint) == 0) then

            call solver%spec_to_grid(divnm, divg)

            u = ug
            v = vg
            p = pg
            htime = time/3600.0_wp

            write( stdout, '(/A)' ) ' steady nonlinear rotated flow:'
            write( stdout, fmt=write_fmt ) &
                ' cycle number              ', ncycle, &
                ' model time in  hours      ', htime, &
                ' time step in seconds      ', dt, &
                ' number of latitudes       ', nlat, &
                ' number of longitudes      ', nlon, &
                ' max wave number           ', ntrunc, &
                ' rotation rate        '     , omega, &
                ' mean height          '     , pzero, &
                ' maximum velocity     '     , uzero, &
                ' tilt angle           '     , alphad

            dvgm = 0.0_wp
            dvmax = 0.0_wp
            dpmax = 0.0_wp
            evmax = 0.0_wp
            epmax = 0.0_wp

            do j=1, nlat
                do i=1, nlon
                    dvgm = max(dvgm, abs(divg(i, j)))
                    dvmax = dvmax+(u(i, j)-uxact(i, j))**2+(v(i, j)-vxact(i, j))**2
                    dpmax = dpmax+(p(i, j)-pxact(i, j))**2
                    evmax = &
                        max(evmax, abs(v(i, j)-vxact(i, j)), abs(u(i, j)-uxact(i, j)))
                    epmax = max(epmax, abs(p(i, j)-pxact(i, j)))
                end do
            end do

            dvmax = sqrt(dvmax/v2max)
            dpmax = sqrt(dpmax/p2max)
            evmax = evmax/vmax
            epmax = epmax/pmax

            write( stdout, fmt='(2(A, 1pe15.6)/, A, 1pe15.6)') &
                ' max error in velocity', evmax, &
                ' max error in geopot. ', epmax, &
                ' l2 error in velocity ', dvmax, &
                ' l2 error in geopot.  ', dpmax, &
                ' maximum divergence   ', dvgm

        end if
        !
        !==> Compute right-hand sides of prognostic eqns
        !
        scrg1 = ug * (vrtg + f)
        scrg2 = vg * (vrtg + f)

        call solver%get_vrtdivspec(scrg1, scrg2, ddivdtnm(:,new), dvrtdtnm(:,new))

        dvrtdtnm(:,new) = -dvrtdtnm(:,new)
        scrg1 = ug*(pg+pzero)
        scrg2 = vg*(pg+pzero)

        call solver%get_vrtdivspec(scrg1, scrg2, scrnm, dpdtnm(:,new))

        dpdtnm(:,new) = -dpdtnm(:,new)
        scrg1 = pg + 0.5_wp * (ug**2+vg**2)

        call solver%grid_to_spec(scrg1, scrnm)


        associate( lap => solver%laplacian_coefficients )

            ddivdtnm(:,new) = ddivdtnm(:,new)-lap*scrnm

        end associate
        !
        !==> Update vrt and div with third-order adams-bashforth
        !
        select case (ncycle)
            !
            !==> forward euler, then 2nd-order adams-bashforth time steps to start
            !
            case (0)
                dvrtdtnm(:,now) = dvrtdtnm(:,new)
                dvrtdtnm(:,old) = dvrtdtnm(:,new)
                ddivdtnm(:,now) = ddivdtnm(:,new)
                ddivdtnm(:,old) = ddivdtnm(:,new)
                dpdtnm(:,now) = dpdtnm(:,new)
                dpdtnm(:,old) = dpdtnm(:,new)
            case (1)
                dvrtdtnm(:,old) = dvrtdtnm(:,new)
                ddivdtnm(:,old) = ddivdtnm(:,new)
                dpdtnm(:,old) = dpdtnm(:,new)
        end select

        vrtnm = vrtnm + dt*( &
            (23.0_wp/12)*dvrtdtnm(:,new) &
            - (16.0_wp/12)*dvrtdtnm(:,now) &
            + (5.0_wp/12)*dvrtdtnm(:,old) )

        divnm = divnm + dt*( &
            (23.0_wp/12)*ddivdtnm(:,new) &
            - (16.0_wp/12)*ddivdtnm(:,now) &
            + (5.0_wp/12)*ddivdtnm(:,old) )

        pnm = pnm + dt*( &
            (23.0_wp/12)*dpdtnm(:,new) &
            - (16.0_wp/12)*dpdtnm(:,now) &
            + (5.0_wp/12)*dpdtnm(:,old) )

        !
        !==> Switch indices
        !
        temp_save_new = new
        temp_save_now = now
        new = old
        now = temp_save_new
        old = temp_save_now

    end do

    !
    !==> Release memory
    !
    call solver%destroy()
    deallocate( write_fmt )

end program shallow
