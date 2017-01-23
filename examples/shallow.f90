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
!     vtses.f90, dives.f90, vrtes.f90, grades.f90, type_SpherepackAux.f90, type_RealPeriodicTransform.f90,
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
program shallow

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        TWO_PI

    use type_ShallowWaterSolver, only: &
        ShallowwaterSolver

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter           :: NLON = 128 ! Number of longitudinal points
    integer(ip), parameter           :: NLAT = NLON/2+1 ! Number of latitudinal points
    integer(ip), parameter           :: NTRUNC = 42
    integer(ip), parameter           :: NMDIM = (NTRUNC+1)*(NTRUNC+2)/2
    integer(ip), parameter           :: nl = 91
    integer(ip)                      :: max_iter, mprint, i, j, ncycle
    integer(ip)                      :: temp_save_new, temp_save_now,old,now,new
    real(wp), dimension(NLON, NLAT)  :: uxact, vxact, pxact, u, v, p, f
    real(wp), dimension(NLON, NLAT)  :: ug, vg, pg, vrtg, divg, scrg1, scrg2
    real(wp)                         :: phlt(nl-2)
    real(wp)                         :: lhat, uhat, aa, uzero, pzero, HALF_PI, radian_unit, omega, alpha_in_degrees
    real(wp)                         :: alpha, fzero, dt, theta, sth
    real(wp)                         :: cth, cthclh, cthslh
    real(wp)                         :: clh, time, that, slh, htime
    real(wp)                         :: v2max, p2max, vmax, pmax
    real(wp), parameter              :: ZERO = 0.0_wp, HALF = 0.5_wp, ONE = 1.0_wp
    complex(wp), dimension(NMDIM)    :: vrtnm, divnm, pnm, scrnm
    complex(wp), dimension(NMDIM, 3) :: dvrtdtnm, ddivdtnm, dpdtnm
    type(ShallowWaterSolver)         :: solver
    character(len=:), allocatable    :: write_fmt

    write( stdout, '(/a)') '     shallow *** TEST RUN *** '

    !  Parameters for test
    HALF_PI = PI/2
    radian_unit = PI/180
    aa = 6.37122e+6_wp ! Earth radius
    omega = 7.292e-5_wp ! Rotation rate
    fzero = 2.0_wp * omega
    uzero = 40.0_wp
    pzero = 2.94e+4_wp
    alpha_in_degrees = 60.0_wp
    alpha = radian_unit * alpha_in_degrees
    dt = 300.0_wp ! Time step in seconds
    max_iter = nint(86400.0_wp * HALF/dt, kind=ip) ! Integration length in days
    mprint = max_iter/10

    !  Allocate memory. Setup spherical harmonic instance
    call solver%create(nlat=NLAT, nlon=NLON, ntrunc=NTRUNC, rsphere=aa)

    allocate( write_fmt, source=&
        '(a, i10, a, f10.2/, a, f10.0, a, i10/, a, i10, '&
        //'a, i10/, a, 1pe15.6, a, 1pe15.6, /a, 1pe15.6, a, 1pe15.6)' )

    ! Compute the derivative of the unrotated geopotential
    ! p as a function of latitude
    block
        real(wp) :: cfn, dlath

        cfn = ONE / (nl - 1)
        dlath = PI / (nl - 1)
        do i=1, size(phlt)
            theta = real(i, kind=wp) * dlath
            uhat = solver%get_initial_velocity(uzero, HALF_PI-theta)
            phlt(i) = cfn * cos(theta) * uhat * (uhat/sin(theta) + aa*fzero)
        end do
    end block

    ! Compute sine transform of the derivative of the geopotential
    ! for the purpose of computing the geopotential by integration
    ! see equation (3.9) in reference [1] above
    call solver%sine_transform(phlt)

    ! Compute the cosine coefficients of the unrotated geopotential
    ! by the formal integration of the sine series representation
    do i=1, size(phlt)
        phlt(i) = -phlt(i)/i
    end do

    !     phlt(i) contains the coefficients in the cosine series
    !     representation of the unrotated geopotential that are used
    !     below to compute the geopotential on the rotated grid.
    !
    !     compute the initial values of  east longitudinal
    !     and latitudinal velocities u and v as well as the
    !     geopotential p and coriolis f on the rotated grid.
    block
        real(wp) :: cosa, sina, lambda_mesh, theta_mesh
        real(wp) :: lambda, cosl, sinl
        real(wp) :: sint, cost

        cosa = cos(alpha)
        sina = sin(alpha)
        lambda_mesh = TWO_PI/NLON
        theta_mesh = PI/(NLAT-1)

        do j=1, NLON
            lambda = real(j - 1, kind=wp) * lambda_mesh
            cosl = cos(lambda)
            sinl = sin(lambda)
            do i=1, NLAT
                !
                !     lambda is longitude, theta is colatitude, and pi/2-theta is
                !     latitude on the rotated grid. lhat and that are longitude
                !     and colatitude on the unrotated grid. see text starting at
                !     equation (3.10)
                !
                theta = real(i - 1, kind=wp) * theta_mesh
                sint = cos(theta)
                cost = sin(theta)
                sth = cosa*sint+sina*cost*cosl
                cthclh = cosa*cost*cosl-sina*sint
                cthslh = cost*sinl
                lhat = solver%atanxy(cthclh, cthslh)
                clh = cos(lhat)
                slh = sin(lhat)
                cth = clh*cthclh+slh*cthslh
                that = solver%atanxy(sth, cth)
                uhat = solver%get_initial_velocity(uzero, HALF_PI-that)
                pxact(j, i) = solver%cosine_transform(that, phlt)
                uxact(j, i) = uhat*(cosa*sinl*slh+cosl*clh)
                vxact(j, i) = uhat*(cosa*cosl*slh*sint-clh*sinl*sint+sina*slh*cost)
                f(j, i) = fzero*sth ! Coriolis
            end do
        end do
    end block

    vmax = ZERO
    pmax = ZERO
    v2max = ZERO
    p2max = ZERO

    do j=1, NLAT
        do i=1, NLON
            v2max = v2max + uxact(i, j)**2 + vxact(i, j)**2
            p2max = p2max + pxact(i, j)**2
            vmax = max(abs(uxact(i, j)), abs(vxact(i, j)), vmax)
            pmax = max(abs(pxact(i, j)), pmax)
        end do
    end do

    ! Initialize first time step
    u = uxact
    v = vxact
    p = pxact
    ug = u
    vg = v
    pg = p

    !  Compute spectral coeffs of initial vrt, div, p
    call solver%get_vrtdivspec(ug, vg, vrtnm, divnm)
    call solver%grid_to_spec(pg, pnm)

    ! Time step loop
    new = 1
    now = 2
    old = 3
    do ncycle=0, max_iter

        time = real(ncycle, kind=wp) * dt

        !  Inverse transform to get vort and phig on grid
        call solver%spec_to_grid(vrtnm, vrtg)
        call solver%spec_to_grid(pnm, pg)

        ! Compute u and v on grid from spectral coeffs of vort and div.
        call solver%get_uv(vrtnm, divnm, ug, vg)

        ! Calculate error statistics
        if (mod(ncycle, mprint) == 0) then

            call solver%spec_to_grid(divnm, divg)

            u = ug
            v = vg
            p = pg
            htime = time/3600.0_wp

            write( stdout, '(/a)' ) ' steady nonlinear rotated flow:'
            write( stdout, fmt=write_fmt ) &
                ' cycle number              ', ncycle, &
                ' model time in  hours      ', htime, &
                ' time step in seconds      ', dt, &
                ' number of latitudes       ', NLAT, &
                ' number of longitudes      ', NLON, &
                ' max wave number           ', NTRUNC, &
                ' rotation rate        '     , omega, &
                ' mean height          '     , pzero, &
                ' maximum velocity     '     , uzero, &
                ' tilt angle           '     , alpha_in_degrees

            block
                real(wp) :: dvgm, dvmax, dpmax, evmax, epmax

                dvgm = ZERO
                dvmax = ZERO
                dpmax = ZERO
                evmax = ZERO
                epmax = ZERO

                do j=1, NLAT
                    do i=1, NLON
                        dvgm = max(dvgm, abs(divg(i, j)))
                        dvmax = dvmax+(u(i, j)-uxact(i, j))**2 + (v(i, j)-vxact(i, j))**2
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

                write( stdout, fmt='(2(a, 1pe15.6)/, a, 1pe15.6)') &
                    ' max error in velocity', evmax, &
                    ' max error in geopot. ', epmax, &
                    ' l2 error in velocity ', dvmax, &
                    ' l2 error in geopot.  ', dpmax, &
                    ' maximum divergence   ', dvgm
            end block
        end if

        !  Compute right-hand sides of prognostic eqns
        scrg1 = ug * (vrtg + f)
        scrg2 = vg * (vrtg + f)

        call solver%get_vrtdivspec(scrg1, scrg2, ddivdtnm(:,new), dvrtdtnm(:,new))

        dvrtdtnm(:,new) = -dvrtdtnm(:,new)
        scrg1 = ug * (pg + pzero)
        scrg2 = vg * (pg + pzero)

        call solver%get_vrtdivspec(scrg1, scrg2, scrnm, dpdtnm(:,new))

        dpdtnm(:,new) = -dpdtnm(:,new)
        scrg1 = pg + HALF * (ug**2 + vg**2)

        call solver%grid_to_spec(scrg1, scrnm)

        associate( lap => solver%laplacian_coefficients )
            ddivdtnm(:,new) = ddivdtnm(:,new) - lap * scrnm
        end associate

        ! Update vrt and div with third-order adams-bashforth
        ! Forward Euler, then 2nd-order Adams-Bashforth time steps to start
        select case (ncycle)
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

        vrtnm = vrtnm + dt * ( &
            (23.0_wp/12) * dvrtdtnm(:,new) &
            - (16.0_wp/12) * dvrtdtnm(:,now) &
            + (5.0_wp/12) * dvrtdtnm(:,old) )

        divnm = divnm + dt*( &
            (23.0_wp/12) * ddivdtnm(:,new) &
            - (16.0_wp/12) * ddivdtnm(:,now) &
            + (5.0_wp/12) * ddivdtnm(:,old) )

        pnm = pnm + dt*( &
            (23.0_wp/12) * dpdtnm(:,new) &
            - (16.0_wp/12) * dpdtnm(:,now) &
            + (5.0_wp/12) * dpdtnm(:,old) )

        !  Switch indices
        temp_save_new = new
        temp_save_now = now
        new = old
        now = temp_save_new
        old = temp_save_now
    end do

    !  Release memory
    call solver%destroy()
    deallocate( write_fmt )

end program shallow
