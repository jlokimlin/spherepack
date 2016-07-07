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
program shallow

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use, intrinsic :: iso_fortran_env, only: &
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
    real (wp)                         :: theta_mesh, lambda, lhat, uhat, aa, uzero, pzero, HALF_PI, radian_unit, omega, alphad
    real (wp)                         :: alpha, fzero, dt, cfn, dlath, theta, sth
    real (wp)                         :: cth, cosa, sina, lambda_mesh, st, ct, cthclh, cthslh
    real (wp)                         :: clh, time, that, sl, slh, evmax, epmax, dvmax, dpmax, htime, dvgm, cl
    real (wp)                         :: v2max, p2max, vmax, pmax
    complex (wp), dimension(nmdim)    :: vrtnm, divnm, pnm, scrnm
    complex (wp), dimension(nmdim, 3) :: dvrtdtnm, ddivdtnm, dpdtnm
    type (ShallowWaterSolver)         :: solver
    character (len=:), allocatable    :: write_fmt
    !--------------------------------------------------------------------------------

    write( stdout, '(/a)') '     shallow *** TEST RUN *** '

    !
    !==> Initialize constants
    !
    HALF_PI = PI/2
    radian_unit = PI/180.0_wp
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
        '(a, i10, a, f10.2/, a, f10.0, a, i10/, a, i10, '&
        //'a, i10/, a, 1pe15.6, a, 1pe15.6, /a, 1pe15.6, a, 1pe15.6)' )

    !
    !==> compute the derivative of the unrotated geopotential
    !    p as a function of latitude
    !
    cfn = 1.0_wp / nlm1
    dlath = PI / nlm1
    do i=1, nlm2
        theta = i*dlath
        sth = sin(theta)
        cth = cos(theta)
        uhat = solver%get_initial_velocity(uzero, HALF_PI-theta)
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
    !     below to compute the geopotential on the rotated grid.
    !
    !     compute the initial values of  east longitudinal
    !     and latitudinal velocities u and v as well as the
    !     geopotential p and coriolis f on the rotated grid.
    !
    cosa = cos(alpha)
    sina = sin(alpha)
    lambda_mesh = (2.0_wp*PI)/nlon
    theta_mesh = PI/(nlat-1)

    do j=1, nlon
        lambda = real(j - 1, kind=wp) * lambda_mesh
        cl = cos(lambda)
        sl = sin(lambda)
        do i=1, nlat
            !
            !     lambda is longitude, theta is colatitude, and pi/2-theta is
            !     latitude on the rotated grid. lhat and that are longitude
            !     and colatitude on the unrotated grid. see text starting at
            !     equation (3.10)
            !
            theta = real(i - 1, kind=wp)*theta_mesh
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
            uhat = solver%get_initial_velocity(uzero, HALF_PI-that)
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

            write( stdout, '(/a)' ) ' steady nonlinear rotated flow:'
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

            write( stdout, fmt='(2(a, 1pe15.6)/, a, 1pe15.6)') &
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
