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
! Purpose:
!
! Program shallow solves the nonlinear shallow-water equations
! on the sphere using spherepack software.
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
!     ncycle        cycle number
!     time          model time in seconds
!     dt            time step
!     lambda        longitude
!     theta         colatitude
!
!   the first dimension of the following two dimensional arrays
!   corresponds to the latitude index with values i=1, ..., nlat
!   where i=1 is the north pole and i=nlat is the south pole.
!   the second dimension is longitude with values j=1, ..., nlon
!   where j=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!     u(i, j)       east longitudinal velocity component at t=time
!     v(i, j)       latitudinal velocity component at t=time
!     p(i, j)       +pzero = geopotential at t=time
!
!     unew(i, j)    east longitudinal velocity component at t=time+dt
!     vnew(i, j)    latitudinal velocity component at t=time+dt
!     pnew(i, j)    +pzero = geopotential at t=time+dt
!
!     uold(i, j)    east longitudinal velocity component at t=time-dt
!     vold(i, j)    latitudinal velocity component at t=time-dt
!     pold(i, j)    +pzero = geopotential at t=time-dt
!
!     divg(i, j)    divergence (d/dtheta (cos(theta) v)
!                                          + du/dlambda)/cos(theta)
!     vort(i, j)    vorticity  (d/dtheta (cos(theta) u)
!                                          - dv/dlambda)/cos(theta)
!
!     ut(i, j)      latitudinal derivative of longitudinal
!                  velocity component
!     vt(i, j)      latitudinal derivative of latitudinal
!                  velocity component
!
!     dudt(i, j)    time derivative of longitudinal velocity component
!     dvdt(i, j)    time derivative of latitudinal  velocity component
!     dpdt(i, j)    time derivative of geopotential
!
!     gpdl(i, j)    first component of the gradient of p(i, j)
!                  the longitudinal derivative of the geopotential
!                  divided by the cosine of the latitude
!
!     gpdt(i, j)    second component of the gradient of p(i, j)
!                  the latitudinal derivative of the geopotential
!
!     uxact(i, j)   the "exact" longitudinal veloctiy component
!     vxact(i, j)   the "exact" latitudinal  veloctiy component
!     uxact(i, j)   the "exact" geopotential
!
!     f(i, j)       the coriolis force on rotated grid
!
!   the following two dimensional arrays are nonzero in the triangle
!   n=1, ..., nlat and m less than or equal to n.
!
!     a(m, n), b(m, n)    spectral coefficients of the geopotential
!
!     br(m, n), bi(m, n)  spectral coefficients of the velocity
!     cr(m, n), ci(m, n)  vector [u(i, j), v(i, j)]
!
!
!     phlt(i)      the coefficients in the cosine series
!                  representation of the unrotated geopotential
!
program shallow

    use spherepack

    ! Explicit typing only
    implicit none

    real(wp) :: aa
    real(wp) :: alpha
    real(wp) :: alphad
    
    real(wp) :: cfn
    
    real(wp) :: clh
    real(wp) :: cost
    real(wp) :: cth
    real(wp) :: cthclh
    real(wp) :: cthslh
    real(wp) :: ctime
    real(wp) :: dlath
    real(wp) :: dpmax
    real(wp) :: dt
    
    real(wp) :: degree_to_radian

    real(wp) :: dvmax, dvgm
    real(wp) :: epmax
    real(wp) :: evmax
    real(wp) :: fzero
    real(wp) :: htime
    integer(ip) :: i
    integer(ip) :: ierror
    integer(ip) :: j
    integer(ip) :: lwsha
    integer(ip) :: lwshs
    integer(ip) :: lwvha
    integer(ip) :: lwvhs
    integer(ip) :: lwvts
    integer(ip) :: ncycle
    integer(ip) :: nl
    real(wp) :: omega
    real(wp) :: p2max
    real(wp) :: pmax
    real(wp) :: pzero
    
    
    real(wp) :: slh
    real(wp) :: sint
    real(wp) :: sth
    real(wp) :: two_dt
    real(wp) :: that
    real(wp) :: theta
    real(wp) :: time
    real(wp) :: uhat
    real(wp) :: uzero
    real(wp) :: v2max
    real(wp) :: vmax

    integer(ip), parameter :: isym = 0, nt = 1, lwork = 40000
    integer(ip), parameter :: idp = 73, jdp = 144, mdab = 73, ndab = 73
    integer(ip), parameter :: itmax = 720, mprint = 72, max_wave_number = 42
    integer(ip), parameter :: nlat = 65, nlon = 128
    real(wp) :: work(lwork)
    real(wp), dimension(idp, jdp)   :: u, v, p, f
    real(wp), dimension(idp, jdp)   :: unew, vnew, pnew, uold, vold, pold
    real(wp), dimension(idp, jdp)   :: uxact, vxact, pxact, divg, vort
    real(wp), dimension(idp, jdp)   :: ut, vt, dudt, dvdt, dpdt, gpdt, gpdl
    real(wp), dimension(mdab, ndab) :: a, b, br, bi, cr, ci
    real(wp)                        :: phlt(361)
    real(wp), allocatable, dimension(:) :: wsha, wshs, wvha, wvhs, wvts
    real(wp), parameter :: ZERO = 0.0_wp, ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    real(wp) :: lambda, lhat

    lwsha = 70928
    lwshs = 70928
    lwvha = 141647
    lwvhs = 141647
    lwvts = 141647

    degree_to_radian = PI/180
    aa = 6.37122e+6_wp
    omega = 7.292e-5_wp
    fzero = TWO * omega
    uzero = 40.0_wp
    pzero = 2.94e4_wp
    alphad = 60.0_wp
    alpha = degree_to_radian * alphad
    dt = 600.0_wp
    two_dt = TWO * dt

    ! Initialize wavetables
    call initialize_shaes(nlat, nlon, wsha, ierror)
    call initialize_shses(nlat, nlon, wshs, ierror)
    call initialize_vhaes(nlat, nlon, wvha, ierror)
    call initialize_vhses(nlat, nlon, wvhs, ierror)
    call initialize_vtses(nlat, nlon, wvts, ierror)

    !
    !
    !     compute the derivative of the unrotated geopotential
    !             p as a function of latitude
    !
    nl = 91
    cfn = ONE/(nl - 1)
    dlath = PI/(nl - 1)
    do i=1, nl - 2
        theta = real(i, kind=wp) * dlath
        sth = sin(theta)
        cth = cos(theta)
        uhat = ui(uzero, HALF_PI-theta)
        phlt(i) = cfn*cth*uhat*(uhat/sth+aa*fzero)
    end do
    !
    !     compute sine transform of the derivative of the geopotential
    !     for the purpose of computing the geopotential by integration
    !     see equation (3.9) in reference [1] above
    !
    call sine_transform(phlt(1:nl-2))
    !
    !     compute the cosine coefficients of the unrotated geopotential
    !     by the formal integration of the sine series representation
    !
    do i=1, nl - 2
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
    block
        real(wp) :: cosa, sina, dlam, dtheta
        real(wp) :: cosl, sinl

        cosa = cos(alpha)
        sina = sin(alpha)
        dtheta = pi/(nlat-1)
        dlam = TWO_PI/nlon
        do j=1, nlon
            lambda = real(j - 1, kind=wp) * dlam
            cosl = cos(lambda)
            sinl = sin(lambda)
            do i=1, nlat
                !
                !     lambda is longitude, theta is colatitude, and pi/2-theta is
                !     latitude on the rotated grid. lhat and that are longitude
                !     and colatitude on the unrotated grid. see text starting at
                !     equation (3.10)
                !
                theta = (i-1)*dtheta
                sint = cos(theta)
                cost = sin(theta)
                sth = cosa*sint+sina*cost*cosl
                cthclh = cosa*cost*cosl-sina*sint
                cthslh = cost*sinl
                lhat = atanxy(cthclh, cthslh)
                clh = cos(lhat)
                slh = sin(lhat)
                cth = clh*cthclh+slh*cthslh
                that = atanxy(sth, cth)
                uhat = ui(uzero, HALF_PI-that)
                pxact(i, j) = cosine(that, nl-2, phlt)
                uxact(i, j) = uhat*(cosa*sinl*slh+cosl*clh)
                vxact(i, j) = uhat*(cosa*cosl*slh*sint-clh*sinl*sint+sina*slh*cost)
                f(i, j) = fzero*sth
            end do
        end do
    end block

    vmax = ZERO
    pmax = ZERO
    v2max = ZERO
    p2max = ZERO
    do  j=1, nlon
        do i=1, nlat
            v2max = v2max+uxact(i, j)**2+vxact(i, j)**2
            p2max = p2max+pxact(i, j)**2
            vmax = max(abs(uxact(i, j)), abs(vxact(i, j)), vmax)
            pmax = max(abs(pxact(i, j)), pmax)
        end do
    end do

    ! Initialize first time step
    u(1: nlat, 1: nlon) = uxact(1: nlat, 1: nlon)
    v(1: nlat, 1: nlon) = vxact(1: nlat, 1: nlon)
    p(1: nlat, 1: nlon) = pxact(1: nlat, 1: nlon)

    time = ZERO
    ctime = ZERO
    ncycle = ZERO
    !
    !     start of the time loop
    !
    !   begin step 1, section 3
    !
    !     analyze the velocity components (u, v)
    !
    do while (ncycle <= itmax)
        call vhaesgo(nlat, nlon, isym, nt, u, v, idp, jdp, br, bi, cr, ci, &
            mdab, ndab, wvha, ierror)
        if (ierror /= 0) write (*, 91) ierror
91      format(' error' i4 ' in vhaes')
        !
        !     truncate spectrum to eliminate aliasing of the
        !     product terms in the shallow-water equations
        !
        call trunc(nlat, max_wave_number, mdab, br, bi)
        call trunc(nlat, max_wave_number, mdab, cr, ci)
        !
        !     resynthesize the velocity components
        !
        call vhsesgo(nlat, nlon, isym, nt, u, v, idp, jdp, br, bi, cr, ci, &
            mdab, ndab, wvhs, lwvhs, work, lwork, ierror)
        if (ierror /= 0) write (*, 92) ierror
92      format(' error' i4 ' in vhses')
        !
        !   begin step 2, section 3
        !
        !     analyze geopotential p
        !
        call shaes(nlat, nlon, isym, nt, p, idp, jdp, a, b, mdab, ndab, &
            wsha, ierror)
        if (ierror /= 0) write (*, 93) ierror
93      format(' error' i4 ' in shaes')
        !
        !     truncate spectrum to eliminate aliasing of the
        !     product terms in the shallow-water equations
        !
        call trunc(nlat, max_wave_number, mdab, a, b)
        !
        !     resynthesize the geopotential p
        !
        call shses(nlat, nlon, isym, nt, p, idp, jdp, a, b, mdab, ndab, &
            wshs, ierror)
        if (ierror /= 0) write (*, 94) ierror
94      format(' error' i4 ' in shses')
        !
        !
        !   begin step 3, section 3
        !
        !     compute the vorticity of the velocity (u, v)
        !
        call vrtes(nlat, nlon, isym, nt, vort, idp, jdp, cr, ci, mdab, ndab, &
            wshs, lwshs, work, lwork, ierror)
        if (ierror /= 0) write (*, 95) ierror
95      format(' error' i4 ' in vrtes')
        !
        !     compute the divergence of the velocity (u, v)
        !
        call dives(nlat, nlon, isym, nt, divg, idp, jdp, br, bi, mdab, ndab, &
            wshs, ierror)
        if (ierror /= 0) write (*, 96) ierror
96      format(' error' i4 ' in dives')
        !
        !   begin step 4, section 3
        !
        !     compute the derivative of the velocity (u, v) with
        !     respect to colatitude theta.
        !
        call vtsesgo(nlat, nlon, isym, nt, ut, vt, idp, jdp, br, bi, cr, ci, &
            mdab, ndab, wvts, lwvts, work, lwork, ierror)
        if (ierror /= 0) write (*, 97) ierror
97      format(' error' i4 ' in vtsesgo')
        !
        !   begin step 5, section 3
        !
        !     compute the gradient of the geopotential p
        !
        call gradesgo(nlat, nlon, isym, nt, gpdl, gpdt, idp, jdp, a, b, mdab, ndab, &
            wvhs, lwvhs, work, lwork, ierror)
        if (ierror /= 0) write (*, 98) ierror
98      format(' error' i4 ' in grades')
        !
        !     compute the time derivatives of the velocity (u, v)
        !     and the geopotential p using the shallow-water
        !     equations (2.8), (2.9), and (2.10), section 3.
        !
        dudt(1: nlat, 1: nlon) = (u(1: nlat, 1: nlon)*(vt(1: nlat, 1: nlon) &
            -divg(1: nlat, 1: nlon))-v(1: nlat, 1: nlon)*ut(1: nlat, 1: nlon) &
            -gpdl(1: nlat, 1: nlon))/aa+f(1: nlat, 1: nlon)*v(1: nlat, 1: nlon)
        dvdt(1: nlat, 1: nlon) = -(u(1: nlat, 1: nlon)*(vort(1: nlat, 1: nlon) &
            +ut(1: nlat, 1: nlon))+v(1: nlat, 1: nlon)*vt(1: nlat, 1: nlon) &
            +gpdt(1: nlat, 1: nlon))/aa-f(1: nlat, 1: nlon)*u(1: nlat, 1: nlon)
        dpdt(1: nlat, 1: nlon) = -((p(1: nlat, 1: nlon)+pzero)*divg(1: nlat, 1: nlon) &
            +v(1: nlat, 1: nlon)*gpdt(1: nlat, 1: nlon) &
            +u(1: nlat, 1: nlon)*gpdl(1: nlat, 1: nlon))/aa

        if (mod(ncycle, mprint) == 0) then
            htime = time/3600.

            write (*, 390) ncycle, htime, dt, nlat, nlon, max_wave_number, omega, pzero, &
                uzero, alphad
390         format(//' steady nonlinear rotated flow, test case 3'/ &
                ' cycle number              ' i10 &
                ' model time in  hours      ' f10.2/ &
                ' time step in seconds      ' f10.0 &
                ' number of latitudes       ' i10/ &
                ' number of longitudes      ' i10 &
                ' max wave number           ' i10/ &
                ' rotation rate        ' 1pe15.6 &
                ' mean height          ' 1pe15.6/ &
                ' maximum velocity     ' 1pe15.6 &
                ' tilt angle                ' f10.2)
            dvgm = ZERO
            dvmax = ZERO
            dpmax = ZERO
            evmax = ZERO
            epmax = ZERO
            do j=1, nlon
                do i=1, nlat
                    dvgm = max(dvgm, abs(divg(i, j)))
                    dvmax = dvmax+(u(i, j)-uxact(i, j))**2+(v(i, j)-vxact(i, j))**2
                    dpmax = dpmax+(p(i, j)-pxact(i, j))**2
                    evmax = max(evmax, abs(v(i, j)-vxact(i, j)), abs(u(i, j)-uxact(i, j)))
                    epmax = max(epmax, abs(p(i, j)-pxact(i, j)))
                end do
            end do

            dvmax = sqrt(dvmax/v2max)
            dpmax = sqrt(dpmax/p2max)
            evmax = evmax/vmax
            epmax = epmax/pmax

            write (*, 391) evmax, epmax, dvmax, dpmax, dvgm
391         format(' max error in velocity' 1pe15.6 &
                ' max error in geopot. ' 1pe15.6/ &
                ' l2 error in velocity ' 1pe15.6 &
                ' l2 error in geopot.  ' 1pe15.6/ &
                ' maximum divergence   ' 1pe15.6)
        end if

        ! Set values at time = -dt to values at time = 0.
        if (ncycle == 0) then
            uold(1: nlat, 1: nlon) = u(1: nlat, 1: nlon)
            vold(1: nlat, 1: nlon) = v(1: nlat, 1: nlon)
            pold(1: nlat, 1: nlon) = p(1: nlat, 1: nlon)
        end if

        ! Compute values at next time level using leap frog
        ! time differencing
        unew(1: nlat, 1: nlon) = uold(1: nlat, 1: nlon)+two_dt*dudt(1: nlat, 1: nlon)
        vnew(1: nlat, 1: nlon) = vold(1: nlat, 1: nlon)+two_dt*dvdt(1: nlat, 1: nlon)
        pnew(1: nlat, 1: nlon) = pold(1: nlat, 1: nlon)+two_dt*dpdt(1: nlat, 1: nlon)

        ! Update values to next time level
        uold(1: nlat, 1: nlon) = u(1: nlat, 1: nlon)
        vold(1: nlat, 1: nlon) = v(1: nlat, 1: nlon)
        pold(1: nlat, 1: nlon) = p(1: nlat, 1: nlon)
        u(1: nlat, 1: nlon) = unew(1: nlat, 1: nlon)
        v(1: nlat, 1: nlon) = vnew(1: nlat, 1: nlon)
        p(1: nlat, 1: nlon) = pnew(1: nlat, 1: nlon)

        ncycle = ncycle + 1
        time = time + dt
    end do

    ! Release memory
    deallocate (wsha)
    deallocate (wshs)
    deallocate (wvha)
    deallocate (wvhs)
    deallocate (wvts)

contains

    subroutine vtsesgo(nlat, nlon, ityp, nt, ut, vt, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvts, lwvts, work, lwork, ierror)

        real(wp), dimension(..) :: br, bi, cr, ci
        integer(ip) :: i
        integer(ip) :: idvw
        integer(ip) :: ierror
        integer(ip) :: ityp
        integer(ip) :: j
        integer(ip) :: jdvw
        integer(ip) :: k
        integer(ip) :: lwork
        integer(ip) :: lwvts
        integer(ip) :: mdab
        integer(ip) :: ndab
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: ut(idvw, jdvw, *), vt(idvw, jdvw, *)
        real(wp) :: wvts(:), work(:)
        !
        !     vtsesgo computes the latitudinal derivatives of the
        !     velocity components using subroutine vtses which
        !     assumes the velocity components are given in terms
        !     of mathematical coordinates

        call vtses(nlat, nlon, ityp, nt, vt, ut, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, ierror)

        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    ut(i, j, k) = -ut(i, j, k)
                end do
            end do
        end do

    end subroutine vtsesgo

    !     computes the initial unrotated longitudinal velocity
    !     see section 3.3.
    pure function ui(amp, thetad) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: amp
        real(wp), intent(in) :: thetad
        real(wp)             :: return_value

        ! Local variables
        real(wp) :: thetae, thetab, x, xe

        thetab = -PI/6
        thetae = HALF_PI
        xe = 3.0e-1_wp

        x = xe * (thetad-thetab)/(thetae-thetab)
        if (x <= ZERO .or. x >= xe) then
            return_value = ZERO
        else
            return_value = amp * exp(-ONE/x-ONE/(xe-x)+4./xe)
        end if

    end function ui

    pure function atanxy(x, y) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: x
        real(wp), intent(in) :: y
        real(wp)             :: return_value

        return_value = ZERO

        if (x == ZERO .and. y == ZERO) return

        return_value = atan2(y, x)

    end function atanxy

    subroutine sine_transform(x)

        ! Dummy arguments
        real(wp), intent(inout) :: x(:)

        ! Local variables
        integer(ip) :: i, j

        associate (n => size(x))
            block
                real(wp) :: w(n), arg

                arg = PI/(n+1)

                do j=1, n
                    w(j) = ZERO
                    do i=1, n
                        w(j) = w(j)+x(i)*sin(real(i*j, kind=wp)*arg)
                    end do
                end do

                x = TWO * w
            end block
        end associate

    end subroutine sine_transform

    real function cosine(theta, n, cf)

        real(wp) :: cf
        integer(ip) :: i
        integer(ip) :: n
        real(wp) :: theta
        !
        !     computes the cosine transform
        !
        dimension cf(n)
        cosine = ZERO
        do i=1, n
            cosine = cosine+cf(i)*cos(i*theta)
        end do

    end function cosine

    subroutine trunc(nm, ms, id, a, b)

        real(wp) :: a(id,*)
        real(wp) :: b(id,*)
        integer(ip) :: id
        integer(ip) :: m
        
        integer(ip) :: ms
        integer(ip) :: n
        integer(ip) :: nm
        !
        !     truncates spectral coefficients so that aliasing
        !     does not occur when computing the spectral representations
        !     of the product terms.
        !

        do n=ms + 2, nm
            do m=1, n
                a(m, n) = ZERO
                b(m, n) = ZERO
            end do
        end do

    end subroutine trunc

    subroutine vhaesgo(nlat, nlon, ityp, nt, u, v, iduv, jduv, &
        br, bi, cr, ci, mdab, ndab, wsav, ierror)

        real(wp) :: bi
        real(wp) :: br
        real(wp) :: ci
        real(wp) :: cr
        integer(ip) :: i
        integer(ip) :: iduv
        integer(ip) :: ierror
        integer(ip) :: ityp
        integer(ip) :: j
        integer(ip) :: jduv
        integer(ip) :: k
        
        
        integer(ip) :: mdab
        integer(ip) :: ndab
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: u
        real(wp) :: v
        real(wp) :: wsav(:)
        dimension u(iduv, jduv, *), v(iduv, jduv, *), br(mdab, ndab, *), &
            bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)

        !     vhaesgo computes the vector harmonic analysis of (u, v) using vhaes which
        !     assumes the velocity components are given in mathematical coordinates
        !
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    v(i, j, k) = -v(i, j, k)
                end do
            end do
        end do

        call vhaes(nlat, nlon, ityp, nt, v, u, iduv, jduv, &
            br, bi, cr, ci, mdab, ndab, wsav, ierror)
        !
        !     restore v
        !
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    v(i, j, k) = -v(i, j, k)
                end do
            end do
        end do

        if (ierror /= 0) return

    end subroutine vhaesgo

    subroutine vhsesgo(nlat, nlon, ityp, nt, u, v, iduv, jduv, &
        br, bi, cr, ci, mdab, ndab, wsav, lwsav, work, lwork, ierror)
        implicit none
        real(wp) :: bi
        real(wp) :: br
        real(wp) :: ci
        real(wp) :: cr
        integer(ip) :: i
        integer(ip) :: iduv
        integer(ip) :: ierror
        integer(ip) :: ityp
        integer(ip) :: j
        integer(ip) :: jduv
        integer(ip) :: k
        integer(ip) :: lwork
        integer(ip) :: lwsav
        integer(ip) :: mdab
        integer(ip) :: ndab
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: u
        real(wp) :: v
        real(wp) :: work(..)
        real(wp) :: wsav(:)
        dimension u(iduv, jduv, *), v(iduv, jduv, *), br(mdab, ndab, *), &
            bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
        !
        !     vhsesgo computes a vector harmonic synthesis in (u, v) using vhses which
        !     assumes the velocity components are given in mathematical coordinates
        !
        call vhses(nlat, nlon, ityp, nt, v, u, iduv, jduv, &
            br, bi, cr, ci, mdab, ndab, wsav, ierror)

        if (ierror /= 0) return

        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    v(i, j, k) = -v(i, j, k)
                end do
            end do
        end do

    end subroutine vhsesgo

    subroutine gradesgo(nlat, nlon, isym, nt, u, v, iduv, jduv, a, b, &
        mdab, ndab, wsav, lwsav, work, lwork, ierror)
        implicit none
        real(wp) :: a
        real(wp) :: b
        integer(ip) :: i
        integer(ip) :: iduv
        integer(ip) :: ierror
        integer(ip) :: isym
        integer(ip) :: j
        integer(ip) :: jduv
        integer(ip) :: k
        integer(ip) :: lwork
        integer(ip) :: lwsav
        integer(ip) :: mdab
        integer(ip) :: ndab
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: u
        real(wp) :: v
        real(wp) :: work
        real(wp) :: wsav
        dimension u(iduv, jduv, nt), v(iduv, jduv, nt)
        dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
        dimension wsav(lwsav), work(lwork)
        !
        !     gradesgo computes the gradient in (u, v) using grades which assumes
        !     the velocity components are given in mathematical coordinates
        !
        call grades(nlat, nlon, isym, nt, v, u, iduv, jduv, a, b, &
            mdab, ndab, wsav, ierror)

        if (ierror /= 0) return

        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    v(i, j, k) = -v(i, j, k)
                end do
            end do
        end do

    end subroutine gradesgo

end program shallow
