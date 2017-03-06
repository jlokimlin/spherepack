!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK                               *
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
! ... file testrvsph.f
!
!     this file contains a test program for trvsph.f
!
! ... required files
!
!     trvsph.f, sphcom.f, hrfft.f, gaqd.f, vhaec.f, vhsec.f, vhagc.f, vhsgc.f
!
! ... description (see documentation in file trvsph.f)
!
!     Subroutine trvsph is used to demonstrate vector data transfer between a coarse
!     ten degree equally spaced grid and a higher resolution T42 Global Spectral
!     Gaussian grid.
!
!     Ten Degree Grid (Mathematical Spherical Coordinates)
!
!     The equally spaced vector data, given in mathematical spherical
!     coordinates,  is stored in 19 X 36 colatitude-longitude arrays (ue, ve).  The
!     colatitudinal values are stored north to south with increasing colatitude
!     subscript values.  ue is the east longitudinal component and ve is the
!     colatitudinal component of the vector.
!
!     Gaussian Grid (Geophysical Spherical Coordinates)
!
!     The T42 Gaussian grid vector data is in geophysical spherical coordinates.
!     It is stored in 128 X 64 longitude-latitude arrays (ug, vg).  Values are
!     stored south to north with increasing latitude subscript value.  ug
!     is the longitudinal component and vg is the latitudinal componenet of the
!     vector data.
!
!     For testing purposes we use an analytic vector field (u, v).  Let t and p
!     be colatitude and longitude and x=sin(t)cos(p), y=sint(t)*sin(p), z=cos(t)
!     be the cartesian coordinates restricted to the sphere in mathematical
!     coordinates. We derive the vector field (u, v) from the stream function
!     S and velocity function P given by
!
!               y    -z           x    z
!          S = e  + e  ,     P = e  + e
!
!     The corresponding vector field has the form
!
!                x          -z          y
!          u = -e sin(p) + e  sin(t) + e cos(t)sin(p)
!
!                y          z          x
!          v = -e cos(p) - e sin(t) + e cos(t)cos(p)
!
!
!     in mathematical spherical coordinates.  Values in geophysical
!     coordinates can be obtained by negating v.
!
!     In the code below, (ue, ve) is set equal to (u, v) and trvsph is used to
!     transfer (ue, ve) to (ug, vg).  (ug, vg) is then compared with (u, v) in
!     geophysical coordinates on the T42 Gaussian grid.  Finally, trvsph is
!     used to transfer (ug, vg) back to (ue, vg) and this is compared with (u, v).
!     Output from executing the test program on separate platforms with 32 bit
!     and 64 bit floating point arithmetic is listed below.  The minimum required
!     saved and unsaved work space lengths were predetermined by an earlier call
!     trvsph with nlone=36, nlate=19, nlong=128, nlatg=64, lsave=0, lwork=0 and printout
!     of lsvmin and lwkmin.
!
!
! **********************************************************************
!
!     OUTPUT FROM EXECUTING CODE IN THIS FILE
!
! **********************************************************************
!
!      EQUALLY SPACED TO GAUSSIAN GRID TRANSFER
!      trvsph input arguments:
!      intl =  0
!      igride(1) = -1   igride(2) =  1
!      nlone =  36   nlate =  19
!      ive =  1
!      igridg(1) =  2   igridg(2) =  0
!      nlong = 128   nlatg =  64
!      ivg =  0
!      lsave =   21814   lwork =   75173
!
!      trvsph output:
!      ier =        0  lsvmin =   21792  lwkmin =   75173
! ***  32 BIT FLOATING POINT ARITHMETIC
!      least squares error in u =  0.307E-06
!      least squares error in v =  0.272E-06
! ***  64 BIT FLOATING POINT ARITHMETIC
!      least squares error in u =  0.841E-12
!      least squares error in v =  0.603E-12
!
!      GAUSSIAN TO EQUALLY SPACED GRID TRANSFER
!      trvsph input arguments:
!      intl =  0
!      igridg(1) =  2   igridg(2) =  0
!      nlong = 128   nlatg =  64
!      ivg =  0
!      igride(1) = -1   igride(2) =  1
!      nlone =  36   nlate =  19
!      ive =  1
!      lsave =   21814   lwork =   75173
!
!      trvsph output:
!      ier =        0  lsvmin =   21814  lwkmin =   75173
! ***  32 BIT FLOATING POINT ARITHMETIC
!      least squares error in u =  0.374E-06
!      least squares error in v =  0.364E-06
! ***  64 BIT FLOATING POINT ARITHMETIC
!      least squares error in u =  0.170E-12
!      least squares error in v =  0.161E-12
!
! **********************************************************************
!
!     END OF OUTPUT ... CODE FOLLOWS
!
! **********************************************************************
program testrvsph

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack

    ! Explicit typing only
    implicit none

    !     set grid sizes with parameter statements
    integer(ip), parameter :: nlatg=64, nlong=128, nlate=19, nlone=36
    !
    !     set predetermined minimum saved and unsaved work space lengths
    !
    integer(ip), parameter :: lwork =  76000
    integer(ip), parameter :: lsave = 21814

    !     set real workspace length
    integer(ip), parameter :: ldwork = 2*nlatg*(nlatg+1)+1
    !
    !     dimension and type data arrays and grid vectors and internal variables
    !
    real(wp) :: ue(nlate, nlone), ve(nlate, nlone)
    real(wp) :: ug(nlong, nlatg), vg(nlong, nlatg)
    real(wp) :: work(lwork), wsave(lsave), thetag(nlatg)
    real(wp) :: dtheta(nlatg), dwts(nlatg), dwork(ldwork)
    integer(ip) :: igride(2), igridg(2), ive, ivg
    real(wp) :: dlate, dlone, dlong, t, p, cosp, sinp, cost, sint, x, y, z
    real(wp) :: erru2, errv2, ex, ey, ez, emz, uee, vee
    integer(ip) :: i, j, ib, intl, error_flag, lsvmin, lwkmin

    ! Set equally spaced grid increments
    dlate = PI/(nlate-1)
    dlone = TWO_PI/nlone
    dlong = TWO_PI/nlong

    ! Set vector data in (ue, ve)
    do  j=1, nlone
        p = (j-1)*dlone
        cosp = cos(p)
        sinp = sin(p)
        do i=1, nlate
            t = (i-1)*dlate
            cost = cos(t)
            sint = sin(t)
            x = sint*cosp
            y = sint*sinp
            z = cost
            ex = exp(x)
            ey = exp(y)
            ez = exp(z)
            emz = exp(-z)
            ve(i, j) = ex*cost*cosp-ey*cosp-ez*sint
            ue(i, j) = -ex*sinp + emz*sint + ey*cost*sinp
        end do
    end do
    !
    !     set initial call flag
    !
    intl = 0
    !
    !     flag (ue, ve) grid as north to south equally spaced
    !
    igride(1) = -1
    !
    !     flag (ue, ve) as nlate by nlone arrays
    !
    igride(2) = 1
    !
    !     flag ve as colatitude component of vector
    ive = 1
    !
    !     flag (ug, vg) as south to north gaussian
    !
    igridg(1) = 2
    !
    !     flag (ug, vg) as nlong by nlatg arrays
    !
    igridg(2) = 0
    !
    !     flag vg as latitude component of vector
    ivg = 0
    !
    !     print trvsph input arguments
    !
    write (stdout, 100) intl, igride(1), igride(2), nlone, nlate, ive, &
        igridg(1), igridg(2), nlong, nlatg, ivg, lsave, lwork, ldwork
100 format(//' EQUALLY SPACED TO GAUSSIAN GRID TRANSFER ' , &
        /' trvsph input arguments: ' , &
        /' intl = ', i2, &
        /' igride(1) = ', i2, 2x, ' igride(2) = ', i2, &
        /' nlone = ', i3, 2x, ' nlate = ', i3, &
        /' ive = ', i2, &
        /' igridg(1) = ', i2, 2x, ' igridg(2) = ', i2, &
        /' nlong = ', i3, 2x, ' nlatg = ', i3, &
        /' ivg = ', i2 &
        /' lsave = ', i7, 2x, ' lwork = ', i7, 2x, ' ldwork = ', i5)
    !
    !     transfer  (ue, ve) to (ug, vg)
    !
    call trvsph(intl, igride, nlone, nlate, ive, ue, ve, igridg, nlong, &
        nlatg, ivg, ug, vg, wsave, lsave, lsvmin, work, lwork, lwkmin, dwork, &
        ldwork, error_flag)
    !
    !     print output arguments
    !
    write (stdout, 200) error_flag, lsvmin, lwkmin
200 format(//' trvsph output: ' &
        / ' error_flag = ', i8, 2x, 'lsvmin = ', i7, 2x, 'lwkmin = ', i7)

    if (error_flag == 0) then
        !
        !     compute nlatg gaussian colatitude points and
        !     set with south to north orientation in thetag
        !
        call compute_gaussian_latitudes_and_weights(nlatg, dtheta, dwts, error_flag)
        do  i=1, nlatg
            ib = nlatg-i+1
            thetag(i) = dtheta(ib)
        end do
          !
          !     compute the least squares error in (ug, vg)
          !     by comparing with exact geophysical vector
          !
        errv2 = 0.0
        erru2 = 0.0
        do  j=1, nlong
            p = (j-1)*dlong
            cosp = cos(p)
            sinp = sin(p)
            do i=1, nlatg
                t = thetag(i)
                cost = cos(t)
                sint = sin(t)
                x = sint*cosp
                y = sint*sinp
                z = cost
                ex = exp(x)
                ey = exp(y)
                ez = exp(z)
                emz = exp(-z)
                vee = -ex*cost*cosp+ey*cosp+ez*sint
                uee = -ex*sinp + emz*sint + ey*cost*sinp
                erru2 = erru2 + (ug(j, i)-uee)**2
                errv2 = errv2 + (vg(j, i)-vee)**2
            end do
        end do
        erru2 = sqrt(erru2/(nlong*nlatg))
        errv2 = sqrt(errv2/(nlong*nlatg))
        write (stdout, 300) erru2, errv2
300     format(' least squares error in u = ', e10.3 &
            /' least squares error in v = ', e10.3)
    end if

    ! Now transfer (ug, vg) back to (ue, ve)
    ue = 0.0_wp
    ve = 0.0_wp

    write (stdout, 101) intl, igridg(1), igridg(2), nlong, nlatg, ivg, &
        igride(1), igride(2), nlone, nlate, ive, lsave, lwork, ldwork
101 format(//' GAUSSIAN TO EQUALLY SPACED GRID TRANSFER ' , &
        /' trvsph input arguments: ' , &
        /' intl = ', i2, &
        /' igridg(1) = ', i2, 2x, ' igridg(2) = ', i2, &
        /' nlong = ', i3, 2x, ' nlatg = ', i3, &
        /' ivg = ', i2, &
        /' igride(1) = ', i2, 2x, ' igride(2) = ', i2, &
        /' nlone = ', i3, 2x, ' nlate = ', i3, &
        /' ive = ', i2 &
        /' lsave = ', i7, 2x, ' lwork = ', i7, 2x, ' ldwork = ', i5)
    call trvsph(intl, igridg, nlong, nlatg, ivg, ug, vg, igride, nlone, nlate, &
        ive, ue, ve, wsave, lsave, lsvmin, work, lwork, lwkmin, dwork, ldwork, error_flag)
    !
    !     print output arguments
    !
    write (stdout, 200) error_flag, lsvmin, lwkmin
    if (error_flag == 0) then
        !
        !     compute the least squares error in (ue, ve)
        !     by comparing with exact mathematical vector
        !
        errv2 = 0.0
        erru2 = 0.0
        do  j=1, nlone
            p = (j-1)*dlone
            cosp = cos(p)
            sinp = sin(p)
            do i=1, nlate
                t = (i-1)*dlate
                cost = cos(t)
                sint = sin(t)
                x = sint*cosp
                y = sint*sinp
                z = cost
                ex = exp(x)
                ey = exp(y)
                ez = exp(z)
                emz = exp(-z)
                vee =  ex*cost*cosp-ey*cosp-ez*sint
                uee = -ex*sinp + emz*sint + ey*cost*sinp
                erru2 = erru2 + (ue(i, j)-uee)**2
                errv2 = errv2 + (ve(i, j)-vee)**2
            end do
        end do
        erru2 = sqrt(erru2/(nlone*nlate))
        errv2 = sqrt(errv2/(nlone*nlate))
        write (stdout, 300) erru2, errv2
    end if

end program testrvsph
