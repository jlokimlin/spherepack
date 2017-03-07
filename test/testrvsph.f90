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
!     trvsph with nlon_reg=36, nlat_reg=19, nlon_gau=128, nlat_gau=64, lsave=0, lwork=0 and printout
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
!      igrid_reg(1) = -1   igrid_reg(2) =  1
!      nlon_reg =  36   nlat_reg =  19
!      ive =  1
!      igrid_gau(1) =  2   igrid_gau(2) =  0
!      nlon_gau = 128   nlat_gau =  64
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
!      igrid_gau(1) =  2   igrid_gau(2) =  0
!      nlon_gau = 128   nlat_gau =  64
!      ivg =  0
!      igrid_reg(1) = -1   igrid_reg(2) =  1
!      nlon_reg =  36   nlat_reg =  19
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
    integer(ip), parameter :: NLAT_GAU=64, NLON_GAU=128, NLAT_REG=19, NLON_REG=36
    !
    !     dimension and type data arrays and grid vectors and internal variables
    !
    real(wp) :: u_reg(NLAT_REG, NLON_REG), v_reg(NLAT_REG, NLON_REG)
    real(wp) :: u_gau(NLON_GAU, NLAT_GAU), v_gau(NLON_GAU, NLAT_GAU)
    real(wp) :: theta_gau(NLAT_GAU)
    real(wp) :: gaussian_latitudes(NLAT_GAU), gaussian_weights(NLAT_GAU)
    integer(ip) :: igrid_reg(2), igrid_gau(2), iv_reg, iv_gau
    real(wp) :: dlat_reg, dlon_reg, dlon_gau
    real(wp) :: theta, phi, cosp, sinp, cost, sint, x, y, z
    real(wp) :: erru2, errv2, ex, ey, ez, emz, uee, vee
    integer(ip) :: i, j, ib, intl, error_flag

    ! Set equally spaced grid increments
    dlat_reg = PI/(NLAT_REG-1)
    dlon_reg = TWO_PI/NLON_REG
    dlon_gau = TWO_PI/NLON_GAU

    ! Set initial call flag
    intl = 0

    ! flag (ue, ve) grid as north to south equally spaced
    igrid_reg(1) = -1

    ! flag (ue, ve) as nlat_reg by nlon_reg arrays
    igrid_reg(2) = 1
    !     flag ve as colatitude component of vector
    iv_reg = 1

    ! flag (ug, vg) as south to north gaussian
    igrid_gau(1) = 2

    ! flag (ug, vg) as nlon_gau by nlat_gau arrays
    igrid_gau(2) = 0

    ! flag vg as latitude component of vector
    iv_gau = 0

    ! Set vector data in (ue, ve)
    do  j=1, NLON_REG
        phi = real(j - 1, kind=wp)*dlon_reg
        cosp = cos(phi)
        sinp = sin(phi)
        do i=1, NLAT_REG
            theta = real(i - 1, kind=wp)*dlat_reg
            cost = cos(theta)
            sint = sin(theta)
            x = sint*cosp
            y = sint*sinp
            z = cost
            ex = exp(x)
            ey = exp(y)
            ez = exp(z)
            emz = exp(-z)
            v_reg(i, j) = ex*cost*cosp-ey*cosp-ez*sint
            u_reg(i, j) = -ex*sinp + emz*sint + ey*cost*sinp
        end do
    end do

    ! print trvsph input arguments
    write (stdout, 100) intl, igrid_reg(1), igrid_reg(2), NLON_REG, NLAT_REG, iv_reg, &
        igrid_gau(1), igrid_gau(2), NLON_GAU, NLAT_GAU, iv_gau
100 format(//' EQUALLY SPACED TO GAUSSIAN GRID TRANSFER ' , &
        /' trvsph input arguments: ' , &
        /' intl = ', i2, &
        /' igrid_reg(1) = ', i2, 2x, ' igrid_reg(2) = ', i2, &
        /' nlon_reg = ', i3, 2x, ' nlat_reg = ', i3, &
        /' ive = ', i2, &
        /' igrid_gau(1) = ', i2, 2x, ' igrid_gau(2) = ', i2, &
        /' nlon_gau = ', i3, 2x, ' nlat_gau = ', i3, &
        /' ivg = ', i2)

    !     transfer  (ue, ve) to (ug, vg)
    call trvsph(intl, igrid_reg, NLON_REG, NLAT_REG, iv_reg, u_reg, v_reg, igrid_gau, NLON_GAU, &
        NLAT_GAU, iv_gau, u_gau, v_gau, error_flag)

    if (error_flag == 0) then
        !
        !     compute nlat_gau gaussian colatitude points and
        !     set with south to north orientation in thetag
        !
        call compute_gaussian_latitudes_and_weights( &
            NLAT_GAU, gaussian_latitudes, gaussian_weights, error_flag)
        do  i=1, NLAT_GAU
            ib = NLAT_GAU-i+1
            theta_gau(i) = gaussian_latitudes(ib)
        end do
          !
          !     compute the least squares error in (ug, vg)
          !     by comparing with exact geophysical vector
          !
        errv2 = 0.0
        erru2 = 0.0
        do  j=1, NLON_GAU
            phi = real(j-1, kind=wp)*dlon_gau
            cosp = cos(phi)
            sinp = sin(phi)
            do i=1, NLAT_GAU
                theta = theta_gau(i)
                cost = cos(theta)
                sint = sin(theta)
                x = sint*cosp
                y = sint*sinp
                z = cost
                ex = exp(x)
                ey = exp(y)
                ez = exp(z)
                emz = exp(-z)
                vee = -ex*cost*cosp+ey*cosp+ez*sint
                uee = -ex*sinp + emz*sint + ey*cost*sinp
                erru2 = erru2 + (u_gau(j, i)-uee)**2
                errv2 = errv2 + (v_gau(j, i)-vee)**2
            end do
        end do
        erru2 = sqrt(erru2/(NLON_GAU*NLAT_GAU))
        errv2 = sqrt(errv2/(NLON_GAU*NLAT_GAU))
        call print_least_squared_error(erru2, errv2)
    end if

    write (stdout, 101) intl, igrid_gau(1), igrid_gau(2), NLON_GAU, NLAT_GAU, iv_gau, &
        igrid_reg(1), igrid_reg(2), NLON_REG, NLAT_REG, iv_reg
101 format(//' GAUSSIAN TO EQUALLY SPACED GRID TRANSFER ' , &
        /' trvsph input arguments: ' , &
        /' intl = ', i2, &
        /' igrid_gau(1) = ', i2, 2x, ' igrid_gau(2) = ', i2, &
        /' nlon_gau = ', i3, 2x, ' nlat_gau = ', i3, &
        /' iv_gau = ', i2, &
        /' igrid_reg(1) = ', i2, 2x, ' igrid_reg(2) = ', i2, &
        /' nlon_reg = ', i3, 2x, ' nlat_reg = ', i3, &
        /' iv_reg = ', i2)

    call trvsph(intl, igrid_gau, NLON_GAU, NLAT_GAU, iv_gau, u_gau, v_gau, igrid_reg, &
        NLON_REG, NLAT_REG, iv_reg, u_reg, v_reg, error_flag)

    if (error_flag == 0) then
        !
        !     compute the least squares error in (ue, ve)
        !     by comparing with exact mathematical vector
        !
        errv2 = 0.0_wp
        erru2 = 0.0_wp
        do  j=1, NLON_REG
            phi = real(j-1, kind=wp)*dlon_reg
            cosp = cos(phi)
            sinp = sin(phi)
            do i=1, NLAT_REG
                theta = real(i-1, kind=wp)*dlat_reg
                cost = cos(theta)
                sint = sin(theta)
                x = sint*cosp
                y = sint*sinp
                z = cost
                ex = exp(x)
                ey = exp(y)
                ez = exp(z)
                emz = exp(-z)
                vee =  ex*cost*cosp-ey*cosp-ez*sint
                uee = -ex*sinp + emz*sint + ey*cost*sinp
                erru2 = erru2 + (u_reg(i, j)-uee)**2
                errv2 = errv2 + (v_reg(i, j)-vee)**2
            end do
        end do
        erru2 = sqrt(erru2/(NLON_REG*NLAT_REG))
        errv2 = sqrt(errv2/(NLON_REG*NLAT_REG))
        call print_least_squared_error(erru2, errv2)
    end if

contains

    subroutine print_least_squared_error(erru2, errv2)

        ! Dummy arguments
        real(wp), intent(in) :: erru2, errv2

        write (stdout, '(2(a,e10.3/))') &
            ' least squares error in u = ', erru2, &
            ' least squares error in v = ', errv2

    end subroutine print_least_squared_error

end program testrvsph
