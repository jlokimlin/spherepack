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
!
! ... file testrssph.f
!
!     this file contains a test program for subroutine trssph
!
! ... author
!
!     John C Adams (1996, NCAR)
!
! ... required files
!
!     subroutine trssph and the entire spherepack library
!
! ... description (see documentation in file trssph.f)
!
!     Subroutine trssph is used to demonstrate data transfer between a coarse
!     ten degree equally spaced grid and a higher resolution T64 Global Spectral
!     Gaussian grid.  The equally spaced data is stored in a 19 X 36 colatitude-
!     longitude array data_reg which runs north to south with increasing subscript
!     values.  The Gaussian grid data is stored in a 192 X 94 longitude-latitude
!     array data_gau which runs south to north with increasing latitude subscript.
!     First trssph is used to transfer data_reg to data_gau.  Then trssph is used to
!     transfer data_gau back to data_reg.
!
!     For testing purposes, data_reg is set equal the analytic function
!
!                       x*y*z
!           f(x, y, z) = e
!
!     in Cartesian coordinates x, y, z restricted to the surface of the sphere.
!     The same function is used to compute error in data_gau after the data transfer
!     with trssph.  Finally this is used to compute error in data_reg after the transfer
!     back with trssph.  Output from executing the test program on machines with
!     32 bit and 64 bit arithmetic is listed below.  The minimum required saved
!     and unsaved work space lengths were predetermined by an earlier call to
!     trssph with nlon_reg=36, nlat_reg=19, nlon_gau=192, nlat_gau=94, lsave=1, lwork=1 and printout
!     of lsvmin and lwkmin.
!
!
! **********************************************************************
!
!     OUTPUT FROM EXECUTING CODE IN THIS FILE
!
! **********************************************************************
!
!     EQUALLY SPACED TO GAUSSIAN GRID TRANSFER
!
!     trssph input parameters:
!     intl =  0
!     igrid_reg(1) = -1   igrid_reg(2) =  1
!     nlon_reg =  36   nlat_reg =  19
!     igrid_gau(1) =  2   igrid_gau(2) =  0
!     nlon_gau = 194   nlat_gau =  92
!     lsave =   22213  lwork =   53347
!
!     trssph output:
!     ier =  0  lsvmin =   22213  lwkmin =   53347
!     *** 32 BIT ARITHMETIC
!     least squares error =  0.201E-06
!     *** 64 BIT ARITHMETIC
!     least squares error =  0.763E-11
!
!
!     GAUSSIAN TO EQUALLY SPACED GRID TRANSFER
!
!     trssph input parameters:
!     intl =  0
!     igrid_gau(1) =  2   igrid_gau(2) =  0
!     nlon_gau = 194   nlat_gau =  92
!     igrid_reg(1) = -1   igrid_reg(2) =  1
!     nlon_reg =  36   nlat_reg =  19
!     lsave =   22213  lwork =   53347
!
!     trssph output:
!     ier =  0  lsvmin =   22213  lwkmin =   53347
!     *** 32 BIT ARITHMETIC
!     least squares error =  0.618E-06
!     *** 64 BIT ARITHMETIC
!     least squares error =  0.547E-11
!
! **********************************************************************
!
!     END OF OUTPUT ... CODE FOLLOWS
!
! **********************************************************************
!
program test_trssph

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack

    ! Explicit typing only
    implicit none

    !Set grid sizes with parameter statements
    integer(ip), parameter :: nlat_gau = 92, nlon_gau = 194
    integer(ip), parameter :: nlat_reg = 19, nlon_reg = 36
    real(wp) :: data_reg(nlat_reg, nlon_reg), data_gau(nlon_gau, nlat_gau)
    real(wp), dimension(nlat_gau) :: colatitudes, gaussian_latitudes, gaussian_weights
    integer(ip) :: igrid_reg(2), igrid_gau(2)
    real(wp)    :: dlat_reg, dlon_reg, dlon_gau
    real(wp)    :: cosp, sinp, cost, sint, xyz, err2, theta, phi, dif
    integer(ip) :: i, j, intl, error_flag

    ! Set equally spaced grid increments
    dlat_reg = PI/(nlat_reg-1)
    dlon_reg = TWO_PI/nlon_reg
    dlon_gau = TWO_PI/nlon_gau

    ! Set given data in data_reg from f(x, y, z)= exp(x*y*z) restricted
    ! to nlat_reg by nlon_reg equally spaced grid on the sphere
    do  j=1, nlon_reg
        phi = real(j-1, kind=wp)*dlon_reg
        cosp = cos(phi)
        sinp = sin(phi)
        do i=1, nlat_reg
            ! Set north to south oriented colatitude point
            theta = real(i-1, kind=wp)*dlat_reg
            cost = cos(theta)
            sint = sin(theta)
            xyz = (sint*(sint*cost*sinp*cosp))
            data_reg(i, j) = exp(xyz)
        end do
    end do

    ! Set initial call flag
    intl = 0

    ! Flag data_reg grid as north to south equally spaced
    igrid_reg(1) = -1

    ! Flag data_reg grid as colatitude by longitude
    igrid_reg(2) = 1

    ! Flag data_gau grid as south to north Gaussian
    igrid_gau(1) = 2

    ! Flag data_gau grid as longitude by latitude
    igrid_gau(2) = 0

    ! Print trssph input parameters
    write (stdout, 100) intl, igrid_reg(1), igrid_reg(2), nlon_reg, nlat_reg, &
        igrid_gau(1), igrid_gau(2), nlon_gau, nlat_gau
100 format(//' EQUALLY SPACED TO GAUSSIAN GRID TRANSFER ' , &
        /' trssph input arguments: ' , &
        /' intl = ', i2, &
        /' igrid_reg(1) = ', i2, 2x, ' igrid_reg(2) = ', i2, &
        /' nlon_reg = ', i3, 2x, ' nlat_reg = ', i3, &
        /' igrid_gau(1) = ', i2, 2x, ' igrid_gau(2) = ', i2, &
        /' nlon_gau = ', i3, 2x, ' nlat_gau = ', i3)

    ! Transfer data from data_reg to data_gau
    call trssph(intl, igrid_reg, nlon_reg, nlat_reg, data_reg, igrid_gau, nlon_gau, &
        nlat_gau, data_gau, error_flag)

    if (error_flag == 0) then
        !
        ! Compute nlat_gau gaussian colatitudinal points
        ! and set in colatitudes with south to north orientation
        ! for computing error in data_gau
        call compute_gaussian_latitudes_and_weights( &
            nlat_gau, gaussian_latitudes, gaussian_weights, error_flag)

        colatitudes = PI - gaussian_latitudes

        ! Compute the least squares error in data_gau
        err2 = 0.0_wp
        do j=1, nlon_gau
            phi = real(j - 1, kind=wp) * dlon_gau
            cosp = cos(phi)
            sinp = sin(phi)
            do i=1, nlat_gau
                theta = colatitudes(i)
                cost = cos(theta)
                sint = sin(theta)
                xyz = (sint*(sint*cost*sinp*cosp))
                dif = abs(data_gau(j, i)-exp(xyz))
                err2 = err2 + dif**2
            end do
        end do
        err2 = sqrt(err2/(nlon_gau*nlat_gau))
        write (stdout, 300) err2
300     format(' least squares error = ', e10.3)
    end if

    ! Set data_reg to zero
    data_reg = 0.0_wp

    write (stdout, 400) intl, igrid_gau(1), igrid_gau(2), nlon_gau, nlat_gau, igrid_reg(1), &
        igrid_reg(2), nlon_reg, nlat_reg
400 format(/' GAUSSIAN TO EQUALLY SPACED GRID TRANSFER ' , &
        /' trssph input arguments: ' , &
        /' intl = ', i2, &
        /' igrid_gau(1) = ', i2, 2x, ' igrid_gau(2) = ', i2, &
        /' nlon_gau = ', i3, 2x, ' nlat_gau = ', i3, &
        /' igrid_reg(1) = ', i2, 2x, ' igrid_reg(2) = ', i2, &
        /' nlon_reg = ', i3, 2x, ' nlat_reg = ', i3)

    ! Transfer data_gau back to data_reg
    call trssph(intl, igrid_gau, nlon_gau, nlat_gau, data_gau, igrid_reg, nlon_reg, &
        nlat_reg, data_reg, error_flag)

    if (error_flag == 0) then

        ! Compute the least squares error in data_reg
        err2 = 0.0
        do j=1, nlon_reg
            phi = real(j-1, kind=wp)*dlon_reg
            cosp = cos(phi)
            sinp = sin(phi)
            do i=1, nlat_reg
                theta = real(i - 1, kind=wp)*dlat_reg
                cost = cos(theta)
                sint = sin(theta)
                xyz = (sint*(sint*cost*sinp*cosp))
                dif = abs(data_reg(i, j)-exp(xyz))
                err2 = err2+dif**2
            end do
        end do
        err2 = sqrt(err2/(nlat_reg*nlon_reg))
        write (stdout, 300) err2
    end if

end program test_trssph
