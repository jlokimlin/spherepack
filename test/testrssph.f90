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
!     longitude array DATAE which runs north to south with increasing subscript
!     values.  The Gaussian grid data is stored in a 192 X 94 longitude-latitude
!     array DATAG which runs south to north with increasing latitude subscript.
!     First trssph is used to transfer DATAE to DATAG.  Then trssph is used to
!     transfer DATAG back to DATAE.
!
!     For testing purposes, DATAE is set equal the analytic function
!
!                       x*y*z
!           f(x,y,z) = e
!
!     in Cartesian coordinates x,y,z restricted to the surface of the sphere.
!     The same function is used to compute error in DATAG after the data transfer
!     with trssph.  Finally this is used to compute error in DATAE after the transfer
!     back with trssph.  Output from executing the test program on machines with
!     32 bit and 64 bit arithmetic is listed below.  The minimum required saved
!     and unsaved work space lengths were predetermined by an earlier call to
!     trssph with nlone=36,nlate=19,nlong=192,nlatg=94,lsave=1,lwork=1 and printout
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
!     igride(1) = -1   igride(2) =  1
!     nlone =  36   nlate =  19
!     igridg(1) =  2   igridg(2) =  0
!     nlong = 194   nlatg =  92
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
!     igridg(1) =  2   igridg(2) =  0
!     nlong = 194   nlatg =  92
!     igride(1) = -1   igride(2) =  1
!     nlone =  36   nlate =  19
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
program testrssph
    implicit none
    !
    !     set grid sizes with parameter statements
    !
    integer :: nnlatg
    integer :: nnlong
    integer :: nnlate
    integer :: nnlone
    integer :: llwork
    integer :: llsave
    integer :: lldwork
    parameter (nnlatg=92, nnlong=194, nnlate=19,nnlone=36)
    !
    !     set predetermined minimum saved and unsaved work space lengths
    !
    parameter (llwork = 53347, llsave = 22213)
    !
    !     set real work space lengt (see shagci.f documentation)
    !
    parameter (lldwork = nnlatg*(nnlatg+4))
    !
    !     dimension and type data arrays and grid vectors and internal variables
    !
    real :: regular_data(nnlate,nnlone)
    real :: gaussian_data(nnlong,nnlatg)
    real :: work(llwork)
    real :: wsave(llsave)
    real :: thetag(nnlatg)
    real :: dwork(lldwork)
    real :: dtheta(nnlatg)
    real :: dwts(nnlatg)
    integer :: igride(2)
    integer :: igridg(2)
    integer :: GAUSSIAN_NLATS
    integer :: GAUSSIAN_NLONS
    integer :: REGULAR_NLATS
    integer :: REGULAR_NLONS
    integer :: lwork
    integer :: lsave
    integer :: ldwork
    real :: pi
    real :: dlate
    real :: dlone
    real :: dlong
    real :: cp
    real :: sp
    real :: ct
    real :: st
    real :: xyz
    real :: err2
    real :: t
    real :: p
    real :: dif
    integer :: i
    integer :: j
    integer :: intl
    integer :: ier
    integer :: lsvmin
    integer :: lwkmin
    !
    !     set grid sizes and dimensions from parameter statements
    !
    GAUSSIAN_NLATS = nnlatg
    GAUSSIAN_NLONS = nnlong
    REGULAR_NLATS = nnlate
    REGULAR_NLONS = nnlone
    lwork = llwork
    ldwork = lldwork
    lsave = llsave
    !
    !     set equally spaced grid increments
    !
    pi = acos(-1.0)
    dlate = pi/(REGULAR_NLATS-1)
    dlone = (pi+pi)/REGULAR_NLONS
    dlong = (pi+pi)/GAUSSIAN_NLONS
    !
    !     set given data in DATAE from f(x,y,z)= exp(x*y*z) restricted
    !     to nlate by nlone equally spaced grid on the sphere
    !
    do  j=1,REGULAR_NLONS
        p = (j-1)*dlone
        cp = cos(p)
        sp = sin(p)
        do i=1,REGULAR_NLATS
            !
            !     set north to south oriented colatitude point
            !
            t = (i-1)*dlate
            ct = cos(t)
            st = sin(t)
            xyz = (st*(st*ct*sp*cp))
            regular_data(i,j) = exp(xyz)
        end do
    end do
    !
    !     set initial call flag
    !
    intl = 0
    !
    !     flag DATAE grid as north to south equally spaced
    !
    igride(1) = -1
    !
    !     flag DATAE grid as colatitude by longitude
    !
    igride(2) = 1
    !
    !     flag DATAG grid as south to north Gaussian
    !
    igridg(1) = 2
    !
    !     flag DATAG grid as longitude by latitude
    !
    igridg(2) = 0
    !
    !     print trssph input parameters
    !
    write(*,100) intl,igride(1),igride(2),REGULAR_NLONS,REGULAR_NLATS, &
        igridg(1),igridg(2),GAUSSIAN_NLONS,GAUSSIAN_NLATS,lsave,lwork,ldwork
100 format(//' EQUALLY SPACED TO GAUSSIAN GRID TRANSFER ' , &
        /' trssph input arguments: ' , &
        /' intl = ',i2, &
        /' igride(1) = ',i2,2x,' igride(2) = ',i2, &
        /' nlone = ',i3,2x,' nlate = ',i3, &
        /' igridg(1) = ',i2,2x,' igridg(2) = ',i2, &
        /' nlong = ',i3,2x,' nlatg = ',i3, &
        /' lsave = ',i7,2x,' lwork = ',i7,2x,' ldwork = ',i5)
    !
    !     transfer data from DATAE to DATAG
    !
    call trssph(intl,igride,REGULAR_NLONS,REGULAR_NLATS,regular_data,igridg,GAUSSIAN_NLONS, &
        GAUSSIAN_NLATS,gaussian_data,wsave,lsave,lsvmin,work,lwork,lwkmin,dwork, &
        ldwork,ier)
    !
    !     print output parameters
    !
    write (*,200) ier, lsvmin, lwkmin
200 format(//' trssph output: ' &
        / ' ier = ', i2,2x, 'lsvmin = ',i7, 2x,'lwkmin = ',i7)
    if (ier == 0) then
        !
        !     compute nlatg gaussian colatitude points using spherepack routine "gaqd"
        !     and set in single precision vector thetag with south to north orientation
        !     for computing error in DATAG
        !
        call gaqd(GAUSSIAN_NLATS,dtheta,dwts,dwork,ldwork,ier)
        do  i=1,GAUSSIAN_NLATS
            thetag(i) = pi-dtheta(i)
        end do
        !
        !     compute the least squares error in DATAG
        !
        err2 = 0.0
        do j=1,GAUSSIAN_NLONS
            p = (j-1)*dlong
            cp = cos(p)
            sp = sin(p)
            do i=1,GAUSSIAN_NLATS
                t = thetag(i)
                ct = cos(t)
                st = sin(t)
                xyz = (st*(st*ct*sp*cp))
                dif = abs(gaussian_data(j,i)-exp(xyz))
                err2 = err2+dif*dif
            end do
        end do
        err2 = sqrt(err2/size(gaussian_data))
        write (6,300) err2
300     format(' least squares error = ',e10.3)
    end if
    !
    !     set DATAE to zero
    !
    do j=1,REGULAR_NLONS
        do i=1,REGULAR_NLATS
            regular_data(i,j) = 0.0
        end do
    end do

    write(*,400) intl,igridg(1),igridg(2),GAUSSIAN_NLONS,GAUSSIAN_NLATS,igride(1), &
        igride(2),REGULAR_NLONS,REGULAR_NLATS,lsave,lwork,ldwork
400 format(/' GAUSSIAN TO EQUALLY SPACED GRID TRANSFER ' , &
        /' trssph input arguments: ' , &
        /' intl = ',i2, &
        /' igridg(1) = ',i2,2x,' igridg(2) = ',i2, &
        /' nlong = ',i3,2x,' nlatg = ',i3, &
        /' igride(1) = ',i2,2x,' igride(2) = ',i2, &
        /' nlone = ',i3,2x,' nlate = ',i3, &
        /' lsave = ',i7,2x,'lwork = ',i7, 2x, 'ldwork = ',i7)
    !
    !     transfer DATAG back to DATAE
    !
    call trssph(intl,igridg,GAUSSIAN_NLONS,GAUSSIAN_NLATS,gaussian_data,igride,REGULAR_NLONS, &
        REGULAR_NLATS,regular_data,wsave,lsave,lsvmin,work,lwork,lwkmin,dwork, &
        ldwork,ier)
    !
    !     print output parameters
    !
    write (*,200) ier, lsvmin, lwkmin
    if (ier == 0) then
        !
        !     compute the least squares error in DATAE
        !
        err2 = 0.0
        do j=1,REGULAR_NLONS
            p = (j-1)*dlone
            cp = cos(p)
            sp = sin(p)
            do i=1,REGULAR_NLATS
                t = (i-1)*dlate
                ct = cos(t)
                st = sin(t)
                xyz = (st*(st*ct*sp*cp))
                dif = abs(regular_data(i,j)-exp(xyz))
                err2 = err2+dif*dif
            end do
        end do
        err2 = sqrt(err2/size(regular_data))
        write (6,300) err2
    end if

end program testrssph
