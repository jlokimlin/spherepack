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
!
! ... file shaes.f
!
!     this file contains code and documentation for subroutines
!     shaes and shaesi
!
! ... files which must be loaded with shaes.f
!
!     sphcom.f, hrfft.f
!
!     subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, 
!    +                 wshaes, lshaes, work, lwork, ierror)
!
!     subroutine shaes performs the spherical harmonic analysis
!     on the array g and stores the result in the arrays a and b.
!     the analysis is performed on an equally spaced grid.  the
!     associated legendre functions are stored rather than recomputed
!     as they are in subroutine shaec.  the analysis is described
!     below at output parameters a, b.
!
!     sphcom.f, hrfft.f
!
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     isym   = 0  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
!                 (see description of g below)
!
!            = 1  g is antisymmetric about the equator. the analysis
!                 is performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the analysis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the analysis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!            = 2  g is symmetric about the equator. the analysis is
!                 performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the analysis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the analysis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!     nt     the number of analyses.  in the program that calls shaes, 
!            the arrays g, a and b can be three dimensional in which
!            case multiple analyses will be performed.  the third
!            index is the analysis index which assumes the values
!            k=1, ..., nt.  for a single analysis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that the arrays g, a and b
!            have only two dimensions.
!
!     g      a two or three dimensional array (see input parameter
!            nt) that contains the discrete function to be analyzed.
!            g(i, j) contains the value of the function at the colatitude
!            point theta(i) = (i-1)*pi/(nlat-1) and longitude point
!            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined
!            above at the input parameter isym.
!
!
!     idg    the first dimension of the array g as it appears in the
!            program that calls shaes.  if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg
!            must be at least nlat/2 if nlat is even or at least
!            (nlat+1)/2 if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shaes.  jdg must be at least nlon.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shaes. mdab must be at least
!            min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shaes. ndab must be at least nlat
!
!     wshaes an array which must be initialized by subroutine shaesi.
!            once initialized, wshaes can be used repeatedly by shaes
!            as long as nlon and nlat remain unchanged.  wshaes must
!            not be altered between calls of shaes.
!
!     lshaes the dimension of the array wshaes as it appears in the
!            program that calls shaes. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshaes must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shaes.  define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if isym is zero then lwork must be at least
!            (nt+1)*nlat*nlon. if isym is not zero then
!            lwork must be at least (nt+1)*l2*nlon.
!
!
!     **************************************************************
!
!     output parameters
!
!     a, b    both a, b are two or three dimensional arrays (see input
!            parameter nt) that contain the spherical harmonic
!            coefficients in the representation of g(i, j) given in the
!            discription of subroutine shses. for isym=0, a(m, n) and
!            b(m, n) are given by the equations listed below. symmetric
!            versions are used when isym is greater than zero.
!
!
!
!     definitions
!
!     1. the normalized associated legendre functions
!
!     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                       *sin(theta)**m/(2**n*factorial(n)) times the
!                       (n+m)th derivative of (x**2-1)**n with respect
!                       to x=cos(theta)
!
!     2. the normalized z functions for m even
!
!     zbar(m, n, theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
!                       the integral from tau = 0 to tau = pi of
!                       cos(k*theta)*cos(k*tau)*pbar(m, n, tau)*sin(tau)
!                       (first and last terms in this sum are divided
!                       by 2)
!
!     3. the normalized z functions for m odd
!
!     zbar(m, n, theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
!                       of the integral from tau = 0 to tau = pi of
!                       sin(k*theta)*sin(k*tau)*pbar(m, n, tau)*sin(tau)
!
!     4. the fourier transform of g(i, j).
!
!     c(m, i)          = 2/nlon times the sum from j=1 to j=nlon
!                       of g(i, j)*cos((m-1)*(j-1)*2*pi/nlon)
!                       (the first and last terms in this sum
!                       are divided by 2)
!
!     s(m, i)          = 2/nlon times the sum from j=2 to j=nlon
!                       of g(i, j)*sin((m-1)*(j-1)*2*pi/nlon)
!
!     5. the maximum (plus one) longitudinal wave number
!
!            mmax = min(nlat, (nlon+2)/2) if nlon is even or
!            mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!     then for m=0, ..., mmax-1 and  n=m, ..., nlat-1  the arrays a, b are
!     given by
!
!     a(m+1, n+1)      = the sum from i=1 to i=nlat of
!                       c(m+1, i)*zbar(m, n, theta(i))
!                       (first and last terms in this sum are
!                       divided by 2)
!
!     b(m+1, n+1)      = the sum from i=1 to i=nlat of
!                       s(m+1, i)*zbar(m, n, theta(i))
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of isym
!            = 4  error in the specification of nt
!            = 5  error in the specification of idg
!            = 6  error in the specification of jdg
!            = 7  error in the specification of mdab
!            = 8  error in the specification of ndab
!            = 9  error in the specification of lshaes
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!     subroutine shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, 
!    +                  ldwork, ierror)
!
!     subroutine shaesi initializes the array wshaes which can then
!     be used repeatedly by subroutine shaes
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     lshaes the dimension of the array wshaes as it appears in the
!            program that calls shaesi. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshaes must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a real   work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shaesi.  define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lwork must be at least
!
!               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
!
!
!     dwork  a real work array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls shaesi.  ldwork must be at least nlat+1
!
!
!     output parameters
!
!     wshaes an array which is initialized for use by subroutine shaes.
!            once initialized, wshaes can be used repeatedly by shaes
!            as long as nlon and nlat remain unchanged.  wshaes must
!            not be altered between calls of shaes.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshaes
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
!
!
subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
    wshaes, lshaes, work, lwork, ierror)
    implicit none
    ! External routines: shaes1
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: nlat
    integer, intent (in)     :: nlon
    integer, intent (in)     :: isym
    integer, intent (in)     :: nt
    real,    intent (in)     :: g(idg, jdg, 1)
    integer, intent (in)     :: idg
    integer, intent (in)     :: jdg
    real,    intent (out)    :: a(mdab, ndab, 1)
    real,    intent (out)    :: b(mdab, ndab, 1)
    integer, intent (in)     :: mdab
    integer, intent (in)     :: ndab
    real,    intent (in out) :: wshaes(1)
    integer, intent (in)     :: lshaes
    real,    intent (in out) :: work(1)
    integer, intent (in)     :: lwork
    integer, intent (out)    :: ierror
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer :: ist, mmax, imid, idz, lzimn, ls, nln
    !----------------------------------------------------------------------

    !
    !==> Set constants
    mmax = min(nlat, nlon/2+1)
    imid = (nlat+1)/2
    idz = (mmax*(nlat+nlat-mmax+1))/2
    lzimn = idz*imid
    ls = nlat

    if (isym > 0) then
        ls = imid
    end if

    nln = nt*ls*nlon

    !
    !==> Check validity of input arguments
    !

    ! Initialize error flag
    ierror = 0

    ! Check case 1
    if (nlat < 3) then
        ierror = 1
        return
    end if

    ! Check case 2
    if (nlon < 4) then
        ierror = 2
        return
    end if

    ! Check case 3
    if (isym < 0 .or. isym > 2) then
        ierror = 3
        return
    end if

    ! Check case 4
    if (nt < 0) then
        ierror = 4
        return
    end if

    ! Check case 5
    if ((isym == 0 .and. idg < nlat) .or. &
        (isym /= 0 .and. idg < (nlat+1)/2)) then
        ierror = 5
        return
    end if

    ! Check case 6
    if (jdg < nlon) then
        ierror = 6
        return
    end if


    ! Check case 7
    if (mdab < mmax) then
        ierror = 7
        return
    end if

    ! Check case 8
    if (ndab < nlat) then
        ierror = 8
        return
    end if

    ! Check case 9
    if (lshaes < lzimn+nlon+15) then
        ierror = 9
        return
    end if

    ! Check case 10
    if (lwork < nln+ls*nlon) then
        ierror = 10
        return
    end if

    ! Set calling argument for analsyis
    select case (isym)
        case (0)
            ist = imid
        case default
            ist = 0
    end select


    call shaes1(nlat, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshaes, idz, &
        ls, nlon, work, work(ist+1), work(nln+1), wshaes(lzimn+1))

end subroutine shaes



subroutine shaes1(nlat, isym, nt, g, idgs, jdgs, a, b, mdab, ndab, &
    z, idz, idg, jdg, ge, go, work, whrfft)
    implicit none
    ! External routines: hrfftf
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: nlat
    integer, intent (in)     :: isym
    integer, intent (in)     :: nt
    real,    intent (in)     :: g(idgs, jdgs, 1)
    integer, intent (in)     :: idgs
    integer, intent (in)     :: jdgs
    real,    intent (in out) :: a(mdab, ndab, 1)
    real,    intent (in out) :: b(mdab, ndab, 1)
    integer, intent (in)     :: mdab
    integer, intent (in)     :: ndab
    real,    intent (in out) :: z(idz, 1)
    integer, intent (in)     :: idz
    integer, intent (in)     :: idg
    integer, intent (in)     :: jdg
    real,    intent (in out) :: ge(idg, jdg, 1)
    real,    intent (in out) :: go(idg, jdg, 1)
    real,    intent (in out) :: work(1)
    real,    intent (in out) :: whrfft(1)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer :: i, j, k, m, mb, ls, mp1, np1, mp2, mdo, ndo
    real    :: fsn, tsn
    integer :: imm1, nlp1, imid, modl, mmax, nlon
    !----------------------------------------------------------------------


    ls = idg
    nlon = jdg
    mmax = min(nlat, nlon/2+1)
    mdo = mmax

    if (mdo+mdo-1 > nlon) then
        mdo = mmax-1
    end if

    nlp1 = nlat+1
    tsn = 2.0/nlon
    fsn = 4.0/nlon
    imid = (nlat+1)/2
    modl = mod(nlat, 2)
    imm1 = imid

    if (modl /= 0) then
        imm1 = imid-1
    end if

    if (isym /= 0) then
        go to 15
    end if

    do k=1, nt
        do i=1, imm1
            do j=1, nlon
                ge(i, j, k) = tsn*(g(i, j, k)+g(nlp1-i, j, k))
                go(i, j, k) = tsn*(g(i, j, k)-g(nlp1-i, j, k))
            end do
        end do
    end do
    go to 30
    15 do k=1, nt
        do i=1, imm1
            do j=1, nlon
                ge(i, j, k) = fsn*g(i, j, k)
            end do
        end do
    end do

    if (isym == 1) then
        go to 27
    end if

30  if (modl == 0) then
        go to 27
    end if

    do k=1, nt
        do j=1, nlon
            ge(imid, j, k) = tsn*g(imid, j, k)
        end do
    end do

    27 do k=1, nt
        call hrfftf(ls, nlon, ge(1, 1, k), ls, whrfft, work)
        if (mod(nlon, 2) /= 0) cycle
        do i=1, ls
            ge(i, nlon, k) = 0.5*ge(i, nlon, k)
        end do
    end do

    do k=1, nt
        do mp1=1, mmax
            do np1=mp1, nlat
                a(mp1, np1, k) = 0.0
                b(mp1, np1, k) = 0.0
            end do
        end do
    end do

    if (isym == 1) then
        go to 145
    end if

    do k=1, nt
        do i=1, imid
            do np1=1, nlat, 2
                a(1, np1, k) = a(1, np1, k)+z(np1, i)*ge(i, 1, k)
            end do
        end do
    end do

    ndo = nlat

    if (mod(nlat, 2) == 0) then
        ndo = nlat-1
    end if

    do mp1=2, mdo
        m = mp1-1
        mb = m*(nlat-1)-(m*(m-1))/2
        do k=1, nt
            do i=1, imid
                do np1=mp1, ndo, 2
                    a(mp1, np1, k) = a(mp1, np1, k)+z(np1+mb, i)*ge(i, 2*mp1-2, k)
                    b(mp1, np1, k) = b(mp1, np1, k)+z(np1+mb, i)*ge(i, 2*mp1-1, k)
                end do
            end do
        end do
    end do

    if (mdo == mmax .or. mmax > ndo) then
        go to 135
    end if

    mb = mdo*(nlat-1)-(mdo*(mdo-1))/2

    do k=1, nt
        do i=1, imid
            do np1=mmax, ndo, 2
                a(mmax, np1, k) = a(mmax, np1, k)+z(np1+mb, i)*ge(i, 2*mmax-2, k)
            end do
        end do
    end do

135 if (isym == 2) then
        return
    end if

    145 do k=1, nt
        do i=1, imm1
            do np1=2, nlat, 2
                a(1, np1, k) = a(1, np1, k)+z(np1, i)*go(i, 1, k)
            end do
        end do
    end do

    ndo = nlat

    if (mod(nlat, 2) /= 0) then
        ndo = nlat-1
    end if

    do mp1=2, mdo
        m = mp1-1
        mp2 = mp1+1
        mb = m*(nlat-1)-(m*(m-1))/2
        do k=1, nt
            do i=1, imm1
                do np1=mp2, ndo, 2
                    a(mp1, np1, k) = a(mp1, np1, k)+z(np1+mb, i)*go(i, 2*mp1-2, k)
                    b(mp1, np1, k) = b(mp1, np1, k)+z(np1+mb, i)*go(i, 2*mp1-1, k)
                end do
            end do
        end do
    end do

    mp2 = mmax+1

    if (mdo == mmax .or. mp2 > ndo) then
        return
    end if

    mb = mdo*(nlat-1)-(mdo*(mdo-1))/2

    do k=1, nt
        do i=1, imm1
            do np1=mp2, ndo, 2
                a(mmax, np1, k) = a(mmax, np1, k)+z(np1+mb, i)*go(i, 2*mmax-2, k)
            end do
        end do
    end do

end subroutine shaes1



subroutine shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, &
    ldwork, ierror)
    implicit none
    ! External routines: sea1, hrffti
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: nlat
    integer, intent (in)     :: nlon
    real,    intent (in out) :: wshaes(*)
    integer, intent (in)     :: lshaes
    real,    intent (in)     :: work(*)
    integer, intent (in)     :: lwork
    real,    intent (in out) :: dwork(*)
    integer, intent (in)     :: ldwork
    integer, intent (out)    :: ierror
    !----------------------------------------------------------------------

    !
    !==> Remarks:
    !    length of wshaes is (l*(l+1)*imid)/2+nlon+15
    !    length of work is 5*l*imid + 3*((l-3)*l+2)/2
    !

    !
    !==> Check validity of input values
    !

    ! Initialize error flag
    ierror = 0

    ! Check case 1
    if (nlat<3) then
        ierror = 1
        return
    end if

    ! Check case 2
    if (nlon<4) then
        ierror = 2
        return
    end if

    associate( &
        mmax => min(nlat, nlon/2+1), &
        imid => (nlat+1)/2 &
        )
        associate( lzimn => (imid*mmax*(nlat+nlat-mmax+1))/2 )

            ! Check case 3
            if (lshaes < lzimn+nlon+15) then
                ierror = 3
                return
            end if

            ! Check case 4
            associate( labc => 3*((mmax-2)*(nlat+nlat-mmax-1))/2 )
                if (lwork < 5*nlat*imid + labc) then
                    ierror = 4
                    return
                end if
            end associate

            ! Check case 5
            if (ldwork < nlat+1) then
                ierror = 5
                return
            end if

            !
            !==> Compute workspace
            !
            associate( &
                iw1 => 3*nlat*imid+1, &
                idz => (mmax*(nlat+nlat-mmax+1))/2 &
                )
                call sea1(nlat, nlon, imid, wshaes, idz, work, work(iw1), dwork)
                call hrffti(nlon, wshaes(lzimn+1))
            end associate
        end associate
    end associate

end subroutine shaesi
