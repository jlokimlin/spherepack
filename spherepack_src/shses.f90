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
! ... file shses.f
!
!     this file contains code and documentation for subroutines
!     shses and shsesi
!
! ... files which must be loaded with shses.f
!
!     sphcom.f, hrfft.f
!
!     subroutine shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, 
!    +                 wshses, lshses, work, lwork, ierror)
!
!     subroutine shses performs the spherical harmonic synthesis
!     on the arrays a and b and stores the result in the array g.
!     the synthesis is performed on an equally spaced grid.  the
!     associated legendre functions are stored rather than recomputed
!     as they are in subroutine shsec.  the synthesis is described
!     below at output parameter g.
!
! *** required files from spherepack2
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
!     isym   = 0  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
!                 (see description of g below)
!
!            = 1  g is antisymmetric about the equator. the synthesis
!                 is performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the synthesis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the synthesis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!            = 2  g is symmetric about the equator. the synthesis is
!                 performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the synthesis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the synthesis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!     nt     the number of syntheses.  in the program that calls shses, 
!            the arrays g, a and b can be three dimensional in which
!            case multiple syntheses will be performed.  the third
!            index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that the arrays g, a and b
!            have only two dimensions.
!
!     idg    the first dimension of the array g as it appears in the
!            program that calls shses.  if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg
!            must be at least nlat/2 if nlat is even or at least
!            (nlat+1)/2 if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shses.  jdg must be at least nlon.
!
!     a, b    two or three dimensional arrays (see the input parameter
!            nt) that contain the coefficients in the spherical harmonic
!            expansion of g(i, j) given below at the definition of the
!            output parameter g.  a(m, n) and b(m, n) are defined for
!            indices m=1, ..., mmax and n=m, ..., nlat where mmax is the
!            maximum (plus one) longitudinal wave number given by
!            mmax = min(nlat, (nlon+2)/2) if nlon is even or
!            mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shses. mdab must be at least
!            min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shses. ndab must be at least nlat
!
!     wshses an array which must be initialized by subroutine shsesi.
!            once initialized, wshses can be used repeatedly by shses
!            as long as nlon and nlat remain unchanged.  wshses must
!            not be altered between calls of shses.
!
!     lshses the dimension of the array wshses as it appears in the
!            program that calls shses. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshses must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shses.  define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if isym is zero then lwork must be at least
!
!               (nt+1)*nlat*nlon
!
!            if isym is nonzero lwork must be at least
!
!               (nt+1)*l2*nlon.
!
!     **************************************************************
!
!     output parameters
!
!     g      a two or three dimensional array (see input parameter
!            nt) that contains the spherical harmonic synthesis of
!            the arrays a and b at the colatitude point theta(i) =
!            (i-1)*pi/(nlat-1) and longitude point phi(j) =
!            (j-1)*2*pi/nlon. the index ranges are defined above at
!            at the input parameter isym.  for isym=0, g(i, j) is
!            given by the the equations listed below.  symmetric
!            versions are used when isym is greater than zero.
!
!     the normalized associated legendre functions are given by
!
!     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                       *sin(theta)**m/(2**n*factorial(n)) times the
!                       (n+m)th derivative of (x**2-1)**n with respect
!                       to x=cos(theta)
!
!     define the maximum (plus one) longitudinal wave number
!     as   mmax = min(nlat, (nlon+2)/2) if nlon is even or
!          mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!     then g(i, j) = the sum from n=0 to n=nlat-1 of
!
!                   .5*pbar(0, n, theta(i))*a(1, n+1)
!
!              plus the sum from m=1 to m=mmax-1 of
!
!                   the sum from n=m to n=nlat-1 of
!
!              pbar(m, n, theta(i))*(a(m+1, n+1)*cos(m*phi(j))
!                                    -b(m+1, n+1)*sin(m*phi(j)))
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
!            = 9  error in the specification of lshses
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!     subroutine shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, 
!    +                  ldwork, ierror)
!
!     subroutine shsesi initializes the array wshses which can then
!     be used repeatedly by subroutine shses.
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
!     lshses the dimension of the array wshses as it appears in the
!            program that calls shsesi. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshses must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a real   work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in
!            the program that calls shsesi.  define
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
!            program that calls shsesi.  ldwork must be at least nlat+1
!
!
!     output parameters
!
!     wshses an array which is initialized for use by subroutine shses.
!            once initialized, wshses can be used repeatedly by shses
!            as long as nlon and nlat remain unchanged.  wshses must
!            not be altered between calls of shses.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshses
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
subroutine shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
    wshses, lshses, work, lwork, ierror)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: nlat
    integer (ip), intent (in)     :: nlon
    integer (ip), intent (in)     :: isym
    integer (ip), intent (in)     :: nt
    real (wp),    intent (out)    :: g(idg, jdg, 1)
    integer (ip), intent (in)     :: idg
    integer (ip), intent (in)     :: jdg
    real (wp),    intent (in)     :: a(mdab, ndab, 1)
    real (wp),    intent (in)     :: b(mdab, ndab, 1)
    integer (ip), intent (in)     :: mdab
    integer (ip), intent (in)     :: ndab
    real (wp),    intent (in out) :: wshses(1)
    integer (ip), intent (in)     :: lshses
    real (wp),    intent (in out) :: work(1)
    integer (ip), intent (in)     :: lwork
    integer (ip), intent (out)    :: ierror
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip) :: ls, nln, ist, imid, mmax, lpimn
    !----------------------------------------------------------------------

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
    if ( &
        (isym == 0 .and. idg < nlat) &
        .or. &
        (isym /= 0 .and. idg < (nlat+1)/2) &
        ) then
        ierror = 5
        return
    end if

    ! Check case 6
    if (jdg < nlon) then
        ierror = 6
        return
    end if


    mmax = min(nlat, nlon/2+1)

    ! Check case 7
    if (mdab < mmax)  then
        ierror = 7
        return
    end if

    ! Check case 8
    if (ndab < nlat) then
        ierror = 8
        return
    end if

    imid = (nlat+1)/2
    lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2

    ! Check case 9
    if (lshses < lpimn+nlon+15) then
        ierror = 9
        return
    end if

    if (isym > 0) then
        ls = imid
    else
        ls = nlat
    end if

    nln = nt*ls*nlon

    ! Check case 10
    if (lwork < nln+ls*nlon) then
        ierror = 10
        return
    end if

    if (isym == 0) then
        ist = imid
    else
        ist = 0
    end if

    call shses1(nlat, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshses, imid, &
        ls, nlon, work, work(ist+1), work(nln+1), wshses(lpimn+1))

end subroutine shses


subroutine shses1(nlat, isym, nt, g, idgs, jdgs, a, b, mdab, ndab, p, imid, &
    idg, jdg, ge, go, work, whrfft)

    ! External subroutines: hrfftb

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: nlat
    integer (ip), intent (in)     :: isym
    integer (ip), intent (in)     :: nt
    real (wp),    intent (in out) :: g(idgs, jdgs, 1)
    integer (ip), intent (in)     :: idgs
    integer (ip), intent (in)     :: jdgs
    real (wp),    intent (in)     :: a(mdab, ndab, 1)
    real (wp),    intent (in)     :: b(mdab, ndab, 1)
    integer (ip), intent (in)     :: mdab
    integer (ip), intent (in)     :: ndab
    real (wp),    intent (in out) :: p(imid, 1)
    integer (ip), intent (in)     :: imid
    integer (ip), intent (in)     :: idg
    integer (ip), intent (in)     :: jdg
    real (wp),    intent (in out) :: ge(idg, jdg, 1)
    real (wp),    intent (in out) :: go(idg, jdg, 1)
    real (wp),    intent (in out) :: work(1)
    real (wp),    intent (in out) :: whrfft(1)
    !----------------------------------------------------------------------
    ! Dictionary: local  variables
    !----------------------------------------------------------------------
    integer (ip) ::  i, j, k, m, mb, mn, ls, mp1, np1, mp2
    integer (ip) ::  mdo, ndo, imm1, nlp1, modl, mmax, nlo, nlon
    !----------------------------------------------------------------------

    ls = idg
    nlon = jdg
    mmax = min(nlat, nlon/2+1)

    if (2*mmax-1 > nlon) then
        mdo = mmax-1
    else
        mdo = mmax
    end if

    nlp1 = nlat+1
    modl = mod(nlat, 2)

    if (modl /= 0) then
        imm1 = imid-1
    else
        imm1 = imid
    end if

    ! Initialize
    ge = 0.0_wp

    if (isym /= 1) then
        do k = 1,nt
            do np1=1, nlat, 2
                ge(1:imid, 1, k) = &
                    ge(1:imid, 1, k) + &
                    a(1, np1, k) * p(1:imid, np1)
            end do
        end do

        if (mod(nlat, 2) == 0) then
            ndo = nlat-1
        else
            ndo = nlat
        end if

        do mp1=2, mdo
            m = mp1-1
            mb = m*(nlat-1)-(m*(m-1))/2
            do np1=mp1, ndo, 2
                mn = mb+np1
                do k=1, nt
                    ge(1:imid, 2*mp1-2, k) = &
                        ge(1:imid, 2*mp1-2, k)+a(mp1, np1, k)*p(1:imid, mn)
                    ge(1:imid, 2*mp1-1, k) = &
                        ge(1:imid, 2*mp1-1, k)+b(mp1, np1, k)*p(1:imid, mn)
                end do
            end do
        end do

        if ( .not.(mdo == mmax .or. mmax > ndo) ) then
            if (isym == 2) then
                mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
                do k = 1, nt
                    do np1=mmax, ndo, 2
                        mn = mb+np1
                        ge(1:imid, 2*mmax-2, k) = &
                            ge(1:imid, 2*mmax-2, k) &
                            + a(mmax, np1, k) * p(1:imid, mn)
                    end do
                end do
            end if

            do k = 1, nt
                do np1=2, nlat, 2
                    go(1:imm1, 1, k)= &
                        go(1:imm1, 1, k) &
                        + a(1, np1, k)*p(1:imm1, np1)
                end do
            end do

            if (mod(nlat, 2) /= 0) then
                ndo = nlat-1
            else
                ndo = nlat
            end if

            do mp1=2, mdo
                mp2 = mp1+1
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                do np1=mp2, ndo, 2
                    mn = mb+np1
                    do k=1, nt
                        do i=1, imm1
                            go(i, 2*mp1-2, k) = &
                                go(i, 2*mp1-2, k) + a(mp1, np1, k) * p(i, mn)
                            go(i, 2*mp1-1, k) = &
                                go(i, 2*mp1-1, k) + b(mp1, np1, k) * p(i, mn)
                        end do
                    end do
                end do
            end do

            mp2 = mmax+1

            if (.not.(mdo == mmax .or. mp2 > ndo) ) then
                mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
                do np1=mp2, ndo, 2
                    mn = mb+np1
                    do k=1, nt
                        go(1:imm1, 2*mmax-2, k) = &
                            go(1:imm1, 2*mmax-2, k) + &
                            a(mmax, np1, k) * p(1:imm1, mn)
                    end do
                end do
            end if
        end if
    end if

    do k=1, nt
        if (mod(nlon, 2) == 0) then
            ge(1:ls, nlon, k) = 2.0_wp * ge(1:ls, nlon, k)
        end if
        call hrfftb(ls, nlon, ge(1, 1, k), ls, whrfft, work)
    end do

    if (isym == 0) then
        do k=1, nt
            do j=1, nlon
                do i=1, imm1
                    g(i, j, k) = 0.5_wp * (ge(i, j, k)+go(i, j, k))
                    g(nlp1-i, j, k) = 0.5_wp * (ge(i, j, k)-go(i, j, k))
                end do
                if (modl == 0) then
                    exit
                end if
                g(imid, j, k) = 0.5_wp * ge(imid, j, k)
            end do
        end do
    else
        do k=1, nt
            do j=1, nlon
                g(1:imid, j, k) = 0.5_wp * ge(1:imid, j, k)
            end do
        end do
    end if

end subroutine shses1



subroutine shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, &
    ldwork, ierror)
    ! External subroutines :: ses1, hrffti
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: nlat
    integer (ip), intent (in)     :: nlon
    real (wp),    intent (in out) :: wshses(*)
    integer (ip), intent (in)     :: lshses
    real (wp),    intent (in out) :: work(*)
    integer (ip), intent (in)     :: lwork
    real (wp),    intent (in out) :: dwork(*)
    integer (wp), intent (in)     :: ldwork
    integer (ip), intent (out)    :: ierror
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip) :: mmax, imid, lpimn, labc, iw1
    !----------------------------------------------------------------------

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

    mmax = min(nlat, nlon/2+1)
    imid = (nlat+1)/2
    lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2

    ! Check case 3
    if (lshses < lpimn+nlon+15) then
        ierror = 3
        return
    end if

    labc = 3*((mmax-2)*(nlat+nlat-mmax-1))/2

    ! Check case 4
    if (lwork < 5*nlat*imid + labc) then
        ierror = 4
        return
    end if

    ! Check case 5
    if (ldwork < nlat+1) then
        ierror = 5
        return
    end if

    iw1 = 3*nlat*imid+1

    call ses1(nlat, nlon, imid, wshses, work, work(iw1), dwork)

    call hrffti(nlon, wshses(lpimn+1))

end subroutine shsesi
