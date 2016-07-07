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
! ... file shaec.f
!
!     this file contains code and documentation for subroutines
!     shaec and shaeci
!
! ... files which must be loaded with shaec.f
!
!     type_SpherepackAux.f, type_HFFTpack.f
!
!     subroutine shaec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, 
!    +                 wshaec, lshaec, work, lwork, ierror)
!
!     subroutine shaec performs the spherical harmonic analysis
!     on the array g and stores the result in the arrays a and b.
!     the analysis is performed on an equally spaced grid.  the
!     associated legendre functions are recomputed rather than stored
!     as they are in subroutine shaes.  the analysis is described
!     below at output parameters a, b.
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
!     nt     the number of analyses.  in the program that calls shaec, 
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
!            program that calls shaec.  if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg
!            must be at least nlat/2 if nlat is even or at least
!            (nlat+1)/2 if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shaec.  jdg must be at least nlon.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shaec. mdab must be at least
!            min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shaec. ndab must be at least nlat
!
!     wshaec an array which must be initialized by subroutine shaeci.
!            once initialized, wshaec can be used repeatedly by shaec
!            as long as nlon and nlat remain unchanged.  wshaec must
!            not be altered between calls of shaec.
!
!     lshaec the dimension of the array wshaec as it appears in the
!            program that calls shaec. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshaec must be at least
!
!            2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shaec. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if isym is zero then lwork must be at least
!
!                    nlat*(nt*nlon+max(3*l2, nlon))
!
!            if isym is not zero then lwork must be at least
!
!                    l2*(nt*nlon+max(3*nlat, nlon))
!
!     **************************************************************
!
!     output parameters
!
!     a, b    both a, b are two or three dimensional arrays (see input
!            parameter nt) that contain the spherical harmonic
!            coefficients in the representation of g(i, j) given in the
!            discription of subroutine shsec. for isym=0, a(m, n) and
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
!
!     then for m=0, ..., mmax-1 and n=m, ..., nlat-1 the arrays a, b
!     are given by
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
!            = 9  error in the specification of lshaec
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!     subroutine shaeci(nlat, nlon, wshaec, lshaec, dwork, ldwork, ierror)
!
!     subroutine shaeci initializes the array wshaec which can then
!     be used repeatedly by subroutine shaec.
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
!     lshaec the dimension of the array wshaec as it appears in the
!            program that calls shaeci. the array wshaec is an output
!            parameter which is described below. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshaec must be at least
!
!            2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
!
!     dwork  a real dwork array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls shaeci.  ldwork  must be at least
!            nlat+1.
!
!
!     output parameters
!
!     wshaec an array which is initialized for use by subroutine shaec.
!            once initialized, wshaec can be used repeatedly by shaec
!            as long as nlon and nlat remain unchanged.  wshaec must
!            not be altered between calls of shaec.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshaec
!            = 4  error in the specification of ldwork
!
!
!
module module_shaec

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: shaec
    public :: shaeci

contains

    subroutine shaec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
        wshaec, lshaec, work, lwork, ierror)
        real (wp) :: a
        real (wp) :: b
        real (wp) :: g
        integer (ip) :: idg
        integer (ip) :: ierror
        integer (ip) :: imid
        integer (ip) :: ist
        integer (ip) :: isym
        integer (ip) :: iw1, jw1, jw2, jw3
        integer (ip) :: jdg
        integer (ip) :: labc
        integer (ip) :: ls
        integer (ip) :: lshaec
        integer (ip) :: lwork
        integer (ip) :: lzz1
        integer (ip) :: mdab
        integer (ip) :: mmax
        integer (ip) :: ndab
        integer (ip) :: nlat
        integer (ip) :: nln
        integer (ip) :: nlon
        integer (ip) :: nt
        real (wp) :: work
        real (wp) :: wshaec
        dimension g(idg, jdg, *), a(mdab, ndab, *), b(mdab, ndab, *), wshaec(lshaec), &
            work(lwork)

        mmax = min(nlat, nlon/2+1)
        imid = (nlat+1)/2
        lzz1 = 2*nlat*imid
        labc = 3*((mmax-2)*(2*nlat-mmax-1))/2

        if (isym > 0) then
            ls = imid
        else
            ls = nlat
        end if

        nln = nt*ls*nlon

        if (isym == 0) then
            ist = imid
        else
            ist = 0
        end if

        ! Check validity of input arguments
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 4) then
            ierror = 2
            return
        else if (isym < 0 .or. isym > 2) then
            ierror = 3
            return
        else if (nt < 0) then
            ierror = 4
            return
        else if ((isym == 0 .and. idg < nlat) &
            .or. &
            (isym /= 0 .and. idg < (nlat+1)/2)) then
            ierror = 5
            return
        else if (jdg < nlon) then
            ierror = 6
            return
        else if (mdab < mmax) then
            ierror = 7
            return
        else if (ndab < nlat) then
            ierror = 8
            return
        else if (lshaec < lzz1+labc+nlon+15) then
            ierror = 9
            return
        else if (lwork < nln+max(ls*nlon, 3*nlat*imid)) then
            ierror = 10
            return
        else
            ierror = 0
        end if

        ! Set workspace pointers
        iw1 = lzz1+labc+1
        jw1 = ist + 1
        jw2 = nln + 1
        jw3 = jw2

        call shaec1(nlat, isym, nt, g, idg, jdg, a, b, mdab, ndab, imid, ls, nlon, &
            work, work(jw1), work(jw2), work(jw3), wshaec, wshaec(iw1))


    end subroutine shaec


    subroutine shaec1(nlat, isym, nt, g, idgs, jdgs, a, b, mdab, ndab, imid, &
        idg, jdg, ge, go, work, zb, wzfin, whrfft)

        real (wp) :: a
        real (wp) :: b
        real (wp) :: fsn
        real (wp) :: g
        real (wp) :: ge
        real (wp) :: go
        integer (ip) :: i
        integer (ip) :: i3
        integer (ip) :: idg
        integer (ip) :: idgs
        integer (ip) :: imid
        integer (ip) :: imm1
        integer (ip) :: isym
        integer (ip) :: j
        integer (ip) :: jdg
        integer (ip) :: jdgs
        integer (ip) :: k
        integer (ip) :: ls
        integer (ip) :: m
        integer (ip) :: mdab
        integer (ip) :: mdo
        integer (ip) :: mmax
        integer (ip) :: modl
        integer (ip) :: mp1
        integer (ip) :: mp2
        integer (ip) :: ndab
        integer (ip) :: ndo
        integer (ip) :: nlat
        integer (ip) :: nlon
        integer (ip) :: nlp1
        integer (ip) :: np1
        integer (ip) :: nt
        real (wp) :: tsn
        real (wp) :: whrfft
        real (wp) :: work
        real (wp) :: wzfin
        real (wp) :: zb
        !
        !     whrfft must have at least nlon+15 locations
        !     wzfin must have 2*l*(nlat+1)/2 + ((l-3)*l+2)/2 locations
        !     zb must have 3*l*(nlat+1)/2 locations
        !     work must have ls*nlon locations
        !
        dimension g(idgs, jdgs, *), a(mdab, ndab, *), b(mdab, ndab, *), &
            ge(idg, jdg, *), go(idg, jdg, *), zb(imid, nlat, 3), wzfin(*), &
            whrfft(*), work(*)

        type (HFFTpack)      :: hfft
        type (SpherepackAux) :: sphere_aux

        ls = idg
        nlon = jdg
        mmax = min(nlat, nlon/2+1)
        mdo = mmax
        if (mdo+mdo-1 > nlon) mdo = mmax-1
        nlp1 = nlat+1
        tsn = 2./nlon
        fsn = 4./nlon
        modl = mod(nlat, 2)
        imm1 = imid
        if (modl /= 0) imm1 = imid-1
        if (isym /= 0) then
            goto 15
        end if

        do k=1, nt
            do i=1, imm1
                do j=1, nlon
                    ge(i, j, k) = tsn*(g(i, j, k)+g(nlp1-i, j, k))
                    go(i, j, k) = tsn*(g(i, j, k)-g(nlp1-i, j, k))
                end do
            end do
        end do
        goto 30

        15 do k=1, nt
            do i=1, imm1
                do j=1, nlon
                    ge(i, j, k) = fsn*g(i, j, k)
                end do
            end do
        end do

        if (isym == 1) then
            goto 27
        end if

30      if (modl == 0) then
            goto 27
        end if

        do k=1, nt
            do j=1, nlon
                ge(imid, j, k) = tsn*g(imid, j, k)
            end do
        end do

        27 do k=1, nt
            call hfft%forward(ls, nlon, ge(1, 1, k), ls, whrfft, work)
            if (mod(nlon, 2) /= 0) exit !goto 35
            do i=1, ls
                ge(i, nlon, k) = .5*ge(i, nlon, k)
            end do
        end do

        do k=1, nt
            do mp1=1, mmax
                do np1=mp1, nlat
                    a(mp1, np1, k) = 0.
                    b(mp1, np1, k) = 0.
                end do
            end do
        end do
        if (isym == 1) then
            goto 145
        end if

        call sphere_aux%zfin(2, nlat, nlon, 0, zb, i3, wzfin)

        do k=1, nt
            do i=1, imid
                do np1=1, nlat, 2
                    a(1, np1, k) = a(1, np1, k)+zb(i, np1, i3)*ge(i, 1, k)
                end do
            end do
        end do

        ndo = nlat

        if (mod(nlat, 2) == 0) then
            ndo = nlat-1
        end if

        do mp1=2, mdo
            m = mp1-1
            call sphere_aux%zfin(2, nlat, nlon, m, zb, i3, wzfin)
            do k=1, nt
                do i=1, imid
                    do np1=mp1, ndo, 2
                        a(mp1, np1, k) = a(mp1, np1, k)+zb(i, np1, i3)*ge(i, 2*mp1-2, k)
                        b(mp1, np1, k) = b(mp1, np1, k)+zb(i, np1, i3)*ge(i, 2*mp1-1, k)
                    end do
                end do
            end do
        end do

        if (mdo == mmax .or. mmax > ndo) then
            goto 135
        end if

        call sphere_aux%zfin(2, nlat, nlon, mdo, zb, i3, wzfin)

        do k=1, nt
            do i=1, imid
                do np1=mmax, ndo, 2
                    a(mmax, np1, k) = a(mmax, np1, k)+zb(i, np1, i3)*ge(i, 2*mmax-2, k)
                end do
            end do
        end do

135     if (isym == 2) then
            return
        end if

145     call sphere_aux%zfin(1, nlat, nlon, 0, zb, i3, wzfin)

        do k=1, nt
            do i=1, imm1
                do np1=2, nlat, 2
                    a(1, np1, k) = a(1, np1, k)+zb(i, np1, i3)*go(i, 1, k)
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
            call sphere_aux%zfin(1, nlat, nlon, m, zb, i3, wzfin)
            do k=1, nt
                do i=1, imm1
                    do np1=mp2, ndo, 2
                        a(mp1, np1, k) = a(mp1, np1, k)+zb(i, np1, i3)*go(i, 2*mp1-2, k)
                        b(mp1, np1, k) = b(mp1, np1, k)+zb(i, np1, i3)*go(i, 2*mp1-1, k)
                    end do
                end do
            end do
        end do

        mp2 = mmax+1

        if (mdo == mmax .or. mp2 > ndo) then
            return
        end if

        call sphere_aux%zfin(1, nlat, nlon, mdo, zb, i3, wzfin)

        do k=1, nt
            do i=1, imm1
                do np1=mp2, ndo, 2
                    a(mmax, np1, k) = a(mmax, np1, k)+zb(i, np1, i3)*go(i, 2*mmax-2, k)
                end do
            end do
        end do

    end subroutine shaec1




    subroutine shaeci(nlat, nlon, wshaec, lshaec, dwork, ldwork, ierror)
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        real (wp),    intent (out) :: wshaec(lshaec)
        integer (ip), intent (in)  :: lshaec
        real (wp),    intent (out) :: dwork(ldwork)
        integer (ip), intent (in)  :: ldwork
        integer (ip), intent (out) :: ierror
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip)         :: imid, iw1, labc, lzz1, mmax
        type (HFFTpack)      :: hfft
        type (SpherepackAux) :: sphere_aux
        !--------------------------------------------------------------

        imid = (nlat+1)/2
        mmax = min(nlat, nlon/2+1)
        lzz1 = 2*nlat*imid
        labc = 3*((mmax-2)*(2*nlat-mmax-1))/2

        ! Check validity of input arguments
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 4) then
            ierror = 2
            return
        else if (lshaec < lzz1+labc+nlon+15) then
            ierror = 3
            return
        else if (ldwork < nlat+1) then
            ierror = 4
            return
        else
            ierror = 0
        end if

        call sphere_aux%zfinit(nlat, nlon, wshaec, dwork)

        ! Set workspace pointer
        iw1 = lzz1+labc+1

        call hfft%initialize(nlon, wshaec(iw1))

    end subroutine shaeci

end module module_shaec
