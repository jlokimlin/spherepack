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
! ... file shagc.f
!
!     this file contains code and documentation for subroutines
!     shagc and shagci
!
! ... files which must be loaded with shagc.f
!
!     type_SpherepackAux.f, type_HFFTpack.f, gaqd.f
!
!
!     subroutine shagc(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, 
!    +                 wshagc, lshagc, work, lwork, ierror)
!
!     subroutine shagc performs the spherical harmonic analysis
!     on the array g and stores the result in the arrays a and b.
!     the analysis is performed on a gaussian grid in colatitude
!     and an equally spaced grid in longitude.  the associated
!     legendre functions are recomputed rather than stored as they
!     are in subroutine shags.  the analysis is described below
!     at output parameters a, b.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are compu
!            in radians in theta(1), ..., theta(nlat) by subroutine gaqd.
!            if nlat is odd the equator will be included as the grid poi
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
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
!     nt     the number of analyses.  in the program that calls shagc, 
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
!            g(i, j) contains the value of the function at the gaussian
!            point theta(i) and longitude point phi(j) = (j-1)*2*pi/nlon
!            the index ranges are defined above at the input parameter
!            isym.
!
!     idg    the first dimension of the array g as it appears in the
!            program that calls shagc. if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg must
!            be at least nlat/2 if nlat is even or at least (nlat+1)/2
!            if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shagc. jdg must be at least nlon.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shagc. mdab must be at least
!            min((nlon+2)/2, nlat) if nlon is even or at least
!            min((nlon+1)/2, nlat) if nlon is odd
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shaec. ndab must be at least nlat
!
!     wshagc an array which must be initialized by subroutine shagci.
!            once initialized, wshagc can be used repeatedly by shagc.
!            as long as nlat and nlon remain unchanged.  wshagc must
!            not be altered between calls of shagc.
!
!     lshagc the dimension of the array wshagc as it appears in the
!            program that calls shagc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshagc must be at least
!
!                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shagc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if isym is zero then lwork must be at least
!
!                      nlat*(nlon*nt+max(3*l2, nlon))
!
!            if isym is not zero then lwork must be at least
!
!                      l2*(nlon*nt+max(3*nlat, nlon))
!
!     **************************************************************
!
!     output parameters
!
!     a, b    both a, b are two or three dimensional arrays (see input
!            parameter nt) that contain the spherical harmonic
!            coefficients in the representation of g(i, j) given in the
!            discription of subroutine shagc. for isym=0, a(m, n) and
!            b(m, n) are given by the equations listed below. symmetric
!            versions are used when isym is greater than zero.
!
!     definitions
!
!     1. the normalized associated legendre functions
!
!     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                       *sin(theta)**m/(2**n*factorial(n)) times the
!                       (n+m)th derivative of (x**2-1)**n with respect
!                       to x=cos(theta).
!
!     2. the fourier transform of g(i, j).
!
!     c(m, i)          = 2/nlon times the sum from j=1 to j=nlon of
!                       g(i, j)*cos((m-1)*(j-1)*2*pi/nlon)
!                       (the first and last terms in this sum
!                       are divided by 2)
!
!     s(m, i)          = 2/nlon times the sum from j=2 to j=nlon of
!                       g(i, j)*sin((m-1)*(j-1)*2*pi/nlon)
!
!
!     3. the gaussian points and weights on the sphere
!        (computed by subroutine gaqd).
!
!        theta(1), ..., theta(nlat) (gaussian pts in radians)
!        wts(1), ..., wts(nlat) (corresponding gaussian weights)
!
!     4. the maximum (plus one) longitudinal wave number
!
!            mmax = min(nlat, (nlon+2)/2) if nlon is even or
!            mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!
!     then for m=0, ..., mmax-1 and n=m, ..., nlat-1 the arrays a, b
!     are given by
!
!     a(m+1, n+1)     =  the sum from i=1 to i=nlat of
!                       c(m+1, i)*wts(i)*pbar(m, n, theta(i))
!
!     b(m+1, n+1)      = the sum from i=1 to nlat of
!                       s(m+1, i)*wts(i)*pbar(m, n, theta(i))
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
!            = 9  error in the specification of lshagc
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!
!     subroutine shagci(nlat, nlon, wshagc, lshagc, dwork, ldwork, ierror)
!
!     subroutine shagci initializes the array wshagc which can then
!     be used repeatedly by subroutines shagc. it precomputes
!     and stores in wshagc quantities such as gaussian weights, 
!     legendre polynomial coefficients, and fft trigonometric tables.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are compu
!            in radians in theta(1), ..., theta(nlat) by subroutine gaqd.
!            if nlat is odd the equator will be included as the grid poi
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     wshagc an array which must be initialized by subroutine shagci.
!            once initialized, wshagc can be used repeatedly by shagc
!            as long as nlat and nlon remain unchanged.  wshagc must
!            not be altered between calls of shagc.
!
!     lshagc the dimension of the array wshagc as it appears in the
!            program that calls shagc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshagc must be at least
!
!                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
!
!     dwork   a real work array that does not have to be saved.
!
!     ldwork  the dimension of the array dwork as it appears in the
!            program that calls shagci. ldwork must be at least
!
!                nlat*(nlat+4)
!
!     output parameter
!
!     wshagc an array which must be initialized before calling shagc or
!            once initialized, wshagc can be used repeatedly by shagc or
!            as long as nlat and nlon remain unchanged.  wshagc must not
!            altered between calls of shagc.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshagc
!            = 4  error in the specification of ldwork
!            = 5  failure in gaqd to compute gaussian points
!                 (due to failure in eigenvalue routine)
!
!
module module_shagc

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    use module_gaqd, only: &
        gaqd

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: shagc
    public :: shagci

contains

    subroutine shagc(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
        wshagc, lshagc, work, lwork, ierror)

        real :: a
        real :: b
        real :: g
        integer :: idg
        integer :: ierror
        integer :: ifft
        integer :: ipmn
        integer :: isym
        integer :: iwts
        integer :: jdg
        integer :: l
        integer :: l1
        integer :: l2
        integer :: lat
        integer :: late
        integer :: lshagc
        integer :: lwork
        integer :: mdab
        integer :: ndab
        integer :: nlat
        integer :: nlon
        integer :: nt
        real :: work
        real :: wshagc
        !     subroutine shagc performs the spherical harmonic analysis on
        !     a gaussian grid on the array(s) in g and returns the coefficients
        !     in array(s) a, b. the necessary legendre polynomials are computed
        !     as needed in this version.
        !
        dimension g(idg, jdg, *), a(mdab, ndab, *), b(mdab, ndab,*), &
            wshagc(lshagc), work(lwork)
        !     check input parameters
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        ierror = 3
        if (isym < 0 .or.isym > 2) return
        ierror = 4
        if (nt < 1) return
        !     set upper limit on m for spherical harmonic basis
        l = min((nlon+2)/2, nlat)
        !     set gaussian point nearest equator pointer
        late = (nlat+mod(nlat, 2))/2
        !     set number of grid points for analysis/synthesis
        lat = nlat
        if (isym /= 0) lat = late
        ierror = 5
        if (idg < lat) return
        ierror = 6
        if (jdg < nlon) return
        ierror = 7
        if (mdab < l) return
        ierror = 8
        if (ndab < nlat) return
        l1 = l
        l2 = late
        ierror = 9
        !     check permanent work space length
        if (lshagc < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
        ierror = 10
        !     check temporary work space length
        if (isym == 0) then
            if (lwork <nlat*(nlon*nt+max(3*l2, nlon))) return
        else
            !     isym.ne.0
            if (lwork <l2*(nlon*nt+max(3*nlat, nlon))) return
        end if
        ierror = 0
        !     starting address for gaussian wts in shigc and fft values
        iwts = 1
        ifft = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+1
        !     set pointers for internal storage of g and legendre polys
        ipmn = lat*nlon*nt+1
        call shagc1(nlat, nlon, l, lat, isym, g, idg, jdg, nt, a, b, mdab, ndab, &
            wshagc, wshagc(iwts), wshagc(ifft), late, work(ipmn), work)

    contains


        subroutine shagc1(nlat, nlon, l, lat, mode, gs, idg, jdg, nt, a, b, mdab, &
            ndab, w, wts, wfft, late, pmn, g)

            real :: a
            real :: b
            real :: g
            real :: gs
            integer :: i
            integer :: idg
            integer :: is
            integer :: j
            integer :: jdg
            integer :: k
            integer :: km
            integer :: l
            integer :: lat
            integer :: late
            integer :: lm1
            integer :: lp1
            integer :: m
            integer :: mdab
            integer :: mode
            integer :: mp1
            integer :: mp2
            integer :: ms
            integer :: ndab
            integer :: nl2
            integer :: nlat
            integer :: nlon
            integer :: np1
            integer :: ns
            integer :: nt
            real :: pmn
            real :: sfn
            real :: t1
            real :: t2
            real :: w
            real :: wfft
            real :: wts
            dimension gs(idg, jdg, nt), a(mdab, ndab, nt), &
                b(mdab, ndab, nt), g(lat, nlon, nt)
            dimension w(*), wts(nlat), wfft(*), pmn(nlat, late, 3)

            type (HFFTpack)      :: hfft
            type (SpherepackAux) :: sphere_aux

            !     set gs array internally in shagc1
            do k=1, nt
                do j=1, nlon
                    do i=1, lat
                        g(i, j, k) = gs(i, j, k)
                    end do
                end do
            end do
            !     do fourier transform
            do k=1, nt
                call hfft%forward(lat, nlon, g(1, 1, k), lat, wfft, pmn)
            end do
            !     scale result
            sfn = 2.0/real(nlon)
            do k=1, nt
                do j=1, nlon
                    do i=1, lat
                        g(i, j, k) = sfn*g(i, j, k)
                    end do
                end do
            end do
            !     compute using gaussian quadrature
            !     a(n, m) = s (ga(theta, m)*pnm(theta)*sin(theta)*dtheta)
            !     b(n, m) = s (gb(theta, m)*pnm(theta)*sin(theta)*dtheta)
            !     here ga, gb are the cos(phi), sin(phi) coefficients of
            !     the fourier expansion of g(theta, phi) in phi.  as a result
            !     of the above fourier transform they are stored in array
            !     g as follows:
            !     for each theta(i) and k= l-1
            !     ga(0), ga(1), gb(1), ga(2), gb(2), ..., ga(k-1), gb(k-1), ga(k)
            !     correspond to (in the case nlon=l+l-2)
            !     g(i, 1), g(i, 2), g(i, 3), g(i, 4), g(i, 5), ..., g(i, 2l-4), g(i, 2l-3), g(i, 2l-
            !     initialize coefficients to zero
            do k=1, nt
                do np1=1, nlat
                    do mp1=1, l
                        a(mp1, np1, k) = 0.0
                        b(mp1, np1, k) = 0.0
                    end do
                end do
            end do
            !     set m+1 limit on b(m+1) calculation
            lm1 = l
            if (nlon == 2*l-2) lm1 = l-1
            if (mode == 0) then
                !     for full sphere (mode=0) and even/odd reduction:
                !     overwrite g(i) with (g(i)+g(nlat-i+1))*wts(i)
                !     overwrite g(nlat-i+1) with (g(i)-g(nlat-i+1))*wts(i)
                nl2 = nlat/2
                do k=1, nt
                    do j=1, nlon
                        do i=1, nl2
                            is = nlat-i+1
                            t1 = g(i, j, k)
                            t2 = g(is, j, k)
                            g(i, j, k) = wts(i)*(t1+t2)
                            g(is, j, k) = wts(i)*(t1-t2)
                        end do
                        !     adjust equator if necessary(nlat odd)
                        if (mod(nlat, 2)/=0) g(late, j, k) = wts(late)*g(late, j, k)
                    end do
                end do
                !     set m = 0 coefficients first
                m = 0
                call sphere_aux%legin(mode, l, nlat, m, w, pmn, km)
                do k=1, nt
                    do i=1, late
                        is = nlat-i+1
                        do np1=1, nlat, 2
                            !     n even
                            a(1, np1, k) = a(1, np1, k)+g(i, 1, k)*pmn(np1, i, km)
                        end do
                        do np1=2, nlat, 2
                            !     n odd
                            a(1, np1, k) = a(1, np1, k)+g(is, 1, k)*pmn(np1, i, km)
                        end do
                    end do
                end do
                !     compute coefficients for which b(m, n) is available
                do mp1=2, lm1
                    m = mp1-1
                    mp2 = m+2
                    !     compute pmn for all i and n=m, ..., l-1
                    call sphere_aux%legin(mode, l, nlat, m, w, pmn, km)
                    do k=1, nt
                        do i=1, late
                            is = nlat-i+1
                            !     n-m even
                            do np1=mp1, nlat, 2
                                a(mp1, np1, k) = a(mp1, np1, k)+g(i, 2*m, k)*pmn(np1, i, km)
                                b(mp1, np1, k) = b(mp1, np1, k)+g(i, 2*m+1, k)*pmn(np1, i, km)
                            end do
                            !     n-m odd
                            do np1=mp2, nlat, 2
                                a(mp1, np1, k) = a(mp1, np1, k)+g(is, 2*m, k)*pmn(np1, i, km)
                                b(mp1, np1, k) = b(mp1, np1, k)+g(is, 2*m+1, k)*pmn(np1, i, km)
                            end do
                        end do
                    end do
                end do
                if (nlon == 2*l-2) then
                    !     compute a(l, np1) coefficients only
                    m = l-1
                    call sphere_aux%legin(mode, l, nlat, m, w, pmn, km)
                    do k=1, nt
                        do i=1, late
                            is = nlat-i+1
                            !     n-m even
                            do np1=l, nlat, 2
                                a(l, np1, k) = a(l, np1, k)+0.5*g(i, nlon, k)*pmn(np1, i, km)
                            end do
                            lp1 = l+1
                            !     n-m odd
                            do np1=lp1, nlat, 2
                                a(l, np1, k) = a(l, np1, k)+0.5*g(is, nlon, k)*pmn(np1, i, km)
                            end do
                        end do
                    end do
                end if
            else
                !     half sphere
                !     overwrite g(i) with wts(i)*(g(i)+g(i)) for i=1, ..., nlate/2
                nl2 = nlat/2
                do  k=1, nt
                    do j=1, nlon
                        do i=1, nl2
                            g(i, j, k) = wts(i)*(g(i, j, k)+g(i, j, k))
                        end do
                        !     adjust equator separately if a grid point
                        if (nl2<late) g(late, j, k) = wts(late)*g(late, j, k)
                    end do
                end do
                !     set m = 0 coefficients first
                m = 0
                call sphere_aux%legin(mode, l, nlat, m, w, pmn, km)
                ms = 1
                if (mode == 1) ms = 2
                do k=1, nt
                    do i=1, late
                        do np1=ms, nlat, 2
                            a(1, np1, k) = a(1, np1, k)+g(i, 1, k)*pmn(np1, i, km)
                        end do
                    end do
                end do
                !     compute coefficients for which b(m, n) is available
                do mp1=2, lm1
                    m = mp1-1
                    ms = mp1
                    if (mode == 1) ms = mp1+1
                    !     compute pmn for all i and n=m, ..., nlat-1
                    call sphere_aux%legin(mode, l, nlat, m, w, pmn, km)
                    do  k=1, nt
                        do  i=1, late
                            do np1=ms, nlat, 2
                                a(mp1, np1, k) = a(mp1, np1, k)+g(i, 2*m, k)*pmn(np1, i, km)
                                b(mp1, np1, k) = b(mp1, np1, k)+g(i, 2*m+1, k)*pmn(np1, i, km)
                            end do
                        end do
                    end do
                end do
                if (nlon==2*l-2) then
                    !     compute coefficient a(l, np1) only
                    m = l-1
                    call sphere_aux%legin(mode, l, nlat, m, w, pmn, km)
                    ns = l
                    if (mode == 1) ns = l+1
                    do k=1, nt
                        do i=1, late
                            do np1=ns, nlat, 2
                                a(l, np1, k) = a(l, np1, k)+0.5*g(i, nlon, k)*pmn(np1, i, km)
                            end do
                        end do
                    end do
                end if
            end if

        end subroutine shagc1


    end subroutine shagc

    subroutine shagci(nlat, nlon, wshagc, lshagc, dwork, ldwork, ierror)
        implicit none
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: i6
        integer :: i7
        integer :: idth
        integer :: idwts
        integer :: ierror
        integer :: iw
        integer :: l
        integer :: l1
        integer :: l2
        integer :: late
        integer :: ldwork
        integer :: lshagc
        integer :: nlat
        integer :: nlon
        real :: wshagc
        !     this subroutine must be called before calling shagc with
        !     fixed nlat, nlon. it precomputes quantites such as the gaussian
        !     points and weights, m=0, m=1 legendre polynomials, recursion
        !     recursion coefficients.
        dimension wshagc(lshagc)
        real dwork(ldwork)
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        !     set triangular truncation limit for spherical harmonic basis
        l = min((nlon+2)/2, nlat)
        !     set equator or nearest point (if excluded) pointer
        late = (nlat+mod(nlat, 2))/2
        l1 = l
        l2 = late
        ierror = 3
        !     check permanent work space length
        if (lshagc < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
        ierror = 4
        if (ldwork<nlat*(nlat+4))return
        ierror = 0
        !     set pointers
        i1 = 1
        i2 = i1+nlat
        i3 = i2+nlat*late
        i4 = i3+nlat*late
        i5 = i4+l*(l-1)/2 +(nlat-l)*(l-1)
        i6 = i5+l*(l-1)/2 +(nlat-l)*(l-1)
        i7 = i6+l*(l-1)/2 +(nlat-l)*(l-1)
        !     set indices in temp work for real gaussian wts and pts
        idth = 1
        idwts = idth+nlat
        iw = idwts+nlat
        call shagci1(nlat, nlon, l, late, wshagc(i1), wshagc(i2), wshagc(i3), &
            wshagc(i4), wshagc(i5), wshagc(i6), wshagc(i7), dwork(idth), &
            dwork(idwts), dwork(iw), ierror)
        if (ierror /= 0) ierror = 5

    contains


        subroutine shagci1(nlat, nlon, l, late, wts, p0n, p1n, abel, bbel, cbel, &
            wfft, dtheta, dwts, work, ier)
            real :: abel
            real :: bbel
            real :: cbel
            integer :: i
            integer :: ier
            integer :: imn
            integer :: imndx
            integer :: indx
            integer :: l
            integer :: late
            integer :: lw
            integer :: m
            integer :: mlim
            integer :: n
            integer :: nlat
            integer :: nlon
            integer :: np1
            real :: p0n
            real :: p1n
            real :: wfft
            real :: wts
            dimension wts(nlat), p0n(nlat, late), p1n(nlat, late), abel(*), bbel(*), &
                cbel(*), wfft(*)
            real pb, dtheta(nlat), dwts(nlat), work(*)

            type (HFFTpack)      :: hfft
            type (SpherepackAux) :: sphere_aux

            !     compute the nlat  gaussian points and weights, the
            !     m=0, 1 legendre polys for gaussian points and all n,
            !     and the legendre recursion coefficients
            !     define index function used in storing
            !     arrays for recursion coefficients (functions of (m, n))
            !     the index function indx(m, n) is defined so that
            !     the pairs (m, n) map to [1, 2, ..., indx(l-1, l-1)] with no
            !     "holes" as m varies from 2 to n and n varies from 2 to l-1.
            !     (m=0, 1 are set from p0n, p1n for all n)
            !     define for 2.le.n.le.l-1
            !indx(m, n) = imn = (n-1)*(n-2)/2+m-1
            !     define index function for l.le.n.le.nlat
            !imndx(m, n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1
            !     preset quantites for fourier transform
            call hfft%initialize(nlon, wfft)
            !     compute real gaussian points and weights
            !     lw = 4*nlat*(nlat+1)+2
            lw = nlat*(nlat+2)
            call gaqd(nlat, dtheta, dwts, work(1), lw, ier)
            if (ier/=0) return
            !     store gaussian weights single precision to save computation
            !     in inner loops in analysis
            do i=1, nlat
                wts(i) = dwts(i)
            end do
            !     initialize p0n, p1n using real dnlfk, dnlft
            do np1=1, nlat
                do i=1, late
                    p0n(np1, i) = 0.0
                    p1n(np1, i) = 0.0
                end do
            end do
            !     compute m=n=0 legendre polynomials for all theta(i)
            np1 = 1
            n = 0
            m = 0
            call sphere_aux%dnlfk(m, n, work)
            do i=1, late
                call sphere_aux%dnlft(m, n, dtheta(i), work, pb)
                p0n(1, i) = pb
            end do
            !     compute p0n, p1n for all theta(i) when n.gt.0
            do np1=2, nlat
                n = np1-1
                m = 0
                call sphere_aux%dnlfk(m, n, work)
                do i=1, late
                    call sphere_aux%dnlft(m, n, dtheta(i), work, pb)
                    p0n(np1, i) = pb
                end do
                !     compute m=1 legendre polynomials for all n and theta(i)
                m = 1
                call sphere_aux%dnlfk(m, n, work)
                do i=1, late
                    call sphere_aux%dnlft(m, n, dtheta(i), work, pb)
                    p1n(np1, i) = pb
                end do
            end do
            !     compute and store swarztrauber recursion coefficients
            !     for 2.le.m.le.n and 2.le.n.le.nlat in abel, bbel, cbel
            do n=2, nlat
                mlim = min(n, l)
                do m=2, mlim
                    if (n >= l) then
                        imn = l*(l-1)/2+(n-l-1)*(l-1)+m-1
                    else
                        imn = (n-1)*(n-2)/2+m-1
                    end if
                    abel(imn)=sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                        real(((2*n-3)*(m+n-1)*(m+n))))
                    bbel(imn)=sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                        real(((2*n-3)*(m+n-1)*(m+n))))
                    cbel(imn)=sqrt(real((n-m+1)*(n-m+2))/ &
                        real(((n+m-1)*(n+m))))
                end do
            end do

        end subroutine shagci1

    end subroutine shagci

end module module_shagc
