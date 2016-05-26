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
! ... file shigc.f
!
!     this file contains code and documentation for subroutine shigc
!
! ... files which must be loaded with shigc.f
!
!     sphcom.f, hrfft.f, gaqd.f
!
!     3/6/98
!
! *** shigc is functionally the same as shagci or shsgci.  It
!     it included in this version of spherepack because 
!     older versions of spherepack call shigc to initialize
!     the saved work space wshigc, for either shagc or shsgc
!
!     subroutine shigc(nlat, nlon, wshigc, lshigc, dwork, ldwork, ierror)
!
!     subroutine shigc initializes the array wshigc which can then
!     be used repeatedly by subroutines shsgc or shagc. it precomputes
!     and stores in wshigc quantities such as gaussian weights, 
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
!     wshigc an array which must be initialized by subroutine shigc.
!            once initialized, wshigc can be used repeatedly by shsgc
!            or shagc as long as nlat and nlon remain unchanged.  wshigc
!            must not be altered between calls of shsgc or shagc.
!
!     lshigc the dimension of the array wshigc as it appears in the
!            program that calls shsgc or shagc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshigc must be at least
!
!                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
!
!     dwork  a real work array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls shigc. ldwork must be at least
!
!               nlat*(nlat+4)
!
!     output parameter
!
!     wshigc an array which must be initialized before calling shsgc or shagc.
!            once initialized, wshigc can be used repeatedly by shsgc or shagc
!            as long as nlat and nlon remain unchanged.  wshigc must not
!            altered between calls of shsgc or shagc
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshigc
!            = 4  error in the specification of ldwork
!            = 5  failure in gaqd to compute gaussian points
!                 (due to failure in eigenvalue routine)
!
!
! ****************************************************************
subroutine shigc(nlat, nlon, wshigc, lshigc, dwork, ldwork, ierror)
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
    integer :: lshigc
    integer :: nlat
    integer :: nlon
    real :: wshigc
    !     this subroutine must be called before calling shsgc/shagc with
    !     fixed nlat, nlon. it precomputes quantites such as the gaussian
    !     points and weights, m=0, m=1 legendre polynomials, recursion
    !     recursion coefficients.
    dimension wshigc(lshigc)
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
    if (lshigc < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
    ierror = 4
    !     if (lwork.lt.4*nlat*(nlat+2)+2) return
    if (ldwork < nlat*(nlat+4)) return
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
    !     idwts = idth+2*nlat
    !     iw = idwts+2*nlat
    idwts = idth+nlat
    iw = idwts+nlat
    call shigc1(nlat, nlon, l, late, wshigc(i1), wshigc(i2), wshigc(i3), &
        wshigc(i4), wshigc(i5), wshigc(i6), wshigc(i7), dwork(idth), &
        dwork(idwts), dwork(iw), ierror)
    if (ierror /= 0) ierror = 5

end subroutine shigc


subroutine shigc1(nlat, nlon, l, late, wts, p0n, p1n, abel, bbel, cbel, &
    wfft, dtheta, dwts, work, ier)
    implicit none
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
    dimension wts(nlat), p0n(nlat, late), p1n(nlat, late), abel(1), bbel(1), &
        cbel(1), wfft(1), dtheta(nlat), dwts(nlat)
    real pb, dtheta, dwts, work(*)
    !     compute the nlat  gaussian points and weights, the
    !     m=0, 1 legendre polys for gaussian points and all n,
    !     and the legendre recursion coefficients
    !     define index function used in storing
    !     arrays for recursion coefficients (functions of (m, n))
    !     the index function indx(m, n) is defined so that
    !     the pairs (m, n) map to [1, 2, ..., indx(l-1, l-1)] with no
    !     "holes" as m varies from 2 to n and n varies from 2 to l-1.
    !     (m=0, 1 are set from p0n, p1n for all n)
    !     define for 2<=n<=l-1
    indx(m, n) = (n-1)*(n-2)/2+m-1
    !     define index function for l<=n<=nlat
    imndx(m, n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1
    !     preset quantites for fourier transform
    call hrffti(nlon, wfft)
    !     compute real gaussian points and weights
    !     lw = 4*nlat*(nlat+1)+2
    lw = nlat*(nlat+2)

    call gaqd(nlat, dtheta, dwts, work, lw, ier)

    if (ier/=0) then
        return
    end if
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
    call dnlfk(m, n, work)

    do i=1, late
        call dnlft(m, n, dtheta(i), work, pb)
        p0n(1, i) = pb
    end do

    !     compute p0n, p1n for all theta(i) when n.gt.0
    do np1=2, nlat
        n = np1-1
        m = 0
        call dnlfk(m, n, work)
        do i=1, late
            call dnlft(m, n, dtheta(i), work, pb)
            p0n(np1, i) = pb
        end do
        !     compute m=1 legendre polynomials for all n and theta(i)
        m = 1
        call dnlfk(m, n, work)
        do i=1, late
            call dnlft(m, n, dtheta(i), work, pb)
            p1n(np1, i) = pb
        end do
    end do
    !     compute and store swarztrauber recursion coefficients
    !     for 2<=m<=n and 2<=n<=nlat in abel, bbel, cbel
    do n=2, nlat
        mlim = min(n, l)
        do m=2, mlim
            imn = indx(m, n)
            if (n >= l) then
                imn = imndx(m, n)
            end if
            abel(imn)=sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                real(((2*n-3)*(m+n-1)*(m+n))))
            bbel(imn)=sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                real(((2*n-3)*(m+n-1)*(m+n))))
            cbel(imn)=sqrt(real((n-m+1)*(n-m+2))/ &
                real(((n+m-1)*(n+m))))
        end do
    end do

end subroutine shigc1
