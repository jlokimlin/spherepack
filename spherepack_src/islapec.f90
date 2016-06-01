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
! ... file islapec.f
!
!     this file includes documentation and code for
!     subroutine islapec         i
!
! ... files which must be loaded with islapec.f
!
!     sphcom.f, hrfft.f, shaec.f, shsec.f
!
!     subroutine islapec(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, 
!    +mdab, ndab, wshsec, lshsec, work, lwork, pertrb, ierror)
!
!     islapec inverts the laplace or helmholz operator on an equally
!     spaced latitudinal grid using o(n**2) storage. given the
!     spherical harmonic coefficients a(m, n) and b(m, n) of the right
!     hand side slap(i, j), islapec computes a solution sf(i, j) to
!     the following helmhotz equation :
!
!           2                2
!     [d(sf(i, j))/dlambda /sint + d(sint*d(sf(i, j))/dtheta)/dtheta]/sint
!
!                   - xlmbda * sf(i, j) = slap(i, j)
!
!      where sf(i, j) is computed at colatitude
!
!                 theta(i) = (i-1)*pi/(nlat-1)
!
!            and longitude
!
!                 lambda(j) = (j-1)*2*pi/nlon
!
!            for i=1, ..., nlat and j=1, ..., nlon.
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
!     nlon   the number of distinct longitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     isym   this parameter should have the same value input to subroutine
!            shaec to compute the coefficients a and b for the scalar field
!            slap.  isym is set as follows:
!
!            = 0  no symmetries exist in slap about the equator. scalar
!                 synthesis is used to compute sf on the entire sphere.
!                 i.e., in the array sf(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon.
!
!           = 1  sf and slap are antisymmetric about the equator. the
!                synthesis used to compute sf is performed on the
!                northern hemisphere only.  if nlat is odd, sf(i, j) is
!                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
!                nlat is even, sf(i, j) is computed for i=1, ..., nlat/2
!                and j=1, ..., nlon.
!
!
!           = 2  sf and slap are symmetric about the equator. the
!                synthesis used to compute sf is performed on the
!                northern hemisphere only.  if nlat is odd, sf(i, j) is
!                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
!                nlat is even, sf(i, j) is computed for i=1, ..., nlat/2
!                and j=1, ..., nlon.
!
!
!   nt       the number of solutions. in the program that calls islapec
!            the arrays sf, a, and b can be three dimensional in which
!            case multiple solutions are computed. the third index
!            is the solution index with values k=1, ..., nt.
!            for a single solution set nt=1. the description of the
!            remaining parameters is simplified by assuming that nt=1
!            and sf, a, b are two dimensional.
!
!   xlmbda   a one dimensional array with nt elements. if xlmbda is
!            is identically zero islapec solves poisson's equation.
!            if xlmbda > 0.0 islapec solves the helmholtz equation.
!            if xlmbda < 0.0 the nonfatal error flag ierror=-1 is
!            returned. negative xlambda could result in a division
!            by zero.
!
!   ids      the first dimension of the array sf as it appears in the
!            program that calls islapec.  if isym = 0 then ids must be at
!            least nlat.  if isym > 0 and nlat is even then ids must be
!            at least nlat/2. if isym > 0 and nlat is odd then ids must
!            be at least (nlat+1)/2.
!
!   jds      the second dimension of the array sf as it appears in the
!            program that calls islapec. jds must be at least nlon.
!
!
!   a, b      two or three dimensional arrays (see input parameter nt)
!            that contain scalar spherical harmonic coefficients
!            of the scalar field slap. a, b must be computed by shaec
!            prior to calling islapec.
!
!
!   mdab     the first dimension of the arrays a and b as it appears
!            in the program that calls islapec.  mdab must be at
!            least min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!   ndab     the second dimension of the arrays a and b as it appears
!            in the program that calls islapec. ndab must be at least
!            least nlat.
!
!            mdab, ndab should have the same values input to shaec to
!            compute the coefficients a and b.
!
!
!   wshsec   an array which must be initialized by subroutine shseci.
!            once initialized, wshsec can be used repeatedly by
!            islapec as long as nlat and nlon  remain unchanged.
!            wshsec must not be altered between calls of islapec.
!
!   lshsec   the dimension of the array wshsec as it appears in the
!            program that calls islapec.  let
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lsave must be greater than or equal to
!
!               2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls islapec. define
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            if isym = 0 let
!
!               lwkmin = nlat*(2*nt*nlon+max(6*l2, nlon)+2*nt*l1+1).
!
!            if isym > 0 let
!
!               lwkmin = l2*(2*nt*nlon+max(6*nlat, nlon))+nlat*(2*nt*l1+1)
!
!
!     then lwork must be greater than or equal to lwkmin (see ierror=10)
!
!     **************************************************************
!
!     output parameters
!
!
!    sf      two or three dimensional arrays (see input parameter nt)
!            that contain the solution to either the helmholtz 
!            (xlmbda>0.0) or poisson's equation. sf(i, j) is computed
!            at colatitude
!
!                 theta(i) = (i-1)*pi/(nlat-1)
!
!            and longitude
!
!                 lambda(j) = (j-1)*2*pi/nlon
!
!            for i=1, ..., nlat and j=1, ..., nlon.
!
!   pertrb  a one dimensional array with nt elements (see input 
!           parameter nt). in the discription that follows we assume 
!           that nt=1. if xlmbda > 0.0 then pertrb=0.0 is always 
!           returned because the helmholtz operator is invertible.
!           if xlmbda = 0.0 then a solution exists only if a(1, 1)
!           is zero. islapec sets a(1, 1) to zero. the resulting
!           solution sf(i, j) solves poisson's equation with
!           pertrb = a(1, 1)/(2.*sqrt(2.)) subtracted from the
!           right side slap(i, j).
!
!
!  ierror   a parameter which flags errors in input parameters as follows:
!
!            =-1  xlmbda is input negative (nonfatal error)
!
!            = 0  no errors detected
!
!            = 1  error in the specification of nlat
!
!            = 2  error in the specification of nlon
!
!            = 3  error in the specification of ityp
!
!            = 4  error in the specification of nt
!
!            = 5  error in the specification of ids
!
!            = 6  error in the specification of jds
!
!            = 7  error in the specification of mdbc
!
!            = 8  error in the specification of ndbc
!
!            = 9  error in the specification of lsave
!
!            = 10 error in the specification of lwork
!
!
! **********************************************************************
!                                                                              
!     end of documentation for islapec
!
! **********************************************************************
!


subroutine islapec(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
    mdab, ndab, wshsec, lshsec, work, lwork, pertrb, ierror)
    implicit none
    real :: a
    real :: b
    integer :: ia
    integer :: ib
    integer :: ids
    integer :: ierror
    integer :: ifn
    integer :: imid
    integer :: isym
    integer :: iwk
    integer :: jds
    integer :: k
    integer :: l1
    integer :: l2
    integer :: ls
    integer :: lshsec
    integer :: lwk
    integer :: lwkmin
    integer :: lwmin
    integer :: lwork
    integer :: mdab
    integer :: mmax
    integer :: mn
    integer :: ndab
    integer :: nlat
    integer :: nln
    integer :: nlon
    integer :: nt
    real :: pertrb
    real :: sf
    real :: work
    real :: wshsec
    real :: xlmbda
    dimension sf(ids, jds, nt), a(mdab, ndab, nt), b(mdab, ndab, nt)
    dimension wshsec(lshsec), work(lwork), pertrb(nt), xlmbda(nt)
    !
    !     check input parameters
    !
    ierror = 1
    if (nlat < 3) return
    ierror = 2
    if (nlon < 4) return
    ierror = 3
    if (isym < 0 .or. isym > 2) return
    ierror = 4
    if (nt < 0) return
    ierror = 5
    imid = (nlat+1)/2
    if ((isym == 0 .and. ids<nlat) .or. &
        (isym>0 .and. ids<imid)) return
    ierror = 6
    if (jds < nlon) return
    ierror = 7
    mmax = min(nlat, nlon/2+1)
    if (mdab < mmax) return
    ierror = 8
    if (ndab < nlat) return
    ierror = 9
    !
    !     set and verify saved work space length
    !
    !
    l1 = min(nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwmin = 2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
    if (lshsec < lwmin) return
    ierror = 10
    !
    !     set and verify unsaved work space length
    !
    ls = nlat
    if (isym > 0) ls = imid
    nln = nt*ls*nlon
    mn = mmax*nlat*nt
    !     lwmin = nln+ls*nlon+2*mn+nlat
    !     if (lwork .lt. lwmin) return
    l2 = (nlat+1)/2
    l1 = min(nlat, nlon/2+1)
    if (isym == 0) then
        lwkmin = nlat*(2*nt*nlon+max(6*l2, nlon)+2*nt*l1+1)
    else
        lwkmin = l2*(2*nt*nlon+max(6*nlat, nlon))+nlat*(2*nt*l1+1)
    end if
    if (lwork < lwkmin) return
    ierror = 0
    !
    !     check sign of xlmbda
    !
    do  k=1, nt
        if (xlmbda(k)<0.0) then
            ierror = -1
        end if
    end do
    !
    !     set work space pointers
    !
    ia = 1
    ib = ia+mn
    ifn = ib+mn
    iwk = ifn+nlat
    lwk = lwork-2*mn-nlat
    call islpec1(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, mdab, ndab, &
        work(ia), work(ib), mmax, work(ifn), wshsec, lshsec, work(iwk), lwk, &
        pertrb, ierror)
    return
end subroutine islapec

subroutine islpec1(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
    mdab, ndab, as, bs, mmax, fnn, wshsec, lshsec, wk, lwk, pertrb, ierror)
    implicit none
    real :: a
    real :: as
    real :: b
    real :: bs
    real :: fn
    real :: fnn
    integer :: ids
    integer :: ierror
    integer :: isym
    integer :: jds
    integer :: k
    integer :: lshsec
    integer :: lwk
    integer :: m
    integer :: mdab
    integer :: mmax
    integer :: n
    integer :: ndab
    integer :: nlat
    integer :: nlon
    integer :: nt
    real :: pertrb
    real :: sf
    real :: wk
    real :: wshsec
    real :: xlmbda
    dimension sf(ids, jds, nt), a(mdab, ndab, nt), b(mdab, ndab, nt)
    dimension as(mmax, nlat, nt), bs(mmax, nlat, nt), fnn(nlat)
    dimension wshsec(lshsec), wk(lwk), pertrb(nt), xlmbda(nt)
    !
    !     set multipliers and preset synthesis coefficients to zero
    !
    do n=1, nlat
        fn = real(n - 1)
        fnn(n) = fn*(fn+1.0)
        do m=1, mmax
            do k=1, nt
                as(m, n, k) = 0.0
                bs(m, n, k) = 0.0
            end do
        end do
    end do

    do k=1, nt
            !
            !     compute synthesis coefficients for xlmbda zero or nonzero
            !
        if (xlmbda(k) == 0.0) then
            do n=2, nlat
                as(1, n, k) = -a(1, n, k)/fnn(n)
                bs(1, n, k) = -b(1, n, k)/fnn(n)
            end do
            do m=2, mmax
                do n=m, nlat
                    as(m, n, k) = -a(m, n, k)/fnn(n)
                    bs(m, n, k) = -b(m, n, k)/fnn(n)
                end do
            end do
        else
            !
            !     xlmbda nonzero so operator invertible unless
            !     -n*(n-1) = xlmbda(k) < 0.0  for some n
            !
            pertrb(k) = 0.0
            do n=1, nlat
                as(1, n, k) = -a(1, n, k)/(fnn(n)+xlmbda(k))
                bs(1, n, k) = -b(1, n, k)/(fnn(n)+xlmbda(k))
            end do
            do m=2, mmax
                do n=m, nlat
                    as(m, n, k) = -a(m, n, k)/(fnn(n)+xlmbda(k))
                    bs(m, n, k) = -b(m, n, k)/(fnn(n)+xlmbda(k))
                end do
            end do
        end if
    end do
    !
    !     synthesize as, bs into sf
    !
    call shsec(nlat, nlon, isym, nt, sf, ids, jds, as, bs, mmax, nlat, &
        wshsec, lshsec, wk, lwk, ierror)
    return
end subroutine islpec1
