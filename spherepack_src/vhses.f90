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
! ... file vhses.f90
!
!     this file contains code and documentation for subroutines
!     vhses and vhsesi
!
! ... files which must be loaded with vhses.f90
!
!     sphcom.f90, hrfft.f90
!
!   
!     subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvhses, lvhses, work, lwork, ierror)
!
!     subroutine vhses performs the vector spherical harmonic synthesis
!     of the arrays br, bi, cr, and ci and stores the result in the
!     arrays v and w. v(i, j) and w(i, j) are the colatitudinal 
!     (measured from the north pole) and east longitudinal components
!     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given below at output parameters v, w.
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
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     ityp   = 0  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon.   
!
!            = 1  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon. the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 2  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 3  v is symmetric and w is antisymmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 4  v is symmetric and w is antisymmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 5  v is symmetric and w is antisymmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 6  v is antisymmetric and w is symmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 7  v is antisymmetric and w is symmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 8  v is antisymmetric and w is symmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!
!     nt     the number of syntheses.  in the program that calls vhses, 
!            the arrays v, w, br, bi, cr, and ci can be three dimensional
!            in which case multiple syntheses will be performed.
!            the third index is the synthesis index which assumes the 
!            values k=1, ..., nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhaes. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhses. jdvw must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given
!            below at the discription of output parameters v and w.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhses. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhses. ndab must be at
!            least nlat.
!
!     wvhses an array which must be initialized by subroutine vhsesi.
!            once initialized, wvhses can be used repeatedly by vhses
!            as long as nlon and nlat remain unchanged.  wvhses must
!            not be altered between calls of vhses.
!
!     lvhses the dimension of the array wvhses as it appears in the
!            program that calls vhses. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhses must be at least
!
!                 l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhses. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!
!                       (2*nt+1)*nlat*nlon
!
!            if ityp .gt. 2 then lwork must be at least
!
!                        (2*nt+1)*l2*nlon 
!
!     **************************************************************
!
!     output parameters
!
!     v, w    two or three dimensional arrays (see input parameter nt)
!            in which the synthesis is stored. v is the colatitudinal
!            component and w is the east longitudinal component. 
!            v(i, j), w(i, j) contain the components at colatitude
!            theta(i) = (i-1)*pi/(nlat-1) and longitude phi(j) =
!            (j-1)*2*pi/nlon. the index ranges are defined above at
!            the input parameter ityp. v and w are computed from the 
!            formulas given below
!
!
!     define
!
!     1.  theta is colatitude and phi is east longitude
!
!     2.  the normalized associated legendre funnctions
!
!         pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
!                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
!                        factorial(n)) times the (n+m)th derivative
!                        of (x**2-1)**n with respect to x=cos(theta)
!
!     3.  vbar(m, n, theta) = the derivative of pbar(m, n, theta) with
!                           respect to theta divided by the square
!                           root of n(n+1).
!
!         vbar(m, n, theta) is more easily computed in the form
!
!         vbar(m, n, theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1, n, theta)
!         -sqrt((n-m)*(n+m+1))*pbar(m+1, n, theta))/(2*sqrt(n*(n+1)))
!
!     4.  wbar(m, n, theta) = m/(sin(theta))*pbar(m, n, theta) divided
!                           by the square root of n(n+1).
!
!         wbar(m, n, theta) is more easily computed in the form
!
!         wbar(m, n, theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1))
!         *pbar(m-1, n-1, theta)+sqrt((n-m)*(n-m-1))*pbar(m+1, n-1, theta))
!         /(2*sqrt(n*(n+1)))
!
!
!    the colatitudnal dependence of the normalized surface vector
!                spherical harmonics are defined by
!
!     5.    bbar(m, n, theta) = (vbar(m, n, theta), i*wbar(m, n, theta))
!
!     6.    cbar(m, n, theta) = (i*wbar(m, n, theta), -vbar(m, n, theta))
!
!
!    the coordinate to index mappings 
!
!     7.   theta(i) = (i-1)*pi/(nlat-1) and phi(j) = (j-1)*2*pi/nlon
!
!    
!     the maximum (plus one) longitudinal wave number
!
!     8.     mmax = min(nlat, nlon/2) if nlon is even or
!            mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!    if we further define the output vector as
!
!     9.    h(i, j) = (v(i, j), w(i, j))
!
!    and the complex coefficients
!
!     10.   b(m, n) = cmplx(br(m+1, n+1), bi(m+1, n+1))
!
!     11.   c(m, n) = cmplx(cr(m+1, n+1), ci(m+1, n+1))
!
!
!    then for i=1, ..., nlat and  j=1, ..., nlon
!
!        the expansion for real h(i, j) takes the form
!
!     h(i, j) = the sum from n=1 to n=nlat-1 of the real part of
!
!         .5*(b(0, n)*bbar(0, n, theta(i))+c(0, n)*cbar(0, n, theta(i)))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
!     n=nlat-1 of the real part of
!
!              b(m, n)*bbar(m, n, theta(i))*exp(i*m*phi(j))
!             +c(m, n)*cbar(m, n, theta(i))*exp(i*m*phi(j))
!
!   *************************************************************
!
!   in terms of real variables this expansion takes the form
!
!             for i=1, ..., nlat and  j=1, ..., nlon
!
!     v(i, j) = the sum from n=1 to n=nlat-1 of
!
!               .5*br(1, n+1)*vbar(0, n, theta(i))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
!     n=nlat-1 of the real part of
!
!       (br(m+1, n+1)*vbar(m, n, theta(i))-ci(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *cos(m*phi(j))
!      -(bi(m+1, n+1)*vbar(m, n, theta(i))+cr(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *sin(m*phi(j))
!
!    and for i=1, ..., nlat and  j=1, ..., nlon
!
!     w(i, j) = the sum from n=1 to n=nlat-1 of
!
!              -.5*cr(1, n+1)*vbar(0, n, theta(i))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
!     n=nlat-1 of the real part of
!
!      -(cr(m+1, n+1)*vbar(m, n, theta(i))+bi(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *cos(m*phi(j))
!      +(ci(m+1, n+1)*vbar(m, n, theta(i))-br(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *sin(m*phi(j))
!
!
!      br(m+1, nlat), bi(m+1, nlat), cr(m+1, nlat), and ci(m+1, nlat) are
!      assumed zero for m even.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of ityp
!            = 4  error in the specification of nt
!            = 5  error in the specification of idvw
!            = 6  error in the specification of jdvw
!            = 7  error in the specification of mdab
!            = 8  error in the specification of ndab
!            = 9  error in the specification of lvhses
!            = 10 error in the specification of lwork
!
! ************************************************************
!
!     subroutine vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, 
!    +                  ldwork, ierror)
!
!     subroutine vhsesi initializes the array wvhses which can then be
!     used repeatedly by subroutine vhses until nlat or nlon is changed.
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
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhses the dimension of the array wvhses as it appears in the
!            program that calls vhses. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhses must be at least
!
!                  l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhses. lwork must be at least
!
!              3*(max(l1-2, 0)*(nlat+nlat-l1-1))/2+5*l2*nlat
!
!     dwork  an unsaved real work space
!
!     ldwork the length of the array dwork as it appears in the
!            program that calls vhsesi.  ldwork must be at least
!            2*(nlat+1)
!
!
!
!     output parameters
!
!     wvhses an array which is initialized for use by subroutine vhses.
!            once initialized, wvhses can be used repeatedly by vhses
!            as long as nlat or nlon remain unchanged.  wvhses must not
!            be altered between calls of vhses.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhses
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
!
subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    mdab, ndab, wvhses, lvhses, work, lwork, ierror)
    implicit none
    real :: bi
    real :: br
    real :: ci
    real :: cr
    integer :: idv
    integer :: idvw
    integer :: idz
    integer :: ierror
    integer :: imid
    integer :: ist
    integer :: ityp
    integer :: iw1
    integer :: iw2
    integer :: iw3
    integer :: iw4
    integer :: jdvw
    integer :: jw1
    integer :: jw2
    integer :: lnl
    integer :: lvhses
    integer :: lwork
    integer :: lzimn
    integer :: mdab
    integer :: mmax
    integer :: ndab
    integer :: nlat
    integer :: nlon
    integer :: nt
    real :: v
    real :: w
    real :: work
    real :: wvhses
    dimension v(idvw, jdvw, 1), w(idvw, jdvw, 1), br(mdab, ndab, 1), &
        bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
        work(1), wvhses(1)

    imid = (nlat+1)/2
    mmax = min(nlat, (nlon+1)/2)
    idz = (mmax*(2*nlat-mmax+1))/2
    lzimn = idz*imid

    if (ityp > 2) then
        idv = imid
    else
        idv = nlat
    end if

    lnl = nt*idv*nlon

    if (ityp <= 2) then
        ist = imid
    else
        ist = 0
    end if

    !
    !==> Check validity of input arguments
    !
    if (nlat < 3) then
        ierror = 1
        return
    else if (nlon < 1) then
        ierror = 2
        return
    else if (ityp < 0 .or. ityp > 8) then
        ierror = 3
        return
    else if (nt < 0) then
        ierror = 4
        return
    else if ((ityp <= 2 .and. idvw < nlat) .or. (ityp > 2 .and. idvw < imid)) then
        ierror = 5
        return
    else if (jdvw < nlon) then
        ierror = 6
        return
    else if (mdab < mmax) then
        ierror = 7
        return
    else if (ndab < nlat) then
        ierror = 8
        return
    else if (lvhses < 2*lzimn+nlon+15) then
        ierror = 9
        return
    else if (lwork < 2*lnl+idv*nlon) then
        ierror = 10
        return
    else
        ierror = 0
    end if

    iw1 = ist+1
    iw2 = lnl+1
    iw3 = iw2+ist
    iw4 = iw2+lnl
    jw1 = lzimn+1
    jw2 = jw1+lzimn

    call vhses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
        br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
        work(iw4), idz, wvhses, wvhses(jw1), wvhses(jw2))

contains

    subroutine vhses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)
        implicit none
        real :: bi
        real :: br
        real :: ci
        real :: cr
        integer :: i
        integer :: idv
        integer :: idvw
        integer :: idz
        integer :: imid
        integer :: imm1
        integer :: ityp

        integer :: j
        integer :: jdvw
        integer :: k
        integer :: m
        integer :: mb
        integer :: mdab
        integer :: mlat
        integer :: mlon
        integer :: mmax
        integer :: mn
        integer :: mp1
        integer :: mp2
        integer :: ndab
        integer :: ndo1
        integer :: ndo2
        integer :: nlat
        integer :: nlon
        integer :: nlp1
        integer :: np1
        integer :: nt
        real :: v
        real :: vb
        real :: ve
        real :: vo
        real :: w
        real :: wb
        real :: we
        real :: wo
        real :: work
        real :: wrfft
        dimension v(idvw, jdvw, 1), w(idvw, jdvw, 1), br(mdab, ndab, 1), &
            bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
            ve(idv, nlon, 1), vo(idv, nlon, 1), we(idv, nlon, 1), &
            wo(idv, nlon, 1), work(1), wrfft(1), &
            vb(imid, 1), wb(imid, 1)

        nlp1 = nlat+1
        mlat = mod(nlat, 2)
        mlon = mod(nlon, 2)
        mmax = min(nlat, (nlon+1)/2)

        if (mlat /= 0) then
            imm1 = imid-1
        else
            imm1 = imid
        end if

        do k=1, nt
            do j=1, nlon
                do i=1, idv
                    ve(i, j, k) = 0.
                    we(i, j, k) = 0.
                end do
            end do
        end do

        if (mlat /= 0) then
            ndo1 = nlat-1
        else
            ndo1 = nlat
        end if

        if (mlat == 0) then
            ndo2 = nlat-1
        else
            ndo2 = nlat
        end if

        case_block: block
            select case (ityp)
                case (0)
                    !
                    !     case ityp=0   no symmetries
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=2, ndo2, 2
                            do i=1, imid
                                ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    do k=1, nt
                        do np1=3, ndo1, 2
                            do i=1, imm1
                                vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then
                            do k=1, nt
                                do np1=mp1, ndo1, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do

                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            -ci(mp1, np1, k)*wb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +cr(mp1, np1, k)*wb(imid, mn)
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -bi(mp1, np1, k)*wb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            +br(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if

                        if (mp2 <= ndo2) then
                            do k=1, nt
                                do np1=mp2, ndo2, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            +br(mp1, np1, k)*vb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +bi(mp1, np1, k)*vb(imid, mn)
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -cr(mp1, np1, k)*vb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            -ci(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (1)
                    !
                    !     case ityp=1   no symmetries,  cr and ci equal zero
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=2, ndo2, 2
                            do i=1, imid
                                ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    do k=1, nt
                        do np1=3, ndo1, 2
                            do i=1, imm1
                                vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1
                        if (mp1 <= ndo1) then
                            do k=1, nt
                                do np1=mp1, ndo1, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -bi(mp1, np1, k)*wb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            +br(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if

                        if (mp2 <= ndo2) then
                            do k=1, nt
                                do np1=mp2, ndo2, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            +br(mp1, np1, k)*vb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +bi(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (2)
                    !
                    !     case ityp=2   no symmetries,  br and bi are equal to zero
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=2, ndo2, 2
                            do i=1, imid
                                we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do np1=3, ndo1, 2
                            do i=1, imm1
                                wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1
                        if (mp1 <= ndo1) then
                            do k=1, nt
                                do np1=mp1, ndo1, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    end do

                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            -ci(mp1, np1, k)*wb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +cr(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if

                        if (mp2 <= ndo2) then
                            do k=1, nt
                                do np1=mp2, ndo2, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -cr(mp1, np1, k)*vb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            -ci(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (3)
                     !
                     !     case ityp=3   v even,  w odd
                     !
                     !     case m = 0
                     !
                    do k=1, nt
                        do np1=2, ndo2, 2
                            do i=1, imid
                                ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do np1=3, ndo1, 2
                            do i=1, imm1
                                wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1
                        if (mp1 <= ndo1) then
                            do k=1, nt
                                do np1=mp1, ndo1, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            -ci(mp1, np1, k)*wb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +cr(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if

                        if (mp2 <= ndo2) then
                            do k=1, nt
                                do np1=mp2, ndo2, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            +br(mp1, np1, k)*vb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +bi(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (4)
                    !
                    !     case ityp=4   v even,  w odd, and both cr and ci equal zero
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=2, ndo2, 2
                            do i=1, imid
                                ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1
                        if (mp2 <= ndo2) then
                            do k=1, nt
                                do np1=mp2, ndo2, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            +br(mp1, np1, k)*vb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +bi(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (5)
                    !
                    !     case ityp=5   v even,  w odd,     br and bi equal zero
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=3, ndo1, 2
                            do i=1, imm1
                                wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then
                            do k=1, nt
                                do np1=mp1, ndo1, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            -ci(mp1, np1, k)*wb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +cr(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (6)
                    !
                    !     case ityp=6   v odd  ,  w even
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=2, ndo2, 2
                            do i=1, imid
                                we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do np1=3, ndo1, 2
                            do i=1, imm1
                                vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1
                        if (mp1 <= ndo1) then
                            do k=1, nt
                                do np1=mp1, ndo1, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -bi(mp1, np1, k)*wb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            +br(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if

                        if (mp2 <= ndo2) then
                            do k=1, nt
                                do np1=mp2, ndo2, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -cr(mp1, np1, k)*vb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            -ci(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (7)
                    !
                    !     case ityp=7   v odd, w even   cr and ci equal zero
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=3, ndo1, 2
                            do i=1, imm1
                                vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1
                        if (mp1 <= ndo1) then
                            do k=1, nt
                                do np1=mp1, ndo1, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -bi(mp1, np1, k)*wb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            +br(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                case (8)
                    !
                    !     case ityp=8   v odd,  w even   br and bi equal zero
                    !
                    !     case m = 0
                    !
                    do k=1, nt
                        do np1=2, ndo2, 2
                            do i=1, imid
                                we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                            end do
                        end do
                    end do
                    !
                    !     case m = 1 through nlat-1
                    !
                    if (mmax < 2) exit case_block

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*(nlat-1)-(m*(m-1))/2
                        mp2 = mp1+1
                        if (mp2 <= ndo2) then
                            do k=1, nt
                                do np1=mp2, ndo2, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    end do
                                    if (mlat /= 0) then
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -cr(mp1, np1, k)*vb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            -ci(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
            end select
        end block case_block

        do k=1, nt
            call hrfftb(idv, nlon, ve(1, 1, k), idv, wrfft, work)
            call hrfftb(idv, nlon, we(1, 1, k), idv, wrfft, work)
        end do

        if (ityp <= 2) then
            do k=1, nt
                do j=1, nlon
                    do i=1, imm1
                        v(i, j, k) = .5*(ve(i, j, k)+vo(i, j, k))
                        w(i, j, k) = .5*(we(i, j, k)+wo(i, j, k))
                        v(nlp1-i, j, k) = .5*(ve(i, j, k)-vo(i, j, k))
                        w(nlp1-i, j, k) = .5*(we(i, j, k)-wo(i, j, k))
                    end do
                end do
            end do
        else
            do k=1, nt
                do j=1, nlon
                    do i=1, imm1
                        v(i, j, k) = .5*ve(i, j, k)
                        w(i, j, k) = .5*we(i, j, k)
                    end do
                end do
            end do
        end if

        if (mlat /= 0) then
            do k=1, nt
                do j=1, nlon
                    v(imid, j, k) = .5*ve(imid, j, k)
                    w(imid, j, k) = .5*we(imid, j, k)
                end do
            end do
        end if

    end subroutine vhses1

end subroutine vhses

subroutine vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, &
    ldwork, ierror)
    implicit none
    integer :: idz
    integer :: ierror
    integer :: imid
    integer :: iw1
    integer :: labc
    integer :: ldwork
    integer :: lvhses
    integer :: lwork
    integer :: lzimn
    integer :: mmax
    integer :: nlat
    integer :: nlon
    real :: work
    real :: wvhses
    dimension wvhses(lvhses), work(lwork)
    real dwork(ldwork)
    ierror = 1
    if (nlat < 3) return
    ierror = 2
    if (nlon < 1) return
    ierror = 3
    mmax = min(nlat, (nlon+1)/2)
    imid = (nlat+1)/2
    lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
    if (lvhses < lzimn+lzimn+nlon+15) return
    ierror = 4
    labc = 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
    if (lwork < 5*nlat*imid+labc) return
    ierror = 5
    if (ldwork < 2*(nlat+1)) return
    ierror = 0
    iw1 = 3*nlat*imid+1
    idz = (mmax*(nlat+nlat-mmax+1))/2

    call ves1(nlat, nlon, imid, wvhses, wvhses(lzimn+1), idz, work, work(iw1), dwork)
    call hrffti(nlon, wvhses(2*lzimn+1))

contains

    subroutine ves1(nlat, nlon, imid, vb, wb, idz, vin, wzvin, dwork)
        implicit none
        integer :: i
        integer :: i3
        integer :: idz
        integer :: imid
        integer :: m
        integer :: mmax
        integer :: mn
        integer :: mp1
        integer :: nlat
        integer :: nlon
        integer :: np1
        real :: vb
        real :: vin
        real :: wb
        real :: wzvin
        dimension vb(imid, *), wb(imid, *), vin(imid, nlat, 3), wzvin(*)
        real dwork(*)

        mmax = min(nlat, (nlon+1)/2)

        call vbinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call vbin (0, nlat, nlon, m, vin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                do i=1, imid
                    vb(i, mn) = vin(i, np1, i3)
                end do
            end do
        end do

        call wbinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call wbin (0, nlat, nlon, m, vin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                do i=1, imid
                    wb(i, mn) = vin(i, np1, i3)
                end do
            end do
        end do

    end subroutine ves1

end subroutine vhsesi
