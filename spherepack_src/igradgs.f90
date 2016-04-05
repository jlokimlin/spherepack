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
! ... file igradgs.f
!
!     this file includes documentation and code for
!     subroutine igradgs         i
!
! ... files which must be loaded with igradgs.f
!
!     sphcom.f, hrfft.f, shsgs.f, vhags.f
!
!     subroutine igradgs(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, 
!    +                   wshsgs, lshsgs, work, lwork, ierror)
!
!     let br, bi, cr, ci be the vector spherical harmonic coefficients
!     precomputed by vhags for a vector field (v, w).  let (v', w') be
!     the irrotational component of (v, w) (i.e., (v', w') is generated
!     by assuming cr, ci are zero and synthesizing br, bi with vhsgs).
!     then subroutine igradgs computes a scalar field sf such that
!
!            gradient(sf) = (v', w').
!
!     i.e., 
!
!            v'(i, j) = d(sf(i, j))/dtheta          (colatitudinal component of
!                                                 the gradient)
!     and
!
!            w'(i, j) = 1/sint*d(sf(i, j))/dlambda  (east longitudinal component
!                                                 of the gradient)
!
!     at the gaussian colatitude theta(i) (see nlat as input parameter)
!     and longitude lambda(j) = (j-1)*2*pi/nlon where sint = sin(theta(i)).
!
!     note:  for an irrotational vector field (v, w), subroutine igradgs
!     computes a scalar field whose gradient is (v, w).  in ay case, 
!     subroutine igradgs "inverts" the gradient subroutine gradgs.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are computed
!            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
!            if nlat is odd the equator will be included as the grid point
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater than
!            3.  the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   a parameter which determines whether the scalar field sf is
!            computed on the full or half sphere as follows:
!
!      = 0
!
!            the symmetries/antsymmetries described in isym=1, 2 below
!            do not exist in (v, w) about the equator.  in this case sf
!            is neither symmetric nor antisymmetric about the equator.
!            sf is computed on the entire sphere.  i.e., in the array
!            sf(i, j) for i=1, ..., nlat and  j=1, ..., nlon
!
!      = 1
!
!            w is antisymmetric and v is symmetric about the equator.
!            in this case sf is antisymmetyric about the equator and
!            is computed for the northern hemisphere only.  i.e., 
!            if nlat is odd sf is computed in the array sf(i, j) for
!            i=1, ..., (nlat+1)/2 and for j=1, ..., nlon.  if nlat is even
!            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!      = 2
!
!            w is symmetric and v is antisymmetric about the equator.
!            in this case sf is symmetyric about the equator and
!            is computed for the northern hemisphere only.  i.e., 
!            if nlat is odd sf is computed in the array sf(i, j) for
!            i=1, ..., (nlat+1)/2 and for j=1, ..., nlon.  if nlat is even
!            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            the arrays br, bi, and sf can be three dimensional corresponding
!            to an indexed multiple vector field (v, w).  in this case, 
!            multiple scalar synthesis will be performed to compute each
!            scalar field.  the third index for br, bi, and sf is the synthesis
!            index which assumes the values k = 1, ..., nt.  for a single
!            synthesis set nt = 1.  the description of the remaining
!            parameters is simplified by assuming that nt=1 or that br, bi, 
!            and sf are two dimensional arrays.
!
!     isf    the first dimension of the array sf as it appears in
!            the program that calls igradgs. if isym = 0 then isf
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then isf must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then isf must be at least (nlat+1)/2.
!
!     jsf    the second dimension of the array sf as it appears in
!            the program that calls igradgs. jsf must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhags.
!     ***    br, bi must be computed by vhags prior to calling igradgs.
!
!     mdb    the first dimension of the arrays br and bi as it appears in
!            the program that calls igradgs (and vhags). mdb must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays br and bi as it appears in
!            the program that calls igradgs (and vhags). ndb must be at
!            least nlat.
!
!
!  wshsgs    an array which must be initialized by subroutine igradgsi
!            (or equivalently by subroutine shsesi).  once initialized, 
!            wshsgs can be used repeatedly by igradgs as long as nlon
!            and nlat remain unchanged.  wshsgs must not be altered
!            between calls of igradgs.
!
!
!  lshsgs    the dimension of the array wshsgs as it appears in the
!            program that calls igradgs. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd.
!
!
!            then lshsgs must be greater than or equal to
!
!               nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls igradgs. define
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            if isym = 0 lwork must be greater than or equal to
!
!               nlat*((nt+1)*nlon+2*nt*l1+1)
!
!            if isym > 0 lwork must be greater than or equal to
!
!               (nt+1)*l2*nlon+nlat*(2*nt*l1+1)
!
!
!
!     **************************************************************
!
!     output parameters
!
!
!     sf    a two or three dimensional array (see input parameter nt) that
!           contain a scalar field whose gradient is the irrotational
!           component of the vector field (v, w).  the vector spherical
!           harmonic coefficients br, bi were precomputed by subroutine
!           vhags.  sf(i, j) is given at the gaussian colatitude theta(i)
!           and longitude lambda(j) = (j-1)*2*pi/nlon.  the index ranges
!           are defined at input parameter isym.
!
!
!  ierror   = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of isf
!           = 6  error in the specification of jsf
!           = 7  error in the specification of mdb
!           = 8  error in the specification of ndb
!           = 9  error in the specification of lshsgs
!           = 10 error in the specification of lwork
!
! **********************************************************************
!   
subroutine igradgs(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, &
    wshsgs, lshsgs, work, lwork, ierror)
    dimension sf(isf, jsf, nt)
    dimension br(mdb, ndb, nt), bi(mdb, ndb, nt)
    dimension wshsgs(lshsgs), work(lwork)
    !
    !     check input parameters
    !
    ierror = 1
    if(nlat < 3) return
    ierror = 2
    if(nlon < 4) return
    ierror = 3
    if(isym<0 .or. isym>2) return
    ierror = 4
    if(nt < 0) return
    ierror = 5
    imid = (nlat+1)/2
    if((isym==0 .and. isf<nlat) .or. &
        (isym/=0 .and. isf<imid)) return
    ierror = 6
    if(jsf < nlon) return
    ierror = 7
    mmax = min(nlat, (nlon+2)/2)
    if(mdb < min(nlat, (nlon+1)/2)) return
    ierror = 8
    if(ndb < nlat) return
    ierror = 9
    !
    !     verify saved work space length
    !
    l2 = (nlat+mod(nlat, 2))/2
    l1 = min((nlon+2)/2, nlat)
    lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    if(lshsgs<lp) return
    ierror = 10
    !
    !     set minimum and verify unsaved work space
    !
    ls = nlat
    if(isym > 0) ls = imid
    nln = nt*ls*nlon
    !
    !     set first dimension for a, b (as requried by shses)
    !
    mab = min(nlat, nlon/2+1)
    mn = mab*nlat*nt
    lwkmin = nln+ls*nlon+2*mn+nlat
    if (lwork < lwkmin) return
    ierror = 0
    !
    !     set work space pointers
    !
    ia = 1
    ib = ia + mn
    is = ib + mn
    iwk = is + nlat
    liwk = lwork-2*mn-nlat

    call igrdgs1(nlat, nlon, isym, nt, sf, isf, jsf, work(ia), work(ib), mab, &
        work(is), mdb, ndb, br, bi, wshsgs, lshsgs, work(iwk), liwk, ierror)

end subroutine igradgs



subroutine igrdgs1(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, &
    sqnn, mdb, ndb, br, bi, wsav, lsav, wk, lwk, ierror)
    dimension sf(isf, jsf, nt)
    dimension br(mdb, ndb, nt), bi(mdb, ndb, nt), sqnn(nlat)
    dimension a(mab, nlat, nt), b(mab, nlat, nt)
    dimension wsav(lsav), wk(lwk)
    !
    !     preset coefficient multiplyers in vector
    !
    do n=2, nlat
        fn = real(n-1)
        sqnn(n) = 1.0/sqrt(fn*(fn + 1.0))
    end do
    !
    !     set upper limit for vector m subscript
    !
    mmax = min(nlat, (nlon+1)/2)
    !
    !     compute multiple scalar field coefficients
    !
    do k=1, nt
        !
        !     preset to 0.0
        !
        do n=1, nlat
            do m=1, mab
                a(m, n, k) = 0.0
                b(m, n, k) = 0.0
            end do
        end do
        !
        !     compute m=0 coefficients
        !
        do n=2, nlat
            a(1, n, k) = br(1, n, k)*sqnn(n)
            b(1, n, k)= bi(1, n, k)*sqnn(n)
        end do
        !
        !     compute m>0 coefficients
        !
        do m=2, mmax
            do n=m, nlat
                a(m, n, k) = sqnn(n)*br(m, n, k)
                b(m, n, k) = sqnn(n)*bi(m, n, k)
            end do
        end do
    end do
    !
    !     scalar sythesize a, b into sf
    !
    call shsgs(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, nlat, wsav, &
        lsav, wk, lwk, ierror)

end subroutine igrdgs1
