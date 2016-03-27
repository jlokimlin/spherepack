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
! ... file igradec.f
!
!     this file includes documentation and code for
!     subroutine igradec         i
!
! ... files which must be loaded with igradec.f
!
!     sphcom.f, hrfft.f, shsec.f, vhaec.f
!
!     subroutine igradec(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, 
!    +                   wshsec, lshsec, work, lwork, ierror)
!
!     let br, bi, cr, ci be the vector spherical harmonic coefficients
!     precomputed by vhaec for a vector field (v, w).  let (v', w') be
!     the irrotational component of (v, w) (i.e., (v', w') is generated
!     by assuming cr, ci are zero and synthesizing br, bi with vhsec).
!     then subroutine igradec computes a scalar field sf such that
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
!     at colatitude
!
!            theta(i) = (i-1)*pi/(nlat-1)
!
!     and longitude
!
!            lambda(j) = (j-1)*2*pi/nlon
!
!     where sint = sin(theta(i)).  required associated legendre polynomials
!     are recomputed rather than stored as they are in subroutine igrades. this
!     saves storage (compare lshsec and lshses in igrades) but increases
!     computational requirements.
!
!     note:  for an irrotational vector field (v, w), subroutine igradec
!     computes a scalar field whose gradient is (v, w).  in ay case, 
!     subroutine igradec "inverts" the gradient subroutine gradec.
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
!            the program that calls igradec. if isym = 0 then isf
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then isf must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then isf must be at least (nlat+1)/2.
!
!     jsf    the second dimension of the array sf as it appears in
!            the program that calls igradec. jsf must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhaec.
!     ***    br, bi must be computed by vhaec prior to calling igradec.
!
!     mdb    the first dimension of the arrays br and bi as it appears in
!            the program that calls igradec (and vhaec). mdb must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays br and bi as it appears in
!            the program that calls igradec (and vhaec). ndb must be at
!            least nlat.
!
!
!  wshsec    an array which must be initialized by subroutine shseci.
!            once initialized, 
!            wshsec can be used repeatedly by igradec as long as nlon
!            and nlat remain unchanged.  wshsec must not be altered
!            between calls of igradec.
!
!
!  lshsec    the dimension of the array wshsec as it appears in the
!            program that calls igradec. define
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
!            then lshsec must be greater than or equal to
!
!               2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls igradec. define
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!               l1 = min(nlat, (nlon+2)/2)     if nlon is even or
!               l1 = min(nlat, (nlon+1)/2)     if nlon is odd
!
!            if isym is zero then lwork must be at least
!
!               nlat*(nt*nlon+max(3*l2, nlon)+2*nt*l1+1)
!
!            if isym is not zero then lwork must be at least
!
!               l2*(nt*nlon+max(3*nlat, nlon))+nlat*(2*nt*l1+1)
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
!           vhaec.  sf(i, j) is given at the gaussian colatitude theta(i)
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
!           = 9  error in the specification of lshsec
!           = 10 error in the specification of lwork
!
! **********************************************************************
!   
subroutine igradec(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, &
    wshsec, lshsec, work, lwork, ierror)
    dimension sf(isf, jsf, nt)
    dimension br(mdb, ndb, nt), bi(mdb, ndb, nt)
    dimension wshsec(lshsec), work(lwork)
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
    imid = (nlat+1)/2
    lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
    !
    !     verify saved work space length
    !
    l1 = min(nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwkmin=2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
    if (lshsec < lwkmin) return
    ierror = 10
    !
    !     set minimum and verify unsaved work space
    !
    ls = nlat
    if(isym > 0) ls = imid
    nln = nt*ls*nlon
    !
    !     set first dimension for a, b (as requried by shsec)
    !
    mab = min(nlat, nlon/2+1)
    mn = mab*nlat*nt
    !     lwkmin = nln+ls*nlon+2*mn+nlat
    if (isym == 0) then
        lwkmin = nlat*(nt*nlon+max(3*l2, nlon)+2*nt*l1+1)
    else
        lwkmin = l2*(nt*nlon+max(3*nlat, nlon))+nlat*(2*nt*l1+1)
    end if
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
    call igrdec1(nlat, nlon, isym, nt, sf, isf, jsf, work(ia), work(ib), mab, &
        work(is), mdb, ndb, br, bi, wshsec, lshsec, work(iwk), liwk, ierror)
    return
end subroutine igradec

subroutine igrdec1(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, &
    sqnn, mdb, ndb, br, bi, wshsec, lshsec, wk, lwk, ierror)
    dimension sf(isf, jsf, nt)
    dimension br(mdb, ndb, nt), bi(mdb, ndb, nt), sqnn(nlat)
    dimension a(mab, nlat, nt), b(mab, nlat, nt)
    dimension wshsec(lshsec), wk(lwk)
    !
    !     preset coefficient multiplyers in vector
    !
    do 1 n=2, nlat
        fn = real(n-1)
        sqnn(n) = 1.0/sqrt(fn*(fn+1.))
1   continue
    !
    !     set upper limit for vector m subscript
    !
    mmax = min(nlat, (nlon+1)/2)
    !
    !     compute multiple scalar field coefficients
    !
    do 2 k=1, nt
        !
        !     preset to 0.0
        !
        do 3 n=1, nlat
            do 4 m=1, mab
                a(m, n, k) = 0.0
                b(m, n, k) = 0.0
4           continue
3       continue
        !
        !     compute m=0 coefficients
        !
        do 5 n=2, nlat
            a(1, n, k) = br(1, n, k)*sqnn(n)
            b(1, n, k)= bi(1, n, k)*sqnn(n)
5       continue
        !
        !     compute m>0 coefficients
        !
        do 6 m=2, mmax
            do 7 n=m, nlat
                a(m, n, k) = sqnn(n)*br(m, n, k)
                b(m, n, k) = sqnn(n)*bi(m, n, k)
7           continue
6       continue
2   continue
    !
    !     scalar sythesize a, b into sf
    !
    call shsec(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, nlat, &
        wshsec, lshsec, wk, lwk, ierror)
    return
end subroutine igrdec1
