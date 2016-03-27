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
! ... file gradgc.f
!
!     this file includes documentation and code for
!     subroutine gradgc         i
!
! ... files which must be loaded with gradgc.f
!
!     sphcom.f, hrfft.f, shagc.f, vhsgc.f
!
!     subroutine gradgc(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, 
!    +                  wvhsgc, lvhsgc, work, lwork, ierror)
!
!     given the scalar spherical harmonic coefficients a and b, precomputed
!     by subroutine shagc for a scalar field sf, subroutine gradgc computes
!     an irrotational vector field (v, w) such that
!
!           gradient(sf) = (v, w).
!
!     v is the colatitudinal and w is the east longitudinal component
!     of the gradient.  i.e., 
!
!            v(i, j) = d(sf(i, j))/dtheta
!
!     and
!
!            w(i, j) = 1/sint*d(sf(i, j))/dlambda
!
!     at the gaussian colatitude point theta(i) (see nlat as input
!     parameter) and longitude lambda(j) = (j-1)*2*pi/nlon where
!     where sint = sin(theta(i)).  required associated legendre polynomials
!     are recomputed rather than stored as they are in subroutine gradgs. this
!     saves storage (compare lsav with lsav in gradgs) but increases
!     computational requirements.
!
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
!     isym   this has the same value as the isym that was input to
!            subroutine shagc to compute the arrays a and b from the
!            scalar field sf.  isym determines whether (v, w) are
!            computed on the full or half sphere as follows:
!
!      = 0
!
!           sf is not symmetric about the equator. in this case
!           the vector field (v, w) is computed on the entire sphere.
!           i.e., in the arrays  v(i, j), w(i, j) for i=1, ..., nlat and
!           j=1, ..., nlon.
!
!      = 1
!
!           sf is antisymmetric about the equator. in this case w is
!           antisymmetric and v is symmetric about the equator. w
!           and v are computed on the northern hemisphere only.  i.e., 
!           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!           and j=1, ..., nlon.  if nlat is even they are computed for
!           i=1, ..., nlat/2 and j=1, ..., nlon.
!
!      = 2
!
!           sf is symmetric about the equator. in this case w is
!           symmetric and v is antisymmetric about the equator. w
!           and v are computed on the northern hemisphere only.  i.e., 
!           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!           and j=1, ..., nlon.  if nlat is even they are computed for
!           i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            the arrays a, b, v, and w can be three dimensional corresponding
!            to an indexed multiple array sf.  in this case, multiple
!            vector synthesis will be performed to compute each vector
!            field.  the third index for a, b, v, and w is the synthesis
!            index which assumes the values k = 1, ..., nt.  for a single
!            synthesis set nt = 1.  the description of the remaining
!            parameters is simplified by assuming that nt=1 or that a, b, v, 
!            and w are two dimensional arrays.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls gradgc. if isym = 0 then idvw
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idvw must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls gradgc. jdvw must be at least nlon.
!
!     a, b    two or three dimensional arrays (see input parameter nt)
!            that contain scalar spherical harmonic coefficients
!            of the scalar field array sf as computed by subroutine shagc.
!     ***    a, b must be computed by shagc prior to calling gradgc.
!
!     mdab   the first dimension of the arrays a and b as it appears in
!            the program that calls gradgc (and shagc). mdab must be at
!            least min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears in
!            the program that calls gradgc (and shagc). ndab must be at
!            least nlat.
!
!
!     wvhsgc an array which must be initialized by subroutine vhsgci.
!            once initialized, 
!            wvhsgc can be used repeatedly by gradgc as long as nlon
!            and nlat remain unchanged.  wvhsgc must not be altered
!            between calls of gradgc.
!
!
!     lvhsgc the dimension of the array wvhsgc as it appears in the
!            program that calls gradgc. Let
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2                  if nlat is even or
!               l2 = (nlat+1)/2              if nlat is odd
!
!            lvhsgc must be at least
!
!               4*nlat*l2+3*max(l1-2, 0)*(nlat+nlat-l1-1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls gradgc. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2                  if nlat is even or
!               l2 = (nlat+1)/2              if nlat is odd
!
!            if isym = 0 then lwork must be at least
!
!                nlat*(2*nt*nlon+max(6*l2, nlon)) + nlat*(2*l1*nt+1)
!
!            if isym = 1 or 2 then lwork must be at least
!
!                l2*(2*nt*nlon+max(6*nlat, nlon)) + nlat*(2*l1*nt+1)
!
!
!
!     **************************************************************
!
!     output parameters
!
!
!     v, w   two or three dimensional arrays (see input parameter nt) that
!           contain an irrotational vector field such that the gradient of
!           the scalar field sf is (v, w).  w(i, j) is the east longitude
!           component and v(i, j) is the colatitudinal component of velocity
!           at gaussian colatitude and longitude lambda(j) = (j-1)*2*pi/nlon
!           the indices for v and w are defined at the input parameter
!           isym.  the vorticity of (v, w) is zero.  note that any nonzero
!           vector field on the sphere will be multiple valued at the poles
!           [reference swarztrauber].
!
!
!  ierror   = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of idvw
!           = 6  error in the specification of jdvw
!           = 7  error in the specification of mdab
!           = 8  error in the specification of ndab
!           = 9  error in the specification of lvhsgc
!           = 10 error in the specification of lwork
! **********************************************************************
!                                                                              
!   
subroutine gradgc(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
    wvhsgc, lvhsgc, work, lwork, ierror)
    dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt)
    dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
    dimension wvhsgc(lvhsgc), work(lwork)
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
    if((isym==0 .and. idvw<nlat) .or. &
        (isym/=0 .and. idvw<imid)) return
    ierror = 6
    if(jdvw < nlon) return
    ierror = 7
    mmax = min(nlat, (nlon+1)/2)
    if(mdab < min(nlat, (nlon+2)/2)) return
    ierror = 8
    if(ndab < nlat) return
    ierror = 9
    !
    !     verify minimum saved work space length
    !
    l1 = min(nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lwmin =   4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
    if (lvhsgc < lwmin) return
    ierror = 10
    !
    !     verify minimum unsaved work space length
    !
    if (isym == 0) then
        lwkmin = nlat*(2*nt*nlon+max(6*l2, nlon)+2*l1*nt+1)
    else
        lwkmin = l2*(2*nt*nlon+max(6*nlat, nlon)) + nlat*(2*l1*nt+1)
    end if
    if (lwork < lwkmin) return

    ierror = 0
    !
    !     set work space pointers
    !
    mn = mmax*nlat*nt
    ibr = 1
    ibi = ibr + mn
    is = ibi + mn
    iwk = is + nlat
    liwk = lwork-2*mn-nlat
    call gradgc1(nlat, nlon, isym, nt, v, w, idvw, jdvw, work(ibr), work(ibi), &
        mmax, work(is), mdab, ndab, a, b, wvhsgc, lvhsgc, work(iwk), liwk, &
        ierror)
    return
end subroutine gradgc
subroutine gradgc1(nlat, nlon, isym, nt, v, w, idvw, jdvw, br, bi, mmax, &
    sqnn, mdab, ndab, a, b, wvhsgc, lvhsgc, wk, lwk, ierror)
    dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt)
    dimension br(mmax, nlat, nt), bi(mmax, nlat, nt), sqnn(nlat)
    dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
    dimension wvhsgc(lvhsgc), wk(lwk)
    !
    !     preset coefficient multiplyers in vector
    !
    do 1 n=2, nlat
        fn = real(n-1)
        sqnn(n) = sqrt(fn*(fn+1.))
1   continue
    !
    !     compute multiple vector fields coefficients
    !
    do 2 k=1, nt
        !
        !     preset br, bi to 0.0
        !
        do 3 n=1, nlat
            do 4 m=1, mmax
                br(m, n, k) = 0.0
                bi(m, n, k) = 0.0
4           continue
3       continue
        !
        !     compute m=0 coefficients
        !
        do 5 n=2, nlat
            br(1, n, k) = sqnn(n)*a(1, n, k)
            bi(1, n, k) = sqnn(n)*b(1, n, k)
5       continue
        !
        !     compute m>0 coefficients
        !
        do 6 m=2, mmax
            do 7 n=m, nlat
                br(m, n, k) = sqnn(n)*a(m, n, k)
                bi(m, n, k) = sqnn(n)*b(m, n, k)
7           continue
6       continue
2   continue
    !
    !     set ityp for irrotational vector synthesis to compute gradient
    !
    if (isym==0) then
        ityp = 1
    else if (isym==1) then
        ityp = 4
    else if (isym==2) then
        ityp = 7
    end if
    !
    !     vector sythesize br, bi into (v, w) (cr, ci are dummy variables)
    !
    call vhsgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mmax, nlat, wvhsgc, lvhsgc, wk, lwk, ierror)
    return
end subroutine gradgc1
