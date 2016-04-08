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
! ... file idivec.f
!
!     this file includes documentation and code for
!     subroutine idivgs          i
!
! ... files which must be loaded with idivgs.f
!
!     sphcom.f, hrfft.f, vhsgs.f, shags.f
!
!
!     subroutine idivgs(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, 
!    +                  wvhsgs, lvhsgs, work, lwork, pertrb, ierror)
!
!     given the scalar spherical harmonic coefficients a and b, precomputed
!     by subroutine shags for a scalar array divg, subroutine idivgs computes
!     an irrotational vector field (v, w) whose divergence is divg - pertrb.
!     w is the east longitude component and v is the colatitudinal component.
!     pertrb is a constant which must be subtracted from divg for (v, w) to
!     exist (see the description of pertrb below).  usually pertrb is zero
!     or small relative to divg.  the vorticity of (v, w) is the zero scalar
!     field.  v(i, j) and w(i, j) are the velocity components at the gaussian
!     colatitude theta(i) (see nlat) and longitude lambda(j)=(j-1)*2*pi/nlon.
!     the
!
!            divergence[v(i, j), w(i, j)]
!
!          = [d(w(i, j)/dlambda + d(sint*v(i, j))/dtheta]/sint
!
!          = divg(i, j) - pertrb
!
!     and
!
!            vorticity(v(i, j), w(i, j))
!
!         =  [dv/dlambda - d(sint*w)/dtheta]/sint
!
!         =  0.0
!
!     where sint = sin(theta(i)).
!
!     input parameters
!
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
!            subroutine shags to compute the arrays a and b from the
!            scalar field divg.  isym determines whether (v, w) are
!            computed on the full or half sphere as follows:
!
!      = 0
!
!           divg is not symmetric about the equator. in this case
!           the vector field (v, w) is computed on the entire sphere.
!           i.e., in the arrays  v(i, j), w(i, j) for i=1, ..., nlat and
!           j=1, ..., nlon.
!
!      = 1
!
!           divg is antisymmetric about the equator. in this case w is
!           antisymmetric and v is symmetric about the equator. w
!           and v are computed on the northern hemisphere only.  i.e., 
!           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!           and j=1, ..., nlon.  if nlat is even they are computed for
!           i=1, ..., nlat/2 and j=1, ..., nlon.
!
!      = 2
!
!           divg is symmetric about the equator. in this case w is
!           symmetric and v is antisymmetric about the equator. w
!           and v are computed on the northern hemisphere only.  i.e., 
!           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!           and j=1, ..., nlon.  if nlat is even they are computed for
!           i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!     nt     nt is the number of divergence and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            the arrays a, b, v, and w can be three dimensional and pertrb
!            can be one dimensional corresponding to an indexed multiple
!            array divg.  in this case, multiple vector synthesis will be
!            performed to compute each vector field.  the third index for
!            a, b, v, w and first for pertrb is the synthesis index which
!            assumes the values k = 1, ..., nt.  for a single synthesis set
!            nt = 1.  the description of the remaining parameters is
!            simplified by assuming that nt=1 or that a, b, v, w are two
!            dimensional and pertrb is a constant.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls idivgs. if isym = 0 then idvw
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idvw must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls idivgs. jdvw must be at least nlon.
!
!     a, b    two or three dimensional arrays (see input parameter nt)
!            that contain scalar spherical harmonic coefficients
!            of the divergence array divg as computed by subroutine shags.
!     ***    a, b must be computed by shags prior to calling idivgs.
!
!     mdab   the first dimension of the arrays a and b as it appears in
!            the program that calls idivgs (and shags). mdab must be at
!            least min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears in
!            the program that calls idivgs (and shags). ndab must be at
!            least nlat.
!
!
!  wvhsgs    an array which must be initialized by subroutine vhsgsi.
!            once initialized, 
!            wvhsgs can be used repeatedly by idivgs as long as nlon
!            and nlat remain unchanged.  wvhsgs must not be altered
!            between calls of idivgs.
!
!
!  lvhsgs    the dimension of the array wvhsgs as it appears in the
!            program that calls idivgs. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhsgs must be at least
!
!               l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls idivgs. define
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
!               (2*nt+1)*nlat*nlon + nlat*(2*nt*l1+1)
!
!            if isym = 1 or 2 then lwork must be at least
!
!               (2*nt+1)*l2*nlon + nlat*(2*nt*l1+1)
!
!     **************************************************************
!
!     output parameters
!
!
!     v, w   two or three dimensional arrays (see input parameter nt) that
!           contain an irrotational vector field whose divergence is
!           divg-pertrb at the guassian colatitude point theta(i) and
!           longitude point lambda(j)=(j-1)*2*pi/nlon.  w is the east
!           longitude component and v is the colatitudinal component.  the
!           indices for w and v are defined at the input parameter isym.
!           the curl or vorticity of (v, w) is the zero vector field.  note
!           that any nonzero vector field on the sphere will be multiple
!           valued at the poles [reference swarztrauber].
!
!   pertrb  a nt dimensional array (see input parameter nt and assume nt=1
!           for the description that follows).  divg - pertrb is a scalar
!           field which can be the divergence of a vector field (v, w).
!           pertrb is related to the scalar harmonic coefficients a, b
!           of divg (computed by shags) by the formula
!
!                pertrb = a(1, 1)/(2.*sqrt(2.))
!
!
!
!           the unperturbed scalar field divg can be the divergence of a
!           vector field only if a(1, 1) is zero.  if a(1, 1) is nonzero
!           (flagged by pertrb nonzero) then subtracting pertrb from
!           divg yields a scalar field for which a(1, 1) is zero.
!
!    ierror = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of idvw
!           = 6  error in the specification of jdvw
!           = 7  error in the specification of mdab
!           = 8  error in the specification of ndab
!           = 9  error in the specification of lvhsgs
!           = 10 error in the specification of lwork
! **********************************************************************
!                                                                              
!   
subroutine idivgs(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
    wvhsgs, lvhsgs, work, lwork, pertrb, ierror)
    dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt), pertrb(nt)
    dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
    dimension wvhsgs(lvhsgs), work(lwork)
    !
    !     check input parameters
    !
    ierror = 1
    if (nlat < 3) return
    ierror = 2
    if (nlon < 4) return
    ierror = 3
    if (isym<0 .or. isym>2) return
    ierror = 4
    if (nt < 0) return
    ierror = 5
    imid = (nlat+1)/2
    if ((isym==0 .and. idvw<nlat) .or. &
        (isym/=0 .and. idvw<imid)) return
    ierror = 6
    if (jdvw < nlon) return
    ierror = 7
    mmax = min(nlat, (nlon+1)/2)
    if (mdab < min(nlat, (nlon+2)/2)) return
    ierror = 8
    if (ndab < nlat) return
    ierror = 9
    idz = (mmax*(nlat+nlat-mmax+1))/2
    lzimn = idz*imid
    l1 = min(nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lwmin = l1*l2*(nlat+nlat-l1+1)+nlon+15
    if (lvhsgs < lwmin) return
    ierror = 10
    !
    !     verify unsaved work space length
    !
    mn = mmax*nlat*nt
    if (isym/=0  .and. lwork < &
        nlat*(2*nt*nlon+max(6*imid, nlon))+2*mn+nlat) return
    if (isym==0  .and. lwork < &
        imid*(2*nt*nlon+max(6*nlat, nlon))+2*mn+nlat) return
    ierror = 0
    !
    !     set work space pointers
    !
    ibr = 1
    ibi = ibr + mn
    is = ibi + mn
    iwk = is + nlat
    liwk = lwork-2*mn-nlat
    call idvgs1(nlat, nlon, isym, nt, v, w, idvw, jdvw, work(ibr), work(ibi), &
        mmax, work(is), mdab, ndab, a, b, wvhsgs, lvhsgs, work(iwk), &
        liwk, pertrb, ierror)

end subroutine idivgs



subroutine idvgs1(nlat, nlon, isym, nt, v, w, idvw, jdvw, br, bi, mmax, &
    sqnn, mdab, ndab, a, b, wsav, lwsav, wk, lwk, pertrb, ierror)
    dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt), pertrb(nt)
    dimension br(mmax, nlat, nt), bi(mmax, nlat, nt), sqnn(nlat)
    dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
    dimension wsav(lwsav), wk(lwk)
    !
    !     preset coefficient multiplyers in vector
    !
    do n=2, nlat
        fn = real(n - 1)
        sqnn(n) = sqrt(fn * (fn + 1.0))
    end do
    !
    !     compute multiple vector fields coefficients
    !
    do k=1, nt
        !
        !     set divergence field perturbation adjustment
        !
        pertrb(k) = a(1, 1, k)/(2.0 * sqrt(2.0))
        !
        !     preset br, bi to 0.0
        !
        do n=1, nlat
            do m=1, mmax
                br(m, n, k) = 0.0
                bi(m, n, k) = 0.0
            end do
        end do
        !
        !     compute m=0 coefficients
        !
        do n=2, nlat
            br(1, n, k) = -a(1, n, k)/sqnn(n)
            bi(1, n, k) = -b(1, n, k)/sqnn(n)
        end do
        !
        !     compute m>0 coefficients
        !
        do m=2, mmax
            do n=m, nlat
                br(m, n, k) = -a(m, n, k)/sqnn(n)
                bi(m, n, k) = -b(m, n, k)/sqnn(n)
            end do
        end do
    end do
    !
    !     set ityp for vector synthesis with curl=0
    !
    select case (isym)
        case (0)
            ityp = 1
        case (1)
            ityp = 4
        case (2)
            ityp = 7
    end select
    !
    !     vector sythesize br, bi into irrotational (v, w)
    !
    call vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mmax, nlat, wsav, lwsav, wk, lwk, ierror)

end subroutine idvgs1

